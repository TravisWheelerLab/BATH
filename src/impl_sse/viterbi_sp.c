/* SSE-accelerated Spliced Viterbi algorithms; full matrix.
 *
 * Probability-space DP (mul/max) with sparse rescaling.  Uses only
 * 3-nt codons (codon_lengths == 1 profile).  The P (split-codon)
 * state is vectorised via an internal pvx[SIGNAL_MEM_SIZE * Q] buffer
 * laid out in stripe order.  M/D/I states are fully vectorised.
 *
 * pvx stores probability-space values (unlike the scalar P_scores which
 * uses log-space).  This eliminates all expf()/logf() calls in the
 * P-state read and reduces the rescaling adjustment to a single
 * vectorised _mm_mul_ps.
 *
 * Contents:
 *   1. p7_Viterbi_SplicedGlobal()
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_splice.h"
#include "impl_sse.h"

/* Local nucleotide sentinel (mirrors MAX_NUC in generic_viterbi_spliced.c) */
#define SP_MAX_NUC 4

/* Probability-space constant for the P->M transition (= exp(log(4.58e-5))). */
#define TSC_P_PROB  4.58e-5f

/*****************************************************************
 * 1. p7_Viterbi_SplicedGlobal()
 *****************************************************************/

/* run_dd_passes()
 * Internal helper: run the D->D propagation loop (first pass + convergence
 * passes) for one DP row.  On entry <dpc> holds the current row, <dcv> is
 * the M*MD carry from the last stripe of the main loop (already right-shifted
 * by the caller), and tfv_dd points to the DD-transition block (om_fs->tfv +
 * 7*Q).  D[stripe 0] is zeroed (k=1 boundary) before the first pass.
 */
static void
run_dd_passes(__m128 *dpc, __m128 dcv, const __m128 *tfv_dd, int Q, int M)
{
  const __m128 *tp;
  __m128 sv, zerov = _mm_setzero_ps();
  int q, j;

  dpc[0 * p7X_NSCELLS_SP + p7X_D] = zerov;   /* D[k=1] boundary = 0 */
  tp = tfv_dd;
  for (q = 0; q < Q; q++)
    {
      __m128 *dp_d = &dpc[q * p7X_NSCELLS_SP + p7X_D];
      *dp_d = _mm_max_ps(dcv, *dp_d);
      dcv   = _mm_mul_ps(*dp_d, *tp); tp++;
    }

  if (M < 100)
    {
      for (j = 1; j < 4; j++)
        {
          dcv = esl_sse_rightshiftz_float(dcv);
          tp  = tfv_dd;
          for (q = 0; q < Q; q++)
            {
              __m128 *dp_d = &dpc[q * p7X_NSCELLS_SP + p7X_D];
              *dp_d = _mm_max_ps(dcv, *dp_d);
              dcv   = _mm_mul_ps(*dp_d, *tp); tp++;
            }
        }
    }
  else
    {
      for (j = 1; j < 4; j++)
        {
          register __m128 cv;
          dcv = esl_sse_rightshiftz_float(dcv);
          tp  = tfv_dd;
          cv  = zerov;
          for (q = 0; q < Q; q++)
            {
              __m128 *dp_d = &dpc[q * p7X_NSCELLS_SP + p7X_D];
              sv    = _mm_max_ps(dcv, *dp_d);
              cv    = _mm_or_ps(cv, _mm_cmpgt_ps(sv, *dp_d));
              *dp_d = sv;
              dcv   = _mm_mul_ps(sv, *tp); tp++;
            }
          if (!_mm_movemask_ps(cv)) break;
        }
    }
}


/* Function:  p7_Viterbi_SplicedGlobal()
 * Synopsis:  SSE-accelerated spliced global Viterbi, 3-nt codons, full matrix.
 *
 * Purpose:   Global Viterbi for a sub-region [i_start..i_end] of a nucleotide
 *            sequence against an optimised 1-codon-length profile <om_fs>,
 *            storing the full DP matrix in <ox>.  Models splice junctions
 *            via a P (split-codon) state backed by the caller-maintained
 *            an internally allocated probability-space <pvx> buffer and
 *            <signal_scores>.
 *
 *            This is the SSE accelerated equivalent of p7_GViterbi_SplicedGlobal().
 *            It requires a profile configured with codon_lengths == 1.
 *            k_start and k_end must equal 1 and om_fs->M respectively (full
 *            model; sub-range striping is not yet supported).
 *
 *            <ox> must be created with
 *              p7_omx_Create_dpf(om_fs->M, L, L, p7X_NSCELLS_SP)
 *            where L = i_end - i_start + 1.
 *
 *            On return the raw (scaled) Viterbi E value is stored in
 *              ox->xmx[L * p7X_NXCELLS + p7X_E]
 *            The log-probability score (nats) can be reconstructed as
 *              ox->totscale + logf(ox->xmx[L*p7X_NXCELLS + p7X_E]
 *                                 * om_fs->xf[p7O_E][p7O_MOVE])
 *
 * Args:      sub_dsq       - digital nucleotide sequence, i_start..i_end
 *            om_fs         - optimised 1-codon profile
 *            ox            - DP matrix (p7X_NSCELLS_SP cells per stripe)
 *            signal_scores - [p7S_SPLICE_SIGNALS] signal log-scores
 *            i_start,i_end - sequence sub-range (1-based absolute coords)
 *            k_start,k_end - model sub-range (must be 1..om_fs->M)
 *            min_intron    - minimum intron length
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile not codon_lengths==1.
 */
int
p7_Viterbi_SplicedGlobal(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_fs, P7_OMX *ox,
                         const float *signal_scores,
                         int i_start, int i_end, int k_start, int k_end, int min_intron)
{
  /* SSE register names mirror viterbi_fs.c conventions */
  register __m128  mpv3, dpv3, ipv3, ppv3;  /* right-shifted prev3 MDIP (k-1 access) */
  register __m128  sv;                       /* scratch vector                         */
  register __m128  msv;                      /* M-state best-path value                */
  register __m128  dcv;                      /* delayed D(i,q+1) carry                 */
  register __m128  xEv;                      /* E-state partial max                    */
  __m128           zerov;                    /* splatted 0.0                            */
  __m128           tsc_p_v;                  /* splatted TSC_P_PROB                     */
  __m128           sig_gtag_v;              /* splatted signal_prob[p7S_GTAG]          */
  __m128           sig_gcag_v;              /* splatted signal_prob[p7S_GCAG]          */
  __m128           sig_atac_v;              /* splatted signal_prob[p7S_ATAC]          */

  __m128  *dpc;                  /* current DP row                              */
  __m128  *dpp3;                 /* DP row i-3                                  */
  __m128  *dpd;                  /* DP row i-min_intron-3 (donor lookback)      */
  __m128  *tp;                   /* transition pointer into om_fs->tfv          */
  __m128  *pvx     = NULL;       /* vectorised P_scores: [SIGNAL_MEM_SIZE*Q]    */
  void    *pvx_mem = NULL;       /* raw allocation for pvx (needs 16-b align)   */

  int      L     = i_end - i_start + 1;
  int      M     = om_fs->M;
  int      Q     = p7O_NQF(M);

  float    xE;                   /* E-state scalar                              */
  float    xB0;                  /* B at row 0 (global entry)                   */

  /* Nucleotide rolling window (mirrors r,s,t,u,v,w,x in scalar) */
  int      r, s, t, u;           /* donor-site rolling window                  */
  int      v, w, x;              /* acceptor/codon rolling window               */

  /* Codon indices for the current row */
  int      C0;                   /* 3-nt unsplit codon index                    */
  int      C1[SP_MAX_NUC+1];     /* 1-nt-before-donor split codon (nuc1=0..4)   */
  int      C2[SP_MAX_NUC+1];     /* 2-nt-before-donor split codon (nuc3=0..4)   */
  int      nuc1, nuc3;

  /* Acceptor / donor state for current row */
  int      acc0, acc1, acc2;
  int      don0, don1, don2;
  int      don_sig;

  int      loop_end;
  int      i, q, r_lane, sig_idx;
  int      sub_i, k_local;

  float    insert_adj;           /* scale correction for I(i,k) from dpp3      */
  float    TMP_SC_prob;          /* scratch probability-space score             */

  union { __m128 v; float p[4]; } u_vec;
  int      status;

  if (om_fs->codon_lengths != 1)
    ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");

  /* Allocate pvx: SIGNAL_MEM_SIZE circular slots x Q stripes.
   * Layout: pvx[signal * Q + q] holds 4-lane prob-space P_scores for
   * stripe q at signal index signal (0..SIGNAL_MEM_SIZE-1).            */
  ESL_ALLOC(pvx_mem, sizeof(__m128) * SIGNAL_MEM_SIZE * Q + 15);
  pvx = (__m128 *) (((unsigned long int) pvx_mem + 15) & (~0xf));

  ox->M              = M;
  ox->L              = L;
  ox->has_own_scales = TRUE;
  ox->totscale       = 0.0f;

  zerov   = _mm_setzero_ps();
  tsc_p_v = _mm_set1_ps(TSC_P_PROB);

  /* Pre-compute probability-space signal scores (used in P-state read). */
  sig_gtag_v = _mm_set1_ps(expf(signal_scores[p7S_GTAG]));
  sig_gcag_v = _mm_set1_ps(expf(signal_scores[p7S_GCAG]));
  sig_atac_v = _mm_set1_ps(expf(signal_scores[p7S_ATAC]));

  /* ------------------------------------------------------------------
   * Reset pvx to 0.0 (probability-space zero = log-space -inf).
   * ------------------------------------------------------------------ */
  for (sig_idx = 0; sig_idx < SIGNAL_MEM_SIZE; sig_idx++)
    for (q = 0; q < Q; q++)
      pvx[sig_idx * Q + q] = zerov;

  /* ------------------------------------------------------------------
   * Zero-initialise rows 0..L of the DP matrix.
   * ------------------------------------------------------------------ */
  for (i = 0; i <= L; i++)
    for (q = 0; q < Q; q++)
      MMO_SP(ox->dpf[i], q) = DMO_SP(ox->dpf[i], q) =
      IMO_SP(ox->dpf[i], q) = PMO_SP(ox->dpf[i], q) = zerov;

  /* ------------------------------------------------------------------
   * Special-state scalars.  Global mode: N(0)=1, B(0)=N->B, rest=0.
   * ------------------------------------------------------------------ */
  xB0 = om_fs->xf[p7O_N][p7O_MOVE];

  for (i = 0; i <= 2; i++)
    {
      ox->xmx[i*p7X_NXCELLS + p7X_SCALE] = 1.0f;
      ox->xmx[i*p7X_NXCELLS + p7X_E]     = 0.0f;
      ox->xmx[i*p7X_NXCELLS + p7X_N]     = (i == 0) ? 1.0f : 0.0f;
      ox->xmx[i*p7X_NXCELLS + p7X_J]     = 0.0f;
      ox->xmx[i*p7X_NXCELLS + p7X_B]     = (i == 0) ? xB0  : 0.0f;
      ox->xmx[i*p7X_NXCELLS + p7X_C]     = 0.0f;
    }

  /* ------------------------------------------------------------------
   * Nucleotide rolling-window prime (rows 1 and 2 advance x twice).
   * ------------------------------------------------------------------ */
  r = s = t = u = v = w = x = -1;
  acc0 = acc1 = acc2 = -1;
  don0 = don1 = don2 = -1;

  /* row 1 advance */
  w = x;  /* = -1 */
  x = (sub_dsq[i_start] < SP_MAX_NUC) ? sub_dsq[i_start] : p7P_MAXCODONS1;

  /* row 2 advance */
  v = w;  /* = -1 */
  w = x;
  x = (sub_dsq[i_start + 1] < SP_MAX_NUC) ? sub_dsq[i_start + 1] : p7P_MAXCODONS1;

  /* ------------------------------------------------------------------
   * Row i=3: global entry at k=1 only.
   * dpp3 = dpf[0] = all zeros; M(3,k>1) = 0; I(3,k) = 0; P(3,k) = 0.
   * M(3,1) = B(0) * rfv[C0][q=0].lane[0].
   * D(3,k) propagates west-to-east from M(3,1).
   * ------------------------------------------------------------------ */
  if (L >= 3)
    {
      /* Advance codon window for row 3 */
      v = w;
      w = x;
      x = (sub_dsq[i_start + 2] < SP_MAX_NUC) ? sub_dsq[i_start + 2] : p7P_MAXCODONS1;

      C0 = p7P_CODON3_FS1(v, w, x);
      C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);

      /* Acceptor tracking for row 3 */
      acc0 = acc1;
      acc1 = acc2;
      if      (v >= SP_MAX_NUC || w >= SP_MAX_NUC) acc2 = -1;
      else if (SIGNAL(v, w) == ACCEPT_AG)           acc2 = ACCEPT_AG;
      else if (SIGNAL(v, w) == ACCEPT_AC)           acc2 = ACCEPT_AC;
      else                                          acc2 = -1;

      dpc = ox->dpf[3];   /* already zeroed */

      /* M(3,1): B(0) * emit(k=1, C0) */
      float b_emit = xB0 * p7_fs_oprofile_FGetEmission(om_fs, 1, C0);
      u_vec.v    = MMO_SP(dpc, 0);
      u_vec.p[0] = b_emit;
      MMO_SP(dpc, 0) = u_vec.v;
      xEv = MMO_SP(dpc, 0);

      /* D(3,2): main-loop carries D[q] = M[q-1]*MD[q-1].
       * Only M[q=0] is non-zero (= {b_emit, 0, 0, 0}), so only D[q=1] gets a
       * contribution: D[q=1] = M[q=0] * MD[q=0].
       * (Q >= 2 always by p7O_NQF definition.)                               */
      DMO_SP(dpc, 1) = _mm_mul_ps(MMO_SP(dpc, 0), *(om_fs->tfv + 4));  /* MD[q=0] */

      /* DD propagation (dcv=0 since no carry from last M stripe) */
      dcv = esl_sse_rightshiftz_float(zerov);
      run_dd_passes(dpc, dcv, om_fs->tfv + 7*Q, Q, M);

      /* Collect D into xEv */
      for (q = 0; q < Q; q++)
        xEv = _mm_max_ps(xEv, DMO_SP(dpc, q));
      esl_sse_hmax_ps(xEv, &xE);

      /* Sparse rescaling */
      if (xE > 1.0e4f)
        {
          float  sfv_f = 1.0f / xE;
          __m128 sfv   = _mm_set1_ps(sfv_f);
          for (q = 0; q < Q; q++)
            {
              MMO_SP(dpc, q) = _mm_mul_ps(MMO_SP(dpc, q), sfv);
              DMO_SP(dpc, q) = _mm_mul_ps(DMO_SP(dpc, q), sfv);
              /* IMO_SP and PMO_SP are zero, skip */
            }
          ox->xmx[3*p7X_NXCELLS + p7X_SCALE] = xE;
          ox->totscale += logf(xE);
          xE = 1.0f;
        }
      else ox->xmx[3*p7X_NXCELLS + p7X_SCALE] = 1.0f;

      ox->xmx[3*p7X_NXCELLS + p7X_E] = xE;
    }

  /* ------------------------------------------------------------------
   * Phase 2: rows i = 4 .. min(L, min_intron+2)
   * M/D/I states only; P = 0.  No global entry (dpp3 non-zero).
   * ------------------------------------------------------------------ */
  loop_end = ESL_MIN(L, min_intron + 2);

  for (i = 4; i <= loop_end; i++)
    {
      sub_i = i_start + i - 1;

      v = w;
      w = x;
      x = (sub_dsq[sub_i] < SP_MAX_NUC) ? sub_dsq[sub_i] : p7P_MAXCODONS1;

      C0 = p7P_CODON3_FS1(v, w, x);
      C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);

      acc0 = acc1;
      acc1 = acc2;
      if      (v >= SP_MAX_NUC || w >= SP_MAX_NUC) acc2 = -1;
      else if (SIGNAL(v, w) == ACCEPT_AG)           acc2 = ACCEPT_AG;
      else if (SIGNAL(v, w) == ACCEPT_AC)           acc2 = ACCEPT_AC;
      else                                          acc2 = -1;

      dpc  = ox->dpf[i];
      dpp3 = ox->dpf[i - 3];

      insert_adj = 1.0f
                 / (ox->xmx[(i-2)*p7X_NXCELLS + p7X_SCALE]
                 *  ox->xmx[(i-1)*p7X_NXCELLS + p7X_SCALE]);
      register __m128 adj_v = _mm_set1_ps(insert_adj);

      mpv3 = esl_sse_rightshiftz_float(MMO_SP(dpp3, Q-1));
      dpv3 = esl_sse_rightshiftz_float(DMO_SP(dpp3, Q-1));
      ipv3 = esl_sse_rightshiftz_float(IMO_SP(dpp3, Q-1));

      tp  = om_fs->tfv;
      dcv = zerov;
      xEv = zerov;

      for (q = 0; q < Q; q++)
        {
          tp++;  /* skip BM */
          sv  =                _mm_mul_ps(mpv3, *tp); tp++;  /* MM */
          sv  = _mm_max_ps(sv, _mm_mul_ps(ipv3, *tp)); tp++; /* IM */
          sv  = _mm_max_ps(sv, _mm_mul_ps(dpv3, *tp)); tp++; /* DM */
          msv = _mm_mul_ps(sv, om_fs->rfv[C0][q]);
          xEv = _mm_max_ps(xEv, msv);

          mpv3 = MMO_SP(dpp3, q);
          dpv3 = DMO_SP(dpp3, q);
          ipv3 = IMO_SP(dpp3, q);

          MMO_SP(dpc, q) = msv;
          DMO_SP(dpc, q) = dcv;
          dcv = _mm_mul_ps(msv, *tp); tp++;   /* MD */

          sv          = _mm_mul_ps(_mm_mul_ps(MMO_SP(dpp3, q), adj_v), *tp); tp++; /* MI */
          IMO_SP(dpc, q) = _mm_max_ps(sv,
                           _mm_mul_ps(_mm_mul_ps(IMO_SP(dpp3, q), adj_v), *tp)); tp++; /* II */

          PMO_SP(dpc, q) = zerov;   /* no P state in Phase 2 */
        }

      dcv = esl_sse_rightshiftz_float(dcv);
      run_dd_passes(dpc, dcv, om_fs->tfv + 7*Q, Q, M);

      for (q = 0; q < Q; q++)
        xEv = _mm_max_ps(xEv, DMO_SP(dpc, q));
      esl_sse_hmax_ps(xEv, &xE);

      if (xE > 1.0e4f)
        {
          float  sfv_f = 1.0f / xE;
          __m128 sfv   = _mm_set1_ps(sfv_f);
          for (q = 0; q < Q; q++)
            {
              MMO_SP(dpc, q) = _mm_mul_ps(MMO_SP(dpc, q), sfv);
              DMO_SP(dpc, q) = _mm_mul_ps(DMO_SP(dpc, q), sfv);
              IMO_SP(dpc, q) = _mm_mul_ps(IMO_SP(dpc, q), sfv);
              PMO_SP(dpc, q) = _mm_mul_ps(PMO_SP(dpc, q), sfv);
            }
          ox->xmx[i*p7X_NXCELLS + p7X_SCALE] = xE;
          ox->totscale += logf(xE);
          xE = 1.0f;
        }
      else ox->xmx[i*p7X_NXCELLS + p7X_SCALE] = 1.0f;

      ox->xmx[i*p7X_NXCELLS + p7X_E] = xE;
    }  /* end Phase 2 */

  /* ------------------------------------------------------------------
   * Initialise donor-lookback rolling window r,s,t,u.
   * Mirrors generic_viterbi_spliced.c lines 231-238.
   * ------------------------------------------------------------------ */
  s = (sub_dsq[i_start]   < SP_MAX_NUC) ? (int)sub_dsq[i_start]   : p7P_MAXCODONS1;
  t = (sub_dsq[i_start+1] < SP_MAX_NUC) ? (int)sub_dsq[i_start+1] : p7P_MAXCODONS1;
  u = (sub_dsq[i_start+2] < SP_MAX_NUC) ? (int)sub_dsq[i_start+2] : p7P_MAXCODONS1;

  /* ------------------------------------------------------------------
   * Phase 3: rows i = min_intron+3 .. L
   * Full M/D/I/P recursion + donor-site k-loop.
   * ------------------------------------------------------------------ */
  for (i = min_intron + 3; i <= L; i++)
    {
      sub_i = i_start + i - 1;

      /* Advance donor-site window */
      r = s; s = t; t = u;
      u = (sub_dsq[sub_i - min_intron + 1] < SP_MAX_NUC)
            ? (int)sub_dsq[sub_i - min_intron + 1] : p7P_MAXCODONS1;

      /* Advance codon window */
      v = w; w = x;
      x = (sub_dsq[sub_i] < SP_MAX_NUC) ? (int)sub_dsq[sub_i] : p7P_MAXCODONS1;

      C0 = p7P_CODON3_FS1(v, w, x);
      C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);

      for (nuc1 = 0; nuc1 < SP_MAX_NUC; nuc1++)
        C1[nuc1] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, w, x), p7P_DEGEN1_C);
      C1[SP_MAX_NUC] = p7P_MINIDX(p7P_CODON3_FS1(p7P_MAXCODONS1, w, x), p7P_DEGEN1_C);

      nuc3 = ESL_MIN(x, SP_MAX_NUC);

      /* Acceptor tracking */
      acc0 = acc1; acc1 = acc2;
      if      (SIGNAL(v, w) == ACCEPT_AG) acc2 = ACCEPT_AG;
      else if (SIGNAL(v, w) == ACCEPT_AC) acc2 = ACCEPT_AC;
      else                                acc2 = -1;

      /* Donor tracking */
      don0 = don1; don1 = don2;
      if      (SIGNAL(t, u) == DONOR_GT) don2 = DONOR_GT;
      else if (SIGNAL(t, u) == DONOR_GC) don2 = DONOR_GC;
      else if (SIGNAL(t, u) == DONOR_AT) don2 = DONOR_AT;
      else                               don2 = -1;

      dpc  = ox->dpf[i];
      dpp3 = ox->dpf[i - 3];
      dpd  = ox->dpf[i - min_intron - 3];

      insert_adj = 1.0f
                 / (ox->xmx[(i-2)*p7X_NXCELLS + p7X_SCALE]
                 *  ox->xmx[(i-1)*p7X_NXCELLS + p7X_SCALE]);
      register __m128 adj_v = _mm_set1_ps(insert_adj);

      mpv3 = esl_sse_rightshiftz_float(MMO_SP(dpp3, Q-1));
      dpv3 = esl_sse_rightshiftz_float(DMO_SP(dpp3, Q-1));
      ipv3 = esl_sse_rightshiftz_float(IMO_SP(dpp3, Q-1));
      ppv3 = esl_sse_rightshiftz_float(PMO_SP(dpp3, Q-1));

      tp  = om_fs->tfv;
      dcv = zerov;
      xEv = zerov;

      for (q = 0; q < Q; q++)
        {
          /* M(i,k) = max(M(i-3,k-1)*MM, I(i-3,k-1)*IM,
           *              D(i-3,k-1)*DM, P(i-3,k-1)*TSC_P) * rfv[C0] */
          tp++;  /* skip BM */
          sv  =                _mm_mul_ps(mpv3, *tp); tp++;  /* MM */
          sv  = _mm_max_ps(sv, _mm_mul_ps(ipv3, *tp)); tp++; /* IM */
          sv  = _mm_max_ps(sv, _mm_mul_ps(dpv3, *tp)); tp++; /* DM */
          sv  = _mm_max_ps(sv, _mm_mul_ps(ppv3, tsc_p_v));   /* P*TSC_P */
          msv = _mm_mul_ps(sv, om_fs->rfv[C0][q]);
          xEv = _mm_max_ps(xEv, msv);

          mpv3 = MMO_SP(dpp3, q);
          dpv3 = DMO_SP(dpp3, q);
          ipv3 = IMO_SP(dpp3, q);
          ppv3 = PMO_SP(dpp3, q);

          MMO_SP(dpc, q) = msv;
          DMO_SP(dpc, q) = dcv;
          dcv = _mm_mul_ps(msv, *tp); tp++;   /* MD */

          sv          = _mm_mul_ps(_mm_mul_ps(MMO_SP(dpp3, q), adj_v), *tp); tp++; /* MI */
          IMO_SP(dpc, q) = _mm_max_ps(sv,
                           _mm_mul_ps(_mm_mul_ps(IMO_SP(dpp3, q), adj_v), *tp)); tp++; /* II */

          /* P state: vectorised probability-space lookup.
           * pvx[signal*Q + q] holds the 4-lane prob-space P_scores for
           * stripe q; signal vectors are pre-multiplied into sig_*_v.
           * Padding lanes (k_local >= M) are zero in pvx, contributing 0. */
          if (acc0 < 0 && acc1 < 0 && acc2 < 0)
            {
              PMO_SP(dpc, q) = zerov;
            }
          else
            {
              __m128 pv        = zerov;
              __m128 emit_c0_v = om_fs->rfv[C0][q];

              /* acc0: SSX0 * signal_prob * emit_C0 */
              if (acc0 == ACCEPT_AG)
                {
                  __m128 sg = _mm_mul_ps(pvx[p7S_GTAG * Q + q], sig_gtag_v);
                  __m128 sc = _mm_mul_ps(pvx[p7S_GCAG * Q + q], sig_gcag_v);
                  pv = _mm_max_ps(pv, _mm_mul_ps(_mm_max_ps(sg, sc), emit_c0_v));
                }
              else if (acc0 == ACCEPT_AC)
                {
                  __m128 sa = _mm_mul_ps(pvx[p7S_ATAC * Q + q], sig_atac_v);
                  pv = _mm_max_ps(pv, _mm_mul_ps(sa, emit_c0_v));
                }

              /* acc1: SSX1(nuc1) * signal_prob * emit_C1[nuc1] */
              if (acc1 == ACCEPT_AG)
                {
                  for (nuc1 = 0; nuc1 <= SP_MAX_NUC; nuc1++)
                    {
                      int base = SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS;
                      __m128 sg = _mm_mul_ps(pvx[(base + p7S_GTAG) * Q + q], sig_gtag_v);
                      __m128 sc = _mm_mul_ps(pvx[(base + p7S_GCAG) * Q + q], sig_gcag_v);
                      pv = _mm_max_ps(pv, _mm_mul_ps(_mm_max_ps(sg, sc), om_fs->rfv[C1[nuc1]][q]));
                    }
                }
              else if (acc1 == ACCEPT_AC)
                {
                  for (nuc1 = 0; nuc1 <= SP_MAX_NUC; nuc1++)
                    {
                      int base = SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS;
                      __m128 sa = _mm_mul_ps(pvx[(base + p7S_ATAC) * Q + q], sig_atac_v);
                      pv = _mm_max_ps(pv, _mm_mul_ps(sa, om_fs->rfv[C1[nuc1]][q]));
                    }
                }

              /* acc2: SSX2(nuc3) * signal_prob (emission baked in at write time) */
              if (acc2 == ACCEPT_AG)
                {
                  int base = SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS;
                  __m128 sg = _mm_mul_ps(pvx[(base + p7S_GTAG) * Q + q], sig_gtag_v);
                  __m128 sc = _mm_mul_ps(pvx[(base + p7S_GCAG) * Q + q], sig_gcag_v);
                  pv = _mm_max_ps(pv, _mm_max_ps(sg, sc));
                }
              else if (acc2 == ACCEPT_AC)
                {
                  int base = SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS;
                  pv = _mm_max_ps(pv, _mm_mul_ps(pvx[(base + p7S_ATAC) * Q + q], sig_atac_v));
                }

              PMO_SP(dpc, q) = pv;
              xEv = _mm_max_ps(xEv, pv);
            }
        } /* end q loop */

      dcv = esl_sse_rightshiftz_float(dcv);
      run_dd_passes(dpc, dcv, om_fs->tfv + 7*Q, Q, M);

      for (q = 0; q < Q; q++)
        xEv = _mm_max_ps(xEv, DMO_SP(dpc, q));
      esl_sse_hmax_ps(xEv, &xE);

      if (xE > 1.0e4f)
        {
          float  sfv_f = 1.0f / xE;
          __m128 sfv   = _mm_set1_ps(sfv_f);
          for (q = 0; q < Q; q++)
            {
              MMO_SP(dpc, q) = _mm_mul_ps(MMO_SP(dpc, q), sfv);
              DMO_SP(dpc, q) = _mm_mul_ps(DMO_SP(dpc, q), sfv);
              IMO_SP(dpc, q) = _mm_mul_ps(IMO_SP(dpc, q), sfv);
              PMO_SP(dpc, q) = _mm_mul_ps(PMO_SP(dpc, q), sfv);
            }
          /* Keep pvx in sync: pvx is probability-space so divide by xE
           * (same as scaling M/D/I/P cells).                             */
          for (sig_idx = 0; sig_idx < SIGNAL_MEM_SIZE; sig_idx++)
            for (q = 0; q < Q; q++)
              pvx[sig_idx * Q + q] = _mm_mul_ps(pvx[sig_idx * Q + q], sfv);
          ox->xmx[i*p7X_NXCELLS + p7X_SCALE] = xE;
          ox->totscale += logf(xE);
          xE = 1.0f;
        }
      else ox->xmx[i*p7X_NXCELLS + p7X_SCALE] = 1.0f;

      ox->xmx[i*p7X_NXCELLS + p7X_E] = xE;

      /* ------------------------------------------------------------------
       * Donor-site k-loop: update pvx for future P-state reads.
       * Reads from dpd = dpf[i - min_intron - 3].
       * TMP_SC_prob = max(M(dpd,k-1), D(dpd,k-1)) — probability-space.
       * pvx stores probability-space; no logf() conversion needed.
       * ------------------------------------------------------------------ */
      if (don2 >= 0)
        {
          /* C2-type: 2 nt before donor; emission baked into pvx at write time. */
          for (nuc3 = 0; nuc3 < SP_MAX_NUC; nuc3++)
            C2[nuc3] = p7P_MINIDX(p7P_CODON3_FS1(r, s, nuc3), p7P_DEGEN1_C);
          C2[SP_MAX_NUC] = p7P_MINIDX(p7P_CODON3_FS1(r, s, p7P_MAXCODONS1), p7P_DEGEN1_C);

          don_sig = (don2 == DONOR_GT) ? p7S_GTAG
                  : (don2 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;

          for (q = 0; q < Q; q++)
            for (r_lane = 0; r_lane < 4; r_lane++)
              {
                k_local = q + r_lane * Q + 1;   /* 1-based */
                if (k_local < 2 || k_local >= M) continue;

                int km1 = k_local - 1;
                int q2  = (km1 - 1) % Q;
                int r2  = (km1 - 1) / Q;

                u_vec.v = MMO_SP(dpd, q2); TMP_SC_prob  = u_vec.p[r2];
                u_vec.v = DMO_SP(dpd, q2); TMP_SC_prob  = ESL_MAX(TMP_SC_prob, u_vec.p[r2]);

                if (TMP_SC_prob <= 0.0f) continue;
                /* Probability-space: no logf needed */

                for (nuc3 = 0; nuc3 <= SP_MAX_NUC; nuc3++)
                  {
                    float new_sc = TMP_SC_prob * p7_fs_oprofile_FGetEmission(om_fs, k_local, C2[nuc3]);
                    int   off    = SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + don_sig;
                    u_vec.v = pvx[off * Q + q];
                    if (new_sc > u_vec.p[r_lane]) u_vec.p[r_lane] = new_sc;
                    pvx[off * Q + q] = u_vec.v;
                  }
              }
        } /* end don2 */

      if (don1 >= 0)
        {
          nuc1 = ESL_MIN(r, SP_MAX_NUC);
          don_sig = (don1 == DONOR_GT) ? p7S_GTAG
                  : (don1 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;

          for (q = 0; q < Q; q++)
            for (r_lane = 0; r_lane < 4; r_lane++)
              {
                k_local = q + r_lane * Q + 1;
                if (k_local < 2 || k_local >= M) continue;

                int km1 = k_local - 1;
                int q2  = (km1 - 1) % Q;
                int r2  = (km1 - 1) / Q;

                u_vec.v = MMO_SP(dpd, q2); TMP_SC_prob  = u_vec.p[r2];
                u_vec.v = DMO_SP(dpd, q2); TMP_SC_prob  = ESL_MAX(TMP_SC_prob, u_vec.p[r2]);

                if (TMP_SC_prob <= 0.0f) continue;

                int off = SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + don_sig;
                u_vec.v = pvx[off * Q + q];
                if (TMP_SC_prob > u_vec.p[r_lane]) u_vec.p[r_lane] = TMP_SC_prob;
                pvx[off * Q + q] = u_vec.v;
              }
        } /* end don1 */

      if (don0 >= 0)
        {
          don_sig = (don0 == DONOR_GT) ? p7S_GTAG
                  : (don0 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;

          for (q = 0; q < Q; q++)
            for (r_lane = 0; r_lane < 4; r_lane++)
              {
                k_local = q + r_lane * Q + 1;
                if (k_local < 2 || k_local >= M) continue;

                int km1 = k_local - 1;
                int q2  = (km1 - 1) % Q;
                int r2  = (km1 - 1) / Q;

                u_vec.v = MMO_SP(dpd, q2); TMP_SC_prob  = u_vec.p[r2];
                u_vec.v = DMO_SP(dpd, q2); TMP_SC_prob  = ESL_MAX(TMP_SC_prob, u_vec.p[r2]);

                if (TMP_SC_prob <= 0.0f) continue;

                u_vec.v = pvx[don_sig * Q + q];
                if (TMP_SC_prob > u_vec.p[r_lane]) u_vec.p[r_lane] = TMP_SC_prob;
                pvx[don_sig * Q + q] = u_vec.v;
              }
        } /* end don0 */

    } /* end Phase 3 */

  /* ------------------------------------------------------------------
   * Final exit: E(L) = max(M(L,M), D(L,M)).
   * M is model position M; in striped layout: stripe q_last, lane r_last.
   * ------------------------------------------------------------------ */
  {
    int q_last = (M - 1) % Q;
    int r_last = (M - 1) / Q;
    __m128 *dpL = ox->dpf[L];

    u_vec.v = MMO_SP(dpL, q_last);
    float mLM = u_vec.p[r_last];
    u_vec.v = DMO_SP(dpL, q_last);
    float dLM = u_vec.p[r_last];
    xE = ESL_MAX(mLM, dLM);

    ox->xmx[L*p7X_NXCELLS + p7X_E] = xE;
    ox->xmx[L*p7X_NXCELLS + p7X_C] = xE * om_fs->xf[p7O_E][p7O_MOVE];
  }

  free(pvx_mem);
  return eslOK;

 ERROR:
  if (pvx_mem) free(pvx_mem);
  return status;
}
/*------------------ end, p7_Viterbi_SplicedGlobal() ---------------*/



/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/
#ifdef p7VITERBI_SP_BENCHMARK
/*
   gcc -g -O3 -msse2 -std=gnu99 -o viterbi_sp_benchmark -I.. -L.. -I../../easel -L../../easel -Dp7VITERBI_SP_BENCHMARK viterbi_sp.c -lhmmer -leasel -lm

   ./viterbi_sp_benchmark <hmmfile>
   ./viterbi_sp_benchmark -c -N100 <hmmfile>   compare scores to generic spliced Viterbi
*/
#include "p7_config.h"

#include <math.h>
#include <stdio.h>

#include "easel.h"
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS benchmark_options[] = {
  /* name    type          default  env  range  toggles reqs incomp  help                                             docgroup */
  { "-h",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",                0 },
  { "-c",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "compare scores to generic implementation (debug)",    0 },
  { "-s",  eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                       0 },
  { "-N",  eslARG_INT,    "100", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                        0 },
  { "-I",  eslARG_INT,    "200", NULL, "n>0", NULL,  NULL, NULL, "length of simulated intron (excl. GT..AG signals)",   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char benchmark_usage[]  = "[-options] <hmmfile>";
static char benchmark_banner[] = "benchmark driver for SSE spliced Viterbi";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go           = p7_CreateDefaultApp(benchmark_options, 1, argc, argv, benchmark_banner, benchmark_usage);
  char           *hmmfile      = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w            = esl_stopwatch_Create();
  ESL_RANDOMNESS *r            = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA        = NULL;
  ESL_ALPHABET   *abcDNA       = esl_alphabet_Create(eslDNA);
  P7_HMMFILE     *hfp          = NULL;
  P7_HMM         *hmm          = NULL;
  P7_BG          *bgAA         = NULL;
  P7_PROFILE     *gm           = NULL;
  P7_FS_PROFILE  *gm_tr        = NULL;
  P7_FS_OPROFILE *om_fs        = NULL;
  ESL_GENCODE    *gcode        = NULL;
  P7_CODONTABLE  *codon_table  = NULL;
  ESL_SQ         *sq           = NULL;
  SPLICE_PIPELINE *pli         = NULL;
  P7_OMX         *ox           = NULL;
  int             N            = esl_opt_GetInteger(go, "-N");
  int             I            = esl_opt_GetInteger(go, "-I");
  int             intron_total = I + 4;
  ESL_DSQ        *dsq          = NULL;
  int             i, j, k, L_amino, L_dna_total;
  int64_t         total_cells;
  float           sc, sc2;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)          != eslOK) p7_Fail("Failed to read HMM");

  gcode       = esl_gencode_Create(abcDNA, abcAA);
  bgAA        = p7_bg_Create(abcAA);
  gm          = p7_profile_Create(hmm->M, abcAA);
  gm_tr       = p7_profile_fs_Create(hmm->M, abcAA, 1);
  codon_table = p7_codontable_Create(gcode);
  sq          = esl_sq_CreateDigital(abcAA);

  p7_ProfileConfig   (hmm, bgAA,        gm,    hmm->M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, hmm->M, p7_UNILOCAL);

  om_fs = p7_fs_oprofile_Create(hmm->M, abcAA, 1);
  p7_fs_oprofile_Convert(gm_tr, om_fs);

  ox  = p7_omx_Create_dpf(hmm->M, hmm->M, hmm->M, p7X_NSCELLS_SP);
  pli = p7_splicepipeline_Create(NULL, hmm->M, hmm->M * 3);
  p7_splicescores_GrowTo(pli->splice_scores, hmm->M);

  /* Baseline: time to generate sequences alone */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, NULL);
      L_amino     = sq->n;
      L_dna_total = L_amino * 3 + intron_total;
      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) p7_Fail("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;
      j = 1;
      for (k = 1; k <= L_amino / 2; k++) { p7_codontable_GetCodon(codon_table, r, sq->dsq[k], dsq + j); j += 3; }
      dsq[j++] = 2;  /* G */
      dsq[j++] = 3;  /* T */
      for (k = 0; k < I; k++) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = 0;  /* A */
      dsq[j++] = 2;  /* G */
      for (k = L_amino / 2 + 1; k <= L_amino; k++) { p7_codontable_GetCodon(codon_table, r, sq->dsq[k], dsq + j); j += 3; }
    }
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark */
  total_cells = 0;
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, NULL);
      L_amino     = sq->n;
      L_dna_total = L_amino * 3 + intron_total;
      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) p7_Fail("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;
      j = 1;
      for (k = 1; k <= L_amino / 2; k++) { p7_codontable_GetCodon(codon_table, r, sq->dsq[k], dsq + j); j += 3; }
      dsq[j++] = 2;  /* G */
      dsq[j++] = 3;  /* T */
      for (k = 0; k < I; k++) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = 0;  /* A */
      dsq[j++] = 2;  /* G */
      for (k = L_amino / 2 + 1; k <= L_amino; k++) { p7_codontable_GetCodon(codon_table, r, sq->dsq[k], dsq + j); j += 3; }

      p7_fs_ReconfigLength(gm_tr, L_dna_total / 3);
      p7_fs_oprofile_ReconfigLength(om_fs, L_dna_total / 3);
      p7_omx_GrowTo_dpf(ox, hmm->M, L_dna_total, L_dna_total);

      if (esl_opt_GetBoolean(go, "-c"))
        {
          p7_gmx_GrowTo(pli->vit, hmm->M, L_dna_total, L_dna_total);
          p7_GViterbi_SplicedGlobal(dsq, gm_tr, pli->vit,
                                    pli->splice_scores->P_scores, pli->splice_scores->signal_scores,
                                    1, L_dna_total, 1, hmm->M, pli->min_intron);
          sc2 = pli->vit->xmx[L_dna_total * p7G_NXCELLS + p7G_C];
        }

      p7_Viterbi_SplicedGlobal(dsq, om_fs, ox,
                               pli->splice_scores->signal_scores,
                               1, L_dna_total, 1, hmm->M, pli->min_intron);

      if (esl_opt_GetBoolean(go, "-c"))
        {
          sc = (ox->xmx[L_dna_total * p7X_NXCELLS + p7X_C] > 0.0f)
               ? ox->totscale + logf(ox->xmx[L_dna_total * p7X_NXCELLS + p7X_C])
               : -eslINFINITY;
          printf("%.4f %.4f\n", sc, sc2);
        }

      total_cells += (int64_t) L_dna_total * hmm->M;
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) total_cells * 1e-6 / bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M          = %d\n",   hmm->M);
  printf("# N          = %d\n",   N);
  printf("# I          = %d\n",   I);
  printf("# %.1f Mc/s\n", Mcs);

  if (dsq != NULL) free(dsq);
  p7_omx_Destroy(ox);
  p7_splicepipeline_Destroy(pli);
  esl_sq_Destroy(sq);
  p7_codontable_Destroy(codon_table);
  p7_fs_oprofile_Destroy(om_fs);
  p7_profile_fs_Destroy(gm_tr);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bgAA);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7VITERBI_SP_BENCHMARK*/
/*----------------- end, benchmark driver -----------------------*/



/*****************************************************************
 * 4. Unit tests.
 *****************************************************************/
#ifdef p7VITERBI_SP_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_viterbi_sp()
 *
 * For N random HMMs, a profile-emitted amino-acid sequence is reverse-
 * translated into a DNA sequence with a simulated GT..AG intron inserted
 * at the midpoint (matching the construction in the generic_viterbi_spliced
 * unit test).  Both the scalar p7_GViterbi_SplicedGlobal() and the SSE
 * p7_Viterbi_SplicedGlobal() are run on the same sequence.
 *
 * The final C-state log-probability is extracted from each:
 *   scalar:  gx->xmx[L * p7G_NXCELLS + p7G_C]                (log-space)
 *   SSE:     ox->totscale + logf(ox->xmx[L * p7X_NXCELLS + p7X_C])
 *
 * Test passes when the two scores agree to within 0.001 nats.
 * If both implementations agree the sequence has no valid path
 * (C = -inf), the sequence is skipped rather than counted as a pass.
 */
static void
utest_viterbi_sp(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
                 ESL_GENCODE *gcode, P7_BG *bgAA, P7_CODONTABLE *codon_table,
                 int M, int N, int intron_len)
{
  char            *msg         = "utest_viterbi_sp failed";
  P7_HMM          *hmm         = NULL;
  P7_PROFILE      *gm          = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE   *gm_tr       = p7_profile_fs_Create(M, abcAA, 1);
  P7_FS_OPROFILE  *om_fs       = p7_fs_oprofile_Create(M, abcAA, 1);
  ESL_SQ          *sq          = esl_sq_CreateDigital(abcAA);
  P7_OMX          *ox          = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS_SP);
  SPLICE_PIPELINE *pli         = NULL;
  ESL_DSQ         *dsq         = NULL;
  int              intron_total = intron_len + 4;   /* GT + intron_len random + AG */
  int              L_amino, L_dna_total;
  float            gsc, osc;
  int              i, j, min_intron;

  if (!gm || !gm_tr || !om_fs || !sq || !ox) esl_fatal(msg);

  /* Pipeline creates pli->vit as a P7_GMX with the splice-specific 4-cell
   * layout (M/D/I/P) needed by p7_GViterbi_SplicedGlobal().              */
  pli = p7_splicepipeline_Create(NULL, M, M * 3);
  min_intron = pli->min_intron;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig   (hmm, bgAA, gm,    M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, M, p7_UNILOCAL);
  p7_fs_oprofile_Convert(gm_tr, om_fs);

  while (N--)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, NULL);
      L_amino     = sq->n;
      L_dna_total = L_amino * 3 + intron_total;

      /* Build the nucleotide sequence: exon1 + GT+intron+AG + exon2 */
      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) esl_fatal("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;

      j = 1;
      for (i = 1; i <= L_amino / 2; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }
      dsq[j++] = 2;  /* G — eslDNA: A=0 C=1 G=2 T=3 */
      dsq[j++] = 3;  /* T */
      for (i = 0; i < intron_len; i++) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = 0;  /* A */
      dsq[j++] = 2;  /* G */
      for (i = L_amino / 2 + 1; i <= L_amino; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }

      p7_fs_ReconfigLength(gm_tr, L_dna_total / 3);
      p7_fs_oprofile_ReconfigLength(om_fs, L_dna_total / 3);

      p7_gmx_GrowTo(pli->vit, M, L_dna_total, L_dna_total);
      p7_omx_GrowTo_dpf(ox, M, L_dna_total, L_dna_total);
      p7_splicescores_GrowTo(pli->splice_scores, M);

      /* --- Scalar Viterbi --- */
      p7_GViterbi_SplicedGlobal(dsq, gm_tr, pli->vit,
                                pli->splice_scores->P_scores,
                                pli->splice_scores->signal_scores,
                                1, L_dna_total, 1, M, min_intron);
      gsc = pli->vit->xmx[L_dna_total * p7G_NXCELLS + p7G_C];

      /* --- SSE Viterbi --- */
      p7_Viterbi_SplicedGlobal(dsq, om_fs, ox,
                               pli->splice_scores->signal_scores,
                               1, L_dna_total, 1, M, min_intron);

      if (ox->xmx[L_dna_total * p7X_NXCELLS + p7X_C] <= 0.0f)
        {
          /* SSE reports no viable path; scalar must agree */
          if (gsc != -eslINFINITY)
            esl_fatal("%s: SSE C=0 but scalar C=%.4f", msg, gsc);
          continue;   /* both -inf; skip score comparison */
        }
      osc = ox->totscale + logf(ox->xmx[L_dna_total * p7X_NXCELLS + p7X_C]);

      if (gsc == -eslINFINITY)
        esl_fatal("%s: scalar C=-inf but SSE C=%.4f", msg, osc);

      if (fabs(gsc - osc) > 0.001f) {
        p7_omx_DumpSP(stdout, ox);
        p7_gmx_Dump(stdout, pli->vit, P7_DEFAULT);
        esl_fatal("%s: scalar C=%.4f  SSE C=%.4f  diff=%.4f nats",
                  msg, gsc, osc, fabs(gsc - osc));
      }
    }

  if (dsq != NULL) free(dsq);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_tr);
  p7_fs_oprofile_Destroy(om_fs);
  p7_omx_Destroy(ox);
  p7_splicepipeline_Destroy(pli);
}
#endif /*p7VITERBI_SP_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/



/*****************************************************************
 * 5. Test driver.
 *****************************************************************/
#ifdef p7VITERBI_SP_TESTDRIVE
/*
  gcc -g -Wall -msse2 -std=gnu99 -o viterbi_sp_utest -I.. -L.. -I../../easel -L../../easel -Dp7VITERBI_SP_TESTDRIVE viterbi_sp.c -lhmmer -leasel -lm
  ./viterbi_sp_utest
*/
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name   type          default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h", eslARG_NONE,   FALSE,  NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",           0 },
  { "-s", eslARG_INT,     "42",  NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M", eslARG_INT,     "72",  NULL, NULL, NULL, NULL, NULL, "size of random models to sample",                0 },
  { "-N", eslARG_INT,     "20",  NULL, NULL, NULL, NULL, NULL, "number of random sequences to sample",           0 },
  { "-I", eslARG_INT,    "500",  NULL, "n>0",NULL, NULL, NULL, "simulated intron length (random nt between GT..AG)", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE p7_Viterbi_SplicedGlobal()";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go    = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r     = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA  = esl_alphabet_Create(eslAMINO);
  ESL_ALPHABET   *abcDNA = esl_alphabet_Create(eslDNA);
  ESL_GENCODE    *gcode  = esl_gencode_Create(abcDNA, abcAA);
  P7_CODONTABLE  *ct     = p7_codontable_Create(gcode);
  P7_BG          *bgAA   = p7_bg_Create(abcAA);
  int             M      = esl_opt_GetInteger(go, "-M");
  int             N      = esl_opt_GetInteger(go, "-N");
  int             I      = esl_opt_GetInteger(go, "-I");

  p7_FLogsumInit();

  utest_viterbi_sp(r, abcAA, abcDNA, gcode, bgAA, ct, M, N, I);

  p7_bg_Destroy(bgAA);
  p7_codontable_Destroy(ct);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  printf("All tests passed.\n");
  return eslOK;
}
#endif /*p7VITERBI_SP_TESTDRIVE*/
/*----------------- end, test driver ----------------------------*/
