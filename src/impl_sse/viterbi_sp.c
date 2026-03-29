/* SSE-accelerated Spliced Viterbi algorithms; full matrix.
 *
 * Probability-space DP (mul/max) with sparse rescaling.  Uses only
 * 3-nt codons (codon_lengths == 1 profile).  The P (split-codon)
 * state is computed per-lane (scalar) because P_scores[][SIGNAL_MEM_SIZE]
 * is indexed by model position, not by stripe.  M/D/I states are fully
 * vectorised.
 *
 * P_scores is kept in log-space (as in the scalar implementation);
 * conversions between log- and prob-space occur at the donor-update
 * and P-state-read boundaries.
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
 *            log-space <P_scores> table and <signal_scores>.
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
 *            P_scores      - [0..M-1][SIGNAL_MEM_SIZE] splice-site log-scores
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
                         float **P_scores, const float *signal_scores,
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

  __m128  *dpc;                  /* current DP row                              */
  __m128  *dpp3;                 /* DP row i-3                                  */
  __m128  *dpd;                  /* DP row i-min_intron-3 (donor lookback)      */
  __m128  *tp;                   /* transition pointer into om_fs->tfv          */

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
  int      i, q, r_lane;
  int      sub_i, k_local;

  float    insert_adj;           /* scale correction for I(i,k) from dpp3      */
  float    TMP_SC_prob;          /* scratch probability-space score             */
  float    TMP_SC_log;           /* scratch log-space score                     */

  /* P-state per-lane scratch */
  float    p_lane[4];
  union { __m128 v; float p[4]; } u_vec;

  if (om_fs->codon_lengths != 1)
    ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");

  ox->M              = M;
  ox->L              = L;
  ox->has_own_scales = TRUE;
  ox->totscale       = 0.0f;

  zerov   = _mm_setzero_ps();
  tsc_p_v = _mm_set1_ps(TSC_P_PROB);

  /* ------------------------------------------------------------------
   * Reset P_scores to -inf.
   * ------------------------------------------------------------------ */
  for (k_local = 0; k_local < M; k_local++)
    for (i = 0; i < SIGNAL_MEM_SIZE; i++)
      P_scores[k_local][i] = -eslINFINITY;

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

          /* P state: per-lane scalar lookup (log->prob conversion) */
          if (acc0 < 0 && acc1 < 0 && acc2 < 0)
            {
              PMO_SP(dpc, q) = zerov;
            }
          else
            {
              for (r_lane = 0; r_lane < 4; r_lane++)
                {
                  k_local = q + r_lane * Q + 1;   /* 1-based model position */
                  p_lane[r_lane] = 0.0f;
                  if (k_local >= M) continue;      /* beyond model end       */

                  float emit_c0 = p7_fs_oprofile_FGetEmission(om_fs, k_local, C0);

                  /* acc0: C0-type (0 nt before donor) */
                  if (acc0 == ACCEPT_AG)
                    {
                      float sg = (P_scores[k_local][p7S_GTAG] > -eslINFINITY)
                                   ? expf(P_scores[k_local][p7S_GTAG] + signal_scores[p7S_GTAG]) : 0.0f;
                      float sc = (P_scores[k_local][p7S_GCAG] > -eslINFINITY)
                                   ? expf(P_scores[k_local][p7S_GCAG] + signal_scores[p7S_GCAG]) : 0.0f;
                      p_lane[r_lane] = ESL_MAX(p_lane[r_lane], ESL_MAX(sg, sc) * emit_c0);
                    }
                  else if (acc0 == ACCEPT_AC)
                    {
                      float sa = (P_scores[k_local][p7S_ATAC] > -eslINFINITY)
                                   ? expf(P_scores[k_local][p7S_ATAC] + signal_scores[p7S_ATAC]) : 0.0f;
                      p_lane[r_lane] = ESL_MAX(p_lane[r_lane], sa * emit_c0);
                    }

                  /* acc1: C1-type (1 nt before donor) */
                  if (acc1 == ACCEPT_AG)
                    {
                      for (nuc1 = 0; nuc1 <= SP_MAX_NUC; nuc1++)
                        {
                          int off = SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS;
                          float sg = (P_scores[k_local][off + p7S_GTAG] > -eslINFINITY)
                                       ? expf(P_scores[k_local][off + p7S_GTAG] + signal_scores[p7S_GTAG]) : 0.0f;
                          float sc = (P_scores[k_local][off + p7S_GCAG] > -eslINFINITY)
                                       ? expf(P_scores[k_local][off + p7S_GCAG] + signal_scores[p7S_GCAG]) : 0.0f;
                          float ec1 = p7_fs_oprofile_FGetEmission(om_fs, k_local, C1[nuc1]);
                          p_lane[r_lane] = ESL_MAX(p_lane[r_lane], ESL_MAX(sg, sc) * ec1);
                        }
                    }
                  else if (acc1 == ACCEPT_AC)
                    {
                      for (nuc1 = 0; nuc1 <= SP_MAX_NUC; nuc1++)
                        {
                          int off = SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS;
                          float sa = (P_scores[k_local][off + p7S_ATAC] > -eslINFINITY)
                                       ? expf(P_scores[k_local][off + p7S_ATAC] + signal_scores[p7S_ATAC]) : 0.0f;
                          float ec1 = p7_fs_oprofile_FGetEmission(om_fs, k_local, C1[nuc1]);
                          p_lane[r_lane] = ESL_MAX(p_lane[r_lane], sa * ec1);
                        }
                    }

                  /* acc2: C2-type (2 nt before donor, nuc3 = x) */
                  if (acc2 == ACCEPT_AG)
                    {
                      int off = SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS;
                      float sg = (P_scores[k_local][off + p7S_GTAG] > -eslINFINITY)
                                   ? expf(P_scores[k_local][off + p7S_GTAG] + signal_scores[p7S_GTAG]) : 0.0f;
                      float sc = (P_scores[k_local][off + p7S_GCAG] > -eslINFINITY)
                                   ? expf(P_scores[k_local][off + p7S_GCAG] + signal_scores[p7S_GCAG]) : 0.0f;
                      p_lane[r_lane] = ESL_MAX(p_lane[r_lane], ESL_MAX(sg, sc));
                    }
                  else if (acc2 == ACCEPT_AC)
                    {
                      int off = SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS;
                      float sa = (P_scores[k_local][off + p7S_ATAC] > -eslINFINITY)
                                   ? expf(P_scores[k_local][off + p7S_ATAC] + signal_scores[p7S_ATAC]) : 0.0f;
                      p_lane[r_lane] = ESL_MAX(p_lane[r_lane], sa);
                    }
                } /* end r_lane */

              PMO_SP(dpc, q) = _mm_set_ps(p_lane[3], p_lane[2], p_lane[1], p_lane[0]);
              xEv = _mm_max_ps(xEv, PMO_SP(dpc, q));
            }
        } /* end q loop */

      dcv = esl_sse_rightshiftz_float(dcv);
      run_dd_passes(dpc, dcv, om_fs->tfv + 7*Q, Q, M);

      for (q = 0; q < Q; q++)
        xEv = _mm_max_ps(xEv, DMO_SP(dpc, q));
      esl_sse_hmax_ps(xEv, &xE);

      if (xE > 1.0e4f)
        {
          float  sfv_f   = 1.0f / xE;
          float  log_adj = -logf(xE);
          __m128 sfv     = _mm_set1_ps(sfv_f);
          for (q = 0; q < Q; q++)
            {
              MMO_SP(dpc, q) = _mm_mul_ps(MMO_SP(dpc, q), sfv);
              DMO_SP(dpc, q) = _mm_mul_ps(DMO_SP(dpc, q), sfv);
              IMO_SP(dpc, q) = _mm_mul_ps(IMO_SP(dpc, q), sfv);
              PMO_SP(dpc, q) = _mm_mul_ps(PMO_SP(dpc, q), sfv);
            }
          /* Keep P_scores in sync with current scale (mirrors IVX rescaling in
           * viterbi_fs.c).  P_scores stores log-space values; dividing the
           * probability scale by xE subtracts logf(xE) from all log entries.  */
          for (k_local = 1; k_local < M; k_local++)
            for (q = 0; q < SIGNAL_MEM_SIZE; q++)
              if (P_scores[k_local][q] > -eslINFINITY)
                P_scores[k_local][q] += log_adj;
          ox->xmx[i*p7X_NXCELLS + p7X_SCALE] = xE;
          ox->totscale += logf(xE);
          xE = 1.0f;
        }
      else ox->xmx[i*p7X_NXCELLS + p7X_SCALE] = 1.0f;

      ox->xmx[i*p7X_NXCELLS + p7X_E] = xE;

      /* ------------------------------------------------------------------
       * Donor-site k-loop: update P_scores for future P-state reads.
       * Reads from dpd = dpf[i - min_intron - 3].
       * TMP_SC_prob = max(M(dpd,k-1), D(dpd,k-1)) — probability-space.
       * P_scores stores log-space; convert with logf() at write time.
       * ------------------------------------------------------------------ */
      if (don2 >= 0)
        {
          /* C2-type: 2 nt before donor; emission baked in. */
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

                /* k-1 position in striped layout */
                int km1 = k_local - 1;
                int q2  = (km1 - 1) % Q;
                int r2  = (km1 - 1) / Q;

                u_vec.v = MMO_SP(dpd, q2); TMP_SC_prob  = u_vec.p[r2];
                u_vec.v = DMO_SP(dpd, q2); TMP_SC_prob  = ESL_MAX(TMP_SC_prob, u_vec.p[r2]);

                if (TMP_SC_prob <= 0.0f) continue;
                TMP_SC_log = logf(TMP_SC_prob);

                for (nuc3 = 0; nuc3 <= SP_MAX_NUC; nuc3++)
                  {
                    float emit = p7_fs_oprofile_FGetEmission(om_fs, k_local, C2[nuc3]);
                    float log_emit = (emit > 0.0f) ? logf(emit) : -eslINFINITY;
                    int   off      = SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + don_sig;
                    float new_sc   = TMP_SC_log + log_emit;
                    if (new_sc > P_scores[k_local][off])
                      P_scores[k_local][off] = new_sc;
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
                TMP_SC_log = logf(TMP_SC_prob);

                int off = SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + don_sig;
                if (TMP_SC_log > P_scores[k_local][off])
                  P_scores[k_local][off] = TMP_SC_log;
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
                TMP_SC_log = logf(TMP_SC_prob);

                if (TMP_SC_log > P_scores[k_local][don_sig])
                  P_scores[k_local][don_sig] = TMP_SC_log;
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

  return eslOK;
}
/*------------------ end, p7_Viterbi_SplicedGlobal() ---------------*/



/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/




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
                               pli->splice_scores->P_scores,
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

      if (fabs(gsc - osc) > 0.001f)
        esl_fatal("%s: scalar C=%.4f  SSE C=%.4f  diff=%.4f nats",
                  msg, gsc, osc, fabs(gsc - osc));
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
