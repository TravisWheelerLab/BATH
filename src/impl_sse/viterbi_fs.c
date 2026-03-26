/* SSE-accelerated frameshift-aware Viterbi algorithm; full matrix.
 *
 * Probability-space DP: _mm_mul_ps along a path, _mm_max_ps to select
 * the best path at each branch point. Sparse rescaling is used to
 * achieve full dynamic range on long sequences. The final score is:
 *
 *   ret_sc = ox->totscale + logf(max_xC * xf[C][MOVE])
 *
 * where max_xC = max(xC(L), xC(L-1)*LOOP, xC(L-2)*LOOP), accounting
 * for the three possible final codon lengths (3, 4, or 5 nt) that can
 * end exactly at position L.
 *
 * Based on p7_Forward_Frameshift() in fwdback_fs.c. Every sum-of-paths
 * (_mm_add_ps combining alternative paths) is replaced by a
 * max-over-paths (_mm_max_ps); multiplications along a single path
 * (_mm_mul_ps) are unchanged.
 *
 * Contents:
 *   1. p7_Viterbi_Frameshift() implementation.
 *   2. p7_Viterbi_Frameshift_Trace() implementation.
 *   3. Benchmark driver.
 *   4. Unit tests.
 *   5. Test driver.
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
#include "impl_sse.h"

/* IVX intermediate matrix access: [slot 0..p7P_5CODONS-1][stripe 0..Q-1] */
#define IVX(slot, q) (ivxf[(slot)*Q + (q)])

/*****************************************************************
 * 1. p7_Viterbi_Frameshift() implementation
 *****************************************************************/

/* Function:  p7_Viterbi_Frameshift()
 * Synopsis:  SSE-accelerated frameshift-aware Viterbi algorithm, 5 codon lengths, full matrix.
 *
 * Purpose:   Full O(ML) Viterbi algorithm for frameshift-aware translated
 *            comparison between a nucleotide sequence and a frameshift-aware
 *            codon HMM, using five codon lengths (1-5 nucleotides), with
 *            SIMD parallelism and sparse rescaling.
 *
 *            Finds the maximum-probability alignment of nucleotide sequence
 *            <dsq> of length <L> to frameshift-aware codon profile <om_fs>,
 *            storing the full DP matrix in <ox> to enable subsequent
 *            traceback.  The Viterbi score (in nats) is returned in <opt_sc>.
 *
 *            The DP is computed in probability space (scaled odds ratios),
 *            mirroring p7_Forward_Frameshift(), with all sum-over-paths
 *            operations replaced by max-over-paths operations.
 *
 *            <ox> must be created with p7_omx_Create_dpf(M, L, L, p7X_NSCELLS).
 *            Scale factors at each row are written to
 *            <ox->xmx[i*p7X_NXCELLS+p7X_SCALE]>.
 *
 * Args:      dsq    - nucleotide sequence, 1..L
 *            L      - length of dsq
 *            om_fs  - optimized frameshift profile (codon_lengths == 5)
 *            ox     - DP matrix, p7_omx_Create_dpf(M, L, L, p7X_NSCELLS)
 *            opt_sc - optRETURN: Viterbi score in nats
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on bad inputs; <eslERANGE> on score underflow.
 */
int
p7_Viterbi_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_OMX *ox, float *opt_sc)
{
  register __m128 mpv1, dpv1, ipv1;     /* right-shifted prev1 MDI for BM/MM/IM/DM        */
  register __m128 sv;                    /* temporary IVX accumulator                       */
  register __m128 msv;                   /* M-state best-path value for current position    */
  register __m128 dcv;                   /* delayed D(i,q+1) carry                          */
  register __m128 xEv;                   /* E-state partial max (horizontal reduce later)   */
  register __m128 xBv1;                  /* splatted B(i-1)                                 */
  __m128   zerov;                        /* splatted 0.0                                    */
  float    xN, xE, xB, xC, xJ;          /* special state scalars                           */
  float    xN_buf[PARSER_ROWS_FWD];
  float    xB_buf[PARSER_ROWS_FWD];
  float    xJ_buf[PARSER_ROWS_FWD];
  float    xC_buf[PARSER_ROWS_FWD];
  int      b, b1, b3;                    /* circular buffer slots: i, i-1, i-3             */
  int      ivx_1, ivx_2, ivx_3, ivx_4, ivx_5;
  __m128  *dpc, *dpp1, *dpp3;
  __m128  *tp;
  __m128  *ivxf     = NULL;
  __m128  *ivxf_mem = NULL;
  int      Q = p7O_NQF(om_fs->M);
  int      i, q, j, r;
  int      c1, c2, c3, c4, c5;
  int      t, u, v, w, x;
  float    insert_adj;                   /* scale correction for I(i,k) from dpf[i-3]      */
  int      status;

  if (om_fs->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");

  ox->M              = om_fs->M;
  ox->L              = L;
  ox->has_own_scales = TRUE;
  ox->totscale       = 0.0;
  zerov = _mm_setzero_ps();

  /* Allocate IVX: p7P_5CODONS circular rows x Q stripes */
  ESL_ALLOC(ivxf_mem, sizeof(__m128) * p7P_5CODONS * Q + 15);
  ivxf = (__m128 *) (((unsigned long int) ivxf_mem + 15) & (~0xf));

  /* Zero-initialize rows 0..L of the full DP matrix */
  for (r = 0; r <= L; r++)
    for (q = 0; q < Q; q++)
      MMO(ox->dpf[r],q) = DMO(ox->dpf[r],q) = IMO(ox->dpf[r],q) = zerov;

  /* Zero-initialize all IVX rows */
  for (r = 0; r < p7P_5CODONS; r++)
    for (q = 0; q < Q; q++)
      IVX(r, q) = zerov;

  /* Initialize special-state circular buffers.
   * N(0)=N(1)=N(2)=1; B(0)=B(1)=B(2)=T_NM; E=J=C=0. */
  for (r = 0; r < PARSER_ROWS_FWD; r++)
    xN_buf[r] = xB_buf[r] = xJ_buf[r] = xC_buf[r] = 0.0f;
  xN_buf[0] = xN_buf[1] = xN_buf[2] = 1.0f;
  xB_buf[0] = xB_buf[1] = xB_buf[2] = om_fs->xf[p7O_N][p7O_MOVE];

  /* Write rows 0, 1, 2 specials to ox->xmx */
  for (r = 0; r < 3; r++)
    {
      ox->xmx[r*p7X_NXCELLS+p7X_SCALE] = 1.0f;
      ox->xmx[r*p7X_NXCELLS+p7X_E]     = 0.0f;
      ox->xmx[r*p7X_NXCELLS+p7X_N]     = 1.0f;
      ox->xmx[r*p7X_NXCELLS+p7X_J]     = 0.0f;
      ox->xmx[r*p7X_NXCELLS+p7X_B]     = om_fs->xf[p7O_N][p7O_MOVE];
      ox->xmx[r*p7X_NXCELLS+p7X_C]     = 0.0f;
    }

  /* Initialize nucleotide rolling window */
  t = u = v = w = p7P_MAXCODONS5;
  if (dsq[1] < p7P_MAXNUC) x = dsq[1]; else x = p7P_MAXCODONS5;

  /*----------------------------------------------------------------
   * Initialization: i=1 (only 1-nt codon c1 available)
   *----------------------------------------------------------------*/
  i = 1;
  c1 = p7P_CODON1_FS5(x); c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);

  ivx_1 = 1;  /* i % p7P_5CODONS */

  dpc  = ox->dpf[1];
  dpp1 = ox->dpf[0];  /* all zeros */

  xBv1 = _mm_set1_ps(xB_buf[0]);
  tp   = om_fs->tfv;
  dcv  = zerov;
  xEv  = zerov;
  mpv1 = dpv1 = ipv1 = zerov;

  for (q = 0; q < Q; q++)
    {
      /* IVX(1,q) = max(B(0)*BM, 0, 0, 0) = B(0)*BM */
      sv  =                _mm_mul_ps(xBv1, *tp); tp++;   /* BM */
      sv  = _mm_max_ps(sv, _mm_mul_ps(mpv1, *tp)); tp++;  /* MM=0 */
      sv  = _mm_max_ps(sv, _mm_mul_ps(ipv1, *tp)); tp++;  /* IM=0 */
      sv  = _mm_max_ps(sv, _mm_mul_ps(dpv1, *tp)); tp++;  /* DM=0 */
      IVX(ivx_1, q) = sv;

      /* M_C0(1,q): 1-nt codon only */
      msv = _mm_mul_ps(sv, om_fs->rfv[c1][q]);
      MMO(dpc, q) = msv;
      xEv = _mm_max_ps(xEv, msv);

      DMO(dpc, q) = dcv;
      dcv = _mm_mul_ps(msv, *tp); tp++;   /* MD */

      /* I(1,q) = 0 (dpp3 all zeros) */
      IMO(dpc, q) = zerov;
      tp += 2;  /* skip MI, II */
    }

  /* DD paths */
  dcv        = esl_sse_rightshiftz_float(dcv);
  DMO(dpc,0) = zerov;
  tp         = om_fs->tfv + 7*Q;
  for (q = 0; q < Q; q++)
    {
      DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
      dcv         = _mm_mul_ps(DMO(dpc, q), *tp); tp++;
    }
  if (om_fs->M < 100)
    {
      for (j = 1; j < 4; j++)
        {
          dcv = esl_sse_rightshiftz_float(dcv);
          tp  = om_fs->tfv + 7*Q;
          for (q = 0; q < Q; q++)
            {
              DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
              dcv         = _mm_mul_ps(dcv, *tp); tp++;
            }
        }
    }
  else
    {
      for (j = 1; j < 4; j++)
        {
          register __m128 cv;
          dcv = esl_sse_rightshiftz_float(dcv);
          tp  = om_fs->tfv + 7*Q;
          cv  = zerov;
          for (q = 0; q < Q; q++)
            {
              sv          = _mm_max_ps(dcv, DMO(dpc, q));
              cv          = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO(dpc, q)));
              DMO(dpc, q) = sv;
              dcv         = _mm_mul_ps(dcv, *tp); tp++;
            }
          if (! _mm_movemask_ps(cv)) break;
        }
    }
  /* Add D to xEv, then reduce to scalar xE */
  for (q = 0; q < Q; q++)
    xEv = _mm_max_ps(xEv, DMO(dpc, q));
  esl_sse_hmax_ps(xEv, &xE);

  xN = 1.0f;
  xJ = xE * om_fs->xf[p7O_E][p7O_LOOP];
  xC = xE * om_fs->xf[p7O_E][p7O_MOVE];
  xB = ESL_MAX(xN * om_fs->xf[p7O_N][p7O_MOVE], xJ * om_fs->xf[p7O_J][p7O_MOVE]);

  /* Sparse rescaling at i=1 */
  if (xE > 1.0e4f)
    {
      float scale_factor = 1.0f / xE;
      xN *= scale_factor; xJ *= scale_factor; xC *= scale_factor; xB *= scale_factor;
      xEv = _mm_set1_ps(scale_factor);
      for (q = 0; q < Q; q++)
        {
          MMO(dpc,q) = _mm_mul_ps(MMO(dpc,q), xEv);
          DMO(dpc,q) = _mm_mul_ps(DMO(dpc,q), xEv);
          IMO(dpc,q) = _mm_mul_ps(IMO(dpc,q), xEv);
        }
      for (r = 0; r < p7P_5CODONS; r++)
        for (q = 0; q < Q; q++)
          IVX(r, q) = _mm_mul_ps(IVX(r, q), xEv);
      for (r = 0; r < PARSER_ROWS_FWD; r++)
        {
          xN_buf[r] *= scale_factor; xB_buf[r] *= scale_factor;
          xJ_buf[r] *= scale_factor; xC_buf[r] *= scale_factor;
        }
      ox->xmx[1*p7X_NXCELLS+p7X_SCALE] = xE;
      ox->totscale += log(xE);
      xE = 1.0f;
    }
  else ox->xmx[1*p7X_NXCELLS+p7X_SCALE] = 1.0f;

  xN_buf[1] = xN; xB_buf[1] = xB; xJ_buf[1] = xJ; xC_buf[1] = xC;
  ox->xmx[1*p7X_NXCELLS+p7X_E] = xE;
  ox->xmx[1*p7X_NXCELLS+p7X_N] = xN;
  ox->xmx[1*p7X_NXCELLS+p7X_J] = xJ;
  ox->xmx[1*p7X_NXCELLS+p7X_B] = xB;
  ox->xmx[1*p7X_NXCELLS+p7X_C] = xC;

  /*----------------------------------------------------------------
   * Initialization: i=2 (1-nt and 2-nt codons available)
   *----------------------------------------------------------------*/
  i = 2;
  t = u = v = p7P_MAXCODONS5;
  w = x;
  if (dsq[2] < p7P_MAXNUC) x = dsq[2]; else x = p7P_MAXCODONS5;

  c1 = p7P_CODON1_FS5(x);    c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
  c2 = p7P_CODON2_FS5(w, x); c2 = p7P_MINIDX(c2, p7P_DEGEN5_QC1);

  ivx_1 = 2;  /* i % p7P_5CODONS */
  ivx_2 = 1;  /* (i-1) % p7P_5CODONS */

  dpc  = ox->dpf[2];
  dpp1 = ox->dpf[1];

  xBv1 = _mm_set1_ps(xB_buf[1]);
  tp   = om_fs->tfv;
  dcv  = zerov;
  xEv  = zerov;

  mpv1 = esl_sse_rightshiftz_float(MMO(dpp1, Q-1));
  dpv1 = esl_sse_rightshiftz_float(DMO(dpp1, Q-1));
  ipv1 = esl_sse_rightshiftz_float(IMO(dpp1, Q-1));

  for (q = 0; q < Q; q++)
    {
      sv  =                _mm_mul_ps(xBv1, *tp); tp++;
      sv  = _mm_max_ps(sv, _mm_mul_ps(mpv1, *tp)); tp++;
      sv  = _mm_max_ps(sv, _mm_mul_ps(ipv1, *tp)); tp++;
      sv  = _mm_max_ps(sv, _mm_mul_ps(dpv1, *tp)); tp++;
      IVX(ivx_1, q) = sv;

      __m128 mc1 = _mm_mul_ps(sv,             om_fs->rfv[c1][q]);
      __m128 mc2 = _mm_mul_ps(IVX(ivx_2, q), om_fs->rfv[c2][q]);
      msv = _mm_max_ps(mc1, mc2);
      xEv = _mm_max_ps(xEv, msv);

      mpv1 = MMO(dpp1, q);
      dpv1 = DMO(dpp1, q);
      ipv1 = IMO(dpp1, q);

      MMO(dpc, q) = msv;
      DMO(dpc, q) = dcv;

      dcv = _mm_mul_ps(msv, *tp); tp++;   /* MD */

      /* I(2,q) = 0 (no dpp3 yet) */
      IMO(dpc, q) = zerov;
      tp += 2;  /* skip MI, II */
    }

  /* DD paths */
  dcv        = esl_sse_rightshiftz_float(dcv);
  DMO(dpc,0) = zerov;
  tp         = om_fs->tfv + 7*Q;
  for (q = 0; q < Q; q++)
    {
      DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
      dcv         = _mm_mul_ps(DMO(dpc, q), *tp); tp++;
    }
  if (om_fs->M < 100)
    {
      for (j = 1; j < 4; j++)
        {
          dcv = esl_sse_rightshiftz_float(dcv);
          tp  = om_fs->tfv + 7*Q;
          for (q = 0; q < Q; q++)
            {
              DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
              dcv         = _mm_mul_ps(dcv, *tp); tp++;
            }
        }
    }
  else
    {
      for (j = 1; j < 4; j++)
        {
          register __m128 cv;
          dcv = esl_sse_rightshiftz_float(dcv);
          tp  = om_fs->tfv + 7*Q;
          cv  = zerov;
          for (q = 0; q < Q; q++)
            {
              sv          = _mm_max_ps(dcv, DMO(dpc, q));
              cv          = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO(dpc, q)));
              DMO(dpc, q) = sv;
              dcv         = _mm_mul_ps(dcv, *tp); tp++;
            }
          if (! _mm_movemask_ps(cv)) break;
        }
    }
  for (q = 0; q < Q; q++)
    xEv = _mm_max_ps(xEv, DMO(dpc, q));
  esl_sse_hmax_ps(xEv, &xE);

  xN = 1.0f;
  xJ = xE * om_fs->xf[p7O_E][p7O_LOOP];
  xC = xE * om_fs->xf[p7O_E][p7O_MOVE];
  xB = ESL_MAX(xN * om_fs->xf[p7O_N][p7O_MOVE], xJ * om_fs->xf[p7O_J][p7O_MOVE]);

  if (xE > 1.0e4f)
    {
      float scale_factor = 1.0f / xE;
      xN *= scale_factor; xJ *= scale_factor; xC *= scale_factor; xB *= scale_factor;
      xEv = _mm_set1_ps(scale_factor);
      for (q = 0; q < Q; q++)
        {
          MMO(dpc,q) = _mm_mul_ps(MMO(dpc,q), xEv);
          DMO(dpc,q) = _mm_mul_ps(DMO(dpc,q), xEv);
          IMO(dpc,q) = _mm_mul_ps(IMO(dpc,q), xEv);
        }
      /* Rows 0 and 1 are committed; do not retroactively scale them.
       * The insert_adj correction handles the scale difference when they
       * are read as dpp3 at i=3 and i=4 respectively. */
      for (r = 0; r < p7P_5CODONS; r++)
        for (q = 0; q < Q; q++)
          IVX(r, q) = _mm_mul_ps(IVX(r, q), xEv);
      for (r = 0; r < PARSER_ROWS_FWD; r++)
        {
          xN_buf[r] *= scale_factor; xB_buf[r] *= scale_factor;
          xJ_buf[r] *= scale_factor; xC_buf[r] *= scale_factor;
        }
      ox->xmx[2*p7X_NXCELLS+p7X_SCALE] = xE;
      ox->totscale += log(xE);
      xE = 1.0f;
    }
  else ox->xmx[2*p7X_NXCELLS+p7X_SCALE] = 1.0f;

  xN_buf[2] = xN; xB_buf[2] = xB; xJ_buf[2] = xJ; xC_buf[2] = xC;
  ox->xmx[2*p7X_NXCELLS+p7X_E] = xE;
  ox->xmx[2*p7X_NXCELLS+p7X_N] = xN;
  ox->xmx[2*p7X_NXCELLS+p7X_J] = xJ;
  ox->xmx[2*p7X_NXCELLS+p7X_B] = xB;
  ox->xmx[2*p7X_NXCELLS+p7X_C] = xC;

  /*----------------------------------------------------------------
   * Main recurrence: i = 3..L
   *
   * M_Cn(i,k) = IVX(i-n+1, k) * Rn(k)              for n = 1..5
   * M_C0(i,k) = max(M_C1, M_C2, M_C3, M_C4, M_C5)  (best codon length)
   * I(i,k)    = max(M(i-3,k)*MI, I(i-3,k)*II)
   * D(i,k)    = max(M(i,k-1)*MD, D(i,k-1)*DD)
   *
   * insert_adj: I(i,k) reads dpf[i-3] which was committed under the
   * cumulative scale at row i-3; running scale is that of row i-1.
   * Correction = 1 / (S[i-2] * S[i-1]).
   *----------------------------------------------------------------*/
  for (i = 3; i <= L; i++)
    {
      t = u; u = v; v = w; w = x;
      if (dsq[i] < p7P_MAXNUC) x = dsq[i]; else x = p7P_MAXCODONS5;

      c1 = p7P_CODON1_FS5(x);             c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
      c2 = p7P_CODON2_FS5(w, x);          c2 = p7P_MINIDX(c2, p7P_DEGEN5_QC1);
      c3 = p7P_CODON3_FS5(v, w, x);       c3 = p7P_MINIDX(c3, p7P_DEGEN5_C);
      c4 = p7P_CODON4_FS5(u, v, w, x);    c4 = p7P_MINIDX(c4, p7P_DEGEN5_QC1);
      c5 = p7P_CODON5_FS5(t, u, v, w, x); c5 = p7P_MINIDX(c5, p7P_DEGEN5_QC2);

      ivx_1 =     i               % p7P_5CODONS;
      ivx_2 = ((i-1) % p7P_5CODONS + p7P_5CODONS) % p7P_5CODONS;
      ivx_3 = ((i-2) % p7P_5CODONS + p7P_5CODONS) % p7P_5CODONS;
      ivx_4 = ((i-3) % p7P_5CODONS + p7P_5CODONS) % p7P_5CODONS;
      ivx_5 = ((i-4) % p7P_5CODONS + p7P_5CODONS) % p7P_5CODONS;

      b  =     i               % PARSER_ROWS_FWD;
      b1 = ((i-1) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD;
      b3 = ((i-3) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD;

      dpc  = ox->dpf[i];
      dpp1 = ox->dpf[i-1];
      dpp3 = ox->dpf[i-3];

      /* insert_adj: dpp3 committed at scale S[i-3]; running scale = S[i-1].
       * Correction = 1.0 / (S[i-2] * S[i-1]). */
      insert_adj = 1.0f
                 / (ox->xmx[(i-2)*p7X_NXCELLS+p7X_SCALE]
                 *  ox->xmx[(i-1)*p7X_NXCELLS+p7X_SCALE]);

      mpv1 = esl_sse_rightshiftz_float(MMO(dpp1, Q-1));
      dpv1 = esl_sse_rightshiftz_float(DMO(dpp1, Q-1));
      ipv1 = esl_sse_rightshiftz_float(IMO(dpp1, Q-1));

      xBv1 = _mm_set1_ps(xB_buf[b1]);
      tp   = om_fs->tfv;
      dcv  = zerov;
      xEv  = zerov;

      register __m128 adj_v = _mm_set1_ps(insert_adj);

      for (q = 0; q < Q; q++)
        {
          /* IVX(i,k) = max(B(i-1)*BM, M(i-1,k-1)*MM, I(i-1,k-1)*IM, D(i-1,k-1)*DM) */
          sv  =                _mm_mul_ps(xBv1, *tp); tp++;
          sv  = _mm_max_ps(sv, _mm_mul_ps(mpv1, *tp)); tp++;
          sv  = _mm_max_ps(sv, _mm_mul_ps(ipv1, *tp)); tp++;
          sv  = _mm_max_ps(sv, _mm_mul_ps(dpv1, *tp)); tp++;
          IVX(ivx_1, q) = sv;

          /* M_C0 = best codon-length match at this model position */
          msv = _mm_max_ps(_mm_max_ps(_mm_max_ps(_mm_mul_ps(sv,             om_fs->rfv[c1][q]),
                                                  _mm_mul_ps(IVX(ivx_2, q), om_fs->rfv[c2][q])),
                                      _mm_max_ps(_mm_mul_ps(IVX(ivx_3, q), om_fs->rfv[c3][q]),
                                                  _mm_mul_ps(IVX(ivx_4, q), om_fs->rfv[c4][q]))),
                                      _mm_mul_ps(IVX(ivx_5, q),             om_fs->rfv[c5][q]));
          xEv = _mm_max_ps(xEv, msv);

          mpv1 = MMO(dpp1, q);
          dpv1 = DMO(dpp1, q);
          ipv1 = IMO(dpp1, q);

          MMO(dpc, q) = msv;
          DMO(dpc, q) = dcv;
          dcv = _mm_mul_ps(msv, *tp); tp++;   /* MD */

          /* I(i,k) = max(M(i-3,k)*MI, I(i-3,k)*II); scale-corrected */
          sv          = _mm_mul_ps(_mm_mul_ps(MMO(dpp3, q), adj_v), *tp); tp++;  /* MI */
          IMO(dpc, q) = _mm_max_ps(sv, _mm_mul_ps(_mm_mul_ps(IMO(dpp3, q), adj_v), *tp)); tp++;  /* II */
        }

      /* DD paths */
      dcv        = esl_sse_rightshiftz_float(dcv);
      DMO(dpc,0) = zerov;
      tp         = om_fs->tfv + 7*Q;
      for (q = 0; q < Q; q++)
        {
          DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
          dcv         = _mm_mul_ps(DMO(dpc, q), *tp); tp++;
        }
      if (om_fs->M < 100)
        {
          for (j = 1; j < 4; j++)
            {
              dcv = esl_sse_rightshiftz_float(dcv);
              tp  = om_fs->tfv + 7*Q;
              for (q = 0; q < Q; q++)
                {
                  DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
                  dcv         = _mm_mul_ps(dcv, *tp); tp++;
                }
            }
        }
      else
        {
          for (j = 1; j < 4; j++)
            {
              register __m128 cv;
              dcv = esl_sse_rightshiftz_float(dcv);
              tp  = om_fs->tfv + 7*Q;
              cv  = zerov;
              for (q = 0; q < Q; q++)
                {
                  sv          = _mm_max_ps(dcv, DMO(dpc, q));
                  cv          = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO(dpc, q)));
                  DMO(dpc, q) = sv;
                  dcv         = _mm_mul_ps(dcv, *tp); tp++;
                }
              if (! _mm_movemask_ps(cv)) break;
            }
        }

      for (q = 0; q < Q; q++)
        xEv = _mm_max_ps(xEv, DMO(dpc, q));
      esl_sse_hmax_ps(xEv, &xE);

      xN = xN_buf[b3] * om_fs->xf[p7O_N][p7O_LOOP];
      xJ = ESL_MAX(xJ_buf[b3] * om_fs->xf[p7O_J][p7O_LOOP], xE * om_fs->xf[p7O_E][p7O_LOOP]);
      xC = ESL_MAX(xC_buf[b3] * om_fs->xf[p7O_C][p7O_LOOP], xE * om_fs->xf[p7O_E][p7O_MOVE]);
      xB = ESL_MAX(xN         * om_fs->xf[p7O_N][p7O_MOVE],  xJ * om_fs->xf[p7O_J][p7O_MOVE]);

      /* Sparse rescaling */
      if (xE > 1.0e4f)
        {
          float scale_factor = 1.0f / xE;
          xN *= scale_factor; xJ *= scale_factor; xC *= scale_factor; xB *= scale_factor;
          xEv = _mm_set1_ps(scale_factor);
          for (q = 0; q < Q; q++)
            {
              MMO(dpc,q) = _mm_mul_ps(MMO(dpc,q), xEv);
              DMO(dpc,q) = _mm_mul_ps(DMO(dpc,q), xEv);
              IMO(dpc,q) = _mm_mul_ps(IMO(dpc,q), xEv);
            }
          for (r = 0; r < p7P_5CODONS; r++)
            for (q = 0; q < Q; q++)
              IVX(r, q) = _mm_mul_ps(IVX(r, q), xEv);
          for (r = 0; r < PARSER_ROWS_FWD; r++)
            {
              xN_buf[r] *= scale_factor; xB_buf[r] *= scale_factor;
              xJ_buf[r] *= scale_factor; xC_buf[r] *= scale_factor;
            }
          ox->xmx[i*p7X_NXCELLS+p7X_SCALE] = xE;
          ox->totscale += log(xE);
          xE = 1.0f;
        }
      else ox->xmx[i*p7X_NXCELLS+p7X_SCALE] = 1.0f;

      xN_buf[b] = xN; xB_buf[b] = xB; xJ_buf[b] = xJ; xC_buf[b] = xC;
      ox->xmx[i*p7X_NXCELLS+p7X_E] = xE;
      ox->xmx[i*p7X_NXCELLS+p7X_N] = xN;
      ox->xmx[i*p7X_NXCELLS+p7X_J] = xJ;
      ox->xmx[i*p7X_NXCELLS+p7X_B] = xB;
      ox->xmx[i*p7X_NXCELLS+p7X_C] = xC;
    } /* end main loop i=3..L */

  /* Final score: max of xC(L), xC(L-1)*LOOP, xC(L-2)*LOOP, then *MOVE.
   * The best-scoring alignment can end with a codon of length 3, 4, or 5 nt. */
  {
    float xCL   = xC_buf[   L    % PARSER_ROWS_FWD];
    float xCLm1 = xC_buf[((L-1) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD];
    float xCLm2 = xC_buf[((L-2) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD];
    float xCtot = ESL_MAX(xCL,
                  ESL_MAX(xCLm1 * om_fs->xf[p7O_C][p7O_LOOP],
                          xCLm2 * om_fs->xf[p7O_C][p7O_LOOP]));

    if      (isnan(xCtot))           ESL_EXCEPTION(eslERANGE, "Viterbi score is NaN");
    else if (isinf(xCtot) == 1)      ESL_EXCEPTION(eslERANGE, "Viterbi score overflow (is infinity)");
    else if (L > 1 && xCtot == 0.0f) {
      if (opt_sc != NULL) *opt_sc = -eslINFINITY;
      free(ivxf_mem);
      return eslERANGE;
    }

    if (opt_sc != NULL)
      *opt_sc = ox->totscale + logf(xCtot * om_fs->xf[p7O_C][p7O_MOVE]);
  }

  free(ivxf_mem);
  return eslOK;

 ERROR:
  if (ivxf_mem) free(ivxf_mem);
  return status;
}
/*------------------ end, p7_Viterbi_Frameshift() ---------------*/



/*****************************************************************
 * 2. p7_Viterbi_Frameshift_Trace() implementation
 *****************************************************************/

/* Forward declarations for vit_select_*_fs() helpers */
static inline int vit_select_m_fs      (const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k);
static inline int vit_select_d_fs      (const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k);
static inline int vit_select_i_fs      (const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k);
static inline int vit_select_n_fs      (int i);
static inline int vit_select_c_fs      (const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i);
static inline int vit_select_j_fs      (const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i);
static inline int vit_select_e_fs      (const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int *ret_k);
static inline int vit_select_b_fs      (const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i);
static inline int vit_select_codon_len_fs(const ESL_DSQ *dsq, int i, int k, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox);


/* Function:  p7_Viterbi_Framshift_Trace()
 * Synopsis:  Traceback from a Viterbi DP matrix (FS 8-cell layout); SSE version.
 *
 * Purpose:   Extract the optimal (Viterbi) alignment path from the DP matrix
 *            <ox> that was filled by <p7_Viterbi_Frameshift()>.  The trace is
 *            written into <tr> (caller provides an already-allocated trace via
 *            <p7_trace_fs_Create()>).
 *
 * Args:      dsq   - digital nucleotide sequence, 1..L (used only for L)
 *            L     - sequence length
 *            om_fs - optimized frameshift profile
 *            ox    - Viterbi DP matrix from p7_Viterbi_Frameshift()
 *            tr    - RETURN: traceback; caller provides initial alloc
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if an impossible state is reached.
 */
int
p7_Viterbi_Frameshift_Trace(const ESL_DSQ *dsq, int L,
                 const P7_FS_OPROFILE *om_fs, const P7_OMX *ox,
                 P7_TRACE *tr)
{
  int i  = L;
  int k  = 0;
  int c  = 0;
  int s0, s1;
  int status;

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i, c)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i, c)) != eslOK) return status;
  s0 = p7T_C;

  while (s0 != p7T_S)
    {
      switch (s0) {
      case p7T_M: s1 = vit_select_m_fs(om_fs, ox, i, k);  k--;    break;
      case p7T_D: s1 = vit_select_d_fs(om_fs, ox, i, k);  k--;    break;
      case p7T_I: s1 = vit_select_i_fs(om_fs, ox, i, k);  i -= 3; break;
      case p7T_N: s1 = vit_select_n_fs(i);                         break;
      case p7T_C: s1 = vit_select_c_fs(om_fs, ox, i);              break;
      case p7T_J: s1 = vit_select_j_fs(om_fs, ox, i);              break;
      case p7T_E: s1 = vit_select_e_fs(om_fs, ox, i, &k);          break;
      case p7T_B: s1 = vit_select_b_fs(om_fs, ox, i);              break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in Viterbi traceback");
      }
      if (s1 == -1) ESL_EXCEPTION(eslEINVAL, "Viterbi traceback choice failed");

      /* For M states: recompute codon length c (1..5) from stored DP cells. */
      if (s1 == p7T_M) {
        c = vit_select_codon_len_fs(dsq, i, k, om_fs, ox);
        if (i - c < 0) s1 = p7T_B;
      } else {
        c = 0;
      }

      if ((status = p7_trace_fs_Append(tr, s1, k, i, c)) != eslOK) return status;

      /* For NCJ states: deferred i decrement (one nucleotide per loop iteration). */
      if ((s1 == p7T_N || s1 == p7T_C || s1 == p7T_J) && s1 == s0) i--;

      s0 = s1;
      i -= c;
    }

  tr->M = om_fs->M;
  tr->L = L;
  return p7_trace_fs_Reverse(tr);
}
/*------------------ end, p7_Viterbi_Frameshift_Trace() --------------------*/


/*****************************************************************
 * 2a. vit_select_*_fs() helper functions
 *****************************************************************/

/* M(i,k) is reached from B(i), M(i,k-1), I(i,k-1), or D(i,k-1).
 * All predecessors are at the same row i.
 */
static inline int
vit_select_m_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m128 *tp    = om_fs->tfv + 7*q;   /* BM, MM, IM, DM, MD, MI, II */
  __m128  xBv   = _mm_set1_ps(ox->xmx[i*p7X_NXCELLS + p7X_B]);
  __m128  mpv, dpv, ipv;
  union { __m128 v; float p[4]; } u;
  float   path[4];
  int     state[4] = { p7T_B, p7T_M, p7T_I, p7T_D };

  if (q > 0) {
    mpv = MMO(ox->dpf[i], q-1);
    dpv = DMO(ox->dpf[i], q-1);
    ipv = IMO(ox->dpf[i], q-1);
  } else {
    mpv = esl_sse_rightshiftz_float(MMO(ox->dpf[i], Q-1));
    dpv = esl_sse_rightshiftz_float(DMO(ox->dpf[i], Q-1));
    ipv = esl_sse_rightshiftz_float(IMO(ox->dpf[i], Q-1));
  }

  u.v = _mm_mul_ps(xBv, *tp); tp++;  path[0] = u.p[r];  /* B  * T_BM */
  u.v = _mm_mul_ps(mpv, *tp); tp++;  path[1] = u.p[r];  /* M' * T_MM */
  u.v = _mm_mul_ps(ipv, *tp); tp++;  path[2] = u.p[r];  /* I' * T_IM */
  u.v = _mm_mul_ps(dpv, *tp);        path[3] = u.p[r];  /* D' * T_DM */
  return state[esl_vec_FArgMax(path, 4)];
}

/* D(i,k) is reached from M(i,k-1) or D(i,k-1). Same row i. */
static inline int
vit_select_d_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m128  mpv, dpv;
  __m128  tmdv, tddv;
  union { __m128 v; float p[4]; } u;
  float   path[2];
  int     state[2] = { p7T_M, p7T_D };

  if (q > 0) {
    mpv  = MMO(ox->dpf[i], q-1);
    dpv  = DMO(ox->dpf[i], q-1);
    tmdv = om_fs->tfv[7*(q-1) + p7O_MD];
    tddv = om_fs->tfv[7*Q + (q-1)];
  } else {
    mpv  = esl_sse_rightshiftz_float(MMO(ox->dpf[i], Q-1));
    dpv  = esl_sse_rightshiftz_float(DMO(ox->dpf[i], Q-1));
    tmdv = esl_sse_rightshiftz_float(om_fs->tfv[7*(Q-1) + p7O_MD]);
    tddv = esl_sse_rightshiftz_float(om_fs->tfv[8*Q-1]);
  }

  u.v = _mm_mul_ps(mpv, tmdv); path[0] = u.p[r];
  u.v = _mm_mul_ps(dpv, tddv); path[1] = u.p[r];
  return state[esl_vec_FArgMax(path, 2)];
}

/* I(i,k) is reached from M(i-3,k) or I(i-3,k). */
static inline int
vit_select_i_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q    = p7O_NQF(ox->M);
  int     q    = (k-1) % Q;
  int     r    = (k-1) / Q;
  __m128  mpv  = MMO(ox->dpf[i-3], q);
  __m128  ipv  = IMO(ox->dpf[i-3], q);
  __m128 *tp   = om_fs->tfv + 7*q + p7O_MI;
  union { __m128 v; float p[4]; } u;
  float   path[2];
  int     state[2] = { p7T_M, p7T_I };

  u.v = _mm_mul_ps(mpv, *tp); tp++;  path[0] = u.p[r];  /* M * T_MI */
  u.v = _mm_mul_ps(ipv, *tp);        path[1] = u.p[r];  /* I * T_II */
  return state[esl_vec_FArgMax(path, 2)];
}

/* N(i) must come from N(i-1) for i>0, or S for i==0. */
static inline int
vit_select_n_fs(int i)
{
  return (i == 0) ? p7T_S : p7T_N;
}

/* C(i) is reached from E(i) or C at one of three codon offsets (i-3, i-2, i-1).
 * Scale-factor adjustment brings all candidates to a common scale before argmax.
 */
static inline int
vit_select_c_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i)
{
  float loop  = om_fs->xf[p7O_C][p7O_LOOP];
  float s_im2, s_im1, s_i;
  float path[4];

  if (i < 4) return p7T_E;

  s_im2 = ox->xmx[(i-2)*p7X_NXCELLS + p7X_SCALE];
  s_im1 = ox->xmx[(i-1)*p7X_NXCELLS + p7X_SCALE];
  s_i   = ox->xmx[ i   *p7X_NXCELLS + p7X_SCALE];

  path[0] = ox->xmx[(i-3)*p7X_NXCELLS + p7X_C] * loop;
  path[1] = ox->xmx[(i-2)*p7X_NXCELLS + p7X_C] * loop  * s_im2;
  path[2] = ox->xmx[(i-1)*p7X_NXCELLS + p7X_C] * loop  * s_im2 * s_im1;
  path[3] = ox->xmx[ i   *p7X_NXCELLS + p7X_E] * om_fs->xf[p7O_E][p7O_MOVE] * s_im2 * s_im1 * s_i;
  return (esl_vec_FArgMax(path, 4) < 3) ? p7T_C : p7T_E;
}

/* J(i) is reached from E(i) or J at one of three codon offsets (i-3, i-2, i-1).
 * Same scale-factor adjustment logic as vit_select_c_fs().
 */
static inline int
vit_select_j_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i)
{
  float loop  = om_fs->xf[p7O_J][p7O_LOOP];
  float s_im2, s_im1, s_i;
  float path[4];

  if (i < 4) return p7T_E;

  s_im2 = ox->xmx[(i-2)*p7X_NXCELLS + p7X_SCALE];
  s_im1 = ox->xmx[(i-1)*p7X_NXCELLS + p7X_SCALE];
  s_i   = ox->xmx[ i   *p7X_NXCELLS + p7X_SCALE];

  path[0] = ox->xmx[(i-3)*p7X_NXCELLS + p7X_J] * loop;
  path[1] = ox->xmx[(i-2)*p7X_NXCELLS + p7X_J] * loop  * s_im2;
  path[2] = ox->xmx[(i-1)*p7X_NXCELLS + p7X_J] * loop  * s_im2 * s_im1;
  path[3] = ox->xmx[ i   *p7X_NXCELLS + p7X_E] * om_fs->xf[p7O_E][p7O_LOOP] * s_im2 * s_im1 * s_i;
  return (esl_vec_FArgMax(path, 4) < 3) ? p7T_J : p7T_E;
}

/* E(i) is reached from any M(i,k=1..M) or D(i,k=2..M).
 * Scan all (q,r) to find the maximum M or D cell; return p7T_M or p7T_D
 * and set *ret_k to the winning model position.
 */
static inline int
vit_select_e_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int *ret_k)
{
  int    Q    = p7O_NQF(ox->M);
  union { __m128 v; float p[4]; } u;
  float  best = -1.0f;
  int    best_k = 1;
  int    best_state = p7T_M;
  int    q, r, k;

  for (q = 0; q < Q; q++)
    {
      u.v = MMO(ox->dpf[i], q);
      for (r = 0; r < 4; r++) {
        k = r*Q + q + 1;
        if (k <= ox->M && u.p[r] > best) { best = u.p[r]; best_k = k; best_state = p7T_M; }
      }
      u.v = DMO(ox->dpf[i], q);
      for (r = 0; r < 4; r++) {
        k = r*Q + q + 1;
        if (k <= ox->M && u.p[r] > best) { best = u.p[r]; best_k = k; best_state = p7T_D; }
      }
    }
  *ret_k = best_k;
  return best_state;
}

/* B(i) is reached from N(i) or J(i). Same row, no scale adjustment needed. */
static inline int
vit_select_b_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i)
{
  float   path[2];
  int     state[2] = { p7T_N, p7T_J };

  path[0] = ox->xmx[i*p7X_NXCELLS + p7X_N] * om_fs->xf[p7O_N][p7O_MOVE];
  path[1] = ox->xmx[i*p7X_NXCELLS + p7X_J] * om_fs->xf[p7O_J][p7O_MOVE];
  return state[esl_vec_FArgMax(path, 2)];
}

/* For M state at (i,k): recompute codon length c (1..5) by reconstructing
 * per-codon IVX scores from stored DP cells, with row-scale corrections.
 *
 * For codon length n, the predecessor DP row is i-n (at cumulative scale
 * S[i-n]).  To compare across n values at a common scale (S[i-1]), multiply
 * each raw score by corr = product(SCALE[r] for r = i-n+1 .. i-1).
 * Built cumulatively: after using corr for n, update corr *= SCALE[i-n].
 */
static inline int
vit_select_codon_len_fs(const ESL_DSQ *dsq, int i, int k,
                        const P7_FS_OPROFILE *om_fs, const P7_OMX *ox)
{
  int    Q     = p7O_NQF(ox->M);
  int    q_k   = (k-1) % Q;   /* stripe for model pos k   */
  int    r_k   = (k-1) / Q;   /* element for model pos k  */
  int    q_km1 = (k > 1) ? (k-2) % Q : 0;
  int    r_km1 = (k > 1) ? (k-2) / Q : 0;
  union { __m128 v; float p[4]; } u;
  float  tbm, tmm, tim, tdm;
  float  bprev, mprev, iprev, dprev;
  float  raw_ivx, rfv_cn;
  float  scores[5];
  float  corr;
  int    x2, w2, v2, u2, t2;
  int    cn, n;

  /* Scalar transition probs for model position k */
  u.v = om_fs->tfv[7*q_k + 0]; tbm = u.p[r_k];
  u.v = om_fs->tfv[7*q_k + 1]; tmm = u.p[r_k];
  u.v = om_fs->tfv[7*q_k + 2]; tim = u.p[r_k];
  u.v = om_fs->tfv[7*q_k + 3]; tdm = u.p[r_k];

  /* Nucleotide indices for codon lengths 1..5 ending at position i */
  x2 = (i >= 1 && dsq[i]   < p7P_MAXNUC) ? (int)dsq[i]   : p7P_MAXCODONS5;
  w2 = (i >= 2 && dsq[i-1] < p7P_MAXNUC) ? (int)dsq[i-1] : p7P_MAXCODONS5;
  v2 = (i >= 3 && dsq[i-2] < p7P_MAXNUC) ? (int)dsq[i-2] : p7P_MAXCODONS5;
  u2 = (i >= 4 && dsq[i-3] < p7P_MAXNUC) ? (int)dsq[i-3] : p7P_MAXCODONS5;
  t2 = (i >= 5 && dsq[i-4] < p7P_MAXNUC) ? (int)dsq[i-4] : p7P_MAXCODONS5;

  corr = 1.0f;
  for (n = 1; n <= 5; n++)
    {
      if (i < n) { scores[n-1] = 0.0f; continue; }

      /* Predecessor B/M/I/D at row i-n */
      bprev = ox->xmx[(i-n)*p7X_NXCELLS + p7X_B];
      if (k > 1) {
        u.v = MMO(ox->dpf[i-n], q_km1); mprev = u.p[r_km1];
        u.v = IMO(ox->dpf[i-n], q_km1); iprev = u.p[r_km1];
        u.v = DMO(ox->dpf[i-n], q_km1); dprev = u.p[r_km1];
      } else {
        mprev = iprev = dprev = 0.0f;
      }

      raw_ivx = ESL_MAX(bprev*tbm, ESL_MAX(mprev*tmm, ESL_MAX(iprev*tim, dprev*tdm)));

      switch (n) {
        case 1: cn = p7P_CODON1_FS5(x2);                   cn = p7P_MINIDX(cn, p7P_DEGEN5_QC2); break;
        case 2: cn = p7P_CODON2_FS5(w2, x2);               cn = p7P_MINIDX(cn, p7P_DEGEN5_QC1); break;
        case 3: cn = p7P_CODON3_FS5(v2, w2, x2);           cn = p7P_MINIDX(cn, p7P_DEGEN5_C);   break;
        case 4: cn = p7P_CODON4_FS5(u2, v2, w2, x2);       cn = p7P_MINIDX(cn, p7P_DEGEN5_QC1); break;
        default:cn = p7P_CODON5_FS5(t2, u2, v2, w2, x2);  cn = p7P_MINIDX(cn, p7P_DEGEN5_QC2); break;
      }
      u.v = om_fs->rfv[cn][q_k]; rfv_cn = u.p[r_k];

      scores[n-1] = raw_ivx * rfv_cn * corr;

      /* Accumulate scale correction: bring row i-n to scale of row i-1 */
      corr *= ox->xmx[(i-n)*p7X_NXCELLS + p7X_SCALE];
    }

  return esl_vec_FArgMax(scores, 5) + 1;  /* returns 1..5 */
}
/*-------------- end, vit_select_*_fs() helpers -----------------*/


/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/
#ifdef p7VITERBI_FS_BENCHMARK
/*
   gcc -g -O3 -msse2 -std=gnu99 -o viterbi_fs_benchmark -I.. -L.. -I../../easel -L../../easel -Dp7VITERBI_FS_BENCHMARK viterbi_fs.c -lhmmer -leasel -lm
   icc  -O3 -static  -o viterbi_fs_benchmark -I.. -L.. -I../../easel -L../../easel -Dp7VITERBI_FS_BENCHMARK viterbi_fs.c -lhmmer -leasel -lm

   ./viterbi_fs_benchmark <hmmfile>            benchmark Viterbi + traceback together
   ./viterbi_fs_benchmark -V <hmmfile>         benchmark Viterbi only
   ./viterbi_fs_benchmark -T <hmmfile>         benchmark traceback only
   ./viterbi_fs_benchmark -c -N100 <hmmfile>   compare scores to generic Viterbi
*/
#include "p7_config.h"

#include <math.h>
#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name    type           default  env  range  toggles reqs incomp  help                                             docgroup */
  { "-h",  eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL,  "show brief help on version and usage",                0 },
  { "-c",  eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL,  "compare scores to generic implementation (debug)",    0 },
  { "-s",  eslARG_INT,     "42", NULL, NULL,   NULL,  NULL, NULL,  "set random number seed to <n>",                       0 },
  { "-L",  eslARG_INT,   "1200", NULL, "n>0",  NULL,  NULL, NULL,  "length of random target DNA seqs (multiple of 3)",    0 },
  { "-N",  eslARG_INT,   "2000", NULL, "n>0",  NULL,  NULL, NULL,  "number of random target seqs",                        0 },
  { "-V",  eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, "-T",  "benchmark Viterbi only (skip traceback)",             0 },
  { "-T",  eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, "-V",  "benchmark traceback only (skip Viterbi timing)",      0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for SSE frameshift Viterbi";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA   = NULL;
  ESL_ALPHABET   *abcDNA  = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bgAA    = NULL;
  P7_BG          *bgDNA   = NULL;
  ESL_GENCODE    *gcode   = NULL;
  P7_FS_PROFILE  *gm_fs5  = NULL;
  P7_FS_OPROFILE *om_fs5  = NULL;
  P7_OMX         *ox      = NULL;
  P7_GMX         *gx      = NULL;   /* only allocated if -c */
  P7_IVX         *iv      = NULL;   /* only allocated if -c */
  P7_TRACE       *tr      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc, sc2;
  double          base_time, bench_time, Mcs;
  int             do_viterbi   = ! esl_opt_GetBoolean(go, "-T");
  int             do_traceback = ! esl_opt_GetBoolean(go, "-V");

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)          != eslOK) p7_Fail("Failed to read HMM");

  abcDNA = esl_alphabet_Create(eslDNA);
  gcode  = esl_gencode_Create(abcDNA, abcAA);
  bgAA   = p7_bg_Create(abcAA);
  bgDNA  = p7_bg_Create(abcDNA);
  p7_bg_SetLength(bgDNA, L);

  gm_fs5 = p7_profile_fs_Create(hmm->M, abcAA, p7P_5CODONS);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_UNILOCAL);

  om_fs5 = p7_fs_oprofile_Create(hmm->M, abcAA, p7P_5CODONS);
  p7_fs_oprofile_Convert(gm_fs5, om_fs5);
  p7_fs_oprofile_ReconfigLength(om_fs5, L/3);

  ox = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS);
  tr = p7_trace_fs_Create();

  if (esl_opt_GetBoolean(go, "-c")) {
    gx = p7_gmx_Create(gm_fs5->M, L, L, p7G_NSCELLS);
    iv = p7_ivx_Create(gm_fs5->M, p7P_5CODONS);
  }

  /* If benchmarking traceback only, pre-fill the matrix once. */
  if (! do_viterbi) {
    esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
    p7_Viterbi_Frameshift(dsq, L, om_fs5, ox, &sc);
  }

  /* Baseline time: just sequence generation */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark time */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);

      if (do_viterbi)
        p7_Viterbi_Frameshift(dsq, L, om_fs5, ox, &sc);

      if (esl_opt_GetBoolean(go, "-c")) {
        p7_GViterbi_Frameshift(dsq, L, gm_fs5, gx, iv, &sc2);
        printf("%.4f %.4f\n", sc, sc2);
        p7_gmx_Reuse(gx);
      }

      if (do_traceback) {
        p7_Viterbi_Frameshift_Trace(dsq, L, om_fs5, ox, tr);
        p7_trace_Reuse(tr);
      }
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm_fs5->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n", gm_fs5->M);
  printf("# L    = %d\n", L);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_trace_fs_Destroy(tr);
  p7_omx_Destroy(ox);
  if (iv) p7_ivx_Destroy(iv);
  if (gx) p7_gmx_Destroy(gx);
  p7_fs_oprofile_Destroy(om_fs5);
  p7_profile_fs_Destroy(gm_fs5);
  p7_bg_Destroy(bgDNA);
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
#endif /*p7VITERBI_FS_BENCHMARK*/
/*----------------- end, benchmark driver -----------------------*/



/*****************************************************************
 * 4. Unit tests.
 *****************************************************************/
#ifdef p7VITERBI_FS_TESTDRIVE
#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_viterbi_fs()
 *
 * For N random HMMs and randomly reverse-translated DNA sequences:
 *   1. SSE and generic Viterbi scores agree (< 0.001 nat tolerance).
 *   2. Each SSE trace is structurally valid.
 *   3. Each generic trace is structurally valid.
 *
 * Direct path comparison is intentionally omitted: when multiple
 * paths tie for the best score, the two argmax implementations may
 * break ties differently and still both be correct.
 */
static void
utest_viterbi_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
                 ESL_GENCODE *gcode, P7_BG *bgAA,
                 P7_CODONTABLE *codon_table, int M, int N)
{
  char           *msg    = "utest_viterbi_fs failed";
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_fs5 = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE *om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq    = NULL;
  P7_TRACE       *tr     = p7_trace_fs_Create();
  P7_TRACE       *trg    = p7_trace_fs_Create();
  P7_IVX         *iv5    = p7_ivx_Create(M, p7P_5CODONS);
  P7_GMX         *gx     = p7_gmx_Create(M, M, M, p7G_NSCELLS);
  P7_OMX         *ox     = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS);
  char            errbuf[eslERRBUFSIZE];
  float           gsc, osc;
  int             curr_L, i, j;

  if (!gm || !gm_fs5 || !om_fs5 || !sq || !tr || !trg || !iv5 || !gx || !ox) esl_fatal(msg);

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig(hmm, bgAA, gm, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, M, p7_LOCAL);
  p7_fs_oprofile_Convert(gm_fs5, om_fs5);

  while (N--)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, NULL);
      curr_L = sq->n * 3;

      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (curr_L + 2))) == NULL) esl_fatal("malloc failed");
      j = 1;
      for (i = 1; i <= (int)sq->n; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }

      p7_fs_oprofile_ReconfigLength(om_fs5, sq->n);
      p7_fs_ReconfigLength(gm_fs5, sq->n);

      p7_gmx_GrowTo(gx, M, curr_L, curr_L);
      p7_ivx_GrowTo(iv5, M, p7P_5CODONS);
      p7_omx_GrowTo_dpf(ox, M, curr_L, curr_L);

      /* Generic Viterbi */
      if (p7_GViterbi_Frameshift(dsq, curr_L, gm_fs5, gx, iv5, &gsc) != eslOK) esl_fatal(msg);

      /* SSE Viterbi */
      { int s = p7_Viterbi_Frameshift(dsq, curr_L, om_fs5, ox, &osc);
        if (s == eslERANGE) continue;
        if (s != eslOK)     esl_fatal(msg); }

      /* Scores must agree */
      if (fabs(gsc - osc) > 0.001)
        esl_fatal("%s: generic vit %.4f != SSE vit %.4f", msg, gsc, osc);

      /* Traces */
      if (p7_GVTrace_Frameshift(dsq, curr_L, gm_fs5, gx, trg)           != eslOK) esl_fatal(msg);
      if (p7_Viterbi_Frameshift_Trace(dsq, curr_L, om_fs5, ox, tr)       != eslOK) esl_fatal(msg);

      if (p7_trace_fs_Validate(tr,  abcDNA, dsq, errbuf) != eslOK)
        esl_fatal("%s: SSE trace invalid: %s",     msg, errbuf);
      if (p7_trace_fs_Validate(trg, abcDNA, dsq, errbuf) != eslOK)
        esl_fatal("%s: generic trace invalid: %s", msg, errbuf);

      p7_trace_Reuse(tr);
      p7_trace_Reuse(trg);
    }

  free(dsq);
  p7_trace_fs_Destroy(tr);
  p7_trace_fs_Destroy(trg);
  p7_omx_Destroy(ox);
  p7_ivx_Destroy(iv5);
  p7_gmx_Destroy(gx);
  p7_hmm_Destroy(hmm);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_fs5);
  p7_fs_oprofile_Destroy(om_fs5);
}
#endif /*p7VITERBI_FS_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/



/*****************************************************************
 * 5. Test driver.
 *****************************************************************/
#ifdef p7VITERBI_FS_TESTDRIVE
/*
  gcc -g -Wall -msse2 -std=gnu99 -o viterbi_fs_utest -I.. -L.. -I../../easel -L../../easel -Dp7VITERBI_FS_TESTDRIVE viterbi_fs.c -lhmmer -leasel -lm
  ./viterbi_fs_utest
*/
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name    type           default env range toggles reqs incomp help                                       docgroup */
  { "-h",  eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",  eslARG_INT,     "42", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M",  eslARG_INT,     "72", NULL, NULL, NULL, NULL, NULL, "size of random models to sample",                0 },
  { "-N",  eslARG_INT,     "20", NULL, NULL, NULL, NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE p7_Viterbi_Frameshift() and p7_Viterbi_Frameshift_Trace()";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA  = esl_alphabet_Create(eslAMINO);
  ESL_ALPHABET   *abcDNA = esl_alphabet_Create(eslDNA);
  ESL_GENCODE    *gcode  = esl_gencode_Create(abcDNA, abcAA);
  P7_CODONTABLE  *ct     = p7_codontable_Create(gcode);
  P7_BG          *bgAA   = p7_bg_Create(abcAA);
  int             M      = esl_opt_GetInteger(go, "-M");
  int             N      = esl_opt_GetInteger(go, "-N");

  p7_FLogsumInit();

  utest_viterbi_fs(r, abcAA, abcDNA, gcode, bgAA, ct, M, N);

  p7_bg_Destroy(bgAA);
  p7_codontable_Destroy(ct);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7VITERBI_FS_TESTDRIVE*/
/*----------------- end, test driver ----------------------------*/
