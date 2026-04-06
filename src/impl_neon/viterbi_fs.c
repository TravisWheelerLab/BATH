/* NEON-accelerated frameshift-aware Viterbi algorithm; full matrix.
 *
 * Log-space DP: vaddq_f32 along a path, vmaxq_f32 to select the best
 * path at each branch point. No rescaling is required. The profile must
 * be in log-space before calling. The final score is returned directly
 * as a log-odds score in nats:
 *
 *   opt_sc = max(xC(L), xC(L-1)+LOOP, xC(L-2)+LOOP) + xf[C][MOVE]
 *
 * where the max accounts for the three possible final codon lengths
 * (3, 4, or 5 nt) that can end exactly at position L.
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

#include <arm_neon.h>

#include "easel.h"
#include "esl_neon.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_neon.h"

/* IVX intermediate matrix access: [slot 0..p7P_5CODONS-1][stripe 0..Q-1] */
#define IVX(slot, q) (ivx[(slot)][(q)])

/*****************************************************************
 * 1. p7_Viterbi_Frameshift() implementation
 *****************************************************************/

/* Function:  p7_Viterbi_Frameshift()
 * Synopsis:  NEON-accelerated frameshift-aware Viterbi algorithm, log-space, 5 codon lengths.
 *
 * Purpose:   Full O(ML) Viterbi algorithm for frameshift-aware translated
 *            comparison between a nucleotide sequence and a frameshift-aware
 *            codon HMM, using five codon lengths (1-5 nucleotides), with
 *            SIMD parallelism.
 *
 *            Finds the maximum-scoring alignment of nucleotide sequence <dsq>
 *            of length <L> to frameshift-aware codon profile <om_fs>, storing
 *            the full DP matrix in <ox> to enable subsequent traceback.  The
 *            Viterbi score (in nats) is returned in <opt_sc>.
 *
 *            The DP is computed in log-space (log-odds scores).  The profile
 *            <om_fs> must be in log space.
 *
 * Args:      dsq    - nucleotide sequence, 1..L
 *            L      - length of dsq
 *            om_fs  - optimized log space frameshift profile
 *            ox     - DP matrix, p7_omx_Create_dpf(M, L, L, p7X_NSCELLS)
 *            opt_sc - optRETURN: Viterbi score in nats
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on bad inputs; <eslERANGE> if no valid path exists.
 */
int
p7_Viterbi_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_OMX *ox, P7_OIVX *ov, float *opt_sc)
{
  register float32x4_t mpv1, dpv1, ipv1;   /* right-shifted prev1 MDI for BM/MM/IM/DM        */
  register float32x4_t sv;                  /* temporary IVX accumulator                       */
  register float32x4_t msv;                 /* M-state best-path value for current position    */
  register float32x4_t dcv;                 /* delayed D(i,q+1) carry                          */
  register float32x4_t xEv;                 /* E-state partial max (horizontal reduce later)   */
  register float32x4_t xBv1;               /* splatted B(i-1)                                 */
  float32x4_t   infv;                       /* splatted -eslINFINITY: log-space "impossible"   */
  float    xN, xE, xB, xC, xJ;             /* special state scalars (log-space)               */
  float    xN_buf[PARSER_ROWS_FWD];
  float    xB_buf[PARSER_ROWS_FWD];
  float    xJ_buf[PARSER_ROWS_FWD];
  float    xC_buf[PARSER_ROWS_FWD];
  int      b, b1, b3;                       /* circular buffer slots: i, i-1, i-3             */
  int      ivx_1, ivx_2, ivx_3, ivx_4, ivx_5;
  float32x4_t  *dpc, *dpp1, *dpp3;
  float32x4_t  *tp;
  float32x4_t **ivx      = NULL;
  int      Q = p7O_NQF(om_fs->M);
  int      i, q, j, r;
  int      c1, c2, c3, c4, c5;
  int      t, u, v, w, x;
  int      status;

  if (om_fs->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");

  ox->M              = om_fs->M;
  ox->L              = L;
  ox->has_own_scales = FALSE;
  ox->totscale       = 0.0;
  infv = vdupq_n_f32(-eslINFINITY);

  ivx = ov->ivx;

  /* Initialize rows 0..L of the full DP matrix to -inf (log-space impossible). */
  for (r = 0; r <= L; r++)
    for (q = 0; q < Q; q++)
      MMO(ox->dpf[r],q) = DMO(ox->dpf[r],q) = IMO(ox->dpf[r],q) = infv;

  /* Initialize all IVX rows to -inf. */
  for (r = 0; r < p7P_5CODONS; r++)
    for (q = 0; q < Q; q++)
      IVX(r, q) = infv;

  /* Initialize special-state circular buffers (log-space).
   * N(0)=N(1)=N(2)=0 (=log 1); B(0)=B(1)=B(2)=log T_NM; E=J=C=-inf. */
  for (r = 0; r < PARSER_ROWS_FWD; r++)
    xN_buf[r] = xB_buf[r] = xJ_buf[r] = xC_buf[r] = -eslINFINITY;
  xN_buf[0] = xN_buf[1] = xN_buf[2] = 0.0f;
  xB_buf[0] = xB_buf[1] = xB_buf[2] = om_fs->xf[p7O_N][p7O_MOVE];

  /* Write rows 0, 1, 2 specials to ox->xmx (log-space; SCALE=0 means no rescaling). */
  for (r = 0; r < 3; r++)
    {
      ox->xmx[r*p7X_NXCELLS+p7X_SCALE] = 0.0f;
      ox->xmx[r*p7X_NXCELLS+p7X_E]     = -eslINFINITY;
      ox->xmx[r*p7X_NXCELLS+p7X_N]     = 0.0f;
      ox->xmx[r*p7X_NXCELLS+p7X_J]     = -eslINFINITY;
      ox->xmx[r*p7X_NXCELLS+p7X_B]     = om_fs->xf[p7O_N][p7O_MOVE];
      ox->xmx[r*p7X_NXCELLS+p7X_C]     = -eslINFINITY;
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
  dpp1 = ox->dpf[0];  /* all -inf */

  xBv1 = vdupq_n_f32(xB_buf[0]);
  tp   = om_fs->tfv;
  dcv  = infv;
  xEv  = infv;
  mpv1 = dpv1 = ipv1 = infv;

  for (q = 0; q < Q; q++)
    {
      /* IVX(1,q) = max(B(0)+BM, -inf, -inf, -inf) = B(0)+BM */
      sv  =                vaddq_f32(xBv1, *tp); tp++;    /* BM */
      sv  = vmaxq_f32(sv,  vaddq_f32(mpv1, *tp)); tp++;   /* MM=-inf */
      sv  = vmaxq_f32(sv,  vaddq_f32(ipv1, *tp)); tp++;   /* IM=-inf */
      sv  = vmaxq_f32(sv,  vaddq_f32(dpv1, *tp)); tp++;   /* DM=-inf */
      IVX(ivx_1, q) = sv;

      /* M_C0(1,q): 1-nt codon only */
      msv = vaddq_f32(sv, om_fs->rfv[c1][q]);
      MMO(dpc, q) = msv;
      xEv = vmaxq_f32(xEv, msv);

      DMO(dpc, q) = dcv;  /* dcv = -inf initially */
      dcv = vaddq_f32(msv, *tp); tp++;    /* MD */

      /* I(1,q) = -inf (dpp3 not yet available) */
      IMO(dpc, q) = infv;
      tp += 2;  /* skip MI, II */
    }

  /* DD paths */
  dcv        = vextq_f32(infv, dcv, 3);
  DMO(dpc,0) = infv;
  tp         = om_fs->tfv + 7*Q;
  for (q = 0; q < Q; q++)
    {
      DMO(dpc, q) = vmaxq_f32(dcv, DMO(dpc, q));
      dcv         = vaddq_f32(DMO(dpc, q), *tp); tp++;
    }
  if (om_fs->M < 100)
    {
      for (j = 1; j < 4; j++)
        {
          dcv = vextq_f32(infv, dcv, 3);
          tp  = om_fs->tfv + 7*Q;
          for (q = 0; q < Q; q++)
            {
              DMO(dpc, q) = vmaxq_f32(dcv, DMO(dpc, q));
              dcv         = vaddq_f32(dcv, *tp); tp++;
            }
        }
    }
  else
    {
      for (j = 1; j < 4; j++)
        {
          uint32x4_t cv;
          dcv = vextq_f32(infv, dcv, 3);
          tp  = om_fs->tfv + 7*Q;
          cv  = vdupq_n_u32(0);
          for (q = 0; q < Q; q++)
            {
              sv          = vmaxq_f32(dcv, DMO(dpc, q));
              cv          = vorrq_u32(cv, vcgtq_f32(sv, DMO(dpc, q)));
              DMO(dpc, q) = sv;
              dcv         = vaddq_f32(dcv, *tp); tp++;
            }
          if (vgetq_lane_u64(vreinterpretq_u64_u32(cv), 0) == 0 &&
              vgetq_lane_u64(vreinterpretq_u64_u32(cv), 1) == 0) break;
        }
    }
  /* Add D to xEv, then reduce to scalar xE */
  for (q = 0; q < Q; q++)
    xEv = vmaxq_f32(xEv, DMO(dpc, q));
  xE = esl_neon_hmax_f32((esl_neon_128f_t){ .f32x4 = xEv });

  xN = 0.0f;
  xJ = xE + om_fs->xf[p7O_E][p7O_LOOP];
  xC = xE + om_fs->xf[p7O_E][p7O_MOVE];
  xB = ESL_MAX(xN + om_fs->xf[p7O_N][p7O_MOVE], xJ + om_fs->xf[p7O_J][p7O_MOVE]);

  xN_buf[1] = xN; xB_buf[1] = xB; xJ_buf[1] = xJ; xC_buf[1] = xC;
  ox->xmx[1*p7X_NXCELLS+p7X_SCALE] = 0.0f;
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

  xBv1 = vdupq_n_f32(xB_buf[1]);
  tp   = om_fs->tfv;
  dcv  = infv;
  xEv  = infv;

  mpv1 = vextq_f32(infv, MMO(dpp1, Q-1), 3);
  dpv1 = vextq_f32(infv, DMO(dpp1, Q-1), 3);
  ipv1 = vextq_f32(infv, IMO(dpp1, Q-1), 3);

  for (q = 0; q < Q; q++)
    {
      sv  =                vaddq_f32(xBv1, *tp); tp++;
      sv  = vmaxq_f32(sv,  vaddq_f32(mpv1, *tp)); tp++;
      sv  = vmaxq_f32(sv,  vaddq_f32(ipv1, *tp)); tp++;
      sv  = vmaxq_f32(sv,  vaddq_f32(dpv1, *tp)); tp++;
      IVX(ivx_1, q) = sv;

      float32x4_t mc1 = vaddq_f32(sv,             om_fs->rfv[c1][q]);
      float32x4_t mc2 = vaddq_f32(IVX(ivx_2, q), om_fs->rfv[c2][q]);
      msv = vmaxq_f32(mc1, mc2);
      xEv = vmaxq_f32(xEv, msv);

      mpv1 = MMO(dpp1, q);
      dpv1 = DMO(dpp1, q);
      ipv1 = IMO(dpp1, q);

      MMO(dpc, q) = msv;
      DMO(dpc, q) = dcv;

      dcv = vaddq_f32(msv, *tp); tp++;    /* MD */

      /* I(2,q) = -inf (dpp3 not yet available) */
      IMO(dpc, q) = infv;
      tp += 2;  /* skip MI, II */
    }

  /* DD paths */
  dcv        = vextq_f32(infv, dcv, 3);
  DMO(dpc,0) = infv;
  tp         = om_fs->tfv + 7*Q;
  for (q = 0; q < Q; q++)
    {
      DMO(dpc, q) = vmaxq_f32(dcv, DMO(dpc, q));
      dcv         = vaddq_f32(DMO(dpc, q), *tp); tp++;
    }
  if (om_fs->M < 100)
    {
      for (j = 1; j < 4; j++)
        {
          dcv = vextq_f32(infv, dcv, 3);
          tp  = om_fs->tfv + 7*Q;
          for (q = 0; q < Q; q++)
            {
              DMO(dpc, q) = vmaxq_f32(dcv, DMO(dpc, q));
              dcv         = vaddq_f32(dcv, *tp); tp++;
            }
        }
    }
  else
    {
      for (j = 1; j < 4; j++)
        {
          uint32x4_t cv;
          dcv = vextq_f32(infv, dcv, 3);
          tp  = om_fs->tfv + 7*Q;
          cv  = vdupq_n_u32(0);
          for (q = 0; q < Q; q++)
            {
              sv          = vmaxq_f32(dcv, DMO(dpc, q));
              cv          = vorrq_u32(cv, vcgtq_f32(sv, DMO(dpc, q)));
              DMO(dpc, q) = sv;
              dcv         = vaddq_f32(dcv, *tp); tp++;
            }
          if (vgetq_lane_u64(vreinterpretq_u64_u32(cv), 0) == 0 &&
              vgetq_lane_u64(vreinterpretq_u64_u32(cv), 1) == 0) break;
        }
    }
  for (q = 0; q < Q; q++)
    xEv = vmaxq_f32(xEv, DMO(dpc, q));
  xE = esl_neon_hmax_f32((esl_neon_128f_t){ .f32x4 = xEv });

  xN = 0.0f;
  xJ = xE + om_fs->xf[p7O_E][p7O_LOOP];
  xC = xE + om_fs->xf[p7O_E][p7O_MOVE];
  xB = ESL_MAX(xN + om_fs->xf[p7O_N][p7O_MOVE], xJ + om_fs->xf[p7O_J][p7O_MOVE]);

  xN_buf[2] = xN; xB_buf[2] = xB; xJ_buf[2] = xJ; xC_buf[2] = xC;
  ox->xmx[2*p7X_NXCELLS+p7X_SCALE] = 0.0f;
  ox->xmx[2*p7X_NXCELLS+p7X_E] = xE;
  ox->xmx[2*p7X_NXCELLS+p7X_N] = xN;
  ox->xmx[2*p7X_NXCELLS+p7X_J] = xJ;
  ox->xmx[2*p7X_NXCELLS+p7X_B] = xB;
  ox->xmx[2*p7X_NXCELLS+p7X_C] = xC;

  /*----------------------------------------------------------------
   * Main recurrence: i = 3..L
   *
   * M_Cn(i,k) = IVX(i-n+1, k) + Rn(k)              for n = 1..5
   * M_C0(i,k) = max(M_C1, M_C2, M_C3, M_C4, M_C5)  (best codon length)
   * I(i,k)    = max(M(i-3,k)+MI, I(i-3,k)+II)
   * D(i,k)    = max(M(i,k-1)+MD, D(i,k-1)+DD)
   *
   * In log-space there is no scaling, so dpp3 is read directly with
   * no insert_adj correction.
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

      mpv1 = vextq_f32(infv, MMO(dpp1, Q-1), 3);
      dpv1 = vextq_f32(infv, DMO(dpp1, Q-1), 3);
      ipv1 = vextq_f32(infv, IMO(dpp1, Q-1), 3);

      xBv1 = vdupq_n_f32(xB_buf[b1]);
      tp   = om_fs->tfv;
      dcv  = infv;
      xEv  = infv;

      for (q = 0; q < Q; q++)
        {
          /* IVX(i,k) = max(B(i-1)+BM, M(i-1,k-1)+MM, I(i-1,k-1)+IM, D(i-1,k-1)+DM) */
          sv  =                vaddq_f32(xBv1, *tp); tp++;
          sv  = vmaxq_f32(sv,  vaddq_f32(mpv1, *tp)); tp++;
          sv  = vmaxq_f32(sv,  vaddq_f32(ipv1, *tp)); tp++;
          sv  = vmaxq_f32(sv,  vaddq_f32(dpv1, *tp)); tp++;
          IVX(ivx_1, q) = sv;

          /* M_C0 = best codon-length match at this model position */
          msv = vmaxq_f32(vmaxq_f32(vmaxq_f32(vaddq_f32(sv,             om_fs->rfv[c1][q]),
                                               vaddq_f32(IVX(ivx_2, q), om_fs->rfv[c2][q])),
                                    vmaxq_f32(vaddq_f32(IVX(ivx_3, q), om_fs->rfv[c3][q]),
                                               vaddq_f32(IVX(ivx_4, q), om_fs->rfv[c4][q]))),
                                    vaddq_f32(IVX(ivx_5, q),            om_fs->rfv[c5][q]));
          xEv = vmaxq_f32(xEv, msv);

          mpv1 = MMO(dpp1, q);
          dpv1 = DMO(dpp1, q);
          ipv1 = IMO(dpp1, q);

          MMO(dpc, q) = msv;
          DMO(dpc, q) = dcv;
          dcv = vaddq_f32(msv, *tp); tp++;    /* MD */

          /* I(i,k) = max(M(i-3,k)+MI, I(i-3,k)+II); no scale correction in log-space */
          sv          = vaddq_f32(MMO(dpp3, q), *tp); tp++;  /* MI */
          IMO(dpc, q) = vmaxq_f32(sv, vaddq_f32(IMO(dpp3, q), *tp)); tp++;  /* II */
        }

      /* DD paths */
      dcv        = vextq_f32(infv, dcv, 3);
      DMO(dpc,0) = infv;
      tp         = om_fs->tfv + 7*Q;
      for (q = 0; q < Q; q++)
        {
          DMO(dpc, q) = vmaxq_f32(dcv, DMO(dpc, q));
          dcv         = vaddq_f32(DMO(dpc, q), *tp); tp++;
        }
      if (om_fs->M < 100)
        {
          for (j = 1; j < 4; j++)
            {
              dcv = vextq_f32(infv, dcv, 3);
              tp  = om_fs->tfv + 7*Q;
              for (q = 0; q < Q; q++)
                {
                  DMO(dpc, q) = vmaxq_f32(dcv, DMO(dpc, q));
                  dcv         = vaddq_f32(dcv, *tp); tp++;
                }
            }
        }
      else
        {
          for (j = 1; j < 4; j++)
            {
              uint32x4_t cv;
              dcv = vextq_f32(infv, dcv, 3);
              tp  = om_fs->tfv + 7*Q;
              cv  = vdupq_n_u32(0);
              for (q = 0; q < Q; q++)
                {
                  sv          = vmaxq_f32(dcv, DMO(dpc, q));
                  cv          = vorrq_u32(cv, vcgtq_f32(sv, DMO(dpc, q)));
                  DMO(dpc, q) = sv;
                  dcv         = vaddq_f32(dcv, *tp); tp++;
                }
              if (vgetq_lane_u64(vreinterpretq_u64_u32(cv), 0) == 0 &&
                  vgetq_lane_u64(vreinterpretq_u64_u32(cv), 1) == 0) break;
            }
        }

      for (q = 0; q < Q; q++)
        xEv = vmaxq_f32(xEv, DMO(dpc, q));
      xE = esl_neon_hmax_f32((esl_neon_128f_t){ .f32x4 = xEv });

      xN = xN_buf[b3] + om_fs->xf[p7O_N][p7O_LOOP];
      xJ = ESL_MAX(xJ_buf[b3] + om_fs->xf[p7O_J][p7O_LOOP], xE + om_fs->xf[p7O_E][p7O_LOOP]);
      xC = ESL_MAX(xC_buf[b3] + om_fs->xf[p7O_C][p7O_LOOP], xE + om_fs->xf[p7O_E][p7O_MOVE]);
      xB = ESL_MAX(xN         + om_fs->xf[p7O_N][p7O_MOVE],  xJ + om_fs->xf[p7O_J][p7O_MOVE]);

      xN_buf[b] = xN; xB_buf[b] = xB; xJ_buf[b] = xJ; xC_buf[b] = xC;
      ox->xmx[i*p7X_NXCELLS+p7X_SCALE] = 0.0f;
      ox->xmx[i*p7X_NXCELLS+p7X_E] = xE;
      ox->xmx[i*p7X_NXCELLS+p7X_N] = xN;
      ox->xmx[i*p7X_NXCELLS+p7X_J] = xJ;
      ox->xmx[i*p7X_NXCELLS+p7X_B] = xB;
      ox->xmx[i*p7X_NXCELLS+p7X_C] = xC;
    } /* end main loop i=3..L */

  /* Final score: max of xC(L), xC(L-1)+LOOP, xC(L-2)+LOOP, then +MOVE.
   * The best-scoring alignment can end with a codon of length 3, 4, or 5 nt. */
  {
    float xCL   = xC_buf[   L    % PARSER_ROWS_FWD];
    float xCLm1 = xC_buf[((L-1) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD];
    float xCLm2 = xC_buf[((L-2) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD];
    float xCtot = ESL_MAX(xCL,
                  ESL_MAX(xCLm1 + om_fs->xf[p7O_C][p7O_LOOP],
                          xCLm2 + om_fs->xf[p7O_C][p7O_LOOP]));

    if      (isnan(xCtot))            ESL_EXCEPTION(eslERANGE, "Viterbi score is NaN");
    else if (isinf(xCtot) == 1)       ESL_EXCEPTION(eslERANGE, "Viterbi score overflow (is +infinity)");
    else if (L > 1 && isinf(xCtot) == -1) {
      if (opt_sc != NULL) *opt_sc = -eslINFINITY;
      return eslERANGE;
    }

    if (opt_sc != NULL)
      *opt_sc = xCtot + om_fs->xf[p7O_C][p7O_MOVE];
  }

  return eslOK;

 ERROR:
  return status;
}
/*------------------ end, p7_Viterbi_Frameshift() ---------------*/



/*****************************************************************
 * 2. p7_Viterbi_Frameshift_Trace() implementation
 *****************************************************************/

/* Function:  p7_Viterbi_Frameshift_Trace()
 * Synopsis:  Traceback from a Viterbi DP matrix (FS layout); NEON version.
 *
 * Purpose:   Extract the optimal (Viterbi) alignment path from the DP matrix
 *            <ox> that was filled by <p7_Viterbi_Frameshift()>.  The trace is
 *            written into <tr> (caller provides an already-allocated trace via
 *            <p7_trace_fs_Create()>).
 *
 * Args:      dsq   - digital nucleotide sequence, 1..L
 *            L     - sequence length
 *            om_fs - optimized frameshift profile
 *            ox    - Viterbi DP matrix from p7_Viterbi_Frameshift()
 *            tr    - RETURN: traceback; caller provides initial alloc
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslFAIL> if an impossible state is reached during traceback.
 */
int
p7_Viterbi_Frameshift_Trace(const ESL_DSQ *dsq, int L,
                            const P7_FS_OPROFILE *om_fs, const P7_OMX *ox,
                            P7_TRACE *tr)
{
  int   Q      = p7O_NQF(om_fs->M);
  int   M      = om_fs->M;
  int   i      = L;
  int   k      = 0;
  int   c      = 0;
  int   prev_c = 0;
  float r_tol  = 1e-5;
  float a_tol  = 1e-4;
  int   sprv, scur;
  float path[4];
  int   state4[4] = { p7T_B, p7T_M, p7T_I, p7T_D };
  float match_codon[5];
  int   status;

#define OMMo(i,k)  ((k)<1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_M, (i), (k)))
#define ODMo(i,k)  ((k)<1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_D, (i), (k)))
#define OIMo(i,k)  ((k)<1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_I, (i), (k)))
#define OXMXo(i,s)  (ox->xmx[(i)*p7X_NXCELLS+(s)])

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i, c)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i, c)) != eslOK) return status;
  sprv = p7T_C;

  while (sprv != p7T_S) {
    switch (sprv) {
    case p7T_C:     /* C(i) comes from C(i-3) or E(i) */
      if (OXMXo(i, p7X_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);

      if      (OXMXo(i, p7X_C) < OXMXo(i-2, p7X_C) || OXMXo(i, p7X_C) < OXMXo(i-1, p7X_C))                             scur = p7T_C;
      else if (esl_FCompare(OXMXo(i, p7X_C), OXMXo(i-3, p7X_C) + om_fs->xf[p7O_C][p7O_LOOP], r_tol, a_tol) == eslOK)  scur = p7T_C;
      else if (esl_FCompare(OXMXo(i, p7X_C), OXMXo(i,   p7X_E) + om_fs->xf[p7O_E][p7O_MOVE], r_tol, a_tol) == eslOK)  scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7T_E:     /* E connects from any M state; k set here */
      if (OXMXo(i, p7X_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);

      if (p7_fs_oprofile_IsLocal(om_fs))
        {
          scur = p7T_M;     /* can't come from D in local mode */
          for (k = M; k >= 1; k--) if (esl_FCompare(OXMXo(i, p7X_E), OMMo(i,k), r_tol, a_tol) == eslOK) break;
          if (k == 0) ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
        }
      else              /* glocal: only M_M or D_M */
        {
          if      (esl_FCompare(OXMXo(i, p7X_E), OMMo(i,M), r_tol, a_tol) == eslOK) { scur = p7T_M; k = M; }
          else if (esl_FCompare(OXMXo(i, p7X_E), ODMo(i,M), r_tol, a_tol) == eslOK) { scur = p7T_D; k = M; }
          else    ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
        }
      break;

    case p7T_M:     /* M(i,k) comes from B/M/I/D at row i-prev_c, pos k-1 */
      {
        union { float32x4_t v; float p[4]; } u;
        int   q_k   = (k-1) % Q;
        int   r_k   = (k-1) / Q;
        int   q_km1 = (k > 1) ? (k-2) % Q : 0;
        int   r_km1 = (k > 1) ? (k-2) / Q : 0;
        float tbm, tmm, tim, tdm;
        float bprev, mprev, iprev, dprev;
        int   ipred = i - prev_c;

        if (OMMo(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k, i);

        u.v = om_fs->tfv[7*q_k + p7O_BM]; tbm = u.p[r_k];
        u.v = om_fs->tfv[7*q_k + p7O_MM]; tmm = u.p[r_k];
        u.v = om_fs->tfv[7*q_k + p7O_IM]; tim = u.p[r_k];
        u.v = om_fs->tfv[7*q_k + p7O_DM]; tdm = u.p[r_k];

        bprev = OXMXo(ipred, p7X_B);
        if (k > 1) {
          u.v = MMO(ox->dpf[ipred], q_km1); mprev = u.p[r_km1];
          u.v = IMO(ox->dpf[ipred], q_km1); iprev = u.p[r_km1];
          u.v = DMO(ox->dpf[ipred], q_km1); dprev = u.p[r_km1];
        } else {
          mprev = iprev = dprev = -eslINFINITY;
        }

        path[0] = bprev + tbm;
        path[1] = mprev + tmm;
        path[2] = iprev + tim;
        path[3] = dprev + tdm;
        scur = state4[esl_vec_FArgMax(path, 4)];
        k--; i -= prev_c;
        break;
      }

    case p7T_D:     /* D(i,k) comes from M(i,k-1) or D(i,k-1) */
      {
        union { float32x4_t v; float p[4]; } u;
        int   q_km1 = (k-2) % Q;
        int   r_km1 = (k-2) / Q;
        float tmd, tdd, mm_km1, dm_km1;

        if (ODMo(i, k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k, i);
        if (k < 2) ESL_EXCEPTION(eslFAIL, "D at k=1,i=%d can't be traced", i);

        u.v = om_fs->tfv[7*q_km1 + p7O_MD]; tmd = u.p[r_km1];
        u.v = om_fs->tfv[7*Q + q_km1];      tdd = u.p[r_km1];
        u.v = MMO(ox->dpf[i], q_km1); mm_km1 = u.p[r_km1];
        u.v = DMO(ox->dpf[i], q_km1); dm_km1 = u.p[r_km1];

        if      (esl_FCompare(ODMo(i,k), mm_km1 + tmd, r_tol, a_tol) == eslOK) scur = p7T_M;
        else if (esl_FCompare(ODMo(i,k), dm_km1 + tdd, r_tol, a_tol) == eslOK) scur = p7T_D;
        else ESL_EXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k, i);
        k--;
        break;
      }

    case p7T_I:     /* I(i,k) comes from M(i-3,k) or I(i-3,k) */
      {
        union { float32x4_t v; float p[4]; } u;
        int   q_k = (k-1) % Q;
        int   r_k = (k-1) / Q;
        float tmi, tii, mprev_i, iprev_i;

        if (OIMo(i, k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k, i);

        u.v = om_fs->tfv[7*q_k + p7O_MI]; tmi = u.p[r_k];
        u.v = om_fs->tfv[7*q_k + p7O_II]; tii = u.p[r_k];
        u.v = MMO(ox->dpf[i-3], q_k); mprev_i = u.p[r_k];
        u.v = IMO(ox->dpf[i-3], q_k); iprev_i = u.p[r_k];

        if      (esl_FCompare(OIMo(i,k), mprev_i + tmi, r_tol, a_tol) == eslOK) scur = p7T_M;
        else if (esl_FCompare(OIMo(i,k), iprev_i + tii, r_tol, a_tol) == eslOK) scur = p7T_I;
        else ESL_EXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k, i);
        i -= 3;
        break;
      }

    case p7T_N:     /* N connects from S or N */
      if (OXMXo(i, p7X_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = (i == 0) ? p7T_S : p7T_N;
      break;

    case p7T_B:     /* B connects from N or J */
      if (OXMXo(i, p7X_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

      if      (esl_FCompare(OXMXo(i, p7X_B), OXMXo(i, p7X_N) + om_fs->xf[p7O_N][p7O_MOVE], r_tol, a_tol) == eslOK) scur = p7T_N;
      else if (esl_FCompare(OXMXo(i, p7X_B), OXMXo(i, p7X_J) + om_fs->xf[p7O_J][p7O_MOVE], r_tol, a_tol) == eslOK) scur = p7T_J;
      else ESL_EXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    case p7T_J:     /* J connects from E(i) or J(i-3) */
      if (OXMXo(i, p7X_J) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible J reached at i=%d", i);

      if      (esl_FCompare(OXMXo(i, p7X_J), OXMXo(i-3, p7X_J) + om_fs->xf[p7O_J][p7O_LOOP], r_tol, a_tol) == eslOK) scur = p7T_J;
      else if (esl_FCompare(OXMXo(i, p7X_J), OXMXo(i,   p7X_E) + om_fs->xf[p7O_E][p7O_LOOP], r_tol, a_tol) == eslOK) scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "J at i=%d couldn't be traced", i);
      break;

    default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over sprv */

    /* For M states: recompute codon length c (1..5) from stored DP cells.
     * At this point k and i have already been decremented to the predecessor
     * position, so we compute the codon length for the M state being appended
     * (scur == p7T_M) at (k, i).
     */
    if (scur == p7T_M) {
      union { float32x4_t v; float p[4]; } u;
      int   q_k   = (k-1) % Q;
      int   r_k   = (k-1) / Q;
      int   q_km1 = (k > 1) ? (k-2) % Q : 0;
      int   r_km1 = (k > 1) ? (k-2) / Q : 0;
      float tbm, tmm, tim, tdm;
      int   x2 = (i >= 1 && dsq[i]   < p7P_MAXNUC) ? (int)dsq[i]   : p7P_MAXCODONS5;
      int   w2 = (i >= 2 && dsq[i-1] < p7P_MAXNUC) ? (int)dsq[i-1] : p7P_MAXCODONS5;
      int   v2 = (i >= 3 && dsq[i-2] < p7P_MAXNUC) ? (int)dsq[i-2] : p7P_MAXCODONS5;
      int   u2 = (i >= 4 && dsq[i-3] < p7P_MAXNUC) ? (int)dsq[i-3] : p7P_MAXCODONS5;
      int   t2 = (i >= 5 && dsq[i-4] < p7P_MAXNUC) ? (int)dsq[i-4] : p7P_MAXCODONS5;
      int   cn, n;
      float bprev, mprev, iprev, dprev, ivx_n;

      u.v = om_fs->tfv[7*q_k + p7O_BM]; tbm = u.p[r_k];
      u.v = om_fs->tfv[7*q_k + p7O_MM]; tmm = u.p[r_k];
      u.v = om_fs->tfv[7*q_k + p7O_IM]; tim = u.p[r_k];
      u.v = om_fs->tfv[7*q_k + p7O_DM]; tdm = u.p[r_k];

      for (n = 1; n <= 5; n++) {
        if (i < n) { match_codon[n-1] = -eslINFINITY; continue; }

        bprev = OXMXo(i-n, p7X_B);
        if (k > 1) {
          u.v = MMO(ox->dpf[i-n], q_km1); mprev = u.p[r_km1];
          u.v = IMO(ox->dpf[i-n], q_km1); iprev = u.p[r_km1];
          u.v = DMO(ox->dpf[i-n], q_km1); dprev = u.p[r_km1];
        } else {
          mprev = iprev = dprev = -eslINFINITY;
        }
        ivx_n = ESL_MAX(bprev + tbm, ESL_MAX(mprev + tmm, ESL_MAX(iprev + tim, dprev + tdm)));

        switch (n) {
          case 1: cn = p7P_CODON1_FS5(x2);                  cn = p7P_MINIDX(cn, p7P_DEGEN5_QC2); break;
          case 2: cn = p7P_CODON2_FS5(w2, x2);              cn = p7P_MINIDX(cn, p7P_DEGEN5_QC1); break;
          case 3: cn = p7P_CODON3_FS5(v2, w2, x2);          cn = p7P_MINIDX(cn, p7P_DEGEN5_C);   break;
          case 4: cn = p7P_CODON4_FS5(u2, v2, w2, x2);      cn = p7P_MINIDX(cn, p7P_DEGEN5_QC1); break;
          default:cn = p7P_CODON5_FS5(t2, u2, v2, w2, x2); cn = p7P_MINIDX(cn, p7P_DEGEN5_QC2); break;
        }
        u.v = om_fs->rfv[cn][q_k];
        match_codon[n-1] = ivx_n + u.p[r_k];
      }
      c = esl_vec_FArgMax(match_codon, 5) + 1;
    } else {
      c = 0;
    }

    if ((status = p7_trace_fs_Append(tr, scur, k, i, c)) != eslOK) return status;

    /* Deferred i decrements for N, C (by 1) and J (by 3). */
    if ((scur == p7T_N || scur == p7T_C) && scur == sprv) i--;
    if ( scur == p7T_J                   && scur == sprv) i -= 3;

    prev_c = c;
    c = 0;
    sprv = scur;
  } /* end traceback, at S state */

#undef OMMo
#undef ODMo
#undef OIMo
#undef OXMXo

  tr->M = om_fs->M;
  tr->L = L;
  return p7_trace_fs_Reverse(tr);
}
/*------------------ end, p7_Viterbi_Frameshift_Trace() --------------------*/



/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/
#ifdef p7VITERBI_FS_BENCHMARK
/*
   gcc -g -O3 -march=armv8-a -std=gnu99 -o viterbi_fs_benchmark -I.. -L.. -I../../easel -L../../easel -Dp7VITERBI_FS_BENCHMARK viterbi_fs.c -lhmmer -leasel -lm

   ./viterbi_fs_benchmark <hmmfile>            benchmark Viterbi + traceback together
   ./viterbi_fs_benchmark -V <hmmfile>         benchmark Viterbi only
   ./viterbi_fs_benchmark -T <hmmfile>         benchmark traceback only
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
#include "impl_neon.h"

static ESL_OPTIONS options[] = {
  /* name    type           default  env  range  toggles reqs incomp  help                                             docgroup */
  { "-h",  eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL,  "show brief help on version and usage",                0 },
  { "-s",  eslARG_INT,     "42", NULL, NULL,   NULL,  NULL, NULL,  "set random number seed to <n>",                       0 },
  { "-L",  eslARG_INT,   "1200", NULL, "n>0",  NULL,  NULL, NULL,  "length of random target DNA seqs (multiple of 3)",    0 },
  { "-N",  eslARG_INT,   "2000", NULL, "n>0",  NULL,  NULL, NULL,  "number of random target seqs",                        0 },
  { "-V",  eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, "-T",  "benchmark Viterbi only (skip traceback)",             0 },
  { "-T",  eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, "-V",  "benchmark traceback only (skip Viterbi timing)",      0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for NEON frameshift Viterbi";

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
  p7_fs_oprofile_Convert_Log(gm_fs5, om_fs5);
  p7_fs_oprofile_ReconfigLength_Log(om_fs5, L/3);

  ox = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS);
  tr = p7_trace_fs_Create();
  P7_OIVX *ov5 = p7_oivx_Create(hmm->M, p7P_5CODONS);

  /* If benchmarking traceback only, pre-fill the matrix once. */
  if (! do_viterbi) {
    esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
    p7_Viterbi_Frameshift(dsq, L, om_fs5, ox, ov5, &sc);
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
        p7_Viterbi_Frameshift(dsq, L, om_fs5, ox, ov5, &sc);

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
  p7_oivx_Destroy(ov5);
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
 *   1. NEON and generic Viterbi scores agree (< 0.001 nat tolerance).
 *   2. Each NEON trace is structurally valid.
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
  P7_OIVX        *ov5    = p7_oivx_Create(M, p7P_5CODONS);
  P7_GMX         *gx     = p7_gmx_Create(M, M, M, p7G_NSCELLS);
  P7_OMX         *ox     = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS);
  char            errbuf[eslERRBUFSIZE];
  float           gsc, osc;
  int             curr_L, i, j;

  if (!gm || !gm_fs5 || !om_fs5 || !sq || !tr || !trg || !iv5 || !ov5 || !gx || !ox) esl_fatal(msg);

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig(hmm, bgAA, gm, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, M, p7_LOCAL);
  p7_fs_oprofile_Convert_Log(gm_fs5, om_fs5);

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

      p7_fs_oprofile_ReconfigLength_Log(om_fs5, sq->n);
      p7_fs_ReconfigLength(gm_fs5, sq->n);

      p7_gmx_GrowTo(gx, M, curr_L, curr_L);
      p7_ivx_GrowTo(iv5, M, p7P_5CODONS);
      p7_oivx_GrowTo(ov5, M, p7P_5CODONS);
      p7_omx_GrowTo_dpf(ox, M, curr_L, curr_L);

      /* Generic Viterbi */
      if (p7_GViterbi_Frameshift(dsq, curr_L, gm_fs5, gx, iv5, &gsc) != eslOK) esl_fatal(msg);

      /* NEON Viterbi */
      { int s = p7_Viterbi_Frameshift(dsq, curr_L, om_fs5, ox, ov5, &osc);
        if (s == eslERANGE) continue;
        if (s != eslOK)     esl_fatal(msg); }

      /* Scores must agree */
      if (fabs(gsc - osc) > 0.001)
        esl_fatal("%s: generic vit %.4f != NEON vit %.4f", msg, gsc, osc);

      /* Traces */
      if (p7_GVTrace_Frameshift(dsq, curr_L, gm_fs5, gx, trg)           != eslOK) esl_fatal(msg);
      if (p7_Viterbi_Frameshift_Trace(dsq, curr_L, om_fs5, ox, tr)       != eslOK) esl_fatal(msg);

      if (p7_trace_fs_Validate(tr,  abcDNA, dsq, errbuf) != eslOK)
        esl_fatal("%s: NEON trace invalid: %s",    msg, errbuf);
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
  p7_oivx_Destroy(ov5);
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
  gcc -g -Wall -march=armv8-a -std=gnu99 -o viterbi_fs_utest -I.. -L.. -I../../easel -L../../easel -Dp7VITERBI_FS_TESTDRIVE viterbi_fs.c -lhmmer -leasel -lm
  ./viterbi_fs_utest
*/
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_neon.h"

static ESL_OPTIONS options[] = {
  /* name    type           default env range toggles reqs incomp help                                       docgroup */
  { "-h",  eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",  eslARG_INT,     "42", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M",  eslARG_INT,     "72", NULL, NULL, NULL, NULL, NULL, "size of random models to sample",                0 },
  { "-N",  eslARG_INT,     "20", NULL, NULL, NULL, NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for NEON p7_Viterbi_Frameshift() and p7_Viterbi_Frameshift_Trace()";

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
