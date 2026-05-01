/* SSE frameshift-aware Viterbi algorithm; full matrix.
 * Ported from impl_sse/viterbi_fs.c for the impl_avx runtime-dispatch model.
 * All public functions carry the _sse suffix; the non-suffixed dispatchers
 * live in viterbi_fs.c.
 *
 * Contents:
 *   1. p7_Viterbi_Frameshift_sse() implementation.
 *   2. p7_Viterbi_Frameshift_Trace_sse() implementation.
 */
#include "p7_config.h"

#ifdef eslENABLE_SSE

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"

/* IVX intermediate matrix access: [slot 0..p7P_5CODONS-1][stripe 0..Q-1] */
#define IVX(slot, q) (ivx[(slot)][(q)])

/*****************************************************************
 * 1. p7_Viterbi_Frameshift_sse() implementation
 *****************************************************************/

/* Function:  p7_Viterbi_Frameshift_sse()
 * Synopsis:  SSE frameshift-aware Viterbi algorithm, log-space, 5 codon lengths.
 *
 * Purpose:   SSE implementation; see p7_Viterbi_Frameshift() in viterbi_fs.c
 *            for full documentation.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on bad inputs; <eslERANGE> if no valid path exists.
 */
int
p7_Viterbi_Frameshift_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_OMX *ox, P7_OIVX *ov, float *opt_sc)
{
  register __m128 mpv1, dpv1, ipv1;     /* right-shifted prev1 MDI for BM/MM/IM/DM        */
  register __m128 sv;                    /* temporary IVX accumulator                       */
  register __m128 msv;                   /* M-state best-path value for current position    */
  register __m128 dcv;                   /* delayed D(i,q+1) carry                          */
  register __m128 xEv;                   /* E-state partial max (horizontal reduce later)   */
  register __m128 xBv1;                  /* splatted B(i-1)                                 */
  __m128   infv;                         /* splatted -eslINFINITY: log-space "impossible"   */
  float    xN, xE, xB, xC, xJ;          /* special state scalars (log-space)               */
  float    xN_buf[PARSER_ROWS_FWD];
  float    xB_buf[PARSER_ROWS_FWD];
  float    xJ_buf[PARSER_ROWS_FWD];
  float    xC_buf[PARSER_ROWS_FWD];
  int      b, b1, b3;                    /* circular buffer slots: i, i-1, i-3             */
  int      ivx_1, ivx_2, ivx_3, ivx_4, ivx_5;
  __m128  *dpc, *dpp1, *dpp3;
  __m128  *tp;
  __m128 **ivx      = NULL;
  int      Q = p7O_NQF(om_fs->M);
  int      i, q, j, r;
  int      c1, c2, c3, c4, c5;
  int      t, u, v, w, x;

  if (om_fs->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");

  ox->M              = om_fs->M;
  ox->L              = L;
  ox->has_own_scales = FALSE;
  ox->totscale       = 0.0;
  infv = _mm_set1_ps(-eslINFINITY);

  ivx = ov->ivx;

  { __m128 *drow = ox->dpf[0];
    for (q = 0; q < Q; q++)
      MMO(drow, q) = DMO(drow, q) = IMO(drow, q) = infv; }

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

  xBv1 = _mm_set1_ps(xB_buf[0]);
  tp   = om_fs->tfv;
  dcv  = infv;
  xEv  = infv;
  mpv1 = dpv1 = ipv1 = infv;

  for (q = 0; q < Q; q++)
    {
      /* IVX(1,q) = max(B(0)+BM, -inf, -inf, -inf) = B(0)+BM */
      sv  =                _mm_add_ps(xBv1, *tp); tp++;   /* BM */
      sv  = _mm_max_ps(sv, _mm_add_ps(mpv1, *tp)); tp++;  /* MM=-inf */
      sv  = _mm_max_ps(sv, _mm_add_ps(ipv1, *tp)); tp++;  /* IM=-inf */
      sv  = _mm_max_ps(sv, _mm_add_ps(dpv1, *tp)); tp++;  /* DM=-inf */
      IVX(ivx_1, q) = sv;

      /* M_C0(1,q): 1-nt codon only */
      msv = _mm_add_ps(sv, om_fs->rfv[c1][q]);
      MMO(dpc, q) = msv;
      xEv = _mm_max_ps(xEv, msv);

      DMO(dpc, q) = dcv;  /* dcv = -inf initially */
      dcv = _mm_add_ps(msv, *tp); tp++;   /* MD */

      /* I(1,q) = -inf (dpp3 not yet available) */
      IMO(dpc, q) = infv;
      tp += 2;  /* skip MI, II */
    }

  /* DD paths */
  dcv        = esl_sse_rightshift_ps(dcv, infv);
  DMO(dpc,0) = infv;
  tp         = om_fs->tfv + 7*Q;
  for (q = 0; q < Q; q++)
    {
      DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
      dcv         = _mm_add_ps(DMO(dpc, q), *tp); tp++;
    }
  if (om_fs->M < 100)
    {
      for (j = 1; j < 4; j++)
        {
          dcv = esl_sse_rightshift_ps(dcv, infv);
          tp  = om_fs->tfv + 7*Q;
          for (q = 0; q < Q; q++)
            {
              DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
              dcv         = _mm_add_ps(dcv, *tp); tp++;
            }
        }
    }
  else
    {
      for (j = 1; j < 4; j++)
        {
          register __m128 cv;
          dcv = esl_sse_rightshift_ps(dcv, infv);
          tp  = om_fs->tfv + 7*Q;
          cv  = _mm_setzero_ps();
          for (q = 0; q < Q; q++)
            {
              sv          = _mm_max_ps(dcv, DMO(dpc, q));
              cv          = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO(dpc, q)));
              DMO(dpc, q) = sv;
              dcv         = _mm_add_ps(dcv, *tp); tp++;
            }
          if (! _mm_movemask_ps(cv)) break;
        }
    }
  /* Add D to xEv, then reduce to scalar xE */
  for (q = 0; q < Q; q++)
    xEv = _mm_max_ps(xEv, DMO(dpc, q));
  esl_sse_hmax_ps(xEv, &xE);

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

  xBv1 = _mm_set1_ps(xB_buf[1]);
  tp   = om_fs->tfv;
  dcv  = infv;
  xEv  = infv;

  mpv1 = esl_sse_rightshift_ps(MMO(dpp1, Q-1), infv);
  dpv1 = esl_sse_rightshift_ps(DMO(dpp1, Q-1), infv);
  ipv1 = esl_sse_rightshift_ps(IMO(dpp1, Q-1), infv);

  for (q = 0; q < Q; q++)
    {
      sv  =                _mm_add_ps(xBv1, *tp); tp++;
      sv  = _mm_max_ps(sv, _mm_add_ps(mpv1, *tp)); tp++;
      sv  = _mm_max_ps(sv, _mm_add_ps(ipv1, *tp)); tp++;
      sv  = _mm_max_ps(sv, _mm_add_ps(dpv1, *tp)); tp++;
      IVX(ivx_1, q) = sv;

      __m128 mc1 = _mm_add_ps(sv,             om_fs->rfv[c1][q]);
      __m128 mc2 = _mm_add_ps(IVX(ivx_2, q), om_fs->rfv[c2][q]);
      msv = _mm_max_ps(mc1, mc2);
      xEv = _mm_max_ps(xEv, msv);

      mpv1 = MMO(dpp1, q);
      dpv1 = DMO(dpp1, q);
      ipv1 = IMO(dpp1, q);

      MMO(dpc, q) = msv;
      DMO(dpc, q) = dcv;

      dcv = _mm_add_ps(msv, *tp); tp++;   /* MD */

      /* I(2,q) = -inf (dpp3 not yet available) */
      IMO(dpc, q) = infv;
      tp += 2;  /* skip MI, II */
    }

  /* DD paths */
  dcv        = esl_sse_rightshift_ps(dcv, infv);
  DMO(dpc,0) = infv;
  tp         = om_fs->tfv + 7*Q;
  for (q = 0; q < Q; q++)
    {
      DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
      dcv         = _mm_add_ps(DMO(dpc, q), *tp); tp++;
    }
  if (om_fs->M < 100)
    {
      for (j = 1; j < 4; j++)
        {
          dcv = esl_sse_rightshift_ps(dcv, infv);
          tp  = om_fs->tfv + 7*Q;
          for (q = 0; q < Q; q++)
            {
              DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
              dcv         = _mm_add_ps(dcv, *tp); tp++;
            }
        }
    }
  else
    {
      for (j = 1; j < 4; j++)
        {
          register __m128 cv;
          dcv = esl_sse_rightshift_ps(dcv, infv);
          tp  = om_fs->tfv + 7*Q;
          cv  = _mm_setzero_ps();
          for (q = 0; q < Q; q++)
            {
              sv          = _mm_max_ps(dcv, DMO(dpc, q));
              cv          = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO(dpc, q)));
              DMO(dpc, q) = sv;
              dcv         = _mm_add_ps(dcv, *tp); tp++;
            }
          if (! _mm_movemask_ps(cv)) break;
        }
    }
  for (q = 0; q < Q; q++)
    xEv = _mm_max_ps(xEv, DMO(dpc, q));
  esl_sse_hmax_ps(xEv, &xE);

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

      mpv1 = esl_sse_rightshift_ps(MMO(dpp1, Q-1), infv);
      dpv1 = esl_sse_rightshift_ps(DMO(dpp1, Q-1), infv);
      ipv1 = esl_sse_rightshift_ps(IMO(dpp1, Q-1), infv);

      xBv1 = _mm_set1_ps(xB_buf[b1]);
      tp   = om_fs->tfv;
      dcv  = infv;
      xEv  = infv;

      for (q = 0; q < Q; q++)
        {
          sv  =                _mm_add_ps(xBv1, *tp); tp++;
          sv  = _mm_max_ps(sv, _mm_add_ps(mpv1, *tp)); tp++;
          sv  = _mm_max_ps(sv, _mm_add_ps(ipv1, *tp)); tp++;
          sv  = _mm_max_ps(sv, _mm_add_ps(dpv1, *tp)); tp++;
          IVX(ivx_1, q) = sv;

          msv = _mm_max_ps(_mm_max_ps(_mm_max_ps(_mm_add_ps(sv,             om_fs->rfv[c1][q]),
                                                  _mm_add_ps(IVX(ivx_2, q), om_fs->rfv[c2][q])),
                                      _mm_max_ps(_mm_add_ps(IVX(ivx_3, q), om_fs->rfv[c3][q]),
                                                  _mm_add_ps(IVX(ivx_4, q), om_fs->rfv[c4][q]))),
                                      _mm_add_ps(IVX(ivx_5, q),             om_fs->rfv[c5][q]));
          xEv = _mm_max_ps(xEv, msv);

          mpv1 = MMO(dpp1, q);
          dpv1 = DMO(dpp1, q);
          ipv1 = IMO(dpp1, q);

          MMO(dpc, q) = msv;
          DMO(dpc, q) = dcv;
          dcv = _mm_add_ps(msv, *tp); tp++;   /* MD */

          sv          = _mm_add_ps(MMO(dpp3, q), *tp); tp++;  /* MI */
          IMO(dpc, q) = _mm_max_ps(sv, _mm_add_ps(IMO(dpp3, q), *tp)); tp++;  /* II */
        }

      /* DD paths */
      dcv        = esl_sse_rightshift_ps(dcv, infv);
      DMO(dpc,0) = infv;
      tp         = om_fs->tfv + 7*Q;
      for (q = 0; q < Q; q++)
        {
          DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
          dcv         = _mm_add_ps(DMO(dpc, q), *tp); tp++;
        }
      if (om_fs->M < 100)
        {
          for (j = 1; j < 4; j++)
            {
              dcv = esl_sse_rightshift_ps(dcv, infv);
              tp  = om_fs->tfv + 7*Q;
              for (q = 0; q < Q; q++)
                {
                  DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
                  dcv         = _mm_add_ps(dcv, *tp); tp++;
                }
            }
        }
      else
        {
          for (j = 1; j < 4; j++)
            {
              register __m128 cv;
              dcv = esl_sse_rightshift_ps(dcv, infv);
              tp  = om_fs->tfv + 7*Q;
              cv  = _mm_setzero_ps();
              for (q = 0; q < Q; q++)
                {
                  sv          = _mm_max_ps(dcv, DMO(dpc, q));
                  cv          = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO(dpc, q)));
                  DMO(dpc, q) = sv;
                  dcv         = _mm_add_ps(dcv, *tp); tp++;
                }
              if (! _mm_movemask_ps(cv)) break;
            }
        }

      for (q = 0; q < Q; q++)
        xEv = _mm_max_ps(xEv, DMO(dpc, q));
      esl_sse_hmax_ps(xEv, &xE);

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

  /* Final score */
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
}
/*------------------ end, p7_Viterbi_Frameshift_sse() ---------------*/



/*****************************************************************
 * 2. p7_Viterbi_Frameshift_Trace_sse() implementation
 *****************************************************************/

/* Function:  p7_Viterbi_Frameshift_Trace_sse()
 * Synopsis:  Traceback from a Viterbi DP matrix (FS layout); SSE version.
 *
 * Purpose:   SSE implementation; see p7_Viterbi_Frameshift_Trace() in viterbi_fs.c
 *            for full documentation.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslFAIL> if an impossible state is reached during traceback.
 */
int
p7_Viterbi_Frameshift_Trace_sse(const ESL_DSQ *dsq, int L,
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
    case p7T_C:
      if (OXMXo(i, p7X_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);

      if      (OXMXo(i, p7X_C) < OXMXo(i-2, p7X_C) || OXMXo(i, p7X_C) < OXMXo(i-1, p7X_C))                             scur = p7T_C;
      else if (esl_FCompare(OXMXo(i, p7X_C), OXMXo(i-3, p7X_C) + om_fs->xf[p7O_C][p7O_LOOP], r_tol, a_tol) == eslOK)  scur = p7T_C;
      else if (esl_FCompare(OXMXo(i, p7X_C), OXMXo(i,   p7X_E) + om_fs->xf[p7O_E][p7O_MOVE], r_tol, a_tol) == eslOK)  scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7T_E:
      if (OXMXo(i, p7X_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);

      if (p7_fs_oprofile_IsLocal(om_fs))
        {
          scur = p7T_M;
          for (k = M; k >= 1; k--) if (esl_FCompare(OXMXo(i, p7X_E), OMMo(i,k), r_tol, a_tol) == eslOK) break;
          if (k == 0) ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
        }
      else
        {
          if      (esl_FCompare(OXMXo(i, p7X_E), OMMo(i,M), r_tol, a_tol) == eslOK) { scur = p7T_M; k = M; }
          else if (esl_FCompare(OXMXo(i, p7X_E), ODMo(i,M), r_tol, a_tol) == eslOK) { scur = p7T_D; k = M; }
          else    ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
        }
      break;

    case p7T_M:
      {
        union { __m128 v; float p[4]; } u;
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

    case p7T_D:
      {
        union { __m128 v; float p[4]; } u;
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

    case p7T_I:
      {
        union { __m128 v; float p[4]; } u;
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

    case p7T_N:
      if (OXMXo(i, p7X_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = (i == 0) ? p7T_S : p7T_N;
      break;

    case p7T_B:
      if (OXMXo(i, p7X_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

      if      (esl_FCompare(OXMXo(i, p7X_B), OXMXo(i, p7X_N) + om_fs->xf[p7O_N][p7O_MOVE], r_tol, a_tol) == eslOK) scur = p7T_N;
      else if (esl_FCompare(OXMXo(i, p7X_B), OXMXo(i, p7X_J) + om_fs->xf[p7O_J][p7O_MOVE], r_tol, a_tol) == eslOK) scur = p7T_J;
      else ESL_EXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    case p7T_J:
      if (OXMXo(i, p7X_J) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible J reached at i=%d", i);

      if      (esl_FCompare(OXMXo(i, p7X_J), OXMXo(i-3, p7X_J) + om_fs->xf[p7O_J][p7O_LOOP], r_tol, a_tol) == eslOK) scur = p7T_J;
      else if (esl_FCompare(OXMXo(i, p7X_J), OXMXo(i,   p7X_E) + om_fs->xf[p7O_E][p7O_LOOP], r_tol, a_tol) == eslOK) scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "J at i=%d couldn't be traced", i);
      break;

    default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over sprv */

    if (scur == p7T_M) {
      union { __m128 v; float p[4]; } u;
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

    if ((scur == p7T_N || scur == p7T_C) && scur == sprv) i--;
    if ( scur == p7T_J                   && scur == sprv) i -= 3;

    prev_c = c;
    c = 0;
    sprv = scur;
  } /* end traceback */

#undef OMMo
#undef ODMo
#undef OIMo
#undef OXMXo

  tr->M = om_fs->M;
  tr->L = L;
  return p7_trace_fs_Reverse(tr);
}
/*------------------ end, p7_Viterbi_Frameshift_Trace_sse() --------------------*/

#endif /* eslENABLE_SSE */
