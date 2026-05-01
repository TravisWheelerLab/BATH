/* SSE implementation of frameshift-aware Forward and Backward algorithms.
 * Ported from impl_sse/fwdback_fs.c for the impl_avx runtime-dispatch model.
 * All public functions carry the _sse suffix; the non-suffixed dispatchers
 * live in fwdback_fs.c.
 */
#include "p7_config.h"

#ifdef eslENABLE_SSE

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>   /* SSE  */
#include <emmintrin.h>   /* SSE2 */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_sse.h"

#include "hmmer.h"
#include "impl_avx.h"

/* IVX intermediate matrix access: [slot 0..p7P_xCODONS-1][stripe 0..Q-1] */
#define IVX(slot, q) (ivx[(slot)][(q)])

/*****************************************************************
 * 1. Forward/Backward
 *****************************************************************/

/* Function:  p7_ForwardParser_Frameshift_3Codons_sse()
 * Synopsis:  SSE frameshift-aware Forward, 3 codon lengths, linear memory.
 */
int
p7_ForwardParser_Frameshift_3Codons_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_OMX *ox, P7_OIVX *ov, float *opt_sc)
{
  register __m128 mpv2, dpv2, ipv2;   /* right-shifted prev2 row MDI for BM/MM/IM/DM   */
  register __m128 sv;                 /* temporary IVX accumulator                      */
  register __m128 msv;                /* M state accumulator for current position       */
  register __m128 dcv;                /* delayed D(i,q+1) = M(i,q)*tMD carry value     */
  register __m128 xEv;                /* E state sum: horizontal reduce to scalar xE    */
  register __m128 xBv2;               /* splatted B(i-2) for B->Mk transitions          */
  __m128   zerov;                     /* splatted 0.0                                   */
  float    xN, xE, xB, xC, xJ;       /* special state values at position i             */
  float    xN_buf[PARSER_ROWS_FWD];
  float    xB_buf[PARSER_ROWS_FWD];
  float    xJ_buf[PARSER_ROWS_FWD];
  float    xC_buf[PARSER_ROWS_FWD];
  int      b, b2, b3;
  int      curr, prev2, prev3;
  int      ivx_2, ivx_3, ivx_4;
  __m128  *dpc, *dpp2, *dpp3;
  __m128  *tp;
  __m128 **ivx = ov->ivx;
  int      Q = p7O_NQF(om_fs->M);
  int      i, q, j, r;
  int      c2, c3, c4;
  int      u, v, w, x;

  if (om_fs->codon_lengths != 3) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 3 codon lengths");

#if eslDEBUGLEVEL > 0
  if (om_fs->M > ox->allocQ4*4)        ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few columns)");
  if (ox->validR < PARSER_ROWS_FWD)    ESL_EXCEPTION(eslEINVAL, "DP matrix needs at least PARSER_ROWS_FWD MDI rows");
  if (L >= ox->allocXR)                ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few X rows)");
  if (! p7_fs_oprofile_IsLocal(om_fs)) ESL_EXCEPTION(eslEINVAL, "Forward implementation assumes local alignment mode");
#endif

  ox->M  = om_fs->M;
  ox->L  = L;
  ox->has_own_scales = TRUE;
  ox->totscale       = 0.0;
  zerov = _mm_setzero_ps();

  for (r = 0; r < PARSER_ROWS_FWD; r++)
    for (q = 0; q < Q; q++)
      MMO(ox->dpf[r],q) = IMO(ox->dpf[r],q) = DMO(ox->dpf[r],q) = zerov;

  for (r = 0; r < p7P_3CODONS; r++)
    for (q = 0; q < Q; q++)
      IVX(r, q) = zerov;

  for (r = 0; r < PARSER_ROWS_FWD; r++)
    xN_buf[r] = xB_buf[r] = xJ_buf[r] = xC_buf[r] = 0.0f;
  xN_buf[0] = xN_buf[1] = 1.0f;
  xB_buf[0] = xB_buf[1] = om_fs->xf[p7O_N][p7O_MOVE];

  ox->xmx[p7X_SCALE] = 1.0;
  ox->xmx[p7X_E]     = 0.0;
  ox->xmx[p7X_N]     = 1.0;
  ox->xmx[p7X_J]     = 0.0;
  ox->xmx[p7X_B]     = om_fs->xf[p7O_N][p7O_MOVE];
  ox->xmx[p7X_C]     = 0.0;
  ox->xmx[1*p7X_NXCELLS+p7X_SCALE] = 1.0;
  ox->xmx[1*p7X_NXCELLS+p7X_E]     = 0.0;
  ox->xmx[1*p7X_NXCELLS+p7X_N]     = 1.0;
  ox->xmx[1*p7X_NXCELLS+p7X_J]     = 0.0;
  ox->xmx[1*p7X_NXCELLS+p7X_B]     = om_fs->xf[p7O_N][p7O_MOVE];
  ox->xmx[1*p7X_NXCELLS+p7X_C]     = 0.0;

  u = v = p7P_MAXCODONS3;
  if (dsq[1] < p7P_MAXNUC) w = dsq[1]; else w = p7P_MAXCODONS3;
  if (dsq[2] < p7P_MAXNUC) x = dsq[2]; else x = p7P_MAXCODONS3;

  /* Initialization: i=2 */
  i = 2;
  c2 = p7P_CODON2_FS3(w, x); c2 = p7P_MINIDX(c2, p7P_DEGEN3_QC1);

  curr  = 2;  prev2 = 0;  prev3 = 3;
  ivx_2 = 2;  ivx_3 = 1;  ivx_4 = 0;

  dpc  = ox->dpf[curr];
  dpp2 = ox->dpf[prev2];
  dpp3 = ox->dpf[prev3];

  xBv2 = _mm_set1_ps(xB_buf[0]);
  tp   = om_fs->tfv;
  dcv  = zerov;
  xEv  = zerov;
  mpv2 = dpv2 = ipv2 = zerov;

  for (q = 0; q < Q; q++)
    {
      sv  =                _mm_mul_ps(xBv2, *tp); tp++;
      sv  = _mm_add_ps(sv, _mm_mul_ps(mpv2, *tp)); tp++;
      sv  = _mm_add_ps(sv, _mm_mul_ps(ipv2, *tp)); tp++;
      sv  = _mm_add_ps(sv, _mm_mul_ps(dpv2, *tp)); tp++;
      IVX(ivx_2, q) = sv;

      msv = _mm_mul_ps(sv, om_fs->rfv[c2][q]);
      xEv = _mm_add_ps(xEv, msv);

      MMO(dpc, q) = msv;
      DMO(dpc, q) = dcv;

      dcv = _mm_mul_ps(msv, *tp); tp++;

      IMO(dpc, q) = zerov;
      tp += 2;
    }

  dcv        = esl_sse_rightshiftz_float(dcv);
  DMO(dpc,0) = zerov;
  tp         = om_fs->tfv + 7*Q;
  for (q = 0; q < Q; q++)
    {
      DMO(dpc, q) = _mm_add_ps(dcv, DMO(dpc, q));
      dcv         = _mm_mul_ps(DMO(dpc, q), *tp); tp++;
    }
  if (om_fs->M < 100)
    {
      for (j = 1; j < 4; j++)
        {
          dcv = esl_sse_rightshiftz_float(dcv);
          tp  = om_fs->tfv + 7*Q;
          for (q = 0; q < Q; q++)
            { DMO(dpc, q) = _mm_add_ps(dcv, DMO(dpc, q)); dcv = _mm_mul_ps(dcv, *tp); tp++; }
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
              sv          = _mm_add_ps(dcv, DMO(dpc, q));
              cv          = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO(dpc, q)));
              DMO(dpc, q) = sv;
              dcv         = _mm_mul_ps(dcv, *tp); tp++;
            }
          if (! _mm_movemask_ps(cv)) break;
        }
    }

  for (q = 0; q < Q; q++) xEv = _mm_add_ps(DMO(dpc, q), xEv);
  xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(0, 3, 2, 1)));
  xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(1, 0, 3, 2)));
  _mm_store_ss(&xE, xEv);

  xN = 1.0f;
  xJ = xE * om_fs->xf[p7O_E][p7O_LOOP];
  xC = xE * om_fs->xf[p7O_E][p7O_MOVE];
  xB = xN * om_fs->xf[p7O_N][p7O_MOVE] + xJ * om_fs->xf[p7O_J][p7O_MOVE];

  if (xE > 1.0e4f)
    {
      float scale_factor = 1.0f / xE;
      xN *= scale_factor; xJ *= scale_factor; xC *= scale_factor; xB *= scale_factor;
      xEv = _mm_set1_ps(scale_factor);
      for (r = 0; r < PARSER_ROWS_FWD; r++)
        for (q = 0; q < Q; q++)
          { MMO(ox->dpf[r],q) = _mm_mul_ps(MMO(ox->dpf[r],q), xEv);
            DMO(ox->dpf[r],q) = _mm_mul_ps(DMO(ox->dpf[r],q), xEv);
            IMO(ox->dpf[r],q) = _mm_mul_ps(IMO(ox->dpf[r],q), xEv); }
      for (r = 0; r < p7P_3CODONS; r++)
        for (q = 0; q < Q; q++)
          IVX(r, q) = _mm_mul_ps(IVX(r, q), xEv);
      for (r = 0; r < PARSER_ROWS_FWD; r++)
        { xN_buf[r] *= scale_factor; xB_buf[r] *= scale_factor;
          xJ_buf[r] *= scale_factor; xC_buf[r] *= scale_factor; }
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

  /* Main recurrence: i = 3..L */
  for (i = 3; i <= L; i++)
    {
      u = v; v = w; w = x;
      if (dsq[i] < p7P_MAXNUC) x = dsq[i]; else x = p7P_MAXCODONS3;

      c2 = p7P_CODON2_FS3(w, x); c2 = p7P_MINIDX(c2, p7P_DEGEN3_QC1);
      c3 = p7P_CODON3_FS3(v, w, x); c3 = p7P_MINIDX(c3, p7P_DEGEN3_C);
      c4 = p7P_CODON4_FS3(u, v, w, x); c4 = p7P_MINIDX(c4, p7P_DEGEN3_QC1);

      curr  =     i               % PARSER_ROWS_FWD;
      prev2 = ((i-2) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD;
      prev3 = ((i-3) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD;

      ivx_2 =     i               % p7P_3CODONS;
      ivx_3 = ((i-1) % p7P_3CODONS + p7P_3CODONS) % p7P_3CODONS;
      ivx_4 = ((i-2) % p7P_3CODONS + p7P_3CODONS) % p7P_3CODONS;

      b  =     i               % PARSER_ROWS_FWD;
      b2 = ((i-2) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD;
      b3 = ((i-3) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD;

      dpc  = ox->dpf[curr];
      dpp2 = ox->dpf[prev2];
      dpp3 = ox->dpf[prev3];

      mpv2 = esl_sse_rightshiftz_float(MMO(dpp2, Q-1));
      dpv2 = esl_sse_rightshiftz_float(DMO(dpp2, Q-1));
      ipv2 = esl_sse_rightshiftz_float(IMO(dpp2, Q-1));

      xBv2 = _mm_set1_ps(xB_buf[b2]);
      tp   = om_fs->tfv;
      dcv  = zerov;
      xEv  = zerov;

      for (q = 0; q < Q; q++)
        {
          sv  =                _mm_mul_ps(xBv2, *tp); tp++;
          sv  = _mm_add_ps(sv, _mm_mul_ps(mpv2, *tp)); tp++;
          sv  = _mm_add_ps(sv, _mm_mul_ps(ipv2, *tp)); tp++;
          sv  = _mm_add_ps(sv, _mm_mul_ps(dpv2, *tp)); tp++;
          IVX(ivx_2, q) = sv;

          msv =                _mm_mul_ps(sv,              om_fs->rfv[c2][q]);
          msv = _mm_add_ps(msv, _mm_mul_ps(IVX(ivx_3, q), om_fs->rfv[c3][q]));
          msv = _mm_add_ps(msv, _mm_mul_ps(IVX(ivx_4, q), om_fs->rfv[c4][q]));
          xEv = _mm_add_ps(xEv, msv);

          mpv2 = MMO(dpp2, q);
          dpv2 = DMO(dpp2, q);
          ipv2 = IMO(dpp2, q);

          MMO(dpc, q) = msv;
          DMO(dpc, q) = dcv;

          dcv = _mm_mul_ps(msv, *tp); tp++;

          sv          =                _mm_mul_ps(MMO(dpp3, q), *tp); tp++;
          IMO(dpc, q) = _mm_add_ps(sv, _mm_mul_ps(IMO(dpp3, q), *tp)); tp++;
        }

      dcv        = esl_sse_rightshiftz_float(dcv);
      DMO(dpc,0) = zerov;
      tp         = om_fs->tfv + 7*Q;
      for (q = 0; q < Q; q++)
        { DMO(dpc, q) = _mm_add_ps(dcv, DMO(dpc, q)); dcv = _mm_mul_ps(DMO(dpc, q), *tp); tp++; }
      if (om_fs->M < 100)
        {
          for (j = 1; j < 4; j++)
            {
              dcv = esl_sse_rightshiftz_float(dcv);
              tp  = om_fs->tfv + 7*Q;
              for (q = 0; q < Q; q++)
                { DMO(dpc, q) = _mm_add_ps(dcv, DMO(dpc, q)); dcv = _mm_mul_ps(dcv, *tp); tp++; }
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
                  sv          = _mm_add_ps(dcv, DMO(dpc, q));
                  cv          = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO(dpc, q)));
                  DMO(dpc, q) = sv;
                  dcv         = _mm_mul_ps(dcv, *tp); tp++;
                }
              if (! _mm_movemask_ps(cv)) break;
            }
        }

      for (q = 0; q < Q; q++) xEv = _mm_add_ps(DMO(dpc, q), xEv);
      xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(0, 3, 2, 1)));
      xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(1, 0, 3, 2)));
      _mm_store_ss(&xE, xEv);

      xN = xN_buf[b3] * om_fs->xf[p7O_N][p7O_LOOP];
      xJ = xJ_buf[b3] * om_fs->xf[p7O_J][p7O_LOOP] + xE * om_fs->xf[p7O_E][p7O_LOOP];
      xC = xC_buf[b3] * om_fs->xf[p7O_C][p7O_LOOP] + xE * om_fs->xf[p7O_E][p7O_MOVE];
      xB = xN         * om_fs->xf[p7O_N][p7O_MOVE]  + xJ * om_fs->xf[p7O_J][p7O_MOVE];

      if (xE > 1.0e4f)
        {
          float scale_factor = 1.0f / xE;
          xN *= scale_factor; xJ *= scale_factor; xC *= scale_factor; xB *= scale_factor;
          xEv = _mm_set1_ps(scale_factor);
          for (r = 0; r < PARSER_ROWS_FWD; r++)
            for (q = 0; q < Q; q++)
              { MMO(ox->dpf[r],q) = _mm_mul_ps(MMO(ox->dpf[r],q), xEv);
                DMO(ox->dpf[r],q) = _mm_mul_ps(DMO(ox->dpf[r],q), xEv);
                IMO(ox->dpf[r],q) = _mm_mul_ps(IMO(ox->dpf[r],q), xEv); }
          for (r = 0; r < p7P_3CODONS; r++)
            for (q = 0; q < Q; q++)
              IVX(r, q) = _mm_mul_ps(IVX(r, q), xEv);
          for (r = 0; r < PARSER_ROWS_FWD; r++)
            { xN_buf[r] *= scale_factor; xB_buf[r] *= scale_factor;
              xJ_buf[r] *= scale_factor; xC_buf[r] *= scale_factor; }
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
    }

  /* Final score */
  {
    float xCL   = xC_buf[   L    % PARSER_ROWS_FWD];
    float xCLm1 = xC_buf[((L-1) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD];
    float xCLm2 = xC_buf[((L-2) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD];
    float xCtot = xCL
                + xCLm1 * om_fs->xf[p7O_C][p7O_LOOP]
                + xCLm2 * om_fs->xf[p7O_C][p7O_LOOP];

    if      (isnan(xCtot))           ESL_EXCEPTION(eslERANGE, "forward score is NaN");
    else if (isinf(xCtot) == 1)      ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");
    else if (L > 2 && xCtot == 0.0f) {
      if (opt_sc != NULL) *opt_sc = -eslINFINITY;
      return eslERANGE;
    }
    if (opt_sc != NULL)
      *opt_sc = ox->totscale + logf(xCtot * om_fs->xf[p7O_C][p7O_MOVE]);
  }

  return eslOK;
}


/* Function:  p7_BackwardParser_Frameshift_3Codons_sse()
 * Synopsis:  SSE frameshift-aware Backward, 3 codon lengths, linear memory.
 */
int
p7_BackwardParser_Frameshift_3Codons_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
  register __m128 sv;
  register __m128 dcv;
  register __m128 ivx_carry;
  __m128   tmmv, timv, tdmv;
  register __m128 xBv;
  register __m128 xEv;
  __m128   zerov;
  float    xN, xE, xB, xC, xJ;
  float    xN_buf[PARSER_ROWS_BWD];
  float    xB_buf[PARSER_ROWS_BWD];
  float    xJ_buf[PARSER_ROWS_BWD];
  float    xC_buf[PARSER_ROWS_BWD];
  int      b, b3;
  int      curr, prev2, prev3, prev4;
  __m128  *dpc, *dpp2, *dpp3, *dpp4;
  __m128  *tp7, *tp_dd, *tp_mm;
  __m128  *ivxf;
  int      Q = p7O_NQF(om_fs->M);
  int      i, q, j, r;
  int      c2, c3, c4;
  int      u, v, w, x;
  float    scale;

  ivxf = ov->ivx[0];

  if (om_fs->codon_lengths != 3)
    ESL_EXCEPTION(eslEINVAL, "profile not allocated for 3 codon lengths");
#if eslDEBUGLEVEL > 0
  if (om_fs->M > bck->allocQ4*4)
    ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few columns)");
  if (bck->validR < PARSER_ROWS_BWD)
    ESL_EXCEPTION(eslEINVAL, "DP matrix needs at least PARSER_ROWS_BWD MDI rows");
  if (L >= bck->allocXR)
    ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few X rows)");
  if (! p7_fs_oprofile_IsLocal(om_fs))
    ESL_EXCEPTION(eslEINVAL, "Backward implementation assumes local alignment mode");
#endif

  bck->M              = om_fs->M;
  bck->L              = L;
  bck->has_own_scales = FALSE;
  bck->totscale       = 0.0;
  zerov = _mm_setzero_ps();

  for (r = 0; r < PARSER_ROWS_BWD; r++)
    for (q = 0; q < Q; q++)
      MMO(bck->dpf[r],q) = IMO(bck->dpf[r],q) = DMO(bck->dpf[r],q) = zerov;

  for (r = 0; r < PARSER_ROWS_BWD; r++)
    xN_buf[r] = xB_buf[r] = xJ_buf[r] = xC_buf[r] = 0.0f;

  /* Initialization: rows L and L-1 */
  for (i = L; i >= L-1; i--)
    {
      b    = i % PARSER_ROWS_BWD;
      curr = b;

      if      (i == L)   xC = om_fs->xf[p7O_C][p7O_MOVE];
      else if (i == L-1) xC = om_fs->xf[p7O_C][p7O_LOOP] * om_fs->xf[p7O_C][p7O_MOVE];

      xN = xB = xJ = 0.0f;
      xE = xC * om_fs->xf[p7O_E][p7O_MOVE];
      xEv = _mm_set1_ps(xE);

      dpc = bck->dpf[curr];
      for (q = 0; q < Q; q++)
        { MMO(dpc,q) = xEv; DMO(dpc,q) = xEv; IMO(dpc,q) = zerov; }

      dcv    = esl_sse_leftshiftz_float(DMO(dpc,0));
      tp_dd  = om_fs->tfv + 8*Q - 1;
      for (q = Q-1; q >= 0; q--)
        { sv = _mm_mul_ps(dcv, *tp_dd); tp_dd--;
          DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), sv); dcv = DMO(dpc,q); }
      dcv = sv;
      for (j = 1; j < 4; j++)
        { dcv   = esl_sse_leftshiftz_float(dcv);
          tp_dd = om_fs->tfv + 8*Q - 1;
          for (q = Q-1; q >= 0; q--)
            { sv = _mm_mul_ps(dcv, *tp_dd); tp_dd--;
              DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), sv);
              dcv = _mm_mul_ps(dcv, *(tp_dd+1)); } }

      dcv  = esl_sse_leftshiftz_float(DMO(dpc,0));
      tp7  = om_fs->tfv + 7*Q - 3;
      for (q = Q-1; q >= 0; q--)
        { MMO(dpc,q) = _mm_add_ps(MMO(dpc,q), _mm_mul_ps(dcv, *tp7)); tp7 -= 7;
          dcv = DMO(dpc,q); }

      scale = fwd->xmx[i*p7X_NXCELLS+p7X_SCALE];
      bck->xmx[i*p7X_NXCELLS+p7X_SCALE] = scale;
      if (scale > 1.0f)
        { float sf = 1.0f / scale;
          xN *= sf; xJ *= sf; xC *= sf; xB *= sf; xE *= sf;
          xEv = _mm_set1_ps(sf);
          for (r = 0; r < PARSER_ROWS_BWD; r++)
            for (q = 0; q < Q; q++)
              { MMO(bck->dpf[r],q) = _mm_mul_ps(MMO(bck->dpf[r],q), xEv);
                DMO(bck->dpf[r],q) = _mm_mul_ps(DMO(bck->dpf[r],q), xEv);
                IMO(bck->dpf[r],q) = _mm_mul_ps(IMO(bck->dpf[r],q), xEv); }
          bck->totscale += log(scale); }

      xN_buf[b] = xN; xB_buf[b] = xB; xJ_buf[b] = xJ; xC_buf[b] = xC;
      bck->xmx[i*p7X_NXCELLS+p7X_E] = xE;
      bck->xmx[i*p7X_NXCELLS+p7X_N] = xN;
      bck->xmx[i*p7X_NXCELLS+p7X_J] = xJ;
      bck->xmx[i*p7X_NXCELLS+p7X_B] = xB;
      bck->xmx[i*p7X_NXCELLS+p7X_C] = xC;
    }

  /* Initialization: row L-2 */
  u = v = p7P_MAXCODONS3;
  if (dsq[L]   < p7P_MAXNUC) w = dsq[L];   else w = p7P_MAXCODONS3;
  if (dsq[L-1] < p7P_MAXNUC) x = dsq[L-1]; else x = p7P_MAXCODONS3;

  i     = L-2;
  b     = i % PARSER_ROWS_BWD;
  curr  = i    % PARSER_ROWS_BWD;
  prev2 = (i+2) % PARSER_ROWS_BWD;

  c2 = p7P_CODON2_FS3(x, w); c2 = p7P_MINIDX(c2, p7P_DEGEN3_QC1);

  dpc  = bck->dpf[curr];
  dpp2 = bck->dpf[prev2];

  for (q = 0; q < Q; q++)
    ivxf[q] = _mm_mul_ps(MMO(dpp2,q), om_fs->rfv[c2][q]);

  xBv = zerov;
  tp7 = om_fs->tfv;
  for (q = 0; q < Q; q++) { xBv = _mm_add_ps(xBv, _mm_mul_ps(ivxf[q], *tp7)); tp7 += 7; }
  xBv = _mm_add_ps(xBv, _mm_shuffle_ps(xBv, xBv, _MM_SHUFFLE(0,3,2,1)));
  xBv = _mm_add_ps(xBv, _mm_shuffle_ps(xBv, xBv, _MM_SHUFFLE(1,0,3,2)));
  _mm_store_ss(&xB, xBv);

  xC = om_fs->xf[p7O_C][p7O_LOOP] * om_fs->xf[p7O_C][p7O_MOVE];
  xJ = xB * om_fs->xf[p7O_J][p7O_MOVE];
  xN = xB * om_fs->xf[p7O_N][p7O_MOVE];
  xE = xJ * om_fs->xf[p7O_E][p7O_LOOP] + xC * om_fs->xf[p7O_E][p7O_MOVE];
  xEv = _mm_set1_ps(xE);

  for (q = 0; q < Q; q++) { MMO(dpc,q) = xEv; DMO(dpc,q) = xEv; IMO(dpc,q) = zerov; }

  ivx_carry = esl_sse_leftshiftz_float(ivxf[0]);
  tmmv = esl_sse_leftshiftz_float(*(om_fs->tfv+1));
  timv = esl_sse_leftshiftz_float(*(om_fs->tfv+2));
  tdmv = esl_sse_leftshiftz_float(*(om_fs->tfv+3));
  tp_mm = om_fs->tfv;
  for (q = Q-1; q >= 0; q--)
    { MMO(dpc,q) = _mm_add_ps(MMO(dpc,q), _mm_mul_ps(ivx_carry, tmmv));
      IMO(dpc,q) = _mm_add_ps(IMO(dpc,q), _mm_mul_ps(ivx_carry, timv));
      DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), _mm_mul_ps(ivx_carry, tdmv));
      ivx_carry  = ivxf[q];
      tp_mm = (tp_mm == om_fs->tfv) ? om_fs->tfv + 7*(Q-1) : tp_mm - 7;
      tmmv = *(tp_mm+1); timv = *(tp_mm+2); tdmv = *(tp_mm+3); }

  dcv   = esl_sse_leftshiftz_float(DMO(dpc,0));
  tp_dd = om_fs->tfv + 8*Q - 1;
  for (q = Q-1; q >= 0; q--)
    { sv = _mm_mul_ps(dcv, *tp_dd); tp_dd--;
      DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), sv); dcv = DMO(dpc,q); }
  dcv = sv;
  for (j = 1; j < 4; j++)
    { dcv = esl_sse_leftshiftz_float(dcv); tp_dd = om_fs->tfv + 8*Q - 1;
      for (q = Q-1; q >= 0; q--)
        { sv = _mm_mul_ps(dcv, *tp_dd); tp_dd--;
          DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), sv);
          dcv = _mm_mul_ps(dcv, *(tp_dd+1)); } }
  dcv = esl_sse_leftshiftz_float(DMO(dpc,0));
  tp7 = om_fs->tfv + 7*Q - 3;
  for (q = Q-1; q >= 0; q--)
    { MMO(dpc,q) = _mm_add_ps(MMO(dpc,q), _mm_mul_ps(dcv, *tp7)); tp7 -= 7;
      dcv = DMO(dpc,q); }

  scale = fwd->xmx[i*p7X_NXCELLS+p7X_SCALE];
  bck->xmx[i*p7X_NXCELLS+p7X_SCALE] = scale;
  if (scale > 1.0f)
    { float sf = 1.0f / scale;
      xN *= sf; xJ *= sf; xC *= sf; xB *= sf; xE *= sf;
      xEv = _mm_set1_ps(sf);
      for (r = 0; r < PARSER_ROWS_BWD; r++)
        for (q = 0; q < Q; q++)
          { MMO(bck->dpf[r],q) = _mm_mul_ps(MMO(bck->dpf[r],q), xEv);
            DMO(bck->dpf[r],q) = _mm_mul_ps(DMO(bck->dpf[r],q), xEv);
            IMO(bck->dpf[r],q) = _mm_mul_ps(IMO(bck->dpf[r],q), xEv); }
      for (r = 0; r < PARSER_ROWS_BWD; r++)
        { xN_buf[r] *= sf; xB_buf[r] *= sf; xJ_buf[r] *= sf; xC_buf[r] *= sf; }
      bck->totscale += log(scale); }

  xN_buf[b] = xN; xB_buf[b] = xB; xJ_buf[b] = xJ; xC_buf[b] = xC;
  bck->xmx[i*p7X_NXCELLS+p7X_E] = xE;
  bck->xmx[i*p7X_NXCELLS+p7X_N] = xN;
  bck->xmx[i*p7X_NXCELLS+p7X_J] = xJ;
  bck->xmx[i*p7X_NXCELLS+p7X_B] = xB;
  bck->xmx[i*p7X_NXCELLS+p7X_C] = xC;

  /* Main recurrence: i = L-3 down to 1 */
  for (i = L-3; i >= 1; i--)
    {
      u = v; v = w; w = x;
      if (dsq[i+1] < p7P_MAXNUC) x = dsq[i+1]; else x = p7P_MAXCODONS3;

      c2 = p7P_CODON2_FS3(x, w); c2 = p7P_MINIDX(c2, p7P_DEGEN3_QC1);
      c3 = p7P_CODON3_FS3(x, w, v); c3 = p7P_MINIDX(c3, p7P_DEGEN3_C);
      c4 = p7P_CODON4_FS3(x, w, v, u); c4 = p7P_MINIDX(c4, p7P_DEGEN3_QC1);

      curr  =  i       % PARSER_ROWS_BWD;
      prev2 = (i+2)    % PARSER_ROWS_BWD;
      prev3 = (i+3)    % PARSER_ROWS_BWD;
      prev4 = (i+4)    % PARSER_ROWS_BWD;
      b     =  i       % PARSER_ROWS_BWD;
      b3    = (i+3)    % PARSER_ROWS_BWD;

      dpc  = bck->dpf[curr];
      dpp2 = bck->dpf[prev2];
      dpp3 = bck->dpf[prev3];
      dpp4 = bck->dpf[prev4];

      for (q = 0; q < Q; q++)
        ivxf[q] = _mm_add_ps(_mm_add_ps(
                    _mm_mul_ps(MMO(dpp2,q), om_fs->rfv[c2][q]),
                    _mm_mul_ps(MMO(dpp3,q), om_fs->rfv[c3][q])),
                    _mm_mul_ps(MMO(dpp4,q), om_fs->rfv[c4][q]));

      xBv = zerov;
      tp7 = om_fs->tfv;
      for (q = 0; q < Q; q++) { xBv = _mm_add_ps(xBv, _mm_mul_ps(ivxf[q], *tp7)); tp7 += 7; }
      xBv = _mm_add_ps(xBv, _mm_shuffle_ps(xBv, xBv, _MM_SHUFFLE(0,3,2,1)));
      xBv = _mm_add_ps(xBv, _mm_shuffle_ps(xBv, xBv, _MM_SHUFFLE(1,0,3,2)));
      _mm_store_ss(&xB, xBv);

      xC = xC_buf[b3] * om_fs->xf[p7O_C][p7O_LOOP];
      xJ = xJ_buf[b3] * om_fs->xf[p7O_J][p7O_LOOP] + xB * om_fs->xf[p7O_J][p7O_MOVE];
      xN = xN_buf[b3] * om_fs->xf[p7O_N][p7O_LOOP] + xB * om_fs->xf[p7O_N][p7O_MOVE];
      xE = xJ * om_fs->xf[p7O_E][p7O_LOOP] + xC * om_fs->xf[p7O_E][p7O_MOVE];
      xEv = _mm_set1_ps(xE);

      for (q = 0; q < Q; q++)
        { MMO(dpc,q) = xEv; DMO(dpc,q) = xEv; IMO(dpc,q) = zerov; }

      tp7 = om_fs->tfv;
      for (q = 0; q < Q; q++)
        { MMO(dpc,q) = _mm_add_ps(MMO(dpc,q), _mm_mul_ps(IMO(dpp3,q), *(tp7+5)));
          IMO(dpc,q) = _mm_add_ps(IMO(dpc,q), _mm_mul_ps(IMO(dpp3,q), *(tp7+6)));
          tp7 += 7; }

      ivx_carry = esl_sse_leftshiftz_float(ivxf[0]);
      tmmv = esl_sse_leftshiftz_float(*(om_fs->tfv+1));
      timv = esl_sse_leftshiftz_float(*(om_fs->tfv+2));
      tdmv = esl_sse_leftshiftz_float(*(om_fs->tfv+3));
      tp_mm = om_fs->tfv;
      for (q = Q-1; q >= 0; q--)
        { MMO(dpc,q) = _mm_add_ps(MMO(dpc,q), _mm_mul_ps(ivx_carry, tmmv));
          IMO(dpc,q) = _mm_add_ps(IMO(dpc,q), _mm_mul_ps(ivx_carry, timv));
          DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), _mm_mul_ps(ivx_carry, tdmv));
          ivx_carry  = ivxf[q];
          tp_mm = (tp_mm == om_fs->tfv) ? om_fs->tfv + 7*(Q-1) : tp_mm - 7;
          tmmv = *(tp_mm+1); timv = *(tp_mm+2); tdmv = *(tp_mm+3); }

      dcv   = esl_sse_leftshiftz_float(DMO(dpc,0));
      tp_dd = om_fs->tfv + 8*Q - 1;
      for (q = Q-1; q >= 0; q--)
        { sv = _mm_mul_ps(dcv, *tp_dd); tp_dd--;
          DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), sv); dcv = DMO(dpc,q); }
      dcv = sv;
      for (j = 1; j < 4; j++)
        { dcv = esl_sse_leftshiftz_float(dcv); tp_dd = om_fs->tfv + 8*Q - 1;
          for (q = Q-1; q >= 0; q--)
            { sv = _mm_mul_ps(dcv, *tp_dd); tp_dd--;
              DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), sv);
              dcv = _mm_mul_ps(dcv, *(tp_dd+1)); } }

      dcv = esl_sse_leftshiftz_float(DMO(dpc,0));
      tp7 = om_fs->tfv + 7*Q - 3;
      for (q = Q-1; q >= 0; q--)
        { MMO(dpc,q) = _mm_add_ps(MMO(dpc,q), _mm_mul_ps(dcv, *tp7)); tp7 -= 7;
          dcv = DMO(dpc,q); }

      if (xB > 1.0e16f) bck->has_own_scales = TRUE;

      if (bck->has_own_scales) scale = (xB > 1.0e4f) ? xB : 1.0f;
      else                     scale = fwd->xmx[i*p7X_NXCELLS+p7X_SCALE];

      bck->xmx[i*p7X_NXCELLS+p7X_SCALE] = scale;
      if (scale > 1.0f)
        { float sf = 1.0f / scale;
          xN *= sf; xJ *= sf; xC *= sf; xB *= sf; xE *= sf;
          xEv = _mm_set1_ps(sf);
          for (r = 0; r < PARSER_ROWS_BWD; r++)
            for (q = 0; q < Q; q++)
              { MMO(bck->dpf[r],q) = _mm_mul_ps(MMO(bck->dpf[r],q), xEv);
                DMO(bck->dpf[r],q) = _mm_mul_ps(DMO(bck->dpf[r],q), xEv);
                IMO(bck->dpf[r],q) = _mm_mul_ps(IMO(bck->dpf[r],q), xEv); }
          for (r = 0; r < PARSER_ROWS_BWD; r++)
            { xN_buf[r] *= sf; xB_buf[r] *= sf; xJ_buf[r] *= sf; xC_buf[r] *= sf; }
          bck->totscale += log(scale); }

      xN_buf[b] = xN; xB_buf[b] = xB; xJ_buf[b] = xJ; xC_buf[b] = xC;
      bck->xmx[i*p7X_NXCELLS+p7X_E] = xE;
      bck->xmx[i*p7X_NXCELLS+p7X_N] = xN;
      bck->xmx[i*p7X_NXCELLS+p7X_J] = xJ;
      bck->xmx[i*p7X_NXCELLS+p7X_B] = xB;
      bck->xmx[i*p7X_NXCELLS+p7X_C] = xC;
    }

  /* Termination: i=0 */
  u = v; v = w; w = x;
  if (dsq[1] < p7P_MAXNUC) x = dsq[1]; else x = p7P_MAXCODONS3;

  c2 = p7P_CODON2_FS3(x, w); c2 = p7P_MINIDX(c2, p7P_DEGEN3_QC1);
  c3 = p7P_CODON3_FS3(x, w, v); c3 = p7P_MINIDX(c3, p7P_DEGEN3_C);
  c4 = p7P_CODON4_FS3(x, w, v, u); c4 = p7P_MINIDX(c4, p7P_DEGEN3_QC1);

  prev2 = 2 % PARSER_ROWS_BWD;
  prev3 = 3 % PARSER_ROWS_BWD;
  prev4 = 4 % PARSER_ROWS_BWD;

  dpp2 = bck->dpf[prev2];
  dpp3 = bck->dpf[prev3];
  dpp4 = bck->dpf[prev4];

  for (q = 0; q < Q; q++)
    ivxf[q] = _mm_add_ps(_mm_add_ps(
                _mm_mul_ps(MMO(dpp2,q), om_fs->rfv[c2][q]),
                _mm_mul_ps(MMO(dpp3,q), om_fs->rfv[c3][q])),
                _mm_mul_ps(MMO(dpp4,q), om_fs->rfv[c4][q]));

  xBv = zerov;
  tp7 = om_fs->tfv;
  for (q = 0; q < Q; q++) { xBv = _mm_add_ps(xBv, _mm_mul_ps(ivxf[q], *tp7)); tp7 += 7; }
  xBv = _mm_add_ps(xBv, _mm_shuffle_ps(xBv, xBv, _MM_SHUFFLE(0,3,2,1)));
  xBv = _mm_add_ps(xBv, _mm_shuffle_ps(xBv, xBv, _MM_SHUFFLE(1,0,3,2)));
  _mm_store_ss(&xB, xBv);

  xN = xN_buf[3 % PARSER_ROWS_BWD] * om_fs->xf[p7O_N][p7O_LOOP]
     + xB                           * om_fs->xf[p7O_N][p7O_MOVE];

  bck->xmx[p7X_B]     = xB;
  bck->xmx[p7X_N]     = xN;
  bck->xmx[p7X_J]     = 0.0f;
  bck->xmx[p7X_C]     = 0.0f;
  bck->xmx[p7X_E]     = 0.0f;
  bck->xmx[p7X_SCALE] = 1.0f;

  {
    float xNtot = xN
                + xN_buf[1 % PARSER_ROWS_BWD]
                + xN_buf[2 % PARSER_ROWS_BWD];

    if      (isnan(xNtot))           ESL_EXCEPTION(eslERANGE, "backward score is NaN");
    else if (isinf(xNtot) == 1)      ESL_EXCEPTION(eslERANGE, "backward score overflow (is infinity)");
    else if (L > 0 && xNtot == 0.0f) {
      if (opt_sc != NULL) *opt_sc = -eslINFINITY;
      return eslERANGE;
    }
    if (opt_sc != NULL)
      *opt_sc = bck->totscale + logf(xNtot);
  }

  return eslOK;
}

/* Function:  p7_ForwardParser_Frameshift_5Codons_sse()
 * Synopsis:  SSE frameshift-aware Forward, 5 codon lengths, linear memory.
 */
int
p7_ForwardParser_Frameshift_5Codons_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_OMX *fwd, P7_OIVX *ov, float *opt_sc)
{
  register __m128 mpv1, dpv1, ipv1;
  register __m128 sv;
  register __m128 msv;
  register __m128 dcv;
  register __m128 xEv;
  register __m128 xBv1;
  __m128   zerov;
  float    xN, xE, xB, xC, xJ;
  float    xN_buf[PARSER_ROWS_FWD];
  float    xB_buf[PARSER_ROWS_FWD];
  float    xJ_buf[PARSER_ROWS_FWD];
  float    xC_buf[PARSER_ROWS_FWD];
  int      b, b1, b3;
  int      curr, prev1, prev3;
  int      ivx_1, ivx_2, ivx_3, ivx_4, ivx_5;
  __m128  *dpc, *dpp1, *dpp3;
  __m128  *tp;
  __m128 **ivx = ov->ivx;
  int      Q = p7O_NQF(om_fs->M);
  int      i, q, j, r;
  int      c1, c2, c3, c4, c5;
  int      t, u, v, w, x;

  if (om_fs->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");

#if eslDEBUGLEVEL > 0
  if (om_fs->M > fwd->allocQ4*4)        ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few columns)");
  if (fwd->validR < PARSER_ROWS_FWD)    ESL_EXCEPTION(eslEINVAL, "DP matrix needs at least PARSER_ROWS_FWD MDI rows");
  if (L >= fwd->allocXR)                ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few X rows)");
  if (! p7_fs_oprofile_IsLocal(om_fs)) ESL_EXCEPTION(eslEINVAL, "Forward implementation assumes local alignment mode");
#endif

  fwd->M              = om_fs->M;
  fwd->L              = L;
  fwd->has_own_scales = TRUE;
  fwd->totscale       = 0.0;
  zerov = _mm_setzero_ps();

  for (r = 0; r < PARSER_ROWS_FWD; r++)
    for (q = 0; q < Q; q++)
      MMO(fwd->dpf[r],q) = IMO(fwd->dpf[r],q) = DMO(fwd->dpf[r],q) = zerov;

  for (r = 0; r < p7P_5CODONS; r++)
    for (q = 0; q < Q; q++)
      IVX(r, q) = zerov;

  for (r = 0; r < PARSER_ROWS_FWD; r++)
    xN_buf[r] = xB_buf[r] = xJ_buf[r] = xC_buf[r] = 0.0f;
  xN_buf[0] = xN_buf[1] = xN_buf[2] = 1.0f;
  xB_buf[0] = xB_buf[1] = xB_buf[2] = om_fs->xf[p7O_N][p7O_MOVE];

  for (r = 0; r < 3; r++)
    {
      fwd->xmx[r*p7X_NXCELLS+p7X_SCALE] = 1.0f;
      fwd->xmx[r*p7X_NXCELLS+p7X_E]     = 0.0f;
      fwd->xmx[r*p7X_NXCELLS+p7X_N]     = 1.0f;
      fwd->xmx[r*p7X_NXCELLS+p7X_J]     = 0.0f;
      fwd->xmx[r*p7X_NXCELLS+p7X_B]     = om_fs->xf[p7O_N][p7O_MOVE];
      fwd->xmx[r*p7X_NXCELLS+p7X_C]     = 0.0f;
    }

  t = u = v = w = p7P_MAXCODONS5;
  if (dsq[1] < p7P_MAXNUC) x = dsq[1]; else x = p7P_MAXCODONS5;

  /* Initialization: i=1 (only 1-nt codon c1) */
  i = 1;
  c1 = p7P_CODON1_FS5(x); c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);

  curr  = 1;  prev1 = 0;  prev3 = 2;
  ivx_1 = 1;

  dpc  = fwd->dpf[curr];
  dpp1 = fwd->dpf[prev1];
  dpp3 = fwd->dpf[prev3];

  xBv1 = _mm_set1_ps(xB_buf[0]);
  tp   = om_fs->tfv;
  dcv  = zerov; xEv = zerov;
  mpv1 = dpv1 = ipv1 = zerov;

  for (q = 0; q < Q; q++)
    {
      sv  =                _mm_mul_ps(xBv1, *tp); tp++;
      sv  = _mm_add_ps(sv, _mm_mul_ps(mpv1, *tp)); tp++;
      sv  = _mm_add_ps(sv, _mm_mul_ps(ipv1, *tp)); tp++;
      sv  = _mm_add_ps(sv, _mm_mul_ps(dpv1, *tp)); tp++;
      IVX(ivx_1, q) = sv;
      msv = _mm_mul_ps(sv, om_fs->rfv[c1][q]);
      xEv = _mm_add_ps(xEv, msv);
      MMO(dpc, q) = msv;
      DMO(dpc, q) = dcv;
      dcv = _mm_mul_ps(msv, *tp); tp++;
      IMO(dpc, q) = zerov;
      tp += 2;
    }

  dcv = esl_sse_rightshiftz_float(dcv);
  DMO(dpc,0) = zerov;
  tp = om_fs->tfv + 7*Q;
  for (q = 0; q < Q; q++)
    { DMO(dpc, q) = _mm_add_ps(dcv, DMO(dpc, q)); dcv = _mm_mul_ps(DMO(dpc, q), *tp); tp++; }
  if (om_fs->M < 100)
    { for (j = 1; j < 4; j++) { dcv = esl_sse_rightshiftz_float(dcv); tp = om_fs->tfv + 7*Q;
        for (q = 0; q < Q; q++) { DMO(dpc,q) = _mm_add_ps(dcv,DMO(dpc,q)); dcv=_mm_mul_ps(dcv,*tp); tp++; } } }
  else
    { for (j = 1; j < 4; j++) { register __m128 cv; dcv=esl_sse_rightshiftz_float(dcv); tp=om_fs->tfv+7*Q; cv=zerov;
        for (q=0;q<Q;q++) { sv=_mm_add_ps(dcv,DMO(dpc,q)); cv=_mm_or_ps(cv,_mm_cmpgt_ps(sv,DMO(dpc,q))); DMO(dpc,q)=sv; dcv=_mm_mul_ps(dcv,*tp); tp++; }
        if (!_mm_movemask_ps(cv)) break; } }

  for (q = 0; q < Q; q++) xEv = _mm_add_ps(DMO(dpc,q), xEv);
  xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv,xEv,_MM_SHUFFLE(0,3,2,1)));
  xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv,xEv,_MM_SHUFFLE(1,0,3,2)));
  _mm_store_ss(&xE, xEv);

  xN = 1.0f;
  xJ = xE * om_fs->xf[p7O_E][p7O_LOOP];
  xC = xE * om_fs->xf[p7O_E][p7O_MOVE];
  xB = xN * om_fs->xf[p7O_N][p7O_MOVE] + xJ * om_fs->xf[p7O_J][p7O_MOVE];

  if (xE > 1.0e4f)
    {
      float scale_factor = 1.0f / xE;
      xN *= scale_factor; xJ *= scale_factor; xC *= scale_factor; xB *= scale_factor;
      xEv = _mm_set1_ps(scale_factor);
      for (r = 0; r < PARSER_ROWS_FWD; r++)
        for (q = 0; q < Q; q++)
          { MMO(fwd->dpf[r],q)=_mm_mul_ps(MMO(fwd->dpf[r],q),xEv);
            DMO(fwd->dpf[r],q)=_mm_mul_ps(DMO(fwd->dpf[r],q),xEv);
            IMO(fwd->dpf[r],q)=_mm_mul_ps(IMO(fwd->dpf[r],q),xEv); }
      for (r = 0; r < p7P_5CODONS; r++)
        for (q = 0; q < Q; q++) IVX(r,q)=_mm_mul_ps(IVX(r,q),xEv);
      for (r = 0; r < PARSER_ROWS_FWD; r++)
        { xN_buf[r]*=scale_factor; xB_buf[r]*=scale_factor; xJ_buf[r]*=scale_factor; xC_buf[r]*=scale_factor; }
      fwd->xmx[1*p7X_NXCELLS+p7X_SCALE] = xE;
      fwd->totscale += log(xE);
      xE = 1.0f;
    }
  else fwd->xmx[1*p7X_NXCELLS+p7X_SCALE] = 1.0f;

  xN_buf[1]=xN; xB_buf[1]=xB; xJ_buf[1]=xJ; xC_buf[1]=xC;
  fwd->xmx[1*p7X_NXCELLS+p7X_E]=xE; fwd->xmx[1*p7X_NXCELLS+p7X_N]=xN;
  fwd->xmx[1*p7X_NXCELLS+p7X_J]=xJ; fwd->xmx[1*p7X_NXCELLS+p7X_B]=xB;
  fwd->xmx[1*p7X_NXCELLS+p7X_C]=xC;

  /* Initialization: i=2 (1-nt and 2-nt codons) */
  i = 2;
  t = u = v = p7P_MAXCODONS5; w = x;
  if (dsq[2] < p7P_MAXNUC) x = dsq[2]; else x = p7P_MAXCODONS5;

  c1 = p7P_CODON1_FS5(x);    c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
  c2 = p7P_CODON2_FS5(w, x); c2 = p7P_MINIDX(c2, p7P_DEGEN5_QC1);

  curr = 2; prev1 = 1; prev3 = 3;
  ivx_1 = 2; ivx_2 = 1;

  dpc  = fwd->dpf[curr];
  dpp1 = fwd->dpf[prev1];
  dpp3 = fwd->dpf[prev3];

  xBv1 = _mm_set1_ps(xB_buf[1]);
  tp   = om_fs->tfv; dcv = zerov; xEv = zerov;
  mpv1 = esl_sse_rightshiftz_float(MMO(dpp1, Q-1));
  dpv1 = esl_sse_rightshiftz_float(DMO(dpp1, Q-1));
  ipv1 = esl_sse_rightshiftz_float(IMO(dpp1, Q-1));

  for (q = 0; q < Q; q++)
    {
      sv  =                _mm_mul_ps(xBv1, *tp); tp++;
      sv  = _mm_add_ps(sv, _mm_mul_ps(mpv1, *tp)); tp++;
      sv  = _mm_add_ps(sv, _mm_mul_ps(ipv1, *tp)); tp++;
      sv  = _mm_add_ps(sv, _mm_mul_ps(dpv1, *tp)); tp++;
      IVX(ivx_1, q) = sv;
      msv =                _mm_mul_ps(sv,              om_fs->rfv[c1][q]);
      msv = _mm_add_ps(msv, _mm_mul_ps(IVX(ivx_2, q), om_fs->rfv[c2][q]));
      xEv = _mm_add_ps(xEv, msv);
      mpv1 = MMO(dpp1,q); dpv1 = DMO(dpp1,q); ipv1 = IMO(dpp1,q);
      MMO(dpc,q) = msv; DMO(dpc,q) = dcv;
      dcv = _mm_mul_ps(msv, *tp); tp++;
      IMO(dpc,q) = zerov; tp += 2;
    }

  dcv = esl_sse_rightshiftz_float(dcv);
  DMO(dpc,0) = zerov; tp = om_fs->tfv + 7*Q;
  for (q = 0; q < Q; q++)
    { DMO(dpc,q) = _mm_add_ps(dcv,DMO(dpc,q)); dcv = _mm_mul_ps(DMO(dpc,q),*tp); tp++; }
  if (om_fs->M < 100)
    { for (j=1;j<4;j++) { dcv=esl_sse_rightshiftz_float(dcv); tp=om_fs->tfv+7*Q;
        for (q=0;q<Q;q++) { DMO(dpc,q)=_mm_add_ps(dcv,DMO(dpc,q)); dcv=_mm_mul_ps(dcv,*tp); tp++; } } }
  else
    { for (j=1;j<4;j++) { register __m128 cv; dcv=esl_sse_rightshiftz_float(dcv); tp=om_fs->tfv+7*Q; cv=zerov;
        for (q=0;q<Q;q++) { sv=_mm_add_ps(dcv,DMO(dpc,q)); cv=_mm_or_ps(cv,_mm_cmpgt_ps(sv,DMO(dpc,q))); DMO(dpc,q)=sv; dcv=_mm_mul_ps(dcv,*tp); tp++; }
        if (!_mm_movemask_ps(cv)) break; } }

  for (q = 0; q < Q; q++) xEv = _mm_add_ps(DMO(dpc,q), xEv);
  xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv,xEv,_MM_SHUFFLE(0,3,2,1)));
  xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv,xEv,_MM_SHUFFLE(1,0,3,2)));
  _mm_store_ss(&xE, xEv);

  xN = 1.0f;
  xJ = xE * om_fs->xf[p7O_E][p7O_LOOP];
  xC = xE * om_fs->xf[p7O_E][p7O_MOVE];
  xB = xN * om_fs->xf[p7O_N][p7O_MOVE] + xJ * om_fs->xf[p7O_J][p7O_MOVE];

  if (xE > 1.0e4f)
    {
      float scale_factor = 1.0f / xE;
      xN*=scale_factor; xJ*=scale_factor; xC*=scale_factor; xB*=scale_factor;
      xEv = _mm_set1_ps(scale_factor);
      for (r=0;r<PARSER_ROWS_FWD;r++) for (q=0;q<Q;q++)
        { MMO(fwd->dpf[r],q)=_mm_mul_ps(MMO(fwd->dpf[r],q),xEv);
          DMO(fwd->dpf[r],q)=_mm_mul_ps(DMO(fwd->dpf[r],q),xEv);
          IMO(fwd->dpf[r],q)=_mm_mul_ps(IMO(fwd->dpf[r],q),xEv); }
      for (r=0;r<p7P_5CODONS;r++) for (q=0;q<Q;q++) IVX(r,q)=_mm_mul_ps(IVX(r,q),xEv);
      for (r=0;r<PARSER_ROWS_FWD;r++)
        { xN_buf[r]*=scale_factor; xB_buf[r]*=scale_factor; xJ_buf[r]*=scale_factor; xC_buf[r]*=scale_factor; }
      fwd->xmx[2*p7X_NXCELLS+p7X_SCALE]=xE; fwd->totscale+=log(xE); xE=1.0f;
    }
  else fwd->xmx[2*p7X_NXCELLS+p7X_SCALE] = 1.0f;

  xN_buf[2]=xN; xB_buf[2]=xB; xJ_buf[2]=xJ; xC_buf[2]=xC;
  fwd->xmx[2*p7X_NXCELLS+p7X_E]=xE; fwd->xmx[2*p7X_NXCELLS+p7X_N]=xN;
  fwd->xmx[2*p7X_NXCELLS+p7X_J]=xJ; fwd->xmx[2*p7X_NXCELLS+p7X_B]=xB;
  fwd->xmx[2*p7X_NXCELLS+p7X_C]=xC;

  /* Main recurrence: i = 3..L */
  for (i = 3; i <= L; i++)
    {
      t=u; u=v; v=w; w=x;
      if (dsq[i] < p7P_MAXNUC) x = dsq[i]; else x = p7P_MAXCODONS5;

      c1 = p7P_CODON1_FS5(x);           c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
      c2 = p7P_CODON2_FS5(w, x);        c2 = p7P_MINIDX(c2, p7P_DEGEN5_QC1);
      c3 = p7P_CODON3_FS5(v, w, x);     c3 = p7P_MINIDX(c3, p7P_DEGEN5_C);
      c4 = p7P_CODON4_FS5(u, v, w, x);  c4 = p7P_MINIDX(c4, p7P_DEGEN5_QC1);
      c5 = p7P_CODON5_FS5(t, u, v, w, x); c5 = p7P_MINIDX(c5, p7P_DEGEN5_QC2);

      curr  =     i               % PARSER_ROWS_FWD;
      prev1 = ((i-1) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD;
      prev3 = ((i-3) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD;

      ivx_1 =     i               % p7P_5CODONS;
      ivx_2 = ((i-1) % p7P_5CODONS + p7P_5CODONS) % p7P_5CODONS;
      ivx_3 = ((i-2) % p7P_5CODONS + p7P_5CODONS) % p7P_5CODONS;
      ivx_4 = ((i-3) % p7P_5CODONS + p7P_5CODONS) % p7P_5CODONS;
      ivx_5 = ((i-4) % p7P_5CODONS + p7P_5CODONS) % p7P_5CODONS;

      b  =     i               % PARSER_ROWS_FWD;
      b1 = ((i-1) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD;
      b3 = ((i-3) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD;

      dpc  = fwd->dpf[curr];
      dpp1 = fwd->dpf[prev1];
      dpp3 = fwd->dpf[prev3];

      mpv1 = esl_sse_rightshiftz_float(MMO(dpp1, Q-1));
      dpv1 = esl_sse_rightshiftz_float(DMO(dpp1, Q-1));
      ipv1 = esl_sse_rightshiftz_float(IMO(dpp1, Q-1));

      xBv1 = _mm_set1_ps(xB_buf[b1]);
      tp = om_fs->tfv; dcv = zerov; xEv = zerov;

      for (q = 0; q < Q; q++)
        {
          sv  =                _mm_mul_ps(xBv1, *tp); tp++;
          sv  = _mm_add_ps(sv, _mm_mul_ps(mpv1, *tp)); tp++;
          sv  = _mm_add_ps(sv, _mm_mul_ps(ipv1, *tp)); tp++;
          sv  = _mm_add_ps(sv, _mm_mul_ps(dpv1, *tp)); tp++;
          IVX(ivx_1, q) = sv;
          msv =                _mm_mul_ps(sv,               om_fs->rfv[c1][q]);
          msv = _mm_add_ps(msv, _mm_mul_ps(IVX(ivx_2, q),  om_fs->rfv[c2][q]));
          msv = _mm_add_ps(msv, _mm_mul_ps(IVX(ivx_3, q),  om_fs->rfv[c3][q]));
          msv = _mm_add_ps(msv, _mm_mul_ps(IVX(ivx_4, q),  om_fs->rfv[c4][q]));
          msv = _mm_add_ps(msv, _mm_mul_ps(IVX(ivx_5, q),  om_fs->rfv[c5][q]));
          xEv = _mm_add_ps(xEv, msv);
          mpv1 = MMO(dpp1,q); dpv1 = DMO(dpp1,q); ipv1 = IMO(dpp1,q);
          MMO(dpc,q) = msv; DMO(dpc,q) = dcv;
          dcv = _mm_mul_ps(msv, *tp); tp++;
          sv =                _mm_mul_ps(MMO(dpp3,q), *tp); tp++;
          IMO(dpc,q) = _mm_add_ps(sv, _mm_mul_ps(IMO(dpp3,q), *tp)); tp++;
        }

      dcv=esl_sse_rightshiftz_float(dcv); DMO(dpc,0)=zerov; tp=om_fs->tfv+7*Q;
      for (q=0;q<Q;q++) { DMO(dpc,q)=_mm_add_ps(dcv,DMO(dpc,q)); dcv=_mm_mul_ps(DMO(dpc,q),*tp); tp++; }
      if (om_fs->M < 100)
        { for (j=1;j<4;j++) { dcv=esl_sse_rightshiftz_float(dcv); tp=om_fs->tfv+7*Q;
            for (q=0;q<Q;q++) { DMO(dpc,q)=_mm_add_ps(dcv,DMO(dpc,q)); dcv=_mm_mul_ps(dcv,*tp); tp++; } } }
      else
        { for (j=1;j<4;j++) { register __m128 cv; dcv=esl_sse_rightshiftz_float(dcv); tp=om_fs->tfv+7*Q; cv=zerov;
            for (q=0;q<Q;q++) { sv=_mm_add_ps(dcv,DMO(dpc,q)); cv=_mm_or_ps(cv,_mm_cmpgt_ps(sv,DMO(dpc,q))); DMO(dpc,q)=sv; dcv=_mm_mul_ps(dcv,*tp); tp++; }
            if (!_mm_movemask_ps(cv)) break; } }

      for (q=0;q<Q;q++) xEv=_mm_add_ps(DMO(dpc,q),xEv);
      xEv=_mm_add_ps(xEv,_mm_shuffle_ps(xEv,xEv,_MM_SHUFFLE(0,3,2,1)));
      xEv=_mm_add_ps(xEv,_mm_shuffle_ps(xEv,xEv,_MM_SHUFFLE(1,0,3,2)));
      _mm_store_ss(&xE, xEv);

      xN = xN_buf[b3] * om_fs->xf[p7O_N][p7O_LOOP];
      xJ = xJ_buf[b3] * om_fs->xf[p7O_J][p7O_LOOP] + xE * om_fs->xf[p7O_E][p7O_LOOP];
      xC = xC_buf[b3] * om_fs->xf[p7O_C][p7O_LOOP] + xE * om_fs->xf[p7O_E][p7O_MOVE];
      xB = xN         * om_fs->xf[p7O_N][p7O_MOVE]  + xJ * om_fs->xf[p7O_J][p7O_MOVE];

      if (xE > 1.0e4f)
        {
          float scale_factor = 1.0f / xE;
          xN*=scale_factor; xJ*=scale_factor; xC*=scale_factor; xB*=scale_factor;
          xEv=_mm_set1_ps(scale_factor);
          for (r=0;r<PARSER_ROWS_FWD;r++) for (q=0;q<Q;q++)
            { MMO(fwd->dpf[r],q)=_mm_mul_ps(MMO(fwd->dpf[r],q),xEv);
              DMO(fwd->dpf[r],q)=_mm_mul_ps(DMO(fwd->dpf[r],q),xEv);
              IMO(fwd->dpf[r],q)=_mm_mul_ps(IMO(fwd->dpf[r],q),xEv); }
          for (r=0;r<p7P_5CODONS;r++) for (q=0;q<Q;q++) IVX(r,q)=_mm_mul_ps(IVX(r,q),xEv);
          for (r=0;r<PARSER_ROWS_FWD;r++)
            { xN_buf[r]*=scale_factor; xB_buf[r]*=scale_factor; xJ_buf[r]*=scale_factor; xC_buf[r]*=scale_factor; }
          fwd->xmx[i*p7X_NXCELLS+p7X_SCALE]=xE; fwd->totscale+=log(xE); xE=1.0f;
        }
      else fwd->xmx[i*p7X_NXCELLS+p7X_SCALE] = 1.0f;

      xN_buf[b]=xN; xB_buf[b]=xB; xJ_buf[b]=xJ; xC_buf[b]=xC;
      fwd->xmx[i*p7X_NXCELLS+p7X_E]=xE; fwd->xmx[i*p7X_NXCELLS+p7X_N]=xN;
      fwd->xmx[i*p7X_NXCELLS+p7X_J]=xJ; fwd->xmx[i*p7X_NXCELLS+p7X_B]=xB;
      fwd->xmx[i*p7X_NXCELLS+p7X_C]=xC;
    }

  /* Final score */
  {
    float xCL   = xC_buf[   L    % PARSER_ROWS_FWD];
    float xCLm1 = xC_buf[((L-1) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD];
    float xCLm2 = xC_buf[((L-2) % PARSER_ROWS_FWD + PARSER_ROWS_FWD) % PARSER_ROWS_FWD];
    float xCtot = xCL + xCLm1 * om_fs->xf[p7O_C][p7O_LOOP] + xCLm2 * om_fs->xf[p7O_C][p7O_LOOP];

    if      (isnan(xCtot))           ESL_EXCEPTION(eslERANGE, "forward score is NaN");
    else if (isinf(xCtot) == 1)      ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");
    else if (L > 1 && xCtot == 0.0f) { if (opt_sc) *opt_sc = -eslINFINITY; return eslERANGE; }
    if (opt_sc) *opt_sc = fwd->totscale + logf(xCtot * om_fs->xf[p7O_C][p7O_MOVE]);
  }

  return eslOK;
}


/* Function:  p7_BackwardParser_Frameshift_5Codons_sse()
 * Synopsis:  SSE frameshift-aware Backward, 5 codon lengths, linear memory.
 */
int
p7_BackwardParser_Frameshift_5Codons_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
  register __m128 sv;
  register __m128 dcv;
  register __m128 ivx_carry;
  __m128   tmmv, timv, tdmv;
  register __m128 xBv;
  register __m128 xEv;
  __m128   zerov;
  float    xN, xE, xB, xC, xJ;
  float    xN_buf[PARSER_ROWS_BWD];
  float    xB_buf[PARSER_ROWS_BWD];
  float    xJ_buf[PARSER_ROWS_BWD];
  float    xC_buf[PARSER_ROWS_BWD];
  int      b, b3;
  int      curr, prev1, prev2, prev3, prev4, prev5;
  __m128  *dpc, *dpp1, *dpp2, *dpp3, *dpp4, *dpp5;
  __m128  *tp7, *tp_dd, *tp_mm;
  __m128  *ivxf;
  int      Q = p7O_NQF(om_fs->M);
  int      i, q, j, r;
  int      c1, c2, c3, c4, c5;
  int      t, u, v, w, x;
  float    scale;

  ivxf = ov->ivx[0];

  if (om_fs->codon_lengths != 5)
    ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");
#if eslDEBUGLEVEL > 0
  if (om_fs->M > bck->allocQ4*4)
    ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few columns)");
  if (bck->validR < PARSER_ROWS_BWD)
    ESL_EXCEPTION(eslEINVAL, "DP matrix needs at least PARSER_ROWS_BWD MDI rows");
  if (L >= bck->allocXR)
    ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few X rows)");
  if (! p7_fs_oprofile_IsLocal(om_fs))
    ESL_EXCEPTION(eslEINVAL, "Backward implementation assumes local alignment mode");
#endif

  bck->M              = om_fs->M;
  bck->L              = L;
  bck->has_own_scales = FALSE;
  bck->totscale       = 0.0;
  zerov = _mm_setzero_ps();

  for (r = 0; r < PARSER_ROWS_BWD; r++)
    for (q = 0; q < Q; q++)
      MMO(bck->dpf[r],q) = IMO(bck->dpf[r],q) = DMO(bck->dpf[r],q) = zerov;

  for (r = 0; r < PARSER_ROWS_BWD; r++)
    xN_buf[r] = xB_buf[r] = xJ_buf[r] = xC_buf[r] = 0.0f;

  xC_buf[(L+1) % PARSER_ROWS_BWD] = om_fs->xf[p7O_C][p7O_MOVE];
  xC_buf[(L+2) % PARSER_ROWS_BWD] = om_fs->xf[p7O_C][p7O_MOVE];

  /* Initialization: row L */
  i = L;  b = i % PARSER_ROWS_BWD;  curr = b;

  xC = om_fs->xf[p7O_C][p7O_MOVE];
  xN = xB = xJ = 0.0f;
  xE = xC * om_fs->xf[p7O_E][p7O_MOVE];
  xEv = _mm_set1_ps(xE);

  dpc = bck->dpf[curr];
  for (q = 0; q < Q; q++)
    { MMO(dpc,q) = xEv; DMO(dpc,q) = xEv; IMO(dpc,q) = zerov; }

  dcv   = esl_sse_leftshiftz_float(DMO(dpc,0));
  tp_dd = om_fs->tfv + 8*Q - 1;
  for (q = Q-1; q >= 0; q--)
    { sv=_mm_mul_ps(dcv,*tp_dd); tp_dd--; DMO(dpc,q)=_mm_add_ps(DMO(dpc,q),sv); dcv=DMO(dpc,q); }
  dcv = sv;
  for (j = 1; j < 4; j++)
    { dcv=esl_sse_leftshiftz_float(dcv); tp_dd=om_fs->tfv+8*Q-1;
      for (q=Q-1;q>=0;q--) { sv=_mm_mul_ps(dcv,*tp_dd); tp_dd--; DMO(dpc,q)=_mm_add_ps(DMO(dpc,q),sv); dcv=_mm_mul_ps(dcv,*(tp_dd+1)); } }

  dcv=esl_sse_leftshiftz_float(DMO(dpc,0)); tp7=om_fs->tfv+7*Q-3;
  for (q=Q-1;q>=0;q--)
    { MMO(dpc,q)=_mm_add_ps(MMO(dpc,q),_mm_mul_ps(dcv,*tp7)); tp7-=7; dcv=DMO(dpc,q); }

  scale = fwd->xmx[L*p7X_NXCELLS+p7X_SCALE];
  bck->xmx[L*p7X_NXCELLS+p7X_SCALE] = scale;
  if (scale > 1.0f)
    { float sf = 1.0f/scale;
      xN*=sf; xJ*=sf; xC*=sf; xB*=sf; xE*=sf;
      xEv=_mm_set1_ps(sf);
      for (r=0;r<PARSER_ROWS_BWD;r++) for (q=0;q<Q;q++)
        { MMO(bck->dpf[r],q)=_mm_mul_ps(MMO(bck->dpf[r],q),xEv);
          DMO(bck->dpf[r],q)=_mm_mul_ps(DMO(bck->dpf[r],q),xEv);
          IMO(bck->dpf[r],q)=_mm_mul_ps(IMO(bck->dpf[r],q),xEv); }
      bck->totscale += log(scale); }

  xN_buf[b]=xN; xB_buf[b]=xB; xJ_buf[b]=xJ; xC_buf[b]=xC;
  bck->xmx[L*p7X_NXCELLS+p7X_E]=xE; bck->xmx[L*p7X_NXCELLS+p7X_N]=xN;
  bck->xmx[L*p7X_NXCELLS+p7X_J]=xJ; bck->xmx[L*p7X_NXCELLS+p7X_B]=xB;
  bck->xmx[L*p7X_NXCELLS+p7X_C]=xC;

  /* Main recurrence: i = L-1 down to 1 */
  t = u = v = w = x = p7P_MAXCODONS5;

  for (i = L-1; i >= 1; i--)
    {
      t=u; u=v; v=w; w=x;
      if (dsq[i+1] < p7P_MAXNUC) x = dsq[i+1]; else x = p7P_MAXCODONS5;

      c1 = p7P_CODON1_FS5(x);           c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
      c2 = p7P_CODON2_FS5(x, w);        c2 = p7P_MINIDX(c2, p7P_DEGEN5_QC1);
      c3 = p7P_CODON3_FS5(x, w, v);     c3 = p7P_MINIDX(c3, p7P_DEGEN5_C);
      c4 = p7P_CODON4_FS5(x, w, v, u);  c4 = p7P_MINIDX(c4, p7P_DEGEN5_QC1);
      c5 = p7P_CODON5_FS5(x, w, v, u, t); c5 = p7P_MINIDX(c5, p7P_DEGEN5_QC2);

      curr  =  i       % PARSER_ROWS_BWD;
      prev1 = (i+1)    % PARSER_ROWS_BWD;
      prev2 = (i+2)    % PARSER_ROWS_BWD;
      prev3 = (i+3)    % PARSER_ROWS_BWD;
      prev4 = (i+4)    % PARSER_ROWS_BWD;
      prev5 = (i+5)    % PARSER_ROWS_BWD;
      b     =  i       % PARSER_ROWS_BWD;
      b3    = (i+3)    % PARSER_ROWS_BWD;

      dpc  = bck->dpf[curr];
      dpp1 = bck->dpf[prev1];
      dpp2 = bck->dpf[prev2];
      dpp3 = bck->dpf[prev3];
      dpp4 = bck->dpf[prev4];
      dpp5 = bck->dpf[prev5];

      for (q = 0; q < Q; q++)
        ivxf[q] = _mm_add_ps(
                    _mm_add_ps(
                      _mm_add_ps(_mm_mul_ps(MMO(dpp1,q), om_fs->rfv[c1][q]),
                                 _mm_mul_ps(MMO(dpp2,q), om_fs->rfv[c2][q])),
                      _mm_add_ps(_mm_mul_ps(MMO(dpp3,q), om_fs->rfv[c3][q]),
                                 _mm_mul_ps(MMO(dpp4,q), om_fs->rfv[c4][q]))),
                    _mm_mul_ps(MMO(dpp5,q), om_fs->rfv[c5][q]));

      xBv = zerov; tp7 = om_fs->tfv;
      for (q=0;q<Q;q++) { xBv=_mm_add_ps(xBv,_mm_mul_ps(ivxf[q],*tp7)); tp7+=7; }
      xBv=_mm_add_ps(xBv,_mm_shuffle_ps(xBv,xBv,_MM_SHUFFLE(0,3,2,1)));
      xBv=_mm_add_ps(xBv,_mm_shuffle_ps(xBv,xBv,_MM_SHUFFLE(1,0,3,2)));
      _mm_store_ss(&xB, xBv);

      xC = xC_buf[b3] * om_fs->xf[p7O_C][p7O_LOOP];
      xJ = xJ_buf[b3] * om_fs->xf[p7O_J][p7O_LOOP] + xB * om_fs->xf[p7O_J][p7O_MOVE];
      xN = xN_buf[b3] * om_fs->xf[p7O_N][p7O_LOOP] + xB * om_fs->xf[p7O_N][p7O_MOVE];
      xE = xJ * om_fs->xf[p7O_E][p7O_LOOP] + xC * om_fs->xf[p7O_E][p7O_MOVE];
      xEv = _mm_set1_ps(xE);

      for (q=0;q<Q;q++) { MMO(dpc,q)=xEv; DMO(dpc,q)=xEv; IMO(dpc,q)=zerov; }

      tp7 = om_fs->tfv;
      for (q=0;q<Q;q++)
        { MMO(dpc,q)=_mm_add_ps(MMO(dpc,q),_mm_mul_ps(IMO(dpp3,q),*(tp7+5)));
          IMO(dpc,q)=_mm_add_ps(IMO(dpc,q),_mm_mul_ps(IMO(dpp3,q),*(tp7+6))); tp7+=7; }

      ivx_carry=esl_sse_leftshiftz_float(ivxf[0]);
      tmmv=esl_sse_leftshiftz_float(*(om_fs->tfv+1));
      timv=esl_sse_leftshiftz_float(*(om_fs->tfv+2));
      tdmv=esl_sse_leftshiftz_float(*(om_fs->tfv+3));
      tp_mm=om_fs->tfv;
      for (q=Q-1;q>=0;q--)
        { MMO(dpc,q)=_mm_add_ps(MMO(dpc,q),_mm_mul_ps(ivx_carry,tmmv));
          IMO(dpc,q)=_mm_add_ps(IMO(dpc,q),_mm_mul_ps(ivx_carry,timv));
          DMO(dpc,q)=_mm_add_ps(DMO(dpc,q),_mm_mul_ps(ivx_carry,tdmv));
          ivx_carry=ivxf[q];
          tp_mm=(tp_mm==om_fs->tfv) ? om_fs->tfv+7*(Q-1) : tp_mm-7;
          tmmv=*(tp_mm+1); timv=*(tp_mm+2); tdmv=*(tp_mm+3); }

      dcv=esl_sse_leftshiftz_float(DMO(dpc,0)); tp_dd=om_fs->tfv+8*Q-1;
      for (q=Q-1;q>=0;q--)
        { sv=_mm_mul_ps(dcv,*tp_dd); tp_dd--; DMO(dpc,q)=_mm_add_ps(DMO(dpc,q),sv); dcv=DMO(dpc,q); }
      dcv=sv;
      for (j=1;j<4;j++)
        { dcv=esl_sse_leftshiftz_float(dcv); tp_dd=om_fs->tfv+8*Q-1;
          for (q=Q-1;q>=0;q--) { sv=_mm_mul_ps(dcv,*tp_dd); tp_dd--; DMO(dpc,q)=_mm_add_ps(DMO(dpc,q),sv); dcv=_mm_mul_ps(dcv,*(tp_dd+1)); } }

      dcv=esl_sse_leftshiftz_float(DMO(dpc,0)); tp7=om_fs->tfv+7*Q-3;
      for (q=Q-1;q>=0;q--)
        { MMO(dpc,q)=_mm_add_ps(MMO(dpc,q),_mm_mul_ps(dcv,*tp7)); tp7-=7; dcv=DMO(dpc,q); }

      if (xB > 1.0e16f) bck->has_own_scales = TRUE;
      if (bck->has_own_scales) scale = (xB > 1.0e4f) ? xB : 1.0f;
      else                     scale = fwd->xmx[i*p7X_NXCELLS+p7X_SCALE];

      bck->xmx[i*p7X_NXCELLS+p7X_SCALE] = scale;
      if (scale > 1.0f)
        { float sf = 1.0f/scale;
          xN*=sf; xJ*=sf; xC*=sf; xB*=sf; xE*=sf;
          xEv=_mm_set1_ps(sf);
          for (r=0;r<PARSER_ROWS_BWD;r++) for (q=0;q<Q;q++)
            { MMO(bck->dpf[r],q)=_mm_mul_ps(MMO(bck->dpf[r],q),xEv);
              DMO(bck->dpf[r],q)=_mm_mul_ps(DMO(bck->dpf[r],q),xEv);
              IMO(bck->dpf[r],q)=_mm_mul_ps(IMO(bck->dpf[r],q),xEv); }
          for (r=0;r<PARSER_ROWS_BWD;r++)
            { xN_buf[r]*=sf; xB_buf[r]*=sf; xJ_buf[r]*=sf; xC_buf[r]*=sf; }
          bck->totscale += log(scale); }

      xN_buf[b]=xN; xB_buf[b]=xB; xJ_buf[b]=xJ; xC_buf[b]=xC;
      bck->xmx[i*p7X_NXCELLS+p7X_E]=xE; bck->xmx[i*p7X_NXCELLS+p7X_N]=xN;
      bck->xmx[i*p7X_NXCELLS+p7X_J]=xJ; bck->xmx[i*p7X_NXCELLS+p7X_B]=xB;
      bck->xmx[i*p7X_NXCELLS+p7X_C]=xC;
    }

  /* Termination: i=0 */
  if (dsq[1] < p7P_MAXNUC) x=dsq[1]; else x=p7P_MAXCODONS5;
  if (L>=2 && dsq[2]<p7P_MAXNUC) w=dsq[2]; else w=p7P_MAXCODONS5;
  if (L>=3 && dsq[3]<p7P_MAXNUC) v=dsq[3]; else v=p7P_MAXCODONS5;
  if (L>=4 && dsq[4]<p7P_MAXNUC) u=dsq[4]; else u=p7P_MAXCODONS5;
  if (L>=5 && dsq[5]<p7P_MAXNUC) t=dsq[5]; else t=p7P_MAXCODONS5;

  c1=p7P_CODON1_FS5(x);            c1=p7P_MINIDX(c1,p7P_DEGEN5_QC2);
  c2=p7P_CODON2_FS5(x,w);          c2=p7P_MINIDX(c2,p7P_DEGEN5_QC1);
  c3=p7P_CODON3_FS5(x,w,v);        c3=p7P_MINIDX(c3,p7P_DEGEN5_C);
  c4=p7P_CODON4_FS5(x,w,v,u);      c4=p7P_MINIDX(c4,p7P_DEGEN5_QC1);
  c5=p7P_CODON5_FS5(x,w,v,u,t);    c5=p7P_MINIDX(c5,p7P_DEGEN5_QC2);

  prev1=1%PARSER_ROWS_BWD; prev2=2%PARSER_ROWS_BWD; prev3=3%PARSER_ROWS_BWD;
  prev4=4%PARSER_ROWS_BWD; prev5=5%PARSER_ROWS_BWD;
  dpp1=bck->dpf[prev1]; dpp2=bck->dpf[prev2]; dpp3=bck->dpf[prev3];
  dpp4=bck->dpf[prev4]; dpp5=bck->dpf[prev5];

  for (q=0;q<Q;q++)
    ivxf[q]=_mm_add_ps(
              _mm_add_ps(
                _mm_add_ps(_mm_mul_ps(MMO(dpp1,q),om_fs->rfv[c1][q]),
                           _mm_mul_ps(MMO(dpp2,q),om_fs->rfv[c2][q])),
                _mm_add_ps(_mm_mul_ps(MMO(dpp3,q),om_fs->rfv[c3][q]),
                           _mm_mul_ps(MMO(dpp4,q),om_fs->rfv[c4][q]))),
              _mm_mul_ps(MMO(dpp5,q),om_fs->rfv[c5][q]));

  xBv=zerov; tp7=om_fs->tfv;
  for (q=0;q<Q;q++) { xBv=_mm_add_ps(xBv,_mm_mul_ps(ivxf[q],*tp7)); tp7+=7; }
  xBv=_mm_add_ps(xBv,_mm_shuffle_ps(xBv,xBv,_MM_SHUFFLE(0,3,2,1)));
  xBv=_mm_add_ps(xBv,_mm_shuffle_ps(xBv,xBv,_MM_SHUFFLE(1,0,3,2)));
  _mm_store_ss(&xB, xBv);

  xN = xN_buf[3%PARSER_ROWS_BWD] * om_fs->xf[p7O_N][p7O_LOOP]
     + xB                         * om_fs->xf[p7O_N][p7O_MOVE];

  bck->xmx[p7X_B]=xB; bck->xmx[p7X_N]=xN; bck->xmx[p7X_J]=0.0f;
  bck->xmx[p7X_C]=0.0f; bck->xmx[p7X_E]=0.0f; bck->xmx[p7X_SCALE]=1.0f;

  {
    float xNtot = xN + xN_buf[1%PARSER_ROWS_BWD] + xN_buf[2%PARSER_ROWS_BWD];

    if      (isnan(xNtot))           ESL_EXCEPTION(eslERANGE, "backward score is NaN");
    else if (isinf(xNtot) == 1)      ESL_EXCEPTION(eslERANGE, "backward score overflow (is infinity)");
    else if (L > 0 && xNtot == 0.0f) { if (opt_sc) *opt_sc = -eslINFINITY; return eslERANGE; }
    if (opt_sc) *opt_sc = bck->totscale + logf(xNtot);
  }

  return eslOK;
}

/* Function:  p7_Forward_Frameshift_sse()
 * Synopsis:  SSE frameshift-aware Forward, 5 codon lengths, full O(ML) matrix.
 */
int
p7_Forward_Frameshift_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_OMX *fwd, P7_OIVX *ov, float *opt_sc)
{
  register __m128 mpv1, dpv1, ipv1;
  register __m128 sv;
  register __m128 msv;
  register __m128 dcv;
  register __m128 xEv;
  register __m128 xBv1;
  __m128   zerov;
  float    xN, xE, xB, xC, xJ;
  float    xN_buf[PARSER_ROWS_FWD];
  float    xB_buf[PARSER_ROWS_FWD];
  float    xJ_buf[PARSER_ROWS_FWD];
  float    xC_buf[PARSER_ROWS_FWD];
  int      b, b1, b3;
  int      ivx_1, ivx_2, ivx_3, ivx_4, ivx_5;
  __m128  *dpc, *dpp1, *dpp3;
  __m128  *tp;
  __m128 **ivx = ov->ivx;
  int      Q = p7O_NQF(om_fs->M);
  int      i, q, j, r;
  int      c1, c2, c3, c4, c5;
  int      t, u, v, w, x;
  float    insert_adj;

  if (om_fs->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");

  fwd->M              = om_fs->M;
  fwd->L              = L;
  fwd->has_own_scales = TRUE;
  fwd->totscale       = 0.0;
  zerov = _mm_setzero_ps();

  { __m128 *drow = fwd->dpf[0];
    for (q = 0; q < Q; q++)
      MMO_FS(drow,q,p7X_FS_C0) = MMO_FS(drow,q,p7X_FS_C1) =
      MMO_FS(drow,q,p7X_FS_C2) = MMO_FS(drow,q,p7X_FS_C3) =
      MMO_FS(drow,q,p7X_FS_C4) = MMO_FS(drow,q,p7X_FS_C5) =
      DMO_FS(drow,q)            = IMO_FS(drow,q)            = zerov; }

  for (r = 0; r < p7P_5CODONS; r++)
    for (q = 0; q < Q; q++)
      IVX(r, q) = zerov;

  for (r = 0; r < PARSER_ROWS_FWD; r++)
    xN_buf[r] = xB_buf[r] = xJ_buf[r] = xC_buf[r] = 0.0f;
  xN_buf[0] = xN_buf[1] = xN_buf[2] = 1.0f;
  xB_buf[0] = xB_buf[1] = xB_buf[2] = om_fs->xf[p7O_N][p7O_MOVE];

  for (r = 0; r < 3; r++)
    {
      fwd->xmx[r*p7X_NXCELLS+p7X_SCALE] = 1.0f;
      fwd->xmx[r*p7X_NXCELLS+p7X_E]     = 0.0f;
      fwd->xmx[r*p7X_NXCELLS+p7X_N]     = 1.0f;
      fwd->xmx[r*p7X_NXCELLS+p7X_J]     = 0.0f;
      fwd->xmx[r*p7X_NXCELLS+p7X_B]     = om_fs->xf[p7O_N][p7O_MOVE];
      fwd->xmx[r*p7X_NXCELLS+p7X_C]     = 0.0f;
    }

  t = u = v = w = p7P_MAXCODONS5;
  if (dsq[1] < p7P_MAXNUC) x = dsq[1]; else x = p7P_MAXCODONS5;

  /* Initialization: i=1 */
  i = 1;
  c1 = p7P_CODON1_FS5(x); c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
  ivx_1 = 1;
  dpc  = fwd->dpf[1];
  dpp1 = fwd->dpf[0];

  xBv1 = _mm_set1_ps(xB_buf[0]);
  tp = om_fs->tfv; dcv = zerov; xEv = zerov;
  mpv1 = dpv1 = ipv1 = zerov;

  for (q = 0; q < Q; q++)
    {
      sv  =                _mm_mul_ps(xBv1, *tp); tp++;
      sv  = _mm_add_ps(sv, _mm_mul_ps(mpv1, *tp)); tp++;
      sv  = _mm_add_ps(sv, _mm_mul_ps(ipv1, *tp)); tp++;
      sv  = _mm_add_ps(sv, _mm_mul_ps(dpv1, *tp)); tp++;
      IVX(ivx_1, q) = sv;
      msv = _mm_mul_ps(sv, om_fs->rfv[c1][q]);
      MMO_FS(dpc,q,p7X_FS_C0) = msv; MMO_FS(dpc,q,p7X_FS_C1) = msv;
      MMO_FS(dpc,q,p7X_FS_C2) = zerov; MMO_FS(dpc,q,p7X_FS_C3) = zerov;
      MMO_FS(dpc,q,p7X_FS_C4) = zerov; MMO_FS(dpc,q,p7X_FS_C5) = zerov;
      xEv = _mm_add_ps(xEv, msv);
      DMO_FS(dpc,q) = dcv;
      dcv = _mm_mul_ps(msv, *tp); tp++;
      IMO_FS(dpc,q) = zerov;
      tp += 2;
    }

  dcv=esl_sse_rightshiftz_float(dcv); DMO_FS(dpc,0)=zerov; tp=om_fs->tfv+7*Q;
  for (q=0;q<Q;q++) { DMO_FS(dpc,q)=_mm_add_ps(dcv,DMO_FS(dpc,q)); dcv=_mm_mul_ps(DMO_FS(dpc,q),*tp); tp++; }
  if (om_fs->M<100)
    { for (j=1;j<4;j++) { dcv=esl_sse_rightshiftz_float(dcv); tp=om_fs->tfv+7*Q;
        for (q=0;q<Q;q++) { DMO_FS(dpc,q)=_mm_add_ps(dcv,DMO_FS(dpc,q)); dcv=_mm_mul_ps(dcv,*tp); tp++; } } }
  else
    { for (j=1;j<4;j++) { register __m128 cv; dcv=esl_sse_rightshiftz_float(dcv); tp=om_fs->tfv+7*Q; cv=zerov;
        for (q=0;q<Q;q++) { sv=_mm_add_ps(dcv,DMO_FS(dpc,q)); cv=_mm_or_ps(cv,_mm_cmpgt_ps(sv,DMO_FS(dpc,q))); DMO_FS(dpc,q)=sv; dcv=_mm_mul_ps(dcv,*tp); tp++; }
        if (!_mm_movemask_ps(cv)) break; } }

  for (q=0;q<Q;q++) xEv=_mm_add_ps(DMO_FS(dpc,q),xEv);
  xEv=_mm_add_ps(xEv,_mm_shuffle_ps(xEv,xEv,_MM_SHUFFLE(0,3,2,1)));
  xEv=_mm_add_ps(xEv,_mm_shuffle_ps(xEv,xEv,_MM_SHUFFLE(1,0,3,2)));
  _mm_store_ss(&xE, xEv);

  xN=1.0f; xJ=xE*om_fs->xf[p7O_E][p7O_LOOP]; xC=xE*om_fs->xf[p7O_E][p7O_MOVE];
  xB=xN*om_fs->xf[p7O_N][p7O_MOVE]+xJ*om_fs->xf[p7O_J][p7O_MOVE];

  if (xE > 1.0e4f)
    {
      float scale_factor = 1.0f / xE;
      xN*=scale_factor; xJ*=scale_factor; xC*=scale_factor; xB*=scale_factor;
      xEv=_mm_set1_ps(scale_factor);
      for (q=0;q<Q;q++)
        { MMO_FS(dpc,q,p7X_FS_C0)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C0),xEv);
          MMO_FS(dpc,q,p7X_FS_C1)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C1),xEv);
          MMO_FS(dpc,q,p7X_FS_C2)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C2),xEv);
          MMO_FS(dpc,q,p7X_FS_C3)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C3),xEv);
          MMO_FS(dpc,q,p7X_FS_C4)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C4),xEv);
          MMO_FS(dpc,q,p7X_FS_C5)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C5),xEv);
          DMO_FS(dpc,q)=_mm_mul_ps(DMO_FS(dpc,q),xEv);
          IMO_FS(dpc,q)=_mm_mul_ps(IMO_FS(dpc,q),xEv); }
      for (r=0;r<p7P_5CODONS;r++) for (q=0;q<Q;q++) IVX(r,q)=_mm_mul_ps(IVX(r,q),xEv);
      for (r=0;r<PARSER_ROWS_FWD;r++)
        { xN_buf[r]*=scale_factor; xB_buf[r]*=scale_factor; xJ_buf[r]*=scale_factor; xC_buf[r]*=scale_factor; }
      fwd->xmx[1*p7X_NXCELLS+p7X_SCALE]=xE; fwd->totscale+=log(xE); xE=1.0f;
    }
  else fwd->xmx[1*p7X_NXCELLS+p7X_SCALE]=1.0f;

  xN_buf[1]=xN; xB_buf[1]=xB; xJ_buf[1]=xJ; xC_buf[1]=xC;
  fwd->xmx[1*p7X_NXCELLS+p7X_E]=xE; fwd->xmx[1*p7X_NXCELLS+p7X_N]=xN;
  fwd->xmx[1*p7X_NXCELLS+p7X_J]=xJ; fwd->xmx[1*p7X_NXCELLS+p7X_B]=xB;
  fwd->xmx[1*p7X_NXCELLS+p7X_C]=xC;

  /* Initialization: i=2 */
  i=2; t=u=v=p7P_MAXCODONS5; w=x;
  if (dsq[2]<p7P_MAXNUC) x=dsq[2]; else x=p7P_MAXCODONS5;

  c1=p7P_CODON1_FS5(x);    c1=p7P_MINIDX(c1,p7P_DEGEN5_QC2);
  c2=p7P_CODON2_FS5(w,x);  c2=p7P_MINIDX(c2,p7P_DEGEN5_QC1);
  ivx_1=2; ivx_2=1;
  dpc=fwd->dpf[2]; dpp1=fwd->dpf[1];

  xBv1=_mm_set1_ps(xB_buf[1]); tp=om_fs->tfv; dcv=zerov; xEv=zerov;
  mpv1=esl_sse_rightshiftz_float(MMO_FS(dpp1,Q-1,p7X_FS_C0));
  dpv1=esl_sse_rightshiftz_float(DMO_FS(dpp1,Q-1));
  ipv1=esl_sse_rightshiftz_float(IMO_FS(dpp1,Q-1));

  for (q=0;q<Q;q++)
    {
      sv =               _mm_mul_ps(xBv1,*tp); tp++;
      sv=_mm_add_ps(sv,  _mm_mul_ps(mpv1,*tp)); tp++;
      sv=_mm_add_ps(sv,  _mm_mul_ps(ipv1,*tp)); tp++;
      sv=_mm_add_ps(sv,  _mm_mul_ps(dpv1,*tp)); tp++;
      IVX(ivx_1,q)=sv;
      __m128 mc1=_mm_mul_ps(sv,             om_fs->rfv[c1][q]);
      __m128 mc2=_mm_mul_ps(IVX(ivx_2,q),  om_fs->rfv[c2][q]);
      msv=_mm_add_ps(mc1,mc2); xEv=_mm_add_ps(xEv,msv);
      mpv1=MMO_FS(dpp1,q,p7X_FS_C0); dpv1=DMO_FS(dpp1,q); ipv1=IMO_FS(dpp1,q);
      MMO_FS(dpc,q,p7X_FS_C0)=msv; MMO_FS(dpc,q,p7X_FS_C1)=mc1; MMO_FS(dpc,q,p7X_FS_C2)=mc2;
      MMO_FS(dpc,q,p7X_FS_C3)=zerov; MMO_FS(dpc,q,p7X_FS_C4)=zerov; MMO_FS(dpc,q,p7X_FS_C5)=zerov;
      DMO_FS(dpc,q)=dcv; dcv=_mm_mul_ps(msv,*tp); tp++;
      IMO_FS(dpc,q)=zerov; tp+=2;
    }

  dcv=esl_sse_rightshiftz_float(dcv); DMO_FS(dpc,0)=zerov; tp=om_fs->tfv+7*Q;
  for (q=0;q<Q;q++) { DMO_FS(dpc,q)=_mm_add_ps(dcv,DMO_FS(dpc,q)); dcv=_mm_mul_ps(DMO_FS(dpc,q),*tp); tp++; }
  if (om_fs->M<100)
    { for (j=1;j<4;j++) { dcv=esl_sse_rightshiftz_float(dcv); tp=om_fs->tfv+7*Q;
        for (q=0;q<Q;q++) { DMO_FS(dpc,q)=_mm_add_ps(dcv,DMO_FS(dpc,q)); dcv=_mm_mul_ps(dcv,*tp); tp++; } } }
  else
    { for (j=1;j<4;j++) { register __m128 cv; dcv=esl_sse_rightshiftz_float(dcv); tp=om_fs->tfv+7*Q; cv=zerov;
        for (q=0;q<Q;q++) { sv=_mm_add_ps(dcv,DMO_FS(dpc,q)); cv=_mm_or_ps(cv,_mm_cmpgt_ps(sv,DMO_FS(dpc,q))); DMO_FS(dpc,q)=sv; dcv=_mm_mul_ps(dcv,*tp); tp++; }
        if (!_mm_movemask_ps(cv)) break; } }

  for (q=0;q<Q;q++) xEv=_mm_add_ps(DMO_FS(dpc,q),xEv);
  xEv=_mm_add_ps(xEv,_mm_shuffle_ps(xEv,xEv,_MM_SHUFFLE(0,3,2,1)));
  xEv=_mm_add_ps(xEv,_mm_shuffle_ps(xEv,xEv,_MM_SHUFFLE(1,0,3,2)));
  _mm_store_ss(&xE,xEv);

  xN=1.0f; xJ=xE*om_fs->xf[p7O_E][p7O_LOOP]; xC=xE*om_fs->xf[p7O_E][p7O_MOVE];
  xB=xN*om_fs->xf[p7O_N][p7O_MOVE]+xJ*om_fs->xf[p7O_J][p7O_MOVE];

  if (xE>1.0e4f)
    {
      float scale_factor=1.0f/xE;
      xN*=scale_factor; xJ*=scale_factor; xC*=scale_factor; xB*=scale_factor;
      xEv=_mm_set1_ps(scale_factor);
      for (q=0;q<Q;q++)
        { MMO_FS(dpc,q,p7X_FS_C0)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C0),xEv);
          MMO_FS(dpc,q,p7X_FS_C1)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C1),xEv);
          MMO_FS(dpc,q,p7X_FS_C2)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C2),xEv);
          MMO_FS(dpc,q,p7X_FS_C3)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C3),xEv);
          MMO_FS(dpc,q,p7X_FS_C4)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C4),xEv);
          MMO_FS(dpc,q,p7X_FS_C5)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C5),xEv);
          DMO_FS(dpc,q)=_mm_mul_ps(DMO_FS(dpc,q),xEv); IMO_FS(dpc,q)=_mm_mul_ps(IMO_FS(dpc,q),xEv); }
      for (r=0;r<p7P_5CODONS;r++) for (q=0;q<Q;q++) IVX(r,q)=_mm_mul_ps(IVX(r,q),xEv);
      for (r=0;r<PARSER_ROWS_FWD;r++)
        { xN_buf[r]*=scale_factor; xB_buf[r]*=scale_factor; xJ_buf[r]*=scale_factor; xC_buf[r]*=scale_factor; }
      fwd->xmx[2*p7X_NXCELLS+p7X_SCALE]=xE; fwd->totscale+=log(xE); xE=1.0f;
    }
  else fwd->xmx[2*p7X_NXCELLS+p7X_SCALE]=1.0f;

  xN_buf[2]=xN; xB_buf[2]=xB; xJ_buf[2]=xJ; xC_buf[2]=xC;
  fwd->xmx[2*p7X_NXCELLS+p7X_E]=xE; fwd->xmx[2*p7X_NXCELLS+p7X_N]=xN;
  fwd->xmx[2*p7X_NXCELLS+p7X_J]=xJ; fwd->xmx[2*p7X_NXCELLS+p7X_B]=xB;
  fwd->xmx[2*p7X_NXCELLS+p7X_C]=xC;

  /* Main recurrence: i = 3..L */
  for (i = 3; i <= L; i++)
    {
      t=u; u=v; v=w; w=x;
      if (dsq[i]<p7P_MAXNUC) x=dsq[i]; else x=p7P_MAXCODONS5;

      c1=p7P_CODON1_FS5(x);           c1=p7P_MINIDX(c1,p7P_DEGEN5_QC2);
      c2=p7P_CODON2_FS5(w,x);         c2=p7P_MINIDX(c2,p7P_DEGEN5_QC1);
      c3=p7P_CODON3_FS5(v,w,x);       c3=p7P_MINIDX(c3,p7P_DEGEN5_C);
      c4=p7P_CODON4_FS5(u,v,w,x);     c4=p7P_MINIDX(c4,p7P_DEGEN5_QC1);
      c5=p7P_CODON5_FS5(t,u,v,w,x);   c5=p7P_MINIDX(c5,p7P_DEGEN5_QC2);

      ivx_1=    i              %p7P_5CODONS;
      ivx_2=((i-1)%p7P_5CODONS+p7P_5CODONS)%p7P_5CODONS;
      ivx_3=((i-2)%p7P_5CODONS+p7P_5CODONS)%p7P_5CODONS;
      ivx_4=((i-3)%p7P_5CODONS+p7P_5CODONS)%p7P_5CODONS;
      ivx_5=((i-4)%p7P_5CODONS+p7P_5CODONS)%p7P_5CODONS;

      b =    i              %PARSER_ROWS_FWD;
      b1=((i-1)%PARSER_ROWS_FWD+PARSER_ROWS_FWD)%PARSER_ROWS_FWD;
      b3=((i-3)%PARSER_ROWS_FWD+PARSER_ROWS_FWD)%PARSER_ROWS_FWD;

      dpc=fwd->dpf[i]; dpp1=fwd->dpf[i-1]; dpp3=fwd->dpf[i-3];

      insert_adj = 1.0f
                 / (fwd->xmx[(i-2)*p7X_NXCELLS+p7X_SCALE]
                 *  fwd->xmx[(i-1)*p7X_NXCELLS+p7X_SCALE]);

      mpv1=esl_sse_rightshiftz_float(MMO_FS(dpp1,Q-1,p7X_FS_C0));
      dpv1=esl_sse_rightshiftz_float(DMO_FS(dpp1,Q-1));
      ipv1=esl_sse_rightshiftz_float(IMO_FS(dpp1,Q-1));

      xBv1=_mm_set1_ps(xB_buf[b1]);
      tp=om_fs->tfv; dcv=zerov; xEv=zerov;

      register __m128 adj_v=_mm_set1_ps(insert_adj);

      for (q=0;q<Q;q++)
        {
          sv =               _mm_mul_ps(xBv1,*tp); tp++;
          sv=_mm_add_ps(sv,  _mm_mul_ps(mpv1,*tp)); tp++;
          sv=_mm_add_ps(sv,  _mm_mul_ps(ipv1,*tp)); tp++;
          sv=_mm_add_ps(sv,  _mm_mul_ps(dpv1,*tp)); tp++;
          IVX(ivx_1,q)=sv;
          __m128 mc1=_mm_mul_ps(sv,              om_fs->rfv[c1][q]);
          __m128 mc2=_mm_mul_ps(IVX(ivx_2,q),   om_fs->rfv[c2][q]);
          __m128 mc3=_mm_mul_ps(IVX(ivx_3,q),   om_fs->rfv[c3][q]);
          __m128 mc4=_mm_mul_ps(IVX(ivx_4,q),   om_fs->rfv[c4][q]);
          __m128 mc5=_mm_mul_ps(IVX(ivx_5,q),   om_fs->rfv[c5][q]);
          msv=_mm_add_ps(_mm_add_ps(_mm_add_ps(mc1,mc2),_mm_add_ps(mc3,mc4)),mc5);
          xEv=_mm_add_ps(xEv,msv);
          mpv1=MMO_FS(dpp1,q,p7X_FS_C0); dpv1=DMO_FS(dpp1,q); ipv1=IMO_FS(dpp1,q);
          MMO_FS(dpc,q,p7X_FS_C0)=msv; MMO_FS(dpc,q,p7X_FS_C1)=mc1; MMO_FS(dpc,q,p7X_FS_C2)=mc2;
          MMO_FS(dpc,q,p7X_FS_C3)=mc3; MMO_FS(dpc,q,p7X_FS_C4)=mc4; MMO_FS(dpc,q,p7X_FS_C5)=mc5;
          DMO_FS(dpc,q)=dcv; dcv=_mm_mul_ps(msv,*tp); tp++;
          sv=_mm_mul_ps(_mm_mul_ps(MMO_FS(dpp3,q,p7X_FS_C0),adj_v),*tp); tp++;
          IMO_FS(dpc,q)=_mm_add_ps(sv,_mm_mul_ps(_mm_mul_ps(IMO_FS(dpp3,q),adj_v),*tp)); tp++;
        }

      dcv=esl_sse_rightshiftz_float(dcv); DMO_FS(dpc,0)=zerov; tp=om_fs->tfv+7*Q;
      for (q=0;q<Q;q++) { DMO_FS(dpc,q)=_mm_add_ps(dcv,DMO_FS(dpc,q)); dcv=_mm_mul_ps(DMO_FS(dpc,q),*tp); tp++; }
      if (om_fs->M<100)
        { for (j=1;j<4;j++) { dcv=esl_sse_rightshiftz_float(dcv); tp=om_fs->tfv+7*Q;
            for (q=0;q<Q;q++) { DMO_FS(dpc,q)=_mm_add_ps(dcv,DMO_FS(dpc,q)); dcv=_mm_mul_ps(dcv,*tp); tp++; } } }
      else
        { for (j=1;j<4;j++) { register __m128 cv; dcv=esl_sse_rightshiftz_float(dcv); tp=om_fs->tfv+7*Q; cv=zerov;
            for (q=0;q<Q;q++) { sv=_mm_add_ps(dcv,DMO_FS(dpc,q)); cv=_mm_or_ps(cv,_mm_cmpgt_ps(sv,DMO_FS(dpc,q))); DMO_FS(dpc,q)=sv; dcv=_mm_mul_ps(dcv,*tp); tp++; }
            if (!_mm_movemask_ps(cv)) break; } }

      for (q=0;q<Q;q++) xEv=_mm_add_ps(DMO_FS(dpc,q),xEv);
      xEv=_mm_add_ps(xEv,_mm_shuffle_ps(xEv,xEv,_MM_SHUFFLE(0,3,2,1)));
      xEv=_mm_add_ps(xEv,_mm_shuffle_ps(xEv,xEv,_MM_SHUFFLE(1,0,3,2)));
      _mm_store_ss(&xE,xEv);

      xN=xN_buf[b3]*om_fs->xf[p7O_N][p7O_LOOP];
      xJ=xJ_buf[b3]*om_fs->xf[p7O_J][p7O_LOOP]+xE*om_fs->xf[p7O_E][p7O_LOOP];
      xC=xC_buf[b3]*om_fs->xf[p7O_C][p7O_LOOP]+xE*om_fs->xf[p7O_E][p7O_MOVE];
      xB=xN*om_fs->xf[p7O_N][p7O_MOVE]+xJ*om_fs->xf[p7O_J][p7O_MOVE];

      if (xE>1.0e4f)
        {
          float scale_factor=1.0f/xE;
          xN*=scale_factor; xJ*=scale_factor; xC*=scale_factor; xB*=scale_factor;
          xEv=_mm_set1_ps(scale_factor);
          for (q=0;q<Q;q++)
            { MMO_FS(dpc,q,p7X_FS_C0)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C0),xEv);
              MMO_FS(dpc,q,p7X_FS_C1)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C1),xEv);
              MMO_FS(dpc,q,p7X_FS_C2)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C2),xEv);
              MMO_FS(dpc,q,p7X_FS_C3)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C3),xEv);
              MMO_FS(dpc,q,p7X_FS_C4)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C4),xEv);
              MMO_FS(dpc,q,p7X_FS_C5)=_mm_mul_ps(MMO_FS(dpc,q,p7X_FS_C5),xEv);
              DMO_FS(dpc,q)=_mm_mul_ps(DMO_FS(dpc,q),xEv); IMO_FS(dpc,q)=_mm_mul_ps(IMO_FS(dpc,q),xEv); }
          for (r=0;r<p7P_5CODONS;r++) for (q=0;q<Q;q++) IVX(r,q)=_mm_mul_ps(IVX(r,q),xEv);
          for (r=0;r<PARSER_ROWS_FWD;r++)
            { xN_buf[r]*=scale_factor; xB_buf[r]*=scale_factor; xJ_buf[r]*=scale_factor; xC_buf[r]*=scale_factor; }
          fwd->xmx[i*p7X_NXCELLS+p7X_SCALE]=xE; fwd->totscale+=log(xE); xE=1.0f;
        }
      else fwd->xmx[i*p7X_NXCELLS+p7X_SCALE]=1.0f;

      xN_buf[b]=xN; xB_buf[b]=xB; xJ_buf[b]=xJ; xC_buf[b]=xC;
      fwd->xmx[i*p7X_NXCELLS+p7X_E]=xE; fwd->xmx[i*p7X_NXCELLS+p7X_N]=xN;
      fwd->xmx[i*p7X_NXCELLS+p7X_J]=xJ; fwd->xmx[i*p7X_NXCELLS+p7X_B]=xB;
      fwd->xmx[i*p7X_NXCELLS+p7X_C]=xC;
    }

  {
    float xCL   = xC_buf[   L    %PARSER_ROWS_FWD];
    float xCLm1 = xC_buf[((L-1)%PARSER_ROWS_FWD+PARSER_ROWS_FWD)%PARSER_ROWS_FWD];
    float xCLm2 = xC_buf[((L-2)%PARSER_ROWS_FWD+PARSER_ROWS_FWD)%PARSER_ROWS_FWD];
    float xCtot = xCL + xCLm1*om_fs->xf[p7O_C][p7O_LOOP] + xCLm2*om_fs->xf[p7O_C][p7O_LOOP];

    if      (isnan(xCtot))           ESL_EXCEPTION(eslERANGE, "forward score is NaN");
    else if (isinf(xCtot)==1)        ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");
    else if (L>1 && xCtot==0.0f)     { if (opt_sc) *opt_sc=-eslINFINITY; return eslERANGE; }
    if (opt_sc) *opt_sc=fwd->totscale+logf(xCtot*om_fs->xf[p7O_C][p7O_MOVE]);
  }

  return eslOK;
}


/* Function:  p7_Backward_Frameshift_sse()
 * Synopsis:  SSE frameshift-aware Backward, 5 codon lengths, full O(ML) matrix.
 */
int
p7_Backward_Frameshift_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
  register __m128 sv;
  register __m128 dcv;
  register __m128 ivx_carry;
  __m128   tmmv, timv, tdmv;
  register __m128 xBv;
  register __m128 xEv;
  __m128   zerov;
  float    xN, xE, xB, xC, xJ;
  float    xN_buf[PARSER_ROWS_BWD];
  float    xB_buf[PARSER_ROWS_BWD];
  float    xJ_buf[PARSER_ROWS_BWD];
  float    xC_buf[PARSER_ROWS_BWD];
  int      b, b3;
  __m128  *dpc, *dpp1, *dpp2, *dpp3, *dpp4, *dpp5;
  __m128  *zero_row = NULL;
  __m128  *zero_mem = NULL;
  __m128  *tp7, *tp_dd, *tp_mm;
  __m128  *ivxf;
  int      Q = p7O_NQF(om_fs->M);
  int      i, q, j, r;
  int      c1, c2, c3, c4, c5;
  int      t, u, v, w, x;
  float    scale;
  float    adj2, adj3, adj4, adj5;
  int      status;

  if (om_fs->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");

  ivxf = ov->ivx[0];

  bck->M=om_fs->M; bck->L=L; bck->has_own_scales=FALSE; bck->totscale=0.0;
  zerov=_mm_setzero_ps();

  ESL_ALLOC(zero_mem, sizeof(__m128)*Q*p7X_NSCELLS+15);
  zero_row=(__m128*)(((unsigned long int)zero_mem+15)&(~0xf));
  for (q=0;q<Q*p7X_NSCELLS;q++) zero_row[q]=zerov;

  for (r=0;r<PARSER_ROWS_BWD;r++)
    xN_buf[r]=xB_buf[r]=xJ_buf[r]=xC_buf[r]=0.0f;

  xC_buf[(L+1)%PARSER_ROWS_BWD]=om_fs->xf[p7O_C][p7O_MOVE];
  xC_buf[(L+2)%PARSER_ROWS_BWD]=om_fs->xf[p7O_C][p7O_MOVE];

  /* Initialization: row L */
  i=L; b=i%PARSER_ROWS_BWD;
  xC=om_fs->xf[p7O_C][p7O_MOVE]; xN=xB=xJ=0.0f;
  xE=xC*om_fs->xf[p7O_E][p7O_MOVE]; xEv=_mm_set1_ps(xE);

  dpc=bck->dpf[L];
  for (q=0;q<Q;q++) { MMO(dpc,q)=xEv; DMO(dpc,q)=xEv; IMO(dpc,q)=zerov; }

  dcv=esl_sse_leftshiftz_float(DMO(dpc,0)); tp_dd=om_fs->tfv+8*Q-1;
  for (q=Q-1;q>=0;q--)
    { sv=_mm_mul_ps(dcv,*tp_dd); tp_dd--; DMO(dpc,q)=_mm_add_ps(DMO(dpc,q),sv); dcv=DMO(dpc,q); }
  dcv=sv;
  for (j=1;j<4;j++)
    { dcv=esl_sse_leftshiftz_float(dcv); tp_dd=om_fs->tfv+8*Q-1;
      for (q=Q-1;q>=0;q--) { sv=_mm_mul_ps(dcv,*tp_dd); tp_dd--; DMO(dpc,q)=_mm_add_ps(DMO(dpc,q),sv); dcv=_mm_mul_ps(dcv,*(tp_dd+1)); } }

  dcv=esl_sse_leftshiftz_float(DMO(dpc,0)); tp7=om_fs->tfv+7*Q-3;
  for (q=Q-1;q>=0;q--) { MMO(dpc,q)=_mm_add_ps(MMO(dpc,q),_mm_mul_ps(dcv,*tp7)); tp7-=7; dcv=DMO(dpc,q); }

  scale=fwd->xmx[L*p7X_NXCELLS+p7X_SCALE]; bck->xmx[L*p7X_NXCELLS+p7X_SCALE]=scale;
  if (scale>1.0f)
    { float sf=1.0f/scale; xN*=sf; xJ*=sf; xC*=sf; xB*=sf; xE*=sf; xEv=_mm_set1_ps(sf);
      for (q=0;q<Q;q++) { MMO(dpc,q)=_mm_mul_ps(MMO(dpc,q),xEv); DMO(dpc,q)=_mm_mul_ps(DMO(dpc,q),xEv); IMO(dpc,q)=_mm_mul_ps(IMO(dpc,q),xEv); }
      bck->totscale+=log(scale); }

  xN_buf[b]=xN; xB_buf[b]=xB; xJ_buf[b]=xJ; xC_buf[b]=xC;
  bck->xmx[L*p7X_NXCELLS+p7X_E]=xE; bck->xmx[L*p7X_NXCELLS+p7X_N]=xN;
  bck->xmx[L*p7X_NXCELLS+p7X_J]=xJ; bck->xmx[L*p7X_NXCELLS+p7X_B]=xB;
  bck->xmx[L*p7X_NXCELLS+p7X_C]=xC;

  /* Main recurrence: i = L-1 down to 1 */
  t=u=v=w=x=p7P_MAXCODONS5;

  for (i=L-1;i>=1;i--)
    {
      t=u; u=v; v=w; w=x;
      if (dsq[i+1]<p7P_MAXNUC) x=dsq[i+1]; else x=p7P_MAXCODONS5;

      c1=p7P_CODON1_FS5(x);              c1=p7P_MINIDX(c1,p7P_DEGEN5_QC2);
      c2=p7P_CODON2_FS5(x,w);            c2=p7P_MINIDX(c2,p7P_DEGEN5_QC1);
      c3=p7P_CODON3_FS5(x,w,v);          c3=p7P_MINIDX(c3,p7P_DEGEN5_C);
      c4=p7P_CODON4_FS5(x,w,v,u);        c4=p7P_MINIDX(c4,p7P_DEGEN5_QC1);
      c5=p7P_CODON5_FS5(x,w,v,u,t);      c5=p7P_MINIDX(c5,p7P_DEGEN5_QC2);

      b=i%PARSER_ROWS_BWD; b3=(i+3)%PARSER_ROWS_BWD;

      dpc=bck->dpf[i];
      dpp1=bck->dpf[i+1];
      dpp2=(i+2<=L)?bck->dpf[i+2]:zero_row;
      dpp3=(i+3<=L)?bck->dpf[i+3]:zero_row;
      dpp4=(i+4<=L)?bck->dpf[i+4]:zero_row;
      dpp5=(i+5<=L)?bck->dpf[i+5]:zero_row;

      adj2=(i+2<=L)?1.0f/fwd->xmx[(i+1)*p7X_NXCELLS+p7X_SCALE]:1.0f;
      adj3=(i+3<=L)?adj2/fwd->xmx[(i+2)*p7X_NXCELLS+p7X_SCALE]:1.0f;
      adj4=(i+4<=L)?adj3/fwd->xmx[(i+3)*p7X_NXCELLS+p7X_SCALE]:1.0f;
      adj5=(i+5<=L)?adj4/fwd->xmx[(i+4)*p7X_NXCELLS+p7X_SCALE]:1.0f;

      register __m128 adj2v=_mm_set1_ps(adj2);
      register __m128 adj3v=_mm_set1_ps(adj3);
      register __m128 adj4v=_mm_set1_ps(adj4);
      register __m128 adj5v=_mm_set1_ps(adj5);

      for (q=0;q<Q;q++)
        ivxf[q]=_mm_add_ps(
                  _mm_add_ps(
                    _mm_add_ps(_mm_mul_ps(MMO(dpp1,q),om_fs->rfv[c1][q]),
                               _mm_mul_ps(_mm_mul_ps(MMO(dpp2,q),adj2v),om_fs->rfv[c2][q])),
                    _mm_add_ps(_mm_mul_ps(_mm_mul_ps(MMO(dpp3,q),adj3v),om_fs->rfv[c3][q]),
                               _mm_mul_ps(_mm_mul_ps(MMO(dpp4,q),adj4v),om_fs->rfv[c4][q]))),
                  _mm_mul_ps(_mm_mul_ps(MMO(dpp5,q),adj5v),om_fs->rfv[c5][q]));

      xBv=zerov; tp7=om_fs->tfv;
      for (q=0;q<Q;q++) { xBv=_mm_add_ps(xBv,_mm_mul_ps(ivxf[q],*tp7)); tp7+=7; }
      xBv=_mm_add_ps(xBv,_mm_shuffle_ps(xBv,xBv,_MM_SHUFFLE(0,3,2,1)));
      xBv=_mm_add_ps(xBv,_mm_shuffle_ps(xBv,xBv,_MM_SHUFFLE(1,0,3,2)));
      _mm_store_ss(&xB,xBv);

      xC=xC_buf[b3]*om_fs->xf[p7O_C][p7O_LOOP];
      xJ=xJ_buf[b3]*om_fs->xf[p7O_J][p7O_LOOP]+xB*om_fs->xf[p7O_J][p7O_MOVE];
      xN=xN_buf[b3]*om_fs->xf[p7O_N][p7O_LOOP]+xB*om_fs->xf[p7O_N][p7O_MOVE];
      xE=xJ*om_fs->xf[p7O_E][p7O_LOOP]+xC*om_fs->xf[p7O_E][p7O_MOVE];
      xEv=_mm_set1_ps(xE);

      for (q=0;q<Q;q++) { MMO(dpc,q)=xEv; DMO(dpc,q)=xEv; IMO(dpc,q)=zerov; }

      tp7=om_fs->tfv;
      for (q=0;q<Q;q++)
        { MMO(dpc,q)=_mm_add_ps(MMO(dpc,q),_mm_mul_ps(_mm_mul_ps(IMO(dpp3,q),adj3v),*(tp7+5)));
          IMO(dpc,q)=_mm_add_ps(IMO(dpc,q),_mm_mul_ps(_mm_mul_ps(IMO(dpp3,q),adj3v),*(tp7+6))); tp7+=7; }

      ivx_carry=esl_sse_leftshiftz_float(ivxf[0]);
      tmmv=esl_sse_leftshiftz_float(*(om_fs->tfv+1));
      timv=esl_sse_leftshiftz_float(*(om_fs->tfv+2));
      tdmv=esl_sse_leftshiftz_float(*(om_fs->tfv+3));
      tp_mm=om_fs->tfv;
      for (q=Q-1;q>=0;q--)
        { MMO(dpc,q)=_mm_add_ps(MMO(dpc,q),_mm_mul_ps(ivx_carry,tmmv));
          IMO(dpc,q)=_mm_add_ps(IMO(dpc,q),_mm_mul_ps(ivx_carry,timv));
          DMO(dpc,q)=_mm_add_ps(DMO(dpc,q),_mm_mul_ps(ivx_carry,tdmv));
          ivx_carry=ivxf[q];
          tp_mm=(tp_mm==om_fs->tfv)?om_fs->tfv+7*(Q-1):tp_mm-7;
          tmmv=*(tp_mm+1); timv=*(tp_mm+2); tdmv=*(tp_mm+3); }

      dcv=esl_sse_leftshiftz_float(DMO(dpc,0)); tp_dd=om_fs->tfv+8*Q-1;
      for (q=Q-1;q>=0;q--)
        { sv=_mm_mul_ps(dcv,*tp_dd); tp_dd--; DMO(dpc,q)=_mm_add_ps(DMO(dpc,q),sv); dcv=DMO(dpc,q); }
      dcv=sv;
      for (j=1;j<4;j++)
        { dcv=esl_sse_leftshiftz_float(dcv); tp_dd=om_fs->tfv+8*Q-1;
          for (q=Q-1;q>=0;q--) { sv=_mm_mul_ps(dcv,*tp_dd); tp_dd--; DMO(dpc,q)=_mm_add_ps(DMO(dpc,q),sv); dcv=_mm_mul_ps(dcv,*(tp_dd+1)); } }

      dcv=esl_sse_leftshiftz_float(DMO(dpc,0)); tp7=om_fs->tfv+7*Q-3;
      for (q=Q-1;q>=0;q--) { MMO(dpc,q)=_mm_add_ps(MMO(dpc,q),_mm_mul_ps(dcv,*tp7)); tp7-=7; dcv=DMO(dpc,q); }

      if (bck->has_own_scales) scale=(xB>1.0e4f)?xB:1.0f;
      else                     scale=fwd->xmx[i*p7X_NXCELLS+p7X_SCALE];
      if (xB>1.0e16f) bck->has_own_scales=TRUE;

      bck->xmx[i*p7X_NXCELLS+p7X_SCALE]=scale;
      if (scale>1.0f)
        { float sf=1.0f/scale; xN*=sf; xJ*=sf; xC*=sf; xB*=sf; xE*=sf; xEv=_mm_set1_ps(sf);
          for (q=0;q<Q;q++) { MMO(dpc,q)=_mm_mul_ps(MMO(dpc,q),xEv); DMO(dpc,q)=_mm_mul_ps(DMO(dpc,q),xEv); IMO(dpc,q)=_mm_mul_ps(IMO(dpc,q),xEv); }
          for (r=0;r<PARSER_ROWS_BWD;r++) { xN_buf[r]*=sf; xB_buf[r]*=sf; xJ_buf[r]*=sf; xC_buf[r]*=sf; }
          bck->totscale+=log(scale); }

      xN_buf[b]=xN; xB_buf[b]=xB; xJ_buf[b]=xJ; xC_buf[b]=xC;
      bck->xmx[i*p7X_NXCELLS+p7X_E]=xE; bck->xmx[i*p7X_NXCELLS+p7X_N]=xN;
      bck->xmx[i*p7X_NXCELLS+p7X_J]=xJ; bck->xmx[i*p7X_NXCELLS+p7X_B]=xB;
      bck->xmx[i*p7X_NXCELLS+p7X_C]=xC;
    }

  /* Termination: i=0 */
  if (dsq[1]<p7P_MAXNUC) x=dsq[1]; else x=p7P_MAXCODONS5;
  if (L>=2 && dsq[2]<p7P_MAXNUC) w=dsq[2]; else w=p7P_MAXCODONS5;
  if (L>=3 && dsq[3]<p7P_MAXNUC) v=dsq[3]; else v=p7P_MAXCODONS5;
  if (L>=4 && dsq[4]<p7P_MAXNUC) u=dsq[4]; else u=p7P_MAXCODONS5;
  if (L>=5 && dsq[5]<p7P_MAXNUC) t=dsq[5]; else t=p7P_MAXCODONS5;

  c1=p7P_CODON1_FS5(x);             c1=p7P_MINIDX(c1,p7P_DEGEN5_QC2);
  c2=p7P_CODON2_FS5(x,w);           c2=p7P_MINIDX(c2,p7P_DEGEN5_QC1);
  c3=p7P_CODON3_FS5(x,w,v);         c3=p7P_MINIDX(c3,p7P_DEGEN5_C);
  c4=p7P_CODON4_FS5(x,w,v,u);       c4=p7P_MINIDX(c4,p7P_DEGEN5_QC1);
  c5=p7P_CODON5_FS5(x,w,v,u,t);     c5=p7P_MINIDX(c5,p7P_DEGEN5_QC2);

  dpp1=bck->dpf[1];
  dpp2=(L>=2)?bck->dpf[2]:zero_row;
  dpp3=(L>=3)?bck->dpf[3]:zero_row;
  dpp4=(L>=4)?bck->dpf[4]:zero_row;
  dpp5=(L>=5)?bck->dpf[5]:zero_row;

  adj2=(L>=2)?1.0f/fwd->xmx[1*p7X_NXCELLS+p7X_SCALE]:1.0f;
  adj3=(L>=3)?adj2/fwd->xmx[2*p7X_NXCELLS+p7X_SCALE]:1.0f;
  adj4=(L>=4)?adj3/fwd->xmx[3*p7X_NXCELLS+p7X_SCALE]:1.0f;
  adj5=(L>=5)?adj4/fwd->xmx[4*p7X_NXCELLS+p7X_SCALE]:1.0f;

  {
    register __m128 adj2v=_mm_set1_ps(adj2);
    register __m128 adj3v=_mm_set1_ps(adj3);
    register __m128 adj4v=_mm_set1_ps(adj4);
    register __m128 adj5v=_mm_set1_ps(adj5);

    for (q=0;q<Q;q++)
      ivxf[q]=_mm_add_ps(
                _mm_add_ps(
                  _mm_add_ps(_mm_mul_ps(MMO(dpp1,q),om_fs->rfv[c1][q]),
                             _mm_mul_ps(_mm_mul_ps(MMO(dpp2,q),adj2v),om_fs->rfv[c2][q])),
                  _mm_add_ps(_mm_mul_ps(_mm_mul_ps(MMO(dpp3,q),adj3v),om_fs->rfv[c3][q]),
                             _mm_mul_ps(_mm_mul_ps(MMO(dpp4,q),adj4v),om_fs->rfv[c4][q]))),
                _mm_mul_ps(_mm_mul_ps(MMO(dpp5,q),adj5v),om_fs->rfv[c5][q]));
  }

  xBv=zerov; tp7=om_fs->tfv;
  for (q=0;q<Q;q++) { xBv=_mm_add_ps(xBv,_mm_mul_ps(ivxf[q],*tp7)); tp7+=7; }
  xBv=_mm_add_ps(xBv,_mm_shuffle_ps(xBv,xBv,_MM_SHUFFLE(0,3,2,1)));
  xBv=_mm_add_ps(xBv,_mm_shuffle_ps(xBv,xBv,_MM_SHUFFLE(1,0,3,2)));
  _mm_store_ss(&xB,xBv);

  xN=xN_buf[3%PARSER_ROWS_BWD]*om_fs->xf[p7O_N][p7O_LOOP]+xB*om_fs->xf[p7O_N][p7O_MOVE];

  bck->xmx[p7X_B]=xB; bck->xmx[p7X_N]=xN; bck->xmx[p7X_J]=0.0f;
  bck->xmx[p7X_C]=0.0f; bck->xmx[p7X_E]=0.0f; bck->xmx[p7X_SCALE]=1.0f;

  {
    float xNtot=xN+xN_buf[1%PARSER_ROWS_BWD]+xN_buf[2%PARSER_ROWS_BWD];

    if      (isnan(xNtot))           ESL_EXCEPTION(eslERANGE, "backward score is NaN");
    else if (isinf(xNtot)==1)        ESL_EXCEPTION(eslERANGE, "backward score overflow (is infinity)");
    else if (L>0 && xNtot==0.0f)     { if (opt_sc) *opt_sc=-eslINFINITY; free(zero_mem); return eslERANGE; }
    if (opt_sc) *opt_sc=bck->totscale+logf(xNtot);
  }

  free(zero_mem);
  return eslOK;

 ERROR:
  if (zero_mem) free(zero_mem);
  return status;
}

#endif /* eslENABLE_SSE */
