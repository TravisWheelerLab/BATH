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

static int backward_frameshift_engine_sse(int do_full, const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, P7_OMX *pp, float *opt_sc);

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

  /* C of rows L+1 and L+2, past the end, read by rows L-1 and L-2 */
  xC_buf[(L+1) % PARSER_ROWS_BWD] = om_fs->xf[p7O_C][p7O_MOVE];
  xC_buf[(L+2) % PARSER_ROWS_BWD] = om_fs->xf[p7O_C][p7O_MOVE];

  /* Initialization: rows L and L-1 */
  for (i = L; i >= L-1; i--)
    {
      b    = i % PARSER_ROWS_BWD;
      curr = b;

      if      (i == L)   xC = om_fs->xf[p7O_C][p7O_MOVE];
      else if (i == L-1) xC = xC_buf[(L+2) % PARSER_ROWS_BWD] * om_fs->xf[p7O_C][p7O_LOOP];

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
          for (r = 0; r < PARSER_ROWS_BWD; r++) { xN_buf[r] *= sf; xB_buf[r] *= sf; xJ_buf[r] *= sf; xC_buf[r] *= sf; }
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

  xC = xC_buf[(L+1) % PARSER_ROWS_BWD] * om_fs->xf[p7O_C][p7O_LOOP];
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
      /* rescale the ring too: it holds C of rows L+1 and L+2, read by rows L-1 and L-2 */
      for (r = 0; r < PARSER_ROWS_BWD; r++) { xN_buf[r] *= sf; xB_buf[r] *= sf; xJ_buf[r] *= sf; xC_buf[r] *= sf; }
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
/* Function:  p7_ForwardRow_Frameshift_sse()
 * Synopsis:  One row of the frameshift Forward algorithm.
 *
 * Purpose:   Forward row <i> of the 5-codon profile <om_fs> against the
 *            DNA sequence <dsq> into <dpc>, from rows i-1 (<dpp1>) and
 *            i-3 (<dpp3>, read for i >= 3), the codon partial sums of
 *            rows i-1..i-4 in <ivx> (<p7P_5CODONS> rows of a P7_OIVX)
 *            and the special states of earlier rows in the ring <rg>.
 *            Writes row <i>'s special states and scale factor into
 *            <xmx>; when the row is rescaled, the log of its scale factor
 *            is added to <*totscale>. Row 0 is the initialization, of
 *            <dpc>, <ivx>, <rg> and row 0 of <xmx>. Same calculation as
 *            each row of <p7_Forward_Frameshift_sse()>, which calls it.
 */
void
p7_ForwardRow_Frameshift_sse(const ESL_DSQ *dsq, const P7_FS_OPROFILE *om_fs, float *xmx, float *totscale, P7_FS_FWDRING *rg,
                       __m128 **ivx, __m128 *dpc, const __m128 *dpp1, const __m128 *dpp3, int i)
{
  register __m128 mpv1, dpv1, ipv1;
  register __m128 sv;
  register __m128 msv;
  register __m128 dcv;
  register __m128 xEv;
  register __m128 xBv1;
  __m128   zerov;
  float    xN, xE, xB, xC, xJ;
  int      b, b1, b3;
  int      ivx_1, ivx_2, ivx_3, ivx_4, ivx_5;
  __m128  *tp;
  int      Q = p7O_NQF(om_fs->M);
  int      q, j, r;
  int      c1, c2, c3, c4, c5;
  int      t, u, v, w, x;
  float    insert_adj;


  zerov = _mm_setzero_ps();

  if (i == 0)
    {
      { __m128 *drow = dpc;
        for (q = 0; q < Q; q++)
          MMO_FS(drow,q,p7X_FS_C0) = MMO_FS(drow,q,p7X_FS_C1) =
          MMO_FS(drow,q,p7X_FS_C2) = MMO_FS(drow,q,p7X_FS_C3) =
          MMO_FS(drow,q,p7X_FS_C4) = MMO_FS(drow,q,p7X_FS_C5) =
          DMO_FS(drow,q)            = IMO_FS(drow,q)            = zerov; }

      for (r = 0; r < p7P_5CODONS; r++)
        for (q = 0; q < Q; q++)
          IVX(r, q) = zerov;

      for (r = 0; r < PARSER_ROWS_FWD; r++)
        rg->xN[r] = rg->xB[r] = rg->xJ[r] = rg->xC[r] = 0.0f;
      rg->xN[0] = rg->xN[1] = rg->xN[2] = 1.0f;
      rg->xB[0] = rg->xB[1] = rg->xB[2] = om_fs->xf[p7O_N][p7O_MOVE];

      xmx[p7X_SCALE] = 1.0f;
      xmx[p7X_E]     = 0.0f;
      xmx[p7X_N]     = 1.0f;
      xmx[p7X_J]     = 0.0f;
      xmx[p7X_B]     = om_fs->xf[p7O_N][p7O_MOVE];
      xmx[p7X_C]     = 0.0f;
      return;
    }

  x = (            dsq[i]   < p7P_MAXNUC) ? dsq[i]   : p7P_MAXCODONS5;
  w = (i >= 2 && dsq[i-1] < p7P_MAXNUC) ? dsq[i-1] : p7P_MAXCODONS5;
  v = (i >= 3 && dsq[i-2] < p7P_MAXNUC) ? dsq[i-2] : p7P_MAXCODONS5;
  u = (i >= 4 && dsq[i-3] < p7P_MAXNUC) ? dsq[i-3] : p7P_MAXCODONS5;
  t = (i >= 5 && dsq[i-4] < p7P_MAXNUC) ? dsq[i-4] : p7P_MAXCODONS5;

  if (i == 1)
    {
    c1 = p7P_CODON1_FS5(x); c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
    ivx_1 = 1;
    xBv1 = _mm_set1_ps(rg->xB[0]);
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
          { rg->xN[r]*=scale_factor; rg->xB[r]*=scale_factor; rg->xJ[r]*=scale_factor; rg->xC[r]*=scale_factor; }
        xmx[1*p7X_NXCELLS+p7X_SCALE]=xE; *totscale+=log(xE); xE=1.0f;
      }
    else xmx[1*p7X_NXCELLS+p7X_SCALE]=1.0f;

    rg->xN[1]=xN; rg->xB[1]=xB; rg->xJ[1]=xJ; rg->xC[1]=xC;
    xmx[1*p7X_NXCELLS+p7X_E]=xE; xmx[1*p7X_NXCELLS+p7X_N]=xN;
    xmx[1*p7X_NXCELLS+p7X_J]=xJ; xmx[1*p7X_NXCELLS+p7X_B]=xB;
    xmx[1*p7X_NXCELLS+p7X_C]=xC;
      return;
    }

  if (i == 2)
    {
    c1=p7P_CODON1_FS5(x);    c1=p7P_MINIDX(c1,p7P_DEGEN5_QC2);
    c2=p7P_CODON2_FS5(w,x);  c2=p7P_MINIDX(c2,p7P_DEGEN5_QC1);
    ivx_1=2; ivx_2=1;

    xBv1=_mm_set1_ps(rg->xB[1]); tp=om_fs->tfv; dcv=zerov; xEv=zerov;
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
          { rg->xN[r]*=scale_factor; rg->xB[r]*=scale_factor; rg->xJ[r]*=scale_factor; rg->xC[r]*=scale_factor; }
        xmx[2*p7X_NXCELLS+p7X_SCALE]=xE; *totscale+=log(xE); xE=1.0f;
      }
    else xmx[2*p7X_NXCELLS+p7X_SCALE]=1.0f;

    rg->xN[2]=xN; rg->xB[2]=xB; rg->xJ[2]=xJ; rg->xC[2]=xC;
    xmx[2*p7X_NXCELLS+p7X_E]=xE; xmx[2*p7X_NXCELLS+p7X_N]=xN;
    xmx[2*p7X_NXCELLS+p7X_J]=xJ; xmx[2*p7X_NXCELLS+p7X_B]=xB;
    xmx[2*p7X_NXCELLS+p7X_C]=xC;
      return;
    }

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

  insert_adj = 1.0f
             / (xmx[(i-2)*p7X_NXCELLS+p7X_SCALE]
             *  xmx[(i-1)*p7X_NXCELLS+p7X_SCALE]);

  mpv1=esl_sse_rightshiftz_float(MMO_FS(dpp1,Q-1,p7X_FS_C0));
  dpv1=esl_sse_rightshiftz_float(DMO_FS(dpp1,Q-1));
  ipv1=esl_sse_rightshiftz_float(IMO_FS(dpp1,Q-1));

  xBv1=_mm_set1_ps(rg->xB[b1]);
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

  xN=rg->xN[b3]*om_fs->xf[p7O_N][p7O_LOOP];
  xJ=rg->xJ[b3]*om_fs->xf[p7O_J][p7O_LOOP]+xE*om_fs->xf[p7O_E][p7O_LOOP];
  xC=rg->xC[b3]*om_fs->xf[p7O_C][p7O_LOOP]+xE*om_fs->xf[p7O_E][p7O_MOVE];
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
        { rg->xN[r]*=scale_factor; rg->xB[r]*=scale_factor; rg->xJ[r]*=scale_factor; rg->xC[r]*=scale_factor; }
      xmx[i*p7X_NXCELLS+p7X_SCALE]=xE; *totscale+=log(xE); xE=1.0f;
    }
  else xmx[i*p7X_NXCELLS+p7X_SCALE]=1.0f;

  rg->xN[b]=xN; rg->xB[b]=xB; rg->xJ[b]=xJ; rg->xC[b]=xC;
  xmx[i*p7X_NXCELLS+p7X_E]=xE; xmx[i*p7X_NXCELLS+p7X_N]=xN;
  xmx[i*p7X_NXCELLS+p7X_J]=xJ; xmx[i*p7X_NXCELLS+p7X_B]=xB;
  xmx[i*p7X_NXCELLS+p7X_C]=xC;
}


int
p7_Forward_Frameshift_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_OMX *fwd, P7_OIVX *ov, float *opt_sc)
{
  P7_FS_FWDRING rg;
  int           i;

  if (om_fs->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");

  fwd->M              = om_fs->M;
  fwd->L              = L;
  fwd->has_own_scales = TRUE;
  fwd->totscale       = 0.0;

  p7_ForwardRow_Frameshift_sse(dsq, om_fs, fwd->xmx, &fwd->totscale, &rg, ov->ivx, fwd->dpf[0], NULL, NULL, 0);
  for (i = 1; i <= L; i++)
    p7_ForwardRow_Frameshift_sse(dsq, om_fs, fwd->xmx, &fwd->totscale, &rg, ov->ivx, fwd->dpf[i], fwd->dpf[i-1], (i >= 3) ? fwd->dpf[i-3] : NULL, i);

  {
    float xCL   = rg.xC[   L    %PARSER_ROWS_FWD];
    float xCLm1 = rg.xC[((L-1)%PARSER_ROWS_FWD+PARSER_ROWS_FWD)%PARSER_ROWS_FWD];
    float xCLm2 = rg.xC[((L-2)%PARSER_ROWS_FWD+PARSER_ROWS_FWD)%PARSER_ROWS_FWD];
    float xCtot = xCL + xCLm1*om_fs->xf[p7O_C][p7O_LOOP] + xCLm2*om_fs->xf[p7O_C][p7O_LOOP];

    if      (isnan(xCtot))           ESL_EXCEPTION(eslERANGE, "forward score is NaN");
    else if (isinf(xCtot)==1)        ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");
    else if (L>1 && xCtot==0.0f)     { if (opt_sc) *opt_sc=-eslINFINITY; return eslERANGE; }
    if (opt_sc) *opt_sc=fwd->totscale+logf(xCtot*om_fs->xf[p7O_C][p7O_MOVE]);
  }

  return eslOK;
}


/* fwd_fs_ring_load()
 * Ring entries before row <b> (> 0) is computed again: rows b-4 .. b-1
 * (from 0), from the special states stored in <xmx>, each multiplied by
 * the inverse scale factors of the rows computed after it up to b-1, in
 * the same order as when the rows were first computed, so the entries
 * are identical.
 */
static void
fwd_fs_ring_load(P7_FS_FWDRING *rg, const float *xmx, int b)
{
  float sf;
  int   j, r, s;

  for (r = 0; r < PARSER_ROWS_FWD; r++)
    rg->xN[r] = rg->xB[r] = rg->xJ[r] = rg->xC[r] = 0.0f;
  for (j = ESL_MAX(0, b-PARSER_ROWS_FWD); j < b; j++)
    {
      s = j % PARSER_ROWS_FWD;
      rg->xN[s] = xmx[j*p7X_NXCELLS+p7X_N];
      rg->xB[s] = xmx[j*p7X_NXCELLS+p7X_B];
      rg->xJ[s] = xmx[j*p7X_NXCELLS+p7X_J];
      rg->xC[s] = xmx[j*p7X_NXCELLS+p7X_C];
      for (r = j+1; r < b; r++)
        if (xmx[r*p7X_NXCELLS+p7X_SCALE] > 1.0f)
          {
            sf = 1.0f / xmx[r*p7X_NXCELLS+p7X_SCALE];
            rg->xN[s] *= sf; rg->xB[s] *= sf; rg->xJ[s] *= sf; rg->xC[s] *= sf;
          }
    }
}


/* forward_ivx_load_sse()
 * Codon partial sums before row <b> (> 0) is computed again: rows
 * b-4 .. b-1 (from 1) into <ivx>, each from the row before it (core
 * cells in <dp>, xB in <xmx>) as p7_ForwardRow_Frameshift_sse() computes
 * it, then multiplied by the inverse scale factors of its own row and
 * the rows after it up to b-1, in the same order as when the rows were
 * first computed, so the values are identical. The other slots are set
 * to 0.
 */
static void
forward_ivx_load_sse(const P7_FS_OPROFILE *om_fs, __m128 * const *dp, const float *xmx, __m128 **ivx, int b)
{
  __m128   zerov = _mm_setzero_ps();
  __m128   mpv, dpv, ipv, xBv, sv, sfv;
  __m128  *tp;
  __m128  *dpp;
  int  Q = p7O_NQF(om_fs->M);
  int  q, r, j, slot;

  for (r = 0; r < p7P_5CODONS; r++)
    for (q = 0; q < Q; q++)
      ivx[r][q] = zerov;

  for (r = ESL_MAX(1, b-4); r < b; r++)
    {
      slot = r % p7P_5CODONS;
      dpp  = dp[r-1];
      xBv  = _mm_set1_ps(xmx[(r-1)*p7X_NXCELLS+p7X_B]);
      mpv  = esl_sse_rightshiftz_float(MMO_FS(dpp,Q-1,p7X_FS_C0));
      dpv  = esl_sse_rightshiftz_float(DMO_FS(dpp,Q-1));
      ipv  = esl_sse_rightshiftz_float(IMO_FS(dpp,Q-1));
      tp   = om_fs->tfv;
      for (q = 0; q < Q; q++)
        {
          sv  =            _mm_mul_ps(xBv, *tp);  tp++;
          sv  = _mm_add_ps(sv,   _mm_mul_ps(mpv, *tp)); tp++;
          sv  = _mm_add_ps(sv,   _mm_mul_ps(ipv, *tp)); tp++;
          sv  = _mm_add_ps(sv,   _mm_mul_ps(dpv, *tp)); tp++;
          tp += 3;
          ivx[slot][q] = sv;
          mpv = MMO_FS(dpp,q,p7X_FS_C0); dpv = DMO_FS(dpp,q); ipv = IMO_FS(dpp,q);
        }
      for (j = r; j < b; j++)
        if (xmx[j*p7X_NXCELLS+p7X_SCALE] > 1.0f)
          {
            sfv = _mm_set1_ps(1.0f / xmx[j*p7X_NXCELLS+p7X_SCALE]);
            for (q = 0; q < Q; q++) ivx[slot][q] = _mm_mul_ps(ivx[slot][q], sfv);
          }
    }
}


/* Function:  p7_ForwardChk_Frameshift_sse()
 * Synopsis:  Frameshift Forward in a checkpointed matrix.
 *
 * Purpose:   The Forward algorithm of <p7_Forward_Frameshift_sse()> for the
 *            5-codon profile <om_fs> against <dsq> of length <L>, computed
 *            in the checkpointed matrix <cx>, with <ov> as workspace. Each
 *            row is computed by <p7_ForwardRow_Frameshift_sse()>, so cells,
 *            special states, scale factors and score are identical to
 *            <p7_Forward_Frameshift_sse()>.
 *
 *            Rows are computed from 0 up to <L>, one block at a time;
 *            checkpoints of the 5 rows before each block are saved, so
 *            <p7_ForwardChkLoad_Frameshift_sse()> can compute any block
 *            again, with identical values. On return the block holding
 *            row <L> is loaded and the special states and scale factors
 *            of all rows are in <cx->xmx>.
 *
 *            The caller creates <cx> with <p7_CHK_FORWARD> (or sets that
 *            direction with <p7_chkomx_SetDirection()>), <p7X_NSCELLS_FS>
 *            vectors per stripe and at least 5 saved rows, and grows it
 *            for <om_fs->M> by <L> with the block size it wants
 *            (<p7_chkomx_GrowTo()>).
 *
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq (>= 1)
 *            om_fs  - 5-codon frameshift profile
 *            cx     - checkpointed DP matrix
 *            ov     - workspace for the profile's M
 *            opt_sc - optRETURN: Forward lod score in nats
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if the score underflows to 0 (L > 1); <*opt_sc>
 *            is then -infinity.
 *
 * Throws:    <eslEINVAL> if <cx> isn't set up for this computation or
 *            <om_fs> isn't a 5-codon profile.
 *            <eslERANGE> if the score is NaN or overflows.
 */
int
p7_ForwardChk_Frameshift_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, float *opt_sc)
{
  P7_FS_FWDRING rg;
  float         totscale = 0.0f;
  float         xCL, xCLm1, xCLm2, xCtot;
  int           b, i;
  int           status;

  if (om_fs->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");
  if (cx->dir != p7_CHK_FORWARD || cx->saveR < 5 || cx->nscells != p7X_NSCELLS_FS)
    ESL_EXCEPTION(eslEINVAL, "matrix not set up for the checkpointed frameshift Forward");
  if (L < 1) ESL_EXCEPTION(eslEINVAL, "sequence length must be >= 1");

  cx->M = om_fs->M;
  cx->L = L;

  for (b = 0; b <= L; b += cx->blockR)
    {
      if ((status = p7_chkomx_StartBlock_sse(cx, b)) != eslOK) return status;
      for (i = b; i <= ESL_MIN(L, b + cx->blockR - 1); i++)
        p7_ForwardRow_Frameshift_sse(dsq, om_fs, cx->xmx, &totscale, &rg, ov->ivx, cx->dpf[i],
                                     (i >= 1) ? cx->dpf[i-1] : NULL, (i >= 3) ? cx->dpf[i-3] : NULL, i);
    }

  xCL   = rg.xC[   L    %PARSER_ROWS_FWD];
  xCLm1 = rg.xC[((L-1)%PARSER_ROWS_FWD+PARSER_ROWS_FWD)%PARSER_ROWS_FWD];
  xCLm2 = rg.xC[((L-2)%PARSER_ROWS_FWD+PARSER_ROWS_FWD)%PARSER_ROWS_FWD];
  xCtot = xCL + xCLm1*om_fs->xf[p7O_C][p7O_LOOP] + xCLm2*om_fs->xf[p7O_C][p7O_LOOP];

  if      (isnan(xCtot))         ESL_EXCEPTION(eslERANGE, "forward score is NaN");
  else if (isinf(xCtot) == 1)    ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");
  else if (L > 1 && xCtot == 0.0f) { if (opt_sc) *opt_sc = -eslINFINITY; return eslERANGE; }
  if (opt_sc != NULL) *opt_sc = totscale + logf(xCtot * om_fs->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}


/* Function:  p7_ForwardChkLoad_Frameshift_sse()
 * Synopsis:  Compute one block of a checkpointed frameshift Forward again.
 *
 * Purpose:   After <p7_ForwardChk_Frameshift_sse()> has filled <cx> for
 *            <om_fs> and <dsq>, compute block <b> (a multiple of
 *            <cx->blockR>) again from the checkpoint before it, with <ov>
 *            as workspace: the codon partial sums of the 4 rows before
 *            the block are computed from the saved rows and the special
 *            states of the rows before it are taken from <cx->xmx>
 *            first. Values are identical to those of the fill. Blocks are
 *            loaded in decreasing order, as a computation that runs
 *            backward through the rows reads them
 *            (<p7_chkomx_LoadBlock()>).
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if block <b> can't be loaded.
 */
int
p7_ForwardChkLoad_Frameshift_sse(const ESL_DSQ *dsq, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, int b)
{
  P7_FS_FWDRING rg;
  float         totscale = 0.0f;   /* unused: the scale factors are already in <cx->xmx> */
  int           i;
  int           status;

  if ((status = p7_chkomx_LoadBlock_sse(cx, b)) != eslOK) return status;
  if (b > 0)
    {
      fwd_fs_ring_load(&rg, cx->xmx, b);
      forward_ivx_load_sse(om_fs, cx->dpf, cx->xmx, ov->ivx, b);
    }
  for (i = b; i <= ESL_MIN(cx->L, b + cx->blockR - 1); i++)
    p7_ForwardRow_Frameshift_sse(dsq, om_fs, cx->xmx, &totscale, &rg, ov->ivx, cx->dpf[i],
                                 (i >= 1) ? cx->dpf[i-1] : NULL, (i >= 3) ? cx->dpf[i-3] : NULL, i);
  return eslOK;
}


/* Function:  p7_Backward_Frameshift_sse()
 * Synopsis:  SSE frameshift-aware Backward, 5 codon lengths, full O(ML) matrix.
 */
int
p7_Backward_Frameshift_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
  return backward_frameshift_engine_sse(TRUE, dsq, L, om_fs, fwd, bck, ov, NULL, opt_sc);
}


/* Function:  p7_BackwardDecoding_Frameshift_sse()
 * Synopsis:  Frameshift-aware Backward in a ring of rows, with posterior decoding.
 *
 * Purpose:   The Backward calculation of p7_Backward_Frameshift_sse(), kept
 *            in a ring of <PARSER_ROWS_BWD> rows of <bck> (plus special
 *            states for all rows), writing the posterior decoding over
 *            the Forward matrix <fwd> as the rows are finished from <L>
 *            down. On return <fwd> holds the posterior decoding matrix
 *            that p7_Decoding_Frameshift_sse() would produce, and <bck> holds
 *            the Backward special states and scale factors for all rows.
 *            The caller allocates <bck> with
 *            <p7_omx_GrowTo_dpf(bck, M, PARSER_ROWS_BWD, L)>.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if the Backward score or the posterior scaling
 *            is out of range, as for p7_Backward_Frameshift_sse() and
 *            p7_Decoding_Frameshift_sse().
 */
int
p7_BackwardDecoding_Frameshift_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
#if eslDEBUGLEVEL > 0
  if (bck->validR < PARSER_ROWS_BWD) ESL_EXCEPTION(eslEINVAL, "DP matrix needs at least PARSER_ROWS_BWD MDI rows");
  if (L >= bck->allocXR)             ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few X rows)");
  if (L >= fwd->validR)              ESL_EXCEPTION(eslEINVAL, "fwd matrix isn't a full matrix (too few MDI rows)");
#endif
  return backward_frameshift_engine_sse(FALSE, dsq, L, om_fs, fwd, bck, ov, fwd, opt_sc);
}


/* backward_row_sse()
 * Row <r> of <bck> for backward_frameshift_engine_sse(): row <r> of a full
 * matrix, or ring row <r> % PARSER_ROWS_BWD. The engine reads rows i+1..i+5
 * while computing row i and never changes a row once it is stored, so
 * the ring holds the same values as the full matrix.
 */
static inline __m128 *
backward_row_sse(const P7_OMX *bck, int do_full, int r)
{
  return bck->dpf[do_full ? r : r % PARSER_ROWS_BWD];
}


/* backward_decode_row_frameshift_sse()
 * For backward_frameshift_engine_sse() with a posterior matrix <pp>, which
 * is the Forward matrix: as soon as Backward row <i> (<dpc>, and its
 * special states in <bck>) is final, write F*B over row <i> of <pp> for
 * the M cells of every codon length and for I, 0 for D, and
 * F(i-3)*B(i)*t_loop for N, J and C (B(i) alone for N when i <= 2). The
 * engine multiplies them by the row scale factors after the sweep,
 * completing the same products, in the same order, as
 * p7_Decoding_Frameshift_sse(). Posterior row <i> reads Forward row <i>
 * and the Forward special states of row <i-3>, which are not overwritten
 * yet when the rows are done from L down.
 */
static inline void
backward_decode_row_frameshift_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *bck, const __m128 *dpc, P7_OMX *pp, int i, int Q)
{
  __m128 *ppv = pp->dpf[i];
  __m128  bM, bI;
  int     q, c;

  for (q = 0; q < Q; q++)
    {
      bM = dpc[q * p7X_NSCELLS + p7X_M];
      bI = dpc[q * p7X_NSCELLS + p7X_I];

      ppv[q * p7X_NSCELLS_FS + p7X_FS_D] = _mm_setzero_ps();
      ppv[q * p7X_NSCELLS_FS + p7X_FS_I] = _mm_mul_ps(ppv[q * p7X_NSCELLS_FS + p7X_FS_I], bI);
      for (c = 0; c < 6; c++)
        ppv[q * p7X_NSCELLS_FS + p7X_FS_M + c] = _mm_mul_ps(ppv[q * p7X_NSCELLS_FS + p7X_FS_M + c], bM);
    }

  pp->xmx[i*p7X_NXCELLS + p7X_E] = 0.0f;
  pp->xmx[i*p7X_NXCELLS + p7X_B] = 0.0f;
  if (i > 2)
    {
      pp->xmx[i*p7X_NXCELLS + p7X_N] = pp->xmx[(i-3)*p7X_NXCELLS + p7X_N] * bck->xmx[i*p7X_NXCELLS + p7X_N] * om_fs->xf[p7O_N][p7O_LOOP];
      pp->xmx[i*p7X_NXCELLS + p7X_J] = pp->xmx[(i-3)*p7X_NXCELLS + p7X_J] * bck->xmx[i*p7X_NXCELLS + p7X_J] * om_fs->xf[p7O_J][p7O_LOOP];
      pp->xmx[i*p7X_NXCELLS + p7X_C] = pp->xmx[(i-3)*p7X_NXCELLS + p7X_C] * bck->xmx[i*p7X_NXCELLS + p7X_C] * om_fs->xf[p7O_C][p7O_LOOP];
    }
  else
    {
      pp->xmx[i*p7X_NXCELLS + p7X_N] = bck->xmx[i*p7X_NXCELLS + p7X_N];
      pp->xmx[i*p7X_NXCELLS + p7X_J] = 0.0f;
      pp->xmx[i*p7X_NXCELLS + p7X_C] = 0.0f;
    }
}


/* Scale factor a Backward row takes, for backward_frameshift_row_L_sse()
 * and backward_frameshift_row_sse():
 *   BWD_SCALE_FORWARD: the Forward matrix's factor for the row, or the
 *     row's own factor once the Backward has turned own scaling on
 *     (p7_Backward_Frameshift_sse());
 *   BWD_SCALE_OWN:     the row's own factor, xB when xB > 1e4, else 1;
 *   BWD_SCALE_STORED:  the factor already stored for the row, when a
 *     row is computed again.
 */
enum bwd_scale_e { BWD_SCALE_FORWARD = 0, BWD_SCALE_OWN = 1, BWD_SCALE_STORED = 2 };

/* BWD_FS_RING
 * The N, B, J and C special states of the last PARSER_ROWS_BWD Backward
 * rows, indexed by row % PARSER_ROWS_BWD. Row i reads those of row i+3.
 * They are kept at the scale of the most recent row: when a row is
 * rescaled, every entry is rescaled with it.
 */
typedef struct {
  float xN[PARSER_ROWS_BWD];
  float xB[PARSER_ROWS_BWD];
  float xJ[PARSER_ROWS_BWD];
  float xC[PARSER_ROWS_BWD];
} BWD_FS_RING;


/* bwd_fs_ring_init()
 * Ring entries before row L is computed: 0, except C for the rows L+1
 * and L+2 past the end of the sequence, which rows L-1 and L-2 read.
 */
static void
bwd_fs_ring_init(BWD_FS_RING *rg, const P7_FS_OPROFILE *om_fs, int L)
{
  int r;
  for (r=0;r<PARSER_ROWS_BWD;r++)
    rg->xN[r]=rg->xB[r]=rg->xJ[r]=rg->xC[r]=0.0f;

  rg->xC[(L+1)%PARSER_ROWS_BWD]=om_fs->xf[p7O_C][p7O_MOVE];
  rg->xC[(L+2)%PARSER_ROWS_BWD]=om_fs->xf[p7O_C][p7O_MOVE];
}


/* backward_frameshift_row_L_sse()
 * Frameshift Backward row <L> into <dpc>; special states and scale factor
 * into row <L> of <xmx> and into the ring <rg>. <fwd_xmx> holds the Forward
 * special states and scale factors (needed only for BWD_SCALE_FORWARD).
 * If <totscale> is non-NULL the log of a scale factor > 1 is added to it.
 */
static void
backward_frameshift_row_L_sse(const P7_FS_OPROFILE *om_fs, const float *fwd_xmx, float *xmx, float *totscale, BWD_FS_RING *rg, __m128 *dpc, int L, int scale_mode)
{
  register __m128 sv;
  register __m128 dcv;
  register __m128 xEv;
  __m128   zerov;
  float    xN, xE, xB, xC, xJ;
  int      b;
  __m128  *tp7, *tp_dd;
  int      Q = p7O_NQF(om_fs->M);
  int      q, j, r;
  float    scale;

  zerov=_mm_setzero_ps();
  b=L%PARSER_ROWS_BWD;
  xC=om_fs->xf[p7O_C][p7O_MOVE]; xN=xB=xJ=0.0f;
  xE=xC*om_fs->xf[p7O_E][p7O_MOVE]; xEv=_mm_set1_ps(xE);
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

  if      (scale_mode == BWD_SCALE_FORWARD) scale=fwd_xmx[L*p7X_NXCELLS+p7X_SCALE];
  else if (scale_mode == BWD_SCALE_OWN)     scale=(xB>1.0e4f)?xB:1.0f;
  else                                      scale=xmx[L*p7X_NXCELLS+p7X_SCALE];
  xmx[L*p7X_NXCELLS+p7X_SCALE]=scale;
  if (scale>1.0f)
    { float sf=1.0f/scale; xN*=sf; xJ*=sf; xC*=sf; xB*=sf; xE*=sf; xEv=_mm_set1_ps(sf);
      for (q=0;q<Q;q++) { MMO(dpc,q)=_mm_mul_ps(MMO(dpc,q),xEv); DMO(dpc,q)=_mm_mul_ps(DMO(dpc,q),xEv); IMO(dpc,q)=_mm_mul_ps(IMO(dpc,q),xEv); }
      /* rescale the ring too: it holds C of rows L+1 and L+2, read by rows L-1 and L-2 */
      for (r = 0; r < PARSER_ROWS_BWD; r++) { rg->xN[r] *= sf; rg->xB[r] *= sf; rg->xJ[r] *= sf; rg->xC[r] *= sf; }
      if (totscale != NULL) *totscale+=log(scale); }

  rg->xN[b]=xN; rg->xB[b]=xB; rg->xJ[b]=xJ; rg->xC[b]=xC;
  xmx[L*p7X_NXCELLS+p7X_E]=xE; xmx[L*p7X_NXCELLS+p7X_N]=xN;
  xmx[L*p7X_NXCELLS+p7X_J]=xJ; xmx[L*p7X_NXCELLS+p7X_B]=xB;
  xmx[L*p7X_NXCELLS+p7X_C]=xC;
}


/* backward_frameshift_row_sse()
 * Frameshift Backward row <i> (1 <= i < L) of a sequence <dsq> of length
 * <L>. <dp>[0] is row <i>, written here; <dp>[d] is row <i>+d for d =
 * 1..5, or a row of zeros past <L>. The look-back rows are corrected
 * by the scale factors stored in <xmx> for the rows between; row <i>+3's
 * special states come from the ring <rg>. Writes row <i>'s special
 * states and scale factor into <xmx> and <rg>, with <ivxf> as workspace.
 * <fwd_xmx> holds the Forward special states and scale factors (needed
 * only for BWD_SCALE_FORWARD). If <totscale> is non-NULL the log of a
 * scale factor > 1 is added to it. <has_own_scales> is the flag that
 * BWD_SCALE_FORWARD turns on once xB exceeds 1e16.
 */
static void
backward_frameshift_row_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const float *fwd_xmx, float *xmx, float *totscale, int *has_own_scales,
                        BWD_FS_RING *rg, __m128 *const *dp, __m128 *ivxf, int i, int scale_mode)
{
  register __m128 sv;
  register __m128 dcv;
  register __m128 ivx_carry;
  __m128   tmmv, timv, tdmv;
  register __m128 xBv;
  register __m128 xEv;
  __m128   zerov;
  float    xN, xE, xB, xC, xJ;
  int      b, b3;
  __m128  *dpc, *dpp1, *dpp2, *dpp3, *dpp4, *dpp5;
  __m128  *tp7, *tp_dd, *tp_mm;
  int      Q = p7O_NQF(om_fs->M);
  int      q, j, r;
  int      c1, c2, c3, c4, c5;
  int      t, u, v, w, x;
  float    scale;
  float    adj2, adj3, adj4, adj5;

  zerov=_mm_setzero_ps();
  x = (            dsq[i+1]<p7P_MAXNUC) ? dsq[i+1] : p7P_MAXCODONS5;
  w = (i+2<=L && dsq[i+2]<p7P_MAXNUC) ? dsq[i+2] : p7P_MAXCODONS5;
  v = (i+3<=L && dsq[i+3]<p7P_MAXNUC) ? dsq[i+3] : p7P_MAXCODONS5;
  u = (i+4<=L && dsq[i+4]<p7P_MAXNUC) ? dsq[i+4] : p7P_MAXCODONS5;
  t = (i+5<=L && dsq[i+5]<p7P_MAXNUC) ? dsq[i+5] : p7P_MAXCODONS5;

  c1=p7P_CODON1_FS5(x);              c1=p7P_MINIDX(c1,p7P_DEGEN5_QC2);
  c2=p7P_CODON2_FS5(x,w);            c2=p7P_MINIDX(c2,p7P_DEGEN5_QC1);
  c3=p7P_CODON3_FS5(x,w,v);          c3=p7P_MINIDX(c3,p7P_DEGEN5_C);
  c4=p7P_CODON4_FS5(x,w,v,u);        c4=p7P_MINIDX(c4,p7P_DEGEN5_QC1);
  c5=p7P_CODON5_FS5(x,w,v,u,t);      c5=p7P_MINIDX(c5,p7P_DEGEN5_QC2);

  b=i%PARSER_ROWS_BWD; b3=(i+3)%PARSER_ROWS_BWD;
  dpc=dp[0]; dpp1=dp[1]; dpp2=dp[2]; dpp3=dp[3]; dpp4=dp[4]; dpp5=dp[5];

  /* Row i+d (d=2..5) was finished before rows i+1..i+d-1 were scaled, so it must also be divided
   * by their scale factors. Use the factors Backward actually applied (xmx), which differ
   * from Forward's once has_own_scales is set; otherwise they are the same values. */
  adj2=(i+2<=L)?1.0f/xmx[(i+1)*p7X_NXCELLS+p7X_SCALE]:1.0f;
  adj3=(i+3<=L)?adj2/xmx[(i+2)*p7X_NXCELLS+p7X_SCALE]:1.0f;
  adj4=(i+4<=L)?adj3/xmx[(i+3)*p7X_NXCELLS+p7X_SCALE]:1.0f;
  adj5=(i+5<=L)?adj4/xmx[(i+4)*p7X_NXCELLS+p7X_SCALE]:1.0f;

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

  xC=rg->xC[b3]*om_fs->xf[p7O_C][p7O_LOOP];
  xJ=rg->xJ[b3]*om_fs->xf[p7O_J][p7O_LOOP]+xB*om_fs->xf[p7O_J][p7O_MOVE];
  xN=rg->xN[b3]*om_fs->xf[p7O_N][p7O_LOOP]+xB*om_fs->xf[p7O_N][p7O_MOVE];
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
  if (scale_mode == BWD_SCALE_FORWARD)
    {
      if (*has_own_scales) scale=(xB>1.0e4f)?xB:1.0f;
      else                 scale=fwd_xmx[i*p7X_NXCELLS+p7X_SCALE];
      if (xB>1.0e16f) *has_own_scales=TRUE;
    }
  else if (scale_mode == BWD_SCALE_OWN) scale=(xB>1.0e4f)?xB:1.0f;
  else                                  scale=xmx[i*p7X_NXCELLS+p7X_SCALE];
  xmx[i*p7X_NXCELLS+p7X_SCALE]=scale;
  if (scale>1.0f)
    { float sf=1.0f/scale; xN*=sf; xJ*=sf; xC*=sf; xB*=sf; xE*=sf; xEv=_mm_set1_ps(sf);
      for (q=0;q<Q;q++) { MMO(dpc,q)=_mm_mul_ps(MMO(dpc,q),xEv); DMO(dpc,q)=_mm_mul_ps(DMO(dpc,q),xEv); IMO(dpc,q)=_mm_mul_ps(IMO(dpc,q),xEv); }
      for (r=0;r<PARSER_ROWS_BWD;r++) { rg->xN[r]*=sf; rg->xB[r]*=sf; rg->xJ[r]*=sf; rg->xC[r]*=sf; }
      if (totscale != NULL) *totscale+=log(scale); }

  rg->xN[b]=xN; rg->xB[b]=xB; rg->xJ[b]=xJ; rg->xC[b]=xC;
  xmx[i*p7X_NXCELLS+p7X_E]=xE; xmx[i*p7X_NXCELLS+p7X_N]=xN;
  xmx[i*p7X_NXCELLS+p7X_J]=xJ; xmx[i*p7X_NXCELLS+p7X_B]=xB;
  xmx[i*p7X_NXCELLS+p7X_C]=xC;
}


/* backward_frameshift_row_0_sse()
 * Frameshift Backward row 0 of a sequence <dsq> of length <L>: only the N
 * and B states, from rows 1..5 in <dp>[1..5] (rows past <L> are rows
 * of zeros) and row 3's N from the ring <rg>, with <ivxf> as workspace.
 * Writes row 0's special states into <xmx> and returns its N value.
 */
static float
backward_frameshift_row_0_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, float *xmx, const BWD_FS_RING *rg, __m128 *const *dp, __m128 *ivxf)
{
  register __m128 xBv;
  __m128   zerov;
  float    xN, xB;
  __m128  *dpp1, *dpp2, *dpp3, *dpp4, *dpp5;
  __m128  *tp7;
  int      Q = p7O_NQF(om_fs->M);
  int      q;
  int      c1, c2, c3, c4, c5;
  int      t, u, v, w, x;
  float    adj2, adj3, adj4, adj5;

  zerov=_mm_setzero_ps();
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
  dpp1=dp[1]; dpp2=dp[2]; dpp3=dp[3]; dpp4=dp[4]; dpp5=dp[5];

  adj2=(L>=2)?1.0f/xmx[1*p7X_NXCELLS+p7X_SCALE]:1.0f;   /* bck's own scale factors, as in the main loop */
  adj3=(L>=3)?adj2/xmx[2*p7X_NXCELLS+p7X_SCALE]:1.0f;
  adj4=(L>=4)?adj3/xmx[3*p7X_NXCELLS+p7X_SCALE]:1.0f;
  adj5=(L>=5)?adj4/xmx[4*p7X_NXCELLS+p7X_SCALE]:1.0f;

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

  xN=rg->xN[3%PARSER_ROWS_BWD]*om_fs->xf[p7O_N][p7O_LOOP]+xB*om_fs->xf[p7O_N][p7O_MOVE];

  xmx[p7X_B]=xB; xmx[p7X_N]=xN; xmx[p7X_J]=0.0f;
  xmx[p7X_C]=0.0f; xmx[p7X_E]=0.0f; xmx[p7X_SCALE]=1.0f;
  return xN;
}


/* backward_frameshift_engine_sse()
 * The Backward calculation of p7_Backward_Frameshift_sse() (<do_full> TRUE,
 * <pp> NULL) and p7_BackwardDecoding_Frameshift_sse() (<do_full> FALSE, rows
 * kept in a ring; <pp> is the Forward matrix, which receives the
 * posterior decoding).
 */
static int
backward_frameshift_engine_sse(int do_full, const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, P7_OMX *pp, float *opt_sc)
{
  __m128   zerov;
  float    xN;
  __m128  *zero_row = NULL;
  __m128  *zero_mem = NULL;
  __m128  *ivxf;
  int      Q = p7O_NQF(om_fs->M);
  int      i, q;
  int      status;
  float   *log_sfwd = NULL;
  float   *log_sbck = NULL;
  BWD_FS_RING rg;
  __m128  *dp[6];
  int      d;

  if (om_fs->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");

  ivxf = ov->ivx[0];

  bck->M=om_fs->M; bck->L=L; bck->has_own_scales=FALSE; bck->totscale=0.0;
  zerov=_mm_setzero_ps();

  ESL_ALLOC(zero_mem, sizeof(__m128)*Q*p7X_NSCELLS+15);
  zero_row=(__m128*)(((unsigned long int)zero_mem+15)&(~0xf));
  for (q=0;q<Q*p7X_NSCELLS;q++) zero_row[q]=zerov;

  bwd_fs_ring_init(&rg, om_fs, L);

  /* Initialization: row L */
  dp[0]=backward_row_sse(bck, do_full, L);
  backward_frameshift_row_L_sse(om_fs, fwd->xmx, bck->xmx, &bck->totscale, &rg, dp[0], L, BWD_SCALE_FORWARD);
  if (pp != NULL) backward_decode_row_frameshift_sse(om_fs, bck, dp[0], pp, L, Q);

  /* Main recurrence: i = L-1 down to 1 */
  for (i=L-1;i>=1;i--)
    {
      for (d=0;d<6;d++) dp[d]=(i+d<=L)?backward_row_sse(bck, do_full, i+d):zero_row;
      backward_frameshift_row_sse(dsq, L, om_fs, fwd->xmx, bck->xmx, &bck->totscale, &bck->has_own_scales, &rg, dp, ivxf, i, BWD_SCALE_FORWARD);
      if (pp != NULL) backward_decode_row_frameshift_sse(om_fs, bck, dp[0], pp, i, Q);
    }

  /* Termination: i=0 */
  dp[0]=NULL;
  for (d=1;d<6;d++) dp[d]=(d<=L)?backward_row_sse(bck, do_full, d):zero_row;
  xN=backward_frameshift_row_0_sse(dsq, L, om_fs, bck->xmx, &rg, dp, ivxf);

  {
    float xNtot=xN+rg.xN[1%PARSER_ROWS_BWD]+rg.xN[2%PARSER_ROWS_BWD];

    if      (isnan(xNtot))           ESL_EXCEPTION(eslERANGE, "backward score is NaN");
    else if (isinf(xNtot)==1)        ESL_EXCEPTION(eslERANGE, "backward score overflow (is infinity)");
    else if (L>0 && xNtot==0.0f)     { if (opt_sc) *opt_sc=-eslINFINITY; free(zero_mem); return eslERANGE; }
    if (opt_sc) *opt_sc=bck->totscale+logf(xNtot);
  }

  /* Posterior decoding: scale each row as p7_Decoding_Frameshift_sse() does */
  if (pp != NULL)
    {
      float  log_inv_Z, factor_mdi, factor_njc, factor_nsmall;
      __m128 scv;

      ESL_ALLOC(log_sfwd, sizeof(float) * (L + 2));
      ESL_ALLOC(log_sbck, sizeof(float) * (L + 2));

      log_sfwd[0] = logf(fwd->xmx[p7X_SCALE]);
      for (i = 1; i <= L; i++)
        log_sfwd[i] = log_sfwd[i-1] + logf(fwd->xmx[i*p7X_NXCELLS + p7X_SCALE]);

      log_sbck[L+1] = 0.0f;
      for (i = L; i >= 0; i--)
        log_sbck[i] = log_sbck[i+1] + logf(bck->xmx[i*p7X_NXCELLS + p7X_SCALE]);

      log_inv_Z = -p7_FLogsum(
                     logf(bck->xmx[p7X_N]) + log_sbck[0],
                     p7_FLogsum(
                       logf(bck->xmx[1*p7X_NXCELLS + p7X_N]) + log_sbck[1],
                       logf(bck->xmx[2*p7X_NXCELLS + p7X_N]) + log_sbck[2]));

      for (q = 0; q < Q * p7X_NSCELLS_FS; q++) pp->dpf[0][q] = zerov;
      pp->xmx[p7X_E] = pp->xmx[p7X_N] = pp->xmx[p7X_J] = pp->xmx[p7X_B] = pp->xmx[p7X_C] = 0.0f;

      for (i = 1; i <= L; i++)
        {
          factor_mdi = expf(log_sfwd[i] + log_sbck[i] + log_inv_Z);
          if (isinf(factor_mdi)) { status = eslERANGE; goto ERROR; }

          scv = _mm_set1_ps(factor_mdi);
          for (q = 0; q < Q * p7X_NSCELLS_FS; q++)
            pp->dpf[i][q] = _mm_mul_ps(pp->dpf[i][q], scv);

          if (i > 2)
            {
              factor_njc = expf(log_sfwd[i-3] + log_sbck[i] + log_inv_Z);
              pp->xmx[i*p7X_NXCELLS + p7X_N] *= factor_njc;
              pp->xmx[i*p7X_NXCELLS + p7X_J] *= factor_njc;
              pp->xmx[i*p7X_NXCELLS + p7X_C] *= factor_njc;
            }
          else
            {
              factor_nsmall = expf(log_sbck[i] + log_inv_Z);
              pp->xmx[i*p7X_NXCELLS + p7X_N] *= factor_nsmall;
            }
          if (isinf(pp->xmx[i*p7X_NXCELLS + p7X_N] + pp->xmx[i*p7X_NXCELLS + p7X_J] + pp->xmx[i*p7X_NXCELLS + p7X_C])) { status = eslERANGE; goto ERROR; }
        }

      free(log_sfwd);
      free(log_sbck);
    }

  free(zero_mem);
  return eslOK;

 ERROR:
  if (zero_mem) free(zero_mem);
  if (log_sfwd) free(log_sfwd);
  if (log_sbck) free(log_sbck);
  return status;
}


/*****************************************************************
 * Checkpointed frameshift Backward, SSE
 *****************************************************************/

/* bwd_fs_ring_load()
 * Ring entries before row <top> (< L) is computed again: rows <top>+1
 * .. <top>+3, the ones rows <top> .. <top>-2 read, from the special
 * states stored in <xmx> (the starting values of bwd_fs_ring_init()
 * for rows past <L>), each multiplied by the inverse scale factors of
 * the rows computed after it down to <top>+1, in the same order as when
 * the rows were first computed, so the entries are identical.
 */
static void
bwd_fs_ring_load(BWD_FS_RING *rg, const P7_FS_OPROFILE *om_fs, const float *xmx, int L, int top)
{
  float sf;
  int   j, r, s;

  bwd_fs_ring_init(rg, om_fs, L);
  for (j = top+1; j <= top+3; j++)
    {
      s = j % PARSER_ROWS_BWD;
      if (j <= L)
        {
          rg->xN[s] = xmx[j*p7X_NXCELLS+p7X_N];
          rg->xB[s] = xmx[j*p7X_NXCELLS+p7X_B];
          rg->xJ[s] = xmx[j*p7X_NXCELLS+p7X_J];
          rg->xC[s] = xmx[j*p7X_NXCELLS+p7X_C];
        }
      for (r = ESL_MIN(j-1, L); r > top; r--)
        if (xmx[r*p7X_NXCELLS+p7X_SCALE] > 1.0f)
          {
            sf = 1.0f / xmx[r*p7X_NXCELLS+p7X_SCALE];
            rg->xN[s] *= sf; rg->xB[s] *= sf; rg->xJ[s] *= sf; rg->xC[s] *= sf;
          }
    }
}


/* backward_frameshift_chk_block_sse()
 * Compute rows <b> .. <b>+<blockR>-1 (up to <cx->L>) of the loaded
 * block <b> of the checkpointed frameshift Backward matrix <cx>, top
 * row first, with scale factors chosen by <scale_mode>. <rg> must hold
 * the ring entries for the block's top row (bwd_fs_ring_init() for the
 * block holding row L, else the state left by the block above or
 * bwd_fs_ring_load()). <zero_row> is a row of zeros for rows past L;
 * <ivxf> is workspace. For block 0 this includes row 0, whose core
 * cells are set to 0. If <totscale> is non-NULL the logs of the scale
 * factors are added to it. Returns the sum of the N values of rows 0,
 * 1 and 2 at row 0's scale (the Backward score before scaling) for
 * block 0, else 0.
 */
static float
backward_frameshift_chk_block_sse(const ESL_DSQ *dsq, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, BWD_FS_RING *rg,
                                  __m128 *zero_row, __m128 *ivxf, int b, float *totscale, int scale_mode)
{
  __m128 *dp[6];
  int     L    = cx->L;
  int     top  = ESL_MIN(L, b + cx->blockR - 1);
  int     Q    = p7O_NQF(om_fs->M);
  float   xN   = 0.0f;
  int     i, d, q;

  for (i = top; i >= ESL_MAX(b, 1); i--)
    {
      if (i == L)
        backward_frameshift_row_L_sse(om_fs, NULL, cx->xmx, totscale, rg, cx->dpf[L], L, scale_mode);
      else
        {
          for (d = 0; d < 6; d++) dp[d] = (i+d <= L) ? cx->dpf[i+d] : zero_row;
          backward_frameshift_row_sse(dsq, L, om_fs, NULL, cx->xmx, totscale, NULL, rg, dp, ivxf, i, scale_mode);
        }
    }

  if (b == 0)
    {
      dp[0] = NULL;
      for (d = 1; d < 6; d++) dp[d] = (d <= L) ? cx->dpf[d] : zero_row;
      xN  = backward_frameshift_row_0_sse(dsq, L, om_fs, cx->xmx, rg, dp, ivxf);
      xN += rg->xN[1%PARSER_ROWS_BWD] + rg->xN[2%PARSER_ROWS_BWD];
      for (q = 0; q < Q * p7X_NSCELLS; q++) cx->dpf[0][q] = _mm_setzero_ps();
    }
  return xN;
}


/* Function:  p7_BackwardChk_Frameshift_sse()
 * Synopsis:  Frameshift Backward in a checkpointed matrix, with its own scale factors.
 *
 * Purpose:   The Backward algorithm of <p7_Backward_Frameshift_sse()> for
 *            the 5-codon profile <om_fs> and sequence <dsq> of length
 *            <L>, computed in the checkpointed matrix <cx>, with <ov>
 *            as workspace.
 *
 *            This version is meant to be run before Forward, and it uses
 *            its own scale factors: each row is scaled by xB when xB >
 *            1e4 (else not scaled), not by the Forward matrix's factor for
 *            the row, so no Forward matrix is needed. The factors are
 *            stored in <cx->xmx> with the special states. Values and
 *            score equal those of <p7_Backward_Frameshift_sse()> up to
 *            rounding.
 *
 *            Rows are computed from <L> down to 0, one block at a time;
 *            checkpoints of the 5 rows after each block are saved, so
 *            <p7_BackwardChkLoad_Frameshift_sse()> can compute any block
 *            again, with identical values. On return block 0 is loaded
 *            and the special states and scale factors of all rows are in
 *            <cx->xmx>. Row 0 holds only special states; its core cells
 *            are set to 0.
 *
 *            The caller creates <cx> with <p7_CHK_BACKWARD>, <p7X_NSCELLS>
 *            vectors per stripe and at least 5 saved rows, and grows it
 *            for <om_fs->M> by <L> with the block size it wants
 *            (<p7_chkomx_GrowTo()>).
 *
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq (>= 1)
 *            om_fs  - 5-codon frameshift profile
 *            cx     - checkpointed DP matrix
 *            ov     - workspace for the profile's M
 *            opt_sc - optRETURN: Backward lod score in nats
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if the score underflows to 0; <*opt_sc> is then
 *            -infinity.
 *
 * Throws:    <eslEINVAL> if <cx> isn't set up for this computation or
 *            <om_fs> isn't a 5-codon profile.
 *            <eslERANGE> if the score is NaN or overflows.
 *            <eslEMEM> on allocation failure.
 */
int
p7_BackwardChk_Frameshift_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, float *opt_sc)
{
  BWD_FS_RING rg;
  __m128     *zero_mem = NULL;
  __m128     *zero_row;
  int         Q        = p7O_NQF(om_fs->M);
  float       totscale = 0.0f;
  float       xNtot    = 0.0f;
  int         b, q;
  int         status;

  if (om_fs->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");
  if (cx->dir != p7_CHK_BACKWARD || cx->saveR < 5 || cx->nscells != p7X_NSCELLS)
    ESL_EXCEPTION(eslEINVAL, "matrix not set up for the checkpointed frameshift Backward");
  if (L < 1) ESL_EXCEPTION(eslEINVAL, "sequence length must be >= 1");

  ESL_ALLOC(zero_mem, sizeof(__m128) * Q * p7X_NSCELLS + 15);
  zero_row = (__m128 *) (((unsigned long int) zero_mem + 15) & (~0xf));
  for (q = 0; q < Q * p7X_NSCELLS; q++) zero_row[q] = _mm_setzero_ps();

  cx->M = om_fs->M;
  cx->L = L;

  bwd_fs_ring_init(&rg, om_fs, L);
  for (b = (L / cx->blockR) * cx->blockR; b >= 0; b -= cx->blockR)
    {
      if ((status = p7_chkomx_StartBlock_sse(cx, b)) != eslOK) goto ERROR;
      xNtot = backward_frameshift_chk_block_sse(dsq, om_fs, cx, &rg, zero_row, ov->ivx[0], b, &totscale, BWD_SCALE_OWN);
    }
  free(zero_mem);

  if      (isnan(xNtot))      ESL_EXCEPTION(eslERANGE, "backward score is NaN");
  else if (xNtot == 0.0f)     { if (opt_sc != NULL) *opt_sc = -eslINFINITY; return eslERANGE; }
  else if (isinf(xNtot) == 1) ESL_EXCEPTION(eslERANGE, "backward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = totscale + logf(xNtot);
  return eslOK;

 ERROR:
  free(zero_mem);
  return status;
}


/* Function:  p7_BackwardChkLoad_Frameshift_sse()
 * Synopsis:  Compute one block of a checkpointed frameshift Backward matrix again.
 *
 * Purpose:   Load block <b> (a multiple of <cx->blockR>) of the Backward
 *            matrix that <p7_BackwardChk_Frameshift_sse()> filled in <cx>
 *            for profile <om_fs> and sequence <dsq>, computing its rows
 *            again from the checkpoint after it with the scale factors
 *            stored in <cx->xmx>, with <ov> as workspace. Afterwards rows
 *            <b> .. <b>+<blockR>-1 (up to <cx->L>) hold the same values as
 *            when <p7_BackwardChk_Frameshift_sse()> computed them. Blocks
 *            can be loaded in any order.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <b> isn't the start of a block.
 *            <eslEMEM> on allocation failure.
 */
int
p7_BackwardChkLoad_Frameshift_sse(const ESL_DSQ *dsq, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, int b)
{
  BWD_FS_RING rg;
  __m128     *zero_mem = NULL;
  __m128     *zero_row;
  int         Q        = p7O_NQF(om_fs->M);
  int         top;
  int         q;
  int         status;

  if ((status = p7_chkomx_LoadBlock_sse(cx, b)) != eslOK) return status;

  ESL_ALLOC(zero_mem, sizeof(__m128) * Q * p7X_NSCELLS + 15);
  zero_row = (__m128 *) (((unsigned long int) zero_mem + 15) & (~0xf));
  for (q = 0; q < Q * p7X_NSCELLS; q++) zero_row[q] = _mm_setzero_ps();

  top = ESL_MIN(cx->L, b + cx->blockR - 1);
  if (top == cx->L) bwd_fs_ring_init(&rg, om_fs, cx->L);
  else              bwd_fs_ring_load(&rg, om_fs, cx->xmx, cx->L, top);

  backward_frameshift_chk_block_sse(dsq, om_fs, cx, &rg, zero_row, ov->ivx[0], b, NULL, BWD_SCALE_STORED);
  free(zero_mem);
  return eslOK;

 ERROR:
  free(zero_mem);
  return status;
}

#endif /* eslENABLE_SSE */
