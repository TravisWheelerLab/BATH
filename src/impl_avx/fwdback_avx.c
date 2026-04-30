/* Forward and Backward algorithms, AVX2 implementation.
 * Ported from fwdback_sse.c using 256-bit vectors.
 *
 * Contents:
 *    1. p7_Forward_avx, p7_ForwardParser_avx, p7_Backward_avx, p7_BackwardParser_avx
 *    2. Static engine implementations (forward_engine_avx, backward_engine_avx)
 */

#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <stdio.h>
#include <math.h>

#include <immintrin.h>

#include "easel.h"
#include "esl_avx.h"

#include "hmmer.h"
#include "impl_avx.h"

static int forward_engine_avx (int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
static int backward_engine_avx(int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);


/*****************************************************************
 * 1. Forward/Backward AVX2 API.
 *****************************************************************/

int
p7_Forward_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *opt_sc)
{
#if eslDEBUGLEVEL > 0
  if (om->M >  ox->allocQ4_avx*8)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few columns)");
  if (L     >= ox->validR)          ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few MDI rows)");
  if (L     >= ox->allocXR)         ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few X rows)");
  if (! p7_oprofile_IsLocal(om))    ESL_EXCEPTION(eslEINVAL, "Forward implementation makes assumptions that only work for local alignment");
#endif
  return forward_engine_avx(TRUE, dsq, L, om, ox, opt_sc);
}

int
p7_ForwardParser_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *opt_sc)
{
#if eslDEBUGLEVEL > 0
  if (om->M >  ox->allocQ4_avx*8)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few columns)");
  if (ox->validR < 1)               ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few MDI rows)");
  if (L     >= ox->allocXR)         ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few X rows)");
  if (! p7_oprofile_IsLocal(om))    ESL_EXCEPTION(eslEINVAL, "Forward implementation makes assumptions that only work for local alignment");
#endif
  return forward_engine_avx(FALSE, dsq, L, om, ox, opt_sc);
}

int
p7_Backward_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc)
{
#if eslDEBUGLEVEL > 0
  if (om->M >  bck->allocQ4_avx*8)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few columns)");
  if (L     >= bck->validR)          ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few MDI rows)");
  if (L     >= bck->allocXR)         ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few X rows)");
  if (L     != fwd->L)               ESL_EXCEPTION(eslEINVAL, "fwd matrix size doesn't agree with length L");
  if (! p7_oprofile_IsLocal(om))     ESL_EXCEPTION(eslEINVAL, "Forward implementation makes assumptions that only work for local alignment");
#endif
  return backward_engine_avx(TRUE, dsq, L, om, fwd, bck, opt_sc);
}

int
p7_BackwardParser_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc)
{
#if eslDEBUGLEVEL > 0
  if (om->M >  bck->allocQ4_avx*8)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few columns)");
  if (bck->validR < 1)               ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few MDI rows)");
  if (L     >= bck->allocXR)         ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few X rows)");
  if (L     != fwd->L)               ESL_EXCEPTION(eslEINVAL, "fwd matrix size doesn't agree with length L");
  if (! p7_oprofile_IsLocal(om))     ESL_EXCEPTION(eslEINVAL, "Forward implementation makes assumptions that only work for local alignment");
#endif
  return backward_engine_avx(FALSE, dsq, L, om, fwd, bck, opt_sc);
}


/*****************************************************************
 * 2. Static engine implementations.
 *****************************************************************/

static int
forward_engine_avx(int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *opt_sc)
{
  register __m256 mpv, dpv, ipv;
  register __m256 sv;
  register __m256 dcv;
  register __m256 xEv;
  register __m256 xBv;
  __m256   zerov;
  float    xN, xE, xB, xC, xJ;
  int i;
  int q;
  int j;
  int Q     = p7O_NQF_AVX(om->M);
  __m256 *dpc = ox->dpf_avx[0];
  __m256 *dpp;
  __m256 *rp;
  __m256 *tp;

  ox->M  = om->M;
  ox->L  = L;
  ox->has_own_scales = TRUE;
  zerov  = _mm256_setzero_ps();
  for (q = 0; q < Q; q++)
    MMO(dpc,q) = IMO(dpc,q) = DMO(dpc,q) = zerov;
  xE = ox->xmx[p7X_E] = 0.0f;
  xN = ox->xmx[p7X_N] = 1.0f;
  xJ = ox->xmx[p7X_J] = 0.0f;
  xB = ox->xmx[p7X_B] = om->xf[p7O_N][p7O_MOVE];
  xC = ox->xmx[p7X_C] = 0.0f;

  ox->xmx[p7X_SCALE] = 1.0f;
  ox->totscale        = 0.0;

  for (i = 1; i <= L; i++) {
    dpp = dpc;
    dpc = ox->dpf_avx[do_full * i];
    rp  = om->rfv_avx[dsq[i]];
    tp  = om->tfv_avx;
    dcv = zerov;
    xEv = zerov;
    xBv = _mm256_set1_ps(xB);

    mpv = esl_avx_rightshiftz_float(MMO(dpp,Q-1));
    dpv = esl_avx_rightshiftz_float(DMO(dpp,Q-1));
    ipv = esl_avx_rightshiftz_float(IMO(dpp,Q-1));

    for (q = 0; q < Q; q++) {
      sv   =                _mm256_mul_ps(xBv, *tp);  tp++;
      sv   = _mm256_add_ps(sv, _mm256_mul_ps(mpv, *tp)); tp++;
      sv   = _mm256_add_ps(sv, _mm256_mul_ps(ipv, *tp)); tp++;
      sv   = _mm256_add_ps(sv, _mm256_mul_ps(dpv, *tp)); tp++;
      sv   = _mm256_mul_ps(sv, *rp);                     rp++;
      xEv  = _mm256_add_ps(xEv, sv);

      mpv = MMO(dpp,q);
      dpv = DMO(dpp,q);
      ipv = IMO(dpp,q);

      MMO(dpc,q) = sv;
      DMO(dpc,q) = dcv;

      dcv = _mm256_mul_ps(sv, *tp); tp++;

      sv         =                _mm256_mul_ps(mpv, *tp);  tp++;
      IMO(dpc,q) = _mm256_add_ps(sv, _mm256_mul_ps(ipv, *tp)); tp++;
    }

    /* DD paths */
    dcv        = esl_avx_rightshiftz_float(dcv);
    DMO(dpc,0) = zerov;
    tp         = om->tfv_avx + 7*Q;
    for (q = 0; q < Q; q++) {
      DMO(dpc,q) = _mm256_add_ps(dcv, DMO(dpc,q));
      dcv        = _mm256_mul_ps(DMO(dpc,q), *tp); tp++;
    }

    if (om->M < 100) {
      for (j = 1; j < 4; j++) {
        dcv = esl_avx_rightshiftz_float(dcv);
        tp  = om->tfv_avx + 7*Q;
        for (q = 0; q < Q; q++) {
          DMO(dpc,q) = _mm256_add_ps(dcv, DMO(dpc,q));
          dcv        = _mm256_mul_ps(dcv, *tp); tp++;
        }
      }
    } else {
      for (j = 1; j < 4; j++) {
        register __m256 cv;
        dcv = esl_avx_rightshiftz_float(dcv);
        tp  = om->tfv_avx + 7*Q;
        cv  = zerov;
        for (q = 0; q < Q; q++) {
          sv         = _mm256_add_ps(dcv, DMO(dpc,q));
          cv         = _mm256_or_ps(cv, _mm256_cmp_ps(sv, DMO(dpc,q), _CMP_GT_OS));
          DMO(dpc,q) = sv;
          dcv        = _mm256_mul_ps(dcv, *tp); tp++;
        }
        if (!_mm256_movemask_ps(cv)) break;
      }
    }

    /* D's contribute to xE in Forward */
    for (q = 0; q < Q; q++) xEv = _mm256_add_ps(DMO(dpc,q), xEv);

    /* Horizontal sum of xEv */
    esl_avx_hsum_ps(xEv, &xE);

    xN = xN * om->xf[p7O_N][p7O_LOOP];
    xC = (xC * om->xf[p7O_C][p7O_LOOP]) + (xE * om->xf[p7O_E][p7O_MOVE]);
    xJ = (xJ * om->xf[p7O_J][p7O_LOOP]) + (xE * om->xf[p7O_E][p7O_LOOP]);
    xB = (xJ * om->xf[p7O_J][p7O_MOVE]) + (xN * om->xf[p7O_N][p7O_MOVE]);

    /* Sparse rescaling */
    if (xE > 1.0e4f) {
      xN  = xN / xE;
      xC  = xC / xE;
      xJ  = xJ / xE;
      xB  = xB / xE;
      xEv = _mm256_set1_ps(1.0f / xE);
      for (q = 0; q < Q; q++) {
        MMO(dpc,q) = _mm256_mul_ps(MMO(dpc,q), xEv);
        DMO(dpc,q) = _mm256_mul_ps(DMO(dpc,q), xEv);
        IMO(dpc,q) = _mm256_mul_ps(IMO(dpc,q), xEv);
      }
      ox->xmx[i*p7X_NXCELLS+p7X_SCALE] = xE;
      ox->totscale += log(xE);
      xE = 1.0f;
    } else {
      ox->xmx[i*p7X_NXCELLS+p7X_SCALE] = 1.0f;
    }

    ox->xmx[i*p7X_NXCELLS+p7X_E] = xE;
    ox->xmx[i*p7X_NXCELLS+p7X_N] = xN;
    ox->xmx[i*p7X_NXCELLS+p7X_J] = xJ;
    ox->xmx[i*p7X_NXCELLS+p7X_B] = xB;
    ox->xmx[i*p7X_NXCELLS+p7X_C] = xC;
  }

  if      (isnan(xC))             ESL_EXCEPTION(eslERANGE, "forward score is NaN");
  else if (L > 0 && xC == 0.0f)  ESL_EXCEPTION(eslERANGE, "forward score underflow (is 0.0)");
  else if (isinf(xC) == 1)        ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = ox->totscale + log(xC * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}


static int
backward_engine_avx(int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc)
{
  register __m256 mpv, ipv, dpv;
  register __m256 mcv, dcv;
  register __m256 tmmv, timv, tdmv;
  register __m256 xBv;
  register __m256 xEv;
  __m256   zerov;
  float    xN, xE, xB, xC, xJ;
  int      i;
  int      q;
  int      Q   = p7O_NQF_AVX(om->M);
  int      j;
  __m256  *dpc;
  __m256  *dpp;
  __m256  *rp;
  __m256  *tp;

  bck->M = om->M;
  bck->L = L;
  bck->has_own_scales = FALSE;
  dpc    = bck->dpf_avx[L * do_full];
  xJ     = 0.0f;
  xB     = 0.0f;
  xN     = 0.0f;
  xC     = om->xf[p7O_C][p7O_MOVE];
  xE     = xC * om->xf[p7O_E][p7O_MOVE];
  xEv    = _mm256_set1_ps(xE);
  zerov  = _mm256_setzero_ps();
  dcv    = zerov;
  for (q = 0; q < Q; q++) MMO(dpc,q) = DMO(dpc,q) = xEv;
  for (q = 0; q < Q; q++) IMO(dpc,q) = zerov;

  /* init L row DD paths */
  tp  = om->tfv_avx + 8*Q - 1;
  dpv = esl_avx_leftshiftz_float(DMO(dpc,Q-1));
  for (q = Q-1; q >= 0; q--) {
    dcv        = _mm256_mul_ps(dpv, *tp);      tp--;
    DMO(dpc,q) = _mm256_add_ps(DMO(dpc,q), dcv);
    dpv        = DMO(dpc,q);
  }
  for (j = 1; j < 4; j++) {
    tp  = om->tfv_avx + 8*Q - 1;
    dcv = esl_avx_leftshiftz_float(dcv);
    for (q = Q-1; q >= 0; q--) {
      dcv        = _mm256_mul_ps(dcv, *tp); tp--;
      DMO(dpc,q) = _mm256_add_ps(DMO(dpc,q), dcv);
    }
  }
  /* MD init */
  tp  = om->tfv_avx + 7*Q - 3;
  dcv = esl_avx_leftshiftz_float(DMO(dpc,0));
  for (q = Q-1; q >= 0; q--) {
    MMO(dpc,q) = _mm256_add_ps(MMO(dpc,q), _mm256_mul_ps(dcv, *tp)); tp -= 7;
    dcv        = DMO(dpc,q);
  }

  /* Sparse rescaling for row L */
  if (fwd->xmx[L*p7X_NXCELLS+p7X_SCALE] > 1.0f) {
    float sc = fwd->xmx[L*p7X_NXCELLS+p7X_SCALE];
    xE /= sc;  xN /= sc;  xC /= sc;  xJ /= sc;  xB /= sc;
    xEv = _mm256_set1_ps(1.0f / sc);
    for (q = 0; q < Q; q++) {
      MMO(dpc,q) = _mm256_mul_ps(MMO(dpc,q), xEv);
      DMO(dpc,q) = _mm256_mul_ps(DMO(dpc,q), xEv);
      IMO(dpc,q) = _mm256_mul_ps(IMO(dpc,q), xEv);
    }
  }
  bck->xmx[L*p7X_NXCELLS+p7X_SCALE] = fwd->xmx[L*p7X_NXCELLS+p7X_SCALE];
  bck->totscale                      = log(bck->xmx[L*p7X_NXCELLS+p7X_SCALE]);

  bck->xmx[L*p7X_NXCELLS+p7X_E] = xE;
  bck->xmx[L*p7X_NXCELLS+p7X_N] = xN;
  bck->xmx[L*p7X_NXCELLS+p7X_J] = xJ;
  bck->xmx[L*p7X_NXCELLS+p7X_B] = xB;
  bck->xmx[L*p7X_NXCELLS+p7X_C] = xC;

  /* Main backward recursion */
  for (i = L-1; i >= 1; i--) {
    dpc = bck->dpf_avx[i     * do_full];
    dpp = bck->dpf_avx[(i+1) * do_full];
    rp  = om->rfv_avx[dsq[i+1]] + Q-1;
    tp  = om->tfv_avx + 7*Q - 1;

    tmmv = esl_avx_leftshiftz_float(om->tfv_avx[1]);
    timv = esl_avx_leftshiftz_float(om->tfv_avx[2]);
    tdmv = esl_avx_leftshiftz_float(om->tfv_avx[3]);

    mpv = _mm256_mul_ps(MMO(dpp,0), om->rfv_avx[dsq[i+1]][0]);
    mpv = esl_avx_leftshiftz_float(mpv);

    xBv = zerov;
    for (q = Q-1; q >= 0; q--) {
      ipv = IMO(dpp,q);
      IMO(dpc,q) = _mm256_add_ps(_mm256_mul_ps(ipv, *tp), _mm256_mul_ps(mpv, timv));   tp--;
      DMO(dpc,q) =                                         _mm256_mul_ps(mpv, tdmv);
      mcv        = _mm256_add_ps(_mm256_mul_ps(ipv, *tp), _mm256_mul_ps(mpv, tmmv));   tp -= 2;

      mpv        = _mm256_mul_ps(MMO(dpp,q), *rp);  rp--;
      MMO(dpc,q) = mcv;

      tdmv = *tp;  tp--;
      timv = *tp;  tp--;
      tmmv = *tp;  tp--;

      xBv = _mm256_add_ps(xBv, _mm256_mul_ps(mpv, *tp)); tp--;
    }

    /* Specials */
    esl_avx_hsum_ps(xBv, &xB);

    xC = xC * om->xf[p7O_C][p7O_LOOP];
    xJ = (xB * om->xf[p7O_J][p7O_MOVE]) + (xJ * om->xf[p7O_J][p7O_LOOP]);
    xN = (xB * om->xf[p7O_N][p7O_MOVE]) + (xN * om->xf[p7O_N][p7O_LOOP]);
    xE = (xC * om->xf[p7O_E][p7O_MOVE]) + (xJ * om->xf[p7O_E][p7O_LOOP]);
    xEv = _mm256_set1_ps(xE);

    /* {MD}->E and D->D */
    tp  = om->tfv_avx + 8*Q - 1;
    dpv = _mm256_add_ps(DMO(dpc,0), xEv);
    dpv = esl_avx_leftshiftz_float(dpv);
    for (q = Q-1; q >= 0; q--) {
      dcv        = _mm256_mul_ps(dpv, *tp); tp--;
      DMO(dpc,q) = _mm256_add_ps(DMO(dpc,q), _mm256_add_ps(dcv, xEv));
      dpv        = DMO(dpc,q);
      MMO(dpc,q) = _mm256_add_ps(MMO(dpc,q), xEv);
    }
    for (j = 1; j < 4; j++) {
      dcv = esl_avx_leftshiftz_float(dcv);
      tp  = om->tfv_avx + 8*Q - 1;
      for (q = Q-1; q >= 0; q--) {
        dcv        = _mm256_mul_ps(dcv, *tp); tp--;
        DMO(dpc,q) = _mm256_add_ps(DMO(dpc,q), dcv);
      }
    }
    /* M->D paths */
    dcv = esl_avx_leftshiftz_float(DMO(dpc,0));
    tp  = om->tfv_avx + 7*Q - 3;
    for (q = Q-1; q >= 0; q--) {
      MMO(dpc,q) = _mm256_add_ps(MMO(dpc,q), _mm256_mul_ps(dcv, *tp)); tp -= 7;
      dcv        = DMO(dpc,q);
    }

    /* Sparse rescaling */
    if (xB > 1.0e16f) bck->has_own_scales = TRUE;
    if (bck->has_own_scales)
      bck->xmx[i*p7X_NXCELLS+p7X_SCALE] = (xB > 1.0e4f) ? xB : 1.0f;
    else
      bck->xmx[i*p7X_NXCELLS+p7X_SCALE] = fwd->xmx[i*p7X_NXCELLS+p7X_SCALE];

    if (bck->xmx[i*p7X_NXCELLS+p7X_SCALE] > 1.0f) {
      float sc = bck->xmx[i*p7X_NXCELLS+p7X_SCALE];
      xE /= sc;  xN /= sc;  xJ /= sc;  xB /= sc;  xC /= sc;
      xBv = _mm256_set1_ps(1.0f / sc);
      for (q = 0; q < Q; q++) {
        MMO(dpc,q) = _mm256_mul_ps(MMO(dpc,q), xBv);
        DMO(dpc,q) = _mm256_mul_ps(DMO(dpc,q), xBv);
        IMO(dpc,q) = _mm256_mul_ps(IMO(dpc,q), xBv);
      }
      bck->totscale += log(sc);
    }

    bck->xmx[i*p7X_NXCELLS+p7X_E] = xE;
    bck->xmx[i*p7X_NXCELLS+p7X_N] = xN;
    bck->xmx[i*p7X_NXCELLS+p7X_J] = xJ;
    bck->xmx[i*p7X_NXCELLS+p7X_B] = xB;
    bck->xmx[i*p7X_NXCELLS+p7X_C] = xC;
  }

  /* Termination at i=0 */
  dpp = bck->dpf_avx[1 * do_full];
  tp  = om->tfv_avx;
  rp  = om->rfv_avx[dsq[1]];
  xBv = zerov;
  for (q = 0; q < Q; q++) {
    mpv = _mm256_mul_ps(MMO(dpp,q), *rp);  rp++;
    mpv = _mm256_mul_ps(mpv,        *tp);  tp += 7;
    xBv = _mm256_add_ps(xBv,        mpv);
  }
  esl_avx_hsum_ps(xBv, &xB);

  xN = (xB * om->xf[p7O_N][p7O_MOVE]) + (xN * om->xf[p7O_N][p7O_LOOP]);

  bck->xmx[p7X_B]     = xB;
  bck->xmx[p7X_C]     = 0.0f;
  bck->xmx[p7X_J]     = 0.0f;
  bck->xmx[p7X_N]     = xN;
  bck->xmx[p7X_E]     = 0.0f;
  bck->xmx[p7X_SCALE] = 1.0f;

  if      (isnan(xN))              ESL_EXCEPTION(eslERANGE, "backward score is NaN");
  else if (L > 0 && xN == 0.0f)   ESL_EXCEPTION(eslERANGE, "backward score underflow (is 0.0)");
  else if (isinf(xN) == 1)         ESL_EXCEPTION(eslERANGE, "backward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = bck->totscale + log(xN);
  return eslOK;
}

#endif /* eslENABLE_AVX */
