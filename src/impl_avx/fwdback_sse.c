/* Forward and Backward algorithms, SSE implementation.
 * Ported from impl_sse/fwdback.c with _sse suffix for runtime dispatch.
 *
 * Contents:
 *    1. p7_Forward_sse, p7_ForwardParser_sse, p7_Backward_sse, p7_BackwardParser_sse
 *    2. Static engine implementations (forward_engine_sse, backward_engine_sse)
 */

#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#ifdef eslENABLE_SSE
#include <xmmintrin.h>
#include <emmintrin.h>
#endif

#include "easel.h"
#include "esl_sse.h"

#include "hmmer.h"
#include "impl_avx.h"

#ifdef eslENABLE_SSE
static int forward_engine_sse (int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
static int backward_engine_sse(int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);
#endif

/*****************************************************************
 * 1. Forward/Backward SSE API.
 *****************************************************************/

int
p7_Forward_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *opt_sc)
{
#ifdef eslENABLE_SSE
#if eslDEBUGLEVEL > 0
  if (om->M >  ox->allocQ4*4)    ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few columns)");
  if (L     >= ox->validR)       ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few MDI rows)");
  if (L     >= ox->allocXR)      ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few X rows)");
  if (! p7_oprofile_IsLocal(om)) ESL_EXCEPTION(eslEINVAL, "Forward implementation makes assumptions that only work for local alignment");
#endif
  return forward_engine_sse(TRUE, dsq, L, om, ox, opt_sc);
#else
  return eslENORESULT;
#endif
}

int
p7_ForwardParser_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *opt_sc)
{
#ifdef eslENABLE_SSE
#if eslDEBUGLEVEL > 0
  if (om->M >  ox->allocQ4*4)    ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few columns)");
  if (ox->validR < 1)            ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few MDI rows)");
  if (L     >= ox->allocXR)      ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few X rows)");
  if (! p7_oprofile_IsLocal(om)) ESL_EXCEPTION(eslEINVAL, "Forward implementation makes assumptions that only work for local alignment");
#endif
  return forward_engine_sse(FALSE, dsq, L, om, ox, opt_sc);
#else
  return eslENORESULT;
#endif
}

int
p7_Backward_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc)
{
#ifdef eslENABLE_SSE
#if eslDEBUGLEVEL > 0
  if (om->M >  bck->allocQ4*4)    ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few columns)");
  if (L     >= bck->validR)       ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few MDI rows)");
  if (L     >= bck->allocXR)      ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few X rows)");
  if (L     != fwd->L)            ESL_EXCEPTION(eslEINVAL, "fwd matrix size doesn't agree with length L");
  if (! p7_oprofile_IsLocal(om))  ESL_EXCEPTION(eslEINVAL, "Forward implementation makes assumptions that only work for local alignment");
#endif
  return backward_engine_sse(TRUE, dsq, L, om, fwd, bck, opt_sc);
#else
  return eslENORESULT;
#endif
}

int
p7_BackwardParser_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc)
{
#ifdef eslENABLE_SSE
#if eslDEBUGLEVEL > 0
  if (om->M >  bck->allocQ4*4)    ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few columns)");
  if (bck->validR < 1)            ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few MDI rows)");
  if (L     >= bck->allocXR)      ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few X rows)");
  if (L     != fwd->L)            ESL_EXCEPTION(eslEINVAL, "fwd matrix size doesn't agree with length L");
  if (! p7_oprofile_IsLocal(om))  ESL_EXCEPTION(eslEINVAL, "Forward implementation makes assumptions that only work for local alignment");
#endif
  return backward_engine_sse(FALSE, dsq, L, om, fwd, bck, opt_sc);
#else
  return eslENORESULT;
#endif
}


/*****************************************************************
 * 2. Static engine implementations.
 *****************************************************************/

#ifdef eslENABLE_SSE

static int
forward_engine_sse(int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *opt_sc)
{
  register __m128 mpv, dpv, ipv;
  register __m128 sv;
  register __m128 dcv;
  register __m128 xEv;
  register __m128 xBv;
  __m128   zerov;
  float    xN, xE, xB, xC, xJ;
  int i;
  int q;
  int j;
  int Q     = p7O_NQF(om->M);
  __m128 *dpc = ox->dpf[0];
  __m128 *dpp;
  __m128 *rp;
  __m128 *tp;

  ox->M  = om->M;
  ox->L  = L;
  ox->has_own_scales = TRUE;
  zerov  = _mm_setzero_ps();
  for (q = 0; q < Q; q++)
    MMO(dpc,q) = IMO(dpc,q) = DMO(dpc,q) = zerov;
  xE = ox->xmx[p7X_E] = 0.0f;
  xN = ox->xmx[p7X_N] = 1.0f;
  xJ = ox->xmx[p7X_J] = 0.0f;
  xB = ox->xmx[p7X_B] = om->xf[p7O_N][p7O_MOVE];
  xC = ox->xmx[p7X_C] = 0.0f;

  ox->xmx[p7X_SCALE] = 1.0f;
  ox->totscale        = 0.0;

#if eslDEBUGLEVEL > 0
  if (ox->debugging) p7_omx_DumpFBRow(ox, TRUE, 0, 9, 5, xE, xN, xJ, xB, xC);
#endif

  for (i = 1; i <= L; i++) {
    dpp = dpc;
    dpc = ox->dpf[do_full * i];
    rp  = om->rfv[dsq[i]];
    tp  = om->tfv;
    dcv = _mm_setzero_ps();
    xEv = _mm_setzero_ps();
    xBv = _mm_set1_ps(xB);

    mpv = esl_sse_rightshiftz_float(MMO(dpp,Q-1));
    dpv = esl_sse_rightshiftz_float(DMO(dpp,Q-1));
    ipv = esl_sse_rightshiftz_float(IMO(dpp,Q-1));

    for (q = 0; q < Q; q++) {
      sv   =                _mm_mul_ps(xBv, *tp);  tp++;
      sv   = _mm_add_ps(sv, _mm_mul_ps(mpv, *tp)); tp++;
      sv   = _mm_add_ps(sv, _mm_mul_ps(ipv, *tp)); tp++;
      sv   = _mm_add_ps(sv, _mm_mul_ps(dpv, *tp)); tp++;
      sv   = _mm_mul_ps(sv, *rp);                  rp++;
      xEv  = _mm_add_ps(xEv, sv);

      mpv = MMO(dpp,q);
      dpv = DMO(dpp,q);
      ipv = IMO(dpp,q);

      MMO(dpc,q) = sv;
      DMO(dpc,q) = dcv;

      dcv = _mm_mul_ps(sv, *tp); tp++;

      sv         =                _mm_mul_ps(mpv, *tp);  tp++;
      IMO(dpc,q) = _mm_add_ps(sv, _mm_mul_ps(ipv, *tp)); tp++;
    }

    /* DD paths */
    dcv        = esl_sse_rightshiftz_float(dcv);
    DMO(dpc,0) = zerov;
    tp         = om->tfv + 7*Q;
    for (q = 0; q < Q; q++) {
      DMO(dpc,q) = _mm_add_ps(dcv, DMO(dpc,q));
      dcv        = _mm_mul_ps(DMO(dpc,q), *tp); tp++;
    }

    if (om->M < 100) {
      for (j = 1; j < 4; j++) {
        dcv = esl_sse_rightshiftz_float(dcv);
        tp  = om->tfv + 7*Q;
        for (q = 0; q < Q; q++) {
          DMO(dpc,q) = _mm_add_ps(dcv, DMO(dpc,q));
          dcv        = _mm_mul_ps(dcv, *tp); tp++;
        }
      }
    } else {
      for (j = 1; j < 4; j++) {
        register __m128 cv;
        dcv = esl_sse_rightshiftz_float(dcv);
        tp  = om->tfv + 7*Q;
        cv  = zerov;
        for (q = 0; q < Q; q++) {
          sv         = _mm_add_ps(dcv, DMO(dpc,q));
          cv         = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO(dpc,q)));
          DMO(dpc,q) = sv;
          dcv        = _mm_mul_ps(dcv, *tp); tp++;
        }
        if (!_mm_movemask_ps(cv)) break;
      }
    }

    /* D's contribute to xE in Forward */
    for (q = 0; q < Q; q++) xEv = _mm_add_ps(DMO(dpc,q), xEv);

    /* Horizontal sum of xEv */
    xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(0, 3, 2, 1)));
    xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(1, 0, 3, 2)));
    _mm_store_ss(&xE, xEv);

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
      xEv = _mm_set1_ps(1.0f / xE);
      for (q = 0; q < Q; q++) {
        MMO(dpc,q) = _mm_mul_ps(MMO(dpc,q), xEv);
        DMO(dpc,q) = _mm_mul_ps(DMO(dpc,q), xEv);
        IMO(dpc,q) = _mm_mul_ps(IMO(dpc,q), xEv);
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

#if eslDEBUGLEVEL > 0
    if (ox->debugging) p7_omx_DumpFBRow(ox, TRUE, i, 9, 5, xE, xN, xJ, xB, xC);
#endif
  }

  if      (isnan(xC))        ESL_EXCEPTION(eslERANGE, "forward score is NaN");
  else if (L > 0 && xC == 0.0f) ESL_EXCEPTION(eslERANGE, "forward score underflow (is 0.0)");
  else if (isinf(xC) == 1)  ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = ox->totscale + log(xC * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}


static int
backward_engine_sse(int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc)
{
  register __m128 mpv, ipv, dpv;
  register __m128 mcv, dcv;
  register __m128 tmmv, timv, tdmv;
  register __m128 xBv;
  register __m128 xEv;
  __m128   zerov;
  float    xN, xE, xB, xC, xJ;
  int      i;
  int      q;
  int      Q   = p7O_NQF(om->M);
  int      j;
  __m128  *dpc;
  __m128  *dpp;
  __m128  *rp;
  __m128  *tp;

  bck->M = om->M;
  bck->L = L;
  bck->has_own_scales = FALSE;
  dpc    = bck->dpf[L * do_full];
  xJ     = 0.0f;
  xB     = 0.0f;
  xN     = 0.0f;
  xC     = om->xf[p7O_C][p7O_MOVE];
  xE     = xC * om->xf[p7O_E][p7O_MOVE];
  xEv    = _mm_set1_ps(xE);
  zerov  = _mm_setzero_ps();
  dcv    = zerov;
  for (q = 0; q < Q; q++) MMO(dpc,q) = DMO(dpc,q) = xEv;
  for (q = 0; q < Q; q++) IMO(dpc,q) = zerov;

  /* init L row DD paths */
  tp  = om->tfv + 8*Q - 1;
  dpv = _mm_move_ss(DMO(dpc,Q-1), zerov);
  dpv = _mm_shuffle_ps(dpv, dpv, _MM_SHUFFLE(0,3,2,1));
  for (q = Q-1; q >= 0; q--) {
    dcv        = _mm_mul_ps(dpv, *tp);      tp--;
    DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), dcv);
    dpv        = DMO(dpc,q);
  }
  for (j = 1; j < 4; j++) {
    tp  = om->tfv + 8*Q - 1;
    dcv = _mm_move_ss(dcv, zerov);
    dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(0,3,2,1));
    for (q = Q-1; q >= 0; q--) {
      dcv        = _mm_mul_ps(dcv, *tp); tp--;
      DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), dcv);
    }
  }
  /* MD init */
  tp  = om->tfv + 7*Q - 3;
  dcv = _mm_move_ss(DMO(dpc,0), zerov);
  dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(0,3,2,1));
  for (q = Q-1; q >= 0; q--) {
    MMO(dpc,q) = _mm_add_ps(MMO(dpc,q), _mm_mul_ps(dcv, *tp)); tp -= 7;
    dcv        = DMO(dpc,q);
  }

  /* Sparse rescaling for row L */
  if (fwd->xmx[L*p7X_NXCELLS+p7X_SCALE] > 1.0f) {
    float sc = fwd->xmx[L*p7X_NXCELLS+p7X_SCALE];
    xE /= sc;  xN /= sc;  xC /= sc;  xJ /= sc;  xB /= sc;
    xEv = _mm_set1_ps(1.0f / sc);
    for (q = 0; q < Q; q++) {
      MMO(dpc,q) = _mm_mul_ps(MMO(dpc,q), xEv);
      DMO(dpc,q) = _mm_mul_ps(DMO(dpc,q), xEv);
      IMO(dpc,q) = _mm_mul_ps(IMO(dpc,q), xEv);
    }
  }
  bck->xmx[L*p7X_NXCELLS+p7X_SCALE] = fwd->xmx[L*p7X_NXCELLS+p7X_SCALE];
  bck->totscale                      = log(bck->xmx[L*p7X_NXCELLS+p7X_SCALE]);

  bck->xmx[L*p7X_NXCELLS+p7X_E] = xE;
  bck->xmx[L*p7X_NXCELLS+p7X_N] = xN;
  bck->xmx[L*p7X_NXCELLS+p7X_J] = xJ;
  bck->xmx[L*p7X_NXCELLS+p7X_B] = xB;
  bck->xmx[L*p7X_NXCELLS+p7X_C] = xC;

#if eslDEBUGLEVEL > 0
  if (bck->debugging) p7_omx_DumpFBRow(bck, TRUE, L, 9, 4, xE, xN, xJ, xB, xC);
#endif

  /* Main backward recursion */
  for (i = L-1; i >= 1; i--) {
    dpc = bck->dpf[i     * do_full];
    dpp = bck->dpf[(i+1) * do_full];
    rp  = om->rfv[dsq[i+1]] + Q-1;
    tp  = om->tfv + 7*Q - 1;

    tmmv = _mm_move_ss(om->tfv[1], zerov); tmmv = _mm_shuffle_ps(tmmv, tmmv, _MM_SHUFFLE(0,3,2,1));
    timv = _mm_move_ss(om->tfv[2], zerov); timv = _mm_shuffle_ps(timv, timv, _MM_SHUFFLE(0,3,2,1));
    tdmv = _mm_move_ss(om->tfv[3], zerov); tdmv = _mm_shuffle_ps(tdmv, tdmv, _MM_SHUFFLE(0,3,2,1));

    mpv = _mm_mul_ps(MMO(dpp,0), om->rfv[dsq[i+1]][0]);
    mpv = _mm_move_ss(mpv, zerov);
    mpv = _mm_shuffle_ps(mpv, mpv, _MM_SHUFFLE(0,3,2,1));

    xBv = zerov;
    for (q = Q-1; q >= 0; q--) {
      ipv = IMO(dpp,q);
      IMO(dpc,q) = _mm_add_ps(_mm_mul_ps(ipv, *tp), _mm_mul_ps(mpv, timv));   tp--;
      DMO(dpc,q) =                                   _mm_mul_ps(mpv, tdmv);
      mcv        = _mm_add_ps(_mm_mul_ps(ipv, *tp), _mm_mul_ps(mpv, tmmv));   tp -= 2;

      mpv        = _mm_mul_ps(MMO(dpp,q), *rp);  rp--;
      MMO(dpc,q) = mcv;

      tdmv = *tp;  tp--;
      timv = *tp;  tp--;
      tmmv = *tp;  tp--;

      xBv = _mm_add_ps(xBv, _mm_mul_ps(mpv, *tp)); tp--;
    }

    /* Specials */
    xBv = _mm_add_ps(xBv, _mm_shuffle_ps(xBv, xBv, _MM_SHUFFLE(0, 3, 2, 1)));
    xBv = _mm_add_ps(xBv, _mm_shuffle_ps(xBv, xBv, _MM_SHUFFLE(1, 0, 3, 2)));
    _mm_store_ss(&xB, xBv);

    xC = xC * om->xf[p7O_C][p7O_LOOP];
    xJ = (xB * om->xf[p7O_J][p7O_MOVE]) + (xJ * om->xf[p7O_J][p7O_LOOP]);
    xN = (xB * om->xf[p7O_N][p7O_MOVE]) + (xN * om->xf[p7O_N][p7O_LOOP]);
    xE = (xC * om->xf[p7O_E][p7O_MOVE]) + (xJ * om->xf[p7O_E][p7O_LOOP]);
    xEv = _mm_set1_ps(xE);

    /* {MD}->E and D->D */
    tp  = om->tfv + 8*Q - 1;
    dpv = _mm_add_ps(DMO(dpc,0), xEv);
    dpv = _mm_move_ss(dpv, zerov);
    dpv = _mm_shuffle_ps(dpv, dpv, _MM_SHUFFLE(0,3,2,1));
    for (q = Q-1; q >= 0; q--) {
      dcv        = _mm_mul_ps(dpv, *tp); tp--;
      DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), _mm_add_ps(dcv, xEv));
      dpv        = DMO(dpc,q);
      MMO(dpc,q) = _mm_add_ps(MMO(dpc,q), xEv);
    }
    for (j = 1; j < 4; j++) {
      dcv = _mm_move_ss(dcv, zerov);
      dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(0,3,2,1));
      tp  = om->tfv + 8*Q - 1;
      for (q = Q-1; q >= 0; q--) {
        dcv        = _mm_mul_ps(dcv, *tp); tp--;
        DMO(dpc,q) = _mm_add_ps(DMO(dpc,q), dcv);
      }
    }
    /* M->D paths */
    dcv = _mm_move_ss(DMO(dpc,0), zerov);
    dcv = _mm_shuffle_ps(dcv, dcv, _MM_SHUFFLE(0,3,2,1));
    tp  = om->tfv + 7*Q - 3;
    for (q = Q-1; q >= 0; q--) {
      MMO(dpc,q) = _mm_add_ps(MMO(dpc,q), _mm_mul_ps(dcv, *tp)); tp -= 7;
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
      xBv = _mm_set1_ps(1.0f / sc);
      for (q = 0; q < Q; q++) {
        MMO(dpc,q) = _mm_mul_ps(MMO(dpc,q), xBv);
        DMO(dpc,q) = _mm_mul_ps(DMO(dpc,q), xBv);
        IMO(dpc,q) = _mm_mul_ps(IMO(dpc,q), xBv);
      }
      bck->totscale += log(sc);
    }

    bck->xmx[i*p7X_NXCELLS+p7X_E] = xE;
    bck->xmx[i*p7X_NXCELLS+p7X_N] = xN;
    bck->xmx[i*p7X_NXCELLS+p7X_J] = xJ;
    bck->xmx[i*p7X_NXCELLS+p7X_B] = xB;
    bck->xmx[i*p7X_NXCELLS+p7X_C] = xC;

#if eslDEBUGLEVEL > 0
    if (bck->debugging) p7_omx_DumpFBRow(bck, TRUE, i, 9, 4, xE, xN, xJ, xB, xC);
#endif
  }

  /* Termination at i=0 */
  dpp = bck->dpf[1 * do_full];
  tp  = om->tfv;
  rp  = om->rfv[dsq[1]];
  xBv = zerov;
  for (q = 0; q < Q; q++) {
    mpv = _mm_mul_ps(MMO(dpp,q), *rp);  rp++;
    mpv = _mm_mul_ps(mpv,        *tp);  tp += 7;
    xBv = _mm_add_ps(xBv,        mpv);
  }
  xBv = _mm_add_ps(xBv, _mm_shuffle_ps(xBv, xBv, _MM_SHUFFLE(0, 3, 2, 1)));
  xBv = _mm_add_ps(xBv, _mm_shuffle_ps(xBv, xBv, _MM_SHUFFLE(1, 0, 3, 2)));
  _mm_store_ss(&xB, xBv);

  xN = (xB * om->xf[p7O_N][p7O_MOVE]) + (xN * om->xf[p7O_N][p7O_LOOP]);

  bck->xmx[p7X_B]     = xB;
  bck->xmx[p7X_C]     = 0.0f;
  bck->xmx[p7X_J]     = 0.0f;
  bck->xmx[p7X_N]     = xN;
  bck->xmx[p7X_E]     = 0.0f;
  bck->xmx[p7X_SCALE] = 1.0f;

#if eslDEBUGLEVEL > 0
  dpc = bck->dpf[0];
  for (q = 0; q < Q; q++)
    MMO(dpc,q) = DMO(dpc,q) = IMO(dpc,q) = zerov;
  if (bck->debugging) p7_omx_DumpFBRow(bck, TRUE, 0, 9, 4,
                                        bck->xmx[p7X_E], bck->xmx[p7X_N],
                                        bck->xmx[p7X_J], bck->xmx[p7X_B],
                                        bck->xmx[p7X_C]);
#endif

  if      (isnan(xN))         ESL_EXCEPTION(eslERANGE, "backward score is NaN");
  else if (L > 0 && xN == 0.0f) ESL_EXCEPTION(eslERANGE, "backward score underflow (is 0.0)");
  else if (isinf(xN) == 1)   ESL_EXCEPTION(eslERANGE, "backward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = bck->totscale + log(xN);
  return eslOK;
}

#endif /* eslENABLE_SSE */
