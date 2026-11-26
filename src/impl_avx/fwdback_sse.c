/* Forward and Backward algorithms, SSE implementation.
 * Ported from impl_sse/fwdback.c with _sse suffix for runtime dispatch.
 *
 * Contents:
 *    1. p7_Forward_sse, p7_ForwardParser_sse, p7_Backward_sse, p7_BackwardParser_sse
 *    2. Static engine implementations (forward_engine_sse, backward_engine_sse)
 *    3. p7_BackwardChk_sse, p7_BackwardChkLoad_sse: checkpointed Backward
 *       p7_ForwardChk_sse, p7_ForwardChkLoad_sse: checkpointed Forward
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

/* Function:  p7_ForwardRow_sse()
 * Synopsis:  One row of the SSE Forward algorithm.
 *
 * Purpose:   Forward row <i> of profile <om> against <dsq> into <dpc>,
 *            from row i-1 in <dpp> (not used for <i> = 0), reading the
 *            special states of row i-1 from <xmx> and writing row <i>'s
 *            special states and scale factor into it. When the row is
 *            rescaled, the log of its scale factor is added to
 *            <*totscale>. Row 0 is the initialization. Same calculation
 *            as each row of <p7_Forward_sse()>, which calls it.
 */
void
p7_ForwardRow_sse(const ESL_DSQ *dsq, const P7_OPROFILE *om, float *xmx, float *totscale, __m128 *dpc, const __m128 *dpp, int i)
{
  register __m128 mpv, dpv, ipv;
  register __m128 sv;
  register __m128 dcv;
  register __m128 xEv;
  register __m128 xBv;
  __m128   zerov = _mm_setzero_ps();
  float    xN, xE, xB, xC, xJ;
  int      q;
  int      j;
  int      Q     = p7O_NQF(om->M);
  __m128  *rp;
  __m128  *tp;

  if (i == 0)
    {
      for (q = 0; q < Q; q++)
        MMO(dpc,q) = IMO(dpc,q) = DMO(dpc,q) = zerov;
      xmx[p7X_E]     = 0.0f;
      xmx[p7X_N]     = 1.0f;
      xmx[p7X_J]     = 0.0f;
      xmx[p7X_B]     = om->xf[p7O_N][p7O_MOVE];
      xmx[p7X_C]     = 0.0f;
      xmx[p7X_SCALE] = 1.0f;
      return;
    }

  xN = xmx[(i-1)*p7X_NXCELLS+p7X_N];
  xJ = xmx[(i-1)*p7X_NXCELLS+p7X_J];
  xB = xmx[(i-1)*p7X_NXCELLS+p7X_B];
  xC = xmx[(i-1)*p7X_NXCELLS+p7X_C];

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
    xmx[i*p7X_NXCELLS+p7X_SCALE] = xE;
    *totscale += log(xE);
    xE = 1.0f;
  } else {
    xmx[i*p7X_NXCELLS+p7X_SCALE] = 1.0f;
  }

  xmx[i*p7X_NXCELLS+p7X_E] = xE;
  xmx[i*p7X_NXCELLS+p7X_N] = xN;
  xmx[i*p7X_NXCELLS+p7X_J] = xJ;
  xmx[i*p7X_NXCELLS+p7X_B] = xB;
  xmx[i*p7X_NXCELLS+p7X_C] = xC;
}


static int
forward_engine_sse(int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *opt_sc)
{
  float    xC;
  int      i;
  __m128  *dpc = ox->dpf[0];
  __m128  *dpp;

  ox->M  = om->M;
  ox->L  = L;
  ox->has_own_scales = TRUE;
  ox->totscale       = 0.0;
  p7_ForwardRow_sse(dsq, om, ox->xmx, &ox->totscale, dpc, NULL, 0);

#if eslDEBUGLEVEL > 0
  if (ox->debugging) p7_omx_DumpFBRow(ox, TRUE, 0, 9, 5, ox->xmx[p7X_E], ox->xmx[p7X_N], ox->xmx[p7X_J], ox->xmx[p7X_B], ox->xmx[p7X_C]);
#endif

  for (i = 1; i <= L; i++) {
    dpp = dpc;
    dpc = ox->dpf[do_full * i];
    p7_ForwardRow_sse(dsq, om, ox->xmx, &ox->totscale, dpc, dpp, i);
#if eslDEBUGLEVEL > 0
    if (ox->debugging) p7_omx_DumpFBRow(ox, TRUE, i, 9, 5, ox->xmx[i*p7X_NXCELLS+p7X_E], ox->xmx[i*p7X_NXCELLS+p7X_N], ox->xmx[i*p7X_NXCELLS+p7X_J], ox->xmx[i*p7X_NXCELLS+p7X_B], ox->xmx[i*p7X_NXCELLS+p7X_C]);
#endif
  }
  xC = ox->xmx[L*p7X_NXCELLS+p7X_C];

  if      (isnan(xC))        ESL_EXCEPTION(eslERANGE, "forward score is NaN");
  else if (L > 0 && xC == 0.0f) ESL_EXCEPTION(eslERANGE, "forward score underflow (is 0.0)");
  else if (isinf(xC) == 1)  ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = ox->totscale + log(xC * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}


/* Scale factor a Backward row takes, for backward_row_L_sse() and
 * backward_row_sse():
 *   BWD_SCALE_FORWARD: the Forward matrix's factor for the row, or the
 *     row's own factor once its xB has exceeded 1e16 (p7_Backward());
 *   BWD_SCALE_OWN:     the row's own factor, xB when xB > 1e4, else 1;
 *   BWD_SCALE_STORED:  the factor already stored for the row, when a
 *     row is computed again.
 */
enum bwd_scale_e { BWD_SCALE_FORWARD = 0, BWD_SCALE_OWN = 1, BWD_SCALE_STORED = 2 };


/* backward_row_L_sse()
 * Backward row <L> into <dpc>, special states and scale factor into
 * row <L> of <xmx>. <fwd_xmx> holds the Forward special states and
 * scale factors (needed only for BWD_SCALE_FORWARD). If <totscale> is
 * non-NULL it is set to the log of the row's scale factor.
 */
static void
backward_row_L_sse(const P7_OPROFILE *om, const float *fwd_xmx, float *xmx, float *totscale, __m128 *dpc, int L, int scale_mode)
{
  register __m128 dpv;
  register __m128 dcv;
  register __m128 xEv;
  __m128 zerov;
  float xN, xE, xB, xC, xJ;
  int q;
  int Q   = p7O_NQF(om->M);
  int j;
  __m128 *tp;
  float    sc;

  zerov = _mm_setzero_ps();
  xJ     = 0.0f;
  xB     = 0.0f;
  xN     = 0.0f;
  xC     = om->xf[p7O_C][p7O_MOVE];
  xE     = xC * om->xf[p7O_E][p7O_MOVE];
  xEv    = _mm_set1_ps(xE);
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
  if      (scale_mode == BWD_SCALE_FORWARD) sc = fwd_xmx[L*p7X_NXCELLS+p7X_SCALE];
  else if (scale_mode == BWD_SCALE_OWN)     sc = (xB > 1.0e4f) ? xB : 1.0f;
  else                                      sc = xmx[L*p7X_NXCELLS+p7X_SCALE];
  if (sc > 1.0f) {
    xE /= sc;  xN /= sc;  xC /= sc;  xJ /= sc;  xB /= sc;
    xEv = _mm_set1_ps(1.0f / sc);
    for (q = 0; q < Q; q++) {
      MMO(dpc,q) = _mm_mul_ps(MMO(dpc,q), xEv);
      DMO(dpc,q) = _mm_mul_ps(DMO(dpc,q), xEv);
      IMO(dpc,q) = _mm_mul_ps(IMO(dpc,q), xEv);
    }
  }
  xmx[L*p7X_NXCELLS+p7X_SCALE] = sc;
  if (totscale != NULL) *totscale = log(xmx[L*p7X_NXCELLS+p7X_SCALE]);

  xmx[L*p7X_NXCELLS+p7X_E] = xE;
  xmx[L*p7X_NXCELLS+p7X_N] = xN;
  xmx[L*p7X_NXCELLS+p7X_J] = xJ;
  xmx[L*p7X_NXCELLS+p7X_B] = xB;
  xmx[L*p7X_NXCELLS+p7X_C] = xC;
}


/* backward_row_sse()
 * Backward row <i> (1 <= i < L) into <dpc> from row <i>+1 in <dpp>,
 * reading the special states of row <i>+1 from <xmx> and writing row
 * <i>'s special states and scale factor into it. <fwd_xmx> holds the
 * Forward special states and scale factors (needed only for
 * BWD_SCALE_FORWARD). If <totscale> is non-NULL the log of a scale
 * factor > 1 is added to it. <has_own_scales> is the flag that
 * BWD_SCALE_FORWARD turns on once xB exceeds 1e16.
 */
static void
backward_row_sse(const ESL_DSQ *dsq, const P7_OPROFILE *om, const float *fwd_xmx, float *xmx, float *totscale, int *has_own_scales,
             __m128 *dpc, const __m128 *dpp, int i, int scale_mode)
{
  register __m128 mpv, ipv, dpv;
  register __m128 mcv, dcv;
  register __m128 tmmv, timv, tdmv;
  register __m128 xBv;
  register __m128 xEv;
  __m128 zerov;
  float xN, xE, xB, xC, xJ;
  int q;
  int Q   = p7O_NQF(om->M);
  int j;
  __m128 *rp;
  __m128 *tp;

  zerov = _mm_setzero_ps();
  xN = xmx[(i+1)*p7X_NXCELLS+p7X_N];
  xJ = xmx[(i+1)*p7X_NXCELLS+p7X_J];
  xC = xmx[(i+1)*p7X_NXCELLS+p7X_C];

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
  if (scale_mode == BWD_SCALE_FORWARD)
    {
      if (xB > 1.0e16f) *has_own_scales = TRUE;
      if (*has_own_scales)
        xmx[i*p7X_NXCELLS+p7X_SCALE] = (xB > 1.0e4f) ? xB : 1.0f;
      else
        xmx[i*p7X_NXCELLS+p7X_SCALE] = fwd_xmx[i*p7X_NXCELLS+p7X_SCALE];
    }
  else if (scale_mode == BWD_SCALE_OWN)
    xmx[i*p7X_NXCELLS+p7X_SCALE] = (xB > 1.0e4f) ? xB : 1.0f;

  if (xmx[i*p7X_NXCELLS+p7X_SCALE] > 1.0f) {
    float sc = xmx[i*p7X_NXCELLS+p7X_SCALE];
    xE /= sc;  xN /= sc;  xJ /= sc;  xB /= sc;  xC /= sc;
    xBv = _mm_set1_ps(1.0f / sc);
    for (q = 0; q < Q; q++) {
      MMO(dpc,q) = _mm_mul_ps(MMO(dpc,q), xBv);
      DMO(dpc,q) = _mm_mul_ps(DMO(dpc,q), xBv);
      IMO(dpc,q) = _mm_mul_ps(IMO(dpc,q), xBv);
    }
    if (totscale != NULL) *totscale += log(sc);
  }

  xmx[i*p7X_NXCELLS+p7X_E] = xE;
  xmx[i*p7X_NXCELLS+p7X_N] = xN;
  xmx[i*p7X_NXCELLS+p7X_J] = xJ;
  xmx[i*p7X_NXCELLS+p7X_B] = xB;
  xmx[i*p7X_NXCELLS+p7X_C] = xC;
}


/* backward_row_0_sse()
 * Backward row 0: only the N and B states, from row 1 in <dpp> and the
 * special states of row 1 in <xmx>. Writes row 0's special states into
 * <xmx> and returns its N value.
 */
static float
backward_row_0_sse(const ESL_DSQ *dsq, const P7_OPROFILE *om, float *xmx, const __m128 *dpp)
{
  register __m128 mpv;
  register __m128 xBv;
  __m128 zerov;
  float xN, xB;
  int q;
  int Q   = p7O_NQF(om->M);
  __m128 *rp;
  __m128 *tp;

  zerov = _mm_setzero_ps();
  xN = xmx[1*p7X_NXCELLS+p7X_N];
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

  xmx[p7X_B]     = xB;
  xmx[p7X_C]     = 0.0f;
  xmx[p7X_J]     = 0.0f;
  xmx[p7X_N]     = xN;
  xmx[p7X_E]     = 0.0f;
  xmx[p7X_SCALE] = 1.0f;
  return xN;
}


static int
backward_engine_sse(int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc)
{
  float xN;
  int i;
  __m128 *dpc;
  __m128 *dpp;
#if eslDEBUGLEVEL > 0
  __m128 zerov = _mm_setzero_ps();
  int q;
  int Q   = p7O_NQF(om->M);
#endif

  bck->M = om->M;
  bck->L = L;
  bck->has_own_scales = FALSE;
  dpc    = bck->dpf[L * do_full];
  backward_row_L_sse(om, fwd->xmx, bck->xmx, &bck->totscale, dpc, L, BWD_SCALE_FORWARD);

#if eslDEBUGLEVEL > 0
  if (bck->debugging) p7_omx_DumpFBRow(bck, TRUE, L, 9, 4, bck->xmx[L*p7X_NXCELLS+p7X_E], bck->xmx[L*p7X_NXCELLS+p7X_N], bck->xmx[L*p7X_NXCELLS+p7X_J], bck->xmx[L*p7X_NXCELLS+p7X_B], bck->xmx[L*p7X_NXCELLS+p7X_C]);
#endif

  /* Main backward recursion */
  for (i = L-1; i >= 1; i--) {
    dpc = bck->dpf[i     * do_full];
    dpp = bck->dpf[(i+1) * do_full];
    backward_row_sse(dsq, om, fwd->xmx, bck->xmx, &bck->totscale, &bck->has_own_scales, dpc, dpp, i, BWD_SCALE_FORWARD);

#if eslDEBUGLEVEL > 0
    if (bck->debugging) p7_omx_DumpFBRow(bck, TRUE, i, 9, 4, bck->xmx[i*p7X_NXCELLS+p7X_E], bck->xmx[i*p7X_NXCELLS+p7X_N], bck->xmx[i*p7X_NXCELLS+p7X_J], bck->xmx[i*p7X_NXCELLS+p7X_B], bck->xmx[i*p7X_NXCELLS+p7X_C]);
#endif
  }

  /* Termination at i=0 */
  dpp = bck->dpf[1 * do_full];
  xN  = backward_row_0_sse(dsq, om, bck->xmx, dpp);

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

/*****************************************************************
 * 3. Checkpointed Backward, SSE.
 *****************************************************************/

/* backward_chk_block_sse()
 * Compute rows <b> .. <b>+<blockR>-1 (up to <cx->L>) of the loaded
 * block <b> of the checkpointed Backward matrix <cx>, top row first,
 * with scale factors chosen by <scale_mode>. For block 0 this includes
 * row 0, whose core cells are set to 0. If <totscale> is non-NULL the
 * logs of the scale factors are added to it. Returns row 0's N value
 * for block 0, else 0.
 */
static float
backward_chk_block_sse(const ESL_DSQ *dsq, const P7_OPROFILE *om, P7_CHKOMX *cx, int b, float *totscale, int scale_mode)
{
  int   L   = cx->L;
  int   top = ESL_MIN(L, b + cx->blockR - 1);
  int   Q   = p7O_NQF(om->M);
  float xN  = 0.0f;
  int   i, q;

  for (i = top; i >= ESL_MAX(b, 1); i--)
    {
      if (i == L) backward_row_L_sse(om, NULL, cx->xmx, totscale, cx->dpf[L], L, scale_mode);
      else        backward_row_sse(dsq, om, NULL, cx->xmx, totscale, NULL, cx->dpf[i], cx->dpf[i+1], i, scale_mode);
    }

  if (b == 0)
    {
      xN = backward_row_0_sse(dsq, om, cx->xmx, cx->dpf[1]);
      for (q = 0; q < Q * p7X_NSCELLS; q++) cx->dpf[0][q] = _mm_setzero_ps();
    }
  return xN;
}


/* Function:  p7_BackwardChk_sse()
 * Synopsis:  Backward in a checkpointed matrix, with its own scale factors.
 *
 * Purpose:   The Backward algorithm of <p7_Backward_sse()> for profile
 *            <om> and sequence <dsq> of length <L>, computed in the
 *            checkpointed matrix <cx>.
 *
 *            This version is meant to be run before Forward, and it uses
 *            its own scale factors: each row is scaled by xB when xB >
 *            1e4 (else not scaled), not by the Forward matrix's factor for
 *            the row, so no Forward matrix is needed. The factors are
 *            stored in <cx->xmx> with the special states. Values and
 *            score equal those of <p7_Backward_sse()> up to rounding.
 *
 *            Rows are computed from <L> down to 0, one block at a time;
 *            checkpoints of the row after each block are saved, so
 *            <p7_BackwardChkLoad_sse()> can compute any block again, with
 *            identical values. On return block 0 is loaded and the special
 *            states and scale factors of all rows are in <cx->xmx>. Row 0
 *            holds only special states; its core cells are set to 0.
 *
 *            The caller creates <cx> with <p7_CHK_BACKWARD>, <p7X_NSCELLS>
 *            vectors per stripe and at least 1 saved row, and grows it for
 *            <om->M> by <L> with the block size it wants
 *            (<p7_chkomx_GrowTo()>).
 *
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq (>= 1)
 *            om     - profile
 *            cx     - checkpointed DP matrix
 *            opt_sc - optRETURN: Backward lod score in nats
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <cx> isn't set up for this computation.
 *            <eslERANGE> if the score is NaN, underflows to 0 or overflows.
 */
int
p7_BackwardChk_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc)
{
  float totscale = 0.0f;
  float xN       = 0.0f;
  int   b;
  int   status;

  if (cx->dir != p7_CHK_BACKWARD || cx->saveR < 1 || cx->nscells != p7X_NSCELLS)
    ESL_EXCEPTION(eslEINVAL, "matrix not set up for the checkpointed Backward");
  if (L < 1) ESL_EXCEPTION(eslEINVAL, "sequence length must be >= 1");

  cx->M = om->M;
  cx->L = L;

  for (b = (L / cx->blockR) * cx->blockR; b >= 0; b -= cx->blockR)
    {
      if ((status = p7_chkomx_StartBlock_sse(cx, b)) != eslOK) return status;
      xN = backward_chk_block_sse(dsq, om, cx, b, &totscale, BWD_SCALE_OWN);
    }

  if      (isnan(xN))      ESL_EXCEPTION(eslERANGE, "backward score is NaN");
  else if (xN == 0.0f)     ESL_EXCEPTION(eslERANGE, "backward score underflow (is 0.0)");
  else if (isinf(xN) == 1) ESL_EXCEPTION(eslERANGE, "backward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = totscale + log(xN);
  return eslOK;
}


/* Function:  p7_BackwardChkLoad_sse()
 * Synopsis:  Compute one block of a checkpointed Backward matrix again.
 *
 * Purpose:   Load block <b> (a multiple of <cx->blockR>) of the Backward
 *            matrix that <p7_BackwardChk_sse()> filled in <cx> for profile
 *            <om> and sequence <dsq>, computing its rows again from the
 *            checkpoint after it with the scale factors stored in
 *            <cx->xmx>. Afterwards rows <b> .. <b>+<blockR>-1 (up to
 *            <cx->L>) hold the same values as when <p7_BackwardChk_sse()>
 *            computed them. Blocks can be loaded in any order.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <b> isn't the start of a block.
 */
int
p7_BackwardChkLoad_sse(const ESL_DSQ *dsq, const P7_OPROFILE *om, P7_CHKOMX *cx, int b)
{
  int status;

  if ((status = p7_chkomx_LoadBlock_sse(cx, b)) != eslOK) return status;
  backward_chk_block_sse(dsq, om, cx, b, NULL, BWD_SCALE_STORED);
  return eslOK;
}

/* Function:  p7_ForwardChk_sse()
 * Synopsis:  Checkpointed Forward, SSE.
 *
 * Purpose:   Forward of profile <om> against <dsq> of length <L> into
 *            the checkpointed matrix <cx>. Each row is computed as in
 *            <p7_Forward_sse()>, by <p7_ForwardRow_sse()>, so cells,
 *            special states, scale factors and score are identical to
 *            it.
 *
 *            Rows are computed from 0 up to <L>, one block at a time;
 *            checkpoints of the rows before each block are saved, so
 *            <p7_ForwardChkLoad_sse()> can compute any block again, with
 *            identical values. On return the block holding row <L> is
 *            loaded and the special states and scale factors of all rows
 *            are in <cx->xmx>.
 *
 *            The caller creates <cx> with <p7_CHK_FORWARD> (or sets that
 *            direction with <p7_chkomx_SetDirection()>), <p7X_NSCELLS>
 *            vectors per stripe and at least 1 saved row, and grows it
 *            for <om->M> by <L> with the block size it wants
 *            (<p7_chkomx_GrowTo()>).
 *
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq (>= 1)
 *            om     - profile
 *            cx     - checkpointed DP matrix
 *            opt_sc - optRETURN: Forward lod score in nats
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <cx> isn't set up for this computation.
 *            <eslERANGE> if the score is NaN, underflows to 0 or overflows.
 */
int
p7_ForwardChk_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc)
{
  float totscale = 0.0f;
  float xC;
  int   b, i;
  int   status;

  if (cx->dir != p7_CHK_FORWARD || cx->saveR < 1 || cx->nscells != p7X_NSCELLS)
    ESL_EXCEPTION(eslEINVAL, "matrix not set up for the checkpointed Forward");
  if (L < 1) ESL_EXCEPTION(eslEINVAL, "sequence length must be >= 1");

  cx->M = om->M;
  cx->L = L;

  for (b = 0; b <= L; b += cx->blockR)
    {
      if ((status = p7_chkomx_StartBlock_sse(cx, b)) != eslOK) return status;
      for (i = b; i <= ESL_MIN(L, b + cx->blockR - 1); i++)
        p7_ForwardRow_sse(dsq, om, cx->xmx, &totscale, cx->dpf[i], (i > 0 ? cx->dpf[i-1] : NULL), i);
    }
  xC = cx->xmx[L*p7X_NXCELLS+p7X_C];

  if      (isnan(xC))      ESL_EXCEPTION(eslERANGE, "forward score is NaN");
  else if (xC == 0.0f)     ESL_EXCEPTION(eslERANGE, "forward score underflow (is 0.0)");
  else if (isinf(xC) == 1) ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = totscale + log(xC * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}


/* Function:  p7_ForwardChkLoad_sse()
 * Synopsis:  Compute one block of a checkpointed Forward again.
 *
 * Purpose:   After <p7_ForwardChk_sse()> has filled <cx>, compute block
 *            <b> (a multiple of <cx->blockR>) again from the checkpoint
 *            before it, loading the block and its saved rows into <cx>.
 *            Values are identical to those of the fill. Blocks are loaded
 *            in decreasing order, as a computation that runs backward
 *            through the rows reads them (<p7_chkomx_LoadBlock()>).
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if block <b> can't be loaded.
 */
int
p7_ForwardChkLoad_sse(const ESL_DSQ *dsq, const P7_OPROFILE *om, P7_CHKOMX *cx, int b)
{
  float totscale = 0.0f;   /* unused: the scale factors are already in <cx->xmx> */
  int   i;
  int   status;

  if ((status = p7_chkomx_LoadBlock_sse(cx, b)) != eslOK) return status;
  for (i = b; i <= ESL_MIN(cx->L, b + cx->blockR - 1); i++)
    p7_ForwardRow_sse(dsq, om, cx->xmx, &totscale, cx->dpf[i], (i > 0 ? cx->dpf[i-1] : NULL), i);
  return eslOK;
}


#endif /* eslENABLE_SSE */
