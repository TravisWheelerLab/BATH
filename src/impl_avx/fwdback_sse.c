/* Forward and Backward algorithms, SSE implementation.
 * Ported from impl_sse/fwdback.c with _sse suffix for runtime dispatch.
 *
 * Contents:
 *    1. p7_Forward_sse, p7_ForwardParser_sse, p7_Backward_sse, p7_BackwardParser_sse
 *    2. Row functions p7_ForwardRow_sse, p7_BackwardRow_sse, and the static
 *       engines (forward_engine_sse, backward_engine_sse)
 *
 * The checkpointed Forward and Backward are in fwdback_chk_sse.c.
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

  for (i = 1; i <= L; i++) {
    dpp = dpc;
    dpc = ox->dpf[do_full * i];
    p7_ForwardRow_sse(dsq, om, ox->xmx, &ox->totscale, dpc, dpp, i);
  }
  xC = ox->xmx[L*p7X_NXCELLS+p7X_C];

  if      (isnan(xC))        ESL_EXCEPTION(eslERANGE, "forward score is NaN");
  else if (L > 0 && xC == 0.0f) ESL_EXCEPTION(eslERANGE, "forward score underflow (is 0.0)");
  else if (isinf(xC) == 1)  ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = ox->totscale + log(xC * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}


/* Function:  p7_BackwardRow_sse()
 * Synopsis:  One row of the Backward algorithm.
 *
 * Purpose:   Backward row <i> (0..L) into <dpc>, from row <i>+1 in <dpp> (not read
 * for <i> = L) and the special states of row <i>+1 in <xmx>; writes row
 * <i>'s special states and scale factor into <xmx>. Row L is the
 * initialization. Row 0 holds only the N and B states: its core cells
 * are set to 0 and its N value is returned; 0 is returned for the other
 * rows. <fwd_xmx> holds the Forward special states and scale factors
 * (needed only for p7_BWD_SCALE_FORWARD). If <totscale> is non-NULL it is
 * set to the log of row L's scale factor, and the log of each later
 * scale factor > 1 is added to it. <has_own_scales> is the flag that
 * p7_BWD_SCALE_FORWARD turns on once xB exceeds 1e16. <dpc> and <dpp> may
 * be the same row (a Backward parser's single row): every cell of row
 * <i>+1 is read before row <i>'s cell overwrites it.
 */
float
p7_BackwardRow_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const float *fwd_xmx, float *xmx, float *totscale, int *has_own_scales,
                   __m128 *dpc, const __m128 *dpp, int i, int scale_mode)
{
  if (i == L)
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
      if      (scale_mode == p7_BWD_SCALE_FORWARD) sc = fwd_xmx[L*p7X_NXCELLS+p7X_SCALE];
      else if (scale_mode == p7_BWD_SCALE_OWN)     sc = (xB > 1.0e4f) ? xB : 1.0f;
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
      return 0.0f;
    }

  if (i == 0)
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

      /* row 0 has no core states: its cells are set to 0 */
      for (q = 0; q < Q * p7X_NSCELLS; q++) dpc[q] = _mm_setzero_ps();
      return xN;
    }

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
    if (scale_mode == p7_BWD_SCALE_FORWARD)
      {
        if (xB > 1.0e16f) *has_own_scales = TRUE;
        if (*has_own_scales)
          xmx[i*p7X_NXCELLS+p7X_SCALE] = (xB > 1.0e4f) ? xB : 1.0f;
        else
          xmx[i*p7X_NXCELLS+p7X_SCALE] = fwd_xmx[i*p7X_NXCELLS+p7X_SCALE];
      }
    else if (scale_mode == p7_BWD_SCALE_OWN)
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
  return 0.0f;
}


static int
backward_engine_sse(int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc)
{
  float xN = 0.0f;
  int   i;

  bck->M = om->M;
  bck->L = L;
  bck->has_own_scales = FALSE;	/* backwards scale factors are *usually* given by <fwd> */

  /* rows L down to 0; row 0 returns the N value */
  for (i = L; i >= 0; i--)
    {
      xN = p7_BackwardRow_sse(dsq, L, om, fwd->xmx, bck->xmx, &bck->totscale, &bck->has_own_scales,
                       bck->dpf[i * do_full], (i < L) ? bck->dpf[(i+1) * do_full] : NULL, i, p7_BWD_SCALE_FORWARD);
    }

  if      (isnan(xN))         ESL_EXCEPTION(eslERANGE, "backward score is NaN");
  else if (L > 0 && xN == 0.0f) ESL_EXCEPTION(eslERANGE, "backward score underflow (is 0.0)");
  else if (isinf(xN) == 1)   ESL_EXCEPTION(eslERANGE, "backward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = bck->totscale + log(xN);

  return eslOK;
}

#endif /* eslENABLE_SSE */
