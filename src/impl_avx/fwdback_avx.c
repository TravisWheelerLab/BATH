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
static int backward_engine_avx(int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, P7_OMX *pp, float *opt_sc);


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
  return backward_engine_avx(TRUE, dsq, L, om, fwd, bck, NULL, opt_sc);
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
  return backward_engine_avx(FALSE, dsq, L, om, fwd, bck, NULL, opt_sc);
}

int
p7_BackwardDecoding_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *fwd, P7_OMX *bck, float *opt_sc)
{
#if eslDEBUGLEVEL > 0
  if (om->M >  bck->allocQ4_avx*8)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few columns)");
  if (bck->validR < 1)               ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few MDI rows)");
  if (L     >= bck->allocXR)         ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (too few X rows)");
  if (L     != fwd->L)               ESL_EXCEPTION(eslEINVAL, "fwd matrix size doesn't agree with length L");
  if (L     >= fwd->validR)          ESL_EXCEPTION(eslEINVAL, "fwd matrix isn't a full matrix (too few MDI rows)");
  if (! p7_oprofile_IsLocal(om))     ESL_EXCEPTION(eslEINVAL, "Forward implementation makes assumptions that only work for local alignment");
#endif
  return backward_engine_avx(FALSE, dsq, L, om, fwd, bck, fwd, opt_sc);
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
      for (j = 1; j < 8; j++) {
        dcv = esl_avx_rightshiftz_float(dcv);
        tp  = om->tfv_avx + 7*Q;
        for (q = 0; q < Q; q++) {
          DMO(dpc,q) = _mm256_add_ps(dcv, DMO(dpc,q));
          dcv        = _mm256_mul_ps(dcv, *tp); tp++;
        }
      }
    } else {
      for (j = 1; j < 8; j++) {
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


/* backward_decode_row_avx()
 * For backward_engine_avx() with a posterior matrix <pp>, which is
 * the Forward matrix <fwd>: as soon as Backward row <i> (<dpc>, and the
 * special states in <bck>) is final, write F*B for M and I over row <i>
 * of <pp>, 0 for D, and F(i-1)*B(i)*t_loop for N, J and C. The engine
 * multiplies each row by its scale product after the sweep, completing
 * the same products, in the same order, as p7_Decoding_avx().
 * Posterior row <i> reads Forward row <i> and the Forward special
 * states of row <i-1>, which are not overwritten yet when the rows are
 * done from L down.
 */
static inline void
backward_decode_row_avx(const P7_OPROFILE *om, const P7_OMX *fwd, const P7_OMX *bck, __m256 *dpc, P7_OMX *pp, int i, int Q)
{
  __m256 *ppv = pp->dpf_avx[i];
  __m256 *fv  = fwd->dpf_avx[i];
  int     q;

  for (q = 0; q < Q; q++)
    {
      MMO(ppv,q) = _mm256_mul_ps(MMO(fv,q), MMO(dpc,q));
      DMO(ppv,q) = _mm256_setzero_ps();
      IMO(ppv,q) = _mm256_mul_ps(IMO(fv,q), IMO(dpc,q));
    }
  pp->xmx[i*p7X_NXCELLS+p7X_E] = 0.0;
  pp->xmx[i*p7X_NXCELLS+p7X_N] = fwd->xmx[(i-1)*p7X_NXCELLS+p7X_N] * bck->xmx[i*p7X_NXCELLS+p7X_N] * om->xf[p7O_N][p7O_LOOP];
  pp->xmx[i*p7X_NXCELLS+p7X_J] = fwd->xmx[(i-1)*p7X_NXCELLS+p7X_J] * bck->xmx[i*p7X_NXCELLS+p7X_J] * om->xf[p7O_J][p7O_LOOP];
  pp->xmx[i*p7X_NXCELLS+p7X_C] = fwd->xmx[(i-1)*p7X_NXCELLS+p7X_C] * bck->xmx[i*p7X_NXCELLS+p7X_C] * om->xf[p7O_C][p7O_LOOP];
  pp->xmx[i*p7X_NXCELLS+p7X_B] = 0.0;
}


/* Scale factor a Backward row takes, for backward_row_L_avx() and
 * backward_row_avx():
 *   BWD_SCALE_FORWARD: the Forward matrix's factor for the row, or the
 *     row's own factor once its xB has exceeded 1e16 (p7_Backward());
 *   BWD_SCALE_OWN:     the row's own factor, xB when xB > 1e4, else 1;
 *   BWD_SCALE_STORED:  the factor already stored for the row, when a
 *     row is computed again.
 */
enum bwd_scale_e { BWD_SCALE_FORWARD = 0, BWD_SCALE_OWN = 1, BWD_SCALE_STORED = 2 };


/* backward_row_L_avx()
 * Backward row <L> into <dpc>, special states and scale factor into
 * row <L> of <xmx>. <fwd_xmx> holds the Forward special states and
 * scale factors (needed only for BWD_SCALE_FORWARD). If <totscale> is
 * non-NULL it is set to the log of the row's scale factor.
 */
static void
backward_row_L_avx(const P7_OPROFILE *om, const float *fwd_xmx, float *xmx, float *totscale, __m256 *dpc, int L, int scale_mode)
{
  register __m256 dpv, dcv, xEv;
  __m256  zerov = _mm256_setzero_ps();
  float   xN, xE, xB, xC, xJ, sc;
  int     Q     = p7O_NQF_AVX(om->M);
  int     q, j;
  __m256 *tp;

  xJ     = 0.0f;
  xB     = 0.0f;
  xN     = 0.0f;
  xC     = om->xf[p7O_C][p7O_MOVE];
  xE     = xC * om->xf[p7O_E][p7O_MOVE];
  xEv    = _mm256_set1_ps(xE);
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
  for (j = 1; j < 8; j++) {
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
  if      (scale_mode == BWD_SCALE_FORWARD) sc = fwd_xmx[L*p7X_NXCELLS+p7X_SCALE];
  else if (scale_mode == BWD_SCALE_OWN)     sc = (xB > 1.0e4f) ? xB : 1.0f;
  else                                      sc = xmx[L*p7X_NXCELLS+p7X_SCALE];
  if (sc > 1.0f) {
    xE /= sc;  xN /= sc;  xC /= sc;  xJ /= sc;  xB /= sc;
    xEv = _mm256_set1_ps(1.0f / sc);
    for (q = 0; q < Q; q++) {
      MMO(dpc,q) = _mm256_mul_ps(MMO(dpc,q), xEv);
      DMO(dpc,q) = _mm256_mul_ps(DMO(dpc,q), xEv);
      IMO(dpc,q) = _mm256_mul_ps(IMO(dpc,q), xEv);
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


/* backward_row_avx()
 * Backward row <i> (1 <= i < L) into <dpc> from row <i>+1 in <dpp>,
 * reading the special states of row <i>+1 from <xmx> and writing row
 * <i>'s special states and scale factor into it. <fwd_xmx> holds the
 * Forward special states and scale factors (needed only for
 * BWD_SCALE_FORWARD). If <totscale> is non-NULL the log of a scale
 * factor > 1 is added to it. <has_own_scales> is the flag that
 * BWD_SCALE_FORWARD turns on once xB exceeds 1e16.
 */
static void
backward_row_avx(const ESL_DSQ *dsq, const P7_OPROFILE *om, const float *fwd_xmx, float *xmx, float *totscale, int *has_own_scales,
                 __m256 *dpc, const __m256 *dpp, int i, int scale_mode)
{
  register __m256 mpv, ipv, dpv;
  register __m256 mcv, dcv;
  register __m256 tmmv, timv, tdmv;
  register __m256 xBv;
  register __m256 xEv;
  __m256  zerov = _mm256_setzero_ps();
  float   xN, xE, xB, xC, xJ, sc;
  int     Q     = p7O_NQF_AVX(om->M);
  int     q, j;
  __m256 *rp;
  __m256 *tp;

  xN = xmx[(i+1)*p7X_NXCELLS+p7X_N];
  xJ = xmx[(i+1)*p7X_NXCELLS+p7X_J];
  xC = xmx[(i+1)*p7X_NXCELLS+p7X_C];

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
  for (j = 1; j < 8; j++) {
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
    sc = xmx[i*p7X_NXCELLS+p7X_SCALE];
    xE /= sc;  xN /= sc;  xJ /= sc;  xB /= sc;  xC /= sc;
    xBv = _mm256_set1_ps(1.0f / sc);
    for (q = 0; q < Q; q++) {
      MMO(dpc,q) = _mm256_mul_ps(MMO(dpc,q), xBv);
      DMO(dpc,q) = _mm256_mul_ps(DMO(dpc,q), xBv);
      IMO(dpc,q) = _mm256_mul_ps(IMO(dpc,q), xBv);
    }
    if (totscale != NULL) *totscale += log(sc);
  }

  xmx[i*p7X_NXCELLS+p7X_E] = xE;
  xmx[i*p7X_NXCELLS+p7X_N] = xN;
  xmx[i*p7X_NXCELLS+p7X_J] = xJ;
  xmx[i*p7X_NXCELLS+p7X_B] = xB;
  xmx[i*p7X_NXCELLS+p7X_C] = xC;
}


/* backward_row_0_avx()
 * Backward row 0: only the N and B states, from row 1 in <dpp> and the
 * special states of row 1 in <xmx>. Writes row 0's special states into
 * <xmx> and returns its N value.
 */
static float
backward_row_0_avx(const ESL_DSQ *dsq, const P7_OPROFILE *om, float *xmx, const __m256 *dpp)
{
  register __m256 mpv, xBv;
  float   xN, xB;
  int     Q  = p7O_NQF_AVX(om->M);
  int     q;
  __m256 *rp;
  __m256 *tp;

  xN  = xmx[1*p7X_NXCELLS+p7X_N];
  tp  = om->tfv_avx;
  rp  = om->rfv_avx[dsq[1]];
  xBv = _mm256_setzero_ps();
  for (q = 0; q < Q; q++) {
    mpv = _mm256_mul_ps(MMO(dpp,q), *rp);  rp++;
    mpv = _mm256_mul_ps(mpv,        *tp);  tp += 7;
    xBv = _mm256_add_ps(xBv,        mpv);
  }
  esl_avx_hsum_ps(xBv, &xB);

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
backward_engine_avx(int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, P7_OMX *pp, float *opt_sc)
{
  __m256   zerov = _mm256_setzero_ps();
  float    xN;
  int      i;
  int      q;
  int      Q   = p7O_NQF_AVX(om->M);
  __m256  *dpc;
  __m256  *dpp;

  bck->M = om->M;
  bck->L = L;
  bck->has_own_scales = FALSE;

  dpc = bck->dpf_avx[L * do_full];
  backward_row_L_avx(om, fwd->xmx, bck->xmx, &bck->totscale, dpc, L, BWD_SCALE_FORWARD);
  if (pp != NULL) backward_decode_row_avx(om, fwd, bck, dpc, pp, L, Q);

  /* Main backward recursion */
  for (i = L-1; i >= 1; i--) {
    dpc = bck->dpf_avx[i     * do_full];
    dpp = bck->dpf_avx[(i+1) * do_full];
    backward_row_avx(dsq, om, fwd->xmx, bck->xmx, &bck->totscale, &bck->has_own_scales, dpc, dpp, i, BWD_SCALE_FORWARD);
    if (pp != NULL) backward_decode_row_avx(om, fwd, bck, dpc, pp, i, Q);
  }

  /* Termination at i=0 */
  dpp = bck->dpf_avx[1 * do_full];
  xN  = backward_row_0_avx(dsq, om, bck->xmx, dpp);

  if      (isnan(xN))              ESL_EXCEPTION(eslERANGE, "backward score is NaN");
  else if (L > 0 && xN == 0.0f)   ESL_EXCEPTION(eslERANGE, "backward score underflow (is 0.0)");
  else if (isinf(xN) == 1)         ESL_EXCEPTION(eslERANGE, "backward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = bck->totscale + log(xN);

  /* Posterior decoding: scale each row by its scale product, as p7_Decoding_avx() does */
  if (pp != NULL)
    {
      float   scaleproduct = 1.0 / bck->xmx[p7X_N];
      __m256  totrv;
      __m256 *ppv;

      pp->M = om->M;
      pp->L = L;

      ppv = pp->dpf_avx[0];
      for (q = 0; q < Q; q++)
        MMO(ppv,q) = DMO(ppv,q) = IMO(ppv,q) = zerov;
      pp->xmx[p7X_E] = 0.0;
      pp->xmx[p7X_N] = 0.0;
      pp->xmx[p7X_J] = 0.0;
      pp->xmx[p7X_C] = 0.0;
      pp->xmx[p7X_B] = 0.0;

      for (i = 1; i <= L; i++)
        {
          ppv   = pp->dpf_avx[i];
          totrv = _mm256_set1_ps(scaleproduct * fwd->xmx[i*p7X_NXCELLS+p7X_SCALE]);
          for (q = 0; q < Q; q++)
            {
              MMO(ppv,q) = _mm256_mul_ps(MMO(ppv,q), totrv);
              IMO(ppv,q) = _mm256_mul_ps(IMO(ppv,q), totrv);
            }
          pp->xmx[i*p7X_NXCELLS+p7X_N] *= scaleproduct;
          pp->xmx[i*p7X_NXCELLS+p7X_J] *= scaleproduct;
          pp->xmx[i*p7X_NXCELLS+p7X_C] *= scaleproduct;

          if (bck->has_own_scales) scaleproduct *= fwd->xmx[i*p7X_NXCELLS+p7X_SCALE] /  bck->xmx[i*p7X_NXCELLS+p7X_SCALE];
        }
      if (isinf(scaleproduct)) return eslERANGE;
    }
  return eslOK;
}

#endif /* eslENABLE_AVX */
