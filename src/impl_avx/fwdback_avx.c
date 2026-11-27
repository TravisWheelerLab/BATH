/* Forward and Backward algorithms, AVX2 implementation.
 * Ported from fwdback_sse.c using 256-bit vectors.
 *
 * Contents:
 *    1. p7_Forward_avx, p7_ForwardParser_avx, p7_Backward_avx, p7_BackwardParser_avx
 *    2. Row functions p7_ForwardRow_avx, p7_BackwardRow_avx, and the static
 *       engines (forward_engine_avx, backward_engine_avx)
 *
 * The checkpointed Forward and Backward are in fwdback_chk_avx.c.
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

/* Function:  p7_ForwardRow_avx()
 * Synopsis:  One row of the AVX Forward algorithm.
 *
 * Purpose:   Forward row <i> of profile <om> against <dsq> into <dpc>,
 *            from row i-1 in <dpp> (not used for <i> = 0), reading the
 *            special states of row i-1 from <xmx> and writing row <i>'s
 *            special states and scale factor into it. When the row is
 *            rescaled, the log of its scale factor is added to
 *            <*totscale>. Row 0 is the initialization. Same calculation
 *            as each row of <p7_Forward_avx()>, which calls it.
 */
void
p7_ForwardRow_avx(const ESL_DSQ *dsq, const P7_OPROFILE *om, float *xmx, float *totscale, __m256 *dpc, const __m256 *dpp, int i)
{
  register __m256 mpv, dpv, ipv;
  register __m256 sv;
  register __m256 dcv;
  register __m256 xEv;
  register __m256 xBv;
  __m256   zerov = _mm256_setzero_ps();
  float    xN, xE, xB, xC, xJ;
  int      q;
  int      j;
  int      Q     = p7O_NQF_AVX(om->M);
  __m256  *rp;
  __m256  *tp;

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
forward_engine_avx(int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *opt_sc)
{
  float    xC;
  int      i;
  __m256  *dpc = ox->dpf_avx[0];
  __m256  *dpp;

  ox->M  = om->M;
  ox->L  = L;
  ox->has_own_scales = TRUE;
  ox->totscale       = 0.0;
  p7_ForwardRow_avx(dsq, om, ox->xmx, &ox->totscale, dpc, NULL, 0);

  for (i = 1; i <= L; i++) {
    dpp = dpc;
    dpc = ox->dpf_avx[do_full * i];
    p7_ForwardRow_avx(dsq, om, ox->xmx, &ox->totscale, dpc, dpp, i);
  }
  xC = ox->xmx[L*p7X_NXCELLS+p7X_C];

  if      (isnan(xC))             ESL_EXCEPTION(eslERANGE, "forward score is NaN");
  else if (L > 0 && xC == 0.0f)  ESL_EXCEPTION(eslERANGE, "forward score underflow (is 0.0)");
  else if (isinf(xC) == 1)        ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = ox->totscale + log(xC * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}


/* Function:  p7_BackwardRow_avx()
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
p7_BackwardRow_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const float *fwd_xmx, float *xmx, float *totscale, int *has_own_scales,
                   __m256 *dpc, const __m256 *dpp, int i, int scale_mode)
{
  if (i == L)
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
      if      (scale_mode == p7_BWD_SCALE_FORWARD) sc = fwd_xmx[L*p7X_NXCELLS+p7X_SCALE];
      else if (scale_mode == p7_BWD_SCALE_OWN)     sc = (xB > 1.0e4f) ? xB : 1.0f;
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
      return 0.0f;
    }

  if (i == 0)
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

      /* row 0 has no core states: its cells are set to 0 */
      for (q = 0; q < Q * p7X_NSCELLS; q++) dpc[q] = _mm256_setzero_ps();
      return xN;
    }

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
  return 0.0f;
}


static int
backward_engine_avx(int do_full, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc)
{
  float    xN = 0.0f;
  int      i;

  bck->M = om->M;
  bck->L = L;
  bck->has_own_scales = FALSE;

  /* rows L down to 0; row 0 returns the N value */
  for (i = L; i >= 0; i--)
    xN = p7_BackwardRow_avx(dsq, L, om, fwd->xmx, bck->xmx, &bck->totscale, &bck->has_own_scales,
                          bck->dpf_avx[i * do_full], (i < L) ? bck->dpf_avx[(i+1) * do_full] : NULL, i, p7_BWD_SCALE_FORWARD);

  if      (isnan(xN))              ESL_EXCEPTION(eslERANGE, "backward score is NaN");
  else if (L > 0 && xN == 0.0f)   ESL_EXCEPTION(eslERANGE, "backward score underflow (is 0.0)");
  else if (isinf(xN) == 1)         ESL_EXCEPTION(eslERANGE, "backward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = bck->totscale + log(xN);

  return eslOK;
}

#endif /* eslENABLE_AVX */
