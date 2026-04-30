/* Viterbi filter SSE implementation.
 * Ported from impl_sse/vitfilter.c with _sse suffix for runtime dispatch.
 *
 * Contents:
 *    1. p7_ViterbiFilter_sse()     — int16 precision Viterbi score
 *    2. p7_ViterbiFilter_BATH_sse() — score + above-threshold window list
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
#include "esl_gumbel.h"

#include "hmmer.h"
#include "impl_avx.h"

/*****************************************************************
 * 1. p7_ViterbiFilter_sse()
 *****************************************************************/

/* Function:  p7_ViterbiFilter_sse()
 * Synopsis:  Viterbi score in int16 precision, SSE path.
 *
 * Purpose:   Calculate an approximation of the Viterbi score for digital
 *            sequence <dsq> of length <L>, using optimized profile <om> and
 *            one-row DP matrix <ox>.  Returns the estimated Viterbi score
 *            (nats) in <*ret_sc>.
 *
 *            Score may overflow on high-scoring sequences (<eslERANGE>).
 *            The profile must be in local alignment mode.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if score overflows int16 range; <*ret_sc> is
 *            <eslINFINITY>.
 *
 * Throws:    <eslEINVAL> if <ox> is too small or profile is not local.
 */
int
p7_ViterbiFilter_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
#ifdef eslENABLE_SSE
  register __m128i mpv, dpv, ipv;
  register __m128i sv;
  register __m128i dcv;
  register __m128i xEv;
  register __m128i xBv;
  register __m128i Dmaxv;
  int16_t  xE, xB, xC, xJ, xN;
  int16_t  Dmax;
  int i;
  int q;
  int Q       = p7O_NQW(om->M);
  __m128i *dp = ox->dpw[0];
  __m128i *rsc;
  __m128i *tsc;
  __m128i  negInfv;

  if (Q > ox->allocQ8)                                 ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  if (om->mode != p7_LOCAL && om->mode != p7_UNILOCAL) ESL_EXCEPTION(eslEINVAL, "Fast filter only works for local alignment");
  ox->M = om->M;

  negInfv = _mm_set1_epi16(-32768);
  negInfv = _mm_srli_si128(negInfv, 14);   /* lowest 2 bytes = -32768, rest 0 */

  for (q = 0; q < Q; q++)
    MMXo(q) = IMXo(q) = DMXo(q) = _mm_set1_epi16(-32768);
  xN = om->base_w;
  xB = xN + om->xw[p7O_N][p7O_MOVE];
  xJ = -32768;
  xC = -32768;
  xE = -32768;

#if eslDEBUGLEVEL > 0
  if (ox->debugging) p7_omx_DumpVFRow(ox, 0, xE, 0, xJ, xB, xC);
#endif

  for (i = 1; i <= L; i++) {
    rsc   = om->rwv[dsq[i]];
    tsc   = om->twv;
    dcv   = _mm_set1_epi16(-32768);
    xEv   = _mm_set1_epi16(-32768);
    Dmaxv = _mm_set1_epi16(-32768);
    xBv   = _mm_set1_epi16(xB);

    mpv = MMXo(Q-1);  mpv = _mm_slli_si128(mpv, 2);  mpv = _mm_or_si128(mpv, negInfv);
    dpv = DMXo(Q-1);  dpv = _mm_slli_si128(dpv, 2);  dpv = _mm_or_si128(dpv, negInfv);
    ipv = IMXo(Q-1);  ipv = _mm_slli_si128(ipv, 2);  ipv = _mm_or_si128(ipv, negInfv);

    for (q = 0; q < Q; q++) {
      sv   =                    _mm_adds_epi16(xBv, *tsc);  tsc++;
      sv   = _mm_max_epi16(sv,  _mm_adds_epi16(mpv, *tsc)); tsc++;
      sv   = _mm_max_epi16(sv,  _mm_adds_epi16(ipv, *tsc)); tsc++;
      sv   = _mm_max_epi16(sv,  _mm_adds_epi16(dpv, *tsc)); tsc++;
      sv   = _mm_adds_epi16(sv, *rsc);                      rsc++;
      xEv  = _mm_max_epi16(xEv, sv);

      mpv     = MMXo(q);
      dpv     = DMXo(q);
      ipv     = IMXo(q);

      MMXo(q) = sv;
      DMXo(q) = dcv;

      dcv   = _mm_adds_epi16(sv, *tsc);  tsc++;
      Dmaxv = _mm_max_epi16(dcv, Dmaxv);

      sv      =                    _mm_adds_epi16(mpv, *tsc); tsc++;
      IMXo(q) = _mm_max_epi16(sv,  _mm_adds_epi16(ipv, *tsc)); tsc++;
    }

    xE = esl_sse_hmax_epi16(xEv);
    if (xE >= 32767) { *ret_sc = eslINFINITY; return eslERANGE; }

    xN = xN +  om->xw[p7O_N][p7O_LOOP];
    xC = ESL_MAX(xC + om->xw[p7O_C][p7O_LOOP], xE + om->xw[p7O_E][p7O_MOVE]);
    xJ = ESL_MAX(xJ + om->xw[p7O_J][p7O_LOOP], xE + om->xw[p7O_E][p7O_LOOP]);
    xB = ESL_MAX(xJ + om->xw[p7O_J][p7O_MOVE], xN + om->xw[p7O_N][p7O_MOVE]);

    /* Lazy F loop */
    Dmax = esl_sse_hmax_epi16(Dmaxv);
    if (Dmax + om->ddbound_w > xB) {
      dcv = _mm_slli_si128(dcv, 2);
      dcv = _mm_or_si128(dcv, negInfv);
      tsc = om->twv + 7*Q;
      for (q = 0; q < Q; q++) {
        DMXo(q) = _mm_max_epi16(dcv, DMXo(q));
        dcv     = _mm_adds_epi16(DMXo(q), *tsc); tsc++;
      }
      do {
        dcv = _mm_slli_si128(dcv, 2);
        dcv = _mm_or_si128(dcv, negInfv);
        tsc = om->twv + 7*Q;
        for (q = 0; q < Q; q++) {
          if (!esl_sse_any_gt_epi16(dcv, DMXo(q))) break;
          DMXo(q) = _mm_max_epi16(dcv, DMXo(q));
          dcv     = _mm_adds_epi16(DMXo(q), *tsc); tsc++;
        }
      } while (q == Q);
    } else {
      dcv = _mm_slli_si128(dcv, 2);
      DMXo(0) = _mm_or_si128(dcv, negInfv);
    }

#if eslDEBUGLEVEL > 0
    if (ox->debugging) p7_omx_DumpVFRow(ox, i, xE, 0, xJ, xB, xC);
#endif
  }

  if (xC > -32768) {
    *ret_sc  = (float) xC + (float) om->xw[p7O_C][p7O_MOVE] - (float) om->base_w;
    *ret_sc /= om->scale_w;
    *ret_sc -= 3.0f;
  } else {
    *ret_sc = -eslINFINITY;
  }
  return eslOK;
#else
  return eslENORESULT;
#endif
}
/*---------------- end, p7_ViterbiFilter_sse() ------------------*/


/*****************************************************************
 * 2. p7_ViterbiFilter_BATH_sse()
 *****************************************************************/

/* Function:  p7_ViterbiFilter_BATH_sse()
 * Synopsis:  Viterbi filter + above-threshold window list, SSE path.
 *
 * Purpose:   Run the int16 Viterbi DP without dp-matrix reset, accumulating
 *            the global Viterbi score in <*ret_sc>.  Additionally, at each
 *            row i where the row-maximum xE >= <sc_thresh>, identify the model
 *            position k_start, extend a forward diagonal in SSV score space,
 *            and append a window to <windowlist>.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if score overflows; <*ret_sc> is <eslINFINITY>.
 *
 * Throws:    <eslEINVAL> if <ox> is too small or profile is not local.
 */
int
p7_ViterbiFilter_BATH_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox,
                           const P7_SCOREDATA *ssvdata, float filtersc, double P,
                           P7_HMM_WINDOWLIST *windowlist, float *ret_sc)
{
#ifdef eslENABLE_SSE
  register __m128i mpv, dpv, ipv;
  register __m128i sv;
  register __m128i dcv;
  register __m128i xEv;
  register __m128i xBv;
  register __m128i Dmaxv;
  __m128i  negInfv;
  int16_t  xE, xB, xC, xJ, xN;
  int16_t  Dmax;
  int      i;
  int      q;
  int      Q        = p7O_NQW(om->M);
  __m128i *dp       = ox->dpw[0];
  __m128i *rsc;
  __m128i *tsc;

  int16_t sc_thresh;
  int     sc_ext_thresh;
  float   invP;
  int     z, k;
  int     skip_until = 0;
  union { __m128i v; int16_t i[8]; } tmp;

  invP      = esl_gumbel_invsurv(P, om->evparam[p7_VMU], om->evparam[p7_VLAMBDA]);
  sc_thresh = (int16_t) ceil(((filtersc + eslCONST_LOG2 * invP + 3.0) * om->scale_w)
              - (float) om->xw[p7O_E][p7O_MOVE] - (float) om->xw[p7O_C][p7O_MOVE]
              + (float) om->base_w);

  invP          = esl_gumbel_invsurv(P, om->evparam[p7_MMU], om->evparam[p7_MLAMBDA]);
  sc_ext_thresh = (int) ceil(((filtersc + eslCONST_LOG2 * invP + 3.0) * om->scale_b)
                  + om->base_b + om->tec_b + om->tjb_b);

  if (Q > ox->allocQ8)                                 ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  if (om->mode != p7_LOCAL && om->mode != p7_UNILOCAL) ESL_EXCEPTION(eslEINVAL, "Fast filter only works for local alignment");
  ox->M = om->M;

  negInfv = _mm_set1_epi16(-32768);
  negInfv = _mm_srli_si128(negInfv, 14);

  for (q = 0; q < Q; q++)
    MMXo(q) = IMXo(q) = DMXo(q) = _mm_set1_epi16(-32768);
  xN = om->base_w;
  xB = xN + om->xw[p7O_N][p7O_MOVE];
  xJ = -32768;
  xC = -32768;
  xE = -32768;

  for (i = 1; i <= L; i++) {
    rsc   = om->rwv[dsq[i]];
    tsc   = om->twv;
    dcv   = _mm_set1_epi16(-32768);
    xEv   = _mm_set1_epi16(-32768);
    Dmaxv = _mm_set1_epi16(-32768);
    xBv   = _mm_set1_epi16(xB);

    mpv = MMXo(Q-1);  mpv = _mm_slli_si128(mpv, 2);  mpv = _mm_or_si128(mpv, negInfv);
    dpv = DMXo(Q-1);  dpv = _mm_slli_si128(dpv, 2);  dpv = _mm_or_si128(dpv, negInfv);
    ipv = IMXo(Q-1);  ipv = _mm_slli_si128(ipv, 2);  ipv = _mm_or_si128(ipv, negInfv);

    for (q = 0; q < Q; q++) {
      sv   =                    _mm_adds_epi16(xBv, *tsc);  tsc++;
      sv   = _mm_max_epi16(sv,  _mm_adds_epi16(mpv, *tsc)); tsc++;
      sv   = _mm_max_epi16(sv,  _mm_adds_epi16(ipv, *tsc)); tsc++;
      sv   = _mm_max_epi16(sv,  _mm_adds_epi16(dpv, *tsc)); tsc++;
      sv   = _mm_adds_epi16(sv, *rsc);                      rsc++;
      xEv  = _mm_max_epi16(xEv, sv);

      mpv     = MMXo(q);
      dpv     = DMXo(q);
      ipv     = IMXo(q);

      MMXo(q) = sv;
      DMXo(q) = dcv;

      dcv   = _mm_adds_epi16(sv, *tsc);  tsc++;
      Dmaxv = _mm_max_epi16(dcv, Dmaxv);

      sv      =                    _mm_adds_epi16(mpv, *tsc); tsc++;
      IMXo(q) = _mm_max_epi16(sv,  _mm_adds_epi16(ipv, *tsc)); tsc++;
    }

    xE = esl_sse_hmax_epi16(xEv);
    if (xE >= 32767) { *ret_sc = eslINFINITY; return eslERANGE; }

    xN = xN +  om->xw[p7O_N][p7O_LOOP];
    xC = ESL_MAX(xC + om->xw[p7O_C][p7O_LOOP], xE + om->xw[p7O_E][p7O_MOVE]);
    xJ = ESL_MAX(xJ + om->xw[p7O_J][p7O_LOOP], xE + om->xw[p7O_E][p7O_LOOP]);
    xB = ESL_MAX(xJ + om->xw[p7O_J][p7O_MOVE], xN + om->xw[p7O_N][p7O_MOVE]);

    if (i > skip_until && xE >= sc_thresh) {
      /* Find k_start: model position where MMXo score == xE */
      int k_start = 0;
      for (q = 0; q < Q && k_start == 0; q++) {
        tmp.v = MMXo(q);
        for (z = 0; z < 8; z++) {
          k = q + Q*z + 1;
          if (tmp.i[z] == xE && k <= om->M) { k_start = k; break; }
        }
      }

      /* Forward diagonal extension in SSV score space */
      int max_k_end     = k_start;
      int max_i_end     = i;
      int sc_ext        = sc_ext_thresh;
      int max_sc_ext    = sc_ext;
      int pos_since_max = 0;
      int kk = k_start + 1;
      int nn = i + 1;
      while (kk <= om->M && nn <= L) {
        sc_ext += om->bias_b - ssvdata->ssv_scores[kk * om->abc->Kp + dsq[nn]];
        if (sc_ext >= max_sc_ext) {
          max_sc_ext    = sc_ext;
          max_k_end     = kk;
          max_i_end     = nn;
          pos_since_max = 0;
        } else {
          if (++pos_since_max == 5) break;
        }
        kk++;
        nn++;
      }

      p7_hmmwindow_new(windowlist, 0, i, max_k_end, max_k_end - k_start + 1, 0.0f, p7_NOCOMPLEMENT, L);
      skip_until = max_i_end;
    }

    /* Lazy F loop */
    Dmax = esl_sse_hmax_epi16(Dmaxv);
    if (Dmax + om->ddbound_w > xB) {
      dcv = _mm_slli_si128(dcv, 2);
      dcv = _mm_or_si128(dcv, negInfv);
      tsc = om->twv + 7*Q;
      for (q = 0; q < Q; q++) {
        DMXo(q) = _mm_max_epi16(dcv, DMXo(q));
        dcv     = _mm_adds_epi16(DMXo(q), *tsc); tsc++;
      }
      do {
        dcv = _mm_slli_si128(dcv, 2);
        dcv = _mm_or_si128(dcv, negInfv);
        tsc = om->twv + 7*Q;
        for (q = 0; q < Q; q++) {
          if (!esl_sse_any_gt_epi16(dcv, DMXo(q))) break;
          DMXo(q) = _mm_max_epi16(dcv, DMXo(q));
          dcv     = _mm_adds_epi16(DMXo(q), *tsc); tsc++;
        }
      } while (q == Q);
    } else {
      dcv = _mm_slli_si128(dcv, 2);
      DMXo(0) = _mm_or_si128(dcv, negInfv);
    }
  }

  if (xC > -32768) {
    *ret_sc  = (float) xC + (float) om->xw[p7O_C][p7O_MOVE] - (float) om->base_w;
    *ret_sc /= om->scale_w;
    *ret_sc -= 3.0f;
  } else {
    *ret_sc = -eslINFINITY;
  }
  return eslOK;
#else
  return eslENORESULT;
#endif
}
/*---------------- end, p7_ViterbiFilter_BATH_sse() -------------*/
