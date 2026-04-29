/* Viterbi filter AVX2 implementation.
 * Ported from vitfilter_sse.c using 256-bit vectors.
 *
 * Contents:
 *    1. p7_ViterbiFilter_avx()      — int16 precision Viterbi score
 *    2. p7_ViterbiFilter_BATH_avx() — score + above-threshold window list
 */

#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <stdio.h>
#include <math.h>

#include <immintrin.h>

#include "easel.h"
#include "esl_avx.h"
#include "esl_gumbel.h"

#include "hmmer.h"
#include "impl_avx.h"


/*****************************************************************
 * 1. p7_ViterbiFilter_avx()
 *****************************************************************/

/* Function:  p7_ViterbiFilter_avx()
 * Synopsis:  Viterbi score in int16 precision, AVX2 path.
 *
 * Purpose:   Calculate an approximation of the Viterbi score for digital
 *            sequence <dsq> of length <L>, using optimized profile <om> and
 *            one-row DP matrix <ox>.  Returns the estimated Viterbi score
 *            (nats) in <*ret_sc>.
 *
 *            Uses 256-bit AVX2 vectors: Q = ceil(M/16) stripes of 16
 *            int16 lanes each.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if score overflows int16 range; <*ret_sc> is
 *            <eslINFINITY>.
 *
 * Throws:    <eslEINVAL> if <ox> is too small or profile is not local.
 */
int
p7_ViterbiFilter_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register __m256i mpv, dpv, ipv;
  register __m256i sv;
  register __m256i dcv;
  register __m256i xEv;
  register __m256i xBv;
  register __m256i Dmaxv;
  int16_t  xE, xB, xC, xJ, xN;
  int16_t  Dmax;
  int i;
  int q;
  int Q        = p7O_NQW_AVX(om->M);
  __m256i *dp  = ox->dpw_avx[0];
  __m256i *rsc;
  __m256i *tsc;
  __m256i  negInfv;

  if (Q > ox->allocQ8_avx)                               ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  if (om->mode != p7_LOCAL && om->mode != p7_UNILOCAL)   ESL_EXCEPTION(eslEINVAL, "Fast filter only works for local alignment");
  ox->M = om->M;

  /* -32768 only at element 0; used to fill right-shifted-in position */
  negInfv = _mm256_set_epi16(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,(int16_t)-32768);

  for (q = 0; q < Q; q++)
    MMXo(q) = IMXo(q) = DMXo(q) = _mm256_set1_epi16(-32768);
  xN = om->base_w;
  xB = xN + om->xw[p7O_N][p7O_MOVE];
  xJ = -32768;
  xC = -32768;
  xE = -32768;

  for (i = 1; i <= L; i++) {
    rsc   = om->rwv_avx[dsq[i]];
    tsc   = om->twv_avx;
    dcv   = _mm256_set1_epi16(-32768);
    xEv   = _mm256_set1_epi16(-32768);
    Dmaxv = _mm256_set1_epi16(-32768);
    xBv   = _mm256_set1_epi16(xB);

    mpv = esl_avx_rightshift_int16(MMXo(Q-1), negInfv);
    dpv = esl_avx_rightshift_int16(DMXo(Q-1), negInfv);
    ipv = esl_avx_rightshift_int16(IMXo(Q-1), negInfv);

    for (q = 0; q < Q; q++) {
      sv   =                     _mm256_adds_epi16(xBv, *tsc);  tsc++;
      sv   = _mm256_max_epi16(sv, _mm256_adds_epi16(mpv, *tsc)); tsc++;
      sv   = _mm256_max_epi16(sv, _mm256_adds_epi16(ipv, *tsc)); tsc++;
      sv   = _mm256_max_epi16(sv, _mm256_adds_epi16(dpv, *tsc)); tsc++;
      sv   = _mm256_adds_epi16(sv, *rsc);                        rsc++;
      xEv  = _mm256_max_epi16(xEv, sv);

      mpv     = MMXo(q);
      dpv     = DMXo(q);
      ipv     = IMXo(q);

      MMXo(q) = sv;
      DMXo(q) = dcv;

      dcv   = _mm256_adds_epi16(sv, *tsc);  tsc++;
      Dmaxv = _mm256_max_epi16(dcv, Dmaxv);

      sv      =                     _mm256_adds_epi16(mpv, *tsc); tsc++;
      IMXo(q) = _mm256_max_epi16(sv, _mm256_adds_epi16(ipv, *tsc)); tsc++;
    }

    xE = esl_avx_hmax_epi16(xEv);
    if (xE >= 32767) { *ret_sc = eslINFINITY; return eslERANGE; }

    xN = xN +  om->xw[p7O_N][p7O_LOOP];
    xC = ESL_MAX(xC + om->xw[p7O_C][p7O_LOOP], xE + om->xw[p7O_E][p7O_MOVE]);
    xJ = ESL_MAX(xJ + om->xw[p7O_J][p7O_LOOP], xE + om->xw[p7O_E][p7O_LOOP]);
    xB = ESL_MAX(xJ + om->xw[p7O_J][p7O_MOVE], xN + om->xw[p7O_N][p7O_MOVE]);

    /* Lazy F loop */
    Dmax = esl_avx_hmax_epi16(Dmaxv);
    if (Dmax + om->ddbound_w > xB) {
      dcv = esl_avx_rightshift_int16(dcv, negInfv);
      tsc = om->twv_avx + 7*Q;
      for (q = 0; q < Q; q++) {
        DMXo(q) = _mm256_max_epi16(dcv, DMXo(q));
        dcv     = _mm256_adds_epi16(DMXo(q), *tsc); tsc++;
      }
      do {
        dcv = esl_avx_rightshift_int16(dcv, negInfv);
        tsc = om->twv_avx + 7*Q;
        for (q = 0; q < Q; q++) {
          if (!esl_avx_any_gt_epi16(dcv, DMXo(q))) break;
          DMXo(q) = _mm256_max_epi16(dcv, DMXo(q));
          dcv     = _mm256_adds_epi16(DMXo(q), *tsc); tsc++;
        }
      } while (q == Q);
    } else {
      DMXo(0) = esl_avx_rightshift_int16(dcv, negInfv);
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
}
/*---------------- end, p7_ViterbiFilter_avx() ------------------*/


/*****************************************************************
 * 2. p7_ViterbiFilter_BATH_avx()
 *****************************************************************/

/* Function:  p7_ViterbiFilter_BATH_avx()
 * Synopsis:  Viterbi filter + above-threshold window list, AVX2 path.
 *
 * Purpose:   Run the int16 Viterbi DP without dp-matrix reset, accumulating
 *            the global Viterbi score in <*ret_sc>.  Additionally, at each
 *            row i where the row-maximum xE >= <sc_thresh>, identify the model
 *            position k_start, extend a forward diagonal in SSV score space,
 *            and append a window to <windowlist>.
 *
 *            Uses 256-bit AVX2 vectors: Q = ceil(M/16) stripes of 16
 *            int16 lanes each.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if score overflows; <*ret_sc> is <eslINFINITY>.
 *
 * Throws:    <eslEINVAL> if <ox> is too small or profile is not local.
 */
int
p7_ViterbiFilter_BATH_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox,
                           const P7_SCOREDATA *ssvdata, float filtersc, double P,
                           P7_HMM_WINDOWLIST *windowlist, float *ret_sc)
{
  register __m256i mpv, dpv, ipv;
  register __m256i sv;
  register __m256i dcv;
  register __m256i xEv;
  register __m256i xBv;
  register __m256i Dmaxv;
  __m256i  negInfv;
  int16_t  xE, xB, xC, xJ, xN;
  int16_t  Dmax;
  int      i;
  int      q;
  int      Q        = p7O_NQW_AVX(om->M);
  __m256i *dp       = ox->dpw_avx[0];
  __m256i *rsc;
  __m256i *tsc;

  int16_t sc_thresh;
  int     sc_ext_thresh;
  float   invP;
  int     z, k;
  int     skip_until = 0;
  union { __m256i v; int16_t i[16]; } tmp;

  invP      = esl_gumbel_invsurv(P, om->evparam[p7_VMU], om->evparam[p7_VLAMBDA]);
  sc_thresh = (int16_t) ceil(((filtersc + eslCONST_LOG2 * invP + 3.0) * om->scale_w)
              - (float) om->xw[p7O_E][p7O_MOVE] - (float) om->xw[p7O_C][p7O_MOVE]
              + (float) om->base_w);

  invP          = esl_gumbel_invsurv(P, om->evparam[p7_MMU], om->evparam[p7_MLAMBDA]);
  sc_ext_thresh = (int) ceil(((filtersc + eslCONST_LOG2 * invP + 3.0) * om->scale_b)
                  + om->base_b + om->tec_b + om->tjb_b);

  if (Q > ox->allocQ8_avx)                               ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  if (om->mode != p7_LOCAL && om->mode != p7_UNILOCAL)   ESL_EXCEPTION(eslEINVAL, "Fast filter only works for local alignment");
  ox->M = om->M;

  negInfv = _mm256_set_epi16(0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,(int16_t)-32768);

  for (q = 0; q < Q; q++)
    MMXo(q) = IMXo(q) = DMXo(q) = _mm256_set1_epi16(-32768);
  xN = om->base_w;
  xB = xN + om->xw[p7O_N][p7O_MOVE];
  xJ = -32768;
  xC = -32768;
  xE = -32768;

  for (i = 1; i <= L; i++) {
    rsc   = om->rwv_avx[dsq[i]];
    tsc   = om->twv_avx;
    dcv   = _mm256_set1_epi16(-32768);
    xEv   = _mm256_set1_epi16(-32768);
    Dmaxv = _mm256_set1_epi16(-32768);
    xBv   = _mm256_set1_epi16(xB);

    mpv = esl_avx_rightshift_int16(MMXo(Q-1), negInfv);
    dpv = esl_avx_rightshift_int16(DMXo(Q-1), negInfv);
    ipv = esl_avx_rightshift_int16(IMXo(Q-1), negInfv);

    for (q = 0; q < Q; q++) {
      sv   =                     _mm256_adds_epi16(xBv, *tsc);  tsc++;
      sv   = _mm256_max_epi16(sv, _mm256_adds_epi16(mpv, *tsc)); tsc++;
      sv   = _mm256_max_epi16(sv, _mm256_adds_epi16(ipv, *tsc)); tsc++;
      sv   = _mm256_max_epi16(sv, _mm256_adds_epi16(dpv, *tsc)); tsc++;
      sv   = _mm256_adds_epi16(sv, *rsc);                        rsc++;
      xEv  = _mm256_max_epi16(xEv, sv);

      mpv     = MMXo(q);
      dpv     = DMXo(q);
      ipv     = IMXo(q);

      MMXo(q) = sv;
      DMXo(q) = dcv;

      dcv   = _mm256_adds_epi16(sv, *tsc);  tsc++;
      Dmaxv = _mm256_max_epi16(dcv, Dmaxv);

      sv      =                     _mm256_adds_epi16(mpv, *tsc); tsc++;
      IMXo(q) = _mm256_max_epi16(sv, _mm256_adds_epi16(ipv, *tsc)); tsc++;
    }

    xE = esl_avx_hmax_epi16(xEv);
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
        for (z = 0; z < 16; z++) {
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
    Dmax = esl_avx_hmax_epi16(Dmaxv);
    if (Dmax + om->ddbound_w > xB) {
      dcv = esl_avx_rightshift_int16(dcv, negInfv);
      tsc = om->twv_avx + 7*Q;
      for (q = 0; q < Q; q++) {
        DMXo(q) = _mm256_max_epi16(dcv, DMXo(q));
        dcv     = _mm256_adds_epi16(DMXo(q), *tsc); tsc++;
      }
      do {
        dcv = esl_avx_rightshift_int16(dcv, negInfv);
        tsc = om->twv_avx + 7*Q;
        for (q = 0; q < Q; q++) {
          if (!esl_avx_any_gt_epi16(dcv, DMXo(q))) break;
          DMXo(q) = _mm256_max_epi16(dcv, DMXo(q));
          dcv     = _mm256_adds_epi16(DMXo(q), *tsc); tsc++;
        }
      } while (q == Q);
    } else {
      DMXo(0) = esl_avx_rightshift_int16(dcv, negInfv);
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
}
/*---------------- end, p7_ViterbiFilter_BATH_avx() -------------*/

#endif /* eslENABLE_AVX */
