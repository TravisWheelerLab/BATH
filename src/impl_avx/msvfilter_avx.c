/* MSV filter AVX2 implementation.
 * Ported from msvfilter_sse.c using 256-bit vectors.
 *
 * Contents:
 *    1. p7_MSVFilter_avx()       — fast byte-precision MSV score
 *    2. p7_SSVFilter_BATH_avx()  — BATH ungapped window finder
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
 * 1. p7_MSVFilter_avx()
 *****************************************************************/

/* Function:  p7_MSVFilter_avx()
 * Synopsis:  MSV score in uint8 precision, AVX2 path.
 *
 * Purpose:   Calculate an approximation of the MSV score for digital
 *            sequence <dsq> of length <L>, using optimized profile <om>
 *            and pre-allocated one-row DP matrix <ox>.  Returns the
 *            estimated MSV score (nats) in <*ret_sc>.
 *
 *            Tries p7_SSVFilter_avx() first; falls back to the full MSV
 *            calculation only when SSV returns <eslENORESULT>.
 *
 *            Uses 256-bit AVX2 vectors: Q = ceil(M/32) stripes of 32
 *            uint8 lanes each.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if the score overflows uint8 range (high-scoring hit).
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_MSVFilter_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register __m256i mpv;
  register __m256i xEv;
  register __m256i xBv;
  register __m256i sv;
  register __m256i biasv;
  uint8_t  xJ;
  int i;
  int q;
  int Q         = p7O_NQB_AVX(om->M);
  __m256i *dp   = ox->dpb_avx[0];
  __m256i *rsc;

  __m256i xJv;
  __m256i tjbmv;
  __m256i tecv;
  __m256i basev;
  __m256i ceilingv;
  __m256i zerov;
  __m256i tempv;

  int cmp;
  int status = eslOK;

  if (Q > ox->allocQ16_avx) ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M = om->M;

  /* Try highly optimized SSV filter first */
  status = p7_SSVFilter_avx(dsq, L, om, ret_sc);
  if (status != eslENORESULT) return status;

  /* Initialization. In offset unsigned arithmetic, -infinity is 0,
   * and 0 is om->base_b. */
  biasv    = _mm256_set1_epi8((int8_t) om->bias_b);
  zerov    = _mm256_setzero_si256();
  for (q = 0; q < Q; q++) dp[q] = zerov;
  xJ       = 0;

  ceilingv = _mm256_cmpeq_epi8(biasv, biasv);
  basev    = _mm256_set1_epi8((int8_t) om->base_b);
  tjbmv    = _mm256_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);
  tecv     = _mm256_set1_epi8((int8_t) om->tec_b);

  xJv = _mm256_subs_epu8(biasv, biasv);
  xBv = _mm256_subs_epu8(basev, tjbmv);

  for (i = 1; i <= L; i++) {
    rsc = om->rbv_avx[dsq[i]];
    xEv = zerov;

    /* Right-shift by 1 byte: shift all 32 bytes right, zero fills. */
    mpv = esl_avx_rightshift_int8(dp[Q-1], zerov);
    for (q = 0; q < Q; q++) {
      sv   = _mm256_max_epu8(mpv, xBv);
      sv   = _mm256_adds_epu8(sv, biasv);
      sv   = _mm256_subs_epu8(sv, *rsc);  rsc++;
      xEv  = _mm256_max_epu8(xEv, sv);

      mpv   = dp[q];
      dp[q] = sv;
    }

    /* Overflow test */
    tempv = _mm256_adds_epu8(xEv, biasv);
    tempv = _mm256_cmpeq_epi8(tempv, ceilingv);
    cmp   = _mm256_movemask_epi8(tempv);

    if (cmp != 0x00000000) {
      *ret_sc = eslINFINITY;
      return eslERANGE;
    }

    /* Horizontal max — broadcast back into xEv */
    xJ  = esl_avx_hmax_epu8(xEv);
    xEv = _mm256_set1_epi8((int8_t) xJ);

    xEv = _mm256_subs_epu8(xEv, tecv);
    xJv = _mm256_max_epu8(xJv, xEv);
    xBv = _mm256_max_epu8(basev, xJv);
    xBv = _mm256_subs_epu8(xBv, tjbmv);
  }

  xJ = esl_avx_hmax_epu8(xJv);

  *ret_sc  = ((float)(xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0f;

  return eslOK;
}
/*------------------ end, p7_MSVFilter_avx() ----------------------*/


/*****************************************************************
 * 2. p7_SSVFilter_BATH_avx()
 *****************************************************************/

/* Function:  p7_SSVFilter_BATH_avx()
 * Synopsis:  BATH ungapped window finder, AVX2 path.
 *
 * Purpose:   SSV scan over <dsq[1..L]> using profile <om> and
 *            one-row DP matrix <ox>.  Appends to <windowlist> every
 *            ungapped diagonal whose score exceeds the threshold
 *            implied by <P>, <bg>, and the MSV EVD parameters.
 *
 *            Uses 256-bit AVX2 vectors: Q = ceil(M/32) stripes of 32
 *            uint8 lanes each.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_SSVFilter_BATH_avx(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox,
                      const P7_SCOREDATA *ssvdata, P7_BG *bg, double P,
                      P7_HMM_WINDOWLIST *windowlist)
{
  register __m256i mpv;
  register __m256i xEv;
  register __m256i xBv;
  register __m256i sv;
  register __m256i biasv;
  int i;
  int q;
  int Q         = p7O_NQB_AVX(om->M);
  __m256i *dp   = ox->dpb_avx[0];
  __m256i *rsc;
  __m256i tjbmv;
  __m256i basev;
  __m256i ceilingv;
  __m256i zerov;
  __m256i tempv;
  int cmp;
  int k, n, end, rem_sc, start, target_end, target_start;
  int max_end, max_sc, sc, pos_since_max;
  float nullsc;
  float ret_sc;

  union { __m256i v; uint8_t b[32]; } u;

  __m256i sc_threshv;
  uint8_t sc_thresh;
  float invP = esl_gumbel_invsurv(P, om->evparam[p7_MMU], om->evparam[p7_MLAMBDA]);

  if (Q > ox->allocQ16_avx) ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M = om->M;

  p7_bg_SetLength(bg, L);
  p7_oprofile_ReconfigMSVLength(om, L);
  p7_bg_NullOne(bg, dsq, L, &nullsc);

  sc_thresh  = (int) ceil(((nullsc + (invP * eslCONST_LOG2) + 3.0) * om->scale_b)
                           + om->base_b + om->tec_b + om->tjb_b);
  sc_threshv = _mm256_set1_epi8((int8_t)(255 - sc_thresh));

  zerov    = _mm256_setzero_si256();
  biasv    = _mm256_set1_epi8((int8_t) om->bias_b);
  ceilingv = _mm256_cmpeq_epi8(biasv, biasv);
  for (q = 0; q < Q; q++) dp[q] = zerov;

  basev = _mm256_set1_epi8((int8_t) om->base_b);
  tjbmv = _mm256_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);
  xBv   = _mm256_subs_epu8(basev, tjbmv);

  for (i = 1; i <= L; i++) {
    rsc = om->rbv_avx[dsq[i]];
    xEv = zerov;

    mpv = esl_avx_rightshift_int8(dp[Q-1], zerov);
    for (q = 0; q < Q; q++) {
      sv   = _mm256_max_epu8(mpv, xBv);
      sv   = _mm256_adds_epu8(sv, biasv);
      sv   = _mm256_subs_epu8(sv, *rsc);  rsc++;
      xEv  = _mm256_max_epu8(xEv, sv);

      mpv   = dp[q];
      dp[q] = sv;
    }

    /* Test if any lane exceeded the p-value threshold */
    tempv = _mm256_adds_epu8(xEv, sc_threshv);
    tempv = _mm256_cmpeq_epi8(tempv, ceilingv);
    cmp   = _mm256_movemask_epi8(tempv);

    if (cmp != 0) {
      /* Find which model position hit threshold */
      end    = -1;
      rem_sc = -1;
      for (q = 0; q < Q; q++) {
        u.v = dp[q];
        for (k = 0; k < 32; k++) {
          if (u.b[k] >= sc_thresh && u.b[k] > rem_sc && (q + Q*k + 1) <= om->M) {
            end    = q + Q*k + 1;
            rem_sc = u.b[k];
          }
        }
        dp[q] = zerov; /* reset so next iter starts from xB */
      }

      /* Recover the diagonal that hit threshold */
      start       = end;
      target_end  = target_start = i;
      sc          = rem_sc;
      while (rem_sc > om->base_b - om->tjb_b - om->tbm_b) {
        rem_sc -= om->bias_b - ssvdata->ssv_scores[start * om->abc->Kp + dsq[target_start]];
        --start;
        --target_start;
      }
      start++;
      target_start++;

      /* Extend diagonal greedily */
      k            = end + 1;
      n            = target_end + 1;
      max_end      = target_end;
      max_sc       = sc;
      pos_since_max = 0;
      while (k < om->M && n <= L) {
        sc += om->bias_b - ssvdata->ssv_scores[k * om->abc->Kp + dsq[n]];
        if (sc >= max_sc) {
          max_sc       = sc;
          max_end      = n;
          pos_since_max = 0;
        } else {
          pos_since_max++;
          if (pos_since_max == 5) break;
        }
        k++;
        n++;
      }

      end        += (max_end - target_end);
      target_end  = max_end;

      ret_sc  = ((float)(max_sc - om->tjb_b) - (float) om->base_b);
      ret_sc /= om->scale_b;
      ret_sc -= 3.0f;

      p7_hmmwindow_new(windowlist,
                       0,
                       target_start,
                       end,
                       end - start + 1,
                       ret_sc,
                       p7_NOCOMPLEMENT,
                       L);

      i = target_end; /* skip forward */
    }
  }
  return eslOK;
}
/*------------------ end, p7_SSVFilter_BATH_avx() ------------------*/

#endif /* eslENABLE_AVX */
