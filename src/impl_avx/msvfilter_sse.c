/* MSV filter SSE implementation.
 * Ported from impl_sse/msvfilter.c with _sse suffix for runtime dispatch.
 *
 * Contents:
 *    1. p7_MSVFilter_sse()    — fast byte-precision MSV score
 *    2. p7_SSVFilter_BATH_sse() — BATH ungapped window finder
 */

#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#ifdef eslENABLE_SSE
#include <xmmintrin.h>
#include <emmintrin.h>
#endif

#include "easel.h"
#include "esl_gumbel.h"

#include "hmmer.h"
#include "impl_avx.h"

/*****************************************************************
 * 1. p7_MSVFilter_sse()
 *****************************************************************/

/* Function:  p7_MSVFilter_sse()
 * Synopsis:  MSV score in uint8 precision, SSE path.
 *
 * Purpose:   Calculate an approximation of the MSV score for digital
 *            sequence <dsq> of length <L>, using optimized profile <om>
 *            and pre-allocated one-row DP matrix <ox>.  Returns the
 *            estimated MSV score (nats) in <*ret_sc>.
 *
 *            Tries p7_SSVFilter_sse() first; falls back to the full MSV
 *            calculation only when SSV returns <eslENORESULT>.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if the score overflows uint8 range (high-scoring hit).
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_MSVFilter_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
#ifdef eslENABLE_SSE
  register __m128i mpv;
  register __m128i xEv;
  register __m128i xBv;
  register __m128i sv;
  register __m128i biasv;
  uint8_t  xJ;
  int i;
  int q;
  int Q        = p7O_NQB(om->M);
  __m128i *dp  = ox->dpb[0];
  __m128i *rsc;

  __m128i xJv;
  __m128i tjbmv;
  __m128i tecv;
  __m128i basev;
  __m128i ceilingv;
  __m128i tempv;

  int cmp;
  int status = eslOK;

  if (Q > ox->allocQ16) ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M = om->M;

  /* Try highly optimized SSV filter first */
  status = p7_SSVFilter_sse(dsq, L, om, ret_sc);
  if (status != eslENORESULT) return status;

  /* Initialization. In offset unsigned arithmetic, -infinity is 0,
   * and 0 is om->base_b. */
  biasv    = _mm_set1_epi8((int8_t) om->bias_b);
  for (q = 0; q < Q; q++) dp[q] = _mm_setzero_si128();
  xJ       = 0;

  ceilingv = _mm_cmpeq_epi8(biasv, biasv);
  basev    = _mm_set1_epi8((int8_t) om->base_b);
  tjbmv    = _mm_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);
  tecv     = _mm_set1_epi8((int8_t) om->tec_b);

  xJv = _mm_subs_epu8(biasv, biasv);
  xBv = _mm_subs_epu8(basev, tjbmv);

#if eslDEBUGLEVEL > 0
  if (ox->debugging) {
    uint8_t xB;
    xB = _mm_extract_epi16(xBv, 0);
    xJ = _mm_extract_epi16(xJv, 0);
    p7_omx_DumpMFRow(ox, 0, 0, 0, xJ, xB, xJ);
  }
#endif

  for (i = 1; i <= L; i++) {
    rsc = om->rbv[dsq[i]];
    xEv = _mm_setzero_si128();

    /* Right-shift by 1 byte (little-endian: left bit shift).
     * Zeros fill automatically — our -infinity. */
    mpv = _mm_slli_si128(dp[Q-1], 1);
    for (q = 0; q < Q; q++) {
      sv   = _mm_max_epu8(mpv, xBv);
      sv   = _mm_adds_epu8(sv, biasv);
      sv   = _mm_subs_epu8(sv, *rsc);  rsc++;
      xEv  = _mm_max_epu8(xEv, sv);

      mpv   = dp[q];
      dp[q] = sv;
    }

    /* Overflow test */
    tempv = _mm_adds_epu8(xEv, biasv);
    tempv = _mm_cmpeq_epi8(tempv, ceilingv);
    cmp   = _mm_movemask_epi8(tempv);

    /* Horizontal max via shuffles */
    tempv = _mm_shuffle_epi32(xEv, _MM_SHUFFLE(2, 3, 0, 1));
    xEv   = _mm_max_epu8(xEv, tempv);
    tempv = _mm_shuffle_epi32(xEv, _MM_SHUFFLE(0, 1, 2, 3));
    xEv   = _mm_max_epu8(xEv, tempv);
    tempv = _mm_shufflelo_epi16(xEv, _MM_SHUFFLE(2, 3, 0, 1));
    xEv   = _mm_max_epu8(xEv, tempv);
    tempv = _mm_srli_si128(xEv, 1);
    xEv   = _mm_max_epu8(xEv, tempv);
    xEv   = _mm_shuffle_epi32(xEv, _MM_SHUFFLE(0, 0, 0, 0));

    if (cmp != 0x0000) {
      *ret_sc = eslINFINITY;
      return eslERANGE;
    }

    xEv = _mm_subs_epu8(xEv, tecv);
    xJv = _mm_max_epu8(xJv, xEv);
    xBv = _mm_max_epu8(basev, xJv);
    xBv = _mm_subs_epu8(xBv, tjbmv);

#if eslDEBUGLEVEL > 0
    if (ox->debugging) {
      uint8_t xB, xE;
      xB = _mm_extract_epi16(xBv, 0);
      xE = _mm_extract_epi16(xEv, 0);
      xJ = _mm_extract_epi16(xJv, 0);
      p7_omx_DumpMFRow(ox, i, xE, 0, xJ, xB, xJ);
    }
#endif
  }

  xJ = (uint8_t) _mm_extract_epi16(xJv, 0);

  *ret_sc  = ((float)(xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0f;

  return eslOK;
#else
  return eslENORESULT;
#endif
}
/*------------------ end, p7_MSVFilter_sse() ----------------------*/


/*****************************************************************
 * 2. p7_SSVFilter_BATH_sse()
 *****************************************************************/

/* Function:  p7_SSVFilter_BATH_sse()
 * Synopsis:  BATH ungapped window finder, SSE path.
 *
 * Purpose:   SSV scan over <dsq[1..L]> using profile <om> and
 *            one-row DP matrix <ox>.  Appends to <windowlist> every
 *            ungapped diagonal whose score exceeds the threshold
 *            implied by <P>, <bg>, and the MSV EVD parameters.
 *
 *            Used by the frameshift (--fs) and splice (--splice) paths
 *            to seed candidate exon windows.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_SSVFilter_BATH_sse(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox,
                      const P7_SCOREDATA *ssvdata, P7_BG *bg, double P,
                      P7_HMM_WINDOWLIST *windowlist)
{
#ifdef eslENABLE_SSE
  register __m128i mpv;
  register __m128i xEv;
  register __m128i xBv;
  register __m128i sv;
  register __m128i biasv;
  int i;
  int q;
  int Q        = p7O_NQB(om->M);
  __m128i *dp  = ox->dpb[0];
  __m128i *rsc;
  __m128i tjbmv;
  __m128i basev;
  __m128i ceilingv;
  __m128i tempv;
  int cmp;
  int k, n, end, rem_sc, start, target_end, target_start;
  int max_end, max_sc, sc, pos_since_max;
  float nullsc;
  float ret_sc;

  union { __m128i v; uint8_t b[16]; } u;

  __m128i sc_threshv;
  uint8_t sc_thresh;
  float invP = esl_gumbel_invsurv(P, om->evparam[p7_MMU], om->evparam[p7_MLAMBDA]);

  if (Q > ox->allocQ16) ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M = om->M;

  p7_bg_SetLength(bg, L);
  p7_oprofile_ReconfigMSVLength(om, L);
  p7_bg_NullOne(bg, dsq, L, &nullsc);

  sc_thresh  = (int) ceil(((nullsc + (invP * eslCONST_LOG2) + 3.0) * om->scale_b)
                           + om->base_b + om->tec_b + om->tjb_b);
  sc_threshv = _mm_set1_epi8((int8_t)(255 - sc_thresh));

  biasv    = _mm_set1_epi8((int8_t) om->bias_b);
  ceilingv = _mm_cmpeq_epi8(biasv, biasv);
  for (q = 0; q < Q; q++) dp[q] = _mm_setzero_si128();

  basev = _mm_set1_epi8((int8_t) om->base_b);
  tjbmv = _mm_set1_epi8((int8_t) om->tjb_b + (int8_t) om->tbm_b);
  xBv   = _mm_subs_epu8(basev, tjbmv);

  for (i = 1; i <= L; i++) {
    rsc = om->rbv[dsq[i]];
    xEv = _mm_setzero_si128();

    mpv = _mm_slli_si128(dp[Q-1], 1);
    for (q = 0; q < Q; q++) {
      sv   = _mm_max_epu8(mpv, xBv);
      sv   = _mm_adds_epu8(sv, biasv);
      sv   = _mm_subs_epu8(sv, *rsc);  rsc++;
      xEv  = _mm_max_epu8(xEv, sv);

      mpv   = dp[q];
      dp[q] = sv;
    }

    /* Test if any lane exceeded the p-value threshold */
    tempv = _mm_adds_epu8(xEv, sc_threshv);
    tempv = _mm_cmpeq_epi8(tempv, ceilingv);
    cmp   = _mm_movemask_epi8(tempv);

    if (cmp != 0) {
      /* Find which model position hit threshold */
      end    = -1;
      rem_sc = -1;
      for (q = 0; q < Q; q++) {
        u.v = dp[q];
        for (k = 0; k < 16; k++) {
          if (u.b[k] >= sc_thresh && u.b[k] > rem_sc && (q + Q*k + 1) <= om->M) {
            end    = q + Q*k + 1;
            rem_sc = u.b[k];
          }
        }
        dp[q] = _mm_set1_epi8(0); /* reset so next iter starts from xB */
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
#else
  return eslENORESULT;
#endif
}
/*------------------ end, p7_SSVFilter_BATH_sse() ------------------*/
