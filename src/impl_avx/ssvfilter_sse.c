/* SSV filter SSE implementation.
 * Ported from impl_sse/ssvfilter.c with _sse suffix for runtime dispatch.
 * See the original file for a detailed explanation of the multi-band macro scheme.
 *
 * Contents:
 *   1. Band calculation helper functions (static, via macros)
 *   2. p7_SSVFilter_sse()
 */

#include "p7_config.h"

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

/* Note that some ifdefs below have to be changed if these values are changed.
   Values chosen based on speed tests: two registers are generally used for
   other purposes, leaving 14 on 64-bit and 6 on 32-bit. */
#ifdef __x86_64__
#define MAX_BANDS 14
#else
#define MAX_BANDS 6
#endif


#define STEP_SINGLE(sv)                         \
  sv   = _mm_subs_epi8(sv, *rsc); rsc++;        \
  xEv  = _mm_max_epu8(xEv, sv);


#define LENGTH_CHECK(label)                     \
  if (i >= L) goto label;


#define NO_CHECK(label)


#define STEP_BANDS_1()                          \
  STEP_SINGLE(sv00)

#define STEP_BANDS_2()                          \
  STEP_BANDS_1()                                \
  STEP_SINGLE(sv01)

#define STEP_BANDS_3()                          \
  STEP_BANDS_2()                                \
  STEP_SINGLE(sv02)

#define STEP_BANDS_4()                          \
  STEP_BANDS_3()                                \
  STEP_SINGLE(sv03)

#define STEP_BANDS_5()                          \
  STEP_BANDS_4()                                \
  STEP_SINGLE(sv04)

#define STEP_BANDS_6()                          \
  STEP_BANDS_5()                                \
  STEP_SINGLE(sv05)

#define STEP_BANDS_7()                          \
  STEP_BANDS_6()                                \
  STEP_SINGLE(sv06)

#define STEP_BANDS_8()                          \
  STEP_BANDS_7()                                \
  STEP_SINGLE(sv07)

#define STEP_BANDS_9()                          \
  STEP_BANDS_8()                                \
  STEP_SINGLE(sv08)

#define STEP_BANDS_10()                         \
  STEP_BANDS_9()                                \
  STEP_SINGLE(sv09)

#define STEP_BANDS_11()                         \
  STEP_BANDS_10()                               \
  STEP_SINGLE(sv10)

#define STEP_BANDS_12()                         \
  STEP_BANDS_11()                               \
  STEP_SINGLE(sv11)

#define STEP_BANDS_13()                         \
  STEP_BANDS_12()                               \
  STEP_SINGLE(sv12)

#define STEP_BANDS_14()                         \
  STEP_BANDS_13()                               \
  STEP_SINGLE(sv13)

#define STEP_BANDS_15()                         \
  STEP_BANDS_14()                               \
  STEP_SINGLE(sv14)

#define STEP_BANDS_16()                         \
  STEP_BANDS_15()                               \
  STEP_SINGLE(sv15)

#define STEP_BANDS_17()                         \
  STEP_BANDS_16()                               \
  STEP_SINGLE(sv16)

#define STEP_BANDS_18()                         \
  STEP_BANDS_17()                               \
  STEP_SINGLE(sv17)


#define CONVERT_STEP(step, length_check, label, sv, pos)        \
  length_check(label)                                           \
  rsc = om->sbv[dsq[i]] + pos;                                 \
  step()                                                        \
  sv = _mm_slli_si128(sv, 1);                                   \
  sv = _mm_or_si128(sv, beginv);                                \
  i++;


#define CONVERT_1(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv00, Q - 1)

#define CONVERT_2(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv01, Q - 2)  \
  CONVERT_1(step, LENGTH_CHECK, label)

#define CONVERT_3(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv02, Q - 3)  \
  CONVERT_2(step, LENGTH_CHECK, label)

#define CONVERT_4(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv03, Q - 4)  \
  CONVERT_3(step, LENGTH_CHECK, label)

#define CONVERT_5(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv04, Q - 5)  \
  CONVERT_4(step, LENGTH_CHECK, label)

#define CONVERT_6(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv05, Q - 6)  \
  CONVERT_5(step, LENGTH_CHECK, label)

#define CONVERT_7(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv06, Q - 7)  \
  CONVERT_6(step, LENGTH_CHECK, label)

#define CONVERT_8(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv07, Q - 8)  \
  CONVERT_7(step, LENGTH_CHECK, label)

#define CONVERT_9(step, LENGTH_CHECK, label)            \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv08, Q - 9)  \
  CONVERT_8(step, LENGTH_CHECK, label)

#define CONVERT_10(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv09, Q - 10) \
  CONVERT_9(step, LENGTH_CHECK, label)

#define CONVERT_11(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv10, Q - 11) \
  CONVERT_10(step, LENGTH_CHECK, label)

#define CONVERT_12(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv11, Q - 12) \
  CONVERT_11(step, LENGTH_CHECK, label)

#define CONVERT_13(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv12, Q - 13) \
  CONVERT_12(step, LENGTH_CHECK, label)

#define CONVERT_14(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv13, Q - 14) \
  CONVERT_13(step, LENGTH_CHECK, label)

#define CONVERT_15(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv14, Q - 15) \
  CONVERT_14(step, LENGTH_CHECK, label)

#define CONVERT_16(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv15, Q - 16) \
  CONVERT_15(step, LENGTH_CHECK, label)

#define CONVERT_17(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv16, Q - 17) \
  CONVERT_16(step, LENGTH_CHECK, label)

#define CONVERT_18(step, LENGTH_CHECK, label)           \
  CONVERT_STEP(step, LENGTH_CHECK, label, sv17, Q - 18) \
  CONVERT_17(step, LENGTH_CHECK, label)


#define RESET_1()                               \
  register __m128i sv00 = beginv;

#define RESET_2()                               \
  RESET_1()                                     \
  register __m128i sv01 = beginv;

#define RESET_3()                               \
  RESET_2()                                     \
  register __m128i sv02 = beginv;

#define RESET_4()                               \
  RESET_3()                                     \
  register __m128i sv03 = beginv;

#define RESET_5()                               \
  RESET_4()                                     \
  register __m128i sv04 = beginv;

#define RESET_6()                               \
  RESET_5()                                     \
  register __m128i sv05 = beginv;

#define RESET_7()                               \
  RESET_6()                                     \
  register __m128i sv06 = beginv;

#define RESET_8()                               \
  RESET_7()                                     \
  register __m128i sv07 = beginv;

#define RESET_9()                               \
  RESET_8()                                     \
  register __m128i sv08 = beginv;

#define RESET_10()                              \
  RESET_9()                                     \
  register __m128i sv09 = beginv;

#define RESET_11()                              \
  RESET_10()                                    \
  register __m128i sv10 = beginv;

#define RESET_12()                              \
  RESET_11()                                    \
  register __m128i sv11 = beginv;

#define RESET_13()                              \
  RESET_12()                                    \
  register __m128i sv12 = beginv;

#define RESET_14()                              \
  RESET_13()                                    \
  register __m128i sv13 = beginv;

#define RESET_15()                              \
  RESET_14()                                    \
  register __m128i sv14 = beginv;

#define RESET_16()                              \
  RESET_15()                                    \
  register __m128i sv15 = beginv;

#define RESET_17()                              \
  RESET_16()                                    \
  register __m128i sv16 = beginv;

#define RESET_18()                              \
  RESET_17()                                    \
  register __m128i sv17 = beginv;


#define CALC(reset, step, convert, width)       \
  int i;                                        \
  int i2;                                       \
  int Q        = p7O_NQB(om->M);                \
  __m128i *rsc;                                 \
                                                \
  int w = width;                                \
                                                \
  dsq++;                                        \
                                                \
  reset()                                       \
                                                \
  for (i = 0; i < L && i < Q - q - w; i++)      \
    {                                           \
      rsc = om->sbv[dsq[i]] + i + q;            \
      step()                                    \
    }                                           \
                                                \
  i = Q - q - w;                                \
  convert(step, LENGTH_CHECK, done1)            \
done1:                                          \
                                                \
 for (i2 = Q - q; i2 < L - Q; i2 += Q)          \
   {                                            \
     for (i = 0; i < Q - w; i++)                \
       {                                        \
         rsc = om->sbv[dsq[i2 + i]] + i;        \
         step()                                 \
       }                                        \
                                                \
     i += i2;                                   \
     convert(step, NO_CHECK, )                  \
   }                                            \
                                                \
 for (i = 0; i2 + i < L && i < Q - w; i++)      \
   {                                            \
     rsc = om->sbv[dsq[i2 + i]] + i;            \
     step()                                     \
   }                                            \
                                                \
 i+=i2;                                         \
 convert(step, LENGTH_CHECK, done2)             \
done2:                                          \
                                                \
 return xEv;


/*****************************************************************
 * 1. Band calculation helper functions (static via macros).
 *****************************************************************/

static __m128i
calc_band_1(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_1, STEP_BANDS_1, CONVERT_1, 1)
}

static __m128i
calc_band_2(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_2, STEP_BANDS_2, CONVERT_2, 2)
}

static __m128i
calc_band_3(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_3, STEP_BANDS_3, CONVERT_3, 3)
}

static __m128i
calc_band_4(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_4, STEP_BANDS_4, CONVERT_4, 4)
}

static __m128i
calc_band_5(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_5, STEP_BANDS_5, CONVERT_5, 5)
}

static __m128i
calc_band_6(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_6, STEP_BANDS_6, CONVERT_6, 6)
}

#if MAX_BANDS > 6
static __m128i
calc_band_7(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_7, STEP_BANDS_7, CONVERT_7, 7)
}

static __m128i
calc_band_8(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_8, STEP_BANDS_8, CONVERT_8, 8)
}

static __m128i
calc_band_9(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_9, STEP_BANDS_9, CONVERT_9, 9)
}

static __m128i
calc_band_10(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_10, STEP_BANDS_10, CONVERT_10, 10)
}

static __m128i
calc_band_11(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_11, STEP_BANDS_11, CONVERT_11, 11)
}

static __m128i
calc_band_12(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_12, STEP_BANDS_12, CONVERT_12, 12)
}

static __m128i
calc_band_13(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_13, STEP_BANDS_13, CONVERT_13, 13)
}

static __m128i
calc_band_14(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_14, STEP_BANDS_14, CONVERT_14, 14)
}
#endif /* MAX_BANDS > 6 */

#if MAX_BANDS > 14
static __m128i
calc_band_15(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_15, STEP_BANDS_15, CONVERT_15, 15)
}

static __m128i
calc_band_16(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_16, STEP_BANDS_16, CONVERT_16, 16)
}

static __m128i
calc_band_17(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_17, STEP_BANDS_17, CONVERT_17, 17)
}

static __m128i
calc_band_18(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, int q, __m128i beginv, register __m128i xEv)
{
  CALC(RESET_18, STEP_BANDS_18, CONVERT_18, 18)
}
#endif /* MAX_BANDS > 14 */


/*****************************************************************
 * 2. p7_SSVFilter_sse()
 *****************************************************************/

static uint8_t
get_xE_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om)
{
  __m128i xEv;
  __m128i beginv;

  int q;
  int Q     = p7O_NQB(om->M);
  int bands;
  int last_q = 0;
  int i;

  __m128i (*fs[MAX_BANDS + 1])(const ESL_DSQ *, int, const P7_OPROFILE *, int, register __m128i, __m128i)
    = { NULL
      , calc_band_1,  calc_band_2,  calc_band_3,  calc_band_4,  calc_band_5,  calc_band_6
#if MAX_BANDS > 6
      , calc_band_7,  calc_band_8,  calc_band_9,  calc_band_10, calc_band_11, calc_band_12, calc_band_13, calc_band_14
#endif
#if MAX_BANDS > 14
      , calc_band_15, calc_band_16, calc_band_17, calc_band_18
#endif
    };

  beginv = _mm_set1_epi8(-128);
  xEv    = beginv;

  bands = (Q + MAX_BANDS - 1) / MAX_BANDS;

  for (i = 0; i < bands; i++) {
    q = (Q * (i + 1)) / bands;
    xEv = fs[q - last_q](dsq, L, om, last_q, beginv, xEv);
    last_q = q;
  }

  return esl_sse_hmax_epu8(xEv);
}


/* Function:  p7_SSVFilter_sse()
 * Synopsis:  SSV filter, SSE path.
 *
 * Purpose:   Calculate an approximation of the highest scoring single
 *            ungapped diagonal (SSV score) for digital sequence <dsq>
 *            of length <L> against profile <om>.  Returns the estimated
 *            score (nats) in <*ret_sc>.
 *
 *            Returns <eslENORESULT> when J-state use cannot be ruled out
 *            or the overflow analysis is inconclusive, signalling the
 *            caller to fall back to the full MSV filter.
 *
 *            Returns <eslERANGE> when the score provably overflows
 *            (high-scoring hit).
 *
 * Returns:   <eslOK> on success.
 *            <eslENORESULT> when the score may be underestimated.
 *            <eslERANGE> on confirmed overflow.
 */
int
p7_SSVFilter_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
  uint16_t xE;
  uint16_t xJ;

  if (om->tjb_b + om->tbm_b + om->tec_b + om->bias_b >= 127)
    return eslENORESULT;

  xE = get_xE_sse(dsq, L, om);

  if (xE >= 255 - om->bias_b) {
    *ret_sc = eslINFINITY;

    if (om->base_b - om->tjb_b - om->tbm_b < 128)
      return eslENORESULT;

    return eslERANGE;
  }

  xE += om->base_b - om->tjb_b - om->tbm_b;
  xE -= 128;

  if (xE >= 255 - om->bias_b) {
    *ret_sc = eslINFINITY;
    return eslERANGE;
  }

  xJ = xE - om->tec_b;

  if (xJ > om->base_b) return eslENORESULT;

  *ret_sc  = ((float)(xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0f;

  return eslOK;
}

#endif /* eslENABLE_SSE */
