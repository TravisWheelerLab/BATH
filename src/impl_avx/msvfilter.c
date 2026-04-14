/* MSV filter dispatcher.
 * Provides the non-suffixed extern functions declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"

/* Function:  p7_MSVFilter()
 *
 * Purpose:   Dispatch MSV filter to the fastest available ISA path.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> on score overflow (high-scoring hit).
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_MSVFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_MSVFilter_avx512(dsq, L, om, ox, ret_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_MSVFilter_sse(dsq, L, om, ox, ret_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_MSVFilter_sse(dsq, L, om, ox, ret_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_SSVFilter_BATH()
 *
 * Purpose:   Dispatch BATH SSV window finder to the fastest available ISA path.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_SSVFilter_BATH(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox,
                  const P7_SCOREDATA *ssvdata, P7_BG *bg, double P,
                  P7_HMM_WINDOWLIST *windowlist)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_SSVFilter_BATH_avx512(dsq, L, om, ox, ssvdata, bg, P, windowlist);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_SSVFilter_BATH_sse(dsq, L, om, ox, ssvdata, bg, P, windowlist);
#endif
#ifdef eslENABLE_SSE
  return p7_SSVFilter_BATH_sse(dsq, L, om, ox, ssvdata, bg, P, windowlist);
#else
  return eslENORESULT;
#endif
}
