/* Null2 biased composition correction dispatcher.
 * Provides the non-suffixed p7_Null2_ByExpectation() and p7_Null2_ByTrace()
 * declared in impl_avx.h.  Delegates to the fastest available ISA.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Null2_ByExpectation()
 *
 * Purpose:   Dispatch null2 estimation (expectation) to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_Null2_ByExpectation(const P7_OPROFILE *om, const P7_OMX *pp, float *null2)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Null2_ByExpectation_avx512(om, pp, null2);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Null2_ByExpectation_avx(om, pp, null2);
#endif
#ifdef eslENABLE_SSE
  return p7_Null2_ByExpectation_sse(om, pp, null2);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_Null2_ByTrace()
 *
 * Purpose:   Dispatch null2 estimation (trace) to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_Null2_ByTrace(const P7_OPROFILE *om, const P7_TRACE *tr, int zstart, int zend,
                 P7_OMX *wrk, float *null2)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Null2_ByTrace_avx512(om, tr, zstart, zend, wrk, null2);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Null2_ByTrace_avx(om, tr, zstart, zend, wrk, null2);
#endif
#ifdef eslENABLE_SSE
  return p7_Null2_ByTrace_sse(om, tr, zstart, zend, wrk, null2);
#else
  return eslENORESULT;
#endif
}
