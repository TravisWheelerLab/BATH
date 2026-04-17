/* SSV filter dispatcher.
 * Provides the non-suffixed p7_SSVFilter() declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"

/* Forward declaration of dispatcher */
static int p7_SSVFilter_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_SSVFilter)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc) = p7_SSVFilter_Dispatcher;

/* Function:  p7_SSVFilter_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for SSV filter.
 *
 * Returns:   <eslOK> on success.
 *            <eslENORESULT> when J-state use cannot be ruled out.
 *            <eslERANGE> on confirmed overflow (high-scoring hit).
 */
static int
p7_SSVFilter_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_SSVFilter = p7_SSVFilter_avx512; return p7_SSVFilter_avx512(dsq, L, om, ret_sc); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_SSVFilter = p7_SSVFilter_avx;    return p7_SSVFilter_avx(dsq, L, om, ret_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_SSVFilter = p7_SSVFilter_sse;
  return p7_SSVFilter_sse(dsq, L, om, ret_sc);
#else
  p7_Die("p7_SSVFilter: no SIMD implementation available");
  return eslENORESULT;
#endif
}
