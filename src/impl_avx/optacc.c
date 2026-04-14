/* Optimal accuracy alignment dispatcher.
 * Provides the non-suffixed p7_OptimalAccuracy() and p7_OATrace() declared
 * in impl_avx.h.  Delegates to the fastest available ISA implementation
 * at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_OptimalAccuracy()
 *
 * Purpose:   Dispatch optimal accuracy DP fill to the fastest available ISA.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_OptimalAccuracy(const P7_OPROFILE *om, const P7_OMX *pp, P7_OMX *ox, float *ret_e)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_OptimalAccuracy_avx512(om, pp, ox, ret_e);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_OptimalAccuracy_avx(om, pp, ox, ret_e);
#endif
#ifdef eslENABLE_SSE
  return p7_OptimalAccuracy_sse(om, pp, ox, ret_e);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_OATrace()
 *
 * Purpose:   Dispatch optimal accuracy traceback to the fastest available ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> if trace is not empty.
 */
int
p7_OATrace(const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, P7_TRACE *tr)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_OATrace_avx512(om, pp, ox, tr);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_OATrace_avx(om, pp, ox, tr);
#endif
#ifdef eslENABLE_SSE
  return p7_OATrace_sse(om, pp, ox, tr);
#else
  return eslENORESULT;
#endif
}
