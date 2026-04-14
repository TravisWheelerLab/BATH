/* Frameshift optimal accuracy dispatcher.
 * Provides the non-suffixed p7_OptimalAccuracy_Frameshift() and
 * p7_OATrace_Frameshift() declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_OptimalAccuracy_Frameshift()
 *
 * Purpose:   Dispatch frameshift optimal accuracy DP fill to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_OptimalAccuracy_Frameshift(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                               P7_OMX *ox, float *ret_e)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_OptimalAccuracy_Frameshift_avx512(om_fs, pp, ox, ret_e);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_OptimalAccuracy_Frameshift_avx(om_fs, pp, ox, ret_e);
#endif
#ifdef eslENABLE_SSE
  return p7_OptimalAccuracy_Frameshift_sse(om_fs, pp, ox, ret_e);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_OATrace_Frameshift()
 *
 * Purpose:   Dispatch frameshift optimal accuracy traceback to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if trace is not empty.
 *            <eslEMEM> on allocation error.
 */
int
p7_OATrace_Frameshift(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                      const P7_OMX *ox, P7_TRACE *tr)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_OATrace_Frameshift_avx512(om_fs, pp, ox, tr);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_OATrace_Frameshift_avx(om_fs, pp, ox, tr);
#endif
#ifdef eslENABLE_SSE
  return p7_OATrace_Frameshift_sse(om_fs, pp, ox, tr);
#else
  return eslENORESULT;
#endif
}
