/* Frameshift null2 dispatcher.
 * Provides the non-suffixed p7_Null2_fs_ByExpectation() declared in
 * impl_avx.h.  Delegates to the fastest available ISA implementation.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Null2_fs_ByExpectation()
 *
 * Purpose:   Dispatch frameshift null2 estimation (expectation) to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if codon_lengths is not 1, 3, or 5.
 */
int
p7_Null2_fs_ByExpectation(const P7_FS_OPROFILE *om_fs, P7_OMX *pp, float *null2)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Null2_fs_ByExpectation_avx512(om_fs, pp, null2);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Null2_fs_ByExpectation_avx(om_fs, pp, null2);
#endif
#ifdef eslENABLE_SSE
  return p7_Null2_fs_ByExpectation_sse(om_fs, pp, null2);
#else
  return eslENORESULT;
#endif
}
