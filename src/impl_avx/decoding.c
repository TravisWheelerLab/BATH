/* Posterior decoding dispatcher.
 * Provides the non-suffixed p7_Decoding() and p7_DomainDecoding() declared
 * in impl_avx.h.  Delegates to the fastest available ISA implementation
 * at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Decoding()
 *
 * Purpose:   Dispatch posterior decoding of residue assignment to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if numeric range is exceeded.
 */
int
p7_Decoding(const P7_OPROFILE *om, const P7_OMX *oxf, P7_OMX *oxb, P7_OMX *pp)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Decoding_avx512(om, oxf, oxb, pp);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Decoding_sse(om, oxf, oxb, pp);
#endif
#ifdef eslENABLE_SSE
  return p7_Decoding_sse(om, oxf, oxb, pp);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_DomainDecoding()
 *
 * Purpose:   Dispatch posterior decoding of domain location to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> on numeric overflow.
 */
int
p7_DomainDecoding(const P7_OPROFILE *om, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_DomainDecoding_avx512(om, oxf, oxb, ddef);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_DomainDecoding_sse(om, oxf, oxb, ddef);
#endif
#ifdef eslENABLE_SSE
  return p7_DomainDecoding_sse(om, oxf, oxb, ddef);
#else
  return eslENORESULT;
#endif
}
