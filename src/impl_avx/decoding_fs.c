/* Frameshift posterior decoding dispatcher.
 * Provides the non-suffixed p7_Decoding_Frameshift() and
 * p7_DomainDecoding_Frameshift() declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Decoding_Frameshift()
 *
 * Purpose:   Dispatch frameshift posterior decoding (in-place overwrite of fwd)
 *            to the fastest available ISA path.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if numeric range is exceeded.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_Decoding_Frameshift(const P7_FS_OPROFILE *om_fs, P7_OMX *fwd, const P7_OMX *bck)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Decoding_Frameshift_avx512(om_fs, fwd, bck);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Decoding_Frameshift_avx(om_fs, fwd, bck);
#endif
#ifdef eslENABLE_SSE
  return p7_Decoding_Frameshift_sse(om_fs, fwd, bck);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_DomainDecoding_Frameshift()
 *
 * Purpose:   Dispatch frameshift domain decoding to the fastest available ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_DomainDecoding_Frameshift(const P7_FS_OPROFILE *om_fs, const P7_OMX *oxf, const P7_OMX *oxb,
                              P7_DOMAINDEF *ddef)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_DomainDecoding_Frameshift_avx512(om_fs, oxf, oxb, ddef);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_DomainDecoding_Frameshift_avx(om_fs, oxf, oxb, ddef);
#endif
#ifdef eslENABLE_SSE
  return p7_DomainDecoding_Frameshift_sse(om_fs, oxf, oxb, ddef);
#else
  return eslENORESULT;
#endif
}
