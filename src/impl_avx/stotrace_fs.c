/* Frameshift stochastic traceback dispatcher.
 * Provides the non-suffixed p7_StochasticTrace_Frameshift() declared in
 * impl_avx.h.  Delegates to the fastest available ISA implementation at
 * runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_StochasticTrace_Frameshift()
 *
 * Purpose:   Dispatch frameshift stochastic traceback to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on various internal failures.
 */
int
p7_StochasticTrace_Frameshift(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L,
                               const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, P7_TRACE *tr)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_StochasticTrace_Frameshift_avx512(rng, dsq, L, om_fs, ox, tr);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_StochasticTrace_Frameshift_avx(rng, dsq, L, om_fs, ox, tr);
#endif
#ifdef eslENABLE_SSE
  return p7_StochasticTrace_Frameshift_sse(rng, dsq, L, om_fs, ox, tr);
#else
  return eslENORESULT;
#endif
}
