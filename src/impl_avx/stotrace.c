/* Stochastic traceback dispatcher.
 * Provides the non-suffixed p7_StochasticTrace() declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_StochasticTrace()
 *
 * Purpose:   Dispatch stochastic traceback to the fastest available ISA path.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> on various internal failures.
 */
int
p7_StochasticTrace(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,
                   const P7_OMX *ox, P7_TRACE *tr)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_StochasticTrace_avx512(rng, dsq, L, om, ox, tr);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_StochasticTrace_sse(rng, dsq, L, om, ox, tr);
#endif
#ifdef eslENABLE_SSE
  return p7_StochasticTrace_sse(rng, dsq, L, om, ox, tr);
#else
  return eslENORESULT;
#endif
}
