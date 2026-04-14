/* Viterbi score (one-row, no traceback) dispatcher.
 * Provides the non-suffixed p7_ViterbiScore() declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_ViterbiScore()
 *
 * Purpose:   Dispatch one-row Viterbi score to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_ViterbiScore(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_ViterbiScore_avx512(dsq, L, om, ox, ret_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_ViterbiScore_sse(dsq, L, om, ox, ret_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_ViterbiScore_sse(dsq, L, om, ox, ret_sc);
#else
  return eslENORESULT;
#endif
}
