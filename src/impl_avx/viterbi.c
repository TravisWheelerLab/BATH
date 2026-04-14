/* Viterbi full-matrix dispatcher.
 * Provides the non-suffixed p7_Viterbi() and p7_Viterbi_Trace() declared in
 * impl_avx.h.  Delegates to the fastest available ISA implementation.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Viterbi()
 *
 * Purpose:   Dispatch full-matrix Viterbi to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_Viterbi(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Viterbi_avx512(dsq, L, om, ox, ret_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Viterbi_avx(dsq, L, om, ox, ret_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_Viterbi_sse(dsq, L, om, ox, ret_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_Viterbi_Trace()
 *
 * Purpose:   Dispatch Viterbi traceback to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if an impossible state is reached.
 */
int
p7_Viterbi_Trace(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *ox,
                 P7_TRACE *tr)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Viterbi_Trace_avx512(dsq, L, om, ox, tr);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Viterbi_Trace_avx(dsq, L, om, ox, tr);
#endif
#ifdef eslENABLE_SSE
  return p7_Viterbi_Trace_sse(dsq, L, om, ox, tr);
#else
  return eslENORESULT;
#endif
}
