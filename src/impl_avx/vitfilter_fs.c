/* Frameshift-aware Viterbi (full matrix) dispatcher.
 * Provides the non-suffixed p7_Viterbi_Frameshift() and
 * p7_Viterbi_Frameshift_Trace() declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Viterbi_Frameshift()
 *
 * Purpose:   Dispatch frameshift-aware full-matrix Viterbi to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on bad inputs; <eslERANGE> if no valid path exists.
 */
int
p7_Viterbi_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                      P7_OMX *ox, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Viterbi_Frameshift_avx512(dsq, L, om_fs, ox, ov, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Viterbi_Frameshift_sse(dsq, L, om_fs, ox, ov, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_Viterbi_Frameshift_sse(dsq, L, om_fs, ox, ov, opt_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_Viterbi_Frameshift_Trace()
 *
 * Purpose:   Dispatch frameshift Viterbi traceback to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslFAIL> if an impossible state is reached during traceback.
 */
int
p7_Viterbi_Frameshift_Trace(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                             const P7_OMX *ox, P7_TRACE *tr)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Viterbi_Frameshift_Trace_avx512(dsq, L, om_fs, ox, tr);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Viterbi_Frameshift_Trace_sse(dsq, L, om_fs, ox, tr);
#endif
#ifdef eslENABLE_SSE
  return p7_Viterbi_Frameshift_Trace_sse(dsq, L, om_fs, ox, tr);
#else
  return eslENORESULT;
#endif
}
