/* Spliced Viterbi algorithms dispatcher.
 * Provides the non-suffixed p7_Viterbi_Spliced() and p7_Viterbi_SplicedTrace()
 * declared in impl_avx.h.  Delegates to the fastest available ISA implementation.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Viterbi_Spliced()
 *
 * Purpose:   Dispatch spliced Viterbi DP fill to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on bad inputs.
 */
int
p7_Viterbi_Spliced(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_OMX *ox,
                   const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                   int i_start, int i_end, int min_intron, int global_start, int global_end)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Viterbi_Spliced_avx512(sub_dsq, om_tr, ox, signal_scores, acc_ov, don_ov, i_start, i_end, min_intron, global_start, global_end);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Viterbi_Spliced_avx(sub_dsq, om_tr, ox, signal_scores, acc_ov, don_ov, i_start, i_end, min_intron, global_start, global_end);
#endif
#ifdef eslENABLE_SSE
  return p7_Viterbi_Spliced_sse(sub_dsq, om_tr, ox, signal_scores, acc_ov, don_ov, i_start, i_end, min_intron, global_start, global_end);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_Viterbi_SplicedTrace()
 *
 * Purpose:   Dispatch spliced Viterbi traceback to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile is not 1-codon-length.
 *            <eslEFAIL>  if traceback fails.
 */
int
p7_Viterbi_SplicedTrace(const ESL_DSQ *sub_dsq, const P7_OMX *ox,
                        const P7_FS_PROFILE *gm_tr, const float *signal_scores,
                        P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                        int min_intron, float *vitsc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Viterbi_SplicedTrace_avx512(sub_dsq, ox, gm_tr, signal_scores, tr, i_start, i_end, k_start, k_end, min_intron, vitsc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Viterbi_SplicedTrace_avx(sub_dsq, ox, gm_tr, signal_scores, tr, i_start, i_end, k_start, k_end, min_intron, vitsc);
#endif
#ifdef eslENABLE_SSE
  return p7_Viterbi_SplicedTrace_sse(sub_dsq, ox, gm_tr, signal_scores, tr, i_start, i_end, k_start, k_end, min_intron, vitsc);
#else
  return eslENORESULT;
#endif
}
