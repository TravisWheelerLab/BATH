/* Viterbi filter dispatcher.
 * Provides the non-suffixed extern functions declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"

/* Function:  p7_ViterbiFilter()
 *
 * Purpose:   Dispatch Viterbi filter to the fastest available ISA path.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if score overflows; <*ret_sc> is <eslINFINITY>.
 * Throws:    <eslEINVAL> if <ox> is too small or profile is not local.
 */
int
p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_ViterbiFilter_avx512(dsq, L, om, ox, ret_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_ViterbiFilter_avx(dsq, L, om, ox, ret_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_ViterbiFilter_sse(dsq, L, om, ox, ret_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_ViterbiFilter_BATH()
 *
 * Purpose:   Dispatch BATH Viterbi filter (score + windows) to the fastest
 *            available ISA path.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if score overflows; <*ret_sc> is <eslINFINITY>.
 * Throws:    <eslEINVAL> if <ox> is too small or profile is not local.
 */
int
p7_ViterbiFilter_BATH(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox,
                      const P7_SCOREDATA *ssvdata, float filtersc, double P,
                      P7_HMM_WINDOWLIST *windowlist, float *ret_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_ViterbiFilter_BATH_avx512(dsq, L, om, ox, ssvdata, filtersc, P, windowlist, ret_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_ViterbiFilter_BATH_avx(dsq, L, om, ox, ssvdata, filtersc, P, windowlist, ret_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_ViterbiFilter_BATH_sse(dsq, L, om, ox, ssvdata, filtersc, P, windowlist, ret_sc);
#else
  return eslENORESULT;
#endif
}
