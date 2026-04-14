/* Forward/Backward dispatcher.
 * Provides the non-suffixed p7_Forward, p7_ForwardParser, p7_Backward,
 * and p7_BackwardParser declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"

/* Function:  p7_Forward()
 *
 * Purpose:   Dispatch Forward algorithm to the fastest available ISA path.
 *
 * Returns:   <eslOK> on success.  <*opt_sc> is the log Forward score in
 *            nats.
 * Throws:    <eslEINVAL> if <ox> allocation is too small, or profile is
 *            not in a local alignment mode.
 *            <eslEMEM> on reallocation failure.
 */
int
p7_Forward(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *fwd, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Forward_avx512(dsq, L, om, fwd, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Forward_sse(dsq, L, om, fwd, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_Forward_sse(dsq, L, om, fwd, opt_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_ForwardParser()
 *
 * Purpose:   Dispatch Forward parser (linear memory) to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
int
p7_ForwardParser(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *fwd, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_ForwardParser_avx512(dsq, L, om, fwd, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_ForwardParser_sse(dsq, L, om, fwd, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_ForwardParser_sse(dsq, L, om, fwd, opt_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_Backward()
 *
 * Purpose:   Dispatch Backward algorithm to the fastest available ISA path.
 *
 * Returns:   <eslOK> on success.  <*opt_sc> is the log Backward score.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
int
p7_Backward(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Backward_avx512(dsq, L, om, fwd, bck, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Backward_sse(dsq, L, om, fwd, bck, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_Backward_sse(dsq, L, om, fwd, bck, opt_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_BackwardParser()
 *
 * Purpose:   Dispatch Backward parser (linear memory) to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
int
p7_BackwardParser(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_BackwardParser_avx512(dsq, L, om, fwd, bck, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_BackwardParser_sse(dsq, L, om, fwd, bck, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_BackwardParser_sse(dsq, L, om, fwd, bck, opt_sc);
#else
  return eslENORESULT;
#endif
}
