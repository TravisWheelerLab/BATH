/* Frameshift Forward/Backward dispatcher.
 * Provides the non-suffixed frameshift Forward/Backward functions declared
 * in impl_avx.h.  Delegates to the fastest available ISA implementation
 * at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_ForwardParser_Frameshift_3Codons()
 *
 * Purpose:   Dispatch 3-codon frameshift Forward parser to the fastest ISA.
 *
 * Returns:   <eslOK> on success.  <*opt_sc> is the log Forward score in nats.
 * Throws:    <eslEINVAL> if allocation is too small.
 *            <eslEMEM> on reallocation failure.
 */
int
p7_ForwardParser_Frameshift_3Codons(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                    P7_OMX *ox, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_ForwardParser_Frameshift_3Codons_avx512(dsq, L, om_fs, ox, ov, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_ForwardParser_Frameshift_3Codons_avx(dsq, L, om_fs, ox, ov, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_ForwardParser_Frameshift_3Codons_sse(dsq, L, om_fs, ox, ov, opt_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_BackwardParser_Frameshift_3Codons()
 *
 * Purpose:   Dispatch 3-codon frameshift Backward parser to the fastest ISA.
 *
 * Returns:   <eslOK> on success.  <*opt_sc> is the log Backward score.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
int
p7_BackwardParser_Frameshift_3Codons(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                     const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_BackwardParser_Frameshift_3Codons_avx512(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_BackwardParser_Frameshift_3Codons_avx(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_BackwardParser_Frameshift_3Codons_sse(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_ForwardParser_Frameshift_5Codons()
 *
 * Purpose:   Dispatch 5-codon frameshift Forward parser to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
int
p7_ForwardParser_Frameshift_5Codons(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                    P7_OMX *ox, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_ForwardParser_Frameshift_5Codons_avx512(dsq, L, om_fs, ox, ov, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_ForwardParser_Frameshift_5Codons_avx(dsq, L, om_fs, ox, ov, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_ForwardParser_Frameshift_5Codons_sse(dsq, L, om_fs, ox, ov, opt_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_BackwardParser_Frameshift_5Codons()
 *
 * Purpose:   Dispatch 5-codon frameshift Backward parser to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
int
p7_BackwardParser_Frameshift_5Codons(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                     const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_BackwardParser_Frameshift_5Codons_avx512(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_BackwardParser_Frameshift_5Codons_avx(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_BackwardParser_Frameshift_5Codons_sse(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_Forward_Frameshift()
 *
 * Purpose:   Dispatch full-matrix 5-codon frameshift Forward to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
int
p7_Forward_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                      P7_OMX *ox, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Forward_Frameshift_avx512(dsq, L, om_fs, ox, ov, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Forward_Frameshift_avx(dsq, L, om_fs, ox, ov, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_Forward_Frameshift_sse(dsq, L, om_fs, ox, ov, opt_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_Backward_Frameshift()
 *
 * Purpose:   Dispatch full-matrix 5-codon frameshift Backward to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
int
p7_Backward_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                       const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Backward_Frameshift_avx512(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Backward_Frameshift_avx(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_Backward_Frameshift_sse(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#else
  return eslENORESULT;
#endif
}
