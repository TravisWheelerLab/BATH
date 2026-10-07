/* SSV filter dispatcher.
 * Provides the non-suffixed p7_SSVFilter() declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"

/* Forward declaration of dispatcher */
static int p7_SSVFilter_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_SSVFilter)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc) = p7_SSVFilter_Dispatcher;

/* Function:  p7_SSVFilter_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for SSV filter.
 *
 * Returns:   <eslOK> on success.
 *            <eslENORESULT> when J-state use cannot be ruled out.
 *            <eslERANGE> on confirmed overflow (high-scoring hit).
 */
static int
p7_SSVFilter_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_SSVFilter = p7_SSVFilter_avx512; return p7_SSVFilter_avx512(dsq, L, om, ret_sc); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_SSVFilter = p7_SSVFilter_avx;    return p7_SSVFilter_avx(dsq, L, om, ret_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_SSVFilter = p7_SSVFilter_sse;
  return p7_SSVFilter_sse(dsq, L, om, ret_sc);
#else
  p7_Die("p7_SSVFilter: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Function:  p7_SSVFilter_FromXE()
 * Synopsis:  The SSV filter's result, from the SSV matrix's highest cell.
 *
 * Purpose:   <xE> is the highest cell of the SSV matrix of a sequence against
 *            <om>, as p7_SSVFilter_OrfBlock() gives it. Return what
 *            p7_SSVFilter() returns for that sequence, and its score in
 *            <*ret_sc>; <om> must be configured for the sequence's length.
 *
 * Returns:   <eslOK>, <eslENORESULT> or <eslERANGE>, as p7_SSVFilter().
 */
int
p7_SSVFilter_FromXE(int xE_in, const P7_OPROFILE *om, float *ret_sc)
{
  uint16_t xE = (uint16_t) xE_in;
  uint16_t xJ;

  if (om->tjb_b + om->tbm_b + om->tec_b + om->bias_b >= 127)
    return eslENORESULT;

  /* Saturation floors every diagonal at the begin score (128), so a
   * max of 128 means no diagonal scored above it and the true best
   * may be lower; let the full MSV filter compute it. */
  if (xE <= 128) return eslENORESULT;

  if (xE >= 255 - om->bias_b) {
    *ret_sc = eslINFINITY;

    if (om->base_b - om->tjb_b - om->tbm_b < 128)
      return eslENORESULT;

    return eslERANGE;
  }

  xE += om->base_b - om->tjb_b - om->tbm_b;
  xE -= 128;

  if (xE >= 255 - om->bias_b) {
    *ret_sc = eslINFINITY;
    return eslERANGE;
  }

  xJ = xE - om->tec_b;

  if (xJ > om->base_b) return eslENORESULT;

  *ret_sc  = ((float)(xJ - om->tjb_b) - (float) om->base_b);
  *ret_sc /= om->scale_b;
  *ret_sc -= 3.0f;

  return eslOK;
}


/* Function:  p7_SSVFilter_OrfBlock()
 * Synopsis:  SSV maximum of every ORF in a block.
 *
 * Purpose:   For each of the <n> ORFs in <orf>, store in <xE[i]> the highest
 *            cell of its SSV matrix against <om>. It does not depend on the
 *            length <om> is configured for; p7_SSVFilter_FromXE() makes the
 *            filter's status and score from it.
 *
 * Returns:   <eslOK> on success; <eslENORESULT> if there is no SIMD
 *            implementation, and <xE> is unset.
 */
int
p7_SSVFilter_OrfBlock(const P7_OPROFILE *om, const ESL_ORF *orf, int n, uint8_t *xE)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_SSVFilter_OrfBlock_avx(om, orf, n, xE); return eslOK; }
#endif
#ifdef eslENABLE_SSE
  p7_SSVFilter_OrfBlock_sse(om, orf, n, xE);
  return eslOK;
#else
  return eslENORESULT;
#endif
}
