/* P7_OIVX SSE implementations.
 * Vectorized intermediate-value matrix for the spliced Viterbi.
 * Ported from impl_sse/p7_oivx.c with _sse suffix for runtime dispatch.
 *
 * Contents:
 *    1. The P7_OIVX SSE lifecycle functions.
 */

#include "p7_config.h"

#include "easel.h"
#include "hmmer.h"
#include "impl_avx.h"

/*****************************************************************
 * 1. The P7_OIVX SSE lifecycle functions.
 *****************************************************************/

/* Function:  p7_oivx_Create_sse()
 *
 * Purpose:   Allocate a <P7_OIVX> for a model of up to <M_hint> nodes
 *            and <C> channels using SSE (128-bit) vectors.  Storage
 *            layout is [C][Q] where Q = p7O_NQF(M_hint).  AVX and
 *            AVX-512 pointer fields are initialised to NULL so the
 *            object can be safely destroyed by any ISA path.
 *
 * Returns:   Pointer to the new <P7_OIVX> on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_OIVX *
p7_oivx_Create_sse(int M_hint, int C)
{
#ifdef eslENABLE_SSE
  P7_OIVX *ov;
  int      nq = p7O_NQF(M_hint);
  int      c;
  int      status;

  ov = NULL;
  ESL_ALLOC(ov, sizeof(P7_OIVX));
  ov->ivx_mem = NULL;
  ov->ivx     = NULL;

  ov->allocM  = M_hint;
  ov->allocC  = C;
  ov->allocQ4 = nq;

#ifdef eslENABLE_AVX
  ov->ivx_mem_avx = NULL;
  ov->ivx_avx     = NULL;
  ov->allocQ4_avx = 0;
#endif
#ifdef eslENABLE_AVX512
  ov->ivx_mem_avx512 = NULL;
  ov->ivx_avx512     = NULL;
  ov->allocQ4_avx512 = 0;
#endif

  ESL_ALLOC(ov->ivx_mem, sizeof(__m128) * C * nq + 15);
  ESL_ALLOC(ov->ivx,     sizeof(__m128 *) * C);

  ov->ivx[0] = (__m128 *) (((unsigned long int) ov->ivx_mem + 15) & (~0xf));
  for (c = 1; c < C; c++)
    ov->ivx[c] = ov->ivx[0] + (c * nq);

  return ov;

 ERROR:
  p7_oivx_Destroy_sse(ov);
  return NULL;
#else
  return NULL;
#endif
}


/* Function:  p7_oivx_GrowTo_sse()
 *
 * Purpose:   Ensure that <ov> (an SSE-allocated <P7_OIVX>) has storage
 *            for at least <M> model nodes and <C> channels.  Returns
 *            immediately if the existing allocation is sufficient.
 *            Otherwise frees and re-allocates the SSE backing memory
 *            and re-establishes the pointer array.
 *
 *            Any data in <ov> must be considered invalid after a
 *            successful reallocation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; <ov> is destroyed.
 */
int
p7_oivx_GrowTo_sse(P7_OIVX *ov, int M, int C)
{
#ifdef eslENABLE_SSE
  int   nq = p7O_NQF(M);
  int   c;
  void *p;
  int   status;

  if (C * nq <= ov->allocC * ov->allocQ4 && C <= ov->allocC) return eslOK;

  if (ov->ivx_mem != NULL) { free(ov->ivx_mem); ov->ivx_mem = NULL; }

  if (C > ov->allocC) {
    ESL_RALLOC(ov->ivx, p, sizeof(__m128 *) * C);
  }

  ESL_ALLOC(ov->ivx_mem, sizeof(__m128) * C * nq + 15);

  ov->ivx[0] = (__m128 *) (((unsigned long int) ov->ivx_mem + 15) & (~0xf));
  for (c = 1; c < C; c++)
    ov->ivx[c] = ov->ivx[0] + (c * nq);

  ov->allocM  = M;
  ov->allocC  = C;
  ov->allocQ4 = nq;

  return eslOK;

 ERROR:
  p7_oivx_Destroy_sse(ov);
  return status;
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_oivx_Destroy_sse()
 *
 * Purpose:   Free a <P7_OIVX> that was allocated by p7_oivx_Create_sse().
 *            Safe to call on NULL.
 */
void
p7_oivx_Destroy_sse(P7_OIVX *ov)
{
  if (ov == NULL) return;

#ifdef eslENABLE_SSE
  if (ov->ivx_mem != NULL) free(ov->ivx_mem);
  if (ov->ivx     != NULL) free(ov->ivx);
#endif

  free(ov);
}
