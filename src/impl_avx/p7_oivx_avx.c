/* p7_oivx_avx.c — AVX2 (256-bit) implementations of P7_OIVX management.
 *
 * Ported from p7_oivx_sse.c with:
 *   - __m128 / p7O_NQF / allocQ4 / ivx / ivx_mem
 *       → __m256 / p7O_NQF_AVX / allocQ4_avx / ivx_avx / ivx_mem_avx
 *   - Memory alignment padding +15/~0xf → +31/~0x1f (32-byte AVX alignment)
 *   - sizeof(__m128) / sizeof(__m128 *) → sizeof(__m256) / sizeof(__m256 *)
 */
#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <immintrin.h>

#include "easel.h"
#include "hmmer.h"
#include "impl_avx.h"


/*****************************************************************
 * 1. The P7_OIVX AVX2 lifecycle functions.
 *****************************************************************/

/* Function:  p7_oivx_Create_avx()
 *
 * Purpose:   Allocate a <P7_OIVX> for a model of up to <M_hint> nodes
 *            and <C> channels using AVX2 (256-bit) vectors.  Storage
 *            layout is [C][Q] where Q = p7O_NQF_AVX(M_hint).  SSE and
 *            AVX-512 pointer fields are initialised to NULL so the
 *            object can be safely destroyed by any ISA path.
 *
 * Returns:   Pointer to the new <P7_OIVX> on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_OIVX *
p7_oivx_Create_avx(int M_hint, int C)
{
  P7_OIVX *ov;
  int      nq = p7O_NQF_AVX(M_hint);
  int      c;
  int      status;

  ov = NULL;
  ESL_ALLOC(ov, sizeof(P7_OIVX));

  /* SSE fields — not used on the AVX path */
  ov->ivx_mem = NULL;
  ov->ivx     = NULL;
  ov->allocQ4 = 0;

  /* AVX fields */
  ov->ivx_mem_avx = NULL;
  ov->ivx_avx     = NULL;
  ov->allocQ4_avx = 0;

#ifdef eslENABLE_AVX512
  ov->ivx_mem_avx512 = NULL;
  ov->ivx_avx512     = NULL;
  ov->allocQ4_avx512 = 0;
#endif

  ov->allocM = M_hint;
  ov->allocC = C;

  ESL_ALLOC(ov->ivx_mem_avx, sizeof(__m256) * C * nq + 31);
  ESL_ALLOC(ov->ivx_avx,     sizeof(__m256 *) * C);

  ov->ivx_avx[0] = (__m256 *) (((unsigned long int) ov->ivx_mem_avx + 31) & (~0x1f));
  for (c = 1; c < C; c++)
    ov->ivx_avx[c] = ov->ivx_avx[0] + (c * nq);

  ov->allocQ4_avx = nq;

  return ov;

 ERROR:
  p7_oivx_Destroy_avx(ov);
  return NULL;
}


/* Function:  p7_oivx_GrowTo_avx()
 *
 * Purpose:   Ensure that <ov> (an AVX2-allocated <P7_OIVX>) has storage
 *            for at least <M> model nodes and <C> channels.  Returns
 *            immediately if the existing allocation is sufficient.
 *            Otherwise frees and re-allocates the AVX backing memory
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
p7_oivx_GrowTo_avx(P7_OIVX *ov, int M, int C)
{
  int   nq = p7O_NQF_AVX(M);
  int   c;
  void *p;
  int   status;

  if (C * nq <= ov->allocC * ov->allocQ4_avx && C <= ov->allocC) return eslOK;

  if (ov->ivx_mem_avx != NULL) { free(ov->ivx_mem_avx); ov->ivx_mem_avx = NULL; }

  if (C > ov->allocC) {
    ESL_RALLOC(ov->ivx_avx, p, sizeof(__m256 *) * C);
  }

  ESL_ALLOC(ov->ivx_mem_avx, sizeof(__m256) * C * nq + 31);

  ov->ivx_avx[0] = (__m256 *) (((unsigned long int) ov->ivx_mem_avx + 31) & (~0x1f));
  for (c = 1; c < C; c++)
    ov->ivx_avx[c] = ov->ivx_avx[0] + (c * nq);

  ov->allocM      = M;
  ov->allocC      = C;
  ov->allocQ4_avx = nq;

  return eslOK;

 ERROR:
  p7_oivx_Destroy_avx(ov);
  return status;
}


/* Function:  p7_oivx_Destroy_avx()
 *
 * Purpose:   Free a <P7_OIVX> that was allocated by p7_oivx_Create_avx().
 *            Safe to call on NULL.
 */
void
p7_oivx_Destroy_avx(P7_OIVX *ov)
{
  if (ov == NULL) return;

  if (ov->ivx_mem_avx != NULL) free(ov->ivx_mem_avx);
  if (ov->ivx_avx     != NULL) free(ov->ivx_avx);

  free(ov);
}

#endif /* eslENABLE_AVX */
