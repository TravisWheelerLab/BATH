/* P7_OIVX:
 * Vectorized intermediate-value matrix for the SSE spliced Viterbi.
 * Mirrors P7_IVX (p7_ivx.c) but stores __m128 vectors in striped
 * [channel][stripe] layout so that the inner k loop can be vectorized.
 *
 * Contents:
 *    1. The P7_OIVX object.
 */

#include "p7_config.h"

#include "easel.h"
#include "hmmer.h"
#include "impl_sse.h"

/*****************************************************************
 * 1. The P7_OIVX object.
 *****************************************************************/

/* Function:  p7_oivx_Create()
 *
 * Purpose:   Allocate a <P7_OIVX> for a model of up to <M_hint> nodes
 *            and <C> channels.  The storage layout is [C][Q] where
 *            Q = p7O_NQF(M_hint); element [c][q] holds the striped
 *            __m128 vector for channel c and stripe q.  Backing memory
 *            is over-allocated by 15 bytes so the aligned pointer can
 *            be computed without a separate posix_memalign call,
 *            consistent with the alignment strategy used throughout
 *            impl_sse.
 *
 * Returns:   Pointer to the new <P7_OIVX> on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_OIVX *
p7_oivx_Create(int M_hint, int C)
{
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

  ESL_ALLOC(ov->ivx_mem, sizeof(__m128) * C * nq + 15);
  ESL_ALLOC(ov->ivx,     sizeof(__m128 *) * C);

  ov->ivx[0] = (__m128 *) (((unsigned long int) ov->ivx_mem + 15) & (~0xf));
  for (c = 1; c < C; c++)
    ov->ivx[c] = ov->ivx[0] + (c * nq);

  return ov;

 ERROR:
  p7_oivx_Destroy(ov);
  return NULL;
}


/* Function:  p7_oivx_GrowTo()
 *
 * Purpose:   Ensure that <ov> has storage for at least <M> model nodes
 *            and <C> channels.  If the existing allocation is sufficient,
 *            returns immediately.  Otherwise frees and re-allocates the
 *            backing memory and re-establishes the pointer array.
 *
 *            Any data in <ov> must be considered invalid after a
 *            successful reallocation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; <ov> is destroyed and
 *            should not be used.
 */
int
p7_oivx_GrowTo(P7_OIVX *ov, int M, int C)
{
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
  p7_oivx_Destroy(ov);
  return status;
}


/* Function:  p7_oivx_Destroy()
 *
 * Purpose:   Free a <P7_OIVX>.  Safe to call on NULL.
 */
void
p7_oivx_Destroy(P7_OIVX *ov)
{
  if (ov == NULL) return;

  if (ov->ivx_mem != NULL) free(ov->ivx_mem);
  if (ov->ivx     != NULL) free(ov->ivx);

  free(ov);
}
