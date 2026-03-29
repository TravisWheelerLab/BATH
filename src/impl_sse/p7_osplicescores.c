/* OSPLICE_SCORES:
 * Vectorized splice-signal scores and P-state score storage for the
 * optimized (SSE) spliced Viterbi algorithms.
 *
 * Mirrors the scalar SPLICE_SCORES in p7_splicescores.c, but stores
 * scores as striped __m128 vectors indexed [slot][stripe] so that the
 * inner k loop in the SSE Viterbi can be fully vectorized.
 *
 * Contents:
 *    1. The OSPLICE_SCORES object.
 */

#include "p7_config.h"

#include <math.h>
#include <string.h>

#include "easel.h"

#include "hmmer.h"
#include "impl_sse.h"
#include "p7_splice.h"


/*****************************************************************
 * 1. The OSPLICE_SCORES structure.
 *****************************************************************/


/* Function:  p7_osplicescores_Create()
 *
 * Purpose:   Allocate an <OSPLICE_SCORES>, broadcast the hardcoded splice-signal
 *            log-probabilities into SSE vectors, and allocate P-state score
 *            storage for a model of up to <M_hint> nodes.
 *
 *            The P_scores array is laid out as [SIGNAL_MEM_SIZE][nq] where
 *            nq = p7O_NQF(M_hint), mirroring the [slot][stripe] indexing used
 *            by the OSSX macros.  Memory is over-allocated by 15 bytes so the
 *            aligned pointer can be computed without a separate posix_memalign
 *            call, consistent with the alignment strategy used throughout
 *            impl_sse.
 *
 * Returns:   Pointer to the new <OSPLICE_SCORES> on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
OSPLICE_SCORES *
p7_osplicescores_Create(int M_hint)
{
  OSPLICE_SCORES *os;
  int nq = p7O_NQF(M_hint);
  int s;
  int status;

  os = NULL;
  ESL_ALLOC(os, sizeof(OSPLICE_SCORES));

  os->allocM       = M_hint;
  os->allocQ4      = nq;
  os->P_scores_mem = NULL;
  os->P_scores     = NULL;

  ESL_ALLOC(os->P_scores_mem, sizeof(__m128) * SIGNAL_MEM_SIZE * nq + 15);
  ESL_ALLOC(os->P_scores,     sizeof(__m128 *) * SIGNAL_MEM_SIZE);

  os->P_scores[0] = (__m128 *) (((unsigned long int) os->P_scores_mem + 15) & (~0xf));
  for (s = 1; s < SIGNAL_MEM_SIZE; s++)
    os->P_scores[s] = os->P_scores[0] + (s * nq);

  /* Broadcast splice-signal log-probabilities to all four SSE lanes.
   * Probabilities from Sheth et al. 2006, matching p7_splicescores_SignalScores(). */
  os->signal_scores[p7S_GTAG] = _mm_set1_ps(logf(0.9921f)); /* GT-AG */
  os->signal_scores[p7S_GCAG] = _mm_set1_ps(logf(0.0073f)); /* GC-AG */
  os->signal_scores[p7S_ATAC] = _mm_set1_ps(logf(0.0006f)); /* AT-AC */

  return os;

 ERROR:
  p7_osplicescores_Destroy(os);
  return NULL;
}


/* Function:  p7_osplicescores_GrowTo()
 *
 * Purpose:   Ensure that <os> has P-state score storage for a model of
 *            at least <M> nodes.  If <os->allocM >= M>, returns immediately.
 *            Otherwise frees the existing backing memory, re-allocates at
 *            the new size, and re-establishes the pointer array.
 *
 *            P_scores is always reset to -inf at the start of each Viterbi
 *            call, so the old contents need not be preserved.
 *
 *            The signal_scores vectors are constants and are left unchanged.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; <os> is destroyed and should
 *            not be used.
 */
int
p7_osplicescores_GrowTo(OSPLICE_SCORES *os, int M)
{
  int nq = p7O_NQF(M);
  int s;
  int status;

  if (M <= os->allocM) return eslOK;

  /* Free old backing memory and re-allocate at the new size.
   * The pointer array (P_scores) is always SIGNAL_MEM_SIZE entries; reuse it. */
  if (os->P_scores_mem != NULL) { free(os->P_scores_mem); os->P_scores_mem = NULL; }

  ESL_ALLOC(os->P_scores_mem, sizeof(__m128) * SIGNAL_MEM_SIZE * nq + 15);

  os->P_scores[0] = (__m128 *) (((unsigned long int) os->P_scores_mem + 15) & (~0xf));
  for (s = 1; s < SIGNAL_MEM_SIZE; s++)
    os->P_scores[s] = os->P_scores[0] + (s * nq);

  os->allocM  = M;
  os->allocQ4 = nq;

  return eslOK;

 ERROR:
  p7_osplicescores_Destroy(os);
  return status;
}


/* Function:  p7_osplicescores_Destroy()
 *
 * Purpose:   Free an <OSPLICE_SCORES>.  Safe to call on NULL.
 */
void
p7_osplicescores_Destroy(OSPLICE_SCORES *os)
{
  if (os == NULL) return;

  if (os->P_scores_mem != NULL) free(os->P_scores_mem);
  if (os->P_scores     != NULL) free(os->P_scores);

  free(os);
}
