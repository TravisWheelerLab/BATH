/* OSPLICE_SCORES:
 * Probability-space splice signal scores for SSE multiply-based DP.
 * Mirrors p7_splicescores.c but values are in probability space
 * (1.0 = valid, 0.0 = invalid; signal probs not log-transformed).
 * P-state accumulation arrays are stored as striped __m128 vectors.
 *
 * Contents:
 *    1. Init helpers (signal/donor/acceptor arrays)
 *    2. The OSPLICE_SCORES object.
 */

#include "p7_config.h"

#include <string.h>

#include <xmmintrin.h>
#include <emmintrin.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_sse.h"
#include "../p7_splice.h"


/*****************************************************************
 * 1. Init helpers
 *****************************************************************/

/* Splice signal probabilities taken from
 * "Comprehensive splice-site analysis using comparative genomics",
 * Nihar Sheth et al., 2006
 */
static void
p7_osplicescores_SignalScores(float *f)
{
  f[p7S_GTAG] = 0.9921f;   /* GT-AG */
  f[p7S_GCAG] = 0.0073f;   /* GC-AG */
  f[p7S_ATAC] = 0.0006f;   /* AT-AC */
}

static void
p7_osplicescores_DonorGT(float *f)
{
  /* G = 2; T = 3; x*4 + y */
  int i;
  for (i = 0; i < 16; i++) f[i] = 0.0f;
  f[11] = 1.0f;   /* G*4 + T = 11 */
}

static void
p7_osplicescores_DonorGC(float *f)
{
  /* G = 2; C = 1; x*4 + y */
  int i;
  for (i = 0; i < 16; i++) f[i] = 0.0f;
  f[9] = 1.0f;    /* G*4 + C = 9 */
}

static void
p7_osplicescores_DonorAT(float *f)
{
  /* A = 0; T = 3; x*4 + y */
  int i;
  for (i = 0; i < 16; i++) f[i] = 0.0f;
  f[3] = 1.0f;    /* A*4 + T = 3 */
}

static void
p7_osplicescores_AcceptorAG(float *f)
{
  /* A = 0; G = 2; x*4 + y */
  int i;
  for (i = 0; i < 16; i++) f[i] = 0.0f;
  f[2] = 1.0f;    /* A*4 + G = 2 */
}

static void
p7_osplicescores_AcceptorAC(float *f)
{
  /* A = 0; C = 1; x*4 + y */
  int i;
  for (i = 0; i < 16; i++) f[i] = 0.0f;
  f[1] = 1.0f;    /* A*4 + C = 1 */
}


/*****************************************************************
 * 2. The OSPLICE_SCORES object.
 *****************************************************************/

/* Function:  p7_osplicescores_Create()
 *
 * Purpose:   Allocates an <OSPLICE_SCORES>, initializes probability-space
 *            hardcoded data, and allocates striped __m128 P-state score
 *            storage for an <M_hint> length model.
 *
 * Returns:   Pointer to new <OSPLICE_SCORES> on success.
 *
 * Throws:    <NULL> on allocation error.
 */
OSPLICE_SCORES *
p7_osplicescores_Create(int M_hint)
{
  OSPLICE_SCORES *ss;
  int Q = p7O_NQF(M_hint);
  int i;
  int status;

  ss = NULL;
  ESL_ALLOC(ss, sizeof(OSPLICE_SCORES));

  ss->allocM     = M_hint;
  ss->allocQ     = Q;
  ss->oscore_raw = NULL;
  ss->oscore_base= NULL;
  ss->oscore     = NULL;

  /* Allocate striped __m128 P-state accumulation: SIGNAL_MEM_SIZE arrays of Q vectors.
   * Use +15 trick for 16-byte alignment. */
  ESL_ALLOC(ss->oscore_raw, sizeof(__m128) * SIGNAL_MEM_SIZE * Q + 15);
  ss->oscore_base = (__m128 *)(((unsigned long int) ss->oscore_raw + 15) & (~0xfUL));

  ESL_ALLOC(ss->oscore, sizeof(__m128 *) * SIGNAL_MEM_SIZE);
  for (i = 0; i < SIGNAL_MEM_SIZE; i++)
    ss->oscore[i] = ss->oscore_base + i * Q;

  ss->signal_scores = NULL;
  ESL_ALLOC(ss->signal_scores, sizeof(float) * p7S_SPLICE_SIGNALS);
  p7_osplicescores_SignalScores(ss->signal_scores);

  ss->donor_GT = ss->donor_GC = ss->donor_AT = NULL;
  ss->acceptor_AG = ss->acceptor_AC = NULL;

  ESL_ALLOC(ss->donor_GT,    sizeof(float) * 16);
  ESL_ALLOC(ss->donor_GC,    sizeof(float) * 16);
  ESL_ALLOC(ss->donor_AT,    sizeof(float) * 16);
  ESL_ALLOC(ss->acceptor_AG, sizeof(float) * 16);
  ESL_ALLOC(ss->acceptor_AC, sizeof(float) * 16);

  p7_osplicescores_DonorGT(ss->donor_GT);
  p7_osplicescores_DonorGC(ss->donor_GC);
  p7_osplicescores_DonorAT(ss->donor_AT);
  p7_osplicescores_AcceptorAG(ss->acceptor_AG);
  p7_osplicescores_AcceptorAC(ss->acceptor_AC);

  return ss;

  ERROR:
    p7_osplicescores_Destroy(ss);
    return NULL;
}


/* Function:  p7_osplicescores_GrowTo()
 * Synopsis:  Grow P-state score storage.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_osplicescores_GrowTo(OSPLICE_SCORES *ss, int M)
{
  int Q = p7O_NQF(M);
  void *new_raw;
  int i;
  int status;

  if (M <= ss->allocM) return eslOK;

  ESL_ALLOC(new_raw, sizeof(__m128) * SIGNAL_MEM_SIZE * Q + 15);
  free(ss->oscore_raw);
  ss->oscore_raw  = new_raw;
  ss->oscore_base = (__m128 *)(((unsigned long int) ss->oscore_raw + 15) & (~0xfUL));

  for (i = 0; i < SIGNAL_MEM_SIZE; i++)
    ss->oscore[i] = ss->oscore_base + i * Q;

  ss->allocM = M;
  ss->allocQ = Q;
  return eslOK;

  ERROR:
    p7_osplicescores_Destroy(ss);
    return status;
}


/* Function:  p7_osplicescores_Destroy()
 *
 * Purpose:   Frees an <OSPLICE_SCORES>.
 */
void
p7_osplicescores_Destroy(OSPLICE_SCORES *ss)
{
  if (ss == NULL) return;

  if (ss->oscore     != NULL) free(ss->oscore);
  if (ss->oscore_raw != NULL) free(ss->oscore_raw);

  if (ss->signal_scores != NULL) free(ss->signal_scores);
  if (ss->donor_GT      != NULL) free(ss->donor_GT);
  if (ss->donor_GC      != NULL) free(ss->donor_GC);
  if (ss->donor_AT      != NULL) free(ss->donor_AT);
  if (ss->acceptor_AG   != NULL) free(ss->acceptor_AG);
  if (ss->acceptor_AC   != NULL) free(ss->acceptor_AC);

  free(ss);
}
