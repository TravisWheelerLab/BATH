/* p7_chkoivx_avx.c — AVX2 (256-bit) implementations of P7_CHKOIVX management.
 */
#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <string.h>

#include <immintrin.h>

#include "easel.h"
#include "hmmer.h"
#include "impl_avx.h"


/* ck_vec(): first vector of channel <c> of checkpoint <k>. Float
 * channels are c = 0..allocCF-1, int channels follow them.
 */
static inline __m256 *
ck_vec(const P7_CHKOIVX *cv, int k, int c)
{
  __m256 *base = (__m256 *) (((unsigned long int) ((char *) cv->ck_mem_avx + 31)) & (~0x1f));
  return base + ((size_t) k * (cv->allocCF + cv->allocCI) + c) * cv->allocQ4_avx;
}

/* set_ivx_ptrs(): point the working int channels into ivx_mem_avx. */
static void
set_ivx_ptrs(P7_CHKOIVX *cv, int CI, int nq)
{
  int c;

  if (CI == 0) return;
  cv->ivx_avx[0] = (__m256i *) (((unsigned long int) cv->ivx_mem_avx + 31) & (~0x1f));
  for (c = 1; c < CI; c++)
    cv->ivx_avx[c] = cv->ivx_avx[0] + (c * nq);
}


/*****************************************************************
 * 1. The P7_CHKOIVX AVX2 lifecycle functions.
 *****************************************************************/

/* Function:  p7_chkoivx_Create_avx()
 *
 * Purpose:   Allocate a <P7_CHKOIVX> for a model of up to <M_hint>
 *            nodes using AVX2 (256-bit) vectors: <CI> working int32
 *            channels, and <K> checkpoints that each hold <CF> float
 *            and <CI> int channels. Channels are p7O_NQF_AVX(M_hint)
 *            vectors. SSE pointer fields are initialized to NULL so
 *            the object can be safely destroyed by any ISA path.
 *
 * Returns:   Pointer to the new <P7_CHKOIVX> on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_CHKOIVX *
p7_chkoivx_Create_avx(int M_hint, int CF, int CI, int K)
{
  P7_CHKOIVX *cv;
  int         nq = p7O_NQF_AVX(M_hint);
  int         status;

  cv = NULL;
  ESL_ALLOC(cv, sizeof(P7_CHKOIVX));

  /* SSE fields — not used on the AVX path */
#ifdef eslENABLE_SSE
  cv->ivx     = NULL;
  cv->ivx_mem = NULL;
  cv->ck_mem  = NULL;
  cv->allocQ4 = 0;
#endif

  /* AVX fields */
  cv->ivx_avx     = NULL;
  cv->ivx_mem_avx = NULL;
  cv->ck_mem_avx  = NULL;
  cv->allocQ4_avx = 0;

  cv->allocM  = M_hint;
  cv->allocCF = CF;
  cv->allocCI = CI;
  cv->allocK  = K;

  ESL_ALLOC(cv->ivx_mem_avx, sizeof(__m256i) * CI * nq + 31);
  ESL_ALLOC(cv->ivx_avx,     sizeof(__m256i *) * ESL_MAX(CI, 1));
  ESL_ALLOC(cv->ck_mem_avx,  sizeof(__m256) * (size_t) K * (CF + CI) * nq + 31);

  set_ivx_ptrs(cv, CI, nq);
  cv->allocQ4_avx = nq;

  return cv;

 ERROR:
  p7_chkoivx_Destroy_avx(cv);
  return NULL;
}


/* Function:  p7_chkoivx_GrowTo_avx()
 *
 * Purpose:   Ensure that <cv> (an AVX2-allocated <P7_CHKOIVX>) has
 *            storage for at least <M> model nodes, <CF> float and <CI>
 *            int channels, and <K> checkpoints. Returns immediately if
 *            the existing allocation is sufficient. Otherwise frees and
 *            re-allocates the AVX backing memory and re-establishes the
 *            pointer array.
 *
 *            Any data in <cv> must be considered invalid after a
 *            successful reallocation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; <cv> is destroyed.
 */
int
p7_chkoivx_GrowTo_avx(P7_CHKOIVX *cv, int M, int CF, int CI, int K)
{
  int   nq = p7O_NQF_AVX(M);
  void *p;
  int   status;

  if (nq <= cv->allocQ4_avx && CF <= cv->allocCF && CI <= cv->allocCI && K <= cv->allocK) return eslOK;

  nq = ESL_MAX(nq, cv->allocQ4_avx);
  CF = ESL_MAX(CF, cv->allocCF);
  CI = ESL_MAX(CI, cv->allocCI);
  K  = ESL_MAX(K,  cv->allocK);

  if (cv->ivx_mem_avx != NULL) { free(cv->ivx_mem_avx); cv->ivx_mem_avx = NULL; }
  if (cv->ck_mem_avx  != NULL) { free(cv->ck_mem_avx);  cv->ck_mem_avx  = NULL; }

  if (CI > cv->allocCI) {
    ESL_RALLOC(cv->ivx_avx, p, sizeof(__m256i *) * CI);
  }

  ESL_ALLOC(cv->ivx_mem_avx, sizeof(__m256i) * CI * nq + 31);
  ESL_ALLOC(cv->ck_mem_avx,  sizeof(__m256) * (size_t) K * (CF + CI) * nq + 31);

  set_ivx_ptrs(cv, CI, nq);

  cv->allocM      = ESL_MAX(M, cv->allocM);
  cv->allocCF     = CF;
  cv->allocCI     = CI;
  cv->allocK      = K;
  cv->allocQ4_avx = nq;

  return eslOK;

 ERROR:
  p7_chkoivx_Destroy_avx(cv);
  return status;
}


/* Function:  p7_chkoivx_Destroy_avx()
 *
 * Purpose:   Free a <P7_CHKOIVX> that was allocated by
 *            p7_chkoivx_Create_avx(). Safe to call on NULL.
 */
void
p7_chkoivx_Destroy_avx(P7_CHKOIVX *cv)
{
  if (cv == NULL) return;

  if (cv->ivx_mem_avx != NULL) free(cv->ivx_mem_avx);
  if (cv->ivx_avx     != NULL) free(cv->ivx_avx);
  if (cv->ck_mem_avx  != NULL) free(cv->ck_mem_avx);

  free(cv);
}


/*****************************************************************
 * 2. Saving and loading checkpoints.
 *****************************************************************/

/* Function:  p7_chkoivx_Save_avx()
 *
 * Purpose:   Copy the first <allocCF> float channels of <ov> and the
 *            working int channels of <cv> into checkpoint <k>, for a
 *            model of length <M>.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <k>, <M> or <ov> don't fit <cv>.
 */
int
p7_chkoivx_Save_avx(P7_CHKOIVX *cv, const P7_OIVX *ov, int M, int k)
{
  size_t nbytes = sizeof(__m256) * p7O_NQF_AVX(M);
  int    c;

  if (k < 0 || k >= cv->allocK || p7O_NQF_AVX(M) > cv->allocQ4_avx || ov->allocC < cv->allocCF || p7O_NQF_AVX(M) > ov->allocQ4_avx)
    ESL_EXCEPTION(eslEINVAL, "checkpoint %d for M=%d doesn't fit", k, M);

  for (c = 0; c < cv->allocCF; c++) memcpy(ck_vec(cv, k, c),               ov->ivx_avx[c], nbytes);
  for (c = 0; c < cv->allocCI; c++) memcpy(ck_vec(cv, k, cv->allocCF + c), cv->ivx_avx[c], nbytes);
  return eslOK;
}


/* Function:  p7_chkoivx_Load_avx()
 *
 * Purpose:   Copy checkpoint <k> back into the first <allocCF> float
 *            channels of <ov> and the working int channels of <cv>, for
 *            a model of length <M>.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <k>, <M> or <ov> don't fit <cv>.
 */
int
p7_chkoivx_Load_avx(P7_CHKOIVX *cv, P7_OIVX *ov, int M, int k)
{
  size_t nbytes = sizeof(__m256) * p7O_NQF_AVX(M);
  int    c;

  if (k < 0 || k >= cv->allocK || p7O_NQF_AVX(M) > cv->allocQ4_avx || ov->allocC < cv->allocCF || p7O_NQF_AVX(M) > ov->allocQ4_avx)
    ESL_EXCEPTION(eslEINVAL, "checkpoint %d for M=%d doesn't fit", k, M);

  for (c = 0; c < cv->allocCF; c++) memcpy(ov->ivx_avx[c], ck_vec(cv, k, c),               nbytes);
  for (c = 0; c < cv->allocCI; c++) memcpy(cv->ivx_avx[c], ck_vec(cv, k, cv->allocCF + c), nbytes);
  return eslOK;
}

#endif /* eslENABLE_AVX */
