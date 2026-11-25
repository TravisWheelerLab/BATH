/* p7_chkomx_avx.c — AVX2 (256-bit) P7_CHKOMX matrix lifecycle.
 *
 * Contents:
 *    1. P7_CHKOMX AVX lifecycle.
 *    2. Loading blocks and checkpoints.
 *    3. Cell access.
 */
#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <string.h>

#include <immintrin.h>

#include "easel.h"

#include "hmmer.h"
#include "impl_avx.h"


/* dp_row(): 32-byte aligned storage slot <r> of <cx->dp_mem>. */
static inline __m256 *
dp_row(const P7_CHKOMX *cx, int r)
{
  __m256 *base = (__m256 *) (((unsigned long int) ((char *) cx->dp_mem + 31)) & (~0x1f));
  return base + (size_t) r * cx->allocQ4_avx * cx->nscells;
}

/* ck_row(): 32-byte aligned row <r> of checkpoint <k> in <cx->ck_mem>. */
static inline __m256 *
ck_row(const P7_CHKOMX *cx, int k, int r)
{
  __m256 *base = (__m256 *) (((unsigned long int) ((char *) cx->ck_mem + 31)) & (~0x1f));
  return base + ((size_t) k * cx->saveR + r) * cx->allocQ4_avx * cx->nscells;
}

/* block_end(): last row of the block starting at row <b>. */
static inline int
block_end(const P7_CHKOMX *cx, int b)
{
  return ESL_MIN(cx->L, b + cx->blockR - 1);
}

/* top_block(): first row of the last block, the one a backward
 * computation starts in. */
static inline int
top_block(const P7_CHKOMX *cx)
{
  return (cx->L / cx->blockR) * cx->blockR;
}

/* set_block_ptrs(): point the rows of block <b>, and the <saveR> rows
 * next to it (before it going forward, after it going backward), at
 * their storage slots. Rows of the previously loaded block are set to
 * NULL first. Uses <cx->L> to bound the rows.
 */
static void
set_block_ptrs(P7_CHKOMX *cx, int b)
{
  int i;

  if (cx->b >= 0)
    {
      if (cx->dir == p7_CHK_FORWARD) for (i = ESL_MAX(0, cx->b - cx->saveR); i <= block_end(cx, cx->b);                              i++) cx->dpf_avx[i] = NULL;
      else                           for (i = cx->b;                          i <= ESL_MIN(cx->L, cx->b + cx->blockR + cx->saveR - 1); i++) cx->dpf_avx[i] = NULL;
    }

  if (cx->dir == p7_CHK_FORWARD)
    {
      for (i = ESL_MAX(0, b - cx->saveR); i < b; i++)  cx->dpf_avx[i] = dp_row(cx, i - (b - cx->saveR));
      for (i = b; i <= block_end(cx, b); i++)          cx->dpf_avx[i] = dp_row(cx, cx->saveR + (i - b));
    }
  else
    {
      for (i = b; i <= block_end(cx, b); i++)                                            cx->dpf_avx[i] = dp_row(cx, i - b);
      for (i = b + cx->blockR; i <= ESL_MIN(cx->L, b + cx->blockR + cx->saveR - 1); i++) cx->dpf_avx[i] = dp_row(cx, i - b);
    }

  cx->b = b;
}


/*****************************************************************
 * 1. P7_CHKOMX AVX lifecycle.
 *****************************************************************/

/* Function:  p7_chkomx_Create_avx()
 *
 * Purpose:   Allocate a <P7_CHKOMX> for models up to size <allocM> and
 *            sequences up to length <allocL>, using AVX2 256-bit
 *            vectors, for a computation that runs through the rows in
 *            direction <dir>: <p7_CHK_FORWARD> (rows 0..L, each row
 *            needing rows before it) or <p7_CHK_BACKWARD> (rows L..0,
 *            each row needing rows after it). Rows are computed in
 *            blocks of <blockR> rows, and each checkpoint saves the
 *            <saveR> rows next to a block that its computation reads:
 *            before it going forward, after it going backward.
 *            <saveR> must be <= <blockR>; GrowTo() can change <blockR>. Each row is
 *            p7O_NQF_AVX(allocM) stripes of <nscells> vectors.  SSE
 *            fields are set to NULL and their allocQ fields to 0.
 *
 * Returns:   pointer to new <P7_CHKOMX> on success.
 * Throws:    NULL on allocation failure or bad arguments.
 */
P7_CHKOMX *
p7_chkomx_Create_avx(int allocM, int allocL, int blockR, int saveR, int nscells, int dir)
{
  P7_CHKOMX *cx = NULL;
  int        i;
  int        status;

  ESL_ALLOC(cx, sizeof(P7_CHKOMX));
  cx->dp_mem  = NULL;
  cx->ck_mem  = NULL;
  cx->x_mem   = NULL;
  cx->xmx     = NULL;
#ifdef eslENABLE_SSE
  cx->dpf     = NULL;
#endif
  cx->allocQ4 = 0;
  cx->dpf_avx = NULL;

  if (blockR < 1 || saveR < 0 || saveR > blockR) ESL_XEXCEPTION(eslEINVAL, "bad block size %d or saved rows %d", blockR, saveR);
  if (dir != p7_CHK_FORWARD && dir != p7_CHK_BACKWARD) ESL_XEXCEPTION(eslEINVAL, "bad direction %d", dir);

  cx->dir         = dir;
  cx->nscells     = nscells;
  cx->blockR      = blockR;
  cx->saveR       = saveR;
  cx->allocQ4_avx = p7O_NQF_AVX(allocM);
  cx->allocR      = saveR + blockR;
  cx->nvec        = (size_t) cx->allocR * cx->allocQ4_avx * nscells;
  ESL_ALLOC(cx->dp_mem, sizeof(__m256) * cx->nvec + 31);

  cx->allocP = allocL + 1;
  ESL_ALLOC(cx->dpf_avx, sizeof(__m256 *) * cx->allocP);
  for (i = 0; i < cx->allocP; i++) cx->dpf_avx[i] = NULL;

  cx->allocK = allocL / blockR;
  cx->nckvec = (size_t) cx->allocK * saveR * cx->allocQ4_avx * nscells;
  ESL_ALLOC(cx->ck_mem, sizeof(__m256) * cx->nckvec + 31);

  cx->allocXR = allocL + 1;
  ESL_ALLOC(cx->x_mem, sizeof(float) * cx->allocXR * p7X_NXCELLS + 15);
  cx->xmx = (float *) (((unsigned long int) ((char *) cx->x_mem + 15)) & (~0xf));

  cx->M  = 0;
  cx->L  = 0;
  cx->b  = -1;
  cx->nK = 0;
  return cx;

 ERROR:
  p7_chkomx_Destroy_avx(cx);
  return NULL;
}


/* Function:  p7_chkomx_GrowTo_avx()
 *
 * Purpose:   Ensure that <cx> (an AVX2-allocated <P7_CHKOMX>) can hold a
 *            comparison of a model of size <allocM> to a sequence of
 *            length <allocL> in blocks of <blockR> rows, reallocating as
 *            needed. The direction, saved rows and cells per stripe stay
 *            as they were set by p7_chkomx_Create_avx() or
 *            p7_chkomx_SetDirection(); <saveR> must be <= <blockR>.
 *
 *            Memory is reallocated only when this computation needs more
 *            than is held, with rows of this model's width, so <cx>
 *            holds the most any one computation has needed.
 *
 *            All row pointers are set to NULL and no block is loaded
 *            afterwards; any stored rows and checkpoints are invalid.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on failure, <eslEINVAL> if <blockR> is smaller
 *            than one row or than <saveR>.
 */
int
p7_chkomx_GrowTo_avx(P7_CHKOMX *cx, int allocM, int allocL, int blockR)
{
  void  *p;
  int    nqf    = p7O_NQF_AVX(allocM);
  int    R      = cx->saveR + blockR;
  int    K      = allocL / blockR;
  size_t nvec   = (size_t) R * nqf * cx->nscells;
  size_t nckvec = (size_t) K * cx->saveR * nqf * cx->nscells;
  int    i;
  int    status;

  if (blockR < 1 || cx->saveR > blockR) ESL_EXCEPTION(eslEINVAL, "bad block size %d for %d saved rows", blockR, cx->saveR);

  /* Storage is reallocated only when this computation needs more than
   * is held, and rows take the width of this model, so the memory held
   * is the most any one computation has needed. */
  if (nvec > cx->nvec) {
    ESL_RALLOC(cx->dp_mem, p, sizeof(__m256) * nvec + 31);
    cx->nvec = nvec;
  }
  if (nckvec > cx->nckvec) {
    ESL_RALLOC(cx->ck_mem, p, sizeof(__m256) * nckvec + 31);
    cx->nckvec = nckvec;
  }
  cx->blockR = blockR;
  cx->allocQ4_avx = nqf;
  cx->allocR = (int) (cx->nvec / ((size_t) nqf * cx->nscells));
  cx->allocK = (cx->saveR > 0 ? (int) (cx->nckvec / ((size_t) cx->saveR * nqf * cx->nscells)) : K);

  if (allocL + 1 > cx->allocP) {
    ESL_RALLOC(cx->dpf_avx, p, sizeof(__m256 *) * (allocL + 1));
    cx->allocP = allocL + 1;
  }
  for (i = 0; i < cx->allocP; i++) cx->dpf_avx[i] = NULL;

  if (allocL + 1 > cx->allocXR) {
    ESL_RALLOC(cx->x_mem, p, sizeof(float) * (allocL + 1) * p7X_NXCELLS + 15);
    cx->allocXR = allocL + 1;
    cx->xmx     = (float *) (((unsigned long int) ((char *) cx->x_mem + 15)) & (~0xf));
  }

  cx->M  = 0;
  cx->L  = 0;
  cx->b  = -1;
  cx->nK = 0;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_chkomx_Destroy_avx()
 *
 * Purpose:   Free a <P7_CHKOMX> allocated by p7_chkomx_Create_avx().
 *            Safe to call on NULL.
 */
void
p7_chkomx_Destroy_avx(P7_CHKOMX *cx)
{
  if (cx == NULL) return;

  if (cx->x_mem   != NULL) free(cx->x_mem);
  if (cx->dp_mem  != NULL) free(cx->dp_mem);
  if (cx->ck_mem  != NULL) free(cx->ck_mem);
  if (cx->dpf_avx != NULL) free(cx->dpf_avx);
  free(cx);
}


/*****************************************************************
 * 2. Loading blocks and checkpoints.
 *****************************************************************/

/* Function:  p7_chkomx_StartBlock_avx()
 *
 * Purpose:   Make block <b> (a multiple of <cx->blockR>, <= <cx->L>) the
 *            loaded block, ready to be computed, for a computation that
 *            runs through the blocks in <cx->dir> order.
 *
 *            Forward: blocks go 0, blockR, 2*blockR, ... For <b> = 0
 *            this just points the block's rows at their storage.
 *            Otherwise block <b> - <blockR> must be the loaded block:
 *            the <saveR> rows before <b> are copied into checkpoint
 *            <b>/<blockR> - 1 and into the slots before the block.
 *
 *            Backward: blocks go from the last one (the block holding
 *            row <L>) down to 0. For the last block this just points
 *            its rows at their storage. Otherwise block <b> + <blockR>
 *            must be the loaded block: its first <saveR> rows (fewer if
 *            it ends at <L> sooner) are copied into checkpoint
 *            <b>/<blockR> and into the slots after the block.
 *
 *            Caller sets <cx->M> and <cx->L> before the first call.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <b> isn't the next block in order.
 */
int
p7_chkomx_StartBlock_avx(P7_CHKOMX *cx, int b)
{
  size_t nbytes = sizeof(__m256) * p7O_NQF_AVX(cx->M) * cx->nscells;
  int    k;
  int    r;

  if (cx->dir == p7_CHK_FORWARD)
    {
      if (b == 0) { set_block_ptrs(cx, 0); cx->nK = 0; return eslOK; }

      if (b % cx->blockR != 0 || b != cx->b + cx->blockR || b > cx->L)
        ESL_EXCEPTION(eslEINVAL, "block %d doesn't follow loaded block %d", b, cx->b);

      k = b / cx->blockR - 1;
      for (r = 0; r < cx->saveR; r++)
        {
          memcpy(ck_row(cx, k, r), cx->dpf_avx[b - cx->saveR + r], nbytes);
          memcpy(dp_row(cx, r),    cx->dpf_avx[b - cx->saveR + r], nbytes);
        }
      set_block_ptrs(cx, b);
      cx->nK = k + 1;
    }
  else
    {
      if (b == top_block(cx)) { set_block_ptrs(cx, b); cx->nK = 0; return eslOK; }

      if (b % cx->blockR != 0 || b != cx->b - cx->blockR || b < 0)
        ESL_EXCEPTION(eslEINVAL, "block %d doesn't follow loaded block %d", b, cx->b);

      k = b / cx->blockR;
      for (r = 0; r < cx->saveR && b + cx->blockR + r <= cx->L; r++)
        {
          memcpy(ck_row(cx, k, r),           cx->dpf_avx[b + cx->blockR + r], nbytes);
          memcpy(dp_row(cx, cx->blockR + r), cx->dpf_avx[b + cx->blockR + r], nbytes);
        }
      set_block_ptrs(cx, b);
      cx->nK++;
    }
  return eslOK;
}


/* Function:  p7_chkomx_LoadBlock_avx()
 *
 * Purpose:   Make block <b> (a multiple of <cx->blockR>, <= <cx->L>) the
 *            loaded block, so it can be computed again, in any order.
 *            The block's checkpoint rows are copied into the slots next
 *            to it: forward, checkpoint <b>/<blockR> - 1 before the
 *            block (none for <b> = 0); backward, checkpoint <b>/<blockR>
 *            after it (none for the last block). The rows of block <b>
 *            are pointed at their storage; their contents are left over
 *            from earlier blocks until the caller computes them.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if the block's checkpoint hasn't been saved.
 */
int
p7_chkomx_LoadBlock_avx(P7_CHKOMX *cx, int b)
{
  size_t nbytes = sizeof(__m256) * p7O_NQF_AVX(cx->M) * cx->nscells;
  int    top    = top_block(cx);
  int    k;
  int    r;

  if (b % cx->blockR != 0 || b < 0 || b > cx->L)
    ESL_EXCEPTION(eslEINVAL, "%d isn't the start of a block", b);

  if (cx->dir == p7_CHK_FORWARD)
    {
      if (b / cx->blockR > cx->nK) ESL_EXCEPTION(eslEINVAL, "no checkpoint for block %d", b);
      if (b > 0)
        {
          k = b / cx->blockR - 1;
          for (r = 0; r < cx->saveR; r++)
            memcpy(dp_row(cx, r), ck_row(cx, k, r), nbytes);
        }
    }
  else
    {
      if (b < top - cx->nK * cx->blockR) ESL_EXCEPTION(eslEINVAL, "no checkpoint for block %d", b);
      if (b < top)
        {
          k = b / cx->blockR;
          for (r = 0; r < cx->saveR && b + cx->blockR + r <= cx->L; r++)
            memcpy(dp_row(cx, cx->blockR + r), ck_row(cx, k, r), nbytes);
        }
    }
  set_block_ptrs(cx, b);
  return eslOK;
}


/*****************************************************************
 * 3. Cell access.
 *****************************************************************/

/* Function:  p7_chkomx_FGetMDI_avx()
 * Synopsis:  Read one float DP cell of a loaded row.
 *
 * Purpose:   Return state <s> (p7X_M, p7X_D or p7X_I) of model position
 *            <k> on row <i>. Row <i> must be loaded.
 */
float
p7_chkomx_FGetMDI_avx(const P7_CHKOMX *cx, int s, int i, int k)
{
  union { __m256 v; float p[8]; } u;
  int Q = p7O_NQF_AVX(cx->M);
  u.v = cx->dpf_avx[i][cx->nscells * ((k-1) % Q) + s];
  return u.p[(k-1)/Q];
}

#endif /* eslENABLE_AVX */
