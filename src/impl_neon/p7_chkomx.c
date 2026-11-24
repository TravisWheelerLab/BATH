/* P7_CHKOMX: checkpointed float dynamic programming matrix (NEON).
 *
 * Contents:
 *    1. P7_CHKOMX NEON lifecycle.
 *    2. Loading blocks and checkpoints.
 *    3. Cell access.
 *    4. Unit tests.
 *    5. Test driver.
 */
#include "p7_config.h"

#include <string.h>

#include <arm_neon.h>

#include "easel.h"

#include "hmmer.h"
#include "impl_neon.h"


/* dp_row(): 16-byte aligned storage slot <r> of <cx->dp_mem>. */
static inline float32x4_t *
dp_row(const P7_CHKOMX *cx, int r)
{
  float32x4_t *base = (float32x4_t *) (((unsigned long int) ((char *) cx->dp_mem + 15)) & (~0xf));
  return base + (size_t) r * cx->allocQ4 * cx->nscells;
}

/* ck_row(): 16-byte aligned row <r> of checkpoint <k> in <cx->ck_mem>. */
static inline float32x4_t *
ck_row(const P7_CHKOMX *cx, int k, int r)
{
  float32x4_t *base = (float32x4_t *) (((unsigned long int) ((char *) cx->ck_mem + 15)) & (~0xf));
  return base + ((size_t) k * cx->saveR + r) * cx->allocQ4 * cx->nscells;
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
      if (cx->dir == p7_CHK_FORWARD) for (i = ESL_MAX(0, cx->b - cx->saveR); i <= block_end(cx, cx->b);                              i++) cx->dpf[i] = NULL;
      else                           for (i = cx->b;                          i <= ESL_MIN(cx->L, cx->b + cx->blockR + cx->saveR - 1); i++) cx->dpf[i] = NULL;
    }

  if (cx->dir == p7_CHK_FORWARD)
    {
      for (i = ESL_MAX(0, b - cx->saveR); i < b; i++)  cx->dpf[i] = dp_row(cx, i - (b - cx->saveR));
      for (i = b; i <= block_end(cx, b); i++)          cx->dpf[i] = dp_row(cx, cx->saveR + (i - b));
    }
  else
    {
      for (i = b; i <= block_end(cx, b); i++)                                            cx->dpf[i] = dp_row(cx, i - b);
      for (i = b + cx->blockR; i <= ESL_MIN(cx->L, b + cx->blockR + cx->saveR - 1); i++) cx->dpf[i] = dp_row(cx, i - b);
    }

  cx->b = b;
}


/*****************************************************************
 * 1. P7_CHKOMX NEON lifecycle.
 *****************************************************************/

/* Function:  p7_chkomx_Create()
 *
 * Purpose:   Allocate a <P7_CHKOMX> for models up to size <allocM> and
 *            sequences up to length <allocL>, using NEON 128-bit
 *            vectors, for a computation that runs through the rows in
 *            direction <dir>: <p7_CHK_FORWARD> (rows 0..L, each row
 *            needing rows before it) or <p7_CHK_BACKWARD> (rows L..0,
 *            each row needing rows after it). Rows are computed in
 *            blocks of <blockR> rows, and each checkpoint saves the
 *            <saveR> rows next to a block that its computation reads:
 *            before it going forward, after it going backward.
 *            <saveR> must be <= <blockR>; GrowTo() can change <blockR>. Each row is
 *            p7O_NQF(allocM) stripes of <nscells> vectors.
 *
 * Returns:   pointer to new <P7_CHKOMX> on success.
 * Throws:    NULL on allocation failure or bad arguments.
 */
P7_CHKOMX *
p7_chkomx_Create(int allocM, int allocL, int blockR, int saveR, int nscells, int dir)
{
  P7_CHKOMX *cx = NULL;
  int        i;
  int        status;

  ESL_ALLOC(cx, sizeof(P7_CHKOMX));
  cx->dp_mem  = NULL;
  cx->ck_mem  = NULL;
  cx->x_mem   = NULL;
  cx->xmx     = NULL;
  cx->dpf     = NULL;

  if (blockR < 1 || saveR < 0 || saveR > blockR) ESL_XEXCEPTION(eslEINVAL, "bad block size %d or saved rows %d", blockR, saveR);
  if (dir != p7_CHK_FORWARD && dir != p7_CHK_BACKWARD) ESL_XEXCEPTION(eslEINVAL, "bad direction %d", dir);

  cx->dir         = dir;
  cx->nscells     = nscells;
  cx->blockR      = blockR;
  cx->saveR       = saveR;
  cx->allocQ4 = p7O_NQF(allocM);
  cx->allocR      = saveR + blockR;
  cx->nvec        = (size_t) cx->allocR * cx->allocQ4 * nscells;
  ESL_ALLOC(cx->dp_mem, sizeof(float32x4_t) * cx->nvec + 15);

  cx->allocP = allocL + 1;
  ESL_ALLOC(cx->dpf, sizeof(float32x4_t *) * cx->allocP);
  for (i = 0; i < cx->allocP; i++) cx->dpf[i] = NULL;

  cx->allocK = allocL / blockR;
  ESL_ALLOC(cx->ck_mem, sizeof(float32x4_t) * (size_t) cx->allocK * saveR * cx->allocQ4 * nscells + 15);

  cx->allocXR = allocL + 1;
  ESL_ALLOC(cx->x_mem, sizeof(float) * cx->allocXR * p7X_NXCELLS + 15);
  cx->xmx = (float *) (((unsigned long int) ((char *) cx->x_mem + 15)) & (~0xf));

  cx->M  = 0;
  cx->L  = 0;
  cx->b  = -1;
  cx->nK = 0;
  return cx;

 ERROR:
  p7_chkomx_Destroy(cx);
  return NULL;
}


/* Function:  p7_chkomx_GrowTo()
 *
 * Purpose:   Ensure that <cx> can hold a
 *            comparison of a model of size <allocM> to a sequence of
 *            length <allocL> in blocks of <blockR> rows, reallocating as
 *            needed. The direction, saved rows and cells per stripe stay
 *            as they were set by p7_chkomx_Create(); <saveR> must be <= <blockR>.
 *
 *            All row pointers are set to NULL and no block is loaded
 *            afterwards; any stored rows and checkpoints are invalid.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on failure, <eslEINVAL> if <blockR> is smaller
 *            than one row or than <saveR>.
 */
int
p7_chkomx_GrowTo(P7_CHKOMX *cx, int allocM, int allocL, int blockR)
{
  void  *p;
  int    nqf  = ESL_MAX(cx->allocQ4, p7O_NQF(allocM));
  int    R    = ESL_MAX(cx->allocR, cx->saveR + blockR);
  int    K    = allocL / blockR;
  size_t nvec = (size_t) R * nqf * cx->nscells;
  int    i;
  int    status;

  if (blockR < 1 || cx->saveR > blockR) ESL_EXCEPTION(eslEINVAL, "bad block size %d for %d saved rows", blockR, cx->saveR);
  cx->blockR = blockR;
  cx->allocR = R;

  if (nvec > cx->nvec) {
    ESL_RALLOC(cx->dp_mem, p, sizeof(float32x4_t) * nvec + 15);
    cx->nvec = nvec;
  }

  if (K > cx->allocK || nqf > cx->allocQ4) {
    ESL_RALLOC(cx->ck_mem, p, sizeof(float32x4_t) * (size_t) ESL_MAX(K, cx->allocK) * cx->saveR * nqf * cx->nscells + 15);
    cx->allocK = ESL_MAX(K, cx->allocK);
  }
  cx->allocQ4 = nqf;

  if (allocL + 1 > cx->allocP) {
    ESL_RALLOC(cx->dpf, p, sizeof(float32x4_t *) * (allocL + 1));
    cx->allocP = allocL + 1;
  }
  for (i = 0; i < cx->allocP; i++) cx->dpf[i] = NULL;

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


/* Function:  p7_chkomx_Destroy()
 *
 * Purpose:   Free a <P7_CHKOMX> allocated by p7_chkomx_Create().
 *            Safe to call on NULL.
 */
void
p7_chkomx_Destroy(P7_CHKOMX *cx)
{
  if (cx == NULL) return;

  if (cx->x_mem   != NULL) free(cx->x_mem);
  if (cx->dp_mem  != NULL) free(cx->dp_mem);
  if (cx->ck_mem  != NULL) free(cx->ck_mem);
  if (cx->dpf != NULL) free(cx->dpf);
  free(cx);
}


/*****************************************************************
 * 2. Loading blocks and checkpoints.
 *****************************************************************/

/* Function:  p7_chkomx_StartBlock()
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
p7_chkomx_StartBlock(P7_CHKOMX *cx, int b)
{
  size_t nbytes = sizeof(float32x4_t) * p7O_NQF(cx->M) * cx->nscells;
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
          memcpy(ck_row(cx, k, r), cx->dpf[b - cx->saveR + r], nbytes);
          memcpy(dp_row(cx, r),    cx->dpf[b - cx->saveR + r], nbytes);
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
          memcpy(ck_row(cx, k, r),           cx->dpf[b + cx->blockR + r], nbytes);
          memcpy(dp_row(cx, cx->blockR + r), cx->dpf[b + cx->blockR + r], nbytes);
        }
      set_block_ptrs(cx, b);
      cx->nK++;
    }
  return eslOK;
}


/* Function:  p7_chkomx_LoadBlock()
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
p7_chkomx_LoadBlock(P7_CHKOMX *cx, int b)
{
  size_t nbytes = sizeof(float32x4_t) * p7O_NQF(cx->M) * cx->nscells;
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

/* Function:  p7_chkomx_FGetMDI()
 * Synopsis:  Read one float DP cell of a loaded row.
 *
 * Purpose:   Return state <s> (p7X_M, p7X_D or p7X_I) of model position
 *            <k> on row <i>. Row <i> must be loaded.
 */
float
p7_chkomx_FGetMDI(const P7_CHKOMX *cx, int s, int i, int k)
{
  union { float32x4_t v; float p[4]; } u;
  int Q = p7O_NQF(cx->M);
  u.v = cx->dpf[i][cx->nscells * ((k-1) % Q) + s];
  return u.p[(k-1)/Q];
}


/*****************************************************************
 * 4. Unit tests.
 *****************************************************************/
#ifdef p7CHKOMX_TESTDRIVE

/* Each row's floats are set to row * 10000 + position in the row, so
 * any row can be checked against what was written to it. <V> is the
 * number of floats per vector for the ISA under test.
 */
typedef struct {
  P7_CHKOMX *(*Create)    (int, int, int, int, int, int);
  int        (*GrowTo)    (P7_CHKOMX *, int, int, int);
  int        (*StartBlock)(P7_CHKOMX *, int);
  int        (*LoadBlock) (P7_CHKOMX *, int);
  float      (*FGetMDI)   (const P7_CHKOMX *, int, int, int);
  void       (*Destroy)   (P7_CHKOMX *);
  int          V;
  const char  *name;
} CHKOMX_ISA;

static float *
row_floats(const P7_CHKOMX *cx, const CHKOMX_ISA *isa, int i)
{
  return (float *) cx->dpf[i];
}

static void
write_row(P7_CHKOMX *cx, const CHKOMX_ISA *isa, int i, int nfloat)
{
  float *p = row_floats(cx, isa, i);
  int    n;

  if (p == NULL) esl_fatal("%s: row %d isn't loaded", isa->name, i);
  for (n = 0; n < nfloat; n++) p[n] = (float) (i * 10000 + n);
}

static void
check_row(const P7_CHKOMX *cx, const CHKOMX_ISA *isa, int i, int nfloat)
{
  float *p = row_floats(cx, isa, i);
  int    n;

  if (p == NULL) esl_fatal("%s: row %d isn't loaded", isa->name, i);
  for (n = 0; n < nfloat; n++)
    if (p[n] != (float) (i * 10000 + n)) esl_fatal("%s: row %d holds wrong values", isa->name, i);
}

/* check_loaded(): the rows of block <b> and its saved rows (before it
 * going forward, after it going backward) are loaded, all others NULL. */
static void
check_loaded(const P7_CHKOMX *cx, const CHKOMX_ISA *isa, int b)
{
  int i, lo, hi;

  if (cx->dir == p7_CHK_FORWARD) { lo = b - cx->saveR; hi = b + cx->blockR - 1;             }
  else                           { lo = b;             hi = b + cx->blockR + cx->saveR - 1; }
  for (i = 0; i <= cx->L; i++)
    {
      int loaded = (i >= lo && i <= hi);
      if (loaded != (row_floats(cx, isa, i) != NULL)) esl_fatal("%s: row %d pointer wrong for block %d", isa->name, i, b);
    }
}

/* utest_blocks():
 * Run forward through all blocks of an <M> by <L> comparison, writing
 * each row as it is "computed", then load the blocks again from last
 * to first, rewriting each block, and check that the saved rows before
 * each block and the block's own rows hold the right values. Also
 * checks p7_chkomx_FGetMDI() against the striped layout.
 */
static void
utest_blocks(P7_CHKOMX *cx, const CHKOMX_ISA *isa, int M, int L)
{
  int Q      = ESL_MAX(2, (M - 1) / isa->V + 1);
  int nfloat = Q * cx->nscells * isa->V;
  int b, i, k, s;

  cx->M = M;
  cx->L = L;

  for (b = 0; b <= L; b += cx->blockR)
    {
      if (isa->StartBlock(cx, b) != eslOK) esl_fatal("%s: StartBlock(%d) failed", isa->name, b);
      check_loaded(cx, isa, b);
      for (i = ESL_MAX(0, b - cx->saveR); i < b; i++) check_row(cx, isa, i, nfloat);
      for (i = b; i <= ESL_MIN(L, b + cx->blockR - 1); i++) write_row(cx, isa, i, nfloat);
    }
  if (cx->nK != L / cx->blockR) esl_fatal("%s: %d checkpoints saved, expected %d", isa->name, cx->nK, L / cx->blockR);

  for (b = (L / cx->blockR) * cx->blockR; b >= 0; b -= cx->blockR)
    {
      if (isa->LoadBlock(cx, b) != eslOK) esl_fatal("%s: LoadBlock(%d) failed", isa->name, b);
      check_loaded(cx, isa, b);
      for (i = ESL_MAX(0, b - cx->saveR); i < b; i++) check_row(cx, isa, i, nfloat);
      for (i = b; i <= ESL_MIN(L, b + cx->blockR - 1); i++) write_row(cx, isa, i, nfloat);
      for (i = b; i <= ESL_MIN(L, b + cx->blockR - 1); i++) check_row(cx, isa, i, nfloat);

      /* FGetMDI reads float lane (k-1)/Q of vector nscells*((k-1)%Q)+s */
      i = b;
      for (k = 1; k <= M; k++)
        for (s = 0; s < cx->nscells; s++)
          {
            int n = (cx->nscells * ((k-1) % Q) + s) * isa->V + (k-1) / Q;
            if (isa->FGetMDI(cx, s, i, k) != (float) (i * 10000 + n)) esl_fatal("%s: FGetMDI wrong at i=%d k=%d s=%d", isa->name, i, k, s);
          }
    }

  /* A block that doesn't follow the loaded one is refused */
  if (L >= 2 * cx->blockR && isa->StartBlock(cx, 2 * cx->blockR) != eslEINVAL) esl_fatal("%s: out-of-order StartBlock accepted", isa->name);
}

/* utest_blocks_backward():
 * Backward direction. Run backward through all blocks, from the one
 * holding row L down to 0, writing each block's rows from its last row
 * down, and check the saved rows after each block. Then load the blocks
 * again from first to last, as a forward sweep reading the rows in
 * increasing order would, rewriting each block, and check the saved rows
 * and the block's own rows.
 */
static void
utest_blocks_backward(P7_CHKOMX *cx, const CHKOMX_ISA *isa, int M, int L)
{
  int Q      = ESL_MAX(2, (M - 1) / isa->V + 1);
  int nfloat = Q * cx->nscells * isa->V;
  int top, b, i;

  cx->M = M;
  cx->L = L;
  top   = (L / cx->blockR) * cx->blockR;

  for (b = top; b >= 0; b -= cx->blockR)
    {
      if (isa->StartBlock(cx, b) != eslOK) esl_fatal("%s: backward StartBlock(%d) failed", isa->name, b);
      check_loaded(cx, isa, b);
      for (i = b + cx->blockR; i <= ESL_MIN(L, b + cx->blockR + cx->saveR - 1); i++) check_row(cx, isa, i, nfloat);
      for (i = ESL_MIN(L, b + cx->blockR - 1); i >= b; i--) write_row(cx, isa, i, nfloat);
    }
  if (cx->nK != L / cx->blockR) esl_fatal("%s: %d checkpoints saved, expected %d", isa->name, cx->nK, L / cx->blockR);

  for (b = 0; b <= top; b += cx->blockR)
    {
      if (isa->LoadBlock(cx, b) != eslOK) esl_fatal("%s: backward LoadBlock(%d) failed", isa->name, b);
      check_loaded(cx, isa, b);
      for (i = b + cx->blockR; i <= ESL_MIN(L, b + cx->blockR + cx->saveR - 1); i++) check_row(cx, isa, i, nfloat);
      for (i = ESL_MIN(L, b + cx->blockR - 1); i >= b; i--) write_row(cx, isa, i, nfloat);
      for (i = b; i <= ESL_MIN(L, b + cx->blockR - 1); i++) check_row(cx, isa, i, nfloat);
    }

  /* A block that doesn't follow the loaded one is refused */
  if (top >= 2 * cx->blockR)
    {
      if (isa->StartBlock(cx, top) != eslOK) esl_fatal("%s: backward StartBlock(top) failed", isa->name);
      if (isa->StartBlock(cx, top - 2 * cx->blockR) != eslEINVAL) esl_fatal("%s: out-of-order backward StartBlock accepted", isa->name);
    }
}

static void
utest_isa(const CHKOMX_ISA *isa)
{
  P7_CHKOMX *cx = NULL;

  if ((cx = isa->Create(13, 37, 8, 3, p7X_NSCELLS, p7_CHK_FORWARD)) == NULL) esl_fatal("%s: Create failed", isa->name);
  utest_blocks(cx, isa, 13, 37);
  utest_blocks(cx, isa, 13, 32);                   /* L a multiple of the block size */
  utest_blocks(cx, isa,  5,  6);                   /* a single block                 */

  if (isa->GrowTo(cx, 70, 101, 8) != eslOK) esl_fatal("%s: GrowTo failed", isa->name);
  utest_blocks(cx, isa, 70, 101);
  if (isa->GrowTo(cx, 20, 40, 8) != eslOK) esl_fatal("%s: GrowTo (smaller) failed", isa->name);
  utest_blocks(cx, isa, 20, 40);
  isa->Destroy(cx);

  if ((cx = isa->Create(7, 50, 8, 8, p7X_NSCELLS, p7_CHK_FORWARD)) == NULL) esl_fatal("%s: Create failed", isa->name);
  utest_blocks(cx, isa, 7, 50);                    /* saveR == blockR                */
  isa->Destroy(cx);

  if ((cx = isa->Create(13, 40, 8, 5, p7X_NSCELLS, p7_CHK_BACKWARD)) == NULL) esl_fatal("%s: backward Create failed", isa->name);
  utest_blocks_backward(cx, isa, 13, 37);          /* last block shorter than saveR  */
  utest_blocks_backward(cx, isa, 13, 40);          /* last block is the single row L */
  utest_blocks_backward(cx, isa, 13, 39);          /* last block one row short       */
  utest_blocks_backward(cx, isa,  5,  6);          /* a single block                 */
  if (isa->GrowTo(cx, 30, 200, 32) != eslOK) esl_fatal("%s: backward GrowTo failed", isa->name);
  utest_blocks_backward(cx, isa, 30, 200);
  if (isa->GrowTo(cx, 30, 200, 5)  != eslOK) esl_fatal("%s: backward GrowTo (block == saveR) failed", isa->name);
  utest_blocks_backward(cx, isa, 30, 200);
  if (isa->GrowTo(cx, 70, 140, 13) != eslOK) esl_fatal("%s: backward GrowTo (bigger block, wider rows) failed", isa->name);
  utest_blocks_backward(cx, isa, 70, 140);
  if (isa->GrowTo(cx, 70, 140, 4)  != eslEINVAL) esl_fatal("%s: GrowTo with blockR < saveR accepted", isa->name);
  isa->Destroy(cx);

  if (isa->Create(7, 50, 4, 5, p7X_NSCELLS, p7_CHK_FORWARD) != NULL) esl_fatal("%s: saveR > blockR accepted", isa->name);
  if (isa->Create(7, 50, 8, 5, p7X_NSCELLS, 2)              != NULL) esl_fatal("%s: bad direction accepted", isa->name);
}

#endif /*p7CHKOMX_TESTDRIVE*/


/*****************************************************************
 * 5. Test driver.
 *****************************************************************/
#ifdef p7CHKOMX_TESTDRIVE
/*
   make chkomx_utest
   ./chkomx_utest
 */
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default   env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the P7_CHKOMX checkpointed matrix";

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);

  esl_exception_SetHandler(&esl_nonfatal_handler);   /* the tests check that bad arguments are refused */

  { CHKOMX_ISA neon = { p7_chkomx_Create, p7_chkomx_GrowTo, p7_chkomx_StartBlock, p7_chkomx_LoadBlock, p7_chkomx_FGetMDI, p7_chkomx_Destroy, 4, "NEON" };
    utest_isa(&neon); }

  esl_getopts_Destroy(go);
  printf("All tests passed.\n");
  return eslOK;
}
#endif /*p7CHKOMX_TESTDRIVE*/
