/* P7_CHKGMX : a checkpointed generic dynamic programming matrix
 *
 * Contents:
 *   1. The <P7_CHKGMX> object
 *   2. Loading blocks and checkpoints
 *   3. Unit tests
 *   4. Test driver
 */
#include "p7_config.h"

#include <math.h>
#include <string.h>

#include "easel.h"

#include "hmmer.h"


/* dp_row(): storage slot <r> of <cx->dp_mem>. */
static inline float *
dp_row(const P7_CHKGMX *cx, int r)
{
  return cx->dp_mem + (ptrdiff_t) r * (ptrdiff_t) cx->allocW * (ptrdiff_t) cx->nscells;
}

/* ck_row(): row <r> of checkpoint <k> in <cx->ck_mem>. */
static inline float *
ck_row(const P7_CHKGMX *cx, int k, int r)
{
  return cx->ck_mem + ((ptrdiff_t) k * cx->saveR + r) * (ptrdiff_t) cx->allocW * (ptrdiff_t) cx->nscells;
}

/* block_end(): last row of the block starting at row <b>. */
static inline int
block_end(const P7_CHKGMX *cx, int b)
{
  return ESL_MIN(cx->L, b + cx->blockR - 1);
}

/* top_block(): first row of the last block, the one a backward
 * computation starts in. */
static inline int
top_block(const P7_CHKGMX *cx)
{
  return (cx->L / cx->blockR) * cx->blockR;
}

/* set_block_ptrs(): point the rows of block <b>, and the <saveR> rows
 * next to it (before it going forward, after it going backward), at
 * their storage slots. Rows of the previously loaded block are set to
 * NULL first. Uses <cx->L> to bound the rows.
 */
static void
set_block_ptrs(P7_CHKGMX *cx, int b)
{
  int i;

  if (cx->b >= 0)
    {
      if (cx->dir == p7_CHK_FORWARD) for (i = ESL_MAX(0, cx->b - cx->saveR); i <= block_end(cx, cx->b);                         i++) cx->dp[i] = NULL;
      else                           for (i = cx->b;                          i <= ESL_MIN(cx->L, cx->b + cx->blockR + cx->saveR - 1); i++) cx->dp[i] = NULL;
    }

  if (cx->dir == p7_CHK_FORWARD)
    {
      for (i = ESL_MAX(0, b - cx->saveR); i < b; i++)  cx->dp[i] = dp_row(cx, i - (b - cx->saveR));
      for (i = b; i <= block_end(cx, b); i++)          cx->dp[i] = dp_row(cx, cx->saveR + (i - b));
    }
  else
    {
      for (i = b; i <= block_end(cx, b); i++)                                         cx->dp[i] = dp_row(cx, i - b);
      for (i = b + cx->blockR; i <= ESL_MIN(cx->L, b + cx->blockR + cx->saveR - 1); i++) cx->dp[i] = dp_row(cx, i - b);
    }

  cx->b = b;
}


/*****************************************************************
 *= 1. The <P7_CHKGMX> object.
 *****************************************************************/

/* Function:  p7_chkgmx_Create()
 * Synopsis:  Allocate a new <P7_CHKGMX>.
 *
 * Purpose:   Allocate a reusable, resizeable <P7_CHKGMX> for models up
 *            to size <allocM> and sequences up to length <allocL>, with
 *            <nscells> cells per model position, for a computation that
 *            runs through the rows in direction <dir>: <p7_CHK_FORWARD>
 *            (rows 0..L, each row needing rows before it) or
 *            <p7_CHK_BACKWARD> (rows L..0, each row needing rows after
 *            it). Rows are computed in blocks of <blockR> rows, and each
 *            checkpoint saves the <saveR> rows next to a block that the
 *            computation of the block reads: before it going forward,
 *            after it going backward. <saveR> must be <= <blockR>.
 *            <p7_chkgmx_GrowTo()> can change <blockR> later.
 *
 * Returns:   a pointer to the new <P7_CHKGMX>.
 *
 * Throws:    <NULL> on allocation error or bad arguments.
 */
P7_CHKGMX *
p7_chkgmx_Create(int allocM, int allocL, int blockR, int saveR, int nscells, int dir)
{
  P7_CHKGMX *cx = NULL;
  int        i;
  int        status;

  ESL_ALLOC(cx, sizeof(P7_CHKGMX));
  cx->dp     = NULL;
  cx->dp_mem = NULL;
  cx->ck_mem = NULL;
  cx->xmx    = NULL;

  if (blockR < 1 || saveR < 0 || saveR > blockR) ESL_XEXCEPTION(eslEINVAL, "bad block size %d or saved rows %d", blockR, saveR);
  if (dir != p7_CHK_FORWARD && dir != p7_CHK_BACKWARD) ESL_XEXCEPTION(eslEINVAL, "bad direction %d", dir);

  cx->dir     = dir;
  cx->nscells = nscells;
  cx->blockR  = blockR;
  cx->saveR   = saveR;
  cx->allocW  = allocM + 1;
  cx->allocR  = saveR + blockR;
  ESL_ALLOC(cx->dp_mem, sizeof(float) * cx->allocR * cx->allocW * nscells);

  cx->allocP = allocL + 1;
  ESL_ALLOC(cx->dp, sizeof(float *) * cx->allocP);
  for (i = 0; i < cx->allocP; i++) cx->dp[i] = NULL;

  cx->allocK = allocL / blockR;
  ESL_ALLOC(cx->ck_mem, sizeof(float) * ESL_MAX(1, (size_t) cx->allocK * saveR * cx->allocW * nscells));

  cx->allocXR = allocL + 1;
  ESL_ALLOC(cx->xmx, sizeof(float) * cx->allocXR * p7G_NXCELLS);

  cx->M  = 0;
  cx->L  = 0;
  cx->b  = -1;
  cx->nK = 0;
  return cx;

 ERROR:
  p7_chkgmx_Destroy(cx);
  return NULL;
}


/* Function:  p7_chkgmx_GrowTo()
 * Synopsis:  Assure that a checkpointed DP matrix is big enough.
 *
 * Purpose:   Assures that <cx> can hold a comparison of a model of size
 *            up to <M> to a sequence of length up to <L> in blocks of
 *            <blockR> rows; reallocates if necessary. The direction,
 *            saved rows and cells per model position stay as they were
 *            set by p7_chkgmx_Create(); <saveR> must be <= <blockR>.
 *
 *            All row pointers are set to NULL and no block is loaded
 *            afterwards; any stored rows and checkpoints are invalid.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure, <eslEINVAL> if <blockR>
 *            is smaller than one row or than <saveR>.
 */
int
p7_chkgmx_GrowTo(P7_CHKGMX *cx, int M, int L, int blockR)
{
  void *p;
  int   W = ESL_MAX(cx->allocW, M + 1);
  int   R = ESL_MAX(cx->allocR, cx->saveR + blockR);
  int   K = L / blockR;
  int   i;
  int   status;

  if (blockR < 1 || cx->saveR > blockR) ESL_EXCEPTION(eslEINVAL, "bad block size %d for %d saved rows", blockR, cx->saveR);
  cx->blockR = blockR;

  if (W > cx->allocW || R > cx->allocR) {
    ESL_RALLOC(cx->dp_mem, p, sizeof(float) * R * W * cx->nscells);
    cx->allocR = R;
  }

  if (K > cx->allocK || W > cx->allocW) {
    ESL_RALLOC(cx->ck_mem, p, sizeof(float) * ESL_MAX(1, (size_t) ESL_MAX(K, cx->allocK) * cx->saveR * W * cx->nscells));
    cx->allocK = ESL_MAX(K, cx->allocK);
  }
  cx->allocW = W;

  if (L + 1 > cx->allocP) {
    ESL_RALLOC(cx->dp, p, sizeof(float *) * (L + 1));
    cx->allocP = L + 1;
  }
  for (i = 0; i < cx->allocP; i++) cx->dp[i] = NULL;

  if (L + 1 > cx->allocXR) {
    ESL_RALLOC(cx->xmx, p, sizeof(float) * (L + 1) * p7G_NXCELLS);
    cx->allocXR = L + 1;
  }

  cx->M  = 0;
  cx->L  = 0;
  cx->b  = -1;
  cx->nK = 0;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_chkgmx_Destroy()
 * Synopsis:  Frees a checkpointed DP matrix.
 *
 * Purpose:   Frees a <P7_CHKGMX>. Safe to call on NULL.
 *
 * Returns:   (void)
 */
void
p7_chkgmx_Destroy(P7_CHKGMX *cx)
{
  if (cx == NULL) return;

  if (cx->dp     != NULL) free(cx->dp);
  if (cx->dp_mem != NULL) free(cx->dp_mem);
  if (cx->ck_mem != NULL) free(cx->ck_mem);
  if (cx->xmx    != NULL) free(cx->xmx);
  free(cx);
  return;
}


/* Function:  p7_chk_BlockRows()
 * Synopsis:  Rows per block for a checkpointed matrix.
 *
 * Purpose:   Return the number of rows per block for a checkpointed
 *            matrix (<P7_CHKGMX> or <P7_CHKOMX>) for a model of length
 *            <M>, a sequence of length <L>, <saveR> rows saved at each
 *            checkpoint and <nscells> floats per cell: one block of all
 *            <L>+1 rows if that block and its saved rows take at most
 *            <p7_CHK_ONEBLOCK_BYTES>, else ceil(sqrt(saveR*L)) rows, which
 *            keeps the block and the checkpoints at about
 *            2*sqrt(saveR*L) rows together.
 */
int
p7_chk_BlockRows(int M, int L, int saveR, int nscells)
{
  double bytes  = (double) (L + 1 + saveR) * (double) M * (double) nscells * sizeof(float);
  int    blockR;

  if (bytes <= p7_CHK_ONEBLOCK_BYTES) return L + 1;
  blockR = (int) ceil(sqrt((double) saveR * (double) L));
  return ESL_MAX(blockR, saveR);
}

/* Function:  p7_chkgmx_SetDirection()
 * Synopsis:  Change the direction of a checkpointed matrix.
 *
 * Purpose:   Set up <cx> for a computation that runs through the rows in
 *            direction <dir> (<p7_CHK_FORWARD> or <p7_CHK_BACKWARD>),
 *            saving <saveR> rows at each checkpoint, so that one
 *            allocation can serve a Forward and a Backward computation
 *            in turn. Afterwards no block is loaded and no checkpoints
 *            are saved, and the caller must call <p7_chkgmx_GrowTo()> before the
 *            next computation; that call also makes room for a larger
 *            <saveR> and checks it against the block size.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on a bad direction or <saveR> < 1.
 */
int
p7_chkgmx_SetDirection(P7_CHKGMX *cx, int dir, int saveR)
{
  if (dir != p7_CHK_FORWARD && dir != p7_CHK_BACKWARD) ESL_EXCEPTION(eslEINVAL, "bad direction %d", dir);
  if (saveR < 1)                                       ESL_EXCEPTION(eslEINVAL, "saveR must be >= 1, not %d", saveR);
  cx->dir   = dir;
  cx->saveR = saveR;
  cx->M     = 0;
  cx->L     = 0;
  cx->b     = -1;
  cx->nK    = 0;
  return eslOK;
}

/*****************************************************************
 *= 2. Loading blocks and checkpoints.
 *****************************************************************/

/* Function:  p7_chkgmx_StartBlock()
 * Synopsis:  Move a computation into its next block.
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
 *            must be the loaded block: the first <saveR> rows of it
 *            (fewer if it ends at <L> sooner) are copied into checkpoint
 *            <b>/<blockR> and into the slots after the block.
 *
 *            Caller sets <cx->M> and <cx->L> before the first call.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <b> isn't the next block in order.
 */
int
p7_chkgmx_StartBlock(P7_CHKGMX *cx, int b)
{
  size_t nbytes = sizeof(float) * (cx->M + 1) * cx->nscells;
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
          memcpy(ck_row(cx, k, r), cx->dp[b - cx->saveR + r], nbytes);
          memcpy(dp_row(cx, r),    cx->dp[b - cx->saveR + r], nbytes);
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
          memcpy(ck_row(cx, k, r),              cx->dp[b + cx->blockR + r], nbytes);
          memcpy(dp_row(cx, cx->blockR + r),    cx->dp[b + cx->blockR + r], nbytes);
        }
      set_block_ptrs(cx, b);
      cx->nK++;
    }
  return eslOK;
}


/* Function:  p7_chkgmx_LoadBlock()
 * Synopsis:  Load a block again from its checkpoint.
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
 *
 * Throws:    <eslEINVAL> if the block's checkpoint hasn't been saved.
 */
int
p7_chkgmx_LoadBlock(P7_CHKGMX *cx, int b)
{
  size_t nbytes = sizeof(float) * (cx->M + 1) * cx->nscells;
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
 * 3. Unit tests
 *****************************************************************/
#ifdef p7CHKGMX_TESTDRIVE

/* Each row's cells are set to row * 10000 + position in the row. */
static void
write_row(P7_CHKGMX *cx, int i)
{
  int n;
  if (cx->dp[i] == NULL) esl_fatal("row %d isn't loaded", i);
  for (n = 0; n < (cx->M + 1) * cx->nscells; n++) cx->dp[i][n] = (float) (i * 10000 + n);
}

static void
check_row(const P7_CHKGMX *cx, int i)
{
  int n;
  if (cx->dp[i] == NULL) esl_fatal("row %d isn't loaded", i);
  for (n = 0; n < (cx->M + 1) * cx->nscells; n++)
    if (cx->dp[i][n] != (float) (i * 10000 + n)) esl_fatal("row %d holds wrong values", i);
}

/* check_loaded(): the rows of block <b> and its saved rows (before it
 * going forward, after it going backward) are loaded, all others NULL. */
static void
check_loaded(const P7_CHKGMX *cx, int b)
{
  int i, lo, hi;

  if (cx->dir == p7_CHK_FORWARD) { lo = b - cx->saveR; hi = b + cx->blockR - 1;              }
  else                           { lo = b;             hi = b + cx->blockR + cx->saveR - 1;  }
  for (i = 0; i <= cx->L; i++)
    if ((i >= lo && i <= hi) != (cx->dp[i] != NULL)) esl_fatal("row %d pointer wrong for block %d", i, b);
}

/* check_macros(): the {MDI}MX() macros read row <i> through cx->dp. */
static void
check_macros(P7_CHKGMX *cx, int i)
{
  float **dp = cx->dp;
  int     k;

  for (k = 0; k <= cx->M; k++)
    if (MMX(i,k) != (float) (i * 10000 + k * p7G_NSCELLS + p7G_M) ||
        DMX(i,k) != (float) (i * 10000 + k * p7G_NSCELLS + p7G_D) ||
        IMX(i,k) != (float) (i * 10000 + k * p7G_NSCELLS + p7G_I))
      esl_fatal("MDI macros wrong at i=%d k=%d", i, k);
}

/* utest_blocks():
 * Forward direction. Run forward through all blocks of an <M> by <L>
 * comparison, writing each row as it is "computed", then load the
 * blocks again from last to first, rewriting each block, and check that
 * the saved rows before each block and the block's own rows hold the
 * right values. Also checks that the {MDI}MX() macros read through
 * cx->dp.
 */
static void
utest_blocks(P7_CHKGMX *cx, int M, int L)
{
  int b, i;

  cx->M = M;
  cx->L = L;

  for (b = 0; b <= L; b += cx->blockR)
    {
      if (p7_chkgmx_StartBlock(cx, b) != eslOK) esl_fatal("StartBlock(%d) failed", b);
      check_loaded(cx, b);
      for (i = ESL_MAX(0, b - cx->saveR); i < b; i++) check_row(cx, i);
      for (i = b; i <= ESL_MIN(L, b + cx->blockR - 1); i++) write_row(cx, i);
    }
  if (cx->nK != L / cx->blockR) esl_fatal("%d checkpoints saved, expected %d", cx->nK, L / cx->blockR);

  for (b = (L / cx->blockR) * cx->blockR; b >= 0; b -= cx->blockR)
    {
      if (p7_chkgmx_LoadBlock(cx, b) != eslOK) esl_fatal("LoadBlock(%d) failed", b);
      check_loaded(cx, b);
      for (i = ESL_MAX(0, b - cx->saveR); i < b; i++) check_row(cx, i);
      for (i = b; i <= ESL_MIN(L, b + cx->blockR - 1); i++) write_row(cx, i);
      for (i = b; i <= ESL_MIN(L, b + cx->blockR - 1); i++) check_row(cx, i);
      check_macros(cx, b);
    }

  /* A block that doesn't follow the loaded one is refused */
  if (L >= 2 * cx->blockR && p7_chkgmx_StartBlock(cx, 2 * cx->blockR) != eslEINVAL) esl_fatal("out-of-order StartBlock accepted");
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
utest_blocks_backward(P7_CHKGMX *cx, int M, int L)
{
  int top, b, i;

  cx->M = M;
  cx->L = L;
  top   = (L / cx->blockR) * cx->blockR;

  for (b = top; b >= 0; b -= cx->blockR)
    {
      if (p7_chkgmx_StartBlock(cx, b) != eslOK) esl_fatal("backward StartBlock(%d) failed", b);
      check_loaded(cx, b);
      for (i = b + cx->blockR; i <= ESL_MIN(L, b + cx->blockR + cx->saveR - 1); i++) check_row(cx, i);
      for (i = ESL_MIN(L, b + cx->blockR - 1); i >= b; i--) write_row(cx, i);
    }
  if (cx->nK != L / cx->blockR) esl_fatal("%d checkpoints saved, expected %d", cx->nK, L / cx->blockR);

  for (b = 0; b <= top; b += cx->blockR)
    {
      if (p7_chkgmx_LoadBlock(cx, b) != eslOK) esl_fatal("backward LoadBlock(%d) failed", b);
      check_loaded(cx, b);
      for (i = b + cx->blockR; i <= ESL_MIN(L, b + cx->blockR + cx->saveR - 1); i++) check_row(cx, i);
      for (i = ESL_MIN(L, b + cx->blockR - 1); i >= b; i--) write_row(cx, i);
      for (i = b; i <= ESL_MIN(L, b + cx->blockR - 1); i++) check_row(cx, i);
      check_macros(cx, b);
    }

  /* A block that doesn't follow the loaded one is refused */
  if (top >= 2 * cx->blockR)
    {
      if (p7_chkgmx_StartBlock(cx, top) != eslOK) esl_fatal("backward StartBlock(top) failed");
      if (p7_chkgmx_StartBlock(cx, top - 2 * cx->blockR) != eslEINVAL) esl_fatal("out-of-order backward StartBlock accepted");
    }
}

static void
utest_chkgmx(void)
{
  P7_CHKGMX *cx = NULL;

  if ((cx = p7_chkgmx_Create(13, 37, 8, 3, p7G_NSCELLS, p7_CHK_FORWARD)) == NULL) esl_fatal("Create failed");
  utest_blocks(cx, 13, 37);
  utest_blocks(cx, 13, 32);                    /* L a multiple of the block size */
  utest_blocks(cx,  5,  6);                    /* a single block                 */

  if (p7_chkgmx_GrowTo(cx, 70, 101, 8)  != eslOK) esl_fatal("GrowTo failed");
  utest_blocks(cx, 70, 101);
  if (p7_chkgmx_GrowTo(cx, 20, 40, 8)   != eslOK) esl_fatal("GrowTo (smaller) failed");
  utest_blocks(cx, 20, 40);
  if (p7_chkgmx_GrowTo(cx, 20, 140, 13) != eslOK) esl_fatal("GrowTo (bigger block) failed");
  utest_blocks(cx, 20, 140);
  if (p7_chkgmx_GrowTo(cx, 20, 140, 4)  != eslOK) esl_fatal("GrowTo (smaller block) failed");
  utest_blocks(cx, 20, 140);
  if (p7_chkgmx_GrowTo(cx, 20, 140, 2)  != eslEINVAL) esl_fatal("GrowTo with blockR < saveR accepted");
  p7_chkgmx_Destroy(cx);

  if ((cx = p7_chkgmx_Create(7, 50, 8, 8, p7G_NSCELLS, p7_CHK_FORWARD)) == NULL) esl_fatal("Create failed");
  utest_blocks(cx, 7, 50);                     /* saveR == blockR                */
  p7_chkgmx_Destroy(cx);

  if ((cx = p7_chkgmx_Create(13, 40, 8, 5, p7G_NSCELLS_FS, p7_CHK_BACKWARD)) == NULL) esl_fatal("backward Create failed");
  utest_blocks_backward(cx, 13, 37);           /* last block shorter than saveR  */
  utest_blocks_backward(cx, 13, 40);           /* last block is the single row L */
  utest_blocks_backward(cx, 13, 39);           /* last block one row short       */
  utest_blocks_backward(cx,  5,  6);           /* a single block                 */
  if (p7_chkgmx_GrowTo(cx, 30, 200, 32) != eslOK) esl_fatal("backward GrowTo failed");
  utest_blocks_backward(cx, 30, 200);
  if (p7_chkgmx_GrowTo(cx, 30, 200, 5)  != eslOK) esl_fatal("backward GrowTo (block == saveR) failed");
  utest_blocks_backward(cx, 30, 200);

  /* the same allocation turned around: forward, then backward with more saved rows */
  if (p7_chkgmx_SetDirection(cx, p7_CHK_FORWARD, 3)  != eslOK) esl_fatal("SetDirection (forward) failed");
  if (p7_chkgmx_GrowTo(cx, 13, 37, 8)                != eslOK) esl_fatal("GrowTo after SetDirection failed");
  utest_blocks(cx, 13, 37);
  if (p7_chkgmx_SetDirection(cx, p7_CHK_BACKWARD, 7) != eslOK) esl_fatal("SetDirection (backward) failed");
  if (p7_chkgmx_GrowTo(cx, 30, 200, 9)               != eslOK) esl_fatal("GrowTo after SetDirection failed");
  utest_blocks_backward(cx, 30, 200);
  if (p7_chkgmx_SetDirection(cx, 2, 3)               != eslEINVAL) esl_fatal("SetDirection with bad direction accepted");
  if (p7_chkgmx_SetDirection(cx, p7_CHK_FORWARD, 0)  != eslEINVAL) esl_fatal("SetDirection with saveR < 1 accepted");
  p7_chkgmx_Destroy(cx);

  if (p7_chkgmx_Create(7, 50, 4, 5, p7G_NSCELLS, p7_CHK_FORWARD)  != NULL) esl_fatal("saveR > blockR accepted");
  if (p7_chkgmx_Create(7, 50, 8, 5, p7G_NSCELLS, 2)               != NULL) esl_fatal("bad direction accepted");
}

#endif /*p7CHKGMX_TESTDRIVE*/


/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7CHKGMX_TESTDRIVE
/*
   gcc -g -Wall -std=gnu99 -o p7_chkgmx_utest -I. -L. -I../easel -L../easel -Dp7CHKGMX_TESTDRIVE p7_chkgmx.c -lhmmer -leasel -lm
   ./p7_chkgmx_utest
 */
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default   env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the P7_CHKGMX checkpointed matrix";

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);

  esl_exception_SetHandler(&esl_nonfatal_handler);   /* the tests check that bad arguments are refused */
  utest_chkgmx();

  esl_getopts_Destroy(go);
  printf("All tests passed.\n");
  return eslOK;
}
#endif /*p7CHKGMX_TESTDRIVE*/
