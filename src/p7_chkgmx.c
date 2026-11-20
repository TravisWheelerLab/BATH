/* P7_CHKGMX : a checkpointed generic dynamic programming matrix
 *
 * Contents:
 *   1. The <P7_CHKGMX> object
 *   2. Loading blocks and checkpoints
 *   3. Unit tests
 *   4. Test driver
 */
#include "p7_config.h"

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

/* set_block_ptrs(): point the rows of block <b>, and the <saveR> rows
 * before it, at their storage slots. Rows of the previously loaded
 * block are set to NULL first. Uses <cx->L> to bound the last block.
 */
static void
set_block_ptrs(P7_CHKGMX *cx, int b)
{
  int i;

  if (cx->b >= 0)
    for (i = ESL_MAX(0, cx->b - cx->saveR); i <= ESL_MIN(cx->L, cx->b + cx->blockR - 1); i++)
      cx->dp[i] = NULL;

  for (i = ESL_MAX(0, b - cx->saveR); i < b; i++)
    cx->dp[i] = dp_row(cx, i - (b - cx->saveR));
  for (i = b; i <= ESL_MIN(cx->L, b + cx->blockR - 1); i++)
    cx->dp[i] = dp_row(cx, cx->saveR + (i - b));

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
 *            <nscells> cells per model position. Rows are computed in
 *            blocks of <blockR> rows, and each checkpoint saves the
 *            <saveR> rows before a block; <saveR> must be <= <blockR>.
 *
 * Returns:   a pointer to the new <P7_CHKGMX>.
 *
 * Throws:    <NULL> on allocation error or bad arguments.
 */
P7_CHKGMX *
p7_chkgmx_Create(int allocM, int allocL, int blockR, int saveR, int nscells)
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
 *            up to <M> to a sequence of length up to <L>; reallocates
 *            if necessary. The block size, saved rows and cells per
 *            model position stay as they were set by p7_chkgmx_Create().
 *
 *            All row pointers are set to NULL and no block is loaded
 *            afterwards; any stored rows and checkpoints are invalid.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_chkgmx_GrowTo(P7_CHKGMX *cx, int M, int L)
{
  void *p;
  int   W = ESL_MAX(cx->allocW, M + 1);
  int   K = L / cx->blockR;
  int   i;
  int   status;

  if (W > cx->allocW) {
    ESL_RALLOC(cx->dp_mem, p, sizeof(float) * cx->allocR * W * cx->nscells);
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


/*****************************************************************
 *= 2. Loading blocks and checkpoints.
 *****************************************************************/

/* Function:  p7_chkgmx_StartBlock()
 * Synopsis:  Move a forward computation into the next block.
 *
 * Purpose:   Make block <b> (a multiple of <cx->blockR>, <= <cx->L>) the
 *            loaded block, for a computation that runs forward through
 *            the rows. For <b> = 0 this just points rows 0..blockR-1
 *            at their storage. Otherwise block <b> - <blockR> must be
 *            the loaded block: the <saveR> rows before <b> are copied
 *            into checkpoint <b>/<blockR> - 1 and into the slots before
 *            the block, and the rows of block <b> are pointed at their
 *            storage, ready to be computed.
 *
 *            Caller sets <cx->M> and <cx->L> before the first call.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <b> isn't the block after the loaded one.
 */
int
p7_chkgmx_StartBlock(P7_CHKGMX *cx, int b)
{
  size_t nbytes = sizeof(float) * (cx->M + 1) * cx->nscells;
  int    k;
  int    r;

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
  return eslOK;
}


/* Function:  p7_chkgmx_LoadBlock()
 * Synopsis:  Load a block again from its checkpoint.
 *
 * Purpose:   Make block <b> (a multiple of <cx->blockR>, <= <cx->L>) the
 *            loaded block, so it can be computed again. For <b> > 0 the
 *            <saveR> rows of checkpoint <b>/<blockR> - 1 are copied into
 *            the slots before the block. The rows of block <b> are
 *            pointed at their storage; their contents are left over
 *            from earlier blocks until the caller computes them.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if checkpoint <b>/<blockR> - 1 hasn't been saved.
 */
int
p7_chkgmx_LoadBlock(P7_CHKGMX *cx, int b)
{
  size_t nbytes = sizeof(float) * (cx->M + 1) * cx->nscells;
  int    k;
  int    r;

  if (b % cx->blockR != 0 || b > cx->L || b / cx->blockR > cx->nK)
    ESL_EXCEPTION(eslEINVAL, "no checkpoint for block %d", b);

  if (b > 0)
    {
      k = b / cx->blockR - 1;
      for (r = 0; r < cx->saveR; r++)
        memcpy(dp_row(cx, r), ck_row(cx, k, r), nbytes);
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

/* check_loaded(): rows b-saveR .. block end are loaded, all others NULL. */
static void
check_loaded(const P7_CHKGMX *cx, int b)
{
  int i;
  for (i = 0; i <= cx->L; i++)
    if ((i >= b - cx->saveR && i < b + cx->blockR) != (cx->dp[i] != NULL)) esl_fatal("row %d pointer wrong for block %d", i, b);
}

/* utest_blocks():
 * Run forward through all blocks of an <M> by <L> comparison, writing
 * each row as it is "computed", then load the blocks again from last
 * to first, rewriting each block, and check that the saved rows before
 * each block and the block's own rows hold the right values. Also
 * checks that the {MDI}MX() macros read through cx->dp.
 */
static void
utest_blocks(P7_CHKGMX *cx, int M, int L)
{
  float **dp = cx->dp;
  int     b, i, k;

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

      i = b;
      for (k = 0; k <= M; k++)
        if (MMX(i,k) != (float) (i * 10000 + k * p7G_NSCELLS + p7G_M) ||
            DMX(i,k) != (float) (i * 10000 + k * p7G_NSCELLS + p7G_D) ||
            IMX(i,k) != (float) (i * 10000 + k * p7G_NSCELLS + p7G_I))
          esl_fatal("MDI macros wrong at i=%d k=%d", i, k);
    }

  /* A block that doesn't follow the loaded one is refused */
  if (L >= 2 * cx->blockR && p7_chkgmx_StartBlock(cx, 2 * cx->blockR) != eslEINVAL) esl_fatal("out-of-order StartBlock accepted");
}

static void
utest_chkgmx(void)
{
  P7_CHKGMX *cx = NULL;

  if ((cx = p7_chkgmx_Create(13, 37, 8, 3, p7G_NSCELLS)) == NULL) esl_fatal("Create failed");
  utest_blocks(cx, 13, 37);
  utest_blocks(cx, 13, 32);                    /* L a multiple of the block size */
  utest_blocks(cx,  5,  6);                    /* a single block                 */

  if (p7_chkgmx_GrowTo(cx, 70, 101) != eslOK) esl_fatal("GrowTo failed");
  utest_blocks(cx, 70, 101);
  if (p7_chkgmx_GrowTo(cx, 20, 40)  != eslOK) esl_fatal("GrowTo (smaller) failed");
  utest_blocks(cx, 20, 40);
  p7_chkgmx_Destroy(cx);

  if ((cx = p7_chkgmx_Create(7, 50, 8, 8, p7G_NSCELLS)) == NULL) esl_fatal("Create failed");
  utest_blocks(cx, 7, 50);                     /* saveR == blockR                */
  p7_chkgmx_Destroy(cx);

  if (p7_chkgmx_Create(7, 50, 4, 5, p7G_NSCELLS) != NULL) esl_fatal("saveR > blockR accepted");
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
