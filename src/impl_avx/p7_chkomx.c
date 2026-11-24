/* P7_CHKOMX dispatcher.
 * Provides the non-suffixed function pointers declared in impl_avx.h.
 * ISA-specific work is delegated to p7_chkomx_sse.c and p7_chkomx_avx.c.
 *
 * Contents:
 *    1. Dispatcher wrappers.
 *    2. Unit tests.
 *    3. Test driver.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/*****************************************************************
 * 1. Dispatcher wrappers.
 *****************************************************************/

static P7_CHKOMX *p7_chkomx_Create_Dispatcher    (int allocM, int allocL, int blockR, int saveR, int nscells, int dir);
static int        p7_chkomx_GrowTo_Dispatcher    (P7_CHKOMX *cx, int allocM, int allocL, int blockR);
static int        p7_chkomx_StartBlock_Dispatcher(P7_CHKOMX *cx, int b);
static int        p7_chkomx_LoadBlock_Dispatcher (P7_CHKOMX *cx, int b);
static float      p7_chkomx_FGetMDI_Dispatcher   (const P7_CHKOMX *cx, int s, int i, int k);
static void       p7_chkomx_Destroy_Dispatcher   (P7_CHKOMX *cx);

P7_CHKOMX *(*p7_chkomx_Create)    (int allocM, int allocL, int blockR, int saveR, int nscells, int dir) = p7_chkomx_Create_Dispatcher;
int        (*p7_chkomx_GrowTo)    (P7_CHKOMX *cx, int allocM, int allocL, int blockR)                      = p7_chkomx_GrowTo_Dispatcher;
int        (*p7_chkomx_StartBlock)(P7_CHKOMX *cx, int b)                                       = p7_chkomx_StartBlock_Dispatcher;
int        (*p7_chkomx_LoadBlock) (P7_CHKOMX *cx, int b)                                       = p7_chkomx_LoadBlock_Dispatcher;
float      (*p7_chkomx_FGetMDI)   (const P7_CHKOMX *cx, int s, int i, int k)                   = p7_chkomx_FGetMDI_Dispatcher;
void       (*p7_chkomx_Destroy)   (P7_CHKOMX *cx)                                              = p7_chkomx_Destroy_Dispatcher;


static P7_CHKOMX *
p7_chkomx_Create_Dispatcher(int allocM, int allocL, int blockR, int saveR, int nscells, int dir)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_chkomx_Create = p7_chkomx_Create_avx; return p7_chkomx_Create_avx(allocM, allocL, blockR, saveR, nscells, dir); }
#endif
#ifdef eslENABLE_SSE
  p7_chkomx_Create = p7_chkomx_Create_sse;
  return p7_chkomx_Create_sse(allocM, allocL, blockR, saveR, nscells, dir);
#else
  p7_Die("p7_chkomx_Create: no SIMD implementation available");
  return NULL;
#endif
}

static int
p7_chkomx_GrowTo_Dispatcher(P7_CHKOMX *cx, int allocM, int allocL, int blockR)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_chkomx_GrowTo = p7_chkomx_GrowTo_avx; return p7_chkomx_GrowTo_avx(cx, allocM, allocL, blockR); }
#endif
#ifdef eslENABLE_SSE
  p7_chkomx_GrowTo = p7_chkomx_GrowTo_sse;
  return p7_chkomx_GrowTo_sse(cx, allocM, allocL, blockR);
#else
  p7_Die("p7_chkomx_GrowTo: no SIMD implementation available");
  return eslENORESULT;
#endif
}

static int
p7_chkomx_StartBlock_Dispatcher(P7_CHKOMX *cx, int b)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_chkomx_StartBlock = p7_chkomx_StartBlock_avx; return p7_chkomx_StartBlock_avx(cx, b); }
#endif
#ifdef eslENABLE_SSE
  p7_chkomx_StartBlock = p7_chkomx_StartBlock_sse;
  return p7_chkomx_StartBlock_sse(cx, b);
#else
  p7_Die("p7_chkomx_StartBlock: no SIMD implementation available");
  return eslENORESULT;
#endif
}

static int
p7_chkomx_LoadBlock_Dispatcher(P7_CHKOMX *cx, int b)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_chkomx_LoadBlock = p7_chkomx_LoadBlock_avx; return p7_chkomx_LoadBlock_avx(cx, b); }
#endif
#ifdef eslENABLE_SSE
  p7_chkomx_LoadBlock = p7_chkomx_LoadBlock_sse;
  return p7_chkomx_LoadBlock_sse(cx, b);
#else
  p7_Die("p7_chkomx_LoadBlock: no SIMD implementation available");
  return eslENORESULT;
#endif
}

static float
p7_chkomx_FGetMDI_Dispatcher(const P7_CHKOMX *cx, int s, int i, int k)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_chkomx_FGetMDI = p7_chkomx_FGetMDI_avx; return p7_chkomx_FGetMDI_avx(cx, s, i, k); }
#endif
#ifdef eslENABLE_SSE
  p7_chkomx_FGetMDI = p7_chkomx_FGetMDI_sse;
  return p7_chkomx_FGetMDI_sse(cx, s, i, k);
#else
  p7_Die("p7_chkomx_FGetMDI: no SIMD implementation available");
  return 0.0f;
#endif
}

static void
p7_chkomx_Destroy_Dispatcher(P7_CHKOMX *cx)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_chkomx_Destroy = p7_chkomx_Destroy_avx; p7_chkomx_Destroy_avx(cx); return; }
#endif
#ifdef eslENABLE_SSE
  p7_chkomx_Destroy = p7_chkomx_Destroy_sse;
  p7_chkomx_Destroy_sse(cx);
#else
  p7_Die("p7_chkomx_Destroy: no SIMD implementation available");
#endif
}


/*****************************************************************
 * 2. Unit tests.
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
#ifdef eslENABLE_AVX
  if (isa->V == 8) return (float *) cx->dpf_avx[i];
#endif
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
 * 3. Test driver.
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

#ifdef eslENABLE_SSE
  { CHKOMX_ISA sse = { p7_chkomx_Create_sse, p7_chkomx_GrowTo_sse, p7_chkomx_StartBlock_sse, p7_chkomx_LoadBlock_sse, p7_chkomx_FGetMDI_sse, p7_chkomx_Destroy_sse, 4, "SSE" };
    utest_isa(&sse); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    { CHKOMX_ISA avx = { p7_chkomx_Create_avx, p7_chkomx_GrowTo_avx, p7_chkomx_StartBlock_avx, p7_chkomx_LoadBlock_avx, p7_chkomx_FGetMDI_avx, p7_chkomx_Destroy_avx, 8, "AVX" };
      utest_isa(&avx); }
#endif

  esl_getopts_Destroy(go);
  printf("All tests passed.\n");
  return eslOK;
}
#endif /*p7CHKOMX_TESTDRIVE*/
