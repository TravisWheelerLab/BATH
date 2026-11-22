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

static P7_CHKOMX *p7_chkomx_Create_Dispatcher    (int allocM, int allocL, int blockR, int saveR, int nscells);
static int        p7_chkomx_GrowTo_Dispatcher    (P7_CHKOMX *cx, int allocM, int allocL);
static int        p7_chkomx_StartBlock_Dispatcher(P7_CHKOMX *cx, int b);
static int        p7_chkomx_LoadBlock_Dispatcher (P7_CHKOMX *cx, int b);
static float      p7_chkomx_FGetMDI_Dispatcher   (const P7_CHKOMX *cx, int s, int i, int k);
static void       p7_chkomx_Destroy_Dispatcher   (P7_CHKOMX *cx);

P7_CHKOMX *(*p7_chkomx_Create)    (int allocM, int allocL, int blockR, int saveR, int nscells) = p7_chkomx_Create_Dispatcher;
int        (*p7_chkomx_GrowTo)    (P7_CHKOMX *cx, int allocM, int allocL)                      = p7_chkomx_GrowTo_Dispatcher;
int        (*p7_chkomx_StartBlock)(P7_CHKOMX *cx, int b)                                       = p7_chkomx_StartBlock_Dispatcher;
int        (*p7_chkomx_LoadBlock) (P7_CHKOMX *cx, int b)                                       = p7_chkomx_LoadBlock_Dispatcher;
float      (*p7_chkomx_FGetMDI)   (const P7_CHKOMX *cx, int s, int i, int k)                   = p7_chkomx_FGetMDI_Dispatcher;
void       (*p7_chkomx_Destroy)   (P7_CHKOMX *cx)                                              = p7_chkomx_Destroy_Dispatcher;


static P7_CHKOMX *
p7_chkomx_Create_Dispatcher(int allocM, int allocL, int blockR, int saveR, int nscells)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_chkomx_Create = p7_chkomx_Create_avx; return p7_chkomx_Create_avx(allocM, allocL, blockR, saveR, nscells); }
#endif
#ifdef eslENABLE_SSE
  p7_chkomx_Create = p7_chkomx_Create_sse;
  return p7_chkomx_Create_sse(allocM, allocL, blockR, saveR, nscells);
#else
  p7_Die("p7_chkomx_Create: no SIMD implementation available");
  return NULL;
#endif
}

static int
p7_chkomx_GrowTo_Dispatcher(P7_CHKOMX *cx, int allocM, int allocL)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_chkomx_GrowTo = p7_chkomx_GrowTo_avx; return p7_chkomx_GrowTo_avx(cx, allocM, allocL); }
#endif
#ifdef eslENABLE_SSE
  p7_chkomx_GrowTo = p7_chkomx_GrowTo_sse;
  return p7_chkomx_GrowTo_sse(cx, allocM, allocL);
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
  P7_CHKOMX *(*Create)    (int, int, int, int, int);
  int        (*GrowTo)    (P7_CHKOMX *, int, int);
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

/* check_loaded(): rows b-saveR .. block end are loaded, all others NULL. */
static void
check_loaded(const P7_CHKOMX *cx, const CHKOMX_ISA *isa, int b)
{
  int i;

  for (i = 0; i <= cx->L; i++)
    {
      int loaded = (i >= b - cx->saveR && i < b + cx->blockR);
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

static void
utest_isa(const CHKOMX_ISA *isa)
{
  P7_CHKOMX *cx = NULL;

  if ((cx = isa->Create(13, 37, 8, 3, p7X_NSCELLS)) == NULL) esl_fatal("%s: Create failed", isa->name);
  utest_blocks(cx, isa, 13, 37);
  utest_blocks(cx, isa, 13, 32);                   /* L a multiple of the block size */
  utest_blocks(cx, isa,  5,  6);                   /* a single block                 */

  if (isa->GrowTo(cx, 70, 101) != eslOK) esl_fatal("%s: GrowTo failed", isa->name);
  utest_blocks(cx, isa, 70, 101);
  if (isa->GrowTo(cx, 20, 40) != eslOK) esl_fatal("%s: GrowTo (smaller) failed", isa->name);
  utest_blocks(cx, isa, 20, 40);
  isa->Destroy(cx);

  if ((cx = isa->Create(7, 50, 8, 8, p7X_NSCELLS)) == NULL) esl_fatal("%s: Create failed", isa->name);
  utest_blocks(cx, isa, 7, 50);                    /* saveR == blockR                */
  isa->Destroy(cx);

  if (isa->Create(7, 50, 4, 5, p7X_NSCELLS) != NULL) esl_fatal("%s: saveR > blockR accepted", isa->name);
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
