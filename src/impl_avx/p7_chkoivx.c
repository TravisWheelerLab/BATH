/* p7_chkoivx.c — runtime dispatch for P7_CHKOIVX management.
 *
 * Contains:
 *   1. Self-patching dispatcher function pointer definitions (extern in impl_avx.h).
 *   2. Unit tests.
 *   3. Test driver.
 *
 * ISA-specific implementations live in:
 *   p7_chkoivx_sse.c — SSE (128-bit) versions, suffix _sse
 *   p7_chkoivx_avx.c — AVX-256 versions, suffix _avx
 *
 * Each pointer is initialized to a static dispatcher.  On the first call,
 * the dispatcher detects the ISA, patches the pointer, and calls the
 * implementation.  All subsequent calls skip dispatch entirely.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/*****************************************************************
 * 1. Self-patching dispatcher function pointer definitions.
 *****************************************************************/

static P7_CHKOIVX *p7_chkoivx_Create_Dispatcher (int M_hint, int CF, int CI, int K);
static int         p7_chkoivx_GrowTo_Dispatcher (P7_CHKOIVX *cv, int M, int CF, int CI, int K);
static int         p7_chkoivx_Save_Dispatcher   (P7_CHKOIVX *cv, const P7_OIVX *ov, int M, int k);
static int         p7_chkoivx_Load_Dispatcher   (P7_CHKOIVX *cv, P7_OIVX *ov, int M, int k);
static void        p7_chkoivx_Destroy_Dispatcher(P7_CHKOIVX *cv);

P7_CHKOIVX *(*p7_chkoivx_Create) (int M_hint, int CF, int CI, int K)               = p7_chkoivx_Create_Dispatcher;
int         (*p7_chkoivx_GrowTo) (P7_CHKOIVX *cv, int M, int CF, int CI, int K)    = p7_chkoivx_GrowTo_Dispatcher;
int         (*p7_chkoivx_Save)   (P7_CHKOIVX *cv, const P7_OIVX *ov, int M, int k) = p7_chkoivx_Save_Dispatcher;
int         (*p7_chkoivx_Load)   (P7_CHKOIVX *cv, P7_OIVX *ov, int M, int k)       = p7_chkoivx_Load_Dispatcher;
void        (*p7_chkoivx_Destroy)(P7_CHKOIVX *cv)                                  = p7_chkoivx_Destroy_Dispatcher;

static P7_CHKOIVX *
p7_chkoivx_Create_Dispatcher(int M_hint, int CF, int CI, int K)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_chkoivx_Create = p7_chkoivx_Create_avx; return p7_chkoivx_Create_avx(M_hint, CF, CI, K); }
#endif
#ifdef eslENABLE_SSE
  p7_chkoivx_Create = p7_chkoivx_Create_sse;
  return p7_chkoivx_Create_sse(M_hint, CF, CI, K);
#else
  p7_Die("p7_chkoivx_Create: no SIMD implementation available");
  return NULL;
#endif
}

static int
p7_chkoivx_GrowTo_Dispatcher(P7_CHKOIVX *cv, int M, int CF, int CI, int K)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_chkoivx_GrowTo = p7_chkoivx_GrowTo_avx; return p7_chkoivx_GrowTo_avx(cv, M, CF, CI, K); }
#endif
#ifdef eslENABLE_SSE
  p7_chkoivx_GrowTo = p7_chkoivx_GrowTo_sse;
  return p7_chkoivx_GrowTo_sse(cv, M, CF, CI, K);
#else
  p7_Die("p7_chkoivx_GrowTo: no SIMD implementation available");
  return eslENORESULT;
#endif
}

static int
p7_chkoivx_Save_Dispatcher(P7_CHKOIVX *cv, const P7_OIVX *ov, int M, int k)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_chkoivx_Save = p7_chkoivx_Save_avx; return p7_chkoivx_Save_avx(cv, ov, M, k); }
#endif
#ifdef eslENABLE_SSE
  p7_chkoivx_Save = p7_chkoivx_Save_sse;
  return p7_chkoivx_Save_sse(cv, ov, M, k);
#else
  p7_Die("p7_chkoivx_Save: no SIMD implementation available");
  return eslENORESULT;
#endif
}

static int
p7_chkoivx_Load_Dispatcher(P7_CHKOIVX *cv, P7_OIVX *ov, int M, int k)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_chkoivx_Load = p7_chkoivx_Load_avx; return p7_chkoivx_Load_avx(cv, ov, M, k); }
#endif
#ifdef eslENABLE_SSE
  p7_chkoivx_Load = p7_chkoivx_Load_sse;
  return p7_chkoivx_Load_sse(cv, ov, M, k);
#else
  p7_Die("p7_chkoivx_Load: no SIMD implementation available");
  return eslENORESULT;
#endif
}

static void
p7_chkoivx_Destroy_Dispatcher(P7_CHKOIVX *cv)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_chkoivx_Destroy = p7_chkoivx_Destroy_avx; p7_chkoivx_Destroy_avx(cv); return; }
#endif
#ifdef eslENABLE_SSE
  p7_chkoivx_Destroy = p7_chkoivx_Destroy_sse;
  p7_chkoivx_Destroy_sse(cv);
#else
  p7_Die("p7_chkoivx_Destroy: no SIMD implementation available");
#endif
}


/*****************************************************************
 * 2. Unit tests.
 *****************************************************************/
#ifdef p7CHKOIVX_TESTDRIVE
#include "p7_splice.h"

typedef struct {
  P7_CHKOIVX *(*Create) (int, int, int, int);
  int         (*GrowTo) (P7_CHKOIVX *, int, int, int, int);
  int         (*Save)   (P7_CHKOIVX *, const P7_OIVX *, int, int);
  int         (*Load)   (P7_CHKOIVX *, P7_OIVX *, int, int);
  void        (*Destroy)(P7_CHKOIVX *);
  P7_OIVX    *(*OCreate)(int, int);
  void        (*ODestroy)(P7_OIVX *);
  int           V;     /* floats (or int32s) per vector */
  const char   *name;
} CHKOIVX_ISA;

static float *
float_channel(const P7_OIVX *ov, const CHKOIVX_ISA *isa, int c)
{
#ifdef eslENABLE_AVX
  if (isa->V == 8) return (float *) ov->ivx_avx[c];
#endif
  return (float *) ov->ivx[c];
}

static int32_t *
int_channel(const P7_CHKOIVX *cv, const CHKOIVX_ISA *isa, int c)
{
#ifdef eslENABLE_AVX
  if (isa->V == 8) return (int32_t *) cv->ivx_avx[c];
#endif
  return (int32_t *) cv->ivx[c];
}

/* Fill the working channels with values tagged by <tag>: floats
 * tag*100000 + c*1000 + n, ints -(tag*100000 + c*1000 + n).
 */
static void
write_channels(P7_OIVX *ov, P7_CHKOIVX *cv, const CHKOIVX_ISA *isa, int CF, int CI, int n, int tag)
{
  int c, i;
  for (c = 0; c < CF; c++) { float   *p = float_channel(ov, isa, c); for (i = 0; i < n; i++) p[i] =  (float) (tag * 100000 + c * 1000 + i); }
  for (c = 0; c < CI; c++) { int32_t *p = int_channel  (cv, isa, c); for (i = 0; i < n; i++) p[i] = -(tag * 100000 + c * 1000 + i); }
}

static void
check_channels(const P7_OIVX *ov, const P7_CHKOIVX *cv, const CHKOIVX_ISA *isa, int CF, int CI, int n, int tag)
{
  int c, i;
  for (c = 0; c < CF; c++) { float   *p = float_channel(ov, isa, c); for (i = 0; i < n; i++) if (p[i] !=  (float) (tag * 100000 + c * 1000 + i)) esl_fatal("%s: float channel %d wrong for checkpoint %d", isa->name, c, tag); }
  for (c = 0; c < CI; c++) { int32_t *p = int_channel  (cv, isa, c); for (i = 0; i < n; i++) if (p[i] != -(tag * 100000 + c * 1000 + i)) esl_fatal("%s: int channel %d wrong for checkpoint %d",   isa->name, c, tag); }
}

/* utest_saveload():
 * Save distinct values into each of <K> checkpoints, overwrite the
 * working channels, then load every checkpoint back (last to first)
 * and check its values.
 */
static void
utest_saveload(P7_OIVX *ov, P7_CHKOIVX *cv, const CHKOIVX_ISA *isa, int M, int CF, int CI, int K)
{
  int n = ESL_MAX(2, (M - 1) / isa->V + 1) * isa->V;
  int k;

  for (k = 0; k < K; k++)
    {
      write_channels(ov, cv, isa, CF, CI, n, k + 1);
      if (isa->Save(cv, ov, M, k) != eslOK) esl_fatal("%s: Save(%d) failed", isa->name, k);
    }
  write_channels(ov, cv, isa, CF, CI, n, 99);
  for (k = K - 1; k >= 0; k--)
    {
      if (isa->Load(cv, ov, M, k) != eslOK) esl_fatal("%s: Load(%d) failed", isa->name, k);
      check_channels(ov, cv, isa, CF, CI, n, k + 1);
    }

  if (isa->Save(cv, ov, M, K) != eslEINVAL) esl_fatal("%s: Save past the last checkpoint accepted", isa->name);
}

static void
utest_isa(const CHKOIVX_ISA *isa)
{
  P7_OIVX    *ov = NULL;
  P7_CHKOIVX *cv = NULL;

  if ((ov = isa->OCreate(13, SIGNAL_MEM_SIZE))                       == NULL) esl_fatal("%s: P7_OIVX Create failed", isa->name);
  if ((cv = isa->Create(13, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, 4))    == NULL) esl_fatal("%s: Create failed", isa->name);
  utest_saveload(ov, cv, isa, 13, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, 4);
  utest_saveload(ov, cv, isa,  5, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, 4);   /* smaller M, same allocation */
  isa->ODestroy(ov);

  if ((ov = isa->OCreate(70, SIGNAL_MEM_SIZE))                                    == NULL) esl_fatal("%s: P7_OIVX Create failed", isa->name);
  if (isa->GrowTo(cv, 70, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, 9)                    != eslOK) esl_fatal("%s: GrowTo failed", isa->name);
  utest_saveload(ov, cv, isa, 70, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, 9);
  isa->ODestroy(ov);
  isa->Destroy(cv);
}

#endif /*p7CHKOIVX_TESTDRIVE*/


/*****************************************************************
 * 3. Test driver.
 *****************************************************************/
#ifdef p7CHKOIVX_TESTDRIVE
/*
   make chkoivx_utest
   ./chkoivx_utest
 */
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default   env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for P7_CHKOIVX checkpoint vectors";

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);

  esl_exception_SetHandler(&esl_nonfatal_handler);   /* the tests check that bad arguments are refused */

#ifdef eslENABLE_SSE
  { CHKOIVX_ISA sse = { p7_chkoivx_Create_sse, p7_chkoivx_GrowTo_sse, p7_chkoivx_Save_sse, p7_chkoivx_Load_sse, p7_chkoivx_Destroy_sse, p7_oivx_Create_sse, p7_oivx_Destroy_sse, 4, "SSE" };
    utest_isa(&sse); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    { CHKOIVX_ISA avx = { p7_chkoivx_Create_avx, p7_chkoivx_GrowTo_avx, p7_chkoivx_Save_avx, p7_chkoivx_Load_avx, p7_chkoivx_Destroy_avx, p7_oivx_Create_avx, p7_oivx_Destroy_avx, 8, "AVX" };
      utest_isa(&avx); }
#endif

  esl_getopts_Destroy(go);
  printf("All tests passed.\n");
  return eslOK;
}
#endif /*p7CHKOIVX_TESTDRIVE*/
