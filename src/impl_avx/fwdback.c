/* Forward/Backward dispatcher.
 * Provides the non-suffixed p7_Forward, p7_ForwardParser, p7_Backward,
 * p7_BackwardParser, p7_BackwardChk, p7_BackwardChkLoad, p7_ForwardChk
 * and p7_ForwardChkLoad declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"
#include <string.h>

/* Forward declaration of dispatcher */
static int p7_Forward_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *fwd, float *opt_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Forward)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *fwd, float *opt_sc) = p7_Forward_Dispatcher;

/* Function:  p7_Forward_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for Forward algorithm.
 *
 * Returns:   <eslOK> on success.  <*opt_sc> is the log Forward score in nats.
 * Throws:    <eslEINVAL> if <ox> allocation is too small, or profile is
 *            not in a local alignment mode.
 *            <eslEMEM> on reallocation failure.
 */
static int
p7_Forward_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *fwd, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_Forward = p7_Forward_avx512; return p7_Forward_avx512(dsq, L, om, fwd, opt_sc); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Forward = p7_Forward_avx;    return p7_Forward_avx(dsq, L, om, fwd, opt_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_Forward = p7_Forward_sse;
  return p7_Forward_sse(dsq, L, om, fwd, opt_sc);
#else
  p7_Die("p7_Forward: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_ForwardParser_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *fwd, float *opt_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_ForwardParser)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *fwd, float *opt_sc) = p7_ForwardParser_Dispatcher;

/* Function:  p7_ForwardParser_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for Forward parser (linear memory).
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
static int
p7_ForwardParser_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *fwd, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_ForwardParser = p7_ForwardParser_avx512; return p7_ForwardParser_avx512(dsq, L, om, fwd, opt_sc); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_ForwardParser = p7_ForwardParser_avx;    return p7_ForwardParser_avx(dsq, L, om, fwd, opt_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_ForwardParser = p7_ForwardParser_sse;
  return p7_ForwardParser_sse(dsq, L, om, fwd, opt_sc);
#else
  p7_Die("p7_ForwardParser: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_Backward_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Backward)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc) = p7_Backward_Dispatcher;

/* Function:  p7_Backward_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for Backward algorithm.
 *
 * Returns:   <eslOK> on success.  <*opt_sc> is the log Backward score.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
static int
p7_Backward_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_Backward = p7_Backward_avx512; return p7_Backward_avx512(dsq, L, om, fwd, bck, opt_sc); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Backward = p7_Backward_avx;    return p7_Backward_avx(dsq, L, om, fwd, bck, opt_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_Backward = p7_Backward_sse;
  return p7_Backward_sse(dsq, L, om, fwd, bck, opt_sc);
#else
  p7_Die("p7_Backward: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_BackwardParser_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_BackwardParser)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc) = p7_BackwardParser_Dispatcher;

/* Function:  p7_BackwardParser_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for Backward parser (linear memory).
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
static int
p7_BackwardParser_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_BackwardParser = p7_BackwardParser_avx512; return p7_BackwardParser_avx512(dsq, L, om, fwd, bck, opt_sc); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_BackwardParser = p7_BackwardParser_avx;    return p7_BackwardParser_avx(dsq, L, om, fwd, bck, opt_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_BackwardParser = p7_BackwardParser_sse;
  return p7_BackwardParser_sse(dsq, L, om, fwd, bck, opt_sc);
#else
  p7_Die("p7_BackwardParser: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_BackwardChk_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_BackwardChk)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc) = p7_BackwardChk_Dispatcher;

/* Function:  p7_BackwardChk_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for the checkpointed Backward with its own scale factors.
 */
static int
p7_BackwardChk_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_BackwardChk = p7_BackwardChk_avx; return p7_BackwardChk_avx(dsq, L, om, cx, opt_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_BackwardChk = p7_BackwardChk_sse;
  return p7_BackwardChk_sse(dsq, L, om, cx, opt_sc);
#else
  p7_Die("p7_BackwardChk: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_BackwardChkLoad_Dispatcher(const ESL_DSQ *dsq, const P7_OPROFILE *om, P7_CHKOMX *cx, int b);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_BackwardChkLoad)(const ESL_DSQ *dsq, const P7_OPROFILE *om, P7_CHKOMX *cx, int b) = p7_BackwardChkLoad_Dispatcher;

/* Function:  p7_BackwardChkLoad_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for loading one block of the checkpointed Backward again.
 */
static int
p7_BackwardChkLoad_Dispatcher(const ESL_DSQ *dsq, const P7_OPROFILE *om, P7_CHKOMX *cx, int b)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_BackwardChkLoad = p7_BackwardChkLoad_avx; return p7_BackwardChkLoad_avx(dsq, om, cx, b); }
#endif
#ifdef eslENABLE_SSE
  p7_BackwardChkLoad = p7_BackwardChkLoad_sse;
  return p7_BackwardChkLoad_sse(dsq, om, cx, b);
#else
  p7_Die("p7_BackwardChkLoad: no SIMD implementation available");
  return eslENORESULT;
#endif
}
/* Forward declarations of dispatchers */
static int p7_ForwardChk_Dispatcher    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc);
static int p7_ForwardChkLoad_Dispatcher(const ESL_DSQ *dsq,        const P7_OPROFILE *om, P7_CHKOMX *cx, int b);

/* Global function pointers, initially pointing at the dispatchers */
int (*p7_ForwardChk)    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc) = p7_ForwardChk_Dispatcher;
int (*p7_ForwardChkLoad)(const ESL_DSQ *dsq,        const P7_OPROFILE *om, P7_CHKOMX *cx, int b)     = p7_ForwardChkLoad_Dispatcher;

/* Function:  p7_ForwardChk_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for the checkpointed Forward.
 *            See p7_ForwardChk_avx().
 */
static int
p7_ForwardChk_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_ForwardChk = p7_ForwardChk_avx; return p7_ForwardChk_avx(dsq, L, om, cx, opt_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_ForwardChk = p7_ForwardChk_sse;
  return p7_ForwardChk_sse(dsq, L, om, cx, opt_sc);
#else
  p7_Die("p7_ForwardChk: no SIMD implementation available");
  return eslENORESULT;
#endif
}

/* Function:  p7_ForwardChkLoad_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for computing one block of a
 *            checkpointed Forward again. See p7_ForwardChkLoad_avx().
 */
static int
p7_ForwardChkLoad_Dispatcher(const ESL_DSQ *dsq, const P7_OPROFILE *om, P7_CHKOMX *cx, int b)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_ForwardChkLoad = p7_ForwardChkLoad_avx; return p7_ForwardChkLoad_avx(dsq, om, cx, b); }
#endif
#ifdef eslENABLE_SSE
  p7_ForwardChkLoad = p7_ForwardChkLoad_sse;
  return p7_ForwardChkLoad_sse(dsq, om, cx, b);
#else
  p7_Die("p7_ForwardChkLoad: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/*****************************************************************
 * 4. Benchmark driver.
 *****************************************************************/
#ifdef p7FWDBACK_BENCHMARK
/* -c, -x options are for debugging and testing: see msvfilter.c for explanation */
/*
   gcc -g -O3 -mavx2 -std=gnu99 -o benchmark-fwdback -I.. -L.. -I../../easel -L../../easel -Dp7FWDBACK_BENCHMARK fwdback.c -lhmmer -leasel -lm
   icc  -O3 -static -o benchmark-fwdback -I.. -L.. -I../../easel -L../../easel -Dp7FWDBACK_BENCHMARK fwdback.c -lhmmer -leasel -lm

   ./benchmark-fwdback <hmmfile>           runs benchmark on both Forward and Backward parser
   ./benchmark-fwdback -c -N100 <hmmfile>  compare scores of AVX to generic impl
   ./benchmark-fwdback -x -N100 <hmmfile>  test that scores match trusted implementation.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_avx.h"
#include <string.h>

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-x", "compare scores to generic implementation (debug)", 0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-c", "equate scores to trusted implementation (debug)",  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-B", "only benchmark Forward",                           0 },
  { "-B",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-F", "only benchmark Backward",                          0 },
  { "-P",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "benchmark parsing version, not full version",      0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for Forward, Backward implementations";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_GMX         *gx      = NULL;
  P7_OMX         *fwd     = NULL;
  P7_OMX         *bck     = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           fsc, bsc;
  float           fsc2, bsc2;
  double          base_time, bench_time, Mcs;

  p7_FLogsumInit();

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);

  if (esl_opt_GetBoolean(go, "-x") && p7_FLogsumError(-0.4, -0.5) > 0.0001)
    p7_Fail("-x here requires p7_Logsum() recompiled in slow exact mode");

  if (esl_opt_GetBoolean(go, "-P")) {
    fwd = p7_omx_Create(gm->M, 0, L);
    bck = p7_omx_Create(gm->M, 0, L);
  } else {
    fwd = p7_omx_Create(gm->M, L, L);
    bck = p7_omx_Create(gm->M, L, L);
  }
  gx  = p7_gmx_Create(gm->M, L, L, p7G_NSCELLS);

  /* Get a baseline time: how long it takes just to generate the sequences */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      if (esl_opt_GetBoolean(go, "-P")) {
	if (! esl_opt_GetBoolean(go, "-B"))  p7_ForwardParser (dsq, L, om,      fwd, &fsc);
	if (! esl_opt_GetBoolean(go, "-F"))  p7_BackwardParser(dsq, L, om, fwd, bck, &bsc);
      } else {
	if (! esl_opt_GetBoolean(go, "-B"))  p7_Forward (dsq, L, om,      fwd, &fsc);
	if (! esl_opt_GetBoolean(go, "-F"))  p7_Backward(dsq, L, om, fwd, bck, &bsc);
      }

      if (esl_opt_GetBoolean(go, "-c") || esl_opt_GetBoolean(go, "-x"))
	{
	  p7_GForward (dsq, L, gm, gx, &fsc2);
	  p7_GBackward(dsq, L, gm, gx, &bsc2);
	  printf("%.4f %.4f %.4f %.4f\n", fsc, bsc, fsc2, bsc2);
	}
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_omx_Destroy(bck);
  p7_omx_Destroy(fwd);
  p7_gmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7FWDBACK_BENCHMARK*/
/*------------------- end, benchmark driver ---------------------*/


/*****************************************************************
 * 5. Unit tests.
 *****************************************************************/
#ifdef p7FWDBACK_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"

/* 
 * compare to GForward() scores.
 */
static void
utest_fwdback(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char        *msg = "forward/backward unit test failed";
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *fwd = p7_omx_Create(M, 0, L);
  P7_OMX      *bck = p7_omx_Create(M, 0, L);
  P7_OMX      *oxf = p7_omx_Create(M, L, L);
  P7_OMX      *oxb = p7_omx_Create(M, L, L);
  P7_GMX      *gx  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  float tolerance;
  float fsc1, fsc2;
  float bsc1, bsc2;
  float generic_sc;

  p7_FLogsumInit();
  if (p7_FLogsumError(-0.4, -0.5) > 0.0001) tolerance = 1.0;  /* weaker test against GForward()   */
  else tolerance = 0.0001;   /* stronger test: FLogsum() is in slow exact mode. */

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_Forward       (dsq, L, om, oxf,      &fsc1);
      p7_Backward      (dsq, L, om, oxf, oxb, &bsc1);
      p7_ForwardParser (dsq, L, om, fwd,      &fsc2);
      p7_BackwardParser(dsq, L, om, fwd, bck, &bsc2);
      p7_GForward      (dsq, L, gm, gx,  &generic_sc);

      /* Forward and Backward scores should agree with high tolerance */
      if (fabs(fsc1-bsc1) > 0.0001)    esl_fatal(msg);
      if (fabs(fsc2-bsc2) > 0.0001)    esl_fatal(msg);
      if (fabs(fsc1-fsc2) > 0.0001)    esl_fatal(msg);

      /* GForward scores should approximate Forward scores, 
       * with tolerance that depends on how logsum.c was compiled
       */
      if (fabs(fsc1-generic_sc) > tolerance) esl_fatal(msg);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(oxb);
  p7_omx_Destroy(oxf);
  p7_omx_Destroy(bck);
  p7_omx_Destroy(fwd);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}

/* utest_backward_chk()
 * Checkpointed Backward with its own scale factors, and checkpointed
 * Forward, on random sequences
 * and the model's consensus between random flanks, up to length <L>:
 *   - a fill in one block and a fill in blocks of 1..8 rows give
 *     identical scores, special states and scale factors;
 *   - every block loaded again, in increasing order as a forward
 *     sweep reads them, is identical to the one-block fill;
 *   - scores and cells (with each matrix's own row scale factors
 *     taken out) agree with p7_Backward() up to rounding.
 *   - the same matrix, turned around by p7_chkomx_SetDirection(),
 *     filled by the checkpointed Forward and each block loaded again
 *     from the last to the first, is identical to p7_Forward().
 * With <require_scaling>, fail if no row was rescaled. With <use_sse>,
 * the _sse versions are called directly, so an AVX CPU also tests the
 * SSE engine.
 */
static void
utest_backward_chk(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N, int require_scaling, int use_sse)
{
  char        *msg   = "checkpointed backward unit test failed";
  P7_HMM      *hmm   = NULL;
  P7_PROFILE  *gm    = NULL;
  P7_OPROFILE *om    = NULL;
  ESL_SQ      *sq    = esl_sq_CreateDigital(abc);
  ESL_DSQ     *dsq   = malloc(sizeof(ESL_DSQ) * (L+2));
  double      *lsf1  = malloc(sizeof(double) * (L+2));
  double      *lsf2  = malloc(sizeof(double) * (L+2));
  P7_OMX      *oxf   = NULL;
  P7_OMX      *oxb   = NULL;
  P7_CHKOMX   *cx1   = NULL;
  P7_CHKOMX   *cx2   = NULL;
  int        (*chk_growto)(P7_CHKOMX *, int, int, int)             = use_sse ? p7_chkomx_GrowTo_sse    : p7_chkomx_GrowTo;
  float      (*chk_get)   (const P7_CHKOMX *, int, int, int)       = use_sse ? p7_chkomx_FGetMDI_sse   : p7_chkomx_FGetMDI;
  float      (*omx_get)   (const P7_OMX *, int, int, int)          = use_sse ? p7_omx_FGetMDI_sse      : p7_omx_FGetMDI;
  int        (*bck_chk)   (const ESL_DSQ *, int, const P7_OPROFILE *, P7_CHKOMX *, float *) = use_sse ? p7_BackwardChk_sse     : p7_BackwardChk;
  int        (*bck_load)  (const ESL_DSQ *, const P7_OPROFILE *, P7_CHKOMX *, int)          = use_sse ? p7_BackwardChkLoad_sse : p7_BackwardChkLoad;
  int        (*fwd_chk)   (const ESL_DSQ *, int, const P7_OPROFILE *, P7_CHKOMX *, float *) = use_sse ? p7_ForwardChk_sse      : p7_ForwardChk;
  int        (*fwd_load)  (const ESL_DSQ *, const P7_OPROFILE *, P7_CHKOMX *, int)          = use_sse ? p7_ForwardChkLoad_sse  : p7_ForwardChkLoad;
  float        fsc, bsc, sc1, sc2;
  float        a, c;
  double       d;
  int          nscaled = 0;
  int          n, Ln, blockR, b, i, k, s, x;

  /* a sampled model with sharpened match emissions and M->M transitions,
   * so its consensus scores high */
  if (p7_hmm_Sample(r, M, abc, &hmm) != eslOK) esl_fatal(msg);
  for (k = 1; k <= M; k++)
    {
      x = esl_vec_FArgMax(hmm->mat[k], abc->K);
      esl_vec_FSet(hmm->mat[k], abc->K, 0.15 / (abc->K - 1));
      hmm->mat[k][x] = 0.85;
      if (k < M) { hmm->t[k][p7H_MM] = 0.9; hmm->t[k][p7H_MI] = 0.05; hmm->t[k][p7H_MD] = 0.05; }
    }
  gm = p7_profile_Create(M, abc);
  if (p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL) != eslOK) esl_fatal(msg);
  if (use_sse)
    {
      om  = p7_oprofile_Create_sse(M, abc);
      p7_oprofile_Convert_sse(gm, om);
      oxf = p7_omx_Create_sse(M, L, L);
      oxb = p7_omx_Create_sse(M, L, L);
      cx1 = p7_chkomx_Create_sse(M, L, L+1, 1, p7X_NSCELLS, p7_CHK_BACKWARD);
      cx2 = p7_chkomx_Create_sse(M, L, 1,   1, p7X_NSCELLS, p7_CHK_BACKWARD);
    }
  else
    {
      om  = p7_oprofile_Create(M, abc);
      p7_oprofile_Convert(gm, om);
      oxf = p7_omx_Create(M, L, L);
      oxb = p7_omx_Create(M, L, L);
      cx1 = p7_chkomx_Create(M, L, L+1, 1, p7X_NSCELLS, p7_CHK_BACKWARD);
      cx2 = p7_chkomx_Create(M, L, 1,   1, p7X_NSCELLS, p7_CHK_BACKWARD);
    }

  for (n = 0; n < N; n++)
    {
      /* alternate random sequences and the model's consensus between
       * random flanks, which scores high enough for rows to be rescaled */
      Ln = 1 + esl_rnd_Roll(r, L);
      esl_rsq_xfIID(r, bg->f, abc->K, Ln, dsq);
      if (n % 2 == 1 && M <= L)
        {
          esl_sq_Reuse(sq);
          p7_emit_SimpleConsensus(hmm, sq);
          Ln = ESL_MAX(Ln, M);
          esl_rsq_xfIID(r, bg->f, abc->K, Ln, dsq);
          memcpy(dsq + 1 + esl_rnd_Roll(r, Ln - M + 1), sq->dsq + 1, sizeof(ESL_DSQ) * M);
        }

      if (use_sse)
        {
          if (p7_Forward_sse (dsq, Ln, om, oxf,      &fsc) != eslOK) esl_fatal(msg);
          if (p7_Backward_sse(dsq, Ln, om, oxf, oxb, &bsc) != eslOK) esl_fatal(msg);
        }
      else
        {
          if (p7_Forward (dsq, Ln, om, oxf,      &fsc) != eslOK) esl_fatal(msg);
          if (p7_Backward(dsq, Ln, om, oxf, oxb, &bsc) != eslOK) esl_fatal(msg);
        }

      blockR = 1 + esl_rnd_Roll(r, 8);
      if (chk_growto(cx1, M, Ln, Ln+1)  != eslOK) esl_fatal(msg);
      if (chk_growto(cx2, M, Ln, blockR) != eslOK) esl_fatal(msg);
      if (bck_chk(dsq, Ln, om, cx1, &sc1) != eslOK) esl_fatal(msg);
      if (bck_chk(dsq, Ln, om, cx2, &sc2) != eslOK) esl_fatal(msg);

      if (sc1 != sc2)                esl_fatal("%s: one-block and blocked scores differ (%f %f)", msg, sc1, sc2);
      if (fabs(sc1 - bsc) > 0.0001 + 1e-5 * fabs(bsc)) esl_fatal("%s: score %f, p7_Backward() %f", msg, sc1, bsc);
      for (i = 0; i <= Ln; i++)
        for (x = 0; x < p7X_NXCELLS; x++)
          if (cx1->xmx[i*p7X_NXCELLS+x] != cx2->xmx[i*p7X_NXCELLS+x]) esl_fatal("%s: special states differ", msg);

      for (b = 0; b <= Ln; b += blockR)
        {
          if (bck_load(dsq, om, cx2, b) != eslOK) esl_fatal(msg);
          for (i = b; i <= ESL_MIN(Ln, b + blockR - 1); i++)
            for (k = 1; k <= M; k++)
              for (s = 0; s < p7X_NSCELLS; s++)
                if (chk_get(cx2, s, i, k) != chk_get(cx1, s, i, k))
                  esl_fatal("%s: loaded block %d differs at row %d, k %d (%g %g)", msg, b, i, k, chk_get(cx2, s, i, k), chk_get(cx1, s, i, k));
        }
      for (i = 0; i <= Ln; i++)
        for (x = 0; x < p7X_NXCELLS; x++)
          if (cx1->xmx[i*p7X_NXCELLS+x] != cx2->xmx[i*p7X_NXCELLS+x]) esl_fatal("%s: special states changed by loading", msg);

      /* the same matrix turned around for the checkpointed Forward, whose
       * fill and every block loaded again, from the last down to the first,
       * are identical to p7_Forward() */
      if (p7_chkomx_SetDirection(cx2, p7_CHK_FORWARD, 1) != eslOK) esl_fatal(msg);
      if (chk_growto(cx2, M, Ln, blockR)                 != eslOK) esl_fatal(msg);
      if (fwd_chk(dsq, Ln, om, cx2, &sc2)                != eslOK) esl_fatal(msg);
      if (sc2 != fsc) esl_fatal("%s: checkpointed Forward score %f, p7_Forward() %f", msg, sc2, fsc);
      for (i = 0; i <= Ln; i++)
        for (x = 0; x < p7X_NXCELLS; x++)
          if (cx2->xmx[i*p7X_NXCELLS+x] != oxf->xmx[i*p7X_NXCELLS+x]) esl_fatal("%s: checkpointed Forward special states differ", msg);
      for (b = (Ln / blockR) * blockR; b >= 0; b -= blockR)
        {
          if (fwd_load(dsq, om, cx2, b) != eslOK) esl_fatal(msg);
          for (i = b; i <= ESL_MIN(Ln, b + blockR - 1); i++)
            for (k = 1; k <= M; k++)
              for (s = 0; s < p7X_NSCELLS; s++)
                if (chk_get(cx2, s, i, k) != omx_get(oxf, s, i, k))
                  esl_fatal("%s: checkpointed Forward block %d differs at row %d, k %d", msg, b, i, k);
        }
      if (p7_chkomx_SetDirection(cx2, p7_CHK_BACKWARD, 1) != eslOK) esl_fatal(msg);

      /* log of the product of each matrix's scale factors for rows i..L */
      lsf1[Ln+1] = lsf2[Ln+1] = 0.;
      for (i = Ln; i >= 1; i--)
        {
          lsf1[i] = lsf1[i+1] + log(oxb->xmx[i*p7X_NXCELLS+p7X_SCALE]);
          lsf2[i] = lsf2[i+1] + log(cx1->xmx[i*p7X_NXCELLS+p7X_SCALE]);
          if (cx1->xmx[i*p7X_NXCELLS+p7X_SCALE] > 1.0f) nscaled++;
        }
      /* Cells below 1e-25 are skipped: terms that fall below the smallest
       * normal float are lost, which can be most of such a cell's value;
       * above 1e-25 the lost terms are negligible. */
      for (i = 1; i <= Ln; i++)
        for (k = 1; k <= M; k++)
          for (s = 0; s < p7X_NSCELLS; s++)
            {
              a = omx_get(oxb, s, i, k);
              c = chk_get(cx1, s, i, k);
              if (a < 1e-25 || c < 1e-25) continue;
              d = (log(a) + lsf1[i]) - (log(c) + lsf2[i]);
              if (fabs(d) > 0.001) esl_fatal("%s: row %d, k %d, state %d: log values differ by %g", msg, i, k, s, d);
            }
    }
  if (require_scaling && nscaled == 0) esl_fatal("%s: no row was rescaled; test is too weak", msg);

  free(dsq);
  free(lsf1);
  free(lsf2);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  if (use_sse)
    {
      p7_omx_Destroy_sse(oxf);
      p7_omx_Destroy_sse(oxb);
      p7_chkomx_Destroy_sse(cx1);
      p7_chkomx_Destroy_sse(cx2);
      p7_oprofile_Destroy_sse(om);
    }
  else
    {
      p7_omx_Destroy(oxf);
      p7_omx_Destroy(oxb);
      p7_chkomx_Destroy(cx1);
      p7_chkomx_Destroy(cx2);
      p7_oprofile_Destroy(om);
    }
}

#if defined(eslENABLE_SSE) && defined(eslENABLE_AVX)
/* utest_sse_vs_avx_fwdback()
 *
 * Run Forward_sse/Backward_sse and Forward_avx/Backward_avx on the same
 * sequences and profile, comparing scores.  SSE uses Q=ceil(M/4) stripes;
 * AVX uses Q=ceil(M/8).  Different summation orders from different stripe
 * widths cause small floating-point differences.  Tolerance is 0.01 nats,
 * tight enough to catch wrong implementations but loose enough to pass
 * expected rounding differences.
 * Skipped silently if AVX is not available at runtime.
 */
static void
utest_sse_vs_avx_fwdback(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char         msg[]  = "utest_sse_vs_avx_fwdback: scores differ";
  P7_HMM      *hmm    = NULL;
  P7_PROFILE  *gm     = NULL;
  P7_OPROFILE *om_tmp = NULL;
  P7_OPROFILE *om_sse = NULL;
  P7_OPROFILE *om_avx = NULL;
  ESL_DSQ     *dsq    = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *oxf_sse = NULL;
  P7_OMX      *oxb_sse = NULL;
  P7_OMX      *oxf_avx = NULL;
  P7_OMX      *oxb_avx = NULL;
  float        fsc_sse, bsc_sse;
  float        fsc_avx, bsc_avx;
  int          n       = N;

  if (!esl_cpu_has_avx()) { free(dsq); return; }

  p7_FLogsumInit();

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om_tmp);
  p7_oprofile_Destroy(om_tmp);

  om_sse = p7_oprofile_Create_sse(M, abc);
  om_avx = p7_oprofile_Create_avx(M, abc);
  p7_oprofile_Convert_sse(gm, om_sse);
  p7_oprofile_Convert_avx(gm, om_avx);

  oxf_sse = p7_omx_Create_sse(M, L, L);
  oxb_sse = p7_omx_Create_sse(M, L, L);
  oxf_avx = p7_omx_Create_avx(M, L, L);
  oxb_avx = p7_omx_Create_avx(M, L, L);

  while (n--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_Forward_sse (dsq, L, om_sse, oxf_sse, &fsc_sse);
      p7_Backward_sse(dsq, L, om_sse, oxf_sse, oxb_sse, &bsc_sse);
      p7_Forward_avx (dsq, L, om_avx, oxf_avx, &fsc_avx);
      p7_Backward_avx(dsq, L, om_avx, oxf_avx, oxb_avx, &bsc_avx);

      if (fabs(fsc_sse - fsc_avx) > 0.01) esl_fatal("%s: Forward sse=%.4f avx=%.4f",  msg, fsc_sse, fsc_avx);
      if (fabs(bsc_sse - bsc_avx) > 0.01) esl_fatal("%s: Backward sse=%.4f avx=%.4f", msg, bsc_sse, bsc_avx);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy_sse(oxf_sse);
  p7_omx_Destroy_sse(oxb_sse);
  p7_omx_Destroy_avx(oxf_avx);
  p7_omx_Destroy_avx(oxb_avx);
  p7_oprofile_Destroy_sse(om_sse);
  p7_oprofile_Destroy_avx(om_avx);
  p7_profile_Destroy(gm);
}
#endif /* eslENABLE_SSE && eslENABLE_AVX */

#endif /*p7FWDBACK_TESTDRIVE*/
/*---------------------- end, unit tests ------------------------*/




/*****************************************************************
 * 6. Test driver
 *****************************************************************/
#ifdef p7FWDBACK_TESTDRIVE
/* 
   gcc -g -Wall -msse2 -std=gnu99 -o fwdback_utest -I.. -L.. -I../../easel -L../../easel -Dp7FWDBACK_TESTDRIVE fwdback.c -lhmmer -leasel -lm
   ./fwdback_utest
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_avx.h"
#include <string.h>

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE Forward, Backward implementations";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = NULL;
  P7_BG          *bg   = NULL;
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");

  impl_Init();

  /* First round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))            == NULL)  esl_fatal("failed to create null model");

  utest_fwdback(r, abc, bg, M, L, N);   /* normal sized models */
  utest_fwdback(r, abc, bg, 1, L, 10);  /* size 1 models       */
  utest_fwdback(r, abc, bg, M, 1, 10);  /* size 1 sequences    */
  utest_backward_chk(r, abc, bg, M, L, N,  TRUE,  FALSE);
  utest_backward_chk(r, abc, bg, 1, L, 10, FALSE, FALSE);
  utest_backward_chk(r, abc, bg, M, 1, 10, FALSE, FALSE);
#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse()) {
    utest_backward_chk(r, abc, bg, M, L, N,  TRUE,  TRUE);
    utest_backward_chk(r, abc, bg, 1, L, 10, FALSE, TRUE);
    utest_backward_chk(r, abc, bg, M, 1, 10, FALSE, TRUE);
  }
#endif
#if defined(eslENABLE_SSE) && defined(eslENABLE_AVX)
  utest_sse_vs_avx_fwdback(r, abc, bg, M, L, N);
  utest_sse_vs_avx_fwdback(r, abc, bg, 1, L, 10);
  utest_sse_vs_avx_fwdback(r, abc, bg, M, 1, 10);
#endif

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  utest_fwdback(r, abc, bg, M, L, N);
  utest_fwdback(r, abc, bg, 1, L, 10);
  utest_fwdback(r, abc, bg, M, 1, 10);
  utest_backward_chk(r, abc, bg, M, L, N,  TRUE,  FALSE);
  utest_backward_chk(r, abc, bg, 1, L, 10, FALSE, FALSE);
  utest_backward_chk(r, abc, bg, M, 1, 10, FALSE, FALSE);
#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse()) {
    utest_backward_chk(r, abc, bg, M, L, N,  TRUE,  TRUE);
    utest_backward_chk(r, abc, bg, 1, L, 10, FALSE, TRUE);
    utest_backward_chk(r, abc, bg, M, 1, 10, FALSE, TRUE);
  }
#endif
#if defined(eslENABLE_SSE) && defined(eslENABLE_AVX)
  utest_sse_vs_avx_fwdback(r, abc, bg, M, L, N);
  utest_sse_vs_avx_fwdback(r, abc, bg, 1, L, 10);
  utest_sse_vs_avx_fwdback(r, abc, bg, M, 1, 10);
#endif

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*p7FWDBACK_TESTDRIVE*/
/*--------------------- end, test driver ------------------------*/



/*****************************************************************
 * 7. Example
 *****************************************************************/
#ifdef p7FWDBACK_EXAMPLE
/* Useful for debugging on small HMMs and sequences.
 * 
 * Compares to GForward().
 * 
   gcc -g -Wall -msse2 -std=gnu99 -o fwdback_example -I.. -L.. -I../../easel -L../../easel -Dp7FWDBACK_EXAMPLE fwdback.c -lhmmer -leasel -lm
   ./fwdback_example <hmmfile> <seqfile>
 */ 
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "impl_avx.h"
#include <string.h>

#define STYLES     "--fs,--sw,--ls,--s"	               /* Exclusive choice for alignment mode     */

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range  toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "output in one line awkable format",                0 },
  { "-P",        eslARG_NONE,   FALSE, NULL, NULL,   NULL,  NULL, NULL, "output in profmark format",                        0 },
  { "--fs",      eslARG_NONE,"default",NULL, NULL, STYLES,  NULL, NULL, "multihit local alignment",                         0 },
  { "--sw",      eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "unihit local alignment",                           0 },
  { "--ls",      eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "multihit glocal alignment",                        0 },
  { "--s",       eslARG_NONE,   FALSE, NULL, NULL, STYLES,  NULL, NULL, "unihit glocal alignment",                          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of Forward/Backward (SSE versions)";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_GMX         *gx      = NULL;
  P7_OMX         *fwd     = NULL;
  P7_OMX         *bck     = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           fraw, braw, nullsc, fsc;
  float           gfraw, gbraw, gfsc;
  double          P, gP;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  /* Open sequence file for reading */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

  /* create default null model, then create and optimize profile */
  bg = p7_bg_Create(abc);               
  p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc); 

  /* Now reconfig the models however we were asked to */
  if      (esl_opt_GetBoolean(go, "--fs"))  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL);
  else if (esl_opt_GetBoolean(go, "--sw"))  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_UNILOCAL);
  else if (esl_opt_GetBoolean(go, "--ls"))  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_GLOCAL);
  else if (esl_opt_GetBoolean(go, "--s"))   p7_ProfileConfig(hmm, bg, gm, sq->n, p7_UNIGLOCAL);

  /* now the optimized profile */
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  /* p7_oprofile_Dump(stdout, om);  */

  /* allocate DP matrices for O(M+L) parsers */
  fwd = p7_omx_Create(gm->M, 0, sq->n);
  bck = p7_omx_Create(gm->M, 0, sq->n);
  gx  = p7_gmx_Create(gm->M,    sq->n, sq->n, p7G_NSCELLS);

  /* allocate DP matrices for O(ML) fills */
  /* fwd = p7_omx_Create(gm->M, sq->n, sq->n); */
  /* bck = p7_omx_Create(gm->M, sq->n, sq->n); */

  /* p7_omx_SetDumpMode(stdout, fwd, TRUE); */     /* makes the fast DP algorithms dump their matrices */
  /* p7_omx_SetDumpMode(stdout, bck, TRUE); */  

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_ReconfigLength(gm,          sq->n);
      p7_bg_SetLength(bg,            sq->n);
      p7_omx_GrowTo(fwd, om->M, 0,   sq->n); 
      p7_omx_GrowTo(bck, om->M, 0,   sq->n); 
      p7_gmx_GrowTo(gx,  gm->M,      sq->n, sq->n);

      p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);
    
      p7_ForwardParser (sq->dsq, sq->n, om,      fwd, &fraw);
      p7_BackwardParser(sq->dsq, sq->n, om, fwd, bck, &braw);

      /* p7_Forward (sq->dsq, sq->n, om,      fwd, &fsc);        printf("forward:              %.2f nats\n", fsc);  */
      /* p7_Backward(sq->dsq, sq->n, om, fwd, bck, &bsc);        printf("backward:             %.2f nats\n", bsc);  */

      /* Comparison to other F/B implementations */
      p7_GForward     (sq->dsq, sq->n, gm, gx,  &gfraw);
      p7_GBackward    (sq->dsq, sq->n, gm, gx,  &gbraw);

      /* p7_gmx_Dump(stdout, gx, p7_DEFAULT);  */

      fsc  =  (fraw-nullsc) / eslCONST_LOG2;
      gfsc = (gfraw-nullsc) / eslCONST_LOG2;
      P  = esl_exp_surv(fsc,   om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
      gP = esl_exp_surv(gfsc,  gm->evparam[p7_FTAU],  gm->evparam[p7_FLAMBDA]);

      if (esl_opt_GetBoolean(go, "-1"))
	{
	  printf("%-30s\t%-20s\t%9.2g\t%6.1f\t%9.2g\t%6.1f\n", sq->name, hmm->name, P, fsc, gP, gfsc);
	}
      else if (esl_opt_GetBoolean(go, "-P"))
	{ /* output suitable for direct use in profmark benchmark postprocessors: */
	  printf("%g\t%.2f\t%s\t%s\n", P, fsc, sq->name, hmm->name);
	}
      else
	{
	  printf("target sequence:      %s\n",        sq->name);
	  printf("fwd filter raw score: %.2f nats\n", fraw);
	  printf("bck filter raw score: %.2f nats\n", braw);
	  printf("null score:           %.2f nats\n", nullsc);
	  printf("per-seq score:        %.2f bits\n", fsc);
	  printf("P-value:              %g\n",        P);
	  printf("GForward raw score:   %.2f nats\n", gfraw);
	  printf("GBackward raw score:  %.2f nats\n", gbraw);
	  printf("GForward seq score:   %.2f bits\n", gfsc);
	  printf("GForward P-value:     %g\n",        gP);
	}

      esl_sq_Reuse(sq);
    }

  /* cleanup */
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_omx_Destroy(bck);
  p7_omx_Destroy(fwd);
  p7_gmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7FWDBACK_EXAMPLE*/

