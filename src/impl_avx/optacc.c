/* Optimal accuracy alignment dispatcher.
 * Provides the non-suffixed p7_OptimalAccuracy() and p7_OATrace() declared
 * in impl_avx.h.  Delegates to the fastest available ISA implementation
 * at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Forward declaration of dispatcher */
static int p7_OptimalAccuracy_Dispatcher(const P7_OPROFILE *om, const P7_OMX *pp, P7_OMX *ox, float *ret_e);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_OptimalAccuracy)(const P7_OPROFILE *om, const P7_OMX *pp, P7_OMX *ox, float *ret_e) = p7_OptimalAccuracy_Dispatcher;

/* Function:  p7_OptimalAccuracy_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for optimal accuracy DP fill.
 *            <ox> may be <pp>: the OA matrix can overwrite the
 *            posterior decoding matrix.
 *
 * Returns:   <eslOK> on success.
 */
static int
p7_OptimalAccuracy_Dispatcher(const P7_OPROFILE *om, const P7_OMX *pp, P7_OMX *ox, float *ret_e)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_OptimalAccuracy = p7_OptimalAccuracy_avx512; return p7_OptimalAccuracy_avx512(om, pp, ox, ret_e); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_OptimalAccuracy = p7_OptimalAccuracy_avx;    return p7_OptimalAccuracy_avx(om, pp, ox, ret_e); }
#endif
#ifdef eslENABLE_SSE
  p7_OptimalAccuracy = p7_OptimalAccuracy_sse;
  return p7_OptimalAccuracy_sse(om, pp, ox, ret_e);
#else
  p7_Die("p7_OptimalAccuracy: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_OATrace_Dispatcher(const P7_OPROFILE *om, const P7_OMX *ox, P7_TRACE *tr);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_OATrace)(const P7_OPROFILE *om, const P7_OMX *ox, P7_TRACE *tr) = p7_OATrace_Dispatcher;

/* Function:  p7_OATrace_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for optimal accuracy traceback.
 *            The traceback reads only the OA matrix <ox>, which may have
 *            been written over the posterior decoding matrix; posterior
 *            probabilities are recovered from <ox>.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> if trace is not empty.
 */
static int
p7_OATrace_Dispatcher(const P7_OPROFILE *om, const P7_OMX *ox, P7_TRACE *tr)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_OATrace = p7_OATrace_avx512; return p7_OATrace_avx512(om, ox, tr); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_OATrace = p7_OATrace_avx;    return p7_OATrace_avx(om, ox, tr); }
#endif
#ifdef eslENABLE_SSE
  p7_OATrace = p7_OATrace_sse;
  return p7_OATrace_sse(om, ox, tr);
#else
  p7_Die("p7_OATrace: no SIMD implementation available");
  return eslENORESULT;
#endif
}
/* Forward declarations of dispatchers */
static int p7_OptimalAccuracySweep_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, P7_OATMX *ox, float *opt_null2, float *opt_fsc, float *ret_oasc);
static int p7_OATraceSweep_Dispatcher(const P7_OPROFILE *om, const P7_OATMX *ox, P7_TRACE *tr);

/* Global function pointers, initially pointing at the dispatchers */
int (*p7_OptimalAccuracySweep)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, P7_OATMX *ox, float *opt_null2, float *opt_fsc, float *ret_oasc) = p7_OptimalAccuracySweep_Dispatcher;
int (*p7_OATraceSweep)(const P7_OPROFILE *om, const P7_OATMX *ox, P7_TRACE *tr) = p7_OATraceSweep_Dispatcher;

/* Function:  p7_OptimalAccuracySweep_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for the OA sweep over a
 *            checkpointed Backward. See p7_OptimalAccuracySweep_avx().
 */
static int
p7_OptimalAccuracySweep_Dispatcher(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, P7_OATMX *ox, float *opt_null2, float *opt_fsc, float *ret_oasc)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_OptimalAccuracySweep = p7_OptimalAccuracySweep_avx; return p7_OptimalAccuracySweep_avx(dsq, L, om, cx, ox, opt_null2, opt_fsc, ret_oasc); }
#endif
#ifdef eslENABLE_SSE
  p7_OptimalAccuracySweep = p7_OptimalAccuracySweep_sse;
  return p7_OptimalAccuracySweep_sse(dsq, L, om, cx, ox, opt_null2, opt_fsc, ret_oasc);
#else
  p7_Die("p7_OptimalAccuracySweep: no SIMD implementation available");
  return eslENORESULT;
#endif
}

/* Function:  p7_OATraceSweep_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for the OA trace from what the OA
 *            sweep stored. See p7_OATraceSweep_avx().
 */
static int
p7_OATraceSweep_Dispatcher(const P7_OPROFILE *om, const P7_OATMX *ox, P7_TRACE *tr)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_OATraceSweep = p7_OATraceSweep_avx; return p7_OATraceSweep_avx(om, ox, tr); }
#endif
#ifdef eslENABLE_SSE
  p7_OATraceSweep = p7_OATraceSweep_sse;
  return p7_OATraceSweep_sse(om, ox, tr);
#else
  p7_Die("p7_OATraceSweep: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/
#ifdef p7OPTACC_BENCHMARK
/*
   gcc -g -O3 -mavx2 -std=gnu99 -o benchmark-optacc -I.. -L.. -I../../easel -L../../easel -Dp7OPTACC_BENCHMARK optacc.c -lhmmer -leasel -lm
   icc  -O3 -static -o benchmark-optacc -I.. -L.. -I../../easel -L../../easel -Dp7OPTACC_BENCHMARK optacc.c -lhmmer -leasel -lm

   ./benchmark-optacc <hmmfile>         runs benchmark on optimal accuracy fill and trace
   ./benchmark-optacc -c -N1 <hmmfile>  compare scores of AVX version to generic impl
   ./benchmark-optacc -x -N1 <hmmfile>  test that scores match trusted implementation.

                    RRM_1 (M=72)       Caudal_act (M=136)     SMC_N (M=1151)
                 -----------------    ------------------     ---------------
   20 Aug 08:     13.11u (110 Mc/s)     23.39u (116 Mc/s)    332.62u (69 Mc/s)

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

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-x", "compare scores to generic implementation (debug)", 0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-c", "equate scores to trusted implementation (debug)",  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  { "--notrace", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark the DP fill stage",                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for optimal accuracy alignment, AVX version";

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
  P7_GMX         *gx1     = NULL;
  P7_GMX         *gx2     = NULL;
  P7_OMX         *ox1     = NULL;
  P7_OMX         *ox2     = NULL;
  P7_TRACE       *tr      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           fsc, bsc, accscore;
  float           fsc_g, bsc_g, accscore_g;
  double          Mcs;

  p7_FLogsumInit();

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);                 p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);    p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);    p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);

  if (esl_opt_GetBoolean(go, "-x") && p7_FLogsumError(-0.4, -0.5) > 0.0001)
    p7_Fail("-x here requires p7_Logsum() recompiled in slow exact mode");

  ox1 = p7_omx_Create(gm->M, L, L);
  ox2 = p7_omx_Create(gm->M, L, L);
  tr  = p7_trace_CreateWithPP();

  esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  p7_Forward (dsq, L, om, ox1,      &fsc);
  p7_Backward(dsq, L, om, ox1, ox2, &bsc);
  p7_Decoding(om, ox1, ox2, ox2);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      p7_OptimalAccuracy(om, ox2, ox1, &accscore);

      if (! esl_opt_GetBoolean(go, "--notrace"))
	{
	  p7_OATrace(om, ox1, tr);
	  p7_trace_Reuse(tr);
	}
    }
  esl_stopwatch_Stop(w);

  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) w->user;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  if (esl_opt_GetBoolean(go, "-c") || esl_opt_GetBoolean(go, "-x") )
    {
      gx1 = p7_gmx_Create(gm->M, L, L, p7G_NSCELLS);
      gx2 = p7_gmx_Create(gm->M, L, L, p7G_NSCELLS);

      p7_GForward (dsq, L, gm, gx1, &fsc_g);
      p7_GBackward(dsq, L, gm, gx2, &bsc_g);
      p7_GDecoding(gm, gx1, gx2, gx2);
      p7_GOptimalAccuracy(gm, gx2, gx1, &accscore_g);

      printf("generic:  fwd=%8.4f  bck=%8.4f  acc=%8.4f\n", fsc_g, bsc_g, accscore_g);
      printf("AVX:      fwd=%8.4f  bck=%8.4f  acc=%8.4f\n", fsc,   bsc,   accscore);

      p7_gmx_Destroy(gx1);
      p7_gmx_Destroy(gx2);
    }

  free(dsq);
  p7_omx_Destroy(ox1);
  p7_omx_Destroy(ox2);
  p7_trace_Destroy(tr);
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
#endif /*p7OPTACC_BENCHMARK*/
/*---------------- end, benchmark driver ------------------------*/




/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7OPTACC_TESTDRIVE
#include <math.h>
#include <string.h>

#include "esl_randomseq.h"
#include "esl_vectorops.h"

#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
/* 
 * 1. Compare accscore to GOptimalAccuracy().
 * 2. Compare trace to GOATrace().
 * 
 * Note: This test is subject to some expected noise and can fail
 * for entirely innocent reasons. Generic Forward/Backward calculations with
 * p7_GForward(), p7_GBackward() use coarse-grain table lookups to sum
 * log probabilities, and sufficient roundoff error can accumulate to
 * change the optimal accuracy traceback, causing this test to fail.
 * So, if optacc_utest fails, before you go looking for bugs, first
 * go to ../logsum.c, change the #ifdef to activate the slow/accurate 
 * version, recompile and rerun optacc_utest. If the failure goes away,
 * you can ignore it.   - SRE, Wed Dec 17 09:45:31 2008
 */
static void
utest_optacc(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char        *msg = "optimal accuracy unit test failed";
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_SQ      *sq  = esl_sq_CreateDigital(abc);
  P7_OMX      *ox1 = p7_omx_Create(M, L, L);
  P7_OMX      *ox2 = p7_omx_Create(M, L, L);
  P7_GMX      *gx1 = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX      *gx2 = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_TRACE    *tr  = p7_trace_CreateWithPP();
  P7_TRACE    *trg = p7_trace_CreateWithPP();
  P7_TRACE    *tro = p7_trace_CreateWithPP();
  float        accscore_o;
  float        fsc, bsc, accscore;
  float        fsc_g, bsc_g, accscore_g, accscore_g2;
  float        pptol = 0.01;
  float        sctol = 0.001;
  float        gtol;

  p7_FLogsumInit();
  gtol = ( (p7_FLogsumError(-0.4, -0.5) > 0.0001) ?  0.1 : 0.001);

  if (p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om)!= eslOK) esl_fatal(msg);
  while (N--)
    {
      if (p7_ProfileEmit(r, hmm, gm, bg, sq, tro)          != eslOK) esl_fatal(msg);

      if (p7_omx_GrowTo(ox1, M, sq->n, sq->n)              != eslOK) esl_fatal(msg);
      if (p7_omx_GrowTo(ox2, M, sq->n, sq->n)              != eslOK) esl_fatal(msg);
      if (p7_gmx_GrowTo(gx1, M, sq->n, sq->n) != eslOK) esl_fatal(msg);
      if (p7_gmx_GrowTo(gx2, M, sq->n, sq->n) != eslOK) esl_fatal(msg);

      if (p7_Forward (sq->dsq, sq->n, om, ox1,      &fsc) != eslOK) esl_fatal(msg);
      if (p7_Backward(sq->dsq, sq->n, om, ox1, ox2, &bsc) != eslOK) esl_fatal(msg);
      if (p7_Decoding(om, ox1, ox2, ox2)                  != eslOK) esl_fatal(msg);
      if (p7_OptimalAccuracy(om, ox2, ox1, &accscore)     != eslOK) esl_fatal(msg);

#if 0
      p7_omx_FDeconvert(ox1, gx1); 
      p7_gmx_Dump(stdout, gx1, p7_DEFAULT); 
      p7_omx_FDeconvert(ox2, gx1); 
      p7_gmx_Dump(stdout, gx1, p7_DEFAULT); 
#endif
      if (p7_OATrace(om, ox1, tr)                    != eslOK) esl_fatal(msg);
      
      if (p7_GForward (sq->dsq, sq->n, gm, gx1, &fsc_g)   != eslOK) esl_fatal(msg);
      if (p7_GBackward(sq->dsq, sq->n, gm, gx2, &bsc_g)   != eslOK) esl_fatal(msg);

#if 0
      p7_gmx_Dump(stdout, gx1, p7_DEFAULT); /* fwd */
      p7_gmx_Dump(stdout, gx2, p7_DEFAULT); /* bck */
#endif

      if (p7_GDecoding(gm, gx1, gx2, gx2)                 != eslOK) esl_fatal(msg);
      if (p7_GOptimalAccuracy(gm, gx2, gx1, &accscore_g)  != eslOK) esl_fatal(msg);
      
#if 0
      p7_gmx_Dump(stdout, gx1, p7_DEFAULT); /* oa */
      p7_gmx_Dump(stdout, gx2, p7_DEFAULT); /* pp */
#endif
      if (p7_GOATrace(gm, gx1, trg)                       != eslOK) esl_fatal(msg);

      if (p7_trace_SetPP(tro, gx2)                        != eslOK) esl_fatal(msg);

      if (esl_opt_GetBoolean(go, "--traces"))
	{
	  p7_trace_Dump(stdout, tro, gm, sq->dsq);
	  p7_trace_Dump(stdout, tr,  gm, sq->dsq);
	  p7_trace_Dump(stdout, trg, gm, sq->dsq);
	}

      if (p7_trace_Validate(tr,  abc, sq->dsq, NULL)      != eslOK) esl_fatal(msg);
      if (p7_trace_Validate(trg, abc, sq->dsq, NULL)      != eslOK) esl_fatal(msg);
      if (p7_trace_Compare(tr, trg, pptol)                != eslOK) esl_fatal(msg);

      accscore_o  = p7_trace_GetExpectedAccuracy(tro); /* according to gx2; see p7_trace_SetPP() call above */
      accscore_g2 = p7_trace_GetExpectedAccuracy(trg);

#if 0
      printf("%f %f %f %f\n", accscore, accscore_g, accscore_g2, accscore_o);
#endif

      if (esl_FCompare_old(fsc,        bsc,         sctol)    != eslOK) esl_fatal(msg);
      if (esl_FCompare_old(fsc_g,      bsc_g,       gtol)     != eslOK) esl_fatal(msg);
      if (esl_FCompare_old(fsc,        fsc_g,       gtol)     != eslOK) esl_fatal(msg);
      if (esl_FCompare_old(accscore,   accscore_g,  gtol)     != eslOK) esl_fatal(msg);
      if (esl_FCompare_old(accscore_g, accscore_g2, gtol)     != eslOK) esl_fatal(msg);
      if (accscore_g2 < accscore_o)                                 esl_fatal(msg);
      /* the above deserves explanation:
       *  - accscore_o is the accuracy of the originally emitted trace, according
       *      to the generic posterior decoding matrix <gx2>. This is a lower bound
       *      on the expected # of accurately aligned residues found by a DP 
       *      optimization.
       *  - accscore is the accuracy found by the fast (vector) code DP implementation.
       *  - accscore_g is the accuracy found by the generic DP implementation.
       *      accscore and accscore_g should be nearly identical,
       *      within tolerance of roundoff error accumulation and
       *      the imprecision of Logsum() tables.
       *  - accscore_g2 is the accuracy of the traceback identified by the generic
       *      DP implementation. It should be identical (within order-of-evaluation
       *      roundoff error) to accscore_g.
       *      
       * the "accscore_g2 < accscore_o" test is carefully contrived.
       * accscore_o is a theoretical lower bound but because of fp error, 
       * accscore and (much more rarely) even accscore_g can exceed accscore_o.
       * accscore_g2, however, is calculated with identical order of evaluation
       * as accscore_o if the optimal trace does turn out to be identical to 
       * the originally emitted trace. It should be extremely unlikely (though
       * not impossible) for accscore_o to exceed accscore_g2. (The DP algorithm
       * would have to identify a trace that was different than the original trace,
       * which the DP algorithm, by order-of-evaluation, assigned higher accuracy,
       * but order-of-evaluation in traceback dependent code assigned lower accuracy.
       * [xref J5/29]
       */

      esl_sq_Reuse(sq);
      p7_trace_Reuse(tr);
      p7_trace_Reuse(trg);
      p7_trace_Reuse(tro);
    }

  p7_trace_Destroy(tro);
  p7_trace_Destroy(trg);
  p7_trace_Destroy(tr);
  p7_gmx_Destroy(gx2);
  p7_gmx_Destroy(gx1);
  p7_omx_Destroy(ox2);
  p7_omx_Destroy(ox1);  
  esl_sq_Destroy(sq);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}


/* utest_inplace()
 * On sequences emitted from a sampled profile, the OA fill written over
 * the posterior matrix (<ox> = <pp>) must give the same OA matrix and
 * score as a fill into a separate matrix, and the traceback from either
 * matrix must be the same. The traceback's posterior probabilities,
 * recovered from the OA matrix by subtraction, must agree with the
 * posterior matrix within <pptol>.
 * With <use_sse>, the _sse versions are called directly, so an AVX CPU
 * also tests the SSE implementation.
 */
static void
utest_inplace(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N, int use_sse)
{
  char        *msg   = "optimal accuracy in-place unit test failed";
  P7_HMM      *hmm   = NULL;
  P7_PROFILE  *gm    = NULL;
  P7_OPROFILE *om    = NULL;
  ESL_SQ      *sq    = esl_sq_CreateDigital(abc);
  P7_OMX      *ox1   = NULL;
  P7_OMX      *ox2   = NULL;
  P7_GMX      *gpp   = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX      *gox1  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX      *gox2  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_TRACE    *tro   = p7_trace_Create();
  P7_TRACE    *tr1   = p7_trace_CreateWithPP();
  P7_TRACE    *tr2   = p7_trace_CreateWithPP();
  float        pptol = 0.001;
  float        oasc1, oasc2;
  int          z;

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  if (use_sse)
    {
      p7_oprofile_Destroy(om);
      om  = p7_oprofile_Create_sse(M, abc);
      p7_oprofile_Convert_sse(gm, om);
      ox1 = p7_omx_Create_sse(M, L, L);
      ox2 = p7_omx_Create_sse(M, L, L);
    }
  else
    {
      ox1 = p7_omx_Create(M, L, L);
      ox2 = p7_omx_Create(M, L, L);
    }
  while (N--)
    {
      if (p7_ProfileEmit(r, hmm, gm, bg, sq, tro)                  != eslOK) esl_fatal(msg);
      if (p7_gmx_GrowTo(gpp,  M, sq->n, sq->n)                     != eslOK) esl_fatal(msg);
      if (p7_gmx_GrowTo(gox1, M, sq->n, sq->n)                     != eslOK) esl_fatal(msg);
      if (p7_gmx_GrowTo(gox2, M, sq->n, sq->n)                     != eslOK) esl_fatal(msg);
      if (use_sse)
        {
          if (p7_omx_GrowTo_sse(ox1, M, sq->n, sq->n)              != eslOK) esl_fatal(msg);
          if (p7_omx_GrowTo_sse(ox2, M, sq->n, sq->n)              != eslOK) esl_fatal(msg);
          if (p7_Forward_sse (sq->dsq, sq->n, om, ox1,      NULL)  != eslOK) esl_fatal(msg);
          if (p7_Backward_sse(sq->dsq, sq->n, om, ox1, ox2, NULL)  != eslOK) esl_fatal(msg);
          if (p7_Decoding_sse(om, ox1, ox2, ox2)                   != eslOK) esl_fatal(msg);
          if (p7_omx_FDeconvert_sse(ox2, gpp)                      != eslOK) esl_fatal(msg);
          if (p7_OptimalAccuracy_sse(om, ox2, ox1, &oasc1)         != eslOK) esl_fatal(msg);
          if (p7_OATrace_sse(om, ox1, tr1)                         != eslOK) esl_fatal(msg);
          if (p7_OptimalAccuracy_sse(om, ox2, ox2, &oasc2)         != eslOK) esl_fatal(msg);
          if (p7_OATrace_sse(om, ox2, tr2)                         != eslOK) esl_fatal(msg);
          if (p7_omx_FDeconvert_sse(ox1, gox1)                     != eslOK) esl_fatal(msg);
          if (p7_omx_FDeconvert_sse(ox2, gox2)                     != eslOK) esl_fatal(msg);
        }
      else
        {
          if (p7_omx_GrowTo(ox1, M, sq->n, sq->n)                  != eslOK) esl_fatal(msg);
          if (p7_omx_GrowTo(ox2, M, sq->n, sq->n)                  != eslOK) esl_fatal(msg);
          if (p7_Forward (sq->dsq, sq->n, om, ox1,      NULL)      != eslOK) esl_fatal(msg);
          if (p7_Backward(sq->dsq, sq->n, om, ox1, ox2, NULL)      != eslOK) esl_fatal(msg);
          if (p7_Decoding(om, ox1, ox2, ox2)                       != eslOK) esl_fatal(msg);
          if (p7_omx_FDeconvert(ox2, gpp)                          != eslOK) esl_fatal(msg);
          if (p7_OptimalAccuracy(om, ox2, ox1, &oasc1)             != eslOK) esl_fatal(msg);
          if (p7_OATrace(om, ox1, tr1)                             != eslOK) esl_fatal(msg);
          if (p7_OptimalAccuracy(om, ox2, ox2, &oasc2)             != eslOK) esl_fatal(msg);
          if (p7_OATrace(om, ox2, tr2)                             != eslOK) esl_fatal(msg);
          if (p7_omx_FDeconvert(ox1, gox1)                         != eslOK) esl_fatal(msg);
          if (p7_omx_FDeconvert(ox2, gox2)                         != eslOK) esl_fatal(msg);
        }

      if (oasc1 != oasc2)                                        esl_fatal("%s: OA scores differ", msg);
      if (p7_gmx_Compare(gox2, gox1, 0.0)              != eslOK) esl_fatal("%s: OA matrices differ", msg);
      if (tr1->N != tr2->N)                                      esl_fatal("%s: traces differ", msg);
      for (z = 0; z < tr1->N; z++)
        if (tr1->st[z] != tr2->st[z] || tr1->k[z] != tr2->k[z] || tr1->i[z] != tr2->i[z] || tr1->pp[z] != tr2->pp[z])
          esl_fatal("%s: traces differ", msg);
      for (z = 0; z < tr1->N; z++)
        {
          if      (tr1->st[z] == p7T_M) { if (fabs(tr1->pp[z] - gpp->dp[tr1->i[z]][tr1->k[z]*p7G_NSCELLS + p7G_M]) > pptol) esl_fatal("%s: M posterior", msg); }
          else if (tr1->st[z] == p7T_I) { if (fabs(tr1->pp[z] - gpp->dp[tr1->i[z]][tr1->k[z]*p7G_NSCELLS + p7G_I]) > pptol) esl_fatal("%s: I posterior", msg); }
          else if ((tr1->st[z] == p7T_N || tr1->st[z] == p7T_C || tr1->st[z] == p7T_J) && tr1->i[z] > 0 && tr1->st[z-1] == tr1->st[z])
            {
              if (fabs(tr1->pp[z] - gpp->xmx[tr1->i[z]*p7G_NXCELLS + (tr1->st[z] == p7T_N ? p7G_N : tr1->st[z] == p7T_C ? p7G_C : p7G_J)]) > pptol) esl_fatal("%s: NCJ posterior", msg);
            }
        }

      esl_sq_Reuse(sq);
      p7_trace_Reuse(tro);
      p7_trace_Reuse(tr1);
      p7_trace_Reuse(tr2);
    }

  p7_trace_Destroy(tro);
  p7_trace_Destroy(tr1);
  p7_trace_Destroy(tr2);
  p7_gmx_Destroy(gpp);
  p7_gmx_Destroy(gox1);
  p7_gmx_Destroy(gox2);
  if (use_sse)
    {
      p7_omx_Destroy_sse(ox1);
      p7_omx_Destroy_sse(ox2);
      p7_oprofile_Destroy_sse(om);
    }
  else
    {
      p7_omx_Destroy(ox1);
      p7_omx_Destroy(ox2);
      p7_oprofile_Destroy(om);
    }
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}
/* utest_sweep()
 * The OA sweep over a checkpointed Backward against the full-matrix
 * path on the same Backward: p7_Forward(), the checkpointed Backward
 * filled in one block and copied into a Backward matrix with its own
 * scale factors, p7_Decoding(), p7_OptimalAccuracy(), p7_OATrace() and
 * p7_Null2_ByExpectation(). The sweep, with Backward blocks of 1..8
 * rows, must give the same Forward and OA scores, null2, states, k and
 * i on the trace, and PP characters equal to those the posterior
 * matrix gives. Sequences alternate random ones and the consensus of a
 * model with sharpened match emissions and M->M transitions between
 * random flanks, so that rows are rescaled. With <use_sse>, the _sse
 * versions are called directly, so an AVX CPU also tests the SSE engine.
 */
static void
utest_sweep(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N, int use_sse)
{
  char        *msg  = "OA sweep unit test failed";
  P7_HMM      *hmm  = NULL;
  P7_PROFILE  *gm   = NULL;
  P7_OPROFILE *om   = NULL;
  ESL_SQ      *sq   = esl_sq_CreateDigital(abc);
  ESL_DSQ     *dsq  = malloc(sizeof(ESL_DSQ) * (L+2));
  float       *n2a  = malloc(sizeof(float) * abc->Kp);
  float       *n2b  = malloc(sizeof(float) * abc->Kp);
  P7_OMX      *fwd  = NULL, *bck = NULL, *pp = NULL, *oa = NULL;
  P7_CHKOMX   *cx1  = NULL, *cx2 = NULL;
  P7_OATMX    *ox   = p7_oatmx_Create(M, L, p7OAT_NPLANES, 1);
  P7_TRACE    *tr1  = p7_trace_CreateWithPP();
  P7_TRACE    *tr2  = p7_trace_CreateWithPP();
  int          Q    = use_sse ? p7O_NQF(M) : p7O_NQF_AVX(M);
  size_t       rowb = (use_sse ? 16 : 32) * 3 * (size_t) Q;
  float        fsc, fsc2, oasc1, oasc2, p;
  int          n, Ln, blockR, i, k, x, z;

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
  p7_emit_SimpleConsensus(hmm, sq);

  if (use_sse)
    {
      om  = p7_oprofile_Create_sse(M, abc);  p7_oprofile_Convert_sse(gm, om);
      fwd = p7_omx_Create_sse(M, L, L); bck = p7_omx_Create_sse(M, L, L); pp = p7_omx_Create_sse(M, L, L); oa = p7_omx_Create_sse(M, L, L);
      cx1 = p7_chkomx_Create_sse(M, L, L+1, 1, p7X_NSCELLS, p7_CHK_BACKWARD);
      cx2 = p7_chkomx_Create_sse(M, L, 1,   1, p7X_NSCELLS, p7_CHK_BACKWARD);
    }
  else
    {
      om  = p7_oprofile_Create(M, abc);  p7_oprofile_Convert(gm, om);
      fwd = p7_omx_Create(M, L, L); bck = p7_omx_Create(M, L, L); pp = p7_omx_Create(M, L, L); oa = p7_omx_Create(M, L, L);
      cx1 = p7_chkomx_Create(M, L, L+1, 1, p7X_NSCELLS, p7_CHK_BACKWARD);
      cx2 = p7_chkomx_Create(M, L, 1,   1, p7X_NSCELLS, p7_CHK_BACKWARD);
    }

  for (n = 0; n < N; n++)
    {
      Ln = 1 + esl_rnd_Roll(r, L);
      esl_rsq_xfIID(r, bg->f, abc->K, Ln, dsq);
      if (n % 2 == 1 && M <= L)
        {
          Ln = ESL_MAX(Ln, M);
          esl_rsq_xfIID(r, bg->f, abc->K, Ln, dsq);
          memcpy(dsq + 1 + esl_rnd_Roll(r, Ln - M + 1), sq->dsq + 1, sizeof(ESL_DSQ) * M);
        }
      blockR = 1 + esl_rnd_Roll(r, 8);
      p7_trace_Reuse(tr1);
      p7_trace_Reuse(tr2);

      if (use_sse)
        {
          p7_oprofile_ReconfigLength_sse(om, Ln);
          if (p7_Forward_sse(dsq, Ln, om, fwd, &fsc)                    != eslOK) esl_fatal(msg);
          if (p7_chkomx_GrowTo_sse(cx1, M, Ln, Ln+1)                    != eslOK) esl_fatal(msg);
          if (p7_BackwardChk_sse(dsq, Ln, om, cx1, NULL)                != eslOK) esl_fatal(msg);
          for (i = 0; i <= Ln; i++) memcpy(bck->dpf[i], cx1->dpf[i], rowb);
        }
      else
        {
          p7_oprofile_ReconfigLength(om, Ln);
          if (p7_Forward(dsq, Ln, om, fwd, &fsc)                        != eslOK) esl_fatal(msg);
          if (p7_chkomx_GrowTo(cx1, M, Ln, Ln+1)                        != eslOK) esl_fatal(msg);
          if (p7_BackwardChk(dsq, Ln, om, cx1, NULL)                    != eslOK) esl_fatal(msg);
          for (i = 0; i <= Ln; i++) memcpy(bck->dpf_avx[i], cx1->dpf_avx[i], rowb);
        }
      memcpy(bck->xmx, cx1->xmx, sizeof(float) * (Ln+1) * p7X_NXCELLS);
      bck->M = M; bck->L = Ln; bck->has_own_scales = TRUE;

      if (use_sse)
        {
          if (p7_Decoding_sse(om, fwd, bck, pp)                         != eslOK) esl_fatal(msg);
          p7_OptimalAccuracy_sse(om, pp, oa, &oasc1);
          if (p7_OATrace_sse(om, oa, tr1)                               != eslOK) esl_fatal(msg);
          if (p7_chkomx_GrowTo_sse(cx2, M, Ln, blockR)                  != eslOK) esl_fatal(msg);
          if (p7_BackwardChk_sse(dsq, Ln, om, cx2, NULL)                != eslOK) esl_fatal(msg);
          if (p7_OptimalAccuracySweep_sse(dsq, Ln, om, cx2, ox, n2b, &fsc2, &oasc2) != eslOK) esl_fatal(msg);
          if (p7_OATraceSweep_sse(om, ox, tr2)                          != eslOK) esl_fatal(msg);
        }
      else
        {
          if (p7_Decoding(om, fwd, bck, pp)                             != eslOK) esl_fatal(msg);
          p7_OptimalAccuracy(om, pp, oa, &oasc1);
          if (p7_OATrace(om, oa, tr1)                                   != eslOK) esl_fatal(msg);
          if (p7_chkomx_GrowTo(cx2, M, Ln, blockR)                      != eslOK) esl_fatal(msg);
          if (p7_BackwardChk(dsq, Ln, om, cx2, NULL)                    != eslOK) esl_fatal(msg);
          if (p7_OptimalAccuracySweep(dsq, Ln, om, cx2, ox, n2b, &fsc2, &oasc2) != eslOK) esl_fatal(msg);
          if (p7_OATraceSweep(om, ox, tr2)                              != eslOK) esl_fatal(msg);
        }

      if (fsc2  != fsc)    esl_fatal("%s: Forward scores differ (%f %f)", msg, fsc2, fsc);
      if (oasc2 != oasc1)  esl_fatal("%s: OA scores differ (%f %f)", msg, oasc2, oasc1);
      if (tr1->N != tr2->N) esl_fatal("%s: trace lengths differ (%d %d)", msg, tr1->N, tr2->N);
      for (z = 0; z < tr1->N; z++)
        {
          if (tr1->st[z] != tr2->st[z] || tr1->k[z] != tr2->k[z] || tr1->i[z] != tr2->i[z])
            esl_fatal("%s: traces differ at %d", msg, z);
          switch (tr1->st[z]) {
          case p7T_M: p = use_sse ? p7_omx_FGetMDI_sse(pp, p7X_M, tr1->i[z], tr1->k[z]) : p7_omx_FGetMDI(pp, p7X_M, tr1->i[z], tr1->k[z]); break;
          case p7T_I: p = use_sse ? p7_omx_FGetMDI_sse(pp, p7X_I, tr1->i[z], tr1->k[z]) : p7_omx_FGetMDI(pp, p7X_I, tr1->i[z], tr1->k[z]); break;
          case p7T_N: case p7T_C: case p7T_J:
            /* an N, C or J emits residue i when the state before it is the same one */
            p = (z > 0 && tr1->st[z-1] == tr1->st[z]) ? pp->xmx[tr1->i[z]*p7X_NXCELLS + (tr1->st[z] == p7T_N ? p7X_N : tr1->st[z] == p7T_C ? p7X_C : p7X_J)] : 0.0;
            break;
          default: p = 0.0; break;
          }
          if (p7_alidisplay_EncodePostProb(tr2->pp[z]) != p7_alidisplay_EncodePostProb(p))
            esl_fatal("%s: PP character differs at %d (state %d, i %d, k %d)", msg, z, tr1->st[z], tr1->i[z], tr1->k[z]);
        }

      /* null2 from the full posterior matrix, which it uses as workspace, so last */
      if (use_sse) p7_Null2_ByExpectation_sse(om, pp, n2a);
      else         p7_Null2_ByExpectation    (om, pp, n2a);
      for (x = 0; x < abc->Kp; x++)
        if (n2b[x] != n2a[x]) esl_fatal("%s: null2 differs for residue %d (%g %g)", msg, x, n2b[x], n2a[x]);
    }

  free(dsq); free(n2a); free(n2b);
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr1);
  p7_trace_Destroy(tr2);
  p7_oatmx_Destroy(ox);
  if (use_sse)
    {
      p7_omx_Destroy_sse(fwd); p7_omx_Destroy_sse(bck); p7_omx_Destroy_sse(pp); p7_omx_Destroy_sse(oa);
      p7_chkomx_Destroy_sse(cx1); p7_chkomx_Destroy_sse(cx2);
      p7_oprofile_Destroy_sse(om);
    }
  else
    {
      p7_omx_Destroy(fwd); p7_omx_Destroy(bck); p7_omx_Destroy(pp); p7_omx_Destroy(oa);
      p7_chkomx_Destroy(cx1); p7_chkomx_Destroy(cx2);
      p7_oprofile_Destroy(om);
    }
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}

#endif /*p7OPTACC_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/




/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7OPTACC_TESTDRIVE

/* Failures in this test are to be expected, if you change the defaults.
 * The default RNG seed of 41 is very carefully chosen! Most seeds will
 * cause this test to fail. (Only 13 and 41 work for seeds 1..50.)
 * 
 * The generic fwd/bck algorithms work in log space and suffer from a
 * small amount of imprecision, not only from the use of FLogsum()'s
 * table-driven approximation, but even (apparently) inherent in log()
 * and exp(). To minimize this, the generic decoding algorithm burns
 * time renormalizing each row. Still, it's a problem. See notes at
 * the header of utest_optacc() for more info.
 * 
 * Another expected failure mode is when a fwd, bck nat score are close to
 * 0.0; FCompare() can evaluate two close-to-zero numbers as "different"
 * even if their absolute diff is negligible. (I suppose I could fix this.)
 */

/* 
   gcc -g -Wall -msse2 -std=gnu99 -o optacc_utest -I.. -L.. -I../../easel -L../../easel -Dp7OPTACC_TESTDRIVE optacc.c -lhmmer -leasel -lm
   ./optacc_utest
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"
#include "impl_avx.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "41", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,     "50", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,     "45", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,     "20", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  { "--traces",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "dump all tracebacks",                            0 },
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

  /* first round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))            == NULL)  esl_fatal("failed to create null model");

  utest_optacc(go, r, abc, bg, M, L, N);   /* normal sized models */
  utest_optacc(go, r, abc, bg, 1, L, 10);  /* size 1 models       */
  utest_optacc(go, r, abc, bg, M, 1, 10);  /* size 1 sequences    */

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  utest_optacc(go, r, abc, bg, M, L, N);   
  utest_optacc(go, r, abc, bg, 1, L, 10);  
  utest_optacc(go, r, abc, bg, M, 1, 10);  
  /* after every utest_optacc() call, which passes only for the default seed's draws */
  utest_inplace(r, abc, bg, M, L, N,  FALSE);
  utest_inplace(r, abc, bg, 1, L, 10, FALSE);
  utest_inplace(r, abc, bg, M, 1, 10, FALSE);
#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse()) {
    utest_inplace(r, abc, bg, M, L, N,  TRUE);
    utest_inplace(r, abc, bg, 1, L, 10, TRUE);
    utest_inplace(r, abc, bg, M, 1, 10, TRUE);
  }
#endif
  utest_sweep(r, abc, bg, M, 200, N,  FALSE);
  utest_sweep(r, abc, bg, 1, 200, 10, FALSE);
#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse()) {
    utest_sweep(r, abc, bg, M, 200, N,  TRUE);
    utest_sweep(r, abc, bg, 1, 200, 10, TRUE);
  }
#endif

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*p7OPTACC_TESTDRIVE*/
/*------------------ end, test driver ---------------------------*/




/*****************************************************************
 * 6. Example
 *****************************************************************/
#ifdef p7OPTACC_EXAMPLE
/* 
   gcc -g -Wall -o optacc_example -Dp7OPTACC_EXAMPLE -I.. -I../../easel -L.. -L../../easel optacc.c -lhmmer -leasel -lm
   ./optacc_example <hmmfile> <seqfile>
*/
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "impl_avx.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-d",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "dump posterior residue decoding matrix",           0 },
  { "-m",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "dump OA matrix",                                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of optimal accuracy alignment, SSE implementation";

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
  P7_OMX         *ox1     = NULL;
  P7_OMX         *ox2     = NULL;
  P7_GMX         *gx      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  P7_TRACE       *tr      = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  char            errbuf[eslERRBUFSIZE];
  float           fsc, bsc;
  float           accscore;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Read in one sequence */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
  if  (esl_sqio_Read(sqfp, sq) != eslOK) p7_Fail("Failed to read sequence");
  esl_sqfile_Close(sqfp);
 
  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);                 p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);    p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL); /* multihit local: H3 default */
  om = p7_oprofile_Create(gm->M, abc);    p7_oprofile_Convert(gm, om);

  /* Allocations */
  ox1 = p7_omx_Create(gm->M, sq->n, sq->n);
  ox2 = p7_omx_Create(gm->M, sq->n, sq->n);
  gx  = p7_gmx_Create(gm->M, sq->n, sq->n, p7G_NSCELLS);
  tr  = p7_trace_CreateWithPP();
  p7_FLogsumInit();

  /* Run Forward, Backward; do OA fill and trace */
  p7_Forward (sq->dsq, sq->n, om, ox1,      &fsc);
  p7_Backward(sq->dsq, sq->n, om, ox1, ox2, &bsc);
  p7_Decoding(om, ox1, ox2, ox2);                   /* <gx2> is now the posterior decoding matrix */
  p7_OptimalAccuracy(om, ox2, ox1, &accscore);	    /* <gx1> is now the OA matrix */
  p7_OATrace(om, ox1, tr);
  
  if (esl_opt_GetBoolean(go, "-d")) { p7_omx_FDeconvert(ox2, gx);  p7_gmx_Dump(stdout, gx, p7_DEFAULT); }
  if (esl_opt_GetBoolean(go, "-m")) { p7_omx_FDeconvert(ox1, gx);  p7_gmx_Dump(stdout, gx, p7_DEFAULT); }

  p7_trace_Dump(stdout, tr, gm, sq->dsq);

  if (p7_trace_Validate(tr, abc, sq->dsq, errbuf) != eslOK) p7_Die("trace fails validation:\n%s\n", errbuf);

  printf("fwd = %.4f nats\n", fsc);
  printf("bck = %.4f nats\n", bsc);
  printf("acc = %.4f (%.2f%%)\n", accscore, accscore * 100. / (float) sq->n);

  /* Cleanup */
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_omx_Destroy(ox1);
  p7_omx_Destroy(ox2);
  p7_gmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7OPTACC_EXAMPLE*/
/*-------------------- end, example -----------------------------*/




