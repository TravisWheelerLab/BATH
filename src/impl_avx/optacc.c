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
static int p7_OATrace_Dispatcher(const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, P7_TRACE *tr);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_OATrace)(const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, P7_TRACE *tr) = p7_OATrace_Dispatcher;

/* Function:  p7_OATrace_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for optimal accuracy traceback.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> if trace is not empty.
 */
static int
p7_OATrace_Dispatcher(const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, P7_TRACE *tr)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_OATrace = p7_OATrace_avx512; return p7_OATrace_avx512(om, pp, ox, tr); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_OATrace = p7_OATrace_avx;    return p7_OATrace_avx(om, pp, ox, tr); }
#endif
#ifdef eslENABLE_SSE
  p7_OATrace = p7_OATrace_sse;
  return p7_OATrace_sse(om, pp, ox, tr);
#else
  p7_Die("p7_OATrace: no SIMD implementation available");
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7OPTACC_TESTDRIVE
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
      if (p7_OATrace(om, ox2, ox1, tr)                    != eslOK) esl_fatal(msg);
      
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
      if (p7_GOATrace(gm, gx2, gx1, trg)                  != eslOK) esl_fatal(msg);

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
  p7_OATrace(om, ox2, ox1, tr);
  
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




