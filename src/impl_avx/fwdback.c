/* Forward/Backward dispatcher.
 * Provides the non-suffixed p7_Forward, p7_ForwardParser, p7_Backward,
 * and p7_BackwardParser declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"

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
/*****************************************************************
 * 5. Unit tests.
 *****************************************************************/
#ifdef p7FWDBACK_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

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

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  utest_fwdback(r, abc, bg, M, L, N);   
  utest_fwdback(r, abc, bg, 1, L, 10);  
  utest_fwdback(r, abc, bg, M, 1, 10);  

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

