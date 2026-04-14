/* MSV filter dispatcher.
 * Provides the non-suffixed extern functions declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"

/* Function:  p7_MSVFilter()
 *
 * Purpose:   Dispatch MSV filter to the fastest available ISA path.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> on score overflow (high-scoring hit).
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_MSVFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_MSVFilter_avx512(dsq, L, om, ox, ret_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_MSVFilter_sse(dsq, L, om, ox, ret_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_MSVFilter_sse(dsq, L, om, ox, ret_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_SSVFilter_BATH()
 *
 * Purpose:   Dispatch BATH SSV window finder to the fastest available ISA path.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_SSVFilter_BATH(const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox,
                  const P7_SCOREDATA *ssvdata, P7_BG *bg, double P,
                  P7_HMM_WINDOWLIST *windowlist)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_SSVFilter_BATH_avx512(dsq, L, om, ox, ssvdata, bg, P, windowlist);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_SSVFilter_BATH_sse(dsq, L, om, ox, ssvdata, bg, P, windowlist);
#endif
#ifdef eslENABLE_SSE
  return p7_SSVFilter_BATH_sse(dsq, L, om, ox, ssvdata, bg, P, windowlist);
#else
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7MSVFILTER_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* 
 * We can check that scores are identical (within machine error) to
 * scores of generic DP with scores rounded the same way.  Do this for
 * a random model of length <M>, for <N> test sequences of length <L>.
 * 
 * We assume that we don't accidentally generate a high-scoring random
 * sequence that overflows MSVFilter()'s limited range.
 * 
 */
static void
utest_msv_filter(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *ox  = p7_omx_Create(M, 0, 0);
  P7_GMX      *gx  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  float sc1, sc2;

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  p7_profile_SameAsMF(om, gm);
#if 0
  p7_oprofile_Dump(stdout, om);              /* dumps the optimized profile */
  p7_omx_SetDumpMode(stdout, ox, TRUE);      /* makes the fast DP algorithms dump their matrices */
#endif

  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_MSVFilter(dsq, L, om, ox, &sc1);
      p7_GViterbi (dsq, L, gm, gx, &sc2);
#if 0
      p7_gmx_Dump(stdout, gx, p7_DEFAULT);   /* dumps a generic DP matrix */
#endif

      sc2 = sc2 / om->scale_b - 3.0f;
      if (fabs(sc1-sc2) > 0.001) esl_fatal("msv filter unit test failed: scores differ (%.2f, %.2f)", sc1, sc2);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7MSVFILTER_TESTDRIVE*/
/*-------------------- end, unit tests --------------------------*/




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7MSVFILTER_TESTDRIVE
/* 
   gcc -g -Wall -msse2 -std=gnu99 -I.. -L.. -I../../easel -L../../easel -o msvfilter_utest -Dp7MSVFILTER_TESTDRIVE msvfilter.c -lhmmer -leasel -lm
   ./msvfilter_utest
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
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the SSE MSVFilter() implementation";

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

  /* First round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))            == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("MSVFilter() tests, DNA\n");
  utest_msv_filter(r, abc, bg, M, L, N);   /* normal sized models */
  utest_msv_filter(r, abc, bg, 1, L, 10);  /* size 1 models       */
  utest_msv_filter(r, abc, bg, M, 1, 10);  /* size 1 sequences    */

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("MSVFilter() tests, protein\n");
  utest_msv_filter(r, abc, bg, M, L, N);   
  utest_msv_filter(r, abc, bg, 1, L, 10);  
  utest_msv_filter(r, abc, bg, M, 1, 10);  

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*VITFILTER_TESTDRIVE*/



/*****************************************************************
 * 5. Example
 *****************************************************************/

#ifdef p7MSVFILTER_EXAMPLE
/* A minimal example.
   Also useful for debugging on small HMMs and sequences.

   gcc -g -Wall -msse2 -std=gnu99 -I.. -L.. -I../../easel -L../../easel -o msvfilter_example -Dp7MSVFILTER_EXAMPLE msvfilter.c -lhmmer -leasel -lm
   ./msvfilter_example <hmmfile> <seqfile>
 */ 
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "impl_avx.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "output in one line awkable format",                0 },
  { "-P",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "output in profmark format",                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of MSV filter algorithm";

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
  P7_OMX         *ox      = NULL;
  P7_GMX         *gx      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           msvraw, nullsc, msvscore;
  float           graw, gscore;
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
  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  /* allocate DP matrices, both a generic and an optimized one */
  ox = p7_omx_Create(gm->M, 0, 0); /* one row version */
  gx = p7_gmx_Create(gm->M, sq->n, sq->n, p7G_NSCELLS);

  /* Useful to place and compile in for debugging: 
     p7_oprofile_Dump(stdout, om);              dumps the optimized profile
     p7_omx_SetDumpMode(stdout, ox, TRUE);      makes the fast DP algorithms dump their matrices
     p7_gmx_Dump(stdout, gx, p7_DEFAULT);       dumps a generic DP matrix
     p7_oprofile_SameMSV(om, gm);
  */
  //p7_oprofile_Dump(stdout, om);
  //p7_omx_SetDumpMode(stdout, ox, TRUE);    

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_ReconfigLength(gm,          sq->n);
      p7_bg_SetLength(bg,            sq->n);
      p7_omx_GrowTo(ox, om->M, 0,    sq->n); 
      p7_gmx_GrowTo(gx, gm->M,       sq->n, sq->n);

      p7_MSVFilter   (sq->dsq, sq->n, om, ox, &msvraw);  
      p7_bg_NullOne  (bg, sq->dsq, sq->n, &nullsc);
      msvscore = (msvraw - nullsc) / eslCONST_LOG2;
      P        = esl_gumbel_surv(msvscore,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);

      p7_GMSV(sq->dsq, sq->n, gm, gx, 2.0, &graw);
      gscore   = (graw - nullsc) / eslCONST_LOG2;
      gP       = esl_gumbel_surv(gscore,  gm->evparam[p7_MMU],  gm->evparam[p7_MLAMBDA]);

      if (esl_opt_GetBoolean(go, "-1"))
	{
	  printf("%-30s  %-20s  %9.2g  %7.2f  %9.2g  %7.2f\n", sq->name, hmm->name, P, msvscore, gP, gscore);
	}
      else if (esl_opt_GetBoolean(go, "-P"))
	{ /* output suitable for direct use in profmark benchmark postprocessors: */
	  printf("%g  %.2f  %s  %s\n", P, msvscore, sq->name, hmm->name);
	}
      else
	{
	  printf("target sequence:      %s\n",        sq->name);
	  printf("msv filter raw score: %.2f nats\n", msvraw);
	  printf("null score:           %.2f nats\n", nullsc);
	  printf("per-seq score:        %.2f bits\n", msvscore);
	  printf("P-value:              %g\n",        P);
	  printf("GMSV raw score:       %.2f nats\n", graw);
	  printf("GSMV per-seq score:   %.2f bits\n", gscore);
	  printf("GSMV P-value:         %g\n",        gP);
	}
      
      esl_sq_Reuse(sq);
    }

  /* cleanup */
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_omx_Destroy(ox);
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
#endif /*p7MSVFILTER_EXAMPLE*/
/*---------------------- end, example ---------------------------*/




