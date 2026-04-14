/* Viterbi filter dispatcher.
 * Provides the non-suffixed extern functions declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"

/* Function:  p7_ViterbiFilter()
 *
 * Purpose:   Dispatch Viterbi filter to the fastest available ISA path.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if score overflows; <*ret_sc> is <eslINFINITY>.
 * Throws:    <eslEINVAL> if <ox> is too small or profile is not local.
 */
int
p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_ViterbiFilter_avx512(dsq, L, om, ox, ret_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_ViterbiFilter_sse(dsq, L, om, ox, ret_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_ViterbiFilter_sse(dsq, L, om, ox, ret_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_ViterbiFilter_BATH()
 *
 * Purpose:   Dispatch BATH Viterbi filter (score + windows) to the fastest
 *            available ISA path.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if score overflows; <*ret_sc> is <eslINFINITY>.
 * Throws:    <eslEINVAL> if <ox> is too small or profile is not local.
 */
int
p7_ViterbiFilter_BATH(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox,
                      const P7_SCOREDATA *ssvdata, float filtersc, double P,
                      P7_HMM_WINDOWLIST *windowlist, float *ret_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_ViterbiFilter_BATH_avx512(dsq, L, om, ox, ssvdata, filtersc, P, windowlist, ret_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_ViterbiFilter_BATH_sse(dsq, L, om, ox, ssvdata, filtersc, P, windowlist, ret_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_ViterbiFilter_BATH_sse(dsq, L, om, ox, ssvdata, filtersc, P, windowlist, ret_sc);
#else
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef p7VITFILTER_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* ViterbiFilter() unit test
 * 
 * We can check that scores are identical (within machine error) to
 * scores of generic DP with scores rounded the same way.  Do this for
 * a random model of length <M>, for <N> test sequences of length <L>.
 * 
 * We assume that we don't accidentally generate a high-scoring random
 * sequence that overflows ViterbiFilter()'s limited range.
 * 
 */
static void
utest_viterbi_filter(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *ox  = p7_omx_Create(M, 0, 0);
  P7_GMX      *gx  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  float sc1, sc2;

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  p7_profile_SameAsVF(om, gm);	/* round and scale the scores in <gm> the same as in <om> */

#if 0
  p7_oprofile_Dump(stdout, om);              // dumps the optimized profile
  p7_omx_SetDumpMode(stdout, ox, TRUE);      // makes the fast DP algorithms dump their matrices
#endif

  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_ViterbiFilter(dsq, L, om, ox, &sc1);
      p7_GViterbi     (dsq, L, gm, gx, &sc2);

#if 0
      p7_gmx_Dump(stdout, gx, p7_DEFAULT);   // dumps a generic DP matrix
#endif
      
      sc2 /= om->scale_w;
      sc2 -= 3.0;

      if (fabs(sc1-sc2) > 0.001) esl_fatal("viterbi filter unit test failed: scores differ (%.2f, %.2f)", sc1, sc2);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}

/* utest_viterbi_filter_bath()
 *
 * Checks two properties of p7_ViterbiFilter_BATH():
 *
 * 1. Score agreement: the returned Viterbi score must match
 *    p7_ViterbiFilter() to within floating-point rounding (0.001 nats),
 *    because both run the same DP without resetting the matrix.
 *
 * 2. Window validity: every window appended to the windowlist must have
 *    n in [1,L], k in [1,M], and length in [1,M].
 */
static void
utest_viterbi_filter_bath(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char              msg[]  = "ViterbiFilter_BATH unit test failed";
  P7_HMM           *hmm   = NULL;
  P7_PROFILE       *gm    = NULL;
  P7_OPROFILE      *om    = NULL;
  P7_SCOREDATA     *data  = NULL;
  ESL_DSQ          *dsq   = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX           *ox    = p7_omx_Create(M, 0, 0);
  P7_HMM_WINDOWLIST wlist;
  float             sc1, sc2;
  float             nullsc;
  int               w;

  p7_hmmwindow_init(&wlist);
  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  data = p7_hmm_ScoreDataCreate(om, NULL);

  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      /* Reference score from standard ViterbiFilter */
      p7_ViterbiFilter(dsq, L, om, ox, &sc1);

      /* Score + windows from ViterbiFilter_BATH; use nullsc as filtersc */
      p7_bg_NullOne(bg, dsq, L, &nullsc);
      wlist.count = 0;
      p7_ViterbiFilter_BATH(dsq, L, om, ox, data, nullsc, 0.001, &wlist, &sc2);

      /* 1. Scores must agree */
      if (sc1 != eslINFINITY && sc2 != eslINFINITY)
        if (fabs(sc1 - sc2) > 0.001)
          esl_fatal("%s: scores differ: ViterbiFilter=%.4f  ViterbiFilter_BATH=%.4f", msg, sc1, sc2);

      /* 2. Every window must have valid boundaries */
      for (w = 0; w < wlist.count; w++)
        {
          if (wlist.windows[w].n      <  1) esl_fatal("%s: window %d n=%d < 1",           msg, w, wlist.windows[w].n);
          if (wlist.windows[w].n      >  L) esl_fatal("%s: window %d n=%d > L=%d",        msg, w, wlist.windows[w].n, L);
          if (wlist.windows[w].k      <  1) esl_fatal("%s: window %d k=%d < 1",           msg, w, wlist.windows[w].k);
          if (wlist.windows[w].k      >  M) esl_fatal("%s: window %d k=%d > M=%d",        msg, w, wlist.windows[w].k, M);
          if (wlist.windows[w].length <  1) esl_fatal("%s: window %d length=%d < 1",      msg, w, wlist.windows[w].length);
          if (wlist.windows[w].length >  M) esl_fatal("%s: window %d length=%d > M=%d",   msg, w, wlist.windows[w].length, M);
        }
    }

  free(wlist.windows);
  free(dsq);
  p7_hmm_ScoreDataDestroy(data);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(ox);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7VITFILTER_TESTDRIVE*/


/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7VITFILTER_TESTDRIVE
/* 
   gcc -g -Wall -msse2 -std=gnu99 -I.. -L.. -I../../easel -L../../easel -o vitfilter_utest -Dp7VITFILTER_TESTDRIVE vitfilter.c -lhmmer -leasel -lm
   ./vitfilter_utest
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"

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
static char banner[] = "test driver for the SSE implementation";

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

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiFilter() tests, DNA\n");
  utest_viterbi_filter(r, abc, bg, M, L, N);
  utest_viterbi_filter(r, abc, bg, 1, L, 10);
  utest_viterbi_filter(r, abc, bg, M, 1, 10);

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiFilter_BATH() tests, DNA\n");
  utest_viterbi_filter_bath(r, abc, bg, M, L, N);
  utest_viterbi_filter_bath(r, abc, bg, 1, L, 10);
  utest_viterbi_filter_bath(r, abc, bg, M, 1, 10);

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiFilter() tests, protein\n");
  utest_viterbi_filter(r, abc, bg, M, L, N);
  utest_viterbi_filter(r, abc, bg, 1, L, 10);
  utest_viterbi_filter(r, abc, bg, M, 1, 10);

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiFilter_BATH() tests, protein\n");
  utest_viterbi_filter_bath(r, abc, bg, M, L, N);
  utest_viterbi_filter_bath(r, abc, bg, 1, L, 10);
  utest_viterbi_filter_bath(r, abc, bg, M, 1, 10);

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*VITFILTER_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/



/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7VITFILTER_EXAMPLE
/* A minimal example.
   Also useful for debugging on small HMMs and sequences.

   gcc -g -Wall -msse2 -std=gnu99 -I.. -L.. -I../../easel -L../../easel -o vitfilter_example -Dp7VITFILTER_EXAMPLE vitfilter.c -lhmmer -leasel -lm
   ./vitfilter_example <hmmfile> <seqfile>
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
static char banner[] = "example of Viterbi filter algorithm";

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
  float           vfraw, nullsc, vfscore;
  float           graw, gscore;
  double          P, gP;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  /* Read in one sequence */
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
  ox = p7_omx_Create(gm->M, 0, sq->n);
  gx = p7_gmx_Create(gm->M, sq->n, sq->n, p7G_NSCELLS);

  /* Useful to place and compile in for debugging: 
     p7_oprofile_Dump(stdout, om);         dumps the optimized profile
     p7_omx_SetDumpMode(ox, TRUE);         makes the fast DP algorithms dump their matrices
     p7_gmx_Dump(stdout, gx, p7_DEFAULT);  dumps a generic DP matrix
  */

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_ReconfigLength(gm,          sq->n);
      p7_bg_SetLength(bg,            sq->n);
      p7_omx_GrowTo(ox, om->M, 0,    sq->n); 
      p7_gmx_GrowTo(gx, gm->M,       sq->n, sq->n);

      p7_ViterbiFilter  (sq->dsq, sq->n, om, ox, &vfraw);
      p7_bg_NullOne (bg, sq->dsq, sq->n, &nullsc);
      vfscore = (vfraw - nullsc) / eslCONST_LOG2;
      P        = esl_gumbel_surv(vfscore,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);

      p7_GViterbi       (sq->dsq, sq->n, gm, gx, &graw); 
      gscore   = (graw - nullsc) / eslCONST_LOG2;
      gP       = esl_gumbel_surv(gscore,  gm->evparam[p7_VMU],  gm->evparam[p7_VLAMBDA]);

      if (esl_opt_GetBoolean(go, "-1"))
	{
	  printf("%-30s\t%-20s\t%9.2g\t%7.2f\t%9.2g\t%7.2f\n", sq->name, hmm->name, P, vfscore, gP, gscore);
	}
      else if (esl_opt_GetBoolean(go, "-P"))
	{ /* output suitable for direct use in profmark benchmark postprocessors: */
	  printf("%g\t%.2f\t%s\t%s\n", P, vfscore, sq->name, hmm->name);
	}
      else
	{
	  printf("target sequence:      %s\n",        sq->name);
	  printf("vit filter raw score: %.2f nats\n", vfraw);
	  printf("null score:           %.2f nats\n", nullsc);
	  printf("per-seq score:        %.2f bits\n", vfscore);
	  printf("P-value:              %g\n",        P);
	  printf("GViterbi raw score:   %.2f nats\n", graw);
	  printf("GViterbi seq score:   %.2f bits\n", gscore);
	  printf("GViterbi P-value:     %g\n",        gP);
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
#endif /*p7VITFILTER_EXAMPLE*/
/*-------------------- end, example -----------------------------*/


