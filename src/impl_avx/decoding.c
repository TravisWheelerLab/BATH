/* Posterior decoding dispatcher.
 * Provides the non-suffixed p7_Decoding() and p7_DomainDecoding() declared
 * in impl_avx.h.  Delegates to the fastest available ISA implementation
 * at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Forward declaration of dispatcher */
static int p7_Decoding_Dispatcher(const P7_OPROFILE *om, const P7_OMX *oxf, P7_OMX *oxb, P7_OMX *pp);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Decoding)(const P7_OPROFILE *om, const P7_OMX *oxf, P7_OMX *oxb, P7_OMX *pp) = p7_Decoding_Dispatcher;

/* Function:  p7_Decoding_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for posterior decoding of residue assignment.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if numeric range is exceeded.
 */
static int
p7_Decoding_Dispatcher(const P7_OPROFILE *om, const P7_OMX *oxf, P7_OMX *oxb, P7_OMX *pp)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_Decoding = p7_Decoding_avx512; return p7_Decoding_avx512(om, oxf, oxb, pp); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Decoding = p7_Decoding_avx;    return p7_Decoding_avx(om, oxf, oxb, pp); }
#endif
#ifdef eslENABLE_SSE
  p7_Decoding = p7_Decoding_sse;
  return p7_Decoding_sse(om, oxf, oxb, pp);
#else
  p7_Die("p7_Decoding: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_DomainDecoding_Dispatcher(const P7_OPROFILE *om, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_DomainDecoding)(const P7_OPROFILE *om, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef) = p7_DomainDecoding_Dispatcher;

/* Function:  p7_DomainDecoding_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for posterior decoding of domain location.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> on numeric overflow.
 */
static int
p7_DomainDecoding_Dispatcher(const P7_OPROFILE *om, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_DomainDecoding = p7_DomainDecoding_avx512; return p7_DomainDecoding_avx512(om, oxf, oxb, ddef); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_DomainDecoding = p7_DomainDecoding_avx;    return p7_DomainDecoding_avx(om, oxf, oxb, ddef); }
#endif
#ifdef eslENABLE_SSE
  p7_DomainDecoding = p7_DomainDecoding_sse;
  return p7_DomainDecoding_sse(om, oxf, oxb, ddef);
#else
  p7_Die("p7_DomainDecoding: no SIMD implementation available");
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
#ifdef p7DECODING_BENCHMARK
/*
   gcc -g -O3 -mavx2 -std=gnu99 -o benchmark-decoding -I.. -L.. -I../../easel -L../../easel -Dp7DECODING_BENCHMARK decoding.c -lhmmer -leasel -lm
   icc  -O3 -static -o benchmark-decoding -I.. -L.. -I../../easel -L../../easel -Dp7DECODING_BENCHMARK decoding.c -lhmmer -leasel -lm
   ./benchmark-decoding <hmmfile>
                    RRM_1 (M=72)       Caudal_act (M=136)     SMC_N (M=1151)
                 -----------------    ------------------     ---------------
   21 Aug 08      3.52u (409 Mc/s)     15.36u (177 Mc/s)     318.78u (72.2 Mc/s)

   The length dependency probably indicates L1 cache missing; because we're
   manipulating 3 matrices at the same time, we can't fit the calculation
   in cache.
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
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for posterior residue decoding, AVX version";

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
  P7_OMX         *fwd     = NULL;
  P7_OMX         *bck     = NULL;
  P7_OMX         *pp      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           fsc, bsc;
  double          Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);                 p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);    p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);    p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);

  fwd = p7_omx_Create(gm->M, L, L);
  bck = p7_omx_Create(gm->M, L, L);
  pp  = p7_omx_Create(gm->M, L, L);

  esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  p7_Forward (dsq, L, om, fwd,      &fsc);
  p7_Backward(dsq, L, om, fwd, bck, &bsc);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    p7_Decoding(om, fwd, bck, pp);
  esl_stopwatch_Stop(w);

  Mcs = (double) N * (double) L * (double) gm->M * 1e-6 / (double) w->user;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bck);
  p7_omx_Destroy(pp);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7DECODING_BENCHMARK*/
/*------------------ end, benchmark driver ----------------------*/


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7DECODING_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* compare results to GDecoding(). */
static void
utest_decoding(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N, float tolerance)
{
  char        *msg  = "decoding unit test failed";
  P7_HMM      *hmm  = NULL;
  P7_PROFILE  *gm   = NULL;
  P7_OPROFILE *om   = NULL;
  ESL_DSQ     *dsq  = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *fwd  = p7_omx_Create(M, L, L);
  P7_OMX      *bck  = p7_omx_Create(M, L, L);
  P7_OMX      *pp   = p7_omx_Create(M, L, L);
  P7_GMX      *gxf  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX      *gxb  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX      *gxp1 = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX      *gxp2 = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  float fsc1, fsc2;
  float bsc1, bsc2;

  if (p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om) != eslOK) esl_fatal(msg);
  while (N--)
    {
      if (esl_rsq_xfIID(r, bg->f, abc->K, L, dsq) != eslOK) esl_fatal(msg);
      if (p7_Forward       (dsq, L, om, fwd,      &fsc1) != eslOK) esl_fatal(msg);
      if (p7_Backward      (dsq, L, om, fwd, bck, &bsc1) != eslOK) esl_fatal(msg);
      if (p7_Decoding(om, fwd, bck, pp)                  != eslOK) esl_fatal(msg);
      if (p7_omx_FDeconvert(pp, gxp1)                    != eslOK) esl_fatal(msg);
      
      if (p7_GForward (dsq, L, gm, gxf, &fsc2)           != eslOK) esl_fatal(msg);
      if (p7_GBackward(dsq, L, gm, gxb, &bsc2)           != eslOK) esl_fatal(msg);
      if (p7_GDecoding(gm, gxf, gxb, gxp2)               != eslOK) esl_fatal(msg);

      // p7_gmx_Dump(stdout, gxp1, p7_DEFAULT);
      // p7_gmx_Dump(stdout, gxp2, p7_DEFAULT);

      if (p7_gmx_Compare(gxp1, gxp2, tolerance)          != eslOK) esl_fatal(msg);
    }

  p7_gmx_Destroy(gxp1);
  p7_gmx_Destroy(gxp2);
  p7_gmx_Destroy(gxf);
  p7_gmx_Destroy(gxb);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bck);
  p7_omx_Destroy(pp);
  free(dsq);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}
#endif /*p7DECODING_TESTDRIVE*/
/*--------------------- end, unit tests -------------------------*/




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7DECODING_TESTDRIVE

/*
  gcc -o decoding_utest -msse2 -g -Wall -I.. -L.. -I../../easel -L../../easel -Dp7DECODING_TESTDRIVE decoding.c -lhmmer -leasel -lm
  ./decoding_utest
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "hmmer.h"
#include "impl_avx.h"

static ESL_OPTIONS options[] = {
  /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  { "-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",                  0},
  { "-s",  eslARG_INT,     "42",  NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",        0 },
  { "-t",  eslARG_REAL,  "0.01",  NULL, NULL, NULL, NULL, NULL, "floating point comparison tolerance",  0 },
  { "-L",  eslARG_INT,     "40",  NULL, NULL, NULL, NULL, NULL, "length of sampled sequences",          0 },
  { "-M",  eslARG_INT,     "40",  NULL, NULL, NULL, NULL, NULL, "length of sampled test profile",       0 },
  { "-N",  eslARG_INT,     "10",  NULL, NULL, NULL, NULL, NULL, "number of sampled test sequences",     0 },
  { 0,0,0,0,0,0,0,0,0,0},
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE posterior decoding";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg   = p7_bg_Create(abc);
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");
  float           tol  = esl_opt_GetReal   (go, "-t");
  
  p7_FLogsumInit();
   impl_Init();

  utest_decoding(r, abc, bg, M, L, N, tol);
  
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);
  return eslOK;
}
#endif /*p7DECODING_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/


