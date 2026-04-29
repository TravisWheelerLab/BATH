/* Null2 biased composition correction dispatcher.
 * Provides the non-suffixed p7_Null2_ByExpectation() and p7_Null2_ByTrace()
 * declared in impl_avx.h.  Delegates to the fastest available ISA.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Forward declaration of dispatcher */
static int p7_Null2_ByExpectation_Dispatcher(const P7_OPROFILE *om, const P7_OMX *pp, float *null2);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Null2_ByExpectation)(const P7_OPROFILE *om, const P7_OMX *pp, float *null2) = p7_Null2_ByExpectation_Dispatcher;

/* Function:  p7_Null2_ByExpectation_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for null2 estimation (expectation).
 *
 * Returns:   <eslOK> on success.
 */
static int
p7_Null2_ByExpectation_Dispatcher(const P7_OPROFILE *om, const P7_OMX *pp, float *null2)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_Null2_ByExpectation = p7_Null2_ByExpectation_avx512; return p7_Null2_ByExpectation_avx512(om, pp, null2); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Null2_ByExpectation = p7_Null2_ByExpectation_avx;    return p7_Null2_ByExpectation_avx(om, pp, null2); }
#endif
#ifdef eslENABLE_SSE
  p7_Null2_ByExpectation = p7_Null2_ByExpectation_sse;
  return p7_Null2_ByExpectation_sse(om, pp, null2);
#else
  p7_Die("p7_Null2_ByExpectation: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_Null2_ByTrace_Dispatcher(const P7_OPROFILE *om, const P7_TRACE *tr, int zstart, int zend,
                                        P7_OMX *wrk, float *null2);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Null2_ByTrace)(const P7_OPROFILE *om, const P7_TRACE *tr, int zstart, int zend,
                        P7_OMX *wrk, float *null2) = p7_Null2_ByTrace_Dispatcher;

/* Function:  p7_Null2_ByTrace_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for null2 estimation (trace).
 *
 * Returns:   <eslOK> on success.
 */
static int
p7_Null2_ByTrace_Dispatcher(const P7_OPROFILE *om, const P7_TRACE *tr, int zstart, int zend,
                             P7_OMX *wrk, float *null2)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_Null2_ByTrace = p7_Null2_ByTrace_avx512; return p7_Null2_ByTrace_avx512(om, tr, zstart, zend, wrk, null2); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Null2_ByTrace = p7_Null2_ByTrace_avx;    return p7_Null2_ByTrace_avx(om, tr, zstart, zend, wrk, null2); }
#endif
#ifdef eslENABLE_SSE
  p7_Null2_ByTrace = p7_Null2_ByTrace_sse;
  return p7_Null2_ByTrace_sse(om, tr, zstart, zend, wrk, null2);
#else
  p7_Die("p7_Null2_ByTrace: no SIMD implementation available");
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7NULL2_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"

/* compare results to GDecoding(). */
static void
utest_null2_expectation(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N, float tolerance)
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
  P7_GMX      *gpp  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  float       *on2  = malloc(sizeof(float) * abc->Kp);
  float       *gn2  = malloc(sizeof(float) * abc->Kp);
  float fsc1, fsc2;
  float bsc1, bsc2;

  if (!gn2 || !on2) esl_fatal(msg);

  if (p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om) != eslOK) esl_fatal(msg);
  while (N--)
    {
      if (esl_rsq_xfIID(r, bg->f, abc->K, L, dsq) != eslOK) esl_fatal(msg);
      if (p7_Forward       (dsq, L, om, fwd,      &fsc1) != eslOK) esl_fatal(msg);
      if (p7_Backward      (dsq, L, om, fwd, bck, &bsc1) != eslOK) esl_fatal(msg);
      if (p7_Decoding(om, fwd, bck, pp)                  != eslOK) esl_fatal(msg);
      if (p7_Null2_ByExpectation(om, pp, on2)            != eslOK) esl_fatal(msg);
      
      if (p7_GForward (dsq, L, gm, gxf, &fsc2)           != eslOK) esl_fatal(msg);
      if (p7_GBackward(dsq, L, gm, gxb, &bsc2)           != eslOK) esl_fatal(msg);
      if (p7_GDecoding(gm, gxf, gxb, gpp)                != eslOK) esl_fatal(msg);
      if (p7_GNull2_ByExpectation(gm, gpp, gn2)          != eslOK) esl_fatal(msg);

      if (esl_vec_FCompare(gn2, on2, abc->Kp, tolerance) != eslOK) esl_fatal(msg);
    }

  p7_gmx_Destroy(gpp);
  p7_gmx_Destroy(gxf);
  p7_gmx_Destroy(gxb);
  p7_omx_Destroy(pp);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bck);
  free(on2);
  free(gn2);
  free(dsq);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}
#endif /*p7NULL2_TESTDRIVE*/
/*--------------------- end, unit tests -------------------------*/




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7NULL2_TESTDRIVE
/* 
   gcc -g -Wall -msse2 -std=gnu99 -o null2_utest -I.. -L.. -I../../easel -L../../easel -Dp7NULL2_TESTDRIVE null2.c -lhmmer -leasel -lm
   ./null2_utest
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
  { "-t",        eslARG_REAL,  "0.01", NULL, NULL,  NULL,  NULL, NULL, "floating point comparison tolerance",            0 },
  { "-L",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,     "72", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,     "10", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE implementation of null2 model";

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

  utest_null2_expectation(r, abc, bg, M, L, N, tol);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);
  return eslOK;
}
#endif /*p7NULL2_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/






