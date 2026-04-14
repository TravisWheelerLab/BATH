/* Viterbi full-matrix dispatcher.
 * Provides the non-suffixed p7_Viterbi() and p7_Viterbi_Trace() declared in
 * impl_avx.h.  Delegates to the fastest available ISA implementation.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Viterbi()
 *
 * Purpose:   Dispatch full-matrix Viterbi to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_Viterbi(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Viterbi_avx512(dsq, L, om, ox, ret_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Viterbi_sse(dsq, L, om, ox, ret_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_Viterbi_sse(dsq, L, om, ox, ret_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_Viterbi_Trace()
 *
 * Purpose:   Dispatch Viterbi traceback to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if an impossible state is reached.
 */
int
p7_Viterbi_Trace(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *ox,
                 P7_TRACE *tr)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Viterbi_Trace_avx512(dsq, L, om, ox, tr);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Viterbi_Trace_sse(dsq, L, om, ox, tr);
#endif
#ifdef eslENABLE_SSE
  return p7_Viterbi_Trace_sse(dsq, L, om, ox, tr);
#else
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef p7VITERBI_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_viterbi()
 *
 * Compare p7_Viterbi() scores to p7_GViterbi() scores.  Differences
 * should be negligible floating-point roundoff only (< 0.001 nats).
 * The optimized profile must be in lspace (p7_oprofile_Logify()).
 */
static void
utest_viterbi(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char        *msg = "SSE p7_Viterbi() unit test failed";
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *ox  = p7_omx_Create(M, L, L);
  P7_GMX      *gx  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_TRACE    *tr  = p7_trace_fs_Create();
  P7_TRACE    *trg = p7_trace_fs_Create();
  char         errbuf[eslERRBUFSIZE];
  float        sc1, sc2, tsc1, tsc2;

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  p7_oprofile_Logify(om);
  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_Viterbi (dsq, L, om, ox, &sc1);
      p7_GViterbi(dsq, L, gm, gx, &sc2);
      if (fabs(sc1 - sc2) > 0.001) esl_fatal("%s: Viterbi scores disagree: %.4f vs %.4f", msg, sc1, sc2);

      p7_Viterbi_Trace(dsq, L, om, ox, tr);
      p7_GTrace(dsq, L, gm, gx, trg);

      if (p7_trace_Validate(tr,  abc, dsq, errbuf) != eslOK) esl_fatal("%s: SSE trace invalid: %s",     msg, errbuf);
      if (p7_trace_Validate(trg, abc, dsq, errbuf) != eslOK) esl_fatal("%s: generic trace invalid: %s", msg, errbuf);

      p7_trace_Score(tr,  dsq, gm, &tsc1);
      p7_trace_Score(trg, dsq, gm, &tsc2);
      if (fabs(tsc1 - tsc2) > 0.001) esl_fatal("%s: trace scores disagree: %.4f vs %.4f", msg, tsc1, tsc2);

      p7_trace_Reuse(tr);
      p7_trace_Reuse(trg);
    }

  free(dsq);
  p7_trace_fs_Destroy(tr);
  p7_trace_fs_Destroy(trg);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7VITERBI_TESTDRIVE*/
/*-------------------- end, unit tests --------------------------*/

/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7VITERBI_TESTDRIVE
/*
   gcc -g -Wall -msse2 -std=gnu99 -I.. -L.. -I../../easel -L../../easel -o viterbi_utest -Dp7VITERBI_TESTDRIVE viterbi.c -lhmmer -leasel -lm
   ./viterbi_utest
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
static char banner[] = "test driver for SSE p7_Viterbi() full-matrix implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r   = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc = NULL;
  P7_BG          *bg  = NULL;
  int             M   = esl_opt_GetInteger(go, "-M");
  int             L   = esl_opt_GetInteger(go, "-L");
  int             N   = esl_opt_GetInteger(go, "-N");

  impl_Init();

  /* First round of tests for DNA alphabets. */
  if ((abc = esl_alphabet_Create(eslDNA))   == NULL) esl_fatal("failed to create alphabet");
  if ((bg  = p7_bg_Create(abc))             == NULL) esl_fatal("failed to create null model");

  utest_viterbi(r, abc, bg, M, L, N);   /* normal sized models  */
  utest_viterbi(r, abc, bg, 1, L, 10);  /* size 1 models        */
  utest_viterbi(r, abc, bg, M, 1, 10);  /* size 1 sequences     */

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets. */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL) esl_fatal("failed to create alphabet");
  if ((bg  = p7_bg_Create(abc))             == NULL) esl_fatal("failed to create null model");

  utest_viterbi(r, abc, bg, M, L, N);
  utest_viterbi(r, abc, bg, 1, L, 10);
  utest_viterbi(r, abc, bg, M, 1, 10);

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*p7VITERBI_TESTDRIVE*/
/*--------------------- end, test driver ------------------------*/


