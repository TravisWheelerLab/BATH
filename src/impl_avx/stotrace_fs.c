/* Frameshift stochastic traceback dispatcher.
 * Provides the non-suffixed p7_StochasticTrace_Frameshift() declared in
 * impl_avx.h.  Delegates to the fastest available ISA implementation at
 * runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Forward declaration of dispatcher */
static int p7_StochasticTrace_Frameshift_Dispatcher(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L,
                                                     const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, P7_TRACE *tr);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_StochasticTrace_Frameshift)(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L,
                                     const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, P7_TRACE *tr) = p7_StochasticTrace_Frameshift_Dispatcher;

/* Function:  p7_StochasticTrace_Frameshift_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for frameshift stochastic traceback.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on various internal failures.
 */
static int
p7_StochasticTrace_Frameshift_Dispatcher(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L,
                                          const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, P7_TRACE *tr)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_StochasticTrace_Frameshift = p7_StochasticTrace_Frameshift_avx512; return p7_StochasticTrace_Frameshift_avx512(rng, dsq, L, om_fs, ox, tr); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_StochasticTrace_Frameshift = p7_StochasticTrace_Frameshift_avx;    return p7_StochasticTrace_Frameshift_avx(rng, dsq, L, om_fs, ox, tr); }
#endif
#ifdef eslENABLE_SSE
  p7_StochasticTrace_Frameshift = p7_StochasticTrace_Frameshift_sse;
  return p7_StochasticTrace_Frameshift_sse(rng, dsq, L, om_fs, ox, tr);
#else
  p7_Die("p7_StochasticTrace_Frameshift: no SIMD implementation available");
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7STOTRACE_FS_TESTDRIVE
#include "esl_getopts.h"
#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_stotrace_fs()
 * Tests:
 *   1. Generic and SSE forward scores agree to tolerance.
 *   2. Each SSE stochastic trace is structurally valid: N >= 3, starts with
 *      S, ends with T.
 *   3. Both the generic and SSE stochastic trace functions complete without
 *      error on every sequence.
 */
static void
utest_stotrace_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
                  ESL_GENCODE *gcode, P7_BG *bgAA,
                  P7_CODONTABLE *codon_table, int M, int N)
{
  char           *msg    = "utest_stotrace_fs failed";
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_fs5 = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE *om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq    = NULL;
  P7_TRACE       *tr     = p7_trace_fs_Create();
  P7_IVX         *iv5    = p7_ivx_Create(M, p7P_5CODONS);
  P7_OIVX        *ov5    = p7_oivx_Create(M, p7P_5CODONS);
  P7_GMX         *gxf    = p7_gmx_Create(M, M, M, p7G_NSCELLS_FS);
  P7_OMX         *fwd    = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS_FS);
  float           gfsc, ofsc;
  float           tolerance;
  int             i, j, curr_L, idx;

  tolerance = (p7_FLogsumError(-0.4, -0.5) > 0.0001) ? 1.0f : 0.001f;

  if (!gm || !gm_fs5 || !om_fs5 || !sq || !tr || !iv5 || !gxf || !fwd) esl_fatal(msg);

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig(hmm, bgAA, gm, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, M, p7_LOCAL);
  p7_fs_oprofile_Convert(gm_fs5, om_fs5);

  while (N--)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, NULL);
      curr_L = sq->n * 3;

      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (curr_L + 2))) == NULL) esl_fatal("malloc failed");
      j = 1;
      for (i = 1; i <= (int)sq->n; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }

      p7_fs_oprofile_ReconfigLength(om_fs5, sq->n);
      p7_fs_ReconfigLength(gm_fs5, sq->n);

      /* Generic forward */
      p7_gmx_GrowTo(gxf, M, curr_L, curr_L);
      p7_ivx_GrowTo(iv5, M, p7P_5CODONS);
      if (p7_GForward_Frameshift(dsq, curr_L, gm_fs5, gxf, iv5, &gfsc) != eslOK) esl_fatal(msg);

      /* SSE forward */
      p7_omx_GrowTo_dpf(fwd, M, curr_L, curr_L);
      { int s = p7_Forward_Frameshift(dsq, curr_L, om_fs5, fwd, ov5, &ofsc);
        if (s == eslERANGE) continue;
        if (s != eslOK)     esl_fatal(msg); }

      /* Forward scores must agree */
      if (esl_FCompare_old(gfsc, ofsc, tolerance) != eslOK)
        esl_fatal("%s: generic fwd %.4f != SSE fwd %.4f (tol %.4f)", msg, gfsc, ofsc, tolerance);

      /* Sample SSE stochastic traces and validate each with p7_trace_fs_Validate */
      for (idx = 0; idx < 10; idx++)
        {
          char errbuf[eslERRBUFSIZE];
          if (p7_StochasticTrace_Frameshift(r, dsq, curr_L, om_fs5, fwd, tr) != eslOK)
            esl_fatal("%s: SSE stochastic trace failed", msg);
          if (p7_trace_fs_Validate(tr, abcDNA, dsq, errbuf) != eslOK)
            esl_fatal("%s: SSE trace invalid: %s", msg, errbuf);
          p7_trace_Reuse(tr);
        }

      /* Verify generic stochastic trace also runs cleanly */
      if (p7_GStochasticTrace_Frameshift(r, dsq, curr_L, gm_fs5, gxf, tr) != eslOK)
        esl_fatal("%s: generic stochastic trace failed", msg);
      p7_trace_Reuse(tr);
    }

  free(dsq);
  p7_trace_fs_Destroy(tr);
  p7_omx_Destroy(fwd);
  p7_ivx_Destroy(iv5);
  p7_gmx_Destroy(gxf);
  p7_oivx_Destroy(ov5);
  p7_hmm_Destroy(hmm);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_fs5);
  p7_fs_oprofile_Destroy(om_fs5);
}
#endif /*p7STOTRACE_FS_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/



/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7STOTRACE_FS_TESTDRIVE
/*
  gcc -g -Wall -msse2 -std=gnu99 -o stotrace_fs_utest -I.. -L.. -I../../easel -L../../easel -Dp7STOTRACE_FS_TESTDRIVE stotrace_fs.c -lhmmer -leasel -lm
  ./stotrace_fs_utest
*/
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_avx.h"

static ESL_OPTIONS options[] = {
  /* name    type           default env range toggles reqs incomp help                                       docgroup */
  { "-h",  eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",  eslARG_INT,     "42", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M",  eslARG_INT,     "72", NULL, NULL, NULL, NULL, NULL, "size of random models to sample",                0 },
  { "-N",  eslARG_INT,     "10", NULL, NULL, NULL, NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE frameshift stochastic traceback";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA  = esl_alphabet_Create(eslAMINO);
  ESL_ALPHABET   *abcDNA = esl_alphabet_Create(eslDNA);
  ESL_GENCODE    *gcode  = esl_gencode_Create(abcDNA, abcAA);
  P7_CODONTABLE  *ct     = p7_codontable_Create(gcode);
  P7_BG          *bgAA   = p7_bg_Create(abcAA);
  int             M      = esl_opt_GetInteger(go, "-M");
  int             N      = esl_opt_GetInteger(go, "-N");

  p7_FLogsumInit();
  impl_Init();

  utest_stotrace_fs(r, abcAA, abcDNA, gcode, bgAA, ct, M, N);

  p7_bg_Destroy(bgAA);
  p7_codontable_Destroy(ct);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcAA);
  esl_alphabet_Destroy(abcDNA);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7STOTRACE_FS_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/


/*****************************************************************
 * 6. Example.
 *****************************************************************/

/*------------------------ end, example -------------------------*/


