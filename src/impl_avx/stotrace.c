/* Stochastic traceback dispatcher.
 * Provides the non-suffixed p7_StochasticTrace() declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Forward declaration of dispatcher */
static int p7_StochasticTrace_Dispatcher(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,
                                         const P7_OMX *ox, P7_TRACE *tr);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_StochasticTrace)(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,
                          const P7_OMX *ox, P7_TRACE *tr) = p7_StochasticTrace_Dispatcher;

/* Function:  p7_StochasticTrace_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for stochastic traceback.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> on various internal failures.
 */
static int
p7_StochasticTrace_Dispatcher(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,
                               const P7_OMX *ox, P7_TRACE *tr)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_StochasticTrace = p7_StochasticTrace_avx512; return p7_StochasticTrace_avx512(rng, dsq, L, om, ox, tr); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_StochasticTrace = p7_StochasticTrace_avx;    return p7_StochasticTrace_avx(rng, dsq, L, om, ox, tr); }
#endif
#ifdef eslENABLE_SSE
  p7_StochasticTrace = p7_StochasticTrace_sse;
  return p7_StochasticTrace_sse(rng, dsq, L, om, ox, tr);
#else
  p7_Die("p7_StochasticTrace: no SIMD implementation available");
  return eslENORESULT;
#endif
}

/* Forward declaration of dispatcher */
static int p7_StochasticTraceChk_Dispatcher(ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,
                                            P7_CHKOMX *cx, P7_TRACE **tr);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_StochasticTraceChk)(ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,
                             P7_CHKOMX *cx, P7_TRACE **tr) = p7_StochasticTraceChk_Dispatcher;

/* Function:  p7_StochasticTraceChk_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for sampling tracebacks from a
 *            checkpointed Forward matrix. See p7_StochasticTraceChk_avx().
 */
static int
p7_StochasticTraceChk_Dispatcher(ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,
                                 P7_CHKOMX *cx, P7_TRACE **tr)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_StochasticTraceChk = p7_StochasticTraceChk_avx; return p7_StochasticTraceChk_avx(rng, ntr, dsq, L, om, cx, tr); }
#endif
#ifdef eslENABLE_SSE
  p7_StochasticTraceChk = p7_StochasticTraceChk_sse;
  return p7_StochasticTraceChk_sse(rng, ntr, dsq, L, om, cx, tr);
#else
  p7_Die("p7_StochasticTraceChk: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7STOTRACE_TESTDRIVE
#include "esl_getopts.h"
#include "esl_randomseq.h"
#include "esl_sq.h"
#include "esl_vectorops.h"
#include <string.h>

/* tests: 
 *   1. each sampled trace must validate.
 *   2. each trace must be <= viterbi trace score
 *   3. in a large # of traces, one is "equal" to the viterbi trace score.
 *      (this of course is stochastic; but it's true for the particular
 *       choice of RNG seed used in tests here.)
 */
static void
utest_stotrace(ESL_GETOPTS *go, ESL_RANDOMNESS *rng, ESL_ALPHABET *abc, P7_PROFILE *gm, P7_OPROFILE *om, ESL_DSQ *dsq, int L, int ntrace)
{
  P7_GMX   *gx  = NULL;
  P7_OMX   *ox  = NULL;
  P7_TRACE *vtr = NULL;
  P7_TRACE *tr  = NULL;
  char      errbuf[eslERRBUFSIZE];
  int       idx;
  float     maxsc = -eslINFINITY;
  float     vsc, sc;

  if ((gx     = p7_gmx_Create(gm->M, L, L, p7G_NSCELLS)) == NULL)  esl_fatal("generic DP matrix creation failed");
  if ((ox     = p7_omx_Create(gm->M, L, L))              == NULL)  esl_fatal("optimized DP matrix create failed");
  if ((tr     = p7_trace_Create())                       == NULL)  esl_fatal("trace creation failed");
  if ((vtr    = p7_trace_Create())                       == NULL)  esl_fatal("trace creation failed");

  if (p7_GViterbi(dsq, L, gm, gx, &vsc) != eslOK) esl_fatal("viterbi failed");
  if (p7_GTrace  (dsq, L, gm, gx, vtr)  != eslOK) esl_fatal("viterbi trace failed");
  if (p7_Forward (dsq, L, om, ox, NULL) != eslOK) esl_fatal("forward failed");

  for (idx = 0; idx < ntrace; idx++)
    {
      if (p7_StochasticTrace(rng, dsq, L, om, ox, tr) != eslOK) esl_fatal("stochastic trace failed");
      if (p7_trace_Validate(tr, abc, dsq, errbuf)     != eslOK) esl_fatal("trace invalid:\n%s", errbuf);
      if (p7_trace_Score(tr, dsq, gm, &sc)            != eslOK) esl_fatal("trace scoring failed"); 

      maxsc = ESL_MAX(sc, maxsc);
      if (sc > vsc + 0.001){	/* need a little tolerance of floating point math here  */
	//p7_trace_Dump(stdout, vtr, gm, dsq);
	//p7_trace_Dump(stdout, tr,  gm, dsq);
	esl_fatal("sampled trace has score > optimal Viterbi path; not possible (%f > %f)", sc, vsc);
      }
      p7_trace_Reuse(tr);
    }
  if (esl_FCompare_old(maxsc, vsc, 0.1) != eslOK) esl_fatal("stochastic trace failed to sample the Viterbi path");
  
  p7_trace_Destroy(tr);
  p7_trace_Destroy(vtr);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
}
/* utest_stotrace_chk()
 * Stochastic traces from a checkpointed Forward matrix against
 * p7_StochasticTrace() on the full Forward matrix. For random sequences
 * and the consensus of a model with sharpened match emissions and M->M
 * transitions between random flanks, <ntr> traces sampled together by
 * p7_StochasticTraceChk() must be identical to those p7_StochasticTrace()
 * samples one at a time with the same per-trace random number
 * generators, with the checkpointed matrix in one block and in blocks of
 * 1..8 rows. With <use_sse>, the _sse versions are called directly, so
 * an AVX CPU also tests the SSE engine.
 */
static void
utest_stotrace_chk(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N, int use_sse)
{
  char            *msg  = "checkpointed stochastic trace unit test failed";
  int              ntr  = 20;
  P7_HMM          *hmm  = NULL;
  P7_PROFILE      *gm   = NULL;
  P7_OPROFILE     *om   = NULL;
  P7_OMX          *ox   = NULL;
  P7_CHKOMX       *cx   = NULL;
  ESL_SQ          *sq   = esl_sq_CreateDigital(abc);
  ESL_DSQ         *dsq  = malloc(sizeof(ESL_DSQ) * (L+2));
  ESL_RANDOMNESS  *rng1 = esl_randomness_Create(1);
  ESL_RANDOMNESS **rng  = malloc(sizeof(ESL_RANDOMNESS *) * ntr);
  P7_TRACE        *tr1  = p7_trace_Create();
  P7_TRACE       **tr2  = malloc(sizeof(P7_TRACE *) * ntr);
  uint32_t        *seed = malloc(sizeof(uint32_t) * ntr);
  int              ntested = 0;
  int              n, Ln, blockR, t, z, k, x;

  for (t = 0; t < ntr; t++) { rng[t] = esl_randomness_Create(1); tr2[t] = p7_trace_Create(); }

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
      om = p7_oprofile_Create_sse(M, abc);  p7_oprofile_Convert_sse(gm, om);
      ox = p7_omx_Create_sse(M, L, L);
      cx = p7_chkomx_Create_sse(M, L, L+1, 1, p7X_NSCELLS, p7_CHK_FORWARD);
    }
  else
    {
      om = p7_oprofile_Create(M, abc);  p7_oprofile_Convert(gm, om);
      ox = p7_omx_Create(M, L, L);
      cx = p7_chkomx_Create(M, L, L+1, 1, p7X_NSCELLS, p7_CHK_FORWARD);
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
      blockR = (n % 3 == 2) ? Ln+1 : 1 + esl_rnd_Roll(r, 8);
      if (use_sse)
        {
          p7_oprofile_ReconfigLength_sse(om, Ln);
          if (p7_Forward_sse(dsq, Ln, om, ox, NULL)       != eslOK) continue;
          if (p7_chkomx_GrowTo_sse(cx, M, Ln, blockR)     != eslOK) esl_fatal(msg);
          if (p7_ForwardChk_sse(dsq, Ln, om, cx, NULL)    != eslOK) esl_fatal(msg);
        }
      else
        {
          p7_oprofile_ReconfigLength(om, Ln);
          if (p7_Forward(dsq, Ln, om, ox, NULL)           != eslOK) continue;
          if (p7_chkomx_GrowTo(cx, M, Ln, blockR)         != eslOK) esl_fatal(msg);
          if (p7_ForwardChk(dsq, Ln, om, cx, NULL)        != eslOK) esl_fatal(msg);
        }
      for (t = 0; t < ntr; t++) { seed[t] = 1 + esl_rnd_Roll(r, 1000000); esl_randomness_Init(rng[t], seed[t]); }
      if (use_sse) { if (p7_StochasticTraceChk_sse(rng, ntr, dsq, Ln, om, cx, tr2) != eslOK) esl_fatal(msg); }
      else         { if (p7_StochasticTraceChk    (rng, ntr, dsq, Ln, om, cx, tr2) != eslOK) esl_fatal(msg); }
      for (t = 0; t < ntr; t++)
        {
          esl_randomness_Init(rng1, seed[t]);
          if (use_sse) { if (p7_StochasticTrace_sse(rng1, dsq, Ln, om, ox, tr1) != eslOK) esl_fatal(msg); }
          else         { if (p7_StochasticTrace    (rng1, dsq, Ln, om, ox, tr1) != eslOK) esl_fatal(msg); }
          if (tr1->N != tr2[t]->N) esl_fatal("%s: trace %d lengths differ (%d %d)", msg, t, tr1->N, tr2[t]->N);
          for (z = 0; z < tr1->N; z++)
            if (tr1->st[z] != tr2[t]->st[z] || tr1->k[z] != tr2[t]->k[z] || tr1->i[z] != tr2[t]->i[z])
              esl_fatal("%s: trace %d differs at %d (block size %d)", msg, t, z, blockR);
          p7_trace_Reuse(tr1);
          p7_trace_Reuse(tr2[t]);
        }
      ntested++;
    }
  if (ntested == 0) esl_fatal("%s: no sequence was tested", msg);

  for (t = 0; t < ntr; t++) { esl_randomness_Destroy(rng[t]); p7_trace_Destroy(tr2[t]); }
  free(rng); free(tr2); free(seed); free(dsq);
  esl_randomness_Destroy(rng1);
  p7_trace_Destroy(tr1);
  esl_sq_Destroy(sq);
  if (use_sse) { p7_omx_Destroy_sse(ox); p7_chkomx_Destroy_sse(cx); p7_oprofile_Destroy_sse(om); }
  else         { p7_omx_Destroy(ox);     p7_chkomx_Destroy(cx);     p7_oprofile_Destroy(om); }
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}
#endif /*p7STOTRACE_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/



/*****************************************************************
 * 5. Test driver 
 *****************************************************************/
#ifdef p7STOTRACE_TESTDRIVE
/* gcc -std=gnu99 -msse2 -g -Wall -o stotrace_utest -Dp7STOTRACE_TESTDRIVE -I.. -L.. -I../../easel -L../../easel stotrace.c -lhmmer -leasel -lm
 */
#include "easel.h"
#include "esl_getopts.h"
#include "esl_randomseq.h"

#include "p7_config.h"
#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "--vv",      eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be very verbose",                                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for stochastic Viterbi traceback (optimized version)";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc    = NULL;
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = NULL;
  P7_OPROFILE    *om     = NULL;
  P7_BG          *bg     = NULL;
  ESL_DSQ        *dsq    = NULL;
  ESL_SQ         *sq     = NULL;
  int             M      = 6;
  int             L      = 10;
  int             ntrace = 1000;

  impl_Init();

  if ((abc = esl_alphabet_Create(eslAMINO))         == NULL)  esl_fatal("failed to create alphabet");
  if (p7_hmm_Sample(r, M, abc, &hmm)                != eslOK) esl_fatal("failed to sample an HMM");
  if ((bg = p7_bg_Create(abc))                      == NULL)  esl_fatal("failed to create null model");
  if ((gm = p7_profile_Create(hmm->M, abc))         == NULL)  esl_fatal("failed to create profile");
  if (p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL)    != eslOK) esl_fatal("failed to config profile");
  if ((om = p7_oprofile_Create(gm->M, abc))         == NULL)  esl_fatal("failed to create optimized profile");
  if (p7_oprofile_Convert(gm, om)                   != eslOK) esl_fatal("failed to convert profile");

  /* Test with randomly generated (iid) sequence */
  if ((dsq = malloc(sizeof(ESL_DSQ) *(L+2)))  == NULL)  esl_fatal("malloc failed");
  if (esl_rsq_xfIID(r, bg->f, abc->K, L, dsq) != eslOK) esl_fatal("seq generation failed");
  utest_stotrace(go, r, abc, gm, om, dsq, L, ntrace);

  /* Test with seq sampled from profile */
  if ((sq = esl_sq_CreateDigital(abc))             == NULL) esl_fatal("sequence allocation failed");
  if (p7_ProfileEmit(r, hmm, gm, bg, sq, NULL)    != eslOK) esl_fatal("profile emission failed");
  utest_stotrace(go, r, abc, gm, om, sq->dsq, sq->n, ntrace);

  /* Checkpointed Forward and traces sampled together */
  utest_stotrace_chk(r, abc, bg, 45, 200, 15, FALSE);
#ifdef eslENABLE_SSE
  utest_stotrace_chk(r, abc, bg, 45, 200, 15, TRUE);
#endif
   
  esl_sq_Destroy(sq);
  free(dsq);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7STOTRACE_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/


/*****************************************************************
 * 6. Example.
 *****************************************************************/
#ifdef p7STOTRACE_EXAMPLE
/* 
   gcc -g -Wall -msse2 -std=gnu99 -o stotrace_example -I.. -L.. -I../../easel -L../../easel -Dp7STOTRACE_EXAMPLE stotrace.c -lhmmer -leasel -lm
   ./stotrace_example <hmmfile> <seqfile>
 */ 

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "impl_avx.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",     0 },
  { "-m",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "dump the DP matrix to stdout",             0 },
  { "-p",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "dump the profile to stdout",               0 },
  { "-t",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "dump the sampled trace to stdout",         0 },
  { "-N",        eslARG_INT,      "1", NULL, NULL,  NULL,  NULL, NULL, "number of traces to sample",               0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of stochastic backtrace (SSE version)";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  ESL_RANDOMNESS *rng     = esl_randomness_CreateFast(0);
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_GMX         *gx      = NULL;
  P7_OMX         *fwd     = NULL;
  P7_TRACE       *tr      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  int             N       = esl_opt_GetInteger(go, "-N");
  int             i;
  float           vsc, fsc, tsc;
  char            errbuf[eslERRBUFSIZE];
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
  if  (esl_sqio_Read(sqfp, sq) != eslOK) p7_Fail("Failed to read sequence");

  /* create default null model, then create and optimize profile */
  bg = p7_bg_Create(abc);                p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);   p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);   p7_oprofile_Convert(gm, om);

  if (esl_opt_GetBoolean(go, "-p")) p7_oprofile_Dump(stdout, om);  

  fwd = p7_omx_Create(gm->M, sq->n, sq->n);
  gx  = p7_gmx_Create(gm->M, sq->n, sq->n, p7G_NSCELLS);
  tr  = p7_trace_Create();

  if (esl_opt_GetBoolean(go, "-m") == TRUE) p7_omx_SetDumpMode(stdout, fwd, TRUE); 
  p7_GViterbi(sq->dsq, sq->n, gm, gx,  &vsc);
  p7_Forward (sq->dsq, sq->n, om, fwd, &fsc);

  for (i = 0; i < N; i++)
    {
      p7_StochasticTrace(rng, sq->dsq, sq->n, om, fwd, tr);
      p7_trace_Score(tr, sq->dsq, gm, &tsc);
  
      if (esl_opt_GetBoolean(go, "-t") == TRUE) p7_trace_Dump(stdout, tr, gm, sq->dsq);
      if (p7_trace_Validate(tr, abc, sq->dsq, errbuf) != eslOK)  p7_Die("trace %d fails validation:\n%s\n", i, errbuf);

      printf("Sampled trace:  %.4f nats\n", tsc);
      p7_trace_Reuse(tr);
    }
  printf("Forward score:  %.4f nats\n", fsc);
  printf("Viterbi score:  %.4f nats\n", vsc);

  /* cleanup */
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_trace_Destroy(tr);
  p7_omx_Destroy(fwd);
  p7_gmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_randomness_Destroy(rng);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7STOTRACE_EXAMPLE*/
/*------------------------ end, example -------------------------*/


