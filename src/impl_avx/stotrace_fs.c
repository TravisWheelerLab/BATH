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


/* Forward declaration of dispatcher */
static int p7_StochasticTraceChk_Frameshift_Dispatcher(ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                                        P7_CHKOMX *cx, P7_OIVX *ov, P7_TRACE **tr);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_StochasticTraceChk_Frameshift)(ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                        P7_CHKOMX *cx, P7_OIVX *ov, P7_TRACE **tr) = p7_StochasticTraceChk_Frameshift_Dispatcher;

/* Function:  p7_StochasticTraceChk_Frameshift_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for frameshift stochastic traceback
 *            from a checkpointed Forward matrix.
 */
static int
p7_StochasticTraceChk_Frameshift_Dispatcher(ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                             P7_CHKOMX *cx, P7_OIVX *ov, P7_TRACE **tr)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_StochasticTraceChk_Frameshift = p7_StochasticTraceChk_Frameshift_avx; return p7_StochasticTraceChk_Frameshift_avx(rng, ntr, dsq, L, om_fs, cx, ov, tr); }
#endif
#ifdef eslENABLE_SSE
  p7_StochasticTraceChk_Frameshift = p7_StochasticTraceChk_Frameshift_sse;
  return p7_StochasticTraceChk_Frameshift_sse(rng, ntr, dsq, L, om_fs, cx, ov, tr);
#else
  p7_Die("p7_StochasticTraceChk_Frameshift: no SIMD implementation available");
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
#include "esl_vectorops.h"
#include "esl_randomseq.h"

/* utest_stotrace_fs()
 * Tests:
 *   1. Generic and SSE forward scores agree to within <tolerance> nats.
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

  tolerance = (p7_FLogsumError(-0.4, -0.5) > 0.0001) ? 0.1f : 0.001f;   /* nats; the generic Forward's logsum table rounds */

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
      if (fabsf(gfsc - ofsc) > tolerance)
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


/* utest_stotrace_posterior_fs()
 * SIMD stochastic traces must follow the posterior distribution. As the
 * generic utest_stotrace_posterior_frameshift(): for the back-translated
 * consensus of a sharpened model between random flanks with two nucleotides
 * deleted and two inserted, and frameshift probability 0.05, the mean number
 * per trace of M states of each codon length and of I states, over <N>
 * traces from p7_StochasticTrace_Frameshift(), must be within 5
 * standard errors (plus 0.01) of the posterior sums from the generic
 * p7_GDecoding_Frameshift(), and every trace must be valid. With <use_sse>,
 * the _sse versions are called directly.
 */
static void
utest_stotrace_posterior_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA,
                            P7_CODONTABLE *codon_table, int M, int N, int use_sse)
{
  char           *msg    = "frameshift stochastic trace posterior unit test failed";
  P7_HMM         *hmm    = NULL;
  P7_BG          *bgDNA  = p7_bg_Create(abcDNA);
  P7_FS_PROFILE  *gm_fs5 = p7_profile_fs_Create(M, abcAA, p7P_5CODONS);
  P7_FS_OPROFILE *om_fs5 = NULL;
  int             L      = 3*M + 120;
  P7_GMX         *fx     = p7_gmx_Create(M, L, L, p7G_NSCELLS_FS);
  P7_GMX         *bx     = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_IVX         *iv     = p7_ivx_Create(M, p7P_5CODONS);
  P7_OMX         *ox     = NULL;
  P7_OIVX        *ov     = NULL;
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq    = malloc(sizeof(ESL_DSQ) * (L+2));
  ESL_DSQ        *tmp    = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_TRACE       *tr     = p7_trace_fs_Create();
  double          expct[6], sum[6], sumsq[6], mean, se;
  char            errbuf[eslERRBUFSIZE];
  int             cnt[6];
  int             i, j, k, x, n, z, c, Ln, start, del1, del2, ins1, ins2;
  int (*fwd_full)(const ESL_DSQ *, int, const P7_FS_OPROFILE *, P7_OMX *, P7_OIVX *, float *)                = use_sse ? p7_Forward_Frameshift_sse : p7_Forward_Frameshift;
  int (*sto_full)(ESL_RANDOMNESS *, const ESL_DSQ *, int, const P7_FS_OPROFILE *, const P7_OMX *, P7_TRACE *) = use_sse ? p7_StochasticTrace_Frameshift_sse : p7_StochasticTrace_Frameshift;

  if (p7_hmm_Sample(r, M, abcAA, &hmm) != eslOK) esl_fatal(msg);
  for (k = 1; k <= M; k++)
    {
      x = esl_vec_FArgMax(hmm->mat[k], abcAA->K);
      esl_vec_FSet(hmm->mat[k], abcAA->K, 0.15 / (abcAA->K - 1));
      hmm->mat[k][x] = 0.85;
      if (k < M) { hmm->t[k][p7H_MM] = 0.9; hmm->t[k][p7H_MI] = 0.05; hmm->t[k][p7H_MD] = 0.05; }
    }
  hmm->fsprob = 0.05;
  p7_emit_SimpleConsensus(hmm, sq);

  Ln = 3*M + 100;
  esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, tmp);
  start = 1 + esl_rnd_Roll(r, Ln - 3*M + 1);
  for (i = 1, j = start; i <= M; i++, j += 3) p7_codontable_GetCodon(codon_table, r, sq->dsq[i], tmp+j);
  del1 = start + M/2;
  ins1 = start + M   + esl_rnd_Roll(r, M/4);
  del2 = start + 3*M/2;
  ins2 = start + 2*M + esl_rnd_Roll(r, M/4);
  for (i = 1, j = 1; i <= Ln; i++)
    {
      if (i == del1 || i == del2) continue;
      if (i == ins1 || i == ins2) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = tmp[i];
    }
  Ln = j - 1;

  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, Ln/3, p7_LOCAL) != eslOK) esl_fatal(msg);
  if (use_sse)
    {
      om_fs5 = p7_fs_oprofile_Create_sse(M, abcAA, 5);
      p7_fs_oprofile_Convert_sse(gm_fs5, om_fs5);
      ox     = p7_omx_Create_dpf_sse(M, L, L, p7X_NSCELLS_FS);
      ov     = p7_oivx_Create_sse(M, p7P_5CODONS);
      if (fwd_full(dsq, Ln, om_fs5, ox, ov, NULL) != eslOK) esl_fatal(msg);
    }
  else
    {
      om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
      p7_fs_oprofile_Convert(gm_fs5, om_fs5);
      ox     = p7_omx_Create_dpf(M, L, L, p7X_NSCELLS_FS);
      ov     = p7_oivx_Create(M, p7P_5CODONS);
      if (fwd_full(dsq, Ln, om_fs5, ox, ov, NULL) != eslOK) esl_fatal(msg);
    }

  for (n = 0; n < 6; n++) { sum[n] = sumsq[n] = 0.; }
  for (n = 0; n < N; n++)
    {
      if (sto_full(r, dsq, Ln, om_fs5, ox, tr)          != eslOK) esl_fatal(msg);
      if (p7_trace_fs_Validate(tr, abcDNA, dsq, errbuf)  != eslOK) esl_fatal("%s: invalid trace: %s", msg, errbuf);
      for (c = 0; c < 6; c++) cnt[c] = 0;
      for (z = 0; z < tr->N; z++)
        {
          if      (tr->st[z] == p7T_M) { if (tr->c[z] < 1 || tr->c[z] > 5) esl_fatal("%s: M with codon length %d", msg, tr->c[z]); cnt[tr->c[z]]++; }
          else if (tr->st[z] == p7T_I) cnt[0]++;
        }
      for (c = 0; c < 6; c++) { sum[c] += cnt[c]; sumsq[c] += (double) cnt[c] * cnt[c]; }
      p7_trace_Reuse(tr);
    }

  if (p7_GForward_Frameshift (dsq, Ln, gm_fs5, fx, iv, NULL) != eslOK) esl_fatal(msg);
  if (p7_GBackward_Frameshift(dsq, Ln, gm_fs5, bx, iv, NULL) != eslOK) esl_fatal(msg);
  if (p7_GDecoding_Frameshift(gm_fs5, fx, bx)               != eslOK) esl_fatal(msg);
  for (c = 0; c < 6; c++) expct[c] = 0.;
  for (i = 1; i <= Ln; i++)
    for (k = 1; k <= M; k++)
      {
        for (c = 1; c <= 5; c++) expct[c] += fx->dp[i][k*p7G_NSCELLS_FS + p7G_M + c];
        expct[0] += fx->dp[i][k*p7G_NSCELLS_FS + p7G_I];
      }
  for (c = 0; c < 6; c++)
    {
      mean = sum[c] / N;
      se   = sqrt(ESL_MAX(0., sumsq[c] / N - mean * mean) / N);
      if (fabs(mean - expct[c]) > 5. * se + 0.01)
        esl_fatal("%s: %s %d: mean count %.4f per trace, posterior sum %.4f (s.e. %.4f)", msg, c ? "M codon length" : "I", c, mean, expct[c], se);
    }

  free(dsq); free(tmp);
  p7_trace_fs_Destroy(tr);
  esl_sq_Destroy(sq);
  p7_ivx_Destroy(iv);
  p7_gmx_Destroy(fx);
  p7_gmx_Destroy(bx);
  p7_profile_fs_Destroy(gm_fs5);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bgDNA);
  if (use_sse)
    {
      p7_omx_Destroy_sse(ox);
      p7_oivx_Destroy_sse(ov);
      p7_fs_oprofile_Destroy_sse(om_fs5);
    }
  else
    {
      p7_omx_Destroy(ox);
      p7_oivx_Destroy(ov);
      p7_fs_oprofile_Destroy(om_fs5);
    }
}

/* utest_stotrace_chk_fs()
 * Stochastic traces from a checkpointed frameshift Forward matrix against
 * p7_StochasticTrace_Frameshift() on the full Forward matrix. For
 * random DNA and the back-translated consensus of a model with sharpened
 * match emissions between random flanks, <ntr> traces sampled together
 * by p7_StochasticTraceChk_Frameshift() must be identical (states,
 * k, i and codon lengths) to those sampled one at a time with the same
 * per-trace random number generators, with the checkpointed matrix in one
 * block and in blocks of 5..12 rows. With <use_sse>, the _sse versions are called directly, so an AVX CPU
 * also tests the SSE engine.
 */
static void
utest_stotrace_chk_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA,
                      P7_CODONTABLE *codon_table, int M, int L, int N, int use_sse)
{
  char            *msg    = "checkpointed frameshift stochastic trace unit test failed";
  int              ntr    = 20;
  P7_HMM          *hmm    = NULL;
  P7_BG           *bgDNA  = p7_bg_Create(abcDNA);
  P7_FS_PROFILE   *gm_fs5 = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE  *om_fs5 = NULL;
  P7_OMX          *fwd    = NULL;
  P7_OIVX         *ov5    = NULL;
  P7_CHKOMX       *cx     = NULL;
  ESL_SQ          *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ         *dsq    = malloc(sizeof(ESL_DSQ) * (L+2));
  ESL_RANDOMNESS  *rng1   = esl_randomness_Create(1);
  ESL_RANDOMNESS **rng    = malloc(sizeof(ESL_RANDOMNESS *) * ntr);
  P7_TRACE        *tr1    = p7_trace_fs_Create();
  P7_TRACE       **tr2    = malloc(sizeof(P7_TRACE *) * ntr);
  uint32_t        *seed   = malloc(sizeof(uint32_t) * ntr);
  int              ntested = 0;
  int              n, Ln, blockR, start, i, j, t, z, k, x;
  int            (*reconfig)  (P7_FS_OPROFILE *, int)                                                         = use_sse ? p7_fs_oprofile_ReconfigLength_sse : p7_fs_oprofile_ReconfigLength;
  int            (*fwd_full)  (const ESL_DSQ *, int, const P7_FS_OPROFILE *, P7_OMX *, P7_OIVX *, float *)     = use_sse ? p7_Forward_Frameshift_sse : p7_Forward_Frameshift;
  int            (*fwd_chk)   (const ESL_DSQ *, int, const P7_FS_OPROFILE *, P7_CHKOMX *, P7_OIVX *, float *)  = use_sse ? p7_ForwardChk_Frameshift_sse : p7_ForwardChk_Frameshift;
  int            (*chk_growto)(P7_CHKOMX *, int, int, int)                                                     = use_sse ? p7_chkomx_GrowTo_sse : p7_chkomx_GrowTo;
  int            (*sto_full)  (ESL_RANDOMNESS *, const ESL_DSQ *, int, const P7_FS_OPROFILE *, const P7_OMX *, P7_TRACE *)               = use_sse ? p7_StochasticTrace_Frameshift_sse : p7_StochasticTrace_Frameshift;
  int            (*sto_chk)   (ESL_RANDOMNESS **, int, const ESL_DSQ *, int, const P7_FS_OPROFILE *, P7_CHKOMX *, P7_OIVX *, P7_TRACE **) = use_sse ? p7_StochasticTraceChk_Frameshift_sse : p7_StochasticTraceChk_Frameshift;

  for (t = 0; t < ntr; t++) { rng[t] = esl_randomness_Create(1); tr2[t] = p7_trace_fs_Create(); }

  if (p7_hmm_Sample(r, M, abcAA, &hmm) != eslOK) esl_fatal(msg);
  for (k = 1; k <= M; k++)
    {
      x = esl_vec_FArgMax(hmm->mat[k], abcAA->K);
      esl_vec_FSet(hmm->mat[k], abcAA->K, 0.15 / (abcAA->K - 1));
      hmm->mat[k][x] = 0.85;
      if (k < M) { hmm->t[k][p7H_MM] = 0.9; hmm->t[k][p7H_MI] = 0.05; hmm->t[k][p7H_MD] = 0.05; }
    }
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_LOCAL) != eslOK) esl_fatal(msg);
  p7_emit_SimpleConsensus(hmm, sq);
  if (use_sse)
    {
      om_fs5 = p7_fs_oprofile_Create_sse(M, abcAA, 5);
      p7_fs_oprofile_Convert_sse(gm_fs5, om_fs5);
      fwd    = p7_omx_Create_dpf_sse(M, L, L, p7X_NSCELLS_FS);
      ov5    = p7_oivx_Create_sse(M, p7P_5CODONS);
      cx     = p7_chkomx_Create_sse(M, L, L+1, 5, p7X_NSCELLS_FS, p7_CHK_FORWARD);
    }
  else
    {
      om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
      p7_fs_oprofile_Convert(gm_fs5, om_fs5);
      fwd    = p7_omx_Create_dpf(M, L, L, p7X_NSCELLS_FS);
      ov5    = p7_oivx_Create(M, p7P_5CODONS);
      cx     = p7_chkomx_Create(M, L, L+1, 5, p7X_NSCELLS_FS, p7_CHK_FORWARD);
    }

  for (n = 0; n < N; n++)
    {
      Ln = 15 + esl_rnd_Roll(r, L - 14);
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, dsq);
      if (n % 2 == 1 && 3*M <= L)
        {
          Ln    = ESL_MAX(Ln, 3*M);
          esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, dsq);
          start = 1 + esl_rnd_Roll(r, Ln - 3*M + 1);
          for (i = 1, j = start; i <= M; i++, j += 3)
            p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq+j);
        }
      blockR = (n % 3 == 2) ? Ln+1 : 5 + esl_rnd_Roll(r, 8);

      reconfig(om_fs5, Ln/3);
      if (fwd_full(dsq, Ln, om_fs5, fwd, ov5, NULL) != eslOK) continue;
      if (chk_growto(cx, M, Ln, blockR)                != eslOK) esl_fatal(msg);
      if (fwd_chk(dsq, Ln, om_fs5, cx, ov5, NULL)      != eslOK) esl_fatal(msg);

      for (t = 0; t < ntr; t++) { seed[t] = 1 + esl_rnd_Roll(r, 1000000); esl_randomness_Init(rng[t], seed[t]); }
      if (sto_chk(rng, ntr, dsq, Ln, om_fs5, cx, ov5, tr2) != eslOK) esl_fatal(msg);
      for (t = 0; t < ntr; t++)
        {
          esl_randomness_Init(rng1, seed[t]);
          if (sto_full(rng1, dsq, Ln, om_fs5, fwd, tr1) != eslOK) esl_fatal(msg);
          if (tr1->N != tr2[t]->N) esl_fatal("%s: trace %d lengths differ (%d %d)", msg, t, tr1->N, tr2[t]->N);
          for (z = 0; z < tr1->N; z++)
            if (tr1->st[z] != tr2[t]->st[z] || tr1->k[z] != tr2[t]->k[z] || tr1->i[z] != tr2[t]->i[z] || tr1->c[z] != tr2[t]->c[z])
              esl_fatal("%s: trace %d differs at %d (block size %d)", msg, t, z, blockR);
          p7_trace_Reuse(tr1);
          p7_trace_Reuse(tr2[t]);
        }
      ntested++;
    }
  if (ntested == 0) esl_fatal("%s: no sequence was tested", msg);

  for (t = 0; t < ntr; t++) { esl_randomness_Destroy(rng[t]); p7_trace_fs_Destroy(tr2[t]); }
  free(rng); free(tr2); free(seed); free(dsq);
  esl_randomness_Destroy(rng1);
  p7_trace_fs_Destroy(tr1);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bgDNA);
  p7_profile_fs_Destroy(gm_fs5);
  if (use_sse)
    {
      p7_omx_Destroy_sse(fwd);
      p7_oivx_Destroy_sse(ov5);
      p7_chkomx_Destroy_sse(cx);
      p7_fs_oprofile_Destroy_sse(om_fs5);
    }
  else
    {
      p7_omx_Destroy(fwd);
      p7_oivx_Destroy(ov5);
      p7_chkomx_Destroy(cx);
      p7_fs_oprofile_Destroy(om_fs5);
    }
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
  utest_stotrace_chk_fs(r, abcAA, abcDNA, gcode, bgAA, ct, M, 3*M+150, N, FALSE);
  utest_stotrace_posterior_fs(r, abcAA, abcDNA, gcode, bgAA, ct, 45, 2000, FALSE);
#ifdef eslENABLE_SSE
  utest_stotrace_chk_fs(r, abcAA, abcDNA, gcode, bgAA, ct, M, 3*M+150, N, TRUE);
  utest_stotrace_posterior_fs(r, abcAA, abcDNA, gcode, bgAA, ct, 45, 2000, TRUE);
#endif

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


