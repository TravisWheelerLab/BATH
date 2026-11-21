/* Spliced Viterbi algorithms dispatcher.
 * Provides the non-suffixed p7_Viterbi_Spliced(), p7_Viterbi_SplicedTraceOld(),
 * p7_Viterbi_SplicedTrace(), p7_Viterbi_SplicedChk() and
 * p7_Viterbi_SplicedTraceChk() declared in impl_avx.h.  Delegates to the fastest available ISA implementation.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Forward declaration of dispatcher */
static int p7_Viterbi_Spliced_Dispatcher(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_OMX *ox,
                                         const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                                         int i_start, int i_end, int min_intron, int global_start, int global_end);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Viterbi_Spliced)(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_OMX *ox,
                          const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                          int i_start, int i_end, int min_intron, int global_start, int global_end) = p7_Viterbi_Spliced_Dispatcher;

/* Function:  p7_Viterbi_Spliced_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for spliced Viterbi DP fill.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on bad inputs.
 */
static int
p7_Viterbi_Spliced_Dispatcher(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_OMX *ox,
                               const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                               int i_start, int i_end, int min_intron, int global_start, int global_end)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_Viterbi_Spliced = p7_Viterbi_Spliced_avx512; return p7_Viterbi_Spliced_avx512(sub_dsq, om_tr, ox, signal_scores, acc_ov, don_ov, i_start, i_end, min_intron, global_start, global_end); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Viterbi_Spliced = p7_Viterbi_Spliced_avx;    return p7_Viterbi_Spliced_avx(sub_dsq, om_tr, ox, signal_scores, acc_ov, don_ov, i_start, i_end, min_intron, global_start, global_end); }
#endif
#ifdef eslENABLE_SSE
  p7_Viterbi_Spliced = p7_Viterbi_Spliced_sse;
  return p7_Viterbi_Spliced_sse(sub_dsq, om_tr, ox, signal_scores, acc_ov, don_ov, i_start, i_end, min_intron, global_start, global_end);
#else
  p7_Die("p7_Viterbi_Spliced: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_Viterbi_SplicedTraceOld_Dispatcher(const ESL_DSQ *sub_dsq, const P7_OMX *ox,
                                               const P7_FS_PROFILE *gm_tr, const float *signal_scores,
                                               P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                                               int min_intron, float *vitsc, float *max_drop);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Viterbi_SplicedTraceOld)(const ESL_DSQ *sub_dsq, const P7_OMX *ox,
                               const P7_FS_PROFILE *gm_tr, const float *signal_scores,
                               P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                               int min_intron, float *vitsc, float *max_drop) = p7_Viterbi_SplicedTraceOld_Dispatcher;

/* Function:  p7_Viterbi_SplicedTraceOld_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for the spliced Viterbi traceback
 *            kept for reference (p7_GViterbi_SplicedTraceOld()).
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile is not 1-codon-length.
 *            <eslFAIL>   if traceback fails.
 */
static int
p7_Viterbi_SplicedTraceOld_Dispatcher(const ESL_DSQ *sub_dsq, const P7_OMX *ox,
                                    const P7_FS_PROFILE *gm_tr, const float *signal_scores,
                                    P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                                    int min_intron, float *vitsc, float *max_drop)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_Viterbi_SplicedTraceOld = p7_Viterbi_SplicedTrace_avx512; return p7_Viterbi_SplicedTrace_avx512(sub_dsq, ox, gm_tr, signal_scores, tr, i_start, i_end, k_start, k_end, min_intron, vitsc, max_drop); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Viterbi_SplicedTraceOld = p7_Viterbi_SplicedTraceOld_avx;    return p7_Viterbi_SplicedTraceOld_avx(sub_dsq, ox, gm_tr, signal_scores, tr, i_start, i_end, k_start, k_end, min_intron, vitsc, max_drop); }
#endif
#ifdef eslENABLE_SSE
  p7_Viterbi_SplicedTraceOld = p7_Viterbi_SplicedTraceOld_sse;
  return p7_Viterbi_SplicedTraceOld_sse(sub_dsq, ox, gm_tr, signal_scores, tr, i_start, i_end, k_start, k_end, min_intron, vitsc, max_drop);
#else
  p7_Die("p7_Viterbi_SplicedTraceOld: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_Viterbi_SplicedTrace_Dispatcher(const ESL_DSQ *sub_dsq, const P7_OMX *ox,
                                                     const P7_FS_PROFILE *gm_tr, const float *signal_scores,
                                                     P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                                                     int min_intron, float *vitsc, float *max_drop);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Viterbi_SplicedTrace)(const ESL_DSQ *sub_dsq, const P7_OMX *ox,
                                     const P7_FS_PROFILE *gm_tr, const float *signal_scores,
                                     P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                                     int min_intron, float *vitsc, float *max_drop) = p7_Viterbi_SplicedTrace_Dispatcher;

/* Function:  p7_Viterbi_SplicedTrace_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for spliced Viterbi traceback.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile is not 1-codon-length.
 *            <eslFAIL>   if traceback fails.
 */
static int
p7_Viterbi_SplicedTrace_Dispatcher(const ESL_DSQ *sub_dsq, const P7_OMX *ox,
                                          const P7_FS_PROFILE *gm_tr, const float *signal_scores,
                                          P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                                          int min_intron, float *vitsc, float *max_drop)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Viterbi_SplicedTrace = p7_Viterbi_SplicedTrace_avx; return p7_Viterbi_SplicedTrace_avx(sub_dsq, ox, gm_tr, signal_scores, tr, i_start, i_end, k_start, k_end, min_intron, vitsc, max_drop); }
#endif
#ifdef eslENABLE_SSE
  p7_Viterbi_SplicedTrace = p7_Viterbi_SplicedTrace_sse;
  return p7_Viterbi_SplicedTrace_sse(sub_dsq, ox, gm_tr, signal_scores, tr, i_start, i_end, k_start, k_end, min_intron, vitsc, max_drop);
#else
  p7_Die("p7_Viterbi_SplicedTrace: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_Viterbi_SplicedChk_Dispatcher(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                                            const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                                            int i_start, int i_end, int min_intron, int global_start, int global_end);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Viterbi_SplicedChk)(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                             const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                             int i_start, int i_end, int min_intron, int global_start, int global_end) = p7_Viterbi_SplicedChk_Dispatcher;

/* Function:  p7_Viterbi_SplicedChk_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for the checkpointed spliced
 *            Viterbi fill.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on bad inputs.
 */
static int
p7_Viterbi_SplicedChk_Dispatcher(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                                 const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                                 int i_start, int i_end, int min_intron, int global_start, int global_end)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Viterbi_SplicedChk = p7_Viterbi_SplicedChk_avx; return p7_Viterbi_SplicedChk_avx(sub_dsq, om_tr, cx, cv, signal_scores, acc_ov, don_ov, i_start, i_end, min_intron, global_start, global_end); }
#endif
#ifdef eslENABLE_SSE
  p7_Viterbi_SplicedChk = p7_Viterbi_SplicedChk_sse;
  return p7_Viterbi_SplicedChk_sse(sub_dsq, om_tr, cx, cv, signal_scores, acc_ov, don_ov, i_start, i_end, min_intron, global_start, global_end);
#else
  p7_Die("p7_Viterbi_SplicedChk: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_Viterbi_SplicedTraceChk_Dispatcher(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                                                 const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                                                 P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                                                 int min_intron, int global_start, int global_end, float *vitsc, float *max_drop);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Viterbi_SplicedTraceChk)(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                                  const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                                  P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                                  int min_intron, int global_start, int global_end, float *vitsc, float *max_drop) = p7_Viterbi_SplicedTraceChk_Dispatcher;

/* Function:  p7_Viterbi_SplicedTraceChk_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for the checkpointed spliced
 *            Viterbi traceback.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile is not 1-codon-length or the matrix
 *                        isn't filled for this comparison.
 *            <eslFAIL>   if traceback fails.
 */
static int
p7_Viterbi_SplicedTraceChk_Dispatcher(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                                      const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                                      P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                                      int min_intron, int global_start, int global_end, float *vitsc, float *max_drop)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Viterbi_SplicedTraceChk = p7_Viterbi_SplicedTraceChk_avx; return p7_Viterbi_SplicedTraceChk_avx(sub_dsq, om_tr, cx, cv, gm_tr, signal_scores, acc_ov, don_ov, tr, i_start, i_end, k_start, k_end, min_intron, global_start, global_end, vitsc, max_drop); }
#endif
#ifdef eslENABLE_SSE
  p7_Viterbi_SplicedTraceChk = p7_Viterbi_SplicedTraceChk_sse;
  return p7_Viterbi_SplicedTraceChk_sse(sub_dsq, om_tr, cx, cv, gm_tr, signal_scores, acc_ov, don_ov, tr, i_start, i_end, k_start, k_end, min_intron, global_start, global_end, vitsc, max_drop);
#else
  p7_Die("p7_Viterbi_SplicedTraceChk: no SIMD implementation available");
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/
#ifdef p7VITERBI_SP_BENCHMARK
/*
  gcc -o benchmark-viterbi_sp -std=gnu99 -g -Wall -mavx2 -I.. -L.. -I../../easel -L../../easel -Dp7VITERBI_SP_BENCHMARK viterbi_sp.c -lhmmer -leasel -lm

   ./benchmark-viterbi_sp <hmmfile>
 */

#include "p7_config.h"

#include <math.h>
#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_avx.h"
#include "p7_splice.h"

static ESL_OPTIONS benchmark_options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-N",        eslARG_INT,    "100", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "-I",        eslARG_INT,    "200", NULL, "n>0", NULL,  NULL, NULL, "length of simulated intron (excl. GT..AG signals)", 0 },
  { "-T",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also benchmark Trace after each DP",   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char benchmark_usage[]  = "[-options] <hmmfile>";
static char benchmark_banner[] = "benchmark driver for spliced Viterbi algorithms";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go           = p7_CreateDefaultApp(benchmark_options, 1, argc, argv, benchmark_banner, benchmark_usage);
  char           *hmmfile      = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w            = esl_stopwatch_Create();
  ESL_RANDOMNESS *r            = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA        = NULL;
  ESL_ALPHABET   *abcDNA       = esl_alphabet_Create(eslDNA);
  P7_HMMFILE     *hfp          = NULL;
  P7_HMM         *hmm          = NULL;
  P7_BG          *bgAA         = NULL;
  P7_PROFILE     *gm           = NULL;
  P7_FS_PROFILE  *gm_tr        = NULL;
  P7_FS_OPROFILE *om_tr        = NULL;
  P7_TRACE       *tr           = NULL;
  ESL_GENCODE    *gcode        = NULL;
  P7_CODONTABLE  *codon_table  = NULL;
  ESL_SQ         *sq           = NULL;
  SPLICE_PIPELINE *pli         = NULL;
  P7_OMX          *ox          = NULL;
  int             do_T         = esl_opt_GetBoolean(go, "-T");
  int             N            = esl_opt_GetInteger(go, "-N");
  int             I            = esl_opt_GetInteger(go, "-I");
  int             intron_total = I + 4;
  ESL_DSQ        *dsq          = NULL;
  int             i, j, k, L_amino, L_dna_total;
  int64_t         total_cells;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)          != eslOK) p7_Fail("Failed to read HMM");

  gcode       = esl_gencode_Create(abcDNA, abcAA);
  bgAA        = p7_bg_Create(abcAA);
  gm          = p7_profile_Create(hmm->M, abcAA);
  gm_tr       = p7_profile_fs_Create(hmm->M, abcAA, 1);
  om_tr       = p7_fs_oprofile_Create(hmm->M, abcAA, 1);
  codon_table = p7_codontable_Create(gcode);
  sq          = esl_sq_CreateDigital(abcAA);

  p7_ProfileConfig   (hmm, bgAA,        gm,    hmm->M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, hmm->M, p7_UNILOCAL);
  p7_fs_oprofile_Convert_Log(gm_tr, om_tr);

  pli = p7_splicepipeline_Create(NULL, hmm->M, hmm->M * 3);
  ox  = p7_omx_Create_dpf(hmm->M, hmm->M * 3, hmm->M * 3, p7X_NSCELLS);

  tr = p7_trace_fs_Create();
  /* Baseline: time to generate sequences alone */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, NULL);
      L_amino     = sq->n;
      L_dna_total = L_amino * 3 + intron_total;
      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) p7_Fail("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;
      j = 1;
      for (k = 1; k <= L_amino / 2; k++) { p7_codontable_GetCodon(codon_table, r, sq->dsq[k], dsq + j); j += 3; }
      dsq[j++] = 2;  /* G */
      dsq[j++] = 3;  /* T */
      for (k = 0; k < I; k++) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = 0;  /* A */
      dsq[j++] = 2;  /* G */
      for (k = L_amino / 2 + 1; k <= L_amino; k++) { p7_codontable_GetCodon(codon_table, r, sq->dsq[k], dsq + j); j += 3; }
    }
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark */
  total_cells = 0;
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, NULL);
      L_amino     = sq->n;
      L_dna_total = L_amino * 3 + intron_total;
      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) p7_Fail("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;
      j = 1;
      for (k = 1; k <= L_amino / 2; k++) { p7_codontable_GetCodon(codon_table, r, sq->dsq[k], dsq + j); j += 3; }
      dsq[j++] = 2;  /* G */
      dsq[j++] = 3;  /* T */
      for (k = 0; k < I; k++) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = 0;  /* A */
      dsq[j++] = 2;  /* G */
      for (k = L_amino / 2 + 1; k <= L_amino; k++) { p7_codontable_GetCodon(codon_table, r, sq->dsq[k], dsq + j); j += 3; }

      p7_fs_ReconfigLength(gm_tr, L_dna_total/3);
      p7_fs_oprofile_ReconfigLength_Log(om_tr, L_dna_total/3);
      p7_omx_GrowTo_dpf(ox, hmm->M, L_dna_total, L_dna_total);
      p7_Viterbi_Spliced(dsq, om_tr, ox, pli->signal_scores, pli->acc_ov, pli->don_ov, 1, L_dna_total, pli->min_intron, TRUE, TRUE);

      if(do_T)
        p7_Viterbi_SplicedTrace(dsq, ox, gm_tr, pli->signal_scores, tr, 1, L_dna_total, 1, hmm->M, pli->min_intron, NULL, NULL);
      p7_trace_Reuse(tr);
      total_cells += (int64_t) L_dna_total * hmm->M;
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) total_cells * 1e-6 / bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M          = %d\n",   hmm->M);
  printf("# N          = %d\n",   N);
  printf("# I          = %d\n",   I);
  printf("# %.1f Mc/s\n", Mcs);

  if (dsq != NULL) free(dsq);
  esl_sq_Destroy(sq);
  p7_codontable_Destroy(codon_table);
  p7_profile_fs_Destroy(gm_tr);
  p7_fs_oprofile_Destroy(om_tr);
  p7_omx_Destroy(ox);
  p7_splicepipeline_Destroy(pli);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bgAA);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  p7_trace_fs_Destroy(tr);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}

#endif /*p7VITERBI_SP_BENCHMARK*/
/*----------------- end, benchmark driver -----------------------*/


/*****************************************************************
 * 4. Unit tests.
 *****************************************************************/
#ifdef p7VITERBI_SP_TESTDRIVE

#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "p7_splice.h"

/* utest_compare_traces()
 *
 * Require two traces of the same matrix to have the same states,
 * positions, splice codon phases and Viterbi scores.
 */
static void
utest_compare_traces(const P7_TRACE *tr1, float sc1, const P7_TRACE *tr2, float sc2, const char *msg)
{
  int z;

  if (tr1->N != tr2->N) esl_fatal("%s: trace lengths %d != %d", msg, tr1->N, tr2->N);
  for (z = 0; z < tr1->N; z++)
    if (tr1->st[z] != tr2->st[z] || tr1->k[z] != tr2->k[z] ||
        tr1->i[z]  != tr2->i[z]  || tr1->c[z] != tr2->c[z])
      esl_fatal("%s: traces differ at position %d", msg, z);
  if (sc1 != sc2) esl_fatal("%s: trace scores %f != %f", msg, sc1, sc2);
}

static void
utest_viterbi_sp(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
              ESL_GENCODE *gcode, P7_BG *bgAA, P7_CODONTABLE *codon_table,
              int M, int N, int intron_len)
{
  char           *msg         = "spliced viterbi unit test failed";
  P7_HMM         *hmm         = NULL;
  P7_FS_PROFILE  *gm_tr       = p7_profile_fs_Create(M, abcAA, 1);
  P7_FS_OPROFILE *om_tr       = p7_fs_oprofile_Create(M, abcAA, 1);
  ESL_SQ         *sq          = esl_sq_CreateDigital(abcAA);
  P7_TRACE       *otr         = p7_trace_fs_Create();
  P7_TRACE       *gtr         = p7_trace_fs_Create();
  P7_TRACE       *rtr         = p7_trace_fs_Create();
  ESL_DSQ        *dsq         = NULL;
  P7_GMX         *gx          = NULL;
  P7_IVX         *acc_iv      = NULL;
  P7_IVX         *don_iv      = NULL;
  SPLICE_PIPELINE *pli        = NULL;
  P7_OMX          *ox         = NULL;
  int             intron_total = intron_len + 4;  /* GT + intron_len random nucs + AG */
  int             L_amino, L_dna_total;
  int             i, j;
  int             k_start, k_end;
  int             sub_M;
  float           final_gC, final_oC;
  float           sc_search, sc_replay;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, M, p7_UNILOCAL);
  p7_fs_oprofile_Convert_Log(gm_tr, om_tr);

  pli = p7_splicepipeline_Create(NULL, M, M * 3);
  ox  = p7_omx_Create_dpf(M, M * 3, M * 3, p7X_NSCELLS);
  gx     = p7_gmx_Create(M, M, M, p7G_NSCELLS);
  acc_iv = p7_ivx_Create(M, SPLICE_ROWS);
  don_iv = p7_ivx_Create(M, SIGNAL_MEM_SIZE);
  float signal_scores[p7S_SPLICE_SIGNALS];
  p7_SignalScores(signal_scores);

  p7_emit_SimpleConsensus(hmm, sq);
  L_amino     = sq->n;

  while (N--)
    {
	  k_start = k_end = 0;
	  while(k_end <= k_start || k_end > L_amino) {
	    k_start = esl_rnd_Roll(r, L_amino/4) + 1;
	    k_end   = esl_rnd_Roll(r, L_amino/2) + k_start + L_amino/4;
	  }
      sub_M = k_end - k_start + 1;
      
      L_dna_total = sub_M * 3 + intron_total;
      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) esl_fatal("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;

      /* Reverse-translate first half of the sequence (exon 1) */
      j = 1;
      for (i = k_start; i <= k_start+sub_M/2; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }
      /* Simulated intron: GT + intron_len random nucleotides + AG
       * eslDNA alphabet: A=0, C=1, G=2, T=3                        */
      dsq[j++] = 2;  /* G */
      dsq[j++] = 3;  /* T */
      for (i = 0; i < intron_len; i++) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = 0;  /* A */
      dsq[j++] = 2;  /* G */
      /* Reverse-translate second half of the sequence (exon 2) */
      for (i = k_start+ sub_M/2 + 1; i <= k_end; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }

      /* --- Score comparison: generic Dummy and SSE Dummy must agree --- */
	  p7_fs_oprofile_SubConvert_Log(gm_tr, om_tr, k_start, k_end);
      p7_fs_oprofile_ReconfigLength_Log(om_tr, L_dna_total/3);
	  p7_omx_GrowTo_dpf(ox, sub_M, L_dna_total, L_dna_total);
	  p7_oivx_GrowTo(pli->acc_ov, sub_M, SPLICE_ROWS);
	  p7_oivx_GrowTo(pli->don_ov, sub_M, SIGNAL_MEM_SIZE);

      p7_Viterbi_Spliced(dsq, om_tr, ox, signal_scores, pli->acc_ov, pli->don_ov, 1, L_dna_total, pli->min_intron, TRUE, TRUE);
      final_oC = ox->xmx[L_dna_total * p7X_NXCELLS + p7X_C];

      p7_fs_ReconfigLength(gm_tr, L_dna_total/3);
      p7_gmx_GrowTo(gx, sub_M, L_dna_total, L_dna_total);
      p7_ivx_GrowTo(acc_iv, sub_M, SPLICE_ROWS);
	  p7_ivx_GrowTo(don_iv, sub_M, SIGNAL_MEM_SIZE);

      p7_GViterbi_Spliced(dsq, gm_tr, gx, acc_iv, don_iv, signal_scores, 1, L_dna_total, k_start, k_end, pli->min_intron, TRUE, TRUE);
      final_gC = gx->xmx[L_dna_total * p7G_NXCELLS + p7G_C];

      if (fabs(final_gC - final_oC) > 0.001) 
        esl_fatal("%s: generic %.4f != SSE %.4f", msg, final_gC, final_oC); 

      if(final_gC == -eslINFINITY) continue;   
      
      p7_Viterbi_SplicedTraceOld      (dsq, ox, gm_tr, signal_scores, otr, 1, L_dna_total, k_start, k_end, pli->min_intron, &sc_search, NULL);
      p7_Viterbi_SplicedTrace(dsq, ox, gm_tr, signal_scores, rtr, 1, L_dna_total, k_start, k_end, pli->min_intron, &sc_replay, NULL);
      p7_GViterbi_SplicedTraceOld(dsq, gm_tr, gx, signal_scores, gtr, 1, L_dna_total, k_start, k_end, pli->min_intron, NULL, NULL);

      p7_trace_Compare(otr, gtr, 0.0);
      utest_compare_traces(otr, sc_search, rtr, sc_replay, msg);
      p7_trace_Reuse(rtr);
        
      p7_trace_Reuse(otr);
      p7_trace_Reuse(gtr);

      p7_Viterbi_Spliced(dsq, om_tr, ox, signal_scores, pli->acc_ov, pli->don_ov, 1, L_dna_total, pli->min_intron, TRUE, FALSE);
      final_oC = ox->xmx[L_dna_total * p7X_NXCELLS + p7X_C];

      p7_GViterbi_Spliced(dsq, gm_tr, gx, acc_iv, don_iv, signal_scores, 1, L_dna_total, k_start, k_end, pli->min_intron, TRUE, FALSE);
      final_gC = gx->xmx[L_dna_total * p7G_NXCELLS + p7G_C];
       
      p7_Viterbi_SplicedTraceOld      (dsq, ox, gm_tr, signal_scores, otr, 1, L_dna_total, k_start, k_end, pli->min_intron, &sc_search, NULL);
      p7_Viterbi_SplicedTrace(dsq, ox, gm_tr, signal_scores, rtr, 1, L_dna_total, k_start, k_end, pli->min_intron, &sc_replay, NULL);
      p7_GViterbi_SplicedTraceOld(dsq, gm_tr, gx, signal_scores, gtr, 1, L_dna_total, k_start, k_end, pli->min_intron, NULL, NULL);

      p7_trace_Compare(otr, gtr, 0.0);
      utest_compare_traces(otr, sc_search, rtr, sc_replay, msg);
      p7_trace_Reuse(rtr);

      p7_trace_Reuse(otr);
      p7_trace_Reuse(gtr);

      p7_Viterbi_Spliced(dsq, om_tr, ox, signal_scores, pli->acc_ov, pli->don_ov, 1, L_dna_total, pli->min_intron, FALSE, TRUE);
      final_oC = ox->xmx[L_dna_total * p7X_NXCELLS + p7X_C];

      p7_GViterbi_Spliced(dsq, gm_tr, gx, acc_iv, don_iv, signal_scores, 1, L_dna_total, k_start, k_end, pli->min_intron, FALSE, TRUE);
      final_gC = gx->xmx[L_dna_total * p7G_NXCELLS + p7G_C];
      
      p7_Viterbi_SplicedTraceOld      (dsq, ox, gm_tr, signal_scores, otr, 1, L_dna_total, k_start, k_end, pli->min_intron, &sc_search, NULL);
      p7_Viterbi_SplicedTrace(dsq, ox, gm_tr, signal_scores, rtr, 1, L_dna_total, k_start, k_end, pli->min_intron, &sc_replay, NULL);
      p7_GViterbi_SplicedTraceOld(dsq, gm_tr, gx, signal_scores, gtr, 1, L_dna_total, k_start, k_end, pli->min_intron, NULL, NULL);

      p7_trace_Compare(otr, gtr, 0.0);
      utest_compare_traces(otr, sc_search, rtr, sc_replay, msg);
      p7_trace_Reuse(rtr);

      p7_trace_Reuse(otr);
      p7_trace_Reuse(gtr);
    }

  if (dsq != NULL) free(dsq);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_profile_fs_Destroy(gm_tr);
  p7_fs_oprofile_Destroy(om_tr);
  p7_omx_Destroy(ox);
  p7_splicepipeline_Destroy(pli);
  p7_gmx_Destroy(gx);
  p7_ivx_Destroy(acc_iv);
  p7_ivx_Destroy(don_iv);
  p7_trace_fs_Destroy(otr);
  p7_trace_fs_Destroy(gtr);
  p7_trace_fs_Destroy(rtr);
}

#if defined(eslENABLE_SSE) && defined(eslENABLE_AVX)
/* utest_sse_vs_avx_viterbi_sp()
 *
 * Run Viterbi_Spliced_sse and Viterbi_Spliced_avx on the same sequences
 * and profile, comparing final C-state scores.  Creates ISA-specific
 * DP matrices directly (no SPLICE_PIPELINE wrapper) so each ISA operates
 * on its own allocated memory.  Tolerance is 0.01 nats.
 * Skipped silently if AVX is not available at runtime.
 */
static void
utest_sse_vs_avx_viterbi_sp(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
                             ESL_GENCODE *gcode, P7_BG *bgAA, P7_CODONTABLE *codon_table,
                             int M, int N, int intron_len)
{
  char            msg[]         = "utest_sse_vs_avx_viterbi_sp failed";
  P7_HMM         *hmm           = NULL;
  P7_FS_PROFILE  *gm_tr         = p7_profile_fs_Create(M, abcAA, 1);
  P7_FS_OPROFILE *om_tr_sse     = NULL;
  P7_FS_OPROFILE *om_tr_avx     = NULL;
  ESL_SQ         *sq            = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq           = NULL;
  P7_OMX         *ox_sse        = NULL;
  P7_OMX         *ox_avx        = NULL;
  P7_OIVX        *acc_ov_sse    = NULL;
  P7_OIVX        *acc_ov_avx    = NULL;
  P7_OIVX        *don_ov_sse    = NULL;
  P7_OIVX        *don_ov_avx    = NULL;
  float           signal_scores[p7S_SPLICE_SIGNALS];
  P7_TRACE       *tr_search     = p7_trace_fs_Create();
  P7_TRACE       *tr_replay     = p7_trace_fs_Create();
  float           sc_sse, sc_avx;
  float           sc_search, sc_replay;
  int             intron_total  = intron_len + 4;
  int             L_amino, L_dna_total;
  int             k_start, k_end, sub_M;
  int             i, j;
  int             min_intron    = 13;
  int             n             = N;

  if (!esl_cpu_has_avx()) goto cleanup;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, M, p7_UNILOCAL);

  om_tr_sse = p7_fs_oprofile_Create_sse(M, abcAA, 1);
  om_tr_avx = p7_fs_oprofile_Create_avx(M, abcAA, 1);
  p7_fs_oprofile_Convert_Log_sse(gm_tr, om_tr_sse);
  p7_fs_oprofile_Convert_Log_avx(gm_tr, om_tr_avx);

  ox_sse     = p7_omx_Create_dpf_sse(M, M * 3, M * 3, p7X_NSCELLS);
  ox_avx     = p7_omx_Create_dpf_avx(M, M * 3, M * 3, p7X_NSCELLS);
  acc_ov_sse = p7_oivx_Create_sse(M, SPLICE_ROWS);
  acc_ov_avx = p7_oivx_Create_avx(M, SPLICE_ROWS);
  don_ov_sse = p7_oivx_Create_sse(M, SIGNAL_MEM_SIZE);
  don_ov_avx = p7_oivx_Create_avx(M, SIGNAL_MEM_SIZE);

  p7_SignalScores(signal_scores);
  p7_emit_SimpleConsensus(hmm, sq);
  L_amino = sq->n;

  while (n--)
    {
      k_start = k_end = 0;
      while (k_end <= k_start || k_end > L_amino) {
        k_start = esl_rnd_Roll(r, L_amino / 4) + 1;
        k_end   = esl_rnd_Roll(r, L_amino / 2) + k_start + L_amino / 4;
      }
      sub_M       = k_end - k_start + 1;
      L_dna_total = sub_M * 3 + intron_total;

      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) esl_fatal("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;

      j = 1;
      for (i = k_start; i <= k_start + sub_M / 2; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }
      dsq[j++] = 2; /* G */
      dsq[j++] = 3; /* T */
      for (i = 0; i < intron_len; i++) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = 0; /* A */
      dsq[j++] = 2; /* G */
      for (i = k_start + sub_M / 2 + 1; i <= k_end; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }

      p7_fs_oprofile_SubConvert_Log_sse(gm_tr, om_tr_sse, k_start, k_end);
      p7_fs_oprofile_SubConvert_Log_avx(gm_tr, om_tr_avx, k_start, k_end);
      p7_fs_oprofile_ReconfigLength_Log_sse(om_tr_sse, L_dna_total / 3);
      p7_fs_oprofile_ReconfigLength_Log_avx(om_tr_avx, L_dna_total / 3);

      p7_omx_GrowTo_dpf_sse(ox_sse, sub_M, L_dna_total, L_dna_total);
      p7_omx_GrowTo_dpf_avx(ox_avx, sub_M, L_dna_total, L_dna_total);
      p7_oivx_GrowTo_sse(acc_ov_sse, sub_M, SPLICE_ROWS);
      p7_oivx_GrowTo_avx(acc_ov_avx, sub_M, SPLICE_ROWS);
      p7_oivx_GrowTo_sse(don_ov_sse, sub_M, SIGNAL_MEM_SIZE);
      p7_oivx_GrowTo_avx(don_ov_avx, sub_M, SIGNAL_MEM_SIZE);

      p7_Viterbi_Spliced_sse(dsq, om_tr_sse, ox_sse, signal_scores, acc_ov_sse, don_ov_sse,
                             1, L_dna_total, min_intron, TRUE, TRUE);
      sc_sse = ox_sse->xmx[L_dna_total * p7X_NXCELLS + p7X_C];

      p7_Viterbi_Spliced_avx(dsq, om_tr_avx, ox_avx, signal_scores, acc_ov_avx, don_ov_avx,
                             1, L_dna_total, min_intron, TRUE, TRUE);
      sc_avx = ox_avx->xmx[L_dna_total * p7X_NXCELLS + p7X_C];

      if (sc_sse == -eslINFINITY && sc_avx == -eslINFINITY) continue;
      if (fabs(sc_sse - sc_avx) > 0.01)
        esl_fatal("%s: sse=%.4f avx=%.4f", msg, sc_sse, sc_avx);

      p7_Viterbi_SplicedTraceOld_sse      (dsq, ox_sse, gm_tr, signal_scores, tr_search, 1, L_dna_total, k_start, k_end, min_intron, &sc_search, NULL);
      p7_Viterbi_SplicedTrace_sse(dsq, ox_sse, gm_tr, signal_scores, tr_replay, 1, L_dna_total, k_start, k_end, min_intron, &sc_replay, NULL);
      utest_compare_traces(tr_search, sc_search, tr_replay, sc_replay, msg);
      p7_trace_Reuse(tr_search);
      p7_trace_Reuse(tr_replay);

      p7_Viterbi_SplicedTraceOld_avx      (dsq, ox_avx, gm_tr, signal_scores, tr_search, 1, L_dna_total, k_start, k_end, min_intron, &sc_search, NULL);
      p7_Viterbi_SplicedTrace_avx(dsq, ox_avx, gm_tr, signal_scores, tr_replay, 1, L_dna_total, k_start, k_end, min_intron, &sc_replay, NULL);
      utest_compare_traces(tr_search, sc_search, tr_replay, sc_replay, msg);
      p7_trace_Reuse(tr_search);
      p7_trace_Reuse(tr_replay);
    }

 cleanup:
  if (dsq)        free(dsq);
  if (ox_sse)     p7_omx_Destroy_sse(ox_sse);
  if (ox_avx)     p7_omx_Destroy_avx(ox_avx);
  if (acc_ov_sse) p7_oivx_Destroy_sse(acc_ov_sse);
  if (acc_ov_avx) p7_oivx_Destroy_avx(acc_ov_avx);
  if (don_ov_sse) p7_oivx_Destroy_sse(don_ov_sse);
  if (don_ov_avx) p7_oivx_Destroy_avx(don_ov_avx);
  if (om_tr_sse)  p7_fs_oprofile_Destroy_sse(om_tr_sse);
  if (om_tr_avx)  p7_fs_oprofile_Destroy_avx(om_tr_avx);
  if (hmm)        p7_hmm_Destroy(hmm);
  p7_trace_fs_Destroy(tr_search);
  p7_trace_fs_Destroy(tr_replay);
  esl_sq_Destroy(sq);
  p7_profile_fs_Destroy(gm_tr);
}
#endif /* eslENABLE_SSE && eslENABLE_AVX */

#ifdef eslENABLE_AVX
/* utest_chk_avx()
 *
 * For sequences with a simulated intron, fill with
 * p7_Viterbi_Spliced_avx() and with p7_Viterbi_SplicedChk_avx() at
 * several block sizes, in all three fill modes, and require the same
 * special states on every row. If the final C state is reachable, trace
 * both with p7_Viterbi_SplicedTrace_avx() and
 * p7_Viterbi_SplicedTraceChk_avx() and require the same trace, Viterbi
 * score and max_drop. Skipped if AVX is not available at runtime.
 */
static void
utest_chk_avx(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_CODONTABLE *codon_table,
              int M, int N, int intron_len)
{
  char            msg[]        = "utest_chk_avx failed";
  P7_HMM         *hmm          = NULL;
  P7_FS_PROFILE  *gm_tr        = p7_profile_fs_Create(M, abcAA, 1);
  P7_FS_OPROFILE *om_tr        = NULL;
  ESL_SQ         *sq           = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq          = NULL;
  P7_OMX         *ox           = NULL;
  P7_OIVX        *acc_ov       = NULL;
  P7_OIVX        *don_ov       = NULL;
  P7_TRACE       *tr_full      = p7_trace_fs_Create();
  P7_TRACE       *tr_chk       = p7_trace_fs_Create();
  float           signal_scores[p7S_SPLICE_SIGNALS];
  int             modes[3][2]  = { { TRUE, TRUE }, { TRUE, FALSE }, { FALSE, TRUE } };
  int             min_intron   = 13;
  int             blocks[3]    = { 13 + 6, 64, SPLICE_CHECKPOINT_ROWS };
  int             intron_total = intron_len + 4;
  int             L_amino, L_dna_total;
  int             k_start, k_end, sub_M;
  int             i, j, md, bi, s;
  int             n = N;

  if (!esl_cpu_has_avx()) goto cleanup;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, M, p7_UNILOCAL);
  om_tr  = p7_fs_oprofile_Create_avx(M, abcAA, 1);
  p7_fs_oprofile_Convert_Log_avx(gm_tr, om_tr);
  ox     = p7_omx_Create_dpf_avx(M, M * 3, M * 3, p7X_NSCELLS);
  acc_ov = p7_oivx_Create_avx(M, SPLICE_ROWS);
  don_ov = p7_oivx_Create_avx(M, SIGNAL_MEM_SIZE);

  p7_SignalScores(signal_scores);
  p7_emit_SimpleConsensus(hmm, sq);
  L_amino = sq->n;

  while (n--)
    {
      k_start = k_end = 0;
      while (k_end <= k_start || k_end > L_amino) {
        k_start = esl_rnd_Roll(r, L_amino / 4) + 1;
        k_end   = esl_rnd_Roll(r, L_amino / 2) + k_start + L_amino / 4;
      }
      sub_M       = k_end - k_start + 1;
      L_dna_total = sub_M * 3 + intron_total;

      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) esl_fatal("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;

      j = 1;
      for (i = k_start; i <= k_start + sub_M / 2; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }
      dsq[j++] = 2; /* G */
      dsq[j++] = 3; /* T */
      for (i = 0; i < intron_len; i++) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = 0; /* A */
      dsq[j++] = 2; /* G */
      for (i = k_start + sub_M / 2 + 1; i <= k_end; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }

      p7_fs_oprofile_SubConvert_Log_avx(gm_tr, om_tr, k_start, k_end);
      p7_fs_oprofile_ReconfigLength_Log_avx(om_tr, L_dna_total / 3);
      p7_fs_ReconfigLength(gm_tr, L_dna_total / 3);
      p7_omx_GrowTo_dpf_avx(ox, sub_M, L_dna_total, L_dna_total);
      p7_oivx_GrowTo_avx(acc_ov, sub_M, SPLICE_ROWS);
      p7_oivx_GrowTo_avx(don_ov, sub_M, SIGNAL_MEM_SIZE);

      for (md = 0; md < 3; md++)
        {
          p7_Viterbi_Spliced_avx(dsq, om_tr, ox, signal_scores, acc_ov, don_ov, 1, L_dna_total, min_intron, modes[md][0], modes[md][1]);

          for (bi = 0; bi < 3; bi++)
            {
              P7_CHKOMX  *cx  = p7_chkomx_Create_avx(sub_M, L_dna_total, blocks[bi], min_intron + 6, p7X_NSCELLS);
              P7_CHKOIVX *cv  = p7_chkoivx_Create_avx(sub_M, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, L_dna_total / blocks[bi]);
              P7_OIVX    *acc = p7_oivx_Create_avx(sub_M, SPLICE_ROWS);
              P7_OIVX    *don = p7_oivx_Create_avx(sub_M, SIGNAL_MEM_SIZE);
              float       sc_full, sc_chk, drop_full, drop_chk;

              if (cx == NULL || cv == NULL || acc == NULL || don == NULL) esl_fatal("%s: allocation failed", msg);
              if (p7_Viterbi_SplicedChk_avx(dsq, om_tr, cx, cv, signal_scores, acc, don, 1, L_dna_total, min_intron, modes[md][0], modes[md][1]) != eslOK)
                esl_fatal("%s: checkpointed fill failed", msg);

              for (i = 0; i <= L_dna_total; i++)
                for (s = 0; s < p7X_NXCELLS; s++)
                  if (cx->xmx[i*p7X_NXCELLS+s] != ox->xmx[i*p7X_NXCELLS+s])
                    esl_fatal("%s: special state %d differs at row %d (block size %d, mode %d)", msg, s, i, blocks[bi], md);

              if (ox->xmx[L_dna_total*p7X_NXCELLS+p7X_C] != -eslINFINITY)
                {
                  if (p7_Viterbi_SplicedTrace_avx(dsq, ox, gm_tr, signal_scores, tr_full, 1, L_dna_total, k_start, k_end, min_intron, &sc_full, &drop_full) != eslOK)
                    esl_fatal("%s: replay trace failed", msg);
                  if (p7_Viterbi_SplicedTraceChk_avx(dsq, om_tr, cx, cv, gm_tr, signal_scores, acc, don, tr_chk, 1, L_dna_total, k_start, k_end,
                                                     min_intron, modes[md][0], modes[md][1], &sc_chk, &drop_chk) != eslOK)
                    esl_fatal("%s: checkpointed trace failed (block size %d, mode %d)", msg, blocks[bi], md);
                  utest_compare_traces(tr_full, sc_full, tr_chk, sc_chk, msg);
                  if (drop_full != drop_chk) esl_fatal("%s: max_drop %f != %f (block size %d, mode %d)", msg, drop_chk, drop_full, blocks[bi], md);
                  p7_trace_Reuse(tr_full);
                  p7_trace_Reuse(tr_chk);
                }

              p7_chkomx_Destroy_avx(cx);
              p7_chkoivx_Destroy_avx(cv);
              p7_oivx_Destroy_avx(acc);
              p7_oivx_Destroy_avx(don);
            }
        }
    }

 cleanup:
  if (dsq)    free(dsq);
  if (ox)     p7_omx_Destroy_avx(ox);
  if (acc_ov) p7_oivx_Destroy_avx(acc_ov);
  if (don_ov) p7_oivx_Destroy_avx(don_ov);
  if (om_tr)  p7_fs_oprofile_Destroy_avx(om_tr);
  if (hmm)    p7_hmm_Destroy(hmm);
  p7_trace_fs_Destroy(tr_full);
  p7_trace_fs_Destroy(tr_chk);
  esl_sq_Destroy(sq);
  p7_profile_fs_Destroy(gm_tr);
}
#endif /* eslENABLE_AVX */

#ifdef eslENABLE_SSE
/* utest_chk_sse()
 *
 * For sequences with a simulated intron, fill with
 * p7_Viterbi_Spliced_sse() and with p7_Viterbi_SplicedChk_sse() at
 * several block sizes, in all three fill modes, and require the same
 * special states on every row. If the final C state is reachable, trace
 * both with p7_Viterbi_SplicedTrace_sse() and
 * p7_Viterbi_SplicedTraceChk_sse() and require the same trace, Viterbi
 * score and max_drop.
 */
static void
utest_chk_sse(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_CODONTABLE *codon_table,
              int M, int N, int intron_len)
{
  char            msg[]        = "utest_chk_sse failed";
  P7_HMM         *hmm          = NULL;
  P7_FS_PROFILE  *gm_tr        = p7_profile_fs_Create(M, abcAA, 1);
  P7_FS_OPROFILE *om_tr        = NULL;
  ESL_SQ         *sq           = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq          = NULL;
  P7_OMX         *ox           = NULL;
  P7_OIVX        *acc_ov       = NULL;
  P7_OIVX        *don_ov       = NULL;
  P7_TRACE       *tr_full      = p7_trace_fs_Create();
  P7_TRACE       *tr_chk       = p7_trace_fs_Create();
  float           signal_scores[p7S_SPLICE_SIGNALS];
  int             modes[3][2]  = { { TRUE, TRUE }, { TRUE, FALSE }, { FALSE, TRUE } };
  int             min_intron   = 13;
  int             blocks[3]    = { 13 + 6, 64, SPLICE_CHECKPOINT_ROWS };
  int             intron_total = intron_len + 4;
  int             L_amino, L_dna_total;
  int             k_start, k_end, sub_M;
  int             i, j, md, bi, s;
  int             n = N;

  if (!esl_cpu_has_sse()) goto cleanup;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, M, p7_UNILOCAL);
  om_tr  = p7_fs_oprofile_Create_sse(M, abcAA, 1);
  p7_fs_oprofile_Convert_Log_sse(gm_tr, om_tr);
  ox     = p7_omx_Create_dpf_sse(M, M * 3, M * 3, p7X_NSCELLS);
  acc_ov = p7_oivx_Create_sse(M, SPLICE_ROWS);
  don_ov = p7_oivx_Create_sse(M, SIGNAL_MEM_SIZE);

  p7_SignalScores(signal_scores);
  p7_emit_SimpleConsensus(hmm, sq);
  L_amino = sq->n;

  while (n--)
    {
      k_start = k_end = 0;
      while (k_end <= k_start || k_end > L_amino) {
        k_start = esl_rnd_Roll(r, L_amino / 4) + 1;
        k_end   = esl_rnd_Roll(r, L_amino / 2) + k_start + L_amino / 4;
      }
      sub_M       = k_end - k_start + 1;
      L_dna_total = sub_M * 3 + intron_total;

      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) esl_fatal("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;

      j = 1;
      for (i = k_start; i <= k_start + sub_M / 2; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }
      dsq[j++] = 2; /* G */
      dsq[j++] = 3; /* T */
      for (i = 0; i < intron_len; i++) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = 0; /* A */
      dsq[j++] = 2; /* G */
      for (i = k_start + sub_M / 2 + 1; i <= k_end; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }

      p7_fs_oprofile_SubConvert_Log_sse(gm_tr, om_tr, k_start, k_end);
      p7_fs_oprofile_ReconfigLength_Log_sse(om_tr, L_dna_total / 3);
      p7_fs_ReconfigLength(gm_tr, L_dna_total / 3);
      p7_omx_GrowTo_dpf_sse(ox, sub_M, L_dna_total, L_dna_total);
      p7_oivx_GrowTo_sse(acc_ov, sub_M, SPLICE_ROWS);
      p7_oivx_GrowTo_sse(don_ov, sub_M, SIGNAL_MEM_SIZE);

      for (md = 0; md < 3; md++)
        {
          p7_Viterbi_Spliced_sse(dsq, om_tr, ox, signal_scores, acc_ov, don_ov, 1, L_dna_total, min_intron, modes[md][0], modes[md][1]);

          for (bi = 0; bi < 3; bi++)
            {
              P7_CHKOMX  *cx  = p7_chkomx_Create_sse(sub_M, L_dna_total, blocks[bi], min_intron + 6, p7X_NSCELLS);
              P7_CHKOIVX *cv  = p7_chkoivx_Create_sse(sub_M, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, L_dna_total / blocks[bi]);
              P7_OIVX    *acc = p7_oivx_Create_sse(sub_M, SPLICE_ROWS);
              P7_OIVX    *don = p7_oivx_Create_sse(sub_M, SIGNAL_MEM_SIZE);
              float       sc_full, sc_chk, drop_full, drop_chk;

              if (cx == NULL || cv == NULL || acc == NULL || don == NULL) esl_fatal("%s: allocation failed", msg);
              if (p7_Viterbi_SplicedChk_sse(dsq, om_tr, cx, cv, signal_scores, acc, don, 1, L_dna_total, min_intron, modes[md][0], modes[md][1]) != eslOK)
                esl_fatal("%s: checkpointed fill failed", msg);

              for (i = 0; i <= L_dna_total; i++)
                for (s = 0; s < p7X_NXCELLS; s++)
                  if (cx->xmx[i*p7X_NXCELLS+s] != ox->xmx[i*p7X_NXCELLS+s])
                    esl_fatal("%s: special state %d differs at row %d (block size %d, mode %d)", msg, s, i, blocks[bi], md);

              if (ox->xmx[L_dna_total*p7X_NXCELLS+p7X_C] != -eslINFINITY)
                {
                  if (p7_Viterbi_SplicedTrace_sse(dsq, ox, gm_tr, signal_scores, tr_full, 1, L_dna_total, k_start, k_end, min_intron, &sc_full, &drop_full) != eslOK)
                    esl_fatal("%s: replay trace failed", msg);
                  if (p7_Viterbi_SplicedTraceChk_sse(dsq, om_tr, cx, cv, gm_tr, signal_scores, acc, don, tr_chk, 1, L_dna_total, k_start, k_end,
                                                     min_intron, modes[md][0], modes[md][1], &sc_chk, &drop_chk) != eslOK)
                    esl_fatal("%s: checkpointed trace failed (block size %d, mode %d)", msg, blocks[bi], md);
                  utest_compare_traces(tr_full, sc_full, tr_chk, sc_chk, msg);
                  if (drop_full != drop_chk) esl_fatal("%s: max_drop %f != %f (block size %d, mode %d)", msg, drop_chk, drop_full, blocks[bi], md);
                  p7_trace_Reuse(tr_full);
                  p7_trace_Reuse(tr_chk);
                }

              p7_chkomx_Destroy_sse(cx);
              p7_chkoivx_Destroy_sse(cv);
              p7_oivx_Destroy_sse(acc);
              p7_oivx_Destroy_sse(don);
            }
        }
    }

 cleanup:
  if (dsq)    free(dsq);
  if (ox)     p7_omx_Destroy_sse(ox);
  if (acc_ov) p7_oivx_Destroy_sse(acc_ov);
  if (don_ov) p7_oivx_Destroy_sse(don_ov);
  if (om_tr)  p7_fs_oprofile_Destroy_sse(om_tr);
  if (hmm)    p7_hmm_Destroy(hmm);
  p7_trace_fs_Destroy(tr_full);
  p7_trace_fs_Destroy(tr_chk);
  esl_sq_Destroy(sq);
  p7_profile_fs_Destroy(gm_tr);
}
#endif /* eslENABLE_SSE */

#endif /*p7VITERBI_SP_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/



/*****************************************************************
 * 5. Test driver.
 *****************************************************************/
#ifdef p7VITERBI_SP_TESTDRIVE
/*
  gcc -g -Wall -msse2 -std=gnu99 -o viterbi_sp_utest -I.. -L.. -I../../easel -L../../easel -Dp7VITERBI_SP_TESTDRIVE viterbi_sp.c -lhmmer -leasel -lm
  ./viterbi_sp_utest
*/
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_avx.h"

static ESL_OPTIONS options[] = {
  /* name    type           default env range toggles reqs incomp help                                       docgroup */
  { "-h",  eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",  eslARG_INT,     "42", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M",  eslARG_INT,    "200", NULL,"n>=5",NULL,NULL, NULL, "size of random models to sample",                0 },
  { "-N",  eslARG_INT,     "20", NULL, NULL, NULL, NULL, NULL, "number of random sequences to sample",           0 },
  { "-I",  eslARG_INT,    "500", NULL,"n>0", NULL, NULL, NULL, "simulated intron length",                        0 }, 
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for p7_Viterbi_Spliced()";

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
  int             I      = esl_opt_GetInteger(go, "-I");

  p7_FLogsumInit();
  impl_Init();

  utest_viterbi_sp(r, abcAA, abcDNA, gcode, bgAA, ct, M, N, I);
#if defined(eslENABLE_SSE) && defined(eslENABLE_AVX)
  utest_sse_vs_avx_viterbi_sp(r, abcAA, abcDNA, gcode, bgAA, ct, M, N, I);
#endif
#ifdef eslENABLE_AVX
  utest_chk_avx(r, abcAA, gcode, bgAA, ct, M, N, I);
#endif
#ifdef eslENABLE_SSE
  utest_chk_sse(r, abcAA, gcode, bgAA, ct, M, N, I);
#endif

  p7_bg_Destroy(bgAA);
  p7_codontable_Destroy(ct);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7VITERBI_SP_TESTDRIVE*/
/*----------------- end, test driver ----------------------------*/
