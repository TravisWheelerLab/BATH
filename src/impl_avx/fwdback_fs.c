/* Frameshift Forward/Backward dispatcher.
 * Provides the non-suffixed frameshift Forward/Backward functions declared
 * in impl_avx.h.  Delegates to the fastest available ISA implementation
 * at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Forward declaration of dispatcher */
static int p7_ForwardParser_Frameshift_3Codons_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                                          P7_OMX *ox, P7_OIVX *ov, float *opt_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_ForwardParser_Frameshift_3Codons)(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                           P7_OMX *ox, P7_OIVX *ov, float *opt_sc) = p7_ForwardParser_Frameshift_3Codons_Dispatcher;

/* Function:  p7_ForwardParser_Frameshift_3Codons_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for 3-codon frameshift Forward parser.
 *
 * Returns:   <eslOK> on success.  <*opt_sc> is the log Forward score in nats.
 * Throws:    <eslEINVAL> if allocation is too small.
 *            <eslEMEM> on reallocation failure.
 */
static int
p7_ForwardParser_Frameshift_3Codons_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                               P7_OMX *ox, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_ForwardParser_Frameshift_3Codons = p7_ForwardParser_Frameshift_3Codons_avx512; return p7_ForwardParser_Frameshift_3Codons_avx512(dsq, L, om_fs, ox, ov, opt_sc); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_ForwardParser_Frameshift_3Codons = p7_ForwardParser_Frameshift_3Codons_avx;    return p7_ForwardParser_Frameshift_3Codons_avx(dsq, L, om_fs, ox, ov, opt_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_ForwardParser_Frameshift_3Codons = p7_ForwardParser_Frameshift_3Codons_sse;
  return p7_ForwardParser_Frameshift_3Codons_sse(dsq, L, om_fs, ox, ov, opt_sc);
#else
  p7_Die("p7_ForwardParser_Frameshift_3Codons: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_BackwardParser_Frameshift_3Codons_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                                           const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_BackwardParser_Frameshift_3Codons)(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                            const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc) = p7_BackwardParser_Frameshift_3Codons_Dispatcher;

/* Function:  p7_BackwardParser_Frameshift_3Codons_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for 3-codon frameshift Backward parser.
 *
 * Returns:   <eslOK> on success.  <*opt_sc> is the log Backward score.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
static int
p7_BackwardParser_Frameshift_3Codons_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                                const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_BackwardParser_Frameshift_3Codons = p7_BackwardParser_Frameshift_3Codons_avx512; return p7_BackwardParser_Frameshift_3Codons_avx512(dsq, L, om_fs, fwd, bck, ov, opt_sc); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_BackwardParser_Frameshift_3Codons = p7_BackwardParser_Frameshift_3Codons_avx;    return p7_BackwardParser_Frameshift_3Codons_avx(dsq, L, om_fs, fwd, bck, ov, opt_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_BackwardParser_Frameshift_3Codons = p7_BackwardParser_Frameshift_3Codons_sse;
  return p7_BackwardParser_Frameshift_3Codons_sse(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#else
  p7_Die("p7_BackwardParser_Frameshift_3Codons: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_ForwardParser_Frameshift_5Codons_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                                          P7_OMX *ox, P7_OIVX *ov, float *opt_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_ForwardParser_Frameshift_5Codons)(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                           P7_OMX *ox, P7_OIVX *ov, float *opt_sc) = p7_ForwardParser_Frameshift_5Codons_Dispatcher;

/* Function:  p7_ForwardParser_Frameshift_5Codons_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for 5-codon frameshift Forward parser.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
static int
p7_ForwardParser_Frameshift_5Codons_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                               P7_OMX *ox, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_ForwardParser_Frameshift_5Codons = p7_ForwardParser_Frameshift_5Codons_avx512; return p7_ForwardParser_Frameshift_5Codons_avx512(dsq, L, om_fs, ox, ov, opt_sc); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_ForwardParser_Frameshift_5Codons = p7_ForwardParser_Frameshift_5Codons_avx;    return p7_ForwardParser_Frameshift_5Codons_avx(dsq, L, om_fs, ox, ov, opt_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_ForwardParser_Frameshift_5Codons = p7_ForwardParser_Frameshift_5Codons_sse;
  return p7_ForwardParser_Frameshift_5Codons_sse(dsq, L, om_fs, ox, ov, opt_sc);
#else
  p7_Die("p7_ForwardParser_Frameshift_5Codons: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_BackwardParser_Frameshift_5Codons_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                                           const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_BackwardParser_Frameshift_5Codons)(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                            const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc) = p7_BackwardParser_Frameshift_5Codons_Dispatcher;

/* Function:  p7_BackwardParser_Frameshift_5Codons_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for 5-codon frameshift Backward parser.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
static int
p7_BackwardParser_Frameshift_5Codons_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                                const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_BackwardParser_Frameshift_5Codons = p7_BackwardParser_Frameshift_5Codons_avx512; return p7_BackwardParser_Frameshift_5Codons_avx512(dsq, L, om_fs, fwd, bck, ov, opt_sc); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_BackwardParser_Frameshift_5Codons = p7_BackwardParser_Frameshift_5Codons_avx;    return p7_BackwardParser_Frameshift_5Codons_avx(dsq, L, om_fs, fwd, bck, ov, opt_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_BackwardParser_Frameshift_5Codons = p7_BackwardParser_Frameshift_5Codons_sse;
  return p7_BackwardParser_Frameshift_5Codons_sse(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#else
  p7_Die("p7_BackwardParser_Frameshift_5Codons: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_Forward_Frameshift_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                            P7_OMX *ox, P7_OIVX *ov, float *opt_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Forward_Frameshift)(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                             P7_OMX *ox, P7_OIVX *ov, float *opt_sc) = p7_Forward_Frameshift_Dispatcher;

/* Function:  p7_Forward_Frameshift_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for full-matrix 5-codon frameshift Forward.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
static int
p7_Forward_Frameshift_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                 P7_OMX *ox, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_Forward_Frameshift = p7_Forward_Frameshift_avx512; return p7_Forward_Frameshift_avx512(dsq, L, om_fs, ox, ov, opt_sc); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Forward_Frameshift = p7_Forward_Frameshift_avx;    return p7_Forward_Frameshift_avx(dsq, L, om_fs, ox, ov, opt_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_Forward_Frameshift = p7_Forward_Frameshift_sse;
  return p7_Forward_Frameshift_sse(dsq, L, om_fs, ox, ov, opt_sc);
#else
  p7_Die("p7_Forward_Frameshift: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_Backward_Frameshift_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                             const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Backward_Frameshift)(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                              const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc) = p7_Backward_Frameshift_Dispatcher;

/* Function:  p7_Backward_Frameshift_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for full-matrix 5-codon frameshift Backward.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
static int
p7_Backward_Frameshift_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                  const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_Backward_Frameshift = p7_Backward_Frameshift_avx512; return p7_Backward_Frameshift_avx512(dsq, L, om_fs, fwd, bck, ov, opt_sc); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Backward_Frameshift = p7_Backward_Frameshift_avx;    return p7_Backward_Frameshift_avx(dsq, L, om_fs, fwd, bck, ov, opt_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_Backward_Frameshift = p7_Backward_Frameshift_sse;
  return p7_Backward_Frameshift_sse(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#else
  p7_Die("p7_Backward_Frameshift: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_BackwardDecoding_Frameshift_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                                     P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_BackwardDecoding_Frameshift)(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                      P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc) = p7_BackwardDecoding_Frameshift_Dispatcher;

/* Function:  p7_BackwardDecoding_Frameshift_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for the 5-codon frameshift Backward
 *            with posterior decoding: Backward kept in a ring of
 *            <PARSER_ROWS_BWD> rows of <bck>, posteriors written over the
 *            full Forward matrix <fwd>, as p7_Decoding_Frameshift() would
 *            produce them.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if the Backward score or the posterior scaling
 *            is out of range.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
static int
p7_BackwardDecoding_Frameshift_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                          P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_BackwardDecoding_Frameshift = p7_BackwardDecoding_Frameshift_avx;    return p7_BackwardDecoding_Frameshift_avx(dsq, L, om_fs, fwd, bck, ov, opt_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_BackwardDecoding_Frameshift = p7_BackwardDecoding_Frameshift_sse;
  return p7_BackwardDecoding_Frameshift_sse(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#else
  p7_Die("p7_BackwardDecoding_Frameshift: no SIMD implementation available");
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
#ifdef p7FWDBACK_FS_BENCHMARK
/*
   gcc -g -O3 -mavx2 -std=gnu99 -o benchmark-fwdback_fs -I.. -L.. -I../../easel -L../../easel -Dp7FWDBACK_FS_BENCHMARK fwdback_fs.c -lhmmer -leasel -lm
   icc -O3 -static -o benchmark-fwdback_fs -I.. -L.. -I../../easel -L../../easel -Dp7FWDBACK_FS_BENCHMARK fwdback_fs.c -lhmmer -leasel -lm
   ./benchmark-fwdback_fs <hmmfile>
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
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,   "1200", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs (nucleotides)",     0 },
  { "-N",        eslARG_INT,   "2000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "-B",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark Backward functions",              0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark Forward functions",               0 },
  { "-P",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark parser (linear memory) versions", 0 },
  { "-U",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark full matrix versions",            0 },
  { "-T",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  "-P", NULL, "only benchmark 3-codon length parser",          0 },
  { "-V",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  "-P", NULL, "only benchmark 5-codon length parser",          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for AVX frameshift Forward/Backward";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA   = NULL;
  ESL_ALPHABET   *abcDNA  = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bgDNA   = NULL;
  P7_BG          *bgAA    = NULL;
  P7_FS_PROFILE  *gm_fs5  = NULL;
  P7_FS_PROFILE  *gm_fs3  = NULL;
  P7_FS_OPROFILE *om_fs5  = NULL;
  P7_FS_OPROFILE *om_fs3  = NULL;
  P7_OMX         *fwd_p   = NULL;  /* parser fwd: linear memory, 3-cell layout  */
  P7_OMX         *bck_p   = NULL;  /* parser bck: linear memory, 3-cell layout  */
  P7_OMX         *fwd     = NULL;  /* full fwd: O(ML), 8-cell FS layout          */
  P7_OMX         *bck     = NULL;  /* full bck: O(ML), 3-cell layout             */
  ESL_GENCODE    *gcode   = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)          != eslOK) p7_Fail("Failed to read HMM");

  abcDNA = esl_alphabet_Create(eslDNA);
  gcode  = esl_gencode_Create(abcDNA, abcAA);
  bgDNA  = p7_bg_Create(abcDNA);
  bgAA   = p7_bg_Create(abcAA);

  gm_fs5 = p7_profile_fs_Create(hmm->M, abcAA, p7P_5CODONS);
  gm_fs3 = p7_profile_fs_Create(hmm->M, abcAA, p7P_3CODONS);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_UNILOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs3, L/3, p7_UNILOCAL);

  om_fs5 = p7_fs_oprofile_Create(hmm->M, abcAA, p7P_5CODONS);
  om_fs3 = p7_fs_oprofile_Create(hmm->M, abcAA, p7P_3CODONS);
  p7_fs_oprofile_Convert(gm_fs5, om_fs5);
  p7_fs_oprofile_Convert(gm_fs3, om_fs3);
  p7_fs_oprofile_ReconfigLength(om_fs5, L/3);
  p7_fs_oprofile_ReconfigLength(om_fs3, L/3);

  fwd_p = p7_omx_Create(hmm->M, PARSER_ROWS_FWD, L);
  bck_p = p7_omx_Create(hmm->M, PARSER_ROWS_BWD, L);
  fwd   = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS_FS);
  bck   = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS);
  P7_OIVX *ov3 = p7_oivx_Create(hmm->M, p7P_3CODONS);
  P7_OIVX *ov5 = p7_oivx_Create(hmm->M, p7P_5CODONS);

  /* Baseline: time to generate sequences alone */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) {
    esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
    if(esl_opt_GetBoolean(go, "-B")) {
      if (! esl_opt_GetBoolean(go, "-U") && ! esl_opt_GetBoolean(go, "-V"))
        p7_ForwardParser_Frameshift_3Codons(dsq, L, om_fs3, fwd_p, ov3, &sc);
      if (! esl_opt_GetBoolean(go, "-U") && ! esl_opt_GetBoolean(go, "-T"))
        p7_ForwardParser_Frameshift_5Codons(dsq, L, om_fs5, fwd_p, ov5, &sc);
      if (! esl_opt_GetBoolean(go, "-P"))
        p7_Forward_Frameshift(dsq, L, om_fs5, fwd, ov5, &sc);
	}
  }
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);

      if (! esl_opt_GetBoolean(go, "-U") && ! esl_opt_GetBoolean(go, "-V"))
        p7_ForwardParser_Frameshift_3Codons(dsq, L, om_fs3, fwd_p, ov3, &sc);
      if (! esl_opt_GetBoolean(go, "-U") && ! esl_opt_GetBoolean(go, "-T"))
        p7_ForwardParser_Frameshift_5Codons(dsq, L, om_fs5, fwd_p, ov5, &sc);
      if (! esl_opt_GetBoolean(go, "-P"))
        p7_Forward_Frameshift(dsq, L, om_fs5, fwd, ov5, &sc);

      if (! esl_opt_GetBoolean(go, "-F")) {
        if (! esl_opt_GetBoolean(go, "-U") && ! esl_opt_GetBoolean(go, "-V"))
          p7_BackwardParser_Frameshift_3Codons(dsq, L, om_fs3, fwd_p, bck_p, ov3, NULL);
        if (! esl_opt_GetBoolean(go, "-U") && ! esl_opt_GetBoolean(go, "-T"))
          p7_BackwardParser_Frameshift_5Codons(dsq, L, om_fs5, fwd_p, bck_p, ov5, NULL);
        if (! esl_opt_GetBoolean(go, "-P"))
          p7_Backward_Frameshift(dsq, L, om_fs5, fwd, bck, ov5, NULL);
      }
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) hmm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   hmm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_omx_Destroy(bck);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bck_p);
  p7_omx_Destroy(fwd_p);
  p7_oivx_Destroy(ov3);
  p7_oivx_Destroy(ov5);
  p7_fs_oprofile_Destroy(om_fs5);
  p7_fs_oprofile_Destroy(om_fs3);
  p7_profile_fs_Destroy(gm_fs5);
  p7_profile_fs_Destroy(gm_fs3);
  p7_bg_Destroy(bgDNA);
  p7_bg_Destroy(bgAA);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7FWDBACK_FS_BENCHMARK*/
/*---------------- end, benchmark driver ----------------*/


/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef p7FWDBACK_FS_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"

/*
 * compare to GForward() scores.
 */
static void
utest_fwdbackfs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA, P7_CODONTABLE *codon_table, int M, int N)
{
  int  i,j;
  int  curr_L;
  char           *msg    = "forward/backward fs unit test failed";
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_fs3 = p7_profile_fs_Create(M, abcAA, 3);
  P7_FS_PROFILE  *gm_fs5 = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE *om_fs3 = p7_fs_oprofile_Create(M, abcAA, 3);
  P7_FS_OPROFILE *om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq    = NULL;
  P7_TRACE       *tr     = p7_trace_Create();
  P7_OMX         *oxf    = p7_omx_Create(M, PARSER_ROWS_FWD, M);
  P7_OMX         *oxb    = p7_omx_Create(M, PARSER_ROWS_BWD, M);
  P7_OMX         *fwd    = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS_FS);
  P7_OMX         *bck    = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS);
  P7_GMX         *fgx    = p7_gmx_Create(M, PARSER_ROWS_FWD, M, p7X_NSCELLS);
  P7_IVX         *iv3    = p7_ivx_Create(M, p7P_3CODONS);
  P7_IVX         *iv5    = p7_ivx_Create(M, p7P_5CODONS);
  P7_OIVX        *ov3    = p7_oivx_Create(M, p7P_3CODONS);
  P7_OIVX        *ov5    = p7_oivx_Create(M, p7P_5CODONS);
  float tolerance, generic_tolerance;
  float fsc3, bsc3;
  float fsc5, bsc5;
  float full_fsc, full_bsc;
  float generic_fsc3;
  float generic_fsc5;

  p7_FLogsumInit();
  if (p7_FLogsumError(-0.4, -0.5) > 0.0001) generic_tolerance = 1.0;  /* weaker test against generic   */
  else generic_tolerance = 0.001;   /* stronger test: FLogsum() is in slow exact mode. */
  tolerance = 0.0001;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig(hmm, bgAA, gm, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs3, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, M, p7_LOCAL);
  p7_fs_oprofile_Convert(gm_fs3, om_fs3);
  p7_fs_oprofile_ReconfigLength(om_fs3, M);
  p7_fs_oprofile_Convert(gm_fs5, om_fs5);
  p7_fs_oprofile_ReconfigLength(om_fs5, M);

  while (N--)
    {

	  p7_ProfileEmit(r, hmm, gm, bgAA, sq, tr);
      curr_L = sq->n*3;

	  if(dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) *(curr_L+2))) == NULL)  esl_fatal("malloc failed");

	  j = 1;
	  for(i = 1; i <= sq->n; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq+j);
        j+=3;
      }

      p7_fs_oprofile_ReconfigLength(om_fs3, sq->n);
      p7_fs_ReconfigLength(gm_fs3, sq->n);
      p7_fs_oprofile_ReconfigLength(om_fs5, sq->n);
      p7_fs_ReconfigLength(gm_fs5, sq->n);

      p7_omx_GrowTo(oxf, M, PARSER_ROWS_FWD, curr_L);
      p7_gmx_GrowTo(fgx, M, PARSER_ROWS_FWD, curr_L);

      /* 3-codon SSE vs scalar */
      if (p7_ForwardParser_Frameshift_3Codons(dsq, curr_L, om_fs3, oxf, ov3, &fsc3) == eslERANGE) continue;
	  p7_GForwardParser_Frameshift_3Codons(dsq, curr_L, gm_fs3, fgx, iv3, &generic_fsc3);

      if (fabs(fsc3-generic_fsc3) > generic_tolerance) esl_fatal(msg);

      p7_omx_GrowTo(oxb, M, PARSER_ROWS_BWD, curr_L);
      if (p7_BackwardParser_Frameshift_3Codons(dsq, curr_L, om_fs3, oxf, oxb, ov3, &bsc3) == eslERANGE) continue;

      if (fabs(fsc3-bsc3) > tolerance) esl_fatal(msg);

      /* 5-codon SSE vs scalar */
      p7_omx_GrowTo(oxf, M, PARSER_ROWS_FWD, curr_L);
      p7_gmx_GrowTo(fgx, M, PARSER_ROWS_FWD, curr_L);

      if (p7_ForwardParser_Frameshift_5Codons(dsq, curr_L, om_fs5, oxf, ov5, &fsc5) == eslERANGE) continue;
      p7_GForwardParser_Frameshift_5Codons(dsq, curr_L, gm_fs5, fgx, iv5, &generic_fsc5);

      if (fabs(fsc5-generic_fsc5) > generic_tolerance) esl_fatal(msg);

      p7_omx_GrowTo(oxb, M, PARSER_ROWS_BWD, curr_L);
      if (p7_BackwardParser_Frameshift_5Codons(dsq, curr_L, om_fs5, oxf, oxb, ov5, &bsc5) == eslERANGE) continue;

      if (fabs(fsc5-bsc5) > tolerance) esl_fatal(msg);

      p7_omx_GrowTo_dpf(fwd, M, curr_L, curr_L);
      if (p7_Forward_Frameshift(dsq, curr_L, om_fs5, fwd, ov5, &full_fsc)            == eslERANGE) continue;

      if (fabs(fsc5-full_fsc) > tolerance) esl_fatal(msg);

      p7_omx_GrowTo_dpf(bck, M, curr_L, curr_L);
      if (p7_Backward_Frameshift(dsq, curr_L, om_fs5, fwd, bck, ov5, &full_bsc)      == eslERANGE) continue;

      if (fabs(bsc5-full_bsc) > tolerance) esl_fatal(msg);
    }

  free(dsq);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(oxf);
  p7_omx_Destroy(oxb);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bck);
  p7_gmx_Destroy(fgx);
  p7_ivx_Destroy(iv3);
  p7_ivx_Destroy(iv5);
  p7_oivx_Destroy(ov3);
  p7_oivx_Destroy(ov5);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_fs3);
  p7_profile_fs_Destroy(gm_fs5);
  p7_fs_oprofile_Destroy(om_fs3);
  p7_fs_oprofile_Destroy(om_fs5);
}

/* utest_backward_decoding_fs()
 * Compare p7_BackwardDecoding_Frameshift() with p7_Backward_Frameshift()
 * followed by p7_Decoding_Frameshift() on sequences emitted from a
 * sampled model: the Backward scores, the Backward special states and
 * scale factors, and the posterior matrices must be identical. With
 * <use_sse>, the _sse versions are called directly, so an AVX CPU also
 * tests the SSE engine.
 */
static void
utest_backward_decoding_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_CODONTABLE *codon_table, int M, int N, int use_sse)
{
  char           *msg    = "backward decoding fs unit test failed";
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_fs5 = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE *om_fs5 = NULL;
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq    = NULL;
  P7_TRACE       *tr     = p7_trace_Create();
  P7_OMX         *fwd1   = NULL;
  P7_OMX         *bck1   = NULL;
  P7_OMX         *fwd2   = NULL;
  P7_OMX         *bck2   = NULL;
  P7_OIVX        *ov5    = NULL;
  P7_GMX         *gxp1   = p7_gmx_Create(M, M, M, p7G_NSCELLS_FS);
  P7_GMX         *gxp2   = p7_gmx_Create(M, M, M, p7G_NSCELLS_FS);
  float           bsc1, bsc2;
  int             s1, s2;
  int             curr_L, i, j, x;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig(hmm, bgAA, gm, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, M, p7_LOCAL);
  if (use_sse)
    {
      om_fs5 = p7_fs_oprofile_Create_sse(M, abcAA, 5);
      p7_fs_oprofile_Convert_sse(gm_fs5, om_fs5);
      fwd1 = p7_omx_Create_dpf_sse(M, M, M, p7X_NSCELLS_FS);
      bck1 = p7_omx_Create_dpf_sse(M, M, M, p7X_NSCELLS);
      fwd2 = p7_omx_Create_dpf_sse(M, M, M, p7X_NSCELLS_FS);
      bck2 = p7_omx_Create_dpf_sse(M, PARSER_ROWS_BWD, M, p7X_NSCELLS);
      ov5  = p7_oivx_Create_sse(M, p7P_5CODONS);
    }
  else
    {
      om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
      p7_fs_oprofile_Convert(gm_fs5, om_fs5);
      fwd1 = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS_FS);
      bck1 = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS);
      fwd2 = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS_FS);
      bck2 = p7_omx_Create_dpf(M, PARSER_ROWS_BWD, M, p7X_NSCELLS);
      ov5  = p7_oivx_Create(M, p7P_5CODONS);
    }

  while (N--)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, tr);
      curr_L = sq->n*3;
      if (curr_L < 15 || curr_L > 30000) { esl_sq_Reuse(sq); p7_trace_Reuse(tr); continue; }

      free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (curr_L+2))) == NULL) esl_fatal("malloc failed");
      for (i = 1, j = 1; i <= sq->n; i++, j += 3)
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq+j);

      p7_gmx_GrowTo(gxp1, M, curr_L, curr_L);
      p7_gmx_GrowTo(gxp2, M, curr_L, curr_L);
      if (use_sse)
        {
          p7_fs_oprofile_ReconfigLength_sse(om_fs5, sq->n);
          p7_omx_GrowTo_dpf_sse(fwd1, M, curr_L, curr_L);
          p7_omx_GrowTo_dpf_sse(bck1, M, curr_L, curr_L);
          p7_omx_GrowTo_dpf_sse(fwd2, M, curr_L, curr_L);
          p7_omx_GrowTo_dpf_sse(bck2, M, PARSER_ROWS_BWD, curr_L);
          if (p7_Forward_Frameshift_sse(dsq, curr_L, om_fs5, fwd1, ov5, NULL) != eslOK) continue;
          if (p7_Forward_Frameshift_sse(dsq, curr_L, om_fs5, fwd2, ov5, NULL) != eslOK) esl_fatal(msg);
          s1 = p7_Backward_Frameshift_sse(dsq, curr_L, om_fs5, fwd1, bck1, ov5, &bsc1);
          if (s1 == eslOK) s1 = p7_Decoding_Frameshift_sse(om_fs5, fwd1, bck1);
          s2 = p7_BackwardDecoding_Frameshift_sse(dsq, curr_L, om_fs5, fwd2, bck2, ov5, &bsc2);
          if (s1 != s2) esl_fatal("%s: return codes differ (%d %d)", msg, s1, s2);
          if (s1 != eslOK) continue;
          p7_omx_FDeconvert_sse(fwd1, gxp1);
          p7_omx_FDeconvert_sse(fwd2, gxp2);
        }
      else
        {
          p7_fs_oprofile_ReconfigLength(om_fs5, sq->n);
          p7_omx_GrowTo_dpf(fwd1, M, curr_L, curr_L);
          p7_omx_GrowTo_dpf(bck1, M, curr_L, curr_L);
          p7_omx_GrowTo_dpf(fwd2, M, curr_L, curr_L);
          p7_omx_GrowTo_dpf(bck2, M, PARSER_ROWS_BWD, curr_L);
          if (p7_Forward_Frameshift(dsq, curr_L, om_fs5, fwd1, ov5, NULL) != eslOK) continue;
          if (p7_Forward_Frameshift(dsq, curr_L, om_fs5, fwd2, ov5, NULL) != eslOK) esl_fatal(msg);
          s1 = p7_Backward_Frameshift(dsq, curr_L, om_fs5, fwd1, bck1, ov5, &bsc1);
          if (s1 == eslOK) s1 = p7_Decoding_Frameshift(om_fs5, fwd1, bck1);
          s2 = p7_BackwardDecoding_Frameshift(dsq, curr_L, om_fs5, fwd2, bck2, ov5, &bsc2);
          if (s1 != s2) esl_fatal("%s: return codes differ (%d %d)", msg, s1, s2);
          if (s1 != eslOK) continue;
          p7_omx_FDeconvert(fwd1, gxp1);
          p7_omx_FDeconvert(fwd2, gxp2);
        }

      if (bsc1 != bsc2)                       esl_fatal("%s: Backward scores differ (%f %f)", msg, bsc1, bsc2);
      if (p7_gmx_Compare(gxp2, gxp1, 0.0) != eslOK) esl_fatal("%s: posterior matrices differ", msg);
      for (i = 0; i <= curr_L; i++)
        for (x = 0; x < p7X_NXCELLS; x++)
          if (bck2->xmx[i*p7X_NXCELLS+x] != bck1->xmx[i*p7X_NXCELLS+x]) esl_fatal("%s: Backward special states differ", msg);

      esl_sq_Reuse(sq);
      p7_trace_Reuse(tr);
    }

  free(dsq);
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_hmm_Destroy(hmm);
  p7_gmx_Destroy(gxp1);
  p7_gmx_Destroy(gxp2);
  if (use_sse)
    {
      p7_omx_Destroy_sse(fwd1);
      p7_omx_Destroy_sse(bck1);
      p7_omx_Destroy_sse(fwd2);
      p7_omx_Destroy_sse(bck2);
      p7_oivx_Destroy_sse(ov5);
      p7_fs_oprofile_Destroy_sse(om_fs5);
    }
  else
    {
      p7_omx_Destroy(fwd1);
      p7_omx_Destroy(bck1);
      p7_omx_Destroy(fwd2);
      p7_omx_Destroy(bck2);
      p7_oivx_Destroy(ov5);
      p7_fs_oprofile_Destroy(om_fs5);
    }
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_fs5);
}

/* rowL_compare_specials()
 * Compare the special states of every row of SIMD Backward matrix <bck>,
 * with its row scale factors taken out, to those of the generic matrix
 * <gx>, for a sequence of length <L>. Fail if any differs by more than
 * <tol> nats.
 */
static void
rowL_compare_specials(const P7_OMX *bck, const P7_GMX *gx, int L, float tol, const char *msg, const char *what)
{
  int    sx[5]  = { p7X_N, p7X_J, p7X_C, p7X_B, p7X_E };
  int    gsx[5] = { p7G_N, p7G_J, p7G_C, p7G_B, p7G_E };
  double lsf    = 0.;
  float  a, g;
  int    i, s;

  for (i = L; i >= 1; i--)
    {
      lsf += log(bck->xmx[i*p7X_NXCELLS+p7X_SCALE]);
      for (s = 0; s < 5; s++)
        {
          a = bck->xmx[i*p7X_NXCELLS+sx[s]];
          g = gx->xmx[i*p7G_NXCELLS+gsx[s]];
          if (g == -eslINFINITY || a < 1e-30) continue;
          if (fabs(log(a) + lsf - g) > tol)
            esl_fatal("%s: %s, row %d of %d, special state %d: %.3f vs generic %.3f", msg, what, i, L, s, log(a) + lsf, g);
        }
    }
}

/* utest_backward_rowL_scale_fs()
 * The frameshift Backward variants (full matrix, 5-codon parser, 3-codon
 * parser) when Forward rescales one of the last rows. Each sequence is
 * the back-translated consensus of a model with sharpened match
 * emissions and M->M transitions, followed by 0 or 1 random nucleotides,
 * so the hit ends at row L or L-1, where Forward often rescales. The
 * special states of every row, with the Backward scale factors taken
 * out, must match the generic versions. Fails if Forward never rescaled
 * row L (full matrix) or row L or L-1 (3-codon parser), so the case is
 * exercised. With <use_sse>, the _sse versions are called directly, so
 * an AVX CPU also tests the SSE engine.
 */
static void
utest_backward_rowL_scale_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA,
                             P7_CODONTABLE *codon_table, int M, int nmodels, int use_sse)
{
  char            *msg    = "frameshift backward last-row scaling unit test failed";
  P7_FS_OPROFILE *(*om_create) (int, const ESL_ALPHABET *, int)            = use_sse ? p7_fs_oprofile_Create_sse         : p7_fs_oprofile_Create;
  int             (*om_convert)(const P7_FS_PROFILE *, P7_FS_OPROFILE *)    = use_sse ? p7_fs_oprofile_Convert_sse        : p7_fs_oprofile_Convert;
  int             (*om_reconf) (P7_FS_OPROFILE *, int)                     = use_sse ? p7_fs_oprofile_ReconfigLength_sse : p7_fs_oprofile_ReconfigLength;
  void            (*om_destroy)(P7_FS_OPROFILE *)                          = use_sse ? p7_fs_oprofile_Destroy_sse        : p7_fs_oprofile_Destroy;
  P7_OMX         *(*ox_create) (int, int, int, int)                        = use_sse ? p7_omx_Create_dpf_sse             : p7_omx_Create_dpf;
  void            (*ox_destroy)(P7_OMX *)                                  = use_sse ? p7_omx_Destroy_sse                : p7_omx_Destroy;
  P7_OIVX        *(*ov_create) (int, int)                                  = use_sse ? p7_oivx_Create_sse                : p7_oivx_Create;
  void            (*ov_destroy)(P7_OIVX *)                                 = use_sse ? p7_oivx_Destroy_sse               : p7_oivx_Destroy;
  int (*fwd_full)(const ESL_DSQ *, int, const P7_FS_OPROFILE *, P7_OMX *, P7_OIVX *, float *)                 = use_sse ? p7_Forward_Frameshift_sse                : p7_Forward_Frameshift;
  int (*bck_full)(const ESL_DSQ *, int, const P7_FS_OPROFILE *, const P7_OMX *, P7_OMX *, P7_OIVX *, float *) = use_sse ? p7_Backward_Frameshift_sse               : p7_Backward_Frameshift;
  int (*fwd_p5)  (const ESL_DSQ *, int, const P7_FS_OPROFILE *, P7_OMX *, P7_OIVX *, float *)                 = use_sse ? p7_ForwardParser_Frameshift_5Codons_sse  : p7_ForwardParser_Frameshift_5Codons;
  int (*bck_p5)  (const ESL_DSQ *, int, const P7_FS_OPROFILE *, const P7_OMX *, P7_OMX *, P7_OIVX *, float *) = use_sse ? p7_BackwardParser_Frameshift_5Codons_sse : p7_BackwardParser_Frameshift_5Codons;
  int (*fwd_p3)  (const ESL_DSQ *, int, const P7_FS_OPROFILE *, P7_OMX *, P7_OIVX *, float *)                 = use_sse ? p7_ForwardParser_Frameshift_3Codons_sse  : p7_ForwardParser_Frameshift_3Codons;
  int (*bck_p3)  (const ESL_DSQ *, int, const P7_FS_OPROFILE *, const P7_OMX *, P7_OMX *, P7_OIVX *, float *) = use_sse ? p7_BackwardParser_Frameshift_3Codons_sse : p7_BackwardParser_Frameshift_3Codons;
  int             Lmax   = 3*M + 1;
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq    = malloc(sizeof(ESL_DSQ) * (Lmax+2));
  P7_FS_PROFILE  *gm5    = p7_profile_fs_Create(M, abcAA, p7P_5CODONS);
  P7_FS_PROFILE  *gm3    = p7_profile_fs_Create(M, abcAA, p7P_3CODONS);
  P7_FS_OPROFILE *om5    = om_create(M, abcAA, p7P_5CODONS);
  P7_FS_OPROFILE *om3    = om_create(M, abcAA, p7P_3CODONS);
  P7_OMX         *fwd    = ox_create(M, Lmax, Lmax, p7X_NSCELLS_FS);
  P7_OMX         *bck    = ox_create(M, Lmax, Lmax, p7X_NSCELLS);
  P7_OMX         *pf     = ox_create(M, PARSER_ROWS_FWD, Lmax, p7G_NSCELLS);
  P7_OMX         *pb     = ox_create(M, PARSER_ROWS_BWD, Lmax, p7G_NSCELLS);
  P7_OIVX        *ov5    = ov_create(M, p7P_5CODONS);
  P7_OIVX        *ov3    = ov_create(M, p7P_3CODONS);
  P7_IVX         *iv5    = p7_ivx_Create(M, p7P_5CODONS);
  P7_IVX         *iv3    = p7_ivx_Create(M, p7P_3CODONS);
  P7_GMX         *gx     = p7_gmx_Create(M, Lmax, Lmax, p7G_NSCELLS);
  P7_HMM         *hmm    = NULL;
  float           tol;
  int             nfull  = 0;
  int             np3    = 0;
  int             n, tail, L, i, j, k, x;

  p7_FLogsumInit();
  tol = (p7_FLogsumError(-0.4, -0.5) > 0.0001) ? 1.0 : 0.01;

  for (n = 0; n < nmodels; n++)
    {
      if (p7_hmm_Sample(r, M, abcAA, &hmm) != eslOK) esl_fatal(msg);
      for (k = 1; k <= M; k++)
        {
          x = esl_vec_FArgMax(hmm->mat[k], abcAA->K);
          esl_vec_FSet(hmm->mat[k], abcAA->K, 0.15 / (abcAA->K - 1));
          hmm->mat[k][x] = 0.85;
          if (k < M) { hmm->t[k][p7H_MM] = 0.9; hmm->t[k][p7H_MI] = 0.05; hmm->t[k][p7H_MD] = 0.05; }
        }
      esl_sq_Reuse(sq);
      p7_emit_SimpleConsensus(hmm, sq);

      for (tail = 0; tail <= 1; tail++)
        {
          L = 3*M + tail;
          esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
          for (i = 1, j = 1; i <= M; i++, j += 3)
            p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq+j);

          if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm5, L/3, p7_LOCAL) != eslOK) esl_fatal(msg);
          if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm3, L/3, p7_LOCAL) != eslOK) esl_fatal(msg);
          om_convert(gm5, om5);  om_reconf(om5, L/3);
          om_convert(gm3, om3);  om_reconf(om3, L/3);

          if (fwd_full(dsq, L, om5, fwd,      ov5, NULL) != eslOK) esl_fatal(msg);
          if (bck_full(dsq, L, om5, fwd, bck, ov5, NULL) != eslOK) esl_fatal(msg);
          if (fwd->xmx[L*p7X_NXCELLS+p7X_SCALE] > 1.0f) nfull++;
          p7_GBackward_Frameshift(dsq, L, gm5, gx, iv5, NULL);
          rowL_compare_specials(bck, gx, L, tol, msg, "full matrix");

          if (fwd_p5(dsq, L, om5, pf,     ov5, NULL) != eslOK) esl_fatal(msg);
          if (bck_p5(dsq, L, om5, pf, pb, ov5, NULL) != eslOK) esl_fatal(msg);
          p7_GBackwardParser_Frameshift_5Codons(dsq, L, gm5, gx, iv5, NULL);
          rowL_compare_specials(pb, gx, L, tol, msg, "5-codon parser");

          if (fwd_p3(dsq, L, om3, pf,     ov3, NULL) != eslOK) esl_fatal(msg);
          if (pf->xmx[L*p7X_NXCELLS+p7X_SCALE] > 1.0f || pf->xmx[(L-1)*p7X_NXCELLS+p7X_SCALE] > 1.0f) np3++;
          if (bck_p3(dsq, L, om3, pf, pb, ov3, NULL) != eslOK) esl_fatal(msg);
          p7_GBackwardParser_Frameshift_3Codons(dsq, L, gm3, gx, iv3, NULL);
          rowL_compare_specials(pb, gx, L, tol, msg, "3-codon parser");
        }
      p7_hmm_Destroy(hmm);
    }
  if (nfull == 0 || np3 == 0) esl_fatal("%s: Forward never rescaled the last rows; test is too weak", msg);

  free(dsq);
  esl_sq_Destroy(sq);
  p7_profile_fs_Destroy(gm5);
  p7_profile_fs_Destroy(gm3);
  om_destroy(om5);
  om_destroy(om3);
  ox_destroy(fwd);
  ox_destroy(bck);
  ox_destroy(pf);
  ox_destroy(pb);
  ov_destroy(ov5);
  ov_destroy(ov3);
  p7_ivx_Destroy(iv5);
  p7_ivx_Destroy(iv3);
  p7_gmx_Destroy(gx);
}

#if defined(eslENABLE_SSE) && defined(eslENABLE_AVX)
/* utest_sse_vs_avx_fwdbackfs()
 *
 * Run the frameshift Forward/Backward variants for both SSE and AVX on the
 * same sequences and profiles, comparing scores with 0.01 nat tolerance.
 * Tests 3-codon parser, 5-codon parser, and full 5-codon Forward/Backward.
 * Skipped silently if AVX is not available at runtime.
 */
static void
utest_sse_vs_avx_fwdbackfs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
                            ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA,
                            P7_CODONTABLE *codon_table, int M, int N)
{
  char            msg[]       = "utest_sse_vs_avx_fwdbackfs failed";
  P7_HMM         *hmm         = NULL;
  P7_PROFILE     *gm          = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_fs3      = p7_profile_fs_Create(M, abcAA, 3);
  P7_FS_PROFILE  *gm_fs5      = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE *om_fs3_sse  = NULL;
  P7_FS_OPROFILE *om_fs3_avx  = NULL;
  P7_FS_OPROFILE *om_fs5_sse  = NULL;
  P7_FS_OPROFILE *om_fs5_avx  = NULL;
  ESL_SQ         *sq          = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq         = NULL;
  P7_OMX         *oxf_sse     = NULL;
  P7_OMX         *oxb_sse     = NULL;
  P7_OMX         *oxf_avx     = NULL;
  P7_OMX         *oxb_avx     = NULL;
  P7_OMX         *fwd_sse     = NULL;
  P7_OMX         *bck_sse     = NULL;
  P7_OMX         *fwd_avx     = NULL;
  P7_OMX         *bck_avx     = NULL;
  P7_OIVX        *ov3_sse     = NULL;
  P7_OIVX        *ov3_avx     = NULL;
  P7_OIVX        *ov5_sse     = NULL;
  P7_OIVX        *ov5_avx     = NULL;
  float           fsc3_sse, fsc3_avx;
  float           bsc3_sse, bsc3_avx;
  float           fsc5_sse, fsc5_avx;
  float           bsc5_sse, bsc5_avx;
  float           full_fsc_sse, full_fsc_avx;
  float           full_bsc_sse, full_bsc_avx;
  int             curr_L, i, j;
  int             n = N;

  if (!esl_cpu_has_avx()) goto cleanup;

  p7_FLogsumInit();

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig(hmm, bgAA, gm, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs3, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, M, p7_LOCAL);

  om_fs3_sse = p7_fs_oprofile_Create_sse(M, abcAA, 3);
  om_fs3_avx = p7_fs_oprofile_Create_avx(M, abcAA, 3);
  om_fs5_sse = p7_fs_oprofile_Create_sse(M, abcAA, 5);
  om_fs5_avx = p7_fs_oprofile_Create_avx(M, abcAA, 5);
  p7_fs_oprofile_Convert_sse(gm_fs3, om_fs3_sse);  p7_fs_oprofile_ReconfigLength_sse(om_fs3_sse, M);
  p7_fs_oprofile_Convert_avx(gm_fs3, om_fs3_avx);  p7_fs_oprofile_ReconfigLength_avx(om_fs3_avx, M);
  p7_fs_oprofile_Convert_sse(gm_fs5, om_fs5_sse);  p7_fs_oprofile_ReconfigLength_sse(om_fs5_sse, M);
  p7_fs_oprofile_Convert_avx(gm_fs5, om_fs5_avx);  p7_fs_oprofile_ReconfigLength_avx(om_fs5_avx, M);

  ov3_sse = p7_oivx_Create_sse(M, p7P_3CODONS);
  ov3_avx = p7_oivx_Create_avx(M, p7P_3CODONS);
  ov5_sse = p7_oivx_Create_sse(M, p7P_5CODONS);
  ov5_avx = p7_oivx_Create_avx(M, p7P_5CODONS);

  oxf_sse = p7_omx_Create_sse(M, PARSER_ROWS_FWD, M);
  oxf_avx = p7_omx_Create_avx(M, PARSER_ROWS_FWD, M);
  oxb_sse = p7_omx_Create_sse(M, PARSER_ROWS_BWD, M);
  oxb_avx = p7_omx_Create_avx(M, PARSER_ROWS_BWD, M);
  fwd_sse = p7_omx_Create_dpf_sse(M, M, M, p7X_NSCELLS_FS);
  fwd_avx = p7_omx_Create_dpf_avx(M, M, M, p7X_NSCELLS_FS);
  bck_sse = p7_omx_Create_dpf_sse(M, M, M, p7X_NSCELLS);
  bck_avx = p7_omx_Create_dpf_avx(M, M, M, p7X_NSCELLS);

  while (n--)
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

      p7_fs_oprofile_ReconfigLength_sse(om_fs3_sse, sq->n);
      p7_fs_oprofile_ReconfigLength_avx(om_fs3_avx, sq->n);
      p7_fs_oprofile_ReconfigLength_sse(om_fs5_sse, sq->n);
      p7_fs_oprofile_ReconfigLength_avx(om_fs5_avx, sq->n);

      p7_omx_GrowTo_sse(oxf_sse, M, PARSER_ROWS_FWD, curr_L);
      p7_omx_GrowTo_avx(oxf_avx, M, PARSER_ROWS_FWD, curr_L);
      p7_omx_GrowTo_sse(oxb_sse, M, PARSER_ROWS_BWD, curr_L);
      p7_omx_GrowTo_avx(oxb_avx, M, PARSER_ROWS_BWD, curr_L);

      /* 3-codon parser */
      if (p7_ForwardParser_Frameshift_3Codons_sse(dsq, curr_L, om_fs3_sse, oxf_sse, ov3_sse, &fsc3_sse) == eslERANGE) continue;
      if (p7_ForwardParser_Frameshift_3Codons_avx(dsq, curr_L, om_fs3_avx, oxf_avx, ov3_avx, &fsc3_avx) == eslERANGE) continue;
      if (fabs(fsc3_sse - fsc3_avx) > 0.01) esl_fatal("%s: 3cod Fwd sse=%.4f avx=%.4f", msg, fsc3_sse, fsc3_avx);

      if (p7_BackwardParser_Frameshift_3Codons_sse(dsq, curr_L, om_fs3_sse, oxf_sse, oxb_sse, ov3_sse, &bsc3_sse) == eslERANGE) continue;
      if (p7_BackwardParser_Frameshift_3Codons_avx(dsq, curr_L, om_fs3_avx, oxf_avx, oxb_avx, ov3_avx, &bsc3_avx) == eslERANGE) continue;
      if (fabs(bsc3_sse - bsc3_avx) > 0.01) esl_fatal("%s: 3cod Bck sse=%.4f avx=%.4f", msg, bsc3_sse, bsc3_avx);

      /* 5-codon parser */
      p7_omx_GrowTo_sse(oxf_sse, M, PARSER_ROWS_FWD, curr_L);
      p7_omx_GrowTo_avx(oxf_avx, M, PARSER_ROWS_FWD, curr_L);
      if (p7_ForwardParser_Frameshift_5Codons_sse(dsq, curr_L, om_fs5_sse, oxf_sse, ov5_sse, &fsc5_sse) == eslERANGE) continue;
      if (p7_ForwardParser_Frameshift_5Codons_avx(dsq, curr_L, om_fs5_avx, oxf_avx, ov5_avx, &fsc5_avx) == eslERANGE) continue;
      if (fabs(fsc5_sse - fsc5_avx) > 0.01) esl_fatal("%s: 5cod Fwd sse=%.4f avx=%.4f", msg, fsc5_sse, fsc5_avx);

      p7_omx_GrowTo_sse(oxb_sse, M, PARSER_ROWS_BWD, curr_L);
      p7_omx_GrowTo_avx(oxb_avx, M, PARSER_ROWS_BWD, curr_L);
      if (p7_BackwardParser_Frameshift_5Codons_sse(dsq, curr_L, om_fs5_sse, oxf_sse, oxb_sse, ov5_sse, &bsc5_sse) == eslERANGE) continue;
      if (p7_BackwardParser_Frameshift_5Codons_avx(dsq, curr_L, om_fs5_avx, oxf_avx, oxb_avx, ov5_avx, &bsc5_avx) == eslERANGE) continue;
      if (fabs(bsc5_sse - bsc5_avx) > 0.01) esl_fatal("%s: 5cod Bck sse=%.4f avx=%.4f", msg, bsc5_sse, bsc5_avx);

      /* Full 5-codon Forward/Backward */
      p7_omx_GrowTo_dpf_sse(fwd_sse, M, curr_L, curr_L);
      p7_omx_GrowTo_dpf_avx(fwd_avx, M, curr_L, curr_L);
      if (p7_Forward_Frameshift_sse(dsq, curr_L, om_fs5_sse, fwd_sse, ov5_sse, &full_fsc_sse) == eslERANGE) continue;
      if (p7_Forward_Frameshift_avx(dsq, curr_L, om_fs5_avx, fwd_avx, ov5_avx, &full_fsc_avx) == eslERANGE) continue;
      if (fabs(full_fsc_sse - full_fsc_avx) > 0.01) esl_fatal("%s: full Fwd sse=%.4f avx=%.4f", msg, full_fsc_sse, full_fsc_avx);

      p7_omx_GrowTo_dpf_sse(bck_sse, M, curr_L, curr_L);
      p7_omx_GrowTo_dpf_avx(bck_avx, M, curr_L, curr_L);
      if (p7_Backward_Frameshift_sse(dsq, curr_L, om_fs5_sse, fwd_sse, bck_sse, ov5_sse, &full_bsc_sse) == eslERANGE) continue;
      if (p7_Backward_Frameshift_avx(dsq, curr_L, om_fs5_avx, fwd_avx, bck_avx, ov5_avx, &full_bsc_avx) == eslERANGE) continue;
      if (fabs(full_bsc_sse - full_bsc_avx) > 0.01) esl_fatal("%s: full Bck sse=%.4f avx=%.4f", msg, full_bsc_sse, full_bsc_avx);
    }

 cleanup:
  if (dsq)        free(dsq);
  if (oxf_sse)    p7_omx_Destroy_sse(oxf_sse);
  if (oxf_avx)    p7_omx_Destroy_avx(oxf_avx);
  if (oxb_sse)    p7_omx_Destroy_sse(oxb_sse);
  if (oxb_avx)    p7_omx_Destroy_avx(oxb_avx);
  if (fwd_sse)    p7_omx_Destroy_sse(fwd_sse);
  if (fwd_avx)    p7_omx_Destroy_avx(fwd_avx);
  if (bck_sse)    p7_omx_Destroy_sse(bck_sse);
  if (bck_avx)    p7_omx_Destroy_avx(bck_avx);
  if (ov3_sse)    p7_oivx_Destroy_sse(ov3_sse);
  if (ov3_avx)    p7_oivx_Destroy_avx(ov3_avx);
  if (ov5_sse)    p7_oivx_Destroy_sse(ov5_sse);
  if (ov5_avx)    p7_oivx_Destroy_avx(ov5_avx);
  if (om_fs3_sse) p7_fs_oprofile_Destroy_sse(om_fs3_sse);
  if (om_fs3_avx) p7_fs_oprofile_Destroy_avx(om_fs3_avx);
  if (om_fs5_sse) p7_fs_oprofile_Destroy_sse(om_fs5_sse);
  if (om_fs5_avx) p7_fs_oprofile_Destroy_avx(om_fs5_avx);
  if (hmm)        p7_hmm_Destroy(hmm);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_fs3);
  p7_profile_fs_Destroy(gm_fs5);
}
#endif /* eslENABLE_SSE && eslENABLE_AVX */

#endif /*p7FWDBACK_FS_TESTDRIVE*/
/*---------------------- end, unit tests ------------------------*/



/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7FWDBACK_FS_TESTDRIVE
/*
   gcc -g -Wall -msse2 -std=gnu99 -o fwdback_fs_utest -I.. -L.. -I../../easel -L../../easel -Dp7FWDBACK_FS_TESTDRIVE fwdback_fs.c -lhmmer -leasel -lm
   ./fwdback_fs_utest
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
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE Forward, Backward implementations";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_GENCODE    *gcode  = NULL;
  ESL_ALPHABET   *abcAA  = NULL;
  ESL_ALPHABET   *abcDNA = NULL;
  P7_BG          *bgAA   = NULL;
  P7_BG          *bgDNA  = NULL;
  P7_CODONTABLE  *ct     = NULL;
  int             M      = esl_opt_GetInteger(go, "-M");
  int             N      = esl_opt_GetInteger(go, "-N");

  impl_Init();

  if ((abcDNA = esl_alphabet_Create(eslDNA))      == NULL)  esl_fatal("failed to create alphabet");
  if ((bgDNA  = p7_bg_Create(abcDNA))             == NULL)  esl_fatal("failed to create null model");
  if ((abcAA  = esl_alphabet_Create(eslAMINO))    == NULL)  esl_fatal("failed to create alphabet");
  if ((bgAA   = p7_bg_Create(abcAA))              == NULL)  esl_fatal("failed to create null model");
  if ((gcode  = esl_gencode_Create(abcDNA,abcAA)) == NULL)  esl_fatal("failed to create gencode");
  if ((ct     = p7_codontable_Create(gcode))      == NULL)  esl_fatal("failed to create codon table");
  
  utest_fwdbackfs(r, abcAA, abcDNA, gcode, bgAA, bgDNA, ct, M, N);
#if defined(eslENABLE_SSE) && defined(eslENABLE_AVX)
  utest_sse_vs_avx_fwdbackfs(r, abcAA, abcDNA, gcode, bgAA, bgDNA, ct, M, N);
#endif
  utest_backward_decoding_fs(r, abcAA, gcode, bgAA, ct, M, N, FALSE);
  utest_backward_rowL_scale_fs(r, abcAA, abcDNA, gcode, bgAA, bgDNA, ct, M, 10, FALSE);
#ifdef eslENABLE_SSE
  if (esl_cpu_has_sse()) utest_backward_decoding_fs(r, abcAA, gcode, bgAA, ct, M, N, TRUE);
  if (esl_cpu_has_sse()) utest_backward_rowL_scale_fs(r, abcAA, abcDNA, gcode, bgAA, bgDNA, ct, M, 10, TRUE);
#endif

  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  p7_bg_Destroy(bgDNA);
  p7_bg_Destroy(bgAA);
  esl_gencode_Destroy(gcode);
  p7_codontable_Destroy(ct);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*p7FWDBACK_FS_TESTDRIVE*/


/*--------------------- end, test driver ------------------------*/



/*****************************************************************
 * 5. Example
 *****************************************************************/

