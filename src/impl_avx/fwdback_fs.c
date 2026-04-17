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

