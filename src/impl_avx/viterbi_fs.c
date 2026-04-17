/* Frameshift Viterbi dispatcher.
 * Provides the non-suffixed p7_Viterbi_Frameshift() and
 * p7_Viterbi_Frameshift_Trace() declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Forward declaration of dispatcher */
static int p7_Viterbi_Frameshift_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_OMX *ox, P7_OIVX *ov, float *opt_sc);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Viterbi_Frameshift)(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_OMX *ox, P7_OIVX *ov, float *opt_sc) = p7_Viterbi_Frameshift_Dispatcher;

/* Function:  p7_Viterbi_Frameshift_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for frameshift-aware Viterbi.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslERANGE> as documented in the _sse implementation.
 */
static int
p7_Viterbi_Frameshift_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_OMX *ox, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_Viterbi_Frameshift = p7_Viterbi_Frameshift_avx512; return p7_Viterbi_Frameshift_avx512(dsq, L, om_fs, ox, ov, opt_sc); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Viterbi_Frameshift = p7_Viterbi_Frameshift_avx;    return p7_Viterbi_Frameshift_avx(dsq, L, om_fs, ox, ov, opt_sc); }
#endif
#ifdef eslENABLE_SSE
  p7_Viterbi_Frameshift = p7_Viterbi_Frameshift_sse;
  return p7_Viterbi_Frameshift_sse(dsq, L, om_fs, ox, ov, opt_sc);
#else
  p7_Die("p7_Viterbi_Frameshift: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_Viterbi_Frameshift_Trace_Dispatcher(const ESL_DSQ *dsq, int L,
                                                   const P7_FS_OPROFILE *om_fs, const P7_OMX *ox,
                                                   P7_TRACE *tr);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Viterbi_Frameshift_Trace)(const ESL_DSQ *dsq, int L,
                                   const P7_FS_OPROFILE *om_fs, const P7_OMX *ox,
                                   P7_TRACE *tr) = p7_Viterbi_Frameshift_Trace_Dispatcher;

/* Function:  p7_Viterbi_Frameshift_Trace_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for frameshift-aware Viterbi traceback.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslFAIL> if an impossible state is reached during traceback.
 */
static int
p7_Viterbi_Frameshift_Trace_Dispatcher(const ESL_DSQ *dsq, int L,
                                        const P7_FS_OPROFILE *om_fs, const P7_OMX *ox,
                                        P7_TRACE *tr)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_Viterbi_Frameshift_Trace = p7_Viterbi_Frameshift_Trace_avx512; return p7_Viterbi_Frameshift_Trace_avx512(dsq, L, om_fs, ox, tr); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Viterbi_Frameshift_Trace = p7_Viterbi_Frameshift_Trace_avx;    return p7_Viterbi_Frameshift_Trace_avx(dsq, L, om_fs, ox, tr); }
#endif
#ifdef eslENABLE_SSE
  p7_Viterbi_Frameshift_Trace = p7_Viterbi_Frameshift_Trace_sse;
  return p7_Viterbi_Frameshift_Trace_sse(dsq, L, om_fs, ox, tr);
#else
  p7_Die("p7_Viterbi_Frameshift_Trace: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef p7VITERBI_FS_TESTDRIVE
#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"

static void
utest_viterbi_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
                 ESL_GENCODE *gcode, P7_BG *bgAA,
                 P7_CODONTABLE *codon_table, int M, int N)
{
  char           *msg    = "utest_viterbi_fs failed";
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_fs5 = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE *om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq    = NULL;
  P7_TRACE       *tr     = p7_trace_fs_Create();
  P7_TRACE       *trg    = p7_trace_fs_Create();
  P7_IVX         *iv5    = p7_ivx_Create(M, p7P_5CODONS);
  P7_OIVX        *ov5    = p7_oivx_Create(M, p7P_5CODONS);
  P7_GMX         *gx     = p7_gmx_Create(M, M, M, p7G_NSCELLS);
  P7_OMX         *ox     = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS);
  char            errbuf[eslERRBUFSIZE];
  float           gsc, osc;
  int             curr_L, i, j;

  if (!gm || !gm_fs5 || !om_fs5 || !sq || !tr || !trg || !iv5 || !ov5 || !gx || !ox) esl_fatal(msg);

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig(hmm, bgAA, gm, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, M, p7_LOCAL);
  p7_fs_oprofile_Convert_Log(gm_fs5, om_fs5);

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

      p7_fs_oprofile_ReconfigLength_Log(om_fs5, sq->n);
      p7_fs_ReconfigLength(gm_fs5, sq->n);

      p7_gmx_GrowTo(gx, M, curr_L, curr_L);
      p7_ivx_GrowTo(iv5, M, p7P_5CODONS);
      p7_oivx_GrowTo(ov5, M, p7P_5CODONS);
      p7_omx_GrowTo_dpf(ox, M, curr_L, curr_L);

      if (p7_GViterbi_Frameshift(dsq, curr_L, gm_fs5, gx, iv5, &gsc) != eslOK) esl_fatal(msg);

      { int s = p7_Viterbi_Frameshift(dsq, curr_L, om_fs5, ox, ov5, &osc);
        if (s == eslERANGE) continue;
        if (s != eslOK)     esl_fatal(msg); }

      if (fabs(gsc - osc) > 0.001)
        esl_fatal("%s: generic vit %.4f != SSE vit %.4f", msg, gsc, osc);

      if (p7_GVTrace_Frameshift(dsq, curr_L, gm_fs5, gx, trg)           != eslOK) esl_fatal(msg);
      if (p7_Viterbi_Frameshift_Trace(dsq, curr_L, om_fs5, ox, tr)       != eslOK) esl_fatal(msg);

      if (p7_trace_fs_Validate(tr,  abcDNA, dsq, errbuf) != eslOK)
        esl_fatal("%s: SSE trace invalid: %s",     msg, errbuf);
      if (p7_trace_fs_Validate(trg, abcDNA, dsq, errbuf) != eslOK)
        esl_fatal("%s: generic trace invalid: %s", msg, errbuf);

      p7_trace_Reuse(tr);
      p7_trace_Reuse(trg);
    }

  free(dsq);
  p7_trace_fs_Destroy(tr);
  p7_trace_fs_Destroy(trg);
  p7_omx_Destroy(ox);
  p7_ivx_Destroy(iv5);
  p7_oivx_Destroy(ov5);
  p7_gmx_Destroy(gx);
  p7_hmm_Destroy(hmm);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_fs5);
  p7_fs_oprofile_Destroy(om_fs5);
}

#if defined(eslENABLE_SSE) && defined(eslENABLE_AVX)
/* utest_sse_vs_avx_viterbi_fs()
 *
 * Run Viterbi_Frameshift_sse and Viterbi_Frameshift_avx on the same sequences
 * and profile, comparing scores.  Tests the 5-codon path only (most common).
 * Tolerance is 0.01 nats to accommodate float rounding from different stripe widths.
 * Skipped silently if AVX is not available at runtime.
 */
static void
utest_sse_vs_avx_viterbi_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
                             ESL_GENCODE *gcode, P7_BG *bgAA,
                             P7_CODONTABLE *codon_table, int M, int N)
{
  char           *msg       = "utest_sse_vs_avx_viterbi_fs failed";
  P7_HMM         *hmm       = NULL;
  P7_PROFILE     *gm        = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_fs5    = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE *om_fs_sse = NULL;
  P7_FS_OPROFILE *om_fs_avx = NULL;
  ESL_SQ         *sq        = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq       = NULL;
  P7_OIVX        *ov_sse    = NULL;
  P7_OIVX        *ov_avx    = NULL;
  P7_OMX         *ox_sse    = NULL;
  P7_OMX         *ox_avx    = NULL;
  float           sc_sse, sc_avx;
  int             curr_L, i, j;
  int             n         = N;

  if (!esl_cpu_has_avx()) goto cleanup;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig(hmm, bgAA, gm, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, M, p7_LOCAL);

  om_fs_sse = p7_fs_oprofile_Create_sse(M, abcAA, 5);
  om_fs_avx = p7_fs_oprofile_Create_avx(M, abcAA, 5);
  p7_fs_oprofile_Convert_Log_sse(gm_fs5, om_fs_sse);
  p7_fs_oprofile_Convert_Log_avx(gm_fs5, om_fs_avx);

  ov_sse = p7_oivx_Create_sse(M, p7P_5CODONS);
  ov_avx = p7_oivx_Create_avx(M, p7P_5CODONS);
  ox_sse = p7_omx_Create_dpf_sse(M, M, M, p7X_NSCELLS);
  ox_avx = p7_omx_Create_dpf_avx(M, M, M, p7X_NSCELLS);

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

      p7_fs_oprofile_ReconfigLength_Log_sse(om_fs_sse, sq->n);
      p7_fs_oprofile_ReconfigLength_Log_avx(om_fs_avx, sq->n);
      p7_omx_GrowTo_dpf_sse(ox_sse, M, curr_L, curr_L);
      p7_omx_GrowTo_dpf_avx(ox_avx, M, curr_L, curr_L);
      p7_oivx_GrowTo_sse(ov_sse, M, p7P_5CODONS);
      p7_oivx_GrowTo_avx(ov_avx, M, p7P_5CODONS);

      { int s_sse = p7_Viterbi_Frameshift_sse(dsq, curr_L, om_fs_sse, ox_sse, ov_sse, &sc_sse);
        int s_avx = p7_Viterbi_Frameshift_avx(dsq, curr_L, om_fs_avx, ox_avx, ov_avx, &sc_avx);
        if (s_sse == eslERANGE || s_avx == eslERANGE) continue;
        if (s_sse != eslOK || s_avx != eslOK) esl_fatal(msg); }

      if (fabs(sc_sse - sc_avx) > 0.01)
        esl_fatal("%s: sse=%.4f avx=%.4f", msg, sc_sse, sc_avx);
    }

 cleanup:
  if (dsq)       free(dsq);
  if (ox_sse)    p7_omx_Destroy_sse(ox_sse);
  if (ox_avx)    p7_omx_Destroy_avx(ox_avx);
  if (ov_sse)    p7_oivx_Destroy_sse(ov_sse);
  if (ov_avx)    p7_oivx_Destroy_avx(ov_avx);
  if (om_fs_sse) p7_fs_oprofile_Destroy_sse(om_fs_sse);
  if (om_fs_avx) p7_fs_oprofile_Destroy_avx(om_fs_avx);
  if (hmm)       p7_hmm_Destroy(hmm);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_fs5);
}
#endif /* eslENABLE_SSE && eslENABLE_AVX */

#endif /*p7VITERBI_FS_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/



/*****************************************************************
 * 4. Test driver.
 *****************************************************************/
#ifdef p7VITERBI_FS_TESTDRIVE
/*
  gcc -g -Wall -msse2 -mavx2 -std=gnu99 -o viterbi_fs_utest -I.. -L.. -I../../easel -L../../easel -Dp7VITERBI_FS_TESTDRIVE viterbi_fs.c -lhmmer -leasel -lm
  ./viterbi_fs_utest
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
  { "-M",  eslARG_INT,     "72", NULL, NULL, NULL, NULL, NULL, "size of random models to sample",                0 },
  { "-N",  eslARG_INT,     "20", NULL, NULL, NULL, NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for p7_Viterbi_Frameshift() and p7_Viterbi_Frameshift_Trace()";

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

  utest_viterbi_fs(r, abcAA, abcDNA, gcode, bgAA, ct, M, N);
#if defined(eslENABLE_SSE) && defined(eslENABLE_AVX)
  utest_sse_vs_avx_viterbi_fs(r, abcAA, abcDNA, gcode, bgAA, ct, M, N);
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
#endif /*p7VITERBI_FS_TESTDRIVE*/
/*----------------- end, test driver ----------------------------*/
