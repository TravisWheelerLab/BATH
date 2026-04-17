/* Frameshift null2 dispatcher.
 * Provides the non-suffixed p7_Null2_fs_ByExpectation() declared in
 * impl_avx.h.  Delegates to the fastest available ISA implementation.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Forward declaration of dispatcher */
static int p7_Null2_fs_ByExpectation_Dispatcher(const P7_FS_OPROFILE *om_fs, P7_OMX *pp, float *null2);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Null2_fs_ByExpectation)(const P7_FS_OPROFILE *om_fs, P7_OMX *pp, float *null2) = p7_Null2_fs_ByExpectation_Dispatcher;

/* Function:  p7_Null2_fs_ByExpectation_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for frameshift null2 estimation (expectation).
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if codon_lengths is not 1, 3, or 5.
 */
static int
p7_Null2_fs_ByExpectation_Dispatcher(const P7_FS_OPROFILE *om_fs, P7_OMX *pp, float *null2)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_Null2_fs_ByExpectation = p7_Null2_fs_ByExpectation_avx512; return p7_Null2_fs_ByExpectation_avx512(om_fs, pp, null2); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Null2_fs_ByExpectation = p7_Null2_fs_ByExpectation_avx;    return p7_Null2_fs_ByExpectation_avx(om_fs, pp, null2); }
#endif
#ifdef eslENABLE_SSE
  p7_Null2_fs_ByExpectation = p7_Null2_fs_ByExpectation_sse;
  return p7_Null2_fs_ByExpectation_sse(om_fs, pp, null2);
#else
  p7_Die("p7_Null2_fs_ByExpectation: no SIMD implementation available");
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7NULL2_FS_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"

/* utest_null2_fs_expectation()
 * Compare p7_Null2_fs_ByExpectation() against the generic
 * p7_GNull2_fs_ByExpectation(). For each sampled sequence, run the full
 * pipeline on both paths and compare the resulting null2 odds vectors.
 *
 * The generic path works in log space (using p7_FLogsum); the SSE path
 * works in linear space.  Tolerance is set wide (0.2) when FLogsum is
 * in its fast approximate mode, tight (0.001) when in exact mode.
 */
static void
utest_null2_fs_expectation(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
                            ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA,
                            P7_CODONTABLE *codon_table, int M, int N)
{
  char           *msg    = "null2_fs_expectation unit test failed";
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_fs5 = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE *om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq    = NULL;
  P7_TRACE       *tr     = p7_trace_Create();
  P7_IVX         *iv5    = p7_ivx_Create(M, p7P_5CODONS);
  P7_OIVX        *ov5    = p7_oivx_Create(M, p7P_5CODONS);
  P7_GMX         *gxf    = p7_gmx_Create(M, M, M, p7G_NSCELLS_FS);
  P7_GMX         *gxb    = p7_gmx_Create(M, M, M, p7G_NSCELLS);
  P7_OMX         *fwd    = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS_FS);
  P7_OMX         *bck    = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS);
  float          *gn2    = malloc(sizeof(float) * abcAA->Kp);
  float          *on2    = malloc(sizeof(float) * abcAA->Kp);
  float           fsc, bsc;
  float           tolerance;
  int             i, j, curr_L;

  tolerance = (p7_FLogsumError(-0.4, -0.5) > 0.0001) ? 0.2f : 0.001f;

  if (!gn2 || !on2) esl_fatal(msg);

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig(hmm, bgAA, gm, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, M, p7_LOCAL);
  p7_fs_oprofile_Convert(gm_fs5, om_fs5);

  while (N--)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, tr);
      curr_L = sq->n * 3;

      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (curr_L + 2))) == NULL) esl_fatal("malloc failed");

      j = 1;
      for (i = 1; i <= sq->n; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }

      p7_fs_oprofile_ReconfigLength(om_fs5, sq->n);
      p7_fs_ReconfigLength(gm_fs5, sq->n);

      /* Generic reference: forward -> backward -> decoding -> null2_fs */
      p7_gmx_GrowTo(gxf, M, curr_L, curr_L);
      p7_gmx_GrowTo(gxb, M, curr_L, curr_L);
      p7_ivx_GrowTo(iv5, M, p7P_5CODONS);

      if (p7_GForward_Frameshift (dsq, curr_L, gm_fs5, gxf, iv5, &fsc) != eslOK) esl_fatal(msg);
      if (p7_GBackward_Frameshift(dsq, curr_L, gm_fs5, gxb, iv5, &bsc) != eslOK) esl_fatal(msg);
      if (p7_GDecoding_Frameshift(gm_fs5, gxf, gxb)                    != eslOK) esl_fatal(msg);
      if (p7_GNull2_fs_ByExpectation(gm_fs5, gxf, gn2)                  != eslOK) esl_fatal(msg);

      /* SSE path under test: forward -> backward -> decoding -> null2_fs */
      p7_omx_GrowTo_dpf(fwd, M, curr_L, curr_L);
      p7_omx_GrowTo_dpf(bck, M, curr_L, curr_L);

      if (p7_Forward_Frameshift (dsq, curr_L, om_fs5, fwd,           ov5, &fsc) == eslERANGE) continue;
      if (p7_Backward_Frameshift(dsq, curr_L, om_fs5, fwd, bck, ov5, &bsc) == eslERANGE) continue;
      if (p7_Decoding_Frameshift(om_fs5, fwd, bck)                    != eslOK) esl_fatal(msg);
      if (p7_Null2_fs_ByExpectation(om_fs5, fwd, on2)                 != eslOK) esl_fatal(msg);

      if (esl_vec_FCompare(gn2, on2, abcAA->Kp, tolerance) != eslOK) esl_fatal(msg);
    }

  free(dsq);
  free(gn2);
  free(on2);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_gmx_Destroy(gxf);
  p7_gmx_Destroy(gxb);
  p7_ivx_Destroy(iv5);
  p7_oivx_Destroy(ov5);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bck);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_fs5);
  p7_fs_oprofile_Destroy(om_fs5);
}
#endif /*p7NULL2_FS_TESTDRIVE*/
/*--------------------- end, unit tests -------------------------*/




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7NULL2_FS_TESTDRIVE
/*
  gcc -g -Wall -msse2 -std=gnu99 -o null2_fs_utest -I.. -L.. -I../../easel -L../../easel -Dp7NULL2_FS_TESTDRIVE null2_fs.c -lhmmer -leasel -lm
  ./null2_fs_utest
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
  { "-M",        eslARG_INT,     "72", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,     "10", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE frameshift null2 implementation";

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
  P7_BG          *bgDNA  = p7_bg_Create(abcDNA);
  int             M      = esl_opt_GetInteger(go, "-M");
  int             N      = esl_opt_GetInteger(go, "-N");

  p7_FLogsumInit();
  impl_Init();

  utest_null2_fs_expectation(r, abcAA, abcDNA, gcode, bgAA, bgDNA, ct, M, N);

  p7_bg_Destroy(bgAA);
  p7_bg_Destroy(bgDNA);
  p7_codontable_Destroy(ct);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcAA);
  esl_alphabet_Destroy(abcDNA);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7NULL2_FS_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/






