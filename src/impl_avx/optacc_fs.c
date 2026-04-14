/* Frameshift optimal accuracy dispatcher.
 * Provides the non-suffixed p7_OptimalAccuracy_Frameshift() and
 * p7_OATrace_Frameshift() declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_OptimalAccuracy_Frameshift()
 *
 * Purpose:   Dispatch frameshift optimal accuracy DP fill to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_OptimalAccuracy_Frameshift(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                               P7_OMX *ox, float *ret_e)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_OptimalAccuracy_Frameshift_avx512(om_fs, pp, ox, ret_e);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_OptimalAccuracy_Frameshift_sse(om_fs, pp, ox, ret_e);
#endif
#ifdef eslENABLE_SSE
  return p7_OptimalAccuracy_Frameshift_sse(om_fs, pp, ox, ret_e);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_OATrace_Frameshift()
 *
 * Purpose:   Dispatch frameshift optimal accuracy traceback to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if trace is not empty.
 *            <eslEMEM> on allocation error.
 */
int
p7_OATrace_Frameshift(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                      const P7_OMX *ox, P7_TRACE *tr)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_OATrace_Frameshift_avx512(om_fs, pp, ox, tr);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_OATrace_Frameshift_sse(om_fs, pp, ox, tr);
#endif
#ifdef eslENABLE_SSE
  return p7_OATrace_Frameshift_sse(om_fs, pp, ox, tr);
#else
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7OPTACC_FS_TESTDRIVE
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_optacc_fs:
 *   Compare p7_OptimalAccuracy_Frameshift and p7_OATrace_Frameshift
 *   against their scalar counterparts.
 *
 *   Strategy: run the full SSE pipeline (Forward->Backward->Decoding->OA->Trace)
 *   and the scalar pipeline (via deconverted pp matrix).  Compare OA scores,
 *   OA DP matrices, and tracebacks for equality within tolerance.
 *
 *   Note: Like utest_optacc in optacc.c, this test can fail for innocent
 *   reasons when numerical ties in the OA matrix lead to different traceback
 *   choices between SSE and scalar implementations.
 */
static void
utest_optacc_fs(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA,
                ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode,
                P7_BG *bgAA, P7_BG *bgDNA, int M, int L, int N)
{
  char            *msg   = "optacc_fs SSE unit test failed";
  P7_HMM          *hmm   = NULL;
  P7_FS_PROFILE   *gm_fs = NULL;
  P7_FS_OPROFILE  *om_fs = NULL;
  P7_GMX          *gx_pp  = NULL;   /* scalar pp (8-cell FS) */
  P7_GMX          *gx_oa  = NULL;   /* scalar OA (3-cell)    */
  P7_GMX          *gx_oa2 = NULL;   /* SSE OA deconverted    */
  P7_OMX          *ox_fwd = NULL;   /* SSE forward/pp (8-cell FS) */
  P7_OMX          *ox_bck = NULL;   /* SSE backward (8-cell FS)   */
  P7_TRACE        *tr_sse = NULL;   /* SSE OA traceback           */
  P7_TRACE        *tr_ref = NULL;   /* scalar OA traceback        */
  ESL_DSQ         *dsq    = NULL;
  float            fsc, bsc, accscore_ref, accscore_sse;
  float            tol    = 0.01f;
  float            pptol  = 0.01f;
  int              seq_L;

  p7_FLogsumInit();

  if (p7_hmm_Sample(r, M, abcAA, &hmm)                              != eslOK) esl_fatal(msg);

  gm_fs = p7_profile_fs_Create(hmm->M, abcAA, p7P_5CODONS);
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs, L, p7_UNILOCAL) != eslOK) esl_fatal(msg);

  om_fs = p7_fs_oprofile_Create(hmm->M, abcAA, p7P_5CODONS);
  if (p7_fs_oprofile_Convert(gm_fs, om_fs)                          != eslOK) esl_fatal(msg);

  seq_L = L;
  dsq   = malloc(sizeof(ESL_DSQ) * (seq_L + 2));

  /* SSE: forward/backward use 8-cell FS layout; OA uses standard 3-cell */
  ox_fwd = p7_omx_Create_dpf(hmm->M, seq_L, seq_L, p7X_NSCELLS_FS);
  ox_bck = p7_omx_Create_dpf(hmm->M, seq_L, seq_L, p7X_NSCELLS);
  P7_OIVX *ov5 = p7_oivx_Create(hmm->M, p7P_5CODONS);

  /* Scalar: pp uses 8-cell FS layout; OA matrices use standard 3-cell */
  gx_pp  = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS_FS);
  gx_oa  = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS);
  gx_oa2 = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS);

  tr_sse = p7_trace_fs_CreateWithPP();
  tr_ref = p7_trace_fs_CreateWithPP();

  while (N--)
    {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, seq_L, dsq);
      
      /* SSE pipeline: Forward -> Backward -> Decoding (in-place into ox_fwd) -> OA -> Trace */
      if (p7_Forward_Frameshift (dsq, seq_L, om_fs,         ox_fwd,      ov5, &fsc) == eslERANGE) continue;
      if (p7_Backward_Frameshift(dsq, seq_L, om_fs, ox_fwd, ox_bck, ov5, &bsc) == eslERANGE) continue;
      if (esl_FCompare_old(fsc, bsc, tol)                                      != eslOK) esl_fatal(msg);
      if (p7_Decoding_Frameshift(om_fs, ox_fwd, ox_bck)                    != eslOK) esl_fatal(msg);
      /* ox_fwd now holds the SSE posterior probability matrix */
      if (p7_OptimalAccuracy_Frameshift(om_fs, ox_fwd, ox_bck, &accscore_sse) != eslOK) esl_fatal(msg);
      if (p7_OATrace_Frameshift(om_fs, ox_fwd, ox_bck, tr_sse)              != eslOK) esl_fatal(msg);

      /* Reference: convert SSE pp to scalar GMX, run scalar OA and trace */
      if (p7_omx_FDeconvert(ox_fwd, gx_pp)                                    != eslOK) esl_fatal(msg);
      if (p7_GOptimalAccuracy_Frameshift(gm_fs, gx_pp, gx_oa, &accscore_ref)   != eslOK) esl_fatal(msg);
      if (p7_GOATrace_Frameshift(gm_fs, gx_pp, gx_oa, tr_ref)                  != eslOK) esl_fatal(msg);

      /* Convert SSE OA to scalar GMX for matrix-level comparison */
      if (p7_omx_FDeconvert(ox_bck, gx_oa2)                                    != eslOK) esl_fatal(msg);
      
      if (esl_FCompare_old(accscore_sse, accscore_ref, tol) != eslOK) esl_fatal(msg);
      if (p7_gmx_Compare(gx_oa, gx_oa2, tol)               != eslOK) esl_fatal(msg);
      if (p7_trace_Compare(tr_sse, tr_ref, pptol)           != eslOK) esl_fatal(msg);
      
      if (esl_opt_GetBoolean(go, "--traces"))
        {
          p7_trace_fs_Dump(stdout, tr_sse, gm_fs, dsq, abcDNA);
          p7_trace_fs_Dump(stdout, tr_ref, gm_fs, dsq, abcDNA);
        }

      p7_trace_Reuse(tr_sse);
      p7_trace_Reuse(tr_ref);
    }

  free(dsq);
  p7_trace_fs_Destroy(tr_sse);
  p7_trace_fs_Destroy(tr_ref);
  p7_omx_Destroy(ox_fwd);
  p7_omx_Destroy(ox_bck);
  p7_oivx_Destroy(ov5);
  p7_gmx_Destroy(gx_pp);
  p7_gmx_Destroy(gx_oa);
  p7_gmx_Destroy(gx_oa2);
  p7_fs_oprofile_Destroy(om_fs);
  p7_profile_fs_Destroy(gm_fs);
  p7_hmm_Destroy(hmm);
}
#endif /*p7OPTACC_FS_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/




/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7OPTACC_FS_TESTDRIVE
/*
   gcc -g -Wall -msse2 -std=gnu99 -o optacc_fs_utest \
       -I.. -L.. -I../../easel -L../../easel \
       -Dp7OPTACC_FS_TESTDRIVE optacc_fs.c -lhmmer -leasel -lm
   ./optacc_fs_utest
*/
#include "p7_config.h"
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "hmmer.h"
#include "impl_avx.h"

static ESL_OPTIONS options[] = {
  { "-h",       eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help",              0 },
  { "-s",       eslARG_INT,    "42", NULL, NULL, NULL, NULL, NULL, "set random number seed",       0 },
  { "-L",       eslARG_INT,   "150", NULL, NULL, NULL, NULL, NULL, "length of random target seqs", 0 },
  { "-M",       eslARG_INT,    "45", NULL, NULL, NULL, NULL, NULL, "size of random models",        0 },
  { "-N",       eslARG_INT,    "20", NULL, NULL, NULL, NULL, NULL, "number of test sequences",     0 },
  { "--traces", eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "dump all tracebacks",          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE frameshift optimal accuracy";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA  = esl_alphabet_Create(eslAMINO);
  ESL_ALPHABET   *abcDNA = esl_alphabet_Create(eslDNA);
  ESL_GENCODE    *gcode  = esl_gencode_Create(abcDNA, abcAA);
  P7_BG          *bgAA   = p7_bg_Create(abcAA);
  P7_BG          *bgDNA  = p7_bg_Create(abcDNA);
  int             M      = esl_opt_GetInteger(go, "-M");
  int             L      = esl_opt_GetInteger(go, "-L");
  int             N      = esl_opt_GetInteger(go, "-N");

  impl_Init();

  utest_optacc_fs(go, r, abcAA, abcDNA, gcode, bgAA, bgDNA, M, L, N);

  p7_bg_Destroy(bgAA);
  p7_bg_Destroy(bgDNA);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcAA);
  esl_alphabet_Destroy(abcDNA);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7OPTACC_FS_TESTDRIVE*/
/*------------------ end, test driver ---------------------------*/




