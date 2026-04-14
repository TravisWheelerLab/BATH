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


/* Function:  p7_ForwardParser_Frameshift_3Codons()
 *
 * Purpose:   Dispatch 3-codon frameshift Forward parser to the fastest ISA.
 *
 * Returns:   <eslOK> on success.  <*opt_sc> is the log Forward score in nats.
 * Throws:    <eslEINVAL> if allocation is too small.
 *            <eslEMEM> on reallocation failure.
 */
int
p7_ForwardParser_Frameshift_3Codons(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                    P7_OMX *ox, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_ForwardParser_Frameshift_3Codons_avx512(dsq, L, om_fs, ox, ov, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_ForwardParser_Frameshift_3Codons_sse(dsq, L, om_fs, ox, ov, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_ForwardParser_Frameshift_3Codons_sse(dsq, L, om_fs, ox, ov, opt_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_BackwardParser_Frameshift_3Codons()
 *
 * Purpose:   Dispatch 3-codon frameshift Backward parser to the fastest ISA.
 *
 * Returns:   <eslOK> on success.  <*opt_sc> is the log Backward score.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
int
p7_BackwardParser_Frameshift_3Codons(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                     const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_BackwardParser_Frameshift_3Codons_avx512(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_BackwardParser_Frameshift_3Codons_sse(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_BackwardParser_Frameshift_3Codons_sse(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_ForwardParser_Frameshift_5Codons()
 *
 * Purpose:   Dispatch 5-codon frameshift Forward parser to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
int
p7_ForwardParser_Frameshift_5Codons(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                    P7_OMX *ox, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_ForwardParser_Frameshift_5Codons_avx512(dsq, L, om_fs, ox, ov, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_ForwardParser_Frameshift_5Codons_sse(dsq, L, om_fs, ox, ov, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_ForwardParser_Frameshift_5Codons_sse(dsq, L, om_fs, ox, ov, opt_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_BackwardParser_Frameshift_5Codons()
 *
 * Purpose:   Dispatch 5-codon frameshift Backward parser to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
int
p7_BackwardParser_Frameshift_5Codons(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                     const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_BackwardParser_Frameshift_5Codons_avx512(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_BackwardParser_Frameshift_5Codons_sse(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_BackwardParser_Frameshift_5Codons_sse(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_Forward_Frameshift()
 *
 * Purpose:   Dispatch full-matrix 5-codon frameshift Forward to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
int
p7_Forward_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                      P7_OMX *ox, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Forward_Frameshift_avx512(dsq, L, om_fs, ox, ov, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Forward_Frameshift_sse(dsq, L, om_fs, ox, ov, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_Forward_Frameshift_sse(dsq, L, om_fs, ox, ov, opt_sc);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_Backward_Frameshift()
 *
 * Purpose:   Dispatch full-matrix 5-codon frameshift Backward to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL>, <eslEMEM> as above.
 */
int
p7_Backward_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                       const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Backward_Frameshift_avx512(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Backward_Frameshift_sse(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#endif
#ifdef eslENABLE_SSE
  return p7_Backward_Frameshift_sse(dsq, L, om_fs, fwd, bck, ov, opt_sc);
#else
  return eslENORESULT;
#endif
}
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

