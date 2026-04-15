/* Frameshift posterior decoding dispatcher.
 * Provides the non-suffixed p7_Decoding_Frameshift() and
 * p7_DomainDecoding_Frameshift() declared in impl_avx.h.
 * Delegates to the fastest available ISA implementation at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Decoding_Frameshift()
 *
 * Purpose:   Dispatch frameshift posterior decoding (in-place overwrite of fwd)
 *            to the fastest available ISA path.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if numeric range is exceeded.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_Decoding_Frameshift(const P7_FS_OPROFILE *om_fs, P7_OMX *fwd, const P7_OMX *bck)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Decoding_Frameshift_avx512(om_fs, fwd, bck);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Decoding_Frameshift_avx(om_fs, fwd, bck);
#endif
#ifdef eslENABLE_SSE
  return p7_Decoding_Frameshift_sse(om_fs, fwd, bck);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_DomainDecoding_Frameshift()
 *
 * Purpose:   Dispatch frameshift domain decoding to the fastest available ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_DomainDecoding_Frameshift(const P7_FS_OPROFILE *om_fs, const P7_OMX *oxf, const P7_OMX *oxb,
                              P7_DOMAINDEF *ddef)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_DomainDecoding_Frameshift_avx512(om_fs, oxf, oxb, ddef);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_DomainDecoding_Frameshift_avx(om_fs, oxf, oxb, ddef);
#endif
#ifdef eslENABLE_SSE
  return p7_DomainDecoding_Frameshift_sse(om_fs, oxf, oxb, ddef);
#else
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7DECODING_FS_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_decoding_fs()
 * Compare p7_Decoding_Frameshift() against the generic p7_Decoding_Frameshift().
 * Both are given the same sequence and profile.  We run the generic full-matrix
 * forward/backward to produce a generic PP matrix (gx1), and the SSE full-matrix
 * forward/backward to produce an SSE PP matrix (fwd).  We compare the special-state
 * posteriors (N, J, C per row) between the two.
 */
static void
utest_decoding_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA, P7_CODONTABLE *codon_table, int M, int N)
{
  int  i, j;
  int  curr_L;
  char           *msg    = "decoding fs unit test failed";
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_fs5 = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE *om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq    = NULL;
  P7_TRACE       *tr     = p7_trace_Create();
  P7_GMX         *gx1    = p7_gmx_Create(M, M, M, p7G_NSCELLS_FS);
  P7_GMX         *gx2    = p7_gmx_Create(M, M, M, p7G_NSCELLS);
  P7_IVX         *iv5    = p7_ivx_Create(M, p7P_5CODONS);
  P7_OIVX        *ov5    = p7_oivx_Create(M, p7P_5CODONS);
  P7_OMX         *fwd    = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS_FS);
  P7_OMX         *bck    = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS);
  P7_GMX         *gxpp   = p7_gmx_Create(M, M, M, p7G_NSCELLS_FS);
  float           fsc, bsc;
  float           tolerance;

  tolerance = (p7_FLogsumError(-0.4, -0.5) > 0.0001) ? 0.2f : 0.001f;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig(hmm, bgAA, gm, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, M, p7_LOCAL);
  p7_fs_oprofile_Convert(gm_fs5, om_fs5);
  p7_fs_oprofile_ReconfigLength(om_fs5, M);

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

      /* Generic reference: p7_Decoding_Frameshift overwrites gx1 with PP  */
      p7_gmx_GrowTo(gx1, M, curr_L, curr_L);
      p7_gmx_GrowTo(gx2, M, curr_L, curr_L);
      p7_gmx_GrowTo(gxpp, M, curr_L, curr_L);
      p7_ivx_GrowTo(iv5, M, p7P_5CODONS);

      p7_GForward_Frameshift (dsq, curr_L, gm_fs5, gx1, iv5, &fsc);
      p7_GBackward_Frameshift(dsq, curr_L, gm_fs5, gx2, iv5, &bsc);
      p7_GDecoding_Frameshift(gm_fs5, gx1, gx2);   /* gx1 is now the PP matrix */

      /*  simd under test: p7_Decoding_Frameshift overwrites fwd with PP  */
      p7_omx_GrowTo_dpf(fwd, M, curr_L, curr_L);
      p7_omx_GrowTo_dpf(bck, M, curr_L, curr_L);

      if (p7_Forward_Frameshift (dsq, curr_L, om_fs5, fwd, ov5, &fsc) == eslERANGE) continue;
      if (p7_Backward_Frameshift(dsq, curr_L, om_fs5, fwd, bck, ov5, &bsc) == eslERANGE) continue;
      p7_Decoding_Frameshift(om_fs5, fwd, bck);  /* fwd is now the PP matrix */

      /* --- Deconvert SSE PP to generic layout and compare --- */
      p7_omx_FDeconvert(fwd, gxpp);
      if (p7_gmx_Compare(gxpp, gx1, tolerance) != eslOK) esl_fatal(msg);
    }

  free(dsq);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_gmx_Destroy(gx1);
  p7_gmx_Destroy(gx2);
  p7_gmx_Destroy(gxpp);
  p7_ivx_Destroy(iv5);
  p7_oivx_Destroy(ov5);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bck);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_fs5);
  p7_fs_oprofile_Destroy(om_fs5);
}

static void
utest_domdef(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA, P7_CODONTABLE *codon_table, int M, int N)
{
  int  i,j;
  int  curr_L;
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_fs3 = p7_profile_fs_Create(M, abcAA, 3);
  P7_FS_OPROFILE *om_fs3 = p7_fs_oprofile_Create(M, abcAA, 3);
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq    = NULL;
  P7_TRACE       *tr     = p7_trace_Create();
  P7_OMX         *fwd    = p7_omx_Create(M, PARSER_ROWS_FWD, M);
  P7_OMX         *bwd    = p7_omx_Create(M, PARSER_ROWS_BWD, M);
  P7_GMX         *fgx    = p7_gmx_Create(M, PARSER_ROWS_FWD, M, p7G_NSCELLS);
  P7_GMX         *bgx    = p7_gmx_Create(M, PARSER_ROWS_BWD, M, p7G_NSCELLS);
  P7_IVX         *iv     = p7_ivx_Create(M, p7P_3CODONS);
  P7_OIVX        *ov3    = p7_oivx_Create(M, p7P_3CODONS);
  P7_DOMAINDEF   *oddef  = NULL; 
  P7_DOMAINDEF   *gddef  = NULL;
  float tolerance;
  int status;

  float fsc3, bsc3;
  float generic_fsc3, generic_bsc3;
  
  if (p7_FLogsumError(-0.4, -0.5) > 0.0001) tolerance = 1.0;  /* weaker test against generic   */
  else tolerance = 0.001;   /* stronger test: FLogsum() is in slow exact mode. */  

  ESL_ALLOC(oddef, sizeof(P7_DOMAINDEF));
  ESL_ALLOC(gddef, sizeof(P7_DOMAINDEF));

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig(hmm, bgAA, gm, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs3, M, p7_LOCAL);
  p7_fs_oprofile_Convert(gm_fs3, om_fs3);
  p7_fs_oprofile_ReconfigLength(om_fs3, M);

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

      p7_omx_GrowTo(fwd, M, PARSER_ROWS_FWD, curr_L);
      p7_omx_GrowTo(bwd, M, PARSER_ROWS_BWD, curr_L);

      if (p7_ForwardParser_Frameshift_3Codons(dsq, curr_L, om_fs3, fwd, ov3, &fsc3)     == eslERANGE) continue;
      if (p7_BackwardParser_Frameshift_3Codons(dsq, curr_L, om_fs3, fwd, bwd, ov3, &generic_fsc3) == eslERANGE) continue;

      oddef->mocc = oddef->btot = oddef->etot = NULL;
      ESL_ALLOC(oddef->mocc, sizeof(float) * (curr_L+1));
      ESL_ALLOC(oddef->btot, sizeof(float) * (curr_L+1));
      ESL_ALLOC(oddef->etot, sizeof(float) * (curr_L+1));
      
      p7_DomainDecoding_Frameshift(om_fs3, fwd, bwd, oddef);

      p7_gmx_GrowTo(fgx, M, PARSER_ROWS_FWD, curr_L);
      p7_gmx_GrowTo(bgx, M, PARSER_ROWS_BWD, curr_L);

      p7_GForwardParser_Frameshift_3Codons(dsq, curr_L, gm_fs3, fgx, iv, &bsc3);
      p7_GBackwardParser_Frameshift_3Codons(dsq, curr_L, gm_fs3, bgx, iv, &generic_bsc3);

      gddef->mocc = gddef->btot = gddef->etot = NULL;
      ESL_ALLOC(gddef->mocc, sizeof(float) * (curr_L+1));
      ESL_ALLOC(gddef->btot, sizeof(float) * (curr_L+1));
      ESL_ALLOC(gddef->etot, sizeof(float) * (curr_L+1));
      p7_GDomainDecoding_Frameshift(gm_fs3, fgx, bgx, gddef);
     
      for(i = 0; i <= curr_L; i++) {
        if(fabs(oddef->mocc[i] - gddef->mocc[i]) > tolerance) esl_fatal("domain def fs unit test failed at mocc");
        if(fabs(oddef->btot[i] - gddef->btot[i]) > tolerance) esl_fatal("domain def fs unit test failed at btot");
        if(fabs(oddef->etot[i] - gddef->etot[i]) > tolerance) esl_fatal("domain def fs unit test failed at etot");
      }

      free(oddef->mocc);
      free(oddef->btot);
      free(oddef->etot);
      free(gddef->mocc);
      free(gddef->btot);
      free(gddef->etot);

    }

  free(dsq);
  free(oddef);
  free(gddef);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bwd);
  p7_gmx_Destroy(fgx);
  p7_gmx_Destroy(bgx);
  p7_ivx_Destroy(iv);
  p7_oivx_Destroy(ov3);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_fs3);
  p7_fs_oprofile_Destroy(om_fs3);

ERROR:
  return;
}
#endif /*p7DECODING_FS_TESTDRIVE*/


/*--------------------- end, unit tests -------------------------*/




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7DECODING_FS_TESTDRIVE
/* 
  gcc -g -Wall -msse2 -std=gnu99 -o decoding_fs_utest -I.. -L.. -I../../easel -L../../easel -Dp7DECODING_FS_TESTDRIVE decoding_fs.c -lhmmer -leasel -lm
 ./decoding_fs_utest
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
static char banner[] = "test driver for SSE frameshift domain definition ";

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
  if ((bgDNA = p7_bg_Create(abcDNA))              == NULL)  esl_fatal("failed to create null model");
  if ((abcAA = esl_alphabet_Create(eslAMINO))     == NULL)  esl_fatal("failed to create alphabet");
  if ((bgAA = p7_bg_Create(abcAA))                == NULL)  esl_fatal("failed to create null model");
  if ((gcode  = esl_gencode_Create(abcDNA,abcAA)) == NULL)  esl_fatal("failed to create gencode");
  if ((ct     = p7_codontable_Create(gcode))      == NULL)  esl_fatal("failed to create codon table");

  impl_Init();
  p7_FLogsumInit();

  utest_domdef      (r, abcAA, abcDNA, gcode, bgAA, bgDNA, ct, M, N);
  utest_decoding_fs (r, abcAA, abcDNA, gcode, bgAA, bgDNA, ct, M, N);

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

#endif /*p7DECODING_FS_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/


