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


/* Forward declaration of dispatcher */
static int p7_Decoding_Frameshift_Dispatcher(const P7_FS_OPROFILE *om_fs, P7_OMX *fwd, const P7_OMX *bck);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Decoding_Frameshift)(const P7_FS_OPROFILE *om_fs, P7_OMX *fwd, const P7_OMX *bck) = p7_Decoding_Frameshift_Dispatcher;

/* Function:  p7_Decoding_Frameshift_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for frameshift posterior decoding.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if numeric range is exceeded.
 * Throws:    <eslEMEM> on allocation failure.
 */
static int
p7_Decoding_Frameshift_Dispatcher(const P7_FS_OPROFILE *om_fs, P7_OMX *fwd, const P7_OMX *bck)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_Decoding_Frameshift = p7_Decoding_Frameshift_avx512; return p7_Decoding_Frameshift_avx512(om_fs, fwd, bck); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Decoding_Frameshift = p7_Decoding_Frameshift_avx;    return p7_Decoding_Frameshift_avx(om_fs, fwd, bck); }
#endif
#ifdef eslENABLE_SSE
  p7_Decoding_Frameshift = p7_Decoding_Frameshift_sse;
  return p7_Decoding_Frameshift_sse(om_fs, fwd, bck);
#else
  p7_Die("p7_Decoding_Frameshift: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_DomainDecoding_Frameshift_Dispatcher(const P7_FS_OPROFILE *om_fs, const P7_OMX *oxf, const P7_OMX *oxb,
                                                    P7_DOMAINDEF *ddef);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_DomainDecoding_Frameshift)(const P7_FS_OPROFILE *om_fs, const P7_OMX *oxf, const P7_OMX *oxb,
                                    P7_DOMAINDEF *ddef) = p7_DomainDecoding_Frameshift_Dispatcher;

/* Function:  p7_DomainDecoding_Frameshift_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for frameshift domain decoding.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
static int
p7_DomainDecoding_Frameshift_Dispatcher(const P7_FS_OPROFILE *om_fs, const P7_OMX *oxf, const P7_OMX *oxb,
                                         P7_DOMAINDEF *ddef)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_DomainDecoding_Frameshift = p7_DomainDecoding_Frameshift_avx512; return p7_DomainDecoding_Frameshift_avx512(om_fs, oxf, oxb, ddef); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_DomainDecoding_Frameshift = p7_DomainDecoding_Frameshift_avx;    return p7_DomainDecoding_Frameshift_avx(om_fs, oxf, oxb, ddef); }
#endif
#ifdef eslENABLE_SSE
  p7_DomainDecoding_Frameshift = p7_DomainDecoding_Frameshift_sse;
  return p7_DomainDecoding_Frameshift_sse(om_fs, oxf, oxb, ddef);
#else
  p7_Die("p7_DomainDecoding_Frameshift: no SIMD implementation available");
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
#ifdef p7DECODING_FS_BENCHMARK
/*
   gcc -g -O3 -mavx2 -std=gnu99 -o benchmark-decoding_fs -I.. -L.. -I../../easel -L../../easel -Dp7DECODING_FS_BENCHMARK decoding_fs.c -lhmmer -leasel -lm
   icc -O3 -static -o benchmark-decoding_fs -I.. -L.. -I../../easel -L../../easel -Dp7DECODING_FS_BENCHMARK decoding_fs.c -lhmmer -leasel -lm
   ./benchmark-decoding_fs <hmmfile>

   Because p7_Decoding_Frameshift() overwrites the forward matrix in-place with
   posterior probabilities, we cannot pre-fill fwd/bck once and loop over
   Decoding alone.  Instead we time a baseline loop (Forward+Backward only) and
   a benchmark loop (Forward+Backward+Decoding), then subtract to isolate the
   decoding cost.
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
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-L",        eslARG_INT,   "1200", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs (nucleotides)",       0 },
  { "-N",        eslARG_INT,   "5000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for posterior residue decoding, AVX frameshift version";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA   = NULL;
  ESL_ALPHABET   *abcDNA  = NULL;
  ESL_GENCODE    *gcode   = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bgAA    = NULL;
  P7_BG          *bgDNA   = NULL;
  P7_FS_PROFILE  *gm_fs5  = NULL;
  P7_FS_OPROFILE *om_fs5  = NULL;
  P7_OMX         *fwd     = NULL;
  P7_OMX         *bck     = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           fsc, bsc;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)          != eslOK) p7_Fail("Failed to read HMM");

  abcDNA = esl_alphabet_Create(eslDNA);
  gcode  = esl_gencode_Create(abcDNA, abcAA);
  bgAA   = p7_bg_Create(abcAA);                 p7_bg_SetLength(bgAA, L/3);
  bgDNA  = p7_bg_Create(abcDNA);                p7_bg_SetLength(bgDNA, L);

  gm_fs5 = p7_profile_fs_Create(hmm->M, abcAA, p7P_5CODONS);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_UNILOCAL);

  om_fs5 = p7_fs_oprofile_Create(hmm->M, abcAA, p7P_5CODONS);
  p7_fs_oprofile_Convert(gm_fs5, om_fs5);
  p7_fs_oprofile_ReconfigLength(om_fs5, L/3);

  fwd = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS_FS);
  bck = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS);
  P7_OIVX *ov5 = p7_oivx_Create(hmm->M, p7P_5CODONS);

  /* Baseline: time Forward+Backward alone */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
      p7_Forward_Frameshift (dsq, L, om_fs5, fwd,      ov5, &fsc);
      p7_Backward_Frameshift(dsq, L, om_fs5, fwd, bck, ov5, &bsc);
    }
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark: Forward+Backward+Decoding */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
      p7_Forward_Frameshift (dsq, L, om_fs5, fwd,      ov5, &fsc);
      p7_Backward_Frameshift(dsq, L, om_fs5, fwd, bck, ov5, &bsc);
      p7_Decoding_Frameshift(om_fs5, fwd, bck);  /* overwrites fwd in-place */
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;

  Mcs = (double) N * (double) L * (double) hmm->M * 1e-6 / bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   hmm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_omx_Destroy(bck);
  p7_omx_Destroy(fwd);
  p7_oivx_Destroy(ov5);
  p7_fs_oprofile_Destroy(om_fs5);
  p7_profile_fs_Destroy(gm_fs5);
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
#endif /*p7DECODING_FS_BENCHMARK*/
/*------------------ end, benchmark driver ----------------------*/


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7DECODING_FS_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_decoding_fs()
 * Compare p7_Decoding_Frameshift() against the generic
 * p7_GDecoding_Frameshift() on <N> random DNA sequences of length <L>,
 * as utest_decoding() does for the standard model. Each sequence is
 * decoded with the generic and the vector full-matrix Forward/Backward,
 * and the two posterior matrices are compared within the relative
 * <tolerance>.
 */
static void
utest_decoding_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA, int M, int L, int N, float tolerance)
{
  char           *msg    = "decoding fs unit test failed";
  P7_HMM         *hmm    = NULL;
  P7_FS_PROFILE  *gm_fs5 = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE *om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
  ESL_DSQ        *dsq    = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_GMX         *gxf    = p7_gmx_Create(M, L, L, p7G_NSCELLS_FS);
  P7_GMX         *gxb    = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX         *gxpp   = p7_gmx_Create(M, L, L, p7G_NSCELLS_FS);
  P7_IVX         *iv5    = p7_ivx_Create(M, p7P_5CODONS);
  P7_OIVX        *ov5    = p7_oivx_Create(M, p7P_5CODONS);
  P7_OMX         *fwd    = p7_omx_Create_dpf(M, L, L, p7X_NSCELLS_FS);
  P7_OMX         *bck    = p7_omx_Create_dpf(M, L, L, p7X_NSCELLS);
  float           fsc, bsc;

  if (p7_hmm_Sample(r, M, abcAA, &hmm)                                  != eslOK) esl_fatal(msg);
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_LOCAL)      != eslOK) esl_fatal(msg);
  if (p7_fs_oprofile_Convert(gm_fs5, om_fs5)                            != eslOK) esl_fatal(msg);
  if (p7_fs_oprofile_ReconfigLength(om_fs5, L/3)                        != eslOK) esl_fatal(msg);

  while (N--)
    {
      if (esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq)                 != eslOK) esl_fatal(msg);

      /* generic: gxf becomes the posterior matrix */
      if (p7_GForward_Frameshift (dsq, L, gm_fs5, gxf, iv5, &fsc)       != eslOK) esl_fatal(msg);
      if (p7_GBackward_Frameshift(dsq, L, gm_fs5, gxb, iv5, &bsc)       != eslOK) esl_fatal(msg);
      if (p7_GDecoding_Frameshift(gm_fs5, gxf, gxb)                     != eslOK) esl_fatal(msg);

      /* vector: fwd becomes the posterior matrix */
      if (p7_Forward_Frameshift (dsq, L, om_fs5, fwd, ov5, &fsc)        != eslOK) esl_fatal(msg);
      if (p7_Backward_Frameshift(dsq, L, om_fs5, fwd, bck, ov5, &bsc)   != eslOK) esl_fatal(msg);
      if (p7_Decoding_Frameshift(om_fs5, fwd, bck)                      != eslOK) esl_fatal(msg);
      if (p7_omx_FDeconvert(fwd, gxpp)                                  != eslOK) esl_fatal(msg);

      if (p7_gmx_Compare(gxpp, gxf, tolerance)                          != eslOK) esl_fatal(msg);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_gmx_Destroy(gxf);
  p7_gmx_Destroy(gxb);
  p7_gmx_Destroy(gxpp);
  p7_ivx_Destroy(iv5);
  p7_oivx_Destroy(ov5);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bck);
  p7_profile_fs_Destroy(gm_fs5);
  p7_fs_oprofile_Destroy(om_fs5);
}

static void
utest_domdef(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA, int M, int L, int N)
{
  char           *msg    = "domain def fs unit test failed";
  P7_HMM         *hmm    = NULL;
  P7_FS_PROFILE  *gm_fs3 = p7_profile_fs_Create(M, abcAA, 3);
  P7_FS_OPROFILE *om_fs3 = p7_fs_oprofile_Create(M, abcAA, 3);
  ESL_DSQ        *dsq    = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX         *fwd    = p7_omx_Create(M, PARSER_ROWS_FWD, L);
  P7_OMX         *bwd    = p7_omx_Create(M, PARSER_ROWS_BWD, L);
  P7_GMX         *fgx    = p7_gmx_Create(M, PARSER_ROWS_FWD, L, p7G_NSCELLS);
  P7_GMX         *bgx    = p7_gmx_Create(M, PARSER_ROWS_BWD, L, p7G_NSCELLS);
  P7_IVX         *iv     = p7_ivx_Create(M, p7P_3CODONS);
  P7_OIVX        *ov3    = p7_oivx_Create(M, p7P_3CODONS);
  P7_DOMAINDEF    oddef, gddef;
  float           tolerance;
  int             i;

  if (p7_FLogsumError(-0.4, -0.5) > 0.0001) tolerance = 1.0;    /* weaker test against generic   */
  else                                      tolerance = 0.001;  /* stronger test: FLogsum() is in slow exact mode. */

  oddef.mocc = malloc(sizeof(float) * (L+1));
  oddef.btot = malloc(sizeof(float) * (L+1));
  oddef.etot = malloc(sizeof(float) * (L+1));
  gddef.mocc = malloc(sizeof(float) * (L+1));
  gddef.btot = malloc(sizeof(float) * (L+1));
  gddef.etot = malloc(sizeof(float) * (L+1));

  if (p7_hmm_Sample(r, M, abcAA, &hmm)                             != eslOK) esl_fatal(msg);
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs3, L/3, p7_LOCAL) != eslOK) esl_fatal(msg);
  if (p7_fs_oprofile_Convert(gm_fs3, om_fs3)                       != eslOK) esl_fatal(msg);
  if (p7_fs_oprofile_ReconfigLength(om_fs3, L/3)                   != eslOK) esl_fatal(msg);

  while (N--)
    {
      if (esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq)                      != eslOK) esl_fatal(msg);

      if (p7_ForwardParser_Frameshift_3Codons (dsq, L, om_fs3, fwd,      ov3, NULL) != eslOK) esl_fatal(msg);
      if (p7_BackwardParser_Frameshift_3Codons(dsq, L, om_fs3, fwd, bwd, ov3, NULL) != eslOK) esl_fatal(msg);
      if (p7_DomainDecoding_Frameshift(om_fs3, fwd, bwd, &oddef)             != eslOK) esl_fatal(msg);

      if (p7_GForwardParser_Frameshift_3Codons (dsq, L, gm_fs3, fgx, iv, NULL) != eslOK) esl_fatal(msg);
      if (p7_GBackwardParser_Frameshift_3Codons(dsq, L, gm_fs3, bgx, iv, NULL) != eslOK) esl_fatal(msg);
      if (p7_GDomainDecoding_Frameshift(gm_fs3, fgx, bgx, &gddef)            != eslOK) esl_fatal(msg);

      for (i = 0; i <= L; i++)
        {
          if (fabs(oddef.mocc[i] - gddef.mocc[i]) > tolerance) esl_fatal("domain def fs unit test failed at mocc");
          if (fabs(oddef.btot[i] - gddef.btot[i]) > tolerance) esl_fatal("domain def fs unit test failed at btot");
          if (fabs(oddef.etot[i] - gddef.etot[i]) > tolerance) esl_fatal("domain def fs unit test failed at etot");
        }
    }

  free(oddef.mocc);  free(oddef.btot);  free(oddef.etot);
  free(gddef.mocc);  free(gddef.btot);  free(gddef.etot);
  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bwd);
  p7_gmx_Destroy(fgx);
  p7_gmx_Destroy(bgx);
  p7_ivx_Destroy(iv);
  p7_oivx_Destroy(ov3);
  p7_profile_fs_Destroy(gm_fs3);
  p7_fs_oprofile_Destroy(om_fs3);
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
  /* name  type         default  env   range togs  reqs  incomp  help                                         docgrp */
  { "-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",                            0 },
  { "-s",  eslARG_INT,     "42",  NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",                  0 },
  { "-t",  eslARG_REAL,  "0.01",  NULL, NULL, NULL, NULL, NULL, "relative tolerance for the decoding test",       0 },
  { "-L",  eslARG_INT,    "120",  NULL, NULL, NULL, NULL, NULL, "length of random DNA seqs (nucleotides)",        0 },
  { "-M",  eslARG_INT,     "40",  NULL, NULL, NULL, NULL, NULL, "length of sampled test profile",                 0 },
  { "-N",  eslARG_INT,     "10",  NULL, NULL, NULL, NULL, NULL, "number of random test sequences",                0 },
  { 0,0,0,0,0,0,0,0,0,0},
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
  int             M      = esl_opt_GetInteger(go, "-M");
  int             N      = esl_opt_GetInteger(go, "-N");
  int             L      = esl_opt_GetInteger(go, "-L");
  float           tol    = esl_opt_GetReal   (go, "-t");

  if ((abcDNA = esl_alphabet_Create(eslDNA))      == NULL)  esl_fatal("failed to create alphabet");
  if ((bgDNA = p7_bg_Create(abcDNA))              == NULL)  esl_fatal("failed to create null model");
  if ((abcAA = esl_alphabet_Create(eslAMINO))     == NULL)  esl_fatal("failed to create alphabet");
  if ((bgAA = p7_bg_Create(abcAA))                == NULL)  esl_fatal("failed to create null model");
  if ((gcode  = esl_gencode_Create(abcDNA,abcAA)) == NULL)  esl_fatal("failed to create gencode");

  p7_FLogsumInit();
  impl_Init();

  utest_domdef      (r, abcAA, abcDNA, gcode, bgAA, bgDNA, M, L, N);
  utest_decoding_fs (r, abcAA, abcDNA, gcode, bgAA, bgDNA, M, L, N, tol);

  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  p7_bg_Destroy(bgDNA);
  p7_bg_Destroy(bgAA);
  esl_gencode_Destroy(gcode);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}

#endif /*p7DECODING_FS_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/


