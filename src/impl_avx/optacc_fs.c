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


/* Forward declaration of dispatcher */
static int p7_OptimalAccuracy_Frameshift_Dispatcher(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                                                    P7_OMX *ox, float *ret_e);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_OptimalAccuracy_Frameshift)(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                                     P7_OMX *ox, float *ret_e) = p7_OptimalAccuracy_Frameshift_Dispatcher;

/* Function:  p7_OptimalAccuracy_Frameshift_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for frameshift optimal accuracy DP fill.
 *
 * Returns:   <eslOK> on success.
 */
static int
p7_OptimalAccuracy_Frameshift_Dispatcher(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                                         P7_OMX *ox, float *ret_e)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_OptimalAccuracy_Frameshift = p7_OptimalAccuracy_Frameshift_avx512; return p7_OptimalAccuracy_Frameshift_avx512(om_fs, pp, ox, ret_e); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_OptimalAccuracy_Frameshift = p7_OptimalAccuracy_Frameshift_avx;    return p7_OptimalAccuracy_Frameshift_avx(om_fs, pp, ox, ret_e); }
#endif
#ifdef eslENABLE_SSE
  p7_OptimalAccuracy_Frameshift = p7_OptimalAccuracy_Frameshift_sse;
  return p7_OptimalAccuracy_Frameshift_sse(om_fs, pp, ox, ret_e);
#else
  p7_Die("p7_OptimalAccuracy_Frameshift: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_OATrace_Frameshift_Dispatcher(const P7_FS_OPROFILE *om_fs, const float *ppx,
                                             const P7_OMX *ox, P7_TRACE *tr);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_OATrace_Frameshift)(const P7_FS_OPROFILE *om_fs, const float *ppx,
                             const P7_OMX *ox, P7_TRACE *tr) = p7_OATrace_Frameshift_Dispatcher;

/* Function:  p7_OATrace_Frameshift_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for frameshift optimal accuracy traceback.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if trace is not empty.
 *            <eslEMEM> on allocation error.
 */
static int
p7_OATrace_Frameshift_Dispatcher(const P7_FS_OPROFILE *om_fs, const float *ppx,
                                  const P7_OMX *ox, P7_TRACE *tr)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_OATrace_Frameshift = p7_OATrace_Frameshift_avx512; return p7_OATrace_Frameshift_avx512(om_fs, ppx, ox, tr); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_OATrace_Frameshift = p7_OATrace_Frameshift_avx;    return p7_OATrace_Frameshift_avx(om_fs, ppx, ox, tr); }
#endif
#ifdef eslENABLE_SSE
  p7_OATrace_Frameshift = p7_OATrace_Frameshift_sse;
  return p7_OATrace_Frameshift_sse(om_fs, ppx, ox, tr);
#else
  p7_Die("p7_OATrace_Frameshift: no SIMD implementation available");
  return eslENORESULT;
#endif
}
/* Forward declarations of dispatchers */
static int p7_OptimalAccuracySweep_Frameshift_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, P7_OATMX *ox, float *opt_null2, float *opt_fsc, float *ret_oasc);
static int p7_OATraceSweep_Frameshift_Dispatcher(const P7_FS_OPROFILE *om_fs, const P7_OATMX *ox, P7_TRACE *tr);

/* Global function pointers, initially pointing at the dispatchers */
int (*p7_OptimalAccuracySweep_Frameshift)(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, P7_OATMX *ox, float *opt_null2, float *opt_fsc, float *ret_oasc) = p7_OptimalAccuracySweep_Frameshift_Dispatcher;
int (*p7_OATraceSweep_Frameshift)(const P7_FS_OPROFILE *om_fs, const P7_OATMX *ox, P7_TRACE *tr) = p7_OATraceSweep_Frameshift_Dispatcher;

/* Function:  p7_OptimalAccuracySweep_Frameshift_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for the frameshift OA sweep over a
 *            checkpointed Backward. See p7_OptimalAccuracySweep_Frameshift_avx().
 */
static int
p7_OptimalAccuracySweep_Frameshift_Dispatcher(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, P7_OATMX *ox, float *opt_null2, float *opt_fsc, float *ret_oasc)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_OptimalAccuracySweep_Frameshift = p7_OptimalAccuracySweep_Frameshift_avx; return p7_OptimalAccuracySweep_Frameshift_avx(dsq, L, om_fs, cx, ov, ox, opt_null2, opt_fsc, ret_oasc); }
#endif
#ifdef eslENABLE_SSE
  p7_OptimalAccuracySweep_Frameshift = p7_OptimalAccuracySweep_Frameshift_sse;
  return p7_OptimalAccuracySweep_Frameshift_sse(dsq, L, om_fs, cx, ov, ox, opt_null2, opt_fsc, ret_oasc);
#else
  p7_Die("p7_OptimalAccuracySweep_Frameshift: no SIMD implementation available");
  return eslENORESULT;
#endif
}

/* Function:  p7_OATraceSweep_Frameshift_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for the frameshift OA trace from what
 *            the OA sweep stored. See p7_OATraceSweep_Frameshift_avx().
 */
static int
p7_OATraceSweep_Frameshift_Dispatcher(const P7_FS_OPROFILE *om_fs, const P7_OATMX *ox, P7_TRACE *tr)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx()) { p7_OATraceSweep_Frameshift = p7_OATraceSweep_Frameshift_avx; return p7_OATraceSweep_Frameshift_avx(om_fs, ox, tr); }
#endif
#ifdef eslENABLE_SSE
  p7_OATraceSweep_Frameshift = p7_OATraceSweep_Frameshift_sse;
  return p7_OATraceSweep_Frameshift_sse(om_fs, ox, tr);
#else
  p7_Die("p7_OATraceSweep_Frameshift: no SIMD implementation available");
  return eslENORESULT;
#endif
}

/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/
#ifdef p7OPTACC_FS_BENCHMARK
/*
   gcc -g -O3 -mavx2 -std=gnu99 -o benchmark-optacc_fs -I.. -L.. -I../../easel -L../../easel -Dp7OPTACC_FS_BENCHMARK optacc_fs.c -lhmmer -leasel -lm
   icc -O3 -static -o benchmark-optacc_fs -I.. -L.. -I../../easel -L../../easel -Dp7OPTACC_FS_BENCHMARK optacc_fs.c -lhmmer -leasel -lm
   ./benchmark-optacc_fs <hmmfile>
   ./benchmark-optacc_fs --notrace <hmmfile>   (benchmark DP fill only)

   Fill pipeline (Forward->Backward->Decoding) is run once to produce a
   posterior probability matrix before the timed loop.  p7_Decoding_Frameshift()
   overwrites the forward matrix in-place with the pp matrix; because pp is
   read-only in p7_OptimalAccuracy_Frameshift() and p7_OATrace_Frameshift(),
   the timed loop measures only OA fill and trace cost.
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
  { "--notrace", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark the DP fill stage",                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for optimal accuracy alignment, AVX frameshift version";

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
  P7_OMX         *ox_fwd  = NULL;   /* forward matrix; overwritten in-place by Decoding -> pp */
  P7_OMX         *ox_bck  = NULL;   /* backward matrix                                  */
  P7_OMX         *ox_oa   = NULL;   /* OA matrix (8-cell FS layout)                     */
  P7_TRACE       *tr      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           fsc, bsc, accscore;
  double          Mcs;

  p7_FLogsumInit();

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

  /* ox_fwd: 8-cell FS layout (Forward input; becomes pp after Decoding)  */
  ox_fwd = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS_FS);
  ox_bck = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS);
  ox_oa  = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS_FS);
  P7_OIVX *ov5 = p7_oivx_Create(hmm->M, p7P_5CODONS);
  tr     = p7_trace_fs_CreateWithPP();

  /* Fill pipeline once to produce pp matrix in ox_fwd */
  esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
  p7_Forward_Frameshift_5Codons (dsq, L, om_fs5, ox_fwd, ov5, &fsc);
  p7_Backward_Frameshift_5Codons(dsq, L, om_fs5, ox_fwd, ox_bck, ov5, &bsc);
  p7_Decoding_Frameshift(om_fs5, ox_fwd, ox_bck);  /* ox_fwd now holds pp */

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      p7_OptimalAccuracy_Frameshift(om_fs5, ox_fwd, ox_oa, &accscore);

      if (! esl_opt_GetBoolean(go, "--notrace"))
        {
          p7_OATrace_Frameshift(om_fs5, ox_fwd->xmx, ox_oa, tr);
          p7_trace_Reuse(tr);
        }
    }
  esl_stopwatch_Stop(w);

  Mcs = (double) N * (double) L * (double) hmm->M * 1e-6 / w->user;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   hmm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_trace_fs_Destroy(tr);
  p7_omx_Destroy(ox_bck);
  p7_omx_Destroy(ox_oa);
  p7_omx_Destroy(ox_fwd);
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
#endif /*p7OPTACC_FS_BENCHMARK*/
/*---------------- end, benchmark driver ------------------------*/


/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7OPTACC_FS_TESTDRIVE
#include <string.h>

#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"

/* oa_fs_compare()
 * Compare two frameshift OA matrices (scalar, or SIMD deconverted): every
 * slot of rows 1..L, the OA slots (C1, C2, D) of row 0, which is all the
 * fill writes there, and the special states. With <tol> 0 the values must
 * be identical.
 */
static int
oa_fs_compare(const P7_GMX *gx1, const P7_GMX *gx2, float tol)
{
  int i, k, s, x;

  if (gx1->M != gx2->M || gx1->L != gx2->L) return eslFAIL;

  for (i = 0; i <= gx1->L; i++)
    {
      for (k = 1; k <= gx1->M; k++)
        for (s = 0; s < p7G_NSCELLS_FS; s++)
          if (i > 0 || s == p7G_D || s == p7G_M + p7G_C1 || s == p7G_M + p7G_C2)
            if (esl_FCompare_old(gx1->dp[i][k*p7G_NSCELLS_FS + s], gx2->dp[i][k*p7G_NSCELLS_FS + s], tol) != eslOK) return eslFAIL;
      for (x = 0; x < p7G_NXCELLS; x++)
        if (esl_FCompare_old(gx1->xmx[i*p7G_NXCELLS + x], gx2->xmx[i*p7G_NXCELLS + x], tol) != eslOK) return eslFAIL;
    }
  return eslOK;
}

/* utest_optacc_fs:
 *   Compare p7_OptimalAccuracy_Frameshift and p7_OATrace_Frameshift
 *   against their scalar counterparts.
 *
 *   Strategy: run the full SSE pipeline (Forward->Backward->Decoding->OA->Trace)
 *   and the scalar pipeline (via deconverted pp matrix).  Compare OA scores,
 *   OA DP matrices, and tracebacks for equality within tolerance.
 *
 *   With <use_sse> TRUE the test runs the SSE engines (the _sse
 *   functions) instead of the dispatched ones, as on a CPU without AVX.
 *
 *   Note: Like utest_optacc in optacc.c, this test can fail for innocent
 *   reasons when numerical ties in the OA matrix lead to different traceback
 *   choices between SSE and scalar implementations.
 */
static void
utest_optacc_fs(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA,
                ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode,
                P7_BG *bgAA, P7_BG *bgDNA, int M, int L, int N, int use_sse)
{
  char            *msg   = (use_sse ? "optacc_fs SSE unit test failed" : "optacc_fs unit test failed");
  P7_HMM          *hmm   = NULL;
  P7_FS_PROFILE   *gm_fs = NULL;
  P7_FS_OPROFILE  *om_fs = NULL;
  P7_GMX          *gx_pp  = NULL;   /* scalar pp (8-cell FS) */
  P7_GMX          *gx_oa  = NULL;   /* scalar OA (8-cell FS) */
  P7_GMX          *gx_oa2 = NULL;   /* SSE OA deconverted (8-cell) */
  P7_OMX          *ox_fwd = NULL;   /* SSE forward/pp (8-cell FS) */
  P7_OMX          *ox_bck = NULL;   /* SSE backward (3-cell)      */
  P7_OMX          *ox_oa  = NULL;   /* SIMD OA (8-cell FS)         */
  P7_TRACE        *tr_sse = NULL;   /* SSE OA traceback           */
  P7_TRACE        *tr_ref = NULL;   /* scalar OA traceback        */
  P7_TRACE        *tr_ip  = NULL;   /* scalar OA traceback, in place */
  float           *ppx    = NULL;   /* copy of the pp special states */
  ESL_DSQ         *dsq    = NULL;
  float            fsc, bsc, accscore_ref, accscore_sse, accscore_ip, accscore_simd_ip;
  float            tol    = 0.01f;
  float            pptol  = 0.01f;
  int              seq_L;

  /* SIMD functions: the SSE engines, or the dispatched ones */
  P7_FS_OPROFILE *(*oprofile_create)(int, const ESL_ALPHABET *, int)                                  = p7_fs_oprofile_Create;
  int             (*oprofile_convert)(const P7_FS_PROFILE *, P7_FS_OPROFILE *)                         = p7_fs_oprofile_Convert;
  void            (*oprofile_destroy)(P7_FS_OPROFILE *)                                                = p7_fs_oprofile_Destroy;
  P7_OMX         *(*omx_create)(int, int, int, int)                                                    = p7_omx_Create_dpf;
  int             (*omx_deconvert)(P7_OMX *, P7_GMX *)                                                 = p7_omx_FDeconvert;
  void            (*omx_destroy)(P7_OMX *)                                                             = p7_omx_Destroy;
  P7_OIVX        *(*oivx_create)(int, int)                                                             = p7_oivx_Create;
  void            (*oivx_destroy)(P7_OIVX *)                                                           = p7_oivx_Destroy;
  int             (*forward)(const ESL_DSQ *, int, const P7_FS_OPROFILE *, P7_OMX *, P7_OIVX *, float *) = p7_Forward_Frameshift_5Codons;
  int             (*backward)(const ESL_DSQ *, int, const P7_FS_OPROFILE *, const P7_OMX *, P7_OMX *, P7_OIVX *, float *) = p7_Backward_Frameshift_5Codons;
  int             (*decoding)(const P7_FS_OPROFILE *, P7_OMX *, const P7_OMX *)                        = p7_Decoding_Frameshift;
  int             (*oa_fill)(const P7_FS_OPROFILE *, const P7_OMX *, P7_OMX *, float *)                = p7_OptimalAccuracy_Frameshift;
  int             (*oa_trace)(const P7_FS_OPROFILE *, const float *, const P7_OMX *, P7_TRACE *)        = p7_OATrace_Frameshift;

#ifdef eslENABLE_SSE
  if (use_sse)
    {
      oprofile_create  = p7_fs_oprofile_Create_sse;
      oprofile_convert = p7_fs_oprofile_Convert_sse;
      oprofile_destroy = p7_fs_oprofile_Destroy_sse;
      omx_create       = p7_omx_Create_dpf_sse;
      omx_deconvert    = p7_omx_FDeconvert_sse;
      omx_destroy      = p7_omx_Destroy_sse;
      oivx_create      = p7_oivx_Create_sse;
      oivx_destroy     = p7_oivx_Destroy_sse;
      forward          = p7_Forward_Frameshift_5Codons_sse;
      backward         = p7_Backward_Frameshift_5Codons_sse;
      decoding         = p7_Decoding_Frameshift_sse;
      oa_fill          = p7_OptimalAccuracy_Frameshift_sse;
      oa_trace         = p7_OATrace_Frameshift_sse;
    }
#endif

  p7_FLogsumInit();

  if (p7_hmm_Sample(r, M, abcAA, &hmm)                              != eslOK) esl_fatal(msg);

  gm_fs = p7_profile_fs_Create(hmm->M, abcAA, p7P_5CODONS);
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs, L, p7_UNILOCAL) != eslOK) esl_fatal(msg);

  om_fs = oprofile_create(hmm->M, abcAA, p7P_5CODONS);
  if (oprofile_convert(gm_fs, om_fs)                          != eslOK) esl_fatal(msg);

  seq_L = L;
  dsq   = malloc(sizeof(ESL_DSQ) * (seq_L + 2));

  /* SIMD: Forward and pp use the 8-cell FS layout, Backward 3-cell, OA 8-cell */
  ox_fwd = omx_create(hmm->M, seq_L, seq_L, p7X_NSCELLS_FS);
  ox_bck = omx_create(hmm->M, seq_L, seq_L, p7X_NSCELLS);
  ox_oa  = omx_create(hmm->M, seq_L, seq_L, p7X_NSCELLS_FS);
  P7_OIVX *ov5 = oivx_create(hmm->M, p7P_5CODONS);

  /* Scalar: pp, OA and the deconverted SIMD OA use the 8-cell FS layout */
  gx_pp  = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS_FS);
  gx_oa  = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS_FS);
  gx_oa2 = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS_FS);

  tr_sse = p7_trace_fs_CreateWithPP();
  tr_ref = p7_trace_fs_CreateWithPP();
  tr_ip  = p7_trace_fs_CreateWithPP();
  ppx    = malloc(sizeof(float) * (seq_L+1) * p7X_NXCELLS);   /* p7X_NXCELLS >= p7G_NXCELLS */

  while (N--)
    {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, seq_L, dsq);
      
      /* SSE pipeline: Forward -> Backward -> Decoding (in-place into ox_fwd) -> OA -> Trace */
      if (forward(dsq, seq_L, om_fs,         ox_fwd,      ov5, &fsc) == eslERANGE) continue;
      if (backward(dsq, seq_L, om_fs, ox_fwd, ox_bck, ov5, &bsc) == eslERANGE) continue;
      if (esl_FCompare_old(fsc, bsc, tol)                                      != eslOK) esl_fatal(msg);
      if (decoding(om_fs, ox_fwd, ox_bck)                    != eslOK) esl_fatal(msg);
      /* ox_fwd now holds the SSE posterior probability matrix */
      if (oa_fill(om_fs, ox_fwd, ox_oa, &accscore_sse)  != eslOK) esl_fatal(msg);
      if (oa_trace(om_fs, ox_fwd->xmx, ox_oa, tr_sse)          != eslOK) esl_fatal(msg);

      /* Reference: convert SSE pp to scalar GMX, run scalar OA and trace */
      if (omx_deconvert(ox_fwd, gx_pp)                                    != eslOK) esl_fatal(msg);
      if (p7_GOptimalAccuracy_Frameshift(gm_fs, gx_pp, gx_oa, &accscore_ref)   != eslOK) esl_fatal(msg);
      if (p7_GOATrace_Frameshift(gm_fs, gx_pp->xmx, gx_oa, tr_ref)             != eslOK) esl_fatal(msg);

      /* In place: the scalar OA written over the posteriors, with a copy of
       * their special states for the trace, matches the separate OA exactly */
      esl_vec_FCopy(gx_pp->xmx, (seq_L+1) * p7G_NXCELLS, ppx);
      if (p7_GOptimalAccuracy_Frameshift(gm_fs, gx_pp, gx_pp, &accscore_ip)  != eslOK) esl_fatal(msg);
      if (p7_GOATrace_Frameshift(gm_fs, ppx, gx_pp, tr_ip)                   != eslOK) esl_fatal(msg);
      if (accscore_ip != accscore_ref)                                                 esl_fatal(msg);
      if (oa_fs_compare(gx_pp, gx_oa, 0.0)                                    != eslOK) esl_fatal(msg);
      if (p7_trace_Compare(tr_ip, tr_ref, 0.0)                               != eslOK) esl_fatal(msg);

      /* Convert SSE OA to scalar GMX for matrix-level comparison */
      if (omx_deconvert(ox_oa, gx_oa2)                                     != eslOK) esl_fatal(msg);
      
      if (esl_FCompare_old(accscore_sse, accscore_ref, tol) != eslOK) esl_fatal(msg);
      if (oa_fs_compare(gx_oa, gx_oa2, tol)                != eslOK) esl_fatal(msg);
      if (p7_trace_Compare(tr_sse, tr_ref, pptol)           != eslOK) esl_fatal(msg);

      /* In place: the SIMD OA written over the posteriors, with a copy of
       * their special states for the trace, matches the separate OA exactly */
      p7_trace_Reuse(tr_ip);
      esl_vec_FCopy(ox_fwd->xmx, (seq_L+1) * p7X_NXCELLS, ppx);
      if (oa_fill(om_fs, ox_fwd, ox_fwd, &accscore_simd_ip) != eslOK) esl_fatal(msg);
      if (oa_trace(om_fs, ppx, ox_fwd, tr_ip)                        != eslOK) esl_fatal(msg);
      if (accscore_simd_ip != accscore_sse)                                                esl_fatal(msg);
      if (omx_deconvert(ox_fwd, gx_pp)                                        != eslOK) esl_fatal(msg);
      if (oa_fs_compare(gx_pp, gx_oa2, 0.0)                                       != eslOK) esl_fatal(msg);
      if (p7_trace_Compare(tr_ip, tr_sse, 0.0)                                   != eslOK) esl_fatal(msg);
      
      if (esl_opt_GetBoolean(go, "--traces"))
        {
          p7_trace_fs_Dump(stdout, tr_sse, gm_fs, dsq, abcDNA);
          p7_trace_fs_Dump(stdout, tr_ref, gm_fs, dsq, abcDNA);
        }

      p7_trace_Reuse(tr_sse);
      p7_trace_Reuse(tr_ref);
      p7_trace_Reuse(tr_ip);
    }

  free(dsq);
  p7_trace_fs_Destroy(tr_sse);
  p7_trace_fs_Destroy(tr_ref);
  p7_trace_fs_Destroy(tr_ip);
  free(ppx);
  omx_destroy(ox_fwd);
  omx_destroy(ox_bck);
  omx_destroy(ox_oa);
  oivx_destroy(ov5);
  p7_gmx_Destroy(gx_pp);
  p7_gmx_Destroy(gx_oa);
  p7_gmx_Destroy(gx_oa2);
  oprofile_destroy(om_fs);
  p7_profile_fs_Destroy(gm_fs);
  p7_hmm_Destroy(hmm);
}
/* pp_cell_fs()
 * Value of cell <s> (frameshift layout, e.g. p7X_FS_M + p7X_FS_C0) of
 * model position <k> on row <i> of <ox>, with <use_sse> choosing the
 * SSE or AVX rows.
 */
static float
pp_cell_fs(const P7_OMX *ox, int i, int k, int s, int use_sse)
{
  int Q = use_sse ? p7O_NQF(ox->M) : p7O_NQF_AVX(ox->M);
  int q = (k-1) % Q;
  int r = (k-1) / Q;
  const float *v = use_sse ? (const float *) (ox->dpf[i] + q * p7X_NSCELLS_FS + s) : (const float *) (ox->dpf_avx[i] + q * p7X_NSCELLS_FS + s);
  return v[r];
}

/* utest_sweep_fs()
 * The frameshift OA sweep over a checkpointed Backward against the
 * full-matrix path on the same Backward: p7_Forward_Frameshift_5Codons(), the
 * checkpointed Backward filled in one block and copied into a Backward
 * matrix, p7_Decoding_Frameshift(), p7_OptimalAccuracy_Frameshift(),
 * p7_OATrace_Frameshift() and p7_Null2_fs_ByExpectation(). The sweep,
 * with Backward blocks of 5..12 rows, must give the same Forward and OA
 * scores, null2, states, k, i and codon lengths on the trace, and PP
 * characters equal to those the posterior matrix gives. Sequences
 * alternate random DNA and the back-translated consensus of a model
 * with sharpened match emissions and M->M transitions between random
 * flanks, so that rows are rescaled. With <use_sse>, the _sse versions
 * are called directly, so an AVX CPU also tests the SSE engine.
 */
static void
utest_sweep_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA,
               int M, int L, int N, int use_sse)
{
  char           *msg  = "frameshift OA sweep unit test failed";
  P7_CODONTABLE  *ct   = p7_codontable_Create(gcode);
  P7_HMM         *hmm  = NULL;
  P7_FS_PROFILE  *gm   = p7_profile_fs_Create(M, abcAA, p7P_5CODONS);
  P7_FS_OPROFILE *om   = NULL;
  P7_OIVX        *ov   = NULL;
  ESL_SQ         *sq   = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq  = malloc(sizeof(ESL_DSQ) * (L+2));
  float          *n2a  = malloc(sizeof(float) * abcAA->Kp);
  float          *n2b  = malloc(sizeof(float) * abcAA->Kp);
  float          *ppx  = malloc(sizeof(float) * (L+1) * p7X_NXCELLS);
  P7_OMX         *fwd  = NULL, *bck = NULL, *oa = NULL;
  P7_CHKOMX      *cx1  = NULL, *cx2 = NULL;
  P7_OATMX       *ox   = p7_oatmx_Create(M, L, p7OAT_NPLANES_FS, 1);
  P7_TRACE       *tr1  = p7_trace_fs_CreateWithPP();
  P7_TRACE       *tr2  = p7_trace_fs_CreateWithPP();
  int             Q    = use_sse ? p7O_NQF(M) : p7O_NQF_AVX(M);
  size_t          rowb = (use_sse ? 16 : 32) * p7X_NSCELLS * (size_t) Q;
  float           fsc, fsc2, oasc1, oasc2, p;
  int             n, Ln, blockR, i, j, k, x, z, ntested = 0;

  if (p7_hmm_Sample(r, M, abcAA, &hmm) != eslOK) esl_fatal(msg);
  for (k = 1; k <= M; k++)
    {
      x = esl_vec_FArgMax(hmm->mat[k], abcAA->K);
      esl_vec_FSet(hmm->mat[k], abcAA->K, 0.15 / (abcAA->K - 1));
      hmm->mat[k][x] = 0.85;
      if (k < M) { hmm->t[k][p7H_MM] = 0.9; hmm->t[k][p7H_MI] = 0.05; hmm->t[k][p7H_MD] = 0.05; }
    }
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm, L/3, p7_LOCAL) != eslOK) esl_fatal(msg);
  p7_emit_SimpleConsensus(hmm, sq);

  if (use_sse)
    {
      om  = p7_fs_oprofile_Create_sse(M, abcAA, p7P_5CODONS);  p7_fs_oprofile_Convert_sse(gm, om);
      ov  = p7_oivx_Create_sse(M, p7P_5CODONS);
      fwd = p7_omx_Create_dpf_sse(M, L, L, p7X_NSCELLS_FS); bck = p7_omx_Create_dpf_sse(M, L, L, p7X_NSCELLS); oa = p7_omx_Create_dpf_sse(M, L, L, p7X_NSCELLS_FS);
      cx1 = p7_chkomx_Create_sse(M, L, L+1, 5, p7X_NSCELLS, p7_CHK_BACKWARD);
      cx2 = p7_chkomx_Create_sse(M, L, 5,   5, p7X_NSCELLS, p7_CHK_BACKWARD);
    }
  else
    {
      om  = p7_fs_oprofile_Create(M, abcAA, p7P_5CODONS);  p7_fs_oprofile_Convert(gm, om);
      ov  = p7_oivx_Create(M, p7P_5CODONS);
      fwd = p7_omx_Create_dpf(M, L, L, p7X_NSCELLS_FS); bck = p7_omx_Create_dpf(M, L, L, p7X_NSCELLS); oa = p7_omx_Create_dpf(M, L, L, p7X_NSCELLS_FS);
      cx1 = p7_chkomx_Create(M, L, L+1, 5, p7X_NSCELLS, p7_CHK_BACKWARD);
      cx2 = p7_chkomx_Create(M, L, 5,   5, p7X_NSCELLS, p7_CHK_BACKWARD);
    }

  for (n = 0; n < N; n++)
    {
      Ln = 15 + esl_rnd_Roll(r, L - 14);
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, dsq);
      if (n % 2 == 1 && 3*M <= L)
        {
          Ln = ESL_MAX(Ln, 3*M);
          esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, dsq);
          for (i = 1, j = 1 + esl_rnd_Roll(r, Ln - 3*M + 1); i <= M; i++, j += 3)
            p7_codontable_GetCodon(ct, r, sq->dsq[i], dsq+j);
        }
      dsq[0] = dsq[Ln+1] = eslDSQ_SENTINEL;
      blockR = 5 + esl_rnd_Roll(r, 8);
      p7_trace_Reuse(tr1);
      p7_trace_Reuse(tr2);

      if (use_sse)
        {
          p7_fs_oprofile_ReconfigLength_sse(om, Ln/3);
          if (p7_Forward_Frameshift_5Codons_sse(dsq, Ln, om, fwd, ov, &fsc)       != eslOK) continue;
          if (p7_chkomx_GrowTo_sse(cx1, M, Ln, Ln+1)                       != eslOK) esl_fatal(msg);
          if (p7_BackwardChk_Frameshift_5Codons_sse(dsq, Ln, om, cx1, ov, NULL)    != eslOK) continue;
          for (i = 0; i <= Ln; i++) memcpy(bck->dpf[i], cx1->dpf[i], rowb);
        }
      else
        {
          p7_fs_oprofile_ReconfigLength(om, Ln/3);
          if (p7_Forward_Frameshift_5Codons(dsq, Ln, om, fwd, ov, &fsc)           != eslOK) continue;
          if (p7_chkomx_GrowTo(cx1, M, Ln, Ln+1)                           != eslOK) esl_fatal(msg);
          if (p7_BackwardChk_Frameshift_5Codons(dsq, Ln, om, cx1, ov, NULL)        != eslOK) continue;
          for (i = 0; i <= Ln; i++) memcpy(bck->dpf_avx[i], cx1->dpf_avx[i], rowb);
        }
      memcpy(bck->xmx, cx1->xmx, sizeof(float) * (Ln+1) * p7X_NXCELLS);
      bck->M = M; bck->L = Ln;

      if (use_sse)
        {
          if (p7_Decoding_Frameshift_sse(om, fwd, bck)                     != eslOK) esl_fatal(msg);
          memcpy(ppx, fwd->xmx, sizeof(float) * (Ln+1) * p7X_NXCELLS);
          p7_OptimalAccuracy_Frameshift_sse(om, fwd, oa, &oasc1);
          if (p7_OATrace_Frameshift_sse(om, ppx, oa, tr1)                  != eslOK) esl_fatal(msg);
          if (p7_chkomx_GrowTo_sse(cx2, M, Ln, blockR)                     != eslOK) esl_fatal(msg);
          if (p7_BackwardChk_Frameshift_5Codons_sse(dsq, Ln, om, cx2, ov, NULL)    != eslOK) esl_fatal(msg);
          if (p7_OptimalAccuracySweep_Frameshift_sse(dsq, Ln, om, cx2, ov, ox, n2b, &fsc2, &oasc2) != eslOK) esl_fatal(msg);
          if (p7_OATraceSweep_Frameshift_sse(om, ox, tr2)                  != eslOK) esl_fatal(msg);
        }
      else
        {
          if (p7_Decoding_Frameshift(om, fwd, bck)                         != eslOK) esl_fatal(msg);
          memcpy(ppx, fwd->xmx, sizeof(float) * (Ln+1) * p7X_NXCELLS);
          p7_OptimalAccuracy_Frameshift(om, fwd, oa, &oasc1);
          if (p7_OATrace_Frameshift(om, ppx, oa, tr1)                      != eslOK) esl_fatal(msg);
          if (p7_chkomx_GrowTo(cx2, M, Ln, blockR)                         != eslOK) esl_fatal(msg);
          if (p7_BackwardChk_Frameshift_5Codons(dsq, Ln, om, cx2, ov, NULL)        != eslOK) esl_fatal(msg);
          if (p7_OptimalAccuracySweep_Frameshift(dsq, Ln, om, cx2, ov, ox, n2b, &fsc2, &oasc2) != eslOK) esl_fatal(msg);
          if (p7_OATraceSweep_Frameshift(om, ox, tr2)                      != eslOK) esl_fatal(msg);
        }

      if (fsc2  != fsc)    esl_fatal("%s: Forward scores differ (%f %f)", msg, fsc2, fsc);
      if (oasc2 != oasc1)  esl_fatal("%s: OA scores differ (%f %f)", msg, oasc2, oasc1);
      if (tr1->N != tr2->N) esl_fatal("%s: trace lengths differ (%d %d)", msg, tr1->N, tr2->N);
      for (z = 0; z < tr1->N; z++)
        {
          if (tr1->st[z] != tr2->st[z] || tr1->k[z] != tr2->k[z] || tr1->i[z] != tr2->i[z] || tr1->c[z] != tr2->c[z])
            esl_fatal("%s: traces differ at %d", msg, z);
          switch (tr1->st[z]) {
          case p7T_M: p = pp_cell_fs(fwd, tr1->i[z], tr1->k[z], p7X_FS_M + p7X_FS_C0, use_sse); break;
          case p7T_I: p = pp_cell_fs(fwd, tr1->i[z], tr1->k[z], p7X_FS_I,              use_sse); break;
          case p7T_N: case p7T_C: case p7T_J:
            /* an N, C or J emits residue i when the state before it is the same one */
            p = (z > 0 && tr1->st[z-1] == tr1->st[z]) ? ppx[tr1->i[z]*p7X_NXCELLS + (tr1->st[z] == p7T_N ? p7X_N : tr1->st[z] == p7T_C ? p7X_C : p7X_J)] : 0.0;
            break;
          default: p = 0.0; break;
          }
          if (p7_alidisplay_EncodePostProb(tr2->pp[z]) != p7_alidisplay_EncodePostProb(p))
            esl_fatal("%s: PP character differs at %d (state %d, i %d, k %d)", msg, z, tr1->st[z], tr1->i[z], tr1->k[z]);
        }

      /* null2 from the full posterior matrix, which it uses as workspace, so last */
      if (use_sse) p7_Null2_fs_ByExpectation_sse(om, fwd, n2a);
      else         p7_Null2_fs_ByExpectation    (om, fwd, n2a);
      for (x = 0; x < abcAA->Kp; x++)
        if (n2b[x] != n2a[x]) esl_fatal("%s: null2 differs for residue %d (%g %g)", msg, x, n2b[x], n2a[x]);
      ntested++;
    }
  if (ntested == 0) esl_fatal("%s: no sequence was tested", msg);

  free(dsq); free(n2a); free(n2b); free(ppx);
  esl_sq_Destroy(sq);
  p7_trace_fs_Destroy(tr1);
  p7_trace_fs_Destroy(tr2);
  p7_oatmx_Destroy(ox);
  if (use_sse)
    {
      p7_omx_Destroy_sse(fwd); p7_omx_Destroy_sse(bck); p7_omx_Destroy_sse(oa);
      p7_chkomx_Destroy_sse(cx1); p7_chkomx_Destroy_sse(cx2);
      p7_oivx_Destroy_sse(ov);
      p7_fs_oprofile_Destroy_sse(om);
    }
  else
    {
      p7_omx_Destroy(fwd); p7_omx_Destroy(bck); p7_omx_Destroy(oa);
      p7_chkomx_Destroy(cx1); p7_chkomx_Destroy(cx2);
      p7_oivx_Destroy(ov);
      p7_fs_oprofile_Destroy(om);
    }
  p7_profile_fs_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_codontable_Destroy(ct);
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

  utest_optacc_fs(go, r, abcAA, abcDNA, gcode, bgAA, bgDNA, M, L, N, FALSE);
#ifdef eslENABLE_SSE
  utest_optacc_fs(go, r, abcAA, abcDNA, gcode, bgAA, bgDNA, M, L, N, TRUE);
#endif
  utest_sweep_fs(r, abcAA, abcDNA, gcode, bgAA, bgDNA, M, 3*M+150, N,  FALSE);
  utest_sweep_fs(r, abcAA, abcDNA, gcode, bgAA, bgDNA, 1, 60,      10, FALSE);
#ifdef eslENABLE_SSE
  utest_sweep_fs(r, abcAA, abcDNA, gcode, bgAA, bgDNA, M, 3*M+150, N,  TRUE);
  utest_sweep_fs(r, abcAA, abcDNA, gcode, bgAA, bgDNA, 1, 60,      10, TRUE);
#endif

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




