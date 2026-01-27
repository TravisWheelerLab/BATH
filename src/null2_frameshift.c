/* "null2" model: biased composition correction
 * 
 * Contents:
 *   1. Null2 estimation algorithms.
 *   2. Benchmark driver.
 *
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

/*****************************************************************
 * 1. Null2 estimation algorithms.
 *****************************************************************/

/* Function:  p7_Null2_fs_ByExpectation()
 * Synopsis:  Calculate null2 model from posterior probabilities.
 *
 * Purpose:   Calculate the "null2" model for the envelope encompassed
 *            by a posterior probability calculation <pp> for frameshift
 *            aware model <gm_fs>. Return the null2 odds amino acide
 *            emission probabilities $\frac{f'{x}}{f{x}}$ in <null2>,
 *            which caller provides as space for at least <alphabet->Kp>
 *            residues.
 *
 *            Make sure that the posterior probability matrix <pp> has
 *            been calculated by the caller for only the envelope; thus
 *            its rows are numbered <1..Ld>, for envelope <ienv..jenv>
 *            of length <Ld=jenv-ienv+1>.
 *
 * Args:      gm_fs - profile, in any mode, target length model set to <L>
 *            pp    - posterior prob matrix, for <gm> against domain envelope <dsq+i-1> (offset)
 *            null2 - RETURN: null2 odds ratios per residue; <0..Kp-1>; caller allocated space
 *
 * Returns:   <eslOK> on success; <null2> contains the null2 scores. The 0
 *            row of <pp> has been used as temp space, and happens to contain
 *            the expected frequency that each M,I,N,C,J state is used in this
 *            <pp> matrix to generate residues.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_Null2_fs_ByExpectation(const P7_FS_PROFILE *gm_fs, P7_GMX *pp, float *null2)
{
  int      M      = gm_fs->M;
  int      Ld     = pp->L;
  float  **dp     = pp->dp;
  float   *xmx    = pp->xmx;
  float    xfactor;
  int      i;			/* over offset envelope dsq positions 1..Ld  */
  int      k;			/* over model M states 1..M, I states 1..M-1 */
  int      x;

  /* Calculate expected # of times that each emitting state was used
   * in generating the Ld residues in this domain.
   * The 0 row in <wrk> is used to hold these numbers.
   */
    
  esl_vec_FCopy(pp->dp[1],            (M+1)*p7G_NSCELLS_FS, pp->dp[0]); 
  esl_vec_FCopy(pp->xmx+p7G_NXCELLS,        p7G_NXCELLS,    pp->xmx);  
  for (i = 2; i <= Ld; i++)
  {
      esl_vec_FAdd(pp->dp[0], pp->dp[i],             (M+1)*p7G_NSCELLS_FS);
      esl_vec_FAdd(pp->xmx,   pp->xmx+i*p7G_NXCELLS,       p7G_NXCELLS); 
  }
 
  /* Convert those expected #'s to log frequencies; these we'll use as
   * the log posterior weights.
   */
  esl_vec_FLog(pp->dp[0], (M+1)*p7G_NSCELLS_FS);
  esl_vec_FLog(pp->xmx,   p7G_NXCELLS);  

  esl_vec_FIncrement(pp->dp[0], (M+1)*p7G_NSCELLS_FS, -log((float)Ld));
  esl_vec_FIncrement(pp->xmx,   p7G_NXCELLS,       -log((float)Ld)); 

  /* Calculate null2's log odds emission probabilities, by taking
   * posterior weighted sum over all emission vectors used in paths
   * explaining the domain.
   * This is dog-slow; a point for future optimization.
   */
  xfactor = XMX_FS(0,p7G_N);
  xfactor = p7_FLogsum(xfactor, XMX_FS(0,p7G_C));
  xfactor = p7_FLogsum(xfactor, XMX_FS(0,p7G_J));
  
  esl_vec_FSet(null2, gm_fs->abc->K, -eslINFINITY);

  for (x = 0; x < gm_fs->abc->K; x++)
  {
      for (k = 1; k < M; k++)
        {
          null2[x] = p7_FLogsum(null2[x], MMX_FS(0,k,p7G_C0) + p7P_MSC_AMINO(gm_fs, k, x));
          null2[x] = p7_FLogsum(null2[x], IMX_FS(0,k));
        }
      null2[x] = p7_FLogsum(null2[x], MMX_FS(0,M,p7G_C0) + p7P_MSC_AMINO(gm_fs, M, x));
      null2[x] = p7_FLogsum(null2[x], xfactor);
    }
    
  esl_vec_FExp (null2, gm_fs->abc->K);
  /* now null2[x] = \frac{f_d(x)}{f_0(x)} for all x in alphabet,
   * 0..K-1, where f_d(x) are the ad hoc "null2" residue frequencies
   * for this envelope.
   */

  /* make valid scores for all degeneracies, by averaging the odds ratios. */
  esl_abc_FAvgScVec(gm_fs->abc, null2); /* does not set gap, nonres, missing  */
  null2[gm_fs->abc->K]    = 1.0;        /* gap character    */
  null2[gm_fs->abc->Kp-2] = 1.0;        /* nonresidue "*"   */
  null2[gm_fs->abc->Kp-1] = 1.0;        /* missing data "~" */

  return eslOK;
}



/*****************************************************************
 * 2. Benchmark driver
 *****************************************************************/
#ifdef p7NULL2_FRAMESHIFT_BENCHMARK
/*
    gcc -g -O2      -o null2_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7NULL2_FRAMESHIFT_BENCHMARK null2_frameshift.c -lhmmer -leasel -lm
   icc -O3 -static -o null2_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7NULL2_FRAMESHIFT_BENCHMARK null2_frameshift.c -lhmmer -leasel -lm
   ./null2_frameshift_benchmark <hmmfile>
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,   "1200", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for posterior residue null2, generic version";

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
  P7_FS_PROFILE  *gm_fs   = NULL;
  P7_GMX         *gx1     = NULL;
  P7_GMX         *gx2     = NULL;
  P7_GMX         *pp      = NULL;
  P7_IVX         *iv      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  float           null2[p7_MAXCODE];
  int             i;
  float           fsc, bsc;
  double          Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)          != eslOK) p7_Fail("Failed to read HMM");

  abcDNA = esl_alphabet_Create(eslDNA);
  gcode  = esl_gencode_Create(abcDNA, abcAA);
  bgDNA  = p7_bg_Create(abcDNA);                p7_bg_SetLength(bgDNA, L);
  bgAA   = p7_bg_Create(abcAA);                 p7_bg_SetLength(bgAA, L/3);
  gm_fs  = p7_profile_fs_Create(hmm->M, abcAA); p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs, L/3, p7_LOCAL);
  gx1    = p7_gmx_fs_Create(gm_fs->M, L, L, p7P_5CODONS);  
  gx2    = p7_gmx_fs_Create(gm_fs->M, L, L, 0);
  pp     = p7_gmx_fs_Create(gm_fs->M, L, L, p7P_5CODONS);
  iv     = p7_ivx_Create(gm_fs->M, p7P_5CODONS);

  esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
  p7_Forward_Frameshift (dsq, gcode, L, gm_fs, gx1, iv, &fsc);
  p7_Backward_Frameshift(dsq, gcode, L, gm_fs, gx2, iv, &bsc);
  p7_Decoding_Frameshift(gm_fs, gx1, gx2, pp);   

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) 
    p7_Null2_fs_ByExpectation(gm_fs, pp, null2);   
  esl_stopwatch_Stop(w);

  Mcs  = (double) N * (double) L * (double) gm_fs->M * 1e-6 / w->user;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n", gm_fs->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_gmx_Destroy(gx1);
  p7_gmx_Destroy(gx2);
  p7_gmx_Destroy(pp);
  p7_ivx_Destroy(iv);
  p7_profile_fs_Destroy(gm_fs);
  p7_bg_Destroy(bgAA);
  p7_bg_Destroy(bgDNA);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abcAA);
  esl_alphabet_Destroy(abcDNA);
  esl_gencode_Destroy(gcode);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7NULL2_FRAMESHIFT_BENCHMARK*/
/*------------------ end, benchmark driver ----------------------*/



