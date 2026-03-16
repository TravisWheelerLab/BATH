/* Posterior decoding algorithms; frameshift versions.
 * 
 * Contents:
 *   1. Posterior decoding algorithms.
 *   2. Benchmark driver.
 */
#include "p7_config.h"

#include <math.h>
#include "easel.h"
#include "hmmer.h"

/*****************************************************************
 * 1. Posterior decoding algorithms.
 *****************************************************************/

/* Function:  p7_Decoding_Frameshift() - BATH
 * Synopsis:  Posterior decoding of residue assignments.
 *
 * Purpose:   Calculates a posterior decoding of the residues and 
 *            of the codons in a DNA target sequence, given framshift 
 *            aware codon profile  <gm_fs5> and filled Forward and 
 *            Backward matrices <fwd>, <bck> for the profile aligned 
 *            to that target sequence. The codon posterior decoding 
 *            is overwriten to the <fwd> matrix.
 *            
 * Args:      gm_fs5 - frameshift codon profile 
 *            fwd   - filled Forward matrix - get overwriten 
 *            bck   - filled Backward matrix
 *
 * Returns:   <eslOK> on success.
 *
 *    
 */
int
p7_GDecoding_Frameshift(const P7_FS_PROFILE *gm_fs5, P7_GMX *fwd, P7_GMX *bck)
{
  float      **dp   = fwd->dp;
  float       *xmx  = fwd->xmx;
  float        N0, N1, N2, N3;
  float        J0, J1, J2, J3;
  float        C0, C1, C2, C3;
  int          L    = fwd->L;
  int          M    = gm_fs5->M;
  int          i,k;
  float overall_sc = p7_FLogsum( bck->xmx[0*p7G_NXCELLS+p7G_N],
                     p7_FLogsum( bck->xmx[1*p7G_NXCELLS+p7G_N],
                                 bck->xmx[2*p7G_NXCELLS+p7G_N]));
  float        denom;
  
  N0 = fwd->xmx[p7G_NXCELLS*0 + p7G_N];
  J0 = fwd->xmx[p7G_NXCELLS*0 + p7G_J];
  C0 = fwd->xmx[p7G_NXCELLS*0 + p7G_C];

  N1 = N2 = N3 = 0.0;
  J1 = J2 = J3 = 0.0;
  C1 = C2 = C3 = 0.0;

  XMX_FS(0, p7G_E) = 0.0; 
  XMX_FS(0, p7G_N) = 0.0; 
  XMX_FS(0, p7G_J) = 0.0;
  XMX_FS(0, p7G_B) = 0.0;
  XMX_FS(0, p7G_C) = 0.0;

  for (k = 0; k <= M; k++) {
    MMX_FS(0,k,p7G_C1) = MMX_FS(0,k,p7G_C2) = MMX_FS(0,k,p7G_C3) = MMX_FS(0,k,p7G_C4) = MMX_FS(0,k,p7G_C5) = 0.0;
    MMX_FS(0,k,p7G_C0) = IMX_FS(0,k) = DMX_FS(0,k) = 0.0;
  }

  for (i = 1; i <= L; i++)
    {
	 N3 = N2; N2 = N1; N1 = N0;
     J3 = J2; J2 = J1; J1 = J0;
	 C3 = C2; C2 = C1; C1 = C0;

     denom = 0.0;
     /* 0th position in model */    
	  MMX_FS(i,0,p7G_C1) = MMX_FS(i,0,p7G_C2) = MMX_FS(i,0,p7G_C3) = MMX_FS(i,0,p7G_C4) = MMX_FS(i,0,p7G_C5) = 0.0;
      MMX_FS(i,0, p7G_C0) = IMX_FS(i,0) = DMX_FS(i,0) = 0.0;

      for (k = 1; k < M; k++)
      {
        MMX_FS(i,k,p7G_C0) = expf(fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C0] + bck->dp[i][k*p7G_NSCELLS + p7G_M] - overall_sc); denom += MMX_FS(i,k, p7G_C0);
        MMX_FS(i,k,p7G_C1) = expf(fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1] + bck->dp[i][k*p7G_NSCELLS + p7G_M] - overall_sc); 
        MMX_FS(i,k,p7G_C2) = expf(fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2] + bck->dp[i][k*p7G_NSCELLS + p7G_M] - overall_sc);
        MMX_FS(i,k,p7G_C3) = expf(fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3] + bck->dp[i][k*p7G_NSCELLS + p7G_M] - overall_sc);
        MMX_FS(i,k,p7G_C4) = expf(fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4] + bck->dp[i][k*p7G_NSCELLS + p7G_M] - overall_sc);
        MMX_FS(i,k,p7G_C5) = expf(fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C5] + bck->dp[i][k*p7G_NSCELLS + p7G_M] - overall_sc);        

        IMX_FS(i,k) = expf(fwd->dp[i][k*p7G_NSCELLS_FS + p7G_I] + bck->dp[i][k*p7G_NSCELLS + p7G_I] - overall_sc); denom += IMX_FS(i,k); 

        DMX_FS(i,k) = 0.0;
      }

      MMX_FS(i,M, p7G_C0) = expf(fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C0] + bck->dp[i][M*p7G_NSCELLS + p7G_M] - overall_sc); denom += MMX_FS(i,M, p7G_C0); 
      MMX_FS(i,M, p7G_C1) = expf(fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C1] + bck->dp[i][M*p7G_NSCELLS + p7G_M] - overall_sc);
      MMX_FS(i,M, p7G_C2) = expf(fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C2] + bck->dp[i][M*p7G_NSCELLS + p7G_M] - overall_sc);
      MMX_FS(i,M, p7G_C3) = expf(fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C3] + bck->dp[i][M*p7G_NSCELLS + p7G_M] - overall_sc);
      MMX_FS(i,M, p7G_C4) = expf(fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C4] + bck->dp[i][M*p7G_NSCELLS + p7G_M] - overall_sc);
      MMX_FS(i,M, p7G_C5) = expf(fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C5] + bck->dp[i][M*p7G_NSCELLS + p7G_M] - overall_sc);
                       
      IMX_FS(i,M) = 0.0; 
      DMX_FS(i,M) = 0.0;

      /* no emission from E or B states */
      XMX_FS(i,p7G_E)  = 0.0;
      XMX_FS(i,p7G_B)  = 0.0; 

      /* probaility from N, J and C states */
	  N0 = fwd->xmx[p7G_NXCELLS*i + p7G_N];
	  J0 = fwd->xmx[p7G_NXCELLS*i + p7G_J];
	  C0 = fwd->xmx[p7G_NXCELLS*i + p7G_C];
      if(i > 2)
      {
         XMX_FS(i,p7G_N) = expf(N3 + bck->xmx[p7G_NXCELLS*i + p7G_N] + gm_fs5->xsc[p7P_N][p7P_LOOP] - overall_sc);

         XMX_FS(i,p7G_C) = expf(C3 + bck->xmx[p7G_NXCELLS*i + p7G_C] + gm_fs5->xsc[p7P_C][p7P_LOOP] - overall_sc);
         
         XMX_FS(i,p7G_J) = expf(J3 + bck->xmx[p7G_NXCELLS*i + p7G_J] + gm_fs5->xsc[p7P_J][p7P_LOOP] - overall_sc);

         denom += XMX_FS(i,p7G_N) + XMX_FS(i,p7G_J) + XMX_FS(i,p7G_C);
      }
      else
      {
         XMX_FS(i,p7G_N) = expf(bck->xmx[p7G_NXCELLS*i + p7G_N] - overall_sc);
         XMX_FS(i,p7G_C) = 0.0;
         XMX_FS(i,p7G_J) = 0.0;
         denom += XMX_FS(i,p7G_N);
      }

        denom = 1.0 / denom;
 
        for (k = 1; k < M; k++) { 
          MMX_FS(i,k, p7G_C0) *= denom; 
          MMX_FS(i,k, p7G_C1) *= denom;
          MMX_FS(i,k, p7G_C2) *= denom;
          MMX_FS(i,k, p7G_C3) *= denom;
          MMX_FS(i,k, p7G_C4) *= denom;
          MMX_FS(i,k, p7G_C5) *= denom;

          IMX_FS(i,k) *= denom;
        }
        MMX_FS(i,M, p7G_C0) *= denom;
        MMX_FS(i,M, p7G_C1) *= denom;
        MMX_FS(i,M, p7G_C2) *= denom;
        MMX_FS(i,M, p7G_C3) *= denom;
        MMX_FS(i,M, p7G_C4) *= denom;
        MMX_FS(i,M, p7G_C5) *= denom;

        XMX_FS(i,p7G_N) *= denom;
        XMX_FS(i,p7G_C) *= denom;
        XMX_FS(i,p7G_J) *= denom;
    }

  return eslOK;
}




/* Function:  p7_DomainDecoding_Frameshift
 * Synopsis:  Posterior decoding of domain location for frameshift aware alignments.
 *
 * Purpose:   The caller has calculated Forward and Backward matrices
 *            <fwd> and <bck> for frameshift codon model <gm_fs5> 
 *            aligned to a target sequence of length <L>.
 * 
 *            We use this information to calculate the posterior
 *            probabilities that we're in a begin state B, end state
 *            E, or any core model state {M,D,I} at each target
 *            sequence position <i = 1..L>.
 * 
 *            <ddef->btot[i]> stores the cumulative expectation
 *            $\sum_1^i$ of the number of i's that were emitted (by an
 *            Mk state) immediately after a B : i.e., the expected
 *            number of times domains have started at or before
 *            position i.
 * 
 *            <ddef->etot[i]> stores the cumulative expectation
 *            $\sum_1^i$ of the number of i's that were emitted (by
 *            an Mk or Dk state) and immediately followed by an end
 *            transition to E : i.e., the expected number of times
 *            domains have ended at or before position i.
 * 
 *            <ddef->mocc[i]> stores the probability that residue i is
 *            emitted by the core model, as opposed to the flanking
 *            N,C,J states : i.e., the probability that i is in a
 *            domain.
 * 
 *            Upon return, each of these arrays has been made, and
 *            <ddef->L> has * been set.
 *
 * Args:      gm_fs5 - profile
 *            fwd   - filled Forward matrix
 *            bck   - filled Backward matrix
 *            ddef  - container for the results.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 * 
 */
int
p7_GDomainDecoding_Frameshift(const P7_FS_PROFILE *gm_fs5, const P7_GMX *fwd, const P7_GMX *bck, P7_DOMAINDEF *ddef)
{
  int   L            = fwd->L;
  float overall_logp = p7_FLogsum( bck->xmx[0*p7G_NXCELLS+p7G_N],
                       p7_FLogsum( bck->xmx[1*p7G_NXCELLS+p7G_N],
                                   bck->xmx[2*p7G_NXCELLS+p7G_N]));
  float njcp; 
  int   i;

  /* First three psoitions are set to 0. They will be included in a domain that starts at i=3 */
  ddef->btot[0] = 0.;
  ddef->btot[1] = 0.;
  ddef->btot[2] = 0.;
  ddef->etot[0] = 0.;
  ddef->etot[1] = 0.;
  ddef->etot[2] = 0.;

  for (i = 3; i <= L; i++)
  {
    /* Accumulate probabiliies at the B and E states - proabilities of transitioning in or out of the core martix */
    ddef->btot[i] = ddef->btot[i-3] + expf(fwd->xmx[(i-3)*p7G_NXCELLS+p7G_B] + bck->xmx[(i-3)*p7G_NXCELLS+p7G_B] - overall_logp); 

    ddef->etot[i] = ddef->etot[i-3] + expf(fwd->xmx[i*p7G_NXCELLS+p7G_E] + bck->xmx[i*p7G_NXCELLS+p7G_E] - overall_logp);;
  }  

  ddef->mocc[0] = 0.;
  ddef->mocc[1] = 0.;
  ddef->mocc[2] = 0.;

  for (i = 3; i < L-1; i++)
  {

    /* Sum poropabilities in the N, J and C states for all codons in which i is present */
    njcp = 0.0;
  
    /*N state */
    njcp += expf(fwd->xmx[(i-3)*p7G_NXCELLS+p7G_N] + bck->xmx[i*p7G_NXCELLS+p7G_N]     + gm_fs5->xsc[p7P_N][p7P_LOOP] - overall_logp);
    njcp += expf(fwd->xmx[(i-2)*p7G_NXCELLS+p7G_N] + bck->xmx[(i+1)*p7G_NXCELLS+p7G_N] + gm_fs5->xsc[p7P_N][p7P_LOOP] - overall_logp);
    njcp += expf(fwd->xmx[(i-1)*p7G_NXCELLS+p7G_N] + bck->xmx[(i+2)*p7G_NXCELLS+p7G_N] + gm_fs5->xsc[p7P_N][p7P_LOOP] - overall_logp);    

    /* C state */
    njcp += expf(fwd->xmx[(i-3)*p7G_NXCELLS+p7G_C] + bck->xmx[i*p7G_NXCELLS+p7G_C]     + gm_fs5->xsc[p7P_C][p7P_LOOP] - overall_logp);
    njcp += expf(fwd->xmx[(i-2)*p7G_NXCELLS+p7G_C] + bck->xmx[(i+1)*p7G_NXCELLS+p7G_C] + gm_fs5->xsc[p7P_C][p7P_LOOP] - overall_logp);
    njcp += expf(fwd->xmx[(i-1)*p7G_NXCELLS+p7G_C] + bck->xmx[(i+2)*p7G_NXCELLS+p7G_C] + gm_fs5->xsc[p7P_C][p7P_LOOP] - overall_logp);

    /* J state */
    njcp += expf(fwd->xmx[(i-3)*p7G_NXCELLS+p7G_J] + bck->xmx[i*p7G_NXCELLS+p7G_J]     + gm_fs5->xsc[p7P_J][p7P_LOOP] - overall_logp);
    njcp += expf(fwd->xmx[(i-2)*p7G_NXCELLS+p7G_J] + bck->xmx[(i+1)*p7G_NXCELLS+p7G_J] + gm_fs5->xsc[p7P_J][p7P_LOOP] - overall_logp);
    njcp += expf(fwd->xmx[(i-1)*p7G_NXCELLS+p7G_J] + bck->xmx[(i+2)*p7G_NXCELLS+p7G_J] + gm_fs5->xsc[p7P_J][p7P_LOOP] - overall_logp);
    
    /* Probability of i emitted by the core model is aaprixmated as 1 - probability of i emited by the specials */
    ddef->mocc[i] = 1. - njcp;
  }
  njcp = 0.0;

  /*N state */
  njcp += expf(fwd->xmx[(L-4)*p7G_NXCELLS+p7G_N] + bck->xmx[(L-1)*p7G_NXCELLS+p7G_N]     + gm_fs5->xsc[p7P_N][p7P_LOOP] - overall_logp);
  njcp += expf(fwd->xmx[(L-3)*p7G_NXCELLS+p7G_N] + bck->xmx[L*p7G_NXCELLS+p7G_N] + gm_fs5->xsc[p7P_N][p7P_LOOP] - overall_logp);

  /* C state */
  njcp += expf(fwd->xmx[(L-4)*p7G_NXCELLS+p7G_C] + bck->xmx[(L-1)*p7G_NXCELLS+p7G_C]     + gm_fs5->xsc[p7P_C][p7P_LOOP] - overall_logp);
  njcp += expf(fwd->xmx[(L-3)*p7G_NXCELLS+p7G_C] + bck->xmx[L*p7G_NXCELLS+p7G_C] + gm_fs5->xsc[p7P_C][p7P_LOOP] - overall_logp);

  /* J state */
  njcp += expf(fwd->xmx[(L-4)*p7G_NXCELLS+p7G_J] + bck->xmx[(L-1)*p7G_NXCELLS+p7G_J]     + gm_fs5->xsc[p7P_J][p7P_LOOP] - overall_logp);
  njcp += expf(fwd->xmx[(L-3)*p7G_NXCELLS+p7G_J] + bck->xmx[L*p7G_NXCELLS+p7G_J] + gm_fs5->xsc[p7P_J][p7P_LOOP] - overall_logp);

  ddef->mocc[L-1] = 1. - njcp;

  i++;
  njcp = 0.0;
  /*N state */
  njcp += expf(fwd->xmx[(L-3)*p7G_NXCELLS+p7G_N] + bck->xmx[L*p7G_NXCELLS+p7G_N] + gm_fs5->xsc[p7P_N][p7P_LOOP] - overall_logp);

  /* C state */
  njcp += expf(fwd->xmx[(L-3)*p7G_NXCELLS+p7G_C] + bck->xmx[L*p7G_NXCELLS+p7G_C] + gm_fs5->xsc[p7P_C][p7P_LOOP] - overall_logp);

  /* J state */
  njcp += expf(fwd->xmx[(L-3)*p7G_NXCELLS+p7G_J] + bck->xmx[L*p7G_NXCELLS+p7G_J] + gm_fs5->xsc[p7P_J][p7P_LOOP] - overall_logp);

  ddef->mocc[L] = 1. - njcp;

  ddef->L = L;
  return eslOK;

}

/*------------------ end, decoding algorithms -------------------*/



/*****************************************************************
 * 2. Benchmark driver
 *****************************************************************/
#ifdef p7GENERIC_DECODING_FRAMESHIFT_BENCHMARK
/*
   gcc -g -O2      -o generic_decoding_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_DECODING_FRAMESHIFT_BENCHMARK generic_decoding_frameshift.c -lhmmer -leasel -lm

   icc -O3 -static -o generic_decoding_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_DECODING_FRAMESHIFT_BENCHMARK generic_decoding_frameshift.c -lhmmer -leasel -lm
   ./generic_benchmark_decoding_frameshift <hmmfile>
   
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
  { "-N",        eslARG_INT,   "5000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for posterior residue decoding, generic version";

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
  P7_FS_PROFILE  *gm_fs5   = NULL;
  P7_GMX         *fwd     = NULL;
  P7_GMX         *bck     = NULL;
  P7_GMX         *pp      = NULL;
  P7_IVX         *iv      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           fsc, bsc;
  double          Mcs;

  p7_FLogsumInit();

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  abcAA  = esl_alphabet_Create(eslAMINO);
  abcDNA = esl_alphabet_Create(eslDNA); 
  bgAA   = p7_bg_Create(abcAA);                  p7_bg_SetLength(bgAA, L/3);
  gcode  = esl_gencode_Create(abcDNA,abcAA);
  gm_fs5 = p7_profile_fs_Create(hmm->M, abcAA, p7P_5CODONS);  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_LOCAL);
  fwd    = p7_gmx_Create(gm_fs5->M, L, L, p7G_NSCELLS_FS);  
  bck    = p7_gmx_Create(gm_fs5->M, L, L, p7G_NSCELLS);
  pp     = p7_gmx_Create(gm_fs5->M, L, L, p7G_NSCELLS_FS);
  iv     = p7_ivx_Create(gm_fs5->M, p7P_5CODONS);
  bgDNA  = p7_bg_Create(abcDNA);                p7_bg_SetLength(bgDNA, L);

  esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
  p7_GForward_Frameshift(dsq, L, gm_fs5, fwd, iv, &fsc);
  p7_GBackward_Frameshift(dsq, L, gm_fs5, bck, iv, &bsc);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) 
    p7_GDecoding_Frameshift(gm_fs5, fwd, bck, pp);
  esl_stopwatch_Stop(w);

  Mcs  = (double) N * (double) L * (double) gm_fs5->M * 1e-6 / w->user;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n", gm_fs5->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_gmx_Destroy(pp);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(bck);
  p7_ivx_Destroy(iv);
  p7_profile_fs_Destroy(gm_fs5);
  p7_bg_Destroy(bgDNA);
  p7_bg_Destroy(bgAA);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcAA);
  esl_alphabet_Destroy(abcDNA);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_DECODING_FRAMESHIFT_BENCHMARK*/
/*------------------ end, benchmark driver ----------------------*/




