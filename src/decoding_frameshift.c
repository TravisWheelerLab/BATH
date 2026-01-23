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
 *            aware codon profile  <gm_fs> and filled Forward and 
 *            Backward matrices <fwd>, <bck> for the profile aligned 
 *            to that target sequence. The nucleotide posterior decoding 
 *            is stored in a DP matrix <pp>, provided by the caller, and 
 *            the codon posterior decoding is overwriten to the <fwd>
 *            matrix.
 *            
 *            Each residue <i> must have been emitted by match state
 *            <1..M>, insert state <1..M-1>, or an NN, CC, or JJ loop
 *            transition.  For <dp = pp->dp>, <xmx = pp->xmx>,
 *            <MMX_FS(i,k,p7G_C0)> is the probability that match <k> 
 *            emitted the residue <i> as the last nucletotide in any 
 *            (quasi)codon while <MMX_FS(i,k,p7G_C1-5)> gives the same 
 *            probability but for each scepific (quasi)codon length; 
 *            <IMX(i,k)> is the probability that insert
 *            <k> emitted <i> as the last nucleotide in a cdodon ending 
 *            in residue <i>; <XMX(i,N)>,<XMX(i,C)>, <XMX(i,J)> are the 
 *            probabilities that residue <i> was the last nucleotide in 
 *            a codon emitted on an N->N, C->C, or J->J transition. 
 *            The frameshift aware implementation requires two seperate 
 *            normalizations. In the <pp> matrix we normalize so that 
 *            the sum over all emssions that include the residue <i> at
 *            any position in in the (quasi)codon. In the <fwd> matrix we 
 *            normalize so that the sum over all emissions which end at 
 *            <i>.The <pp> posteriors are used for alignment creation
 *            and bias scoring, and the <fwd> posteriors are used for 
 *            the alignemtn post-porb display line. 
 *
 * Args:      gm_fs - frameshift codon profile 
 *            fwd   - filled Forward matrix 
 *            bck   - filled Backward matrix
 *            pp    - RESULT: posterior decoding matrix.
 *
 * Returns:   <eslOK> on success.
 *
 *    
 */
int
p7_Decoding_Frameshift(const P7_FS_PROFILE *gm_fs, const P7_GMX *fwd, P7_GMX *bck, P7_GMX *pp)
{
  float      **dp   = pp->dp;
  float       *xmx  = pp->xmx;
  int          L    = fwd->L;
  int          M    = gm_fs->M;
  int          i,k;
  float        overall_sc = fwd->xmx[p7G_NXCELLS*L + p7G_C] + gm_fs->xsc[p7P_C][p7P_MOVE];
  float        denom, bias_denom;
  float        back_sc; 
  
  pp->M = M;
  pp->L = L;

  XMX_FS(0, p7G_N) = -eslINFINITY;
  
  XMX_FS(0, p7G_B) = -eslINFINITY;
  XMX_FS(0, p7G_E) = -eslINFINITY;
  XMX_FS(0, p7G_J) = -eslINFINITY;    
  XMX_FS(0, p7G_C) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_FS(0,k, p7G_C0) = MMX_FS(0,k, p7G_C1) = MMX_FS(0,k, p7G_C2) = MMX_FS(0,k, p7G_C3) = 
    MMX_FS(0,k, p7G_C4) = MMX_FS(0,k, p7G_C5) = IMX_FS(0,k) = DMX_FS(0,k) = -eslINFINITY;
 
  for (i = 1; i <= L; i++)
    {

     /* 0th position in model */    
      MMX_FS(i,0, p7G_C0) = MMX_FS(i,0, p7G_C1) = MMX_FS(i,0, p7G_C2) = MMX_FS(i,0, p7G_C3) = 
      MMX_FS(i,0, p7G_C4) = MMX_FS(i,0, p7G_C5) = IMX_FS(i,0) = DMX_FS(i,0) = -eslINFINITY;

      for (k = 1; k < M; k++)
      {
        /* probabilty from backward matrix */
        back_sc = bck->dp[i][k*p7G_NSCELLS + p7G_M] - overall_sc;  
         /* probability fom all 5 codons and their sum in the forward matrix */
         MMX_FS(i,k, p7G_C0) = fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C0] + back_sc;    
         MMX_FS(i,k, p7G_C1) = fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1] + back_sc;
         MMX_FS(i,k, p7G_C2) = fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2] + back_sc;
         MMX_FS(i,k, p7G_C3) = fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3] + back_sc;
         MMX_FS(i,k, p7G_C4) = fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4] + back_sc;
         MMX_FS(i,k, p7G_C5) = fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C5] + back_sc; 
         
         /* insert state probablity */
         IMX_FS(i,k) = fwd->dp[i][k*p7G_NSCELLS_FS + p7G_I] + bck->dp[i][k*p7G_NSCELLS + p7G_I] - overall_sc;

         /* no emision from delete state */
         DMX_FS(i,k) = -eslINFINITY;
      }

      /* final model position - no insert state */
      back_sc = bck->dp[i][M*p7G_NSCELLS + p7G_M] - overall_sc;
      
      MMX_FS(i,M, p7G_C0) = fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C0] + back_sc; 
      MMX_FS(i,M, p7G_C1) = fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C1] + back_sc; 
      MMX_FS(i,M, p7G_C2) = fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C2] + back_sc;
      MMX_FS(i,M, p7G_C3) = fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C3] + back_sc; 
      MMX_FS(i,M, p7G_C4) = fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C4] + back_sc; 
      MMX_FS(i,M, p7G_C5) = fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C5] + back_sc; 
      
      IMX_FS(i,M) = -eslINFINITY; 
      DMX_FS(i,M) = -eslINFINITY;

      /* no emission from E or B states */
      XMX_FS(i,p7G_E)   = -eslINFINITY;
      XMX(i,p7G_B)     = -eslINFINITY;  

      /* probaility from N, J and C states */
       if(i > 2)
       {
         XMX_FS(i,p7G_N) = fwd->xmx[p7G_NXCELLS*(i-3) + p7G_N] +  gm_fs->xsc[p7P_N][p7P_LOOP] +
                                bck->xmx[p7G_NXCELLS*i + p7G_N]     -  overall_sc;
         XMX_FS(i,p7G_C) = fwd->xmx[p7G_NXCELLS*(i-3) + p7G_C] +  gm_fs->xsc[p7P_C][p7P_LOOP] +
                                bck->xmx[p7G_NXCELLS*i + p7G_C]     -  overall_sc;
         
         XMX_FS(i,p7G_J) = fwd->xmx[p7G_NXCELLS*(i-3) + p7G_J] +  gm_fs->xsc[p7P_J][p7P_LOOP] +
                                bck->xmx[p7G_NXCELLS*i + p7G_J]     -  overall_sc;
       }
       else
       {
         XMX_FS(i,p7G_N) = bck->xmx[p7G_NXCELLS*i + p7G_N]     -  overall_sc;
         XMX_FS(i,p7G_C) = -eslINFINITY;
         XMX_FS(i,p7G_J) = -eslINFINITY;
        }
    }

  /* denom normailzes i for all codons in which i may be present 
   * while bias_denom only normailizes i for all codons ending in i*/
   for (i = 1; i <= L ; i++) {
    
    denom = -eslINFINITY;
    for (k = 1; k < M; k++) {  
      denom = p7_FLogsum(MMX_FS(i,k,p7G_C0), denom); 
      denom = p7_FLogsum(IMX_FS(i,k), denom);
    }
    denom = p7_FLogsum(MMX_FS(i,M,p7G_C0), denom);       
    
    denom = p7_FLogsum(XMX_FS(i,p7G_N), denom);
    denom = p7_FLogsum(XMX_FS(i,p7G_J), denom);
    denom = p7_FLogsum(XMX_FS(i,p7G_C), denom);

    bias_denom = -1*denom;

    for (k = 1; k < M; k++) {
      if(i < L)
      {
        denom = p7_FLogsum(MMX_FS(i+1,k,p7G_C5), denom);
        denom = p7_FLogsum(MMX_FS(i+1,k,p7G_C4), denom);
        denom = p7_FLogsum(MMX_FS(i+1,k,p7G_C3), denom);
        denom = p7_FLogsum(MMX_FS(i+1,k,p7G_C2), denom);
        denom = p7_FLogsum(IMX_FS(i+1,k) , denom);
      }
 
      if(i < L-1)
      {
        denom = p7_FLogsum(MMX_FS(i+2,k,p7G_C5), denom);
        denom = p7_FLogsum(MMX_FS(i+2,k,p7G_C4), denom);
        denom = p7_FLogsum(MMX_FS(i+2,k,p7G_C3), denom);
        denom = p7_FLogsum(IMX_FS(i+2,k), denom);
      } 

      if(i < L-2)
      {
        denom = p7_FLogsum(MMX_FS(i+3,k,p7G_C4), denom);
        denom = p7_FLogsum(MMX_FS(i+3,k,p7G_C5), denom);
      }    

      if(i < L-3)
        denom = p7_FLogsum(MMX_FS(i+4,k,p7G_C5), denom);
    }

    if(i < L)
    {
      denom = p7_FLogsum(MMX_FS(i+1,M,p7G_C5), denom);
      denom = p7_FLogsum(MMX_FS(i+1,M,p7G_C4), denom);
      denom = p7_FLogsum(MMX_FS(i+1,M,p7G_C3), denom);
      denom = p7_FLogsum(MMX_FS(i+1,M,p7G_C2), denom);
      denom = p7_FLogsum(XMX_FS(i+1,p7G_N), denom);
      denom = p7_FLogsum(XMX_FS(i+1,p7G_J), denom);
      denom = p7_FLogsum(XMX_FS(i+1,p7G_C), denom);
    }
 
    if(i < L-1)
    {
      denom = p7_FLogsum(MMX_FS(i+2,M,p7G_C5), denom);
      denom = p7_FLogsum(MMX_FS(i+2,M,p7G_C4), denom);
      denom = p7_FLogsum(MMX_FS(i+2,M,p7G_C3), denom);
      denom = p7_FLogsum(XMX_FS(i+2,p7G_N), denom);
      denom = p7_FLogsum(XMX_FS(i+2,p7G_J), denom);
      denom = p7_FLogsum(XMX_FS(i+2,p7G_C), denom);
    }
 
    if(i < L-2)
    {
      denom = p7_FLogsum(MMX_FS(i+3,M,p7G_C4), denom);
      denom = p7_FLogsum(MMX_FS(i+3,M,p7G_C5), denom);
    }

    if(i < L-3)
      denom = p7_FLogsum(MMX_FS(i+4,M,p7G_C5), denom);
    
    denom = -1*denom;
 
    for (k = 1; k < M; k++) {  
      fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C0] = expf(MMX_FS(i,k,p7G_C0) + bias_denom);    
      fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1] = expf(MMX_FS(i,k,p7G_C1) + bias_denom);
      fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2] = expf(MMX_FS(i,k,p7G_C2) + bias_denom);
      fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3] = expf(MMX_FS(i,k,p7G_C3) + bias_denom);
      fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4] = expf(MMX_FS(i,k,p7G_C4) + bias_denom);
      fwd->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C5] = expf(MMX_FS(i,k,p7G_C5) + bias_denom);
      fwd->dp[i][k*p7G_NSCELLS_FS + p7G_I]          = expf(IMX_FS(i,k)        + bias_denom);
      MMX_FS(i,k,p7G_C1) = MMX_FS(i,k,p7G_C1) + denom; 
      MMX_FS(i,k,p7G_C2) = MMX_FS(i,k,p7G_C2) + denom; 
      MMX_FS(i,k,p7G_C3) = MMX_FS(i,k,p7G_C3) + denom; 
      MMX_FS(i,k,p7G_C4) = MMX_FS(i,k,p7G_C4) + denom; 
      MMX_FS(i,k,p7G_C5) = MMX_FS(i,k,p7G_C5) + denom; 
      MMX_FS(i,k,p7G_C0) = MMX_FS(i,k,p7G_C0) + denom; 
      IMX_FS(i,k)        = IMX_FS(i,k)        + denom;
    }

    fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C0] = expf(MMX_FS(i,M,p7G_C0) + bias_denom); 
    fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C1] = expf(MMX_FS(i,M,p7G_C1) + bias_denom);
    fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C2] = expf(MMX_FS(i,M,p7G_C2) + bias_denom);
    fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C3] = expf(MMX_FS(i,M,p7G_C3) + bias_denom);
    fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C4] = expf(MMX_FS(i,M,p7G_C4) + bias_denom);
    fwd->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C5] = expf(MMX_FS(i,M,p7G_C5) + bias_denom);
    fwd->xmx[p7G_NXCELLS*i + p7G_N]               = expf(XMX_FS(i,p7G_N)    + bias_denom);
    fwd->xmx[p7G_NXCELLS*i + p7G_J]               = expf(XMX_FS(i,p7G_J)    + bias_denom);
    fwd->xmx[p7G_NXCELLS*i + p7G_C]               = expf(XMX_FS(i,p7G_C)    + bias_denom);

    MMX_FS(i,M,p7G_C1) =  MMX_FS(i,M,p7G_C1) + denom; 
    MMX_FS(i,M,p7G_C2) =  MMX_FS(i,M,p7G_C2) + denom; 
    MMX_FS(i,M,p7G_C3) =  MMX_FS(i,M,p7G_C3) + denom; 
    MMX_FS(i,M,p7G_C4) =  MMX_FS(i,M,p7G_C4) + denom; 
    MMX_FS(i,M,p7G_C5) =  MMX_FS(i,M,p7G_C5) + denom; 
    MMX_FS(i,M,p7G_C0) =  MMX_FS(i,M,p7G_C0) + denom; 
    XMX_FS(i,p7G_N) = XMX_FS(i,p7G_N) + denom;
    XMX_FS(i,p7G_J) = XMX_FS(i,p7G_J) + denom;
    XMX_FS(i,p7G_C) = XMX_FS(i,p7G_C) + denom;
 
  } 

  return eslOK;
}

/* Function:  p7_DomainDecoding_Frameshift
 * Synopsis:  Posterior decoding of domain location.
 *
 * Purpose:   The caller has calculated Forward and Backward matrices
 *            <fwd> and <bck> for frameshift codon model <gm_fs> 
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
 * Args:      gm_fs - profile
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
p7_DomainDecoding_Frameshift(const P7_FS_PROFILE *gm_fs, const P7_GMX *fwd, const P7_GMX *bck, P7_DOMAINDEF *ddef)
{
  int   L            = fwd->L;
  float overall_logp = p7_FLogsum( fwd->xmx[(L)*p7G_NXCELLS+p7G_C],
                       p7_FLogsum( fwd->xmx[(L-1)*p7G_NXCELLS+p7G_C],
                                   fwd->xmx[(L-2)*p7G_NXCELLS+p7G_C])) +
                                   gm_fs->xsc[p7P_C][p7P_MOVE];
  float njcp; 
  int   i;

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

    /* Sum poropabilities in the N, J and C states */
    njcp = 0.0;
    njcp += expf(fwd->xmx[(i-3)*p7G_NXCELLS+p7G_N] + bck->xmx[i*p7G_NXCELLS+p7G_N] + gm_fs->xsc[p7P_N][p7P_LOOP] - overall_logp);
    njcp += expf(fwd->xmx[(i-3)*p7G_NXCELLS+p7G_J] + bck->xmx[i*p7G_NXCELLS+p7G_J] + gm_fs->xsc[p7P_J][p7P_LOOP] - overall_logp);
    njcp += expf(fwd->xmx[(i-3)*p7G_NXCELLS+p7G_C] + bck->xmx[i*p7G_NXCELLS+p7G_C] + gm_fs->xsc[p7P_C][p7P_LOOP] - overall_logp);
    
    njcp += expf(fwd->xmx[(i-2)*p7G_NXCELLS+p7G_N] + bck->xmx[(i+1)*p7G_NXCELLS+p7G_N] + gm_fs->xsc[p7P_N][p7P_LOOP] - overall_logp);  
    njcp += expf(fwd->xmx[(i-2)*p7G_NXCELLS+p7G_J] + bck->xmx[(i+1)*p7G_NXCELLS+p7G_J] + gm_fs->xsc[p7P_J][p7P_LOOP] - overall_logp);
    njcp += expf(fwd->xmx[(i-2)*p7G_NXCELLS+p7G_C] + bck->xmx[(i+1)*p7G_NXCELLS+p7G_C] + gm_fs->xsc[p7P_C][p7P_LOOP] - overall_logp);
  
    njcp += expf(fwd->xmx[(i-1)*p7G_NXCELLS+p7G_N] + bck->xmx[(i+2)*p7G_NXCELLS+p7G_N] + gm_fs->xsc[p7P_N][p7P_LOOP] - overall_logp);
    njcp += expf(fwd->xmx[(i-1)*p7G_NXCELLS+p7G_J] + bck->xmx[(i+2)*p7G_NXCELLS+p7G_J] + gm_fs->xsc[p7P_J][p7P_LOOP] - overall_logp);
    njcp += expf(fwd->xmx[(i-1)*p7G_NXCELLS+p7G_C] + bck->xmx[(i+2)*p7G_NXCELLS+p7G_C] + gm_fs->xsc[p7P_C][p7P_LOOP] - overall_logp);
    
    /* Probability of i emitted by the core model is aaprixmated as 1 - probability of i emited by the specials */
    ddef->mocc[i] = 1. - njcp;
  }
  njcp = 0.0;
  njcp += expf(fwd->xmx[(i-3)*p7G_NXCELLS+p7G_N] + bck->xmx[i*p7G_NXCELLS+p7G_N] + gm_fs->xsc[p7P_N][p7P_LOOP] - overall_logp);
  njcp += expf(fwd->xmx[(i-3)*p7G_NXCELLS+p7G_J] + bck->xmx[i*p7G_NXCELLS+p7G_J] + gm_fs->xsc[p7P_J][p7P_LOOP] - overall_logp);
  njcp += expf(fwd->xmx[(i-3)*p7G_NXCELLS+p7G_C] + bck->xmx[i*p7G_NXCELLS+p7G_C] + gm_fs->xsc[p7P_C][p7P_LOOP] - overall_logp);

  njcp += expf(fwd->xmx[(i-2)*p7G_NXCELLS+p7G_N] + bck->xmx[(i+1)*p7G_NXCELLS+p7G_N] + gm_fs->xsc[p7P_N][p7P_LOOP] - overall_logp);
  njcp += expf(fwd->xmx[(i-2)*p7G_NXCELLS+p7G_J] + bck->xmx[(i+1)*p7G_NXCELLS+p7G_J] + gm_fs->xsc[p7P_J][p7P_LOOP] - overall_logp);
  njcp += expf(fwd->xmx[(i-2)*p7G_NXCELLS+p7G_C] + bck->xmx[(i+1)*p7G_NXCELLS+p7G_C] + gm_fs->xsc[p7P_C][p7P_LOOP] - overall_logp);
  ddef->mocc[L-1] = 1. - njcp;

  njcp = 0.0;
  njcp += expf(fwd->xmx[(i-3)*p7G_NXCELLS+p7G_N] + bck->xmx[i*p7G_NXCELLS+p7G_N] + gm_fs->xsc[p7P_N][p7P_LOOP] - overall_logp);
  njcp += expf(fwd->xmx[(i-3)*p7G_NXCELLS+p7G_J] + bck->xmx[i*p7G_NXCELLS+p7G_J] + gm_fs->xsc[p7P_J][p7P_LOOP] - overall_logp);
  njcp += expf(fwd->xmx[(i-3)*p7G_NXCELLS+p7G_C] + bck->xmx[i*p7G_NXCELLS+p7G_C] + gm_fs->xsc[p7P_C][p7P_LOOP] - overall_logp); 
  ddef->mocc[L] = 1. - njcp;

  ddef->L = L;
  return eslOK;
}

/*------------------ end, decoding algorithms -------------------*/



/*****************************************************************
 * 2. Benchmark driver
 *****************************************************************/
#ifdef p7DECODING_FRAMESHIFT_BENCHMARK
/*
   gcc -g -O2      -o decoding_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7DECODING_FRAMESHIFT_BENCHMARK decoding_frameshift.c -lhmmer -leasel -lm

   icc -O3 -static -o decoding_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7DECODING_FRAMESHIFT_BENCHMARK decoding_frameshift.c -lhmmer -leasel -lm
   ./benchmark_decoding_frameshift <hmmfile>
   
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
  P7_FS_PROFILE  *gm_fs   = NULL;
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

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  abcAA  = esl_alphabet_Create(eslAMINO);
  abcDNA = esl_alphabet_Create(eslDNA); 
  bgAA  = p7_bg_Create(abcAA);                  p7_bg_SetLength(bgAA, L/3);
  gcode = esl_gencode_Create(abcDNA,abcAA);
  gm_fs = p7_profile_fs_Create(hmm->M, abcAA);  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs, L/3, p7_LOCAL);
  fwd   = p7_gmx_fs_Create(gm_fs->M, L, L, p7P_5CODONS);  
  bck   = p7_gmx_fs_Create(gm_fs->M, L, L, 0);
  pp    = p7_gmx_fs_Create(gm_fs->M, L, L, p7P_5CODONS);
  iv    = p7_ivx_Create(gm_fs->M, p7P_5CODONS);
  bgDNA = p7_bg_Create(abcDNA);

  esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
  p7_Forward_Frameshift(dsq, gcode, L, gm_fs, fwd, iv, &fsc);
  p7_Backward_Frameshift(dsq, gcode, L, gm_fs, bck, iv, &bsc);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) 
    p7_Decoding_Frameshift(gm_fs, fwd, bck, pp);
  esl_stopwatch_Stop(w);

  Mcs  = (double) N * (double) L * (double) gm_fs->M * 1e-6 / w->user;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n", gm_fs->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_gmx_Destroy(pp);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(bck);
  p7_ivx_Destroy(iv);
  p7_profile_fs_Destroy(gm_fs);
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
#endif /*p7DECODING_FRAMESHIFT_BENCHMARK*/
/*------------------ end, benchmark driver ----------------------*/




