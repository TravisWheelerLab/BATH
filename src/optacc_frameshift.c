/* Optimal accuracy alignment; frameshift version.
 * 
 * Contents:
 *   1. Optimal alignment accuracy fill.
 *   2. Optimal alignment accuracy traceback.
 *   3. Benchmark driver
 *   4. Unit tests
 *   5. Test driver
 * 
 */
#include "p7_config.h"

#include <float.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

/*****************************************************************
 * 1. Optimal alignment fill and traceback.
 *****************************************************************/

#define TSCDELTA(s,k) ( (tsc[(k) * p7P_NTRANS + (s)] == -eslINFINITY) ? -eslINFINITY : 0.0)
#define MVX(i,k,c) (max_val[(k)*p7P_FULL_CODONS+L3-(i)+(c)])
enum p7e_optacc_etrace {
  p7E_C1 = 0,
  p7E_C2 = 1,
  p7E_C3 = 2,
  p7E_C4 = 3,
  p7E_C5 = 4,
  p7E_D  = 5,
  p7E_E  = 6,
};
#define p7E_EXIT 7
#define p7E_MATCH 5

/* The TSCDELTA is used to make impossible paths impossible in the
 * optimal accuracy decoding algorithm;
 */


/* Function:  p7_OptimalAccuracy_Frameshift()
 * Synopsis:  Optimal accuracy decoding: fill. 
 *
 * Purpose:   Calculates the fill step of the optimal accuracy decoding
 *            algorithm.
 *            
 *            Caller provides the posterior decoding matrix <pp>,
 *            which was calculated by Forward/Backward on a target sequence
 *            of length <L> using the query model <gmi_fs>.
 *            
 *            Caller also provides a DP matrix <gx>, allocated for the
 *            <gm->M> by <pp->L> comparison. The routine fills this in
 *            with OA scores.
 *
 *            For the frameshift version, porobabilities are smaller due to
 *            being normaalized across multile frames, so they are kept in 
 *            log space 
 *            
 * Args:      gm_fs - query profile      
 *            pp    - posterior decoding matrix created by <p7_GPosteriorDecoding()>
 *            gx    - RESULT: caller provided DP matrix for <gm->M> by <L> 
 *            ret_e - RETURN: expected number of correctly decoded positions 
 *
 * Returns:   <eslOK> on success, and <*ret_e> contains the final OA
 *            score, which is the expected number of correctly decoded
 *            positions in the target sequence (up to <L>).
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_OptimalAccuracy_Frameshift(const P7_FS_PROFILE *gm_fs, const P7_GMX *pp, P7_GMX *gx, float *ret_e)
{
  int          L    = pp->L;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float const *tsc  = gm_fs->tsc;
  int          i,k;
  int          M    = gm_fs->M;
  float        esc  = p7_fs_profile_IsLocal(gm_fs) ? 1.0 : 0.0;
  float        nn, jj, cc;
  float        nb, jb, ej, ec;
  float        max1, max2, max3, max4, max5;

  nn = ((gm_fs->xsc[p7P_N][p7P_LOOP] == -eslINFINITY) ? -eslINFINITY : 0.0);
  jj = ((gm_fs->xsc[p7P_J][p7P_LOOP] == -eslINFINITY) ? -eslINFINITY : 0.0);
  cc = ((gm_fs->xsc[p7P_C][p7P_LOOP] == -eslINFINITY) ? -eslINFINITY : 0.0);
  nb = ((gm_fs->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? -eslINFINITY : 0.0);
  jb = ((gm_fs->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? -eslINFINITY : 0.0);
  ej = ((gm_fs->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) ? -eslINFINITY : 0.0); 
  ec = ((gm_fs->xsc[p7P_E][p7P_MOVE] == -eslINFINITY) ? -eslINFINITY : 0.0);
 
  /* Initialization of the zero row (i=0; no residues to account for.  */
  XMX(0,p7G_N) = 0.;                                                /* S->N, p=1            */
  XMX(0,p7G_B) = 0.;                                                /* S->N->B, no N-tail   */
  XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) = -eslINFINITY;        /* need seq to get here */
  for (k = 0; k <= M; k++)
    MMX(0,k) = IMX(0,k) = DMX(0,k) = -eslINFINITY;                  /* need seq to get here */

  /* Initialization of row 1 */
  MMX(1,0) = IMX(1,0) = DMX(1,0) = XMX(1,p7G_E) = -eslINFINITY;
  for (k = 1; k < M; k++) {
    MMX(1,k)  = TSCDELTA(p7P_BM, k-1) + pp->dp[1][k*p7G_NSCELLS_FS + p7G_M + p7G_C1];
    IMX(1,k)  = -eslINFINITY;  
    DMX(1,k)  = ESL_MAX( TSCDELTA(p7P_MD, k-1) + MMX(1,k-1),
                         TSCDELTA(p7P_DD, k-1) + DMX(1,k-1));

    XMX(1,p7G_E) = ESL_MAX(XMX(1,p7G_E), esc * MMX(1,k));
  }
  
  MMX(1,M)  = TSCDELTA(p7P_BM, M-1) + pp->dp[1][M*p7G_NSCELLS_FS + p7G_M + p7G_C1];
  IMX(1,M)  = -eslINFINITY;
  DMX(1,M)  = ESL_MAX( TSCDELTA(p7P_MD, M-1) + MMX(1,M-1),
                       TSCDELTA(p7P_DD, M-1) + DMX(1,M-1));

  XMX(1,p7G_E) = ESL_MAX(XMX(1,p7G_E), ESL_MAX(MMX(1,M), DMX(1, M)));

  XMX(1,p7G_J) = ej + XMX(1,p7G_E);
  XMX(1,p7G_C) = ec + XMX(1,p7G_E);
  XMX(1,p7G_N) = nn + pp->xmx[1*p7G_NXCELLS + p7G_N];
  XMX(1,p7G_B) = ESL_MAX(nb + XMX(1,p7G_N), jb + XMX(1,p7G_J));
  
  /* Initialization of row 2 */
  MMX(2,0) = IMX(2,0) = DMX(2,0) = XMX(2,p7G_E) = -eslINFINITY;  
  for (k = 1; k < M; k++) {
    max1 = ESL_MAX(TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(1,k-1),   pp->dp[2][k*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
           ESL_MAX(TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(1,k-1),   pp->dp[2][k*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
                   TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(1,p7G_B), pp->dp[2][k*p7G_NSCELLS_FS + p7G_M + p7G_C1])));

    MMX(2,k)  = ESL_MAX( max1, TSCDELTA(p7P_BM, k-1) + pp->dp[2][k*p7G_NSCELLS_FS + p7G_M + p7G_C2]);
    IMX(2,k)  = -eslINFINITY;  
    DMX(2,k)  = ESL_MAX( TSCDELTA(p7P_MD, k-1) + MMX(2,k-1),
                         TSCDELTA(p7P_DD, k-1) + DMX(2,k-1));
           
    XMX(2,p7G_E) = ESL_MAX(XMX(2,p7G_E), esc * MMX(2,k));
  }

  max1 = ESL_MAX(TSCDELTA(p7P_MM, M-1) + p7_FLogsum(MMX(2,M-1),   pp->dp[2][M*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
         ESL_MAX(TSCDELTA(p7P_DM, M-1) + p7_FLogsum(DMX(2,M-1),   pp->dp[2][M*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
                 TSCDELTA(p7P_BM, M-1) + p7_FLogsum(XMX(2,p7G_B), pp->dp[2][M*p7G_NSCELLS_FS + p7G_M + p7G_C1])));

  MMX(2,M)  = ESL_MAX( max1, TSCDELTA(p7P_BM, M-1) + pp->dp[2][M*p7G_NSCELLS_FS + p7G_M + p7G_C2]);
  IMX(2,M)  = -eslINFINITY;
  DMX(2,M)  = ESL_MAX( TSCDELTA(p7P_MD, M-1) + MMX(2,M-1),
                       TSCDELTA(p7P_DD, M-1) + DMX(2,M-1));

  XMX(2,p7G_E) = ESL_MAX(XMX(2,p7G_E), ESL_MAX(MMX(2,M), DMX(2, M)));

  XMX(2,p7G_J) = ej + XMX(2,p7G_E);
  XMX(2,p7G_C) = ec + XMX(2,p7G_E);
  XMX(2,p7G_N) = nn + pp->xmx[2*p7G_NXCELLS + p7G_N];
  XMX(2,p7G_B) = ESL_MAX(nb + XMX(2,p7G_N), jb + XMX(2,p7G_J));

  /* Initialization of row 3-4 */
  for (i = 3; i < 5; i++)
  {
    MMX(i,0) = IMX(i,0) = DMX(i,0) = XMX(i,p7G_E) = -eslINFINITY;
    for (k = 1; k < M; k++) {
      max1 = ESL_MAX( TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(i-1,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
             ESL_MAX( TSCDELTA(p7P_IM, k-1) + p7_FLogsum(IMX(i-1,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
             ESL_MAX( TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(i-1,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
                      TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(i-1,p7G_B), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1]))));   

      max2 = ESL_MAX( TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(i-2,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2]),
             ESL_MAX( TSCDELTA(p7P_IM, k-1) + p7_FLogsum(IMX(i-2,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2]),
             ESL_MAX( TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(i-2,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2]),
                      TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(i-2,p7G_B), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2]))));

      max3 = ESL_MAX( TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(i-3,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3]),
             ESL_MAX( TSCDELTA(p7P_IM, k-1) + p7_FLogsum(IMX(i-3,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3]),
             ESL_MAX( TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(i-3,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3]),
                      TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(i-3,p7G_B), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3]))));
   
      if(i == 4) 
        max4 = ESL_MAX( TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(i-4,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4]),
               ESL_MAX( TSCDELTA(p7P_IM, k-1) + p7_FLogsum(IMX(i-4,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4]),
               ESL_MAX( TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(i-4,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4]),
                        TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(i-4,p7G_B), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4]))));
      else max4 = -eslINFINITY;

      MMX(i,k)  = ESL_MAX( max1, ESL_MAX( max2, ESL_MAX( max3, max4)));
      IMX(i,k)  = ESL_MAX( TSCDELTA(p7P_MI, k) + p7_FLogsum(MMX(i-3,k), pp->dp[i][k*p7G_NSCELLS_FS + p7G_I]),
                           TSCDELTA(p7P_II, k) + p7_FLogsum(IMX(i-3,k), pp->dp[i][k*p7G_NSCELLS_FS + p7G_I]));    
      DMX(i,k)  = ESL_MAX( TSCDELTA(p7P_MD, k-1) + MMX(i,k-1),
                           TSCDELTA(p7P_DD, k-1) + DMX(i,k-1)); 

      XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), esc * MMX(i,k));    
    }

    max1 = ESL_MAX( TSCDELTA(p7P_MM, M-1) + p7_FLogsum(MMX(i-1,M-1),   pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
           ESL_MAX( TSCDELTA(p7P_IM, M-1) + p7_FLogsum(IMX(i-1,M-1),   pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
           ESL_MAX( TSCDELTA(p7P_DM, M-1) + p7_FLogsum(DMX(i-1,M-1),   pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
                    TSCDELTA(p7P_BM, M-1) + p7_FLogsum(XMX(i-1,p7G_B), pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C1]))));

    max2 = ESL_MAX( TSCDELTA(p7P_MM, M-1) + p7_FLogsum(MMX(i-2,M-1),   pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C2]),
           ESL_MAX( TSCDELTA(p7P_IM, M-1) + p7_FLogsum(IMX(i-2,M-1),   pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C2]),
           ESL_MAX( TSCDELTA(p7P_DM, M-1) + p7_FLogsum(DMX(i-2,M-1),   pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C2]),
                    TSCDELTA(p7P_BM, M-1) + p7_FLogsum(XMX(i-2,p7G_B), pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C2]))));

    max3 = ESL_MAX( TSCDELTA(p7P_MM, M-1) + p7_FLogsum(MMX(i-3,M-1),   pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C3]),
           ESL_MAX( TSCDELTA(p7P_IM, M-1) + p7_FLogsum(IMX(i-3,M-1),   pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C3]),
           ESL_MAX( TSCDELTA(p7P_DM, M-1) + p7_FLogsum(DMX(i-3,M-1),   pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C3]),
                    TSCDELTA(p7P_BM, M-1) + p7_FLogsum(XMX(i-3,p7G_B), pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C3]))));

    if(i == 4)
      max4 = ESL_MAX( TSCDELTA(p7P_MM, M-1) + p7_FLogsum(MMX(i-4,M-1),   pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C4]),
             ESL_MAX( TSCDELTA(p7P_IM, M-1) + p7_FLogsum(IMX(i-4,M-1),   pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C4]),
             ESL_MAX( TSCDELTA(p7P_DM, M-1) + p7_FLogsum(DMX(i-4,M-1),   pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C4]),
                      TSCDELTA(p7P_BM, M-1) + p7_FLogsum(XMX(i-4,p7G_B), pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C4]))));
    else max4 = -eslINFINITY;   
  
    MMX(i,M)  = ESL_MAX( max1, ESL_MAX( max2, ESL_MAX( max3, max4)));
    IMX(i,M)  = -eslINFINITY; 
    DMX(i,M)  = ESL_MAX( TSCDELTA(p7P_MD, M-1) + MMX(i,M-1),
                         TSCDELTA(p7P_DD, M-1) + DMX(i,M-1));

    XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), ESL_MAX(MMX(i,M), DMX(i, M))); 
        
    XMX(i, p7G_J) = ESL_MAX( jj + p7_FLogsum(XMX(i-3,p7G_J), pp->xmx[i*p7G_NXCELLS + p7G_J]),
                             ej + XMX(i,  p7G_E));
  
    XMX(i,p7G_C) = ESL_MAX( cc + p7_FLogsum(XMX(i-3,p7G_C),  pp->xmx[i*p7G_NXCELLS + p7G_C]),
                            ec + XMX(i,  p7G_E));
  
    XMX(i,p7G_N) = nn +  p7_FLogsum(XMX(i-3,p7G_N), pp->xmx[i*p7G_NXCELLS + p7G_N]);  
     
    XMX(i,p7G_B) = ESL_MAX( nb + XMX(i,  p7G_N), jb + XMX(i,  p7G_J)); 

  }
 
  for (i = 5; i <= L; i++) {
    MMX(i,0) = IMX(i,0) = DMX(i,0) = XMX(i,p7G_E) = -eslINFINITY;
   
    for (k = 1; k < M; k++) {

      max1 = ESL_MAX( TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(i-1,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
             ESL_MAX( TSCDELTA(p7P_IM, k-1) + p7_FLogsum(IMX(i-1,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
             ESL_MAX( TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(i-1,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
                      TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(i-1,p7G_B), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1]))));

      max2 = ESL_MAX( TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(i-2,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2]),
             ESL_MAX( TSCDELTA(p7P_IM, k-1) + p7_FLogsum(IMX(i-2,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2]),
             ESL_MAX( TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(i-2,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2]),
                      TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(i-2,p7G_B), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2]))));

      max3 = ESL_MAX( TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(i-3,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3]),
             ESL_MAX( TSCDELTA(p7P_IM, k-1) + p7_FLogsum(IMX(i-3,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3]),
             ESL_MAX( TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(i-3,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3]),
                      TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(i-3,p7G_B), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3]))));

      max4 = ESL_MAX( TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(i-4,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4]),
             ESL_MAX( TSCDELTA(p7P_IM, k-1) + p7_FLogsum(IMX(i-4,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4]),
             ESL_MAX( TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(i-4,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4]),
                      TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(i-4,p7G_B), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4]))));

      max5 = ESL_MAX( TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(i-5,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C5]),
             ESL_MAX( TSCDELTA(p7P_IM, k-1) + p7_FLogsum(IMX(i-5,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C5]),
             ESL_MAX( TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(i-5,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C5]),
                      TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(i-5,p7G_B), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C5]))));

      MMX(i,k)  = ESL_MAX( max1, ESL_MAX( max2, ESL_MAX( max3, ESL_MAX(max4, max5))));

      IMX(i,k)  = ESL_MAX( TSCDELTA(p7P_MI, k) + p7_FLogsum(MMX(i-3,k), pp->dp[i][k*p7G_NSCELLS_FS + p7G_I]),
                           TSCDELTA(p7P_II, k) + p7_FLogsum(IMX(i-3,k), pp->dp[i][k*p7G_NSCELLS_FS + p7G_I]));
      DMX(i,k)  = ESL_MAX( TSCDELTA(p7P_MD, k-1) + MMX(i,k-1), TSCDELTA(p7P_DD, k-1) + DMX(i,k-1));

      XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), esc * MMX(i,k));  
    } 

    max1 = ESL_MAX( TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(i-1,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
           ESL_MAX( TSCDELTA(p7P_IM, k-1) + p7_FLogsum(IMX(i-1,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
           ESL_MAX( TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(i-1,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1]),
                    TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(i-1,p7G_B), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1]))));

    max2 = ESL_MAX( TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(i-2,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2]),
           ESL_MAX( TSCDELTA(p7P_IM, k-1) + p7_FLogsum(IMX(i-2,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2]),
           ESL_MAX( TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(i-2,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2]),
                    TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(i-2,p7G_B), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2]))));

    max3 = ESL_MAX( TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(i-3,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3]),
           ESL_MAX( TSCDELTA(p7P_IM, k-1) + p7_FLogsum(IMX(i-3,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3]),
           ESL_MAX( TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(i-3,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3]),
                    TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(i-3,p7G_B), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3]))));

    max4 = ESL_MAX( TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(i-4,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4]),
           ESL_MAX( TSCDELTA(p7P_IM, k-1) + p7_FLogsum(IMX(i-4,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4]),
           ESL_MAX( TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(i-4,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4]),
                    TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(i-4,p7G_B), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4]))));

    max5 = ESL_MAX( TSCDELTA(p7P_MM, k-1) + p7_FLogsum(MMX(i-5,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C5]),
           ESL_MAX( TSCDELTA(p7P_IM, k-1) + p7_FLogsum(IMX(i-5,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C5]),
           ESL_MAX( TSCDELTA(p7P_DM, k-1) + p7_FLogsum(DMX(i-5,k-1),   pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C5]),
                    TSCDELTA(p7P_BM, k-1) + p7_FLogsum(XMX(i-5,p7G_B), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C5]))));
  
    
    MMX(i,M)     = ESL_MAX( max1, ESL_MAX( max2, ESL_MAX( max3, ESL_MAX( max4, max5))));
    IMX(i,M)     = -eslINFINITY;
    DMX(i,M)     = ESL_MAX( TSCDELTA(p7P_MD, M-1) + MMX(i,M-1), TSCDELTA(p7P_DD, M-1) + DMX(i,M-1));

    XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), ESL_MAX(MMX(i,M), DMX(i, M)));

    XMX(i, p7G_J) = ESL_MAX( jj + p7_FLogsum(XMX(i-3,p7G_J), pp->xmx[i*p7G_NXCELLS + p7G_J]),
                             ej + XMX(i,  p7G_E));

    XMX(i,p7G_C) = ESL_MAX( cc + p7_FLogsum(XMX(i-3,p7G_C),  pp->xmx[i*p7G_NXCELLS + p7G_C]),
                            ec + XMX(i,  p7G_E));

    XMX(i,p7G_N) = nn +  p7_FLogsum(XMX(i-3,p7G_N), pp->xmx[i*p7G_NXCELLS + p7G_N]);

    XMX(i,p7G_B) = ESL_MAX( nb + XMX(i,  p7G_N), jb + XMX(i,  p7G_J));
                            
  }
  *ret_e = p7_FLogsum( XMX(L    ,p7G_C),
           p7_FLogsum( XMX(L-1  ,p7G_C),
                       XMX(L-2  ,p7G_C)));

  return eslOK;
}

    
/*---------------------- end, oa fill ---------------------------*/

/*****************************************************************
 * 2. Optimal alignment accuracy, traceback
 *****************************************************************/

static inline float get_postprob(const P7_GMX *pp, int scur, int sprv, int k, int i);
static inline int select_m(const P7_FS_PROFILE *gm_fs,                   const P7_GMX *gx, int c, int i, int k);
static inline int select_d(const P7_FS_PROFILE *gm_fs,                   const P7_GMX *gx,        int i, int k);
static inline int select_i(const P7_FS_PROFILE *gm_fs,                   const P7_GMX *gx,        int i, int k);
static inline int select_n(int i);
static inline int select_c(const P7_FS_PROFILE *gm_fs, const P7_GMX *pp, const P7_GMX *gx,        int i);
static inline int select_j(const P7_FS_PROFILE *gm_fs, const P7_GMX *pp, const P7_GMX *gx,        int i);
static inline int select_e(const P7_FS_PROFILE *gm_fs,                   const P7_GMX *gx,        int i, int *ret_k);
static inline int select_b(const P7_FS_PROFILE *gm_fs,                   const P7_GMX *gx,        int i);


/* Function:  p7_GOATrace()
 * Synopsis:  Optimal accuracy decoding: traceback.
 * Incept:    SRE, Fri Feb 29 12:59:11 2008 [Janelia]
 *
 * Purpose:   The traceback stage of the optimal accuracy decoding algorithm
 *            \citep{Kall05}.
 *            
 *            Caller provides the OA DP matrix <gx> that was just
 *            calculated by <p7_GOptimalAccuracy()>, as well as the
 *            posterior decoding matrix <pp>, which was calculated by
 *            Forward/Backward on a target sequence of length <L>
 *            using the query model <gm>.
 *            
 *            Caller provides an empty traceback structure <tr> to
 *            hold the result, allocated to hold optional posterior
 *            probability annotation on residues (with
 *            <p7_trace_CreateWithPP()>, generally).  This will be
 *            internally reallocated as needed for larger traces.
 *
 * Args:      gm    - query profile      
 *            pp    - posterior decoding matrix created by <p7_PosteriorDecoding()>
 *            gx    - OA DP matrix calculated by  <p7_OptimalAccuracyDP()>
 *            tr    - RESULT: OA traceback, allocated with posterior probs
 *
 * Returns:   <eslOK> on success, and <tr> contains the OA traceback.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_OATrace_Frameshift(const P7_FS_PROFILE *gm_fs, const P7_GMX *pp, const P7_GMX *gx, const P7_GMX *probs, P7_TRACE *tr)
{
  int           i   = gx->L;  /* position in seq (1..L)         */
  int           k   = 0;  /* position in model (1..M)       */
  //int           arg = 0;
  ESL_DSQ       c = 0;
  float        postprob;
  int          sprv, scur;
  int          status;
  float match_codon[5];

#if eslDEBUGLEVEL > 0
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif
  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i, c)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i, c)) != eslOK) return status;

  sprv = p7T_C;
  while (sprv != p7T_S) 
    { 
     switch (sprv) {
      case p7T_M: scur = select_m(gm_fs,     gx, c, i,  k);          k--;  break;
      case p7T_D: scur = select_d(gm_fs,     gx,    i,  k);          k--;  break;
      case p7T_I: scur = select_i(gm_fs,     gx,    i,  k); i -= 3;        break;
      case p7T_N: scur = select_n(                  i);                    break;
      case p7T_C: scur = select_c(gm_fs, pp, gx,    i);                    break;
      case p7T_J: scur = select_j(gm_fs, pp, gx,    i);                    break;
      case p7T_E: scur = select_e(gm_fs,     gx,    i, &k);                break;
      case p7T_B: scur = select_b(gm_fs,     gx,    i);                    break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }
      if (scur == -1) ESL_EXCEPTION(eslEINVAL, "OA traceback choice failed");
      if(scur == p7T_M) 
      {
        match_codon[0] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1];
        match_codon[1] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2];
        match_codon[2] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3];
        match_codon[3] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4];
        match_codon[4] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C5]; 
       
        c = esl_vec_FArgMax(match_codon, 5) + 1;
      }
      else c = 0;
     
      postprob = get_postprob(probs, scur, sprv, k, i); 

      if ((status = p7_trace_fs_AppendWithPP(tr, scur, k, i, c, postprob)) != eslOK) return status;

      /* For NCJ, we had to defer i decrement. */
      if ( (scur == p7T_N || scur == p7T_C || scur == p7T_J) && scur == sprv) i--;
      sprv = scur;
      i-=c;
    }
  tr->M = gm_fs->M;
  tr->L = gx->L;
  return p7_trace_fs_Reverse(tr);
}

static inline float
get_postprob(const P7_GMX *pp, int scur, int sprv, int k, int i)
{
  float **dp  = pp->dp;
  float  *xmx = pp->xmx;
  switch (scur) {
  case p7T_M: return MMX_FS(i,k,p7G_C0);
  case p7T_I: return IMX_FS(i,k);
  case p7T_N: if (sprv == scur) return XMX_FS(i,p7G_N);
  case p7T_C: if (sprv == scur) return XMX_FS(i,p7G_C); 
  case p7T_J: if (sprv == scur) return XMX_FS(i,p7G_J); 
  default:    return 0.0;
  }
}

static inline int
select_m(const P7_FS_PROFILE *gm_fs, const P7_GMX *gx, int c, int i, int k)
{
  float      **dp   = gx->dp;  /* so {MDI}MX() macros work       */
  float       *xmx  = gx->xmx; /* so XMX() macro works           */
  float const *tsc  = gm_fs->tsc;  /* so TSCDELTA() macro works */
  float path[4];
  int   state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };
 
  path[0] = TSCDELTA(p7P_MM, k-1) + MMX(i,k-1);
  path[1] = TSCDELTA(p7P_IM, k-1) + IMX(i,k-1);
  path[2] = TSCDELTA(p7P_DM, k-1) + DMX(i,k-1);
  path[3] = TSCDELTA(p7P_BM, k-1) + XMX(i,p7G_B);

  return state[esl_vec_FArgMax(path, 4)];

}

static inline int
select_d(const P7_FS_PROFILE *gm_fs, const P7_GMX *gx, int i, int k)
{
  float      **dp   = gx->dp;  /* so {MDI}MX() macros work       */
  float const *tsc  = gm_fs->tsc;  /* so TSCDELTA() macro works */
  float        path[2];

  path[0] = TSCDELTA(p7P_MD, k-1) + MMX(i, k-1);
  path[1] = TSCDELTA(p7P_DD, k-1) + DMX(i, k-1);
  return ((path[0] >= path[1]) ? p7T_M : p7T_D);
}
  
static inline int
select_i(const P7_FS_PROFILE *gm_fs, const P7_GMX *gx, int i, int k)
{
  float      **dp   = gx->dp;  /* so {MDI}MX() macros work       */
  float const *tsc  = gm_fs->tsc;  /* so TSCDELTA() macro works */
  float        path[2];

  path[0] = TSCDELTA(p7P_MI, k) + MMX(i-3,k);
  path[1] = TSCDELTA(p7P_II, k) + IMX(i-3,k);

  return ((path[0] >= path[1]) ? p7T_M : p7T_I);
}

static inline int
select_n(int i)
{
  return ((i==0) ? p7T_S : p7T_N);
}

static inline int
select_c(const P7_FS_PROFILE *gm_fs, const P7_GMX *pp, const P7_GMX *gx, int i)
{
  float  t1   =  ( (gm_fs->xsc[p7P_C][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  float  t2   =  ( (gm_fs->xsc[p7P_E][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  float *xmx  = gx->xmx;  /* so XMX() macro works           */
  float  path[4];
  int    state[4] = { p7T_C, p7T_C, p7T_C, p7T_E };

  /* If we have gotten all the way up to the last possible codon enter the model */
  if(i < 4)                                                           return p7T_E;

  /* Possible paths include any of the three C state frames or moving to the E state at the current index */
  path[0] = t1 * expf(p7_FLogsum(XMX(i-3, p7G_C), pp->xmx[i*p7G_NXCELLS + p7G_C]));
  if(i < gx->L)
    path[1] = t1 * expf(p7_FLogsum(XMX(i-2, p7G_C), pp->xmx[(i+1)*p7G_NXCELLS + p7G_C]));
  else
    path[1] = FLT_MIN;
  if(i < gx->L-1)
    path[2] = t1 * expf(p7_FLogsum(XMX(i-1, p7G_C), pp->xmx[(i+2)*p7G_NXCELLS + p7G_C]));
  else
    path[2] = FLT_MIN;
  path[3] = t2 *  expf(XMX(i,p7G_E));
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
select_j(const P7_FS_PROFILE *gm_fs, const P7_GMX *pp, const P7_GMX *gx, int i)
{
  float  t1   = ( (gm_fs->xsc[p7P_J][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  float  t2   = ( (gm_fs->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  float *xmx  = gx->xmx;  /* so XMX() macro works           */
  float  path[2];
  int    state[2] = { p7T_J, p7T_E };
 
  if(i <= 5) return p7T_E;

  path[0] = t1 * expf(p7_FLogsum(XMX(i,p7G_J), pp->xmx[i*p7G_NXCELLS + p7G_J]));
  path[1] = t2 * expf(XMX(i,p7G_E));
  return state[esl_vec_FArgMax(path, 2)];
}

static inline int
select_e(const P7_FS_PROFILE *gm_fs, const P7_GMX *gx, int i, int *ret_k)
{
  float **dp   = gx->dp;  /* so {MDI}MX() macros work       */
  float   max  = -eslINFINITY;
  int     smax = -1;    /* will be returned as "error code" if no max found */
  int     kmax = -1;
  int     k;

  if (! p7_fs_profile_IsLocal(gm_fs)) /* glocal/global is easier */
    {
      *ret_k = gm_fs->M;
      return ((expf(MMX(i,gm_fs->M)) >= expf(DMX(i,gm_fs->M))) ? p7T_M : p7T_D);
    }

  for (k = 1; k <= gm_fs->M; k++)
    {
      if (expf(MMX(i,k)) >  max) { max = expf(MMX(i,k)); smax = p7T_M; kmax = k; }
      if (expf(DMX(i,k)) >  max) { max = expf(DMX(i,k)); smax = p7T_D; kmax = k; }
    }
  *ret_k = kmax;
  return smax;
}

  
static inline int
select_b(const P7_FS_PROFILE *gm_fs, const P7_GMX *gx, int i)
{
  float t1 = ( (gm_fs->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  float t2 = ( (gm_fs->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  float *xmx  = gx->xmx;  /* so XMX() macro works           */
  float path[2];
  
  path[0] = t1 * expf(XMX(i, p7G_N));
  path[1] = t2 * expf(XMX(i, p7G_J));
  return  ((path[0] > path[1]) ? p7T_N : p7T_J);

}



/*------------------------ end, oa traceback --------------------*/




/*****************************************************************
 * 3. Benchmark driver
 *****************************************************************/
#ifdef p7GENERIC_OPTACC_BENCHMARK
/*
   icc -O3 -static -o generic_optacc_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_OPTACC_BENCHMARK generic_optacc.c -lhmmer -leasel -lm
   ./benchmark-generic-optacc <hmmfile>
                   RRM_1 (M=72)       Caudal_act (M=136)      SMC_N (M=1151)
                 -----------------    ------------------     -------------------
   20 Aug 08:    67.96u (21.2 Mc/s)   128.14u (21.2 Mc/s)    1091.90u (21.1 Mc/s)
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
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,   "5000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "--notrace", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark the DP fill stage",               0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for optimal accuracy alignment, generic version";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMX         *gx1     = NULL;
  P7_GMX         *gx2     = NULL;
  P7_TRACE       *tr      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           fsc, bsc, accscore;
  double          Mcs;

  p7_FLogsumInit();
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL);
  gx1 = p7_gmx_Create(gm->M, L);
  gx2 = p7_gmx_Create(gm->M, L);
  tr  = p7_trace_CreateWithPP();

  esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  p7_GForward (dsq, L, gm, gx1, &fsc);
  p7_GBackward(dsq, L, gm, gx2, &bsc);
  p7_GDecoding(gm, gx1, gx2, gx2);                   /* <gx2> is now the posterior decoding matrix */

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      p7_GOptimalAccuracy(gm, gx2, gx1, &accscore);       /* <gx1> is now the OA matrix */

      if (! esl_opt_GetBoolean(go, "--notrace"))
  {
    p7_GOATrace(gm, gx2, gx1, tr);
    p7_trace_Reuse(tr);
  }
    }
  esl_stopwatch_Stop(w);
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / w->user;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n", gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_trace_Destroy(tr);
  p7_gmx_Destroy(gx1);
  p7_gmx_Destroy(gx2);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_OPTACC_BENCHMARK*/
/*---------------- end, benchmark driver ------------------------*/


/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7GENERIC_OPTACC_TESTDRIVE

#endif /*p7GENERIC_OPTACC_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/


/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7GENERIC_OPTACC_TESTDRIVE

#endif /*p7GENERIC_OPTACC_TESTDRIVE*/
/*------------------ end, test driver ---------------------------*/




/*****************************************************************
 * 6. Example
 *****************************************************************/
#ifdef p7GENERIC_OPTACC_EXAMPLE
/* 
   gcc -g -Wall -o generic_optacc_example -Dp7GENERIC_OPTACC_EXAMPLE -I. -I../easel -L. -L../easel generic_optacc.c -lhmmer -leasel -lm
   ./generic_optacc_example <hmmfile> <seqfile>
*/
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-d",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "dump posterior residue decoding matrix",           0 },
  { "-m",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "dump OA matrix",                                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of optimal accuracy alignment, generic implementation";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_GMX         *gx1     = NULL;
  P7_GMX         *gx2     = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  P7_TRACE       *tr      = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  char            errbuf[eslERRBUFSIZE];
  float           fsc, bsc, vsc;
  float           accscore;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);
 
  /* Read in one sequence */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
  if  (esl_sqio_Read(sqfp, sq) != eslOK) p7_Fail("Failed to read sequence");
  esl_sqfile_Close(sqfp);
 
  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL); /* multihit local: H3 default */
  
  /* Allocations */
  gx1 = p7_gmx_Create(gm->M, sq->n);
  gx2 = p7_gmx_Create(gm->M, sq->n);
  tr  = p7_trace_CreateWithPP();

  p7_FLogsumInit();
  /* Run Forward, Backward; do OA fill and trace */
  p7_GForward (sq->dsq, sq->n, gm, gx1, &fsc);
  p7_GBackward(sq->dsq, sq->n, gm, gx2, &bsc);
  p7_GDecoding(gm, gx1, gx2, gx2);                   /* <gx2> is now the posterior decoding matrix */
  p7_GOptimalAccuracy(gm, gx2, gx1, &accscore);       /* <gx1> is now the OA matrix */
  p7_GOATrace(gm, gx2, gx1, tr);

  if (esl_opt_GetBoolean(go, "-d")) p7_gmx_Dump(stdout, gx2, p7_DEFAULT);
  if (esl_opt_GetBoolean(go, "-m")) p7_gmx_Dump(stdout, gx1, p7_DEFAULT);

  p7_trace_Dump(stdout, tr, gm, sq->dsq);
  if (p7_trace_Validate(tr, abc, sq->dsq, errbuf) != eslOK) p7_Die("trace fails validation:\n%s\n", errbuf);

  printf("fwd = %.4f nats\n", fsc);
  printf("bck = %.4f nats\n", bsc);
  printf("acc = %.4f (%.2f%%)\n", accscore, accscore * 100. / (float) sq->n);

  p7_trace_Reuse(tr);

  p7_GViterbi(sq->dsq, sq->n, gm, gx1, &vsc);
  p7_GTrace  (sq->dsq, sq->n, gm, gx1, tr);
  p7_trace_SetPP(tr, gx2);
  p7_trace_Dump(stdout, tr, gm, sq->dsq);

  printf("vit = %.4f nats\n", vsc);
  printf("acc = %.4f\n", p7_trace_GetExpectedAccuracy(tr));

  /* Cleanup */
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_gmx_Destroy(gx1);
  p7_gmx_Destroy(gx2);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_OPTACC_EXAMPLE*/
/*-------------------- end, example -----------------------------*/




