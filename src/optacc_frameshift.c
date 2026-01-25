/* Optimal accuracy alignment; frameshift version.
 * 
 * Contents:
 *   1. Optimal alignment accuracy fill.
 *   2. Optimal alignment accuracy traceback.
 *   3. Benchmark driver
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
 *            For the frameshift version, probabilities are smaller due to
 *            being normalized across multiple frames, so they are kept in 
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
  int          i,k,c;
  int          M    = gm_fs->M;
  float        esc  = p7_fs_profile_IsLocal(gm_fs) ? 1.0 : 0.0;
  float        nn, jj, cc;
  float        nb, jb, ej, ec;
  float        max1, max2, max3, max4, max5;
  float        codon[5];

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

    codon[0] = pp->dp[2][k*p7G_NSCELLS_FS + p7G_M + p7G_C1];
    codon[1] = pp->dp[2][k*p7G_NSCELLS_FS + p7G_M + p7G_C2];

    c = esl_vec_FArgMax(codon, 2) + 1;

    MMX(2,k)  = ESL_MAX( TSCDELTA(p7P_MM, k-1) + MMX(2-c,k-1),
                ESL_MAX( TSCDELTA(p7P_IM, k-1) + IMX(2-c,k-1),
                ESL_MAX( TSCDELTA(p7P_DM, k-1) + DMX(2-c,k-1),
                         TSCDELTA(p7P_BM, k-1) + XMX(2-c,p7G_B))));

    /* c+2 turns the codon length into the correct matrix index */
    MMX(2,k)  = (MMX(2,k)  != -eslINFINITY ?  p7_FLogsum(MMX(2,k), pp->dp[2][k*p7G_NSCELLS_FS + p7G_M + (c+2)]) : -eslINFINITY);
    IMX(2,k)  = -eslINFINITY;  
    DMX(2,k)  = ESL_MAX( TSCDELTA(p7P_MD, k-1) + MMX(2,k-1),
                         TSCDELTA(p7P_DD, k-1) + DMX(2,k-1));
           
    XMX(2,p7G_E) = ESL_MAX(XMX(2,p7G_E), esc * MMX(2,k));
  }

  codon[0] = pp->dp[2][M*p7G_NSCELLS_FS + p7G_M + p7G_C1];
  codon[1] = pp->dp[2][M*p7G_NSCELLS_FS + p7G_M + p7G_C2];

  c = esl_vec_FArgMax(codon, 2) + 1;

  MMX(2,M)  = ESL_MAX( TSCDELTA(p7P_MM, M-1) + MMX(2-c,M-1),
              ESL_MAX( TSCDELTA(p7P_IM, M-1) + IMX(2-c,M-1),
              ESL_MAX( TSCDELTA(p7P_DM, M-1) + DMX(2-c,M-1),
                       TSCDELTA(p7P_BM, M-1) + XMX(2-c,p7G_B))));

  /* c+2 turns the codon length into the correct matrix index */
  MMX(2,M)  = (MMX(2,M)  != -eslINFINITY ?  p7_FLogsum(MMX(2,M), pp->dp[2][M*p7G_NSCELLS_FS + p7G_M + (c+2)]) : -eslINFINITY);
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

      /* First determine the codon length */
      codon[0] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1];
      codon[1] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2];
      codon[2] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3];
      codon[3] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4];
      
      /* Use i to determine which codon lengths are permisable */
      c = esl_vec_FArgMax(codon, i) + 1;

      MMX(i,k)  = ESL_MAX( TSCDELTA(p7P_MM, k-1) + MMX(i-c,k-1),
                  ESL_MAX( TSCDELTA(p7P_IM, k-1) + IMX(i-c,k-1),
                  ESL_MAX( TSCDELTA(p7P_DM, k-1) + DMX(i-c,k-1),
                           TSCDELTA(p7P_BM, k-1) + XMX(i-c,p7G_B))));

      /* c+2 turns the codon length into the correct matrix index */
      MMX(i,k)  = (MMX(i,k)  != -eslINFINITY ?  p7_FLogsum(MMX(i,k), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + (c+2)]) : -eslINFINITY);
 
      IMX(i,k)  = ESL_MAX( TSCDELTA(p7P_MI, k) + MMX(i-3,k),
                           TSCDELTA(p7P_II, k) + IMX(i-3,k));
      IMX(i,k)  = (IMX(i,k) != -eslINFINITY ? p7_FLogsum(IMX(i,k), pp->dp[i][k*p7G_NSCELLS_FS + p7G_I]) : -eslINFINITY);
      DMX(i,k)  = ESL_MAX( TSCDELTA(p7P_MD, k-1) + MMX(i,k-1),
                           TSCDELTA(p7P_DD, k-1) + DMX(i,k-1)); 

      XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), esc * MMX(i,k));    
    }

    /* First determine the codon length */
    codon[0] = pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C1];
    codon[1] = pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C2];
    codon[2] = pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C3];
    codon[3] = pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C4];
    codon[4] = pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C5];

    c = esl_vec_FArgMax(codon, 5) + 1;

    MMX(i,M)  = ESL_MAX( TSCDELTA(p7P_MM, M-1) + MMX(i-c,M-1),
                ESL_MAX( TSCDELTA(p7P_IM, M-1) + IMX(i-c,M-1),
                ESL_MAX( TSCDELTA(p7P_DM, M-1) + DMX(i-c,M-1),
                         TSCDELTA(p7P_BM, M-1) + XMX(i-c,p7G_B))));

    /* c+2 turns the codon length into the correct matrix index */
    MMX(i,M)  = (MMX(i,M)  != -eslINFINITY ?  p7_FLogsum(MMX(i,M), pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + (c+2)]) : -eslINFINITY);

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

      /* First determine the codon length */
      codon[0] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1];
      codon[1] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C2];
      codon[2] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C3];
      codon[3] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C4];
      codon[4] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C5];
      
      c = esl_vec_FArgMax(codon, 5) + 1; 

      MMX(i,k)  = ESL_MAX( TSCDELTA(p7P_MM, k-1) + MMX(i-c,k-1),
                  ESL_MAX( TSCDELTA(p7P_IM, k-1) + IMX(i-c,k-1),
                  ESL_MAX( TSCDELTA(p7P_DM, k-1) + DMX(i-c,k-1),
                           TSCDELTA(p7P_BM, k-1) + XMX(i-c,p7G_B))));

      /* c+2 turns the codon length into the correct matrix index */
      MMX(i,k)  = (MMX(i,k)  != -eslINFINITY ?  p7_FLogsum(MMX(i,k), pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + (c+2)]) : -eslINFINITY);

      IMX(i,k)  = ESL_MAX( TSCDELTA(p7P_MI, k) + MMX(i-3,k), 
                           TSCDELTA(p7P_II, k) + IMX(i-3,k));

	  IMX(i,k)  = (IMX(i,k) != -eslINFINITY ? p7_FLogsum(IMX(i,k), pp->dp[i][k*p7G_NSCELLS_FS + p7G_I]) : -eslINFINITY);

      DMX(i,k)  = ESL_MAX( TSCDELTA(p7P_MD, k-1) + MMX(i,k-1), TSCDELTA(p7P_DD, k-1) + DMX(i,k-1));

      XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), esc * MMX(i,k));  
    } 

    /* First determine the codon length */
    codon[0] = pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C1];
    codon[1] = pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C2];
    codon[2] = pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C3];
    codon[3] = pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C4];
    codon[4] = pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + p7G_C5];

    c = esl_vec_FArgMax(codon, 5) + 1;

    MMX(i,M)  = ESL_MAX( TSCDELTA(p7P_MM, M-1) + MMX(i-c,M-1),
                ESL_MAX( TSCDELTA(p7P_IM, M-1) + IMX(i-c,M-1),
                ESL_MAX( TSCDELTA(p7P_DM, M-1) + DMX(i-c,M-1),
                         TSCDELTA(p7P_BM, M-1) + XMX(i-c,p7G_B))));

    /* c+2 turns the codon length into the correct matrix index */
    MMX(i,M)  = (MMX(i,M)  != -eslINFINITY ?  p7_FLogsum(MMX(i,M), pp->dp[i][M*p7G_NSCELLS_FS + p7G_M + (c+2)]) : -eslINFINITY);

    IMX(i,M)     = -eslINFINITY;
    DMX(i,M)     = ESL_MAX( TSCDELTA(p7P_MD, M-1) + MMX(i,M-1), TSCDELTA(p7P_DD, M-1) + DMX(i,M-1));

    XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), ESL_MAX(MMX(i,M), DMX(i, M)));

    XMX(i,p7G_J) = ESL_MAX( jj + p7_FLogsum(XMX(i-3,p7G_J), pp->xmx[i*p7G_NXCELLS + p7G_J]),
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
static inline int select_m(const P7_FS_PROFILE *gm_fs,                   const P7_GMX *gx, int i, int k);
static inline int select_d(const P7_FS_PROFILE *gm_fs,                   const P7_GMX *gx, int i, int k);
static inline int select_i(const P7_FS_PROFILE *gm_fs,                   const P7_GMX *gx, int i, int k);
static inline int select_n(int i);
static inline int select_c(const P7_FS_PROFILE *gm_fs, const P7_GMX *pp, const P7_GMX *gx, int i);
static inline int select_j(const P7_FS_PROFILE *gm_fs, const P7_GMX *pp, const P7_GMX *gx, int i);
static inline int select_e(const P7_FS_PROFILE *gm_fs,                   const P7_GMX *gx, int i, int *ret_k);
static inline int select_b(const P7_FS_PROFILE *gm_fs,                   const P7_GMX *gx, int i);


/* Function:  p7_OATrace_Frameshift()
 * Synopsis:  Optimal accuracy decoding: traceback.
 *
 * Purpose:   The traceback stage of the frameshift aware optimal 
 *            accuracy decoding algorithm
 *            
 *            Caller provides the OA DP matrix <gx> that was just
 *            calculated by <p7_OptimalAccuracy_Frameshift()>, as 
 *            well as the posterior decoding matricies <pp>, and 
 *            <probs> which were calculated by Forward/Backward on 
 *            a target sequence of length <L> using the query model 
 *           <gm_fs>.
 *            
 *            Caller provides an empty traceback structure <tr> to
 *            hold the result, allocated to hold optional posterior
 *            probability annotation on residues (with
 *            <p7_trace_fs_CreateWithPP()>).  This will be
 *            internally reallocated as needed for larger traces.
 *
 * Args:      gm_fs - query profile      
 *            pp    - posterior decoding (i normalized accross all codons containing i)
 *            gx    - OA DP matrix calculated by  <p7_OptimalAccuracyDP()>
 *            probs -  posterior decoding (i normalized accross all codons ending at i)
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
      case p7T_M: scur = select_m(gm_fs,     gx, i,  k);          k--;  break;
      case p7T_D: scur = select_d(gm_fs,     gx, i,  k);          k--;  break;
      case p7T_I: scur = select_i(gm_fs,     gx, i,  k); i -= 3;        break;
      case p7T_N: scur = select_n(               i);                    break;
      case p7T_C: scur = select_c(gm_fs, pp, gx, i);                    break;
      case p7T_J: scur = select_j(gm_fs, pp, gx, i);                    break;
      case p7T_E: scur = select_e(gm_fs,     gx, i, &k);                break;
      case p7T_B: scur = select_b(gm_fs,     gx, i);                    break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }
      if (scur == -1) ESL_EXCEPTION(eslEINVAL, "OA traceback choice failed");

      /* To prevent the trace from containing excessive frameshifts we defer to the 
       * posterior probabilities matrix in determining match state emissions lenghts */
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
select_m(const P7_FS_PROFILE *gm_fs, const P7_GMX *gx, int i, int k)
{
  float      **dp   = gx->dp;  /* so {MDI}MX() macros work       */
  float       *xmx  = gx->xmx; /* so XMX() macro works           */
  float const *tsc  = gm_fs->tsc;  /* so TSCDELTA() macro works */
  float path[4];
  int   state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };
  
  /* Codon length has already been determined and i moved to the pervious state */
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
#ifdef p7OPTACC_FRAMESHIFT_BENCHMARK
/*
    gcc -g -O2      -o optacc_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7OPTACC_FRAMESHIFT_BENCHMARK optacc_frameshift.c -lhmmer -leasel -lm

   icc -O3 -static -o optacc_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7OPTACC_FRAMESHIFT_BENCHMARK optacc_frameshift.c -lhmmer -leasel -lm
   ./optacc_frameshift_benchmark <hmmfile>
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
static char banner[] = "benchmark driver for optimal accuracy alignment, frameshift version";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcDNA  = NULL;
  ESL_ALPHABET   *abcAA   = NULL;
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
  P7_TRACE       *tr      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           fsc, bsc, accscore;
  double          Mcs;

  p7_FLogsumInit();
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  abcDNA = esl_alphabet_Create(eslDNA);
  gcode = esl_gencode_Create(abcDNA, abcAA);
  bgAA  = p7_bg_Create(abcAA);
  bgDNA = p7_bg_Create(abcDNA);
  p7_bg_SetLength(bgAA, L/3);
  p7_bg_SetLength(bgDNA, L);
  gm_fs = p7_profile_fs_Create(hmm->M, abcAA);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs, L, p7_UNILOCAL);
  gx1 = p7_gmx_fs_Create(gm_fs->M, L, L, p7P_5CODONS);
  gx2 = p7_gmx_fs_Create(gm_fs->M, L, L, 0);
  pp  = p7_gmx_fs_Create(gm_fs->M, L, L, p7P_5CODONS);
  iv  = p7_ivx_Create(gm_fs->M, p7P_5CODONS);
  tr  = p7_trace_fs_CreateWithPP();

  esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
  p7_Forward_Frameshift(dsq, gcode, L, gm_fs, gx1, iv, &fsc);
  p7_Backward_Frameshift(dsq, gcode, L, gm_fs, gx2, iv, &bsc);
  p7_Decoding_Frameshift(gm_fs, gx1, gx2, pp);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
  {
    p7_OptimalAccuracy_Frameshift(gm_fs, pp, gx2, &accscore);
    if (! esl_opt_GetBoolean(go, "--notrace"))
    {
      p7_OATrace_Frameshift(gm_fs, pp, gx2, gx1, tr);
      p7_trace_Reuse(tr);
    }  

  }
   
  esl_stopwatch_Stop(w);
  Mcs        = (double) N * (double) L * (double) gm_fs->M * 1e-6 / w->user;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n", gm_fs->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_trace_fs_Destroy(tr);
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
#endif /*p7OPTACC_FRAMESHIFT_BENCHMARK*/
/*---------------- end, benchmark driver ------------------------*/


