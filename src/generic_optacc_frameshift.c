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
#define TSCDELTA(s,k) ( (tsc[(k) * p7P_NTRANS + (s)] == -eslINFINITY) ? FLT_MIN : 1.0)
/* The TSCDELTA is used to make impossible paths impossible in the
 * optimal accuracy decoding algorithm;
 */

/* The OA matrix uses the frameshift (p7G_NSCELLS_FS) cell layout, so the
 * fill can write it over the posterior matrix. Each filled cell holds
 * the OA scores of M, I and D in the C1, C2 and D slots, the codon length
 * the trace takes for an M there in C3, and the posteriors pp_C0 and
 * pp_I, which the trace reports, in C0 and I. C4 and C5 keep pp_C4 and
 * pp_C5, which nothing reads after the fill.
 */
#define MMX_OA(i,k)   (dp[(i)][(k) * p7G_NSCELLS_FS + p7G_M + p7G_C1])
#define IMX_OA(i,k)   (dp[(i)][(k) * p7G_NSCELLS_FS + p7G_M + p7G_C2])
#define DMX_OA(i,k)   (dp[(i)][(k) * p7G_NSCELLS_FS + p7G_D])
#define CODON_OA(i,k) (dp[(i)][(k) * p7G_NSCELLS_FS + p7G_M + p7G_C3])

#define PPM_FS(k,c)   (pp[(k) * p7G_NSCELLS_FS + p7G_M + (c)])
#define PPI_FS(k)     (pp[(k) * p7G_NSCELLS_FS + p7G_I])

/* oa_fs_max()
 * For p7_GOptimalAccuracyRow_Frameshift(): the best OA path into an M
 * at (i,k) through codon length <c>, from row i-c (cells at k-1 and B).
 */
static inline float
oa_fs_max(const float *tsc, float **dp, const float *xmx, const float *pp, int i, int k, int c)
{
  return ESL_MAX( TSCDELTA(p7P_MM, k-1) * (MMX_OA(i-c,k-1)   + PPM_FS(k, p7G_C0 + c)),
         ESL_MAX( TSCDELTA(p7P_IM, k-1) * (IMX_OA(i-c,k-1)   + PPM_FS(k, p7G_C0 + c)),
         ESL_MAX( TSCDELTA(p7P_DM, k-1) * (DMX_OA(i-c,k-1)   + PPM_FS(k, p7G_C0 + c)),
                  TSCDELTA(p7P_BM, k-1) * (XMX(i-c,p7G_B)    + PPM_FS(k, p7G_C0 + c)))));
}

/* Function:  p7_GOptimalAccuracyRow_Frameshift()
 * Synopsis:  One row of the optimal accuracy DP fill, frameshift model.
 *
 * Purpose:   OA row <i> (0..L) of the 5-codon profile <gm_fs5> into
 *            the M, I and D OA cells of <dp[i]> (the C1, C2 and D slots
 *            of the <p7G_NSCELLS_FS> layout) and row <i> of the special
 *            states <xmx>, from OA rows i-1..i-5 and the posterior row
 *            <pp> (one row, [k*p7G_NSCELLS_FS+s]) with its special
 *            states <ppx>. Row 0 is the initialization (<pp>, <ppx>
 *            unused). <dp[i]>, <xmx> row i may be the posterior row
 *            itself: each posterior is read before the OA score at its
 *            place is written, and the C0, C3, C4, C5 and I slots are
 *            not written. Used by <p7_GOptimalAccuracy_Frameshift()> for
 *            each row and by the OA sweep
 *            <p7_GOptimalAccuracySweep_Frameshift()>.
 */
void
p7_GOptimalAccuracyRow_Frameshift(const P7_FS_PROFILE *gm_fs5, const float *pp, const float *ppx, float **dp, float *xmx, int i)
{
  float const *tsc  = gm_fs5->tsc;
  int          M    = gm_fs5->M;
  float        esc  = p7_fs_profile_IsLocal(gm_fs5) ? 1.0 : 0.0;
  float        nn   = ((gm_fs5->xsc[p7P_N][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  float        jj   = ((gm_fs5->xsc[p7P_J][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  float        cc   = ((gm_fs5->xsc[p7P_C][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  float        nb   = ((gm_fs5->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  float        jb   = ((gm_fs5->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  float        ej   = ((gm_fs5->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  float        ec   = ((gm_fs5->xsc[p7P_E][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  float        max1, max2, max3, max4, max5;
  int          k;

  if (i == 0)
    {
      XMX(0,p7G_N) = 0.;                                         /* S->N, p=1            */
      XMX(0,p7G_B) = 0.;                                         /* S->N->B, no N-tail   */
      XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) = -eslINFINITY; /* need seq to get here */
      for (k = 0; k <= M; k++)
        MMX_OA(0,k) = IMX_OA(0,k) = DMX_OA(0,k) = -eslINFINITY;  /* need seq to get here */
      return;
    }

  MMX_OA(i,0) = IMX_OA(i,0) = DMX_OA(i,0) = XMX(i,p7G_E) = -eslINFINITY;

  if (i == 1)
    {
      for (k = 1; k < M; k++) {
        MMX_OA(1,k)  = TSCDELTA(p7P_BM,k-1) * PPM_FS(k,p7G_C1);
        IMX_OA(1,k)  = -eslINFINITY;
        DMX_OA(1,k)  = ESL_MAX( TSCDELTA(p7P_MD, k-1) * MMX_OA(1,k-1),
                                TSCDELTA(p7P_DD, k-1) * DMX_OA(1,k-1));
        XMX(1,p7G_E) = ESL_MAX(XMX(1,p7G_E), esc * MMX_OA(1,k));
      }
      MMX_OA(1,M) = TSCDELTA(p7P_BM, M-1) * PPM_FS(M,p7G_C1);
      IMX_OA(1,M)  = -eslINFINITY;
      DMX_OA(1,M)  = ESL_MAX( TSCDELTA(p7P_MD, M-1) * MMX_OA(1,M-1),
                              TSCDELTA(p7P_DD, M-1) * DMX_OA(1,M-1));
      XMX(1,p7G_E) = ESL_MAX(XMX(1,p7G_E), ESL_MAX(MMX_OA(1,M), DMX_OA(1, M)));

      XMX(1,p7G_J) = ej * XMX(1,p7G_E);
      XMX(1,p7G_C) = ec * XMX(1,p7G_E);
      XMX(1,p7G_N) = nn * ppx[p7G_N];
      XMX(1,p7G_B) = ESL_MAX(nb * XMX(1,p7G_N), jb * XMX(1,p7G_J));
      return;
    }

  if (i == 2)
    {
      for (k = 1; k <= M; k++) {
        max1 = oa_fs_max(tsc, dp, xmx, pp, 2, k, 1);
        max2 = TSCDELTA(p7P_BM, k-1) * (XMX(0,p7G_B) + PPM_FS(k,p7G_C2));
        MMX_OA(2,k)  = ESL_MAX( max1, max2);
        IMX_OA(2,k)  = -eslINFINITY;
        DMX_OA(2,k)  = ESL_MAX( TSCDELTA(p7P_MD, k-1) * MMX_OA(2,k-1),
                                TSCDELTA(p7P_DD, k-1) * DMX_OA(2,k-1));
        if (k < M) XMX(2,p7G_E) = ESL_MAX(XMX(2,p7G_E), esc * MMX_OA(2,k));
      }
      XMX(2,p7G_E) = ESL_MAX(XMX(2,p7G_E), ESL_MAX(MMX_OA(2,M), DMX_OA(2, M)));

      XMX(2,p7G_J) = ej * XMX(2,p7G_E);
      XMX(2,p7G_C) = ec * XMX(2,p7G_E);
      XMX(2,p7G_N) = nn * ppx[p7G_N];
      XMX(2,p7G_B) = ESL_MAX(nb * XMX(2,p7G_N), jb * XMX(2,p7G_J));
      return;
    }

  for (k = 1; k <= M; k++)
    {
      max1 = oa_fs_max(tsc, dp, xmx, pp, i, k, 1);
      max2 = oa_fs_max(tsc, dp, xmx, pp, i, k, 2);
      max3 = oa_fs_max(tsc, dp, xmx, pp, i, k, 3);
      if (i >= 5)
        {
          max4 = oa_fs_max(tsc, dp, xmx, pp, i, k, 4);
          max5 = oa_fs_max(tsc, dp, xmx, pp, i, k, 5);
          MMX_OA(i,k) = ESL_MAX(ESL_MAX(max1, max2),
                                ESL_MAX(ESL_MAX(max3, max4), max5));
        }
      else
        {
          max4 = (i == 4) ? TSCDELTA(p7P_BM, k-1) * (XMX(0,p7G_B) + PPM_FS(k,p7G_C4)) : -eslINFINITY;
          MMX_OA(i,k) = ESL_MAX( ESL_MAX( max1, max2),
                                 ESL_MAX( max3, max4));
        }

      if (k < M)
        IMX_OA(i,k) = ESL_MAX( TSCDELTA(p7P_MI, k) * (MMX_OA(i-3,k) + PPI_FS(k)),
                               TSCDELTA(p7P_II, k) * (IMX_OA(i-3,k) + PPI_FS(k)));
      else
        IMX_OA(i,k) = -eslINFINITY;

      DMX_OA(i,k) = ESL_MAX( TSCDELTA(p7P_MD, k-1) * MMX_OA(i,k-1),
                             TSCDELTA(p7P_DD, k-1) * DMX_OA(i,k-1));

      if (k < M) XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), esc * MMX_OA(i,k));
    }
  XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), ESL_MAX(MMX_OA(i,M), DMX_OA(i, M)));

  XMX(i,p7G_J) = ESL_MAX( jj * (XMX(i-3,p7G_J) + ppx[p7G_J]),
                          ej *  XMX(i,  p7G_E));
  XMX(i,p7G_C) = ESL_MAX( cc * (XMX(i-3,p7G_C) + ppx[p7G_C]),
                          ec *  XMX(i,  p7G_E));
  XMX(i,p7G_N) = nn * (XMX(i-3,p7G_N) + ppx[p7G_N]);
  XMX(i,p7G_B) = ESL_MAX( nb * XMX(i,  p7G_N), jb * XMX(i,  p7G_J));
}


/* Function:  p7_GOptimalAccuracy_Frameshift()
 * Synopsis:  Optimal accuracy decoding for frameshift aware alignments: fill. 
 *
 * Purpose:   Calculates the fill step of the optimal accuracy decoding
 *            algorithm.
 *            
 *            Caller provides the posterior decoding matrix <pp>,
 *            which was calculated by Forward/Backward on a target sequence
 *            of length <L> using the query model <gm_fs5>.
 *            
 *            Caller also provides a DP matrix <gx> with the frameshift
 *            cell layout (<p7G_NSCELLS_FS>), allocated for the
 *            <gm_fs5->M> by <pp->L> comparison, which may be <pp>
 *            itself. Each row's OA scores of M, I and D (C1, C2 and D
 *            slots) and special states are computed by
 *            <p7_GOptimalAccuracyRow_Frameshift()>; each cell also gets
 *            the codon length that <p7_GOATrace_Frameshift()> takes for
 *            an M there (C3; the argmax of pp_C1..C5, the first on
 *            ties, found before the row is computed) and the posteriors
 *            pp_C0 and pp_I (C0 and I). A cell's posteriors are read
 *            before the cell is written, and the special states of a
 *            row before they are overwritten, so the OA matrix can be
 *            written over <pp>.
 *            The trace also needs the posteriors of the special states;
 *            when <gx> is <pp>, the caller copies <pp->xmx> before the
 *            fill.
 *
 * Args:      gm_fs5 - query profile      
 *            pp    - posterior decoding matrix created by <p7_GDecoding_Frameshift()>
 *            gx    - RESULT: caller provided DP matrix for <gm_fs5->M> by <L>; may be <pp> 
 *            ret_e - RETURN: expected number of correctly decoded positions 
 *
 * Returns:   <eslOK> on success, and <*ret_e> contains the final OA
 *            score, which is the expected number of correctly decoded
 *            positions in the target sequence (up to <L>).
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_GOptimalAccuracy_Frameshift(const P7_FS_PROFILE *gm_fs5, const P7_GMX *pp, P7_GMX *gx, float *ret_e)
{
  int     L     = pp->L;
  int     M     = gm_fs5->M;
  float  *xmx   = gx->xmx;
  float  *cbest = NULL;     /* codon length (1..5) an M at (i,k) takes in the trace */
  float   codon[p7P_5CODONS];
  float  *pc, *gc;
  int     i, k, c;
  int     status;

  gx->M = M;
  gx->L = L;
  ESL_ALLOC(cbest, sizeof(float) * (M+1));

  p7_GOptimalAccuracyRow_Frameshift(gm_fs5, NULL, NULL, gx->dp, xmx, 0);
  for (i = 1; i <= L; i++)
    {
      for (k = 1; k <= M; k++)
        {
          for (c = 0; c < p7P_5CODONS; c++)
            codon[c] = pp->dp[i][k*p7G_NSCELLS_FS + p7G_M + p7G_C1 + c];
          cbest[k] = (float) (esl_vec_FArgMax(codon, p7P_5CODONS) + 1);
        }

      p7_GOptimalAccuracyRow_Frameshift(gm_fs5, pp->dp[i], pp->xmx + i*p7G_NXCELLS, gx->dp, xmx, i);

      /* when <gx> is <pp> these copies leave the values where they are */
      for (k = 1; k <= M; k++)
        {
          pc = pp->dp[i] + k*p7G_NSCELLS_FS;
          gc = gx->dp[i] + k*p7G_NSCELLS_FS;
          gc[p7G_M + p7G_C0] = pc[p7G_M + p7G_C0];
          gc[p7G_M + p7G_C4] = pc[p7G_M + p7G_C4];
          gc[p7G_M + p7G_C5] = pc[p7G_M + p7G_C5];
          gc[p7G_I]          = pc[p7G_I];
          gc[p7G_M + p7G_C3] = cbest[k];
        }
    }

  *ret_e = XMX(L,  p7G_C) +
           XMX(L-1,p7G_C) +
           XMX(L-2,p7G_C);

  free(cbest);
  return eslOK;

 ERROR:
  return status;
}

   
/*---------------------- end, oa fill ---------------------------*/

/*****************************************************************
 * 2. Optimal alignment accuracy, traceback
 *****************************************************************/

static inline float get_postprob(const P7_GMX *gx, const float *ppx, int scur, int sprv, int k, int i);
static inline int select_m(const P7_FS_PROFILE *gm_fs5, const P7_GMX *gx, int i, int k);
static inline int select_d(const P7_FS_PROFILE *gm_fs5, const P7_GMX *gx, int i, int k);
static inline int select_i(const P7_FS_PROFILE *gm_fs5, const P7_GMX *gx, int i, int k);
static inline int select_n(int i);
static inline int select_c(const P7_FS_PROFILE *gm_fs5, const float *ppx, const P7_GMX *gx, int i);
static inline int select_j(const P7_FS_PROFILE *gm_fs5, const float *ppx, const P7_GMX *gx, int i);
static inline int select_e(const P7_FS_PROFILE *gm_fs5, const P7_GMX *gx, int i, int *ret_k);
static inline int select_b(const P7_FS_PROFILE *gm_fs5, const P7_GMX *gx, int i);
static inline int select_codon(const P7_GMX *gx, int i, int k);

/* Function:  p7_GOATrace_Frameshift()
 * Synopsis:  Optimal accuracy decoding: traceback.
 *
 * Purpose:   The traceback stage of the frameshift aware optimal 
 *            accuracy decoding algorithm
 *            
 *            Caller provides the OA DP matrix <gx> that was just
 *            calculated by <p7_GOptimalAccuracy_Frameshift()> for a
 *            target sequence of length <L> and the query model
 *            <gm_fs5>, and <ppx>, the posteriors of the special states
 *            (<pp->xmx> of the posterior decoding matrix, copied before
 *            the fill when the fill wrote <gx> over <pp>). The codon
 *            length of each M, and the posteriors of M and I, are in
 *            the cells of <gx>.
 *            
 *            Caller provides an empty traceback structure <tr> to
 *            hold the result, allocated to hold optional posterior
 *            probability annotation on residues (with
 *            <p7_trace_fs_CreateWithPP()>).  This will be
 *            internally reallocated as needed for larger traces.
 *
 * Args:      gm_fs5 - query profile      
 *            ppx   - posteriors of the special states, <p7G_NXCELLS> per row 0..L
 *            gx    - OA DP matrix calculated by  <p7_GOptimalAccuracy_Frameshift()>
 *            tr    - RESULT: OA traceback, allocated with posterior probs
 *
 * Returns:   <eslOK> on success, and <tr> contains the OA traceback.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_GOATrace_Frameshift(const P7_FS_PROFILE *gm_fs5, const float *ppx, const P7_GMX *gx, P7_TRACE *tr)
{
  int          i   = gx->L; /* position in seq (1..L)         */
  int          k   = 0;     /* position in model (1..M)       */
  int          c   = 0;
  float        postprob;
  int          sprv, scur;
  int          status;

#if eslDEBUGLEVEL > 0
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif

  if ((status = p7_trace_fs_AppendWithPP(tr, p7T_T, k, i, c, 0.0f)) != eslOK) return status;
  if ((status = p7_trace_fs_AppendWithPP(tr, p7T_C, k, i, c, 0.0f)) != eslOK) return status;

  sprv = p7T_C;
  while (sprv != p7T_S)
    {

     switch (sprv) {
      case p7T_M: scur = select_m(gm_fs5,     gx, i,  k); k--;  break;
      case p7T_D: scur = select_d(gm_fs5,     gx, i,  k); k--;  break;
      case p7T_I: scur = select_i(gm_fs5,     gx, i,  k); i-=3; break;
      case p7T_N: scur = select_n(                i);           break;
      case p7T_C: scur = select_c(gm_fs5, ppx, gx, i);          break;
      case p7T_J: scur = select_j(gm_fs5, ppx, gx, i);          break;
      case p7T_E: scur = select_e(gm_fs5,     gx, i, &k);       break;
      case p7T_B: scur = select_b(gm_fs5,     gx, i);           break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }
      if (scur == -1) ESL_EXCEPTION(eslEINVAL, "OA traceback choice failed %s", gm_fs5->name);

      postprob = get_postprob(gx, ppx, scur, sprv, k, i);

      if(scur == p7T_M) c = select_codon(gx, i, k);
      else              c = 0;

      if ((status = p7_trace_fs_AppendWithPP(tr, scur, k, i, c, postprob)) != eslOK) return status;

      /* For NCJ, we had to defer i decrement. */
      if ( (scur == p7T_N || scur == p7T_C || scur == p7T_J) && scur == sprv) i--;
      sprv = scur;
      i -= c;
    }
  tr->M = gm_fs5->M;
  tr->L = gx->L;
  return p7_trace_fs_Reverse(tr);
}


/* get_postprob()
 * Posterior probability of state <scur> at (<i>,<k>): pp_C0 and pp_I as
 * the fill kept them in the cells of <gx>, and the special states'
 * posteriors from <ppx>.
 */
static inline float
get_postprob(const P7_GMX *gx, const float *ppx, int scur, int sprv, int k, int i)
{
  float **dp  = gx->dp;
  
  switch (scur) {
  case p7T_M: return MMX_FS(i,k,p7G_C0);
  case p7T_I: return IMX_FS(i,k);
  case p7T_N: if (sprv == scur) return ppx[i*p7G_NXCELLS + p7G_N];
  case p7T_C: if (sprv == scur) return ppx[i*p7G_NXCELLS + p7G_C]; 
  case p7T_J: if (sprv == scur) return ppx[i*p7G_NXCELLS + p7G_J]; 
  default:    return 0.0;
  }
}

static inline int
select_m(const P7_FS_PROFILE *gm_fs5, const P7_GMX *gx, int i, int k)
{
  float      **dp   = gx->dp;  /* so {MDI}MX() macros work       */
  float       *xmx  = gx->xmx; /* so XMX() macro works           */
  float const *tsc  = gm_fs5->tsc;  /* so TSCDELTA() macro works */
  float path[4];
  int   state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };

  path[0] = TSCDELTA(p7P_MM, k-1) * MMX_OA(i,k-1);
  path[1] = TSCDELTA(p7P_IM, k-1) * IMX_OA(i,k-1);
  path[2] = TSCDELTA(p7P_DM, k-1) * DMX_OA(i,k-1); 
  path[3] = TSCDELTA(p7P_BM, k-1) * XMX(i,p7G_B);

  return state[esl_vec_FArgMax(path, 4)];

}

static inline int
select_d(const P7_FS_PROFILE *gm_fs5, const P7_GMX *gx, int i, int k)
{
  float      **dp   = gx->dp;  /* so {MDI}MX() macros work       */
  float const *tsc  = gm_fs5->tsc;  /* so TSCDELTA() macro works */
  float        path[2];

  path[0] = TSCDELTA(p7P_MD, k-1) * MMX_OA(i, k-1);
  path[1] = TSCDELTA(p7P_DD, k-1) * DMX_OA(i, k-1);
  return ((path[0] >= path[1]) ? p7T_M : p7T_D);
}
 
  
static inline int
select_i(const P7_FS_PROFILE *gm_fs5, const P7_GMX *gx, int i, int k)
{
  float      **dp   = gx->dp;  /* so {MDI}MX() macros work       */
  float const *tsc  = gm_fs5->tsc;  /* so TSCDELTA() macro works */
  float        path[2];

  path[0] = TSCDELTA(p7P_MI, k) * MMX_OA(i-3,k);
  path[1] = TSCDELTA(p7P_II, k) * IMX_OA(i-3,k);

  return ((path[0] >= path[1]) ? p7T_M : p7T_I);
}


static inline int
select_n(int i)
{
  return ((i==0) ? p7T_S : p7T_N);
}

static inline int
select_c(const P7_FS_PROFILE *gm_fs5, const float *ppx, const P7_GMX *gx, int i)
{
  float  t1   =  ( (gm_fs5->xsc[p7P_C][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  float  t2   =  ( (gm_fs5->xsc[p7P_E][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  float *xmx  = gx->xmx;  /* so XMX() macro works           */
  float  path[4];
  int    state[4] = { p7T_C, p7T_C, p7T_C, p7T_E };

  /* If we have gotten all the way up to the last possible codon enter the model */
  if(i < 4)                                                           return p7T_E;

  /* Possible paths include any of the three C state frames or moving to the E state at the current index */
  path[0] = t1 * (XMX(i-3, p7G_C) + ppx[i*p7G_NXCELLS + p7G_C]);

  if(i < gx->L) path[1] = t1 * (XMX(i-2, p7G_C) + ppx[(i+1)*p7G_NXCELLS + p7G_C]);
  else          path[1] = FLT_MIN;

  if(i < gx->L-1) path[2] = t1 * (XMX(i-1, p7G_C) + ppx[(i+2)*p7G_NXCELLS + p7G_C]);
  else            path[2] = FLT_MIN;

  path[3] = t2 * XMX(i,p7G_E);

  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
select_j(const P7_FS_PROFILE *gm_fs5, const float *ppx, const P7_GMX *gx, int i)
{
  float  t1   = ( (gm_fs5->xsc[p7P_J][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  float  t2   = ( (gm_fs5->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  float *xmx  = gx->xmx;  /* so XMX() macro works           */
  float  path[2];
  int    state[2] = { p7T_J, p7T_E };

  if(i <= 5) return p7T_E;

  path[0] = t1 * (XMX(i,p7G_J) + ppx[i*p7G_NXCELLS + p7G_J]);
  path[1] = t2 *  XMX(i,p7G_E);
  return state[esl_vec_FArgMax(path, 2)];
}

static inline int
select_e(const P7_FS_PROFILE *gm_fs5, const P7_GMX *gx, int i, int *ret_k)
{
  float **dp   = gx->dp;  /* so {MDI}MX() macros work       */
  float   max  = -eslINFINITY;
  int     smax = -1;    /* will be returned as "error code" if no max found */
  int     kmax = -1;
  int     k;

  if (! p7_fs_profile_IsLocal(gm_fs5)) /* glocal/global is easier */
    {
      *ret_k = gm_fs5->M;
      return (MMX_OA(i,gm_fs5->M) >= DMX_OA(i,gm_fs5->M)) ? p7T_M : p7T_D;
    }

  for (k = 1; k <= gm_fs5->M; k++)
    {
      if (MMX_OA(i,k) >  max) { max = MMX_OA(i,k); smax = p7T_M; kmax = k; }
      if (DMX_OA(i,k) >  max) { max = DMX_OA(i,k); smax = p7T_D; kmax = k; }
    }
  *ret_k = kmax;

  return smax;
}


static inline int
select_b(const P7_FS_PROFILE *gm_fs5, const P7_GMX *gx, int i)
{
  float t1 = ( (gm_fs5->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  float t2 = ( (gm_fs5->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  float *xmx  = gx->xmx;  /* so XMX() macro works           */
  float path[2];

  path[0] = t1 * XMX(i, p7G_N);
  path[1] = t2 * XMX(i, p7G_J);
  return  ((path[0] > path[1]) ? p7T_N : p7T_J);

}

/* select_codon()
 * Codon length of the M at (<i>,<k>), which the fill stored in the cell:
 * the argmax of the posteriors pp_C1..C5.
 */
static inline int
select_codon(const P7_GMX *gx, int i, int k)
{
  float **dp   = gx->dp;

  return (int) CODON_OA(i,k);
}

/*------------------------ end, oa traceback --------------------*/




/*****************************************************************
 * 3. Benchmark driver
 *****************************************************************/
#ifdef p7GENERIC_OPTACC_FRAMESHIFT_BENCHMARK
/*
   gcc -g -O3      -o generic_optacc_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_OPTACC_FRAMESHIFT_BENCHMARK generic_optacc_frameshift.c -lhmmer -leasel -lm

   icc -O3 -static -o generic_optacc_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_OPTACC_FRAMESHIFT_BENCHMARK generic_optacc_frameshift.c -lhmmer -leasel -lm
   ./generic_optacc_frameshift_benchmark <hmmfile>
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
  P7_FS_PROFILE  *gm_fs5   = NULL;
  P7_GMX         *gx1     = NULL;
  P7_GMX         *gx2     = NULL;
  P7_GMX         *gbck    = NULL;
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
  gm_fs5 = p7_profile_fs_Create(hmm->M, abcAA, p7P_5CODONS);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L, p7_UNILOCAL);
  gx1 = p7_gmx_Create(gm_fs5->M, L, L, p7G_NSCELLS_FS);
  gx2 = p7_gmx_Create(gm_fs5->M, L, L, p7G_NSCELLS_FS);
  gbck = p7_gmx_Create(gm_fs5->M, L, L, p7G_NSCELLS);
  iv  = p7_ivx_Create(gm_fs5->M, p7P_5CODONS);
  tr  = p7_trace_fs_CreateWithPP();

  esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
  p7_GForward_Frameshift_5Codons(dsq, L, gm_fs5, gx1, iv, &fsc);
  p7_GBackward_Frameshift_5Codons(dsq, L, gm_fs5, gbck, iv, &bsc);
  p7_GDecoding_Frameshift(gm_fs5, gx1, gbck);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
  {
    p7_GOptimalAccuracy_Frameshift(gm_fs5, gx1, gx2, &accscore);
    if (! esl_opt_GetBoolean(go, "--notrace"))
    {
      p7_GOATrace_Frameshift(gm_fs5, gx1->xmx, gx2, tr);
      p7_trace_Reuse(tr);
    }  

  }
   
  esl_stopwatch_Stop(w);
  Mcs        = (double) N * (double) L * (double) gm_fs5->M * 1e-6 / w->user;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n", gm_fs5->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_trace_fs_Destroy(tr);
  p7_gmx_Destroy(gx1);
  p7_gmx_Destroy(gx2);
  p7_gmx_Destroy(gbck);
  p7_ivx_Destroy(iv);
  p7_profile_fs_Destroy(gm_fs5);
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
#endif /*p7GENERIC_OPTACC_FRAMESHIFT_BENCHMARK*/
/*---------------- end, benchmark driver ------------------------*/


