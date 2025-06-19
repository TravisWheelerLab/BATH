/* Global Viterbi traceback; stadard and frameshift versions. */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"

#include "hmmer.h"

/* Function: p7_global_Trace()
 * 
 * Purpose:  Traceback of a global Viterbi matrix: retrieval 
 *           of optimum alignment.
 *           
 * Args:     dsq    - digital sequence aligned to, 1..L 
 *           L      - length of <dsq>
 *           gm     - profile
 *           mx     - Viterbi matrix to trace, L x M
 *           tr     - storage for the recovered traceback.
 *           
 * Return:   <eslOK> on success.
 *           <eslFAIL> if even the optimal path has zero probability;
 *           in this case, the trace is set blank (<tr->N = 0>).
 *
 * Note:     Care is taken to evaluate the prev+tsc+emission
 *           calculations in exactly the same order that Viterbi did
 *           them, lest you get numerical problems with
 *           a+b+c = d; d-c != a+b because d,c are nearly equal.
 *           (This bug appeared in dev: xref J1/121.)
 */
int
p7_global_Trace(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, const P7_GMX *gx, P7_TRACE *tr)
{
  int          i   = L;		/* position in seq (1..L)         */
  int          k   = 0;		/* position in model (1..M)       */
  int          M   = gm->M;
  float      **dp  = gx->dp;	/* so {MDI}MX() macros work       */
  float       *xmx = gx->xmx;	/* so XMX() macro works           */
  float        tol = 1e-5;	/* floating point "equality" test */
  float const *tsc = gm->tsc;
  int     sprv, scur;		/* previous, current state in trace */
  int     status;

#if eslDEBUGLEVEL > 0
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif

  if ((status = p7_trace_Append(tr, p7T_T, k, i)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_C, k, i)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_E, k, i)) != eslOK) return status;
  sprv = p7T_E;
  while (sprv != p7T_S) {
    float const *rsc = (i>0 ? gm->rsc[dsq[i]] : NULL);
    switch (sprv) {
    case p7T_E:		/* E connects from any M, I or D state. k set here */
      if (XMX(i, p7G_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);

      if      (esl_FCompare_old(XMX(i, p7G_E), MMX(i,M), tol) == eslOK)  scur = p7T_M; 
      else if (esl_FCompare_old(XMX(i, p7G_E), DMX(i,M), tol) == eslOK)  scur = p7T_D; 
      else if (esl_FCompare_old(XMX(i, p7G_E), IMX(i,M), tol) == eslOK)  scur = p7T_I; 

      k = M;
      break;      

    case p7T_M:			/* M connects from i-1,k-1, or B */
      if (MMX(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);

      if      (esl_FCompare_old(MMX(i,k), XMX(i-1,p7G_B)                    + MSC(k), tol) == eslOK) scur = p7T_B;
      else if (esl_FCompare_old(MMX(i,k), MMX(i-1,k-1)   + TSC(p7P_MM, k-1) + MSC(k), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(MMX(i,k), IMX(i-1,k-1)   + TSC(p7P_IM, k-1) + MSC(k), tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare_old(MMX(i,k), DMX(i-1,k-1)   + TSC(p7P_DM, k-1) + MSC(k), tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
      k--; i--;
      break;

    case p7T_D:			/* D connects from M,D at i,k-1 , or B*/
      if (DMX(i, k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);
      
      if      (esl_FCompare_old(DMX(i,k), XMX(i, p7G_B),                     tol) == eslOK) scur = p7T_B;
      else if (esl_FCompare_old(DMX(i,k), MMX(i, k-1)    + TSC(p7P_MD, k-1), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(DMX(i,k), DMX(i, k-1)    + TSC(p7P_DD, k-1), tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      k--;
      break;

    case p7T_I:			/* I connects from M,I at i-1,k, or B*/
      if (IMX(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);
      if      (esl_FCompare_old(IMX(i,k), XMX(i-1,p7G_B),                  tol) == eslOK) scur = p7T_B;
      else if (esl_FCompare_old(IMX(i,k), MMX(i-1,k)     + TSC(p7P_MI, k), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(IMX(i,k), IMX(i-1,k)     + TSC(p7P_II, k), tol) == eslOK) scur = p7T_I;
      else ESL_EXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      i--;
      break;

    case p7T_N:			/* N connects from S, N */
      if (XMX(i, p7G_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = ( (i == 0) ? p7T_S : p7T_N);
      break;

    case p7T_B:			/* B connects from N, J */
      if (XMX(i,p7G_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);
      scur = p7T_N;
      break;

    default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over statetype[tpos-1] */

    /* Append this state and the current i,k to be explained to the growing trace */
    if ((status = p7_trace_Append(tr, scur, k, i)) != eslOK) return status;

    /* For NCJ, we had to defer i decrement. */
    if ( (scur == p7T_N || scur == p7T_J || scur == p7T_C) && scur == sprv) i--;

    sprv = scur;
  } /* end traceback, at S state */

  tr->M = gm->M;
  tr->L = L;
  return p7_trace_Reverse(tr);
}

/* Function: p7_fs_global_Trace()
 * 
 * Purpose:  Traceback of a global Viterbi matrix: retrieval 
 *           of optimum alignment.
 *           
 * Args:     dsq    - digital sequence aligned to, 1..L 
 *           L      - length of <dsq>
 *           gm     - profile
 *           mx     - Viterbi matrix to trace, L x M
 *           tr     - storage for the recovered traceback.
 *           
 * Return:   <eslOK> on success.
 *           <eslFAIL> if even the optimal path has zero probability;
 *           in this case, the trace is set blank (<tr->N = 0>).
 *
 * Note:     Care is taken to evaluate the prev+tsc+emission
 *           calculations in exactly the same order that Viterbi did
 *           them, lest you get numerical problems with
 *           a+b+c = d; d-c != a+b because d,c are nearly equal.
 *           (This bug appeared in dev: xref J1/121.)
 */
int
p7_fs_global_Trace(const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs, const P7_GMX *gx, P7_TRACE *tr)
{
  int          i   = L;		/* position in seq (1..L)         */
  int          k   = 0;		/* position in model (1..M)       */
  int          c   = 0;
  int          t;
  int          M   = gm_fs->M;
  float      **dp  = gx->dp;	/* so {MDI}MX() macros work       */
  float       *xmx = gx->xmx;	/* so XMX() macro works           */
  float        tol = 1e-5;	/* floating point "equality" test */
  float        match_codon[5];
  float        match_trans[3];
  float const *tsc = gm_fs->tsc;
  int     sprv, scur;		/* previous, current state in trace */
  int     status;

#if eslDEBUGLEVEL > 0
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i, c)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i, c)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_E, k, i, c)) != eslOK) return status;
  sprv = p7T_E;
  while (sprv != p7T_S) {
    switch (sprv) {
    case p7T_E:		/* E connects from any M, I or D state. k set here */
      if (XMX(i, p7G_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
      
      if      (esl_FCompare_old(XMX(i, p7G_E), MMX_FS(i,M,p7G_C0), tol) == eslOK)  scur = p7T_M; 
      else if (esl_FCompare_old(XMX(i, p7G_E), DMX_FS(i,M),        tol) == eslOK)  scur = p7T_D; 
      else if (esl_FCompare_old(XMX(i, p7G_E), IMX_FS(i,M),        tol) == eslOK)  scur = p7T_I; 
     
      if(scur == p7T_M) {
        match_codon[0] = MMX_FS(i,M,p7G_C1);
        match_codon[1] = MMX_FS(i,M,p7G_C2);
        match_codon[2] = MMX_FS(i,M,p7G_C3);
        match_codon[3] = MMX_FS(i,M,p7G_C4);
        match_codon[4] = MMX_FS(i,M,p7G_C5);
  
       c = esl_vec_FArgMax(match_codon, 5) + 1;
      }
      else if (scur == p7T_I) c = 3;
      else                   c = 0;
      k = M; 
      break;      

    case p7T_M:			/* M connects from i-1,k-1, or B */
      if (MMX_FS(i,k,p7G_C0) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);

      if(i == c && k == 1) scur = p7T_B;
      else {
        match_trans[0] = MMX_FS(i-c,k-1,p7G_C0) + TSC(p7P_MM, k-1); 
        match_trans[1] = IMX_FS(i-c,k-1)        + TSC(p7P_IM, k-1);
        match_trans[2] = DMX_FS(i-c,k-1)        + TSC(p7P_IM, k-1);

        t = esl_vec_FArgMax(match_trans, 3);
        if      ( t==0 ) scur = p7T_M;
        else if ( t==1 ) scur = p7T_I;
        else if ( t==2 ) scur = p7T_D;
      }
      k--; i-=c;

      if(scur == p7T_M) {
        match_codon[0] = MMX_FS(i,k,p7G_C1);
        match_codon[1] = MMX_FS(i,k,p7G_C2);
        match_codon[2] = MMX_FS(i,k,p7G_C3);
        match_codon[3] = MMX_FS(i,k,p7G_C4);
        match_codon[4] = MMX_FS(i,k,p7G_C5);

       c = esl_vec_FArgMax(match_codon, 5) + 1;
      }
      else if (scur == p7T_I) c = 3;
      else                    c = 0;

      break;

    case p7T_D:			/* D connects from M,D at i,k-1 , or B*/
      if (DMX_FS(i, k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);
      
      if(i == 0 && k == 1) scur = p7T_B;
      else if (esl_FCompare_old(DMX_FS(i,k), MMX_FS(i, k-1, p7G_C0) + TSC(p7P_MD, k-1), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(DMX_FS(i,k), DMX_FS(i, k-1)         + TSC(p7P_DD, k-1), tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      k--; 
      
      if(scur == p7T_M) {
        match_codon[0] = MMX_FS(i,k,p7G_C1);
        match_codon[1] = MMX_FS(i,k,p7G_C2);
        match_codon[2] = MMX_FS(i,k,p7G_C3);
        match_codon[3] = MMX_FS(i,k,p7G_C4);
        match_codon[4] = MMX_FS(i,k,p7G_C5);

       c = esl_vec_FArgMax(match_codon, 5) + 1;
      }
      else                    c = 0;

      break;

    case p7T_I:			/* I connects from M,I at i-1,k, or B*/
      if (IMX_FS(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);

      if (i == 3 && k == 0) scur = p7T_B;
      else if (esl_FCompare_old(IMX_FS(i,k), MMX_FS(i-3,k,p7G_C0)  + TSC(p7P_MI, k), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(IMX_FS(i,k), IMX_FS(i-3,k)         + TSC(p7P_II, k), tol) == eslOK) scur = p7T_I;
      else ESL_EXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      i-=3; 

      if(scur == p7T_M) {
        match_codon[0] = MMX_FS(i,k,p7G_C1);
        match_codon[1] = MMX_FS(i,k,p7G_C2);
        match_codon[2] = MMX_FS(i,k,p7G_C3);
        match_codon[3] = MMX_FS(i,k,p7G_C4);
        match_codon[4] = MMX_FS(i,k,p7G_C5);

       c = esl_vec_FArgMax(match_codon, 5) + 1;
      }
      else if (scur == p7T_I) c = 3;
      else                    c = 0;

      break;

    case p7T_N:			/* N connects from S, N */
      if (XMX(i, p7G_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = ( (i == 0) ? p7T_S : p7T_N);
      c = 0;
      break;

    case p7T_B:			/* B connects from N, J */
      if (XMX(i,p7G_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);
      scur = p7T_N;
      c = 0;
      break;

    default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over statetype[tpos-1] */

    /* Append this state and the current i,k to be explained to the growing trace */
    if ((status = p7_trace_fs_Append(tr, scur, k, i, c)) != eslOK) return status;
  
    /* For NCJ, we had to defer i decrement. */
    if ( (scur == p7T_N || scur == p7T_J || scur == p7T_C) && scur == sprv) i--;

    sprv = scur;
  } /* end traceback, at S state */

  tr->M = gm_fs->M;
  tr->L = L;
  return p7_trace_fs_Reverse(tr);
}


