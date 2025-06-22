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

int
p7_fs_semiglobal_VTrace(const ESL_DSQ *sub_dsq, int L, const ESL_GENCODE *gcode, const P7_FS_PROFILE *sub_gm, const P7_GMX *gx, P7_TRACE *tr)
{
  int          i   = L;     /* position in seq (1..L)         */
  int          k   = 0;     /* position in model (1..M)       */
  int          M   = sub_gm->M;
  int          t,u,v,w,x;
  int          c, c1, c2, c3, c4, c5;
  int          sprv,scur;
  float        emit1, emit2, emit3, emit4, emit5;
  float      **dp  = gx->dp;    /* so {MDI}MX() macros work       */
  float       *xmx = gx->xmx;   /* so XMX() macro works           */
  float        tol = 1e-5;  /* floating point "equality" test */
  float const *tsc = sub_gm->tsc; 
  int     status;

#if eslDEBUGLEVEL > 0
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i, 0)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i, 0)) != eslOK) return status;

  sprv = p7T_C;
 
  while (sprv != p7T_S) {

    switch (sprv) {
    case p7T_C:     /* C(i) comes from C(i-1) or E(i) */
      if   (XMX(i,p7G_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
  
      if      (XMX(i, p7G_C) + XMX(i-2, p7G_C) || XMX(i, p7G_C) < XMX(i-1, p7G_C))                         scur = p7T_C; 
      if      (esl_FCompare_old(XMX(i, p7G_C), XMX(i-3, p7G_C) + sub_gm->xsc[p7P_C][p7P_LOOP], tol) == eslOK)  scur = p7T_C;
      else if (esl_FCompare_old(XMX(i, p7G_C), XMX(i,   p7G_E) + sub_gm->xsc[p7P_E][p7P_MOVE], tol) == eslOK)  scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
 
      break;

    case p7T_E:     /* E connects from any M state. k set here */
      if (XMX(i, p7G_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);

      if      (esl_FCompare_old(XMX(i, p7G_E), MMX_SP(i,M), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(XMX(i, p7G_E), DMX_SP(i,M), tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
 
      k = M;
      break;

    case p7T_M:         /* M connects from i-1,k-1, or B */
      if (MMX_SP(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);
 
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i]))   x = sub_dsq[i];
      else                                                  x = p7P_MAXCODONS;      
      if(i > 1) {
        if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i-1])) w = sub_dsq[i-1];
        else                                                  w = p7P_MAXCODONS;
      } 
      if(i > 2) {
        if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i-2])) v = sub_dsq[i-2];
        else                                                  v = p7P_MAXCODONS;
      }
      if(i > 3) {
        if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i-3])) u = sub_dsq[i-3];
        else                                                  u = p7P_MAXCODONS;
      }
      if(i > 4) {
        if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i-4])) t = sub_dsq[i-4];
        else                                                  t = p7P_MAXCODONS;
      }
     
      /* find correct index for looking up scores of codons and quasicodons */
      c1 = p7P_CODON1(x);
      c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

      c2 = p7P_CODON2(w, x);
      c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

      c3 = p7P_CODON3(v, w, x);
      c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

      c4 = p7P_CODON4(u, v, w, x);
      c4 = p7P_MINIDX(c4, p7P_DEGEN_QC1);

      c5 = p7P_CODON5(t, u, v, w, x);
      c5 = p7P_MINIDX(c5, p7P_DEGEN_QC2);     

      emit1 = p7P_MSC_CODON(sub_gm, k, c1);
      emit2 = p7P_MSC_CODON(sub_gm, k, c2);
      emit3 = p7P_MSC_CODON(sub_gm, k, c3);
      emit4 = p7P_MSC_CODON(sub_gm, k, c4);
      emit5 = p7P_MSC_CODON(sub_gm, k, c5);
  
      
      if      (esl_FCompare_old(MMX_SP(i,k), MMX_SP(i-3, k-1) + TSC(p7P_MM, k-1) + emit3, tol) == eslOK) { scur = p7T_M; c = 3; }
      else if (esl_FCompare_old(MMX_SP(i,k), IMX_SP(i-3, k-1) + TSC(p7P_IM, k-1) + emit3, tol) == eslOK) { scur = p7T_I; c = 3; }
      else if (esl_FCompare_old(MMX_SP(i,k), DMX_SP(i-3, k-1) + TSC(p7P_DM, k-1) + emit3, tol) == eslOK) { scur = p7T_D; c = 3; }     
      else if (esl_FCompare_old(MMX_SP(i,k), XMX(i-3, p7G_B)  + TSC(p7P_BM, k-1) + emit3, tol) == eslOK) { scur = p7T_B; c = 3; }     
      else if (esl_FCompare_old(MMX_SP(i,k), MMX_SP(i-1, k-1) + TSC(p7P_MM, k-1) + emit1, tol) == eslOK) { scur = p7T_M; c = 1; }
      else if (esl_FCompare_old(MMX_SP(i,k), IMX_SP(i-1, k-1) + TSC(p7P_IM, k-1) + emit1, tol) == eslOK) { scur = p7T_I; c = 1; }
      else if (esl_FCompare_old(MMX_SP(i,k), DMX_SP(i-1, k-1) + TSC(p7P_DM, k-1) + emit1, tol) == eslOK) { scur = p7T_D; c = 1; }
      else if (esl_FCompare_old(MMX_SP(i,k), XMX(i-1, p7G_B)  + TSC(p7P_BM, k-1) + emit1, tol) == eslOK) { scur = p7T_B; c = 1; }
      else if (esl_FCompare_old(MMX_SP(i,k), MMX_SP(i-2, k-1) + TSC(p7P_MM, k-1) + emit2, tol) == eslOK) { scur = p7T_M; c = 2; }
      else if (esl_FCompare_old(MMX_SP(i,k), IMX_SP(i-2, k-1) + TSC(p7P_IM, k-1) + emit2, tol) == eslOK) { scur = p7T_I; c = 2; }
      else if (esl_FCompare_old(MMX_SP(i,k), DMX_SP(i-2, k-1) + TSC(p7P_DM, k-1) + emit2, tol) == eslOK) { scur = p7T_D; c = 2; }
      else if (esl_FCompare_old(MMX_SP(i,k), XMX(i-2, p7G_B)  + TSC(p7P_BM, k-1) + emit2, tol) == eslOK) { scur = p7T_B; c = 2; }
      else if (esl_FCompare_old(MMX_SP(i,k), MMX_SP(i-4, k-1) + TSC(p7P_MM, k-1) + emit4, tol) == eslOK) { scur = p7T_M; c = 4; }
      else if (esl_FCompare_old(MMX_SP(i,k), IMX_SP(i-4, k-1) + TSC(p7P_IM, k-1) + emit4, tol) == eslOK) { scur = p7T_I; c = 4; }
      else if (esl_FCompare_old(MMX_SP(i,k), DMX_SP(i-4, k-1) + TSC(p7P_DM, k-1) + emit4, tol) == eslOK) { scur = p7T_D; c = 4; }
      else if (esl_FCompare_old(MMX_SP(i,k), XMX(i-4, p7G_B)  + TSC(p7P_BM, k-1) + emit4, tol) == eslOK) { scur = p7T_B; c = 4; }
      else if (esl_FCompare_old(MMX_SP(i,k), MMX_SP(i-5, k-1) + TSC(p7P_MM, k-1) + emit5, tol) == eslOK) { scur = p7T_M; c = 5; }
      else if (esl_FCompare_old(MMX_SP(i,k), IMX_SP(i-5, k-1) + TSC(p7P_IM, k-1) + emit5, tol) == eslOK) { scur = p7T_I; c = 5; }
      else if (esl_FCompare_old(MMX_SP(i,k), DMX_SP(i-5, k-1) + TSC(p7P_DM, k-1) + emit5, tol) == eslOK) { scur = p7T_D; c = 5; }
      else if (esl_FCompare_old(MMX_SP(i,k), XMX(i-5, p7G_B)  + TSC(p7P_BM, k-1) + emit5, tol) == eslOK) { scur = p7T_B; c = 5; }
      else ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
     
      /* pervious M state has to be appended here so we can know the codon length */
      if ((status = p7_trace_fs_Append(tr, sprv, k, i, c)) != eslOK) return status;
 
      k--; i-=c;
      break;

    case p7T_D:         /* D connects from M,D at i,k-1 */
      if (DMX_SP(i, k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);
      
      if      (esl_FCompare_old(DMX_SP(i,k), MMX_SP(i, k-1) + TSC(p7P_MD, k-1), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(DMX_SP(i,k), DMX_SP(i, k-1) + TSC(p7P_DD, k-1), tol) == eslOK) scur = p7T_D;
      else if (esl_FCompare_old(DMX_SP(i,k), XMX(i, p7G_B)  + TSC(p7P_BM, k-1), tol) == eslOK) scur = p7T_B;
      else ESL_EXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      k--;
      break;

    case p7T_I:         /* I connects from M,I at i-1,k*/
      if (IMX_SP(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);

      if      (esl_FCompare_old(IMX_SP(i,k), MMX_SP(i-3,k) + TSC(p7P_MI, k), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(IMX_SP(i,k), IMX_SP(i-3,k) + TSC(p7P_II, k), tol) == eslOK) scur = p7T_I;
      else ESL_EXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      i-=3;
      break;

    case p7T_N:         /* N connects from S, N */
      if (XMX(i, p7G_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = ( (i == 0) ? p7T_S : p7T_N);
      break;

    case p7T_B:         /* B connects from N, J */
      if (XMX(i,p7G_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

      if      (esl_FCompare_old(XMX(i,p7G_B), XMX(i, p7G_N) + sub_gm->xsc[p7P_N][p7P_MOVE], tol) == eslOK) scur = p7T_N;
      else  ESL_EXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over statetype[tpos-1] */

    if(scur != p7T_M) { if ((status = p7_trace_fs_Append(tr, scur, k, i, 0)) != eslOK) return status; }
  
    /* For NCJ, we had to defer i decrement. */
    if ( (scur == p7T_N || scur == p7T_C) && scur == sprv) i--; 
   
    sprv = scur;
  } /* end traceback, at S state */

  tr->M = sub_gm->M;
  tr->L = L;

  return p7_trace_fs_Reverse(tr);
}

