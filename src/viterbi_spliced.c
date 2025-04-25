/* Frameshift aware Viterbi algorithm.*/
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

/* Function:  p7_sp_trans_semiglobal_Viterbi
 * Synopsis:  The Viterbi algorithm for aligning codons to aminos, 
 *            global with regards to the hmm, with splice states 
 * 
 * Purpose:   The Viterbi dynamic programming algorithm for aligning 
 *            DNA to proteins forcing all HMM positions to be used 
 *            and using R and P states to allow introns. 
 *
 *            This specialized version of Viterbi is for finding small 
 *            or otherwise hard to find exons between two existing hits.
 *
 *            Given a digital nucleotide sub-sequence <sub_dsq> of 
 *            length <L>, a codon sub-profile <sub_gm>, and DP matrix 
 *            <gx> allocated for at least <L> by <sub_gm->M> cells; 
 *            calculate the maximum scoring path by semi-global spliced
 *            Viterbi; return the Viterbi matrix is in <gx>.
 *            
 *            The caller may then retrieve the Viterbi path by calling
 *            <p7_sp_trans_semiglobal_VTrace()>.
 *           
 * Args:      sub_dsq - sub-sequence in digitized form, 1..L
 *            gcode   - genetic code for translation
 *            L       - length of dsq
 *            sub_gm  - sub-profile 
 *            gx      - DP matrix with room for an MxL alignment
 *           
 * Return:   <eslOK> on success.
 */
int
p7_sp_trans_semiglobal_Viterbi(const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *sub_gm, P7_GMX *gx)
{
  float const *tsc  = sub_gm->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  int          M    = sub_gm->M;
  int          i,k;
  int          c3;
  int          v, w, x;

  /* Initialization of the zero row.  */
  XMX(0,p7G_N) = 0.;                                  /* S->N, p=1            */
  XMX(0,p7G_B) = sub_gm->xsc[p7P_N][p7P_MOVE];            /* S->N->B, no N-tail   */ 
  XMX(0,p7G_E) = XMX(0,p7G_J) = XMX(0,p7G_C) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_SP(0,k) = IMX_SP(0,k) = RMX_SP(0,k) = PMX_SP(0,k) = -eslINFINITY;

  DMX_SP(0,0) = -eslINFINITY;
  DMX_SP(0,1) = XMX(0,p7G_B) + TSC(p7P_BM,0); 
  for (k = 2; k <= M; k++)
    DMX_SP(0,k) = DMX_SP(0,k-1) + TSC(p7P_DD,k-1);

  /*Special cases for the first 2 rows */
  v = w = x = -1;
  for(i = 1; i <= 2; i++)
  {
    v = w;
    w = x;
  
    /* if new nucleotide is not A,C,G, or T set it to placeholder vlaue */
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                            x = p7P_MAXCODONS;

    XMX(i,p7G_N) =  XMX(i-1,p7G_N) + sub_gm->xsc[p7P_N][p7P_LOOP];
    XMX(i,p7G_B) =  XMX(i,p7G_N)   + sub_gm->xsc[p7P_N][p7P_MOVE]; 

    for (k = 0; k <= M; k++) 
      MMX_SP(i,k) = IMX_SP(i,k) = RMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;

    DMX_SP(i,0) = -eslINFINITY;
    DMX_SP(i,1) = XMX(i,p7G_B) + TSC(p7P_BM,0);
    for (k = 2; k <= M; k++)
      DMX_SP(i,k) = DMX_SP(i,k-1) + TSC(p7P_DD,k-1);        

    XMX(i,p7G_E) = XMX(i,p7G_J) = XMX(i,p7G_C) = -eslINFINITY;

  } 

  /* R & P state create special cases for rows 3-11 */
  for(i = 3; i <= ESL_MIN(11,L); i++)
  {
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder vlaue */
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                                x = p7P_MAXCODONS;
    
    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    XMX(i,p7G_N) = XMX(i-1,p7G_N) + sub_gm->xsc[p7P_N][p7P_LOOP];
    XMX(i,p7G_B) = XMX(i,p7G_N)   + sub_gm->xsc[p7P_N][p7P_MOVE];

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = RMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;
  
    MMX_SP(i,1) = XMX(i-3,p7G_B) + TSC(p7P_BM,0) + p7P_MSC_CODON(sub_gm, 1, c3);  //B->M
    
    IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,1),
                          IMX_SP(i-3,1) + TSC(p7P_II,1));

    DMX_SP(i,1) = XMX(i,p7G_B) + TSC(p7P_BM,0);
    
	RMX_SP(i,1) = PMX_SP(i,1) = -eslINFINITY;

    for (k = 2; k < M; k++) {
      
      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,k-1),
                            DMX_SP(i-3,k-1) + TSC(p7P_DM,k-1))) + p7P_MSC_CODON(sub_gm, k, c3);

      IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,k),
                            IMX_SP(i-3,k) + TSC(p7P_II,k));

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,k-1));

      RMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;
    }        
 
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,M-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,M-1),
                          DMX_SP(i-3,M-1) + TSC(p7P_DM,M-1))) + p7P_MSC_CODON(sub_gm, M, c3);

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,M-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,M-1));

    RMX_SP(i,M) = PMX_SP(i,M) = -eslINFINITY; 
 
    XMX(i,p7G_E) = ESL_MAX(MMX_SP(i,M), DMX_SP(i,M));
    XMX(i,p7G_C) = ESL_MAX(XMX(i-1,p7G_C) + sub_gm->xsc[p7P_C][p7P_LOOP],
                           XMX(i,p7G_E)   + sub_gm->xsc[p7P_E][p7P_MOVE]);
    
    XMX(i,p7G_J) = -eslINFINITY;    

  }
  
  /* Main DP recursion */
  for (i = 12; i <= L; i++) 
  {
    v = w;
    w = x;

    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                            x = p7P_MAXCODONS;

    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    XMX(i,p7G_N) = XMX(i-1,p7G_N) + sub_gm->xsc[p7P_N][p7P_LOOP];
    XMX(i,p7G_B) = XMX(i,p7G_N)   + sub_gm->xsc[p7P_N][p7P_MOVE];

    
    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = RMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    MMX_SP(i,1) = XMX(i-3,p7G_B) + TSC(p7P_BM,0) + p7P_MSC_CODON(sub_gm, 1, c3); //B->M

    IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,1),
                          IMX_SP(i-3,1) + TSC(p7P_II,1));

    DMX_SP(i,1) = XMX(i,p7G_B) + TSC(p7P_BM,0);
    
    RMX_SP(i,1) = PMX_SP(i,1) = -eslINFINITY;

    for (k = 2; k < M; k++) {

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,k-1),
                    ESL_MAX(DMX_SP(i-3,k-1) + TSC(p7P_DM,k-1),
                            PMX_SP(i-3,k-1) + TSC(p7P_BM,k-1)))) + p7P_MSC_CODON(sub_gm, k, c3);

      IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,k),
                            IMX_SP(i-3,k) + TSC(p7P_II,k));

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,k-1),
                    ESL_MAX(DMX_SP(i,k-1) + TSC(p7P_DD,k-1),
                            PMX_SP(i,k-1) + TSC(p7P_BM,k-1)));

      RMX_SP(i,k) = ESL_MAX(MMX_SP(i-12,k), DMX_SP(i-12,k)) + (sub_gm->xsc[p7P_J][p7P_LOOP] * 12);
      PMX_SP(i,k) = ESL_MAX(RMX_SP(i-1,k), PMX_SP(i-1,k))   + sub_gm->xsc[p7P_J][p7P_LOOP];
    }

    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,M-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,M-1),
                  ESL_MAX(DMX_SP(i-3,M-1) + TSC(p7P_DM,M-1),
                          PMX_SP(i-3,M-1) + TSC(p7P_BM,M-1)))) + p7P_MSC_CODON(sub_gm, M, c3);

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,M-1),
                  ESL_MAX(DMX_SP(i,M-1) + TSC(p7P_DD,M-1),
                          PMX_SP(i,M-1) + TSC(p7P_BM,M-1)));

    RMX_SP(i,M) =  PMX_SP(i,M) = -eslINFINITY;

    XMX(i,p7G_E) = ESL_MAX(MMX_SP(i,M), DMX_SP(i,M));
    XMX(i,p7G_C) = ESL_MAX(XMX(i-1,p7G_C) + sub_gm->xsc[p7P_C][p7P_LOOP],
                           XMX(i,p7G_E)   + sub_gm->xsc[p7P_E][p7P_MOVE]);

    XMX(i,p7G_J) = -eslINFINITY;    

  }
  
  gx->M = sub_gm->M;
  gx->L = L;

  //p7_gmx_sp_Dump(stdout, gx, p7_DEFAULT);
  return eslOK;

}

int
p7_sp_trans_local_Viterbi(const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *sub_gm, P7_GMX *gx)
{

  float const *tsc  = sub_gm->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  int          M    = sub_gm->M;
  int          i,k;
  int          c3;
  int          v, w, x;

  /* Initialization of the zero row.  */
  XMX(0,p7G_N) = 0;                                           /* S->N, p=1            */
  XMX(0,p7G_B) = sub_gm->xsc[p7P_N][p7P_MOVE];                    /* S->N->B, no N-tail   */
  XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) = -eslINFINITY;  /* need seq to get here */
  for (k = 0; k <= sub_gm->M; k++)
    MMX_SP(0,k) = IMX_SP(0,k) = DMX_SP(0,k) = RMX_SP(0,k) = PMX_SP(0,k) = -eslINFINITY;            /* need seq to get here */

 /*Special cases for the first 2 rows */
  v = w = x = -1;
  for(i = 1; i <= 2; i++)
  {
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder vlaue */
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                            x = p7P_MAXCODONS;

    XMX(i,p7G_N) =  XMX(i-1,p7G_N) + sub_gm->xsc[p7P_N][p7P_LOOP];
    XMX(i,p7G_B) =  XMX(i,p7G_N)   + sub_gm->xsc[p7P_N][p7P_MOVE];

    for (k = 0; k <= M; k++)
      MMX_SP(i,k) = IMX_SP(i,k) = DMX_SP(i,k) = RMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;

    XMX(i,p7G_E) = XMX(i,p7G_J) = XMX(i,p7G_C) = -eslINFINITY;

  } 

  /* R & P state create special cases for rows 3-11 */
  for(i = 3; i <= ESL_MIN(11,L); i++)
  {
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder vlaue */
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                                x = p7P_MAXCODONS;

    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    XMX(i,p7G_N) = XMX(i-1,p7G_N) + sub_gm->xsc[p7P_N][p7P_LOOP];
    XMX(i,p7G_B) = XMX(i,p7G_N)   + sub_gm->xsc[p7P_N][p7P_MOVE];

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = RMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    XMX(i,p7G_E) = -eslINFINITY;
    for (k = 1; k <= M; k++) {

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,k-1),
                    ESL_MAX(DMX_SP(i-3,k-1) + TSC(p7P_DM,k-1),
                            XMX(i-3,p7G_B)  + TSC(p7P_BM,k-1)))) + p7P_MSC_CODON(sub_gm, k, c3);

      IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,k),
                            IMX_SP(i-3,k) + TSC(p7P_II,k));

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,k-1));

      RMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;

      XMX(i,p7G_E) = ESL_MAX(MMX_SP(i,k), XMX(i,p7G_E)); 
    }


    XMX(i,p7G_C) = ESL_MAX(XMX(i-1,p7G_C) + sub_gm->xsc[p7P_C][p7P_LOOP],
                           XMX(i,p7G_E)   + sub_gm->xsc[p7P_E][p7P_MOVE]);

    XMX(i,p7G_J) = -eslINFINITY;

  }
 
  /* Main DP recursion */
  for (i = 12; i <= L; i++)
  {
    v = w;
    w = x;

    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                            x = p7P_MAXCODONS;

    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    XMX(i,p7G_N) = XMX(i-1,p7G_N) + sub_gm->xsc[p7P_N][p7P_LOOP];
    XMX(i,p7G_B) = XMX(i,p7G_N)   + sub_gm->xsc[p7P_N][p7P_MOVE];


    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = RMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    XMX(i,p7G_E) = -eslINFINITY;
    for (k = 1; k <= M; k++) {

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,k-1),
                    ESL_MAX(DMX_SP(i-3,k-1) + TSC(p7P_DM,k-1),
                    ESL_MAX(PMX_SP(i-3,k-1) + TSC(p7P_BM,k-1),
                            XMX(i-3,p7G_B)  + TSC(p7P_BM,k-1))))) + p7P_MSC_CODON(sub_gm, k, c3);

      IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,k),
                            IMX_SP(i-3,k) + TSC(p7P_II,k));

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,k-1),
                    ESL_MAX(DMX_SP(i,k-1) + TSC(p7P_DD,k-1),
                            PMX_SP(i,k-1) + TSC(p7P_BM,k-1)));

      RMX_SP(i,k) = ESL_MAX(MMX_SP(i-12,k), DMX_SP(i-12,k)) + (sub_gm->xsc[p7P_J][p7P_LOOP] * 12);
      PMX_SP(i,k) = ESL_MAX(RMX_SP(i-1,k), PMX_SP(i-1,k))   + sub_gm->xsc[p7P_J][p7P_LOOP];

    
    }

    XMX(i,p7G_C) = ESL_MAX(XMX(i-1,p7G_C) + sub_gm->xsc[p7P_C][p7P_LOOP],
                           XMX(i,p7G_E)   + sub_gm->xsc[p7P_E][p7P_MOVE]);

    XMX(i,p7G_J) = -eslINFINITY;

  }
  
  gx->M = sub_gm->M;
  gx->L = L;
  return eslOK;
}



int
p7_sp_trans_semiglobal_VTrace(const ESL_DSQ *sub_dsq, int L, const ESL_GENCODE *gcode, const P7_FS_PROFILE *sub_gm, const P7_GMX *gx, P7_TRACE *tr)
{
  int          i   = L;     /* position in seq (1..L)         */
  int          k   = 0;     /* position in model (1..M)       */
  int          M   = sub_gm->M;
  int          v,w,x;
  int          c, c3;
  int          sprv,scur;
  float        emit;
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
   
      if      (esl_FCompare_old(XMX(i, p7G_C), XMX(i-1, p7G_C) + sub_gm->xsc[p7P_C][p7P_LOOP], tol) == eslOK)  scur = p7T_C;
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
     
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i-2])) v = sub_dsq[i-2];
      else                                              v = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i-1])) w = sub_dsq[i-1];
      else                                              w = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i]))   x = sub_dsq[i];
      else                                              x = p7P_MAXCODONS; 

      c3 = p7P_CODON3(v, w, x);

      emit = p7P_MSC_CODON(sub_gm, k, c3);

      if      (esl_FCompare_old(MMX_SP(i,k), MMX_SP(i-3, k-1) + TSC(p7P_MM, k-1) + emit, tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(MMX_SP(i,k), IMX_SP(i-3, k-1) + TSC(p7P_IM, k-1) + emit, tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare_old(MMX_SP(i,k), DMX_SP(i-3, k-1) + TSC(p7P_DM, k-1) + emit, tol) == eslOK) scur = p7T_D;      
      else if (esl_FCompare_old(MMX_SP(i,k), PMX_SP(i-3, k-1) + TSC(p7P_BM, k-1) + emit, tol) == eslOK) scur = p7T_P;
      else if (esl_FCompare_old(MMX_SP(i,k), XMX(i-3, p7G_B)  + TSC(p7P_BM, k-1) + emit, tol) == eslOK) scur = p7T_B;      
      else ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
      
      k--; i-=3;
      break;

    case p7T_D:         /* D connects from M,D at i,k-1 */
      if (DMX_SP(i, k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);
      
      if      (esl_FCompare_old(DMX_SP(i,k), MMX_SP(i, k-1) + TSC(p7P_MD, k-1), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(DMX_SP(i,k), DMX_SP(i, k-1) + TSC(p7P_DD, k-1), tol) == eslOK) scur = p7T_D;
      else if (esl_FCompare_old(DMX_SP(i,k), PMX_SP(i, k-1) + TSC(p7P_BM, k-1), tol) == eslOK) scur = p7T_P;
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

    case p7T_P:
      if (PMX_SP(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible P reached at k=%d,i=%d", k,i);

      if      (esl_FCompare_old(PMX_SP(i,k), PMX_SP(i-1, k) + sub_gm->xsc[p7P_J][p7P_LOOP], tol) == eslOK) scur = p7T_P;
      else if (esl_FCompare_old(PMX_SP(i,k), RMX_SP(i-1, k) + sub_gm->xsc[p7P_J][p7P_LOOP], tol) == eslOK) scur = p7T_R;
      else ESL_EXCEPTION(eslFAIL, "P at k=%d,i=%d couldn't be traced", k,i);  

      i--;
      break;

    case p7T_R:
      if (RMX_SP(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible R reached at k=%d,i=%d", k,i); 

      if      (esl_FCompare_old(RMX_SP(i,k), MMX_SP(i-12, k) + (sub_gm->xsc[p7P_J][p7P_LOOP]*12), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(RMX_SP(i,k), DMX_SP(i-12, k) + (sub_gm->xsc[p7P_J][p7P_LOOP]*12), tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "R at k=%d,i=%d couldn't be traced", k,i);

      i-=12;
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

    c = (scur == p7T_M ? 3:0);
    if ((status = p7_trace_fs_Append(tr, scur, k, i, c)) != eslOK) return status;
  
    /* For NCJ, we had to defer i decrement. */
    if ( (scur == p7T_N || scur == p7T_J || scur == p7T_C) && scur == sprv) i--; 
   
    sprv = scur;
  } /* end traceback, at S state */

  tr->M = sub_gm->M;
  tr->L = L;

  return p7_trace_fs_Reverse(tr);
}

int
p7_sp_trans_local_VTrace(const ESL_DSQ *sub_dsq, int L, const ESL_GENCODE *gcode, const P7_FS_PROFILE *sub_gm, const P7_GMX *gx, P7_TRACE *tr)
{
  int          i   = L;     /* position in seq (1..L)         */
  int          k   = 0;     /* position in model (1..M)       */
  int          M   = sub_gm->M;
  int          v,w,x;
  int          c, c3;
  int          sprv,scur;
  float        emit;
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
   
      if      (esl_FCompare_old(XMX(i, p7G_C), XMX(i-1, p7G_C) + sub_gm->xsc[p7P_C][p7P_LOOP], tol) == eslOK)  scur = p7T_C;
      else if (esl_FCompare_old(XMX(i, p7G_C), XMX(i,   p7G_E) + sub_gm->xsc[p7P_E][p7P_MOVE], tol) == eslOK)  scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
 
      break;

    case p7T_E:     /* E connects from any M state. k set here */
      if (XMX(i, p7G_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
  
      for (k = M; k >= 1; k--) if (esl_FCompare_old(XMX(i, p7G_E), MMX_SP(i,k), tol) == eslOK) break;
      if (k == 0) ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      scur = p7T_M; 
      break;

    case p7T_M:         /* M connects from i-1,k-1, or B */
      if (MMX_SP(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);
     
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i-2])) v = sub_dsq[i-2];
      else                                              v = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i-1])) w = sub_dsq[i-1];
      else                                              w = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i]))   x = sub_dsq[i];
      else                                              x = p7P_MAXCODONS; 

      c3 = p7P_CODON3(v, w, x);

      emit = p7P_MSC_CODON(sub_gm, k, c3);

      if      (esl_FCompare_old(MMX_SP(i,k), MMX_SP(i-3, k-1) + TSC(p7P_MM, k-1) + emit, tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(MMX_SP(i,k), IMX_SP(i-3, k-1) + TSC(p7P_IM, k-1) + emit, tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare_old(MMX_SP(i,k), DMX_SP(i-3, k-1) + TSC(p7P_DM, k-1) + emit, tol) == eslOK) scur = p7T_D;      
      else if (esl_FCompare_old(MMX_SP(i,k), PMX_SP(i-3, k-1) + TSC(p7P_BM, k-1) + emit, tol) == eslOK) scur = p7T_P;
      else if (esl_FCompare_old(MMX_SP(i,k), XMX(i-3, p7G_B)  + TSC(p7P_BM, k-1) + emit, tol) == eslOK) scur = p7T_B;      
      else ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
      
      k--; i-=3;
      break;

    case p7T_D:         /* D connects from M,D at i,k-1 */
      if (DMX_SP(i, k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);
      
      if      (esl_FCompare_old(DMX_SP(i,k), MMX_SP(i, k-1) + TSC(p7P_MD, k-1), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(DMX_SP(i,k), DMX_SP(i, k-1) + TSC(p7P_DD, k-1), tol) == eslOK) scur = p7T_D;
      else if (esl_FCompare_old(DMX_SP(i,k), PMX_SP(i, k-1) + TSC(p7P_BM, k-1), tol) == eslOK) scur = p7T_P;
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

    case p7T_P:
      if (PMX_SP(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible P reached at k=%d,i=%d", k,i);

      if      (esl_FCompare_old(PMX_SP(i,k), PMX_SP(i-1, k) + sub_gm->xsc[p7P_J][p7P_LOOP], tol) == eslOK) scur = p7T_P;
      else if (esl_FCompare_old(PMX_SP(i,k), RMX_SP(i-1, k) + sub_gm->xsc[p7P_J][p7P_LOOP], tol) == eslOK) scur = p7T_R;
      else ESL_EXCEPTION(eslFAIL, "P at k=%d,i=%d couldn't be traced", k,i);  

      i--;
      break;

    case p7T_R:
      if (RMX_SP(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible R reached at k=%d,i=%d", k,i); 

      if      (esl_FCompare_old(RMX_SP(i,k), MMX_SP(i-12, k) + (sub_gm->xsc[p7P_J][p7P_LOOP]*12), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(RMX_SP(i,k), DMX_SP(i-12, k) + (sub_gm->xsc[p7P_J][p7P_LOOP]*12), tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "R at k=%d,i=%d couldn't be traced", k,i);

      i-=12;
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

    c = (scur == p7T_M ? 3:0);
    if ((status = p7_trace_fs_Append(tr, scur, k, i, c)) != eslOK) return status;
  
    /* For NCJ, we had to defer i decrement. */
    if ( (scur == p7T_N || scur == p7T_J || scur == p7T_C) && scur == sprv) i--; 

    sprv = scur;
  } /* end traceback, at S state */

  tr->M = sub_gm->M;
  tr->L = L;

  return p7_trace_fs_Reverse(tr);
}


