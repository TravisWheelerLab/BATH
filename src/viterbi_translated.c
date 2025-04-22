/* Frameshift aware Viterbi algorithm.*/
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

/* Function:  p7_fs_Viterbi()
 * Synopsis:  The Viterbi algorithm, with frameshift awareness.
 * Incept:    SRE, Tue Jan 30 10:50:53 2007 [Einstein's, St. Louis]
 * 
 * Purpose:   The Viterbi dynamic programming algorithm for aligning 
 *             DNA to proteins with frameshift awareness. 
 *
 *            Given a digital nucleotide sequence <dsq> of length <L>, 
 *            a frameshift aware codon profile <gm_fs>, and DP matrix 
 *            <gx> allocated for at least <L> by <gm_fs->M> cells; 
 *            calculate the maximum scoring path by fremashift-aware 
 *            Viterbi; return the Viterbi score in <ret_sc>, and the
 *            Viterbi matrix is in <gx>.
 *            
 *            The caller may then retrieve the Viterbi path by calling
 *            <p7_fs_VTrace()>.
 *           
 *            The Viterbi lod score is returned in nats. The caller
 *            needs to subtract a null model lod score, then convert
 *            to bits.
 *           
 * Args:      dsq    - sequence in digitized form, 1..L
 *            gcode  - genetic code for translation
 *            L      - length of dsq
 *            gm_fs  - profile. 
 *            gx     - DP matrix with room for an MxL alignment
 *            opt_sc - optRETURN: Viterbi lod score in nats
 *           
 * Return:   <eslOK> on success.
 */
int
p7_trans_Viterbi(const ESL_DSQ *dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx, float *opt_sc)
{
  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  int          M    = gm_fs->M;
  int          i,k;
  int          c3;
  int          v, w, x;
  float        esc  = p7_fs_profile_IsLocal(gm_fs) ? 0 : -eslINFINITY;

  /* Initialization of the zero row.  */
  XMX(0,p7G_N) = 0.;                                  /* S->N, p=1            */
  XMX(0,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE];         /* S->N->B, no N-tail   */
  XMX(0,p7G_E) = XMX(0,p7G_J) = XMX(0,p7G_C) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX(0,k) = IMX(0,k) = DMX(0,k) = -eslINFINITY;

  /*Special cases for the first 2 rows */
  v = w = x = -1;
  for(i = 1; i < 3; i++)
  {
    v = w;
    w = x;
  
    /* if new nucleotide is not A,C,G, or T set it to placeholder vlaue */
    if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[i])) x = dsq[i];
    else                                            x = p7P_MAXCODONS;

    for (k = 0; k <= M; k++) 
      MMX(i,k) = IMX(i,k) = DMX(i,k) = -eslINFINITY;
    
    XMX(i,p7G_E) = XMX(i,p7G_J) = XMX(i,p7G_C) = -eslINFINITY;

    XMX(i,p7G_N) = 0.;
    XMX(i,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE]; 
  } 

  /* Main DP recursion */
  for (i = 3; i <= L; i++) 
  {
    MMX(i,0) = IMX(i,0) = DMX(i,0) = -eslINFINITY;

    XMX(i, p7G_E) = -eslINFINITY;

    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder vlaue */
    if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[i])) x = dsq[i];
    else                                            x = p7P_MAXCODONS;

    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    for (k = 1; k < M; k++)
    {


      MMX(i,k) = ESL_MAX(MMX(i-3,k-1)   + TSC(p7P_MM,k-1),
                        ESL_MAX(IMX(i-3,k-1)   + TSC(p7P_IM,k-1),
                        ESL_MAX(DMX(i-3,k-1)   + TSC(p7P_DM,k-1),
                                XMX(i-3,p7G_B) + TSC(p7P_BM,k-1)))) + p7P_MSC_CODON(gm_fs, k, c3);

      IMX(i,k) = ESL_MAX(MMX(i-3,k) + TSC(p7P_MI,k),
                         IMX(i-3,k)  + TSC(p7P_II,k));

      DMX(i,k) = ESL_MAX(MMX(i,k-1) + TSC(p7P_MD,k-1),
                         DMX(i,k-1) + TSC(p7P_DD,k-1));

      XMX(i,p7G_E) = ESL_MAX(MMX(i,k) + esc, XMX(i,p7G_E));
    }

    MMX(i,M) = ESL_MAX(MMX(i-3,M-1)   + TSC(p7P_MM,M-1),
               ESL_MAX(IMX(i-3,M-1)   + TSC(p7P_IM,M-1),
               ESL_MAX(DMX(i-3,M-1)   + TSC(p7P_DM,M-1),
                       XMX(i-3,p7G_B) + TSC(p7P_BM,M-1)))) + p7P_MSC_CODON(gm_fs, M, c3);

    IMX(i,M) = -eslINFINITY;

    DMX(i,M) = ESL_MAX(MMX(i,M-1) + TSC(p7P_MD,M-1),
                       DMX(i,M-1) + TSC(p7P_DD,M-1));

    XMX(i,p7G_E) = ESL_MAX(MMX(i,M), XMX(i,p7G_E));

    XMX(i,p7G_J) = ESL_MAX(XMX(i-1,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP],
                           XMX(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_LOOP]);

    XMX(i,p7G_C) = ESL_MAX(XMX(i-1,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                           XMX(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);

    XMX(i,p7G_N) =         XMX(i-1,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];

    XMX(i,p7G_B) = ESL_MAX(XMX(i,p7G_N) + gm_fs->xsc[p7P_N][p7P_MOVE],
                           XMX(i,p7G_J) + gm_fs->xsc[p7P_J][p7P_MOVE]);
  }
  
  /* T state (not stored) */
  if (opt_sc != NULL) *opt_sc = XMX(L,p7G_C) + gm_fs->xsc[p7P_C][p7P_MOVE];

  gx->M = gm_fs->M;
  gx->L = L;

  return eslOK;

}





int
p7_trans_VTrace(const ESL_DSQ *dsq, int L, const ESL_GENCODE *gcode, const P7_FS_PROFILE *gm_fs, const P7_GMX *gx, P7_TRACE *tr)
{
  int          i   = L;     /* position in seq (1..L)         */
  int          k   = 0;     /* position in model (1..M)       */
  int          M   = gm_fs->M;
  int          v,w,x;
  int          c, c3;
  int          sprv,scur;
  float        emit;
  float      **dp  = gx->dp;    /* so {MDI}MX() macros work       */
  float       *xmx = gx->xmx;   /* so XMX() macro works           */
  float        tol = 1e-5;  /* floating point "equality" test */
  float const *tsc = gm_fs->tsc; 
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
   
      if      (esl_FCompare_old(XMX(i, p7G_C), XMX(i-1, p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP], tol) == eslOK)  scur = p7T_C;
      else if (esl_FCompare_old(XMX(i, p7G_C), XMX(i,   p7G_E) + gm_fs->xsc[p7P_E][p7P_MOVE], tol) == eslOK)  scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
 
      break;

    case p7T_E:     /* E connects from any M state. k set here */
      if (XMX(i, p7G_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);

      scur = p7T_M;     /* can't come from D, in a *local* Viterbi trace. */
      for (k = M; k >= 1; k--) if (esl_FCompare_old(XMX(i, p7G_E), MMX(i,k), tol) == eslOK) break;
      if (k == 0) ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      break;

    case p7T_M:         /* M connects from i-1,k-1, or B */
      if (MMX(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);
     
      if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[i-2])) v = dsq[i-2];
      else                                              v = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[i-1])) w = dsq[i-1];
      else                                              w = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[i]))   x = dsq[i];
      else                                              x = p7P_MAXCODONS; 

      c3 = p7P_CODON3(v, w, x);

      emit = p7P_MSC_CODON(gm_fs, k, c3);

      if      (esl_FCompare_old(MMX(i,k), MMX(i-3, k-1)   + TSC(p7P_MM, k-1) + emit, tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(MMX(i,k), IMX(i-3, k-1)   + TSC(p7P_IM, k-1) + emit, tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare_old(MMX(i,k), DMX(i-3, k-1)   + TSC(p7P_DM, k-1) + emit, tol) == eslOK) scur = p7T_D;      
      else if (esl_FCompare_old(MMX(i,k), XMX(i-3, p7G_B) + TSC(p7P_BM, k-1) + emit, tol) == eslOK) scur = p7T_B;      
      else ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
      
      k--; i-=3;
      break;

    case p7T_D:         /* D connects from M,D at i,k-1 */
      if (DMX(i, k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);

      if      (esl_FCompare_old(DMX(i,k), MMX(i, k-1) + TSC(p7P_MD, k-1), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(DMX(i,k), DMX(i, k-1) + TSC(p7P_DD, k-1), tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      k--;
      break;

    case p7T_I:         /* I connects from M,I at i-1,k*/
      if (IMX(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);

      if      (esl_FCompare_old(IMX(i,k), MMX(i-3,k) + TSC(p7P_MI, k), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(IMX(i,k), IMX(i-3,k) + TSC(p7P_II, k), tol) == eslOK) scur = p7T_I;
      else ESL_EXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      i-=3;
      break;

    case p7T_N:         /* N connects from S, N */
      if (XMX(i, p7G_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = ( (i == 0) ? p7T_S : p7T_N);
      break;

    case p7T_B:         /* B connects from N, J */
      if (XMX(i,p7G_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

      if      (esl_FCompare_old(XMX(i,p7G_B), XMX(i, p7G_N) + gm_fs->xsc[p7P_N][p7P_MOVE], tol) == eslOK) scur = p7T_N;
      else if (esl_FCompare_old(XMX(i,p7G_B), XMX(i, p7G_J) + gm_fs->xsc[p7P_J][p7P_MOVE], tol) == eslOK) scur = p7T_J;
      else  ESL_EXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    case p7T_J:         /* J connects from E(i) or J(i-1) */
      if (XMX(i,p7G_J) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible J reached at i=%d", i);
    
      if      (esl_FCompare_old(XMX(i,p7G_J), XMX(i-1,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP], tol) == eslOK) scur = p7T_J;
      else if (esl_FCompare_old(XMX(i,p7G_J), XMX(i,  p7G_E) + gm_fs->xsc[p7P_E][p7P_LOOP], tol) == eslOK) scur = p7T_E;
      else  ESL_EXCEPTION(eslFAIL, "J at i=%d couldn't be traced", i);
      
      break;

    default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over statetype[tpos-1] */

    c = (scur == p7T_M ? 3:0);
    if ((status = p7_trace_fs_Append(tr, scur, k, i, c)) != eslOK) return status;
  
    /* For NCJ, we had to defer i decrement. */
    if ( (scur == p7T_N || scur == p7T_J || scur == p7T_C) && scur == sprv) i--; 
    
    sprv = scur;
  } /* end traceback, at S state */

  tr->M = gm_fs->M;
  tr->L = L;

  return p7_trace_fs_Reverse(tr);
}


