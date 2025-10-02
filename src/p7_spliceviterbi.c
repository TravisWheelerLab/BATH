/* Spliced Viterbi algorithm and trace back - translated and fraemshift versions. */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_splice.h"

#define IVX(i,k,c) (iv[((k)*p7P_CODONS)+L+3-(i)+(c)])
#define TSC_R -14
#define TSC_P -14


int
p7_spliceviterbi_translated_semiglobal(SPLICE_SITE_IDX *splice_signals, float *signal_scores, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx)
{
  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float        P_tmp;
  int          M    = gm_fs->M;
  int          i,k,s;
  int          c3, c2, c1;
  int          t, u, v, w, x, y, z;
  int          donor_idx;
  int          AGXXX, ACXXX;
  int          AGXX,  ACXX;
  int          AGX,   ACX;
  
  
  /* Initialization of the zero row.  */
  XMX(0,p7G_N) = 0.;                                  /* S->N, p=1            */
  XMX(0,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE];            /* S->N->B, no N-tail   */ 
  XMX(0,p7G_E) = XMX(0,p7G_J) = XMX(0,p7G_C) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_SP(0,k) = IMX_SP(0,k) = PMX_SP(0,k) = -eslINFINITY;

  DMX_SP(0,0) = -eslINFINITY;
  DMX_SP(0,1) = XMX(0,p7G_B) + TSC(p7P_BM,0); 
  for (k = 2; k <= M; k++)
    DMX_SP(0,k) = DMX_SP(0,k-1) + TSC(p7P_DD,k-1);

  /*Special cases for the first 2 rows */
  t = u = v = w = x = -1;
  for(i = 1; i <= 2; i++)
  {
    w = x;
  
    /* if new nucleotide is not A,C,G, or T set it to placeholder vlaue */
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                                x = p7P_MAXCODONS;

    XMX(i,p7G_N) =  0; 
    XMX(i,p7G_B) =  gm_fs->xsc[p7P_N][p7P_MOVE]; 

    for (k = 0; k <= M; k++) 
      MMX_SP(i,k) = IMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;

    DMX_SP(i,0) = -eslINFINITY;
    DMX_SP(i,1) = XMX(i,p7G_B) + TSC(p7P_BM,0);
    for (k = 2; k <= M; k++)
      DMX_SP(i,k) = DMX_SP(i,k-1) + TSC(p7P_DD,k-1);        

    XMX(i,p7G_E) = XMX(i,p7G_J) = XMX(i,p7G_C) = -eslINFINITY;

  } 

  AGXXX = ACXXX = AGXX = ACXX = AGX = ACX = FALSE; 
  /* Main DP recursion */
  for (i = 3; i <= L; i++) 
  {
    AGXXX = AGXX;
    ACXXX = ACXX;
    AGXX  = AGX;
    ACXX  = ACX;

    t = u;
    u = v;
    v = w;
    w = x;

    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                            x = p7P_MAXCODONS;

    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    /* Check if we at at a potential donor site and record index*/
    if( w == 2 ) {
      if( x == 3) { //GT
        p7_splicesiteidx_Grow(splice_signals, L, p7S_GTAG);
        splice_signals->index[p7S_GTAG][splice_signals->N[p7S_GTAG]] = i;
        splice_signals->N[p7S_GTAG]++;
      }
      else if( x == 2) { //GC
        p7_splicesiteidx_Grow(splice_signals, L, p7S_GCAG);
        splice_signals->index[p7S_GCAG][splice_signals->N[p7S_GCAG]] = i;
        splice_signals->N[p7S_GCAG]++;
      }
    }
    else if( w == 0 && x == 3 ) { //AT
       p7_splicesiteidx_Grow(splice_signals, L, p7S_ATAC);
       splice_signals->index[p7S_ATAC][splice_signals->N[p7S_ATAC]] = i;
       splice_signals->N[p7S_ATAC]++;
    } 

    /* Check if we at a potential acceptor site */
    AGX = ACX = FALSE;
    if ( v == 0  && w == 2) AGX = TRUE;
    else                    AGX = FALSE;
   
    if(  v == 0  && w == 1) ACX = TRUE;
    else                    ACX = FALSE;

    XMX(i,p7G_N) = XMX(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];
    XMX(i,p7G_B) = XMX(i,p7G_N)   + gm_fs->xsc[p7P_N][p7P_MOVE];

    
    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    MMX_SP(i,1) = XMX(i-3,p7G_B) + TSC(p7P_BM,0) + p7P_MSC_CODON(gm_fs, 1, c3); //B->M

    if(p7P_MSC_CODON(gm_fs, 1, c3) == -eslINFINITY) IMX_SP(i,1) = -eslINFINITY;
    else                                             IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,1),
                                                     IMX_SP(i-3,1) + TSC(p7P_II,1));

    DMX_SP(i,1) = XMX(i,p7G_B) + TSC(p7P_BM,0);
    
    PMX_SP(i,1) = -eslINFINITY;

    for (k = 2; k < M; k++) {

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,k-1),
                    ESL_MAX(DMX_SP(i-3,k-1) + TSC(p7P_DM,k-1),
                            PMX_SP(i-3,k-1) + TSC_P)))         + p7P_MSC_CODON(gm_fs, k, c3);

      if(p7P_MSC_CODON(gm_fs, k, c3) == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;
      else                                             IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,k),
                                                                             IMX_SP(i-3,k) + TSC(p7P_II,k));

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,k-1));

      PMX_SP(i,k) = -eslINFINITY;
      /*Check for acceptor sites and then match to correct donor sites to create splice-codons */
      if ( AGXXX ) {
         P_tmp = ESL_MAX(MMX_SP(donor_idx-2,k-1) , DMX_SP(donor_idx-2,k-1));
         /* Check all possible GTAGXXX */
         for(s = 0; s < splice_signals->N[p7S_GTAG]; s++) {
           donor_idx = splice_signals->index[p7S_GTAG][s];
           if(i - (donor_idx-2) > MIN_INTRON_LENG) {
             PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k),
                                   P_tmp + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c3));
           }
         }
           
         /* Check all possible GCAGXXX */
         for(s = 0; s < splice_signals->N[p7S_GCAG]; s++) {
           donor_idx = splice_signals->index[p7S_GCAG][s];
           if(i - (donor_idx-2) > MIN_INTRON_LENG) {
             PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k),
                                   P_tmp + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c3));
           }
         }
      }
      else if ( ACXXX ) {
        P_tmp = ESL_MAX(MMX_SP(donor_idx-2,k-1) , DMX_SP(donor_idx-2,k-1)); 
        /* Check all possible ATACXXX */
        for(s = 0; s < splice_signals->N[p7S_ATAC]; s++) {
          donor_idx = splice_signals->index[p7S_ATAC][s];
 
          if(i - (donor_idx-2) > MIN_INTRON_LENG) {
            PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k),
                          P_tmp + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c3));
          }
        }
      }
      
      if ( AGXX ) {
         P_tmp = ESL_MAX(MMX_SP(donor_idx-3,k-1) , DMX_SP(donor_idx-3,k-1)); 
         /* Check all possible XGTAGXX */
         for(s = 0; s < splice_signals->N[p7S_GTAG]; s++) {
           donor_idx = splice_signals->index[p7S_GTAG][s];
           if(i - (donor_idx-3) > MIN_INTRON_LENG) {
             /* Create codon with one nuc from before the donor and two from after the acceptor */
             if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-2])) z = sub_dsq[donor_idx-2];
             else                                                          z = p7P_MAXCODONS;

             c2 = p7P_CODON3(z, w, x);
             c2 = p7P_MINIDX(c2, p7P_DEGEN_C);  
          
             PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k),
                                   P_tmp + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c2));
           }
         }
        
         /* Check all possible XGCAGXX */
         for(s = 0; s < splice_signals->N[p7S_GCAG]; s++) {
           donor_idx = splice_signals->index[p7S_GCAG][s];
           if(i - (donor_idx-3) > MIN_INTRON_LENG) {
             /* Create codon with one nuc from before the donor and two from after the acceptor */
             if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-2])) z = sub_dsq[donor_idx-2];
             else                                                          z = p7P_MAXCODONS;

             c2 = p7P_CODON3(z, w, x);
             c2 = p7P_MINIDX(c2, p7P_DEGEN_C);

             PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k),
                                   P_tmp + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c2));
           }
         }
      }
      else if (ACXX ) {
        P_tmp = ESL_MAX(MMX_SP(donor_idx-3,k-1) , DMX_SP(donor_idx-3,k-1));
        /* Check all possible ATACXXX */
        for(s = 0; s < splice_signals->N[p7S_ATAC]; s++) {
          donor_idx = splice_signals->index[p7S_ATAC][s];
 
          if(i - (donor_idx-3) > MIN_INTRON_LENG) {
           /* Create codon with one nuc from before the donor and two from after the acceptor */
             if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-2])) z = sub_dsq[donor_idx-2];
             else                                                          z = p7P_MAXCODONS;

             c2 = p7P_CODON3(z, w, x);
             c2 = p7P_MINIDX(c2, p7P_DEGEN_C);

             PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k),
                                   P_tmp + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c2));
          }
        }
      }

      if ( AGX && i > 4) {
         P_tmp = ESL_MAX(MMX_SP(donor_idx-4,k-1) , DMX_SP(donor_idx-4,k-1)); 
         /* Check all possible XXGTAGX */
         for(s = 0; s < splice_signals->N[p7S_GTAG]; s++) {
           donor_idx = splice_signals->index[p7S_GTAG][s];
           if(i - (donor_idx-4) > MIN_INTRON_LENG) {
             /* Create codon with one nuc from before the donor and two from after the acceptor */
             if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-3])) y = sub_dsq[donor_idx-3];
             else                                                          y = p7P_MAXCODONS;

             if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-2])) z = sub_dsq[donor_idx-2];
             else                                                          z = p7P_MAXCODONS;

             c1 = p7P_CODON3(y, z, x);
             c1 = p7P_MINIDX(c1, p7P_DEGEN_C);  
          
             PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k),
                                   P_tmp + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c1));
           }
         }
        
         /* Check all possible XXGCAGX */
         for(s = 0; s < splice_signals->N[p7S_GCAG]; s++) {
           donor_idx = splice_signals->index[p7S_GCAG][s];
           if(i - (donor_idx-4) > MIN_INTRON_LENG) {
             /* Create codon with one nuc from before the donor and two from after the acceptor */
             if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-3])) y = sub_dsq[donor_idx-3];
             else                                                          y = p7P_MAXCODONS;

             if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-2])) z = sub_dsq[donor_idx-2];
             else                                                          z = p7P_MAXCODONS;

             c1 = p7P_CODON3(y, z, x);
             c1 = p7P_MINIDX(c1, p7P_DEGEN_C); 
             
             PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k),
                                   P_tmp + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c1));
           }
         }
      }
      else if (ACX && i > 4) {
        P_tmp = ESL_MAX(MMX_SP(donor_idx-4,k-1) , DMX_SP(donor_idx-4,k-1));
        /* Check all possible XXATACX */
        for(s = 0; s < splice_signals->N[p7S_ATAC]; s++) {
          donor_idx = splice_signals->index[p7S_ATAC][s];
 
          if(i - (donor_idx-4) > MIN_INTRON_LENG) {
           /* Create codon with one nuc from before the donor and two from after the acceptor */
             if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-2])) z = sub_dsq[donor_idx-2];
             else                                                          z = p7P_MAXCODONS;

             c2 = p7P_CODON3(z, w, x);
             c2 = p7P_MINIDX(c2, p7P_DEGEN_C);

             PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k),
                                   P_tmp + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c1));
          }
        }
      } 

    } 

    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,M-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,M-1),
                  ESL_MAX(DMX_SP(i-3,M-1) + TSC(p7P_DM,M-1),
                          PMX_SP(i-3,M-1) + TSC_P)))         + p7P_MSC_CODON(gm_fs, M, c3);

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,M-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,M-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX(i,p7G_E) = ESL_MAX(MMX_SP(i,M), DMX_SP(i,M));
    XMX(i,p7G_C) = ESL_MAX(XMX(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                           XMX(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);

    XMX(i,p7G_J) = -eslINFINITY;    

  }
  
  gx->M = gm_fs->M;
  gx->L = L;

  //p7_gmx_sp_Dump(stdout, gx, p7_DEFAULT);
  return eslOK;

}


int
p7_splicevitebi_translated_semiglobal_trace(SPLICE_SITE_IDX *splice_signals, float *signal_scores, const ESL_DSQ *sub_dsq, int L, const ESL_GENCODE *gcode, const P7_FS_PROFILE *sub_gm, const P7_GMX *gx, P7_TRACE *tr)
{
  int          i   = L;     /* position in seq (1..L)         */
  int          k   = 0;     /* position in model (1..M)       */
  int          M   = sub_gm->M;
  int          t,u,v,w,x,y,z;
  int          s, donor_idx;
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
  
      if      (XMX(i, p7G_C) + XMX(i-2, p7G_C) || XMX(i, p7G_C) < XMX(i-1, p7G_C))                        scur = p7T_C; 
      if      (esl_FCompare_old(XMX(i, p7G_C), XMX(i-3, p7G_C) + sub_gm->xsc[p7P_C][p7P_LOOP], tol) == eslOK) scur = p7T_C;
      else if (esl_FCompare_old(XMX(i, p7G_C), XMX(i,   p7G_E) + sub_gm->xsc[p7P_E][p7P_MOVE], tol) == eslOK) scur = p7T_E;
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
      c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

      emit = p7P_MSC_CODON(sub_gm, k, c3);

      if      (esl_FCompare_old(MMX_SP(i,k), MMX_SP(i-3, k-1) + TSC(p7P_MM, k-1) + emit, tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(MMX_SP(i,k), IMX_SP(i-3, k-1) + TSC(p7P_IM, k-1) + emit, tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare_old(MMX_SP(i,k), DMX_SP(i-3, k-1) + TSC(p7P_DM, k-1) + emit, tol) == eslOK) scur = p7T_D;      
      else if (esl_FCompare_old(MMX_SP(i,k), PMX_SP(i-3, k-1) + TSC_P            + emit, tol) == eslOK) scur = p7T_P;
      else if (esl_FCompare_old(MMX_SP(i,k), XMX(i-3, p7G_B)  + TSC(p7P_BM, k-1) + emit, tol) == eslOK) scur = p7T_B;      
      else ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
      
      k--; i-=3;
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

    case p7T_P:
      if (PMX_SP(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible P reached at k=%d,i=%d", k,i);

      if(i > 4 && esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i-4])) t = sub_dsq[i-4];
      else                                              t = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i-3])) u = sub_dsq[i-3];
      else                                              u = p7P_MAXCODONS; 
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i-2])) v = sub_dsq[i-2];
      else                                              v = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i-1])) w = sub_dsq[i-1];
      else                                              w = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i]))   x = sub_dsq[i];
      else

      if(t == 0 && u == 2) { //AGXXX 
        for(s = 0; s < splice_signals->N[p7S_GTAG]; s++) {
          donor_idx = splice_signals->index[p7S_GTAG][s];
          if(i - (donor_idx-2) > MIN_INTRON_LENG) {
            c3 = p7P_CODON3(v, w, x);
            c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

            emit = p7P_MSC_CODON(sub_gm, k, c3); 
            
            if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-2, k-1) + signal_scores[p7S_GTAG] + emit, tol) == eslOK) { scur = p7T_M; i = donor_idx-2; c = 0; break; }
            else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-2, k-1) + signal_scores[p7S_GTAG] + emit, tol) == eslOK) { scur = p7T_D; i = donor_idx-2; c = 0; break; }
          }
        }              
        if(scur != p7T_P) break;

        for(s = 0; s < splice_signals->N[p7S_GCAG]; s++) {
          donor_idx = splice_signals->index[p7S_GCAG][s];
          if(i - (donor_idx-2) > MIN_INTRON_LENG) {
            c3 = p7P_CODON3(v, w, x);
            c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

            emit = p7P_MSC_CODON(sub_gm, k, c3);

            if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-2, k-1) + signal_scores[p7S_GCAG] + emit, tol) == eslOK) { scur = p7T_M; i = donor_idx-2; c = 0; break; }
            else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-2, k-1) + signal_scores[p7S_GCAG] + emit, tol) == eslOK) { scur = p7T_D; i = donor_idx-2; c = 0; break; }
          }
        }
        if(scur != p7T_P) break;
      }
      else if(t == 0 && u == 1) { //ACXXX 

        for(s = 0; s < splice_signals->N[p7S_ATAC]; s++) {
          donor_idx = splice_signals->index[p7S_ATAC][s];
          if(i - (donor_idx-2) > MIN_INTRON_LENG) {
            c3 = p7P_CODON3(v, w, x);
            c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

            emit = p7P_MSC_CODON(sub_gm, k, c3);

            if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-2, k-1) + signal_scores[p7S_ATAC] + emit, tol) == eslOK) { scur = p7T_M; i = donor_idx-2; c = 0; break; }
            else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-2, k-1) + signal_scores[p7S_ATAC] + emit, tol) == eslOK) { scur = p7T_D; i = donor_idx-2; c = 0; break; }
          }
        }
        if(scur != p7T_P) break;
      }
   
      if(u == 0 && v == 2) { //XAGXX 
        for(s = 0; s < splice_signals->N[p7S_GTAG]; s++) {
          donor_idx = splice_signals->index[p7S_GTAG][s];
          if(i - (donor_idx-3) > MIN_INTRON_LENG) {
            if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-3])) z = sub_dsq[donor_idx-3];
            else                                                          z = p7P_MAXCODONS;     
            c3 = p7P_CODON3(z, w, x);
            c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

            emit = p7P_MSC_CODON(sub_gm, k, c3);

            if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-3, k-1) + signal_scores[p7S_GTAG] + emit, tol) == eslOK) { scur = p7T_M; i = donor_idx-3; c = 1; break; }
            else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-3, k-1) + signal_scores[p7S_GTAG] + emit, tol) == eslOK) { scur = p7T_D; i = donor_idx-3; c = 1; break; }
          }
        } 
        if(scur != p7T_P) break;

        for(s = 0; s < splice_signals->N[p7S_GCAG]; s++) {
          donor_idx = splice_signals->index[p7S_GCAG][s];
          if(i - (donor_idx-3) > MIN_INTRON_LENG) {
            if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-3])) z = sub_dsq[donor_idx-3];
            else                                                          z = p7P_MAXCODONS;     
            c3 = p7P_CODON3(z, w, x);
            c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

            emit = p7P_MSC_CODON(sub_gm, k, c3);

            if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-3, k-1) + signal_scores[p7S_GCAG] + emit, tol) == eslOK) { scur = p7T_M; i = donor_idx-3; c = 1; break; }
            else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-3, k-1) + signal_scores[p7S_GCAG] + emit, tol) == eslOK) { scur = p7T_D; i = donor_idx-3; c = 1; break; }
          }
        }
        if(scur != p7T_P) break;
      }
      else if(u == 0 && v == 1) { //XACXX 
         for(s = 0; s < splice_signals->N[p7S_ATAC]; s++) {
          donor_idx = splice_signals->index[p7S_ATAC][s];
          if(i - (donor_idx-3) > MIN_INTRON_LENG) {
            if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-3])) z = sub_dsq[donor_idx-3];
            else                                                          z = p7P_MAXCODONS;
            c3 = p7P_CODON3(z, w, x);
            c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

            emit = p7P_MSC_CODON(sub_gm, k, c3);

            if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-3, k-1) + signal_scores[p7S_ATAC] + emit, tol) == eslOK) { scur = p7T_M; i = donor_idx-3; c = 1; break; }
            else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-3, k-1) + signal_scores[p7S_ATAC] + emit, tol) == eslOK) { scur = p7T_D; i = donor_idx-3; c = 1; break; }
          }
        }
        if(scur != p7T_P) break;
      }

      if(v == 0 && w == 2 && i > 4) { //XXAGX
        for(s = 0; s < splice_signals->N[p7S_GTAG]; s++) {
          donor_idx = splice_signals->index[p7S_GTAG][s];
          if(i - (donor_idx-4) > MIN_INTRON_LENG) {
            if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-4])) y = sub_dsq[donor_idx-4];
            else                                                          y = p7P_MAXCODONS;
            if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-3])) z = sub_dsq[donor_idx-3];
            else                                                          z = p7P_MAXCODONS;
            c3 = p7P_CODON3(y, z, x);
            c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

            emit = p7P_MSC_CODON(sub_gm, k, c3);

            if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-4, k-1) + signal_scores[p7S_GTAG] + emit, tol) == eslOK) { scur = p7T_M; i = donor_idx-4; c = 2; break; }
            else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-4, k-1) + signal_scores[p7S_GTAG] + emit, tol) == eslOK) { scur = p7T_D; i = donor_idx-4; c = 2; break; }
          }
        }  
        if(scur != p7T_P) break;
 
        for(s = 0; s < splice_signals->N[p7S_GCAG]; s++) {
          donor_idx = splice_signals->index[p7S_GCAG][s];
          if(i - (donor_idx-4) > MIN_INTRON_LENG) {
            if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-4])) y = sub_dsq[donor_idx-4];
            else                                                          y = p7P_MAXCODONS;
            if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-3])) z = sub_dsq[donor_idx-3];
            else                                                          z = p7P_MAXCODONS;
            c3 = p7P_CODON3(y, z, x);
            c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

            emit = p7P_MSC_CODON(sub_gm, k, c3);

            if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-4, k-1) + signal_scores[p7S_GCAG] + emit, tol) == eslOK) { scur = p7T_M; i = donor_idx-4; c = 2; break; }
            else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-4, k-1) + signal_scores[p7S_GCAG] + emit, tol) == eslOK) { scur = p7T_D; i = donor_idx-4; c = 2; break; }
          }
        }
        if(scur != p7T_P) break; 
      }
      else if(v == 0 && w == 1) { //XXACX
        for(s = 0; s < splice_signals->N[p7S_ATAC]; s++) {
          donor_idx = splice_signals->index[p7S_ATAC][s];
          if(i - (donor_idx-4) > MIN_INTRON_LENG) {
            if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-4])) y = sub_dsq[donor_idx-4];
            else                                                          y = p7P_MAXCODONS;
            if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-3])) z = sub_dsq[donor_idx-3];
            else                                                          z = p7P_MAXCODONS;
            c3 = p7P_CODON3(y, z, x);
            c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

            emit = p7P_MSC_CODON(sub_gm, k, c3);

            if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-4, k-1) + signal_scores[p7S_ATAC] + emit, tol) == eslOK) { scur = p7T_M; i = donor_idx-4; c = 2; break; }
            else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-4, k-1) + signal_scores[p7S_ATAC] + emit, tol) == eslOK) { scur = p7T_D; i = donor_idx-4; c = 2; break; }
          }
        }
        if(scur != p7T_P) break;
      }
      if(scur == p7T_P)  ESL_EXCEPTION(eslFAIL, "P at k=%d,i=%d couldn't be traced", k,i);  

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

    if     (scur == p7T_M) c = 3;
    else if(scur != p7T_P) c = 0;
    
    if ((status = p7_trace_fs_Append(tr, scur, k, i, c)) != eslOK) return status;
  
    /* For NCJ, we had to defer i decrement. */
    if ( (scur == p7T_N || scur == p7T_C) && scur == sprv) i--; 
   
    sprv = scur;
  } /* end traceback, at S state */

  tr->M = sub_gm->M;
  tr->L = L;

  return p7_trace_fs_Reverse(tr);
}





