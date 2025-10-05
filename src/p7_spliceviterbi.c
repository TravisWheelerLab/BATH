/* Spliced Viterbi algorithm and trace back - translated and fraemshift versions. */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_splice.h"

#define TSC_P -14

int
p7_spliceviterbi_translated_semiglobal(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx)
{
  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float      **score = pli->sig_idx->score;
  float       *signal_scores = pli->signal_scores;
  int        **index = pli->sig_idx->index;
  int        **lookback = pli->sig_idx->lookback;
  int          M    = gm_fs->M;
  int          i,k,s;
  int          c3, c2, c1;
  int          u, v, w, x;
  int          AGXXX, ACXXX;
  int          AGXX, ACXX;
  int          AGX,  ACX;
  int          GT, GC, AT;
  int          nuc1, nuc2;
  float        TMP_SC; 
 
  /*Initialize the signal index */
  for(i = 0; i <= L; i++) {
    for(k = 0; k < M; k++) {
      lookback[i][k] = -1;
    }
  }

  for(k = 0; k <  M; k++) {
    for(s = 0; s < SIGNAL_MEM_SIZE; s++) {
      index[k][s] = -1;
      score[k][s] = -eslINFINITY;
    }
  }

  /* Initialization of the zero row.  */
  XMX(0,p7G_N) = 0.;                                     /* S->N, p=1            */
  XMX(0,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE];            /* S->N->B, no N-tail   */ 
  XMX(0,p7G_E) = XMX(0,p7G_J) = XMX(0,p7G_C) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_SP(0,k) = IMX_SP(0,k) = PMX_SP(0,k) = -eslINFINITY;

  DMX_SP(0,0) = -eslINFINITY;
  DMX_SP(0,1) = XMX(0,p7G_B) + TSC(p7P_BM,0); 
  for (k = 2; k <= M; k++)
    DMX_SP(0,k) = DMX_SP(0,k-1) + TSC(p7P_DD,k-1);

  /*Special cases for the first 2 rows */
  u = v = w = x = -1;
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

  /*Special cases for the rows 3-6*/
  for (i = 3; i <= 6; i++) 
  {
    u = v;
    v = w;
    w = x;

    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                                x = p7P_MAXCODONS;

    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    XMX(i,p7G_N) = XMX(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];
    XMX(i,p7G_B) = XMX(i,p7G_N)   + gm_fs->xsc[p7P_N][p7P_MOVE];


    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    MMX_SP(i,1) = XMX(i-3,p7G_B) + TSC(p7P_BM,0) + p7P_MSC_CODON(gm_fs, 1, c3); //B->M

    if(p7P_MSC_CODON(gm_fs, 1, c3) == -eslINFINITY)  IMX_SP(i,1) = -eslINFINITY;
    else                                             IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,1),
                                                     IMX_SP(i-3,1) + TSC(p7P_II,1));

    DMX_SP(i,1) = XMX(i,p7G_B) + TSC(p7P_BM,0);

    PMX_SP(i,1) = -eslINFINITY;

    for (k = 2; k < M; k++) {
      
      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,k-1),
                            DMX_SP(i-3,k-1) + TSC(p7P_DM,k-1))) + p7P_MSC_CODON(gm_fs, k, c3);
      
      if(p7P_MSC_CODON(gm_fs, k, c3) == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;
      else                                            IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,k),
                                                                            IMX_SP(i-3,k) + TSC(p7P_II,k));

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,k-1));

      PMX_SP(i,k) = -eslINFINITY;
    } 

    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,M-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,M-1),
                          DMX_SP(i-3,M-1) + TSC(p7P_DM,M-1)))+ p7P_MSC_CODON(gm_fs, M, c3);

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,M-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,M-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX(i,p7G_E) = ESL_MAX(MMX_SP(i,M), DMX_SP(i,M));
    XMX(i,p7G_C) = ESL_MAX(XMX(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                           XMX(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);

    XMX(i,p7G_J) = -eslINFINITY;

  }


  AGXXX = ACXXX = AGXX = ACXX = AGX = ACX  = FALSE;  
  GT = GC = AT = FALSE; 
  /* Main DP recursion */
  for (i = 7; i <= L; i++) 
  {
    /* get nucleotides and codon */
    u = v;
    v = w;
    w = x;

    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                            x = p7P_MAXCODONS;
    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);
 
    /* get acceptor site */
    AGXXX = AGXX;
    ACXXX = ACXX;
    AGXX  = AGX;
    ACXX  = ACX;
    if(SIGNAL(v, w) == ACCEPT_AG) AGX = TRUE;
    else                          AGX = FALSE;
    
    if(SIGNAL(v, w) == ACCEPT_AC) ACX = TRUE;
    else                          ACX = FALSE;

    /* get donor site */
    if(SIGNAL(w, x) == DONOR_GT) GT = TRUE;
    else                         GT = FALSE;
   
    if(SIGNAL(w, x) == DONOR_GC) GC = TRUE;
    else                         GC = FALSE;
 
    if(SIGNAL(w, x) == DONOR_AT) AT = TRUE;
    else                         AT = FALSE;
 
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



      /* For each valid splice signal pair, codon type, and emissions 
       * (total of 63, but max 34 per sequence position) find the 
       * maximum score of the P state */
      PMX_SP(i,k) = -eslINFINITY;
      if(AGXXX) {
      
        TMP_SC = SSX0(k, p7S_GTAG) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c3);
        if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX0(k, p7S_GTAG);
        }

        TMP_SC = SSX0(k, p7S_GCAG) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c3);
        if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX0(k, p7S_GCAG);
        }
      }
      else if(ACXXX) {
        TMP_SC = SSX0(k, p7S_ATAC) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c3);
        if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX0(k, p7S_ATAC);
        }
      }

      if(AGXX) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(nuc1, w, x);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);
       
          TMP_SC = SSX1(k, p7S_GTAG, nuc1) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c2); 
          if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX1(k, p7S_GTAG, nuc1);
          }

          TMP_SC = SSX1(k, p7S_GCAG, nuc1) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c2);
          if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX1(k, p7S_GCAG, nuc1);
          }
        }
      }
      else if(ACXX) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(nuc1, w, x);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);

          TMP_SC = SSX1(k, p7S_ATAC, nuc1) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c2);
          if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX1(k, p7S_ATAC, nuc1);
          }
        }
      }
    
      if(AGX) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(nuc1, nuc2, x);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);

            TMP_SC = SSX2(k, p7S_GTAG, nuc1, nuc2) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c1);
            if(TMP_SC > PMX_SP(i,k)) {
              PMX_SP(i,k) = TMP_SC;
              lookback[i][k] = SIX2(k, p7S_GTAG, nuc1, nuc2);
            }

            TMP_SC = SSX2(k, p7S_GCAG, nuc1, nuc2) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c1);
            if(TMP_SC > PMX_SP(i,k)) {
              PMX_SP(i,k) = TMP_SC;
              lookback[i][k] = SIX2(k, p7S_GCAG, nuc1, nuc2);
            }
          }
        }
      }
      else if(ACX) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(nuc1, nuc2, x);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);

            TMP_SC = SSX2(k, p7S_ATAC, nuc1, nuc2) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c1);
            if(TMP_SC > PMX_SP(i,k)) {
              PMX_SP(i,k) = TMP_SC;
              lookback[i][k] = SIX2(k, p7S_ATAC, nuc1, nuc2);
            }
          }
        }
      }
      
      
      /* Record best scoring donor sites*/
      if(GT) {
        
        TMP_SC = ESL_MAX(MMX_SP(i-4,k-1), DMX_SP(i-4,k-1)); 
        if(u < 4 && v < 4 && TMP_SC > SSX2(k, p7S_GTAG, u, v)) {
          SSX2(k, p7S_GTAG, u, v) = TMP_SC;
          SIX2(k, p7S_GTAG, u, v) = i;
        }
 
        TMP_SC = ESL_MAX(MMX_SP(i-3,k-1), DMX_SP(i-3,k-1));
        if(v < 4 && TMP_SC > SSX1(k, p7S_GTAG, v)) {
          SSX1(k, p7S_GTAG, v) = TMP_SC;
          SIX1(k, p7S_GTAG, v) = i;
        } 

        TMP_SC = ESL_MAX(MMX_SP(i-2,k-1), DMX_SP(i-2,k-1));
        if(TMP_SC > SSX0(k, p7S_GTAG)) {
          SSX0(k, p7S_GTAG) = TMP_SC;
          SIX0(k, p7S_GTAG) = i;
        }
         
      }
      else if(GC) {
        
        TMP_SC = ESL_MAX(MMX_SP(i-4,k-1), DMX_SP(i-4,k-1));
        if(u < 4 && v < 4 && TMP_SC > SSX2(k, p7S_GCAG, u, v)) {
          SSX2(k, p7S_GCAG, u, v) = TMP_SC;
          SIX2(k, p7S_GCAG, u, v) = i;
        }

        TMP_SC = ESL_MAX(MMX_SP(i-3,k-1), DMX_SP(i-3,k-1));
        if(v < 4 && TMP_SC > SSX1(k, p7S_GCAG, v)) {
          SSX1(k, p7S_GCAG, v) = TMP_SC;
          SIX1(k, p7S_GCAG, v) = i;
        }

        TMP_SC = ESL_MAX(MMX_SP(i-2,k-1), DMX_SP(i-2,k-1));
        if(TMP_SC > SSX0(k, p7S_GCAG)) {
          SSX0(k, p7S_GCAG) = TMP_SC;
          SIX0(k, p7S_GCAG) = i;
        }
        

      }
      else if(AT) {

        TMP_SC = ESL_MAX(MMX_SP(i-4,k-1), DMX_SP(i-4,k-1));
        if(u < 4 && v < 4 && TMP_SC > SSX2(k, p7S_ATAC, u, v)) {
          SSX2(k, p7S_ATAC, u, v) = TMP_SC;
          SIX2(k, p7S_ATAC, u, v) = i;
        }

        TMP_SC = ESL_MAX(MMX_SP(i-3,k-1), DMX_SP(i-3,k-1));
        if(v < 4 && TMP_SC > SSX1(k, p7S_ATAC, v)) {
          SSX1(k, p7S_ATAC, v) = TMP_SC;
          SIX1(k, p7S_ATAC, v) = i;
        }

        TMP_SC = ESL_MAX(MMX_SP(i-2,k-1), DMX_SP(i-2,k-1));
        if(TMP_SC > SSX0(k, p7S_ATAC)) {
          SSX0(k, p7S_ATAC) = TMP_SC;
          SIX0(k, p7S_ATAC) = i;
        }
      }
    } //end loop over M 

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

  } // end loop over L
  
  gx->M = gm_fs->M;
  gx->L = L;

  //p7_gmx_sp_Dump(stdout, gx, p7_DEFAULT);
  return eslOK;

}


int
p7_splicevitebi_translated_semiglobal_trace(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, int L, const ESL_GENCODE *gcode, const P7_FS_PROFILE *sub_gm, const P7_GMX *gx, P7_TRACE *tr)
{
  int          i   = L;     /* position in seq (1..L)         */
  int          k   = 0;     /* position in model (1..M)       */
  int          M   = sub_gm->M;
  int          t,u,v,w,x;
  int          donor_idx;
  int          c, c0, c1, c2, c3;
  int          sprv,scur,snxt;
  float        emit, emit0, emit1, emit2;
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

  snxt = p7T_X; 
  sprv = p7T_C;
  
 
  while (sprv != p7T_S) {
    switch (sprv) {
    case p7T_C:     /* C(i) comes from C(i-1) or E(i) */
      if   (XMX(i,p7G_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
  
      if      (XMX(i, p7G_C) < XMX(i-2, p7G_C) || XMX(i, p7G_C) < XMX(i-1, p7G_C))                        scur = p7T_C; 
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

      scur = snxt;
      k--; i = donor_idx - c - 2;
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


    /* Need to find which type of splice codon we have */
    if (scur == p7T_P) {
      if (PMX_SP(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible P reached at k=%d,i=%d", k,i);

      donor_idx = pli->sig_idx->lookback[i][k];
      if(donor_idx == -1) ESL_EXCEPTION(eslFAIL, "P at k=%d,i=%d not in lookback", k,i);     

      /* nucleotides upstream of donor */
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-3])) t = sub_dsq[donor_idx-3];
      else                                                          t = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[donor_idx-2])) u = sub_dsq[donor_idx-2];
      else                                                          u = p7P_MAXCODONS;

      /* nucleotides upstream of acceptor */
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i-2])) v = sub_dsq[i-2];
      else                                                  v = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i-1])) w = sub_dsq[i-1];
      else                                                  w = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i]))   x = sub_dsq[i];
      else                                                  x = p7P_MAXCODONS; 
     
      c0 = p7P_CODON3(v, w, x);
      c0 = p7P_MINIDX(c0, p7P_DEGEN_C);
      emit0 = p7P_MSC_CODON(sub_gm, k, c0);

      c1 = p7P_CODON3(u, w, x);
      c1 = p7P_MINIDX(c1, p7P_DEGEN_C);
      emit1 = p7P_MSC_CODON(sub_gm, k, c1);

      c2 = p7P_CODON3(t, u, x);
      c2 = p7P_MINIDX(c2, p7P_DEGEN_C);
      emit2 = p7P_MSC_CODON(sub_gm, k, c2);
     
      if(SIGNAL(sub_dsq[donor_idx-1], sub_dsq[donor_idx]) == DONOR_GT) {
        
        if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-2,k-1) + pli->signal_scores[p7S_GTAG] + emit0, tol) == eslOK) { snxt = p7T_M; c = 0; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-2,k-1) + pli->signal_scores[p7S_GTAG] + emit0, tol) == eslOK) { snxt = p7T_D; c = 0; }
        else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-3,k-1) + pli->signal_scores[p7S_GTAG] + emit1, tol) == eslOK) { snxt = p7T_M; c = 1; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-3,k-1) + pli->signal_scores[p7S_GTAG] + emit1, tol) == eslOK) { snxt = p7T_D; c = 1; }
        else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-4,k-1) + pli->signal_scores[p7S_GTAG] + emit2, tol) == eslOK) { snxt = p7T_M; c = 2; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-4,k-1) + pli->signal_scores[p7S_GTAG] + emit2, tol) == eslOK) { snxt = p7T_D; c = 2; }
        else ESL_EXCEPTION(eslFAIL, "P at k=%d,i=%d couldn't be traced", k,i);
      }
      else if (SIGNAL(sub_dsq[donor_idx-1], sub_dsq[donor_idx]) == DONOR_GC) {
        if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-2,k-1) + pli->signal_scores[p7S_GCAG] + emit0, tol) == eslOK) { snxt = p7T_M; c = 0; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-2,k-1) + pli->signal_scores[p7S_GCAG] + emit0, tol) == eslOK) { snxt = p7T_D; c = 0; }
        else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-3,k-1) + pli->signal_scores[p7S_GCAG] + emit1, tol) == eslOK) { snxt = p7T_M; c = 1; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-3,k-1) + pli->signal_scores[p7S_GCAG] + emit1, tol) == eslOK) { snxt = p7T_D; c = 1; }
        else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-4,k-1) + pli->signal_scores[p7S_GCAG] + emit2, tol) == eslOK) { snxt = p7T_M; c = 2; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-4,k-1) + pli->signal_scores[p7S_GCAG] + emit2, tol) == eslOK) { snxt = p7T_D; c = 2; }
        else ESL_EXCEPTION(eslFAIL, "P at k=%d,i=%d couldn't be traced", k,i); 
      }
      else if (SIGNAL(sub_dsq[donor_idx-1], sub_dsq[donor_idx]) == DONOR_AT) {
        if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-2,k-1) + pli->signal_scores[p7S_ATAC] + emit0, tol) == eslOK) { snxt = p7T_M; c = 0; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-2,k-1) + pli->signal_scores[p7S_ATAC] + emit0, tol) == eslOK) { snxt = p7T_D; c = 0; }
        else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-3,k-1) + pli->signal_scores[p7S_ATAC] + emit1, tol) == eslOK) { snxt = p7T_M; c = 1; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-3,k-1) + pli->signal_scores[p7S_ATAC] + emit1, tol) == eslOK) { snxt = p7T_D; c = 1; }
        else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-4,k-1) + pli->signal_scores[p7S_ATAC] + emit2, tol) == eslOK) { snxt = p7T_M; c = 2; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-4,k-1) + pli->signal_scores[p7S_ATAC] + emit2, tol) == eslOK) { snxt = p7T_D; c = 2; }
        else ESL_EXCEPTION(eslFAIL, "P at k=%d,i=%d couldn't be traced", k,i); 
      }
      else ESL_EXCEPTION(eslFAIL, "P at k=%d,i=%d not a donor site", k,donor_idx); 
    } 

    if      (scur == p7T_M) c = 3;
    else if (scur != p7T_P) c = 0;

    if ((status = p7_trace_fs_Append(tr, scur, k, i, c)) != eslOK) return status;
  
    /* For NCJ, we had to defer i decrement. */
    if ( (scur == p7T_N || scur == p7T_C) && scur == sprv) i--; 
   
    sprv = scur;
  } /* end traceback, at S state */

  tr->M = sub_gm->M;
  tr->L = L;

  return p7_trace_fs_Reverse(tr);
}



