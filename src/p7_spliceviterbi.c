/* Spliced Viterbi algorithm and trace back - translated and fraemshift versions. */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_splice.h"

#define TSC_P -10


// Global on model from position 1 to M
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
  int          t, u, v, w, x;
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
  XMX_SP(0,p7G_N) = 0.;                                     /* S->N, p=1            */
  XMX_SP(0,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE];            /* S->N->B, no N-tail   */ 
  XMX_SP(0,p7G_E) = XMX_SP(0,p7G_C) = XMX_SP(0,p7G_J) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_SP(0,k) = DMX_SP(0,k) = IMX_SP(0,k) = PMX_SP(0,k) = -eslINFINITY;

  /*Special cases for the first 2 rows */
  t = u = v = w = x = -1;
  for(i = 1; i <= 2; i++)
  {
    w = x;
  
    /* if new nucleotide is not A,C,G, or T set it to placeholder vlaue */
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                                x = p7P_MAXCODONS;

    XMX_SP(i,p7G_N) =  0; 
    XMX_SP(i,p7G_B) =  gm_fs->xsc[p7P_N][p7P_MOVE]; 
    XMX_SP(i,p7G_J) = -eslINFINITY;
 
    for (k = 0; k <= M; k++) 
      MMX_SP(i,k) = DMX_SP(i,k) = IMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = -eslINFINITY;

  } 

  /*Special cases for the rows 3-MIN_INTRON_LENG+4*/
  AGXXX = ACXXX = AGXX = ACXX = AGX = ACX  = FALSE;  
  for (i = 3; i <= MIN_INTRON_LENG+4; i++) 
  {
    v = w;
    w = x;

    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                                x = p7P_MAXCODONS;

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

    XMX_SP(i,p7G_N) = XMX_SP(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];
    XMX_SP(i,p7G_B) = XMX_SP(i,p7G_N)   + gm_fs->xsc[p7P_N][p7P_MOVE];
    XMX_SP(i,p7G_J) = -eslINFINITY;

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    MMX_SP(i,1) = XMX_SP(i-3,p7G_B) + p7P_MSC_CODON(gm_fs, 1, c3); //B->M

    if(p7P_MSC_CODON(gm_fs, 1, c3) == -eslINFINITY)  IMX_SP(i,1) = -eslINFINITY;
    else                                             IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,1),
                                                     IMX_SP(i-3,1) + TSC(p7P_II,1));
    DMX_SP(i,1) = PMX_SP(i,1) = -eslINFINITY;

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

    XMX_SP(i,p7G_E) = ESL_MAX(MMX_SP(i,M), DMX_SP(i,M));
    XMX_SP(i,p7G_C) = ESL_MAX(XMX_SP(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                           XMX_SP(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);

  }

  GT = GC = AT = FALSE; 
  /* Main DP recursion */
  for (i = MIN_INTRON_LENG+5; i <= L; i++) 
  {
    /* get nucleotides and codon */
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
    t = sub_dsq[i-MIN_INTRON_LENG-3];
    u = sub_dsq[i-MIN_INTRON_LENG-2]; 
    if(SIGNAL(sub_dsq[i-MIN_INTRON_LENG-1], sub_dsq[i-MIN_INTRON_LENG]) == DONOR_GT) GT = TRUE;
    else                         GT = FALSE;
 
    if(SIGNAL(sub_dsq[i-MIN_INTRON_LENG-1], sub_dsq[i-MIN_INTRON_LENG]) == DONOR_GC) GC = TRUE;
    else                         GC = FALSE;
 
    if(SIGNAL(sub_dsq[i-MIN_INTRON_LENG-1], sub_dsq[i-MIN_INTRON_LENG]) == DONOR_AT) AT = TRUE;
    else                         AT = FALSE;

    XMX_SP(i,p7G_N) = XMX_SP(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];
    XMX_SP(i,p7G_B) = XMX_SP(i,p7G_N)   + gm_fs->xsc[p7P_N][p7P_MOVE];
    XMX_SP(i,p7G_J) = -eslINFINITY;   
 
    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    MMX_SP(i,1) = XMX_SP(i-3,p7G_B) + p7P_MSC_CODON(gm_fs, 1, c3); //B->M

    if(p7P_MSC_CODON(gm_fs, 1, c3) == -eslINFINITY) IMX_SP(i,1) = -eslINFINITY;
    else                                            IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,1),
                                                    IMX_SP(i-3,1) + TSC(p7P_II,1));
    
    DMX_SP(i,1) = PMX_SP(i,1) = -eslINFINITY;

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
    } 

    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,M-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,M-1),
                  ESL_MAX(DMX_SP(i-3,M-1) + TSC(p7P_DM,M-1),
                          PMX_SP(i-3,M-1) + TSC_P)))         + p7P_MSC_CODON(gm_fs, M, c3);

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,M-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,M-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = ESL_MAX(MMX_SP(i,M), DMX_SP(i,M));
    XMX_SP(i,p7G_C) = ESL_MAX(XMX_SP(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                              XMX_SP(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);
    
    if(AGXXX) {
      for (k = 2; k < M; k++) {      
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
    }
    else if(ACXXX) {
      for (k = 2; k < M; k++) {
        TMP_SC = SSX0(k, p7S_ATAC) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c3);
        if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX0(k, p7S_ATAC);
        }
      }
    }

    if(AGXX) {
      for (k = 2; k < M; k++) {
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
    }
    else if(ACXX) {
      for (k = 2; k < M; k++) {
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
    }
    if(AGX) {
      for (k = 2; k < M; k++) {
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
    }
    else if(ACX) {
      for (k = 2; k < M; k++) {
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
    }
      
    /* Record best scoring donor sites*/
    if(GT) {
      for (k = 2; k < M; k++) {  
        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-4,k-1), DMX_SP(i-MIN_INTRON_LENG-4,k-1)); 
        if(t < 4 && u < 4 && TMP_SC > SSX2(k, p7S_GTAG, t, u)) {
          SSX2(k, p7S_GTAG, t, u) = TMP_SC;
          SIX2(k, p7S_GTAG, t, u) = i-MIN_INTRON_LENG;
        }
 
        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-3,k-1), DMX_SP(i-MIN_INTRON_LENG-3,k-1));
        if(u < 4 && TMP_SC > SSX1(k, p7S_GTAG, u)) {
          SSX1(k, p7S_GTAG, u) = TMP_SC;
          SIX1(k, p7S_GTAG, u) = i-MIN_INTRON_LENG;
        } 

        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-2,k-1), DMX_SP(i-MIN_INTRON_LENG-2,k-1));
        if(TMP_SC > SSX0(k, p7S_GTAG)) {
          SSX0(k, p7S_GTAG) = TMP_SC;
          SIX0(k, p7S_GTAG) = i-MIN_INTRON_LENG;
        }
      }
    }
    else if(GC) {
      for (k = 2; k < M; k++) {    
        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-4,k-1), DMX_SP(i-MIN_INTRON_LENG-4,k-1));
        if(t < 4 && u < 4 && TMP_SC > SSX2(k, p7S_GCAG, t, u)) {
          SSX2(k, p7S_GCAG, t, u) = TMP_SC;
          SIX2(k, p7S_GCAG, t, u) = i-MIN_INTRON_LENG;
        }

        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-3,k-1), DMX_SP(i-MIN_INTRON_LENG-3,k-1));
        if(u < 4 && TMP_SC > SSX1(k, p7S_GCAG, u)) {
          SSX1(k, p7S_GCAG, u) = TMP_SC;
          SIX1(k, p7S_GCAG, u) = i-MIN_INTRON_LENG;
        }

        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-2,k-1), DMX_SP(i-MIN_INTRON_LENG-2,k-1));
        if(TMP_SC > SSX0(k, p7S_GCAG)) {
          SSX0(k, p7S_GCAG) = TMP_SC;
          SIX0(k, p7S_GCAG) = i-MIN_INTRON_LENG;
        }
      }
    }
    else if(AT) {
      for (k = 2; k < M; k++) {
        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-4,k-1), DMX_SP(i-MIN_INTRON_LENG-4,k-1));
        if(t < 4 && u < 4 && TMP_SC > SSX2(k, p7S_ATAC, t, u)) {
          SSX2(k, p7S_ATAC, t, u) = TMP_SC;
          SIX2(k, p7S_ATAC, t, u) = i-MIN_INTRON_LENG;
        }

        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-3,k-1), DMX_SP(i-MIN_INTRON_LENG-3,k-1));
        if(u < 4 && TMP_SC > SSX1(k, p7S_ATAC, u)) {
          SSX1(k, p7S_ATAC, u) = TMP_SC;
          SIX1(k, p7S_ATAC, u) = i-MIN_INTRON_LENG;
        }

        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-2,k-1), DMX_SP(i-MIN_INTRON_LENG-2,k-1));
        if(TMP_SC > SSX0(k, p7S_ATAC)) {
          SSX0(k, p7S_ATAC) = TMP_SC;
          SIX0(k, p7S_ATAC) = i-MIN_INTRON_LENG;
        }
      }
    }
    
  } // end loop over L
  
  gx->M = gm_fs->M;
  gx->L = L;
 // printf("gx->L %d\n", gx->L );
 // p7_gmx_sp_Dump(stdout, gx, p7_DEFAULT);
  return eslOK;

}

int
p7_spliceviterbi_translated_semiglobal2(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx, int* remove_idx, int remove_cnt)
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
  int          r_i;
  int          c3, c2, c1;
  int          t, u, v, w, x;
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
  XMX_SP(0,p7G_N) = 0.;                                     /* S->N, p=1            */
  XMX_SP(0,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE];            /* S->N->B, no N-tail   */ 
  XMX_SP(0,p7G_E) = XMX_SP(0,p7G_C) = XMX_SP(0,p7G_J) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_SP(0,k) = DMX_SP(0,k) = IMX_SP(0,k) = PMX_SP(0,k) = -eslINFINITY;

  /*Special cases for the first 2 rows */
  t = u = v = w = x = -1;
  for(i = 1; i <= 2; i++)
  {
    w = x;
  
    /* if new nucleotide is not A,C,G, or T set it to placeholder vlaue */
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                                x = p7P_MAXCODONS;

    XMX_SP(i,p7G_N) =  0; 
    XMX_SP(i,p7G_B) =  gm_fs->xsc[p7P_N][p7P_MOVE]; 
    XMX_SP(i,p7G_J) = -eslINFINITY;
 
    for (k = 0; k <= M; k++) 
      MMX_SP(i,k) = DMX_SP(i,k) = IMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = -eslINFINITY;

  } 
  r_i = 0;
  /*Special cases for the rows 3-MIN_INTRON_LENG+4*/
  AGXXX = ACXXX = AGXX = ACXX = AGX = ACX  = FALSE;  
  for (i = 3; i <= MIN_INTRON_LENG+4; i++) 
  {
    v = w;
    w = x;

    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                                x = p7P_MAXCODONS;

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

    XMX_SP(i,p7G_N) = XMX_SP(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];
    XMX_SP(i,p7G_B) = XMX_SP(i,p7G_N)   + gm_fs->xsc[p7P_N][p7P_MOVE];
    XMX_SP(i,p7G_J) = -eslINFINITY;

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;
    
    MMX_SP(i,1) = XMX_SP(i-3,p7G_B) + p7P_MSC_CODON(gm_fs, 1, c3); //B->M

    if(p7P_MSC_CODON(gm_fs, 1, c3) == -eslINFINITY)  IMX_SP(i,1) = -eslINFINITY;
    else                                             IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,1),
                                                     IMX_SP(i-3,1) + TSC(p7P_II,1));
    DMX_SP(i,1) = PMX_SP(i,1) = -eslINFINITY;

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
  
      MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,M-1),
                    ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,M-1),
                            DMX_SP(i-3,M-1) + TSC(p7P_DM,M-1)))+ p7P_MSC_CODON(gm_fs, M, c3);
  
      IMX_SP(i,M) = -eslINFINITY;
  
      DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,M-1),
                            DMX_SP(i,M-1) + TSC(p7P_DD,M-1));
  
      PMX_SP(i,M) = -eslINFINITY;
    }
 
    XMX_SP(i,p7G_E) = ESL_MAX(MMX_SP(i,M), DMX_SP(i,M));
    XMX_SP(i,p7G_C) = ESL_MAX(XMX_SP(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                              XMX_SP(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);

    /* create 3 frames where only special states or the P state can emit */
    if(r_i < remove_cnt         &&  
       (i == remove_idx[r_i] - 1 ||
        i == remove_idx[r_i]     ||
        i == remove_idx[r_i] + 1)) {
        
      for(k = 1; k <= M; k++)
        MMX_SP(i,k) = IMX_SP(i,k) = DMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;
             
      if(i == remove_idx[r_i] + 1) r_i++;
        
    }

  }

  GT = GC = AT = FALSE; 
  /* Main DP recursion */
  for (i = MIN_INTRON_LENG+5; i <= L; i++) 
  {
    /* get nucleotides and codon */
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
    t = sub_dsq[i-MIN_INTRON_LENG-3];
    u = sub_dsq[i-MIN_INTRON_LENG-2]; 
    if(SIGNAL(sub_dsq[i-MIN_INTRON_LENG-1], sub_dsq[i-MIN_INTRON_LENG]) == DONOR_GT) GT = TRUE;
    else                         GT = FALSE;
 
    if(SIGNAL(sub_dsq[i-MIN_INTRON_LENG-1], sub_dsq[i-MIN_INTRON_LENG]) == DONOR_GC) GC = TRUE;
    else                         GC = FALSE;
 
    if(SIGNAL(sub_dsq[i-MIN_INTRON_LENG-1], sub_dsq[i-MIN_INTRON_LENG]) == DONOR_AT) AT = TRUE;
    else                         AT = FALSE;

    XMX_SP(i,p7G_N) = XMX_SP(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];
    XMX_SP(i,p7G_B) = XMX_SP(i,p7G_N)   + gm_fs->xsc[p7P_N][p7P_MOVE];
    XMX_SP(i,p7G_J) = -eslINFINITY;   
 
    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;


    MMX_SP(i,1) = XMX_SP(i-3,p7G_B) + p7P_MSC_CODON(gm_fs, 1, c3); //B->M

    if(p7P_MSC_CODON(gm_fs, 1, c3) == -eslINFINITY) IMX_SP(i,1) = -eslINFINITY;
    else                                            IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,1),
                                                    IMX_SP(i-3,1) + TSC(p7P_II,1));
    
    DMX_SP(i,1) = PMX_SP(i,1) = -eslINFINITY;

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
    } 

    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,M-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,M-1),
                  ESL_MAX(DMX_SP(i-3,M-1) + TSC(p7P_DM,M-1),
                          PMX_SP(i-3,M-1) + TSC_P)))         + p7P_MSC_CODON(gm_fs, M, c3);

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,M-1),
                        DMX_SP(i,M-1) + TSC(p7P_DD,M-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = ESL_MAX(MMX_SP(i,M), DMX_SP(i,M));
    XMX_SP(i,p7G_C) = ESL_MAX(XMX_SP(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                           XMX_SP(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);
    if(AGXXX) {
      for (k = 2; k < M; k++) {      
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
    }
    else if(ACXXX) {
      for (k = 2; k < M; k++) {
        TMP_SC = SSX0(k, p7S_ATAC) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c3);
        if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX0(k, p7S_ATAC);
        }
      }
    }

    if(AGXX) {
      for (k = 2; k < M; k++) {
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
    }
    else if(ACXX) {
      for (k = 2; k < M; k++) {
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
    }
    if(AGX) {
      for (k = 2; k < M; k++) {
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
    }
    else if(ACX) {
      for (k = 2; k < M; k++) {
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
    }

    /* Record best scoring donor sites*/
    if(GT) {
      for (k = 2; k < M; k++) {  
        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-4,k-1), DMX_SP(i-MIN_INTRON_LENG-4,k-1)); 
        if(t < 4 && u < 4 && TMP_SC > SSX2(k, p7S_GTAG, t, u)) {
          SSX2(k, p7S_GTAG, t, u) = TMP_SC;
          SIX2(k, p7S_GTAG, t, u) = i-MIN_INTRON_LENG;
        }
 
        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-3,k-1), DMX_SP(i-MIN_INTRON_LENG-3,k-1));
        if(u < 4 && TMP_SC > SSX1(k, p7S_GTAG, u)) {
          SSX1(k, p7S_GTAG, u) = TMP_SC;
          SIX1(k, p7S_GTAG, u) = i-MIN_INTRON_LENG;
        } 

        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-2,k-1), DMX_SP(i-MIN_INTRON_LENG-2,k-1));
        if(TMP_SC > SSX0(k, p7S_GTAG)) {
          SSX0(k, p7S_GTAG) = TMP_SC;
          SIX0(k, p7S_GTAG) = i-MIN_INTRON_LENG;
        }
      }
    }
    else if(GC) {
      for (k = 2; k < M; k++) {    
        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-4,k-1), DMX_SP(i-MIN_INTRON_LENG-4,k-1));
        if(t < 4 && u < 4 && TMP_SC > SSX2(k, p7S_GCAG, t, u)) {
          SSX2(k, p7S_GCAG, t, u) = TMP_SC;
          SIX2(k, p7S_GCAG, t, u) = i-MIN_INTRON_LENG;
        }

        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-3,k-1), DMX_SP(i-MIN_INTRON_LENG-3,k-1));
        if(u < 4 && TMP_SC > SSX1(k, p7S_GCAG, u)) {
          SSX1(k, p7S_GCAG, u) = TMP_SC;
          SIX1(k, p7S_GCAG, u) = i-MIN_INTRON_LENG;
        }

        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-2,k-1), DMX_SP(i-MIN_INTRON_LENG-2,k-1));
        if(TMP_SC > SSX0(k, p7S_GCAG)) {
          SSX0(k, p7S_GCAG) = TMP_SC;
          SIX0(k, p7S_GCAG) = i-MIN_INTRON_LENG;
        }
      }
    }
    else if(AT) {
      for (k = 2; k < M; k++) {
        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-4,k-1), DMX_SP(i-MIN_INTRON_LENG-4,k-1));
        if(t < 4 && u < 4 && TMP_SC > SSX2(k, p7S_ATAC, t, u)) {
          SSX2(k, p7S_ATAC, t, u) = TMP_SC;
          SIX2(k, p7S_ATAC, t, u) = i-MIN_INTRON_LENG;
        }

        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-3,k-1), DMX_SP(i-MIN_INTRON_LENG-3,k-1));
        if(u < 4 && TMP_SC > SSX1(k, p7S_ATAC, u)) {
          SSX1(k, p7S_ATAC, u) = TMP_SC;
          SIX1(k, p7S_ATAC, u) = i-MIN_INTRON_LENG;
        }

        TMP_SC = ESL_MAX(MMX_SP(i-MIN_INTRON_LENG-2,k-1), DMX_SP(i-MIN_INTRON_LENG-2,k-1));
        if(TMP_SC > SSX0(k, p7S_ATAC)) {
          SSX0(k, p7S_ATAC) = TMP_SC;
          SIX0(k, p7S_ATAC) = i-MIN_INTRON_LENG;
        }
      }
    }
 
    /* create 3 frames where only special states or the intron portion of a P state can emit */
    if(r_i < remove_cnt         &&  
       i >= remove_idx[r_i] - 2 &&
       i <= remove_idx[r_i] + 2) {
      for(k = 1; k <= M; k++)
        MMX_SP(i,k) = IMX_SP(i,k) = DMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;
             
      if(i == remove_idx[r_i] + 2) r_i++;
       //printf("r_i %d i %d \n", r_i, i); 
    }
   
  } // end loop over L
  
  gx->M = gm_fs->M;
  gx->L = L;

//p7_gmx_sp_Dump(stdout, gx, p7_DEFAULT);
  return eslOK;

}

/* MAX lookback is i-MIN_INTRON_LENG-4 */
//Global on both ends of model (must enter at 1 and exit at M) 
//Used to find long introns and determine right end termination on both sequence and model 
int
p7_spliceviterbi_parser_semiglobal(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx)
{
  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float      **score = pli->sig_idx->score;
  float       *signal_scores = pli->signal_scores;
  float       *parser_scores = pli->sig_idx->parser_scores;
  int        **index = pli->sig_idx->index;
  int        **lookback = pli->sig_idx->lookback;
  int         *parser_index = pli->sig_idx->parser_index;
  int          M    = gm_fs->M;
  int          i,k,s;
  int          c3, c2, c1;
  int          t, u, v, w, x;
  int          AGXXX, ACXXX;
  int          AGXX, ACXX;
  int          AGX,  ACX;
  int          nuc1, nuc2;
  int          curr_i, prev_i;
  int          donor_i;
  float        TMP_SC; 
 
  /*Initialize the signal index */
  for(k = 0; k < M; k++) 
    lookback[0][k] = -1;
  
  for(i = 0; i <= L; i++) { 
    PI(i,p7S_PI) = PI(i,p7S_PK) = PI(i,p7S_MK) = -1; 
    PS(i,p7S_P) =  PS(i,p7S_M) = -eslINFINITY;  
  }

  for(k = 0; k <  M; k++) {
    for(s = 0; s < SIGNAL_MEM_SIZE; s++) {
      index[k][s] = -1;
      score[k][s] = -eslINFINITY;
    }
  }

  /* Initialization of the zero row.  */
  XMX_SP(0,p7G_N) = 0.;                                     /* S->N, p=1            */
  XMX_SP(0,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE];            /* S->N->B, no N-tail   */ 
  XMX_SP(0,p7G_E) = XMX_SP(0,p7G_C) = XMX_SP(0,p7G_J) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_SP(0,k) = DMX_SP(0,k) = IMX_SP(0,k) = PMX_SP(0,k) = -eslINFINITY;
  
  //p7_gmx_sp_DumpHeader(stdout, gx, 0, gm_fs->M, p7_DEFAULT);
  //p7_gmx_sp_DumpRow(stdout, gx, 0, 0, 0, gm_fs->M, p7_DEFAULT);
  /*Special cases for the first 2 rows */
  t = u = v = w = x = -1;
  for(i = 1; i <= 2; i++)
  {
    w = x;
  
    /* if new nucleotide is not A,C,G, or T set it to placeholder vlaue */
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                                x = p7P_MAXCODONS;

    XMX_SP(i,p7G_N) = 0; 
    XMX_SP(i,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE]; 

    curr_i = i % (MIN_INTRON_LENG+5);
    for (k = 0; k <= M; k++) 
      MMX_SP(curr_i,k) = DMX_SP(curr_i,k) = IMX_SP(curr_i,k) = PMX_SP(curr_i,k) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;

  // p7_gmx_sp_DumpRow(stdout, gx, i, curr_i, 0, gm_fs->M, p7_DEFAULT);
  } 

  /*Special cases for the rows 3-MIN_INTRON_LENG+4 (no P state)*/
  for (i = 3; i <= MIN_INTRON_LENG+4; i++) 
  {
    v = w;
    w = x;

    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                                x = p7P_MAXCODONS;

    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    XMX_SP(i,p7G_N) = XMX_SP(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];
    XMX_SP(i,p7G_B) = XMX_SP(i,p7G_N)   + gm_fs->xsc[p7P_N][p7P_MOVE];

    curr_i = i % (MIN_INTRON_LENG+5);
    prev_i = (i-3) % (MIN_INTRON_LENG+5);

    MMX_SP(curr_i,0) = IMX_SP(curr_i,0) = DMX_SP(curr_i,0) = PMX_SP(curr_i,0) = -eslINFINITY;
    
    MMX_SP(curr_i,1) = XMX_SP(i-3,p7G_B) + p7P_MSC_CODON(gm_fs, 1, c3);

    if(MMX_SP(curr_i,1) > PS(i,p7S_M)) {
      PS(i,p7S_M)  = MMX_SP(curr_i,1);
      PI(i,p7S_MK) = 1;
    }

    if(p7P_MSC_CODON(gm_fs, 1, c3) == -eslINFINITY)
      IMX_SP(curr_i,1) = -eslINFINITY;
    else
      IMX_SP(curr_i,1) = ESL_MAX(MMX_SP(prev_i,1) + TSC(p7P_MI,1),
                                 IMX_SP(prev_i,1) + TSC(p7P_II,1));
    
    DMX_SP(curr_i,1) = PMX_SP(curr_i,1) = -eslINFINITY; 

    for (k = 2; k < M; k++) {
      
      MMX_SP(curr_i,k) = ESL_MAX(MMX_SP(prev_i,k-1) + TSC(p7P_MM,k-1),
                         ESL_MAX(IMX_SP(prev_i,k-1) + TSC(p7P_IM,k-1),
                                 DMX_SP(prev_i,k-1) + TSC(p7P_DM,k-1))) + p7P_MSC_CODON(gm_fs, k, c3);

      if(MMX_SP(curr_i,k) > PS(i,p7S_M)) {
        PS(i,p7S_M)  = MMX_SP(curr_i,k);
        PI(i,p7S_MK) = k;
      }      

      if(p7P_MSC_CODON(gm_fs, k, c3) == -eslINFINITY) 
        IMX_SP(curr_i,k) = -eslINFINITY;
      else                                            
        IMX_SP(curr_i,k) = ESL_MAX(MMX_SP(prev_i,k) + TSC(p7P_MI,k),
                                   IMX_SP(prev_i,k) + TSC(p7P_II,k));

      DMX_SP(curr_i,k) = ESL_MAX(MMX_SP(curr_i,k-1) + TSC(p7P_MD,k-1),
                                 DMX_SP(curr_i,k-1) + TSC(p7P_DD,k-1));

      PMX_SP(curr_i,k) = -eslINFINITY;
    } 

    MMX_SP(curr_i,M) = ESL_MAX(MMX_SP(prev_i,M-1) + TSC(p7P_MM,M-1),
                       ESL_MAX(IMX_SP(prev_i,M-1) + TSC(p7P_IM,M-1),
                               DMX_SP(prev_i,M-1) + TSC(p7P_DM,M-1))) + p7P_MSC_CODON(gm_fs, M, c3);
    
    if(MMX_SP(curr_i,M) > PS(i,p7S_M)) {
      PS(i,p7S_M)  = MMX_SP(curr_i,M);
      PI(i,p7S_MK) = M;
    }

    IMX_SP(curr_i,M) = -eslINFINITY;

    DMX_SP(curr_i,M) = ESL_MAX(MMX_SP(curr_i,M-1) + TSC(p7P_MD,M-1),
                               DMX_SP(curr_i,M-1) + TSC(p7P_DD,M-1));

    PMX_SP(curr_i,M) = -eslINFINITY;
 
    XMX_SP(i,p7G_E) = ESL_MAX(MMX_SP(curr_i,M), DMX_SP(curr_i,M)); 
    
    XMX_SP(i,p7G_C) = ESL_MAX(XMX_SP(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                              XMX_SP(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);
    
    XMX_SP(i,p7G_J) = -eslINFINITY; 
    //p7_gmx_sp_DumpRow(stdout, gx, i, curr_i, 0, gm_fs->M, p7_DEFAULT);
  }
 
  AGXXX = ACXXX = AGXX = ACXX = AGX = ACX  = FALSE;  
  /* Main DP recursion */
  for (i = MIN_INTRON_LENG+5; i <= L; i++) 
  {
    /* get nucleotides and codon */
    v = w;
    w = x;

    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                            x = p7P_MAXCODONS;
    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);
 
    XMX_SP(i,p7G_N) = XMX_SP(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];
    XMX_SP(i,p7G_B) = XMX_SP(i,p7G_N)   + gm_fs->xsc[p7P_N][p7P_MOVE];

    curr_i = i % (MIN_INTRON_LENG+5);
    prev_i = (i-3) % (MIN_INTRON_LENG+5);
   
    MMX_SP(curr_i,0) = IMX_SP(curr_i,0) = DMX_SP(curr_i,0) = PMX_SP(curr_i,0) = -eslINFINITY;

    MMX_SP(curr_i,1) = XMX_SP(i-3,p7G_B) + p7P_MSC_CODON(gm_fs, 1, c3);

    if(MMX_SP(curr_i,1) > PS(i,p7S_M)) {
      PS(i,p7S_M)  = MMX_SP(curr_i,1);
      PI(i,p7S_MK) = 1;
    }

    if(p7P_MSC_CODON(gm_fs, 1, c3) == -eslINFINITY)
      IMX_SP(curr_i,1) = -eslINFINITY;
    else
      IMX_SP(curr_i,1) = ESL_MAX(MMX_SP(prev_i,1) + TSC(p7P_MI,1),
                                 IMX_SP(prev_i,1) + TSC(p7P_II,1));
    
    DMX_SP(curr_i,1) = PMX_SP(curr_i,1) = -eslINFINITY; 

    for (k = 2; k < M; k++) {

      MMX_SP(curr_i,k) = ESL_MAX(MMX_SP(prev_i,k-1) + TSC(p7P_MM,k-1),
                         ESL_MAX(IMX_SP(prev_i,k-1) + TSC(p7P_IM,k-1),
                         ESL_MAX(DMX_SP(prev_i,k-1) + TSC(p7P_DM,k-1),
                                 PMX_SP(prev_i,k-1) + TSC_P))) + p7P_MSC_CODON(gm_fs, k, c3);

      if(MMX_SP(curr_i,k) > PS(i,p7S_M)) {
        PS(i,p7S_M)  = MMX_SP(curr_i,k);
        PI(i,p7S_MK) = k;
      }     

      if(p7P_MSC_CODON(gm_fs, k, c3) == -eslINFINITY) 
        IMX_SP(curr_i,k) = -eslINFINITY;
      else                                            
        IMX_SP(curr_i,k) = ESL_MAX(MMX_SP(prev_i,k) + TSC(p7P_MI,k),
                                   IMX_SP(prev_i,k) + TSC(p7P_II,k));

      DMX_SP(curr_i,k) = ESL_MAX(MMX_SP(curr_i,k-1) + TSC(p7P_MD,k-1),
                                 DMX_SP(curr_i,k-1) + TSC(p7P_DD,k-1));

      PMX_SP(curr_i,k) = -eslINFINITY;
    }  
     
    MMX_SP(curr_i,M) = ESL_MAX(MMX_SP(prev_i,M-1) + TSC(p7P_MM,M-1),
                       ESL_MAX(IMX_SP(prev_i,M-1) + TSC(p7P_IM,M-1),
                       ESL_MAX(DMX_SP(prev_i,M-1) + TSC(p7P_DM,M-1),
                               PMX_SP(prev_i,M-1) + TSC_P)))         + p7P_MSC_CODON(gm_fs, M, c3);

    if(MMX_SP(curr_i,M) > PS(i,p7S_M)) {
        PS(i,p7S_M)  = MMX_SP(curr_i,M);
        PI(i,p7S_MK) = M;
    }

    IMX_SP(curr_i,M) = -eslINFINITY;

    DMX_SP(curr_i,M) = ESL_MAX(MMX_SP(curr_i,M-1) + TSC(p7P_MD,M-1),
                               DMX_SP(curr_i,M-1) + TSC(p7P_DD,M-1));

    PMX_SP(curr_i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = ESL_MAX(MMX_SP(curr_i,M), DMX_SP(curr_i,M));

    XMX_SP(i,p7G_C) = ESL_MAX(XMX_SP(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                              XMX_SP(i,  p7G_E) + gm_fs->xsc[p7P_E][p7P_MOVE]);

    XMX_SP(i,p7G_J) = -eslINFINITY; 

    /* get acceptor site */
    AGXXX = AGXX;
    AGXX  = AGX;
    if(SIGNAL(v, w) == ACCEPT_AG) AGX = TRUE;
    else                          AGX = FALSE;
 
    ACXXX = ACXX;
    ACXX  = ACX;
    if(SIGNAL(v, w) == ACCEPT_AC) ACX = TRUE;
    else                          ACX = FALSE;

    /* Sperate K loops for all acccptor and donor sites to prevent conditional in the inner loops */
    if(AGXXX) {
      for (k = 2; k < M; k++) { 
        TMP_SC = SSX0(k, p7S_GTAG) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c3);
        if(TMP_SC > PMX_SP(curr_i,k)) {
          PMX_SP(curr_i,k) = TMP_SC;
          lookback[0][k] = SIX0(k, p7S_GTAG);
        }
  
        TMP_SC = SSX0(k, p7S_GCAG) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c3);
        if(TMP_SC > PMX_SP(curr_i,k)) {
          PMX_SP(curr_i,k) = TMP_SC;
          lookback[0][k] = SIX0(k, p7S_GCAG);
        }
      }
    }
    else if(ACXXX) {
      for (k = 2; k < M; k++) {
        TMP_SC = SSX0(k, p7S_ATAC) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c3);
        if(TMP_SC > PMX_SP(curr_i,k)) {
          PMX_SP(curr_i,k) = TMP_SC;
          lookback[0][k] = SIX0(k, p7S_ATAC);
        }
      }
    }

    if(AGXX) {
      for (k = 2; k < M; k++) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(nuc1, w, x);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);
       
          TMP_SC = SSX1(k, p7S_GTAG, nuc1) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c2); 
          if(TMP_SC > PMX_SP(curr_i,k)) {
            PMX_SP(curr_i,k) = TMP_SC;
            lookback[0][k] = SIX1(k, p7S_GTAG, nuc1);
          }
  
          TMP_SC = SSX1(k, p7S_GCAG, nuc1) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c2);
          if(TMP_SC > PMX_SP(curr_i,k)) {
            PMX_SP(curr_i,k) = TMP_SC;
            lookback[0][k] = SIX1(k, p7S_GCAG, nuc1);
          }
        }
      }
    }
    else if(ACXX) {
      for (k = 2; k < M; k++) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(nuc1, w, x);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);
  
          TMP_SC = SSX1(k, p7S_ATAC, nuc1) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c2);
          if(TMP_SC > PMX_SP(curr_i,k)) {
            PMX_SP(curr_i,k) = TMP_SC;
            lookback[0][k] = SIX1(k, p7S_ATAC, nuc1);
          }
        }
      }
    }
  
    if(AGX) {
      for (k = 2; k < M; k++) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(nuc1, nuc2, x);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);
  
            TMP_SC = SSX2(k, p7S_GTAG, nuc1, nuc2) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c1);
            if(TMP_SC > PMX_SP(curr_i,k)) {
              PMX_SP(curr_i,k) = TMP_SC;
              lookback[0][k] = SIX2(k, p7S_GTAG, nuc1, nuc2);
            }
  
            TMP_SC = SSX2(k, p7S_GCAG, nuc1, nuc2) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c1);
            if(TMP_SC > PMX_SP(curr_i,k)) {
              PMX_SP(curr_i,k) = TMP_SC;
              lookback[0][k] = SIX2(k, p7S_GCAG, nuc1, nuc2);
            }
          }
        }
      }
    }
    else if(ACX) {
      for (k = 2; k < M; k++) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(nuc1, nuc2, x);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);
  
            TMP_SC = SSX2(k, p7S_ATAC, nuc1, nuc2) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c1);
            if(TMP_SC > PMX_SP(curr_i,k)) {
              PMX_SP(curr_i,k) = TMP_SC;
              lookback[0][k] = SIX2(k, p7S_ATAC, nuc1, nuc2);
            }
          }
        }
      }
    }

    for (k = 2; k < M; k++) { 
      if(PMX_SP(curr_i,k) > PS(i,p7S_P)) {
        PS(i,p7S_P)  = PMX_SP(curr_i,k);
        PI(i,p7S_PI) = lookback[0][k]; 
        PI(i,p7S_PK) = k;
      }
    }

    /* Record best scoring donor sites*/
    if(SIGNAL(sub_dsq[i-MIN_INTRON_LENG-1], sub_dsq[i-MIN_INTRON_LENG]) == DONOR_GT) {
      t = sub_dsq[i-MIN_INTRON_LENG-3];
      u = sub_dsq[i-MIN_INTRON_LENG-2];
      if(t < 4 && u < 4) {
        for (k = 2; k < M; k++) {
          donor_i = (i-MIN_INTRON_LENG-4) % (MIN_INTRON_LENG+5); 
          TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1)); 
          if(TMP_SC > SSX2(k, p7S_GTAG, t, u)) {
            SSX2(k, p7S_GTAG, t, u) = TMP_SC;
            SIX2(k, p7S_GTAG, t, u) = i-MIN_INTRON_LENG;
          }
        }
      }
      
      if(u < 4) {
        for (k = 2; k < M; k++) {
          donor_i = (i-MIN_INTRON_LENG-3) % (MIN_INTRON_LENG+5);
          TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
          if(TMP_SC > SSX1(k, p7S_GTAG, u)) {
            SSX1(k, p7S_GTAG, u) = TMP_SC;
            SIX1(k, p7S_GTAG, u) = i-MIN_INTRON_LENG; 
          } 
        }
      }
      
      for (k = 2; k < M; k++) {
        donor_i = (i-MIN_INTRON_LENG-2) % (MIN_INTRON_LENG+5);
        TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
        if(TMP_SC > SSX0(k, p7S_GTAG)) {
          SSX0(k, p7S_GTAG) = TMP_SC;
          SIX0(k, p7S_GTAG) = i-MIN_INTRON_LENG;
        }
      } 
    }
    else if(SIGNAL(sub_dsq[i-MIN_INTRON_LENG-1], sub_dsq[i-MIN_INTRON_LENG]) == DONOR_GC) {
      t = sub_dsq[i-MIN_INTRON_LENG-3];
      u = sub_dsq[i-MIN_INTRON_LENG-2];
      if(t < 4 && u < 4) {
        for (k = 2; k < M; k++) {
          donor_i = (i-MIN_INTRON_LENG-4) % (MIN_INTRON_LENG+5);
          TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
          if(TMP_SC > SSX2(k, p7S_GCAG, t, u)) {
            SSX2(k, p7S_GCAG, t, u) = TMP_SC;
            SIX2(k, p7S_GCAG, t, u) = i-MIN_INTRON_LENG;
          }
        }
      }

      if(u < 4) {
        for (k = 2; k < M; k++) {
          donor_i = (i-MIN_INTRON_LENG-3) % (MIN_INTRON_LENG+5);
          TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
          if(TMP_SC > SSX1(k, p7S_GCAG, u)) {
            SSX1(k, p7S_GCAG, u) = TMP_SC;
            SIX1(k, p7S_GCAG, u) = i-MIN_INTRON_LENG;
          }
        }
      }

      for (k = 2; k < M; k++) {
        donor_i = (i-MIN_INTRON_LENG-2) % (MIN_INTRON_LENG+5);
        TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
        if(TMP_SC > SSX0(k, p7S_GCAG)) {
          SSX0(k, p7S_GCAG) = TMP_SC;
          SIX0(k, p7S_GCAG) = i-MIN_INTRON_LENG;
        }
      }    
    }
    else if(SIGNAL(sub_dsq[i-MIN_INTRON_LENG-1], sub_dsq[i-MIN_INTRON_LENG]) == DONOR_AT) {
      t = sub_dsq[i-MIN_INTRON_LENG-3];
      u = sub_dsq[i-MIN_INTRON_LENG-2];
      if(t < 4 && u < 4) {
        for (k = 2; k < M; k++) {
          donor_i = (i-MIN_INTRON_LENG-4) % (MIN_INTRON_LENG+5);
          TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
          if(TMP_SC > SSX2(k, p7S_ATAC, t, u)) {
            SSX2(k, p7S_ATAC, t, u) = TMP_SC;
            SIX2(k, p7S_ATAC, t, u) = i-MIN_INTRON_LENG;
          }
        }
      }

      if(u < 4) {
        for (k = 2; k < M; k++) {
          donor_i = (i-MIN_INTRON_LENG-3) % (MIN_INTRON_LENG+5);
          TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
          if(TMP_SC > SSX1(k, p7S_ATAC, u)) {
            SSX1(k, p7S_ATAC, u) = TMP_SC;
            SIX1(k, p7S_ATAC, u) = i-MIN_INTRON_LENG;
          }
        }
      }

      for (k = 2; k < M; k++) {
        donor_i = (i-MIN_INTRON_LENG-2) % (MIN_INTRON_LENG+5);
        TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
        if(TMP_SC > SSX0(k, p7S_ATAC)) {
          SSX0(k, p7S_ATAC) = TMP_SC;
          SIX0(k, p7S_ATAC) = i-MIN_INTRON_LENG;
        }
      }
    }

        //p7_gmx_sp_DumpRow(stdout, gx, i, curr_i, 0, gm_fs->M, p7_DEFAULT);
  } // end loop over L
 
  gx->M = gm_fs->M;
  gx->L = L;

  return eslOK;

}



/* MAX lookback is i-MIN_INTRON_LENG-4 */
//Global on left end of model (must enter at 1) 
//Used to find long introns and determine right end termination on both sequence and model 
int
p7_spliceviterbi_rightparser_semiglobal(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx)
{
  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float      **score = pli->sig_idx->score;
  float       *signal_scores = pli->signal_scores;
  float       *parser_scores = pli->sig_idx->parser_scores;
  int        **index = pli->sig_idx->index;
  int        **lookback = pli->sig_idx->lookback;
  int         *parser_index = pli->sig_idx->parser_index;
  int          M    = gm_fs->M;
  int          i,k,s;
  int          c3, c2, c1;
  int          t, u, v, w, x;
  int          AGXXX, ACXXX;
  int          AGXX, ACXX;
  int          AGX,  ACX;
  int          nuc1, nuc2;
  int          curr_i, prev_i;
  int          donor_i;
  float        TMP_SC; 
 
  /*Initialize the signal index */
  for(k = 0; k < M; k++) 
    lookback[0][k] = -1;
  
  for(i = 0; i <= L; i++) { 
    PI(i,p7S_PI) = PI(i,p7S_PK) = PI(i,p7S_MK) = -1;
    PS(i,p7S_P) =  PS(i,p7S_M) = -eslINFINITY;
  }

  for(k = 0; k <  M; k++) {
    for(s = 0; s < SIGNAL_MEM_SIZE; s++) {
      index[k][s] = -1;
      score[k][s] = -eslINFINITY;
    }
  }

  /* Initialization of the zero row.  */
  XMX_SP(0,p7G_N) = 0.;                                     /* S->N, p=1            */
  XMX_SP(0,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE];            /* S->N->B, no N-tail   */ 
  XMX_SP(0,p7G_E) = XMX_SP(0,p7G_C) = XMX_SP(0,p7G_J) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_SP(0,k) = DMX_SP(0,k) = IMX_SP(0,k) = PMX_SP(0,k) = -eslINFINITY;
  
  //p7_gmx_sp_DumpHeader(stdout, gx, 0, gm_fs->M, p7_DEFAULT);
  //p7_gmx_sp_DumpRow(stdout, gx, 0, 0, 0, gm_fs->M, p7_DEFAULT);
  /*Special cases for the first 2 rows */
  t = u = v = w = x = -1;
  for(i = 1; i <= 2; i++)
  {
    w = x;
  
    /* if new nucleotide is not A,C,G, or T set it to placeholder vlaue */
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                                x = p7P_MAXCODONS;

    XMX_SP(i,p7G_N) = 0; 
    XMX_SP(i,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE]; 

    curr_i = i % (MIN_INTRON_LENG+5);
    for (k = 0; k <= M; k++) 
      MMX_SP(curr_i,k) = DMX_SP(curr_i,k) = IMX_SP(curr_i,k) = PMX_SP(curr_i,k) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;

   //p7_gmx_sp_DumpRow(stdout, gx, i, curr_i, 0, gm_fs->M, p7_DEFAULT);
  } 

  /*Special cases for the rows 3-MIN_INTRON_LENG+4 (no P state)*/
  for (i = 3; i <= MIN_INTRON_LENG+4; i++) 
  {
    v = w;
    w = x;

    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                                x = p7P_MAXCODONS;

    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    XMX_SP(i,p7G_N) = XMX_SP(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];
    XMX_SP(i,p7G_B) = XMX_SP(i,p7G_N)   + gm_fs->xsc[p7P_N][p7P_MOVE];

    curr_i = i % (MIN_INTRON_LENG+5);
    prev_i = (i-3) % (MIN_INTRON_LENG+5);

    MMX_SP(curr_i,0) = IMX_SP(curr_i,0) = DMX_SP(curr_i,0) = PMX_SP(curr_i,0) = -eslINFINITY;
    
    MMX_SP(curr_i,1) = XMX_SP(i-3,p7G_B) + p7P_MSC_CODON(gm_fs, 1, c3);

    if(MMX_SP(curr_i,1) > PS(i,p7S_M)) {
      PS(i,p7S_M)  = MMX_SP(curr_i,1);
      PI(i,p7S_MK) = 1;
    }

    if(p7P_MSC_CODON(gm_fs, 1, c3) == -eslINFINITY)
      IMX_SP(curr_i,1) = -eslINFINITY;
    else
      IMX_SP(curr_i,1) = ESL_MAX(MMX_SP(prev_i,1) + TSC(p7P_MI,1),
                                 IMX_SP(prev_i,1) + TSC(p7P_II,1));
    
    DMX_SP(curr_i,1) = PMX_SP(curr_i,1) = -eslINFINITY; 

    XMX_SP(i,p7G_E) = MMX_SP(curr_i,1); 

    for (k = 2; k < M; k++) {
      
      MMX_SP(curr_i,k) = ESL_MAX(MMX_SP(prev_i,k-1) + TSC(p7P_MM,k-1),
                         ESL_MAX(IMX_SP(prev_i,k-1) + TSC(p7P_IM,k-1),
                                 DMX_SP(prev_i,k-1) + TSC(p7P_DM,k-1))) + p7P_MSC_CODON(gm_fs, k, c3);

      if(MMX_SP(curr_i,k) > PS(i,p7S_M)) {
        PS(i,p7S_M)  = MMX_SP(curr_i,k);
        PI(i,p7S_MK) = k;
      }      

      if(p7P_MSC_CODON(gm_fs, k, c3) == -eslINFINITY) 
        IMX_SP(curr_i,k) = -eslINFINITY;
      else                                            
        IMX_SP(curr_i,k) = ESL_MAX(MMX_SP(prev_i,k) + TSC(p7P_MI,k),
                                   IMX_SP(prev_i,k) + TSC(p7P_II,k));

      DMX_SP(curr_i,k) = ESL_MAX(MMX_SP(curr_i,k-1) + TSC(p7P_MD,k-1),
                                 DMX_SP(curr_i,k-1) + TSC(p7P_DD,k-1));

      PMX_SP(curr_i,k) = -eslINFINITY;

      XMX_SP(i,p7G_E) = ESL_MAX(XMX_SP(i,p7G_E),
                        ESL_MAX(MMX_SP(curr_i,k), DMX_SP(curr_i,k)));
    } 

    MMX_SP(curr_i,M) = ESL_MAX(MMX_SP(prev_i,M-1) + TSC(p7P_MM,M-1),
                       ESL_MAX(IMX_SP(prev_i,M-1) + TSC(p7P_IM,M-1),
                               DMX_SP(prev_i,M-1) + TSC(p7P_DM,M-1))) + p7P_MSC_CODON(gm_fs, M, c3);

    if(MMX_SP(curr_i,M) > PS(i,p7S_M)) {
      PS(i,p7S_M)  = MMX_SP(curr_i,M);
      PI(i,p7S_MK) = M;
    }

    IMX_SP(curr_i,M) = -eslINFINITY;

    DMX_SP(curr_i,M) = ESL_MAX(MMX_SP(curr_i,M-1) + TSC(p7P_MD,M-1),
                               DMX_SP(curr_i,M-1) + TSC(p7P_DD,M-1));

    PMX_SP(curr_i,M) = -eslINFINITY;
 
    XMX_SP(i,p7G_E) = ESL_MAX(XMX_SP(i,p7G_E),
                      ESL_MAX(MMX_SP(curr_i,M), DMX_SP(curr_i,M))); 
    
    XMX_SP(i,p7G_C) = ESL_MAX(XMX_SP(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                              XMX_SP(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);
    
    XMX_SP(i,p7G_J) = -eslINFINITY; 
//    p7_gmx_sp_DumpRow(stdout, gx, i, curr_i, 0, gm_fs->M, p7_DEFAULT);
  }
 
  AGXXX = ACXXX = AGXX = ACXX = AGX = ACX  = FALSE;  
  /* Main DP recursion */
  for (i = MIN_INTRON_LENG+5; i <= L; i++) 
  {
    /* get nucleotides and codon */
    v = w;
    w = x;

    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
    else                                            x = p7P_MAXCODONS;
    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);
 
    XMX_SP(i,p7G_N) = XMX_SP(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];
    XMX_SP(i,p7G_B) = XMX_SP(i,p7G_N)   + gm_fs->xsc[p7P_N][p7P_MOVE];

    curr_i = i % (MIN_INTRON_LENG+5);
    prev_i = (i-3) % (MIN_INTRON_LENG+5);
   
    MMX_SP(curr_i,0) = IMX_SP(curr_i,0) = DMX_SP(curr_i,0) = PMX_SP(curr_i,0) = -eslINFINITY;

    MMX_SP(curr_i,1) = XMX_SP(i-3,p7G_B) + p7P_MSC_CODON(gm_fs, 1, c3);

    if(MMX_SP(curr_i,1) > PS(i,p7S_M)) {
      PS(i,p7S_M)  = MMX_SP(curr_i,1);
      PI(i,p7S_MK) = 1;
    }

    if(p7P_MSC_CODON(gm_fs, 1, c3) == -eslINFINITY)
      IMX_SP(curr_i,1) = -eslINFINITY;
    else
      IMX_SP(curr_i,1) = ESL_MAX(MMX_SP(prev_i,1) + TSC(p7P_MI,1),
                                 IMX_SP(prev_i,1) + TSC(p7P_II,1));
    
    DMX_SP(curr_i,1) = PMX_SP(curr_i,1) = -eslINFINITY; 

    XMX_SP(i,p7G_E) = MMX_SP(curr_i,1); 

    for (k = 2; k < M; k++) {

      MMX_SP(curr_i,k) = ESL_MAX(MMX_SP(prev_i,k-1) + TSC(p7P_MM,k-1),
                         ESL_MAX(IMX_SP(prev_i,k-1) + TSC(p7P_IM,k-1),
                         ESL_MAX(DMX_SP(prev_i,k-1) + TSC(p7P_DM,k-1),
                                 PMX_SP(prev_i,k-1) + TSC_P))) + p7P_MSC_CODON(gm_fs, k, c3);
    
      if(MMX_SP(curr_i,k) > PS(i,p7S_M)) {
        PS(i,p7S_M)  = MMX_SP(curr_i,k);
        PI(i,p7S_MK) = k;
      } 

      if(p7P_MSC_CODON(gm_fs, k, c3) == -eslINFINITY) 
        IMX_SP(curr_i,k) = -eslINFINITY;
      else                                            
        IMX_SP(curr_i,k) = ESL_MAX(MMX_SP(prev_i,k) + TSC(p7P_MI,k),
                                   IMX_SP(prev_i,k) + TSC(p7P_II,k));

      DMX_SP(curr_i,k) = ESL_MAX(MMX_SP(curr_i,k-1) + TSC(p7P_MD,k-1),
                                 DMX_SP(curr_i,k-1) + TSC(p7P_DD,k-1));

      PMX_SP(curr_i,k) = -eslINFINITY;

      XMX_SP(i,p7G_E) = ESL_MAX(XMX_SP(i,p7G_E),
                        ESL_MAX(MMX_SP(curr_i,k), DMX_SP(curr_i,k)));
    }  
     
    MMX_SP(curr_i,M) = ESL_MAX(MMX_SP(prev_i,M-1) + TSC(p7P_MM,M-1),
                       ESL_MAX(IMX_SP(prev_i,M-1) + TSC(p7P_IM,M-1),
                       ESL_MAX(DMX_SP(prev_i,M-1) + TSC(p7P_DM,M-1),
                               PMX_SP(prev_i,M-1) + TSC_P)))         + p7P_MSC_CODON(gm_fs, M, c3);

    if(MMX_SP(curr_i,M) > PS(i,p7S_M)) {
      PS(i,p7S_M)  = MMX_SP(curr_i,M);
      PI(i,p7S_MK) = M;
    }

    IMX_SP(curr_i,M) = -eslINFINITY;

    DMX_SP(curr_i,M) = ESL_MAX(MMX_SP(curr_i,M-1) + TSC(p7P_MD,M-1),
                               DMX_SP(curr_i,M-1) + TSC(p7P_DD,M-1));

    PMX_SP(curr_i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = ESL_MAX(XMX_SP(i,p7G_E),
                      ESL_MAX(MMX_SP(curr_i,M), DMX_SP(curr_i,M)));

    XMX_SP(i,p7G_C) = ESL_MAX(XMX_SP(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                              XMX_SP(i,  p7G_E) + gm_fs->xsc[p7P_E][p7P_MOVE]);

    XMX_SP(i,p7G_J) = -eslINFINITY; 

    /* get acceptor site */
    AGXXX = AGXX;
    ACXXX = ACXX;
    if(SIGNAL(v, w) == ACCEPT_AG) AGX = TRUE;
    else                          AGX = FALSE;
 
    AGXX  = AGX;
    ACXX  = ACX;
    if(SIGNAL(v, w) == ACCEPT_AC) ACX = TRUE;
    else                          ACX = FALSE;

    
    /* Sperate K loops for all acccptor and donor sites to prevent conditional in the inner loops */
    if(AGXXX) {
      for (k = 2; k < M; k++) { 
        TMP_SC = SSX0(k, p7S_GTAG) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c3);
        if(TMP_SC > PMX_SP(curr_i,k)) {
          PMX_SP(curr_i,k) = TMP_SC;
          lookback[0][k] = SIX0(k, p7S_GTAG);
        }
  
        TMP_SC = SSX0(k, p7S_GCAG) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c3);
        if(TMP_SC > PMX_SP(curr_i,k)) {
          PMX_SP(curr_i,k) = TMP_SC;
          lookback[0][k] = SIX0(k, p7S_GCAG);
        }
      }
    }
    else if(ACXXX) {
      for (k = 2; k < M; k++) {
        TMP_SC = SSX0(k, p7S_ATAC) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c3);
        if(TMP_SC > PMX_SP(curr_i,k)) {
          PMX_SP(curr_i,k) = TMP_SC;
          lookback[0][k] = SIX0(k, p7S_ATAC);
        }
      }
    }

    if(AGXX) {
      for (k = 2; k < M; k++) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(nuc1, w, x);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);
       
          TMP_SC = SSX1(k, p7S_GTAG, nuc1) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c2); 
          if(TMP_SC > PMX_SP(curr_i,k)) {
            PMX_SP(curr_i,k) = TMP_SC;
            lookback[0][k] = SIX1(k, p7S_GTAG, nuc1);
          }
  
          TMP_SC = SSX1(k, p7S_GCAG, nuc1) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c2);
          if(TMP_SC > PMX_SP(curr_i,k)) {
            PMX_SP(curr_i,k) = TMP_SC;
            lookback[0][k] = SIX1(k, p7S_GCAG, nuc1);
          }
        }
      }
    }
    else if(ACXX) {
      for (k = 2; k < M; k++) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(nuc1, w, x);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);
  
          TMP_SC = SSX1(k, p7S_ATAC, nuc1) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c2);
          if(TMP_SC > PMX_SP(curr_i,k)) {
            PMX_SP(curr_i,k) = TMP_SC;
            lookback[0][k] = SIX1(k, p7S_ATAC, nuc1);
          }
        }
      }
    }
  
    if(AGX) {
      for (k = 2; k < M; k++) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(nuc1, nuc2, x);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);
  
            TMP_SC = SSX2(k, p7S_GTAG, nuc1, nuc2) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c1);
            if(TMP_SC > PMX_SP(curr_i,k)) {
              PMX_SP(curr_i,k) = TMP_SC;
              lookback[0][k] = SIX2(k, p7S_GTAG, nuc1, nuc2);
            }
  
            TMP_SC = SSX2(k, p7S_GCAG, nuc1, nuc2) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c1);
            if(TMP_SC > PMX_SP(curr_i,k)) {
              PMX_SP(curr_i,k) = TMP_SC;
              lookback[0][k] = SIX2(k, p7S_GCAG, nuc1, nuc2);
            }
          }
        }
      }
    }
    else if(ACX) {
      for (k = 2; k < M; k++) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(nuc1, nuc2, x);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);
  
            TMP_SC = SSX2(k, p7S_ATAC, nuc1, nuc2) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c1);
            if(TMP_SC > PMX_SP(curr_i,k)) {
              PMX_SP(curr_i,k) = TMP_SC;
              lookback[0][k] = SIX2(k, p7S_ATAC, nuc1, nuc2);
            }
          }
        }
      }
    }

    for (k = 2; k < M; k++) {
      if(PMX_SP(curr_i,k) > PS(i,p7S_P)) {
        PS(i,p7S_P)  = PMX_SP(curr_i,k);
        PI(i,p7S_PI) = lookback[0][k];
        PI(i,p7S_PK) = k;
      }
    }

    /* Record best scoring donor sites*/
    if(SIGNAL(sub_dsq[i-MIN_INTRON_LENG-1], sub_dsq[i-MIN_INTRON_LENG]) == DONOR_GT) {
      t = sub_dsq[i-MIN_INTRON_LENG-3];
      u = sub_dsq[i-MIN_INTRON_LENG-2];
      if(t < 4 && u < 4) {
        for (k = 2; k < M; k++) {
          donor_i = (i-MIN_INTRON_LENG-4) % (MIN_INTRON_LENG+5); 
          TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1)); 
          if(TMP_SC > SSX2(k, p7S_GTAG, t, u)) {
            SSX2(k, p7S_GTAG, t, u) = TMP_SC;
            SIX2(k, p7S_GTAG, t, u) = i-MIN_INTRON_LENG;
          }
        }
      }
      
      if(u < 4) {
        for (k = 2; k < M; k++) {
          donor_i = (i-MIN_INTRON_LENG-3) % (MIN_INTRON_LENG+5);
          TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
          if(TMP_SC > SSX1(k, p7S_GTAG, u)) {
            SSX1(k, p7S_GTAG, u) = TMP_SC;
            SIX1(k, p7S_GTAG, u) = i-MIN_INTRON_LENG; 
          } 
        }
      }
      
      for (k = 2; k < M; k++) {
        donor_i = (i-MIN_INTRON_LENG-2) % (MIN_INTRON_LENG+5);
        TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
        if(TMP_SC > SSX0(k, p7S_GTAG)) {
          SSX0(k, p7S_GTAG) = TMP_SC;
          SIX0(k, p7S_GTAG) = i-MIN_INTRON_LENG;
        }
      } 
    }
    else if(SIGNAL(sub_dsq[i-MIN_INTRON_LENG-1], sub_dsq[i-MIN_INTRON_LENG]) == DONOR_GC) {
      t = sub_dsq[i-MIN_INTRON_LENG-3];
      u = sub_dsq[i-MIN_INTRON_LENG-2];
      if(t < 4 && u < 4) {
        for (k = 2; k < M; k++) {
          donor_i = (i-MIN_INTRON_LENG-4) % (MIN_INTRON_LENG+5);
          TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
          if(TMP_SC > SSX2(k, p7S_GCAG, t, u)) {
            SSX2(k, p7S_GCAG, t, u) = TMP_SC;
            SIX2(k, p7S_GCAG, t, u) = i-MIN_INTRON_LENG;
          }
        }
      }

      if(u < 4) {
        for (k = 2; k < M; k++) {
          donor_i = (i-MIN_INTRON_LENG-3) % (MIN_INTRON_LENG+5);
          TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
          if(TMP_SC > SSX1(k, p7S_GCAG, u)) {
            SSX1(k, p7S_GCAG, u) = TMP_SC;
            SIX1(k, p7S_GCAG, u) = i-MIN_INTRON_LENG;
          }
        }
      }

      for (k = 2; k < M; k++) {
        donor_i = (i-MIN_INTRON_LENG-2) % (MIN_INTRON_LENG+5);
        TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
        if(TMP_SC > SSX0(k, p7S_GCAG)) {
          SSX0(k, p7S_GCAG) = TMP_SC;
          SIX0(k, p7S_GCAG) = i-MIN_INTRON_LENG;
        }
      }    
    }
    else if(SIGNAL(sub_dsq[i-MIN_INTRON_LENG-1], sub_dsq[i-MIN_INTRON_LENG]) == DONOR_AT) {
      t = sub_dsq[i-MIN_INTRON_LENG-3];
      u = sub_dsq[i-MIN_INTRON_LENG-2];
      if(t < 4 && u < 4) {
        for (k = 2; k < M; k++) {
          donor_i = (i-MIN_INTRON_LENG-4) % (MIN_INTRON_LENG+5);
          TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
          if(TMP_SC > SSX2(k, p7S_ATAC, t, u)) {
            SSX2(k, p7S_ATAC, t, u) = TMP_SC;
            SIX2(k, p7S_ATAC, t, u) = i-MIN_INTRON_LENG;
          }
        }
      }

      if(u < 4) {
        for (k = 2; k < M; k++) {
          donor_i = (i-MIN_INTRON_LENG-3) % (MIN_INTRON_LENG+5);
          TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
          if(TMP_SC > SSX1(k, p7S_ATAC, u)) {
            SSX1(k, p7S_ATAC, u) = TMP_SC;
            SIX1(k, p7S_ATAC, u) = i-MIN_INTRON_LENG;
          }
        }
      }

      for (k = 2; k < M; k++) {
        donor_i = (i-MIN_INTRON_LENG-2) % (MIN_INTRON_LENG+5);
        TMP_SC = ESL_MAX(MMX_SP(donor_i,k-1), DMX_SP(donor_i,k-1));
        if(TMP_SC > SSX0(k, p7S_ATAC)) {
          SSX0(k, p7S_ATAC) = TMP_SC;
          SIX0(k, p7S_ATAC) = i-MIN_INTRON_LENG;
        }
      }
    }

    //    p7_gmx_sp_DumpRow(stdout, gx, i, curr_i, 0, gm_fs->M, p7_DEFAULT);
  } // end loop over L
 
  gx->M = gm_fs->M;
  gx->L = L;

  return eslOK;

}

/* MAX lookback is i+MIN_INTRON_LENG+7 */
//Global on right end of model (must enter at M)
////Used to find long introns and determine left end termination on both sequence and model
int
p7_spliceviterbi_leftparser_semiglobal(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx)
{
  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float      **score = pli->sig_idx->score;
  float       *signal_scores = pli->signal_scores;
  float       *parser_scores = pli->sig_idx->parser_scores;
  int        **index = pli->sig_idx->index;
  int        **lookback = pli->sig_idx->lookback;
  int         *parser_index = pli->sig_idx->parser_index;
  int          M    = gm_fs->M;
  int          i,k,s;
  int          v, w, x, y, z;
  int          c3, c2, c1;
  int          nuc1, nuc2;
  int          XGT, XXGT, XXXGT;
  int          XGC, XXGC, XXXGC;
  int          XAT, XXAT, XXXAT;
  int          curr_i, prev_i;
  float        TMP_SC; 


  /*Initialize the signal index */
  for(k = 0; k < M; k++) 
    lookback[0][k] = -1;
  
  for(i = 0; i <= L; i++) { 
    PI(i,p7S_PI) = PI(i,p7S_PK) = PI(i,p7S_MK) = -1;
    PS(i,p7S_P) =  PS(i,p7S_M) = -eslINFINITY;
  }

  for(k = 0; k <  M; k++) {
    for(s = 0; s < SIGNAL_MEM_SIZE; s++) {
      index[k][s] = -1;
      score[k][s] = -eslINFINITY;
    }
  }
//p7_gmx_sp_DumpHeader(stdout, gx, 0, gm_fs->M, p7_DEFAULT);
  /* Initialization of the L to L-1 rows.  */
  for(i = L; i >= L-1; i--) {
    curr_i = i % (MIN_INTRON_LENG+8);
    XMX_SP(i,p7G_J) = XMX_SP(i,p7G_B) = XMX_SP(i,p7G_N) = -eslINFINITY; 
    XMX_SP(i,p7G_C) = 0.; 
    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE];
    for (k = M; k >= 0; k--) {
      MMX_SP(curr_i,k) = IMX_SP(curr_i,k) = PMX_SP(curr_i,k) = DMX_SP(curr_i,k) = -eslINFINITY;
    }
//p7_gmx_sp_DumpRow(stdout, gx, i, curr_i, 0, gm_fs->M, p7_DEFAULT);
  }
   
 
  /* Initialization of the L-2 row because we have no L+1 to transtion from  */
  if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[L-2])) v = sub_dsq[L-2];
  else                                              v = p7P_MAXCODONS;
  if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[L-1])) w = sub_dsq[L-1];
  else                                              w = p7P_MAXCODONS;
  if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[L]))   x = sub_dsq[L];
  else                                              x = p7P_MAXCODONS;

  c3 = p7P_CODON3(v, w, x);
  c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

  curr_i = (L-2) % (MIN_INTRON_LENG+8);

  XMX_SP(L-2,p7G_J) = XMX_SP(L-2,p7G_B) = XMX_SP(L-2,p7G_N) = -eslINFINITY;
  XMX_SP(L-2,p7G_C) = 0.;  
  XMX_SP(L-2,p7G_E) = XMX_SP(L-2,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE];

  MMX_SP(curr_i,M) = gm_fs->xsc[p7P_E][p7P_MOVE] + p7P_MSC_CODON(gm_fs, M, c3);
  IMX_SP(curr_i,M) = DMX_SP(curr_i,M) = PMX_SP(curr_i,M) = -eslINFINITY;
  for (k = M-1; k >= 1; k--) {
    MMX_SP(curr_i,k) = IMX_SP(curr_i,k) = PMX_SP(curr_i,k) = -eslINFINITY; 
    DMX_SP(curr_i,k) = ESL_MAX(MMX_SP(curr_i,k+1) + TSC(p7P_DM,k),
                               DMX_SP(curr_i,k+1) + TSC(p7P_DD,k));
  }
  MMX_SP(curr_i,0) = IMX_SP(curr_i,0) = DMX_SP(curr_i,0) = PMX_SP(curr_i,0) = -eslINFINITY;

//p7_gmx_sp_DumpRow(stdout, gx, L-2, curr_i, 0, gm_fs->M, p7_DEFAULT);
  /* No P state for rows L-3 to L-MIN_INTRON_LENG-6 */
  for (i = L-3; i >= ESL_MAX(L-MIN_INTRON_LENG-6, 1); i--)
  {

    x = w;
    w = v;
  
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) v = sub_dsq[i];
    else                                                v = p7P_MAXCODONS;

    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    curr_i = i % (MIN_INTRON_LENG+8);
    prev_i = (i+3) % (MIN_INTRON_LENG+8);

    XMX_SP(i,p7G_J) = -eslINFINITY;
    XMX_SP(i,p7G_C) = XMX_SP(i+3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP];
    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C)   + gm_fs->xsc[p7P_E][p7P_MOVE];
       
    MMX_SP(curr_i,M) = XMX_SP(i+3,p7G_E) + p7P_MSC_CODON(gm_fs, M, c3);
    IMX_SP(curr_i,M) = DMX_SP(curr_i,M) = PMX_SP(curr_i,M)  = -eslINFINITY;

    XMX_SP(i,p7G_B) = MMX_SP(curr_i,M);

    if(MMX_SP(curr_i,M) > PS(i,p7S_M)) {
      PS(i,p7S_M)  = MMX_SP(curr_i,M);
      PI(i,p7S_MK) = M;
    }

    for (k = M-1; k >= 1; k--) {
      MMX_SP(curr_i,k) = ESL_MAX(MMX_SP(prev_i,k+1) + TSC(p7P_MM,k),
                         ESL_MAX(IMX_SP(prev_i,k)   + TSC(p7P_MI,k),
                                 DMX_SP(prev_i,k+1) + TSC(p7P_MD,k))) + p7P_MSC_CODON(gm_fs, k, c3);

      if(p7P_MSC_CODON(gm_fs, k, c3) == -eslINFINITY)
        IMX_SP(curr_i,k) = -eslINFINITY;
      else
        IMX_SP(curr_i,k) = ESL_MAX(MMX_SP(prev_i,k+1) + TSC(p7P_IM,k),  
                                   IMX_SP(prev_i,k)   + TSC(p7P_II,k));

      DMX_SP(curr_i,k) = ESL_MAX(MMX_SP(curr_i,k+1) + TSC(p7P_DM,k), 
                                 DMX_SP(curr_i,k+1)   + TSC(p7P_DD,k)); 

      PMX_SP(curr_i,k) = -eslINFINITY;

      XMX_SP(i,p7G_B) = ESL_MAX(XMX_SP(i,p7G_B), MMX_SP(curr_i,k));

      if(MMX_SP(curr_i,k) > PS(i,p7S_M)) {
        PS(i,p7S_M)  = MMX_SP(curr_i,k);
        PI(i,p7S_MK) = k;
      }  
      
    }
    MMX_SP(curr_i,0) = IMX_SP(curr_i,0) = DMX_SP(curr_i,0) = PMX_SP(curr_i,0) = -eslINFINITY;  

    XMX_SP(i,p7G_N) = ESL_MAX(XMX_SP(i+3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP],
                              XMX_SP(i,p7G_B)   + gm_fs->xsc[p7P_N][p7P_MOVE]);
//p7_gmx_sp_DumpRow(stdout, gx, i, curr_i, 0, gm_fs->M, p7_DEFAULT);
  } 


  /* Main recursion */
  XXXGT = XXGT = XGT = FALSE;
  XXXGC = XXGC = XGC = FALSE;
  XXXAT = XXAT = XAT = FALSE;
  for (i = L-MIN_INTRON_LENG-7; i >= 1; i--)
  {
    x = w;
    w = v;

    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) v = sub_dsq[i];
    else                                                v = p7P_MAXCODONS;

    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    curr_i = i % (MIN_INTRON_LENG+8);
    prev_i = (i+3) % (MIN_INTRON_LENG+8);

    XMX_SP(i,p7G_J) = -eslINFINITY;
    XMX_SP(i,p7G_C) = XMX_SP(i+3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP];
    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C)   + gm_fs->xsc[p7P_E][p7P_MOVE];

    MMX_SP(curr_i,M) = XMX_SP(i+3,p7G_E) +  p7P_MSC_CODON(gm_fs, M, c3);
    IMX_SP(curr_i,M) = DMX_SP(curr_i,M) = PMX_SP(curr_i,M) = -eslINFINITY;

    XMX_SP(i,p7G_B) = MMX_SP(curr_i,M);

    if(MMX_SP(curr_i,M) > PS(i,p7S_M)) {
      PS(i,p7S_M)  = MMX_SP(curr_i,M);
      PI(i,p7S_MK) = M;
    }

    for (k = M-1; k >= 1; k--) {
      MMX_SP(curr_i,k) = ESL_MAX(MMX_SP(prev_i,k+1) + TSC(p7P_MM,k),
                         ESL_MAX(IMX_SP(prev_i,k)   + TSC(p7P_MI,k),
                         ESL_MAX(DMX_SP(prev_i,k+1) + TSC(p7P_MD,k),
                                 PMX_SP(prev_i,k+1)))) + p7P_MSC_CODON(gm_fs, k, c3);

      if(p7P_MSC_CODON(gm_fs, k, c3) == -eslINFINITY)
        IMX_SP(curr_i,k) = -eslINFINITY;
      else
        IMX_SP(curr_i,k) = ESL_MAX(MMX_SP(prev_i,k+1) + TSC(p7P_IM,k),
                                   IMX_SP(prev_i,k)   + TSC(p7P_II,k));

      DMX_SP(curr_i,k) = ESL_MAX(MMX_SP(curr_i,k+1) + TSC(p7P_DM,k),
                         ESL_MAX(DMX_SP(curr_i,k+1) + TSC(p7P_DD,k),
                                 PMX_SP(curr_i,k+1))); 

      PMX_SP(curr_i,k) = -eslINFINITY;

      XMX_SP(i,p7G_B) = ESL_MAX(XMX_SP(i,p7G_B), MMX_SP(curr_i,k));

      if(MMX_SP(curr_i,k) > PS(i,p7S_M)) {
        PS(i,p7S_M)  = MMX_SP(curr_i,k);
        PI(i,p7S_MK) = k;
      }
    }
    MMX_SP(curr_i,0) = IMX_SP(curr_i,0) = DMX_SP(curr_i,0) = PMX_SP(curr_i,0) = -eslINFINITY;

    XMX_SP(i,p7G_N) = ESL_MAX(XMX_SP(i+3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP],
                              XMX_SP(i,p7G_B)   + gm_fs->xsc[p7P_N][p7P_MOVE]);

    /* Check for donor */
    XXXGT = XXGT;
    XXGT  = XGT;
    if(SIGNAL(w,x) == DONOR_GT) XGT = TRUE;
    else                        XGT = FALSE;
  
    XXXGC = XXGC;
    XXGC  = XGC;
    if(SIGNAL(w,x) == DONOR_GC) XGC = TRUE;
    else                        XGC = FALSE;

    XXXAT = XXAT;
    XXAT  = XAT;
    if(SIGNAL(w,x) == DONOR_AT) XAT = TRUE;
    else                        XAT = FALSE;


    if(XXXGT) {
      for(k = M-1; k >= 2; k--) {
        TMP_SC = SSX0(k, p7S_GTAG) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c3);
        if(TMP_SC > PMX_SP(curr_i,k)){
          PMX_SP(curr_i,k) = TMP_SC;
          lookback[0][k] = SIX0(k, p7S_GTAG);
        }
      }
    } 
    else if(XXXGC) {
      for(k = M-1; k >= 2; k--) {
        TMP_SC = SSX0(k, p7S_GCAG) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c3);
        if(TMP_SC > PMX_SP(curr_i,k)){
          PMX_SP(curr_i,k) = TMP_SC;
          lookback[0][k] = SIX0(k, p7S_GCAG);
        }
      }  
    }
    else if(XXXAT) {
      for(k = M-1; k >= 2; k--) {
        TMP_SC = SSX0(k, p7S_ATAC) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c3);
        if(TMP_SC > PMX_SP(curr_i,k)){
          PMX_SP(curr_i,k) = TMP_SC;
          lookback[0][k] = SIX0(k, p7S_ATAC);
        }
      }
    } 

    if(XXGT) {
      for(k = M-1; k >= 2; k--) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(v, w, nuc1);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);
          TMP_SC = SSX1(k, p7S_GTAG, nuc1) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c2);
          if(TMP_SC > PMX_SP(curr_i,k)) {
            PMX_SP(curr_i,k) = TMP_SC;
            lookback[0][k] = SIX1(k, p7S_GTAG, nuc1);
          }
        }
      }
    }
    else if(XXGC) {
      for(k = M-1; k >= 2; k--) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(v, w, nuc1);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);
          TMP_SC = SSX1(k, p7S_GCAG, nuc1) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c2);
          if(TMP_SC > PMX_SP(curr_i,k)) {
            PMX_SP(curr_i,k) = TMP_SC;
            lookback[0][k] = SIX1(k, p7S_GCAG, nuc1);
          }
        }
      }
    }
    else if(XXAT) {
      for(k = M-1; k >= 2; k--) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(v, w, nuc1);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);
          TMP_SC = SSX1(k, p7S_ATAC, nuc1) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c2);
          if(TMP_SC > PMX_SP(curr_i,k)) {
            PMX_SP(curr_i,k) = TMP_SC;
            lookback[0][k] = SIX1(k, p7S_ATAC, nuc1);
          }
        }
      }
    }

    if(XGT) {
      for(k = M-1; k >= 2; k--) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(v, nuc1, nuc2);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);
            TMP_SC = SSX2(k, p7S_GTAG, nuc1, nuc2) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, k, c1);
            if(TMP_SC > PMX_SP(curr_i,k)) {
              PMX_SP(curr_i,k) = TMP_SC;
              lookback[0][k] = SIX2(k, p7S_GTAG, nuc1, nuc2);
            }
          }
        }
      }
    }
    else if(XGC) {
      for(k = M-1; k >= 2; k--) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(v, nuc1, nuc2);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);
            TMP_SC = SSX2(k, p7S_GCAG, nuc1, nuc2) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, k, c1);
            if(TMP_SC > PMX_SP(curr_i,k)) {
              PMX_SP(curr_i,k) = TMP_SC;
              lookback[0][k] = SIX2(k, p7S_GCAG, nuc1, nuc2);
            }
          }
        }
      }
    }
    else if(XAT) {
      for(k = M-1; k >= 2; k--) {
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(v, nuc1, nuc2);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);
            TMP_SC = SSX2(k, p7S_ATAC, nuc1, nuc2) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, k, c1);
            if(TMP_SC > PMX_SP(curr_i,k)) {
              PMX_SP(curr_i,k) = TMP_SC;
              lookback[0][k] = SIX2(k, p7S_ATAC, nuc1, nuc2);
            }
          }
        }
      }
    }

    for (k = 2; k < M; k++) {
      if(PMX_SP(curr_i,k) > PS(i,p7S_P)) {
        PS(i,p7S_P)  = PMX_SP(curr_i,k);
        PI(i,p7S_PI) = lookback[0][k];
        PI(i,p7S_PK) = k;
      }
    }

    /*Check for acceptor */
    if(SIGNAL(sub_dsq[i+MIN_INTRON_LENG+1],sub_dsq[i+MIN_INTRON_LENG+2]) == ACCEPT_AG) {
      y = sub_dsq[i+MIN_INTRON_LENG+3];
      z = sub_dsq[i+MIN_INTRON_LENG+4];
      if(y < 4 && z < 4) {
        for(k = M-1; k >= 2; k--) {
          TMP_SC = MMX_SP((i+MIN_INTRON_LENG+7) % (MIN_INTRON_LENG+8),k+1) + TSC_P;
          if(TMP_SC > SSX2(k, p7S_GTAG, y, z)) {
            SSX2(k, p7S_GTAG, y, z) = TMP_SC;
            SIX2(k, p7S_GTAG, y, z) = i+MIN_INTRON_LENG+2;
          }
          if(TMP_SC > SSX2(k, p7S_GCAG, y, z)) {
            SSX2(k, p7S_GCAG, y, z) = TMP_SC;
            SIX2(k, p7S_GCAG, y, z) = i+MIN_INTRON_LENG+2;
          }
        }
      }

      if(y < 4) {
        for(k = M-1; k >= 2; k--) {
          TMP_SC = MMX_SP((i+MIN_INTRON_LENG+6) % (MIN_INTRON_LENG+8),k+1) + TSC_P;
          if(TMP_SC > SSX1(k, p7S_GTAG, y)) {
            SSX1(k, p7S_GTAG, y) = TMP_SC;
            SIX1(k, p7S_GTAG, y) = i+MIN_INTRON_LENG+2;
          }
          if(TMP_SC > SSX1(k, p7S_GCAG, y)) {
            SSX1(k, p7S_GCAG, y) = TMP_SC;
            SIX1(k, p7S_GCAG, y) = i+MIN_INTRON_LENG+2;
          }
        }
      }
     
      for(k = M-1; k >= 2; k--) {
        TMP_SC = MMX_SP((i+MIN_INTRON_LENG+5) % (MIN_INTRON_LENG+8),k+1) + TSC_P;
        if(TMP_SC > SSX0(k, p7S_GTAG)) {
          SSX0(k, p7S_GTAG) = TMP_SC;
          SIX0(k, p7S_GTAG) = i+MIN_INTRON_LENG+2;
        }
        if(TMP_SC > SSX0(k, p7S_GCAG)) {
          SSX0(k, p7S_GCAG) = TMP_SC;
          SIX0(k, p7S_GCAG) = i+MIN_INTRON_LENG+2;
        }
      }

    } 
    else if(SIGNAL(sub_dsq[i+MIN_INTRON_LENG+1],sub_dsq[i+MIN_INTRON_LENG+2]) == ACCEPT_AC) {
      y = sub_dsq[i+MIN_INTRON_LENG+3];
      z = sub_dsq[i+MIN_INTRON_LENG+4];
      if(y < 4 && z < 4) {
        for(k = M-1; k >= 2; k--) {
          TMP_SC = MMX_SP((i+MIN_INTRON_LENG+7) % (MIN_INTRON_LENG+8),k+1) + TSC_P;
          if(TMP_SC > SSX2(k, p7S_ATAC, y, z)) {
            SSX2(k, p7S_ATAC, y, z) = TMP_SC; 
            SIX2(k, p7S_ATAC, y, z) = i+MIN_INTRON_LENG+2;
          }
        }
      }

      if(y < 4) {
        for(k = M-1; k >= 2; k--) {
          TMP_SC = MMX_SP((i+MIN_INTRON_LENG+6) % (MIN_INTRON_LENG+8),k+1) + TSC_P;
          if(TMP_SC > SSX1(k, p7S_ATAC, y)) {
            SSX1(k, p7S_ATAC, y) = TMP_SC;
            SIX1(k, p7S_ATAC, y) = i+MIN_INTRON_LENG+2;
          }
        }
      }

      for(k = M-1; k >= 2; k--) {
        TMP_SC = MMX_SP((i+MIN_INTRON_LENG+5) % (MIN_INTRON_LENG+8),k+1) + TSC_P;
        if(TMP_SC > SSX0(k, p7S_ATAC)) {
          SSX0(k, p7S_ATAC) = TMP_SC; 
          SIX0(k, p7S_ATAC) = i+MIN_INTRON_LENG+2;
        }
      }
    }
//p7_gmx_sp_DumpRow(stdout, gx, i, curr_i, 0, gm_fs->M, p7_DEFAULT);
  } //End L loop 

  /* Empty 0 row */
  XMX_SP(0,p7G_C) = XMX_SP(0,p7G_J) = XMX_SP(0,p7G_E) = XMX_SP(0,p7G_B) = XMX_SP(0,p7G_N) = -eslINFINITY;
 
  gx->M = gm_fs->M;
  gx->L = L;
//p7_gmx_sp_Dump(stdout, gx, p7_DEFAULT);
  return eslOK;

}

int
p7_splicevitebi_exon_definition(SPLICE_PIPELINE *pli, SPLICE_PATH *path, P7_GMX *gx, P7_FS_PROFILE *sub_gm, ESL_SQ *sub_seq, int **exon_idx, int *idx_size) 
{

  int i,j,k,s;
  int tmp_j;
  //int step;
  int true_i, true_j;
  int donor;
  int i_size;
  int a_size;
  float tol = 1e-5;
  float *xmx = gx->xmx;
  float *parser_scores = pli->sig_idx->parser_scores;
  int   *parser_index  = pli->sig_idx->parser_index; 
  int   *e_idx = NULL;
  int status;

  ESL_ALLOC(e_idx, sizeof(int) * path->path_len * 4);
  a_size = path->path_len * 2;

  /* First find end of alignment */ 
  for(j = gx->L; j > MIN_INTRON_LENG; j--) {
    if   (XMX_SP(j,p7G_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", j);
    if   (XMX_SP(j, p7G_C) < XMX_SP(j-2, p7G_C) || XMX_SP(j, p7G_C) < XMX_SP(j-1, p7G_C)) continue;  
   //printf("j %d XMX_SP(j, p7G_C) %f XMX_SP(j, p7G_E) %f sub_gm->xsc[p7P_E][p7P_MOVE] %f  \n", j, XMX_SP(j, p7G_C), XMX_SP(j, p7G_E), sub_gm->xsc[p7P_E][p7P_MOVE]); 
    if (esl_FCompare_old(XMX_SP(j, p7G_C), XMX_SP(j, p7G_E) + sub_gm->xsc[p7P_E][p7P_MOVE], tol) == eslOK) break;
  }
  i_size = 0;

  /* Find which step in path we are entering at */  
  if(path->revcomp) true_j = sub_seq->n - j + sub_seq->end;
  else              true_j = sub_seq->start + j - 1;

/*
  for(s = path->path_len-1; s >= 0; s--) {
    if(path->revcomp){
      if(path->iali[s] < true_j) {
        e_idx[(i_size+1)*2]   = sub_seq->n - path->iali[s] + sub_seq->end;
        e_idx[(i_size+1)*2+1] = sub_seq->n - path->jali[s] + sub_seq->end;
        i_size++;

        if((i_size+2) * 2 > a_size) {
          ESL_REALLOC(e_idx, sizeof(int) * (i_size+2) * 4);
          a_size = (i_size+2) * 2;
        }
      }
      //Extend the last exon included in parser tarce 
      else if(path->jali[s] < true_j) {
        true_j = path->jali[s];
        j = sub_seq->n - true_j + sub_seq->end;
      }
    }
    else {
      //Add an additional exon after parser exits model 

      if(path->iali[s] > true_j) {
        e_idx[i_size*2]   = path->iali[s] - sub_seq->start + 1;
        e_idx[i_size*2+1] = path->jali[s] - sub_seq->start + 1;
        i_size++;
        
        if((i_size+2) * 2 > a_size) {
          ESL_REALLOC(e_idx, sizeof(int) * (i_size+2) * 4);
          a_size = (i_size+2) * 2;
        }        
      }
      //Extend the last exon included in parser tarce 
      else if(path->jali[s] > true_j) {
        true_j = path->jali[s];
        j = true_j - sub_seq->start + 1;
      }
    }
  }
*/
  //printf("gx->L %d\n", gx->L);
  //printf("true_j %d\n", true_j);
  //printf("j %d\n", j);
  e_idx[i_size*2 + 1]   = j;
    
  i = j;
  for( i = j; i > MIN_INTRON_LENG; i--) {

    /*If the max P if greater than the max M */
    if(PS(i,p7S_P) > PS(i,p7S_M)) {
      k = PI(i,p7S_PK);
      donor = PI(i,p7S_PI);
          
      /* If the donor site for the max P is more than MAX_INTRON_INCL away */
      if(i - donor > MAX_INTRON_INCL) {

        /* Check if maxP at j is greater than all max M up to donor*/
        j = i;
        while(j > donor && PS(i,p7S_P) > PS(j,p7S_M)) j--;

        if(j == donor) {
/*
          closest_s = path->path_len;
          closest_j = MAX_INTRON_INCL;
          for(s = 0; s < path->path_len; s++) {
            if(path->revcomp) { 
              iali = sub_seq->n + sub_seq->end - path->iali[s];
              jali = sub_seq->n + sub_seq->end - path->jali[s];
            }
            else {
              iali = path->iali[s] - sub_seq->start + 1;
              jali = path->jali[s] - sub_seq->start + 1; 
            }   
          
            if(abs(jali - donor) < closest_j) {
              closest_j = abs(jali - donor);
              closest_s = s;
            }    
          }
         
          if(s < path->path_len) {
            closest_k = closest_j * 2 / 3;
              
          }

*/

           //printf("acceptor %d donor %d\n", i, donor);
          /* Grow exon index if needed */
          if((i_size+2) * 2 > a_size) {
            ESL_REALLOC(e_idx, sizeof(int) * (i_size+2) * 4);
            a_size = (i_size+2) * 2;
          }

         /* start of previous exon region (from last to first) */
          e_idx[i_size*2] = i;
          /* start of next exon region (from last to first) */
          e_idx[(i_size+1)*2+1] = j; 
          
        //printf("i %d j %d k %d\n", i, j, k);  
          if(path->revcomp) {
            true_i = sub_seq->n - i + sub_seq->end;
            true_j = sub_seq->n - j + sub_seq->end;
            //printf("true_i %d true_j %d\n", true_i, true_j);
            /* Extend into intro region from upstream*/
            for(s = path->path_len-1; s >= 0; s--) {
              /* If there is a step in the path that ends after the previous exon region
               * starts but extends further downstream, extend the exon region to match */
               if(path->jali[s] <= true_i && path->iali[s] > true_i) {
                 e_idx[i_size*2] = sub_seq->n - path->iali[s] + sub_seq->end;
                 true_i = path->iali[s];
               }
            }
    
            //printf("true_i %d true_j %d\n", true_i, true_j);
            /* Extend into intro region from downstream*/
            tmp_j = j;
            //printf("tmp_j %d \n", tmp_j);
            for(s = 0; s < path->path_len; s++) {
               /* If there is a step in the path that starts before the next exon region
               * ends but extends further upstream, extend the exon region to match */
              if(path->iali[s] >= true_j && path->jali[s] < true_j) {
                e_idx[(i_size+1)*2+1] = sub_seq->n - path->jali[s] + sub_seq->end;
                tmp_j  = e_idx[(i_size+1)*2+1];
                true_j = path->jali[s];
              }
            }
            
            //printf("true_i %d true_j %d\n", true_i, true_j);
            /* Add new exon region */
            for(s = path->path_len-1; s >= 0; s--) {
              /* If there is a step in the path that lies between the pervious and the
               * next exon regions add a new exon region */
              if(path->jali[s] > true_i && path->iali[s] < true_j) {
                e_idx[(i_size+1)*2]   = sub_seq->n - path->iali[s] + sub_seq->end;
                e_idx[(i_size+1)*2+1] = sub_seq->n - path->jali[s] + sub_seq->end;
//printf("added e_idx[(i_size+1)*2] %d e_idx[(i_size+1)*2+1] %d\n", e_idx[(i_size+1)*2], e_idx[(i_size+1)*2+1]);
                i_size++;

                if((i_size+2) * 2 > a_size) {
                  ESL_REALLOC(e_idx, sizeof(int) * (i_size+2) * 4);
                  a_size = (i_size+2) * 2;
                }
                e_idx[(i_size+1)*2+1] = tmp_j;     
              }
            } 
             
            //printf("true_i %d true_j %d\n", true_i, true_j); 
          }
          else {
            true_i = sub_seq->start + i - 1;
            true_j = sub_seq->start + j - 1;

            /* Extend into intro region from upstream*/
            for(s = path->path_len-1; s >= 0; s--) {
              /* If there is a step in the path that ends after the previous exon region
               * starts but extends further downstream, extend the exon region to match */
              if(path->jali[s] >= true_i && path->iali[s] < true_i) {
                e_idx[i_size*2] = path->iali[s] - sub_seq->start + 1;
                true_i = path->iali[s];
              }
            }

            //printf("i %d j %d\n", e_idx[i_size*2], e_idx[(i_size+1)*2+1]);
            /* Extend into intron region from downstream*/
            tmp_j = j;
            for(s = 0; s < path->path_len; s++) {
              /* If there is a step in the path that starts before the next exon region
               * ends but extends further upstream, extend the exon region to match */
              if(path->iali[s] <= true_j && path->jali[s] > true_j) {
                e_idx[(i_size+1)*2+1] = path->jali[s] - sub_seq->start + 1;
                tmp_j                 = path->jali[s] - sub_seq->start + 1;
                true_j = path->jali[s];
              }
            }
           
            //printf("i %d j %d\n", e_idx[i_size*2], e_idx[(i_size+1)*2+1]);
            /* Add new exon region */
            for(s = path->path_len-1; s >= 0; s--) {
              /* If there is a step in the path that lies between the pervious and the
               * next exon regions add a new exon region */
              if(path->jali[s] < true_i && path->iali[s] > true_j) {
                e_idx[(i_size+1)*2]   = path->iali[s] - sub_seq->start + 1;
                e_idx[(i_size+1)*2+1] = path->jali[s] - sub_seq->start + 1;
                i_size++;
            
                if((i_size+2) * 2 > a_size) {
                  ESL_REALLOC(e_idx, sizeof(int) * (i_size+2) * 4);
                  a_size = (i_size+2) * 2;
                }
                e_idx[(i_size+1)*2+1] = tmp_j;
              }
            } 
           //printf("i %d j %d\n", e_idx[i_size*2], e_idx[(i_size+1)*2+1]);
          }
                  
          i_size++;
          i = j;
        }
      }
    }
  } 
  e_idx[i_size*2] = 1;
  i_size++;
 
  /* ensure exon regions are at least one nucletoide long and do not overlap */
  //printf("e_idx[0] %d e_idx[1] %d\n", e_idx[0], e_idx[1]);
  if(e_idx[0] > e_idx[1]) e_idx[0] = e_idx[1];
  for(i = 1; i < i_size; i++) {
  //  printf("e_idx[%d] %d e_idx[%d] %d\n", i*2, e_idx[i*2], i*2+1, e_idx[i*2+1]);
    if(e_idx[i*2+1] >= e_idx[(i-1)*2]) e_idx[i*2+1] = e_idx[(i-1)*2] - 1; 
    if(e_idx[i*2] > e_idx[i*2+1])      e_idx[i*2]   = e_idx[i*2+1];
  }

  *exon_idx = e_idx;
  *idx_size   = i_size;

  return eslOK;

  ERROR:
    if(e_idx != NULL) free(e_idx);
    return status;
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
      if   (XMX_SP(i,p7G_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
  
      if      (XMX_SP(i, p7G_C) < XMX_SP(i-2, p7G_C) || XMX_SP(i, p7G_C) < XMX_SP(i-1, p7G_C))                        scur = p7T_C; 
      if      (esl_FCompare_old(XMX_SP(i, p7G_C), XMX_SP(i-3, p7G_C) + sub_gm->xsc[p7P_C][p7P_LOOP], tol) == eslOK) scur = p7T_C;
      else if (esl_FCompare_old(XMX_SP(i, p7G_C), XMX_SP(i,   p7G_E) + sub_gm->xsc[p7P_E][p7P_MOVE], tol) == eslOK) scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);

      break;

    case p7T_E:     /* E connects from any M state. k set here */
      if (XMX_SP(i, p7G_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);

      if      (esl_FCompare_old(XMX_SP(i, p7G_E), MMX_SP(i,M), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(XMX_SP(i, p7G_E), DMX_SP(i,M), tol) == eslOK) scur = p7T_D;
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

      if(k == 1) {
        if (esl_FCompare_old(MMX_SP(i,k), XMX_SP(i-3, p7G_B) + emit, tol) == eslOK) scur = p7T_B;
        else ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
      } 
      else {
        
        if      (esl_FCompare_old(MMX_SP(i,k), MMX_SP(i-3, k-1) + TSC(p7P_MM, k-1) + emit, tol) == eslOK) scur = p7T_M;
        else if (esl_FCompare_old(MMX_SP(i,k), IMX_SP(i-3, k-1) + TSC(p7P_IM, k-1) + emit, tol) == eslOK) scur = p7T_I;
        else if (esl_FCompare_old(MMX_SP(i,k), DMX_SP(i-3, k-1) + TSC(p7P_DM, k-1) + emit, tol) == eslOK) scur = p7T_D;      
        else if (esl_FCompare_old(MMX_SP(i,k), PMX_SP(i-3, k-1) + TSC_P            + emit, tol) == eslOK) scur = p7T_P;
        else ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
         
      }

      k--; i-=3;
      break;

    case p7T_D:         /* D connects from M,D at i,k-1 */
      if (DMX_SP(i, k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);
      
      if      (esl_FCompare_old(DMX_SP(i,k), MMX_SP(i, k-1) + TSC(p7P_MD, k-1), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(DMX_SP(i,k), DMX_SP(i, k-1) + TSC(p7P_DD, k-1), tol) == eslOK) scur = p7T_D;
      else if (esl_FCompare_old(DMX_SP(i,k), XMX_SP(i, p7G_B)  + TSC(p7P_BM, k-1), tol) == eslOK) scur = p7T_B;
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
      //printf("i %d\n", i);
      break;

    case p7T_N:         /* N connects from S, N */
      if (XMX_SP(i, p7G_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = ( (i == 0) ? p7T_S : p7T_N);
      break;

    case p7T_B:         /* B connects from N, J */
      if (XMX_SP(i,p7G_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

      if      (esl_FCompare_old(XMX_SP(i,p7G_B), XMX_SP(i, p7G_N) + sub_gm->xsc[p7P_N][p7P_MOVE], tol) == eslOK) scur = p7T_N;
      else  ESL_EXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over statetype[tpos-1] */


    /* Need to find which type of splice codon we have */
    if (scur == p7T_P) {
      if (PMX_SP(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible P reached at k=%d,i=%d", k,i);

      donor_idx = pli->sig_idx->lookback[i][k];
      //printf("donor_idx %d\n", donor_idx);  
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
        //printf("snxt %d c %d\n", snxt, c);
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





