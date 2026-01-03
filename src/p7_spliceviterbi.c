/* Spliced Viterbi algorithm and trace back */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_splice.h"

#define TSC_P -10.0


int
p7_spliceviterbi_translated_semiglobal_extendup(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, const P7_FS_PROFILE *gm_fs, P7_GMX *gx, int i_start, int i_end, int k_start, int k_end)
{
  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float      **score = pli->sig_idx->score;
  float       *signal_scores = pli->signal_scores;
  int        **index = pli->sig_idx->index;
  int        **lookback = pli->sig_idx->lookback;
  int          L = (i_end - i_start + 1);
  int          M = (k_end - k_start + 1);
  int          i,k,s;
  int          c3, c2, c1;
  int          t, u, v, w, x;
  int          AGXXX, ACXXX;
  int          AGXX, ACXX;
  int          AGX,  ACX;
  int          GT, GC, AT;
  int          nuc1, nuc2;
  int          loop_end;
  int          sub_i,sub_k;
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
    sub_i = i_start + i - 1;
    /* if new nucleotide is not A,C,G, or T set it to placeholder value */
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i])) x = sub_dsq[sub_i];
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
  loop_end = ESL_MIN(L, MIN_INTRON_LENG+4);
  for (i = 3; i <= loop_end; i++)
  {
    v = w;
    w = x;
    sub_i = i_start + i - 1;
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i])) x = sub_dsq[sub_i];
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

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    XMX_SP(i,p7G_E) = -eslINFINITY;
    for (k = 1; k < M; k++) {
      sub_k = k_start + k -1;

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,sub_k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,sub_k-1),
                    ESL_MAX(DMX_SP(i-3,k-1) + TSC(p7P_DM,sub_k-1),
                            XMX_SP(i-3,p7G_B)))) + p7P_MSC_CODON(gm_fs, sub_k, c3);

      if(p7P_MSC_CODON(gm_fs, sub_k, c3) == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;
      else                                                IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,sub_k),
                                                                                IMX_SP(i-3,k) + TSC(p7P_II,sub_k));

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,sub_k-1));

      PMX_SP(i,k) = -eslINFINITY;
    }

    sub_k = k_start + M -1;
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,sub_k-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,sub_k-1),
                          DMX_SP(i-3,M-1) + TSC(p7P_DM,sub_k-1)))+ p7P_MSC_CODON(gm_fs, sub_k, c3);

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,sub_k-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,sub_k-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;

  }

  GT = GC = AT = FALSE;
  /* Main DP recursion */
  for (i = MIN_INTRON_LENG+5; i <= L; i++)
  {
    /* get nucleotides and codon */
    v = w;
    w = x;

    sub_i = i_start + i - 1;
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i])) x = sub_dsq[sub_i];
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
    t = sub_dsq[sub_i-MIN_INTRON_LENG-3];
    u = sub_dsq[sub_i-MIN_INTRON_LENG-2];
    if(SIGNAL(sub_dsq[sub_i-MIN_INTRON_LENG-1], sub_dsq[sub_i-MIN_INTRON_LENG]) == DONOR_GT) GT = TRUE;
    else                         GT = FALSE;

    if(SIGNAL(sub_dsq[sub_i-MIN_INTRON_LENG-1], sub_dsq[sub_i-MIN_INTRON_LENG]) == DONOR_GC) GC = TRUE;
    else                         GC = FALSE;

    if(SIGNAL(sub_dsq[sub_i-MIN_INTRON_LENG-1], sub_dsq[sub_i-MIN_INTRON_LENG]) == DONOR_AT) AT = TRUE;
    else                         AT = FALSE;

    XMX_SP(i,p7G_N) = XMX_SP(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];
    XMX_SP(i,p7G_B) = XMX_SP(i,p7G_N)   + gm_fs->xsc[p7P_N][p7P_MOVE];

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    for (k = 1; k < M; k++) {

      sub_k = k_start + k -1;

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,sub_k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,sub_k-1),
                    ESL_MAX(DMX_SP(i-3,k-1) + TSC(p7P_DM,sub_k-1),
                    ESL_MAX(XMX_SP(i-3,p7G_B),
                            PMX_SP(i-3,k-1) + TSC_P)))) + p7P_MSC_CODON(gm_fs, sub_k, c3);
      
      if(p7P_MSC_CODON(gm_fs, sub_k, c3) == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;
      else                                                IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,sub_k),
                                                                                IMX_SP(i-3,k) + TSC(p7P_II,sub_k));

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,sub_k-1));


      PMX_SP(i,k) = -eslINFINITY;
    }

    sub_k = k_start + M -1;
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,sub_k-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,sub_k-1),
                  ESL_MAX(DMX_SP(i-3,M-1) + TSC(p7P_DM,sub_k-1),
                          PMX_SP(i-3,M-1) + TSC_P)))         + p7P_MSC_CODON(gm_fs, sub_k, c3);

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,sub_k-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,sub_k-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;

    if(AGXXX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        TMP_SC = SSX0(k, p7S_GTAG) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, sub_k, c3);
        if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX0(k, p7S_GTAG);
        }

        TMP_SC = SSX0(k, p7S_GCAG) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, sub_k, c3);
        if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX0(k, p7S_GCAG);
        }
      }
    }
    else if(ACXXX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        TMP_SC = SSX0(k, p7S_ATAC) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, sub_k, c3);
        if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX0(k, p7S_ATAC);
        }
      }
    }

    if(AGXX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(nuc1, w, x);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);

          TMP_SC = SSX1(k, p7S_GTAG, nuc1) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, sub_k, c2);

          if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX1(k, p7S_GTAG, nuc1);
          }

          TMP_SC = SSX1(k, p7S_GCAG, nuc1) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, sub_k, c2);
          if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX1(k, p7S_GCAG, nuc1);
          }
        }
      }
    }
    else if(ACXX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(nuc1, w, x);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);

          TMP_SC = SSX1(k, p7S_ATAC, nuc1) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, sub_k, c2);
          if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX1(k, p7S_ATAC, nuc1);
          }
        }
      }
    }
    if(AGX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(nuc1, nuc2, x);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);

            TMP_SC = SSX2(k, p7S_GTAG, nuc1, nuc2) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, sub_k, c1);
            if(TMP_SC > PMX_SP(i,k)) {
              PMX_SP(i,k) = TMP_SC;
              lookback[i][k] = SIX2(k, p7S_GTAG, nuc1, nuc2);
            }

            TMP_SC = SSX2(k, p7S_GCAG, nuc1, nuc2) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, sub_k, c1);
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
        sub_k = k_start + k -1;
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(nuc1, nuc2, x);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);

            TMP_SC = SSX2(k, p7S_ATAC, nuc1, nuc2) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, sub_k, c1);
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

  /*exit from last i and k */
  XMX_SP(L,p7G_E) = ESL_MAX(MMX_SP(L,M), DMX_SP(L,M));
  XMX_SP(L,p7G_C) = XMX_SP(L,p7G_E) + gm_fs->xsc[p7P_E][p7P_MOVE];

  gx->M = M;
  gx->L = L;

  return eslOK;

}


int
p7_spliceviterbi_translated_semiglobal(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, const P7_FS_PROFILE *gm_fs, P7_GMX *gx, int i_start, int i_end, int k_start, int k_end)
{
  
  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float      **score = pli->sig_idx->score;
  float       *signal_scores = pli->signal_scores;
  int        **index = pli->sig_idx->index;
  int        **lookback = pli->sig_idx->lookback;
  int          L = (i_end - i_start + 1);
  int          M = (k_end - k_start + 1);
  int          i,k,s;
  int          c3, c2, c1;
  int          t, u, v, w, x;
  int          AGXXX, ACXXX;
  int          AGXX, ACXX;
  int          AGX,  ACX;
  int          GT, GC, AT;
  int          nuc1, nuc2;
  int          loop_end;
  int          sub_i,sub_k;
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
    sub_i = i_start + i - 1;
    /* if new nucleotide is not A,C,G, or T set it to placeholder value */
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i])) x = sub_dsq[sub_i];
    else                                                x = p7P_MAXCODONS;

    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    for (k = 0; k <= M; k++)
      MMX_SP(i,k) = DMX_SP(i,k) = IMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;
  }

  /*Special cases for the rows 3-MIN_INTRON_LENG+4*/
  AGXXX = ACXXX = AGXX = ACXX = AGX = ACX  = FALSE;
  loop_end = ESL_MIN(L, MIN_INTRON_LENG+4);
  for (i = 3; i <= loop_end; i++)
  {
    v = w;
    w = x;
    sub_i = i_start + i - 1;
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i])) x = sub_dsq[sub_i];
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

    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    sub_k = k_start;

    MMX_SP(i,1) = XMX_SP(i-3,p7G_B) + p7P_MSC_CODON(gm_fs, sub_k, c3);

    if(p7P_MSC_CODON(gm_fs, sub_k, c3) == -eslINFINITY)
      IMX_SP(i,1) = -eslINFINITY;
    else
      IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,sub_k),
                            IMX_SP(i-3,1) + TSC(p7P_II,sub_k));

    DMX_SP(i,1) = PMX_SP(i,1) = -eslINFINITY;

    for (k = 2; k < M; k++) {
      sub_k = k_start + k -1;

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,sub_k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,sub_k-1),
                            DMX_SP(i-3,k-1) + TSC(p7P_DM,sub_k-1))) + p7P_MSC_CODON(gm_fs, sub_k, c3);
       
      if(p7P_MSC_CODON(gm_fs, sub_k, c3) == -eslINFINITY) 
        IMX_SP(i,k) = -eslINFINITY;
      else  
        IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,sub_k),
                              IMX_SP(i-3,k) + TSC(p7P_II,sub_k));

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,sub_k-1));

      PMX_SP(i,k) = -eslINFINITY;
    }

    sub_k = k_start + M -1;
    
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,sub_k-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,sub_k-1),
                          DMX_SP(i-3,M-1) + TSC(p7P_DM,sub_k-1)))+ p7P_MSC_CODON(gm_fs, sub_k, c3);

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,sub_k-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,sub_k-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;

  }

  GT = GC = AT = FALSE;
  /* Main DP recursion */
  for (i = MIN_INTRON_LENG+5; i <= L; i++)
  {
    /* get nucleotides and codon */
    v = w;
    w = x;

    sub_i = i_start + i - 1;
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i])) x = sub_dsq[sub_i];
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
    t = sub_dsq[sub_i-MIN_INTRON_LENG-3];
    u = sub_dsq[sub_i-MIN_INTRON_LENG-2];
    if(SIGNAL(sub_dsq[sub_i-MIN_INTRON_LENG-1], sub_dsq[sub_i-MIN_INTRON_LENG]) == DONOR_GT) GT = TRUE;
    else                         GT = FALSE;

    if(SIGNAL(sub_dsq[sub_i-MIN_INTRON_LENG-1], sub_dsq[sub_i-MIN_INTRON_LENG]) == DONOR_GC) GC = TRUE;
    else                         GC = FALSE;

    if(SIGNAL(sub_dsq[sub_i-MIN_INTRON_LENG-1], sub_dsq[sub_i-MIN_INTRON_LENG]) == DONOR_AT) AT = TRUE;
    else                         AT = FALSE;

    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    for (k = 1; k < M; k++) {

      sub_k = k_start + k -1;

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,sub_k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,sub_k-1),
                    ESL_MAX(DMX_SP(i-3,k-1) + TSC(p7P_DM,sub_k-1),
                            PMX_SP(i-3,k-1) + TSC_P))) + p7P_MSC_CODON(gm_fs, sub_k, c3);
      
      if(p7P_MSC_CODON(gm_fs, sub_k, c3) == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;
      else                                                IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,sub_k),
                                                                                IMX_SP(i-3,k) + TSC(p7P_II,sub_k));

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,sub_k-1));


      PMX_SP(i,k) = -eslINFINITY;
    }

    sub_k = k_start + M -1;
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,sub_k-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,sub_k-1),
                  ESL_MAX(DMX_SP(i-3,M-1) + TSC(p7P_DM,sub_k-1),
                          PMX_SP(i-3,M-1) + TSC_P)))         + p7P_MSC_CODON(gm_fs, sub_k, c3);
    
    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,sub_k-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,sub_k-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;

    if(AGXXX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        TMP_SC = SSX0(k, p7S_GTAG) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, sub_k, c3);
        if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX0(k, p7S_GTAG);
        }

        TMP_SC = SSX0(k, p7S_GCAG) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, sub_k, c3);
        if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX0(k, p7S_GCAG);
        }
      }
    }
    else if(ACXXX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        TMP_SC = SSX0(k, p7S_ATAC) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, sub_k, c3);
        if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX0(k, p7S_ATAC);
        }
      }
    }

    if(AGXX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(nuc1, w, x);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);
          TMP_SC = SSX1(k, p7S_GTAG, nuc1) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, sub_k, c2);

          if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX1(k, p7S_GTAG, nuc1);
          }

          TMP_SC = SSX1(k, p7S_GCAG, nuc1) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, sub_k, c2);
          if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX1(k, p7S_GCAG, nuc1);
          }
        }
      }
    }
    else if(ACXX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(nuc1, w, x);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);

          TMP_SC = SSX1(k, p7S_ATAC, nuc1) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, sub_k, c2);
          if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX1(k, p7S_ATAC, nuc1);
          }
        }
      }
    }
    if(AGX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(nuc1, nuc2, x);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);

            TMP_SC = SSX2(k, p7S_GTAG, nuc1, nuc2) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, sub_k, c1);
            if(TMP_SC > PMX_SP(i,k)) {
              PMX_SP(i,k) = TMP_SC;
              lookback[i][k] = SIX2(k, p7S_GTAG, nuc1, nuc2);
            }

            TMP_SC = SSX2(k, p7S_GCAG, nuc1, nuc2) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, sub_k, c1);
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
        sub_k = k_start + k -1;
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(nuc1, nuc2, x);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);

            TMP_SC = SSX2(k, p7S_ATAC, nuc1, nuc2) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, sub_k, c1);
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

  /*exit from last i and k */
  XMX_SP(L,p7G_E) = ESL_MAX(MMX_SP(L,M), DMX_SP(L,M));
  XMX_SP(L,p7G_C) = XMX_SP(L,p7G_E) + gm_fs->xsc[p7P_E][p7P_MOVE];

  gx->M = M;
  gx->L = L;
//p7_gmx_sp_Dump(stdout, gx, p7_DEFAULT);
  return eslOK;
}


int
p7_spliceviterbi_translated_semiglobal_extenddown(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, const P7_FS_PROFILE *gm_fs, P7_GMX *gx, int i_start, int i_end, int k_start, int k_end)
{
  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float      **score = pli->sig_idx->score;
  float       *signal_scores = pli->signal_scores;
  int        **index = pli->sig_idx->index;
  int        **lookback = pli->sig_idx->lookback;
  int          L = (i_end - i_start + 1);
  int          M = (k_end - k_start + 1);
  int          i,k,s;
  int          c3, c2, c1;
  int          t, u, v, w, x;
  int          AGXXX, ACXXX;
  int          AGXX, ACXX;
  int          AGX,  ACX;
  int          GT, GC, AT;
  int          nuc1, nuc2;
  int          loop_end;
  int          sub_i,sub_k;
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
    sub_i = i_start + i - 1;
    /* if new nucleotide is not A,C,G, or T set it to placeholder value */
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i])) x = sub_dsq[sub_i];
    else                                                x = p7P_MAXCODONS;

    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    for (k = 0; k <= M; k++)
      MMX_SP(i,k) = DMX_SP(i,k) = IMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;
  }

  /*Special cases for the rows 3-MIN_INTRON_LENG+4*/
  AGXXX = ACXXX = AGXX = ACXX = AGX = ACX  = FALSE;
  loop_end = ESL_MIN(L, MIN_INTRON_LENG+4);
  for (i = 3; i <= loop_end; i++)
  {
    v = w;
    w = x;
    sub_i = i_start + i - 1;
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i])) x = sub_dsq[sub_i];
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

    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    XMX_SP(i,p7G_E) = -eslINFINITY;
 
    sub_k = k_start;

    MMX_SP(i,1) = XMX_SP(i-3,p7G_B) + p7P_MSC_CODON(gm_fs, sub_k, c3); 
    
    if(p7P_MSC_CODON(gm_fs, sub_k, c3) == -eslINFINITY) 
      IMX_SP(i,1) = -eslINFINITY;
    else                                                
      IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,sub_k),
                            IMX_SP(i-3,1) + TSC(p7P_II,sub_k));

    DMX_SP(i,1) = PMX_SP(i,1) = -eslINFINITY; 

    XMX_SP(i,p7G_E) = ESL_MAX(XMX_SP(i,p7G_E), MMX_SP(i,1));
    
    for (k = 2; k < M; k++) {
    
      sub_k = k_start + k -1;  

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,sub_k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,sub_k-1),
                            DMX_SP(i-3,k-1) + TSC(p7P_DM,sub_k-1))) + p7P_MSC_CODON(gm_fs, sub_k, c3);
      if(p7P_MSC_CODON(gm_fs, sub_k, c3) == -eslINFINITY) 
        IMX_SP(i,k) = -eslINFINITY;
      else                         
        IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,sub_k),
                              IMX_SP(i-3,k) + TSC(p7P_II,sub_k));

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,sub_k-1));

      PMX_SP(i,k) = -eslINFINITY;

      XMX_SP(i,p7G_E) = ESL_MAX(XMX_SP(i,p7G_E), 
                        ESL_MAX(MMX_SP(i,k), DMX_SP(i,k)));
    }

    sub_k = k_start + M -1;
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,sub_k-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,sub_k-1),
                          DMX_SP(i-3,M-1) + TSC(p7P_DM,sub_k-1)))+ p7P_MSC_CODON(gm_fs, sub_k, c3);

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,sub_k-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,sub_k-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = ESL_MAX(XMX_SP(i,p7G_E),
                      ESL_MAX(MMX_SP(i,M), DMX_SP(i,M)));

    XMX_SP(i,p7G_C) = ESL_MAX(XMX_SP(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                              XMX_SP(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]); 

    XMX_SP(i,p7G_J) = -eslINFINITY;

  }

  GT = GC = AT = FALSE;
  /* Main DP recursion */
  for (i = MIN_INTRON_LENG+5; i <= L; i++)
  {
    /* get nucleotides and codon */
    v = w;
    w = x;

    sub_i = i_start + i - 1;
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i])) x = sub_dsq[sub_i];
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
    t = sub_dsq[sub_i-MIN_INTRON_LENG-3];
    u = sub_dsq[sub_i-MIN_INTRON_LENG-2];
    if(SIGNAL(sub_dsq[sub_i-MIN_INTRON_LENG-1], sub_dsq[sub_i-MIN_INTRON_LENG]) == DONOR_GT) GT = TRUE;
    else                         GT = FALSE;

    if(SIGNAL(sub_dsq[sub_i-MIN_INTRON_LENG-1], sub_dsq[sub_i-MIN_INTRON_LENG]) == DONOR_GC) GC = TRUE;
    else                         GC = FALSE;

    if(SIGNAL(sub_dsq[sub_i-MIN_INTRON_LENG-1], sub_dsq[sub_i-MIN_INTRON_LENG]) == DONOR_AT) AT = TRUE;
    else                         AT = FALSE;

    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    XMX_SP(i,p7G_E) = -eslINFINITY;

    for (k = 1; k < M; k++) {

      sub_k = k_start + k -1;

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,sub_k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,sub_k-1),
                    ESL_MAX(DMX_SP(i-3,k-1) + TSC(p7P_DM,sub_k-1),
                            PMX_SP(i-3,k-1) + TSC_P))) + p7P_MSC_CODON(gm_fs, sub_k, c3);
       
      if(p7P_MSC_CODON(gm_fs, sub_k, c3) == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;
      else                                                IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,sub_k),
                                                                                IMX_SP(i-3,k) + TSC(p7P_II,sub_k));

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,sub_k-1));

      PMX_SP(i,k) = -eslINFINITY;

      XMX_SP(i,p7G_E) = ESL_MAX(XMX_SP(i,p7G_E),
                        ESL_MAX(MMX_SP(i,k), DMX_SP(i,k)));
    }

    sub_k = k_start + M -1;
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,sub_k-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,sub_k-1),
                  ESL_MAX(DMX_SP(i-3,M-1) + TSC(p7P_DM,sub_k-1),
                          PMX_SP(i-3,M-1) + TSC_P)))         + p7P_MSC_CODON(gm_fs, sub_k, c3);

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,sub_k-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,sub_k-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = ESL_MAX(XMX_SP(i,p7G_E),
                      ESL_MAX(MMX_SP(i,M), DMX_SP(i,M)));

    XMX_SP(i,p7G_C) = ESL_MAX(XMX_SP(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                              XMX_SP(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);

    XMX_SP(i,p7G_J) = -eslINFINITY;

    if(AGXXX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        TMP_SC = SSX0(k, p7S_GTAG) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, sub_k, c3);
        if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX0(k, p7S_GTAG);
        }

        TMP_SC = SSX0(k, p7S_GCAG) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, sub_k, c3);
        if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX0(k, p7S_GCAG);
        }
      }
    }
    else if(ACXXX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        TMP_SC = SSX0(k, p7S_ATAC) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, sub_k, c3);
        if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX0(k, p7S_ATAC);
        }
      }
    }

    if(AGXX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(nuc1, w, x);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);

          TMP_SC = SSX1(k, p7S_GTAG, nuc1) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, sub_k, c2);

          if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX1(k, p7S_GTAG, nuc1);
          }

          TMP_SC = SSX1(k, p7S_GCAG, nuc1) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, sub_k, c2);
          if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX1(k, p7S_GCAG, nuc1);
          }
        }
      }
    }
    else if(ACXX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          c2 = p7P_CODON3(nuc1, w, x);
          c2 = p7P_MINIDX(c2, p7P_DEGEN_C);

          TMP_SC = SSX1(k, p7S_ATAC, nuc1) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, sub_k, c2);
          if(TMP_SC > PMX_SP(i,k)) {
            PMX_SP(i,k) = TMP_SC;
            lookback[i][k] = SIX1(k, p7S_ATAC, nuc1);
          }
        }
      }
    }
    if(AGX) {
      for (k = 2; k < M; k++) {
        sub_k = k_start + k -1;
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(nuc1, nuc2, x);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);

            TMP_SC = SSX2(k, p7S_GTAG, nuc1, nuc2) + signal_scores[p7S_GTAG] + p7P_MSC_CODON(gm_fs, sub_k, c1);
            if(TMP_SC > PMX_SP(i,k)) {
              PMX_SP(i,k) = TMP_SC;
              lookback[i][k] = SIX2(k, p7S_GTAG, nuc1, nuc2);
            }

            TMP_SC = SSX2(k, p7S_GCAG, nuc1, nuc2) + signal_scores[p7S_GCAG] + p7P_MSC_CODON(gm_fs, sub_k, c1);
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
        sub_k = k_start + k -1;
        for(nuc1 = 0; nuc1 < 4; nuc1++) {
          for(nuc2 = 0; nuc2 < 4; nuc2++) {
            c1 = p7P_CODON3(nuc1, nuc2, x);
            c1 = p7P_MINIDX(c1, p7P_DEGEN_C);

            TMP_SC = SSX2(k, p7S_ATAC, nuc1, nuc2) + signal_scores[p7S_ATAC] + p7P_MSC_CODON(gm_fs, sub_k, c1);
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

  gx->M = M;
  gx->L = L;

  return eslOK;

}


int
p7_splicevitebi_translated_semiglobal_trace(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, const P7_FS_PROFILE *gm_fs, const P7_GMX *gx, P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end, float *ali_score)
{
  int          M   = k_end - k_start + 1;
  int          L   = i_end - i_start + 1;
  int          i   = L;     /* position in seq (1..L)         */
  int          k   = 0;     /* position in model (1..M)       */
  int          t,u,v,w,x;
  int          donor_idx;
  int          c, c0, c1, c2, c3;
  int          sprv,scur,snxt;
  int          sub_i, sub_d, sub_k;
  float        emit, emit0, emit1, emit2;
  float      **dp  = gx->dp;    /* so {MDI}MX() macros work       */
  float       *xmx = gx->xmx;   /* so XMX() macro works           */
  float        tol = 1e-5;  /* floating point "equality" test */
  float const *tsc = gm_fs->tsc; 
  int     status;

#if eslDEBUGLEVEL > 0
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i+i_start-1, 0)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i+i_start-1, 0)) != eslOK) return status;

  snxt = p7T_X; 
  sprv = p7T_C;
  
 
  while (sprv != p7T_S) {
    switch (sprv) {
    case p7T_C:     /* C(i) comes from C(i-1) or E(i) */

      if      (XMX_SP(i, p7G_C) < XMX_SP(i-2, p7G_C) || XMX_SP(i, p7G_C) < XMX_SP(i-1, p7G_C))                     scur = p7T_C; 
      else if (XMX_SP(i,p7G_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
      else if (esl_FCompare_old(XMX_SP(i, p7G_C), XMX_SP(i-3, p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP], tol) == eslOK) scur = p7T_C;
      else if (esl_FCompare_old(XMX_SP(i, p7G_C), XMX_SP(i,   p7G_E) + gm_fs->xsc[p7P_E][p7P_MOVE], tol) == eslOK) scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);

      break;

    case p7T_E:     /* E connects from any M state. k set here */
      if (XMX_SP(i, p7G_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
      *ali_score = XMX_SP(i, p7G_E);
      for (k = M; k >= 1; k--) {
        if (esl_FCompare_old(XMX_SP(i, p7G_E), MMX_SP(i,k), tol) == eslOK) { scur = p7T_M; break; }
        if (esl_FCompare_old(XMX_SP(i, p7G_E), DMX_SP(i,k), tol) == eslOK) { scur = p7T_D; break; }
      }
      if (k == 0) ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
 
      break;

    case p7T_M:         /* M connects from i-1,k-1, or B */
      if (MMX_SP(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);
     
      sub_i = i_start + i - 1; 
      
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i-2])) v = sub_dsq[sub_i-2];
      else                                              v = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i-1])) w = sub_dsq[sub_i-1];
      else                                              w = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i]))   x = sub_dsq[sub_i];
      else                                              x = p7P_MAXCODONS; 

      c3 = p7P_CODON3(v, w, x);
      c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

      sub_k  = k_start + k - 1;
      emit = p7P_MSC_CODON(gm_fs, sub_k, c3);
      
      if      (esl_FCompare_old(MMX_SP(i,k), MMX_SP(i-3, k-1) + TSC(p7P_MM, sub_k-1) + emit, tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(MMX_SP(i,k), IMX_SP(i-3, k-1) + TSC(p7P_IM, sub_k-1) + emit, tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare_old(MMX_SP(i,k), DMX_SP(i-3, k-1) + TSC(p7P_DM, sub_k-1) + emit, tol) == eslOK) scur = p7T_D;      
      else if (esl_FCompare_old(MMX_SP(i,k), XMX_SP(i-3, p7G_B)                      + emit, tol) == eslOK) scur = p7T_B;
      else if (esl_FCompare_old(MMX_SP(i,k), PMX_SP(i-3, k-1) + TSC_P                + emit, tol) == eslOK) scur = p7T_P;
      else ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);  
      
      k--; i-=3;
      break;

    case p7T_D:         /* D connects from M,D at i,k-1 */
      if (DMX_SP(i, k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);
     
      sub_k  = k_start + k - 1; 
      if      (esl_FCompare_old(DMX_SP(i,k), MMX_SP(i, k-1) + TSC(p7P_MD, sub_k-1), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(DMX_SP(i,k), DMX_SP(i, k-1) + TSC(p7P_DD, sub_k-1), tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      k--;
      break;

    case p7T_I:         /* I connects from M,I at i-1,k*/
      if (IMX_SP(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);

      sub_k  = k_start + k - 1;
      if      (esl_FCompare_old(IMX_SP(i,k), MMX_SP(i-3,k) + TSC(p7P_MI, sub_k), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(IMX_SP(i,k), IMX_SP(i-3,k) + TSC(p7P_II, sub_k), tol) == eslOK) scur = p7T_I;
      else ESL_EXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      i-=3;
      break;

    case p7T_P:
      *ali_score -= TSC_P;
      scur = snxt;
      k--; i = donor_idx - c - 2;
      break;

    case p7T_N:         /* N connects from S, N */
      if (XMX_SP(i, p7G_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = ( (i == 0) ? p7T_S : p7T_N);
      break;

    case p7T_B:         /* B connects from N, J */
      if (XMX_SP(i,p7G_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);
      *ali_score -= XMX_SP(i,p7G_B);
      if      (esl_FCompare_old(XMX_SP(i,p7G_B), XMX_SP(i, p7G_N) + gm_fs->xsc[p7P_N][p7P_MOVE], tol) == eslOK) scur = p7T_N;
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
      sub_d = i_start + donor_idx - 1;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_d-3])) t = sub_dsq[sub_d-3];
      else                                                          t = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_d-2])) u = sub_dsq[sub_d-2];
      else                                                          u = p7P_MAXCODONS;

      sub_i = i_start + i - 1;
      /* nucleotides upstream of acceptor */
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i-2])) v = sub_dsq[sub_i-2];
      else                                                  v = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i-1])) w = sub_dsq[sub_i-1];
      else                                                  w = p7P_MAXCODONS;
      if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[sub_i]))   x = sub_dsq[sub_i];
      else                                                  x = p7P_MAXCODONS; 
     
      sub_k = k_start + k - 1;
      c0 = p7P_CODON3(v, w, x);
      c0 = p7P_MINIDX(c0, p7P_DEGEN_C);
      emit0 = p7P_MSC_CODON(gm_fs, sub_k, c0);

      c1 = p7P_CODON3(u, w, x);
      c1 = p7P_MINIDX(c1, p7P_DEGEN_C);
      emit1 = p7P_MSC_CODON(gm_fs, sub_k, c1);

      c2 = p7P_CODON3(t, u, x);
      c2 = p7P_MINIDX(c2, p7P_DEGEN_C);
      emit2 = p7P_MSC_CODON(gm_fs, sub_k, c2);
            
      if(SIGNAL(sub_dsq[sub_d-1], sub_dsq[sub_d]) == DONOR_GT) {
        
        if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-2,k-1) + pli->signal_scores[p7S_GTAG] + emit0, tol) == eslOK) { snxt = p7T_M; c = 0; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-2,k-1) + pli->signal_scores[p7S_GTAG] + emit0, tol) == eslOK) { snxt = p7T_D; c = 0; }
        else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-3,k-1) + pli->signal_scores[p7S_GTAG] + emit1, tol) == eslOK) { snxt = p7T_M; c = 1; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-3,k-1) + pli->signal_scores[p7S_GTAG] + emit1, tol) == eslOK) { snxt = p7T_D; c = 1; }
        else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-4,k-1) + pli->signal_scores[p7S_GTAG] + emit2, tol) == eslOK) { snxt = p7T_M; c = 2; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-4,k-1) + pli->signal_scores[p7S_GTAG] + emit2, tol) == eslOK) { snxt = p7T_D; c = 2; }
        else ESL_EXCEPTION(eslFAIL, "P at k=%d,i=%d couldn't be traced", k,i);
      }
      else if (SIGNAL(sub_dsq[sub_d-1], sub_dsq[sub_d]) == DONOR_GC) {
        if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-2,k-1) + pli->signal_scores[p7S_GCAG] + emit0, tol) == eslOK) { snxt = p7T_M; c = 0; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-2,k-1) + pli->signal_scores[p7S_GCAG] + emit0, tol) == eslOK) { snxt = p7T_D; c = 0; }
        else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-3,k-1) + pli->signal_scores[p7S_GCAG] + emit1, tol) == eslOK) { snxt = p7T_M; c = 1; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-3,k-1) + pli->signal_scores[p7S_GCAG] + emit1, tol) == eslOK) { snxt = p7T_D; c = 1; }
        else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(donor_idx-4,k-1) + pli->signal_scores[p7S_GCAG] + emit2, tol) == eslOK) { snxt = p7T_M; c = 2; }
        else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(donor_idx-4,k-1) + pli->signal_scores[p7S_GCAG] + emit2, tol) == eslOK) { snxt = p7T_D; c = 2; }
        else ESL_EXCEPTION(eslFAIL, "P at k=%d,i=%d couldn't be traced", k,i); 
      }
      else if (SIGNAL(sub_dsq[sub_d-1], sub_dsq[sub_d]) == DONOR_AT) {
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

    if ((status = p7_trace_fs_Append(tr, scur, k_start+k-1, i_start+i-1, c)) != eslOK) return status;
    
    /* For NC, we had to defer i decrement. */
    if ( (scur == p7T_N || scur == p7T_C) && scur == sprv) i--; 
   
    sprv = scur;
  } /* end traceback, at S state */

  tr->M = M;
  tr->L = L;

  return p7_trace_fs_Reverse(tr);
}



