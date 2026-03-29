/* Spliced Viterbi algorithms and traceback.
 *
 * Three DP variants are provided, covering the three ways a sub-sequence
 * region can be aligned relative to the splice graph:
 *
 *   Global     -- global entry (k=1) and exit (k=M)
 *   ExtendDown -- global entry, semi-global exit (any k,i)
 *   ExtendUp   -- semi-global entry (any k,i), global exit
 *
 * All three use log-space probability DP and share a common P-state
 * mechanism: donor-site scores are accumulated in the SPLICE_SCORES
 * arrays (SSX macros) during the forward pass, and acceptor-site
 * scores from those arrays are read back to compute P(i,k) at each
 * position.
 *
 * Contents:
 *    1. p7_GViterbi_SplicedGlobal()
 *    2. p7_GViterbi_SplicedExtendDown()
 *    3. p7_GViterbi_SplicedExtendUp()
 *    4. p7_GViterbi_SplicedTrace()
 *    5. Benchmark driver.
 *    6. Unit tests.
 *    7. Test driver.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_splice.h"

#define TSC_P log(4.58e-5)

/* Function:  p7_GViterbi_SplicedGlobal()
 * Synopsis:  Fully global translated spliced Viterbi algorithm
 *
 * Purpose:   For finding the maxiumum scoring splice site between two 
 *            exons. Algins globally from poistion <i_start> to <i_end> 
 *            on the <sub_dsq> and from <k_start> to <k_end> on the 
 *            <gm_tr>. The DP matrix <gx> must include room for the 
 *            standard core model stats <M, I, D> as well as the splice 
 *            state <P>. The <P> state acts as a modiifed <M> state 
 *            that emits a codon that is made of either two nucleotides 
 *            from before the donor site and one from after the acceptor 
 *            site <C2>, one nucleotide from from before the donor site 
 *            and two from after the acceptor <C1>, or from the three 
 *            after the acceptor <C0>. 
 *
 *            Potetnial donor sites scores are recorded in the  <score> 
 *            matrix (63*M) via the SSX macro. Splice siginal scores are
 *            sotored in <signal_scores> array and <acceptor_..> and 
 *            <donor_..> arrays return 0 for valid sites and -inf for 
 *            all others. 
 *
 * Args:      sub_dsq       - nucleotide sequence 
 *            gm_tr         - a codon profile.
 *            gx            - DP matrix with room for an MxL alignment
 *            P_scores      - storage for max P state scores
 *            signal_scores - scores for splice site signals
 *            i_start       - start poition on the <sub_dsq>
 *            i_end         - end poition on the <sub_dsq>
 *            k_start       - start poition on the <gm_tr>
 *            k_end         - end poition on the <gm_tr>
 *            min_intron    - minimum intron length
 *
 * Return:    <eslOK> on success.
 */
int
p7_GViterbi_SplicedGlobal(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, P7_GMX *gx, float **P_scores, const float *signal_scores, int i_start, int i_end, int k_start, int k_end, int min_intron)
{
  
  float const *tsc  = gm_tr->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  int          L = (i_end - i_start + 1);
  int          M = (k_end - k_start + 1);
  int          C2[5] = {0, 0, 0, 0, 0}; 
  int          C1[5] = {0, 0, 0, 0, 0};
  int          C0;
  int          i,k;
  int          r, s, t, u;
  int          v, w, x;
  int          nuc1, nuc3;
  int          loop_end;
  int          sub_i,sub_k;
  float        TMP_SC;
  int          acc0, acc1, acc2;
  int          don0, don1, don2;
  int          don_sig; 
  float const *rsc_c0;
  float const *rsc_c1[5];

  if(gm_tr->codon_lengths != 1) ESL_EXCEPTION(eslEINVAL, "proflie not allocated for 1 codon length");

  /* Note on variable names 
   * ..0 are for splice codons where zero nucleortides come from before the donor,
   * ..1 are for splice codons where one nucleortide comes from before the donor,
   * ..2 are for splice codons where two nucleortides come from before the donor,
   */

  acc0 = acc1 = acc2 = -1;
  don0 = don1 = don2 = -1;

  /* Reset splice site score storage */
  for(k = 0; k < M; k++) {
    for(i = 0; i < SIGNAL_MEM_SIZE; i++)
      P_scores[k][i] = -eslINFINITY;
  }

  /* Initialization of the zero row.  */
  XMX_SP(0,p7G_N) = 0.;                                     /* S->N, p=1            */
  XMX_SP(0,p7G_B) = gm_tr->xsc[p7P_N][p7P_MOVE];            /* S->N->B, no N-tail   */
  XMX_SP(0,p7G_E) = XMX_SP(0,p7G_C) = XMX_SP(0,p7G_J) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_SP(0,k) = DMX_SP(0,k) = IMX_SP(0,k) = PMX_SP(0,k) = -eslINFINITY;

  /*Special cases for the first 2 rows (no codons yet)*/
  v = w = x = -1;
  for(i = 1; i <= 2; i++)
  {
    w = x;
    sub_i = i_start + i - 1;
    /* if new nucleotide is not A,C,G, or T set it to placeholder value */
    if(sub_dsq[sub_i] < p7P_MAXNUC) x = sub_dsq[sub_i];
    else                         x = p7P_MAXCODONS1;

    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    for (k = 0; k <= M; k++)
      MMX_SP(i,k) = DMX_SP(i,k) = IMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;
  }

  /*Special cases for the rows 3-min_intron+2 (no donor site look backs)*/
  loop_end = ESL_MIN(L, min_intron+2);
  for (i = 3; i <= loop_end; i++)
  {
    v = w;
    w = x;
    sub_i = i_start + i - 1;
    if(sub_dsq[sub_i] < p7P_MAXNUC) x = sub_dsq[sub_i];
    else                         x = p7P_MAXCODONS1; 

    C0 = p7P_CODON3_FS1(v, w, x);
    C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);
    rsc_c0 = gm_tr->rsc[C0];

	/* get acceptor site */
    acc0 = acc1;
    acc1 = acc2;
    if      (v >= p7P_MAXNUC || w >= p7P_MAXNUC) acc2 = -1;
    else if (SIGNAL(v,w) == ACCEPT_AG)     acc2 = ACCEPT_AG;
    else if (SIGNAL(v,w) == ACCEPT_AC)     acc2 = ACCEPT_AC;
    else                                   acc2 = -1;

    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

	/* Global model entry at i = 3, k = 1 */
    sub_k = k_start;
	if(i == 3) MMX_SP(i,1) = XMX_SP(i-3,p7G_B) + rsc_c0[sub_k];
	else       MMX_SP(i,1) = -eslINFINITY;

	IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,sub_k),
                          IMX_SP(i-3,1) + TSC(p7P_II,sub_k));

	/* No I state for stop codons */
	if(rsc_c0[sub_k] == -eslINFINITY) IMX_SP(i,1) = -eslINFINITY;

    DMX_SP(i,1) = PMX_SP(i,1) = -eslINFINITY;

    for (k = 2; k < M; k++) {
      sub_k = k_start + k -1;

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,sub_k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,sub_k-1),
                            DMX_SP(i-3,k-1) + TSC(p7P_DM,sub_k-1))) + rsc_c0[sub_k];
       
	  IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,sub_k),
                            IMX_SP(i-3,k) + TSC(p7P_II,sub_k));

	  /* No I state for stop codons */
      if(rsc_c0[sub_k] == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,sub_k-1));

      PMX_SP(i,k) = -eslINFINITY;
    }

    sub_k = k_start + M -1;
    
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,sub_k-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,sub_k-1),
                          DMX_SP(i-3,M-1) + TSC(p7P_DM,sub_k-1)))+ rsc_c0[sub_k];

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,sub_k-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,sub_k-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;

  }


  /* Example Indexing for r, s, t, u:
   * i_start = 5;
   * min_intron = 13;
   * s = sub_dsq[5];
   * t = sub_dsq[6];
   * u = sub_dsq[7];
   * ...
   * i = min_intron+3;        // 16; 
   * sub_i = i_start + i - 1; // 20;
   * r = s;                   // sub_dsq[5];
   * s = t;                   // sub_dsq[6];
   * t = u;                   // sub_dsq[7];
   * u = sub_dsq[sub_i-min_intron+1]; // sub_dsq[8];
   */

  if(sub_dsq[i_start] < p7P_MAXNUC)   s = sub_dsq[i_start];
  else                             s = p7P_MAXCODONS1;

  if(sub_dsq[i_start+1] < p7P_MAXNUC) t = sub_dsq[i_start+1];
  else                             t = p7P_MAXCODONS1;

  if(sub_dsq[i_start+2] < p7P_MAXNUC) u = sub_dsq[i_start+2];
  else                             u = p7P_MAXCODONS1;

  /* Main DP recursion */
  for (i = min_intron+3; i <= L; i++)
  {
    /* get nucleotides and codon */
    sub_i = i_start + i - 1;

	r = s;
    s = t;
    t = u;
	if(sub_dsq[sub_i-min_intron+1] < p7P_MAXNUC) u = sub_dsq[sub_i-min_intron+1];
	else                                      u = p7P_MAXCODONS1;
 
    v = w;
    w = x;

    if(sub_dsq[sub_i] < p7P_MAXNUC) x = sub_dsq[sub_i];
    else                         x = p7P_MAXCODONS1;

	/* Codon indexing for unsplit codon */
    C0 = p7P_CODON3_FS1(v, w, x);
    C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);
    rsc_c0 = gm_tr->rsc[C0];

	/* Codon indexing for split codons */
    for (nuc1 = 0; nuc1 < p7P_MAXNUC; nuc1++)
      C1[nuc1] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, w, x), p7P_DEGEN1_C);
    C1[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(p7P_MAXCODONS1, w, x), p7P_DEGEN1_C);

    for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++)
      rsc_c1[nuc1] = gm_tr->rsc[C1[nuc1]];

    nuc3 = ESL_MIN(x, p7P_MAXNUC);

    /* get acceptor sites */
    acc0 = acc1;
    acc1 = acc2;
    if      (SIGNAL(v,w) == ACCEPT_AG) acc2 = ACCEPT_AG;
    else if (SIGNAL(v,w) == ACCEPT_AC) acc2 = ACCEPT_AC;
    else                               acc2 = -1;

	/* get donor sites */
	don0 = don1;
    don1 = don2;
    if      (SIGNAL(t,u) == DONOR_GT) don2 = DONOR_GT;
    else if (SIGNAL(t,u) == DONOR_GC) don2 = DONOR_GC;
    else if (SIGNAL(t,u) == DONOR_AT) don2 = DONOR_AT;
    else                              don2 = -1;

    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    for (k = 1; k < M; k++) {

      sub_k = k_start + k -1;

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,sub_k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,sub_k-1),
                    ESL_MAX(DMX_SP(i-3,k-1) + TSC(p7P_DM,sub_k-1),
                            PMX_SP(i-3,k-1) + TSC_P))) + rsc_c0[sub_k];
      
	  IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,sub_k),
                            IMX_SP(i-3,k) + TSC(p7P_II,sub_k));

	  /* No I state for stop codons */
	  if(rsc_c0[sub_k] == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,sub_k-1));

      /* P state scoring */
      PMX_SP(i,k) = -eslINFINITY;

      if (acc0 >= 0 || acc1 >= 0 || acc2 >= 0) {
        if (acc0 == ACCEPT_AG) {
          TMP_SC = ESL_MAX(SSX0(k, p7S_GTAG) + signal_scores[p7S_GTAG],
                           SSX0(k, p7S_GCAG) + signal_scores[p7S_GCAG]) + rsc_c0[sub_k];
          PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
        } else if (acc0 == ACCEPT_AC) {
          TMP_SC = SSX0(k, p7S_ATAC) + signal_scores[p7S_ATAC] + rsc_c0[sub_k];
          PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
        }

        if (acc1 == ACCEPT_AG) {
          for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
            TMP_SC = ESL_MAX(SSX1(k, p7S_GTAG, nuc1) + signal_scores[p7S_GTAG],
                             SSX1(k, p7S_GCAG, nuc1) + signal_scores[p7S_GCAG]) + rsc_c1[nuc1][sub_k];
            PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
          }
        } else if (acc1 == ACCEPT_AC) {
          for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
            TMP_SC = SSX1(k, p7S_ATAC, nuc1) + signal_scores[p7S_ATAC] + rsc_c1[nuc1][sub_k];
            PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
          }
        }

        if (acc2 == ACCEPT_AG) {
          TMP_SC = ESL_MAX(SSX2(k, p7S_GTAG, nuc3) + signal_scores[p7S_GTAG],
                           SSX2(k, p7S_GCAG, nuc3) + signal_scores[p7S_GCAG]);
          PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
        } else if (acc2 == ACCEPT_AC) {
          TMP_SC = SSX2(k, p7S_ATAC, nuc3) + signal_scores[p7S_ATAC];
          PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
        }
      }
    }
	  
    sub_k = k_start + M -1;
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,sub_k-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,sub_k-1),
                  ESL_MAX(DMX_SP(i-3,M-1) + TSC(p7P_DM,sub_k-1),
                          PMX_SP(i-3,M-1) + TSC_P)))         + rsc_c0[sub_k];
    
    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,sub_k-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,sub_k-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;

	/* Separate k loop for donor sites */
	if (don2 >= 0) {
	 
      /* Get codon indicies for all C2 type split codons */
      for(nuc3 = 0; nuc3 < p7P_MAXNUC; nuc3++)
        C2[nuc3] = p7P_MINIDX(p7P_CODON3_FS1(r, s, nuc3), p7P_DEGEN1_C);
      C2[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(r, s, p7P_MAXCODONS1), p7P_DEGEN1_C);

      don_sig = (don2 == DONOR_GT) ? p7S_GTAG : (don2 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
      for (k = 2; k < M; k++) {
        sub_k  = k_start + k - 1;
        TMP_SC = ESL_MAX(MMX_SP(i-min_intron-3,k-1), DMX_SP(i-min_intron-3,k-1));
        for(nuc3 = 0; nuc3 <= p7P_MAXNUC; nuc3++)
          SSX2(k,don_sig,nuc3) = ESL_MAX(SSX2(k,don_sig,nuc3), TMP_SC + gm_tr->rsc[C2[nuc3]][sub_k]);
      }
    }

	if(don1  >= 0) {
	  nuc1 = ESL_MIN(r, p7P_MAXNUC); 
	  don_sig = (don1 == DONOR_GT) ? p7S_GTAG : (don1 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
	  for (k = 2; k < M; k++) {
	    sub_k  = k_start + k - 1;
		TMP_SC = ESL_MAX(MMX_SP(i-min_intron-3,k-1), DMX_SP(i-min_intron-3,k-1));
		SSX1(k, don_sig, nuc1) = ESL_MAX(SSX1(k, don_sig, nuc1), TMP_SC);
	  }
	}

	if (don0 >= 0) {
      don_sig = (don0 == DONOR_GT) ? p7S_GTAG : (don0 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
      for (k = 2; k < M; k++) {
        sub_k  = k_start + k - 1;
        TMP_SC = ESL_MAX(MMX_SP(i-min_intron-3,k-1), DMX_SP(i-min_intron-3,k-1));
        SSX0(k, don_sig) = ESL_MAX(SSX0(k, don_sig), TMP_SC);
	  }
	}
   
  } // end loop over L

  /*exit from last i and k */
  XMX_SP(L,p7G_E) = ESL_MAX(MMX_SP(L,M), DMX_SP(L,M));
  XMX_SP(L,p7G_C) = XMX_SP(L,p7G_E) + gm_tr->xsc[p7P_E][p7P_MOVE];

  gx->M = M;
  gx->L = L;

  return eslOK;
}

int
p7_GViterbi_SplicedGlobal_Dummy(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, P7_GMX *gx, float **P_scores, const float *signal_scores, int i_start, int i_end, int k_start, int k_end, int min_intron)
{
  
  float const *tsc  = gm_tr->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  int          L = (i_end - i_start + 1);
  int          M = (k_end - k_start + 1);
  int          C0;
  int          i,k;
  int          v, w, x;
  int          loop_end;
  int          sub_i,sub_k;
  float const *rsc_c0;

  if(gm_tr->codon_lengths != 1) ESL_EXCEPTION(eslEINVAL, "proflie not allocated for 1 codon length");

  /* Note on variable names 
   * ..0 are for splice codons where zero nucleortides come from before the donor,
   * ..1 are for splice codons where one nucleortide comes from before the donor,
   * ..2 are for splice codons where two nucleortides come from before the donor,
   */

  /* Reset splice site score storage */
  for(k = 0; k < M; k++) {
    for(i = 0; i < SIGNAL_MEM_SIZE; i++)
      P_scores[k][i] = -eslINFINITY;
  }

  /* Initialization of the zero row.  */
  XMX_SP(0,p7G_N) = 0.;                                     /* S->N, p=1            */
  XMX_SP(0,p7G_B) = gm_tr->xsc[p7P_N][p7P_MOVE];            /* S->N->B, no N-tail   */
  XMX_SP(0,p7G_E) = XMX_SP(0,p7G_C) = XMX_SP(0,p7G_J) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_SP(0,k) = DMX_SP(0,k) = IMX_SP(0,k) = PMX_SP(0,k) = -eslINFINITY;

  /*Special cases for the first 2 rows (no codons yet)*/
  v = w = x = -1;
  for(i = 1; i <= 2; i++)
  {
    w = x;
    sub_i = i_start + i - 1;
    /* if new nucleotide is not A,C,G, or T set it to placeholder value */
    if(sub_dsq[sub_i] < p7P_MAXNUC) x = sub_dsq[sub_i];
    else                         x = p7P_MAXCODONS1;

    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    for (k = 0; k <= M; k++)
      MMX_SP(i,k) = DMX_SP(i,k) = IMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;
  }

  /*Special cases for the rows 3-min_intron+2 (no donor site look backs)*/
  loop_end = ESL_MIN(L, min_intron+2);
  for (i = 3; i <= loop_end; i++)
  {
    v = w;
    w = x;
    sub_i = i_start + i - 1;
    if(sub_dsq[sub_i] < p7P_MAXNUC) x = sub_dsq[sub_i];
    else                         x = p7P_MAXCODONS1; 

    C0 = p7P_CODON3_FS1(v, w, x);
    C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);
    rsc_c0 = gm_tr->rsc[C0];

    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

	/* Global model entry at i = 3, k = 1 */
    sub_k = k_start;
	if(i == 3) MMX_SP(i,1) = XMX_SP(i-3,p7G_B) + rsc_c0[sub_k];
	else       MMX_SP(i,1) = -eslINFINITY;

	IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,sub_k),
                          IMX_SP(i-3,1) + TSC(p7P_II,sub_k));

	/* No I state for stop codons */
	if(rsc_c0[sub_k] == -eslINFINITY) IMX_SP(i,1) = -eslINFINITY;

    DMX_SP(i,1) = PMX_SP(i,1) = -eslINFINITY;

    for (k = 2; k < M; k++) {
      sub_k = k_start + k -1;

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,sub_k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,sub_k-1),
                            DMX_SP(i-3,k-1) + TSC(p7P_DM,sub_k-1))) + rsc_c0[sub_k];
       
	  IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,sub_k),
                            IMX_SP(i-3,k) + TSC(p7P_II,sub_k));

	  /* No I state for stop codons */
      if(rsc_c0[sub_k] == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,sub_k-1));

      PMX_SP(i,k) = -eslINFINITY;
    }

    sub_k = k_start + M -1;
    
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,sub_k-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,sub_k-1),
                          DMX_SP(i-3,M-1) + TSC(p7P_DM,sub_k-1)))+ rsc_c0[sub_k];

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,sub_k-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,sub_k-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;

  }


  /* Main DP recursion */
  for (i = min_intron+3; i <= L; i++)
  {
    /* get nucleotides and codon */
    sub_i = i_start + i - 1;

    v = w;
    w = x;

    if(sub_dsq[sub_i] < p7P_MAXNUC) x = sub_dsq[sub_i];
    else                         x = p7P_MAXCODONS1;

	/* Codon indexing for unsplit codon */
    C0 = p7P_CODON3_FS1(v, w, x);
    C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);
    rsc_c0 = gm_tr->rsc[C0];

    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    for (k = 1; k < M; k++) {

      sub_k = k_start + k -1;

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,sub_k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,sub_k-1),
                    ESL_MAX(DMX_SP(i-3,k-1) + TSC(p7P_DM,sub_k-1),
                            PMX_SP(i-3,k-1) + TSC_P))) + rsc_c0[sub_k];
      
	  IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,sub_k),
                            IMX_SP(i-3,k) + TSC(p7P_II,sub_k));

	  /* No I state for stop codons */
	  if(rsc_c0[sub_k] == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,sub_k-1));

      /* P state scoring */
      PMX_SP(i,k) = -eslINFINITY;

    }
	  
    sub_k = k_start + M -1;
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,sub_k-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,sub_k-1),
                  ESL_MAX(DMX_SP(i-3,M-1) + TSC(p7P_DM,sub_k-1),
                          PMX_SP(i-3,M-1) + TSC_P)))         + rsc_c0[sub_k];
    
    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,sub_k-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,sub_k-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;
   
  } // end loop over L

  /*exit from last i and k */
  XMX_SP(L,p7G_E) = ESL_MAX(MMX_SP(L,M), DMX_SP(L,M));
  XMX_SP(L,p7G_C) = XMX_SP(L,p7G_E) + gm_tr->xsc[p7P_E][p7P_MOVE];

  gx->M = M;
  gx->L = L;

  return eslOK;
}



/* Function:  p7_GViterbi_SplicedExtendDown()
 * Synopsis:  Semi global translated spliced Viterbi algorithm (global start)
 *
 * Purpose:   For finding the maxiumum scoring splice site between an 
 *            exons and a downstream region that may contiain one or 
 *            more other exons. Aligns globally from the start postions 
 *            <i_start> on the <sub_dsq> and <k_start> on the <gm_tr>, 
 *            but can exit at any point at ro before <i_end> on the 
 *            <sub_dsq> and <k_end> on the <gm_tr>. The DP matrix <gx> 
 *            must include room for the standard core model stats <M,I,D> 
 *            as well as the splice state <P>. The <P> state acts as a modiifed <M> state
 *            that emits a codon that is made of either two nucleotides
 *            from before the donor site and one from after the acceptor
 *            site <C2>, one nucleotide from from before the donor site
 *            and two from after the acceptor <C1>, or from the three
 *            after the acceptor <C0>.
 *
 *            Potetnial donor sites scores are recorded in the  <score>
 *            matrix (63*M) via the SSX macro. Splice siginal scores are
 *            sotored in <signal_scores> array and <acceptor_..> and
 *            <donor_..> arrays return 0 for valid sites and -inf for
 *            all others. 
 *             
 * Args:      sub_dsq       - nucleotide sequence 
 *            gm_tr         - a codon profile.
 *            gx            - DP matrix with room for an MxL alignment
 *            P_scores      - storage for max P state scores
 *            signal_scores - scores for splice site signals
 *            i_start       - start poition on the <sub_dsq>
 *            i_end         - end poition on the <sub_dsq>
 *            k_start       - start poition on the <gm_tr>
 *            k_end         - end poition on the <gm_tr>
 *            min_intron    - minimum intron length
 *
 * Return:    <eslOK> on success.
 */
int
p7_GViterbi_SplicedExtendDown(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, P7_GMX *gx, float **P_scores, const float *signal_scores, int i_start, int i_end, int k_start, int k_end, int min_intron)
{
  float const *tsc  = gm_tr->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  int          L = (i_end - i_start + 1);
  int          M = (k_end - k_start + 1);
  int          C2[5] = {0, 0, 0, 0, 0};
  int          C1[5] = {0, 0, 0, 0, 0};
  int          C0;
  int          i,k;
  int          r, s, t, u;
  int          v, w, x;
  int          nuc1, nuc3;
  int          loop_end;
  int          sub_i,sub_k;
  float        TMP_SC;
  int          acc0, acc1, acc2;     /* acceptor signal at codon offsets 0,1,2: ACCEPT_AG, ACCEPT_AC, or -1 */
  int          don0, don1, don2; /* donor signal at codon offsets 0,1,2: DONOR_GT/GC/AT, or -1 */
  int          don_sig;
  float const *rsc_c0;
  float const *rsc_c1[5];

  if(gm_tr->codon_lengths != 1) ESL_EXCEPTION(eslEINVAL, "proflie not allocated for 1 codon length");

  /* Note on variable names
   * ..0 are for splice codons where zero nucleortides come from before the donor,
   * ..1 are for splice codons where one nucleortide comes from before the donor,
   * ..2 are for splice codons where two nucleortides come from before the donor,
   */

  acc0 = acc1 = acc2 = -1;
  don0 = don1 = don2 = -1;

  /* Reset splice site score storage */
  for(k = 0; k < M; k++) {
    for(i = 0; i < SIGNAL_MEM_SIZE; i++)
      P_scores[k][i] = -eslINFINITY;
  }

  /* Initialization of the zero row.  */
  XMX_SP(0,p7G_N) = 0.;                                     /* S->N, p=1            */
  XMX_SP(0,p7G_B) = gm_tr->xsc[p7P_N][p7P_MOVE];            /* S->N->B, no N-tail   */
  XMX_SP(0,p7G_E) = XMX_SP(0,p7G_C) = XMX_SP(0,p7G_J) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_SP(0,k) = DMX_SP(0,k) = IMX_SP(0,k) = PMX_SP(0,k) = -eslINFINITY;

  /*Special cases for the first 2 rows (no codons yet) */
  t = u = v = w = x = -1;
  for(i = 1; i <= 2; i++)
  {
    w = x;
    sub_i = i_start + i - 1;
    if(sub_dsq[sub_i] < p7P_MAXNUC) x = sub_dsq[sub_i];
    else                         x = p7P_MAXCODONS1;

    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    for (k = 0; k <= M; k++)
      MMX_SP(i,k) = DMX_SP(i,k) = IMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;
  }

  /*Special cases for the rows 3-min_intron+2 (no donor site look backs) */
  loop_end = ESL_MIN(L, min_intron+2);
  for (i = 3; i <= loop_end; i++)
  {
    v = w;
    w = x;
    sub_i = i_start + i - 1;
    if(sub_dsq[sub_i] < p7P_MAXNUC) x = sub_dsq[sub_i];
    else                         x = p7P_MAXCODONS1;

    C0 = p7P_CODON3_FS1(v, w, x);
    C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);
    rsc_c0 = gm_tr->rsc[C0];

    /* get acceptor site */
    acc0 = acc1;
    acc1 = acc2;
    if      (v >= p7P_MAXNUC || w >= p7P_MAXNUC) acc2 = -1;
    else if (SIGNAL(v,w) == ACCEPT_AG)     acc2 = ACCEPT_AG;
    else if (SIGNAL(v,w) == ACCEPT_AC)     acc2 = ACCEPT_AC;
    else                                   acc2 = -1;

    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    XMX_SP(i,p7G_E) = -eslINFINITY;

    sub_k = k_start;

    /* Global entry at i = 3, k = 1 */
	if (i == 3) MMX_SP(i,1) = XMX_SP(i-3,p7G_B) + rsc_c0[sub_k];
    else        MMX_SP(i,1) = -eslINFINITY;

    IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,sub_k),
                          IMX_SP(i-3,1) + TSC(p7P_II,sub_k));

    /* No I state for stop codons */
    if(rsc_c0[sub_k] == -eslINFINITY) IMX_SP(i,1) = -eslINFINITY;

    DMX_SP(i,1) = PMX_SP(i,1) = -eslINFINITY;

    XMX_SP(i,p7G_E) = ESL_MAX(XMX_SP(i,p7G_E), MMX_SP(i,1));

    for (k = 2; k < M; k++) {
      sub_k = k_start + k -1;

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,sub_k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,sub_k-1),
                            DMX_SP(i-3,k-1) + TSC(p7P_DM,sub_k-1))) + rsc_c0[sub_k];

      IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,sub_k),
                            IMX_SP(i-3,k) + TSC(p7P_II,sub_k));

      /* No I state for stop codons */
      if(rsc_c0[sub_k] == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,sub_k-1));

      PMX_SP(i,k) = -eslINFINITY;

      XMX_SP(i,p7G_E) = ESL_MAX(XMX_SP(i,p7G_E),
                        ESL_MAX(MMX_SP(i,k), DMX_SP(i,k)));
    }

    sub_k = k_start + M -1;
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,sub_k-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,sub_k-1),
                          DMX_SP(i-3,M-1) + TSC(p7P_DM,sub_k-1))) + rsc_c0[sub_k];

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,sub_k-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,sub_k-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = ESL_MAX(XMX_SP(i,p7G_E),
                      ESL_MAX(MMX_SP(i,M), DMX_SP(i,M)));

    XMX_SP(i,p7G_C) = ESL_MAX(XMX_SP(i-3,p7G_C) + gm_tr->xsc[p7P_C][p7P_LOOP],
                              XMX_SP(i,p7G_E)   + gm_tr->xsc[p7P_E][p7P_MOVE]);

    XMX_SP(i,p7G_J) = -eslINFINITY;

  }

  if(sub_dsq[i_start] < p7P_MAXNUC)   s = sub_dsq[i_start];
  else                             s = p7P_MAXCODONS1;

  if(sub_dsq[i_start+1] < p7P_MAXNUC) t = sub_dsq[i_start+1];
  else                             t = p7P_MAXCODONS1;

  if(sub_dsq[i_start+2] < p7P_MAXNUC) u = sub_dsq[i_start+2];
  else                             u = p7P_MAXCODONS1;

  /* Main DP recursion */
  for (i = min_intron+3; i <= L; i++)
  {
    sub_i = i_start + i - 1;

    r = s;
    s = t;
    t = u;
    if(sub_dsq[sub_i-min_intron+1] < p7P_MAXNUC) u = sub_dsq[sub_i-min_intron+1];
    else                                      u = p7P_MAXCODONS1;

    v = w;
    w = x;

    if(sub_dsq[sub_i] < p7P_MAXNUC) x = sub_dsq[sub_i];
    else                         x = p7P_MAXCODONS1;

    /* Codon indexing */
    C0 = p7P_CODON3_FS1(v, w, x);
    C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);

    for (nuc1 = 0; nuc1 < p7P_MAXNUC; nuc1++)
      C1[nuc1] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, w, x), p7P_DEGEN1_C);
    C1[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(p7P_MAXCODONS1, w, x), p7P_DEGEN1_C);

    nuc3 = ESL_MIN(x, p7P_MAXNUC);

    rsc_c0 = gm_tr->rsc[C0];
    for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++)
      rsc_c1[nuc1] = gm_tr->rsc[C1[nuc1]];

    /* get acceptor sites */
    acc0 = acc1;
    acc1 = acc2;
    if      (SIGNAL(v,w) == ACCEPT_AG) acc2 = ACCEPT_AG;
    else if (SIGNAL(v,w) == ACCEPT_AC) acc2 = ACCEPT_AC;
    else                               acc2 = -1;

    /* get donor sites */
    don0 = don1;
    don1 = don2;
    if      (SIGNAL(t,u) == DONOR_GT) don2 = DONOR_GT;
    else if (SIGNAL(t,u) == DONOR_GC) don2 = DONOR_GC;
    else if (SIGNAL(t,u) == DONOR_AT) don2 = DONOR_AT;
    else                              don2 = -1;


    XMX_SP(i,p7G_N) =  XMX_SP(i,p7G_B) =  -eslINFINITY;

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    XMX_SP(i,p7G_E) = -eslINFINITY;

    for (k = 1; k < M; k++) {

      sub_k = k_start + k -1;

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,sub_k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,sub_k-1),
                    ESL_MAX(DMX_SP(i-3,k-1) + TSC(p7P_DM,sub_k-1),
                            PMX_SP(i-3,k-1) + TSC_P))) + rsc_c0[sub_k];

      IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,sub_k),
                            IMX_SP(i-3,k) + TSC(p7P_II,sub_k));

      /* No I state for stop codons */
      if(rsc_c0[sub_k] == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,sub_k-1));

      /* P state scoring */
      PMX_SP(i,k) = -eslINFINITY;

      if (acc0 >= 0 || acc1 >= 0 || acc2 >= 0) {
        if (acc0 == ACCEPT_AG) {
          TMP_SC = ESL_MAX(SSX0(k, p7S_GTAG) + signal_scores[p7S_GTAG],
                           SSX0(k, p7S_GCAG) + signal_scores[p7S_GCAG]) + rsc_c0[sub_k];
          PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
        } else if (acc0 == ACCEPT_AC) {
          TMP_SC = SSX0(k, p7S_ATAC) + signal_scores[p7S_ATAC] + rsc_c0[sub_k];
          PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
        }
        if (acc1 == ACCEPT_AG) {
          for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
            TMP_SC = ESL_MAX(SSX1(k, p7S_GTAG, nuc1) + signal_scores[p7S_GTAG],
                             SSX1(k, p7S_GCAG, nuc1) + signal_scores[p7S_GCAG]) + rsc_c1[nuc1][sub_k];
            PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
          }
        } else if (acc1 == ACCEPT_AC) {
          for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
            TMP_SC = SSX1(k, p7S_ATAC, nuc1) + signal_scores[p7S_ATAC] + rsc_c1[nuc1][sub_k];
            PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
          }
        }
        if (acc2 == ACCEPT_AG) {
          TMP_SC = ESL_MAX(SSX2(k, p7S_GTAG, nuc3) + signal_scores[p7S_GTAG],
                           SSX2(k, p7S_GCAG, nuc3) + signal_scores[p7S_GCAG]);
          PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
        } else if (acc2 == ACCEPT_AC) {
          TMP_SC = SSX2(k, p7S_ATAC, nuc3) + signal_scores[p7S_ATAC];
          PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
        }
      }

      XMX_SP(i,p7G_E) = ESL_MAX(XMX_SP(i,p7G_E),
                        ESL_MAX(MMX_SP(i,k), DMX_SP(i,k)));

    } // end k loop

    /* Separate k loop for donor sites */
    if (don2 >= 0) {

      /* Get codon indices for all C2 type split codons */
      for(nuc3 = 0; nuc3 < p7P_MAXNUC; nuc3++)
        C2[nuc3] = p7P_MINIDX(p7P_CODON3_FS1(r, s, nuc3), p7P_DEGEN1_C);
      C2[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(r, s, p7P_MAXCODONS1), p7P_DEGEN1_C);

      don_sig = (don2 == DONOR_GT) ? p7S_GTAG : (don2 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
      for (k = 2; k < M; k++) {
        sub_k  = k_start + k - 1;
        TMP_SC = ESL_MAX(MMX_SP(i-min_intron-3,k-1), DMX_SP(i-min_intron-3,k-1));
        for(nuc3 = 0; nuc3 <= p7P_MAXNUC; nuc3++)
          SSX2(k,don_sig,nuc3) = ESL_MAX(SSX2(k,don_sig,nuc3), TMP_SC + gm_tr->rsc[C2[nuc3]][sub_k]);
      }
    }

    if (don1 >= 0) {
      nuc1 = ESL_MIN(r, p7P_MAXNUC);
      don_sig = (don1 == DONOR_GT) ? p7S_GTAG : (don1 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
      for (k = 2; k < M; k++) {
        sub_k  = k_start + k - 1;
        TMP_SC = ESL_MAX(MMX_SP(i-min_intron-3,k-1), DMX_SP(i-min_intron-3,k-1));
        SSX1(k, don_sig, nuc1) = ESL_MAX(SSX1(k, don_sig, nuc1), TMP_SC);
      }
    }

    if (don0 >= 0) {
      don_sig = (don0 == DONOR_GT) ? p7S_GTAG : (don0 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
      for (k = 2; k < M; k++) {
        sub_k  = k_start + k - 1;
        TMP_SC = ESL_MAX(MMX_SP(i-min_intron-3,k-1), DMX_SP(i-min_intron-3,k-1));
        SSX0(k, don_sig) = ESL_MAX(SSX0(k, don_sig), TMP_SC);
      }
    }


    sub_k = k_start + M -1;
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,sub_k-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,sub_k-1),
                  ESL_MAX(DMX_SP(i-3,M-1) + TSC(p7P_DM,sub_k-1),
                          PMX_SP(i-3,M-1) + TSC_P)))         + rsc_c0[sub_k];

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,sub_k-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,sub_k-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = ESL_MAX(XMX_SP(i,p7G_E),
                      ESL_MAX(MMX_SP(i,M), DMX_SP(i,M)));

    XMX_SP(i,p7G_C) = ESL_MAX(XMX_SP(i-3,p7G_C) + gm_tr->xsc[p7P_C][p7P_LOOP],
                              XMX_SP(i,p7G_E)   + gm_tr->xsc[p7P_E][p7P_MOVE]);

    XMX_SP(i,p7G_J) = -eslINFINITY;

  } // end loop over L

  gx->M = M;
  gx->L = L;

  return eslOK;

}


/* Function:  p7_GViterbi_SplicedExtendUp()
 * Synopsis:  Semi global translated spliced Viterbi algorithm (global start)
 *
 * Purpose:   For finding the maxiumum scoring splice site between an
 *            exons and an upstream region that may contiain one or
 *            more other exons. Can enter the model at any point at or 
 *            after <i_start> on the <sub_dsq> and <k_start> on the 
 *            <gm_tr>, but must exit globally ar at <i_end> on the
 *            <sub_dsq> and <k_end> on the <gm_tr>. The DP matrix <gx>
 *            must include room for the standard core model stats <M,I,D>
 *            as well as the splice state <P>. The <P> state acts as a modiifed <M> state
 *            that emits a codon that is made of either two nucleotides
 *            from before the donor site and one from after the acceptor
 *            site <C2>, one nucleotide from from before the donor site
 *            and two from after the acceptor <C1>, or from the three
 *            after the acceptor <C0>.
 *
 *            Potetnial donor sites scores are recorded in the  <score>
 *            matrix (63*M) via the SSX macro. Splice siginal scores are
 *            sotored in <signal_scores> array and <acceptor_..> and
 *            <donor_..> arrays return 0 for valid sites and -inf for
 *            all others. 
 *
 * Args:      sub_dsq       - nucleotide sequence
 *            gm_tr         - a codon profile.
 *            gx            - DP matrix with room for an MxL alignment
 *            P_scores      - storage for max P state scores
 *            signal_scores - scores for splice site signals
 *            i_start       - start poition on the <sub_dsq>
 *            i_end         - end poition on the <sub_dsq>
 *            k_start       - start poition on the <gm_tr>
 *            k_end         - end poition on the <gm_tr>
 *            min_intron    - minimum intron length
 *
 * Return:    <eslOK> on success.
 */
int
p7_GViterbi_SplicedExtendUp(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, P7_GMX *gx, float **P_scores, const float *signal_scores, int i_start, int i_end, int k_start, int k_end, int min_intron)
{
  float const *tsc  = gm_tr->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  int          L = (i_end - i_start + 1);
  int          M = (k_end - k_start + 1);
  int          C2[5] = {0, 0, 0, 0, 0};
  int          C1[5] = {0, 0, 0, 0, 0};
  int          C0;
  int          i,k;
  int          r, s, t, u;
  int          v, w, x;
  int          nuc1, nuc3;
  int          loop_end;
  int          sub_i,sub_k;
  float        TMP_SC;
  int          acc0, acc1, acc2;     /* acceptor signal at codon offsets 0,1,2: ACCEPT_AG, ACCEPT_AC, or -1 */
  int          don0, don1, don2; /* donor signal at codon offsets 0,1,2: DONOR_GT/GC/AT, or -1 */
  int          don_sig;
  float const *rsc_c0;
  float const *rsc_c1[5];

  if(gm_tr->codon_lengths != 1) ESL_EXCEPTION(eslEINVAL, "proflie not allocated for 1 codon length");

  /* Note on variable names
   * ..0 are for splice codons where zero nucleortides come from before the donor,
   * ..1 are for splice codons where one nucleortide comes from before the donor,
   * ..2 are for splice codons where two nucleortides come from before the donor,
   */

  acc0 = acc1 = acc2 = -1;
  don0 = don1 = don2 = -1;

  /* Reset splice site score storage */
  for(k = 0; k < M; k++) {
    for(i = 0; i < SIGNAL_MEM_SIZE; i++)
      P_scores[k][i] = -eslINFINITY;
  }

  /* Initialization of the zero row.  */
  XMX_SP(0,p7G_N) = 0.;                                     /* S->N, p=1            */
  XMX_SP(0,p7G_B) = gm_tr->xsc[p7P_N][p7P_MOVE];            /* S->N->B, no N-tail   */
  XMX_SP(0,p7G_E) = XMX_SP(0,p7G_C) = XMX_SP(0,p7G_J) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_SP(0,k) = DMX_SP(0,k) = IMX_SP(0,k) = PMX_SP(0,k) = -eslINFINITY;

  /*Special cases for the first 2 rows */
  t = u = v = w = x = -1;
  for(i = 1; i <= 2; i++)
  {
    w = x;
    sub_i = i_start + i - 1;
    if(sub_dsq[sub_i] < p7P_MAXNUC) x = sub_dsq[sub_i];
    else                         x = p7P_MAXCODONS1;

    XMX_SP(i,p7G_N) =  0;
    XMX_SP(i,p7G_B) =  gm_tr->xsc[p7P_N][p7P_MOVE];
    XMX_SP(i,p7G_J) = -eslINFINITY;

    for (k = 0; k <= M; k++)
      MMX_SP(i,k) = DMX_SP(i,k) = IMX_SP(i,k) = PMX_SP(i,k) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = -eslINFINITY;
  }

  /*Special cases for the rows 3-min_intron+2 (no donor site look backs) */
  loop_end = ESL_MIN(L, min_intron+2);
  for (i = 3; i <= loop_end; i++)
  {
    v = w;
    w = x;
    sub_i = i_start + i - 1;
    if(sub_dsq[sub_i] < p7P_MAXNUC) x = sub_dsq[sub_i];
    else                         x = p7P_MAXCODONS1;

    C0 = p7P_CODON3_FS1(v, w, x);
    C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);
    rsc_c0 = gm_tr->rsc[C0];

    /* get acceptor site */
    acc0 = acc1;
    acc1 = acc2;
    if      (v >= p7P_MAXNUC || w >= p7P_MAXNUC) acc2 = -1;
    else if (SIGNAL(v,w) == ACCEPT_AG)     acc2 = ACCEPT_AG;
    else if (SIGNAL(v,w) == ACCEPT_AC)     acc2 = ACCEPT_AC;
    else                                   acc2 = -1;

    XMX_SP(i,p7G_N) = XMX_SP(i-3,p7G_N) + gm_tr->xsc[p7P_N][p7P_LOOP];
    XMX_SP(i,p7G_B) = XMX_SP(i,p7G_N)   + gm_tr->xsc[p7P_N][p7P_MOVE];

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    for (k = 1; k < M; k++) {
      sub_k = k_start + k -1;

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,sub_k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,sub_k-1),
                    ESL_MAX(DMX_SP(i-3,k-1) + TSC(p7P_DM,sub_k-1),
                            XMX_SP(i-3,p7G_B)))) + rsc_c0[sub_k];

      IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,sub_k),
                            IMX_SP(i-3,k) + TSC(p7P_II,sub_k));

      /* No I state for stop codons */
      if(rsc_c0[sub_k] == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,sub_k-1));

      PMX_SP(i,k) = -eslINFINITY;
    }

    sub_k = k_start + M -1;
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,sub_k-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,sub_k-1),
                  ESL_MAX(DMX_SP(i-3,M-1) + TSC(p7P_DM,sub_k-1),
                          XMX_SP(i-3,p7G_B)))) + rsc_c0[sub_k];

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,sub_k-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,sub_k-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;

  }

  if(sub_dsq[i_start] < p7P_MAXNUC)   s = sub_dsq[i_start];
  else                             s = p7P_MAXCODONS1;

  if(sub_dsq[i_start+1] < p7P_MAXNUC) t = sub_dsq[i_start+1];
  else                             t = p7P_MAXCODONS1;

  if(sub_dsq[i_start+2] < p7P_MAXNUC) u = sub_dsq[i_start+2];
  else                             u = p7P_MAXCODONS1;

  /* Main DP recursion */
  for (i = min_intron+3; i <= L; i++)
  {
    sub_i = i_start + i - 1;

    r = s;
    s = t;
    t = u;
    if(sub_dsq[sub_i-min_intron+1] < p7P_MAXNUC) u = sub_dsq[sub_i-min_intron+1];
    else                                      u = p7P_MAXCODONS1;

    v = w;
    w = x;

    if(sub_dsq[sub_i] < p7P_MAXNUC) x = sub_dsq[sub_i];
    else                         x = p7P_MAXCODONS1;

    /* Codon indexing */
    C0 = p7P_CODON3_FS1(v, w, x);
    C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);

    for (nuc1 = 0; nuc1 < p7P_MAXNUC; nuc1++)
      C1[nuc1] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, w, x), p7P_DEGEN1_C);
    C1[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(p7P_MAXCODONS1, w, x), p7P_DEGEN1_C);

    nuc3 = ESL_MIN(x, p7P_MAXNUC);

    rsc_c0 = gm_tr->rsc[C0];
    for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++)
      rsc_c1[nuc1] = gm_tr->rsc[C1[nuc1]];

    /* get acceptor sites */
    acc0 = acc1;
    acc1 = acc2;
    if      (SIGNAL(v,w) == ACCEPT_AG) acc2 = ACCEPT_AG;
    else if (SIGNAL(v,w) == ACCEPT_AC) acc2 = ACCEPT_AC;
    else                               acc2 = -1;

    /* get donor sites */
    don0 = don1;
    don1 = don2;
    if      (SIGNAL(t,u) == DONOR_GT) don2 = DONOR_GT;
    else if (SIGNAL(t,u) == DONOR_GC) don2 = DONOR_GC;
    else if (SIGNAL(t,u) == DONOR_AT) don2 = DONOR_AT;
    else                              don2 = -1;


    XMX_SP(i,p7G_N) = XMX_SP(i-3,p7G_N) + gm_tr->xsc[p7P_N][p7P_LOOP];
    XMX_SP(i,p7G_B) = XMX_SP(i,p7G_N)   + gm_tr->xsc[p7P_N][p7P_MOVE];

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = PMX_SP(i,0) = -eslINFINITY;

    for (k = 1; k < M; k++) {

      sub_k = k_start + k -1;

      MMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k-1) + TSC(p7P_MM,sub_k-1),
                    ESL_MAX(IMX_SP(i-3,k-1) + TSC(p7P_IM,sub_k-1),
                    ESL_MAX(DMX_SP(i-3,k-1) + TSC(p7P_DM,sub_k-1),
                    ESL_MAX(XMX_SP(i-3,p7G_B),
                            PMX_SP(i-3,k-1) + TSC_P)))) + rsc_c0[sub_k];

      IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,sub_k),
                            IMX_SP(i-3,k) + TSC(p7P_II,sub_k));

      /* No I state for stop codons */
      if(rsc_c0[sub_k] == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,sub_k-1));

      /* P state scoring */
      PMX_SP(i,k) = -eslINFINITY;

      if (acc0 >= 0 || acc1 >= 0 || acc2 >= 0) {
        if (acc0 == ACCEPT_AG) {
          TMP_SC = ESL_MAX(SSX0(k, p7S_GTAG) + signal_scores[p7S_GTAG],
                           SSX0(k, p7S_GCAG) + signal_scores[p7S_GCAG]) + rsc_c0[sub_k];
          PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
        } else if (acc0 == ACCEPT_AC) {
          TMP_SC = SSX0(k, p7S_ATAC) + signal_scores[p7S_ATAC] + rsc_c0[sub_k];
          PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
        }
        if (acc1 == ACCEPT_AG) {
          for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
            TMP_SC = ESL_MAX(SSX1(k, p7S_GTAG, nuc1) + signal_scores[p7S_GTAG],
                             SSX1(k, p7S_GCAG, nuc1) + signal_scores[p7S_GCAG]) + rsc_c1[nuc1][sub_k];
            PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
          }
        } else if (acc1 == ACCEPT_AC) {
          for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
            TMP_SC = SSX1(k, p7S_ATAC, nuc1) + signal_scores[p7S_ATAC] + rsc_c1[nuc1][sub_k];
            PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
          }
        }
        if (acc2 == ACCEPT_AG) {
          TMP_SC = ESL_MAX(SSX2(k, p7S_GTAG, nuc3) + signal_scores[p7S_GTAG],
                           SSX2(k, p7S_GCAG, nuc3) + signal_scores[p7S_GCAG]);
          PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
        } else if (acc2 == ACCEPT_AC) {
          TMP_SC = SSX2(k, p7S_ATAC, nuc3) + signal_scores[p7S_ATAC];
          PMX_SP(i,k) = ESL_MAX(PMX_SP(i,k), TMP_SC);
        }
      }

    } // end k loop

    /* Separate k loop for donor sites */
    if (don2 >= 0) {

      /* Get codon indices for all C2 type split codons */
      for(nuc3 = 0; nuc3 < p7P_MAXNUC; nuc3++)
        C2[nuc3] = p7P_MINIDX(p7P_CODON3_FS1(r, s, nuc3), p7P_DEGEN1_C);
      C2[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(r, s, p7P_MAXCODONS1), p7P_DEGEN1_C);

      don_sig = (don2 == DONOR_GT) ? p7S_GTAG : (don2 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
      for (k = 2; k < M; k++) {
        sub_k  = k_start + k - 1;
        TMP_SC = ESL_MAX(MMX_SP(i-min_intron-3,k-1), DMX_SP(i-min_intron-3,k-1));
        for(nuc3 = 0; nuc3 <= p7P_MAXNUC; nuc3++)
          SSX2(k,don_sig,nuc3) = ESL_MAX(SSX2(k,don_sig,nuc3), TMP_SC + gm_tr->rsc[C2[nuc3]][sub_k]);
      }
    }

    if (don1 >= 0) {
      nuc1 = ESL_MIN(r, p7P_MAXNUC);
      don_sig = (don1 == DONOR_GT) ? p7S_GTAG : (don1 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
      for (k = 2; k < M; k++) {
        sub_k  = k_start + k - 1;
        TMP_SC = ESL_MAX(MMX_SP(i-min_intron-3,k-1), DMX_SP(i-min_intron-3,k-1));
        SSX1(k, don_sig, nuc1) = ESL_MAX(SSX1(k, don_sig, nuc1), TMP_SC);
      }
    }

    if (don0 >= 0) {
      don_sig = (don0 == DONOR_GT) ? p7S_GTAG : (don0 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
      for (k = 2; k < M; k++) {
        sub_k  = k_start + k - 1;
        TMP_SC = ESL_MAX(MMX_SP(i-min_intron-3,k-1), DMX_SP(i-min_intron-3,k-1));
        SSX0(k, don_sig) = ESL_MAX(SSX0(k, don_sig), TMP_SC);
      }
    }

    sub_k = k_start + M -1;
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i-3,M-1) + TSC(p7P_MM,sub_k-1),
                  ESL_MAX(IMX_SP(i-3,M-1) + TSC(p7P_IM,sub_k-1),
                  ESL_MAX(DMX_SP(i-3,M-1) + TSC(p7P_DM,sub_k-1),
                  ESL_MAX(XMX_SP(i-3,p7G_B),
                          PMX_SP(i-3,M-1) + TSC_P))))         + rsc_c0[sub_k];

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,sub_k-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,sub_k-1));

    PMX_SP(i,M) = -eslINFINITY;

    XMX_SP(i,p7G_E) = XMX_SP(i,p7G_C) = XMX_SP(i,p7G_J) = -eslINFINITY;

  } // end loop over L

  /*Global exit from last i and k */
  XMX_SP(L,p7G_E) = ESL_MAX(MMX_SP(L,M), DMX_SP(L,M));
  XMX_SP(L,p7G_C) = XMX_SP(L,p7G_E) + gm_tr->xsc[p7P_E][p7P_MOVE];

  gx->M = M;
  gx->L = L;

  return eslOK;

}


/* Function:  p7_GViterbi_SplicedTrace()
 * Synopsis:  Create a trace fot any of the translted spliced vitebi algorithms
 *
 * Purpose:   Create a trace that includes the exons(s) and any splice site(s) 
 *            from the filled translated spliced viterbi matrix <gx> and the 
 *            donor site index matrix <lookback>. Returns the filled trance <tr> 
 *
 * Args:      sub_dsq       - nucleotide sequence
 *            gm_tr         - a codon profile.
 *            gx            - filled spliced viterbi DP matrix
 *            signal_scores - scores for splice site signals 
 *            tr            - trace to fil
 *            i_start       - start poition on the <sub_dsq>
 *            i_end         - end poition on the <sub_dsq>
 *            k_start       - start poition on the <gm_tr>
 *            k_end         - end poition on the <gm_tr>
 *            min_intron    - minimum intron length
 *
 * Return:    <eslOK> on success. 
 *            <eslFAIL> if even the optimal path has zero probability;
 *            in this case, the trace is set blank (<tr->N = 0>).
 *
 */            
int
p7_GViterbi_SplicedTrace(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, const P7_GMX *gx, const float *signal_scores, P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end, int min_intron)
{

  float const *tsc = gm_tr->tsc;
  float      **dp  = gx->dp;              /* so {MDI}MX() macros work       */
  float       *xmx = gx->xmx;             /* so XMX() macro works           */
  float        tol = 1e-5;                /* floating point "equality" test */
  int          M   = k_end - k_start + 1; /* sub model length               */
  int          L   = i_end - i_start + 1; /* sub seq length                 */
  int          i   = L;                   /* position in seq (1..L)         */
  int          k   = 0;                   /* position in model (1..M)       */
  int          j;
  int          t,u,v,w,x;
  int          c, c0, c1, c2, c3;
  int          sprv,scur,snxt;
  int          sub_i, sub_d, sub_k;
  float        emit, emit0, emit1, emit2;
  int          status;

  if(gm_tr->codon_lengths != 1) ESL_EXCEPTION(eslEINVAL, "proflie not allocated for 1 codon length");

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
      else if (esl_FCompare_old(XMX_SP(i, p7G_C), XMX_SP(i-3, p7G_C) + gm_tr->xsc[p7P_C][p7P_LOOP], tol) == eslOK) scur = p7T_C;
      else if (esl_FCompare_old(XMX_SP(i, p7G_C), XMX_SP(i,   p7G_E) + gm_tr->xsc[p7P_E][p7P_MOVE], tol) == eslOK) scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);

      break;

    case p7T_E:     /* E connects from any M state. k set here */
      if (XMX_SP(i, p7G_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
      for (k = M; k >= 1; k--) {
        if (esl_FCompare_old(XMX_SP(i, p7G_E), MMX_SP(i,k), tol) == eslOK) { scur = p7T_M; break; }
        if (esl_FCompare_old(XMX_SP(i, p7G_E), DMX_SP(i,k), tol) == eslOK) { scur = p7T_D; break; }
      }
      if (k == 0) ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
 
      break;

    case p7T_M:         /* M connects from i-1,k-1, or B */
      if (MMX_SP(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);
     
      sub_i = i_start + i - 1; 
      
      if(sub_dsq[sub_i-2] < p7P_MAXNUC) v = sub_dsq[sub_i-2];
      else                           v = p7P_MAXCODONS1;

      if(sub_dsq[sub_i-1] < p7P_MAXNUC) w = sub_dsq[sub_i-1];
      else                           w = p7P_MAXCODONS1;

      if(sub_dsq[sub_i]   < p7P_MAXNUC) x = sub_dsq[sub_i];
      else                           x = p7P_MAXCODONS1;

      c3 = p7P_CODON3_FS1(v, w, x);
      c3 = p7P_MINIDX(c3, p7P_DEGEN1_C);

      sub_k  = k_start + k - 1;
      emit = p7P_MSC_CODON(gm_tr, sub_k, c3);
      
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
      scur = snxt;
      k--; i = j - c - 2;
      break;

    case p7T_N:         /* N connects from S, N */
      if (XMX_SP(i, p7G_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = ( (i == 0) ? p7T_S : p7T_N);
      break;

    case p7T_B:         /* B connects from N, J */
      if (XMX_SP(i,p7G_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);
      if      (esl_FCompare_old(XMX_SP(i,p7G_B), XMX_SP(i, p7G_N) + gm_tr->xsc[p7P_N][p7P_MOVE], tol) == eslOK) scur = p7T_N;
      else  ESL_EXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

      default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over statetype[tpos-1] */


    /* Find donor site */
    if (scur == p7T_P) {
      if (PMX_SP(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible P reached at k=%d,i=%d", k,i);

	  sub_i = i_start + i - 1;
	  sub_k = k_start + k - 1;

	  if(sub_dsq[sub_i-2] < p7P_MAXNUC) v = sub_dsq[sub_i-2];
      else                           v = p7P_MAXCODONS1;

	  if(sub_dsq[sub_i-1] < p7P_MAXNUC) w = sub_dsq[sub_i-1];
      else                           w = p7P_MAXCODONS1;

	  if(sub_dsq[sub_i]   < p7P_MAXNUC) x = sub_dsq[sub_i];
      else                           x = p7P_MAXCODONS1;


      c0 = p7P_CODON3_FS1(v, w, x);
      c0 = p7P_MINIDX(c0, p7P_DEGEN1_C);
      emit0 = p7P_MSC_CODON(gm_tr, sub_k, c0);

	  for (j = i-min_intron+1; j > 0; j--)
	  { 
	    sub_d = i_start + j - 1;
        
		if(sub_dsq[sub_d-3]   < p7P_MAXNUC) t = sub_dsq[sub_d-3];
        else                           t = p7P_MAXCODONS1;

		if(sub_dsq[sub_d-2] < p7P_MAXNUC) u = sub_dsq[sub_d-2];
        else                           u = p7P_MAXCODONS1;

        c1 = p7P_CODON3_FS1(u, w, x);
        c1 = p7P_MINIDX(c1, p7P_DEGEN1_C);
        emit1 = p7P_MSC_CODON(gm_tr, sub_k, c1);

        c2 = p7P_CODON3_FS1(t, u, x);
        c2 = p7P_MINIDX(c2, p7P_DEGEN1_C);
        emit2 = p7P_MSC_CODON(gm_tr, sub_k, c2);

        if(SIGNAL(sub_dsq[sub_d-1], sub_dsq[sub_d]) == DONOR_GT) {
          if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(j-2,k-1) + signal_scores[p7S_GTAG] + emit0, tol) == eslOK) { snxt = p7T_M; c = 0; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(j-2,k-1) + signal_scores[p7S_GTAG] + emit0, tol) == eslOK) { snxt = p7T_D; c = 0; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(j-3,k-1) + signal_scores[p7S_GTAG] + emit1, tol) == eslOK) { snxt = p7T_M; c = 1; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(j-3,k-1) + signal_scores[p7S_GTAG] + emit1, tol) == eslOK) { snxt = p7T_D; c = 1; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(j-4,k-1) + signal_scores[p7S_GTAG] + emit2, tol) == eslOK) { snxt = p7T_M; c = 2; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(j-4,k-1) + signal_scores[p7S_GTAG] + emit2, tol) == eslOK) { snxt = p7T_D; c = 2; break; }
        }
        else if (SIGNAL(sub_dsq[sub_d-1], sub_dsq[sub_d]) == DONOR_GC) {
          if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(j-2,k-1) + signal_scores[p7S_GCAG] + emit0, tol) == eslOK) { snxt = p7T_M; c = 0; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(j-2,k-1) + signal_scores[p7S_GCAG] + emit0, tol) == eslOK) { snxt = p7T_D; c = 0; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(j-3,k-1) + signal_scores[p7S_GCAG] + emit1, tol) == eslOK) { snxt = p7T_M; c = 1; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(j-3,k-1) + signal_scores[p7S_GCAG] + emit1, tol) == eslOK) { snxt = p7T_D; c = 1; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(j-4,k-1) + signal_scores[p7S_GCAG] + emit2, tol) == eslOK) { snxt = p7T_M; c = 2; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(j-4,k-1) + signal_scores[p7S_GCAG] + emit2, tol) == eslOK) { snxt = p7T_D; c = 2; break; }
        }
        else if (SIGNAL(sub_dsq[sub_d-1], sub_dsq[sub_d]) == DONOR_AT) {
          if      (esl_FCompare_old(PMX_SP(i,k), MMX_SP(j-2,k-1) + signal_scores[p7S_ATAC] + emit0, tol) == eslOK) { snxt = p7T_M; c = 0; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(j-2,k-1) + signal_scores[p7S_ATAC] + emit0, tol) == eslOK) { snxt = p7T_D; c = 0; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(j-3,k-1) + signal_scores[p7S_ATAC] + emit1, tol) == eslOK) { snxt = p7T_M; c = 1; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(j-3,k-1) + signal_scores[p7S_ATAC] + emit1, tol) == eslOK) { snxt = p7T_D; c = 1; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), MMX_SP(j-4,k-1) + signal_scores[p7S_ATAC] + emit2, tol) == eslOK) { snxt = p7T_M; c = 2; break; }
          else if (esl_FCompare_old(PMX_SP(i,k), DMX_SP(j-4,k-1) + signal_scores[p7S_ATAC] + emit2, tol) == eslOK) { snxt = p7T_D; c = 2; break; }
        }
      }
	  if (j == 0) ESL_EXCEPTION(eslFAIL, "P at k=%d,i=%d couldn't be traced", k,i);
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
/*----------------- end, p7_GViterbi_SplicedTrace() ------*/


/*****************************************************************
 * 5. Benchmark driver.
 *****************************************************************/
#ifdef p7GENERIC_VITERBI_SPLICED_BENCHMARK
/*
   gcc -g -O3 -Wall -std=gnu99 -o generic_viterbi_spliced_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_VITERBI_SPLICED_BENCHMARK generic_viterbi_spliced.c -lhmmer -leasel -lm
   ./generic_viterbi_spliced_benchmark <hmmfile>
 */
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

static ESL_OPTIONS benchmark_options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-N",        eslARG_INT,    "100", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "-I",        eslARG_INT,    "200", NULL, "n>0", NULL,  NULL, NULL, "length of simulated intron (excl. GT..AG signals)", 0 },
  { "-G",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark Global",                0 },
  { "-D",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark ExtendDown",  0 },
  { "-U",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark ExtendUp",    0 },
  { "-T",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also benchmark Trace after each DP",   0 },
  { "--dummy",   eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark dummy",   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char benchmark_usage[]  = "[-options] <hmmfile>";
static char benchmark_banner[] = "benchmark driver for spliced Viterbi DP algorithms";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go           = p7_CreateDefaultApp(benchmark_options, 1, argc, argv, benchmark_banner, benchmark_usage);
  char           *hmmfile      = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w            = esl_stopwatch_Create();
  ESL_RANDOMNESS *r            = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA        = NULL;
  ESL_ALPHABET   *abcDNA       = esl_alphabet_Create(eslDNA);
  P7_HMMFILE     *hfp          = NULL;
  P7_HMM         *hmm          = NULL;
  P7_BG          *bgAA         = NULL;
  P7_PROFILE     *gm           = NULL;
  P7_FS_PROFILE  *gm_tr        = NULL;
  ESL_GENCODE    *gcode        = NULL;
  P7_CODONTABLE  *codon_table  = NULL;
  ESL_SQ         *sq           = NULL;
  SPLICE_PIPELINE *pli         = NULL;
  P7_TRACE        *tr          = NULL;
  int             N            = esl_opt_GetInteger(go, "-N");
  int             I            = esl_opt_GetInteger(go, "-I");
  int             intron_total = I + 4;
  int             do_G         = (esl_opt_GetBoolean(go, "-G") || (!esl_opt_GetBoolean(go, "-D") && !esl_opt_GetBoolean(go, "-U") && !esl_opt_GetBoolean(go, "--dummy")));
  int             do_D         = (esl_opt_GetBoolean(go, "-D") || (!esl_opt_GetBoolean(go, "-G") && !esl_opt_GetBoolean(go, "-U") && !esl_opt_GetBoolean(go, "--dummy")));
  int             do_U         = (esl_opt_GetBoolean(go, "-U") || (!esl_opt_GetBoolean(go, "-G") && !esl_opt_GetBoolean(go, "-D") && !esl_opt_GetBoolean(go, "--dummy")));
  int             do_T         = esl_opt_GetBoolean(go, "-T");
  int             do_dummy     = esl_opt_GetBoolean(go, "--dummy");
  ESL_DSQ        *dsq          = NULL;
  int             i, j, k, L_amino, L_dna_total;
  int64_t         total_cells;
  float           final_C;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)          != eslOK) p7_Fail("Failed to read HMM");

  gcode       = esl_gencode_Create(abcDNA, abcAA);
  bgAA        = p7_bg_Create(abcAA);
  gm          = p7_profile_Create(hmm->M, abcAA);
  gm_tr       = p7_profile_fs_Create(hmm->M, abcAA, 1);
  codon_table = p7_codontable_Create(gcode);
  sq          = esl_sq_CreateDigital(abcAA);

  p7_ProfileConfig   (hmm, bgAA,        gm,    hmm->M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, hmm->M, p7_UNILOCAL);

  pli = p7_splicepipeline_Create(NULL, hmm->M, hmm->M * 3);
  p7_splicescores_GrowTo(pli->splice_scores, hmm->M);
  tr = p7_trace_fs_Create();

  /* Baseline: time to generate sequences alone */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, NULL);
      L_amino     = sq->n;
      L_dna_total = L_amino * 3 + intron_total;
      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) p7_Fail("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;
      j = 1;
      for (k = 1; k <= L_amino / 2; k++) { p7_codontable_GetCodon(codon_table, r, sq->dsq[k], dsq + j); j += 3; }
      dsq[j++] = 2;  /* G */
      dsq[j++] = 3;  /* T */
      for (k = 0; k < I; k++) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = 0;  /* A */
      dsq[j++] = 2;  /* G */
      for (k = L_amino / 2 + 1; k <= L_amino; k++) { p7_codontable_GetCodon(codon_table, r, sq->dsq[k], dsq + j); j += 3; }
    }
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark */
  total_cells = 0;
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, NULL);
      L_amino     = sq->n;
      L_dna_total = L_amino * 3 + intron_total;
      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) p7_Fail("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;
      j = 1;
      for (k = 1; k <= L_amino / 2; k++) { p7_codontable_GetCodon(codon_table, r, sq->dsq[k], dsq + j); j += 3; }
      dsq[j++] = 2;  /* G */
      dsq[j++] = 3;  /* T */
      for (k = 0; k < I; k++) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = 0;  /* A */
      dsq[j++] = 2;  /* G */
      for (k = L_amino / 2 + 1; k <= L_amino; k++) { p7_codontable_GetCodon(codon_table, r, sq->dsq[k], dsq + j); j += 3; }

      p7_fs_ReconfigLength(gm_tr, L_dna_total / 3);
      p7_gmx_GrowTo(pli->vit, hmm->M, L_dna_total, L_dna_total);

	  if (do_G) {
        p7_GViterbi_SplicedGlobal(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, 1, hmm->M, pli->min_intron);
	    final_C = pli->vit->xmx[L_dna_total * p7G_NXCELLS + p7G_C];
        if (final_C != -eslINFINITY && do_T) {
		  p7_GViterbi_SplicedTrace(dsq, gm_tr, pli->vit, pli->splice_scores->signal_scores, tr, 1, L_dna_total, 1, hmm->M, pli->min_intron);
          p7_trace_Reuse(tr);
		}
	  }
      if (do_D) {
	    p7_GViterbi_SplicedExtendDown(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, 1, hmm->M, pli->min_intron);
        if (do_T) { 
		  p7_GViterbi_SplicedTrace(dsq, gm_tr, pli->vit, pli->splice_scores->signal_scores, tr, 1, L_dna_total, 1, hmm->M, pli->min_intron);
          p7_trace_Reuse(tr);
		}
	  }
      if (do_U) { 
	    p7_GViterbi_SplicedExtendUp(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, 1, hmm->M, pli->min_intron);
        if (do_T) {
		  p7_GViterbi_SplicedTrace(dsq, gm_tr, pli->vit, pli->splice_scores->signal_scores, tr, 1, L_dna_total, 1, hmm->M, pli->min_intron);
          p7_trace_Reuse(tr);
        }
	  }
      if (do_dummy) {
        p7_GViterbi_SplicedGlobal_Dummy(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, 1, hmm->M, pli->min_intron);
      } 
      total_cells += (int64_t) L_dna_total * hmm->M;
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) total_cells * 1e-6 / bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M          = %d\n",   hmm->M);
  printf("# N          = %d\n",   N);
  printf("# I          = %d\n",   I);
  printf("# %.1f Mc/s\n", Mcs);

  if (dsq != NULL) free(dsq);
  p7_trace_fs_Destroy(tr);
  p7_splicepipeline_Destroy(pli);
  esl_sq_Destroy(sq);
  p7_codontable_Destroy(codon_table);
  p7_profile_fs_Destroy(gm_tr);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bgAA);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_VITERBI_SPLICED_BENCHMARK*/
/*---------------- end, benchmark driver ----------------*/


/*****************************************************************
 * 6. Unit tests.
 *****************************************************************/
#ifdef p7GENERIC_VITERBI_SPLICED_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_viterbi():
 *
 * A simulated intron (GT + <intron_len> random nucleotides + AG) is inserted
 * at the midpoint of a profile-emitted, reverse-translated DNA sequence,
 * splitting it into two exons.  For each sequence, three tests are run on
 * the full spliced sequence (both exons + intron):
 *
 * 1. Global + Trace: the trace must contain at least one
 *    P state (p7T_P), confirming that the splice junction was detected.
 *
 * 2. ExtendDown: the final C state must be reachable
 *    (not -eslINFINITY).
 *
 * 3. ExtendUp: the final C state must be reachable
 *    (not -eslINFINITY).
 */
static void
utest_viterbi(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
              ESL_GENCODE *gcode, P7_BG *bgAA, P7_CODONTABLE *codon_table,
              int M, int N, int intron_len)
{
  char           *msg         = "splice viterbi unit test failed";
  P7_HMM         *hmm         = NULL;
  P7_PROFILE     *gm          = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_tr       = p7_profile_fs_Create(M, abcAA, 1);
  ESL_SQ         *sq          = esl_sq_CreateDigital(abcAA);
  P7_TRACE       *tr          = p7_trace_fs_Create();
  ESL_DSQ        *dsq         = NULL;
  SPLICE_PIPELINE *pli        = NULL;
  int             intron_total = intron_len + 4;  /* GT + intron_len random nucs + AG */
  int             L_amino, L_dna_total;
  int             i, j, n_p;
  float           final_C;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig   (hmm, bgAA, gm,    M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, M, p7_UNILOCAL);

  pli = p7_splicepipeline_Create(NULL, M, M * 3);

  while (N--)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, NULL);
      L_amino     = sq->n;
      L_dna_total = L_amino * 3 + intron_total;

      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) esl_fatal("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;

      /* Reverse-translate first half of the sequence (exon 1) */
      j = 1;
      for (i = 1; i <= L_amino / 2; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }
      /* Simulated intron: GT + intron_len random nucleotides + AG
       * eslDNA alphabet: A=0, C=1, G=2, T=3                        */
      dsq[j++] = 2;  /* G */
      dsq[j++] = 3;  /* T */
      for (i = 0; i < intron_len; i++) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = 0;  /* A */
      dsq[j++] = 2;  /* G */
      /* Reverse-translate second half of the sequence (exon 2) */
      for (i = L_amino / 2 + 1; i <= L_amino; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }

      p7_fs_ReconfigLength(gm_tr, L_dna_total/3);
      p7_gmx_GrowTo(pli->vit, M, L_dna_total, L_dna_total);
      p7_splicescores_GrowTo(pli->splice_scores, M);

      /* --- Test 1: Global + Trace; trace must contain >= 1 P state --- */
      p7_GViterbi_SplicedGlobal(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, 1, hmm->M, pli->min_intron);
      final_C = pli->vit->xmx[L_dna_total * p7G_NXCELLS + p7G_C];
	  if (final_C != -eslINFINITY) {
        p7_GViterbi_SplicedTrace(dsq, gm_tr, pli->vit, pli->splice_scores->signal_scores, tr, 1, L_dna_total, 1, M, pli->min_intron);
        n_p = 0;
        for (i = 0; i < tr->N; i++)
          if (tr->st[i] == p7T_P) n_p++;
        if (n_p < 1) esl_fatal(msg);
        p7_trace_Reuse(tr);
	  }

      /* --- Test 2: ExtendDown; final C must be reachable --- */
      p7_GViterbi_SplicedExtendDown(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, 1, hmm->M, pli->min_intron);
      final_C = pli->vit->xmx[L_dna_total * p7G_NXCELLS + p7G_C];
      if (final_C == -eslINFINITY) esl_fatal(msg);

      /* --- Test 3: ExtendUp; final C must be reachable --- */
      p7_GViterbi_SplicedExtendUp(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, 1, hmm->M, pli->min_intron);
      final_C = pli->vit->xmx[L_dna_total * p7G_NXCELLS + p7G_C];
      if (final_C == -eslINFINITY) esl_fatal(msg);
    }

  if (dsq != NULL) free(dsq);
  esl_sq_Destroy(sq);
  p7_trace_fs_Destroy(tr);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_tr);
  p7_splicepipeline_Destroy(pli);
}

/* utest_splice():
 *
 * A simulated intron (GT + random nucleotides + AG) is inserted at the 
 * midpoint of a profile-emitted consensus, reverse-translated DNA sequence,
 * splitting it into two exons. Intron is inserted at different postions, 
 * with and without degenerate nucleotides. Test ability to detect correct 
 * splice site.
 *
 */
static void
utest_splice(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
              ESL_GENCODE *gcode, P7_BG *bgAA, P7_CODONTABLE *codon_table)
{
  char           *msg         = "splice viterbi unit test failed";
  int             M           = 200;
  P7_HMM         *hmm         = NULL;
  P7_PROFILE     *gm          = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_tr       = p7_profile_fs_Create(M, abcAA, 1);
  ESL_SQ         *sq          = esl_sq_CreateDigital(abcAA);
  P7_TRACE       *tr          = p7_trace_fs_Create();
  ESL_DSQ        *dsq         = NULL;
  SPLICE_PIPELINE *pli        = NULL;
  int             intron_len  = 500;  
  int             L_amino, L_dna_total;
  int             n, i, j, c;
  int             v, x;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig   (hmm, bgAA, gm,    M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, M, p7_UNILOCAL);

  pli = p7_splicepipeline_Create(NULL, M, M * 3);

  p7_emit_SimpleConsensus(hmm, sq);
  L_amino     = sq->n;
  L_dna_total = L_amino * 3 + intron_len + 4;

  p7_fs_ReconfigLength(gm_tr, L_dna_total/3);
  p7_gmx_GrowTo(pli->vit, M, L_dna_total, L_dna_total);
  p7_splicescores_GrowTo(pli->splice_scores, M);

  for(n = 0; n < 3; n++) {
    if (dsq != NULL) free(dsq);  
    if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) esl_fatal("malloc failed");
    dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;

    j = 1;
	for (i = 1; i <= L_amino / 2; i++) {
      p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
      j += 3;
    }

	v = dsq[j-2];
	x = dsq[j-1];

    j -= n;
  
    dsq[j++] = 2;  /* G */
    dsq[j++] = 3;  /* T */
    for (i = 0; i < intron_len; i++) dsq[j++] = esl_rnd_Roll(r, 4);
    dsq[j++] = 0;  /* A */
    dsq[j++] = 2;  /* G */	

	if(n == 2)
	  dsq[j++] = v;
	if(n >= 1) 
	  dsq[j++] = x;
   
    /* Reverse-translate second half of the sequence (exon 2) */
    for (i = L_amino / 2 + 1; i <= L_amino; i++) {
      p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
      j += 3;
    } 

    p7_GViterbi_SplicedGlobal(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, 1, M, pli->min_intron);
	p7_GViterbi_SplicedTrace(dsq, gm_tr, pli->vit, pli->splice_scores->signal_scores, tr, 1, L_dna_total, 1, M, pli->min_intron);

	c = -1;
	for (i = 0; i < tr->N; i++)
	  if (tr->st[i] == p7T_P) c = tr->c[i];
    
    if( n == 0 && c != 0) esl_fatal(msg);
    if( n == 1 && c != 2) esl_fatal(msg);
    if( n == 2 && c != 1) esl_fatal(msg);

	p7_trace_Reuse(tr);
    p7_GViterbi_SplicedExtendDown(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, 1, M, pli->min_intron);
	p7_GViterbi_SplicedTrace(dsq, gm_tr, pli->vit, pli->splice_scores->signal_scores, tr, 1, L_dna_total, 1, M, pli->min_intron);

	c = -1;
	for (i = 0; i < tr->N; i++)
      if (tr->st[i] == p7T_P) c = tr->c[i];
    
    if( n == 0 && c != 0) esl_fatal(msg);
    if( n == 1 && c != 2) esl_fatal(msg);
    if( n == 2 && c != 1) esl_fatal(msg);
	
	p7_trace_Reuse(tr);
	p7_GViterbi_SplicedExtendUp(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, 1, M, pli->min_intron);
    p7_GViterbi_SplicedTrace(dsq, gm_tr, pli->vit, pli->splice_scores->signal_scores, tr, 1, L_dna_total, 1, M, pli->min_intron);

	c = -1;
	for (i = 0; i < tr->N; i++)
      if (tr->st[i] == p7T_P) c = tr->c[i];
    
    if( n == 0 && c != 0) esl_fatal(msg);
    if( n == 1 && c != 2) esl_fatal(msg);
    if( n == 2 && c != 1) esl_fatal(msg);

	p7_trace_Reuse(tr);
	/* Test handling of degerenate nucleotides */ 
    if (dsq != NULL) free(dsq);  
    if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) esl_fatal("malloc failed");
    dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;

    j = 1;
	for (i = 1; i <= L_amino / 2; i++) {
      p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
      j += 3;
    }

	dsq[j-3] = abcDNA->Kp-3;
	dsq[j-2] = abcDNA->Kp-3;
	dsq[j-1] = abcDNA->Kp-3;

	v = dsq[j-2];
	x = dsq[j-1];

    j -= n;
  
    dsq[j++] = 2;  /* G */
    dsq[j++] = 3;  /* T */
    for (i = 0; i < intron_len; i++) dsq[j++] = esl_rnd_Roll(r, 4);
    dsq[j++] = 0;  /* A */
    dsq[j++] = 2;  /* G */	

	if(n == 2)
	  dsq[j++] = v;
	if(n >= 1) 
	  dsq[j++] = x;
   
    /* Reverse-translate second half of the sequence (exon 2) */
    for (i = L_amino / 2 + 1; i <= L_amino; i++) {
      p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
      j += 3;
    } 
    
    p7_GViterbi_SplicedGlobal(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, 1, M, pli->min_intron);
	p7_GViterbi_SplicedTrace(dsq, gm_tr, pli->vit, pli->splice_scores->signal_scores, tr, 1, L_dna_total, 1, M, pli->min_intron);

	c = -1;
	for (i = 0; i < tr->N; i++)
	  if (tr->st[i] == p7T_P) c = tr->c[i];

    if( n == 0 && c != 0) esl_fatal(msg);
    if( n == 1 && c != 2) esl_fatal(msg);
    if( n == 2 && c != 1) esl_fatal(msg);

	p7_trace_Reuse(tr);
    p7_GViterbi_SplicedExtendDown(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, 1, M, pli->min_intron);
	p7_GViterbi_SplicedTrace(dsq, gm_tr, pli->vit, pli->splice_scores->signal_scores, tr, 1, L_dna_total, 1, M, pli->min_intron);

	c = -1;
	for (i = 0; i < tr->N; i++)
      if (tr->st[i] == p7T_P) c = tr->c[i];

    if( n == 0 && c != 0) esl_fatal(msg);
    if( n == 1 && c != 2) esl_fatal(msg);
    if( n == 2 && c != 1) esl_fatal(msg);  
	
	p7_trace_Reuse(tr);
	p7_GViterbi_SplicedExtendUp(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, 1, M, pli->min_intron);
    p7_GViterbi_SplicedTrace(dsq, gm_tr, pli->vit, pli->splice_scores->signal_scores, tr, 1, L_dna_total, 1, M, pli->min_intron);

	c = -1;
	for (i = 0; i < tr->N; i++)
      if (tr->st[i] == p7T_P) c = tr->c[i];

    if( n == 0 && c != 0) esl_fatal(msg);
    if( n == 1 && c != 2) esl_fatal(msg);
    if( n == 2 && c != 1) esl_fatal(msg);

	p7_trace_Reuse(tr);
  }
  
  if (dsq != NULL) free(dsq);
  esl_sq_Destroy(sq);
  p7_trace_fs_Destroy(tr);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_tr);
  p7_splicepipeline_Destroy(pli);
}

#endif /*p7GENERIC_VITERBI_SPLICED_TESTDRIVE*/
/*---------------------- end, unit tests ------------------------*/


/*****************************************************************
 * 7. Test driver.
 *****************************************************************/
#ifdef p7GENERIC_VITERBI_SPLICED_TESTDRIVE
/*
   gcc -g -Wall -std=gnu99 -o generic_viterbi_spliced_utest -I. -L. -I../easel -L../easel -Dp7GENERIC_VITERBI_SPLICED_TESTDRIVE generic_viterbi_spliced.c -lhmmer -leasel -lm
   ./generic_viterbi_spliced_utest
 */
#include "esl_gencode.h"
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default   env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42",  NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M",        eslARG_INT,     "100", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,     "50",  NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  { "-I",        eslARG_INT,    "500",  NULL, "n>0", NULL,  NULL, NULL, "simulated intron length (random nucleotides between GT..AG)", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for spliced Viterbi DP algorithms";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA  = esl_alphabet_Create(eslAMINO);
  ESL_ALPHABET   *abcDNA = esl_alphabet_Create(eslDNA);
  P7_BG          *bgAA   = p7_bg_Create(abcAA);
  ESL_GENCODE    *gcode  = esl_gencode_Create(abcDNA, abcAA);
  P7_CODONTABLE  *ct     = p7_codontable_Create(gcode);
  int             M      = esl_opt_GetInteger(go, "-M");
  int             N      = esl_opt_GetInteger(go, "-N");
  int             I      = esl_opt_GetInteger(go, "-I");

  utest_viterbi(r, abcAA, abcDNA, gcode, bgAA, ct, M, N, I);
  utest_splice(r, abcAA, abcDNA, gcode, bgAA, ct);

  esl_alphabet_Destroy(abcAA);
  esl_alphabet_Destroy(abcDNA);
  p7_bg_Destroy(bgAA);
  esl_gencode_Destroy(gcode);
  p7_codontable_Destroy(ct);
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  printf("All tests passed.\n");
  return eslOK;
}
#endif /*p7GENERIC_VITERBI_SPLICED_TESTDRIVE*/
/*--------------------- end, test driver ------------------------*/


