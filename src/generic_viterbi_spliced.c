/* Spliced Viterbi algorithms and traceback.
 *
 * Produce a spliced alignment with splice site detection via the P state
 *
 * Contents:
 *    1. p7_GViterbi_Spliced()
 *       p7_GViterbi_SplicedChk()
 *    2. p7_GViterbi_SplicedTrace()
 *       p7_GViterbi_SplicedTraceReplay()
 *       p7_GViterbi_SplicedTraceChk()
 *    3. Benchmark driver.
 *    4. Unit tests.
 *    5. Test driver.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_splice.h"

#define PVX4(i,k) (acc_ivx[((k)*SPLICE_ROWS) + (i)])

/* nuc_at(): nucleotide at row <row> of the sub sequence starting at
 * <i_start>, with anything other than A,C,G,T as p7P_MAXCODONS1. */
static inline int
nuc_at(const ESL_DSQ *sub_dsq, int i_start, int row)
{
  ESL_DSQ x = sub_dsq[i_start + row - 1];
  return (x < p7P_MAXNUC) ? x : p7P_MAXCODONS1;
}

/* donor_signal(): splice signal index (p7S_GTAG, p7S_GCAG or p7S_ATAC)
 * of a donor dinucleotide, or -1. */
static inline int
donor_signal(int nuc1, int nuc2)
{
  switch (SIGNAL(nuc1, nuc2)) {
    case DONOR_GT: return p7S_GTAG;
    case DONOR_GC: return p7S_GCAG;
    case DONOR_AT: return p7S_ATAC;
    default:       return -1;
  }
}

/* acceptor_signal(): ACCEPT_AG or ACCEPT_AC for an acceptor
 * dinucleotide, or -1. */
static inline int
acceptor_signal(int nuc1, int nuc2)
{
  switch (SIGNAL(nuc1, nuc2)) {
    case ACCEPT_AG: return ACCEPT_AG;
    case ACCEPT_AC: return ACCEPT_AC;
    default:        return -1;
  }
}

/* acceptor_at(): acceptor signal (ACCEPT_AG or ACCEPT_AC) that the fill
 * records at row <i>, from rows i-2 and i-1, or -1. */
static inline int
acceptor_at(const ESL_DSQ *sub_dsq, int i_start, int i)
{
  if (i < 3) return -1;
  return acceptor_signal(nuc_at(sub_dsq, i_start, i-2), nuc_at(sub_dsq, i_start, i-1));
}

/* donor_at(): splice signal of the donor that the fill's donor update
 * at row <iu> handles as don2 (the dinucleotide starting at row
 * <iu>-<min_intron>), or -1. The fill updates donors from row
 * <min_intron>+3 on.
 */
static inline int
donor_at(const ESL_DSQ *sub_dsq, int i_start, int min_intron, int iu)
{
  if (iu < min_intron + 3) return -1;
  return donor_signal(nuc_at(sub_dsq, i_start, iu - min_intron), nuc_at(sub_dsq, i_start, iu - min_intron + 1));
}

/* Function:  p7_GViterbi_Spliced()
 * Synopsis:  Translated spliced Viterbi algorithm
 *
 * Purpose:   For finding the maximum scoring splice site between two
 *            or more exons. Aligns from position <i_start> to <i_end>
 *            on the <sub_dsq> and from <k_start> to <k_end> on the
 *            <gm_tr>, in either fully global, semi-global (one end only)
 *            or local mode. The DP matrix <gx> must include room for the
 *            standard core model states <M, I, D>. The splice <P> state
 *            is only stored temporarily in the P7_IVX buffer. The
 *            <P> state acts as a modified <M> state, emitting a codon
 *            that is made of either two nucleotides from before the
 *            donor site and one from after the acceptor site <C2>, one
 *            nucleotide from before the donor site and two from
 *            after the acceptor <C1>, or from the three after the
 *            acceptor <C0>.
 *
 *            Potential donor site scores are recorded in the <P_scores>
 *            matrix (33*M) via the SSX macro. Splice signal scores are
 *            stored in <signal_scores> array and <acceptor_..> and
 *            <donor_..> arrays return 0 for valid sites and -inf for 
 *            all others. 
 *
 * Args:      sub_dsq       - nucleotide sequence 
 *            gm_tr         - a codon profile.
 *            gx            - DP matrix with room for an MxL alignment
 *            acc_iv        - a buffer for the P state at the acceptor site
 *            don_iv        - a buffer for the P state at donor site
 *            signal_scores - array of splice site signal scores
 *            i_start       - start position on the <sub_dsq>
 *            i_end         - end position on the <sub_dsq>
 *            k_start       - start position on the <gm_tr>
 *            k_end         - end position on the <gm_tr>
 *            min_intron    - minimum intron length
 *            global_start  - bool value controlling global vs local entry
 *            global_end    - bool value controlling global vs local exit
 *
 * Return:    <eslOK> on success.
 */
int
p7_GViterbi_Spliced(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, P7_GMX *gx, P7_IVX *acc_iv, P7_IVX *don_iv, float *signal_scores, int i_start, int i_end, int k_start, int k_end, int min_intron, int global_start, int global_end)
{
  
  float const *tsc  = gm_tr->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float       *acc_ivx  = acc_iv->ivx;
  float       *don_ivx  = don_iv->ivx;
  float        entry = (global_start ? -eslINFINITY : 0.f);
  float        exit  = (global_end   ? -eslINFINITY : 0.f);
  int          L = (i_end - i_start + 1);
  int          M = (k_end - k_start + 1);
  int          C0;
  int          C1[5] = {0, 0, 0, 0, 0};
  int          C2[5] = {0, 0, 0, 0, 0};
  int          i,k;
  int          r, s, t, u;
  int          v, w, x;
  int          nuc1, nuc3;
  int          loop_end;
  int          sub_k;
  int          pv_i, pv_pi;
  int          acc0, acc1, acc2;
  int          don0, don1, don2;
  int          don_sig;
  float        TMP_SC;
  float const *rsc_c0;
  float const *rsc_c1[5];
  
  if(gm_tr->codon_lengths != 1) ESL_EXCEPTION(eslEINVAL, "proflie not allocated for 1 codon length");

  /* Note on variable names 
   * ..0 are for splice codons where zero nucleotides come from before the donor,
   * ..1 are for splice codons where one nucleotide comes from before the donor,
   * ..2 are for splice codons where two nucleotides come from before the donor,
   */

  acc0 = acc1 = acc2 = -1;
  don0 = don1 = don2 = -1;

  /* Reset donor site score buffer */
  for(k = 0; k < M; k++) {
    for(i = 0; i < SIGNAL_MEM_SIZE; i++)
     don_ivx[k*SIGNAL_MEM_SIZE+i] = -eslINFINITY;
  }

  /* Reset acceptor site score buffer */
  for (i = 0; i < SPLICE_ROWS; i++) {
    for(k = 0; k <= M; k++) {
      acc_ivx[k*SPLICE_ROWS+i] = -eslINFINITY;
    }
  }

  /* Initialization of the zero row.  */
  XMX(0,p7G_N) = 0.;                                     /* S->N, p=1            */
  XMX(0,p7G_B) = gm_tr->xsc[p7P_N][p7P_MOVE];            /* S->N->B, no N-tail   */
  XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX(0,k) = DMX(0,k) = IMX(0,k) = -eslINFINITY;

  /*Special cases for the first 2 rows (no codons yet)*/
  v = w = x = -1;
  for(i = 1; i <= 2; i++)
  {
    w = x;
    /* if new nucleotide is not A,C,G, or T set it to placeholder value */
    x = nuc_at(sub_dsq, i_start, i);

    XMX(i,p7G_N) =  0;
    XMX(i,p7G_B) =  gm_tr->xsc[p7P_N][p7P_MOVE];

    for (k = 0; k <= M; k++)
      MMX(i,k) = DMX(i,k) = IMX(i,k) = -eslINFINITY;

    XMX(i,p7G_E) = XMX(i,p7G_C) = XMX(i,p7G_J) = -eslINFINITY;
  }

  /*Special cases for the rows 3-min_intron+2 (no donor site look backs)*/
  loop_end = ESL_MIN(L, min_intron+2);
  for (i = 3; i <= loop_end; i++)
  {
    v = w;
    w = x;
    x = nuc_at(sub_dsq, i_start, i);

    C0 = p7P_CODON3_FS1(v, w, x);
    C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);
    rsc_c0 = gm_tr->rsc[C0];

    /* get acceptor site */
    acc0 = acc1;
    acc1 = acc2;
    acc2 = acceptor_signal(v, w);
    
	if(!global_start) {
      XMX(i,p7G_N) = XMX(i-3,p7G_N) + gm_tr->xsc[p7P_N][p7P_LOOP];
      XMX(i,p7G_B) = XMX(i,p7G_N)   + gm_tr->xsc[p7P_N][p7P_MOVE];
    }
	else XMX(i,p7G_N) = XMX(i,p7G_B) = -eslINFINITY;

    MMX(i,0) = IMX(i,0) = DMX(i,0) = -eslINFINITY;

    if(global_start && i != 3)
      MMX(i,1) = -eslINFINITY;
    else
      MMX(i,1) = XMX(i-3,p7G_B) + rsc_c0[k_start];

    IMX(i,1) = ESL_MAX(MMX(i-3,1) + TSC(p7P_MI,k_start),
                       IMX(i-3,1) + TSC(p7P_II,k_start));

    /* No I state for stop codons */
    if(rsc_c0[k_start] == -eslINFINITY) IMX(i,1) = -eslINFINITY;

    DMX(i,1) = -eslINFINITY;

    XMX(i,p7G_E) = MMX(i,1);

    for (k = 2; k < M; k++) {
      sub_k = k_start + k -1;

      MMX(i,k) = ESL_MAX(MMX(i-3,k-1)   + TSC(p7P_MM,sub_k-1),
                 ESL_MAX(IMX(i-3,k-1)   + TSC(p7P_IM,sub_k-1),
                 ESL_MAX(DMX(i-3,k-1)   + TSC(p7P_DM,sub_k-1),
                         XMX(i-3,p7G_B) + entry))) + rsc_c0[sub_k]; // no B->M transitions for exons. 
       
	  IMX(i,k) = ESL_MAX(MMX(i-3,k) + TSC(p7P_MI,sub_k),
                         IMX(i-3,k) + TSC(p7P_II,sub_k));

	  /* No I state for stop codons */
      if(rsc_c0[sub_k] == -eslINFINITY) IMX(i,k) = -eslINFINITY;

      DMX(i,k) = ESL_MAX(MMX(i,k-1) + TSC(p7P_MD,sub_k-1),
                         DMX(i,k-1) + TSC(p7P_DD,sub_k-1));
 
      XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E),
                     ESL_MAX(MMX(i,k),
                             DMX(i,k))) + exit;
                             
    }

    sub_k = k_start + M -1;
    
    MMX(i,M) = ESL_MAX(MMX(i-3,M-1)   + TSC(p7P_MM,sub_k-1),
               ESL_MAX(IMX(i-3,M-1)   + TSC(p7P_IM,sub_k-1),
               ESL_MAX(DMX(i-3,M-1)   + TSC(p7P_DM,sub_k-1),
                       XMX(i-3,p7G_B) + entry))) + rsc_c0[sub_k];

    IMX(i,M) = -eslINFINITY;

    DMX(i,M) = ESL_MAX(MMX(i,M-1) + TSC(p7P_MD,sub_k-1),
                       DMX(i,M-1) + TSC(p7P_DD,sub_k-1));

    if(!global_end) {
      XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E),
                     ESL_MAX(MMX(i,M),
                             DMX(i,M)));

     XMX(i,p7G_C) = ESL_MAX(XMX(i-3,p7G_C) + gm_tr->xsc[p7P_C][p7P_LOOP],
                            XMX(i,p7G_E)   + gm_tr->xsc[p7P_E][p7P_MOVE]);
   }
   else XMX(i,p7G_C) = -eslINFINITY;

   XMX(i,p7G_J) = -eslINFINITY;

  }

  /* Example Indexing for r, s, t, u (row n is sub_dsq[i_start+n-1]):
   * i_start = 5;
   * min_intron = 13;
   * s = row 1;               // sub_dsq[5];
   * t = row 2;               // sub_dsq[6];
   * u = row 3;               // sub_dsq[7];
   * ...
   * i = min_intron+3;        // 16;
   * r = s;                   // sub_dsq[5];
   * s = t;                   // sub_dsq[6];
   * t = u;                   // sub_dsq[7];
   * u = row i-min_intron+1;  // sub_dsq[8];
   */

  s = nuc_at(sub_dsq, i_start, 1);

  t = nuc_at(sub_dsq, i_start, 2);

  u = nuc_at(sub_dsq, i_start, 3);

  /* Main DP recursion */
  for (i = min_intron+3; i <= L; i++)
  {
    /* get nucleotides and codon */
    pv_i  =  i    % SPLICE_ROWS;
    pv_pi = (i-3) % SPLICE_ROWS; 

    r = s;
    s = t;
    t = u;
    u = nuc_at(sub_dsq, i_start, i-min_intron+1);

    v = w;
    w = x;

    x = nuc_at(sub_dsq, i_start, i);

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
    acc2 = acceptor_signal(v, w);

    /* get donor sites */
    don0 = don1;
    don1 = don2;
    don2 = donor_signal(t, u);

	if(!global_start) {
      XMX(i,p7G_N) = XMX(i-3,p7G_N) + gm_tr->xsc[p7P_N][p7P_LOOP];
      XMX(i,p7G_B) = XMX(i,p7G_N)   + gm_tr->xsc[p7P_N][p7P_MOVE];
    }
	else XMX(i,p7G_N) = XMX(i,p7G_B) = -eslINFINITY;

    MMX(i,0) = IMX(i,0) = DMX(i,0) = -eslINFINITY;

    XMX(i,p7G_E) = -eslINFINITY;

    for (k = 1; k < M; k++) {

      sub_k = k_start + k -1;

      PVX4(pv_i,k) = -eslINFINITY;
      if (acc0 >= 0 || acc1 >= 0 || acc2 >= 0) {
        
        if (acc0 == ACCEPT_AG) {
          TMP_SC = ESL_MAX(SSX0(k, p7S_GTAG) + signal_scores[p7S_GTAG],
                           SSX0(k, p7S_GCAG) + signal_scores[p7S_GCAG]) + rsc_c0[sub_k];
          PVX4(pv_i,k) = ESL_MAX(PVX4(pv_i,k), TMP_SC);
        } else if (acc0 == ACCEPT_AC) {
          TMP_SC = SSX0(k, p7S_ATAC) + signal_scores[p7S_ATAC] + rsc_c0[sub_k];
          PVX4(pv_i,k) = ESL_MAX(PVX4(pv_i,k), TMP_SC);
        }

        if (acc1 == ACCEPT_AG) {
          for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
            TMP_SC = ESL_MAX(SSX1(k, p7S_GTAG, nuc1) + signal_scores[p7S_GTAG],
                             SSX1(k, p7S_GCAG, nuc1) + signal_scores[p7S_GCAG]) + rsc_c1[nuc1][sub_k];

            PVX4(pv_i,k) = ESL_MAX(PVX4(pv_i,k), TMP_SC);
 
          }
        } else if (acc1 == ACCEPT_AC) {
          for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
            TMP_SC = SSX1(k, p7S_ATAC, nuc1) + signal_scores[p7S_ATAC] + rsc_c1[nuc1][sub_k];
            PVX4(pv_i,k) = ESL_MAX(PVX4(pv_i,k), TMP_SC);
          }
        }

        if (acc2 == ACCEPT_AG) {
          TMP_SC = ESL_MAX(SSX2(k, p7S_GTAG, nuc3) + signal_scores[p7S_GTAG],
                           SSX2(k, p7S_GCAG, nuc3) + signal_scores[p7S_GCAG]);
          PVX4(pv_i,k) = ESL_MAX(PVX4(pv_i,k), TMP_SC);

        } else if (acc2 == ACCEPT_AC) {
          TMP_SC = SSX2(k, p7S_ATAC, nuc3) + signal_scores[p7S_ATAC];
          PVX4(pv_i,k) = ESL_MAX(PVX4(pv_i,k), TMP_SC);
        }
      }

      MMX(i,k) = ESL_MAX(MMX(i-3,k-1)    + TSC(p7P_MM,sub_k-1),
                 ESL_MAX(IMX(i-3,k-1)    + TSC(p7P_IM,sub_k-1),
                 ESL_MAX(DMX(i-3,k-1)    + TSC(p7P_DM,sub_k-1),
                 ESL_MAX(XMX(i-3,p7G_B)  + entry,        
                         PVX4(pv_pi,k-1) + TSC_P)))) + rsc_c0[sub_k];

	  IMX(i,k) = ESL_MAX(MMX(i-3,k) + TSC(p7P_MI,sub_k),
                         IMX(i-3,k) + TSC(p7P_II,sub_k));

	  /* No I state for stop codons */
	  if(rsc_c0[sub_k] == -eslINFINITY) IMX(i,k) = -eslINFINITY;

      DMX(i,k) = ESL_MAX(MMX(i,k-1) + TSC(p7P_MD,sub_k-1),
                            DMX(i,k-1) + TSC(p7P_DD,sub_k-1));

      XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E),
                     ESL_MAX(MMX(i,k),
                             DMX(i,k))) + exit;
    }
	  
    sub_k = k_start + M -1;
    MMX(i,M) = ESL_MAX(MMX(i-3,M-1)    + TSC(p7P_MM,sub_k-1),
               ESL_MAX(IMX(i-3,M-1)    + TSC(p7P_IM,sub_k-1),
               ESL_MAX(DMX(i-3,M-1)    + TSC(p7P_DM,sub_k-1),
               ESL_MAX(XMX(i-3,p7G_B)  + entry,
                       PVX4(pv_pi,M-1) + TSC_P))))+ rsc_c0[sub_k];
                           
    IMX(i,M) = -eslINFINITY;

    DMX(i,M) = ESL_MAX(MMX(i,M-1) + TSC(p7P_MD,sub_k-1),
                       DMX(i,M-1) + TSC(p7P_DD,sub_k-1));

	if(!global_end) {
      XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E),
                     ESL_MAX(MMX(i,M),
                             DMX(i,M)));

      XMX(i,p7G_C) = ESL_MAX(XMX(i-3,p7G_C) + gm_tr->xsc[p7P_C][p7P_LOOP],
                             XMX(i,p7G_E)   + gm_tr->xsc[p7P_E][p7P_MOVE]);
    }
	else XMX(i,p7G_C) = -eslINFINITY;
	   
    XMX(i,p7G_J) = -eslINFINITY;

   	/* Separate k loop for donor sites */
	if (don2 >= 0) {
	 
      /* Get codon indices for all C2 type split codons */
      for(nuc3 = 0; nuc3 < p7P_MAXNUC; nuc3++)
        C2[nuc3] = p7P_MINIDX(p7P_CODON3_FS1(r, s, nuc3), p7P_DEGEN1_C);
      C2[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(r, s, p7P_MAXCODONS1), p7P_DEGEN1_C);

      don_sig = don2;
      for (k = 2; k < M; k++) {
        sub_k  = k_start + k - 1;
        TMP_SC = ESL_MAX(MMX(i-min_intron-3,k-1), DMX(i-min_intron-3,k-1));
        for(nuc3 = 0; nuc3 <= p7P_MAXNUC; nuc3++) 
          SSX2(k,don_sig,nuc3) = ESL_MAX(SSX2(k,don_sig,nuc3), TMP_SC + gm_tr->rsc[C2[nuc3]][sub_k]);
      }
    }

	if(don1  >= 0) {
	  nuc1 = ESL_MIN(r, p7P_MAXNUC); 
	  don_sig = don1;
	  for (k = 2; k < M; k++) {
	    sub_k  = k_start + k - 1;
		TMP_SC = ESL_MAX(MMX(i-min_intron-3,k-1), DMX(i-min_intron-3,k-1));
		SSX1(k, don_sig, nuc1) = ESL_MAX(SSX1(k, don_sig, nuc1), TMP_SC);
	  }
	}

	if (don0 >= 0) {
      don_sig = don0;
      for (k = 2; k < M; k++) {
        sub_k  = k_start + k - 1;
        TMP_SC = ESL_MAX(MMX(i-min_intron-3,k-1), DMX(i-min_intron-3,k-1));
        SSX0(k, don_sig) = ESL_MAX(SSX0(k, don_sig), TMP_SC);
	  }
	}

  } // end loop over L
   
  if(global_end) { 
    XMX(L,p7G_E) = ESL_MAX(MMX(L,M), DMX(L,M));
    XMX(L,p7G_C) = XMX(L,p7G_E) + gm_tr->xsc[p7P_E][p7P_MOVE];
  }

  gx->M = M;
  gx->L = L;

  return eslOK;
}



/*****************************************************************
 * Checkpointed spliced Viterbi: p7_GViterbi_SplicedChk() and
 * p7_GViterbi_SplicedTraceChk()
 *
 * The fill computes the same rows as p7_GViterbi_Spliced(), but into a
 * P7_CHKGMX that holds one block of rows at a time. At the start of row
 * b-3 for each block start b, the donor scores (don_iv) and the rows
 * they came from (don_idx, the working ints of the P7_CHKIVX) are saved
 * as checkpoint b/blockR - 1; at row b the P7_CHKGMX saves the saveR
 * rows before b. To compute block b again, the trace restores both,
 * recomputes the acceptor P values and donor updates of rows b-3..b-1
 * (their M, I, D rows are among the saved rows), and then computes the
 * rows of the block in full. The same row functions do both, so the
 * rows come out identical.
 *
 * saveR must be at least min_intron+6: row i reads row i-3, and its
 * donor update reads row i-min_intron-3.
 *****************************************************************/

/* Arguments shared by the checkpointed fill, the trace and the row
 * functions they both use. */
typedef struct {
  const ESL_DSQ       *sub_dsq;
  const P7_FS_PROFILE *gm_tr;
  P7_CHKGMX           *cx;
  P7_CHKIVX           *cv;
  P7_IVX              *acc_iv;
  P7_IVX              *don_iv;
  const float         *signal_scores;
  int                  i_start;
  int                  k_start;
  int                  M;
  int                  L;
  int                  min_intron;
  int                  global_start;
  int                  global_end;
} SPLICED_CHK;

/* spliced_init():
 * Clear the donor scores, donor rows and acceptor P values, and set rows
 * 0..2, as at the start of p7_GViterbi_Spliced(). Block 0 must be loaded.
 */
static void
spliced_init(const SPLICED_CHK *sc)
{
  const P7_FS_PROFILE *gm_tr   = sc->gm_tr;
  float              **dp      = sc->cx->dp;
  float               *xmx     = sc->cx->xmx;
  float               *acc_ivx = sc->acc_iv->ivx;
  float               *don_ivx = sc->don_iv->ivx;
  int                 *don_idx = sc->cv->ivx;
  int                  M       = sc->M;
  int                  i, k;

  for (k = 0; k < M; k++)
    for (i = 0; i < SIGNAL_MEM_SIZE; i++) {
      don_ivx[k*SIGNAL_MEM_SIZE+i] = -eslINFINITY;
      don_idx[k*SIGNAL_MEM_SIZE+i] = -1;
    }

  for (i = 0; i < SPLICE_ROWS; i++)
    for (k = 0; k <= M; k++)
      acc_ivx[k*SPLICE_ROWS+i] = -eslINFINITY;

  XMX(0,p7G_N) = 0.;
  XMX(0,p7G_B) = gm_tr->xsc[p7P_N][p7P_MOVE];
  XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX(0,k) = DMX(0,k) = IMX(0,k) = -eslINFINITY;

  for (i = 1; i <= ESL_MIN(2, sc->L); i++) {
    XMX(i,p7G_N) =  0;
    XMX(i,p7G_B) =  gm_tr->xsc[p7P_N][p7P_MOVE];
    for (k = 0; k <= M; k++)
      MMX(i,k) = DMX(i,k) = IMX(i,k) = -eslINFINITY;
    XMX(i,p7G_E) = XMX(i,p7G_C) = XMX(i,p7G_J) = -eslINFINITY;
  }
}

/* spliced_acc_row():
 * P state values of row <i> (acceptor sites), from the donor scores as
 * they are at the start of row <i>. Rows before min_intron+3 have none.
 */
static void
spliced_acc_row(const SPLICED_CHK *sc, int i)
{
  const P7_FS_PROFILE *gm_tr         = sc->gm_tr;
  const float         *signal_scores = sc->signal_scores;
  float               *acc_ivx       = sc->acc_iv->ivx;
  float               *don_ivx       = sc->don_iv->ivx;
  int                  M             = sc->M;
  int                  pv_i          = i % SPLICE_ROWS;
  int                  acc0, acc1, acc2;
  int                  v, w, x;
  int                  C0;
  int                  C1[p7P_MAXNUC+1];
  int                  nuc1, nuc3;
  int                  k, sub_k;
  float                TMP_SC;
  float const         *rsc_c0;
  float const         *rsc_c1[p7P_MAXNUC+1];

  for (k = 1; k < M; k++) PVX4(pv_i,k) = -eslINFINITY;
  if (i < sc->min_intron + 3) return;

  acc2 = acceptor_at(sc->sub_dsq, sc->i_start, i);
  acc1 = acceptor_at(sc->sub_dsq, sc->i_start, i-1);
  acc0 = acceptor_at(sc->sub_dsq, sc->i_start, i-2);
  if (acc0 < 0 && acc1 < 0 && acc2 < 0) return;

  v = nuc_at(sc->sub_dsq, sc->i_start, i-2);
  w = nuc_at(sc->sub_dsq, sc->i_start, i-1);
  x = nuc_at(sc->sub_dsq, sc->i_start, i);

  C0     = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);
  rsc_c0 = gm_tr->rsc[C0];
  for (nuc1 = 0; nuc1 < p7P_MAXNUC; nuc1++)
    C1[nuc1] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, w, x), p7P_DEGEN1_C);
  C1[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(p7P_MAXCODONS1, w, x), p7P_DEGEN1_C);
  for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++)
    rsc_c1[nuc1] = gm_tr->rsc[C1[nuc1]];
  nuc3 = ESL_MIN(x, p7P_MAXNUC);

  for (k = 1; k < M; k++) {
    sub_k = sc->k_start + k - 1;

    if (acc0 == ACCEPT_AG) {
      TMP_SC = ESL_MAX(SSX0(k, p7S_GTAG) + signal_scores[p7S_GTAG],
                       SSX0(k, p7S_GCAG) + signal_scores[p7S_GCAG]) + rsc_c0[sub_k];
      PVX4(pv_i,k) = ESL_MAX(PVX4(pv_i,k), TMP_SC);
    } else if (acc0 == ACCEPT_AC) {
      TMP_SC = SSX0(k, p7S_ATAC) + signal_scores[p7S_ATAC] + rsc_c0[sub_k];
      PVX4(pv_i,k) = ESL_MAX(PVX4(pv_i,k), TMP_SC);
    }

    if (acc1 == ACCEPT_AG) {
      for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
        TMP_SC = ESL_MAX(SSX1(k, p7S_GTAG, nuc1) + signal_scores[p7S_GTAG],
                         SSX1(k, p7S_GCAG, nuc1) + signal_scores[p7S_GCAG]) + rsc_c1[nuc1][sub_k];
        PVX4(pv_i,k) = ESL_MAX(PVX4(pv_i,k), TMP_SC);
      }
    } else if (acc1 == ACCEPT_AC) {
      for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
        TMP_SC = SSX1(k, p7S_ATAC, nuc1) + signal_scores[p7S_ATAC] + rsc_c1[nuc1][sub_k];
        PVX4(pv_i,k) = ESL_MAX(PVX4(pv_i,k), TMP_SC);
      }
    }

    if (acc2 == ACCEPT_AG) {
      TMP_SC = ESL_MAX(SSX2(k, p7S_GTAG, nuc3) + signal_scores[p7S_GTAG],
                       SSX2(k, p7S_GCAG, nuc3) + signal_scores[p7S_GCAG]);
      PVX4(pv_i,k) = ESL_MAX(PVX4(pv_i,k), TMP_SC);
    } else if (acc2 == ACCEPT_AC) {
      TMP_SC = SSX2(k, p7S_ATAC, nuc3) + signal_scores[p7S_ATAC];
      PVX4(pv_i,k) = ESL_MAX(PVX4(pv_i,k), TMP_SC);
    }
  }
}

/* spliced_mdi_row():
 * M, I, D and special states of row <i> (i >= 3), as in
 * p7_GViterbi_Spliced(). Rows before min_intron+3 have no P state.
 */
static void
spliced_mdi_row(const SPLICED_CHK *sc, int i)
{
  const P7_FS_PROFILE *gm_tr   = sc->gm_tr;
  float const         *tsc     = gm_tr->tsc;
  float              **dp      = sc->cx->dp;
  float               *xmx     = sc->cx->xmx;
  float               *acc_ivx = sc->acc_iv->ivx;
  float                entry   = (sc->global_start ? -eslINFINITY : 0.f);
  float                exit    = (sc->global_end   ? -eslINFINITY : 0.f);
  int                  M       = sc->M;
  int                  k_start = sc->k_start;
  int                  pv_pi   = (i-3) % SPLICE_ROWS;
  int                  v, w, x;
  int                  C0;
  int                  k, sub_k;
  float const         *rsc_c0;

  v = nuc_at(sc->sub_dsq, sc->i_start, i-2);
  w = nuc_at(sc->sub_dsq, sc->i_start, i-1);
  x = nuc_at(sc->sub_dsq, sc->i_start, i);
  C0     = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);
  rsc_c0 = gm_tr->rsc[C0];

  if (!sc->global_start) {
    XMX(i,p7G_N) = XMX(i-3,p7G_N) + gm_tr->xsc[p7P_N][p7P_LOOP];
    XMX(i,p7G_B) = XMX(i,p7G_N)   + gm_tr->xsc[p7P_N][p7P_MOVE];
  }
  else XMX(i,p7G_N) = XMX(i,p7G_B) = -eslINFINITY;

  MMX(i,0) = IMX(i,0) = DMX(i,0) = -eslINFINITY;

  if (i <= sc->min_intron + 2)
    {
      if (sc->global_start && i != 3)
        MMX(i,1) = -eslINFINITY;
      else
        MMX(i,1) = XMX(i-3,p7G_B) + rsc_c0[k_start];

      IMX(i,1) = ESL_MAX(MMX(i-3,1) + TSC(p7P_MI,k_start),
                         IMX(i-3,1) + TSC(p7P_II,k_start));
      if (rsc_c0[k_start] == -eslINFINITY) IMX(i,1) = -eslINFINITY;

      DMX(i,1) = -eslINFINITY;

      XMX(i,p7G_E) = MMX(i,1);

      for (k = 2; k < M; k++) {
        sub_k = k_start + k - 1;

        MMX(i,k) = ESL_MAX(MMX(i-3,k-1)   + TSC(p7P_MM,sub_k-1),
                   ESL_MAX(IMX(i-3,k-1)   + TSC(p7P_IM,sub_k-1),
                   ESL_MAX(DMX(i-3,k-1)   + TSC(p7P_DM,sub_k-1),
                           XMX(i-3,p7G_B) + entry))) + rsc_c0[sub_k];

        IMX(i,k) = ESL_MAX(MMX(i-3,k) + TSC(p7P_MI,sub_k),
                           IMX(i-3,k) + TSC(p7P_II,sub_k));
        if (rsc_c0[sub_k] == -eslINFINITY) IMX(i,k) = -eslINFINITY;

        DMX(i,k) = ESL_MAX(MMX(i,k-1) + TSC(p7P_MD,sub_k-1),
                           DMX(i,k-1) + TSC(p7P_DD,sub_k-1));

        XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E),
                       ESL_MAX(MMX(i,k),
                               DMX(i,k))) + exit;
      }

      sub_k = k_start + M - 1;
      MMX(i,M) = ESL_MAX(MMX(i-3,M-1)   + TSC(p7P_MM,sub_k-1),
                 ESL_MAX(IMX(i-3,M-1)   + TSC(p7P_IM,sub_k-1),
                 ESL_MAX(DMX(i-3,M-1)   + TSC(p7P_DM,sub_k-1),
                         XMX(i-3,p7G_B) + entry))) + rsc_c0[sub_k];
    }
  else
    {
      XMX(i,p7G_E) = -eslINFINITY;

      for (k = 1; k < M; k++) {
        sub_k = k_start + k - 1;

        MMX(i,k) = ESL_MAX(MMX(i-3,k-1)    + TSC(p7P_MM,sub_k-1),
                   ESL_MAX(IMX(i-3,k-1)    + TSC(p7P_IM,sub_k-1),
                   ESL_MAX(DMX(i-3,k-1)    + TSC(p7P_DM,sub_k-1),
                   ESL_MAX(XMX(i-3,p7G_B)  + entry,
                           PVX4(pv_pi,k-1) + TSC_P)))) + rsc_c0[sub_k];

        IMX(i,k) = ESL_MAX(MMX(i-3,k) + TSC(p7P_MI,sub_k),
                           IMX(i-3,k) + TSC(p7P_II,sub_k));
        if (rsc_c0[sub_k] == -eslINFINITY) IMX(i,k) = -eslINFINITY;

        DMX(i,k) = ESL_MAX(MMX(i,k-1) + TSC(p7P_MD,sub_k-1),
                           DMX(i,k-1) + TSC(p7P_DD,sub_k-1));

        XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E),
                       ESL_MAX(MMX(i,k),
                               DMX(i,k))) + exit;
      }

      sub_k = k_start + M - 1;
      MMX(i,M) = ESL_MAX(MMX(i-3,M-1)    + TSC(p7P_MM,sub_k-1),
                 ESL_MAX(IMX(i-3,M-1)    + TSC(p7P_IM,sub_k-1),
                 ESL_MAX(DMX(i-3,M-1)    + TSC(p7P_DM,sub_k-1),
                 ESL_MAX(XMX(i-3,p7G_B)  + entry,
                         PVX4(pv_pi,M-1) + TSC_P)))) + rsc_c0[sub_k];
    }

  IMX(i,M) = -eslINFINITY;

  DMX(i,M) = ESL_MAX(MMX(i,M-1) + TSC(p7P_MD,sub_k-1),
                     DMX(i,M-1) + TSC(p7P_DD,sub_k-1));

  if (!sc->global_end) {
    XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E),
                   ESL_MAX(MMX(i,M),
                           DMX(i,M)));

    XMX(i,p7G_C) = ESL_MAX(XMX(i-3,p7G_C) + gm_tr->xsc[p7P_C][p7P_LOOP],
                           XMX(i,p7G_E)   + gm_tr->xsc[p7P_E][p7P_MOVE]);
  }
  else XMX(i,p7G_C) = -eslINFINITY;

  XMX(i,p7G_J) = -eslINFINITY;
}

/* spliced_donor_row():
 * Donor updates of row <i>: fold the look-back row i-min_intron-3 into
 * the donor scores (don_iv) at donor sites, and record in don_idx the
 * row each score came from. Of equal scores the later row is kept.
 */
static void
spliced_donor_row(const SPLICED_CHK *sc, int i)
{
  const P7_FS_PROFILE *gm_tr   = sc->gm_tr;
  float              **dp      = sc->cx->dp;
  float               *don_ivx = sc->don_iv->ivx;
  int                 *don_idx = sc->cv->ivx;
  int                  M       = sc->M;
  int                  lb      = i - sc->min_intron - 3;
  int                  don0, don1, don2;
  int                  r, s;
  int                  C2[p7P_MAXNUC+1];
  int                  nuc1, nuc3;
  int                  k, sub_k;
  float                TMP_SC, NEW_SC;

  if (i < sc->min_intron + 3) return;

  don2 = donor_at(sc->sub_dsq, sc->i_start, sc->min_intron, i);
  don1 = donor_at(sc->sub_dsq, sc->i_start, sc->min_intron, i-1);
  don0 = donor_at(sc->sub_dsq, sc->i_start, sc->min_intron, i-2);
  if (don0 < 0 && don1 < 0 && don2 < 0) return;

  r = nuc_at(sc->sub_dsq, sc->i_start, i - sc->min_intron - 2);
  s = nuc_at(sc->sub_dsq, sc->i_start, i - sc->min_intron - 1);

  if (don2 >= 0) {
    for (nuc3 = 0; nuc3 < p7P_MAXNUC; nuc3++)
      C2[nuc3] = p7P_MINIDX(p7P_CODON3_FS1(r, s, nuc3), p7P_DEGEN1_C);
    C2[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(r, s, p7P_MAXCODONS1), p7P_DEGEN1_C);

    for (k = 2; k < M; k++) {
      sub_k  = sc->k_start + k - 1;
      TMP_SC = ESL_MAX(MMX(lb,k-1), DMX(lb,k-1));
      for (nuc3 = 0; nuc3 <= p7P_MAXNUC; nuc3++) {
        NEW_SC = TMP_SC + gm_tr->rsc[C2[nuc3]][sub_k];
        if (NEW_SC >= SSX2(k,don2,nuc3)) { SSX2(k,don2,nuc3) = NEW_SC; DIX2(k,don2,nuc3) = lb; }
      }
    }
  }

  if (don1 >= 0) {
    nuc1 = ESL_MIN(r, p7P_MAXNUC);
    for (k = 2; k < M; k++) {
      TMP_SC = ESL_MAX(MMX(lb,k-1), DMX(lb,k-1));
      if (TMP_SC >= SSX1(k,don1,nuc1)) { SSX1(k,don1,nuc1) = TMP_SC; DIX1(k,don1,nuc1) = lb; }
    }
  }

  if (don0 >= 0) {
    for (k = 2; k < M; k++) {
      TMP_SC = ESL_MAX(MMX(lb,k-1), DMX(lb,k-1));
      if (TMP_SC >= SSX0(k,don0)) { SSX0(k,don0) = TMP_SC; DIX0(k,don0) = lb; }
    }
  }
}

/* spliced_rows():
 * Compute rows <from>..<to>: for each row the acceptor P values, then
 * (if <full>) the M, I, D and special states, then the donor updates.
 * With <full> FALSE, only the P values and donor updates are done, for
 * rows whose M, I, D values are already loaded.
 */
static void
spliced_rows(const SPLICED_CHK *sc, int from, int to, int full)
{
  int i;

  for (i = from; i <= to; i++) {
    spliced_acc_row(sc, i);
    if (full) spliced_mdi_row(sc, i);
    spliced_donor_row(sc, i);
  }
}

/* spliced_global_end(): the global exit at row L, as at the end of
 * p7_GViterbi_Spliced(). Row L must be loaded. */
static void
spliced_global_end(const SPLICED_CHK *sc)
{
  float **dp  = sc->cx->dp;
  float  *xmx = sc->cx->xmx;
  int     L   = sc->L;
  int     M   = sc->M;

  XMX(L,p7G_E) = ESL_MAX(MMX(L,M), DMX(L,M));
  XMX(L,p7G_C) = XMX(L,p7G_E) + sc->gm_tr->xsc[p7P_E][p7P_MOVE];
}

/* spliced_check_sizes(): make sure <cx>, <cv> fit an M by L comparison
 * with this <min_intron>. */
static int
spliced_check_sizes(const SPLICED_CHK *sc)
{
  if (sc->gm_tr->codon_lengths != 1)                          ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");
  if (sc->cx->saveR < sc->min_intron + 6)                     ESL_EXCEPTION(eslEINVAL, "checkpoints save %d rows; need %d", sc->cx->saveR, sc->min_intron + 6);
  if (sc->cx->allocW < sc->M + 1 || sc->cx->allocP < sc->L + 1 || sc->cx->allocXR < sc->L + 1) ESL_EXCEPTION(eslEINVAL, "checkpointed matrix too small");
  if (sc->cv->allocCF != SIGNAL_MEM_SIZE || sc->cv->allocCI != SIGNAL_MEM_SIZE) ESL_EXCEPTION(eslEINVAL, "checkpoint values need %d channels", SIGNAL_MEM_SIZE);
  if (sc->cv->allocM < sc->M || sc->cv->allocK < sc->L / sc->cx->blockR)         ESL_EXCEPTION(eslEINVAL, "checkpoint values too small");
  return eslOK;
}


/* Function:  p7_GViterbi_SplicedChk()
 * Synopsis:  Translated spliced Viterbi into a checkpointed matrix
 *
 * Purpose:   Same as p7_GViterbi_Spliced(), but the rows go into the
 *            checkpointed matrix <cx>, with checkpoints of the donor
 *            scores and their rows in <cv>, so that
 *            p7_GViterbi_SplicedTraceChk() can compute any block of rows
 *            again. On return the last block is loaded, and the special
 *            states of every row are in <cx->xmx>.
 *
 *            <cx> must be grown to <k_end>-<k_start>+1 by
 *            <i_end>-<i_start>+1 and save at least <min_intron>+6 rows
 *            per checkpoint. <cv> must have SIGNAL_MEM_SIZE float and
 *            int channels, room for the model length, and (L /
 *            <cx->blockR>) checkpoints.
 *
 * Args:      as p7_GViterbi_Spliced(), with <cx> and <cv> in place of <gx>
 *
 * Return:    <eslOK> on success.
 *            <eslEINVAL> if <cx> or <cv> doesn't fit the comparison.
 */
int
p7_GViterbi_SplicedChk(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, P7_CHKGMX *cx, P7_CHKIVX *cv, P7_IVX *acc_iv, P7_IVX *don_iv, const float *signal_scores, int i_start, int i_end, int k_start, int k_end, int min_intron, int global_start, int global_end)
{
  SPLICED_CHK sc = { sub_dsq, gm_tr, cx, cv, acc_iv, don_iv, signal_scores, i_start, k_start,
                     k_end - k_start + 1, i_end - i_start + 1, min_intron, global_start, global_end };
  int         blockR = cx->blockR;
  int         b, nb, end;
  int         status;

  if ((status = spliced_check_sizes(&sc)) != eslOK) return status;

  cx->M = sc.M;
  cx->L = sc.L;
  if ((status = p7_chkgmx_StartBlock(cx, 0)) != eslOK) return status;
  spliced_init(&sc);

  for (b = 0; b <= sc.L; b += blockR)
    {
      if (b > 0 && (status = p7_chkgmx_StartBlock(cx, b)) != eslOK) return status;
      end = ESL_MIN(sc.L, b + blockR - 1);
      nb  = b + blockR;
      if (nb <= sc.L)
        {   /* save the donor state for block nb at the start of row nb-3 */
          spliced_rows(&sc, ESL_MAX(3, b), nb - 4, TRUE);
          if ((status = p7_chkivx_Save(cv, don_iv, sc.M, nb / blockR - 1)) != eslOK) return status;
          spliced_rows(&sc, nb - 3, end, TRUE);
        }
      else spliced_rows(&sc, ESL_MAX(3, b), end, TRUE);
    }

  if (global_end) spliced_global_end(&sc);
  return eslOK;
}


/* Function:  p7_GViterbi_SplicedTrace()
 * Synopsis:  Create a trace for spliced Viterbi 
 *
 * Purpose:   Create a trace that includes the exon(s) and any splice 
 *            site(s) from the filled translated spliced viterbi matrix 
 *            <gx> and calculate the Viterbi score of the spliced sequence.
 *            Returns the filled trace <tr> and optional Viterbi score 
 *            <vitsc>.
 *
 * Args:      sub_dsq       - nucleotide sequence
 *            gm_tr         - a codon profile.
 *            gx            - filled spliced viterbi DP matrix
 *            signal_scores - scores for splice site signals 
 *            tr            - trace to fill
 *            i_start       - start position on the <sub_dsq>
 *            i_end         - end position on the <sub_dsq>
 *            k_start       - start position on the <gm_tr>
 *            k_end         - end position on the <gm_tr>
 *            min_intron    - minimum intron length
 *            vitsc         - optional return Viterbi score
 *            max_drop      - optional returned x-drop signal
 *
 * Return:    <eslOK> on success.
 *            <eslFAIL> if even the optimal path has zero probability;
 *            in this case, the trace is set blank (<tr->N = 0>).
 *
 */
int
p7_GViterbi_SplicedTrace(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, const P7_GMX *gx, const float *signal_scores, P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end, int min_intron, float *vitsc, float *max_drop)
{

  float const *tsc   = gm_tr->tsc;
  float      **dp    = gx->dp;              /* so {MDI}MX() macros work                              */
  float       *xmx   = gx->xmx;             /* so XMX() macro works                                  */
  float        r_tol = 1e-5;                /* relative tolerance for floating point "equality" test */
  float        a_tol = 1e-4;                /* absolute tolerance for floating point "equality" test */
  int          M     = k_end - k_start + 1; /* sub model length                                      */
  int          L     = i_end - i_start + 1; /* sub seq length                                        */
  int          i     = L;                   /* position in seq (1..L)                                */
  int          k     = 0;                   /* position in model (1..M)                              */
  int          j, z;
  int          t,u,v,w,x;
  int          d1, d2;
  int          c, c0, c1, c2, c3;
  int          sprv,scur;
  int          sub_i, sub_k;
  int          donor_i;
  int          acc[3] = {0, 0, 0};
  int          don_sig;
  int          P_cnt, P_idx;
  float        P_state;
  float        emit, emit0, emit1, emit2;
  float        vsc;
  float        curr_P_score;
  float        curr_M_score;
  float        max_M_score;
  float        curr_max_drop;
  float       *P_scores;
  int          status;

  P_scores = NULL;
  P_cnt    = 0;
  ESL_ALLOC(P_scores, sizeof(float) * M);

  if(gm_tr->codon_lengths != 1) ESL_EXCEPTION(eslEINVAL, "proflie not allocated for 1 codon length");

#if eslDEBUGLEVEL > 0
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif

  vsc = XMX(i, p7G_C) + gm_tr->xsc[p7P_C][p7P_MOVE]; 

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i+i_start-1, 0)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i+i_start-1, 0)) != eslOK) return status;
  
  sprv = p7T_C;
  
  while (sprv != p7T_S) {
    switch (sprv) {
    case p7T_C:     /* C(i) comes from C(i-1) or E(i) */
       
      if      (XMX(i, p7G_C) < XMX(i-2, p7G_C) || XMX(i, p7G_C) < XMX(i-1, p7G_C))                     scur = p7T_C; 
      else if (XMX(i,p7G_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
      else if (esl_FCompare(XMX(i, p7G_C), XMX(i-3, p7G_C) + gm_tr->xsc[p7P_C][p7P_LOOP], r_tol, a_tol) == eslOK) scur = p7T_C;
      else if (esl_FCompare(XMX(i, p7G_C), XMX(i,   p7G_E) + gm_tr->xsc[p7P_E][p7P_MOVE], r_tol, a_tol) == eslOK) scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);

      break;

    case p7T_E:     /* E connects from any M state. k set here */
      if (XMX(i, p7G_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
      for (k = M; k >= 1; k--) {
        if (esl_FCompare(XMX(i, p7G_E), MMX(i,k), r_tol, a_tol) == eslOK) { scur = p7T_M; break; }
        if (esl_FCompare(XMX(i, p7G_E), DMX(i,k), r_tol, a_tol) == eslOK) { scur = p7T_D; break; }
      }
      if (k == 0) ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
 
      break;

    case p7T_M:         /* M connects from i-1,k-1, or B */
      if (MMX(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);
     
      sub_i = i_start + i - 1; 
      
      v = nuc_at(sub_dsq, i_start, i-2);

      w = nuc_at(sub_dsq, i_start, i-1);

      x = nuc_at(sub_dsq, i_start, i);

      c3 = p7P_CODON3_FS1(v, w, x);
      c3 = p7P_MINIDX(c3, p7P_DEGEN1_C);

      sub_k  = k_start + k - 1;
      emit = p7P_MSC_CODON(gm_tr, sub_k, c3);
      
      if      (esl_FCompare(MMX(i,k), MMX(i-3, k-1) + TSC(p7P_MM, sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(MMX(i,k), IMX(i-3, k-1) + TSC(p7P_IM, sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare(MMX(i,k), DMX(i-3, k-1) + TSC(p7P_DM, sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_D;      
      else if (esl_FCompare(MMX(i,k), XMX(i-3, p7G_B)                      + emit, r_tol, a_tol) == eslOK) scur = p7T_B;
      else {    /* The P state is not stored - must be recreated */
        if (i < min_intron+7) ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);

        vsc -= TSC_P;

        /* detect which acceptor sites are available */
        acc[0] = 0;
        if(SIGNAL(sub_dsq[sub_i-7], sub_dsq[sub_i-6]) == ACCEPT_AG) acc[0] = 1;
        if(SIGNAL(sub_dsq[sub_i-7], sub_dsq[sub_i-6]) == ACCEPT_AC) acc[0] = 2; 
       
        acc[1] = 0;
        if(SIGNAL(sub_dsq[sub_i-6], sub_dsq[sub_i-5]) == ACCEPT_AG) acc[1] = 1;
        if(SIGNAL(sub_dsq[sub_i-6], sub_dsq[sub_i-5]) == ACCEPT_AC) acc[1] = 2;

        acc[2] = 0;
        if(SIGNAL(sub_dsq[sub_i-5], sub_dsq[sub_i-4]) == ACCEPT_AG) acc[2] = 1;
        if(SIGNAL(sub_dsq[sub_i-5], sub_dsq[sub_i-4]) == ACCEPT_AC) acc[2] = 2;

        if(acc[0] == 0 && acc[1] == 0 && acc[2] == 0) ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
        
        /* Scan for donor sites */
        for(j = 0; j < i-min_intron-4; j++) {
        
          /* anything other than A,C,G,T gets a placeholder, so it can't form a donor signal */
          d1 = (sub_dsq[sub_i-min_intron-j-1] < p7P_MAXNUC) ? (int)sub_dsq[sub_i-min_intron-j-1] : p7P_MAXCODONS1;
          d2 = (sub_dsq[sub_i-min_intron-j]   < p7P_MAXNUC) ? (int)sub_dsq[sub_i-min_intron-j]   : p7P_MAXCODONS1;
          don_sig = -1;
          if      (SIGNAL(d1, d2) == DONOR_GT) don_sig = p7S_GTAG;
          else if (SIGNAL(d1, d2) == DONOR_GC) don_sig = p7S_GCAG;
          else if (SIGNAL(d1, d2) == DONOR_AT) don_sig = p7S_ATAC;
        
          if(don_sig != -1) {
            /* If we have a donor site get emissions */
            
            t = nuc_at(sub_dsq, i_start, i-min_intron-j-3);

            u = nuc_at(sub_dsq, i_start, i-min_intron-j-2);
            
            v = nuc_at(sub_dsq, i_start, i-5);

            w = nuc_at(sub_dsq, i_start, i-4);
            
            x = nuc_at(sub_dsq, i_start, i-3);

            c2 = p7P_CODON3_FS1(t, u, x);
            c2 = p7P_MINIDX(c2, p7P_DEGEN1_C);
            emit2 = p7P_MSC_CODON(gm_tr, sub_k-1, c2);

            c1 = p7P_CODON3_FS1(u, w, x);
            c1 = p7P_MINIDX(c1, p7P_DEGEN1_C);
            emit1 = p7P_MSC_CODON(gm_tr, sub_k-1, c1);

            c0 = p7P_CODON3_FS1(v, w, x);
            c0 = p7P_MINIDX(c0, p7P_DEGEN1_C);
            emit0 = p7P_MSC_CODON(gm_tr, sub_k-1, c0);

            /* Pair donor and acceptor sites and compare scores */
            if(don_sig == p7S_GTAG || don_sig == p7S_GCAG) {
              if (acc[2] == 1) {
                P_state = ESL_MAX(MMX(i-min_intron-j-4,k-2), DMX(i-min_intron-j-4,k-2)) + signal_scores[don_sig] + emit2;
				if(esl_FCompare(MMX(i,k), P_state + TSC_P + emit, r_tol, a_tol) == eslOK) { scur = p7T_P; c = 2; donor_i = i-min_intron-j-4; vsc -= signal_scores[don_sig]; break; }
              }
              if (acc[1] == 1) {
                P_state = ESL_MAX(MMX(i-min_intron-j-3,k-2), DMX(i-min_intron-j-3,k-2)) + signal_scores[don_sig] + emit1;
                if(esl_FCompare(MMX(i,k), P_state + TSC_P + emit, r_tol, a_tol) == eslOK) { scur = p7T_P; c = 1; donor_i = i-min_intron-j-3;  vsc -= signal_scores[don_sig]; break; }
              }
              if (acc[0] == 1) {
                P_state = ESL_MAX(MMX(i-min_intron-j-2,k-2), DMX(i-min_intron-j-2,k-2)) + signal_scores[don_sig] + emit0;
                if(esl_FCompare(MMX(i,k), P_state + TSC_P + emit, r_tol, a_tol) == eslOK) { scur = p7T_P; c = 0; donor_i = i-min_intron-j-2;  vsc -= signal_scores[don_sig]; break; }
              }
            }
            else if(don_sig == p7S_ATAC) {
              if (acc[2] == 2) {
                P_state = ESL_MAX(MMX(i-min_intron-j-4,k-2), DMX(i-min_intron-j-4,k-2)) + signal_scores[don_sig] + emit2;
                if(esl_FCompare(MMX(i,k), P_state + TSC_P + emit, r_tol, a_tol) == eslOK) { scur = p7T_P; c = 2; donor_i = i-min_intron-j-4;  vsc -= signal_scores[don_sig]; break; }
              }
              if (acc[1] == 2) {
                P_state = ESL_MAX(MMX(i-min_intron-j-3,k-2), DMX(i-min_intron-j-3,k-2)) + signal_scores[don_sig] + emit1;
                if(esl_FCompare(MMX(i,k), P_state + TSC_P + emit, r_tol, a_tol) == eslOK) { scur = p7T_P; c = 1; donor_i = i-min_intron-j-3;  vsc -= signal_scores[don_sig]; break; }
              }
              if (acc[0] == 2) {
                P_state = ESL_MAX(MMX(i-min_intron-j-2,k-2), DMX(i-min_intron-j-2,k-2)) + signal_scores[don_sig] + emit0;
                if(esl_FCompare(MMX(i,k), P_state + TSC_P + emit, r_tol, a_tol) == eslOK) { scur = p7T_P; c = 0; donor_i = i-min_intron-j-2;  vsc -= signal_scores[don_sig]; break; }
              }
            }
          }
        }
        if(scur != p7T_P) ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);

        /* Track P state scores for x_drop calculation */
        P_scores[P_cnt++] = TSC_P + signal_scores[don_sig];
      }

      k--; i-=3;
      break;

    case p7T_D:         /* D connects from M,D at i,k-1 */
      if (DMX(i, k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);
     
      sub_k  = k_start + k - 1; 
      if      (esl_FCompare(DMX(i,k), MMX(i, k-1) + TSC(p7P_MD, sub_k-1), r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(DMX(i,k), DMX(i, k-1) + TSC(p7P_DD, sub_k-1), r_tol, a_tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      k--;
      break;

    case p7T_I:         /* I connects from M,I at i-1,k*/
      if (IMX(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);

      sub_k  = k_start + k - 1;
      if      (esl_FCompare(IMX(i,k), MMX(i-3,k) + TSC(p7P_MI, sub_k), r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(IMX(i,k), IMX(i-3,k) + TSC(p7P_II, sub_k), r_tol, a_tol) == eslOK) scur = p7T_I;
      else ESL_EXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      i-=3;
      break;

    case p7T_P:
      if (MMX(donor_i,k-1) > DMX(donor_i,k-1)) scur = p7T_M;
      else                                     scur = p7T_D;

      k--; i = donor_i; 
      break;

    case p7T_N:         /* N connects from S, N */
      if (XMX(i, p7G_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = ( (i == 0) ? p7T_S : p7T_N);
      break;

    case p7T_B:         /* B connects from N, J */
      vsc += p7P_TSC(gm_tr, k, p7P_BM); 
      
      if (XMX(i,p7G_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);
      if      (esl_FCompare(XMX(i,p7G_B), XMX(i, p7G_N) + gm_tr->xsc[p7P_N][p7P_MOVE], r_tol, a_tol) == eslOK) scur = p7T_N;
      else  ESL_EXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

      default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over statetype[tpos-1] */

    if      (scur == p7T_M) c = 3;
    else if (scur != p7T_P) c = 0;

    if ((status = p7_trace_fs_Append(tr, scur, k_start+k-1, i_start+i-1, c)) != eslOK) return status;
    
    /* For NC, we had to defer i decrement. */
    if ( (scur == p7T_N || scur == p7T_C) && scur == sprv) i--; 
   
    sprv = scur;
  } /* end traceback, at S state */

  tr->M = M;
  tr->L = L;

  if(vitsc != NULL) *vitsc = vsc;

  status = p7_trace_fs_Reverse(tr);

  /* When max_drop != NULL calculate the maximum M state score drop found
   * in the alignment for comparison to an x-drop threshold. */
  if (status == eslOK && max_drop != NULL) {
    P_idx         = P_cnt - 1;
    curr_P_score  = 0.0f;
    max_M_score   = -eslINFINITY;
    curr_max_drop = -eslINFINITY;

    for (z = 0; z < tr->N; z++) {
      /* Track the current P state score contribution so that it can be
       * subtracted out from M state scores to ensure the curr_max_drop
       * only looks for drops in M state scores, not splice penalties */
      if (tr->st[z] == p7T_P) {
        curr_P_score += P_scores[P_idx--];
      }
      else if (tr->st[z] == p7T_M) {
        i = tr->i[z] - i_start + 1;
        k = tr->k[z] - k_start + 1;
        curr_M_score = MMX(i, k) - curr_P_score;
        if (max_M_score != -eslINFINITY)
          curr_max_drop = ESL_MAX(curr_max_drop, max_M_score - curr_M_score);
        max_M_score = ESL_MAX(max_M_score, curr_M_score);
      }
    }
    *max_drop = curr_max_drop;
  }

  free(P_scores);

  return status;

 ERROR:
  free(P_scores);
  return status;
}


/* The donor behind a P state score, as chosen by replay_donors() */
typedef struct {
  float sc;     /* P state score                                       */
  int   row;    /* donor row: the last row before the split codon      */
  int   c;      /* nucleotides of the split codon before the donor     */
  int   sig;    /* splice signal                                       */
} DONOR_CHOICE;

/* keep_donor(): replace <best> with the given donor if it scores higher,
 * or scores the same and is nearer the acceptor. */
static inline void
keep_donor(DONOR_CHOICE *best, float sc, int row, int c, int sig)
{
  if (sc > best->sc || (sc == best->sc && row > best->row)) {
    best->sc  = sc;
    best->row = row;
    best->c   = c;
    best->sig = sig;
  }
}

/* Donor slot scores for one model column, laid out as SSX0/1/2 in the
 * fill, with the row each score came from. */
typedef struct {
  float sc0[p7S_SPLICE_SIGNALS];
  float sc1[p7P_MAXNUC+1][p7S_SPLICE_SIGNALS];
  float sc2[p7P_MAXNUC+1][p7S_SPLICE_SIGNALS];
  int   row0[p7S_SPLICE_SIGNALS];
  int   row1[p7P_MAXNUC+1][p7S_SPLICE_SIGNALS];
  int   row2[p7P_MAXNUC+1][p7S_SPLICE_SIGNALS];
} DONOR_SLOTS;

/* replay_donors_from():
 *
 * Find the donor behind the P state score that the spliced Viterbi fill
 * built at acceptor row <acc_i> and model column <kk>. M(acc_i+3, kk+1)
 * enters through that P state.
 *
 * The fill keeps only the best score of each donor slot, not the row it
 * came from, and later donors overwrite it. Starting from the slot
 * scores and rows in <ds> as they were at the start of row <s0>, this
 * repeats the fill's donor updates for column <kk> alone up to
 * <acc_i>, reading each look-back row through <dp>, and records each
 * slot's row along with its score. Of equal scores the later (nearer)
 * donor is kept. The acceptor is then scored as in the fill, again
 * keeping the nearer donor of equal scores. <ds> is modified.
 *
 * Returns <eslOK> with the P state score in <ret_sc>, the donor row
 * (the last row before the split codon) in <ret_donor_i>, the number of
 * nucleotides of the split codon before the donor in <ret_c>, and the
 * splice signal in <ret_don_sig>. Returns <eslFAIL> if no donor reaches
 * this acceptor.
 */
static int
replay_donors_from(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, float **dp, const float *signal_scores,
                   int i_start, int k_start, int min_intron, int s0, DONOR_SLOTS *ds, int acc_i, int kk,
                   float *ret_sc, int *ret_donor_i, int *ret_c, int *ret_don_sig)
{
  int          sub_kk = k_start + kk - 1;
  int          don0, don1, don2;
  int          acc0, acc1, acc2;
  int          r, s, v, w, x;
  int          nuc, sig, lb, iu, C;
  float        tsc, sc;
  DONOR_CHOICE best = { -eslINFINITY, -1, -1, -1 };

  /* The fill updates donors only for columns 2..M-1 and only from row
   * min_intron+3 on, and scores acceptors from that row on. */
  if (kk < 2 || acc_i < min_intron + 3) return eslFAIL;

  /* Donor updates, as in the fill. At row iu the donor dinucleotide
   * that starts at row iu-min_intron is don2, and the look-back row lb
   * is the same for don2, don1 and don0. */
  for (iu = s0; iu < acc_i; iu++) {
    don2 = donor_at(sub_dsq, i_start, min_intron, iu);
    don1 = donor_at(sub_dsq, i_start, min_intron, iu - 1);
    don0 = donor_at(sub_dsq, i_start, min_intron, iu - 2);
    if (don0 < 0 && don1 < 0 && don2 < 0) continue;

    r   = nuc_at(sub_dsq, i_start, iu - min_intron - 2);
    lb  = iu - min_intron - 3;
    tsc = ESL_MAX(MMX(lb, kk-1), DMX(lb, kk-1));

    if (don2 >= 0) {
      s = nuc_at(sub_dsq, i_start, iu - min_intron - 1);
      for (nuc = 0; nuc <= p7P_MAXNUC; nuc++) {
        C  = p7P_MINIDX(p7P_CODON3_FS1(r, s, (nuc < p7P_MAXNUC ? nuc : p7P_MAXCODONS1)), p7P_DEGEN1_C);
        sc = tsc + gm_tr->rsc[C][sub_kk];
        if (sc >= ds->sc2[nuc][don2]) { ds->sc2[nuc][don2] = sc; ds->row2[nuc][don2] = lb; }
      }
    }
    if (don1 >= 0) {
      nuc = ESL_MIN(r, p7P_MAXNUC);
      if (tsc >= ds->sc1[nuc][don1]) { ds->sc1[nuc][don1] = tsc; ds->row1[nuc][don1] = lb; }
    }
    if (don0 >= 0) {
      if (tsc >= ds->sc0[don0]) { ds->sc0[don0] = tsc; ds->row0[don0] = lb; }
    }
  }

  /* Acceptor at row acc_i, as in the fill's P state */
  v = nuc_at(sub_dsq, i_start, acc_i - 2);
  w = nuc_at(sub_dsq, i_start, acc_i - 1);
  x = nuc_at(sub_dsq, i_start, acc_i);
  acc2 = acceptor_at(sub_dsq, i_start, acc_i);
  acc1 = acceptor_at(sub_dsq, i_start, acc_i - 1);
  acc0 = acceptor_at(sub_dsq, i_start, acc_i - 2);

  if (acc0 >= 0) {
    C = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);
    for (sig = 0; sig < p7S_SPLICE_SIGNALS; sig++) {
      if ((sig == p7S_ATAC) != (acc0 == ACCEPT_AC)) continue;
      sc = ds->sc0[sig] + signal_scores[sig] + gm_tr->rsc[C][sub_kk];
      keep_donor(&best, sc, ds->row0[sig], 0, sig);
    }
  }
  if (acc1 >= 0) {
    for (nuc = 0; nuc <= p7P_MAXNUC; nuc++) {
      C = p7P_MINIDX(p7P_CODON3_FS1((nuc < p7P_MAXNUC ? nuc : p7P_MAXCODONS1), w, x), p7P_DEGEN1_C);
      for (sig = 0; sig < p7S_SPLICE_SIGNALS; sig++) {
        if ((sig == p7S_ATAC) != (acc1 == ACCEPT_AC)) continue;
        sc = ds->sc1[nuc][sig] + signal_scores[sig] + gm_tr->rsc[C][sub_kk];
        keep_donor(&best, sc, ds->row1[nuc][sig], 1, sig);
      }
    }
  }
  if (acc2 >= 0) {
    nuc = ESL_MIN(x, p7P_MAXNUC);
    for (sig = 0; sig < p7S_SPLICE_SIGNALS; sig++) {
      if ((sig == p7S_ATAC) != (acc2 == ACCEPT_AC)) continue;
      sc = ds->sc2[nuc][sig] + signal_scores[sig];
      keep_donor(&best, sc, ds->row2[nuc][sig], 2, sig);
    }
  }
  if (best.sc == -eslINFINITY) return eslFAIL;

  *ret_sc      = best.sc;
  *ret_donor_i = best.row;
  *ret_c       = best.c;
  *ret_don_sig = best.sig;
  return eslOK;
}

/* donor_slots_clear(): no donors seen yet. */
static void
donor_slots_clear(DONOR_SLOTS *ds)
{
  int nuc, sig;

  for (sig = 0; sig < p7S_SPLICE_SIGNALS; sig++) {
    ds->sc0[sig] = -eslINFINITY; ds->row0[sig] = -1;
    for (nuc = 0; nuc <= p7P_MAXNUC; nuc++) {
      ds->sc1[nuc][sig] = ds->sc2[nuc][sig] = -eslINFINITY;
      ds->row1[nuc][sig] = ds->row2[nuc][sig] = -1;
    }
  }
}

/* replay_donors():
 *
 * replay_donors_from() for the full matrix <gx> of p7_GViterbi_Spliced(),
 * replaying from the first row that has donor updates.
 */
static int
replay_donors(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, const P7_GMX *gx, const float *signal_scores,
              int i_start, int k_start, int min_intron, int acc_i, int kk,
              float *ret_sc, int *ret_donor_i, int *ret_c, int *ret_don_sig)
{
  DONOR_SLOTS ds;

  donor_slots_clear(&ds);
  return replay_donors_from(sub_dsq, gm_tr, gx->dp, signal_scores, i_start, k_start, min_intron, min_intron + 3, &ds,
                            acc_i, kk, ret_sc, ret_donor_i, ret_c, ret_don_sig);
}


/* Function:  p7_GViterbi_SplicedTraceReplay()
 * Synopsis:  Spliced Viterbi traceback that replays donor updates
 *
 * Purpose:   Same as p7_GViterbi_SplicedTrace(), except for how a P
 *            state's donor is found. Instead of searching back through
 *            earlier rows for a donor whose score matches, it replays
 *            the fill's donor updates for the one model column up to
 *            the acceptor (see replay_donors()), which gives the donor
 *            that produced the P state score. Of equal scores the
 *            nearer donor is used.
 *
 * Args:      as p7_GViterbi_SplicedTrace()
 *
 * Return:    <eslOK> on success.
 *            <eslFAIL> if even the optimal path has zero probability;
 *            in this case, the trace is set blank (<tr->N = 0>).
 */
int
p7_GViterbi_SplicedTraceReplay(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, const P7_GMX *gx, const float *signal_scores, P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end, int min_intron, float *vitsc, float *max_drop)
{

  float const *tsc   = gm_tr->tsc;
  float      **dp    = gx->dp;              /* so {MDI}MX() macros work                              */
  float       *xmx   = gx->xmx;             /* so XMX() macro works                                  */
  float        r_tol = 1e-5;                /* relative tolerance for floating point "equality" test */
  float        a_tol = 1e-4;                /* absolute tolerance for floating point "equality" test */
  int          M     = k_end - k_start + 1; /* sub model length                                      */
  int          L     = i_end - i_start + 1; /* sub seq length                                        */
  int          i     = L;                   /* position in seq (1..L)                                */
  int          k     = 0;                   /* position in model (1..M)                              */
  int          z;
  int          v,w,x;
  int          c, c3;
  int          sprv,scur;
  int          sub_k;
  int          donor_i;
  int          don_sig;
  int          P_cnt, P_idx;
  float        P_state;
  float        emit;
  float        vsc;
  float        curr_P_score;
  float        curr_M_score;
  float        max_M_score;
  float        curr_max_drop;
  float       *P_scores;
  int          status;

  P_scores = NULL;
  P_cnt    = 0;
  ESL_ALLOC(P_scores, sizeof(float) * M);

  if(gm_tr->codon_lengths != 1) ESL_EXCEPTION(eslEINVAL, "proflie not allocated for 1 codon length");

#if eslDEBUGLEVEL > 0
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif

  vsc = XMX(i, p7G_C) + gm_tr->xsc[p7P_C][p7P_MOVE];

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i+i_start-1, 0)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i+i_start-1, 0)) != eslOK) return status;

  sprv = p7T_C;

  while (sprv != p7T_S) {
    switch (sprv) {
    case p7T_C:     /* C(i) comes from C(i-1) or E(i) */

      if      (XMX(i, p7G_C) < XMX(i-2, p7G_C) || XMX(i, p7G_C) < XMX(i-1, p7G_C))                     scur = p7T_C;
      else if (XMX(i,p7G_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
      else if (esl_FCompare(XMX(i, p7G_C), XMX(i-3, p7G_C) + gm_tr->xsc[p7P_C][p7P_LOOP], r_tol, a_tol) == eslOK) scur = p7T_C;
      else if (esl_FCompare(XMX(i, p7G_C), XMX(i,   p7G_E) + gm_tr->xsc[p7P_E][p7P_MOVE], r_tol, a_tol) == eslOK) scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);

      break;

    case p7T_E:     /* E connects from any M state. k set here */
      if (XMX(i, p7G_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
      for (k = M; k >= 1; k--) {
        if (esl_FCompare(XMX(i, p7G_E), MMX(i,k), r_tol, a_tol) == eslOK) { scur = p7T_M; break; }
        if (esl_FCompare(XMX(i, p7G_E), DMX(i,k), r_tol, a_tol) == eslOK) { scur = p7T_D; break; }
      }
      if (k == 0) ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);

      break;

    case p7T_M:         /* M connects from i-1,k-1, B, or a P state at i-3,k-1 */
      if (MMX(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);

      v = nuc_at(sub_dsq, i_start, i-2);

      w = nuc_at(sub_dsq, i_start, i-1);

      x = nuc_at(sub_dsq, i_start, i);

      c3 = p7P_CODON3_FS1(v, w, x);
      c3 = p7P_MINIDX(c3, p7P_DEGEN1_C);

      sub_k  = k_start + k - 1;
      emit = p7P_MSC_CODON(gm_tr, sub_k, c3);

      if      (esl_FCompare(MMX(i,k), MMX(i-3, k-1) + TSC(p7P_MM, sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(MMX(i,k), IMX(i-3, k-1) + TSC(p7P_IM, sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare(MMX(i,k), DMX(i-3, k-1) + TSC(p7P_DM, sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_D;
      else if (esl_FCompare(MMX(i,k), XMX(i-3, p7G_B)                      + emit, r_tol, a_tol) == eslOK) scur = p7T_B;
      else {    /* The P state is not stored - recreate it from the donor updates */
        if (replay_donors(sub_dsq, gm_tr, gx, signal_scores, i_start, k_start, min_intron, i-3, k-1, &P_state, &donor_i, &c, &don_sig) != eslOK)
          ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
        if (esl_FCompare(MMX(i,k), P_state + TSC_P + emit, r_tol, a_tol) != eslOK)
          ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);

        scur = p7T_P;
        vsc -= TSC_P;
        vsc -= signal_scores[don_sig];

        /* Track P state scores for x_drop calculation */
        P_scores[P_cnt++] = TSC_P + signal_scores[don_sig];
      }

      k--; i-=3;
      break;

    case p7T_D:         /* D connects from M,D at i,k-1 */
      if (DMX(i, k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);

      sub_k  = k_start + k - 1;
      if      (esl_FCompare(DMX(i,k), MMX(i, k-1) + TSC(p7P_MD, sub_k-1), r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(DMX(i,k), DMX(i, k-1) + TSC(p7P_DD, sub_k-1), r_tol, a_tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      k--;
      break;

    case p7T_I:         /* I connects from M,I at i-1,k*/
      if (IMX(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);

      sub_k  = k_start + k - 1;
      if      (esl_FCompare(IMX(i,k), MMX(i-3,k) + TSC(p7P_MI, sub_k), r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(IMX(i,k), IMX(i-3,k) + TSC(p7P_II, sub_k), r_tol, a_tol) == eslOK) scur = p7T_I;
      else ESL_EXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      i-=3;
      break;

    case p7T_P:
      if (MMX(donor_i,k-1) > DMX(donor_i,k-1)) scur = p7T_M;
      else                                     scur = p7T_D;

      k--; i = donor_i;
      break;

    case p7T_N:         /* N connects from S, N */
      if (XMX(i, p7G_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = ( (i == 0) ? p7T_S : p7T_N);
      break;

    case p7T_B:         /* B connects from N, J */
      vsc += p7P_TSC(gm_tr, k, p7P_BM);

      if (XMX(i,p7G_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);
      if      (esl_FCompare(XMX(i,p7G_B), XMX(i, p7G_N) + gm_tr->xsc[p7P_N][p7P_MOVE], r_tol, a_tol) == eslOK) scur = p7T_N;
      else  ESL_EXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

      default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over statetype[tpos-1] */

    if      (scur == p7T_M) c = 3;
    else if (scur != p7T_P) c = 0;

    if ((status = p7_trace_fs_Append(tr, scur, k_start+k-1, i_start+i-1, c)) != eslOK) return status;

    /* For NC, we had to defer i decrement. */
    if ( (scur == p7T_N || scur == p7T_C) && scur == sprv) i--;

    sprv = scur;
  } /* end traceback, at S state */

  tr->M = M;
  tr->L = L;

  if(vitsc != NULL) *vitsc = vsc;

  status = p7_trace_fs_Reverse(tr);

  /* When max_drop != NULL calculate the maximum M state score drop found
   * in the alignment for comparison to an x-drop threshold. */
  if (status == eslOK && max_drop != NULL) {
    P_idx         = P_cnt - 1;
    curr_P_score  = 0.0f;
    max_M_score   = -eslINFINITY;
    curr_max_drop = -eslINFINITY;

    for (z = 0; z < tr->N; z++) {
      /* Track the current P state score contribution so that it can be
       * subtracted out from M state scores to ensure the curr_max_drop
       * only looks for drops in M state scores, not splice penalties */
      if (tr->st[z] == p7T_P) {
        curr_P_score += P_scores[P_idx--];
      }
      else if (tr->st[z] == p7T_M) {
        i = tr->i[z] - i_start + 1;
        k = tr->k[z] - k_start + 1;
        curr_M_score = MMX(i, k) - curr_P_score;
        if (max_M_score != -eslINFINITY)
          curr_max_drop = ESL_MAX(curr_max_drop, max_M_score - curr_M_score);
        max_M_score = ESL_MAX(max_M_score, curr_M_score);
      }
    }
    *max_drop = curr_max_drop;
  }

  free(P_scores);

  return status;

 ERROR:
  free(P_scores);
  return status;
}

/* spliced_load_block():
 * Make block <b> the loaded block of <sc->cx> and compute its rows
 * again, as the fill did.
 */
static int
spliced_load_block(const SPLICED_CHK *sc, int b)
{
  int end = ESL_MIN(sc->L, b + sc->cx->blockR - 1);
  int status;

  if ((status = p7_chkgmx_LoadBlock(sc->cx, b)) != eslOK) return status;

  if (b == 0)
    {
      spliced_init(sc);
      spliced_rows(sc, 3, end, TRUE);
    }
  else
    {
      if ((status = p7_chkivx_Load(sc->cv, sc->don_iv, sc->M, b / sc->cx->blockR - 1)) != eslOK) return status;
      spliced_rows(sc, b - 3, b - 1, FALSE);
      spliced_rows(sc, b, end, TRUE);
    }

  if (sc->global_end && end == sc->L) spliced_global_end(sc);
  return eslOK;
}

/* spliced_need_row():
 * Make sure row <i> and the rows it reads are loaded: if <i> is before
 * the loaded block, load the block that contains it.
 */
static int
spliced_need_row(const SPLICED_CHK *sc, int i)
{
  int b = (i / sc->cx->blockR) * sc->cx->blockR;

  if (b == sc->cx->b) return eslOK;
  return spliced_load_block(sc, b);
}

/* replay_donors_chk():
 * replay_donors_from() for the loaded block of <sc->cx>, replaying from
 * the start of that block's checkpoint (row b-3), or from the first row
 * with donor updates for block 0. <acc_i> must be >= b-3.
 */
static int
replay_donors_chk(const SPLICED_CHK *sc, int acc_i, int kk, float *ret_sc, int *ret_donor_i, int *ret_c, int *ret_don_sig)
{
  P7_CHKIVX  *cv = sc->cv;
  int         b  = sc->cx->b;
  DONOR_SLOTS ds;
  int         s0;
  int         ck;
  int         nuc, sig;
  size_t      base;

  if (b == 0)
    {
      donor_slots_clear(&ds);
      s0 = sc->min_intron + 3;
    }
  else
    {   /* checkpoint values of column kk, in P7_IVX layout [kk*SIGNAL_MEM_SIZE + slot] */
      ck   = b / sc->cx->blockR - 1;
      base = (size_t) kk * SIGNAL_MEM_SIZE;
      for (sig = 0; sig < p7S_SPLICE_SIGNALS; sig++) {
        ds.sc0[sig]  = cv->ckf[(size_t) ck * (cv->allocM+1) * cv->allocCF + base + sig];
        ds.row0[sig] = cv->cki[(size_t) ck * (cv->allocM+1) * cv->allocCI + base + sig];
        for (nuc = 0; nuc <= p7P_MAXNUC; nuc++) {
          ds.sc1[nuc][sig]  = cv->ckf[(size_t) ck * (cv->allocM+1) * cv->allocCF + base + SPLICE_OFFSET_1 + nuc * p7S_SPLICE_SIGNALS + sig];
          ds.row1[nuc][sig] = cv->cki[(size_t) ck * (cv->allocM+1) * cv->allocCI + base + SPLICE_OFFSET_1 + nuc * p7S_SPLICE_SIGNALS + sig];
          ds.sc2[nuc][sig]  = cv->ckf[(size_t) ck * (cv->allocM+1) * cv->allocCF + base + SPLICE_OFFSET_2 + nuc * p7S_SPLICE_SIGNALS + sig];
          ds.row2[nuc][sig] = cv->cki[(size_t) ck * (cv->allocM+1) * cv->allocCI + base + SPLICE_OFFSET_2 + nuc * p7S_SPLICE_SIGNALS + sig];
        }
      }
      s0 = b - 3;
    }

  return replay_donors_from(sc->sub_dsq, sc->gm_tr, sc->cx->dp, sc->signal_scores, sc->i_start, sc->k_start, sc->min_intron,
                            s0, &ds, acc_i, kk, ret_sc, ret_donor_i, ret_c, ret_don_sig);
}


/* Function:  p7_GViterbi_SplicedTraceChk()
 * Synopsis:  Spliced Viterbi traceback from a checkpointed matrix
 *
 * Purpose:   Same as p7_GViterbi_SplicedTraceReplay(), for a matrix
 *            filled by p7_GViterbi_SplicedChk() with the same arguments.
 *            The trace runs from the last block back to the first,
 *            computing each block it passes through again from its
 *            checkpoint, and skipping blocks that lie inside introns.
 *            Each P state's donor comes from replaying donor updates
 *            from the checkpoint of the block that holds the acceptor.
 *            The M state scores on the trace are recorded as it goes,
 *            for <max_drop>.
 *
 * Args:      as p7_GViterbi_SplicedTraceReplay(), with <cx> and <cv> in
 *            place of <gx>, plus <acc_iv>, <don_iv>, <global_start> and
 *            <global_end> as passed to p7_GViterbi_SplicedChk(), for
 *            computing blocks again.
 *
 * Return:    <eslOK> on success.
 *            <eslFAIL> if even the optimal path has zero probability;
 *            in this case, the trace is set blank (<tr->N = 0>).
 */
int
p7_GViterbi_SplicedTraceChk(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, P7_CHKGMX *cx, P7_CHKIVX *cv, P7_IVX *acc_iv, P7_IVX *don_iv, const float *signal_scores, P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end, int min_intron, int global_start, int global_end, float *vitsc, float *max_drop)
{
  SPLICED_CHK  sc    = { sub_dsq, gm_tr, cx, cv, acc_iv, don_iv, signal_scores, i_start, k_start,
                         k_end - k_start + 1, i_end - i_start + 1, min_intron, global_start, global_end };
  float const *tsc   = gm_tr->tsc;
  float      **dp    = cx->dp;              /* so {MDI}MX() macros work                              */
  float       *xmx   = cx->xmx;             /* so XMX() macro works                                  */
  float        r_tol = 1e-5;                /* relative tolerance for floating point "equality" test */
  float        a_tol = 1e-4;                /* absolute tolerance for floating point "equality" test */
  int          M     = k_end - k_start + 1; /* sub model length                                      */
  int          L     = i_end - i_start + 1; /* sub seq length                                        */
  int          i     = L;                   /* position in seq (1..L)                                */
  int          k     = 0;                   /* position in model (1..M)                              */
  int          z;
  int          v,w,x;
  int          c, c3;
  int          sprv,scur;
  int          sub_k;
  int          donor_i;
  int          don_sig;
  int          P_cnt, P_idx;
  float        P_state;
  float        emit;
  float        vsc;
  float        curr_P_score;
  float        curr_M_score;
  float        max_M_score;
  float        curr_max_drop;
  float       *P_scores;
  float       *M_scores;                    /* M state score of each trace position, in trace order  */
  int          M_alloc;
  void        *p;
  int          status;

  P_scores = NULL;
  M_scores = NULL;
  P_cnt    = 0;
  M_alloc  = 256;
  ESL_ALLOC(P_scores, sizeof(float) * M);
  ESL_ALLOC(M_scores, sizeof(float) * M_alloc);

  if(gm_tr->codon_lengths != 1) ESL_XEXCEPTION(eslEINVAL, "proflie not allocated for 1 codon length");
  if (cx->M != M || cx->L != L || cx->b != (L / cx->blockR) * cx->blockR) ESL_XEXCEPTION(eslEINVAL, "matrix isn't filled for this comparison");

#if eslDEBUGLEVEL > 0
  if (tr->N != 0) ESL_XEXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif

  vsc = XMX(i, p7G_C) + gm_tr->xsc[p7P_C][p7P_MOVE];

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i+i_start-1, 0)) != eslOK) goto ERROR;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i+i_start-1, 0)) != eslOK) goto ERROR;
  M_scores[0] = M_scores[1] = 0.;

  sprv = p7T_C;

  while (sprv != p7T_S) {
    switch (sprv) {
    case p7T_C:     /* C(i) comes from C(i-1) or E(i) */

      if      (XMX(i, p7G_C) < XMX(i-2, p7G_C) || XMX(i, p7G_C) < XMX(i-1, p7G_C))                     scur = p7T_C;
      else if (XMX(i,p7G_C) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
      else if (esl_FCompare(XMX(i, p7G_C), XMX(i-3, p7G_C) + gm_tr->xsc[p7P_C][p7P_LOOP], r_tol, a_tol) == eslOK) scur = p7T_C;
      else if (esl_FCompare(XMX(i, p7G_C), XMX(i,   p7G_E) + gm_tr->xsc[p7P_E][p7P_MOVE], r_tol, a_tol) == eslOK) scur = p7T_E;
      else ESL_XEXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);

      break;

    case p7T_E:     /* E connects from any M state. k set here */
      if (XMX(i, p7G_E) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
      if ((status = spliced_need_row(&sc, i)) != eslOK) goto ERROR;
      for (k = M; k >= 1; k--) {
        if (esl_FCompare(XMX(i, p7G_E), MMX(i,k), r_tol, a_tol) == eslOK) { scur = p7T_M; break; }
        if (esl_FCompare(XMX(i, p7G_E), DMX(i,k), r_tol, a_tol) == eslOK) { scur = p7T_D; break; }
      }
      if (k == 0) ESL_XEXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);

      break;

    case p7T_M:         /* M connects from i-1,k-1, B, or a P state at i-3,k-1 */
      if ((status = spliced_need_row(&sc, i)) != eslOK) goto ERROR;
      if (MMX(i,k) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);

      v = nuc_at(sub_dsq, i_start, i-2);

      w = nuc_at(sub_dsq, i_start, i-1);

      x = nuc_at(sub_dsq, i_start, i);

      c3 = p7P_CODON3_FS1(v, w, x);
      c3 = p7P_MINIDX(c3, p7P_DEGEN1_C);

      sub_k  = k_start + k - 1;
      emit = p7P_MSC_CODON(gm_tr, sub_k, c3);

      if      (esl_FCompare(MMX(i,k), MMX(i-3, k-1) + TSC(p7P_MM, sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(MMX(i,k), IMX(i-3, k-1) + TSC(p7P_IM, sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare(MMX(i,k), DMX(i-3, k-1) + TSC(p7P_DM, sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_D;
      else if (esl_FCompare(MMX(i,k), XMX(i-3, p7G_B)                      + emit, r_tol, a_tol) == eslOK) scur = p7T_B;
      else {    /* The P state is not stored - recreate it from the donor updates */
        if (replay_donors_chk(&sc, i-3, k-1, &P_state, &donor_i, &c, &don_sig) != eslOK)
          ESL_XEXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
        if (esl_FCompare(MMX(i,k), P_state + TSC_P + emit, r_tol, a_tol) != eslOK)
          ESL_XEXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);

        scur = p7T_P;
        vsc -= TSC_P;
        vsc -= signal_scores[don_sig];

        /* Track P state scores for x_drop calculation */
        P_scores[P_cnt++] = TSC_P + signal_scores[don_sig];
      }

      k--; i-=3;
      break;

    case p7T_D:         /* D connects from M,D at i,k-1 */
      if ((status = spliced_need_row(&sc, i)) != eslOK) goto ERROR;
      if (DMX(i, k) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);

      sub_k  = k_start + k - 1;
      if      (esl_FCompare(DMX(i,k), MMX(i, k-1) + TSC(p7P_MD, sub_k-1), r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(DMX(i,k), DMX(i, k-1) + TSC(p7P_DD, sub_k-1), r_tol, a_tol) == eslOK) scur = p7T_D;
      else ESL_XEXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      k--;
      break;

    case p7T_I:         /* I connects from M,I at i-1,k*/
      if ((status = spliced_need_row(&sc, i)) != eslOK) goto ERROR;
      if (IMX(i,k) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);

      sub_k  = k_start + k - 1;
      if      (esl_FCompare(IMX(i,k), MMX(i-3,k) + TSC(p7P_MI, sub_k), r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(IMX(i,k), IMX(i-3,k) + TSC(p7P_II, sub_k), r_tol, a_tol) == eslOK) scur = p7T_I;
      else ESL_XEXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      i-=3;
      break;

    case p7T_P:
      if ((status = spliced_need_row(&sc, donor_i)) != eslOK) goto ERROR;
      if (MMX(donor_i,k-1) > DMX(donor_i,k-1)) scur = p7T_M;
      else                                     scur = p7T_D;

      k--; i = donor_i;
      break;

    case p7T_N:         /* N connects from S, N */
      if (XMX(i, p7G_N) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = ( (i == 0) ? p7T_S : p7T_N);
      break;

    case p7T_B:         /* B connects from N, J */
      vsc += p7P_TSC(gm_tr, k, p7P_BM);

      if (XMX(i,p7G_B) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible B reached at i=%d", i);
      if      (esl_FCompare(XMX(i,p7G_B), XMX(i, p7G_N) + gm_tr->xsc[p7P_N][p7P_MOVE], r_tol, a_tol) == eslOK) scur = p7T_N;
      else  ESL_XEXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

      default: ESL_XEXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over statetype[tpos-1] */

    if      (scur == p7T_M) c = 3;
    else if (scur != p7T_P) c = 0;

    if ((status = p7_trace_fs_Append(tr, scur, k_start+k-1, i_start+i-1, c)) != eslOK) goto ERROR;

    /* Record the M state's score now, while its row is loaded */
    if (tr->N > M_alloc) {
      M_alloc *= 2;
      ESL_RALLOC(M_scores, p, sizeof(float) * M_alloc);
    }
    M_scores[tr->N-1] = (scur == p7T_M ? MMX(i,k) : 0.);

    /* For NC, we had to defer i decrement. */
    if ( (scur == p7T_N || scur == p7T_C) && scur == sprv) i--;

    sprv = scur;
  } /* end traceback, at S state */

  tr->M = M;
  tr->L = L;

  if(vitsc != NULL) *vitsc = vsc;

  status = p7_trace_fs_Reverse(tr);

  /* When max_drop != NULL calculate the maximum M state score drop found
   * in the alignment for comparison to an x-drop threshold. Trace
   * position z was appended as position N-1-z. */
  if (status == eslOK && max_drop != NULL) {
    P_idx         = P_cnt - 1;
    curr_P_score  = 0.0f;
    max_M_score   = -eslINFINITY;
    curr_max_drop = -eslINFINITY;

    for (z = 0; z < tr->N; z++) {
      /* Track the current P state score contribution so that it can be
       * subtracted out from M state scores to ensure the curr_max_drop
       * only looks for drops in M state scores, not splice penalties */
      if (tr->st[z] == p7T_P) {
        curr_P_score += P_scores[P_idx--];
      }
      else if (tr->st[z] == p7T_M) {
        curr_M_score = M_scores[tr->N - 1 - z] - curr_P_score;
        if (max_M_score != -eslINFINITY)
          curr_max_drop = ESL_MAX(curr_max_drop, max_M_score - curr_M_score);
        max_M_score = ESL_MAX(max_M_score, curr_M_score);
      }
    }
    *max_drop = curr_max_drop;
  }

  free(P_scores);
  free(M_scores);
  return status;

 ERROR:
  free(P_scores);
  free(M_scores);
  return status;
}

/*****************************************************************
 * 3. Benchmark driver.
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
  { "-T",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also benchmark Trace after each DP",   0 },
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
  P7_GMX         *gx           = NULL;
  P7_IVX         *acc_iv       = NULL;
  P7_IVX         *don_iv       = NULL;
  P7_TRACE        *tr          = NULL;
  int             N            = esl_opt_GetInteger(go, "-N");
  int             I            = esl_opt_GetInteger(go, "-I");
  int             intron_total = I + 4;
  int             do_T         = esl_opt_GetBoolean(go, "-T");
  int             min_intron   = 13;
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

  gx     = p7_gmx_Create(hmm->M, hmm->M, hmm->M, p7G_NSCELLS);
  acc_iv = p7_ivx_Create(hmm->M, SPLICE_ROWS);
  don_iv = p7_ivx_Create(hmm->M, SIGNAL_MEM_SIZE); 
  float signal_scores[p7S_SPLICE_SIGNALS];
  p7_SignalScores(signal_scores);
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
      p7_gmx_GrowTo(gx, hmm->M, L_dna_total, L_dna_total);

      p7_GViterbi_Spliced(dsq, gm_tr, gx, acc_iv, don_iv, signal_scores, 1, L_dna_total, 1, hmm->M, min_intron, TRUE, TRUE);
      final_C = gx->xmx[L_dna_total * p7G_NXCELLS + p7G_C];
      if (do_T && final_C != -eslINFINITY) {
        p7_GViterbi_SplicedTrace(dsq, gm_tr, gx, signal_scores, tr, 1, L_dna_total, 1, hmm->M, min_intron, NULL, NULL);
        p7_trace_Reuse(tr);
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
  p7_gmx_Destroy(gx);
  p7_ivx_Destroy(acc_iv);
  p7_ivx_Destroy(don_iv);
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
 * 4. Unit tests.
 *****************************************************************/
#ifdef p7GENERIC_VITERBI_SPLICED_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_compare_traces():
 *
 * Trace the filled matrix <gx> with p7_GViterbi_SplicedTraceReplay()
 * into <tr_replay> and require the same states, positions, splice codon
 * phases and Viterbi score as <tr_search>, which the caller traced with
 * p7_GViterbi_SplicedTrace(). Leaves both traces empty for reuse.
 */
static void
utest_compare_traces(ESL_DSQ *dsq, P7_FS_PROFILE *gm_tr, P7_GMX *gx, float *signal_scores,
                     P7_TRACE *tr_search, P7_TRACE *tr_replay, int L, int M, int min_intron, char *msg)
{
  float sc_search, sc_replay;
  int   z;

  p7_trace_Reuse(tr_search);
  p7_GViterbi_SplicedTrace      (dsq, gm_tr, gx, signal_scores, tr_search, 1, L, 1, M, min_intron, &sc_search, NULL);
  p7_GViterbi_SplicedTraceReplay(dsq, gm_tr, gx, signal_scores, tr_replay, 1, L, 1, M, min_intron, &sc_replay, NULL);

  if (tr_search->N != tr_replay->N) esl_fatal("%s: replay trace length %d != search trace length %d", msg, tr_replay->N, tr_search->N);
  for (z = 0; z < tr_search->N; z++)
    if (tr_search->st[z] != tr_replay->st[z] || tr_search->k[z] != tr_replay->k[z] ||
        tr_search->i[z]  != tr_replay->i[z]  || tr_search->c[z] != tr_replay->c[z])
      esl_fatal("%s: replay and search traces differ at position %d", msg, z);
  if (sc_search != sc_replay) esl_fatal("%s: replay score %f != search score %f", msg, sc_replay, sc_search);

  p7_trace_Reuse(tr_search);
  p7_trace_Reuse(tr_replay);
}

/* utest_compare_chk():
 *
 * Fill <dsq> with p7_GViterbi_SplicedChk() into a checkpointed matrix
 * with blocks of <blockR> rows, and require the same special states on
 * every row as the full matrix <gx>, which the caller filled with
 * p7_GViterbi_Spliced() and the same <global_start>, <global_end>. If
 * the final C state is reachable, trace both with
 * p7_GViterbi_SplicedTraceReplay() and p7_GViterbi_SplicedTraceChk()
 * and require the same trace, Viterbi score and max_drop.
 */
static void
utest_compare_chk(ESL_DSQ *dsq, P7_FS_PROFILE *gm_tr, P7_GMX *gx, float *signal_scores, int L, int M, int min_intron,
                  int global_start, int global_end, int blockR, char *msg)
{
  P7_CHKGMX *cx     = p7_chkgmx_Create(M, L, blockR, min_intron + 6, p7G_NSCELLS);
  P7_CHKIVX *cv     = p7_chkivx_Create(M, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, L / blockR);
  P7_IVX    *acc_iv = p7_ivx_Create(M, SPLICE_ROWS);
  P7_IVX    *don_iv = p7_ivx_Create(M, SIGNAL_MEM_SIZE);
  P7_TRACE  *tr_full = p7_trace_fs_Create();
  P7_TRACE  *tr_chk  = p7_trace_fs_Create();
  float      sc_full, sc_chk, drop_full, drop_chk;
  int        i, s, z;

  if (cx == NULL || cv == NULL || acc_iv == NULL || don_iv == NULL) esl_fatal("%s: allocation failed", msg);

  if (p7_GViterbi_SplicedChk(dsq, gm_tr, cx, cv, acc_iv, don_iv, signal_scores, 1, L, 1, M, min_intron, global_start, global_end) != eslOK)
    esl_fatal("%s: checkpointed fill failed", msg);

  for (i = 0; i <= L; i++)
    for (s = 0; s < p7G_NXCELLS; s++)
      if (cx->xmx[i*p7G_NXCELLS+s] != gx->xmx[i*p7G_NXCELLS+s])
        esl_fatal("%s: checkpointed fill differs at row %d, special state %d (block size %d)", msg, i, s, blockR);

  if (gx->xmx[L*p7G_NXCELLS+p7G_C] != -eslINFINITY)
    {
      if (p7_GViterbi_SplicedTraceReplay(dsq, gm_tr, gx, signal_scores, tr_full, 1, L, 1, M, min_intron, &sc_full, &drop_full) != eslOK)
        esl_fatal("%s: replay trace failed", msg);
      if (p7_GViterbi_SplicedTraceChk(dsq, gm_tr, cx, cv, acc_iv, don_iv, signal_scores, tr_chk, 1, L, 1, M, min_intron, global_start, global_end, &sc_chk, &drop_chk) != eslOK)
        esl_fatal("%s: checkpointed trace failed (block size %d)", msg, blockR);

      if (tr_full->N != tr_chk->N) esl_fatal("%s: checkpointed trace length %d != %d (block size %d)", msg, tr_chk->N, tr_full->N, blockR);
      for (z = 0; z < tr_full->N; z++)
        if (tr_full->st[z] != tr_chk->st[z] || tr_full->k[z] != tr_chk->k[z] ||
            tr_full->i[z]  != tr_chk->i[z]  || tr_full->c[z] != tr_chk->c[z])
          esl_fatal("%s: checkpointed trace differs at position %d (block size %d)", msg, z, blockR);
      if (sc_full   != sc_chk)   esl_fatal("%s: checkpointed trace score %f != %f", msg, sc_chk, sc_full);
      if (drop_full != drop_chk) esl_fatal("%s: checkpointed max_drop %f != %f",    msg, drop_chk, drop_full);
    }

  p7_trace_fs_Destroy(tr_full);
  p7_trace_fs_Destroy(tr_chk);
  p7_ivx_Destroy(acc_iv);
  p7_ivx_Destroy(don_iv);
  p7_chkivx_Destroy(cv);
  p7_chkgmx_Destroy(cx);
}

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
  P7_TRACE       *tr_NoP      = p7_trace_fs_Create();
  ESL_DSQ        *dsq         = NULL;
  P7_GMX         *gx          = NULL;
  P7_IVX         *acc_iv      = NULL;
  P7_IVX         *don_iv      = NULL;
  int             intron_total = intron_len + 4;  /* GT + intron_len random nucs + AG */
  int             min_intron   = 13;
  int             L_amino, L_dna_total;
  int             i, j, n_p;
  float           final_C;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig   (hmm, bgAA, gm,    M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, M, p7_UNILOCAL);

  gx     = p7_gmx_Create(hmm->M, hmm->M, hmm->M, p7G_NSCELLS);
  acc_iv = p7_ivx_Create(hmm->M, SPLICE_ROWS);
  don_iv = p7_ivx_Create(hmm->M, SIGNAL_MEM_SIZE);
  float signal_scores[p7S_SPLICE_SIGNALS];
  p7_SignalScores(signal_scores);

  while (N--)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, NULL);
      L_amino     = sq->n;
      if(L_amino < 4) continue;
      L_dna_total = L_amino * 3 + intron_total;
      
      /* Test can fail if there is a non-spliced path to L,M */
      if(L_dna_total %3 == 0) { intron_len ++; intron_total++; L_dna_total ++; }

      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) esl_fatal("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;

      /* Reverse-translate first half of the sequence (exon 1) */
      j = 1;
      for (i = 1; i <= L_amino / 2 ; i++) {
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
      p7_gmx_GrowTo(gx, M, L_dna_total, L_dna_total);

      p7_GViterbi_Spliced(dsq, gm_tr, gx, acc_iv, don_iv, signal_scores, 1, L_dna_total, 1, hmm->M, min_intron, TRUE, TRUE);
      final_C = gx->xmx[L_dna_total * p7G_NXCELLS + p7G_C];
      if (final_C != -eslINFINITY) {
        p7_GViterbi_SplicedTrace(dsq, gm_tr, gx, signal_scores, tr_NoP, 1, L_dna_total, 1, M, min_intron, NULL, NULL);
        n_p = 0;
        for (i = 0; i < tr_NoP->N; i++)
          if (tr_NoP->st[i] == p7T_P) n_p++;
        if (n_p < 1) esl_fatal(msg);
        utest_compare_traces(dsq, gm_tr, gx, signal_scores, tr_NoP, tr, L_dna_total, M, min_intron, msg);
      }
      utest_compare_chk(dsq, gm_tr, gx, signal_scores, L_dna_total, M, min_intron, TRUE, TRUE, min_intron + 6, msg);
      utest_compare_chk(dsq, gm_tr, gx, signal_scores, L_dna_total, M, min_intron, TRUE, TRUE, 64, msg);
      utest_compare_chk(dsq, gm_tr, gx, signal_scores, L_dna_total, M, min_intron, TRUE, TRUE, SPLICE_CHECKPOINT_ROWS, msg);

      /* --- Test 2: ExtendDown; final C must be reachable --- */
      p7_GViterbi_Spliced(dsq, gm_tr, gx, acc_iv, don_iv, signal_scores, 1, L_dna_total, 1, hmm->M, min_intron, TRUE, FALSE);
      final_C = gx->xmx[L_dna_total * p7G_NXCELLS + p7G_C];
      if (final_C == -eslINFINITY) esl_fatal(msg);
      utest_compare_traces(dsq, gm_tr, gx, signal_scores, tr_NoP, tr, L_dna_total, M, min_intron, msg);
      utest_compare_chk(dsq, gm_tr, gx, signal_scores, L_dna_total, M, min_intron, TRUE, FALSE, min_intron + 6, msg);
      utest_compare_chk(dsq, gm_tr, gx, signal_scores, L_dna_total, M, min_intron, TRUE, FALSE, 64, msg);

      /* --- Test 3: ExtendUp; final C must be reachable --- */
      p7_GViterbi_Spliced(dsq, gm_tr, gx, acc_iv, don_iv, signal_scores, 1, L_dna_total, 1, hmm->M, min_intron, FALSE, TRUE);
      final_C = gx->xmx[L_dna_total * p7G_NXCELLS + p7G_C];
      if (final_C == -eslINFINITY) esl_fatal(msg);
      utest_compare_traces(dsq, gm_tr, gx, signal_scores, tr_NoP, tr, L_dna_total, M, min_intron, msg);
      utest_compare_chk(dsq, gm_tr, gx, signal_scores, L_dna_total, M, min_intron, FALSE, TRUE, min_intron + 6, msg);
      utest_compare_chk(dsq, gm_tr, gx, signal_scores, L_dna_total, M, min_intron, FALSE, TRUE, 64, msg);
    }

  if (dsq != NULL) free(dsq);
  esl_sq_Destroy(sq);
  p7_trace_fs_Destroy(tr);
  p7_trace_fs_Destroy(tr_NoP);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_tr);
  p7_gmx_Destroy(gx);
  p7_ivx_Destroy(acc_iv);
  p7_ivx_Destroy(don_iv);

}

#endif /*p7GENERIC_VITERBI_SPLICED_TESTDRIVE*/
/*---------------------- end, unit tests ------------------------*/


/*****************************************************************
 * 5. Test driver.
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


