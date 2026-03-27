/* SSE-accelerated spliced Viterbi algorithm.
 *
 * Probability-space DP: _mm_mul_ps along a path, _mm_max_ps to select
 * the best path at each branch point. Sparse rescaling is used to
 * achieve full dynamic range on long sequences.
 *
 * The spliced Viterbi algorithm handles introns via a P (splice) state.
 * Donor site scores are accumulated into the OSPLICE_SCORES arrays
 * (OSS macros, striped __m128 vectors) during the forward pass.
 * Acceptor site scores are read back from those arrays to compute
 * P(i,k) at each sequence position.
 *
 * Differences from p7_Viterbi_Frameshift():
 *   - Single codon length (3-nt only): no IVX circular buffer.
 *     M(i,q) predecessors come from dpf[i-3] directly.
 *   - Global alignment: B enters only at i=3, k=1; E exits only at i=L, k=M.
 *   - P state stored in a 3-row circular buffer (pvbuf) alongside the
 *     main dp matrix.  P(i-3, k-1) is read from pvbuf via the same
 *     right-shift trick used for M/I/D.
 *   - Donor accumulation loop fills OSS arrays from dpf[i-min_intron-3].
 *
 * Contents:
 *   1. p7_Viterbi_SplicedGlobal()
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_sse.h"
#include "p7_splice.h"

#define MAX_NUC 4   /* same sentinel as generic: nucleotides 0..3 are valid */


/*****************************************************************
 * 1. p7_Viterbi_SplicedGlobal()
 *****************************************************************/

/* Function:  p7_Viterbi_SplicedGlobal()
 * Synopsis:  SSE-accelerated fully-global translated spliced Viterbi.
 *
 * Purpose:   SSE probability-space equivalent of
 *            p7_GViterbi_spliced_TranslatedGlobal().  Aligns globally
 *            (i_start..i_end on the nucleotide sequence, k_start..k_end
 *            on the profile) finding the maximum-scoring alignment that
 *            may cross one or more GT-AG/GC-AG/AT-AC splice junctions.
 *
 *            Donor site scores are accumulated into <oss> during the
 *            forward pass; acceptor lookups from <oss> compute the P
 *            (splice) state at each cell.  P state values are stored in
 *            a 3-row circular buffer (pvbuf) for use as predecessors
 *            three rows later.
 *
 *            The DP matrix <ox> must have been pre-grown to hold at
 *            least M model positions and L = i_end-i_start+1 sequence
 *            rows using p7_gmx_GrowTo().  The profile <om_fs> must have
 *            codon_lengths == 1 (single 3-nt codon length).
 *
 *            The final alignment score (in nats) is accessible as
 *              ox->totscale + logf(xC)
 *            where xC is stored in ox->xmx[L*p7X_NXCELLS + p7X_C].
 *
 * Args:      oss        - probability-space splice scores (OSS arrays)
 *            sub_dsq    - nucleotide sequence (1-indexed, i_start..i_end valid)
 *            om_fs      - optimized frameshift profile (codon_lengths == 1)
 *            ox         - DP matrix, pre-grown to (M, L, L, p7X_NSCELLS)
 *            i_start    - first nucleotide position (1-indexed in sub_dsq)
 *            i_end      - last  nucleotide position
 *            k_start    - first model position (1-indexed in om_fs)
 *            k_end      - last  model position
 *            min_intron - minimum intron length (nucleotides)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if om_fs->codon_lengths != 1.
 *            <eslEMEM>   on allocation failure.
 */
int
p7_Viterbi_SplicedGlobal(OSPLICE_SCORES *oss,
                         const ESL_DSQ *sub_dsq,
                         const P7_FS_OPROFILE *om_fs,
                         P7_OMX *ox,
                         int i_start, int i_end,
                         int k_start, int k_end,
                         int min_intron)
{
  int             Q          = p7O_NQF(om_fs->M);
  int             L          = i_end - i_start + 1;
  int             M          = k_end - k_start + 1;

  /* SSE registers */
  register __m128 mpv3, dpv3, ipv3, ppv3; /* right-shifted i-3 predecessors for M update */
  register __m128 sv;
  register __m128 msv;
  register __m128 pv;
  register __m128 dcv;
  register __m128 xEv;
  __m128          zerov      = _mm_setzero_ps();
  __m128          tsc_pv     = _mm_set1_ps(4.58e-5f);  /* P→M: exp(log(4.58e-5)) */
  __m128          sig_gtag_v, sig_gcag_v, sig_atac_v;
  __m128         *dpc, *dpp3;
  __m128         *tp;

  /* P-state circular buffer: 3 rows × Q stripes */
  __m128         *pvbuf     = NULL;
  __m128         *pvbuf_mem = NULL;
  __m128         *ppc;        /* current P-state row:     pvbuf[(i   %3)*Q] */
  __m128         *ppp3_row;   /* P-state row i-3:         pvbuf[((i-3)%3)*Q] */

  /* Scale correction for I(i,k) reading from dpp3 committed at S[i-3] */
  float           insert_adj;
  __m128          adj_v;

  /* Special state scalars */
  float           xB0;   /* N→B probability: B(row 0) */
  float           xE, xC;

  /* Nucleotide rolling windows */
  int             v, w, x;         /* codon window: C0 = codon(v,w,x) ending at i */
  int             r, s, t, u;      /* donor window: donor signal at (r,s) etc.     */
  int             acc0, acc1, acc2; /* acceptor signals delayed 0/1/2 rows          */
  int             donor0, donor1, donor2;
  int             C0, C1[4], C2[16];
  int             nuc1, nuc2;

  int             i, q, j, row, sub_i;
  int             status;

  if (om_fs->codon_lengths != 1)
    ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");

  ox->M              = M;
  ox->L              = L;
  ox->has_own_scales = TRUE;
  ox->totscale       = 0.0;

  /* Allocate P-state circular buffer */
  ESL_ALLOC(pvbuf_mem, sizeof(__m128) * 3 * Q + 15);
  pvbuf = (__m128 *)(((unsigned long int) pvbuf_mem + 15) & (~0xfUL));
  for (row = 0; row < 3 * Q; row++) pvbuf[row] = zerov;

  /* Zero-initialize DP matrix rows 0..L; set all scale factors to 1 */
  for (row = 0; row <= L; row++) {
    for (q = 0; q < Q; q++)
      MMO(ox->dpf[row], q) = DMO(ox->dpf[row], q) = IMO(ox->dpf[row], q) = zerov;
    ox->xmx[row * p7X_NXCELLS + p7X_SCALE] = 1.0f;
    ox->xmx[row * p7X_NXCELLS + p7X_E]     = 0.0f;
    ox->xmx[row * p7X_NXCELLS + p7X_C]     = 0.0f;
  }

  /* Reset OSS P-state accumulators to 0 (probability space: 0.0 = -inf) */
  for (q = 0; q < Q; q++)
    for (j = 0; j < SIGNAL_MEM_SIZE; j++)
      oss->oscore[q][j] = zerov;

  xB0        = om_fs->xf[p7O_N][p7O_MOVE];  /* N→B transition probability */
  sig_gtag_v = _mm_set1_ps(oss->signal_scores[p7S_GTAG]);
  sig_gcag_v = _mm_set1_ps(oss->signal_scores[p7S_GCAG]);
  sig_atac_v = _mm_set1_ps(oss->signal_scores[p7S_ATAC]);

  /* Initialize codon nucleotide window (v,w,x) and acceptor state */
  v = w = x = p7P_MAXNUC;
  acc0 = acc1 = acc2 = -1;

  /* Initialize donor nucleotide window (r,s,t,u): primed with positions 1,2,3 */
  r = p7P_MAXNUC;
  s = (i_start     <= i_end) ? (int) sub_dsq[i_start]     : p7P_MAXNUC;
  t = (i_start + 1 <= i_end) ? (int) sub_dsq[i_start + 1] : p7P_MAXNUC;
  u = (i_start + 2 <= i_end) ? (int) sub_dsq[i_start + 2] : p7P_MAXNUC;
  donor0 = donor1 = donor2 = -1;

  /*----------------------------------------------------------------
   * Main DP loop: i = 3..L
   *
   * Rows 0..2 are all zeros (no full 3-nt codons yet) and remain so.
   *----------------------------------------------------------------*/
  for (i = 3; i <= L; i++) {

    /* Advance codon window: v←w, w←x, x←sub_dsq[i] */
    sub_i = i_start + i - 1;
    v = w; w = x;
    x = (sub_dsq[sub_i] < MAX_NUC) ? (int) sub_dsq[sub_i] : p7P_MAXNUC;

    /* 3-nt codon index for current row */
    C0 = p7P_CODON3_FS1(v, w, x);
    C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);

    /* Splice codon indices for C1 (1 pre-donor nt) and C2 (2 pre-donor nt) */
    for (nuc1 = 0; nuc1 < 4; nuc1++) {
      C1[nuc1] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, w, x), p7P_DEGEN1_C);
      for (nuc2 = 0; nuc2 < 4; nuc2++)
        C2[nuc1*4 + nuc2] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, nuc2, x), p7P_DEGEN1_C);
    }

    /* Acceptor signal at (v,w): delay 0/1/2 rows for C0/C1/C2 splice codons */
    acc0 = acc1;
    acc1 = acc2;
    if      (v >= MAX_NUC || w >= MAX_NUC) acc2 = -1;
    else if (SIGNAL(v, w) == ACCEPT_AG)    acc2 = ACCEPT_AG;
    else if (SIGNAL(v, w) == ACCEPT_AC)    acc2 = ACCEPT_AC;
    else                                   acc2 = -1;

    /* Advance donor window and determine donor signals (only once we can look back) */
    if (i >= min_intron + 3) {
      r = s; s = t; t = u;
      u = (int) sub_dsq[i_start + i - min_intron];  /* nucleotide min_intron back from i */

      donor0 = donor1 = donor2 = -1;
      if (!(r >= MAX_NUC || s >= MAX_NUC)) {
        if      (SIGNAL(r, s) == DONOR_GT) donor0 = DONOR_GT;
        else if (SIGNAL(r, s) == DONOR_GC) donor0 = DONOR_GC;
        else if (SIGNAL(r, s) == DONOR_AT) donor0 = DONOR_AT;
        if (t < MAX_NUC) {
          if      (SIGNAL(s, t) == DONOR_GT) donor1 = DONOR_GT;
          else if (SIGNAL(s, t) == DONOR_GC) donor1 = DONOR_GC;
          else if (SIGNAL(s, t) == DONOR_AT) donor1 = DONOR_AT;
          if (u < MAX_NUC) {
            if      (SIGNAL(t, u) == DONOR_GT) donor2 = DONOR_GT;
            else if (SIGNAL(t, u) == DONOR_GC) donor2 = DONOR_GC;
            else if (SIGNAL(t, u) == DONOR_AT) donor2 = DONOR_AT;
          }
        }
      }
    }

    dpc      = ox->dpf[i];
    dpp3     = ox->dpf[i - 3];
    ppc      = pvbuf + ( i      % 3) * Q;
    ppp3_row = pvbuf + ((i - 3) % 3) * Q;

    /* Scale correction: I(i,k) reads dpp3 committed at cumulative scale S[i-3];
     * current running scale is S[i-1].  Correction = 1/(S[i-2]*S[i-1]). */
    insert_adj = 1.0f
               / (ox->xmx[(i-2) * p7X_NXCELLS + p7X_SCALE]
               *  ox->xmx[(i-1) * p7X_NXCELLS + p7X_SCALE]);
    adj_v = _mm_set1_ps(insert_adj);

    /* Right-shifted row-i-3 predecessors for M(i,k) ← (M/I/D/P)(i-3, k-1) */
    mpv3 = esl_sse_rightshiftz_float(MMO(dpp3,     Q-1));
    dpv3 = esl_sse_rightshiftz_float(DMO(dpp3,     Q-1));
    ipv3 = esl_sse_rightshiftz_float(IMO(dpp3,     Q-1));
    ppv3 = esl_sse_rightshiftz_float(ppp3_row[Q-1]);

    dcv  = zerov;
    xEv  = zerov;
    tp   = om_fs->tfv;

    /*-- Inner stripe loop: compute M, I, D, P for all model positions in row i --*/
    for (q = 0; q < Q; q++) {

      /* sv = max( M(i-3,k-1)*MM, I(i-3,k-1)*IM, D(i-3,k-1)*DM, P(i-3,k-1)*PM )
       * (BM transition skipped: global entry injected separately at i=3) */
      tp++;                                                              /* skip BM  */
      sv  =             _mm_mul_ps(mpv3, *tp); tp++;                    /* MM       */
      sv  = _mm_max_ps(sv, _mm_mul_ps(ipv3, *tp)); tp++;                /* IM       */
      sv  = _mm_max_ps(sv, _mm_mul_ps(dpv3, *tp)); tp++;                /* DM       */
      sv  = _mm_max_ps(sv, _mm_mul_ps(ppv3, tsc_pv));                   /* PM       */

      msv = _mm_mul_ps(sv, om_fs->rfv[C0][q]);
      xEv = _mm_max_ps(xEv, msv);

      /* Update right-shifted predecessors for next stripe (k→k+1) */
      mpv3 = MMO(dpp3, q);
      dpv3 = DMO(dpp3, q);
      ipv3 = IMO(dpp3, q);
      ppv3 = ppp3_row[q];

      MMO(dpc, q) = msv;
      DMO(dpc, q) = dcv;
      dcv = _mm_mul_ps(msv, *tp); tp++;                                  /* MD       */

      /* I(i,k) = max( M(i-3,k)*MI, I(i-3,k)*II ); mpv3/ipv3 now hold same-stripe values */
      sv          = _mm_mul_ps(_mm_mul_ps(mpv3, adj_v), *tp); tp++;     /* MI       */
      IMO(dpc, q) = _mm_max_ps(sv, _mm_mul_ps(_mm_mul_ps(ipv3, adj_v), *tp)); tp++;  /* II */

      /* P(i,k): acceptor-site lookup from accumulated donor scores in oss */
      pv = zerov;
      if (acc0 >= 0 || acc1 >= 0 || acc2 >= 0) {

        if (acc0 == ACCEPT_AG) {
          sv = _mm_max_ps(_mm_mul_ps(OSS0(oss, q, p7S_GTAG), sig_gtag_v),
                          _mm_mul_ps(OSS0(oss, q, p7S_GCAG), sig_gcag_v));
          pv = _mm_max_ps(pv, _mm_mul_ps(sv, om_fs->rfv[C0][q]));
        } else if (acc0 == ACCEPT_AC) {
          sv = _mm_mul_ps(OSS0(oss, q, p7S_ATAC), sig_atac_v);
          pv = _mm_max_ps(pv, _mm_mul_ps(sv, om_fs->rfv[C0][q]));
        }

        if (acc1 == ACCEPT_AG) {
          for (nuc1 = 0; nuc1 < 4; nuc1++) {
            sv = _mm_max_ps(_mm_mul_ps(OSS1(oss, q, p7S_GTAG, nuc1), sig_gtag_v),
                            _mm_mul_ps(OSS1(oss, q, p7S_GCAG, nuc1), sig_gcag_v));
            pv = _mm_max_ps(pv, _mm_mul_ps(sv, om_fs->rfv[C1[nuc1]][q]));
          }
        } else if (acc1 == ACCEPT_AC) {
          for (nuc1 = 0; nuc1 < 4; nuc1++) {
            sv = _mm_mul_ps(OSS1(oss, q, p7S_ATAC, nuc1), sig_atac_v);
            pv = _mm_max_ps(pv, _mm_mul_ps(sv, om_fs->rfv[C1[nuc1]][q]));
          }
        }

        if (acc2 == ACCEPT_AG) {
          for (nuc1 = 0; nuc1 < 4; nuc1++)
            for (nuc2 = 0; nuc2 < 4; nuc2++) {
              sv = _mm_max_ps(_mm_mul_ps(OSS2(oss, q, p7S_GTAG, nuc1, nuc2), sig_gtag_v),
                              _mm_mul_ps(OSS2(oss, q, p7S_GCAG, nuc1, nuc2), sig_gcag_v));
              pv = _mm_max_ps(pv, _mm_mul_ps(sv, om_fs->rfv[C2[nuc1*4+nuc2]][q]));
            }
        } else if (acc2 == ACCEPT_AC) {
          for (nuc1 = 0; nuc1 < 4; nuc1++)
            for (nuc2 = 0; nuc2 < 4; nuc2++) {
              sv = _mm_mul_ps(OSS2(oss, q, p7S_ATAC, nuc1, nuc2), sig_atac_v);
              pv = _mm_max_ps(pv, _mm_mul_ps(sv, om_fs->rfv[C2[nuc1*4+nuc2]][q]));
            }
        }
      }
      ppc[q] = pv;  /* Store P(i,k) for use as predecessor at row i+3 */

    } /* end inner q loop */

    /* Global entry at i=3: inject B→M(k=1) into element r=0 of stripe q=0.
     * k=1 is at (q=0, element r=0).  BM transition = 1.0 for global entry,
     * so M(3,k=1) = xB0 * rfv[C0][q=0][r=0]. */
    if (i == 3) {
      union { __m128 v; float p[4]; } u_entry, u_rfv;
      u_rfv.v      = om_fs->rfv[C0][0];
      u_entry.v    = MMO(dpc, 0);          /* currently zerov */
      u_entry.p[0] = ESL_MAX(u_entry.p[0], xB0 * u_rfv.p[0]);
      MMO(dpc, 0)  = u_entry.v;
    }

    /* DD paths: propagate D(i,k) ← max(M(i,k-1)*MD, D(i,k-1)*DD) along k */
    dcv        = esl_sse_rightshiftz_float(dcv);
    DMO(dpc,0) = zerov;
    tp         = om_fs->tfv + 7*Q;
    for (q = 0; q < Q; q++) {
      DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
      dcv         = _mm_mul_ps(DMO(dpc, q), *tp); tp++;
    }
    if (om_fs->M < 100) {
      for (j = 1; j < 4; j++) {
        dcv = esl_sse_rightshiftz_float(dcv);
        tp  = om_fs->tfv + 7*Q;
        for (q = 0; q < Q; q++) {
          DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
          dcv         = _mm_mul_ps(dcv, *tp); tp++;
        }
      }
    } else {
      for (j = 1; j < 4; j++) {
        register __m128 cv;
        dcv = esl_sse_rightshiftz_float(dcv);
        tp  = om_fs->tfv + 7*Q;
        cv  = zerov;
        for (q = 0; q < Q; q++) {
          sv          = _mm_max_ps(dcv, DMO(dpc, q));
          cv          = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO(dpc, q)));
          DMO(dpc, q) = sv;
          dcv         = _mm_mul_ps(dcv, *tp); tp++;
        }
        if (! _mm_movemask_ps(cv)) break;
      }
    }

    /* Donor accumulation: record max(M,D) at row i-min_intron-3 into OSS arrays.
     * This stores potential donor-site scores for acceptor lookup min_intron+3
     * rows later.  First possible P target is k=2, so zero out k=1
     * (element r=0 of stripe q=0) before accumulating. */
    if (i >= min_intron + 3 && (donor0 >= 0 || donor1 >= 0 || donor2 >= 0)) {
      __m128 *dpd = ox->dpf[i - min_intron - 3];
      for (q = 0; q < Q; q++) {
        __m128 tmp_sc = _mm_max_ps(MMO(dpd, q), DMO(dpd, q));
        if (q == 0) {
          /* Mask out k=1 (element r=0): P state requires k >= 2 */
          union { __m128 v; float p[4]; } u_tmp;
          u_tmp.v    = tmp_sc;
          u_tmp.p[0] = 0.0f;
          tmp_sc     = u_tmp.v;
        }
        if      (donor2 == DONOR_GT) OSS2(oss,q,p7S_GTAG,r,s) = _mm_max_ps(OSS2(oss,q,p7S_GTAG,r,s), tmp_sc);
        else if (donor2 == DONOR_GC) OSS2(oss,q,p7S_GCAG,r,s) = _mm_max_ps(OSS2(oss,q,p7S_GCAG,r,s), tmp_sc);
        else if (donor2 == DONOR_AT) OSS2(oss,q,p7S_ATAC,r,s) = _mm_max_ps(OSS2(oss,q,p7S_ATAC,r,s), tmp_sc);

        if      (donor1 == DONOR_GT) OSS1(oss,q,p7S_GTAG,r) = _mm_max_ps(OSS1(oss,q,p7S_GTAG,r), tmp_sc);
        else if (donor1 == DONOR_GC) OSS1(oss,q,p7S_GCAG,r) = _mm_max_ps(OSS1(oss,q,p7S_GCAG,r), tmp_sc);
        else if (donor1 == DONOR_AT) OSS1(oss,q,p7S_ATAC,r) = _mm_max_ps(OSS1(oss,q,p7S_ATAC,r), tmp_sc);

        if      (donor0 == DONOR_GT) OSS0(oss,q,p7S_GTAG) = _mm_max_ps(OSS0(oss,q,p7S_GTAG), tmp_sc);
        else if (donor0 == DONOR_GC) OSS0(oss,q,p7S_GCAG) = _mm_max_ps(OSS0(oss,q,p7S_GCAG), tmp_sc);
        else if (donor0 == DONOR_AT) OSS0(oss,q,p7S_ATAC) = _mm_max_ps(OSS0(oss,q,p7S_ATAC), tmp_sc);
      }
    }

    /* Accumulate D into xEv for rescaling trigger */
    for (q = 0; q < Q; q++)
      xEv = _mm_max_ps(xEv, DMO(dpc, q));
    esl_sse_hmax_ps(xEv, &xE);

    /* Sparse rescaling: scale row when max value exceeds threshold */
    if (xE > 1.0e4f) {
      float   scale_factor = 1.0f / xE;
      __m128  scale_v      = _mm_set1_ps(scale_factor);
      for (q = 0; q < Q; q++) {
        MMO(dpc, q) = _mm_mul_ps(MMO(dpc, q), scale_v);
        DMO(dpc, q) = _mm_mul_ps(DMO(dpc, q), scale_v);
        IMO(dpc, q) = _mm_mul_ps(IMO(dpc, q), scale_v);
        ppc[q]      = _mm_mul_ps(ppc[q],      scale_v);
      }
      ox->xmx[i * p7X_NXCELLS + p7X_SCALE] = xE;
      ox->totscale += logf(xE);
      /* Note: xmx[i][SCALE] was initialized to 1.0 in init loop above;
       * only written here if actual rescaling occurs. */
    }

  } /* end main loop i = 3..L */

  /* Global exit: E = max(M(L,k=M), D(L,k=M)) → C = E * T_EC
   * k=M is at stripe q_M = (M-1)%Q, element r_M = (M-1)/Q. */
  {
    int   q_M = (M - 1) % Q;
    int   r_M = (M - 1) / Q;
    union { __m128 v; float p[4]; } u_m, u_d;
    u_m.v = MMO(ox->dpf[L], q_M);
    u_d.v = DMO(ox->dpf[L], q_M);
    xE    = ESL_MAX(u_m.p[r_M], u_d.p[r_M]);
    xC    = xE * om_fs->xf[p7O_E][p7O_MOVE];
    ox->xmx[L * p7X_NXCELLS + p7X_E] = xE;
    ox->xmx[L * p7X_NXCELLS + p7X_C] = xC;
  }

  free(pvbuf_mem);
  return eslOK;

 ERROR:
  if (pvbuf_mem) free(pvbuf_mem);
  return status;
}
/*-------------- end, p7_Viterbi_SplicedGlobal() --------------*/
