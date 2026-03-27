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
 *   - P state stored in the OMX matrix using PMO_SP (p7X_NSCELLS_SP=4
 *     cells per stripe: M, D, I, P).  P(i-3, k-1) is read from the
 *     stored row via the same right-shift trick used for M/I/D.
 *   - Donor accumulation loop fills OSS arrays from dpf[i-min_intron-3].
 *
 * Contents:
 *   1. p7_Viterbi_SplicedGlobal()
 *   2. p7_Viterbi_SplicedTrace()
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Benchmark driver.
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
 *            (splice) state at each cell.  P state values are stored
 *            directly in the OMX matrix via PMO_SP() (the matrix must
 *            have been created with p7X_NSCELLS_SP cells per stripe).
 *
 *            The DP matrix <ox> must have been pre-grown to hold at
 *            least M model positions and L = i_end-i_start+1 sequence
 *            rows using p7_omx_GrowTo_dpf() with p7X_NSCELLS_SP.
 *            The profile <om_fs> must have codon_lengths == 1.
 *
 *            The final alignment score (in nats) is accessible as
 *              ox->totscale + logf(xC)
 *            where xC is stored in ox->xmx[L*p7X_NXCELLS + p7X_C].
 *
 * Args:      oss        - probability-space splice scores (OSS arrays)
 *            sub_dsq    - nucleotide sequence (1-indexed, i_start..i_end valid)
 *            om_fs      - optimized frameshift profile (codon_lengths == 1)
 *            ox         - DP matrix, pre-grown to (M, L, L, p7X_NSCELLS_SP)
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

  if (om_fs->codon_lengths != 1)
    ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");

  ox->M              = M;
  ox->L              = L;
  ox->has_own_scales = TRUE;
  ox->totscale       = 0.0;

  /* Zero-initialize DP matrix rows 0..L; set all scale factors to 1 */
  for (row = 0; row <= L; row++) {
    for (q = 0; q < Q; q++)
      MMO_SP(ox->dpf[row], q) = DMO_SP(ox->dpf[row], q) = IMO_SP(ox->dpf[row], q) = PMO_SP(ox->dpf[row], q) = zerov;
    ox->xmx[row * p7X_NXCELLS + p7X_SCALE] = 1.0f;
    ox->xmx[row * p7X_NXCELLS + p7X_E]     = 0.0f;
    ox->xmx[row * p7X_NXCELLS + p7X_N]     = 0.0f;
    ox->xmx[row * p7X_NXCELLS + p7X_J]     = 0.0f;
    ox->xmx[row * p7X_NXCELLS + p7X_B]     = 0.0f;
    ox->xmx[row * p7X_NXCELLS + p7X_C]     = 0.0f;
  }

  /* Reset OSS P-state accumulators to 0 (probability space: 0.0 = -inf) */
  for (q = 0; q < Q; q++)
    for (j = 0; j < SIGNAL_MEM_SIZE; j++)
      oss->oscore[q][j] = zerov;

  xB0        = om_fs->xf[p7O_N][p7O_MOVE];  /* N→B transition probability */
  ox->xmx[0 * p7X_NXCELLS + p7X_N] = 1.0f; /* S->N, p=1; store for dump/debug */
  ox->xmx[0 * p7X_NXCELLS + p7X_B] = xB0;  /* store for dump/debug */
  sig_gtag_v = _mm_set1_ps(oss->signal_scores[p7S_GTAG]);
  sig_gcag_v = _mm_set1_ps(oss->signal_scores[p7S_GCAG]);
  sig_atac_v = _mm_set1_ps(oss->signal_scores[p7S_ATAC]);

  /* Initialize codon nucleotide window (v,w,x) and acceptor state.
   * Pre-prime w and x with the first two nucleotides so that at i=3,
   * after v←w, w←x, x←sub_dsq[i_start+2], the first codon is valid.
   * This mirrors the generic's two-row pre-loop (i=1,2). */
  v   = p7P_MAXNUC;
  w   = ((int)sub_dsq[i_start]   < MAX_NUC) ? (int)sub_dsq[i_start]   : p7P_MAXNUC;
  x   = ((int)sub_dsq[i_start+1] < MAX_NUC) ? (int)sub_dsq[i_start+1] : p7P_MAXNUC;
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

    dpc  = ox->dpf[i];
    dpp3 = ox->dpf[i - 3];

    /* Right-shifted row-i-3 predecessors for M(i,k) ← (M/I/D/P)(i-3, k-1) */
    mpv3 = esl_sse_rightshiftz_float(MMO_SP(dpp3, Q-1));
    dpv3 = esl_sse_rightshiftz_float(DMO_SP(dpp3, Q-1));
    ipv3 = esl_sse_rightshiftz_float(IMO_SP(dpp3, Q-1));
    ppv3 = esl_sse_rightshiftz_float(PMO_SP(dpp3, Q-1));

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
      mpv3 = MMO_SP(dpp3, q);
      dpv3 = DMO_SP(dpp3, q);
      ipv3 = IMO_SP(dpp3, q);
      ppv3 = PMO_SP(dpp3, q);

      MMO_SP(dpc, q) = msv;
      DMO_SP(dpc, q) = dcv;
      dcv = _mm_mul_ps(msv, *tp); tp++;                                  /* MD       */

      /* I(i,k) = max( M(i-3,k)*MI, I(i-3,k)*II ); mpv3/ipv3 now hold same-stripe values */
      sv             = _mm_mul_ps(mpv3, *tp); tp++;  /* MI */
      IMO_SP(dpc, q) = _mm_max_ps(sv, _mm_mul_ps(ipv3, *tp)); tp++;  /* II */
      /* Zero I at stop-codon positions (matches generic: if rsc=-inf then I=-inf) */
      IMO_SP(dpc, q) = _mm_and_ps(IMO_SP(dpc, q), _mm_cmpgt_ps(om_fs->rfv[C0][q], zerov));

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
      PMO_SP(dpc, q) = pv;  /* Store P(i,k) in OMX for use as predecessor at row i+3 */

    } /* end inner q loop */

    /* Global entry at i=3: inject B→M(k=1) into element r=0 of stripe q=0.
     * k=1 is at (q=0, element r=0).  BM transition = 1.0 for global entry,
     * so M(3,k=1) = xB0 * rfv[C0][q=0][r=0]. */
    if (i == 3) {
      union { __m128 v; float p[4]; } u_entry, u_rfv;
      u_rfv.v      = om_fs->rfv[C0][0];
      u_entry.v    = MMO_SP(dpc, 0);        /* currently zerov */
      u_entry.p[0] = ESL_MAX(u_entry.p[0], xB0 * u_rfv.p[0]);
      MMO_SP(dpc, 0) = u_entry.v;
    }

    /* DD paths: propagate D(i,k) ← max(M(i,k-1)*MD, D(i,k-1)*DD) along k */
    dcv           = esl_sse_rightshiftz_float(dcv);
    DMO_SP(dpc,0) = zerov;
    tp            = om_fs->tfv + 7*Q;
    for (q = 0; q < Q; q++) {
      DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
      dcv            = _mm_mul_ps(DMO_SP(dpc, q), *tp); tp++;
    }
    if (om_fs->M < 100) {
      for (j = 1; j < 4; j++) {
        dcv = esl_sse_rightshiftz_float(dcv);
        tp  = om_fs->tfv + 7*Q;
        for (q = 0; q < Q; q++) {
          DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
          dcv            = _mm_mul_ps(dcv, *tp); tp++;
        }
      }
    } else {
      for (j = 1; j < 4; j++) {
        register __m128 cv;
        dcv = esl_sse_rightshiftz_float(dcv);
        tp  = om_fs->tfv + 7*Q;
        cv  = zerov;
        for (q = 0; q < Q; q++) {
          sv             = _mm_max_ps(dcv, DMO_SP(dpc, q));
          cv             = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO_SP(dpc, q)));
          DMO_SP(dpc, q) = sv;
          dcv            = _mm_mul_ps(dcv, *tp); tp++;
        }
        if (! _mm_movemask_ps(cv)) break;
      }
    }

    /* At i=3, the DD fixup above used dcv=0 (dpp3=row 0, all zero), so
     * D(3,k≥2) from the global-entry M(3,k=1) was not propagated.
     * Re-run DD propagation seeded by M(3,k=1)*MD. */
    if (i == 3) {
      dcv = _mm_mul_ps(MMO_SP(dpc, 0), om_fs->tfv[p7O_MD]);  /* M(k)*MD for q=0 → D(k+1) in q=1 */
      tp  = om_fs->tfv + 7*Q + 1;                             /* DD for q=1..Q-1 */
      for (q = 1; q < Q; q++) {
        DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
        dcv = _mm_mul_ps(DMO_SP(dpc, q), *tp); tp++;
      }
      if (om_fs->M < 100) {
        for (j = 1; j < 4; j++) {
          dcv = esl_sse_rightshiftz_float(dcv);
          tp  = om_fs->tfv + 7*Q;
          for (q = 0; q < Q; q++) {
            DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
            dcv = _mm_mul_ps(DMO_SP(dpc, q), *tp); tp++;
          }
        }
      } else {
        for (j = 1; j < 4; j++) {
          register __m128 cv;
          dcv = esl_sse_rightshiftz_float(dcv);
          tp  = om_fs->tfv + 7*Q;
          cv  = zerov;
          for (q = 0; q < Q; q++) {
            sv             = _mm_max_ps(dcv, DMO_SP(dpc, q));
            cv             = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO_SP(dpc, q)));
            DMO_SP(dpc, q) = sv;
            dcv            = _mm_mul_ps(dcv, *tp); tp++;
          }
          if (! _mm_movemask_ps(cv)) break;
        }
      }
    }

    /* Donor accumulation: record max(M(j,k-1),D(j,k-1)) at row j=i-min_intron-3
     * into OSS[q] for each stripe q.  The generic stores max(M,D) at k-1 under
     * index k, so we right-shift the donor vector by one stripe position using the
     * same rolling trick used for M predecessors.  The rightshiftz inserts 0 at
     * stripe 0 element 0, masking k=1 (P state requires k >= 2) for free. */
    if (i >= min_intron + 3 && (donor0 >= 0 || donor1 >= 0 || donor2 >= 0)) {
      __m128 *dpd      = ox->dpf[i - min_intron - 3];
      __m128  donor_prev = esl_sse_rightshiftz_float(_mm_max_ps(MMO_SP(dpd, Q-1), DMO_SP(dpd, Q-1)));
      for (q = 0; q < Q; q++) {
        __m128 tmp_sc = donor_prev;
        donor_prev    = _mm_max_ps(MMO_SP(dpd, q), DMO_SP(dpd, q));

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
      xEv = _mm_max_ps(xEv, DMO_SP(dpc, q));
    esl_sse_hmax_ps(xEv, &xE);

    /* Sparse rescaling: scale row when max value exceeds threshold */
    if (xE > 1.0e4f) {
      float   scale_factor = 1.0f / xE;
      __m128  scale_v      = _mm_set1_ps(scale_factor);
      for (q = 0; q < Q; q++) {
        MMO_SP(dpc, q) = _mm_mul_ps(MMO_SP(dpc, q), scale_v);
        DMO_SP(dpc, q) = _mm_mul_ps(DMO_SP(dpc, q), scale_v);
        IMO_SP(dpc, q) = _mm_mul_ps(IMO_SP(dpc, q), scale_v);
        PMO_SP(dpc, q) = _mm_mul_ps(PMO_SP(dpc, q), scale_v);
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
    u_m.v = MMO_SP(ox->dpf[L], q_M);
    u_d.v = DMO_SP(ox->dpf[L], q_M);
    xE    = ESL_MAX(u_m.p[r_M], u_d.p[r_M]);
    xC    = xE * om_fs->xf[p7O_E][p7O_MOVE];
    ox->xmx[L * p7X_NXCELLS + p7X_E] = xE;
    ox->xmx[L * p7X_NXCELLS + p7X_C] = xC;
  }

  return eslOK;
}
/*-------------- end, p7_Viterbi_SplicedGlobal() --------------*/


/*****************************************************************
 * 2. p7_Viterbi_SplicedTrace()
 *****************************************************************/

/* Forward declarations for vit_select_*_sp() helpers */
static inline float get_msc_sp(const P7_OMX *ox, int i, int k, int Q);
static inline float get_dsc_sp(const P7_OMX *ox, int i, int k, int Q);
static inline float get_isc_sp(const P7_OMX *ox, int i, int k, int Q);
static inline float get_psc_sp(const P7_OMX *ox, int i, int k, int Q);
static inline float get_rsc_sp(const P7_FS_OPROFILE *om_fs, int cn, int k, int Q);
static inline float get_tsc_sp(const P7_FS_OPROFILE *om_fs, int t_off, int k, int Q);
static inline float get_dd_sp (const P7_FS_OPROFILE *om_fs, int k, int Q);
static inline int   vit_select_m_sp(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int M_loc, int i, int k, int Q);
static inline int   vit_select_d_sp(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k, int Q);
static inline int   vit_select_i_sp(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k, int Q);
static inline int   vit_select_e_sp(const P7_OMX *ox, int i, int M_loc, int Q, int *ret_k);


/*****************************************************************
 * 2a. vit_select_*_sp() helper functions
 *****************************************************************/

/* Scalar M/D/I/P cell extractors: strip r=(k-1)/Q from stripe q=(k-1)%Q.
 * Must use the _SP macros because ox was created with p7X_NSCELLS_SP stride. */
static inline float
get_msc_sp(const P7_OMX *ox, int i, int k, int Q)
{ union { __m128 v; float p[4]; } u; u.v = MMO_SP(ox->dpf[i], (k-1)%Q); return u.p[(k-1)/Q]; }

static inline float
get_dsc_sp(const P7_OMX *ox, int i, int k, int Q)
{ union { __m128 v; float p[4]; } u; u.v = DMO_SP(ox->dpf[i], (k-1)%Q); return u.p[(k-1)/Q]; }

static inline float
get_isc_sp(const P7_OMX *ox, int i, int k, int Q)
{ union { __m128 v; float p[4]; } u; u.v = IMO_SP(ox->dpf[i], (k-1)%Q); return u.p[(k-1)/Q]; }

static inline float
get_psc_sp(const P7_OMX *ox, int i, int k, int Q)
{ union { __m128 v; float p[4]; } u; u.v = PMO_SP(ox->dpf[i], (k-1)%Q); return u.p[(k-1)/Q]; }

/* Emission probability for codon index cn at model position k */
static inline float
get_rsc_sp(const P7_FS_OPROFILE *om_fs, int cn, int k, int Q)
{ union { __m128 v; float p[4]; } u; u.v = om_fs->rfv[cn][(k-1)%Q]; return u.p[(k-1)/Q]; }

/* Main transition: t_off 0=BM 1=MM 2=IM 3=DM 4=MD 5=MI 6=II; layout tfv[7*q + t_off] */
static inline float
get_tsc_sp(const P7_FS_OPROFILE *om_fs, int t_off, int k, int Q)
{ union { __m128 v; float p[4]; } u; u.v = om_fs->tfv[7 * ((k-1)%Q) + t_off]; return u.p[(k-1)/Q]; }

/* DD transition: layout tfv[7*Q + q] */
static inline float
get_dd_sp(const P7_FS_OPROFILE *om_fs, int k, int Q)
{ union { __m128 v; float p[4]; } u; u.v = om_fs->tfv[7*Q + (k-1)%Q]; return u.p[(k-1)/Q]; }


/* M(i,k) predecessors: M/I/D/P from row i-3 at k-1, or B (global entry at i=3, k=1).
 * P state is read directly from PMO_SP in the OMX matrix. */
static inline int
vit_select_m_sp(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox,
                int M_loc, int i, int k, int Q)
{
  float path[5];
  int   state[5] = { p7T_M, p7T_I, p7T_D, p7T_P, p7T_B };

  if (k > 1) {
    path[0] = get_msc_sp(ox, i-3, k-1, Q) * get_tsc_sp(om_fs, 1, k, Q);  /* M * MM */
    path[1] = get_isc_sp(ox, i-3, k-1, Q) * get_tsc_sp(om_fs, 2, k, Q);  /* I * IM */
    path[2] = get_dsc_sp(ox, i-3, k-1, Q) * get_tsc_sp(om_fs, 3, k, Q);  /* D * DM */
    path[3] = get_psc_sp(ox, i-3, k-1, Q) * 4.58e-5f;                     /* P * PM */
  } else {
    path[0] = path[1] = path[2] = path[3] = 0.0f;
  }
  /* Global entry B→M only at i=3, k=1. */
  path[4] = (k == 1 && i == 3) ? om_fs->xf[p7O_N][p7O_MOVE] : 0.0f;

  return state[esl_vec_FArgMax(path, 5)];
}

/* D(i,k) predecessors M(i,k-1) and D(i,k-1) are at the same row; no scale correction. */
static inline int
vit_select_d_sp(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k, int Q)
{
  float path[2];
  int   state[2] = { p7T_M, p7T_D };

  if (k > 1) {
    path[0] = get_msc_sp(ox, i, k-1, Q) * get_tsc_sp(om_fs, 4, k, Q);  /* M * MD */
    path[1] = get_dsc_sp(ox, i, k-1, Q) * get_dd_sp (om_fs, k, Q);     /* D * DD */
  } else {
    path[0] = path[1] = 0.0f;
  }
  return state[esl_vec_FArgMax(path, 2)];
}

/* I(i,k) predecessors M(i-3,k) and I(i-3,k) are both at row i-3; no relative correction. */
static inline int
vit_select_i_sp(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k, int Q)
{
  float path[2];
  int   state[2] = { p7T_M, p7T_I };

  path[0] = get_msc_sp(ox, i-3, k, Q) * get_tsc_sp(om_fs, 5, k, Q);  /* M * MI */
  path[1] = get_isc_sp(ox, i-3, k, Q) * get_tsc_sp(om_fs, 6, k, Q);  /* I * II */
  return state[esl_vec_FArgMax(path, 2)];
}

/* Global exit: E is reached only from M(L,M) or D(L,M). */
static inline int
vit_select_e_sp(const P7_OMX *ox, int i, int M_loc, int Q, int *ret_k)
{
  float mval = get_msc_sp(ox, i, M_loc, Q);
  float dval = get_dsc_sp(ox, i, M_loc, Q);
  *ret_k = M_loc;
  return (mval >= dval) ? p7T_M : p7T_D;
}
/*-------------- end, vit_select_*_sp() helpers -----------------*/


/* Function:  p7_Viterbi_SplicedTrace()
 * Synopsis:  Traceback from a spliced Viterbi DP matrix; SSE probability-space version.
 *
 * Purpose:   Extract the optimal alignment path from the DP matrix <ox> filled
 *            by p7_Viterbi_SplicedGlobal().  P-state predecessors of M cells
 *            are read directly from PMO_SP() in the OMX matrix (the matrix
 *            must have been created with p7X_NSCELLS_SP cells per stripe).
 *            Donor site search uses per-row cumulative scale factors stored
 *            in ox->xmx to normalize candidates across different DP rows
 *            before argmax comparison.
 *            <signal_scores> is an array of p7S_SPLICE_SIGNALS probabilities
 *            (GT-AG, GC-AG, AT-AC) as filled by p7_osplicescores_Create().
 *
 * Args:      sub_dsq       - nucleotide sequence (1-indexed)
 *            om_fs         - optimized frameshift profile (codon_lengths == 1)
 *            ox            - DP matrix from p7_Viterbi_SplicedGlobal()
 *            signal_scores - splice signal probabilities [p7S_GTAG/GCAG/ATAC]
 *            tr            - RETURN: traceback (caller provides empty trace)
 *            i_start       - first nucleotide position
 *            i_end         - last  nucleotide position
 *            k_start       - first model position (must be 1 for current impl.)
 *            k_end         - last  model position
 *            min_intron    - minimum intron length
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if traceback reaches an impossible state.
 *            <eslEMEM>   on allocation failure.
 */
int
p7_Viterbi_SplicedTrace(const ESL_DSQ *sub_dsq,
                        const P7_FS_OPROFILE *om_fs,
                        const P7_OMX *ox,
                        const float *signal_scores,
                        P7_TRACE *tr,
                        int i_start, int i_end,
                        int k_start, int k_end,
                        int min_intron)
{
  int     Q         = p7O_NQF(om_fs->M);
  int     L         = i_end - i_start + 1;
  int     M_loc     = k_end - k_start + 1;
  int     i         = L;
  int     k         = 0;
  int     c         = 0;
  int     j         = 0;     /* donor DP row, set during P donor search       */
  int     snxt      = p7T_X; /* state following donor, set during P search    */
  int     s0, s1;
  int     status;

  /* Per-row log-cumulative scale for donor-site candidate normalization.
   * logcumscale[r] = sum_{j=1}^{r} log( scale_factor[j] ).
   * Ratio cumscale(a)/cumscale(b) = exp( logcumscale[a] - logcumscale[b] ). */
  float  *logcumscale = NULL;
  int     r;
  ESL_ALLOC(logcumscale, (L + 2) * sizeof(float));
  logcumscale[0] = 0.0f;
  for (r = 1; r <= L; r++)
    logcumscale[r] = logcumscale[r-1] + logf(ox->xmx[r * p7X_NXCELLS + p7X_SCALE]);

  if (ox->xmx[L * p7X_NXCELLS + p7X_C] == 0.0f)
    ESL_EXCEPTION(eslFAIL, "impossible C reached at L=%d: alignment score is -inf", L);

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i, 0)) != eslOK) goto ERROR;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i, 0)) != eslOK) goto ERROR;
  s0 = p7T_C;

  while (s0 != p7T_S) {
    switch (s0) {

    case p7T_C:                       /* global: C always comes directly from E */
      s1 = p7T_E;
      break;

    case p7T_E:                       /* global exit only from M(L,M) or D(L,M) */
      s1 = vit_select_e_sp(ox, i, M_loc, Q, &k);
      break;

    case p7T_M:
      s1 = vit_select_m_sp(om_fs, ox, M_loc, i, k, Q);
      k--; i -= 3;
      break;

    case p7T_D:
      s1 = vit_select_d_sp(om_fs, ox, i, k, Q);
      k--;
      break;

    case p7T_I:
      s1 = vit_select_i_sp(om_fs, ox, i, k, Q);
      i -= 3;
      break;

    case p7T_P:                       /* P: jump to donor exon predecessor */
      s1 = snxt;
      k--;
      i  = j - c - 2;
      break;

    case p7T_B:                       /* global: B always comes from N */
      s1 = p7T_N;
      break;

    case p7T_N:
      s1 = (i == 0) ? p7T_S : p7T_N;
      break;

    default: ESL_EXCEPTION(eslEINVAL, "bogus state in Viterbi_SplicedTrace");
    }

    /* Donor site search: when s1==P, (i,k) is already the P state position
     * (updated by case p7T_M above).  Scan backward over j to find which
     * donor signal and codon offset c gave the best score into P(i,k). */
    if (s1 == p7T_P) {
      int   sub_i  = i_start + i - 1;
      int   v      = (sub_i >= 2 && sub_dsq[sub_i-2] < MAX_NUC) ? (int)sub_dsq[sub_i-2] : p7P_MAXCODONS1;
      int   w      = (sub_i >= 1 && sub_dsq[sub_i-1] < MAX_NUC) ? (int)sub_dsq[sub_i-1] : p7P_MAXCODONS1;
      int   x      = (              sub_dsq[sub_i  ] < MAX_NUC) ? (int)sub_dsq[sub_i  ] : p7P_MAXCODONS1;
      int   c0_idx = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);
      float emit0  = get_rsc_sp(om_fs, c0_idx, k, Q);

      float best   = -1.0f;
      int   best_j = -1;
      int   jj, cc;
      snxt = p7T_M;

      for (jj = i - min_intron + 1; jj > 3; jj--) {
        int   sub_d  = i_start + jj - 1;
        int   dn0    = (int)sub_dsq[sub_d - 1];
        int   dn1    = (int)sub_dsq[sub_d    ];
        float sig;

        if (dn0 >= MAX_NUC || dn1 >= MAX_NUC) continue;
        if      (SIGNAL(dn0, dn1) == DONOR_GT) sig = signal_scores[p7S_GTAG];
        else if (SIGNAL(dn0, dn1) == DONOR_GC) sig = signal_scores[p7S_GCAG];
        else if (SIGNAL(dn0, dn1) == DONOR_AT) sig = signal_scores[p7S_ATAC];
        else continue;

        /* Nucleotides just before the donor for c=1 (1 nt) and c=2 (2 nt) variants */
        int tu = (sub_d >= 2 && sub_dsq[sub_d-2] < MAX_NUC) ? (int)sub_dsq[sub_d-2] : p7P_MAXCODONS1;
        int tt = (sub_d >= 3 && sub_dsq[sub_d-3] < MAX_NUC) ? (int)sub_dsq[sub_d-3] : p7P_MAXCODONS1;
        int c1_idx = p7P_MINIDX(p7P_CODON3_FS1(tu, w, x), p7P_DEGEN1_C);
        int c2_idx = p7P_MINIDX(p7P_CODON3_FS1(tt, tu, x), p7P_DEGEN1_C);
        float emit1 = get_rsc_sp(om_fs, c1_idx, k, Q);
        float emit2 = get_rsc_sp(om_fs, c2_idx, k, Q);

        /* Candidates at three codon offsets; normalize to scale of row i via
         * adj = exp( logcumscale[row] - logcumscale[i] ). */
        float adj, cand;
        float mval, dval;

        /* c=0: predecessor DP row = jj-2 */
        if (jj >= 2) {
          adj  = expf(logcumscale[jj-2] - logcumscale[i]);
          mval = get_msc_sp(ox, jj-2, k-1, Q) * sig * emit0 * adj;
          dval = get_dsc_sp(ox, jj-2, k-1, Q) * sig * emit0 * adj;
          cand = ESL_MAX(mval, dval);
          if (cand > best) { best = cand; best_j = jj; c = 0; snxt = (mval >= dval) ? p7T_M : p7T_D; }
        }
        /* c=1: predecessor DP row = jj-3 */
        if (jj >= 3) {
          adj  = expf(logcumscale[jj-3] - logcumscale[i]);
          mval = get_msc_sp(ox, jj-3, k-1, Q) * sig * emit1 * adj;
          dval = get_dsc_sp(ox, jj-3, k-1, Q) * sig * emit1 * adj;
          cand = ESL_MAX(mval, dval);
          if (cand > best) { best = cand; best_j = jj; c = 1; snxt = (mval >= dval) ? p7T_M : p7T_D; }
        }
        /* c=2: predecessor DP row = jj-4 */
        if (jj >= 4) {
          adj  = expf(logcumscale[jj-4] - logcumscale[i]);
          mval = get_msc_sp(ox, jj-4, k-1, Q) * sig * emit2 * adj;
          dval = get_dsc_sp(ox, jj-4, k-1, Q) * sig * emit2 * adj;
          cand = ESL_MAX(mval, dval);
          if (cand > best) { best = cand; best_j = jj; c = 2; snxt = (mval >= dval) ? p7T_M : p7T_D; }
        }
      }
      if (best_j == -1) ESL_EXCEPTION(eslEINVAL, "P at k=%d,i=%d couldn't be traced", k, i);
      j = best_j;
    } /* end donor search */

    /* Codon length: M states always use c=3; P state uses c from donor search. */
    if      (s1 == p7T_M) c = 3;
    else if (s1 != p7T_P) c = 0;

    if ((status = p7_trace_fs_Append(tr, s1, k_start + k - 1, i_start + i - 1, c)) != eslOK) goto ERROR;

    /* N/C: deferred i decrement for self-transitions */
    if ((s1 == p7T_N || s1 == p7T_C) && s1 == s0) i--;

    s0 = s1;
  } /* end traceback loop */

  tr->M = M_loc;
  tr->L = L;

  free(logcumscale);
  return p7_trace_fs_Reverse(tr);

 ERROR:
  if (logcumscale) free(logcumscale);
  return status;
}
/*-------------- end, p7_Viterbi_SplicedTrace() --------------*/


/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef p7VITERBI_SP_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_spliced():
 *
 * Build the same protein-emitted, reverse-translated DNA sequence with a
 * simulated GT..AG intron that the generic unit test uses.  For each sequence:
 *
 * 1. Run p7_GViterbi_spliced_TranslatedGlobal (generic) and read the final
 *    C-state log-probability score.
 * 2. Run p7_Viterbi_SplicedGlobal (SSE) and convert its probability-space
 *    C-state score to log-probability; compare against the generic score.
 * 3. Run p7_Viterbi_SplicedTrace (SSE) and verify the trace contains at
 *    least one P state, confirming the splice junction was detected.
 */
static void
utest_spliced(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
              ESL_GENCODE *gcode, P7_BG *bgAA, P7_CODONTABLE *codon_table,
              int M, int N, int intron_len)
{
  char           *msg          = "viterbi_sp utest_spliced failed";
  P7_HMM         *hmm          = NULL;
  P7_PROFILE     *gm           = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_tr        = p7_profile_fs_Create(M, abcAA, 1);
  P7_FS_OPROFILE *om_fs        = p7_fs_oprofile_Create(M, abcAA, 1);
  ESL_SQ         *sq           = esl_sq_CreateDigital(abcAA);
  P7_TRACE       *tr_sse       = p7_trace_fs_Create();
  ESL_DSQ        *dsq          = NULL;
  SPLICE_PIPELINE *pli         = NULL;
  P7_OMX         *ox           = p7_omx_Create_dpf(M, M*3, M*3, p7X_NSCELLS_SP);
  OSPLICE_SCORES *oss          = p7_osplicescores_Create(M);
  int             intron_total = intron_len + 4;   /* GT + intron_len + AG */
  int             L_amino, L_dna_total;
  int             i, j, n_p;
  float           ref_sc, sse_sc;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig   (hmm, bgAA, gm,    M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, M, p7_UNILOCAL);
  p7_fs_oprofile_Convert(gm_tr, om_fs);

  pli = p7_splicepipeline_Create(NULL, M, M * 3);

  while (N--)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, NULL);
      L_amino     = sq->n;
      L_dna_total = L_amino * 3 + intron_total;

      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) esl_fatal(msg);
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;

      /* Reverse-translate first exon */
      j = 1;
      for (i = 1; i <= L_amino / 2; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }
      /* Simulated intron: GT + intron_len random nucleotides + AG */
      dsq[j++] = 2;  /* G */
      dsq[j++] = 3;  /* T */
      for (i = 0; i < intron_len; i++) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = 0;  /* A */
      dsq[j++] = 2;  /* G */
      /* Reverse-translate second exon */
      for (i = L_amino / 2 + 1; i <= L_amino; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }

      /* Reconfigure lengths */
      p7_fs_ReconfigLength(gm_tr, L_dna_total / 3);
      p7_fs_oprofile_ReconfigLength(om_fs, L_dna_total / 3);

      /* Grow matrices */
      p7_gmx_GrowTo(pli->vit, M, L_dna_total, L_dna_total);
      p7_omx_GrowTo_dpf(ox, M, L_dna_total, L_dna_total);
      p7_splicescores_GrowTo(pli->splice_scores, M);
      p7_osplicescores_GrowTo(oss, M);

      /* --- Test 1: compare final C-state score (generic log-prob vs SSE) --- */
      p7_GViterbi_spliced_TranslatedGlobal(pli, dsq, gm_tr, pli->vit, 1, L_dna_total, 1, M);
      ref_sc = pli->vit->xmx[L_dna_total * p7G_NXCELLS + p7G_C];

      p7_Viterbi_SplicedGlobal(oss, dsq, om_fs, ox, 1, L_dna_total, 1, M, pli->min_intron);
      sse_sc = ox->totscale + logf(ox->xmx[L_dna_total * p7X_NXCELLS + p7X_C]);
      printf("N %d sse_sc %f ref_sc %f\n", N, sse_sc, ref_sc);
      if (fabsf(sse_sc - ref_sc) > 0.001f) {
        printf("Score mismatch: dumping SSE matrix (logified):\n");
        p7_omx_DumpSP(stdout, ox, 0, L_dna_total, 0, M, TRUE);
        printf("Generic matrix:\n");
        p7_gmx_DumpWindow(stdout, pli->vit, 0, L_dna_total, 0, M, p7_DEFAULT);
        esl_fatal(msg);
      }

      /* --- Test 2: SSE trace must contain at least one P state ---
       * Skip if score is -inf (impossible alignment: no valid trace exists). */
      if (!isinf(sse_sc)) {
        p7_trace_Reuse(tr_sse);
        p7_Viterbi_SplicedTrace(dsq, om_fs, ox, oss->signal_scores,
                                tr_sse, 1, L_dna_total, 1, M, pli->min_intron);
        n_p = 0;
        for (i = 0; i < tr_sse->N; i++)
          if (tr_sse->st[i] == p7T_P) n_p++;
        if (n_p < 1) esl_fatal(msg);
      }
    }

  if (dsq != NULL) free(dsq);
  esl_sq_Destroy(sq);
  p7_trace_fs_Destroy(tr_sse);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_tr);
  p7_fs_oprofile_Destroy(om_fs);
  p7_omx_Destroy(ox);
  p7_osplicescores_Destroy(oss);
  p7_splicepipeline_Destroy(pli);
}
#endif /*p7VITERBI_SP_TESTDRIVE*/
/*---------------------- end, unit tests ------------------------*/


/*****************************************************************
 * 4. Test driver.
 *****************************************************************/
#ifdef p7VITERBI_SP_TESTDRIVE
/*
   gcc -g -Wall -std=gnu99 -o viterbi_sp_utest -I. -I.. -L.. -I../../easel -L../../easel -Dp7VITERBI_SP_TESTDRIVE viterbi_sp.c -lhmmer -leasel -lm
   ./viterbi_sp_utest
 */
#include "esl_gencode.h"
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name    type          default env range toggles reqs incomp help                                  docgroup */
  { "-h", eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",        0 },
  { "-s", eslARG_INT,    "42", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",               0 },
  { "-M", eslARG_INT,   "100", NULL, NULL, NULL, NULL, NULL, "size of random models to sample",             0 },
  { "-N", eslARG_INT,    "20", NULL, NULL, NULL, NULL, NULL, "number of random sequences to sample",        0 },
  { "-I", eslARG_INT,   "500", NULL, "n>0",NULL, NULL, NULL, "simulated intron length (random nucs between GT..AG)", 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE spliced Viterbi: score comparison and P-state detection";

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

  utest_spliced(r, abcAA, abcDNA, gcode, bgAA, ct, M, N, I);

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
#endif /*p7VITERBI_SP_TESTDRIVE*/
/*--------------------- end, test driver ------------------------*/


/*****************************************************************
 * 5. Benchmark driver.
 *****************************************************************/
#ifdef p7VITERBI_SP_BENCHMARK
/*
   gcc -g -O3 -Wall -std=gnu99 -msse2 -o viterbi_sp_benchmark -I. -I.. -I../.. -I../../easel -L../.. -L../../easel -Dp7VITERBI_SP_BENCHMARK viterbi_sp.c -lhmmer -leasel -lm
   ./viterbi_sp_benchmark <hmmfile>
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
  { "-T",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "also benchmark SplicedTrace after each DP",      0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char benchmark_usage[]  = "[-options] <hmmfile>";
static char benchmark_banner[] = "benchmark driver for SSE spliced Viterbi DP";

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
  P7_FS_OPROFILE *om_fs        = NULL;
  ESL_GENCODE    *gcode        = NULL;
  P7_CODONTABLE  *codon_table  = NULL;
  ESL_SQ         *sq           = NULL;
  SPLICE_PIPELINE *pli         = NULL;
  P7_OMX         *ox           = NULL;
  OSPLICE_SCORES *oss          = NULL;
  P7_TRACE       *tr           = NULL;
  int             N            = esl_opt_GetInteger(go, "-N");
  int             I            = esl_opt_GetInteger(go, "-I");
  int             intron_total = I + 4;
  int             do_T         = esl_opt_GetBoolean(go, "-T");
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
  om_fs       = p7_fs_oprofile_Create(hmm->M, abcAA, 1);
  codon_table = p7_codontable_Create(gcode);
  sq          = esl_sq_CreateDigital(abcAA);

  p7_ProfileConfig   (hmm, bgAA,        gm,    hmm->M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, hmm->M, p7_UNILOCAL);
  p7_fs_oprofile_Convert(gm_tr, om_fs);

  pli = p7_splicepipeline_Create(NULL, hmm->M, hmm->M * 3);
  ox  = p7_omx_Create_dpf(hmm->M, hmm->M * 3, hmm->M * 3, p7X_NSCELLS_SP);
  oss = p7_osplicescores_Create(hmm->M);
  tr  = p7_trace_fs_Create();

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
      p7_fs_oprofile_ReconfigLength(om_fs, L_dna_total / 3);
      p7_omx_GrowTo_dpf(ox, hmm->M, L_dna_total, L_dna_total);
      p7_osplicescores_GrowTo(oss, hmm->M);

      p7_Viterbi_SplicedGlobal(oss, dsq, om_fs, ox, 1, L_dna_total, 1, hmm->M, pli->min_intron);
      final_C = ox->xmx[L_dna_total * p7X_NXCELLS + p7X_C];
      if (final_C != 0.0f && do_T) {
        p7_Viterbi_SplicedTrace(dsq, om_fs, ox, oss->signal_scores,
                                tr, 1, L_dna_total, 1, hmm->M, pli->min_intron);
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
  p7_osplicescores_Destroy(oss);
  p7_omx_Destroy(ox);
  p7_splicepipeline_Destroy(pli);
  esl_sq_Destroy(sq);
  p7_codontable_Destroy(codon_table);
  p7_fs_oprofile_Destroy(om_fs);
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
#endif /*p7VITERBI_SP_BENCHMARK*/
/*---------------- end, benchmark driver ----------------*/

