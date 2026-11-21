/* Checkpointed spliced Viterbi (SSE).
 *
 * Contents:
 *    1. p7_Viterbi_SplicedChk()
 *    2. p7_Viterbi_SplicedTraceChk()
 *    3. Unit tests.
 *    4. Test driver.
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
#include "p7_splice.h"
#include "impl_sse.h"

/* row_get(): one float DP cell from row <i> of a striped matrix given by
 * its row pointers <dp>: state <s> of model position <k>. */
static inline float
row_get(__m128 **dp, int s, int i, int k, int Q)
{
  union { __m128 v; float p[4]; } u;
  u.v = dp[i][p7X_NSCELLS * ((k-1) % Q) + s];
  return u.p[(k-1)/Q];
}

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


/*****************************************************************
 * 1. p7_Viterbi_SplicedChk()
 *
 * The fill computes the same rows as p7_Viterbi_Spliced(), but into
 * a P7_CHKOMX that holds one block of rows at a time. At the start of
 * row b-3 for each block start b, the donor scores (don_ov) and the rows
 * they came from (don_idx, the working int vectors of the P7_CHKOIVX)
 * are saved as checkpoint b/blockR - 1; at row b the P7_CHKOMX saves the
 * saveR rows before b. To compute block b again, the trace restores
 * both, recomputes the acceptor P values and donor updates of rows
 * b-3..b-1 (their M, I, D rows are among the saved rows), and then
 * computes the rows of the block in full. Both use
 * viterbi_spliced_row(), so the rows come out identical.
 *
 * saveR must be at least min_intron+6: row i reads row i-3, and its
 * donor update reads row i-min_intron-3.
 *****************************************************************/

/* viterbi_spliced_init()
 *
 * Clear the donor scores, donor rows and acceptor P values and set rows
 * 0..2, as at the start of p7_Viterbi_Spliced(). With <do_xmx>, also
 * set the special states of every row. Inlined, because both the fill
 * and the trace (when it computes block 0 again) start this way. Block 0
 * must be loaded, and <cx->M>, <cx->L> set.
 */
static inline void
viterbi_spliced_init(const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                     int global_start, int do_xmx)
{
  __m128  **acc_ovx = acc_ov->ivx;
  __m128  **don_ovx = don_ov->ivx;
  __m128i **don_idx = cv->ivx;
  __m128    infv    = _mm_set1_ps(-eslINFINITY);
  __m128i   nonev   = _mm_set1_epi32(-1);
  float     xB      = om_tr->xf[p7O_N][p7O_MOVE];
  int       Q       = p7O_NQF(om_tr->M);
  int       L       = cx->L;
  int       ri, q, j;

  for (ri = 0; ri < SPLICE_ROWS; ri++)
    for (q = 0; q < Q; q++)
      acc_ovx[ri][q] = infv;

  for (j = 0; j < SIGNAL_MEM_SIZE; j++)
    for (q = 0; q < Q; q++) {
      don_ovx[j][q] = infv;
      don_idx[j][q] = nonev;
    }

  for (ri = 0; ri <= 2 && ri <= L; ri++)
    for (q = 0; q < Q; q++)
      MMO(cx->dpf[ri], q) = DMO(cx->dpf[ri], q) = IMO(cx->dpf[ri], q) = infv;

  if (do_xmx)
    {
      for (ri = 0; ri <= L; ri++) {
        cx->xmx[ri * p7X_NXCELLS + p7X_SCALE] = 0.0f;
        cx->xmx[ri * p7X_NXCELLS + p7X_E]     = -eslINFINITY;
        cx->xmx[ri * p7X_NXCELLS + p7X_N]     = -eslINFINITY;
        cx->xmx[ri * p7X_NXCELLS + p7X_J]     = -eslINFINITY;
        cx->xmx[ri * p7X_NXCELLS + p7X_B]     = -eslINFINITY;
        cx->xmx[ri * p7X_NXCELLS + p7X_C]     = -eslINFINITY;
      }
      cx->xmx[0 * p7X_NXCELLS + p7X_B] = xB;
      if (!global_start) cx->xmx[0 * p7X_NXCELLS + p7X_N] = 0.f;

      if (L >= 1) { cx->xmx[1 * p7X_NXCELLS + p7X_N] = 0.f; cx->xmx[1 * p7X_NXCELLS + p7X_B] = xB; }
      if (L >= 2) { cx->xmx[2 * p7X_NXCELLS + p7X_N] = 0.f; cx->xmx[2 * p7X_NXCELLS + p7X_B] = xB; }
    }
}

/* viterbi_spliced_row()
 *
 * Compute row <i> (i >= 3) of the checkpointed spliced Viterbi matrix
 * <cx>, as p7_Viterbi_Spliced() computes it: the acceptor P values
 * of row i, then (if <do_mdi>) its M, I, D and special states, then its
 * donor updates, which fold the look-back row i-min_intron-3 into the
 * donor scores and record in don_idx the row each score came from (of
 * equal scores, the later row). With <do_mdi> FALSE, only the P values
 * and donor updates are done, for a row whose M, I, D values are already
 * loaded. Inlined, because the fill and the trace both compute rows.
 */
static inline void
viterbi_spliced_row(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                    const float *signal_scores, int i_start, int min_intron, int global_start, int global_end, int i, int do_mdi)
{
  register __m128 mpv, dpv, ipv;
  register __m128 sv, msv, dcv;
  __m128   infv      = _mm_set1_ps(-eslINFINITY);
  __m128   tsc_p_vec = _mm_set1_ps(TSC_P);
  __m128   sig_gtag  = _mm_set1_ps(signal_scores[p7S_GTAG]);
  __m128   sig_gcag  = _mm_set1_ps(signal_scores[p7S_GCAG]);
  __m128   sig_atac  = _mm_set1_ps(signal_scores[p7S_ATAC]);
  __m128   b_entry;
  __m128   ppv;
  __m128   lb_mpv, lb_dpv, lb_tsv;
  __m128   new_sv, mask, lbv;
  __m128  *dpc, *dpp3, *dpp_lb;
  __m128  *tp;
  __m128 **acc_ovx   = acc_ov->ivx;
  __m128 **don_ovx   = don_ov->ivx;
  __m128i **don_idx  = cv->ivx;
  float    xB        = om_tr->xf[p7O_N][p7O_MOVE];
  float    xE;
  int      Q         = p7O_NQF(om_tr->M);
  int      pv_i      =  i    % SPLICE_ROWS;
  int      pv_pi     = (i-3) % SPLICE_ROWS;
  int      lb        = i - min_intron - 3;
  int      acc0, acc1, acc2;
  int      don0, don1, don2;
  int      r, s, v, w, x;
  int      C0;
  int      C1[p7P_MAXNUC + 1];
  int      C2[p7P_MAXNUC + 1];
  int      nuc1, nuc3;
  int      slot;
  int      q, j;

  v = nuc_at(sub_dsq, i_start, i-2);
  w = nuc_at(sub_dsq, i_start, i-1);
  x = nuc_at(sub_dsq, i_start, i);

  C0 = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);

  /* Acceptor P values of row i, from the donor scores as they are at
   * the start of row i. */
  acc2 = acceptor_at(sub_dsq, i_start, i);
  acc1 = acceptor_at(sub_dsq, i_start, i-1);
  acc0 = acceptor_at(sub_dsq, i_start, i-2);

  for (nuc1 = 0; nuc1 < p7P_MAXNUC; nuc1++)
    C1[nuc1] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, w, x), p7P_DEGEN1_C);
  C1[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(p7P_MAXCODONS1, w, x), p7P_DEGEN1_C);
  nuc3 = ESL_MIN(x, p7P_MAXNUC);

  for (q = 0; q < Q; q++) {
    __m128 psv = infv;
    if (acc0 == ACCEPT_AG) {
      psv = _mm_max_ps(psv, _mm_add_ps(
          _mm_max_ps(_mm_add_ps(don_ovx[p7S_GTAG][q], sig_gtag),
                        _mm_add_ps(don_ovx[p7S_GCAG][q], sig_gcag)),
          om_tr->rfv[C0][q]));
    } else if (acc0 == ACCEPT_AC) {
      psv = _mm_max_ps(psv, _mm_add_ps(
          _mm_add_ps(don_ovx[p7S_ATAC][q], sig_atac),
          om_tr->rfv[C0][q]));
    }
    if (acc1 == ACCEPT_AG) {
      for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
        psv = _mm_max_ps(psv, _mm_add_ps(
            _mm_max_ps(_mm_add_ps(don_ovx[SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + p7S_GTAG][q], sig_gtag),
                          _mm_add_ps(don_ovx[SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + p7S_GCAG][q], sig_gcag)),
            om_tr->rfv[C1[nuc1]][q]));
      }
    } else if (acc1 == ACCEPT_AC) {
      for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
        psv = _mm_max_ps(psv, _mm_add_ps(
            _mm_add_ps(don_ovx[SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + p7S_ATAC][q], sig_atac),
            om_tr->rfv[C1[nuc1]][q]));
      }
    }
    if (acc2 == ACCEPT_AG) {
      psv = _mm_max_ps(psv,
          _mm_max_ps(_mm_add_ps(don_ovx[SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + p7S_GTAG][q], sig_gtag),
                        _mm_add_ps(don_ovx[SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + p7S_GCAG][q], sig_gcag)));
    } else if (acc2 == ACCEPT_AC) {
      psv = _mm_max_ps(psv,
          _mm_add_ps(don_ovx[SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + p7S_ATAC][q], sig_atac));
    }
    acc_ovx[pv_i][q] = psv;
  }

  /* M, I, D and special states of row i */
  if (do_mdi)
    {
      dpc  = cx->dpf[i];
      dpp3 = cx->dpf[i - 3];

      if (!global_start) {
        float xN_im3 = cx->xmx[(i-3) * p7X_NXCELLS + p7X_N];
        float xN_i   = xN_im3 + om_tr->xf[p7O_N][p7O_LOOP];
        float xB_i   = xN_i   + om_tr->xf[p7O_N][p7O_MOVE];
        cx->xmx[i * p7X_NXCELLS + p7X_N] = xN_i;
        cx->xmx[i * p7X_NXCELLS + p7X_B] = xB_i;
      }

      mpv = esl_sse_rightshift_ps(MMO(dpp3, Q - 1), infv);
      dpv = esl_sse_rightshift_ps(DMO(dpp3, Q - 1), infv);
      ipv = esl_sse_rightshift_ps(IMO(dpp3, Q - 1), infv);

      ppv = esl_sse_rightshift_ps(acc_ovx[pv_pi][Q - 1], infv);

      tp  = om_tr->tfv;
      dcv = infv;

      if (global_start) {
        if (i == 3) {
          union { __m128 v; float p[4]; } u2;
          u2.v    = infv;
          u2.p[0] = xB;
          b_entry = u2.v;
        } else {
          b_entry = infv;
        }
      } else {
        b_entry = _mm_set1_ps(cx->xmx[(i-3) * p7X_NXCELLS + p7X_B]);
      }

      for (q = 0; q < Q; q++) {
        __m128 ppv_next = acc_ovx[pv_pi][q];

        tp++;                                                    /* skip BM */
        sv  =                  _mm_add_ps(mpv, *tp); tp++;   /* MM */
        sv  = _mm_max_ps(sv, _mm_add_ps(ipv, *tp)); tp++; /* IM */
        sv  = _mm_max_ps(sv, _mm_add_ps(dpv, *tp)); tp++; /* DM */
        sv  = _mm_max_ps(sv, _mm_add_ps(ppv, tsc_p_vec)); /* PM */
        sv  = _mm_max_ps(sv, b_entry);
        if (global_start) b_entry = infv;

        mpv = MMO(dpp3, q);
        dpv = DMO(dpp3, q);
        ipv = IMO(dpp3, q);

        msv = _mm_add_ps(sv, om_tr->rfv[C0][q]);
        MMO(dpc, q) = msv;
        DMO(dpc, q) = dcv;
        dcv = _mm_add_ps(msv, *tp); tp++;                    /* MD */

        sv  = _mm_add_ps(mpv, *tp); tp++;                    /* MI */
        sv  = _mm_max_ps(sv, _mm_add_ps(ipv, *tp)); tp++; /* II */
        { /* stop-codon mask */
          __m128 stop_mask = _mm_cmpeq_ps(om_tr->rfv[C0][q], infv);
          sv = _mm_or_ps(_mm_andnot_ps(stop_mask, sv), _mm_and_ps(stop_mask, infv));
        }
        IMO(dpc, q) = sv;

        ppv = ppv_next;
      }

      /* DD sweep */
      dcv = esl_sse_rightshift_ps(dcv, infv);
      DMO(dpc, 0) = infv;
      tp = om_tr->tfv + 7 * Q;
      for (q = 0; q < Q; q++) {
        DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
        dcv = _mm_add_ps(DMO(dpc, q), *tp); tp++;
      }
      if (om_tr->M < 100) {
        for (j = 1; j < 4; j++) {
          dcv = esl_sse_rightshift_ps(dcv, infv);
          tp  = om_tr->tfv + 7 * Q;
          for (q = 0; q < Q; q++) {
            DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
            dcv = _mm_add_ps(dcv, *tp); tp++;
          }
        }
      } else {
        for (j = 1; j < 4; j++) {
          register __m128 cv_any;
          dcv    = esl_sse_rightshift_ps(dcv, infv);
          tp     = om_tr->tfv + 7 * Q;
          cv_any = _mm_setzero_ps();
          for (q = 0; q < Q; q++) {
            sv             = _mm_max_ps(dcv, DMO(dpc, q));
            cv_any         = _mm_or_ps(cv_any, _mm_cmpgt_ps(sv, DMO(dpc, q)));
            DMO(dpc, q)    = sv;
            dcv            = _mm_add_ps(dcv, *tp); tp++;
          }
          if (!_mm_movemask_ps(cv_any)) break;
        }
      }

      /* Semi-global exit */
      if (!global_end) {
        __m128 xEv = infv;
        for (q = 0; q < Q; q++)
          xEv = _mm_max_ps(xEv, _mm_max_ps(MMO(dpc, q), DMO(dpc, q)));
        xEv = _mm_max_ps(xEv, _mm_movehl_ps(xEv, xEv));
        xEv = _mm_max_ps(xEv, _mm_shuffle_ps(xEv, xEv, 0x55));
        { union { __m128 v; float p[4]; } uE; uE.v = xEv; xE = uE.p[0]; }
        cx->xmx[i * p7X_NXCELLS + p7X_E] = xE;
        cx->xmx[i * p7X_NXCELLS + p7X_C] = ESL_MAX(
            cx->xmx[(i-3) * p7X_NXCELLS + p7X_C] + om_tr->xf[p7O_C][p7O_LOOP],
            xE + om_tr->xf[p7O_E][p7O_MOVE]);
      }
    }

  /* Donor updates of row i. Where the new score is at least the kept
   * one, it replaces it and row lb replaces the kept row. */
  if (i < min_intron + 3) return;

  don2 = donor_at(sub_dsq, i_start, min_intron, i);
  don1 = donor_at(sub_dsq, i_start, min_intron, i-1);
  don0 = donor_at(sub_dsq, i_start, min_intron, i-2);
  if (don0 < 0 && don1 < 0 && don2 < 0) return;

  r      = nuc_at(sub_dsq, i_start, i - min_intron - 2);
  s      = nuc_at(sub_dsq, i_start, i - min_intron - 1);
  dpp_lb = cx->dpf[lb];
  lbv    = _mm_castsi128_ps(_mm_set1_epi32(lb));

  if (don2 >= 0) {
    for (nuc3 = 0; nuc3 < p7P_MAXNUC; nuc3++)
      C2[nuc3] = p7P_MINIDX(p7P_CODON3_FS1(r, s, nuc3), p7P_DEGEN1_C);
    C2[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(r, s, p7P_MAXCODONS1), p7P_DEGEN1_C);

    lb_mpv = esl_sse_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
    lb_dpv = esl_sse_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
    for (q = 0; q < Q; q++) {
      lb_tsv = _mm_max_ps(lb_mpv, lb_dpv);
      for (nuc3 = 0; nuc3 <= p7P_MAXNUC; nuc3++) {
        slot   = SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + don2;
        new_sv = _mm_add_ps(lb_tsv, om_tr->rfv[C2[nuc3]][q]);
        mask   = _mm_cmpge_ps(new_sv, don_ovx[slot][q]);
        don_ovx[slot][q] = _mm_or_ps(_mm_and_ps(mask, new_sv), _mm_andnot_ps(mask, don_ovx[slot][q]));
        don_idx[slot][q] = _mm_castps_si128(_mm_or_ps(_mm_and_ps(mask, lbv), _mm_andnot_ps(mask, _mm_castsi128_ps(don_idx[slot][q]))));
      }
      lb_mpv = MMO(dpp_lb, q);
      lb_dpv = DMO(dpp_lb, q);
    }
  }

  if (don1 >= 0) {
    nuc1   = ESL_MIN(r, p7P_MAXNUC);
    slot   = SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + don1;
    lb_mpv = esl_sse_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
    lb_dpv = esl_sse_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
    for (q = 0; q < Q; q++) {
      new_sv = _mm_max_ps(lb_mpv, lb_dpv);
      mask   = _mm_cmpge_ps(new_sv, don_ovx[slot][q]);
      don_ovx[slot][q] = _mm_or_ps(_mm_and_ps(mask, new_sv), _mm_andnot_ps(mask, don_ovx[slot][q]));
      don_idx[slot][q] = _mm_castps_si128(_mm_or_ps(_mm_and_ps(mask, lbv), _mm_andnot_ps(mask, _mm_castsi128_ps(don_idx[slot][q]))));
      lb_mpv = MMO(dpp_lb, q);
      lb_dpv = DMO(dpp_lb, q);
    }
  }

  if (don0 >= 0) {
    slot   = don0;
    lb_mpv = esl_sse_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
    lb_dpv = esl_sse_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
    for (q = 0; q < Q; q++) {
      new_sv = _mm_max_ps(lb_mpv, lb_dpv);
      mask   = _mm_cmpge_ps(new_sv, don_ovx[slot][q]);
      don_ovx[slot][q] = _mm_or_ps(_mm_and_ps(mask, new_sv), _mm_andnot_ps(mask, don_ovx[slot][q]));
      don_idx[slot][q] = _mm_castps_si128(_mm_or_ps(_mm_and_ps(mask, lbv), _mm_andnot_ps(mask, _mm_castsi128_ps(don_idx[slot][q]))));
      lb_mpv = MMO(dpp_lb, q);
      lb_dpv = DMO(dpp_lb, q);
    }
  }
}


/* Function:  p7_Viterbi_SplicedChk()
 *
 * Purpose:   Checkpointed spliced Viterbi fill: same as
 *            p7_Viterbi_Spliced(), but the rows go into the
 *            checkpointed matrix <cx>, with checkpoints of the donor
 *            scores and their rows in <cv>, so that
 *            p7_Viterbi_SplicedTraceChk() can compute any block of
 *            rows again. On return the last block is loaded, and the
 *            special states of every row are in <cx->xmx>.
 *
 *            <cx> must be grown to the model length of <om_tr> by
 *            <i_end>-<i_start>+1 and save at least <min_intron>+6 rows
 *            per checkpoint. <cv> must have SIGNAL_MEM_SIZE float and int
 *            channels, room for the model length, and (L / <cx->blockR>)
 *            checkpoints.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <cx>, <cv> or the vectors don't fit.
 */
int
p7_Viterbi_SplicedChk(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                          const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                          int i_start, int i_end, int min_intron, int global_start, int global_end)
{
  int   M      = om_tr->M;
  int   L      = i_end - i_start + 1;
  int   Q      = p7O_NQF(M);
  int   blockR = cx->blockR;
  int   b, nb, end, i;
  int   status;

  if (om_tr->codon_lengths != 1)                                          ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");
  if (cx->nscells != p7X_NSCELLS)                                         ESL_EXCEPTION(eslEINVAL, "DP matrix must have p7X_NSCELLS cells per stripe position");
  if (cx->saveR < min_intron + 6)                                         ESL_EXCEPTION(eslEINVAL, "checkpoints save %d rows; need %d", cx->saveR, min_intron + 6);
  if (cx->allocQ4 < Q || cx->allocP < L + 1 || cx->allocXR < L + 1)   ESL_EXCEPTION(eslEINVAL, "checkpointed matrix too small");
  if (cv->allocCF != SIGNAL_MEM_SIZE || cv->allocCI != SIGNAL_MEM_SIZE)   ESL_EXCEPTION(eslEINVAL, "checkpoint vectors need %d channels", SIGNAL_MEM_SIZE);
  if (cv->allocQ4 < Q || cv->allocK < L / blockR)                     ESL_EXCEPTION(eslEINVAL, "checkpoint vectors too small");
  if (don_ov->allocC < SIGNAL_MEM_SIZE || don_ov->allocQ4 < Q)        ESL_EXCEPTION(eslEINVAL, "donor vectors too small");
  if (acc_ov->allocC < SPLICE_ROWS     || acc_ov->allocQ4 < Q)        ESL_EXCEPTION(eslEINVAL, "acceptor vectors too small");

  cx->M = M;
  cx->L = L;
  if ((status = p7_chkomx_StartBlock(cx, 0)) != eslOK) return status;
  viterbi_spliced_init(om_tr, cx, cv, acc_ov, don_ov, global_start, TRUE);

  for (b = 0; b <= L; b += blockR)
    {
      if (b > 0 && (status = p7_chkomx_StartBlock(cx, b)) != eslOK) return status;
      end = ESL_MIN(L, b + blockR - 1);
      nb  = b + blockR;

      for (i = ESL_MAX(3, b); i <= end; i++)
        {
          /* save the donor state for block nb at the start of row nb-3 */
          if (i == nb - 3 && nb <= L && (status = p7_chkoivx_Save(cv, don_ov, M, nb / blockR - 1)) != eslOK) return status;
          viterbi_spliced_row(sub_dsq, om_tr, cx, cv, acc_ov, don_ov, signal_scores, i_start, min_intron, global_start, global_end, i, TRUE);
        }
    }

  if (global_end) {
    union { __m128 v; float p[4]; } um, ud;
    int   qM = (M - 1) % Q;
    int   rM = (M - 1) / Q;
    float xE;
    um.v = MMO(cx->dpf[L], qM);
    ud.v = DMO(cx->dpf[L], qM);
    xE   = ESL_MAX(um.p[rM], ud.p[rM]);
    cx->xmx[L * p7X_NXCELLS + p7X_E] = xE;
    cx->xmx[L * p7X_NXCELLS + p7X_C] = xE + om_tr->xf[p7O_E][p7O_MOVE];
  }

  return eslOK;
}


/*****************************************************************
 * 2. p7_Viterbi_SplicedTraceChk()
 *****************************************************************/

/* The donor behind a P state score, as chosen by replay_donors_from() */
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

/* Donor slot scores for one model column, laid out as the don_ov
 * channels, with the row each score came from. */
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
 * Starting from the slot scores and rows in <ds> as they were at the
 * start of row <s0>, this repeats the fill's donor updates for column
 * <kk> alone up to <acc_i>, reading each look-back row through <dp>,
 * and records each slot's row along with its score. Of equal scores the
 * later (nearer) donor is kept. The acceptor is then scored as in the
 * fill, again keeping the nearer donor of equal scores. <ds> is
 * modified.
 *
 * Returns <eslOK> with the P state score in <ret_sc>, the donor row
 * (the last row before the split codon) in <ret_donor_i>, the number of
 * nucleotides of the split codon before the donor in <ret_c>, and the
 * splice signal in <ret_don_sig>. Returns <eslFAIL> if no donor reaches
 * this acceptor.
 */
static int
replay_donors_from(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, __m128 **dp, int Q, const float *signal_scores,
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

  if (kk < 2 || acc_i < min_intron + 3) return eslFAIL;

  for (iu = s0; iu < acc_i; iu++) {
    don2 = donor_at(sub_dsq, i_start, min_intron, iu);
    don1 = donor_at(sub_dsq, i_start, min_intron, iu - 1);
    don0 = donor_at(sub_dsq, i_start, min_intron, iu - 2);
    if (don0 < 0 && don1 < 0 && don2 < 0) continue;

    r   = nuc_at(sub_dsq, i_start, iu - min_intron - 2);
    lb  = iu - min_intron - 3;
    tsc = ESL_MAX(row_get(dp, p7X_M, lb, kk-1, Q), row_get(dp, p7X_D, lb, kk-1, Q));

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


/* Function:  p7_Viterbi_SplicedTraceChk()
 *
 * Purpose:   Traceback of a matrix filled by
 *            p7_Viterbi_SplicedChk() with the same arguments, giving
 *            the same trace as p7_Viterbi_SplicedTrace() gives for a
 *            full matrix. The trace runs from the last block back to the
 *            first, computing each block it passes through again from
 *            its checkpoint and skipping blocks that lie inside introns.
 *            Each P state's donor comes from replaying donor updates from
 *            the checkpoint of the block that holds the acceptor. The M
 *            state scores on the trace are recorded as it goes, for
 *            <max_drop>.
 *
 *            <om_tr>, <acc_ov>, <don_ov>, <global_start> and <global_end>
 *            are as passed to p7_Viterbi_SplicedChk(), for computing
 *            blocks again.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if the profile is not 1-codon-length or <cx>
 *                        isn't filled for this comparison.
 *            <eslFAIL>   if traceback fails.
 */
int
p7_Viterbi_SplicedTraceChk(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                               const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                               P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                               int min_intron, int global_start, int global_end, float *vitsc, float *max_drop)
{
  float const *tsc    = gm_tr->tsc;
  float        r_tol  = 1e-5;
  float        a_tol  = 1e-4;
  int          M      = k_end - k_start + 1;
  int          L      = i_end - i_start + 1;
  int          Q      = p7O_NQF(om_tr->M);
  int          blockR = cx->blockR;
  int          i      = L;
  int          k      = 0;
  int          b, end, i2, s0, nuc, sig, slot, q, rl;
  int          z;
  int          v, w, x;
  int          c, c3;
  int          sprv, scur;
  int          sub_k;
  int          donor_i;
  int          don_sig;
  int          P_cnt, P_idx;
  float        P_state;
  float        emit;
  float        vsc = 0.0f;
  float        curr_P_score;
  float        curr_M_score;
  float        max_M_score;
  float        curr_max_drop;
  float       *P_scores;
  float       *M_scores;     /* M state score of each trace position, in trace order */
  int          M_alloc;
  float        slot_sc[SIGNAL_MEM_SIZE];
  int          slot_row[SIGNAL_MEM_SIZE];
  DONOR_SLOTS  ds;
  __m128      *ckv;
  void        *p;
  int          status;

  P_scores = NULL;
  M_scores = NULL;
  P_cnt    = 0;
  M_alloc  = 256;
  ESL_ALLOC(P_scores, sizeof(float) * M);
  ESL_ALLOC(M_scores, sizeof(float) * M_alloc);

#define OMMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : p7_chkomx_FGetMDI(cx, p7X_M, (ii), (kk)))
#define ODMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : p7_chkomx_FGetMDI(cx, p7X_D, (ii), (kk)))
#define OIMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : p7_chkomx_FGetMDI(cx, p7X_I, (ii), (kk)))
#define OXMXo(ii,s)  (cx->xmx[(ii) * p7X_NXCELLS + (s)])

  if (gm_tr->codon_lengths != 1) ESL_XEXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");
  if (cx->M != M || cx->L != L || cx->b != (L / blockR) * blockR) ESL_XEXCEPTION(eslEINVAL, "matrix isn't filled for this comparison");

  vsc = OXMXo(i, p7X_C) + gm_tr->xsc[p7P_C][p7P_MOVE];

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i+i_start-1, 0)) != eslOK) goto ERROR;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i+i_start-1, 0)) != eslOK) goto ERROR;
  M_scores[0] = M_scores[1] = 0.;

  sprv = p7T_C;

  while (sprv != p7T_S) {

    /* E, M, D and I read row i (and i-3), and P reads row donor_i. If
     * that row is before the loaded block, load its block and compute it
     * again: from the start for block 0, otherwise from the block's
     * checkpoint, redoing the P values and donor updates of rows b-3..b-1
     * and then the block's rows in full. */
    if (sprv == p7T_E || sprv == p7T_M || sprv == p7T_D || sprv == p7T_I || sprv == p7T_P)
      {
        b = (((sprv == p7T_P) ? donor_i : i) / blockR) * blockR;
        if (b != cx->b)
          {
            if ((status = p7_chkomx_LoadBlock(cx, b)) != eslOK) goto ERROR;
            end = ESL_MIN(L, b + blockR - 1);
            if (b == 0)
              {
                viterbi_spliced_init(om_tr, cx, cv, acc_ov, don_ov, global_start, FALSE);
                for (i2 = 3; i2 <= end; i2++)
                  viterbi_spliced_row(sub_dsq, om_tr, cx, cv, acc_ov, don_ov, signal_scores, i_start, min_intron, global_start, global_end, i2, TRUE);
              }
            else
              {
                if ((status = p7_chkoivx_Load(cv, don_ov, M, b / blockR - 1)) != eslOK) goto ERROR;
                for (i2 = b - 3; i2 < b; i2++)
                  viterbi_spliced_row(sub_dsq, om_tr, cx, cv, acc_ov, don_ov, signal_scores, i_start, min_intron, global_start, global_end, i2, FALSE);
                for (i2 = b; i2 <= end; i2++)
                  viterbi_spliced_row(sub_dsq, om_tr, cx, cv, acc_ov, don_ov, signal_scores, i_start, min_intron, global_start, global_end, i2, TRUE);
              }
            if (global_end && end == L) {
              union { __m128 v; float p[4]; } um, ud;
              um.v = MMO(cx->dpf[L], (M - 1) % Q);
              ud.v = DMO(cx->dpf[L], (M - 1) % Q);
              OXMXo(L, p7X_E) = ESL_MAX(um.p[(M - 1) / Q], ud.p[(M - 1) / Q]);
              OXMXo(L, p7X_C) = OXMXo(L, p7X_E) + om_tr->xf[p7O_E][p7O_MOVE];
            }
          }
      }

    switch (sprv) {

    case p7T_C:
      if      (OXMXo(i, p7X_C) < OXMXo(i-2, p7X_C) || OXMXo(i, p7X_C) < OXMXo(i-1, p7X_C))                           scur = p7T_C;
      else if (OXMXo(i, p7X_C) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
      else if (esl_FCompare(OXMXo(i, p7X_C), OXMXo(i-3, p7X_C) + gm_tr->xsc[p7P_C][p7P_LOOP], r_tol, a_tol) == eslOK) scur = p7T_C;
      else if (esl_FCompare(OXMXo(i, p7X_C), OXMXo(i,   p7X_E) + gm_tr->xsc[p7P_E][p7P_MOVE], r_tol, a_tol) == eslOK) scur = p7T_E;
      else ESL_XEXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7T_E:
      if (OXMXo(i, p7X_E) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
      for (k = M; k >= 1; k--) {
        if (esl_FCompare(OXMXo(i, p7X_E), OMMo(i,k), r_tol, a_tol) == eslOK) { scur = p7T_M; break; }
        if (esl_FCompare(OXMXo(i, p7X_E), ODMo(i,k), r_tol, a_tol) == eslOK) { scur = p7T_D; break; }
      }
      if (k == 0) ESL_XEXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      break;

    case p7T_M:
      if (OMMo(i,k) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);

      v     = nuc_at(sub_dsq, i_start, i-2);
      w     = nuc_at(sub_dsq, i_start, i-1);
      x     = nuc_at(sub_dsq, i_start, i);
      c3    = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);
      sub_k = k_start + k - 1;
      emit  = p7P_MSC_CODON(gm_tr, sub_k, c3);

      scur = -1;
      if      (esl_FCompare(OMMo(i,k), OMMo(i-3,k-1) + TSC(p7P_MM,sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(OMMo(i,k), OIMo(i-3,k-1) + TSC(p7P_IM,sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare(OMMo(i,k), ODMo(i-3,k-1) + TSC(p7P_DM,sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_D;
      else if (esl_FCompare(OMMo(i,k), OXMXo(i-3, p7X_B)                   + emit, r_tol, a_tol) == eslOK) scur = p7T_B;
      else {    /* The P state is not stored - recreate it from the donor updates */

        /* Replay from the loaded block's checkpoint (the start of row
         * b-3, which is at or before the acceptor row i-3), or from the
         * first row with donor updates for block 0. In the checkpoint,
         * channel c's value for column k-1 is lane (k-2)/Q of vector
         * (k-2)%Q; float channels come first, then the int channels. */
        if (cx->b == 0)
          {
            donor_slots_clear(&ds);
            s0 = min_intron + 3;
          }
        else
          {
            ckv = (__m128 *) (((unsigned long int) ((char *) cv->ck_mem + 15)) & (~0xf));
            ckv += (size_t) (cx->b / blockR - 1) * (cv->allocCF + cv->allocCI) * cv->allocQ4;
            q  = (k-2) % Q;
            rl = (k-2) / Q;
            for (slot = 0; slot < SIGNAL_MEM_SIZE; slot++) {
              union { __m128 v; float p[4]; int32_t n[4]; } u;
              u.v = ckv[(size_t) slot * cv->allocQ4 + q];                       slot_sc[slot]  = u.p[rl];
              u.v = ckv[(size_t) (cv->allocCF + slot) * cv->allocQ4 + q];      slot_row[slot] = u.n[rl];
            }
            for (sig = 0; sig < p7S_SPLICE_SIGNALS; sig++) {
              ds.sc0[sig]  = slot_sc[sig];
              ds.row0[sig] = slot_row[sig];
              for (nuc = 0; nuc <= p7P_MAXNUC; nuc++) {
                ds.sc1[nuc][sig]  = slot_sc [SPLICE_OFFSET_1 + nuc * p7S_SPLICE_SIGNALS + sig];
                ds.row1[nuc][sig] = slot_row[SPLICE_OFFSET_1 + nuc * p7S_SPLICE_SIGNALS + sig];
                ds.sc2[nuc][sig]  = slot_sc [SPLICE_OFFSET_2 + nuc * p7S_SPLICE_SIGNALS + sig];
                ds.row2[nuc][sig] = slot_row[SPLICE_OFFSET_2 + nuc * p7S_SPLICE_SIGNALS + sig];
              }
            }
            s0 = cx->b - 3;
          }

        if (replay_donors_from(sub_dsq, gm_tr, cx->dpf, Q, signal_scores, i_start, k_start, min_intron, s0, &ds, i-3, k-1, &P_state, &donor_i, &c, &don_sig) != eslOK)
          ESL_XEXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
        if (esl_FCompare(OMMo(i,k), P_state + TSC_P + emit, r_tol, a_tol) != eslOK)
          ESL_XEXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);

        scur = p7T_P;
        vsc -= TSC_P;
        vsc -= signal_scores[don_sig];

        /* Track P state scores for x_drop calculation */
        P_scores[P_cnt++] = TSC_P + signal_scores[don_sig];
      }

      k--; i -= 3;
      break;

    case p7T_D:
      if (ODMo(i,k) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);
      sub_k = k_start + k - 1;
      if      (esl_FCompare(ODMo(i,k), OMMo(i,k-1) + TSC(p7P_MD,sub_k-1), r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(ODMo(i,k), ODMo(i,k-1) + TSC(p7P_DD,sub_k-1), r_tol, a_tol) == eslOK) scur = p7T_D;
      else ESL_XEXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      k--;
      break;

    case p7T_I:
      if (OIMo(i,k) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);
      sub_k = k_start + k - 1;
      if      (esl_FCompare(OIMo(i,k), OMMo(i-3,k) + TSC(p7P_MI,sub_k), r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(OIMo(i,k), OIMo(i-3,k) + TSC(p7P_II,sub_k), r_tol, a_tol) == eslOK) scur = p7T_I;
      else ESL_XEXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      i -= 3;
      break;

    case p7T_P:
      scur = (OMMo(donor_i,k-1) >= ODMo(donor_i,k-1)) ? p7T_M : p7T_D;
      k--; i = donor_i;
      break;

    case p7T_N:
      if (i > 0 && OXMXo(i, p7X_N) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = (i == 0) ? p7T_S : p7T_N;
      break;

    case p7T_B:
      vsc += p7P_TSC(gm_tr, k, p7P_BM);
      if (OXMXo(i, p7X_B) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible B reached at i=%d", i);
      if (i == 0) { scur = p7T_N; break; }
      if (esl_FCompare(OXMXo(i, p7X_B), OXMXo(i, p7X_N) + gm_tr->xsc[p7P_N][p7P_MOVE], r_tol, a_tol) == eslOK) scur = p7T_N;
      else ESL_XEXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    default: ESL_XEXCEPTION(eslFAIL, "bogus state in traceback");
    }

    if      (scur == p7T_M) c = 3;
    else if (scur != p7T_P) c = 0;

    if ((status = p7_trace_fs_Append(tr, scur, k_start+k-1, i_start+i-1, c)) != eslOK) goto ERROR;

    /* Record the M state's score now, while its row is loaded */
    if (tr->N > M_alloc) {
      M_alloc *= 2;
      ESL_RALLOC(M_scores, p, sizeof(float) * M_alloc);
    }
    M_scores[tr->N-1] = (scur == p7T_M ? OMMo(i, k) : 0.);

    if ((scur == p7T_N || scur == p7T_C) && scur == sprv) i--;

    sprv = scur;
  }

  tr->M = M;
  tr->L = L;

  if (vitsc != NULL) *vitsc = vsc;

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

#undef OMMo
#undef ODMo
#undef OIMo
#undef OXMXo

  free(P_scores);
  free(M_scores);
  return status;

 ERROR:
  free(P_scores);
  free(M_scores);
  return status;
}


/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef p7VITERBI_SP_CHK_TESTDRIVE

#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_compare_traces()
 *
 * Require two traces of the same matrix to have the same states,
 * positions, splice codon phases and Viterbi scores.
 */
static void
utest_compare_traces(const P7_TRACE *tr1, float sc1, const P7_TRACE *tr2, float sc2, const char *msg)
{
  int z;

  if (tr1->N != tr2->N) esl_fatal("%s: trace lengths %d != %d", msg, tr1->N, tr2->N);
  for (z = 0; z < tr1->N; z++)
    if (tr1->st[z] != tr2->st[z] || tr1->k[z] != tr2->k[z] ||
        tr1->i[z]  != tr2->i[z]  || tr1->c[z] != tr2->c[z])
      esl_fatal("%s: traces differ at position %d", msg, z);
  if (sc1 != sc2) esl_fatal("%s: trace scores %f != %f", msg, sc1, sc2);
}

/* utest_chk()
 *
 * For sequences with a simulated intron, fill with
 * p7_Viterbi_Spliced() and with p7_Viterbi_SplicedChk() at
 * several block sizes, in all three fill modes, and require the same
 * special states on every row. If the final C state is reachable, trace
 * both with p7_Viterbi_SplicedTrace() and
 * p7_Viterbi_SplicedTraceChk() and require the same trace, Viterbi
 * score and max_drop.
 */
static void
utest_chk(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_CODONTABLE *codon_table,
              int M, int N, int intron_len)
{
  char            msg[]        = "utest_chk failed";
  P7_HMM         *hmm          = NULL;
  P7_FS_PROFILE  *gm_tr        = p7_profile_fs_Create(M, abcAA, 1);
  P7_FS_OPROFILE *om_tr        = NULL;
  ESL_SQ         *sq           = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq          = NULL;
  P7_OMX         *ox           = NULL;
  P7_OIVX        *acc_ov       = NULL;
  P7_OIVX        *don_ov       = NULL;
  P7_TRACE       *tr_full      = p7_trace_fs_Create();
  P7_TRACE       *tr_chk       = p7_trace_fs_Create();
  float           signal_scores[p7S_SPLICE_SIGNALS];
  int             modes[3][2]  = { { TRUE, TRUE }, { TRUE, FALSE }, { FALSE, TRUE } };
  int             min_intron   = 13;
  int             blocks[3]    = { 13 + 6, 64, SPLICE_CHECKPOINT_ROWS };
  int             intron_total = intron_len + 4;
  int             L_amino, L_dna_total;
  int             k_start, k_end, sub_M;
  int             i, j, md, bi, s;
  int             n = N;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, M, p7_UNILOCAL);
  om_tr  = p7_fs_oprofile_Create(M, abcAA, 1);
  p7_fs_oprofile_Convert_Log(gm_tr, om_tr);
  ox     = p7_omx_Create_dpf(M, M * 3, M * 3, p7X_NSCELLS);
  acc_ov = p7_oivx_Create(M, SPLICE_ROWS);
  don_ov = p7_oivx_Create(M, SIGNAL_MEM_SIZE);

  p7_SignalScores(signal_scores);
  p7_emit_SimpleConsensus(hmm, sq);
  L_amino = sq->n;

  while (n--)
    {
      k_start = k_end = 0;
      while (k_end <= k_start || k_end > L_amino) {
        k_start = esl_rnd_Roll(r, L_amino / 4) + 1;
        k_end   = esl_rnd_Roll(r, L_amino / 2) + k_start + L_amino / 4;
      }
      sub_M       = k_end - k_start + 1;
      L_dna_total = sub_M * 3 + intron_total;

      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) esl_fatal("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;

      j = 1;
      for (i = k_start; i <= k_start + sub_M / 2; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }
      dsq[j++] = 2; /* G */
      dsq[j++] = 3; /* T */
      for (i = 0; i < intron_len; i++) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = 0; /* A */
      dsq[j++] = 2; /* G */
      for (i = k_start + sub_M / 2 + 1; i <= k_end; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }

      p7_fs_oprofile_SubConvert_Log(gm_tr, om_tr, k_start, k_end);
      p7_fs_oprofile_ReconfigLength_Log(om_tr, L_dna_total / 3);
      p7_fs_ReconfigLength(gm_tr, L_dna_total / 3);
      p7_omx_GrowTo_dpf(ox, sub_M, L_dna_total, L_dna_total);
      p7_oivx_GrowTo(acc_ov, sub_M, SPLICE_ROWS);
      p7_oivx_GrowTo(don_ov, sub_M, SIGNAL_MEM_SIZE);

      for (md = 0; md < 3; md++)
        {
          p7_Viterbi_Spliced(dsq, om_tr, ox, signal_scores, acc_ov, don_ov, 1, L_dna_total, min_intron, modes[md][0], modes[md][1]);

          for (bi = 0; bi < 3; bi++)
            {
              P7_CHKOMX  *cx  = p7_chkomx_Create(sub_M, L_dna_total, blocks[bi], min_intron + 6, p7X_NSCELLS);
              P7_CHKOIVX *cv  = p7_chkoivx_Create(sub_M, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, L_dna_total / blocks[bi]);
              P7_OIVX    *acc = p7_oivx_Create(sub_M, SPLICE_ROWS);
              P7_OIVX    *don = p7_oivx_Create(sub_M, SIGNAL_MEM_SIZE);
              float       sc_full, sc_chk, drop_full, drop_chk;

              if (cx == NULL || cv == NULL || acc == NULL || don == NULL) esl_fatal("%s: allocation failed", msg);
              if (p7_Viterbi_SplicedChk(dsq, om_tr, cx, cv, signal_scores, acc, don, 1, L_dna_total, min_intron, modes[md][0], modes[md][1]) != eslOK)
                esl_fatal("%s: checkpointed fill failed", msg);

              for (i = 0; i <= L_dna_total; i++)
                for (s = 0; s < p7X_NXCELLS; s++)
                  if (cx->xmx[i*p7X_NXCELLS+s] != ox->xmx[i*p7X_NXCELLS+s])
                    esl_fatal("%s: special state %d differs at row %d (block size %d, mode %d)", msg, s, i, blocks[bi], md);

              if (ox->xmx[L_dna_total*p7X_NXCELLS+p7X_C] != -eslINFINITY)
                {
                  if (p7_Viterbi_SplicedTrace(dsq, ox, gm_tr, signal_scores, tr_full, 1, L_dna_total, k_start, k_end, min_intron, &sc_full, &drop_full) != eslOK)
                    esl_fatal("%s: replay trace failed", msg);
                  if (p7_Viterbi_SplicedTraceChk(dsq, om_tr, cx, cv, gm_tr, signal_scores, acc, don, tr_chk, 1, L_dna_total, k_start, k_end,
                                                     min_intron, modes[md][0], modes[md][1], &sc_chk, &drop_chk) != eslOK)
                    esl_fatal("%s: checkpointed trace failed (block size %d, mode %d)", msg, blocks[bi], md);
                  utest_compare_traces(tr_full, sc_full, tr_chk, sc_chk, msg);
                  if (drop_full != drop_chk) esl_fatal("%s: max_drop %f != %f (block size %d, mode %d)", msg, drop_chk, drop_full, blocks[bi], md);
                  p7_trace_Reuse(tr_full);
                  p7_trace_Reuse(tr_chk);
                }

              p7_chkomx_Destroy(cx);
              p7_chkoivx_Destroy(cv);
              p7_oivx_Destroy(acc);
              p7_oivx_Destroy(don);
            }
        }
    }

 cleanup:
  if (dsq)    free(dsq);
  if (ox)     p7_omx_Destroy(ox);
  if (acc_ov) p7_oivx_Destroy(acc_ov);
  if (don_ov) p7_oivx_Destroy(don_ov);
  if (om_tr)  p7_fs_oprofile_Destroy(om_tr);
  if (hmm)    p7_hmm_Destroy(hmm);
  p7_trace_fs_Destroy(tr_full);
  p7_trace_fs_Destroy(tr_chk);
  esl_sq_Destroy(sq);
  p7_profile_fs_Destroy(gm_tr);
}

#endif /*p7VITERBI_SP_CHK_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/


/*****************************************************************
 * 4. Test driver.
 *****************************************************************/
#ifdef p7VITERBI_SP_CHK_TESTDRIVE
/*
  make viterbi_sp_chk_utest
  ./viterbi_sp_chk_utest
*/
#include "esl_alphabet.h"
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name    type           default env range toggles reqs incomp help                                       docgroup */
  { "-h",  eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",  eslARG_INT,     "42", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M",  eslARG_INT,    "200", NULL,"n>=5",NULL,NULL, NULL, "size of random models to sample",                0 },
  { "-N",  eslARG_INT,     "20", NULL, NULL, NULL, NULL, NULL, "number of random sequences to sample",           0 },
  { "-I",  eslARG_INT,    "500", NULL,"n>0", NULL, NULL, NULL, "simulated intron length",                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the checkpointed spliced Viterbi";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA  = esl_alphabet_Create(eslAMINO);
  ESL_ALPHABET   *abcDNA = esl_alphabet_Create(eslDNA);
  ESL_GENCODE    *gcode  = esl_gencode_Create(abcDNA, abcAA);
  P7_CODONTABLE  *ct     = p7_codontable_Create(gcode);
  P7_BG          *bgAA   = p7_bg_Create(abcAA);
  int             M      = esl_opt_GetInteger(go, "-M");
  int             N      = esl_opt_GetInteger(go, "-N");
  int             I      = esl_opt_GetInteger(go, "-I");

  p7_FLogsumInit();
  impl_Init();

  utest_chk(r, abcAA, gcode, bgAA, ct, M, N, I);

  p7_bg_Destroy(bgAA);
  p7_codontable_Destroy(ct);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  printf("All tests passed.\n");
  return eslOK;
}
#endif /*p7VITERBI_SP_CHK_TESTDRIVE*/
/*----------------- end, test driver ----------------------------*/
