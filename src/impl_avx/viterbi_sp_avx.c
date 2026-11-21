/* Spliced Viterbi algorithms (full matrix); AVX2 implementations for impl_avx dispatch.
 * These are the _avx-suffixed versions called by the runtime dispatcher in viterbi_sp.c.
 */

#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <stdio.h>
#include <math.h>

#include <immintrin.h>		/* AVX2 */

#include "easel.h"
#include "esl_avx.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_splice.h"
#include "impl_avx.h"

/* avx_rightshift_ps: right-shift a __m256 by one float lane, filling lane 0 with infv. */
static inline __m256
avx_rightshift_ps(__m256 v, __m256 infv)
{
  return _mm256_blend_ps(esl_avx_rightshiftz_float(v), infv, 0x01);
}

/* avx_hmax_ps: horizontal max across all 8 lanes of a __m256, return scalar. */
static inline void
avx_hmax_ps(__m256 v, float *ret_max)
{
  __m128 a = _mm256_extractf128_ps(v, 0);
  __m128 b = _mm256_extractf128_ps(v, 1);
  a = _mm_max_ps(a, b);
  a = _mm_max_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(2,3,0,1)));
  a = _mm_max_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(1,0,3,2)));
  _mm_store_ss(ret_max, a);
}

/* avx_FGetMDI: extract a single float from the AVX dp matrix at row i, model position k. */
static inline float
avx_FGetMDI(const P7_OMX *ox, int s, int i, int k, int Q)
{
  union { __m256 v; float p[8]; } u;
  u.v = ox->dpf_avx[i][p7X_NSCELLS * ((k-1) % Q) + s];
  return u.p[(k-1)/Q];
}

/* row_get(): one float DP cell from row <i> of a striped matrix given by
 * its row pointers <dp>: state <s> of model position <k>. */
static inline float
row_get(__m256 **dp, int s, int i, int k, int Q)
{
  union { __m256 v; float p[8]; } u;
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


/* Function:  p7_Viterbi_Spliced_avx()
 *
 * Purpose:   AVX2 implementation of translated spliced Viterbi.
 *            See p7_Viterbi_Spliced() documentation for details.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on bad inputs.
 */
int
p7_Viterbi_Spliced_avx(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_OMX *ox,
                        const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                        int i_start, int i_end, int min_intron, int global_start, int global_end)
{
  register __m256 mpv, dpv, ipv;
  register __m256 sv, msv, dcv;
  __m256   infv, b_entry;
  float    xB, xE;
  __m256  *dpc, *dpp3;
  __m256  *tp;
  __m256 **acc_ovx;
  __m256   ppv;
  __m256   tsc_p_vec;
  int      Q  = p7O_NQF_AVX(om_tr->M);
  int      L  = i_end - i_start + 1;
  int      C0;
  int      C1[p7P_MAXNUC + 1];
  int      C2[p7P_MAXNUC + 1];
  int      v, w, x;
  int      r, s, t, u;
  int      nuc1, nuc3;
  int      acc0, acc1, acc2;
  int      don0, don1, don2;
  int      don_sig;
  int      pv_i;
  int      pv_pi;
  int      i, q, j, ri, sub_i;
  __m256 **don_ovx      = don_ov->ivx_avx;
  __m256   sig_gtag, sig_gcag, sig_atac;

  if (om_tr->codon_lengths != 1)
    ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");
  if (ox->nscells != p7X_NSCELLS)
    ESL_EXCEPTION(eslEINVAL, "DP matrix must have p7X_NSCELLS cells per stripe position");

  ox->M              = om_tr->M;
  ox->L              = L;
  ox->has_own_scales = FALSE;
  ox->totscale       = 0.0f;
  infv               = _mm256_set1_ps(-eslINFINITY);
  tsc_p_vec          = _mm256_set1_ps(TSC_P);
  sig_gtag           = _mm256_set1_ps(signal_scores[p7S_GTAG]);
  sig_gcag           = _mm256_set1_ps(signal_scores[p7S_GCAG]);
  sig_atac           = _mm256_set1_ps(signal_scores[p7S_ATAC]);

  acc_ovx = acc_ov->ivx_avx;
  for (ri = 0; ri < SPLICE_ROWS; ri++)
    for (q = 0; q < Q; q++)
      acc_ovx[ri][q] = infv;

  for (j = 0; j < SIGNAL_MEM_SIZE; j++)
    for (q = 0; q < Q; q++)
      don_ovx[j][q] = infv;

  { __m256 *drow;
    for (ri = 0; ri <= 2 && ri <= L; ri++) {
      drow = ox->dpf_avx[ri];
      for (q = 0; q < Q; q++)
        MMO(drow, q) = DMO(drow, q) = IMO(drow, q) = infv;
    }
  }

  for (ri = 0; ri <= L; ri++) {
    ox->xmx[ri * p7X_NXCELLS + p7X_SCALE] = 0.0f;
    ox->xmx[ri * p7X_NXCELLS + p7X_E]     = -eslINFINITY;
    ox->xmx[ri * p7X_NXCELLS + p7X_N]     = -eslINFINITY;
    ox->xmx[ri * p7X_NXCELLS + p7X_J]     = -eslINFINITY;
    ox->xmx[ri * p7X_NXCELLS + p7X_B]     = -eslINFINITY;
    ox->xmx[ri * p7X_NXCELLS + p7X_C]     = -eslINFINITY;
  }
  xB = om_tr->xf[p7O_N][p7O_MOVE];
  ox->xmx[0 * p7X_NXCELLS + p7X_B] = xB;
  if (!global_start) ox->xmx[0 * p7X_NXCELLS + p7X_N] = 0.f;

  if (L >= 1) { ox->xmx[1 * p7X_NXCELLS + p7X_N] = 0.f; ox->xmx[1 * p7X_NXCELLS + p7X_B] = xB; }
  if (L >= 2) { ox->xmx[2 * p7X_NXCELLS + p7X_N] = 0.f; ox->xmx[2 * p7X_NXCELLS + p7X_B] = xB; }

  acc0 = acc1 = acc2 = -1;
  don0 = don1 = don2 = -1;
  v = w = p7P_MAXCODONS1;
  x = (sub_dsq[i_start] < p7P_MAXNUC) ? sub_dsq[i_start] : p7P_MAXCODONS1;
  if (L >= 2) {
    w = x;
    x = (sub_dsq[i_start + 1] < p7P_MAXNUC) ? sub_dsq[i_start + 1] : p7P_MAXCODONS1;
  }

  r = p7P_MAXCODONS1;
  s = (sub_dsq[i_start]     < p7P_MAXNUC) ? sub_dsq[i_start]     : p7P_MAXCODONS1;
  t = (sub_dsq[i_start + 1] < p7P_MAXNUC) ? sub_dsq[i_start + 1] : p7P_MAXCODONS1;
  u = (sub_dsq[i_start + 2] < p7P_MAXNUC) ? sub_dsq[i_start + 2] : p7P_MAXCODONS1;

  for (i = 3; i <= L; i++) {
    sub_i = i_start + i - 1;

    if (i >= min_intron + 3) {
      r = s; s = t; t = u;
      u = (sub_dsq[sub_i - min_intron + 1] < p7P_MAXNUC)
              ? sub_dsq[sub_i - min_intron + 1] : p7P_MAXCODONS1;
    }

    v = w;
    w = x;
    x = (sub_dsq[sub_i] < p7P_MAXNUC) ? sub_dsq[sub_i] : p7P_MAXCODONS1;

    C0 = p7P_CODON3_FS1(v, w, x);
    C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);

    for (nuc1 = 0; nuc1 < p7P_MAXNUC; nuc1++)
      C1[nuc1] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, w, x), p7P_DEGEN1_C);
    C1[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(p7P_MAXCODONS1, w, x), p7P_DEGEN1_C);
    nuc3 = ESL_MIN(x, p7P_MAXNUC);

    acc0 = acc1;
    acc1 = acc2;
    if      (SIGNAL(v, w) == ACCEPT_AG) acc2 = ACCEPT_AG;
    else if (SIGNAL(v, w) == ACCEPT_AC) acc2 = ACCEPT_AC;
    else                                acc2 = -1;

    if (i >= min_intron + 3) {
      don0 = don1; don1 = don2;
      if      (SIGNAL(t, u) == DONOR_GT) don2 = DONOR_GT;
      else if (SIGNAL(t, u) == DONOR_GC) don2 = DONOR_GC;
      else if (SIGNAL(t, u) == DONOR_AT) don2 = DONOR_AT;
      else                               don2 = -1;
    }

    pv_i  =  i    % 4;
    pv_pi = (i-3) % 4;
    dpc  = ox->dpf_avx[i];
    dpp3 = ox->dpf_avx[i - 3];

    if (!global_start) {
      float xN_im3 = ox->xmx[(i-3) * p7X_NXCELLS + p7X_N];
      float xN_i   = xN_im3 + om_tr->xf[p7O_N][p7O_LOOP];
      float xB_i   = xN_i   + om_tr->xf[p7O_N][p7O_MOVE];
      ox->xmx[i * p7X_NXCELLS + p7X_N] = xN_i;
      ox->xmx[i * p7X_NXCELLS + p7X_B] = xB_i;
    }

    mpv = avx_rightshift_ps(MMO(dpp3, Q - 1), infv);
    dpv = avx_rightshift_ps(DMO(dpp3, Q - 1), infv);
    ipv = avx_rightshift_ps(IMO(dpp3, Q - 1), infv);

    /* Pre-loop: compute P(i, k) for all stripes into acc_ovx. */
    for (q = 0; q < Q; q++) {
      __m256 psv = infv;
      if (acc0 == ACCEPT_AG) {
        psv = _mm256_max_ps(psv, _mm256_add_ps(
            _mm256_max_ps(_mm256_add_ps(don_ovx[p7S_GTAG][q], sig_gtag),
                          _mm256_add_ps(don_ovx[p7S_GCAG][q], sig_gcag)),
            om_tr->rfv_avx[C0][q]));
      } else if (acc0 == ACCEPT_AC) {
        psv = _mm256_max_ps(psv, _mm256_add_ps(
            _mm256_add_ps(don_ovx[p7S_ATAC][q], sig_atac),
            om_tr->rfv_avx[C0][q]));
      }
      if (acc1 == ACCEPT_AG) {
        for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
          psv = _mm256_max_ps(psv, _mm256_add_ps(
              _mm256_max_ps(_mm256_add_ps(don_ovx[SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + p7S_GTAG][q], sig_gtag),
                            _mm256_add_ps(don_ovx[SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + p7S_GCAG][q], sig_gcag)),
              om_tr->rfv_avx[C1[nuc1]][q]));
        }
      } else if (acc1 == ACCEPT_AC) {
        for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
          psv = _mm256_max_ps(psv, _mm256_add_ps(
              _mm256_add_ps(don_ovx[SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + p7S_ATAC][q], sig_atac),
              om_tr->rfv_avx[C1[nuc1]][q]));
        }
      }
      if (acc2 == ACCEPT_AG) {
        psv = _mm256_max_ps(psv,
            _mm256_max_ps(_mm256_add_ps(don_ovx[SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + p7S_GTAG][q], sig_gtag),
                          _mm256_add_ps(don_ovx[SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + p7S_GCAG][q], sig_gcag)));
      } else if (acc2 == ACCEPT_AC) {
        psv = _mm256_max_ps(psv,
            _mm256_add_ps(don_ovx[SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + p7S_ATAC][q], sig_atac));
      }
      acc_ovx[pv_i][q] = psv;
    }

    ppv = avx_rightshift_ps(acc_ovx[pv_pi][Q - 1], infv);

    tp  = om_tr->tfv_avx;
    dcv = infv;

    if (global_start) {
      if (i == 3) {
        union { __m256 v; float p[8]; } u2;
        u2.v    = infv;
        u2.p[0] = xB;
        b_entry = u2.v;
      } else {
        b_entry = infv;
      }
    } else {
      b_entry = _mm256_set1_ps(ox->xmx[(i-3) * p7X_NXCELLS + p7X_B]);
    }

    for (q = 0; q < Q; q++) {
      __m256 ppv_next = acc_ovx[pv_pi][q];

      tp++;                                                    /* skip BM */
      sv  =                  _mm256_add_ps(mpv, *tp); tp++;   /* MM */
      sv  = _mm256_max_ps(sv, _mm256_add_ps(ipv, *tp)); tp++; /* IM */
      sv  = _mm256_max_ps(sv, _mm256_add_ps(dpv, *tp)); tp++; /* DM */
      sv  = _mm256_max_ps(sv, _mm256_add_ps(ppv, tsc_p_vec)); /* PM */
      sv  = _mm256_max_ps(sv, b_entry);
      if (global_start) b_entry = infv;

      mpv = MMO(dpp3, q);
      dpv = DMO(dpp3, q);
      ipv = IMO(dpp3, q);

      msv = _mm256_add_ps(sv, om_tr->rfv_avx[C0][q]);
      MMO(dpc, q) = msv;
      DMO(dpc, q) = dcv;
      dcv = _mm256_add_ps(msv, *tp); tp++;                    /* MD */

      sv  = _mm256_add_ps(mpv, *tp); tp++;                    /* MI */
      sv  = _mm256_max_ps(sv, _mm256_add_ps(ipv, *tp)); tp++; /* II */
      { /* stop-codon mask */
        __m256 stop_mask = _mm256_cmp_ps(om_tr->rfv_avx[C0][q], infv, _CMP_EQ_OQ);
        sv = _mm256_or_ps(_mm256_andnot_ps(stop_mask, sv), _mm256_and_ps(stop_mask, infv));
      }
      IMO(dpc, q) = sv;

      ppv = ppv_next;
    }

    /* DD sweep */
    dcv = avx_rightshift_ps(dcv, infv);
    DMO(dpc, 0) = infv;
    tp = om_tr->tfv_avx + 7 * Q;
    for (q = 0; q < Q; q++) {
      DMO(dpc, q) = _mm256_max_ps(dcv, DMO(dpc, q));
      dcv = _mm256_add_ps(DMO(dpc, q), *tp); tp++;
    }
    if (om_tr->M < 100) {
      for (j = 1; j < 8; j++) {
        dcv = avx_rightshift_ps(dcv, infv);
        tp  = om_tr->tfv_avx + 7 * Q;
        for (q = 0; q < Q; q++) {
          DMO(dpc, q) = _mm256_max_ps(dcv, DMO(dpc, q));
          dcv = _mm256_add_ps(dcv, *tp); tp++;
        }
      }
    } else {
      for (j = 1; j < 8; j++) {
        register __m256 cv;
        dcv = avx_rightshift_ps(dcv, infv);
        tp  = om_tr->tfv_avx + 7 * Q;
        cv  = _mm256_setzero_ps();
        for (q = 0; q < Q; q++) {
          sv             = _mm256_max_ps(dcv, DMO(dpc, q));
          cv             = _mm256_or_ps(cv, _mm256_cmp_ps(sv, DMO(dpc, q), _CMP_GT_OQ));
          DMO(dpc, q)    = sv;
          dcv            = _mm256_add_ps(dcv, *tp); tp++;
        }
        if (!_mm256_movemask_ps(cv)) break;
      }
    }

    /* Donor site q-loops. */
    if (don2 >= 0 || don1 >= 0 || don0 >= 0) {
      __m256 *dpp_lb = ox->dpf_avx[i - min_intron - 3];
      __m256  lb_mpv, lb_dpv, lb_tsv;

      if (don2 >= 0) {
        for (nuc3 = 0; nuc3 < p7P_MAXNUC; nuc3++)
          C2[nuc3] = p7P_MINIDX(p7P_CODON3_FS1(r, s, nuc3), p7P_DEGEN1_C);
        C2[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(r, s, p7P_MAXCODONS1), p7P_DEGEN1_C);

        don_sig = (don2 == DONOR_GT) ? p7S_GTAG : (don2 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
        lb_mpv  = avx_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
        lb_dpv  = avx_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
        for (q = 0; q < Q; q++) {
          lb_tsv = _mm256_max_ps(lb_mpv, lb_dpv);
          for (nuc3 = 0; nuc3 <= p7P_MAXNUC; nuc3++) {
            int slot = SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + don_sig;
            don_ovx[slot][q] = _mm256_max_ps(don_ovx[slot][q],
                                              _mm256_add_ps(lb_tsv, om_tr->rfv_avx[C2[nuc3]][q]));
          }
          lb_mpv = MMO(dpp_lb, q);
          lb_dpv = DMO(dpp_lb, q);
        }
      }

      if (don1 >= 0) {
        nuc1    = ESL_MIN(r, p7P_MAXNUC);
        don_sig = (don1 == DONOR_GT) ? p7S_GTAG : (don1 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
        {
          int slot = SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + don_sig;
          lb_mpv = avx_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
          lb_dpv = avx_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
          for (q = 0; q < Q; q++) {
            lb_tsv = _mm256_max_ps(lb_mpv, lb_dpv);
            don_ovx[slot][q] = _mm256_max_ps(don_ovx[slot][q], lb_tsv);
            lb_mpv = MMO(dpp_lb, q);
            lb_dpv = DMO(dpp_lb, q);
          }
        }
      }

      if (don0 >= 0) {
        don_sig = (don0 == DONOR_GT) ? p7S_GTAG : (don0 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
        lb_mpv  = avx_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
        lb_dpv  = avx_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
        for (q = 0; q < Q; q++) {
          lb_tsv = _mm256_max_ps(lb_mpv, lb_dpv);
          don_ovx[don_sig][q] = _mm256_max_ps(don_ovx[don_sig][q], lb_tsv);
          lb_mpv = MMO(dpp_lb, q);
          lb_dpv = DMO(dpp_lb, q);
        }
      }
    }

    /* Semi-global exit */
    if (!global_end) {
      __m256 xEv = infv;
      for (q = 0; q < Q; q++)
        xEv = _mm256_max_ps(xEv, _mm256_max_ps(MMO(dpc, q), DMO(dpc, q)));
      xE = -eslINFINITY;
      avx_hmax_ps(xEv, &xE);
      ox->xmx[i * p7X_NXCELLS + p7X_E] = xE;
      ox->xmx[i * p7X_NXCELLS + p7X_C] = ESL_MAX(
          ox->xmx[(i-3) * p7X_NXCELLS + p7X_C] + om_tr->xf[p7O_C][p7O_LOOP],
          xE + om_tr->xf[p7O_E][p7O_MOVE]);
    }
  } /* end main loop i = 3..L */

  if (global_end) {
    union { __m256 v; float p[8]; } um, ud;
    int qM = (om_tr->M - 1) % Q;
    int rM = (om_tr->M - 1) / Q;
    um.v = MMO(ox->dpf_avx[L], qM);
    ud.v = DMO(ox->dpf_avx[L], qM);
    xE   = ESL_MAX(um.p[rM], ud.p[rM]);
    ox->xmx[L * p7X_NXCELLS + p7X_E] = xE;
    ox->xmx[L * p7X_NXCELLS + p7X_C] = xE + om_tr->xf[p7O_E][p7O_MOVE];
  }

  return eslOK;
}


/*****************************************************************
 * Checkpointed spliced Viterbi: p7_Viterbi_SplicedChk_avx() and
 * p7_Viterbi_SplicedTraceChk_avx()
 *
 * The fill computes the same rows as p7_Viterbi_Spliced_avx(), but into
 * a P7_CHKOMX that holds one block of rows at a time. At the start of
 * row b-3 for each block start b, the donor scores (don_ov) and the rows
 * they came from (don_idx, the working int vectors of the P7_CHKOIVX)
 * are saved as checkpoint b/blockR - 1; at row b the P7_CHKOMX saves the
 * saveR rows before b. To compute block b again, the trace restores
 * both, recomputes the acceptor P values and donor updates of rows
 * b-3..b-1 (their M, I, D rows are among the saved rows), and then
 * computes the rows of the block in full. The same row functions do
 * both, so the rows come out identical.
 *
 * saveR must be at least min_intron+6: row i reads row i-3, and its
 * donor update reads row i-min_intron-3.
 *****************************************************************/

/* Arguments shared by the checkpointed fill, the trace and the row
 * functions they both use. */
typedef struct {
  const ESL_DSQ        *sub_dsq;
  const P7_FS_OPROFILE *om_tr;
  P7_CHKOMX            *cx;
  P7_CHKOIVX           *cv;
  P7_OIVX              *acc_ov;
  P7_OIVX              *don_ov;
  const float          *signal_scores;
  int                   i_start;
  int                   M;
  int                   L;
  int                   Q;
  int                   min_intron;
  int                   global_start;
  int                   global_end;
} SPLICED_CHK;

/* spliced_init():
 * Clear the donor scores, donor rows and acceptor P values and set rows
 * 0..2, as at the start of p7_Viterbi_Spliced_avx(). With <do_xmx>, also
 * set the special states of every row. Block 0 must be loaded.
 */
static void
spliced_init(const SPLICED_CHK *sc, int do_xmx)
{
  P7_CHKOMX *cx      = sc->cx;
  __m256   **acc_ovx = sc->acc_ov->ivx_avx;
  __m256   **don_ovx = sc->don_ov->ivx_avx;
  __m256i  **don_idx = sc->cv->ivx_avx;
  __m256     infv    = _mm256_set1_ps(-eslINFINITY);
  __m256i    nonev   = _mm256_set1_epi32(-1);
  float      xB      = sc->om_tr->xf[p7O_N][p7O_MOVE];
  int        L       = sc->L;
  int        ri, q, j;

  for (ri = 0; ri < SPLICE_ROWS; ri++)
    for (q = 0; q < sc->Q; q++)
      acc_ovx[ri][q] = infv;

  for (j = 0; j < SIGNAL_MEM_SIZE; j++)
    for (q = 0; q < sc->Q; q++) {
      don_ovx[j][q] = infv;
      don_idx[j][q] = nonev;
    }

  for (ri = 0; ri <= 2 && ri <= L; ri++)
    for (q = 0; q < sc->Q; q++)
      MMO(cx->dpf_avx[ri], q) = DMO(cx->dpf_avx[ri], q) = IMO(cx->dpf_avx[ri], q) = infv;

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
      if (!sc->global_start) cx->xmx[0 * p7X_NXCELLS + p7X_N] = 0.f;

      if (L >= 1) { cx->xmx[1 * p7X_NXCELLS + p7X_N] = 0.f; cx->xmx[1 * p7X_NXCELLS + p7X_B] = xB; }
      if (L >= 2) { cx->xmx[2 * p7X_NXCELLS + p7X_N] = 0.f; cx->xmx[2 * p7X_NXCELLS + p7X_B] = xB; }
    }
}

/* spliced_acc_row():
 * P state values of row <i> (i >= 3), from the donor scores as they are
 * at the start of row <i>.
 */
static void
spliced_acc_row(const SPLICED_CHK *sc, int i)
{
  const P7_FS_OPROFILE *om_tr   = sc->om_tr;
  __m256              **acc_ovx = sc->acc_ov->ivx_avx;
  __m256              **don_ovx = sc->don_ov->ivx_avx;
  __m256                infv    = _mm256_set1_ps(-eslINFINITY);
  __m256                sig_gtag = _mm256_set1_ps(sc->signal_scores[p7S_GTAG]);
  __m256                sig_gcag = _mm256_set1_ps(sc->signal_scores[p7S_GCAG]);
  __m256                sig_atac = _mm256_set1_ps(sc->signal_scores[p7S_ATAC]);
  int                   pv_i    = i % SPLICE_ROWS;
  int                   acc0, acc1, acc2;
  int                   v, w, x;
  int                   C0;
  int                   C1[p7P_MAXNUC + 1];
  int                   nuc1, nuc3;
  int                   q;

  acc2 = acceptor_at(sc->sub_dsq, sc->i_start, i);
  acc1 = acceptor_at(sc->sub_dsq, sc->i_start, i-1);
  acc0 = acceptor_at(sc->sub_dsq, sc->i_start, i-2);

  v = nuc_at(sc->sub_dsq, sc->i_start, i-2);
  w = nuc_at(sc->sub_dsq, sc->i_start, i-1);
  x = nuc_at(sc->sub_dsq, sc->i_start, i);

  C0 = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);
  for (nuc1 = 0; nuc1 < p7P_MAXNUC; nuc1++)
    C1[nuc1] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, w, x), p7P_DEGEN1_C);
  C1[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(p7P_MAXCODONS1, w, x), p7P_DEGEN1_C);
  nuc3 = ESL_MIN(x, p7P_MAXNUC);

  for (q = 0; q < sc->Q; q++) {
    __m256 psv = infv;
    if (acc0 == ACCEPT_AG) {
      psv = _mm256_max_ps(psv, _mm256_add_ps(
          _mm256_max_ps(_mm256_add_ps(don_ovx[p7S_GTAG][q], sig_gtag),
                        _mm256_add_ps(don_ovx[p7S_GCAG][q], sig_gcag)),
          om_tr->rfv_avx[C0][q]));
    } else if (acc0 == ACCEPT_AC) {
      psv = _mm256_max_ps(psv, _mm256_add_ps(
          _mm256_add_ps(don_ovx[p7S_ATAC][q], sig_atac),
          om_tr->rfv_avx[C0][q]));
    }
    if (acc1 == ACCEPT_AG) {
      for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
        psv = _mm256_max_ps(psv, _mm256_add_ps(
            _mm256_max_ps(_mm256_add_ps(don_ovx[SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + p7S_GTAG][q], sig_gtag),
                          _mm256_add_ps(don_ovx[SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + p7S_GCAG][q], sig_gcag)),
            om_tr->rfv_avx[C1[nuc1]][q]));
      }
    } else if (acc1 == ACCEPT_AC) {
      for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
        psv = _mm256_max_ps(psv, _mm256_add_ps(
            _mm256_add_ps(don_ovx[SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + p7S_ATAC][q], sig_atac),
            om_tr->rfv_avx[C1[nuc1]][q]));
      }
    }
    if (acc2 == ACCEPT_AG) {
      psv = _mm256_max_ps(psv,
          _mm256_max_ps(_mm256_add_ps(don_ovx[SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + p7S_GTAG][q], sig_gtag),
                        _mm256_add_ps(don_ovx[SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + p7S_GCAG][q], sig_gcag)));
    } else if (acc2 == ACCEPT_AC) {
      psv = _mm256_max_ps(psv,
          _mm256_add_ps(don_ovx[SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + p7S_ATAC][q], sig_atac));
    }
    acc_ovx[pv_i][q] = psv;
  }
}

/* spliced_mdi_row():
 * M, I, D and special states of row <i> (i >= 3), as in
 * p7_Viterbi_Spliced_avx().
 */
static void
spliced_mdi_row(const SPLICED_CHK *sc, int i)
{
  const P7_FS_OPROFILE *om_tr   = sc->om_tr;
  P7_CHKOMX            *cx      = sc->cx;
  __m256              **acc_ovx = sc->acc_ov->ivx_avx;
  register __m256       mpv, dpv, ipv;
  register __m256       sv, msv, dcv;
  __m256                infv      = _mm256_set1_ps(-eslINFINITY);
  __m256                tsc_p_vec = _mm256_set1_ps(TSC_P);
  __m256                b_entry;
  __m256                ppv;
  __m256               *dpc  = cx->dpf_avx[i];
  __m256               *dpp3 = cx->dpf_avx[i-3];
  __m256               *tp;
  float                 xB   = om_tr->xf[p7O_N][p7O_MOVE];
  float                 xE;
  int                   Q     = sc->Q;
  int                   pv_pi = (i-3) % SPLICE_ROWS;
  int                   v, w, x;
  int                   C0;
  int                   q, j;

  v = nuc_at(sc->sub_dsq, sc->i_start, i-2);
  w = nuc_at(sc->sub_dsq, sc->i_start, i-1);
  x = nuc_at(sc->sub_dsq, sc->i_start, i);
  C0 = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);

  if (!sc->global_start) {
    float xN_im3 = cx->xmx[(i-3) * p7X_NXCELLS + p7X_N];
    float xN_i   = xN_im3 + om_tr->xf[p7O_N][p7O_LOOP];
    float xB_i   = xN_i   + om_tr->xf[p7O_N][p7O_MOVE];
    cx->xmx[i * p7X_NXCELLS + p7X_N] = xN_i;
    cx->xmx[i * p7X_NXCELLS + p7X_B] = xB_i;
  }

  mpv = avx_rightshift_ps(MMO(dpp3, Q - 1), infv);
  dpv = avx_rightshift_ps(DMO(dpp3, Q - 1), infv);
  ipv = avx_rightshift_ps(IMO(dpp3, Q - 1), infv);

  ppv = avx_rightshift_ps(acc_ovx[pv_pi][Q - 1], infv);

  tp  = om_tr->tfv_avx;
  dcv = infv;

  if (sc->global_start) {
    if (i == 3) {
      union { __m256 v; float p[8]; } u2;
      u2.v    = infv;
      u2.p[0] = xB;
      b_entry = u2.v;
    } else {
      b_entry = infv;
    }
  } else {
    b_entry = _mm256_set1_ps(cx->xmx[(i-3) * p7X_NXCELLS + p7X_B]);
  }

  for (q = 0; q < Q; q++) {
    __m256 ppv_next = acc_ovx[pv_pi][q];

    tp++;                                                    /* skip BM */
    sv  =                  _mm256_add_ps(mpv, *tp); tp++;   /* MM */
    sv  = _mm256_max_ps(sv, _mm256_add_ps(ipv, *tp)); tp++; /* IM */
    sv  = _mm256_max_ps(sv, _mm256_add_ps(dpv, *tp)); tp++; /* DM */
    sv  = _mm256_max_ps(sv, _mm256_add_ps(ppv, tsc_p_vec)); /* PM */
    sv  = _mm256_max_ps(sv, b_entry);
    if (sc->global_start) b_entry = infv;

    mpv = MMO(dpp3, q);
    dpv = DMO(dpp3, q);
    ipv = IMO(dpp3, q);

    msv = _mm256_add_ps(sv, om_tr->rfv_avx[C0][q]);
    MMO(dpc, q) = msv;
    DMO(dpc, q) = dcv;
    dcv = _mm256_add_ps(msv, *tp); tp++;                    /* MD */

    sv  = _mm256_add_ps(mpv, *tp); tp++;                    /* MI */
    sv  = _mm256_max_ps(sv, _mm256_add_ps(ipv, *tp)); tp++; /* II */
    { /* stop-codon mask */
      __m256 stop_mask = _mm256_cmp_ps(om_tr->rfv_avx[C0][q], infv, _CMP_EQ_OQ);
      sv = _mm256_or_ps(_mm256_andnot_ps(stop_mask, sv), _mm256_and_ps(stop_mask, infv));
    }
    IMO(dpc, q) = sv;

    ppv = ppv_next;
  }

  /* DD sweep */
  dcv = avx_rightshift_ps(dcv, infv);
  DMO(dpc, 0) = infv;
  tp = om_tr->tfv_avx + 7 * Q;
  for (q = 0; q < Q; q++) {
    DMO(dpc, q) = _mm256_max_ps(dcv, DMO(dpc, q));
    dcv = _mm256_add_ps(DMO(dpc, q), *tp); tp++;
  }
  if (om_tr->M < 100) {
    for (j = 1; j < 8; j++) {
      dcv = avx_rightshift_ps(dcv, infv);
      tp  = om_tr->tfv_avx + 7 * Q;
      for (q = 0; q < Q; q++) {
        DMO(dpc, q) = _mm256_max_ps(dcv, DMO(dpc, q));
        dcv = _mm256_add_ps(dcv, *tp); tp++;
      }
    }
  } else {
    for (j = 1; j < 8; j++) {
      register __m256 cv;
      dcv = avx_rightshift_ps(dcv, infv);
      tp  = om_tr->tfv_avx + 7 * Q;
      cv  = _mm256_setzero_ps();
      for (q = 0; q < Q; q++) {
        sv             = _mm256_max_ps(dcv, DMO(dpc, q));
        cv             = _mm256_or_ps(cv, _mm256_cmp_ps(sv, DMO(dpc, q), _CMP_GT_OQ));
        DMO(dpc, q)    = sv;
        dcv            = _mm256_add_ps(dcv, *tp); tp++;
      }
      if (!_mm256_movemask_ps(cv)) break;
    }
  }

  /* Semi-global exit */
  if (!sc->global_end) {
    __m256 xEv = infv;
    for (q = 0; q < Q; q++)
      xEv = _mm256_max_ps(xEv, _mm256_max_ps(MMO(dpc, q), DMO(dpc, q)));
    xE = -eslINFINITY;
    avx_hmax_ps(xEv, &xE);
    cx->xmx[i * p7X_NXCELLS + p7X_E] = xE;
    cx->xmx[i * p7X_NXCELLS + p7X_C] = ESL_MAX(
        cx->xmx[(i-3) * p7X_NXCELLS + p7X_C] + om_tr->xf[p7O_C][p7O_LOOP],
        xE + om_tr->xf[p7O_E][p7O_MOVE]);
  }
}

/* keep_donor_vec(): where <new_sc> is at least <*sc>, take <new_sc> and
 * record row <lbv> in <*idx>. Of equal scores the later row is kept. */
static inline void
keep_donor_vec(__m256 *sc, __m256i *idx, __m256 new_sc, __m256i lbv)
{
  __m256 mask = _mm256_cmp_ps(new_sc, *sc, _CMP_GE_OQ);
  *sc  = _mm256_blendv_ps(*sc, new_sc, mask);
  *idx = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(*idx), _mm256_castsi256_ps(lbv), mask));
}

/* spliced_donor_row():
 * Donor updates of row <i>: fold the look-back row i-min_intron-3 into
 * the donor scores (don_ov) at donor sites, and record in don_idx the
 * row each score came from.
 */
static void
spliced_donor_row(const SPLICED_CHK *sc, int i)
{
  const P7_FS_OPROFILE *om_tr   = sc->om_tr;
  __m256              **don_ovx = sc->don_ov->ivx_avx;
  __m256i             **don_idx = sc->cv->ivx_avx;
  __m256                infv    = _mm256_set1_ps(-eslINFINITY);
  __m256               *dpp_lb;
  __m256                lb_mpv, lb_dpv, lb_tsv;
  __m256i               lbv;
  int                   Q  = sc->Q;
  int                   lb = i - sc->min_intron - 3;
  int                   don0, don1, don2;
  int                   r, s;
  int                   C2[p7P_MAXNUC + 1];
  int                   nuc1, nuc3;
  int                   slot;
  int                   q;

  if (i < sc->min_intron + 3) return;

  don2 = donor_at(sc->sub_dsq, sc->i_start, sc->min_intron, i);
  don1 = donor_at(sc->sub_dsq, sc->i_start, sc->min_intron, i-1);
  don0 = donor_at(sc->sub_dsq, sc->i_start, sc->min_intron, i-2);
  if (don0 < 0 && don1 < 0 && don2 < 0) return;

  r      = nuc_at(sc->sub_dsq, sc->i_start, i - sc->min_intron - 2);
  s      = nuc_at(sc->sub_dsq, sc->i_start, i - sc->min_intron - 1);
  dpp_lb = sc->cx->dpf_avx[lb];
  lbv    = _mm256_set1_epi32(lb);

  if (don2 >= 0) {
    for (nuc3 = 0; nuc3 < p7P_MAXNUC; nuc3++)
      C2[nuc3] = p7P_MINIDX(p7P_CODON3_FS1(r, s, nuc3), p7P_DEGEN1_C);
    C2[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(r, s, p7P_MAXCODONS1), p7P_DEGEN1_C);

    lb_mpv = avx_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
    lb_dpv = avx_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
    for (q = 0; q < Q; q++) {
      lb_tsv = _mm256_max_ps(lb_mpv, lb_dpv);
      for (nuc3 = 0; nuc3 <= p7P_MAXNUC; nuc3++) {
        slot = SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + don2;
        keep_donor_vec(&don_ovx[slot][q], &don_idx[slot][q], _mm256_add_ps(lb_tsv, om_tr->rfv_avx[C2[nuc3]][q]), lbv);
      }
      lb_mpv = MMO(dpp_lb, q);
      lb_dpv = DMO(dpp_lb, q);
    }
  }

  if (don1 >= 0) {
    nuc1   = ESL_MIN(r, p7P_MAXNUC);
    slot   = SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + don1;
    lb_mpv = avx_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
    lb_dpv = avx_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
    for (q = 0; q < Q; q++) {
      lb_tsv = _mm256_max_ps(lb_mpv, lb_dpv);
      keep_donor_vec(&don_ovx[slot][q], &don_idx[slot][q], lb_tsv, lbv);
      lb_mpv = MMO(dpp_lb, q);
      lb_dpv = DMO(dpp_lb, q);
    }
  }

  if (don0 >= 0) {
    lb_mpv = avx_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
    lb_dpv = avx_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
    for (q = 0; q < Q; q++) {
      lb_tsv = _mm256_max_ps(lb_mpv, lb_dpv);
      keep_donor_vec(&don_ovx[don0][q], &don_idx[don0][q], lb_tsv, lbv);
      lb_mpv = MMO(dpp_lb, q);
      lb_dpv = DMO(dpp_lb, q);
    }
  }
}

/* spliced_rows():
 * Compute rows <from>..<to> (all >= 3): for each row the acceptor P
 * values, then (if <full>) the M, I, D and special states, then the
 * donor updates. With <full> FALSE, only the P values and donor updates
 * are done, for rows whose M, I, D values are already loaded.
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
 * p7_Viterbi_Spliced_avx(). Row L must be loaded. */
static void
spliced_global_end(const SPLICED_CHK *sc)
{
  union { __m256 v; float p[8]; } um, ud;
  int   qM = (sc->M - 1) % sc->Q;
  int   rM = (sc->M - 1) / sc->Q;
  float xE;

  um.v = MMO(sc->cx->dpf_avx[sc->L], qM);
  ud.v = DMO(sc->cx->dpf_avx[sc->L], qM);
  xE   = ESL_MAX(um.p[rM], ud.p[rM]);
  sc->cx->xmx[sc->L * p7X_NXCELLS + p7X_E] = xE;
  sc->cx->xmx[sc->L * p7X_NXCELLS + p7X_C] = xE + sc->om_tr->xf[p7O_E][p7O_MOVE];
}

/* spliced_check_sizes(): make sure the matrices and vectors fit an M by
 * L comparison with this <min_intron>. */
static int
spliced_check_sizes(const SPLICED_CHK *sc)
{
  if (sc->om_tr->codon_lengths != 1)                              ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");
  if (sc->cx->nscells != p7X_NSCELLS)                             ESL_EXCEPTION(eslEINVAL, "DP matrix must have p7X_NSCELLS cells per stripe position");
  if (sc->cx->saveR < sc->min_intron + 6)                         ESL_EXCEPTION(eslEINVAL, "checkpoints save %d rows; need %d", sc->cx->saveR, sc->min_intron + 6);
  if (sc->cx->allocQ4_avx < sc->Q || sc->cx->allocP < sc->L + 1 || sc->cx->allocXR < sc->L + 1) ESL_EXCEPTION(eslEINVAL, "checkpointed matrix too small");
  if (sc->cv->allocCF != SIGNAL_MEM_SIZE || sc->cv->allocCI != SIGNAL_MEM_SIZE) ESL_EXCEPTION(eslEINVAL, "checkpoint vectors need %d channels", SIGNAL_MEM_SIZE);
  if (sc->cv->allocQ4_avx < sc->Q || sc->cv->allocK < sc->L / sc->cx->blockR)   ESL_EXCEPTION(eslEINVAL, "checkpoint vectors too small");
  if (sc->don_ov->allocC < SIGNAL_MEM_SIZE || sc->don_ov->allocQ4_avx < sc->Q)  ESL_EXCEPTION(eslEINVAL, "donor vectors too small");
  if (sc->acc_ov->allocC < SPLICE_ROWS     || sc->acc_ov->allocQ4_avx < sc->Q)  ESL_EXCEPTION(eslEINVAL, "acceptor vectors too small");
  return eslOK;
}


/* Function:  p7_Viterbi_SplicedChk_avx()
 *
 * Purpose:   AVX2 checkpointed spliced Viterbi fill: same as
 *            p7_Viterbi_Spliced_avx(), but the rows go into the
 *            checkpointed matrix <cx>, with checkpoints of the donor
 *            scores and their rows in <cv>, so that
 *            p7_Viterbi_SplicedTraceChk_avx() can compute any block of
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
p7_Viterbi_SplicedChk_avx(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                          const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                          int i_start, int i_end, int min_intron, int global_start, int global_end)
{
  SPLICED_CHK sc = { sub_dsq, om_tr, cx, cv, acc_ov, don_ov, signal_scores, i_start,
                     om_tr->M, i_end - i_start + 1, p7O_NQF_AVX(om_tr->M), min_intron, global_start, global_end };
  int         blockR = cx->blockR;
  int         b, nb, end;
  int         status;

  if ((status = spliced_check_sizes(&sc)) != eslOK) return status;

  cx->M = sc.M;
  cx->L = sc.L;
  if ((status = p7_chkomx_StartBlock_avx(cx, 0)) != eslOK) return status;
  spliced_init(&sc, TRUE);

  for (b = 0; b <= sc.L; b += blockR)
    {
      if (b > 0 && (status = p7_chkomx_StartBlock_avx(cx, b)) != eslOK) return status;
      end = ESL_MIN(sc.L, b + blockR - 1);
      nb  = b + blockR;
      if (nb <= sc.L)
        {   /* save the donor state for block nb at the start of row nb-3 */
          spliced_rows(&sc, ESL_MAX(3, b), nb - 4, TRUE);
          if ((status = p7_chkoivx_Save_avx(cv, don_ov, sc.M, nb / blockR - 1)) != eslOK) return status;
          spliced_rows(&sc, nb - 3, end, TRUE);
        }
      else spliced_rows(&sc, ESL_MAX(3, b), end, TRUE);
    }

  if (global_end) spliced_global_end(&sc);
  return eslOK;
}


/* Function:  p7_Viterbi_SplicedTrace_avx()
 *
 * Purpose:   AVX2 implementation of spliced Viterbi traceback.
 *            See p7_Viterbi_SplicedTrace() documentation for details.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile is not 1-codon-length.
 *            <eslFAIL>   if traceback fails.
 */
int
p7_Viterbi_SplicedTrace_avx(const ESL_DSQ *sub_dsq, const P7_OMX *ox,
                              const P7_FS_PROFILE *gm_tr, const float *signal_scores,
                              P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                              int min_intron, float *vitsc, float *max_drop)
{
  float const *tsc    = gm_tr->tsc;
  float        r_tol  = 1e-5;
  float        a_tol  = 1e-4;
  int          M      = k_end - k_start + 1;
  int          L      = i_end - i_start + 1;
  int          Q      = p7O_NQF_AVX(ox->M);
  int          i      = L;
  int          k      = 0;
  int          j, z;
  int          t, u, v, w, x;
  int          d1, d2;
  int          c, c0, c1, c2, c3;
  int          sprv, scur;
  int          sub_i, sub_k;
  int          donor_i;
  int          acc[3] = {0, 0, 0};
  int          don_sig;
  int          P_cnt, P_idx;
  float        P_state;
  float        emit, emit0, emit1, emit2;
  float        vsc = 0.0f;
  float        curr_P_score;
  float        curr_M_score;
  float        max_M_score;
  float        curr_max_drop;
  float       *P_scores;
  int          status;

  P_scores = NULL;
  P_cnt    = 0;
  ESL_ALLOC(P_scores, sizeof(float) * M);

#define OMMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : avx_FGetMDI(ox, p7X_M, (ii), (kk), Q))
#define ODMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : avx_FGetMDI(ox, p7X_D, (ii), (kk), Q))
#define OIMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : avx_FGetMDI(ox, p7X_I, (ii), (kk), Q))
#define OXMXo(ii,s)  (ox->xmx[(ii) * p7X_NXCELLS + (s)])

  if (gm_tr->codon_lengths != 1) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");

  vsc = OXMXo(i, p7X_C) + gm_tr->xsc[p7P_C][p7P_MOVE];

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i+i_start-1, 0)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i+i_start-1, 0)) != eslOK) return status;

  sprv = p7T_C;

  while (sprv != p7T_S) {
    switch (sprv) {

    case p7T_C:
      if      (OXMXo(i, p7X_C) < OXMXo(i-2, p7X_C) || OXMXo(i, p7X_C) < OXMXo(i-1, p7X_C))                           scur = p7T_C;
      else if (OXMXo(i, p7X_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
      else if (esl_FCompare(OXMXo(i, p7X_C), OXMXo(i-3, p7X_C) + gm_tr->xsc[p7P_C][p7P_LOOP], r_tol, a_tol) == eslOK) scur = p7T_C;
      else if (esl_FCompare(OXMXo(i, p7X_C), OXMXo(i,   p7X_E) + gm_tr->xsc[p7P_E][p7P_MOVE], r_tol, a_tol) == eslOK) scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7T_E:
      if (OXMXo(i, p7X_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
      for (k = M; k >= 1; k--) {
        if (esl_FCompare(OXMXo(i, p7X_E), OMMo(i,k), r_tol, a_tol) == eslOK) { scur = p7T_M; break; }
        if (esl_FCompare(OXMXo(i, p7X_E), ODMo(i,k), r_tol, a_tol) == eslOK) { scur = p7T_D; break; }
      }
      if (k == 0) ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      break;

    case p7T_M:
      if (OMMo(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);

      sub_i = i_start + i - 1;
      v     = (sub_dsq[sub_i-2] < p7P_MAXNUC) ? (int)sub_dsq[sub_i-2] : p7P_MAXCODONS1;
      w     = (sub_dsq[sub_i-1] < p7P_MAXNUC) ? (int)sub_dsq[sub_i-1] : p7P_MAXCODONS1;
      x     = (sub_dsq[sub_i]   < p7P_MAXNUC) ? (int)sub_dsq[sub_i]   : p7P_MAXCODONS1;
      c3    = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);
      sub_k = k_start + k - 1;
      emit  = p7P_MSC_CODON(gm_tr, sub_k, c3);

      scur = -1;
      if      (esl_FCompare(OMMo(i,k), OMMo(i-3,k-1) + TSC(p7P_MM,sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(OMMo(i,k), OIMo(i-3,k-1) + TSC(p7P_IM,sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare(OMMo(i,k), ODMo(i-3,k-1) + TSC(p7P_DM,sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_D;
      else if (esl_FCompare(OMMo(i,k), OXMXo(i-3, p7X_B)                   + emit, r_tol, a_tol) == eslOK) scur = p7T_B;
      else {
        if (i < min_intron + 7) ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);

        vsc -= TSC_P;
        acc[0] = 0;
        if (SIGNAL(sub_dsq[sub_i-7], sub_dsq[sub_i-6]) == ACCEPT_AG) acc[0] = 1;
        if (SIGNAL(sub_dsq[sub_i-7], sub_dsq[sub_i-6]) == ACCEPT_AC) acc[0] = 2;
        acc[1] = 0;
        if (SIGNAL(sub_dsq[sub_i-6], sub_dsq[sub_i-5]) == ACCEPT_AG) acc[1] = 1;
        if (SIGNAL(sub_dsq[sub_i-6], sub_dsq[sub_i-5]) == ACCEPT_AC) acc[1] = 2;
        acc[2] = 0;
        if (SIGNAL(sub_dsq[sub_i-5], sub_dsq[sub_i-4]) == ACCEPT_AG) acc[2] = 1;
        if (SIGNAL(sub_dsq[sub_i-5], sub_dsq[sub_i-4]) == ACCEPT_AC) acc[2] = 2;

        if (acc[0] == 0 && acc[1] == 0 && acc[2] == 0)
          ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced: no acceptor site", k,i);

        for (j = 0; j < i - min_intron - 4; j++) {
          /* anything other than A,C,G,T gets a placeholder, so it can't form a donor signal */
          d1 = (sub_dsq[sub_i-min_intron-j-1] < p7P_MAXNUC) ? (int)sub_dsq[sub_i-min_intron-j-1] : p7P_MAXCODONS1;
          d2 = (sub_dsq[sub_i-min_intron-j]   < p7P_MAXNUC) ? (int)sub_dsq[sub_i-min_intron-j]   : p7P_MAXCODONS1;
          don_sig = -1;
          if      (SIGNAL(d1, d2) == DONOR_GT) don_sig = p7S_GTAG;
          else if (SIGNAL(d1, d2) == DONOR_GC) don_sig = p7S_GCAG;
          else if (SIGNAL(d1, d2) == DONOR_AT) don_sig = p7S_ATAC;

          if (don_sig != -1) {
            t = (sub_dsq[sub_i-min_intron-j-3] < p7P_MAXNUC) ? (int)sub_dsq[sub_i-min_intron-j-3] : p7P_MAXCODONS1;
            u = (sub_dsq[sub_i-min_intron-j-2] < p7P_MAXNUC) ? (int)sub_dsq[sub_i-min_intron-j-2] : p7P_MAXCODONS1;
            v = (sub_dsq[sub_i-5]              < p7P_MAXNUC) ? (int)sub_dsq[sub_i-5]              : p7P_MAXCODONS1;
            w = (sub_dsq[sub_i-4]              < p7P_MAXNUC) ? (int)sub_dsq[sub_i-4]              : p7P_MAXCODONS1;
            x = (sub_dsq[sub_i-3]              < p7P_MAXNUC) ? (int)sub_dsq[sub_i-3]              : p7P_MAXCODONS1;

            c2    = p7P_MINIDX(p7P_CODON3_FS1(t, u, x), p7P_DEGEN1_C);
            emit2 = p7P_MSC_CODON(gm_tr, sub_k-1, c2);
            c1    = p7P_MINIDX(p7P_CODON3_FS1(u, w, x), p7P_DEGEN1_C);
            emit1 = p7P_MSC_CODON(gm_tr, sub_k-1, c1);
            c0    = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);
            emit0 = p7P_MSC_CODON(gm_tr, sub_k-1, c0);

            if (don_sig == p7S_GTAG || don_sig == p7S_GCAG) {
              if (acc[2] == 1) {
                P_state = ESL_MAX(OMMo(i-min_intron-j-4,k-2), ODMo(i-min_intron-j-4,k-2)) + signal_scores[don_sig] + emit2;
                if (esl_FCompare(OMMo(i,k), P_state + TSC_P + emit, r_tol, a_tol) == eslOK) { scur = p7T_P; c = 2; donor_i = i-min_intron-j-4; vsc -= signal_scores[don_sig]; break; }
              }
              if (acc[1] == 1) {
                P_state = ESL_MAX(OMMo(i-min_intron-j-3,k-2), ODMo(i-min_intron-j-3,k-2)) + signal_scores[don_sig] + emit1;
                if (esl_FCompare(OMMo(i,k), P_state + TSC_P + emit, r_tol, a_tol) == eslOK) { scur = p7T_P; c = 1; donor_i = i-min_intron-j-3; vsc -= signal_scores[don_sig]; break; }
              }
              if (acc[0] == 1) {
                P_state = ESL_MAX(OMMo(i-min_intron-j-2,k-2), ODMo(i-min_intron-j-2,k-2)) + signal_scores[don_sig] + emit0;
                if (esl_FCompare(OMMo(i,k), P_state + TSC_P + emit, r_tol, a_tol) == eslOK) { scur = p7T_P; c = 0; donor_i = i-min_intron-j-2; vsc -= signal_scores[don_sig]; break; }
              }
            } else { /* p7S_ATAC */
              if (acc[2] == 2) {
                P_state = ESL_MAX(OMMo(i-min_intron-j-4,k-2), ODMo(i-min_intron-j-4,k-2)) + signal_scores[don_sig] + emit2;
                if (esl_FCompare(OMMo(i,k), P_state + TSC_P + emit, r_tol, a_tol) == eslOK) { scur = p7T_P; c = 2; donor_i = i-min_intron-j-4; vsc -= signal_scores[don_sig]; break; }
              }
              if (acc[1] == 2) {
                P_state = ESL_MAX(OMMo(i-min_intron-j-3,k-2), ODMo(i-min_intron-j-3,k-2)) + signal_scores[don_sig] + emit1;
                if (esl_FCompare(OMMo(i,k), P_state + TSC_P + emit, r_tol, a_tol) == eslOK) { scur = p7T_P; c = 1; donor_i = i-min_intron-j-3; vsc -= signal_scores[don_sig]; break; }
              }
              if (acc[0] == 2) {
                P_state = ESL_MAX(OMMo(i-min_intron-j-2,k-2), ODMo(i-min_intron-j-2,k-2)) + signal_scores[don_sig] + emit0;
                if (esl_FCompare(OMMo(i,k), P_state + TSC_P + emit, r_tol, a_tol) == eslOK) { scur = p7T_P; c = 0; donor_i = i-min_intron-j-2; vsc -= signal_scores[don_sig]; break; }
              }
            }
          }
        }
        if (scur == -1) ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);

        /* Track P state scores for x_drop calculation */
        P_scores[P_cnt++] = TSC_P + signal_scores[don_sig];
      }

      k--; i -= 3;
      break;

    case p7T_D:
      if (ODMo(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);
      sub_k = k_start + k - 1;
      if      (esl_FCompare(ODMo(i,k), OMMo(i,k-1) + TSC(p7P_MD,sub_k-1), r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(ODMo(i,k), ODMo(i,k-1) + TSC(p7P_DD,sub_k-1), r_tol, a_tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      k--;
      break;

    case p7T_I:
      if (OIMo(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);
      sub_k = k_start + k - 1;
      if      (esl_FCompare(OIMo(i,k), OMMo(i-3,k) + TSC(p7P_MI,sub_k), r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(OIMo(i,k), OIMo(i-3,k) + TSC(p7P_II,sub_k), r_tol, a_tol) == eslOK) scur = p7T_I;
      else ESL_EXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      i -= 3;
      break;

    case p7T_P:
      scur = (OMMo(donor_i,k-1) >= ODMo(donor_i,k-1)) ? p7T_M : p7T_D;
      k--; i = donor_i;
      break;

    case p7T_N:
      if (i > 0 && OXMXo(i, p7X_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = (i == 0) ? p7T_S : p7T_N;
      break;

    case p7T_B:
      vsc += p7P_TSC(gm_tr, k, p7P_BM);
      if (OXMXo(i, p7X_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);
      if (i == 0) { scur = p7T_N; break; }
      if (esl_FCompare(OXMXo(i, p7X_B), OXMXo(i, p7X_N) + gm_tr->xsc[p7P_N][p7P_MOVE], r_tol, a_tol) == eslOK) scur = p7T_N;
      else ESL_EXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    }

    if      (scur == p7T_M) c = 3;
    else if (scur != p7T_P) c = 0;

    if ((status = p7_trace_fs_Append(tr, scur, k_start+k-1, i_start+i-1, c)) != eslOK) return status;

    if ((scur == p7T_N || scur == p7T_C) && scur == sprv) i--;

    sprv = scur;
  }

  tr->M = M;
  tr->L = L;

  if (vitsc != NULL) *vitsc = vsc;

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
        curr_M_score = OMMo(i, k) - curr_P_score;
        if (max_M_score != -eslINFINITY)
          curr_max_drop = ESL_MAX(curr_max_drop, max_M_score - curr_M_score);
        max_M_score = ESL_MAX(max_M_score, curr_M_score);
      }
    }
    *max_drop = curr_max_drop;
  }

  free(P_scores);

#undef OMMo
#undef ODMo
#undef OIMo
#undef OXMXo

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
replay_donors_from(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, __m256 **dp, int Q, const float *signal_scores,
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
 * replay_donors_from() for the full matrix <ox> of p7_Viterbi_Spliced_avx(),
 * replaying from the first row that has donor updates.
 */
static int
replay_donors(const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, const P7_OMX *ox, const float *signal_scores,
              int i_start, int k_start, int min_intron, int acc_i, int kk,
              float *ret_sc, int *ret_donor_i, int *ret_c, int *ret_don_sig)
{
  DONOR_SLOTS ds;

  donor_slots_clear(&ds);
  return replay_donors_from(sub_dsq, gm_tr, ox->dpf_avx, p7O_NQF_AVX(ox->M), signal_scores, i_start, k_start, min_intron, min_intron + 3, &ds,
                            acc_i, kk, ret_sc, ret_donor_i, ret_c, ret_don_sig);
}


/* Function:  p7_Viterbi_SplicedTraceReplay_avx()
 *
 * Purpose:   AVX2 implementation of the spliced Viterbi traceback that
 *            replays donor updates to find each P state's donor.
 *            See p7_GViterbi_SplicedTrace() for details.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile is not 1-codon-length.
 *            <eslFAIL>   if traceback fails.
 */
int
p7_Viterbi_SplicedTraceReplay_avx(const ESL_DSQ *sub_dsq, const P7_OMX *ox,
                              const P7_FS_PROFILE *gm_tr, const float *signal_scores,
                              P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                              int min_intron, float *vitsc, float *max_drop)
{
  float const *tsc    = gm_tr->tsc;
  float        r_tol  = 1e-5;
  float        a_tol  = 1e-4;
  int          M      = k_end - k_start + 1;
  int          L      = i_end - i_start + 1;
  int          Q      = p7O_NQF_AVX(ox->M);
  int          i      = L;
  int          k      = 0;
  int          z;
  int          v, w, x;
  int          c, c3;
  int          sprv, scur;
  int          sub_i, sub_k;
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
  int          status;

  P_scores = NULL;
  P_cnt    = 0;
  ESL_ALLOC(P_scores, sizeof(float) * M);

#define OMMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : avx_FGetMDI(ox, p7X_M, (ii), (kk), Q))
#define ODMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : avx_FGetMDI(ox, p7X_D, (ii), (kk), Q))
#define OIMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : avx_FGetMDI(ox, p7X_I, (ii), (kk), Q))
#define OXMXo(ii,s)  (ox->xmx[(ii) * p7X_NXCELLS + (s)])

  if (gm_tr->codon_lengths != 1) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");

  vsc = OXMXo(i, p7X_C) + gm_tr->xsc[p7P_C][p7P_MOVE];

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i+i_start-1, 0)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i+i_start-1, 0)) != eslOK) return status;

  sprv = p7T_C;

  while (sprv != p7T_S) {
    switch (sprv) {

    case p7T_C:
      if      (OXMXo(i, p7X_C) < OXMXo(i-2, p7X_C) || OXMXo(i, p7X_C) < OXMXo(i-1, p7X_C))                           scur = p7T_C;
      else if (OXMXo(i, p7X_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
      else if (esl_FCompare(OXMXo(i, p7X_C), OXMXo(i-3, p7X_C) + gm_tr->xsc[p7P_C][p7P_LOOP], r_tol, a_tol) == eslOK) scur = p7T_C;
      else if (esl_FCompare(OXMXo(i, p7X_C), OXMXo(i,   p7X_E) + gm_tr->xsc[p7P_E][p7P_MOVE], r_tol, a_tol) == eslOK) scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7T_E:
      if (OXMXo(i, p7X_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
      for (k = M; k >= 1; k--) {
        if (esl_FCompare(OXMXo(i, p7X_E), OMMo(i,k), r_tol, a_tol) == eslOK) { scur = p7T_M; break; }
        if (esl_FCompare(OXMXo(i, p7X_E), ODMo(i,k), r_tol, a_tol) == eslOK) { scur = p7T_D; break; }
      }
      if (k == 0) ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      break;

    case p7T_M:
      if (OMMo(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);

      sub_i = i_start + i - 1;
      v     = (sub_dsq[sub_i-2] < p7P_MAXNUC) ? (int)sub_dsq[sub_i-2] : p7P_MAXCODONS1;
      w     = (sub_dsq[sub_i-1] < p7P_MAXNUC) ? (int)sub_dsq[sub_i-1] : p7P_MAXCODONS1;
      x     = (sub_dsq[sub_i]   < p7P_MAXNUC) ? (int)sub_dsq[sub_i]   : p7P_MAXCODONS1;
      c3    = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);
      sub_k = k_start + k - 1;
      emit  = p7P_MSC_CODON(gm_tr, sub_k, c3);

      scur = -1;
      if      (esl_FCompare(OMMo(i,k), OMMo(i-3,k-1) + TSC(p7P_MM,sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(OMMo(i,k), OIMo(i-3,k-1) + TSC(p7P_IM,sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare(OMMo(i,k), ODMo(i-3,k-1) + TSC(p7P_DM,sub_k-1) + emit, r_tol, a_tol) == eslOK) scur = p7T_D;
      else if (esl_FCompare(OMMo(i,k), OXMXo(i-3, p7X_B)                   + emit, r_tol, a_tol) == eslOK) scur = p7T_B;
      else {    /* The P state is not stored - recreate it from the donor updates */
        if (replay_donors(sub_dsq, gm_tr, ox, signal_scores, i_start, k_start, min_intron, i-3, k-1, &P_state, &donor_i, &c, &don_sig) != eslOK)
          ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);
        if (esl_FCompare(OMMo(i,k), P_state + TSC_P + emit, r_tol, a_tol) != eslOK)
          ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);

        scur = p7T_P;
        vsc -= TSC_P;
        vsc -= signal_scores[don_sig];

        /* Track P state scores for x_drop calculation */
        P_scores[P_cnt++] = TSC_P + signal_scores[don_sig];
      }

      k--; i -= 3;
      break;

    case p7T_D:
      if (ODMo(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);
      sub_k = k_start + k - 1;
      if      (esl_FCompare(ODMo(i,k), OMMo(i,k-1) + TSC(p7P_MD,sub_k-1), r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(ODMo(i,k), ODMo(i,k-1) + TSC(p7P_DD,sub_k-1), r_tol, a_tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      k--;
      break;

    case p7T_I:
      if (OIMo(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);
      sub_k = k_start + k - 1;
      if      (esl_FCompare(OIMo(i,k), OMMo(i-3,k) + TSC(p7P_MI,sub_k), r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(OIMo(i,k), OIMo(i-3,k) + TSC(p7P_II,sub_k), r_tol, a_tol) == eslOK) scur = p7T_I;
      else ESL_EXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      i -= 3;
      break;

    case p7T_P:
      scur = (OMMo(donor_i,k-1) >= ODMo(donor_i,k-1)) ? p7T_M : p7T_D;
      k--; i = donor_i;
      break;

    case p7T_N:
      if (i > 0 && OXMXo(i, p7X_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = (i == 0) ? p7T_S : p7T_N;
      break;

    case p7T_B:
      vsc += p7P_TSC(gm_tr, k, p7P_BM);
      if (OXMXo(i, p7X_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);
      if (i == 0) { scur = p7T_N; break; }
      if (esl_FCompare(OXMXo(i, p7X_B), OXMXo(i, p7X_N) + gm_tr->xsc[p7P_N][p7P_MOVE], r_tol, a_tol) == eslOK) scur = p7T_N;
      else ESL_EXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    }

    if      (scur == p7T_M) c = 3;
    else if (scur != p7T_P) c = 0;

    if ((status = p7_trace_fs_Append(tr, scur, k_start+k-1, i_start+i-1, c)) != eslOK) return status;

    if ((scur == p7T_N || scur == p7T_C) && scur == sprv) i--;

    sprv = scur;
  }

  tr->M = M;
  tr->L = L;

  if (vitsc != NULL) *vitsc = vsc;

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
        curr_M_score = OMMo(i, k) - curr_P_score;
        if (max_M_score != -eslINFINITY)
          curr_max_drop = ESL_MAX(curr_max_drop, max_M_score - curr_M_score);
        max_M_score = ESL_MAX(max_M_score, curr_M_score);
      }
    }
    *max_drop = curr_max_drop;
  }

  free(P_scores);

#undef OMMo
#undef ODMo
#undef OIMo
#undef OXMXo

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

  if ((status = p7_chkomx_LoadBlock_avx(sc->cx, b)) != eslOK) return status;

  if (b == 0)
    {
      spliced_init(sc, FALSE);
      spliced_rows(sc, 3, end, TRUE);
    }
  else
    {
      if ((status = p7_chkoivx_Load_avx(sc->cv, sc->don_ov, sc->M, b / sc->cx->blockR - 1)) != eslOK) return status;
      spliced_rows(sc, b - 3, b - 1, FALSE);
      spliced_rows(sc, b, end, TRUE);
    }

  if (sc->global_end && end == sc->L) spliced_global_end(sc);
  return eslOK;
}

/* spliced_need_row():
 * Make sure row <i> and the rows it reads are loaded: if <i> isn't in
 * the loaded block, load the block that contains it.
 */
static int
spliced_need_row(const SPLICED_CHK *sc, int i)
{
  int b = (i / sc->cx->blockR) * sc->cx->blockR;

  if (b == sc->cx->b) return eslOK;
  return spliced_load_block(sc, b);
}

/* ck_lane(): value of model position <kk> in channel <c> of checkpoint
 * <ck> of <cv>, as float <ret_f> (float channels) or int <ret_i> (int
 * channels, c >= allocCF). Layout as in P7_CHKOIVX: [ck][c][stripe].
 */
static inline void
ck_lane(const P7_CHKOIVX *cv, int ck, int c, int kk, int Q, float *ret_f, int *ret_i)
{
  __m256 *base = (__m256 *) (((unsigned long int) ((char *) cv->ck_mem_avx + 31)) & (~0x1f));
  __m256 *vec  = base + ((size_t) ck * (cv->allocCF + cv->allocCI) + c) * cv->allocQ4_avx + (kk-1) % Q;
  union { __m256 v; float p[8]; int32_t n[8]; } u;

  u.v = *vec;
  if (ret_f) *ret_f = u.p[(kk-1)/Q];
  if (ret_i) *ret_i = u.n[(kk-1)/Q];
}

/* replay_donors_chk():
 * replay_donors_from() for the loaded block of <sc->cx>, replaying from
 * the start of that block's checkpoint (row b-3), or from the first row
 * with donor updates for block 0. <acc_i> must be >= b-3.
 */
static int
replay_donors_chk(const SPLICED_CHK *sc, const P7_FS_PROFILE *gm_tr, int k_start, int acc_i, int kk,
                  float *ret_sc, int *ret_donor_i, int *ret_c, int *ret_don_sig)
{
  const P7_CHKOIVX *cv = sc->cv;
  int               b  = sc->cx->b;
  DONOR_SLOTS       ds;
  int               s0;
  int               ck;
  int               nuc, sig, slot;

  if (b == 0)
    {
      donor_slots_clear(&ds);
      s0 = sc->min_intron + 3;
    }
  else
    {
      ck = b / sc->cx->blockR - 1;
      for (sig = 0; sig < p7S_SPLICE_SIGNALS; sig++) {
        ck_lane(cv, ck, sig,                  kk, sc->Q, &ds.sc0[sig], NULL);
        ck_lane(cv, ck, cv->allocCF + sig,    kk, sc->Q, NULL, &ds.row0[sig]);
        for (nuc = 0; nuc <= p7P_MAXNUC; nuc++) {
          slot = SPLICE_OFFSET_1 + nuc * p7S_SPLICE_SIGNALS + sig;
          ck_lane(cv, ck, slot,               kk, sc->Q, &ds.sc1[nuc][sig], NULL);
          ck_lane(cv, ck, cv->allocCF + slot, kk, sc->Q, NULL, &ds.row1[nuc][sig]);
          slot = SPLICE_OFFSET_2 + nuc * p7S_SPLICE_SIGNALS + sig;
          ck_lane(cv, ck, slot,               kk, sc->Q, &ds.sc2[nuc][sig], NULL);
          ck_lane(cv, ck, cv->allocCF + slot, kk, sc->Q, NULL, &ds.row2[nuc][sig]);
        }
      }
      s0 = b - 3;
    }

  return replay_donors_from(sc->sub_dsq, gm_tr, sc->cx->dpf_avx, sc->Q, sc->signal_scores, sc->i_start, k_start, sc->min_intron,
                            s0, &ds, acc_i, kk, ret_sc, ret_donor_i, ret_c, ret_don_sig);
}


/* Function:  p7_Viterbi_SplicedTraceChk_avx()
 *
 * Purpose:   AVX2 traceback of a matrix filled by
 *            p7_Viterbi_SplicedChk_avx() with the same arguments. Same
 *            as p7_Viterbi_SplicedTraceReplay_avx(), except that the
 *            trace runs from the last block back to the first, computing
 *            each block it passes through again from its checkpoint and
 *            skipping blocks that lie inside introns. Each P state's
 *            donor comes from replaying donor updates from the checkpoint
 *            of the block that holds the acceptor. The M state scores on
 *            the trace are recorded as it goes, for <max_drop>.
 *
 *            <om_tr>, <acc_ov>, <don_ov>, <global_start> and <global_end>
 *            are as passed to p7_Viterbi_SplicedChk_avx(), for computing
 *            blocks again.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if the profile is not 1-codon-length or <cx>
 *                        isn't filled for this comparison.
 *            <eslFAIL>   if traceback fails.
 */
int
p7_Viterbi_SplicedTraceChk_avx(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                               const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                               P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                               int min_intron, int global_start, int global_end, float *vitsc, float *max_drop)
{
  SPLICED_CHK  sc     = { sub_dsq, om_tr, cx, cv, acc_ov, don_ov, signal_scores, i_start,
                          om_tr->M, i_end - i_start + 1, p7O_NQF_AVX(om_tr->M), min_intron, global_start, global_end };
  float const *tsc    = gm_tr->tsc;
  float        r_tol  = 1e-5;
  float        a_tol  = 1e-4;
  int          M      = k_end - k_start + 1;
  int          L      = i_end - i_start + 1;
  int          i      = L;
  int          k      = 0;
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
  void        *p;
  int          status;

  P_scores = NULL;
  M_scores = NULL;
  P_cnt    = 0;
  M_alloc  = 256;
  ESL_ALLOC(P_scores, sizeof(float) * M);
  ESL_ALLOC(M_scores, sizeof(float) * M_alloc);

#define OMMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : p7_chkomx_FGetMDI_avx(cx, p7X_M, (ii), (kk)))
#define ODMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : p7_chkomx_FGetMDI_avx(cx, p7X_D, (ii), (kk)))
#define OIMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : p7_chkomx_FGetMDI_avx(cx, p7X_I, (ii), (kk)))
#define OXMXo(ii,s)  (cx->xmx[(ii) * p7X_NXCELLS + (s)])

  if (gm_tr->codon_lengths != 1) ESL_XEXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");
  if (cx->M != M || cx->L != L || cx->b != (L / cx->blockR) * cx->blockR) ESL_XEXCEPTION(eslEINVAL, "matrix isn't filled for this comparison");

  vsc = OXMXo(i, p7X_C) + gm_tr->xsc[p7P_C][p7P_MOVE];

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i+i_start-1, 0)) != eslOK) goto ERROR;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i+i_start-1, 0)) != eslOK) goto ERROR;
  M_scores[0] = M_scores[1] = 0.;

  sprv = p7T_C;

  while (sprv != p7T_S) {
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
      if ((status = spliced_need_row(&sc, i)) != eslOK) goto ERROR;
      for (k = M; k >= 1; k--) {
        if (esl_FCompare(OXMXo(i, p7X_E), OMMo(i,k), r_tol, a_tol) == eslOK) { scur = p7T_M; break; }
        if (esl_FCompare(OXMXo(i, p7X_E), ODMo(i,k), r_tol, a_tol) == eslOK) { scur = p7T_D; break; }
      }
      if (k == 0) ESL_XEXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      break;

    case p7T_M:
      if ((status = spliced_need_row(&sc, i)) != eslOK) goto ERROR;
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
        if (replay_donors_chk(&sc, gm_tr, k_start, i-3, k-1, &P_state, &donor_i, &c, &don_sig) != eslOK)
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
      if ((status = spliced_need_row(&sc, i)) != eslOK) goto ERROR;
      if (ODMo(i,k) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);
      sub_k = k_start + k - 1;
      if      (esl_FCompare(ODMo(i,k), OMMo(i,k-1) + TSC(p7P_MD,sub_k-1), r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(ODMo(i,k), ODMo(i,k-1) + TSC(p7P_DD,sub_k-1), r_tol, a_tol) == eslOK) scur = p7T_D;
      else ESL_XEXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      k--;
      break;

    case p7T_I:
      if ((status = spliced_need_row(&sc, i)) != eslOK) goto ERROR;
      if (OIMo(i,k) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);
      sub_k = k_start + k - 1;
      if      (esl_FCompare(OIMo(i,k), OMMo(i-3,k) + TSC(p7P_MI,sub_k), r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(OIMo(i,k), OIMo(i-3,k) + TSC(p7P_II,sub_k), r_tol, a_tol) == eslOK) scur = p7T_I;
      else ESL_XEXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      i -= 3;
      break;

    case p7T_P:
      if ((status = spliced_need_row(&sc, donor_i)) != eslOK) goto ERROR;
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

#undef TSC_P

#endif /* eslENABLE_AVX */
