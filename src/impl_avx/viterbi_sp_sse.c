/* Spliced Viterbi algorithms (full matrix); SSE implementations for impl_avx dispatch.
 * These are the _sse-suffixed versions called by the runtime dispatcher in viterbi_sp.c.
 */

#include "p7_config.h"

#ifdef eslENABLE_SSE

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>
#include <emmintrin.h>

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_splice.h"
#include "impl_avx.h"

/* Log-space constant for the P->M transition (= logf(4.58e-5)). */
#define TSC_P logf(4.58e-5f)


/* Function:  p7_Viterbi_Spliced_sse()
 *
 * Purpose:   SSE implementation of translated spliced Viterbi.
 *            See p7_Viterbi_Spliced() documentation for details.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on bad inputs.
 */
int
p7_Viterbi_Spliced_sse(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_OMX *ox,
                       const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                       int i_start, int i_end, int min_intron, int global_start, int global_end)
{
  register __m128 mpv, dpv, ipv;
  register __m128 sv, msv, dcv;
  __m128   infv, b_entry;
  float    xB, xE;
  __m128  *dpc, *dpp3;
  __m128  *tp;
  __m128 **acc_ovx;
  __m128   ppv;
  __m128   tsc_p_vec;
  int      Q  = p7O_NQF(om_tr->M);
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
  __m128 **don_ovx      = don_ov->ivx;
  __m128   sig_gtag, sig_gcag, sig_atac;

  if (om_tr->codon_lengths != 1)
    ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");
  if (ox->nscells != p7X_NSCELLS)
    ESL_EXCEPTION(eslEINVAL, "DP matrix must have p7X_NSCELLS cells per stripe position");

  ox->M              = om_tr->M;
  ox->L              = L;
  ox->has_own_scales = FALSE;
  ox->totscale       = 0.0f;
  infv               = _mm_set1_ps(-eslINFINITY);
  tsc_p_vec          = _mm_set1_ps(TSC_P);
  sig_gtag           = _mm_set1_ps(signal_scores[p7S_GTAG]);
  sig_gcag           = _mm_set1_ps(signal_scores[p7S_GCAG]);
  sig_atac           = _mm_set1_ps(signal_scores[p7S_ATAC]);

  acc_ovx = acc_ov->ivx;
  for (ri = 0; ri < SPLICE_ROWS; ri++)
    for (q = 0; q < Q; q++)
      acc_ovx[ri][q] = infv;

  for (j = 0; j < SIGNAL_MEM_SIZE; j++)
    for (q = 0; q < Q; q++)
      don_ovx[j][q] = infv;

  for (ri = 0; ri <= L; ri++)
    for (q = 0; q < Q; q++)
      MMO(ox->dpf[ri], q) = DMO(ox->dpf[ri], q) = IMO(ox->dpf[ri], q) = infv;

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
    dpc  = ox->dpf[i];
    dpp3 = ox->dpf[i - 3];

    if (!global_start) {
      float xN_im3 = ox->xmx[(i-3) * p7X_NXCELLS + p7X_N];
      float xN_i   = xN_im3 + om_tr->xf[p7O_N][p7O_LOOP];
      float xB_i   = xN_i   + om_tr->xf[p7O_N][p7O_MOVE];
      ox->xmx[i * p7X_NXCELLS + p7X_N] = xN_i;
      ox->xmx[i * p7X_NXCELLS + p7X_B] = xB_i;
    }

    mpv = esl_sse_rightshift_ps(MMO(dpp3, Q - 1), infv);
    dpv = esl_sse_rightshift_ps(DMO(dpp3, Q - 1), infv);
    ipv = esl_sse_rightshift_ps(IMO(dpp3, Q - 1), infv);

    /* Pre-loop: compute P(i, k) for all stripes into acc_ovx. */
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
      b_entry = _mm_set1_ps(ox->xmx[(i-3) * p7X_NXCELLS + p7X_B]);
    }

    for (q = 0; q < Q; q++) {
      __m128 ppv_next = acc_ovx[pv_pi][q];

      tp++;                                               /* skip BM */
      sv  =                _mm_add_ps(mpv, *tp); tp++;   /* MM */
      sv  = _mm_max_ps(sv, _mm_add_ps(ipv, *tp)); tp++;  /* IM */
      sv  = _mm_max_ps(sv, _mm_add_ps(dpv, *tp)); tp++;  /* DM */
      sv  = _mm_max_ps(sv, _mm_add_ps(ppv, tsc_p_vec));  /* PM */
      sv  = _mm_max_ps(sv, b_entry);
      if (global_start) b_entry = infv;

      mpv = MMO(dpp3, q);
      dpv = DMO(dpp3, q);
      ipv = IMO(dpp3, q);

      msv = _mm_add_ps(sv, om_tr->rfv[C0][q]);
      MMO(dpc, q) = msv;
      DMO(dpc, q) = dcv;
      dcv = _mm_add_ps(msv, *tp); tp++;                  /* MD */

      sv  = _mm_add_ps(mpv, *tp); tp++;                  /* MI */
      sv  = _mm_max_ps(sv, _mm_add_ps(ipv, *tp)); tp++;  /* II */
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
        register __m128 cv;
        dcv = esl_sse_rightshift_ps(dcv, infv);
        tp  = om_tr->tfv + 7 * Q;
        cv  = _mm_setzero_ps();
        for (q = 0; q < Q; q++) {
          sv             = _mm_max_ps(dcv, DMO(dpc, q));
          cv             = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO(dpc, q)));
          DMO(dpc, q)    = sv;
          dcv            = _mm_add_ps(dcv, *tp); tp++;
        }
        if (!_mm_movemask_ps(cv)) break;
      }
    }

    /* Donor site q-loops. */
    if (don2 >= 0 || don1 >= 0 || don0 >= 0) {
      __m128 *dpp_lb = ox->dpf[i - min_intron - 3];
      __m128  lb_mpv, lb_dpv, lb_tsv;

      if (don2 >= 0) {
        for (nuc3 = 0; nuc3 < p7P_MAXNUC; nuc3++)
          C2[nuc3] = p7P_MINIDX(p7P_CODON3_FS1(r, s, nuc3), p7P_DEGEN1_C);
        C2[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(r, s, p7P_MAXCODONS1), p7P_DEGEN1_C);

        don_sig = (don2 == DONOR_GT) ? p7S_GTAG : (don2 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
        lb_mpv  = esl_sse_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
        lb_dpv  = esl_sse_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
        for (q = 0; q < Q; q++) {
          lb_tsv = _mm_max_ps(lb_mpv, lb_dpv);
          for (nuc3 = 0; nuc3 <= p7P_MAXNUC; nuc3++) {
            int slot = SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + don_sig;
            don_ovx[slot][q] = _mm_max_ps(don_ovx[slot][q],
                                           _mm_add_ps(lb_tsv, om_tr->rfv[C2[nuc3]][q]));
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
          lb_mpv = esl_sse_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
          lb_dpv = esl_sse_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
          for (q = 0; q < Q; q++) {
            lb_tsv = _mm_max_ps(lb_mpv, lb_dpv);
            don_ovx[slot][q] = _mm_max_ps(don_ovx[slot][q], lb_tsv);
            lb_mpv = MMO(dpp_lb, q);
            lb_dpv = DMO(dpp_lb, q);
          }
        }
      }

      if (don0 >= 0) {
        don_sig = (don0 == DONOR_GT) ? p7S_GTAG : (don0 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
        lb_mpv  = esl_sse_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
        lb_dpv  = esl_sse_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
        for (q = 0; q < Q; q++) {
          lb_tsv = _mm_max_ps(lb_mpv, lb_dpv);
          don_ovx[don_sig][q] = _mm_max_ps(don_ovx[don_sig][q], lb_tsv);
          lb_mpv = MMO(dpp_lb, q);
          lb_dpv = DMO(dpp_lb, q);
        }
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
      ox->xmx[i * p7X_NXCELLS + p7X_E] = xE;
      ox->xmx[i * p7X_NXCELLS + p7X_C] = ESL_MAX(
          ox->xmx[(i-3) * p7X_NXCELLS + p7X_C] + om_tr->xf[p7O_C][p7O_LOOP],
          xE + om_tr->xf[p7O_E][p7O_MOVE]);
    }
  } /* end main loop i = 3..L */

  if (global_end) {
    union { __m128 v; float p[4]; } um, ud;
    int qM = (om_tr->M - 1) % Q;
    int rM = (om_tr->M - 1) / Q;
    um.v = MMO(ox->dpf[L], qM);
    ud.v = DMO(ox->dpf[L], qM);
    xE   = ESL_MAX(um.p[rM], ud.p[rM]);
    ox->xmx[L * p7X_NXCELLS + p7X_E] = xE;
    ox->xmx[L * p7X_NXCELLS + p7X_C] = xE + om_tr->xf[p7O_E][p7O_MOVE];
  }

  return eslOK;
}


/* Function:  p7_Viterbi_SplicedTrace_sse()
 *
 * Purpose:   SSE implementation of spliced Viterbi traceback.
 *            See p7_Viterbi_SplicedTrace() documentation for details.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile is not 1-codon-length.
 *            <eslEFAIL>  if traceback fails.
 */
int
p7_Viterbi_SplicedTrace_sse(const ESL_DSQ *sub_dsq, const P7_OMX *ox,
                             const P7_FS_PROFILE *gm_tr, const float *signal_scores,
                             P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                             int min_intron, float *vitsc)
{
  float const *tsc    = gm_tr->tsc;
  float        r_tol  = 1e-5;
  float        a_tol  = 1e-4;
  int          M      = k_end - k_start + 1;
  int          L      = i_end - i_start + 1;
  int          i      = L;
  int          k      = 0;
  int          j;
  int          t, u, v, w, x;
  int          c, c0, c1, c2, c3;
  int          sprv, scur;
  int          sub_i, sub_k;
  int          donor_i;
  int          acc[3] = {0, 0, 0};
  int          don_sig;
  float        P_state;
  float        emit, emit0, emit1, emit2;
  float        vsc = 0.0f;
  int          status;

#define OMMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_M, (ii), (kk)))
#define ODMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_D, (ii), (kk)))
#define OIMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_I, (ii), (kk)))
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
          don_sig = -1;
          if      (SIGNAL(sub_dsq[sub_i-min_intron-j-1], sub_dsq[sub_i-min_intron-j]) == DONOR_GT) don_sig = p7S_GTAG;
          else if (SIGNAL(sub_dsq[sub_i-min_intron-j-1], sub_dsq[sub_i-min_intron-j]) == DONOR_GC) don_sig = p7S_GCAG;
          else if (SIGNAL(sub_dsq[sub_i-min_intron-j-1], sub_dsq[sub_i-min_intron-j]) == DONOR_AT) don_sig = p7S_ATAC;

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

#undef OMMo
#undef ODMo
#undef OIMo
#undef OXMXo

  return p7_trace_fs_Reverse(tr);
}

#undef TSC_P

#endif /* eslENABLE_SSE */
