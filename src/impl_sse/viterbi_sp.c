/* SSE-accelerated Spliced Viterbi algorithms; full matrix.
 *
 * Contents:
 *   1. p7_Viterbi_SplicedGlobal()
 *   2. p7_Viterbi_SplicedTrace_NoP()
 *   3. Benchmark driver.
 *   4. Unit tests.
 *   5. Test driver.
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

/* Log-space constant for the P->M transition (= logf(4.58e-5)). */
#define TSC_P_LOG  logf(4.58e-5f)

/*****************************************************************
 * 1. p7_Viterbi_SplicedGlobal()
 *****************************************************************/

/* Function:  p7_Viterbi_SplicedGlobal()
 * Synopsis:  SSE-accelerated fully global spliced Viterbi algorithm.
 *
 * Purpose:   SSE implementation of p7_GViterbi_SplicedGlobal_Dummy: a
 *            fully global translated Viterbi over a sub-region of the
 *            nucleotide sequence using a sub-region-striped codon
 *            profile.  The P state is always -eslINFINITY (Dummy
 *            behavior); donor/acceptor site tracking and OSPLICE_SCORES
 *            P_scores writes are reserved for the full implementation.
 *
 *            The profile <om_tr> must already be in log space and
 *            re-striped for the sub-region [k_start..k_end] via
 *            p7_fs_oprofile_SubConvert_Log().  om_tr->M is the
 *            sub-region length; i_start/i_end select the sub-sequence
 *            on <sub_dsq>.
 *
 *            The DP matrix <ox> must be allocated with
 *            p7_omx_Create_dpf(M, L, L, p7X_NSCELLS_SP) where
 *            L = i_end - i_start + 1.
 *
 *            On return, ox->xmx[L * p7X_NXCELLS + p7X_C] holds the
 *            global Viterbi score (log-odds, nats), and the full DP
 *            matrix is available in ox->dpf[0..L] for traceback.
 *
 * Args:      sub_dsq   - nucleotide subsequence, 1-based indexing,
 *                        positions 1..i_end (caller provides full dsq;
 *                        i_start is the first position of interest)
 *            om_tr     - log-space sub-region-striped FS oprofile
 *            ox        - full DP matrix, allocated for splice layout
 *            os        - vectorized splice-site score storage
 *            i_start   - first nucleotide position (1-based in sub_dsq)
 *            i_end     - last  nucleotide position (1-based in sub_dsq)
 *            min_intron - minimum intron length (reserved; not used in Dummy)
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile is not 1-codon-length or matrix
 *            nscells is wrong.
 */
int
p7_Viterbi_SplicedGlobal(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_OMX *ox, OSPLICE_SCORES *os, int i_start, int i_end, int min_intron)
{
  register __m128 mpv, dpv, ipv;  /* rightshifted M, D, I from dpp3, representing (i-3, k-1) */
  register __m128 sv;              /* pre-emission best-predecessor accumulator                */
  register __m128 msv;             /* M(i,k): post-emission match score                       */
  register __m128 dcv;             /* delayed D(i,q+1) carry for the D-propagation sweep      */
  __m128   infv;                   /* splatted -eslINFINITY                                   */
  __m128   b_entry;                /* B->M(k=1) global entry: non-inf only at q=0 of i=3      */
  float    xB;                     /* log B(0) = log T(N->B)                                  */
  float    xE;                     /* E-state scalar for final row                            */
  __m128  *dpc, *dpp3;             /* current and i-3 row pointers                            */
  __m128  *tp;                     /* rolling pointer through om_tr->tfv                      */
  int      Q  = p7O_NQF(om_tr->M); /* number of float SIMD stripes                           */
  int      L  = i_end - i_start + 1;
  int      C0;                     /* 3-nt codon index (1-codon-length profile)               */
  int      v, w, x;               /* nucleotide rolling window: v=i-2, w=i-1, x=i            */
  int      i, q, j, r;
  int      sub_i;

  if (om_tr->codon_lengths != 1)
    ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");
  if (ox->nscells != p7X_NSCELLS_SP)
    ESL_EXCEPTION(eslEINVAL, "DP matrix must have p7X_NSCELLS_SP cells per stripe position");

  ox->M              = om_tr->M;
  ox->L              = L;
  ox->has_own_scales = FALSE;
  ox->totscale       = 0.0f;
  infv               = _mm_set1_ps(-eslINFINITY);

  /* Reset OSPLICE_SCORES P_scores to -inf.
   * In the Dummy, P(i,k) is always -inf so P_scores is never updated. */
  for (j = 0; j < SIGNAL_MEM_SIZE; j++)
    for (q = 0; q < Q; q++)
      os->P_scores[j][q] = infv;

  /* Initialize all DP rows 0..L to -inf. */
  for (r = 0; r <= L; r++)
    for (q = 0; q < Q; q++)
      MMO_SP(ox->dpf[r], q) = DMO_SP(ox->dpf[r], q) =
      IMO_SP(ox->dpf[r], q) = PMO_SP(ox->dpf[r], q) = infv;

  /* Initialize special-state (xmx) rows to -inf; set B(0) = log T(N->B). */
  for (r = 0; r <= L; r++) {
    ox->xmx[r * p7X_NXCELLS + p7X_SCALE] = 0.0f;
    ox->xmx[r * p7X_NXCELLS + p7X_E]     = -eslINFINITY;
    ox->xmx[r * p7X_NXCELLS + p7X_N]     = -eslINFINITY;
    ox->xmx[r * p7X_NXCELLS + p7X_J]     = -eslINFINITY;
    ox->xmx[r * p7X_NXCELLS + p7X_B]     = -eslINFINITY;
    ox->xmx[r * p7X_NXCELLS + p7X_C]     = -eslINFINITY;
  }
  xB = om_tr->xf[p7O_N][p7O_MOVE];           /* log-space B(0): N->B transition */
  ox->xmx[0 * p7X_NXCELLS + p7X_B] = xB;

  /* Initialize nucleotide rolling window for rows 1..2.
   * After this setup: w = nuc at i_start, x = nuc at i_start+1.
   * (Rows 1 and 2 remain all -inf since no complete 3-nt codon is available.) */
  v = w = p7P_MAXCODONS1;   /* degenerate sentinel: -1 equivalent */
  x = (sub_dsq[i_start] < p7P_MAXNUC) ? sub_dsq[i_start] : p7P_MAXCODONS1;
  if (L >= 2) {
    w = x;
    x = (sub_dsq[i_start + 1] < p7P_MAXNUC) ? sub_dsq[i_start + 1] : p7P_MAXCODONS1;
  }

  /*----------------------------------------------------------------
   * Main DP recurrence: i = 3..L
   *
   * M(i,k) = max(M(i-3,k-1)*MM, I(i-3,k-1)*IM, D(i-3,k-1)*DM,
   *              B(i-3) [global entry: non-inf only at i=3, k=1])
   *          + emission(C0, k)
   * I(i,k) = max(M(i-3,k)*MI, I(i-3,k)*II)
   * D(i,k) = max(M(i,k-1)*MD, D(i,k-1)*DD)     [via delayed-D sweep]
   * P(i,k) = -inf                                [Dummy: always impossible]
   *
   * Rows 1 and 2 remain all -inf (no 3-nt codon available).
   *----------------------------------------------------------------*/
  for (i = 3; i <= L; i++) {
    sub_i = i_start + i - 1;
    v     = w;
    w     = x;
    x     = (sub_dsq[sub_i] < p7P_MAXNUC) ? sub_dsq[sub_i] : p7P_MAXCODONS1;

    C0 = p7P_CODON3_FS1(v, w, x);
    C0 = p7P_MINIDX(C0, p7P_DEGEN1_C);

    dpc  = ox->dpf[i];
    dpp3 = ox->dpf[i - 3];

    /* Rightshift dpp3's last stripe to seed the k-1 carry for the first stripe. */
    mpv = esl_sse_rightshift_ps(MMO_SP(dpp3, Q - 1), infv);
    dpv = esl_sse_rightshift_ps(DMO_SP(dpp3, Q - 1), infv);
    ipv = esl_sse_rightshift_ps(IMO_SP(dpp3, Q - 1), infv);

    tp  = om_tr->tfv;
    dcv = infv;

    /* Global B->M(k=1) entry: contributes only at i=3, lane 0 of q=0. */
    if (i == 3) {
      union { __m128 v; float p[4]; } u;
      u.v    = infv;
      u.p[0] = xB;   /* lane 0 = model position k=1 */
      b_entry = u.v;
    } else {
      b_entry = infv;
    }

    for (q = 0; q < Q; q++) {
      tp++;                                              /* skip BM transition slot (global entry handled via b_entry) */
      sv  =                _mm_add_ps(mpv, *tp); tp++;  /* MM: M(i-3, k-1) * MM */
      sv  = _mm_max_ps(sv, _mm_add_ps(ipv, *tp)); tp++; /* IM */
      sv  = _mm_max_ps(sv, _mm_add_ps(dpv, *tp)); tp++; /* DM */
      sv  = _mm_max_ps(sv, b_entry);                    /* B global entry (non-inf only at q=0, i=3) */
      b_entry = infv;                                   /* clear for subsequent stripes */

      /* Advance k-1 carry to current stripe (also used for I computation below). */
      mpv = MMO_SP(dpp3, q);
      dpv = DMO_SP(dpp3, q);
      ipv = IMO_SP(dpp3, q);

      msv = _mm_add_ps(sv, om_tr->rfv[C0][q]);          /* M(i,k) = best predecessor + emission */

      MMO_SP(dpc, q) = msv;
      DMO_SP(dpc, q) = dcv;                             /* delayed D from the previous stripe's M->D */
      dcv = _mm_add_ps(msv, *tp); tp++;                 /* MD: stage M->D carry for next stripe */

      /* I(i,k) = max(M(i-3,k)*MI, I(i-3,k)*II).  Uses same-k values from dpp3.
       * Stop codons have emission -inf; I state is biologically impossible there. */
      sv             = _mm_add_ps(mpv, *tp); tp++;      /* MI */
      sv             = _mm_max_ps(sv, _mm_add_ps(ipv, *tp)); tp++;  /* II */
      { /* mask out I at stop-codon positions (emission == -inf) */
        __m128 stop_mask = _mm_cmpeq_ps(om_tr->rfv[C0][q], infv);
        sv = _mm_or_ps(_mm_andnot_ps(stop_mask, sv), _mm_and_ps(stop_mask, infv));
      }
      IMO_SP(dpc, q) = sv;

      PMO_SP(dpc, q) = infv;                            /* Dummy: P always -inf */
    }

    /*------------------------------------------------------------
     * DD propagation sweep.
     * The delayed-D scheme requires a right-shifted pass after the
     * main loop to propagate D(i,k) = max(M(i,k-1)*MD, D(i,k-1)*DD)
     * along the full model.  Three extra passes handle the longest
     * DD chain of length 4Q that can remain after a single sweep.
     *------------------------------------------------------------*/
    dcv          = esl_sse_rightshift_ps(dcv, infv);
    DMO_SP(dpc, 0) = infv;          /* D(i,k=1) boundary: no predecessor at k=0 */
    tp = om_tr->tfv + 7 * Q;        /* DD transitions follow the 7 main transitions */
    for (q = 0; q < Q; q++) {
      DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
      dcv = _mm_add_ps(DMO_SP(dpc, q), *tp); tp++;
    }
    if (om_tr->M < 100) {
      /* Small model: fixed 3 extra passes (no early exit needed). */
      for (j = 1; j < 4; j++) {
        dcv = esl_sse_rightshift_ps(dcv, infv);
        tp  = om_tr->tfv + 7 * Q;
        for (q = 0; q < Q; q++) {
          DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
          dcv = _mm_add_ps(dcv, *tp); tp++;
        }
      }
    } else {
      /* Larger model: early-exit when no D value improved. */
      for (j = 1; j < 4; j++) {
        register __m128 cv;
        dcv = esl_sse_rightshift_ps(dcv, infv);
        tp  = om_tr->tfv + 7 * Q;
        cv  = _mm_setzero_ps();
        for (q = 0; q < Q; q++) {
          sv             = _mm_max_ps(dcv, DMO_SP(dpc, q));
          cv             = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO_SP(dpc, q)));
          DMO_SP(dpc, q) = sv;
          dcv            = _mm_add_ps(dcv, *tp); tp++;
        }
        if (!_mm_movemask_ps(cv)) break;
      }
    }
  } /* end main loop i = 3..L */

  /*----------------------------------------------------------------
   * Final score: global model exit at (L, M).
   * E(L) = max(M(L, M), D(L, M)).
   * C(L) = E(L) + xf[E][MOVE].
   *
   * In the global model, exit is forced at the last model position k=M.
   * We extract M(L,M) and D(L,M) directly from the DP matrix.
   *----------------------------------------------------------------*/
  {
    union { __m128 v; float p[4]; } um, ud;
    int qM = (om_tr->M - 1) % Q;   /* stripe containing k=M */
    int rM = (om_tr->M - 1) / Q;   /* lane within that stripe */
    um.v = MMO_SP(ox->dpf[L], qM);
    ud.v = DMO_SP(ox->dpf[L], qM);
    xE   = ESL_MAX(um.p[rM], ud.p[rM]);
  }
  ox->xmx[L * p7X_NXCELLS + p7X_E] = xE;
  ox->xmx[L * p7X_NXCELLS + p7X_C] = xE + om_tr->xf[p7O_E][p7O_MOVE];

  return eslOK;
}
/*------------------ end, p7_Viterbi_SplicedGlobal() ---------------*/


/*****************************************************************
 * 2b. p7_Viterbi_SplicedGlobal_NoP()
 *
 * SSE counterpart of p7_GViterbi_SplicedGlobal_DummyNoP.
 * Uses the standard 3-cell (M,D,I) layout with a 3-slot circular
 * P-state buffer pvx[slot*Q + q] (slot = i%3) instead of a P column
 * in the main DP matrix.
 *
 * Mirrors p7_GViterbi_SplicedGlobal_DummyNoP fully:
 *   - Acceptor-site q-loop writes to pvx using os->P_scores.
 *   - Donor-site q-loops (after the DD sweep) write to os->P_scores
 *     using the DP values at row i-min_intron-3.
 *
 * The pvx buffer and don tracking are allocated and managed internally.
 *****************************************************************/

int
p7_Viterbi_SplicedGlobal_NoP(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr,
                              P7_OMX *ox, OSPLICE_SCORES *os,
                              int i_start, int i_end, int min_intron)
{
  register __m128 mpv, dpv, ipv;
  register __m128 sv, msv, dcv;
  __m128   infv, b_entry;
  float    xB, xE;
  __m128  *dpc, *dpp3;
  __m128  *tp;
  __m128  *pvx;                    /* 4-slot circular P buffer: pvx[slot*Q + q] */
  __m128   ppv;                    /* P(i, k-1) carry for current row            */
  __m128   tsc_p_vec;              /* broadcast TSC_P_LOG                        */
  int      Q  = p7O_NQF(om_tr->M);
  int      L  = i_end - i_start + 1;
  int      C0;
  int      C1[p7P_MAXNUC + 1];    /* split-codon indices for acc1 acceptor case  */
  int      C2[p7P_MAXNUC + 1];    /* split-codon indices for don2 donor case     */
  int      v, w, x;
  int      r, s, t, u;            /* nucleotide lookback window for donor sites  */
  int      nuc1, nuc3;
  int      acc0, acc1, acc2;      /* acceptor site type for rows i, i-1, i-2    */
  int      don0, don1, don2;      /* donor site type for rows i, i-1, i-2       */
  int      don_sig;               /* p7S_GTAG / p7S_GCAG / p7S_ATAC             */
  int      pv_i;                  /* i % 4: current-row slot index into pvx     */
  int      pv_pi;                 /* (i-3) % 4: i-3 row slot index into pvx    */
  int      i, q, j, ri, sub_i;
  int      status;

  if (om_tr->codon_lengths != 1)
    ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");
  if (ox->nscells != p7X_NSCELLS)
    ESL_EXCEPTION(eslEINVAL, "DP matrix must have p7X_NSCELLS cells per stripe position");

  ox->M              = om_tr->M;
  ox->L              = L;
  ox->has_own_scales = FALSE;
  ox->totscale       = 0.0f;
  infv               = _mm_set1_ps(-eslINFINITY);
  tsc_p_vec          = _mm_set1_ps(TSC_P_LOG);

  pvx = NULL;
  ESL_ALLOC(pvx, sizeof(__m128) * 4 * Q);
  for (ri = 0; ri < 4 * Q; ri++) pvx[ri] = infv;

  /* Reset os->P_scores to -inf; donor writes will accumulate into it. */
  for (j = 0; j < SIGNAL_MEM_SIZE; j++)
    for (q = 0; q < Q; q++)
      os->P_scores[j][q] = infv;

  /* Initialize all DP rows 0..L to -inf. */
  for (ri = 0; ri <= L; ri++)
    for (q = 0; q < Q; q++)
      MMO(ox->dpf[ri], q) = DMO(ox->dpf[ri], q) = IMO(ox->dpf[ri], q) = infv;

  /* Initialize special-state rows; set B(0) = log T(N->B). */
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

  /* Nucleotide rolling window and acceptor/donor site tracker init. */
  acc0 = acc1 = acc2 = -1;
  don0 = don1 = don2 = -1;
  v = w = p7P_MAXCODONS1;
  x = (sub_dsq[i_start] < p7P_MAXNUC) ? sub_dsq[i_start] : p7P_MAXCODONS1;
  if (L >= 2) {
    w = x;
    x = (sub_dsq[i_start + 1] < p7P_MAXNUC) ? sub_dsq[i_start + 1] : p7P_MAXCODONS1;
  }

  /* Donor lookback window init: s/t/u prime the rolling window so that
   * at the first main iteration (i = min_intron+3) the shift r=s; s=t; t=u
   * yields r=sub_dsq[i_start], s=sub_dsq[i_start+1], t=sub_dsq[i_start+2]. */
  r = p7P_MAXCODONS1;
  s = (sub_dsq[i_start]     < p7P_MAXNUC) ? sub_dsq[i_start]     : p7P_MAXCODONS1;
  t = (sub_dsq[i_start + 1] < p7P_MAXNUC) ? sub_dsq[i_start + 1] : p7P_MAXCODONS1;
  u = (sub_dsq[i_start + 2] < p7P_MAXNUC) ? sub_dsq[i_start + 2] : p7P_MAXCODONS1;

  for (i = 3; i <= L; i++) {
    sub_i = i_start + i - 1;

    /* Donor lookback window: only valid when i >= min_intron+3 (avoids
     * accessing sub_dsq before i_start).  Guards don update too. */
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

    /* Split-codon indices for acc1 acceptor case; nuc3 for acc2 SSX2 lookups */
    for (nuc1 = 0; nuc1 < p7P_MAXNUC; nuc1++)
      C1[nuc1] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, w, x), p7P_DEGEN1_C);
    C1[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(p7P_MAXCODONS1, w, x), p7P_DEGEN1_C);
    nuc3 = ESL_MIN(x, p7P_MAXNUC);

    /* Shift acceptor site window: acc0 fires 2 rows ago (SSX0), acc1 one row ago (SSX1),
     * acc2 is this row's acceptor dinucleotide (SSX2). */
    acc0 = acc1;
    acc1 = acc2;
    if      (SIGNAL(v, w) == ACCEPT_AG) acc2 = ACCEPT_AG;
    else if (SIGNAL(v, w) == ACCEPT_AC) acc2 = ACCEPT_AC;
    else                                acc2 = -1;

    /* Shift donor site window (guarded): don2 is this row's donor dinucleotide,
     * don1 one row ago (1 pre-donor nuc = r), don0 two rows ago (0 pre-donor nucs). */
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

    mpv = esl_sse_rightshift_ps(MMO(dpp3, Q - 1), infv);
    dpv = esl_sse_rightshift_ps(DMO(dpp3, Q - 1), infv);
    ipv = esl_sse_rightshift_ps(IMO(dpp3, Q - 1), infv);
    /* Pre-loop: compute P(i, k) for all stripes and write to pvx before
     * the main MM/IM/DM/PM loop reads them.  This ensures the PM transition
     * uses P(i, k-1) from the current row rather than the stale P(i-3, k-1)
     * that would otherwise sit in the circular buffer slot pv_i. */
    for (q = 0; q < Q; q++) {
      __m128 psv = infv;
      if (acc0 == ACCEPT_AG) {
        psv = _mm_max_ps(psv, _mm_add_ps(
            _mm_max_ps(_mm_add_ps(os->P_scores[p7S_GTAG][q], os->signal_scores[p7S_GTAG]),
                       _mm_add_ps(os->P_scores[p7S_GCAG][q], os->signal_scores[p7S_GCAG])),
            om_tr->rfv[C0][q]));
      } else if (acc0 == ACCEPT_AC) {
        psv = _mm_max_ps(psv, _mm_add_ps(
            _mm_add_ps(os->P_scores[p7S_ATAC][q], os->signal_scores[p7S_ATAC]),
            om_tr->rfv[C0][q]));
      }
      if (acc1 == ACCEPT_AG) {
        for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
          psv = _mm_max_ps(psv, _mm_add_ps(
              _mm_max_ps(_mm_add_ps(os->P_scores[SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + p7S_GTAG][q], os->signal_scores[p7S_GTAG]),
                         _mm_add_ps(os->P_scores[SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + p7S_GCAG][q], os->signal_scores[p7S_GCAG])),
              om_tr->rfv[C1[nuc1]][q]));
        }
      } else if (acc1 == ACCEPT_AC) {
        for (nuc1 = 0; nuc1 <= p7P_MAXNUC; nuc1++) {
          psv = _mm_max_ps(psv, _mm_add_ps(
              _mm_add_ps(os->P_scores[SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + p7S_ATAC][q], os->signal_scores[p7S_ATAC]),
              om_tr->rfv[C1[nuc1]][q]));
        }
      }
      if (acc2 == ACCEPT_AG) {
        psv = _mm_max_ps(psv,
            _mm_max_ps(_mm_add_ps(os->P_scores[SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + p7S_GTAG][q], os->signal_scores[p7S_GTAG]),
                       _mm_add_ps(os->P_scores[SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + p7S_GCAG][q], os->signal_scores[p7S_GCAG])));
      } else if (acc2 == ACCEPT_AC) {
        psv = _mm_max_ps(psv,
            _mm_add_ps(os->P_scores[SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + p7S_ATAC][q], os->signal_scores[p7S_ATAC]));
      }
      pvx[pv_i * Q + q] = psv;
    }

    /* Seed ppv from P(i-3, *) in pvx; PM transition uses P(i-3, k-1). */
    ppv = esl_sse_rightshift_ps(pvx[pv_pi * Q + Q - 1], infv);  /* P(i-3, k-1) seed */

    tp  = om_tr->tfv;
    dcv = infv;

    if (i == 3) {
      union { __m128 v; float p[4]; } u;
      u.v    = infv;
      u.p[0] = xB;
      b_entry = u.v;
    } else {
      b_entry = infv;
    }

    for (q = 0; q < Q; q++) {
      __m128 ppv_next = pvx[pv_pi * Q + q];   /* P(i-3, k); carry as P(i-3, k-1) to next stripe */

      tp++;                                              /* skip BM */
      sv  =                _mm_add_ps(mpv, *tp); tp++;  /* MM */
      sv  = _mm_max_ps(sv, _mm_add_ps(ipv, *tp)); tp++; /* IM */
      sv  = _mm_max_ps(sv, _mm_add_ps(dpv, *tp)); tp++; /* DM */
      sv  = _mm_max_ps(sv, _mm_add_ps(ppv, tsc_p_vec)); /* PM */
      sv  = _mm_max_ps(sv, b_entry);
      b_entry = infv;

      mpv = MMO(dpp3, q);
      dpv = DMO(dpp3, q);
      ipv = IMO(dpp3, q);

      msv = _mm_add_ps(sv, om_tr->rfv[C0][q]);
      MMO(dpc, q) = msv;
      DMO(dpc, q) = dcv;
      dcv = _mm_add_ps(msv, *tp); tp++;                 /* MD */

      sv  = _mm_add_ps(mpv, *tp); tp++;                 /* MI */
      sv  = _mm_max_ps(sv, _mm_add_ps(ipv, *tp)); tp++; /* II */
      { /* stop-codon mask: I = -inf where emission == -inf */
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

    /*------------------------------------------------------------
     * Donor site q-loops.
     * Run after the DD sweep so MMO/DMO(dpc) are final.
     * Read DP values at row i-min_intron-3 (the last complete exon
     * codon before the donor site), accumulate into os->P_scores.
     *
     * don2: 2 pre-donor nucs (r, s); C2 codon emission baked in.
     * don1: 1 pre-donor nuc  (r);    emission added at acceptor time.
     * don0: 0 pre-donor nucs;        emission added at acceptor time.
     *------------------------------------------------------------*/
    if (don2 >= 0 || don1 >= 0 || don0 >= 0) {
      __m128 *dpp_lb = ox->dpf[i - min_intron - 3];
      __m128  lb_mpv, lb_dpv, lb_tsv;

      if (don2 >= 0) {
        /* C2[nuc3] = codon(r, s, nuc3) for all post-acceptor nuc3 values.
         * Emission is baked in; nuc3 is reused as loop variable (safe: acceptor
         * block finished using it as ESL_MIN(x, p7P_MAXNUC)). */
        for (nuc3 = 0; nuc3 < p7P_MAXNUC; nuc3++)
          C2[nuc3] = p7P_MINIDX(p7P_CODON3_FS1(r, s, nuc3), p7P_DEGEN1_C);
        C2[p7P_MAXNUC] = p7P_MINIDX(p7P_CODON3_FS1(r, s, p7P_MAXCODONS1), p7P_DEGEN1_C);

        don_sig = (don2 == DONOR_GT) ? p7S_GTAG : (don2 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
        lb_mpv  = esl_sse_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
        lb_dpv  = esl_sse_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
        for (q = 0; q < Q; q++) {
          lb_tsv = _mm_max_ps(lb_mpv, lb_dpv);           /* TMP_SC: best M/D at k-1 */
          for (nuc3 = 0; nuc3 <= p7P_MAXNUC; nuc3++) {
            int slot = SPLICE_OFFSET_2 + nuc3 * p7S_SPLICE_SIGNALS + don_sig;
            os->P_scores[slot][q] = _mm_max_ps(os->P_scores[slot][q],
                                               _mm_add_ps(lb_tsv, om_tr->rfv[C2[nuc3]][q]));
          }
          lb_mpv = MMO(dpp_lb, q);
          lb_dpv = DMO(dpp_lb, q);
        }
      }

      if (don1 >= 0) {
        /* nuc1 = pre-donor nucleotide (r); emission added at acceptor time. */
        nuc1    = ESL_MIN(r, p7P_MAXNUC);
        don_sig = (don1 == DONOR_GT) ? p7S_GTAG : (don1 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
        {
          int slot = SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + don_sig;
          lb_mpv = esl_sse_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
          lb_dpv = esl_sse_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
          for (q = 0; q < Q; q++) {
            lb_tsv = _mm_max_ps(lb_mpv, lb_dpv);
            os->P_scores[slot][q] = _mm_max_ps(os->P_scores[slot][q], lb_tsv);
            lb_mpv = MMO(dpp_lb, q);
            lb_dpv = DMO(dpp_lb, q);
          }
        }
      }

      if (don0 >= 0) {
        /* 0 pre-donor nucs; emission added at acceptor time. */
        don_sig = (don0 == DONOR_GT) ? p7S_GTAG : (don0 == DONOR_GC) ? p7S_GCAG : p7S_ATAC;
        lb_mpv  = esl_sse_rightshift_ps(MMO(dpp_lb, Q - 1), infv);
        lb_dpv  = esl_sse_rightshift_ps(DMO(dpp_lb, Q - 1), infv);
        for (q = 0; q < Q; q++) {
          lb_tsv = _mm_max_ps(lb_mpv, lb_dpv);
          os->P_scores[don_sig][q] = _mm_max_ps(os->P_scores[don_sig][q], lb_tsv);
          lb_mpv = MMO(dpp_lb, q);
          lb_dpv = DMO(dpp_lb, q);
        }
      }
    }
  } /* end main loop i = 3..L */

  {
    union { __m128 v; float p[4]; } um, ud;
    int qM = (om_tr->M - 1) % Q;
    int rM = (om_tr->M - 1) / Q;
    um.v = MMO(ox->dpf[L], qM);
    ud.v = DMO(ox->dpf[L], qM);
    xE   = ESL_MAX(um.p[rM], ud.p[rM]);
  }
  ox->xmx[L * p7X_NXCELLS + p7X_E] = xE;
  ox->xmx[L * p7X_NXCELLS + p7X_C] = xE + om_tr->xf[p7O_E][p7O_MOVE];

  free(pvx);
  return eslOK;

 ERROR:
  if (pvx != NULL) free(pvx);
  return status;
}
/*--------------- end, p7_Viterbi_SplicedGlobal_NoP() -----------*/



/*****************************************************************
 * 2. p7_Viterbi_SplicedTrace_NoP()
 *****************************************************************/

/* Function:  p7_Viterbi_SplicedTrace_NoP()
 * Synopsis:  Traceback through a p7_Viterbi_SplicedGlobal_NoP() DP matrix.
 *
 * Purpose:   Given a filled DP matrix <ox> from p7_Viterbi_SplicedGlobal_NoP(),
 *            trace the optimal path and record it in <tr>.
 *
 *            The P state is not stored in the DP matrix and must be
 *            reconstructed by scanning donor and acceptor dinucleotides in
 *            <sub_dsq>, mirroring the scalar p7_GViterbi_SplicedTrace_NoP().
 *            Emission and transition scores for the reconstruction are read
 *            from the scalar profile <gm_tr>.
 *
 * Args:      sub_dsq      - nucleotide subsequence, 1-based
 *            ox           - filled DP matrix from p7_Viterbi_SplicedGlobal_NoP()
 *            gm_tr        - scalar FS profile (transition/emission scores)
 *            signal_scores - log-probabilities for GT-AG, GC-AG, AT-AC signals
 *            tr           - allocated, empty trace to fill
 *            i_start      - first nucleotide position (1-based in sub_dsq)
 *            i_end        - last  nucleotide position
 *            k_start      - first model position (1-based in gm_tr)
 *            k_end        - last  model position
 *            min_intron   - minimum intron length
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile is not 1-codon-length.
 *            <eslEFAIL>  if traceback fails to identify a predecessor state.
 */
int
p7_Viterbi_SplicedTrace_NoP(const ESL_DSQ *sub_dsq, const P7_OMX *ox,
                             const P7_FS_PROFILE *gm_tr, const float *signal_scores,
                             P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end, int min_intron)
{
  float const *tsc    = gm_tr->tsc;          /* so TSC() macro works                        */
  float        tol    = 1e-5;
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
  int          status;

  /* Safe scalar accessors for the striped SSE DP matrix.
   * Guard k<1 to avoid negative-modulo undefined behaviour at model boundaries. */
#define OMMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_M, (ii), (kk)))
#define ODMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_D, (ii), (kk)))
#define OIMo(ii,kk)  ((kk) < 1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_I, (ii), (kk)))
#define OXMXo(ii,s)  (ox->xmx[(ii) * p7X_NXCELLS + (s)])

  if (gm_tr->codon_lengths != 1) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");
#if eslDEBUGLEVEL > 0
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i+i_start-1, 0)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i+i_start-1, 0)) != eslOK) return status;

  sprv = p7T_C;

  while (sprv != p7T_S) {
    switch (sprv) {

    case p7T_C:
      if      (OXMXo(i, p7X_C) == -eslINFINITY)
        ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
      else if (esl_FCompare_old(OXMXo(i, p7X_C), OXMXo(i-3, p7X_C) + gm_tr->xsc[p7P_C][p7P_LOOP], tol) == eslOK)
        scur = p7T_C;
      else if (esl_FCompare_old(OXMXo(i, p7X_C), OXMXo(i,   p7X_E) + gm_tr->xsc[p7P_E][p7P_MOVE], tol) == eslOK)
        scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7T_E:
      if (OXMXo(i, p7X_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
      for (k = M; k >= 1; k--) {
        if (esl_FCompare_old(OXMXo(i, p7X_E), OMMo(i,k), tol) == eslOK) { scur = p7T_M; break; }
        if (esl_FCompare_old(OXMXo(i, p7X_E), ODMo(i,k), tol) == eslOK) { scur = p7T_D; break; }
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
      if      (esl_FCompare_old(OMMo(i,k), OMMo(i-3,k-1) + TSC(p7P_MM,sub_k-1) + emit, tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(OMMo(i,k), OIMo(i-3,k-1) + TSC(p7P_IM,sub_k-1) + emit, tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare_old(OMMo(i,k), ODMo(i-3,k-1) + TSC(p7P_DM,sub_k-1) + emit, tol) == eslOK) scur = p7T_D;
      else if (esl_FCompare_old(OMMo(i,k), OXMXo(i-3, p7X_B)                   + emit, tol) == eslOK) scur = p7T_B;
      else {
        /* P state is not stored; reconstruct by scanning donor/acceptor sites. */
        if (i < min_intron + 7) ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k,i);

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

        for (j = 0; j < i - min_intron - 4 && scur == -1; j++) {
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
              if (scur == -1 && acc[2] == 1) {
                P_state = ESL_MAX(OMMo(i-min_intron-j-4,k-2), ODMo(i-min_intron-j-4,k-2)) + signal_scores[don_sig] + emit2;
                if (esl_FCompare_old(OMMo(i,k), P_state + TSC_P_LOG + emit, tol) == eslOK) { scur = p7T_P; c = 2; donor_i = i-min_intron-j-4; }
              }
              if (scur == -1 && acc[1] == 1) {
                P_state = ESL_MAX(OMMo(i-min_intron-j-3,k-2), ODMo(i-min_intron-j-3,k-2)) + signal_scores[don_sig] + emit1;
                if (esl_FCompare_old(OMMo(i,k), P_state + TSC_P_LOG + emit, tol) == eslOK) { scur = p7T_P; c = 1; donor_i = i-min_intron-j-3; }
              }
              if (scur == -1 && acc[0] == 1) {
                P_state = ESL_MAX(OMMo(i-min_intron-j-2,k-2), ODMo(i-min_intron-j-2,k-2)) + signal_scores[don_sig] + emit0;
                if (esl_FCompare_old(OMMo(i,k), P_state + TSC_P_LOG + emit, tol) == eslOK) { scur = p7T_P; c = 0; donor_i = i-min_intron-j-2; }
              }
            } else { /* p7S_ATAC */
              if (scur == -1 && acc[2] == 2) {
                P_state = ESL_MAX(OMMo(i-min_intron-j-4,k-2), ODMo(i-min_intron-j-4,k-2)) + signal_scores[don_sig] + emit2;
                if (esl_FCompare_old(OMMo(i,k), P_state + TSC_P_LOG + emit, tol) == eslOK) { scur = p7T_P; c = 2; donor_i = i-min_intron-j-4; }
              }
              if (scur == -1 && acc[1] == 2) {
                P_state = ESL_MAX(OMMo(i-min_intron-j-3,k-2), ODMo(i-min_intron-j-3,k-2)) + signal_scores[don_sig] + emit1;
                if (esl_FCompare_old(OMMo(i,k), P_state + TSC_P_LOG + emit, tol) == eslOK) { scur = p7T_P; c = 1; donor_i = i-min_intron-j-3; }
              }
              if (scur == -1 && acc[0] == 2) {
                P_state = ESL_MAX(OMMo(i-min_intron-j-2,k-2), ODMo(i-min_intron-j-2,k-2)) + signal_scores[don_sig] + emit0;
                if (esl_FCompare_old(OMMo(i,k), P_state + TSC_P_LOG + emit, tol) == eslOK) { scur = p7T_P; c = 0; donor_i = i-min_intron-j-2; }
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
      if      (esl_FCompare_old(ODMo(i,k), OMMo(i,k-1) + TSC(p7P_MD,sub_k-1), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(ODMo(i,k), ODMo(i,k-1) + TSC(p7P_DD,sub_k-1), tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      k--;
      break;

    case p7T_I:
      if (OIMo(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);
      sub_k = k_start + k - 1;
      if      (esl_FCompare_old(OIMo(i,k), OMMo(i-3,k) + TSC(p7P_MI,sub_k), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(OIMo(i,k), OIMo(i-3,k) + TSC(p7P_II,sub_k), tol) == eslOK) scur = p7T_I;
      else ESL_EXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      i -= 3;
      break;

    case p7T_P:
      scur = (OMMo(donor_i,k-1) >= ODMo(donor_i,k-1)) ? p7T_M : p7T_D;
      k--; i = donor_i;
      break;

    case p7T_N:
      /* N is not stored in the SSE fill; the only path is S->N(0)->B(0). */
      scur = (i == 0) ? p7T_S : p7T_N;
      break;

    case p7T_B:
      /* B is only at i=0 and always came from N->B in this global model. */
      scur = p7T_N;
      break;

    default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch */

    if      (scur == p7T_M) c = 3;
    else if (scur != p7T_P) c = 0;

    if ((status = p7_trace_fs_Append(tr, scur, k_start+k-1, i_start+i-1, c)) != eslOK) return status;

    /* For N and C: i decrement is deferred to handle self-loops correctly. */
    if ((scur == p7T_N || scur == p7T_C) && scur == sprv) i--;

    sprv = scur;
  } /* end traceback, at S state */

  tr->M = M;
  tr->L = L;

#undef OMMo
#undef ODMo
#undef OIMo
#undef OXMXo

  return p7_trace_fs_Reverse(tr);
}
/*------------- end, p7_Viterbi_SplicedTrace_NoP() --------------*/



/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/
#ifdef p7VITERBI_SP_BENCHMARK
/*
 gcc -g -O3 -Wall -std=gnu99 -o viterbi_sp_benchmark -I. -I.. -I../../easel -L. -L.. -L../../easel -Dp7VITERBI_SP_BENCHMARK viterbi_sp.c -lhmmer -leasel -lm
 
   ./viterbi_sp_benchmark <hmmfile>
 */

#include "p7_config.h"

#include <math.h>
#include <stdio.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_sse.h"

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
static char benchmark_banner[] = "benchmark driver for spliced Viterbi algorithms";

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
  P7_FS_OPROFILE *om_tr        = NULL;
  P7_OMX         *ox           = NULL;
  P7_TRACE       *tr           = NULL;
  ESL_GENCODE    *gcode        = NULL;
  P7_CODONTABLE  *codon_table  = NULL;
  ESL_SQ         *sq           = NULL;
  SPLICE_PIPELINE *pli         = NULL;
  OSPLICE_SCORES *oss          = NULL;
  int             do_T         = esl_opt_GetBoolean(go, "-T"); 
  int             N            = esl_opt_GetInteger(go, "-N");
  int             I            = esl_opt_GetInteger(go, "-I");
  int             intron_total = I + 4;
  ESL_DSQ        *dsq          = NULL;
  int             i, j, k, L_amino, L_dna_total;
  int64_t         total_cells;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)          != eslOK) p7_Fail("Failed to read HMM");

  gcode       = esl_gencode_Create(abcDNA, abcAA);
  bgAA        = p7_bg_Create(abcAA);
  gm          = p7_profile_Create(hmm->M, abcAA);
  gm_tr       = p7_profile_fs_Create(hmm->M, abcAA, 1);
  om_tr       = p7_fs_oprofile_Create(hmm->M, abcAA, 1); 
  codon_table = p7_codontable_Create(gcode);
  sq          = esl_sq_CreateDigital(abcAA);

  p7_ProfileConfig   (hmm, bgAA,        gm,    hmm->M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, hmm->M, p7_UNILOCAL);
  p7_fs_oprofile_Convert_Log(gm_tr, om_tr);

  pli = p7_splicepipeline_Create(NULL, hmm->M, hmm->M * 3);
  ox = p7_omx_Create_dpf(hmm->M, hmm->M, hmm->M, p7X_NSCELLS); 
  oss = p7_osplicescores_Create(hmm->M);

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

      p7_fs_ReconfigLength(gm_tr, L_dna_total/3);
      p7_fs_oprofile_ReconfigLength_Log(om_tr, L_dna_total/3);
      p7_omx_GrowTo_dpf(ox, hmm->M, L_dna_total, L_dna_total);
      p7_Viterbi_SplicedGlobal_NoP(dsq, om_tr, ox, oss, 1, L_dna_total, 13);

      if(do_T)
        p7_Viterbi_SplicedTrace_NoP(dsq, ox, gm_tr, pli->splice_scores->signal_scores, tr, 1, L_dna_total, 1, hmm->M, pli->min_intron); 
      
      p7_trace_Reuse(tr); 
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
  esl_sq_Destroy(sq);
  p7_codontable_Destroy(codon_table);
  p7_profile_fs_Destroy(gm_tr);
  p7_fs_oprofile_Destroy(om_tr);
  p7_omx_Destroy(ox);
  p7_splicepipeline_Destroy(pli);
  p7_osplicescores_Destroy(oss);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bgAA);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  p7_trace_fs_Destroy(tr);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}

#endif /*p7GENERIC_VITERBI_SPLICED_BENCHMARK*/
/*----------------- end, benchmark driver -----------------------*/



/*****************************************************************
 * 4. Unit tests.
 *****************************************************************/
#ifdef p7VITERBI_SP_TESTDRIVE

#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"

static void
utest_viterbi_sp(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
              ESL_GENCODE *gcode, P7_BG *bgAA, P7_CODONTABLE *codon_table,
              int M, int N, int intron_len)
{
  char           *msg         = "spliced viterbi unit test failed";
  P7_HMM         *hmm         = NULL;
  P7_FS_PROFILE  *gm_tr       = p7_profile_fs_Create(M, abcAA, 1);
  P7_FS_OPROFILE *om_tr       = p7_fs_oprofile_Create(M, abcAA, 1);
  ESL_SQ         *sq          = esl_sq_CreateDigital(abcAA);
  P7_OMX         *ox          = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS_SP);
  P7_OMX         *ox_NoP      = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS);
  P7_GMX         *gx_NoP      = p7_gmx_Create(M, M, M, p7G_NSCELLS);
  P7_IVX         *iv_nop      = p7_ivx_Create(M, SPLICE_ROWS);
  P7_TRACE       *otr         = p7_trace_fs_Create();
  P7_TRACE       *gtr         = p7_trace_fs_Create();
  ESL_DSQ        *dsq         = NULL;
  SPLICE_PIPELINE *pli        = NULL;
  OSPLICE_SCORES  *oss        = NULL;
  int             intron_total = intron_len + 4;  /* GT + intron_len random nucs + AG */
  int             L_amino, L_dna_total;
  int             i, j;
  int             k_start, k_end;
  int             sub_M;
  float           final_gC, final_oC;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, M, p7_UNILOCAL);
  p7_fs_oprofile_Convert_Log(gm_tr, om_tr);

  pli = p7_splicepipeline_Create(NULL, M, M * 3);
  oss = p7_osplicescores_Create(M); 

  p7_emit_SimpleConsensus(hmm, sq);
  L_amino     = sq->n;

  while (N--)
    {
	  k_start = k_end = 0;
	  while(k_end <= k_start || k_end > L_amino) {
	    k_start = esl_rnd_Roll(r, L_amino/4) + 1;
	    k_end   = esl_rnd_Roll(r, L_amino/2) + k_start + L_amino/4;
	  }
      sub_M = k_end - k_start + 1;
      
      L_dna_total = sub_M * 3 + intron_total;
      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna_total + 2))) == NULL) esl_fatal("malloc failed");
      dsq[0] = dsq[L_dna_total + 1] = eslDSQ_SENTINEL;

      /* Reverse-translate first half of the sequence (exon 1) */
      j = 1;
      for (i = k_start; i <= k_start+sub_M/2; i++) {
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
      for (i = k_start+ sub_M/2 + 1; i <= k_end; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }
      
      p7_fs_ReconfigLength(gm_tr, L_dna_total/3);
      p7_gmx_GrowTo(pli->vit, sub_M, L_dna_total, L_dna_total);
      p7_splicescores_GrowTo(pli->splice_scores, sub_M);

      /* --- Score comparison: generic Dummy and SSE Dummy must agree --- */
//      p7_GViterbi_SplicedGlobal_Dummy(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, k_start, k_end, pli->min_intron);
//      final_gC = pli->vit->xmx[L_dna_total * p7G_NXCELLS + p7G_C];

	  p7_fs_oprofile_SubConvert_Log(gm_tr, om_tr, k_start, k_end);
      p7_fs_oprofile_ReconfigLength_Log(om_tr, L_dna_total/3);
	  p7_omx_GrowTo_dpf(ox, sub_M, L_dna_total, L_dna_total);
	  p7_osplicescores_GrowTo(oss, sub_M);

//      p7_Viterbi_SplicedGlobal(dsq, om_tr, ox, oss, 1, L_dna_total , pli->min_intron);
//      final_oC = ox->xmx[L_dna_total * p7X_NXCELLS + p7X_C];

	  /* Scores must agree */
//      if (fabs(final_gC - final_oC) > 0.001)
//        esl_fatal("%s: generic %.4f != SSE %.4f", msg, final_gC, final_oC);

      p7_gmx_GrowTo(gx_NoP, sub_M, L_dna_total, L_dna_total);
      p7_ivx_GrowTo(iv_nop, sub_M, SPLICE_ROWS);
      p7_GViterbi_SplicedGlobal_NoP(dsq, gm_tr, gx_NoP, iv_nop, pli->splice_scores, 1, L_dna_total, k_start, k_end, pli->min_intron);
      final_gC = gx_NoP->xmx[L_dna_total * p7G_NXCELLS + p7G_C];

      p7_omx_GrowTo_dpf(ox_NoP, sub_M, L_dna_total, L_dna_total);
      p7_osplicescores_GrowTo(oss, sub_M);
      p7_Viterbi_SplicedGlobal_NoP(dsq, om_tr, ox_NoP, oss, 1, L_dna_total, pli->min_intron);
      final_oC = ox_NoP->xmx[L_dna_total * p7X_NXCELLS + p7X_C];
     
      if (fabs(final_gC - final_oC) > 0.001) 
        esl_fatal("%s: generic %.4f != SSE %.4f", msg, final_gC, final_oC); 

      if(final_gC == -eslINFINITY) continue;

      p7_Viterbi_SplicedTrace_NoP(dsq, ox_NoP, gm_tr, pli->splice_scores->signal_scores, otr, 1, L_dna_total, k_start, k_end, pli->min_intron);
      p7_GViterbi_SplicedTrace_NoP(dsq, gm_tr, gx_NoP, pli->splice_scores->signal_scores, gtr, 1, L_dna_total, k_start, k_end, pli->min_intron); 

      p7_trace_Compare(otr, gtr, 0.0);
        
      p7_trace_Reuse(otr);
      p7_trace_Reuse(gtr);
    }

  if (dsq != NULL) free(dsq);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_profile_fs_Destroy(gm_tr);
  p7_fs_oprofile_Destroy(om_tr);
  p7_splicepipeline_Destroy(pli);
  p7_omx_Destroy(ox);
  p7_omx_Destroy(ox_NoP);
  p7_gmx_Destroy(gx_NoP);
  p7_ivx_Destroy(iv_nop);
  p7_trace_fs_Destroy(otr);
  p7_trace_fs_Destroy(gtr);
  p7_osplicescores_Destroy(oss);
}



#endif /*p7VITERBI_SP_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/



/*****************************************************************
 * 5. Test driver.
 *****************************************************************/
#ifdef p7VITERBI_SP_TESTDRIVE
/*
  gcc -g -Wall -msse2 -std=gnu99 -o viterbi_sp_utest -I.. -L.. -I../../easel -L../../easel -Dp7VITERBI_SP_TESTDRIVE viterbi_sp.c -lhmmer -leasel -lm
  ./viterbi_sp_utest
*/
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_sse.h"

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
static char banner[] = "test driver for p7_Viterbi_SplicedGlobal()";

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

  utest_viterbi_sp(r, abcAA, abcDNA, gcode, bgAA, ct, M, N, I);

  p7_bg_Destroy(bgAA);
  p7_codontable_Destroy(ct);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7VITERBI_SP_TESTDRIVE*/
/*----------------- end, test driver ----------------------------*/
