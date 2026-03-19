/* SSE-accelerated spliced Viterbi algorithm.
 *
 * Probability-space DP: multiply (_mm_mul_ps) replaces log-space add;
 * _mm_max_ps replaces ESL_MAX; 0.0 maps to -infinity. Sub-range
 * k_start..k_end emissions and transitions are re-striped into local
 * tables at function entry. Donor-site lookback scores accumulate in
 * OSPLICE_SCORES striped probability-space vectors.
 *
 * Contents:
 *    1. Local table helpers (re-stripe emissions/transitions for sub-range).
 *    2. p7_ospliceviterbi_TranslatedGlobal()
 *    3. p7_ospliceviterbi_TranslatedSemiGlobalExtendDown()
 *    4. p7_ospliceviterbi_TranslatedSemiGlobalExtendUp()
 *    5. Benchmark driver.
 *    6. Unit tests.
 *    7. Test driver.
 */

#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>
#include <emmintrin.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sse.h"

#include "hmmer.h"
#include "impl_sse.h"
#include "../p7_splice.h"

/* P→M transition in probability space: exp(-10.0) */
#define TSC_P_PROB  4.58e-5f

/* tfv layout: per-stripe slots 0..6 = BM,MM,IM,DM,MD,MI,II; DD at 7*Q_full+q */
#define TFV_MM  1
#define TFV_IM  2
#define TFV_DM  3
#define TFV_MD  4
#define TFV_MI  5
#define TFV_II  6
#define TFV_DD  7   /* index into our local table: 7 * Q + q */


/*****************************************************************
 * 1. Local table helpers
 *
 * The full-profile tfv/rfv use Q_full = p7O_NQF(om_fs->M) stripes.
 * We re-pack the sub-range k_start..k_end into Q = p7O_NQF(M_local)
 * local stripes so the inner loop uses simple 0-based indexing.
 *****************************************************************/

/* build_local_rfv():
 * Allocate and fill a [p7P_MAXCODONS1 * Q] emission table.
 * local_rfv[c * Q + q] = emission vector for codon c at local stripe q.
 * *raw_out must be free()'d by caller.
 */
static __m128 *
build_local_rfv(const P7_FS_OPROFILE *om_fs,
                int k_start, int M_local, int Q,
                void **raw_out)
{
  int    Q_full = p7O_NQF(om_fs->M);
  __m128 *base;
  int    c, q_l, r_l, k_l, k_g, q_g, r_g;
  int    status;

  *raw_out = NULL;
  ESL_ALLOC(*raw_out, sizeof(__m128) * p7P_MAXCODONS1 * Q + 15);
  base = (__m128 *)(((unsigned long int) *raw_out + 15) & (~0xfUL));

  for (c = 0; c < p7P_MAXCODONS1; c++) {
    __m128 *rfv = base + c * Q;
    for (q_l = 0; q_l < Q; q_l++) {
      union { __m128 v; float p[4]; } em, tv;
      em.v = _mm_setzero_ps();
      for (r_l = 0; r_l < 4; r_l++) {
        k_l = q_l + r_l * Q + 1;
        if (k_l > M_local) break;
        k_g = k_start + k_l - 1;
        q_g = (k_g - 1) % Q_full;
        r_g = (k_g - 1) / Q_full;
        tv.v       = om_fs->rfv[c][q_g];
        em.p[r_l]  = tv.p[r_g];
      }
      rfv[q_l] = em.v;
    }
  }
  return base;

  ERROR:
    return NULL;
}


/* build_local_tfv():
 * Allocate and fill an [8 * Q] transition table.
 * Slots 1..6 (TFV_MM..TFV_II): base[t * Q + q].
 * Slot 7 (TFV_DD): base[7 * Q + q].
 * *raw_out must be free()'d by caller.
 */
static __m128 *
build_local_tfv(const P7_FS_OPROFILE *om_fs,
                int k_start, int M_local, int Q,
                void **raw_out)
{
  int    Q_full = p7O_NQF(om_fs->M);
  __m128 *base;
  int    q_l, r_l, k_l, k_g, q_g, r_g, t;
  int    status;

  *raw_out = NULL;
  ESL_ALLOC(*raw_out, sizeof(__m128) * 8 * Q + 15);
  base = (__m128 *)(((unsigned long int) *raw_out + 15) & (~0xfUL));

  /* Zero everything (slots 0 = BM is unused but should be 0) */
  for (t = 0; t < 8; t++)
    for (q_l = 0; q_l < Q; q_l++)
      base[t * Q + q_l] = _mm_setzero_ps();

  for (q_l = 0; q_l < Q; q_l++) {
    union { __m128 v; float p[4]; } u[8];
    for (t = 0; t < 8; t++) u[t].v = _mm_setzero_ps();

    for (r_l = 0; r_l < 4; r_l++) {
      union { __m128 v; float p[4]; } tv;
      k_l = q_l + r_l * Q + 1;
      if (k_l > M_local) break;
      k_g = k_start + k_l - 1;
      q_g = (k_g - 1) % Q_full;
      r_g = (k_g - 1) / Q_full;

      for (t = TFV_MM; t <= TFV_II; t++) {
        tv.v         = om_fs->tfv[7 * q_g + t];
        u[t].p[r_l]  = tv.p[r_g];
      }
      /* DD is stored at 7*Q_full + q_g in the full profile */
      tv.v              = om_fs->tfv[7 * Q_full + q_g];
      u[TFV_DD].p[r_l]  = tv.p[r_g];
    }

    for (t = TFV_MM; t <= TFV_II; t++)
      base[t * Q + q_l] = u[t].v;
    base[TFV_DD * Q + q_l] = u[TFV_DD].v;
  }
  return base;

  ERROR:
    return NULL;
}

/* Accessor macros for local tables */
#define LRFV(c, q)   (local_rfv[(c) * Q + (q)])
#define LTFV(t, q)   (local_tfv[(t) * Q + (q)])


/*****************************************************************
 * 2. p7_ospliceviterbi_TranslatedGlobal()
 *****************************************************************/

/* Function:  p7_ospliceviterbi_TranslatedGlobal()
 * Synopsis:  SSE-accelerated fully global translated spliced Viterbi.
 *
 * Purpose:   SSE equivalent of p7_spliceviterbi_TranslatedGlobal().
 *            Uses probability-space DP (multiply/max) with P7_FS_OPROFILE
 *            emissions (codon_lengths == 1) and P7_OMX DP matrix
 *            (nscells = p7X_NSCELLS_SP = 4; validR >= L+1).
 *
 *            Donor-site lookback scores accumulate in <oss->oscore>
 *            (striped __m128 vectors, probability space).
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile not codon_lengths==1, <eslEMEM> on alloc.
 */
int
p7_ospliceviterbi_TranslatedGlobal(SPLICE_PIPELINE *pli, OSPLICE_SCORES *oss,
                                   const ESL_DSQ *sub_dsq,
                                   const P7_FS_OPROFILE *om_fs,
                                   P7_OMX *ox,
                                   int i_start, int i_end,
                                   int k_start, int k_end)
{
  int      M         = k_end - k_start + 1;
  int      L         = i_end - i_start + 1;
  int      Q         = p7O_NQF(M);
  int      r_M       = (M - 1) / Q;   /* lane index of position k=M in the last stripe */
  int      min_intron = pli->min_intron;

  __m128  *local_rfv  = NULL;
  __m128  *local_tfv  = NULL;
  void    *rfv_raw    = NULL;
  void    *tfv_raw    = NULL;

  __m128  *dpc, *dpp3, *dplb;

  __m128   zerov     = _mm_setzero_ps();
  __m128   TSC_P_v   = _mm_set1_ps(TSC_P_PROB);

  /* Carry variables */
  __m128   mpv3, dpv3, ipv3, ppv3;   /* right-shifted row i-3, for k-1 transition  */
  __m128   tsv_m, tsv_d;              /* donor lookback right-shifted carry          */
  __m128   msv, isv, dcv;
  __m128   sig_GTAG_v, sig_GCAG_v, sig_ATAC_v;

  /* Acceptor/donor validity (probability space: 1.0 valid, 0.0 invalid) */
  int      acc0, acc1, acc2;           /* ACCEPT_AG, ACCEPT_AC, or -1 */
  int      donor0, donor1, donor2;       /* DONOR_GT/GC/AT, or -1       */
  int      do_pstate, do_donor;
  __m128   AG0v, AG1v, AG2v, AC0v, AC1v, AC2v;
  __m128   GT0v, GC0v, AT0v, GT1v, GC1v, AT1v, GT2v, GC2v, AT2v;

  /* Codon indices */
  int      C0, C1[4], C2[16];

  /* Nucleotide windows:
   *   v, w, x: current codon (positions i-2, i-1, i)
   *   r, s, t, u: donor-lookback window (positions i-min_intron-2..i-min_intron+1)
   */
  int      v, w, x;
  int      r, s, t, u;
  int      nuc1, nuc2;

  float    xB_0;
  int      i, q, j, loop_end, sub_i;
  int      status;

  if (om_fs->codon_lengths != 1)
    ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");

  /* Build local emission and transition tables for sub-range k_start..k_end */
  local_rfv = build_local_rfv(om_fs, k_start, M, Q, &rfv_raw);
  local_tfv = build_local_tfv(om_fs, k_start, M, Q, &tfv_raw);
  if (!local_rfv || !local_tfv) { status = eslEMEM; goto ERROR; }

  sig_GTAG_v = _mm_set1_ps(oss->signal_scores[p7S_GTAG]);
  sig_GCAG_v = _mm_set1_ps(oss->signal_scores[p7S_GCAG]);
  sig_ATAC_v = _mm_set1_ps(oss->signal_scores[p7S_ATAC]);

  /* Reset OSS P-state accumulation tables to 0.0 (= -inf in log space) */
  for (q = 0; q < Q; q++)
    for (j = 0; j < SIGNAL_MEM_SIZE; j++)
      oss->oscore[q][j] = zerov;

  /* B(0) in probability space */
  xB_0 = om_fs->xf[p7O_N][p7O_MOVE];

  /* Zero rows 0, 1, 2 */
  for (i = 0; i <= 2; i++) {
    dpc = ox->dpf[i];
    for (q = 0; q < Q; q++)
      MMO_SP(dpc,q) = DMO_SP(dpc,q) = IMO_SP(dpc,q) = PMO_SP(dpc,q) = zerov;
  }

  /* Initialize acceptor/donor rolling window to impossible */
  acc0 = acc1 = acc2 = -1;

  /* Initialize codon nucleotide window to placeholder (-1 → clamped to degenerate) */
  v = w = x = -1;

  /* ---------------------------------------------------------------
   * Rows 1-2: update w,x for the codon window (no DP — rows zeroed).
   * --------------------------------------------------------------- */
  for (i = 1; i <= 2; i++) {
    w = x;
    sub_i = i_start + i - 1;
    x = (sub_dsq[sub_i] < 4) ? (int)sub_dsq[sub_i] : p7P_MAXCODONS1;
  }

  /* ---------------------------------------------------------------
   * Rows 3 .. min_intron+2: M/I/D recursion only.
   * No P state (OSS arrays are 0) and no donor lookback yet.
   * --------------------------------------------------------------- */
  loop_end = ESL_MIN(L, min_intron + 2);

  for (i = 3; i <= loop_end; i++) {
    v = w;  w = x;
    sub_i = i_start + i - 1;
    x = (sub_dsq[sub_i] < 4) ? (int)sub_dsq[sub_i] : p7P_MAXCODONS1;

    C0 = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);
    __m128 *rfv_c0 = local_rfv + C0 * Q;

    dpc  = ox->dpf[i];
    dpp3 = ox->dpf[i - 3];

    mpv3 = esl_sse_rightshiftz_float(MMO_SP(dpp3, Q-1));
    dpv3 = esl_sse_rightshiftz_float(DMO_SP(dpp3, Q-1));
    ipv3 = esl_sse_rightshiftz_float(IMO_SP(dpp3, Q-1));
    ppv3 = esl_sse_rightshiftz_float(PMO_SP(dpp3, Q-1));
    dcv  = zerov;

    for (q = 0; q < Q; q++) {
      /* M(i,k) = max(M(i-3,k-1)*MM, I(i-3,k-1)*IM, D(i-3,k-1)*DM, P(i-3,k-1)*TSC_P)
       *          * emission(k, C0) */
      msv = _mm_max_ps(_mm_mul_ps(mpv3, LTFV(TFV_MM, q)),
            _mm_max_ps(_mm_mul_ps(ipv3, LTFV(TFV_IM, q)),
            _mm_max_ps(_mm_mul_ps(dpv3, LTFV(TFV_DM, q)),
                       _mm_mul_ps(ppv3, TSC_P_v))));
      msv = _mm_mul_ps(msv, rfv_c0[q]);

      /* Update k-1 carry */
      mpv3 = MMO_SP(dpp3, q);
      dpv3 = DMO_SP(dpp3, q);
      ipv3 = IMO_SP(dpp3, q);
      ppv3 = PMO_SP(dpp3, q);

      MMO_SP(dpc, q) = msv;
      DMO_SP(dpc, q) = dcv;               /* delayed MD carry from stripe q-1 */
      dcv             = _mm_mul_ps(msv, LTFV(TFV_MD, q));

      /* I(i,k) = max(M(i-3,k)*MI, I(i-3,k)*II); zero at stop codons */
      isv = _mm_max_ps(_mm_mul_ps(MMO_SP(dpp3, q), LTFV(TFV_MI, q)),
                       _mm_mul_ps(IMO_SP(dpp3, q), LTFV(TFV_II, q)));
      IMO_SP(dpc, q) = _mm_and_ps(isv, _mm_cmpgt_ps(rfv_c0[q], zerov));

      PMO_SP(dpc, q) = zerov;
    }

    /* Global entry at i=3, k=1 (local position 1 = stripe q=0, lane r=0):
     * M(3,1) = B(0) * emission(k_start, C0).
     * All other positions got 0 from the zero dpp3 row. */
    if (i == 3) {
      union { __m128 v; float p[4]; } em, mc, md0, d1;
      em.v    = rfv_c0[0];
      mc.v    = MMO_SP(dpc, 0);
      mc.p[0] = xB_0 * em.p[0];
      MMO_SP(dpc, 0) = mc.v;
      /* The M loop computed dcv=0 from M(3,1)=0 (before this override), so
       * DMO_SP(dpc,1).lane[0] = 0, which would leave D(3,2..M)=0 in SSE.
       * Inject the correct M(3,1)->D(3,2) carry into DMO_SP(dpc,1).lane[0]. */
      md0.v   = LTFV(TFV_MD, 0);
      d1.v    = DMO_SP(dpc, 1);
      d1.p[0] = mc.p[0] * md0.p[0];
      DMO_SP(dpc, 1) = d1.v;
    }

    /* D-state serialization:
     * After the q loop, dcv holds M(i,Q-1)*MD.  DMO_SP(dpc,q) currently stores
     * the delayed MD carry from stripe q-1, i.e., DMO_SP(dpc,q)=dcv_q-1.
     * We finalize D by: for each q, D(q) = max(dcv_from_prev, D(q)); carry DD. */
    dcv = esl_sse_rightshiftz_float(dcv);
    DMO_SP(dpc, 0) = zerov;
    for (q = 0; q < Q; q++) {
      DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
      dcv = _mm_mul_ps(DMO_SP(dpc, q), LTFV(TFV_DD, q));
    }
    for (j = 1; j < 4; j++) {
      dcv = esl_sse_rightshiftz_float(dcv);
      for (q = 0; q < Q; q++) {
        DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
        dcv = _mm_mul_ps(DMO_SP(dpc, q), LTFV(TFV_DD, q));
      }
    }
  }  /* end early rows */


  /* ---------------------------------------------------------------
   * Initialize donor-lookback nucleotide window (r,s,t,u) for the
   * main DP.  At the first main-loop iteration i = min_intron+3,
   * after "r=s; s=t; t=u; u=sub_dsq[sub_i-min_intron+1]" we need:
   *   r = sub_dsq[i_start], s = sub_dsq[i_start+1],
   *   t = sub_dsq[i_start+2], u = sub_dsq[i_start+3].
   * Pre-load the three prior values:
   * --------------------------------------------------------------- */
  s = (int)sub_dsq[i_start];
  t = (int)sub_dsq[i_start + 1];
  u = (int)sub_dsq[i_start + 2];
  r = 0;  /* will be overwritten by r=s at loop start */

  /* ---------------------------------------------------------------
   * Main DP: rows min_intron+3 .. L.
   * Full M/I/D/P update + donor lookback SSX write.
   * --------------------------------------------------------------- */
  for (i = min_intron + 3; i <= L; i++) {
    __m128 pv, TMP_v;

    sub_i = i_start + i - 1;

    /* Advance donor-lookback window */
    r = s;  s = t;  t = u;
    u = (sub_dsq[sub_i - min_intron + 1] < 4)
        ? (int)sub_dsq[sub_i - min_intron + 1]
        : p7P_MAXCODONS1;

    /* Advance codon window */
    v = w;  w = x;
    x = (sub_dsq[sub_i] < 4) ? (int)sub_dsq[sub_i] : p7P_MAXCODONS1;

    C0 = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);

    __m128 *rfv_c0 = local_rfv + C0 * Q;
    __m128 *rfv_c1[4];
    __m128 *rfv_c2[16];

    /* Shift acceptor window */
    acc0 = acc1;  acc1 = acc2;
    if (v < 0 || v >= 4 || w < 0 || w >= 4)      acc2 = -1;
    else if (SIGNAL(v,w) == ACCEPT_AG)             acc2 = ACCEPT_AG;
    else if (SIGNAL(v,w) == ACCEPT_AC)             acc2 = ACCEPT_AC;
    else                                           acc2 = -1;

    /* Donor validity for all three codon classes */
    donor0 = donor1 = donor2 = -1;
    if (!(r < 0 || r >= 4 || s < 0 || s >= 4)) {
      if      (SIGNAL(r,s) == DONOR_GT) donor0 = DONOR_GT;
      else if (SIGNAL(r,s) == DONOR_GC) donor0 = DONOR_GC;
      else if (SIGNAL(r,s) == DONOR_AT) donor0 = DONOR_AT;
      if (!(t < 0 || t >= 4)) {
        if      (SIGNAL(s,t) == DONOR_GT) donor1 = DONOR_GT;
        else if (SIGNAL(s,t) == DONOR_GC) donor1 = DONOR_GC;
        else if (SIGNAL(s,t) == DONOR_AT) donor1 = DONOR_AT;
        if (!(u < 0 || u >= 4)) {
          if      (SIGNAL(t,u) == DONOR_GT) donor2 = DONOR_GT;
          else if (SIGNAL(t,u) == DONOR_GC) donor2 = DONOR_GC;
          else if (SIGNAL(t,u) == DONOR_AT) donor2 = DONOR_AT;
        }
      }
    }

    do_pstate = (acc0 >= 0 || acc1 >= 0 || acc2 >= 0);
    do_donor  = (donor0 >= 0 || donor1 >= 0 || donor2 >= 0);

    if (do_pstate) {
      for (nuc1 = 0; nuc1 < 4; nuc1++) {
        C1[nuc1] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, w, x), p7P_DEGEN1_C);
        rfv_c1[nuc1] = local_rfv + C1[nuc1] * Q;
        for (nuc2 = 0; nuc2 < 4; nuc2++) {
          C2[nuc1*4+nuc2] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, nuc2, x), p7P_DEGEN1_C);
          rfv_c2[nuc1*4+nuc2] = local_rfv + C2[nuc1*4+nuc2] * Q;
        }
      }
      AG0v = _mm_set1_ps(acc0 == ACCEPT_AG ? 1.0f : 0.0f);
      AG1v = _mm_set1_ps(acc1 == ACCEPT_AG ? 1.0f : 0.0f);
      AG2v = _mm_set1_ps(acc2 == ACCEPT_AG ? 1.0f : 0.0f);
      AC0v = _mm_set1_ps(acc0 == ACCEPT_AC ? 1.0f : 0.0f);
      AC1v = _mm_set1_ps(acc1 == ACCEPT_AC ? 1.0f : 0.0f);
      AC2v = _mm_set1_ps(acc2 == ACCEPT_AC ? 1.0f : 0.0f);
    }
    if (do_donor) {
      GT0v = _mm_set1_ps(donor0 == DONOR_GT ? 1.0f : 0.0f);
      GC0v = _mm_set1_ps(donor0 == DONOR_GC ? 1.0f : 0.0f);
      AT0v = _mm_set1_ps(donor0 == DONOR_AT ? 1.0f : 0.0f);
      GT1v = _mm_set1_ps(donor1 == DONOR_GT ? 1.0f : 0.0f);
      GC1v = _mm_set1_ps(donor1 == DONOR_GC ? 1.0f : 0.0f);
      AT1v = _mm_set1_ps(donor1 == DONOR_AT ? 1.0f : 0.0f);
      GT2v = _mm_set1_ps(donor2 == DONOR_GT ? 1.0f : 0.0f);
      GC2v = _mm_set1_ps(donor2 == DONOR_GC ? 1.0f : 0.0f);
      AT2v = _mm_set1_ps(donor2 == DONOR_AT ? 1.0f : 0.0f);
    }

    dpc  = ox->dpf[i];
    dpp3 = ox->dpf[i - 3];
    dplb = ox->dpf[i - min_intron - 3];

    /* Initialize right-shifted carries for k-1 transitions */
    mpv3  = esl_sse_rightshiftz_float(MMO_SP(dpp3, Q-1));
    dpv3  = esl_sse_rightshiftz_float(DMO_SP(dpp3, Q-1));
    ipv3  = esl_sse_rightshiftz_float(IMO_SP(dpp3, Q-1));
    ppv3  = esl_sse_rightshiftz_float(PMO_SP(dpp3, Q-1));
    tsv_m = esl_sse_rightshiftz_float(MMO_SP(dplb, Q-1));
    tsv_d = esl_sse_rightshiftz_float(DMO_SP(dplb, Q-1));
    dcv   = zerov;

    /* Main stripe loop: k = 1 .. M-1 (q = 0..Q-2) and k = M (q = Q-1).
     * The last stripe (k=M) skips I state exit. */
    for (q = 0; q < Q; q++) {
      int is_last = (q == Q - 1);

      /* ---- M state ---- */
      msv = _mm_max_ps(_mm_mul_ps(mpv3, LTFV(TFV_MM, q)),
            _mm_max_ps(_mm_mul_ps(ipv3, LTFV(TFV_IM, q)),
            _mm_max_ps(_mm_mul_ps(dpv3, LTFV(TFV_DM, q)),
                       _mm_mul_ps(ppv3, TSC_P_v))));
      msv = _mm_mul_ps(msv, rfv_c0[q]);

      mpv3 = MMO_SP(dpp3, q);
      dpv3 = DMO_SP(dpp3, q);
      ipv3 = IMO_SP(dpp3, q);
      ppv3 = PMO_SP(dpp3, q);

      MMO_SP(dpc, q) = msv;
      DMO_SP(dpc, q) = dcv;
      dcv = _mm_mul_ps(msv, LTFV(TFV_MD, q));

      /* ---- I state (zero at k=M only, not entire last stripe) ---- */
      isv = _mm_max_ps(_mm_mul_ps(MMO_SP(dpp3, q), LTFV(TFV_MI, q)),
                       _mm_mul_ps(IMO_SP(dpp3, q), LTFV(TFV_II, q)));
      if (is_last) {
        union { __m128 v; float p[4]; } ii;
        ii.v = _mm_and_ps(isv, _mm_cmpgt_ps(rfv_c0[q], zerov));
        ii.p[r_M] = 0.0f;  /* zero only position k=M; other lanes (k<M) remain valid */
        IMO_SP(dpc, q) = ii.v;
      } else {
        IMO_SP(dpc, q) = _mm_and_ps(isv, _mm_cmpgt_ps(rfv_c0[q], zerov));
      }

      /* ---- P state ---- */
      PMO_SP(dpc, q) = zerov;
      if (do_pstate) {
        pv = zerov;

        /* C0: full 3-nt from acceptor */
        TMP_v = _mm_max_ps(_mm_mul_ps(OSS0(oss,q,p7S_GTAG), _mm_mul_ps(AG0v, sig_GTAG_v)),
                _mm_max_ps(_mm_mul_ps(OSS0(oss,q,p7S_GCAG), _mm_mul_ps(AG0v, sig_GCAG_v)),
                           _mm_mul_ps(OSS0(oss,q,p7S_ATAC), _mm_mul_ps(AC0v, sig_ATAC_v))));
        pv = _mm_max_ps(pv, _mm_mul_ps(TMP_v, rfv_c0[q]));

        /* C1: 1 nuc from before donor + 2 from acceptor */
        for (nuc1 = 0; nuc1 < 4; nuc1++) {
          TMP_v = _mm_max_ps(_mm_mul_ps(OSS1(oss,q,p7S_GTAG,nuc1), _mm_mul_ps(AG1v, sig_GTAG_v)),
                  _mm_max_ps(_mm_mul_ps(OSS1(oss,q,p7S_GCAG,nuc1), _mm_mul_ps(AG1v, sig_GCAG_v)),
                             _mm_mul_ps(OSS1(oss,q,p7S_ATAC,nuc1), _mm_mul_ps(AC1v, sig_ATAC_v))));
          pv = _mm_max_ps(pv, _mm_mul_ps(TMP_v, rfv_c1[nuc1][q]));
        }

        /* C2: 2 nucs from before donor + 1 from acceptor */
        for (nuc1 = 0; nuc1 < 4; nuc1++) {
          for (nuc2 = 0; nuc2 < 4; nuc2++) {
            TMP_v = _mm_max_ps(_mm_mul_ps(OSS2(oss,q,p7S_GTAG,nuc1,nuc2), _mm_mul_ps(AG2v, sig_GTAG_v)),
                    _mm_max_ps(_mm_mul_ps(OSS2(oss,q,p7S_GCAG,nuc1,nuc2), _mm_mul_ps(AG2v, sig_GCAG_v)),
                               _mm_mul_ps(OSS2(oss,q,p7S_ATAC,nuc1,nuc2), _mm_mul_ps(AC2v, sig_ATAC_v))));
            pv = _mm_max_ps(pv, _mm_mul_ps(TMP_v, rfv_c2[nuc1*4+nuc2][q]));
          }
        }
        PMO_SP(dpc, q) = pv;
      }

      /* ---- Donor lookback: update OSS at position k from (lb, k-1) ---- */
      if (do_donor) {
        TMP_v = _mm_max_ps(tsv_m, tsv_d);

        OSS0(oss,q,p7S_GTAG) = _mm_max_ps(OSS0(oss,q,p7S_GTAG), _mm_mul_ps(TMP_v, GT0v));
        OSS0(oss,q,p7S_GCAG) = _mm_max_ps(OSS0(oss,q,p7S_GCAG), _mm_mul_ps(TMP_v, GC0v));
        OSS0(oss,q,p7S_ATAC) = _mm_max_ps(OSS0(oss,q,p7S_ATAC), _mm_mul_ps(TMP_v, AT0v));

        OSS1(oss,q,p7S_GTAG,r) = _mm_max_ps(OSS1(oss,q,p7S_GTAG,r), _mm_mul_ps(TMP_v, GT1v));
        OSS1(oss,q,p7S_GCAG,r) = _mm_max_ps(OSS1(oss,q,p7S_GCAG,r), _mm_mul_ps(TMP_v, GC1v));
        OSS1(oss,q,p7S_ATAC,r) = _mm_max_ps(OSS1(oss,q,p7S_ATAC,r), _mm_mul_ps(TMP_v, AT1v));

        OSS2(oss,q,p7S_GTAG,r,s) = _mm_max_ps(OSS2(oss,q,p7S_GTAG,r,s), _mm_mul_ps(TMP_v, GT2v));
        OSS2(oss,q,p7S_GCAG,r,s) = _mm_max_ps(OSS2(oss,q,p7S_GCAG,r,s), _mm_mul_ps(TMP_v, GC2v));
        OSS2(oss,q,p7S_ATAC,r,s) = _mm_max_ps(OSS2(oss,q,p7S_ATAC,r,s), _mm_mul_ps(TMP_v, AT2v));
      }

      /* Advance donor lookback carry */
      tsv_m = MMO_SP(dplb, q);
      tsv_d = DMO_SP(dplb, q);
    }

    /* D-state serialization */
    dcv = esl_sse_rightshiftz_float(dcv);
    DMO_SP(dpc, 0) = zerov;
    for (q = 0; q < Q; q++) {
      DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
      dcv = _mm_mul_ps(DMO_SP(dpc, q), LTFV(TFV_DD, q));
    }
    for (j = 1; j < 4; j++) {
      dcv = esl_sse_rightshiftz_float(dcv);
      for (q = 0; q < Q; q++) {
        DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
        dcv = _mm_mul_ps(DMO_SP(dpc, q), LTFV(TFV_DD, q));
      }
    }
  }  /* end main DP loop */


  /* ---------------------------------------------------------------
   * Exit: E(L) = max(M(L,M), D(L,M)).
   * Position k=M is at local stripe q_M = (M-1) % Q, lane r_M = (M-1) / Q.
   * --------------------------------------------------------------- */
  {
    int   q_M = (M - 1) % Q;
    int   r_M = (M - 1) / Q;
    float e_M, e_D, xE;
    union { __m128 v; float p[4]; } u;

    dpc = ox->dpf[L];

    u.v = MMO_SP(dpc, q_M);  e_M = u.p[r_M];
    u.v = DMO_SP(dpc, q_M);  e_D = u.p[r_M];
    xE  = ESL_MAX(e_M, e_D);

    ox->xmx[L * p7X_NXCELLS + p7X_E] = xE;
    ox->xmx[L * p7X_NXCELLS + p7X_C] = xE * om_fs->xf[p7O_E][p7O_MOVE];
  }

  ox->M = M;
  ox->L = L;

  free(rfv_raw);
  free(tfv_raw);
  return eslOK;

  ERROR:
    if (rfv_raw) free(rfv_raw);
    if (tfv_raw) free(tfv_raw);
    return status;
}
/*-------------- end, p7_ospliceviterbi_TranslatedGlobal() ------*/


/*****************************************************************
 * 3. p7_ospliceviterbi_TranslatedSemiGlobalExtendDown()
 *****************************************************************/

/* Function:  p7_ospliceviterbi_TranslatedSemiGlobalExtendDown()
 * Synopsis:  SSE-accelerated semi-global (global entry, semi-global exit) translated spliced Viterbi.
 *
 * Purpose:   SSE equivalent of p7_spliceviterbi_TranslatedSemiGlobalExtendDown().
 *            Same DP recurrence as p7_ospliceviterbi_TranslatedGlobal() but exits
 *            semi-globally: at each row i, E(i) = max over k=1..M of max(M(i,k), D(i,k)),
 *            and C(i) = max(C(i-3) * CC_loop, E(i) * EC_move).  The final score is C(L).
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile not codon_lengths==1, <eslEMEM> on alloc.
 */
int
p7_ospliceviterbi_TranslatedSemiGlobalExtendDown(SPLICE_PIPELINE *pli, OSPLICE_SCORES *oss,
                                                  const ESL_DSQ *sub_dsq,
                                                  const P7_FS_OPROFILE *om_fs,
                                                  P7_OMX *ox,
                                                  int i_start, int i_end,
                                                  int k_start, int k_end)
{
  int      M          = k_end - k_start + 1;
  int      L          = i_end - i_start + 1;
  int      Q          = p7O_NQF(M);
  int      r_M        = (M - 1) / Q;
  int      min_intron = pli->min_intron;

  __m128  *local_rfv  = NULL;
  __m128  *local_tfv  = NULL;
  void    *rfv_raw    = NULL;
  void    *tfv_raw    = NULL;
  __m128  *valid_mask = NULL;
  void    *mask_raw   = NULL;

  __m128  *dpc, *dpp3, *dplb;

  __m128   zerov     = _mm_setzero_ps();
  __m128   TSC_P_v   = _mm_set1_ps(TSC_P_PROB);

  __m128   mpv3, dpv3, ipv3, ppv3;
  __m128   tsv_m, tsv_d;
  __m128   msv, isv, dcv;
  __m128   sig_GTAG_v, sig_GCAG_v, sig_ATAC_v;

  int      acc0, acc1, acc2;           /* ACCEPT_AG, ACCEPT_AC, or -1 */
  int      donor0, donor1, donor2;       /* DONOR_GT/GC/AT, or -1       */
  int      do_pstate, do_donor;
  __m128   AG0v, AG1v, AG2v, AC0v, AC1v, AC2v;
  __m128   GT0v, GC0v, AT0v, GT1v, GC1v, AT1v, GT2v, GC2v, AT2v;

  int      C0, C1[4], C2[16];
  int      v, w, x;
  int      r, s, t, u;
  int      nuc1, nuc2;

  float    xB_0;
  float    CC_loop, EC_move;
  int      i, q, j, loop_end, sub_i;
  int      status;

  if (om_fs->codon_lengths != 1)
    ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");

  local_rfv = build_local_rfv(om_fs, k_start, M, Q, &rfv_raw);
  local_tfv = build_local_tfv(om_fs, k_start, M, Q, &tfv_raw);
  if (!local_rfv || !local_tfv) { status = eslEMEM; goto ERROR; }

  /* Valid-lane bitmask per stripe: all-1s for k<=M, all-0s for k>M.
   * Used to exclude invalid lanes (k>M in last stripe) from the E horizontal max. */
  ESL_ALLOC(mask_raw, sizeof(__m128) * Q + 15);
  valid_mask = (__m128 *)(((unsigned long int) mask_raw + 15) & (~0xfUL));
  {
    int rl;
    for (q = 0; q < Q; q++) {
      union { __m128 v; float p[4]; } vm;
      vm.v = _mm_setzero_ps();
      for (rl = 0; rl < 4; rl++)
        if (q + rl * Q + 1 <= M) vm.p[rl] = 1.0f;
      valid_mask[q] = _mm_cmpgt_ps(vm.v, zerov);
    }
  }

  sig_GTAG_v = _mm_set1_ps(oss->signal_scores[p7S_GTAG]);
  sig_GCAG_v = _mm_set1_ps(oss->signal_scores[p7S_GCAG]);
  sig_ATAC_v = _mm_set1_ps(oss->signal_scores[p7S_ATAC]);

  CC_loop = om_fs->xf[p7O_C][p7O_LOOP];
  EC_move = om_fs->xf[p7O_E][p7O_MOVE];

  for (q = 0; q < Q; q++)
    for (j = 0; j < SIGNAL_MEM_SIZE; j++)
      oss->oscore[q][j] = zerov;

  xB_0 = om_fs->xf[p7O_N][p7O_MOVE];

  /* Zero rows 0, 1, 2: DP cells and C state */
  for (i = 0; i <= 2; i++) {
    dpc = ox->dpf[i];
    for (q = 0; q < Q; q++)
      MMO_SP(dpc,q) = DMO_SP(dpc,q) = IMO_SP(dpc,q) = PMO_SP(dpc,q) = zerov;
    ox->xmx[i * p7X_NXCELLS + p7X_E] = 0.0f;
    ox->xmx[i * p7X_NXCELLS + p7X_C] = 0.0f;
  }

  acc0 = acc1 = acc2 = -1;
  v = w = x = -1;

  for (i = 1; i <= 2; i++) {
    w = x;
    sub_i = i_start + i - 1;
    x = (sub_dsq[sub_i] < 4) ? (int)sub_dsq[sub_i] : p7P_MAXCODONS1;
  }

  /* ---------------------------------------------------------------
   * Early rows 3..min_intron+2: M/I/D only, no donor lookback (P=0).
   * Accumulate E and C at each row.
   * --------------------------------------------------------------- */
  loop_end = ESL_MIN(L, min_intron + 2);

  for (i = 3; i <= loop_end; i++) {
    v = w;  w = x;
    sub_i = i_start + i - 1;
    x = (sub_dsq[sub_i] < 4) ? (int)sub_dsq[sub_i] : p7P_MAXCODONS1;

    C0 = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);
    __m128 *rfv_c0 = local_rfv + C0 * Q;

    /* Acceptor window: needed for correct P state when main loop fires early */
    acc0 = acc1;  acc1 = acc2;
    if (v < 0 || v >= 4 || w < 0 || w >= 4)      acc2 = -1;
    else if (SIGNAL(v,w) == ACCEPT_AG)             acc2 = ACCEPT_AG;
    else if (SIGNAL(v,w) == ACCEPT_AC)             acc2 = ACCEPT_AC;
    else                                           acc2 = -1;

    dpc  = ox->dpf[i];
    dpp3 = ox->dpf[i - 3];

    mpv3 = esl_sse_rightshiftz_float(MMO_SP(dpp3, Q-1));
    dpv3 = esl_sse_rightshiftz_float(DMO_SP(dpp3, Q-1));
    ipv3 = esl_sse_rightshiftz_float(IMO_SP(dpp3, Q-1));
    ppv3 = esl_sse_rightshiftz_float(PMO_SP(dpp3, Q-1));
    dcv  = zerov;

    for (q = 0; q < Q; q++) {
      int is_last = (q == Q - 1);

      msv = _mm_max_ps(_mm_mul_ps(mpv3, LTFV(TFV_MM, q)),
            _mm_max_ps(_mm_mul_ps(ipv3, LTFV(TFV_IM, q)),
            _mm_max_ps(_mm_mul_ps(dpv3, LTFV(TFV_DM, q)),
                       _mm_mul_ps(ppv3, TSC_P_v))));
      msv = _mm_mul_ps(msv, rfv_c0[q]);

      mpv3 = MMO_SP(dpp3, q);
      dpv3 = DMO_SP(dpp3, q);
      ipv3 = IMO_SP(dpp3, q);
      ppv3 = PMO_SP(dpp3, q);

      MMO_SP(dpc, q) = msv;
      DMO_SP(dpc, q) = dcv;
      dcv             = _mm_mul_ps(msv, LTFV(TFV_MD, q));

      isv = _mm_max_ps(_mm_mul_ps(MMO_SP(dpp3, q), LTFV(TFV_MI, q)),
                       _mm_mul_ps(IMO_SP(dpp3, q), LTFV(TFV_II, q)));
      if (is_last) {
        union { __m128 v; float p[4]; } ii;
        ii.v = _mm_and_ps(isv, _mm_cmpgt_ps(rfv_c0[q], zerov));
        ii.p[r_M] = 0.0f;
        IMO_SP(dpc, q) = ii.v;
      } else {
        IMO_SP(dpc, q) = _mm_and_ps(isv, _mm_cmpgt_ps(rfv_c0[q], zerov));
      }

      PMO_SP(dpc, q) = zerov;
    }

    if (i == 3) {
      union { __m128 v; float p[4]; } em, mc, md0, d1;
      em.v    = rfv_c0[0];
      mc.v    = MMO_SP(dpc, 0);
      mc.p[0] = xB_0 * em.p[0];
      MMO_SP(dpc, 0) = mc.v;
      md0.v   = LTFV(TFV_MD, 0);
      d1.v    = DMO_SP(dpc, 1);
      d1.p[0] = mc.p[0] * md0.p[0];
      DMO_SP(dpc, 1) = d1.v;
    }

    /* D-state serialization */
    dcv = esl_sse_rightshiftz_float(dcv);
    DMO_SP(dpc, 0) = zerov;
    for (q = 0; q < Q; q++) {
      DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
      dcv = _mm_mul_ps(DMO_SP(dpc, q), LTFV(TFV_DD, q));
    }
    for (j = 1; j < 4; j++) {
      dcv = esl_sse_rightshiftz_float(dcv);
      for (q = 0; q < Q; q++) {
        DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
        dcv = _mm_mul_ps(DMO_SP(dpc, q), LTFV(TFV_DD, q));
      }
    }

    /* Semi-global exit: E(i) = max over k=1..M of max(M,D); C(i) = max(C(i-3)*CC, E(i)*EC) */
    {
      __m128 emax_v = zerov;
      union { __m128 v; float p[4]; } em;
      float xE, xC;
      for (q = 0; q < Q; q++)
        emax_v = _mm_max_ps(emax_v, _mm_and_ps(_mm_max_ps(MMO_SP(dpc, q), DMO_SP(dpc, q)), valid_mask[q]));
      em.v = emax_v;
      xE   = ESL_MAX(ESL_MAX(em.p[0], em.p[1]), ESL_MAX(em.p[2], em.p[3]));
      xC   = ESL_MAX(ox->xmx[(i - 3) * p7X_NXCELLS + p7X_C] * CC_loop, xE * EC_move);
      ox->xmx[i * p7X_NXCELLS + p7X_E] = xE;
      ox->xmx[i * p7X_NXCELLS + p7X_C] = xC;
    }
  }  /* end early rows */


  /* ---------------------------------------------------------------
   * Initialize donor-lookback nucleotide window for the main loop.
   * --------------------------------------------------------------- */
  s = (int)sub_dsq[i_start];
  t = (int)sub_dsq[i_start + 1];
  u = (int)sub_dsq[i_start + 2];
  r = 0;

  /* ---------------------------------------------------------------
   * Main DP: rows min_intron+3 .. L.
   * Full M/I/D/P + donor lookback OSS write + semi-global E/C exit.
   * --------------------------------------------------------------- */
  for (i = min_intron + 3; i <= L; i++) {
    __m128 pv, TMP_v;

    sub_i = i_start + i - 1;

    r = s;  s = t;  t = u;
    u = (sub_dsq[sub_i - min_intron + 1] < 4)
        ? (int)sub_dsq[sub_i - min_intron + 1]
        : p7P_MAXCODONS1;

    v = w;  w = x;
    x = (sub_dsq[sub_i] < 4) ? (int)sub_dsq[sub_i] : p7P_MAXCODONS1;

    C0 = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);

    __m128 *rfv_c0 = local_rfv + C0 * Q;
    __m128 *rfv_c1[4];
    __m128 *rfv_c2[16];

    acc0 = acc1;  acc1 = acc2;
    if (v < 0 || v >= 4 || w < 0 || w >= 4)      acc2 = -1;
    else if (SIGNAL(v,w) == ACCEPT_AG)             acc2 = ACCEPT_AG;
    else if (SIGNAL(v,w) == ACCEPT_AC)             acc2 = ACCEPT_AC;
    else                                           acc2 = -1;

    donor0 = donor1 = donor2 = -1;
    if (!(r < 0 || r >= 4 || s < 0 || s >= 4)) {
      if      (SIGNAL(r,s) == DONOR_GT) donor0 = DONOR_GT;
      else if (SIGNAL(r,s) == DONOR_GC) donor0 = DONOR_GC;
      else if (SIGNAL(r,s) == DONOR_AT) donor0 = DONOR_AT;
      if (!(t < 0 || t >= 4)) {
        if      (SIGNAL(s,t) == DONOR_GT) donor1 = DONOR_GT;
        else if (SIGNAL(s,t) == DONOR_GC) donor1 = DONOR_GC;
        else if (SIGNAL(s,t) == DONOR_AT) donor1 = DONOR_AT;
        if (!(u < 0 || u >= 4)) {
          if      (SIGNAL(t,u) == DONOR_GT) donor2 = DONOR_GT;
          else if (SIGNAL(t,u) == DONOR_GC) donor2 = DONOR_GC;
          else if (SIGNAL(t,u) == DONOR_AT) donor2 = DONOR_AT;
        }
      }
    }

    do_pstate = (acc0 >= 0 || acc1 >= 0 || acc2 >= 0);
    do_donor  = (donor0 >= 0 || donor1 >= 0 || donor2 >= 0);

    if (do_pstate) {
      for (nuc1 = 0; nuc1 < 4; nuc1++) {
        C1[nuc1] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, w, x), p7P_DEGEN1_C);
        rfv_c1[nuc1] = local_rfv + C1[nuc1] * Q;
        for (nuc2 = 0; nuc2 < 4; nuc2++) {
          C2[nuc1*4+nuc2] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, nuc2, x), p7P_DEGEN1_C);
          rfv_c2[nuc1*4+nuc2] = local_rfv + C2[nuc1*4+nuc2] * Q;
        }
      }
      AG0v = _mm_set1_ps(acc0 == ACCEPT_AG ? 1.0f : 0.0f);
      AG1v = _mm_set1_ps(acc1 == ACCEPT_AG ? 1.0f : 0.0f);
      AG2v = _mm_set1_ps(acc2 == ACCEPT_AG ? 1.0f : 0.0f);
      AC0v = _mm_set1_ps(acc0 == ACCEPT_AC ? 1.0f : 0.0f);
      AC1v = _mm_set1_ps(acc1 == ACCEPT_AC ? 1.0f : 0.0f);
      AC2v = _mm_set1_ps(acc2 == ACCEPT_AC ? 1.0f : 0.0f);
    }
    if (do_donor) {
      GT0v = _mm_set1_ps(donor0 == DONOR_GT ? 1.0f : 0.0f);
      GC0v = _mm_set1_ps(donor0 == DONOR_GC ? 1.0f : 0.0f);
      AT0v = _mm_set1_ps(donor0 == DONOR_AT ? 1.0f : 0.0f);
      GT1v = _mm_set1_ps(donor1 == DONOR_GT ? 1.0f : 0.0f);
      GC1v = _mm_set1_ps(donor1 == DONOR_GC ? 1.0f : 0.0f);
      AT1v = _mm_set1_ps(donor1 == DONOR_AT ? 1.0f : 0.0f);
      GT2v = _mm_set1_ps(donor2 == DONOR_GT ? 1.0f : 0.0f);
      GC2v = _mm_set1_ps(donor2 == DONOR_GC ? 1.0f : 0.0f);
      AT2v = _mm_set1_ps(donor2 == DONOR_AT ? 1.0f : 0.0f);
    }

    dpc  = ox->dpf[i];
    dpp3 = ox->dpf[i - 3];
    dplb = ox->dpf[i - min_intron - 3];

    mpv3  = esl_sse_rightshiftz_float(MMO_SP(dpp3, Q-1));
    dpv3  = esl_sse_rightshiftz_float(DMO_SP(dpp3, Q-1));
    ipv3  = esl_sse_rightshiftz_float(IMO_SP(dpp3, Q-1));
    ppv3  = esl_sse_rightshiftz_float(PMO_SP(dpp3, Q-1));
    tsv_m = esl_sse_rightshiftz_float(MMO_SP(dplb, Q-1));
    tsv_d = esl_sse_rightshiftz_float(DMO_SP(dplb, Q-1));
    dcv   = zerov;

    for (q = 0; q < Q; q++) {
      int is_last = (q == Q - 1);

      /* ---- M state ---- */
      msv = _mm_max_ps(_mm_mul_ps(mpv3, LTFV(TFV_MM, q)),
            _mm_max_ps(_mm_mul_ps(ipv3, LTFV(TFV_IM, q)),
            _mm_max_ps(_mm_mul_ps(dpv3, LTFV(TFV_DM, q)),
                       _mm_mul_ps(ppv3, TSC_P_v))));
      msv = _mm_mul_ps(msv, rfv_c0[q]);

      mpv3 = MMO_SP(dpp3, q);
      dpv3 = DMO_SP(dpp3, q);
      ipv3 = IMO_SP(dpp3, q);
      ppv3 = PMO_SP(dpp3, q);

      MMO_SP(dpc, q) = msv;
      DMO_SP(dpc, q) = dcv;
      dcv = _mm_mul_ps(msv, LTFV(TFV_MD, q));

      /* ---- I state (zero at k=M only, not entire last stripe) ---- */
      isv = _mm_max_ps(_mm_mul_ps(MMO_SP(dpp3, q), LTFV(TFV_MI, q)),
                       _mm_mul_ps(IMO_SP(dpp3, q), LTFV(TFV_II, q)));
      if (is_last) {
        union { __m128 v; float p[4]; } ii;
        ii.v = _mm_and_ps(isv, _mm_cmpgt_ps(rfv_c0[q], zerov));
        ii.p[r_M] = 0.0f;
        IMO_SP(dpc, q) = ii.v;
      } else {
        IMO_SP(dpc, q) = _mm_and_ps(isv, _mm_cmpgt_ps(rfv_c0[q], zerov));
      }

      /* ---- P state ---- */
      PMO_SP(dpc, q) = zerov;
      if (do_pstate) {
        pv = zerov;

        TMP_v = _mm_max_ps(_mm_mul_ps(OSS0(oss,q,p7S_GTAG), _mm_mul_ps(AG0v, sig_GTAG_v)),
                _mm_max_ps(_mm_mul_ps(OSS0(oss,q,p7S_GCAG), _mm_mul_ps(AG0v, sig_GCAG_v)),
                           _mm_mul_ps(OSS0(oss,q,p7S_ATAC), _mm_mul_ps(AC0v, sig_ATAC_v))));
        pv = _mm_max_ps(pv, _mm_mul_ps(TMP_v, rfv_c0[q]));

        for (nuc1 = 0; nuc1 < 4; nuc1++) {
          TMP_v = _mm_max_ps(_mm_mul_ps(OSS1(oss,q,p7S_GTAG,nuc1), _mm_mul_ps(AG1v, sig_GTAG_v)),
                  _mm_max_ps(_mm_mul_ps(OSS1(oss,q,p7S_GCAG,nuc1), _mm_mul_ps(AG1v, sig_GCAG_v)),
                             _mm_mul_ps(OSS1(oss,q,p7S_ATAC,nuc1), _mm_mul_ps(AC1v, sig_ATAC_v))));
          pv = _mm_max_ps(pv, _mm_mul_ps(TMP_v, rfv_c1[nuc1][q]));
        }

        for (nuc1 = 0; nuc1 < 4; nuc1++) {
          for (nuc2 = 0; nuc2 < 4; nuc2++) {
            TMP_v = _mm_max_ps(_mm_mul_ps(OSS2(oss,q,p7S_GTAG,nuc1,nuc2), _mm_mul_ps(AG2v, sig_GTAG_v)),
                    _mm_max_ps(_mm_mul_ps(OSS2(oss,q,p7S_GCAG,nuc1,nuc2), _mm_mul_ps(AG2v, sig_GCAG_v)),
                               _mm_mul_ps(OSS2(oss,q,p7S_ATAC,nuc1,nuc2), _mm_mul_ps(AC2v, sig_ATAC_v))));
            pv = _mm_max_ps(pv, _mm_mul_ps(TMP_v, rfv_c2[nuc1*4+nuc2][q]));
          }
        }
        PMO_SP(dpc, q) = pv;
      }

      /* ---- Donor lookback: update OSS at position k from (lb, k-1) ---- */
      if (do_donor) {
        TMP_v = _mm_max_ps(tsv_m, tsv_d);

        OSS0(oss,q,p7S_GTAG) = _mm_max_ps(OSS0(oss,q,p7S_GTAG), _mm_mul_ps(TMP_v, GT0v));
        OSS0(oss,q,p7S_GCAG) = _mm_max_ps(OSS0(oss,q,p7S_GCAG), _mm_mul_ps(TMP_v, GC0v));
        OSS0(oss,q,p7S_ATAC) = _mm_max_ps(OSS0(oss,q,p7S_ATAC), _mm_mul_ps(TMP_v, AT0v));

        OSS1(oss,q,p7S_GTAG,r) = _mm_max_ps(OSS1(oss,q,p7S_GTAG,r), _mm_mul_ps(TMP_v, GT1v));
        OSS1(oss,q,p7S_GCAG,r) = _mm_max_ps(OSS1(oss,q,p7S_GCAG,r), _mm_mul_ps(TMP_v, GC1v));
        OSS1(oss,q,p7S_ATAC,r) = _mm_max_ps(OSS1(oss,q,p7S_ATAC,r), _mm_mul_ps(TMP_v, AT1v));

        OSS2(oss,q,p7S_GTAG,r,s) = _mm_max_ps(OSS2(oss,q,p7S_GTAG,r,s), _mm_mul_ps(TMP_v, GT2v));
        OSS2(oss,q,p7S_GCAG,r,s) = _mm_max_ps(OSS2(oss,q,p7S_GCAG,r,s), _mm_mul_ps(TMP_v, GC2v));
        OSS2(oss,q,p7S_ATAC,r,s) = _mm_max_ps(OSS2(oss,q,p7S_ATAC,r,s), _mm_mul_ps(TMP_v, AT2v));
      }

      tsv_m = MMO_SP(dplb, q);
      tsv_d = DMO_SP(dplb, q);
    }

    /* D-state serialization */
    dcv = esl_sse_rightshiftz_float(dcv);
    DMO_SP(dpc, 0) = zerov;
    for (q = 0; q < Q; q++) {
      DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
      dcv = _mm_mul_ps(DMO_SP(dpc, q), LTFV(TFV_DD, q));
    }
    for (j = 1; j < 4; j++) {
      dcv = esl_sse_rightshiftz_float(dcv);
      for (q = 0; q < Q; q++) {
        DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
        dcv = _mm_mul_ps(DMO_SP(dpc, q), LTFV(TFV_DD, q));
      }
    }

    /* Semi-global exit: E(i) and C(i) */
    {
      __m128 emax_v = zerov;
      union { __m128 v; float p[4]; } em;
      float xE, xC;
      for (q = 0; q < Q; q++)
        emax_v = _mm_max_ps(emax_v, _mm_and_ps(_mm_max_ps(MMO_SP(dpc, q), DMO_SP(dpc, q)), valid_mask[q]));
      em.v = emax_v;
      xE   = ESL_MAX(ESL_MAX(em.p[0], em.p[1]), ESL_MAX(em.p[2], em.p[3]));
      xC   = ESL_MAX(ox->xmx[(i - 3) * p7X_NXCELLS + p7X_C] * CC_loop, xE * EC_move);
      ox->xmx[i * p7X_NXCELLS + p7X_E] = xE;
      ox->xmx[i * p7X_NXCELLS + p7X_C] = xC;
    }
  }  /* end main DP loop */

  ox->M = M;
  ox->L = L;

  free(mask_raw);
  free(rfv_raw);
  free(tfv_raw);
  return eslOK;

  ERROR:
    if (mask_raw) free(mask_raw);
    if (rfv_raw)  free(rfv_raw);
    if (tfv_raw)  free(tfv_raw);
    return status;
}
/*--- end, p7_ospliceviterbi_TranslatedSemiGlobalExtendDown() ---*/


/*****************************************************************
 * 4. p7_ospliceviterbi_TranslatedSemiGlobalExtendUp()
 *****************************************************************/

/* Function:  p7_ospliceviterbi_TranslatedSemiGlobalExtendUp()
 * Synopsis:  SSE-accelerated semi-global (semi-global entry, global exit) translated spliced Viterbi.
 *
 * Purpose:   SSE equivalent of p7_spliceviterbi_TranslatedSemiGlobalExtendUp().
 *            Differs from p7_ospliceviterbi_TranslatedGlobal() in that entry is
 *            semi-global: at every row i the N state advances (N(i) = N(i-3)*NN_loop)
 *            and B(i) = N(i)*NB_move is broadcast into the M recurrence for all k, so
 *            the model can enter at any k.  Exit is global: E(L) = max(M(L,M), D(L,M)),
 *            C(L) = E(L) * EC_move.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile not codon_lengths==1, <eslEMEM> on alloc.
 */
int
p7_ospliceviterbi_TranslatedSemiGlobalExtendUp(SPLICE_PIPELINE *pli, OSPLICE_SCORES *oss,
                                                const ESL_DSQ *sub_dsq,
                                                const P7_FS_OPROFILE *om_fs,
                                                P7_OMX *ox,
                                                int i_start, int i_end,
                                                int k_start, int k_end)
{
  int      M          = k_end - k_start + 1;
  int      L          = i_end - i_start + 1;
  int      Q          = p7O_NQF(M);
  int      r_M        = (M - 1) / Q;
  int      min_intron = pli->min_intron;

  __m128  *local_rfv  = NULL;
  __m128  *local_tfv  = NULL;
  void    *rfv_raw    = NULL;
  void    *tfv_raw    = NULL;

  __m128  *dpc, *dpp3, *dplb;

  __m128   zerov     = _mm_setzero_ps();
  __m128   TSC_P_v   = _mm_set1_ps(TSC_P_PROB);

  __m128   mpv3, dpv3, ipv3, ppv3;
  __m128   tsv_m, tsv_d;
  __m128   msv, isv, dcv;
  __m128   sig_GTAG_v, sig_GCAG_v, sig_ATAC_v;

  int      acc0, acc1, acc2;           /* ACCEPT_AG, ACCEPT_AC, or -1 */
  int      donor0, donor1, donor2;       /* DONOR_GT/GC/AT, or -1       */
  int      do_pstate, do_donor;
  __m128   AG0v, AG1v, AG2v, AC0v, AC1v, AC2v;
  __m128   GT0v, GC0v, AT0v, GT1v, GC1v, AT1v, GT2v, GC2v, AT2v;

  int      C0, C1[4], C2[16];
  int      v, w, x;
  int      r, s, t, u;
  int      nuc1, nuc2;

  float    NN_loop, NB_move, EC_move;
  int      i, q, j, loop_end, sub_i;
  int      status;

  if (om_fs->codon_lengths != 1)
    ESL_EXCEPTION(eslEINVAL, "profile not allocated for 1 codon length");

  local_rfv = build_local_rfv(om_fs, k_start, M, Q, &rfv_raw);
  local_tfv = build_local_tfv(om_fs, k_start, M, Q, &tfv_raw);
  if (!local_rfv || !local_tfv) { status = eslEMEM; goto ERROR; }

  sig_GTAG_v = _mm_set1_ps(oss->signal_scores[p7S_GTAG]);
  sig_GCAG_v = _mm_set1_ps(oss->signal_scores[p7S_GCAG]);
  sig_ATAC_v = _mm_set1_ps(oss->signal_scores[p7S_ATAC]);

  NN_loop = om_fs->xf[p7O_N][p7O_LOOP];
  NB_move = om_fs->xf[p7O_N][p7O_MOVE];
  EC_move = om_fs->xf[p7O_E][p7O_MOVE];

  for (q = 0; q < Q; q++)
    for (j = 0; j < SIGNAL_MEM_SIZE; j++)
      oss->oscore[q][j] = zerov;

  /* Zero rows 0, 1, 2: DP cells.
   * N and B for semi-global entry: N(0..2)=1.0 (exp(0)), B(0..2)=NB_move. */
  for (i = 0; i <= 2; i++) {
    dpc = ox->dpf[i];
    for (q = 0; q < Q; q++)
      MMO_SP(dpc,q) = DMO_SP(dpc,q) = IMO_SP(dpc,q) = PMO_SP(dpc,q) = zerov;
    ox->xmx[i * p7X_NXCELLS + p7X_N] = 1.0f;
    ox->xmx[i * p7X_NXCELLS + p7X_B] = NB_move;
    ox->xmx[i * p7X_NXCELLS + p7X_E] = 0.0f;
    ox->xmx[i * p7X_NXCELLS + p7X_C] = 0.0f;
  }

  acc0 = acc1 = acc2 = -1;
  v = w = x = -1;

  for (i = 1; i <= 2; i++) {
    w = x;
    sub_i = i_start + i - 1;
    x = (sub_dsq[sub_i] < 4) ? (int)sub_dsq[sub_i] : p7P_MAXCODONS1;
  }

  /* ---------------------------------------------------------------
   * Early rows 3..min_intron+2: M/I/D only, no donor lookback (P=0).
   * N/B advance by codon stride; B(i-3) broadcast into M recurrence.
   * --------------------------------------------------------------- */
  loop_end = ESL_MIN(L, min_intron + 2);

  for (i = 3; i <= loop_end; i++) {
    v = w;  w = x;
    sub_i = i_start + i - 1;
    x = (sub_dsq[sub_i] < 4) ? (int)sub_dsq[sub_i] : p7P_MAXCODONS1;

    C0 = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);
    __m128 *rfv_c0 = local_rfv + C0 * Q;

    /* Acceptor window update */
    acc0 = acc1;  acc1 = acc2;
    if (v < 0 || v >= 4 || w < 0 || w >= 4)      acc2 = -1;
    else if (SIGNAL(v,w) == ACCEPT_AG)             acc2 = ACCEPT_AG;
    else if (SIGNAL(v,w) == ACCEPT_AC)             acc2 = ACCEPT_AC;
    else                                           acc2 = -1;

    /* N/B: codon-stride recurrence */
    {
      float xN = ox->xmx[(i - 3) * p7X_NXCELLS + p7X_N] * NN_loop;
      ox->xmx[i * p7X_NXCELLS + p7X_N] = xN;
      ox->xmx[i * p7X_NXCELLS + p7X_B] = xN * NB_move;
    }

    dpc  = ox->dpf[i];
    dpp3 = ox->dpf[i - 3];

    /* B(i-3) broadcast for semi-global entry */
    __m128 xBv = _mm_set1_ps(ox->xmx[(i - 3) * p7X_NXCELLS + p7X_B]);

    mpv3 = esl_sse_rightshiftz_float(MMO_SP(dpp3, Q-1));
    dpv3 = esl_sse_rightshiftz_float(DMO_SP(dpp3, Q-1));
    ipv3 = esl_sse_rightshiftz_float(IMO_SP(dpp3, Q-1));
    ppv3 = esl_sse_rightshiftz_float(PMO_SP(dpp3, Q-1));
    dcv  = zerov;

    for (q = 0; q < Q; q++) {
      int is_last = (q == Q - 1);

      msv = _mm_max_ps(_mm_mul_ps(mpv3, LTFV(TFV_MM, q)),
            _mm_max_ps(_mm_mul_ps(ipv3, LTFV(TFV_IM, q)),
            _mm_max_ps(_mm_mul_ps(dpv3, LTFV(TFV_DM, q)),
            _mm_max_ps(xBv,
                       _mm_mul_ps(ppv3, TSC_P_v)))));
      msv = _mm_mul_ps(msv, rfv_c0[q]);

      mpv3 = MMO_SP(dpp3, q);
      dpv3 = DMO_SP(dpp3, q);
      ipv3 = IMO_SP(dpp3, q);
      ppv3 = PMO_SP(dpp3, q);

      MMO_SP(dpc, q) = msv;
      DMO_SP(dpc, q) = dcv;
      dcv             = _mm_mul_ps(msv, LTFV(TFV_MD, q));

      isv = _mm_max_ps(_mm_mul_ps(MMO_SP(dpp3, q), LTFV(TFV_MI, q)),
                       _mm_mul_ps(IMO_SP(dpp3, q), LTFV(TFV_II, q)));
      if (is_last) {
        union { __m128 v; float p[4]; } ii;
        ii.v = _mm_and_ps(isv, _mm_cmpgt_ps(rfv_c0[q], zerov));
        ii.p[r_M] = 0.0f;
        IMO_SP(dpc, q) = ii.v;
      } else {
        IMO_SP(dpc, q) = _mm_and_ps(isv, _mm_cmpgt_ps(rfv_c0[q], zerov));
      }

      PMO_SP(dpc, q) = zerov;
    }

    /* D-state serialization */
    dcv = esl_sse_rightshiftz_float(dcv);
    DMO_SP(dpc, 0) = zerov;
    for (q = 0; q < Q; q++) {
      DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
      dcv = _mm_mul_ps(DMO_SP(dpc, q), LTFV(TFV_DD, q));
    }
    for (j = 1; j < 4; j++) {
      dcv = esl_sse_rightshiftz_float(dcv);
      for (q = 0; q < Q; q++) {
        DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
        dcv = _mm_mul_ps(DMO_SP(dpc, q), LTFV(TFV_DD, q));
      }
    }
  }  /* end early rows */


  /* ---------------------------------------------------------------
   * Initialize donor-lookback nucleotide window for the main loop.
   * --------------------------------------------------------------- */
  s = (int)sub_dsq[i_start];
  t = (int)sub_dsq[i_start + 1];
  u = (int)sub_dsq[i_start + 2];
  r = 0;

  /* ---------------------------------------------------------------
   * Main DP: rows min_intron+3 .. L.
   * Full M/I/D/P + donor lookback OSS write.
   * N/B advance by codon stride; B(i-3) broadcast into M recurrence.
   * --------------------------------------------------------------- */
  for (i = min_intron + 3; i <= L; i++) {
    __m128 pv, TMP_v;

    sub_i = i_start + i - 1;

    r = s;  s = t;  t = u;
    u = (sub_dsq[sub_i - min_intron + 1] < 4)
        ? (int)sub_dsq[sub_i - min_intron + 1]
        : p7P_MAXCODONS1;

    v = w;  w = x;
    x = (sub_dsq[sub_i] < 4) ? (int)sub_dsq[sub_i] : p7P_MAXCODONS1;

    C0 = p7P_MINIDX(p7P_CODON3_FS1(v, w, x), p7P_DEGEN1_C);

    __m128 *rfv_c0 = local_rfv + C0 * Q;
    __m128 *rfv_c1[4];
    __m128 *rfv_c2[16];

    acc0 = acc1;  acc1 = acc2;
    if (v < 0 || v >= 4 || w < 0 || w >= 4)      acc2 = -1;
    else if (SIGNAL(v,w) == ACCEPT_AG)             acc2 = ACCEPT_AG;
    else if (SIGNAL(v,w) == ACCEPT_AC)             acc2 = ACCEPT_AC;
    else                                           acc2 = -1;

    donor0 = donor1 = donor2 = -1;
    if (!(r < 0 || r >= 4 || s < 0 || s >= 4)) {
      if      (SIGNAL(r,s) == DONOR_GT) donor0 = DONOR_GT;
      else if (SIGNAL(r,s) == DONOR_GC) donor0 = DONOR_GC;
      else if (SIGNAL(r,s) == DONOR_AT) donor0 = DONOR_AT;
      if (!(t < 0 || t >= 4)) {
        if      (SIGNAL(s,t) == DONOR_GT) donor1 = DONOR_GT;
        else if (SIGNAL(s,t) == DONOR_GC) donor1 = DONOR_GC;
        else if (SIGNAL(s,t) == DONOR_AT) donor1 = DONOR_AT;
        if (!(u < 0 || u >= 4)) {
          if      (SIGNAL(t,u) == DONOR_GT) donor2 = DONOR_GT;
          else if (SIGNAL(t,u) == DONOR_GC) donor2 = DONOR_GC;
          else if (SIGNAL(t,u) == DONOR_AT) donor2 = DONOR_AT;
        }
      }
    }

    do_pstate = (acc0 >= 0 || acc1 >= 0 || acc2 >= 0);
    do_donor  = (donor0 >= 0 || donor1 >= 0 || donor2 >= 0);

    if (do_pstate) {
      for (nuc1 = 0; nuc1 < 4; nuc1++) {
        C1[nuc1] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, w, x), p7P_DEGEN1_C);
        rfv_c1[nuc1] = local_rfv + C1[nuc1] * Q;
        for (nuc2 = 0; nuc2 < 4; nuc2++) {
          C2[nuc1*4+nuc2] = p7P_MINIDX(p7P_CODON3_FS1(nuc1, nuc2, x), p7P_DEGEN1_C);
          rfv_c2[nuc1*4+nuc2] = local_rfv + C2[nuc1*4+nuc2] * Q;
        }
      }
      AG0v = _mm_set1_ps(acc0 == ACCEPT_AG ? 1.0f : 0.0f);
      AG1v = _mm_set1_ps(acc1 == ACCEPT_AG ? 1.0f : 0.0f);
      AG2v = _mm_set1_ps(acc2 == ACCEPT_AG ? 1.0f : 0.0f);
      AC0v = _mm_set1_ps(acc0 == ACCEPT_AC ? 1.0f : 0.0f);
      AC1v = _mm_set1_ps(acc1 == ACCEPT_AC ? 1.0f : 0.0f);
      AC2v = _mm_set1_ps(acc2 == ACCEPT_AC ? 1.0f : 0.0f);
    }
    if (do_donor) {
      GT0v = _mm_set1_ps(donor0 == DONOR_GT ? 1.0f : 0.0f);
      GC0v = _mm_set1_ps(donor0 == DONOR_GC ? 1.0f : 0.0f);
      AT0v = _mm_set1_ps(donor0 == DONOR_AT ? 1.0f : 0.0f);
      GT1v = _mm_set1_ps(donor1 == DONOR_GT ? 1.0f : 0.0f);
      GC1v = _mm_set1_ps(donor1 == DONOR_GC ? 1.0f : 0.0f);
      AT1v = _mm_set1_ps(donor1 == DONOR_AT ? 1.0f : 0.0f);
      GT2v = _mm_set1_ps(donor2 == DONOR_GT ? 1.0f : 0.0f);
      GC2v = _mm_set1_ps(donor2 == DONOR_GC ? 1.0f : 0.0f);
      AT2v = _mm_set1_ps(donor2 == DONOR_AT ? 1.0f : 0.0f);
    }

    /* N/B: codon-stride recurrence */
    {
      float xN = ox->xmx[(i - 3) * p7X_NXCELLS + p7X_N] * NN_loop;
      ox->xmx[i * p7X_NXCELLS + p7X_N] = xN;
      ox->xmx[i * p7X_NXCELLS + p7X_B] = xN * NB_move;
    }

    dpc  = ox->dpf[i];
    dpp3 = ox->dpf[i - 3];
    dplb = ox->dpf[i - min_intron - 3];

    /* B(i-3) broadcast for semi-global entry */
    __m128 xBv = _mm_set1_ps(ox->xmx[(i - 3) * p7X_NXCELLS + p7X_B]);

    mpv3  = esl_sse_rightshiftz_float(MMO_SP(dpp3, Q-1));
    dpv3  = esl_sse_rightshiftz_float(DMO_SP(dpp3, Q-1));
    ipv3  = esl_sse_rightshiftz_float(IMO_SP(dpp3, Q-1));
    ppv3  = esl_sse_rightshiftz_float(PMO_SP(dpp3, Q-1));
    tsv_m = esl_sse_rightshiftz_float(MMO_SP(dplb, Q-1));
    tsv_d = esl_sse_rightshiftz_float(DMO_SP(dplb, Q-1));
    dcv   = zerov;

    for (q = 0; q < Q; q++) {
      int is_last = (q == Q - 1);

      /* ---- M state (semi-global entry: B(i-3) broadcast to all k) ---- */
      msv = _mm_max_ps(_mm_mul_ps(mpv3, LTFV(TFV_MM, q)),
            _mm_max_ps(_mm_mul_ps(ipv3, LTFV(TFV_IM, q)),
            _mm_max_ps(_mm_mul_ps(dpv3, LTFV(TFV_DM, q)),
            _mm_max_ps(xBv,
                       _mm_mul_ps(ppv3, TSC_P_v)))));
      msv = _mm_mul_ps(msv, rfv_c0[q]);

      mpv3 = MMO_SP(dpp3, q);
      dpv3 = DMO_SP(dpp3, q);
      ipv3 = IMO_SP(dpp3, q);
      ppv3 = PMO_SP(dpp3, q);

      MMO_SP(dpc, q) = msv;
      DMO_SP(dpc, q) = dcv;
      dcv = _mm_mul_ps(msv, LTFV(TFV_MD, q));

      /* ---- I state (zero at k=M only) ---- */
      isv = _mm_max_ps(_mm_mul_ps(MMO_SP(dpp3, q), LTFV(TFV_MI, q)),
                       _mm_mul_ps(IMO_SP(dpp3, q), LTFV(TFV_II, q)));
      if (is_last) {
        union { __m128 v; float p[4]; } ii;
        ii.v = _mm_and_ps(isv, _mm_cmpgt_ps(rfv_c0[q], zerov));
        ii.p[r_M] = 0.0f;
        IMO_SP(dpc, q) = ii.v;
      } else {
        IMO_SP(dpc, q) = _mm_and_ps(isv, _mm_cmpgt_ps(rfv_c0[q], zerov));
      }

      /* ---- P state ---- */
      PMO_SP(dpc, q) = zerov;
      if (do_pstate) {
        pv = zerov;

        TMP_v = _mm_max_ps(_mm_mul_ps(OSS0(oss,q,p7S_GTAG), _mm_mul_ps(AG0v, sig_GTAG_v)),
                _mm_max_ps(_mm_mul_ps(OSS0(oss,q,p7S_GCAG), _mm_mul_ps(AG0v, sig_GCAG_v)),
                           _mm_mul_ps(OSS0(oss,q,p7S_ATAC), _mm_mul_ps(AC0v, sig_ATAC_v))));
        pv = _mm_max_ps(pv, _mm_mul_ps(TMP_v, rfv_c0[q]));

        for (nuc1 = 0; nuc1 < 4; nuc1++) {
          TMP_v = _mm_max_ps(_mm_mul_ps(OSS1(oss,q,p7S_GTAG,nuc1), _mm_mul_ps(AG1v, sig_GTAG_v)),
                  _mm_max_ps(_mm_mul_ps(OSS1(oss,q,p7S_GCAG,nuc1), _mm_mul_ps(AG1v, sig_GCAG_v)),
                             _mm_mul_ps(OSS1(oss,q,p7S_ATAC,nuc1), _mm_mul_ps(AC1v, sig_ATAC_v))));
          pv = _mm_max_ps(pv, _mm_mul_ps(TMP_v, rfv_c1[nuc1][q]));
        }

        for (nuc1 = 0; nuc1 < 4; nuc1++) {
          for (nuc2 = 0; nuc2 < 4; nuc2++) {
            TMP_v = _mm_max_ps(_mm_mul_ps(OSS2(oss,q,p7S_GTAG,nuc1,nuc2), _mm_mul_ps(AG2v, sig_GTAG_v)),
                    _mm_max_ps(_mm_mul_ps(OSS2(oss,q,p7S_GCAG,nuc1,nuc2), _mm_mul_ps(AG2v, sig_GCAG_v)),
                               _mm_mul_ps(OSS2(oss,q,p7S_ATAC,nuc1,nuc2), _mm_mul_ps(AC2v, sig_ATAC_v))));
            pv = _mm_max_ps(pv, _mm_mul_ps(TMP_v, rfv_c2[nuc1*4+nuc2][q]));
          }
        }
        PMO_SP(dpc, q) = pv;
      }

      /* ---- Donor lookback: update OSS at position k from (lb, k-1) ---- */
      if (do_donor) {
        TMP_v = _mm_max_ps(tsv_m, tsv_d);

        OSS0(oss,q,p7S_GTAG) = _mm_max_ps(OSS0(oss,q,p7S_GTAG), _mm_mul_ps(TMP_v, GT0v));
        OSS0(oss,q,p7S_GCAG) = _mm_max_ps(OSS0(oss,q,p7S_GCAG), _mm_mul_ps(TMP_v, GC0v));
        OSS0(oss,q,p7S_ATAC) = _mm_max_ps(OSS0(oss,q,p7S_ATAC), _mm_mul_ps(TMP_v, AT0v));

        OSS1(oss,q,p7S_GTAG,r) = _mm_max_ps(OSS1(oss,q,p7S_GTAG,r), _mm_mul_ps(TMP_v, GT1v));
        OSS1(oss,q,p7S_GCAG,r) = _mm_max_ps(OSS1(oss,q,p7S_GCAG,r), _mm_mul_ps(TMP_v, GC1v));
        OSS1(oss,q,p7S_ATAC,r) = _mm_max_ps(OSS1(oss,q,p7S_ATAC,r), _mm_mul_ps(TMP_v, AT1v));

        OSS2(oss,q,p7S_GTAG,r,s) = _mm_max_ps(OSS2(oss,q,p7S_GTAG,r,s), _mm_mul_ps(TMP_v, GT2v));
        OSS2(oss,q,p7S_GCAG,r,s) = _mm_max_ps(OSS2(oss,q,p7S_GCAG,r,s), _mm_mul_ps(TMP_v, GC2v));
        OSS2(oss,q,p7S_ATAC,r,s) = _mm_max_ps(OSS2(oss,q,p7S_ATAC,r,s), _mm_mul_ps(TMP_v, AT2v));
      }

      tsv_m = MMO_SP(dplb, q);
      tsv_d = DMO_SP(dplb, q);
    }

    /* D-state serialization */
    dcv = esl_sse_rightshiftz_float(dcv);
    DMO_SP(dpc, 0) = zerov;
    for (q = 0; q < Q; q++) {
      DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
      dcv = _mm_mul_ps(DMO_SP(dpc, q), LTFV(TFV_DD, q));
    }
    for (j = 1; j < 4; j++) {
      dcv = esl_sse_rightshiftz_float(dcv);
      for (q = 0; q < Q; q++) {
        DMO_SP(dpc, q) = _mm_max_ps(dcv, DMO_SP(dpc, q));
        dcv = _mm_mul_ps(DMO_SP(dpc, q), LTFV(TFV_DD, q));
      }
    }
  }  /* end main DP loop */

  /* ---------------------------------------------------------------
   * Global exit: E(L) = max(M(L,M), D(L,M)); C(L) = E(L) * EC_move.
   * --------------------------------------------------------------- */
  {
    int   q_M = (M - 1) % Q;
    float e_M, e_D, xE;
    union { __m128 v; float p[4]; } u;

    dpc = ox->dpf[L];
    u.v = MMO_SP(dpc, q_M);  e_M = u.p[r_M];
    u.v = DMO_SP(dpc, q_M);  e_D = u.p[r_M];
    xE  = ESL_MAX(e_M, e_D);

    ox->xmx[L * p7X_NXCELLS + p7X_E] = xE;
    ox->xmx[L * p7X_NXCELLS + p7X_C] = xE * EC_move;
  }

  ox->M = M;
  ox->L = L;

  free(rfv_raw);
  free(tfv_raw);
  return eslOK;

  ERROR:
    if (rfv_raw) free(rfv_raw);
    if (tfv_raw) free(tfv_raw);
    return status;
}
/*---- end, p7_ospliceviterbi_TranslatedSemiGlobalExtendUp() ----*/


/*****************************************************************
 * 5. Benchmark driver.
 *****************************************************************/
#ifdef p7SPLICED_VITERBI_BENCHMARK
/*
   gcc -g -O3 -Wall -msse2 -std=gnu99 -o spliced_viterbi_benchmark -I.. -L.. -I../../easel -L../../easel -Dp7SPLICED_VITERBI_BENCHMARK spliced_viterbi.c -lhmmer -leasel -lm
   ./spliced_viterbi_benchmark <hmmfile>
 */
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

static ESL_OPTIONS benchmark_options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,   "1200", NULL, "n>0", NULL,  NULL, NULL, "length of random target DNA seqs",               0 },
  { "-N",        eslARG_INT,    "100", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "-G",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark TranslatedGlobal",                0 },
  { "-D",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark TranslatedSemiGlobalExtendDown",  0 },
  { "-U",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark TranslatedSemiGlobalExtendUp",    0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char benchmark_usage[]  = "[-options] <hmmfile>";
static char benchmark_banner[] = "benchmark driver for SSE spliced Viterbi DP algorithms";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(benchmark_options, 1, argc, argv, benchmark_banner, benchmark_usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA   = NULL;
  ESL_ALPHABET   *abcDNA  = esl_alphabet_Create(eslDNA);
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bgAA    = NULL;
  P7_BG          *bgDNA   = p7_bg_Create(abcDNA);
  P7_FS_PROFILE  *gm_tr   = NULL;
  P7_FS_OPROFILE *om_fs   = NULL;
  ESL_GENCODE    *gcode   = NULL;
  SPLICE_PIPELINE *pli    = NULL;
  P7_OMX         *ox      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  int             do_G    = (esl_opt_GetBoolean(go, "-G") || (!esl_opt_GetBoolean(go, "-D") && !esl_opt_GetBoolean(go, "-U")));
  int             do_D    = (esl_opt_GetBoolean(go, "-D") || (!esl_opt_GetBoolean(go, "-G") && !esl_opt_GetBoolean(go, "-U")));
  int             do_U    = (esl_opt_GetBoolean(go, "-U") || (!esl_opt_GetBoolean(go, "-G") && !esl_opt_GetBoolean(go, "-D")));
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L + 2));
  int             i;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)          != eslOK) p7_Fail("Failed to read HMM");

  gcode  = esl_gencode_Create(abcDNA, abcAA);
  bgAA   = p7_bg_Create(abcAA);
  gm_tr  = p7_profile_fs_Create(hmm->M, abcAA, 1);
  om_fs  = p7_fs_oprofile_Create(hmm->M, abcAA, 1);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, L / 3, p7_UNILOCAL);
  p7_fs_ReconfigLength(gm_tr, L / 3);
  p7_fs_oprofile_Convert(gm_tr, om_fs);
  p7_fs_oprofile_ReconfigLength(om_fs, L / 3);

  pli = p7_splicepipeline_Create(NULL, hmm->M, L);
  ox  = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS_SP);
  p7_omx_GrowTo_dpf      (ox,            hmm->M, L, L);
  p7_osplicescores_GrowTo(pli->osplice_scores, hmm->M);

  /* Baseline: time to generate sequences alone */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
      if (do_G) p7_ospliceviterbi_TranslatedGlobal             (pli, pli->osplice_scores, dsq, om_fs, ox, 1, L, 1, hmm->M);
      if (do_D) p7_ospliceviterbi_TranslatedSemiGlobalExtendDown(pli, pli->osplice_scores, dsq, om_fs, ox, 1, L, 1, hmm->M);
      if (do_U) p7_ospliceviterbi_TranslatedSemiGlobalExtendUp  (pli, pli->osplice_scores, dsq, om_fs, ox, 1, L, 1, hmm->M);
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) hmm->M * 1e-6 / bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   hmm->M);
  printf("# L    = %d\n",   L);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_omx_Destroy(ox);
  p7_splicepipeline_Destroy(pli);
  p7_fs_oprofile_Destroy(om_fs);
  p7_profile_fs_Destroy(gm_tr);
  p7_bg_Destroy(bgAA);
  p7_bg_Destroy(bgDNA);
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
#endif /*p7SPLICED_VITERBI_BENCHMARK*/
/*---------------- end, benchmark driver ----------------*/


/*****************************************************************
 * 6. Unit tests.
 *****************************************************************/
#ifdef p7SPLICED_VITERBI_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_sviterbi():
 * Compare all three SSE spliced-Viterbi implementations against their scalar
 * (log-space) counterparts using the same profile-emitted sequence each round.
 * Sequence generation, profile setup, and matrix allocation are done once.
 * Per-sequence loop runs all three scalar/SSE pairs and checks agreement.
 */
static void
utest_sviterbi(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
          ESL_GENCODE *gcode, P7_BG *bgAA, P7_CODONTABLE *codon_table,
          int M, int N)
{
  char           *msg_g  = "spliced viterbi global unit test failed";
  char           *msg_d  = "spliced viterbi extend_down unit test failed";
  char           *msg_u  = "spliced viterbi extend_up unit test failed";
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_tr  = p7_profile_fs_Create(M, abcAA, 1);
  P7_FS_OPROFILE *om_fs  = p7_fs_oprofile_Create(M, abcAA, 1);
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  P7_TRACE       *tr     = p7_trace_Create();
  ESL_DSQ        *dsq    = NULL;
  SPLICE_PIPELINE *pli   = NULL;
  P7_OMX         *ox     = NULL;
  int             L_dna, L_amino;
  int             i, j;
  float           scalar_E, sse_E;
  float           scalar_C, sse_C;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig   (hmm, bgAA, gm,    M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, M, p7_UNILOCAL);
  p7_fs_oprofile_Convert(gm_tr, om_fs);

  pli = p7_splicepipeline_Create(NULL, M, M * 3);
  ox  = p7_omx_Create_dpf(M, M * 3, M * 3, p7X_NSCELLS_SP);

  while (N--)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, tr);
      L_amino = sq->n;
      L_dna   = L_amino * 3;

      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (L_dna + 2))) == NULL) esl_fatal("malloc failed");
      dsq[0] = dsq[L_dna + 1] = eslDSQ_SENTINEL;

      j = 1;
      for (i = 1; i <= sq->n; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }

      p7_fs_ReconfigLength         (gm_tr, L_amino);
      p7_fs_oprofile_ReconfigLength(om_fs,  L_amino);

      p7_gmx_GrowTo          (pli->vit, M, L_dna, L_dna);
      p7_omx_GrowTo_dpf      (ox,       M, L_dna, L_dna);
      p7_splicescores_GrowTo (pli->splice_scores,  M);
      p7_osplicescores_GrowTo(pli->osplice_scores, M);

      /* --- Global --- */
      p7_spliceviterbi_TranslatedGlobal(pli, dsq, gm_tr, pli->vit, 1, L_dna, 1, M);
      scalar_E = pli->vit->xmx[L_dna * p7G_NXCELLS + p7G_E];

      p7_ospliceviterbi_TranslatedGlobal(pli, pli->osplice_scores, dsq, om_fs, ox, 1, L_dna, 1, M);
      sse_E = ox->xmx[L_dna * p7X_NXCELLS + p7X_E];

      if (scalar_E <= -1e30f) {
        if (sse_E > 0.0f) esl_fatal(msg_g);
      } else {
        if (fabsf(expf(scalar_E) - sse_E) > 0.001f) esl_fatal(msg_g);
      }

      /* --- ExtendDown --- */
      p7_spliceviterbi_TranslatedSemiGlobalExtendDown(pli, dsq, gm_tr, pli->vit, 1, L_dna, 1, M);
      scalar_C = pli->vit->xmx[L_dna * p7G_NXCELLS + p7G_C];

      p7_ospliceviterbi_TranslatedSemiGlobalExtendDown(pli, pli->osplice_scores, dsq, om_fs, ox, 1, L_dna, 1, M);
      sse_C = ox->xmx[L_dna * p7X_NXCELLS + p7X_C];

      if (scalar_C <= -1e30f) {
        if (sse_C > 0.0f) esl_fatal(msg_d);
      } else {
        if (fabsf(expf(scalar_C) - sse_C) > 0.001f) esl_fatal(msg_d);
      }

      /* --- ExtendUp --- */
      p7_spliceviterbi_TranslatedSemiGlobalExtendUp(pli, dsq, gm_tr, pli->vit, 1, L_dna, 1, M);
      scalar_C = pli->vit->xmx[L_dna * p7G_NXCELLS + p7G_C];

      p7_ospliceviterbi_TranslatedSemiGlobalExtendUp(pli, pli->osplice_scores, dsq, om_fs, ox, 1, L_dna, 1, M);
      sse_C = ox->xmx[L_dna * p7X_NXCELLS + p7X_C];

      if (scalar_C <= -1e30f) {
        if (sse_C > 0.0f) esl_fatal(msg_u);
      } else {
        if (fabsf(expf(scalar_C) - sse_C) > 0.001f) esl_fatal(msg_u);
      }
    }

  if (dsq != NULL) free(dsq);
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_hmm_Destroy(hmm);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_tr);
  p7_fs_oprofile_Destroy(om_fs);
  p7_splicepipeline_Destroy(pli);
  p7_omx_Destroy(ox);
}
#endif /*p7SPLICED_VITERBI_TESTDRIVE*/
/*---------------------- end, unit tests ------------------------*/


/*****************************************************************
 * 7. Test driver.
 *****************************************************************/
#ifdef p7SPLICED_VITERBI_TESTDRIVE
/*
   gcc -g -Wall -msse2 -std=gnu99 -o spliced_viterbi_utest -I.. -L.. -I../../easel -L../../easel -Dp7SPLICED_VITERBI_TESTDRIVE spliced_viterbi.c -lhmmer -leasel -lm
   ./spliced_viterbi_utest
 */
#include "esl_gencode.h"
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M",        eslARG_INT,     "20", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,     "20", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE spliced Viterbi";

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

  utest_sviterbi(r, abcAA, abcDNA, gcode, bgAA, ct, M, N);

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
#endif /*p7SPLICED_VITERBI_TESTDRIVE*/
/*--------------------- end, test driver ------------------------*/
