/* SSE-accelerated Spliced Viterbi algorithms; full matrix.
 *
 * Contents:
 *   1. p7_Viterbi_SplicedGlobal()
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
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

      /* I(i,k) = max(M(i-3,k)*MI, I(i-3,k)*II).  Uses same-k values from dpp3. */
      sv             = _mm_add_ps(mpv, *tp); tp++;      /* MI */
      IMO_SP(dpc, q) = _mm_max_ps(sv, _mm_add_ps(ipv, *tp)); tp++;  /* II */

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
 * 3. Benchmark driver.
 *****************************************************************/



/*----------------- end, benchmark driver -----------------------*/



/*****************************************************************
 * 4. Unit tests.
 *****************************************************************/
#ifdef p7VITERBI_SP_TESTDRIVE
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

  while (N--)
    {
	  p7_emit_SimpleConsensus(hmm, sq);
      L_amino     = sq->n;

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
      p7_GViterbi_SplicedGlobal_Dummy(dsq, gm_tr, pli->vit, pli->splice_scores->P_scores, pli->splice_scores->signal_scores, 1, L_dna_total, k_start, k_end, pli->min_intron);
      final_gC = pli->vit->xmx[L_dna_total * p7G_NXCELLS + p7G_C];

	  p7_fs_oprofile_SubConvert_Log(gm_tr, om_tr, k_start, k_end);
      p7_fs_oprofile_ReconfigLength_Log(om_tr, L_dna_total/3);
	  p7_omx_GrowTo_dpf(ox, sub_M, L_dna_total, L_dna_total);
	  p7_osplicescores_GrowTo(oss, sub_M);

      p7_Viterbi_SplicedGlobal(dsq, om_tr, ox, oss, 1, L_dna_total , pli->min_intron);
      final_oC = ox->xmx[L_dna_total * p7X_NXCELLS + p7X_C];

	  /* Scores must agree */
      if (fabs(final_gC - final_oC) > 0.001)
        esl_fatal("%s: generic %.4f != SSE %.4f", msg, final_gC, final_oC);
    }

  if (dsq != NULL) free(dsq);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_profile_fs_Destroy(gm_tr);
  p7_fs_oprofile_Destroy(om_tr); 
  p7_splicepipeline_Destroy(pli);
  p7_omx_Destroy(ox);
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
  { "-M",  eslARG_INT,    "200", NULL,"n>=50",NULL,NULL, NULL, "size of random models to sample",                0 },
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
