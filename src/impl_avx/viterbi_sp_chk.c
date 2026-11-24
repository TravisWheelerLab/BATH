/* Checkpointed spliced Viterbi dispatcher.
 * Provides the non-suffixed p7_Viterbi_SplicedChk() and
 * p7_Viterbi_SplicedTraceChk() declared in impl_avx.h.  Delegates to the
 * fastest available ISA implementation (viterbi_sp_chk_avx.c,
 * viterbi_sp_chk_sse.c).
 *
 * Contents:
 *    1. Dispatchers.
 *    2. Unit tests.
 *    3. Test driver.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/*****************************************************************
 * 1. Dispatchers.
 *****************************************************************/

/* Forward declaration of dispatcher */
static int p7_Viterbi_SplicedChk_Dispatcher(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                                            const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                                            int i_start, int i_end, int min_intron, int global_start, int global_end);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Viterbi_SplicedChk)(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                             const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                             int i_start, int i_end, int min_intron, int global_start, int global_end) = p7_Viterbi_SplicedChk_Dispatcher;

/* Function:  p7_Viterbi_SplicedChk_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for the checkpointed spliced
 *            Viterbi fill.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on bad inputs.
 */
static int
p7_Viterbi_SplicedChk_Dispatcher(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                                 const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                                 int i_start, int i_end, int min_intron, int global_start, int global_end)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Viterbi_SplicedChk = p7_Viterbi_SplicedChk_avx; return p7_Viterbi_SplicedChk_avx(sub_dsq, om_tr, cx, cv, signal_scores, acc_ov, don_ov, i_start, i_end, min_intron, global_start, global_end); }
#endif
#ifdef eslENABLE_SSE
  p7_Viterbi_SplicedChk = p7_Viterbi_SplicedChk_sse;
  return p7_Viterbi_SplicedChk_sse(sub_dsq, om_tr, cx, cv, signal_scores, acc_ov, don_ov, i_start, i_end, min_intron, global_start, global_end);
#else
  p7_Die("p7_Viterbi_SplicedChk: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/* Forward declaration of dispatcher */
static int p7_Viterbi_SplicedTraceChk_Dispatcher(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                                                 const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                                                 P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                                                 int min_intron, int global_start, int global_end, float *vitsc, float *max_drop);

/* Global function pointer, initially pointing at the dispatcher */
int (*p7_Viterbi_SplicedTraceChk)(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                                  const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                                  P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                                  int min_intron, int global_start, int global_end, float *vitsc, float *max_drop) = p7_Viterbi_SplicedTraceChk_Dispatcher;

/* Function:  p7_Viterbi_SplicedTraceChk_Dispatcher()
 *
 * Purpose:   Self-patching dispatcher for the checkpointed spliced
 *            Viterbi traceback.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile is not 1-codon-length or the matrix
 *                        isn't filled for this comparison.
 *            <eslFAIL>   if traceback fails.
 */
static int
p7_Viterbi_SplicedTraceChk_Dispatcher(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv,
                                      const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                                      P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                                      int min_intron, int global_start, int global_end, float *vitsc, float *max_drop)
{
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_Viterbi_SplicedTraceChk = p7_Viterbi_SplicedTraceChk_avx; return p7_Viterbi_SplicedTraceChk_avx(sub_dsq, om_tr, cx, cv, gm_tr, signal_scores, acc_ov, don_ov, tr, i_start, i_end, k_start, k_end, min_intron, global_start, global_end, vitsc, max_drop); }
#endif
#ifdef eslENABLE_SSE
  p7_Viterbi_SplicedTraceChk = p7_Viterbi_SplicedTraceChk_sse;
  return p7_Viterbi_SplicedTraceChk_sse(sub_dsq, om_tr, cx, cv, gm_tr, signal_scores, acc_ov, don_ov, tr, i_start, i_end, k_start, k_end, min_intron, global_start, global_end, vitsc, max_drop);
#else
  p7_Die("p7_Viterbi_SplicedTraceChk: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/*****************************************************************
 * 2. Unit tests.
 *****************************************************************/
#ifdef p7VITERBI_SP_CHK_TESTDRIVE

#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "p7_splice.h"

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

#ifdef eslENABLE_AVX
/* utest_chk_avx()
 *
 * For sequences with a simulated intron, fill with
 * p7_Viterbi_Spliced_avx() and with p7_Viterbi_SplicedChk_avx() at
 * several block sizes, in all three fill modes, and require the same
 * special states on every row. If the final C state is reachable, trace
 * both with p7_Viterbi_SplicedTrace_avx() and
 * p7_Viterbi_SplicedTraceChk_avx() and require the same trace, Viterbi
 * score and max_drop. Skipped if AVX is not available at runtime.
 */
static void
utest_chk_avx(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_CODONTABLE *codon_table,
              int M, int N, int intron_len)
{
  char            msg[]        = "utest_chk_avx failed";
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

  if (!esl_cpu_has_avx()) goto cleanup;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, M, p7_UNILOCAL);
  om_tr  = p7_fs_oprofile_Create_avx(M, abcAA, 1);
  p7_fs_oprofile_Convert_Log_avx(gm_tr, om_tr);
  ox     = p7_omx_Create_dpf_avx(M, M * 3, M * 3, p7X_NSCELLS);
  acc_ov = p7_oivx_Create_avx(M, SPLICE_ROWS);
  don_ov = p7_oivx_Create_avx(M, SIGNAL_MEM_SIZE);

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

      p7_fs_oprofile_SubConvert_Log_avx(gm_tr, om_tr, k_start, k_end);
      p7_fs_oprofile_ReconfigLength_Log_avx(om_tr, L_dna_total / 3);
      p7_fs_ReconfigLength(gm_tr, L_dna_total / 3);
      p7_omx_GrowTo_dpf_avx(ox, sub_M, L_dna_total, L_dna_total);
      p7_oivx_GrowTo_avx(acc_ov, sub_M, SPLICE_ROWS);
      p7_oivx_GrowTo_avx(don_ov, sub_M, SIGNAL_MEM_SIZE);

      for (md = 0; md < 3; md++)
        {
          p7_Viterbi_Spliced_avx(dsq, om_tr, ox, signal_scores, acc_ov, don_ov, 1, L_dna_total, min_intron, modes[md][0], modes[md][1]);

          for (bi = 0; bi < 3; bi++)
            {
              P7_CHKOMX  *cx  = p7_chkomx_Create_avx(sub_M, L_dna_total, blocks[bi], min_intron + 6, p7X_NSCELLS, p7_CHK_FORWARD);
              P7_CHKOIVX *cv  = p7_chkoivx_Create_avx(sub_M, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, L_dna_total / blocks[bi]);
              P7_OIVX    *acc = p7_oivx_Create_avx(sub_M, SPLICE_ROWS);
              P7_OIVX    *don = p7_oivx_Create_avx(sub_M, SIGNAL_MEM_SIZE);
              float       sc_full, sc_chk, drop_full, drop_chk;

              if (cx == NULL || cv == NULL || acc == NULL || don == NULL) esl_fatal("%s: allocation failed", msg);
              if (p7_Viterbi_SplicedChk_avx(dsq, om_tr, cx, cv, signal_scores, acc, don, 1, L_dna_total, min_intron, modes[md][0], modes[md][1]) != eslOK)
                esl_fatal("%s: checkpointed fill failed", msg);

              for (i = 0; i <= L_dna_total; i++)
                for (s = 0; s < p7X_NXCELLS; s++)
                  if (cx->xmx[i*p7X_NXCELLS+s] != ox->xmx[i*p7X_NXCELLS+s])
                    esl_fatal("%s: special state %d differs at row %d (block size %d, mode %d)", msg, s, i, blocks[bi], md);

              if (ox->xmx[L_dna_total*p7X_NXCELLS+p7X_C] != -eslINFINITY)
                {
                  if (p7_Viterbi_SplicedTrace_avx(dsq, ox, gm_tr, signal_scores, tr_full, 1, L_dna_total, k_start, k_end, min_intron, &sc_full, &drop_full) != eslOK)
                    esl_fatal("%s: replay trace failed", msg);
                  if (p7_Viterbi_SplicedTraceChk_avx(dsq, om_tr, cx, cv, gm_tr, signal_scores, acc, don, tr_chk, 1, L_dna_total, k_start, k_end,
                                                     min_intron, modes[md][0], modes[md][1], &sc_chk, &drop_chk) != eslOK)
                    esl_fatal("%s: checkpointed trace failed (block size %d, mode %d)", msg, blocks[bi], md);
                  utest_compare_traces(tr_full, sc_full, tr_chk, sc_chk, msg);
                  if (drop_full != drop_chk) esl_fatal("%s: max_drop %f != %f (block size %d, mode %d)", msg, drop_chk, drop_full, blocks[bi], md);
                  p7_trace_Reuse(tr_full);
                  p7_trace_Reuse(tr_chk);
                }

              p7_chkomx_Destroy_avx(cx);
              p7_chkoivx_Destroy_avx(cv);
              p7_oivx_Destroy_avx(acc);
              p7_oivx_Destroy_avx(don);
            }
        }
    }

 cleanup:
  if (dsq)    free(dsq);
  if (ox)     p7_omx_Destroy_avx(ox);
  if (acc_ov) p7_oivx_Destroy_avx(acc_ov);
  if (don_ov) p7_oivx_Destroy_avx(don_ov);
  if (om_tr)  p7_fs_oprofile_Destroy_avx(om_tr);
  if (hmm)    p7_hmm_Destroy(hmm);
  p7_trace_fs_Destroy(tr_full);
  p7_trace_fs_Destroy(tr_chk);
  esl_sq_Destroy(sq);
  p7_profile_fs_Destroy(gm_tr);
}
#endif /* eslENABLE_AVX */

#ifdef eslENABLE_SSE
/* utest_chk_sse()
 *
 * For sequences with a simulated intron, fill with
 * p7_Viterbi_Spliced_sse() and with p7_Viterbi_SplicedChk_sse() at
 * several block sizes, in all three fill modes, and require the same
 * special states on every row. If the final C state is reachable, trace
 * both with p7_Viterbi_SplicedTrace_sse() and
 * p7_Viterbi_SplicedTraceChk_sse() and require the same trace, Viterbi
 * score and max_drop.
 */
static void
utest_chk_sse(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_CODONTABLE *codon_table,
              int M, int N, int intron_len)
{
  char            msg[]        = "utest_chk_sse failed";
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

  if (!esl_cpu_has_sse()) goto cleanup;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_tr, M, p7_UNILOCAL);
  om_tr  = p7_fs_oprofile_Create_sse(M, abcAA, 1);
  p7_fs_oprofile_Convert_Log_sse(gm_tr, om_tr);
  ox     = p7_omx_Create_dpf_sse(M, M * 3, M * 3, p7X_NSCELLS);
  acc_ov = p7_oivx_Create_sse(M, SPLICE_ROWS);
  don_ov = p7_oivx_Create_sse(M, SIGNAL_MEM_SIZE);

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

      p7_fs_oprofile_SubConvert_Log_sse(gm_tr, om_tr, k_start, k_end);
      p7_fs_oprofile_ReconfigLength_Log_sse(om_tr, L_dna_total / 3);
      p7_fs_ReconfigLength(gm_tr, L_dna_total / 3);
      p7_omx_GrowTo_dpf_sse(ox, sub_M, L_dna_total, L_dna_total);
      p7_oivx_GrowTo_sse(acc_ov, sub_M, SPLICE_ROWS);
      p7_oivx_GrowTo_sse(don_ov, sub_M, SIGNAL_MEM_SIZE);

      for (md = 0; md < 3; md++)
        {
          p7_Viterbi_Spliced_sse(dsq, om_tr, ox, signal_scores, acc_ov, don_ov, 1, L_dna_total, min_intron, modes[md][0], modes[md][1]);

          for (bi = 0; bi < 3; bi++)
            {
              P7_CHKOMX  *cx  = p7_chkomx_Create_sse(sub_M, L_dna_total, blocks[bi], min_intron + 6, p7X_NSCELLS, p7_CHK_FORWARD);
              P7_CHKOIVX *cv  = p7_chkoivx_Create_sse(sub_M, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, L_dna_total / blocks[bi]);
              P7_OIVX    *acc = p7_oivx_Create_sse(sub_M, SPLICE_ROWS);
              P7_OIVX    *don = p7_oivx_Create_sse(sub_M, SIGNAL_MEM_SIZE);
              float       sc_full, sc_chk, drop_full, drop_chk;

              if (cx == NULL || cv == NULL || acc == NULL || don == NULL) esl_fatal("%s: allocation failed", msg);
              if (p7_Viterbi_SplicedChk_sse(dsq, om_tr, cx, cv, signal_scores, acc, don, 1, L_dna_total, min_intron, modes[md][0], modes[md][1]) != eslOK)
                esl_fatal("%s: checkpointed fill failed", msg);

              for (i = 0; i <= L_dna_total; i++)
                for (s = 0; s < p7X_NXCELLS; s++)
                  if (cx->xmx[i*p7X_NXCELLS+s] != ox->xmx[i*p7X_NXCELLS+s])
                    esl_fatal("%s: special state %d differs at row %d (block size %d, mode %d)", msg, s, i, blocks[bi], md);

              if (ox->xmx[L_dna_total*p7X_NXCELLS+p7X_C] != -eslINFINITY)
                {
                  if (p7_Viterbi_SplicedTrace_sse(dsq, ox, gm_tr, signal_scores, tr_full, 1, L_dna_total, k_start, k_end, min_intron, &sc_full, &drop_full) != eslOK)
                    esl_fatal("%s: replay trace failed", msg);
                  if (p7_Viterbi_SplicedTraceChk_sse(dsq, om_tr, cx, cv, gm_tr, signal_scores, acc, don, tr_chk, 1, L_dna_total, k_start, k_end,
                                                     min_intron, modes[md][0], modes[md][1], &sc_chk, &drop_chk) != eslOK)
                    esl_fatal("%s: checkpointed trace failed (block size %d, mode %d)", msg, blocks[bi], md);
                  utest_compare_traces(tr_full, sc_full, tr_chk, sc_chk, msg);
                  if (drop_full != drop_chk) esl_fatal("%s: max_drop %f != %f (block size %d, mode %d)", msg, drop_chk, drop_full, blocks[bi], md);
                  p7_trace_Reuse(tr_full);
                  p7_trace_Reuse(tr_chk);
                }

              p7_chkomx_Destroy_sse(cx);
              p7_chkoivx_Destroy_sse(cv);
              p7_oivx_Destroy_sse(acc);
              p7_oivx_Destroy_sse(don);
            }
        }
    }

 cleanup:
  if (dsq)    free(dsq);
  if (ox)     p7_omx_Destroy_sse(ox);
  if (acc_ov) p7_oivx_Destroy_sse(acc_ov);
  if (don_ov) p7_oivx_Destroy_sse(don_ov);
  if (om_tr)  p7_fs_oprofile_Destroy_sse(om_tr);
  if (hmm)    p7_hmm_Destroy(hmm);
  p7_trace_fs_Destroy(tr_full);
  p7_trace_fs_Destroy(tr_chk);
  esl_sq_Destroy(sq);
  p7_profile_fs_Destroy(gm_tr);
}
#endif /* eslENABLE_SSE */

#endif /*p7VITERBI_SP_CHK_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/


/*****************************************************************
 * 3. Test driver.
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

#ifdef eslENABLE_AVX
  utest_chk_avx(r, abcAA, gcode, bgAA, ct, M, N, I);
#endif
#ifdef eslENABLE_SSE
  utest_chk_sse(r, abcAA, gcode, bgAA, ct, M, N, I);
#endif

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
