/* Spliced Viterbi algorithms dispatcher.
 * Provides the non-suffixed p7_Viterbi_Spliced() and p7_Viterbi_SplicedTrace()
 * declared in impl_avx.h.  Delegates to the fastest available ISA implementation.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Viterbi_Spliced()
 *
 * Purpose:   Dispatch spliced Viterbi DP fill to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on bad inputs.
 */
int
p7_Viterbi_Spliced(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_OMX *ox,
                   const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov,
                   int i_start, int i_end, int min_intron, int global_start, int global_end)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Viterbi_Spliced_avx512(sub_dsq, om_tr, ox, signal_scores, acc_ov, don_ov, i_start, i_end, min_intron, global_start, global_end);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Viterbi_Spliced_sse(sub_dsq, om_tr, ox, signal_scores, acc_ov, don_ov, i_start, i_end, min_intron, global_start, global_end);
#endif
#ifdef eslENABLE_SSE
  return p7_Viterbi_Spliced_sse(sub_dsq, om_tr, ox, signal_scores, acc_ov, don_ov, i_start, i_end, min_intron, global_start, global_end);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_Viterbi_SplicedTrace()
 *
 * Purpose:   Dispatch spliced Viterbi traceback to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if profile is not 1-codon-length.
 *            <eslEFAIL>  if traceback fails.
 */
int
p7_Viterbi_SplicedTrace(const ESL_DSQ *sub_dsq, const P7_OMX *ox,
                        const P7_FS_PROFILE *gm_tr, const float *signal_scores,
                        P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end,
                        int min_intron, float *vitsc)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Viterbi_SplicedTrace_avx512(sub_dsq, ox, gm_tr, signal_scores, tr, i_start, i_end, k_start, k_end, min_intron, vitsc);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Viterbi_SplicedTrace_sse(sub_dsq, ox, gm_tr, signal_scores, tr, i_start, i_end, k_start, k_end, min_intron, vitsc);
#endif
#ifdef eslENABLE_SSE
  return p7_Viterbi_SplicedTrace_sse(sub_dsq, ox, gm_tr, signal_scores, tr, i_start, i_end, k_start, k_end, min_intron, vitsc);
#else
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 4. Unit tests.
 *****************************************************************/
#ifdef p7VITERBI_SP_TESTDRIVE

#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "p7_splice.h"

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
  P7_TRACE       *otr         = p7_trace_fs_Create();
  P7_TRACE       *gtr         = p7_trace_fs_Create();
  ESL_DSQ        *dsq         = NULL;
  P7_GMX         *gx          = NULL;
  P7_IVX         *acc_iv      = NULL;
  P7_IVX         *don_iv      = NULL;
  SPLICE_PIPELINE *pli        = NULL;
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
  gx     = p7_gmx_Create(M, M, M, p7G_NSCELLS);
  acc_iv = p7_ivx_Create(M, SPLICE_ROWS);
  don_iv = p7_ivx_Create(M, SIGNAL_MEM_SIZE);
  float signal_scores[p7S_SPLICE_SIGNALS];
  p7_SignalScores(signal_scores);

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

      /* --- Score comparison: generic Dummy and SSE Dummy must agree --- */
	  p7_fs_oprofile_SubConvert_Log(gm_tr, om_tr, k_start, k_end);
      p7_fs_oprofile_ReconfigLength_Log(om_tr, L_dna_total/3);
	  p7_omx_GrowTo_dpf(pli->vit, sub_M, L_dna_total, L_dna_total);
	  p7_oivx_GrowTo(pli->acc_ov, sub_M, SPLICE_ROWS);
	  p7_oivx_GrowTo(pli->don_ov, sub_M, SIGNAL_MEM_SIZE);

      p7_Viterbi_Spliced(dsq, om_tr, pli->vit, signal_scores, pli->acc_ov, pli->don_ov, 1, L_dna_total, pli->min_intron, TRUE, TRUE);
      final_oC = pli->vit->xmx[L_dna_total * p7X_NXCELLS + p7X_C];

      p7_fs_ReconfigLength(gm_tr, L_dna_total/3);
      p7_gmx_GrowTo(gx, sub_M, L_dna_total, L_dna_total);
      p7_ivx_GrowTo(acc_iv, sub_M, SPLICE_ROWS);
	  p7_ivx_GrowTo(don_iv, sub_M, SIGNAL_MEM_SIZE);

      p7_GViterbi_Spliced(dsq, gm_tr, gx, acc_iv, don_iv, signal_scores, 1, L_dna_total, k_start, k_end, pli->min_intron, TRUE, TRUE);
      final_gC = gx->xmx[L_dna_total * p7G_NXCELLS + p7G_C];

      if (fabs(final_gC - final_oC) > 0.001) 
        esl_fatal("%s: generic %.4f != SSE %.4f", msg, final_gC, final_oC); 

      if(final_gC == -eslINFINITY) continue;   
      
      p7_Viterbi_SplicedTrace(dsq, pli->vit, gm_tr, signal_scores, otr, 1, L_dna_total, k_start, k_end, pli->min_intron, NULL);
      p7_GViterbi_SplicedTrace(dsq, gm_tr, gx, signal_scores, gtr, 1, L_dna_total, k_start, k_end, pli->min_intron, NULL); 

      p7_trace_Compare(otr, gtr, 0.0);
        
      p7_trace_Reuse(otr);
      p7_trace_Reuse(gtr);

      p7_Viterbi_Spliced(dsq, om_tr, pli->vit, signal_scores, pli->acc_ov, pli->don_ov, 1, L_dna_total, pli->min_intron, TRUE, FALSE);
      final_oC = pli->vit->xmx[L_dna_total * p7X_NXCELLS + p7X_C];

      p7_GViterbi_Spliced(dsq, gm_tr, gx, acc_iv, don_iv, signal_scores, 1, L_dna_total, k_start, k_end, pli->min_intron, TRUE, FALSE);
      final_gC = gx->xmx[L_dna_total * p7G_NXCELLS + p7G_C];
       
      p7_Viterbi_SplicedTrace(dsq, pli->vit, gm_tr, signal_scores, otr, 1, L_dna_total, k_start, k_end, pli->min_intron, NULL);
      p7_GViterbi_SplicedTrace(dsq, gm_tr, gx, signal_scores, gtr, 1, L_dna_total, k_start, k_end, pli->min_intron, NULL);

      p7_trace_Compare(otr, gtr, 0.0);

      p7_trace_Reuse(otr);
      p7_trace_Reuse(gtr);

      p7_Viterbi_Spliced(dsq, om_tr, pli->vit, signal_scores, pli->acc_ov, pli->don_ov, 1, L_dna_total, pli->min_intron, FALSE, TRUE);
      final_oC = pli->vit->xmx[L_dna_total * p7X_NXCELLS + p7X_C];

      p7_GViterbi_Spliced(dsq, gm_tr, gx, acc_iv, don_iv, signal_scores, 1, L_dna_total, k_start, k_end, pli->min_intron, FALSE, TRUE);
      final_gC = gx->xmx[L_dna_total * p7G_NXCELLS + p7G_C];
      
      p7_Viterbi_SplicedTrace(dsq, pli->vit, gm_tr, signal_scores, otr, 1, L_dna_total, k_start, k_end, pli->min_intron, NULL);
      p7_GViterbi_SplicedTrace(dsq, gm_tr, gx, signal_scores, gtr, 1, L_dna_total, k_start, k_end, pli->min_intron, NULL);

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
  p7_gmx_Destroy(gx);
  p7_ivx_Destroy(acc_iv);
  p7_ivx_Destroy(don_iv);
  p7_trace_fs_Destroy(otr);
  p7_trace_fs_Destroy(gtr);
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
#include "impl_avx.h"

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
