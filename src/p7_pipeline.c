/* H3's accelerated seq/profile comparison pipeline
 * Modified for use by BATH
 *
 * Contents:
 *   1. P7_PIPELINE: allocation, initialization, destruction
 *   2. Pipeline API
 */

#include "p7_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h> 

#include "easel.h"
#include "esl_exponential.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"
#include "esl_gencode.h"
#include "hmmer.h"
#include "p7_splice.h"

/* Struct used to pass a collection of useful temporary objects around
 * within the Frameshift functions (BATH)
 */
typedef struct {
  ESL_SQ           *tmpseq;     // - a new or reused digital sequence object used for p7_alidisplay_Create() call
  P7_OMX          **oxf_holder; // - a temporary list of forward parser matricies for ORFs
  float            *fwdsc;
  double           *P_orf;      // - a temporary list or forwrad P values for ORFs
} P7_PIPELINE_OBJS;


/*****************************************************************
 * 1. The P7_PIPELINE object: allocation, initialization, destruction.
 *****************************************************************/

/* Function:  p7_pipeline_Create_BATH() 
 * Synopsis:  Create a new BATH comparison pipeline.
 *
 * Purpose:   Given an application configuration structure <go>
 *            containing certain standardized options (described
 *            below), and some initial guesses at the model size <M_hint>
 *            and sequence length <L_hint> that will be processed,
 *            create new pipeline object. Only <p7_SEARCH_SEQS> mode 
 *            (searching one model against a sequence database) is 
 *            currently supported.
 *
 *            In search mode, we would generally know the length of
 *            our query profile exactly, and would pass <om->M> as <M_hint>.
 *            Targets will come in various sizes as we read them,
 *            and the pipeline will resize any necessary objects as
 *            needed, so the other (unknown) length is only an
 *            initial allocation.
 *            
 *            The configuration <go> must include settings for the 
 *            following options:
 *            
 *            || option      ||            description                    || usually  ||
 *            | --noali      |  don't output alignments (smaller output)   |   FALSE   |
 *            | -E           |  report hits <= this E-value threshold      |    10.0   |
 *            | -T           |  report hits >= this bit score threshold    |    NULL   |
 *            | -Z           |  set initial hit search space size          |    NULL   |
 *            | --domZ       |  set domain search space size               |    NULL   |
 *            | --domE       |  report domains <= this E-value threshold   |    10.0   |
 *            | --domT       |  report domains <= this bit score threshold |    NULL   |
 *            | --incE       |  include hits <= this E-value threshold     |    0.01   |
 *            | --incT       |  include hits >= this bit score threshold   |    NULL   |
 *            | --incdomE    |  include domains <= this E-value threshold  |    0.01   |
 *            | --incdomT    |  include domains <= this score threshold    |    NULL   |
 *            | --cut_ga     |  model-specific thresholding using GA       |   FALSE   |
 *            | --cut_nc     |  model-specific thresholding using NC       |   FALSE   |
 *            | --cut_tc     |  model-specific thresholding using TC       |   FALSE   |
 *            | --max        |  turn all heuristic filters off             |   FALSE   |
 *            | --F1         |  Stage 1 (MSV) thresh: promote hits P <= F1 |    0.02   |
 *            | --F2         |  Stage 2 (Vit) thresh: promote hits P <= F2 |    1e-3   |
 *            | --F3         |  Stage 3 (Fwd) thresh: promote hits P <= F3 |    1e-5   |
 *            | --F4         |  Stage 3 (FS-Fwd) thresh: promote hits P <= F4 |    5e-4   |
 *            | --nobias     |  turn OFF composition bias filter HMM       |   FALSE   |
 *            | --nonull2    |  turn OFF biased comp score correction      |   FALSE   |
 *            | --seed       |  RNG seed (0=use arbitrary seed)            |      42   |
 *            | --acc        |  prefer accessions over names in output     |   FALSE   |
 *
 *            As a special case, if <go> is <NULL>, defaults are set as above.
 *            This shortcut is used in simplifying test programs and the like.
 *            
 * Returns:   ptr to new <P7_PIPELINE> object on success. Caller frees this
 *            with <p7_pipeline_Destroy()>.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_PIPELINE *
p7_pipeline_Create_BATH(ESL_GETOPTS *go, int M_hint, int L_hint, enum p7_pipemodes_e mode)
{
  P7_PIPELINE *pli  = NULL;
  int          seed = (go ? esl_opt_GetInteger(go, "--seed") : 42);
  int          status;

  ESL_ALLOC(pli, sizeof(P7_PIPELINE));

  pli->do_alignment_score_calc = (go ? esl_opt_IsUsed(go, "--splice") : 0);

   /* Set Alignment Mode */ 
  pli->spliced =  (go ? esl_opt_IsUsed(go, "--splice") : 0); 
  pli->fs_pipe  = (go ? (esl_opt_IsUsed(go, "--fs") || esl_opt_IsUsed(go, "--fsonly")) : 0); 
  pli->std_pipe = (go ? !esl_opt_IsUsed(go, "--fsonly") : 1);

  /* Create sparce memeory forward and backward optimized matricies for use in the 
   * non-frameshift pipeline branch
   */
   if ((pli->oxf  = p7_omx_Create(M_hint, 0, L_hint)) == NULL) goto ERROR;
   if ((pli->oxb  = p7_omx_Create(M_hint, 0, L_hint)) == NULL) goto ERROR;

   /* Create full memeory forward and backward optimized matricies for use in the
   * non-frameshift pipeline branch
   */
   if ((pli->fwd  = p7_omx_Create(M_hint, L_hint, L_hint)) == NULL) goto ERROR;
   if ((pli->bck  = p7_omx_Create(M_hint, L_hint, L_hint)) == NULL) goto ERROR;

  /* Create sparce memeory forward and backward generic frameshift aware matricies 
   * for use in the frameshift pipeline filters
   */  
    if ((pli->oxf_fs  = p7_omx_Create_dpf(M_hint, PARSER_ROWS_FWD, L_hint, p7G_NSCELLS)) == NULL) goto ERROR;
    if ((pli->oxb_fs  = p7_omx_Create_dpf(M_hint, PARSER_ROWS_BWD, L_hint, p7G_NSCELLS)) == NULL) goto ERROR;
    if ((pli->ov3     = p7_oivx_Create(M_hint, p7P_3CODONS))                              == NULL) goto ERROR;
    if ((pli->ov5     = p7_oivx_Create(M_hint, p7P_5CODONS))                              == NULL) goto ERROR;

  /* Create full memeory forward, backward and posterior generic frameshift aware matricies
   * for use in the frameshift pipeline alignment
   */ 
   if ((pli->fwd_fs = p7_omx_Create_dpf(M_hint, L_hint, L_hint, p7G_NSCELLS_FS)) == NULL) goto ERROR;
   if ((pli->bck_fs = p7_omx_Create_dpf(M_hint, L_hint, L_hint, p7G_NSCELLS))    == NULL) goto ERROR;

  /* Normally, we reinitialize the RNG to the original seed every time we're
   * about to collect a stochastic trace ensemble. This eliminates run-to-run
   * variability. As a special case, if seed==0, we choose an arbitrary one-time 
   * seed: time() sets the seed, and we turn off the reinitialization.
   */
   pli->r                  =  esl_randomness_CreateFast(seed);
   pli->do_reseeding       = (seed == 0) ? FALSE : TRUE;
   pli->ddef               = p7_domaindef_Create_BATH(pli->r, go);
   pli->ddef->do_reseeding = pli->do_reseeding;

   /* Configure reporting thresholds */
   pli->by_E            = TRUE;
   pli->E               = (go ? esl_opt_GetReal(go, "-E") : 10.0);
   pli->T               = 0.0;
   pli->dom_by_E        = TRUE;
   pli->domE            = (go ? esl_opt_GetReal(go, "--domE") : 10.0);
   pli->domT            = 0.0;
   pli->use_bit_cutoffs = FALSE;
   if (go && esl_opt_IsOn(go, "-T")) 
   {
     pli->T    = esl_opt_GetReal(go, "-T");  
     pli->by_E = FALSE;
   }
   if (go && esl_opt_IsOn(go, "--domT")) 
   {
     pli->domT     = esl_opt_GetReal(go, "--domT"); 
     pli->dom_by_E = FALSE;
   }

   /* Configure inclusion thresholds */
   pli->inc_by_E           = TRUE;
   pli->incE               = (go ? esl_opt_GetReal(go, "--incE") : 0.01);
   pli->incT               = 0.0;
   pli->incdom_by_E        = TRUE;
   pli->incdomE            = (go ? esl_opt_GetReal(go, "--incdomE") : 0.01);
   pli->incdomT            = 0.0;
   if (go && esl_opt_IsOn(go, "--incT")) 
   {
     pli->incT     = esl_opt_GetReal(go, "--incT"); 
     pli->inc_by_E = FALSE;
   } 
   if (go && esl_opt_IsOn(go, "--incdomT")) 
   {
     pli->incdomT     = esl_opt_GetReal(go, "--incdomT"); 
     pli->incdom_by_E = FALSE;
   }

   /* Configure search space sizes for E value calculations  */
   pli->Z       = pli->domZ       = 0.0;
   pli->Z_setby = pli->domZ_setby = p7_ZSETBY_NTARGETS;
   if (go && esl_opt_IsOn(go, "-Z")) 
   {
     pli->Z_setby = p7_ZSETBY_OPTION;
     pli->Z       = esl_opt_GetReal(go, "-Z");
   }
   if (go && esl_opt_IsOn(go, "--domZ")) 
   {
     pli->domZ_setby = p7_ZSETBY_OPTION;
     pli->domZ       = esl_opt_GetReal(go, "--domZ");
   }

   /* Configure acceleration pipeline thresholds */
   pli->do_max        = FALSE;
   pli->do_biasfilter = TRUE;
   pli->do_null2      = TRUE;
   pli->F1     = ((go && esl_opt_IsOn(go, "--F1")) ? ESL_MIN(1.0, esl_opt_GetReal(go, "--F1")) : 0.02);
   pli->F2     = (go ? ESL_MIN(1.0, esl_opt_GetReal(go, "--F2")) : 1e-3);
   pli->F3     = (go ? ESL_MIN(1.0, esl_opt_GetReal(go, "--F3")) : 1e-5);
   pli->F4     = (go ? ESL_MIN(1.0, esl_opt_GetReal(go, "--F4")) : 5e-4);

   if (go && esl_opt_GetBoolean(go, "--max")) 
   {
    pli->do_max        = TRUE;
    pli->do_biasfilter = FALSE;
    
    pli->F1 = pli->F2 = pli->F3 = pli->F4 = 1.0;

   }
   if (go && esl_opt_GetBoolean(go, "--nonull2")) pli->do_null2      = FALSE;
   if (go && esl_opt_GetBoolean(go, "--nobias"))  pli->do_biasfilter = FALSE;
  
   /* Accounting as we collect results */
   pli->nmodels         = 0;
   pli->nseqs           = 0;
   pli->nres            = 0;
   pli->nnodes          = 0;
   pli->n_past_msv      = 0;
   pli->n_past_bias     = 0;
   pli->n_past_vit      = 0;
   pli->n_past_fwd      = 0;
   pli->pos_past_msv    = 0;
   pli->pos_past_bias   = 0;
   pli->pos_past_vit    = 0;
   pli->pos_past_fwd    = 0;
   pli->mode            = mode;
   pli->show_accessions = (go && esl_opt_GetBoolean(go, "--acc")   ? TRUE  : FALSE);
   pli->show_alignments = (go && esl_opt_GetBoolean(go, "--noali") ? FALSE : TRUE);
   pli->show_translated_sequence = (go && esl_opt_GetBoolean(go, "--notrans") ? FALSE : TRUE); /* TRUE to display translated DNA sequence in alignment display for bathsearch */
   pli->show_frameline = (go && esl_opt_GetBoolean(go, "--frameline") ? TRUE : FALSE); /* TRUE to display the frame of each codon in alignment display for bathsearch */
   pli->show_cigar     = (go && esl_opt_GetBoolean(go, "--cigar") ? TRUE : FALSE); /* TRUE to alignment CIGAR string int tabular output for bathsearch */
   pli->hfp             = NULL;
   pli->errbuf[0]       = '\0';

   return pli;

ERROR:
  p7_pipeline_Destroy_BATH(pli);
  return NULL;
}

/* Function:  p7_pipeline_Reuse_BATH() 
 * Synopsis:  Reuse a BATH for next target.
 *
 */
int
p7_pipeline_Reuse_BATH(P7_PIPELINE *pli)
{
  p7_omx_Reuse(pli->oxf_fs);
  p7_omx_Reuse(pli->oxb_fs);
  p7_omx_Reuse(pli->fwd_fs);
  p7_omx_Reuse(pli->bck_fs);
  p7_omx_Reuse(pli->oxf);
  p7_omx_Reuse(pli->oxb);
  p7_omx_Reuse(pli->fwd);
  p7_omx_Reuse(pli->bck);
  p7_domaindef_Reuse(pli->ddef);
  return eslOK;
}

/* Function:  p7_pipeline_Destroy_BATH() 
 * Synopsis:  Free a BATH <P7_PIPELINE> object.
 *
 */
void
p7_pipeline_Destroy_BATH(P7_PIPELINE *pli)
{
  if (pli == NULL) return;
 
  p7_omx_Destroy(pli->oxf_fs);
  p7_omx_Destroy(pli->oxb_fs);
  p7_omx_Destroy(pli->fwd_fs);
  p7_omx_Destroy(pli->bck_fs);
  p7_oivx_Destroy(pli->ov3);
  p7_oivx_Destroy(pli->ov5);
  p7_omx_Destroy(pli->oxf);
  p7_omx_Destroy(pli->oxb);
  p7_omx_Destroy(pli->fwd);
  p7_omx_Destroy(pli->bck);
  esl_randomness_Destroy(pli->r);
  p7_domaindef_Destroy_BATH(pli->ddef);
  free(pli);
}

/*---------------- end, P7_PIPELINE object ----------------------*/


/*****************************************************************
 * 2. The pipeline API.
 *****************************************************************/

/* Function:  p7_pli_ExtendAndMergeWindows_BATH
 * Synopsis:  Turns a list of ORF ssv diagonals into DNA windows, and 
 *            merges overlapping windows.
 *
 * Purpose:   Accepts a <orf_block> of and creates a set of SSV 
 *            diagonals. For each ORFs best (higest scoreing)
 *            diagonal, extends those windows based on a combination 
 *            of the max_length value from <om> and the prefix and 
 *            suffix lengths stored in <data>, and converts them to
 *            DNA coordinates. Then merges (in place) windows that 
 *            overlap by more than <pct_overlap> percent, ensuring 
 *            that windows stay within the bounds of 1..<L>.
 *
 * Returns:   <eslOK>
 */
int
p7_pli_ExtendAndMergeWindows_BATH(P7_PIPELINE *pli, ESL_SQ_BLOCK *orf_block, ESL_SQ *dnasq, P7_OPROFILE *om, P7_BG *bg, const P7_SCOREDATA *data, P7_HMM_WINDOWLIST *windowlist, float pct_overlap, int complementarity)
{ 

  int i, f;
  int best_window;
  float best_score;
  P7_HMM_WINDOWLIST    tmp_windowlist;
  P7_HMM_WINDOW        *prev_window = NULL;
  P7_HMM_WINDOW        *curr_window = NULL;
  int64_t              window_start;
  int64_t              window_end;
  int32_t              window_len;
  int64_t              overlap_start;
  int64_t              overlap_end;
  int32_t              overlap_len;
  int                  new_hit_cnt = 0;
  ESL_SQ *curr_orf;

  p7_hmmwindow_init(&tmp_windowlist);

  for(f = 0; f < orf_block->count; f++)
  {
    curr_orf = &(orf_block->list[f]);
  
    p7_oprofile_ReconfigLength(om, curr_orf->n);
    p7_omx_GrowTo(pli->oxf, om->M, 0, curr_orf->n);    /* expand the one-row omx if needed */ 

    p7_SSVFilter_BATH(curr_orf->dsq, curr_orf->n, om, pli->oxf, data, bg, pli->F1, &tmp_windowlist); 

    /* If the orf fails to produce a window use the full ORF coords aligned to the center of the model. */
    if(tmp_windowlist.count == 0) {
      if(curr_orf->n >= om->M) {
        window_start = (curr_orf->n - om->M) / 2 + 1;
        window_end   = om->M;
        window_len   = om->M;
      }
      else {
        window_start = 1;
        window_end   = om->M - ((om->M - curr_orf->n) / 2);
        window_len   = curr_orf->n;
      }
      p7_hmmwindow_new(&tmp_windowlist, 0, window_start, window_end, window_len, 0.0, 0, curr_orf->n);
    }

    best_window = 0;
    best_score  = tmp_windowlist.windows[0].score;
    for(i = 1; i < tmp_windowlist.count; i++) {
      if(tmp_windowlist.windows[i].score > best_score) {
        best_window = i;
        best_score  = tmp_windowlist.windows[i].score;
      }
    }

    curr_window = tmp_windowlist.windows+best_window;

    /* Extend ORF coords */
    window_start = curr_window->n -                       (om->max_length * (0.1 + data->prefix_lengths[curr_window->k - curr_window->length + 1])) + 1;
    window_end   = curr_window->n + curr_window->length + (om->max_length * (0.1 + data->suffix_lengths[curr_window->k])) - 2;

    window_start = ESL_MIN(0,           window_start); //move start to at least the begining of the ORF
    window_end   = ESL_MAX(curr_orf->n, window_end);   // move end to at least the end of the ORF    

    /* Convert to DNA coords */
    if(complementarity) {
      window_start   = ESL_MAX(1,        (dnasq->n - curr_orf->start + 1) + (window_start * 3));
      window_end     = ESL_MIN(dnasq->n, (dnasq->n - curr_orf->start + 1) + (window_end * 3));
    }
    else {
      window_start = ESL_MAX(1,        curr_orf->start + (window_start * 3));
      window_end   = ESL_MIN(dnasq->n, curr_orf->start + (window_end * 3));
    }
    
    p7_hmmwindow_new(windowlist, 0, window_start, curr_window->k, window_end-window_start+1, 0.0, complementarity, dnasq->n);    
    curr_orf->idx = windowlist->count - 1; // keep track of which window ORFs belong to
    tmp_windowlist.count = 0;
  }

  if (tmp_windowlist.windows != NULL) free (tmp_windowlist.windows);
  
  if( windowlist->count == 0) return eslOK;

  p7_hmmwindow_SortByStart(windowlist);
  new_hit_cnt = 0;

  /* merge overlapping windows, compressing list in place. */
  for (i=1; i<windowlist->count; i++) {
    prev_window = windowlist->windows+new_hit_cnt;
    curr_window = windowlist->windows+i;

    overlap_start = ESL_MAX(prev_window->n, curr_window->n);
    overlap_end   = ESL_MIN(prev_window->n+prev_window->length-1, curr_window->n+curr_window->length-1);
    overlap_len   = overlap_end - overlap_start + 1;

    window_start  = ESL_MIN(prev_window->n, curr_window->n);
    window_end    = ESL_MAX(prev_window->n+prev_window->length-1, curr_window->n+curr_window->length-1);
    window_len    = window_end - window_start + 1;


    if(((float)(overlap_len)/ESL_MIN(prev_window->length, curr_window->length) > pct_overlap) &&
      window_len < ( 2 * (om->max_length * 3)))
    {
      prev_window->n      = window_start;
      prev_window->length = window_len;
    } 
    else {
      new_hit_cnt++;
      windowlist->windows[new_hit_cnt] = windowlist->windows[i];
    }
    orf_block->list[i].idx = new_hit_cnt;
  }
  windowlist->count = new_hit_cnt+1;

  return eslOK;
}

static int
p7_pli_ComputeLocalCompo(const P7_SCOREDATA *data, const P7_OPROFILE *om, const P7_BG *bg, int k_start, int k_end, float *compo)
{
  int   k, x;
  int   k_len;
  float log_odds;

  /* Enforce minumum window length of 20 */
  k_len = k_end - k_start + 1;
  if(k_len < 20) {
    k_start -= (20-k_len)/2;
    k_end   += (20-k_len)/2;
  }

  /* Clamp to valid model range */
  k_start = ESL_MAX(1,     k_start);
  k_end   = ESL_MIN(om->M, k_end);

  /* Accumulate emission probabilities over the merged, padded span */
  esl_vec_FSet(compo, om->abc->K, 0.0f);
  for (k = k_start; k <= k_end; k++) {
    for (x = 0; x < om->abc->K; x++) {
      /* ssv_scores[k*Kp+x] is a byte emission cost: lower = more favored.
       * (base_b - cost) / scale_b recovers the log-odds emission score.
       * Multiplying by bg->f[x] gives the approximate emission probability. */
      log_odds = ((float)om->base_b - (float)data->ssv_scores[k * om->abc->Kp + x]) / om->scale_b;
      compo[x] += bg->f[x] * expf(log_odds);
    }
  }

  esl_vec_FNorm(compo, om->abc->K);
  return eslOK;
}


int
p7_pli_BuildDNAWindows(P7_PIPELINE *pli, ESL_SQ_BLOCK *orf_block, ESL_SQ *dnasq, P7_OPROFILE *om, P7_BG *bg, const P7_SCOREDATA *data, P7_HMM_WINDOWLIST *windowlist, float pct_overlap, P7_PIPELINE_OBJS *pli_tmp, P7_HMM_WINDOWLIST *hit_windows, int complementarity)
{

  int i, f, w;
  int best_window_idx;
  float best_score;
  P7_HMM_WINDOW        fallback_window;
  P7_HMM_WINDOW        *prev_window = NULL;
  P7_HMM_WINDOW        *curr_window = NULL;
  int64_t              window_start;
  int64_t              window_end;
  int32_t              window_len;
  int64_t              overlap_start;
  int64_t              overlap_end;
  int32_t              overlap_len;
  int                  new_hit_cnt = 0;
  ESL_SQ *curr_orf;

  for(f = 0; f < orf_block->count; f++)
  {
    if(pli_tmp->P_orf[f] > pli->F4) continue;

    curr_orf = &(orf_block->list[f]);

    /* Find the best precomputed window for this ORF (id == f).
     * Windows were produced by p7_ViterbiFilter_BATH or p7_SSVFilter_BATH
     * earlier in the pipeline and stamped with the ORF index in id. */
    best_window_idx = -1;
    best_score      = -eslINFINITY;
    for(w = 0; w < hit_windows->count; w++) {
      if(hit_windows->windows[w].id != f) continue;
      if(hit_windows->windows[w].score > best_score ||
         (hit_windows->windows[w].score == best_score &&
          hit_windows->windows[w].length > (best_window_idx >= 0 ? hit_windows->windows[best_window_idx].length : 0))) {
        best_score      = hit_windows->windows[w].score;
        best_window_idx = w;
      }
    }

    if(best_window_idx >= 0) {
      curr_window = &hit_windows->windows[best_window_idx];
    } else {
      /* No precomputed window for this ORF — fall back to center of model. */
      if(curr_orf->n >= om->M) {
        fallback_window.n      = (curr_orf->n - om->M) / 2 + 1;
        fallback_window.k      = om->M;
        fallback_window.length = om->M;
      } else {
        fallback_window.n      = 1;
        fallback_window.k      = om->M - ((om->M - curr_orf->n) / 2);
        fallback_window.length = curr_orf->n;
      }
      curr_window = &fallback_window;
    }

    /* Extend ORF coords using prefix/suffix lengths */
    window_start = curr_window->n -                       (om->max_length * (0.1 + data->prefix_lengths[curr_window->k - curr_window->length + 1])) + 1;
    window_end   = curr_window->n + curr_window->length + (om->max_length * (0.1 + data->suffix_lengths[curr_window->k])) - 2;

    window_start = ESL_MIN(0,           window_start); /* move start to at least the beginning of the ORF */
    window_end   = ESL_MAX(curr_orf->n, window_end);   /* move end   to at least the end of the ORF      */

    /* Convert to DNA coords */
    if(complementarity) {
      window_start = ESL_MAX(1,        (dnasq->n - curr_orf->start + 1) + (window_start * 3));
      window_end   = ESL_MIN(dnasq->n, (dnasq->n - curr_orf->start + 1) + (window_end * 3));
    }
    else {
      window_start = ESL_MAX(1,        curr_orf->start + (window_start * 3));
      window_end   = ESL_MIN(dnasq->n, curr_orf->start + (window_end * 3));
    }

    p7_hmmwindow_new(windowlist, 0, window_start, curr_window->k, window_end-window_start+1, 0.0, complementarity, dnasq->n);
    curr_orf->idx = windowlist->count - 1; /* keep track of which window ORFs belong to */
  }
  
  if( windowlist->count == 0) return eslOK;

  p7_hmmwindow_SortByStart(windowlist);
  new_hit_cnt = 0;

  /* merge overlapping windows, compressing list in place. */
  for (i=1; i<windowlist->count; i++) {
    prev_window = windowlist->windows+new_hit_cnt;
    curr_window = windowlist->windows+i;

    overlap_start = ESL_MAX(prev_window->n, curr_window->n);
    overlap_end   = ESL_MIN(prev_window->n+prev_window->length-1, curr_window->n+curr_window->length-1);
    overlap_len   = overlap_end - overlap_start + 1;

    window_start  = ESL_MIN(prev_window->n, curr_window->n);
    window_end    = ESL_MAX(prev_window->n+prev_window->length-1, curr_window->n+curr_window->length-1);
    window_len    = window_end - window_start + 1;


    if(((float)(overlap_len)/ESL_MIN(prev_window->length, curr_window->length) > pct_overlap) &&
      window_len < ( 2 * (om->max_length * 3)))
    {
      prev_window->n      = window_start;
      prev_window->length = window_len;
    } 
    else {
      new_hit_cnt++;
      windowlist->windows[new_hit_cnt] = windowlist->windows[i];
    }
    orf_block->list[i].idx = new_hit_cnt;
  }
  windowlist->count = new_hit_cnt+1;

  return eslOK;
}



/* Function:  p7_pli_TargetReportable
 * Synopsis:  Returns TRUE if target score meets reporting threshold.
 *
 * Purpose:   Returns <TRUE> if the bit score <score> and/or 
 *            log P-value <lnP> meet per-target reporting thresholds 
 *            for the processing pipeline.
 */
int
p7_pli_TargetReportable(P7_PIPELINE *pli, float score, double lnP)
{
  if      (  pli->by_E && exp(lnP) <= pli->E) return TRUE; 
  else if (! pli->by_E && score    >= pli->T) return TRUE;

  return FALSE;
}

/* Function:  p7_pli_TargetIncludable()
 * Synopsis:  Returns TRUE if target score meets inclusion threshold.
 */
int
p7_pli_TargetIncludable(P7_PIPELINE *pli, float score, double lnP)
{
  if      (  pli->inc_by_E && exp(lnP) <= pli->incE) return TRUE;
  else if (! pli->inc_by_E   && score  >= pli->incT) return TRUE;

  return FALSE;
}

/* Function:  p7_pli_NewModel()
 * Synopsis:  Prepare pipeline for a new model (target or query)
 *
 * Purpose:   Caller has a new model <om>. Prepare the pipeline <pli>
 *            to receive this model as either a query or a target.
 *
 *            If the "experimental" bias filter HMM is in use, this
 *            call resets it to use the new model's composition. This
 *            overwrites the bias filter HMM's expected length! You
 *            need to call <p7_bg_SetLength()> after a <NewModel()> call.
 *            (Failure to do this is bug #h85, 14 Dec 10.)
 *
 *            The pipeline may alter the null model <bg> in a model-specific
 *            way (if we're using a composition bias filter HMM in the
 *            pipeline).
 *
 * Returns:   <eslOK> on success.
 * 
 *            <eslEINVAL> if pipeline expects to be able to use a
 *            model's bit score thresholds, but this model does not
 *            have the appropriate ones set.
 */
int
p7_pli_NewModel(P7_PIPELINE *pli, const P7_OPROFILE *om, P7_BG *bg)
{
  int status = eslOK;

  pli->nmodels++;
  pli->nnodes += om->M;
  if (pli->Z_setby == p7_ZSETBY_NTARGETS && pli->mode == p7_SCAN_MODELS) pli->Z = pli->nmodels;

  if (pli->do_biasfilter) p7_bg_SetFilter(bg, om->M, om->compo);

  if (pli->mode == p7_SEARCH_SEQS)
    status = p7_pli_NewModelThresholds(pli, om);

  pli->W = om->max_length;

  return status;
}

/* Function:  p7_pli_NewModelThresholds()
 * Synopsis:  Set reporting and inclusion bit score thresholds on a new model.
 *
 * Purpose:   Set the bit score thresholds on a new model, if we're 
 *            using Pfam GA, TC, or NC cutoffs for reporting or
 *            inclusion.
 *            
 *            In a "search" pipeline, this only needs to be done once
 *            per query model, so <p7_pli_NewModelThresholds()> gets 
 *            called by <p7_pli_NewModel()>.
 *            
 *            In a "scan" pipeline, this needs to be called for each
 *            model, and it needs to be called after
 *            <p7_oprofile_ReadRest()>, because that's when the bit
 *            score thresholds get read.
 *
 * Returns:   <eslOK> on success. 
 *            
 *            <eslEINVAL> if pipeline expects to be able to use a
 *            model's bit score thresholds, but this model does not
 *            have the appropriate ones set.
 *
 * Xref:      Written to fix bug #h60.
 */
int
p7_pli_NewModelThresholds(P7_PIPELINE *pli, const P7_OPROFILE *om)
{

  if (pli->use_bit_cutoffs)
  {
    if (pli->use_bit_cutoffs == p7H_GA)
    {
      if (om->cutoff[p7_GA1] == p7_CUTOFF_UNSET)
        ESL_FAIL(eslEINVAL, pli->errbuf, "GA bit thresholds unavailable on model %s\n", om->name);
      pli->T    = pli->incT    = om->cutoff[p7_GA1];
      pli->domT = pli->incdomT = om->cutoff[p7_GA2];
    }
    else if  (pli->use_bit_cutoffs == p7H_TC)
    {
      if (om->cutoff[p7_TC1] == p7_CUTOFF_UNSET)
        ESL_FAIL(eslEINVAL, pli->errbuf, "TC bit thresholds unavailable on model %s\n", om->name);
      pli->T    = pli->incT    = om->cutoff[p7_TC1];
      pli->domT = pli->incdomT = om->cutoff[p7_TC2];
    }
    else if (pli->use_bit_cutoffs == p7H_NC)
    {
      if (om->cutoff[p7_NC1] == p7_CUTOFF_UNSET)
        ESL_FAIL(eslEINVAL, pli->errbuf, "NC bit thresholds unavailable on model %s\n", om->name);
      pli->T    = pli->incT    = om->cutoff[p7_NC1];
      pli->domT = pli->incdomT = om->cutoff[p7_NC2];
    }
  }

  return eslOK;
}

/* Function:  p7_pli_NewSeq()
 * Synopsis:  Prepare pipeline for a new sequence (target or query)
 *
 * Purpose:   Caller has a new sequence <sq>. Prepare the pipeline <pli>
 *            to receive this model as either a query or a target.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_pli_NewSeq(P7_PIPELINE *pli, const ESL_SQ *sq)
{

  pli->nseqs++;
  pli->nres += sq->n;
  if (pli->Z_setby == p7_ZSETBY_NTARGETS && pli->mode == p7_SEARCH_SEQS) pli->Z = pli->nseqs;
  return eslOK;
}

/* Function:  p7_pipeline_Merge()
 * Synopsis:  Merge the pipeline statistics
 *
 * Purpose:   Caller has a new model <om>. Prepare the pipeline <pli>
 *            to receive this model as either a query or a target.
 *
 *            The pipeline may alter the null model <bg> in a model-specific
 *            way (if we're using a composition bias filter HMM in the
 *            pipeline).
 *
 * Returns:   <eslOK> on success.
 * 
 *            <eslEINVAL> if pipeline expects to be able to use a
 *            model's bit score thresholds, but this model does not
 *            have the appropriate ones set.
 */
int
p7_pipeline_Merge(P7_PIPELINE *p1, P7_PIPELINE *p2)
{
  /* if we are searching a sequence database, we need to keep track of the
   * number of sequences and residues processed.
   */
  if (p1->mode == p7_SEARCH_SEQS)
    {
      p1->nseqs   += p2->nseqs;
      p1->nres    += p2->nres;
    }
  else
    {
      p1->nmodels += p2->nmodels;
      p1->nnodes  += p2->nnodes;
    }

  p1->n_past_msv  += p2->n_past_msv;
  p1->n_past_bias += p2->n_past_bias;
  p1->n_past_vit  += p2->n_past_vit;
  p1->n_past_fwd  += p2->n_past_fwd;
  p1->n_output    += p2->n_output;

  p1->pos_past_msv  += p2->pos_past_msv;
  p1->pos_past_bias += p2->pos_past_bias;
  p1->pos_past_vit  += p2->pos_past_vit;
  p1->pos_past_fwd  += p2->pos_past_fwd;
  p1->pos_output    += p2->pos_output;

  if (p1->Z_setby == p7_ZSETBY_NTARGETS) p1->Z += p2->Z;
  else                                   p1->Z =  p2->Z;

  return eslOK;
}

/* Function:  p7_pli_computeAliScores()
 * Synopsis:  Compute per-position scores for a BATH alignment 
 *
 * Purpose:   Compute per-position (Viterbi) scores for a BATH 
 *            alignment for use in splcing algorithms 
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_pli_computeAliScores_BATH(P7_DOMAIN *dom, P7_TRACE *tr, const ESL_SQ *seq, const P7_FS_PROFILE *gm_fs)
{

  int i, k, c, n;
  int codon_idx;
  int amino;
  int indel;
  int z1, z2;
  ESL_DSQ *nuc_dsq;
  int status;

  nuc_dsq   = seq->dsq;
  for(z1 = 0;       z1 < tr->N; z1++) if(tr->st[z1] == p7T_M) break;
  for(z2 = tr->N-1; z2 >= 0;    z2--) if(tr->st[z2] == p7T_M) break;

  dom->per_pos_len = z2 - z1 + 1;

  if(dom->scores_per_pos == NULL)
    ESL_ALLOC( dom->scores_per_pos, sizeof(float) * dom->per_pos_len);
  else
    ESL_REALLOC( dom->scores_per_pos, sizeof(float) * dom->per_pos_len);

  if(dom->k_per_pos == NULL)
    ESL_ALLOC( dom->k_per_pos, sizeof(int) * dom->per_pos_len);
  else
    ESL_REALLOC( dom->k_per_pos, sizeof(int) * dom->per_pos_len);

  
  for (n=0; n<dom->per_pos_len; n++)  dom->scores_per_pos[n] = 0.0;

  n  = 0;
  while (z1<=z2) {
    i = tr->i[z1];
    c = tr->c[z1];
    k = tr->k[z1]; 
     
    if (tr->st[z1] == p7T_M) {
      if(c == 1) {
        if(nuc_dsq[i] < p7P_MAXNUC)
          codon_idx = p7P_CODON1_FS5(nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN5_QC2;
        tr->fs++;
      }
      else if(c == 2) {
        if(nuc_dsq[i-1] < p7P_MAXNUC &&
         nuc_dsq[i] < p7P_MAXNUC)
          codon_idx = p7P_CODON2_FS5(nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN5_QC1;
        tr->fs++;
      }
      else if(c == 3) {
        if(nuc_dsq[i-2] < p7P_MAXNUC &&
         nuc_dsq[i-1] < p7P_MAXNUC &&
         nuc_dsq[i] < p7P_MAXNUC)
          codon_idx = p7P_CODON3_FS5(nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN5_C;
        
        indel = p7P_INDEL(gm_fs, k, codon_idx);
        /* record stop codon */
        if(indel == p7P_XXx || indel == p7P_XxX || indel == p7P_xXX) tr->fs++;
      }
      else if(c == 4) {
        if(nuc_dsq[i-3] < p7P_MAXNUC &&
         nuc_dsq[i-2] < p7P_MAXNUC &&
         nuc_dsq[i-1] < p7P_MAXNUC &&
         nuc_dsq[i] < p7P_MAXNUC)
          codon_idx = p7P_CODON4_FS5(nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN5_QC1;
        tr->fs++;
      }
      else if(c == 5) {
        if(nuc_dsq[i-4] < p7P_MAXNUC &&
         nuc_dsq[i-3] < p7P_MAXNUC &&
         nuc_dsq[i-2] < p7P_MAXNUC &&
         nuc_dsq[i-1] < p7P_MAXNUC &&
         nuc_dsq[i] < p7P_MAXNUC)
          codon_idx = p7P_CODON5_FS5(nuc_dsq[i-4], nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN5_QC2;
        tr->fs++;
      }

      amino = p7P_AMINO(gm_fs, k, codon_idx);
      if(gm_fs->codon_lengths == 5)
        dom->scores_per_pos[n] = p7P_MSC_AMINO5(gm_fs, k, amino);
      else if(gm_fs->codon_lengths == 3)
        dom->scores_per_pos[n] = p7P_MSC_AMINO3(gm_fs, k, amino);
      else if(gm_fs->codon_lengths == 1)
        dom->scores_per_pos[n] = p7P_MSC_AMINO1(gm_fs, k, amino);

      if (tr->st[z1-1] == p7T_I)
        dom->scores_per_pos[n] += p7P_TSC(gm_fs, k-1, p7P_IM);
      else if (tr->st[z1-1] == p7T_D)
        dom->scores_per_pos[n] += p7P_TSC(gm_fs, k-1, p7P_DM);

      dom->k_per_pos[n] = k;

      k++; z1++; n++;

      while(z1 < z2 && tr->st[z1] == p7T_M) {
        c = tr->c[z1];
        i = tr->i[z1];

        if(c == 1) {
          if(nuc_dsq[i] < p7P_MAXNUC)
            codon_idx = p7P_CODON1_FS5(nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN5_QC2;
          tr->fs++;
        }
        else if(c == 2) {
          if(nuc_dsq[i-1] < p7P_MAXNUC &&
           nuc_dsq[i] < p7P_MAXNUC)
            codon_idx = p7P_CODON2_FS5(nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN5_QC1;
          tr->fs++;
        }
        else if(c == 3) {
          if(nuc_dsq[i-2] < p7P_MAXNUC &&
           nuc_dsq[i-1] < p7P_MAXNUC &&
           nuc_dsq[i] < p7P_MAXNUC)
            codon_idx = p7P_CODON3_FS5(nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN5_C;
          indel = p7P_INDEL(gm_fs, k, codon_idx);
          if(indel == p7P_XXx || indel == p7P_XxX || indel == p7P_xXX) tr->fs++;
        }
        else if(c == 4) {
          if(nuc_dsq[i-3] < p7P_MAXNUC &&
           nuc_dsq[i-2] < p7P_MAXNUC &&
           nuc_dsq[i-1] < p7P_MAXNUC &&
           nuc_dsq[i] < p7P_MAXNUC)
            codon_idx = p7P_CODON4_FS5(nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN5_QC1;
          tr->fs++;
        }
        else if(c == 5) {
          if(nuc_dsq[i-4] < p7P_MAXNUC &&
           nuc_dsq[i-3] < p7P_MAXNUC &&
           nuc_dsq[i-2] < p7P_MAXNUC &&
           nuc_dsq[i-1] < p7P_MAXNUC &&
           nuc_dsq[i] < p7P_MAXNUC)
            codon_idx = p7P_CODON5_FS5(nuc_dsq[i-4], nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN5_QC2;
          tr->fs++;
        }

        amino = p7P_AMINO(gm_fs, k, codon_idx);
        if(gm_fs->codon_lengths == 5)
          dom->scores_per_pos[n] = p7P_MSC_AMINO5(gm_fs, k, amino);
        else if(gm_fs->codon_lengths == 3)
          dom->scores_per_pos[n] = p7P_MSC_AMINO3(gm_fs, k, amino);
        else if(gm_fs->codon_lengths == 1)
          dom->scores_per_pos[n] = p7P_MSC_AMINO1(gm_fs, k, amino);

        dom->scores_per_pos[n] += p7P_TSC(gm_fs, k-1, p7P_MM);
        dom->k_per_pos[n] = k;
        k++; z1++; n++;
      }
    }
    else if (tr->st[z1] == p7T_I) {
      
      dom->scores_per_pos[n] = p7P_TSC(gm_fs, k, p7P_MI);
      dom->k_per_pos[n] = k;
      z1++; n++;
      while (z1 < z2 && tr->st[z1] == p7T_I) {
        dom->scores_per_pos[n] = p7P_TSC(gm_fs, k, p7P_II);
        dom->k_per_pos[n] = k;
        z1++; n++;
      }
    }
    else if (tr->st[z1] == p7T_D) {
      dom->scores_per_pos[n] = p7P_TSC(gm_fs, k-1, p7P_MD);
      dom->k_per_pos[n] = k;
      k++; z1++; n++;
      while (z1 < z2 && tr->st[z1] == p7T_D)  {
        dom->scores_per_pos[n] = p7P_TSC(gm_fs, k-1, p7P_DD);
        dom->k_per_pos[n] = k;
        k++; z1++; n++;
      }
    }
    else ESL_XEXCEPTION(eslFAIL, "Impossible state from p7_pli_computeAliScores_BATH()");
  }

  dom->aliscore = 0.0;
  for (n=0; n<dom->per_pos_len; n++)  dom->aliscore += dom->scores_per_pos[n];

  return eslOK;

  ERROR:
    return status; 

}

/* Function:  p7_pli_postDomainDef_Frameshift_BATH()  
 * Synopsis:  the part of the BATH search Pipeline downstream
 *            of Domain Definition
 *
 * Purpose:   This is called by p7_pli_Frameshift(), and runs 
 *            the post-Domain Definition part of the frameshift 
 *            aware branch of the BATH pipeline. It consists of 
 *            running various bookkeeping and sanity checks on hits 
 *
 * Args:      pli             - the main pipeline object
 *            gm_fs5           - fs-aware codon profile (query)
 *            bg              - background model
 *            hitlist         - pointer to hit storage bin
 *            seqidx          - the id # of the DNA sequence  
 *            window_start    - the starting position of the DNA window
 *            dnasq           - the target dna sequence
 *            complementarity - boolean; is the passed window sourced from a complementary sequence block
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 */
static int 
p7_pli_postDomainDef_Frameshift_BATH(P7_PIPELINE *pli, P7_FS_PROFILE *gm_fs5, P7_BG *bg, P7_TOPHITS *hitlist, int64_t seqidx, int window_start, ESL_SQ *dnasq, ESL_SQ *windowsq, int complementarity)
{

  int              d;
  int              ali_len;
  int              env_len;
  int              status;
  int              tmp_i;
  float            bitscore;
  float            dom_bias; 
  float            dom_score;
  float            nullsc;
  double           dom_lnP;
  P7_DOMAIN       *dom        = NULL;      /* convenience variable, ptr to current domain  */
  P7_HIT          *hit        = NULL;      /* ptr to the current hit output data           */

  for (d = 0; d < pli->ddef->ndom; d++)
  {      
    dom = pli->ddef->dcl + d;
    
    ali_len = dom->jali - dom->iali + 1;
    bitscore = dom->envsc;

    if (ali_len < 12)   
    {// anything less than this is a funny byproduct of the Forward score passing a very low threshold, but no reliable alignment existing that supports it
      if(dom->scores_per_pos != NULL) free(dom->scores_per_pos);
      if(dom->k_per_pos != NULL)      free(dom->k_per_pos);
      p7_trace_fs_Destroy(dom->tr);
      continue;
    }
    tmp_i = dom->ienv;
    
    env_len = dom->jenv - dom->ienv + 1;
    /* map alignment and envelope coodinates to orignal DNA target sequence */
    if (!complementarity)
    { 
      dom->ienv       = dnasq->start + window_start + dom->ienv - 2;
      dom->jenv       = dnasq->start + window_start + dom->jenv - 2;
      dom->iali       = dnasq->start + window_start + dom->iali - 2;
      dom->jali       = dnasq->start + window_start + dom->jali - 2;
    }
    else
    {    
      dom->ienv       = dnasq->start - (window_start + dom->ienv) + 2;
      dom->jenv       = dnasq->start - (window_start + dom->jenv) + 2;
      dom->iali       = dnasq->start - (window_start + dom->iali) + 2;
      dom->jali       = dnasq->start - (window_start + dom->jali) + 2;  
    }

    /* Adjust score from env_len to max window length. Note that the loop and move 
     * costs are calculated based on amino lengths but paid per nucleotide*/
    bitscore = dom->envsc;
    bitscore -= 2 * log(2. / ((env_len/3.)+2));
    bitscore += 2 * log(2. / (gm_fs5->max_length+2));
    bitscore -= ((env_len-ali_len)/3.)                              * log((float) (env_len/3.) / (float) ((env_len/3.)+2));
    bitscore += ((ESL_MAX(env_len,gm_fs5->max_length*3)-ali_len)/3.) * log((float) gm_fs5->max_length / (float) (gm_fs5->max_length+2));

    /* Bias calculation and adjustments to Forward score */
    if (pli->do_null2)
      dom_bias = p7_FLogsum(0.0, log(bg->omega) + dom->domcorrection);
    else
      dom_bias = 0.0; 

    p7_bg_SetLength(bg, ESL_MAX(env_len/3,gm_fs5->max_length));
    p7_bg_fs_NullOne  (bg, dnasq->dsq, ESL_MAX(env_len/3,gm_fs5->max_length), &nullsc);
    dom_score  = (bitscore - (nullsc + dom_bias))  / eslCONST_LOG2;
     
    /* P-vaule calculation */
    dom_lnP   = esl_exp_logsurv(dom_score, gm_fs5->evparam[p7_FTAUFS5], gm_fs5->evparam[p7_FLAMBDA]);
     
    /* Check if hit passes the e-value cutoff based on the current
     * residue count. This prevents hits from accumulating and using
     * excessive memmory. */

    pli->Z = (float)pli->nres / (float)gm_fs5->max_length;
    if (pli->inc_by_E ? (exp(dom_lnP) * pli->Z <= pli->E) :  dom_score >= pli->T) 
    { 

      dom->ad = p7_alidisplay_fs_Create(dom->tr, 0, gm_fs5, windowsq, pli->show_cigar);
      dom->ad->sqfrom = dom->iali;
      dom->ad->sqto   = dom->jali;
      dom->ad->L      = dnasq->L;   
    
      p7_tophits_CreateNextHit(hitlist, &hit);
      hit->ndom        = 1;
      hit->best_domain = 0;

      hit->window_length = gm_fs5->max_length;
      hit->target_len = dnasq->n;

      hit->seqidx = seqidx;
      if(!complementarity)
        hit->subseq_start = dom->ienv - tmp_i + 1;
      else
        hit->subseq_start = dom->ienv + tmp_i - 1;

      ESL_ALLOC(hit->dcl, sizeof(P7_DOMAIN) );
      hit->dcl[0] = pli->ddef->dcl[d];
      
      hit->pre_score = bitscore  / eslCONST_LOG2;
      hit->pre_lnP   = esl_exp_logsurv (hit->pre_score,  gm_fs5->evparam[p7_FTAUFS5], gm_fs5->evparam[p7_FLAMBDA]);

      hit->dcl[0].dombias  = dom_bias;
      
      hit->sum_score  = hit->score  = hit->dcl[0].bitscore = dom_score;
      hit->sum_lnP    = hit->lnP    = hit->dcl[0].lnP  = dom_lnP;

      hit->sortkey    = pli->inc_by_E ? -dom_lnP : dom_score; // per-seq output sorts on bit score if inclusion is by score
    
      hit->frameshift = TRUE;
      if (pli->mode == p7_SEARCH_SEQS)
      {
        if (                       (status  = esl_strdup(dnasq->name, -1, &(hit->name)))  != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
        if (dnasq->acc[0]  != '\0' && (status  = esl_strdup(dnasq->acc,  -1, &(hit->acc)))   != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
        if (dnasq->desc[0] != '\0' && (status  = esl_strdup(dnasq->desc, -1, &(hit->desc)))  != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
      } else {
        if ((status  = esl_strdup(gm_fs5->name, -1, &(hit->name)))  != eslOK) esl_fatal("allocation failure");
        if ((status  = esl_strdup(gm_fs5->acc,  -1, &(hit->acc)))   != eslOK) esl_fatal("allocation failure");
        if ((status  = esl_strdup(gm_fs5->desc, -1, &(hit->desc)))  != eslOK) esl_fatal("allocation failure");
      }

    }
    else  
    {
      if(dom->scores_per_pos != NULL) free(dom->scores_per_pos);
      if(dom->k_per_pos != NULL) free(dom->k_per_pos);
      p7_trace_fs_Destroy(dom->tr);
    }
  }

  free(pli->ddef->dcl);
  pli->ddef->dcl = NULL;
  p7_domaindef_Reuse(pli->ddef);

  return eslOK;

ERROR:
  ESL_EXCEPTION(eslEMEM, "Error in Frameshift pipeline\n");

}

/* Function:  p7_pli_postDomainDef_BATH() 
 * Synopsis:  the part of the BATH search Pipeline downstream
 *            of Domain Definition
 *
 * Purpose:   This is called by p7_Pipeline_BATH(), or p7_pli_Frameshift() 
 *            and runs the post-Domain Definition part of non-frameshift 
 *            aware branch of the BATH pipeline. It consists of running 
 *            various bookkeeping and sanity checks on hits 
 *
 * Args:      pli             - the main pipeline object
 *            om              - optimized protien profile (query)  
 *            bg              - background model
 *            hitlist         - pointer to hit storage bin
 *            seqidx          - the id # of the DNA target sequence 
 *            window_start    - the starting position of the DNA window
 *            orfsq           - the target ORF  
 *            dnasq           - the target dna sequence
 *            complementarity - boolean; is the passed window sourced from a complementary sequence block
 *            nullsc          - domain nullsc
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 */

static int 
p7_pli_postDomainDef_BATH(P7_PIPELINE *pli, P7_OPROFILE *om, P7_BG *bg, P7_TOPHITS *hitlist, int64_t seqidx, int window_start, ESL_SQ *orfsq, ESL_SQ *dnasq, ESL_SQ *windowsq, int complementarity)
{

  int              d;
  int              ali_len, env_len;
  int              status;
  int              tmp_i;
  float            bitscore;
  float            dom_score;
  float            dom_bias;
  float            dom_lnP;
  float            nullsc;
  P7_DOMAIN       *dom        = NULL;      /* convenience variable, ptr to current domain  */
  P7_HIT          *hit        = NULL;      /* ptr to the current hit output data           */

  for (d = 0; d < pli->ddef->ndom; d++)
  {      
    dom = pli->ddef->dcl + d;
    
    env_len = dom->jenv - dom->ienv + 1;
    ali_len = (dom->jali - dom->iali + 1) / 3;   

    if (ali_len < 4) 
    {  // anything less than this is a funny byproduct of the Forward score passing a very low threshold, but no reliable alignment existing that supports it
      if(dom->scores_per_pos != NULL) free(dom->scores_per_pos);
      if(dom->k_per_pos != NULL) free(dom->k_per_pos); 
      p7_trace_fs_Destroy(dom->tr);
      continue; 
    }
    tmp_i = dom->ienv;

    /* map alignment and envelope coodinates to orignal DNA target sequence */ 
    if (!complementarity)
    { 
      dom->ienv       = dnasq->start + orfsq->start + dom->ienv*3 - 4; //minus an extra 2 to get to start of codon 
      dom->jenv       = dnasq->start + orfsq->start + dom->jenv*3 - 2;
     
      dom->iali       = dnasq->start + window_start + dom->iali - 2;
      dom->jali       = dnasq->start + window_start + dom->jali - 2;
    }
    else
    {
      dom->ienv       = dnasq->end + orfsq->start - dom->ienv*3 + 2;  //plus 2 to get to end of codon
      dom->jenv       = dnasq->end + orfsq->start - dom->jenv*3; 

      dom->jali       = dnasq->start - (window_start + dom->jali) + 2;
      dom->iali       = dnasq->start - (window_start + dom->iali) + 2; 
    }


    /* Adjust score from env_len to max window length */ 
    bitscore = dom->envsc;
    bitscore -= 2 * log(2. / (env_len+2));
    bitscore += 2 * log(2. / (om->max_length+2));
    bitscore -= (env_len-ali_len)          * log((float) env_len / (float) (env_len+2));  
    bitscore += (om->max_length-ali_len) * log((float) om->max_length / (float) (om->max_length+2)); 
   
    /* Bias calculation and adjustments to Forward score */
    if (pli->do_null2)
      dom_bias = p7_FLogsum(0.0, log(bg->omega) + dom->domcorrection);
    else
      dom_bias = 0.0;

 
     p7_bg_SetLength(bg, om->max_length);
     p7_bg_NullOne  (bg, orfsq->dsq, om->max_length, &nullsc);
     dom_score =  (bitscore - (nullsc + dom_bias)) / eslCONST_LOG2;
     
     /* p-value calculations */
     dom_lnP   = esl_exp_logsurv(dom_score, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
   
     /* Check if hit passes the e-value cutoff based on the current residue count.
     * This prevents unreportable hits from accumulating and using excessive memmory.
     * For spliced alignment also keep all hits with a final P-value below the MSV cuttoff. */
     pli->Z = (float)pli->nres / (float)om->max_length;
     if ((pli->spliced && ((pli->inc_by_E ? (exp(dom_lnP) * pli->Z <= pli->E) :  dom_score >= pli->T) || exp(dom_lnP) < pli->F3)) ||
        (!pli->spliced &&  (pli->inc_by_E ? (exp(dom_lnP) * pli->Z <= pli->E) :  dom_score >= pli->T)))
     {
      
       dom->ad = p7_alidisplay_nonfs_Create(dom->tr, 0, om, windowsq, orfsq, dom->tr->sqfrom[0], pli->show_cigar);
       dom->ad->exon_cnt = 1;
       dom->ad->sqfrom   = dom->iali;
       dom->ad->sqto     = dom->jali;
       dom->ad->L        = dnasq->L; 
        
       /* Add hits to hitlist and check if they are reprotable*/   
       p7_tophits_CreateNextHit(hitlist, &hit);
    
       hit->ndom        = 1;
       hit->best_domain = 0;
       hit->window_length = orfsq->n;
       hit->target_len = dnasq->n;
       hit->seqidx = seqidx;
       if(!complementarity)
         hit->subseq_start = dom->ienv - (orfsq->start - windowsq->start + (tmp_i*3)) + 3;
       else
          hit->subseq_start = dom->ienv + (dnasq->n - orfsq->start + 1) - windowsq->start + (tmp_i*3) - 3;

       ESL_ALLOC(hit->dcl, sizeof(P7_DOMAIN) );
       hit->dcl[0] = pli->ddef->dcl[d];

       hit->pre_score = bitscore  / eslCONST_LOG2;
       hit->pre_lnP   = esl_exp_logsurv (hit->pre_score,  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
       
       hit->dcl[0].dombias  = dom_bias;
   
       hit->sum_score  = hit->score  = hit->dcl[0].bitscore = dom_score;
       hit->sum_lnP    = hit->lnP    = hit->dcl[0].lnP  = dom_lnP;

       hit->sortkey    = pli->inc_by_E ? -dom_lnP : dom_score; // per-seq output sorts on bit score if inclusion is by score

       if (                       (status  = esl_strdup(dnasq->name, -1, &(hit->name)))  != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
       if (dnasq->acc[0]  != '\0' && (status  = esl_strdup(dnasq->acc,  -1, &(hit->acc)))   != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
       if (dnasq->desc[0] != '\0' && (status  = esl_strdup(dnasq->desc, -1, &(hit->desc)))  != eslOK) ESL_EXCEPTION(eslEMEM, "allocation failure");
    } 
    else { //delete unused P7_ALIDSPLAY
      if(dom->scores_per_pos != NULL) free(dom->scores_per_pos);
      if(dom->k_per_pos != NULL) free(dom->k_per_pos);
      p7_trace_fs_Destroy(dom->tr);
    }
  }

  free(pli->ddef->dcl);
  pli->ddef->dcl = NULL;
  p7_domaindef_Reuse(pli->ddef);

return eslOK;

ERROR:
  ESL_EXCEPTION(eslEMEM, "Error in nonFrameshift pipeline\n");
}



/* Function:  p7_pli_Frameshift()
 * Synopsis:  the part of the BATH search Pipeline downstream
 *            of p7_ForwardParser() for framshift search
 *
 * Purpose:   This is called by p7_Pipeline_BATH(), and starts the
 *            frameshift aware pipeline. It consists of creating 
 *            DNA windows for all the ORFs that passed the standard
 *            Forward filter, running the frameshift aware
 *            Forward filter on those windows, and comparring the
 *            standard andframeshift aware Forward p-values. The
 *            lower p-value will dctate which pipeline, standard 
 *            translation or frameshift aware, will be used for the
 *            remained of the pipeline. 
 *
 * Args:      pli             - the main pipeline object
 *            om              - optimized protien profile 
 *            om_fs3          - optimized 3 codon length frameshift profile
 *            om_fs5          - optimized 5 codon length frameshift profile
 *            gm_fs5          - non-optimized 5 codon length frameshift profile 
 *            bg              - background model
 *            hitlist         - pointer to hit storage bin
 *            seqidx          - the id # of the DNA sequence
 *            orf_block       - collection of ORFs translated form <dnasq>
 *            dnasq           - the target DNA sequence
 *            gcode           - genetic code information for codon translation
 *            pli_tmp         - frameshift pipeline object for temporary data
 *            hit_windows     - ORF ungapped aligment windows for DNA window building 
 *            complementarity - boolean; is the passed window sourced from a complementary sequence block
 *
 * Returns:   <eslOK>  
 *
 */
static int
p7_pli_Frameshift(P7_PIPELINE *pli, P7_OPROFILE *om, P7_FS_OPROFILE *om_fs3, P7_FS_OPROFILE *om_fs5, P7_FS_PROFILE *gm_fs5, P7_SCOREDATA *data, P7_BG *bg, P7_TOPHITS *hitlist, int64_t seqidx, ESL_SQ_BLOCK *orf_block, ESL_SQ *dnasq, ESL_GENCODE *gcode, P7_PIPELINE_OBJS *pli_tmp, P7_HMM_WINDOWLIST *hit_windows, int complementarity)
{

  int              i, w, h;
  int              status;
  int              k_min, k_max;
  int              last_window_cnt;
  int              orf_cnt;
  int64_t          orf_start, orf_end;
  int64_t          window_start, window_end;
  float            fwdsc;                      /* framshift forward scores                               */
  float            nullsc;                     /* null score for DNA window            */
  float            filtersc;                   /* global filtersc for DNA window          */
  float            local_filtersc;             /* local filtersc for DNA window          */
  float            seqscore;                   /* the corrected per-seq bit score              */
  float            tot_orfsc;                  /* summed score for all ORFs in current DNA window */
  float            local_compo[p7_MAXABET];    /* Local model composition from windows    */
  double           P_fs;                       /* P-value of frameshift forward for window*/
  double           P_null;                     /* P-value of frameshift forward for window w/o bias adjustments*/
  double           P_tot;                      /* P-value of summed forward score for all ORFs */
  double           P_min;                      /* lowest p-value produced by an ORF */
  ESL_SQ          *orfsq;
  P7_HMM_WINDOWLIST fwd_windowlist;
  P7_HMM_WINDOW     *dna_window;

  p7_hmmwindow_init(&fwd_windowlist);

  /* Build windows from ORF's that pass F4 */
  p7_pli_BuildDNAWindows(pli, orf_block, dnasq, om, bg, data, &fwd_windowlist, 0., pli_tmp, hit_windows, complementarity);
  
  for(w = 0; w < fwd_windowlist.count; w++) {

    dna_window = &(fwd_windowlist.windows[w]); 

    window_start = complementarity ? dnasq->start - (dna_window->n + dna_window->length) : dnasq->start + dna_window->n - 1; 
    window_end   = complementarity ? dnasq->start - dna_window->n + 1 : window_start + dna_window->length - 1;

    pli_tmp->tmpseq->L     = dna_window->length;
    pli_tmp->tmpseq->n     = dna_window->length;
    pli_tmp->tmpseq->start = dna_window->n;
    pli_tmp->tmpseq->end   = dna_window->n + dna_window->length - 1; 
    pli_tmp->tmpseq->dsq   = dnasq->dsq + dna_window->n - 1;;
  
    orf_cnt   = 0;
    tot_orfsc = -eslINFINITY;
    P_tot     = eslINFINITY;
    P_min     = eslINFINITY;
 
    k_min = om->M;
    k_max = 0;
    last_window_cnt = 0;

    /* Get ORF P values for comparision */ 
    for(i = 0; i < orf_block->count; i++) {
      if(pli_tmp->P_orf[i] > pli->F4) continue;

      orfsq = &(orf_block->list[i]);

      if(complementarity) {
        orf_start =  dnasq->start - (dnasq->n - orfsq->end   + 1) + 1;
        orf_end   =  dnasq->start - (dnasq->n - orfsq->start + 1) + 1;
      } else {    
        orf_start = dnasq->start + orfsq->start - 1;
        orf_end   = dnasq->start + orfsq->end   - 1;
      } 

      /* Only process ORF if it in inside the current window */ 
      if(orf_start >= window_start && orf_end <= window_end) {    
        orfsq->idx = w;
        P_min      = ESL_MIN(P_min, pli_tmp->P_orf[i]);  
        tot_orfsc  = p7_FLogsum(tot_orfsc, pli_tmp->fwdsc[i]); 
        orf_cnt++;        

        h = last_window_cnt; 
        while(h < hit_windows->count && hit_windows->windows[h].id != i) h++;
        if (h < hit_windows->count) {
          while(h < hit_windows->count && hit_windows->windows[h].id == i) {
            k_min = ESL_MIN(k_min, hit_windows->windows[h].k - hit_windows->windows[h].length + 1);
            k_max = ESL_MAX(k_max, hit_windows->windows[h].k);
            h++;
          }
          last_window_cnt = h;
        }
      }
    }

    
    P_tot = esl_exp_surv(tot_orfsc / eslCONST_LOG2,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);

    /* Run Frameshift Forward on window */
    p7_bg_SetLength(bg, dna_window->length/3);
    p7_bg_fs_NullOne(bg, pli_tmp->tmpseq->dsq, dna_window->length/3, &nullsc);
        
    if (pli->do_biasfilter) {
      p7_bg_fs_FilterScore(bg, pli_tmp->tmpseq->dsq, pli_tmp->tmpseq->n, gcode, &filtersc);
      if(k_min <= k_max) { // We have ORf windows for this DNA window 
        p7_pli_ComputeLocalCompo(data, om, bg, k_min, k_max, local_compo); 
        p7_bg_SetFilter(bg, om->M, local_compo);
        p7_bg_SetLength(bg, dna_window->length/3);
        p7_bg_fs_FilterScore(bg, pli_tmp->tmpseq->dsq, pli_tmp->tmpseq->n, gcode, &local_filtersc);
        if(local_filtersc > filtersc) filtersc = local_filtersc;
        p7_bg_SetFilter(bg, om->M, om->compo);
        p7_bg_SetLength(bg, dna_window->length/3);
      }
    }
    else filtersc = nullsc;

    p7_omx_GrowTo_dpf(pli->oxf_fs, om->M, PARSER_ROWS_FWD, dna_window->length);
    p7_oivx_GrowTo(pli->ov3, om_fs3->M, p7P_3CODONS);
    p7_fs_oprofile_ReconfigLength(om_fs3, dna_window->length/3);    

    status = p7_ForwardParser_Frameshift_3Codons(pli_tmp->tmpseq->dsq, dna_window->length, om_fs3, pli->oxf_fs, pli->ov3, &fwdsc);

    seqscore = (fwdsc-filtersc) / eslCONST_LOG2;
    P_fs = esl_exp_surv(seqscore,  om_fs3->evparam[p7_FTAUFS3],  om_fs3->evparam[p7_FLAMBDA]);
    P_null = esl_exp_surv((fwdsc-nullsc)/eslCONST_LOG2,  om_fs3->evparam[p7_FTAUFS3],  om_fs3->evparam[p7_FLAMBDA]);

   
    if(!pli->std_pipe) P_tot = 1.0; // for --fsonly
  
    /* Compare window P value to ORF P-values */ 
    /* Use frameshift pip;ine uf frameshift P-value meets threshold ( <= pli->F3 ) AND one of the following constions is met
       1. frameshift P-value (no bais adjustment) is less than the summed ORF P-value
       2. frameshift P-value (no bais adjustment) is eq to the summed ORF P-value and there are 2 or more ORFs for this window
       3. The best ORF P-value does not meet threshold. 
    */
    if(P_fs <= pli->F3 && (P_null < P_tot || (P_null == P_tot && orf_cnt > 1) || P_min > pli->F3)) { 
      
      /* Window p-value won - continue frameshift pipeline */  
      pli->pos_past_fwd += dna_window->length; 
      p7_omx_GrowTo_dpf(pli->oxb_fs, om->M, PARSER_ROWS_BWD, dna_window->length);
      status = p7_BackwardParser_Frameshift_3Codons(pli_tmp->tmpseq->dsq, dna_window->length, om_fs3, pli->oxf_fs, pli->oxb_fs, pli->ov3, NULL);
      if (status == eslERANGE) continue; /* backward underflow; skip domain definition for this window */
  
      status = p7_domaindef_ByPosteriorHeuristics_Frameshift_BATH(pli, pli_tmp->tmpseq, om_fs5, gm_fs5, bg, gcode);
      if (pli->ddef->nregions == 0 || pli->ddef->nenvelopes == 0) continue;
  
      /* Send any hits from the Frameshift aware pipeline to be further processed */ 
      p7_pli_postDomainDef_Frameshift_BATH(pli, gm_fs5, bg, hitlist, seqidx, dna_window->n, dnasq, pli_tmp->tmpseq, complementarity);
  
    } 
    else if (pli->std_pipe) { 
    
      for(i = 0; i < orf_block->count; i++) {
       
        orfsq = &(orf_block->list[i]);

        if(orfsq->idx != w)                continue; // This ORF does not overlap with this window
        if(pli_tmp->P_orf[i] > pli->F3)    continue; // This ORF did not pass Forward
        if(pli_tmp->oxf_holder[i] == NULL) continue; // This ORF has already been aligned
        
        pli->pos_past_fwd += orfsq->n * 3;

        p7_oprofile_ReconfigLength(om, orfsq->n);
        p7_omx_GrowTo(pli->oxb, om->M, 0, orfsq->n);     
          
        p7_BackwardParser(orfsq->dsq, orfsq->n, om, pli_tmp->oxf_holder[i], pli->oxb, NULL);
          
        status = p7_domaindef_ByPosteriorHeuristics_BATH(orfsq, pli_tmp->tmpseq, dnasq->n, om, gm_fs5, pli_tmp->oxf_holder[i], pli->oxb, pli->fwd, pli->bck, pli->ddef);
  
        if (status != eslOK) ESL_FAIL(status, pli->errbuf, "domain definition workflow failure"); /* eslERANGE can happen */
        if (pli->ddef->nregions == 0 || pli->ddef->nenvelopes == 0) {
          p7_omx_Destroy(pli_tmp->oxf_holder[i]);
          pli_tmp->oxf_holder[i] = NULL;
          continue;
        }        
  
        /* Send any hits from the standard pipeline to be further processed */   
        p7_pli_postDomainDef_BATH(pli, om, bg, hitlist, seqidx, dna_window->n, orfsq, dnasq, pli_tmp->tmpseq, complementarity);

        p7_omx_Destroy(pli_tmp->oxf_holder[i]);
        pli_tmp->oxf_holder[i] = NULL;
      }  
    } 
  }

  for(i = 0; i < orf_block->count; i++) { 
    p7_omx_Destroy(pli_tmp->oxf_holder[i]);
    pli_tmp->oxf_holder[i] = NULL;
  }
  if (fwd_windowlist.windows != NULL) free (fwd_windowlist.windows); 
    
  return eslOK;
}




/* Function:  p7_Pipeline_BATH()
 * Synopsis:  Sequence to profile comparison pipeline for 
 *            frameshift aware translated search - bathsearch.
 *
 * Purpose:   Run translated search pipeline to compare a protien 
 *            profile <gm/om> against a DNA sequence <sq>. For the 
 *            first stages of the pipeline (MSV, bias and viterbi 
 *            filters) each DNA strand is translated into ORFs in 
 *            all 3 frames and these are compared directly to an 
 *            optimized protien profile <om>. For the forward filter
 *            both an ORF to <om> and a DNA window to frameshift 
 *            aware codon model <gm_fs5> and a comparison is preformed. 
 *            Which ever Forward filter produces the lower p-value  
 *            determines which target and query form are used for the 
 *            remainder of the pipeline. If a significant hit is 
 *            found, information about it is added to the <hitlist>. 
 *            The pipeline accumulates bean-counting information about 
 *            how many comparisons flow through the pipeline while it's 
 *            active.
 *            
 * Returns:   <eslOK> on success. If a significant hit is obtained,
 *            its information is added to the growing <hitlist>. 
 *            
 *            <eslEINVAL> if (in a scan pipeline) we're supposed to
 *            set GA/TC/NC bit score thresholds but the model doesn't
 *            have any.
 *            
 *            <eslERANGE> on numerical overflow errors in the
 *            optimized vector implementations; particularly in
 *            posterior decoding. I don't believe this is possible for
 *            multihit local models, but I'm set up to catch it
 *            anyway. We may emit a warning to the user, but cleanly
 *            skip the problematic sequence and continue.
 * 
 * Args:      pli             - the main pipeline object
 *            om              - optimized protein profile (query)
 *            gm              - generic protein profile (query)
 *            gm_fs5           - generic fs-aware codon profile (query)
 *            data            - for picking window edges based on 
 *                              maximum prefix/suffix extensions
 *            bg              - background model
 *            hitlist         - pointer to hit storage bin (already 
 *                              allocated)
 *            seqidx          - the id # of the sequence from which 
 *                              the current window was extracted
 *            dnasq           - digital sequence of the DNA window
 *            orf_block       - collection of ORFs translated form <dnasq>
 *            gcode           - genetic code information for codon translation
 *            complementarity - is <sq> from the top strand 
 *                        (p7_NOCOMPLEMENT), or bottom strand 
 *                        (P7_COMPLEMENT)
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 * Xref:      J4/25.
 */
int
p7_Pipeline_BATH(P7_PIPELINE *pli, P7_OPROFILE *om, P7_FS_OPROFILE *om_fs3, P7_FS_OPROFILE *om_fs5, P7_FS_PROFILE *gm_fs5, P7_SCOREDATA *data, P7_BG *bg, P7_TOPHITS *hitlist, int64_t seqidx, ESL_SQ *dnasq, ESL_SQ_BLOCK *orf_block, ESL_GENCODE *gcode,P7_HMM_WINDOWLIST *hit_windows, int complementarity)
{

  int     i, w;
  int     status;
  int     k_min, k_max;
  int     old_window_cnt;
  int64_t orf_start, orf_end; 
  float   local_compo[p7_MAXABET]; /* Local model composition from windows    */
  float   nullsc;                  /* null model score                        */
  float   usc;                     /* msv score                               */
  float   vfsc;                    /* viterbi score                           */
  float   fwdsc;                   /* forward score                           */
  float   seqsc;                   /* null corrected bit score                */
  float   filtersc;                /* bias and null score                     */
  float   local_filtersc;          /* bias and null score using local_compo   */
  double  P;                       /* p-value holder                          */
  ESL_SQ           *orfsq;         /* ORF sequence                            */
  P7_HMM_WINDOW    *window;             
  P7_PIPELINE_OBJS *pli_tmp;   

  if (dnasq->n < 15) return eslOK;         //DNA to short
  if (orf_block->count == 0) return eslOK; //No ORFS translated

  pli_tmp = NULL;
  ESL_ALLOC(pli_tmp, sizeof(P7_PIPELINE_OBJS));
  pli_tmp->tmpseq     = NULL;
  pli_tmp->oxf_holder = NULL;
  pli_tmp->P_orf      = NULL;
  pli_tmp->fwdsc      = NULL;

  ESL_ALLOC(pli_tmp->fwdsc,      sizeof(float)    * orf_block->count);
  ESL_ALLOC(pli_tmp->P_orf,      sizeof(double)   * orf_block->count);
  ESL_ALLOC(pli_tmp->oxf_holder, sizeof(P7_OMX *) * orf_block->count);

  for(i = 0; i < orf_block->count; i++) {
      pli_tmp->oxf_holder[i] = NULL;
      pli_tmp->fwdsc[i] = -eslINFINITY;
      pli_tmp->P_orf[i] = 1.0;
  }
  
  pli_tmp->tmpseq = esl_sq_CreateDigital(dnasq->abc);
  free (pli_tmp->tmpseq->dsq); 
  if ((status = esl_sq_SetName     (pli_tmp->tmpseq, dnasq->name))   != eslOK) goto ERROR;
  if ((status = esl_sq_SetSource   (pli_tmp->tmpseq, dnasq->source)) != eslOK) goto ERROR;
  if ((status = esl_sq_SetAccession(pli_tmp->tmpseq, dnasq->acc))    != eslOK) goto ERROR;
  if ((status = esl_sq_SetDesc     (pli_tmp->tmpseq, dnasq->desc))   != eslOK) goto ERROR;

  for (i = 0; i < orf_block->count; ++i)
  { 
    orfsq = &(orf_block->list[i]);
    if (   (orfsq->start < orfsq->end    &&  orfsq->end < dnasq->C )  ||
           (orfsq->end < orfsq->start    &&  orfsq->start < dnasq->C ) )
        continue; /* don't bother with an orf that showed up completely within a previous window */

    if(orfsq->n > 0) 
    {
      vfsc = -eslINFINITY;

      p7_bg_SetLength(bg, orfsq->n);
      p7_oprofile_ReconfigLength(om, orfsq->n);
      p7_bg_NullOne  (bg, orfsq->dsq, orfsq->n, &nullsc);

      p7_omx_GrowTo(pli->oxf, om->M, 0, orfsq->n);    /* expand the one-row omx if needed */
      
      p7_MSVFilter(orfsq->dsq, orfsq->n, om, pli->oxf, &usc);
      seqsc = (usc - nullsc) / eslCONST_LOG2;
      P = esl_gumbel_surv( seqsc,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
      if (P > pli->F1 ) continue;
    
      pli->pos_past_msv  += orfsq->n * 3; 
      
      /* biased composition HMM filtering */
      if (pli->do_biasfilter)
      {
        p7_bg_FilterScore(bg, orfsq->dsq, orfsq->n, &filtersc);
        seqsc = (usc - filtersc) / eslCONST_LOG2;
        P = esl_gumbel_surv(seqsc,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
        if (P > pli->F1) continue;
      }  else filtersc = nullsc;

      pli->pos_past_bias += orfsq->n * 3;    

      old_window_cnt = hit_windows->count;
      /* Viterbi filer on ORF */
      if (P > pli->F2)
      {
        p7_ViterbiFilter_BATH(orfsq->dsq, orfsq->n, om, pli->oxf, data, filtersc, pli->F2, hit_windows, &vfsc);
        seqsc = (vfsc-filtersc) / eslCONST_LOG2;
        P  = esl_gumbel_surv(seqsc,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
        if (P > pli->F2) { hit_windows->count = old_window_cnt; continue; }
      }
      else  // get windows
        p7_SSVFilter_BATH(orfsq->dsq, orfsq->n, om, pli->oxf, data, bg, pli->F1, hit_windows);

      for(w = old_window_cnt; w < hit_windows->count; w++) 
        hit_windows->windows[w].id = i; // record ORF index for DNA windows
      
      pli->pos_past_vit  += orfsq->n * 3;      

      if (pli->do_biasfilter && old_window_cnt < hit_windows->count) {
        /* Find the min and max hmm postions for all the windows from this ORF  */
        k_max  = hit_windows->windows[old_window_cnt].k;
        k_min  = k_max - hit_windows->windows[old_window_cnt].length + 1;
        for (w = old_window_cnt + 1; w < hit_windows->count; w++) {
          k_max = ESL_MAX(k_max,   hit_windows->windows[w].k);
          k_min = ESL_MIN(k_min, hit_windows->windows[w].k - hit_windows->windows[w].length + 1);
        }

        p7_pli_ComputeLocalCompo(data, om, bg, k_min, k_max, local_compo);
        p7_bg_SetFilter(bg, om->M, local_compo);
        p7_bg_SetLength(bg, orfsq->n);
        p7_bg_FilterScore(bg, orfsq->dsq, orfsq->n, &local_filtersc);
    
        if(local_filtersc > filtersc) {
          filtersc = local_filtersc;
          if(vfsc == -eslINFINITY) { // Viterbi not run 
            seqsc = (usc - filtersc) / eslCONST_LOG2;
            P = esl_gumbel_surv(seqsc, om->evparam[p7_MMU], om->evparam[p7_MLAMBDA]);
            if (P > pli->F2) {
              p7_ViterbiFilter(orfsq->dsq, orfsq->n, om, pli->oxf, &vfsc);
              seqsc = (vfsc - filtersc) / eslCONST_LOG2;
              P = esl_gumbel_surv(seqsc, om->evparam[p7_VMU], om->evparam[p7_VLAMBDA]);
              if (P > pli->F2) { hit_windows->count = old_window_cnt; continue; }
            }
          } 
          else {
            seqsc = (vfsc - filtersc) / eslCONST_LOG2;
            P = esl_gumbel_surv(seqsc, om->evparam[p7_VMU], om->evparam[p7_VLAMBDA]);
            if (P > pli->F2) { hit_windows->count = old_window_cnt; continue; }
          }
        }
        p7_bg_SetFilter(bg, om->M, om->compo);
        p7_bg_SetLength(bg, orfsq->n);
      } 

      if(!pli->fs_pipe) {
        if(pli->spliced) {
          for(w = old_window_cnt; w < hit_windows->count; w++) {
            window = &(hit_windows->windows[w]);

            window->id = seqidx;
            window->complementarity = complementarity;
            if(complementarity)
              window->n = dnasq->end + orfsq->start - ((window->n+window->length-1)*3);
            else
              window->n = dnasq->start + orfsq->start + (window->n*3) - 4;
            window->length *= 3;
          }
        }

        p7_ForwardParser(orfsq->dsq, orfsq->n, om, pli->oxf, &fwdsc);
        seqsc = (fwdsc-filtersc) / eslCONST_LOG2;
        P = esl_exp_surv(seqsc,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
        if (P > pli->F3) continue;

        if(pli->spliced) {
          for(w = old_window_cnt; w < hit_windows->count; w++) 
            hit_windows->windows[w].pass_forward = TRUE;     
        }

        /*Get ORF DNA seqeunce for alignment */
        if(complementarity) {
          orf_start = dnasq->n - orfsq->start + 1;
          orf_end   = dnasq->n - orfsq->end + 1;
        }
        else {
          orf_start = orfsq->start;
          orf_end   = orfsq->end;
        }

        pli_tmp->tmpseq->start = orf_start;
        pli_tmp->tmpseq->end   = orf_end;
        pli_tmp->tmpseq->dsq   = dnasq->dsq + orf_start - 1;
        pli_tmp->tmpseq->n     = orf_end - orf_start + 1;
        pli_tmp->tmpseq->L     = pli_tmp->tmpseq->n;

        pli->pos_past_fwd += orfsq->n * 3;
        p7_omx_GrowTo(pli->oxb, om->M, 0, orfsq->n);

        p7_BackwardParser(orfsq->dsq, orfsq->n, om, pli->oxf, pli->oxb, NULL);

        status = p7_domaindef_ByPosteriorHeuristics_BATH(orfsq, pli_tmp->tmpseq, dnasq->n, om, gm_fs5, pli->oxf, pli->oxb, pli->fwd, pli->bck, pli->ddef);
        if (status != eslOK) ESL_FAIL(status, pli->errbuf, "domain definition workflow failure"); /* eslERANGE can happen */
        if (pli->ddef->nregions   == 0)  continue; /* score passed threshold but there's no discrete domains here     */
        if (pli->ddef->nenvelopes == 0)  continue; /* rarer: region was found, stochastic clustered, no envelope found*/

        p7_pli_postDomainDef_BATH(pli, om, bg, hitlist, seqidx, orf_start, orfsq, dnasq, pli_tmp->tmpseq, complementarity);  
      }
      /* Frameshift pipeline F4 filter */
      else { 

        /* Save the Forward Martix for each ORF so we do not have to rerun
         * Forward in the event that the standard pipeline is selected */
        if ((pli_tmp->oxf_holder[i] = p7_omx_Create_dpf(om->M, 0, orfsq->n, p7G_NSCELLS)) == NULL) goto ERROR;
        p7_ForwardParser(orfsq->dsq, orfsq->n, om, pli_tmp->oxf_holder[i], &fwdsc);    
 
        seqsc = (fwdsc-filtersc) / eslCONST_LOG2;
        pli_tmp->P_orf[i] = esl_exp_surv(seqsc, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
        pli_tmp->fwdsc[i] = fwdsc-nullsc; 

        if(pli_tmp->P_orf[i] > pli->F4) {
          p7_omx_Destroy(pli_tmp->oxf_holder[i]);
          pli_tmp->oxf_holder[i] = NULL;
        }
      }
    }
  }

  if(pli->fs_pipe)  { 
    p7_pli_Frameshift(pli, om, om_fs3, om_fs5, gm_fs5, data, bg, hitlist, seqidx, orf_block, dnasq, gcode, pli_tmp, hit_windows, complementarity);
  }

  /* clean up */ 
  if (pli_tmp != NULL) 
  {
    pli_tmp->tmpseq->dsq = NULL;
    if (pli_tmp->tmpseq     != NULL) esl_sq_Destroy(pli_tmp->tmpseq);
    if (pli_tmp->oxf_holder != NULL) free(pli_tmp->oxf_holder);
    if (pli_tmp->P_orf      != NULL) free(pli_tmp->P_orf);
    if (pli_tmp->fwdsc      != NULL) free(pli_tmp->fwdsc);
    free(pli_tmp);
  }

  return eslOK;

ERROR:
  if (pli_tmp != NULL)
  {
    if (pli_tmp->tmpseq     != NULL) esl_sq_Destroy(pli_tmp->tmpseq);
    if (pli_tmp->oxf_holder != NULL) free(pli_tmp->oxf_holder);
    if (pli_tmp->P_orf      != NULL) free(pli_tmp->P_orf);
    if (pli_tmp->fwdsc      != NULL) free(pli_tmp->fwdsc);
    free(pli_tmp);
  }

  return status;
}

/* Function:  p7_pli_Statistics()
 * Synopsis:  Final statistics output from a processing pipeline.
 *
 * Purpose:   Print a standardized report of the internal statistics of
 *            a finished processing pipeline <pli> to stream <ofp>.
 *            
 *            If stopped, non-<NULL> stopwatch <w> is provided for a
 *            stopwatch that was timing the pipeline, then the report
 *            includes timing information.
 *
 * Returns:   <eslOK> on success.
 */
  int
p7_pli_Statistics(FILE *ofp, P7_PIPELINE *pli, ESL_STOPWATCH *w)
{

  fprintf(ofp, "Internal pipeline statistics summary:\n");
  fprintf(ofp, "-------------------------------------\n");
  if (pli->mode == p7_SEARCH_SEQS) {
    fprintf(ofp, "Query model(s):              %15" PRId64 "  (%" PRId64 " nodes)\n",     pli->nmodels, pli->nnodes);
    fprintf(ofp, "Target %-12s          %15" PRId64 "  (%" PRId64 " residues searched)\n", "sequence(s):", pli->nseqs,   pli->nres);
    
  } else {
    fprintf(ofp, "Query %-12s           %15" PRId64 "  (%" PRId64 " residues searched)\n", "sequence(s):", pli->nseqs,   pli->nres);
    fprintf(ofp, "Target model(s):             %15" PRId64 "  (%" PRId64 " nodes)\n",     pli->nmodels, pli->nnodes);
    
  }

  fprintf(ofp, "Residues passing SSV filter: %15" PRId64 "  (%.3g); expected (%.3g)\n",
      pli->pos_past_msv,
      (double)pli->pos_past_msv / (pli->nres*pli->nmodels) ,
      pli->F1);

  fprintf(ofp, "Residues passing bias filter:%15" PRId64 "  (%.3g); expected (%.3g)\n",
      pli->pos_past_bias,
      (double)pli->pos_past_bias / (pli->nres*pli->nmodels) ,
      pli->F1);

  fprintf(ofp, "Residues passing Vit filter: %15" PRId64 "  (%.3g); expected (%.3g)\n",
      pli->pos_past_vit,
      (double)pli->pos_past_vit / (pli->nres*pli->nmodels) ,
      pli->F2);

  fprintf(ofp, "Residues passing Fwd filter: %15" PRId64 "  (%.3g); expected (%.3g)\n",
      pli->pos_past_fwd,
      (double)pli->pos_past_fwd / (pli->nres*pli->nmodels) ,
      pli->F3);

  fprintf(ofp, "Total number of hits:        %15d  (%.3g)\n",
      (int)pli->n_output,
      (double)pli->pos_output / (pli->nres*pli->nmodels) );

  if (w != NULL) {
    esl_stopwatch_Display(ofp, w, "# CPU time: ");
    fprintf(ofp, "# Mc/sec: %.2f\n", 
        (double) pli->nres * (double) pli->nnodes / (w->elapsed * 1.0e6));
  }

  return eslOK;
}
/*------------------- end, pipeline API -------------------------*/



