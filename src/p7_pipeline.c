/* BATH's seq/profile comparison pipeline (modifies from HMMER3)
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
#include "hmmer.h"
#include "esl_gencode.h"

/* Struct used to pass a collection of useful temporary objects around*/
typedef struct {
  ESL_SQ           *tmpseq; // - a new or reused digital sequence object used for p7_alidisplay_Create() call
  P7_OMX          **oxf_holder; // - a tmeporary list of forward parser matricies for ORFs
} P7_PIPELINE_BATH_OBJS;

/* Struct used to keep track of the # and length of ORFs passing filters */
typedef struct {
  int64_t          *orf_starts;
  int64_t          *orf_ends;
  int               orf_cnt;
} P7_ORF_COORDS;

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
 *            | --fs         |  use frameshift aware algorithms            |   FALSE   |
 *            | --fs         |  use ONLY frameshift aware algorithms       |   FALSE   |
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
 *            with <p7_pipeline_Destroy_BATH()>.
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

  pli->do_alignment_score_calc = 0;

  /* Set Frameshift Mode */
  pli->frameshift = TRUE;
  pli->is_translated = FALSE; 
  pli->fs_pipe  = (go ? (esl_opt_IsUsed(go, "--fs") || esl_opt_IsUsed(go, "--fsonly")) : 0); 
  pli->std_pipe = (go ? !esl_opt_IsUsed(go, "--fsonly") : 1);

  /* Create forward and backward optimized matricies for use in the 
   * non-frameshift pipeline branch
   */
   if ((pli->oxf  = p7_omx_Create(M_hint, 0, L_hint)) == NULL) goto ERROR;
   if ((pli->oxb  = p7_omx_Create(M_hint, 0, L_hint)) == NULL) goto ERROR;

  /* Create sparce memeory forward and backward generic frameshift aware matricies 
   * for use in the frameshift pipeline filters
   */  
   if ((pli->gxf  = p7_gmx_fs_Create(M_hint, PARSER_ROWS_FWD, L_hint, 0)) == NULL) goto ERROR;
   if ((pli->gxb  = p7_gmx_fs_Create(M_hint, PARSER_ROWS_BWD, L_hint, 0)) == NULL) goto ERROR;   

  /* Create full memeory forward and backward generic frameshift aware matricies
   * for use in the frameshift pipeline alignment
   */ 
   if ((pli->gfwd = p7_gmx_fs_Create(M_hint, L_hint, L_hint, p7P_5CODONS)) == NULL) goto ERROR;
   if ((pli->gbck = p7_gmx_fs_Create(M_hint, L_hint, L_hint,  0))         == NULL) goto ERROR;

  /* Create intermediate values matrix */
   if ((pli->iv  = p7_ivx_Create(M_hint, p7P_3CODONS)) == NULL) goto ERROR;

  /* Normally, we reinitialize the RNG to the original seed every time we're
   * about to collect a stochastic trace ensemble. This eliminates run-to-run
   * variability. As a special case, if seed==0, we choose an arbitrary one-time 
   * seed: time() sets the seed, and we turn off the reinitialization.
   */
   pli->r                  =  esl_randomness_CreateFast(seed);
   pli->do_reseeding       = (seed == 0) ? FALSE : TRUE;
   pli->ddef               = p7_domaindef_fs_Create(pli->r, go);
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
   pli->B1     = (go ? esl_opt_GetInteger(go, "--B1") : 100);
   pli->B2     = (go ? esl_opt_GetInteger(go, "--B2") : 240);
   pli->B3     = (go ? esl_opt_GetInteger(go, "--B3") : 1000);

   if (go && esl_opt_GetBoolean(go, "--max")) 
   {
    pli->do_max        = TRUE;
    pli->do_biasfilter = FALSE;

    pli->F2 = pli->F3 = pli->F4 = .0;
    pli->F1 = 1.0; // need to set some threshold for F1 even on long targets. Should this be tighter?
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
 * Synopsis:  Reuse a BATH pipeline for next target.
 *
 */
int
p7_pipeline_Reuse_BATH(P7_PIPELINE *pli)
{
  p7_gmx_Reuse(pli->gxf);
  p7_gmx_Reuse(pli->gxb);
  p7_gmx_Reuse(pli->gfwd);
  p7_gmx_Reuse(pli->gbck);
  p7_omx_Reuse(pli->oxf);
  p7_omx_Reuse(pli->oxb);
  p7_domaindef_Reuse(pli->ddef);
  return eslOK;
}

/* Function:  p7_pipeline_Destroy_BATH() -
 * Synopsis:  Free a BATH <P7_PIPELINE> object.
 *
 */
void
p7_pipeline_Destroy_BATH(P7_PIPELINE *pli)
{
  if (pli == NULL) return;
 
  p7_gmx_Destroy(pli->gxf);
  p7_gmx_Destroy(pli->gxb);
  p7_gmx_Destroy(pli->gfwd);
  p7_gmx_Destroy(pli->gbck);
  p7_omx_Destroy(pli->oxf);
  p7_omx_Destroy(pli->oxb);
  p7_ivx_Destroy(pli->iv);
  esl_randomness_Destroy(pli->r);
  p7_domaindef_fs_Destroy(pli->ddef);
  free(pli);
}

/*---------------- end, P7_PIPELINE object ----------------------*/


/*****************************************************************
 * 2. The pipeline API.
 *****************************************************************/

/* Function:  p7_pli_ExtendAndMergeORFs -  
 * Synopsis:  Creates a DNA window around the coodinated of an ORF 
 *
 * Purpose:   Accepts a <orf_block> of ORFs, extends the DNA coordinates 
 *            based on a combination of the max_length value from <om> 
 *            and the prefix and suffix lengths stored in <data>, then 
 *            merges (in place) ORFs whose DNA coordinates overlap by 
 *            more than <pct_overlap> percent, ensuring that coordinates
 *            stay within the bounds of 1..<dna_sq->L>.
 *
 * Returns:   <eslOK>
 */
int
p7_pli_ExtendAndMergeORFs (P7_PIPELINE *pli, ESL_SQ_BLOCK *orf_block, ESL_SQ *dna_sq, P7_PROFILE *gm, P7_BG *bg, const P7_SCOREDATA *data, P7_HMM_WINDOWLIST *windowlist, float pct_overlap, int complementarity, int64_t seqidx, int32_t *k_coords_list, int32_t *m_coords_list) { 

  int            i;
  int            new_hit_cnt;

  int32_t        i_coords, j_coords;          /* original ORF hit coords */
  int32_t        ext_i_coords, ext_j_coords;  /* extended ORF hit coords */
  int32_t        k_coords, m_coords;          /* original HMM hit coords */
  
  int64_t        window_start, window_end;
  int64_t        overlap_len;
  int64_t        max_window_start, min_window_end;
  int64_t        min_window_start, max_window_end;
 
  float          vsc;
  float          filtersc;
  float          nullsc;
  float          seq_score;
  double         P; 
  
  ESL_SQ        *curr_orf      = NULL;
  P7_HMM_WINDOW *prev_window   = NULL;
  P7_HMM_WINDOW *curr_window   = NULL;
  
  P7_GMX        *vgx           = NULL;
  P7_TRACE      *vtr           = NULL;

  if (orf_block->count == 0)
    return eslOK;
  vgx = p7_gmx_Create(gm->M, 100);
  vtr = p7_trace_Create();

  /* extend each ORF's DNA coordinates based on the model's max length*/
  
  for(i = 0; i < orf_block->count; i++)
  {
    curr_orf = &(orf_block->list[i]);
    p7_gmx_GrowTo(vgx, gm->M, curr_orf->n); 
    p7_ReconfigLength(gm, curr_orf->n);
    
    p7_GViterbi(curr_orf->dsq, curr_orf->n, gm, vgx, &vsc);
    p7_GTrace(curr_orf->dsq, curr_orf->n, gm, vgx, vtr); 
    p7_trace_GetDomainCoords(vtr, 0, &i_coords, &j_coords, &k_coords, &m_coords);

	k_coords_list[i] = k_coords;
    m_coords_list[i] = m_coords;

	/* Rescore bias based on viterbi alignment */
    if (pli->do_biasfilter)
    {
       p7_bg_SetLength(bg, curr_orf->n);
       p7_bg_FilterScore(bg, curr_orf->dsq+i_coords-1, (j_coords-i_coords+1), &filtersc);

       /* Subtract out alignment length nullsc and add orf length null score back in */
       nullsc = (float) (j_coords-i_coords+1) * logf(bg->p1) + logf(1.-bg->p1);
       filtersc -= nullsc;
       nullsc = (float) curr_orf->n * logf(bg->p1) + logf(1.-bg->p1);
       filtersc += nullsc;

       seq_score = (vsc - filtersc) / eslCONST_LOG2;
       P = esl_gumbel_surv(seq_score, gm->evparam[p7_VMU], gm->evparam[p7_VLAMBDA]);
       if (P > pli->F2) {
         p7_gmx_Reuse(vgx);
         p7_trace_Reuse(vtr);
         continue;
       }
    }

    ext_i_coords = ESL_MIN(0,           (i_coords - (gm->max_length * (0.1 + data->prefix_lengths[k_coords]))-1)); //negeative numbers
    ext_j_coords = ESL_MAX(curr_orf->n, (j_coords + (gm->max_length * (0.1 + data->suffix_lengths[m_coords]))+1)); //positive numbers
     
    if(complementarity == p7_NOCOMPLEMENT)
    {
      window_start = ESL_MAX(1,         curr_orf->start + (ext_i_coords * 3));
      window_end   = ESL_MIN(dna_sq->n, curr_orf->start + (ext_j_coords * 3));
     
    }
    else
    {
      window_start   = ESL_MAX(1,         (dna_sq->n - curr_orf->start + 1) + (ext_i_coords * 3));
      window_end     = ESL_MIN(dna_sq->n, (dna_sq->n - curr_orf->start + 1) + (ext_j_coords * 3)); 
    }
    
    p7_hmmwindow_new(windowlist, 0, window_start, window_start-1, k_coords, window_end-window_start+1, 0.0, complementarity, dna_sq->n);
    p7_gmx_Reuse(vgx);
    p7_trace_Reuse(vtr);
  }


  if(windowlist->count == 0) {
    p7_gmx_Destroy(vgx);
    p7_trace_Destroy(vtr);
    return eslOK;
  }

  p7_hmmwindow_SortByStart(windowlist); 
  new_hit_cnt = 0;
   
   /* merge overlapping windows, compressing list in place. */
   
   for (i=1; i<windowlist->count; i++) {
    prev_window = windowlist->windows+new_hit_cnt;
    curr_window = windowlist->windows+i;
    
    max_window_start =  ESL_MAX(prev_window->n, curr_window->n);
    min_window_end   = ESL_MIN(prev_window->n+prev_window->length-1, curr_window->n+curr_window->length-1);
    overlap_len   = min_window_end - max_window_start + 1;
    
    if (  prev_window->complementarity == curr_window->complementarity &&
          (float)(overlap_len)/ESL_MIN(prev_window->length, curr_window->length) > pct_overlap )
    {

      min_window_start        = ESL_MIN(prev_window->n, curr_window->n);
      max_window_end          = ESL_MAX(prev_window->n+prev_window->length-1, curr_window->n+curr_window->length-1);
      /* If length of merged window would not be too long then merge windows */

      if((max_window_end -  min_window_start + 1) < (2 * (gm->max_length * 3))) 
      {
        prev_window->fm_n  -= (prev_window->n - min_window_start);
        prev_window->n      = min_window_start;
        prev_window->length = max_window_end - min_window_start + 1;
      } else {
        new_hit_cnt++;
        windowlist->windows[new_hit_cnt] = windowlist->windows[i];
      }
    } else {
      new_hit_cnt++;
      windowlist->windows[new_hit_cnt] = windowlist->windows[i];
    }
  }
  
  windowlist->count = new_hit_cnt+1;

  p7_gmx_Destroy(vgx);
  p7_trace_Destroy(vtr);
  return eslOK;
}

/* Function:  p7_pli_etPosPast_BATH
 * Synopsis:  Counts DNA positions passing MSV, Viterbi or Forward filter as ORFs 
 *
 * Purpose:    Uses ORF start and end positons on the DNA sequence to 
 *             calculate the number of unique DNA positions passing a 
 *             filter in the frameshift pipeline. Assumes that ORF 
 *             coordinates are in order of start position.             
 */
int
p7_pli_GetPosPast_BATH(P7_ORF_COORDS *coords)
{
  int            i;
  int            pos_cnt = 0;
  int            new_orf_cnt = 0;
  int64_t        prev_start, curr_start;
  int64_t        prev_end, curr_end;
  int64_t        min_start, max_start;
  int64_t        min_end, max_end;
  int64_t        overlap, max_length; 

  if(coords->orf_cnt == 0) return pos_cnt;

  for(i = 1; i < coords->orf_cnt; i++)
  {
    prev_start = coords->orf_starts[new_orf_cnt];
    curr_start = coords->orf_starts[i];
    prev_end =   coords->orf_ends[new_orf_cnt];
    curr_end =   coords->orf_ends[i];
     
    max_start        = ESL_MAX(prev_start, curr_start);
    min_end          = ESL_MIN(prev_end, curr_end);
    min_start        = ESL_MIN(prev_start, curr_start);
    max_end          = ESL_MAX(prev_end, curr_end);
       
    overlap          = min_end - max_start + 1;
    max_length       = max_end - min_start + 1;
   
    if ( (float) overlap / max_length > 0.)
    {
      coords->orf_starts[new_orf_cnt] = min_start;
      coords->orf_ends[new_orf_cnt]   = max_end;
    } else new_orf_cnt++;
  }

  new_orf_cnt++; 
  for(i = 0; i < new_orf_cnt; i++)
    pos_cnt += coords->orf_ends[i] - coords->orf_starts[i] + 1;
  
  return pos_cnt;
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
  if      (  pli->by_E )
    { 
      if      ( pli->frameshift ) { if (exp(lnP) <= pli->E) return TRUE; }  
      else if ( exp(lnP) * pli->Z <= pli->E) return TRUE;
    }
  else if (! pli->by_E   && score         >= pli->T) return TRUE;

  return FALSE;
}

/* Function:  p7_pli_TargetIncludable()
 * Synopsis:  Returns TRUE if target score meets inclusion threshold.
 */
int
p7_pli_TargetIncludable(P7_PIPELINE *pli, float score, double lnP)
{
  if      (  pli->inc_by_E )
    {
      if      ( pli->frameshift ) {
        if (exp(lnP) <= pli->incE) return TRUE;
      }
      else if ( exp(lnP) * pli->Z <= pli->incE) return TRUE;
    }

  else if (! pli->inc_by_E   && score         >= pli->incT) return TRUE;

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

/* Function:  p7_pipeline_Merge_BATH()
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
p7_pipeline_Merge_BATH(P7_PIPELINE *p1, P7_PIPELINE *p2)
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

  if (p1->Z_setby == p7_ZSETBY_NTARGETS)
    {
      if(p2->frameshift)
        p1->Z += p2->Z;
      else 
        p1->Z += (p1->mode == p7_SCAN_MODELS) ? p2->nmodels : p2->nseqs;
    }
  else
    {
      p1->Z = p2->Z;
    }

  return eslOK;
}

/* Function:  p7_pli_computeAliScores()
 * Synopsis:  Compute per-position scores for the alignment for a domain
 *
 * Purpose:   Compute per-position (Viterbi) scores for the alignment for a domain,
 *            for the purpose of optionally printing these scores out in association
 *            with each alignment. Such scores can, for example, be used to detangle
 *            overlapping alignments (from different models)
 *
 * Args:      dom             - domain with the alignment for which we wish to compute scores
 *            seq             - sequence in which domain resides
 *            data            - contains model's emission and transition values in unstriped form
 *            K               - alphabet size
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
static int
p7_pli_computeAliScores(P7_DOMAIN *dom, ESL_DSQ *seq, const P7_SCOREDATA *data, int K)
{
  int status;
  int i, j, k;
  float sc;

  //Compute score contribution of each position in the alignment to the overall Viterbi score
  ESL_ALLOC( dom->scores_per_pos, sizeof(float) * dom->ad->N );
  for (i=0; i<dom->ad->N; i++)  dom->scores_per_pos[i] = 0.0;
  i = dom->iali - 1;        //sequence position
  j = dom->ad->hmmfrom - 1; //model position
  k = 0;
  while ( k<dom->ad->N) {
    if (dom->ad->model[k] != '.' && dom->ad->aseq[k] != '-') { //match
      i++;  j++;
      // Including the MM cost is a hack. The cost of getting to/from this match
      // state does matter, but an IM or DM transition would improperly deflate
      // the score of this column, so just give MM. That amount is offset out of
      // the score shown for preceding indels
      dom->scores_per_pos[k] = data->fwd_scores[K * j + seq[i]]
                             +  (j==1 ? 0 : log(data->fwd_transitions[p7O_MM][j]) );
      k++;
    } else if (dom->ad->model[k] == '.' ) { // insert
      //spin through the insert, accumulating cost;  only assign to final column in gap
      dom->scores_per_pos[k] = -eslINFINITY;

      sc = log(data->fwd_transitions[p7O_MI][j]);
      i++; k++;
      while (k<dom->ad->N && dom->ad->model[k] == '.') { //extend insert
        dom->scores_per_pos[k] = -eslINFINITY;
        sc += log(data->fwd_transitions[p7O_II][j]);
        i++; k++;
      }
      sc += log(data->fwd_transitions[p7O_IM][j+1]) - log(data->fwd_transitions[p7O_MM][j+1]);
      dom->scores_per_pos[k-1] = sc;

    } else if (dom->ad->aseq[k] == '-' ) { // delete
      dom->scores_per_pos[k] = -eslINFINITY;
      sc = log(data->fwd_transitions[p7O_MD][j]);
      j++; k++;
      while (k<dom->ad->N && dom->ad->aseq[k] == '-')  { //extend delete
        dom->scores_per_pos[k] = -eslINFINITY;
        sc += log(data->fwd_transitions[p7O_DD][j]);
        j++; k++;
      }
      sc += log(data->fwd_transitions[p7O_DM][j+1]) - log(data->fwd_transitions[p7O_MM][j+1]);
      dom->scores_per_pos[k-1] = sc;
    }
  }

  return eslOK;

ERROR:
  return eslEMEM;

}

/* Function:  p7_pli_postDomainDef_Frameshift_BATH()  
 * Synopsis:  the part of the BATH search Pipeline downstream
 *            of Domain Definition
 *
 * Purpose:   This is called by p7_pli_postViterbi_Frameshift_BATH(), 
 *            and runs the post-Domain Definition part of the 
 *            frameshift aware branch of the BATH pipeline. It 
 *            consists of running various bookkeeping and sanity 
 *            checks on hits 
 *
 * Args:      pli             - the main pipeline object
 *            gm_fs           - fs-aware codon profile (query)
 *            bg              - background model
 *            hitlist         - pointer to hit storage bin
 *            seqidx          - the id # of the target sequence from which 
 *                              the ORFs were translated
 *            window_start    - the starting position of the DNA window
 *            dnasq           - the target dna sequence
 *            complementarity - boolean; is the passed window sourced from a complementary sequence block
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 */
static int 
p7_pli_postDomainDef_Frameshift_BATH(P7_PIPELINE *pli, P7_FS_PROFILE *gm_fs, P7_BG *bg, P7_TOPHITS *hitlist, 
                              int64_t seqidx, int window_start, ESL_SQ *dnasq, int complementarity 
)
{

  int              d;
  int              ali_len;
  int              env_len;
  int              status;
  
  float            bitscore;
  float            dom_bias; 
  float            dom_score;
  float            nullsc;
  double           dom_lnP;
  P7_DOMAIN       *dom        = NULL;      /* convenience variable, ptr to current domain  */
  P7_HIT          *hit        = NULL;      /* ptr to the current hit output data                           */

  for (d = 0; d < pli->ddef->ndom; d++)
  {      
    dom = pli->ddef->dcl + d;
  
    ali_len = dom->jali - dom->iali + 1;
    bitscore = dom->envsc;
    
    if (ali_len < 12)   
    {// anything less than this is a funny byproduct of the Forward score passing a very low threshold, but no reliable alignment existing that supports it
      p7_alidisplay_Destroy(dom->ad);
      p7_trace_fs_Destroy(dom->tr);
      continue;
    }
 	
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

    dom->ad->sqfrom = dom->iali;
    dom->ad->sqto   = dom->jali;

    /* Adjust score from env_len to max window length. Note that the loop and move 
     * costs are calculated based on amino lengths but paid per nucleotide*/
    bitscore = dom->envsc;
    bitscore -= 2 * log(2. / ((env_len/3.)+2));
    bitscore += 2 * log(2. / (gm_fs->max_length+2));
    bitscore -= ((env_len-ali_len)/3)                              * log((float) (env_len/3.) / (float) ((env_len/3.)+2));
    bitscore += ((ESL_MAX(env_len,gm_fs->max_length*3)-ali_len)/3) * log((float) gm_fs->max_length / (float) (gm_fs->max_length+2));

    /* Bias calculation and adjustments to Forward score */
    if (pli->do_null2)
      dom_bias = p7_FLogsum(0.0, log(bg->omega) + dom->domcorrection);
    else
      dom_bias = 0.0; 

    p7_bg_SetLength(bg, ESL_MAX(env_len/3,gm_fs->max_length));
    p7_bg_fs_NullOne  (bg, dnasq->dsq, ESL_MAX(env_len/3,gm_fs->max_length), &nullsc);
    dom_score  = (bitscore - (nullsc + dom_bias))  / eslCONST_LOG2;
     
    /* P-vaule calculation */	
    dom_lnP   = esl_exp_logsurv(dom_score, gm_fs->evparam[p7_FTAUFS5], gm_fs->evparam[p7_FLAMBDA]);
     
    /* Check if hit passes the e-value cutoff based on the current
     * residue count. This prevents hits from accumulating and using
     * excessive memmory. */
    pli->Z = (float)pli->nres / (float)gm_fs->max_length;
    if ( exp(dom_lnP) * pli->Z <= pli->E ) 
    { 
   
      p7_tophits_CreateNextHit(hitlist, &hit);

      hit->ndom        = 1;
      hit->best_domain = 0;

      hit->window_length = gm_fs->max_length;
      hit->seqidx = seqidx;
      hit->subseq_start = dnasq->start;

      ESL_ALLOC(hit->dcl, sizeof(P7_DOMAIN) );
      hit->dcl[0] = pli->ddef->dcl[d];

      hit->dcl[0].ad->L = 0;     
	
      hit->pre_score = bitscore  / eslCONST_LOG2;
      hit->pre_lnP   = esl_exp_logsurv (hit->pre_score,  gm_fs->evparam[p7_FTAUFS5], gm_fs->evparam[p7_FLAMBDA]);

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
        if ((status  = esl_strdup(gm_fs->name, -1, &(hit->name)))  != eslOK) esl_fatal("allocation failure");
        if ((status  = esl_strdup(gm_fs->acc,  -1, &(hit->acc)))   != eslOK) esl_fatal("allocation failure");
        if ((status  = esl_strdup(gm_fs->desc, -1, &(hit->desc)))  != eslOK) esl_fatal("allocation failure");
      }
    }
    else  //delete unused P7_ALIDSPLAY and P7_TRACE
    {
      p7_alidisplay_Destroy(dom->ad);
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
 * Purpose:   This is called by p7_pli_postViterbi_Frameshift_BATH(), 
 *            or p7_pli_postViterbi_BATH() and runs the post-Domain 
 *            Definition part of non-frameshift aware branch of the 
 *            BATH pipeline. It consists of running various bookkeeping 
 *            and sanity checks on hits 
 *
 * Args:      pli             - the main pipeline object
 *            om              - optimized protien profile (query)  
 *            bg              - background model
 *            hitlist         - pointer to hit storage bin
 *            seqidx          - the id # of the target sequence from which 
 *                              the ORFs were translated
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
p7_pli_postDomainDef_BATH(P7_PIPELINE *pli, P7_OPROFILE *om, P7_BG *bg, P7_TOPHITS *hitlist, int64_t seqidx, int window_start, ESL_SQ *orfsq, ESL_SQ *dnasq, int complementarity, float nullsc)
{

  int              d;
  int              ali_len, env_len;
  int              status;
  float            bitscore;
  float            dom_score;
  float            dom_bias;
  float            dom_lnP;
  P7_DOMAIN       *dom        = NULL;      /* convenience variable, ptr to current domain  */
  P7_HIT          *hit        = NULL;      /* ptr to the current hit output data           */
 
  for (d = 0; d < pli->ddef->ndom; d++)
  {      
    dom = pli->ddef->dcl + d;

    env_len = dom->jenv - dom->ienv + 1;	
    ali_len = (dom->jali - dom->iali + 1) / 3;   
    if (ali_len < 4) 
    {  // anything less than this is a funny byproduct of the Forward score passing a very low threshold, but no reliable alignment existing that supports it
      p7_alidisplay_Destroy(dom->ad);
      p7_trace_fs_Destroy(dom->tr);
      continue; 
    }

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

    dom->ad->sqfrom = dom->iali;
    dom->ad->sqto   = dom->jali;

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

     /* To prevent the accumultion of excessive low quailty hits when 
      * filters are turned off we need to begin weeding out those hits 
      * now. To do this we estimate Z based on crruent target residue 
      * count. This will allways be an understimation so we don't risk 
      * thowing away good hits.  The ture Z is calcualted at the end 
      * by p7_tophits_ComputeBathEvalues() */
     pli->Z = (float)pli->nres / (float)om->max_length;
     /* Check if hit passes the e-value cutoff based on the current
      * residue count. This prevents hits from accumulating and using
      * excessive memmory. */
     if ( exp(dom_lnP) * pli->Z <= pli->E ) 
     { 
       /* Add hits to hitlist and check if they are reprotable*/   
       p7_tophits_CreateNextHit(hitlist, &hit);
    
       hit->ndom        = 1;
       hit->best_domain = 0;
       hit->window_length = orfsq->n;
       hit->seqidx = seqidx;
       hit->subseq_start = orfsq->start;

       ESL_ALLOC(hit->dcl, sizeof(P7_DOMAIN) );
       hit->dcl[0] = pli->ddef->dcl[d];

       hit->dcl[0].ad->L = 0;

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
        p7_alidisplay_Destroy(dom->ad);
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

/* Function:  p7_pli_postViterbi_Frameshift_BATH()
 * Synopsis:  the part of the BATH search Pipeline downstream
 *            of the Viterbi filter for framshift search
 *
 * Purpose:   This is called by p7_Pipeline_Frameshift(), and runs the
 *            post-Viterbi part of the frameshift aware pipeline. It 
 *            consists of running Forward filters on ORFs and protien 
 *            profiles <om> and on corresponding DNA windows and fs-
 *            aware codon models <gm-fs>. Whichever target-query pair 
 *            produces the lowest p-value (that also passes the F3 
 *            threshhold) will be run through the backward filter and 
 *            passed to the remained of the apporptirate branch of the 
 *            pipeline.  
 *
 * Args:      pli             - the main pipeline object
 *            om              - optimized protien profile (query)
 *            gm              - non-optimized protien profile (query)
 *            gm_fs           - fs-aware codon profile (query)
 *            bg              - background model
 *            hitlist         - pointer to hit storage bin
 *            seqidx          - the id # of the target sequence from which 
 *                              the ORFs were translated
 *            dna_window      - a window obeject with the start and length 
 *                              of the window and all associated orfs
 *            dnasq           - the target dna sequence
 *            wrk             - workstate for translating codons
 *            gcode           - genetic code information for codon translation
 *            pli_tmp         - frameshift pipeline object for use in domain definition 
 *            complementarity - boolean; is the passed window sourced from a complementary sequence block
 * Returns:   <eslOK> on success. If a significant hit is obtained,
 *            its information is added to the growing <hitlist>.
 *
 *            <eslERANGE> on numerical overflow errors in the optimized 
 *            vector implementations; particularly in non-frameshift
 *            posterior decoding. I don't believe this is possible for
 *            multihit local models, but I'm set up to catch it
 *            anyway. We may emit a warning to the user, but cleanly
 *            skip the problematic sequence and continue.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 */
static int
p7_pli_postViterbi_Frameshift_BATH(P7_PIPELINE *pli, P7_OPROFILE *om, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, P7_BG *bg, P7_TOPHITS *hitlist,  
                              int64_t seqidx, P7_HMM_WINDOW *dna_window, ESL_SQ_BLOCK *orf_block, ESL_SQ *dnasq, ESL_GENCODE_WORKSTATE *wrk, ESL_GENCODE *gcode,
                             P7_PIPELINE_BATH_OBJS *pli_tmp, int complementarity, int32_t *k_coords_list, int32_t *m_coords_list
)
{

  int              f;
  int              status;
  int64_t          window_start, window_end;
  int64_t          orf_start, orf_end; 
  int32_t          prev_k, prev_m;
  ESL_DSQ         *subseq;                     /* current DNA window holder                    */ 
  ESL_SQ          *curr_orf;                   /* current ORF holder                           */
  float            fwdsc_fs, fwdsc_orf;        /* forward scores                               */
  float            nullsc_orf;                 /* ORF null score for forward filter            */
  float            filtersc_fs, filtersc_orf;  /* total filterscs for forward filters          */
  float            seqscore_fs, seqscore_orf;  /* the corrected per-seq bit score              */
  float 	   tot_orf_sc;                 /* summed score for all ORFs in current DNA window */
  double           P_fs;                       /* P-value of frameshift forward for window*/
  double           P_fs_nobias;                /* P-value of frameshift forward for window w/o bias adjustments*/
  double           tot_orf_P;                  /* P-value of summed forward score for all ORFs */
  double           min_P_orf;                  /* lowest p-value produced by an ORF */
  double          *P_orf;                      /* list of standard forward P-values for each ORf*/

  pli_tmp->oxf_holder = NULL;

  subseq = dnasq->dsq + dna_window->n - 1;
  
  /* Get true coords */
  window_start = complementarity ? dnasq->start - (dna_window->n + dna_window->length) : dnasq->start + dna_window->n - 1; 
  window_end   = complementarity ? dnasq->start - dna_window->n + 1 : window_start + dna_window->length - 1;

  /*set up seq object for domaindef function*/
  if ((status = esl_sq_SetName     (pli_tmp->tmpseq, dnasq->name))   != eslOK) goto ERROR;
  if ((status = esl_sq_SetSource   (pli_tmp->tmpseq, dnasq->source)) != eslOK) goto ERROR;
  if ((status = esl_sq_SetAccession(pli_tmp->tmpseq, dnasq->acc))    != eslOK) goto ERROR;
  if ((status = esl_sq_SetDesc     (pli_tmp->tmpseq, dnasq->desc))   != eslOK) goto ERROR;

  pli_tmp->tmpseq->L = dna_window->length;
  pli_tmp->tmpseq->n = dna_window->length;
  pli_tmp->tmpseq->start = dna_window->n;
  pli_tmp->tmpseq->end = dna_window->n + dna_window->length - 1; 
  pli_tmp->tmpseq->dsq = subseq;
  
  tot_orf_sc = eslINFINITY;
  tot_orf_P  = eslINFINITY;
  min_P_orf  = eslINFINITY;
  P_orf = NULL;
  /* If this search is using the standard translation pipeline
   *  (user did not specify --fsonly) than run the standard
   *  Foward on every ORF that is within the current window */ 
  if(pli->std_pipe) {

    tot_orf_sc = -eslINFINITY;

    ESL_ALLOC(P_orf, sizeof(double) * orf_block->count);
    ESL_ALLOC(pli_tmp->oxf_holder, sizeof(P7_OMX *) * orf_block->count);

    prev_k = 0;
    prev_m = 0;

   for(f = 0; f < orf_block->count; f++) {
     curr_orf = &(orf_block->list[f]);
     pli_tmp->oxf_holder[f] = NULL;

     /* Convert to true coords */
     if(complementarity) {
       orf_start =  dnasq->start - (dnasq->n - curr_orf->end   + 1) + 1;
       orf_end   =  dnasq->start - (dnasq->n - curr_orf->start + 1) + 1;
     } else {    
       orf_start = dnasq->start + curr_orf->start - 1;
       orf_end   = dnasq->start + curr_orf->end   - 1;
     } 

     /* Only process ORF if it belongs to the current window */ 
     if(orf_start >= window_start && orf_end <= window_end) { 
	   
       p7_bg_SetLength(bg, curr_orf->n);
       p7_bg_NullOne  (bg, curr_orf->dsq, curr_orf->n, &nullsc_orf);
         
       if (pli->do_biasfilter)
         p7_bg_FilterScore(bg, curr_orf->dsq, curr_orf->n, &filtersc_orf);
       else filtersc_orf = nullsc_orf;

       /* Save the Forward Martix for each ORF so we do not have to rerun 
        * Forward in the event that the standard pipeline is selected */
       p7_oprofile_ReconfigLength(om, curr_orf->n);
       if ((pli_tmp->oxf_holder[f] = p7_omx_Create(om->M, 0, curr_orf->n)) == NULL) goto ERROR ;
       p7_ForwardParser(curr_orf->dsq, curr_orf->n, om, pli_tmp->oxf_holder[f], &fwdsc_orf);
       
       /* Find the individual p-value (with bias) of each ORF in 
        * the window and store it. Also find the minimum p-value 
        * (with bias) of all ORFs in the window to test it at 
        * least one ORF passed the Forward filter */ 
       seqscore_orf = (fwdsc_orf-filtersc_orf) / eslCONST_LOG2;
       P_orf[f] = esl_exp_surv(seqscore_orf,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]); 
       min_P_orf = ESL_MIN(min_P_orf, P_orf[f]);  
        
       /* Sum the scores of all ORFs in the window (without bias). 
        * If there are multiple ORFs, and thier Viterbi trace 
        * coordinates are within 20 Amino Acids, assume these are 
        * exons and remove addtional core model entry cost to 
        * simulate the use of an intron state rather than a J state*/  
         
       if(prev_m > 0 && prev_k < m_coords_list[f] && abs(k_coords_list[f] - prev_m) < 15) {   
         tot_orf_sc =  p7_FLogsum(tot_orf_sc, (fwdsc_orf - gm->tsc[p7P_NTRANS+p7P_BM]));
       } 
       else  { tot_orf_sc =  p7_FLogsum(tot_orf_sc, fwdsc_orf); }
    
        prev_k = k_coords_list[f];
        prev_m = m_coords_list[f];
      }
    }
    tot_orf_P = esl_exp_surv(tot_orf_sc / eslCONST_LOG2,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
  } 


  P_fs        = eslINFINITY;
  P_fs_nobias = eslINFINITY;

  /*If this search is using the frameshift aware pipeline 
   * (user did not specify --nofs) than run Frameshift 
   * Forward on full Window and save score and P value.*/

  if(pli->fs_pipe && (!pli->std_pipe || min_P_orf <= pli->F4)) {
    p7_bg_SetLength(bg, dna_window->length/3);
    p7_bg_fs_FilterScore(bg, pli_tmp->tmpseq, wrk, gcode, pli->do_biasfilter, &filtersc_fs);

    p7_gmx_fs_GrowTo(pli->gxf, gm_fs->M, PARSER_ROWS_FWD, dna_window->length, 0);
	p7_ivx_GrowTo(pli->iv, gm_fs->M, p7P_3CODONS);
    p7_fs_ReconfigLength(gm_fs, dna_window->length/3);
    
    p7_ForwardParser_Frameshift_3Codons(subseq, gcode, dna_window->length, gm_fs, pli->gxf, pli->iv, &fwdsc_fs);
    seqscore_fs = (fwdsc_fs-filtersc_fs) / eslCONST_LOG2;
    P_fs = esl_exp_surv(seqscore_fs,  gm_fs->evparam[p7_FTAUFS3],  gm_fs->evparam[p7_FLAMBDA]);
    P_fs_nobias = esl_exp_surv(fwdsc_fs/eslCONST_LOG2,  gm_fs->evparam[p7_FTAUFS3],  gm_fs->evparam[p7_FLAMBDA]); 
  }

    /* Compare Pvalues to select either the standard or the frameshift pipeline
   * If the DNA window passed frameshift forward AND produced a lower P-value 
   * than the sumed Forward score of the orfs used to costruct that window 
   * then we proceed with the frameshift pipeline
   */
  
  if(P_fs <= pli->F3 && (P_fs_nobias < tot_orf_P || min_P_orf > pli->F3)) { 
    
    pli->pos_past_fwd += dna_window->length; 
    p7_gmx_fs_GrowTo(pli->gxb, gm_fs->M, PARSER_ROWS_BWD, dna_window->length, 0);
    p7_BackwardParser_Frameshift_3Codons(subseq, gcode, dna_window->length, gm_fs, pli->gxb, pli->iv, NULL);
 
    status = p7_domaindef_ByPosteriorHeuristics_Frameshift(pli_tmp->tmpseq, gm, gm_fs,
           pli->gxf, pli->gxb, pli->gfwd, pli->gbck, pli->iv, pli->ddef, bg, gcode,
           dna_window->n, pli->do_biasfilter);
    if (status != eslOK) ESL_FAIL(status, pli->errbuf, "domain definition workflow failure"); 
    if (pli->ddef->nregions == 0)  return eslOK; /* score passed threshold but there's no discrete domains here     */
    if (pli->ddef->nenvelopes ==   0)  return eslOK; /* rarer: region was found, stochastic clustered, no envelope found*/
   
    /* Send any hits from the Frameshift aware pipeline to be further processed */ 
    p7_pli_postDomainDef_Frameshift_BATH(pli, gm_fs, bg, hitlist, seqidx, dna_window->n, dnasq, complementarity);

  } 

  /* If the DNA window did NOT pass pass frameshift forward OR did NOT produced a 
   * lower P-value than the sumed Forward score of the ORFs used to costruct that 
   * window then we check each of those ORFs individually to to see if they pass 
   * the standard Forward filter and it they do proceed with the standard pipeline
   */
   else if (pli->std_pipe) { 
  
     for(f = 0; f < orf_block->count; f++) {	
      curr_orf = &(orf_block->list[f]);
      if(complementarity) {
         orf_start =  dnasq->start - (dnasq->n - curr_orf->end   + 1) + 1;
         orf_end   =  dnasq->start - (dnasq->n - curr_orf->start + 1) + 1;
       } else {    
         orf_start = dnasq->start + curr_orf->start - 1;
         orf_end   = dnasq->start + curr_orf->end   - 1;
       } 

      /* Ensure current ORF is within the current window and that it passed  the Forward filter */
      if(orf_start >= window_start && orf_end <= window_end && P_orf[f] <= pli->F3) { 
        pli->pos_past_fwd += curr_orf->n * 3;
        p7_oprofile_ReconfigLength(om, curr_orf->n);
        p7_omx_GrowTo(pli->oxb, om->M, 0, curr_orf->n);     
        
        p7_BackwardParser(curr_orf->dsq, curr_orf->n, om, pli_tmp->oxf_holder[f], pli->oxb, NULL);
        
        status = p7_domaindef_ByPosteriorHeuristics_nonFrameshift(curr_orf, pli_tmp->tmpseq, dnasq->n, gcode, om, gm, gm_fs, pli_tmp->oxf_holder[f], pli->oxf, pli->oxb, pli->ddef, bg);
        if (status != eslOK) ESL_FAIL(status, pli->errbuf, "domain definition workflow failure"); /* eslERANGE can happen */
        if (pli->ddef->nregions   == 0)  continue; /* score passed threshold but there's no discrete domains here     */
        if (pli->ddef->nenvelopes == 0)  continue; /* rarer: region was found, stochastic clustered, no envelope found*/
        
        /* Send any hits from the standard pipeline to be further processed */   
        p7_pli_postDomainDef_BATH(pli, om, bg, hitlist, seqidx, dna_window->n, curr_orf, dnasq, complementarity, nullsc_orf);
      }
    }  
  } 

  /* clean up */  
  if(P_orf != NULL)  
    free(P_orf);
   
  if(pli_tmp->oxf_holder != NULL) 
  {
    for(f = 0; f < orf_block->count; f++) 
      p7_omx_Destroy(pli_tmp->oxf_holder[f]);
    free(pli_tmp->oxf_holder);
  } 
    
  return eslOK;

ERROR:
  ESL_EXCEPTION(eslEMEM, "Error in Frameshift pipeline\n");

  if(P_orf != NULL)
    free(P_orf);

  if(pli_tmp->oxf_holder != NULL) 
  {
    for(f = 0; f < orf_block->count; f++) 
      p7_omx_Destroy(pli_tmp->oxf_holder[f]);
    free(pli_tmp->oxf_holder);
  } 

}



/* Function:  p7_pli_postViterbi_BATH()
 * Synopsis:  the part of the BATH search Pipeline downstream
 *            of the Viterbi filter for non-frameshift search
 *
 * Purpose:   This is called by p7_Pipeline_Frameshift(), and runs the
 *            post-Viterbi part of the frameshift aware pipeline. It 
 *            consists of running Forward filters on ORFs and protien 
 *            profiles <om> and on corresponding DNA windows and fs-
 *            aware codon models <gm-fs>. Whichever target-query pair 
 *            produces the lowest p-value (that also passes the F3 
 *            threshhold) will be run through the backward filter and 
 *            passed to the remained of the apporptirate branch of the 
 *            pipeline.  
 *
 * Args:      pli             - the main pipeline object
 *            om              - optimized protien profile (query)
 *            gm              - non-optimized protien profile (query)
 *            gm_fs           - fs-aware codon profile (query)
 *            bg              - background model
 *            hitlist         - pointer to hit storage bin
 *            seqidx          - the id # of the target sequence from which 
 *                              the ORFs were translated
 *            dna_window      - a window obeject with the start and length 
 *                              of the window and all associated orfs
 *            dnasq           - the target dna sequence
 *            wrk             - workstate for translating codons
 *            gcode           - genetic code information for codon translation
 *            pli_tmp         - frameshift pipeline object for use in domain definition 
 *            complementarity - boolean; is the passed window sourced from a complementary sequence block
 * Returns:   <eslOK> on success. If a significant hit is obtained,
 *            its information is added to the growing <hitlist>.
 *
 *            <eslERANGE> on numerical overflow errors in the optimized 
 *            vector implementations; particularly in non-frameshift
 *            posterior decoding. I don't believe this is possible for
 *            multihit local models, but I'm set up to catch it
 *            anyway. We may emit a warning to the user, but cleanly
 *            skip the problematic sequence and continue.
 *
 *
 */
static int
p7_pli_postViterbi_BATH(P7_PIPELINE *pli, P7_OPROFILE *om, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, P7_BG *bg, P7_TOPHITS *hitlist, int64_t seqidx, ESL_SQ_BLOCK *orf_block, ESL_SQ *dnasq, ESL_GENCODE_WORKSTATE *wrk, ESL_GENCODE *gcode, P7_PIPELINE_BATH_OBJS *pli_tmp, int complementarity
)
{

  int              f;
  int              status;
  int64_t          orf_start, orf_end;
  float            fwdsc;                      /* forward scores                               */
  float            nullsc;                     /* ORF null score for forward filter            */
  float            filtersc;                   /* total filterscs for forward filters          */
  float            seqscore;                   /* the corrected per-seq bit score              */
  double           P;                          /* lowest p-value produced by an ORF */
  ESL_SQ          *curr_orf;                   /* current ORF holder                           */


  /*set up seq object for domaindef function*/
  if ((status = esl_sq_SetName     (pli_tmp->tmpseq, dnasq->name))   != eslOK) goto ERROR;
  if ((status = esl_sq_SetSource   (pli_tmp->tmpseq, dnasq->source)) != eslOK) goto ERROR;
  if ((status = esl_sq_SetAccession(pli_tmp->tmpseq, dnasq->acc))    != eslOK) goto ERROR;
  if ((status = esl_sq_SetDesc     (pli_tmp->tmpseq, dnasq->desc))   != eslOK) goto ERROR;

    /* Loop through ORFs */
  for(f = 0; f < orf_block->count; f++) {
    curr_orf = &(orf_block->list[f]);

	if(complementarity) {
	  orf_start = dnasq->n - curr_orf->start + 1;
	  orf_end   = dnasq->n - curr_orf->end + 1;
	}
	else {
	  orf_start = curr_orf->start;
	  orf_end   = curr_orf->end;
	}

    pli_tmp->tmpseq->start = orf_start; 
    pli_tmp->tmpseq->end   = orf_end;
    pli_tmp->tmpseq->dsq   = dnasq->dsq + orf_start - 1; 
    pli_tmp->tmpseq->n     = orf_end - orf_start + 1;
    pli_tmp->tmpseq->L     = pli_tmp->tmpseq->n; 

    p7_bg_SetLength(bg, curr_orf->n);
    p7_bg_NullOne  (bg, curr_orf->dsq, curr_orf->n, &nullsc);

    if (pli->do_biasfilter)
      p7_bg_FilterScore(bg, curr_orf->dsq, curr_orf->n, &filtersc);
    else filtersc = nullsc;
 
    p7_oprofile_ReconfigLength(om, curr_orf->n);
    p7_omx_GrowTo(pli->oxf, om->M, 0, curr_orf->n);
    p7_ForwardParser(curr_orf->dsq, curr_orf->n, om, pli->oxf, &fwdsc);
    seqscore = (fwdsc-filtersc) / eslCONST_LOG2;
    P = esl_exp_surv(seqscore,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);

    if (P > pli->F3) continue; 

    pli->pos_past_fwd += curr_orf->n * 3;
    p7_omx_GrowTo(pli->oxb, om->M, 0, curr_orf->n);

    p7_BackwardParser(curr_orf->dsq, curr_orf->n, om, pli->oxf, pli->oxb, NULL);

    status = p7_domaindef_ByPosteriorHeuristics_nonFrameshift(curr_orf, pli_tmp->tmpseq, dnasq->n, gcode, om, gm, gm_fs, pli->oxf, pli->oxf, pli->oxb, pli->ddef, bg);
    if (status != eslOK) ESL_FAIL(status, pli->errbuf, "domain definition workflow failure"); /* eslERANGE can happen */
    if (pli->ddef->nregions   == 0)  continue; /* score passed threshold but there's no discrete domains here     */
    if (pli->ddef->nenvelopes == 0)  continue; /* rarer: region was found, stochastic clustered, no envelope found*/

    p7_pli_postDomainDef_BATH(pli, om, bg, hitlist, seqidx, orf_start, curr_orf, dnasq, complementarity, nullsc);
  }
    

  return eslOK;

ERROR:
  return status;

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
 *            aware codon model <gm_fs> and a comparison is preformed. 
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
 *            gm_fs           - generic fs-aware codon profile (query)
 *            data            - for picking window edges based on 
 *                              maximum prefix/suffix extensions
 *            bg              - background model
 *            hitlist         - pointer to hit storage bin (already 
 *                              allocated)
 *            seqidx          - the id # of the sequence from which 
 *                              the current window was extracted
 *            dnasq           - digital sequence of the DNA window
 *            orf_block       - collection of ORFs translated form <dnasq>
 *            wrk             - codon translation workstate
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
p7_Pipeline_BATH(P7_PIPELINE *pli, P7_OPROFILE *om, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, P7_SCOREDATA *data, P7_BG *bg, P7_TOPHITS *hitlist, int64_t seqidx, ESL_SQ *dnasq, ESL_SQ_BLOCK *orf_block, ESL_GENCODE_WORKSTATE *wrk, ESL_GENCODE *gcode, int complementarity)
{

  int                i;
  int                status, wstatus;
  float              nullsc;              /* null model score                        */
  float              usc;                 /* msv score                               */
  float              vfsc;                /* viterbi score                           */
  float              seq_score;           /* null corrected bit score                */
  float              filtersc;            /* bias and null score                     */
  double             P;                   /* p-value holder                          */
  int                window_len;          /* length of DNA window                    */
  int                min_length;          /* minimum number of nucs passing a filter */
  int32_t           *k_coords_list, *m_coords_list; /* ORF Viterbi trace HMM coords            */
  ESL_SQ            *orfsq;               /* ORF sequence                            */
  ESL_SQ_BLOCK      *post_vit_orf_block;  /* block of ORFs that pass viterbi         */
  P7_HMM_WINDOWLIST  post_vit_windowlist; /* list of windows from ORFs that pass viterbi */
  P7_ORF_COORDS     *msv_coords, *bias_coords, *vit_coords;  /* number of nucleotieds passing filters */
  P7_PIPELINE_BATH_OBJS *pli_tmp;   

  if (dnasq->n < 15) return eslOK;         //DNA to short
  if (orf_block->count == 0) return eslOK; //No ORFS translated

  post_vit_orf_block = NULL;
  post_vit_orf_block = esl_sq_CreateDigitalBlock(orf_block->listSize, om->abc);
  post_vit_windowlist.windows = NULL;

  k_coords_list = NULL; 
  m_coords_list = NULL;
 
  msv_coords = NULL;
  bias_coords = NULL;
  vit_coords = NULL;

  pli_tmp = NULL;
  ESL_ALLOC(pli_tmp, sizeof(P7_PIPELINE_BATH_OBJS));
  pli_tmp->tmpseq = NULL;

  ESL_ALLOC(msv_coords, sizeof(P7_ORF_COORDS));
  ESL_ALLOC(msv_coords->orf_starts, sizeof(int64_t) *  orf_block->count);
  ESL_ALLOC(msv_coords->orf_ends, sizeof(int64_t) *  orf_block->count);
  msv_coords->orf_cnt = 0;

  ESL_ALLOC(bias_coords, sizeof(P7_ORF_COORDS));
  ESL_ALLOC(bias_coords->orf_starts, sizeof(int64_t) *  orf_block->count);
  ESL_ALLOC(bias_coords->orf_ends, sizeof(int64_t) *  orf_block->count);
  bias_coords->orf_cnt = 0;

  ESL_ALLOC(vit_coords, sizeof(P7_ORF_COORDS));
  ESL_ALLOC(vit_coords->orf_starts, sizeof(int64_t) *  orf_block->count);
  ESL_ALLOC(vit_coords->orf_ends, sizeof(int64_t) *  orf_block->count);
  vit_coords->orf_cnt = 0;
  
  for (i = 0; i < orf_block->count; ++i)
  { 
    orfsq = &(orf_block->list[i]);
   
    if (   (orfsq->start < orfsq->end    &&  orfsq->end < dnasq->C )  ||
           (orfsq->end < orfsq->start    &&  orfsq->start < dnasq->C ) )
        continue; /* don't bother with an orf that showed up completely within a previous window */

    /* use the name, accession, and description from the DNA sequence and
     not from the ORF which is generated by gencode and only for internal use */
    if ((wstatus = esl_sq_SetAccession(orfsq, dnasq->acc))    != eslOK)  ESL_EXCEPTION_SYS(eslEWRITE, "Set query sequence accession failed");
    if ((wstatus = esl_sq_SetDesc     (orfsq, dnasq->desc))   != eslOK)  ESL_EXCEPTION_SYS(eslEWRITE, "Set query sequence description failed");
    
    if(orfsq->n > 0) 
    {
      p7_bg_SetLength(bg, orfsq->n);
      p7_oprofile_ReconfigLength(om, orfsq->n);
      p7_bg_NullOne  (bg, orfsq->dsq, orfsq->n, &nullsc);
	
      p7_omx_GrowTo(pli->oxf, om->M, 0, orfsq->n);    /* expand the one-row omx if needed */
	
      /* MSV Filter on ORF */
      p7_MSVFilter(orfsq->dsq, orfsq->n, om, pli->oxf, &usc);
      seq_score = (usc - nullsc) / eslCONST_LOG2;
      P = esl_gumbel_surv( seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
      if (P > pli->F1 ) continue;    

      msv_coords->orf_starts[msv_coords->orf_cnt] = ESL_MIN(orfsq->start, orfsq->end);
      msv_coords->orf_ends[msv_coords->orf_cnt] =   ESL_MAX(orfsq->start, orfsq->end);
      msv_coords->orf_cnt++;
       
      /* biased composition HMM filtering */
      if (pli->do_biasfilter)
      {
        p7_bg_FilterScore(bg, orfsq->dsq, orfsq->n, &filtersc);

        seq_score = (usc - filtersc) / eslCONST_LOG2;
        P = esl_gumbel_surv(seq_score,  om->evparam[p7_MMU],  om->evparam[p7_MLAMBDA]);
        if (P > pli->F1) continue;
      }  else filtersc = nullsc;

      bias_coords->orf_starts[bias_coords->orf_cnt] = ESL_MIN(orfsq->start, orfsq->end);
      bias_coords->orf_ends[bias_coords->orf_cnt] =   ESL_MAX(orfsq->start, orfsq->end);
      bias_coords->orf_cnt++;

      /* Viterbi filer on ORF */
      if (P > pli->F2)
      {
        p7_ViterbiFilter(orfsq->dsq, orfsq->n, om, pli->oxf, &vfsc);
        seq_score = (vfsc-filtersc) / eslCONST_LOG2;
        P  = esl_gumbel_surv(seq_score,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);
        if (P > pli->F2) continue;
      }

      vit_coords->orf_starts[vit_coords->orf_cnt] = ESL_MIN(orfsq->start, orfsq->end);
      vit_coords->orf_ends[vit_coords->orf_cnt] =   ESL_MAX(orfsq->start, orfsq->end);
      vit_coords->orf_cnt++;
      
      /* Collect all ORFs which passed Viterbi Filter */
      esl_sq_Copy(orfsq, &(post_vit_orf_block->list[post_vit_orf_block->count]));
      post_vit_orf_block->count++;
    }
  }

  min_length = ESL_MIN(dnasq->n, om->max_length * 3);
  pli->pos_past_msv += ESL_MAX(p7_pli_GetPosPast_BATH(msv_coords), min_length);
  pli->pos_past_bias += ESL_MAX(p7_pli_GetPosPast_BATH(bias_coords), min_length);
  pli->pos_past_vit += ESL_MAX(p7_pli_GetPosPast_BATH(vit_coords), min_length);

  pli_tmp->tmpseq = esl_sq_CreateDigital(dnasq->abc);
  free (pli_tmp->tmpseq->dsq); 
 
  /* For frameshift search - create DNA widnows for ORFs that pass viterbi */
  if(pli->fs_pipe) {
    if (data->prefix_lengths == NULL)  //otherwise, already filled in
      p7_hmm_ScoreDataComputeRest(om, data);
  
    /* convert block of ORFs that passed Viterbi into collection of non-overlapping DNA windows */
    p7_hmmwindow_init(&post_vit_windowlist);  
    
    if(post_vit_orf_block->count > 0) 
    {
      ESL_ALLOC(k_coords_list, sizeof(int32_t) * post_vit_orf_block->count); 
      ESL_ALLOC(m_coords_list, sizeof(int32_t) * post_vit_orf_block->count);
    }
  
    p7_pli_ExtendAndMergeORFs (pli, post_vit_orf_block, dnasq, gm, bg, data, &post_vit_windowlist, 0., complementarity, seqidx, k_coords_list, m_coords_list); 
  
    /* Send ORFs and protien models along with DNA windows and fs-aware coddon models to Forward filters */
    for(i = 0; i < post_vit_windowlist.count; i++)
    {
      window_len   = post_vit_windowlist.windows[i].length; 
      if (window_len < 15) continue;
      p7_pli_postViterbi_Frameshift_BATH(pli, om, gm, gm_fs, bg, hitlist, seqidx, &(post_vit_windowlist.windows[i]), post_vit_orf_block, dnasq, wrk, gcode, pli_tmp, complementarity, k_coords_list, m_coords_list);
    }
  }
  else {
    
    p7_pli_postViterbi_BATH(pli, om, gm, gm_fs, bg, hitlist, seqidx, post_vit_orf_block, dnasq, wrk, gcode, pli_tmp, complementarity);
  }


  /* clean up */ 
  if ( msv_coords != NULL)
  {
    if (msv_coords->orf_starts != NULL) free(msv_coords->orf_starts);
    if (msv_coords->orf_ends != NULL) free(msv_coords->orf_ends);
    free(msv_coords);
  }
  if ( bias_coords != NULL)
  {
    if (bias_coords->orf_starts != NULL) free(bias_coords->orf_starts);
    if (bias_coords->orf_ends != NULL) free(bias_coords->orf_ends);
    free(bias_coords);
  }
  if ( vit_coords != NULL)
  {
    if (vit_coords->orf_starts != NULL) free(vit_coords->orf_starts);
    if (vit_coords->orf_ends != NULL) free(vit_coords->orf_ends);
    free(vit_coords);
  }
  if ( post_vit_orf_block != NULL) esl_sq_DestroyBlock(post_vit_orf_block); 
  pli_tmp->tmpseq->dsq = NULL;
  if (pli_tmp != NULL) 
  {
    if (pli_tmp->tmpseq != NULL)  esl_sq_Destroy(pli_tmp->tmpseq);
    free(pli_tmp);
  }
  if (post_vit_windowlist.windows != NULL) free (post_vit_windowlist.windows); 
  if(k_coords_list != NULL) free(k_coords_list);
  if(m_coords_list != NULL) free(m_coords_list);

  return eslOK;

ERROR:
  if ( msv_coords != NULL)
  {
    if (msv_coords->orf_starts != NULL) free(msv_coords->orf_starts);
    if (msv_coords->orf_ends != NULL) free(msv_coords->orf_ends);
    free(msv_coords);
  }

  if ( bias_coords != NULL)
  {
    if (bias_coords->orf_starts != NULL) free(bias_coords->orf_starts);
    if (bias_coords->orf_ends != NULL) free(bias_coords->orf_ends);
    free(bias_coords);
  }

  if ( vit_coords != NULL)
  {
    if (vit_coords->orf_starts != NULL) free(vit_coords->orf_starts);
    if (vit_coords->orf_ends != NULL) free(vit_coords->orf_ends);
    free(vit_coords);
  }

  if (pli_tmp != NULL)
  {
    if (pli_tmp->tmpseq != NULL)  esl_sq_Destroy(pli_tmp->tmpseq);
    free(pli_tmp);
  }

  if ( post_vit_orf_block != NULL) esl_sq_DestroyBlock(post_vit_orf_block); 
  if (post_vit_windowlist.windows != NULL) free (post_vit_windowlist.windows);
  if(k_coords_list != NULL) free(k_coords_list);
  if(m_coords_list != NULL) free(m_coords_list);
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
  double ntargets; 

  fprintf(ofp, "Internal pipeline statistics summary:\n");
  fprintf(ofp, "-------------------------------------\n");
  if (pli->mode == p7_SEARCH_SEQS) {
    fprintf(ofp, "Query model(s):              %15" PRId64 "  (%" PRId64 " nodes)\n",     pli->nmodels, pli->nnodes);
    fprintf(ofp, "Target %-12s          %15" PRId64 "  (%" PRId64 " residues searched)\n", pli->is_translated?"orf(s):":"sequence(s):", pli->nseqs,   pli->nres);
    ntargets = pli->nseqs;
  } else {
    fprintf(ofp, "Query %-12s           %15" PRId64 "  (%" PRId64 " residues searched)\n", pli->is_translated?"orf(s):":"sequence(s):", pli->nseqs,   pli->nres);
    fprintf(ofp, "Target model(s):             %15" PRId64 "  (%" PRId64 " nodes)\n",     pli->nmodels, pli->nnodes);
    ntargets = pli->nmodels;
  }

  if (pli->frameshift) { 
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

  } else { // typical case output

    fprintf(ofp, "Passed MSV filter:           %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n",
        pli->n_past_msv,
        (double) pli->n_past_msv / ntargets,
        pli->F1 * ntargets,
        pli->F1);

    fprintf(ofp, "Passed bias filter:          %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n",
        pli->n_past_bias,
        (double) pli->n_past_bias / ntargets,
        pli->F1 * ntargets,
        pli->F1);

    fprintf(ofp, "Passed Vit filter:           %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n",
        pli->n_past_vit,
        (double) pli->n_past_vit / ntargets,
        pli->F2 * ntargets,
        pli->F2);

    fprintf(ofp, "Passed Fwd filter:           %15" PRId64 "  (%.6g); expected %.1f (%.6g)\n",
        pli->n_past_fwd,
        (double) pli->n_past_fwd / ntargets,
        pli->F3 * ntargets,
        pli->F3);

    fprintf(ofp, "Initial search space (Z):    %15.0f  %s\n", pli->Z,    pli->Z_setby    == p7_ZSETBY_OPTION ? "[as set by --Z on cmdline]"    : "[actual number of targets]");
    fprintf(ofp, "Domain search space  (domZ): %15.0f  %s\n", pli->domZ, pli->domZ_setby == p7_ZSETBY_OPTION ? "[as set by --domZ on cmdline]" : "[number of targets reported over threshold]");
  }

  if (w != NULL) {
    esl_stopwatch_Display(ofp, w, "# CPU time: ");
    fprintf(ofp, "# Mc/sec: %.2f\n", 
        (double) pli->nres * (double) pli->nnodes / (w->elapsed * 1.0e6));
  }

  return eslOK;
}
/*------------------- end, pipeline API -------------------------*/



