/* SPLICE_PIPELINE: spliced seq/profile comparison pipeline
 *
 * Contents:
 *    1. The SPLICE_PIPELINE object.
 *    2. The SPLICE_SITE_IDX object.
 *
 */


#include "p7_config.h"

#include <string.h>

#include "easel.h"
#include "esl_vectorops.h"
#include "esl_gumbel.h"
#include "esl_exponential.h"

#include "hmmer.h"
#include "p7_splice.h"

/*****************************************************************
 * 1. The SPLICE_PIPELINE structure.
 *****************************************************************/

/* Splice singal probabilities taken from
 * "Comprehensive splice-site analysis using comparative genomics",
 * Nihar Sheth et al., 2006
 */
void
p7_splicepipeline_SignalScores(float *f)
{
  f[0] = log(0.9919);     /* GT-AG */
  f[1] = log(0.0073);     /* GC-AG */
  f[2] = log(0.0006);     /* AT-AC */
  return;
}

/* Function:  p7_splicepipeline_Create()
 *
 * Purpose:   Allocates a splice pipeline and set all relevant 
 *            data from the configuration structure <go>. 
 *            Allocate <M_hint> by <L_hint> alignment matricies.
 *
 * Returns:   a pointer to the new <SPLICE_PIPELINE> structure
 *            on success.
 *
 * Throws:    <NULL> on allocation error.
 */
SPLICE_PIPELINE* 
p7_splicepipeline_Create(const ESL_GETOPTS *go, int M_hint, int L_hint)
{
  SPLICE_PIPELINE *pli;
  int              status;
 
  pli = NULL;
  ESL_ALLOC(pli, sizeof(SPLICE_PIPELINE));

  pli->frameshift = TRUE;
  pli->long_targets = FALSE;

  if (go && esl_opt_GetBoolean(go, "--nonull2")) pli->do_null2 = FALSE;
  else                                           pli->do_null2 = TRUE;           
  
  pli->by_E            = TRUE;
  pli->E               = (go ? esl_opt_GetReal(go, "-E") : 10.0);
  pli->T               = 0.0;

  if (go && esl_opt_IsOn(go, "-T")) {
    pli->T    = esl_opt_GetReal(go, "-T");
    pli->by_E = FALSE;
  }

  pli->inc_by_E           = TRUE;
  pli->incE               = (go ? esl_opt_GetReal(go, "--incE") : 0.01);
  pli->incT               = 0.0;
  
  if (go && esl_opt_IsOn(go, "--incT")) {
    pli->incT     = esl_opt_GetReal(go, "--incT");
    pli->inc_by_E = FALSE;
  }

  pli->Z = 0.;

  pli->F1     = ((go && esl_opt_IsOn(go, "--F1")) ? ESL_MIN(1.0, esl_opt_GetReal(go, "--F1")) : 0.02);
  pli->F2     = (go ? ESL_MIN(1.0, esl_opt_GetReal(go, "--F2")) : 1e-3);
  pli->F3     = (go ? ESL_MIN(1.0, esl_opt_GetReal(go, "--F3")) : 1e-5);
  pli->S1     = (go ? ESL_MIN(1.0, esl_opt_GetReal(go, "--S1")) : 0.05);
  pli->S2     = (go ? ESL_MIN(1.0, esl_opt_GetReal(go, "--S2")) : 0.001); 
 
  if (go && esl_opt_GetBoolean(go, "--nobias"))  pli->do_biasfilter = FALSE;
  else                                            pli->do_biasfilter = TRUE;
  
  if (go && esl_opt_GetBoolean(go, "--max")) {
    pli->do_biasfilter = FALSE;
    pli->F1 = pli->F2 = pli->F3 = 1.0;
    pli->S1 = 0.3; // Must set some threshold for SSV windows
    pli->S2 = 1.0;
  }
  
  pli->signal_scores = NULL; 
  ESL_ALLOC(pli->signal_scores, sizeof(float) * p7S_SPLICE_SIGNALS);
  p7_splicepipeline_SignalScores(pli->signal_scores);  

  pli->nuc_sq   = NULL;
  pli->amino_sq = NULL;

  pli->orig_nuc_idx = NULL;

  pli->fwd = NULL;
  pli->bwd = NULL;
  if ((pli->fwd = p7_omx_Create(M_hint, L_hint, L_hint)) == NULL) goto ERROR;
  if ((pli->bwd = p7_omx_Create(M_hint, L_hint, L_hint)) == NULL) goto ERROR;

  pli->gfwd = NULL;
  pli->gbwd = NULL;

  if ((pli->gfwd = p7_gmx_fs_Create(M_hint, L_hint, L_hint, p7P_CODONS)) == NULL) goto ERROR;
  if ((pli->gbwd = p7_gmx_fs_Create(M_hint, L_hint, L_hint, 0         )) == NULL) goto ERROR;

  pli->vit = NULL;
  if ((pli->vit = p7_gmx_sp_Create(M_hint, L_hint, L_hint)) == NULL) goto ERROR;
  
  pli->sig_idx = p7_splicepipline_CreateIndex(M_hint, L_hint, L_hint);

  pli->bg = NULL;

  pli->hit = NULL;
  
  return pli;

  ERROR:
    p7_splicepipeline_Destroy(pli);
    return NULL;

}

/* Function:  p7_splicepipeline_Reuse()
 *
 * Purpose:   Reuse <pli> for next splice path
 * 
 */
void 
p7_splicepipeline_Reuse(SPLICE_PIPELINE *pli)
{

  esl_sq_Destroy(pli->nuc_sq);
  esl_sq_Destroy(pli->amino_sq);
  pli->nuc_sq = NULL;
  pli->amino_sq = NULL; 
 
  p7_omx_Reuse(pli->fwd);
  p7_omx_Reuse(pli->bwd);

  p7_gmx_Reuse(pli->gfwd);
  p7_gmx_Reuse(pli->gbwd);

  p7_gmx_Reuse(pli->vit);

  if(pli->orig_nuc_idx != NULL) free(pli->orig_nuc_idx);
  pli->orig_nuc_idx = NULL;  

  if(pli->hit != NULL && pli->hit->dcl != NULL) { 
    p7_alidisplay_Destroy(pli->hit->dcl->ad);
    p7_trace_splice_Destroy(pli->hit->dcl->tr);
    if(pli->hit->dcl->scores_per_pos != NULL)
      free(pli->hit->dcl->scores_per_pos);
  }
  if(pli->hit != NULL) {
    p7_hit_Destroy(pli->hit);
  }

  pli->hit = NULL;

 return;
  
}

/* Function: p7_splicepipeline_Destroy()
 *
 * Purpose:  Frees a <SPLICE_PIPELINE>
 */
void 
p7_splicepipeline_Destroy(SPLICE_PIPELINE *pli)
{

  if(pli == NULL) return;

  esl_sq_Destroy(pli->nuc_sq);
  esl_sq_Destroy(pli->amino_sq);

  if(pli->signal_scores != NULL) free(pli->signal_scores); 
  if(pli->orig_nuc_idx != NULL) free(pli->orig_nuc_idx);

  p7_omx_Destroy(pli->fwd);
  p7_omx_Destroy(pli->bwd);

  p7_gmx_Destroy(pli->gfwd);
  p7_gmx_Destroy(pli->gbwd);

  p7_gmx_Destroy(pli->vit);

  p7_splicepipeline_DestroyIndex(pli->sig_idx);

  p7_bg_Destroy(pli->bg);

  if(pli->hit != NULL && pli->hit->dcl != NULL) {
    p7_alidisplay_Destroy(pli->hit->dcl->ad);
    p7_trace_splice_Destroy(pli->hit->dcl->tr);
  }
  
  p7_hit_Destroy(pli->hit);

  free(pli);

 return;
  
}



/*****************************************************************
 * 2. The SPLICE_SITE_IDX structure.
 *****************************************************************/

/* Function:  p7_splicepipline_CreateIndex() 
 * Synopsis:  Allocates a splice site idx.
 *
 * Purpose:   Allocates a new <SPLICE_SITE_IDX> 
 *
 * Returns:   <SPLICE_SITE_IDX>  on success.
 *
 * Throws:    <NULL> on allocation error. 
 */
SPLICE_SITE_IDX*
p7_splicepipline_CreateIndex(int M_hint, int L_hint, int Lx_hint)
{
  int i,k;
  SPLICE_SITE_IDX *signal_sites;
  int status;

  signal_sites = NULL;
  ESL_ALLOC(signal_sites, sizeof(SPLICE_SITE_IDX));

  ESL_ALLOC(signal_sites->index_mem, sizeof(int)  * M_hint * SIGNAL_MEM_SIZE);
  ESL_ALLOC(signal_sites->index,     sizeof(int*) * M_hint);

  ESL_ALLOC(signal_sites->score_mem, sizeof(float)  * M_hint * SIGNAL_MEM_SIZE);
  ESL_ALLOC(signal_sites->score,     sizeof(float*) * M_hint);

  for(k = 0; k < M_hint; k++) {
    signal_sites->index[k] = signal_sites->index_mem + (k * SIGNAL_MEM_SIZE);
    signal_sites->score[k] = signal_sites->score_mem + (k * SIGNAL_MEM_SIZE);  
  }

  ESL_ALLOC(signal_sites->lookback_mem, sizeof(int)  * L_hint * M_hint);
  ESL_ALLOC(signal_sites->lookback,     sizeof(int*) * L_hint);

  for(i = 0; i < L_hint; i++) 
    signal_sites->lookback[i] = signal_sites->lookback_mem + (i * M_hint);
 
  ESL_ALLOC(signal_sites->parser_index,  sizeof(int)   * Lx_hint * p7S_PARSE_INDEX);
  ESL_ALLOC(signal_sites->parser_scores, sizeof(float) * Lx_hint * p7S_PARSE_SCORE);
   
  signal_sites->alloc_M  = M_hint;
  signal_sites->alloc_L  = L_hint;
  signal_sites->alloc_Lx = Lx_hint;

  return signal_sites;

  ERROR:
    p7_splicepipeline_DestroyIndex(signal_sites);
    return NULL;

}

/* Function:  p7_splicepipline_GrowIndex() 
 * Synopsis:  Grow a splice site idx.
 *
 * Purpose:   Allocates a larger <SPLICE_SITE_IDX>
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_splicepipline_GrowIndex(SPLICE_SITE_IDX *signal_sites, int M, int L, int mode)
{

  int i, k;
  int status;

  if(mode == ALIGNMENT_MODE && M <= signal_sites->alloc_M && L+1 <= signal_sites->alloc_L) return eslOK;
  if(mode == PARSER_MODE && M <= signal_sites->alloc_M && L+1 <= signal_sites->alloc_Lx) return eslOK;

  if( M > signal_sites->alloc_M) {

    ESL_REALLOC(signal_sites->index_mem, sizeof(int)  * M * SIGNAL_MEM_SIZE);
    ESL_REALLOC(signal_sites->index,     sizeof(int*) * M);

    ESL_REALLOC(signal_sites->score_mem, sizeof(float)  * M * SIGNAL_MEM_SIZE);
    ESL_REALLOC(signal_sites->score,     sizeof(float*) * M);  

    for(k = 0; k < M; k++) {
      signal_sites->index[k] = signal_sites->index_mem + (k * SIGNAL_MEM_SIZE);
      signal_sites->score[k] = signal_sites->score_mem + (k * SIGNAL_MEM_SIZE);
    }
  }

  if(mode == ALIGNMENT_MODE ) {

    if( (L+1) > signal_sites->alloc_L && M > signal_sites->alloc_M )  {
      ESL_REALLOC(signal_sites->lookback_mem, sizeof(int)  * (L+1) * M);
      ESL_REALLOC(signal_sites->lookback,     sizeof(int*) * (L+1));
   
      for(i = 0; i < L+1; i++)
        signal_sites->lookback[i] = signal_sites->lookback_mem + (i * M);
    }
    else if( (L+1) > signal_sites->alloc_L) {
      ESL_REALLOC(signal_sites->lookback_mem, sizeof(int)  * (L+1) * signal_sites->alloc_M);    
      ESL_REALLOC(signal_sites->lookback,     sizeof(int*) * (L+1));
  
      for(i = 0; i < L+1; i++)
        signal_sites->lookback[i] = signal_sites->lookback_mem + (i * signal_sites->alloc_M);
    }
    else if( M > signal_sites->alloc_M )  {
      ESL_REALLOC(signal_sites->lookback_mem, sizeof(int)  * signal_sites->alloc_L * M);
  
      for(i = 0; i < signal_sites->alloc_L; i++)
        signal_sites->lookback[i] = signal_sites->lookback_mem + (i * M);
    }
  
    signal_sites->alloc_M = ESL_MAX(M, signal_sites->alloc_M);
    signal_sites->alloc_L = ESL_MAX((L+1), signal_sites->alloc_L);
  
  }
  else {

    if(M > signal_sites->alloc_M) {
      ESL_REALLOC(signal_sites->lookback_mem, sizeof(int) * signal_sites->alloc_L * M);

      for(i = 0; i < signal_sites->alloc_L; i++)
        signal_sites->lookback[i] = signal_sites->lookback_mem + (i * M);
    }
    if( (L+1) > signal_sites->alloc_Lx) {
      ESL_REALLOC(signal_sites->parser_index,  sizeof(int)   * (L+1) * p7S_PARSE_INDEX);
      ESL_REALLOC(signal_sites->parser_scores, sizeof(float) * (L+1) * p7S_PARSE_SCORE);
      signal_sites->alloc_Lx = L+1;
    }
  }

  return eslOK;

  ERROR:
    p7_splicepipeline_DestroyIndex(signal_sites);
    return status;
}


/* Function:  p7_splicesiteidx_Destroy()
 *
* Purpose:  Frees a <SPLICE_SITE_IDX>
 */
void
p7_splicepipeline_DestroyIndex(SPLICE_SITE_IDX *signal_sites)
{

  if(signal_sites == NULL) return;

  if(signal_sites->index     != NULL) free(signal_sites->index);
  if(signal_sites->index_mem != NULL) free(signal_sites->index_mem);

  if(signal_sites->score     != NULL) free(signal_sites->score);
  if(signal_sites->score_mem != NULL) free(signal_sites->score_mem);

  if(signal_sites->lookback     != NULL) free(signal_sites->lookback);
  if(signal_sites->lookback_mem != NULL) free(signal_sites->lookback_mem);

  if(signal_sites->parser_index  != NULL) free(signal_sites->parser_index);
  if(signal_sites->parser_scores != NULL) free(signal_sites->parser_scores);

  free(signal_sites);
  signal_sites = NULL;

  return;
}
