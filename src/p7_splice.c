/* Create spliced alignments from P7_TOPHITS by
 *   1. Definning coords of target ranges in which spliceable hits exist
 *   2. Build splice graph of hits in each target range
 *   3. Find best path(s) through the splice graph(s)
 *   3. Fill gaps in splice graph(s) by searching for missing exons
 *   4. Realign spliced exons to model and compare score to unspliced scores
 *   5. If spliced alignent is more probable repalce unspliced hits
 *
 * Contents:
 *    1. Macros, structs and struct related functions
 *    2. Main function - SpliceHits() 
 *    3. Internal Routines
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


/* Splice singal probabilities taken from
 *  * "Comprehensive splice-site analysis using comparative genomics",
 *   * Nihar Sheth et al., 2006 */
void
p7_splice_SignalScores(float *f)
{
  f[0] = log(0.9919);     /* GT-AG */
  f[1] = log(0.0074);     /* GC-AG */
  f[2] = log(0.0007);     /* AT-AC */
  f[3] = -eslINFINITY; //log(0.0002);     /* OTHER */
  return;
}

/* Create a TARGET_RANGE object with room for nalloc P7_HIT pointers*/
TARGET_RANGE *
target_range_create(int nalloc)
{
  int status;

  TARGET_RANGE *target_range = NULL;
  ESL_ALLOC(target_range, sizeof(TARGET_RANGE));

  target_range->seqname = NULL;

  target_range->th = NULL;
  ESL_ALLOC(target_range->th, sizeof(P7_TOPHITS));
  target_range->th->is_sorted_by_seqidx = FALSE;

  target_range->th->hit    = NULL;
  ESL_ALLOC(target_range->th->hit,  nalloc * sizeof(P7_HIT*));

  target_range->reportable = NULL;
  ESL_ALLOC(target_range->reportable, nalloc * sizeof(int*));

  target_range->orig_hit_idx = NULL;
  ESL_ALLOC(target_range->orig_hit_idx, nalloc * sizeof(int*));

  target_range->in_target_range = NULL;
  ESL_ALLOC(target_range->in_target_range, nalloc * sizeof(int*));

  target_range->th->N      = 0;
  target_range->th->Nalloc = nalloc;
 
  target_range->orig_N = 0; 

  return target_range;

  ERROR:
    target_range_destroy(target_range);
    return NULL;
}

/* Double size of TARGET_RANGE object */
int
target_range_grow(TARGET_RANGE *target_range)
{
  P7_TOPHITS *th;
  int        status;
 
  th = target_range->th;
  if(th->N == th->Nalloc) {
    th->Nalloc *= 2;
    ESL_REALLOC(th->hit,                       sizeof(P7_HIT *) * th->Nalloc);
    ESL_REALLOC(target_range->in_target_range, sizeof(int)      * th->Nalloc); 
  }

  return eslOK;

  ERROR:
    target_range_destroy(target_range);
    return status;
}

/* Free a TARGET_RANGE object */
void
target_range_destroy(TARGET_RANGE *target_range)
{

  int i;
 
  if (target_range == NULL) return;

  target_range->seqname = NULL;

  for (i = target_range->orig_N; i < target_range->th->N; i++) {
    p7_trace_splice_Destroy(target_range->th->hit[i]->dcl->tr);
    p7_alidisplay_Destroy(target_range->th->hit[i]->dcl->ad);
    free(target_range->th->hit[i]->dcl->scores_per_pos);
    p7_hit_Destroy(target_range->th->hit[i]);
  }

  if (target_range->th->hit != NULL)
    free(target_range->th->hit); 
 
  if (target_range->th != NULL)
    free(target_range->th);

  if(target_range->reportable != NULL)
    free(target_range->reportable);

  if(target_range->orig_hit_idx != NULL)
    free(target_range->orig_hit_idx); 

  if(target_range->in_target_range != NULL)
    free(target_range->in_target_range);

  free(target_range);

  return;

}


/* Create a SPLICE_EDGE object */
SPLICE_EDGE *
splice_edge_create(void)
{

  int status;
  SPLICE_EDGE *edge;

  edge = NULL;
  ESL_ALLOC(edge, sizeof(SPLICE_EDGE));

  edge->frameshift   = FALSE;

  edge->splice_score = -eslINFINITY;
  edge->signal_score = -eslINFINITY;

  edge->prev = NULL;
  edge->next = NULL;

  return edge;

  ERROR:
    if (edge != NULL) free(edge);
    return NULL;
}


/* Create a SPLICE_GRAPH object */
SPLICE_GRAPH *
splice_graph_create(void)
{
  int status;
  SPLICE_GRAPH *graph;

  ESL_ALLOC(graph, sizeof(SPLICE_GRAPH));

  graph->th = NULL;
  ESL_ALLOC(graph->th, sizeof(P7_TOPHITS));
  graph->th->is_sorted_by_seqidx = FALSE;

  graph->nalloc         = 0;
  graph->num_nodes      = 0;
  graph->num_edges      = 0;
  graph->orig_N         = 0;
  graph->split_N        = 0;
  graph->th->N          = 0; 

  graph->reportable     = NULL;
  graph->orig_hit_idx   = NULL;
  graph->in_graph       = NULL;
  
  graph->out_edge_cnt  = NULL;
  graph->in_edge_cnt   = NULL;
  graph->best_out_edge = NULL;
  graph->best_in_edge  = NULL;

  graph->path_scores   = NULL;
  graph->hit_scores    = NULL;

  graph->th->hit       = NULL;
  graph->ds_nodes      = NULL;
  graph->edges         = NULL;

  return graph;

ERROR:
    splice_graph_destroy(graph);
    return NULL;

}

/* Free a SPLICE_GRAPH and all it nodes and edges */
void splice_graph_destroy
(SPLICE_GRAPH *graph)
{

  int i;
  SPLICE_NODE    *tmp_node;
  SPLICE_EDGE    *tmp_edge;

  if (graph == NULL) return;

  if(graph->reportable      != NULL) free(graph->reportable);
  if(graph->orig_hit_idx    != NULL) free(graph->orig_hit_idx);
  if(graph->in_graph        != NULL) free(graph->in_graph);

  if(graph->path_scores     != NULL) free(graph->path_scores);
  if(graph->hit_scores      != NULL) free(graph->hit_scores);

  if(graph->out_edge_cnt    != NULL) free(graph->out_edge_cnt);
  if(graph->in_edge_cnt     != NULL) free(graph->in_edge_cnt);
  if(graph->best_out_edge   != NULL) free(graph->best_out_edge);
  if(graph->best_in_edge    != NULL) free(graph->best_in_edge);

 
  for (i = graph->orig_N; i < graph->th->N; i++) {
    p7_trace_splice_Destroy(graph->th->hit[i]->dcl->tr);
    p7_alidisplay_Destroy(graph->th->hit[i]->dcl->ad);
    free(graph->th->hit[i]->dcl->scores_per_pos);
    p7_hit_Destroy(graph->th->hit[i]);
  }

  if (graph->th->hit != NULL) free(graph->th->hit);
  if (graph->th      != NULL) free(graph->th);

  for(i = 0; i < graph->num_nodes; i++) {
    while(graph->ds_nodes[i] != NULL) {
      tmp_node = graph->ds_nodes[i];
      graph->ds_nodes[i] = tmp_node->next;
      free(tmp_node);
    }
    while(graph->edges[i] != NULL) {
      tmp_edge = graph->edges[i];
      graph->edges[i] = tmp_edge->next;
      
      free(tmp_edge);
    }
  }

  if(graph->ds_nodes != NULL) free(graph->ds_nodes);
  if(graph->edges    != NULL) free(graph->edges);

  graph->seqname = NULL;

  free(graph);
  graph = NULL;

  return;
}

/* Alloacte room for num_nodes nodes and their edges in a SPLICE_GRAPH */
int
splice_graph_create_nodes(SPLICE_GRAPH *graph, int num_nodes)
{

  int i;
  int status;

  if(graph == NULL)  graph = splice_graph_create();

  graph->nalloc         = num_nodes;
  graph->num_nodes      = num_nodes;

  //Allocate array of hit scores
  ESL_ALLOC(graph->hit_scores,        sizeof(float)        * graph->nalloc);

  // Allocate arrays of path scores
  ESL_ALLOC(graph->path_scores,    sizeof(float)        * graph->nalloc);

  // Allocate arrays for keeping track of which nodes have edges ans which edge is best
  ESL_ALLOC(graph->out_edge_cnt,   sizeof(int)          * graph->nalloc);
  ESL_ALLOC(graph->in_edge_cnt,    sizeof(int)          * graph->nalloc);
  ESL_ALLOC(graph->best_out_edge,  sizeof(int)          * graph->nalloc);
  ESL_ALLOC(graph->best_in_edge,   sizeof(int)          * graph->nalloc);

  ESL_ALLOC(graph->ds_nodes,           sizeof(SPLICE_NODE*) * graph->nalloc);
  for(i = 0; i < graph->num_nodes; i++)
    graph->ds_nodes[i] = NULL;

  ESL_ALLOC(graph->edges,             sizeof(SPLICE_EDGE*) * graph->nalloc);
  for(i = 0; i < graph->num_nodes; i++)
    graph->edges[i] = NULL;



  return eslOK;

  ERROR:
    splice_graph_destroy(graph);
    return status;
}

int
splice_graph_grow(SPLICE_GRAPH *graph)
{
  int status;

  if(graph->num_nodes == graph->nalloc) {

    graph->nalloc *= 2;

    ESL_REALLOC(graph->hit_scores,      sizeof(float)        * graph->nalloc);
    ESL_REALLOC(graph->path_scores,     sizeof(float)        * graph->nalloc); 
    ESL_REALLOC(graph->out_edge_cnt,    sizeof(int)          * graph->nalloc);
    ESL_REALLOC(graph->in_edge_cnt,     sizeof(int)          * graph->nalloc);
    ESL_REALLOC(graph->best_out_edge,   sizeof(int)          * graph->nalloc);
    ESL_REALLOC(graph->best_in_edge,    sizeof(int)          * graph->nalloc);

    ESL_REALLOC(graph->ds_nodes,        sizeof(SPLICE_NODE*) * graph->nalloc);
    ESL_REALLOC(graph->edges,           sizeof(SPLICE_EDGE*) * graph->nalloc);
  }

  return eslOK;

  ERROR:
    splice_graph_destroy(graph);
    return status;

}



SPLICE_PATH*
splice_path_create(int path_len)
{

  SPLICE_PATH *path;
  int          status;

  path = NULL;
  ESL_ALLOC(path, sizeof(SPLICE_PATH));

  path->alloc_len = path_len*2; 
  path->path_len  = path_len;

  path->node_id = NULL;
  ESL_ALLOC(path->node_id,                        sizeof(int)*(path_len*2));
  path->split = NULL;
  ESL_ALLOC(path->split,                          sizeof(int)*(path_len*2));

  path->upstream_spliced_amino_end     = NULL;
  path->downstream_spliced_amino_start = NULL;
  ESL_ALLOC(path->upstream_spliced_amino_end,     sizeof(int)*(path_len*2+1));
  ESL_ALLOC(path->downstream_spliced_amino_start, sizeof(int)*(path_len*2+1));

  path->upstream_spliced_nuc_end     = NULL;
  path->downstream_spliced_nuc_start = NULL;
  ESL_ALLOC(path->upstream_spliced_nuc_end,       sizeof(int)*(path_len*2+1));
  ESL_ALLOC(path->downstream_spliced_nuc_start,   sizeof(int)*(path_len*2+1));

  path->hit_scores = NULL;
  ESL_ALLOC(path->hit_scores, sizeof(float)*path_len*2);

  path->edge_scores = NULL;
  ESL_ALLOC(path->edge_scores, sizeof(float)*path_len*2);

  path->signal_scores = NULL;
  ESL_ALLOC(path->signal_scores, sizeof(float)*path_len*2);

  path->hits = NULL;
  ESL_ALLOC(path->hits,          sizeof(P7_HIT*)*path_len*2);

  return path;

  ERROR:
    splice_path_destroy(path);
    return NULL;
}

int
splice_path_grow(SPLICE_PATH *path)
{

  int status;
  
  if(path->path_len < path->alloc_len) return eslOK;

  path->alloc_len *= 2;

  ESL_REALLOC(path->node_id,                        sizeof(int)     * path->alloc_len);
  ESL_REALLOC(path->split,                          sizeof(int)     * path->alloc_len);
  ESL_REALLOC(path->upstream_spliced_amino_end,     sizeof(int)     * (path->alloc_len+1));
  ESL_REALLOC(path->downstream_spliced_amino_start, sizeof(int)     * (path->alloc_len+1));
  ESL_REALLOC(path->upstream_spliced_nuc_end,       sizeof(int)     * (path->alloc_len+1));
  ESL_REALLOC(path->downstream_spliced_nuc_start,   sizeof(int)     * (path->alloc_len+1));
  ESL_REALLOC(path->hit_scores,                     sizeof(float)   * path->alloc_len);
  ESL_REALLOC(path->edge_scores,                    sizeof(float)   * path->alloc_len);
  ESL_REALLOC(path->signal_scores,                  sizeof(float)   * path->alloc_len);
  ESL_REALLOC(path->hits,                           sizeof(P7_HIT*) * path->alloc_len);   

  return eslOK;

  ERROR:
    splice_path_destroy(path);
    return status;

}


void
splice_path_destroy(SPLICE_PATH *path)
{

   int i;

   if(path == NULL) return;

   /* Destroy split hits */
   for(i = 0; i < path->path_len; i++) {
     if(path->split[i]) {
       
       p7_alidisplay_Destroy(path->hits[i]->dcl->ad);
       p7_trace_fs_Destroy(path->hits[i]->dcl->tr);
       free(path->hits[i]->dcl->scores_per_pos);
       p7_hit_Destroy(path->hits[i]);       
     }
   }

   if(path->node_id != NULL) free(path->node_id);
   if(path->split != NULL)   free(path->split);

   if(path->upstream_spliced_amino_end     != NULL)
     free(path->upstream_spliced_amino_end);
   if(path->downstream_spliced_amino_start != NULL)
     free(path->downstream_spliced_amino_start);

   if(path->upstream_spliced_nuc_end       != NULL)
     free(path->upstream_spliced_nuc_end);
   if(path->downstream_spliced_nuc_start   != NULL)
     free(path->downstream_spliced_nuc_start);

   if(path->hit_scores    != NULL) free(path->hit_scores);
   if(path->edge_scores   != NULL) free(path->edge_scores);
   if(path->signal_scores != NULL) free(path->signal_scores);
   if(path->hits          != NULL) free(path->hits);
   if(path                != NULL) free(path);

   return;
}


SPLICE_PIPELINE* 
splice_pipeline_create(const ESL_GETOPTS *go, int M_hint, int L_hint)
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
  
  if (go && esl_opt_GetBoolean(go, "--nobias"))  pli->do_biasfilter = FALSE;
  else                                            pli->do_biasfilter = TRUE;
  
  if (go && esl_opt_GetBoolean(go, "--max")) {
    pli->do_biasfilter = FALSE;
    pli->F1 = pli->F2 = pli->F3 = 1.0;
  }
   

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

  pli->bg = NULL;

  pli->hit = NULL;
  
  return pli;

  ERROR:
    splice_pipeline_destroy(pli);
    return NULL;

}


void 
splice_pipeline_reuse(SPLICE_PIPELINE *pli)
{

  esl_sq_Destroy(pli->nuc_sq);
  esl_sq_Destroy(pli->amino_sq);
  pli->nuc_sq = NULL;
  pli->amino_sq = NULL; 
 
  p7_omx_Reuse(pli->fwd);
  p7_omx_Reuse(pli->bwd);

  p7_gmx_Reuse(pli->gfwd);
  p7_gmx_Reuse(pli->gbwd);

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



void 
splice_pipeline_destroy(SPLICE_PIPELINE *pli)
{

  if(pli == NULL) return;

  esl_sq_Destroy(pli->nuc_sq);
  esl_sq_Destroy(pli->amino_sq);
  
  if(pli->orig_nuc_idx != NULL) free(pli->orig_nuc_idx);

  p7_omx_Destroy(pli->fwd);
  p7_omx_Destroy(pli->bwd);

  p7_gmx_Destroy(pli->gfwd);
  p7_gmx_Destroy(pli->gbwd);

  p7_bg_Destroy(pli->bg);

  if(pli->hit != NULL && pli->hit->dcl != NULL) {
    p7_alidisplay_Destroy(pli->hit->dcl->ad);
    p7_trace_splice_Destroy(pli->hit->dcl->tr);
  }
  
  p7_hit_Destroy(pli->hit);

  free(pli);

 return;
  
}



/*  Function: p7_splice_SpliceHits
 *  Synopsis: SPLASH's splcing pipeline
 *
 *  Purpose : Run the splicing pipeline on a collections of hits 
 *            <tophits> between a protein query model <gm> and a 
 *            nucleotide sequence from the target file <seq_file>.
 *
 * Returns:   <eslOK> on success. If hits are successfully splice
 *            the new spliced alignements, scores, etc will be
 *            included in the <tophits> to be reported.
 *         
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_splice_SpliceHits(P7_TOPHITS *tophits, P7_HMM *hmm, P7_OPROFILE *om, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, P7_SCOREDATA *scoredata, ESL_GETOPTS *go, ESL_GENCODE *gcode, ESL_SQFILE *seq_file, FILE *ofp, int64_t db_nuc_cnt)
{

  int              i;
  int              hit_cnt;
  int              range_cnt;
  int              success;
  int              num_hits_processed; 
  int              prev_num_hits_processed;
  int              is_sorted_by_sortkey;
  int              revcomp, prev_revcomp;
  int              frameshift;
  int              seq_min, seq_max;
  int64_t          seqidx, prev_seqidx;
  int             *hits_processed;
  int             *hit_spliced;
  int64_t         *range_bound_mins;
  int64_t         *range_bound_maxs;
  SPLICE_PIPELINE *pli;
  TARGET_RANGE    *target_range;
  SPLICE_GRAPH   *graph;
  SPLICE_PATH     *path;
  ESL_SQ          *path_seq;
  P7_HIT          *hit; 
  int              status;

  target_range     = NULL;
  graph           = NULL;
  path             = NULL;
  hits_processed   = NULL;
  hit_spliced      = NULL;
  range_bound_mins = NULL;
  range_bound_maxs = NULL;
  path_seq         = NULL;
 
  /* Sort hits by sequence, strand and nucleotide positon. 
   * Note if P7_TOPHITS was perviously sorted by sortkey 
   * so we can return it to it's original state */  
  is_sorted_by_sortkey = tophits->is_sorted_by_sortkey;
  if ((status = p7_tophits_SortBySeqidxAndAlipos(tophits)) != eslOK) goto ERROR;
     
  pli = splice_pipeline_create(go, om->M, om->M * 3);
  pli->bg = p7_bg_Create(om->abc);
  if (pli->do_biasfilter) p7_bg_SetFilter(pli->bg, om->M, om->compo);  

  /* Arrays to keep track of hit status */
  ESL_ALLOC(hits_processed,   tophits->N*sizeof(int));
  ESL_ALLOC(hit_spliced,      tophits->N*sizeof(int));

  /* Check for any hits that are duplicates or have low quailty alignment and mark them as processed */
  hit_cnt = 0;
  num_hits_processed = 0;
  for(i = 0; i < tophits->N; i++) {
    tophits->hit[i]->dcl->ad->exon_cnt = 1;
    if ((tophits->hit[i]->flags & p7_IS_DUPLICATE)) {
      hits_processed[i] = 1;
      num_hits_processed++;
    }
    else { 
      hits_processed[i] = 0;
      hit_cnt++;
    }  
    hit_spliced[i]   = 0;
  }
  
  /* Arrays for keeping track of allowable target range mins 
   * and maxes set by previous ranges or splicings */
  ESL_ALLOC(range_bound_mins, hit_cnt * sizeof(int64_t));
  ESL_ALLOC(range_bound_maxs, hit_cnt * sizeof(int64_t));
  
  //prev_target_range = NULL;
  prev_num_hits_processed = -1;
  range_cnt = 0;

  prev_revcomp = -1;
  prev_seqidx  = -1;
  printf("\nQuery %s LENG %d \n",  gm->name, gm->M);
   fflush(stdout);
  /* loop through until all hits have been processed */
  while(num_hits_processed < tophits->N) {
    
    if(prev_num_hits_processed == num_hits_processed)
      ESL_XEXCEPTION(eslFAIL, "p7_splice_SpliceHits : loop failed to process hits");
    prev_num_hits_processed = num_hits_processed;
  
    /* Find the first unprocessed hit (sorting will ensure that hits 
     * will be on same sequence and strand untill all are processed) */
    i = 0;
    while (i < tophits->N && hits_processed[i]) i++;
    hit = tophits->hit[i];
    seqidx = hit->seqidx;
    revcomp = 1;
    if (hit->dcl->iali < hit->dcl->jali)
      revcomp = 0;

    /* If we have a new sequence and strand we build a new target range and graph */
    if(seqidx != prev_seqidx || revcomp != prev_revcomp) {   

      /* unless this is the first target_range destroy all the old data */
      if( prev_seqidx != -1) {
        target_range_destroy(target_range);
        target_range = NULL;
  
        splice_graph_destroy(graph);
        graph = NULL;
        range_cnt = 0;
      }   

      //if ((graph = splice_graph_create()) == NULL) goto ERROR;

      //add_hits_to_graph(graph, tophits, hits_processed, seqidx, revcomp);

      if ((target_range = target_range_create(tophits->N)) == NULL) goto ERROR;
      build_target_range(target_range, tophits, hits_processed, seqidx, revcomp);


      printf("\nQuery %s Target %s strand %c\n", gm->name, target_range->seqname, (target_range->revcomp ? '-' : '+'));
      fflush(stdout);

     //target_range_dump(stdout, target_range, TRUE);
 
      // if ((graph = splice_graph_create()) == NULL) goto ERROR; 
      if ((graph = splice_graph_create()) == NULL) goto ERROR; 
      
      if ((status = fill_graph_with_nodes(graph, target_range, gm, hmm, pli->bg, gcode, seq_file)) != eslOK) goto ERROR;

      //target_range_dump(stdout, target_range, TRUE);
      
      if((status = longest_path_upstream(target_range, graph)) != eslOK) goto ERROR;

      if ((status = fill_holes_in_graph(target_range, graph, gm, hmm, pli->bg, gcode, seq_file)) != eslOK) goto ERROR;
      target_range_dump(stdout, target_range, TRUE);
       //ds_graph_dump(stdout, target_range, graph); 
     //graph_dump(stdout, target_range, graph, FALSE);
      //check_for_bypasses(target_range, graph);

        
    }

  //  target_range_dump(stdout, target_range, TRUE); 
//graph_dump(stdout, target_range, graph, FALSE);    
 
    enforce_range_bounds(graph, target_range->th, range_bound_mins, range_bound_maxs, range_cnt);

    //graph_dump(stdout, target_range, graph, FALSE);

    reset_edge_counts(target_range, graph);

  //graph_dump(stdout, target_range, graph, FALSE);
    
    path = evaluate_paths(target_range, graph);
   //graph_dump(stdout, target_range, graph, FALSE); 
    path_dump(stdout, path);
   
    seq_min = ESL_MIN(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]);
    seq_max = ESL_MAX(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]);
    path_seq = get_sub_sequence(seq_file, target_range->seqname, seq_min, seq_max, path->revcomp);

    
    split_hits_in_path2 (target_range, graph, path, gm, hmm, pli->bg, gcode, path_seq);
    path_dump(stdout, path);

    if(path->path_len > 1) {
      frameshift = FALSE;
      for(i = 0; i < path->path_len; i++) 
        if(path->hits[i]->dcl->tr->fs != 0) frameshift = TRUE;

      if(!frameshift)
        splice_the_path(target_range,path, pli, tophits, om, gm, gcode, path_seq, db_nuc_cnt, &frameshift, &success);
          
      if(frameshift)
        splice_the_path_frameshift (target_range, path, pli, tophits, gm_fs, om, gcode, path_seq, db_nuc_cnt, &success);
    }
    else
      success = FALSE;    
   
    /* Reset the range_bounds around the spliced hit and set all hits from 
     * the target_range that are outside this new range to unprocessed */
    if (success) {
     
      range_bound_mins[range_cnt] = ESL_MIN(pli->hit->dcl->iali, pli->hit->dcl->jali);
      range_bound_maxs[range_cnt] = ESL_MAX(pli->hit->dcl->iali, pli->hit->dcl->jali);
  
      release_hits_from_target_range(target_range, hits_processed, &num_hits_processed, range_bound_mins[range_cnt], range_bound_maxs[range_cnt]);           

      pli->hit->dcl = NULL;
    } 
	else {
      
      /* If the path failed but there are other hits in the graph that were not in the path 
       * we want to release them to be considered in another target range */
      range_bound_mins[range_cnt] = ESL_MIN(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]);
      range_bound_maxs[range_cnt] = ESL_MAX(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]);
 
      release_hits_from_target_range(target_range, hits_processed, &num_hits_processed, range_bound_mins[range_cnt], range_bound_maxs[range_cnt]);  
 
    }
    
    range_cnt++; 
    prev_seqidx = seqidx;
    prev_revcomp = revcomp; 
  
    if(path  != NULL) splice_path_destroy(path);
    if(path_seq != NULL) esl_sq_Destroy(path_seq);
    path_seq = NULL;
    splice_pipeline_reuse(pli);

  }
 
  /* Leave only footprints */
  if (is_sorted_by_sortkey)
    if ((status = p7_tophits_SortBySortkey(tophits)) != eslOK) goto ERROR;   

  if(pli          != NULL) splice_pipeline_destroy(pli);
  if(target_range != NULL) target_range_destroy(target_range);
  if(graph        != NULL) splice_graph_destroy(graph);

  if (hits_processed   != NULL) free(hits_processed);
  if (hit_spliced      != NULL) free(hit_spliced);
  if (range_bound_mins != NULL) free(range_bound_mins);
  if (range_bound_maxs != NULL) free(range_bound_maxs);


  return status;

  ERROR:
    if(pli          != NULL) splice_pipeline_destroy(pli);
    if(target_range != NULL) target_range_destroy(target_range);
    if(graph        != NULL) splice_graph_destroy(graph);
    if(path         != NULL) splice_path_destroy(path);
    if (hits_processed   != NULL) free(hits_processed);
    if (hit_spliced      != NULL) free(hit_spliced); 
    if (range_bound_mins != NULL) free(range_bound_mins);
    if (range_bound_maxs != NULL) free(range_bound_maxs);
    return status;
}


int
is_upstream(SPLICE_GRAPH* graph, int up_node, int down_node) {

  SPLICE_NODE *tmp_node;

  tmp_node = graph->ds_nodes[up_node];

  while(tmp_node != NULL && tmp_node->node_id != down_node)
    tmp_node = tmp_node->next;

  if(tmp_node == NULL) return FALSE;
  else                 return TRUE;

}

SPLICE_NODE*
get_node(SPLICE_GRAPH* graph, int up_node, int down_node) {

  SPLICE_NODE *node;

  node = graph->ds_nodes[up_node];

  while(node != NULL && node->node_id != down_node)
    node = node->next;

  return node;
  
}

int
edge_exists(SPLICE_GRAPH* graph, int up_node, int down_node) {

  SPLICE_EDGE *tmp_edge;

  tmp_edge = graph->edges[up_node];

  while(tmp_edge != NULL && tmp_edge->downstream_node_id != down_node)
    tmp_edge = tmp_edge->next;

  if(tmp_edge == NULL) return FALSE;
  else                 return TRUE;

}

SPLICE_EDGE*
get_edge(SPLICE_GRAPH* graph, int up_node, int down_node) {

  SPLICE_EDGE *edge;

  edge = graph->edges[up_node];

  while(edge != NULL && edge->downstream_node_id != down_node)
    edge = edge->next;

  return edge;

}

 /* Function: build_target_range
 *
 *  Synopsis: Defines a target range for splicing
 * 
 *  Purpose:  Identify a range of coodinates on a single sequence and 
 *            strand in which there are hits from P7_TOPHITS which may
 *            be splice compatible. 
 *
 *  Returns:  a pointer to the new <TARGET_RANGE> object.
 *  
 *  Throws:   <NULL> on allocation failure.
 */
int *
build_target_range (TARGET_RANGE *target_range, const P7_TOPHITS *tophits, int *hits_processed, int64_t seqidx, int revcomp)
{

  int           i;
  P7_HIT       *curr_hit;

  /*Add all hits from current sequence and strand to TARGET_RANGE */
  for (i = 0; i < tophits->N; i++) {

    curr_hit = tophits->hit[i];

    /* skip hits that are on a different sequence */
    if (curr_hit->seqidx != seqidx) continue;

    /* skip hits that are on a different strand */
    if (revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
    if ((!revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;
  
    target_range->reportable[i] = FALSE;
    if( (!hits_processed[i]) && (curr_hit->flags & p7_IS_REPORTED)) {
      target_range->reportable[i] = TRUE;
    }

    target_range->th->hit[target_range->th->N] = curr_hit;
    target_range->orig_hit_idx[target_range->th->N] = i;

    if( hits_processed[i] ) 
      target_range->in_target_range[target_range->th->N] = FALSE;
    else
      target_range->in_target_range[target_range->th->N] = TRUE;

    target_range->th->N++;

  }

  target_range->orig_N  = target_range->th->N;
  target_range->revcomp = revcomp;
  target_range->seqidx  = seqidx;
  target_range->seqname = target_range->th->hit[0]->name;

  return eslOK;
 
}




/* Function: get_sub_sequence
 *
 *  Synopsis: Fetch a target subsequence based on <target_range> coords
 *
 *  Purpose:  Fetch the sub-sequence dffined by <target_range> from the
 *            sequence file <seq_file>.  An ssi index for the <seq_file>
 *            must exist. 
 *
 *  Returns:  a pointer to an <ESL_SQ> representing the sub-sequence.
 *
 *  Throws:   <NULL> on allocation failure.
 */
ESL_SQ* 
get_sub_sequence(const ESL_SQFILE *seq_file, char* seqname, int64_t seq_min, int64_t seq_max, int revcomp)
{

  ESL_SQ     *target_seq;
  ESL_SQFILE *tmp_file;

  /* Open seq file and ssi */
  esl_sqfile_Open(seq_file->filename,seq_file->format,NULL,&tmp_file);
  esl_sqfile_OpenSSI(tmp_file,NULL);

  /* Get basic sequence info */
  target_seq = esl_sq_Create();
  esl_sqio_FetchInfo(tmp_file, seqname, target_seq);

  target_seq->start = seq_min;
  target_seq->end   = seq_max;
  target_seq->abc   = seq_file->abc;

  /* Make sure target range coords did not extend too far */
  if (target_seq->start < 1)
    target_seq->start = 1; 

  if (target_seq->end > target_seq->L)
    target_seq->end = target_seq->L;

  /* Fetch target range sequcene */
  if (esl_sqio_FetchSubseq(tmp_file,target_seq->name,target_seq->start,target_seq->end,target_seq) != eslOK) 
    esl_fatal(esl_sqfile_GetErrorBuf(tmp_file));

  esl_sqfile_Close(tmp_file);

  esl_sq_SetName(target_seq, seqname);
 
  if(revcomp)
   esl_sq_ReverseComplement(target_seq);

  esl_sq_Digitize(target_seq->abc, target_seq);

  return target_seq;
}


/* Function: release_hits_from_target_range
 *
 *  Synopsis: Mark as unprocessed some hits in a TARGET_RANGE 
 *
 *  Purpose: after the hits in a TARGET_RANGE have been spliced release 
 *           (mark as unprocessed) any hits in the TARGET_RANGE that 
 *           are outside the bounds of the spliced hit.
 *
 *  Returns: eslOK.
 *
 */

int
release_hits_from_target_range(const TARGET_RANGE *target_range, int *hits_processed, int *num_hits_processed, int range_bound_min, int range_bound_max)
{

  int        i;
  int        hit_min, hit_max;
  int        overlap_min, overlap_max;
  int        overlap_len;
  int        hit_idx;
  int        num_hits;
  P7_TOPHITS *th;
  P7_HIT     *hit;

  th = target_range->th;
  num_hits = *num_hits_processed;
  
  for(i = 0; i < th->N; i++) {
    hit = th->hit[i];

    hit_min = ESL_MIN(hit->dcl->iali, hit->dcl->jali);
    hit_max = ESL_MAX(hit->dcl->iali, hit->dcl->jali);
    overlap_min = ESL_MAX(hit_min, range_bound_min);
    overlap_max = ESL_MIN(hit_max, range_bound_max);
    overlap_len = overlap_max - overlap_min + 1;
    
    if(overlap_len > 0) {
      target_range->in_target_range[i] = FALSE;
      if (i < target_range->orig_N) {
        hit_idx =  target_range->orig_hit_idx[i];
        if(!hits_processed[hit_idx]) {
          hits_processed[hit_idx] = 1;
          num_hits++;
        }
      }
    }

  }
  
  *num_hits_processed = num_hits;
  return eslOK;

}          





/*  Function: fill_graph_with_nodes
 *
 *  Synopsis: Use hits in <th> to populate a SPLICE_GRAPH with nodes
 *
 *  Purpose:  Allocate space for <th->N> nodes in <graph> and populate
 *            with basic info on each node.  
 *
 *  Returns:  <eslOK> on success.
 *
 *  Throws:   <eslEMEM> on allocation failure.
 */
int
fill_graph_with_nodes(SPLICE_GRAPH *graph, TARGET_RANGE *target_range, const P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, ESL_GENCODE *gcode, const ESL_SQFILE *seq_file)
{

  int          up,down;
  int          intron_len;
  int          up_min, up_max;
  int          down_min, down_max;
  int          seq_min, seq_max;
  P7_TOPHITS  *th;
  SPLICE_NODE *new_node;
  SPLICE_NODE *tmp_node;
  SPLICE_EDGE *edge;
  ESL_SQ      *target_seq;
  int          status;
 
  target_seq = NULL;
  
  th = target_range->th;

  graph->revcomp = target_range->revcomp;

  /* Initialize nodes and add to graph */
  status = splice_graph_create_nodes(graph, th->N);

  if (status != eslOK) return status;

  /* Set downstreamness of nodes */
  for(up = 0; up < th->N; up++) {

    if(!target_range->in_target_range[up]) continue;
    for(down = 0; down < th->N; down++) {

      if(!target_range->in_target_range[down]) continue;
    
      /* determine if up hit is upstream of down hit */
      if(up == down) continue;

      /*If the hmm coords of up strart after down starts and end after down ends then down is not dowstream of up */
      if(th->hit[up]->dcl->ihmm >= th->hit[down]->dcl->ihmm &&
         th->hit[up]->dcl->jhmm >= th->hit[down]->dcl->jhmm)
        continue;

     /* I the seq coords of down don't start after the seq coords of up end then down is not dowstream of up */
      if (( graph->revcomp  && th->hit[down]->dcl->iali > th->hit[up]->dcl->jali) ||
         ((!graph->revcomp) && th->hit[up]->dcl->jali > th->hit[down]->dcl->iali))
         continue;

      ESL_ALLOC(new_node, sizeof(SPLICE_NODE));
      new_node->node_id = down;
      new_node->gap_checked  = FALSE;
      new_node->path_exists  = FALSE;
      new_node->path_checked = FALSE;
      new_node->next = NULL;

      tmp_node = graph->ds_nodes[up];
      if(tmp_node == NULL)
        graph->ds_nodes[up] = new_node;
      else {
        while(tmp_node->next != NULL) {
          tmp_node = tmp_node->next;
        }
        tmp_node->next = new_node;
      }

    }

    graph->hit_scores[up]  = th->hit[up]->dcl->aliscore;
    graph->path_scores[up] = -eslINFINITY;
    graph->out_edge_cnt[up]  = 0;
    graph->in_edge_cnt[up]   = 0;
    graph->best_out_edge[up] = -1;
    graph->best_in_edge[up]  = -1;
  }

  /* Connect with edges */
  for(up = 0; up < th->N; up++) {
    if(!target_range->in_target_range[up]) continue;
    for(down = 0; down < th->N; down++) {
      if(!target_range->in_target_range[down]) continue;
      /* check than up is upstream of down */
      if(!is_upstream(graph,up,down)) continue;
      
      /* check that hmm and ali coords are spliceable */
      if(th->hit[up]->dcl->jhmm + MAX_AMINO_EXT < th->hit[down]->dcl->ihmm)
        continue;
     
      intron_len = abs(th->hit[down]->dcl->iali - th->hit[up]->dcl->jali) + 1;
      if(intron_len < 1 || intron_len > MAX_INTRON_LEN)
        continue;
      
      /* get min and max coords and fetch sub sequence */      
      up_min = ESL_MIN(th->hit[up]->dcl->iali, th->hit[up]->dcl->jali);
      up_max = ESL_MAX(th->hit[up]->dcl->iali, th->hit[up]->dcl->jali);
 
      down_min = ESL_MIN(th->hit[down]->dcl->iali, th->hit[down]->dcl->jali);
      down_max = ESL_MAX(th->hit[down]->dcl->iali, th->hit[down]->dcl->jali);

      seq_min = ESL_MIN(up_min, down_min);
      seq_max = ESL_MAX(up_max, down_max);

      target_seq = get_sub_sequence(seq_file, target_range->seqname, seq_min, seq_max, graph->revcomp);
      
      /* Find edge */
      edge = connect_nodes_with_edges(th->hit[up], th->hit[down], gm, hmm, bg, gcode, target_seq, graph->revcomp);
      if(edge != NULL) {
        
        edge->upstream_node_id    = up;
        edge->downstream_node_id  = down;
        edge->target_seq_start   = target_seq->start;
        edge->target_seq_end     = target_seq->end;
        edge->target_seq_n       = target_seq->n;
        add_edge_to_graph(graph, edge);
      }
      if ( target_seq != NULL) esl_sq_Destroy(target_seq);
      target_seq = NULL;
    }
  }

  if ( target_seq != NULL) esl_sq_Destroy(target_seq);
  return eslOK;

  ERROR:
    if ( target_seq != NULL) esl_sq_Destroy(target_seq);
    return status;
}

/*  Function: connect_nodes_with_edges
 *
 *  Synopsis: Find and score all viable edges between nodes (hits) in <graph>
 *
 *  Purpose:  For each pair of nodes (hits in <th>) see if they are splice 
 *            compatiable. If so, score the splicing. If there is a viable 
 *            splicing add edge to graph. Edges are directed downstream 
 *            (from nodes with lower hmm coords to nodes with higher hmm coords).
 *
 *  Returns:  <eslOK> on success.
 *
 *  Throws:   <eslEMEM> on allocation failure.
 */

SPLICE_EDGE*
connect_nodes_with_edges(const P7_HIT *upstream_hit, const P7_HIT *downstream_hit, const P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, const ESL_GENCODE *gcode, const ESL_SQ *target_seq, int revcomp)
{

  int          up_amino_start, up_amino_end;
  int          down_amino_start, down_amino_end;
  int          intron_len;
  int          overlap;
  int          num_ext_aminos;
  int          min_hit_len;
  SPLICE_EDGE *edge;
  int          status;

  up_amino_start = upstream_hit->dcl->ihmm;
  up_amino_end   = upstream_hit->dcl->jhmm;

  down_amino_start = downstream_hit->dcl->ihmm;
  down_amino_end   = downstream_hit->dcl->jhmm;

  /* Is the downstream node close enough to the upstream edge in hmm coords to be the next exon  */
  if (up_amino_end + MAX_AMINO_EXT < down_amino_start)
    return NULL;

  /* Is the downstream node close enough to the upstream edge in seq coords to be the next exon  */
  intron_len = abs(downstream_hit->dcl->iali - upstream_hit->dcl->jali) + 1;
  if (intron_len < 1 || intron_len > MAX_INTRON_LEN)
    return NULL;

  /* If the hmm coodrs completly overlap the shorter hit, disregard */
  min_hit_len = ESL_MIN(up_amino_end - up_amino_start + 1, down_amino_end - down_amino_start + 1);
  if(up_amino_end - down_amino_start + 1 > min_hit_len*0.75)
    return NULL;

  /* Ensure overlap region is at lesat MIN_AMINO_OVERLAP long plus any gap extention positions in the hits */
  overlap = up_amino_end - down_amino_start + 1;
  if(overlap >= MIN_AMINO_OVERLAP)  // hits naturally have sufficent overlap
    num_ext_aminos = 0;
  else if( overlap > 0)             //hits overlap nuterally bu still need to be extended
    num_ext_aminos = MIN_AMINO_OVERLAP - overlap;
  else                              // hits do not overlap - gap extention + MIN_AMINO_OVERLAP needed
    num_ext_aminos = (down_amino_start - upstream_hit->dcl->jhmm)*2 + MIN_AMINO_OVERLAP;

  /* If the nodes have reached this point we will give them an edge*/
  edge = splice_edge_create();

  edge->overlap_amino_start = down_amino_start;
  edge->overlap_amino_end   = up_amino_end;

  /* If the hits do not overlap by at least MIN_AMINO_OVERLAP hmm positions, extend them */
  if (num_ext_aminos > 0) {
    num_ext_aminos = (num_ext_aminos+1)/2;
    edge->overlap_amino_start -= num_ext_aminos;
    edge->overlap_amino_end   += num_ext_aminos;
  }

  if(edge->overlap_amino_end > down_amino_end)   edge->overlap_amino_end = down_amino_end;
  if(edge->overlap_amino_start < up_amino_start) edge->overlap_amino_start = up_amino_start;
  
  get_overlap_nuc_coords(edge, upstream_hit->dcl, downstream_hit->dcl, target_seq, revcomp);
  
  /* Addd extra nucleotides for splice sites */
  edge->upstream_nuc_end     = ESL_MIN(edge->upstream_nuc_end + 2, target_seq->n);
  edge->downstream_nuc_start = ESL_MAX(edge->downstream_nuc_start - 2, 1);
  
  edge->target_seq_start = target_seq->start;
  edge->target_seq_end   = target_seq->end;
  edge->target_seq_n     = target_seq->n; 

  
  if ((status = find_optimal_splice_site (edge, upstream_hit->dcl, downstream_hit->dcl, gm, hmm, bg, gcode, target_seq, FALSE)) != eslOK) goto ERROR;

  if(edge->splice_score == -eslINFINITY) {
     free(edge);
    
     return NULL;
 
  }

  if(revcomp) {
    edge->upstream_spliced_nuc_end     = target_seq->n - edge->upstream_spliced_nuc_end     + target_seq->end;
    edge->downstream_spliced_nuc_start = target_seq->n - edge->downstream_spliced_nuc_start + target_seq->end;
  }
  else {
    edge->upstream_spliced_nuc_end     = edge->upstream_spliced_nuc_end     + target_seq->start - 1;
    edge->downstream_spliced_nuc_start = edge->downstream_spliced_nuc_start + target_seq->start - 1;
  }

  
  return edge;

  ERROR:
    if(edge != NULL) free(edge);
    return NULL;
}






/*  Function: get_overlap_nuc_coords 
 *
 *  Synopsis: Get nucleotide and trace coords for an overlap region
 *
 *  Purpose:  For two hits that are splice compaitable, get the range 
 *            of nucleotide and P7_TRACE indicies for both the upstream 
 *            and downstream hit that correspond to the range of an 
 *            amino position overlap.
 */
void 
get_overlap_nuc_coords (SPLICE_EDGE *edge, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const ESL_SQ *target_seq, int revcomp)
{

  int       z1,z2;
  int       strand;
  int       extention_nuc_count;
  int       curr_hmm_pos;
  P7_TRACE *up_trace;
  P7_TRACE *down_trace;
  strand = (revcomp? -1 : 1);

  up_trace = upstream->tr;
  down_trace = downstream->tr;

 /***************UPSTREAM*******************/

 /* Get number of nucleotides that need to be added to the end of the
  * upstream edge nucleotide range if the end of the amino range
  * had to be extended to meet the mimimum overlap */
  extention_nuc_count = (edge->overlap_amino_end - upstream->jhmm) * strand * 3;
  edge->upstream_nuc_end = upstream->jali + extention_nuc_count;
  
  /* Set edge->upstream_nucl_start to the begining of the upstream
   * nulceotides that correpond to overlapping aminos */
  edge->upstream_nuc_start = upstream->jali + strand;
 
  for (z2 = up_trace->N-1 ; z2 >= 0; z2--) if (up_trace->st[z2] == p7T_M) break;

  curr_hmm_pos = upstream->jhmm;
  edge->upstream_trace_end = z2;
  
  while(curr_hmm_pos >= edge->overlap_amino_start) {
	
    if      (up_trace->st[z2] == p7T_M) {
      edge->upstream_nuc_start -= up_trace->c[z2] * strand;
      if(up_trace->c[z2] != 3) edge->frameshift = TRUE;
      curr_hmm_pos--;
    }
    else if (up_trace->st[z2] == p7T_I)
      edge->upstream_nuc_start -= 3 * strand;
    else if (up_trace->st[z2] == p7T_D)
      curr_hmm_pos--;
    else
      break;
    z2--;

  }
  edge->upstream_trace_start = z2+1;

  if(revcomp) {
    edge->upstream_nuc_start = target_seq->start - edge->upstream_nuc_start + 1;
    edge->upstream_nuc_end   = target_seq->start - edge->upstream_nuc_end   + 1;
  }
  else {
    edge->upstream_nuc_start = edge->upstream_nuc_start - target_seq->start + 1;
    edge->upstream_nuc_end   = edge->upstream_nuc_end   - target_seq->start + 1;
  }

/***************DOWNSTREAM*******************/

 /* Get number of nucleotides that need to be added to the start of the
  * downstream edge nucleotide range if the end of the amino range
  * had to be extended to meet the mimimum overlap */
  extention_nuc_count = (downstream->ihmm - edge->overlap_amino_start) * strand * 3;
  edge->downstream_nuc_start = downstream->iali - extention_nuc_count;


 /* Set edge->downstream_nucl_end to the end of the downstream
  * nulceotides that correpond to overlapping aminos */
  edge->downstream_nuc_end   = downstream->iali - strand;
    
  for (z1 = 0; z1 < down_trace->N; z1++) if (down_trace->st[z1] == p7T_M) break;
  
  curr_hmm_pos = downstream->ihmm;
  edge->downstream_trace_start = z1;
  while(curr_hmm_pos <= edge->overlap_amino_end) {
    if (down_trace->st[z1] == p7T_M) {
      edge->downstream_nuc_end += down_trace->c[z1] * strand;
      if(down_trace->c[z1] != 3) edge->frameshift = TRUE; 
      curr_hmm_pos++;
    }
    else if (down_trace->st[z1] == p7T_I)
      edge->downstream_nuc_end += 3 * strand;
    else if (down_trace->st[z1] == p7T_D) 
      curr_hmm_pos++;
	else
      break;
    z1++;
  }
  edge->downstream_trace_end = z1;

  if(revcomp) {
    edge->downstream_nuc_start = target_seq->start - edge->downstream_nuc_start + 1;
    edge->downstream_nuc_end   = target_seq->start - edge->downstream_nuc_end   + 1;
  }
  else {
    edge->downstream_nuc_start = edge->downstream_nuc_start - target_seq->start + 1;
    edge->downstream_nuc_end   = edge->downstream_nuc_end   - target_seq->start + 1;
  }

 return;
}

/*  Function: get_split_nuc_coords 
 *
 *  Synopsis: Get nucleotide and trace coords for an overlap region
 *
 *  Purpose:  For two hits that are splice compaitable, get the range 
 *            of nucleotide and P7_TRACE indicies for both the upstream 
 *            and downstream hit that correspond to the range of an 
 *            amino position overlap.
 */
void 
get_split_nuc_coords (SPLICE_EDGE *edge, const P7_DOMAIN *split_hit, int revcomp)
{

  int       z1,z2;
  int       strand;
  int       curr_hmm_pos;
  P7_TRACE *split_trace;
  strand = (revcomp? -1 : 1);

  split_trace = split_hit->tr;

  edge->upstream_nuc_start = split_hit->iali;

  for (z1 = 0; z1 < split_trace->N; z1++)  if (split_trace->st[z1] == p7T_M) break; 

  curr_hmm_pos = split_hit->ihmm;   

  while(curr_hmm_pos < edge->overlap_amino_start) {
    if (split_trace->st[z1] == p7T_M) {
      edge->upstream_nuc_start += split_trace->c[z1] * strand;
      curr_hmm_pos++;
    }
    else if (split_trace->st[z1] == p7T_I)
      edge->upstream_nuc_start += 3 * strand;
    else if (split_trace->st[z1] == p7T_D)
      curr_hmm_pos++;
    else
      break;
    z1++;
  }
  edge->downstream_trace_start = z1;
  for(z2 = z1-1; z2 < split_trace->N; z2++) {
     if(split_trace->k[z2] > edge->overlap_amino_end) break; 
     if(split_trace->st[z2] == p7T_M && split_trace->c[z2] != 3) 
       edge->frameshift = TRUE;
  } 

  edge->upstream_nuc_end = split_hit->jali;

  for (z2 = split_trace->N-1 ; z2 >= 0; z2--) if (split_trace->st[z2] == p7T_M) break; 
  
  curr_hmm_pos = split_hit->jhmm;

  while(curr_hmm_pos > edge->overlap_amino_end) {

    if      (split_trace->st[z2] == p7T_M) {
      edge->upstream_nuc_end -= split_trace->c[z2] * strand;
      curr_hmm_pos--;
    }
    else if (split_trace->st[z2] == p7T_I)
      edge->upstream_nuc_end -= 3 * strand;
    else if (split_trace->st[z2] == p7T_D)
      curr_hmm_pos--;
    else
      break;
    z2--;

  } 
  edge->upstream_trace_end = z2;
  if(revcomp) {
    edge->upstream_nuc_start = edge->target_seq_start - edge->upstream_nuc_start + 1;
    edge->upstream_nuc_end   = edge->target_seq_start - edge->upstream_nuc_end   + 1;
  }
  else {
    edge->upstream_nuc_start = edge->upstream_nuc_start - edge->target_seq_start + 1;
    edge->upstream_nuc_end   = edge->upstream_nuc_end   - edge->target_seq_start + 1;
  }

  edge->downstream_nuc_start = edge->upstream_nuc_start;
  edge->downstream_nuc_end   = edge->upstream_nuc_end;
  
 return;
}



/*  Function: find_optimal_splice_site
 *
 *  Synopsis: Find the best splice site ffor two hits
 *
 *  Purpose:  Given an overlap region for two splice compatible hits find the best
 *            scoring slice site and store all relvant information to theat splice 
 *            site int the <edge> object.
 *
 *  Returns:  <eslOK> on success.
 *
 *  Throws:   <eslEMEM> on allocation failure.
 */
int
find_optimal_splice_site (SPLICE_EDGE *edge, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, const ESL_GENCODE *gcode, const ESL_SQ *target_seq, int split)
{

  int            n,z;
  int            unp,dnp;
  float          lost_sc;
  float         *signal_scores;
  float         *tp_0;
  float         *tp_M;
  P7_HMM        *sub_hmm;
  P7_PROFILE    *sub_model;
  P7_FS_PROFILE *sub_fs_model;
  int            status;

  /* Get the summed ali score for the overlap region covered by the existing alignments */
  lost_sc = 0.;
  z = edge->upstream_trace_end;
  n = z - (upstream->tr->tfrom[0]+1);
  while(upstream->tr->k[z] >= edge->overlap_amino_start) {
    lost_sc += upstream->scores_per_pos[n];
    n--;
    z--;
  }

  /*split hits only have one hit to lose score from */
  if(!split ) {
    z = edge->downstream_trace_start;
    n = 0;
  
    while(downstream->tr->k[z] && downstream->tr->k[z] <= edge->overlap_amino_end) {
      lost_sc += downstream->scores_per_pos[n];
      n++;
      z++;
    }
  }

  /*Get a submodel that covers the overlap region */
  sub_hmm    = extract_sub_hmm(hmm, edge->overlap_amino_start, edge->overlap_amino_end);
  sub_model = NULL;
  sub_fs_model = NULL;

  if(edge->frameshift) {
    sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
    ESL_REALLOC(sub_fs_model->tsc,  sizeof(float)   * (sub_hmm->M+1) * p7P_NTRANS);
    p7_ProfileConfig_fs(sub_hmm, bg, gcode, sub_fs_model, 100, p7_UNIGLOBAL);
    tp_0 = sub_fs_model->tsc;
    tp_M = sub_fs_model->tsc + sub_fs_model->M * p7P_NTRANS;

    /* Add extra transtion position for insetions at begining and end of alignments*/
    tp_0[p7P_MM] = p7P_TSC(gm, edge->overlap_amino_start-1, p7P_MM);
    tp_0[p7P_MI] = p7P_TSC(gm, edge->overlap_amino_start-1, p7P_MI);
    tp_0[p7P_MD] = p7P_TSC(gm, edge->overlap_amino_start-1, p7P_MD);
    tp_0[p7P_IM] = p7P_TSC(gm, edge->overlap_amino_start-1, p7P_IM);
    tp_0[p7P_II] = p7P_TSC(gm, edge->overlap_amino_start-1, p7P_II);
    tp_0[p7P_DM] = p7P_TSC(gm, edge->overlap_amino_start-1, p7P_DM);
    tp_0[p7P_DD] = p7P_TSC(gm, edge->overlap_amino_start-1, p7P_DD);
  
  
    tp_M[p7P_MM] = p7P_TSC(gm, edge->overlap_amino_end, p7P_MM);
    tp_M[p7P_MI] = p7P_TSC(gm, edge->overlap_amino_end, p7P_MI);
    tp_M[p7P_MD] = p7P_TSC(gm, edge->overlap_amino_end, p7P_MD);
    tp_M[p7P_IM] = p7P_TSC(gm, edge->overlap_amino_end, p7P_IM);
    tp_M[p7P_II] = p7P_TSC(gm, edge->overlap_amino_end, p7P_II);
    tp_M[p7P_DM] = p7P_TSC(gm, edge->overlap_amino_end, p7P_DM);
    tp_M[p7P_DD] = p7P_TSC(gm, edge->overlap_amino_end, p7P_DD);
  }
  
  sub_model  = p7_profile_Create (sub_hmm->M, sub_hmm->abc);
  ESL_REALLOC(sub_model->tsc,  sizeof(float)   * (sub_hmm->M+1) * p7P_NTRANS);
  p7_ProfileConfig(sub_hmm, bg, sub_model, 100, p7_UNIGLOBAL);
  tp_0 = sub_model->tsc;
  tp_M = sub_model->tsc + sub_model->M * p7P_NTRANS; 
  

  /* Add extra transtion position for insetions at begining and end of alignments*/
  tp_0[p7P_MM] = p7P_TSC(gm, edge->overlap_amino_start-1, p7P_MM);
  tp_0[p7P_MI] = p7P_TSC(gm, edge->overlap_amino_start-1, p7P_MI);
  tp_0[p7P_MD] = p7P_TSC(gm, edge->overlap_amino_start-1, p7P_MD);
  tp_0[p7P_IM] = p7P_TSC(gm, edge->overlap_amino_start-1, p7P_IM);
  tp_0[p7P_II] = p7P_TSC(gm, edge->overlap_amino_start-1, p7P_II);
  tp_0[p7P_DM] = p7P_TSC(gm, edge->overlap_amino_start-1, p7P_DM);
  tp_0[p7P_DD] = p7P_TSC(gm, edge->overlap_amino_start-1, p7P_DD);


  tp_M[p7P_MM] = p7P_TSC(gm, edge->overlap_amino_end, p7P_MM);
  tp_M[p7P_MI] = p7P_TSC(gm, edge->overlap_amino_end, p7P_MI);
  tp_M[p7P_MD] = p7P_TSC(gm, edge->overlap_amino_end, p7P_MD);
  tp_M[p7P_IM] = p7P_TSC(gm, edge->overlap_amino_end, p7P_IM);
  tp_M[p7P_II] = p7P_TSC(gm, edge->overlap_amino_end, p7P_II);
  tp_M[p7P_DM] = p7P_TSC(gm, edge->overlap_amino_end, p7P_DM);
  tp_M[p7P_DD] = p7P_TSC(gm, edge->overlap_amino_end, p7P_DD);
  
  signal_scores = NULL;
  ESL_ALLOC(signal_scores, sizeof(float) * p7S_SPLICE_SIGNALS);
  p7_splice_SignalScores(signal_scores);
 
  /* Scan the overlap nucleotides for splice signals */
  for(unp = edge->upstream_nuc_start; unp < edge->upstream_nuc_end; unp++) {
    /*GT-AG*/
    if(target_seq->dsq[unp] == 2 && target_seq->dsq[unp+1] == 3) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {
        if(target_seq->dsq[dnp-1] == 0 && target_seq->dsq[dnp] == 2) {
          if (dnp - unp > MIN_INTRON_LEN)
            select_splice_option(edge, gm, sub_model, sub_fs_model, gcode, target_seq, signal_scores[p7S_GTAG], unp-1, dnp+1);
        }
      }
    }
    /*GC-AG*/
    if(target_seq->dsq[unp] == 2 && target_seq->dsq[unp+1] == 1) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {
        if(target_seq->dsq[dnp-1] == 0 && target_seq->dsq[dnp] == 2) {
          if (dnp - unp > MIN_INTRON_LEN)
            select_splice_option(edge, gm, sub_model, sub_fs_model, gcode, target_seq, signal_scores[p7S_GCAG], unp-1, dnp+1);
        }
      }
    }
    /*AT-AC*/
    if(target_seq->dsq[unp] == 0 && target_seq->dsq[unp+1] == 3) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {  
        if(target_seq->dsq[dnp-1] == 0 && target_seq->dsq[dnp] == 1) {
          if (dnp - unp > MIN_INTRON_LEN)
            select_splice_option(edge, gm, sub_model, sub_fs_model, gcode, target_seq, signal_scores[p7S_ATAC], unp-1, dnp+1);
        }
      }
    }
  }
  
  if(edge->splice_score != -eslINFINITY)
    edge->splice_score -= lost_sc;
  
  if(signal_scores != NULL) free(signal_scores);
  if(sub_hmm       != NULL) p7_hmm_Destroy(sub_hmm);
  if(sub_model     != NULL) p7_profile_Destroy(sub_model); 
  if(sub_fs_model  != NULL) p7_profile_fs_Destroy(sub_fs_model);

  return eslOK;

  ERROR:
    if(signal_scores != NULL) free(signal_scores);
    if(sub_hmm       != NULL) p7_hmm_Destroy(sub_hmm);
    if(sub_model     != NULL) p7_profile_Destroy(sub_model);
    if(sub_fs_model  != NULL) p7_profile_fs_Destroy(sub_fs_model);
    return status;

}






/*  Function: select_splice_option
 *
 *  Synopsis: Find the best splice sitefor two hits at a particular model postion
 *
 *  Purpose:  Given a model posions in the overlap for two splice compatible hits 
 *            and the surronuding upstream and downstream nucleotides find the best
 *            scoring slice site.
 *
 *  Returns:  <eslOK> on success.
 *
 *  Throws:   <eslEMEM> on allocation failure.
 */
int 
select_splice_option (SPLICE_EDGE *edge, const P7_PROFILE *gm, P7_PROFILE *sub_model, P7_FS_PROFILE *sub_fs_model, const ESL_GENCODE *gcode, const ESL_SQ *target_seq, float signal_score, int up_nuc_pos, int down_nuc_pos)
{

  int         i,z;
  int         nuc_seq_idx;
  int         nuc_seq_len;
  int         amino_len;
  int         amino;
  int         upstream_nuc_cnt;
  int         upstream_amino_cnt;
  int         upstream_amino_end;
  float       sum_ali_sc;
  float       vitsc;
  float       overlap_sc;
  ESL_DSQ    *nuc_dsq;
  ESL_DSQ    *amino_dsq;
  P7_GMX     *vit_mx;
  P7_TRACE   *tr;
  int         status;

  nuc_dsq   = NULL;
  amino_dsq = NULL;
  vit_mx    = NULL;
  tr        = NULL;
 
  /*Get the overlap nucleotides that correspond to the current splice signal*/
  nuc_seq_len  = up_nuc_pos - edge->upstream_nuc_start + 1;
  nuc_seq_len += edge->downstream_nuc_end - down_nuc_pos + 1;

  /* If the splice signal is right at the overlap boundries all amino acid positions are deletions */
  if(nuc_seq_len == 0) {
    sum_ali_sc = 0; 
    for(i = edge->overlap_amino_start; i < edge->overlap_amino_end; i++)
      sum_ali_sc += gm->tsc[i * p7P_NTRANS + p7H_DD];

    overlap_sc = sum_ali_sc + signal_score;
    
    if (overlap_sc > edge->splice_score) {
      edge->signal_score = signal_score;
      edge->splice_score = overlap_sc;
      edge->upstream_spliced_nuc_end = up_nuc_pos;
      edge->downstream_spliced_nuc_start = down_nuc_pos;
      edge->upstream_spliced_amino_end = edge->overlap_amino_start - 1;
      edge->downstream_spliced_amino_start = edge->overlap_amino_start;
    }
    return eslOK;
  }

  
  /* When not using the frameshift algorithm we can only 
   * allow mod three length nucelotide sequences  */
  
  if((!edge->frameshift) && nuc_seq_len % 3) return eslOK;

  ESL_ALLOC(nuc_dsq,   sizeof(ESL_DSQ) * (nuc_seq_len+2));
  nuc_seq_idx = 0;
  nuc_dsq[nuc_seq_idx] = eslDSQ_SENTINEL;
  nuc_seq_idx++; 
  for(i = edge->upstream_nuc_start; i <= up_nuc_pos; i++) {
    nuc_dsq[nuc_seq_idx] = target_seq->dsq[i];
    nuc_seq_idx++;
  }
  upstream_nuc_cnt = nuc_seq_idx-1;
  upstream_amino_cnt = (nuc_seq_idx-1) / 3;
  if((nuc_seq_idx-1) % 3)  upstream_amino_cnt++;

  for (i = down_nuc_pos; i <= edge->downstream_nuc_end; i++) {
    nuc_dsq[nuc_seq_idx] = target_seq->dsq[i];
    nuc_seq_idx++;
  }
  nuc_dsq[nuc_seq_idx] = eslDSQ_SENTINEL;

  if(nuc_seq_len % 3) {
     vit_mx = p7_gmx_fs_Create(sub_fs_model->M, nuc_seq_len, nuc_seq_len, p7P_CODONS); 
     p7_fs_ReconfigLength(sub_fs_model, nuc_seq_len);
     p7_fs_global_Viterbi(nuc_dsq, gcode, nuc_seq_len, sub_fs_model, vit_mx, &vitsc);
  }
  else {
    /* Translate overalp nucleotides to amino sequence */
    amino_len = nuc_seq_len / 3;
    ESL_ALLOC(amino_dsq, sizeof(ESL_DSQ) * (amino_len+2));

    amino_dsq[0] = eslDSQ_SENTINEL;
    nuc_seq_idx = 1;
    for(i = 1; i <= amino_len; i++) {
      amino = esl_gencode_GetTranslation(gcode,&nuc_dsq[nuc_seq_idx]);
      amino_dsq[i] = amino;
      nuc_seq_idx+=3;
    }
    amino_dsq[amino_len+1] = eslDSQ_SENTINEL;

    /* Align translated overlap amino acids to submodel */
    vit_mx = p7_gmx_Create(sub_model->M, amino_len);
    p7_ReconfigLength(sub_model, amino_len); 
    p7_global_Viterbi(amino_dsq, amino_len, sub_model, vit_mx, &vitsc);
  }
  
    /* vitsc from global alignment is the same as ali_score */ 
  if (vitsc != -eslINFINITY) {
    overlap_sc = vitsc + signal_score;
    
    if (overlap_sc > edge->splice_score) {
      /*use trace to find splice point*/
     
      if(nuc_seq_len % 3) {
        tr = p7_trace_fs_Create();
        p7_fs_global_Trace(nuc_dsq, nuc_seq_len, sub_fs_model, vit_mx, tr);

        z = 0;
        while(z < tr->N) {
          if(tr->i[z] >= upstream_nuc_cnt) {
            upstream_amino_end = tr->k[z] + edge->overlap_amino_start - 1;
            z = tr->N;
          } 
          z++;
        }
      }
      else {
        tr = p7_trace_Create();
        p7_global_Trace(amino_dsq, amino_len, sub_model, vit_mx, tr); 
       
        z = 0;
        while(z < tr->N) {
          if(tr->i[z] == upstream_amino_cnt) {
            upstream_amino_end = tr->k[z] + edge->overlap_amino_start - 1;
            z = tr->N;
          }
          z++;
        }
      }
      edge->signal_score = signal_score;
      edge->splice_score = overlap_sc;
      edge->upstream_spliced_nuc_end = up_nuc_pos;
      edge->downstream_spliced_nuc_start = down_nuc_pos;
      edge->upstream_spliced_amino_end = upstream_amino_end;
      edge->downstream_spliced_amino_start = upstream_amino_end + 1;
    }
  }
  /* Case where all translated aminos are stop codons - treat them as insertions and all hmm positions as deletions */
  else {
    sum_ali_sc = gm->tsc[(edge->overlap_amino_start-1) * p7P_NTRANS + p7H_MI];
    for(i = 1; i < amino_len; i++)
      sum_ali_sc += gm->tsc[(edge->overlap_amino_start-1) * p7P_NTRANS + p7H_II];
    sum_ali_sc += gm->tsc[(edge->overlap_amino_start-1) * p7P_NTRANS + p7H_DD];
    for(i = edge->overlap_amino_start; i < edge->overlap_amino_end; i++)
      sum_ali_sc += gm->tsc[i * p7P_NTRANS + p7H_DD];

    overlap_sc = sum_ali_sc + signal_score;
    
    if (overlap_sc > edge->splice_score) {
      edge->signal_score = signal_score;
      edge->splice_score = overlap_sc;
      edge->upstream_spliced_nuc_end = up_nuc_pos;
      edge->downstream_spliced_nuc_start = down_nuc_pos;
      edge->upstream_spliced_amino_end = edge->overlap_amino_start - 1;
      edge->downstream_spliced_amino_start = edge->overlap_amino_start;
    }
  }

  if(nuc_dsq   != NULL) free(nuc_dsq);
  if(amino_dsq != NULL) free(amino_dsq);
  if(vit_mx    != NULL) p7_gmx_Destroy(vit_mx);
  if(nuc_seq_len % 3) { if(tr != NULL) p7_trace_fs_Destroy(tr); }
  else                 { if(tr != NULL) p7_trace_Destroy(tr); }
  return eslOK;

  ERROR:
    if(nuc_dsq   != NULL) free(nuc_dsq);
    if(amino_dsq != NULL) free(amino_dsq);
    if(vit_mx    != NULL) p7_gmx_Destroy(vit_mx);
    if(nuc_seq_len % 3) { if(tr != NULL) p7_trace_fs_Destroy(tr); }
    else                 { if(tr != NULL) p7_trace_Destroy(tr); } 
    return status;
}




int
add_edge_to_graph(SPLICE_GRAPH *graph, SPLICE_EDGE *edge)
{

  int up_node_id;
  int down_node_id;
  SPLICE_NODE *tmp_node;
  SPLICE_EDGE *tmp_edge;

  up_node_id   = edge->upstream_node_id;
  down_node_id = edge->downstream_node_id;

  graph->out_edge_cnt[up_node_id]++;
  graph->in_edge_cnt[down_node_id]++;

  /* Find the end of the up nodes list of edges and append the new edge */
  tmp_edge = graph->edges[up_node_id];
  if(tmp_edge == NULL)
    graph->edges[up_node_id] = edge;
  else {
    while(tmp_edge->next != NULL)
      tmp_edge = tmp_edge->next;

    tmp_edge->next = edge;
    edge->prev = tmp_edge;
  }

  /* Set the path exists and path checked valuse to TRUE */
  tmp_node = get_node(graph, up_node_id, down_node_id);
  tmp_node->path_exists  = TRUE;
  tmp_node->path_checked = TRUE;

  graph->num_edges++;
  
  return eslOK;
}


int
fill_holes_in_graph(TARGET_RANGE *target_range, SPLICE_GRAPH *graph, const P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, ESL_GENCODE *gcode, const ESL_SQFILE *seq_file)
{
  int h1, h2;
  int up,down;
  int loop_start, loop_end, loop_inc;
  int seq_min, seq_max;
  int down_hit_min, down_hit_max;
  int hmm_gap_len;
  int seq_gap_len;
  int num_hits;
  int num_in_gap;
  float downstream_loss;
  int *in_gap;
  P7_TOPHITS *th;
  P7_HIT     **top_ten;
  P7_HIT     *up_hit;
  P7_HIT     *down_hit;
  SPLICE_NODE *tmp_node;
  SPLICE_GAP *gap;
  ESL_SQ          *target_seq;
  int         status;

  th = target_range->th;
  
  in_gap = NULL;
  ESL_ALLOC(in_gap, target_range->orig_N * sizeof(int));
  
  target_seq = NULL;


  /* varaibles to controll direction of loops depending on strand */
  if(graph->revcomp) {
    loop_start = target_range->orig_N-1;
    loop_end = -1;
    loop_inc = -1;
  }
  else {
    loop_start = 0; 
    loop_end = target_range->orig_N;
    loop_inc = 1;
  }
   
  /* Find the largest region between two hits in a target range that is less than MAX_GAP_RANGE and pull out 
   * a sub sequence for that region. Use that sub sequence to fill any gaps and to find all possible edges. */ 
  up = loop_start;
  while(up != loop_end) {
    
    if(!target_range->in_target_range[up]) { up += loop_inc; continue; }

    //TODO 
    /* only check gaps if the current downstream path score of up (path_score -hit_score) is less than the path scrore of down*/
    for(down = 0; down < target_range->orig_N; down++) in_gap[down] = FALSE;

    in_gap[up] = TRUE;
    num_in_gap = 1;
    
    up_hit = th->hit[up];

    seq_min = ESL_MIN(up_hit->dcl->iali, up_hit->dcl->jali);   
    seq_max = ESL_MAX(up_hit->dcl->iali, up_hit->dcl->jali);

    /*Find potetnital lost score for breaking the current best downstream path of up*/
    downstream_loss = graph->path_scores[up] - up_hit->dcl->aliscore;

    /* find sub sequence coords */
    down = loop_start;
    while (down != loop_end)  {
     
      
      if ( !target_range->in_target_range[down] || down == up) { down += loop_inc; continue;}
      

      /* Detemine id down is downstream of up */
      tmp_node = get_node(graph,up,down);
 
      if(tmp_node != NULL) {

        if( tmp_node->path_exists || path_exists(graph, up, down)) { 
          tmp_node->path_exists  = TRUE; 
          tmp_node->path_checked = TRUE;
          down += loop_inc; 
          continue;
        }

        if(!tmp_node->gap_checked) {
 
          down_hit = th->hit[down];

          if(graph->revcomp)
            seq_gap_len = up_hit->dcl->jali - down_hit->dcl->iali - 1;
          else
            seq_gap_len = down_hit->dcl->iali - up_hit->dcl->jali - 1;

          if(seq_gap_len < MAX_GAP_RANGE) {
            /* make sure connecting to down is worth losing any current downstream path of up */
            
            if(graph->path_scores[down] > downstream_loss) {    
                      
              in_gap[down] = TRUE;
              num_in_gap++;

              down_hit_min = ESL_MIN(down_hit->dcl->iali, down_hit->dcl->jali);
              down_hit_max = ESL_MAX(down_hit->dcl->iali, down_hit->dcl->jali);

              seq_min = ESL_MIN(down_hit_min, seq_min);
              seq_max = ESL_MAX(down_hit_max, seq_max);
            }
          }         
        }
      }    
      down += loop_inc;
    } 
   
    
    if ( num_in_gap <= 1 ) { up += loop_inc; continue; } 
    
    if(target_seq != NULL) esl_sq_Destroy(target_seq);
    target_seq =  get_sub_sequence(seq_file, target_range->seqname, seq_min, seq_max, target_range->revcomp);
    
    /* For each pair of hits in the current sub sequence range check if they have an 
     * edge or if there is a gap between them that could be filled with new hits */
    h1 = loop_start;
    while (h1 != loop_end) {
      
      if ( !in_gap[h1] ) { h1 += loop_inc; continue; }
        
      up_hit = th->hit[h1];
      h2 = loop_start;
      while(h2 != loop_end) {
        if ( !in_gap[h2] || h1 == h2) { h2 += loop_inc; continue; }
        tmp_node = get_node(graph, h1, h2); 
         
        if(tmp_node != NULL) { 
          /* If we have not yet checked for a path then check for one */
          if( !tmp_node->path_checked) {
            if(path_exists(graph, h1, h2)) tmp_node->path_exists = TRUE;
          }

          if( !tmp_node->gap_checked && !tmp_node->path_exists ) {
            tmp_node->gap_checked = TRUE;
            down_hit = th->hit[h2];
        
            num_hits = 0;
            hmm_gap_len = down_hit->dcl->ihmm - up_hit->dcl->jhmm - 1;
            
            if(graph->revcomp)
              seq_gap_len = up_hit->dcl->jali - down_hit->dcl->iali - 1;
            else
              seq_gap_len = down_hit->dcl->iali - up_hit->dcl->jali - 1;
                    
            if(hmm_gap_len > 0 && seq_gap_len >= hmm_gap_len*3) { 
              
              if ((gap = find_the_gap(th, gm, target_seq, target_range->orig_N, h1, h2, graph->revcomp)) == NULL) goto ERROR;
              //if((top_ten = align_the_gap(target_range, gm, hmm, bg, target_seq, gcode, gap, &num_hits, graph->revcomp)) == NULL) goto ERROR;
              
              if((top_ten = align_the_gap2(target_range, gap, hmm, bg, target_seq, gcode, graph->revcomp, &num_hits, TRUE)) == NULL) goto ERROR; 
              if(num_hits)
                if ((status = bridge_the_gap(target_range, graph, top_ten, gm, hmm, bg, gcode, target_seq, h1, h2, num_hits)) != eslOK) goto ERROR;

              if(gap     != NULL) free(gap);
              if(top_ten != NULL) free(top_ten);
            }           
          }
        }
        h2 += loop_inc;        
      } 
      h1 += loop_inc;
    }
    up += loop_inc;
  }

/* 
  for(up = 0; up < graph->num_nodes; up++) {
  
    if(!target_range->in_target_range[up]) continue; 

    for(down = 0; down < graph->num_nodes; down++) {  
      if ( !target_range->in_target_range[down] || down == up) continue;
      if(up >= target_range->orig_N && down >= target_range->orig_N) continue;
 
      if(is_upstream(graph, up, down) && !edge_exists(graph, up. down)) {
        

      }
      
    }
  }
*/

  if ( in_gap     != NULL) free(in_gap); 
  if ( target_seq != NULL) esl_sq_Destroy(target_seq);
  return eslOK;

  ERROR:
    if ( in_gap     != NULL) free(in_gap);
    if ( target_seq != NULL) esl_sq_Destroy(target_seq);
    if ( gap        != NULL) free(gap);
    if ( top_ten    != NULL) free(top_ten);
    return status;
}






int 
reset_edge_counts(TARGET_RANGE *target_range, SPLICE_GRAPH *graph) {

  int up;
  SPLICE_EDGE *tmp_edge;

  for (up = 0; up < graph->num_nodes; up++) {
    graph->out_edge_cnt[up] = 0;
    graph->in_edge_cnt[up] = 0;
  }

  for (up = 0; up < graph->num_nodes; up++) { 
    tmp_edge = graph->edges[up];
    while(tmp_edge != NULL) {
      if(target_range->in_target_range[tmp_edge->downstream_node_id]) { 
        graph->out_edge_cnt[up]++;
        graph->in_edge_cnt[tmp_edge->downstream_node_id]++;
      }
      tmp_edge = tmp_edge->next;
    }
  }

  return eslOK;
}





int
enforce_range_bounds(SPLICE_GRAPH *graph, P7_TOPHITS *th, int64_t* range_bound_mins, int64_t* range_bound_maxs, int range_cnt)
{
  
  int     i;
  int     up, down;
  int64_t up_hit_min, up_hit_max;
  int64_t down_hit_min, down_hit_max;
  P7_HIT *up_hit;
  P7_HIT *down_hit;
  SPLICE_EDGE *tmp_edge;
  SPLICE_EDGE *prev_edge;
  SPLICE_EDGE *next_edge;

  if (range_cnt == 0 ) return eslOK;

  for(up = 0; up < th->N; up++) {
    for(down = 0; down < th->N; down++) {
      if( up == down) continue;

      tmp_edge = get_edge(graph, up, down);
      if(tmp_edge == NULL) continue;
      prev_edge = tmp_edge->prev;      

      up_hit   = th->hit[up];
      down_hit = th->hit[down];

      up_hit_min   = ESL_MIN(up_hit->dcl->iali, up_hit->dcl->jali);
      up_hit_max   = ESL_MAX(up_hit->dcl->iali, up_hit->dcl->jali);
      down_hit_min = ESL_MIN(down_hit->dcl->iali, down_hit->dcl->jali);
      down_hit_max = ESL_MAX(down_hit->dcl->iali, down_hit->dcl->jali);

      for(i = 0; i < range_cnt; i++) {
        
        /* If the upstream hit ends after the range bound begins */
        if(up_hit_max >= range_bound_mins[i]) {
          /* If the upstream hit starts before the range bound ends or
             the downstream hit strts before the range bound ends*/
          if (up_hit_min <= range_bound_maxs[i] || down_hit_min <= range_bound_maxs[i]) {
             
            next_edge = tmp_edge->next;
            if(next_edge != NULL) next_edge->prev  = prev_edge;
            if(prev_edge == NULL) graph->edges[up] = next_edge; 
            else                  prev_edge->next  = next_edge; 
            free(tmp_edge);
            tmp_edge = NULL;
            graph->out_edge_cnt[up]--;
            graph->in_edge_cnt[down]--;
            graph->num_edges--;
            i = range_cnt;
          }
        }
        else if( down_hit_max >= range_bound_mins[i]) {
          
          next_edge = tmp_edge->next;
          if(next_edge != NULL) next_edge->prev  = prev_edge;
          if(prev_edge == NULL) graph->edges[up] = next_edge;
          else                  prev_edge->next  = next_edge;
          free(tmp_edge);
          tmp_edge = NULL;
          graph->out_edge_cnt[up]--;
          graph->in_edge_cnt[down]--;
          graph->num_edges--;
          i = range_cnt; 
        }
      }
    }
  }

  return eslOK;
}



SPLICE_PATH*
evaluate_paths (TARGET_RANGE *target_range, SPLICE_GRAPH *graph)
{

  int          i;
  int          path_len;
  int          start_node;
  int          prev_node;
  int          curr_node;
  int          next_node;
  int          contains_orig;
  int          step_cnt;
  float        best_start_score;
  SPLICE_EDGE *in_edge;
  SPLICE_EDGE *out_edge;
  SPLICE_EDGE *prev_in_edge;
  SPLICE_EDGE *prev_out_edge;
  SPLICE_EDGE *next_in_edge;
  SPLICE_EDGE *next_out_edge;
  SPLICE_PATH *path;
  P7_TOPHITS  *th;
  int         status;

  th = target_range->th;

  contains_orig = FALSE;
 
  /* Find best scoreing paths */ 
  if((status = longest_path_upstream(target_range, graph)) != eslOK) goto ERROR;

  //graph_dump(stdout, target_range, graph, FALSE);
  while(!contains_orig) { 
    /* Find the best place to start our path */ 
    best_start_score = -eslINFINITY;
    start_node  = -1;
    for (i = 0; i < graph->num_nodes; i++) {
 
      if(graph->path_scores[i] > best_start_score) {
        best_start_score = graph->path_scores[i];
   	    start_node  = i;
      }
    } 
    
    if(start_node < 0)  ESL_XEXCEPTION(eslFAIL, "Failed to find path in splice graph");
   
    if(start_node < target_range->orig_N)
      contains_orig = TRUE;
      
    /* Get first edge in path (if any) and ckeck that it is not backward */ 
    curr_node = start_node;
    
    path_len  = 1;

    prev_in_edge = NULL;
    in_edge = NULL;
    prev_out_edge = NULL;
    out_edge = NULL;
    
    while(graph->best_out_edge[curr_node] >= 0) {
      if(curr_node < target_range->orig_N)
        contains_orig = TRUE;
     
      next_node = graph->best_out_edge[curr_node];
      
      out_edge = get_edge(graph, curr_node, next_node);
      if(out_edge == NULL) ESL_XEXCEPTION(eslFAIL, "Edge does not exist");
      prev_out_edge = out_edge->prev;

      /* Check if edges create backward node (node whose upstream splce occures after its downstream splice */
      if( in_edge != NULL) {
        if((in_edge->downstream_spliced_amino_start >= out_edge->upstream_spliced_amino_end) ||
           ((graph->revcomp && out_edge->upstream_spliced_nuc_end >= in_edge->downstream_spliced_nuc_start) ||
           (!graph->revcomp && in_edge->downstream_spliced_nuc_start >= out_edge->upstream_spliced_nuc_end))) {
           
          //Remove whichever edge has the lowest socre
          if(out_edge->splice_score < in_edge->splice_score) {
            next_out_edge = out_edge->next;
            if(next_out_edge != NULL) next_out_edge->prev     = prev_out_edge;
            if(prev_out_edge == NULL) graph->edges[curr_node] = next_out_edge;
			else                      prev_out_edge->next     = next_out_edge;
            free(out_edge);

            graph->out_edge_cnt[curr_node]--;
            graph->in_edge_cnt[next_node]--;
            graph->num_edges--;
          }
          else {
            next_in_edge = in_edge->next;
            if(next_in_edge != NULL) next_in_edge->prev = prev_in_edge;
            if(prev_in_edge == NULL) graph->edges[prev_node] = next_in_edge;
            else                     prev_in_edge->next      = next_in_edge; 
            free(in_edge);
           
            graph->out_edge_cnt[prev_node]--;
            graph->in_edge_cnt[curr_node]--;
            graph->num_edges--;
          }
          path = evaluate_paths (target_range, graph);
          return path;
        }
      }
      prev_node = curr_node;
      curr_node = next_node;
      prev_in_edge = prev_out_edge;
      in_edge = out_edge;
      path_len++;
    }

    if(curr_node < target_range->orig_N)
      contains_orig = TRUE;

    if(!contains_orig) 
      graph->path_scores[start_node] = -eslINFINITY; 
    
  }    
  
  /* Once a viable path has been found build the path */
  if ((path = splice_path_create(path_len)) == NULL) goto ERROR; 

  path->revcomp = graph->revcomp;

  path->hit_scores[0]    = th->hit[start_node]->dcl->aliscore;
  path->edge_scores[0]   = 0.;
  path->signal_scores[0] = 0.;

  path->node_id[0] = start_node;
  path->split[0]   = FALSE;
  path->hits[0]    = th->hit[start_node];
  
  path->downstream_spliced_amino_start[0] = path->hits[0]->dcl->ihmm;
  path->downstream_spliced_nuc_start[0]   = path->hits[0]->dcl->iali;
  
  curr_node = start_node;
  step_cnt = 1;
  while (step_cnt < path_len) {
    next_node = graph->best_out_edge[curr_node]; 
    out_edge = get_edge(graph, curr_node, next_node);

    path->hit_scores[step_cnt]    = th->hit[next_node]->dcl->aliscore;
    path->edge_scores[step_cnt]   = out_edge->splice_score;
    path->signal_scores[step_cnt] = out_edge->signal_score;
 
    path->node_id[step_cnt] = next_node;
    path->split[step_cnt]   = FALSE;
    path->hits[step_cnt]    = th->hit[next_node];
 
    path->upstream_spliced_amino_end[step_cnt]     = out_edge->upstream_spliced_amino_end;
    path->downstream_spliced_amino_start[step_cnt] = out_edge->downstream_spliced_amino_start;
    path->upstream_spliced_nuc_end[step_cnt]       = out_edge->upstream_spliced_nuc_end;
    path->downstream_spliced_nuc_start[step_cnt]   = out_edge->downstream_spliced_nuc_start;
 
    curr_node = next_node;
    step_cnt++;
  }

  path->upstream_spliced_amino_end[step_cnt] = path->hits[step_cnt-1]->dcl->jhmm;  
  path->upstream_spliced_nuc_end[step_cnt]   = path->hits[step_cnt-1]->dcl->jali;
  
  return path; 
 
  ERROR:
    if(path != NULL) splice_path_destroy(path);
    return NULL;
}


int
longest_path_upstream (TARGET_RANGE *target_range, SPLICE_GRAPH *graph)
{

  int   i;
  int up, down;
  int   stack_size;
  float step_score;
  int   *visited;
  int   *stack;
  SPLICE_EDGE *edge;
  SPLICE_EDGE *tmp_edge;
  int    status;

  for(i = 0; i < graph->num_nodes; i++) {
    if(target_range->in_target_range[i]) 
      graph->path_scores[i]   = graph->hit_scores[i];
    else
      graph->path_scores[i]   = -eslINFINITY;
    graph->best_out_edge[i] = -1;
  }

  /* Append source node downstream of all nodes with no doutgoing edges*/
  if((status = splice_graph_grow(graph)) != eslOK) goto ERROR;
  graph->hit_scores[graph->num_nodes]  = 0.;
  graph->edges[graph->num_nodes] = NULL;
  for(up = 0; up < graph->num_nodes; up++) {

	if(!target_range->in_target_range[up]) continue;

    if(graph->out_edge_cnt[up]) continue; 	
    
    edge = splice_edge_create();
    edge->upstream_node_id = up;
    edge->downstream_node_id = graph->num_nodes;
    edge->splice_score = 0.;
      
    graph->edges[up] = edge;
    edge = NULL;
  } 
  graph->path_scores[graph->num_nodes] = 0.;
  graph->num_nodes++;


  ESL_ALLOC(visited, sizeof(int) * graph->num_nodes);
  esl_vec_ISet(visited,   graph->num_nodes, 0); 
 
  ESL_ALLOC(stack, sizeof(int) * graph->num_nodes);

  stack_size = 0;
  
  for(i = 0; i < graph->num_nodes; i++) {
	if (i < graph->num_nodes-1 && !target_range->in_target_range[i]) continue;
	if(!visited[i]) {
	  topological_sort_upstream(target_range, graph, visited, stack, &stack_size, i);
    }
  } 

   while(stack_size > 0) {
    /*pop top of stack */
	down = stack[stack_size-1];
	stack_size--; 
    if(down == graph->num_nodes-1) continue;
   
	/* Find nodes with ougoing edge to down*/
    for(up = 0; up < graph->num_nodes-1; up++) {
      if (!target_range->in_target_range[up]) continue;
  
      tmp_edge = get_edge(graph, up, down); 
      if(tmp_edge != NULL) {
        step_score = graph->hit_scores[up] + tmp_edge->splice_score + graph->path_scores[down]; 
        if(graph->path_scores[up] <= step_score) {
          graph->path_scores[up]   = step_score;
          graph->best_out_edge[up] = down;
        }        
      }
	}
  }  

  /* Erase source node */
  graph->num_nodes--;

  for(up = 0; up < graph->num_nodes; up++) {
    if(!target_range->in_target_range[up]) continue;
    if(graph->out_edge_cnt[up]) continue;
 
    if(graph->edges[up] != NULL) {
     
      edge = graph->edges[up];
      free(edge);
      graph->edges[up] = NULL;   
    }
  }

  if(visited  != NULL) free(visited);
  if(stack    != NULL) free(stack);
  return eslOK;
 
  ERROR:
	if(visited  != NULL) free(visited);
	if(stack    != NULL) free(stack);
	return status;    
}




int
topological_sort_upstream(TARGET_RANGE *target_range, SPLICE_GRAPH *graph, int *visited, int *stack, int *stack_size, int node)
{

  int i;

  visited[node] = TRUE;

  for(i = 0; i < graph->num_nodes; i++) {
    if(i < graph->num_nodes-1 && !target_range->in_target_range[i]) continue; 
    
    if(visited[i]) continue;

    if(edge_exists(graph, i, node)) 
      topological_sort_upstream(target_range, graph, visited, stack, stack_size, i);
    
  }
  
  stack[*stack_size] = node;
  *stack_size += 1;
  
  return eslOK;

}



int
split_hits_in_path2 (TARGET_RANGE *target_range, SPLICE_GRAPH *graph, SPLICE_PATH *path, const P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, const ESL_GENCODE *gcode, ESL_SQ *path_seq)
{

  int i, h, p, z;
  int z1, z2;
  int I_cnt;
  int low_pp_cnt;
  int splitable;
  int rev;
  int num_hits;
  P7_HIT *curr_hit;
  P7_HIT **split_hits;
  P7_TRACE *tr;
  SPLICE_GAP *gap;
  SPLICE_EDGE *first_edge;
  SPLICE_EDGE *last_edge;
  SPLICE_EDGE *prev_edge;
  SPLICE_EDGE *curr_edge;
  SPLICE_EDGE **split_edges;
  int     status;

  first_edge  = NULL;
  last_edge   = NULL;
  prev_edge   = NULL;
  curr_edge   = NULL;
  split_edges = NULL;

  gap = NULL;
  ESL_ALLOC(gap, sizeof(SPLICE_GAP));

  rev = 1;
  if(path->revcomp) rev = -1;

  p = 0;
  while(p < path->path_len) {
    if(path->node_id[p] >= target_range->orig_N) { p++; continue; }
   
    curr_hit = path->hits[p];
  
    tr = curr_hit->dcl->tr;

    I_cnt = 0;
    low_pp_cnt = 0;
    splitable = FALSE;
	/* find the position in trace that coorespond to the path splice boudries*/
    for (z1 = 0;       z1 < tr->N; z1++) if ( tr->st[z1] == p7T_M ) break; 
    
    while (tr->k[z1] < path->downstream_spliced_amino_start[p]) z1++;
    
    for (z2 = tr->N-1; z2 >= 0;    z2--) { if ( tr->st[z2] == p7T_M && tr->k[z2] <= path->upstream_spliced_amino_end[p+1])   break; }    

	z = z1; 
    while (z < z2) {
      if (tr->st[z] == p7T_I) {
        I_cnt++;
      }
      else if(tr->st[z] == p7T_M && tr->pp[z] < 0.45) {
        low_pp_cnt++;
      }
      else if (I_cnt + low_pp_cnt > 10) {
        splitable = TRUE;
        break;
      } 
      else {
        I_cnt = 0;
        low_pp_cnt = 0;
      }
      z++;
    }
    
    if(splitable) {
       
      gap->hmm_start = path->downstream_spliced_amino_start[p];
      gap->hmm_end   = path->upstream_spliced_amino_end[p+1];

      if(p != 0)  gap->seq_start = path->upstream_spliced_nuc_end[p] + (MIN_INTRON_LEN * rev);
      else        gap->seq_start = path_seq->start; 
 
      if(p != path->path_len -1) gap->seq_end = path->downstream_spliced_nuc_start[p+1] - (MIN_INTRON_LEN * rev);
      else                       gap->seq_end = path_seq->end; 
      
     // gap->seq_start = path->downstream_spliced_nuc_start[p];
     // gap->seq_end   = path->upstream_spliced_nuc_end[p+1];

      split_hits = align_the_gap2(target_range, gap, hmm, bg, path_seq, gcode, path->revcomp, &num_hits, FALSE); 
      //if(p == 4) { p7_alidisplay_Dump(stdout, split_hits[num_hits-1]->dcl->ad);  p7_alidisplay_Dump(stdout, path->hits[p+1]->dcl->ad); }
               
      if(num_hits <= 1) { 
   
        p7_alidisplay_Destroy(split_hits[0]->dcl->ad);
        p7_trace_fs_Destroy(split_hits[0]->dcl->tr);
        free(split_hits[0]->dcl->scores_per_pos);
        p7_hit_Destroy(split_hits[0]);
         
        free(split_hits);
        split_hits = NULL;

        p++; 
        continue; 
 
      }

      /*Ensure set of spliced hits is full length - i.e. can be spliced 
       * to the original hits upstream and downstream path neighbors */
      first_edge = last_edge = NULL;
      if (p!=0)                first_edge = connect_nodes_with_edges(path->hits[p-1],        split_hits[0],   gm, hmm, bg, gcode, path_seq, path->revcomp);  
       
      if (p!=path->path_len-1) last_edge  = connect_nodes_with_edges(split_hits[num_hits-1], path->hits[p+1], gm, hmm, bg, gcode, path_seq, path->revcomp);
        
      /* If we failed to find an edge to the upstream and downstream hits use the original edges positions*/ 
      if(p!=0 && first_edge == NULL) {
         first_edge = splice_edge_create();
         first_edge->upstream_spliced_amino_end     = path->upstream_spliced_amino_end[p];
         first_edge->downstream_spliced_amino_start = path->downstream_spliced_amino_start[p];
         first_edge->upstream_spliced_nuc_end       = path->upstream_spliced_nuc_end[p];
         first_edge->downstream_spliced_nuc_start   = path->downstream_spliced_nuc_start[p];

         first_edge->splice_score = path->edge_scores[p-1]; 
         first_edge->signal_score = path->signal_scores[p-1];         
      }

      if(p!=path->path_len-1 && last_edge == NULL) { 
         last_edge = splice_edge_create();
         last_edge->upstream_spliced_amino_end     = path->upstream_spliced_amino_end[p+1];
         last_edge->downstream_spliced_amino_start = path->downstream_spliced_amino_start[p+1];
         last_edge->upstream_spliced_nuc_end       = path->upstream_spliced_nuc_end[p+1];
         last_edge->downstream_spliced_nuc_start   = path->downstream_spliced_nuc_start[p+1];

         last_edge->splice_score = path->edge_scores[p];
         last_edge->signal_score = path->signal_scores[p];      
      }
       
      /* if first or last edge is backward - abandon ship */
      if(first_edge != NULL) {
        if(first_edge->upstream_spliced_amino_end <= path->downstream_spliced_amino_start[p-1] || 
          ((path->revcomp    && first_edge->upstream_spliced_nuc_end >= path->downstream_spliced_nuc_start[p-1]) ||
          ((!path->revcomp)  && first_edge->upstream_spliced_nuc_end <= path->downstream_spliced_nuc_start[p-1])))
        {
          if(first_edge != NULL) free(first_edge);
          first_edge = NULL;
          if(last_edge != NULL) free(last_edge);
          last_edge = NULL;

          for (i = 0; i < num_hits; i++) {
            p7_alidisplay_Destroy(split_hits[i]->dcl->ad);
            p7_trace_fs_Destroy(split_hits[i]->dcl->tr);
            free(split_hits[i]->dcl->scores_per_pos);
            p7_hit_Destroy(split_hits[i]);
          }
          free(split_hits);
          split_hits = NULL;
          printf("First Edge Backward %d\n", p+1);
          p++;
          continue;
 
        }
      }
   
      if(last_edge != NULL) {
        if(last_edge->downstream_spliced_amino_start >= path->upstream_spliced_amino_end[p+2] ||
          ((path->revcomp    && last_edge->downstream_spliced_nuc_start <= path->upstream_spliced_nuc_end[p+2]) ||
          ((!path->revcomp)  && last_edge->downstream_spliced_nuc_start   >= path->upstream_spliced_nuc_end[p+2])))
        {
          if(first_edge != NULL) free(first_edge);
          first_edge = NULL;
          if(last_edge != NULL) free(last_edge);
          last_edge = NULL;
         
           for (i = 0; i < num_hits; i++) {
            p7_alidisplay_Destroy(split_hits[i]->dcl->ad);
            p7_trace_fs_Destroy(split_hits[i]->dcl->tr);
            free(split_hits[i]->dcl->scores_per_pos);
            p7_hit_Destroy(split_hits[i]);
          }
          free(split_hits);
          split_hits = NULL;
          printf("Last Edge Backward %d\n", p+1);
          p++;
          continue;

        }
      }

      ESL_ALLOC(split_edges, sizeof(SPLICE_EDGE*) * (num_hits-1));
      for(i = 0; i < num_hits-1; i++) split_edges[i] = NULL;

      prev_edge = first_edge;
      for (h = 0; h < num_hits-1; h++) {
 
        curr_edge = NULL;
        curr_edge = splice_edge_create();
       
        curr_edge->overlap_amino_start = ESL_MAX(split_hits[h]->dcl->jhmm   - MIN_AMINO_OVERLAP, (prev_edge == NULL ? split_hits[h]->dcl->ihmm   : prev_edge->downstream_spliced_amino_start+1));
        curr_edge->overlap_amino_end   = ESL_MIN(split_hits[h+1]->dcl->ihmm + MIN_AMINO_OVERLAP, (last_edge == NULL ? split_hits[h+1]->dcl->jhmm : last_edge->upstream_spliced_amino_end-1));

        get_overlap_nuc_coords(curr_edge, split_hits[h]->dcl, split_hits[h+1]->dcl, path_seq, path->revcomp);         

        /* Make sure overlaps do not extend beyond existing splice boundries */
        if(path->revcomp) {
          if(prev_edge != NULL) curr_edge->upstream_nuc_start   = ESL_MAX(curr_edge->upstream_nuc_start,   path_seq->n + path_seq->end - prev_edge->downstream_spliced_nuc_start);
          if(last_edge != NULL) curr_edge->downstream_nuc_end   = ESL_MIN(curr_edge->downstream_nuc_end,   path_seq->n + path_seq->end - last_edge->upstream_spliced_nuc_end);
        }
        else {
          if(prev_edge != NULL) curr_edge->upstream_nuc_start   = ESL_MAX(curr_edge->upstream_nuc_start,   prev_edge->downstream_spliced_nuc_start - path_seq->start + 1);
          if(last_edge != NULL) curr_edge->downstream_nuc_end   = ESL_MIN(curr_edge->downstream_nuc_end,   last_edge->upstream_spliced_nuc_end     - path_seq->start + 1);
        } 
         
        /* If we encounter a backwards edge abandon ship */
        if(curr_edge->upstream_nuc_end <= curr_edge->upstream_nuc_start || curr_edge->downstream_nuc_end <= curr_edge->downstream_nuc_start) {

          if(first_edge != NULL) free(first_edge);
          if(last_edge != NULL) free(last_edge);
          free(curr_edge);
 
          for(i = 0; i < num_hits-1; i++) {
            if(split_edges[i] != NULL) free(split_edges[i]);
          }
          free(split_edges);
          split_edges = NULL;
     
          for (i = 0; i < num_hits; i++) {
            p7_alidisplay_Destroy(split_hits[i]->dcl->ad);
            p7_trace_fs_Destroy(split_hits[i]->dcl->tr);
            free(split_hits[i]->dcl->scores_per_pos);
            p7_hit_Destroy(split_hits[i]);
          }
          free(split_hits); 
          split_hits = NULL;
          printf("Backward Edge %d\n", p+1);
          num_hits = 0;
          break;           
        }
        else {
        
          if ((status = find_optimal_splice_site (curr_edge, split_hits[h]->dcl, split_hits[h+1]->dcl, gm, hmm, bg, gcode, path_seq, TRUE)) != eslOK) goto ERROR;
           
          /* If splice site could not be found - abandon ship */        
          if(curr_edge->splice_score == -eslINFINITY) {
            if(first_edge != NULL) free(first_edge);
            if(last_edge != NULL) free(last_edge);
            free(curr_edge);

            for(i = 0; i < num_hits-1; i++) {
              if(split_edges[i] != NULL) free(split_edges[i]);
            }
            free(split_edges);
            split_edges = NULL;

            for (i = 0; i < num_hits; i++) {
              p7_alidisplay_Destroy(split_hits[i]->dcl->ad);
              p7_trace_fs_Destroy(split_hits[i]->dcl->tr);
              free(split_hits[i]->dcl->scores_per_pos);
              p7_hit_Destroy(split_hits[i]);
            }
            free(split_hits);
            split_hits = NULL;

            printf("No Splice Found %d\n", p+1);
            num_hits = 0;
            break;
          }   
          
          if(path->revcomp) {
            curr_edge->upstream_spliced_nuc_end     = path_seq->n - curr_edge->upstream_spliced_nuc_end     + path_seq->end;
            curr_edge->downstream_spliced_nuc_start = path_seq->n - curr_edge->downstream_spliced_nuc_start + path_seq->end;
          }
          else {
            curr_edge->upstream_spliced_nuc_end     = curr_edge->upstream_spliced_nuc_end     + path_seq->start - 1;
            curr_edge->downstream_spliced_nuc_start = curr_edge->downstream_spliced_nuc_start + path_seq->start - 1;
          } 
  
          curr_edge->upstream_node_id   = path->node_id[p];
          curr_edge->downstream_node_id = path->node_id[p];
          curr_edge->target_seq_start   = path_seq->start;
          curr_edge->target_seq_end     = path_seq->end;
          curr_edge->target_seq_n       = path_seq->n;
          split_edges[h] = curr_edge;
          prev_edge = curr_edge;
        }

      }
   
      if(num_hits == 0) { p++; continue; }

      /* Reset splice boudries to upstream and dowstream neighbors */
      if(p != 0) {
        path->upstream_spliced_amino_end[p]     = first_edge->upstream_spliced_amino_end;
        path->downstream_spliced_amino_start[p] = first_edge->downstream_spliced_amino_start;
        path->upstream_spliced_nuc_end[p]       = first_edge->upstream_spliced_nuc_end;
        path->downstream_spliced_nuc_start[p]   = first_edge->downstream_spliced_nuc_start;
      }
      if(first_edge != NULL) free(first_edge);

      if(p != path->path_len-1) {
        path->upstream_spliced_amino_end[p+1]     = last_edge->upstream_spliced_amino_end;
        path->downstream_spliced_amino_start[p+1] = last_edge->downstream_spliced_amino_start;
        path->upstream_spliced_nuc_end[p+1]       = last_edge->upstream_spliced_nuc_end;
        path->downstream_spliced_nuc_start[p+1]   = last_edge->downstream_spliced_nuc_start;
      }
      if(last_edge != NULL) free(last_edge);
 
      /* Transfer frameshifts and stops so it goes to the right alignment algorithm */
      split_hits[0]->dcl->ad->frameshifts = path->hits[p]->dcl->ad->frameshifts;
      split_hits[0]->dcl->ad->stops       = path->hits[p]->dcl->ad->stops;      

      for (i = 0; i < num_hits-1; i++) 
        if ((status = splice_path_add_hit(path, split_edges[i], split_hits[i], split_hits[i+1], p+i)) != eslOK) goto ERROR;
    
      free(split_hits);
      split_hits = NULL;
  
      for(i = 0; i < num_hits-1; i++) {
        if(split_edges[i] != NULL) free(split_edges[i]);
      }
      free(split_edges);
      split_edges = NULL;

      p += (num_hits-1);
    }
    //path_dump(stdout, path); 
    p++;
  }

  if (gap != NULL) free(gap);
  return eslOK;

  ERROR:
    if(split_hits != NULL) free(split_hits);
    if(split_edges != NULL) free(split_edges);
    if (gap != NULL) free(gap);
    return status;
}


int
splice_path_add_hit(SPLICE_PATH *path, SPLICE_EDGE *edge, P7_HIT *up_hit, P7_HIT *down_hit, int path_step)
{

  int i;

  splice_path_grow(path);

  path->path_len++;
  path->split[path_step] = TRUE;

  /* shift path steps down to make room for split hit */
  for (i = path->path_len; i > path_step; i--) {
    path->upstream_spliced_amino_end[i]     = path->upstream_spliced_amino_end[i-1];
    path->downstream_spliced_amino_start[i] = path->downstream_spliced_amino_start[i-1];
    path->upstream_spliced_nuc_end[i]       = path->upstream_spliced_nuc_end[i-1];
    path->downstream_spliced_nuc_start[i]   = path->downstream_spliced_nuc_start[i-1];

    if(i < path->path_len) {
      path->node_id[i]                      = path->node_id[i-1];
      path->split[i]                        = path->split[i-1];
      path->hit_scores[i]                   = path->hit_scores[i-1];
      path->edge_scores[i]                  = path->edge_scores[i-1];
      path->signal_scores[i]                = path->signal_scores[i-1];
      path->hits[i]                         = path->hits[i-1];
    }
  }

  /*Replace upstream hit */
  path->hits[path_step]       = up_hit;
  path->hit_scores[path_step] = up_hit->dcl->aliscore;  
  path->split[path_step]      = TRUE;

  /*Replace downstream hit */
  path->hits[path_step+1]       = down_hit;
  path->hit_scores[path_step+1] = down_hit->dcl->aliscore;
  path->split[path_step+1]      = TRUE;
  
  /* Set new splice boudries */
  path->upstream_spliced_amino_end[path_step+1]     = edge->upstream_spliced_amino_end;
  path->downstream_spliced_amino_start[path_step+1] = edge->downstream_spliced_amino_start;
  path->upstream_spliced_nuc_end[path_step+1]       = edge->upstream_spliced_nuc_end;
  path->downstream_spliced_nuc_start[path_step+1]   = edge->downstream_spliced_nuc_start;

  path->edge_scores[path_step+1]   = edge->splice_score;
  path->signal_scores[path_step+1] = edge->signal_score;
  
  return eslOK;

}

float
get_partial_ali_score (P7_HIT *hit, int start_k, int end_k) {

  int z;
  int z1, z2;
  int z_offset;
  float score;
  P7_TRACE *tr;
 
  tr = hit->dcl->tr;
 
  for (z1 = 0;       z1 < tr->N; z1++) if (tr->st[z1] == p7T_M) break;
  z_offset = z1;
   
  for (z1 = tr->N-1; z1 >= 0;    z1--) if (tr->k[z1] == start_k) break;
  for (z2 = 0;       z2 < tr->N; z2++) if (tr->k[z2] == end_k)   break;

  score = 0;
  z = z1;
  while(z <= z2) {
    score += hit->dcl->scores_per_pos[z-z_offset];
    z++;
  }
  
  return score;

}


int
splice_the_path (TARGET_RANGE *target_range, SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_TOPHITS *orig_tophits, P7_OPROFILE *om, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *path_seq, int64_t db_nuc_cnt, int *frameshift, int* success) 
{

  int          i;
  int          pos, seq_pos;
  int          exon;
  int          shift;
  int          seq_idx;
  int          amino_len;
  int          amino;
  int          orf_len;
  int          replace_node;
  int          remove_node;
  int          contains_orig;
  int          path_seq_len;
  int64_t      ali_L;
  float        dom_bias;
  float        nullsc;
  float        dom_score;
  double       dom_lnP;
  int         *nuc_index;
  ESL_DSQ     *nuc_dsq;
  ESL_DSQ     *amino_dsq;
  ESL_SQ      *nuc_seq;
  ESL_SQ      *amino_seq;  
  P7_HIT      *replace_hit;
  P7_HIT      *remove_hit;
  P7_TOPHITS  *th;
  int          status;

  th = target_range->th;

  nuc_index = NULL;
  nuc_dsq   = NULL;

  path_seq_len = 0;
  for (i = 0; i < path->path_len; i++)
    path_seq_len += abs(path->upstream_spliced_nuc_end[i+1] - path->downstream_spliced_nuc_start[i]) + 1;

  if(path_seq_len % 3 != 0) { 
    *frameshift = TRUE;
    *success    = FALSE;
    return eslOK;
  }
  
  ESL_ALLOC(nuc_index, sizeof(int64_t) * (path_seq_len+2));  
  ESL_ALLOC(nuc_dsq,   sizeof(ESL_DSQ) * (path_seq_len+2));

  /* Copy spliced nucleotides into single sequence and track their original indicies */
  nuc_index[0] = -1;
  nuc_dsq[0]   = eslDSQ_SENTINEL;
  seq_idx   = 1;
  
  for (i = 0; i < path->path_len; i++) {
    if (path->revcomp) {
      for (pos = path->downstream_spliced_nuc_start[i]; pos >= path->upstream_spliced_nuc_end[i+1]; pos--) {
        seq_pos = path_seq->n - pos + path_seq->end;
        
        nuc_index[seq_idx] = seq_pos;
        nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;
      }
    }
    else {
      for (pos = path->downstream_spliced_nuc_start[i]; pos <= path->upstream_spliced_nuc_end[i+1]; pos++) {
        seq_pos = pos - path_seq->start + 1;      

        nuc_index[seq_idx] = seq_pos;
        nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;    
      }
    }
  }
  
  nuc_index[seq_idx] = -1;
  nuc_dsq[seq_idx]   = eslDSQ_SENTINEL;


  /* Translate spliced nucleotide sequence */
  amino_len = path_seq_len/3;
  ESL_ALLOC(amino_dsq, sizeof(ESL_DSQ) * (amino_len+2));

  amino_dsq[0] = eslDSQ_SENTINEL;
   
  seq_idx = 1;

  for(i = 1; i <= amino_len; i++) {
   
    amino = esl_gencode_GetTranslation(gcode,&nuc_dsq[seq_idx]);

    if(esl_abc_XIsNonresidue(gm->abc, amino)) {
      *frameshift = TRUE;
      *success    = FALSE;
      free(nuc_index);
      free(nuc_dsq);
      free(amino_dsq);
      
      return eslOK;
    }

    amino_dsq[i] = amino;
    seq_idx+=3;
  }
  amino_dsq[amino_len+1] = eslDSQ_SENTINEL;
  
  /* Create ESL_SQs from ESL_DSQs */
  amino_seq   = esl_sq_CreateDigitalFrom(gcode->nt_abc, path_seq->name, amino_dsq, amino_len,    NULL,NULL,NULL); 
  nuc_seq     = esl_sq_CreateDigitalFrom(gcode->aa_abc, path_seq->name, nuc_dsq,   path_seq_len, NULL,NULL,NULL);
  
  pli->nuc_sq       = nuc_seq;
  pli->amino_sq     = amino_seq;   
  pli->orig_nuc_idx = nuc_index;

  /* Algin the splices Amino sequence */
  status = align_spliced_path(pli, om, gm, path_seq, gcode);
 
  /* Alignment failed */
  if(pli->hit == NULL || pli->hit->dcl->ad->exon_cnt == 1) { 

    if(nuc_dsq   != NULL) free(nuc_dsq);
    if(amino_dsq != NULL) free(amino_dsq);
 
    *success = FALSE; 
     return eslOK; 
  }
  
  /* adjust all coords in hit and path */
  if(path->revcomp) { 
    pli->hit->dcl->ad->sqfrom += 2; 
    pli->hit->dcl->ienv = path_seq->n - pli->orig_nuc_idx[1]              + path_seq->end; 
    pli->hit->dcl->jenv = path_seq->n - pli->orig_nuc_idx[pli->nuc_sq->n] + path_seq->end; 
  }  
  else {
    pli->hit->dcl->ienv = pli->orig_nuc_idx[1]              + path_seq->start -1;
    pli->hit->dcl->jenv = pli->orig_nuc_idx[pli->nuc_sq->n] + path_seq->start -1;
  }
  
  /* Adjust spliced hit score from amino_len to om->max_length */ 
  dom_score  = pli->hit->dcl->envsc;
  orf_len = pli->hit->dcl->ad->orfto - pli->hit->dcl->ad->orffrom + 1;
  dom_score -= 2 * log(2. / (amino_len+2));
  dom_score += 2 * log(2. / (om->max_length+2));
  dom_score -= (amino_len-orf_len)      * log((float) (amino_len) / (float) (amino_len+2));
  dom_score += (om->max_length-orf_len) * log((float) om->max_length / (float) (om->max_length+2));
  
  
  /* Bias calculation and adjustments */
  if(pli->do_null2)
    dom_bias = p7_FLogsum(0.0, log(pli->bg->omega) + pli->hit->dcl->domcorrection);
  else
    dom_bias = 0.;

  p7_bg_SetLength(pli->bg, om->max_length);
  p7_bg_NullOne  (pli->bg, pli->amino_sq->dsq, om->max_length, &nullsc);
  dom_score = (dom_score - (nullsc + dom_bias))  / eslCONST_LOG2; 

  /* Add splice signal penalties */
  for(i = 0; i < path->path_len; i++) 
     dom_score += path->signal_scores[i];

  dom_lnP   = esl_exp_logsurv(dom_score, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
  
  /* E-value adjusment */
  dom_lnP += log((float)db_nuc_cnt / (float)om->max_length);
  
  if ((pli->by_E && exp(dom_lnP) <= pli->E) || ((!pli->by_E) && dom_score >= pli->T)) {
    
    *success = TRUE;

    if ( path->path_len >  pli->hit->dcl->ad->exon_cnt) {
      /* Shift the path to start at the first hit that was inculded in the alignment 
       * and end at the last hit that was included in the alignment */ 
      if(path->revcomp) {
        for(shift = 0; shift < path->path_len; shift++) {
          if(path->upstream_spliced_nuc_end[shift+1] <= pli->hit->dcl->iali) break;
        }
      } 
      else {
        for(shift = 0; shift < path->path_len; shift++) {
           if(path->upstream_spliced_nuc_end[shift+1] >= pli->hit->dcl->iali) break;
        } 
      }
     
      /* Ensure path will still contains an original hit after shifting */
      contains_orig = FALSE;
      for(exon = shift; exon < pli->hit->dcl->ad->exon_cnt ; exon ++ ) {
        if(path->node_id[exon] < target_range->orig_N) contains_orig = TRUE;  
      } 

      if(!contains_orig) {
        *success = FALSE;
        if(nuc_dsq   != NULL) free(nuc_dsq);
        if(amino_dsq != NULL) free(amino_dsq);
        return eslOK;
      }

      /* Shift path to start at frist hits that is in alignment */
      path->path_len = pli->hit->dcl->ad->exon_cnt;
      
      for(exon = 0; exon < path->path_len; exon++) {
        path->node_id[exon]                        = path->node_id[shift+exon]; 
        path->upstream_spliced_amino_end[exon]     = path->upstream_spliced_amino_end[shift+exon];
        path->downstream_spliced_amino_start[exon] = path->downstream_spliced_amino_start[shift+exon];
        path->upstream_spliced_nuc_end[exon]       = path->upstream_spliced_nuc_end[shift+exon];   
        path->downstream_spliced_nuc_start[exon]   = path->downstream_spliced_nuc_start[shift+exon];
        path->hit_scores[exon]                     = path->hit_scores[shift+exon];
        path->edge_scores[exon]                    = path->edge_scores[shift+exon];
        path->hit_scores[exon]                     = path->hit_scores[shift+exon];
        path->edge_scores[exon]                    = path->edge_scores[shift+exon];
        path->signal_scores[exon]                  = path->signal_scores[shift+exon];
        path->hits[exon]                           = path->hits[shift+exon];
        path->split[exon]                          = path->split[shift+exon];
      }
      path->downstream_spliced_nuc_start[0]   = pli->hit->dcl->iali;
      path->downstream_spliced_amino_start[0] = pli->hit->dcl->ihmm; 
  
      for(exon = 1; exon < path->path_len; exon++) {
        if (path->downstream_spliced_nuc_start[exon] > pli->hit->dcl->jali)
          break;
      }
     
      path->path_len =  pli->hit->dcl->ad->exon_cnt;
      
      path->upstream_spliced_nuc_end[exon]   = pli->hit->dcl->jali;
      path->upstream_spliced_amino_end[exon] = pli->hit->dcl->jhmm;   
    }  

    /* Replace the first hit in path with the spliced hit 
     * and set all other hits in path to unreportable */ 
    i = 0;
    while( path->node_id[i] >= target_range->orig_N) {
      pli->hit->dcl->ad->exon_orig[i] = FALSE; 
      pli->hit->dcl->ad->exon_split[i] = path->split[i];
      i++; 
    }

    pli->hit->dcl->ad->exon_orig[i] = TRUE; 
    pli->hit->dcl->ad->exon_split[i] = path->split[i];
    replace_node = path->node_id[i];
    replace_hit = th->hit[replace_node]; 
	ali_L = replace_hit->dcl->ad->L;
    p7_domain_Destroy(replace_hit->dcl);

    replace_hit->dcl        = pli->hit->dcl;
	replace_hit->dcl->ad->L = ali_L;
    replace_hit->frameshift = FALSE;
 
    replace_hit->flags =  0;
    replace_hit->flags |= p7_IS_REPORTED; 
    replace_hit->flags |= p7_IS_INCLUDED;
    replace_hit->nreported = 1;
    replace_hit->nincluded = 1;

    replace_hit->dcl->bitscore    = dom_score;
    replace_hit->dcl->lnP         = dom_lnP;
    replace_hit->dcl->dombias     = dom_bias;
    replace_hit->dcl->is_reported = TRUE;
    replace_hit->dcl->is_included = TRUE; 
    
    replace_hit->pre_score = pli->hit->dcl->envsc  / eslCONST_LOG2;
    replace_hit->pre_lnP   = esl_exp_logsurv (replace_hit->pre_score, om->evparam[p7_FTAUFS], om->evparam[p7_FLAMBDA]);
    
    replace_hit->sum_score  = replace_hit->score  = dom_score; 
    replace_hit->sum_lnP    = replace_hit->lnP    = dom_lnP;
  
    replace_hit->sortkey    = pli->inc_by_E ? -dom_lnP : dom_score;
  
    /* Set all other hits in alignment to unreportable */ 
    i++;
    for(  ; i < path->path_len; i++) {
      if(path->node_id[i] >= target_range->orig_N)
        pli->hit->dcl->ad->exon_orig[i] = FALSE;
      else
        pli->hit->dcl->ad->exon_orig[i] = TRUE;

      pli->hit->dcl->ad->exon_split[i] = path->split[i]; 
      /* If the replace node hes bee split we neeed to 
       * make sure we don't set is as unreportable */
      remove_node = path->node_id[i];
      if(remove_node == replace_node)
         continue;

      remove_hit = th->hit[remove_node];

      if(remove_hit->flags & p7_IS_REPORTED ) {
        orig_tophits->nreported--;
        remove_hit->flags &= ~p7_IS_REPORTED;
        remove_hit->dcl->is_reported = FALSE;
      }
      if((remove_hit->flags & p7_IS_INCLUDED)) {
        orig_tophits->nincluded--;     
        remove_hit->flags &= ~p7_IS_INCLUDED;
        remove_hit->dcl->is_included = FALSE;
      }
    }
    
  } 
  else *success = FALSE; 

  if(nuc_dsq   != NULL) free(nuc_dsq);
  if(amino_dsq != NULL) free(amino_dsq);  
  
  return eslOK;
 
  ERROR:
    if(nuc_index != NULL) free(nuc_index);
    if(nuc_dsq   != NULL) free(nuc_dsq);
    if(amino_dsq != NULL) free(amino_dsq);
    return status;

}

int
splice_the_path_frameshift (TARGET_RANGE *target_range, SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_TOPHITS *orig_tophits, P7_FS_PROFILE *gm_fs, P7_OPROFILE *om, ESL_GENCODE *gcode, ESL_SQ *path_seq, int64_t db_nuc_cnt, int* success)
{

  int          i;
  int          pos, seq_pos;
  int          exon;
  int          shift;
  int          seq_idx;
  int          ali_len;
  int          replace_node;
  int          remove_node;
  int          contains_orig;
  int          path_seq_len;
  int64_t      ali_L;
  float        dom_bias;
  float        nullsc;
  float        dom_score;
  double       dom_lnP;
  int         *nuc_index;
  ESL_DSQ     *nuc_dsq;
  ESL_SQ      *nuc_seq;
  P7_HIT      *replace_hit;
  P7_HIT      *remove_hit;
  P7_TOPHITS  *th;
  int          status;

  th = target_range->th;

  nuc_index = NULL;
  nuc_dsq   = NULL;

  path_seq_len = 0;
  for (i = 0; i < path->path_len; i++)
    path_seq_len += abs(path->upstream_spliced_nuc_end[i+1] - path->downstream_spliced_nuc_start[i]) + 1;
  
  ESL_ALLOC(nuc_index, sizeof(int64_t) * (path_seq_len+2));
  ESL_ALLOC(nuc_dsq,   sizeof(ESL_DSQ) * (path_seq_len+2));
 
  /* Copy spliced nucleotides into single sequence and track their original indicies */
  nuc_index[0] = -1;
  nuc_dsq[0]   = eslDSQ_SENTINEL;
  seq_idx   = 1;
   
  for (i = 0; i < path->path_len; i++) {
    if (path->revcomp) {
      for (pos = path->downstream_spliced_nuc_start[i]; pos >= path->upstream_spliced_nuc_end[i+1]; pos--) {
        seq_pos = path_seq->n - pos + path_seq->end;

        nuc_index[seq_idx] = seq_pos;
        nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;
      }
    }
    else {
      for (pos = path->downstream_spliced_nuc_start[i]; pos <= path->upstream_spliced_nuc_end[i+1]; pos++) {
        seq_pos = pos - path_seq->start + 1;

        nuc_index[seq_idx] = seq_pos;
        nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;
      }
    } 
  }

  nuc_index[seq_idx] = -1;
  nuc_dsq[seq_idx]   = eslDSQ_SENTINEL;
 
  nuc_seq = esl_sq_CreateDigitalFrom(gcode->aa_abc, path_seq->name, nuc_dsq, path_seq_len, NULL,NULL,NULL);

  pli->nuc_sq       = nuc_seq;
  pli->amino_sq     = NULL;
  pli->orig_nuc_idx = nuc_index;

  /* Algin the spliced Amino sequence */
  status = align_spliced_path_frameshift(pli, gm_fs, om, path_seq, gcode);

  /* Alignment failed */
  if(pli->hit == NULL ) {
    if(nuc_dsq   != NULL) free(nuc_dsq);
    *success = FALSE;
     return eslOK;
  }

  /* adjust all coords in hit and path */
  if(path->revcomp) {
    pli->hit->dcl->ad->sqfrom += 2;
    pli->hit->dcl->ienv = path_seq->n - pli->orig_nuc_idx[1]              + path_seq->end;
    pli->hit->dcl->jenv = path_seq->n - pli->orig_nuc_idx[pli->nuc_sq->n] + path_seq->end;
  }
  else {
    pli->hit->dcl->ienv = pli->orig_nuc_idx[1]              + path_seq->start -1;
    pli->hit->dcl->jenv = pli->orig_nuc_idx[pli->nuc_sq->n] + path_seq->start -1;
  }
  
  
  ali_len = ESL_MAX(pli->hit->dcl->iali, pli->hit->dcl->jali) - ESL_MIN(pli->hit->dcl->iali, pli->hit->dcl->jali) + 1;
  /* Adjust spliced hit score from nuc_len to gm_fs->max_length*3 */
  dom_score  = pli->hit->dcl->envsc;
  dom_score -= 2 * log(2. / ((pli->nuc_sq->n/3.)+2));
  dom_score += 2 * log(2. / (gm_fs->max_length+2));
  dom_score -= (pli->nuc_sq->n-ali_len)      * log((float) (pli->nuc_sq->n/3.) / (float) ((pli->nuc_sq->n/3.)+2));
  dom_score += (ESL_MAX(pli->nuc_sq->n, gm_fs->max_length*3)-ali_len) * log((float) gm_fs->max_length / (float) (gm_fs->max_length+2));

  /* Bias calculation and adjustments */
  if(pli->do_null2)
    dom_bias = p7_FLogsum(0.0, log(pli->bg->omega) + pli->hit->dcl->domcorrection);
  else
    dom_bias = 0.;

  p7_bg_SetLength(pli->bg, gm_fs->max_length*3);
  p7_bg_NullOne  (pli->bg, pli->nuc_sq->dsq, gm_fs->max_length*3, &nullsc);
  dom_score = (dom_score - (nullsc + dom_bias))  / eslCONST_LOG2;

  /* Add splice signal penalties */
  pli->hit->dcl->jali = pli->hit->dcl->ad->sqto;
  for(i = 0; i < path->path_len; i++)
    dom_score += path->signal_scores[i];

  dom_lnP   = esl_exp_logsurv(dom_score, gm_fs->evparam[p7_FTAUFS], gm_fs->evparam[p7_FLAMBDA]);

  /* E-value adjusment */
  dom_lnP += log((float)db_nuc_cnt / (float)gm_fs->max_length);

  /* If any of the original hits in the path have a lower eval (or higher score)
   * than the spliced alignment then do not report the spliced alignment */
  for(i = 0; i < path->path_len; i++) {
    if(path->node_id[i] < target_range->orig_N) {
      if(pli->by_E && path->hits[i]->dcl->lnP < dom_lnP)
        dom_lnP = eslINFINITY;
      else if((!pli->by_E) && path->hits[i]->dcl->bitscore > dom_score)
        dom_score = -eslINFINITY;
    }
  }

  if ((pli->by_E && exp(dom_lnP) <= pli->E) || ((!pli->by_E) && dom_score >= pli->T)) {

    *success = TRUE;

    if ( path->path_len > pli->hit->dcl->ad->exon_cnt) {
      /* Shift the path to start at the first hit that was inculded in the alignment
       * and end at the last hit that was included in the alignment */
      if(path->revcomp) {
        for(shift = 0; shift < path->path_len; shift++) {
          if(path->upstream_spliced_nuc_end[shift+1] <= pli->hit->dcl->iali) break;
        }
      } 
      else {
        for(shift = 0; shift < path->path_len; shift++) {
           if(path->upstream_spliced_nuc_end[shift+1] >= pli->hit->dcl->iali) break;
        } 
      }
      /* Ensure path will still contains an original hit after shifting */
      contains_orig = FALSE;
      for(exon = shift; exon < pli->hit->dcl->ad->exon_cnt ; exon ++ ) {
        if(path->node_id[exon] < target_range->orig_N) contains_orig = TRUE;  
      } 

      if(!contains_orig) {
        *success = FALSE;
        if(nuc_dsq   != NULL) free(nuc_dsq);
        return eslOK;
      }

      /* Shift path to start at frist hit that is in alignment */
      path->path_len = pli->hit->dcl->ad->exon_cnt;
  
      for(exon = 0; exon < path->path_len; exon++) {
        path->node_id[exon]                        = path->node_id[shift+exon];
        path->upstream_spliced_amino_end[exon]     = path->upstream_spliced_amino_end[shift+exon];
        path->downstream_spliced_amino_start[exon] = path->downstream_spliced_amino_start[shift+exon];
        path->upstream_spliced_nuc_end[exon]       = path->upstream_spliced_nuc_end[shift+exon];
        path->downstream_spliced_nuc_start[exon]   = path->downstream_spliced_nuc_start[shift+exon];
        path->hit_scores[exon]                     = path->hit_scores[shift+exon];
        path->edge_scores[exon]                    = path->edge_scores[shift+exon];
        path->signal_scores[exon]                  = path->signal_scores[shift+exon];
        path->hits[exon]                           = path->hits[shift+exon];
        path->split[exon]                          = path->split[shift+exon];
      }
      path->downstream_spliced_nuc_start[0]   = pli->hit->dcl->iali;
      path->downstream_spliced_amino_start[0] = pli->hit->dcl->ihmm;
      
      for(exon = 1; exon < path->path_len; exon++) {
        if (path->downstream_spliced_nuc_start[exon] > pli->hit->dcl->jali)
          break;
      }
  
      path->upstream_spliced_nuc_end[exon]   = pli->hit->dcl->jali;
      path->upstream_spliced_amino_end[exon] = pli->hit->dcl->jhmm;
    }

    /* Replace the first hit in path with the spliced hit
     * and set all other hits in path to unreportable */
    i = 0;
    while( path->node_id[i] >= target_range->orig_N) {

      pli->hit->dcl->ad->exon_orig[i] = FALSE;
      pli->hit->dcl->ad->exon_split[i] = path->split[i];
      i++;
    }

    pli->hit->dcl->ad->exon_orig[i] = TRUE;
    pli->hit->dcl->ad->exon_split[i] = path->split[i];
    
    replace_node = path->node_id[i];
    replace_hit  = th->hit[replace_node];

    ali_L = replace_hit->dcl->ad->L;
    p7_domain_Destroy(replace_hit->dcl);

    replace_hit->dcl        = pli->hit->dcl;
    replace_hit->dcl->ad->L = ali_L;
    replace_hit->frameshift = TRUE;

    replace_hit->flags =  0;
    replace_hit->flags |= p7_IS_REPORTED;
    replace_hit->flags |= p7_IS_INCLUDED;
    replace_hit->nreported = 1;
    replace_hit->nincluded = 1;

    replace_hit->dcl->bitscore    = dom_score;
    replace_hit->dcl->lnP         = dom_lnP;
    replace_hit->dcl->dombias     = dom_bias;
    replace_hit->dcl->is_reported = TRUE;
    replace_hit->dcl->is_included = TRUE;

    replace_hit->pre_score = pli->hit->dcl->envsc  / eslCONST_LOG2;
    replace_hit->pre_lnP   = esl_exp_logsurv (replace_hit->pre_score, gm_fs->evparam[p7_FTAUFS], gm_fs->evparam[p7_FLAMBDA]);

    replace_hit->sum_score  = replace_hit->score  = dom_score;
    replace_hit->sum_lnP    = replace_hit->lnP    = dom_lnP;

    replace_hit->sortkey    = pli->inc_by_E ? -dom_lnP : dom_score;

    /* Set all other hits in alignment to unreportable */
    i++;
    for(  ; i < path->path_len; i++) {

      if(i < pli->hit->dcl->ad->exon_cnt) {
        if(path->node_id[i] >= target_range->orig_N)
          pli->hit->dcl->ad->exon_orig[i] = FALSE;
        else
          pli->hit->dcl->ad->exon_orig[i] = TRUE;
      }
      pli->hit->dcl->ad->exon_split[i] = path->split[i];
      remove_node = path->node_id[i];
      if(remove_node == replace_node)
        continue;
      remove_hit = th->hit[remove_node];

      if(remove_hit->flags & p7_IS_REPORTED ) {
        orig_tophits->nreported--;
        remove_hit->flags &= ~p7_IS_REPORTED;
        remove_hit->dcl->is_reported = FALSE;
      }
      if((remove_hit->flags & p7_IS_INCLUDED)) {
        orig_tophits->nincluded--;
        remove_hit->flags &= ~p7_IS_INCLUDED;
        remove_hit->dcl->is_included = FALSE;
      }
    }
  }
  else  *success = FALSE; 

  if(nuc_dsq   != NULL) free(nuc_dsq);

  return eslOK;

  ERROR:
    if(nuc_index != NULL) free(nuc_index);
    if(nuc_dsq   != NULL) free(nuc_dsq);
    return status;

}

int 
align_spliced_path (SPLICE_PIPELINE *pli, P7_OPROFILE *om, P7_PROFILE *gm, ESL_SQ *target_seq, ESL_GENCODE *gcode) 
{

  int       i;
  int       splice_cnt;
  float     filtersc;
  float     envsc; 
  float     seq_score;
  float     oasc;
  float     P;
  float     domcorrection;  
  float    null2[p7_MAXCODE];
  P7_HIT   *hit;
  P7_TRACE *tr;
  int       status;
  
  hit          = p7_hit_Create_empty();
  hit->dcl     = p7_domain_Create_empty();
  hit->dcl->tr = NULL;
  hit->dcl->scores_per_pos = NULL;
  tr = p7_trace_CreateWithPP();
 
  p7_oprofile_ReconfigUnihit(om, pli->amino_sq->n);
  p7_omx_GrowTo(pli->fwd, om->M, pli->amino_sq->n, pli->amino_sq->n);
  p7_omx_GrowTo(pli->bwd, om->M, pli->amino_sq->n, pli->amino_sq->n);

  p7_bg_SetLength(pli->bg, pli->amino_sq->n);
  if (pli->do_biasfilter)
    p7_bg_FilterScore(pli->bg, pli->amino_sq->dsq, pli->amino_sq->n, &filtersc);
  else
    p7_bg_NullOne  (pli->bg, pli->amino_sq->dsq, pli->amino_sq->n, &filtersc);

  p7_Forward (pli->amino_sq->dsq, pli->amino_sq->n, om, pli->fwd, &envsc);

  seq_score = (envsc-filtersc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
  
  if (P > pli->F3) {
    p7_hit_Destroy(hit);
    p7_trace_Destroy(tr);   
    return eslOK;
  }
  p7_Backward(pli->amino_sq->dsq, pli->amino_sq->n, om, pli->fwd, pli->bwd, NULL);  

  if((p7_Decoding(om, pli->fwd, pli->bwd, pli->bwd)) == eslERANGE) 
    ESL_XEXCEPTION(eslFAIL, "p7_Decoding failed with eslERANGE");

  p7_OptimalAccuracy(om, pli->bwd, pli->fwd, &oasc);
  p7_OATrace        (om, pli->bwd, pli->fwd, tr);
  
  p7_trace_Index(tr);
   
  p7_splice_compute_ali_scores(hit->dcl, tr, pli->amino_sq->dsq, gm, gm->abc->Kp);
 //p7_trace_Dump(stdout, tr, NULL, NULL); 
  if((hit->dcl->tr = p7_trace_splice_Convert(tr, pli->orig_nuc_idx, &splice_cnt)) == NULL) goto ERROR; 
  //p7_trace_Dump(stdout, hit->dcl->tr, NULL, NULL);
  if((hit->dcl->ad = p7_alidisplay_splice_Create(hit->dcl->tr, 0, om, target_seq, pli->amino_sq, hit->dcl->scores_per_pos, tr->sqfrom[0], splice_cnt)) == NULL) goto ERROR; 

  p7_Null2_ByExpectation(om, pli->bwd, null2);
  domcorrection = 0.;
  for (i = 1; i <= pli->amino_sq->n; i++)
    domcorrection += logf(null2[pli->amino_sq->dsq[i]]);
 
  hit->dcl->domcorrection = domcorrection;
 
  hit->dcl->ihmm = hit->dcl->ad->hmmfrom;
  hit->dcl->jhmm = hit->dcl->ad->hmmto;

  /* Convert sqfrom, sqto to full sequence coords */
  if(target_seq->start < target_seq->end) {
    hit->dcl->ad->sqfrom  = hit->dcl->ad->sqfrom + target_seq->start - 1;
    hit->dcl->ad->sqto    = hit->dcl->ad->sqto   + target_seq->start - 1;
  } else {
    hit->dcl->ad->sqto    = target_seq->n - hit->dcl->ad->sqto   + target_seq->end;
    hit->dcl->ad->sqfrom  = target_seq->n - hit->dcl->ad->sqfrom + target_seq->end;
  }

  hit->dcl->iali = hit->dcl->ad->sqfrom;
  hit->dcl->jali = hit->dcl->ad->sqto;
 
  hit->dcl->envsc = envsc;
  hit->dcl->oasc  = oasc;
  hit->dcl->dombias       = 0.0;
  hit->dcl->bitscore      = 0.0;
  hit->dcl->lnP           = 0.0;
  hit->dcl->is_reported   = FALSE;  
  hit->dcl->is_included   = FALSE;

  pli->hit = hit;

  p7_trace_Destroy(tr);
  return eslOK;

  ERROR:
    p7_trace_Destroy(tr);
    p7_hit_Destroy(hit);
    return status;
}

int
align_spliced_path_frameshift (SPLICE_PIPELINE *pli, P7_FS_PROFILE *gm_fs, P7_OPROFILE *om, ESL_SQ *target_seq, ESL_GENCODE *gcode)
{

  int       i, z;
  int       t, u, v, w, x;
  int       splice_cnt;
  int       codon_idx;
  float     nullsc;
  float     envsc;
  float     seq_score;
  float     oasc;
  float     P;
  float     domcorrection;
  float    null2[p7_MAXCODE];
  P7_GMX   *gxppfs;
  P7_HIT   *hit;
  P7_TRACE *tr;
  int       status;

  tr = p7_trace_fs_CreateWithPP();

  p7_fs_ReconfigUnihit(gm_fs, pli->nuc_sq->n);
  p7_gmx_fs_GrowTo(pli->gfwd, gm_fs->M, pli->nuc_sq->n, pli->nuc_sq->n, p7P_CODONS);
  p7_gmx_fs_GrowTo(pli->gbwd, gm_fs->M, pli->nuc_sq->n, pli->nuc_sq->n, p7P_CODONS);

  p7_bg_SetLength(pli->bg, pli->nuc_sq->n);
  p7_bg_NullOne  (pli->bg, pli->nuc_sq->dsq, pli->nuc_sq->n, &nullsc);

  
  p7_Forward_Frameshift(pli->nuc_sq->dsq, gcode, pli->nuc_sq->n, gm_fs, pli->gfwd, &envsc);

  seq_score = (envsc-nullsc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score, gm_fs->evparam[p7_FTAUFS],  gm_fs->evparam[p7_FLAMBDA]);

  if (P > pli->F3) {
    
    p7_trace_fs_Destroy(tr);
    return eslOK;
  }

  p7_Backward_Frameshift(pli->nuc_sq->dsq, gcode, pli->nuc_sq->n, gm_fs, pli->gbwd, NULL);

  if ((gxppfs = p7_gmx_fs_Create(gm_fs->M, pli->nuc_sq->n, pli->nuc_sq->n, p7P_CODONS)) == NULL) { status = eslFAIL; goto ERROR; }
  p7_Decoding_Frameshift(gm_fs, pli->gfwd, pli->gbwd, gxppfs);

  p7_OptimalAccuracy_Frameshift(gm_fs, gxppfs, pli->gbwd, &oasc);
  p7_OATrace_Frameshift(gm_fs, gxppfs, pli->gbwd, pli->gfwd, tr);   /* <tr>'s seq coords are offset by i-1, rel to orig dsq */

  p7_trace_Index(tr);

  hit          = p7_hit_Create_empty();
  hit->dcl     = p7_domain_Create_empty();
  hit->dcl->tr = NULL;

  p7_splice_compute_ali_scores_fs(hit->dcl, tr, pli->nuc_sq->dsq, gm_fs, target_seq->abc);

  if((hit->dcl->tr = p7_trace_splice_fs_Convert(tr, pli->orig_nuc_idx, &splice_cnt)) == NULL) { status = eslFAIL; goto ERROR; }
  
  if((hit->dcl->ad = p7_alidisplay_splice_fs_Create(hit->dcl->tr, 0, gm_fs, target_seq, pli->nuc_sq->dsq, gcode, hit->dcl->scores_per_pos, pli->orig_nuc_idx, tr->sqfrom[0], splice_cnt)) == NULL) { status = eslFAIL; goto ERROR; }

  p7_Null2_fs_ByExpectation(gm_fs, pli->gfwd, null2);
  domcorrection = 0.;
  t = u = v = w = x = -1;
  z = 0;
  i = 1;
  
  while(i <= pli->nuc_sq->n) {
    if(esl_abc_XIsCanonical(target_seq->abc, pli->nuc_sq->dsq[i])) x = pli->nuc_sq->dsq[i];
    else                                                           x = p7P_MAXCODONS;
  
    switch (tr->st[z]) {
      case p7T_N:
      case p7T_C:
      case p7T_J:  if(tr->i[z] == i && i > 2) i++;
                   z++;   break;
      case p7T_X:
      case p7T_S:
      case p7T_B:
      case p7T_E:
      case p7T_T:
      case p7T_D:  z++;   break;
      case p7T_M:  if(tr->i[z] == i)
                   {
                     if     (tr->c[z] == 1) { codon_idx = p7P_CODON1(x);             codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_QC2); }
                     else if(tr->c[z] == 2) { codon_idx = p7P_CODON2(w, x);          codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_QC1); }
                     else if(tr->c[z] == 3) { codon_idx = p7P_CODON3(v, w, x);       codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_C); }
                     else if(tr->c[z] == 4) { codon_idx = p7P_CODON4(u, v, w, x);    codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_QC1); }
                     else if(tr->c[z] == 5) { codon_idx = p7P_CODON5(t, u, v, w, x); codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_QC2); }
                     domcorrection += logf(null2[p7P_AMINO(gm_fs, tr->k[z], codon_idx)]);
                     z++;
                   }
                   i++;  break;
      case p7T_I:  if(tr->i[z] == i)
                   {
                     codon_idx = p7P_CODON3(v, w, x);
                     codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_C);
                     domcorrection += logf(null2[p7P_AMINO(gm_fs, tr->k[z], codon_idx)]);
                     z++;
                   }
                   i++;  break;
    }
    t = u;
    u = w;
    v = w;
    w = x;
  }
  
  hit->dcl->domcorrection = domcorrection;

  hit->dcl->ihmm = hit->dcl->ad->hmmfrom;
  hit->dcl->jhmm = hit->dcl->ad->hmmto;
  
  /* Convert sqfrom, sqto to full sequence coords */
  if(target_seq->start < target_seq->end) {
    hit->dcl->ad->sqfrom  = hit->dcl->ad->sqfrom + target_seq->start - 1;
    hit->dcl->ad->sqto    = hit->dcl->ad->sqto   + target_seq->start - 1;
  } else {
    hit->dcl->ad->sqto    = target_seq->n - hit->dcl->ad->sqto   + target_seq->end;
    hit->dcl->ad->sqfrom  = target_seq->n - hit->dcl->ad->sqfrom + target_seq->end;
  }
  hit->dcl->iali = hit->dcl->ad->sqfrom;
  hit->dcl->jali = hit->dcl->ad->sqto;

  hit->dcl->envsc = envsc;
  hit->dcl->oasc  = oasc;
  hit->dcl->dombias       = 0.0;
  hit->dcl->bitscore      = 0.0;
  hit->dcl->lnP           = 0.0;
  hit->dcl->is_reported   = FALSE;
  hit->dcl->is_included   = FALSE;

  pli->hit = hit;


  p7_trace_fs_Destroy(tr);
  p7_gmx_Destroy(gxppfs);
  return eslOK;

  ERROR:
    p7_trace_fs_Destroy(tr);
    p7_gmx_Destroy(gxppfs);
    p7_hit_Destroy(hit);
    return status;
}

int
path_exists (SPLICE_GRAPH *graph, int upstream_node, int downstream_node)
{

  int   i;
  int   current;
  int  *visited;
  int status;

  if(edge_exists(graph, upstream_node, downstream_node)) return TRUE;

  visited = NULL;
  ESL_ALLOC(visited, graph->num_nodes * sizeof(int));
  esl_vec_ISet(visited,   graph->num_nodes, 0);

  current = upstream_node;
  for (i = 0; i < graph->num_nodes; i++) {
    if(!visited[i]) {
      if(edge_exists(graph, current, i)) {
        if(path_finder(graph, i, downstream_node, visited)) {
          if ( visited != NULL) free(visited);
          return TRUE;
        }
      }
    }
  }

  if ( visited != NULL) free(visited);

  return FALSE;
  
  ERROR:
    if ( visited != NULL) free(visited);
    return status;
}

int 
path_finder (SPLICE_GRAPH *graph, int upstream_node, int downstream_node, int *visited) {

  int   i;
  int   current;

  if(edge_exists(graph, upstream_node, downstream_node)) return TRUE;

  current = upstream_node;

  for (i = 0; i < graph->num_nodes; i++) {
    if(!visited[i]) {
      if(edge_exists(graph, current, i)) {
        if(path_finder(graph, i, downstream_node, visited)) return TRUE;
      }
    }
  }
  
  return FALSE;
}

SPLICE_GAP*
find_the_gap (P7_TOPHITS *th, const P7_PROFILE *gm, ESL_SQ *target_seq, int orig_N, int upstream_node, int downstream_node, int revcomp) 
{

  int    hmm_ext;
  P7_HIT *ds_hit;
  P7_HIT *us_hit;
  SPLICE_GAP *gap;
  int     status; 

  gap = NULL;
  ESL_ALLOC(gap, sizeof(SPLICE_GAP));

  /* Looking in between exons */
  us_hit = th->hit[upstream_node];
  ds_hit = th->hit[downstream_node];

  if(us_hit->dcl->jhmm < ds_hit->dcl->ihmm)
    hmm_ext = 2;
  else
    hmm_ext = (us_hit->dcl->jhmm - ds_hit->dcl->ihmm) / 2 + 2;
  
  gap->hmm_start = ESL_MAX(us_hit->dcl->jhmm-hmm_ext, 1);
  gap->hmm_end   = ESL_MIN(ds_hit->dcl->ihmm+hmm_ext, gm->M);
  
  if(revcomp) {
    gap->seq_start = ESL_MIN(us_hit->dcl->jali+(hmm_ext*3+2), target_seq->start);
    gap->seq_end   = ESL_MAX(ds_hit->dcl->iali-(hmm_ext*3),   target_seq->end);
  }
  else {
    gap->seq_start = ESL_MAX(us_hit->dcl->jali-(hmm_ext*3+2), target_seq->start);
    gap->seq_end   = ESL_MIN(ds_hit->dcl->iali+(hmm_ext*3), target_seq->end);
  }
  
  gap->upstream_node   = upstream_node;
  gap->downstream_node = downstream_node;
  
  return gap;

  ERROR:
    if(gap != NULL) free(gap);
	return NULL;

}

P7_HIT**
align_the_gap(TARGET_RANGE *target_range, const P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, ESL_SQ *target_seq, ESL_GENCODE *gcode, SPLICE_GAP *gap, int *num_hits, int revcomp)
{
 
  int         i, j, h, z;
  int         z1,z2;
  int         orf_min;
  int         seq_sq_len;
  int         num_saved;
  int         low_idx;
  int         orf_seq_min, orf_seq_max;
  int         hit_seq_min, hit_seq_max;
  int         duplicate;
  float       vitsc;
  float       low_savesc;
  char         *spoofline;
  P7_HMM       *sub_hmm;
  P7_PROFILE   *sub_model;
  P7_OPROFILE  *sub_omodel;
  P7_TOPHITS   *th;
  P7_GMX       *vit_mx;
  P7_HIT       *orf_hit;
  P7_HIT      **top_ten; 
  P7_SCOREDATA *scoredata;
  ESL_DSQ      *sub_dsq;
  ESL_SQ       *sub_sq; 
  ESL_SQ       *orf_sq;
  ESL_GETOPTS  *go;
  ESL_GENCODE_WORKSTATE *wrk;
  int status;

  seq_sq_len =  abs(gap->seq_end - gap->seq_start) + 1;
  if(seq_sq_len < 3) return NULL; 

  th = target_range-> th;
 
  top_ten = NULL;
  ESL_ALLOC(top_ten, sizeof(P7_HIT*) * 10);
  
  sub_hmm    = extract_sub_hmm(hmm, gap->hmm_start, gap->hmm_end);  
  sub_model  = p7_profile_Create (sub_hmm->M, sub_hmm->abc); 
  sub_omodel = p7_oprofile_Create(sub_hmm->M, sub_hmm->abc );
  p7_ProfileConfig(sub_hmm, bg, sub_model, 100, p7_LOCAL);
  p7_oprofile_Convert(sub_model, sub_omodel);

  scoredata  = p7_hmm_ScoreDataCreate(sub_omodel, NULL);
  p7_hmm_ScoreDataComputeRest(sub_omodel, scoredata);
  
  sub_dsq = extract_sub_seq(target_seq, gap->seq_start, gap->seq_end, revcomp);
  sub_sq = esl_sq_CreateDigitalFrom(target_seq->abc, target_seq->name, sub_dsq, seq_sq_len, NULL, NULL, NULL); 
  sub_sq->start = gap->seq_start;
  sub_sq->end   = gap->seq_end; 
  
  orf_min = gap->hmm_end - gap->hmm_start + 1;
  orf_min = ESL_MIN(orf_min-12, 5);

  /* Create spoof getopts to hold flags checked by translation routines */
  spoofline = NULL;
  ESL_ALLOC(spoofline, sizeof(char)* 50); 
  sprintf(spoofline, "getopts -l %d", orf_min);
  go = esl_getopts_Create(Translation_Options);
  if (esl_opt_ProcessSpoof(go, spoofline) != eslOK) esl_fatal("errmsg");
 
  /* Create gencode workstate for use in translting orfs */
  wrk = esl_gencode_WorkstateCreate(go, gcode);
  wrk->orf_block = esl_sq_CreateDigitalBlock(1000, sub_model->abc);

  /* Translate orfs */
  esl_gencode_ProcessStart(gcode, wrk, sub_sq);
  esl_gencode_ProcessPiece(gcode, wrk, sub_sq);
  esl_gencode_ProcessEnd(wrk, sub_sq); 

  num_saved  = 0;
  low_savesc = 0; 
  vit_mx = p7_gmx_Create(sub_model->M,1024);

  for (i = 0; i < wrk->orf_block->count; i++)
  {

    orf_sq = &(wrk->orf_block->list[i]);
    duplicate = FALSE;

    if(revcomp) {
      orf_seq_min = sub_sq->start - sub_sq->n + orf_sq->end;
      orf_seq_max = sub_sq->start - sub_sq->n + orf_sq->start;
    }
    else {
      orf_seq_min = orf_sq->start + sub_sq->start - 1;
      orf_seq_max = orf_sq->end   + sub_sq->start - 1;
    }
    
    for (h = 0 ; h < target_range->orig_N; h++) {
      hit_seq_min = ESL_MIN(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
      hit_seq_max = ESL_MAX(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali); 
      if( hit_seq_min <= orf_seq_min && hit_seq_max >= orf_seq_max) {
        if(th->hit[h]->dcl->ihmm >= gap->hmm_start && th->hit[h]->dcl->jhmm <= gap->hmm_end) {
          duplicate = TRUE;    
          break; 
        }
      }
    }      

    if(duplicate) continue;
    
    p7_ReconfigUnihit(sub_model, orf_sq->n);
    p7_gmx_GrowTo(vit_mx, sub_model->M, orf_sq->n);
    p7_GViterbi(orf_sq->dsq, orf_sq->n, sub_model, vit_mx, &vitsc);

    orf_hit          = p7_hit_Create_empty();
    orf_hit->dcl     = p7_domain_Create_empty();
    orf_hit->dcl->tr = p7_trace_Create();
    
    p7_GTrace(orf_sq->dsq, orf_sq->n, sub_model, vit_mx, orf_hit->dcl->tr);
    p7_trace_Index(orf_hit->dcl->tr);

    for (z1 = orf_hit->dcl->tr->tfrom[0]; z1 < orf_hit->dcl->tr->N; z1++) if (orf_hit->dcl->tr->st[z1] == p7T_M) break;
    for (z2 = orf_hit->dcl->tr->tto[0]; z2 >= 0 ; z2--)                   if (orf_hit->dcl->tr->st[z2] == p7T_M) break;

    ESL_ALLOC(orf_hit->dcl->tr->c,  sizeof(int)  * orf_hit->dcl->tr->nalloc); 
	for (z = z1; z <= z2; z++) {
      if (orf_hit->dcl->tr->st[z] == p7T_M)
        orf_hit->dcl->tr->c[z] = 3;
      else
        orf_hit->dcl->tr->c[z] = 0;
    }

    orf_hit->dcl->ihmm =  orf_hit->dcl->tr->k[z1] + gap->hmm_start - 1;
    orf_hit->dcl->jhmm =  orf_hit->dcl->tr->k[z2] + gap->hmm_start - 1; 
    
    if(revcomp) {
      orf_hit->dcl->iali = sub_sq->start - sub_sq->n + orf_sq->start - 3*(orf_hit->dcl->tr->i[z1]-1);
      orf_hit->dcl->jali = sub_sq->start - sub_sq->n + orf_sq->start - 3*(orf_hit->dcl->tr->i[z2])+1;
    }
    else {
      orf_hit->dcl->iali = sub_sq->start + orf_sq->start + (orf_hit->dcl->tr->i[z1]*3-2) - 2;
      orf_hit->dcl->jali = sub_sq->start + orf_sq->start + (orf_hit->dcl->tr->i[z2]*3)   - 2;
    }



    orf_hit->dcl->ad = p7_alidisplay_Create(orf_hit->dcl->tr, 0, sub_omodel, orf_sq, NULL);

    p7_splice_compute_ali_scores(orf_hit->dcl, orf_hit->dcl->tr, orf_sq->dsq, sub_model, gm->abc->Kp);

    orf_seq_min = ESL_MIN(orf_hit->dcl->iali, orf_hit->dcl->jali);
    orf_seq_max = ESL_MAX(orf_hit->dcl->iali, orf_hit->dcl->jali);

    for (h = 0 ; h < target_range->orig_N; h++) {
     hit_seq_min = ESL_MIN(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
     hit_seq_max = ESL_MAX(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
     if( hit_seq_min <= orf_seq_min && hit_seq_max >= orf_seq_max) {
       if(orf_hit->dcl->ihmm >= th->hit[h]->dcl->ihmm && orf_hit->dcl->jhmm <= th->hit[h]->dcl->jhmm) {
         orf_hit->dcl->aliscore = -eslINFINITY;
         break;
       }
     }
    }
  
    if(orf_hit->dcl->aliscore > low_savesc) {
        
      for (z = z1; z <= z2; z++)
        orf_hit->dcl->tr->k[z] += gap->hmm_start - 1;

      if(num_saved < 10) {
        top_ten[num_saved] = orf_hit;
      }
      else {
        low_savesc =  top_ten[0]->dcl->aliscore;
        low_idx = 0;
        for(j = 1; j < 10; j++) {
          if(top_ten[j]->dcl->aliscore < low_savesc) {
            low_savesc  = top_ten[j]->dcl->aliscore;
            low_idx = j;
          }
        }
        p7_alidisplay_Destroy(top_ten[low_idx]->dcl->ad);
        free(top_ten[low_idx]->dcl->tr->c);
        p7_trace_Destroy(top_ten[low_idx]->dcl->tr);
        free(top_ten[low_idx]->dcl->scores_per_pos);
        p7_hit_Destroy(top_ten[low_idx]); 
       
        top_ten[low_idx] = orf_hit;   
      }
      num_saved++;
      if(low_savesc == -eslINFINITY || orf_hit->dcl->aliscore < low_savesc) low_savesc = orf_hit->dcl->aliscore;
    } 
    else {
      p7_alidisplay_Destroy(orf_hit->dcl->ad);
      free(orf_hit->dcl->tr->c);
      p7_trace_Destroy(orf_hit->dcl->tr);
      free(orf_hit->dcl->scores_per_pos);
      p7_hit_Destroy(orf_hit);
      
    }
  }

  *num_hits = ESL_MIN(num_saved, 10);
 
  p7_hmm_Destroy(sub_hmm);
  p7_profile_Destroy(sub_model);
  p7_oprofile_Destroy(sub_omodel); 
  p7_hmm_ScoreDataDestroy(scoredata);
  esl_sq_Destroy(sub_sq);
  esl_getopts_Destroy(go);
  esl_gencode_WorkstateDestroy(wrk);   
  p7_gmx_Destroy(vit_mx);
  free(spoofline);
  free(sub_dsq);

  return top_ten;


  ERROR:
   if(spoofline != NULL) free(spoofline);
   if(top_ten != NULL) free(top_ten);
   return NULL;

}

P7_HIT**
align_the_gap2(TARGET_RANGE *target_range, SPLICE_GAP *gap, const P7_HMM *hmm, const P7_BG *bg, ESL_SQ *target_seq, const ESL_GENCODE *gcode, int revcomp, int *num_hits, int remove_dups)
{
 
  int         i, h, y, z;
  int         z1, z2;
  int         c;
  int         gap_len;
  int         intron_cnt;
  int         hit_cnt;
  int         start_new;
  int         ihmm, jhmm;
  int         iali, jali;
  int         duplicate;
  int         exon_min, exon_max;
  int         hit_min, hit_max;
  P7_HMM       *sub_hmm;
  P7_FS_PROFILE *sub_fs_model;
  P7_GMX       *vit_mx;
  P7_HIT       *new_hit;
  P7_HIT      **ret_hits; 
  P7_TRACE     *tr; 
  P7_TOPHITS   *th;
  ESL_DSQ      *sub_dsq;
  ESL_SQ       *sub_sq; 
  int status;

  gap_len =  abs(gap->seq_end - gap->seq_start) + 1;
  if(gap_len < 3) return NULL; 

  th = target_range->th;

  /* Build sub model */ 
  sub_hmm    = extract_sub_hmm(hmm, gap->hmm_start, gap->hmm_end);  
  /* Turn off framshifts */
  sub_hmm->fs   = 0.;
  sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
  p7_ProfileConfig_fs(sub_hmm, bg, gcode, sub_fs_model, gap_len*3, p7_GLOBAL); //len*3 to get single nucleotide loops

  /* Get sub sequence */
  sub_dsq = extract_sub_seq(target_seq, gap->seq_start, gap->seq_end, revcomp);
  sub_sq = esl_sq_CreateDigitalFrom(target_seq->abc, target_seq->name, sub_dsq, gap_len, NULL, NULL, NULL); 
  sub_sq->start = gap->seq_start;
  sub_sq->end   = gap->seq_end; 
  vit_mx = p7_gmx_fs_Create(sub_fs_model->M, gap_len, gap_len, 2);
  tr = p7_trace_fs_Create();
  
  p7_sp_trans_semiglobal_Viterbi(sub_dsq, gcode, gap_len, sub_fs_model, vit_mx);
  //if(gap->seq_start == 13153881 ) p7_gmx_sp_Dump(stdout, vit_mx, p7_DEFAULT);
  p7_sp_trans_semiglobal_VTrace(sub_dsq, gap_len, gcode, sub_fs_model, vit_mx, tr); 
//  p7_trace_fs_Dump(stdout, tr, NULL, NULL); 
 
  /* Find number of introns in trace */
  intron_cnt = 0; 
  for(z = 0; z < tr->N; z++)
    if(tr->st[z] == p7T_R) intron_cnt++;

  ret_hits = NULL;
  ESL_ALLOC(ret_hits, sizeof(P7_HIT*) * (intron_cnt+1));
 
  /* Find first M state - start of first hit */
  for(z1 = 0; z1 < tr->N; z1++) if(tr->st[z1] == p7T_M) break;
  
  /* Find last M state state - end of last hit */
  for(z2 = tr->N-1; z1 >= 0; z2--) if(tr->st[z2] == p7T_M) break;

  hit_cnt = 0;
  start_new = TRUE;
  z = z1;
  while(z < z2) {

    if(start_new) {
 
      /* Save z value - currently set to fist M state in exon */
      y = z;

      /*Find end of exon */
      while(tr->st[z] != p7T_R && tr->st[z] != p7T_E) z++;
     
     /*Trace back to last M state of exon*/
      while(tr->st[z] != p7T_M) z--;

      /* Get exon coords */
      ihmm = tr->k[y] + gap->hmm_start - 1;  
      jhmm = tr->k[z] + gap->hmm_start - 1;
  
      if(revcomp) {
        iali = sub_sq->n - tr->i[y] + sub_sq->end + 2;
        jali = sub_sq->n - tr->i[z] + sub_sq->end;
      }
      else {
        iali = sub_sq->start + tr->i[y] - 3;
        jali = sub_sq->start + tr->i[z] - 1;
      }
       
      if(remove_dups) {
        /* check if this exon is a duplicate of an original hit */
        duplicate = FALSE;
  
        exon_min = ESL_MIN(iali, jali);
        exon_max = ESL_MAX(iali, jali);
  
        for (h = 0 ; h < target_range->orig_N; h++) {
          /* If hit is not in gap skip it */
          if(th->hit[h]->dcl->jhmm < gap->hmm_start || th->hit[h]->dcl->ihmm > gap->hmm_end) continue;
  
          hit_min = ESL_MIN(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
          hit_max = ESL_MAX(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);       
          /* Is exon entirley overlaped by hit on the sequence */ 
          if( hit_min <= exon_min && hit_max >= exon_max) {
            duplicate = TRUE;
            break;
          }
        }  
  
        if(duplicate) { z++; start_new = FALSE; continue; }
      }
      /* Create new hit and  set ihmm and iali coords*/
      new_hit          = p7_hit_Create_empty();
      
      new_hit->dcl     = p7_domain_Create_empty();
      new_hit->dcl->tr = p7_trace_fs_Create();
    
      new_hit->dcl->ihmm = ihmm;
      new_hit->dcl->jhmm = jhmm;
      new_hit->dcl->iali = iali;
      new_hit->dcl->jali = jali;
   
      /* Append starting special states */
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_S , 0, 0, 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_N , 0, tr->i[y]-3, 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_B , 0, tr->i[y]-3, 0);
      
      /* Append all states between first and last M state */
      for(i = y; i <= z; i++) {
        c = (tr->st[i] == p7T_M ? 3 : 0);
        p7_trace_fs_Append(new_hit->dcl->tr, tr->st[i] , tr->k[i], tr->i[i], c);
      }

      /* Append ending special states */
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_E, tr->k[z], tr->i[z], 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_C, 0, tr->i[z], 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_T , 0, 0, 0);
       //p7_trace_fs_Dump(stdout, new_hit->dcl->tr, NULL, NULL); 
      new_hit->dcl->ad = p7_alidisplay_fs_Create(new_hit->dcl->tr, 0, sub_fs_model, sub_sq, gcode);
      new_hit->dcl->scores_per_pos = NULL;
      p7_splice_compute_ali_scores_fs(new_hit->dcl, new_hit->dcl->tr, sub_dsq, sub_fs_model, sub_sq->abc);

      ret_hits[hit_cnt] = new_hit;
      hit_cnt++;
 
      start_new = FALSE;
    }
 
    z++;
    if(tr->st[z] == p7T_M) start_new = TRUE;     
     
  }
 
  *num_hits = hit_cnt;
 
  p7_hmm_Destroy(sub_hmm);
  p7_profile_fs_Destroy(sub_fs_model);
  esl_sq_Destroy(sub_sq);
  p7_gmx_Destroy(vit_mx);
  p7_trace_fs_Destroy(tr);
  free(sub_dsq);

  return ret_hits;


  ERROR:
   if(ret_hits != NULL) free(ret_hits);
   return NULL;

}


int
p7_splice_compute_ali_scores(P7_DOMAIN *dom, P7_TRACE *tr, ESL_DSQ *amino_dsq, const P7_PROFILE *gm, int K)
{

  int status;
  int i, k, n;
  int z1, z2;
  int N;

  if(tr->ndom == 0) p7_trace_Index(tr);
  
  N = tr->tto[0] - tr->tfrom[0] - 1;
  ESL_ALLOC( dom->scores_per_pos, sizeof(float) * N);
  for (n=0; n<N; n++)  dom->scores_per_pos[n] = 0.0;  

  i = tr->sqfrom[0];
  k = tr->hmmfrom[0];
  z1 = tr->tfrom[0]+1;
  z2 = tr->tto[0];
  n = 0;

  while (z1 < z2) {
    if (tr->st[z1] == p7T_M) {
      dom->scores_per_pos[n] = p7P_MSC(gm, k, amino_dsq[i]);
      if (tr->st[z1-1] == p7T_I) 
        dom->scores_per_pos[n] += p7P_TSC(gm, k-1, p7P_IM);
      else if (tr->st[z1-1] == p7T_D)
        dom->scores_per_pos[n] += p7P_TSC(gm, k-1, p7P_DM); 
      i++; k++; z1++; n++;
 
      while(z1 < z2 && tr->st[z1] == p7T_M) {
        dom->scores_per_pos[n] =  p7P_MSC(gm, k,   amino_dsq[i]);
        dom->scores_per_pos[n] += p7P_TSC(gm, k-1, p7P_MM);
        i++; k++; z1++; n++;
      }

    }
    else if (tr->st[z1] == p7T_I) {
      dom->scores_per_pos[n] = p7P_TSC(gm, k, p7P_MI);
      i++; z1++; n++;
      while (z1 < z2 && tr->st[z1] == p7T_I) {
        dom->scores_per_pos[n] = p7P_TSC(gm, k, p7P_II);
        i++; z1++; n++;
      }
     }
     else if (tr->st[z1] == p7T_D) {
       dom->scores_per_pos[n] = p7P_TSC(gm, k-1, p7P_MD);
       k++; z1++; n++;
       while (z1 < z2 && tr->st[z1] == p7T_D)  {
         dom->scores_per_pos[n] = p7P_TSC(gm, k-1, p7P_DD);
         k++; z1++; n++;
       } 
     }
     else ESL_XEXCEPTION(eslFAIL, "Impossible state from compute_ali_scores()");
  }

  dom->aliscore = 0.0;
  for (n=0; n<N; n++) dom->aliscore += dom->scores_per_pos[n];

  return eslOK;

  ERROR:
    return status;
}

int
p7_splice_compute_ali_scores_fs(P7_DOMAIN *dom, P7_TRACE *tr, ESL_DSQ *nuc_dsq, P7_FS_PROFILE *gm_fs, const ESL_ALPHABET *abc)
{

  int status;
  int i, k, c, n;
  int codon_idx;
  int z1, z2;
  int N;
 
  if(tr->ndom == 0) p7_trace_Index(tr);

  N = tr->tto[0] - tr->tfrom[0] - 1;
  ESL_ALLOC( dom->scores_per_pos, sizeof(float) * N);
  for (n=0; n<N; n++)  dom->scores_per_pos[n] = 0.0; 

  k  = tr->k[tr->tfrom[0]+1];
  z1 = tr->tfrom[0]+1;
  z2 = tr->tto[0];
  n  = 0;

  while (z1<z2) {
    
    if (tr->st[z1] == p7T_M) {
      i = tr->i[z1];
      c = tr->c[z1];
      if(c == 1) {
        if(esl_abc_XIsCanonical(abc, nuc_dsq[i]))
          codon_idx = p7P_CODON1(nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_QC2;
      }
      else if(c == 2) {
        if(esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i]))
          codon_idx = p7P_CODON2(nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_QC1;
      }
      else if(c == 3) {
        if(esl_abc_XIsCanonical(abc, nuc_dsq[i-2]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i]))
          codon_idx = p7P_CODON3(nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_C;
      }  
      else if(c == 4) {
        if(esl_abc_XIsCanonical(abc, nuc_dsq[i-3]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-2]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i]))
          codon_idx = p7P_CODON4(nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_QC1;
      }
      else if(c == 5) {
        if(esl_abc_XIsCanonical(abc, nuc_dsq[i-4]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-3]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-2]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i]))
          codon_idx = p7P_CODON5(nuc_dsq[i-4], nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_QC2;
      }
      
      dom->scores_per_pos[n] = p7P_MSC_CODON(gm_fs, k, codon_idx);
      if (tr->st[z1-1] == p7T_I)
        dom->scores_per_pos[n] += p7P_TSC(gm_fs, k-1, p7P_IM);
      else if (tr->st[z1-1] == p7T_D)
        dom->scores_per_pos[n] += p7P_TSC(gm_fs, k-1, p7P_DM);
      k++; z1++; n++;

      while(z1 < z2 && tr->st[z1] == p7T_M) {
        c = tr->c[z1];
        i = tr->i[z1];       
    
        if(c == 1) {
          if(esl_abc_XIsCanonical(abc, nuc_dsq[i]))
            codon_idx = p7P_CODON1(nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_QC2;
        }
        else if(c == 2) {
          if(esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i]))
            codon_idx = p7P_CODON2(nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_QC1;
        }
        else if(c == 3) {
          if(esl_abc_XIsCanonical(abc, nuc_dsq[i-2]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i]))
            codon_idx = p7P_CODON3(nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_C;
        }
        else if(c == 4) {
          if(esl_abc_XIsCanonical(abc, nuc_dsq[i-3]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i-2]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i]))
            codon_idx = p7P_CODON4(nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_QC1;
        }
        else if(c == 5) {
          if(esl_abc_XIsCanonical(abc, nuc_dsq[i-4]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i-3]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i-2]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i]))
            codon_idx = p7P_CODON5(nuc_dsq[i-4], nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_QC2;
        }
        
        dom->scores_per_pos[n] = p7P_MSC_CODON(gm_fs, k, codon_idx);
        dom->scores_per_pos[n] += p7P_TSC(gm_fs, k-1, p7P_MM);
        k++; z1++; n++;
      }
    }
    else if (tr->st[z1] == p7T_I) {
      dom->scores_per_pos[n] = p7P_TSC(gm_fs, k, p7P_MI);
      z1++; n++;
      while (z1 < z2 && tr->st[z1] == p7T_I) {
        dom->scores_per_pos[n] = p7P_TSC(gm_fs, k, p7P_II);
        z1++; n++;
      }
    }
    else if (tr->st[z1] == p7T_D) {
      dom->scores_per_pos[n] = p7P_TSC(gm_fs, k-1, p7P_MD);
      k++; z1++; n++;
      while (z1 < z2 && tr->st[z1] == p7T_D)  {
        dom->scores_per_pos[n] = p7P_TSC(gm_fs, k-1, p7P_DD);
        k++; z1++; n++;
      }
    }
    else ESL_XEXCEPTION(eslFAIL, "Impossible state from compute_ali_scores()");
  }

  dom->aliscore = 0.0;
  for (n=0; n<N; n++)  dom->aliscore += dom->scores_per_pos[n];  
 
  return eslOK;

  ERROR:
    return status;
}


P7_HMM*
extract_sub_hmm (const P7_HMM *hmm, int start, int end) 
{

  int i, k;
  int M;
  P7_HMM *sub_hmm;
  int status;

  M = end - start + 1; 
  sub_hmm = NULL;
  sub_hmm = p7_hmm_CreateShell();
  sub_hmm->flags = hmm->flags;  
  
  p7_hmm_CreateBody(sub_hmm, M, hmm->abc);

  if ((status = esl_strdup(hmm->name,   -1, &(sub_hmm->name)))   != eslOK) goto ERROR;
  if ((status = esl_strdup(hmm->acc,    -1, &(sub_hmm->acc)))    != eslOK) goto ERROR;
  if ((status = esl_strdup(hmm->desc,   -1, &(sub_hmm->desc)))   != eslOK) goto ERROR;

  sub_hmm->nseq       = hmm->nseq;
  sub_hmm->eff_nseq   = hmm->eff_nseq;
  sub_hmm->max_length = hmm->max_length;
  sub_hmm->checksum   = hmm->checksum;
  sub_hmm->offset     = hmm->offset;
  sub_hmm->fs         = hmm->fs;

  for (i = 0; i < p7_NEVPARAM; i++) sub_hmm->evparam[i] = hmm->evparam[i];
  for (i = 0; i < p7_NCUTOFFS; i++) sub_hmm->cutoff[i]  = hmm->cutoff[i];
  for (i = 0; i < p7_MAXABET;  i++) sub_hmm->compo[i]   = hmm->compo[i];


  /*Copy zeroth positions */  
  esl_vec_FCopy(hmm->t[0],   p7H_NTRANSITIONS, sub_hmm->t[0]);
  esl_vec_FCopy(hmm->mat[0], hmm->abc->K,      sub_hmm->mat[0]);
  esl_vec_FCopy(hmm->ins[0], hmm->abc->K,      sub_hmm->ins[0]);

  if (hmm->flags & p7H_RF)    sub_hmm->rf[0]        = hmm->rf[0];  
  if (hmm->flags & p7H_MMASK) sub_hmm->mm[0]        = hmm->mm[0];
  if (hmm->flags & p7H_CONS)  sub_hmm->consensus[0] = hmm->consensus[0];
  if (hmm->flags & p7H_CS)    sub_hmm->cs[0]        = hmm->cs[0];

  k = 1;
  for (i = start; i <= end; i++) {
    esl_vec_FCopy(hmm->t[i],   p7H_NTRANSITIONS, sub_hmm->t[k]);
    esl_vec_FCopy(hmm->mat[i], hmm->abc->K,      sub_hmm->mat[k]);
    esl_vec_FCopy(hmm->ins[i], hmm->abc->K,      sub_hmm->ins[k]); 
    if (hmm->flags & p7H_RF)    sub_hmm->rf[k]        = hmm->rf[i];
    if (hmm->flags & p7H_MMASK) sub_hmm->mm[k]        = hmm->mm[i];
    if (hmm->flags & p7H_CONS)  sub_hmm->consensus[k] = hmm->consensus[i];
    if (hmm->flags & p7H_CS)    sub_hmm->cs[k]        = hmm->cs[i];
    k++;
  }
  
  if (hmm->flags & p7H_RF)    sub_hmm->rf[M+1]        = hmm->rf[hmm->M+1]; 
  if (hmm->flags & p7H_MMASK) sub_hmm->mm[M+1]        = hmm->mm[hmm->M+1];
  if (hmm->flags & p7H_CONS)  sub_hmm->consensus[M+1] = hmm->consensus[hmm->M+1];
  if (hmm->flags & p7H_CS)    sub_hmm->cs[M+1]        = hmm->cs[hmm->M+1];

  /* No delete atart at postion M */
  sub_hmm->t[M][p7H_MD] = 0.0;
  sub_hmm->t[M][p7H_DD] = 0.0; 

  
  return sub_hmm;

  ERROR:
    if(sub_hmm != NULL) p7_hmm_Destroy(sub_hmm);
    return NULL;
}



ESL_DSQ*
extract_sub_seq(ESL_SQ *target_seq, int start, int end, int revcomp)
{

  int i,j;
  int L;
  int cp_start, cp_end;
  ESL_DSQ *sub_dsq;
  int status;

  L = abs(end - start) + 1;
  sub_dsq = NULL;
  ESL_ALLOC(sub_dsq, sizeof(ESL_DSQ) * (L+2));

  sub_dsq[0] = eslDSQ_SENTINEL;
 
  if(revcomp) {

    cp_start = target_seq->start - start + 1;
    cp_end   = target_seq->start - end   + 1;

    j = 1;
    for(i = cp_start; i <= cp_end; i++) {
      sub_dsq[j] = target_seq->dsq[i];
      j++;
    }
  }
  else { 
    
    cp_start = start - target_seq->start + 1;
    cp_end   = end   - target_seq->start + 1; 

    j = 1;
    for(i = cp_start; i <= cp_end; i++) {
      sub_dsq[j] = target_seq->dsq[i];
      j++;
    }  
    
  }

  return sub_dsq;

  ERROR:
   if(sub_dsq != NULL) free(sub_dsq);
   return NULL;
}

int
add_missing_node_to_graph(SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_HIT *hit, int M)
{

  int   up,down;
  SPLICE_NODE *new_node;
  SPLICE_NODE *tmp_node;
  int   status;

  /* Make sure we have room for our new node */
  if((status = splice_graph_grow(graph)) != eslOK) goto ERROR;

  graph->hit_scores[graph->num_nodes]     = hit->dcl->aliscore; //This contiains the summed alisc
  graph->path_scores[graph->num_nodes] = -eslINFINITY;
 
  graph->out_edge_cnt[graph->num_nodes]  = 0;
  graph->in_edge_cnt[graph->num_nodes]   = 0;
  graph->best_out_edge[graph->num_nodes] = -1;
  graph->best_in_edge[graph->num_nodes]  = -1;

  graph->ds_nodes[graph->num_nodes] = NULL;
  graph->edges[graph->num_nodes]    = NULL;  

  /*Set new node as downstream of all existing upstream nodes */
  for(up  = 0; up < graph->num_nodes; up++) {
     
    if(th->hit[up]->dcl->ihmm >= th->hit[graph->num_nodes]->dcl->ihmm &&
       th->hit[up]->dcl->jhmm >= th->hit[graph->num_nodes]->dcl->jhmm)
        continue; 

    if (( graph->revcomp  && th->hit[graph->num_nodes]->dcl->jali > th->hit[up]->dcl->iali) ||
         ((!graph->revcomp) && th->hit[graph->num_nodes]->dcl->jali < th->hit[up]->dcl->iali))
         continue;

    ESL_ALLOC(new_node, sizeof(SPLICE_NODE));
    new_node->node_id = graph->num_nodes;
    new_node->gap_checked = FALSE;
    new_node->edge_exists = FALSE;
    new_node->next = NULL;   
 
    tmp_node = graph->ds_nodes[up];
    if(tmp_node == NULL)
      graph->ds_nodes[up] = new_node;
    else {
      while(tmp_node->next != NULL)
        tmp_node = tmp_node->next;

      tmp_node->next = new_node;
    }
  }
  
  /* Set existing nodes that are downstream of new node as downstream */
  for(down  = 0; down < graph->num_nodes; down++) {

    if(th->hit[graph->num_nodes]->dcl->ihmm >= th->hit[down]->dcl->ihmm &&
       th->hit[graph->num_nodes]->dcl->jhmm >= th->hit[down]->dcl->jhmm)
      continue;

    if (( graph->revcomp  && th->hit[down]->dcl->jali > th->hit[graph->num_nodes]->dcl->iali) ||
       ((!graph->revcomp) && th->hit[down]->dcl->jali < th->hit[graph->num_nodes]->dcl->iali))
       continue;

    ESL_ALLOC(new_node, sizeof(SPLICE_NODE));
    new_node->node_id = down;
    new_node->gap_checked = FALSE;
    new_node->edge_exists = FALSE;
    new_node->next = NULL;

    tmp_node = graph->ds_nodes[up];
    if(tmp_node == NULL)
      graph->ds_nodes[up] = new_node;
    else {
      while(tmp_node->next != NULL)
        tmp_node = tmp_node->next;

      tmp_node->next = new_node;
    }
  }
  graph->num_nodes++;

  return eslOK; 

  ERROR:
    return status;
}


int
add_missed_hit_to_target_range(TARGET_RANGE *target_range, P7_HIT *hit)
{

  P7_TOPHITS *th;
  int        status;

  th = target_range->th;

  if((status = target_range_grow(target_range)) != eslOK) goto ERROR;
 
  target_range->in_target_range[th->N] = TRUE; 
  th->hit[th->N] = hit;
  th->N++;

  return eslOK;

  ERROR:
    return status;
}

int
bridge_the_gap(TARGET_RANGE *target_range, SPLICE_GRAPH *graph, P7_HIT **top_ten, const P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, const ESL_GENCODE *gcode, const ESL_SQ *target_seq, int up_gap, int down_gap, int num_hits)
{
  int         i;
  int         up, down;
  int         curr_num_connected;
  int         next_num_connected;
  int         *curr_connected_hits;
  int         *next_connected_hits;
  int         *has_edge;
  int         *node_id;
  SPLICE_EDGE *edge;
  P7_TOPHITS  *th;
  int         status;

  th = target_range->th;

  curr_connected_hits = NULL;
  next_connected_hits = NULL;
  ESL_ALLOC(curr_connected_hits, sizeof(int) * num_hits); 
  ESL_ALLOC(next_connected_hits, sizeof(int) * num_hits);

  has_edge = NULL;
  ESL_ALLOC(has_edge, sizeof(int) * num_hits);

  node_id = NULL;
  ESL_ALLOC(node_id, sizeof(int) * num_hits);
 
  curr_num_connected = 0;
  next_num_connected = 0;
 
  for(i = 0; i < num_hits; i++) {
    curr_connected_hits[i] = FALSE;
    next_connected_hits[i] = FALSE;
    has_edge[i]            = FALSE;
    node_id[i]             = -1;
  }
 
  edge = NULL;
   
  /* Connect new hits upstream of original gap downstream hit */
  for(up = 0; up < num_hits; up++) {
    
    edge = connect_nodes_with_edges(top_ten[up], th->hit[down_gap], gm, hmm, bg, gcode, target_seq, graph->revcomp);
    if(edge != NULL) {
      
      edge->upstream_node_id   = th->N;
      edge->downstream_node_id = down_gap;
      edge->target_seq_start   = target_seq->start;
      edge->target_seq_end     = target_seq->end;
      edge->target_seq_n       = target_seq->n;
        
      if(!has_edge[up]) {
        node_id[up] = th->N;
        if ((status = add_missed_hit_to_target_range(target_range, top_ten[up])) != eslOK) goto ERROR;
        if ((status = add_missing_node_to_graph(graph, th, top_ten[up], gm->M)) != eslOK) goto ERROR;
      }             
      add_edge_to_graph(graph, edge);
      has_edge[up] = TRUE;
      curr_connected_hits[up] = TRUE;
      curr_num_connected++;
    }
  }

  /* Connect new hits to other new hits connected to original hits */
  for(up = 0; up < num_hits; up++) {
    for(down = 0; down < num_hits; down++) {
      if(!curr_connected_hits[down] || up == down) continue;
       
      /* check that up is upstream of down */
      if(top_ten[up]->dcl->ihmm >= top_ten[down]->dcl->ihmm &&
         top_ten[up]->dcl->jhmm >= top_ten[down]->dcl->jhmm) 
        continue;

      if (( graph->revcomp  && top_ten[down]->dcl->jali > top_ten[up]->dcl->iali) ||
         ((!graph->revcomp) && top_ten[down]->dcl->jali < top_ten[up]->dcl->iali)) 
        continue;
      
      edge = connect_nodes_with_edges(top_ten[up], top_ten[down], gm, hmm, bg, gcode, target_seq, graph->revcomp);

      if(edge != NULL) {
        edge->upstream_node_id   = (has_edge[up] ? node_id[up] : th->N);
        edge->downstream_node_id = node_id[down];
        edge->target_seq_start   = target_seq->start;
        edge->target_seq_end     = target_seq->end;
        edge->target_seq_n       = target_seq->n;

        if(!has_edge[up]) {
          node_id[up] = th->N;
          if ((status = add_missed_hit_to_target_range(target_range, top_ten[up])) != eslOK) goto ERROR;
          if ((status = add_missing_node_to_graph(graph, th, top_ten[up], gm->M)) != eslOK) goto ERROR;
        }
        add_edge_to_graph(graph, edge);
        has_edge[up] = TRUE;
        next_connected_hits[up] = TRUE;
        next_num_connected++;
      }
    }
  }
   
  free(curr_connected_hits);
  curr_connected_hits = next_connected_hits;
  curr_num_connected  = next_num_connected;
  next_connected_hits = NULL;
  ESL_ALLOC(next_connected_hits, sizeof(int) * num_hits);
  for(i = 0; i < num_hits; i++)
    next_connected_hits[i] = FALSE;
  next_num_connected = 0;

  /* Continue connecting new hit upstream until we run out of connections */ 
  while(curr_num_connected > 0) {
   
    for(up = 0; up < num_hits; up++) {
      for(down = 0; down < num_hits; down++) {
        if(!curr_connected_hits[down] || up == down) continue;

        /* check that up is upstream of down */
        if(top_ten[up]->dcl->ihmm >= top_ten[down]->dcl->ihmm &&
           top_ten[up]->dcl->jhmm >= top_ten[down]->dcl->jhmm)
          continue;

        if (( graph->revcomp  && top_ten[down]->dcl->jali > top_ten[up]->dcl->iali) ||
           ((!graph->revcomp) && top_ten[down]->dcl->jali < top_ten[up]->dcl->iali)) 
          continue;     
         
        edge = connect_nodes_with_edges(top_ten[up], top_ten[down], gm, hmm, bg, gcode, target_seq, graph->revcomp);
 
        if(edge != NULL) {
          edge->upstream_node_id   = (has_edge[up] ? node_id[up] : th->N);
          edge->downstream_node_id = node_id[down];
          
          edge->target_seq_start   = target_seq->start;
          edge->target_seq_end     = target_seq->end;
          edge->target_seq_n       = target_seq->n;

          if(!has_edge[up]) {
            node_id[up] = th->N;
            if ((status = add_missed_hit_to_target_range(target_range, top_ten[up])) != eslOK) goto ERROR;
            if ((status = add_missing_node_to_graph(graph, th, top_ten[up], gm->M)) != eslOK) goto ERROR;
          }
          add_edge_to_graph(graph, edge);
          has_edge[up] = TRUE;
          next_connected_hits[up] = TRUE;
          next_num_connected++;
        }
      }
    }
    free(curr_connected_hits);
    curr_connected_hits = next_connected_hits;
    curr_num_connected  = next_num_connected;
    next_connected_hits = NULL;
    ESL_ALLOC(next_connected_hits, sizeof(int) * num_hits);
    for(i = 0; i < num_hits; i++)
      next_connected_hits[i] = FALSE;
    next_num_connected = 0;
  } 

  for(i = 0; i < num_hits; i++) {
    curr_connected_hits[i] = FALSE;
    next_connected_hits[i] = FALSE;
  }

  /* Connect new hits downstream of original gap upstream hit * */
  for(down = 0; down < num_hits; down++) { 
    
    edge = connect_nodes_with_edges(th->hit[up_gap], top_ten[down], gm, hmm, bg, gcode, target_seq, graph->revcomp);
    if(edge != NULL) {
        
      edge->upstream_node_id   = up_gap;
      edge->downstream_node_id = (has_edge[down] ? node_id[down] : th->N);;
      edge->target_seq_start   = target_seq->start;
      edge->target_seq_end     = target_seq->end;
      edge->target_seq_n       = target_seq->n;
   
      if(!has_edge[down]) {
   
        node_id[down] = th->N;
         if ((status = add_missed_hit_to_target_range(target_range, top_ten[down])) != eslOK) goto ERROR;
         if ((status = add_missing_node_to_graph(graph, th, top_ten[down], gm->M)) != eslOK) goto ERROR;
      }
      
      add_edge_to_graph(graph, edge);
      has_edge[down] = TRUE;
      curr_connected_hits[down] = TRUE;
      curr_num_connected++;
    }
  }
  
 /* Connect new hits that connected to original hits upstream of any non-original hits */
  for(up = 0; up < num_hits; up++) {
    for(down = 0; down < num_hits; down++) {
      if(!curr_connected_hits[up] || up == down) continue;
   
      /* check that up is upstream of down */
      if(top_ten[up]->dcl->ihmm >= top_ten[down]->dcl->ihmm &&
         top_ten[up]->dcl->jhmm >= top_ten[down]->dcl->jhmm)
        continue;

      if (( graph->revcomp  && top_ten[down]->dcl->jali > top_ten[up]->dcl->iali) ||
         ((!graph->revcomp) && top_ten[down]->dcl->jali < top_ten[up]->dcl->iali))
        continue;
     
      if(edge_exists(graph, node_id[up], (has_edge[down] ? node_id[down] : th->N))) continue;
 
      edge = connect_nodes_with_edges(top_ten[up], top_ten[down], gm, hmm, bg, gcode, target_seq, graph->revcomp);    
      if(edge != NULL) {
        edge->upstream_node_id   = node_id[up];
        edge->downstream_node_id = (has_edge[down] ? node_id[down] : th->N);
        edge->target_seq_start   = target_seq->start;
        edge->target_seq_end     = target_seq->end;
        edge->target_seq_n       = target_seq->n;

        if(!has_edge[down]) {
          node_id[down] = th->N;
          if ((status = add_missed_hit_to_target_range(target_range, top_ten[down])) != eslOK) goto ERROR;
          if ((status = add_missing_node_to_graph(graph, th, top_ten[down], gm->M)) != eslOK) goto ERROR;
        }
        add_edge_to_graph(graph, edge);
        has_edge[down] = TRUE;
        curr_connected_hits[down] = TRUE;
        curr_num_connected++;
      }
    }
  }
   
  free(curr_connected_hits);
  curr_connected_hits = next_connected_hits;
  curr_num_connected  = next_num_connected;
  next_connected_hits = NULL;
  ESL_ALLOC(next_connected_hits, sizeof(int) * num_hits);
  for(i = 0; i < num_hits; i++)
    next_connected_hits[i] = FALSE;
  next_num_connected = 0;
  
  /* Continue connecting non-original hits downstream until we run out of connections */ 
  while(curr_num_connected > 0) {
    for(up = 0; up < num_hits; up++) {
      for(down = 0; down < num_hits; down++) {
        if(!curr_connected_hits[up] || up == down) continue;
     
        /* check that up is upstream of down */
        if(top_ten[up]->dcl->ihmm >= top_ten[down]->dcl->ihmm &&
           top_ten[up]->dcl->jhmm >= top_ten[down]->dcl->jhmm)
          continue;
  
        if (( graph->revcomp  && top_ten[down]->dcl->jali > top_ten[up]->dcl->iali) ||
           ((!graph->revcomp) && top_ten[down]->dcl->jali < top_ten[up]->dcl->iali))
          continue;

        if(edge_exists(graph, node_id[up], (has_edge[down] ? node_id[down] : th->N))) continue;
  
        edge = connect_nodes_with_edges(top_ten[up], top_ten[down], gm, hmm, bg, gcode, target_seq, graph->revcomp);    
        if(edge != NULL) {
          
          edge->upstream_node_id   = node_id[up];
          edge->downstream_node_id = (has_edge[down] ? node_id[down] : th->N);;
          edge->target_seq_start   = target_seq->start;
          edge->target_seq_end     = target_seq->end;
          edge->target_seq_n       = target_seq->n;

          if(!has_edge[down]) {
            node_id[down] = th->N;
            if ((status = add_missed_hit_to_target_range(target_range, top_ten[down])) != eslOK) goto ERROR;
            if ((status = add_missing_node_to_graph(graph, th, top_ten[down], gm->M)) != eslOK) goto ERROR;
          }
          add_edge_to_graph(graph, edge);
          has_edge[down] = TRUE;
          curr_connected_hits[down] = TRUE;
          curr_num_connected++;
        }
      }
    }    
    free(curr_connected_hits);
    curr_connected_hits = next_connected_hits;
    curr_num_connected  = next_num_connected;
    next_connected_hits = NULL;
    ESL_ALLOC(next_connected_hits, sizeof(int) * num_hits);
    for(i = 0; i < num_hits; i++) 
      next_connected_hits[i] = FALSE;
    next_num_connected = 0;
  }
  

  for (i = 0; i < num_hits; i++) {
    if(!has_edge[i]) {
      p7_alidisplay_Destroy(top_ten[i]->dcl->ad);
      free(top_ten[i]->dcl->tr->c);
      p7_trace_Destroy(top_ten[i]->dcl->tr);
      free(top_ten[i]->dcl->scores_per_pos);
      p7_hit_Destroy(top_ten[i]); 
    } 
  }   

  if(curr_connected_hits != NULL) free(curr_connected_hits);
  if(next_connected_hits != NULL) free(next_connected_hits);
  if(has_edge            != NULL) free(has_edge);
  if(node_id             != NULL) free(node_id);

  return eslOK; 

  ERROR:
    if(curr_connected_hits != NULL) free(curr_connected_hits);
    if(next_connected_hits != NULL) free(next_connected_hits);
    if(has_edge            != NULL) free(has_edge);
    if(node_id             != NULL) free(node_id); 
    return status;
}


void 
target_range_dump(FILE *fp, TARGET_RANGE *target_range, int print_hits) {

  int i;
  int hit_cnt;
  P7_HIT *hit = NULL;

  if (target_range == NULL) { fprintf(fp, " [ target range is NULL ]\n"); return; }

  hit_cnt = 0;
  for (i = 0; i < target_range->th->N; i++)
    if(target_range->in_target_range[i]) hit_cnt++;

  fprintf(fp, " Target Range Sequence Name      : %s\n", target_range->seqname);
  fprintf(fp, " Target Range Reverse Complement : %s\n", (target_range->revcomp ? "YES" : "NO") );
  fprintf(fp, " Target Range Hit Count          : %d\n", hit_cnt);
  fprintf(fp, "\n");

  if(print_hits) {
    fprintf(fp, "   %4s %6s %9s %12s %12s %10s %10s\n",
            "Hit", "hmm_from", "hmm_to", "seq_from", "seq_to", "score", "reportable");  

    for(i = 0; i < target_range->th->N; i++) {

      if(!target_range->in_target_range[i])
        continue;

      hit = target_range->th->hit[i];     

      if ( i < target_range->orig_N)
   	    fprintf(fp, "   %4d %6d %9d %12" PRId64 " %12" PRId64 " %10.2f %10s\n", 
                i+1, hit->dcl->ihmm, hit->dcl->jhmm, hit->dcl->iali, hit->dcl->jali, hit->dcl->aliscore, (target_range->reportable[i] ? "YES" : "NO"));	 
      else
        fprintf(fp, "   %4d  %6d  %9d %12" PRId64 " %12" PRId64 " %10.2f\n",
                i+1, hit->dcl->ihmm, hit->dcl->jhmm, hit->dcl->iali, hit->dcl->jali, hit->dcl->aliscore);
   }
  }
  
  fprintf(fp, "\n"); 

  return;
}


void
graph_dump(FILE *fp, TARGET_RANGE *target_range, SPLICE_GRAPH *graph, int print_edges) 
{


  int          i,j;
  int          nuc_end, nuc_start;
  int          num_nodes;
  int          in_range_cnt;
  float        edge_scores[graph->num_nodes];
  SPLICE_EDGE *tmp_edge;

  if (graph == NULL) { fprintf(fp, " [ graph is NULL ]\n"); return; }
  num_nodes = esl_vec_ISum(target_range->in_target_range, graph->num_nodes);  

  fprintf(fp, " SPLICE_GRAPH\n");
  fprintf(fp, " Number of Nodes  : %d\n", num_nodes); 
  fprintf(fp, " Number of Edges  : %d\n", graph->num_edges); 

  if(graph->num_edges > 0) {
    fprintf(fp, "\n Splice Score matrix \n");
    fprintf(fp, "     ");
    for(i = 0; i < num_nodes/2; i++ )
      fprintf(fp, "%8c", ' ');
    fprintf(fp, "  DOWN  \n");

    fprintf(fp, "     ");
    for(i = 0; i < graph->num_nodes; i++ ) {
      if(target_range->in_target_range[i])
        fprintf(fp, "%8d", i+1);
    }
    fprintf(fp, "\n");

    in_range_cnt = 0;
    for(i = 0; i < graph->num_nodes; i++ ) {

      if(!target_range->in_target_range[i]) continue;
      if (in_range_cnt == num_nodes/2)
        fprintf(fp, "UP");
      else
        fprintf(fp, "  ");
      in_range_cnt++;
      for(j = 0; j < graph->num_nodes; j++ ) 
        edge_scores[j] = -eslINFINITY;

      fprintf(fp, "%7d", i+1);


      tmp_edge = graph->edges[i];
      while(tmp_edge != NULL) {
        edge_scores[tmp_edge->downstream_node_id] = tmp_edge->splice_score;
        tmp_edge = tmp_edge->next;
      }

      for(j = 0; j < graph->num_nodes; j++ ) { 
        if(target_range->in_target_range[j])
          fprintf(fp, "%8.2f", edge_scores[j]);
      }
      fprintf(fp, "\n");
    }    
 }

  fprintf(fp, "\n Hit Scores  \n");
  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!target_range->in_target_range[i]) continue;
    fprintf(fp, "%8d", i+1);
  }
  fprintf(fp, "\n");

  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!target_range->in_target_range[i]) continue; 
    fprintf(fp, "%8.2f", graph->hit_scores[i]);
  }
  fprintf(fp, "\n");

  fprintf(fp, "\n Upstream Path Scores  \n");
  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!target_range->in_target_range[i]) continue;
    fprintf(fp, "%8d", i+1);
  }
  fprintf(fp, "\n");

  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!target_range->in_target_range[i]) continue;
    fprintf(fp, "%8.2f", graph->path_scores[i]);
  }
  fprintf(fp, "\n");

  fprintf(fp, "\n Downstream Path Scores  \n");
  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!target_range->in_target_range[i]) continue;
    fprintf(fp, "%8d", i+1);
  }
  fprintf(fp, "\n");

   fprintf(fp, "\n Best Edge  \n");
  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!target_range->in_target_range[i]) continue;
    fprintf(fp, "%8d", i+1);
  }
  fprintf(fp, "\n");

  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!target_range->in_target_range[i]) continue;
    fprintf(fp, "%8d", graph->best_out_edge[i]+1);
  }
  fprintf(fp, "\n");

 
  if (print_edges) {

    fprintf(fp, "\n Edge Data  \n\n");
    for(i = 0; i < graph->num_nodes; i++ ) {
      if(!target_range->in_target_range[i]) continue;
      tmp_edge = graph->edges[i];
      while(tmp_edge != NULL) {
        nuc_end =   tmp_edge->upstream_spliced_nuc_end;
        nuc_start = tmp_edge->downstream_spliced_nuc_start;
        fprintf(fp, "    Edge from Upstream Node %d to Downstream Node %d\n", tmp_edge->upstream_node_id+1, tmp_edge->downstream_node_id+1);
        fprintf(fp, "                                   %s   %s\n", "Amino", "Nucleotide");
        fprintf(fp, "      Upsteam Node End Coords:     %5d  %10d\n", tmp_edge->upstream_spliced_amino_end, nuc_end);
        fprintf(fp, "      Downsteam Node Start Coords: %5d  %10d\n", tmp_edge->downstream_spliced_amino_start, nuc_start);
        fprintf(fp, "\n");
        tmp_edge = tmp_edge->next;
      }
    }
  }
  fprintf(fp, "\n");
 

  return;

}

void
ds_graph_dump(FILE *fp, TARGET_RANGE *target_range, SPLICE_GRAPH *graph) 
{


  int          i,j;
  int          num_nodes;
  int          in_range_cnt;

  if (graph == NULL) { fprintf(fp, " [ graph is NULL ]\n"); return; }
  num_nodes = esl_vec_ISum(target_range->in_target_range, graph->num_nodes);  

  fprintf(fp, " DOWNSTREAM_GRAPH\n");
  fprintf(fp, " Number of Nodes  : %d\n", num_nodes); 

  fprintf(fp, "     ");
  for(i = 0; i < num_nodes/2; i++ )
    fprintf(fp, "%8c", ' ');
  fprintf(fp, "  DOWN  \n");

  fprintf(fp, "     ");
  for(i = 0; i < graph->num_nodes; i++ ) {
    if(target_range->in_target_range[i])
      fprintf(fp, "%8d", i+1);
  }
  fprintf(fp, "\n");

  in_range_cnt = 0;
  for(i = 0; i < graph->num_nodes; i++ ) {

    if(!target_range->in_target_range[i]) continue;
    if (in_range_cnt == num_nodes/2)
      fprintf(fp, "UP");
    else
      fprintf(fp, "  ");
    in_range_cnt++;

    fprintf(fp, "%7d", i+1);
      
    for(j = 0; j < graph->num_nodes; j++ ) { 
      if(!target_range->in_target_range[j]) continue;

      if(is_upstream(graph, i, j)) fprintf(fp, "%8s", "YES");
      else                         fprintf(fp, "%8s", "NO");
        
    }
    fprintf(fp, "\n");
  }    
 

  fprintf(fp, "\n Hit Scores  \n");
  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!target_range->in_target_range[i]) continue;
    fprintf(fp, "%8d", i+1);
  }
  fprintf(fp, "\n");

  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!target_range->in_target_range[i]) continue; 
    fprintf(fp, "%8.2f", graph->hit_scores[i]);
  }
  fprintf(fp, "\n");
  
  return;

}


void
path_dump(FILE *fp, SPLICE_PATH *path)
{
  
  int i;
  
  fprintf(fp, "  Path Length  %d\n", path->path_len);
  fprintf(fp, "  %4s %4s %9s %9s %10s %10s %9s %10s \n", "Step", "Node", "hmm_start", "hmm_end", "seq_start", "seq_end", "hit_score", "edge_score");
  for(i = 0; i < path->path_len; i++) {
    fprintf(fp, "  %4d %4d %9d %9d %10d %10d %9.2f %10.2f\n", i+1, path->node_id[i]+1,
      path->downstream_spliced_amino_start[i], path->upstream_spliced_amino_end[i+1], path->downstream_spliced_nuc_start[i], path->upstream_spliced_nuc_end[i+1], path->hit_scores[i], path->edge_scores[i]);
  }
  
  fprintf(fp, "\n");

  return;
}
