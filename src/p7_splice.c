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

  target_range->orig_hit_idx = NULL;
  ESL_ALLOC(target_range->orig_hit_idx, nalloc * sizeof(int*));

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
     ESL_REALLOC(th->hit, sizeof(P7_HIT *) * th->Nalloc);
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

  if(target_range->orig_hit_idx != NULL)
    free(target_range->orig_hit_idx); 

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

  edge->splice_score = -eslINFINITY;
  edge->signal_score = -eslINFINITY;

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

  graph->num_nodes      = 0;
  graph->orig_num_nodes = 0;
  graph->num_edges      = 0;

  graph->best_path_length = 0;
  graph->best_path_start  = 0;
  graph->best_path_end    = 0;

  graph->out_edge_cnt    = NULL; 
  graph->in_edge_cnt     = NULL;
  graph->orig_out_edge   = NULL;
  graph->orig_in_edge    = NULL;
  graph->best_out_edge   = NULL;
  graph->best_in_edge    = NULL;

  graph->is_upstream     = NULL;
  graph->is_upstream_mem = NULL;

  graph->edge_id         = NULL;
  graph->edge_id_mem     = NULL;

  graph->path_scores     = NULL;
  graph->hit_scores      = NULL;
  graph->edge_scores     = NULL;
  graph->edge_scores_mem = NULL;

  graph->edges           = NULL;

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

  if (graph == NULL) return;

  if(graph->is_upstream     != NULL) free(graph->is_upstream);
  if(graph->is_upstream_mem != NULL) free(graph->is_upstream_mem);
  if(graph->edge_id         != NULL) free(graph->edge_id);
  if(graph->edge_id_mem     != NULL) free(graph->edge_id_mem);
  if(graph->edge_scores     != NULL) free(graph->edge_scores);
  if(graph->edge_scores_mem != NULL) free(graph->edge_scores_mem);
  if(graph->path_scores     != NULL) free(graph->path_scores); 
  if(graph->hit_scores      != NULL) free(graph->hit_scores);

  if(graph->out_edge_cnt    != NULL) free(graph->out_edge_cnt);
  if(graph->in_edge_cnt     != NULL) free(graph->in_edge_cnt);
  if(graph->orig_out_edge   != NULL) free(graph->orig_out_edge);
  if(graph->orig_in_edge    != NULL) free(graph->orig_in_edge);
  if(graph->best_out_edge   != NULL) free(graph->best_out_edge);
  if(graph->best_in_edge    != NULL) free(graph->best_in_edge);
 
  for(i = 0; i < graph->num_edges; i++)
    free(graph->edges[i]);  

  if(graph->edges           != NULL) free(graph->edges);
  
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

  graph->nalloc         = num_nodes*2;
  graph->num_nodes      = num_nodes;
  graph->orig_num_nodes = num_nodes;

  /* Allocate adjacency matrix for splice scores */ 
  ESL_ALLOC(graph->edge_scores,       sizeof(float*)       * graph->nalloc);
  ESL_ALLOC(graph->edge_scores_mem,   sizeof(float)        * (graph->nalloc * graph->nalloc));
 
  for (i = 0; i < graph->nalloc; i++)
    graph->edge_scores[i] = graph->edge_scores_mem + (ptrdiff_t) i * (ptrdiff_t) graph->nalloc;

  /* Allocate array of hit scores */
  ESL_ALLOC(graph->hit_scores,        sizeof(float)        * graph->nalloc);

  /* Allocate arrays of path scores */
  ESL_ALLOC(graph->path_scores,    sizeof(float)        * graph->nalloc);

  /* Allocate arrays for keeping track of which nodes have edges ans which edge is best*/
  ESL_ALLOC(graph->out_edge_cnt,   sizeof(int)          * graph->nalloc);
  ESL_ALLOC(graph->in_edge_cnt,    sizeof(int)          * graph->nalloc);
  ESL_ALLOC(graph->orig_out_edge,  sizeof(int)          * graph->nalloc);
  ESL_ALLOC(graph->orig_in_edge,   sizeof(int)          * graph->nalloc);
  ESL_ALLOC(graph->best_out_edge,  sizeof(int)          * graph->nalloc);
  ESL_ALLOC(graph->best_in_edge,   sizeof(int)          * graph->nalloc);

  /* Allocate adjacency matrix to streo the relative upstream or downstream postions of any two nodes */
  ESL_ALLOC(graph->is_upstream,       sizeof(int*)         * graph->nalloc);
  ESL_ALLOC(graph->is_upstream_mem,   sizeof(int)          * (graph->nalloc * graph->nalloc));
 
  for (i = 0; i < graph->nalloc; i++)
    graph->is_upstream[i] = graph->is_upstream_mem + (ptrdiff_t) i * (ptrdiff_t) graph->nalloc; 

  /* Allocate adjacency matrix for edge index in graph->edges */
  ESL_ALLOC(graph->edge_id,           sizeof(int*)         * graph->nalloc);
  ESL_ALLOC(graph->edge_id_mem,       sizeof(int)          * (graph->nalloc * graph->nalloc));

  for (i = 0; i < graph->nalloc; i++)
    graph->edge_id[i] = graph->edge_id_mem + (ptrdiff_t) i * (ptrdiff_t) graph->nalloc;

  /* Allocate array for SPLICE_EDGE objects */
  ESL_ALLOC(graph->edges,             sizeof(SPLICE_EDGE*) * (graph->nalloc * graph->nalloc));
  
  return eslOK;

  ERROR:
	splice_graph_destroy(graph);
    return status;
}

int
splice_graph_grow(SPLICE_GRAPH *graph)
{
  int     i,j;
  int   **is_upstream;
  int    *is_upstream_mem;
  int   **edge_id;
  int    *edge_id_mem;
  float **edge_scores;
  float  *edge_scores_mem;
  int status;
  if(graph->num_nodes == graph->nalloc) {


    graph->nalloc *= 2;
 
    ESL_ALLOC(is_upstream,     sizeof(int*)   * graph->nalloc);
    ESL_ALLOC(is_upstream_mem, sizeof(int)    * (graph->nalloc * graph->nalloc));   
 
    for (i = 0; i < graph->nalloc; i++)
      is_upstream[i] = is_upstream_mem + (ptrdiff_t) i * (ptrdiff_t) graph->nalloc;

    ESL_ALLOC(edge_scores,     sizeof(float*) * graph->nalloc);
    ESL_ALLOC(edge_scores_mem, sizeof(float)  * (graph->nalloc * graph->nalloc));

    for (i = 0; i < graph->nalloc; i++)
       edge_scores[i] = edge_scores_mem + (ptrdiff_t) i * (ptrdiff_t) graph->nalloc;

    ESL_ALLOC(edge_id,          sizeof(int*)  * graph->nalloc);
    ESL_ALLOC(edge_id_mem,      sizeof(int)   * (graph->nalloc * graph->nalloc));

    for (i = 0; i < graph->nalloc; i++)
      edge_id[i] = edge_id_mem + (ptrdiff_t) i * (ptrdiff_t) graph->nalloc;

    for(i = 0; i < graph->num_nodes; i++) {
      for(j = 0; j < graph->num_nodes; j++) {
        is_upstream[i][j] = graph->is_upstream[i][j];
        edge_scores[i][j] = graph->edge_scores[i][j];
        edge_id[i][j]     = graph->edge_id[i][j]; 
      }
    }

    if(graph->is_upstream     != NULL) free(graph->is_upstream);
    if(graph->is_upstream_mem != NULL) free(graph->is_upstream_mem);
    if(graph->edge_id         != NULL) free(graph->edge_id);
    if(graph->edge_id_mem     != NULL) free(graph->edge_id_mem);
    if(graph->edge_scores     != NULL) free(graph->edge_scores);
    if(graph->edge_scores_mem != NULL) free(graph->edge_scores_mem); 

    graph->is_upstream_mem = is_upstream_mem;
    graph->is_upstream     = is_upstream;
    graph->edge_scores_mem = edge_scores_mem;
    graph->edge_scores     = edge_scores;
    graph->edge_id_mem     = edge_id_mem;
    graph->edge_id         = edge_id;

    ESL_REALLOC(graph->hit_scores,      sizeof(float)        * graph->nalloc);
    ESL_REALLOC(graph->path_scores,     sizeof(float)        * graph->nalloc); 

    ESL_REALLOC(graph->out_edge_cnt,    sizeof(int)          * graph->nalloc);
    ESL_REALLOC(graph->in_edge_cnt,     sizeof(int)          * graph->nalloc);
    ESL_REALLOC(graph->orig_out_edge,   sizeof(int)          * graph->nalloc);
    ESL_REALLOC(graph->orig_in_edge,    sizeof(int)          * graph->nalloc);
    ESL_REALLOC(graph->best_out_edge,   sizeof(int)          * graph->nalloc);
    ESL_REALLOC(graph->best_in_edge,    sizeof(int)          * graph->nalloc);

    ESL_REALLOC(graph->edges,            sizeof(SPLICE_EDGE*) * (graph->nalloc * graph->nalloc));
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
 
  path->path_len = path_len;
  path->seq_len  = 0;

  path->node_id = NULL;
  ESL_ALLOC(path->node_id,                        sizeof(int)*(path_len));

  path->missing = NULL;
  ESL_ALLOC(path->missing,                        sizeof(int)*(path_len));

  path->upstream_spliced_amino_end     = NULL;
  path->downstream_spliced_amino_start = NULL;
  ESL_ALLOC(path->upstream_spliced_amino_end,     sizeof(int)*(path_len+1));
  ESL_ALLOC(path->downstream_spliced_amino_start, sizeof(int)*(path_len+1));

  path->upstream_spliced_nuc_end     = NULL;
  path->downstream_spliced_nuc_start = NULL;
  ESL_ALLOC(path->upstream_spliced_nuc_end,       sizeof(int)*(path_len+1));
  ESL_ALLOC(path->downstream_spliced_nuc_start,   sizeof(int)*(path_len+1));

  path->hit_scores    = NULL;
  path->signal_scores = NULL;
  ESL_ALLOC(path->hit_scores,    sizeof(float)*path_len);
  ESL_ALLOC(path->signal_scores, sizeof(float)*path_len);

  path->hits = NULL;
  ESL_ALLOC(path->hits,          sizeof(P7_HIT*)*path_len);

  return path;

  ERROR:
    splice_path_destroy(path);
    return NULL;
}


int
splice_path_split_hit(SPLICE_PATH *path, SPLICE_EDGE *edge, int split_id)
{
  int i;
  int status;

  path->path_len++;
  path->split_hits++;

  /* grow path to insert split hit */
  ESL_REALLOC(path->node_id,                        sizeof(int)     * path->path_len);
  ESL_REALLOC(path->missing,                        sizeof(int)     * path->path_len);
  ESL_REALLOC(path->upstream_spliced_amino_end,     sizeof(int)     * (path->path_len+1));
  ESL_REALLOC(path->downstream_spliced_amino_start, sizeof(int)     * (path->path_len+1));
  ESL_REALLOC(path->upstream_spliced_nuc_end,       sizeof(int)     * (path->path_len+1));
  ESL_REALLOC(path->downstream_spliced_nuc_start,   sizeof(int)     * (path->path_len+1));
  ESL_REALLOC(path->signal_scores,                  sizeof(float)   * path->path_len);
  ESL_REALLOC(path->hits,                           sizeof(P7_HIT*) * path->path_len);

  /* shift path steps down to make room for split hit */
  for (i = path->path_len; i > split_id; i--) {
    path->upstream_spliced_amino_end[i]     = path->upstream_spliced_amino_end[i-1];
    path->downstream_spliced_amino_start[i] = path->downstream_spliced_amino_start[i-1];
    path->upstream_spliced_nuc_end[i]       = path->upstream_spliced_nuc_end[i-1];
    path->downstream_spliced_nuc_start[i]   = path->downstream_spliced_nuc_start[i-1];

    if(i < path->path_len) {
      path->node_id[i]                      = path->node_id[i-1];
      path->missing[i]                      = path->missing[i-1];
      path->signal_scores[i]                = path->signal_scores[i-1];
      path->hits[i]                         = path->hits[i-1];
    }
  }

  /* Set new splice boudries */
  path->upstream_spliced_amino_end[split_id+1]     = edge->upstream_spliced_amino_end;
  path->downstream_spliced_amino_start[split_id+1] = edge->downstream_spliced_amino_start;
  path->upstream_spliced_nuc_end[split_id+1]       = edge->upstream_spliced_nuc_end;
  path->downstream_spliced_nuc_start[split_id+1]   = edge->downstream_spliced_nuc_start;
  path->signal_scores[split_id]                    = edge->signal_score;

  /* Remove new intron length from path sequence length */
  path->seq_len -= (edge->downstream_spliced_nuc_start - edge->upstream_spliced_nuc_end - 1);

  return eslOK;

   ERROR:
    splice_path_destroy(path);
    return status;

}

void
splice_path_destroy(SPLICE_PATH *path)
{

   if(path == NULL) return;

   if(path->node_id                        != NULL)
     free(path->node_id);

   if(path->missing                        != NULL)
     free(path->missing);

   if(path->upstream_spliced_amino_end     != NULL)
     free(path->upstream_spliced_amino_end);
   if(path->downstream_spliced_amino_start != NULL)
     free(path->downstream_spliced_amino_start);

   if(path->upstream_spliced_nuc_end       != NULL)
     free(path->upstream_spliced_nuc_end);
   if(path->downstream_spliced_nuc_start   != NULL)
     free(path->downstream_spliced_nuc_start);

   if(path->hit_scores    != NULL) free(path->hit_scores);
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

  if(pli->orig_nuc_idx != NULL) free(pli->orig_nuc_idx);
  pli->orig_nuc_idx = NULL;  

  if(pli->hit != NULL && pli->hit->dcl != NULL) { 
    p7_alidisplay_Destroy(pli->hit->dcl->ad);
    p7_trace_splice_Destroy(pli->hit->dcl->tr);
  }
  if(pli->hit != NULL)
    p7_hit_Destroy(pli->hit);
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

  int              i,j;
  int              hit_cnt;
  int              range_cnt;
  int              success;
  int              num_hits_processed; 
  int              prev_num_hits_processed;
  int              is_sorted_by_sortkey;
  int              revcomp;
  int64_t          seqidx;
  int             *hits_processed;
  int             *hit_spliced;
  int64_t         *range_bound_mins;
  int64_t         *range_bound_maxs;
  SPLICE_PIPELINE *pli;
  TARGET_RANGE    *curr_target_range;
  TARGET_RANGE    *prev_target_range;
  SPLICE_GRAPH    *graph;
  SPLICE_EDGE     *edge;
  SPLICE_PATH     *path;
  P7_HIT          *hit; 
  ESL_SQ          *target_seq;
  int              status;

  hits_processed    = NULL;
  hit_spliced      = NULL;
  range_bound_mins = NULL;
  range_bound_maxs = NULL;
 
  /* Sort hits by sequence, strand and nucleotide positon. 
   * Note if P7_TOPHITS was perviously sorted by sortkey 
   * so we can retur it to it's original state */  
  is_sorted_by_sortkey = tophits->is_sorted_by_sortkey;
  if ((status = p7_tophits_SortBySeqidxAndAlipos(tophits)) != eslOK) goto ERROR;
     
  pli = splice_pipeline_create(go, om->M, om->M * 3);
  pli->bg = p7_bg_Create(om->abc);
  if (pli->do_biasfilter) p7_bg_SetFilter(pli->bg, om->M, om->compo);  

  /* Arrays to keep track of hit status */
  ESL_ALLOC(hits_processed, tophits->N*sizeof(int));
  ESL_ALLOC(hit_spliced,   tophits->N*sizeof(int));

  /* Check for any hits that are duplicates and mark them as processed */
  hit_cnt = 0;
  num_hits_processed = 0;
  for(i = 0; i < tophits->N; i++) {
    tophits->hit[i]->dcl->ad->exon_cnt = 1;
    if (tophits->hit[i]->flags & p7_IS_DUPLICATE) {
      hits_processed[i] = 1;
      num_hits_processed++;
    }
    else { 
      hits_processed[i] = 0;
      hit_cnt++;
    }  
    hit_spliced[i]   = 0;
  }
  
  //p7_tophits_AliScores(stdout, om->name, tophits);
  /* Arrays for keeping track of allowable target range mins 
   * and maxes set by previous ranges or splicings */
  ESL_ALLOC(range_bound_mins, hit_cnt * sizeof(int64_t));
  ESL_ALLOC(range_bound_maxs, hit_cnt * sizeof(int64_t));
  
  prev_target_range = NULL;
  prev_num_hits_processed = -1;
  /* loop through until all hits have been processed */
  range_cnt = 0;
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

    /* If the curr_hit is on a different sequences or strand than the prev_target_range reset teh range bounds */
    if (prev_target_range != NULL) {
      if(prev_target_range->seqidx != seqidx || prev_target_range->revcomp != revcomp) 
        range_cnt = 0;
    }

    /* Build a target range from a set of hits. */
    curr_target_range = build_target_range(tophits, hits_processed, &num_hits_processed, range_bound_mins, range_bound_maxs, range_cnt, seqidx, revcomp);
   
    if(curr_target_range->th->N == 0) {
      target_range_destroy(curr_target_range); 
      continue;
    }
  //printf("Target %s strand %c range %d to %d\n", curr_target_range->seqname, (curr_target_range->revcomp ? '-' : '+'), curr_target_range->start, curr_target_range->end);

    range_cnt++;
     //target_range_dump(stdout, curr_target_range, TRUE);
    /* Fetch a sub-sequence thar corresponds to the target range */
    target_seq = get_target_range_sequence(curr_target_range, seq_file); 
   //printf("Sequence n %d L %d start %d end %d\n", target_seq->n, target_seq->L, target_seq->start, target_seq->end);

    /* Build graph from available hits */
    graph = splice_graph_create();   
    graph->revcomp = curr_target_range->revcomp;
       
    if ((status = fill_graph_with_nodes(graph, curr_target_range->th, gm)) != eslOK) goto ERROR;
 
    /* Connect the nodes */
    for(i = 0; i < curr_target_range->th->N; i++) {

      for(j = 0; j < curr_target_range->th->N; j++) {
        if (!graph->is_upstream[i][j]) {
          graph->edge_scores[i][j] = -eslINFINITY;
          graph->edge_id[i][j]     = -1;
          continue;
        }
	 
        edge = connect_nodes_with_edges(curr_target_range->th->hit[i], curr_target_range->th->hit[j], gm, hmm, pli->bg, gcode, target_seq, graph->revcomp);
        if(edge == NULL) {
          graph->edge_scores[i][j] = -eslINFINITY;
          graph->edge_id[i][j]     = -1; 
        }
        else {
          if(ESL_MIN(graph->hit_scores[i], graph->hit_scores[j]) + edge->splice_score > 0) {
            edge->upstream_node_id   = i;
            edge->downstream_node_id = j;
            graph->orig_out_edge[i]  = TRUE;
            graph->orig_in_edge[j]   = TRUE;
            add_edge_to_graph(graph, edge);
          }
          else {
            graph->edge_scores[i][j] = -eslINFINITY;
            graph->edge_id[i][j]     = -1;
            free(edge);
          }
        }
      }
    }

    
   //graph_dump(stdout, graph, target_seq, TRUE); 

    fill_holes_in_graph(graph,  curr_target_range, gm, hmm, pli->bg, target_seq, gcode); 

   //graph_dump(stdout, graph, target_seq, TRUE);
    check_for_loops(graph, curr_target_range->th);
    path = evaluate_paths(graph, curr_target_range->th, target_seq, curr_target_range->orig_N);
    split_hits_in_path(graph, path, gm, hmm, pli->bg, gcode, target_seq, curr_target_range->orig_N); 
    //target_range_dump(stdout, curr_target_range, TRUE);
    if(path->path_len > 1)
      splice_path(graph, path, pli, om, scoredata, target_seq, gcode, db_nuc_cnt, curr_target_range->orig_N, &success);
    else
      success = FALSE;    

    /* Reset the range_bounds around the spliced hit and set all hits from 
     * the target_range that are outside this new range to unprocessed */
    if (success) {
       
      tophits->nreported -= path->path_len-(path->split_hits+1);
      tophits->nreported -= path->path_len-(path->split_hits+1);
 
      range_bound_mins[range_cnt-1] = ESL_MIN(pli->hit->dcl->iali, pli->hit->dcl->jali);
      range_bound_maxs[range_cnt-1] = ESL_MAX(pli->hit->dcl->iali, pli->hit->dcl->jali);

      
      release_hits_from_target_range(curr_target_range, hits_processed, &num_hits_processed, range_bound_mins[range_cnt-1], range_bound_maxs[range_cnt-1]);           

      pli->hit->dcl = NULL;
    } 
	  else {
      
      /* If the path failed but there are other hits in the graph that were not in the path 
       * we want to release them to be considered in another target range */
      if (graph->revcomp) {     
        range_bound_mins[range_cnt-1] = target_seq->n - path->upstream_spliced_nuc_end[path->path_len] + target_seq->end;
        range_bound_maxs[range_cnt-1] = target_seq->n - path->downstream_spliced_nuc_start[0] + target_seq->end;
      }
      else {
        range_bound_mins[range_cnt-1] = path->downstream_spliced_nuc_start[0] + target_seq->start - 1;
        range_bound_maxs[range_cnt-1] = path->upstream_spliced_nuc_end[path->path_len] + target_seq->start - 1;;
      }
      release_hits_from_target_range(curr_target_range, hits_processed, &num_hits_processed, range_bound_mins[range_cnt-1], range_bound_maxs[range_cnt-1]);  

    }
   
   
    target_range_destroy(prev_target_range);
    prev_target_range = curr_target_range; 
    splice_path_destroy(path); 
    splice_pipeline_reuse(pli);
    splice_graph_destroy(graph);
    esl_sq_Destroy(target_seq);

  }

  /* Leave only footprints */
  if (is_sorted_by_sortkey)
    if ((status = p7_tophits_SortBySortkey(tophits)) != eslOK) goto ERROR;   

  splice_pipeline_destroy(pli);
  target_range_destroy(prev_target_range);

  if(hits_processed    != NULL) free(hits_processed);
  if(hit_spliced      != NULL) free(hit_spliced);
  if(range_bound_mins != NULL) free(range_bound_mins);
  if(range_bound_maxs != NULL) free(range_bound_maxs);
  return status;

  ERROR:
    splice_pipeline_destroy(pli);
    target_range_destroy(prev_target_range);
    splice_graph_destroy(graph);
    splice_path_destroy(path);
    esl_sq_Destroy(target_seq);
    if(hits_processed    != NULL) free(hits_processed);
    if(hit_spliced      != NULL) free(hit_spliced); 
    if(range_bound_mins != NULL) free(range_bound_mins);
    if(range_bound_maxs != NULL) free(range_bound_maxs);
    return status;
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
TARGET_RANGE *
build_target_range (const P7_TOPHITS *tophits, int *hits_processed, int *num_hits_processed, int64_t *range_bound_mins, int64_t *range_bound_maxs, int range_num, int64_t seqidx, int revcomp)
{

  int           i;
  int           num_scores;
  int           max_sc_idx;
  int           seed_hit_idx;
  int           curr_hit_idx;
  int           add_hit;
  int           hits_processed_cnt;
  int           seed_hmm_from, seed_hmm_to;
  int64_t       seed_ali_min, seed_ali_max;
  int64_t       curr_ali_min, curr_ali_max;
  int64_t       lower_bound, upper_bound;
  int64_t       curr_range_min, curr_range_max;
  int64_t       new_range_min, new_range_max;
  int          *hit_scores_idx;
  float        *hit_scores;
  TARGET_RANGE *target_range;
  P7_HIT       *curr_hit;
  P7_HIT       *seed_hit;
  int           status;

  hits_processed_cnt = *num_hits_processed;
  hit_scores     = NULL;
  hit_scores_idx = NULL;
 
  /* Arrays for hit score sorting and tracking*/
  ESL_ALLOC(hit_scores    , (tophits->N - hits_processed_cnt) * sizeof(float));
  ESL_ALLOC(hit_scores_idx, (tophits->N - hits_processed_cnt) * sizeof(int));


  /* Add scores from all unprocessed hits ont the current sequence and strand to hit_scores */
  num_scores = 0;
  for (i = 0; i < tophits->N; i++) {

    if (hits_processed[i]) continue; 

    /* Check if hit is on the same sequnce and strand as first unprocessed hit */
    curr_hit = tophits->hit[i];
    if (curr_hit->seqidx != seqidx) continue;
    
    if (revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
    if ((!revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue; 

    /* Add score to socres arrray and keep track of tophits index */
    hit_scores[num_scores]     = curr_hit->score;
    hit_scores_idx[num_scores] = i;
    num_scores++; 
  }  
 
  /* Now that we know the maximum number of hits on the current sequence 
   * and srand we can initialize the TARGET_RANGE */
  target_range = target_range_create(num_scores);
 
  /* Get the index of the Maximum score in the hit score array */
  max_sc_idx = esl_vec_FArgMax(hit_scores,  num_scores);

  /* Get Maximum scoring hit to act at seed hit for target range and mark hit as processed */
  seed_hit_idx = hit_scores_idx[max_sc_idx];
  seed_hit = tophits->hit[seed_hit_idx];
  hits_processed[seed_hit_idx] = 1;
  hits_processed_cnt++; 
 
  /* If the seed hit is bellow the reporting threshold then we are done with this sequence 
   * and strand. Mark any remaining hits as processed and return the empty target range */
  if(!(seed_hit->flags & p7_IS_REPORTED)) {
    for(i = 0; i < tophits->N; i++) {
      if(tophits->hit[i]->seqidx == seqidx && hits_processed[i] == 0) {
        if((revcomp && tophits->hit[i]->dcl->jali <  tophits->hit[i]->dcl->iali) ||
          (!revcomp && tophits->hit[i]->dcl->iali <  tophits->hit[i]->dcl->jali)) {
            hits_processed[i] = 1;
            hits_processed_cnt++;
        }
      }
    }
    if(hit_scores     != NULL) free(hit_scores);
    if(hit_scores_idx != NULL) free(hit_scores_idx);
    *num_hits_processed = hits_processed_cnt;
    return target_range;
  }

  seed_ali_min =  ESL_MIN(seed_hit->dcl->iali, seed_hit->dcl->jali);
  seed_ali_max =  ESL_MAX(seed_hit->dcl->iali, seed_hit->dcl->jali);

  lower_bound  = 0;
  upper_bound  = seed_hit->dcl->ad->L;

  /* Get any upper and lower bounds set by previous ranges or splicings */
  for( i = 0; i < range_num; i++) {

     /* If seed hit exists within range bound return empty target range */
    if((seed_ali_min >= range_bound_mins[i] && seed_ali_min <= range_bound_maxs[i]) ||
       (seed_ali_max >= range_bound_mins[i] && seed_ali_max <= range_bound_maxs[i])) {
      if(hit_scores     != NULL) free(hit_scores);
      if(hit_scores_idx != NULL) free(hit_scores_idx);
      *num_hits_processed = hits_processed_cnt;
      return target_range;
    }
    //printf("range_bound_mins[i] %d range_bound_maxs[i] %d\n", range_bound_mins[i], range_bound_maxs[i]);
    //printf("seed_ali_min %d seed_ali_max %d\n", seed_ali_min, seed_ali_max);
    if(range_bound_maxs[i] < seed_ali_min && range_bound_maxs[i] > lower_bound)
      lower_bound = range_bound_maxs[i];

    if(range_bound_mins[i] > seed_ali_max && range_bound_mins[i] < upper_bound)
      upper_bound = range_bound_mins[i];
  }
 // printf("lower_bound %d  upper_bound %d \n", lower_bound, upper_bound);
  /* Set target range values and add seed hit to P7_TOPHITS */ 
  target_range->revcomp         = revcomp;
  target_range->seqidx          = seed_hit->seqidx;
  target_range->seqname         = seed_hit->name;
  target_range->th->hit[0]      = seed_hit;
  target_range->orig_hit_idx[0] = seed_hit_idx;
  target_range->th->N++;

  seed_hmm_from = seed_hit->dcl->ihmm;
  seed_hmm_to   = seed_hit->dcl->jhmm;

  
  curr_range_min = seed_ali_min; 
  curr_range_max = seed_ali_max;

  for (i = 0; i < num_scores; i++) {

    /* skip the seed hit */
    if (i == max_sc_idx) continue;

   /* Get a hit that has already been check for sequence and strand
    * compatiability with the seed hit and confirmed not processed*/
    curr_hit_idx = hit_scores_idx[i];
    curr_hit = tophits->hit[curr_hit_idx];

    curr_ali_min = ESL_MIN(curr_hit->dcl->iali, curr_hit->dcl->jali);
    curr_ali_max = ESL_MAX(curr_hit->dcl->iali, curr_hit->dcl->jali);
    
    /* Check if current hit is outside upper and lower bounds set by previous ranges or splicings */
    if ( curr_ali_min < lower_bound || curr_ali_max > upper_bound) continue;

    /* Check if the current hit is within the maximum distance from the seed hit */
    if ( curr_ali_min < (seed_ali_min - MAX_TARGET_RANGE_EXT)) continue;
    if ( curr_ali_max > (seed_ali_max + MAX_TARGET_RANGE_EXT)) continue;


    add_hit = FALSE;
    /* Is current hit upstream of seed hit in hmm positions */
    if ( curr_hit->dcl->ihmm < seed_hmm_from ) {

      /* Is current hit upstream of seed hit in seq positions */
      if ( (!revcomp) && curr_ali_min < seed_ali_min) {
        add_hit = TRUE;
        curr_range_min = ESL_MIN(curr_range_min, curr_hit->dcl->iali);
      }
      else if ( revcomp && curr_ali_max > seed_ali_max) {
        add_hit = TRUE;
        curr_range_max = ESL_MAX(curr_range_max, curr_ali_max);
      }
    }
    /* Is current hit downstream of seed hit in hmm positions */
    if ( curr_hit->dcl->jhmm > seed_hmm_to ) {
      /* Is current hit downstream of seed hit in seq positions */
      if ( (!revcomp) && curr_ali_max > seed_ali_max) {
        add_hit = TRUE;
        curr_range_max = ESL_MAX(curr_range_max, curr_ali_max);
      }
      else if ( revcomp && curr_ali_min < seed_ali_min) {
        add_hit = TRUE;
        curr_range_min = ESL_MIN(curr_range_min, curr_ali_min);
      }
    }

    /* If hit is upstream  or donwstream on both hmm and seq positions add it to target range */
    if(add_hit) {
      target_range->th->hit[target_range->th->N] = curr_hit;
      target_range->orig_hit_idx[target_range->th->N] = curr_hit_idx;
      target_range->th->N++;
      hits_processed[curr_hit_idx] = 1;
      hits_processed_cnt++;
    }
  }
 

  new_range_min = curr_range_min;
  new_range_max = curr_range_max;

  /* Add any hits that fall into the current range coords to the target range */
  for (i = 0; i < num_scores; i++) {

    curr_hit_idx = hit_scores_idx[i];
    curr_hit = tophits->hit[curr_hit_idx];

    if (hits_processed[curr_hit_idx]) continue;

    curr_ali_min = ESL_MIN(curr_hit->dcl->iali, curr_hit->dcl->jali);
    curr_ali_max = ESL_MAX(curr_hit->dcl->iali, curr_hit->dcl->jali);   

    if( (curr_ali_min < curr_range_max && curr_ali_max > curr_range_max) ||
        (curr_ali_max > curr_range_min && curr_ali_min < curr_range_max)) {
      
        target_range->th->hit[target_range->th->N] = curr_hit;
        target_range->orig_hit_idx[target_range->th->N] = curr_hit_idx;
        target_range->th->N++;
        hits_processed[curr_hit_idx] = 1;
        hits_processed_cnt++;
  
       new_range_min = ESL_MIN(new_range_min, curr_ali_min);
       new_range_max = ESL_MAX(new_range_max, curr_ali_max);
    }
  } 

  target_range->orig_N =  target_range->th->N;

  range_bound_mins[range_num] = new_range_min;
  range_bound_maxs[range_num] = new_range_max;
  
  /* Extend the target range to make room for possible missing terminal exons */  
  target_range->start = ESL_MAX(new_range_min-TERM_RANGE_EXT, lower_bound);
  target_range->end   = ESL_MIN(new_range_max+TERM_RANGE_EXT, upper_bound); 

  *num_hits_processed = hits_processed_cnt;

  if(hit_scores     != NULL) free(hit_scores);
  if(hit_scores_idx != NULL) free(hit_scores_idx);

  return target_range;

  ERROR:
    if(hit_scores     != NULL) free(hit_scores);
    if(hit_scores_idx != NULL) free(hit_scores_idx);
    return NULL;
}


/* Function: get_target_range_sequence
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
get_target_range_sequence(const TARGET_RANGE *target_range, const ESL_SQFILE *seq_file)
{
  
  ESL_SQ     *target_seq;
  ESL_SQFILE *tmp_file;

  /* Open seq file and ssi */
  esl_sqfile_Open(seq_file->filename,seq_file->format,NULL,&tmp_file);
  esl_sqfile_OpenSSI(tmp_file,NULL);

  /* Get basic sequence info */
  target_seq = esl_sq_Create();
  esl_sqio_FetchInfo(tmp_file, target_range->seqname, target_seq);

  target_seq->start = target_range->start;
  target_seq->end   = target_range->end;
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

  esl_sq_SetName(target_seq, target_range->seqname);
 
  if(target_range->revcomp)
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
  int        hit_idx;
  int        released;
  P7_TOPHITS *th;
 
  th = target_range->th;
  released = 0;
        
  if (target_range->revcomp) {     

    for(i = 0; i < target_range->orig_N; i++) {
      if      (th->hit[i]->dcl->jali > range_bound_max) {
        hit_idx =  target_range->orig_hit_idx[i];
        hits_processed[hit_idx] = 0;
        released++;
      }
      else if (th->hit[i]->dcl->iali < range_bound_min) {
        hit_idx =  target_range->orig_hit_idx[i];
        hits_processed[hit_idx] = 0;
        released++; 
      }
    }
  } 
  else {

    for(i = 0; i < target_range->orig_N; i++) {
      if      (th->hit[i]->dcl->iali > range_bound_max) {
        hit_idx =  target_range->orig_hit_idx[i];
        hits_processed[hit_idx] = 0;
        released++; 
      }
      else if (th->hit[i]->dcl->jali < range_bound_min) { 
        hit_idx =  target_range->orig_hit_idx[i];
        hits_processed[hit_idx] = 0;
        released++; 
      }
    }
  } 
 
  *num_hits_processed -= released;
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
fill_graph_with_nodes(SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm)
{

  int     i,j;
  float   sum_ali_sc;
  int     status;

  status = splice_graph_create_nodes(graph, th->N);

  if (status != eslOK) return status; 

  for(i = 0; i < th->N; i++) {
    for(j = 0; j < th->N; j++) {
      if(i == j) {
        graph->is_upstream[i][j] = FALSE;
        continue;
      }
       
      if(th->hit[i]->dcl->ihmm >= th->hit[j]->dcl->ihmm &&
         th->hit[i]->dcl->jhmm >= th->hit[j]->dcl->jhmm) {
        graph->is_upstream[i][j] = FALSE;
        continue;
      }

      if (( graph->revcomp  && th->hit[j]->dcl->jali > th->hit[i]->dcl->iali) ||
         ((!graph->revcomp) && th->hit[j]->dcl->jali < th->hit[i]->dcl->iali)) {
         graph->is_upstream[i][j] = FALSE;
         continue;
      } 
      
      graph->is_upstream[i][j] = TRUE; 
    }
  
    sum_ali_sc = 0.;
    
    for(j = 0; j < th->hit[i]->dcl->ad->N; j++)  
      sum_ali_sc += th->hit[i]->dcl->scores_per_pos[j];

    graph->hit_scores[i]  = sum_ali_sc;
    graph->path_scores[i] = -eslINFINITY;

	graph->out_edge_cnt[i]  = 0;
    graph->in_edge_cnt[i]   = 0;
    graph->orig_out_edge[i] = 0;
    graph->orig_in_edge[i]  = 0;
    graph->best_out_edge[i] = -1;
    graph->best_in_edge[i]  = -1;   
  }

  return eslOK; 

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
connect_nodes_with_edges(P7_HIT *upstream_hit, P7_HIT *downstream_hit, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, ESL_GENCODE *gcode, ESL_SQ *target_seq, int revcomp)
{

  int          up_amino_start, up_amino_end;
  int          down_amino_start, down_amino_end;
  int          overlap;
  int          num_ext_aminos;
  SPLICE_EDGE *edge;
  int          status;

  up_amino_start = upstream_hit->dcl->ihmm;
  up_amino_end   = upstream_hit->dcl->jhmm;

  down_amino_start = downstream_hit->dcl->ihmm;
  down_amino_end   = downstream_hit->dcl->jhmm;

  /* Is the downstream node close enough to the upstream edge to be the next exon  */
  if (up_amino_end + MAX_AMINO_EXT < down_amino_start)
    return NULL;

  /* If the overlap between the nodes is greater than MAX_AMINO_OVERLAP disregard */
  if(down_amino_start + MAX_AMINO_OVERLAP < up_amino_end)
    return NULL;

  overlap = ESL_MAX(MIN_AMINO_OVERLAP, down_amino_start - up_amino_end + 2);
  num_ext_aminos = overlap - (up_amino_end - down_amino_start + 1);

  if (num_ext_aminos > 0) {
     if(( revcomp  && (upstream_hit->dcl->jali - downstream_hit->dcl->iali) < num_ext_aminos*3) ||
       ((!revcomp) && (downstream_hit->dcl->iali - upstream_hit->dcl->jali) < num_ext_aminos*3))
       return NULL;
  }

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
  if(revcomp) {
    edge->upstream_nuc_end     += 2;
    edge->downstream_nuc_start -= 2;
  }
  else {
    edge->upstream_nuc_end     += 2;
    edge->downstream_nuc_start -= 2;
  }

  if(edge->downstream_nuc_start - edge->upstream_nuc_end > MAX_INTRON_LEN) {
    free(edge);
    return NULL;
  }

  if ((status = find_optimal_splice_site (edge, upstream_hit->dcl, downstream_hit->dcl, gm, hmm, bg, gcode, target_seq)) != eslOK) goto ERROR;

  if(edge->splice_score == -eslINFINITY) {
     free(edge);
     return NULL;
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
get_overlap_nuc_coords (SPLICE_EDGE *edge, P7_DOMAIN *upstream, P7_DOMAIN *downstream, ESL_SQ *target_seq, int revcomp)
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
  edge->upstream_trace_end = z2+1;
  while(curr_hmm_pos >= edge->overlap_amino_start) {
	
    if      (up_trace->st[z2] == p7T_M) {
      edge->upstream_nuc_start -= up_trace->c[z2] * strand;
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
find_optimal_splice_site (SPLICE_EDGE *edge, P7_DOMAIN *upstream, P7_DOMAIN *downstream, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, ESL_GENCODE *gcode, ESL_SQ *target_seq)
{

  int         n,z;
  int         unp,dnp;
  float       lost_sc;
  float      *signal_scores;
  P7_HMM     *sub_hmm;
  P7_PROFILE *sub_model;
  int         status;

  /* Get the summed ali score for the overlap region covered by the existing alignments */
  lost_sc = 0.;
  z = edge->upstream_trace_end;
  n = upstream->tr->tto[0] - upstream->tr->tfrom[0] - 2;
  while(upstream->tr->k[z] >= edge->overlap_amino_start) {
    lost_sc += upstream->scores_per_pos[n];
    n--;
    z--;
  }

  z = edge->downstream_trace_start;
  n = 0;
  while(downstream->tr->k[z] && downstream->tr->k[z] <= edge->overlap_amino_end) {
    lost_sc += downstream->scores_per_pos[n];
    n++;
    z++;
  }

  /*Get a submodel that covers the overlap region */
  sub_hmm    = extract_sub_hmm(hmm, edge->overlap_amino_start, edge->overlap_amino_end);
  sub_model  = p7_profile_Create (sub_hmm->M, sub_hmm->abc);
  p7_ProfileConfig(sub_hmm, bg, sub_model, 100, p7_UNIGLOCAL);

  /* Initialize splice signal score array (in bitscore) */
  signal_scores = NULL;
  ESL_ALLOC(signal_scores, sizeof(float) * p7S_SPLICE_SIGNALS);
  p7_splice_SignalScores(signal_scores);

  /* Scan the overlap nucleotides for splice signals */
  for(unp = edge->upstream_nuc_start; unp < edge->upstream_nuc_end; unp++) {

    /*GT-AG*/
    if(target_seq->dsq[unp] == 2 && target_seq->dsq[unp+1] == 3) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {
        if(target_seq->dsq[dnp-1] == 0 && target_seq->dsq[dnp] == 2) {
          //printf("p7S_GTAG\n");
          if (dnp - unp > MIN_INTRON_LEN)
            select_splice_option(edge, gm, sub_model, gcode, target_seq, signal_scores[p7S_GTAG], unp-1, dnp+1);
        }
      }
    }
    /*GC-AG*/
    if(target_seq->dsq[unp] == 2 && target_seq->dsq[unp+1] == 1) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {
        if(target_seq->dsq[dnp-1] == 0 && target_seq->dsq[dnp] == 2) {
          //printf("p7S_GCAG\n");
          if (dnp - unp > MIN_INTRON_LEN)
                   select_splice_option(edge, gm, sub_model, gcode, target_seq, signal_scores[p7S_GCAG], unp-1, dnp+1);
        }
      }
    }
    /*AT-AC*/
    if(target_seq->dsq[unp] == 0 && target_seq->dsq[unp+1] == 3) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {
        if(target_seq->dsq[dnp-1] == 0 && target_seq->dsq[dnp] == 1) {
          //printf("p7S_ATAC\n");
          if (dnp - unp > MIN_INTRON_LEN)
                   select_splice_option(edge, gm, sub_model, gcode, target_seq, signal_scores[p7S_ATAC], unp-1, dnp+1);
        }
      }
    }
  }

  if(edge->splice_score != -eslINFINITY)
    edge->splice_score -= lost_sc;

  if(signal_scores != NULL) free(signal_scores);
  if(sub_hmm       != NULL) p7_hmm_Destroy(sub_hmm);
  if(sub_model     != NULL) p7_profile_Destroy(sub_model);

  return eslOK;

  ERROR:
    if(signal_scores != NULL) free(signal_scores);
    if(sub_hmm       != NULL) p7_hmm_Destroy(sub_hmm);
    if(sub_model     != NULL) p7_profile_Destroy(sub_model);
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
select_splice_option (SPLICE_EDGE *edge, P7_PROFILE *gm, P7_PROFILE *sub_model, ESL_GENCODE *gcode, ESL_SQ *target_seq, float signal_score, int up_nuc_pos, int down_nuc_pos)
{

  int         i,z;
  int         nuc_seq_idx;
  int         nuc_seq_len;
  int         amino_len;
  int         amino;
  int         N_cnt, C_cnt;
  int         last_state;
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

  /* The viterbi alignment does not curently allow frameshifts */
  if(nuc_seq_len % 3) return eslOK;

  /* If the splice signal is right at the overlap boundries all amino acid positions are deletions */
  if(nuc_seq_len == 0) {
    sum_ali_sc = gm->tsc[(edge->overlap_amino_start-1) * p7P_NTRANS + p7H_MD];
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

  ESL_ALLOC(nuc_dsq,   sizeof(ESL_DSQ) * nuc_seq_len);
  nuc_seq_idx = 0;
  for(i = edge->upstream_nuc_start; i <= up_nuc_pos; i++) {
    nuc_dsq[nuc_seq_idx] = target_seq->dsq[i];
    nuc_seq_idx++;
  }
  upstream_amino_cnt = (nuc_seq_idx+1) / 3;
  if((nuc_seq_idx+1) % 3 > 1)  upstream_amino_cnt++;
  for(i = down_nuc_pos; i <= edge->downstream_nuc_end; i++) {
    nuc_dsq[nuc_seq_idx] = target_seq->dsq[i];
    nuc_seq_idx++;
  }

  /* Translate overalp nucleotides to amino sequence */
  amino_len = nuc_seq_len / 3;
  ESL_ALLOC(amino_dsq, sizeof(ESL_DSQ) * (amino_len+2));

  amino_dsq[0] = eslDSQ_SENTINEL;
  nuc_seq_idx = 0;
  for(i = 1; i <= amino_len; i++) {
    amino = esl_gencode_GetTranslation(gcode,&nuc_dsq[nuc_seq_idx]);
    amino_dsq[i] = amino;
    nuc_seq_idx+=3;
  }
  amino_dsq[amino_len+1] = eslDSQ_SENTINEL;

  /* Align translated overlap amino acids to submodel */
  vit_mx = p7_gmx_Create(sub_model->M, 20);
  tr = p7_trace_Create();
  p7_gmx_GrowTo(vit_mx, sub_model->M, amino_len);
  p7_GViterbi(amino_dsq, amino_len, sub_model, vit_mx, &vitsc);

  if (vitsc != -eslINFINITY) {
    p7_GTrace(amino_dsq, amino_len, sub_model, vit_mx, tr);
    p7_trace_Index(tr);

    N_cnt = C_cnt = 0;
    sum_ali_sc = 0.;
    /* Get ali score form trace. Treat N and C states as insert states */
    for(z = 0; z < tr->N; z++) {
      switch(tr->st[z]) {
        case p7T_N: N_cnt++;
                    if      (N_cnt == 2) sum_ali_sc += gm->tsc[(edge->overlap_amino_start-1) * p7P_NTRANS + p7H_MI];
                    else if (N_cnt > 2)  sum_ali_sc += gm->tsc[(edge->overlap_amino_start-1) * p7P_NTRANS + p7H_II];
                    if(tr->i[z] == upstream_amino_cnt) upstream_amino_end = tr->k[z] + edge->overlap_amino_start - 1;
                    break;
        case p7T_C: C_cnt++;
                    if      (C_cnt == 2) sum_ali_sc += gm->tsc[edge->overlap_amino_end * p7P_NTRANS + p7H_MI];
                    else if (C_cnt > 2)  sum_ali_sc += gm->tsc[edge->overlap_amino_end * p7P_NTRANS + p7H_II];
                    if(tr->i[z] == upstream_amino_cnt) upstream_amino_end = tr->k[z] + edge->overlap_amino_start - 1;
                    break;
        case p7T_M: sum_ali_sc += (tr->k[z] == 1 ? 0. : sub_model->tsc[(tr->k[z]-1) * p7P_NTRANS + p7H_MM]);
                    sum_ali_sc += sub_model->rsc[amino_dsq[tr->i[z]]][tr->k[z] * p7P_NR + p7P_MSC];
                    if(tr->i[z] == upstream_amino_cnt) upstream_amino_end = tr->k[z] + edge->overlap_amino_start - 1;
                    N_cnt = 0;
                    break;
        case p7T_I: if (last_state == p7T_M) sum_ali_sc += sub_model->tsc[(tr->k[z]-1) * p7P_NTRANS + p7H_MI];
                    else                     sum_ali_sc += sub_model->tsc[(tr->k[z]-1) * p7P_NTRANS + p7H_II];
                    if(tr->i[z] == upstream_amino_cnt) upstream_amino_end = tr->k[z] + edge->overlap_amino_start - 1;
                    break;
        case p7T_D: if (last_state == p7T_M) sum_ali_sc += sub_model->tsc[(tr->k[z]-1) * p7P_NTRANS + p7H_MD];
                    else                     sum_ali_sc += sub_model->tsc[(tr->k[z]-1) * p7P_NTRANS + p7H_DD];
                    break;
        case p7T_S:
        case p7T_B:
        case p7T_E:
        case p7T_T: break;
        default:    ESL_EXCEPTION(eslEINVAL, "no such state at find_optimal_splice_site");
      }
      last_state = tr->st[z];
    }

    overlap_sc = sum_ali_sc + signal_score;
    if (overlap_sc > edge->splice_score) {
      edge->signal_score = signal_score;
      edge->splice_score = overlap_sc;
      edge->upstream_spliced_nuc_end = up_nuc_pos;
      edge->downstream_spliced_nuc_start = down_nuc_pos;
      edge->upstream_spliced_amino_end = upstream_amino_end;
      edge->downstream_spliced_amino_start = upstream_amino_end + 1;
    }
  }
  /* Case where all transalted aminos are stop codons - treat them as instertions and all hmm positions as deletions */
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
  if(tr        != NULL) p7_trace_Destroy(tr);
  return eslOK;

  ERROR:
    if(nuc_dsq   != NULL) free(nuc_dsq);
    if(amino_dsq != NULL) free(amino_dsq);
    if(vit_mx    != NULL) p7_gmx_Destroy(vit_mx);
    if(tr        != NULL) p7_trace_Destroy(tr);
    return status;
}



/*  Function: ali_score_at_postion
 *
 *  Synopsis: Get alignment score fo a particular model postions, amino emission, and state transition.
 *
 *  Purpose:  Given an amino acid, model poistions and states, calculate an alignment score that 
 *            is the sum of the emssiosn score (for match states) and the transtions score. 
 *
 *  Returns:  The score.
 *
 */

float 
ali_score_at_postion ( P7_PROFILE *gm, int amino, int model_pos, int trans_pos, int prev_state, int curr_state)
{

  int   transition;
  float transition_score;
  float emission_score;

  /* Emission score */
  if(curr_state == p7T_M)
    emission_score = p7P_MSC(gm,model_pos,amino);
  else
    emission_score = 0.;

  if (isinf(emission_score))
    emission_score = 0.;

  /* Transtion Score */
  transition = p7P_NTRANS;
  if      (curr_state == p7T_M) {
    if      (prev_state == p7T_M)  transition = p7P_MM;
    else if (prev_state == p7T_I)  transition = p7P_IM;
    else if (prev_state == p7T_D)  transition = p7P_DM;
  }
  else if (curr_state == p7T_I) {
    if      (prev_state == p7T_M)  transition = p7P_MI;
    else if (prev_state == p7T_I)  transition = p7P_II;
  }
  else if (curr_state == p7T_D) {
    if      (prev_state == p7T_M)  transition = p7P_MD;
    else if (prev_state == p7T_D)  transition = p7P_DD;
  }


  if (transition != p7P_NTRANS)
    transition_score = p7P_TSC(gm,trans_pos,transition);
  else
    transition_score = 0.; //special case for first or last potition in overlap

  return transition_score + emission_score;

}

int
add_edge_to_graph(SPLICE_GRAPH *graph, SPLICE_EDGE *edge)
{

  int up_node_id;
  int down_node_id;

  up_node_id   = edge->upstream_node_id;
  down_node_id = edge->downstream_node_id; 
 
  graph->out_edge_cnt[up_node_id]++; 
  graph->in_edge_cnt[down_node_id]++;
  graph->edge_scores[up_node_id][down_node_id] = edge->splice_score;
  
  graph->edges[graph->num_edges]               = edge;
  graph->edge_id[up_node_id][down_node_id]     = graph->num_edges;
  graph->num_edges++; 
 
  return eslOK; 
}

int
fill_holes_in_graph(SPLICE_GRAPH *graph, TARGET_RANGE *target_range, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, ESL_SQ *target_seq, ESL_GENCODE *gcode)
{
  int up,down;
  int hit;
  int prev_N;
  int hmm_gap_len;
  int seq_gap_len;
  int new_edge;
  int num_hits;
  int duplicate;
  int dup_hits;
  P7_TOPHITS *th;
  P7_HIT     **top_ten;
  SPLICE_GAP *gap;

  th = target_range->th;

  for(up = 0; up < target_range->orig_N; up++) {
 
    for(down = 0; down < target_range->orig_N; down++) {
      num_hits = 0;
      if(!graph->out_edge_cnt[up] && !graph->in_edge_cnt[down] && graph->is_upstream[up][down])  { 
        hmm_gap_len = th->hit[down]->dcl->ihmm - th->hit[up]->dcl->jhmm - 1;
        if(graph->revcomp)
          seq_gap_len = th->hit[up]->dcl->jali - th->hit[down]->dcl->iali - 1;
        else
          seq_gap_len = th->hit[down]->dcl->iali - th->hit[up]->dcl->jali - 1;

        if(seq_gap_len < hmm_gap_len*3) continue;

       // printf("\nup %d down %d\n", up+1, down+1);
        gap = find_the_gap(graph, th, gm, target_seq, target_range->orig_N, up, down);
      //printf("seq start %d end %d hmm start %d end %d\n", gap->seq_start, gap->seq_end, gap->hmm_start, gap->hmm_end); 
       
        top_ten = align_the_gap(graph, th, gm, hmm, bg, target_seq, gcode, gap, &num_hits);
 
        prev_N = th->N;
        dup_hits = 0;
        for(hit = 0; hit < num_hits; hit++) {
          /* Add all new hits to target range and graph */
          add_missed_hit_to_target_range(target_range, top_ten[hit], &duplicate);
          if(!duplicate)
            add_missing_node_to_graph(graph, th, top_ten[hit], gm->M);
          dup_hits += duplicate;
        }
        num_hits -= dup_hits; 		

        if(num_hits) 
          new_edge = bridge_the_gap(graph, th, gm, hmm, bg, gcode, target_seq, prev_N, target_range->orig_N);         
       
        if(gap     != NULL) free(gap);
        if(top_ten != NULL) free(top_ten); 

        if(new_edge) {


        }
      }
    }
  }
  

  return eslOK;
}


int
check_for_loops (SPLICE_GRAPH *graph, P7_TOPHITS *th)
{

  int i,j;
  //int    *visited;
  //int    *rec_stack;
  //int     status; 

  for(i = 0; i < graph->num_nodes; i++) {
    for(j = 0; j < graph->num_nodes; j++) {
      if(i == j) continue;
      if(graph->edge_scores[i][j] != -eslINFINITY && graph->edge_scores[j][i] != -eslINFINITY) {
        if(graph->edge_scores[i][j] < graph->edge_scores[j][i]) { 
          graph->edge_scores[i][j] = -eslINFINITY;
          graph->out_edge_cnt[i]--;
          graph->in_edge_cnt[j]--;
        }
        else { 
          graph->edge_scores[j][i] = -eslINFINITY;
          graph->out_edge_cnt[j]--;
          graph->in_edge_cnt[i]--;
        }
      }
    }
  }

  return eslOK;
  
}



SPLICE_PATH*
evaluate_paths (SPLICE_GRAPH *graph, P7_TOPHITS *th, ESL_SQ *target_seq, int orig_N)
{

  int          i;
  int          path_len;
  int          start_node;
  int          prev_node;
  int          curr_node;
  int          next_node;
  int          contains_orig;
  int          path_ihmm, path_jhmm;
  int          edge_id;
  int          step_cnt;
  float        best_start_score;
  SPLICE_EDGE *in_edge;
  SPLICE_EDGE *out_edge;
  SPLICE_PATH *path;

  contains_orig = FALSE;

  /* Make sure we clear the field before we begin */
  for(i = 0; i < graph->num_nodes; i++) {
    graph->path_scores[i]   = -eslINFINITY;
    graph->best_out_edge[i] = -1;
  }

  /* Find best scoreing paths */ 
  longest_path_upstream(graph);

  while(!contains_orig) { 
    /* Find the best place to start our path */ 
    best_start_score = -eslINFINITY;
    start_node  = -1;
    for (i = 0; i < graph->num_nodes; i++) {
      if((!graph->in_edge_cnt[i]) && graph->path_scores[i] > best_start_score) {
        best_start_score = graph->path_scores[i];
  	  start_node  = i;
      }
    } 
  
    if(start_node < orig_N)
      contains_orig = TRUE;
      
    /* Get first edge in path (if any) and ckeck that it is not backward */ 
    curr_node = start_node;
    path_len  = 1;
    
    if (graph->out_edge_cnt[start_node]) {
      curr_node = graph->best_out_edge[start_node];
      edge_id   = graph->edge_id[start_node][curr_node];
      in_edge   = graph->edges[edge_id];
      prev_node = start_node;
      path_len++;
   
      path_ihmm = th->hit[prev_node]->dcl->ihmm; 
      if(path_ihmm >= in_edge->upstream_spliced_amino_end) {
     
        graph->edge_scores[prev_node][curr_node] = -eslINFINITY;
        graph->out_edge_cnt[prev_node]--;
        graph->in_edge_cnt[curr_node]--;
        /* If we found a backward node we need to start over */
        path = evaluate_paths (graph, th, target_seq, orig_N);
        return path;
      }
    }
  
    /* Get all other edged in path (if any) and ckeck that they are not backward */
    while(graph->out_edge_cnt[curr_node]) {
      if(curr_node < orig_N)
        contains_orig = TRUE;
  
      next_node = graph->best_out_edge[curr_node];
      edge_id   = graph->edge_id[curr_node][next_node];
      out_edge  = graph->edges[edge_id];
      /* If this is a backward node (in edge ends after out edge begins) 
       * delete the edge with the lower splice score */
     
      if((in_edge->downstream_spliced_amino_start >= out_edge->upstream_spliced_amino_end) || 
         (in_edge->downstream_spliced_nuc_start   >= out_edge->upstream_spliced_nuc_end)) { 
        
        if(out_edge->splice_score < in_edge->splice_score) {
      
          graph->edge_scores[curr_node][next_node] = -eslINFINITY;
          graph->out_edge_cnt[curr_node]--;
          graph->in_edge_cnt[next_node]--;
        }
        else {
          graph->edge_scores[prev_node][curr_node] = -eslINFINITY;
          graph->out_edge_cnt[prev_node]--;
          graph->in_edge_cnt[curr_node]--; 
        }
        /* If we found a backward node we need to start over */
        path = evaluate_paths (graph, th, target_seq, orig_N);
        return path; 
      }
      in_edge   = out_edge;
      prev_node = curr_node;
      curr_node = next_node;
      path_len++;
    }
  
    if(curr_node < orig_N)
      contains_orig = TRUE; 
  
    /* Check that last edge (if any) is not backward */
    if(path_len > 1) {
      path_jhmm = th->hit[curr_node]->dcl->jhmm;
      if(in_edge->downstream_spliced_amino_start >= path_jhmm) {
         
        graph->edge_scores[prev_node][curr_node] = -eslINFINITY;
        graph->out_edge_cnt[prev_node]--;
        graph->in_edge_cnt[curr_node]--;
        /* If we found a backward node we need to start over */
        path = evaluate_paths (graph, th, target_seq, orig_N);
        return path;
      }
    }

    if(!contains_orig) 
      graph->path_scores[start_node] = -eslINFINITY; 
  }    

  graph->best_path_start  = start_node;
  graph->best_path_end    = curr_node;
  graph->best_path_length = path_len;

  path = splice_path_create(path_len);
  path->revcomp = graph->revcomp;

  path->signal_scores[0] = 0.;

  path->node_id[0] = start_node;
  path->hits[0]    = th->hit[start_node];

  if(start_node < orig_N)
    path->missing[0] = FALSE;
  else
    path->missing[0] = TRUE;

  path->downstream_spliced_amino_start[0] = path->hits[0]->dcl->ihmm;
 
  /* Nuc coords need to be on the reverse strand */
  if(graph->revcomp) 
    path->downstream_spliced_nuc_start[0] = target_seq->start - path->hits[0]->dcl->iali + 1;
  else               
    path->downstream_spliced_nuc_start[0] = path->hits[0]->dcl->iali - target_seq->start + 1; 
  
  curr_node = start_node;
  step_cnt = 1;
  while (step_cnt < path_len) {
    next_node = graph->best_out_edge[curr_node]; 
    edge_id   = graph->edge_id[curr_node][next_node];
    out_edge  = graph->edges[edge_id];

    path->signal_scores[step_cnt] = out_edge->signal_score;
 
    path->node_id[step_cnt] = next_node;
    path->hits[step_cnt]    = th->hit[next_node];
 
    if(next_node < orig_N)
      path->missing[step_cnt] = FALSE;
    else
      path->missing[step_cnt] = TRUE;
 
    path->upstream_spliced_amino_end[step_cnt]     = out_edge->upstream_spliced_amino_end;
    path->downstream_spliced_amino_start[step_cnt] = out_edge->downstream_spliced_amino_start;
  
    path->upstream_spliced_nuc_end[step_cnt]     = out_edge->upstream_spliced_nuc_end;
    path->downstream_spliced_nuc_start[step_cnt] = out_edge->downstream_spliced_nuc_start;
     
    path->seq_len = path->seq_len + (path->upstream_spliced_nuc_end[step_cnt] - path->downstream_spliced_nuc_start[step_cnt-1] + 1); 
    
    curr_node = next_node;
    step_cnt++;
  }

  path->upstream_spliced_amino_end[step_cnt] = path->hits[step_cnt-1]->dcl->jhmm;  
  if(graph->revcomp)
    path->upstream_spliced_nuc_end[step_cnt] = target_seq->start - path->hits[step_cnt-1]->dcl->jali + 1;
  else
    path->upstream_spliced_nuc_end[step_cnt] = path->hits[step_cnt-1]->dcl->jali - target_seq->start + 1;  
  
  path->seq_len = path->seq_len + (path->upstream_spliced_nuc_end[step_cnt] - path->downstream_spliced_nuc_start[step_cnt-1] + 1);

//  path_dump(stdout, path, target_seq); 
  return path; 
}

int
longest_path_upstream (SPLICE_GRAPH *graph)
{

  int   i;
  int   node;
  int   stack_size;
  float step_score;
  int   *visited;
  int   *stack;
  int    status;

  /* Append source node */
  if((status = splice_graph_grow(graph)) != eslOK) goto ERROR;
  graph->hit_scores[graph->num_nodes]  = 0.;
  graph->edge_scores[graph->num_nodes][graph->num_nodes]  = -eslINFINITY;
  for(i = 0; i < graph->num_nodes; i++) {
    graph->edge_scores[graph->num_nodes][i] = -eslINFINITY;
    if(!graph->out_edge_cnt[i])
      graph->edge_scores[i][graph->num_nodes] = 0.;
    else 
      graph->edge_scores[i][graph->num_nodes] = -eslINFINITY; 
  } 
  graph->path_scores[graph->num_nodes] = 0.;
  graph->num_nodes++;


  ESL_ALLOC(visited, sizeof(int) * graph->num_nodes);
  esl_vec_ISet(visited,   graph->num_nodes, 0); 
 
  ESL_ALLOC(stack, sizeof(int) * graph->num_nodes);

  stack_size = 0;
  for(i = 0; i < graph->num_nodes; i++) {
    if(!visited[i])
      topological_sort_upstream(graph, visited, stack, &stack_size, i);
  } 


  while(stack_size > 0) {

    node = stack[stack_size-1];
    stack_size--;
    
    if(graph->path_scores[node] != -eslINFINITY) {
      for(i = 0; i < graph->num_nodes; i++) {
        if(graph->edge_scores[i][node] != -eslINFINITY) {
          step_score = graph->path_scores[node] + graph->edge_scores[i][node] + graph->hit_scores[i];
          if(graph->path_scores[i] <= step_score) {
            graph->path_scores[i]   = step_score;
            graph->best_out_edge[i] = node;
          }
        }
      }
    }
  }  

  /* Erase source node */
  graph->num_nodes--;
  if(visited  != NULL) free(visited);
  if(stack    != NULL) free(stack);
  return eslOK;
 
  ERROR:
    if(visited  != NULL) free(visited);
    if(stack    != NULL) free(stack);
    return status;    
}




int
topological_sort_upstream(SPLICE_GRAPH *graph, int *visited, int *stack, int *stack_size, int node)
{

  int i;
  
  visited[node] = TRUE;

  for(i = 0; i < graph->num_nodes; i++) {
    if(graph->edge_scores[i][node] != -eslINFINITY) {
      if(!visited[i])
         topological_sort_upstream(graph, visited, stack, stack_size, i);
    }
  }
  
  stack[*stack_size] = node;
  *stack_size += 1;
  
  return eslOK;

}


int
split_hits_in_path (SPLICE_GRAPH *graph, SPLICE_PATH *path, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, ESL_GENCODE *gcode, ESL_SQ *target_seq, int orig_N)
{

  int i, z;
  int z1, z2;
  int I_cnt;
  int I_start, I_end;
  P7_HIT *curr_hit;
  P7_TRACE *tr;
  SPLICE_EDGE *edge;
  int status;


  for(i = 0; i < path->path_len; i++) {

    curr_hit = path->hits[i];
    tr = curr_hit->dcl->tr;

    I_cnt = 0;

    /* find the position in trace that coorespond to the path splice boudries*/
    for (z1 = 0;       z1 < tr->N; z1++) if (tr->k[z1] >= path->downstream_spliced_amino_start[i]) break;
    for (z2 = tr->N-1; z2 >= 0;    z2--) if (tr->k[z2] <= path->upstream_spliced_amino_end[i+1]) break;

    for (z = z1; z < z2; z++) {
      if(tr->st[z] == p7T_I) {
        if(I_cnt == 0) I_start = z;
        I_cnt++;
      }
      else {
        if(I_cnt > 9) {
          I_end = z-1;

          edge = NULL;
          edge = splice_edge_create();

          edge->overlap_amino_start = ESL_MAX(path->downstream_spliced_amino_start[i], tr->k[I_start]-MIN_AMINO_OVERLAP);
          edge->overlap_amino_end   = ESL_MIN(path->upstream_spliced_amino_end[i+1],   tr->k[I_end]+MIN_AMINO_OVERLAP);

          get_overlap_nuc_coords(edge, curr_hit->dcl, curr_hit->dcl, target_seq, path->revcomp);

          /* Add extra nucleotides for splice sites */
          if(path->revcomp) {
            edge->upstream_nuc_end     += 2;
            edge->downstream_nuc_start -= 2;
          }
          else {
            edge->upstream_nuc_end     += 2;
            edge->downstream_nuc_start -= 2;
          }


          if ((status = find_optimal_splice_site (edge, curr_hit->dcl, curr_hit->dcl, gm, hmm, bg, gcode, target_seq)) != eslOK) goto ERROR;

          /* Make sure we are not attempting to replice the same insertion */
          if(edge->upstream_spliced_nuc_end <= path->downstream_spliced_nuc_start[i]) {
            if (edge != NULL) free(edge);
            I_cnt = 0;
            continue;
          }

          if (graph->hit_scores[path->node_id[i]] + edge->splice_score > 0) {
            if ((status = splice_path_split_hit(path, edge, i)) != eslOK) goto ERROR;
            if (edge != NULL) free(edge);
            I_cnt = 0;
            break;
          }

          if (edge != NULL) free(edge);

        }

        I_cnt = 0;
      }

    }
  }
  return eslOK;

  ERROR:
    if(edge != NULL) free(edge);
    return status;
}



int
splice_path (SPLICE_GRAPH *graph, SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_OPROFILE *om, P7_SCOREDATA *scoredata, ESL_SQ *target_seq, ESL_GENCODE *gcode, int64_t db_nuc_cnt, int orig_N, int* success) 
{

  int          i;
  int          pos;
  int          exon;
  int          shift;
  int          seq_idx;
  int          amino_len;
  int          amino;
  int          env_len;
  int          orf_len;
  int          replace_node;
  int          remove_node;
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
  int          status;

  nuc_index = NULL;
  nuc_dsq   = NULL;

  ESL_ALLOC(nuc_index, sizeof(int64_t) * (path->seq_len+2));  
  ESL_ALLOC(nuc_dsq,   sizeof(ESL_DSQ) * (path->seq_len+2));

  /* Copy spliced nucleotides into single sequence and track their original indicies */
  nuc_index[0] = -1;
  nuc_dsq[0]   = eslDSQ_SENTINEL;
  seq_idx   = 1;
  for (i = 0; i < path->path_len; i++) {
    for (pos = path->downstream_spliced_nuc_start[i]; pos <= path->upstream_spliced_nuc_end[i+1]; pos++) {
      nuc_index[seq_idx] = pos;
      nuc_dsq[seq_idx]   = target_seq->dsq[pos];
      seq_idx++;    
    }
  }

  nuc_index[seq_idx] = -1;
  nuc_dsq[seq_idx]   = eslDSQ_SENTINEL;


  /* Translate spliced nucleotide sequence */
  amino_len = path->seq_len/3;
  ESL_ALLOC(amino_dsq, sizeof(ESL_DSQ) * (amino_len+2));

  amino_dsq[0] = eslDSQ_SENTINEL;
   
  seq_idx = 1;

  for(i = 1; i <= amino_len; i++) {
    amino = esl_gencode_GetTranslation(gcode,&nuc_dsq[seq_idx]);
    amino_dsq[i] = amino;
    seq_idx+=3;
  }
  amino_dsq[amino_len+1] = eslDSQ_SENTINEL;
 
  /* Create ESL_SQs from ESL_DSQs */
  amino_seq   = esl_sq_CreateDigitalFrom(gcode->nt_abc, target_seq->name, amino_dsq, amino_len,      NULL,NULL,NULL); 
  nuc_seq     = esl_sq_CreateDigitalFrom(gcode->aa_abc, target_seq->name, nuc_dsq,   path->seq_len,  NULL,NULL,NULL);

  pli->nuc_sq       = nuc_seq;
  pli->amino_sq     = amino_seq;   
  pli->orig_nuc_idx = nuc_index;

  /* Algin the splices Amino sequence */
  status = align_spliced_path(pli, om, scoredata, target_seq, gcode);
 
  /* Alignment failed */
  if(pli->hit == NULL ) { 
 
    if(nuc_dsq   != NULL) free(nuc_dsq);
    if(amino_dsq != NULL) free(amino_dsq);
 
    *success = FALSE; 
     return eslOK; 
  }
 
  /* adjust all coords in hit and path */
  if(path->revcomp) { 
    pli->hit->dcl->ad->sqfrom += 2; 
    pli->hit->dcl->ienv = target_seq->n - pli->orig_nuc_idx[1]              + target_seq->end; 
    pli->hit->dcl->jenv = target_seq->n - pli->orig_nuc_idx[pli->nuc_sq->n] + target_seq->end; 
    env_len = (pli->hit->dcl->ienv - pli->hit->dcl->jenv + 1)/3.; 
  }  
  else {
    pli->hit->dcl->ienv = pli->orig_nuc_idx[1]              + target_seq->start -1;
    pli->hit->dcl->jenv = pli->orig_nuc_idx[pli->nuc_sq->n] + target_seq->start -1;
    env_len = (pli->hit->dcl->jenv - pli->hit->dcl->ienv + 1)/3.;
  }
  env_len = ESL_MAX(env_len, om->max_length);
  //printf(" pli->hit->dcl->iali %d  pli->hit->dcl->jali %d \n", pli->hit->dcl->iali,  pli->hit->dcl->jali);
  if ( path->path_len > 1) {
    /* Shift the path to start at the first hit that was inculded in the alignment 
     * and end at the last hit that was included in the alignment */ 
    for(exon = 0; exon < path->path_len; exon++) {
      if(path->upstream_spliced_nuc_end[exon+1] > pli->hit->dcl->iali)  
        break;
    }
    /* Shift path to start at frist hits that is in alignment */
    path->path_len -= exon;
    
    for(shift = 0; shift < exon; shift++) {
      path->node_id[shift]                        = path->node_id[shift+exon]; 
      path->upstream_spliced_amino_end[shift]     = path->upstream_spliced_amino_end[shift+exon];
      path->downstream_spliced_amino_start[shift] = path->downstream_spliced_amino_start[shift+exon];
      path->upstream_spliced_nuc_end[shift]       = path->upstream_spliced_nuc_end[shift+exon];   
      path->downstream_spliced_nuc_start[shift]   = path->downstream_spliced_nuc_start[shift+exon];
      path->hit_scores[shift]                     = path->hit_scores[shift+exon];
      path->signal_scores[shift]                  = path->signal_scores[shift+exon];
      path->hits[shift]                           = path->hits[shift+exon];
      path->missing[shift]                        = path->missing[shift+exon];
    }
    path->downstream_spliced_nuc_start[0]   = pli->hit->dcl->iali;
    path->downstream_spliced_amino_start[0] = pli->hit->dcl->ihmm; 

    for(exon = 1; exon < path->path_len; exon++) {
      if (path->downstream_spliced_nuc_start[exon] > pli->hit->dcl->jali)
        break;
    }
   
    path->path_len -= (path->path_len-exon);
    
    path->upstream_spliced_nuc_end[exon]   = pli->hit->dcl->jali;
    path->upstream_spliced_amino_end[exon] = pli->hit->dcl->jhmm;   
  }  

  pli->hit->dcl->iali = pli->hit->dcl->ad->sqfrom;
  pli->hit->dcl->jali = pli->hit->dcl->ad->sqto;


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
    /* Replace the first hit in path with the spliced hit 
     * and set all other hits in path to unreportable */ 
    i = 0;
    while( path->node_id[i] >= orig_N) {
      pli->hit->dcl->ad->exon_orig[i] = FALSE; 
      i++; 
    }

    pli->hit->dcl->ad->exon_orig[i] = TRUE;
    replace_hit = path->hits[i]; 
    replace_node = path->node_id[i];
	ali_L = replace_hit->dcl->ad->L;
    p7_domain_Destroy(replace_hit->dcl);

    replace_hit->dcl        = pli->hit->dcl;
	replace_hit->dcl->ad->L = ali_L;
 
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
  
    replace_hit->frameshift = FALSE;

    /* Set all other hits in alignment to unreportable */ 
    i++;
    for(  ; i < path->path_len; i++) {
      if(path->node_id[i] >= orig_N)
        pli->hit->dcl->ad->exon_orig[i] = FALSE;
      else
        pli->hit->dcl->ad->exon_orig[i] = TRUE;
     
      /* If the replace node hes bee split we neeed to 
       * make sure we don't set is as unreportable */
      remove_node = path->node_id[i];
      if(remove_node == replace_node)
         continue;

      remove_hit = path->hits[i];

      if(remove_hit->flags & p7_IS_REPORTED ) {
        remove_hit->flags &= ~p7_IS_REPORTED;
        remove_hit->dcl->is_reported = FALSE;
      }
      if((remove_hit->flags & p7_IS_INCLUDED)) {
        remove_hit->flags &= ~p7_IS_INCLUDED;
        remove_hit->dcl->is_included = FALSE;
      }
    }
    
  } 
  else *success = FALSE; //printf("EVAL fail \n"); }

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
align_spliced_path (SPLICE_PIPELINE *pli, P7_OPROFILE *om, P7_SCOREDATA *scoredata, ESL_SQ *target_seq, ESL_GENCODE *gcode) 
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
  
  if (P > pli->F3) return eslOK;

  p7_Backward(pli->amino_sq->dsq, pli->amino_sq->n, om, pli->fwd, pli->bwd, NULL);  

  if((p7_Decoding(om, pli->fwd, pli->bwd, pli->bwd)) == eslERANGE) 
    ESL_XEXCEPTION(eslFAIL, "p7_Decoding failed with eslERANGE");

  p7_OptimalAccuracy(om, pli->bwd, pli->fwd, &oasc);
  p7_OATrace        (om, pli->bwd, pli->fwd, tr);
  
  p7_trace_Index(tr);
   
  if (scoredata->prefix_lengths == NULL)
    p7_hmm_ScoreDataComputeRest(om, scoredata);
 
  compute_ali_scores(hit->dcl, tr, pli->amino_sq->dsq, scoredata, om->abc->Kp);
 
  if((hit->dcl->tr = p7_trace_splice_Convert(tr, pli->orig_nuc_idx, &splice_cnt)) == NULL) goto ERROR; 
  
  if((hit->dcl->ad = p7_alidisplay_splice_Create(hit->dcl->tr, 0, om, target_seq, pli->amino_sq, hit->dcl->scores_per_pos, tr->sqfrom[0], splice_cnt)) == NULL) goto ERROR; 
  
  p7_Null2_ByExpectation(om, pli->bwd, null2);
  domcorrection = 0.;
  for (i = 0; i < pli->amino_sq->n; i++)
    domcorrection += logf(null2[pli->amino_sq->dsq[i]]);
 
  hit->dcl->domcorrection = domcorrection;
 
  hit->dcl->ihmm = hit->dcl->ad->hmmfrom;
  hit->dcl->jhmm = hit->dcl->ad->hmmto;
  /*Keep iali, jali in sub-sequence teaget_seq coords for now */
  hit->dcl->iali = hit->dcl->ad->sqfrom;
  hit->dcl->jali = hit->dcl->ad->sqto;

  /* Convert sqfrom, sqto to full sequence coords */
  if(target_seq->start < target_seq->end) {
    hit->dcl->ad->sqfrom  = hit->dcl->ad->sqfrom + target_seq->start - 1;
    hit->dcl->ad->sqto    = hit->dcl->ad->sqto   + target_seq->start - 1;
  } else {
    hit->dcl->ad->sqto    = target_seq->n - hit->dcl->ad->sqto   + target_seq->end;
    hit->dcl->ad->sqfrom  = target_seq->n - hit->dcl->ad->sqfrom + target_seq->end;
  }
 
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

SPLICE_GAP*
find_the_gap (SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm, ESL_SQ *target_seq, int orig_N, int upstream_node, int downstream_node) 
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
  //printf("us_hit->dcl->ihmm %d us_hit->dcl->jhmm %d us_hit->dcl->iali %d us_hit->dcl->jali %d\n", us_hit->dcl->ihmm, us_hit->dcl->jhmm, us_hit->dcl->iali, us_hit->dcl->jali);
  //printf("ds_hit->dcl->ihmm %d ds_hit->dcl->jhmm %d ds_hit->dcl->iali %d ds_hit->dcl->jali %d\n", ds_hit->dcl->ihmm, ds_hit->dcl->jhmm, ds_hit->dcl->iali, ds_hit->dcl->jali);
  gap->hmm_start = ESL_MAX(us_hit->dcl->jhmm-hmm_ext, 1);
  gap->hmm_end   = ESL_MIN(ds_hit->dcl->ihmm+hmm_ext, gm->M);
  
  if(graph->revcomp) {
    gap->seq_start = ESL_MIN(us_hit->dcl->jali+(hmm_ext*3+2), target_seq->start);
    gap->seq_end   = ESL_MAX(ds_hit->dcl->iali-(hmm_ext*3),   target_seq->end);
  }
  else {
    gap->seq_start = ESL_MAX(us_hit->dcl->jali-(hmm_ext*3+2), target_seq->start);
    gap->seq_end   = ESL_MIN(ds_hit->dcl->iali+(hmm_ext*3), target_seq->end);
  }
  //printf("revcomp %d hmm_start %d hmm_end %d seq_start %d seq_end %d\n", graph->revcomp, gap->hmm_start, gap->hmm_end, gap->seq_start, gap->seq_end);
  gap->upstream_node   = upstream_node;
  gap->downstream_node = downstream_node;
  
  return gap;

  ERROR:
    if(gap != NULL) free(gap);
	return NULL;

}

SPLICE_GAP*
terminal_gap (SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm, ESL_SQ *target_seq, int orig_N, int final_node, int path_node, int look_upstream) 
{

  P7_HIT *path_hit;
  P7_HIT *final_hit;
  SPLICE_GAP *gap;
  int     status; 

  gap = NULL;
  ESL_ALLOC(gap, sizeof(SPLICE_GAP));

  path_hit  = th->hit[path_node];
  final_hit = th->hit[final_node];

  if(look_upstream) {
  
    gap->hmm_start = 1; 
    gap->hmm_end   = ESL_MIN(path_hit->dcl->ihmm+2, gm->M);
  
    if(graph->revcomp) {
      gap->seq_start = ESL_MIN(final_hit->dcl->iali+TERM_RANGE_EXT, target_seq->start);
      gap->seq_end   = ESL_MAX(final_hit->dcl->iali-6,   target_seq->end);
    }
    else {
      gap->seq_start = ESL_MAX(final_hit->dcl->iali-TERM_RANGE_EXT, target_seq->start);
      gap->seq_end   = ESL_MIN(final_hit->dcl->iali+6, target_seq->end);
    }
  }
  else {
    gap->hmm_start = ESL_MAX(path_hit->dcl->jhmm-2, 1);
    gap->hmm_end   = gm->M; 
 
    if(graph->revcomp) {
      gap->seq_start = ESL_MIN(final_hit->dcl->iali+6, target_seq->start);
      gap->seq_end   = ESL_MAX(final_hit->dcl->iali-TERM_RANGE_EXT,   target_seq->end);
    }
    else {
      gap->seq_start = ESL_MAX(final_hit->dcl->jali-6, target_seq->start);
      gap->seq_end   = ESL_MIN(final_hit->dcl->jali+TERM_RANGE_EXT, target_seq->end);
    }
  }
  gap->upstream_node   = path_node;
  gap->downstream_node = path_node;
  
  return gap;

  ERROR:
    if(gap != NULL) free(gap);
	return NULL;

}






P7_HIT**
align_the_gap(SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, ESL_SQ *target_seq, ESL_GENCODE *gcode, SPLICE_GAP *gap, int *num_hits)
{
 
  int         i, j, z;
  int         z1,z2;
  int         orf_min;
  int         seq_sq_len;
  int         num_saved;
  int         low_idx;
  float       vitsc;
  float       seqsc;
  float       low_savesc;
  char         *spoofline;
  P7_HMM       *sub_hmm;
  P7_PROFILE   *sub_model;
  P7_OPROFILE  *sub_omodel;
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

  top_ten = NULL;
  ESL_ALLOC(top_ten, sizeof(P7_HIT*) * 10);
  
  sub_hmm    = extract_sub_hmm(hmm, gap->hmm_start, gap->hmm_end);  
  sub_model  = p7_profile_Create (sub_hmm->M, sub_hmm->abc); 
  sub_omodel = p7_oprofile_Create(sub_hmm->M, sub_hmm->abc );
  p7_ProfileConfig(sub_hmm, bg, sub_model, 100, p7_LOCAL);
  p7_oprofile_Convert(sub_model, sub_omodel);

  scoredata  = p7_hmm_ScoreDataCreate(sub_omodel, NULL);
  p7_hmm_ScoreDataComputeRest(sub_omodel, scoredata);

  sub_dsq = extract_sub_seq(target_seq, gap->seq_start, gap->seq_end, graph->revcomp);
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

    if(graph->revcomp) {
      orf_hit->dcl->iali = sub_sq->start - sub_sq->n + orf_sq->start - 3*(orf_hit->dcl->tr->i[z1]-1);
      orf_hit->dcl->jali = sub_sq->start - sub_sq->n + orf_sq->start - 3*(orf_hit->dcl->tr->i[z2])+1;
    }
    else {
      orf_hit->dcl->iali = sub_sq->start + orf_sq->start + (orf_hit->dcl->tr->i[z1]*3-2) - 2;
      orf_hit->dcl->jali = sub_sq->start + orf_sq->start + (orf_hit->dcl->tr->i[z2]*3)   - 2;
    }
    
    orf_hit->dcl->ad = p7_alidisplay_Create(orf_hit->dcl->tr, 0, sub_omodel, orf_sq, NULL);
  
    p7_pli_computeAliScores(orf_hit->dcl, orf_sq->dsq, scoredata, sub_omodel->abc->Kp, TRUE); 

    seqsc = 0.;
    for(j = 0; j < orf_hit->dcl->ad->N; j++)
      seqsc += orf_hit->dcl->scores_per_pos[j];

    if(seqsc > low_savesc) {

      orf_hit->score = seqsc;

      if(num_saved < 10) 
        top_ten[num_saved] = orf_hit;
      else {
        low_savesc =  top_ten[0]->score;
        low_idx = 0;
        for(j = 1; j < 10; j++) {
          if(top_ten[j]->score < low_savesc) {
             low_savesc  = top_ten[j]->score;
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
      if(low_savesc == -eslINFINITY || seqsc < low_savesc) low_savesc = seqsc;
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


int
compute_ali_scores(P7_DOMAIN *dom, P7_TRACE *tr, ESL_DSQ *amino_dsq, const P7_SCOREDATA *data, int K)
{

  int status;
  int i, j, k;
  int z;
  int N;

  N = tr->tto[0] - tr->tfrom[0] - 1;
  ESL_ALLOC( dom->scores_per_pos, sizeof(float) * N);

  for (i=0; i<N; i++)  dom->scores_per_pos[i] = 0.0;
  i = tr->sqfrom[0] - 1;
  j = tr->hmmfrom[0] - 1;
  k = 0;
  z = tr->tfrom[0]+1;
  while (k<N) {
    if (tr->st[z] == p7T_M) {
      i++;  j++;
      dom->scores_per_pos[k] = data->fwd_scores[K * j + amino_dsq[i]]
                             +  (j==1 ? 0 : log(data->fwd_transitions[p7O_MM][j]) );
      k++; z++;
    }
    else if (tr->st[z] == p7T_I) {
      dom->scores_per_pos[k] = log(data->fwd_transitions[p7O_MI][j]);
      i++; k++; z++;
      while (k<N && tr->st[z] == p7T_I) {
        dom->scores_per_pos[k] = log(data->fwd_transitions[p7O_II][j]);
        i++; k++; z++;
      }
     }
     else if (tr->st[z] == p7T_D) {
       dom->scores_per_pos[k] = log(data->fwd_transitions[p7O_DD][j]);
       j++; k++; z++;
       while (k<N && tr->st[z] == p7T_D)  {
         dom->scores_per_pos[k] = log(data->fwd_transitions[p7O_DD][j]);
         j++; k++; z++;
       } 
     }
     else ESL_XEXCEPTION(eslFAIL, "Impossible state from compute_ali_scores()");
  }

  return eslOK;

  ERROR:
    return status;
}


P7_HMM*
extract_sub_hmm (P7_HMM *hmm, int start, int end) 
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

  int   i;
  int   status;

  /* Make sure we have room for our new node */
  if((status = splice_graph_grow(graph)) != eslOK) goto ERROR;

  graph->hit_scores[graph->num_nodes]     = hit->score; //This contiains the summed alisc
  graph->path_scores[graph->num_nodes] = -eslINFINITY;
 
  graph->out_edge_cnt[graph->num_nodes]  = 0;
  graph->in_edge_cnt[graph->num_nodes]   = 0;
  graph->orig_out_edge[graph->num_nodes] = 0;
  graph->orig_in_edge[graph->num_nodes]  = 0;
  graph->best_out_edge[graph->num_nodes] = -1;
  graph->best_in_edge[graph->num_nodes]  = -1;

  for(i = 0; i < graph->num_nodes; i++) {
    graph->edge_scores[i][graph->num_nodes] = -eslINFINITY;
    graph->edge_scores[graph->num_nodes][i] = -eslINFINITY;
    graph->edge_id[i][graph->num_nodes]     = -1;
    graph->edge_id[graph->num_nodes][i]     = -1;
  }

  graph->edge_scores[graph->num_nodes][graph->num_nodes] = -eslINFINITY;
  graph->edge_id[graph->num_nodes][graph->num_nodes]     = -1;
  graph->is_upstream[graph->num_nodes][graph->num_nodes] = FALSE;
  
  for(i  = 0; i < graph->num_nodes; i++) {
    graph->is_upstream[i][graph->num_nodes] = TRUE; 
    graph->is_upstream[graph->num_nodes][i] = TRUE;
    
    if(th->hit[i]->dcl->ihmm >= hit->dcl->ihmm &&
       th->hit[i]->dcl->jhmm >= hit->dcl->jhmm) {
      graph->is_upstream[i][graph->num_nodes] = FALSE;
    }

    if(hit->dcl->ihmm >= th->hit[i]->dcl->ihmm &&
       hit->dcl->jhmm >= th->hit[i]->dcl->jhmm) {
      graph->is_upstream[graph->num_nodes][i] = FALSE;
    }

    if (( graph->revcomp  && hit->dcl->jali > th->hit[i]->dcl->iali) ||
       ((!graph->revcomp) && hit->dcl->jali < th->hit[i]->dcl->iali)) {
       graph->is_upstream[i][graph->num_nodes] = FALSE;
    } 

    if (( graph->revcomp  && th->hit[i]->dcl->jali > hit->dcl->iali) ||
       ((!graph->revcomp) && th->hit[i]->dcl->jali < hit->dcl->iali)) {
       graph->is_upstream[graph->num_nodes][i] = FALSE;
    }
    //printf("i %d graph->num_nodes %d i>n %d n>i %d\n", i+1, graph->num_nodes+1,  graph->is_upstream[i][graph->num_nodes], graph->is_upstream[graph->num_nodes][i]); 
  }
    

  graph->num_nodes++;

  return eslOK; 

  ERROR:
    return status;
}

int
add_missed_hit_to_target_range(TARGET_RANGE *target_range, P7_HIT *hit, int *duplicate)
{

  int        h;
  int        start, end;
  P7_TOPHITS *th;
  int        status;

  th = target_range->th;

  *duplicate = FALSE;
  for(h = 0; h < target_range->th->N; h++) {
    if(hit->dcl->ihmm >= th->hit[h]->dcl->ihmm && hit->dcl->jhmm <= th->hit[h]->dcl->jhmm) {
      if ((target_range->revcomp   && (hit->dcl->iali <= th->hit[h]->dcl->iali && hit->dcl->jali >= th->hit[h]->dcl->jali)) ||
         ((!target_range->revcomp) && (hit->dcl->iali >= th->hit[h]->dcl->iali && hit->dcl->jali <= th->hit[h]->dcl->jali)))   {
        *duplicate = TRUE;        
        return eslOK;
      }
    }
  }
  
  if((status = target_range_grow(target_range)) != eslOK) goto ERROR;
  
  th->hit[th->N] = hit;
  th->N++;

  start = ESL_MIN(hit->dcl->iali, hit->dcl->jali);
  end   = ESL_MAX(hit->dcl->iali, hit->dcl->jali);

  target_range->start = ESL_MIN(start, target_range->start);
  target_range->end   = ESL_MAX(end, target_range->end);

  return eslOK;

  ERROR:
    return status;
}


int
bridge_the_gap(SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, ESL_GENCODE *gcode, ESL_SQ *target_seq, int prev_N, int orig_N)
{
  int         i;
  int         up, down;
  int         new_hits;
  int         new_edge;
  int         curr_num_connected;
  int         next_num_connected;
  int         *curr_connected_hits;
  int         *next_connected_hits;
  SPLICE_EDGE *edge;
  int         status;

  new_hits = th->N - orig_N;

  curr_connected_hits = NULL;
  next_connected_hits = NULL;
  ESL_ALLOC(curr_connected_hits, sizeof(int) * new_hits); 
  ESL_ALLOC(next_connected_hits, sizeof(int) * new_hits);
  curr_num_connected = 0;
  next_num_connected = 0;
 
  for(i = 0; i < new_hits; i++) {
    curr_connected_hits[i] = FALSE;
    next_connected_hits[i] = FALSE;
  }
 
  edge = NULL;
  new_edge = FALSE;
 

  /* Connect new hits upstream of original hits */
  for(up = prev_N; up < th->N; up++) {
    for(down = 0; down < orig_N; down++) {

      if(graph->is_upstream[up][down]) {
        edge = connect_nodes_with_edges(th->hit[up], th->hit[down], gm, hmm, bg, gcode, target_seq, graph->revcomp);
        if(edge != NULL) {
          if(ESL_MIN(graph->hit_scores[up], graph->hit_scores[down]) + edge->splice_score > 0) {
            new_edge = TRUE;
            edge->upstream_node_id   = up;
            edge->downstream_node_id = down;
            add_edge_to_graph(graph, edge);

            curr_connected_hits[up-prev_N] = TRUE;
            curr_num_connected++;
          }
          else free(edge); 
        }
      }   
    }
  }  

  /* Connect new hits that connected to original hits downstream of any non-original hits */
  for(up = orig_N; up < th->N; up++) {
    for(down = prev_N; down < th->N; down++) {
      if(curr_connected_hits[down-prev_N] && graph->is_upstream[up][down] && graph->edge_id[up][down] == -1) {

        edge = connect_nodes_with_edges(th->hit[up], th->hit[down], gm, hmm, bg, gcode, target_seq, graph->revcomp);

        if(edge != NULL) {
          if(ESL_MIN(graph->hit_scores[up], graph->hit_scores[down])  + edge->splice_score > 0) {
            edge->upstream_node_id   = up;
            edge->downstream_node_id = down;
            add_edge_to_graph(graph, edge);
            next_connected_hits[up-orig_N] = TRUE;
            next_num_connected++;
          }
          else free(edge);
        }
      }
    }
  }
   
  free(curr_connected_hits);
  curr_connected_hits = next_connected_hits;
  curr_num_connected  = next_num_connected;
  next_connected_hits = NULL;
  ESL_ALLOC(next_connected_hits, sizeof(int) * new_hits);
  for(i = 0; i < new_hits; i++)
    next_connected_hits[i] = FALSE;
  next_num_connected = 0;

  /* Continue connecting non-original his upstream until we run out of connections */ 
  while(curr_num_connected > 0) {
   
    for(up = orig_N; up < th->N; up++) {
      for(down = orig_N; down < th->N; down++) {
        if(curr_connected_hits[down-orig_N] && graph->is_upstream[up][down] && graph->edge_id[up][down] == -1) {
    
          edge = connect_nodes_with_edges(th->hit[up], th->hit[down], gm, hmm, bg, gcode, target_seq, graph->revcomp);
          if(edge != NULL) {
            if(ESL_MIN(graph->hit_scores[up], graph->hit_scores[down])  + edge->splice_score > 0) {
              edge->upstream_node_id   = up;
              edge->downstream_node_id = down;
              add_edge_to_graph(graph, edge);
              next_connected_hits[up-orig_N] = TRUE;
              next_num_connected++;
            }
            else free(edge);
          }
        }
      }
    }
    free(curr_connected_hits);
    curr_connected_hits = next_connected_hits;
    curr_num_connected  = next_num_connected;
    next_connected_hits = NULL;
    ESL_ALLOC(next_connected_hits, sizeof(int) * new_hits);
    for(i = 0; i < new_hits; i++)
      next_connected_hits[i] = FALSE;
    next_num_connected = 0;
  } 

  for(i = 0; i < new_hits; i++) {
    curr_connected_hits[i] = FALSE;
    next_connected_hits[i] = FALSE;
  }

  /* Connect new hits downstream of original hits */
  for(up = 0; up < orig_N; up++) {
    for(down = prev_N; down < th->N; down++) { 
      if(graph->is_upstream[up][down]) {
        //printf("up %d down %d\n", up+1, down+1);        
        
        edge = connect_nodes_with_edges(th->hit[up], th->hit[down], gm, hmm, bg, gcode, target_seq, graph->revcomp);
        if(edge != NULL) {
          // printf("graph->hit_scores[up] %f graph->hit_scores[down] %f edge->splice_score %f\n", graph->hit_scores[up], graph->hit_scores[down], edge->splice_score);
          if(ESL_MIN(graph->hit_scores[up], graph->hit_scores[down]) + edge->splice_score > 0) {

          //printf("up %d down %d\n", up+1, down+1); 
            new_edge = TRUE;
            edge->upstream_node_id   = up;
            edge->downstream_node_id = down;
            add_edge_to_graph(graph, edge);

            curr_connected_hits[down-prev_N] = TRUE;
            curr_num_connected++;
          }
          else free(edge);
        }
      }
    }
  }

  
 /* Connect new hits that connected to original hits upstream of any non-original hits */
  for(up = prev_N; up < th->N; up++) {
    for(down = orig_N; down < th->N; down++) {
      if(curr_connected_hits[up-prev_N] && graph->is_upstream[up][down] && graph->edge_id[up][down] == -1) {
	
        edge = connect_nodes_with_edges(th->hit[up], th->hit[down], gm, hmm, bg, gcode, target_seq, graph->revcomp);

        if(edge != NULL) {
          if(ESL_MIN(graph->hit_scores[up], graph->hit_scores[down])  + edge->splice_score > 0) {
            edge->upstream_node_id   = up;
            edge->downstream_node_id = down;
            add_edge_to_graph(graph, edge);
            next_connected_hits[down-orig_N] = TRUE;
            next_num_connected++;
          }
          else free(edge);
        }
      }
    }
  }
   
  free(curr_connected_hits);
  curr_connected_hits = next_connected_hits;
  curr_num_connected  = next_num_connected;
  next_connected_hits = NULL;
  ESL_ALLOC(next_connected_hits, sizeof(int) * new_hits);
  for(i = 0; i < new_hits; i++)
    next_connected_hits[i] = FALSE;
  next_num_connected = 0;
  
  /* Continue connecting non-original hits downstream until we run out of connections */ 
  while(curr_num_connected > 0) {
  
    for(up = orig_N; up < th->N; up++) {
      for(down = orig_N; down < th->N; down++) {
        if(curr_connected_hits[up-orig_N] && graph->is_upstream[up][down] && graph->edge_id[up][down] == -1) {
			
          edge = connect_nodes_with_edges(th->hit[up], th->hit[down], gm, hmm, bg, gcode, target_seq, graph->revcomp);
          if(edge != NULL) {
            if(ESL_MIN(graph->hit_scores[up], graph->hit_scores[down]) + edge->splice_score > 0) {
              edge->upstream_node_id   = up;
              edge->downstream_node_id = down;
              add_edge_to_graph(graph, edge);
              next_connected_hits[down-orig_N] = TRUE;
              next_num_connected++;
            }
            else free(edge);
          }
        }
      }
    }    
    free(curr_connected_hits);
    curr_connected_hits = next_connected_hits;
    curr_num_connected  = next_num_connected;
    next_connected_hits = NULL;
    ESL_ALLOC(next_connected_hits, sizeof(int) * new_hits);
    for(i = 0; i < new_hits; i++) 
      next_connected_hits[i] = FALSE;
    next_num_connected = 0;
  }
    
  if(curr_connected_hits != NULL) free(curr_connected_hits);
  if(next_connected_hits != NULL) free(next_connected_hits);

  return new_edge; 

  ERROR:
    if(curr_connected_hits != NULL) free(curr_connected_hits);
    if(next_connected_hits != NULL) free(next_connected_hits);
    return FALSE;
}


void 
target_range_dump(FILE *fp, TARGET_RANGE *target_range, int print_hits) {

  int i;
  P7_HIT *hit = NULL;

  if (target_range == NULL) { fprintf(fp, " [ target range is NULL ]\n"); return; }

  fprintf(fp, " Target Range Sequence Name      : %s\n", target_range->seqname);
  fprintf(fp, " Target Range Reverse Complement : %s\n", (target_range->revcomp ? "YES" : "NO") );
  fprintf(fp, " Target Range Sequence Start     : %" PRId64 "\n", target_range->start);
  fprintf(fp, " Target Range Sequence End       : %" PRId64 "\n", target_range->end);
  fprintf(fp, " Target Range Hit Count          : %ld\n", target_range->th->N);

  if(print_hits) {
    fprintf(fp, "\n   Hit   hmm_from   hmm_to   seq_from     seq_to   score\n");

    for(i = 0; i < target_range->th->N; i++) {

      hit = target_range->th->hit[i];     

   	  fprintf(fp, "   %d  %6d  %9d %9" PRId64 " %12" PRId64 " %10.2f\n", 
                      i+1, hit->dcl->ihmm, hit->dcl->jhmm, hit->dcl->iali, hit->dcl->jali, hit->score);	 
   }
  }
  
  fprintf(fp, "\n"); 

  return;
}

void
graph_dump(FILE *fp, SPLICE_GRAPH *graph, ESL_SQ *target_seq, int print_edges) 
{


  int          i,j;
  int          edge_id;
  int          nuc_end, nuc_start;
  SPLICE_EDGE *edge;

  if (graph == NULL) { fprintf(fp, " [ graph is NULL ]\n"); return; }
  
  fprintf(fp, " SPLICE_GRAPH\n");
  fprintf(fp, " Number of Nodes  : %d\n", graph->num_nodes);
  fprintf(fp, " Number of Edges  : %d\n", graph->num_edges); 

  if(graph->num_edges > 0) {
    fprintf(fp, "\n Splice Score matrix \n");
    fprintf(fp, "     ");
    for(i = 0; i < graph->num_nodes; i++ ) 
      fprintf(fp, "%8d", i+1);
  
    fprintf(fp, "\n");

    for(i = 0; i < graph->num_nodes; i++ ) {

      fprintf(fp, "%5d", i+1);
      for(j = 0; j < graph->num_nodes; j++ ) {   
        fprintf(fp, "%8.2f", graph->edge_scores[i][j]);
      }
      fprintf(fp, "\n");
    }    
 }

  fprintf(fp, "\n Hit Scores  \n");
  for(i = 0; i < graph->num_nodes; i++ ) 
    fprintf(fp, "%8d", i+1);
  
  fprintf(fp, "\n");

  for(i = 0; i < graph->num_nodes; i++ ) 
    fprintf(fp, "%8.2f", graph->hit_scores[i]);
  
  fprintf(fp, "\n");

  fprintf(fp, "\n Upstream Path Scores  \n");
  for(i = 0; i < graph->num_nodes; i++ )
    fprintf(fp, "%8d", i+1);

  fprintf(fp, "\n");

  for(i = 0; i < graph->num_nodes; i++ )
    fprintf(fp, "%8.2f", graph->path_scores[i]);

  fprintf(fp, "\n");

  fprintf(fp, "\n Downstream Path Scores  \n");
  for(i = 0; i < graph->num_nodes; i++ )
    fprintf(fp, "%8d", i+1);

  fprintf(fp, "\n");

   fprintf(fp, "\n Best Edge  \n");
  for(i = 0; i < graph->num_nodes; i++ )
    fprintf(fp, "%8d", i+1);

  fprintf(fp, "\n");

  for(i = 0; i < graph->num_nodes; i++ )
    fprintf(fp, "%8d", graph->best_out_edge[i]+1);

  fprintf(fp, "\n");

 
  if (print_edges) {

    fprintf(fp, "\n Edge Data  \n\n");
    for(i = 0; i < graph->num_nodes; i++ ) {
      for(j = 0; j < graph->num_nodes; j++ ) {
        if(graph->edge_scores[i][j] != -eslINFINITY){
          edge_id =  graph->edge_id[i][j];
          edge = graph->edges[edge_id]; 
          nuc_end =   edge->upstream_spliced_nuc_end;
          nuc_start = edge->downstream_spliced_nuc_start;
          if(graph->revcomp)  {
            nuc_end   = target_seq->n - nuc_end   + target_seq->end;
            nuc_start = target_seq->n - nuc_start + target_seq->end;
           }

           fprintf(fp, "    Edge from Upstream Node %d to Downstream Node %d\n", i+1, j+1);
           fprintf(fp, "                                   %s   %s\n", "Amino", "Nucleotide");
           fprintf(fp, "      Upsteam Node End Coords:     %5d  %10d\n", edge->upstream_spliced_amino_end, nuc_end); 
           fprintf(fp, "      Downsteam Node Start Coords: %5d  %10d\n", edge->downstream_spliced_amino_start, nuc_start);
           fprintf(fp, "\n");
        }
      }
    }
  }
  fprintf(fp, "\n");
 

  return;

}

void
path_dump(FILE *fp, SPLICE_PATH *path, ESL_SQ *target_seq)
{
  
  int i;
  int nuc_start, nuc_end;
  fprintf(fp, "  Path Length  %d\n", path->path_len);
  
  for(i = 0; i < path->path_len; i++) {
    if(path->revcomp) {
       nuc_start = target_seq->n - path->downstream_spliced_nuc_start[i] + target_seq->end;
       nuc_end   = target_seq->n - path->upstream_spliced_nuc_end[i+1] + target_seq->end;
    }
    else {
      nuc_start = path->downstream_spliced_nuc_start[i] + target_seq->start - 1;
      nuc_end   =  path->upstream_spliced_nuc_end[i+1] + target_seq->start - 1;
    }
    fprintf(fp, "  Step %4d Node %4d hmm coords: %5d %5d  seq coords: %10d %10d \n", i+1, path->node_id[i]+1,
      path->downstream_spliced_amino_start[i], path->upstream_spliced_amino_end[i+1], nuc_start, nuc_end);
  }
  
  fprintf(fp, "\n");

  return;
}
