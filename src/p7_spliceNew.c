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

static int TEST_MODE  = 0;



static TARGET_RANGE* target_range_create(int nalloc);
static void target_range_destroy(TARGET_RANGE *target_range);
static SPLICE_EDGE* splice_edge_create(void);
static SPLICE_GRAPH* splice_graph_create(void);
static void splice_graph_destroy (SPLICE_GRAPH* graph);
static int splice_graph_create_nodes(SPLICE_GRAPH *graph, int num_nodes);
static SPLICE_PIPELINE* splice_pipeline_create(const ESL_GETOPTS *go, int M_hint, int L_hint);
static void splice_pipeline_destroy(SPLICE_PIPELINE *pli);

static TARGET_RANGE* build_target_range (TARGET_RANGE *prev_target_range, const P7_TOPHITS *tophits, int64_t *range_bound_mins, int64_t *range_bound_maxs, int range_num,int *hits_processed, int *num_hits_processed);
static ESL_SQ* get_target_range_sequence(TARGET_RANGE *target_range, ESL_SQFILE *seq_file);
static int fill_graph_with_nodes(SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm); 
static int connect_nodes_with_edges(SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *target_seq); 
static void get_overlap_nuc_coords (SPLICE_EDGE *edge, P7_DOMAIN *upstream, P7_DOMAIN *downstream, ESL_SQ *target_seq, int revcomp); 
static int find_optimal_splice_site (SPLICE_EDGE *edge, P7_DOMAIN *upstream, P7_DOMAIN *downstream, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *target_seq, int revcomp);
static float ali_score_at_postion (P7_PROFILE *gm, int amino, int model_pos, int trans_pos, int prev_state, int curr_state);
static int select_splice_option (SPLICE_EDGE *edge, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *target_seq, float *splice_scores, float up_score, float down_score, int model_pos, int up_nuc_pos, int down_nuc_pos, int up_state, int down_state);
static int  evaluate_paths (SPLICE_GRAPH *graph);
static void push_score_upstream (SPLICE_GRAPH *graph, int node_id);
static int splice_path (SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, P7_TOPHITS *th, P7_OPROFILE *om, P7_PROFILE *gm, ESL_SQ *target_seq, ESL_GENCODE *gcode, int start_node, int end_node, int path_len, int *success);
static int align_spliced_path (SPLICE_PIPELINE *pli, P7_OPROFILE *om, ESL_SQ *target_seq, ESL_GENCODE *gcode);
static P7_TRACE* splice_trace_convert(P7_TRACE *orig_tr, int *orig_nuc_idx, int *splice_cnt);
static P7_ALIDISPLAY* create_spliced_alidisplay(const P7_TRACE *tr, int which, const P7_OPROFILE *om, const ESL_SQ *target_seq, const ESL_SQ *amino_sq, int amino_pos, int splice_cnt);

static void target_range_dump(FILE *fp, TARGET_RANGE *target_range, int print_hits);
static void graph_dump(FILE *fp, SPLICE_GRAPH *graph, P7_TOPHITS *th, ESL_SQ *target_seq, int print_edges);

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

  target_range->th->hit    = NULL;
  ESL_ALLOC(target_range->th->hit,  nalloc * sizeof(P7_HIT*));

  target_range->orig_hit_idx = NULL;
  ESL_ALLOC(target_range->orig_hit_idx, nalloc * sizeof(int*));

  target_range->th->unsrt  = NULL;
  target_range->th->N      = 0;
  target_range->th->Nalloc = nalloc;

  return target_range;

  ERROR:
    target_range_destroy(target_range);
    return NULL;
}


/* Free a TARGET_RANGE object */
void
target_range_destroy(TARGET_RANGE *target_range)
{

  if (target_range == NULL) return;

  target_range->seqname = NULL;

  if (target_range->th != NULL)
    p7_tophits_Destroy(target_range->th);

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

  graph->num_nodes  = 0;
  graph->num_edges  = 0;
  graph->num_n_term = 0;
  graph->num_c_term = 0;

  graph->has_full_path         = 0;
  graph->best_full_path_length = 0;
  graph->best_full_path_start  = 0;
  graph->best_full_path_end    = 0;
  graph->best_full_path_score  = 0.0;

  graph->has_out_edge    = NULL; 
  graph->has_in_edge     = NULL;
  graph->best_out_edge   = NULL;
  graph->best_in_edge    = NULL;
  graph->is_n_terminal   = NULL;
  graph->is_c_terminal   = NULL; 

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

  if(graph->edge_id         != NULL) free(graph->edge_id);
  if(graph->edge_id_mem     != NULL) free(graph->edge_id_mem);
  if(graph->edge_scores     != NULL) free(graph->edge_scores);
  if(graph->edge_scores_mem != NULL) free(graph->edge_scores_mem);
  if(graph->path_scores     != NULL) free(graph->path_scores); 
  if(graph->hit_scores      != NULL) free(graph->hit_scores);

  if(graph->has_out_edge    != NULL) free(graph->has_out_edge);
  if(graph->has_in_edge     != NULL) free(graph->has_in_edge);
  if(graph->best_out_edge   != NULL) free(graph->best_out_edge);
  if(graph->best_in_edge    != NULL) free(graph->best_in_edge);
  if(graph->is_n_terminal   != NULL) free(graph->is_n_terminal);
  if(graph->is_c_terminal   != NULL) free(graph->is_c_terminal);
 
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

  /* Allocate adjacency matrix for splice scores */ 
  ESL_ALLOC(graph->edge_scores,       sizeof(float*)       * num_nodes);
  ESL_ALLOC(graph->edge_scores_mem,   sizeof(float)        * (num_nodes * num_nodes));
 
  for (i = 0; i < num_nodes; i++)
    graph->edge_scores[i] = graph->edge_scores_mem + (ptrdiff_t) i * (ptrdiff_t) num_nodes;

  /* Allocate array of hit scores */
  ESL_ALLOC(graph->hit_scores,        sizeof(float)        * num_nodes);

  /* Allocate array of path scores */
  ESL_ALLOC(graph->path_scores,        sizeof(float)        * num_nodes);

  /* Allocate arrays for keeping track of which nodes have edges ans which edge is best*/
  ESL_ALLOC(graph->has_out_edge,   sizeof(int)          * num_nodes);
  ESL_ALLOC(graph->has_in_edge,    sizeof(int)          * num_nodes);
  ESL_ALLOC(graph->best_out_edge,  sizeof(int)          * num_nodes);
  ESL_ALLOC(graph->best_in_edge,   sizeof(int)          * num_nodes);

  /* Allocate arrays of n and c terminal status of each node */
  ESL_ALLOC(graph->is_n_terminal,     sizeof(int)          * num_nodes);
  ESL_ALLOC(graph->is_c_terminal,     sizeof(int)          * num_nodes);

  /* Allocate adjacency matrix for edge index in graph->edges */
  ESL_ALLOC(graph->edge_id,           sizeof(int*)         * num_nodes);
  ESL_ALLOC(graph->edge_id_mem,       sizeof(int)          * (num_nodes * num_nodes));

  for (i = 0; i < num_nodes; i++)
    graph->edge_id[i] = graph->edge_id_mem + (ptrdiff_t) i * (ptrdiff_t) num_nodes;

  /* Allocate array for SPLICE_EDGE objects */
  ESL_ALLOC(graph->edges,             sizeof(SPLICE_EDGE*) * (num_nodes * num_nodes));
  
  
  graph->num_nodes  = num_nodes;

  return eslOK;

  ERROR:
	splice_graph_destroy(graph);
    return status;
}


SPLICE_PIPELINE* 
splice_pipeline_create(const ESL_GETOPTS *go, int M_hint, int L_hint)
{
  SPLICE_PIPELINE *pli;
  int              status;
 
  pli = NULL;
  ESL_ALLOC(pli, sizeof(SPLICE_PIPELINE));


  if (go && esl_opt_GetBoolean(go, "--nonull2")) pli->do_null2 = FALSE;
  else                                           pli->do_null2 = TRUE;           
  
  if (go && esl_opt_IsOn(go, "--incT")) pli->inc_by_E = FALSE;
  else                                  pli->inc_by_E = TRUE;

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
p7_splice_SpliceHits(P7_TOPHITS *tophits, P7_OPROFILE *om, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, ESL_GETOPTS *go, ESL_GENCODE *gcode, ESL_SQFILE *seq_file, FILE *ofp)
{

  int              exon_set_name_id = 0;
  int              i;
  int              hit_cnt;
  int              range_cnt;
  int              success;
  int              spliced_hit_idx;
  int              num_hits_processed;   
  int              is_sorted_by_sortkey;
  int             *hit_processed;
  int             *hit_spliced;
  int64_t         *range_bound_mins;
  int64_t         *range_bound_maxs;
  SPLICE_PIPELINE *pli;
  TARGET_RANGE    *curr_target_range;
  TARGET_RANGE    *prev_target_range;
  SPLICE_GRAPH    *graph;
  ESL_SQ          *target_seq;
  int              status;

  hit_processed    = NULL;
  hit_spliced      = NULL;
  range_bound_mins = NULL;
  range_bound_maxs = NULL;
 

  ESL_STOPWATCH * Timer = esl_stopwatch_Create();
   esl_stopwatch_Start(Timer);
 
  /* Sort hits by sequence, strand and nucleotide positon. 
   * Note if P7_TOPHITS was perviously sorted by sortkey 
   * so we can retur it to it's original state */  
  is_sorted_by_sortkey = tophits->is_sorted_by_sortkey;
  if ((status = p7_tophits_SortBySeqidxAndAlipos(tophits)) != eslOK) goto ERROR;
     
  pli = splice_pipeline_create(go, om->M, om->M * 3);
  pli->bg = p7_bg_Create(om->abc);
  if (pli->do_biasfilter) p7_bg_SetFilter(pli->bg, om->M, om->compo);  

  /* Arrays to keep track of hit status */
  ESL_ALLOC(hit_processed, tophits->N*sizeof(int));
  ESL_ALLOC(hit_spliced,   tophits->N*sizeof(int));

  /* check for any hits that are duplicates or below the 
   * reporting threshold and mark them as processed */
  hit_cnt = 0;
  num_hits_processed = 0;
  for(i = 0; i < tophits->N; i++) {
    if ( !(tophits->hit[i]->flags & p7_IS_REPORTED)) {
      hit_processed[i] = 1;
      num_hits_processed++;
    }
    else { 
      hit_processed[i] = 0;
      hit_cnt++;
    }  
    hit_spliced[i]   = 0;
  }

  /* Arrays for keeping track of allowable target range mins 
   * and maxes set by previous ranges or splicings */
  ESL_ALLOC(range_bound_mins, hit_cnt * sizeof(int64_t));
  ESL_ALLOC(range_bound_maxs, hit_cnt * sizeof(int64_t));

  prev_target_range = NULL;

  /* loop through until all hits have been processed */
  range_cnt = 0;
  while(num_hits_processed < tophits->N) {

    /* build a target range from a set of hits */
    curr_target_range = build_target_range(prev_target_range, tophits, range_bound_mins, range_bound_maxs, range_cnt, hit_processed, &num_hits_processed);

  //printf("Target %s strand %c range %d to %d\n", curr_target_range->seqname, (curr_target_range->revcomp ? '-' : '+'), curr_target_range->start, curr_target_range->end);

    /* Record target range boundries to limit range of other 
     * target ranges on the same sequence and strand */
    range_bound_mins[range_cnt] = curr_target_range->start;
    range_bound_maxs[range_cnt] = curr_target_range->end;
    range_cnt++;

    /* Extend the target range start and end coords to 
     * make room for possible missing terminal exons */  
    curr_target_range->start -= TERM_RANGE_EXT;
    curr_target_range->end   += TERM_RANGE_EXT;  

    /* Fetch a sub-sequence thar corresponds to the target range */
    target_seq = get_target_range_sequence(curr_target_range, seq_file); 
 //printf("Sequence n %d L %d start %d end %d\n", target_seq->n, target_seq->L, target_seq->start, target_seq->end);

    /* Build graph from available hits */
    graph = splice_graph_create();   
    graph->revcomp = curr_target_range->revcomp;
    
    if ((status = fill_graph_with_nodes(graph, curr_target_range->th, gm)) != eslOK) goto ERROR;

    if ((status = connect_nodes_with_edges(graph, curr_target_range->th, gm, gcode, target_seq)) != eslOK) goto ERROR;


    /* Score the paths in the graph.  If there are full mdoel length paths present splice them */
    if ((status = evaluate_paths(graph)) != eslOK) goto ERROR;

     //graph_dump(stdout, graph, curr_target_range->th, target_seq, TRUE); 
    if(graph->has_full_path) {
  
      splice_path(graph, pli, curr_target_range->th, om, gm, target_seq, gcode, graph->best_full_path_start, graph->best_full_path_end, graph->best_full_path_length, &success);
     /* Reset the range_bounds around the spliced hit and set all hits from 
       * the target_range that are outside this new range to unprocessed */
      
      if (success) {
        if (graph->revcomp) {     

          range_bound_mins[range_cnt-1] = pli->hit->dcl->jali;
          range_bound_maxs[range_cnt-1] = pli->hit->dcl->iali;

          for(i = 0; i < curr_target_range->th->N; i++) {
            if      (curr_target_range->th->hit[i]->dcl->jali > range_bound_maxs[range_cnt-1]) {
              spliced_hit_idx =  curr_target_range->orig_hit_idx[i];
              hit_processed[spliced_hit_idx] = 0;
              num_hits_processed--;
            }
            else if (curr_target_range->th->hit[i]->dcl->iali < range_bound_mins[range_cnt-1]) {
              spliced_hit_idx =  curr_target_range->orig_hit_idx[i];
              hit_processed[spliced_hit_idx] = 0;
              num_hits_processed--;
            }
          }
          /* For Rev Comp we might need to reset the splices hit to processed 
           * because the range_bound coords are on the reverse strand*/ 
          spliced_hit_idx =  curr_target_range->orig_hit_idx[graph->best_full_path_start];
          if(hit_processed[spliced_hit_idx] == 0 ) {
            hit_processed[spliced_hit_idx] = 1;
            num_hits_processed++;
          }
        } 
        else {
          range_bound_mins[range_cnt-1] = pli->hit->dcl->iali;
          range_bound_maxs[range_cnt-1] = pli->hit->dcl->jali;
        
          for(i = 0; i < curr_target_range->th->N; i++) {
            if      (curr_target_range->th->hit[i]->dcl->iali > range_bound_maxs[range_cnt-1]) {
              spliced_hit_idx =  curr_target_range->orig_hit_idx[i];
              hit_processed[spliced_hit_idx] = 0;
              num_hits_processed--;
            }
            else if (curr_target_range->th->hit[i]->dcl->jali < range_bound_mins[range_cnt-1]) { 
              spliced_hit_idx =  curr_target_range->orig_hit_idx[i];
              hit_processed[spliced_hit_idx] = 0;
              num_hits_processed--;
            }
          }

        }
        if(TEST_MODE) {
   
          P7_PIPELINE *tmp_pli = p7_pipeline_fs_Create(go, 100, 300, p7_SEARCH_SEQS);
          dump_splash_header(graph, target_seq, exon_set_name_id+1, graph->best_full_path_start, graph->best_full_path_length, pli->hit->dcl->ihmm, pli->hit->dcl->jhmm, pli->hit->dcl->iali, pli->hit->dcl->jali, om, curr_target_range->th, pli->hit->dcl->ad, tmp_pli, ofp, 150);    
          p7_pipeline_fs_Destroy(tmp_pli); 

        }
        exon_set_name_id++;
      
        pli->hit->dcl = NULL;
      } 
      splice_pipeline_reuse(pli);
    }


    
    target_range_destroy(prev_target_range);
    prev_target_range = curr_target_range; 

    splice_graph_destroy(graph);
    esl_sq_Destroy(target_seq);
     
  }

  esl_stopwatch_Stop(Timer);
  if(TEST_MODE) {
    fprintf(stderr,"  Time spent splashing\n  : ");
    esl_stopwatch_Display(stderr,Timer,NULL);
    fprintf(stderr,"\n\n");
  }
  esl_stopwatch_Destroy(Timer);

  /* Leave only footprints */
  if (is_sorted_by_sortkey)
    if ((status = p7_tophits_SortBySortkey(tophits)) != eslOK) goto ERROR;   

  splice_pipeline_destroy(pli);
  target_range_destroy(prev_target_range);

  if(hit_processed    != NULL) free(hit_processed);
  if(hit_spliced      != NULL) free(hit_spliced);
  if(range_bound_mins != NULL) free(range_bound_mins);
  if(range_bound_maxs != NULL) free(range_bound_maxs);
  return status;

  ERROR:
    splice_pipeline_destroy(pli);
    target_range_destroy(prev_target_range);
    splice_graph_destroy(graph);
    esl_sq_Destroy(target_seq);
    if(hit_processed    != NULL) free(hit_processed);
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
build_target_range (TARGET_RANGE *prev_target_range,  const P7_TOPHITS *tophits, int64_t *range_bound_mins, int64_t *range_bound_maxs, int range_num, int *hits_processed, int *num_hits_processed)
{

  int           i;
  int           seqidx;
  int           revcomp;
  int           num_scores;
  int           max_sc_idx;
  int           seed_hit_idx;
  int           curr_hit_idx;
  int           add_hit;
  int           hits_processed_cnt;
  int           seed_hmm_from, seed_hmm_to;
  int64_t       seed_seq_from, seed_seq_to;
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

 
  /* Find first unprocessed hit */
  i = 0;
  while (i < tophits->N && hits_processed[i]) i++;

  /* no hits left to proccess */
  if (i == tophits->N) return NULL;

 
  /* Get sequnce and strand of first unprocessed hit */
  curr_hit = tophits->hit[i];
  seqidx = curr_hit->seqidx;
  if (curr_hit->dcl->iali < curr_hit->dcl->jali)  
    revcomp = 0;
  else
    revcomp = 1;

  /* Check if the previous target range is on the same sequenece and 
   * strand the first unprocessed hit. If so we will keep and use 
   * the range bound values. If not we will reset those values. */ 
  if (prev_target_range != NULL) {
    if(prev_target_range->seqidx != seqidx || prev_target_range->revcomp != revcomp) {
      for (i = 0; i < range_num; i++) {
        range_bound_mins[i] = 0;     
        range_bound_maxs[i] = curr_hit->dcl->ad->L;
      }
      range_num = 0;
    }
  }


  /* Arrays for hit score sorting and tracking*/
  ESL_ALLOC(hit_scores    , (tophits->N - hits_processed_cnt) * sizeof(float));
  ESL_ALLOC(hit_scores_idx, (tophits->N - hits_processed_cnt) * sizeof(int));


  /* Add scores from all hits from the same sequence and 
   * strand as tha first hit to hit_scores */
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

  /* Set target range values and add seed hit to P7_TOPHITS */ 
  target_range->revcomp   = revcomp;
  target_range->seqidx     = seed_hit->seqidx;
  target_range->seqname    = seed_hit->name;
  target_range->th->hit[0] = seed_hit;
  target_range->orig_hit_idx[0] = seed_hit_idx;
  target_range->th->N++;

  seed_seq_from = seed_hit->dcl->iali;
  seed_seq_to   = seed_hit->dcl->jali;

  seed_hmm_from = seed_hit->dcl->ihmm;
  seed_hmm_to   = seed_hit->dcl->jhmm;

  lower_bound = 0;
  upper_bound =  seed_hit->dcl->ad->L;

  curr_range_min = ESL_MIN(seed_seq_from, seed_seq_to);
  curr_range_max = ESL_MAX(seed_seq_from, seed_seq_to);

  seed_ali_min = curr_range_min;
  seed_ali_max = curr_range_max;
  /* Get any upper and lower bounds set by previous ranges or splicings */
  for( i = 0; i < range_num; i++) {
    if(range_bound_maxs[i] < seed_ali_min && range_bound_maxs[i] > lower_bound)
      lower_bound = range_bound_maxs[i];

    if(range_bound_mins[i] > seed_ali_max && range_bound_mins[i] < upper_bound)
      upper_bound = range_bound_mins[i];
  }

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
 
  target_range->start = new_range_min;
  target_range->end   = new_range_max;

  *num_hits_processed = hits_processed_cnt;

  if(hit_scores     != NULL) free(hit_scores);
  if(hit_scores_idx != NULL) free(hit_scores_idx);

  //target_range_dump(stdout, target_range, TRUE);
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
get_target_range_sequence(TARGET_RANGE *target_range, ESL_SQFILE *seq_file)
{
  
  ESL_SQ *target_seq;
  ESL_SQFILE *tmp_file;

  /* Open seq file and ssi */
  esl_sqfile_Open(seq_file->filename,seq_file->format,NULL,&tmp_file);
  esl_sqfile_OpenSSI(tmp_file,NULL);

  /* Get basic sequence info */
  target_seq = esl_sq_Create();
  esl_sqio_FetchInfo(tmp_file,target_range->seqname,target_seq);

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

  int          i;
  int          status;

  status = splice_graph_create_nodes(graph, th->N);

  if (status != eslOK) return status; 

  for(i = 0; i < th->N; i++) {

   graph->hit_scores[i]  = th->hit[i]->score;
   graph->path_scores[i] = -eslINFINITY;

    if(th->hit[i]->dcl->ihmm == 1) {
      graph->is_n_terminal[i] = 1;
      graph->num_n_term++;
    }
    else
      graph->is_n_terminal[i] = 0;

    if(th->hit[i]->dcl->jhmm == gm->M) {
      graph->is_c_terminal[i] = 1;
      graph->num_c_term++;
    }
    else
      graph->is_c_terminal[i] = 0;
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
int
connect_nodes_with_edges(SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *target_seq)
{
    
  int          i,j;
  int          up_amino_start, up_amino_end;
  int          down_amino_start, down_amino_end;
  int          up_nucl_end;
  int          down_nucl_start;
  int          num_ext_aminos;  
  P7_HIT      *upstream_hit;
  P7_HIT      *downstream_hit;
  SPLICE_EDGE *edge;
  int          status;

  for(i = 0; i < th->N; i++) {

    graph->has_out_edge[i]  = 0;
    graph->has_in_edge[i]   = 0;
    graph->best_out_edge[i] = -1;
    graph->best_in_edge[i]  = -1;

    for(j = 0; j < th->N; j++) {

      /* No self egde */
      if (i == j) { 
        graph->edge_scores[i][j] = -eslINFINITY; 
        continue; 
      }

      upstream_hit   = th->hit[i];
      downstream_hit = th->hit[j];
    
      up_amino_start = upstream_hit->dcl->ihmm;
      up_amino_end   = upstream_hit->dcl->jhmm;

      down_amino_start = downstream_hit->dcl->ihmm;
      down_amino_end   = downstream_hit->dcl->jhmm;
      
      /* Is the downstream node actually downstream */
      if (up_amino_start > down_amino_start || up_amino_end > down_amino_end) {
        graph->edge_scores[i][j] = -eslINFINITY; 
        graph->edge_id[i][j]     = -1;
        continue;
      }
  
      /* Is the downstream node close enough to the upstream edge to be the next exon  */
      if (up_amino_end + MAX_AMINO_EXT < down_amino_start) {
        graph->edge_scores[i][j] = -eslINFINITY;
        graph->edge_id[i][j]     = -1;
        continue;
      }
 
      up_nucl_end     = upstream_hit->dcl->jali;
      down_nucl_start = downstream_hit->dcl->iali;

      /* Is there a gap of at least MIN_INTRON_LEN between the nucleotide positions */ 
      if (( graph->revcomp   && (down_nucl_start + MIN_INTRON_LEN > up_nucl_end)) ||
          ((!graph->revcomp) && (down_nucl_start - MIN_INTRON_LEN < up_nucl_end))) {
        graph->edge_scores[i][j] = -eslINFINITY;
        graph->edge_id[i][j]     = -1;
        continue;      
      }

     /* Is there a gap of no more than MAX_INTRON_LEN between the nucleotide positions */ 
     if (( graph->revcomp   && (down_nucl_start + MAX_INTRON_LEN < up_nucl_end)) ||
         ((!graph->revcomp) && (down_nucl_start - MAX_INTRON_LEN > up_nucl_end))) {
        graph->edge_scores[i][j] = -eslINFINITY;
        graph->edge_id[i][j]     = -1;
        continue;
      }      

      /* If the nodes have reached this point we will give them an edge*/
      edge = splice_edge_create();
 
      edge->overlap_amino_start = down_amino_start;
      edge->overlap_amino_end   = up_amino_end;
   
      /* If the hits do not overlap by at least MIN_AMINO_OVERLAP hmm positions, extend them */
      num_ext_aminos = MIN_AMINO_OVERLAP - (edge->overlap_amino_end - edge->overlap_amino_start + 1);
      if (num_ext_aminos > 0) {
        num_ext_aminos = (num_ext_aminos+1)/2;
        edge->overlap_amino_start -= num_ext_aminos;
        edge->overlap_amino_end   += num_ext_aminos;
      }

      get_overlap_nuc_coords(edge, upstream_hit->dcl, downstream_hit->dcl, target_seq, graph->revcomp);
         
      /* Addd extra nucleotides for splice sites */
      if(graph->revcomp) {
        edge->upstream_nuc_end     -= 2;
        edge->downstream_nuc_start += 2;   
      }
      else {
        edge->upstream_nuc_end     += 2;
        edge->downstream_nuc_start -= 2;
      }
      //printf("\nup hmmfrom %d hmmto %d seqfrom %d seqto %d \n", upstream_hit->dcl->ihmm, upstream_hit->dcl->jhmm, upstream_hit->dcl->jali, upstream_hit->dcl->jali);

      //printf("down hmmfrom %d hmmto %d seqfrom %d seqto %d \n", downstream_hit->dcl->ihmm, downstream_hit->dcl->jhmm, downstream_hit->dcl->jali, downstream_hit->dcl->jali);
     
      if ((status = find_optimal_splice_site (edge, upstream_hit->dcl, downstream_hit->dcl, gm, gcode, target_seq, graph->revcomp)) != eslOK) goto ERROR;


      graph->edge_scores[i][j] = edge->splice_score;
      if(edge->splice_score == -eslINFINITY) {
         graph->edge_id[i][j] = -1;
         free(edge);
      }
      else {
        graph->has_out_edge[i] = 1;
        graph->has_in_edge[j]  = 1;
        graph->edges[graph->num_edges] = edge;
        graph->edge_id[i][j] = graph->num_edges;     
        graph->num_edges++;
      }
    }
  }

  return eslOK;

  ERROR:
    if(edge != NULL) free(edge);
    return status;
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
      edge->upstream_nuc_start -= 3 * strand;
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

    if (down_trace->st[z1] != p7T_I) 
      curr_hmm_pos++;

    if (down_trace->st[z1] != p7T_D)
      edge->downstream_nuc_end += 3 * strand;

    if (down_trace->st[z1] > p7T_I) break;
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
find_optimal_splice_site (SPLICE_EDGE *edge, P7_DOMAIN *upstream, P7_DOMAIN *downstream, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *target_seq, int revcomp)
{

  int       z1,z2;
  int       us_idx, ds_idx;
  int       us_start, ds_start;
  int       us_alloc, ds_alloc;
  int       curr_model_pos;
  int       curr_nuc_pos;
  int       amino; 
  int       prev_state;
  float     baseline_score;
  int      *up_model_pos;
  int      *down_model_pos;
  int      *up_nuc_pos;
  int      *down_nuc_pos;
  int      *up_states;
  int      *down_states;
  float    *up_scores;
  float    *down_scores;
  float    *splice_scores;
  P7_TRACE *up_trace;
  P7_TRACE *down_trace;
  int       status;

  /***************************** UPSTREAM **********************************/

  up_model_pos = NULL;
  up_nuc_pos  = NULL;
  up_states   = NULL;
  up_scores   = NULL;

  up_trace = upstream->tr;

  /* Allocate enough room for all trace positions in the original overlap plus the overlap extention */
  us_alloc = (edge->upstream_trace_end - edge->upstream_trace_start + 1) + (edge->overlap_amino_end - upstream->jhmm);

  ESL_ALLOC(up_model_pos, us_alloc*sizeof(int));
  ESL_ALLOC(up_nuc_pos,  us_alloc*sizeof(int));
  ESL_ALLOC(up_states,   us_alloc*sizeof(int));
  ESL_ALLOC(up_scores,   us_alloc*sizeof(float));

  /* Zeroth position is place holder for transtioning into first postion */
   up_nuc_pos[0]  = -1;
   up_model_pos[0] = -1;
   up_states[0]   = p7T_X;
   up_scores[0]   = 0.;
   
  /* Calculate a Viterbi-like score for the upstream hit in the overlap region */
  curr_model_pos = edge->overlap_amino_start;
  curr_nuc_pos = edge->upstream_nuc_start;
  prev_state     = p7T_X;
  us_idx         = 1;

  /* Overlap region that is part of original hit */
  z2 = edge->upstream_trace_start;
  while (z2 < edge->upstream_trace_end) {

    up_nuc_pos[us_idx]  = curr_nuc_pos;
    up_model_pos[us_idx] = curr_model_pos;
    up_states[us_idx]   = up_trace->st[z2];
     
    if(up_states[us_idx] == p7T_M)
      amino = esl_gencode_GetTranslation(gcode, &target_seq->dsq[curr_nuc_pos]);
     
    up_scores[us_idx] = up_scores[us_idx-1];
    up_scores[us_idx] += ali_score_at_postion(gm, amino, curr_model_pos, curr_model_pos-1, prev_state, up_states[us_idx]);

   if(up_states[us_idx] != p7T_D) 
      curr_nuc_pos += 3;
    if(up_states[us_idx] != p7T_I)
      curr_model_pos++;

    prev_state = up_states[us_idx];
    us_idx++;
    z2++;
  }

  baseline_score = up_scores[us_idx-1];

  /*Overlap region that is part of extention */
  while (curr_model_pos <= edge->overlap_amino_end) {

    up_nuc_pos[us_idx]  = curr_nuc_pos;
    up_model_pos[us_idx] = curr_model_pos;
    up_states[us_idx]   = p7T_M;

    amino = esl_gencode_GetTranslation(gcode, &target_seq->dsq[curr_nuc_pos]);
    
    up_scores[us_idx]  = up_scores[us_idx-1];
    up_scores[us_idx] += ali_score_at_postion(gm, amino, curr_model_pos, curr_model_pos-1,prev_state, p7T_M); 

    curr_nuc_pos  += 3;
    curr_model_pos++;

    prev_state = p7T_M;
    us_idx++;

  }
  
  /***************************** DONWSTREAM **********************************/
 
  down_model_pos = NULL;
  down_nuc_pos  = NULL;
  down_states   = NULL;
  down_scores   = NULL;

  down_trace = downstream->tr;

   /* Allocate enough room for all trace positions in the original overlap plus the overlap extention */ 
  ds_alloc = (edge->downstream_trace_end - edge->downstream_trace_start + 1) + (downstream->ihmm - edge->overlap_amino_start + 1);

  ESL_ALLOC(down_model_pos, ds_alloc*sizeof(int));
  ESL_ALLOC(down_nuc_pos,  ds_alloc*sizeof(int));
  ESL_ALLOC(down_states,   ds_alloc*sizeof(int));
  ESL_ALLOC(down_scores,   ds_alloc*sizeof(float));

  /* Last position is place holder for transtioning into Last-1 postion */
  down_nuc_pos[ds_alloc-1]  = -1;
  down_model_pos[ds_alloc-1] = -1;
  down_states[ds_alloc-1]   = p7T_X;
  down_scores[ds_alloc-1]   = 0.; 

  /* Calculate a Viterbi-like score for the donwstream hit in the overlap region */
  curr_model_pos = edge->overlap_amino_end;
  curr_nuc_pos   = edge->downstream_nuc_end-2;
  prev_state     = p7T_X;
  ds_idx         = ds_alloc-2;
 
  z1 = edge->downstream_trace_end;
  while (z1 > edge->downstream_trace_start) {
   
    down_nuc_pos[ds_idx]   = curr_nuc_pos-2;
    down_model_pos[ds_idx] = curr_model_pos;
    down_states[ds_idx]    = down_trace->st[z1];

    if(down_states[ds_idx] == p7T_M)
      amino = esl_gencode_GetTranslation(gcode, &target_seq->dsq[curr_nuc_pos]);

    down_scores[ds_idx] = down_scores[ds_idx+1];
    down_scores[ds_idx] += ali_score_at_postion(gm, amino, curr_model_pos, curr_model_pos, down_states[ds_idx], prev_state);

    if(down_states[ds_idx] != p7T_D)
      curr_nuc_pos -= 3;
    if(down_states[ds_idx] != p7T_I)
      curr_model_pos--;

    prev_state = down_states[ds_idx];
    ds_idx--;
    z1--; 

  }

  baseline_score += down_scores[ds_idx+1];

  /*Overlap region that is part of extention */
  while (curr_model_pos >= edge->overlap_amino_start) {
    
    down_nuc_pos[ds_idx]  = curr_nuc_pos - 2;
    down_model_pos[ds_idx] = curr_model_pos;
    down_states[ds_idx]   = p7T_M;

    amino = esl_gencode_GetTranslation(gcode, &target_seq->dsq[curr_nuc_pos]);

    down_scores[ds_idx]  = down_scores[ds_idx+1];
    down_scores[ds_idx] += ali_score_at_postion(gm, amino, curr_model_pos, curr_model_pos, p7T_M, prev_state);

    curr_nuc_pos  -= 3;
    curr_model_pos--;

    prev_state = p7T_M;
    ds_idx--; 
  }

  /* Initialize splice signal score array (in bitscore) */
  splice_scores = NULL;
  ESL_ALLOC(splice_scores, sizeof(float) * p7S_SPLICE_SIGNALS);
  p7_splice_SignalScores(splice_scores);

  /* Loop through possible splice sites and fine highest scoring splicing */
  ds_start = ds_idx+1;
  us_start = 1;  
  for (curr_model_pos = edge->overlap_amino_start; curr_model_pos <=  edge->overlap_amino_end; curr_model_pos++)
  {
    while (up_model_pos[us_start] < curr_model_pos) us_start++;
    while (down_model_pos[ds_start] < curr_model_pos) ds_start++;
   
    us_idx = us_start;
    while(us_idx < us_alloc && up_model_pos[us_idx] == curr_model_pos) { 
      ds_idx = ds_start;
      while(ds_idx < ds_alloc && down_model_pos[ds_idx] == curr_model_pos) {
        if ((status = select_splice_option(edge, gm, gcode, target_seq, splice_scores, up_scores[us_idx-1], down_scores[ds_idx+1], curr_model_pos, up_nuc_pos[us_idx], down_nuc_pos[ds_idx], up_states[us_idx-1], down_states[ds_idx+1])) != eslOK) goto ERROR; 

        ds_idx++;
      }
      us_idx++;
    }

    if (us_idx >= us_alloc || ds_idx >= ds_alloc) break;
    us_start = us_idx;
    ds_start = ds_idx;
  }
 
  edge->splice_score -= baseline_score;

  if(up_model_pos != NULL) free(up_model_pos);
  if(up_nuc_pos   != NULL) free(up_nuc_pos);
  if(up_states    != NULL) free(up_states);
  if(up_scores    != NULL) free(up_scores);

  if(down_model_pos != NULL) free(down_model_pos);
  if(down_nuc_pos   != NULL) free(down_nuc_pos);
  if(down_states    != NULL) free(down_states);
  if(down_scores    != NULL) free(down_scores);

  if(splice_scores != NULL) free(splice_scores);

  return eslOK;

  ERROR:
    if(up_model_pos != NULL) free(up_model_pos);
    if(up_nuc_pos   != NULL) free(up_nuc_pos);
    if(up_states    != NULL) free(up_states);
    if(up_scores    != NULL) free(up_scores);

    if(down_model_pos != NULL) free(down_model_pos);
    if(down_nuc_pos   != NULL) free(down_nuc_pos);
    if(down_states    != NULL) free(down_states);
    if(down_scores    != NULL) free(down_scores);

    if(splice_scores != NULL) free(splice_scores);
 
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
select_splice_option (SPLICE_EDGE *edge, P7_PROFILE *gm,  ESL_GENCODE *gcode, ESL_SQ *target_seq, float *splice_scores, float up_score, float down_score, int model_pos, int up_nuc_pos, int down_nuc_pos, int up_state, int down_state)
{

  int      amino;
  int      best_opt;
  int      donor_one;
  int      donor_two;
  int      accept_one;
  int      accept_two;
  float    opt_score;
  float    trans_score;
  float    splice_score;
  ESL_DSQ *codon;
  int      status;

  codon = NULL;
  ESL_ALLOC(codon, sizeof(ESL_DSQ) * 3);

  best_opt      = -1;
  splice_score = -eslINFINITY;

  /* option 0 |ABCxx|...yy| */
  codon[0] = target_seq->dsq[up_nuc_pos];
  codon[1] = target_seq->dsq[up_nuc_pos+1];
  codon[2] = target_seq->dsq[up_nuc_pos+2];

  amino     = esl_gencode_GetTranslation(gcode,&codon[0]);
  if (amino < gcode->aa_abc->K) {

    opt_score = p7P_MSC(gm,model_pos,amino);

    donor_one  = target_seq->dsq[up_nuc_pos+3];
    donor_two  = target_seq->dsq[up_nuc_pos+4];
    accept_one = target_seq->dsq[down_nuc_pos+3];
    accept_two = target_seq->dsq[down_nuc_pos+4];

    /* GT-AG */
    if      (donor_one == 2 && donor_two == 3 && accept_one == 0 && accept_two == 2) 
      opt_score += splice_scores[p7S_GTAG];
    /* GC-AG */
    else if (donor_one == 2 && donor_two == 1 && accept_one == 0 && accept_two == 2) 
      opt_score += splice_scores[p7S_GCAG];
    /* AT-AC */
    else if (donor_one == 0 && donor_two == 1 && accept_one == 0 && accept_two == 1)
      opt_score += splice_scores[p7S_ATAC];
    /* OTHER */
    else
      opt_score += splice_scores[p7S_OTHER];  

    best_opt     = 0;
    splice_score = opt_score;
  }

  /* option 1 |ABxx.|..yyF| */
  codon[2] = target_seq->dsq[down_nuc_pos+4];

  amino = esl_gencode_GetTranslation(gcode,&codon[0]);
 
  if (amino < gcode->aa_abc->K) {

    opt_score = p7P_MSC(gm,model_pos,amino);

    donor_one  = target_seq->dsq[up_nuc_pos+2];
    donor_two  = target_seq->dsq[up_nuc_pos+3];
    accept_one = target_seq->dsq[down_nuc_pos+2];
    accept_two = target_seq->dsq[down_nuc_pos+3];

   /* GT-AG */
    if      (donor_one == 2 && donor_two == 3 && accept_one == 0 && accept_two == 2)
      opt_score += splice_scores[p7S_GTAG];
    /* GC-AG */
    else if (donor_one == 2 && donor_two == 1 && accept_one == 0 && accept_two == 2)
      opt_score += splice_scores[p7S_GCAG];
    /* AT-AC */
    else if (donor_one == 0 && donor_two == 1 && accept_one == 0 && accept_two == 1)
      opt_score += splice_scores[p7S_ATAC];
    /* OTHER */
    else
      opt_score += splice_scores[p7S_OTHER];
  
    if(opt_score > splice_score) {
      best_opt = 1;
      splice_score = opt_score;
    }   
  }

  /* option 2 |Axx..|.yyEF|  */
  codon[1] = target_seq->dsq[down_nuc_pos+3];

  amino = esl_gencode_GetTranslation(gcode,&codon[0]);

  if (amino < gcode->aa_abc->K) {

    opt_score = p7P_MSC(gm,model_pos,amino);
  
    donor_one  = target_seq->dsq[up_nuc_pos+1];
    donor_two  = target_seq->dsq[up_nuc_pos+2];
    accept_one = target_seq->dsq[down_nuc_pos+1];
    accept_two = target_seq->dsq[down_nuc_pos+2];

    /* GT-AG */
    if      (donor_one == 2 && donor_two == 3 && accept_one == 0 && accept_two == 2)
      opt_score += splice_scores[p7S_GTAG];
    /* GC-AG */
    else if (donor_one == 2 && donor_two == 1 && accept_one == 0 && accept_two == 2)
      opt_score += splice_scores[p7S_GCAG];
    /* AT-AC */
    else if (donor_one == 0 && donor_two == 1 && accept_one == 0 && accept_two == 1)
      opt_score += splice_scores[p7S_ATAC];
    /* OTHER */
    else
      opt_score += splice_scores[p7S_OTHER];

    if(opt_score > splice_score) {
      best_opt = 2;
      splice_score = opt_score;
    }
  }

  /* option 3 |xx...|yyDEF| */
  codon[0] = target_seq->dsq[down_nuc_pos+2];

  amino = esl_gencode_GetTranslation(gcode,&codon[0]);

  if (amino < gcode->aa_abc->K) { 

    opt_score = p7P_MSC(gm,model_pos,amino);

    donor_one  = target_seq->dsq[up_nuc_pos];
    donor_two  = target_seq->dsq[up_nuc_pos+1];
    accept_one = target_seq->dsq[down_nuc_pos];
    accept_two = target_seq->dsq[down_nuc_pos+1];

    /* GT-AG */
    if      (donor_one == 2 && donor_two == 3 && accept_one == 0 && accept_two == 2)
      opt_score += splice_scores[p7S_GTAG];
    /* GC-AG */
    else if (donor_one == 2 && donor_two == 1 && accept_one == 0 && accept_two == 2)
      opt_score += splice_scores[p7S_GCAG];
    /* AT-AC */
    else if (donor_one == 0 && donor_two == 1 && accept_one == 0 && accept_two == 1)
      opt_score += splice_scores[p7S_ATAC];
    /* OTHER */
    else
      opt_score += splice_scores[p7S_OTHER];

    if(opt_score > splice_score) {
      best_opt = 3;
      splice_score = opt_score;
    }
  }

  if      (up_state   == p7T_M) trans_score = p7P_TSC(gm,model_pos-1,p7P_MM);
  else if (up_state   == p7T_D) trans_score = p7P_TSC(gm,model_pos-1,p7P_DM);
  else if (up_state   == p7T_I) trans_score = p7P_TSC(gm,model_pos-1,p7P_IM);
  else                          trans_score = 0.;

  if      (down_state == p7T_M) trans_score += p7P_TSC(gm,model_pos,p7P_MM);
  else if (down_state == p7T_D) trans_score += p7P_TSC(gm,model_pos,p7P_MD);
  else if (down_state == p7T_I) trans_score += p7P_TSC(gm,model_pos,p7P_MI);
  else                          trans_score += 0.;

  splice_score += (trans_score + up_score + down_score); 

  if (splice_score > edge->splice_score) {
  
    edge->splice_score = splice_score;
    
    edge->upstream_spliced_nuc_end     = up_nuc_pos   + 2 - best_opt;
    edge->downstream_spliced_nuc_start = down_nuc_pos + 5 - best_opt;
    
    if( best_opt < 2) {
      edge-> upstream_spliced_amino_end    = model_pos;
      edge->downstream_spliced_amino_start = model_pos + 1;
    }
    else {
      edge-> upstream_spliced_amino_end    = model_pos - 1;
      edge->downstream_spliced_amino_start = model_pos;
    }
  }

  if(codon != NULL) free(codon);

  return eslOK;

  ERROR:
    if(codon != NULL) free(codon);
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
  int   status;

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


  if(model_pos > gm->M) ESL_XEXCEPTION(eslFAIL, "impossible model position reached %d for model of length %d", model_pos, gm->M);
  if (transition != p7P_NTRANS)
    transition_score = p7P_TSC(gm,trans_pos,transition);
  else
    transition_score = 0.; //special case for first or last potition in overlap
//  printf("transition_score %f emission_score %f\n", transition_score, emission_score);
  return transition_score + emission_score;

  ERROR:
    return -eslINFINITY;

}



int
evaluate_paths (SPLICE_GRAPH *graph)
{

  int          i;
  int          path_len;
  int          start_node;
  int          prev_node;
  int          curr_node;
  int          next_node;
  int          edge_id;
  int          path_is_broken;
  int          path_num;
  int          max_sc_idx;
  int         *full_path_starts;
  float       *full_path_scores;
  SPLICE_EDGE *in_edge;
  SPLICE_EDGE *out_edge;
  int          status;


  /* Find best scoreing paths */ 
  for (i = 0; i < graph->num_nodes; i++) {
    if( ! graph->has_out_edge[i] ) {
       graph->path_scores[i] = graph->hit_scores[i];
       push_score_upstream(graph, i);
    }
  }

  
  /* Find best scoring full length paths */
  full_path_starts  = NULL;
  full_path_scores = NULL;

  /* No full paths */
  if(graph->num_n_term == 0 || graph->num_c_term == 0)  
     return eslOK;

  path_is_broken = FALSE;
  path_len       = 0;
  /* Only one possible full path. 
   *   1. Count the length of the path 
   *   2. Check for beckward nodes
   *   3. Check if path is full lenghth */
  if(graph->num_n_term == 1) {
     start_node = 0;
     while(!graph->is_n_terminal[start_node]) start_node++;
     curr_node = start_node;
     path_len  = 1;
     if(graph->has_out_edge[start_node]) {
       curr_node = graph->best_out_edge[start_node];
       edge_id   = graph->edge_id[start_node][curr_node];
       in_edge   = graph->edges[edge_id];
       prev_node = start_node;
       path_len++;
     }
     while(graph->has_out_edge[curr_node]) {
       next_node = graph->best_out_edge[curr_node];
       edge_id   = graph->edge_id[curr_node][next_node];
       out_edge  = graph->edges[edge_id];
       /* If this is a backward node (in edge ends after out edge begins) delete the out edge */
      if(in_edge->downstream_spliced_amino_start > out_edge->upstream_spliced_amino_end){ 
        
        if(out_edge->splice_score < in_edge->splice_score) {

          graph->edge_scores[curr_node][next_node] = -eslINFINITY;
          graph->has_out_edge[curr_node] = 0;
          for(i = 0; i < graph->num_nodes; i++) {
            if(graph->edge_scores[curr_node][i] != -eslINFINITY)
              graph->has_out_edge[curr_node] = 1;
            }
        }
        else {

          graph->edge_scores[prev_node][curr_node] = -eslINFINITY;
          graph->has_out_edge[prev_node] = 0;
          for(i = 0; i < graph->num_nodes; i++) {
            if(graph->edge_scores[prev_node][i] != -eslINFINITY)
              graph->has_out_edge[prev_node] = 1;
          }
        }
        path_is_broken = TRUE;
        break;   
      }
      in_edge   = out_edge;
      prev_node = curr_node;
      curr_node = next_node;
      path_len++;
    } 
   
    /* If we found a backward node we need to start over */ 
    if(path_is_broken) { 
      for(i = 0; i < graph->num_nodes; i++) {
        graph->path_scores[i]   = -eslINFINITY;
        graph->best_out_edge[i] = -1;
      }
      evaluate_paths (graph);
      return eslOK;
    }
    /* Check if path is full length */
    else if(path_len > 1 && graph->is_c_terminal[curr_node]) {
        graph->has_full_path         = 1;
        graph->best_full_path_start  = start_node;
        graph->best_full_path_end    = curr_node;
        graph->best_full_path_length = path_len;
        graph->best_full_path_score  = graph->path_scores[start_node];   
    }
    return eslOK;
  }

  /* Multiple possible full paths - find with highest scoring */
  if (graph->num_n_term > 1) {
    ESL_ALLOC(full_path_starts, sizeof(int)   * graph->num_n_term);
    ESL_ALLOC(full_path_scores, sizeof(float) * graph->num_n_term);
 
    path_num = 0;
    for (i = 0; i < graph->num_nodes; i++) {
      if(graph->is_n_terminal[i]) {
         full_path_scores[path_num] = graph->path_scores[i];
         full_path_starts[path_num]  = i;
         path_num++;
      }
    }

    while((!graph->has_full_path) &&  path_num > 0) { 
      max_sc_idx = esl_vec_FArgMax(full_path_scores, graph->num_n_term); 
      path_len   = 1;
      start_node = full_path_starts[max_sc_idx];
      curr_node = start_node; 
      if(graph->has_out_edge[start_node]) {
        curr_node = graph->best_out_edge[start_node];
        edge_id   = graph->edge_id[start_node][curr_node];
        in_edge   = graph->edges[edge_id];
        prev_node = start_node;
        path_len++;
      }
      while(graph->has_out_edge[curr_node]) {
        next_node = graph->best_out_edge[curr_node];
        edge_id   = graph->edge_id[curr_node][next_node];
        out_edge  = graph->edges[edge_id];

        /* If this is a backward node (in edge ends after out edge begins) 
         * delete the edge with the lower splice score*/
        if(in_edge->downstream_spliced_amino_start > out_edge->upstream_spliced_amino_end){
          if(out_edge->splice_score < in_edge->splice_score) {

            graph->edge_scores[curr_node][next_node] = -eslINFINITY;
            graph->has_out_edge[curr_node] = 0;
            for(i = 0; i < graph->num_nodes; i++) {
              if(graph->edge_scores[curr_node][i] != -eslINFINITY)
                graph->has_out_edge[curr_node] = 1;
            }
          } 
          else {

            graph->edge_scores[prev_node][curr_node] = -eslINFINITY;
            graph->has_out_edge[prev_node] = 0;
            for(i = 0; i < graph->num_nodes; i++) {
              if(graph->edge_scores[prev_node][i] != -eslINFINITY)
                graph->has_out_edge[prev_node] = 1;         
            }
          }
          path_is_broken = TRUE;
          break;
        } 
        in_edge   = out_edge;
        prev_node = curr_node;
        curr_node = next_node;
        path_len++;
      }
     
      /* If we found a backward node we need to start over */ 
      if(path_is_broken) {
        if(full_path_starts != NULL) free(full_path_starts);
        if(full_path_scores != NULL) free(full_path_scores);

        for(i = 0; i < graph->num_nodes; i++) {
          graph->path_scores[i]   = -eslINFINITY;
          graph->best_out_edge[i] = -1;
        }
        evaluate_paths (graph);
        return eslOK;
      }
      /* Check if path is full length */
      else if(path_len > 1 && graph->is_c_terminal[curr_node]) {
        graph->has_full_path         = 1;
        graph->best_full_path_start  = start_node;
        graph->best_full_path_end    = curr_node;
        graph->best_full_path_length = path_len;
        graph->best_full_path_score  = graph->path_scores[start_node]; 
      }
      else
        path_num--;
    }  
  }     
 
    
  if(full_path_starts != NULL) free(full_path_starts);
  if(full_path_scores != NULL) free(full_path_scores);
  return eslOK;
 
  ERROR:
   if(full_path_starts != NULL) free(full_path_starts);
   if(full_path_scores != NULL) free(full_path_scores);
   return status;
 
}



void
push_score_upstream (SPLICE_GRAPH *graph, int node_id)
{

  int up_node, down_node;
  float step_score;

  for (up_node = 0 ; up_node < graph->num_nodes; up_node++) {
    /* If there is a edge between these nodes scheck if pushing the 
     * score from the downstream node will increase teh score of the 
      * best best path to the upstream node */ 
    if (graph->edge_scores[up_node][node_id] != -eslINFINITY) {
      step_score  = graph->edge_scores[up_node][node_id];
      step_score += graph->path_scores[node_id];
      step_score += graph->hit_scores[up_node]; 
            
      if (step_score > graph->path_scores[up_node]) {
        graph->path_scores[up_node] = step_score;
        graph->best_out_edge[up_node] = node_id;
      }
      push_score_upstream(graph, up_node);
    }
  }
  return;
}






int
splice_path (SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, P7_TOPHITS *th, P7_OPROFILE *om, P7_PROFILE *gm, ESL_SQ *target_seq, ESL_GENCODE *gcode, int start_node, int end_node, int path_len, int* success) 
{

  int          i;
  int          pos;
  int          seq_idx;
  int          exon_cnt;
  int          curr_node;
  int          next_node;
  int          curr_edge;
  int          seq_len;
  int          amino_len;
  int          amino;
  int         *start_pos;
  int         *end_pos;
  int         *nuc_index;
  float        sum_hit_sc;
  float        dom_bias;
  float        nullsc;
  float        dom_score;
  float        dom_lnP;
  ESL_DSQ     *nuc_dsq;
  ESL_DSQ     *amino_dsq;
  ESL_SQ      *nuc_seq;
  ESL_SQ      *amino_seq;  
  P7_HIT      *node_hit;
  SPLICE_EDGE *edge;
  int          status;

  char          *alphaAmino;
  char          *alphaDNA;
  alphaAmino = om->abc->sym;
  alphaDNA   = target_seq->abc->sym;
  
  start_pos = NULL;
  end_pos   = NULL;
  nuc_index = NULL;
  nuc_dsq   = NULL;

  ESL_ALLOC(start_pos, sizeof(int64_t) * path_len);
  ESL_ALLOC(end_pos,   sizeof(int64_t) * path_len);

  /* Get start and end positions of first exon in path */
  node_hit     = th->hit[start_node];
  start_pos[0] = node_hit->dcl->iali;
  if(graph->revcomp) start_pos[0] = target_seq->start - start_pos[0] + 1;
  else               start_pos[0] = start_pos[0] - target_seq->start + 1; 

  /* Get summed hit score for all hits in path */
  sum_hit_sc = graph->hit_scores[start_node];

  next_node  = graph->best_out_edge[start_node];
  curr_edge  = graph->edge_id[start_node][next_node];

  sum_hit_sc += graph->hit_scores[next_node]; 

  edge = graph->edges[curr_edge];

  end_pos[0] = edge->upstream_spliced_nuc_end; 
  
  curr_node = next_node;
  exon_cnt  = 1;
   
  /* Get start and end positions of middle exons in path */
  while (exon_cnt < (path_len-1)) {

    start_pos[exon_cnt] = edge->downstream_spliced_nuc_start;
     
    next_node  = graph->best_out_edge[curr_node];
    curr_edge = graph->edge_id[curr_node][next_node];

    sum_hit_sc += graph->hit_scores[next_node];

    edge = graph->edges[curr_edge];

    end_pos[exon_cnt] = edge->upstream_spliced_nuc_end;
  
    curr_node = next_node;
    exon_cnt++;
  }

  /* Get start and end positions of last exon in path */
  start_pos[exon_cnt] = edge->downstream_spliced_nuc_start;

  node_hit = th->hit[next_node];
  end_pos[exon_cnt] = node_hit->dcl->jali;
  if(graph->revcomp) end_pos[exon_cnt] = target_seq->start - end_pos[exon_cnt] + 1;
  else               end_pos[exon_cnt] = end_pos[exon_cnt] - target_seq->start + 1;
 
  /* Get length of spliced nucleotide sequence */
  seq_len = 0;
  for (i = 0; i < path_len; i++)  
    seq_len += (end_pos[i] - start_pos[i] + 1);
  
  ESL_ALLOC(nuc_index, sizeof(int64_t) * (seq_len+2));  
  ESL_ALLOC(nuc_dsq,   sizeof(ESL_DSQ) * (seq_len+2));

  /* Copy spliced nucleotides into single sequence and track their original indicies */
  nuc_index[0] = -1;
  nuc_dsq[0]   = eslDSQ_SENTINEL;
  seq_idx   = 1;


  for (i = 0; i < path_len; i++) {
    for (pos = start_pos[i]; pos <= end_pos[i]; pos++) {
      nuc_index[seq_idx] = pos;
      nuc_dsq[seq_idx]   = target_seq->dsq[pos];
      seq_idx++;    
    }
  }
 
  nuc_index[seq_idx] = -1;
  nuc_dsq[seq_idx]   = eslDSQ_SENTINEL;


  /* Translate spliced nucleotide sequence */
  amino_len = seq_len/3;
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
  amino_seq   = esl_sq_CreateDigitalFrom(gcode->nt_abc, target_seq->name, amino_dsq, amino_len,NULL,NULL,NULL); 
  nuc_seq     = esl_sq_CreateDigitalFrom(gcode->aa_abc, target_seq->name, nuc_dsq,   seq_len,  NULL,NULL,NULL);

  pli->nuc_sq       = nuc_seq;
  pli->amino_sq     = amino_seq;   
  pli->orig_nuc_idx = nuc_index;

  /* Algin the splices Amino sequence */
  if ((status = align_spliced_path(pli, om, target_seq, gcode)) != eslOK) goto ERROR;
  
  //Alignment failed
  if(pli->hit == NULL) {
    *success = FALSE; 
     return eslOK; 
  }

  if(graph->revcomp) {
    pli->hit->dcl->iali       += 2;
    pli->hit->dcl->ad->sqfrom += 2; 
    pli->hit->dcl->ienv = target_seq->n - pli->orig_nuc_idx[1]              + target_seq->end; 
    pli->hit->dcl->jenv = target_seq->n - pli->orig_nuc_idx[pli->nuc_sq->n] + target_seq->end; 
  }  
  else {
    pli->hit->dcl->ienv = pli->orig_nuc_idx[1]              + target_seq->start -1;
    pli->hit->dcl->jenv = pli->orig_nuc_idx[pli->nuc_sq->n] + target_seq->start -1;
  }
 
  /* Bias calculation and adjustments to Forward score */
  if(pli->do_null2)
    dom_bias = p7_FLogsum(0.0, log(pli->bg->omega) + pli->hit->dcl->domcorrection);
  else
    dom_bias = 0.;

  p7_bg_SetLength(pli->bg, pli->amino_sq->n);
  p7_bg_NullOne  (pli->bg, pli->amino_sq->dsq, pli->amino_sq->n, &nullsc);
  dom_score = (pli->hit->dcl->envsc - (nullsc + dom_bias))  / eslCONST_LOG2;
  
  /* Subrtact extra B->M trasition scores from summed hit score */
  sum_hit_sc -= gm->tsc[p7P_NTRANS+p7P_BM] * (path_len-1); 
  
  /* If the spliced hit score is higher than the summed hit score, replace 
   * the first domain from th first hit in the path with the domain from 
   * the spliced hit and set all other hits as being bellow the inculsion 
   * threshhold so they are not redported */ 
  
  if(dom_score > sum_hit_sc) {

    dom_lnP   = esl_exp_logsurv(dom_score, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

    pli->hit->dcl->bitscore    = dom_score;
    pli->hit->dcl->lnP         = dom_lnP;
    pli->hit->dcl->dombias     = dom_bias;
    pli->hit->dcl->is_reported = TRUE;
    pli->hit->dcl->is_included = TRUE; 
    node_hit = th->hit[start_node];
     
    p7_domain_Destroy(node_hit->dcl);
    node_hit->dcl = pli->hit->dcl;
    
    node_hit->pre_score = node_hit->dcl->envsc  / eslCONST_LOG2;
    node_hit->pre_lnP   = esl_exp_logsurv (node_hit->pre_score, om->evparam[p7_FTAUFS], om->evparam[p7_FLAMBDA]);
    
    node_hit->sum_score  = node_hit->score  = dom_score; 
    node_hit->sum_lnP    = node_hit->lnP    = dom_lnP;

    node_hit->sortkey    = pli->inc_by_E ? -dom_lnP : dom_score;

    node_hit->frameshift = FALSE;

    next_node = graph->best_out_edge[start_node];
    node_hit  = th->hit[next_node];
   
    node_hit->flags &= ~p7_IS_REPORTED;
    node_hit->flags &= ~p7_IS_INCLUDED;
    node_hit->dcl->is_reported = FALSE;
    node_hit->dcl->is_included = FALSE; 

    curr_node = next_node;
    exon_cnt  = 1;
    while (exon_cnt < (path_len-1)) {
      next_node  = graph->best_out_edge[curr_node];
      node_hit  = th->hit[next_node];

      node_hit->flags &= ~p7_IS_REPORTED;
      node_hit->flags &= ~p7_IS_INCLUDED;
      node_hit->dcl->is_reported = FALSE;
      node_hit->dcl->is_included = FALSE;

      curr_node = next_node;
      exon_cnt++;
    }

    *success = TRUE;
  }
  else *success = FALSE;
  
  if(start_pos != NULL) free(start_pos);
  if(end_pos   != NULL) free(end_pos);
  if(nuc_dsq   != NULL) free(nuc_dsq);
  if(amino_dsq != NULL) free(amino_dsq);  

  return eslOK;
 
  ERROR:
    if(start_pos != NULL) free(start_pos);
    if(end_pos   != NULL) free(end_pos);
    if(nuc_index != NULL) free(nuc_index);
    if(nuc_dsq   != NULL) free(nuc_dsq);
    if(amino_dsq != NULL) free(amino_dsq);
    return status;

}


int 
align_spliced_path (SPLICE_PIPELINE *pli, P7_OPROFILE *om, ESL_SQ *target_seq, ESL_GENCODE *gcode) 
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

  p7_oprofile_ReconfigUnihit(om, pli->amino_sq->L);
  p7_omx_GrowTo(pli->fwd, om->M, pli->amino_sq->L, pli->amino_sq->L);
  p7_omx_GrowTo(pli->bwd, om->M, pli->amino_sq->L, pli->amino_sq->L);

  p7_bg_SetLength(pli->bg, pli->amino_sq->n);
  if (pli->do_biasfilter)
    p7_bg_FilterScore(pli->bg, pli->amino_sq->dsq, pli->amino_sq->n, &filtersc);
  else
    p7_bg_NullOne  (pli->bg, pli->amino_sq->dsq, pli->amino_sq->n, &filtersc);

  p7_Forward (pli->amino_sq->dsq, pli->amino_sq->L, om,           pli->fwd, &envsc);

  seq_score = (envsc-filtersc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score,  om->evparam[p7_FTAU],  om->evparam[p7_FLAMBDA]);
  
  if (P > pli->F3) ESL_XEXCEPTION(eslFAIL, "align_spliced_path FAIL");

  p7_Backward(pli->amino_sq->dsq, pli->amino_sq->L, om, pli->fwd, pli->bwd, NULL);  

  status = p7_Decoding(om, pli->fwd, pli->bwd, pli->bwd);
  if (status == eslERANGE) return eslFAIL; 
   
  p7_OptimalAccuracy(om, pli->bwd, pli->fwd, &oasc);
  p7_OATrace        (om, pli->bwd, pli->fwd, tr);

  p7_trace_Index(tr);

  //p7_trace_Dump(stdout, tr, NULL,NULL); 
  if((hit->dcl->tr = splice_trace_convert(tr, pli->orig_nuc_idx, &splice_cnt)) == NULL) goto ERROR; 

  if(splice_cnt == 0) return eslOK;

  hit->dcl->ad = create_spliced_alidisplay(hit->dcl->tr, 0, om, target_seq, pli->amino_sq, tr->sqfrom[0], splice_cnt);   
  
  p7_Null2_ByExpectation(om, pli->bwd, null2);
  domcorrection = 0.;
  for (i = 0; i < pli->amino_sq->n; i++)
    domcorrection += logf(null2[pli->amino_sq->dsq[i]]);

  hit->dcl->domcorrection = domcorrection;

  hit->dcl->ihmm = hit->dcl->ad->hmmfrom;
  hit->dcl->jhmm = hit->dcl->ad->hmmto;
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
    return status;
}


P7_TRACE *
splice_trace_convert(P7_TRACE *orig_tr, int *orig_nuc_idx, int *splice_cnt)
{

  int       z;
  int       curr_nuc_idx;
  int       prev_nuc_idx;
  int       sp_cnt;
  P7_TRACE *new_tr;
  
  if ((new_tr = p7_trace_splice_CreateWithPP()) == NULL) goto ERROR;
  sp_cnt = 0; 

  z = 0;
  prev_nuc_idx = orig_nuc_idx[(orig_tr->i[z]*3)];
  while (z < orig_tr->N) {
    curr_nuc_idx = orig_nuc_idx[(orig_tr->i[z]*3)];
    switch(orig_tr->st[z]) {
      case p7T_N: p7_trace_splice_AppendWithPP(new_tr, p7T_N, orig_tr->k[z], curr_nuc_idx, 3, -1, orig_tr->pp[z]); break;
      case p7T_C: p7_trace_splice_AppendWithPP(new_tr, p7T_C, orig_tr->k[z], curr_nuc_idx, 3, -1, orig_tr->pp[z]); break;
      case p7T_J: p7_trace_splice_AppendWithPP(new_tr, p7T_J, orig_tr->k[z], curr_nuc_idx, 3, -1, orig_tr->pp[z]); break;
      case p7T_M:
        /* Check if the sequence was spliced at this M position. If so, 
         * then determine the splice option and use the p7T_MS state*/  
        if (prev_nuc_idx > 1 && curr_nuc_idx > prev_nuc_idx+3) {
          if     (orig_nuc_idx[(orig_tr->i[z]*3)-2] - prev_nuc_idx > 1) {
            new_tr->sp[new_tr->N-1] = p7S_xxyyABC;
            p7_trace_splice_AppendWithPP(new_tr, p7T_R, 0,             prev_nuc_idx+2, 0, p7S_xxyyABC, 0.);
            p7_trace_splice_AppendWithPP(new_tr, p7T_P, 0,             curr_nuc_idx-4, 0, p7S_xxyyABC, 0.);
            p7_trace_splice_AppendWithPP(new_tr, p7T_A, orig_tr->k[z], curr_nuc_idx,   3, p7S_xxyyABC, orig_tr->pp[z]);
            sp_cnt++;
          }
          else if(orig_nuc_idx[(orig_tr->i[z]*3)-1] - orig_nuc_idx[(orig_tr->i[z]*3)-2] > 1) {
            p7_trace_splice_AppendWithPP(new_tr, p7T_R, 0,             prev_nuc_idx+3, 0, p7S_AxxyyBC, 0.);
            p7_trace_splice_AppendWithPP(new_tr, p7T_P, 0,             curr_nuc_idx-3, 0, p7S_AxxyyBC, 0.);
            p7_trace_splice_AppendWithPP(new_tr, p7T_A, orig_tr->k[z], curr_nuc_idx,   3, p7S_AxxyyBC, orig_tr->pp[z]);
            sp_cnt++;
          }
          else if(orig_nuc_idx[(orig_tr->i[z]*3)]   - orig_nuc_idx[(orig_tr->i[z]*3)-1] > 1) {
            p7_trace_splice_AppendWithPP(new_tr, p7T_R, orig_tr->k[z], prev_nuc_idx+4, 3, p7S_ABxxyyC, orig_tr->pp[z]);
            p7_trace_splice_AppendWithPP(new_tr, p7T_P, 0,             curr_nuc_idx-2, 0, p7S_ABxxyyC, 0.);
            p7_trace_splice_AppendWithPP(new_tr, p7T_A, 0,             curr_nuc_idx,   0, p7S_ABxxyyC, 0.);
            sp_cnt++;
          }
        }
        else { 
          p7_trace_splice_AppendWithPP(new_tr, p7T_M, orig_tr->k[z], curr_nuc_idx, 3, -1, orig_tr->pp[z]); 
         }
        break;
      case p7T_I: 
        /* Check if the sequence was spliced at this I position. If so,
         * then determine the splice option and use the p7T_IS state */
        if(prev_nuc_idx > 1 && curr_nuc_idx > prev_nuc_idx+3) {
          if     (orig_nuc_idx[(orig_tr->i[z]*3)-2] - prev_nuc_idx > 1) {
            new_tr->sp[new_tr->N-1] = p7S_xxyyABC;
            p7_trace_splice_AppendWithPP(new_tr, p7T_R, 0,             prev_nuc_idx+2, 0, p7S_xxyyABC, 0.);
            p7_trace_splice_AppendWithPP(new_tr, p7T_P, 0,             curr_nuc_idx-4, 0, p7S_xxyyABC, 0.);
            p7_trace_splice_AppendWithPP(new_tr, p7T_A, orig_tr->k[z], curr_nuc_idx,   3, p7S_xxyyABC, orig_tr->pp[z]);
            sp_cnt++;
          }
          else if(orig_nuc_idx[(orig_tr->i[z]*3)-1] - orig_nuc_idx[(orig_tr->i[z]*3)-2] > 1) {
            p7_trace_splice_AppendWithPP(new_tr, p7T_R, 0,             prev_nuc_idx+3, 0, p7S_AxxyyBC, 0.);
            p7_trace_splice_AppendWithPP(new_tr, p7T_P, 0,             curr_nuc_idx-3, 0, p7S_AxxyyBC, 0.);
            p7_trace_splice_AppendWithPP(new_tr, p7T_A, orig_tr->k[z], curr_nuc_idx,   3, p7S_AxxyyBC, orig_tr->pp[z]);
            sp_cnt++;
          }
          else if(orig_nuc_idx[(orig_tr->i[z]*3)]   - orig_nuc_idx[(orig_tr->i[z]*3)-1] > 1) {
            p7_trace_splice_AppendWithPP(new_tr, p7T_R, orig_tr->k[z], prev_nuc_idx+4, 3, p7S_ABxxyyC, orig_tr->pp[z]);
            p7_trace_splice_AppendWithPP(new_tr, p7T_P, 0,             curr_nuc_idx-2, 0, p7S_ABxxyyC, 0.);
            p7_trace_splice_AppendWithPP(new_tr, p7T_A, 0,             curr_nuc_idx,   0, p7S_ABxxyyC, 0.);
            sp_cnt++;
          } 
        }
        else
          p7_trace_splice_AppendWithPP(new_tr, p7T_I, orig_tr->k[z], curr_nuc_idx, 0, -1, orig_tr->pp[z]); 
        break; 
      case p7T_D: p7_trace_splice_AppendWithPP(new_tr, p7T_D, orig_tr->k[z], curr_nuc_idx, 0, -1, orig_tr->pp[z]); break;
      case p7T_X: p7_trace_splice_AppendWithPP(new_tr, p7T_X, orig_tr->k[z], curr_nuc_idx, 0, -1, orig_tr->pp[z]); break;
      case p7T_S: p7_trace_splice_AppendWithPP(new_tr, p7T_S, orig_tr->k[z], curr_nuc_idx, 0, -1, orig_tr->pp[z]); break;
      case p7T_B: p7_trace_splice_AppendWithPP(new_tr, p7T_B, orig_tr->k[z], curr_nuc_idx, 0, -1, orig_tr->pp[z]); break;
      case p7T_E: p7_trace_splice_AppendWithPP(new_tr, p7T_E, orig_tr->k[z], curr_nuc_idx, 0, -1, orig_tr->pp[z]); break;
      case p7T_T: p7_trace_splice_AppendWithPP(new_tr, p7T_T, orig_tr->k[z], curr_nuc_idx, 0, -1, orig_tr->pp[z]); break;
      default:    return NULL; 
    }   
    prev_nuc_idx = curr_nuc_idx;
    z++;
  } 
  *splice_cnt = sp_cnt; 
  return new_tr;

  ERROR:
    return NULL;

}



/* Function:  create_spliced_alidisplay() 
 * Synopsis:  Create an alignment display, from trace and oprofile.
 *
 * Purpose:   Creates and returns a BATH formated and spliced alignment 
 *            display for domain number <which> in traceback <tr>, 
 *            where the traceback corresponds to an alignment of 
 *            optimized profile <om> to digital sequence <dsq>, and the 
 *            unique name of that target sequence <dsq> is <sqname>. 
 *            The <which> index starts at 0.
 *
 *            It will be a little faster if the trace is indexed with
 *            <p7_trace_Index()> first. The number of domains is then
 *            in <tr->ndom>. If the caller wants to create alidisplays
 *            for all of these, it would loop <which> from
 *            <0..tr->ndom-1>.
 *
 *            However, even without an index, the routine will work fine.
 *
 * Args:      tr           - traceback
 *            which        - domain number, 0..tr->ndom-1
 *            om           - optimized profile (query)
 *            target_seq   - digital nucleotide sequence (unspliced target)
 *            amino_sq     - digital amino sequence (translation of spliced nucleotides 
 *            orig_nuc_idx - array of indicies in <target_seq> that correspond to the spliced <tr->i> indicies 
 *            amino_pos    - first position in the alignmant for <amino_sq>
 *            exon_cnt     - the total number of spliced exons - 1. 
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <NULL> on allocation failure, or if something's internally corrupt
 *            in the data.
 */
P7_ALIDISPLAY *
create_spliced_alidisplay(const P7_TRACE *tr, int which, const P7_OPROFILE *om, const ESL_SQ *target_seq, const ESL_SQ *amino_sq, int amino_pos, int splice_cnt)
{


  int            z1, z2;
  int            n, pos;
  int            x,y,z;
  int            k,s,c,p; 
  int64_t        i;
  int            sq_namelen,  sq_acclen,  sq_desclen;
  int            hmm_namelen, hmm_acclen, hmm_desclen;
  int            orf_namelen;
  char          *alphaAmino;
  char          *alphaDNA;  
  P7_ALIDISPLAY *ad;        
  int            status;


  alphaAmino = om->abc->sym;
  alphaDNA   = target_seq->abc->sym;
  /* First figure out which piece of the trace (from first match to last match)
   * we're going to represent, and how big it is.
   */

  if (tr->ndom > 0) {       /* if we have an index, this is a little faster: */
    for (z1 = tr->tfrom[which]; z1 < tr->N; z1++) if (tr->st[z1] == p7T_M) break;  /* find next M state    */
    if (z1 == tr->N) return NULL;                                                  /* no M? corrupt trace    */
    for (z2 = tr->tto[which];   z2 >= 0 ;   z2--) if (tr->st[z2] == p7T_M) break;  /* find prev M state      */
    if (z2 == -1) return NULL;                                                     /* no M? corrupt trace    */
  } else {          /* without an index, we can still do it fine:    */

    for (z1 = 0; which >= 0 && z1 < tr->N; z1++) if (tr->st[z1] == p7T_B) which--; /* find the right B state */
    if (z1 == tr->N) return NULL;                                                  /* no such domain <which> */
    for (; z1 < tr->N; z1++) if (tr->st[z1] == p7T_M) break;                       /* find next M state      */
    if (z1 == tr->N) return NULL;                                                  /* no M? corrupt trace    */
    for (z2 = z1; z2 < tr->N; z2++) if (tr->st[z2] == p7T_E) break;                /* find the next E state  */
    for (; z2 >= 0;    z2--) if (tr->st[z2] == p7T_M) break;                       /* find prev M state      */
    if (z2 == -1) return NULL;                                                     /* no M? corrupt trace    */
  }

  /* Allocate enough space for each M, D, and I trace position plus up to 3 for every exon splicing */
  n = (z2-z1+2) * 3;                        /* model, mline, aseq mandatory         */
  n += 5*(z2-z1+1)+1;                       /* nucleotide sequence                  */
  n += (z2-z1+2);                           /* codon lengths                        */
  if (om->rf[0]  != 0)    n += (z2-z1+2);   /* optional reference line              */
  if (om->mm[0]  != 0)    n += (z2-z1+2);   /* optional reference line              */
  if (om->cs[0]  != 0)    n += (z2-z1+2);   /* optional structure line              */
  if (tr->pp     != NULL) n += (z2-z1+2);   /* optional posterior prob line         */
  hmm_namelen = strlen(om->name);                           n += hmm_namelen + 1;
  hmm_acclen  = (om->acc  != NULL ? strlen(om->acc)  : 0);  n += hmm_acclen  + 1;
  hmm_desclen = (om->desc != NULL ? strlen(om->desc) : 0);  n += hmm_desclen + 1;
  sq_namelen  = strlen(target_seq->name);                   n += sq_namelen  + 1;
  sq_acclen   = strlen(target_seq->acc);                    n += sq_acclen   + 1; /* sq->acc is "\0" when unset */
  sq_desclen  = strlen(target_seq->desc);                   n += sq_desclen  + 1; /* same for desc              */
  orf_namelen = strlen(target_seq->orfid);                  n += orf_namelen + 1; /* same for orfname          */
  
  ad = NULL;
  ESL_ALLOC(ad, sizeof(P7_ALIDISPLAY));
  ad->mem = NULL;

  pos = 0;
  ad->memsize = sizeof(char) * n;
  ESL_ALLOC(ad->mem, ad->memsize);

  if (om->rf[0]  != 0) { ad->rfline = ad->mem + pos; pos += z2-z1+2; } else { ad->rfline = NULL; }
  ad->mmline = NULL;
  if (om->cs[0]  != 0) { ad->csline = ad->mem + pos; pos += z2-z1+2; } else { ad->csline = NULL; }
  ad->model   = ad->mem + pos;  pos += z2-z1+2;
  ad->mline   = ad->mem + pos;  pos += z2-z1+2;
  ad->aseq    = ad->mem + pos;  pos += z2-z1+2;
  ad->codon   = ad->mem + pos;  pos += z2-z1+2;
  ad->ntseq  = ad->mem + pos;  pos += 5*(z2-z1+1)+1; /* for the nucleotide sequence there will be 5 times as many bytes */
  if (tr->pp != NULL)  { ad->ppline = ad->mem + pos;  pos += z2-z1+2;} else { ad->ppline = NULL; }
  ad->hmmname = ad->mem + pos;  pos += hmm_namelen +1;
  ad->hmmacc  = ad->mem + pos;  pos += hmm_acclen +1;
  ad->hmmdesc = ad->mem + pos;  pos += hmm_desclen +1;
  ad->sqname  = ad->mem + pos;  pos += sq_namelen +1;
  ad->sqacc   = ad->mem + pos;  pos += sq_acclen +1;
  ad->sqdesc  = ad->mem + pos;  pos += sq_desclen +1;
  ad->orfname = ad->mem + pos;  pos += orf_namelen +1;

  strcpy(ad->hmmname, om->name);
  if (om->acc  != NULL) strcpy(ad->hmmacc,  om->acc);  else ad->hmmacc[0]  = 0;
  if (om->desc != NULL) strcpy(ad->hmmdesc, om->desc); else ad->hmmdesc[0] = 0;

  strcpy(ad->sqname,  target_seq->name);
  strcpy(ad->sqacc,   target_seq->acc);
  strcpy(ad->sqdesc,  target_seq->desc);
  strcpy(ad->orfname, target_seq->orfid);

  ESL_ALLOC(ad->exon_starts, sizeof(int64_t) * splice_cnt);

  /* Determine hit coords */
  ad->hmmfrom     = tr->k[z1];
  ad->hmmto       = tr->k[z2];
  tr->hmmfrom[0]  = tr->k[z1];   //TODO: Delete
  tr->hmmto[0]    = tr->k[z2];   //TODO: Delete
  ad->M           = om->M;
  ad->frameshifts = 0;
  ad->stops       = 0;

  if(target_seq->start < target_seq->end) {
    ad->sqfrom  = tr->i[z1] + target_seq->start - 1 - (tr->c[z1] - 1);
    ad->sqto    = tr->i[z2] + target_seq->start - 1;
  } else {
    ad->sqto    = target_seq->n - tr->i[z2] + target_seq->end;
    ad->sqfrom  = target_seq->n - tr->i[z1] + target_seq->end;
  }
  
  ad->L       = target_seq->n;

  /* optional rf line */
  if (ad->rfline != NULL) {
    for (z = z1; z <= z2; z++) {
      if(tr->st[z] == p7T_R || tr->st[z] == p7T_A) ad->rfline[z-z1] = ' ';
      else if(tr->st[z] == p7T_I)                  ad->rfline[z-z1] = '.';
      else                                         ad->rfline[z-z1] = om->rf[tr->k[z]];
    }
    ad->rfline[z-z1] = '\0';
  }
  /* optional mm line */
  if (ad->mmline != NULL) {
    for (z = z1; z <= z2; z++) {
      if(tr->st[z] == p7T_R || tr->st[z] == p7T_A) ad->mmline[z-z1] = ' ';
      else if(tr->st[z] == p7T_I)                  ad->mmline[z-z1] = '.';
      else                                         ad->mmline[z-z1] = om->mm[tr->k[z]];
    }
    ad->mmline[z-z1] = '\0';
  }
  /* optional cs line */
  if (ad->csline != NULL) {
    for (z = z1; z <= z2; z++) {
      if(tr->st[z] == p7T_R || tr->st[z] == p7T_A) ad->csline[z-z1] = ' ';
      else if(tr->st[z] == p7T_I)                  ad->csline[z-z1] = '.';
      else                                         ad->csline[z-z1] = om->cs[tr->k[z]];
    }
    ad->csline[z-z1] = '\0';
  }

  if (ad->ppline != NULL) {
    for (z = z1; z <= z2; z++) {
      if      (tr->st[z] == p7T_D)                              ad->ppline[z-z1] = '.';
      else if (tr->st[z] == p7T_P)                              ad->ppline[z-z1] = ' ';
      else if (tr->st[z] == p7T_R && tr->sp[z] == p7S_ABxxyyC)  ad->ppline[z-z1] = p7_alidisplay_EncodePostProb(tr->pp[z]); 
      else if (tr->st[z] == p7T_R)                              ad->ppline[z-z1] = ' ';
      else if (tr->st[z] == p7T_A && tr->sp[z] == p7S_ABxxyyC)  ad->ppline[z-z1] = ' ';
      else if (tr->st[z] == p7T_A)                              ad->ppline[z-z1] = p7_alidisplay_EncodePostProb(tr->pp[z]); 
      else                                                      ad->ppline[z-z1] = p7_alidisplay_EncodePostProb(tr->pp[z]);
    }
    ad->ppline[z-z1] = '\0';
  }

  /* There are three ways that the splice signals can apprear in our aligment        */
  /* xx is the doner splice signal and yy is the acceptor splice signal              */
  /* A,B, and C are the nucleotides of the codon surrounding the splice signals      */
  /* N represents the nucleotides present in the last state beofre the splice signal */
  /* $ is used to let p7_alidisplay_Print() know hoew to display the splice signals  */
  /*                 p7T_M      p7T_R      p7T_P      p7T_A                          */
  /* p7S_xxyyABC    " NNNx"    "x    "    "$$  y"    "yABC"                          */
  /* p7S_AxxyyBC    " NNN "    " Axx "    "$ $  "    "yyBC "                         */
  /* p7S_ABxxyyC    " NNN "    " ABxx"    "$  $ "    " yyC "                         */


  /* mandatory three alignment display lines: model, mline, aseq */
  x = y = 0;
  for (z = z1; z <= z2; z++)
  {
    k = tr->k[z];
    i = tr->i[z];
    s = tr->st[z];
    p = tr->sp[z];

    switch (s) {
      case p7T_M:
        ad->codon[y] = 3;
        ad->model[z-z1] = om->consensus[k];
        ad->aseq[z-z1]  = toupper(alphaAmino[amino_sq->dsq[amino_pos]]);
        if      (amino_sq->dsq[amino_pos] == esl_abc_DigitizeSymbol(om->abc, om->consensus[k]))
          ad->mline[z-z1] = ad->model[z-z1];
        else if (p7_oprofile_FGetEmission(om, k, amino_sq->dsq[amino_pos]) > 1.0)
          ad->mline[z-z1] = '+'; /* >1 not >0; om has odds ratios, not scores */
        else
          ad->mline[z-z1] = ' ';
        amino_pos++;
      
        ad->ntseq [5*(z-z1)]   = ' ';
        ad->ntseq [5*(z-z1)+1] = toupper(alphaDNA[target_seq->dsq[i-2]]);
        ad->ntseq [5*(z-z1)+2] = toupper(alphaDNA[target_seq->dsq[i-1]]);
        ad->ntseq [5*(z-z1)+3] = toupper(alphaDNA[target_seq->dsq[i]]);
        if (p == p7S_xxyyABC ) {
          ad->ntseq [5*(z-z1)+4] = tolower(alphaDNA[target_seq->dsq[i+1]]);
          ad->codon[y]++;
        }
        else
          ad->ntseq [5*(z-z1)+4] = ' ';

               break;

      case p7T_I:
        ad->codon[y] = 3;
        ad->model [z-z1] = '.';
        ad->aseq[z-z1] = toupper(alphaAmino[amino_sq->dsq[amino_pos]]);
        ad->mline [z-z1] = ' ';
        amino_pos++;

        ad->ntseq [5*(z-z1)]   = ' ';
        ad->ntseq [5*(z-z1)+1] = toupper(alphaDNA[target_seq->dsq[i-2]]);
        ad->ntseq [5*(z-z1)+2] = toupper(alphaDNA[target_seq->dsq[i-1]]);
        ad->ntseq [5*(z-z1)+3] = toupper(alphaDNA[target_seq->dsq[i]]);
        if (p == p7S_xxyyABC ) {
          ad->ntseq [5*(z-z1)+4] = tolower(alphaDNA[target_seq->dsq[i+1]]);
          ad->codon[y]++;
        }
        else
          ad->ntseq [5*(z-z1)+4] = ' ';      
  
        break;

      case p7T_D:
        ad->codon[y] = 0;
        ad->model [z-z1] = om->consensus[k];
        ad->mline [z-z1] = ' ';
        ad->aseq  [z-z1] = '-';
      
        ad->ntseq [5*(z-z1)] = ' ';
        ad->ntseq [5*(z-z1)+1] = '-';
        ad->ntseq [5*(z-z1)+2] = '-';
        ad->ntseq [5*(z-z1)+3] = '-';
        if (p == p7S_xxyyABC ) {
          ad->ntseq [5*(z-z1)+4] = tolower(alphaDNA[target_seq->dsq[i+1]]);
          ad->codon[y]++;
        }
        else
          ad->ntseq [5*(z-z1)+4] = ' ';
        break;

      case p7T_R:
                
        if(p == p7S_xxyyABC) {
          ad->model [z-z1] = ' ';
          ad->mline [z-z1] = ' ';
          ad->aseq  [z-z1] = ' ';
          ad->codon[y] = 1;

          ad->ntseq [5*(z-z1)]   = tolower(alphaDNA[target_seq->dsq[i]]);
          ad->ntseq [5*(z-z1)+1] = ' ';
          ad->ntseq [5*(z-z1)+2] = ' '; 
          ad->ntseq [5*(z-z1)+3] = ' ';
          ad->ntseq [5*(z-z1)+4] = ' ';
        }
        else if(p == p7S_AxxyyBC) {
          ad->model [z-z1] = ' ';
          ad->mline [z-z1] = ' ';
          ad->aseq  [z-z1] = ' ';
          ad->codon[y] = 3;

          ad->ntseq [5*(z-z1)]   = ' '; 
          ad->ntseq [5*(z-z1)+1] = toupper(alphaDNA[target_seq->dsq[i-2]]); 
          ad->ntseq [5*(z-z1)+2] = tolower(alphaDNA[target_seq->dsq[i-1]]);
          ad->ntseq [5*(z-z1)+3] = tolower(alphaDNA[target_seq->dsq[i]]);
          ad->ntseq [5*(z-z1)+4] = ' ';
        }
        else if(p == p7S_ABxxyyC) {
          ad->model[z-z1] = om->consensus[k];
          ad->aseq[z-z1]  = toupper(alphaAmino[amino_sq->dsq[amino_pos]]);
          if      (amino_sq->dsq[amino_pos] == esl_abc_DigitizeSymbol(om->abc, om->consensus[k]))
            ad->mline[z-z1] = ad->model[z-z1];
          else if (p7_oprofile_FGetEmission(om, k, amino_sq->dsq[amino_pos]) > 1.0)
            ad->mline[z-z1] = '+'; /* >1 not >0; om has odds ratios, not scores */
          else
           ad->mline[z-z1] = ' ';
          amino_pos++;
          ad->codon[y] = 4;

          ad->ntseq [5*(z-z1)]   = ' ';
          ad->ntseq [5*(z-z1)+1] = toupper(alphaDNA[target_seq->dsq[i-3]]);
          ad->ntseq [5*(z-z1)+2] = toupper(alphaDNA[target_seq->dsq[i-2]]);
          ad->ntseq [5*(z-z1)+3] = tolower(alphaDNA[target_seq->dsq[i-1]]);
          ad->ntseq [5*(z-z1)+4] = tolower(alphaDNA[target_seq->dsq[i]]);
        }
        break;
      case p7T_P:
        ad->model [z-z1] = ' ';
        ad->mline [z-z1] = ' ';
        ad->aseq  [z-z1] = ' ';
        if(target_seq->start < target_seq->end)
          ad->exon_starts[x] = i + target_seq->start - 1; 
        else
          ad->exon_starts[x] = target_seq->n-i+target_seq->end;
        x++;
       
        if(p == p7S_xxyyABC) {
          ad->codon[y] = 1;
          ad->ntseq [5*(z-z1)]   = '$';
          ad->ntseq [5*(z-z1)+1] = '$';
          ad->ntseq [5*(z-z1)+2] = ' ';
          ad->ntseq [5*(z-z1)+3] = ' ';
          ad->ntseq [5*(z-z1)+4] = tolower(alphaDNA[target_seq->dsq[i]]);;
        }
        else if(p == p7S_AxxyyBC) {
          ad->codon[y] = 0;
          ad->ntseq [5*(z-z1)]   = '$';
          ad->ntseq [5*(z-z1)+1] = ' ';
          ad->ntseq [5*(z-z1)+2] = '$';
          ad->ntseq [5*(z-z1)+3] = ' ';
          ad->ntseq [5*(z-z1)+4] = ' ';
        }
        else if(p == p7S_ABxxyyC) {
          ad->codon[y] = 0;
          ad->ntseq [5*(z-z1)]   = '$';
          ad->ntseq [5*(z-z1)+1] = ' ';
          ad->ntseq [5*(z-z1)+2] = ' ';
          ad->ntseq [5*(z-z1)+3] = '$';
          ad->ntseq [5*(z-z1)+4] = ' ';
        }
        break;
      case p7T_A:
        if(p == p7S_xxyyABC) {
          ad->model[z-z1] = om->consensus[k];
          ad->aseq[z-z1]  = toupper(alphaAmino[amino_sq->dsq[amino_pos]]);
          if      (amino_sq->dsq[amino_pos] == esl_abc_DigitizeSymbol(om->abc, om->consensus[k]))
            ad->mline[z-z1] = ad->model[z-z1];
          else if (p7_oprofile_FGetEmission(om, k, amino_sq->dsq[amino_pos]) > 1.0)
            ad->mline[z-z1] = '+'; /* >1 not >0; om has odds ratios, not scores */
          else
           ad->mline[z-z1] = ' ';
          amino_pos++;       
          
          ad->codon[y] = 4;
          ad->ntseq [5*(z-z1)]   = tolower(alphaDNA[target_seq->dsq[i-3]]); 
          ad->ntseq [5*(z-z1)+1] = toupper(alphaDNA[target_seq->dsq[i-2]]);
          ad->ntseq [5*(z-z1)+2] = toupper(alphaDNA[target_seq->dsq[i-1]]);
          ad->ntseq [5*(z-z1)+3] = toupper(alphaDNA[target_seq->dsq[i]]);
          ad->ntseq [5*(z-z1)+4] = ' ';
        }
        else if(p == p7S_AxxyyBC) {
          ad->model[z-z1] = om->consensus[k];
          ad->aseq[z-z1]  = toupper(alphaAmino[amino_sq->dsq[amino_pos]]);
          if      (amino_sq->dsq[amino_pos] == esl_abc_DigitizeSymbol(om->abc, om->consensus[k]))
            ad->mline[z-z1] = ad->model[z-z1];
          else if (p7_oprofile_FGetEmission(om, k, amino_sq->dsq[amino_pos]) > 1.0)
            ad->mline[z-z1] = '+'; /* >1 not >0; om has odds ratios, not scores */
          else
           ad->mline[z-z1] = ' ';
          amino_pos++;
           
          ad->codon[y] = 4;
          ad->ntseq [5*(z-z1)]   = tolower(alphaDNA[target_seq->dsq[i-3]]);
          ad->ntseq [5*(z-z1)+1] = tolower(alphaDNA[target_seq->dsq[i-2]]);
          ad->ntseq [5*(z-z1)+2] = toupper(alphaDNA[target_seq->dsq[i-1]]);
          ad->ntseq [5*(z-z1)+3] = toupper(alphaDNA[target_seq->dsq[i]]);
          ad->ntseq [5*(z-z1)+4] = ' ';
        }
        else if(p == p7S_ABxxyyC) {
          ad->model [z-z1] = ' ';
          ad->mline [z-z1] = ' ';
          ad->aseq  [z-z1] = ' '; 

          ad->codon[y] = 3;
          ad->ntseq [5*(z-z1)]   = ' ';
          ad->ntseq [5*(z-z1)+1] = tolower(alphaDNA[target_seq->dsq[i-2]]);
          ad->ntseq [5*(z-z1)+2] = tolower(alphaDNA[target_seq->dsq[i-1]]);
          ad->ntseq [5*(z-z1)+3] =  toupper(alphaDNA[target_seq->dsq[i]]);
          ad->ntseq [5*(z-z1)+4] = ' ';
        }
        break;
      default: ESL_XEXCEPTION(eslEINVAL, "invalid state in trace: not M,D,I");
    } 
    y++;
  }

  ad->model [z2-z1+1] = '\0';
  ad->mline [z2-z1+1] = '\0';
  ad->aseq  [z2-z1+1] = '\0';
  ad->ntseq  [5*(z2-z1+1)] = '\0';
  ad->N = z2-z1+1;

  return ad;

  ERROR:
    p7_alidisplay_Destroy(ad);
    return NULL;
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
    fprintf(fp, "\n   Hit   hmm_from   hmm_to   seq_from     seq_to\n");

    for(i = 0; i < target_range->th->N; i++) {

      hit = target_range->th->hit[i];     

   	  fprintf(fp, "   %d  %6d  %9d %9" PRId64 " %12" PRId64 "\n", 
                      i+1, hit->dcl->ihmm, hit->dcl->jhmm, hit->dcl->iali,            hit->dcl->jali);	 
   }
  }
  
  fprintf(fp, "\n"); 

  return;
}

void
graph_dump(FILE *fp, SPLICE_GRAPH *graph, P7_TOPHITS *th, ESL_SQ *target_seq, int print_edges) 
{


  int          i,j;
  int          edge_id;
  int          nuc_end, nuc_start;
  int          amino_start, amino_end;
  int          curr_node, next_node;
  int          curr_edge;
  P7_HIT      *node_hit;
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

    fprintf(fp, "\n Hit Scores  \n");
    for(i = 0; i < graph->num_nodes; i++ ) 
      fprintf(fp, "%8d", i+1);
    
    fprintf(fp, "\n");

    for(i = 0; i < graph->num_nodes; i++ ) 
      fprintf(fp, "%8.2f", graph->hit_scores[i]);
    
    fprintf(fp, "\n");

    fprintf(fp, "\n Path Scores  \n");
    for(i = 0; i < graph->num_nodes; i++ )
      fprintf(fp, "%8d", i+1);

    fprintf(fp, "\n");

    for(i = 0; i < graph->num_nodes; i++ )
      fprintf(fp, "%8.2f", graph->path_scores[i]);

    fprintf(fp, "\n");

     fprintf(fp, "\n Best Edge  \n");
    for(i = 0; i < graph->num_nodes; i++ )
      fprintf(fp, "%8d", i+1);

    fprintf(fp, "\n");

    for(i = 0; i < graph->num_nodes; i++ )
      fprintf(fp, "%8d", graph->best_out_edge[i]+1);

    fprintf(fp, "\n");

  } 

  fprintf(fp, " Graph has full length path : %s\n", (graph->has_full_path ? "YES" : "NO"));
  
  if(graph->has_full_path) {
    fprintf(fp, "\n Full Path  \n");
   
    curr_node   = graph->best_full_path_start; 
    node_hit    = th->hit[curr_node];
    nuc_start   = node_hit->dcl->iali;
    amino_start = node_hit->dcl->ihmm;
    next_node   = graph->best_out_edge[curr_node];
    curr_edge   = graph->edge_id[curr_node][next_node];
    edge        = graph->edges[curr_edge];
    nuc_end     = edge->upstream_spliced_nuc_end;
    amino_end   = edge->upstream_spliced_amino_end;

    fprintf(fp, "START -> node %3d hmm coords : %5d %5d  seq coords : %10d %10d -> ", curr_node+1, amino_start, amino_end, nuc_start, nuc_end);

    for(i = 1; i < graph->best_full_path_length; i++) {

      nuc_start   = edge->downstream_spliced_nuc_start;
      amino_start = edge->downstream_spliced_amino_start;
      next_node   = graph->best_out_edge[curr_node];
      curr_edge   = graph->edge_id[curr_node][next_node];
      edge        = graph->edges[curr_edge];
      nuc_end     = edge->upstream_spliced_nuc_end;
      amino_end   = edge->upstream_spliced_amino_end;
      curr_node   = next_node;

      fprintf(fp, "\n         node %3d hmm coords : %5d %5d  seq coords : %10d %10d ->  ", curr_node+1, amino_start, amino_end, nuc_start, nuc_end);
    } 
    fprintf(fp, "END\n");
  }


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
           fprintf(fp, "      Upsteam Node End Coords:     %5d  %10" PRId64 "\n", edge->upstream_spliced_amino_end, nuc_end); 
           fprintf(fp, "      Downsteam Node Start Coords: %5d  %10" PRId64 "\n", edge->downstream_spliced_amino_start, nuc_start);
           fprintf(fp, "\n");
        }
      }
    }
  }
  fprintf(fp, "\n");
 

  return;

}
