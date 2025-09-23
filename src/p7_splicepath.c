/* SPLICE_PATH: a set of spliced hits 
 *
 * Contents:
 *    1. The SPLICE_PATH object.
 *    2. Path Finding Algorithms
 *    3. Debugging tools.
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

static int longest_path_upstream (SPLICE_GRAPH *graph);
static int topological_sort_upstream(SPLICE_GRAPH *graph, int *visited, int *stack, int *stack_size, int node);
static SPLICE_PATH* check_bypass(SPLICE_GRAPH *graph, SPLICE_PATH *path);
static SPLICE_PATH* check_bypass_unspliced(SPLICE_GRAPH *graph, SPLICE_PATH *path);
static SPLICE_PATH* check_bypass_extend(SPLICE_GRAPH *orig_graph, SPLICE_GRAPH *extend_graph, SPLICE_PATH *path);
static float get_sub_path_score(SPLICE_GRAPH *graph, int source_node, int termination_node);
static int hit_between(P7_DOMAIN *upstream, P7_DOMAIN *middle, P7_DOMAIN *downstream, int revcomp);
static int has_out_edge(SPLICE_GRAPH *graph, int node_id);

/*****************************************************************
 * 1. The SPLICE_PATH structure.
 *****************************************************************/

/* Function:  p7_splicepath_Create()
 *
 * Purpose:   Allocates a splice path with room for 2x <path_len> 
 *            hits.
 *
 * Returns:   a pointer to the new <SPLICE_PATH> structure
 *            on success.
 *
 * Throws:    <NULL> on allocation error.
 */
SPLICE_PATH*
p7_splicepath_Create(int path_len)
{

  SPLICE_PATH *path;
  int          status;

  path = NULL;
  ESL_ALLOC(path, sizeof(SPLICE_PATH));

  path->frameshift = FALSE;

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
    p7_splicepath_Destroy(path);
    return NULL;
}

/* Function:  p7_splicepath_Grow()
 * Synopsis:  Reallocates a larger splice path, if needed.
 *
 * Purpose:   If <SPLICE_PATH> cannot hold another hit,
 *            doubles the internal allocation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 */
int
p7_splicepath_Grow(SPLICE_PATH *path)
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
    p7_splicepath_Destroy(path);
    return status;

}

/* Function:  p7_splicepath_Insert()
 * Synopsis:  Insert a new path hit between two existing hits. 
 *
 * Purpose:   Add new <P7_HIT> to existing <SPLICE_PATH> one 
 *            step downstream of the provided upstream idx
 *
 * Returns:   <eslOK>.
 *
 */
int
p7_splicepath_Insert(SPLICE_PATH *path, P7_HIT *new_hit, float edge_score, int step)
{

  int s;
  
  p7_splicepath_Grow(path); 

  for(s = path->path_len; s > step; s--) {

    path->hits[s]        = path->hits[s-1];
    path->hit_scores[s]  = path->hit_scores[s-1];
    path->edge_scores[s] = path->edge_scores[s-1];
    path->node_id[s]     = path->node_id[s-1];

    path->downstream_spliced_amino_start[s] = path->downstream_spliced_amino_start[s-1];
    path->upstream_spliced_amino_end[s+1]     = path->upstream_spliced_amino_end[s];
    path->downstream_spliced_nuc_start[s]   = path->downstream_spliced_nuc_start[s-1];
    path->upstream_spliced_nuc_end[s+1]       = path->upstream_spliced_nuc_end[s];

  }

  path->path_len++;
 
  path->hits[step] = new_hit;
  
  path->hit_scores[step]    = new_hit->dcl->aliscore;
  path->edge_scores[step]   = edge_score;
  path->signal_scores[step] = 0.; 

  path->node_id[step] = -1;
  path->split[step]   = FALSE;
  
  path->downstream_spliced_amino_start[step] = new_hit->dcl->ihmm;
  path->downstream_spliced_nuc_start[step]   = new_hit->dcl->iali;

  path->upstream_spliced_amino_end[step+1]     = new_hit->dcl->jhmm;
  path->upstream_spliced_nuc_end[step+1]       = new_hit->dcl->jali;

  return eslOK;
 
}

/* Function: p7_splicepath_Destroy()
 *
 * Purpose:  Frees a <SPLICE_PATH>
 */
void
p7_splicepath_Destroy(SPLICE_PATH *path)
{

   if(path == NULL) return;

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

/*****************************************************************
 * 3. Path Finding Algorithms
 *****************************************************************/

SPLICE_PATH*
p7_splicepath_GetBestPath_Unspliced(SPLICE_GRAPH *graph)
{

  int          i;
  int          path_len;
  int          start_node;
  int          curr_node;
  int          next_node;
  int          step_cnt;
  int          contains_orig;
  float        best_start_score;
  SPLICE_EDGE *out_edge;
  SPLICE_PATH *path;
  P7_TOPHITS  *th;
  int         status;

  th = graph->th;
  contains_orig = FALSE;

  /* Find best scoreing paths */
  if((status = longest_path_upstream(graph)) != eslOK) goto ERROR;

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
  
    /* We have run out of paths */ 
    if(start_node < 0) return NULL; 
 
    curr_node = start_node;
  
    path_len  = 1;
  
    out_edge = NULL;
    
    /* Check each set of niehgbring edges to ensure they are not "backwards" 
       (upstream splice is not downstream of downstream splice) */
    while(graph->best_out_edge[curr_node] >= 0) {

      if(curr_node < graph->orig_N) contains_orig = TRUE;

      next_node = graph->best_out_edge[curr_node];
      
      out_edge = p7_splicegraph_GetEdge(graph, curr_node, next_node);
      if(out_edge == NULL || out_edge->splice_score == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "Edge does not exist");
  
      curr_node = next_node;
      path_len++;
    } 
  
    if(curr_node < graph->orig_N) contains_orig = TRUE;  

    if(!contains_orig) graph->path_scores[start_node] = -eslINFINITY;
  }

  /* Once a viable path has been found build the path */
  if ((path = p7_splicepath_Create(path_len)) == NULL) goto ERROR;
    
  path->revcomp = graph->revcomp;
  
  path->hit_scores[0]    = th->hit[start_node]->dcl->aliscore;
  path->edge_scores[0]   = 0.;
  path->signal_scores[0] = 0.;

  path->node_id[0] = start_node;
  
  path->split[0]   = FALSE;
  path->hits[0]    = th->hit[start_node];

  path->downstream_spliced_amino_start[0] = path->hits[0]->dcl->ihmm;
  path->downstream_spliced_nuc_start[0]   = path->hits[0]->dcl->iali;

  if(th->hit[start_node]->dcl->tr->fs) path->frameshift = TRUE; //{ printf("start_node %d\n" ,start_node); fflush(stdout); path->frameshift = TRUE; }
  
  curr_node = start_node;
  step_cnt = 1;
  while (step_cnt < path_len) {
    next_node = graph->best_out_edge[curr_node];
    out_edge = p7_splicegraph_GetEdge(graph, curr_node, next_node);

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

    if(th->hit[next_node]->dcl->tr->fs) path->frameshift = TRUE; //{ printf("next_node %d\n" ,next_node); fflush(stdout); path->frameshift = TRUE; }    

    curr_node = next_node;
    step_cnt++;
  }

  path->upstream_spliced_amino_end[step_cnt] = path->hits[step_cnt-1]->dcl->jhmm;
  path->upstream_spliced_nuc_end[step_cnt]   = path->hits[step_cnt-1]->dcl->jali;

  path = check_bypass_unspliced(graph, path);

  return path;

  ERROR:
    if(path != NULL) p7_splicepath_Destroy(path);
    return NULL; 

}

SPLICE_PATH*
p7_splicepath_GetBestPath_Extension(SPLICE_GRAPH *orig_graph, SPLICE_GRAPH *extend_graph)
{

  int          i;
  int          path_len;
  int          start_node;
  int          curr_node;
  int          next_node;
  int          step_cnt;
  int          contains_orig;
  float        best_start_score;
  SPLICE_EDGE *out_edge;
  SPLICE_PATH *path;
  P7_TOPHITS  *th;
  int         status;

  th = extend_graph->th;
  contains_orig = FALSE;

  /* Find best scoreing paths */
  if((status = longest_path_upstream(extend_graph)) != eslOK) goto ERROR;

  while(!contains_orig) {
    /* Find the best place to start our path */
    best_start_score = -eslINFINITY;
    start_node  = -1;
    for (i = 0; i < extend_graph->num_nodes; i++) {
      if(extend_graph->path_scores[i] > best_start_score) {
        best_start_score = extend_graph->path_scores[i];
        start_node  = i;
      }
    }
  
    /* We have run out of paths */ 
    if(start_node < 0) return NULL; 
 
    curr_node = start_node;
  
    path_len  = 1;
  
    out_edge = NULL;
    
    /* Check each set of niehgbring edges to ensure they are not "backwards" 
       (upstream splice is not downstream of downstream splice) */
    while(extend_graph->best_out_edge[curr_node] >= 0) {

      if(curr_node < extend_graph->orig_N) contains_orig = TRUE;

      next_node = extend_graph->best_out_edge[curr_node];
      
      out_edge = p7_splicegraph_GetEdge(extend_graph, curr_node, next_node);
      if(out_edge == NULL || out_edge->splice_score == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "Edge does not exist");
  
      curr_node = next_node;
      path_len++;
    } 
  
    if(curr_node < extend_graph->orig_N) contains_orig = TRUE;  

    if(!contains_orig) extend_graph->path_scores[start_node] = -eslINFINITY;
  }

  /* Once a viable path has been found build the path */
  if ((path = p7_splicepath_Create(path_len)) == NULL) goto ERROR;
    
  path->revcomp = extend_graph->revcomp;
  
  path->hit_scores[0]    = th->hit[start_node]->dcl->aliscore;
  path->edge_scores[0]   = 0.;
  path->signal_scores[0] = 0.;

  path->node_id[0] = start_node;
  
  path->split[0]   = FALSE;
  path->hits[0]    = th->hit[start_node];

  path->downstream_spliced_amino_start[0] = path->hits[0]->dcl->ihmm;
  path->downstream_spliced_nuc_start[0]   = path->hits[0]->dcl->iali;

  if(th->hit[start_node]->dcl->tr->fs) path->frameshift = TRUE; //{ printf("start_node %d\n" ,start_node); fflush(stdout); path->frameshift = TRUE; }
  
  curr_node = start_node;
  step_cnt = 1;
  while (step_cnt < path_len) {
    next_node = extend_graph->best_out_edge[curr_node];
    out_edge = p7_splicegraph_GetEdge(extend_graph, curr_node, next_node);

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

    if(th->hit[next_node]->dcl->tr->fs) path->frameshift = TRUE; //{ printf("next_node %d\n" ,next_node); fflush(stdout); path->frameshift = TRUE; }    

    curr_node = next_node;
    step_cnt++;
  }

  path->upstream_spliced_amino_end[step_cnt] = path->hits[step_cnt-1]->dcl->jhmm;
  path->upstream_spliced_nuc_end[step_cnt]   = path->hits[step_cnt-1]->dcl->jali;

  path = check_bypass_extend(orig_graph, extend_graph, path);

  return path;

  ERROR:
    if(path != NULL) p7_splicepath_Destroy(path);
    return NULL; 

}




SPLICE_PATH*
p7_splicepath_GetBestPath(SPLICE_GRAPH *graph)
{

  int          i, s;
  int          path_len;
  int          start_node;
  int          curr_node;
  int          next_node;
  int          contains_orig;
  int          step_cnt;
  float        best_start_score;
  SPLICE_EDGE *in_edge;
  SPLICE_EDGE *out_edge;
  SPLICE_PATH *path;
  P7_TOPHITS  *th;
  int         status;

  th = graph->th;
  contains_orig = FALSE;

  /* Find best scoreing paths */
  if((status = longest_path_upstream(graph)) != eslOK) goto ERROR;

  /* Find the best scoring path in graph that contains and original or split hit */
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
    
    /* We have run out of paths */ 
    if(start_node < 0) return NULL; 

    curr_node = start_node;

    path_len  = 1;

    in_edge = NULL;
    out_edge = NULL;
    
    /* Check each set of niehgbring edges to ensure they are not "backwards" 
       (upstream splice is not downstream of downstream splice) */
    while(graph->best_out_edge[curr_node] >= 0) {

      if(graph->split_orig_id[curr_node] >= 0) contains_orig = TRUE;

      next_node = graph->best_out_edge[curr_node];

      out_edge = p7_splicegraph_GetEdge(graph, curr_node, next_node);
      if(out_edge == NULL || out_edge->splice_score == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "Edge does not exist");

      if( in_edge != NULL) {
        /* Prevent exons with less than one codon */
        if((in_edge->downstream_spliced_amino_start >= out_edge->upstream_spliced_amino_end) ||
           ((graph->revcomp && out_edge->upstream_spliced_nuc_end >= in_edge->downstream_spliced_nuc_start-1) ||
           (!graph->revcomp && in_edge->downstream_spliced_nuc_start >= out_edge->upstream_spliced_nuc_end-1))) {

          /* Remove whichever edge has the lowest score */
          if(out_edge->splice_score < in_edge->splice_score) {
            out_edge->splice_score = -eslINFINITY;
            out_edge->downstream_node_id = -1;
            graph->tot_edges--;
          }
          else {
            in_edge->splice_score = -eslINFINITY;
            in_edge->downstream_node_id = -1;
            graph->tot_edges--;
          }
          path = p7_splicepath_GetBestPath(graph);
          return path;
        }
      }
      curr_node = next_node;
      in_edge = out_edge;
      path_len++;
    } 
  
    if(graph->split_orig_id[curr_node] >= 0) contains_orig = TRUE;

    if(!contains_orig) graph->path_scores[start_node] = -eslINFINITY;

  }

  /* Once a viable path has been found build the path */
  if ((path = p7_splicepath_Create(path_len)) == NULL) goto ERROR;
    
  path->revcomp = graph->revcomp;
  
  path->hit_scores[0]    = th->hit[start_node]->dcl->aliscore;
  path->edge_scores[0]   = 0.;
  path->signal_scores[0] = 0.;

  path->node_id[0] = start_node;
  path->split[0]   = FALSE;
  path->hits[0]    = th->hit[start_node];

  path->downstream_spliced_amino_start[0] = path->hits[0]->dcl->ihmm;
  path->downstream_spliced_nuc_start[0]   = path->hits[0]->dcl->iali;

  if(th->hit[start_node]->dcl->tr->fs) path->frameshift = TRUE;

  curr_node = start_node;
  step_cnt = 1;
  while (step_cnt < path_len) {
    next_node = graph->best_out_edge[curr_node];
    out_edge = p7_splicegraph_GetEdge(graph, curr_node, next_node);

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

    if(th->hit[next_node]->dcl->tr->fs) path->frameshift = TRUE;     

    curr_node = next_node;
    step_cnt++;
  }

  path->upstream_spliced_amino_end[step_cnt] = path->hits[step_cnt-1]->dcl->jhmm;
  path->upstream_spliced_nuc_end[step_cnt]   = path->hits[step_cnt-1]->dcl->jali;

  /* Find and set splits */
  if(path_len > 1 && graph->split_orig_id[path->node_id[0]] >= 0 &&
     graph->split_orig_id[path->node_id[0]] == graph->split_orig_id[path->node_id[1]])
    path->split[0]   = TRUE;

  for(s = 1; s < path_len-1; s++) {
    if(graph->split_orig_id[path->node_id[s]] < 0) continue;
    if(graph->split_orig_id[path->node_id[s]] == graph->split_orig_id[path->node_id[s-1]])
      path->split[s]   = TRUE;
    if(graph->split_orig_id[path->node_id[s]] == graph->split_orig_id[path->node_id[s+1]])
      path->split[s]   = TRUE;
  }
  if(path_len > 1 && graph->split_orig_id[path->node_id[path_len-1]] >= 0 &&
     graph->split_orig_id[path->node_id[path_len-1]] == graph->split_orig_id[path->node_id[path_len-2]])
    path->split[path_len-1]   = TRUE;

  path = check_bypass(graph, path);

  return path;

  ERROR:
    if(path != NULL) p7_splicepath_Destroy(path);
    return NULL; 

}





SPLICE_PATH*
check_bypass(SPLICE_GRAPH *graph, SPLICE_PATH *path) 
{

  int   p, h;
  int   next, step;
  int   up_path, down_path;
  int   reconnects;
  float max_bypass_score;
  float sub_path_score;
  P7_TOPHITS  *th;
  SPLICE_EDGE *tmp_edge;

  th = graph->th; 

  /* Check each edge in a path to see if it baypasses a node with downstream edges will be severed by this path. */ 
  for(p = 0; p < path->path_len-1; p ++) {
    up_path   = path->node_id[p];
    down_path = path->node_id[p+1];
             
    tmp_edge = p7_splicegraph_GetEdge(graph, up_path, down_path);
    
    if(tmp_edge->bypass_checked) continue;
    max_bypass_score = 0.;

    for(h = 0; h < graph->num_nodes; h++) {
      
      if(!graph->node_in_graph[h]) continue;
      if(graph->split_orig_id[h] == -1) continue;
      if(h == up_path || h == down_path) continue;
      
      /* Is hit h bypassed by the edge between up_path and down_path */
      if(!hit_between(th->hit[up_path]->dcl, th->hit[h]->dcl, th->hit[down_path]->dcl, graph->revcomp)) continue;
     
      /* If hit h has no downstream edges than we assume the bypass is correct */
      if(graph->best_out_edge[h] == -1) continue;

      /* If best path downstrem of hit h eventually reconnects with the current path we assume the bypass is correct */
      //reconnects = FALSE;
      //next = h;
      //while(graph->best_out_edge[next] != -1) {
      //  for(step = p+1; step < path->path_len; step++) 
      //    if(graph->best_out_edge[next] == path->node_id[step]) reconnects = TRUE;
      //  next = graph->best_out_edge[next];
      //}
      //if (reconnects) continue;

      /* If hit h is conected to both up_path and down_path then the bypass is correct */
      if(p7_splicegraph_PathExists(graph, up_path, h) && p7_splicegraph_PathExists(graph, h, down_path))  continue; 

      sub_path_score = get_sub_path_score(graph, h, down_path);
      
      max_bypass_score = ESL_MAX(max_bypass_score, sub_path_score);
    }
   
    if(max_bypass_score > 0) {
      tmp_edge->splice_score -= max_bypass_score;
      tmp_edge->bypass_checked = TRUE;
      p7_splicepath_Destroy(path);
      path =  p7_splicepath_GetBestPath(graph);    
      return path;
    }
  }

  return path;
}

SPLICE_PATH*
check_bypass_unspliced(SPLICE_GRAPH *graph, SPLICE_PATH *path) 
{

  int   p, h;
  int   next, step;
  int   up_path, down_path;
  int   reconnects;
  float max_bypass_score;
  float sub_path_score;
  P7_TOPHITS  *th;
  SPLICE_EDGE *tmp_edge;

  th = graph->th; 

  /* Check each edge in a path to see if it baypasses a node with downstream edges will be severed by this path. */ 
  for(p = 0; p < path->path_len-1; p ++) {
    up_path   = path->node_id[p];
    down_path = path->node_id[p+1];
             
    tmp_edge = p7_splicegraph_GetEdge(graph, up_path, down_path);
    
    if(tmp_edge->bypass_checked) continue;
    max_bypass_score = 0.;

    for(h = 0; h < graph->num_nodes; h++) {
      
      if(!graph->node_in_graph[h]) continue;
      if(graph->split_orig_id[h] == -1) continue;
      if(h == up_path || h == down_path) continue;
      
      /* Is hit h bypassed by the edge between up_path and down_path */
      if(!hit_between(th->hit[up_path]->dcl, th->hit[h]->dcl, th->hit[down_path]->dcl, graph->revcomp)) continue;
     
      /* If hit h has no downstream edges than we assume the bypass is correct */
      if(graph->best_out_edge[h] == -1) continue;

      /* If hit h is conected to both up_path and down_path then the bypass is correct */
      if(p7_splicegraph_PathExists(graph, up_path, h) && p7_splicegraph_PathExists(graph, h, down_path))  continue; 

      sub_path_score = get_sub_path_score(graph, h, down_path);
      
      max_bypass_score = ESL_MAX(max_bypass_score, sub_path_score);
    }
   
    if(max_bypass_score > 0) {
      tmp_edge->splice_score -= max_bypass_score;
      tmp_edge->bypass_checked = TRUE;
      p7_splicepath_Destroy(path);
      path =  p7_splicepath_GetBestPath_Unspliced(graph);    
      return path;
    }
  }

  return path;
}

SPLICE_PATH*
check_bypass_extend(SPLICE_GRAPH *orig_graph, SPLICE_GRAPH *extend_graph, SPLICE_PATH *path) 
{

  int   p, h;
  int   next, step;
  int   up, down;
  int   up_path, down_path;
  int   reconnects;
  float max_bypass_score;
  float sub_path_score;
  P7_TOPHITS  *orig_th;
  P7_TOPHITS  *extend_th;
  SPLICE_EDGE *tmp_edge;

  orig_th = orig_graph->th; 
  extend_th = extend_graph->th;

  /* Check each edge in a path to see if it baypasses a node with downstream edges will be severed by this path. */ 
  for(p = 0; p < path->path_len-1; p ++) {
    up_path   = path->node_id[p];
    down_path = path->node_id[p+1];
    up        = extend_graph->orig_hit_idx[up_path];
    down      = extend_graph->orig_hit_idx[down_path];  
 
    /* We are only concerned with the edges that connect to the original node */
    if(up == -1 && down == -1) continue;
    tmp_edge = p7_splicegraph_GetEdge(extend_graph, up_path, down_path);
    
    max_bypass_score = 0.;
    
    for(h = 0; h < orig_graph->orig_N; h++) {
      
      if(h == up || h == down) continue;
      
      /* Is hit h bypassed by the edge between up_path and down_path */
      if(!hit_between(path->hits[p]->dcl, orig_th->hit[h]->dcl, path->hits[p+1]->dcl, orig_graph->revcomp)) continue;
      
      if     (up == -1   && p7_splice_HitUpstream(orig_th->hit[h]->dcl, path->hits[p+1]->dcl, path->revcomp)) 
        max_bypass_score += ESL_MAX(0.0, orig_th->hit[h]->dcl->aliscore);
      else if(down == -1 && p7_splice_HitUpstream(path->hits[p]->dcl, orig_th->hit[h]->dcl, path->revcomp)) 
        max_bypass_score += ESL_MAX(0.0, orig_th->hit[h]->dcl->aliscore);
    }
   
    if(max_bypass_score > 0) {
      tmp_edge->splice_score -= max_bypass_score;
      tmp_edge->bypass_checked = TRUE;
      p7_splicepath_Destroy(path);
      path =  p7_splicepath_GetBestPath_Extension(orig_graph, extend_graph);    
      return path;
    }
  }

  return path;
}



int
longest_path_upstream (SPLICE_GRAPH *graph)
{

  int         i;
  int         up, down;
  int         stack_size;
  float       step_score;
  int         *visited;
  int         *stack;
  SPLICE_EDGE *edge;
  SPLICE_EDGE *tmp_edge;
  int          status;

  /* Reset path and edge data */ 
  for(i = 0; i < graph->num_nodes; i++) {
    if(graph->node_in_graph[i])
      graph->path_scores[i]   = graph->ali_scores[i];
    else
      graph->path_scores[i]   = -eslINFINITY;
    graph->best_out_edge[i] = -1;
  }

  /* Append source node downstream of all nodes with no outgoing edges*/
  if((status = p7_splicegraph_Grow(graph)) != eslOK) goto ERROR;
  graph->ali_scores[graph->num_nodes]  = 0.;
  graph->edges[graph->num_nodes] = NULL;
  for(up = 0; up < graph->num_nodes; up++) {

    if(!graph->node_in_graph[up])     continue;
    if(has_out_edge(graph, up)) continue;
    edge = p7_splicegraph_AddEdge(graph, up, graph->num_nodes);

    edge->splice_score = 0.;
    graph->num_edges[graph->num_nodes] = 0;

  }
  graph->path_scores[graph->num_nodes] = 0.;
  graph->num_nodes++;

  ESL_ALLOC(visited, sizeof(int) * graph->num_nodes);
  esl_vec_ISet(visited,   graph->num_nodes, 0);

  ESL_ALLOC(stack, sizeof(int) * graph->num_nodes);

  stack_size = 0;

  for(i = 0; i < graph->num_nodes; i++) {
    if (i < graph->num_nodes-1 && !graph->node_in_graph[i]) continue;
    if(!visited[i]) {
      topological_sort_upstream(graph, visited, stack, &stack_size, i);
    }
  }
   
   while(stack_size > 0) {
    /*pop top of stack */
    down = stack[stack_size-1];
    stack_size--;
    if(down == graph->num_nodes-1) continue;

    /* Find nodes with ougoing edge to down*/
    for(up = 0; up < graph->num_nodes-1; up++) {
      if (!graph->node_in_graph[up]) continue;
         
      tmp_edge = p7_splicegraph_GetEdge(graph, up, down);
      if(tmp_edge != NULL && tmp_edge->splice_score != -eslINFINITY) {
        step_score = graph->ali_scores[up] + tmp_edge->splice_score + graph->path_scores[down];
        if(graph->path_scores[up] <= step_score) {
          graph->path_scores[up]   = step_score;
          if(down != graph->num_nodes-1)
            graph->best_out_edge[up] = down;
        }
      }
    }
  }
  
  /* Erase source node */
  graph->num_nodes--;

  for(up = 0; up < graph->num_nodes; up++) {
    if(!graph->node_in_graph[up]) continue;
    for(i = 0; i < graph->num_edges[up]; i++) { 
      if(graph->edges[up][i].downstream_node_id == graph->num_nodes) {
        graph->edges[up][i].downstream_node_id = -1;
        i = graph->num_edges[up];
      }
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
topological_sort_upstream(SPLICE_GRAPH *graph, int *visited, int *stack, int *stack_size, int node)
{

  int i;

  visited[node] = TRUE;

  for(i = 0; i < graph->num_nodes; i++) {
    if(i < graph->num_nodes-1 && !graph->node_in_graph[i]) continue;

    if(visited[i]) continue;

    if(p7_splicegraph_EdgeExists(graph, i, node))
      topological_sort_upstream(graph, visited, stack, stack_size, i);

  }

  stack[*stack_size] = node;
  *stack_size += 1;

  return eslOK;

}

int
hit_between(P7_DOMAIN *upstream, P7_DOMAIN *middle, P7_DOMAIN *downstream, int revcomp)
{

  if (( revcomp  && upstream->jali <= downstream->iali) ||
     ((!revcomp) && upstream->jali >= downstream->iali))
    return FALSE;

  if (( revcomp  && upstream->jali <= middle->iali) ||
     ((!revcomp) && upstream->jali >= middle->iali))
    return FALSE;

  if (( revcomp  && middle->jali <= downstream->iali) ||
     ((!revcomp) && middle->jali >= downstream->iali))
    return FALSE;

  return TRUE;

}


int
has_out_edge(SPLICE_GRAPH *graph, int node_id)
{

  int i;
  
  for(i = 0; i < graph->num_edges[node_id]; i++) {
    if(graph->edges[node_id][i].downstream_node_id >= 0)
      return TRUE;
  }
  return FALSE;
}

float 
get_sub_path_score(SPLICE_GRAPH *graph, int source_node, int termination_node) 
{

  int    curr_node;
  int    prev_node;
  float  sub_path_score;
  P7_HIT *source_hit;
  P7_HIT *term_hit;
  P7_HIT *curr_hit;
  SPLICE_EDGE *tmp_edge;  

  source_hit = graph->th->hit[source_node]; 
  term_hit   = graph->th->hit[termination_node];
  sub_path_score = graph->path_scores[source_node];
  prev_node = source_node;
  curr_node = graph->best_out_edge[source_node];
  while(curr_node >= 0) {
    curr_hit = graph->th->hit[curr_node];
    if(!hit_between(source_hit->dcl, curr_hit->dcl, term_hit->dcl, graph->revcomp)) {
              
      sub_path_score -= graph->path_scores[curr_node];
      tmp_edge = p7_splicegraph_GetEdge(graph, prev_node, curr_node);
      sub_path_score -= tmp_edge->splice_score; 
      break;
    }
    prev_node = curr_node;
    curr_node = graph->best_out_edge[curr_node];
  }
  
  return sub_path_score;

}
/*****************************************************************
 * 2. Debugging tools.
 *****************************************************************/


/* Function:  p7_splicepath_Check()
 *
 * Purpose: Checks that path sequence coordinates 
 *          are all in uptream to downstream order 
 *
 */
int
p7_splicepath_Check(SPLICE_PATH *path) 
{
  int s1, s2;

  if(path->revcomp) {
    for(s1 = 0; s1 < path->path_len; s1++)
    {
      for(s2 = s1+1; s2 < path->path_len; s2++)
      {
        if(path->hits[s2]->dcl->iali > path->hits[s1]->dcl->iali)
          return FALSE;
      }
    }
  }
  else {

    for(s1 = 0; s1 < path->path_len; s1++) 
    {
      for(s2 = s1+1; s2 < path->path_len; s2++)
      {
        if(path->hits[s2]->dcl->iali < path->hits[s1]->dcl->iali)
          return FALSE;
      }
    } 
  }
 
  return TRUE;
  
}


/* Function:  p7_splicepath_Dump()
 *
 * Purpose: Dumps splice coords and score data for 
 *          each hit in <path> 
 *
 */
void
p7_splicepath_Dump(FILE *fp, SPLICE_PATH *path)
{

  int i;

  if(path == NULL) return;

  fprintf(fp, "  Path Length  %d\n", path->path_len);
  fprintf(fp, "  %4s %4s %9s %9s %10s %10s %9s %10s \n", "Step", "Node", "hmm_start", "hmm_end", "seq_start", "seq_end", "hit_score", "edge_score");
  for(i = 0; i < path->path_len; i++) {
    fprintf(fp, "  %4d %4d %9d %9d %10d %10d %9.2f %10.2f\n", i+1, path->node_id[i]+1,
      path->downstream_spliced_amino_start[i], path->upstream_spliced_amino_end[i+1], path->downstream_spliced_nuc_start[i], path->upstream_spliced_nuc_end[i+1], path->hit_scores[i], path->edge_scores[i]);
  }

  fprintf(fp, "\n");

  return;
}
