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
  ESL_ALLOC(path->node_id, sizeof(int)*path->alloc_len);

  path->ihmm = NULL;
  path->jhmm = NULL;
  ESL_ALLOC(path->ihmm,    sizeof(int)*path->alloc_len);
  ESL_ALLOC(path->jhmm,    sizeof(int)*path->alloc_len);

  path->iali = NULL;
  path->jali = NULL;
  ESL_ALLOC(path->iali,    sizeof(int64_t)*path->alloc_len);
  ESL_ALLOC(path->jali,    sizeof(int64_t)*path->alloc_len);

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

  ESL_REALLOC(path->node_id, sizeof(int)     * path->alloc_len);
  ESL_REALLOC(path->ihmm,    sizeof(int)     * path->alloc_len);
  ESL_REALLOC(path->jhmm,    sizeof(int)     * path->alloc_len);
  ESL_REALLOC(path->iali,    sizeof(int64_t) * path->alloc_len);
  ESL_REALLOC(path->jali,    sizeof(int64_t) * path->alloc_len);

  return eslOK;

  ERROR:
    p7_splicepath_Destroy(path);
    return status;

}


SPLICE_PATH*
p7_splicepath_Clone(SPLICE_PATH *path)
{
  int s;
  SPLICE_PATH *ret_path;

  ret_path = p7_splicepath_Create(path->path_len);

  ret_path->revcomp    = path->revcomp;
  ret_path->frameshift = path->frameshift;

  for(s = 0; s < path->path_len; s++) {
    ret_path->node_id[s] = path->node_id[s];
    
    ret_path->ihmm[s] = path->ihmm[s];
    ret_path->jhmm[s] = path->jhmm[s];
    ret_path->iali[s] = path->iali[s];
    ret_path->jali[s] = path->jali[s];
  }

  return ret_path;
}
/* Function:  p7_splicepath_Insert()
 * Synopsis:  Insert a new step to a  path
 *
 * Purpose:   Add new <step> idx top existing <SPLICE_PATH> 
 *
 * Returns:   <eslOK>.
 *
 */
int
p7_splicepath_Insert(SPLICE_PATH *path, int step)
{

  int s;
  
  p7_splicepath_Grow(path); 

  for(s = path->path_len; s > step; s--) {

    path->node_id[s] = path->node_id[s-1];
    path->ihmm[s]    = path->ihmm[s-1];
    path->jhmm[s]    = path->jhmm[s-1];
    path->iali[s]    = path->iali[s-1];
    path->jali[s]    = path->jali[s-1];

  }

  path->path_len++;
 
  path->node_id[step] = -1;
  path->ihmm[step]    = -1;
  path->jhmm[step]    = -1;
  path->iali[step]    = -1;
  path->jali[step]    = -1; 

  return eslOK;
 
}


/* Function:  p7_splicepath_Remove()
 * Synopsis:  Remove a step from a path.
 *
 * Purpose:   Remove the <step> inx from <SPLICE_PATH>
 *
 * Returns:   <eslOK>.
 */
int
p7_splicepath_Remove(SPLICE_PATH *path, int step)
{

  int s;

  for(s = step; s < path->path_len-1; s++) {

    path->node_id[s] = path->node_id[s+1];
    path->ihmm[s]    = path->ihmm[s+1];
    path->jhmm[s]    = path->jhmm[s+1];
    path->iali[s]    = path->iali[s+1];
    path->jali[s]    = path->jali[s+1];   
  }

  path->path_len--;
  
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
   if(path->ihmm    != NULL) free(path->ihmm);
   if(path->jhmm    != NULL) free(path->jhmm);
   if(path->iali    != NULL) free(path->iali);
   if(path->jali    != NULL) free(path->jali); 
   if(path          != NULL) free(path);

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
  
  path->node_id[0] = start_node;
  
  path->ihmm[0] = th->hit[start_node]->dcl->ihmm;
  path->iali[0] = th->hit[start_node]->dcl->iali;

  if(th->hit[start_node]->dcl->tr->fs) path->frameshift = TRUE; 
  
  curr_node = start_node;
  step_cnt = 1;
  while (step_cnt < path_len) {
    next_node = graph->best_out_edge[curr_node];
    out_edge = p7_splicegraph_GetEdge(graph, curr_node, next_node);

    path->node_id[step_cnt] = next_node;
    
    path->jhmm[step_cnt-1] = out_edge->upstream_spliced_amino_end;
    path->ihmm[step_cnt]   = out_edge->downstream_spliced_amino_start;
    path->jali[step_cnt-1] = out_edge->upstream_spliced_nuc_end;
    path->iali[step_cnt]   = out_edge->downstream_spliced_nuc_start;

    if(th->hit[next_node]->dcl->tr->fs) path->frameshift = TRUE; 

    curr_node = next_node;
    step_cnt++;
  }
  
  path->jhmm[step_cnt-1] = th->hit[curr_node]->dcl->jhmm;
  path->jali[step_cnt-1] = th->hit[curr_node]->dcl->jali;
 
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
  
  path->node_id[0] = start_node;
  
  path->ihmm[0] = th->hit[start_node]->dcl->ihmm;
  path->iali[0] = th->hit[start_node]->dcl->iali;

  if(th->hit[start_node]->dcl->tr->fs) path->frameshift = TRUE; 
  
  curr_node = start_node;
  step_cnt = 1;
  while (step_cnt < path_len) {
    next_node = extend_graph->best_out_edge[curr_node];
    out_edge = p7_splicegraph_GetEdge(extend_graph, curr_node, next_node);

    path->node_id[step_cnt] = next_node;
    
    path->jhmm[step_cnt-1] = out_edge->upstream_spliced_amino_end;
    path->ihmm[step_cnt]   = out_edge->downstream_spliced_amino_start;
    path->jali[step_cnt-1] = out_edge->upstream_spliced_nuc_end;
    path->iali[step_cnt]   = out_edge->downstream_spliced_nuc_start;

    if(th->hit[next_node]->dcl->tr->fs) path->frameshift = TRUE; 

    curr_node = next_node;
    step_cnt++;
  }

  path->jhmm[step_cnt-1] =  th->hit[curr_node]->dcl->jhmm; 
  path->jali[step_cnt-1] =  th->hit[curr_node]->dcl->jali;

  path = check_bypass_extend(orig_graph, extend_graph, path);

  return path;

  ERROR:
    if(path != NULL) p7_splicepath_Destroy(path);
    return NULL; 

}




SPLICE_PATH*
check_bypass_unspliced(SPLICE_GRAPH *graph, SPLICE_PATH *path) 
{

  int   p, h;
  int   up_path, down_path;
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
  int   up, down;
  int   up_path, down_path;
  float max_bypass_score;
  P7_TOPHITS  *orig_th;
  SPLICE_EDGE *tmp_edge;

  orig_th = orig_graph->th; 

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
      if(!hit_between(extend_graph->th->hit[path->node_id[p]]->dcl, orig_th->hit[h]->dcl, extend_graph->th->hit[path->node_id[p+1]]->dcl, orig_graph->revcomp)) continue;
      
      if     (up == -1   && p7_splice_HitUpstream(orig_th->hit[h]->dcl, extend_graph->th->hit[path->node_id[p+1]]->dcl, path->revcomp)) 
        max_bypass_score += ESL_MAX(0.0, orig_th->hit[h]->dcl->aliscore);
      else if(down == -1 && p7_splice_HitUpstream(extend_graph->th->hit[path->node_id[p]]->dcl, orig_th->hit[h]->dcl, path->revcomp)) 
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




void
p7_splicepath_Dump(FILE *fp, SPLICE_PATH *path)
{

  int i;

  if(path == NULL) return;

  fprintf(fp, "  Path Length  %d\n", path->path_len);
  fprintf(fp, "  %4s %4s %9s %9s %10s %10s \n", "Step", "Node", "hmm_start", "hmm_end", "seq_start", "seq_end");
  for(i = 0; i < path->path_len; i++) {
    fprintf(fp, "  %4d %4d %9d %9d %10ld %10ld \n", i+1, path->node_id[i]+1,
      path->ihmm[i], path->jhmm[i], path->iali[i], path->jali[i]);
  }

  fprintf(fp, "\n");

  return;
}
