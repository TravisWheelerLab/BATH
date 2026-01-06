/* SPLICE_GRAPH: graph structure with nodes (hits) and edges (splice positions and scores)
 *
 * Contents:
 *    1. The SPLICE_GRAPH object.
 *    2. Access routines.
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


static int path_finder (SPLICE_GRAPH *graph, int upstream_node, int downstream_node, int *visited); 


/*****************************************************************
 * 1. The SPLICE_GRAPH structure.
 *****************************************************************/

/* Function:  p7_splicegraph_Create()
 * Synopsis:  Allocates a splice graph with no nodes or edges.
 *
 * Purpose:   Allocates a new <SPLICE_GRAPH>. Does NOT allocate 
 *            any memory for nodes or edges.  Sets reverse 
 *            complement startus to <revcomp>
 *
 * Returns:   a pointer to the new <SPLICE_GRAPH> structure 
 *            on success.
 *
 * Throws:    <NULL> on allocation error.
 */
SPLICE_GRAPH*
p7_splicegraph_Create()
{
  int status;
  SPLICE_GRAPH *graph;

  ESL_ALLOC(graph, sizeof(SPLICE_GRAPH));

  graph->th = NULL;
  ESL_ALLOC(graph->th, sizeof(P7_TOPHITS));
  graph->th->is_sorted_by_seqidx = FALSE;

  graph->nalloc         = 0;
  graph->num_nodes      = 0;
  graph->anchor_N         = 0;
  graph->th->N          = 0; 

  graph->reportable     = NULL;
  graph->orig_hit_idx   = NULL;
  graph->split_orig_id  = NULL;
  graph->node_in_graph  = NULL;
  graph->tmp_node       = NULL;
  
  graph->best_out_edge = NULL;

  graph->path_scores   = NULL;
  graph->ali_scores    = NULL;

  graph->th->hit       = NULL;
  graph->edges         = NULL;
  graph->num_edges     = NULL;
  return graph;

ERROR:
    p7_splicegraph_Destroy(graph);
    return NULL;

}

/* Function:  p7_splicegraph_CreateNodes()
 * Synopsis:  Allocates memory for <num_nodes> nodes in a 
 *            splice graph.
 *
 * Purpose:   Allocates memory for <num_nodes> nodes and
 *            all accompanying data. If <SPLICE_GRAPH> has 
 *            not previously been allocated, calls 
 *            p7_splicegraph_Create().
 *  
 * Returns:  <eslOK> on success. 
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_splicegraph_CreateNodes(SPLICE_GRAPH *graph, int num_nodes)
{

  int i;
  int status;

  graph->nalloc         = num_nodes*2;
  graph->num_nodes      = 0;

  ESL_ALLOC(graph->reportable,    sizeof(int)  * graph->nalloc);
  ESL_ALLOC(graph->orig_hit_idx,  sizeof(int)  * graph->nalloc);
  ESL_ALLOC(graph->split_orig_id, sizeof(int)  * graph->nalloc);
  ESL_ALLOC(graph->node_in_graph, sizeof(int)  * graph->nalloc);
  ESL_ALLOC(graph->tmp_node,      sizeof(int)  * graph->nalloc);
 
  ESL_ALLOC(graph->best_out_edge, sizeof(int)   * graph->nalloc);

  ESL_ALLOC(graph->ali_scores,    sizeof(float) * graph->nalloc);
  ESL_ALLOC(graph->path_scores,   sizeof(float) * graph->nalloc);

  ESL_ALLOC(graph->th->hit,   sizeof(P7_HIT*)      * graph->nalloc);  
  ESL_ALLOC(graph->edges,     sizeof(SPLICE_EDGE*) * graph->nalloc);
  ESL_ALLOC(graph->edge_mem,  sizeof(int)          * graph->nalloc);
  ESL_ALLOC(graph->num_edges, sizeof(int)          * graph->nalloc);

  for(i = 0; i < graph->num_nodes; i++) {
    ESL_ALLOC(graph->edges[i], sizeof(SPLICE_EDGE) * EDGE_ALLOC);
    graph->edge_mem[i]  = EDGE_ALLOC;
    graph->num_edges[i] = 0;
  }
  
  return eslOK;

  ERROR:
    p7_splicegraph_Destroy(graph);
    return status;
}


/* Function:  p7_splicegraph_Grow()
 * Synopsis:  Reallocates a larger splice graph, if needed.
 *
 * Purpose:   If <SPLICE_GRAPH> cannot hold another node, 
 *            doubles the internal allocation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *   
 */
int
p7_splicegraph_Grow(SPLICE_GRAPH *graph)
{
  int status;

  if(graph->num_nodes == graph->nalloc) {

    graph->nalloc *= 2;

    ESL_REALLOC(graph->node_in_graph, sizeof(int)   * graph->nalloc);
    ESL_REALLOC(graph->tmp_node,      sizeof(int)   * graph->nalloc);
    ESL_REALLOC(graph->reportable,    sizeof(int)   * graph->nalloc);
    ESL_REALLOC(graph->orig_hit_idx,  sizeof(int)   * graph->nalloc);    
    ESL_REALLOC(graph->split_orig_id, sizeof(int)   * graph->nalloc);
    ESL_REALLOC(graph->best_out_edge, sizeof(int)   * graph->nalloc);

    ESL_REALLOC(graph->ali_scores,    sizeof(float) * graph->nalloc);
    ESL_REALLOC(graph->path_scores,   sizeof(float) * graph->nalloc); 

    ESL_REALLOC(graph->th->hit,  sizeof(P7_HIT*)      * graph->nalloc);
    ESL_REALLOC(graph->edges,    sizeof(SPLICE_EDGE*) * graph->nalloc);

    ESL_REALLOC(graph->edge_mem,  sizeof(int)         * graph->nalloc);
    ESL_REALLOC(graph->num_edges, sizeof(int)         * graph->nalloc);   
  }
     
  return eslOK;

  ERROR:
    p7_splicegraph_Destroy(graph);
    return status;

}


/* Function:  p7_splicegraph_Destroy()
 *
 * Purpose:  Frees a <SPLICE_GRAPH>
 */
void 
p7_splicegraph_Destroy(SPLICE_GRAPH *graph)
{

  int i;

  if (graph == NULL) return;

  if(graph->reportable    != NULL) free(graph->reportable);
  if(graph->orig_hit_idx  != NULL) free(graph->orig_hit_idx);
  if(graph->split_orig_id != NULL) free(graph->split_orig_id);
  if(graph->node_in_graph != NULL) free(graph->node_in_graph);
 
  if(graph->path_scores   != NULL) free(graph->path_scores);
  if(graph->ali_scores    != NULL) free(graph->ali_scores);

  if(graph->best_out_edge != NULL) free(graph->best_out_edge);

  for(i = 0; i < graph->num_nodes; i++) {
    if(graph->tmp_node[i]) {
      p7_trace_fs_Destroy(graph->th->hit[i]->dcl->tr);
      p7_hit_Destroy(graph->th->hit[i]);
    }   
  }

  if(graph->tmp_node != NULL) free(graph->tmp_node);

  if (graph->th->hit != NULL) free(graph->th->hit);
  if (graph->th      != NULL) free(graph->th);

  for(i = 0; i < graph->num_nodes; i++) 
   if(graph->edges[i] != NULL) free(graph->edges[i]);

  if(graph->edges     != NULL) free(graph->edges);
  if(graph->edge_mem  != NULL) free(graph->edge_mem);
  if(graph->num_edges != NULL) free(graph->num_edges);

  graph->seqname = NULL;

  free(graph);
  graph = NULL;

  return;
}


/* Function:  p7_splicegraph_AddNode()
 * Synopsis:  create new node and append to splice graph
 *
 * Purpose:   Ask the splice graph object <graph> to do any 
 *            necessary internal allocation to add a new node. 
 *            Use hit dats to determine the upstream/downstream 
 *            postion of the new node relative to all existing 
 *            nodes. 
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_splicegraph_AddNode(SPLICE_GRAPH *graph, P7_HIT *hit)
{

  P7_TOPHITS  *th;
  int status;

  th = graph->th;

  /* Make sure we have room for our new node */
  if((status = p7_splicegraph_Grow(graph)) != eslOK) goto ERROR;

  th->hit[graph->num_nodes] = hit;
  th->N++;

  graph->node_in_graph[graph->num_nodes] = TRUE;
  graph->tmp_node[graph->num_nodes]      = FALSE;

  graph->ali_scores[graph->num_nodes]  = hit->dcl->aliscore;
  graph->path_scores[graph->num_nodes] = -eslINFINITY;
 
  graph->best_out_edge[graph->num_nodes] = -1;
  graph->split_orig_id[graph->num_nodes] = -1;

  ESL_ALLOC(graph->edges[graph->num_nodes], sizeof(SPLICE_EDGE) * EDGE_ALLOC);
  graph->edge_mem[graph->num_nodes]  = EDGE_ALLOC;
  graph->num_edges[graph->num_nodes] = 0;

 
  graph->num_nodes++;

  return eslOK; 

  ERROR:
    return status;
}


/* Function:  p7_splicegraph_AddEdge()
 * Synopsis:  get the next edge from up_node
 *
 * Purpose:   Retrieve a potiner to the next available 
 *            SPLICE_EDGE from up_node and points it 
 *            to down_node. Allocate more room for 
 *            edges from up_node, if needed. 
 *
 * Returns:   <SPLICE_EDGE*> on success.
 *
 */
SPLICE_EDGE*
p7_splicegraph_AddEdge(SPLICE_GRAPH *graph, int up_node, int down_node)
{

  SPLICE_EDGE *ret_edge;
  int status;

  if(graph->num_edges[up_node] == graph->edge_mem[up_node]) {
    graph->edge_mem[up_node] *= 2;
    ESL_REALLOC(graph->edges[up_node], sizeof(SPLICE_EDGE) * graph->edge_mem[up_node]);
  }
  
  ret_edge = &(graph->edges[up_node][graph->num_edges[up_node]]);
  graph->num_edges[up_node]++; 

  ret_edge->i_start = -1;
  ret_edge->k_start = -1; 

  ret_edge->i_end = -1;
  ret_edge->k_end = -1;

  ret_edge->next_i_start = -1;
  ret_edge->next_k_start = -1;

  ret_edge->upstream_node_id   = up_node;
  ret_edge->downstream_node_id = down_node; 

  ret_edge->jump_edge = FALSE;

  ret_edge->splice_score = 0.;

  return ret_edge;

  ERROR:
    return NULL;

}



/*****************************************************************
 * 2. Access routines
 *****************************************************************/


/* Function:  p7_splicegraph_EdgeExists()
 *
 * Purpose:   Determine if <up_node> has a splice edge to 
 *            <down_node>
 *
 * Returns:   TRUE if <up_node> has an edge to  <down_node>,
 *            FALSE otherwise. 
 *            
 */
int
p7_splicegraph_EdgeExists(SPLICE_GRAPH* graph, int up_node, int down_node) 
{

  int i;
  
  for(i = 0; i < graph->num_edges[up_node]; i++) {
    if(graph->edges[up_node][i].downstream_node_id == down_node) {  
       if(graph->edges[up_node][i].splice_score != -eslINFINITY) return TRUE;
       else return FALSE;
    }
  }

  return FALSE;

}

/* Function:  p7_splicegraph_GetEdge()
 *
 * Purpose:   Fetch SPLICE_EDGE from <up_node> to <down_node> 
 *
 * Returns:   SPLICE_EDGE* if found, NULL otherwise.
 *
 */
SPLICE_EDGE*
p7_splicegraph_GetEdge(SPLICE_GRAPH* graph, int up_node, int down_node) 
{

  int i;
  SPLICE_EDGE *ret_edge;

  for(i = 0; i < graph->num_edges[up_node]; i++) {
    if(graph->edges[up_node][i].downstream_node_id == down_node) {
      ret_edge = &(graph->edges[up_node][i]);
      return ret_edge;
    } 
  }
  
  return NULL;

}


int
p7_splicegraph_AliScoreEdge(SPLICE_EDGE *edge, const P7_PROFILE *gm, const P7_DOMAIN *upstream_dom, const P7_DOMAIN *downstream_dom)
{

  int p, s, z;
  int overlap_start;
  int overlap_end;
  int overlap_len;
  int z1, z2;
  float min_lost_sc;
  float curr_lost_sc;
  float upstream_lost;
  float downstream_lost;
  float *spp;
  float *upstream_suffix_sum;
  float *downstream_prefix_sum;
  P7_TRACE    *tr;
  int status;

  upstream_suffix_sum   = NULL;
  downstream_prefix_sum = NULL;

  /* return if no there is no hmm overlap */
  if(downstream_dom->ihmm > upstream_dom->jhmm)  return eslOK;

  overlap_start = ESL_MAX(upstream_dom->ihmm, downstream_dom->ihmm);
  overlap_end   = ESL_MIN(upstream_dom->jhmm, downstream_dom->jhmm);
  if(overlap_start == upstream_dom->ihmm)   overlap_start++;
  if(overlap_end   == downstream_dom->jhmm) overlap_end--;

  overlap_len   = overlap_end - overlap_start + 1;

  if(overlap_len < 1) {
    edge->splice_score = -eslINFINITY;
    return eslOK;
  }

  ESL_ALLOC(upstream_suffix_sum,   sizeof(float) * overlap_len);
  ESL_ALLOC(downstream_prefix_sum, sizeof(float) * overlap_len);

  esl_vec_FSet(upstream_suffix_sum,   overlap_len, 0.0);
  esl_vec_FSet(downstream_prefix_sum, overlap_len, 0.0);

  /* Fill the upstream array */
  tr  = upstream_dom->tr;
  spp = upstream_dom->scores_per_pos;

  z1 = tr->tfrom[0];
  z2 = tr->tto[0];
  while(tr->st[z2] != p7T_M) z2--;

  /* Move to the end of the overlap in the trace and the scores_per_pos array */
  p = z2 - z1 - 2;
  z = z2;
  while(tr->k[z] > overlap_end) { z--; p--;}

  s = overlap_len - 1;
  while(tr->k[z]  >= overlap_start) {
    upstream_suffix_sum[s] += spp[p];
    z--;
    p--;
    if(tr->k[z] < tr->k[z+1]) s--;
  }

  for(s = overlap_len - 2; s >= 0; s--)
    upstream_suffix_sum[s] += upstream_suffix_sum[s+1];

  upstream_lost = 0.;
  /*Calculate any lost score from the end of the upstream hit */
  if(upstream_dom->jhmm > overlap_end) {
    z = z1+1;
    p = 0;
    while(tr->k[z] <= overlap_end) { z++; p++; }
    while(z < z2) {
      upstream_lost += spp[p];
      z++;
      p++;
    }
  }

  /* Fill the downstream array */
  tr  = downstream_dom->tr;
  spp = downstream_dom->scores_per_pos;

  z1 = tr->tfrom[0];
  z2 = tr->tto[0];

 /* Move to the start of the overlap in the trace and the scores_per_pos array */
  p = 0;
  z = z1+1;
  
  while(tr->k[z] < overlap_start) { z++; p++;}

  s = 0;
  while(tr->k[z] <= overlap_end && z < z2) {
    downstream_prefix_sum[s] += spp[p];
    z++;
    p++;
    if(tr->k[z] > tr->k[z-1]) s++;
  }

  for(s = 1; s < overlap_len; s++)
    downstream_prefix_sum[s] +=  downstream_prefix_sum[s-1];

  downstream_lost = 0.;
  /*Calculate any lost score from the start of the downstream hit */
  if(downstream_dom->ihmm < overlap_start) {
    z = z2-1;
    p = z2 - z1 - 2;
    while(tr->k[z] >= overlap_start) { z--; p--; }
    while(z > z1) {
      downstream_lost += spp[p];
      z--;
      p--;
    }
  }


  /* Find the minimiom score loss to eliminate the overlap*/
  min_lost_sc = upstream_suffix_sum[0];
  for(s = 1; s < overlap_len; s++) {
    curr_lost_sc = upstream_suffix_sum[s] + downstream_prefix_sum[s-1];
    min_lost_sc = ESL_MIN(min_lost_sc, curr_lost_sc);
  }
  min_lost_sc = ESL_MIN(min_lost_sc, downstream_prefix_sum[overlap_len-1]);

  edge->splice_score -= (min_lost_sc + upstream_lost + downstream_lost);

  if(upstream_suffix_sum   != NULL) free(upstream_suffix_sum);
  if(downstream_prefix_sum != NULL) free(downstream_prefix_sum);

  return eslOK;

  ERROR:
    if(upstream_suffix_sum   != NULL) free(upstream_suffix_sum);
    if(downstream_prefix_sum != NULL) free(downstream_prefix_sum);
    return status;

}


/* Function:  p7_splicegraph_PathExists()
 *
 * Purpose:   Determine if <down_node> is reachaeble from 
 *            <up_node> by splice edges
 *
 * Returns:   TRUE if <down_node> is reachable from <up_node>,
 *            FALSE otherwise.
 *
 */
int
p7_splicegraph_PathExists (SPLICE_GRAPH *graph, int up_node, int down_node)
{

  int   i;
  int   current;
  int  *visited;
  int status;

  if(p7_splicegraph_EdgeExists(graph, up_node, down_node)) return TRUE;

  visited = NULL;
  ESL_ALLOC(visited, graph->num_nodes * sizeof(int));
  esl_vec_ISet(visited,   graph->num_nodes, 0);

  current = up_node;
  visited[current] = 1;
  
  for (i = 0; i < graph->num_nodes; i++) {
    if(!visited[i]) {
      if(p7_splicegraph_EdgeExists(graph, current, i)) {
        if(path_finder(graph, i, down_node, visited)) {
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

/* Function:  path_finder()
 *
 * Purpose:   helper function for p7_splicegraph_PathExists() 
 */
int 
path_finder (SPLICE_GRAPH *graph, int upstream_node, int downstream_node, int *visited) {

  int   i;
  int   current;

  if(p7_splicegraph_EdgeExists(graph, upstream_node, downstream_node)) return TRUE;

  current = upstream_node;

  for (i = 0; i < graph->num_nodes; i++) {
    if(!visited[i]) {
      if(p7_splicegraph_EdgeExists(graph, current, i)) {
        visited[i] = 1;
        if(path_finder(graph, i, downstream_node, visited)) return TRUE;
      }
    }
  }
  
  return FALSE;
}


/*****************************************************************
 * 3. Debugging tools.
 *****************************************************************/


/* Function:  p7_splicegraph_DumpHits()
 *
 * Purpose: Dumps alignment coords and score of hits in 
 *         <graph> 
 */
void
p7_splicegraph_DumpHits(FILE *fp, SPLICE_GRAPH *graph)
{

  int i;
  int hit_cnt;
  P7_HIT *hit = NULL;

  if (graph == NULL) { fprintf(fp, " [ target range is NULL ]\n"); return; }

  hit_cnt = 0;
  for (i = 0; i < graph->th->N; i++) {
    if(graph->node_in_graph[i]) {
      hit_cnt++;
    }
  }

  fprintf(fp, " Graph Sequence Name      : %s\n", graph->seqname);
  fprintf(fp, " Graph Reverse Complement : %s\n", (graph->revcomp ? "YES" : "NO") );
  fprintf(fp, " Graph Total Hit Count    : %d\n", hit_cnt);
  fprintf(fp, "\n");

  fprintf(fp, "   %4s %6s %9s %12s %12s %10s %10s %10s %10s\n",
          "Hit", "hmm_from", "hmm_to", "seq_from", "seq_to", "score", "p-value", "reportable", "split_orig");

  for(i = 0; i < graph->th->N; i++) {

    if(!graph->node_in_graph[i])
      continue;

    hit = graph->th->hit[i];

    if ( i < graph->anchor_N)
      fprintf(fp, "   %4d %6d %9d %12" PRId64 " %12" PRId64 " %10.2f %10f %10s\n",
              i+1, hit->dcl->ihmm, hit->dcl->jhmm, hit->dcl->iali, hit->dcl->jali, hit->dcl->aliscore, exp(hit->sum_lnP), (graph->reportable[i] ? "YES" : "NO"));
    else if (graph->split_orig_id[i] >= 0)
      fprintf(fp, "   %4d %6d %9d %12" PRId64 " %12" PRId64 " %10.2f %10s %10s %10d\n",
              i+1, hit->dcl->ihmm, hit->dcl->jhmm, hit->dcl->iali, hit->dcl->jali, hit->dcl->aliscore, " ", (graph->reportable[graph->split_orig_id[i]] ? "YES" : "NO"), graph->split_orig_id[i]+1);
    else
      fprintf(fp, "   %4d %6d %9d %12" PRId64 " %12" PRId64 " %10.2f\n",
              i+1, hit->dcl->ihmm, hit->dcl->jhmm, hit->dcl->iali, hit->dcl->jali, hit->dcl->aliscore);
  }

  fprintf(fp, "\n");

  return;
}


/* Function:  p7_splicegraph_DumpEdges()
 *
 * Purpose: Dumps splice cooridnates of esch edge.
 */
void
p7_splicegraph_DumpEdges(FILE *fp, SPLICE_GRAPH *graph) 
{


  int          i,j;
  int          nuc_end, nuc_start;
  SPLICE_EDGE *tmp_edge;

  if (graph == NULL) { fprintf(fp, " [ graph is NULL ]\n"); return; }

  fprintf(fp, "\n Edge Data  \n\n");
  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!graph->node_in_graph[i]) continue;
    for(j = 0; j < graph->num_edges[i]; j++) {
      tmp_edge = &(graph->edges[i][j]);
    
      nuc_end =   tmp_edge->upstream_nuc_end;
      nuc_start = tmp_edge->downstream_nuc_start;
      fprintf(fp, "    Edge from Upstream Node %d to Downstream Node %d\n", tmp_edge->upstream_node_id+1, tmp_edge->downstream_node_id+1);
      fprintf(fp, "                                   %s   %s\n", "Amino", "Nucleotide");
      fprintf(fp, "      Upsteam Node End Coords:     %5d  %10d\n", tmp_edge->upstream_amino_end, nuc_end);
      fprintf(fp, "      Downsteam Node Start Coords: %5d  %10d\n", tmp_edge->downstream_amino_start, nuc_start);
      fprintf(fp, "\n");
   
    }
  }
  fprintf(fp, "\n");
 

  return;

}



/* Function:  p7_splicegraph_DumpGraph()
 *
 * Purpose: Dumps matrix of graph scores. If no edge exits 
 *          prints "-inf". Also Dumps ali and path scores 
 *          and best outgoing edge data. 
 */
void
p7_splicegraph_DumpGraph(FILE *fp, SPLICE_GRAPH *graph) 
{


  int          i,j;
  int          num_nodes;
  int          in_range_cnt;
  float        edge_scores[graph->num_nodes];
  SPLICE_EDGE *tmp_edge;

  if (graph == NULL) { fprintf(fp, " [ graph is NULL ]\n"); return; }
  num_nodes = esl_vec_ISum(graph->node_in_graph, graph->num_nodes);  

  fprintf(fp, " SPLICE_GRAPH\n");
  fprintf(fp, " Number of Nodes  : %d\n", num_nodes); 


  fprintf(fp, "\n Splice Score matrix \n");
  fprintf(fp, "     ");
  for(i = 0; i < num_nodes/2; i++ )
    fprintf(fp, "%8c", ' ');
  fprintf(fp, "  DOWN  \n");

  fprintf(fp, "     ");
  for(i = 0; i < graph->num_nodes; i++ ) {
    if(graph->node_in_graph[i])
      fprintf(fp, "%8d", i+1);
  }
  fprintf(fp, "\n");

  in_range_cnt = 0;
  for(i = 0; i < graph->num_nodes; i++ ) {

    if(!graph->node_in_graph[i]) continue;
    if (in_range_cnt == num_nodes/2)
      fprintf(fp, "UP");
    else
      fprintf(fp, "  ");
    in_range_cnt++;
    for(j = 0; j < graph->num_nodes; j++ ) 
      edge_scores[j] = -eslINFINITY;

    fprintf(fp, "%7d", i+1);

    for(j = 0; j < graph->num_edges[i]; j++) {
      tmp_edge = &(graph->edges[i][j]);
      if(tmp_edge->downstream_node_id >= 0)
      edge_scores[tmp_edge->downstream_node_id] = tmp_edge->splice_score;
    }

    for(j = 0; j < graph->num_nodes; j++ ) { 
      if(graph->node_in_graph[j])
        fprintf(fp, "%8.2f", edge_scores[j]);
    }
    fprintf(fp, "\n");
  }    
 

  fprintf(fp, "\n Alignment Scores  \n");
  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!graph->node_in_graph[i]) continue;
    fprintf(fp, "%8d", i+1);
  }
  fprintf(fp, "\n");

  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!graph->node_in_graph[i]) continue; 
    fprintf(fp, "%8.2f", graph->ali_scores[i]);
  }
  fprintf(fp, "\n");

  fprintf(fp, "\n Path Scores  \n");
  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!graph->node_in_graph[i]) continue;
    fprintf(fp, "%8d", i+1);
  }
  fprintf(fp, "\n");

  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!graph->node_in_graph[i]) continue;
    fprintf(fp, "%8.2f", graph->path_scores[i]);
  }
  fprintf(fp, "\n");

  fprintf(fp, "\n Best Edge  \n");
  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!graph->node_in_graph[i]) continue;
    fprintf(fp, "%8d", i+1);
  }
  fprintf(fp, "\n");

  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!graph->node_in_graph[i]) continue;
    fprintf(fp, "%8d", graph->best_out_edge[i]+1);
  }
  fprintf(fp, "\n");

  return;

}


