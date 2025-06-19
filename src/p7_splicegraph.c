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
  graph->num_edges      = 0;
  graph->orig_N         = 0;
  graph->split_N        = 0;
  graph->recover_N      = 0;
  graph->th->N          = 0; 

  graph->reportable     = NULL;
  graph->orig_hit_idx   = NULL;
  graph->split_orig_id  = NULL;
  graph->node_in_graph   = NULL;
  
  graph->edge_in_graph = NULL;
  graph->best_out_edge = NULL;
  graph->best_in_edge  = NULL;

  graph->path_scores   = NULL;
  graph->ali_scores    = NULL;

  graph->th->hit       = NULL;
  graph->edges         = NULL;

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

  if(graph == NULL)  graph = p7_splicegraph_Create();

  graph->nalloc         = num_nodes*2;
  graph->num_nodes      = 0;

  ESL_ALLOC(graph->reportable,    sizeof(int)  * graph->nalloc);
  ESL_ALLOC(graph->orig_hit_idx,  sizeof(int)  * graph->nalloc);
  ESL_ALLOC(graph->split_orig_id, sizeof(int)  * graph->nalloc);
  ESL_ALLOC(graph->node_in_graph,  sizeof(int)  * graph->nalloc);
 
  ESL_ALLOC(graph->edge_in_graph, sizeof(int)  * graph->nalloc); 
  ESL_ALLOC(graph->best_out_edge, sizeof(int)   * graph->nalloc);
  ESL_ALLOC(graph->best_in_edge,  sizeof(int)   * graph->nalloc);

  ESL_ALLOC(graph->ali_scores,    sizeof(float) * graph->nalloc);
  ESL_ALLOC(graph->path_scores,   sizeof(float) * graph->nalloc);

  ESL_ALLOC(graph->th->hit,  sizeof(P7_HIT*)      * graph->nalloc);  
  ESL_ALLOC(graph->edges,    sizeof(SPLICE_EDGE*) * graph->nalloc);

  for(i = 0; i < graph->num_nodes; i++) 
    graph->edges[i]   = NULL;
  
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

    ESL_REALLOC(graph->node_in_graph,  sizeof(int)   * graph->nalloc);
    ESL_REALLOC(graph->reportable,    sizeof(int)   * graph->nalloc);
    ESL_REALLOC(graph->orig_hit_idx,  sizeof(int)   * graph->nalloc);    
    ESL_REALLOC(graph->split_orig_id, sizeof(int)   * graph->nalloc);

    ESL_REALLOC(graph->edge_in_graph, sizeof(int)   * graph->nalloc);
    ESL_REALLOC(graph->best_out_edge, sizeof(int)   * graph->nalloc);
    ESL_REALLOC(graph->best_in_edge,  sizeof(int)   * graph->nalloc);

    ESL_REALLOC(graph->ali_scores,    sizeof(float) * graph->nalloc);
    ESL_REALLOC(graph->path_scores,   sizeof(float) * graph->nalloc); 

    ESL_REALLOC(graph->th->hit,  sizeof(P7_HIT*)      * graph->nalloc);
    ESL_REALLOC(graph->edges,    sizeof(SPLICE_EDGE*) * graph->nalloc);
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
  SPLICE_EDGE    *tmp_edge;

  if (graph == NULL) return;

  if(graph->reportable    != NULL) free(graph->reportable);
  if(graph->orig_hit_idx  != NULL) free(graph->orig_hit_idx);
  if(graph->split_orig_id != NULL) free(graph->split_orig_id);
  if(graph->node_in_graph  != NULL) free(graph->node_in_graph);

  if(graph->path_scores   != NULL) free(graph->path_scores);
  if(graph->ali_scores    != NULL) free(graph->ali_scores);

  if(graph->edge_in_graph != NULL) free(graph->edge_in_graph);
  if(graph->best_out_edge != NULL) free(graph->best_out_edge);
  if(graph->best_in_edge  != NULL) free(graph->best_in_edge);

 
  for (i = graph->orig_N; i < graph->th->N; i++) {
    p7_trace_splice_Destroy(graph->th->hit[i]->dcl->tr);
    free(graph->th->hit[i]->dcl->scores_per_pos);
    p7_hit_Destroy(graph->th->hit[i]);
  }

  if (graph->th->hit != NULL) free(graph->th->hit);
  if (graph->th      != NULL) free(graph->th);

  for(i = 0; i < graph->num_nodes; i++) {
    while(graph->edges[i] != NULL) {
      tmp_edge = graph->edges[i];
      graph->edges[i] = tmp_edge->next;
      free(tmp_edge);
    }
  }

  if(graph->edges   != NULL) free(graph->edges);

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

  int h;
  P7_TOPHITS  *th;
  int status;

  th = graph->th;

  /* Make sure we have room for our new node */
  if((status = p7_splicegraph_Grow(graph)) != eslOK) goto ERROR;

  th->hit[graph->num_nodes] = hit;
  th->N++;

  graph->node_in_graph[graph->num_nodes]  = TRUE;

  graph->ali_scores[graph->num_nodes]  = hit->dcl->aliscore; //This contiains the summed alisc
  graph->path_scores[graph->num_nodes] = -eslINFINITY;
 
  graph->best_out_edge[graph->num_nodes] = -1;
  graph->best_in_edge[graph->num_nodes]  = -1;

  graph->edges[graph->num_nodes]   = NULL;  

  graph->split_orig_id[graph->num_nodes]   = -1;
 
  graph->num_nodes++;

  return eslOK; 

  ERROR:
    return status;
}


/* Function:  p7_splicegraph_AddEdge()
 * Synopsis:  append an edge to a splice graph
 *
 * Purpose:   Takes a splice edge <edge> that has been 
 *            allocated and aissigned an unpstream and 
 *            downstream node, and adds it to a splice 
 *            graph <graph> 
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslFAIL> if assigned nodes either do not 
 *            exist or do not have correct upstream/
 *            downstream orientation.
 */
int
p7_splicegraph_AddEdge(SPLICE_GRAPH *graph, SPLICE_EDGE *edge)
{

  int up_node;
  int down_node;
  int status;

  SPLICE_EDGE *tmp_edge;

  up_node   = edge->upstream_node_id;
  down_node = edge->downstream_node_id;

  if(up_node >= graph->num_nodes || down_node >= graph->num_nodes) goto ERROR; 

  /* Find the end of the up nodes list of edges and append the new edge */
  tmp_edge = graph->edges[up_node];
  if(tmp_edge == NULL)
    graph->edges[up_node] = edge;
  else {
    while(tmp_edge->next != NULL)
      tmp_edge = tmp_edge->next;

    tmp_edge->next = edge;
    edge->prev = tmp_edge;
  }

  graph->num_edges++;
  
  return eslOK;

  ERROR:
    ESL_XEXCEPTION(eslFAIL, "Edge appended with impossible node assignment.");

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

  SPLICE_EDGE *tmp_edge;

  if(up_node >= graph->num_nodes || down_node >= graph->num_nodes) return FALSE;
 
  tmp_edge = graph->edges[up_node];

  while(tmp_edge != NULL && tmp_edge->downstream_node_id != down_node)
    tmp_edge = tmp_edge->next;

  if     (tmp_edge == NULL)                       return FALSE;
  else if(tmp_edge->splice_score == -eslINFINITY) return FALSE;
  else                                            return TRUE;

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

  SPLICE_EDGE *edge;

  edge = graph->edges[up_node];

  while(edge != NULL && edge->downstream_node_id != down_node)
    edge = edge->next;

  return edge;

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

    if ( i < graph->orig_N)
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


  int          i;
  int          nuc_end, nuc_start;
  SPLICE_EDGE *tmp_edge;

  if (graph == NULL) { fprintf(fp, " [ graph is NULL ]\n"); return; }

  fprintf(fp, "\n Edge Data  \n\n");
  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!graph->node_in_graph[i]) continue;
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
  fprintf(fp, " Number of Edges  : %d\n", graph->num_edges); 

  if(graph->num_edges > 0) {
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


      tmp_edge = graph->edges[i];
      while(tmp_edge != NULL) {
        edge_scores[tmp_edge->downstream_node_id] = tmp_edge->splice_score;
        tmp_edge = tmp_edge->next;
      }

      for(j = 0; j < graph->num_nodes; j++ ) { 
        if(graph->node_in_graph[j])
          fprintf(fp, "%8.2f", edge_scores[j]);
      }
      fprintf(fp, "\n");
    }    
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

  fprintf(fp, "\n Upstream Path Scores  \n");
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

  fprintf(fp, "\n Downstream Path Scores  \n");
  for(i = 0; i < graph->num_nodes; i++ ) {
    if(!graph->node_in_graph[i]) continue;
    fprintf(fp, "%8d", i+1);
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


