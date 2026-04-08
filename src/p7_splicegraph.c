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


static int  path_finder      (SPLICE_GRAPH *graph, int upstream_node, int downstream_node, int *visited);
static int* find_components  (SPLICE_GRAPH *graph, int *K_out);


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

  graph->orig_hit_idx   = NULL;
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

  ESL_ALLOC(graph->orig_hit_idx,  sizeof(int)  * graph->nalloc);
  ESL_ALLOC(graph->node_in_graph, sizeof(int)  * graph->nalloc);
  ESL_ALLOC(graph->tmp_node,      sizeof(int)  * graph->nalloc);
 
  ESL_ALLOC(graph->best_out_edge, sizeof(int)   * graph->nalloc);

  ESL_ALLOC(graph->ali_scores,    sizeof(float) * graph->nalloc);
  ESL_ALLOC(graph->path_scores,   sizeof(float) * graph->nalloc);

  ESL_ALLOC(graph->th->hit,   sizeof(P7_HIT*)      * graph->nalloc);  
  ESL_ALLOC(graph->edges,     sizeof(SPLICE_EDGE*) * graph->nalloc);
  ESL_ALLOC(graph->edge_mem,  sizeof(int)          * graph->nalloc);
  ESL_ALLOC(graph->num_edges, sizeof(int)          * graph->nalloc);

  ESL_ALLOC(graph->in_nodes,     sizeof(int*) * graph->nalloc);
  ESL_ALLOC(graph->in_node_mem,  sizeof(int)  * graph->nalloc);
  ESL_ALLOC(graph->num_in_nodes, sizeof(int)  * graph->nalloc);

  for(i = 0; i < graph->nalloc; i++) {
    graph->edges[i]       = NULL;
    graph->edge_mem[i]    = 0;
    graph->num_edges[i]   = 0;
    graph->in_nodes[i]    = NULL;
    graph->in_node_mem[i] = 0;
    graph->num_in_nodes[i]= 0;
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
  int i;
  int old_alloc;
  int status;

  if(graph->num_nodes == graph->nalloc) {

    old_alloc = graph->nalloc;
    graph->nalloc *= 2;

    ESL_REALLOC(graph->node_in_graph, sizeof(int)   * graph->nalloc);
    ESL_REALLOC(graph->tmp_node,      sizeof(int)   * graph->nalloc);
    ESL_REALLOC(graph->orig_hit_idx,  sizeof(int)   * graph->nalloc);    
    ESL_REALLOC(graph->best_out_edge, sizeof(int)   * graph->nalloc);

    ESL_REALLOC(graph->ali_scores,    sizeof(float) * graph->nalloc);
    ESL_REALLOC(graph->path_scores,   sizeof(float) * graph->nalloc); 

    ESL_REALLOC(graph->th->hit,  sizeof(P7_HIT*)      * graph->nalloc);
    ESL_REALLOC(graph->edges,    sizeof(SPLICE_EDGE*) * graph->nalloc);

    ESL_REALLOC(graph->edge_mem,  sizeof(int)  * graph->nalloc);
    ESL_REALLOC(graph->num_edges, sizeof(int)  * graph->nalloc);

    ESL_REALLOC(graph->in_nodes,     sizeof(int*) * graph->nalloc);
    ESL_REALLOC(graph->in_node_mem,  sizeof(int)  * graph->nalloc);
    ESL_REALLOC(graph->num_in_nodes, sizeof(int)  * graph->nalloc);

    for(i = old_alloc; i < graph->nalloc; i++) {
      graph->edges[i]        = NULL;
      graph->edge_mem[i]     = 0;
      graph->num_edges[i]    = 0;
      graph->in_nodes[i]     = NULL;
      graph->in_node_mem[i]  = 0;
      graph->num_in_nodes[i] = 0;
    }

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

  if(graph->orig_hit_idx  != NULL) free(graph->orig_hit_idx);
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

  for(i = 0; i < graph->nalloc; i++)
    if(graph->edges[i] != NULL) free(graph->edges[i]);

  if(graph->edges     != NULL) free(graph->edges);
  if(graph->edge_mem  != NULL) free(graph->edge_mem);
  if(graph->num_edges != NULL) free(graph->num_edges);

  for(i = 0; i < graph->nalloc; i++)
    if(graph->in_nodes[i] != NULL) free(graph->in_nodes[i]);

  if(graph->in_nodes     != NULL) free(graph->in_nodes);
  if(graph->in_node_mem  != NULL) free(graph->in_node_mem);
  if(graph->num_in_nodes != NULL) free(graph->num_in_nodes);

  graph->seqname = NULL;

  free(graph);
  graph = NULL;

  return;
}


/* Function: find_components()  [static]
 *
 * Purpose:  Find weakly-connected components in <graph> using BFS
 *           over both outgoing and incoming edges.  Only nodes where
 *           node_in_graph[i] is TRUE are visited.
 *
 *           Returns an allocated int array of length graph->num_nodes
 *           where labels[i] is the component index (0..K-1) for
 *           active nodes and -1 for inactive nodes.  Sets *K_out to
 *           the total number of components found.
 *
 * Returns:  Pointer to labels array on success, NULL on allocation failure.
 */
static int*
find_components(SPLICE_GRAPH *graph, int *K_out)
{
  int  i, j;
  int  K   = 0;
  int  cur, nb;
  int  head, tail;
  int *labels = NULL;
  int *queue  = NULL;
  int  status;

  ESL_ALLOC(labels, sizeof(int) * graph->num_nodes);
  ESL_ALLOC(queue,  sizeof(int) * graph->num_nodes);

  for(i = 0; i < graph->num_nodes; i++) labels[i] = -1;

  for(i = 0; i < graph->num_nodes; i++) {
    if(!graph->node_in_graph[i]) continue;
    if(labels[i] != -1)          continue;  /* already labelled */

    /* BFS from node i to label component K */
    head = tail = 0;
    queue[tail++] = i;
    labels[i]     = K;

    while(head < tail) {
      cur = queue[head++];

      /* Follow outgoing edges */
      for(j = 0; j < graph->num_edges[cur]; j++) {
        nb = graph->edges[cur][j].downstream_node_id;
        if(nb < 0 || nb >= graph->num_nodes) continue;
        if(!graph->node_in_graph[nb])        continue;
        if(labels[nb] != -1)                 continue;
        labels[nb]    = K;
        queue[tail++] = nb;
      }

      /* Follow incoming edges */
      for(j = 0; j < graph->num_in_nodes[cur]; j++) {
        nb = graph->in_nodes[cur][j];
        if(nb < 0 || nb >= graph->num_nodes) continue;
        if(!graph->node_in_graph[nb])        continue;
        if(labels[nb] != -1)                 continue;
        labels[nb]    = K;
        queue[tail++] = nb;
      }
    }
    K++;
  }

  free(queue);
  *K_out = K;
  return labels;

  ERROR:
    if(labels) free(labels);
    if(queue)  free(queue);
    return NULL;
}


/* Function: p7_splicegraph_Split()
 *
 * Purpose:  Split <graph> into weakly-connected subgraphs.
 *
 *           If the graph is already connected (K==1), returns a
 *           1-element array containing the original <graph> pointer
 *           unchanged; *num_subgraphs is set to 1.  The caller should
 *           destroy the original graph normally and free the returned
 *           array pointer.
 *
 *           If K>1, returns an array of K independent subgraphs, one
 *           per component.  Each subgraph owns deep copies of all its
 *           data arrays.  Nodes not belonging to a component have
 *           node_in_graph[i]=FALSE in that subgraph.  The caller should
 *           call p7_splicegraph_DestroySubgraph() on each element and
 *           then free() the array.  The original <graph> should still
 *           be destroyed normally by the outer code.
 *
 * Returns:  Pointer to SPLICE_GRAPH* array on success, NULL on failure.
 */
SPLICE_GRAPH**
p7_splicegraph_Split(SPLICE_GRAPH *graph, int *num_subgraphs)
{
  int           i, k, K;
  int          *labels     = NULL;
  SPLICE_GRAPH *sub        = NULL;
  SPLICE_GRAPH **subgraphs = NULL;
  int           status;

  labels = find_components(graph, &K);
  if(labels == NULL) goto ERROR;

  /* K==1: graph is already connected, return it unchanged */
  if(K == 1) {
    ESL_ALLOC(subgraphs, sizeof(SPLICE_GRAPH*));
    subgraphs[0]   = graph;
    *num_subgraphs = 1;
    free(labels);
    return subgraphs;
  }

  /* K>1: deep-copy each component into an independent subgraph */
  ESL_ALLOC(subgraphs, sizeof(SPLICE_GRAPH*) * K);
  for(k = 0; k < K; k++) subgraphs[k] = NULL;

  for(k = 0; k < K; k++) {

    ESL_ALLOC(sub, sizeof(SPLICE_GRAPH));

    /* Initialize all pointers to NULL so DestroySubgraph is safe on partial alloc */
    sub->node_in_graph = NULL;
    sub->tmp_node      = NULL;
    sub->orig_hit_idx  = NULL;
    sub->best_out_edge = NULL;
    sub->ali_scores    = NULL;
    sub->path_scores   = NULL;
    sub->th            = NULL;
    sub->edges         = NULL;
    sub->edge_mem      = NULL;
    sub->num_edges     = NULL;
    sub->in_nodes      = NULL;
    sub->in_node_mem   = NULL;
    sub->num_in_nodes  = NULL;

    subgraphs[k] = sub;   /* store early so ERROR cleanup works */

    sub->nalloc    = graph->nalloc;
    sub->num_nodes = graph->num_nodes;
    sub->anchor_N  = graph->anchor_N;
    sub->revcomp   = graph->revcomp;
    sub->seqidx    = graph->seqidx;
    sub->seqname   = graph->seqname;  /* shared pointer, not owned */

    ESL_ALLOC(sub->node_in_graph, sizeof(int)   * sub->nalloc);
    ESL_ALLOC(sub->tmp_node,      sizeof(int)   * sub->nalloc);
    ESL_ALLOC(sub->orig_hit_idx,  sizeof(int)   * sub->nalloc);
    ESL_ALLOC(sub->best_out_edge, sizeof(int)   * sub->nalloc);
    ESL_ALLOC(sub->ali_scores,    sizeof(float) * sub->nalloc);
    ESL_ALLOC(sub->path_scores,   sizeof(float) * sub->nalloc);

    ESL_ALLOC(sub->th,     sizeof(P7_TOPHITS));
    sub->th->N                   = 0;
    sub->th->is_sorted_by_seqidx = FALSE;
    sub->th->hit                 = NULL;
    ESL_ALLOC(sub->th->hit, sizeof(P7_HIT*) * sub->nalloc);

    ESL_ALLOC(sub->edges,     sizeof(SPLICE_EDGE*) * sub->nalloc);
    ESL_ALLOC(sub->edge_mem,  sizeof(int)          * sub->nalloc);
    ESL_ALLOC(sub->num_edges, sizeof(int)          * sub->nalloc);

    ESL_ALLOC(sub->in_nodes,     sizeof(int*) * sub->nalloc);
    ESL_ALLOC(sub->in_node_mem,  sizeof(int)  * sub->nalloc);
    ESL_ALLOC(sub->num_in_nodes, sizeof(int)  * sub->nalloc);

    /* Pre-NULL the per-node pointer arrays so partial-fill is safe on error */
    for(i = 0; i < sub->nalloc; i++) {
      sub->edges[i]    = NULL;
      sub->in_nodes[i] = NULL;
    }

    for(i = 0; i < sub->nalloc; i++) {

      if(i < graph->num_nodes && graph->node_in_graph[i] && labels[i] == k) {

        sub->node_in_graph[i] = TRUE;
        sub->tmp_node[i]      = graph->tmp_node[i];
        sub->orig_hit_idx[i]  = graph->orig_hit_idx[i];
        sub->best_out_edge[i] = graph->best_out_edge[i];
        sub->ali_scores[i]    = graph->ali_scores[i];
        sub->path_scores[i]   = -eslINFINITY;
        sub->th->hit[i]       = graph->th->hit[i];
        sub->th->N++;

        /* Deep copy outgoing edges */
        sub->num_edges[i] = graph->num_edges[i];
        sub->edge_mem[i]  = graph->edge_mem[i];
        if(graph->edges[i] != NULL) {
          ESL_ALLOC(sub->edges[i], sizeof(SPLICE_EDGE) * sub->edge_mem[i]);
          memcpy(sub->edges[i], graph->edges[i], sizeof(SPLICE_EDGE) * sub->num_edges[i]);
        }

        /* Deep copy incoming edge list */
        sub->num_in_nodes[i] = graph->num_in_nodes[i];
        sub->in_node_mem[i]  = graph->in_node_mem[i];
        if(graph->in_nodes[i] != NULL) {
          ESL_ALLOC(sub->in_nodes[i], sizeof(int) * sub->in_node_mem[i]);
          memcpy(sub->in_nodes[i], graph->in_nodes[i], sizeof(int) * sub->num_in_nodes[i]);
        }

      } else {

        sub->node_in_graph[i] = FALSE;
        sub->tmp_node[i]      = FALSE;
        sub->orig_hit_idx[i]  = -1;
        sub->best_out_edge[i] = -1;
        sub->ali_scores[i]    = -eslINFINITY;
        sub->path_scores[i]   = -eslINFINITY;
        sub->th->hit[i]       = NULL;
        sub->num_edges[i]     = 0;
        sub->edge_mem[i]      = 0;
        sub->num_in_nodes[i]  = 0;
        sub->in_node_mem[i]   = 0;
      }
    }
  }

  free(labels);
  *num_subgraphs = K;
  return subgraphs;

  ERROR:
    if(labels)    free(labels);
    if(subgraphs) {
      for(k = 0; k < K; k++)
        if(subgraphs[k] != NULL) p7_splicegraph_DestroySubgraph(subgraphs[k]);
      free(subgraphs);
    }
    return NULL;
}


/* Function: p7_splicegraph_DestroySubgraph()
 *
 * Purpose:  Frees a subgraph created by p7_splicegraph_Split() for K>1.
 *           All arrays are independently owned and are freed here.
 *           The seqname pointer is not owned and is not freed.
 *           Hit objects with tmp_node[i]==TRUE (created during splicing)
 *           are freed; all other hit objects are not freed.
 */
void
p7_splicegraph_DestroySubgraph(SPLICE_GRAPH *graph)
{
  int i;

  if(graph == NULL) return;

  if(graph->orig_hit_idx  != NULL) free(graph->orig_hit_idx);
  if(graph->node_in_graph != NULL) free(graph->node_in_graph);

  if(graph->path_scores   != NULL) free(graph->path_scores);
  if(graph->ali_scores    != NULL) free(graph->ali_scores);

  if(graph->best_out_edge != NULL) free(graph->best_out_edge);

  if(graph->tmp_node != NULL && graph->th != NULL && graph->th->hit != NULL) {
    for(i = 0; i < graph->num_nodes; i++) {
      if(graph->tmp_node[i]) {
        p7_trace_fs_Destroy(graph->th->hit[i]->dcl->tr);
        p7_hit_Destroy(graph->th->hit[i]);
      }
    }
  }

  if(graph->tmp_node != NULL) free(graph->tmp_node);

  if(graph->th != NULL) {
    if(graph->th->hit != NULL) free(graph->th->hit);
    free(graph->th);
  }

  if(graph->edges != NULL) {
    for(i = 0; i < graph->nalloc; i++)
      if(graph->edges[i] != NULL) free(graph->edges[i]);
    free(graph->edges);
  }
  if(graph->edge_mem  != NULL) free(graph->edge_mem);
  if(graph->num_edges != NULL) free(graph->num_edges);

  if(graph->in_nodes != NULL) {
    for(i = 0; i < graph->nalloc; i++)
      if(graph->in_nodes[i] != NULL) free(graph->in_nodes[i]);
    free(graph->in_nodes);
  }
  if(graph->in_node_mem  != NULL) free(graph->in_node_mem);
  if(graph->num_in_nodes != NULL) free(graph->num_in_nodes);

  graph->seqname = NULL;

  free(graph);
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
 
  if(graph->edges[up_node] == NULL) {
    ESL_ALLOC(graph->edges[up_node], sizeof(SPLICE_EDGE) * EDGE_ALLOC);
    graph->edge_mem[up_node]  = EDGE_ALLOC;
    graph->num_edges[up_node] = 0;   
  }

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

  ret_edge->jump_edge  = FALSE;
  ret_edge->edge_score = 0.;

  /* Add up_node to the incoming node list of down_node */
  if(graph->in_nodes[down_node] == NULL) {
    ESL_ALLOC(graph->in_nodes[down_node], sizeof(int) * EDGE_ALLOC);
    graph->in_node_mem[down_node]  = EDGE_ALLOC;
    graph->num_in_nodes[down_node] = 0;
  }
  if(graph->num_in_nodes[down_node] == graph->in_node_mem[down_node]) {
    graph->in_node_mem[down_node] *= 2;
    ESL_REALLOC(graph->in_nodes[down_node], sizeof(int) * graph->in_node_mem[down_node]);
  }
  graph->in_nodes[down_node][graph->num_in_nodes[down_node]++] = up_node;

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
       if(graph->edges[up_node][i].edge_score != -eslINFINITY) return TRUE;
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
p7_splicegraph_AliScoreEdge(SPLICE_EDGE *edge, const P7_DOMAIN *upstream_dom, const P7_DOMAIN *downstream_dom)
{

  int p, s;
  int overlap_start;
  int overlap_end;
  int overlap_len;
  int last_k;
  float min_lost_sc;
  float curr_lost_sc;
  float upstream_lost;
  float downstream_lost;
  int   *kpp;
  float *spp;
  float *upstream_suffix_sum;
  float *downstream_prefix_sum;
  int status;

  upstream_suffix_sum   = NULL;
  downstream_prefix_sum = NULL;

  /* return if no there is no hmm overlap */
  if(downstream_dom->ihmm > upstream_dom->jhmm)  return eslOK;

  overlap_start = ESL_MAX(upstream_dom->ihmm, downstream_dom->ihmm);
  overlap_end   = ESL_MIN(upstream_dom->jhmm, downstream_dom->jhmm);
  
  overlap_len   = overlap_end - overlap_start + 1;

  if(overlap_len < 1) {
    edge->edge_score = -eslINFINITY;
    return eslOK;
  }

  ESL_ALLOC(upstream_suffix_sum,   sizeof(float) * overlap_len);
  ESL_ALLOC(downstream_prefix_sum, sizeof(float) * overlap_len);

  esl_vec_FSet(upstream_suffix_sum,   overlap_len, 0.0);
  esl_vec_FSet(downstream_prefix_sum, overlap_len, 0.0);

  /* Fill the upstream array */
  spp = upstream_dom->scores_per_pos;
  kpp = upstream_dom->k_per_pos;

  /* Move to the end of the overlap in the scores_per_pos array */
  p = upstream_dom->per_pos_len-1;
  while(p >= 0 && kpp[p] != overlap_end) p--;
  if(p == -1) ESL_XEXCEPTION(eslFAIL, "Edge Scoring Failed"); 

  /* Add all scores in the scores_per_pos for each k to the same index in the upstream_suffix_sum array */
  last_k = overlap_end;
  s = overlap_len - 1;
  upstream_suffix_sum[s] += spp[p];
  p--;
  while(p >=0 && kpp[p] >= overlap_start) {
    if(kpp[p] != last_k) s--;
    last_k = kpp[p]; 

    upstream_suffix_sum[s] += spp[p];
    p--;
  }
  
  /* Create running sum of suffix scores */
  for(s = overlap_len - 2; s >= 0; s--)
    upstream_suffix_sum[s] += upstream_suffix_sum[s+1];

  /* In the case that the upstream hit extends past the end fo the downstream hit, calculate any lost score from the end of the upstream hit */
  upstream_lost = 0.;
  if(upstream_dom->jhmm > overlap_end) {
    p = upstream_dom->per_pos_len-1;
    while(kpp[p] > overlap_end) {
      upstream_lost += spp[p];
      p--;
    } 
  }

  /* Fill the downstream array */
  spp = downstream_dom->scores_per_pos;
  kpp = downstream_dom->k_per_pos;

  /* Move to the start of the overlap in the scores_per_pos array */
  p = 0;
  while(p < downstream_dom->per_pos_len && kpp[p] != overlap_start) p++;
  if(p == downstream_dom->per_pos_len) ESL_XEXCEPTION(eslFAIL, "Edge Scoring Failed");

  /* Add all scores in the scores_per_pos for each k to the same index in the downstream_prefix_sum array */
  last_k = overlap_start;
  s = 0;
  downstream_prefix_sum[s] += spp[p];
  p++;
  while(p < downstream_dom->per_pos_len && kpp[p] <= overlap_end) {
    if(kpp[p] != last_k) s++;
    last_k = kpp[p];

    downstream_prefix_sum[s] += spp[p];
    p++;
  }

  /* Create running sum of prefix scores */
  for(s = 1; s < overlap_len; s++)
    downstream_prefix_sum[s] +=  downstream_prefix_sum[s-1];

  /* In the case that the downstream hit extends past the start fo the upstream hit, calculate any lost score from the start of the downstream hit */
  downstream_lost = 0.;
  if(downstream_dom->ihmm < overlap_start) {
    p = 0;
    while(kpp[p] < overlap_start) {
      downstream_lost += spp[p];
      p++;
    }
  }


  /* Find the minimum score loss to eliminate the overlap*/
  /*start with all positions belonging to the downstream hit - not allowed if upstream hit heas no positions before the overlap. */
  if(upstream_dom->ihmm == overlap_start) min_lost_sc = eslINFINITY;
  else                                    min_lost_sc = upstream_suffix_sum[0]; 
  for(s = 1; s < overlap_len; s++) {
    /* at each step add another overlap postion to the upstream hit */
    curr_lost_sc = upstream_suffix_sum[s] + downstream_prefix_sum[s-1]; 
    min_lost_sc = ESL_MIN(min_lost_sc, curr_lost_sc);
  }
  /* end with all positions belonging to the upstream hit - not allowed if downstream hit has no postions after the overlap */
  if(downstream_dom->jhmm > overlap_end) min_lost_sc = ESL_MIN(min_lost_sc, downstream_prefix_sum[overlap_len-1]);

  edge->edge_score -= (min_lost_sc + upstream_lost + downstream_lost);

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


int
p7_splicegraph_NodeOverlap(SPLICE_GRAPH *graph, int node_id, SPLICE_PATH *path, int step_id) {

  int     overlap_hmm_start;
  int     overlap_hmm_end;
  int64_t overlap_seq_start;
  int64_t overlap_seq_end;

  overlap_hmm_start = ESL_MAX(graph->th->hit[node_id]->dcl->ihmm, path->ihmm[step_id]);
  overlap_hmm_end   = ESL_MIN(graph->th->hit[node_id]->dcl->jhmm, path->jhmm[step_id]);
  
  if(overlap_hmm_end - overlap_hmm_start + 1 <= 0) return FALSE;

  if(graph->revcomp) {
     overlap_seq_start = ESL_MAX(graph->th->hit[node_id]->dcl->jali, path->jali[step_id]);
     overlap_seq_end   = ESL_MIN(graph->th->hit[node_id]->dcl->iali, path->iali[step_id]);
  }
  else {
    overlap_seq_start = ESL_MAX(graph->th->hit[node_id]->dcl->iali, path->iali[step_id]);
    overlap_seq_end   = ESL_MIN(graph->th->hit[node_id]->dcl->jali, path->jali[step_id]);
  }
  
  if(overlap_seq_end - overlap_seq_start + 1 <= 0) return FALSE;

  return TRUE;

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

  fprintf(fp, "   %4s %6s %9s %12s %12s %10s %10s\n",
          "Hit", "hmm_from", "hmm_to", "seq_from", "seq_to", "score", "p-value");

  for(i = 0; i < graph->th->N; i++) {

    if(!graph->node_in_graph[i]) continue;

    hit = graph->th->hit[i];

    if ( i < graph->anchor_N)
      fprintf(fp, "   %4d %6d %9d %12" PRId64 " %12" PRId64 " %10.2f %10f\n",
              i+1, hit->dcl->ihmm, hit->dcl->jhmm, hit->dcl->iali, hit->dcl->jali, hit->dcl->aliscore, exp(hit->sum_lnP));
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
      edge_scores[tmp_edge->downstream_node_id] = tmp_edge->edge_score;
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


