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

#include "esl_gumbel.h"
#include "easel.h"
#include "hmmer.h"


/************************************************************************
 * 1. Macros, sturcts and functions for allocation, initialization, destruction.
 ************************************************************************/

static float SSSCORE[2]      = {-0.7,0.0}; // Non-canon vs canon splice site
static float EDGE_FAIL_SCORE = -14773.0;   // Makes me thirsty for a latte!

static char  AMINO_CHARS[21] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-'};
static char    DNA_CHARS[ 6] = {'A','C','G','T','-','N'};
static char LC_DNA_CHARS[ 6] = {'a','c','g','t','-','n'};


//static int MAX_TARGET_RANGE_EXT = 3000000;

// How many amino acids are we willing to extend to bridge two hits?  
// How many overlapping aminos do we require to perform bridging?
static int MAX_AMINO_EXT     = 6;
static int MIN_AMINO_OVERLAP = 6;


// When we're performing sub-model search, what's the maximum size
// range we're willing to consider?
static int MAX_SUB_MODEL_RANGE =    100;
static int MAX_SUB_NUCL_RANGE  = 100000;


int intMax (int a, int b) { if (a>b) return a; return b; }
int intMin (int a, int b) { if (a<b) return a; return b; }


/* Struct used to store the coodinates of ranges on a target strand 
 * and pointers to the hits inside that range.                      */
typedef struct _target_range {

  int64_t      seqidx;            /* target squence id       */
  int64_t      start;             /* range start position    */
  int64_t      end;               /* range end postion       */
  int          complementarity;   /* reverse complementarity */
  char        *seqname;          /* target sequnce name     */
  P7_TOPHITS  *th;                /* hits in range           */

} TARGET_RANGE;


/* Struct used to store pointers to TARGET_RANGE objects */
typedef struct _target_range_set {

  int            N;              /* number of target ranges in set */
  int            Nalloc;         /* current allocation size        */
  TARGET_RANGE **target_range;   /* array of target ranges         */

} TARGET_RANGE_SET;


/* Struct used to store all information relevant to splice 
 * graph edge between upstream and downstream hits */
typedef struct _splice_edge {

  const ESL_ALPHABET * ntalpha;

  int amino_start;
  int amino_end;

  int upstream_hit_id;
  int upstream_nucl_start;
  int upstream_nucl_end;
  int upstream_ext_len;
  int upstream_trace_start;
  int upstream_trace_end;

  P7_TOPHITS    * UpstreamTopHits;
  P7_ALIDISPLAY * UpstreamDisplay;  //Delete
  P7_TRACE      * UpstreamTrace;
  ESL_DSQ       * UpstreamNucls;   // Possibly delete and replace with indicies to original target_seq

  int downstream_hit_id;
  int downstream_nucl_start;
  int downstream_nucl_end;
  int downstream_ext_len;
  int downstream_trace_start;
  int downstream_trace_end;

  P7_TOPHITS    * DownstreamTopHits;
  P7_ALIDISPLAY * DownstreamDisplay;  // Delete
  P7_TRACE      * DownstreamTrace;
  ESL_DSQ       * DownstreamNucls;    // Possibly delete and replace with indicies to original target_seq

  int   upstream_exon_terminus;
  int downstream_exon_terminus;

  int   upstream_spliced_nucl_end;
  int downstream_spliced_nucl_start;

  float score_density;
  float score;

} SPLICE_EDGE;


typedef struct _splice_node {


  int node_id;


  int was_missed;
  int hit_id;

  int in_edge_cap;
  int num_in_edges;
  int best_in_edge; // Relative to the following arrays
  SPLICE_EDGE ** InEdges;
  struct _splice_node ** UpstreamNodes;


  int out_edge_cap;
  int num_out_edges;
  SPLICE_EDGE ** OutEdges;
  struct _splice_node ** DownstreamNodes;


  int is_n_terminal;
  int is_c_terminal;


  float hit_score;
  float cumulative_score;


} SPLICE_NODE;


/*struct that conatins splice nodes and edges as well as 
 * other information to find optimal splice paths */
typedef struct _splice_graph {

  P7_TOPHITS  *    TopHits;
  P7_TOPHITS  * MissedHits; // This is a fake P7_TOPHITS!  Don't be fooled!

  P7_PROFILE  *  Model;
  P7_OPROFILE * OModel;

  int revcomp;

  // Because we build based on the SPLICE_EDGE datastructure,
  // it's helpful to be able to find each node ID by way of
  // these lookups.
  int * TH_HitDomToNodeID; //    TopHits
  int * MH_HitToNodeID;    // MissedHits


  int num_nodes;
  int num_edges;
  SPLICE_NODE ** Nodes; // NOTE: these go from [1..num_nodes]


  int   num_n_term;
  int   num_c_term;
  int * CTermNodeIDs; // When we generate this, we sort by cumulative score

  int   has_full_path;
  int   best_full_path_length;
  int   best_full_path_start;
  int   best_full_path_end;
  float best_full_path_score;


} SPLICE_GRAPH;


/* Free a TARGET_RANGE object */
void
Target_Range_Destroy(TARGET_RANGE *tr)
{

  if (tr == NULL) return;
  
  tr->seqname = NULL;
  if (tr->th != NULL)
    p7_tophits_Destroy(tr->th);

  free(tr);

  return;

}


/* Create a TARGET_RANGE object with room for nalloc P7_HIT pointers*/
TARGET_RANGE *
Target_Range_Create(int nalloc)
{
  int status;

  TARGET_RANGE *tr = NULL;
  ESL_ALLOC(tr, sizeof(TARGET_RANGE));

  tr->seqname = NULL;

  tr->th = NULL;
  ESL_ALLOC(tr->th, sizeof(P7_TOPHITS));
 
  tr->th->hit    = NULL;
  ESL_ALLOC(tr->th->hit,  nalloc * sizeof(P7_HIT *));

  tr->th->unsrt  = NULL;
  tr->th->N      = 0;
  tr->th->Nalloc = nalloc;

  return tr;

  ERROR:
    Target_Range_Destroy(tr);
    return NULL;   
}

/* free a TARGET_RANGE_SET object */
void
Target_Range_Set_Destroy(TARGET_RANGE_SET *trs)
{

  int i;

  if (trs == NULL) return;
  if (trs->target_range != NULL)
  {
    for(i=0; i< trs->N; i++) 
      Target_Range_Destroy(trs->target_range[i]);
    
    free(trs->target_range);
  }

  free(trs);

  return;

}

/* Create a TARGET_RANGE_SET object */
TARGET_RANGE_SET *
Target_Range_Set_Create(void)
{
  int status;

  TARGET_RANGE_SET *trs = NULL;
  ESL_ALLOC(trs, sizeof(TARGET_RANGE_SET));
   
  trs->target_range = NULL;
  ESL_ALLOC(trs->target_range, 100 * sizeof(TARGET_RANGE*));

  trs->N      = 0;
  trs->Nalloc = 100;

  return trs;

  ERROR:
    Target_Range_Set_Destroy(trs);
    return NULL;
}

/* Double the size of a TARGET_RANGE_SET object */
int
Target_Range_Set_Grow(TARGET_RANGE_SET *trs)
{
  void *p;
  int Nalloc = trs->Nalloc * 2;
  int status;

  if (trs->N < trs->Nalloc) return eslOK;
 
  ESL_RALLOC(trs->target_range, p, Nalloc * sizeof(TARGET_RANGE*));
  trs->Nalloc = Nalloc;

  return eslOK;

  ERROR:
    return eslEMEM;

}

/* Add new TARGET_RANGE pointer to TARGET_RANGE_SET */
int
Target_Range_Set_CreateNextRange(TARGET_RANGE_SET *trs, TARGET_RANGE **ret_range, int nalloc)
{
  int status;
  
  if ((status = Target_Range_Set_Grow(trs)) != eslOK) goto ERROR;  

  *ret_range = Target_Range_Create(nalloc);
  trs->target_range[trs->N] = *ret_range;
  trs->N++;

  return eslOK;

  ERROR:
    ret_range = NULL;
    return status;
}

/* Free a SPLICE_EDGE pointer */
void SPLICE_EDGE_Destroy
(SPLICE_EDGE * edge)
{
  
  free(edge->UpstreamNucls);
  free(edge->DownstreamNucls);
  free(edge);
}

/* Free a SPLICE_GRAPH and all it nodes and edges */
void SPLICE_GRAPH_Destroy
(SPLICE_GRAPH * Graph)
{
  
  int node_id, in_edge_index; //, out_edge_index;

  // Because we don't want to accidentally double-free any
  // of the SPLICE_EDGE structs, we'll put each node in
  // charge of its incoming DOs.
  if (Graph->Nodes) {
    
    SPLICE_NODE * Node;
     int tmp = 0;
    for (node_id=1; node_id<=Graph->num_nodes; node_id++) {
      Node = Graph->Nodes[node_id];
      free(Node->DownstreamNodes);
      free(Node->UpstreamNodes);
      //printf("Node->num_in_edges %d\n", Node->num_in_edges);
      for (in_edge_index = 0; in_edge_index < Node->num_in_edges; in_edge_index++) { 
        SPLICE_EDGE_Destroy(Node->InEdges[in_edge_index]);
        tmp++;
       }

      free(Node->InEdges);
      free(Node->OutEdges);
      free(Node);

    }
    free(Graph->Nodes);

    if (Graph->CTermNodeIDs)
      free(Graph->CTermNodeIDs);
  }

  int hit_id;
  if (Graph->TH_HitDomToNodeID) {
    free(Graph->TH_HitDomToNodeID);
  
    if (Graph->MH_HitToNodeID) 
      free(Graph->MH_HitToNodeID);
  }

  // We need to manually destroy the fake 'MissedHits'
  if (Graph->MissedHits) {

    for (hit_id=0; hit_id<Graph->MissedHits->N; hit_id++) {

      // Luckily, the AD is for real (generally) so we should
      // be able to use the generic destructor 
      if(Graph->MissedHits->hit[hit_id]->dcl->ad)
        if(Graph->MissedHits->hit[hit_id]->dcl->ad->ntseq) 
          free(Graph->MissedHits->hit[hit_id]->dcl->ad->ntseq);
      p7_alidisplay_Destroy(Graph->MissedHits->hit[hit_id]->dcl->ad);
      p7_trace_fs_Destroy(Graph->MissedHits->hit[hit_id]->dcl->tr);
      free(Graph->MissedHits->hit[hit_id]->dcl);
      free(Graph->MissedHits->hit[hit_id]);

    }
    free(Graph->MissedHits->hit);
    free(Graph->MissedHits);

  }


  free(Graph);
  Graph = NULL;

}

void SPLICE_NODE_Destroy
(SPLICE_NODE * Node)
{

  int in_edge_id;
  for (in_edge_id = 0; in_edge_id < Node->num_in_edges; in_edge_id++) {
    if (Node->InEdges[in_edge_id]) {
      SPLICE_EDGE_Destroy(Node->InEdges[in_edge_id]);
      Node->InEdges[in_edge_id] = NULL;
    }
  }

  int out_edge_id;
  for (out_edge_id = 0; out_edge_id < Node->num_out_edges; out_edge_id++) {
    if (Node->OutEdges[out_edge_id]) {
      SPLICE_EDGE_Destroy(Node->OutEdges[out_edge_id]);
      Node->OutEdges[out_edge_id] = NULL;
    }
  }

  free(Node->InEdges);
  free(Node->OutEdges);
  free(Node);
  Node = NULL;

}


SPLICE_NODE * SpliceNode_Create
(
  SPLICE_GRAPH * Graph,
  int node_id, 
  int hit_id, 
  int was_missed
)
{

  SPLICE_NODE *NewNode;
  int status;
  
  ESL_ALLOC(NewNode, sizeof(SPLICE_NODE));
  
  NewNode->node_id = node_id;
  NewNode->hit_id  = hit_id;

  NewNode->out_edge_cap    = 10;
  NewNode->num_out_edges   =  0;
  ESL_ALLOC(NewNode->OutEdges,        NewNode->out_edge_cap*sizeof(SPLICE_EDGE *));
  ESL_ALLOC(NewNode->DownstreamNodes, NewNode->out_edge_cap*sizeof(SPLICE_NODE *));

  NewNode->in_edge_cap    = 10;
  NewNode->num_in_edges   =  0;
  NewNode->best_in_edge   = -1;
  ESL_ALLOC(NewNode->InEdges,         NewNode->in_edge_cap*sizeof(SPLICE_EDGE *));
  ESL_ALLOC(NewNode->UpstreamNodes,   NewNode->in_edge_cap*sizeof(SPLICE_NODE *));

  NewNode->was_missed = was_missed;

  if (was_missed) {
  
    NewNode->is_n_terminal = 0;
    if (Graph->MissedHits->hit[hit_id]->dcl->tr->hmmfrom[0] == 1)
      NewNode->is_n_terminal = 1;

    NewNode->is_c_terminal = 0;
    if (Graph->MissedHits->hit[hit_id]->dcl->tr->hmmto[0] == Graph->Model->M)
      NewNode->is_c_terminal = 1;

    NewNode->hit_score = Graph->MissedHits->hit[hit_id]->dcl->bitscore;

  } else {
  
    NewNode->is_n_terminal = 0;
    if ((&(Graph->TopHits->hit[hit_id]->dcl[0]))->tr->hmmfrom[0] == 1)
      NewNode->is_n_terminal = 1;

    NewNode->is_c_terminal = 0;
    if ((&(Graph->TopHits->hit[hit_id]->dcl[0]))->tr->hmmto[0] == Graph->Model->M)
      NewNode->is_c_terminal = 1;

    NewNode->hit_score = Graph->TopHits->hit[hit_id]->dcl[0].bitscore;

  }

  // Initialize this score to a recognizable and impossibly low value
  NewNode->cumulative_score = EDGE_FAIL_SCORE;

  return NewNode;


  ERROR:
   SPLICE_NODE_Destroy(NewNode);
   return NULL;
    
}

int
SpliceNode_Grow(SPLICE_NODE *node)
{
 
  int i;
  SPLICE_EDGE **NewOutEdges;
  SPLICE_EDGE **NewInEdges;
  SPLICE_NODE **NewDSNodes;
  SPLICE_NODE **NewUSNodes;
  int status;

  if (node->num_out_edges == node->out_edge_cap) {
    node->out_edge_cap *= 2;
  
    ESL_ALLOC(NewOutEdges, node->out_edge_cap * sizeof(SPLICE_EDGE *));
    ESL_ALLOC(NewDSNodes,  node->out_edge_cap * sizeof(SPLICE_NODE *));

    for (i=0; i<node->num_out_edges; i++) {
      NewOutEdges[i] = node->OutEdges[i];
      NewDSNodes[i]  = node->DownstreamNodes[i];
    }

    free(node->OutEdges);
    free(node->DownstreamNodes);

    node->OutEdges = NewOutEdges;
    node->DownstreamNodes = NewDSNodes;
  }

  if (node->num_in_edges == node->in_edge_cap) {

    node->in_edge_cap *= 2;

    ESL_ALLOC(NewInEdges, node->in_edge_cap * sizeof(SPLICE_EDGE *));
    ESL_ALLOC(NewUSNodes, node->in_edge_cap * sizeof(SPLICE_NODE *));

    for (i=0; i<node->num_in_edges; i++) {
      NewInEdges[i] = node->InEdges[i];
      NewUSNodes[i] = node->UpstreamNodes[i];
    }

    free(node->InEdges);
    free(node->UpstreamNodes);

    node->InEdges = NewInEdges;
    node->UpstreamNodes = NewUSNodes;
  }

  return eslOK;

  ERROR:
    return eslEMEM;

}


SPLICE_GRAPH *
Splice_Graph_Create(void)
{
  int status;
  SPLICE_GRAPH *graph;
  
  ESL_ALLOC(graph, sizeof(SPLICE_GRAPH));

  graph->num_nodes  = 0;
  graph->num_edges  = 0;
  graph->num_n_term = 0;
  graph->num_c_term = 0;

  graph->has_full_path        = 0;
  graph->best_full_path_start = 0;
  graph->best_full_path_end   = 0;
  graph->best_full_path_score = 0.0;

  graph->TopHits           = NULL;
  graph->MissedHits        = NULL;
  graph->Model             = NULL;
  graph->OModel            = NULL;
  graph->TH_HitDomToNodeID = NULL;
  graph->MH_HitToNodeID    = NULL;
  graph->Nodes             = NULL;
  graph->CTermNodeIDs      = NULL;

  return graph;

ERROR:
    SPLICE_GRAPH_Destroy(graph);
    return NULL;


}
/* Internal functions */
static TARGET_RANGE_SET* build_target_ranges       (const P7_TOPHITS *th);
static void              get_target_range_coords   (TARGET_RANGE_SET* targets, P7_HIT **hits, int *sort_score_idx, const int num_hits, const int complementarity);
static ESL_SQ* get_target_range_sequence (ESL_SQFILE * genomic_seq_file, TARGET_RANGE *target_range);

static SPLICE_GRAPH* build_splice_graph ( TARGET_RANGE *target_range, ESL_SQ *target_seq, P7_TOPHITS *TopHits, P7_PROFILE *gm, P7_OPROFILE *om, ESL_GENCODE *gcode);
static SPLICE_EDGE ** create_splice_edges ( TARGET_RANGE *target_range, ESL_SQ *target_seq, P7_PROFILE *gm, ESL_GENCODE *gcode, int *num_splice_edges);
static int are_hits_splice_comaptible (P7_DOMAIN *Upstream, P7_DOMAIN *Downstream, int complementarity);
static void sketch_splice_edge (SPLICE_EDGE *Edge, ESL_SQ *target_seq, P7_PROFILE *gm, ESL_GENCODE *gcode, int complementarity);
static void splice_edge_overlap (SPLICE_EDGE *Overlap, P7_PROFILE *gm, ESL_GENCODE *gcode);
static void get_nuc_coords_form_amino_overlap (SPLICE_EDGE * Edge, int complementarity);
static ESL_DSQ* grab_seq_range (ESL_SQ* target_seq, int range_start, int range_end, int complementarity);
static float find_optimal_splice_site ( SPLICE_EDGE *Overlap, P7_PROFILE *gm, ESL_GENCODE *gcode, int *upstream_splice_index, int *downstream_splice_index, int *split_amino_model_index, int *codon_split_option);
static float ali_score_at_postion ( P7_PROFILE *gm, int amino, int model_pos, int trans_pos, int prev_state, int curr_state);
static void fill_out_graph_structure (SPLICE_GRAPH *Graph, SPLICE_EDGE **SpliceEdges, int num_splice_edges, int complementarity);
static void add_missing_exons( SPLICE_GRAPH *Graph, ESL_SQ *target_seq, ESL_GENCODE *gcode);
static P7_TOPHITS* seek_missing_exons( SPLICE_GRAPH *Graph, ESL_SQ *target_seq, ESL_GENCODE *gcode);
static int * get_bounded_search_regions( SPLICE_GRAPH *Graph, ESL_SQ *target_seq);
static int * identify_term_search_regions( SPLICE_GRAPH *Graph, ESL_SQ *target_seq, int *NoInEdgeNodes, int num_no_in_edge, int *NoOutEdgeNodes, int num_no_out_edge);
static void run_model_on_exon_sets (SPLICE_GRAPH *Graph, ESL_SQ *target_seq, ESL_GENCODE *gcode, ESL_GETOPTS *go, FILE *ofp, int textw); 
static ESL_DSQ * grab_exon_coord_set_nucls (int *ExonCoordSet, ESL_SQ *target_seq, int complementarity, int *coding_region_len);
static void report_spliced_tophits (SPLICE_GRAPH *Graph, P7_TOPHITS *ExonSetTopHits,  P7_PIPELINE *ExonSetPipeline, ESL_SQ *target_seq, int *ExonCoordSet, int *exon_set_name_id, FILE *ofp, int textw); 
static void print_spliced_alignment (P7_ALIDISPLAY *AD, ESL_SQ *target_seq, int *ExonCoordSet, int exon_set_name_id, FILE *ofp, int textw);
//static void print_exon (EXON_DISPLAY_INFO *EDI, ESL_SQ *target_seq, int *ad_nucl_read_pos, int *ad_amino_read_pos, int *codon_pos);
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                           BEGIN SPLICING STUFF
 * 
 *  Catalogue of Ships
 *  ==================
 *
 *  Fleet 0: Structs
 *
 *  + SPLICE_EDGE    : Information related to how two HMMER hits are spliced together (effectively a splice edge)
 *  + SPLICE_NODE       : A node in the splice graph (corresponds to a HMMER hit)
 *  + SPLICE_GRAPH      : The container for all of our splice nodes and supporting data
 *  + EXON_DISPLAY_INFO : Content used for printing the final exon alignments
 *
 *
 *  Fleet 1: Data Preparation / Management
 *
 *  + FloatHighLowSortIndex : Generate an index that sorts an array of floats in descending order
 *  + FloatLowHighSortIndex : Generate an index that sorts an array of floats in  ascending order
 *  + GetMinAndMaxCoords    : [DEPRECATED] Get the minimum and maximum nucleotide coordinates in a P7_TOPHITS set
 *  + SetTargetSeqRange     : Given a P7_TOPHITS object, predict what the true coding region is for a protein
 *  + get_target_range_sequence      : Extract the nucleotides corresponding to the range chosen by 'SetTargetSeqRange'
 *  + GrabNuclRange         : Get the nucleotide sequence of a window within the range extracted by 'get_target_range_sequence'
 *  
 *
 *
 *  Fleet 2: Determining How to Splice Pairs of Hits
 *
 *  + SelectSpliceOpt              : Given a model position and nucleotides, determine an optimal splice junction
 *  + GetContestedUpstreamNucls    : Pull a nucleotide sequence corresponding to a candidate 5' splice site
 *  + GetContestedDownstreamNucls  : Pull a nucleotide sequence corresponding to a candidate 3' splice site
 *  + AminoScoreAtPosition         : ASSUMING THAT WE'RE MOVING THROUGH THE MODEL, determine a positional score
 *  + FindOptimalSpliceSite        : Examine a SPLICE_EDGE to determine an optimal splicing of two P7_HITs
 *  + SpliceOverlappingDomains     : Do a bit of bookkeeping around 'FindOptimalSpliceSite'
 *  + GetNuclRangesFromAminoCoords : Given two P7_HITs that we're considering splicing, determine their nucleotide coordinate ranges
 *  + SketchSpliceEdge             : Perform any necessary extension to make two P7_HITs splice-able and call 'SpliceOverlappingDomains'
 *  + HitsAreSpliceComaptible      : Determine whether two hits are positioned correctly in the model/target to be spliced
 *  + ExcessiveGapContent          : Determine whether a hit is too gappy to be a strong initial candidate for splicing (can be recovered in sub-model search)
 *  + OutsideSearchArea            : Determine whether a hit is outside of the area selected in 'SetTargetSeqRange'
 *  + create_splice_edges      : Given a P7_TOPHITS, identify all pairs of hits that might be splice-able
 *
 *
 *
 *  Fleet 3: Producing a Splice Graph According to the Spliced Hit Pairs
 *
 *  + InitSpliceNode        : Allocate and initialize a node in our splice graph
 *  + ConnectNodesByEdge    : Given two nodes implicated in a satisfactorily spliced SPLICE_EDGE, draw an edge in the splice graph
 *  + GatherCTermNodes      : Derive a list of all nodes in the graph whose corresponding hits terminate at the C-terminal end of the model
 *  + EdgeWouldEraseNode    : Check if the best path through a node would (effectively) skip that node
 *  + PullUpCumulativeScore : Recursively determine the best cumulative score leading up to each node in the graph
 *  + EvaluatePaths         : Prepare for and call 'PullUpCumulativeScore' so we can determine the best path through the graph
 *  + FillOutGraphStructure : Initialize the splice graph and call 'EvaluatePaths' to produce our complete splice graph
 *  + FindBestFullPath      : Given a complete splice graph produced by 'FillOutGraphStructure' find an optimal path from model positions 1..M
 *  + build_splice_graph      : Top-level function for splice graph construction.  Calls 'create_splice_edges', 'FillOutGraphStructure', and 'FindBestFullPath'
 *
 *
 *
 *  Fleet 4: Filling in Gaps in the Splice Graph (Sub-Model Search)
 *
 *  + ExtractSubProfile         : Given a P7_PROFILE, extract a portion of that model for windowed Viterbi
 *  + NodesAreDCCCompatible     : Could a pair of nodes be upstream/downstream partners in a "disconnected component"?
 *  + IdentifyDisConnComponents : Determine regions where sub-model search might fill in "disconnected components"
 *  + IdentifyMidSearchRegions  : Determine regions where possible homologous overextension might benefit from sub-model search
 *  + IdentifyTermSearchRegions : Determine regions where we might need sub-model search to recover terminal exons
 *  + GetBoundedSearchRegions   : Call each of the 'Identify___' functions and organize sub-model search
 *  + SelectFinalSubHits        : Given a collection of sub-model hits, determine which ones we like as candidate exons
 *  + AminoToIndex              : [PROBABLY REDUNDANT WITH SOMETHING] Go from an amino acid character to a numeric index
 *  + ComputeRoughAliScore      : Use 'AminoScoreAtPosition' to get a quick sense of the score a Viterbi alignment
 *  + FindSubHits               : Given ranges on the model and the nucleotide target, look for high-scoring missed exons using Viterbi
 *  + IntegrateMissedHits       : Add the missed exons chosen by 'SelectFinalSubHits' to the splice graph
 *  + SeekMissingExons          : Manage the search for missed exons, calls 'GetBoundedSearchRegions' and 'FindSubHits'
 *  + AddMissingExonsToGraph    : Manage filling in holes in the splice graph, calls 'SeekMissingExons' and 'IntegrateMissedHits' (among others)
 *
 *
 *
 *  Fleet 5: Producing the Final "Exon Set" Alignments
 *
 *  + GetExonSetFromEndNode     : Given a node in the splice graph, trace an optimal path to that node and record the implicated model and target coordinates
 *  + FindComponentBestEnd      : Determine which node in a connected component has the best cumulative score
 *  + TranslateExonSetNucls     : Translate a nucleotide sequence to amino acids (to serve as the final target for Forward/Backward/Stuff)
 *  + GrabExonCoordSetNucls     : For a set of exon coordinates, extract their single (spliced) nucleotide sequence
 *  + GetSplicedExonCoordSets   : Uses 'FindComponentBestEnd' and 'GetExonSetFromEndNode' to determine which paths through the graph we want to report
 *  + RightAlignStr             : Given a string, return a string with a specific length where the input string's characters are right-aligned
 *  + IntToCharArr              : Convert an integer to a string
 *  + GetALSBlockLengths        : Determine the length of strings needed for aesthetically pleasing alignment output
 *  + PrintExon                 : Output a single exon from a spliced exon set
 *  + PrintSplicedAlignment     : Output a collection of spliced exons
 *  + ReportSplicedTopHits      : Manage output
 *  + CheckHitMatchesExonCoords : Make sure that the results of the P7 pipeline match the exon set we fed into it
 *  + ExonSetCleanup            : [SHOULD BE UNNECESSARY] Checks for exons that move backwards throught the model
 *  + RunModelOnExonSets        : Calls 'GetSplicedExonCoordSets', runs the P7 pipeline on the resulting translations, and calls 'ReportSplicedTopHits'
 *
 *
 *
 *  FLAGSHIP: SpliceHits : Calls 'get_target_range_sequence' to acquire a general coding region on the nucleotide target
 *                       : Calls 'build_splice_graph' to build a first pass of the splice graph using the results from translated HMMER
 *                       : Calls 'AddMissingExonsToGraph' to identify and patch any potential holes in the splice graph (missing exons)
 *                       : Calls 'RunModelOnExonSets' to generate the final spliced alignments of sets of exons to the model
 *
 *
 */



// Before we get to the fun stuff, let's just set up some
// bureaucratic stuff to make debugging relatively (hopefully)
// painless
static int ALEX_MODE  = 1; // Print some extra metadata around hits
static int DEBUGGING1 = 0; // Print debugging output?
static int DEBUGGING2 = 0; // Print debugging output?

static time_t INIT_SECONDS;

// Ever want to know what function you're in, and how deep it
// is (roughly)? Well, wonder no more!
int FUNCTION_DEPTH = 0;
void DEBUG_OUT (const char * message, const int func_depth_change) {

  time_t curr_seconds;
  curr_seconds = time(NULL);

  if (func_depth_change > 0) 
    FUNCTION_DEPTH += func_depth_change;
  
  fprintf(stderr,"  SplDebug:");
  int debug_depth;
  for (debug_depth=0; debug_depth<FUNCTION_DEPTH; debug_depth++) 
    fprintf(stderr,"  ");
  fprintf(stderr,"%s Time: %ld\n",message, (curr_seconds - INIT_SECONDS));
  fflush(stderr);
 
  if (func_depth_change < 0) 
    FUNCTION_DEPTH += func_depth_change;

}



//////////////////////////////////////////////////////////////////////////////////////////////
//
//  Struct: TARGET_SET
//
//  Desc. : 
//
typedef struct _target_set {

  int num_target_seqs;

  char    ** TargetSeqNames; // Also (individually) borrowed pointers!
  int64_t  * TargetSeqIdx;
  int64_t  * TargetStarts;
  int64_t  * TargetEnds;
  int      * TargetIsRevcomp;

} TARGET_SET;









//////////////////////////////////////////////////////////////////////////////////////////////
//
//  Struct: EXON_DISPLAY_INFO
//
//  Desc. :
//
typedef struct _exon_display_info {

  FILE * ofp;
  int textw;
  P7_ALIDISPLAY * AD;

  // These are all right-aligned, equal-length strings.
  // The 'NameBlank' is just an empty string of the same length.
  char * HMMName;
  char * TransName;
  char * NuclName;
  char * NameBlank;
  int    name_str_len;

  // How many chars would we need to represent the largest
  // model position / genome coordinate as a string?
  char * CoordBlank;
  int    coord_str_len;

  // The relevant data from the 'ExonCoordSet'
  int exon_set_id;
  int exon_id;
  int Nterm;
  int Cterm;

  int hmm_start;
  int hmm_end;
  int nucl_start;
  int nucl_end;
  int revcomp;

} EXON_DISPLAY_INFO;









/* DEBUGGING FUNCTION: DumpNode  */
//
// NOTE: Because this function depends on a 'best (up/down)stream node'
//       it will implode if run before 'EvaluatePaths'.
//       Yes, there is an obvious way to fix this fact, but I'm trying
//       to steamroll through some debugging as I'm leaving this note.
//
void DumpNode (SPLICE_NODE * Node)
{
  fprintf(stderr,"   |\n");
  fprintf(stderr,"   |       NODE %d\n",Node->node_id);
  fprintf(stderr,"   |     ,---------------------------------------------------\n");
  fprintf(stderr,"   |     |\n");
  fprintf(stderr,"   |     |  Source in P7_TOPHITS  :  Hit %d\n", Node->hit_id);
  fprintf(stderr,"   |     |\n");
  fprintf(stderr,"   |     |  Score of Source Hit      :  %f\n",Node->hit_score);
  fprintf(stderr,"   |     |  Score of Path Up To Node :  %f\n",Node->cumulative_score);
  fprintf(stderr,"   |     |\n");
  if (Node->num_in_edges > 0) {
    fprintf(stderr,"   |     |  Number of Incoming Edges :  %d\n",Node->num_in_edges);
    if (Node->best_in_edge != -1)
      fprintf(stderr,"   |     |  Best Upstream Node       :  Node %d\n",Node->UpstreamNodes[Node->best_in_edge]->node_id);
  } else if (Node->is_n_terminal) {
    fprintf(stderr,"   |     |  * N-TERMINAL NODE\n");
  } else {
    fprintf(stderr,"   |     |  - No Incoming Edges\n");
  }
  fprintf(stderr,"   |     |\n");
  if (Node->num_out_edges > 0) {
    fprintf(stderr,"   |     |  Number of Outgoing Edges :  %d\n",Node->num_out_edges);
  } else if (Node->is_c_terminal) {
    fprintf(stderr,"   |     |  * C-TERMINAL NODE\n");
  } else {
    fprintf(stderr,"   |     |  - No Outgoing Edges\n");
  }
  fprintf(stderr,"   |     |\n");
  fprintf(stderr,"   |     '\n");
  fprintf(stderr,"   |\n");
  fflush(stderr);
}
/* DEBUGGING FUNCTION: DumpGraph */
void DumpGraph(SPLICE_GRAPH * Graph)
{
  fprintf(stderr,"\n\n");
  fprintf(stderr,"     SPLICE GRAPH\n");
  fprintf(stderr,"   +=========================================================+\n");
  fprintf(stderr,"   |\n");
  fprintf(stderr,"   |  Total Number of Nodes      : %d\n",Graph->num_nodes);
  fprintf(stderr,"   |  Total Number of Edges      : %d\n",Graph->num_edges);
  fprintf(stderr,"   |\n");
  fprintf(stderr,"   |  Number of N-terminal nodes : %d\n",Graph->num_n_term);
  fprintf(stderr,"   |  Number of C-terminal nodes : %d\n",Graph->num_c_term);
  fprintf(stderr,"   |\n");
  if (Graph->has_full_path)
    fprintf(stderr,"   |  * This graph has at least one full path through the HMM!\n");
  fprintf(stderr,"   |\n");
  fprintf(stderr,"   |\n");
  int node_id;
  for (node_id=1; node_id<=Graph->num_nodes; node_id++) 
    DumpNode(Graph->Nodes[node_id]);
  fprintf(stderr,"   |\n");
  fprintf(stderr,"   |\n");
  fprintf(stderr,"   +=========================================================+\n");
  fprintf(stderr,"\n\n\n");
  fflush(stderr);
}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function:  TARGET_SET_Destroy
 *
 */
void TARGET_SET_Destroy
(TARGET_SET * TS)
{
  free(TS->TargetSeqNames);
  free(TS->TargetSeqIdx);
  free(TS->TargetStarts);
  free(TS->TargetEnds);
  free(TS->TargetIsRevcomp);
  free(TS);
}
















/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: FloatHighLowSortIndex
 *
 *  Desc. : Generate a sort index for an array of floats (in decreasing value)
 *          using a basic mergesort implementation.
 *
 *  Inputs:  1. Vals     : The floats that we want to produce a sorting for.
 *           2. num_vals : The number of values in the array.
 *
 *  Output:  An array of indices corresponding to a sorting of 'Vals'.
 *
 */
int * FloatHighLowSortIndex
(float * Vals, int num_vals)
{
  if (DEBUGGING1) DEBUG_OUT("Starting 'FloatHighLowSortIndex'",1);

  int * Write = malloc(num_vals * sizeof(int));
  int * Read  = malloc(num_vals * sizeof(int));
  int * Tmp;

  int i;
  for (i=0; i<num_vals; i++) {
    Write[i] = i;
    Read[i]  = i;
  }

  // Of course I'm doing a merge sort
  int writer,left_reader,left_break,right_reader,right_break;
  int ms_block_size = 1;
  while (ms_block_size < num_vals) {

    writer = 0;
    while (writer+ms_block_size < num_vals) {

      left_reader = writer;
      left_break  = writer + ms_block_size;

      right_reader = left_break;
      right_break  = right_reader + ms_block_size;
      if (right_break > num_vals)
        right_break = num_vals;

      while (left_reader < left_break && right_reader < right_break) {
        if (Vals[Read[left_reader]] > Vals[Read[right_reader]]) {
          Write[writer] = Read[left_reader];
          writer++;
          left_reader++;
        }
        else { 
          Write[writer] = Read[right_reader];
          writer++;
          right_reader++;
        }
      }
      while ( left_reader <  left_break) {
        Write[writer] = Read[left_reader];
        writer++;
        left_reader++;
      }
      while (right_reader < right_break) {
        Write[writer] = Read[right_reader];
        writer++;
        right_reader++;
      }

    }

    while (writer < num_vals) {
      Write[writer] = Read[writer];
      writer++;
    }

    // Flip 'em
    Tmp = Read;
    Read = Write;
    Write = Tmp;

    // Again! Again!
    ms_block_size *= 2;

  }

  if (DEBUGGING1) DEBUG_OUT("'FloatHighLowSortIndex' Complete",-1);

  free(Write);
  return Read;

}
/*  Function: FloatLowHighSortIndex
 *
 *  Desc. : Invert the sort index produced by 'FloatHighLowSortIndex'
 *
 */
int * FloatLowHighSortIndex
(float * Vals, int num_vals)
{
  if (DEBUGGING1) DEBUG_OUT("Starting 'FloatLowHighSortIndex'",1);
  int * SortIndex = FloatHighLowSortIndex(Vals,num_vals);
  int i;
  for (i=0; i<num_vals/2; i++) {
    int tmp = SortIndex[i];
    SortIndex[i] = SortIndex[(num_vals-1)-i];
    SortIndex[(num_vals-1)-i] = tmp;
  }
  if (DEBUGGING1) DEBUG_OUT("'FloatLowHighSortIndex' Complete",-1);
  return SortIndex;
}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: get_taget_range_coords
 *
 *  Desc.:  Define the target ranges that exist on a single strand of one target quequence 
 *         
 *
 */
void
get_target_range_coords(TARGET_RANGE_SET *target_set, P7_HIT **hits, int *sort_score_idx, const int num_hits, const int complementarity)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'get_target_range_coords'",1);

  int           i, j;
  int           seed_idx;
  int           seed_seq_from, seed_seq_to;
  int64_t       seed_hmm_from, seed_hmm_to;
  int64_t       curr_range_min, curr_range_max;
  int64_t       new_range_min, new_range_max;
  int64_t       lower_bound, upper_bound;
  int64_t      *prev_range_min = NULL;
  int64_t      *prev_range_max = NULL;
  P7_HIT       *seed_hit;
  TARGET_RANGE *target_range;
  int           status;  

  ESL_ALLOC(prev_range_min, num_hits * sizeof(int64_t));
  ESL_ALLOC(prev_range_max, num_hits * sizeof(int64_t));
  
  int num_ranges = 0;
  for (i = 0; i < num_hits; i++) {
  
    target_range = NULL;

    seed_idx = sort_score_idx[i];  //index of highest scoreing hit from which we will build a target range
    seed_hit = hits[seed_idx];

    if (seed_hit->in_target_range) continue;  // skip any hits already assigned to a target range
    if ( !(seed_hit->flags & p7_IS_REPORTED)) continue; // skip any hits below the e-value threshold
     
    /* We have a seed hit for a potential target range 
     * 1) create traget range and
     * 2) add seed hit                                 */
    Target_Range_Set_CreateNextRange(target_set, &target_range, num_hits); 
    target_range->seqname = seed_hit->name;
    target_range->seqidx = seed_hit->seqidx;
    target_range->complementarity = complementarity;   
    target_range->th->hit[0] = seed_hit;
    target_range->th->N++;
    seed_hit->in_target_range = TRUE;

    /* Extend the target range from the seed hit */
    if ( !complementarity) {   // target is on forward strand
      
      seed_seq_from = seed_hit->dcl[0].iali;
      seed_seq_to   = seed_hit->dcl[0].jali;

      seed_hmm_from = seed_hit->dcl[0].ad->hmmfrom;
      seed_hmm_to   = seed_hit->dcl[0].ad->hmmto;
      
      /* Get bounds set by pervious ranges to prevent overlapping ranges */
      lower_bound = 0;
      upper_bound =  seed_hit->dcl[0].ad->L;
      for(j = 0; j < num_ranges; j++) {
         if(prev_range_max[j] < seed_seq_from && prev_range_max[j] > lower_bound)
           lower_bound = prev_range_max[j];

         if(prev_range_min[j] > seed_seq_to && prev_range_min[j] < upper_bound) 
           upper_bound = prev_range_min[j];        
      }
      
      /* loop through hits - if they are compatiable with seed hit, 
       * add them to target range and keep track of min and max cooords */  
      curr_range_min = seed_seq_from;
      curr_range_max = seed_seq_to;      
      for (j = 0; j < num_hits; j++) {
       
        if (hits[j]->in_target_range) continue;                                                     // skip any hits already assigned to a target range

        if ( !(hits[j]->flags & p7_IS_REPORTED)) continue;                                          // skip any hits below the e-value threshold

        if ( hits[j]->dcl[0].iali < lower_bound || hits[j]->dcl[0].jali > upper_bound) continue;  // skip any hits outside bounds

        /* If hit is downstream of seed in both hmm and seq coord, add it to target range */          
        if( hits[j]->dcl[0].ad->hmmto > seed_hmm_to && hits[j]->dcl[0].jali > seed_seq_to) {

          target_range->th->hit[target_range->th->N] = hits[j];
          target_range->th->N++;
          hits[j]->in_target_range = TRUE;

          if (hits[j]->dcl[0].jali > curr_range_max)
            curr_range_max = hits[j]->dcl[0].jali;     
        }

        /* If hit is upstream of seed in both hmm and seq coord, add it to target range */ 
        if( hits[j]->dcl[0].ad->hmmfrom < seed_hmm_from && hits[j]->dcl[0].iali < seed_seq_from) {
             
            target_range->th->hit[target_range->th->N] = hits[j];
            target_range->th->N++; 
            hits[j]->in_target_range = TRUE;

            if(hits[j]->dcl[0].iali < curr_range_min)
              curr_range_min = hits[j]->dcl[0].iali;                  
        }        
      }  

     /* Add any hits that overlap with current taget range
      * and get new range min and max coords  */
      new_range_min = curr_range_min;
      new_range_max = curr_range_max;
      for (j = 0; j < num_hits; j++) {

        if (hits[j]->in_target_range) continue;  // skip any hits already assigned to a target range
        if ( !(hits[j]->flags & p7_IS_REPORTED)) continue; // skip any hits below the e-value threshold

        if ( (hits[j]->dcl[0].iali >= curr_range_min && hits[j]->dcl[0].iali <= curr_range_max) ||
             (hits[j]->dcl[0].jali <= curr_range_max && hits[j]->dcl[0].jali >= curr_range_min)) {

          target_range->th->hit[target_range->th->N] = hits[j];
          target_range->th->N++;
          hits[j]->in_target_range = TRUE;

          if(hits[j]->dcl[0].iali < new_range_min)
            new_range_min = hits[j]->dcl[0].iali;

          if(hits[j]->dcl[0].jali > new_range_max)
            new_range_max = hits[j]->dcl[0].jali;
      
        }
      }
    }
    else {  // reverse strand

      seed_seq_from = seed_hit->dcl[0].jali;
      seed_seq_to   = seed_hit->dcl[0].iali;

      seed_hmm_from = seed_hit->dcl[0].ad->hmmfrom;
      seed_hmm_to   = seed_hit->dcl[0].ad->hmmto;   
      //printf("seed_hmm_from %d seed_hmm_to %d seed_seq_from %d seed_seq_to %d\n", seed_hmm_from, seed_hmm_to, seed_seq_from, seed_seq_to);
      /* Get bounds set by pervious ranges to prevent overlapping ranges */
      lower_bound = 0;
      upper_bound =  seed_hit->dcl[0].ad->L;
      for(j = 0; j < num_ranges; j++) {
         if(prev_range_max[j] < seed_seq_from && prev_range_max[j] > lower_bound)
           lower_bound = prev_range_max[j];

         if(prev_range_min[j] > seed_seq_to && prev_range_min[j] < upper_bound) 
           upper_bound = prev_range_min[j];        
      }
     //printf("lower %d upper %d\n", lower_bound, upper_bound); 
      /* loop through hits - if they are compatiable with seed hit,
       * add them to target range and keep track of min and max cooords */
      curr_range_min = seed_seq_from;
      curr_range_max = seed_seq_to;
      for (j = 0; j < num_hits; j++) {
        if (hits[j]->in_target_range) continue;                                                     // skip any hits already assigned to a target range
        if ( !(hits[j]->flags & p7_IS_REPORTED)) continue;                                          // skip any hits below the e-value threshold
        if ( hits[j]->dcl[0].jali < lower_bound || hits[j]->dcl[0].iali > upper_bound) continue;  // skip any hits outside bounds

        if( hits[j]->dcl[0].ad->hmmto > seed_hmm_to && hits[j]->dcl[0].jali < seed_seq_from) {
      
          target_range->th->hit[target_range->th->N] = hits[j];
          target_range->th->N++;
          hits[j]->in_target_range = TRUE;
           
         if(hits[j]->dcl[0].jali < curr_range_min)
           curr_range_min = hits[j]->dcl[0].jali;

        }
        
        if( hits[j]->dcl[0].ad->hmmfrom < seed_hmm_from && hits[j]->dcl[0].iali > seed_seq_to) {
       
            target_range->th->hit[target_range->th->N] = hits[j];
            target_range->th->N++; 
            hits[j]->in_target_range = TRUE;
         
           if(hits[j]->dcl[0].iali > curr_range_max)
             curr_range_max = hits[j]->dcl[0].iali;              
        }        
      }
      
      /* Add any hits that overlap with current taget range 
       * and get new range min and max coords  */     
      new_range_min = curr_range_min;
      new_range_max = curr_range_max;
      for (j = 0; j < num_hits; j++) {

        if (hits[j]->in_target_range) continue;  // skip any hits already assigned to a target range
        if ( !(hits[j]->flags & p7_IS_REPORTED)) continue; // skip any hits below the e-value threshold
     
        if ( (hits[j]->dcl[0].jali >= curr_range_min && hits[j]->dcl[0].jali <= curr_range_max) ||
             (hits[j]->dcl[0].iali <= curr_range_max && hits[j]->dcl[0].iali >= curr_range_min)) {
      
          target_range->th->hit[target_range->th->N] = hits[j];
          target_range->th->N++;
          hits[j]->in_target_range = TRUE;

          if(hits[j]->dcl[0].iali < new_range_min)
            new_range_min = hits[j]->dcl[0].iali;

          if(hits[j]->dcl[0].jali > new_range_max)
            new_range_max = hits[j]->dcl[0].jali;
             
        }
      }
      
    }

    target_range->start = new_range_min;
    target_range->end   = new_range_max;
    prev_range_min[num_ranges] = new_range_min;
    prev_range_max[num_ranges] = new_range_max;
    target_range = NULL;
    num_ranges++;
  }

  seed_hit = NULL;
  if(prev_range_min != NULL) free(prev_range_min);
  if(prev_range_max != NULL) free(prev_range_max);
  if (DEBUGGING1) DEBUG_OUT("'get_target_range_coords' Complete",-1);

  return;

  ERROR:
    if(prev_range_min != NULL) free(prev_range_min);
    if(prev_range_max != NULL) free(prev_range_max);
    return;

}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: build_target_ranges
 *
 *  Desc.: Scan P7_TOPHITS to define a set of ranges on the target(s) that 
 *         bound possible splice regions. 
 *
 */
TARGET_RANGE_SET * 
build_target_ranges(const P7_TOPHITS * TopHits)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'build_target_ranges'",1);

  
  int               i;
  int               prev_comp, curr_comp;
  int               num_hits;               // counts number ot hits on each target
  int               prev_target_first_hit;  // keeps tack of which hit is the first for each target
  int              *HitScoreSort;
  float            *HitScores;
  TARGET_RANGE_SET *Targets;
  int               status;

  /* Create Target Range Set */
  Targets = Target_Range_Set_Create();

  /* Get score of first hit */  
  ESL_ALLOC(HitScores, TopHits->N * sizeof(float));
  HitScores[0] = TopHits->hit[0]->dcl[0].envsc;

  /* Get complementarity of the first hit */
  prev_comp = 0;
  if (TopHits->hit[0]->dcl[0].iali > TopHits->hit[0]->dcl[0].jali)
     prev_comp = 1;

  num_hits    = 1;
  prev_target_first_hit = 0;
  
  for(i = 1; i < TopHits->N; i++) {
    
    /* Get complementarity of current hit */
    curr_comp = 0;
    if (TopHits->hit[i]->dcl[0].iali > TopHits->hit[i]->dcl[0].jali) 
      curr_comp = 1;     
  
    /* If the current hit is on a different sequence or strand 
     * than the previous hit then this is a new target. Process
     * the old target and then begin a new one */
    if (TopHits->hit[i]->seqidx != TopHits->hit[i-1]->seqidx || curr_comp != prev_comp) {

      /* Sort scores of all hits on previous target from high to low */
      HitScoreSort = FloatHighLowSortIndex(HitScores,num_hits);   
      
      get_target_range_coords(Targets, &TopHits->hit[prev_target_first_hit], HitScoreSort, num_hits, prev_comp);           
 
      HitScores[0] =  TopHits->hit[i]->dcl[0].envsc;   
      
      prev_target_first_hit = i;
      num_hits = 0;
      free(HitScoreSort); 
    } 
    else {
      /* Get scores off hits on current target */
      HitScores[num_hits] =  TopHits->hit[i]->dcl[0].envsc;
    }
    prev_comp = curr_comp;
    num_hits++;
  }

  // Process last target  
  if(num_hits) {
    HitScoreSort = FloatHighLowSortIndex(HitScores,num_hits);
    get_target_range_coords(Targets, &TopHits->hit[prev_target_first_hit], HitScoreSort, num_hits, prev_comp);
    free(HitScoreSort);

  }
  if (DEBUGGING1) DEBUG_OUT("'build_target_ranges' Complete",-1);

  free(HitScores);

  return Targets;

  ERROR:
     free(HitScores);
     return NULL;    
}












/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: get_target_range_sequence
 *
 *  Desc. :  Given a sequence file and a target range, extract a sub-region of the genome
 *           that coresponds to that target range 
 *
 *  Inputs:  1. genomic_seq_file : An ESL_SQFILE struct representing the input genom(e/ic sequence) file.
 *           2. taget_range      : A TARGET_RANGE struct with seqidx, comlementarity, and coords of the 
                                   sub-sequnce to be extracted
 *
 *
 */
ESL_SQ * get_target_range_sequence
(ESL_SQFILE * genomic_seq_file, TARGET_RANGE *target_range)
{


  int         fetch_err_code;
  ESL_SQ     *target_seq;
  ESL_SQFILE *seq_file;

  if (DEBUGGING1) DEBUG_OUT("Starting 'get_target_range_sequence'",1);

  /* Open seq file and ssi */
  esl_sqfile_Open(genomic_seq_file->filename,genomic_seq_file->format,NULL,&seq_file);
  esl_sqfile_OpenSSI(seq_file,NULL);

  /* Get basic sequence info */
  target_seq = esl_sq_Create();
  esl_sqio_FetchInfo(seq_file,target_range->seqname,target_seq);

  target_seq->start = target_range->start;
  target_seq->end   = target_range->end;
  target_seq->abc   = genomic_seq_file->abc;

  /* Extend target range search area and check that it stays in bound of sequence */ 
  target_seq->start -= MAX_SUB_NUCL_RANGE;
  if (target_seq->start < 1)
    target_seq->start = 1;
  
  target_seq->end += MAX_SUB_NUCL_RANGE;
  if (target_seq->end > target_seq->L)
    target_seq->end = target_seq->L;

  
  /* Fetch target range sequcene */
  fetch_err_code    = esl_sqio_FetchSubseq(seq_file,target_seq->name,target_seq->start,target_seq->end,target_seq);
  esl_sq_SetName(target_seq, target_range->seqname);

  esl_sqfile_Close(seq_file);

  if (fetch_err_code != eslOK) {
    fprintf(stderr,"\n  ERROR: Failed to fetch target subsequence (is there an .ssi index for the sequence file?)\n");
    fprintf(stderr,"         Requested search area: %s:%ld..%ld\n\n",target_seq->name,target_seq->start,target_seq->end);
    exit(1);
  }

  /* If target range is on the reverse strand, reverse the seq */
  if(target_range->complementarity) 
   esl_sq_ReverseComplement(target_seq);
   
  esl_sq_Digitize(target_seq->abc, target_seq);

  if (DEBUGGING1) DEBUG_OUT("'get_target_range_sequence' Complete",-1);

  return target_seq;

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: grab_seq_range
 *
 *  Desc. :  Extract a sub-sequence from our (already a sub-)sequence stored
 *           in a ESL_SQ.
 *
 *           Importantly, the 'start' and 'end' coordinate values are with respect
 *           to the full sequence of which TargetNuclSeq is a part.  This is because
 *           we're expecting to be working with the coordinates listed in an ALIDISPLAY,
 *           which are w.r.t. the full input sequence.
 *
 *  Inputs:  1. TargetNuclSeq : The nucleotide sub-sequence that we want to
 *                              extract a (further sub-)sequence from.
 *           2.         start : The (inclusive) start coordinate of the sub-sequence we
 *                              want to extract.
 *           3.           end : The (inclusive) end coordinate of the sub-sequence we
 *                              want to extract.
 *
 *  Output:  An ESL_DSQ containing the numeric nucleotide codes for the requested
 *           sub-sequence.
 *
 *           TODO: Add catches for out-of-bounds requests (these *shouldn't* happen)
 *
 */
ESL_DSQ * grab_seq_range
(ESL_SQ* target_seq, int range_start, int range_end, int complementarity)
{

  int i;
  int range_len;
  int64_t read_index;
  ESL_DSQ *digital_subseq;
  int status;

  range_len = abs(range_end - range_start) + 1;
  //printf("start %d end %d seq_start %d seq_end %d\n", range_start, range_end, target_seq->start, target_seq->end);
  ESL_ALLOC(digital_subseq, sizeof(ESL_DSQ) * (range_len+2)); 
  digital_subseq[0] = eslDSQ_SENTINEL;
  //printf("%d %d Dig[i] %d\n", range_start, 0, digital_subseq[0]); 
  // Keep in mind that DSQs are [1..n]
  if (!complementarity) {
     read_index = range_start - target_seq->start +1;

    for (i=1; i<=range_len; i++) {

      digital_subseq[i] = target_seq->dsq[read_index];
    //  printf("%d %d Dig[i] %d\n", read_index, i, digital_subseq[i]);
      read_index++;
    }
  } else {
     
     read_index = target_seq->start - range_start + 1;
     for (i=1; i<=range_len; i++) {
      digital_subseq[i] = target_seq->dsq[read_index];
     // printf("%d %d Dig[i] %d\n", read_index, i, digital_subseq[i]);
      read_index++;
    }
  }

  digital_subseq[range_len+1] = eslDSQ_SENTINEL;
 // printf("%d %d Dig[i] %d\n", read_index, range_len+1, digital_subseq[range_len+1]);

  return digital_subseq;

  ERROR:
    free(digital_subseq);
    return NULL;
}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: SelectSpliceOpt
 *
 */
int SelectSpliceOpt
(
  SPLICE_EDGE * Overlap,
  ESL_DSQ * UN,
  ESL_DSQ * DN,
  int model_pos,
  P7_PROFILE * gm,
  ESL_GENCODE * gcode,
  float * splice_score
)
{

  int status;
  int best_opt  = -1;
  *splice_score = EDGE_FAIL_SCORE;


  ESL_DSQ * Codon;
  esl_abc_CreateDsq(Overlap->ntalpha,"AAA",&Codon);

  // 0. |ABCxx|...yy|
  Codon[1] = UN[1];
  Codon[2] = UN[2];
  Codon[3] = UN[3];

  int amino_index = esl_gencode_GetTranslation(gcode,&Codon[1]);
  float opt_score;
  if (amino_index >= 0 && amino_index <= 20) {

    opt_score = p7P_MSC(gm,model_pos,amino_index);

    if (UN[4] == 2 && UN[5] == 3) opt_score += SSSCORE[1];
    else                          opt_score += SSSCORE[0];

    if (DN[4] == 0 && DN[5] == 2) opt_score += SSSCORE[1];
    else                          opt_score += SSSCORE[0];

    // I don't understand how this happens, but it does.
    if (!isinf(opt_score)) {
      *splice_score = opt_score;
      best_opt = 0;
    }
    else if(opt_score == eslINFINITY)
      ESL_XEXCEPTION(eslFAIL, "impossible infinite opt score");

  }

  // 1. |ABxx.|..yyF|
  Codon[3] = DN[5];
  
  amino_index = esl_gencode_GetTranslation(gcode,&Codon[1]);
  if (amino_index >= 0 && amino_index <= 20) {

    opt_score = p7P_MSC(gm,model_pos,amino_index);

    if (UN[3] == 2 && UN[4] == 3) opt_score += SSSCORE[1];
    else                          opt_score += SSSCORE[0];

    if (DN[3] == 0 && DN[4] == 2) opt_score += SSSCORE[1];
    else                          opt_score += SSSCORE[0];

    if (!isinf(opt_score) && opt_score > *splice_score) {
      *splice_score = opt_score;
      best_opt = 1;
    }
    else if(opt_score == eslINFINITY)
      ESL_XEXCEPTION(eslFAIL, "impossible infinite opt score");
  }

  // 2. |Axx..|.yyEF|
  Codon[2] = DN[4];

  amino_index = esl_gencode_GetTranslation(gcode,&Codon[1]);
  if (amino_index >= 0 && amino_index <= 20) {

    opt_score = p7P_MSC(gm,model_pos,amino_index);

    if (UN[2] == 2 && UN[3] == 3) opt_score += SSSCORE[1];
    else                          opt_score += SSSCORE[0];

    if (DN[2] == 0 && DN[3] == 2) opt_score += SSSCORE[1];
    else                          opt_score += SSSCORE[0];

    if (!isinf(opt_score) && opt_score > *splice_score) {
      *splice_score = opt_score;
      best_opt = 2;
    }
    else if(opt_score == eslINFINITY)
      ESL_XEXCEPTION(eslFAIL, "impossible infinite opt score");
  }

  // 3. |xx...|yyDEF|
  Codon[1] = DN[3];

  amino_index = esl_gencode_GetTranslation(gcode,&Codon[1]);
  if (amino_index >= 0 && amino_index <= 20) {

    opt_score = p7P_MSC(gm,model_pos,amino_index);

    if (UN[1] == 2 && UN[2] == 3) opt_score += SSSCORE[1];
    else                          opt_score += SSSCORE[0];

    if (DN[1] == 0 && DN[2] == 2) opt_score += SSSCORE[1];
    else                          opt_score += SSSCORE[0];

    if (!isinf(opt_score) && opt_score > *splice_score) {
      *splice_score = opt_score;
      best_opt = 3;
    }
    else if(opt_score == eslINFINITY)
      ESL_XEXCEPTION(eslFAIL, "impossible infinite opt score");
  }


  free(Codon);

  return best_opt;

  ERROR:
	return 0;
}












/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetContestedUpstreamNucls
 *
 */
void GetContestedUpstreamNucls
(ESL_DSQ * NuclSeq, int read_pos, ESL_DSQ * UN)
{
  int status;
  int write_pos = 1;
  while (write_pos<6) {
	if (NuclSeq[read_pos] == eslDSQ_SENTINEL)
      ESL_XEXCEPTION(eslFAIL, "impossible nucleotide position reached %d in GetContestedUpstreamNucls()", read_pos);
    else  
      UN[write_pos++] = NuclSeq[read_pos];
    read_pos++;
  }
  ERROR:
	return;
}




/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetContestedDownstreamNucls
 *
 */
void GetContestedDownstreamNucls
(ESL_DSQ * NuclSeq, int read_pos, ESL_DSQ * DN)
{
  int status;
  int write_pos = 5;
  //printf("%d%d%d%d%d\n", NuclSeq[read_pos-4], NuclSeq[read_pos-3], NuclSeq[read_pos-2], NuclSeq[read_pos-1], NuclSeq[read_pos]);
  while (write_pos) {
    if (NuclSeq[read_pos] == eslDSQ_SENTINEL) 
      ESL_XEXCEPTION(eslFAIL, "impossible nucleotide position reached %d in GetContestedDownstreamNucls()", read_pos); 
    else
      DN[write_pos--] = NuclSeq[read_pos];
    read_pos--;
  }
  ERROR:
	return;
}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: AminoScoreAtPosition
 *
 */
float AminoScoreAtPosition
(
  P7_PROFILE * gm, 
  int amino_index, 
  int model_pos,
  P7_ALIDISPLAY * AD, // If AD == NULL we're in the extension (all matches!)
  int display_pos,
  int * state
)
{

  float transition_score;
  float emission_score;
  
  if (amino_index == 27) {
    if      (*state == p7P_MM) { *state = p7P_MI; }
    else if (*state == p7P_MI) { *state = p7P_II; }
    else if (*state == p7P_IM) { *state = p7P_MM; }
    else if (*state != p7P_II) { printf("incorrect state %d, should be %d\n", *state, p7P_II); }
    // Otherwise, we *should* be staying II
  } else if (AD && AD->model[display_pos] == '.') {
    if      (*state == p7P_MM) { *state = p7P_MD; }
    else if (*state == p7P_MD) { *state = p7P_DD; }
    else if (*state == p7P_DM) { *state = p7P_MM; }
    else if (*state != p7P_DD) { printf("incorrect state %d, should be %d\n", *state, p7P_DD); }
    // Otherwise, we *should* be staying DD
  } else {
    if      (*state == p7P_MI) { *state = p7P_IM; }
    else if (*state == p7P_II) { *state = p7P_IM; }
    else if (*state == p7P_MD) { *state = p7P_DM; }
    else if (*state == p7P_DD) { *state = p7P_DM; }
    else                       { *state = p7P_MM; }
  }
  
  if (model_pos < gm->M) transition_score = p7P_TSC(gm,model_pos,*state);
  else                   transition_score = 0.;
   
   emission_score = p7P_MSC(gm,model_pos,amino_index);
   if (emission_score == -eslINFINITY) emission_score = 0;

  return transition_score + emission_score;


}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: ali_score_at_postion
 *
 */


float ali_score_at_postion
(
  P7_PROFILE * gm, 
  int amino, 
  int model_pos,
  int trans_pos,
  int prev_state,
  int curr_state
)
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

  return transition_score + emission_score;

  ERROR:
    return -eslINFINITY;

}





/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: find_optimal_splice_site
 *
 *  Desc. :
 *
 *  Inputs:  1.                 Overlap :
 *           2.                      gm : The straightforward profile for the protein / family.
 *           3.                   gcode : An ESL_GENCODE struct (mainly used for translation).
 *           4.   upstream_splice_index :
 *           5. downstream_splice_index :
 *           6. split_amino_model_index :
 *
 *  Output:  A float representing the total score of the optimal splicing of the two hits.
 *
 */
float find_optimal_splice_site
(
  SPLICE_EDGE * Overlap,
  P7_PROFILE     * gm,
  ESL_GENCODE    * gcode,
  int *   upstream_splice_index,
  int * downstream_splice_index,
  int * split_amino_model_index,
  int * codon_split_option
)
{
  int     z1,z2,i;
  int     us_start, ds_start;
  int     us_pos, ds_pos;
  int     us_arr_cap, ds_arr_cap; // malloc sizes
  int     us_trans_len, ds_trans_len; //counted of number of columns in overlap regions
  int     us_pre_ext_end_pos, ds_pre_ext_end_pos;  //positiotn of firt or last column that is part of orignal hit
  int     nucl_read_pos;
  int     model_pos;
  int     amino_index;
  int     prev_state;
  int     next_state;
  int     splice_option; 
  int     optimal_us_pos;  
  int     optimal_ds_pos;   
  int     optimal_model_pos; 
  int     optimal_splice_opt; 

  float   baseline_score;
  float   splice_score;
  float   transition_score;
  float   sum_score;
  float   optimal_score;  

  int    *USModelPos;
  int    *USNuclPos;
  int    *USGaps;
  int    *USStates;
  int    *DSModelPos; 
  int    *DSNuclPos;  
  int    *DSGaps;    
  int    *DSStates;

  float    *USScores;
  float    *DSScores;   

  P7_TRACE *up_trace;
  P7_TRACE *down_trace;

  ESL_DSQ  *UN;
  ESL_DSQ  *DN;

  int status;

  if (DEBUGGING1) DEBUG_OUT("Starting 'find_optimal_splice_site'",1);

  up_trace   = Overlap->UpstreamTrace;
  down_trace = Overlap->DownstreamTrace;
  
   //printf("\nOverlap->amino_start %d Overlap->amino_end %d\n", Overlap->amino_start, Overlap->amino_end);
   //printf("up_trace->hmmfrom[0] %d up_trace->hmmto[0] %d down_trace->hmmfrom[0] %d down_trace->hmmto[0] %d\n", up_trace->hmmfrom[0], up_trace->hmmto[0],  down_trace->hmmfrom[0], down_trace->hmmto[0]);
   //printf("Overlap->downstream_nucl_start %d Overlap->downstream_nucl_end %d\n", Overlap->downstream_nucl_start, Overlap->downstream_nucl_end);
   //printf("Overlap->upstream_nucl_start %d Overlap->upstream_nucl_end %d\n", Overlap->upstream_nucl_start, Overlap->upstream_nucl_end); 

/*************UPSTREAM****************/

  /* Allocation size for upstream overlap data equals the trace 
   * length of the hit overlap plus the extention length  */
  us_arr_cap = (Overlap->upstream_trace_end - Overlap->upstream_trace_start + 1) + (Overlap->amino_end - up_trace->hmmto[0]);
 
  /* Allocate the arrays we'll use to capture the
   * upstream candidates for splice site splitting */
  ESL_ALLOC(USModelPos, us_arr_cap*sizeof(int));
  ESL_ALLOC(USNuclPos,  us_arr_cap*sizeof(int));
  ESL_ALLOC(USStates,   us_arr_cap*sizeof(int));
  ESL_ALLOC(USGaps,     us_arr_cap*sizeof(int));
  ESL_ALLOC(USScores,   us_arr_cap*sizeof(float));


  /* Zeroth position is place holder for transtioting into first postion */
  USScores[0] = 0.;
  USStates[0]  = p7T_X;
  USModelPos[0] = 0;
  USGaps[0] = 0;

  model_pos     = Overlap->amino_start;
  z2            = Overlap->upstream_trace_start;
  prev_state    = p7T_X; //placeholder to signify no transtion into first overlap postion for ali_score_at_postion()
  nucl_read_pos = 1;
  us_trans_len  = 1; 

  while (z2 < Overlap->upstream_trace_end) {

    USNuclPos[us_trans_len]  = nucl_read_pos;
    USModelPos[us_trans_len] = model_pos;
    USStates[us_trans_len]   = up_trace->st[z2];

    if      (up_trace->st[z2] == p7T_M) {
      amino_index = esl_gencode_GetTranslation(gcode,&(Overlap->UpstreamNucls[nucl_read_pos]));
      nucl_read_pos += 3;   
      USGaps[us_trans_len] = 0;
    }
    else if (up_trace->st[z2] == p7T_I) {
      nucl_read_pos += 3;     
      USGaps[us_trans_len] = 0;
    }
    else if (up_trace->st[z2] == p7T_D) 
      USGaps[us_trans_len] = 1;
    
    USScores[us_trans_len] = ali_score_at_postion(gm, amino_index, model_pos, model_pos-1, prev_state, up_trace->st[z2]);

    if (!(up_trace->st[z2] == p7T_I)) 
      model_pos++;

    prev_state = up_trace->st[z2];
    us_trans_len++;
    z2++;
  }
  us_pre_ext_end_pos = us_trans_len-1; // For scoring purposes
  
  // Incorporate the extension
  while (model_pos <= Overlap->amino_end) {
  
    USNuclPos[us_trans_len] = nucl_read_pos;
    USModelPos[us_trans_len] = model_pos;
    USStates[us_trans_len]  = p7T_M;

    amino_index = esl_gencode_GetTranslation(gcode,&(Overlap->UpstreamNucls[nucl_read_pos]));
    nucl_read_pos += 3;
    USGaps[us_trans_len] = 0;

    USScores[us_trans_len]   = ali_score_at_postion(gm, amino_index, model_pos, model_pos-1,prev_state, p7T_M);
  
    prev_state = p7T_M;
    model_pos++;
    us_trans_len++;
  }

 /*************DOWNSTREAM****************/
 
  ds_arr_cap = (Overlap->downstream_trace_end - Overlap->downstream_trace_start + 1) + (down_trace->hmmfrom[0] - Overlap->amino_start + 1);

  // Time to malloc!
  ESL_ALLOC(DSModelPos, ds_arr_cap*sizeof(int));
  ESL_ALLOC(DSNuclPos,  ds_arr_cap*sizeof(int));
  ESL_ALLOC(DSStates,   ds_arr_cap*sizeof(int));
  ESL_ALLOC(DSGaps,     ds_arr_cap*sizeof(int));
  ESL_ALLOC(DSScores,   ds_arr_cap*sizeof(float));

  ds_trans_len = 0;
  nucl_read_pos    = 3; // Because we sneak in 2 extra nucleotides for splice signal
  model_pos        = Overlap->amino_start; // plus one to start atfer the first splice site
  
  while (model_pos < down_trace->hmmfrom[0]) {
 
    amino_index = esl_gencode_GetTranslation(gcode,&(Overlap->DownstreamNucls[nucl_read_pos]));
    nucl_read_pos += 3;
    DSGaps[ds_trans_len] = 0;

    DSNuclPos[ds_trans_len] = nucl_read_pos-1;
    DSModelPos[ds_trans_len] = model_pos; // minus one to recored the location of the splice that would result in the downstream hit startig here
    DSStates[ds_trans_len]    = p7T_M;
    DSScores[ds_trans_len]   = ali_score_at_postion(gm, amino_index, model_pos, model_pos, p7T_M, p7T_M);

    model_pos++;
    ds_trans_len++;
  }

  z1 = Overlap->downstream_trace_start;
  while (z1 < Overlap->downstream_trace_end) {
    
    next_state = down_trace->st[z1+1];
    
    if      (down_trace->st[z1] == p7T_M) {
      amino_index = esl_gencode_GetTranslation(gcode,&(Overlap->DownstreamNucls[nucl_read_pos]));
      nucl_read_pos += 3;
      DSGaps[ds_trans_len] = 0;
    }
    else if (down_trace->st[z1] == p7T_I) {
      nucl_read_pos += 3;
      DSGaps[ds_trans_len] = 0;
    }
    else if (down_trace->st[z1] == p7T_D) 
      DSGaps[ds_trans_len] = 1;
    
    DSNuclPos[ds_trans_len]  = nucl_read_pos-1;
    DSModelPos[ds_trans_len] = model_pos;
    DSStates[ds_trans_len]   = down_trace->st[z1];
    DSScores[ds_trans_len]   = ali_score_at_postion(gm, amino_index, model_pos, model_pos, down_trace->st[z1], next_state);
    
    if (!(down_trace->st[z1] == p7T_I)) 
      model_pos++;
     
    ds_trans_len++;
    z1++;
  }

  /* Final downstream state is at Overlap->amino_end+1 so that a 
   * splice at Overlap->amino_end has somewhere to transition to */
  DSScores[ds_trans_len] = 0.; 
  DSStates[ds_trans_len] = p7T_X;
  DSModelPos[ds_trans_len] = model_pos;


  /* sum scores and gaps */
  for (i=1; i<us_trans_len; i++) {
    USScores[i] += USScores[i-1];
    USGaps[i] += USGaps[i-1];
  }

  for (i=ds_trans_len-2; i>=0; i--) 
    DSScores[i] += DSScores[i+1];

  for (i=1; i<ds_trans_len; i++)
    DSGaps[i] += DSGaps[i-1];
   

  /*Calculate a baseline score of the sum of just the hit 
   * overlap reigons  of the two hits */
  ds_pre_ext_end_pos =  down_trace->hmmfrom[0] - Overlap->amino_start;
  baseline_score     = DSScores[ds_pre_ext_end_pos];
  baseline_score    += USScores[us_pre_ext_end_pos];


  // Arrays we'll use for splice site evaluation
  ESL_ALLOC(UN, 6*sizeof(ESL_DSQ));
  ESL_ALLOC(DN, 6*sizeof(ESL_DSQ));

  optimal_us_pos     = 0;
  optimal_ds_pos     = 0;
  optimal_model_pos  = 0;
  optimal_splice_opt = 0;
  optimal_score      = EDGE_FAIL_SCORE;
  
  us_start = 1;
  ds_start = 0;
  for (model_pos = Overlap->amino_start; model_pos <= Overlap->amino_end; model_pos++) {

    while (USModelPos[us_start] < model_pos) us_start++;
    while (DSModelPos[ds_start] < model_pos) ds_start++;

    us_pos = us_start;
    while (us_pos < us_trans_len && USModelPos[us_pos] == model_pos) {
      ds_pos = ds_start;
 		
      while (ds_pos < ds_trans_len && DSModelPos[ds_pos] == model_pos) {
        // Pull in the "contested" nucleotides for this position
        GetContestedUpstreamNucls(Overlap->UpstreamNucls,USNuclPos[us_pos],UN);
        GetContestedDownstreamNucls(Overlap->DownstreamNucls,DSNuclPos[ds_pos],DN);

        splice_option = SelectSpliceOpt(Overlap,UN,DN,model_pos,gm,gcode,&splice_score);
  
        //Calculation scores for upstream and downtream trnsitions into the M state of the splice site
        if      (model_pos == 1)             transition_score = 0;
        else if (USStates[us_pos-1] == p7T_M) transition_score = p7P_TSC(gm,model_pos-1,p7P_MM);
        else if (USStates[us_pos-1] == p7T_I) transition_score = p7P_TSC(gm,model_pos-1,p7P_IM);
        else if (USStates[us_pos-1] == p7T_D) transition_score = p7P_TSC(gm,model_pos-1,p7P_DM);
        else                                 transition_score = 0.; 

        if      (DSStates[ds_pos+1] == p7T_M) transition_score += p7P_TSC(gm,model_pos,p7P_MM);
        else if (DSStates[ds_pos+1] == p7T_I) transition_score += p7P_TSC(gm,model_pos,p7P_MI);
        else if (DSStates[ds_pos+1] == p7T_D) transition_score += p7P_TSC(gm,model_pos,p7P_MD);
        
        sum_score = USScores[us_pos-1] + DSScores[ds_pos+1] + splice_score + transition_score;
   //     printf("model_pos %d sum_score %f SScores[us_pos-1] %f DSScores[ds_pos+1] %f splice_score %f transition_score %f\n", model_pos, sum_score, USScores[us_pos-1], DSScores[ds_pos+1], splice_score, transition_score);
        
        if (!isinf(sum_score) && sum_score > optimal_score) {
          optimal_score      = sum_score;
          optimal_us_pos     = us_pos-1;
          optimal_ds_pos     = ds_pos;
          optimal_model_pos  = model_pos;
          optimal_splice_opt = splice_option;
        }
        ds_pos++;
      }
      us_pos++;
    }

    if (us_pos >= us_trans_len || ds_pos >= ds_trans_len)
      break;

    us_start = us_pos;
    ds_start = ds_pos;
  }
 


  /* Clean Up */ 
  free(USModelPos);
  free(USNuclPos);
  free(USScores);
  free(USStates);
  
  free(DSModelPos);
  free(DSNuclPos);
  free(DSScores);
  free(DSStates);

  free(UN);
  free(DN);


  if (optimal_score == EDGE_FAIL_SCORE) {
    free(USGaps);
    free(DSGaps);
    if (DEBUGGING1) DEBUG_OUT("'find_optimal_splice_site' Complete (but sad)",-1);
    return EDGE_FAIL_SCORE;
  }


  // Currently, the optimal positions share an amino,
  // which we're going to want to consider methods
  // for splitting.
  //
  // For that reason, what we return as the splice
  // indices (relative to the *nucleotide* sequences)
  // is the last / first "safe" nucelotide (inside the
  // exon).
  //
  *upstream_splice_index   = 3 * (optimal_us_pos - USGaps[optimal_us_pos]);
  *downstream_splice_index = 3 * (optimal_ds_pos - DSGaps[optimal_ds_pos]) + 4;
  *split_amino_model_index = optimal_model_pos;
  *codon_split_option      = optimal_splice_opt;

  free(USGaps);
  free(DSGaps);


  if (DEBUGGING1) DEBUG_OUT("'find_optimal_splice_site' Complete",-1);

  //printf("optimal_score %f baseline_score %f\n", optimal_score, baseline_score);
  return optimal_score-baseline_score;

  ERROR:
   return EDGE_FAIL_SCORE;
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: splice_edge_overlap
 *
 *  Desc. :
 *
 *  Inputs:  1.       Overlap :
 *           2. TargetNuclSeq :
 *           3.            gm : The straightforward profile for the protein / family.
 *           4.         gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:
 *
 */
void splice_edge_overlap
(
  SPLICE_EDGE * Overlap,
  P7_PROFILE     * gm, 
  ESL_GENCODE    * gcode
)
{

  float splice_score;
  if (DEBUGGING1) DEBUG_OUT("Starting 'splice_edge_overlap'",1);


  // These splice indices are the terminal nucleotides *within* the
  // uncontested coding region, relative to the nucleotide sequences
  // in Overlap
  int   upstream_splice_index;
  int downstream_splice_index;
  int split_amino_model_index;
  int codon_split_option;
  
  splice_score = find_optimal_splice_site(Overlap,gm,gcode,&upstream_splice_index,&downstream_splice_index,&split_amino_model_index,&codon_split_option);

  // This would be really bizarre (and worth checking to see if it
  // ever actually happens...)
  if (codon_split_option == -1) {
    Overlap->score         = EDGE_FAIL_SCORE;
    Overlap->score_density = EDGE_FAIL_SCORE;
    if (DEBUGGING1) DEBUG_OUT("'splice_edge_overlap' Complete (BUT WITH TERRIBLE OPTIONS?!)",-1);
    return;
  }

  if (Overlap->upstream_nucl_start < Overlap->upstream_nucl_end) {
    Overlap->upstream_spliced_nucl_end = Overlap->upstream_nucl_start + (upstream_splice_index + (3 - codon_split_option)) - 1;
    Overlap->downstream_spliced_nucl_start = Overlap->downstream_nucl_start + (downstream_splice_index - codon_split_option) - 1;
  } else {
    Overlap->upstream_spliced_nucl_end = Overlap->upstream_nucl_start - (upstream_splice_index + (3 - codon_split_option)) + 1;
    Overlap->downstream_spliced_nucl_start = Overlap->downstream_nucl_start - (downstream_splice_index - codon_split_option) + 1;
  }


  if (codon_split_option < 2) {
    Overlap->upstream_exon_terminus   = split_amino_model_index;
    Overlap->downstream_exon_terminus = split_amino_model_index + 1;
  } else {
    Overlap->upstream_exon_terminus   = split_amino_model_index - 1;
    Overlap->downstream_exon_terminus = split_amino_model_index;    
  }

  //printf("Overlap->upstream_spliced_nucl_end %d Overlap->downstream_spliced_nucl_start %d Overlap->upstream_exon_terminus %d Overlap->downstream_exon_terminus %d\n",  Overlap->upstream_spliced_nucl_end, Overlap->downstream_spliced_nucl_start, Overlap->upstream_exon_terminus, Overlap->downstream_exon_terminus );

  Overlap->score = splice_score;
  Overlap->score_density = splice_score / (1 + Overlap->amino_end - Overlap->amino_start);


  if (DEBUGGING1) DEBUG_OUT("'splice_edge_overlap' Complete",-1);

}




/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: get_nuc_coords_form_amino_overlap
 *
 *  Desc. :
 *
 *  Inputs:  1. Edge :
 *
 *  Output:
 *
 */
void get_nuc_coords_form_amino_overlap 
(SPLICE_EDGE * Edge, int complementarity)
{
;
   int z1,z2;
   int strand;
   int extention_nuc_count;
   int curr_hmm_pos;
   P7_TRACE      * up_trace;
   P7_TRACE      * down_trace;
   
 
  up_trace = Edge->UpstreamTrace;
  down_trace = Edge->DownstreamTrace; 
  
  strand = (complementarity? -1 : 1);

/***************UPSTREAM*******************/

 /* Get number of nucleotides that need to be added to the end of the 
  * upstream edge nucleotide range if the end of the amino range 
  * had to be extended to meet the mimimum overlap */ 
  extention_nuc_count = (Edge->amino_end - up_trace->hmmto[0]) * strand * 3;
  Edge->upstream_nucl_end = up_trace->sqto[0] + extention_nuc_count;  
 

  /* Set Edge->upstream_nucl_start to the begining of the upstream 
   * nulceotides that correpond to overlapping aminos */ 
  Edge->upstream_nucl_start = up_trace->sqto[0] + strand;

/*
  if(up_trace->ndom > 1) {
    for(tr_dom = 0; tr_dom < up_trace->ndom; tr_dom++) {
      for (z2 = up_trace->tfrom[tr_dom]; z2 >= 0; z2--) 
        if (up_trace->st[z2] == p7T_M) break;
 
      if(z2 != -1) break;
    }
  }
  else */
    for (z2 = up_trace->N-1 ; z2 >= 0; z2--) if (up_trace->st[z2] == p7T_M) break; 
  
  curr_hmm_pos = up_trace->hmmto[0];
  Edge->upstream_trace_end = z2+1;
  while(curr_hmm_pos >= Edge->amino_start) {
    if      (up_trace->st[z2] == p7T_M) {
      Edge->upstream_nucl_start -= 3 * strand;
      curr_hmm_pos--;
    } 
    else if (up_trace->st[z2] == p7T_I) 
      Edge->upstream_nucl_start -= 3 * strand;
    else if (up_trace->st[z2] == p7T_D)
      curr_hmm_pos--;
    else
      break;
    z2--; 
 
  }
  Edge->upstream_trace_start = z2+1;
  
/***************DOWNSTREAM*******************/

/* Get number of nucleotides that need to be added to the start of the
  * downstream edge nucleotide range if the end of the amino range
  * had to be extended to meet the mimimum overlap */
  extention_nuc_count = (down_trace->hmmfrom[0] - Edge->amino_start) * strand * 3;
  Edge->downstream_nucl_start = down_trace->sqfrom[0] - extention_nuc_count;
  
  

 /* Set Edge->downstream_nucl_end to the end of the downstream
   * nulceotides that correpond to overlapping aminos */
  Edge->downstream_nucl_end   = down_trace->sqfrom[0] - strand;
/*
  if(down_trace->ndom > 1) {
    printf("TRACE %d\n", down_trace->ndom);
    for(tr_dom = 0; tr_dom < down_trace->ndom; tr_dom++) {
      for (z1 = down_trace->tfrom[tr_dom]; z1 < down_trace->N; z1++)
        if (down_trace->st[z1] == p7T_M) break;

      if(z1 != down_trace->N) break;
    }
  }
  else*/
    for (z1 = 0; z1 < down_trace->N; z1++) if (down_trace->st[z1] == p7T_M) break;

  curr_hmm_pos = down_trace->hmmfrom[0];
  Edge->downstream_trace_start = z1;
  while(curr_hmm_pos <= Edge->amino_end) {

    if      (down_trace->st[z1] == p7T_M) {
      Edge->downstream_nucl_end += 3 * strand;
      curr_hmm_pos++;
    }
    else if (down_trace->st[z1] == p7T_I)
      Edge->downstream_nucl_end += 3 * strand;
    else if (down_trace->st[z1] == p7T_D)
      curr_hmm_pos++;
    else
      break;
    z1++;
  }
  Edge->downstream_trace_end = z1;

  /* I don't know what this is for but I think its allways zero */
  Edge->downstream_ext_len = 0; 
}





/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: sketch_splice_edge
 *
 *  Desc. :
 *
 *  Inputs:  1.          Edge :
 *           2. TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           3.            gm : The straightforward profile for the protein / family.
 *           4.         gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:
 *
 */
void sketch_splice_edge
(
  SPLICE_EDGE * Edge,
  ESL_SQ       *target_seq,
  P7_PROFILE     * gm,
  ESL_GENCODE    * gcode,
  int              complementarity
)
{

  int num_ext_aminos;
  P7_TRACE *up_trace;
  P7_TRACE *down_trace;
  if (DEBUGGING1) DEBUG_OUT("Starting 'sketch_splice_edge'",1);


  up_trace   = Edge->UpstreamTrace;
  down_trace = Edge->DownstreamTrace;

  Edge->amino_start = down_trace->hmmfrom[0];
  Edge->amino_end   = up_trace->hmmto[0];  

  /* splicing requires a minimum hmm postion overlap between upstream
   * and downstram hits. If hits to do not meet this minimum we will
   * extend the amino_start and amino_end postions of the Edge to 
   * create the minimum overlap */
  num_ext_aminos = MIN_AMINO_OVERLAP - (1 + Edge->amino_end - Edge->amino_start);
  if (num_ext_aminos > 0) {
    num_ext_aminos = (num_ext_aminos+1)/2;
    Edge->amino_start -= num_ext_aminos;
    Edge->amino_end   += num_ext_aminos;
  }

  /* Make sure extention did not push past the bounds of the hits*/
  Edge->amino_start = ESL_MAX(Edge->amino_start, up_trace->hmmfrom[0]);
  Edge->amino_end   = ESL_MIN(Edge->amino_end,   down_trace->hmmto[0]);

  /* Get the upstream and downstream sequence coords 
   * to cover the hmm overlap region */
   get_nuc_coords_form_amino_overlap(Edge, complementarity);


  // Grab them nucleotides!
  // Note that in order to ensure we can consider *every*
  // position in the overlap, we need to grab 2 extra
  // nucleotides to serve as splice signal candidates.
  /* Create ESL_DSQs fo the overalp sequence coords found by get_nuc_coords_form_amino_overlap() */
  if (!complementarity) {
    Edge->UpstreamNucls   = grab_seq_range(target_seq, Edge->upstream_nucl_start, Edge->upstream_nucl_end+2, complementarity); 
    Edge->DownstreamNucls = grab_seq_range(target_seq, Edge->downstream_nucl_start-2, Edge->downstream_nucl_end, complementarity);
  } else {
    Edge->UpstreamNucls   = grab_seq_range(target_seq, Edge->upstream_nucl_start, Edge->upstream_nucl_end-2, complementarity); 
    Edge->DownstreamNucls = grab_seq_range(target_seq, Edge->downstream_nucl_start+2, Edge->downstream_nucl_end, complementarity);
  }

  // Finish off by adding this friendly little pointer
  Edge->ntalpha = target_seq->abc;


  // Big ol' DEBUGGING dump
  if (DEBUGGING2 && 1) {
    fprintf(stderr,"\n");
    fprintf(stderr,"  Overlap  Nucl. Range:   Upstream : %d ... %d\n",  Edge->upstream_nucl_start,  Edge->upstream_nucl_end);
    fprintf(stderr,"                                   : ");
    int i;
    for (i=1; i<=abs(Edge->upstream_nucl_start-Edge->upstream_nucl_end)+1; i++)
      fprintf(stderr,"%c",DNA_CHARS[Edge->UpstreamNucls[i]]);
    fprintf(stderr,"\n");
    fprintf(stderr,"                      : Downstream : %d ... %d\n",Edge->downstream_nucl_start,Edge->downstream_nucl_end);
    fprintf(stderr,"                                   : ");
    for (i=3; i<=abs(Edge->downstream_nucl_start-Edge->downstream_nucl_end)+3; i++)
      fprintf(stderr,"%c",DNA_CHARS[Edge->DownstreamNucls[i]]);
    fprintf(stderr,"\n");
    fprintf(stderr,"                      : Model Pos.s: %d..%d\n",Edge->amino_start,Edge->amino_end);
    fprintf(stderr,"\n\n");
  }



  splice_edge_overlap(Edge,gm,gcode);


  if (DEBUGGING1) DEBUG_OUT("'sketch_splice_edge' Complete",-1);

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: HitsAreSpliceComaptible
 *
 *  Desc. :
 *
 *  Inputs:  1.   Upstream :
 *           2. Downstream :
 *
 *  Output:
 *
 */
int HitsAreSpliceCompatible
(P7_ALIDISPLAY * Upstream, P7_ALIDISPLAY * Downstream)
{

  int     up_amino_start, up_amino_end;
  int     down_amino_start, down_amino_end;
  int64_t up_nucl_start, up_nucl_end;
  int64_t down_nucl_start, down_nucl_end; 
  if (DEBUGGING1) DEBUG_OUT("Starting 'HitsAreSpliceCompatible'",1);

    // Start by checking if we either have amino acid
  // overlap, or are close enough to consider extending
  up_amino_start = Upstream->hmmfrom;
  up_amino_end   = Upstream->hmmto;

  down_amino_start = Downstream->hmmfrom;
  down_amino_end   = Downstream->hmmto;

  // If the upstream ain't upstream, then obviously we can't treat
  // these as splice-compatible!
  //
  // NOTE: We'll allow these to overlap fully to accommodate cases
  //       where there's an optimal splice junction in the middle
  //       that fixes overextension
  //
  if (up_amino_start > down_amino_start || up_amino_end > down_amino_end) {
    if (DEBUGGING1) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
    return 0;
  }

  // Do we have overlap OR sufficient proximity to consider
  // extending?
  if (up_amino_end + MAX_AMINO_EXT < down_amino_start) {
    if (DEBUGGING1) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
    return 0;
  }

  // Fantastic!  The amino acid coordinates support splice
  // compatibility!  Now it's just time to confirm that
  // the nucleotides also look good.

  up_nucl_start = Upstream->sqfrom;
  up_nucl_end   = Upstream->sqto;
  
  int revcomp1 = 0;
  if (up_nucl_start > up_nucl_end)
    revcomp1 = 1;


  down_nucl_start = Downstream->sqfrom;
  down_nucl_end   = Downstream->sqto;
 
  int revcomp2 = 0;
  if (down_nucl_start > down_nucl_end)
    revcomp2 = 1;


  if (revcomp1 != revcomp2) {
    if (DEBUGGING1) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
    return 0;
  }


  // We want to make sure that these aren't unrealistically
  // close together on the genome...
  if (revcomp1) {

    if (down_nucl_start + (3 * MAX_AMINO_EXT) >= up_nucl_end) {
      if (DEBUGGING1) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
      return 0;
    }

  } else {

    if (down_nucl_start - (3 * MAX_AMINO_EXT) <= up_nucl_end) {
      if (DEBUGGING1) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
      return 0;
    }

  }


  if (DEBUGGING1) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);


  // Looks like we've got a viable upstream / downstream pair!
  return 1;

}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: are_hits_splice_comaptible
 *
 *  Desc. :
 *
 *  Inputs:  1.   Upstream :
 *           2. Downstream :
 *
 *  Output:
 *
 */
int are_hits_splice_comaptible
(P7_DOMAIN *Upstream, P7_DOMAIN *Downstream, int complementarity)
{

  int     up_amino_start, up_amino_end;
  int     down_amino_start, down_amino_end;
  int64_t up_nucl_end;
  int64_t down_nucl_start; 
  if (DEBUGGING1) DEBUG_OUT("Starting 'are_hits_splice_comaptible'",1);

  up_amino_start = Upstream->tr->hmmfrom[0];
  up_amino_end   = Upstream->tr->hmmto[0];

  down_amino_start = Downstream->tr->hmmfrom[0];
  down_amino_end   = Downstream->tr->hmmto[0];

  up_nucl_end   = Upstream->tr->sqto[0];
 
  down_nucl_start = Downstream->tr->sqfrom[0];
  
  /* If the upstream hit starts before or ends before the downstream 
   * hit on hmm coords they might be splice compatible */
   //printf("up_amino_start %d down_amino_start %d \n", up_amino_start, down_amino_start);
   //printf("up_amino_end %d 
  if (up_amino_start <= down_amino_start && up_amino_end <= down_amino_end) {
     
    /* If the end of the upstream hit hmm end is not less than 
    * MAX_AMINO_EXT from the downstream hit hmm start then
    * they might be splice compatible */  
    if (up_amino_end + MAX_AMINO_EXT >= down_amino_start) { 
    
      /* If the gap between the upstream and downstream hit seq coords 
       * is at lest 3*MAX_AMINO_EXT, the gap could by an intron, and 
       * the hits might be splice compatible */
   
      if (( complementarity && (down_nucl_start + (3 * MAX_AMINO_EXT) < up_nucl_end)) ||
          (!complementarity && (down_nucl_start - (3 * MAX_AMINO_EXT) > up_nucl_end))) {
  
        if (DEBUGGING1) DEBUG_OUT("'are_hits_splice_comaptible' Complete",-1);
        return 1;
      }
    }
  }

 
  if (DEBUGGING1) DEBUG_OUT("'are_hits_splice_comaptible' Complete",-1);
  /* If you've made it here the hits are not splice compatible */
  return 0;

 }








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: ExcessiveGapContent
 *
 */
int ExcessiveGapContent
(P7_ALIDISPLAY * AD)
{

  int num_gaps = 0;
  int contiguous_gap_cutoff = 15;

  int ad_pos;
  int num_contiguous_gaps = 0;
  for (ad_pos=0; ad_pos<AD->N; ad_pos++) {

    if (AD->model[ad_pos] == '.' || AD->aseq[ad_pos] == '-') {

      num_gaps++;
    
      num_contiguous_gaps++;
      if (num_contiguous_gaps >= contiguous_gap_cutoff)
        return 1;

    } else {

      num_contiguous_gaps = 0;

    }

  }


  // 25% cutoff
  if (4*num_gaps > AD->N)
    return 1;


  return 0;


}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: excessive_gap_content
 *
 */
int excessive_gap_content
(P7_TRACE * TR)
{

  int z1,z2;
  int num_gaps;
  int num_contiguous_gaps;

  for (z1 = 0; z1 < TR->N; z1++) if (TR->st[z1] == p7T_M) break;

  num_gaps = 0;
  num_contiguous_gaps = 0;
  for (z2 = z1 ; z2 < TR->N; z2++) {

    if(TR->st[z2] == p7T_D || TR->st[z2] == p7T_I) {
        num_gaps++;
        num_contiguous_gaps++;
        if (num_contiguous_gaps >= 15) 
          return 1;        
    } 
    else
      num_contiguous_gaps = 0;

  }

  
  // 25% cutoff
  if (4*num_gaps > (z2-z1))
    return 1;


  return 0;

}













/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: create_splice_edges
 *
 *  Desc. :
 *
 *  Inputs:  1.          TopHits :
 *           2.    TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           3.               gm : The straightforward profile for the protein / family.
 *           4.            gcode : An ESL_GENCODE struct (mainly used for translation).
 *           5. num_splice_edges :
 *
 *  Output:
 *
 */
SPLICE_EDGE ** create_splice_edges
(
  TARGET_RANGE *target_range,
  ESL_SQ       *target_seq,
  P7_PROFILE  * gm,
  ESL_GENCODE * gcode,
  int * num_splice_edges
)
{

  int num_hits;
  int edge_capacity;
  int num_edges;
  int   upstream_hit_id;
  int downstream_hit_id;
  SPLICE_EDGE * Edge;
  SPLICE_EDGE ** SpliceEdges;
  P7_HIT * UpstreamHit;
  P7_HIT * DownstreamHit;
  P7_TRACE *up_trace;
  P7_TRACE *down_trace;

  if (DEBUGGING1) DEBUG_OUT("Starting 'create_splice_edges'",1);


  num_hits = target_range->th->N;

  edge_capacity = 2 * num_hits;
  SpliceEdges = (SPLICE_EDGE **)malloc(edge_capacity * sizeof(SPLICE_EDGE *));

  num_edges = 0;
  for (upstream_hit_id = 0; upstream_hit_id < num_hits; upstream_hit_id++) {
   
    UpstreamHit  = target_range->th->hit[upstream_hit_id];
    up_trace = UpstreamHit->dcl[0].tr;
  
    if (excessive_gap_content(up_trace)) continue;
   
    // For each hit, gather all of the indices of other hits that
    // could potentially be downstream exons.
    for (downstream_hit_id=0; downstream_hit_id < num_hits; downstream_hit_id++) {

      if (upstream_hit_id == downstream_hit_id)  continue;

      DownstreamHit  = target_range->th->hit[downstream_hit_id];
      down_trace = DownstreamHit->dcl[0].tr;

      if (excessive_gap_content(down_trace)) continue;
     
      if(are_hits_splice_comaptible(&UpstreamHit->dcl[0], &DownstreamHit->dcl[0], target_range->complementarity)) { 

        // MUST WE RESIZE?!
        if (num_edges == edge_capacity) {

          edge_capacity *= 2;

          SPLICE_EDGE ** MoreSpliceEdges = (SPLICE_EDGE **)malloc(edge_capacity * sizeof(SPLICE_EDGE *));
          int edge_id;
          for (edge_id=0; edge_id<num_edges; edge_id++)
            MoreSpliceEdges[edge_id] = SpliceEdges[edge_id];

          free(SpliceEdges);
          SpliceEdges = MoreSpliceEdges;

        }
 
        Edge = (SPLICE_EDGE *)malloc(sizeof(SPLICE_EDGE));

        Edge->upstream_hit_id   = upstream_hit_id;
        Edge->downstream_hit_id = downstream_hit_id;

        Edge->UpstreamTopHits   = target_range->th;
        Edge->UpstreamDisplay   = NULL; //UpstreamDisplay;
        Edge->UpstreamTrace     = UpstreamHit->dcl[0].tr; 
        Edge->DownstreamTopHits = target_range->th;
        Edge->DownstreamDisplay = NULL; //DownstreamDisplay;
        Edge->DownstreamTrace   = DownstreamHit->dcl[0].tr;

        sketch_splice_edge(Edge, target_seq, gm, gcode, target_range->complementarity);
        
        if (Edge->score != EDGE_FAIL_SCORE) {
          SpliceEdges[num_edges] = Edge;
          num_edges++;
          
        }
        else 
          free(Edge);
      }
    }
  }

  if (DEBUGGING1) DEBUG_OUT("'create_splice_edges' Complete",-1);


  *num_splice_edges = num_edges;
  return SpliceEdges;

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: ConnectNodesByEdge
 *
 *  Desc. :
 *
 *  Inputs:  1.  Edge :
 *           2. Graph :
 *
 *  Output:
 *
 */
void ConnectNodesByEdge
(SPLICE_EDGE * Edge, SPLICE_GRAPH * Graph)
{

  int upstream_node_id;
  int downstream_node_id;
  SPLICE_NODE *UpstreamNode;
  SPLICE_NODE * DownstreamNode;

  if (DEBUGGING1) DEBUG_OUT("Starting 'ConnectNodesByEdge'",1);

  /* Get Upstream */
  if (Edge->UpstreamTopHits == Graph->TopHits) {
    upstream_node_id = Graph->TH_HitDomToNodeID[Edge->upstream_hit_id];
  } else {
    upstream_node_id = Graph->MH_HitToNodeID[Edge->upstream_hit_id];
  }
  UpstreamNode = Graph->Nodes[upstream_node_id];

  /* Get Downstream */
  if (Edge->DownstreamTopHits == Graph->TopHits) {
    downstream_node_id = Graph->TH_HitDomToNodeID[Edge->downstream_hit_id];
  } else {
    downstream_node_id = Graph->MH_HitToNodeID[Edge->downstream_hit_id];
  }
  DownstreamNode = Graph->Nodes[downstream_node_id];

  /* Connect Upstream */
  UpstreamNode->OutEdges[UpstreamNode->num_out_edges] = Edge;
  UpstreamNode->DownstreamNodes[UpstreamNode->num_out_edges] = DownstreamNode;
  UpstreamNode->num_out_edges += 1;

  /* Connect Downstream */
  DownstreamNode->InEdges[DownstreamNode->num_in_edges] = Edge;
  DownstreamNode->UpstreamNodes[DownstreamNode->num_in_edges]  = UpstreamNode;
  DownstreamNode->num_in_edges += 1;

  // Resize?
  SpliceNode_Grow(UpstreamNode);
  SpliceNode_Grow(DownstreamNode);
 
  // Node connected... by an edge!
  Graph->num_edges += 1;

  if (DEBUGGING1) DEBUG_OUT("'ConnectNodesByEdge' Complete",-1);
}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GatherCTermNodes
 *
 *  Desc. :
 *
 *  Inputs:  1. Graph :
 *
 *  Output:
 *
 */
void GatherCTermNodes
(SPLICE_GRAPH * Graph)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'GatherCTermNodes'",1);

  // We'll re-count C terminal hits, just to be certain
  int num_c_term = 0;

  int   * CTermIDs    = malloc(Graph->num_nodes * sizeof(int  ));
  float * CTermScores = malloc(Graph->num_nodes * sizeof(float));

  int node_id;
  for (node_id=1; node_id<=Graph->num_nodes; node_id++) {
    if (Graph->Nodes[node_id]->is_c_terminal) {
      CTermIDs[num_c_term] = node_id;
      CTermScores[num_c_term] = Graph->Nodes[node_id]->cumulative_score;
      num_c_term++;
    }
  }

  int * CTermScoreSort = FloatHighLowSortIndex(CTermScores,num_c_term);
  
  if (Graph->CTermNodeIDs)
    free(Graph->CTermNodeIDs);

  Graph->CTermNodeIDs = malloc(num_c_term * sizeof(int));
  Graph->num_c_term   = num_c_term;
  int c_term_index;
  for (c_term_index=0; c_term_index<num_c_term; c_term_index++)
    Graph->CTermNodeIDs[c_term_index] = CTermIDs[CTermScoreSort[c_term_index]];

  free(CTermIDs);
  free(CTermScores);
  free(CTermScoreSort);

  if (DEBUGGING1) DEBUG_OUT("'GatherCTermNodes' Complete",-1);

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: EdgeWouldEraseNode
 *
 *  Desc. :  This is a catch for the *rare* case that using a path through
 *           the graph could result in a node being erased because the
 *           end of the best edge into it is beyond the start of the best
 *           edge out of it.
 *
 *           When not caught, this creates an ExonCoordinate entry where
 *           hmm_from > hmm_to.
 *
 */
int EdgeWouldEraseNode
(SPLICE_NODE * Node3, int n3_in_edge_index)
{

  // Who's the candidate upstream node for Node3?
  SPLICE_NODE    * Node2  = Node3->UpstreamNodes[n3_in_edge_index];
  SPLICE_EDGE * DO_2_3 = Node3->InEdges[n3_in_edge_index];

  // If Node2 doesn't have a best input edge, there's
  // no risk of overrunning it
  if (Node2->best_in_edge == -1)
    return 0;


  // Who's the best upstream node for Node2?
  // Note that this has already been determined due to the
  // recursive nature of 'PullUpCumulativeScore'
  int n2_in_edge_index    = Node2->best_in_edge;
  SPLICE_EDGE * DO_1_2 = Node2->InEdges[n2_in_edge_index];
  
  // Would the edge connecting nodes 1 and 2 place the
  // end of node 1 after the start of node 3 in the edge
  // between nodes 2 and 3?
  // We'll actually be even ruder and say that this exon
  // would need to contribute at least 2 aminos of its own
  // (we lose 1 to needing a strict less than, hence 3) <- terrible, awful communication
  if (DO_1_2->upstream_exon_terminus >= DO_2_3->downstream_exon_terminus - 3) {

    // We'll actually want to sever this connection, for
    // the purposes of being able to cleanly determine
    // connected components
    int edge_idx;
    Node3->num_in_edges -= 1;
      
    if (n3_in_edge_index <= Node3->num_in_edges) {
      /* edges are freed by InEdges pointers so n3_in_edge_index edge must be freed here */
      SPLICE_EDGE_Destroy(Node3->InEdges[n3_in_edge_index]);
      for(edge_idx = n3_in_edge_index; edge_idx < Node3->num_in_edges; edge_idx++){
        Node3->InEdges[edge_idx]       = Node3->InEdges[edge_idx+1];
        Node3->UpstreamNodes[edge_idx] = Node3->UpstreamNodes[edge_idx+1];
      }
    }

    Node2->num_out_edges -= 1;
    int n2_out_edge_index;
    for (n2_out_edge_index = 0; n2_out_edge_index < Node2->num_out_edges; n2_out_edge_index++) {

      if (Node2->DownstreamNodes[n2_out_edge_index] == Node3) {
        Node2->OutEdges[n2_out_edge_index] = Node2->OutEdges[Node2->num_out_edges];
        Node2->DownstreamNodes[n2_out_edge_index] = Node2->DownstreamNodes[Node2->num_out_edges];
        break;

      }

    }

    return 1;

  }


  // Nope! These two look fine!
  return 0;

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: PullUpCumulativeScore
 *
 */
float PullUpCumulativeScore
(SPLICE_NODE * Node)
{

  int          in_edge_index;
  float        in_edge_score;

  if (DEBUGGING1) DEBUG_OUT("Starting 'PullUpCumulativeScore'",1);

  if (Node->cumulative_score != EDGE_FAIL_SCORE) {
    if (DEBUGGING1) DEBUG_OUT("'PullUpCumulativeScore' Complete",-1);
    return Node->cumulative_score;
  }
  

  // Since 'EdgeWouldEraseNode' can remove edges we'll need to
  // use a 'while' loop instead of a 'for' loop
  in_edge_index = 0;
  while (in_edge_index < Node->num_in_edges) {

    in_edge_score = Node->InEdges[in_edge_index]->score;
    in_edge_score += PullUpCumulativeScore(Node->UpstreamNodes[in_edge_index]);

    if (EdgeWouldEraseNode(Node,in_edge_index)) 
      continue;

    if (in_edge_score > Node->cumulative_score) {
      Node->cumulative_score = in_edge_score;
      Node->best_in_edge     = in_edge_index;
    }

    in_edge_index++;
  }

  if (Node->cumulative_score == EDGE_FAIL_SCORE)
    Node->cumulative_score  = Node->hit_score;
  else
    Node->cumulative_score += Node->hit_score;


  if (DEBUGGING1) DEBUG_OUT("'PullUpCumulativeScore' Complete",-1);


  return Node->cumulative_score;

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: EvaluatePaths
 *
 *  Desc. :
 *
 *  Inputs:  1. Graph :
 *
 *  Output:
 *
 */
void EvaluatePaths
(SPLICE_GRAPH * Graph)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'EvaluatePaths'",1);
 

  // We'll reset all of our pathfinding, in case we're
  // running this after adding missing exons to the graph.
  //
  Graph->has_full_path         = 0;
  Graph->best_full_path_start  = 0;
  Graph->best_full_path_end    = 0;
  Graph->best_full_path_length = 0;
  Graph->best_full_path_score  = EDGE_FAIL_SCORE;
  int node_id;
  for (node_id=1; node_id<=Graph->num_nodes; node_id++) {
    Graph->Nodes[node_id]->best_in_edge     = -1;
    Graph->Nodes[node_id]->cumulative_score = EDGE_FAIL_SCORE;
  }


  // For each node without an outgoing edge, recursively
  // draw score up through the graph.
  //
  for (node_id=1; node_id<=Graph->num_nodes; node_id++) {
    if (Graph->Nodes[node_id]->num_out_edges == 0) 
      PullUpCumulativeScore(Graph->Nodes[node_id]);
  }


  // Once we have all of our cumulative scores pulled up,
  // we can collect and sort the C-terminal nodes by
  // their cumulative scores.
  //
  GatherCTermNodes(Graph);


  // We'll count the number of N-terminal nodes (just for fun)
  Graph->num_n_term = 0;
  for (node_id=1; node_id<=Graph->num_nodes; node_id++)
    if (Graph->Nodes[node_id]->is_n_terminal)
      Graph->num_n_term += 1;


  if (DEBUGGING1) DEBUG_OUT("'EvaluatePaths' Complete",-1);

}





/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: fill_out_graph_structure
 *
 *  Desc. :
 *
 *  Inputs:  1.            Graph :
 *           2.      SpliceEdges :
 *           3. num_splice_edges :
 *
 *  Output:
 *
 */
void fill_out_graph_structure
(
  SPLICE_GRAPH    * Graph, 
  SPLICE_EDGE ** SpliceEdges, 
  int num_splice_edges,
  int complementarity
)
{

  int hit_id;
  int node_id;
  int edge_id;
  uint64_t num_hits;
  int status;

  if (DEBUGGING1) DEBUG_OUT("Starting 'fill_out_graph_structure'",1);

  // We'll want this lookup table to be able to go from
  // SPLICE_EDGE content
  num_hits = Graph->TopHits->N;
  
  // Allocate space for the maximum number of domains
  ESL_ALLOC(Graph->Nodes, (num_hits+1)*sizeof(SPLICE_NODE *));
   ESL_ALLOC(Graph->TH_HitDomToNodeID, (num_hits+1) * sizeof(int));
  node_id = 0;
  for (hit_id=0; hit_id<num_hits; hit_id++) {

    Graph->TH_HitDomToNodeID[hit_id] = ++node_id;
    Graph->Nodes[node_id] = SpliceNode_Create(Graph,node_id,hit_id,0);
  }
  
  Graph->num_nodes = node_id;
  Graph->revcomp = complementarity;

  for (edge_id=0; edge_id<num_splice_edges; edge_id++) {
    if (SpliceEdges[edge_id] != NULL) {
      ConnectNodesByEdge(SpliceEdges[edge_id],Graph);
    }
  }

  EvaluatePaths(Graph);

  if (DEBUGGING1) DEBUG_OUT("'fill_out_graph_structure' Complete",-1);

  return;

  ERROR:
	SPLICE_GRAPH_Destroy(Graph);
    return;

}







/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: FindBestFullPath
 *
 *  Desc. :
 *
 *  Inputs:  1. Graph :
 *
 *  Output:
 *
 */
void FindBestFullPath
(SPLICE_GRAPH * Graph)
{
  
  if (DEBUGGING1) DEBUG_OUT("Starting 'FindBestFullPath'",1);

  // NOTE: Because CTermNodeIDs is sorted by cumulative score,
  //       we can break the first time we find a full path.
  SPLICE_NODE * Walker;
  int c_term_index;
  for (c_term_index=0; c_term_index<Graph->num_c_term; c_term_index++) {
    
    Walker = Graph->Nodes[Graph->CTermNodeIDs[c_term_index]];
    P7_HIT *hit = Graph->TopHits->hit[Walker->hit_id];
    
    Graph->best_full_path_length = 1;
    while (Walker->best_in_edge != -1) {
      Walker = Walker->UpstreamNodes[Walker->best_in_edge];
      hit = Graph->TopHits->hit[Walker->hit_id];
   
      Graph->best_full_path_length += 1;
    }

    if (Walker->is_n_terminal) {
      Graph->has_full_path        = 1;
      Graph->best_full_path_start = Walker->node_id;
      Graph->best_full_path_end   = Graph->CTermNodeIDs[c_term_index];
      Graph->best_full_path_score = Graph->Nodes[Graph->best_full_path_end]->cumulative_score;
      if (DEBUGGING1) DEBUG_OUT("'FindBestFullPath' Complete",-1);
      return;
    } else {
      // RESET!
      Graph->best_full_path_length = 0;
    }

  }

  if (DEBUGGING1) DEBUG_OUT("'FindBestFullPath' Complete",-1);

}







/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: build_splice_graph
 *
 *  Desc. :
 *
 *  Inputs:  1.       TopHits :
 *           2. TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           3.            gm : The straightforward profile for the protein / family.
 *           4.            om : The optimized profile (assumed to be built on 'gm').
 *           5.         gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:
 *
 */
SPLICE_GRAPH * build_splice_graph
(
  TARGET_RANGE *target_range,
  ESL_SQ       *target_seq,
  P7_TOPHITS  * TopHits, 
  P7_PROFILE  * gm,
  P7_OPROFILE * om,
  ESL_GENCODE * gcode
)
{

  int num_splice_edges;
  SPLICE_EDGE ** SpliceEdges;
  SPLICE_GRAPH * Graph;  

  if (DEBUGGING1) DEBUG_OUT("Starting 'build_splice_graph'",1);
  
  // We'll just make an unordered list of our splice edges for now
  num_splice_edges = 0;
  SpliceEdges = create_splice_edges(target_range, target_seq,gm,gcode,&num_splice_edges);

  Graph = Splice_Graph_Create();

  Graph->TopHits = target_range->th;
  Graph->Model   = gm;
  Graph->OModel  = om;

  // Build that stinky graph!
  fill_out_graph_structure (Graph, SpliceEdges, num_splice_edges, target_range->complementarity);

  free(SpliceEdges);

  if (DEBUGGING1) DEBUG_OUT("'build_splice_graph' Complete",-1);

  return Graph;

}






/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  DEBUGGING Function: TestSubModel
 *
 */
void TestSubModel
(P7_PROFILE * Model, P7_OPROFILE * OModel, const char * TargetStr)
{

  P7_GMX   * ViterbiMatrix = p7_gmx_Create(Model->M,1024);
  P7_TRACE * Trace         = p7_trace_Create();

  ESL_SQ * TargetSeq = esl_sq_CreateFrom("target",TargetStr,NULL,NULL,NULL);
  esl_sq_Digitize(Model->abc,TargetSeq);

  float viterbi_score;
  p7_GViterbi(TargetSeq->dsq,strlen(TargetStr),Model,ViterbiMatrix,&viterbi_score);
  p7_GTrace(TargetSeq->dsq,strlen(TargetStr),Model,ViterbiMatrix,Trace);
  p7_trace_Index(Trace);


  P7_ALIDISPLAY * AD = p7_alidisplay_Create(Trace,0,OModel,TargetSeq,NULL);

  fprintf(stdout,"\n\n");
  fprintf(stdout," -- Sub-Model Test\n");
  fprintf(stdout,"    Hit to Positions : %d..%d\n",AD->hmmfrom,AD->hmmto);
  fprintf(stdout,"    Viterbi Score    : %f\n",viterbi_score);
  fprintf(stdout,"\n\nTRACE DUMP\n\n");
  p7_trace_Dump(stdout,Trace,Model,TargetSeq->dsq);
  fprintf(stdout,"\n\n(%s)\n\n",Model->consensus);

  // HARD DEBUGGING
  if (0) {
    fprintf(stdout,"--- Program terminating (TestSubModel kill active) ---\n\n");
    exit(21);
  }

  p7_alidisplay_Destroy(AD);
  esl_sq_Destroy(TargetSeq);
  p7_gmx_Destroy(ViterbiMatrix);
  p7_trace_Destroy(Trace);

}







/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: ExtractSubProfile
 *
 *  Desc. :
 *
 *  Inputs:  1.     FullModel :
 *           2. hmm_start_pos :
 *           3.   hmm_end_pos :
 *
 *  Output:
 *
 */
P7_PROFILE * ExtractSubProfile
(
  P7_PROFILE * FullModel, 
  int hmm_start_pos, 
  int hmm_end_pos
)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'ExtractSubProfile'",1);


  //int fullM = FullModel->M;
  int  subM = 1 + hmm_end_pos - hmm_start_pos;
  P7_PROFILE * SubModel = p7_profile_Create(subM,FullModel->abc);


  SubModel->name = malloc(10*sizeof(char));
  strcpy(SubModel->name,"SubModel\0");


  // I'm just going to declare all of these iterators here.
  int tt,nr,i,j;
  int sub_model_pos,residue_id;


  // 1. TRANSITION SCORES
  //
  for (tt = 0; tt < p7P_NTRANS; tt++) { // "transition type"

    for (sub_model_pos = 0; sub_model_pos < subM; sub_model_pos++)
      SubModel->tsc[sub_model_pos * p7P_NTRANS + tt] = FullModel->tsc[(sub_model_pos-1 + hmm_start_pos) * p7P_NTRANS + tt];

  }


  // 2. EMISSION SCORES
  //
  for (nr = 0; nr < p7P_NR; nr++) {
    for (residue_id = 0; residue_id < FullModel->abc->Kp; residue_id++) {

      // Position 0 is a special little baby
      SubModel->rsc[residue_id][nr] = FullModel->rsc[residue_id][nr];
    
      for (sub_model_pos = 1; sub_model_pos <= subM; sub_model_pos++)
        SubModel->rsc[residue_id][p7P_NR * sub_model_pos + nr] = FullModel->rsc[residue_id][p7P_NR * (sub_model_pos-1 + hmm_start_pos) + nr];
    
    }
  }


  // 3. SPECIAL STATES
  //
  for (i=0; i<p7P_NXSTATES; i++) {
    for (j=0; j<p7P_NXTRANS; j++)
      SubModel->xsc[i][j] = FullModel->xsc[i][j];
  }


  // 4. CONSENSUS SEQUENCE
  //
  SubModel->consensus[0] = FullModel->consensus[0];
  for (sub_model_pos = 1; sub_model_pos <= subM; sub_model_pos++)
    SubModel->consensus[sub_model_pos] = FullModel->consensus[sub_model_pos+hmm_start_pos-1];
  SubModel->consensus[subM+1] = 0;

  //if (SubModel->cs) {
    //free(SubModel->cs);
    //SubModel->cs = NULL;
  //}


  // 5. The REST!
  //
  SubModel->mode       = FullModel->mode;
  SubModel->M          = subM;
  SubModel->L          = FullModel->L;
  SubModel->max_length = FullModel->max_length;
  SubModel->nj         = FullModel->nj;
  SubModel->roff       = -1;
  SubModel->eoff       = -1;
  for (i=0; i< p7_NOFFSETS; i++) SubModel->offs[i]    = FullModel->offs[i];
  for (i=0; i< p7_NEVPARAM; i++) SubModel->evparam[i] = FullModel->evparam[i];
  for (i=0; i< p7_NCUTOFFS; i++) SubModel->cutoff[i]  = FullModel->cutoff[i];
  for (i=0; i< p7_MAXABET;  i++) SubModel->compo[i]   = FullModel->compo[i];



  if (DEBUGGING2 && hmm_start_pos == 1 && hmm_end_pos == FullModel->M) {
    fprintf(stderr,"\n  Full Model Copy Validation:  ");
    if (p7_profile_Compare(SubModel,FullModel,0.0) != eslOK) { fprintf(stderr,"Failed\n\n");   }
    else                                                     { fprintf(stderr,"Success!\n\n"); }
    fflush(stderr);
  }


  if (DEBUGGING1) DEBUG_OUT("'ExtractSubProfile' Complete",-1);


  return SubModel;

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: NodesAreDCCCompatible
 *
 *  Desc. :
 *
 *  Inputs:  1.      upstream_hmm_to :
 *           2.     upstream_nucl_to :
 *           3.  downstream_hmm_from :
 *           4. downstream_nucl_from :
 *
 *  Output:
 *
 */
int NodesAreDCCCompatible
(
  SPLICE_GRAPH * Graph,
  SPLICE_NODE  * UpstreamNode,
  SPLICE_NODE  * DownstreamNode
)
{


  // Pull the relevant data from the upstream and downstream nodes
  P7_DOMAIN * USDom     = &(Graph->TopHits->hit[UpstreamNode->hit_id]->dcl[0]);
  int upstream_hmm_to   = USDom->ad->hmmto;
  int upstream_nucl_to  = USDom->ad->sqto;


  P7_DOMAIN * DSDom        = &(Graph->TopHits->hit[DownstreamNode->hit_id]->dcl[0]);
  int downstream_hmm_from  = DSDom->ad->hmmfrom;
  int downstream_nucl_from = DSDom->ad->sqfrom;


  // Is the "upstream" hit even really upstream?
  int hmm_dist = 1 + downstream_hmm_from - upstream_hmm_to;
  if (hmm_dist <= 0)
    return 0;


  // I mean *REALLY* upstream?
  int nucl_dist = 1;
  if (Graph->revcomp) nucl_dist += upstream_nucl_to - downstream_nucl_from;
  else                nucl_dist += downstream_nucl_from - upstream_nucl_to;
  if (nucl_dist <= 0)
    return 0;


  // Hmmm, I suppose you're oriented correctly...
  // But are you close enough?!
  if (hmm_dist <= MAX_SUB_MODEL_RANGE && nucl_dist <= MAX_SUB_NUCL_RANGE)
    return 1;


  // So close, but not quite ready to work together!
  return 0;

}











/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: IdentifyDisConnComponents
 *
 */
int * IdentifyDisConnComponents
(
  SPLICE_GRAPH * Graph,
  int * NoInEdgeNodes,
  int   num_no_in_edge,
  int * NoOutEdgeNodes,
  int   num_no_out_edge
)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'IdentifyDisConnComponents'",1);


  int i,j,x;


  // This could *surely* be optimized, but I think it's going to
  // be reasonably fast to do an all-versus-all sorta thang...
  int dcc_cap = intMax(num_no_out_edge,num_no_in_edge) + 1;
  int * DisConnCompOuts = malloc(dcc_cap * sizeof(int));
  int * DisConnCompIns  = malloc(dcc_cap * sizeof(int));
  for (i=0; i<dcc_cap; i++) {
    DisConnCompOuts[i] = 0;
    DisConnCompIns[i]  = 0;
  }


  int dcc_ids_issued   = 0; // Note that this is more like "num_dcc_ids issued"
  int num_live_dcc_ids = 0;
  int * LiveDCCIDs = malloc(dcc_cap * sizeof(int));
  for (i=0; i<num_no_out_edge; i++) {


    // Grab the next node without any outgoing edges
    SPLICE_NODE * UpstreamNode = Graph->Nodes[NoOutEdgeNodes[i]];

    
    int dcc_id    = dcc_ids_issued+1; // What's our component ID?
    int connected = 0;                // Did we connect to anyone?

    
    // Scan through the nodes without incoming edges and
    // see who's available for friendship
    for (j=0; j<num_no_in_edge; j++) {


      // Has this node already heard the good news?
      if (DisConnCompIns[j] == dcc_id)
        continue;


      SPLICE_NODE * DownstreamNode = Graph->Nodes[NoInEdgeNodes[j]];


      if (!NodesAreDCCCompatible(Graph,UpstreamNode,DownstreamNode))
        continue;


      // Compatibility achieved!
      // These might be redundantly set multiple times, but that's life!
      DisConnCompOuts[i] = dcc_id;
      connected = 1;


      // This DCC ID is live, baby!
      LiveDCCIDs[num_live_dcc_ids++] = dcc_id;

      // Shouldn't ever happen...
      if (num_live_dcc_ids >= dcc_cap) {
        dcc_cap *= 2;
        int * TmpLive = malloc(dcc_cap*sizeof(int));
        int y;
        for (y=0; y<num_live_dcc_ids; y++)
          TmpLive[y] = LiveDCCIDs[y];
        free(LiveDCCIDs);
        LiveDCCIDs = TmpLive;
      }


      // Has this downstream node already made friends?
      // If so, they're our friends now!
      //
      // This is the lazy way to do this comparison and adjustment,
      // but it's 
      //
      if (DisConnCompIns[j]) {

        int dcc_id_to_replace = DisConnCompIns[j];

        for (x=0; x<num_no_out_edge; x++) {
          if (DisConnCompOuts[x] == dcc_id_to_replace)
            DisConnCompOuts[x] = dcc_id;
        }

        for (x=0; x<num_no_in_edge; x++) {
          if (DisConnCompIns[x] == dcc_id_to_replace) {
            DisConnCompIns[x] = dcc_id;
          }
        }

        // The replaced ID is no longer alive :'(
        num_live_dcc_ids--;
        for (x=0; x<num_live_dcc_ids; x++) {
          if (LiveDCCIDs[x] == dcc_id_to_replace) {
            LiveDCCIDs[x] = dcc_id;
          }
        }


      } else {

        // Just the one friend... for now!
        DisConnCompIns[j] = dcc_id;
      
      }

    }


    // Did we register a new dcc_id?
    //
    // NOTE that it's theoretically possible DCCs could be merged,
    //   so we could end up with an overcount of the number of DCCs.
    //
    if (connected)
      dcc_ids_issued++;

  }



  // Now we can go through all of the hits in each of our
  // "disconnected components" and determine a maximal 
  // search area
  int * DCCSearchRegions = malloc((1 + 4 * num_live_dcc_ids) * sizeof(int));


  // Once again, there is plenty of room for optimization,
  // but I really doubt this is going to be a bottleneck
  // in Splash...
  P7_DOMAIN * DomPtr;
  int meta_dcc_id;
  int final_num_dccs = 0; // In case any search regions are too small
  for (meta_dcc_id = 0; meta_dcc_id < num_live_dcc_ids; meta_dcc_id++) {


    int dcc_id = LiveDCCIDs[meta_dcc_id];


    // 
    // 1. What should we set as the bounds for our search region, based
    //    on the nodes without outgoing edges?
    //

    int dcc_hmm_start  = -1;
    int dcc_nucl_start = -1;

    for (i=0; i<num_no_out_edge; i++) {

      if (DisConnCompOuts[i] != dcc_id)
        continue;

      int node_id = NoOutEdgeNodes[i];
      DomPtr = &(Graph->TopHits->hit[Graph->Nodes[node_id]->hit_id]->dcl[0]);


      // NOTE that the START of the search region will be defined by the
      // minimal END position in the model.
      // It looks like these variable names are at odds, but they aren't,
      // I swear!
      if (dcc_hmm_start == -1 || dcc_hmm_start > DomPtr->ad->hmmto)
        dcc_hmm_start = DomPtr->ad->hmmto;


      if (Graph->revcomp) {
        if (dcc_nucl_start == -1 || dcc_nucl_start < DomPtr->ad->sqto)
          dcc_nucl_start = DomPtr->ad->sqto;
      } else {
        if (dcc_nucl_start == -1 || dcc_nucl_start > DomPtr->ad->sqto)
          dcc_nucl_start = DomPtr->ad->sqto;
      }

    }



    // 
    // 2. What should we set as the bounds for our search region, based
    //    on the nodes without incoming edges?
    //

    int dcc_hmm_end  = -1;
    int dcc_nucl_end = -1;

    for (i=0; i<num_no_in_edge; i++) {

      if (DisConnCompIns[i] != dcc_id)
        continue;

      int node_id = NoInEdgeNodes[i];
      DomPtr = &(Graph->TopHits->hit[Graph->Nodes[node_id]->hit_id]->dcl[0]);


      if (dcc_hmm_end == -1 || dcc_hmm_end < DomPtr->ad->hmmfrom)
        dcc_hmm_end = DomPtr->ad->hmmfrom;


      if (Graph->revcomp) {
        if (dcc_nucl_end == -1 || dcc_nucl_end > DomPtr->ad->sqfrom)
          dcc_nucl_end = DomPtr->ad->sqfrom;
      } else {
        if (dcc_nucl_end == -1 || dcc_nucl_end < DomPtr->ad->sqfrom)
          dcc_nucl_end = DomPtr->ad->sqfrom;
      }

    }



    //
    // 3. Log it! (if it's a long enough stretch to even consider...)
    //
    if (abs(dcc_nucl_end - dcc_nucl_start) > 50) {

      DCCSearchRegions[4*final_num_dccs + 1] = dcc_hmm_start;
      DCCSearchRegions[4*final_num_dccs + 2] = dcc_hmm_end;
      DCCSearchRegions[4*final_num_dccs + 3] = dcc_nucl_start;
      DCCSearchRegions[4*final_num_dccs + 4] = dcc_nucl_end;

      final_num_dccs++;
    
    }


  }
  free(DisConnCompOuts);
  free(DisConnCompIns);
  free(LiveDCCIDs);



  if (DEBUGGING1) DEBUG_OUT("'IdentifyDisConnComponents' Complete",-1);


  DCCSearchRegions[0] = final_num_dccs;

  return DCCSearchRegions;

}











/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: IdentifyMidSearchRegions
 *
 */
int * IdentifyMidSearchRegions
(SPLICE_GRAPH * Graph)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'IdentifyMidSearchRegions'",1);


  int   max_mid_searches = Graph->num_nodes; // Just for mem tracking
  int * MidSearchRegions = malloc((1 + 4 * max_mid_searches) * sizeof(int));
  int   num_mid_searches = 0;
  int   upstream_strong_hmmto;
  int   downstream_strong_hmmfrom;
  int   upstream_nucl_end; 
  int   downstream_nucl_start;
  int   ad_pos, contiguous_stars;
  int   us_node_id, ds_index; //, ds_node_id, mid_index;
  int   min_hmm_gap;
  int * USHMMto   = NULL;
  int * DSHMMfrom = NULL;  

  for (us_node_id=1; us_node_id<=Graph->num_nodes; us_node_id++) {

    SPLICE_NODE * USNode = Graph->Nodes[us_node_id];
    
    if(USHMMto   != NULL) { free(USHMMto);   USHMMto = NULL;}  
    if(DSHMMfrom != NULL) { free(DSHMMfrom); DSHMMfrom = NULL;}
    USHMMto   = malloc(USNode->num_out_edges * sizeof(int));    
    DSHMMfrom = malloc(USNode->num_out_edges * sizeof(int));
  
    for (ds_index=0; ds_index<USNode->num_out_edges; ds_index++) {

      SPLICE_EDGE * DO   = USNode->OutEdges[ds_index];
      P7_ALIDISPLAY  * USAD = DO->UpstreamDisplay;
      P7_ALIDISPLAY  * DSAD = DO->DownstreamDisplay;

      upstream_strong_hmmto    = USAD->hmmto;
      downstream_strong_hmmfrom  = DSAD->hmmfrom; 

      ad_pos = USAD->N;
      contiguous_stars = 0;
      while (ad_pos > 0 && contiguous_stars < 2) {

        ad_pos--;

        if (USAD->ppline[ad_pos] == '*') 
          contiguous_stars++;
        else
          contiguous_stars=0;

        if (USAD->model[ad_pos] != '.')
          upstream_strong_hmmto--;

      }


      ad_pos = 0;
      contiguous_stars = 0;
      while (ad_pos < DSAD->N && contiguous_stars < 2) {

        if (DSAD->ppline[ad_pos] == '*')
          contiguous_stars++;
        else
          contiguous_stars=0;

        if (DSAD->model[ad_pos] != '.')
          downstream_strong_hmmfrom++;

        ad_pos++;

      }


      // We don't want to include the contiguous star
      // positions, so let's nip 'em
      upstream_strong_hmmto     += 2;
      downstream_strong_hmmfrom -= 2;


      // Are the hits strong enough not to warrant a little peek
      // at the intermediate area?
      USHMMto[ds_index]   = upstream_strong_hmmto;       
      DSHMMfrom[ds_index] = downstream_strong_hmmfrom;

    }

    if(USNode->num_out_edges > 0) {
      min_hmm_gap = (DSHMMfrom[0] - USHMMto[0]);
      for (ds_index=1; ds_index<USNode->num_out_edges; ds_index++) {
        if(min_hmm_gap > (DSHMMfrom[ds_index] - USHMMto[ds_index]))
           min_hmm_gap = (DSHMMfrom[ds_index] - USHMMto[ds_index]);  
      }

      if(min_hmm_gap <= 6) continue; 
    }

    for (ds_index=0; ds_index<USNode->num_out_edges; ds_index++) {
  
     if (DSHMMfrom[ds_index] - USHMMto[ds_index] <= 6)
        continue;
     
      SPLICE_EDGE * DO   = USNode->OutEdges[ds_index];
      P7_ALIDISPLAY  * USAD = DO->UpstreamDisplay;
      P7_ALIDISPLAY  * DSAD = DO->DownstreamDisplay;
      

      // Before anything else, is there enough space between the hits
      // that we'd be well-served looking at it?
      upstream_nucl_end   = (int)(USAD->sqto);
      downstream_nucl_start = (int)(DSAD->sqfrom);
      
      if (abs(downstream_nucl_start - upstream_nucl_end) < 50)
        continue;

      // Cool!  Why *not* give it a look?
      MidSearchRegions[4*num_mid_searches+1] = USHMMto[ds_index];
      MidSearchRegions[4*num_mid_searches+2] = DSHMMfrom[ds_index];
      MidSearchRegions[4*num_mid_searches+3] = upstream_nucl_end;
      MidSearchRegions[4*num_mid_searches+4] = downstream_nucl_start;
      num_mid_searches++;


      // Not sure how this would ever happen...
      if (num_mid_searches >= max_mid_searches) {
        max_mid_searches *= 2;
        int * TmpMSR = malloc((1 + 4 * max_mid_searches) * sizeof(int));
        int i;
        for (i=1; i<=4*num_mid_searches; i++)
          TmpMSR[i] = MidSearchRegions[i];
        free(MidSearchRegions);
        MidSearchRegions = TmpMSR;
      }


    }
 
  }

 if(DSHMMfrom != NULL) free(DSHMMfrom); 
 if(USHMMto   != NULL) free(USHMMto);   


  if (DEBUGGING1) DEBUG_OUT("'IdentifyMidSearchRegions' Complete",-1);

  

  MidSearchRegions[0] = num_mid_searches;

  return MidSearchRegions;

}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: identify_mid_search_regions
 *
 */
int * identify_mid_search_regions
(SPLICE_GRAPH * Graph)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'identify_mid_search_regions'",1);

  int   z1, z2;
  int   high_pp_cnt;
  int   max_mid_searches = Graph->num_nodes; // Just for mem tracking
  int * MidSearchRegions = malloc((1 + 4 * max_mid_searches) * sizeof(int));
  int   num_mid_searches = 0;
  int   upstream_strong_hmmto;
  int   downstream_strong_hmmfrom;
  int64_t   upstream_nucl_end; 
  int64_t   downstream_nucl_start;
  int   us_node_id, ds_index; //, ds_node_id, mid_index;
  int   min_hmm_gap;
  int * USHMMto   = NULL;
  int * DSHMMfrom = NULL;  
  float pp;
  P7_TRACE * up_trace;
  P7_TRACE * down_trace;

  for (us_node_id=1; us_node_id<=Graph->num_nodes; us_node_id++) {

    SPLICE_NODE * USNode = Graph->Nodes[us_node_id];
    
    if(USHMMto   != NULL) { free(USHMMto);   USHMMto = NULL;}  
    if(DSHMMfrom != NULL) { free(DSHMMfrom); DSHMMfrom = NULL;}
    USHMMto   = malloc(USNode->num_out_edges * sizeof(int));    
    DSHMMfrom = malloc(USNode->num_out_edges * sizeof(int));
  
    for (ds_index=0; ds_index<USNode->num_out_edges; ds_index++) {

      SPLICE_EDGE * DO   = USNode->OutEdges[ds_index];
      up_trace   = DO->UpstreamTrace;
      down_trace = DO->DownstreamTrace;


      upstream_strong_hmmto    = up_trace->hmmto[0]; 
      downstream_strong_hmmfrom  = down_trace->hmmfrom[0]; 

      
      for (z2 = up_trace->N-1 ; z2 >= 0; z2--) if (up_trace->st[z2] == p7T_M) break;
      high_pp_cnt = 0;
      while(z2 >= 0 && high_pp_cnt < 2) {

        pp = up_trace->pp[z2];
 
        if (pp + 0.05 >= 1.0) 
          high_pp_cnt++;
        else
          high_pp_cnt = 0;
      
        if (up_trace->st[z2] != p7T_I)
          upstream_strong_hmmto--;

        z2--;            
                
      }

     for (z1 = 0; z1 < down_trace->N; z1++) if (down_trace->st[z1] == p7T_M) break;
      
     high_pp_cnt = 0;
     while(z1 <down_trace->N && high_pp_cnt < 2) {  

       pp = down_trace->pp[z1];

       if (pp + 0.05 >= 1.0)
          high_pp_cnt++;
        else
          high_pp_cnt = 0;

       if (down_trace->st[z1] != p7T_I)
         downstream_strong_hmmfrom++;

       z1++;
     }
     
      // We don't want to include the high post 
      // prob positions, so let's nip 'em
      upstream_strong_hmmto     += 2;
      downstream_strong_hmmfrom -= 2;


      // Are the hits strong enough not to warrant a little peek
      // at the intermediate area?
      USHMMto[ds_index]   = upstream_strong_hmmto;       
      DSHMMfrom[ds_index] = downstream_strong_hmmfrom;

    }

    if(USNode->num_out_edges > 0) {
      min_hmm_gap = (DSHMMfrom[0] - USHMMto[0]);
      for (ds_index=1; ds_index<USNode->num_out_edges; ds_index++) {
        if(min_hmm_gap > (DSHMMfrom[ds_index] - USHMMto[ds_index]))
           min_hmm_gap = (DSHMMfrom[ds_index] - USHMMto[ds_index]);  
      }

      if(min_hmm_gap <= 6) continue; 
    }

    for (ds_index=0; ds_index<USNode->num_out_edges; ds_index++) {
  
     if (DSHMMfrom[ds_index] - USHMMto[ds_index] <= 6)
        continue;
     
      SPLICE_EDGE * DO   = USNode->OutEdges[ds_index];
      up_trace   = DO->UpstreamTrace;
      down_trace = DO->DownstreamTrace; 

      // Before anything else, is there enough space between the hits
      // that we'd be well-served looking at it?
      upstream_nucl_end   = up_trace->sqto[0];
      downstream_nucl_start = down_trace->sqfrom[0]; 
      
      if (abs(downstream_nucl_start - upstream_nucl_end) < 50)
        continue;

      // Cool!  Why *not* give it a look?
      MidSearchRegions[4*num_mid_searches+1] = USHMMto[ds_index];
      MidSearchRegions[4*num_mid_searches+2] = DSHMMfrom[ds_index];
      MidSearchRegions[4*num_mid_searches+3] = upstream_nucl_end;
      MidSearchRegions[4*num_mid_searches+4] = downstream_nucl_start;
      num_mid_searches++;


      // Not sure how this would ever happen...
      if (num_mid_searches >= max_mid_searches) {
        max_mid_searches *= 2;
        int * TmpMSR = malloc((1 + 4 * max_mid_searches) * sizeof(int));
        int i;
        for (i=1; i<=4*num_mid_searches; i++)
          TmpMSR[i] = MidSearchRegions[i];
        free(MidSearchRegions);
        MidSearchRegions = TmpMSR;
      }


    }
 
  }

 if(DSHMMfrom != NULL) free(DSHMMfrom); 
 if(USHMMto   != NULL) free(USHMMto);   


  if (DEBUGGING1) DEBUG_OUT("'identify_mid_search_regions' Complete",-1);

  

  MidSearchRegions[0] = num_mid_searches;

  return MidSearchRegions;

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: identify_term_search_regions
 *
 */
int * identify_term_search_regions
(
  SPLICE_GRAPH * Graph,
  ESL_SQ       * target_seq,
  int * NoInEdgeNodes,
  int   num_no_in_edge,
  int * NoOutEdgeNodes,
  int   num_no_out_edge
)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'identify_term_search_regions'",1);


  int i;


  // At most, there will be two search regions (N-,C-terminal).
  //
  int * TermSearchRegions = malloc(9 * sizeof(int));
  int   num_term_searches = 0;


  // Start with the N-terminal scan
  P7_DOMAIN * DomPtr;
  if (Graph->num_n_term == 0) {

    int upstreamest_nucl_pos  = -1;
    int upstreamest_model_pos = Graph->Model->M;

    for (i=0; i<num_no_in_edge; i++) {

      int node_id = NoInEdgeNodes[i];

      DomPtr = &(Graph->TopHits->hit[Graph->Nodes[node_id]->hit_id]->dcl[0]);
      if (DomPtr->ad->hmmfrom < upstreamest_model_pos) {
        upstreamest_model_pos = DomPtr->ad->hmmfrom;
        upstreamest_nucl_pos  = DomPtr->ad->sqfrom;
      }

    }

    // Could just use literals, but in case I change
    // something I'll variable-ize this.
    if (upstreamest_model_pos > 5 && upstreamest_model_pos < 50 && upstreamest_nucl_pos != -1) {

      TermSearchRegions[4*num_term_searches+1] = 1;
      TermSearchRegions[4*num_term_searches+2] = upstreamest_model_pos;
      TermSearchRegions[4*num_term_searches+4] = upstreamest_nucl_pos;
      if (Graph->revcomp) {
        TermSearchRegions[4*num_term_searches+3] = intMin(upstreamest_nucl_pos+10000,(int)(target_seq->start));
      } else {
        TermSearchRegions[4*num_term_searches+3] = intMax(upstreamest_nucl_pos-10000,(int)(target_seq->start));
      }

      // If this search region is too small, skip it
      if (abs(TermSearchRegions[4*num_term_searches+4]-TermSearchRegions[4*num_term_searches+3]) > 50
          && TermSearchRegions[4*num_term_searches+2] > 3)
        num_term_searches++;

    }


  }



  // Any interest in a C-terminal search region?
  if (Graph->num_c_term == 0) {

    int downstreamest_nucl_pos  = -1;
    int downstreamest_model_pos =  1;

    for (i=0; i<num_no_out_edge; i++) {

      int node_id = NoOutEdgeNodes[i];

      DomPtr = &(Graph->TopHits->hit[Graph->Nodes[node_id]->hit_id]->dcl[0]);
      if (DomPtr->ad->hmmto > downstreamest_model_pos) {
        downstreamest_model_pos = DomPtr->ad->hmmto;
        downstreamest_nucl_pos  = DomPtr->ad->sqto;
      }

    }

    int c_term_gap_size = Graph->Model->M - downstreamest_model_pos;
    if (c_term_gap_size > 5 && c_term_gap_size < 50 && downstreamest_nucl_pos != -1) {

      TermSearchRegions[4*num_term_searches+1] = downstreamest_model_pos;
      TermSearchRegions[4*num_term_searches+2] = Graph->Model->M;
      TermSearchRegions[4*num_term_searches+3] = downstreamest_nucl_pos;
      
      if (Graph->revcomp) {
        TermSearchRegions[4*num_term_searches+4] =  intMax(downstreamest_nucl_pos-10000,(int)(target_seq->end));
      } else {
        TermSearchRegions[4*num_term_searches+4] = intMin(downstreamest_nucl_pos+10000,(int)(target_seq->end));
      }

      // If this search region is too small, skip it
      if (abs(TermSearchRegions[4*num_term_searches+4] - downstreamest_nucl_pos) > 50)
        num_term_searches++;

    }

  }



  if (DEBUGGING1) DEBUG_OUT("'identify_term_search_regions' Complete",-1);


  TermSearchRegions[0] = num_term_searches;

  return TermSearchRegions;

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: get_bounded_search_regions
 *
 *  Desc. :
 *
 *  Inputs:
 *
 *  Output:
 *
 */
int * get_bounded_search_regions
(SPLICE_GRAPH * Graph, ESL_SQ *target_seq)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'get_bounded_search_regions'",1);


  // Who all doesn't have an incoming / outgoing edge?
  // (not including "genuine" terminal nodes, of course!)
  int * NoOutEdgeNodes = malloc(Graph->num_nodes * sizeof(int));
  int * NoInEdgeNodes  = malloc(Graph->num_nodes * sizeof(int));
  int num_no_out_edge  = 0;
  int num_no_in_edge   = 0;
  int dest_idx, src_idx;
  int node_id;
  int mid_search1, mid_search2;
  int hmm_gap_start1, hmm_gap_end1;
  int seq_gap_start1, seq_gap_end1;
  int hmm_gap_start2, hmm_gap_end2;
  int seq_gap_start2, seq_gap_end2;
  

  for (node_id = 1; node_id <= Graph->num_nodes; node_id++) {

    SPLICE_NODE * Node = Graph->Nodes[node_id];

    if (Node->num_out_edges == 0 && Node->is_c_terminal == 0)
      NoOutEdgeNodes[num_no_out_edge++] = node_id;

    if (Node->num_in_edges == 0 && Node->is_n_terminal == 0)
      NoInEdgeNodes[num_no_in_edge++] = node_id;

  }

  int *  DCCSearchRegions = IdentifyDisConnComponents(Graph,NoInEdgeNodes,num_no_in_edge,NoOutEdgeNodes,num_no_out_edge);
   int *  MidSearchRegions =  identify_mid_search_regions(Graph); 
//IdentifyMidSearchRegions(Graph); // identify_mid_search_regions(Graph);
                           
  int * TermSearchRegions = identify_term_search_regions(Graph,target_seq,NoInEdgeNodes,num_no_in_edge,NoOutEdgeNodes,num_no_out_edge);

  /*Look for overlaps in mid search regions and combine overlapping regions */ 
  int curr_regions = MidSearchRegions[0];
  int prev_regions = curr_regions + 1;
  while(curr_regions < prev_regions) {
    prev_regions = curr_regions; 
    for(mid_search1 = 0; mid_search1 < MidSearchRegions[0]; mid_search1++) {
      hmm_gap_start1 = MidSearchRegions[4*mid_search1+1];
      hmm_gap_end1   = MidSearchRegions[4*mid_search1+2];
      seq_gap_start1 = ESL_MIN(MidSearchRegions[4*mid_search1+3], MidSearchRegions[4*mid_search1+4]);
      seq_gap_end1   = ESL_MAX(MidSearchRegions[4*mid_search1+3], MidSearchRegions[4*mid_search1+4]);
      for(mid_search2 = 0; mid_search2 < MidSearchRegions[0]; mid_search2++) {
        if(mid_search1 == mid_search2) continue;
        hmm_gap_start2 = MidSearchRegions[4*mid_search2+1];
        hmm_gap_end2   = MidSearchRegions[4*mid_search2+2];
        seq_gap_start2 = ESL_MIN(MidSearchRegions[4*mid_search2+3], MidSearchRegions[4*mid_search2+4]);
        seq_gap_end2   = ESL_MAX(MidSearchRegions[4*mid_search2+3], MidSearchRegions[4*mid_search2+4]);
         
        /* If there is an overlap combine the gap reions, remove old gap regions, and add new one */
        int max_hmm_gap = ESL_MAX(hmm_gap_end1, hmm_gap_end2) - ESL_MIN(hmm_gap_start1, hmm_gap_start2);
        int max_seq_gap = ESL_MAX(seq_gap_end1, seq_gap_end2) - ESL_MIN(seq_gap_start1, seq_gap_start2); 
        if(((hmm_gap_end1 - hmm_gap_start1) + (hmm_gap_end2 - hmm_gap_start2) > max_hmm_gap) &&
           ((seq_gap_end1 - seq_gap_start1) + (seq_gap_end2 - seq_gap_start2) > max_seq_gap)) {
             hmm_gap_start2 = ESL_MIN(hmm_gap_start1, hmm_gap_start2);
             hmm_gap_end2   = ESL_MAX(hmm_gap_end1,     hmm_gap_end2);       
             seq_gap_start2 = ESL_MIN(seq_gap_start1, seq_gap_start2);
             seq_gap_end2   = ESL_MAX(seq_gap_end1,     seq_gap_end2);
  
           int * TmpMSR = malloc((1 + 4 * MidSearchRegions[0]) * sizeof(int));
           dest_idx = 0;
           for (src_idx=0; src_idx<MidSearchRegions[0]; src_idx++){
              if(src_idx == mid_search2) { //This is where we will put the new gap region
                TmpMSR[4*dest_idx+1] = hmm_gap_start2;
                TmpMSR[4*dest_idx+2] = hmm_gap_end2;
                TmpMSR[4*dest_idx+3] = (Graph->revcomp) ? seq_gap_end2   : seq_gap_start2;
                TmpMSR[4*dest_idx+4] = (Graph->revcomp) ? seq_gap_start2 : seq_gap_end2;
              }
              else if(src_idx == mid_search1) { //Remove old gap region
                 dest_idx--;            
              }
              else {
                TmpMSR[4*dest_idx+1] = MidSearchRegions[4*src_idx+1];
                TmpMSR[4*dest_idx+2] = MidSearchRegions[4*src_idx+2];
                TmpMSR[4*dest_idx+3] = MidSearchRegions[4*src_idx+3];
                TmpMSR[4*dest_idx+4] = MidSearchRegions[4*src_idx+4];
              }
              dest_idx++;
          }
          free(MidSearchRegions);
          TmpMSR[0] = dest_idx;
          MidSearchRegions = TmpMSR;
          curr_regions = dest_idx;
        }
      }
    }
  }
 
  // I cast thee out!
  free(NoOutEdgeNodes);
  free(NoInEdgeNodes);
  
  int   final_num_searches    = DCCSearchRegions[0] + MidSearchRegions[0] + TermSearchRegions[0];

  int * SearchRegionAggregate = malloc((1 + 4*final_num_searches)*sizeof(int));

  int x = 0;
  int i;
  for (i=0; i<DCCSearchRegions[0]; i++) {
    SearchRegionAggregate[4*x+1] = DCCSearchRegions[4*i+1];
    SearchRegionAggregate[4*x+2] = DCCSearchRegions[4*i+2];
    SearchRegionAggregate[4*x+3] = DCCSearchRegions[4*i+3];
    SearchRegionAggregate[4*x+4] = DCCSearchRegions[4*i+4];
    x++;
  }
  for (i=0; i<MidSearchRegions[0]; i++) {
    SearchRegionAggregate[4*x+1] = MidSearchRegions[4*i+1];
    SearchRegionAggregate[4*x+2] = MidSearchRegions[4*i+2];
    SearchRegionAggregate[4*x+3] = MidSearchRegions[4*i+3];
    SearchRegionAggregate[4*x+4] = MidSearchRegions[4*i+4];
    x++;
  }
  for (i=0; i<TermSearchRegions[0]; i++) {
    SearchRegionAggregate[4*x+1] = TermSearchRegions[4*i+1];
    SearchRegionAggregate[4*x+2] = TermSearchRegions[4*i+2];
    SearchRegionAggregate[4*x+3] = TermSearchRegions[4*i+3];
    SearchRegionAggregate[4*x+4] = TermSearchRegions[4*i+4];
    x++;
  }
  free(DCCSearchRegions);
  free(MidSearchRegions);
  free(TermSearchRegions);


  if (DEBUGGING1) DEBUG_OUT("'get_bounded_search_regions' Complete",-1);


  SearchRegionAggregate[0] = final_num_searches;

  return SearchRegionAggregate;

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: SelectFinalSubHits
 *
 *  Desc. :
 *
 *  Inputs:  1.          SubHitADs :
 *           2.       SubHitScores :
 *           3.       num_sub_hits :
 *           4.        sub_hmm_len :
 *           5.          hmm_start :
 *           6. final_num_sub_hits :
 *
 *  Output:
 *
 */
P7_DOMAIN ** SelectFinalSubHits
(
  P7_ALIDISPLAY ** SubHitADs,
  P7_TRACE ** SubHitTRs,
  float * SubHitScores,
  int     num_sub_hits,
  int     sub_hmm_len,
  int     hmm_start,
  int   * final_num_sub_hits
)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'SelectFinalSubHits'",1);

  // Get a sorting of the hits by their scores
  // NOTE that these are the uncorrected scores
  //  so they aren't useful for external purposes,
  //  but can still rank our hits *internally*
  int * ADSort = FloatHighLowSortIndex(SubHitScores,num_sub_hits);


  // Now that we have all of our hits, let's find which
  // one(s) will give us coverage of the part of the pHMM
  // we need covered.
  int * Coverage = malloc(sub_hmm_len*sizeof(int));
  int model_pos;
  for (model_pos=0; model_pos<sub_hmm_len; model_pos++) 
    Coverage[model_pos] = 0;

  *final_num_sub_hits = 0;
  int sort_id;
  for (sort_id=0; sort_id<num_sub_hits; sort_id++) {

        
    int sub_hit_id = ADSort[sort_id];
    P7_ALIDISPLAY * AD = SubHitADs[sub_hit_id];
    P7_TRACE      * TR = SubHitTRs[sub_hit_id];
    // Do we get any fresh coverage out of this hit?
    int added_coverage = 0;
    for (model_pos = TR->hmmfrom[0] - hmm_start; model_pos <= TR->hmmto[0] - hmm_start; model_pos++) {
      if (Coverage[model_pos] == 0) {
        Coverage[model_pos] = 1;
        added_coverage++;
      }
    }
    
    // No coverage?! Away with you, filth!
   if (added_coverage < MIN_AMINO_OVERLAP) {
      
      p7_alidisplay_Destroy(AD);
      p7_trace_fs_Destroy(TR);
      SubHitADs[sub_hit_id] = NULL;
      SubHitTRs[sub_hit_id] = NULL;
      continue;
    }


    // THERE WE GO!
    *final_num_sub_hits += 1;


  }
  free(Coverage);
  free(ADSort);



  P7_DOMAIN ** FinalSubHits = malloc(*final_num_sub_hits * sizeof(P7_DOMAIN *));
  *final_num_sub_hits = 0;
  int sub_hit_id;
  for (sub_hit_id = 0; sub_hit_id < num_sub_hits; sub_hit_id++) {

    if (SubHitTRs[sub_hit_id] == NULL)
      continue;

    FinalSubHits[*final_num_sub_hits]           = p7_domain_Create_empty();
    FinalSubHits[*final_num_sub_hits]->ad       = SubHitADs[sub_hit_id];
    FinalSubHits[*final_num_sub_hits]->tr       = SubHitTRs[sub_hit_id];
    FinalSubHits[*final_num_sub_hits]->bitscore = SubHitScores[sub_hit_id];

    *final_num_sub_hits += 1;

  }


  if (DEBUGGING1) DEBUG_OUT("'SelectFinalSubHits' Complete",-1);


  return FinalSubHits;

}












/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: AminoToIndex (I'm sorry, this sucks, but I'm doing it)
 *
 */
int AminoToIndex
(char amino)
{

  // Sanity check / uppercase-ification
  if      (amino > 96 && amino < 123) { amino -= 32; }
  else if (amino < 65 || amino > 90 ) { return   27; }

  switch (amino) {
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'D':
    return 2;
  case 'E':
    return 3;
  case 'F':
    return 4;
  case 'G':
    return 5;
  case 'H':
    return 6;
  case 'I':
    return 7;
  case 'K':
    return 8;
  case 'L':
    return 9;
  case 'M':
    return 10;
  case 'N':
    return 11;
  case 'P':
    return 12;
  case 'Q':
    return 13;
  case 'R':
    return 14;
  case 'S':
    return 15;
  case 'T':
    return 16;
  case 'V':
    return 17;
  case 'W':
    return 18;
  case 'Y':
    return 19;
  case 'U':
    return 20;
  case 'O':
    return 21;
  default:
    return 27;
  }  

  return 27;

}




/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: ComputeRoughAliScore
 *
 */
float ComputeRoughAliScore
(P7_ALIDISPLAY * AD, P7_PROFILE * gm)
{

  int model_pos = AD->hmmfrom;
  int state     = p7P_BM;

  float score = 0.0;

  int ad_pos,amino_index;
  for (ad_pos=0; ad_pos<AD->N; ad_pos++) {

    amino_index = AminoToIndex(AD->aseq[ad_pos]);

    score += AminoScoreAtPosition(gm,amino_index,model_pos,AD,ad_pos,&state);

    if (AD->model[ad_pos] != '.')
      model_pos++;

  }

  return score;

}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: compute_rough_score_from_trace
 *
 */
float compute_rough_score_from_trace
(P7_ALIDISPLAY * AD, P7_PROFILE * gm)
{

  int model_pos = AD->hmmfrom;
  int state     = p7P_BM;

  float score = 0.0;

  int ad_pos,amino_index;
  for (ad_pos=0; ad_pos<AD->N; ad_pos++) {

    amino_index = AminoToIndex(AD->aseq[ad_pos]);

    score += AminoScoreAtPosition(gm,amino_index,model_pos,AD,ad_pos,&state);

    if (AD->model[ad_pos] != '.')
      model_pos++;

  }

  return score;

}





static ESL_OPTIONS options[] = {
  /* name     type         default env_var  range toggles req  incompat help                  docgroup */
 { "--crick", eslARG_NONE, FALSE,  NULL,    NULL, NULL,   NULL, NULL,   "only translate top strand",        99 },
 { "--watson",eslARG_NONE, FALSE,  NULL,    NULL, NULL,   NULL, NULL,   "only translate bottom strand",     99 },
 { "-l",      eslARG_INT,  "5",    NULL,    NULL, NULL,   NULL, NULL,   "minimum ORF length",               99 },
 { "-m",      eslARG_NONE, FALSE,  NULL,    NULL, NULL,   NULL,"-M",    "ORFs must initiate with AUG only", 99 },
 { "-M",      eslARG_NONE, FALSE,  NULL,    NULL, NULL,   NULL,"-m",    "ORFs must start with allowed initiation codon", 99 },

 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* create a set of ORFs for each DNA target sequence */
static int
translate_sequences(ESL_GENCODE *gcode, ESL_GENCODE_WORKSTATE *wrk, ESL_SQ *sq)
{
      esl_gencode_ProcessStart(gcode, wrk, sq);
      esl_gencode_ProcessPiece(gcode, wrk, sq);
      esl_gencode_ProcessEnd(wrk, sq);

  return eslOK;
}


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: find_sub_hits
 *
 *  Desc. :
 *
 *  Inputs:  1.              Graph :
 *           2.      TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           3.       SearchRegion :
 *           4.              gcode : An ESL_GENCODE struct (mainly used for translation).
 *           5. final_num_sub_hits :
 *
 *  Output:
 *
 *  NOTE: This could easily be improved by someone with a keener
 *        understanding of the HMMER internals, but for the purposes
 *        of plugging holes in a splice graph I think this is
 *        adequate... maybe...
 *
 */
P7_DOMAIN ** find_sub_hits
(
  SPLICE_GRAPH * Graph,
  ESL_SQ       * target_seq,
  int          * SearchRegion,
  ESL_GENCODE  * gcode,
  int          * final_num_sub_hits
)
{

   int i;
   int z1,z2;
   int hmm_start, hmm_end;
   int nucl_start, nucl_end;
   int sub_hmm_len;
   int nucl_seq_len;
   int num_sub_hits;
   int max_sub_hits;
   int sub_hit_id;
   int block_size;
   float viterbi_score;
   float nullsc;
   float seq_score;
   double P;
   float * SubHitScores;
   float * NewSubHitScores; 
   char *spoofline;
   P7_PROFILE  *  SubModel;
   P7_OPROFILE * OSubModel;
   ESL_DSQ * SubNucls;
   ESL_SQ * sub_seq;
   ESL_SQ   * ORFAminoSeq;
   P7_BG * bg;
   P7_GMX   * ViterbiMatrix;
   P7_TRACE * Trace;
   P7_ALIDISPLAY ** SubHitADs; 
   P7_TRACE      ** SubHitTRs;
   ESL_GETOPTS *go;  
   ESL_GENCODE_WORKSTATE *wrk;

 
  int submodel_create_err;
  int vit_err_code;
  int gtrace_err_code;


  if (DEBUGGING1) DEBUG_OUT("Starting 'find_sub_hits'",1);


  hmm_start  = SearchRegion[0];
  hmm_end    = SearchRegion[1];
  nucl_start = SearchRegion[2];
  nucl_end   = SearchRegion[3];
 
  // We'll pull in a couple extra HMM positions, since
  // we're interested in having overlaps anyways
  hmm_start = intMax(hmm_start-2,1);
  hmm_end   = intMin(hmm_end  +2,Graph->Model->M);

  sub_hmm_len = 1 + hmm_end - hmm_start;


  // Generate a sub-model and an optimized sub-model
  // for the part of the pHMM we need to fill in
  SubModel = ExtractSubProfile(Graph->Model,hmm_start,hmm_end);
  OSubModel = p7_oprofile_Create(SubModel->M,SubModel->abc);
  submodel_create_err = p7_oprofile_Convert(SubModel,OSubModel);
  
  if (submodel_create_err != eslOK) {
     if (DEBUGGING1) DEBUG_OUT("'find_sub_hits' Failed (could not covert submodel)",-1);
     return NULL;
  }

  // Required for alidisplay generation
  free(OSubModel->name);
  OSubModel->name = malloc(5*sizeof(char));
  strcpy(OSubModel->name,"OSM\0");


  
  // DEBUGGING
  //
  // Note that this is primarily intended to work when there's
  // a particular (tester-known) sequence that we want to confirm
  // would score well under the model.
  //
  //TestSubModel(SubModel,OSubModel,"PPNPSLMSIFRK\0");



  /* Grab the nucleotides we're searching our 
   * sub-model against and create sub sequence */
  SubNucls = grab_seq_range(target_seq,nucl_start,nucl_end,Graph->revcomp);
  nucl_seq_len   = abs(nucl_end-nucl_start)+1;
  sub_seq = esl_sq_CreateDigitalFrom(target_seq->abc, target_seq->name, SubNucls, nucl_seq_len, NULL, NULL, NULL); 

  /* Create spoof getopts to hold flags checked by translation routines */
  spoofline = "getopts --crick -l 5";
  go = esl_getopts_Create(options);
  if (esl_opt_ProcessSpoof(go, spoofline) != eslOK) esl_fatal("errmsg");

  /* Create gencode workstate for use in translting orfs */
  wrk = esl_gencode_WorkstateCreate(go, gcode);
  block_size = (int)(nucl_seq_len/5);  
  wrk->orf_block = esl_sq_CreateDigitalBlock(block_size, Graph->Model->abc);

  /* Translate sub sequnce into set of ORFs */
  translate_sequences(gcode, wrk, sub_seq);

  // Create a whole mess of objects that we'll need to get our
  // our sub-model alignments to fit into an alidisplay...
  ORFAminoSeq   = NULL;
  ViterbiMatrix = p7_gmx_Create(SubModel->M,1024);
  Trace =  p7_trace_Create();    
  bg = p7_bg_fs_Create(Graph->Model->abc);

  // Finally, the array where we'll record all of our incredible triumphs!
  num_sub_hits =  0;
  max_sub_hits = 10;
  SubHitADs = malloc(max_sub_hits * sizeof(P7_ALIDISPLAY *));
  SubHitTRs = malloc(max_sub_hits * sizeof(P7_TRACE *));
  SubHitScores = malloc(max_sub_hits * sizeof(float));

   //printf("wrk->orf_block->count %d\n", wrk->orf_block->count);
   for (i = 0; i < wrk->orf_block->count; i++)
   {

     ORFAminoSeq = &(wrk->orf_block->list[i]);
     /* get null sc */       
     p7_bg_SetLength(bg, ORFAminoSeq->n);
     p7_bg_NullOne  (bg, ORFAminoSeq->dsq, ORFAminoSeq->n, &nullsc); 
     
     /* Configure model */
     p7_ReconfigUnihit(SubModel, ORFAminoSeq->n); 

     p7_gmx_GrowTo(ViterbiMatrix, SubModel->M, ORFAminoSeq->n);
     /* Run Viterbit */
     vit_err_code  = p7_GViterbi(ORFAminoSeq->dsq,ORFAminoSeq->n,SubModel,ViterbiMatrix,&viterbi_score);
     if (vit_err_code != eslOK) {
       fprintf(stderr,"\n  ERROR (find_sub_hits): Failed while running Viterbi\n\n");
     }
     /* Check if Viterbi score passes a threshold */  
     seq_score = (viterbi_score - nullsc) /  eslCONST_LOG2;
     
     P  = esl_gumbel_surv(seq_score, SubModel->evparam[p7_VMU], SubModel->evparam[p7_VLAMBDA]); 
     if(P > 0.1) continue;            


     /* Get the alignment trace */
     gtrace_err_code  = p7_GTrace(ORFAminoSeq->dsq,ORFAminoSeq->n,SubModel,ViterbiMatrix,Trace);
     if (gtrace_err_code != eslOK) {
       fprintf(stderr,"\n  ERROR (find_sub_hits): Failed while generating a generic P7_TRACE for the Viterbi matrix\n\n");
     }


     for (z1 = 0; z1 < Trace->N; z1++) if (Trace->st[z1] == p7T_M) break;
     for (z2 = Trace->N-1 ; z2 >= 0; z2--) if (Trace->st[z2] == p7T_M) break;
     
     if(z2 < z1) { p7_trace_Reuse(Trace); continue; }

     Trace->hmmfrom[0] = Trace->k[z1];
     Trace->hmmto[0]   = Trace->k[z2];  
     Trace->sqfrom[0]  = Trace->i[z1]; 
     Trace->sqto[0]    = Trace->i[z2];

     Trace->hmmfrom[0] += hmm_start - 1;
     Trace->hmmto[0]   += hmm_start - 1;

     if (Graph->revcomp) {
       Trace->sqfrom[0] = nucl_start - ORFAminoSeq->start - ((Trace->sqfrom[0]-2)*3) + 1;
       Trace->sqto[0]   = nucl_start - ORFAminoSeq->start - ((Trace->sqto[0])*3) + 2;
     } else {
       Trace->sqfrom[0] = ((Trace->sqfrom[0]-3)*3) + nucl_start + ORFAminoSeq->start - 1;
       Trace->sqto[0]   = ((Trace->sqto[0])*3) + nucl_start + ORFAminoSeq->start - 2;
     }
     
   
     // Welcome to the list, fella!
     SubHitADs[num_sub_hits] = NULL; //AD;
     
     SubHitScores[num_sub_hits] = seq_score;
     SubHitTRs[num_sub_hits] = p7_trace_fs_Clone(Trace);
     num_sub_hits++;

     // DEBUGGING
     //fprintf(stdout,"PASSED! %d..%d (%f): %s\n",AD->hmmfrom,AD->hmmto,viterbi_score,AD->aseq);
     //p7_trace_Dump(stdout, Trace, SubModel, ORFAminoSeq->dsq);


     // Resize?
     if (num_sub_hits == max_sub_hits) {
        max_sub_hits *= 2;

        P7_ALIDISPLAY ** NewSubHitADs = malloc(max_sub_hits*sizeof(P7_ALIDISPLAY *));
        P7_TRACE      ** NewSubHitTRs = malloc(max_sub_hits*sizeof(P7_TRACE *));
        NewSubHitScores = malloc(max_sub_hits*sizeof(float));
              
        for (sub_hit_id=0; sub_hit_id<num_sub_hits; sub_hit_id++) {
          NewSubHitADs[sub_hit_id]    = SubHitADs[sub_hit_id];
          NewSubHitTRs[sub_hit_id]    = SubHitTRs[sub_hit_id];
          NewSubHitScores[sub_hit_id] = SubHitScores[sub_hit_id];
        }
        free(SubHitADs);    SubHitADs    = NewSubHitADs;
        free(SubHitTRs);    SubHitTRs    = NewSubHitTRs;
        free(SubHitScores); SubHitScores = NewSubHitScores;
           
     }
     // Resize over -- back to our regularly-scheduled programming!




     p7_trace_Reuse(Trace);
     p7_gmx_Reuse(ViterbiMatrix);

  }
  // It takes a lot of fun to need this much cleanup ;)

  esl_sq_Destroy(sub_seq);

  esl_getopts_Destroy(go);
  if(wrk->orf_block != NULL) {
    esl_sq_DestroyBlock(wrk->orf_block);
    wrk->orf_block = NULL;
  } 
  esl_gencode_WorkstateDestroy(wrk);
  p7_bg_Destroy(bg);

  p7_profile_Destroy(SubModel);
  p7_oprofile_Destroy(OSubModel);
  p7_gmx_Destroy(ViterbiMatrix);
  p7_trace_Destroy(Trace);


  // Quick check: Did we find *anything* worth considering?
  if (num_sub_hits == 0) {
    free(SubHitADs);
    free(SubHitTRs);
    free(SubHitScores);
    free(SubNucls);
    if (DEBUGGING1) DEBUG_OUT("'find_sub_hits' Complete (albeit, no hits)",-1);
    *final_num_sub_hits = 0;
    return NULL;
  }


  // Reduce down to just the hits we're really excited about
  
  P7_DOMAIN ** FinalSubHits = SelectFinalSubHits(SubHitADs,SubHitTRs,SubHitScores,num_sub_hits,sub_hmm_len,hmm_start,final_num_sub_hits);
  free(SubHitADs);
  free(SubHitTRs);
  free(SubHitScores);
  free(SubNucls);

  if (DEBUGGING1) DEBUG_OUT("'find_sub_hits' Complete",-1);


  // Happy days!
  return FinalSubHits;

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: IntegrateMissedHits
 *
 *  Desc. :  The goal of this function is to create a new array Graph->Nodes
 *           that now includes any of the little hits that we found to patch
 *           in holes in the original graph.
 *
 *  Inputs:  1.          Graph :
 *           2. NewSpliceEdges :
 *           3.  num_new_edges :
 *
 *  Output:
 *
 */
void IntegrateMissedHits
(SPLICE_GRAPH * Graph, SPLICE_EDGE ** NewSpliceEdges, int num_new_edges)
{
  if (DEBUGGING1) DEBUG_OUT("Starting 'IntegrateMissedHits'",1);


  // First off, we need to copy over all of the existing nodes
  int new_num_nodes = Graph->num_nodes + Graph->MissedHits->N;

  SPLICE_NODE ** NewNodeArray = (SPLICE_NODE **)malloc((new_num_nodes+1)*sizeof(SPLICE_NODE *));

  int node_id;  
  for (node_id = 1; node_id <= Graph->num_nodes; node_id++)
    NewNodeArray[node_id] = Graph->Nodes[node_id];
  for (node_id = Graph->num_nodes+1; node_id <= new_num_nodes; node_id++)
    NewNodeArray[node_id] = NULL;

  free(Graph->Nodes);
  Graph->Nodes = NewNodeArray;


  // Now we can actually integrate the new hits!
  node_id = Graph->num_nodes;
  Graph->num_nodes = new_num_nodes;

  Graph->MH_HitToNodeID = malloc(Graph->MissedHits->N * sizeof(int));
  int missed_hit_id;
  for (missed_hit_id = 0; missed_hit_id < Graph->MissedHits->N; missed_hit_id++) {

    Graph->MH_HitToNodeID[missed_hit_id] = ++node_id;

    Graph->Nodes[node_id] = SpliceNode_Create(Graph,node_id,missed_hit_id,1);

  }


  int new_edge_id;
  for (new_edge_id=0; new_edge_id<num_new_edges; new_edge_id++) {
    if (NewSpliceEdges[new_edge_id] != NULL) {
      ConnectNodesByEdge(NewSpliceEdges[new_edge_id],Graph);
    }
  }

  EvaluatePaths(Graph);


  if (DEBUGGING1) DEBUG_OUT("'IntegrateMissedHits' Complete",-1);

}







/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: seek_missing_exons
 *
 *  Desc. :
 *
 *  Inputs:  1.         Graph :
 *           2. TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           3.         gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:
 *
 */
P7_TOPHITS * seek_missing_exons
(
  SPLICE_GRAPH * Graph, 
  ESL_SQ       * target_seq,
  ESL_GENCODE  * gcode
)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'seek_missing_exons'",1);


  // Grab the sub-regions of our conceptual DP zone that we want
  // to search for missing exons.
  //
  // This 'SearchRegionAggregate' is indexed in groups consisting
  // of: the starting position in the model, the ending position in
  // the model, the starting coordinate on the genome, and the ending
  // coordinate on the genome.
  //
  int * SearchRegionAggregate = get_bounded_search_regions(Graph,target_seq);
  int   num_search_regions    = SearchRegionAggregate[0];

  if (SearchRegionAggregate == NULL) {
    if (DEBUGGING1) DEBUG_OUT("'seek_missing_exons' Complete (none found)",-1);
    return NULL;
  }


  if (DEBUGGING2) {
    fprintf(stderr,"\n  Num Search Regions: %d\n",num_search_regions);
    int i;
    for (i=0; i<num_search_regions; i++) {
      fprintf(stderr,"    - Search Region : %d\n",i+1);
      fprintf(stderr,"      HMM Positions : %d..%d\n",SearchRegionAggregate[i*4 + 1],SearchRegionAggregate[i*4 + 2]);
      fprintf(stderr,"      Nucl. Coord.s : %d..%d\n",SearchRegionAggregate[i*4 + 3],SearchRegionAggregate[i*4 + 4]);
    }
    fprintf(stderr,"\n");
  }


  // Now we can iterate over our list of search regions and,
  // for each:
  //
  //   1. Extract the appropriate sub-region of the model
  //   2. Pull the appropriate portion of the target sequence
  //   3. Search the sub-model against the sub-target,
  //        using the Viterbi algorithm
  //
  //   and...
  //
  //   4. Integrate any new hits into the graph!
  //
  int num_sub_hits = 0;
  int sub_hits_capacity = 20;
  P7_DOMAIN ** SubHits = malloc(sub_hits_capacity*sizeof(P7_DOMAIN *));
  int search_region_id,sub_hit_id,new_sub_hit_id; // iterators
  
  for (search_region_id = 0; search_region_id < num_search_regions; search_region_id++) {

    int num_new_sub_hits = 0;
    P7_DOMAIN ** NewSubHits = find_sub_hits(Graph,target_seq,&SearchRegionAggregate[(4*search_region_id)+1],gcode,&num_new_sub_hits);
   
    if (num_new_sub_hits == 0)
      continue;


    // New sub-hit(s) alert!
    if (DEBUGGING2) {
      fprintf(stderr,"%d additional sub-model hits discovered\n",num_new_sub_hits);
      int i;
      for (i=0; i<num_new_sub_hits; i++) {
        fprintf(stderr,"+ Sub-Hit %d\n",i+1);
        fprintf(stderr,"  Model Range: %d..%d\n",NewSubHits[i]->ad->hmmfrom,NewSubHits[i]->ad->hmmto);
        fprintf(stderr,"  Nucl. Range: %d..%d\n",(int)(NewSubHits[i]->ad->sqfrom),(int)(NewSubHits[i]->ad->sqto));
      }
    }


    // Resize?
    if (num_sub_hits + num_new_sub_hits > sub_hits_capacity) {
      sub_hits_capacity = intMax(sub_hits_capacity*2,num_sub_hits+num_new_sub_hits);
      P7_DOMAIN ** MoreSubHits = malloc(sub_hits_capacity * sizeof(P7_DOMAIN *));
      for (sub_hit_id = 0; sub_hit_id < num_sub_hits; sub_hit_id++)
        MoreSubHits[sub_hit_id] = SubHits[sub_hit_id];
      free(SubHits);
      SubHits = MoreSubHits;
    }


    // Pop dem shrimps on the barbie, Guv! 
    // (as they say way up there in Argentina)
    for (new_sub_hit_id = 0; new_sub_hit_id < num_new_sub_hits; new_sub_hit_id++)
      SubHits[num_sub_hits++] = NewSubHits[new_sub_hit_id];

    
    free(NewSubHits); // clear the pointer

  }
  free(SearchRegionAggregate);


  if (num_sub_hits == 0) {
    free(SubHits);
    return NULL;
  }


  //
  //  NOTE:  From here, what we do is aggregate our new 'SubHits' into
  //         a 'P7_TOPHITS' datastructure.
  //
  //         IMPORTANTLY, we're faking these in order to take advantage
  //         of some of the functions that are already available.
  //         DO NOT ASSUME *ANY* HMMER FUNCTIONS WILL WORK WITH THIS DATASTRUCTURE!
  //
  P7_TOPHITS * MissingHits = malloc(sizeof(P7_TOPHITS));
  MissingHits->N   = (uint64_t)num_sub_hits;
  MissingHits->hit = malloc(num_sub_hits * sizeof(P7_HIT *));

  for (sub_hit_id=0; sub_hit_id<num_sub_hits; sub_hit_id++) {
  
    MissingHits->hit[sub_hit_id]        = malloc(sizeof(P7_HIT));
    MissingHits->hit[sub_hit_id]->name  = target_seq->name;
    MissingHits->hit[sub_hit_id]->ndom  = 1;
    MissingHits->hit[sub_hit_id]->dcl   = SubHits[sub_hit_id];
    MissingHits->hit[sub_hit_id]->score = SubHits[sub_hit_id]->bitscore;

  }
  free(SubHits);
 
  if (DEBUGGING1) DEBUG_OUT("'seek_missing_exons' Complete",-1);


  return MissingHits;

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: add_missing_exons
 *
 *  Desc. :
 *
 *  Inputs:  1.         Graph :
 *           2. TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           3.         gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:
 *
 */
void add_missing_exons
(
  SPLICE_GRAPH * Graph,
  ESL_SQ       * target_seq,
  ESL_GENCODE  * gcode
)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'add_missing_exons'",1);

  
  // As I note in 'SeekMissingExons' (but is worth emphasizing),
  // this datastructure is basically a cheap way to pass around
  // the new hits that we find.
  //
  // MUCH OF THE STANDARD P7_HIT, P7_DOMAIN, and P7_ALIDISPLAY DATA
  // IS NOT CONTAINED IN THESE!  BE CAREFUL BEFORE USING HMMER INTERNAL
  // FUNCTIONS TO INTERROGATE THEM!
  //
  Graph->MissedHits = seek_missing_exons(Graph,target_seq,gcode);

  if (Graph->MissedHits == NULL) {
    if (DEBUGGING1) DEBUG_OUT("'add_missing_exons' Complete (no missed hits added)",-1);
    return;
  }


  // This is basically the same code from 'create_splice_edges,'
  // but reworked to integrate new exons into the splice graph.
  //
  // It probably wouldn't be hard to convert the two to a single function,
  // but the annoying thing is that there's the whole dereferencing thing
  // to get *real* P7_DOMAINs from *actual* P7_HITs (but not my goofy ones)
  //
  int new_splice_edge_cap = 4 * Graph->MissedHits->N;
  SPLICE_EDGE ** NewSpliceEdges = (SPLICE_EDGE **)malloc(new_splice_edge_cap*sizeof(SPLICE_EDGE *));

  int num_missing_hits = (int)(Graph->MissedHits->N);
  int num_new_edges    = 0;
  int missed_hit_id,node_id;
printf("num_missing_hits %d/n", num_missing_hits); 
  for (missed_hit_id = 0; missed_hit_id < num_missing_hits; missed_hit_id++) {

    // Lucky us!  Cheap and easy access!
    P7_ALIDISPLAY * MissedAD = Graph->MissedHits->hit[missed_hit_id]->dcl->ad;
    P7_TRACE * missed_trace = Graph->MissedHits->hit[missed_hit_id]->dcl->tr;
    P7_DOMAIN *missed_domain = Graph->MissedHits->hit[missed_hit_id]->dcl;


    for (node_id = 1; node_id <= Graph->num_nodes; node_id++) {

      int node_hit_id = Graph->Nodes[node_id]->hit_id;

      P7_ALIDISPLAY * NodeAD = (&Graph->TopHits->hit[node_hit_id]->dcl[0])->ad;
      P7_TRACE *node_trace = (&Graph->TopHits->hit[node_hit_id]->dcl[0])->tr;
      P7_DOMAIN *node_domain = Graph->TopHits->hit[node_hit_id]->dcl;
      if (excessive_gap_content(node_trace)) continue;

      // Because order matters for 'HitsAreSpliceCompatible' we
      // need to have a catch for either possibility.
      //if (HitsAreSpliceCompatible(NodeAD,MissedAD)) {
      if (are_hits_splice_comaptible(node_domain, missed_domain, Graph->revcomp)) { 
        NewSpliceEdges[num_new_edges] = (SPLICE_EDGE *)malloc(sizeof(SPLICE_EDGE));

        SPLICE_EDGE * Edge         = NewSpliceEdges[num_new_edges];

        Edge->upstream_hit_id   = node_hit_id;
        Edge->downstream_hit_id = missed_hit_id;

        Edge->UpstreamTopHits   = Graph->TopHits;
        Edge->UpstreamDisplay   = NodeAD;
        Edge->UpstreamTrace     = node_trace; 
        Edge->DownstreamTopHits = Graph->MissedHits;
        Edge->DownstreamDisplay = MissedAD;
        Edge->DownstreamTrace   = missed_trace; 

        num_new_edges++;

    //  } else if (HitsAreSpliceCompatible(MissedAD,NodeAD)) {
      } else if (are_hits_splice_comaptible(missed_domain, node_domain, Graph->revcomp)) {  

         NewSpliceEdges[num_new_edges] = (SPLICE_EDGE *)malloc(sizeof(SPLICE_EDGE));
        SPLICE_EDGE * Edge         = NewSpliceEdges[num_new_edges];

        Edge->upstream_hit_id   = missed_hit_id;
        Edge->downstream_hit_id = node_hit_id;

        Edge->UpstreamTopHits   = Graph->MissedHits;
        Edge->UpstreamDisplay   = MissedAD;
        Edge->UpstreamTrace     = missed_trace;
        Edge->DownstreamTopHits = Graph->TopHits;
        Edge->DownstreamDisplay = NodeAD;
        Edge->DownstreamTrace   = node_trace;

        num_new_edges++;
      }


      // Time to resize?
      if (num_new_edges == new_splice_edge_cap) {
        new_splice_edge_cap *= 2;
       
        SPLICE_EDGE ** MoreNewEdges = (SPLICE_EDGE **)malloc(new_splice_edge_cap*sizeof(SPLICE_EDGE *));
        int i;
        for (i=0; i<num_new_edges; i++)
          MoreNewEdges[i] = NewSpliceEdges[i];

        free(NewSpliceEdges);
        NewSpliceEdges = MoreNewEdges;

      }

    }

    // It's possible we might want to splice two missed exons together
    // (see A1BG human isoform 1)
    int mhi2; // "missed_hit_id_2"
    for (mhi2=missed_hit_id+1; mhi2<num_missing_hits; mhi2++) {

      P7_ALIDISPLAY * MAD2 = Graph->MissedHits->hit[mhi2]->dcl->ad;
      P7_TRACE *missed_trace2 = Graph->MissedHits->hit[mhi2]->dcl->tr;
      P7_DOMAIN *missed_domain2 = Graph->MissedHits->hit[mhi2]->dcl;

      //if (HitsAreSpliceCompatible(MissedAD,MAD2)) {
     if (are_hits_splice_comaptible(missed_domain, missed_domain2, Graph->revcomp)) {   
        NewSpliceEdges[num_new_edges] = (SPLICE_EDGE *)malloc(sizeof(SPLICE_EDGE));
        SPLICE_EDGE * Edge         = NewSpliceEdges[num_new_edges];

        Edge->upstream_hit_id   = missed_hit_id;
        Edge->downstream_hit_id = mhi2;

        Edge->UpstreamTopHits   = Graph->MissedHits;
        Edge->UpstreamDisplay   = MissedAD;
        Edge->UpstreamTrace     = missed_trace; 
        Edge->DownstreamTopHits = Graph->MissedHits;
        Edge->DownstreamDisplay = MAD2;
        Edge->DownstreamTrace  = missed_trace2; 

        num_new_edges++;

     // } else if (HitsAreSpliceCompatible(MAD2,MissedAD)) {
        
        } else if (are_hits_splice_comaptible(missed_domain2, missed_domain, Graph->revcomp)) {
        NewSpliceEdges[num_new_edges] = (SPLICE_EDGE *)malloc(sizeof(SPLICE_EDGE));
        SPLICE_EDGE * Edge         = NewSpliceEdges[num_new_edges];

        Edge->upstream_hit_id   = mhi2;
        Edge->downstream_hit_id = missed_hit_id;

        Edge->UpstreamTopHits   = Graph->MissedHits;
        Edge->UpstreamDisplay   = MAD2;
        Edge->UpstreamTrace     = missed_trace2;
        Edge->DownstreamTopHits = Graph->MissedHits;
        Edge->DownstreamDisplay = MissedAD;
        Edge->DownstreamTrace   = missed_trace; 

        num_new_edges++;

      }

      // Time to resize?
      if (num_new_edges == new_splice_edge_cap) {

        new_splice_edge_cap *= 2;

        SPLICE_EDGE ** MoreNewEdges = (SPLICE_EDGE **)malloc(new_splice_edge_cap*sizeof(SPLICE_EDGE *));
        int i;
        for (i=0; i<num_new_edges; i++)
          MoreNewEdges[i] = NewSpliceEdges[i];

        free(NewSpliceEdges);
        NewSpliceEdges = MoreNewEdges;

      }

    }

  }

  // Back at it again!
  int splice_edge_id;
  for (splice_edge_id = 0; splice_edge_id < num_new_edges; splice_edge_id++) {

    sketch_splice_edge(NewSpliceEdges[splice_edge_id], target_seq,Graph->Model,gcode, Graph->revcomp);


    if (NewSpliceEdges[splice_edge_id]->score == EDGE_FAIL_SCORE) {
      free(NewSpliceEdges[splice_edge_id]);
      NewSpliceEdges[splice_edge_id] = NULL;
    }

  }

  //
  //  If I'm not mistaken.... IT'S PARTY TIME!!!!
  //
  if (num_new_edges)
    IntegrateMissedHits(Graph,NewSpliceEdges,num_new_edges);

  free(NewSpliceEdges);

  if (DEBUGGING1) DEBUG_OUT("'add_missing_exons' Complete",-1);


}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetExonSetFromEndNode
 *
 *  Desc. :
 *
 *  Inputs:  1.         Graph :
 *           2. start_node_id :
 *
 *  Output:
 *
 */
int * GetExonSetFromEndNode
(SPLICE_GRAPH * Graph, int end_node_id)
{
  
  if (DEBUGGING1) DEBUG_OUT("Starting 'GetExonSetFromEndNode'",1);


  // In case we use *every* node
  // BestPathCoords[0] is the number of nodes along the
  // path (Min:1,Max:num_nodes)
  //
  // Each "entry" consists of
  //    [1.] The start position in the genome
  //    [2.] The start position in the profile
  //    [3.] The   end position in the genome
  //    [4.] The   end position in the profile
  //    [5.] The corresponding node's ID
  //
  int * BestPathCoords = malloc((5*Graph->num_nodes + 1)*sizeof(int));
  int   nodes_in_path  = 0;


  SPLICE_NODE * Node = Graph->Nodes[end_node_id];
  BestPathCoords[5]  = Node->node_id;

  // C-Terminal is special
  if (Node->was_missed) {
    BestPathCoords[3] = Graph->MissedHits->hit[Node->hit_id]->dcl->tr->sqto[0];
    BestPathCoords[4] = Graph->MissedHits->hit[Node->hit_id]->dcl->tr->hmmto[0];
  } else {
    BestPathCoords[3] = (&Graph->TopHits->hit[Node->hit_id]->dcl[0])->tr->sqto[0];
    BestPathCoords[4] = (&Graph->TopHits->hit[Node->hit_id]->dcl[0])->tr->hmmto[0];
  }


  // Iterate until we're at the N-terminal node
  SPLICE_NODE * USNode;
  while (Node->best_in_edge != -1) {


    BestPathCoords[nodes_in_path*5 + 1] = Node->InEdges[Node->best_in_edge]->downstream_spliced_nucl_start;
    BestPathCoords[nodes_in_path*5 + 2] = Node->InEdges[Node->best_in_edge]->downstream_exon_terminus;


    // We need to know which edge fed into the current node
    // from its upstream friend.
    USNode = Node->UpstreamNodes[Node->best_in_edge];

    int out_edge_index = 0;
    while (USNode->DownstreamNodes[out_edge_index] != Node)
      out_edge_index++;


    // Swell!
    nodes_in_path++;
    Node = USNode;
    BestPathCoords[nodes_in_path*5 + 5] = Node->node_id;

    BestPathCoords[nodes_in_path*5 + 3] = Node->OutEdges[out_edge_index]->upstream_spliced_nucl_end;
    BestPathCoords[nodes_in_path*5 + 4] = Node->OutEdges[out_edge_index]->upstream_exon_terminus;


  }


  // Complete the list with the N-terminal node's info
  if (Node->was_missed) {
    BestPathCoords[nodes_in_path*5 + 1] = Graph->MissedHits->hit[Node->hit_id]->dcl->tr->sqfrom[0];
    BestPathCoords[nodes_in_path*5 + 2] = Graph->MissedHits->hit[Node->hit_id]->dcl->tr->hmmfrom[0];
  } else {
    BestPathCoords[nodes_in_path*5 + 1] = (&Graph->TopHits->hit[Node->hit_id]->dcl[0])->tr->sqfrom[0];
    BestPathCoords[nodes_in_path*5 + 2] = (&Graph->TopHits->hit[Node->hit_id]->dcl[0])->tr->hmmfrom[0];
  }
  nodes_in_path++;



  // Because we started with the C-terminal exon, these coordinates
  // are in reverse order -- flip 'em!
  int *ExonCoordSet = malloc((5*nodes_in_path + 1)*sizeof(int));
  
  int read_id = nodes_in_path - 1;
  int exon_id,i;
  for (exon_id = 0; exon_id < nodes_in_path; exon_id++) {
  
    for (i=1; i<=5; i++)
      ExonCoordSet[5*exon_id+i] = BestPathCoords[5*read_id+i];
  
    read_id--;
  
  }
  ExonCoordSet[0] = nodes_in_path;
  free(BestPathCoords);


  if (DEBUGGING1) DEBUG_OUT("'GetExonSetFromEndNode' Complete",-1);


  return ExonCoordSet;

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: FindComponentBestEnd
 *
 *  Desc. :
 *
 *  Inputs:  1.            Node :
 *           2.    ComponentIDs :
 *           3.    component_id :
 *           4.   best_comp_end :
 *           5. best_comp_score :
 *
 *  Output:
 *
 */
void FindComponentBestEnd
(
  SPLICE_NODE * Node,
  int   * ComponentIDs,
  int     component_id,
  int   * comp_best_end,
  float * comp_best_score
)
{
  if (DEBUGGING1) DEBUG_OUT("Starting 'FindComponentBestEnd'",1);


  if (ComponentIDs[Node->node_id]) {
    if (DEBUGGING1) DEBUG_OUT("'FindComponentBestEnd' Complete",-1);
    return;
  }

  ComponentIDs[Node->node_id] = component_id;

  if (Node->num_out_edges == 0)
 
  if (Node->num_out_edges == 0 && Node->cumulative_score > *comp_best_score) {
    *comp_best_end   = Node->node_id;
    *comp_best_score = Node->cumulative_score;
  }

  int in_edge_index;
  for (in_edge_index=0; in_edge_index<Node->num_in_edges; in_edge_index++)
    FindComponentBestEnd(Node->UpstreamNodes[in_edge_index],ComponentIDs,component_id,comp_best_end,comp_best_score);

  int out_edge_index;
  for (out_edge_index=0; out_edge_index<Node->num_out_edges; out_edge_index++)
    FindComponentBestEnd(Node->DownstreamNodes[out_edge_index],ComponentIDs,component_id,comp_best_end,comp_best_score);

  if (DEBUGGING1) DEBUG_OUT("'FindComponentBestEnd' Complete",-1);

}











/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: TranslateExonSetNucls
 *
 *  Desc. :
 *
 *  Inputs:  1.      ExonSetNucls :
 *           2. coding_region_len :
 *           3.             gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:
 *
 */
ESL_DSQ * TranslateExonSetNucls
(ESL_DSQ * ExonSetNucls, int coding_region_len, ESL_GENCODE * gcode)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'TranslateExonSetNucls'",1);

  int translation_len = coding_region_len/3;
  ESL_DSQ * ExonSetTrans = malloc((2 + translation_len) * sizeof(ESL_DSQ));

  int trans_index;
  for (trans_index=1; trans_index<=translation_len; trans_index++) 
      ExonSetTrans[trans_index] = esl_gencode_GetTranslation(gcode,&(ExonSetNucls[3*trans_index-2]));


  if (DEBUGGING1) DEBUG_OUT("'TranslateExonSetNucls' Complete",-1);


  return ExonSetTrans;

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: grab_exon_coord_set_nucls
 *
 *  Desc. :
 *
 *  Inputs:  1.      ExonCoordSet :
 *           2.     TargetNuclSeq :
 *           3. coding_region_len :
 *
 *  Output:
 *
 */
ESL_DSQ * grab_exon_coord_set_nucls
(int * ExonCoordSet, ESL_SQ *target_seq, int complementarity, int * coding_region_len)
{

  int num_exons = ExonCoordSet[0];

  int num_nucls = 0;
  int exon_id,nucl_id;
  for (exon_id = 0; exon_id < num_exons; exon_id++) 
    num_nucls += abs(ExonCoordSet[exon_id*5 + 1] - ExonCoordSet[exon_id*5 + 3]) + 1;
  
  ESL_DSQ * ExonSetNucls = malloc((num_nucls+2) * sizeof(ESL_DSQ));
  *coding_region_len = num_nucls;

  int nucl_placer = 1;
  for (exon_id = 0; exon_id < num_exons; exon_id++) {

    int range_start = ExonCoordSet[exon_id*5 + 1];
    int range_end   = ExonCoordSet[exon_id*5 + 3];
    ESL_DSQ *ExonNucls = grab_seq_range(target_seq, range_start, range_end, complementarity);  

    for (nucl_id = 1; nucl_id <= abs(range_end-range_start)+1; nucl_id++)
      ExonSetNucls[nucl_placer++] = ExonNucls[nucl_id];

    free(ExonNucls);
  }

  return ExonSetNucls;
}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  DEBUGGING Function: DumpExonSetSequence
 *
 *  Desc. :
 *
 *  Inputs:  1. ExonCoordSets :
 *           2. num_exon_sets :
 *           3. TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           4.         gcode : An ESL_GENCODE struct (mainly used for translation).
 *
 *  Output:
 *
 */
void DumpExonSets
(int ** ExonCoordSets, int num_exon_sets, ESL_SQ *target_seq, ESL_GENCODE * gcode, int complementarity)
{
  int i,exon_set_id;

  fprintf(stderr,"\n\n+");
  for (i=0; i<60; i++)
    fprintf(stderr,"=");
  fprintf(stderr,"+\n\n");


  
  for (exon_set_id=0; exon_set_id<num_exon_sets; exon_set_id++) {


    int * ExonCoords = ExonCoordSets[exon_set_id];
    int   num_exons  = ExonCoords[0];


    int num_nucls;
    ESL_DSQ * NuclSeq  = grab_exon_coord_set_nucls (ExonCoordSets[exon_set_id], target_seq, complementarity, &num_nucls);
    ESL_DSQ * TransSeq = TranslateExonSetNucls(NuclSeq,num_nucls,gcode);
    


    fprintf(stderr,">ExonSet__%d/%d:Nucls__",exon_set_id+1,num_exon_sets);
    for (i=0; i<num_exons; i++) {
      if (i) fprintf(stderr,",");
      fprintf(stderr,"%d-%d",ExonCoords[i*5 + 1],ExonCoords[i*5 + 3]);
    }


    fprintf(stderr,":Aminos__");
    for (i=0; i<num_exons; i++) {
      if (i) fprintf(stderr,",");
      fprintf(stderr,"%d-%d",ExonCoords[i*5 + 2],ExonCoords[i*5 + 4]);
    }


    int line_length = 60;
    for (i=1; i<=num_nucls; i++) {
      if (i % line_length == 1)
        fprintf(stderr,"\n");
      fprintf(stderr,"%c",DNA_CHARS[NuclSeq[i]]);
    }
    fprintf(stderr,"\n");


    fprintf(stderr,">Translation");
    for (i=1; i<=num_nucls/3; i++) {
      if (i % line_length == 1)
        fprintf(stderr,"\n");
      fprintf(stderr,"%c",AMINO_CHARS[TransSeq[i]]);
    }
    fprintf(stderr,"\n");


    free(NuclSeq);
    free(TransSeq);

    fprintf(stderr,"\n");

  }
  fprintf(stderr,"\n");


  fprintf(stderr,"+");
  for (i=0; i<60; i++)
    fprintf(stderr,"=");
  fprintf(stderr,"+\n\n\n");

}













/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: ExonSetCleanup
 *
 */
int ** ExonSetCleanup
(SPLICE_GRAPH * Graph, int * InputCoordSet, int * num_split_sets)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'ExonSetCleanup'",1);


  int num_exons = InputCoordSet[0];


  // Because we can't immediately guarantee reading frame preservation,
  // when we find an ugly exon (coordinates that don't look like what
  // we expect) we just pluck it out AND recover the start/end coordinates
  // for the adjacent exons from the original hits
  //
  // Begin by marking extranasty bois with '-1' in the +5 (node id).
  // We do a full first pass with this marking so that we don't
  // end up checking against a bad downstream node in the second pass.
  //
  int exon_id;
  for (exon_id = 0; exon_id < num_exons; exon_id++) {

    int exon_hmm_start = InputCoordSet[5*exon_id+2];
    int exon_hmm_end   = InputCoordSet[5*exon_id+4];

    // Non-positive model movement
    if (exon_hmm_end - exon_hmm_start <= 0) {
      InputCoordSet[5*exon_id+5] = -1;
      continue;
    }

    // While we're at it -- pure insanity?
    if (exon_hmm_start <= 0 || exon_hmm_end <= 0) {
      InputCoordSet[5*exon_id+5] = -1;
      continue;
    }

  }


  // Now we can check to make sure each node's coordinates make sense
  // with its up- and downstream friends (individually)
  //
  // NOTE: These catches shouldn't ever be necessary, but since they've
  //       occurred to me as conceivable I'm just going to implement them
  //       out of an abundance of paranoia.
  //
  for (exon_id = 1; exon_id < num_exons-1; exon_id++) {

    if (InputCoordSet[5*exon_id+5] == -1)
      continue;


    int exon_hmm_start = InputCoordSet[5*exon_id+2];
    int exon_hmm_end   = InputCoordSet[5*exon_id+4];


    if (InputCoordSet[5*(exon_id-1)+5 != -1]) {

      int upstream_hmm_end = InputCoordSet[5*(exon_id-1)+4];

      if (upstream_hmm_end > exon_hmm_start) {
        InputCoordSet[5*exon_id+5] = -1;
        continue;
      }

    }


    if (InputCoordSet[5*(exon_id+1)+5] != -1) {

      int downstream_hmm_start = InputCoordSet[5*(exon_id+1)+2];

      if (downstream_hmm_start < exon_hmm_end) {
        InputCoordSet[5*exon_id+5] = -1;
        continue;
      }

    }

  }


  // Last catch -- do the upstream and downstream "swallow" this
  // node?  NOTE that we should have already precluded this, but
  // I refuse to doubt that the computer is conspiring against me.
  //
  for (exon_id = 1; exon_id < num_exons-1; exon_id++) {

    if (InputCoordSet[ 5*(exon_id-1)+5] == -1 
      || InputCoordSet[5* exon_id   +5] == -1 
      || InputCoordSet[5*(exon_id+1)+5] == -1)
      continue;


    int   upstream_hmm_end   = InputCoordSet[5*(exon_id-1)+4];
    int downstream_hmm_start = InputCoordSet[5*(exon_id+1)+2];

    if (upstream_hmm_end+1 >= downstream_hmm_start)
      InputCoordSet[5*exon_id+5] = -1;

  }


  // Now that we've caught all of our problem exons, we'll excise
  // them from the coordinate set, splitting the coordinate set into
  // non-problematic pieces.
  //
  // We'll start by counting how many pieces there will be
  //
  int num_sub_sets = 1;
  for (exon_id = 0; exon_id < num_exons; exon_id++) {
    if (InputCoordSet[5*exon_id+5] == -1)
      num_sub_sets++;
  }


  int ** SplitCoordSets = malloc(num_sub_sets * sizeof(int *));


  // This is what I expect to happen 99.9% of the time...
  //
  // Note that we don't just return '&InputCoordSet' because
  // that would require additional bookkeeping to avoid double-frees.
  // A little wasteful, but *whatever*
  //
  if (num_sub_sets == 1) {
    
    SplitCoordSets[0] = malloc((5*num_exons+1)*sizeof(int));

    int i;
    for (i=0; i<=5*num_exons; i++)
      SplitCoordSets[0][i] = InputCoordSet[i];


    if (DEBUGGING1) DEBUG_OUT("'ExonSetCleanup' Complete",-1);


    *num_split_sets = 1;
    return SplitCoordSets;


  }


  // BUMMER!
  //
  // NOTE that 'num_sub_hits' is the maximum number of subhits.
  // If we only remove terminal exons, then we might still have only
  // one subset of exons to return.
  //
  *num_split_sets = 0;
  int scanner_exon_id = 0;
  while (scanner_exon_id < num_exons) {

    while (scanner_exon_id < num_exons && InputCoordSet[5*scanner_exon_id+5] == -1)
      scanner_exon_id++;
    
    if (scanner_exon_id == num_exons)
      break;

    int start_exon_id = scanner_exon_id;

    while (scanner_exon_id < num_exons && InputCoordSet[5*scanner_exon_id+5] != -1) 
      scanner_exon_id++;

    int end_exon_id = scanner_exon_id-1;


    // If we aren't starting with the first exon, we need to recover the
    // start bounds of the actual hit (as opposed to the now-jettisoned splice
    // coordinates)
    if (start_exon_id > 0) {


      SPLICE_NODE * Node = Graph->Nodes[InputCoordSet[5*start_exon_id+5]];
      P7_TRACE * TR;
      if (Node->was_missed)
        TR = Graph->MissedHits->hit[Node->hit_id]->dcl->tr;
      else
        TR = (&Graph->TopHits->hit[Node->hit_id]->dcl[0])->tr;


      int exon_hmm_start  = TR->hmmfrom[0];
      int exon_nucl_start = TR->sqfrom[0];

      InputCoordSet[5*start_exon_id+1] = exon_nucl_start;
      InputCoordSet[5*start_exon_id+2] = exon_hmm_start;


    }


    // Similarly, if we aren't ending with the last exon we need to go back
    // to what the hit tells us.
    if (end_exon_id < num_exons-1) {


      SPLICE_NODE * Node = Graph->Nodes[InputCoordSet[5*end_exon_id+5]];
      P7_TRACE *TR;
      if (Node->was_missed) {
        TR = Graph->MissedHits->hit[Node->hit_id]->dcl->tr;
      }
      else {
        TR = (&Graph->TopHits->hit[Node->hit_id]->dcl[0])->tr;
      }

     
   
      int exon_hmm_end  = TR->hmmto[0];
      int exon_nucl_end = TR->sqto[0];

      InputCoordSet[5*end_exon_id+3] = exon_nucl_end;
      InputCoordSet[5*end_exon_id+4] = exon_hmm_end;


    }


    // Record this exon set!
    
    int num_sub_exons = end_exon_id - start_exon_id + 1;
    SplitCoordSets[*num_split_sets] = malloc((5*num_sub_exons+1)*sizeof(int));

    SplitCoordSets[*num_split_sets][0] = num_sub_exons;

    int sub_exon_id = 0;
    for (exon_id = start_exon_id; exon_id <= end_exon_id; exon_id++) {
      SplitCoordSets[*num_split_sets][5*sub_exon_id+1] = InputCoordSet[5*exon_id+1]; 
      SplitCoordSets[*num_split_sets][5*sub_exon_id+2] = InputCoordSet[5*exon_id+2]; 
      SplitCoordSets[*num_split_sets][5*sub_exon_id+3] = InputCoordSet[5*exon_id+3]; 
      SplitCoordSets[*num_split_sets][5*sub_exon_id+4] = InputCoordSet[5*exon_id+4]; 
      SplitCoordSets[*num_split_sets][5*sub_exon_id+5] = InputCoordSet[5*exon_id+5]; 
      sub_exon_id++;
    }

    *num_split_sets += 1;


  } 




  return SplitCoordSets;


}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetSplicedExonCoordSets
 *
 *  Desc. :
 *
 *  Inputs:  1.         Graph :
 *           2. num_exon_sets :
 *
 *  Output:
 *
 */
int ** GetSplicedExonCoordSets
(SPLICE_GRAPH * Graph, int * num_coord_sets)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'GetSplicedExonCoordSets'",1);


  // How many connected components do we have?
  // (In the case that we have a full path through
  //  the graph, we just report that one component.)
  int num_conn_comps = 0;
  int * EndNodes = malloc(Graph->num_nodes*sizeof(int));
  int node_id;

  if (Graph->has_full_path) {
    EndNodes[0]  = Graph->best_full_path_end;
    num_conn_comps = 1;

  } else {

    int * ComponentIDs = malloc((Graph->num_nodes+1)*sizeof(int));
    for (node_id=1; node_id<=Graph->num_nodes; node_id++)
      ComponentIDs[node_id] = 0;
    
    for (node_id=1; node_id<=Graph->num_nodes; node_id++) {

      if (!ComponentIDs[node_id]) {

        int   comp_best_end   = 0;
        float comp_best_score = 0.0;
        FindComponentBestEnd(Graph->Nodes[node_id],ComponentIDs,num_conn_comps+1,&comp_best_end,&comp_best_score);

        // There's some way that we can get a node to not be properly tagged
        // so that it *thinks* its connected component hasn't been considered.
        // This results in 'FindComponentBestStart' returning without setting
        // 'best_comp_start'
        // This *is* a bug, but for now I'm just patching over it
        if (comp_best_end != 0)
          EndNodes[num_conn_comps++] = comp_best_end;

      }

    }

    free(ComponentIDs);

  }



  // Cool!  Now let's grab the spliced coordinates
  //
  int ecs_cap = num_conn_comps * 2;
  int ** ExonCoordSets = malloc(ecs_cap * sizeof(int *));

  *num_coord_sets = 0;
  int conn_comp_id, sub_set_id;
  for (conn_comp_id = 0; conn_comp_id < num_conn_comps; conn_comp_id++) {


    int  * InitCoordSet   = GetExonSetFromEndNode(Graph,EndNodes[conn_comp_id]);
    int    num_split_sets = 0;
    int ** SplitCoordSets = ExonSetCleanup(Graph,InitCoordSet,&num_split_sets);


    for (sub_set_id = 0; sub_set_id < num_split_sets; sub_set_id++) {


      // Copy the next subset into ExonCoordSets

      int * SubCoordSet = SplitCoordSets[sub_set_id];
      ExonCoordSets[*num_coord_sets] = malloc((5*SubCoordSet[0]+1) * sizeof(int));
      int x;
      for (x=0; x<=5*SubCoordSet[0]; x++)
        ExonCoordSets[*num_coord_sets][x] = SubCoordSet[x];
      *num_coord_sets += 1;

      free(SubCoordSet);


      // Resize?
      if (*num_coord_sets >= ecs_cap) {
        ecs_cap *= 2;
        int ** TmpECS = malloc(ecs_cap * sizeof(int *));
        for (x=0; x<*num_coord_sets; x++)
          TmpECS[x] = ExonCoordSets[x];
        free(ExonCoordSets); 
        ExonCoordSets = TmpECS;
      }


    }


    // Free up and move forward!
    free(SplitCoordSets);
    free(InitCoordSet);

  }

  free(EndNodes);


  if (DEBUGGING1) DEBUG_OUT("'GetSplicedExonCoordSets' Complete",-1);


  // Easy!
  return ExonCoordSets;

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: RightAlignStr
 *
 */
char * RightAlignStr
(char * InStr, int out_str_len)
{
  int in_str_len = strlen(InStr);
  
  if (in_str_len > out_str_len) {
    fprintf(stderr,"\n  RAS ERROR:  Input string length (%d,'%s') longer than requested output string length (%d)\n\n",in_str_len,InStr,out_str_len);
    return InStr;
  }  

  char * OutStr = malloc((out_str_len+1) * sizeof(char));
  OutStr[out_str_len] = 0;

  int writer = out_str_len-1;
  int i;
  for (i=in_str_len-1; i>=0; i--)
    OutStr[writer--] = InStr[i];

  while (writer >= 0)
    OutStr[writer--] = ' ';

  return OutStr;

}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: IntToCharArr
 *
 */
char * IntToCharArr
(int in_val, int * char_arr_len)
{

  *char_arr_len = 1;

  int is_negative = 0;
  if (in_val < 0) {
    in_val *= -1;
    *char_arr_len = 2;
    is_negative = 1;
  }

  int in_val_copy = in_val;
  while (in_val_copy > 9) {
    in_val_copy /= 10;
    *char_arr_len += 1;
  }

  char * CharArr = malloc((*char_arr_len + 1) * sizeof(char));
  CharArr[*char_arr_len] = 0;
  
  int writer  = *char_arr_len - 1;
  in_val_copy = in_val;
  while (writer >= is_negative) {
    int mod_10 = in_val_copy % 10;
    in_val_copy /= 10;
    CharArr[writer--] = (char)(48 + mod_10);
  }
  if (is_negative) CharArr[0] = '-';

  return CharArr;

}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetALSBlockLengths // 'Alignment Left Side' (deprecated initialism)
 *
 */
void GetALSBlockLengths
(
  P7_ALIDISPLAY * AD, 
  int * block_A_len,  
  int * ExonCoordSet, 
  int * block_B_len
)
{
  
  *block_A_len = intMax(8,intMax(intMax(strlen(AD->hmmname),strlen(AD->sqname)),strlen(AD->orfname)));

  int max_coord = 0;
  int exon_id;
  for (exon_id=0; exon_id<ExonCoordSet[0]; exon_id++) {
    int set_max = intMax(intMax(ExonCoordSet[(exon_id*5)+1],ExonCoordSet[(exon_id*5)+3]),intMax(ExonCoordSet[(exon_id*5)+2],ExonCoordSet[(exon_id*5)+4]));
    max_coord = intMax(max_coord,set_max);
  }

  char * max_coord_str = IntToCharArr(max_coord,block_B_len);
  free(max_coord_str);

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: print_exon
 *
 */
void print_exon
(
  EXON_DISPLAY_INFO * EDI, 
  ESL_SQ     *target_seq,
  int * ad_nucl_read_pos,
  int * ad_amino_read_pos,
  int * codon_pos
)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'print_exon'",1);


  FILE          * ofp = EDI->ofp;
  P7_ALIDISPLAY * AD  = EDI->AD;
  
  ESL_DSQ * LeftDinucls  = NULL;
  ESL_DSQ * RightDinucls = NULL;
  int ali_nucl_start = EDI->nucl_start;
  if (EDI->revcomp) {

    if (!EDI->Nterm) {
      LeftDinucls = grab_seq_range(target_seq, EDI->nucl_start+2, EDI->nucl_start+1, EDI->revcomp);
      ali_nucl_start += 2;
    }

    if (!EDI->Cterm) {
      RightDinucls = grab_seq_range(target_seq, EDI->nucl_end-1, EDI->nucl_end-2, EDI->revcomp);
    }

  } else {

    if (!EDI->Nterm) {
      LeftDinucls = grab_seq_range(target_seq, EDI->nucl_start-2, EDI->nucl_start-1, EDI->revcomp);
      ali_nucl_start -= 2;
    }

    if (!EDI->Cterm) {
      RightDinucls = grab_seq_range(target_seq, EDI->nucl_end+1, EDI->nucl_end+2, EDI->revcomp);
    }

  }


  // This is a rare event (all observations during testing summed up to
  // 20 cases, all in rat (example: ARX)), but there might be Ns infiltrating
  // our dinucleotide sequences, so we need to catch those
  if (LeftDinucls) {
    if ( LeftDinucls[1] > 3)  LeftDinucls[1] = 5;
    if ( LeftDinucls[2] > 3)  LeftDinucls[2] = 5;
  }
  if (RightDinucls) {
    if (RightDinucls[1] > 3) RightDinucls[1] = 5;
    if (RightDinucls[2] > 3) RightDinucls[2] = 5;
  }



  // Prep for each row in the alignment.  Extra length in case of gaps.
  int exon_ali_str_alloc = 2 * (abs(EDI->nucl_end - EDI->nucl_start) + 1);
  char * ModelAli   = malloc(exon_ali_str_alloc * sizeof(char));
  char * QualityAli = malloc(exon_ali_str_alloc * sizeof(char));
  char * TransAli   = malloc(exon_ali_str_alloc * sizeof(char));
  char * NuclAli    = malloc(exon_ali_str_alloc * sizeof(char));
  char * PPAli      = malloc(exon_ali_str_alloc * sizeof(char));
  
  // Left splicing dinucleotides
  int write_pos = 0;
  if (!EDI->Nterm) {
    ModelAli[0]   = ' ';
    ModelAli[1]   = ' ';
    QualityAli[0] = ' ';
    QualityAli[1] = ' ';
    TransAli[0]   = ' ';
    TransAli[1]   = ' ';
    NuclAli[0]    = LC_DNA_CHARS[LeftDinucls[1]];
    NuclAli[1]    = LC_DNA_CHARS[LeftDinucls[2]];
    PPAli[0]      = '|';
    PPAli[1]      = '|';
    write_pos = 2;
  }



  int current_nucl_pos = EDI->nucl_start;
  while ((EDI->revcomp && current_nucl_pos >= EDI->nucl_end)
         || (!EDI->revcomp && current_nucl_pos <= EDI->nucl_end)) {

    NuclAli[write_pos] = AD->ntseq[*ad_nucl_read_pos];
    if (NuclAli[write_pos] != '-') {
      if (EDI->revcomp) {
        current_nucl_pos--;
      } else {
        current_nucl_pos++;
      }
    }
    *ad_nucl_read_pos += 1;


    if (*codon_pos == 1) {
      ModelAli[write_pos]   = AD->model[*ad_amino_read_pos];
      QualityAli[write_pos] = AD->mline[*ad_amino_read_pos];
      TransAli[write_pos]   = AD->aseq[*ad_amino_read_pos];
      PPAli[write_pos]      = AD->ppline[*ad_amino_read_pos];
      *ad_amino_read_pos   += 1;
    } else {
      ModelAli[write_pos]   = ' ';
      QualityAli[write_pos] = ' ';
      TransAli[write_pos]   = ' ';
      PPAli[write_pos]      = ' ';
    }
    write_pos++;


    *codon_pos += 1;
    if (*codon_pos == 3)
      *codon_pos = 0;



    // Similar to the issue of parsing "original" ALIDISPLAYs
    // in 'FindOptimalSpliceSite,' if we're re-allocating
    // it probably means we've ended up with some uncomfortably
    // gappy alignments...
    if (write_pos + 5 >= exon_ali_str_alloc) {
      exon_ali_str_alloc *= 2;
      char * TmpModelAli    = malloc(exon_ali_str_alloc*sizeof(char));
      char * TmpQualityAli  = malloc(exon_ali_str_alloc*sizeof(char));
      char * TmpTransAli    = malloc(exon_ali_str_alloc*sizeof(char));
      char * TmpNuclAli     = malloc(exon_ali_str_alloc*sizeof(char));
      char * TmpPPAli       = malloc(exon_ali_str_alloc*sizeof(char));
      int copy_index;
      for (copy_index=0; copy_index<write_pos; copy_index++) {
        TmpModelAli[copy_index]   = ModelAli[copy_index];
        TmpQualityAli[copy_index] = QualityAli[copy_index];
        TmpTransAli[copy_index]   = TransAli[copy_index];
        TmpNuclAli[copy_index]    = NuclAli[copy_index];
        TmpPPAli[copy_index]      = PPAli[copy_index];
      }
      free(ModelAli);   ModelAli   = TmpModelAli;
      free(QualityAli); QualityAli = TmpQualityAli;
      free(TransAli);   TransAli   = TmpTransAli;
      free(NuclAli);    NuclAli    = TmpNuclAli;
      free(PPAli);      PPAli      = TmpPPAli;
    }

  }


  // Right splicing dinucleotides
  if (!EDI->Cterm) {
    ModelAli[write_pos]     = ' ';
    ModelAli[write_pos+1]   = ' ';
    QualityAli[write_pos]   = ' ';
    QualityAli[write_pos+1] = ' ';
    TransAli[write_pos]     = ' ';
    TransAli[write_pos+1]   = ' ';
    NuclAli[write_pos]      = LC_DNA_CHARS[RightDinucls[1]];
    NuclAli[write_pos+1]    = LC_DNA_CHARS[RightDinucls[2]];
    PPAli[write_pos]        = '|';
    PPAli[write_pos+1]      = '|';
    write_pos += 2;
  }

  // Make the strings *very* printable! (for debugging, mainly...)
  ModelAli[write_pos]   = 0;
  QualityAli[write_pos] = 0;
  TransAli[write_pos]   = 0;
  NuclAli[write_pos]    = 0;
  PPAli[write_pos]      = 0;


  int ali_len = write_pos; // Just for clarity of reading


  // Just wipe these now, since they're stored in the alignment
  if ( LeftDinucls) free( LeftDinucls);
  if (RightDinucls) free(RightDinucls);



  // We're officially well-positioned to print out this exon!
  // Before we start printing the alingment, let's give just
  // a kiss of metadata
  if (ALEX_MODE)
    fprintf(ofp,"\n  %s %s [ Exon Set %d / Exon %d ]\n\n",EDI->NameBlank,EDI->CoordBlank,EDI->exon_set_id,EDI->exon_id);


  // I'm making the executive decision that this is what the
  // Translated Nucleotide String is titled, and nobody can stop
  // me!
  int exon_id_len;
  char * exon_id_str = IntToCharArr(EDI->exon_id,&exon_id_len);
  char * TNSName = malloc((exon_id_len + 6) * sizeof(char));
  TNSName[0] = 'e';
  TNSName[1] = 'x';
  TNSName[2] = 'o';
  TNSName[3] = 'n';
  TNSName[4] = ' ';
  int i;
  for (i=0; i<exon_id_len; i++)
    TNSName[5+i] = exon_id_str[i];
  TNSName[5+exon_id_len] = 0;
  char * FormattedTNSName = RightAlignStr(TNSName,EDI->name_str_len);


  int formatted_int_len; // We just need a pointer for 'IntToCharArr'
  char * CharredInt;
  char * FormattedInt;

  // If they're pranking us, prank 'em right back!
  if (EDI->textw < 1) EDI->textw = ali_len;

  int adj_textw = EDI->textw - (EDI->name_str_len + 2*EDI->coord_str_len + 4);

  int model_pos = EDI->hmm_start;
  int nucl_pos  = ali_nucl_start;
  int line_read_start,read_pos;
  for (line_read_start = 0; line_read_start < ali_len; line_read_start += adj_textw) {

    
    int line_read_end = line_read_start + adj_textw;
    if (line_read_end > ali_len)
      line_read_end = ali_len;


    //////////////////////////////////
    //
    // 1. Model
    //
    CharredInt   = IntToCharArr(model_pos,&formatted_int_len);
    FormattedInt = RightAlignStr(CharredInt,EDI->coord_str_len);

    fprintf(ofp,"  %s %s ",EDI->HMMName,FormattedInt);
    free(CharredInt);
    free(FormattedInt);

    int init_model_pos = model_pos;
    for (read_pos = line_read_start; read_pos < line_read_end; read_pos++) {
      fprintf(ofp,"%c",ModelAli[read_pos]);
      if (ModelAli[read_pos] != ' ' && ModelAli[read_pos] != '.') {
        model_pos++;
      }
    }
    if (model_pos == init_model_pos) {
      fprintf(ofp," %d\n",model_pos);
    } else {
      fprintf(ofp," %d\n",model_pos-1);
    }



    //////////////////////////////////
    //
    // 2. Quality stuff
    //
    fprintf(ofp,"  %s %s ",EDI->NameBlank,EDI->CoordBlank);
    for (read_pos = line_read_start; read_pos < line_read_end; read_pos++)
      fprintf(ofp,"%c",QualityAli[read_pos]);
    fprintf(ofp,"\n");



    //////////////////////////////////
    //
    // 3. Translation
    //
    // OLD: fprintf(ofp,"  %s %s ",EDI->TransName,EDI->CoordBlank);
    fprintf(ofp,"  %s %s ",FormattedTNSName,EDI->CoordBlank);
    for (read_pos = line_read_start; read_pos < line_read_end; read_pos++)
      fprintf(ofp,"%c",TransAli[read_pos]);
    fprintf(ofp,"\n");



    //////////////////////////////////
    //
    // 4. Nucleotides
    //
    CharredInt   = IntToCharArr(nucl_pos,&formatted_int_len);
    FormattedInt = RightAlignStr(CharredInt,EDI->coord_str_len);

    fprintf(ofp,"  %s %s ",EDI->NuclName,FormattedInt);
    free(CharredInt);
    free(FormattedInt);

    int init_nucl_pos = nucl_pos;
    for (read_pos = line_read_start; read_pos < line_read_end; read_pos++) {
      fprintf(ofp,"%c",NuclAli[read_pos]);
      if (NuclAli[read_pos] != '-') {
        if (EDI->revcomp) {
          nucl_pos--;
        } else {
          nucl_pos++;
        }
      }
    }
    if (nucl_pos == init_nucl_pos) {
      fprintf(ofp," %d\n",nucl_pos);
    } else if (EDI->revcomp) {
      fprintf(ofp," %d\n",nucl_pos+1);
    } else {
      fprintf(ofp," %d\n",nucl_pos-1);
    }



    //////////////////////////////////
    //
    // 5. PP Time!
    //
    fprintf(ofp,"  %s %s ",EDI->NameBlank,EDI->CoordBlank);
    for (read_pos = line_read_start; read_pos < line_read_end; read_pos++)
      fprintf(ofp,"%c",PPAli[read_pos]);
    fprintf(ofp,"\n");



    //////////////////////////////////
    //
    // That's another buncha rows!
    //
    fprintf(ofp,"\n");

  }



  free(ModelAli);
  free(QualityAli);
  free(TransAli);
  free(NuclAli);
  free(PPAli);
  free(exon_id_str);
  free(TNSName);
  free(FormattedTNSName);

  if (DEBUGGING1) DEBUG_OUT("'print_exon' Complete",-1);


}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: print_spliced_alignment
 *
 *  Desc. :
 *
 *  Inputs:  
 *
 *  Output:
 *
 */
void print_spliced_alignment
(
  P7_ALIDISPLAY * AD, 
  ESL_SQ        *target_seq,
  int           * ExonCoordSet, 
  int             exon_set_name_id, 
  FILE          * ofp, 
  int             textw
)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'print_spliced_alignment'",1);

  int i;

  EXON_DISPLAY_INFO * EDI = malloc(sizeof(EXON_DISPLAY_INFO));
  EDI->ofp         = ofp;
  EDI->textw       = textw;
  EDI->AD          = AD;
  EDI->exon_set_id = exon_set_name_id;


  GetALSBlockLengths(EDI->AD,&(EDI->name_str_len),ExonCoordSet,&(EDI->coord_str_len));

  EDI->HMMName   = RightAlignStr(AD->hmmname,EDI->name_str_len);
  EDI->TransName = RightAlignStr(AD->orfname,EDI->name_str_len);
  EDI->NuclName  = RightAlignStr(AD->sqname ,EDI->name_str_len);
  EDI->NameBlank = RightAlignStr(" "        ,EDI->name_str_len);


  EDI->CoordBlank = malloc((EDI->coord_str_len+1)*sizeof(char));
  for (i=0; i<EDI->coord_str_len; i++)
    EDI->CoordBlank[i] = ' ';
  EDI->CoordBlank[EDI->coord_str_len] = 0;


  EDI->revcomp = 0;
  if (ExonCoordSet[1] > ExonCoordSet[3])
    EDI->revcomp = 1;
  
  int codon_pos         = 0;
  int ad_nucl_read_pos  = 0;
  int ad_amino_read_pos = 0;

  int num_exons = ExonCoordSet[0];
  int exon_id;
  for (exon_id = 0; exon_id < num_exons; exon_id++) {

    EDI->exon_id     = exon_id + 1;
    EDI->hmm_start   = ExonCoordSet[(5*exon_id)+2];
    EDI->hmm_end     = ExonCoordSet[(5*exon_id)+4];
    EDI->nucl_start  = ExonCoordSet[(5*exon_id)+1];
    EDI->nucl_end    = ExonCoordSet[(5*exon_id)+3];

    if (exon_id == 0) EDI->Nterm = 1;
    else              EDI->Nterm = 0;

    if (exon_id == num_exons-1) EDI->Cterm = 1;
    else                        EDI->Cterm = 0;

    print_exon (EDI, target_seq,&ad_nucl_read_pos,&ad_amino_read_pos,&codon_pos);
  }


  free(EDI->HMMName);
  free(EDI->TransName);
  free(EDI->NuclName);
  free(EDI->NameBlank);
  free(EDI->CoordBlank);
  free(EDI);


  if (DEBUGGING1) DEBUG_OUT("'print_spliced_alignment' Complete",-1);

}







/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: dump_splash_header 
 *
 */
void dump_splash_header
(
  SPLICE_GRAPH * Graph,
  ESL_SQ     *target_seq,
  int exon_set_name_id,
  int * ExonCoordSet, 
  FILE * ofp, 
  int textw
)
{

  int exon_id, i;

  int num_exons   = ExonCoordSet[0];
  int model_start = ExonCoordSet[2];
  int model_end   = ExonCoordSet[4 + 5*(num_exons-1)];
  int nucl_start  = ExonCoordSet[1];
  int nucl_end    = ExonCoordSet[3 + 5*(num_exons-1)];


  int full_coverage = 0;
  if (model_start == 1 && model_end == Graph->Model->M) 
    full_coverage = 1;


  int num_found_exons = 0;
  for (exon_id=0; exon_id<num_exons; exon_id++) {
    int node_id = ExonCoordSet[5 + 5*exon_id];
    if (Graph->Nodes[node_id]->was_missed)
      num_found_exons++;
  }


  fprintf(ofp,"\n\n+");
  for (i=0; i<textw-2; i++)
    fprintf(ofp,"=");
  fprintf(ofp,"+\n");
  fprintf(ofp,"|\n");
  fprintf(ofp,"| splash - spliced alignment of some hits\n");
  fprintf(ofp,"|\n");
  fprintf(ofp,"| = Exon Set %d (%d exons)\n",exon_set_name_id,num_exons);
  fprintf(ofp,"| = Model Positions %d..%d",model_start,model_end);
  if (full_coverage) 
    fprintf(ofp,"  (* Full Model)");
  fprintf(ofp,"\n");
  fprintf(ofp,"| = Target Seq Name %s\n",target_seq->name);
  fprintf(ofp,"| = Nucleotide Coords %d..%d\n",nucl_start,nucl_end);
  for (exon_id=0; exon_id<num_exons; exon_id++)
    fprintf(ofp,"| = Exon %d: %d..%d / %d..%d\n",exon_id+1,ExonCoordSet[5*exon_id+2],ExonCoordSet[5*exon_id+4],ExonCoordSet[5*exon_id+1],ExonCoordSet[5*exon_id+3]);
  if (num_found_exons) 
    fprintf(ofp,"| + Includes Missed Exons\n");
  fprintf(ofp,"|\n");
  fprintf(ofp,":\n");


}








void DumpSplashFooter
(FILE * ofp, int textw)
{
  fprintf(ofp,":\n");
  fprintf(ofp,"|\n");
  fprintf(ofp,"+");
  int i;
  for (i=0; i<textw-2; i++)
    fprintf(ofp,"=");
  fprintf(ofp,"+\n");
  fprintf(ofp,"\n\n"); 
}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: DetermineHitExonCoords
 *
 *  Desc. :
 *
 */
int * DetermineHitExonCoords
(P7_ALIDISPLAY * AD, int * ExonCoords)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'DetermineHitExonCoords'",1);


  // As a sanity check, we want to make sure that the pipeline
  // hit to the same place as where we were hoping (I think
  // this is mainly an issue in highly repetitive sequence)
  //
  int orig_num_exons   = ExonCoords[0];
  int orig_model_start = ExonCoords[2];
  int orig_model_end   = ExonCoords[5*(orig_num_exons-1)+4];
  int orig_nucl_start  = ExonCoords[1];
  int orig_nucl_end    = ExonCoords[5*(orig_num_exons-1)+3];
  if (AD->hmmto <= orig_model_start || AD->hmmfrom >= orig_model_end)
    return NULL;

  int revcomp = 0;
  if (orig_nucl_start > orig_nucl_end)
    revcomp = 1;


  /*
  if (DEBUGGING) {
    fprintf(stdout,"\n\n\n");
    fprintf(stdout,"Original Input: %d Exons\n",ExonCoords[0]);
    int x;
    for (x=0; x<ExonCoords[0]; x++)
      fprintf(stdout,"----> %d..%d / %d..%d\n",ExonCoords[5*x+1],ExonCoords[5*x+3],ExonCoords[5*x+2],ExonCoords[5*x+4]);
    fprintf(stdout,"\n");
    fprintf(stdout,"AD Range: %d..%d\n",AD->hmmfrom,AD->hmmto);
    fprintf(stdout,"\n");
  }
  */


  // Advance through the exon coordinates until we've hit the
  // model position that kicks off the alignment that the
  // HMMER pipeline wants us to use.
  //
  int model_pos = orig_model_start;
  int nucl_pos  = orig_nucl_start;


  // I think this is a bug in the pipeline, but there are *very rare*
  // situations where I've observed we end up doing some serious gapping
  // that pushes the model beyond what we're equipped to handle with our
  // coordinates.
  //
  int map_model_cap = orig_model_end;


  // We want to track the center nucleotide of the codon for this
  // portion
  if (revcomp) nucl_pos--;
  else         nucl_pos++;

  int exon_id = 0;
  while (model_pos < AD->hmmfrom) {

    if (revcomp) {
        
      nucl_pos -= 3;

      // Jump to the next exon?
      if (nucl_pos <= ExonCoords[5*exon_id+3]) {

        // We need to make sure we account for whatever the splice offset is
        int nucl_offset = ExonCoords[5*exon_id+3] - nucl_pos;

        exon_id++;
        nucl_pos = ExonCoords[5*exon_id+1] - nucl_offset;

      }

    } else {

      nucl_pos += 3;

      if (nucl_pos >= ExonCoords[5*exon_id+3]) {

        int nucl_offset = nucl_pos - ExonCoords[5*exon_id+3]; 

        exon_id++;
        nucl_pos = ExonCoords[5*exon_id+1] + nucl_offset;

      }
        
    }

    model_pos++;

  }


  // Initialize our output coordinate set to have as many exons
  // as the input coordinate set (ideally, all we're going to do
  // is create a copy of 'ExonCoords'...)
  int * HitExonCoords = malloc((1 + 5 * orig_num_exons) * sizeof(int));

  // Because we're tracking the middle nucleotide, we'll need to make
  // a little adjustment
  if (revcomp) HitExonCoords[1] = nucl_pos+1;
  else         HitExonCoords[1] = nucl_pos-1;
  HitExonCoords[2] = model_pos;

  int final_num_exons = 0;
  int ali_pos = 0;
  while (model_pos < AD->hmmto && model_pos < map_model_cap) {


    if (AD->model[ali_pos] != '.')
      model_pos++;
    
    if (AD->aseq[ali_pos] != '-') {
      if (revcomp) nucl_pos -= 3;
      else         nucl_pos += 3;
    }

    ali_pos++;


    // Did we just pass into a new exon?
    if (model_pos > ExonCoords[5*exon_id+4]) {

      HitExonCoords[5*final_num_exons+3] = ExonCoords[5*exon_id+3];
      HitExonCoords[5*final_num_exons+4] = model_pos-1;
      HitExonCoords[5*final_num_exons+5] = ExonCoords[5*exon_id+5];

      if (revcomp) {
        int nucl_offset = (ExonCoords[5*exon_id+3] - nucl_pos) - 1;
        exon_id++;
        nucl_pos = ExonCoords[5*exon_id+1] - nucl_offset;
      } else {
        int nucl_offset = (nucl_pos - ExonCoords[5*exon_id+3]) - 1;
        exon_id++;
        nucl_pos = ExonCoords[5*exon_id+1] + nucl_offset;
      }

      final_num_exons++;
      HitExonCoords[5*final_num_exons+1] = ExonCoords[5*exon_id+1];
      HitExonCoords[5*final_num_exons+2] = model_pos;

    }


  }


  if (revcomp) nucl_pos--;
  else         nucl_pos++;
  HitExonCoords[5*final_num_exons+3] = nucl_pos;
  HitExonCoords[5*final_num_exons+4] = model_pos;
  HitExonCoords[5*final_num_exons+5] = ExonCoords[5*exon_id+5];

  HitExonCoords[0] = final_num_exons+1;


  if (DEBUGGING1) DEBUG_OUT("'DetermineHitExonCoords' Complete",-1);


  return HitExonCoords;

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: report_spliced_tophits
 *
 *  Desc. :
 *
 *  Inputs:  1.  ExonSetTopHits :
 *           2. ExonSetPipeline :
 *           3.    ExonCoordSet :
 *           4.             ofp : An open file pointer (specificially, the target for output).
 *           5.           textw : The desired line width for output.
 *
 *  Output:
 *
 */
void report_spliced_tophits
(
  SPLICE_GRAPH * Graph,
  P7_TOPHITS   * ExonSetTopHits, 
  P7_PIPELINE  * ExonSetPipeline, 
  ESL_SQ       *target_seq,
  int          * ExonCoordSet,
  int          * exon_set_name_id,
  FILE         * ofp,
  int            textw
)
{


  // FOR NOW: I'm just going to use this printing style until
  //          we're getting to this point with *every* test case,
  //          and then I'll get the 'ExonCoordSet' data integrated. 
  //
  p7_tophits_SortBySeqidxAndAlipos(ExonSetTopHits);
  p7_tophits_RemoveDuplicates(ExonSetTopHits,ExonSetPipeline->use_bit_cutoffs);
  p7_tophits_SortBySortkey(ExonSetTopHits);
  p7_tophits_Threshold(ExonSetTopHits,ExonSetPipeline);


  // Give the standard metadata
  //
  //fprintf(ofp,"\n\n");
  //p7_tophits_Targets(ofp, ExonSetTopHits, ExonSetPipeline, textw); 
  //fprintf(ofp,"\n\n");
  

  if (DEBUGGING2) p7_tophits_Domains(ofp, ExonSetTopHits, ExonSetPipeline, textw);


  int hit_id;
  for (hit_id = 0; hit_id < (int)(ExonSetTopHits->N); hit_id++) {
      int dom_id = ExonSetTopHits->hit[hit_id]->best_domain;      
      P7_ALIDISPLAY * AD  = (&ExonSetTopHits->hit[hit_id]->dcl[dom_id])->ad;
      P7_TRACE * TR       = (&ExonSetTopHits->hit[hit_id]->dcl[dom_id])->tr;
      if (excessive_gap_content(TR)) continue;

      int * HitExonCoords = DetermineHitExonCoords(AD,ExonCoordSet);

      // If the pipeline hit to a totally different place from
      // where we were expecting, this is the 'bail' code
      if (HitExonCoords == NULL)
        continue;

      *exon_set_name_id += 1;


      // DEBUGGING
      /*
      fprintf(stdout,"\n");
      fprintf(stdout,"Num Exons: %d\n",HitExonCoords[0]);
      int x;
      for (x=0; x<HitExonCoords[0]; x++)
        fprintf(stdout,"  -> %d..%d / %d..%d\n",HitExonCoords[5*x+1],HitExonCoords[5*x+3],HitExonCoords[5*x+2],HitExonCoords[5*x+4]);
      fprintf(stdout,"\n");
      */
      // DEBUGGING


      if (HitExonCoords[0] > 1 )
      {
        // If Alex is running this, he probably wants to
        // make the output as cluttered as possible (and
        // we all love that about him... right?).
        //
        if (ALEX_MODE) dump_splash_header(Graph,target_seq,*exon_set_name_id,HitExonCoords,ofp,textw);

        print_spliced_alignment(AD,target_seq,HitExonCoords,*exon_set_name_id,ofp,textw);


        // You thought Alex was done making a mess of things?! HA!
        if (ALEX_MODE) DumpSplashFooter(ofp,textw);
      }

      free(HitExonCoords);

    
  }

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: run_model_on_exon_sets
 *
 *  Desc. :  This is the culmination of all of our work!
 *           Given a splice graph, use the spliced nucleotide coordinates for each
 *           connected component (ideally there will only be one...) to pull a final
 *           nucleotide sequence and re-run the model on that (spliced) sequence.
 *
 *  Inputs:  1.         Graph : The final SPLICE_GRAPH struct built on the set of unspliced hits.
 *           2. TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *           3.         gcode : An ESL_GENCODE struct (mainly used for translation).
 *           4.            go : An ESL_GETOPTS struct (in case of options... duh).
 *           5.           ofp : An open file pointer (specificially, the target for output).
 *           6.         textw : The desired line width for output.
 *
 *  Output:  Nothing is returned.  The results of running the model against the spliced
 *           nucleotide sequence are printed to ofp.
 *
 */
void run_model_on_exon_sets
(
  SPLICE_GRAPH * Graph, 
  ESL_SQ       *target_seq,
  ESL_GENCODE  * gcode, 
  ESL_GETOPTS  * go,
  FILE         * ofp,
  int            textw
)
{

  if (DEBUGGING1) DEBUG_OUT("Starting 'run_model_on_exon_sets'",1);

  int i;
  int num_coord_sets;
  int ** ExonCoordSets = GetSplicedExonCoordSets(Graph,&num_coord_sets);


  if (DEBUGGING2) DumpExonSets(ExonCoordSets,num_coord_sets,target_seq,gcode, Graph->revcomp);


  // It's better to be re-using these than destroying
  // and re-allocating* NuclSeq           = NULL; 
  ESL_SQ      * AminoSeq          = NULL; 
  P7_PIPELINE * ExonSetPipeline   = NULL;
  P7_TOPHITS  * ExonSetTopHits    = NULL;
  P7_OPROFILE * ExonSetOModel     = p7_oprofile_Clone(Graph->OModel);
  P7_BG       * ExonSetBackground = p7_bg_Create(Graph->OModel->abc);


  int coord_set_id;
  int exon_set_id = 0;
  for (coord_set_id = 0; coord_set_id < num_coord_sets; coord_set_id++) {

    // If there's only one exon in this set of exons, we'll
    // skip reporting it (maybe have this be a user option?)
    //
    // UNLESS: This is a (more or less) full-model hit
    //
    if (ExonCoordSets[coord_set_id][0] <= 1
        && (ExonCoordSets[coord_set_id][2] > 6
            || ExonCoordSets[coord_set_id][4] <= Graph->Model->M - 6))
      continue;

    // Grab the nucleotides for this set of exons and
    // convert them to a textized ESL_SQ
    //
    int coding_region_len;
    ESL_DSQ *ExonSetNucls =  grab_exon_coord_set_nucls(ExonCoordSets[coord_set_id], target_seq, Graph->revcomp, &coding_region_len);
    ESL_SQ *NuclSeq      = esl_sq_CreateDigitalFrom(target_seq->abc,"Exon Set",ExonSetNucls,(int64_t)coding_region_len,NULL,NULL,NULL);
    NuclSeq->idx = coord_set_id+1;
    esl_sq_Textize(NuclSeq);


    // Translate the nucleotide sequence and convert to
    // a digitized ESL_SQ
    //
    int trans_len = coding_region_len / 3;
    ESL_DSQ * ExonSetTrans = TranslateExonSetNucls(ExonSetNucls,coding_region_len,gcode);
    AminoSeq     = esl_sq_CreateDigitalFrom(Graph->OModel->abc,target_seq->name,ExonSetTrans,(int64_t)trans_len,NULL,NULL,NULL);
    AminoSeq->idx = coord_set_id+1;
    strcpy(AminoSeq->orfid,"exon");
    
    // Prep a P7_PIPELINE datastructure and all of the other
    // friends that we need in order to produce our full evaluation
    // of the set of exons as a coding region for this protein model
    //
    // NOTE: The default 'p7_pipeline_Create' in BATH is causing
    //       option-based errors, and 'p7_pipeline_fs_Create' requires
    //       a significant number of pieces of data that I don't have.
    //
    //       This function is a *heavily* simplified version of the
    //       standard 'p7_pipeline_Create' function.
    //
    P7_PIPELINE * ExonSetPipeline  = p7_pipeline_splash_Create(go,Graph->OModel->M,coding_region_len,FALSE,p7_SEARCH_SEQS);
    ExonSetPipeline->is_translated = TRUE;
    ExonSetPipeline->strands       = p7_STRAND_TOPONLY;
    ExonSetPipeline->block_length  = coding_region_len;


    int pipeline_create_err = p7_pli_NewModel(ExonSetPipeline,ExonSetOModel,ExonSetBackground);
    if (pipeline_create_err == eslEINVAL) 
      p7_Fail(ExonSetPipeline->errbuf);

    p7_pli_NewSeq(ExonSetPipeline,AminoSeq);
    p7_bg_SetLength(ExonSetBackground,AminoSeq->n);
    p7_oprofile_ReconfigLength(ExonSetOModel,AminoSeq->n);



    // Create a P7_TOPHITS datastructure to capture the results
    // for this set of exons and run that rowdy ol' p7_Pipeline!
    
    P7_TOPHITS * ExonSetTopHits = p7_tophits_Create();
    int pipeline_execute_err  = p7_Pipeline(ExonSetPipeline,ExonSetOModel,ExonSetBackground,AminoSeq,NuclSeq,ExonSetTopHits,NULL);
    if (pipeline_execute_err != eslOK) {
      fprintf(stderr,"\n  * PIPELINE FAILURE DURING EXON SET RE-ALIGNMENT *\n\n");
      exit(101);
    }

    // If we were successful in our search, report the hit(s)
    // we produced!
    //
    if (ExonSetTopHits->N) 
      report_spliced_tophits(Graph,ExonSetTopHits,ExonSetPipeline,target_seq,ExonCoordSets[coord_set_id],&exon_set_id,ofp,textw);



    // DESTRUCTION AND REBIRTH!
    free(ExonSetNucls);
    free(ExonSetTrans);
    esl_sq_Destroy(NuclSeq);
    esl_sq_Destroy(AminoSeq);
    p7_tophits_Destroy(ExonSetTopHits);
    p7_pipeline_Destroy(ExonSetPipeline);
    free(ExonCoordSets[coord_set_id]);
    ExonCoordSets[coord_set_id] = NULL;
  }


  // That's all we needed!  Good work, team!
  if (ExonSetTopHits) p7_tophits_Destroy(ExonSetTopHits);
  if (ExonSetPipeline) p7_pipeline_Destroy(ExonSetPipeline);
  p7_oprofile_Destroy(ExonSetOModel);
  p7_bg_Destroy(ExonSetBackground);
  for(i = 0; i < num_coord_sets; i++)
    if(ExonCoordSets[i] != NULL) free(ExonCoordSets[i]);
  free(ExonCoordSets);
  if (DEBUGGING1) DEBUG_OUT("'run_model_on_exon_sets' Complete",-1);

}












/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: SpliceHits
 *
 *  Desc. :  This is the top-level function for spliced hmmsearcht
 *
 *  Inputs:  1.        TopHits : A set of hits generated by hmmsearcht.
 *                               These serve as the initial collection of likely exons that
 *                               we try to find ways to stitch together into full-model hits.
 *           2. GenomicSeqFile : An ESL_SQFILE struct holding the genomic file that contains
 *                               all of the target sequences for search.
 *           3.             gm : The straightforward profile for the protein / family.

 *           6.             go : An ESL_GETOPTS struct (in case of options... duh).
 *           7.            ofp : An open file pointer (specificially, the target for output).
 *           8.          textw : The desired line width for output.
 *
c *  Output:  Nothing is returned, but if we find a way to splice together some of the
 *           input hits they're written out to the ofp.
 *
 */
void SpliceHits
(
  P7_TOPHITS  * TopHits,
  ESL_SQFILE  * GenomicSeqFile,
  P7_PROFILE  * gm,
  P7_OPROFILE * om,
  ESL_GENCODE * gcode,
  ESL_GETOPTS * go,
  FILE        * ofp,
  int           textw
)
{


  int               target_set_id;
  TARGET_RANGE_SET *target_range_set;
  TARGET_RANGE     *curr_target_range;
  ESL_SQ           *target_seq;
  SPLICE_GRAPH     *Graph;
//  P7_TOPHITS       *full_path_hits;
//  SPLICE_GRAPH     *full_path_graph;
//  SPLICE_NODE      *full_path_node;
 
  if (DEBUGGING1) DEBUG_OUT("Starting 'SpliceHits'",1);

    INIT_SECONDS = time(NULL);

  // Start the timer!  You're on the clock, splash!
  //
  ESL_STOPWATCH * Timer = esl_stopwatch_Create();
  esl_stopwatch_Start(Timer);


  // Very first thing we want to do is make sure that our hits are
  // organized by the 'Seqidx' (the target genomic sequence) and position
  // within that file.
  //
  // NOTE: I still need to ensure I'm taking advantage of this sorting!
  //
  p7_tophits_SortBySeqidxAndAlipos(TopHits);
  
  // We'll iterate over search regions (max-2MB ranges of chromosomes)
  
  target_range_set = build_target_ranges(TopHits);

  for (target_set_id = 0; target_set_id < target_range_set->N; target_set_id++) {
  
   //printf("Target %s strand %c range %d to %d\n", target_range_set->target_range[target_set_id]->seqname, (target_range_set->target_range[target_set_id]->complementarity ? '-' : '+'), target_range_set->target_range[target_set_id]->start, target_range_set->target_range[target_set_id]->end); 
  
     curr_target_range = target_range_set->target_range[target_set_id]; 

    // Given that our hits are organized by target sequence, we can
    // be a bit more efficient in our file reading by only pulling
    // target sequences as they change (wrt the upstream hit)  
    target_seq = get_target_range_sequence(GenomicSeqFile, curr_target_range); 

    // This function encapsulates a *ton* of the work we do.
    // In short, take the collection of unspliced hits and build
    // a splice graph representing all (reasonable) ways of splicing
    // them.
    //
    Graph = build_splice_graph(curr_target_range, target_seq, TopHits,gm,om,gcode);
  

    // Evaluate the graph for any holes (or possible holes)
    // worth plugging
    //

    add_missing_exons(Graph,target_seq,gcode);

    // If we're debugging, it might be useful to get a quick
    // picture of what our final splice graph looks like.
    //
    if (DEBUGGING2) DumpGraph(Graph);


    FindBestFullPath(Graph); 
    // Re-run the model on the extracted nucleotide sequence(s)
    // for each connected component of the graph.
    run_model_on_exon_sets(Graph,target_seq,gcode,go,ofp,textw);


    // CLEANUP
    SPLICE_GRAPH_Destroy(Graph);
    esl_sq_Destroy(target_seq);
  }

  // return tophits to sortkey sorting 
  p7_tophits_SortBySortkey(TopHits);

  // More cleanup!
  //TARGET_SET_Destroy(TargetSet);
  Target_Range_Set_Destroy(target_range_set);

  if (DEBUGGING1) {
    DEBUG_OUT("'SpliceHits' Complete",-1);
    fprintf(stderr,"\n\n");
  }



  // Very last thing we'll do -- how long were we working on
  // the splash-related stuff?
  esl_stopwatch_Stop(Timer);
  if (ALEX_MODE) {
    fprintf(stderr,"  Time spent splashing\n  : ");
    esl_stopwatch_Display(stderr,Timer,NULL);
    fprintf(stderr,"\n\n");
  }
  esl_stopwatch_Destroy(Timer);


}


