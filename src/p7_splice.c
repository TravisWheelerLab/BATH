/* P7_SPLICE: implementation of ranked list of top-scoring hits
 *
 * Contents:
 *    1. The P7_TOPHITS object.
 *    2. Standard (human-readable) output of pipeline results.
 *    3. Tabular (parsable) output of pipeline results.
 *    4. Benchmark driver.
 *    5. Test driver.
 */
#include "p7_config.h"

#include <string.h>

#include "easel.h"
#include "hmmer.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                           BEGIN SPLICING STUFF
 * 
 *  Catalogue of Ships
 *  ==================
 *
 *  Fleet 0: Structs
 *
 *  + TARGET_SEQ        : The chunk of genomic sequence we're searching around in
 *  + DOMAIN_OVERLAP    : Information related to how two HMMER hits are spliced together (effectively a splice edge)
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
 *  + GetTargetNuclSeq      : Extract the nucleotides corresponding to the range chosen by 'SetTargetSeqRange'
 *  + GrabNuclRange         : Get the nucleotide sequence of a window within the range extracted by 'GetTargetNuclSeq'
 *  + DetermineNuclType     : [REDUNDANT WITH EASEL] Determine if a nucleotide sequence is DNA or RNA (defaults to DNA)
 *  
 *
 *
 *  Fleet 2: Determining How to Splice Pairs of Hits
 *
 *  + SelectSpliceOpt              : Given a model position and nucleotides, determine an optimal splice junction
 *  + GetContestedUpstreamNucls    : Pull a nucleotide sequence corresponding to a candidate 5' splice site
 *  + GetContestedDownstreamNucls  : Pull a nucleotide sequence corresponding to a candidate 3' splice site
 *  + AminoScoreAtPosition         : ASSUMING THAT WE'RE MOVING THROUGH THE MODEL, determine a positional score
 *  + FindOptimalSpliceSite        : Examine a DOMAIN_OVERLAP to determine an optimal splicing of two P7_HITs
 *  + SpliceOverlappingDomains     : Do a bit of bookkeeping around 'FindOptimalSpliceSite'
 *  + GetNuclRangesFromAminoCoords : Given two P7_HITs that we're considering splicing, determine their nucleotide coordinate ranges
 *  + SketchSpliceEdge             : Perform any necessary extension to make two P7_HITs splice-able and call 'SpliceOverlappingDomains'
 *  + HitsAreSpliceComaptible      : Determine whether two hits are positioned correctly in the model/target to be spliced
 *  + ExcessiveGapContent          : Determine whether a hit is too gappy to be a strong initial candidate for splicing (can be recovered in sub-model search)
 *  + OutsideSearchArea            : Determine whether a hit is outside of the area selected in 'SetTargetSeqRange'
 *  + GatherViableSpliceEdges      : Given a P7_TOPHITS, identify all pairs of hits that might be splice-able
 *
 *
 *
 *  Fleet 3: Producing a Splice Graph According to the Spliced Hit Pairs
 *
 *  + InitSpliceNode        : Allocate and initialize a node in our splice graph
 *  + ConnectNodesByEdge    : Given two nodes implicated in a satisfactorily spliced DOMAIN_OVERLAP, draw an edge in the splice graph
 *  + GatherCTermNodes      : Derive a list of all nodes in the graph whose corresponding hits terminate at the C-terminal end of the model
 *  + EdgeWouldEraseNode    : Check if the best path through a node would (effectively) skip that node
 *  + PullUpCumulativeScore : Recursively determine the best cumulative score leading up to each node in the graph
 *  + EvaluatePaths         : Prepare for and call 'PullUpCumulativeScore' so we can determine the best path through the graph
 *  + FillOutGraphStructure : Initialize the splice graph and call 'EvaluatePaths' to produce our complete splice graph
 *  + FindBestFullPath      : Given a complete splice graph produced by 'FillOutGraphStructure' find an optimal path from model positions 1..M
 *  + BuildSpliceGraph      : Top-level function for splice graph construction.  Calls 'GatherViableSpliceEdges', 'FillOutGraphStructure', and 'FindBestFullPath'
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
 *  FLAGSHIP: SpliceHits : Calls 'GetTargetNuclSeq' to acquire a general coding region on the nucleotide target
 *                       : Calls 'BuildSpliceGraph' to build a first pass of the splice graph using the results from translated HMMER
 *                       : Calls 'AddMissingExonsToGraph' to identify and patch any potential holes in the splice graph (missing exons)
 *                       : Calls 'RunModelOnExonSets' to generate the final spliced alignments of sets of exons to the model
 *
 *
 */



// Before we get to the fun stuff, let's just set up some
// bureaucratic stuff to make debugging relatively (hopefully)
// painless
static int ALEX_MODE = 1; // Print some extra metadata around hits
static int DEBUGGING = 0; // Print debugging output?


// Ever want to know what function you're in, and how deep it
// is (roughly)? Well, wonder no more!
int FUNCTION_DEPTH = 0;
void DEBUG_OUT (const char * message, const int func_depth_change) {

  if (func_depth_change > 0) 
    FUNCTION_DEPTH += func_depth_change;
  
  fprintf(stderr,"  SplDebug:");
  int debug_depth;
  for (debug_depth=0; debug_depth<FUNCTION_DEPTH; debug_depth++) 
    fprintf(stderr,"  ");
  fprintf(stderr,"%s\n",message);
  fflush(stderr);
  
  if (func_depth_change < 0) 
    FUNCTION_DEPTH += func_depth_change;

}


static float SSSCORE[2]      = {-0.7,0.0}; // Non-canon vs canon splice site
static float EDGE_FAIL_SCORE = -14773.0;   // Makes me thirsty for a latte!

static char  AMINO_CHARS[21] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','-'};
static char    DNA_CHARS[ 6] = {'A','C','G','T','-','N'};
static char LC_DNA_CHARS[ 6] = {'a','c','g','t','-','n'};
static char    RNA_CHARS[ 6] = {'A','C','G','U','-','N'};
static char LC_RNA_CHARS[ 6] = {'a','c','g','u','-','n'};


// How many chromosome ranges are we willing to try splicing into?
static int MAX_TARGET_REGIONS = 5;


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




//////////////////////////////////////////////////////////////////////////////////////////////
//
//  Struct: TARGET_SET
//
//  Desc. : 
//
typedef struct _target_set {

  int num_target_seqs;

  char    ** TargetSeqNames; // Also (individually) borrowed pointers!
  int64_t  * TargetStarts;
  int64_t  * TargetEnds;
  int      * TargetIsRevcomp;

} TARGET_SET;




//////////////////////////////////////////////////////////////////////////////////////////////
//
//  Struct: TARGET_SEQ
//
//  Desc. : 
//
typedef struct _target_seq {
  
  ESL_SQ  * esl_sq; // For freeing, we'll need this pointer
  ESL_DSQ * Seq;
  const ESL_ALPHABET * abc;

  // What's the chromosome (or, more generally, sequence name?)
  // NOTE that this is a borrowed pointer, so we don't want to free it!
  char * SeqName;


  // NOTE that start is always less than end (even if revcomp)
  int64_t start;
  int64_t end;


  // This is primarily used to make sure that 'in bounds' also
  // means the right complementarity
  int revcomp;


} TARGET_SEQ;




//////////////////////////////////////////////////////////////////////////////////////////////
//
//  Struct: DOMAIN_OVERLAP
//
//  Desc. :
//
typedef struct _domain_overlap {

  const ESL_ALPHABET * ntalpha;


  int amino_start;
  int amino_end;


  int upstream_hit_id;
  int upstream_dom_id;
  int upstream_nucl_start;
  int upstream_nucl_end;
  int upstream_ext_len;
  int upstream_disp_start;

  P7_TOPHITS    * UpstreamTopHits;
  P7_ALIDISPLAY * UpstreamDisplay;
  ESL_DSQ       * UpstreamNucls;
  

  int downstream_hit_id;
  int downstream_dom_id;
  int downstream_nucl_start;
  int downstream_nucl_end;
  int downstream_ext_len;
  int downstream_disp_end;

  P7_TOPHITS    * DownstreamTopHits;
  P7_ALIDISPLAY * DownstreamDisplay;
  ESL_DSQ       * DownstreamNucls;


  int   upstream_exon_terminus;
  int downstream_exon_terminus;

  int   upstream_spliced_nucl_end;
  int downstream_spliced_nucl_start;


  float score_density;
  float score;

} DOMAIN_OVERLAP;







//////////////////////////////////////////////////////////////////////////////////////////////
//
//  Struct: SPLICE_NODE
//
//  Desc. :
//
typedef struct _splice_node {


  int node_id;


  int was_missed;
  int hit_id;
  int dom_id;


  int in_edge_cap;
  int num_in_edges;
  int best_in_edge; // Relative to the following arrays
  DOMAIN_OVERLAP ** InEdges;
  struct _splice_node ** UpstreamNodes;


  int out_edge_cap;
  int num_out_edges;
  DOMAIN_OVERLAP ** OutEdges;
  struct _splice_node ** DownstreamNodes;


  int is_n_terminal;
  int is_c_terminal;


  float hit_score;
  float cumulative_score;


} SPLICE_NODE;







//////////////////////////////////////////////////////////////////////////////////////////////
//
//  Struct: SPLICE_GRAPH
//
//  Desc. :
//
typedef struct _splice_graph {


  P7_TOPHITS  *    TopHits;
  P7_TOPHITS  * MissedHits; // This is a fake P7_TOPHITS!  Don't be fooled!

  P7_PROFILE  *  Model;
  P7_OPROFILE * OModel;

  int revcomp;


  // Because we build based on the DOMAIN_OVERLAP datastructure,
  // it's helpful to be able to find each node ID by way of
  // these lookups.
  int ** TH_HitDomToNodeID; //    TopHits
  int  * MH_HitToNodeID;    // MissedHits


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
  fprintf(stderr,"   |     |  Source in P7_TOPHITS  :  Hit %d, Domain %d\n", Node->hit_id, Node->dom_id);
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
  free(TS->TargetStarts);
  free(TS->TargetEnds);
  free(TS->TargetIsRevcomp);
  free(TS);
}





/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function:  TARGET_SEQ_Destroy
 *
 */
void TARGET_SEQ_Destroy
(TARGET_SEQ * TS)
{
  esl_sq_Destroy(TS->esl_sq); // This takes care of 'TS->Seq' too
  free(TS);
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function:  DOMAIN_OVERLAP_Destroy
 *
 */
void DOMAIN_OVERLAP_Destroy
(DOMAIN_OVERLAP * DO)
{
  free(DO->UpstreamNucls);
  free(DO->DownstreamNucls);
  free(DO);
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function:  SPLICE_NODE_Destroy
 *
 */
void SPLICE_NODE_Destroy
(SPLICE_NODE * Node)
{

  int in_edge_id;
  for (in_edge_id = 0; in_edge_id < Node->num_in_edges; in_edge_id++) {
    if (Node->InEdges[in_edge_id]) {
      DOMAIN_OVERLAP_Destroy(Node->InEdges[in_edge_id]);
      Node->InEdges[in_edge_id] = NULL;
    }
  }

  int out_edge_id;
  for (out_edge_id = 0; out_edge_id < Node->num_out_edges; out_edge_id++) {
    if (Node->OutEdges[out_edge_id]) {
      DOMAIN_OVERLAP_Destroy(Node->OutEdges[out_edge_id]);
      Node->OutEdges[out_edge_id] = NULL;
    }
  }

  free(Node->InEdges);
  free(Node->OutEdges);
  free(Node);
  Node = NULL;

}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function:  SPLICE_GRAPH_Destroy
 *
 */
void SPLICE_GRAPH_Destroy
(SPLICE_GRAPH * Graph)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'SPLICE_GRAPH_Destroy'",1);


  // Because we don't want to accidentally double-free any
  // of the DOMAIN_OVERLAP structs, we'll put each node in
  // charge of its incoming DOs.
  if (Graph->Nodes) {
    
    int node_id, in_edge_index;
    SPLICE_NODE * Node;
    for (node_id=1; node_id<=Graph->num_nodes; node_id++) {
      Node = Graph->Nodes[node_id];
      free(Node->DownstreamNodes);
      free(Node->UpstreamNodes);
      for (in_edge_index = 0; in_edge_index < Node->num_in_edges; in_edge_index++) 
        DOMAIN_OVERLAP_Destroy(Node->InEdges[in_edge_index]);
     
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
  
    for (hit_id=0; hit_id<Graph->TopHits->N; hit_id++)
      free(Graph->TH_HitDomToNodeID[hit_id]);
    free(Graph->TH_HitDomToNodeID);
  
    if (Graph->MH_HitToNodeID) 
      free(Graph->MH_HitToNodeID);
  
  }


  // We need to manually destroy the fake 'MissedHits'
  if (Graph->MissedHits) {

    for (hit_id=0; hit_id<Graph->MissedHits->N; hit_id++) {

      // Luckily, the AD is for real (generally) so we should
      // be able to use the generic destructor
      p7_alidisplay_Destroy(Graph->MissedHits->hit[hit_id]->dcl->ad);
      free(Graph->MissedHits->hit[hit_id]->dcl);
      free(Graph->MissedHits->hit[hit_id]);

    }
    free(Graph->MissedHits);

  }


  free(Graph);
  Graph = NULL;


  if (DEBUGGING) DEBUG_OUT("'SPLICE_GRAPH_Destroy' Complete",-1);

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
  if (DEBUGGING) DEBUG_OUT("Starting 'FloatHighLowSortIndex'",1);

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

  if (DEBUGGING) DEBUG_OUT("'FloatHighLowSortIndex' Complete",-1);

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
  if (DEBUGGING) DEBUG_OUT("Starting 'FloatLowHighSortIndex'",1);
  int * SortIndex = FloatHighLowSortIndex(Vals,num_vals);
  int i;
  for (i=0; i<num_vals/2; i++) {
    int tmp = SortIndex[i];
    SortIndex[i] = SortIndex[(num_vals-1)-i];
    SortIndex[(num_vals-1)-i] = tmp;
  }
  if (DEBUGGING) DEBUG_OUT("'FloatLowHighSortIndex' Complete",-1);
  return SortIndex;
}











/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetMinAndMaxCoords
 *
 *  Desc. :  Iterate over a collection of hits (from non-spliced hmmsearcht)
 *           and identify the minimum and maximum coordinates of hits to the
 *           genomic sequence.
 *
 *  Inputs:  1. TopHits       :
 *           2. TargetNuclSeq : The sub-sequence of the target sequence wherein all hits reside.
 *
 *  Output:
 *
 */
void GetMinAndMaxCoords
(P7_TOPHITS * TopHits, TARGET_SEQ * TargetNuclSeq)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'GetMinAndMaxCoords'",1);

  // First, let's figure out if we're revcomp
  int revcomp = 0;
  int64_t min,max;

  P7_HIT    * Hit;
  P7_DOMAIN * Dom;
  Dom = &(TopHits->hit[0]->dcl[0]);
  if (Dom->ad->sqfrom > Dom->ad->sqto) {

    min = Dom->ad->sqto;
    max = Dom->ad->sqfrom;
    revcomp = 1;

  } else {

    min = Dom->ad->sqfrom;
    max = Dom->ad->sqto;

  }


  int hit_id;
  int dom_id;
  if (revcomp) {

    for (hit_id = 0; hit_id < (int)(TopHits->N); hit_id++) {

      Hit = TopHits->hit[hit_id];
      for (dom_id = 0; dom_id < Hit->ndom; dom_id++) {

        Dom = &(Hit->dcl[dom_id]);

        if (Dom->ad->sqfrom > max)
          max = Dom->ad->sqfrom;

        if (Dom->ad->sqto < min)
          min = Dom->ad->sqto;

      }

    }

  } else {

    for (hit_id = 0; hit_id < (int)(TopHits->N); hit_id++) {

      Hit = TopHits->hit[hit_id];
      for (dom_id = 0; dom_id < Hit->ndom; dom_id++) {

        Dom = &(Hit->dcl[dom_id]);

        if (Dom->ad->sqto > max)
          max = Dom->ad->sqto;

        if (Dom->ad->sqfrom < min)
          min = Dom->ad->sqfrom;

      }

    }

  }

  TargetNuclSeq->start = min;
  TargetNuclSeq->end   = max;

  if (DEBUGGING) DEBUG_OUT("'GetMinAndMaxCoords' Complete",-1);

}













/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: SelectTargetRanges
 *
 *  Desc.: Scan 'TopHits' to determine a bounded search range that we're
 *         inclined to consider the coding region for the input protein.
 *
 */
TARGET_SET * SelectTargetRanges
(P7_TOPHITS * TopHits)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'SelectTargetRanges'",1);


  int hit_id, dom_id, i;


  // Because we go through
  int num_hits = (int)(TopHits->N);
  float * HitScores = malloc(num_hits * sizeof(float));
  for (hit_id = 0; hit_id < num_hits; hit_id++) {
    dom_id = TopHits->hit[hit_id]->best_domain;
    HitScores[hit_id] = (&TopHits->hit[hit_id]->dcl[dom_id])->envsc;
  }
  int * HitScoreSort = FloatHighLowSortIndex(HitScores,num_hits);


  // Initialize our TARGET_SETS struct
  TARGET_SET * TargetSet     = malloc(sizeof(TARGET_SET));
  TargetSet->TargetSeqNames  = malloc(MAX_TARGET_REGIONS * sizeof(char *));
  TargetSet->TargetStarts    = malloc(MAX_TARGET_REGIONS * sizeof(int64_t));
  TargetSet->TargetEnds      = malloc(MAX_TARGET_REGIONS * sizeof(int64_t));
  TargetSet->TargetIsRevcomp = malloc(MAX_TARGET_REGIONS * sizeof(int));



  // We'll need to check each target sequence's length to make sure we don't
  // let our range exceed the length of the chromosome.



  int num_targets = 0;

  int sort_id;
  for (sort_id = 0; sort_id < num_hits; sort_id++) {


    int base_hit_id = HitScoreSort[sort_id];

    char * CandidateSeqName     = TopHits->hit[base_hit_id]->name;
    int candidate_best_dom      = TopHits->hit[base_hit_id]->best_domain;
    P7_ALIDISPLAY * CandidateAD = (&TopHits->hit[base_hit_id]->dcl[candidate_best_dom])->ad;

    int candidate_revcomp = 0;
    if (CandidateAD->sqfrom > CandidateAD->sqto)
      candidate_revcomp = 1;


    int new_target = 1;
    for (i=0; i<num_targets; i++) {
      if (!strcmp(CandidateSeqName,TargetSet->TargetSeqNames[i]) && candidate_revcomp == TargetSet->TargetIsRevcomp[i]) {
        new_target = 0;
        break;
      }
    }

    if (!new_target) 
      continue;


    // Now that we know this is a new chromosome, define a coding region!
    // First, what we'll do is define a hard maximum search area defined
    // by the highest-scoring hit to this sequence.
    int64_t min_coord, max_coord;
    if (candidate_revcomp) {
      min_coord = CandidateAD->sqto;
      max_coord = CandidateAD->sqfrom;
    } else {
      min_coord = CandidateAD->sqfrom;
      max_coord = CandidateAD->sqto;
    }

    int64_t min_cap = min_coord - 1000000;
    int64_t max_cap = max_coord + 1000000;


    // The min_cap and max_cap now define the absolute furthest we're
    // willing to go for our search region, but ideally we can shrink
    // down to a much tighter zone
    for (hit_id = 0; hit_id < num_hits; hit_id++) {


      if (strcmp(CandidateSeqName,TopHits->hit[hit_id]->name))
        continue;


      P7_ALIDISPLAY * AD = (&TopHits->hit[hit_id]->dcl[0])->ad; 

      if ( candidate_revcomp && AD->sqfrom < AD->sqto) continue;
      if (!candidate_revcomp && AD->sqfrom > AD->sqto) continue;


      for (dom_id = 0; dom_id < TopHits->hit[hit_id]->ndom; dom_id++) {
      
          
        AD = (&TopHits->hit[hit_id]->dcl[dom_id])->ad;

        // New minimum?
        // We could revcomp check, but I don't know if it's any faster...
        if (AD->sqto < min_coord && AD->sqto > min_cap)
          min_coord = AD->sqto;
          
        if (AD->sqfrom < min_coord && AD->sqfrom > min_cap)
          min_coord = AD->sqfrom;
          

        // New maximum?
        if (AD->sqto > max_coord && AD->sqto < max_cap)
          max_coord = AD->sqto;
          
        if (AD->sqfrom > max_coord && AD->sqfrom < max_cap)
          max_coord = AD->sqfrom;

      }

    }


    TargetSet->TargetSeqNames[num_targets]  = CandidateSeqName;
    TargetSet->TargetStarts[num_targets]    = min_coord;
    TargetSet->TargetEnds[num_targets]      = max_coord;
    TargetSet->TargetIsRevcomp[num_targets] = candidate_revcomp;

    num_targets++;

    if (num_targets == MAX_TARGET_REGIONS)
      break;


  }


  TargetSet->num_target_seqs = num_targets;


  if (DEBUGGING) DEBUG_OUT("'SelectTargetRanges' Complete",-1);

  free(HitScores);
  free(HitScoreSort);

  return TargetSet;


}












/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetTargetNuclSeq
 *
 *  Desc. :  Given a sequence file and a collection of translated search hits, extract a
 *           sub-region of the genome that contains the all of the nucleotides implicated in
 *           the set of hits.
 *
 *           This assumes (correctly!) that all of the hits are to the same nucleotide sequence.
 *
 *  Inputs:  1. GenomicSeqFile : An ESL_SQFILE struct representing the input genom(e/ic sequence) file.
 *           2.        TopHits : A collection of hits achieved by "standard" hmmsearcht (unspliced)
 *
 *  Output:  A TARGET_SEQ struct, containing the nucleotide subsequence within which all of the
 *           unspliced hmmsearcht hits reside.
 *
 */
TARGET_SEQ * GetTargetNuclSeq
(ESL_SQFILE * GenomicSeqFile, TARGET_SET * TargetSet, int target_set_id)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'GetTargetNuclSeq'",1);


  TARGET_SEQ * TargetNuclSeq = (TARGET_SEQ *)malloc(sizeof(TARGET_SEQ));


  TargetNuclSeq->abc     = GenomicSeqFile->abc;
  TargetNuclSeq->SeqName = TargetSet->TargetSeqNames[target_set_id];
  TargetNuclSeq->start   = TargetSet->TargetStarts[target_set_id];
  TargetNuclSeq->end     = TargetSet->TargetEnds[target_set_id];
  TargetNuclSeq->revcomp = TargetSet->TargetIsRevcomp[target_set_id];


  ESL_SQFILE * TmpSeqFile;
  esl_sqfile_Open(GenomicSeqFile->filename,GenomicSeqFile->format,NULL,&TmpSeqFile);
  esl_sqfile_OpenSSI(TmpSeqFile,NULL);

  ESL_SQ * SeqInfo = esl_sq_Create();
  esl_sqio_FetchInfo(TmpSeqFile,TargetNuclSeq->SeqName,SeqInfo);


  // In case there's a terminal search region we need to consider,
  // pull in a bit of extra sequence
  TargetNuclSeq->start -= MAX_SUB_NUCL_RANGE;
  if (TargetNuclSeq->start < 1)
    TargetNuclSeq->start = 1;
  
  TargetNuclSeq->end += MAX_SUB_NUCL_RANGE;
  if (TargetNuclSeq->end > SeqInfo->L)
    TargetNuclSeq->end = SeqInfo->L;


  TargetNuclSeq->esl_sq = esl_sq_CreateDigital(TargetNuclSeq->abc);
  int fetch_err_code    = esl_sqio_FetchSubseq(TmpSeqFile,TargetNuclSeq->SeqName,TargetNuclSeq->start,TargetNuclSeq->end,TargetNuclSeq->esl_sq);


  esl_sqfile_Close(TmpSeqFile);
  esl_sq_Destroy(SeqInfo);


  if (fetch_err_code != eslOK) {
    fprintf(stderr,"\n  ERROR: Failed to fetch target subsequence (is there an .ssi index for the sequence file?)\n");
    fprintf(stderr,"         Requested search area: %s:%ld..%ld\n\n",TargetNuclSeq->SeqName,TargetNuclSeq->start,TargetNuclSeq->end);
    exit(1);
  }


  TargetNuclSeq->Seq = TargetNuclSeq->esl_sq->dsq;


  if (DEBUGGING) DEBUG_OUT("'GetTargetNuclSeq' Complete",-1);


  return TargetNuclSeq;


}











/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GrabNuclRange
 *
 *  Desc. :  Extract a sub-sequence from our (already a sub-)sequence stored
 *           in a TARGET_SEQ datastructure.
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
ESL_DSQ * GrabNuclRange
(TARGET_SEQ * TargetNuclSeq, int start, int end)
{

  int len = abs(end - start) + 1;

  char * Seq = malloc((len+1) * sizeof(char));

  // Keep in mind that DSQs are [1..n]
  int read_index = start - (int)(TargetNuclSeq->start) + 1;

  if (start < end) {
  
    int i;
    for (i=0; i<len; i++)
      Seq[i] = DNA_CHARS[TargetNuclSeq->Seq[read_index++]];

  } else {

    int i;
    for (i=0; i<len; i++) {
      int fwd_nucl_code = TargetNuclSeq->Seq[read_index--];
      if (fwd_nucl_code < 4) Seq[i] = DNA_CHARS[3 - fwd_nucl_code];
      else                   Seq[i] ='N';
    }

  }
  Seq[len] = 0;


  ESL_DSQ * NuclSubseq;
  esl_abc_CreateDsq(TargetNuclSeq->abc,Seq,&NuclSubseq);


  free(Seq);


  return NuclSubseq;

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: DetermineNuclType
 *
 *  Desc. :  Given a P7_ALIDISPLAY, determine the alphabet type of the nucleotide sequence.
 *
 *  Inputs:  1. AliDisplay : A P7_ALIDISPLAY object, assumed to be from hmmsearcht.
 *
 *  Output:  The easel code for the alphabet type of the nucleotide sequence (DNA or RNA).
 *
 */
int DetermineNuclType 
(P7_ALIDISPLAY * AliDisplay)
{
  int i;
  for (i=0; i<strlen(AliDisplay->ntseq); i++) {
    if (AliDisplay->ntseq[i] == 'T') return eslDNA;
    if (AliDisplay->ntseq[i] == 'U') return eslRNA;
  }
  return 0;
}















/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: SelectSpliceOpt
 *
 */
int SelectSpliceOpt
(
  DOMAIN_OVERLAP * Overlap,
  ESL_DSQ * UN,
  ESL_DSQ * DN,
  int model_pos,
  P7_PROFILE * gm,
  ESL_GENCODE * gcode,
  float * splice_score
)
{

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

  }


  free(Codon);

  return best_opt;

}












/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetContestedUpstreamNucls
 *
 */
void GetContestedUpstreamNucls
(ESL_DSQ * NuclSeq, int read_pos, ESL_DSQ * UN)
{
  int write_pos = 1;
  while (write_pos<6) {
    if (NuclSeq[read_pos] >= 0 && NuclSeq[read_pos] <= 3)
      UN[write_pos++] = NuclSeq[read_pos];
    read_pos++;
  }
}




/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetContestedDownstreamNucls
 *
 */
void GetContestedDownstreamNucls
(ESL_DSQ * NuclSeq, int read_pos, ESL_DSQ * DN)
{
  int write_pos = 5;
  while (write_pos) {
    if (NuclSeq[read_pos] >= 0 && NuclSeq[read_pos] <= 3)
      DN[write_pos--] = NuclSeq[read_pos];
    read_pos--;
  }
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


  if (amino_index > 20) {
    if      (*state == p7P_MM) { *state = p7P_MD; }
    else if (*state == p7P_MD) { *state = p7P_DD; }
    // Otherwise, we *should* be staying DD
  } else if (AD && AD->ntseq[3*display_pos] == '-') {
    if      (*state == p7P_MM) { *state = p7P_MI; }
    else if (*state == p7P_MI) { *state = p7P_II; }
    // Otherwise, we *should* be staying II
  } else {
    if      (*state == p7P_MI) { *state = p7P_IM; }
    else if (*state == p7P_II) { *state = p7P_IM; }
    else if (*state == p7P_MD) { *state = p7P_DM; }
    else if (*state == p7P_DD) { *state = p7P_DM; }
    else                       { *state = p7P_MM; }
  }
  float transition_score = p7P_TSC(gm,model_pos,*state);

  
  float emission_score = 0.0;
  if (*state != p7P_MD && *state != p7P_DD) {
    if (*state == p7P_MI || *state == p7P_II)
      emission_score = p7P_ISC(gm,model_pos,amino_index);
    else
      emission_score = p7P_MSC(gm,model_pos,amino_index);
  }


  return transition_score + emission_score;


}







/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: FindOptimalSpliceSite
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
float FindOptimalSpliceSite
(
  DOMAIN_OVERLAP * Overlap,
  P7_PROFILE     * gm,
  ESL_GENCODE    * gcode,
  int *   upstream_splice_index,
  int * downstream_splice_index,
  int * split_amino_model_index,
  int * codon_split_option
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'FindOptimalSpliceSite'",1);



  // We'll want to track whether we're approaching each new
  // position from having been in a match / deletion / insert
  // (yes, there are also other special states, but those
  // states are for nerds)
  int model_state;


  //
  //  UPSTREAM 
  //


  int upstream_nucl_cnt = abs(Overlap->upstream_nucl_start - Overlap->upstream_nucl_end) + 1;


  // Allocate the arrays we'll use to capture the
  // upstream candidates for splice site splitting
  int     us_arr_cap = upstream_nucl_cnt/3 + 10; // IDK
  int   * USTrans    = malloc(us_arr_cap*sizeof(int));
  int   * USModelPos = malloc(us_arr_cap*sizeof(int));
  int   * USNuclPos  = malloc(us_arr_cap*sizeof(int));
  float * USScores   = malloc(us_arr_cap*sizeof(float));


  // We need to know how many gaps we've encountered for the
  // final position calculation
  int * USGaps = malloc(us_arr_cap * sizeof(int));


  // We'll just assume we matched our way in
  model_state = p7P_MM;

  int us_trans_len   = 0; // Without messy indels, comes out to upstream_nucl_cnt/3
  int nucl_read_pos  = 1;
  int model_pos      = Overlap->amino_start;
  int display_pos    = Overlap->upstream_disp_start;
  while (model_pos <= Overlap->UpstreamDisplay->hmmto) {


    USNuclPos[us_trans_len] = nucl_read_pos;


    int amino_index = 27;
    if (Overlap->UpstreamDisplay->aseq[display_pos] != '-') {

      amino_index = esl_gencode_GetTranslation(gcode,&(Overlap->UpstreamNucls[nucl_read_pos]));
      nucl_read_pos += 3;

      USGaps[us_trans_len] = 0;

    } else {

      USGaps[us_trans_len] = 1;

    }


    USTrans[us_trans_len]    = amino_index;
    USModelPos[us_trans_len] = model_pos;
    USScores[us_trans_len]   = AminoScoreAtPosition(gm,amino_index,model_pos,Overlap->UpstreamDisplay,display_pos,&model_state);


    if (Overlap->UpstreamDisplay->aseq[display_pos] != '.')
      model_pos++;
    display_pos++;
    us_trans_len++;


    // NOTE: this probably means significant gappiness
    //       (*barely* within threshold not to be rejected
    //       by ExcessiveGapContent)
    if (us_trans_len == us_arr_cap) {
      us_arr_cap *= 2;
      int   * TmpTrans    = malloc(us_arr_cap*sizeof(int));
      int   * TmpModelPos = malloc(us_arr_cap*sizeof(int));
      int   * TmpNuclPos  = malloc(us_arr_cap*sizeof(int));
      float * TmpScores   = malloc(us_arr_cap*sizeof(float));
      int   * TmpGaps     = malloc(us_arr_cap*sizeof(int));
      int copy_index;
      for (copy_index=0; copy_index<us_trans_len; copy_index++) {
        TmpTrans[copy_index]    = USTrans[copy_index];
        TmpModelPos[copy_index] = USModelPos[copy_index];
        TmpNuclPos[copy_index]  = USNuclPos[copy_index];
        TmpScores[copy_index]   = USScores[copy_index];
        TmpGaps[copy_index]     = USGaps[copy_index];
      }
      free(USTrans);    USTrans    = TmpTrans;
      free(USModelPos); USModelPos = TmpModelPos;
      free(USNuclPos);  USNuclPos  = TmpNuclPos;
      free(USScores);   USScores   = TmpScores;
      free(USGaps);     USGaps     = TmpGaps;
    }


  }
  int us_pre_ext_end_pos = us_trans_len-1; // For scoring purposes




  // Incorporate the extension
  while (model_pos <= Overlap->amino_end) {


    USNuclPos[us_trans_len] = nucl_read_pos;


    int amino_index = esl_gencode_GetTranslation(gcode,&(Overlap->UpstreamNucls[nucl_read_pos]));
    nucl_read_pos += 3;


    // Necessarily not a gap
    USGaps[us_trans_len] = 0;


    USTrans[us_trans_len]    = amino_index;
    USModelPos[us_trans_len] = model_pos;
    USScores[us_trans_len]   = AminoScoreAtPosition(gm,amino_index,model_pos,NULL,0,&model_state);


    model_pos++;
    us_trans_len++;


    // Resize?
    if (us_trans_len == us_arr_cap) {
      us_arr_cap *= 2;
      int   * TmpTrans    = malloc(us_arr_cap*sizeof(int));
      int   * TmpModelPos = malloc(us_arr_cap*sizeof(int));
      int   * TmpNuclPos  = malloc(us_arr_cap*sizeof(int));
      float * TmpScores   = malloc(us_arr_cap*sizeof(float));
      int   * TmpGaps     = malloc(us_arr_cap*sizeof(int));
      int copy_index;
      for (copy_index=0; copy_index<us_trans_len; copy_index++) {
        TmpTrans[copy_index]    = USTrans[copy_index];
        TmpModelPos[copy_index] = USModelPos[copy_index];
        TmpNuclPos[copy_index]  = USNuclPos[copy_index];
        TmpScores[copy_index]   = USScores[copy_index];
        TmpGaps[copy_index]     = USGaps[copy_index];
      }
      free(USTrans);    USTrans    = TmpTrans;
      free(USModelPos); USModelPos = TmpModelPos;
      free(USNuclPos);  USNuclPos  = TmpNuclPos;
      free(USScores);   USScores   = TmpScores;
      free(USGaps);     USGaps     = TmpGaps;
    }

  }



  //
  //  DOWNSTREAM
  //


  int downstream_nucl_cnt = abs(Overlap->downstream_nucl_start - Overlap->downstream_nucl_end) + 1;


  // Time to malloc!
  int     ds_arr_cap = downstream_nucl_cnt/3 + 10; // Still DK
  int   * DSTrans    = malloc(ds_arr_cap*sizeof(int));
  int   * DSModelPos = malloc(ds_arr_cap*sizeof(int));
  int   * DSNuclPos  = malloc(ds_arr_cap*sizeof(int));
  float * DSScores   = malloc(ds_arr_cap*sizeof(float));
  int   * DSGaps     = malloc(ds_arr_cap*sizeof(int));


  model_state = p7P_MM;

  int ds_trans_len = 0;
  nucl_read_pos    = 3; // Because we sneak in 2 extra nucleotides for splice signal
  model_pos        = Overlap->amino_start;
  while (ds_trans_len < Overlap->downstream_ext_len) {


    DSNuclPos[ds_trans_len] = nucl_read_pos+2;


    int amino_index = esl_gencode_GetTranslation(gcode,&(Overlap->DownstreamNucls[nucl_read_pos]));
    nucl_read_pos += 3;


    // Again, extensions cannot be gaps
    DSGaps[ds_trans_len] = 0;


    DSTrans[ds_trans_len]    = amino_index;
    DSModelPos[ds_trans_len] = model_pos;
    DSScores[ds_trans_len]   = AminoScoreAtPosition(gm,amino_index,model_pos,NULL,0,&model_state);


    model_pos++;
    ds_trans_len++;


    // Resize?
    // (Again, this probably means a gappy hit...)
    if (ds_trans_len == ds_arr_cap) {
      ds_arr_cap *= 2;
      int   * TmpTrans    = malloc(ds_arr_cap*sizeof(int));
      int   * TmpModelPos = malloc(ds_arr_cap*sizeof(int));
      int   * TmpNuclPos  = malloc(ds_arr_cap*sizeof(int));
      float * TmpScores   = malloc(ds_arr_cap*sizeof(float));
      int   * TmpGaps     = malloc(ds_arr_cap*sizeof(int));
      int copy_index;
      for (copy_index=0; copy_index<ds_trans_len; copy_index++) {
        TmpTrans[copy_index]    = DSTrans[copy_index];
        TmpModelPos[copy_index] = DSModelPos[copy_index];
        TmpNuclPos[copy_index]  = DSNuclPos[copy_index];
        TmpScores[copy_index]   = DSScores[copy_index];
        TmpGaps[copy_index]     = DSGaps[copy_index];
      }
      free(DSTrans);    DSTrans    = TmpTrans;
      free(DSModelPos); DSModelPos = TmpModelPos;
      free(DSNuclPos);  DSNuclPos  = TmpNuclPos;
      free(DSScores);   DSScores   = TmpScores;
      free(DSGaps);     DSGaps     = TmpGaps;
    }

  }



  display_pos = 0;
  while (model_pos <= Overlap->amino_end) {


    DSNuclPos[ds_trans_len] = nucl_read_pos+2;


    int amino_index = 27;
    if (Overlap->DownstreamDisplay->aseq[display_pos] != '-') {

      amino_index = esl_gencode_GetTranslation(gcode,&(Overlap->DownstreamNucls[nucl_read_pos]));
      nucl_read_pos += 3;

      DSGaps[ds_trans_len] = 0;

    } else {

      DSGaps[ds_trans_len] = 1;

    }


    DSTrans[ds_trans_len]    = amino_index;
    DSModelPos[ds_trans_len] = model_pos;
    DSScores[ds_trans_len]   = AminoScoreAtPosition(gm,amino_index,model_pos,Overlap->DownstreamDisplay,display_pos,&model_state);


    if (Overlap->DownstreamDisplay->aseq[display_pos] != '.')
      model_pos++;
    display_pos++;
    ds_trans_len++;


    // Resize?
    // (Again, this probably means a gappy hit...)
    if (ds_trans_len == ds_arr_cap) {
      ds_arr_cap *= 2;
      int   * TmpTrans    = malloc(ds_arr_cap*sizeof(int));
      int   * TmpModelPos = malloc(ds_arr_cap*sizeof(int));
      int   * TmpNuclPos  = malloc(ds_arr_cap*sizeof(int));
      float * TmpScores   = malloc(ds_arr_cap*sizeof(float));
      int   * TmpGaps     = malloc(ds_arr_cap*sizeof(int));
      int copy_index;
      for (copy_index=0; copy_index<ds_trans_len; copy_index++) {
        TmpTrans[copy_index]    = DSTrans[copy_index];
        TmpModelPos[copy_index] = DSModelPos[copy_index];
        TmpNuclPos[copy_index]  = DSNuclPos[copy_index];
        TmpScores[copy_index]   = DSScores[copy_index];
        TmpGaps[copy_index]     = DSGaps[copy_index];
      }
      free(DSTrans);    DSTrans    = TmpTrans;
      free(DSModelPos); DSModelPos = TmpModelPos;
      free(DSNuclPos);  DSNuclPos  = TmpNuclPos;
      free(DSScores);   DSScores   = TmpScores;
      free(DSGaps);     DSGaps     = TmpGaps;
    }

  }



  // We're really just interested in the sum scores on each side,
  // so let's switch over to that.
  // Similarly, we want the sum count of gaps up to each position.
  //
  int i;
  for (i=1; i<us_trans_len; i++) {
    USScores[i] += USScores[i-1];
    USGaps[i] += USGaps[i-1];
  }

  for (i=ds_trans_len-2; i>=0; i--)
    DSScores[i] += DSScores[i+1];

  // SO, this may seem weird, but because we count to the first position
  // in the downstream hit from the left in 'SpliceOverlappingDomains'
  // we're going to be counting gaps from the left.
  // 
  for (i=1; i<ds_trans_len; i++)
    DSGaps[i] += DSGaps[i-1];


  // In order to determine the score difference created by the
  // splice, we'll find the sum score at each position right
  // before the extension, and those become the baseline to
  // remove for determining the contribution of our cut.
  //
  float baseline_score = DSScores[Overlap->downstream_ext_len];
  baseline_score      += USScores[us_pre_ext_end_pos];



  // NOTE: As I note down below, too, there is currently a known
  //   bug where infinite scores can occur.  Eventually, this should
  //   be investigated, but for now I'm just patching it out.
  if (isinf(baseline_score)) {
    free(USTrans);
    free(USModelPos);
    free(USNuclPos);
    free(USScores);
    free(USGaps);
    free(DSTrans);
    free(DSModelPos);
    free(DSNuclPos);
    free(DSScores);
    free(DSGaps);
    if (DEBUGGING) DEBUG_OUT("'FindOptimalSpliceSite' Complete (BUT DUE TO INFINITE SCORE ERROR)",-1);
    return EDGE_FAIL_SCORE;
  }



  // What position in the model are we splitting on?
  int   optimal_us_pos     = 0;
  int   optimal_ds_pos     = 0;
  int   optimal_model_pos  = 0;
  int   optimal_splice_opt = 0;
  float optimal_score      = EDGE_FAIL_SCORE;


  // Arrays we'll use for splice site evaluation
  ESL_DSQ * UN = malloc(6*sizeof(ESL_DSQ));
  ESL_DSQ * DN = malloc(6*sizeof(ESL_DSQ));


  int us_start = 0;
  int ds_start = 0;
  for (model_pos = Overlap->amino_start; model_pos <= Overlap->amino_end; model_pos++) {


    while (USModelPos[us_start] < model_pos) us_start++;
    while (DSModelPos[ds_start] < model_pos) ds_start++;


    int us_pos = us_start;
    int ds_pos;
    while (us_pos < us_trans_len && USModelPos[us_pos] == model_pos) {


      ds_pos = ds_start;


      while (ds_pos < ds_trans_len && DSModelPos[ds_pos] == model_pos) {



        // Pull in the "contested" nucleotides for this position
        GetContestedUpstreamNucls(Overlap->UpstreamNucls,USNuclPos[us_pos],UN);
        GetContestedDownstreamNucls(Overlap->DownstreamNucls,DSNuclPos[ds_pos],DN);

  
        float splice_score;
        int splice_option = SelectSpliceOpt(Overlap,UN,DN,model_pos,gm,gcode,&splice_score);


        float sum_score = USScores[us_pos] + DSScores[ds_pos] + splice_score;


        // NOTE: I'm not sure what causes scores of 'inf' but
        //   this does seem to happen in rare cases, so for now
        //   I'm just going to try to force it not to create
        //   downstream chaos.
        //
        if (!isinf(sum_score) && sum_score > optimal_score) {
          optimal_score      = sum_score;
          optimal_us_pos     = us_pos;
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


  free(USTrans);
  free(USModelPos);
  free(USNuclPos);
  free(USScores);
  
  free(DSTrans);
  free(DSModelPos);
  free(DSNuclPos);
  free(DSScores);

  free(UN);
  free(DN);



  if (optimal_score == EDGE_FAIL_SCORE) {
    free(USGaps);
    free(DSGaps);
    if (DEBUGGING) DEBUG_OUT("'FindOptimalSpliceSite' Complete (but sad)",-1);
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


  if (DEBUGGING) DEBUG_OUT("'FindOptimalSpliceSite' Complete",-1);


  return optimal_score-baseline_score;

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: SpliceOverlappingDomains
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
void SpliceOverlappingDomains
(
  DOMAIN_OVERLAP * Overlap,
  TARGET_SEQ     * TargetNuclSeq,
  P7_PROFILE     * gm, 
  ESL_GENCODE    * gcode
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'SpliceOverlappingDomains'",1);


  // These splice indices are the terminal nucleotides *within* the
  // uncontested coding region, relative to the nucleotide sequences
  // in Overlap
  int   upstream_splice_index;
  int downstream_splice_index;
  int split_amino_model_index;
  int codon_split_option;
  float splice_score = FindOptimalSpliceSite(Overlap,gm,gcode,&upstream_splice_index,&downstream_splice_index,&split_amino_model_index,&codon_split_option);



  // This would be really bizarre (and worth checking to see if it
  // ever actually happens...)
  if (codon_split_option == -1) {
    Overlap->score         = EDGE_FAIL_SCORE;
    Overlap->score_density = EDGE_FAIL_SCORE;
    if (DEBUGGING) DEBUG_OUT("'SpliceOverlappingDomains' Complete (BUT WITH TERRIBLE OPTIONS?!)",-1);
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


  Overlap->score = splice_score;
  Overlap->score_density = splice_score / (1 + Overlap->amino_end - Overlap->amino_start);


  if (DEBUGGING) DEBUG_OUT("'SpliceOverlappingDomains' Complete",-1);

}






/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetNuclRangesFromAminoCoords
 *
 *  Desc. :
 *
 *  Inputs:  1. Edge :
 *
 *  Output:
 *
 */
void GetNuclRangesFromAminoCoords 
(DOMAIN_OVERLAP * Edge)
{

  P7_ALIDISPLAY *   UpDisp = Edge->UpstreamDisplay;
  P7_ALIDISPLAY * DownDisp = Edge->DownstreamDisplay;


  int strand = 1;
  if (UpDisp->sqfrom > UpDisp->sqto)
    strand = -1;


  //
  // UPSTREAM
  //

  Edge->upstream_disp_start = UpDisp->N - 1;
  int disp_amino = UpDisp->hmmto;

  Edge->upstream_nucl_end   = UpDisp->sqto + (3 * strand * (Edge->amino_end - disp_amino));
  Edge->upstream_nucl_start = UpDisp->sqto + strand;

  while (Edge->upstream_disp_start >= 0 && disp_amino >= Edge->amino_start) {

    Edge->upstream_nucl_start -= 3 * strand;
    disp_amino--;

    // If we were presumptuous, undo the previous work
    if (UpDisp->model[Edge->upstream_disp_start] == '.') // Insertion relative to pHMM
      disp_amino++;
    if (UpDisp->aseq[Edge->upstream_disp_start] == '-') // Insertion relative to genome
      Edge->upstream_nucl_start += 3 * strand;

    Edge->upstream_disp_start -= 1;

  }
  disp_amino++;
  Edge->upstream_ext_len     = disp_amino - Edge->amino_start;
  Edge->upstream_nucl_start -= 3 * strand * Edge->upstream_ext_len;
  
  // We'll have overstepped by one
  Edge->upstream_disp_start += 1;



  //
  // DOWNSTREAM
  //

  Edge->downstream_disp_end = 0;
  disp_amino = DownDisp->hmmfrom;

  Edge->downstream_nucl_start = DownDisp->sqfrom - (3 * strand * (disp_amino - Edge->amino_start));
  Edge->downstream_nucl_end   = DownDisp->sqfrom - strand;

  while (Edge->downstream_disp_end < DownDisp->N && disp_amino <= Edge->amino_end) {

    Edge->downstream_nucl_end += 3 * strand;
    disp_amino++;

    // Presumptuous? I hardly knumptuous!
    if (DownDisp->model[Edge->downstream_disp_end] == '.') // Insertion relative to pHMM
      disp_amino--;
    if (DownDisp->aseq[Edge->downstream_disp_end] == '-') // Insertion relative to genome
      Edge->downstream_nucl_end -= 3 * strand;

    Edge->downstream_disp_end += 1;

  }
  disp_amino--;
  Edge->downstream_ext_len   = Edge->amino_end - disp_amino;
  Edge->downstream_nucl_end += 3 * strand * Edge->downstream_ext_len;

  // We'll have overstepped by one
  Edge->upstream_disp_start -= 1;


}





/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: SketchSpliceEdge
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
void SketchSpliceEdge
(
  DOMAIN_OVERLAP * Edge,
  TARGET_SEQ     * TargetNuclSeq,
  P7_PROFILE     * gm,
  ESL_GENCODE    * gcode
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'SketchSpliceEdge'",1);


  P7_ALIDISPLAY *   UpDisp = Edge->UpstreamDisplay;
  P7_ALIDISPLAY * DownDisp = Edge->DownstreamDisplay;


  Edge->amino_start = DownDisp->hmmfrom;
  Edge->amino_end   =   UpDisp->hmmto;


  // Do we need to extend beyond the bounds of these
  // hits to have the required number of overlapping
  // amino acids?
  int num_ext_aminos = MIN_AMINO_OVERLAP - (1 + Edge->amino_end - Edge->amino_start);
  if (num_ext_aminos > 0) {

    // Extending the overlap region means going in
    // both directions, so we'll extend each side by
    // half of what's required
    num_ext_aminos = (num_ext_aminos+1)/2;

    Edge->amino_start -= num_ext_aminos;
    Edge->amino_end   += num_ext_aminos;

  }


  // Now we can do the work of finding the (indel-aware)
  // upstream_start and downstream_end coordinates.
  GetNuclRangesFromAminoCoords(Edge);


  // Grab them nucleotides!
  // Note that in order to ensure we can consider *every*
  // position in the overlap, we need to grab 2 extra
  // nucleotides to serve as splice signal candidates.
  if (Edge->upstream_nucl_start < Edge->upstream_nucl_end) {
    Edge->UpstreamNucls   = GrabNuclRange(TargetNuclSeq,Edge->upstream_nucl_start,Edge->upstream_nucl_end+2);
    Edge->DownstreamNucls = GrabNuclRange(TargetNuclSeq,Edge->downstream_nucl_start-2,Edge->downstream_nucl_end);
  } else {
    Edge->UpstreamNucls   = GrabNuclRange(TargetNuclSeq,Edge->upstream_nucl_start,Edge->upstream_nucl_end-2);
    Edge->DownstreamNucls = GrabNuclRange(TargetNuclSeq,Edge->downstream_nucl_start+2,Edge->downstream_nucl_end);
  }


  // Finish off by adding this friendly little pointer
  Edge->ntalpha = TargetNuclSeq->abc;


  // Big ol' DEBUGGING dump
  if (DEBUGGING && 1) {
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



  SpliceOverlappingDomains(Edge,TargetNuclSeq,gm,gcode);


  if (DEBUGGING) DEBUG_OUT("'SketchSpliceEdge' Complete",-1);

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

   
  if (DEBUGGING) DEBUG_OUT("Starting 'HitsAreSpliceCompatible'",1);

  
  // Start by checking if we either have amino acid
  // overlap, or are close enough to consider extending
  int amino_start_1 = Upstream->hmmfrom;
  int amino_end_1   = Upstream->hmmto;

  int amino_start_2 = Downstream->hmmfrom;
  int amino_end_2   = Downstream->hmmto;



  // If the upstream ain't upstream, then obviously we can't treat
  // these as splice-compatible!
  //
  // NOTE: We'll allow these to overlap fully to accommodate cases
  //       where there's an optimal splice junction in the middle
  //       that fixes overextension
  //
  if (amino_start_1 > amino_start_2 || amino_end_1 > amino_end_2) {
    if (DEBUGGING) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
    return 0;
  }


  // Do we have overlap OR sufficient proximity to consider
  // extending?
  if (amino_end_1 + MAX_AMINO_EXT < amino_start_2) {
    if (DEBUGGING) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
    return 0;
  }


  // Fantastic!  The amino acid coordinates support splice
  // compatibility!  Now it's just time to confirm that
  // the nucleotides also look good.

  int  nucl_start_1 = Upstream->sqfrom;
  int  nucl_end_1   = Upstream->sqto;
  
  int revcomp1 = 0;
  if (nucl_start_1 > nucl_end_1)
    revcomp1 = 1;


  int nucl_start_2 = Downstream->sqfrom;
  int nucl_end_2   = Downstream->sqto;
 
  int revcomp2 = 0;
  if (nucl_start_2 > nucl_end_2)
    revcomp2 = 1;


  if (revcomp1 != revcomp2) {
    if (DEBUGGING) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
    return 0;
  }


  // We want to make sure that these aren't unrealistically
  // close together on the genome...
  if (revcomp1) {

    if (nucl_start_2 + (3 * MAX_AMINO_EXT) >= nucl_end_1) {
      if (DEBUGGING) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
      return 0;
    }

  } else {

    if (nucl_start_2 - (3 * MAX_AMINO_EXT) <= nucl_end_1) {
      if (DEBUGGING) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);
      return 0;
    }

  }


  if (DEBUGGING) DEBUG_OUT("'HitsAreSpliceCompatible' Complete",-1);


  // Looks like we've got a viable upstream / downstream pair!
  return 1;

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
 *  Function: OutsideSearchArea
 *
 */
int OutsideSearchArea
(P7_HIT * Hit, TARGET_SEQ * TargetNuclSeq)
{
  
  if (strcmp(TargetNuclSeq->SeqName,Hit->name))
    return 1;

  int dom_id = Hit->best_domain;

  int64_t sqfrom = (&Hit->dcl[dom_id])->ad->sqfrom;
  int64_t sqto = (&Hit->dcl[dom_id])->ad->sqto;

  if ( TargetNuclSeq->revcomp && sqfrom < sqto) return 1;
  if (!TargetNuclSeq->revcomp && sqfrom > sqto) return 1;

  if (sqfrom < TargetNuclSeq->start) return 1;
  if (sqfrom > TargetNuclSeq->end  ) return 1;

  if (sqto < TargetNuclSeq->start) return 1;
  if (sqto > TargetNuclSeq->end  ) return 1;

  // Nice work coloring inside the lines, hit!
  return 0;

}











/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GatherViableSpliceEdges
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
DOMAIN_OVERLAP ** GatherViableSpliceEdges
(
  P7_TOPHITS  * TopHits,
  TARGET_SEQ  * TargetNuclSeq,
  P7_PROFILE  * gm,
  ESL_GENCODE * gcode,
  int * num_splice_edges
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'GatherViableSpliceEdges'",1);


  int num_hits = (int)(TopHits->N);

  int edge_capacity = 2 * num_hits;
  DOMAIN_OVERLAP ** SpliceEdges = (DOMAIN_OVERLAP **)malloc(edge_capacity * sizeof(DOMAIN_OVERLAP *));

  int num_edges = 0;
  int   upstream_hit_id;
  int   upstream_dom_id;
  int downstream_dom_id;
  int downstream_hit_id;
  for (upstream_hit_id = 0; upstream_hit_id < num_hits; upstream_hit_id++) {


    P7_HIT * UpstreamHit  = TopHits->hit[upstream_hit_id];
    int num_upstream_doms = UpstreamHit->ndom;


    // Is this hit inside the valid search area?
    if (OutsideSearchArea(UpstreamHit,TargetNuclSeq))
      continue;


    // For each hit, gather all of the indices of other hits that
    // could potentially be downstream exons.
    for (downstream_hit_id=0; downstream_hit_id < (int)(TopHits->N); downstream_hit_id++) {


      P7_HIT * DownstreamHit  = TopHits->hit[downstream_hit_id];
      int num_downstream_doms = DownstreamHit->ndom;


      // Is the downstream hit chill?
      if (OutsideSearchArea(DownstreamHit,TargetNuclSeq))
        continue;


      for (upstream_dom_id = 0; upstream_dom_id < num_upstream_doms; upstream_dom_id++) {


        P7_ALIDISPLAY * UpstreamDisplay = (&UpstreamHit->dcl[upstream_dom_id])->ad;


        if (ExcessiveGapContent(UpstreamDisplay))       
          continue;


        for (downstream_dom_id = 0; downstream_dom_id < num_downstream_doms; downstream_dom_id++) {

          // NO SELF-SPLICING, YOU ABSOLUTE MANIAC!
          if (upstream_hit_id == downstream_hit_id && upstream_dom_id == downstream_dom_id)
            continue;


          P7_ALIDISPLAY * DownstreamDisplay = (&DownstreamHit->dcl[downstream_dom_id])->ad;


          if (ExcessiveGapContent(DownstreamDisplay))
            continue;


          if (HitsAreSpliceCompatible(UpstreamDisplay,DownstreamDisplay)) {

            // MUST WE RESIZE?!
            if (num_edges == edge_capacity) {

              edge_capacity *= 2;

              DOMAIN_OVERLAP ** MoreSpliceEdges = (DOMAIN_OVERLAP **)malloc(edge_capacity * sizeof(DOMAIN_OVERLAP *));
              int edge_id;
              for (edge_id=0; edge_id<num_edges; edge_id++)
                MoreSpliceEdges[edge_id] = SpliceEdges[edge_id];

              free(SpliceEdges);
              SpliceEdges = MoreSpliceEdges;

            }


            // Record that splice compatibility!
            SpliceEdges[num_edges] = (DOMAIN_OVERLAP *)malloc(sizeof(DOMAIN_OVERLAP));
            DOMAIN_OVERLAP * Edge  = SpliceEdges[num_edges];

            Edge->upstream_hit_id   = upstream_hit_id;
            Edge->upstream_dom_id   = upstream_dom_id;
            Edge->downstream_hit_id = downstream_hit_id;
            Edge->downstream_dom_id = downstream_dom_id;

            Edge->UpstreamTopHits   = TopHits;
            Edge->UpstreamDisplay   = UpstreamDisplay;
            Edge->DownstreamTopHits = TopHits;
            Edge->DownstreamDisplay = DownstreamDisplay;

            num_edges++;

          }

        }

      }

    }

  }


  //
  //  Now that we have our splice edges, we can more fully
  //  sketch out how they connect!
  //
  //
  //  This *could* be part of the above loop, but what's the
  //  rush?
  //


  // We'll run through all of our paired domains and actually
  // splice 'em up (or at least try our best to)!
  int splice_edge_id;
  for (splice_edge_id = 0; splice_edge_id < num_edges; splice_edge_id++) {

    SketchSpliceEdge(SpliceEdges[splice_edge_id],TargetNuclSeq,gm,gcode);

    // If we failed to find a reasonable splice site, then we'll
    // just rip this edge outta consideration.
    if (SpliceEdges[splice_edge_id]->score == EDGE_FAIL_SCORE) {
      free(SpliceEdges[splice_edge_id]);
      SpliceEdges[splice_edge_id] = NULL;
    }

  }



  if (DEBUGGING) DEBUG_OUT("'GatherViableSpliceEdges' Complete",-1);


  *num_splice_edges = num_edges;
  return SpliceEdges;

}











/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: InitSpliceNode
 *
 *  Desc. :
 *
 *  Inputs:  1.      Graph :
 *           2.    node_id :
 *           3.     hit_id :
 *           4.     dom_id :
 *           5. was_missed :
 *
 *  Output:
 *
 */
SPLICE_NODE * InitSpliceNode
(
  SPLICE_GRAPH * Graph,
  int node_id, 
  int hit_id, 
  int dom_id,
  int was_missed
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'InitSpliceNode'",1);
  

  SPLICE_NODE * NewNode = (SPLICE_NODE *)malloc(sizeof(SPLICE_NODE));

  NewNode->node_id = node_id;
  NewNode->hit_id  = hit_id;
  NewNode->dom_id  = dom_id;


  NewNode->out_edge_cap    = 10;
  NewNode->num_out_edges   =  0;
  NewNode->OutEdges        = (DOMAIN_OVERLAP **)malloc(NewNode->out_edge_cap*sizeof(DOMAIN_OVERLAP *));
  NewNode->DownstreamNodes = (SPLICE_NODE    **)malloc(NewNode->out_edge_cap*sizeof(SPLICE_NODE    *));

  NewNode->in_edge_cap    = 10;
  NewNode->num_in_edges   =  0;
  NewNode->best_in_edge   = -1;
  NewNode->InEdges        = (DOMAIN_OVERLAP **)malloc(NewNode->in_edge_cap*sizeof(DOMAIN_OVERLAP *));
  NewNode->UpstreamNodes  = (SPLICE_NODE    **)malloc(NewNode->in_edge_cap*sizeof(SPLICE_NODE    *));


  NewNode->was_missed = was_missed;


  if (was_missed) {
  
    NewNode->is_n_terminal = 0;
    if (Graph->MissedHits->hit[hit_id]->dcl->ad->hmmfrom == 1)
      NewNode->is_n_terminal = 1;

    NewNode->is_c_terminal = 0;
    if (Graph->MissedHits->hit[hit_id]->dcl->ad->hmmto == Graph->Model->M)
      NewNode->is_c_terminal = 1;

    NewNode->hit_score = Graph->MissedHits->hit[hit_id]->dcl->bitscore;

  } else {
  
    NewNode->is_n_terminal = 0;
    if ((&(Graph->TopHits->hit[hit_id]->dcl[dom_id]))->ad->hmmfrom == 1)
      NewNode->is_n_terminal = 1;

    NewNode->is_c_terminal = 0;
    if ((&(Graph->TopHits->hit[hit_id]->dcl[dom_id]))->ad->hmmto == Graph->Model->M)
      NewNode->is_c_terminal = 1;

    NewNode->hit_score = Graph->TopHits->hit[hit_id]->dcl[dom_id].bitscore;

  }


  // Initialize this score to a recognizable and impossibly low value
  NewNode->cumulative_score = EDGE_FAIL_SCORE;


  if (DEBUGGING) DEBUG_OUT("'InitSpliceNode' Complete",-1);


  return NewNode;

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
(DOMAIN_OVERLAP * Edge, SPLICE_GRAPH * Graph)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'ConnectNodesByEdge'",1);


  int upstream_node_id;
  if (Edge->UpstreamTopHits == Graph->TopHits) {
    upstream_node_id = Graph->TH_HitDomToNodeID[Edge->upstream_hit_id][Edge->upstream_dom_id];
  } else {
    upstream_node_id = Graph->MH_HitToNodeID[Edge->upstream_hit_id];
  }
  SPLICE_NODE * UpstreamNode = Graph->Nodes[upstream_node_id];


  int downstream_node_id;
  if (Edge->DownstreamTopHits == Graph->TopHits) {
    downstream_node_id = Graph->TH_HitDomToNodeID[Edge->downstream_hit_id][Edge->downstream_dom_id];
  } else {
    downstream_node_id = Graph->MH_HitToNodeID[Edge->downstream_hit_id];
  }
  SPLICE_NODE * DownstreamNode = Graph->Nodes[downstream_node_id];


  UpstreamNode->OutEdges[UpstreamNode->num_out_edges] = Edge;
  UpstreamNode->DownstreamNodes[UpstreamNode->num_out_edges] = DownstreamNode;
  UpstreamNode->num_out_edges += 1;


  DownstreamNode->InEdges[DownstreamNode->num_in_edges] = Edge;
  DownstreamNode->UpstreamNodes[DownstreamNode->num_in_edges]  = UpstreamNode;
  DownstreamNode->num_in_edges += 1;



  // Resize?
  int edge_id;
  if (UpstreamNode->num_out_edges == UpstreamNode->out_edge_cap) {

    UpstreamNode->out_edge_cap *= 2;

    DOMAIN_OVERLAP ** NewOutEdges = (DOMAIN_OVERLAP **)malloc(UpstreamNode->out_edge_cap*sizeof(DOMAIN_OVERLAP *)); 
    SPLICE_NODE    ** NewDSNodes  = (SPLICE_NODE    **)malloc(UpstreamNode->out_edge_cap*sizeof(SPLICE_NODE    *));

    for (edge_id=0; edge_id<UpstreamNode->num_out_edges; edge_id++) {
      NewOutEdges[edge_id] = UpstreamNode->OutEdges[edge_id];
      NewDSNodes[edge_id]  = UpstreamNode->DownstreamNodes[edge_id];
    }

    free(UpstreamNode->OutEdges);
    free(UpstreamNode->DownstreamNodes);

    UpstreamNode->OutEdges = NewOutEdges;
    UpstreamNode->DownstreamNodes = NewDSNodes;

  }
  // Resize?
  if (DownstreamNode->num_in_edges == DownstreamNode->in_edge_cap) {

    DownstreamNode->in_edge_cap *= 2;

    DOMAIN_OVERLAP ** NewInEdges = (DOMAIN_OVERLAP **)malloc(DownstreamNode->in_edge_cap*sizeof(DOMAIN_OVERLAP *)); 
    SPLICE_NODE    ** NewUSNodes = (SPLICE_NODE    **)malloc(DownstreamNode->in_edge_cap*sizeof(SPLICE_NODE    *));

    for (edge_id=0; edge_id<DownstreamNode->num_in_edges; edge_id++) {
      NewInEdges[edge_id] = DownstreamNode->InEdges[edge_id];
      NewUSNodes[edge_id] = DownstreamNode->UpstreamNodes[edge_id];
    }

    free(DownstreamNode->InEdges);
    free(DownstreamNode->UpstreamNodes);

    DownstreamNode->InEdges = NewInEdges;
    DownstreamNode->UpstreamNodes = NewUSNodes;

  }



  // Node connected... by an edge!
  Graph->num_edges += 1;


  if (DEBUGGING) DEBUG_OUT("'ConnectNodesByEdge' Complete",-1);

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

  if (DEBUGGING) DEBUG_OUT("Starting 'GatherCTermNodes'",1);

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

  if (DEBUGGING) DEBUG_OUT("'GatherCTermNodes' Complete",-1);

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
  DOMAIN_OVERLAP * DO_2_3 = Node3->InEdges[n3_in_edge_index];


  // If Node2 doesn't have a best input edge, there's
  // no risk of overrunning it
  if (Node2->best_in_edge == -1)
    return 0;


  // Who's the best upstream node for Node2?
  // Note that this has already been determined due to the
  // recursive nature of 'PullUpCumulativeScore'
  int n2_in_edge_index    = Node2->best_in_edge;
  DOMAIN_OVERLAP * DO_1_2 = Node2->InEdges[n2_in_edge_index];

  
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
    Node3->num_in_edges -= 1;
    if (n3_in_edge_index < Node3->num_in_edges) {
      Node3->InEdges[n3_in_edge_index] = Node3->InEdges[Node3->num_in_edges];
      Node3->UpstreamNodes[n3_in_edge_index] = Node3->UpstreamNodes[Node3->num_in_edges];
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

  if (DEBUGGING) DEBUG_OUT("Starting 'PullUpCumulativeScore'",1);

  if (Node->cumulative_score != EDGE_FAIL_SCORE) {
    if (DEBUGGING) DEBUG_OUT("'PullUpCumulativeScore' Complete",-1);
    return Node->cumulative_score;
  }


  // Since 'EdgeWouldEraseNode' can remove edges we'll need to
  // use a 'while' loop instead of a 'for' loop
  int in_edge_index = 0;
  while (in_edge_index < Node->num_in_edges) {

    float in_edge_score = Node->InEdges[in_edge_index]->score;
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


  if (DEBUGGING) DEBUG_OUT("'PullUpCumulativeScore' Complete",-1);


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

  if (DEBUGGING) DEBUG_OUT("Starting 'EvaluatePaths'",1);


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


  if (DEBUGGING) DEBUG_OUT("'EvaluatePaths' Complete",-1);

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: FillOutGraphStructure
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
void FillOutGraphStructure
(
  SPLICE_GRAPH    * Graph, 
  TARGET_SEQ      * TargetNuclSeq,
  DOMAIN_OVERLAP ** SpliceEdges, 
  int num_splice_edges
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'FillOutGraphStructure'",1);

   int i;
  // We'll want to be a little careful, just because this is our flagship
  // datastructure...
  int node_id;
  if (Graph->Nodes != NULL) {
    for (node_id=1; node_id<=Graph->num_nodes; node_id++)
      if (Graph->Nodes[node_id] != NULL)
        free(Graph->Nodes[node_id]);
    free(Graph->Nodes);
  }


  // We'll want this lookup table to be able to go from
  // DOMAIN_OVERLAP content
  int num_hits = (int)(Graph->TopHits->N);
  int max_doms = 0;
  int sum_doms = 0;
  int hit_id;
  int dom_id;
  for (hit_id=0; hit_id<num_hits; hit_id++) {
    int hit_doms = Graph->TopHits->hit[hit_id]->ndom;
    if (hit_doms > max_doms)
      max_doms = hit_doms;
    sum_doms += hit_doms;
  }

  // Allocate space for the maximum number of domains
  Graph->Nodes = (SPLICE_NODE **)malloc((sum_doms+1)*sizeof(SPLICE_NODE *));


  Graph->TH_HitDomToNodeID = malloc(num_hits * sizeof(int *));

  
  node_id = 0;
  
  for (hit_id=0; hit_id<num_hits; hit_id++) {

    Graph->TH_HitDomToNodeID[hit_id] = malloc(max_doms * sizeof(int));

    // 1. I don't think we ever actually use this matrix...
    // 2. We need to confirm that this hit isn't out of bounds.
    // 3. Even if it is out of bounds, we'll just have zeros
    int out_of_bounds = OutsideSearchArea(Graph->TopHits->hit[hit_id],TargetNuclSeq);
    
    for (dom_id=0; dom_id<max_doms; dom_id++) {

      Graph->TH_HitDomToNodeID[hit_id][dom_id] = 0;

      if (!out_of_bounds && dom_id < Graph->TopHits->hit[hit_id]->ndom) {

        Graph->TH_HitDomToNodeID[hit_id][dom_id] = ++node_id;
        
        Graph->Nodes[node_id] = InitSpliceNode(Graph,node_id,hit_id,dom_id,0);
      }
    }
  }
  Graph->num_nodes = node_id;


  // Are your hits to the reverse complement of the query nucleotide seq?
  //
  // Even though this is redundant, there are some instances where we just
  // pass along the graph without the TargetNuclSeq (although, really, TNS
  // should be part of Graph, rather than a separate entity...)
  //
  Graph->revcomp = TargetNuclSeq->revcomp;


  int edge_id;
  for (edge_id=0; edge_id<num_splice_edges; edge_id++) {
    if (SpliceEdges[edge_id] != NULL)
      ConnectNodesByEdge(SpliceEdges[edge_id],Graph);
  }


  EvaluatePaths(Graph);

  
  if (DEBUGGING) DEBUG_OUT("'FillOutGraphStructure' Complete",-1);

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
  
  if (DEBUGGING) DEBUG_OUT("Starting 'FindBestFullPath'",1);


  // NOTE: Because CTermNodeIDs is sorted by cumulative score,
  //       we can break the first time we find a full path.
  SPLICE_NODE * Walker;
  int c_term_index;
  for (c_term_index=0; c_term_index<Graph->num_c_term; c_term_index++) {

    Walker = Graph->Nodes[Graph->CTermNodeIDs[c_term_index]];

    Graph->best_full_path_length = 1;
    while (Walker->best_in_edge != -1) {
      Walker = Walker->UpstreamNodes[Walker->best_in_edge];
      Graph->best_full_path_length += 1;
    }

    if (Walker->is_n_terminal) {
      Graph->has_full_path        = 1;
      Graph->best_full_path_start = Walker->node_id;
      Graph->best_full_path_end   = Graph->CTermNodeIDs[c_term_index];
      Graph->best_full_path_score = Graph->Nodes[Graph->best_full_path_end]->cumulative_score;
      if (DEBUGGING) DEBUG_OUT("'FindBestFullPath' Complete",-1);
      return;
    } else {
      // RESET!
      Graph->best_full_path_length = 0;
    }

  }

  if (DEBUGGING) DEBUG_OUT("'FindBestFullPath' Complete",-1);

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: BuildSpliceGraph
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
SPLICE_GRAPH * BuildSpliceGraph
(
  P7_TOPHITS  * TopHits, 
  TARGET_SEQ  * TargetNuclSeq,
  P7_PROFILE  * gm,
  P7_OPROFILE * om,
  ESL_GENCODE * gcode
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'BuildSpliceGraph'",1);

  // We'll just make an unordered list of our splice edges for now
  int num_splice_edges = 0;
  DOMAIN_OVERLAP ** SpliceEdges = GatherViableSpliceEdges(TopHits,TargetNuclSeq,gm,gcode,&num_splice_edges);


  SPLICE_GRAPH * Graph = (SPLICE_GRAPH *)malloc(sizeof(SPLICE_GRAPH));


  Graph->TopHits    = TopHits;
  Graph->MissedHits = NULL;
  
  Graph->Model  = gm;
  Graph->OModel = om;

  Graph->TH_HitDomToNodeID = NULL;
  Graph->MH_HitToNodeID    = NULL;

  Graph->Nodes        = NULL;
  Graph->CTermNodeIDs = NULL;


  // Initialize our basic metadata
  Graph->num_nodes  = 0;
  Graph->num_edges  = 0;
  Graph->num_n_term = 0;
  Graph->num_c_term = 0;


  // Eventually, we'll need to know if we have a full
  // path through this graph.
  Graph->has_full_path        = 0;
  Graph->best_full_path_start = 0;
  Graph->best_full_path_end   = 0;
  Graph->best_full_path_score = 0.0;


  // Build that stinky graph!
  FillOutGraphStructure(Graph,TargetNuclSeq,SpliceEdges,num_splice_edges);
  FindBestFullPath(Graph);

  free(SpliceEdges);

  if (DEBUGGING) DEBUG_OUT("'BuildSpliceGraph' Complete",-1);

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

  if (DEBUGGING) DEBUG_OUT("Starting 'ExtractSubProfile'",1);


  int fullM = FullModel->M;
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



  if (DEBUGGING && hmm_start_pos == 1 && hmm_end_pos == FullModel->M) {
    fprintf(stderr,"\n  Full Model Copy Validation:  ");
    if (p7_profile_Compare(SubModel,FullModel,0.0) != eslOK) { fprintf(stderr,"Failed\n\n");   }
    else                                                     { fprintf(stderr,"Success!\n\n"); }
    fflush(stderr);
  }


  if (DEBUGGING) DEBUG_OUT("'ExtractSubProfile' Complete",-1);


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
  P7_DOMAIN * USDom     = &(Graph->TopHits->hit[UpstreamNode->hit_id]->dcl[UpstreamNode->dom_id]);
  int upstream_hmm_to   = USDom->ad->hmmto;
  int upstream_nucl_to  = USDom->ad->sqto;


  P7_DOMAIN * DSDom        = &(Graph->TopHits->hit[DownstreamNode->hit_id]->dcl[DownstreamNode->dom_id]);
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

  if (DEBUGGING) DEBUG_OUT("Starting 'IdentifyDisConnComponents'",1);


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
      DomPtr = &(Graph->TopHits->hit[Graph->Nodes[node_id]->hit_id]->dcl[Graph->Nodes[node_id]->dom_id]);


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
      DomPtr = &(Graph->TopHits->hit[Graph->Nodes[node_id]->hit_id]->dcl[Graph->Nodes[node_id]->dom_id]);


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



  if (DEBUGGING) DEBUG_OUT("'IdentifyDisConnComponents' Complete",-1);


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

  if (DEBUGGING) DEBUG_OUT("Starting 'IdentifyMidSearchRegions'",1);


  int   max_mid_searches = Graph->num_nodes; // Just for mem tracking
  int * MidSearchRegions = malloc((1 + 4 * max_mid_searches) * sizeof(int));
  int   num_mid_searches = 0;


  int us_node_id, ds_index, ds_node_id, ad_pos, contiguous_stars;
  for (us_node_id=1; us_node_id<=Graph->num_nodes; us_node_id++) {


    SPLICE_NODE * USNode = Graph->Nodes[us_node_id];


    for (ds_index=0; ds_index<USNode->num_out_edges; ds_index++) {

      
      DOMAIN_OVERLAP * DO   = USNode->OutEdges[ds_index];
      P7_ALIDISPLAY  * USAD = DO->UpstreamDisplay;
      P7_ALIDISPLAY  * DSAD = DO->DownstreamDisplay;


      // Before anything else, is there enough space between the hits
      // that we'd be well-served looking at it?
      int   upstream_nucl_end   = (int)(USAD->sqto);
      int downstream_nucl_start = (int)(DSAD->sqfrom);
      if (abs(downstream_nucl_start - upstream_nucl_end) < 50)
        continue;

      
      // What are the pHMM positions at which we've observed two '*'s
      // in the posterior probability output?
      int   upstream_strong_hmmto    = USAD->hmmto;
      int downstream_strong_hmmfrom  = DSAD->hmmfrom;


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
      if (downstream_strong_hmmfrom - upstream_strong_hmmto <= 6)
        continue;


      // Cool!  Why *not* give it a look?
      MidSearchRegions[4*num_mid_searches+1] =   upstream_strong_hmmto;
      MidSearchRegions[4*num_mid_searches+2] = downstream_strong_hmmfrom;
      MidSearchRegions[4*num_mid_searches+3] =   upstream_nucl_end;
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



  if (DEBUGGING) DEBUG_OUT("'IdentifyMidSearchRegions' Complete",-1);


  MidSearchRegions[0] = num_mid_searches;

  return MidSearchRegions;

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: IdentifyTermSearchRegions
 *
 */
int * IdentifyTermSearchRegions
(
  SPLICE_GRAPH * Graph,
  TARGET_SEQ   * TargetNuclSeq,
  int * NoInEdgeNodes,
  int   num_no_in_edge,
  int * NoOutEdgeNodes,
  int   num_no_out_edge
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'IdentifyTermSearchRegions'",1);


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

      DomPtr = &(Graph->TopHits->hit[Graph->Nodes[node_id]->hit_id]->dcl[Graph->Nodes[node_id]->dom_id]);
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
        TermSearchRegions[4*num_term_searches+3] = intMin(upstreamest_nucl_pos+10000,(int)(TargetNuclSeq->end));
      } else {
        TermSearchRegions[4*num_term_searches+3] = intMax(upstreamest_nucl_pos-10000,(int)(TargetNuclSeq->start));
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

      DomPtr = &(Graph->TopHits->hit[Graph->Nodes[node_id]->hit_id]->dcl[Graph->Nodes[node_id]->dom_id]);
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
        TermSearchRegions[4*num_term_searches+4] = intMax(downstreamest_nucl_pos-10000,(int)(TargetNuclSeq->start));
      } else {
        TermSearchRegions[4*num_term_searches+4] = intMin(downstreamest_nucl_pos+10000,(int)(TargetNuclSeq->end));
      }

      // If this search region is too small, skip it
      if (abs(TermSearchRegions[4*num_term_searches+4] - downstreamest_nucl_pos) > 50)
        num_term_searches++;

    }

  }



  if (DEBUGGING) DEBUG_OUT("'IdentifyTermSearchRegions' Complete",-1);


  TermSearchRegions[0] = num_term_searches;

  return TermSearchRegions;

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GetBoundedSearchRegions
 *
 *  Desc. :
 *
 *  Inputs:
 *
 *  Output:
 *
 */
int * GetBoundedSearchRegions
(SPLICE_GRAPH * Graph, TARGET_SEQ * TargetNuclSeq)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'GetBoundedSearchRegions'",1);


  // Who all doesn't have an incoming / outgoing edge?
  // (not including "genuine" terminal nodes, of course!)
  int * NoOutEdgeNodes = malloc(Graph->num_nodes * sizeof(int));
  int * NoInEdgeNodes  = malloc(Graph->num_nodes * sizeof(int));
  int num_no_out_edge  = 0;
  int num_no_in_edge   = 0;

  int node_id;
  for (node_id = 1; node_id <= Graph->num_nodes; node_id++) {

    SPLICE_NODE * Node = Graph->Nodes[node_id];

    if (Node->num_out_edges == 0 && Node->is_c_terminal == 0)
      NoOutEdgeNodes[num_no_out_edge++] = node_id;

    if (Node->num_in_edges == 0 && Node->is_n_terminal == 0)
      NoInEdgeNodes[num_no_in_edge++] = node_id;

  }


  int *  DCCSearchRegions = IdentifyDisConnComponents(Graph,NoInEdgeNodes,num_no_in_edge,NoOutEdgeNodes,num_no_out_edge);
  int *  MidSearchRegions = IdentifyMidSearchRegions(Graph);
  int * TermSearchRegions = IdentifyTermSearchRegions(Graph,TargetNuclSeq,NoInEdgeNodes,num_no_in_edge,NoOutEdgeNodes,num_no_out_edge);


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


  if (DEBUGGING) DEBUG_OUT("'GetBoundedSearchRegions' Complete",-1);


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
  float * SubHitScores,
  int     num_sub_hits,
  int     sub_hmm_len,
  int     hmm_start,
  int   * final_num_sub_hits
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'SelectFinalSubHits'",1);


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


    // Do we get any fresh coverage out of this hit?
    int added_coverage = 0;
    for (model_pos = AD->hmmfrom - hmm_start; model_pos <= AD->hmmto - hmm_start; model_pos++) {
      if (Coverage[model_pos] == 0) {
        Coverage[model_pos] = 1;
        added_coverage++;
      }
    }


    // No coverage?! Away with you, filth!
    if (added_coverage < MIN_AMINO_OVERLAP) {
      p7_alidisplay_Destroy(AD);
      SubHitADs[sub_hit_id] = NULL;
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

    if (SubHitADs[sub_hit_id] == NULL)
      continue;

    FinalSubHits[*final_num_sub_hits]           = p7_domain_Create_empty();
    FinalSubHits[*final_num_sub_hits]->ad       = SubHitADs[sub_hit_id];
    FinalSubHits[*final_num_sub_hits]->bitscore = SubHitScores[sub_hit_id];

    *final_num_sub_hits += 1;

  }



  if (DEBUGGING) DEBUG_OUT("'SelectFinalSubHits' Complete",-1);


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
 *  Function: FindSubHits
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
P7_DOMAIN ** FindSubHits
(
  SPLICE_GRAPH * Graph,
  TARGET_SEQ   * TargetNuclSeq,
  int          * SearchRegion,
  ESL_GENCODE  * gcode,
  int          * final_num_sub_hits
)
{
  if (DEBUGGING) DEBUG_OUT("Starting 'FindSubHits'",1);


  int hmm_start  = SearchRegion[0];
  int hmm_end    = SearchRegion[1];
  int nucl_start = SearchRegion[2];
  int nucl_end   = SearchRegion[3];


  // We'll pull in a couple extra HMM positions, since
  // we're interested in having overlaps anyways
  hmm_start = intMax(hmm_start-2,1);
  hmm_end   = intMin(hmm_end  +2,Graph->Model->M);

  int sub_hmm_len = 1 + hmm_end - hmm_start;


  // Generate a sub-model and an optimized sub-model
  // for the part of the pHMM we need to fill in
  P7_PROFILE  *  SubModel = ExtractSubProfile(Graph->Model,hmm_start,hmm_end);
  P7_OPROFILE * OSubModel = p7_oprofile_Create(SubModel->M,SubModel->abc);
  int submodel_create_err = p7_oprofile_Convert(SubModel,OSubModel);


  // Required for alidisplay generation
  OSubModel->name = malloc(5*sizeof(char));
  strcpy(OSubModel->name,"OSM\0");


  
  // DEBUGGING
  //
  // Note that this is primarily intended to work when there's
  // a particular (tester-known) sequence that we want to confirm
  // would score well under the model.
  //
  //TestSubModel(SubModel,OSubModel,"PPNPSLMSIFRK\0");



  // Grab the nucleotides we're searching our sub-model against
  ESL_DSQ * SubNucls = GrabNuclRange(TargetNuclSeq,nucl_start,nucl_end);
  int nucl_seq_len   = abs(nucl_end-nucl_start)+1;


  // Create a whole mess of objects that we'll need to get our
  // our sub-model alignments to fit into an alidisplay...
  ESL_SQ   * ORFAminoSeq   = NULL;
  P7_GMX   * ViterbiMatrix = p7_gmx_Create(SubModel->M,1024);
  P7_TRACE * Trace         = p7_trace_Create();
  float viterbi_score;


  // Finally, the array where we'll record all of our incredible triumphs!
  int num_sub_hits =  0;
  int max_sub_hits = 10;
  P7_ALIDISPLAY ** SubHitADs = malloc(max_sub_hits * sizeof(P7_ALIDISPLAY *));
  float * SubHitScores = malloc(max_sub_hits * sizeof(float));


  // We'll build up char arrays of each ORF (and its coding nucleotides)
  // that the get converted into the ESL_SQs ORFAminos and ORFNucls
  int orf_str_cap = 400;
  char * AminoStr = malloc(    orf_str_cap * sizeof(char));
  char * NuclStr  = malloc(3 * orf_str_cap * sizeof(char));


  // Loop over each full reading frame
  int max_orf_len = 500;
  int frame,i;
  for (frame=0; frame<3; frame++) {

    int orf_len = 0;
    int frame_start = 1 + frame;

    // As we walk through the reading frame,
    // we'll build up ORFs and search them
    // against our sub-model
    int frame_end = frame_start;
    while (frame_end+2 < nucl_seq_len) {


      // Translate and add to the current frame
      int next_amino_index = esl_gencode_GetTranslation(gcode,&(SubNucls[frame_end]));
      if (next_amino_index < 21) {

        AminoStr[   orf_len] = AMINO_CHARS[next_amino_index];
        NuclStr[  3*orf_len] = DNA_CHARS[SubNucls[frame_end  ]];
        NuclStr[1+3*orf_len] = DNA_CHARS[SubNucls[frame_end+1]];
        NuclStr[2+3*orf_len] = DNA_CHARS[SubNucls[frame_end+2]];

        orf_len++;

        // Resize?
        if (orf_len == orf_str_cap-1) {

          orf_str_cap *= 2;

          char * NewAminoStr = malloc(  orf_str_cap*sizeof(char));
          char * NewNuclStr  = malloc(3*orf_str_cap*sizeof(char));
          for (i=0; i<  orf_len; i++) NewAminoStr[i] = AminoStr[i];
          for (i=0; i<3*orf_len; i++) NewNuclStr[i]  = NuclStr[i];

          free(AminoStr);
          AminoStr = NewAminoStr;

          free(NuclStr);
          NuclStr = NewNuclStr;

        }

      }
      frame_end += 3;


      // Have we hit a stop codon? Are we at the end of the road?
      if (next_amino_index >= 21 || frame_end+2 >= nucl_seq_len || orf_len == max_orf_len) {


        int long_orf_catch = 0;
        if (orf_len == max_orf_len)
          long_orf_catch = 1;


        // Is this ORF long enough to be worth considering as an exon?
        if (orf_len >= 5) {


          // 0-out the ends of the strings
          AminoStr[ orf_len] = 0;
          NuclStr[3*orf_len] = 0;


          ORFAminoSeq = esl_sq_CreateFrom("orf",AminoStr,NULL,NULL,NULL);
          int dsq_err_code  = esl_sq_Digitize(SubModel->abc,ORFAminoSeq);
          if (dsq_err_code != eslOK) {
            fprintf(stderr,"\n  ERROR (FindSubHits): Failed while digitizing ORF amino sequence\n\n");
          }


          int vit_err_code  = p7_GViterbi(ORFAminoSeq->dsq,orf_len,SubModel,ViterbiMatrix,&viterbi_score);
          if (vit_err_code != eslOK) {
            fprintf(stderr,"\n  ERROR (FindSubHits): Failed while running Viterbi\n\n");
          }
          

          // Swell!
          int gtrace_err_code  = p7_GTrace(ORFAminoSeq->dsq,orf_len,SubModel,ViterbiMatrix,Trace);
          if (gtrace_err_code != eslOK) {
            fprintf(stderr,"\n  ERROR (FindSubHits): Failed while generating a generic P7_TRACE for the Viterbi matrix\n\n");
          }


          int trace_index_err_code  = p7_trace_Index(Trace);
          if (trace_index_err_code != eslOK) {
            fprintf(stderr,"\n  ERROR (FindSubHits): Failed while indexing a P7_TRACE\n\n");
          }


          // For each of the domains in our P7_TRACE, see if
          // they look like solid missed exons
          int trace_dom;
          for (trace_dom = 0; trace_dom < Trace->ndom; trace_dom++) {


            P7_ALIDISPLAY * AD = p7_alidisplay_Create(Trace,0,OSubModel,ORFAminoSeq,NULL);
            if (AD == NULL) {
              fprintf(stderr,"\n  ERROR (FindSubHits): Failed while generating P7_ALIDISPLAY\n\n");
            }


            // If this is too short, we'll skip it
            int ali_len = AD->hmmto - AD->hmmfrom + 1;
            if (ali_len < 6) {
              p7_alidisplay_Destroy(AD);
              continue;
            }


            // Get the actual model range
            AD->hmmfrom += hmm_start - 1;
            AD->hmmto   += hmm_start - 1;


            // Add the nucleotide sequence to the ALI_DISPLAY
            AD->ntseq = malloc((3 * AD->N + 1) * sizeof(char));
            
            int nucl_read_pos = 3 * (AD->sqfrom - 1);
            for (i=0; i<AD->N; i++) {
              if (AD->aseq[i] == '-') {
                AD->ntseq[3*i  ] = '-';
                AD->ntseq[3*i+1] = '-';
                AD->ntseq[3*i+2] = '-';
              } else {
                AD->ntseq[3*i  ] = NuclStr[nucl_read_pos++];
                AD->ntseq[3*i+1] = NuclStr[nucl_read_pos++];
                AD->ntseq[3*i+2] = NuclStr[nucl_read_pos++];
              }
            }
            AD->ntseq[3*AD->N] = 0;


            // Compute the nucleotide offsets (within the extracted sequence)
            // Obviously, some of these operations could be combined, but
            // for the sake of making this intelligible I'm going to eat
            // the extra operations.
            //
            int nucl_sq_start = (frame_start + 3 * (AD->sqfrom - 1)    )-1;
            int nucl_sq_end   = (frame_start + 3 * (AD->sqto   - 1) + 2)-1;

            // Convert to the full-sequence coordinates
            if (Graph->revcomp) {
              AD->sqfrom = nucl_start - nucl_sq_start;
              AD->sqto   = nucl_start - nucl_sq_end;
            } else {
              AD->sqfrom = nucl_start + nucl_sq_start;
              AD->sqto   = nucl_start + nucl_sq_end;
            }


            // Hopefully this isn't too cruel, but if your score wasn't
            // very good we'll actually pull our support...
            SubHitScores[num_sub_hits] = ComputeRoughAliScore(AD,Graph->Model);
            float score_cutoff         = 2.0 * (float)ali_len / 3.0;
            if (SubHitScores[num_sub_hits] < score_cutoff || ExcessiveGapContent(AD)) {
              p7_alidisplay_Destroy(AD);
              continue;            
            }


            // Welcome to the list, fella!
            SubHitADs[num_sub_hits] = AD;
            num_sub_hits++;


            // DEBUGGING
            //fprintf(stdout,"PASSED! %d..%d (%f): %s\n",AD->hmmfrom,AD->hmmto,viterbi_score,AD->aseq);
            //p7_trace_Dump(stdout, Trace, SubModel, ORFAminoSeq->dsq);



            // Resize?
            if (num_sub_hits == max_sub_hits) {
              
              max_sub_hits *= 2;

              P7_ALIDISPLAY ** NewSubHitADs = malloc(max_sub_hits*sizeof(P7_ALIDISPLAY *));
              float * NewSubHitScores = malloc(max_sub_hits*sizeof(float));
              
              int sub_hit_id;
              for (sub_hit_id=0; sub_hit_id<num_sub_hits; sub_hit_id++) {
                NewSubHitADs[sub_hit_id]    = SubHitADs[sub_hit_id];
                NewSubHitScores[sub_hit_id] = SubHitScores[sub_hit_id];
              }
              free(SubHitADs);    SubHitADs    = NewSubHitADs;
              free(SubHitScores); SubHitScores = NewSubHitScores;
              
            }
            // Resize over -- back to our regularly-scheduled programming!


          }


          esl_sq_Reuse(ORFAminoSeq);
          p7_trace_Reuse(Trace);
          p7_gmx_Reuse(ViterbiMatrix);


        }


        // Reset and advance!
        // (but be mindful if there was a long ORF!)
        if (long_orf_catch) {

          int reuse_len      = 100;
          int amino_read_pos = max_orf_len - reuse_len;
          
          for (i=0; i<reuse_len; i++) {
            AminoStr[i]    = AminoStr[amino_read_pos];
            NuclStr[3*i  ] = NuclStr[3*amino_read_pos  ];
            NuclStr[3*i+1] = NuclStr[3*amino_read_pos+1];
            NuclStr[3*i+2] = NuclStr[3*amino_read_pos+2];
            amino_read_pos++;
          }

          orf_len = reuse_len;

        } else {

          orf_len     = 0;

        }

        frame_start = frame_end;


      }

    }

  }


  // It takes a lot of fun to need this much cleanup ;)
  free(AminoStr);
  free(NuclStr);


  if (ORFAminoSeq) esl_sq_Destroy(ORFAminoSeq);
  p7_profile_Destroy(SubModel);
  p7_oprofile_Destroy(OSubModel);
  p7_gmx_Destroy(ViterbiMatrix);
  p7_trace_Destroy(Trace);


  // Quick check: Did we find *anything* worth considering?
  if (num_sub_hits == 0) {
    free(SubHitADs);
    free(SubHitScores);
    if (DEBUGGING) DEBUG_OUT("'FindSubHits' Complete (albeit, no hits)",-1);
    *final_num_sub_hits = 0;
    return NULL;
  }


  // Reduce down to just the hits we're really excited about
  P7_DOMAIN ** FinalSubHits = SelectFinalSubHits(SubHitADs,SubHitScores,num_sub_hits,sub_hmm_len,hmm_start,final_num_sub_hits);


  free(SubHitADs);
  free(SubHitScores);


  if (DEBUGGING) DEBUG_OUT("'FindSubHits' Complete",-1);


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
(SPLICE_GRAPH * Graph, DOMAIN_OVERLAP ** NewSpliceEdges, int num_new_edges)
{
  if (DEBUGGING) DEBUG_OUT("Starting 'IntegrateMissedHits'",1);


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

    Graph->Nodes[node_id] = InitSpliceNode(Graph,node_id,missed_hit_id,0,1);

  }


  int new_edge_id;
  for (new_edge_id=0; new_edge_id<num_new_edges; new_edge_id++) {
    if (NewSpliceEdges[new_edge_id] != NULL)
      ConnectNodesByEdge(NewSpliceEdges[new_edge_id],Graph);
  }


  EvaluatePaths(Graph);


  if (DEBUGGING) DEBUG_OUT("'IntegrateMissedHits' Complete",-1);

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: SeekMissingExons
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
P7_TOPHITS * SeekMissingExons
(
  SPLICE_GRAPH * Graph, 
  TARGET_SEQ   * TargetNuclSeq, 
  ESL_GENCODE  * gcode
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'SeekMissingExons'",1);


  // Grab the sub-regions of our conceptual DP zone that we want
  // to search for missing exons.
  //
  // This 'SearchRegionAggregate' is indexed in groups consisting
  // of: the starting position in the model, the ending position in
  // the model, the starting coordinate on the genome, and the ending
  // coordinate on the genome.
  //
  int * SearchRegionAggregate = GetBoundedSearchRegions(Graph,TargetNuclSeq);
  int   num_search_regions    = SearchRegionAggregate[0];

  if (SearchRegionAggregate == NULL) {
    if (DEBUGGING) DEBUG_OUT("'SeekMissingExons' Complete (none found)",-1);
    return NULL;
  }


  if (DEBUGGING) {
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
    P7_DOMAIN ** NewSubHits = FindSubHits(Graph,TargetNuclSeq,&SearchRegionAggregate[(4*search_region_id)+1],gcode,&num_new_sub_hits);


    if (num_new_sub_hits == 0)
      continue;


    // New sub-hit(s) alert!
    if (DEBUGGING) {
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
    MissingHits->hit[sub_hit_id]->name  = TargetNuclSeq->SeqName;
    MissingHits->hit[sub_hit_id]->ndom  = 1;
    MissingHits->hit[sub_hit_id]->dcl   = SubHits[sub_hit_id];
    MissingHits->hit[sub_hit_id]->score = SubHits[sub_hit_id]->bitscore;

  }
  free(SubHits);

 
  if (DEBUGGING) DEBUG_OUT("'SeekMissingExons' Complete",-1);


  return MissingHits;

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: AddMissingExonsToGraph
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
void AddMissingExonsToGraph
(
  SPLICE_GRAPH * Graph,
  TARGET_SEQ   * TargetNuclSeq,
  ESL_GENCODE  * gcode
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'AddMissingExonsToGraph'",1);

  

  // As I note in 'SeekMissingExons' (but is worth emphasizing),
  // this datastructure is basically a cheap way to pass around
  // the new hits that we find.
  //
  // MUCH OF THE STANDARD P7_HIT, P7_DOMAIN, and P7_ALIDISPLAY DATA
  // IS NOT CONTAINED IN THESE!  BE CAREFUL BEFORE USING HMMER INTERNAL
  // FUNCTIONS TO INTERROGATE THEM!
  //
  Graph->MissedHits = SeekMissingExons(Graph,TargetNuclSeq,gcode);

  
  if (Graph->MissedHits == NULL) {
    if (DEBUGGING) DEBUG_OUT("'AddMissingExonsToGraph' Complete (no missed hits added)",-1);
    return;
  }


  // This is basically the same code from 'GatherViableSpliceEdges,'
  // but reworked to integrate new exons into the splice graph.
  //
  // It probably wouldn't be hard to convert the two to a single function,
  // but the annoying thing is that there's the whole dereferencing thing
  // to get *real* P7_DOMAINs from *actual* P7_HITs (but not my goofy ones)
  //
  int new_splice_edge_cap = 4 * Graph->MissedHits->N;
  DOMAIN_OVERLAP ** NewSpliceEdges = (DOMAIN_OVERLAP **)malloc(new_splice_edge_cap*sizeof(DOMAIN_OVERLAP *));

  int num_missing_hits = (int)(Graph->MissedHits->N);
  int num_new_edges    = 0;
  int missed_hit_id,node_id;
  for (missed_hit_id = 0; missed_hit_id < num_missing_hits; missed_hit_id++) {


    // Lucky us!  Cheap and easy access!
    P7_ALIDISPLAY * MissedAD = Graph->MissedHits->hit[missed_hit_id]->dcl->ad;


    for (node_id = 1; node_id <= Graph->num_nodes; node_id++) {


      int node_hit_id = Graph->Nodes[node_id]->hit_id;
      int node_dom_id = Graph->Nodes[node_id]->dom_id;


      P7_ALIDISPLAY * NodeAD = (&Graph->TopHits->hit[node_hit_id]->dcl[node_dom_id])->ad;
      if (ExcessiveGapContent(NodeAD))
        continue;


      // Because order matters for 'HitsAreSpliceCompatible' we
      // need to have a catch for either possibility.
      if (HitsAreSpliceCompatible(NodeAD,MissedAD)) {

        NewSpliceEdges[num_new_edges] = (DOMAIN_OVERLAP *)malloc(sizeof(DOMAIN_OVERLAP));
        DOMAIN_OVERLAP * Edge         = NewSpliceEdges[num_new_edges];

        Edge->upstream_hit_id   = node_hit_id;
        Edge->upstream_dom_id   = node_dom_id;
        Edge->downstream_hit_id = missed_hit_id;
        Edge->downstream_dom_id = 0;

        Edge->UpstreamTopHits   = Graph->TopHits;
        Edge->UpstreamDisplay   = NodeAD;
        Edge->DownstreamTopHits = Graph->MissedHits;
        Edge->DownstreamDisplay = MissedAD;

        num_new_edges++;

      } else if (HitsAreSpliceCompatible(MissedAD,NodeAD)) {

        NewSpliceEdges[num_new_edges] = (DOMAIN_OVERLAP *)malloc(sizeof(DOMAIN_OVERLAP));
        DOMAIN_OVERLAP * Edge         = NewSpliceEdges[num_new_edges];

        Edge->upstream_hit_id   = missed_hit_id;
        Edge->upstream_dom_id   = 0;
        Edge->downstream_hit_id = node_hit_id;
        Edge->downstream_dom_id = node_dom_id;

        Edge->UpstreamTopHits   = Graph->MissedHits;
        Edge->UpstreamDisplay   = MissedAD;
        Edge->DownstreamTopHits = Graph->TopHits;
        Edge->DownstreamDisplay = NodeAD;

        num_new_edges++;

      }


      // Time to resize?
      if (num_new_edges == new_splice_edge_cap) {

        new_splice_edge_cap *= 2;

        DOMAIN_OVERLAP ** MoreNewEdges = (DOMAIN_OVERLAP **)malloc(new_splice_edge_cap*sizeof(DOMAIN_OVERLAP *));
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

      if (HitsAreSpliceCompatible(MissedAD,MAD2)) {

        NewSpliceEdges[num_new_edges] = (DOMAIN_OVERLAP *)malloc(sizeof(DOMAIN_OVERLAP));
        DOMAIN_OVERLAP * Edge         = NewSpliceEdges[num_new_edges];

        Edge->upstream_hit_id   = missed_hit_id;
        Edge->upstream_dom_id   = 0;
        Edge->downstream_hit_id = mhi2;
        Edge->downstream_dom_id = 0;

        Edge->UpstreamTopHits   = Graph->MissedHits;
        Edge->UpstreamDisplay   = MissedAD;
        Edge->DownstreamTopHits = Graph->MissedHits;
        Edge->DownstreamDisplay = MAD2;

        num_new_edges++;

      } else if (HitsAreSpliceCompatible(MAD2,MissedAD)) {

        NewSpliceEdges[num_new_edges] = (DOMAIN_OVERLAP *)malloc(sizeof(DOMAIN_OVERLAP));
        DOMAIN_OVERLAP * Edge         = NewSpliceEdges[num_new_edges];

        Edge->upstream_hit_id   = mhi2;
        Edge->upstream_dom_id   = 0;
        Edge->downstream_hit_id = missed_hit_id;
        Edge->downstream_dom_id = 0;

        Edge->UpstreamTopHits   = Graph->MissedHits;
        Edge->UpstreamDisplay   = MAD2;
        Edge->DownstreamTopHits = Graph->MissedHits;
        Edge->DownstreamDisplay = MissedAD;

        num_new_edges++;

      }

      // Time to resize?
      if (num_new_edges == new_splice_edge_cap) {

        new_splice_edge_cap *= 2;

        DOMAIN_OVERLAP ** MoreNewEdges = (DOMAIN_OVERLAP **)malloc(new_splice_edge_cap*sizeof(DOMAIN_OVERLAP *));
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
  
    SketchSpliceEdge(NewSpliceEdges[splice_edge_id],TargetNuclSeq,Graph->Model,gcode);
  
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


  if (DEBUGGING) DEBUG_OUT("'AddMissingExonsToGraph' Complete",-1);


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
  
  if (DEBUGGING) DEBUG_OUT("Starting 'GetExonSetFromEndNode'",1);


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
    BestPathCoords[3] = Graph->MissedHits->hit[Node->hit_id]->dcl->ad->sqto;
    BestPathCoords[4] = Graph->MissedHits->hit[Node->hit_id]->dcl->ad->hmmto;
  } else {
    BestPathCoords[3] = (&Graph->TopHits->hit[Node->hit_id]->dcl[Node->dom_id])->ad->sqto;
    BestPathCoords[4] = (&Graph->TopHits->hit[Node->hit_id]->dcl[Node->dom_id])->ad->hmmto;
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
    BestPathCoords[nodes_in_path*5 + 1] = Graph->MissedHits->hit[Node->hit_id]->dcl->ad->sqfrom;
    BestPathCoords[nodes_in_path*5 + 2] = Graph->MissedHits->hit[Node->hit_id]->dcl->ad->hmmfrom;
  } else {
    BestPathCoords[nodes_in_path*5 + 1] = (&Graph->TopHits->hit[Node->hit_id]->dcl[Node->dom_id])->ad->sqfrom;
    BestPathCoords[nodes_in_path*5 + 2] = (&Graph->TopHits->hit[Node->hit_id]->dcl[Node->dom_id])->ad->hmmfrom;
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


  if (DEBUGGING) DEBUG_OUT("'GetExonSetFromEndNode' Complete",-1);


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
  if (DEBUGGING) DEBUG_OUT("Starting 'FindComponentBestEnd'",1);


  if (ComponentIDs[Node->node_id]) {
    if (DEBUGGING) DEBUG_OUT("'FindComponentBestEnd' Complete",-1);
    return;
  }

  ComponentIDs[Node->node_id] = component_id;

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


  if (DEBUGGING) DEBUG_OUT("'FindComponentBestEnd' Complete",-1);

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

  if (DEBUGGING) DEBUG_OUT("Starting 'TranslateExonSetNucls'",1);


  int translation_len = coding_region_len/3;
  ESL_DSQ * ExonSetTrans = malloc((1 + translation_len) * sizeof(ESL_DSQ));

  int trans_index;
  for (trans_index=1; trans_index<=translation_len; trans_index++) 
      ExonSetTrans[trans_index] = esl_gencode_GetTranslation(gcode,&(ExonSetNucls[3*trans_index-2]));


  if (DEBUGGING) DEBUG_OUT("'TranslateExonSetNucls' Complete",-1);


  return ExonSetTrans;

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: GrabExonCoordSetNucls
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
ESL_DSQ * GrabExonCoordSetNucls
(int * ExonCoordSet, TARGET_SEQ * TargetNuclSeq, int * coding_region_len)
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
    ESL_DSQ * ExonNucls = GrabNuclRange(TargetNuclSeq,range_start,range_end);


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
(int ** ExonCoordSets, int num_exon_sets, TARGET_SEQ * TargetNuclSeq, ESL_GENCODE * gcode)
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
    ESL_DSQ * NuclSeq  = GrabExonCoordSetNucls(ExonCoordSets[exon_set_id],TargetNuclSeq,&num_nucls);
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
 *  Function: CheckTerminusProximity
 *
 *  Desc. : If we're within spitting distance of the ends of the model,
 *          just extend.
 *
 *          NOTE: It seems like in some cases the fwd/bck pipeline is
 *                dead-set on cutting off the final couple of aminos
 *                (Examples: A1BG:{human1,human2,rat1}).
 *
 *                It seems that this check *can* get us the right
 *                mapping coordinates for a full-model hit, but only 
 *                as output from "DumpExonSets," since I can't force
 *                the model to not cut off those last aminos
 *                (without getting into stuff that I'm *not* messing with).
 *
 *                Because having this on might risk creating discrepancies
 *                between our bookkeeping and the real coordinates of the
 *                alignment produced by running the full pipeline,
 *                I'm going to turn it off.
 *
 *                BUT I'll leave it here in case it's ever useful. 
 *
 */
void CheckTerminusProximity
(
  int * ExonCoordSet, 
  TARGET_SEQ * TargetNuclSeq, 
  SPLICE_GRAPH * Graph,
  ESL_GENCODE * gcode
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'CheckTerminusProximity'",1);

  // The formal definition of "spitting distance"
  int max_ext_dist = 5;


  int num_exons = ExonCoordSet[0];

  int nucl_start = ExonCoordSet[1];
  int hmm_start  = ExonCoordSet[2];

  int nucl_end = ExonCoordSet[5*(num_exons-1)+3];
  int hmm_end  = ExonCoordSet[5*(num_exons-1)+4];


  int revcomp = 0;
  if (nucl_start > nucl_end)
    revcomp = 1;


  // Might we consider extending to the model's N-terminus?
  int ext_amino_index;
  if (hmm_start <= max_ext_dist && hmm_start != 1) {

    int ext_len = hmm_start-1;

    int ext_nucl_end = nucl_start;
    int ext_nucl_start;
    if (revcomp) {
      ext_nucl_start = ext_nucl_end + 3*ext_len;
      ext_nucl_end++;
    } else {
      ext_nucl_start = ext_nucl_end - 3*ext_len;
      ext_nucl_end--;
    }

    ESL_DSQ * NExtSeq = GrabNuclRange(TargetNuclSeq,ext_nucl_start,ext_nucl_end);
    int unusual_codon = 0;
    for (ext_amino_index=0; ext_amino_index<ext_len; ext_amino_index++) {
      int n_ext_amino = esl_gencode_GetTranslation(gcode,&NExtSeq[3*ext_amino_index+1]);
      if (n_ext_amino < 0 || n_ext_amino > 20) {
        unusual_codon = 1;
        break;
      }
    }
    free(NExtSeq);

    if (!unusual_codon) {
      ExonCoordSet[1] = ext_nucl_start;
      ExonCoordSet[2] = 1;
    }

  }



  // How about that C-terminus?
  if (Graph->Model->M - hmm_end <= max_ext_dist && hmm_end != Graph->Model->M) {

    int ext_len = Graph->Model->M - hmm_end;

    int ext_nucl_start = nucl_end;
    int ext_nucl_end;
    if (revcomp) {
      ext_nucl_end = ext_nucl_start - 3*ext_len;
      ext_nucl_start--;
    } else {
      ext_nucl_end = ext_nucl_start + 3*ext_len;
      ext_nucl_start++;
    }

    ESL_DSQ * CExtSeq = GrabNuclRange(TargetNuclSeq,ext_nucl_start,ext_nucl_end);
    int unusual_codon = 0;
    for (ext_amino_index=0; ext_amino_index<ext_len; ext_amino_index++) {
      int c_ext_amino = esl_gencode_GetTranslation(gcode,&CExtSeq[3*ext_amino_index+1]);
      if (c_ext_amino < 0 || c_ext_amino > 20) {
        unusual_codon = 1;
        break;
      }
    }
    free(CExtSeq);

    if (!unusual_codon) {
      ExonCoordSet[5*(num_exons-1)+3] = ext_nucl_end;
      ExonCoordSet[5*(num_exons-1)+4] = Graph->Model->M;
    }

  }


  if (DEBUGGING) DEBUG_OUT("'CheckTerminusProximity' Complete",-1);

}








/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: ExonSetCleanup
 *
 */
int ** ExonSetCleanup
(SPLICE_GRAPH * Graph, int * InputCoordSet, int * num_split_sets)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'ExonSetCleanup'",1);


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


    if (DEBUGGING) DEBUG_OUT("'ExonSetCleanup' Complete",-1);


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


    while (InputCoordSet[5*scanner_exon_id+5] == -1 && scanner_exon_id < num_exons)
      scanner_exon_id++;

    if (scanner_exon_id == num_exons)
      break;

    int start_exon_id = scanner_exon_id;


    while (InputCoordSet[5*scanner_exon_id+5] != -1 && scanner_exon_id < num_exons)
      scanner_exon_id++;

    int end_exon_id = scanner_exon_id-1;


    // If we aren't starting with the first exon, we need to recover the
    // start bounds of the actual hit (as opposed to the now-jettisoned splice
    // coordinates)
    if (start_exon_id > 0) {


      SPLICE_NODE * Node = Graph->Nodes[InputCoordSet[5*start_exon_id+5]];
      P7_ALIDISPLAY * AD;
      if (Node->was_missed)
        AD = Graph->MissedHits->hit[Node->hit_id]->dcl->ad;
      else
        AD = (&Graph->TopHits->hit[Node->hit_id]->dcl[Node->dom_id])->ad;


      int exon_hmm_start  = AD->hmmfrom;
      int exon_nucl_start = AD->sqfrom;

      InputCoordSet[5*start_exon_id+1] = exon_nucl_start;
      InputCoordSet[5*start_exon_id+2] = exon_hmm_start;


    }


    // Similarly, if we aren't ending with the last exon we need to go back
    // to what the hit tells us.
    if (end_exon_id < num_exons-1) {


      SPLICE_NODE * Node = Graph->Nodes[InputCoordSet[5*end_exon_id+5]];
      P7_ALIDISPLAY * AD;
      if (Node->was_missed)
        AD = Graph->MissedHits->hit[Node->hit_id]->dcl->ad;
      else
        AD = (&Graph->TopHits->hit[Node->hit_id]->dcl[Node->dom_id])->ad;


      int exon_hmm_end  = AD->hmmto;
      int exon_nucl_end = AD->sqto;

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


  if (DEBUGGING) DEBUG_OUT("'ExonSetCleanup' Complete",-1);


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

  if (DEBUGGING) DEBUG_OUT("Starting 'GetSplicedExonCoordSets'",1);


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


  if (DEBUGGING) DEBUG_OUT("'GetSplicedExonCoordSets' Complete",-1);


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
 *  Function: PrintExon
 *
 */
void PrintExon
(
  EXON_DISPLAY_INFO * EDI, 
  TARGET_SEQ * TargetNuclSeq, 
  int * ad_nucl_read_pos,
  int * ad_amino_read_pos,
  int * codon_pos
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'PrintExon'",1);


  FILE          * ofp = EDI->ofp;
  P7_ALIDISPLAY * AD  = EDI->AD;
  
  ESL_DSQ * LeftDinucls  = NULL;
  ESL_DSQ * RightDinucls = NULL;
  int ali_nucl_start = EDI->nucl_start;
  if (EDI->revcomp) {

    if (!EDI->Nterm) {
      LeftDinucls = GrabNuclRange(TargetNuclSeq,EDI->nucl_start+2,EDI->nucl_start+1);
      ali_nucl_start += 2;
    }

    if (!EDI->Cterm) {
      RightDinucls = GrabNuclRange(TargetNuclSeq,EDI->nucl_end-1,EDI->nucl_end-2);
    }

  } else {

    if (!EDI->Nterm) {
      LeftDinucls = GrabNuclRange(TargetNuclSeq,EDI->nucl_start-2,EDI->nucl_start-1);
      ali_nucl_start -= 2;
    }

    if (!EDI->Cterm) {
      RightDinucls = GrabNuclRange(TargetNuclSeq,EDI->nucl_end+1,EDI->nucl_end+2);
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
  int exon_ali_str_alloc = 2 * abs(EDI->nucl_end - EDI->nucl_start);
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

  if (DEBUGGING) DEBUG_OUT("'PrintExon' Complete",-1);


}





/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: PrintSplicedAlignment
 *
 *  Desc. :
 *
 *  Inputs:  
 *
 *  Output:
 *
 */
void PrintSplicedAlignment
(
  P7_ALIDISPLAY * AD, 
  TARGET_SEQ    * TargetNuclSeq, 
  int           * ExonCoordSet, 
  int             exon_set_name_id, 
  FILE          * ofp, 
  int             textw
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'PrintSplicedAlignment'",1);

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

    PrintExon(EDI,TargetNuclSeq,&ad_nucl_read_pos,&ad_amino_read_pos,&codon_pos);

  }


  free(EDI->HMMName);
  free(EDI->TransName);
  free(EDI->NuclName);
  free(EDI->NameBlank);
  free(EDI->CoordBlank);
  free(EDI);


  if (DEBUGGING) DEBUG_OUT("'PrintSplicedAlignment' Complete",-1);

}










/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: DumpSplashHeader / DumpSplashFooter
 *
 */
void DumpSplashHeader
(
  SPLICE_GRAPH * Graph,
  TARGET_SEQ * TargetNuclSeq,
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
  fprintf(ofp,"| = Target Seq Name %s\n",TargetNuclSeq->SeqName);
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

  if (DEBUGGING) DEBUG_OUT("Starting 'DetermineHitExonCoords'",1);


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


  if (DEBUGGING) DEBUG_OUT("'DetermineHitExonCoords' Complete",-1);


  return HitExonCoords;

}









/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: ReportSplicedTopHits
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
void ReportSplicedTopHits
(
  SPLICE_GRAPH * Graph,
  P7_TOPHITS   * ExonSetTopHits, 
  P7_PIPELINE  * ExonSetPipeline, 
  TARGET_SEQ   * TargetNuclSeq, 
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
  

  if (DEBUGGING) p7_tophits_Domains(ofp, ExonSetTopHits, ExonSetPipeline, textw);


  int hit_id, dom_id;
  for (hit_id = 0; hit_id < (int)(ExonSetTopHits->N); hit_id++) {
    for (dom_id = 0; dom_id < ExonSetTopHits->hit[hit_id]->ndom; dom_id++) {


      P7_ALIDISPLAY * AD  = (&ExonSetTopHits->hit[hit_id]->dcl[dom_id])->ad;
      if (ExcessiveGapContent(AD))
        continue;


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



      // If Alex is running this, he probably wants to
      // make the output as cluttered as possible (and
      // we all love that about him... right?).
      //
      if (ALEX_MODE) DumpSplashHeader(Graph,TargetNuclSeq,*exon_set_name_id,HitExonCoords,ofp,textw);


      PrintSplicedAlignment(AD,TargetNuclSeq,HitExonCoords,*exon_set_name_id,ofp,textw);


      // You thought Alex was done making a mess of things?! HA!
      if (ALEX_MODE) DumpSplashFooter(ofp,textw);


      free(HitExonCoords);

    }
  }

}











/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *  Function: RunModelOnExonSets
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
void RunModelOnExonSets
(
  SPLICE_GRAPH * Graph, 
  TARGET_SEQ   * TargetNuclSeq, 
  ESL_GENCODE  * gcode, 
  ESL_GETOPTS  * go,
  FILE         * ofp,
  int            textw
)
{

  if (DEBUGGING) DEBUG_OUT("Starting 'RunModelOnExonSets'",1);


  int num_coord_sets;
  int ** ExonCoordSets = GetSplicedExonCoordSets(Graph,&num_coord_sets);


  if (DEBUGGING) DumpExonSets(ExonCoordSets,num_coord_sets,TargetNuclSeq,gcode);


  // It's better to be re-using these than destroying
  // and re-allocating every time
  ESL_SQ      * NuclSeq           = NULL; 
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
    ESL_DSQ * ExonSetNucls = GrabExonCoordSetNucls(ExonCoordSets[coord_set_id],TargetNuclSeq,&coding_region_len);
    ESL_SQ  * NuclSeq      = esl_sq_CreateDigitalFrom(TargetNuclSeq->abc,"Exon Set",ExonSetNucls,(int64_t)coding_region_len,NULL,NULL,NULL);
    NuclSeq->idx = coord_set_id+1;
    esl_sq_Textize(NuclSeq);


    // Translate the nucleotide sequence and convert to
    // a digitized ESL_SQ
    //
    int trans_len = coding_region_len / 3;
    ESL_DSQ * ExonSetTrans = TranslateExonSetNucls(ExonSetNucls,coding_region_len,gcode);
    ESL_SQ  * AminoSeq     = esl_sq_CreateDigitalFrom(Graph->OModel->abc,TargetNuclSeq->SeqName,ExonSetTrans,(int64_t)trans_len,NULL,NULL,NULL);
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
    //
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
      ReportSplicedTopHits(Graph,ExonSetTopHits,ExonSetPipeline,TargetNuclSeq,ExonCoordSets[coord_set_id],&exon_set_id,ofp,textw);



    // DESTRUCTION AND REBIRTH!
    free(ExonSetNucls);
    free(ExonSetTrans);
    esl_sq_Destroy(NuclSeq);
    esl_sq_Destroy(AminoSeq);
    p7_tophits_Destroy(ExonSetTopHits);
    p7_pipeline_Destroy(ExonSetPipeline);
    free(ExonCoordSets[coord_set_id]);

  }


  // That's all we needed!  Good work, team!
  if (ExonSetTopHits) p7_tophits_Destroy(ExonSetTopHits);
  if (ExonSetPipeline) p7_pipeline_Destroy(ExonSetPipeline);
  p7_oprofile_Destroy(ExonSetOModel);
  p7_bg_Destroy(ExonSetBackground);
  free(ExonCoordSets);
  if (DEBUGGING) DEBUG_OUT("'RunModelOnExonSets' Complete",-1);

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
 *           4.             om : The optimized profile (assumed to be built on 'gm').
 *           5.          gcode : An ESL_GENCODE struct (mainly used for translation).
 *           6.             go : An ESL_GETOPTS struct (in case of options... duh).
 *           7.            ofp : An open file pointer (specificially, the target for output).
 *           8.          textw : The desired line width for output.
 *
 *  Output:  Nothing is returned, but if we find a way to splice together some of the
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

  if (DEBUGGING) DEBUG_OUT("Starting 'SpliceHits'",1);


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
  TARGET_SET * TargetSet = SelectTargetRanges(TopHits);



  int target_set_id;
  for (target_set_id = 0; target_set_id < TargetSet->num_target_seqs; target_set_id++) {



    // Given that our hits are organized by target sequence, we can
    // be a bit more efficient in our file reading by only pulling
    // target sequences as they change (wrt the upstream hit)
    //
    TARGET_SEQ * TargetNuclSeq = GetTargetNuclSeq(GenomicSeqFile,TargetSet,target_set_id);



    // This function encapsulates a *ton* of the work we do.
    // In short, take the collection of unspliced hits and build
    // a splice graph representing all (reasonable) ways of splicing
    // them.
    //
    SPLICE_GRAPH * Graph = BuildSpliceGraph(TopHits,TargetNuclSeq,gm,om,gcode);



    // Evaluate the graph for any holes (or possible holes)
    // worth plugging
    //
    AddMissingExonsToGraph(Graph,TargetNuclSeq,gcode);



    // If we're debugging, it might be useful to get a quick
    // picture of what our final splice graph looks like.
    //
    if (DEBUGGING) DumpGraph(Graph);



    // Re-run the model on the extracted nucleotide sequence(s)
    // for each connected component of the graph.
    //
    RunModelOnExonSets(Graph,TargetNuclSeq,gcode,go,ofp,textw);



    // CLEANUP
    SPLICE_GRAPH_Destroy(Graph);
    TARGET_SEQ_Destroy(TargetNuclSeq);


  }


  // More cleanup!
  TARGET_SET_Destroy(TargetSet);



  if (DEBUGGING) {
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
/*
 *                                                                                              END SPLICING STUFF
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */




