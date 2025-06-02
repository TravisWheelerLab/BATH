/* Structs and MACROs for bilding a splice graph */

#include "p7_config.h"

#include <string.h>

#include "easel.h"
#include "hmmer.h"

typedef struct _splice_edge {

  int frameshift;
  int bypass_checked;

  int upstream_node_id;
  int downstream_node_id;

  int overlap_amino_start;
  int overlap_amino_end;

  int upstream_trace_start;
  int upstream_trace_end;

  int upstream_nuc_start;    //sub-seq coords
  int upstream_nuc_end;      //sub-seq coords

  int downstream_trace_start;
  int downstream_trace_end;

  int downstream_nuc_start;  //sub-seq coords
  int downstream_nuc_end;    //sub-seq coords
  
  int upstream_ext_len;
  int downstream_ext_len;

  int upstream_spliced_amino_end;      
  int downstream_spliced_amino_start; 

  int upstream_spliced_nuc_end;      //true seq coords
  int downstream_spliced_nuc_start;  //true seq coords

  int target_seq_start;
  int target_seq_end;
  int target_seq_n;

  float splice_score;
  float signal_score;

  struct _splice_edge *prev;
  struct _splice_edge *next;

} SPLICE_EDGE;


typedef struct _splice_graph {

  /* Graph size intfo */
  int nalloc;
  int num_nodes;
  int num_edges;
  int orig_N;
  int split_N;
  int recover_N; 
 
  /* Target sequence info */
  int          revcomp;
  int64_t      seqidx;      
  char        *seqname;    

  /* Hits and hit info */
  int         *node_in_graph;        /* Is the hit part of the current graph */
  int         *reportable;      /* For orignal hits, do they pass the repoting threshold */
  int         *orig_hit_idx;    /* index of hits in original P7_TOPHITS  */
  int         *split_orig_id;
  

  /*Edge info */  
  int   *edge_in_graph;
  int   *best_out_edge;
  int   *best_in_edge;
  
  /* Scores */
  float *path_scores;  //Path score pulled upstream
  float *ali_scores;

  P7_TOPHITS  *th;  
  SPLICE_EDGE **edges;

} SPLICE_GRAPH;



typedef struct _splice_path {

  int revcomp;
  int frameshift;

  int alloc_len;
  int path_len;
  int seq_len;

  int *node_id;
  int *split;

  int *upstream_spliced_amino_end;
  int *downstream_spliced_amino_start;

  int *upstream_spliced_nuc_end;
  int *downstream_spliced_nuc_start;

  float *hit_scores;
  float *edge_scores;   // edge_scores[i] is score for incoming edge to hits[i]
  float *signal_scores;
  
  P7_HIT **hits;

} SPLICE_PATH;


typedef struct _splice_pipeline 
{
  
  int      by_E;
  int      inc_by_E;
  int      do_null2;
  int      do_biasfilter;
  int      long_targets;
  int      frameshift;

  double   E;
  double   T;
  double   Z;
  double   F1;  
  double   F2;
  double   F3;
  double   incE;
  double   incT;

  float   *signal_scores;

  int     *orig_nuc_idx;

  ESL_SQ  *nuc_sq;
  ESL_SQ  *amino_sq; 

  P7_OMX  *fwd;
  P7_OMX  *bwd;  

  P7_GMX  *gfwd;
  P7_GMX  *gbwd;

  P7_BG   *bg;

  P7_HIT  *hit;
 
} SPLICE_PIPELINE;

typedef struct _splice_gap
{

  int up_node;
  int down_node;  

  int hmm_min;
  int hmm_max;
  int hmm_len;
  int seq_min;
  int seq_max;
  int seq_len;

} SPLICE_GAP;



#define MAX_GAP_RANGE             50000       
#define MAX_INTRON_SHORT          10000 
#define MAX_INTRON_LONG           100000
#define TERMINAL_EXT              10000
#define MIN_INTRON_LEN            13
#define MAX_AMINO_EXT             50
#define MIN_AMINO_OVERLAP         10
#define MAX_AMINO_OVERLAP         20
#define MIN_TERMINAL_LEN          5


/* Indices of p7_splice_SignalScores */
enum p7s_splice_signals_e {
  p7S_GTAG  = 0,
  p7S_GCAG  = 1,
  p7S_ATAC  = 2,
  p7S_OTHER = 3
};
#define p7S_SPLICE_SIGNALS 4

/* p7_spliceedge.c */
extern SPLICE_EDGE* p7_spliceedge_Create(void);
extern SPLICE_EDGE* p7_spliceedge_ConnectHits(SPLICE_PIPELINE *pli, const P7_DOMAIN *upstream_domain, const P7_DOMAIN *downstream_domain, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq, int revcomp);
extern SPLICE_EDGE* p7_spliceedge_ConnectSplits(SPLICE_PIPELINE *pli, const P7_DOMAIN *upstream_domain, const P7_DOMAIN *downstream_domain, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq, int frameshift, int revcomp);


/* p7_splicegap.c */
extern SPLICE_GAP* p7_splicegap_Create(void);
extern P7_HIT** p7_splicegap_AlignGap(SPLICE_GRAPH *graph, SPLICE_GAP *gap, const P7_HMM *hmm, const P7_BG *bg, const ESL_GENCODE *gcode, const ESL_SQFILE *seq_file, int *num_hits);

/* p7_splicegraph.c */
extern SPLICE_GRAPH* p7_splicegraph_Create(void);
extern int p7_splicegraph_CreateNodes(SPLICE_GRAPH *graph, int num_nodes);
extern int p7_splicegraph_Grow(SPLICE_GRAPH *graph);
extern void p7_splicegraph_Destroy(SPLICE_GRAPH *graph);
extern int p7_splicegraph_AddNode(SPLICE_GRAPH *graph, P7_HIT *hit);
extern int p7_splicegraph_AddEdge(SPLICE_GRAPH *graph, SPLICE_EDGE *edge);
extern int p7_splicegraph_EdgeExists(SPLICE_GRAPH* graph, int up_node, int down_node);
extern SPLICE_EDGE* p7_splicegraph_GetEdge(SPLICE_GRAPH* graph, int up_node, int down_node);
extern int p7_splicegraph_PathExists (SPLICE_GRAPH *graph, int upstream_node, int downstream_node);
extern void p7_splicegraph_DumpHits(FILE *fp, SPLICE_GRAPH *graph);
extern void p7_splicegraph_DumpEdges(FILE *fp, SPLICE_GRAPH *graph);
extern void p7_splicegraph_DumpGraph(FILE *fp, SPLICE_GRAPH *graph, int show_j);

/* p7_splicehits.c */
extern SPLICE_SAVED_HITS* p7_splicehits_CreateSavedHits(void);
extern int p7_splicehits_GrowSavedHits(SPLICE_SAVED_HITS *saved_hits);
extern void p7_splicehits_DestroySavedHits(SPLICE_SAVED_HITS *saved_hits);
extern int p7_splicehits_CreateNext(SPLICE_SAVED_HITS *saved_hits, SPLICE_HIT_INFO **ret_info);
extern int p7_splicehits_SortSavedHits(SPLICE_SAVED_HITS *sh);
extern int p7_splicehits_MergeSavedHits(SPLICE_SAVED_HITS *sh1, SPLICE_SAVED_HITS *sh2);
extern int p7_splicehits_FindNodes(SPLICE_GRAPH *graph, SPLICE_SAVED_HITS *sh, int first, int last);
extern int p7_splice_RemoveDuplicates(SPLICE_SAVED_HITS *sh);
extern void p7_splicehits_Dump(FILE *fp, SPLICE_SAVED_HITS *sh);
 
/* p7_splicepath.c */
extern SPLICE_PATH* p7_splicepath_Create(int path_len);
extern int p7_splicepath_Grow(SPLICE_PATH *path);
extern void p7_splicepath_Destroy(SPLICE_PATH *path);
extern SPLICE_PATH* p7_splicepath_GetBestPath(SPLICE_GRAPH *graph);
extern void p7_splicepath_Dump(FILE *fp, SPLICE_PATH *path);


/* p7_splicepipeline.c */
extern SPLICE_PIPELINE* p7_splicepipeline_Create(const ESL_GETOPTS *go, int M_hint, int L_hint);
extern void p7_splicepipeline_Reuse(SPLICE_PIPELINE *pli);
extern void p7_splicepipeline_Destroy(SPLICE_PIPELINE *pli); 


/* p7_splice.c */
extern int p7_splice_AddOriginals(SPLICE_GRAPH *graph, const P7_TOPHITS *tophits, int *hits_processed, int64_t seqidx);
extern int p7_splice_SplitHits(SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, const P7_HMM *hmm, P7_FS_PROFILE *gm_fs, ESL_GENCODE *gcode, const ESL_SQFILE *seq_file, int *hits_to_process);
extern int p7_splice_RecoverHits(SPLICE_GRAPH *graph, SPLICE_SAVED_HITS *saved_hits, SPLICE_PIPELINE *pli, P7_HMM *hmm, ESL_GENCODE *gcode, ESL_SQFILE *seq_file, int first, int last);
extern int p7_splice_ConnectGraph(SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQFILE *seq_file);
extern int p7_splice_FillGaps(SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQFILE *seq_file);
extern int p7_splice_MergePaths(SPLICE_GRAPH *graph, SPLICE_PATH **path_accumulator, SPLICE_PIPELINE *pli, P7_HMM *hmm, ESL_GENCODE *gcode, ESL_SQFILE *seq_file, int *num_paths);
extern int p7_splice_FindTerminals(SPLICE_GRAPH *graph, SPLICE_PATH **path_accumulator, SPLICE_PIPELINE *pli, P7_HMM *hmm, ESL_GENCODE *gcode, ESL_SQFILE *seq_file, int *num_paths);
extern int p7_splice_AlignPath(SPLICE_GRAPH *graph, SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_TOPHITS *tophits, P7_OPROFILE *om, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *path_seq, int64_t db_nuc_cnt, float fs_prob);
extern int p7_splice_AlignFrameshiftPath(SPLICE_GRAPH *graph, SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_TOPHITS *tophits, P7_FS_PROFILE *gm_fs, ESL_GENCODE *gcode, ESL_SQ *path_seq, int64_t db_nuc_cnt);
extern int p7_splice_ReleaseHits(SPLICE_GRAPH *graph, int *hits_processed, int *num_hits_processed, int range_bound_min, int range_bound_maxs);
extern int p7_splice_EnforceRangeBounds(SPLICE_GRAPH *graph, int64_t bound_min, int64_t bound_max);
extern ESL_SQ* p7_splice_GetSubSequence(const ESL_SQFILE *seq_file, char* seqname, int64_t seq_min, int64_t seq_max, int revcomp);
extern P7_HMM* p7_splice_GetSubHMM (const P7_HMM *hmm, int start, int end);
