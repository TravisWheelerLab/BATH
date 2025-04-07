/* Structs and MACROs for bilding a splice graph */

#include "p7_config.h"

#include <string.h>

#include "easel.h"
#include "hmmer.h"


typedef struct _target_range {

  int          revcomp;      /* reverse complementarity               */

  int64_t      start;        /* range start position                  */
  int64_t      end;          /* range end postion                     */

  int64_t      seqidx;       /* target squence id                     */
  char        *seqname;      /* target sequnce name                   */

  int          seed_hit_idx;    /* <th> index of the seed hit            */
  int         *in_target_range; /* Is the hit part of the current target range */
  int         *reportable;      /* For orignal hits, do they pass the repoting threshold */
  int          orig_N;          /* the number of original hits in the th */ 
  int         *orig_hit_idx;    /* index of hits in original P7_TOPHITS  */

  P7_TOPHITS  *th;             /* all hits on current sequence and strand */                        

} TARGET_RANGE;

typedef struct _splice_downstream_nodes {

  int node_id;
  int gap_checked;
  int edge_exists;
  int path_exists;

  struct _splice_downstream_nodes *next;

} SPLICE_DS_NODE;


typedef struct _splice_edge {

  int frameshift;

  int upstream_node_id;
  int downstream_node_id;

  int overlap_amino_start;
  int overlap_amino_end;

  int upstream_trace_start;
  int upstream_trace_end;

  int upstream_nuc_start;
  int upstream_nuc_end;

  int downstream_trace_start;
  int downstream_trace_end;

  int downstream_nuc_start;
  int downstream_nuc_end;
  
  int upstream_ext_len;
  int downstream_ext_len;

  int upstream_spliced_amino_end;
  int downstream_spliced_amino_start;

  int upstream_spliced_nuc_end;
  int downstream_spliced_nuc_start;

  int target_seq_start;
  int target_seq_end;
  int target_seq_n;

  float splice_score;
  float signal_score;

  struct _splice_edge *next;

} SPLICE_EDGE;


typedef struct _splice_graph {

  int nalloc;  

  int revcomp;

  int num_nodes;
  int num_edges;
  int orig_num_nodes;
 
  int   best_path_length;
  int   best_path_start;
  int   best_path_end;

  int   *out_edge_cnt;
  int   *in_edge_cnt;
  int   *best_out_edge;
  int   *best_in_edge;

  float *path_scores;  //Path score pulled upstream
  float *hit_scores;

  SPLICE_DS_NODE **ds_nodes;
  SPLICE_EDGE    **edges;

} SPLICE_GRAPH;



typedef struct _splice_path {

  int revcomp;

  int path_len;
  int seq_len;
  int split_hits;

  int *node_id;

  int *missing;

  int *upstream_spliced_amino_end;
  int *downstream_spliced_amino_start;

  int *upstream_spliced_nuc_end;
  int *downstream_spliced_nuc_start;
 
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
  int upstream_node;
  int downstream_node;

  int hmm_start;
  int hmm_end;
  int seq_start;
  int seq_end;

} SPLICE_GAP;

#define MAX_GAP_RANGE             200000    //2X10^5 
#define MAX_INTRON_LEN            50000      //5x10^4 
#define MIN_INTRON_LEN            10
#define MAX_AMINO_EXT             10
#define MIN_AMINO_OVERLAP         6
#define MAX_AMINO_OVERLAP         12

/* Indices of p7_splice_SignalScores */
enum p7s_splice_signals_e {
  p7S_GTAG  = 0,
  p7S_GCAG  = 1,
  p7S_ATAC  = 2,
  p7S_OTHER = 3
};
#define p7S_SPLICE_SIGNALS 4


static ESL_OPTIONS Translation_Options[] = {
  /* name     type         default env_var  range toggles req  incompat help                  docgroup */
 { "--watson",eslARG_NONE, FALSE,   NULL,    NULL, NULL,   NULL, NULL,   "only translate top strand",        99 },
 { "--crick", eslARG_NONE, FALSE,  NULL,    NULL, NULL,   NULL, NULL,   "only translate bottom strand",     99 },
 { "-l",      eslARG_INT,  "5",    NULL,    NULL, NULL,   NULL, NULL,   "minimum ORF length",               99 },
 { "-m",      eslARG_NONE, FALSE,  NULL,    NULL, NULL,   NULL,"-M",    "ORFs must initiate with AUG only", 99 },
 { "-M",      eslARG_NONE, FALSE,  NULL,    NULL, NULL,   NULL,"-m",    "ORFs must start with allowed initiation codon", 99 },

 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* Allocation and Destruction */
extern TARGET_RANGE* target_range_create(int nalloc);
extern int target_range_grow(TARGET_RANGE *target_range);
extern void target_range_destroy(TARGET_RANGE *target_range);
extern SPLICE_EDGE* splice_edge_create(void);
extern SPLICE_GRAPH* splice_graph_create(void);
extern void splice_graph_destroy (SPLICE_GRAPH* graph);
extern SPLICE_PATH* splice_path_create(int path_len);
extern int splice_path_split_hit(SPLICE_PATH *path, SPLICE_EDGE *edge, int split_id);
extern void splice_path_destroy(SPLICE_PATH *path);
extern SPLICE_PIPELINE* splice_pipeline_create(const ESL_GETOPTS *go, int M_hint, int L_hint);
extern void splice_pipeline_destroy(SPLICE_PIPELINE *pli);


/* Target Range  */
extern int* build_target_range (TARGET_RANGE *target_range, const P7_TOPHITS *tophits, int *hits_processed, int64_t seqidx, int revcomp, int *reportable); 
extern ESL_SQ* get_target_range_sequence(const TARGET_RANGE *target_range, const ESL_SQFILE *seq_file);
extern ESL_SQ* get_sub_sequence(const ESL_SQFILE *seq_file, char* seqname, int64_t seq_min, int64_t seq_max, int revcomp);
extern int release_hits_from_target_range(const TARGET_RANGE *target_range, int *hit_processed, int *num_hits_proccesed, int range_bound_min, int range_bound_max);


/* Initial Splice Graph */
extern int fill_graph_with_nodes(SPLICE_GRAPH *graph, TARGET_RANGE *target_range, P7_PROFILE *gm);
extern SPLICE_EDGE* connect_nodes_with_edges(const P7_HIT *upstream_hit, const P7_HIT *downstream_hit, const P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, const ESL_GENCODE *gcode, const ESL_SQ *target_seq, int revcomp);

extern void get_overlap_nuc_coords (SPLICE_EDGE *edge, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const ESL_SQ *target_seq, int revcomp);
extern int find_optimal_splice_site (SPLICE_EDGE *edge, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, const ESL_GENCODE *gcode, const ESL_SQ *target_seq);
extern float ali_score_at_postion (P7_PROFILE *gm, int amino, int model_pos, int trans_pos, int prev_state, int curr_state);
extern int select_splice_option (SPLICE_EDGE *edge, const P7_PROFILE *gm, P7_PROFILE *sub_model, P7_FS_PROFILE *sub_fs_model, const ESL_GENCODE *gcode, const ESL_SQ *target_seq, float signal_score, int up_nuc_pos, int down_nuc_pos);
extern int add_edge_to_graph(SPLICE_GRAPH *graph, SPLICE_EDGE *edge);
extern void p7_splice_SignalScores(float *f);

/* Missing Exons */
extern int fill_holes_in_graph(TARGET_RANGE *target_range, SPLICE_GRAPH *graph, const P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, ESL_GENCODE *gcode, const ESL_SQFILE *seq_file, int** gap_checked);
extern int path_exists (SPLICE_GRAPH *graph, int upstream_node, int downstream_node, int *visited);
extern SPLICE_GAP* find_the_gap (P7_TOPHITS *th, const P7_PROFILE *gm, ESL_SQ *target_seq, int orig_N, int upstream_node, int downstream_node, int revcomp);
extern P7_HIT** align_the_gap(TARGET_RANGE *target_range, const P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, ESL_SQ *target_seq, ESL_GENCODE *gcode, SPLICE_GAP *gap, int *num_hits, int revcomp);
extern P7_HMM* extract_sub_hmm (const P7_HMM *hmm, int start, int end);
extern ESL_DSQ* extract_sub_seq(ESL_SQ *target_seq, int start, int end, int revcomp);
extern int add_missed_hit_to_target_range(TARGET_RANGE *target_range, P7_HIT *hit);
extern int bridge_the_gap(TARGET_RANGE *target_range, SPLICE_GRAPH *graph, P7_HIT **top_ten, const P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, const ESL_GENCODE *gcode, const ESL_SQ *target_seq, int up_gap, int down_gap, int num_hits);


/* Splice Path */
extern int reset_edge_counts(TARGET_RANGE *target_range, SPLICE_GRAPH *graph);
extern int check_for_loops (SPLICE_GRAPH *graph, P7_TOPHITS *th);
extern int enforce_range_bounds(SPLICE_GRAPH *graph, P7_TOPHITS *th, int64_t* range_bound_mins, int64_t* range_bound_maxs, int range_cnt);
extern SPLICE_PATH* evaluate_paths (TARGET_RANGE *target_range, SPLICE_GRAPH *graph);
extern int longest_path_upstream (TARGET_RANGE *target_range, SPLICE_GRAPH *graph);
extern int topological_sort_upstream (TARGET_RANGE *target_range, SPLICE_GRAPH *graph, int *visited, int *stack, int *stack_size, int node);
extern int split_hits_in_path (SPLICE_GRAPH *graph, SPLICE_PATH *path, const P7_PROFILE *gm, const P7_HMM *hmm, const P7_BG *bg, const ESL_SQFILE *seq_file, const ESL_GENCODE *gcode, char *seqname, int orig_N);

/* Spliced Hit Processing */
extern int splice_path (SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_TOPHITS *orig_tophits, P7_OPROFILE *om, P7_PROFILE *gm, const ESL_SQFILE *seq_file, ESL_GENCODE *gcode, int64_t db_nuc_cnt, char *seqname, int orig_N, int* success); 
extern int splice_path_frameshift (SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_FS_PROFILE *gm_fs, P7_OPROFILE *om, P7_SCOREDATA *scoredata, const ESL_SQFILE *seq_file, ESL_GENCODE *gcode, int64_t db_nuc_cnt, char* seqname, int orig_N, int* success); 
extern int align_spliced_path (SPLICE_PIPELINE *pli, P7_OPROFILE *om, P7_PROFILE *gm, ESL_SQ *target_seq, ESL_GENCODE *gcode);
extern int align_spliced_path_frameshift (SPLICE_PIPELINE *pli, P7_FS_PROFILE *gm_fs, P7_OPROFILE *om, P7_SCOREDATA *scoredata, ESL_SQ *target_seq, ESL_GENCODE *gcode);
extern int p7_splice_compute_ali_scores(P7_DOMAIN *dom, P7_TRACE *tr, ESL_DSQ *amino_dsq, const P7_PROFILE *gm, int K);
extern int p7_splice_compute_ali_scores_fs(P7_DOMAIN *dom, P7_TRACE *tr, ESL_DSQ *nuc_dsq, P7_FS_PROFILE *gm_fs, const ESL_ALPHABET *abc);

/* Debug Dumps */
extern void target_range_dump(FILE *fp, TARGET_RANGE *target_range, int print_hits);
extern void graph_dump(FILE *fp, SPLICE_GRAPH *graph, int print_edges);
extern void path_dump(FILE *fp, SPLICE_PATH *path);



