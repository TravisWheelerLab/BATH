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

  int          orig_N;       /* the number of original hits in the th */ 
  int         *orig_hit_idx; /* index of hits in original P7_TOPHITS  */
  P7_TOPHITS  *th;           /* hits in range                         */

} TARGET_RANGE;


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

  float splice_score;
  float signal_score;

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
  int   *orig_out_edge;
  int   *orig_in_edge;
  int   *best_out_edge;
  int   *best_in_edge;

  int **is_upstream;     //[upstream][downstream]
  int  *is_upstream_mem;

  int **edge_id;      //[upstream][downstream]
  int *edge_id_mem;

  float *path_scores;  //Path score pulled upstream
  float *hit_scores;

  float **edge_scores; //[upstream][downstream]
  float *edge_scores_mem;

  SPLICE_EDGE **edges;

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

#define MAX_TARGET_RANGE_EXT      1000000    //1x10^6 
#define MAX_INTRON_LEN            50000      //5x10^4 
#define MIN_INTRON_LEN            10
#define MAX_AMINO_EXT             10
#define MIN_AMINO_OVERLAP         6
#define MAX_AMINO_OVERLAP         12

/* Splice singal probabilities taken from 
 * "Comprehensive splice-site analysis using comparative genomics", 
 * Nihar Sheth et al., 2006 */
extern void
p7_splice_SignalScores(float *f)
{
  f[0] = log(0.9919);     /* GT-AG */
  f[1] = log(0.0073);     /* GC-AG */
  f[2] = log(0.0006);     /* AT-AC */
  f[3] = log(0.0002);     /* OTHER */
  return;
}

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
extern int splice_graph_grow(SPLICE_GRAPH *graph);
extern int splice_graph_create_nodes(SPLICE_GRAPH *graph, int num_nodes);
extern SPLICE_PATH* splice_path_create(int path_len);
extern int splice_path_split_hit(SPLICE_PATH *path, SPLICE_EDGE *edge, int split_id);
extern void splice_path_destroy(SPLICE_PATH *path);
extern SPLICE_PIPELINE* splice_pipeline_create(const ESL_GETOPTS *go, int M_hint, int L_hint);
extern void splice_pipeline_destroy(SPLICE_PIPELINE *pli);


/* Target Range  */
extern TARGET_RANGE* build_target_range(const P7_TOPHITS *tophits, int *hits_processed, int *num_hits_processed, int64_t *range_bound_mins, int64_t *range_bound_maxs, int range_num, int64_t seqidx, int revcomp);
extern ESL_SQ* get_target_range_sequence(const TARGET_RANGE *target_range, const ESL_SQFILE *seq_file);
extern int release_hits_from_target_range(const TARGET_RANGE *target_range, int *hit_processed, int *num_hits_proccesed, int range_bound_min, int range_bound_max);


/* Initial Splice Graph */
extern int fill_graph_with_nodes(SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm);
extern SPLICE_EDGE* connect_nodes_with_edges(P7_HIT *upstream_hit, P7_HIT *downstream_hit, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, ESL_GENCODE *gcode, ESL_SQ *target_seq, int revcomp);

extern void get_overlap_nuc_coords (SPLICE_EDGE *edge, P7_DOMAIN *upstream, P7_DOMAIN *downstream, ESL_SQ *target_seq, int revcomp);
extern int find_optimal_splice_site (SPLICE_EDGE *edge, P7_DOMAIN *upstream, P7_DOMAIN *downstream, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, ESL_GENCODE *gcode, ESL_SQ *target_seq);
extern float ali_score_at_postion (P7_PROFILE *gm, int amino, int model_pos, int trans_pos, int prev_state, int curr_state);
extern int select_splice_option (SPLICE_EDGE *edge, P7_PROFILE *gm, P7_PROFILE *sub_model, ESL_GENCODE *gcode, ESL_SQ *target_seq, float signal_score, int up_nuc_pos, int down_nuc_pos);
extern int add_edge_to_graph(SPLICE_GRAPH *graph, SPLICE_EDGE *edge);


/* Missing Exons */
extern int fill_holes_in_graph(SPLICE_GRAPH *graph, TARGET_RANGE *target_range, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, ESL_SQ *target_seq, ESL_GENCODE *gcode);
extern SPLICE_GAP* find_the_gap (SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm, ESL_SQ *target_seq, int orig_N, int upstream_node, int downstream_node);
extern SPLICE_GAP* terminal_gap (SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm, ESL_SQ *target_seq, int orig_N, int final_node, int path_node, int look_upstream);
extern P7_HIT** align_the_gap(SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, ESL_SQ *target_seq, ESL_GENCODE *gcode, SPLICE_GAP *gap, int *num_hits);
extern P7_HMM* extract_sub_hmm (P7_HMM *hmm, int start, int end);
extern ESL_DSQ* extract_sub_seq(ESL_SQ *target_seq, int start, int end, int revcomp);
extern int add_missing_node_to_graph(SPLICE_GRAPH *graph, P7_TOPHITS * th, P7_HIT *hit, int M);
extern int add_missed_hit_to_target_range(TARGET_RANGE *target_range, P7_HIT *hit, int *duplicate);
extern int bridge_the_gap(SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, ESL_GENCODE *gcode, ESL_SQ *target_seq, int prev_N, int orig_N);


/* Splice Path */
extern int check_for_loops (SPLICE_GRAPH *graph, P7_TOPHITS *th);
extern SPLICE_PATH* evaluate_paths (SPLICE_GRAPH *graph, P7_TOPHITS *th, ESL_SQ *target_seq, int orig_N);
extern int longest_path_upstream (SPLICE_GRAPH *graph);
extern int topological_sort_upstream(SPLICE_GRAPH *graph, int *visited, int *stack, int *stack_size, int node);
extern int split_hits_in_path (SPLICE_GRAPH *graph, SPLICE_PATH *path, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, ESL_GENCODE *gcode, ESL_SQ *target_seq, int orig_N);

/* Spliced Hit Processing */
extern int splice_path (SPLICE_GRAPH *graph, SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_TOPHITS *orig_tophits, P7_OPROFILE *om, P7_SCOREDATA *scoredata, ESL_SQ *target_seq, ESL_GENCODE *gcode, int64_t db_nuc_cnt, int orig_N, int *success);
extern int splice_path_frameshift (SPLICE_GRAPH *graph, SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_FS_PROFILE *gm_fs, P7_OPROFILE *om, P7_SCOREDATA *scoredata, ESL_SQ *target_seq, ESL_GENCODE *gcode, int64_t db_nuc_cnt, int orig_N, int* success);
extern int align_spliced_path (SPLICE_PIPELINE *pli, P7_OPROFILE *om, P7_SCOREDATA *scoredata, ESL_SQ *target_seq, ESL_GENCODE *gcode);
extern int align_spliced_path_frameshift (SPLICE_PIPELINE *pli, P7_FS_PROFILE *gm_fs, P7_OPROFILE *om, P7_SCOREDATA *scoredata, ESL_SQ *target_seq, ESL_GENCODE *gcode);
extern int compute_ali_scores(P7_DOMAIN *dom, P7_TRACE *tr, ESL_DSQ *amino_dsq, const P7_SCOREDATA *data, int K);
extern int compute_ali_scores_fs(P7_DOMAIN *dom, P7_TRACE *tr, ESL_DSQ *nuc_dsq, const P7_SCOREDATA *data, P7_FS_PROFILE *gm_fs, const ESL_ALPHABET *abc);

/* Debug Dumps */
extern void target_range_dump(FILE *fp, TARGET_RANGE *target_range, int print_hits);
extern void graph_dump(FILE *fp, SPLICE_GRAPH *graph, ESL_SQ *target_seq, int print_edges);
extern void path_dump(FILE *fp, SPLICE_PATH *path, ESL_SQ *target_seq);


