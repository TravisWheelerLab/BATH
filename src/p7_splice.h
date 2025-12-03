/* Structs and MACROs for bilding a splice graph */

#include "p7_config.h"

#include <string.h>
#ifdef HMMER_THREADS
#include <pthread.h>
#endif /*HMMER_THREADS*/

#include "easel.h"
#include "hmmer.h"

typedef struct _splice_edge {

  int frameshift;
  int bypass_checked;
  int jump_edge;

  int upstream_node_id;
  int downstream_node_id;

  int upstream_amino_end;      
  int downstream_amino_start; 

  int upstream_nuc_end;    
  int downstream_nuc_start; 

  float splice_score;
  float signal_score;


} SPLICE_EDGE;


typedef struct _splice_graph {

  /* Graph size intfo */
  int nalloc;
  int num_nodes;
  int tot_edges;
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
  int         *num_edges;
  int         *edge_mem;

} SPLICE_GRAPH;




typedef struct _splice_path {

  int revcomp;
  int frameshift;

  int alloc_len;
  int path_len;

  int *node_id;
  int *path_id;

  int *ihmm;
  int *jhmm;
  
  int64_t *iali;
  int64_t *jali;

} SPLICE_PATH;



typedef struct _splice_site_idx
{
  int **index;
  int *index_mem;

  float **score;
  float *score_mem; 

  int **lookback;
  int *lookback_mem;

  int   *parser_index;
  float *parser_scores;

  int alloc_M;
  int alloc_L;
  int alloc_Lx;

} SPLICE_SITE_IDX;



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
  double   S1;
  double   S2;
  double   incE;
  double   incT;

  float   *signal_scores;

  int64_t *orig_nuc_idx;

  ESL_SQ  *nuc_sq;
  ESL_SQ  *amino_sq; 

  P7_OMX  *fwd;
  P7_OMX  *bwd;  

  P7_GMX  *gfwd;
  P7_GMX  *gbwd;

  P7_GMX  *vit;
  P7_GMX  *spx;

  P7_BG   *bg;

  P7_HIT  *hit;

  SPLICE_SITE_IDX *sig_idx;  
  

} SPLICE_PIPELINE;


typedef struct _splice_info
{
  SPLICE_GRAPH      **graphs;     /* array of graphs to process               */
  SPLICE_GRAPH       *graph;      /* graph to splice                          */
  SPLICE_PIPELINE    *pli;        /* work pipeline                            */
  P7_HMM             *hmm;        /* query hmm                                */
  P7_OPROFILE        *om;         /* optimized query profile                  */
  P7_PROFILE         *gm;         /* non-optimized query profile              */
  P7_FS_PROFILE      *gm_fs;      /* non optimized frameshift query profile   */
  P7_TOPHITS         *tophits;    /* original tophits                         */
  P7_TOPHITS         *seeds;      /* seed hits from SSV                       */
  ESL_GENCODE        *gcode;      /* used for translation                     */
  ESL_SQFILE         *seq_file;   /* target sequence file                     */
  int64_t             db_nuc_cnt; /* sequence database size for e-values      */
  int                 num_graphs; /* total number of graphs                   */
  int                 thread_id;  /* ID of this thread                        */
  int                *graph_idx;  /* current graph index                      */
#ifdef HMMER_THREADS
  pthread_mutex_t    *mutex;      /* mutex for thread synchronization         */
#endif /*HMMER_THREADS*/
} SPLICE_WORKER_INFO;



/* MACROS for the SPLICE_SITE_IDX */


#define SIX0(k, signal)             (index[(k)][(signal)])                                                       //xxxxXXX
#define SIX1(k, signal, nuc1)       (index[(k)][SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + signal])           //XxxxxXX
#define SIX2(k, signal, nuc1, nuc2) (index[(k)][SPLICE_OFFSET_2 + (4*nuc1+nuc2) * p7S_SPLICE_SIGNALS + signal]) //XXxxxxX
 
#define SSX0(k, signal)             (score[(k)][signal])                                                        //xxxxXXX
#define SSX1(k, signal, nuc1)       (score[(k)][SPLICE_OFFSET_1 + nuc1 * p7S_SPLICE_SIGNALS + signal])          //XxxxxXX
#define SSX2(k, signal, nuc1, nuc2) (score[(k)][SPLICE_OFFSET_2 + (4*nuc1+nuc2) * p7S_SPLICE_SIGNALS + signal]) //XXxxxxX

#define SIGNAL(nuc1, nuc2)           (4*nuc1+nuc2)


#define DONOR_GT                  11
#define DONOR_GC                  9
#define DONOR_AT                  3
#define ACCEPT_AG                 2
#define ACCEPT_AC                 1

#define SIGNAL_MEM_SIZE           63       /* total storage = M * SIGNAL_MEM_SIZE */ 
#define SPLICE_OFFSET_1           3        /* start of XxxxxXX codons  */
#define SPLICE_OFFSET_2           15       /* start of XXxxxxX codons  */

#define EDGE_ALLOC                10       /*minimum alloc space for edges for  each node */
#define MAX_INTRON_LENG           200000   /*maximum intron length */
#define MAX_INTRON_EXT            10000    /*maximum extension distance */
#define MIN_INTRON_LENG           13       /*minimum intor length */
#define MAX_INTRON_INCL           1500     /*maximum length on intron to be included in spliced Viterbi search */
#define MAX_AMINO_GAP             100      /*max amino gap the spliced viterbi algoritm will try to bridge without searching the full intron */
#define MAX_SP_AMINO_GAP          10       /*maximum amino gap for spliced edges */
#define MAX_USP_AMINO_GAP         25       /*maximum amino gap fpr unspliced edges */ 
#define MAX_EXT_AMINO_GAP         25       /*maximum amino gap fpr extention edges */
#define MIN_AMINO_OVERLAP         10       /*amino acid splice site window */   
#define ALIGNMENT_EXT             10       /*amino extetion at start and end of path seq for final alignment */

/* SPLICE_SITE_IDX modes */
#define PARSER_MODE    0
#define ALIGNMENT_MODE 1

enum p7_confirm_e {
  p7_CONFIRM_START = 0,
  p7_CONFIRM_END   = 1,
  p7_CONFIRM_NONE  = 2
};

/* Indices of p7_splice_SignalScores */
enum p7s_splice_signals_e {
  p7S_GTAG  = 0,
  p7S_GCAG  = 1,
  p7S_ATAC  = 2,
};
#define p7S_SPLICE_SIGNALS 3

enum p7s_splice_codons_e {
  xxxxXXX  = 0,
  XxxxxXX  = 1,
  XXxxxxX  = 2,
};
#define p7S_SPLICE_CODONS 3

/* p7_spliceedge.c */
extern SPLICE_EDGE* p7_spliceedge_Create(void);
extern int p7_spliceedge_AliScoreEdge(SPLICE_EDGE *edge, const P7_PROFILE *gm, const P7_DOMAIN *upstream_dom, const P7_DOMAIN *downstream_dom);

/* p7_splicegraph.c */
extern SPLICE_GRAPH* p7_splicegraph_Create(void);
extern int p7_splicegraph_CreateNodes(SPLICE_GRAPH *graph, int num_nodes);
extern int p7_splicegraph_Grow(SPLICE_GRAPH *graph);
extern void p7_splicegraph_Destroy(SPLICE_GRAPH *graph);
extern int p7_splicegraph_AddNode(SPLICE_GRAPH *graph, P7_HIT *hit);
extern SPLICE_EDGE* p7_splicegraph_AddEdge(SPLICE_GRAPH *graph, int up_node, int down_node);
extern int p7_splicegraph_EdgeExists(SPLICE_GRAPH* graph, int up_node, int down_node);
extern SPLICE_EDGE* p7_splicegraph_GetEdge(SPLICE_GRAPH* graph, int up_node, int down_node);
extern int p7_splicegraph_PathExists (SPLICE_GRAPH *graph, int upstream_node, int downstream_node);
extern int p7_splicegraph_RemoveDuplicates(SPLICE_GRAPH *graph);
extern void p7_splicegraph_DumpHits(FILE *fp, SPLICE_GRAPH *graph);
extern void p7_splicegraph_DumpEdges(FILE *fp, SPLICE_GRAPH *graph);
extern void p7_splicegraph_DumpGraph(FILE *fp, SPLICE_GRAPH *graph);

/* p7_splicehits.c */
extern SPLICE_SAVED_HITS* p7_splicehits_CreateSavedHits(void);
extern int p7_splicehits_GrowSavedHits(SPLICE_SAVED_HITS *saved_hits);
extern void p7_splicehits_DestroySavedHits(SPLICE_SAVED_HITS *saved_hits);
extern int p7_splicehits_CreateNext(SPLICE_SAVED_HITS *saved_hits, SPLICE_HIT_INFO **ret_info);
extern int p7_splicehits_SortSavedHits(SPLICE_SAVED_HITS *sh);
extern int p7_splicehits_MergeSavedHits(SPLICE_SAVED_HITS *sh1, SPLICE_SAVED_HITS *sh2);
extern int p7_splicehits_AssignNodes(SPLICE_GRAPH *graph, SPLICE_SAVED_HITS *sh, int first, int last);
extern int p7_splicehits_RemoveDuplicates(SPLICE_SAVED_HITS *sh, P7_TOPHITS *th, double F3);
extern P7_TOPHITS* p7_splicehits_GetSeedHits(SPLICE_SAVED_HITS *sh, const P7_TOPHITS *th, P7_HMM *hmm, P7_FS_PROFILE *gm_fs, ESL_SQFILE *seq_file, ESL_GENCODE *gcode, double F3); 
extern void p7_splicehits_Dump(FILE *fp, SPLICE_SAVED_HITS *sh);
 
/* p7_splicepath.c */
extern SPLICE_PATH* p7_splicepath_Create(int path_len);
extern int p7_splicepath_Grow(SPLICE_PATH *path);
extern SPLICE_PATH* p7_splicepath_Clone(SPLICE_PATH *path);
extern int p7_splicepath_Insert(SPLICE_PATH *path, int step);
extern int p7_splicepath_Remove(SPLICE_PATH *path, int step);
extern void p7_splicepath_Destroy(SPLICE_PATH *path);
extern SPLICE_PATH* p7_splicepath_GetBestPath_Unspliced(SPLICE_GRAPH *graph);
extern SPLICE_PATH* p7_splicepath_GetBestPath_Extension(SPLICE_GRAPH *orig_graph, SPLICE_GRAPH *extend_graph);
extern void p7_splicepath_Dump(FILE *fp, SPLICE_PATH *path);
extern void p7_splicepath_DumpScores(FILE *fp, SPLICE_PATH *path, SPLICE_GRAPH *graph);

/* p7_splicepipeline.c */
extern SPLICE_PIPELINE* p7_splicepipeline_Create(const ESL_GETOPTS *go, int M_hint, int L_hint);
extern void p7_splicepipeline_Reuse(SPLICE_PIPELINE *pli);
extern void p7_splicepipeline_Destroy(SPLICE_PIPELINE *pli); 
extern SPLICE_SITE_IDX* p7_splicepipline_CreateIndex(int M_hint, int L_hint, int Lx_hint);
extern int p7_splicepipline_GrowIndex(SPLICE_SITE_IDX *signal_sites, int M, int L, int Lx);
extern void p7_splicepipeline_DestroyIndex(SPLICE_SITE_IDX *signal_sites);


/* p7_spliceviterbi.c */
extern int p7_spliceviterbi_parser_semiglobal(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx);
extern int p7_spliceviterbi_rightparser_semiglobal(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx);
extern int p7_spliceviterbi_leftparser_semiglobal(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx);
extern int p7_spliceviterbi_translated_semiglobal(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx);
extern int p7_splicevitebi_translated_semiglobal_trace(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, int L, const ESL_GENCODE *gcode, const P7_FS_PROFILE *sub_gm, const P7_GMX *gx, P7_TRACE *tr); 
extern int p7_splicevitebi_exon_definition(SPLICE_PIPELINE *pli, SPLICE_PATH *path, P7_GMX *gx, P7_FS_PROFILE *sub_gm, ESL_SQ *sub_seq, int **exon_idx, int *idx_size);
extern int p7_spliceviterbi_lefttranslated_semiglobal(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx);
extern int p7_splicevitebi_lefttranslated_semiglobal_trace(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, int L, const ESL_GENCODE *gcode, const P7_FS_PROFILE *sub_gm, const P7_GMX *gx, P7_TRACE *tr);

/* p7_splice.c */
extern int p7_splice_SpliceGraph(SPLICE_WORKER_INFO *info);
extern int p7_splice_AddOriginals(SPLICE_WORKER_INFO *info, SPLICE_GRAPH *graph, const P7_TOPHITS *tophits);
extern int p7_splice_AddSeeds(SPLICE_GRAPH *graph, const P7_TOPHITS *seed_hits);
extern int p7_splice_ExtendPath(P7_TOPHITS *seed_hits, P7_PROFILE *gm, SPLICE_PATH *path, SPLICE_GRAPH *graph);
extern int p7_splice_CreateUnsplicedEdges(SPLICE_GRAPH *graph, P7_PROFILE *gm);
extern int p7_splice_CreateExtensionEdges(SPLICE_GRAPH *graph, P7_PROFILE *gm);
extern SPLICE_PATH* p7_splice_AlignExons(SPLICE_PIPELINE *pli, P7_HMM *sub_hmm, const P7_FS_PROFILE *gm_fs, P7_BG *bg, ESL_SQ *ali_seq, const ESL_GENCODE *gcode, int removed_start, int removed_end);
extern int p7_splice_RemoveHits(SPLICE_GRAPH *graph, int range_bound_min, int range_bound_max);
extern int p7_splice_EnforceRangeBounds(SPLICE_GRAPH *graph, int64_t bound_min, int64_t bound_max);
extern int p7_splice_HitUpstream(P7_DOMAIN *upstream, P7_DOMAIN *downstream, int revcomp);
extern SPLICE_PATH* p7_splice_FindExons(SPLICE_WORKER_INFO *info, SPLICE_PATH *path, ESL_SQ *path_seq);
extern int p7_splice_AlignPath(SPLICE_GRAPH *graph, SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_TOPHITS *tophits, P7_OPROFILE *om, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *path_seq, int64_t db_nuc_cnt, float fs_prob, SPLICE_WORKER_INFO *info, int *success);
//extern int p7_splice_AlignFrameshiftPath(SPLICE_GRAPH *graph, SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_TOPHITS *tophits, P7_FS_PROFILE *gm_fs, ESL_GENCODE *gcode, ESL_SQ *path_seq, int64_t db_nuc_cnt, SPLICE_WORKER_INFO *info);
extern ESL_SQ* p7_splice_GetSubSequence(const ESL_SQFILE *seq_file, char* seqname, int64_t seq_min, int64_t seq_max, int revcomp, SPLICE_WORKER_INFO *info);
extern ESL_SQ* p7_splice_GetSplicedSequence(SPLICE_PATH *path, ESL_SQ *path_seq, int up_node, int down_node, int intron_include, int **remove_idx, int64_t **nuc_index);
extern ESL_SQ* p7_splice_GetSplicedSequence2(SPLICE_PATH *path, ESL_SQ *path_seq, int *exon_idx, int exon_cnt, int *remove_idx, int *remove_cnt, int64_t **nuc_index); 
extern P7_HMM* p7_splice_GetSubHMM (const P7_HMM *hmm, int start, int end);



