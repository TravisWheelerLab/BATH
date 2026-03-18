/* Structs and MACROs for bilding a splice graph */

#include "p7_config.h"

#include <string.h>
#ifdef HMMER_THREADS
#include <pthread.h>
#endif /*HMMER_THREADS*/

#include "easel.h"
#include "hmmer.h"


typedef struct _splice_edge {

  int jump_edge;

  int upstream_node_id;
  int downstream_node_id;

  /* Spliced coordinates */
  int upstream_amino_end;      
  int downstream_amino_start; 

  int upstream_nuc_end;    
  int downstream_nuc_start; 

  /* unspliced starts of upstream node */
  int i_start;  //unspliced start of upstream node
  int k_start;  

  /* unspliced ends of upstream node */
  int i_end;
  int k_end;

  /* unspliced starts of downstream node */
  int next_i_start;
  int next_k_start;

  float edge_score;


} SPLICE_EDGE;


typedef struct _splice_graph {

  /* Graph size intfo */
  int nalloc;
  int num_nodes;
  int anchor_N;
 
  /* Target sequence info */
  int          revcomp;
  int64_t      seqidx;      
  char        *seqname;    

  /* Node info */
  int         *node_in_graph;   /* Is the hit part of the current graph */
  int         *tmp_node;        /* New nodes found durring splicing */
  int         *orig_hit_idx;    /* index of hits in original P7_TOPHITS  */

  /*Edge info */  
  int   *best_out_edge;
  
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
  int *extension;

  int *ihmm;
  int *jhmm;
  
  int64_t *iali;
  int64_t *jali;

  float *aliscore;

} SPLICE_PATH;


typedef struct _splice_scores
{
  int allocM;

  float **score;
  float *score_mem;

  float *signal_scores;

  float *acceptor_AG;
  float *acceptor_AC;

  float *donor_GT;
  float *donor_GC;
  float *donor_AT;

} SPLICE_SCORES;


typedef struct _splice_pipeline 
{
  
  int by_E;
  int inc_by_E;
  int do_null2;
  int do_biasfilter;
  int show_cigar; 

  int min_intron;
  int max_intron;
  int max_extend;

  double E;
  double T;
  double Z;
  double F1;  
  double F2;
  double F3;
  double incE;
  double incT;

  int64_t *orig_nuc_idx;

  ESL_SQ  *nuc_sq;
  ESL_SQ  *amino_sq; 

  P7_OMX  *fwd;
  P7_OMX  *bwd;  
  P7_OMX  *pp;

  P7_GMX  *gfwd;
  P7_GMX  *gbwd;
  P7_GMX  *gpp;

  P7_GMX  *vit;

  P7_BG   *bg;

  P7_HIT  *hit;

  SPLICE_SCORES *splice_scores;

} SPLICE_PIPELINE;

typedef struct _splice_bounds
{

  int     N;
  int     allocN;

  int     *bound_hmm_mins;
  int     *bound_hmm_maxs;

  int64_t *bound_seq_mins;
  int64_t *bound_seq_maxs;

} SPLICE_BOUNDS;

typedef struct _splice_info
{
  SPLICE_GRAPH      **graphs;     /* array of graphs to process               */
  SPLICE_GRAPH       *graph;      /* graph to splice                          */
  SPLICE_PIPELINE    *pli;        /* work pipeline                            */
  P7_OPROFILE        *om;         /* optimized query profile                  */
  P7_PROFILE         *gm;         /* non-optimized query profile              */
  P7_FS_PROFILE      *gm_tr;      /* non optimized translated query profile   */
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



/* MACROS for the P state score storage */
#define SSX0(k, signal)             (score[(signal)][(k)])                                                            //xxxxXXX
#define SSX1(k, signal, nuc1)       (score[SPLICE_OFFSET_1 + (nuc1) * p7S_SPLICE_SIGNALS + (signal)][(k)])            //XxxxxXX
#define SSX2(k, signal, nuc1, nuc2) (score[SPLICE_OFFSET_2 + (4*(nuc1)+(nuc2)) * p7S_SPLICE_SIGNALS + (signal)][(k)]) //XXxxxxX

#define SIGNAL(nuc1, nuc2)           (4*(nuc1)+(nuc2))

#define DONOR_GT                  11
#define DONOR_GC                  9
#define DONOR_AT                  3
#define ACCEPT_AG                 2
#define ACCEPT_AC                 1

#define SIGNAL_MEM_SIZE           63       /* total storage = M * SIGNAL_MEM_SIZE */ 
#define SPLICE_OFFSET_1           3        /* start of XxxxxXX codons  */
#define SPLICE_OFFSET_2           15       /* start of XXxxxxX codons  */

#define EDGE_ALLOC                10       /*minimum alloc space for edges for  each node           */
#define MAX_INTRON_EXT            10000    /*maximum extension distance                             */
#define MAX_AMINO_GAP             1500     /*maximum amino gap for an edge                          */
#define ALIGNMENT_EXT             30       /*extention at start and end of final alignment sequence */

/* Indices of p7_splice_SignalScores */
enum p7s_splice_signals_e {
  p7S_GTAG  = 0,
  p7S_GCAG  = 1,
  p7S_ATAC  = 2,
};
#define p7S_SPLICE_SIGNALS 3

/* p7_splicebounds.c */
extern SPLICE_BOUNDS* p7_splicebounds_Create(int allocN);
extern int p7_splicebounds_GorwTo(SPLICE_BOUNDS *bounds, int allocN);
extern void p7_splicebounds_Destroy(SPLICE_BOUNDS *bounds);
extern int p7_splicebounds_Add(SPLICE_BOUNDS *bounds, int64_t seq_min, int64_t seq_max, int hmm_min, int hmm_max);

/* p7_splicegraph.c */
extern SPLICE_GRAPH* p7_splicegraph_Create(void);
extern int p7_splicegraph_CreateNodes(SPLICE_GRAPH *graph, int num_nodes);
extern int p7_splicegraph_Grow(SPLICE_GRAPH *graph);
extern void p7_splicegraph_Destroy(SPLICE_GRAPH *graph);
extern int p7_splicegraph_AddNode(SPLICE_GRAPH *graph, P7_HIT *hit);
extern SPLICE_EDGE* p7_splicegraph_AddEdge(SPLICE_GRAPH *graph, int up_node, int down_node);
extern int p7_splicegraph_EdgeExists(SPLICE_GRAPH* graph, int up_node, int down_node);
extern SPLICE_EDGE* p7_splicegraph_GetEdge(SPLICE_GRAPH* graph, int up_node, int down_node);
extern int p7_splicegraph_AliScoreEdge(SPLICE_EDGE *edge, const P7_DOMAIN *upstream_dom, const P7_DOMAIN *downstream_dom);
extern int p7_splicegraph_PathExists (SPLICE_GRAPH *graph, int upstream_node, int downstream_node);
extern int p7_splicegraph_NodeOverlap(SPLICE_GRAPH *graph, int node_id, SPLICE_PATH *path, int step_id);
extern void p7_splicegraph_DumpHits(FILE *fp, SPLICE_GRAPH *graph);
extern void p7_splicegraph_DumpEdges(FILE *fp, SPLICE_GRAPH *graph);
extern void p7_splicegraph_DumpGraph(FILE *fp, SPLICE_GRAPH *graph);

/* p7_splicepath.c */
extern SPLICE_PATH* p7_splicepath_Create(int path_len);
extern int p7_splicepath_Grow(SPLICE_PATH *path);
extern SPLICE_PATH* p7_splicepath_Clone(SPLICE_PATH *path);
extern int p7_splicepath_Insert(SPLICE_PATH *path, int step);
extern int p7_splicepath_Remove(SPLICE_PATH *path, int step);
extern void p7_splicepath_Destroy(SPLICE_PATH *path);
extern SPLICE_PATH* p7_splicepath_GetBestPath(SPLICE_GRAPH *graph, int extend_up, int extend_down);
extern void p7_splicepath_Dump(FILE *fp, SPLICE_PATH *path);
extern void p7_splicepath_DumpScores(FILE *fp, SPLICE_PATH *path, SPLICE_GRAPH *graph);

/* p7_splicepipeline.c */
extern SPLICE_PIPELINE* p7_splicepipeline_Create(const ESL_GETOPTS *go, int M_hint, int L_hint);
extern void p7_splicepipeline_Reuse(SPLICE_PIPELINE *pli);
extern void p7_splicepipeline_Destroy(SPLICE_PIPELINE *pli); 

/* p7_splicescores.c */
SPLICE_SCORES* p7_splicescores_Create(int M_hint);
extern int p7_splicescores_GrowTo(SPLICE_SCORES* splice_scores, int M);
extern void p7_splicescores_Destroy(SPLICE_SCORES* splice_scores);

/* p7_spliceviterbi.c */
extern int p7_spliceviterbi_TranslatedGlobal(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, P7_GMX *gx, int i_start, int i_end, int k_start, int k_end);
extern int p7_spliceviterbi_TranslatedSemiGlobalExtendUp(SPLICE_PIPELINE *pli, const ESL_DSQ *path_dsq, const P7_FS_PROFILE *gm_tr, P7_GMX *gx, int i_start, int i_end, int k_start, int k_end);
extern int p7_spliceviterbi_TranslatedSemiGlobalExtendDown(SPLICE_PIPELINE *pli, const ESL_DSQ *path_dsq, const P7_FS_PROFILE *gm_tr, P7_GMX *gx, int i_start, int i_end, int k_start, int k_end);
extern int p7_spliceviterbi_TranslatedTrace(SPLICE_PIPELINE *pli, const ESL_DSQ *sub_dsq, const P7_FS_PROFILE *gm_tr, const P7_GMX *gx, P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end);


/* p7_splice.c */
extern int p7_splice_SpliceGraph(SPLICE_WORKER_INFO *info);
extern int p7_splice_AddAnchors(SPLICE_WORKER_INFO *info, SPLICE_GRAPH *graph, const P7_TOPHITS *tophits);
extern int p7_splice_AddSeeds(SPLICE_WORKER_INFO *info, SPLICE_GRAPH *graph, const P7_TOPHITS *seed_hits);
extern int p7_splice_ExtendPath(SPLICE_PIPELINE *pli, P7_TOPHITS *seed_hits, SPLICE_PATH *path, SPLICE_PATH *spliced_path, SPLICE_GRAPH *graph, SPLICE_BOUNDS *bounds);
extern int p7_splice_CreateUnsplicedEdges(SPLICE_PIPELINE *pli, SPLICE_GRAPH *graph, P7_FS_PROFILE *gm_tr);
extern int p7_splice_CreateExtensionEdges(SPLICE_PIPELINE *pli, SPLICE_GRAPH *orig_graph, SPLICE_GRAPH *extension_graph);
extern SPLICE_PATH* p7_splice_AlignExons(SPLICE_WORKER_INFO *info, SPLICE_PATH *orig_path, ESL_SQ *path_seq, int down, int i_start, int i_end, int k_start, int k_end, int *next_i_start, int *next_k_start);
extern SPLICE_PATH* p7_splice_AlignExtendDown(SPLICE_WORKER_INFO *info, SPLICE_PATH *spliced_path, ESL_SQ *path_seq, int s_end, int i_start, int i_end, int k_start, int k_end, int *next_i_end, int *next_k_end);
extern SPLICE_PATH* p7_splice_AlignExtendUp(SPLICE_WORKER_INFO *info, SPLICE_PATH *spliced_path, ESL_SQ *path_seq, int s_start, int i_start, int i_end, int k_start, int k_end);
extern SPLICE_PATH* p7_splice_AlignSingle(SPLICE_WORKER_INFO *info, SPLICE_PATH *spliced_path, ESL_SQ *path_seq, int i_start, int i_end, int k_start, int k_end);
extern int p7_splice_EnforceBounds(SPLICE_GRAPH *graph, int64_t bound_min, int64_t bound_max);
extern int p7_splice_HitUpstream(P7_DOMAIN *upstream, P7_DOMAIN *downstream, int revcomp);
extern int p7_splice_HitBetween(P7_DOMAIN *up, P7_DOMAIN *mid, P7_DOMAIN *down, int revcomp);
extern SPLICE_PATH* p7_splice_SpliceExons(SPLICE_WORKER_INFO *info, SPLICE_PATH *path, ESL_SQ *path_seq);
extern int p7_splice_SpliceExtensions(SPLICE_WORKER_INFO *info, SPLICE_PATH *path, ESL_SQ *path_seq);
extern int p7_splice_SpliceSingle(SPLICE_WORKER_INFO *info, SPLICE_PATH *path, ESL_SQ *path_seq);
extern int p7_splice_AlignSplicedPath(SPLICE_WORKER_INFO *info, SPLICE_PATH *orig_path, SPLICE_PATH *spliced_path, ESL_SQ *path_seq, int *success);
extern int p7_splice_CreateSplicedSequnce(SPLICE_WORKER_INFO *info, SPLICE_PATH *spliced_path, ESL_SQ *path_seq);
extern int p7_splice_AlignSplicedSequence(SPLICE_WORKER_INFO *info, SPLICE_PATH *spliced_path, ESL_SQ *path_seq);
extern int p7_splice_FixDecodingErrors(SPLICE_GRAPH *graph, SPLICE_PATH *spliced_path, P7_ALIDISPLAY *ad, ESL_SQ *path_seq);
extern int p7_splice_ScoreExons(SPLICE_PIPELINE *pli, P7_TRACE *tr, P7_ALIDISPLAY *ad, P7_OPROFILE *om, int do_pp);
extern ESL_SQ* p7_splice_GetSubSequence(const ESL_SQFILE *seq_file, char* seqname, int64_t seq_min, int64_t seq_max, int revcomp, SPLICE_WORKER_INFO *info);



