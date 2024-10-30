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
 
  int   has_full_path;
  int   best_path_length;
  int   best_path_start;
  int   best_path_end;

  int   *out_edge_cnt;
  int   *in_edge_cnt;
  int   *orig_out_edge;
  int   *orig_in_edge;
  int   *best_out_edge;
  int   *best_in_edge;

  int   num_n_term;
  int   num_c_term;
  int   *is_n_terminal;
  int   *is_c_terminal;
 
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

  int *node_id;

  int *missing;

  int *upstream_spliced_amino_end;
  int *downstream_spliced_amino_start;

  int *upstream_spliced_nuc_end;
  int *downstream_spliced_nuc_start;
 
  int *checked_upstream;
  int *checked_downstream;
 
  float *hit_scores;
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

#define MAX_TARGET_RANGE_EXT      3000000    //3x10^6 
#define TERM_RANGE_EXT            100000     //1x10^5
#define MAX_INTRON_LEN            10000     //1x10^4 
#define MIN_INTRON_LEN            10
#define MAX_AMINO_EXT             6

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
extern void target_range_destroy(TARGET_RANGE *target_range);
extern SPLICE_EDGE* splice_edge_create(void);
extern SPLICE_GRAPH* splice_graph_create(void);
extern void splice_graph_destroy (SPLICE_GRAPH* graph);
extern int splice_graph_create_nodes(SPLICE_GRAPH *graph, int num_nodes);
extern SPLICE_PATH* splice_path_create(int path_len);
extern void splice_path_destroy(SPLICE_PATH *path);
extern SPLICE_PIPELINE* splice_pipeline_create(const ESL_GETOPTS *go, int M_hint, int L_hint);
extern void splice_pipeline_destroy(SPLICE_PIPELINE *pli);


/* Target Range  */
extern TARGET_RANGE* build_target_range(const P7_TOPHITS *tophits, int *hits_processed, int *num_hits_processed, int64_t *range_bound_mins, int64_t *range_bound_maxs, int range_num, int64_t seqidx, int revcomp);
extern ESL_SQ* get_target_range_sequence(const TARGET_RANGE *target_range, const ESL_SQFILE *seq_file);
extern int release_hits_from_target_range(const TARGET_RANGE *target_range, int *hit_processed, int *num_hits_proccesed, int range_bound_min, int range_bound_max);


/* Initial Splice Graph */
extern int fill_graph_with_nodes(SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm);
extern SPLICE_EDGE* connect_nodes_with_edges(P7_HIT *upstream_hit, P7_HIT *downstream_hit, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *target_seq, int revcomp);
extern void get_overlap_nuc_coords (SPLICE_EDGE *edge, P7_DOMAIN *upstream, P7_DOMAIN *downstream, ESL_SQ *target_seq, int revcomp);
extern int find_optimal_splice_site (SPLICE_EDGE *edge, P7_DOMAIN *upstream, P7_DOMAIN *downstream, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *target_seq, int revcomp);
extern float ali_score_at_postion (P7_PROFILE *gm, int amino, int model_pos, int trans_pos, int prev_state, int curr_state);
extern int select_splice_option (SPLICE_EDGE *edge, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *target_seq, float *splice_scores, float up_score, float down_score, int model_pos, int up_nuc_pos, int down_nuc_pos, int up_state, int down_state);
extern int add_edge_to_graph(SPLICE_GRAPH *graph, SPLICE_EDGE *edge);


/* Missing Exons */
extern int fill_holes_in_graph(SPLICE_GRAPH *graph, TARGET_RANGE *target_range, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, ESL_SQ *target_seq, ESL_GENCODE *gcode);
extern int extend_path(SPLICE_GRAPH *graph, SPLICE_PATH *path, TARGET_RANGE *target_range, P7_HMM *hmm, P7_PROFILE *gm, P7_BG *bg, ESL_GENCODE *gcode, ESL_SQ *target_seq);
extern SPLICE_GAP* find_the_gap (SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm, ESL_SQ *target_seq, int orig_N, int upstream_node, int downstream_node);
extern SPLICE_GAP* terminal_gap (SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm, ESL_SQ *target_seq, int orig_N, int final_node, int path_node, int look_upstream);
extern P7_HIT** align_the_gap(SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm, P7_HMM *hmm, P7_BG *bg, ESL_SQ *target_seq, ESL_GENCODE *gcode, SPLICE_GAP *gap, int *num_hits);
extern P7_HMM* extract_sub_hmm (P7_HMM *hmm, int start, int end);
extern ESL_DSQ* extract_sub_seq(ESL_SQ *target_seq, int start, int end, int revcomp);
extern int add_missing_node_to_graph(SPLICE_GRAPH *graph, P7_TOPHITS * th, P7_HIT *hit, int M);
extern int add_missed_hit_to_target_range(TARGET_RANGE *target_range, P7_HIT *hit);
extern int bridge_the_gap(SPLICE_GRAPH *graph, P7_TOPHITS *th, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *target_seq, int prev_N, int orig_N);


/* Splice Path */
extern int check_for_loops (SPLICE_GRAPH *graph, P7_TOPHITS *th);
extern SPLICE_PATH* evaluate_paths (SPLICE_GRAPH *graph, P7_TOPHITS *th, ESL_SQ *target_seq, int orig_N);
extern int longest_path_upstream (SPLICE_GRAPH *graph);
extern int topological_sort_upstream(SPLICE_GRAPH *graph, int *visited, int *stack, int *stack_size, int node);


/* Spliced Hit Processing */
extern int splice_path (SPLICE_GRAPH *graph, SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_OPROFILE *om, ESL_SQ *target_seq, ESL_GENCODE *gcode, int64_t db_nuc_cnt, int orig_N, int *success);
extern int align_spliced_path (SPLICE_PIPELINE *pli, P7_OPROFILE *om, ESL_SQ *target_seq, ESL_GENCODE *gcode);


/* Debug Dumps */
extern void target_range_dump(FILE *fp, TARGET_RANGE *target_range, int print_hits);
extern void graph_dump(FILE *fp, SPLICE_GRAPH *graph, ESL_SQ *target_seq, int print_edges);
extern void path_dump(FILE *fp, SPLICE_PATH *path, ESL_SQ *target_seq);



/* BELLOW THIS LINE ARE FUNCTIONS FOR TESTING*/
static int
integer_textwidth2(long n)
{
  int w = (n < 0)? 1 : 0;
  while (n != 0) { n /= 10; w++; }
  return w;
}

int
diplay_frame2(int nuc_from, int nuc_to)
{
  int frame = 0;
  if (nuc_from < nuc_to)
  {
    frame = (nuc_to+1) % 3;
    if (frame == 0)
      frame = 3;
  }
  else
  {
    frame = (-1)*(nuc_to % 3);
    if (frame == 0)
      frame = -3;
  }
  return frame;
}

int
Splice_Test_Print(FILE *fp, P7_ALIDISPLAY *ad, int min_aliwidth, int linewidth, P7_PIPELINE *pli, int set_num)
{
  char *buf          = NULL;
  char *show_hmmname = NULL;
  char *show_seqname = NULL;
  int  *frameline    = NULL;
  int   namewidth, coordwidth;
  int   max_aliwidth, cur_aliwidth;
  int   frame;
  int   pos;
  int   status;
  int   ni, nk;
  int   x,y,z;
  long  i1,i2;
  long  c1;
  int   k1,k2;
  int   npos;
  int   i,j;
  int   show_accessions;
  int   exon_cnt;
  int   spliced_ali, is_splice_line;
  int   splice_pos, old_splice_pos;
  int   splice_opt, old_splice_opt;
  char  *exon_name;

  show_accessions = pli->show_accessions;


  /* implement the --acc option for preferring accessions over names in output  */
  show_hmmname = (show_accessions && ad->hmmacc[0] != '\0') ? ad->hmmacc : ad->hmmname;
  show_seqname = (show_accessions && ad->sqacc[0]  != '\0') ? ad->sqacc  : ad->sqname;

  /* dynamically size the output lines */
  namewidth  = ESL_MAX(strlen(show_hmmname), strlen(show_seqname));

  namewidth  = ESL_MAX(namewidth, 8);

  coordwidth = ESL_MAX(
           ESL_MAX(integer_textwidth2(ad->hmmfrom), integer_textwidth2(ad->hmmto)),
           ESL_MAX(integer_textwidth2(ad->sqfrom), integer_textwidth2(ad->sqto)));

  max_aliwidth   = (linewidth > 0) ? linewidth - namewidth - 2*coordwidth - 5 : ad->N;

  if (max_aliwidth < ad->N && max_aliwidth < min_aliwidth) max_aliwidth = min_aliwidth; /* at least, regardless of some silly linewidth setting */

  max_aliwidth -= 2; /* two extra spaces required for printing splice sites */
  max_aliwidth /= 5; /* divide by 5 if printing codons horizontally */

  ESL_ALLOC(buf, sizeof(char) * (max_aliwidth+1));
  buf[max_aliwidth] = 0;

  if(pli->show_frameline) ESL_ALLOC(frameline, sizeof(int) * (max_aliwidth+1));

   /* Break the alignment into multiple blocks of width max_aliwidth for printing */
  i1 = ad->sqfrom;
  if(ad->sqfrom < ad->sqto) i2 = i1 - 1;
  else                      i2 = i1 + 1;

  k1 = ad->hmmfrom;

  ESL_ALLOC(exon_name, sizeof(char)*10);
  if(ad->exon_cnt)  spliced_ali = TRUE;
  else              spliced_ali = FALSE;

  old_splice_opt = splice_opt = 0;
  old_splice_pos = splice_pos = 0;

  exon_cnt = 1;
  pos = 0;
  while(pos < ad->N)
  {

    if (pos > 0) { if (fprintf(fp, "\n") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed"); } /* blank line betweeen blocks */
    cur_aliwidth = max_aliwidth;
    is_splice_line = FALSE;
    if( spliced_ali ) {
      /*Determine if the line we are about to print conatins a splice boundry */
      splice_opt = 0;
      for (z = pos; z < pos + max_aliwidth + 1 && z < ad->N; z++) {
        if (ad->ntseq[z*5] == 0) break;
        /* '$' in the first nucleotide postion denotes a splice site */
        if      (ad->ntseq[z*5]   == '$') {
          is_splice_line = TRUE;
          /* '$' in the 2nd, 3rd, or 4th nucleotide positions
           * denote the exact location of the splice signals */
          if      (ad->ntseq[z*5+1] == '$') splice_opt = 1;
          else if (ad->ntseq[z*5+2] == '$') splice_opt = 2;
          else if (ad->ntseq[z*5+3] == '$') splice_opt = 3;
          x = z;
          break;
        }
      }

      /* If a splice boundry was detected reset the width of the
      * line we are printing to stop at the splice boundry */
      if(is_splice_line) {
        cur_aliwidth    = x-pos;
        splice_pos      = x*5;
      }
    }
    /* Special case were type one splice happens at end of line */
    if (spliced_ali &&  (!is_splice_line) && ad->ntseq[z*5] == '$' && ad->ntseq[z*5+1] == '$') {
         is_splice_line = TRUE;
         splice_opt = 4;
         splice_pos = z*5;
    }
    ni = nk = 0;
    for (z = pos; z < pos + cur_aliwidth && z < ad->N; z++)
    {
      if (ad->model[z] != '.' && ad->model[z] != ' ') nk++; /* k advances except on insert states */
      if (ad->aseq[z]  != '-')                        ni++; /* i advances except on delete states */
    }
    k2 = k1+nk-1;

    /* Optional CS Line */
    if (ad->csline != NULL)
    {
      if (fprintf(fp, "  %*s ", namewidth+coordwidth+1, " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      /* Print space at begining of line to make room for splice signal */
      if (fprintf(fp,  " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      for (i = 0; i < cur_aliwidth; i++)
      {
        if (ad->csline[pos+i] == 0) break;
        if (fprintf(fp,  "  %c  ", ad->csline[pos+i]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      }
      if (fprintf(fp,  " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }

    /* Optional RF Line */
    if (ad->rfline != NULL)
    {
      if (fprintf(fp, "  %*s ", namewidth+coordwidth+1, " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      if (fprintf(fp,  " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      for (i = 0; i < cur_aliwidth; i++)
      {
        if (ad->rfline[pos+i] == 0) break;
        if (fprintf(fp, "  %c  ",  ad->rfline[pos+i]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
     }
     if (fprintf(fp,  " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
     if (fprintf(fp, "RF\n") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }

 /*  Mandatory Model Line */
    if (fprintf(fp, "  %*s %*d ", namewidth,  show_hmmname, coordwidth, k1) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    if (fprintf(fp,  " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    for (i = 0; i < cur_aliwidth; i++)
    {
      if (ad->model[pos+i] == 0) break;
      if (fprintf(fp, "  %c  ", ad->model[pos+i]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }
    if (fprintf(fp,  " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    if (fprintf(fp, " %-*d\n", coordwidth, k2) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");

    /* Mandatory Match Line */
    if (fprintf(fp, "  %*s ", namewidth+coordwidth+1, " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    if (fprintf(fp,  " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    for (i = 0; i < cur_aliwidth; i++)
    {
      if (ad->mline[pos+i] == 0) break;
      if (fprintf(fp,  "  %c  ", ad->mline[pos+i]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }
    if (fprintf(fp,  " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    if (fprintf(fp, "\n") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");

    /* Print Mandatory Translation Line */
    if(spliced_ali ) {
      if (sprintf(exon_name, "%s %d", "exon", exon_cnt) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      if (fprintf(fp, "  %*s", namewidth, exon_name)      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      if (fprintf(fp, " %*s ", coordwidth, "") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }
    else {
      if (fprintf(fp, "  %*s ", namewidth+coordwidth+1, " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }
    if (fprintf(fp,  " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    for (i = 0, y = 0; i < cur_aliwidth; i++, y++)
    {
      if (ad->aseq[pos+i] == 0) break;
      if (fprintf(fp, "  %c  ", ad->aseq[pos+i]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }
    if (fprintf(fp,  " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    if (fprintf(fp, "\n") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");

    /* Print Mandatory Target Line */
    if (fprintf(fp, "  %*s", namewidth, show_seqname)    < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    if (ni > 0) { if (fprintf(fp, " %*ld ", coordwidth, i1) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed"); }
    else        { if (fprintf(fp, " %*s ", coordwidth, "-") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed"); }

    /* Print Acceptor Site Signal or Space */
    if (old_splice_opt  == 1 || old_splice_opt  == 4) {
      if (fprintf(fp, "%c", ad->ntseq[old_splice_pos+4]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      i2++;
    }
    else {
      if (fprintf(fp,  " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }
    /* Print ntseq */
    npos = pos * 5;
    for (i = 0, j = 0; j < cur_aliwidth; i+=5, j++)
    {
      if (ad->ntseq[npos+i] ==  0)  break;
      if (ad->ntseq[npos+i] == '$') break;
      if (fprintf(fp, "%c%c%c%c%c", ad->ntseq[npos+i],ad->ntseq[npos+i+1],ad->ntseq[npos+i+2],ad->ntseq[npos+i+3],ad->ntseq[npos+i+4]) < 0)
        ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");

        if(ad->sqfrom < ad->sqto)       { c1 = i2;     i2 += (ad->codon[pos+j] == 6 ? 3 : ad->codon[pos+j]); }
        else                            { c1 = i2 - 1; i2 -= (ad->codon[pos+j] == 6 ? 3 : ad->codon[pos+j]); }

      /* Get info for optional frame line */
        if(frameline != NULL)
        {
          if(ad->codon[pos+j] == 0 || ad->codon[pos+j] == 6 ) frame = 0;
          else                      frame = diplay_frame2(c1, i2);
        frameline[j] = frame;
      }
    }
    if (splice_opt == 4) {
      if (fprintf(fp, "%c", ad->ntseq[npos+i]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }

    if (ni > 0) { if (fprintf(fp, " %-*ld\n", coordwidth, i2)  < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed"); }
    else        { if (fprintf(fp, " %*s\n",   coordwidth, "-") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed"); }


    /* Print Optional Frame line */
    if(frameline != NULL)
    {
      if (fprintf(fp, "  %*s ", namewidth+coordwidth+1, "")  < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      /* Print space at begining of line to make room for splice signal */
      if (fprintf(fp,  "  ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      for (i = 0; i < cur_aliwidth; i++)
      {
          if (ad->aseq[pos+i] == 0) break;
          if(frameline[i] > 0)            { if (fprintf(fp, "  %d  ", frameline[i]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed"); }
          else if(frameline[i] < 0)       { if (fprintf(fp, " %d  ", frameline[i]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed"); }
        else if(ad->codon[pos+i] == 6) { if (fprintf(fp, "  %d  ", frameline[i]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed"); }
          else                            { if (fprintf(fp, "  .  ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed"); }
      }
      if (fprintf(fp,  " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      if (fprintf(fp, " FRAME\n")  < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }

    /* Print Mandatory Post Prob Line */
    if (fprintf(fp, "  %*s ", namewidth+coordwidth+1, "")  < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    /* Print space or '|' for splice signal*/
    if (old_splice_opt  == 1 || old_splice_opt  == 4) {
      if (fprintf(fp, "|") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }
    else {
       if (fprintf(fp, " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }

    if (old_splice_opt == 1 || old_splice_opt  == 4) {
      if(cur_aliwidth > 2) {
        if (fprintf(fp,  "| %c  ", ad->ppline[pos]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      }
      else {
        if (fprintf(fp,  "| %c |", ad->ppline[pos]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      }
      i = 1;
    }
    else if (old_splice_opt == 2) {
      if (fprintf(fp,  "||%c  ", ad->ppline[pos]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      i = 1;
    }
    else if (old_splice_opt == 3) {
      if (fprintf(fp,  " ||  ")                   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      i = 1;
    }
    else i = 0;

    for (  ; i < cur_aliwidth-2; i++)
    {
      if (ad->ppline[pos+i] == 0) break;
      if (fprintf(fp,  "  %c  ", ad->ppline[pos+i])  < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }
      
    if (splice_opt == 1) {
      if(cur_aliwidth > 2) {
        if (fprintf(fp,  "  %c |", ad->ppline[pos+i])   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      }
      if (fprintf(fp,  "|    ")                       < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }
    else if (splice_opt == 2) {
      if (fprintf(fp,  "  %c  ", ad->ppline[pos+i])   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      if (fprintf(fp,  "  || ")                       < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }
    else if (splice_opt == 3) {
      if (fprintf(fp,  "  %c  ", ad->ppline[pos+i])   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      if (fprintf(fp,  "  %c||", ad->ppline[pos+i+1]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }
    else if (splice_opt == 4) {
        if (fprintf(fp,  "  %c  ", ad->ppline[pos+i])   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
        if (fprintf(fp,  "  %c |", ad->ppline[pos+i+1]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
        if (fprintf(fp,  "|")                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    }
    else {
      for (  ; i < cur_aliwidth; i++)
      {
        if (ad->ppline[pos+i] == 0) break;
        if (fprintf(fp,  "  %c  ", ad->ppline[pos+i])  < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      }

    }
    if (fprintf(fp,  " ") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
    if (fprintf(fp, " PP\n")  < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");

    k1 += nk;
    if   (ad->sqfrom < ad->sqto)  i1 = i2 + 1;
    else                          i1 = i2 - 1;

    pos += cur_aliwidth;
    if(is_splice_line) {
      pos++;                          // pass over splice pos with '$'
      if (splice_opt == 4) pos++;
      i1 = ad->exon_seq_starts[exon_cnt]; //get the nuclueotide state coord for next exon
      if   (ad->sqfrom < ad->sqto)  i2 = i1 - 1;
      else                          i2 = i1 + 1;
      if (fprintf(fp, "\n")  < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");
      if (fprintf(fp, "[ Exon Set %d / Exon %d ]\n", set_num ,exon_cnt)  < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "alignment display write failed");

      exon_cnt++;

    }

    old_splice_pos = splice_pos;
    old_splice_opt = splice_opt;
  }


  fflush(fp);
  free(buf);
  free(exon_name);
  return eslOK;

  ERROR:
    if (buf) free(buf);
    if(exon_name) free(exon_name);
    return status;
}


int dump_splash_header(SPLICE_PATH *path, ESL_SQ *target_seq,int exon_set_name_id,P7_OPROFILE *om, P7_ALIDISPLAY *ad, P7_PIPELINE *tmp_pli, FILE * ofp,int textw)
{
  int status;
  int exon_id, i;
  int num_exons;
  int model_start;
  int model_end;
  int nucl_start;
  int nucl_end;
  int exon_cnt;
  int full_coverage;
  int num_found_exons;
  int *ExonCoordSet;

 
  ESL_ALLOC(ExonCoordSet, sizeof(int)*path->path_len*5);


  num_found_exons = 0;
  exon_cnt  = 0;
  while (exon_cnt < path->path_len) {

    if (path->missing[exon_cnt]) num_found_exons++;

    /*hmm start */
    ExonCoordSet[5*exon_cnt+2] = path->downstream_spliced_amino_start[exon_cnt];

    /* seq start */
    if(path->revcomp)
      ExonCoordSet[5*exon_cnt+1] = target_seq->n - path->downstream_spliced_nuc_start[exon_cnt] + target_seq->end;
    else
      ExonCoordSet[5*exon_cnt+1] = path->downstream_spliced_nuc_start[exon_cnt] + target_seq->start - 1;

    /* hmm end */
    ExonCoordSet[5*exon_cnt+4] = path->upstream_spliced_amino_end[exon_cnt+1];

    /* seq end*/
    if(path->revcomp)
      ExonCoordSet[5*exon_cnt+3] = target_seq->n - path->upstream_spliced_nuc_end[exon_cnt+1] + target_seq->end;
    else
      ExonCoordSet[5*exon_cnt+3] = path->upstream_spliced_nuc_end[exon_cnt+1] + target_seq->start - 1;
    
   
    exon_cnt++;
  }

  num_exons   = path->path_len;
  model_start = ExonCoordSet[2];
  model_end   = ExonCoordSet[5*(exon_cnt-1)+4];
  nucl_start  = ExonCoordSet[1];
  nucl_end    = ExonCoordSet[5*(exon_cnt-1)+3];

  
  full_coverage = 0;
  if (model_start == 1 && model_end == om->M)
    full_coverage = 1;


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
    fprintf(ofp,"| = Exon %d: %d..%d / %d..%d\n",exon_id+1, ExonCoordSet[5*exon_id+2], ExonCoordSet[5*exon_id+4], ExonCoordSet[5*exon_id+1], ExonCoordSet[5*exon_id+3]);
if (num_found_exons)
    fprintf(ofp,"| + Includes Missed Exons\n");

  fprintf(ofp,"|\n");
 fprintf(ofp,":\n");

  Splice_Test_Print(ofp, ad, 100, 150, tmp_pli, exon_set_name_id);

  fprintf(ofp,":\n");
  fprintf(ofp,"|\n");
  fprintf(ofp,"+");
  for (i=0; i<textw-2; i++)
    fprintf(ofp,"=");
  fprintf(ofp,"+\n");
  fprintf(ofp,"\n\n");

  free(ExonCoordSet);
  return eslOK;

  ERROR:
   return status;
}



