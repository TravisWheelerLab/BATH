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

#include "easel.h"
#include "esl_vectorops.h"
#include "esl_gumbel.h"
#include "esl_exponential.h"

#include "hmmer.h"
#include "p7_splice.h"

static ESL_SQ* get_sub_sequence(const ESL_SQFILE *seq_file, char* seqname, int64_t seq_min, int64_t seq_max, int revcomp);
static P7_HMM* get_sub_hmm(const P7_HMM *hmm, int start, int end);
static P7_HIT** split_hit(P7_HIT *hit, const P7_HMM *hmm, const P7_BG *bg, ESL_SQ *hit_seq, const ESL_GENCODE *gcode, int revcomp, int *num_hits);
static int hit_upstream(P7_DOMAIN *upstream, P7_DOMAIN *downstream, int revcomp);
static int hits_spliceable(P7_DOMAIN *upstream, P7_DOMAIN *downstream);
static SPLICE_EDGE* connect_nodes_with_edge(SPLICE_PIPELINE *pli, const P7_DOMAIN *upstream_domain, const P7_DOMAIN *downstream_domain, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq, int revcomp);
static void get_overlap_nuc_coords (SPLICE_EDGE *edge, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const ESL_SQ *splice_seq, int revcomp);
static int find_optimal_splice_site (SPLICE_EDGE *edge, SPLICE_PIPELINE *pli, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq);
static int select_splice_option (SPLICE_EDGE *edge, P7_PROFILE *sub_model, P7_FS_PROFILE *sub_fs_model, const ESL_GENCODE *gcode, const ESL_SQ *target_seq, float signal_score, int up_nuc_pos, int down_nuc_pos);
static P7_HIT** align_the_gap(SPLICE_GRAPH *graph, const P7_HMM *hmm, const P7_BG *bg, ESL_SQ *gap_seq, const ESL_GENCODE *gcode, int seq_i, int gap_len, int hmm_start, int hmm_end, int *num_hits);
static int bridge_the_gap(SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, P7_HIT **gap_hits, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *gap_seq, int num_hits);


/*  Function: p7_splice_SpliceHits
 *  Synopsis: SPLASH's splcing pipeline
 *
 *  Purpose : Run the splicing pipeline on a collections of hits 
 *            <tophits> between a protein query model <gm> and a 
 *            nucleotide sequence from the target file <seq_file>.
 *
 * Returns:   <eslOK> on success. If hits are successfully splice
 *            the new spliced alignements, scores, etc will be
 *            included in the <tophits> to be reported.
 *         
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_splice_SpliceHits(P7_TOPHITS *tophits, SPLICE_SAVED_HITS *saved_hits, P7_HMM *hmm, P7_OPROFILE *om, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, P7_SCOREDATA *scoredata, ESL_GETOPTS *go, ESL_GENCODE *gcode, ESL_SQFILE *seq_file, FILE *ofp, int64_t db_nuc_cnt)
{

  int              i;
  int              hit_cnt;
  int              range_cnt;
  int              success;
  int              num_hits_processed; 
  int              prev_num_hits_processed;
  int              is_sorted_by_sortkey;
  int              revcomp, prev_revcomp;
  int              frameshift;
  int              first, last;
  int64_t          seqidx, prev_seqidx;
  int             *hits_processed;
  int             *hit_spliced;
  int64_t         *range_bound_mins;
  int64_t         *range_bound_maxs;
  SPLICE_PIPELINE *pli;
  SPLICE_GRAPH   *graph;
  SPLICE_PATH     *path;
  P7_HIT          *hit; 
  int              status;

  graph           = NULL;
  path             = NULL;
  hits_processed   = NULL;
  hit_spliced      = NULL;
  range_bound_mins = NULL;
  range_bound_maxs = NULL;

  p7_splice_RemoveDuplicates(saved_hits);

  /* Sort hits by sequence, strand and nucleotide positon. 
   * Note if P7_TOPHITS was perviously sorted by sortkey 
   * so we can return it to it's original state */  
  is_sorted_by_sortkey = tophits->is_sorted_by_sortkey;
  if ((status = p7_tophits_SortBySeqidxAndAlipos(tophits)) != eslOK) goto ERROR;
     
  pli = p7_splicepipeline_Create(go, om->M, om->M * 3);
  pli->bg = p7_bg_Create(om->abc);
  if (pli->do_biasfilter) p7_bg_SetFilter(pli->bg, om->M, om->compo);  

  /* Arrays to keep track of hit status */
  ESL_ALLOC(hits_processed,   tophits->N*sizeof(int));
  ESL_ALLOC(hit_spliced,      tophits->N*sizeof(int));

  /* Check for any hits that are duplicates or have low quailty alignment and mark them as processed */
  hit_cnt = 0;
  num_hits_processed = 0;
  for(i = 0; i < tophits->N; i++) {
    tophits->hit[i]->dcl->ad->exon_cnt = 1;
    if ((tophits->hit[i]->flags & p7_IS_DUPLICATE)) {
      hits_processed[i] = 1;
      num_hits_processed++;
    }
    else { 
      hits_processed[i] = 0;
      hit_cnt++;
    }  
    hit_spliced[i]   = 0;
  }
  
  /* Arrays for keeping track of allowable target range mins 
   * and maxes set by previous ranges or splicings */
  ESL_ALLOC(range_bound_mins, hit_cnt * sizeof(int64_t));
  ESL_ALLOC(range_bound_maxs, hit_cnt * sizeof(int64_t));
  
  prev_num_hits_processed = -1;
  range_cnt = 0;

  prev_revcomp = -1;
  prev_seqidx  = -1;
  printf("\nQuery %s LENG %d \n",  gm->name, gm->M);
   fflush(stdout);
  /* loop through until all hits have been processed */
  while(num_hits_processed < tophits->N) {
    
    if(prev_num_hits_processed == num_hits_processed)
      ESL_XEXCEPTION(eslFAIL, "p7_splice_SpliceHits : loop failed to process hits");
    prev_num_hits_processed = num_hits_processed;
  
    /* Find the first unprocessed hit (sorting will ensure that hits 
     * will be on same sequence and strand untill all are processed) */
    i = 0;
    while (i < tophits->N && hits_processed[i]) i++;
    hit = tophits->hit[i];
    seqidx = hit->seqidx;
    revcomp = 1;
    if (hit->dcl->iali < hit->dcl->jali)
      revcomp = 0;

    /* If we have a new sequence and strand we build a new target range and graph */
    if(seqidx != prev_seqidx || revcomp != prev_revcomp) {   
      /* Find the first and last index of the current seqidx and strand in the saved hits list */
      first = 0;
      while(first < saved_hits->N && (saved_hits->srt[first]->seqidx != seqidx || saved_hits->srt[first]->strand != revcomp)) first++;
      last = first;
      while(last < saved_hits->N && saved_hits->srt[last]->seqidx == seqidx && saved_hits->srt[last]->strand == revcomp) last++;
      last--;
 
      /* unless this is the first graph destroy all the old data */
      if( prev_seqidx != -1) {
  
        p7_splicegraph_Destroy(graph);
        graph = NULL;
        range_cnt = 0;
      }   

      if ((graph = p7_splicegraph_Create()) == NULL) goto ERROR;
      graph->revcomp = revcomp;

      if((p7_splice_AddOriginals(graph, tophits, hits_processed, seqidx)) != eslOK) goto ERROR;      

      p7_splicehits_FindNodes(graph, saved_hits, first, last);
      
      printf("\nQuery %s Target %s strand %c\n", gm->name, graph->seqname, (graph->revcomp ? '-' : '+'));
      fflush(stdout);

      p7_splicehits_Dump(stdout, saved_hits);

      if((p7_splice_SplitHits(graph, hmm, pli->bg, gcode, seq_file)) != eslOK) goto ERROR;
      p7_splicegraph_DumpHits(stdout, graph);        

      p7_splice_ConnectGraph(graph, pli, hmm, gcode, seq_file);

      p7_splice_FillGaps(graph, pli, hmm, gcode, seq_file);
      p7_splicegraph_DumpHits(stdout, graph);
      p7_splicegraph_DumpEdges(stdout, graph, FALSE);
    }
    
    //TODO
    num_hits_processed = tophits->N; 
    p7_splicepipeline_Reuse(pli);
  }
 
  /* Leave only footprints */
  if (is_sorted_by_sortkey)
    if ((status = p7_tophits_SortBySortkey(tophits)) != eslOK) goto ERROR;   

  if(pli   != NULL) p7_splicepipeline_Destroy(pli);
  if(graph != NULL) p7_splicegraph_Destroy(graph);
  

  if (hits_processed   != NULL) free(hits_processed);
  if (hit_spliced      != NULL) free(hit_spliced);
  if (range_bound_mins != NULL) free(range_bound_mins);
  if (range_bound_maxs != NULL) free(range_bound_maxs);


  return status;

  ERROR:
    if(pli   != NULL) p7_splicepipeline_Destroy(pli);
    if(graph != NULL) p7_splicegraph_Destroy(graph);

    if (hits_processed   != NULL) free(hits_processed);
    if (hit_spliced      != NULL) free(hit_spliced); 
    if (range_bound_mins != NULL) free(range_bound_mins);
    if (range_bound_maxs != NULL) free(range_bound_maxs);
    return status;
}

int 
p7_splice_AddOriginals(SPLICE_GRAPH *graph, const P7_TOPHITS *tophits, int *hits_processed, int64_t seqidx)
{

  int     i;
  int     hit_cnt;
  P7_HIT *curr_hit;

  /* Get a count of the hits that will be added to the graph */
  hit_cnt = 0;
  for (i = 0; i < tophits->N; i++) {
    curr_hit = tophits->hit[i];
    if (hits_processed[i]) continue;
    if (curr_hit->seqidx != seqidx) continue;
    if (graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
    if ((!graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;
    
    hit_cnt++;
  }

  p7_splicegraph_CreateNodes(graph, hit_cnt);    
  
  /*Add all hits from current sequence and strand to TARGET_RANGE */
  for (i = 0; i < tophits->N; i++) {

    curr_hit = tophits->hit[i];

    if (hits_processed[i]) continue;

    if (curr_hit->seqidx != seqidx) continue;

    if (graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
    if ((!graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;
  
    if(curr_hit->flags & p7_IS_REPORTED) graph->reportable[i] = TRUE;
    else                                 graph->reportable[i] = FALSE;

    graph->in_graph[graph->th->N]     = TRUE;
    graph->orig_hit_idx[graph->th->N] = i;
    
    p7_splicegraph_AddNode(graph, curr_hit); 
  }

  graph->orig_N  = graph->th->N;
  graph->seqidx  = seqidx;
  graph->seqname = graph->th->hit[0]->name;

  return eslOK;
 
}

int
p7_splice_SplitHits(SPLICE_GRAPH *graph, const P7_HMM *hmm, const P7_BG *bg, ESL_GENCODE *gcode, const ESL_SQFILE *seq_file) {

  int h, s, z;
  int z1, z2;
  int I_cnt;
  int low_pp_cnt;
  int splitable;
  int seq_min, seq_max;
  int num_hits;
  P7_HIT     *curr_hit;
  P7_TRACE   *tr;
  P7_TOPHITS *th;
  P7_HIT    **split_hits; 
  ESL_SQ     *hit_seq;

  th = graph->th;
 
  graph->split_N = graph->orig_N;
  for (h = 0; h < graph->orig_N; h++) {

    curr_hit = th->hit[h]; 

    /*Check if hit looks like it cointins and intron */
    tr = curr_hit->dcl->tr;

    I_cnt = 0;
    low_pp_cnt = 0;
    splitable = FALSE;
    
    for (z1 = 0;       z1 < tr->N; z1++) if ( tr->st[z1] == p7T_M ) break;
    for (z2 = tr->N-1; z2 >= 0;    z2--) if ( tr->st[z2] == p7T_M ) break;

    z = z1;
    while (z < z2) {
      if (tr->st[z] == p7T_I) {
        I_cnt++;
      }
      else if(tr->st[z] == p7T_M && tr->pp[z] < 0.45) {
        low_pp_cnt++;
      }
      else if (I_cnt + low_pp_cnt > 10) {
        splitable = TRUE;
        break;
      }
      else {
        I_cnt = 0;
        low_pp_cnt = 0;
      }
      z++;
    }

    if(splitable) {
    
      seq_min = ESL_MIN(curr_hit->dcl->iali, curr_hit->dcl->jali);
      seq_max = ESL_MAX(curr_hit->dcl->iali, curr_hit->dcl->jali);
      hit_seq = get_sub_sequence(seq_file, graph->seqname, seq_min, seq_max, graph->revcomp);      
                
      split_hits = split_hit(curr_hit, hmm, bg, hit_seq, gcode, graph->revcomp, &num_hits);      
    
      if(num_hits == 1) { 
        
        p7_alidisplay_Destroy(split_hits[0]->dcl->ad);
        p7_trace_fs_Destroy(split_hits[0]->dcl->tr);
        free(split_hits[0]->dcl->scores_per_pos);
        p7_hit_Destroy(split_hits[0]);

        free(split_hits);
        esl_sq_Destroy(hit_seq);
        continue;
      }
 
      /* If the hit was split, set the original hit as not in graph and add new split hits to graph */ 
      graph->in_graph[h] = FALSE;
      
      for(s = 0; s < num_hits; s++) {
     
        p7_splicegraph_Grow(graph);
        graph->in_graph[graph->th->N]     = TRUE;
        graph->orig_hit_idx[graph->th->N] = graph->orig_hit_idx[h];

        p7_splicegraph_AddNode(graph, split_hits[s]);
        
        graph->split_N++; 
        
      } 
      
      free(split_hits); 
      esl_sq_Destroy(hit_seq);
      
    }
    
    
  }
  

  return eslOK;
}



P7_HIT**
split_hit(P7_HIT *hit, const P7_HMM *hmm, const P7_BG *bg, ESL_SQ *hit_seq, const ESL_GENCODE *gcode, int revcomp, int *num_hits)
{
 
  int         i, y, z;
  int         z1, z2;
  int         c;
  int         intron_cnt;
  int         hit_cnt;
  int         start_new;
  int         ihmm, jhmm;
  int         iali, jali;
  P7_HMM       *sub_hmm;
  P7_FS_PROFILE *sub_fs_model;
  P7_GMX       *vit_mx;
  P7_HIT       *new_hit;
  P7_HIT      **ret_hits; 
  P7_TRACE     *tr; 
  int status;

  sub_hmm     = get_sub_hmm(hmm, hit->dcl->ihmm, hit->dcl->jhmm);
  sub_hmm->fs = 0.;
  
  sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
  p7_ProfileConfig_fs(sub_hmm, bg, gcode, sub_fs_model, hit_seq->n*3, p7_GLOBAL); //hit_seq->n*3 to get single nucleotide loops
  
  vit_mx = p7_gmx_fs_Create(sub_fs_model->M, hit_seq->n, hit_seq->n, 2);
  tr = p7_trace_fs_Create();
  
  p7_sp_trans_semiglobal_Viterbi(hit_seq->dsq, gcode, hit_seq->n, sub_fs_model, vit_mx);
  p7_sp_trans_semiglobal_VTrace(hit_seq->dsq, hit_seq->n, gcode, sub_fs_model, vit_mx, tr); 
 
  /* Find number of introns in trace */
  intron_cnt = 0; 
  for(z = 0; z < tr->N; z++)
    if(tr->st[z] == p7T_R) intron_cnt++;

  ret_hits = NULL;
  ESL_ALLOC(ret_hits, sizeof(P7_HIT*) * (intron_cnt+1));
 
  /* Find first M state - start of first hit */
  for(z1 = 0; z1 < tr->N; z1++) if(tr->st[z1] == p7T_M) break;
  
  /* Find last M state state - end of last hit */
  for(z2 = tr->N-1; z1 >= 0; z2--) if(tr->st[z2] == p7T_M) break;

  hit_cnt = 0;
  start_new = TRUE;
  z = z1;
  while(z < z2) {

    if(start_new) {
 
      /* Save z value - currently set to fist M state in exon */
      y = z;

      /*Find end of exon */
      while(tr->st[z] != p7T_R && tr->st[z] != p7T_E) z++;
     
     /*Trace back to last M state of exon*/
      while(tr->st[z] != p7T_M) z--;

      /* Get exon coords */
      ihmm = tr->k[y] + hit->dcl->ihmm - 1;  
      jhmm = tr->k[z] + hit->dcl->ihmm - 1;
  
      if(revcomp) {
        iali = hit_seq->n - tr->i[y] + hit_seq->end + 2;
        jali = hit_seq->n - tr->i[z] + hit_seq->end;
      }
      else {
        iali = hit_seq->start + tr->i[y] - 3;
        jali = hit_seq->start + tr->i[z] - 1;
      }
       
      /* Create new hit and  set ihmm and iali coords*/
      new_hit          = p7_hit_Create_empty();
      
      new_hit->dcl     = p7_domain_Create_empty();
      new_hit->dcl->tr = p7_trace_fs_Create();
    
      new_hit->dcl->ihmm = ihmm;
      new_hit->dcl->jhmm = jhmm;
      new_hit->dcl->iali = iali;
      new_hit->dcl->jali = jali;
   
      /* Append starting special states */
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_S , 0, 0, 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_N , 0, tr->i[y]-3, 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_B , 0, tr->i[y]-3, 0);
      
      /* Append all states between first and last M state */
      for(i = y; i <= z; i++) {
        c = (tr->st[i] == p7T_M ? 3 : 0);
        p7_trace_fs_Append(new_hit->dcl->tr, tr->st[i] , tr->k[i], tr->i[i], c);
      }

      /* Append ending special states */
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_E, tr->k[z], tr->i[z], 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_C, 0, tr->i[z], 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_T , 0, 0, 0);

      new_hit->dcl->ad = p7_alidisplay_fs_Create(new_hit->dcl->tr, 0, sub_fs_model, hit_seq, gcode);
      new_hit->dcl->scores_per_pos = NULL;
      p7_splice_ComputeAliScores_fs(new_hit->dcl, new_hit->dcl->tr, hit_seq->dsq, sub_fs_model, hit_seq->abc);

      ret_hits[hit_cnt] = new_hit;
      hit_cnt++;

      
      start_new = FALSE;
    }
 
    z++;
    if(tr->st[z] == p7T_M) start_new = TRUE;     
     
  }
 
  *num_hits = hit_cnt;
 
  p7_hmm_Destroy(sub_hmm);
  p7_profile_fs_Destroy(sub_fs_model);
  p7_gmx_Destroy(vit_mx);
  p7_trace_fs_Destroy(tr);

  return ret_hits;


  ERROR:
   if(ret_hits != NULL) free(ret_hits);
   return NULL;

}

int
p7_splice_ConnectGraph(SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQFILE *seq_file) 
{
  
  int up, down;
  int up_min, up_max;
  int down_min, down_max;
  int seq_min, seq_max;
  P7_TOPHITS *th;
  ESL_SQ *splice_seq;
  SPLICE_EDGE *edge;

  th = graph->th;

  for(up = 0; up < graph->num_nodes; up++) {
    if(!graph->in_graph[up]) continue;
    for(down = 0; down < graph->num_nodes; down++) {
      if(!graph->in_graph[down]) continue;

      if(!hit_upstream(th->hit[up]->dcl, th->hit[down]->dcl, graph->revcomp)) continue;
      if(p7_splicegraph_EdgeExists(graph,up,down)) continue;

      /* check that hmm and ali coords are spliceable */
      if(!hits_spliceable(th->hit[up]->dcl, th->hit[down]->dcl)) continue;

      /* get min and max coords and fetch sub sequence of the potential splice region*/
      up_min   = ESL_MIN(th->hit[up]->dcl->iali, th->hit[up]->dcl->jali);
      up_max   = ESL_MAX(th->hit[up]->dcl->iali, th->hit[up]->dcl->jali);

      down_min = ESL_MIN(th->hit[down]->dcl->iali, th->hit[down]->dcl->jali);
      down_max = ESL_MAX(th->hit[down]->dcl->iali, th->hit[down]->dcl->jali);

      seq_min  = ESL_MIN(up_min, down_min);
      seq_max  = ESL_MAX(up_max, down_max);

      splice_seq = get_sub_sequence(seq_file, graph->seqname, seq_min, seq_max, graph->revcomp);

      edge = connect_nodes_with_edge(pli, th->hit[up]->dcl, th->hit[down]->dcl, hmm, gcode, splice_seq, graph->revcomp);
 
      if(edge != NULL) {
      
        edge->upstream_node_id   = up;
        edge->downstream_node_id = down;
               
        p7_splicegraph_AddEdge(graph, edge);
      }
      esl_sq_Destroy(splice_seq);   
    }

  }

  return eslOK;

}

int
hit_upstream(P7_DOMAIN *upstream, P7_DOMAIN *downstream, int revcomp) {

  if(upstream->ihmm >= downstream->ihmm && upstream->jhmm >= downstream->jhmm)
    return FALSE;

  if (( revcomp  && downstream->jali > upstream->iali) ||
     ((!revcomp) && downstream->jali < upstream->iali))
    return FALSE;

  return TRUE;
  
}

int
hits_spliceable(P7_DOMAIN *upstream, P7_DOMAIN *downstream) {
  
 
  /* Are nodes close enough on hmm coords */ 
  if(upstream->jhmm + MAX_AMINO_EXT < downstream->ihmm) return FALSE;
 
  /* Are nodes close enough on seq coords */
  if(abs(downstream->iali - upstream->jali) + 1 > MAX_INTRON_LEN) return FALSE;

  /* Do nodes overlap by too much */
  if(upstream->jhmm - downstream->ihmm + 1 > MAX_AMINO_OVERLAP) return FALSE;

  return TRUE;
  
}

SPLICE_EDGE*
connect_nodes_with_edge(SPLICE_PIPELINE *pli, const P7_DOMAIN *upstream_domain, const P7_DOMAIN *downstream_domain, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq, int revcomp)
{

  int          up_amino_start, up_amino_end;
  int          down_amino_start, down_amino_end;
  int          overlap;
  int          num_ext_aminos;
  SPLICE_EDGE *edge;
  int          status;

  up_amino_start = upstream_domain->ihmm;
  up_amino_end   = upstream_domain->jhmm;

  down_amino_start = downstream_domain->ihmm;
  down_amino_end   = downstream_domain->jhmm;

  /* Ensure overlap region is at least MIN_AMINO_OVERLAP long plus any gap extention positions */
  overlap = up_amino_end - down_amino_start + 1;
  if(overlap >= MIN_AMINO_OVERLAP)  // hits naturally have sufficent overlap
    num_ext_aminos = 0;
  else if( overlap > 0)             //hits overlap nuterally bu still need to be extended
    num_ext_aminos = MIN_AMINO_OVERLAP - overlap;
  else                              // hits do not overlap - gap extention + MIN_AMINO_OVERLAP needed
    num_ext_aminos = (down_amino_start - upstream_domain->jhmm)*2 + MIN_AMINO_OVERLAP;

  /* If the nodes have reached this point we will give them an edge*/
  edge = p7_splicegraph_CreateEdge();

  edge->overlap_amino_start = down_amino_start;
  edge->overlap_amino_end   = up_amino_end;

  /* If the hits do not overlap by at least MIN_AMINO_OVERLAP hmm positions, extend them */
  if (num_ext_aminos > 0) {
    num_ext_aminos = (num_ext_aminos+1)/2;
    edge->overlap_amino_start -= num_ext_aminos;
    edge->overlap_amino_end   += num_ext_aminos;
  }

  if(edge->overlap_amino_end > down_amino_end)   edge->overlap_amino_end = down_amino_end;
  if(edge->overlap_amino_start < up_amino_start) edge->overlap_amino_start = up_amino_start;

  get_overlap_nuc_coords(edge, upstream_domain, downstream_domain, splice_seq, revcomp);

  /* Add extra nucleotides for splice sites */
  edge->upstream_nuc_end     = ESL_MIN(edge->upstream_nuc_end + 2, splice_seq->n);
  edge->downstream_nuc_start = ESL_MAX(edge->downstream_nuc_start - 2, 1);

  if ((status = find_optimal_splice_site (edge, pli, upstream_domain, downstream_domain, hmm, gcode, splice_seq)) != eslOK) goto ERROR;
  
  if(edge->splice_score == -eslINFINITY) {
     free(edge);
     return NULL;
  }

  if(revcomp) {
    edge->upstream_spliced_nuc_end     = splice_seq->n - edge->upstream_spliced_nuc_end     + splice_seq->end;
    edge->downstream_spliced_nuc_start = splice_seq->n - edge->downstream_spliced_nuc_start + splice_seq->end;
  }
  else {
    edge->upstream_spliced_nuc_end     = edge->upstream_spliced_nuc_end     + splice_seq->start - 1;
    edge->downstream_spliced_nuc_start = edge->downstream_spliced_nuc_start + splice_seq->start - 1;
  }


  return edge;

  ERROR:
    if(edge != NULL) free(edge);
    return NULL;
}


void
get_overlap_nuc_coords (SPLICE_EDGE *edge, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const ESL_SQ *splice_seq, int revcomp)
{

  int       z1,z2;
  int       strand;
  int       extention_nuc_count;
  int       curr_hmm_pos;
  P7_TRACE *up_trace;
  P7_TRACE *down_trace;
  strand = (revcomp? -1 : 1);

  up_trace = upstream->tr;
  down_trace = downstream->tr;

 /***************UPSTREAM*******************/

 /* Get number of nucleotides that need to be added to the end of the
  * upstream edge nucleotide range if the end of the amino range
  * had to be extended to meet the mimimum overlap */
  extention_nuc_count = (edge->overlap_amino_end - upstream->jhmm) * strand * 3;
  edge->upstream_nuc_end = upstream->jali + extention_nuc_count;

  /* Set edge->upstream_nucl_start to the begining of the upstream
   * nulceotides that correpond to overlapping aminos */
  edge->upstream_nuc_start = upstream->jali + strand;

  for (z2 = up_trace->N-1 ; z2 >= 0; z2--) if (up_trace->st[z2] == p7T_M) break;

  curr_hmm_pos = upstream->jhmm;
  edge->upstream_trace_end = z2;

  while(curr_hmm_pos >= edge->overlap_amino_start) {

    if      (up_trace->st[z2] == p7T_M) {
      edge->upstream_nuc_start -= up_trace->c[z2] * strand;
      if(up_trace->c[z2] != 3) edge->frameshift = TRUE;
      curr_hmm_pos--;
    }
    else if (up_trace->st[z2] == p7T_I)
      edge->upstream_nuc_start -= 3 * strand;
    else if (up_trace->st[z2] == p7T_D)
      curr_hmm_pos--;
    else
      break;
    z2--;

  }
  edge->upstream_trace_start = z2+1;

  if(revcomp) {
    edge->upstream_nuc_start = splice_seq->start - edge->upstream_nuc_start + 1;
    edge->upstream_nuc_end   = splice_seq->start - edge->upstream_nuc_end   + 1;
  }
  else {
    edge->upstream_nuc_start = edge->upstream_nuc_start - splice_seq->start + 1;
    edge->upstream_nuc_end   = edge->upstream_nuc_end   - splice_seq->start + 1;
  }

/***************DOWNSTREAM*******************/

 /* Get number of nucleotides that need to be added to the start of the
  * downstream edge nucleotide range if the end of the amino range
  * had to be extended to meet the mimimum overlap */
  extention_nuc_count = (downstream->ihmm - edge->overlap_amino_start) * strand * 3;
  edge->downstream_nuc_start = downstream->iali - extention_nuc_count;


 /* Set edge->downstream_nucl_end to the end of the downstream
  * nulceotides that correpond to overlapping aminos */
  edge->downstream_nuc_end   = downstream->iali - strand;

  for (z1 = 0; z1 < down_trace->N; z1++) if (down_trace->st[z1] == p7T_M) break;

  curr_hmm_pos = downstream->ihmm;
  edge->downstream_trace_start = z1;
  while(curr_hmm_pos <= edge->overlap_amino_end) {
    if (down_trace->st[z1] == p7T_M) {
      edge->downstream_nuc_end += down_trace->c[z1] * strand;
      if(down_trace->c[z1] != 3) edge->frameshift = TRUE;
      curr_hmm_pos++;
    }
    else if (down_trace->st[z1] == p7T_I)
      edge->downstream_nuc_end += 3 * strand;
    else if (down_trace->st[z1] == p7T_D)
      curr_hmm_pos++;
    else
      break;
    z1++;
  }
  edge->downstream_trace_end = z1;

  if(revcomp) {
    edge->downstream_nuc_start = splice_seq->start - edge->downstream_nuc_start + 1;
    edge->downstream_nuc_end   = splice_seq->start - edge->downstream_nuc_end   + 1;
  }
  else {
    edge->downstream_nuc_start = edge->downstream_nuc_start - splice_seq->start + 1;
    edge->downstream_nuc_end   = edge->downstream_nuc_end   - splice_seq->start + 1;
  }

 return;
}


int
find_optimal_splice_site (SPLICE_EDGE *edge, SPLICE_PIPELINE *pli, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq)
{

  int            n,z;
  int            unp,dnp;
  float          lost_sc;
  float         *tp_0;
  float         *tp_M;
  P7_HMM        *sub_hmm;
  P7_PROFILE    *sub_model;
  P7_FS_PROFILE *sub_fs_model;
  int            status;

  /* Get the summed ali score for the overlap region covered by the existing alignments */
  lost_sc = 0.;
  z = edge->upstream_trace_end;
  n = z - (upstream->tr->tfrom[0]+1);
  while(upstream->tr->k[z] >= edge->overlap_amino_start) {
    lost_sc += upstream->scores_per_pos[n];
    n--;
    z--;
  }

  z = edge->downstream_trace_start;
  n = 0;

  while(downstream->tr->k[z] && downstream->tr->k[z] <= edge->overlap_amino_end) {
    lost_sc += downstream->scores_per_pos[n];
    n++;
    z++;
  }

  /*Get a submodel that covers the overlap region */
  sub_hmm    = get_sub_hmm(hmm, edge->overlap_amino_start, edge->overlap_amino_end);
  sub_model = NULL;
  sub_fs_model = NULL;

  if(edge->frameshift) {
    sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
    ESL_REALLOC(sub_fs_model->tsc,  sizeof(float)   * (sub_hmm->M+1) * p7P_NTRANS);
    p7_ProfileConfig_fs(sub_hmm, pli->bg, gcode, sub_fs_model, 100, p7_UNIGLOBAL);
    tp_0 = sub_fs_model->tsc;
    tp_M = sub_fs_model->tsc + sub_fs_model->M * p7P_NTRANS;

    /* Add extra transtion position for insetions at begining and end of alignments*/
    tp_0[p7P_MM] = log(hmm->t[edge->overlap_amino_start-1][p7H_MM]);
    tp_0[p7P_MI] = log(hmm->t[edge->overlap_amino_start-1][p7H_MI]); 
    tp_0[p7P_MD] = log(hmm->t[edge->overlap_amino_start-1][p7H_MD]);
    tp_0[p7P_IM] = log(hmm->t[edge->overlap_amino_start-1][p7H_IM]);
    tp_0[p7P_II] = log(hmm->t[edge->overlap_amino_start-1][p7H_II]);
    tp_0[p7P_DM] = log(hmm->t[edge->overlap_amino_start-1][p7H_DM]);
    tp_0[p7P_DD] = log(hmm->t[edge->overlap_amino_start-1][p7H_DD]);


    tp_M[p7P_MM] = log(hmm->t[edge->overlap_amino_end][p7H_MM]); 
    tp_M[p7P_MI] = log(hmm->t[edge->overlap_amino_end][p7H_MI]);
    tp_M[p7P_MD] = log(hmm->t[edge->overlap_amino_end][p7H_MD]);
    tp_M[p7P_IM] = log(hmm->t[edge->overlap_amino_end][p7H_IM]);
    tp_M[p7P_II] = log(hmm->t[edge->overlap_amino_end][p7H_II]);
    tp_M[p7P_DM] = log(hmm->t[edge->overlap_amino_end][p7H_DM]);
    tp_M[p7P_DD] = log(hmm->t[edge->overlap_amino_end][p7H_DD]); 
  }

  sub_model  = p7_profile_Create (sub_hmm->M, sub_hmm->abc);
  ESL_REALLOC(sub_model->tsc,  sizeof(float)   * (sub_hmm->M+1) * p7P_NTRANS);
  p7_ProfileConfig(sub_hmm, pli->bg, sub_model, 100, p7_UNIGLOBAL);
  tp_0 = sub_model->tsc;
  tp_M = sub_model->tsc + sub_model->M * p7P_NTRANS;

  /* Add extra transtion position for insetions at begining and end of alignments*/
  tp_0[p7P_MM] = log(hmm->t[edge->overlap_amino_start-1][p7H_MM]); 
  tp_0[p7P_MI] = log(hmm->t[edge->overlap_amino_start-1][p7H_MI]);
  tp_0[p7P_MD] = log(hmm->t[edge->overlap_amino_start-1][p7H_MD]);
  tp_0[p7P_IM] = log(hmm->t[edge->overlap_amino_start-1][p7H_IM]);
  tp_0[p7P_II] = log(hmm->t[edge->overlap_amino_start-1][p7H_II]);
  tp_0[p7P_DM] = log(hmm->t[edge->overlap_amino_start-1][p7H_DM]);
  tp_0[p7P_DD] = log(hmm->t[edge->overlap_amino_start-1][p7H_DD]);


  tp_M[p7P_MM] = log(hmm->t[edge->overlap_amino_end][p7H_MM]); 
  tp_M[p7P_MI] = log(hmm->t[edge->overlap_amino_end][p7H_MI]);
  tp_M[p7P_MD] = log(hmm->t[edge->overlap_amino_end][p7H_MD]);
  tp_M[p7P_IM] = log(hmm->t[edge->overlap_amino_end][p7H_IM]);
  tp_M[p7P_II] = log(hmm->t[edge->overlap_amino_end][p7H_II]);
  tp_M[p7P_DM] = log(hmm->t[edge->overlap_amino_end][p7H_DM]);
  tp_M[p7P_DD] = log(hmm->t[edge->overlap_amino_end][p7H_DD]);

  /* Scan the overlap nucleotides for splice signals */
  for(unp = edge->upstream_nuc_start; unp < edge->upstream_nuc_end; unp++) {
    /*GT-AG*/
    if(splice_seq->dsq[unp] == 2 && splice_seq->dsq[unp+1] == 3) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {
        if(splice_seq->dsq[dnp-1] == 0 && splice_seq->dsq[dnp] == 2) {
          if (dnp - unp > MIN_INTRON_LEN)
            select_splice_option(edge, sub_model, sub_fs_model, gcode, splice_seq, pli->signal_scores[p7S_GTAG], unp-1, dnp+1);
        }
      }
    }
    /*GC-AG*/
    if(splice_seq->dsq[unp] == 2 && splice_seq->dsq[unp+1] == 1) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {
        if(splice_seq->dsq[dnp-1] == 0 && splice_seq->dsq[dnp] == 2) {
          if (dnp - unp > MIN_INTRON_LEN)
            select_splice_option(edge, sub_model, sub_fs_model, gcode, splice_seq, pli->signal_scores[p7S_GCAG], unp-1, dnp+1);
        }
      }
    }
    /*AT-AC*/
    if(splice_seq->dsq[unp] == 0 && splice_seq->dsq[unp+1] == 3) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {
        if(splice_seq->dsq[dnp-1] == 0 && splice_seq->dsq[dnp] == 1) {
          if (dnp - unp > MIN_INTRON_LEN)
            select_splice_option(edge, sub_model, sub_fs_model, gcode, splice_seq, pli->signal_scores[p7S_ATAC], unp-1, dnp+1);
        }
      }
    }
  }

  if(edge->splice_score != -eslINFINITY)
    edge->splice_score -= lost_sc;

  if(sub_hmm       != NULL) p7_hmm_Destroy(sub_hmm);
  if(sub_model     != NULL) p7_profile_Destroy(sub_model);
  if(sub_fs_model  != NULL) p7_profile_fs_Destroy(sub_fs_model);

  return eslOK;

  ERROR:
    if(sub_hmm       != NULL) p7_hmm_Destroy(sub_hmm);
    if(sub_model     != NULL) p7_profile_Destroy(sub_model);
    if(sub_fs_model  != NULL) p7_profile_fs_Destroy(sub_fs_model);
    return status;

}


int
select_splice_option (SPLICE_EDGE *edge, P7_PROFILE *sub_model, P7_FS_PROFILE *sub_fs_model, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq, float signal_score, int up_nuc_pos, int down_nuc_pos)
{

  int         i,k,z;
  int         nuc_seq_idx;
  int         nuc_seq_len;
  int         amino_len;
  int         amino;
  int         upstream_nuc_cnt;
  int         upstream_amino_cnt;
  int         upstream_amino_end;
  float       sum_ali_sc;
  float       vitsc;
  float       overlap_sc;
  ESL_DSQ    *nuc_dsq;
  ESL_DSQ    *amino_dsq;
  P7_GMX     *vit_mx;
  P7_TRACE   *tr;
  int         status;

  nuc_dsq   = NULL;
  amino_dsq = NULL;
  vit_mx    = NULL;
  tr        = NULL;

  /*Get the overlap nucleotides that correspond to the current splice signal*/
  nuc_seq_len  = up_nuc_pos - edge->upstream_nuc_start + 1;
  nuc_seq_len += edge->downstream_nuc_end - down_nuc_pos + 1;

  /* If the splice signal is right at the overlap boundries all amino acid positions are deletions */
  if(nuc_seq_len == 0) {
    sum_ali_sc = 0;
    for(k = 0; k < sub_model->M; k++)
      sum_ali_sc += sub_model->tsc[k * p7P_NTRANS + p7H_DD];

    overlap_sc = sum_ali_sc + signal_score;

    if (overlap_sc > edge->splice_score) {
      edge->signal_score = signal_score;
      edge->splice_score = overlap_sc;
      edge->upstream_spliced_nuc_end = up_nuc_pos;
      edge->downstream_spliced_nuc_start = down_nuc_pos;
      edge->upstream_spliced_amino_end = edge->overlap_amino_start - 1;
      edge->downstream_spliced_amino_start = edge->overlap_amino_start;
    }
    return eslOK;
  }


  /* When not using the frameshift algorithm we can only
   * allow mod three length nucelotide sequences  */

  if((!edge->frameshift) && nuc_seq_len % 3) return eslOK;

  ESL_ALLOC(nuc_dsq,   sizeof(ESL_DSQ) * (nuc_seq_len+2));
  nuc_seq_idx = 0;
  nuc_dsq[nuc_seq_idx] = eslDSQ_SENTINEL;
  nuc_seq_idx++;
  for(i = edge->upstream_nuc_start; i <= up_nuc_pos; i++) {
    nuc_dsq[nuc_seq_idx] = splice_seq->dsq[i];
    nuc_seq_idx++;
  }
  upstream_nuc_cnt = nuc_seq_idx-1;
  upstream_amino_cnt = (nuc_seq_idx-1) / 3;
  if((nuc_seq_idx-1) % 3)  upstream_amino_cnt++;

  for (i = down_nuc_pos; i <= edge->downstream_nuc_end; i++) {
    nuc_dsq[nuc_seq_idx] = splice_seq->dsq[i];
    nuc_seq_idx++;
  }
  nuc_dsq[nuc_seq_idx] = eslDSQ_SENTINEL;

  if(nuc_seq_len % 3) {
     vit_mx = p7_gmx_fs_Create(sub_fs_model->M, nuc_seq_len, nuc_seq_len, p7P_CODONS);
     p7_fs_ReconfigLength(sub_fs_model, nuc_seq_len);
     p7_fs_global_Viterbi(nuc_dsq, gcode, nuc_seq_len, sub_fs_model, vit_mx, &vitsc);
  }
  else {
    /* Translate overalp nucleotides to amino sequence */
    amino_len = nuc_seq_len / 3;
    ESL_ALLOC(amino_dsq, sizeof(ESL_DSQ) * (amino_len+2));

    amino_dsq[0] = eslDSQ_SENTINEL;
    nuc_seq_idx = 1;
    for(i = 1; i <= amino_len; i++) {
      amino = esl_gencode_GetTranslation(gcode,&nuc_dsq[nuc_seq_idx]);
      amino_dsq[i] = amino;
      nuc_seq_idx+=3;
    }
    amino_dsq[amino_len+1] = eslDSQ_SENTINEL;

    /* Align translated overlap amino acids to submodel */
    vit_mx = p7_gmx_Create(sub_model->M, amino_len);
    p7_ReconfigLength(sub_model, amino_len);
    p7_global_Viterbi(amino_dsq, amino_len, sub_model, vit_mx, &vitsc);
  }

    /* vitsc from global alignment is the same as ali_score */
  if (vitsc != -eslINFINITY) {
    overlap_sc = vitsc + signal_score;

    if (overlap_sc > edge->splice_score) {
      /*use trace to find splice point*/

      if(nuc_seq_len % 3) {
        tr = p7_trace_fs_Create();
        p7_fs_global_Trace(nuc_dsq, nuc_seq_len, sub_fs_model, vit_mx, tr);

        z = 0;
        while(z < tr->N) {
          if(tr->i[z] >= upstream_nuc_cnt) {
            upstream_amino_end = tr->k[z] + edge->overlap_amino_start - 1;
            z = tr->N;
          }
          z++;
        }
      }
      else {
        tr = p7_trace_Create();
        p7_global_Trace(amino_dsq, amino_len, sub_model, vit_mx, tr);

        z = 0;
        while(z < tr->N) {
          if(tr->i[z] == upstream_amino_cnt) {
            upstream_amino_end = tr->k[z] + edge->overlap_amino_start - 1;
            z = tr->N;
          }
          z++;
        }
      }
      edge->signal_score = signal_score;
      edge->splice_score = overlap_sc;
      edge->upstream_spliced_nuc_end = up_nuc_pos;
      edge->downstream_spliced_nuc_start = down_nuc_pos;
      edge->upstream_spliced_amino_end = upstream_amino_end;
      edge->downstream_spliced_amino_start = upstream_amino_end + 1;
    }
  }
  /* Case where all translated aminos are stop codons - treat them as insertions and all hmm positions as deletions */
  else {
    for(i = 0; i < amino_len; i++)
      sum_ali_sc += sub_model->tsc[p7H_II];
    for(k = 0; k < sub_model->M; k++)
      sum_ali_sc += sub_model->tsc[k * p7P_NTRANS + p7H_DD];

    overlap_sc = sum_ali_sc + signal_score;

    if (overlap_sc > edge->splice_score) {
      edge->signal_score = signal_score;
      edge->splice_score = overlap_sc;
      edge->upstream_spliced_nuc_end = up_nuc_pos;
      edge->downstream_spliced_nuc_start = down_nuc_pos;
      edge->upstream_spliced_amino_end = edge->overlap_amino_start - 1;
      edge->downstream_spliced_amino_start = edge->overlap_amino_start;
    }
  }

  if(nuc_dsq   != NULL) free(nuc_dsq);
  if(amino_dsq != NULL) free(amino_dsq);
  if(vit_mx    != NULL) p7_gmx_Destroy(vit_mx);
  if(nuc_seq_len % 3) { if(tr != NULL) p7_trace_fs_Destroy(tr); }
  else                 { if(tr != NULL) p7_trace_Destroy(tr); }
  return eslOK;

  ERROR:
    if(nuc_dsq   != NULL) free(nuc_dsq);
    if(amino_dsq != NULL) free(amino_dsq);
    if(vit_mx    != NULL) p7_gmx_Destroy(vit_mx);
    if(nuc_seq_len % 3) { if(tr != NULL) p7_trace_fs_Destroy(tr); }
    else                 { if(tr != NULL) p7_trace_Destroy(tr); }
    return status;
}




int
p7_splice_FillGaps(SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQFILE *seq_file)
{

  int up, down;
  int i, j;
  int hmm_gap_len;
  int seq_gap_len;
  int hmm_overlap_min, hmm_overlap_max;
  int hmm_overlap_len;
  int seq_overlap_min, seq_overlap_max;
  int seq_overlap_len; 
  int seq_min, seq_max;
  int seq_i;
  int hit_cnt;
  int num_gaps;
  int gap_alloc; 
  int        *gap_merged;
  P7_TOPHITS *th;
  P7_HIT     **gap_hits;
  ESL_SQ     *gap_seq;
  SPLICE_GAP **gaps;
  SPLICE_GAP *new_gap; 
  int        status;   

  th = graph->th;

  num_gaps  = 0;
  gap_alloc = th->N;
  ESL_ALLOC(gaps, sizeof(SPLICE_GAP*) * gap_alloc);

  /* Find the total set of possible gaps */  
  for(up = 0; up < th->N; up++) {
    if(!graph->in_graph[up]) continue;

    for(down = 0; down < th->N; down++) {
      if(!graph->in_graph[down] || up == down) continue;
      
      if(!hit_upstream(th->hit[up]->dcl, th->hit[down]->dcl, graph->revcomp)) continue;
      if(p7_splicegraph_PathExists(graph, up, down))  continue;

      hmm_gap_len = th->hit[down]->dcl->ihmm - th->hit[up]->dcl->jhmm - 1;
      if(hmm_gap_len < 10) continue;
 
      if(graph->revcomp)
        seq_gap_len = th->hit[up]->dcl->jali - th->hit[down]->dcl->iali - 1;
      else  
        seq_gap_len = th->hit[down]->dcl->iali - th->hit[up]->dcl->jali - 1;
      if(seq_gap_len < 1 || seq_gap_len > MAX_GAP_RANGE) continue;
      
      
      if(num_gaps == gap_alloc) {
        gap_alloc*=2;
        ESL_REALLOC(gaps, sizeof(SPLICE_GAP*) * gap_alloc);
      }

      ESL_ALLOC(new_gap, sizeof(SPLICE_GAP));
      new_gap->upstream_node   = up;
      new_gap->downstream_node = down;
      new_gap->hmm_min = th->hit[up]->dcl->jhmm   + 1; 
      new_gap->hmm_max = th->hit[down]->dcl->ihmm - 1;  
      new_gap->seq_min = (graph->revcomp? th->hit[down]->dcl->iali + 1 : th->hit[up]->dcl->jali   + 1);
      new_gap->seq_max = (graph->revcomp? th->hit[up]->dcl->jali   - 1 : th->hit[down]->dcl->iali - 1);
      num_gaps++;

      gaps[num_gaps-1] = new_gap;
    }
  }    

  if(num_gaps == 0) {
    free(gaps);
    return eslOK;
  }

  ESL_ALLOC(gap_merged, sizeof(int) * num_gaps);
  esl_vec_ISet(gap_merged, num_gaps, 0);
  ESL_ALLOC(new_gap, sizeof(SPLICE_GAP));

  /* Merge overlapping gaps */
  for(i = 0; i < num_gaps; i++) {
    if(gap_merged[i]) continue;

    new_gap->upstream_node   = gaps[i]->upstream_node;
    new_gap->downstream_node = gaps[i]->downstream_node;

    new_gap->hmm_min = gaps[i]->hmm_min;
    new_gap->hmm_max = gaps[i]->hmm_max; 
    new_gap->seq_min = gaps[i]->seq_min;
    new_gap->seq_max = gaps[i]->seq_max;

    gap_merged[i] = TRUE;
    for(j = 0; j < num_gaps; j++) {
      if(gap_merged[j]) continue;
      
      hmm_overlap_min = ESL_MAX(new_gap->hmm_min, gaps[j]->hmm_min);
      hmm_overlap_max = ESL_MIN(new_gap->hmm_max, gaps[j]->hmm_max);
      
      hmm_overlap_len = hmm_overlap_max - hmm_overlap_min + 1;
      if(hmm_overlap_len < 1) continue;

      seq_overlap_min = ESL_MAX(new_gap->seq_min, gaps[j]->seq_min);
      seq_overlap_max = ESL_MIN(new_gap->seq_max, gaps[j]->seq_max);
 
      seq_overlap_len = seq_overlap_max - seq_overlap_min + 1;
      if(seq_overlap_len < 1 || seq_overlap_len > MAX_GAP_RANGE) continue;
      
      new_gap->hmm_min = ESL_MIN(new_gap->hmm_min, gaps[j]->hmm_min);
      new_gap->hmm_max = ESL_MAX(new_gap->hmm_max, gaps[j]->hmm_max);      
      new_gap->seq_min = ESL_MIN(new_gap->seq_min, gaps[j]->seq_min);      
      new_gap->seq_max = ESL_MAX(new_gap->seq_max, gaps[j]->seq_max); 

      if(graph->revcomp) {
        if(th->hit[gaps[j]->upstream_node]->dcl->iali   > th->hit[new_gap->upstream_node]->dcl->iali)
          new_gap->upstream_node = gaps[j]->upstream_node;
        if(th->hit[gaps[j]->downstream_node]->dcl->jali < th->hit[new_gap->downstream_node]->dcl->jali)
          new_gap->downstream_node = gaps[j]->downstream_node;
      }
      else {
        if(th->hit[gaps[j]->upstream_node]->dcl->iali   < th->hit[new_gap->upstream_node]->dcl->iali)
          new_gap->upstream_node = gaps[j]->upstream_node;
        if(th->hit[gaps[j]->downstream_node]->dcl->jali > th->hit[new_gap->downstream_node]->dcl->jali) 
          new_gap->downstream_node = gaps[j]->downstream_node;
      }

      gap_merged[j] = TRUE;
    }

    /* align merged gap */
    seq_min = ESL_MIN(th->hit[new_gap->upstream_node]->dcl->iali, th->hit[new_gap->downstream_node]->dcl->jali);
    seq_max = ESL_MAX(th->hit[new_gap->upstream_node]->dcl->iali, th->hit[new_gap->downstream_node]->dcl->jali);
    gap_seq = get_sub_sequence(seq_file, graph->seqname, seq_min, seq_max, graph->revcomp);  

    if(graph->revcomp) seq_i = seq_max - th->hit[new_gap->upstream_node]->dcl->jali + 1;
    else               seq_i = th->hit[new_gap->upstream_node]->dcl->jali - seq_min + 1;

    seq_gap_len = new_gap->seq_max - new_gap->seq_min + 1;    
    printf("seq_min %d seq_max %d hmm_min %d hmm_max %d\n", new_gap->seq_min, new_gap->seq_max, new_gap->hmm_min, new_gap->hmm_max);
    gap_hits = align_the_gap(graph, hmm, pli->bg, gap_seq, gcode, seq_i, seq_gap_len, new_gap->hmm_min, new_gap->hmm_max, &hit_cnt);
   
    if(hit_cnt) bridge_the_gap(graph, pli, gap_hits, hmm, gcode, gap_seq, hit_cnt);

    free(gap_hits);
    esl_sq_Destroy(gap_seq);
  } 
  
  for (i = 0; i < num_gaps; i++) free(gaps[i]);
  if(new_gap    != NULL) free(new_gap); 
  if(gaps       != NULL) free(gaps);
  if(gap_merged != NULL) free(gap_merged);
  return eslOK;

  ERROR:
    for (i = 0; i < num_gaps; i++) free(gaps[i]);
    if(new_gap    != NULL) free(new_gap);
    if(gaps       != NULL) free(gaps);
    if(gap_merged != NULL) free(gap_merged);
    return status; 
}



P7_HIT**
align_the_gap(SPLICE_GRAPH *graph, const P7_HMM *hmm, const P7_BG *bg, ESL_SQ *gap_seq, const ESL_GENCODE *gcode, int seq_i, int gap_len, int hmm_start, int hmm_end, int *num_hits) {

  int h, i, y, z;
  int z1, z2;
  int intron_cnt;
  int hit_cnt;
  int start_new;
  int ihmm, jhmm;
  int iali, jali;
  int codon;
  int duplicate;
  int exon_min, exon_max;
  int hit_min, hit_max;
  P7_TOPHITS    *th;
  P7_HMM        *sub_hmm;
  P7_FS_PROFILE *sub_fs_model;
  P7_GMX        *vit_mx;
  P7_TRACE      *tr;
  P7_HIT      **gap_hits;
  P7_HIT       *new_hit; 
  int           status;

  th = graph->th;

  sub_hmm = get_sub_hmm(hmm, hmm_start, hmm_end);
  sub_hmm->fs = 0.;

  sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
  p7_ProfileConfig_fs(sub_hmm, bg, gcode, sub_fs_model, gap_len*3, p7_GLOBAL); //len*3 to get single nucleotide loops

  vit_mx = p7_gmx_fs_Create(sub_fs_model->M, gap_len, gap_len, p7P_SPLICE);
  tr = p7_trace_fs_Create();

  p7_sp_trans_semiglobal_Viterbi(gap_seq->dsq+seq_i, gcode, gap_len, sub_fs_model, vit_mx);
  p7_sp_trans_semiglobal_VTrace(gap_seq->dsq+seq_i, gap_len, gcode, sub_fs_model, vit_mx, tr);

  /* Find number of introns in trace */
  intron_cnt = 0;
  for(z = 0; z < tr->N; z++)
    if(tr->st[z] == p7T_R) intron_cnt++;

  gap_hits = NULL;
  ESL_ALLOC(gap_hits, sizeof(P7_HIT*) * (intron_cnt+1));

  /* Find first M state - start of first hit */
  for(z1 = 0; z1 < tr->N; z1++) if(tr->st[z1] == p7T_M) break;

  /* Find last M state state - end of last hit */
  for(z2 = tr->N-1; z1 >= 0; z2--) if(tr->st[z2] == p7T_M) break;

  hit_cnt = 0;
  start_new = TRUE;
  z = z1;

  while(z < z2) {

    if(start_new) {

      /* Save z value - currently set to fist M state in exon */
      y = z;

      /*Find end of exon */
      while(tr->st[z] != p7T_R && tr->st[z] != p7T_E) z++;

     /*Trace back to last M state of exon*/
      while(tr->st[z] != p7T_M) z--;

      /* Get exon coords */
      ihmm = tr->k[y] + hmm_start - 1;
      jhmm = tr->k[z] + hmm_start - 1;

      if(graph->revcomp) {
        iali = gap_seq->n - tr->i[y] + gap_seq->end - seq_i + 2;
        jali = gap_seq->n - tr->i[z] + gap_seq->end - seq_i;
      }
      else {
        iali = gap_seq->start + seq_i + tr->i[y] - 3;
        jali = gap_seq->start + seq_i + tr->i[z] - 1;
      }

      /* check if this exon is a duplicate of an original hit */
      duplicate = FALSE;

      exon_min = ESL_MIN(iali, jali);
      exon_max = ESL_MAX(iali, jali);

      for (h = 0 ; h < graph->orig_N; h++) {
        /* If hit is not in gap skip it */
        if(th->hit[h]->dcl->jhmm < hmm_start || th->hit[h]->dcl->ihmm > hmm_end) continue;

        hit_min = ESL_MIN(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
        hit_max = ESL_MAX(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
        /* Is exon entirley overlaped by hit on the sequence */
        if( hit_min <= exon_min && hit_max >= exon_max) {
          duplicate = TRUE;
          break;
        }
      }

      if(duplicate) { z++; start_new = FALSE; continue; }

      /* Create new hit and  set ihmm and iali coords*/
      new_hit          = p7_hit_Create_empty();

      new_hit->dcl     = p7_domain_Create_empty();
      new_hit->dcl->tr = p7_trace_fs_Create();

      new_hit->dcl->ihmm = ihmm;
      new_hit->dcl->jhmm = jhmm;
      new_hit->dcl->iali = iali;
      new_hit->dcl->jali = jali;
      printf("ihmm %d jhmm %d iali %d jali %d\n", ihmm, jhmm, iali, jali);
      /* Append starting special states */
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_S , 0, 0, 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_N , 0, tr->i[y]-3, 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_B , 0, tr->i[y]-3, 0);

      /* Append all states between first and last M state */
      for(i = y; i <= z; i++) {
        codon = (tr->st[i] == p7T_M ? 3 : 0);
        p7_trace_fs_Append(new_hit->dcl->tr, tr->st[i] , tr->k[i], tr->i[i], codon);
      }

      /* Append ending special states */
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_E, tr->k[z], tr->i[z], 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_C, 0, tr->i[z], 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_T , 0, 0, 0);

      new_hit->dcl->scores_per_pos = NULL;
      p7_splice_ComputeAliScores_fs(new_hit->dcl, new_hit->dcl->tr, gap_seq->dsq+seq_i, sub_fs_model, gap_seq->abc);

      gap_hits[hit_cnt] = new_hit;
      hit_cnt++;

      start_new = FALSE;
    }

    z++;
    if(tr->st[z] == p7T_M) start_new = TRUE;

  }  

  *num_hits = hit_cnt;

  p7_hmm_Destroy(sub_hmm);
  p7_profile_fs_Destroy(sub_fs_model);
  p7_gmx_Destroy(vit_mx);
  p7_trace_fs_Destroy(tr);

  return gap_hits;

  ERROR:
    p7_hmm_Destroy(sub_hmm);
    p7_profile_fs_Destroy(sub_fs_model);
    p7_gmx_Destroy(vit_mx);
    p7_trace_fs_Destroy(tr);
    if(gap_hits != NULL) free(gap_hits);
    return NULL;
  
}

int
bridge_the_gap(SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, P7_HIT **gap_hits, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *gap_seq, int num_hits)
{
  int         i;
  int         up, down;
  int         curr_num_connected;
  int         next_num_connected;
  int         *curr_connected_hits;
  int         *next_connected_hits;
  int         *has_edge;
  int         *node_id;
  SPLICE_EDGE *edge;
  P7_TOPHITS  *th;
  int         status;

  th = graph->th;

  curr_connected_hits = NULL;
  next_connected_hits = NULL;
  ESL_ALLOC(curr_connected_hits, sizeof(int) * num_hits);
  ESL_ALLOC(next_connected_hits, sizeof(int) * num_hits);

  node_id = NULL;
  ESL_ALLOC(node_id, sizeof(int) * num_hits);

  has_edge = NULL;
  ESL_ALLOC(has_edge, sizeof(int) * num_hits);

  curr_num_connected = 0;
  next_num_connected = 0;

  for(i = 0; i < num_hits; i++) {
    curr_connected_hits[i] = FALSE;
    next_connected_hits[i] = FALSE;
    has_edge[i]            = FALSE;
    node_id[i]             = -1;
  }

  edge = NULL;

  /* Connect new hits upstream of original downstream hits */
  for(up = 0; up < num_hits; up++) {
    for(down = 0; down < graph->split_N; down++) { 
      if(!graph->in_graph[down]) continue;
      if(!hit_upstream(gap_hits[up]->dcl, th->hit[down]->dcl, graph->revcomp)) continue;
      if(!hits_spliceable(gap_hits[up]->dcl, th->hit[down]->dcl)) continue;
      printf("up %d down %d \n", up+1, down +1);
      edge = connect_nodes_with_edge(pli, gap_hits[up]->dcl, th->hit[down]->dcl, hmm, gcode, gap_seq, graph->revcomp);
      if(edge != NULL) {

        edge->upstream_node_id   = (has_edge[up] ? node_id[up] : th->N);
        edge->downstream_node_id = down;

        if(!has_edge[up]) {
          node_id[up] = th->N;
          p7_splicegraph_AddNode(graph, gap_hits[up]);
        }
        p7_splicegraph_AddEdge(graph, edge);
        has_edge[up] = TRUE;
        curr_connected_hits[up] = TRUE;
        curr_num_connected++;
      }
    }
  }

  /* Continue connecting new hits upstream until we run out of connections */
  while(curr_num_connected > 0) {

    for(up = 0; up < num_hits; up++) {
      for(down = 0; down < num_hits; down++) {
        if(!curr_connected_hits[down] || up == down) continue;
        if(has_edge[up] && p7_splicegraph_EdgeExists(graph, node_id[up], node_id[down])) continue;
        if(!hit_upstream(gap_hits[up]->dcl, gap_hits[down]->dcl, graph->revcomp)) continue;
        if(!hits_spliceable(gap_hits[up]->dcl, gap_hits[down]->dcl)) continue;

        edge = connect_nodes_with_edge(pli, gap_hits[up]->dcl, gap_hits[down]->dcl, hmm, gcode, gap_seq, graph->revcomp);

        if(edge != NULL) {
          edge->upstream_node_id   = (has_edge[up] ? node_id[up] : th->N);
          edge->downstream_node_id = node_id[down];

          if(!has_edge[up]) {
            node_id[up] = th->N;
            p7_splicegraph_AddNode(graph, gap_hits[up]); 
          }
          p7_splicegraph_AddEdge(graph, edge); 
          has_edge[up] = TRUE;
          next_connected_hits[up] = TRUE;
          next_num_connected++;
        }
      }
    }
   
    free(curr_connected_hits);
    curr_connected_hits = next_connected_hits;
    curr_num_connected  = next_num_connected;
    next_connected_hits = NULL;
    ESL_ALLOC(next_connected_hits, sizeof(int) * num_hits);
    for(i = 0; i < num_hits; i++)
      next_connected_hits[i] = FALSE;
    next_num_connected = 0;
  }


  /* Connect new hits downstream of original gap upstream hit * */
  for(up = 0; up < graph->split_N; up++) {
    for(down = 0; down < num_hits; down++) {
      if(!graph->in_graph[up]) continue;
      if(!hit_upstream(th->hit[up]->dcl, gap_hits[down]->dcl, graph->revcomp)) continue;
      if(!hits_spliceable(th->hit[up]->dcl, gap_hits[down]->dcl)) continue;
  
      edge = connect_nodes_with_edge(pli, th->hit[up]->dcl, gap_hits[down]->dcl, hmm, gcode, gap_seq, graph->revcomp);
      if(edge != NULL) {
  
        edge->upstream_node_id   = up;
        edge->downstream_node_id = (has_edge[down] ? node_id[down] : th->N);;
  
        if(!has_edge[down]) {
          node_id[down] = th->N;
          p7_splicegraph_AddNode(graph, gap_hits[down]);
        }
        p7_splicegraph_AddEdge(graph, edge);
        has_edge[down] = TRUE;
        curr_connected_hits[down] = TRUE;
        curr_num_connected++;
      }
    }
  }
   
  /* Continue connecting non-original hits downstream until we run out of connections */
  while(curr_num_connected > 0) {
    for(up = 0; up < num_hits; up++) {
      for(down = 0; down < num_hits; down++) {
        if(!curr_connected_hits[up] || up == down) continue; 
        if(has_edge[down] && p7_splicegraph_EdgeExists(graph, node_id[up], node_id[down])) continue;
        if(!hit_upstream(gap_hits[up]->dcl, gap_hits[down]->dcl, graph->revcomp)) continue;
        if(!hits_spliceable(gap_hits[up]->dcl, gap_hits[down]->dcl)) continue;

        edge = connect_nodes_with_edge(pli, gap_hits[up]->dcl, gap_hits[down]->dcl, hmm, gcode, gap_seq, graph->revcomp);
        if(edge != NULL) {

          edge->upstream_node_id   = node_id[up];
          edge->downstream_node_id = (has_edge[down] ? node_id[down] : th->N);;

          if(!has_edge[down]) {
            node_id[down] = th->N;
            p7_splicegraph_AddNode(graph, gap_hits[down]);
          }
          p7_splicegraph_AddEdge(graph, edge);         
          has_edge[down] = TRUE;
          curr_connected_hits[down] = TRUE;
          curr_num_connected++;
        }
      }
    }
    free(curr_connected_hits);
    curr_connected_hits = next_connected_hits;
    curr_num_connected  = next_num_connected;
    next_connected_hits = NULL;
    ESL_ALLOC(next_connected_hits, sizeof(int) * num_hits);
    for(i = 0; i < num_hits; i++)
      next_connected_hits[i] = FALSE;
    next_num_connected = 0;
  }

  for (i = 0; i < num_hits; i++) {
    if(!has_edge[i]) {
      p7_trace_fs_Destroy(gap_hits[i]->dcl->tr);
      free(gap_hits[i]->dcl->scores_per_pos);
      p7_hit_Destroy(gap_hits[i]);
    }
  }

  if(curr_connected_hits != NULL) free(curr_connected_hits);
  if(next_connected_hits != NULL) free(next_connected_hits);
  if(has_edge            != NULL) free(has_edge);
  if(node_id             != NULL) free(node_id);

  return eslOK;

  ERROR:
    if(curr_connected_hits != NULL) free(curr_connected_hits);
    if(next_connected_hits != NULL) free(next_connected_hits);
    if(has_edge            != NULL) free(has_edge);
    if(node_id             != NULL) free(node_id);
    return status; 
}
/* Function: get_sub_sequence
 *
 *  Synopsis: Fetch a target subsequence based on <target_range> coords
 *
 *  Purpose:  Fetch the sub-sequence dffined by <target_range> from the
 *            sequence file <seq_file>.  An ssi index for the <seq_file>
 *            must exist. 
 *
 *  Returns:  a pointer to an <ESL_SQ> representing the sub-sequence.
 *
 *  Throws:   <NULL> on allocation failure.
 */
ESL_SQ* 
get_sub_sequence(const ESL_SQFILE *seq_file, char* seqname, int64_t seq_min, int64_t seq_max, int revcomp)
{

  ESL_SQ     *target_seq;
  ESL_SQFILE *tmp_file;

  /* Open seq file and ssi */
  esl_sqfile_Open(seq_file->filename,seq_file->format,NULL,&tmp_file);
  esl_sqfile_OpenSSI(tmp_file,NULL);

  /* Get basic sequence info */
  target_seq = esl_sq_Create();
  esl_sqio_FetchInfo(tmp_file, seqname, target_seq);

  target_seq->start = seq_min;
  target_seq->end   = seq_max;
  target_seq->abc   = seq_file->abc;

  /* Make sure target range coords did not extend too far */
  if (target_seq->start < 1)
    target_seq->start = 1; 

  if (target_seq->end > target_seq->L)
    target_seq->end = target_seq->L;

  /* Fetch target range sequcene */
  if (esl_sqio_FetchSubseq(tmp_file,target_seq->name,target_seq->start,target_seq->end,target_seq) != eslOK) 
    esl_fatal(esl_sqfile_GetErrorBuf(tmp_file));

  esl_sqfile_Close(tmp_file);

  esl_sq_SetName(target_seq, seqname);
 
  if(revcomp)
   esl_sq_ReverseComplement(target_seq);

  esl_sq_Digitize(target_seq->abc, target_seq);

  return target_seq;
}

P7_HMM*
get_sub_hmm (const P7_HMM *hmm, int start, int end)
{

  int i, k;
  int M;
  P7_HMM *sub_hmm;
  int status;

  M = end - start + 1;
  sub_hmm = NULL;
  sub_hmm = p7_hmm_CreateShell();
  sub_hmm->flags = hmm->flags;

  p7_hmm_CreateBody(sub_hmm, M, hmm->abc);

  if ((status = esl_strdup(hmm->name,   -1, &(sub_hmm->name)))   != eslOK) goto ERROR;
  if ((status = esl_strdup(hmm->acc,    -1, &(sub_hmm->acc)))    != eslOK) goto ERROR;
  if ((status = esl_strdup(hmm->desc,   -1, &(sub_hmm->desc)))   != eslOK) goto ERROR;

  sub_hmm->nseq       = hmm->nseq;
  sub_hmm->eff_nseq   = hmm->eff_nseq;
  sub_hmm->max_length = hmm->max_length;
  sub_hmm->checksum   = hmm->checksum;
  sub_hmm->offset     = hmm->offset;
  sub_hmm->fs         = hmm->fs;

  for (i = 0; i < p7_NEVPARAM; i++) sub_hmm->evparam[i] = hmm->evparam[i];
  for (i = 0; i < p7_NCUTOFFS; i++) sub_hmm->cutoff[i]  = hmm->cutoff[i];
  for (i = 0; i < p7_MAXABET;  i++) sub_hmm->compo[i]   = hmm->compo[i];


  /*Copy zeroth positions */
  esl_vec_FCopy(hmm->t[0],   p7H_NTRANSITIONS, sub_hmm->t[0]);
  esl_vec_FCopy(hmm->mat[0], hmm->abc->K,      sub_hmm->mat[0]);
  esl_vec_FCopy(hmm->ins[0], hmm->abc->K,      sub_hmm->ins[0]);

  if (hmm->flags & p7H_RF)    sub_hmm->rf[0]        = hmm->rf[0];
  if (hmm->flags & p7H_MMASK) sub_hmm->mm[0]        = hmm->mm[0];
  if (hmm->flags & p7H_CONS)  sub_hmm->consensus[0] = hmm->consensus[0];
  if (hmm->flags & p7H_CS)    sub_hmm->cs[0]        = hmm->cs[0];

  k = 1;
  for (i = start; i <= end; i++) {
    esl_vec_FCopy(hmm->t[i],   p7H_NTRANSITIONS, sub_hmm->t[k]);
    esl_vec_FCopy(hmm->mat[i], hmm->abc->K,      sub_hmm->mat[k]);
    esl_vec_FCopy(hmm->ins[i], hmm->abc->K,      sub_hmm->ins[k]);
    if (hmm->flags & p7H_RF)    sub_hmm->rf[k]        = hmm->rf[i];
    if (hmm->flags & p7H_MMASK) sub_hmm->mm[k]        = hmm->mm[i];
    if (hmm->flags & p7H_CONS)  sub_hmm->consensus[k] = hmm->consensus[i];
    if (hmm->flags & p7H_CS)    sub_hmm->cs[k]        = hmm->cs[i];
    k++;
  }

  if (hmm->flags & p7H_RF)    sub_hmm->rf[M+1]        = hmm->rf[hmm->M+1];
  if (hmm->flags & p7H_MMASK) sub_hmm->mm[M+1]        = hmm->mm[hmm->M+1];
  if (hmm->flags & p7H_CONS)  sub_hmm->consensus[M+1] = hmm->consensus[hmm->M+1];
  if (hmm->flags & p7H_CS)    sub_hmm->cs[M+1]        = hmm->cs[hmm->M+1];

  /* No delete atart at postion M */
  sub_hmm->t[M][p7H_MD] = 0.0;
  sub_hmm->t[M][p7H_DD] = 0.0;


  return sub_hmm;

  ERROR:
    if(sub_hmm != NULL) p7_hmm_Destroy(sub_hmm);
    return NULL;
}



int
p7_splice_ComputeAliScores(P7_DOMAIN *dom, P7_TRACE *tr, ESL_DSQ *amino_dsq, const P7_PROFILE *gm)
{

  int status;
  int i, k, n;
  int z1, z2;
  int N;

  if(tr->ndom == 0) p7_trace_Index(tr);
  
  N = tr->tto[0] - tr->tfrom[0] - 1;
  if(dom->scores_per_pos == NULL)
    ESL_ALLOC( dom->scores_per_pos, sizeof(float) * N);
  else
    ESL_REALLOC( dom->scores_per_pos, sizeof(float) * N);
  for (n=0; n<N; n++)  dom->scores_per_pos[n] = 0.0;  

  i = tr->sqfrom[0];
  k = tr->hmmfrom[0];
  z1 = tr->tfrom[0]+1;
  z2 = tr->tto[0];
  n = 0;

  while (z1 < z2) {
    if (tr->st[z1] == p7T_M) {
      dom->scores_per_pos[n] = p7P_MSC(gm, k, amino_dsq[i]);
      if (tr->st[z1-1] == p7T_I) 
        dom->scores_per_pos[n] += p7P_TSC(gm, k-1, p7P_IM);
      else if (tr->st[z1-1] == p7T_D)
        dom->scores_per_pos[n] += p7P_TSC(gm, k-1, p7P_DM); 
      i++; k++; z1++; n++;
 
      while(z1 < z2 && tr->st[z1] == p7T_M) {
        dom->scores_per_pos[n] =  p7P_MSC(gm, k,   amino_dsq[i]);
        dom->scores_per_pos[n] += p7P_TSC(gm, k-1, p7P_MM);
        i++; k++; z1++; n++;
      }

    }
    else if (tr->st[z1] == p7T_I) {
      dom->scores_per_pos[n] = p7P_TSC(gm, k, p7P_MI);
      i++; z1++; n++;
      while (z1 < z2 && tr->st[z1] == p7T_I) {
        dom->scores_per_pos[n] = p7P_TSC(gm, k, p7P_II);
        i++; z1++; n++;
      }
     }
     else if (tr->st[z1] == p7T_D) {
       dom->scores_per_pos[n] = p7P_TSC(gm, k-1, p7P_MD);
       k++; z1++; n++;
       while (z1 < z2 && tr->st[z1] == p7T_D)  {
         dom->scores_per_pos[n] = p7P_TSC(gm, k-1, p7P_DD);
         k++; z1++; n++;
       } 
     }
     else ESL_XEXCEPTION(eslFAIL, "Impossible state from compute_ali_scores()");
  }

  dom->aliscore = 0.0;
  for (n=0; n<N; n++) dom->aliscore += dom->scores_per_pos[n];

  return eslOK;

  ERROR:
    return status;
}




int
p7_splice_ComputeAliScores_fs(P7_DOMAIN *dom, P7_TRACE *tr, ESL_DSQ *nuc_dsq, P7_FS_PROFILE *gm_fs, const ESL_ALPHABET *abc)
{

  int status;
  int i, k, c, n;
  int codon_idx;
  int z1, z2;
  int N;
 
  if(tr->ndom == 0) p7_trace_Index(tr);

  N = tr->tto[0] - tr->tfrom[0] - 1;
  if(dom->scores_per_pos == NULL)
    ESL_ALLOC( dom->scores_per_pos, sizeof(float) * N);
  else
    ESL_REALLOC( dom->scores_per_pos, sizeof(float) * N);
  for (n=0; n<N; n++)  dom->scores_per_pos[n] = 0.0; 

  k  = tr->k[tr->tfrom[0]+1];
  z1 = tr->tfrom[0]+1;
  z2 = tr->tto[0];
  n  = 0;

  while (z1<z2) {
    
    if (tr->st[z1] == p7T_M) {
      i = tr->i[z1];
      c = tr->c[z1];
      if(c == 1) {
        if(esl_abc_XIsCanonical(abc, nuc_dsq[i]))
          codon_idx = p7P_CODON1(nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_QC2;
      }
      else if(c == 2) {
        if(esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i]))
          codon_idx = p7P_CODON2(nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_QC1;
      }
      else if(c == 3) {
        if(esl_abc_XIsCanonical(abc, nuc_dsq[i-2]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i]))
          codon_idx = p7P_CODON3(nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_C;
      }  
      else if(c == 4) {
        if(esl_abc_XIsCanonical(abc, nuc_dsq[i-3]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-2]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i]))
          codon_idx = p7P_CODON4(nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_QC1;
      }
      else if(c == 5) {
        if(esl_abc_XIsCanonical(abc, nuc_dsq[i-4]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-3]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-2]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i]))
          codon_idx = p7P_CODON5(nuc_dsq[i-4], nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_QC2;
      }
      
      dom->scores_per_pos[n] = p7P_MSC_CODON(gm_fs, k, codon_idx);
      if (tr->st[z1-1] == p7T_I)
        dom->scores_per_pos[n] += p7P_TSC(gm_fs, k-1, p7P_IM);
      else if (tr->st[z1-1] == p7T_D)
        dom->scores_per_pos[n] += p7P_TSC(gm_fs, k-1, p7P_DM);
      k++; z1++; n++;

      while(z1 < z2 && tr->st[z1] == p7T_M) {
        c = tr->c[z1];
        i = tr->i[z1];       
    
        if(c == 1) {
          if(esl_abc_XIsCanonical(abc, nuc_dsq[i]))
            codon_idx = p7P_CODON1(nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_QC2;
        }
        else if(c == 2) {
          if(esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i]))
            codon_idx = p7P_CODON2(nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_QC1;
        }
        else if(c == 3) {
          if(esl_abc_XIsCanonical(abc, nuc_dsq[i-2]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i]))
            codon_idx = p7P_CODON3(nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_C;
        }
        else if(c == 4) {
          if(esl_abc_XIsCanonical(abc, nuc_dsq[i-3]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i-2]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i]))
            codon_idx = p7P_CODON4(nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_QC1;
        }
        else if(c == 5) {
          if(esl_abc_XIsCanonical(abc, nuc_dsq[i-4]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i-3]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i-2]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
           esl_abc_XIsCanonical(abc, nuc_dsq[i]))
            codon_idx = p7P_CODON5(nuc_dsq[i-4], nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_QC2;
        }
        
        dom->scores_per_pos[n] = p7P_MSC_CODON(gm_fs, k, codon_idx);
        dom->scores_per_pos[n] += p7P_TSC(gm_fs, k-1, p7P_MM);
        k++; z1++; n++;
      }
    }
    else if (tr->st[z1] == p7T_I) {
      dom->scores_per_pos[n] = p7P_TSC(gm_fs, k, p7P_MI);
      z1++; n++;
      while (z1 < z2 && tr->st[z1] == p7T_I) {
        dom->scores_per_pos[n] = p7P_TSC(gm_fs, k, p7P_II);
        z1++; n++;
      }
    }
    else if (tr->st[z1] == p7T_D) {
      dom->scores_per_pos[n] = p7P_TSC(gm_fs, k-1, p7P_MD);
      k++; z1++; n++;
      while (z1 < z2 && tr->st[z1] == p7T_D)  {
        dom->scores_per_pos[n] = p7P_TSC(gm_fs, k-1, p7P_DD);
        k++; z1++; n++;
      }
    }
    else ESL_XEXCEPTION(eslFAIL, "Impossible state from compute_ali_scores()");
  }

  dom->aliscore = 0.0;
  for (n=0; n<N; n++)  dom->aliscore += dom->scores_per_pos[n];  
 
  return eslOK;

  ERROR:
    return status;
}



