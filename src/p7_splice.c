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
static P7_HIT** split_hit(P7_HIT *hit, const P7_HMM *hmm, const P7_BG *bg, ESL_SQ *hit_seq, const ESL_GENCODE *gcode, int revcomp, int *num_hits);
static int hit_upstream(P7_DOMAIN *upstream, P7_DOMAIN *downstream, int revcomp);
static int hits_spliceable(P7_DOMAIN *upstream, P7_DOMAIN *downstream);
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

      printf("\nQuery %s Target %s strand %c\n", gm->name, graph->seqname, (graph->revcomp ? '-' : '+'));
      fflush(stdout);

      p7_splicegraph_DumpHits(stdout, graph);

      if((p7_splice_SplitHits(graph, hmm, pli->bg, gcode, seq_file)) != eslOK) goto ERROR;
      p7_splicegraph_DumpHits(stdout, graph);        


      p7_splice_RecoverHits(graph, saved_hits, pli, hmm, gcode, seq_file, first, last);     
      p7_splicegraph_DumpHits(stdout, graph);
      

      p7_splice_FillGaps(graph, pli, hmm, gcode, seq_file);
      p7_splicegraph_DumpHits(stdout, graph);

      p7_splice_ConnectGraph(graph, pli, hmm, gcode, seq_file);
      p7_splicegraph_DumpEdges(stdout, graph, FALSE);

      p7_splice_RemoveDisconnected(graph, hits_processed, &num_hits_processed);  
      //use aliscore computations to set the fs and stop codon cnts in the trace.
    }

    //path = p7_splicepath_GetBestPath(graph);   
    //enforce_range_bounds(graph, target_range->th, range_bound_mins, range_bound_maxs, range_cnt);

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
p7_splice_SplitHits(SPLICE_GRAPH *graph, const P7_HMM *hmm, const P7_BG *bg, ESL_GENCODE *gcode, const ESL_SQFILE *seq_file) 
{

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

  sub_hmm     = p7_splice_GetSubHMM(hmm, hit->dcl->ihmm, hit->dcl->jhmm);
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

      edge = p7_spliceedge_ConnectHits(pli, th->hit[up]->dcl, th->hit[down]->dcl, hmm, gcode, splice_seq, graph->revcomp);
 
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
p7_splice_RemoveDisconnected(SPLICE_GRAPH *graph, int *hits_processed, int *num_hits_processed) 
{

  int i, j;
  int connected;
  int split_in_graph;
  P7_TOPHITS *th;

  th = graph->th;

  /* If any hit is not connected to a different original or split hit, set as not in graph */
  for(i = 0; i < th->N; i++) {
    if(!graph->in_graph[i]) continue;
    connected = FALSE;
    for(j = 0; j < graph->split_N; j++) {
      if(!graph->in_graph[j] || i == j) continue;

      if(hit_upstream(th->hit[i]->dcl, th->hit[j]->dcl, graph->revcomp)) {
        if(p7_splicegraph_PathExists(graph,i,j)) {
          connected = TRUE;
          break;
        }
      }
      else if(hit_upstream(th->hit[j]->dcl, th->hit[i]->dcl, graph->revcomp)) {
        if(p7_splicegraph_PathExists(graph,j,i)) {
          connected = TRUE;
          break;
        }
      }
    }
    if(!connected) 
      graph->in_graph[i] = FALSE;
  }

  /* For each original or split hit not in graph, check if any split hits with 
   * the same orig_hit_idx are still in graph and if not set hit as processed */
  for(i = 0; i < graph->split_N; i++) {
    if(graph->in_graph[i] || hits_processed[graph->orig_hit_idx[i]]) continue;
    split_in_graph = FALSE;
    for(j = graph->orig_N; j < graph->split_N; j++) {
      if(graph->orig_hit_idx[i] == graph->orig_hit_idx[j] && graph->in_graph[j]) {
        split_in_graph = TRUE;
        break;
      } 
    }
    if(!split_in_graph) {
      hits_processed[graph->orig_hit_idx[i]] = TRUE;
      *num_hits_processed++;
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

int
p7_splice_RecoverHits(SPLICE_GRAPH *graph, SPLICE_SAVED_HITS *sh, SPLICE_PIPELINE *pli, P7_HMM *hmm, ESL_GENCODE *gcode, ESL_SQFILE *seq_file, int first, int last) {

  int i, j;
  int z1, z2;
  int seq_min, seq_max;
  P7_DOMAIN     *tmp_dom;
  P7_HIT        *new_hit;
  P7_TOPHITS    *th;
  P7_HMM        *sub_hmm;
  P7_FS_PROFILE *sub_fs_model;
  P7_GMX        *vit_mx;
  ESL_SQ        *hit_seq; 

  th = graph->th;
  graph->recover_N = graph->split_N;

  p7_splicehits_FindNodes(graph, sh, first, last);
  p7_splicehits_Dump(stdout, sh);

  tmp_dom = p7_domain_Create_empty();
  for(i = first; i <= last; i++) {
    
    if(sh->srt[i]->duplicate) continue;
    if(sh->srt[i]->node_id >= 0) continue;
  
    tmp_dom->ihmm = sh->srt[i]->hmm_start;
    tmp_dom->jhmm = sh->srt[i]->hmm_end;    
    tmp_dom->iali = sh->srt[i]->seq_start;    
    tmp_dom->jali = sh->srt[i]->seq_end;

    /* Check if hit info is for a hit that is spliceable upstream or downstream of any hits int the graph*/
    for(j = 0; j < th->N; j++) {
      if(!graph->in_graph[j]) continue;
      if(!((hit_upstream(tmp_dom, th->hit[j]->dcl, graph->revcomp) && hits_spliceable(tmp_dom, th->hit[j]->dcl)) ||
           (hit_upstream(th->hit[j]->dcl, tmp_dom, graph->revcomp) && hits_spliceable(th->hit[j]->dcl, tmp_dom)))) continue;

      new_hit          = p7_hit_Create_empty();
      new_hit->dcl     = p7_domain_Create_empty();
      new_hit->dcl->tr = p7_trace_fs_Create();

      seq_min = ESL_MIN(tmp_dom->iali, tmp_dom->jali) - 3;
      seq_max = ESL_MAX(tmp_dom->iali, tmp_dom->jali) + 3; 
      hit_seq = get_sub_sequence(seq_file, graph->seqname, seq_min, seq_max, graph->revcomp);

      sub_hmm = p7_splice_GetSubHMM(hmm, tmp_dom->ihmm-1, tmp_dom->jhmm+1);
      sub_hmm->fs = 0.;
      
      sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
      p7_ProfileConfig_fs(sub_hmm, pli->bg, gcode, sub_fs_model, hit_seq->n*3, p7_UNILOCAL); //hit_seq->n*3 to get single nucleotide loops    
            
      vit_mx = p7_gmx_fs_Create(sub_fs_model->M, hit_seq->n, hit_seq->n, 0);

      p7_trans_Viterbi(hit_seq->dsq, gcode, hit_seq->n, sub_fs_model, vit_mx, NULL); 
      p7_trans_VTrace(hit_seq->dsq, hit_seq->n, gcode, sub_fs_model, vit_mx, new_hit->dcl->tr);

      for(z1 = 0; z1 < new_hit->dcl->tr->N; z1++)    if(new_hit->dcl->tr->st[z1] == p7T_M) break;
      for(z2 = new_hit->dcl->tr->N-1; z2 >= 0; z2--) if(new_hit->dcl->tr->st[z2] == p7T_M) break;
      
      new_hit->dcl->ihmm = new_hit->dcl->tr->k[z1] + tmp_dom->ihmm - 2;
      new_hit->dcl->jhmm = new_hit->dcl->tr->k[z2] + tmp_dom->ihmm - 2; 
      new_hit->dcl->iali = new_hit->dcl->tr->i[z1] + tmp_dom->iali - 6;
      new_hit->dcl->jali = new_hit->dcl->tr->i[z2] + tmp_dom->jali - 4;

      p7_splice_ComputeAliScores_fs(new_hit->dcl, new_hit->dcl->tr, hit_seq->dsq, sub_fs_model, hit_seq->abc);      

      p7_splicegraph_AddNode(graph, new_hit);

      graph->recover_N++;

      p7_hmm_Destroy(sub_hmm);
      p7_profile_fs_Destroy(sub_fs_model);
      p7_gmx_Destroy(vit_mx);
      esl_sq_Destroy(hit_seq);
      break;
    } 
    
  }     
    
  free(tmp_dom);
  return eslOK;
}

int
p7_splice_FillGaps(SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQFILE *seq_file)
{

  int up, down;
  int i, j, h;
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
   
    for(h = 0; h < hit_cnt; h++) 
      p7_splicegraph_AddNode(graph, gap_hits[h]); 
    //if(hit_cnt) bridge_the_gap(graph, pli, gap_hits, hmm, gcode, gap_seq, hit_cnt);

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
  int new_hit_min, new_hit_max;
  int old_hit_min, old_hit_max;
  int new_hit_len, old_hit_len;
  int hmm_overlap_min, hmm_overlap_max;
  int seq_overlap_min, seq_overlap_max;
  int hmm_overlap_len, seq_overlap_len;
  P7_TOPHITS    *th;
  P7_HMM        *sub_hmm;
  P7_FS_PROFILE *sub_fs_model;
  P7_GMX        *vit_mx;
  P7_TRACE      *tr;
  P7_HIT      **gap_hits;
  P7_HIT       *new_hit; 
  int           status;

  th = graph->th;

  sub_hmm = p7_splice_GetSubHMM(hmm, hmm_start, hmm_end);
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

      /* check if this new hit is a duplicate */
      duplicate = FALSE;

      new_hit_min = ESL_MIN(iali, jali);
      new_hit_max = ESL_MAX(iali, jali);
      new_hit_len = new_hit_max - new_hit_min + 1;

      for (h = 0 ; h < th->N; h++) {
        /* If hit is not in gap skip it */
        if(th->hit[h]->dcl->jhmm < hmm_start || th->hit[h]->dcl->ihmm > hmm_end) continue;

        hmm_overlap_min = ESL_MAX(ihmm, th->hit[h]->dcl->ihmm);
        hmm_overlap_max = ESL_MIN(jhmm, th->hit[h]->dcl->jhmm);
        hmm_overlap_len = hmm_overlap_max - hmm_overlap_min + 1;

        if(hmm_overlap_len < 1) continue;

        old_hit_min = ESL_MIN(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
        old_hit_max = ESL_MAX(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
        old_hit_len = old_hit_max - old_hit_min + 1;

        seq_overlap_min = ESL_MAX(new_hit_min, old_hit_min);
        seq_overlap_max = ESL_MIN(new_hit_max, old_hit_max);
        seq_overlap_len = seq_overlap_max - seq_overlap_min + 1; 

        if(seq_overlap_len >= new_hit_len * 0.95 || 
           seq_overlap_len >= old_hit_len * 0.95) {
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
      edge = p7_spliceedge_ConnectHits(pli, gap_hits[up]->dcl, th->hit[down]->dcl, hmm, gcode, gap_seq, graph->revcomp);
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

        edge = p7_spliceedge_ConnectHits(pli, gap_hits[up]->dcl, gap_hits[down]->dcl, hmm, gcode, gap_seq, graph->revcomp);

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
  
      edge = p7_spliceedge_ConnectHits(pli, th->hit[up]->dcl, gap_hits[down]->dcl, hmm, gcode, gap_seq, graph->revcomp);
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

        edge = p7_spliceedge_ConnectHits(pli, gap_hits[up]->dcl, gap_hits[down]->dcl, hmm, gcode, gap_seq, graph->revcomp);
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
p7_splice_GetSubHMM (const P7_HMM *hmm, int start, int end)
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



