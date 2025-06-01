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

static P7_HIT** split_hit(P7_HIT *hit, P7_HMM *sub_hmm, const P7_BG *bg, ESL_SQ *hit_seq, const ESL_GENCODE *gcode, int revcomp, int *num_hits);
static SPLICE_EDGE** connect_split(SPLICE_PIPELINE *pli, P7_HIT **split_hits, const P7_HMM *hmm, P7_FS_PROFILE *gm_fs, ESL_GENCODE *gcode, ESL_SQ *hit_seq, int revcomp, int frameshift, int *num_hits);
static int align_spliced_path (SPLICE_PIPELINE *pli, P7_OPROFILE *om, P7_PROFILE *gm, ESL_SQ *target_seq, ESL_GENCODE *gcode, float fs_prob);
static int align_spliced_path_frameshift (SPLICE_PIPELINE *pli, P7_FS_PROFILE *gm_fs, ESL_SQ *target_seq, ESL_GENCODE *gcode);
static int hit_upstream(P7_DOMAIN *upstream, P7_DOMAIN *downstream, int revcomp);
static int hits_spliceable_short(P7_DOMAIN *upstream, P7_DOMAIN *downstream);
static int hits_spliceable_long(P7_DOMAIN *upstream, P7_DOMAIN *downstream);

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
p7_splice_SpliceHits(P7_TOPHITS *tophits, SPLICE_SAVED_HITS *saved_hits, P7_HMM *hmm, P7_OPROFILE *om, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, P7_SCOREDATA *scoredata, ESL_GETOPTS *go, ESL_GENCODE *gcode, ESL_SQFILE *seq_file, FILE *ofp, int64_t db_nuc_cnt, int fs_pipe, int std_pipe)
{

  int              i, p;
  int              hit_cnt;
  int              success;
  int              hit_to_process;
  int              num_hits_processed; 
  int              prev_num_hits_processed;
  int              is_sorted_by_sortkey;
  int              revcomp, prev_revcomp;
  int              frameshift;
  int              first, last;
  int              seq_min, seq_max;
  int              num_paths;
  int64_t          seqidx, prev_seqidx;
  int64_t          bound_min, bound_max;
  int             *hits_processed;
  SPLICE_PIPELINE *pli;
  SPLICE_GRAPH    *graph;
  SPLICE_PATH     *path;
  SPLICE_PATH     **path_accumulator;
  P7_HIT          *seed_hit; 
  ESL_SQ          *path_seq;
  ESL_SQ          *ali_seq;
  int              status;

  seed_hit         = NULL;
  graph            = NULL;
  path             = NULL;
  hits_processed   = NULL;
  path_accumulator = NULL;

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

  /* Check for any hits that are duplicates or have a p-value >= 1 and mark them as processed */
  hit_cnt = 0;
  num_hits_processed = 0;
  for(i = 0; i < tophits->N; i++) {
    if ((tophits->hit[i]->flags & p7_IS_DUPLICATE)) {
      hits_processed[i] = 1;
      num_hits_processed++;
    }
    else if(!(tophits->hit[i]->flags & p7_IS_REPORTED) && exp(tophits->hit[i]->sum_lnP) >= 1.0) {
      hits_processed[i] = 1;
      num_hits_processed++;
    }
    else { 
      hits_processed[i] = 0;
      hit_cnt++;
    }  
  }
  
  prev_num_hits_processed = -1;

  prev_revcomp = -1;
  prev_seqidx  = -1;
  printf("\nQuery %s LENG %d \n",  gm->name, gm->M);
   fflush(stdout);

  /* loop through until all hits have been processed */
  hit_to_process = tophits->N; 
  
  /*find the firt unprocessed hit */
  if(num_hits_processed != hit_to_process) {
    i = 0;
    while (i < tophits->N && hits_processed[i]) i++;

    seed_hit = tophits->hit[i];
    seqidx = seed_hit->seqidx;
    revcomp = 1;
    if (seed_hit->dcl->iali < seed_hit->dcl->jali)
      revcomp = 0;
  }

  while(num_hits_processed < hit_to_process) {
    
    if(prev_num_hits_processed == num_hits_processed)
      ESL_XEXCEPTION(eslFAIL, "p7_splice_SpliceHits : loop failed to process hits");
    prev_num_hits_processed = num_hits_processed;

    /* If we have a new sequence and strand we build a new graph */
    if(seqidx != prev_seqidx || revcomp != prev_revcomp) {   
      /* Find the first and last index of the current seqidx and strand in the saved hits list */
      first = 0;
      while(first < saved_hits->N && (saved_hits->srt[first]->seqidx != seqidx || saved_hits->srt[first]->strand != revcomp)) first++;
      last = first;
      while(last < saved_hits->N && saved_hits->srt[last]->seqidx == seqidx && saved_hits->srt[last]->strand == revcomp) last++;
      last--;

      if ((graph = p7_splicegraph_Create()) == NULL) goto ERROR;
      graph->revcomp = revcomp;

      if((p7_splice_AddOriginals(graph, tophits, hits_processed, seqidx)) != eslOK) goto ERROR;      

      printf("\nQuery %s Target %s strand %c\n", gm->name, graph->seqname, (graph->revcomp ? '-' : '+'));
      fflush(stdout);

printf("ORIGINAL\n");
p7_splicegraph_DumpHits(stdout, graph);
//fflush(stdout);

     if((p7_splice_SplitHits(graph, pli, hmm, gm_fs, gcode, seq_file, &hit_to_process)) != eslOK) goto ERROR;

printf("SPLIT\n");
p7_splicegraph_DumpHits(stdout, graph);        
//fflush(stdout);

      p7_splice_RecoverHits(graph, saved_hits, pli, hmm, gcode, seq_file, first, last);     
    printf("RECOVER\n");
   p7_splicegraph_DumpHits(stdout, graph);
//    fflush(stdout);  

      p7_splice_ConnectGraph(graph, pli, hmm, gcode, seq_file);
      p7_splice_FillGaps(graph, pli, hmm, gcode, seq_file);

    printf("GAPS\n");
    p7_splicegraph_DumpHits(stdout, graph);
//    fflush(stdout);

      p7_splice_ConnectGraph(graph, pli, hmm, gcode, seq_file);

     //p7_splicegraph_DumpGraph(stdout, graph, FALSE);
 //    p7_splicegraph_DumpEdges(stdout, graph);
      if(esl_vec_ISum(graph->node_in_graph, graph->num_nodes) == 0) {
        prev_seqidx = seqidx;
        prev_revcomp = revcomp;
        continue;
      }

      /* Since every path must inculde an original hit, orig_N is alos the maximum number of possible paths */
      ESL_ALLOC(path_accumulator, sizeof(SPLICE_PATH*) * graph->orig_N);
      num_paths = 0;
    }
    
    /* Find the highest scoring path in graph */
    path = p7_splicepath_GetBestPath(graph);   
    p7_splicepath_Dump(stdout, path);
    //p7_splicegraph_DumpGraph(stdout, graph, FALSE);

    /* Add path to path_accumulator */
    path_accumulator[num_paths] = path;
    num_paths++;

    /* Mark hits in path bounds as processed and break any edges that cross path bounds */
    bound_min = ESL_MIN(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]);
    bound_max = ESL_MAX(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]);

    p7_splice_ReleaseHits(graph, hits_processed, &num_hits_processed, bound_min, bound_max);
    p7_splice_EnforceRangeBounds(graph, bound_min, bound_max);
    
    /*Check if next seed hit is on the current sequence and strand. Sorting ensures that the 
     * seed hit will be on current sequence and strand until all those hits have been process */
    prev_seqidx = seqidx;
    prev_revcomp = revcomp;

    if(num_hits_processed < hit_to_process) {
      i = 0;
      while (i < tophits->N && hits_processed[i]) i++;

      seed_hit = tophits->hit[i];
      seqidx = seed_hit->seqidx;
      revcomp = 1;
      if (seed_hit->dcl->iali < seed_hit->dcl->jali)
        revcomp = 0; 
    } 

    /* If all hits are processed trigger aligmnent for final path set */
    else prev_seqidx = -1;  

    /* If we have processed all the hits on the previous sequence and strand, align the paths */
    if(seqidx != prev_seqidx || revcomp != prev_revcomp) {   
      p7_splice_MergePaths(graph, path_accumulator, pli, hmm, gcode, seq_file, &num_paths);

      for(p = 0; p < num_paths; p++) {
 
        path = path_accumulator[p];
 //     p7_splicepath_Dump(stdout, path);  
        if(path->path_len > 1) {
          seq_min = ESL_MIN(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]);
          seq_max = ESL_MAX(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]);
        
          path_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, path->revcomp);
        
          frameshift = FALSE;      
          if(!std_pipe) frameshift = TRUE;
          for(i = 0; i < path->path_len; i++) {
            if(path->hits[i]->dcl->tr->fs != 0) frameshift = TRUE;
          }
        
          if(!frameshift)
            p7_splice_AlignPath(graph, path, pli, tophits, om, gm, gcode, path_seq, db_nuc_cnt, gm_fs->fs, &frameshift, &success);
         
          if(frameshift)
            p7_splice_AlignFrameshiftPath (graph, path, pli, tophits, gm_fs, gcode, path_seq, db_nuc_cnt, &success);
           
          esl_sq_Destroy(path_seq);
        }
        else success = FALSE; 
     
        p7_splicepipeline_Reuse(pli);
        p7_splicepath_Destroy(path);
      }
      free(path_accumulator);
      path_accumulator = NULL;
      p7_splicegraph_Destroy(graph);
      graph = NULL; 
    }
  }
 
  /* Build any missing alignments */
  for ( i = 0; i < tophits->N; i++) {
    if(tophits->hit[i]->dcl->ad != NULL) continue;
    if(tophits->hit[i]->flags & p7_IS_REPORTED ) {
      if(tophits->hit[i]->dcl->iali < tophits->hit[i]->dcl->jali) {
        seq_min = tophits->hit[i]->subseq_start;
        seq_max = tophits->hit[i]->dcl->jali;
        revcomp = 0;
      }
      else {
        seq_min = tophits->hit[i]->dcl->jali;
        seq_max = tophits->hit[i]->subseq_start;
        revcomp = 1;
      }
      ali_seq = p7_splice_GetSubSequence(seq_file, tophits->hit[i]->name, seq_min, seq_max, revcomp);
      tophits->hit[i]->dcl->ad = p7_alidisplay_fs_Create(tophits->hit[i]->dcl->tr, 0, gm_fs, ali_seq, gcode);
      tophits->hit[i]->dcl->ad->exon_cnt = 1;
      tophits->hit[i]->dcl->ad->sqfrom = tophits->hit[i]->dcl->iali;
      tophits->hit[i]->dcl->ad->sqto   = tophits->hit[i]->dcl->jali;  
      esl_sq_Destroy(ali_seq);
    }
  }
  /* Leave only footprints */
  if (is_sorted_by_sortkey)
    if ((status = p7_tophits_SortBySortkey(tophits)) != eslOK) goto ERROR;   

  if(pli   != NULL) p7_splicepipeline_Destroy(pli);
  if(graph != NULL) p7_splicegraph_Destroy(graph);

  //if (path_accumulator != NULL) free(path_accumulator);  
  if (hits_processed   != NULL) free(hits_processed);

  return status;

  ERROR:
    if(pli   != NULL) p7_splicepipeline_Destroy(pli);
    if(graph != NULL) p7_splicegraph_Destroy(graph);

    if (hits_processed   != NULL) free(hits_processed);
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
  
  /*Add all hits from current sequence and strand to grapgh*/
  for (i = 0; i < tophits->N; i++) {

    curr_hit = tophits->hit[i];

    if (hits_processed[i]) continue;

    if (curr_hit->seqidx != seqidx) continue;

    if (graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
    if ((!graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;
 
    if(curr_hit->flags & p7_IS_REPORTED) graph->reportable[graph->th->N] = TRUE;
    else                                 graph->reportable[graph->th->N] = FALSE;

    graph->node_in_graph[graph->th->N]     = TRUE;
    graph->orig_hit_idx[graph->th->N] = i;
    
    p7_splicegraph_AddNode(graph, curr_hit); 
  }

  graph->orig_N  = graph->th->N;
  graph->seqidx  = seqidx;
  graph->seqname = graph->th->hit[0]->name;

  return eslOK;
 
}

int
p7_splice_SplitHits(SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, const P7_HMM *hmm, P7_FS_PROFILE *gm_fs, ESL_GENCODE *gcode, const ESL_SQFILE *seq_file, int *hits_to_process) 
{

  int h, s, z;
  int z1, z2;
  int I_cnt;
  int low_pp_cnt;
  int splitable;
  int seq_min, seq_max;
  int num_hits;
  float tmp_fs;
  P7_HIT      *curr_hit;
  P7_TRACE    *tr;
  P7_TOPHITS  *th;
  P7_HMM      *sub_hmm;
  P7_HIT      **split_hits; 
  SPLICE_EDGE **split_edges;
  ESL_SQ      *hit_seq;

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
      
      split_hits = NULL;
      split_edges = NULL;
    
      seq_min = ESL_MIN(curr_hit->dcl->iali, curr_hit->dcl->jali);
      seq_max = ESL_MAX(curr_hit->dcl->iali, curr_hit->dcl->jali);
      
      hit_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, graph->revcomp);      

      sub_hmm     = p7_splice_GetSubHMM(hmm, curr_hit->dcl->ihmm, curr_hit->dcl->jhmm);
      tmp_fs = sub_hmm->fs;
      sub_hmm->fs = 0.;    
 
      split_hits = split_hit(curr_hit, sub_hmm, pli->bg, hit_seq, gcode, graph->revcomp, &num_hits);      
      
      if(num_hits > 1)
        split_edges = connect_split(pli, split_hits, hmm, gm_fs, gcode, hit_seq, graph->revcomp, FALSE, &num_hits); 
      
      if( num_hits == 1) {
        p7_trace_fs_Destroy(split_hits[0]->dcl->tr);
        free(split_hits[0]->dcl->scores_per_pos);
        p7_hit_Destroy(split_hits[0]); 
        if(split_hits != NULL) free(split_hits);
      }
      else if(num_hits == 0) {
        if(split_hits != NULL)  free(split_hits);    
      }      

      if (num_hits > 1) {
         
        /* If the hit was split, set the original hit as not in graph and add new split hits to graph */ 
        graph->node_in_graph[h] = FALSE;
        *hits_to_process += num_hits-1;   
   
        /* Add first split hit */
        p7_splicegraph_Grow(graph);
        graph->node_in_graph[graph->th->N]     = TRUE;
        graph->orig_hit_idx[graph->th->N] = graph->orig_hit_idx[h];
        p7_splicegraph_AddNode(graph, split_hits[0]);
        graph->split_orig_id[graph->split_N] = h;
        graph->split_N++; 
  
        /*Add all other splits hits and edges */
        for(s = 1; s < num_hits; s++) {
          p7_splicegraph_Grow(graph);
          graph->node_in_graph[graph->th->N]     = TRUE;
          graph->orig_hit_idx[graph->th->N] = graph->orig_hit_idx[h];
  
          split_edges[s-1]->upstream_node_id   = graph->th->N-1;
          split_edges[s-1]->downstream_node_id = graph->th->N; 
          p7_splicegraph_AddNode(graph, split_hits[s]);
          p7_splicegraph_AddEdge(graph, split_edges[s-1]); 
          graph->split_orig_id[graph->split_N] = h; 
          graph->split_N++; 
        } 
 
        free(split_hits); 
        free(split_edges);
        split_hits = NULL;
        split_edges = NULL;
      }

      /*If the hit has frameshifts we will split again with framshifts enabled */
      if(curr_hit->dcl->tr->fs) {
        sub_hmm->fs = tmp_fs;
        num_hits = 0;
 
        split_hits = split_hit(curr_hit, sub_hmm, pli->bg, hit_seq, gcode, graph->revcomp, &num_hits);      
         
        if(num_hits > 1)
          split_edges = connect_split(pli, split_hits, hmm, gm_fs, gcode, hit_seq, graph->revcomp, TRUE, &num_hits); 
  

        if(num_hits == 1) { 
          
          p7_trace_fs_Destroy(split_hits[0]->dcl->tr);
          free(split_hits[0]->dcl->scores_per_pos);
          p7_hit_Destroy(split_hits[0]);
          if(split_hits != NULL) free(split_hits);
        }
        else if(num_hits == 0) {
          if(split_hits != NULL)  free(split_hits);
        }

        if (num_hits > 1) {
          
          /* If the hit was split, set the original hit as not in graph and add new split hits to graph */ 
          graph->node_in_graph[h] = FALSE;
          *hits_to_process += num_hits-1;   
     
          /* Add first split hit */
          p7_splicegraph_Grow(graph);
          graph->node_in_graph[graph->th->N]     = TRUE;
          graph->orig_hit_idx[graph->th->N] = graph->orig_hit_idx[h];
          p7_splicegraph_AddNode(graph, split_hits[0]);
          graph->split_orig_id[graph->split_N] = h;
          graph->split_N++; 
    
          /*Add all other splits hits and edges */
          for(s = 1; s < num_hits; s++) {
            p7_splicegraph_Grow(graph);
            graph->node_in_graph[graph->th->N]     = TRUE;
            graph->orig_hit_idx[graph->th->N] = graph->orig_hit_idx[h];
    
            split_edges[s-1]->upstream_node_id   = graph->th->N-1;
            split_edges[s-1]->downstream_node_id = graph->th->N; 
            p7_splicegraph_AddNode(graph, split_hits[s]);
            p7_splicegraph_AddEdge(graph, split_edges[s-1]); 
            graph->split_orig_id[graph->split_N] = h; 
            graph->split_N++; 
          } 
   
          free(split_hits); 
          free(split_edges);
          split_hits = NULL;
          split_edges = NULL;
        }
      }

      esl_sq_Destroy(hit_seq);           
      p7_hmm_Destroy(sub_hmm);
    }
  }

  return eslOK;

}



P7_HIT**
split_hit(P7_HIT *hit, P7_HMM *sub_hmm, const P7_BG *bg, ESL_SQ *hit_seq, const ESL_GENCODE *gcode, int revcomp, int *num_hits)
{
  int         i, y, z;
  int         z1, z2;
  int         intron_cnt;
  int         hit_cnt;
  int         start_new;
  int         ihmm, jhmm;
  int         iali, jali;
  P7_FS_PROFILE *sub_fs_model;
  P7_GMX       *vit_mx;
  P7_HIT       *new_hit;
  P7_HIT      **ret_hits; 
  P7_TRACE     *tr; 
  int status;

  sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
  p7_ProfileConfig_fs(sub_hmm, bg, gcode, sub_fs_model, hit_seq->n, p7_GLOBAL); 
  
  vit_mx = p7_gmx_fs_Create(sub_fs_model->M, hit_seq->n, hit_seq->n, p7P_SPLICE);
  tr = p7_trace_fs_Create();
  
  p7_sp_fs_semiglobal_Viterbi(hit_seq->dsq, gcode, hit_seq->n, sub_fs_model, vit_mx);
  p7_sp_fs_semiglobal_VTrace(hit_seq->dsq, hit_seq->n, gcode, sub_fs_model, vit_mx, tr); 
  
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
      for(i = y; i <= z; i++) 
        p7_trace_fs_Append(new_hit->dcl->tr, tr->st[i], tr->k[i], tr->i[i], tr->c[i]);

      /* Append ending special states */
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_E, tr->k[z], tr->i[z], 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_C, 0, tr->i[z], 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_T, 0, 0, 0);
      
      new_hit->dcl->ad             = NULL; 
      new_hit->dcl->scores_per_pos = NULL;
      p7_splice_ComputeAliScores_fs(new_hit->dcl, new_hit->dcl->tr, hit_seq->dsq, sub_fs_model, hit_seq->abc);

      /* Adjust k postions to full HMM */
      for(i = 0; i < new_hit->dcl->tr->N; i++) {
        if(new_hit->dcl->tr->k[i] != 0)
          new_hit->dcl->tr->k[i] += hit->dcl->ihmm - 1;
      }

      ret_hits[hit_cnt] = new_hit;
      hit_cnt++;

      
      start_new = FALSE;
    }
 
    z++;
    if(tr->st[z] == p7T_M) start_new = TRUE;     
     
  }
 
  *num_hits = hit_cnt;
  
  p7_profile_fs_Destroy(sub_fs_model);
  p7_gmx_Destroy(vit_mx);
  p7_trace_fs_Destroy(tr);

  return ret_hits;


  ERROR:
   if(ret_hits != NULL) free(ret_hits);
   return NULL;
 
}

SPLICE_EDGE**
connect_split(SPLICE_PIPELINE *pli, P7_HIT **split_hits, const P7_HMM *hmm, P7_FS_PROFILE *gm_fs, ESL_GENCODE *gcode, ESL_SQ *hit_seq, int revcomp, int frameshift, int *num_hits) 
{

  int s, x;
  int hit_cnt;
  int edge_found;
  int seq_i, seq_j;
  int seq_len;
  int hmm_i, hmm_j;
  int tmp_ihmm, tmp_iali;
  float *tp_0;
  float *tp_M;
  P7_HMM        *sub_hmm; 
  P7_FS_PROFILE *sub_fs_model;
  P7_GMX        *vit_mx; 
  P7_TRACE      *tr;
  SPLICE_EDGE   *edge;
  SPLICE_EDGE   **split_edges;
  int status;

  hit_cnt = *num_hits;

  split_edges = NULL;
  ESL_ALLOC(split_edges, sizeof(SPLICE_EDGE*) * (hit_cnt-1));
  
  /* Splice all split hits */
  for(s = 1; s < hit_cnt; s++) {
    edge = p7_spliceedge_ConnectSplits(pli, split_hits[s-1]->dcl, split_hits[s]->dcl, hmm, gcode, hit_seq, frameshift, revcomp);        
    split_edges[s-1] = edge;
  }
   
  /* Check if any edges were found */
  edge_found = FALSE;
  for(s = 0; s < hit_cnt-1; s++) 
    if(split_edges[s] != NULL) edge_found = TRUE;
      
  /* If no edges found return null */
  if(!edge_found) {
    for(s = 0; s < hit_cnt; s++) {
      p7_trace_fs_Destroy(split_hits[s]->dcl->tr);
      free(split_hits[s]->dcl->scores_per_pos);
      p7_hit_Destroy(split_hits[s]);
    }
    *num_hits = 0;
    free(split_edges);
    return NULL;
  }
  
  
  /* If we have at least one edge, check if any edge are missing. 
   * If so realign the split hits */
  s = 0;
  while(s < hit_cnt-1) {
    if(split_edges[s] == NULL) {
    
      sub_hmm = p7_splice_GetSubHMM(hmm, split_hits[s]->dcl->ihmm, split_hits[s+1]->dcl->jhmm);  
      seq_i = split_hits[s]->dcl->tr->sqfrom[0];
      seq_j = split_hits[s+1]->dcl->tr->sqto[0]; 
      seq_len = seq_j - seq_i + 1;

      sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
      ESL_REALLOC(sub_fs_model->tsc,  sizeof(float)   * (sub_hmm->M+1) * p7P_NTRANS);
      p7_ProfileConfig_fs(sub_hmm, pli->bg, gcode, sub_fs_model, seq_len, p7_UNIGLOBAL);
      tp_0 = sub_fs_model->tsc;
      tp_M = sub_fs_model->tsc + sub_fs_model->M * p7P_NTRANS;
     
      
      hmm_i = ESL_MAX(split_hits[s]->dcl->ihmm-1, 1);
      hmm_j = ESL_MIN(split_hits[s+1]->dcl->jhmm, hmm->M-1);
      /* Add extra transtion position for insetions at begining and end of alignments*/
      tp_0[p7P_MM] = log(hmm->t[hmm_i][p7H_MM]);
      tp_0[p7P_MI] = log(hmm->t[hmm_i][p7H_MI]);
      tp_0[p7P_MD] = log(hmm->t[hmm_i][p7H_MD]);
      tp_0[p7P_IM] = log(hmm->t[hmm_i][p7H_IM]);
      tp_0[p7P_II] = log(hmm->t[hmm_i][p7H_II]);
      tp_0[p7P_DM] = log(hmm->t[hmm_i][p7H_DM]);
      tp_0[p7P_DD] = log(hmm->t[hmm_i][p7H_DD]);
  
      tp_M[p7P_MM] = log(hmm->t[hmm_j][p7H_MM]);
      tp_M[p7P_MI] = log(hmm->t[hmm_j][p7H_MI]);
      tp_M[p7P_MD] = log(hmm->t[hmm_j][p7H_MD]);
      tp_M[p7P_IM] = log(hmm->t[hmm_j][p7H_IM]);
      tp_M[p7P_II] = log(hmm->t[hmm_j][p7H_II]);
      tp_M[p7P_DM] = log(hmm->t[hmm_j][p7H_DM]);
      tp_M[p7P_DD] = log(hmm->t[hmm_j][p7H_DD]);

      vit_mx = p7_gmx_fs_Create(sub_fs_model->M, seq_len, seq_len, p7P_CODONS); 
      p7_fs_global_Viterbi(hit_seq->dsq+seq_i-1, gcode, seq_len, sub_fs_model, vit_mx, NULL);     
      tr = p7_trace_fs_Create();
      p7_fs_global_Trace(hit_seq->dsq+seq_i-1, seq_len, sub_fs_model, vit_mx, tr);
      p7_trace_fs_Index(tr);
     
      /*replace upstream hits trace with new trace */
      p7_trace_fs_Destroy(split_hits[s]->dcl->tr);
      split_hits[s]->dcl->tr = tr;
 
      /* Set new hmm and ali coords */   
      tmp_ihmm = split_hits[s]->dcl->ihmm;
      split_hits[s]->dcl->ihmm = tmp_ihmm + tr->hmmfrom[0] - 1;
      split_hits[s]->dcl->jhmm = tmp_ihmm + tr->hmmto[0]   - 1;
      
      tmp_iali = split_hits[s]->dcl->iali;
      if(revcomp) {
        split_hits[s]->dcl->iali = tmp_iali - tr->sqfrom[0]  + 1;
        split_hits[s]->dcl->jali = tmp_iali - tr->sqto[0]    + 1;
      }
      else {
        split_hits[s]->dcl->iali = tmp_iali + tr->sqfrom[0]  - 1;
        split_hits[s]->dcl->jali = tmp_iali + tr->sqto[0]    - 1; 
      }

      /* Calculate new ali score */
      free(split_hits[s]->dcl->scores_per_pos);
      split_hits[s]->dcl->scores_per_pos = NULL;
      p7_splice_ComputeAliScores_fs(split_hits[s]->dcl, tr, hit_seq->dsq+seq_i-1, sub_fs_model, hit_seq->abc); 
      
      /*Destroy downstream hit */
      p7_trace_fs_Destroy(split_hits[s+1]->dcl->tr);
      free(split_hits[s+1]->dcl->scores_per_pos);
      p7_hit_Destroy(split_hits[s+1]);

      for(x = s+1; x < hit_cnt-1; x++)
        split_hits[x] = split_hits[x+1];

      for(x = s; x < hit_cnt-2; x++)
        split_edges[x] = split_edges[x+1];
      hit_cnt--;    
  
      p7_hmm_Destroy(sub_hmm);
      p7_profile_fs_Destroy(sub_fs_model);
      p7_gmx_Destroy(vit_mx);       
    } 
    else s++;
  }   

  *num_hits = hit_cnt;
  return split_edges;

  ERROR:
    p7_hmm_Destroy(sub_hmm);
    p7_profile_fs_Destroy(sub_fs_model);
    p7_gmx_Destroy(vit_mx);
    if(split_edges != NULL) free(split_edges);
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
    if(!graph->node_in_graph[up]) continue;
    for(down = 0; down < graph->num_nodes; down++) {
      if(!graph->node_in_graph[down]) continue;
      
      if(!hit_upstream(th->hit[up]->dcl, th->hit[down]->dcl, graph->revcomp)) continue;
   
      if(p7_splicegraph_EdgeExists(graph,up,down)) continue;
    
      /* check that hmm and ali coords are spliceable */
      if(!hits_spliceable_short(th->hit[up]->dcl, th->hit[down]->dcl)) continue;
      
      /* get min and max coords and fetch sub sequence of the potential splice region*/
      up_min   = ESL_MIN(th->hit[up]->dcl->iali, th->hit[up]->dcl->jali);
      up_max   = ESL_MAX(th->hit[up]->dcl->iali, th->hit[up]->dcl->jali);

      down_min = ESL_MIN(th->hit[down]->dcl->iali, th->hit[down]->dcl->jali);
      down_max = ESL_MAX(th->hit[down]->dcl->iali, th->hit[down]->dcl->jali);

      seq_min  = ESL_MIN(up_min, down_min);
      seq_max  = ESL_MAX(up_max, down_max);
      splice_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, graph->revcomp);
   printf("up %d down %d\n", up+1, down+1);      
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
hit_upstream(P7_DOMAIN *upstream, P7_DOMAIN *downstream, int revcomp) {

  if(upstream->ihmm > downstream->ihmm || upstream->jhmm > downstream->jhmm)
    return FALSE;
    
 
  if (( revcomp  && upstream->jali <= downstream->iali) ||
     ((!revcomp) && upstream->jali >= downstream->iali))
    return FALSE;

  return TRUE;
  
}


int
hits_spliceable_short(P7_DOMAIN *upstream, P7_DOMAIN *downstream) {
  

  /* Are nodes close enough on hmm coords */ 
  if(upstream->jhmm + MAX_AMINO_EXT < downstream->ihmm) return FALSE;

  /* Are nodes close enough on seq coords */
  if(abs(downstream->iali - upstream->jali) > MAX_INTRON_SHORT) return FALSE;

  return TRUE;
  
}


int
hits_spliceable_long(P7_DOMAIN *upstream, P7_DOMAIN *downstream) {


  /* Are nodes close enough on hmm coords */
  if(upstream->jhmm + MAX_AMINO_EXT < downstream->ihmm) return FALSE;

  /* Are nodes close enough on seq coords to have been spliced with short intron */
  if(abs(downstream->iali - upstream->jali) <= MAX_INTRON_SHORT) return FALSE;

  /* Are nodes close enough on seq coords to bw spliced with long intron*/
  if(abs(downstream->iali - upstream->jali) > MAX_INTRON_LONG) return FALSE;

  return TRUE;

}

int
p7_splice_RecoverHits2(SPLICE_GRAPH *graph, SPLICE_SAVED_HITS *sh, SPLICE_PIPELINE *pli, P7_HMM *hmm, ESL_GENCODE *gcode, ESL_SQFILE *seq_file, int first, int last) 
{

  int i, j;
  int z1, z2;
  int seq_min, seq_max;
  int recovered;
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

  tmp_dom = p7_domain_Create_empty();

  recovered = FALSE;
  for(i = first; i <= last; i++) {
    
    if(sh->srt[i]->duplicate) continue;
    if(sh->srt[i]->node_id >= 0) continue;
  
    tmp_dom->ihmm = sh->srt[i]->hmm_start;
    tmp_dom->jhmm = sh->srt[i]->hmm_end;    
    tmp_dom->iali = sh->srt[i]->seq_start;    
    tmp_dom->jali = sh->srt[i]->seq_end;
    
    /* Check if hit info is for a hit that is spliceable upstream or downstream of any hits int the graph*/
    for(j = 0; j < graph->split_N; j++) {
      if(!graph->node_in_graph[j]) continue;
      if(!((hit_upstream(tmp_dom, th->hit[j]->dcl, graph->revcomp) && hits_spliceable_short(tmp_dom, th->hit[j]->dcl)) ||
           (hit_upstream(th->hit[j]->dcl, tmp_dom, graph->revcomp) && hits_spliceable_short(th->hit[j]->dcl, tmp_dom)))) continue;

      new_hit          = p7_hit_Create_empty();
      new_hit->dcl     = p7_domain_Create_empty();
      new_hit->dcl->tr = p7_trace_fs_Create();

      seq_min = ESL_MIN(tmp_dom->iali, tmp_dom->jali);
      seq_max = ESL_MAX(tmp_dom->iali, tmp_dom->jali); 
   
      hit_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, graph->revcomp);

      sub_hmm = p7_splice_GetSubHMM(hmm, tmp_dom->ihmm, tmp_dom->jhmm);
      sub_hmm->fs = 0.;
      
      sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
      p7_ProfileConfig_fs(sub_hmm, pli->bg, gcode, sub_fs_model, hit_seq->n, p7_UNILOCAL); 
            
      vit_mx = p7_gmx_fs_Create(sub_fs_model->M, hit_seq->n, hit_seq->n, 0);
      p7_trans_Viterbi(hit_seq->dsq, gcode, hit_seq->n, sub_fs_model, vit_mx, NULL); 
      p7_trans_VTrace(hit_seq->dsq, hit_seq->n, gcode, sub_fs_model, vit_mx, new_hit->dcl->tr);

      for(z1 = 0; z1 < new_hit->dcl->tr->N; z1++)    if(new_hit->dcl->tr->st[z1] == p7T_M) break;
      for(z2 = new_hit->dcl->tr->N-1; z2 >= 0; z2--) if(new_hit->dcl->tr->st[z2] == p7T_M) break;
      
      new_hit->dcl->ihmm = tmp_dom->ihmm + new_hit->dcl->tr->k[z1] - 1;
      new_hit->dcl->jhmm = tmp_dom->ihmm + new_hit->dcl->tr->k[z2] - 1; 
      new_hit->dcl->iali = tmp_dom->iali + new_hit->dcl->tr->i[z1] - 3;
      if(graph->revcomp) 
        new_hit->dcl->jali = tmp_dom->iali - new_hit->dcl->tr->i[z2] + 1;
      else
        new_hit->dcl->jali = tmp_dom->iali + new_hit->dcl->tr->i[z2] - 1;

      p7_splice_ComputeAliScores_fs(new_hit->dcl, new_hit->dcl->tr, hit_seq->dsq, sub_fs_model, hit_seq->abc);      

      p7_splicegraph_AddNode(graph, new_hit);
      sh->srt[i]->node_id = graph->recover_N;
      graph->recover_N++;
      recovered = TRUE;     

      p7_hmm_Destroy(sub_hmm);
      p7_profile_fs_Destroy(sub_fs_model);
      p7_gmx_Destroy(vit_mx);
      esl_sq_Destroy(hit_seq);
      break;
    } 
  
  }     
 
  while (recovered) {

    recovered = FALSE;
    for(i = first; i <= last; i++) {

      if(sh->srt[i]->duplicate) continue;
      if(sh->srt[i]->node_id >= 0) continue;
  
      tmp_dom->ihmm = sh->srt[i]->hmm_start;
      tmp_dom->jhmm = sh->srt[i]->hmm_end;
      tmp_dom->iali = sh->srt[i]->seq_start;
      tmp_dom->jali = sh->srt[i]->seq_end;
  
      /* Check if hit info is for a hit that is spliceable upstream or downstream of any hits int the graph*/
      for(j = graph->split_N; j < graph->recover_N; j++) {
        if(!graph->node_in_graph[j]) continue;
        if(!((hit_upstream(tmp_dom, th->hit[j]->dcl, graph->revcomp) && hits_spliceable_short(tmp_dom, th->hit[j]->dcl)) ||
             (hit_upstream(th->hit[j]->dcl, tmp_dom, graph->revcomp) && hits_spliceable_short(th->hit[j]->dcl, tmp_dom)))) continue;
  
        new_hit          = p7_hit_Create_empty();
        new_hit->dcl     = p7_domain_Create_empty();
        new_hit->dcl->tr = p7_trace_fs_Create();
  
        seq_min = ESL_MIN(tmp_dom->iali, tmp_dom->jali);
        seq_max = ESL_MAX(tmp_dom->iali, tmp_dom->jali);
        hit_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, graph->revcomp);
  
        sub_hmm = p7_splice_GetSubHMM(hmm, tmp_dom->ihmm, tmp_dom->jhmm);
        sub_hmm->fs = 0.;
  
        sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
        p7_ProfileConfig_fs(sub_hmm, pli->bg, gcode, sub_fs_model, hit_seq->n, p7_UNILOCAL); 
  
        vit_mx = p7_gmx_fs_Create(sub_fs_model->M, hit_seq->n, hit_seq->n, 0);
  
        p7_trans_Viterbi(hit_seq->dsq, gcode, hit_seq->n, sub_fs_model, vit_mx, NULL);
        p7_trans_VTrace(hit_seq->dsq, hit_seq->n, gcode, sub_fs_model, vit_mx, new_hit->dcl->tr);
  
        for(z1 = 0; z1 < new_hit->dcl->tr->N; z1++)    if(new_hit->dcl->tr->st[z1] == p7T_M) break;
        for(z2 = new_hit->dcl->tr->N-1; z2 >= 0; z2--) if(new_hit->dcl->tr->st[z2] == p7T_M) break;
  
        new_hit->dcl->ihmm = tmp_dom->ihmm + new_hit->dcl->tr->k[z1] - 1;
        new_hit->dcl->jhmm = tmp_dom->ihmm + new_hit->dcl->tr->k[z2] - 1;
        new_hit->dcl->iali = tmp_dom->iali + new_hit->dcl->tr->i[z1] - 3;
        if(graph->revcomp)
          new_hit->dcl->jali = tmp_dom->iali - new_hit->dcl->tr->i[z2] + 1;
        else
          new_hit->dcl->jali = tmp_dom->iali + new_hit->dcl->tr->i[z2] - 1;
  
        p7_splice_ComputeAliScores_fs(new_hit->dcl, new_hit->dcl->tr, hit_seq->dsq, sub_fs_model, hit_seq->abc);
  
        p7_splicegraph_AddNode(graph, new_hit);
 
        sh->srt[i]->node_id = graph->recover_N; 
        graph->recover_N++;
        recovered = TRUE;
  
        p7_hmm_Destroy(sub_hmm);
        p7_profile_fs_Destroy(sub_fs_model);
        p7_gmx_Destroy(vit_mx);
        esl_sq_Destroy(hit_seq);
        break;
      }  
    }
  }

    
  free(tmp_dom);
  return eslOK;
}

int
p7_splice_RecoverHits(SPLICE_GRAPH *graph, SPLICE_SAVED_HITS *sh, SPLICE_PIPELINE *pli, P7_HMM *hmm, ESL_GENCODE *gcode, ESL_SQFILE *seq_file, int first, int last) 
{

  int i, j, s;
  int z1, z2;
  int gap_len;
  int seq_min, seq_max;
  int amino_len;
  int nuc_seq_idx;
  int amino;
  P7_DOMAIN     *tmp_dom;
  P7_HIT        *new_hit;
  P7_TOPHITS    *th;
  P7_HMM        *sub_hmm;
  P7_PROFILE    *sub_model;
  P7_GMX        *vit_mx;
  ESL_SQ        *hit_nuc_seq; 
  ESL_DSQ       *hit_aa_dsq;
  int status;


  th = graph->th;
  graph->recover_N = graph->split_N;
   
  //p7_splicehits_Dump(stdout, sh);
  p7_splicehits_FindNodes(graph, sh, first, last);

  tmp_dom = p7_domain_Create_empty();

  for(i = first; i <= last; i++) {
    
    if(sh->srt[i]->duplicate) continue;
    if(sh->srt[i]->node_id >= 0) continue;
    
    tmp_dom->ihmm = sh->srt[i]->hmm_start;
    tmp_dom->jhmm = sh->srt[i]->hmm_end;    
    tmp_dom->iali = sh->srt[i]->seq_start;    
    tmp_dom->jali = sh->srt[i]->seq_end;

    /* Check if hit info is for a hit that is spliceable upstream or downstream of any hits int the graph*/
    for(j = 0; j < graph->split_N; j++) {
      if(!graph->node_in_graph[j]) continue;
      if(hit_upstream(tmp_dom, th->hit[j]->dcl, graph->revcomp)) {
        if(graph->revcomp) gap_len = tmp_dom->jali - th->hit[j]->dcl->iali - 1;
        else               gap_len = th->hit[j]->dcl->iali - tmp_dom->jali - 1;
        if(gap_len > MAX_GAP_RANGE) continue;
      }
      else if(hit_upstream(th->hit[j]->dcl, tmp_dom, graph->revcomp)) {
        if(graph->revcomp) gap_len = th->hit[j]->dcl->jali - tmp_dom->iali - 1;
        else               gap_len = tmp_dom->iali - th->hit[j]->dcl->jali - 1;
        if(gap_len > MAX_GAP_RANGE) continue;
      }
      else continue;
      
      new_hit          = p7_hit_Create_empty();
      new_hit->dcl     = p7_domain_Create_empty();
      new_hit->dcl->tr = p7_trace_fs_Create();

      seq_min = ESL_MIN(tmp_dom->iali, tmp_dom->jali);
      seq_max = ESL_MAX(tmp_dom->iali, tmp_dom->jali); 
      hit_nuc_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, graph->revcomp);

      amino_len = hit_nuc_seq->n/3;
      ESL_ALLOC(hit_aa_dsq, sizeof(ESL_DSQ) * (amino_len+2));
      hit_aa_dsq[0] = eslDSQ_SENTINEL;

      nuc_seq_idx = 1;

      for(s = 1; s <= amino_len; s++) {
        amino = esl_gencode_GetTranslation(gcode,&hit_nuc_seq->dsq[nuc_seq_idx]);
        hit_aa_dsq[s] = amino;
        nuc_seq_idx+=3;
      }
      hit_aa_dsq[amino_len+1] = eslDSQ_SENTINEL;    

      sub_hmm = p7_splice_GetSubHMM(hmm, tmp_dom->ihmm, tmp_dom->jhmm);
      sub_hmm->fs = 0.;
      
      sub_model = p7_profile_Create(sub_hmm->M, sub_hmm->abc);
      p7_ProfileConfig(sub_hmm, pli->bg, sub_model, amino_len, p7_UNILOCAL);
      //sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
      //p7_ProfileConfig_fs(sub_hmm, pli->bg, gcode, sub_fs_model, hit_seq->n, p7_UNILOCAL); 
          
      //vit_mx = p7_gmx_fs_Create(sub_fs_model->M, hit_seq->n, hit_seq->n, 0);
      vit_mx = p7_gmx_Create(sub_model->M, amino_len);
      //p7_trans_Viterbi(hit_seq->dsq, gcode, hit_seq->n, sub_fs_model, vit_mx, NULL); 
     // p7_trans_VTrace(hit_seq->dsq, hit_seq->n, gcode, sub_fs_model, vit_mx, new_hit->dcl->tr);
      p7_GViterbi(hit_aa_dsq, amino_len, sub_model, vit_mx, NULL);
      p7_GTrace(hit_aa_dsq, amino_len, sub_model, vit_mx, new_hit->dcl->tr);
      p7_splice_ComputeAliScores(new_hit->dcl, new_hit->dcl->tr, hit_aa_dsq, sub_model, hmm->fs);      
      
      if(graph->revcomp)
        p7_trace_fs_Convert(new_hit->dcl->tr, hit_nuc_seq->L - hit_nuc_seq->start + 1, 1);
      else
        p7_trace_fs_Convert(new_hit->dcl->tr, hit_nuc_seq->start, 1);

      for(z1 = 0; z1 < new_hit->dcl->tr->N; z1++)    if(new_hit->dcl->tr->st[z1] == p7T_M) break;
      for(z2 = new_hit->dcl->tr->N-1; z2 >= 0; z2--) if(new_hit->dcl->tr->st[z2] == p7T_M) break;
       
      new_hit->dcl->ihmm = tmp_dom->ihmm + new_hit->dcl->tr->k[z1] - 1;
      new_hit->dcl->jhmm = tmp_dom->ihmm + new_hit->dcl->tr->k[z2] - 1; 
      if(graph->revcomp) {
        new_hit->dcl->iali = new_hit->dcl->tr->i[z2]-2;
        new_hit->dcl->jali = new_hit->dcl->tr->i[z1];
      }
      else {
        new_hit->dcl->iali = new_hit->dcl->tr->i[z1] - 2;
        new_hit->dcl->jali = new_hit->dcl->tr->i[z2];
      }
      
      p7_splicegraph_AddNode(graph, new_hit);
      
      sh->srt[i]->node_id = graph->recover_N;
      graph->recover_N++;

      p7_hmm_Destroy(sub_hmm);
      //p7_profile_fs_Destroy(sub_fs_model);
      p7_profile_Destroy(sub_model);
      p7_gmx_Destroy(vit_mx);
      esl_sq_Destroy(hit_nuc_seq);
      free(hit_aa_dsq);
      break;
    } 
  
  }     
 
  //p7_splicehits_Dump(stdout, sh);
  free(tmp_dom);
  return eslOK;

  ERROR:
	free(tmp_dom);
	free(hit_aa_dsq);
	return status;
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
  int new_seq_min, new_seq_max;
  int new_seq_len;
  int new_hmm_min, new_hmm_max;
  int new_hmm_len;
  int hit_cnt;
  int num_gaps;
  int gap_alloc; 
  int        *gap_merged;
  P7_TOPHITS *th;
  P7_HIT     **gap_hits;
  SPLICE_GAP **gaps;
  SPLICE_GAP *new_gap; 
  int        status;   

  th = graph->th;

  num_gaps  = 0;
  gap_alloc = th->N;
  ESL_ALLOC(gaps, sizeof(SPLICE_GAP*) * gap_alloc);

  /* Find the total set of possible gaps */  
  for(up = 0; up < th->N; up++) {
    if(!graph->node_in_graph[up]) continue;

    for(down = 0; down < th->N; down++) {
      if(!graph->node_in_graph[down] || up == down) continue;
      
      if(!hit_upstream(th->hit[up]->dcl, th->hit[down]->dcl, graph->revcomp)) continue;
      
      if(p7_splicegraph_PathExists(graph, up, down))  continue;
      
      hmm_gap_len = th->hit[down]->dcl->ihmm - th->hit[up]->dcl->jhmm - 1;
      if(hmm_gap_len < MIN_AMINO_OVERLAP) continue;
      
      if(graph->revcomp)
        seq_gap_len = th->hit[up]->dcl->jali - th->hit[down]->dcl->iali - 1;
      else  
        seq_gap_len = th->hit[down]->dcl->iali - th->hit[up]->dcl->jali - 1;
      if(seq_gap_len < 1 || seq_gap_len > MAX_GAP_RANGE) continue;
      
      if(num_gaps == gap_alloc) {
        gap_alloc*=2;
        ESL_REALLOC(gaps, sizeof(SPLICE_GAP*) * gap_alloc);
      }

      new_gap = p7_splicegap_Create();
      new_gap->up_node   = up;
      new_gap->down_node = down;
      
      new_gap->hmm_min = th->hit[up]->dcl->jhmm   + 1; 
      new_gap->hmm_max = th->hit[down]->dcl->ihmm - 1;  
      new_gap->hmm_len = new_gap->hmm_max - new_gap->hmm_min + 1;
      new_gap->seq_min = (graph->revcomp? th->hit[down]->dcl->iali + 1 : th->hit[up]->dcl->jali   + 1);
      new_gap->seq_max = (graph->revcomp? th->hit[up]->dcl->jali   - 1 : th->hit[down]->dcl->iali - 1);
      new_gap->seq_len = new_gap->seq_max - new_gap->seq_min + 1;
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
  
  new_gap = p7_splicegap_Create(); 
   
  /* Merge overlapping gaps */
  for(i = 0; i < num_gaps; i++) {
    if(gap_merged[i]) continue;

    /* Ensure seq gap is long engough to accomidate hmm gap */
    if(gaps[i]->seq_len < gaps[i]->hmm_len * 3) {
      gap_merged[i] = TRUE;
      continue; 
    }
    new_gap->hmm_min = gaps[i]->hmm_min;
    new_gap->hmm_max = gaps[i]->hmm_max; 
    new_gap->seq_min = gaps[i]->seq_min;
    new_gap->seq_max = gaps[i]->seq_max;

    gap_merged[i] = TRUE;
    for(j = 0; j < num_gaps; j++) {
      if(gap_merged[j]) continue;
      /* Ensure seq gap is long engough to accomidate hmm gap */
      if(gaps[j]->seq_len < gaps[j]->hmm_len * 3) {
        gap_merged[j] = TRUE;
        continue;      
      }
   
      /* Ensure gaps are upstream/downstream compatible */
      if(gaps[i]->up_node != gaps[j]->up_node) {
        if(!(hit_upstream(th->hit[gaps[i]->up_node]->dcl, th->hit[gaps[j]->up_node]->dcl, graph->revcomp)) &&
           !(hit_upstream(th->hit[gaps[j]->up_node]->dcl, th->hit[gaps[i]->up_node]->dcl, graph->revcomp))) continue;
      } 
  
      if(gaps[i]->down_node != gaps[j]->down_node) {
        if(!(hit_upstream(th->hit[gaps[i]->down_node]->dcl, th->hit[gaps[j]->down_node]->dcl, graph->revcomp)) &&
           !(hit_upstream(th->hit[gaps[j]->down_node]->dcl, th->hit[gaps[i]->down_node]->dcl, graph->revcomp))) continue;
      }
 
      hmm_overlap_min = ESL_MAX(new_gap->hmm_min, gaps[j]->hmm_min);
      hmm_overlap_max = ESL_MIN(new_gap->hmm_max, gaps[j]->hmm_max);
      
      hmm_overlap_len = hmm_overlap_max - hmm_overlap_min + 1;
      if(hmm_overlap_len < 1) continue;

      seq_overlap_min = ESL_MAX(new_gap->seq_min, gaps[j]->seq_min);
      seq_overlap_max = ESL_MIN(new_gap->seq_max, gaps[j]->seq_max);
 
      seq_overlap_len = seq_overlap_max - seq_overlap_min + 1;
      if(seq_overlap_len < 1) continue;
      
      new_seq_min = ESL_MIN(new_gap->seq_min, gaps[j]->seq_min);
      new_seq_max = ESL_MAX(new_gap->seq_max, gaps[j]->seq_max);
      new_seq_len = new_seq_max - new_seq_min + 1;
      if(new_seq_len > MAX_GAP_RANGE) continue;
      
      new_hmm_min = ESL_MIN(new_gap->hmm_min, gaps[j]->hmm_min);
      new_hmm_max = ESL_MAX(new_gap->hmm_max, gaps[j]->hmm_max);
      new_hmm_len = new_hmm_max -new_hmm_min + 1;
      
      if(new_seq_len < new_hmm_len *3) continue;
   
      new_gap->hmm_min = new_hmm_min; 
      new_gap->hmm_max = new_hmm_max;
      new_gap->hmm_len = new_hmm_len; 
      new_gap->seq_min = new_seq_min; 
      new_gap->seq_max = new_seq_max;
      new_gap->seq_len = new_seq_len;
      gap_merged[j] = TRUE;
    }
    
    /* align merged gap */
    gap_hits = p7_splicegap_AlignGap(graph, new_gap, hmm, pli->bg, gcode, seq_file, &hit_cnt);
    
    for(h = 0; h < hit_cnt; h++) {
      p7_splicegraph_AddNode(graph, gap_hits[h]); 
    }
    free(gap_hits);
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


int
p7_splice_MergePaths(SPLICE_GRAPH *graph, SPLICE_PATH **path_accumulator, SPLICE_PIPELINE *pli, P7_HMM *hmm, ESL_GENCODE *gcode, ESL_SQFILE *seq_file, int *num_paths) 
{


  int          i, p;
  int          up, down;
  int          step;
  int          path_cnt;
  int          new_edge;
  int          up_min, up_max;
  int          down_min, down_max;
  int          seq_min, seq_max;
  int64_t      bound_min, bound_max;
  SPLICE_PATH *up_path;
  SPLICE_PATH *down_path;
  SPLICE_PATH *path;
  SPLICE_EDGE *edge;
  ESL_SQ      *splice_seq;

  path_cnt = *num_paths;

  new_edge = FALSE;

  /* See is last node in one path can be spliced to first node in another path using a longer maxiumum intron length */
  for(up = 0; up < path_cnt; up++) {
    up_path = path_accumulator[up];

    for(down = 0; down < path_cnt; down++) {
      if(up == down) continue;
      down_path = path_accumulator[down];
      
      if(hit_upstream(up_path->hits[up_path->path_len-1]->dcl, down_path->hits[0]->dcl, graph->revcomp) && 
         hits_spliceable_long(up_path->hits[up_path->path_len-1]->dcl, down_path->hits[0]->dcl)) {
         
        up_min   = ESL_MIN(up_path->hits[up_path->path_len-1]->dcl->iali, up_path->hits[up_path->path_len-1]->dcl->jali);
        up_max   = ESL_MAX(up_path->hits[up_path->path_len-1]->dcl->iali, up_path->hits[up_path->path_len-1]->dcl->jali);
        down_min = ESL_MIN(down_path->hits[0]->dcl->iali, down_path->hits[0]->dcl->jali);
        down_max = ESL_MAX(down_path->hits[0]->dcl->iali, down_path->hits[0]->dcl->jali); 

        seq_min = ESL_MIN(up_min, down_min);
        seq_max =  ESL_MAX(up_max, down_max);
        splice_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, graph->revcomp);
       
        edge = p7_spliceedge_ConnectHits(pli, up_path->hits[up_path->path_len-1]->dcl, down_path->hits[0]->dcl, hmm, gcode, splice_seq, graph->revcomp); 

        if(edge != NULL) {
          new_edge = TRUE;
     
          edge->upstream_node_id   = up_path->node_id[up_path->path_len-1];
          edge->downstream_node_id = down_path->node_id[0];

          p7_splicegraph_AddEdge(graph, edge);

        }
        esl_sq_Destroy(splice_seq);
      } 
    }
  }
  
  /* If we have no new edges we are done */
  if(!new_edge) return eslOK;

  esl_vec_ISet(graph->node_in_graph, graph->num_nodes, 0);

  /* Restore all hits and edges from paths with new edges */
  for(p = 0; p < path_cnt; p++) { 
    for(step = 0; step < path_accumulator[p]->path_len; step++) {
      graph->node_in_graph[path_accumulator[p]->node_id[step]] = TRUE;
        
      if(step > 0) {
        edge = p7_splicegraph_GetEdge(graph, path_accumulator[p]->node_id[step-1], path_accumulator[p]->node_id[step]);
        edge->splice_score = path_accumulator[p]->edge_scores[step];
      }
    } 
  }

  /*Destroy old paths */
  for(p = 0; p < path_cnt; p++) { 
      p7_splicepath_Destroy(path_accumulator[p]);
  }

 //p7_splicegraph_DumpGraph(stdout, graph, FALSE);
  path_cnt = 0; 
  while(esl_vec_ISum(graph->node_in_graph, graph->split_N)) {
    path = p7_splicepath_GetBestPath(graph);
 //  p7_splicepath_Dump(stdout,path); 
    path_accumulator[path_cnt] = path; 
    path_cnt++;

    for(step = 0; step < path->path_len; step++) { 
     graph->node_in_graph[path->node_id[step]] = FALSE;
     if(path->node_id[step] >= graph->orig_N && path->node_id[step] < graph->split_N) {
       for(i = 0; i < graph->split_N; i++) {
         if(graph->orig_hit_idx[path->node_id[step]] == graph->orig_hit_idx[i]) {
           graph->node_in_graph[i] = FALSE;
         }
       }
     }
   }

    bound_min = ESL_MIN(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]);
    bound_max = ESL_MAX(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]);
    p7_splice_EnforceRangeBounds(graph, bound_min, bound_max);
  
  }
 
  *num_paths = path_cnt;

  return eslOK;
}


int
p7_splice_AlignPath(SPLICE_GRAPH *graph, SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_TOPHITS *tophits, P7_OPROFILE *om, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *path_seq, int64_t db_nuc_cnt, float fs_prob, int *frameshift, int *success)
{

  int          i;
  int          pos, seq_pos;
  int          exon;
  int          shift;
  int          seq_idx;
  int          amino_len;
  int          amino;
  int          orf_len;
  int          remove_node;
  int          replace_node;
  int          contains_orig;
  int          path_seq_len;
  float        dom_bias;
  float        nullsc;
  float        dom_score;
  double       dom_lnP;
  int         *nuc_index;
  ESL_DSQ     *nuc_dsq;
  ESL_DSQ     *amino_dsq;
  ESL_SQ      *nuc_seq;
  ESL_SQ      *amino_seq;
  P7_HIT      *replace_hit;
  P7_HIT      *remove_hit;
  int          status;

  nuc_index = NULL;
  nuc_dsq   = NULL;
  amino_dsq = NULL;

  path_seq_len = 0;
  for (i = 0; i < path->path_len; i++)
    path_seq_len += abs(path->upstream_spliced_nuc_end[i+1] - path->downstream_spliced_nuc_start[i]) + 1;

  /* If the spliced seqeunce length is non-mod 3 send to farmshift alignment */
  if(path_seq_len % 3 != 0) {
    *frameshift = TRUE;
    *success    = FALSE;
    return eslOK;
  }

  ESL_ALLOC(nuc_index, sizeof(int64_t) * (path_seq_len+2));
  ESL_ALLOC(nuc_dsq,   sizeof(ESL_DSQ) * (path_seq_len+2));

  /* Copy spliced nucleotides into single sequence and track their original indicies */
  nuc_index[0] = -1;
  nuc_dsq[0]   = eslDSQ_SENTINEL;
  seq_idx   = 1;

  for (i = 0; i < path->path_len; i++) {
    if (path->revcomp) {
      for (pos = path->downstream_spliced_nuc_start[i]; pos >= path->upstream_spliced_nuc_end[i+1]; pos--) {
        seq_pos = path_seq->n - pos + path_seq->end;

        nuc_index[seq_idx] = seq_pos;
        nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;
      }
    }
    else {
      for (pos = path->downstream_spliced_nuc_start[i]; pos <= path->upstream_spliced_nuc_end[i+1]; pos++) {
        seq_pos = pos - path_seq->start + 1;

        nuc_index[seq_idx] = seq_pos;
        nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;
      }
    }
  }

  nuc_index[seq_idx] = -1;
  nuc_dsq[seq_idx]   = eslDSQ_SENTINEL;

  
  /* Translate spliced nucleotide sequence */
  amino_len = path_seq_len/3;
  ESL_ALLOC(amino_dsq, sizeof(ESL_DSQ) * (amino_len+2));

  amino_dsq[0] = eslDSQ_SENTINEL;

  seq_idx = 1;

  for(i = 1; i <= amino_len; i++) {

    amino = esl_gencode_GetTranslation(gcode,&nuc_dsq[seq_idx]);

    if(esl_abc_XIsNonresidue(gm->abc, amino)) {
      *frameshift = TRUE;
      *success    = FALSE;
      free(nuc_index);
      free(nuc_dsq);
      free(amino_dsq);

      return eslOK;
    }

    amino_dsq[i] = amino;
    seq_idx+=3;
  }
  amino_dsq[amino_len+1] = eslDSQ_SENTINEL;

  /* Create ESL_SQs from ESL_DSQs */
  amino_seq   = esl_sq_CreateDigitalFrom(gcode->nt_abc, path_seq->name, amino_dsq, amino_len,    NULL,NULL,NULL);
  nuc_seq     = esl_sq_CreateDigitalFrom(gcode->aa_abc, path_seq->name, nuc_dsq,   path_seq_len, NULL,NULL,NULL);

  pli->nuc_sq       = nuc_seq;
  pli->amino_sq     = amino_seq;
  pli->orig_nuc_idx = nuc_index;

  /* Algin the splices Amino sequence */
  status = align_spliced_path(pli, om, gm, path_seq, gcode, fs_prob);

   /* Alignment failed */
  if(pli->hit == NULL || pli->hit->dcl->ad->exon_cnt == 1) {

    if(nuc_dsq   != NULL) free(nuc_dsq);
    if(amino_dsq != NULL) free(amino_dsq);
     
    *success = FALSE;
     return eslOK;
  }

  /* adjust all coords in hit and path */
  if(path->revcomp) {
    pli->hit->dcl->ad->sqfrom += 2;
    pli->hit->dcl->ienv = path_seq->n - pli->orig_nuc_idx[1]              + path_seq->end;
    pli->hit->dcl->jenv = path_seq->n - pli->orig_nuc_idx[pli->nuc_sq->n] + path_seq->end;
  }
  else {
    pli->hit->dcl->ienv = pli->orig_nuc_idx[1]              + path_seq->start -1;
    pli->hit->dcl->jenv = pli->orig_nuc_idx[pli->nuc_sq->n] + path_seq->start -1;
  }

  /* Adjust spliced hit score from amino_len to om->max_length */
  dom_score  = pli->hit->dcl->envsc;
  orf_len = pli->hit->dcl->ad->orfto - pli->hit->dcl->ad->orffrom + 1;
  dom_score -= 2 * log(2. / (amino_len+2));
  dom_score += 2 * log(2. / (om->max_length+2));
  dom_score -= (amino_len-orf_len)      * log((float) (amino_len) / (float) (amino_len+2));
  dom_score += (om->max_length-orf_len) * log((float) om->max_length / (float) (om->max_length+2));

   /* Bias calculation and adjustments */
  if(pli->do_null2)
    dom_bias = p7_FLogsum(0.0, log(pli->bg->omega) + pli->hit->dcl->domcorrection);
  else
    dom_bias = 0.;

  p7_bg_SetLength(pli->bg, om->max_length);
  p7_bg_NullOne  (pli->bg, pli->amino_sq->dsq, om->max_length, &nullsc);
  dom_score = (dom_score - (nullsc + dom_bias))  / eslCONST_LOG2;

  /* Add splice signal penalties */
  for(i = 0; i < path->path_len; i++)
     dom_score += path->signal_scores[i];

  dom_lnP   = esl_exp_logsurv(dom_score, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

  /* E-value adjusment */
  dom_lnP += log((float)db_nuc_cnt / (float)om->max_length);

  if ((pli->by_E && exp(dom_lnP) <= pli->E) || ((!pli->by_E) && dom_score >= pli->T)) {

    *success = TRUE;

    if ( path->path_len >  pli->hit->dcl->ad->exon_cnt) {
      /* Shift the path to start at the first hit that was inculded in the alignment
       * and end at the last hit that was included in the alignment */
      if(path->revcomp) {
        for(shift = 0; shift < path->path_len; shift++) {
          if(path->upstream_spliced_nuc_end[shift+1] <= pli->hit->dcl->iali) break;
        }
      }
      else {
        for(shift = 0; shift < path->path_len; shift++) {
           if(path->upstream_spliced_nuc_end[shift+1] >= pli->hit->dcl->iali) break;
        }
      }

      /* Ensure path will still contains an original hit after shifting */
      contains_orig = FALSE;
      for(exon = shift; exon < pli->hit->dcl->ad->exon_cnt ; exon ++ ) {
        if(path->node_id[exon] < graph->split_N) contains_orig = TRUE;
      }

      if(!contains_orig) {
        *success = FALSE;
        if(nuc_dsq   != NULL) free(nuc_dsq);
        if(amino_dsq != NULL) free(amino_dsq);
        
        return eslOK;
      }

      /* Shift path to start at frist hits that is in alignment */
      path->path_len = pli->hit->dcl->ad->exon_cnt;

      for(exon = 0; exon < path->path_len; exon++) {
        path->node_id[exon]                        = path->node_id[shift+exon];
        path->upstream_spliced_amino_end[exon]     = path->upstream_spliced_amino_end[shift+exon];
        path->downstream_spliced_amino_start[exon] = path->downstream_spliced_amino_start[shift+exon];
        path->upstream_spliced_nuc_end[exon]       = path->upstream_spliced_nuc_end[shift+exon];
        path->downstream_spliced_nuc_start[exon]   = path->downstream_spliced_nuc_start[shift+exon];
        path->hit_scores[exon]                     = path->hit_scores[shift+exon];
        path->edge_scores[exon]                    = path->edge_scores[shift+exon];
        path->hit_scores[exon]                     = path->hit_scores[shift+exon];
        path->edge_scores[exon]                    = path->edge_scores[shift+exon];
        path->signal_scores[exon]                  = path->signal_scores[shift+exon];
        path->hits[exon]                           = path->hits[shift+exon];
        path->split[exon]                          = path->split[shift+exon];
      }
      path->downstream_spliced_nuc_start[0]   = pli->hit->dcl->iali;
      path->downstream_spliced_amino_start[0] = pli->hit->dcl->ihmm;

      for(exon = 1; exon < path->path_len; exon++) {
        if (path->downstream_spliced_nuc_start[exon] > pli->hit->dcl->jali)
          break;
      }
 
      path->path_len =  pli->hit->dcl->ad->exon_cnt;

      path->upstream_spliced_nuc_end[exon]   = pli->hit->dcl->jali;
      path->upstream_spliced_amino_end[exon] = pli->hit->dcl->jhmm;
    }

    /* Find the first original hit in path to copy info*/
    i = 0;
    while( path->node_id[i] >= graph->split_N) {
      pli->hit->dcl->ad->exon_orig[i] = FALSE;
      pli->hit->dcl->ad->exon_split[i] = path->split[i];
      i++;
    }
    pli->hit->dcl->ad->exon_orig[i] = TRUE;
    pli->hit->dcl->ad->exon_split[i] = path->split[i];

    replace_node = path->node_id[i];
    replace_hit = tophits->hit[graph->orig_hit_idx[replace_node]];    
    p7_domain_Destroy(replace_hit->dcl);

    replace_hit->dcl        = pli->hit->dcl;
    replace_hit->frameshift = FALSE;

    replace_hit->flags =  0;
    replace_hit->flags |= p7_IS_REPORTED;
    replace_hit->flags |= p7_IS_INCLUDED;
    replace_hit->nreported = 1;
    replace_hit->nincluded = 1;

    replace_hit->dcl->bitscore    = dom_score;
    replace_hit->dcl->lnP         = dom_lnP;
    replace_hit->dcl->dombias     = dom_bias;
    replace_hit->dcl->is_reported = TRUE;
    replace_hit->dcl->is_included = TRUE;

    replace_hit->pre_score = pli->hit->dcl->envsc  / eslCONST_LOG2;
    replace_hit->pre_lnP   = esl_exp_logsurv (replace_hit->pre_score, om->evparam[p7_FTAUFS], om->evparam[p7_FLAMBDA]);

    replace_hit->sum_score  = replace_hit->score  = dom_score;
    replace_hit->sum_lnP    = replace_hit->lnP    = dom_lnP;

    replace_hit->sortkey    = pli->inc_by_E ? -dom_lnP : dom_score;
      
    /* Set all other original hits in alignment to unreportable */
    i++;
    for(  ; i < path->path_len; i++) {
      remove_node = path->node_id[i];
      if(remove_node >= graph->split_N) {
        pli->hit->dcl->ad->exon_orig[i] = FALSE;
        pli->hit->dcl->ad->exon_split[i] = path->split[i];
      }
      else {
        pli->hit->dcl->ad->exon_orig[i] = TRUE;
        pli->hit->dcl->ad->exon_split[i] = path->split[i];
        if(graph->orig_hit_idx[remove_node] == graph->orig_hit_idx[replace_node])
          continue;   
        remove_hit = tophits->hit[graph->orig_hit_idx[remove_node]];

        if(remove_hit->flags & p7_IS_REPORTED ) {
          tophits->nreported--;
          remove_hit->flags &= ~p7_IS_REPORTED;
          remove_hit->dcl->is_reported = FALSE;
        }
        if((remove_hit->flags & p7_IS_INCLUDED)) {
          tophits->nincluded--;
          remove_hit->flags &= ~p7_IS_INCLUDED;
          remove_hit->dcl->is_included = FALSE;
        }
      }
    }
    pli->hit->dcl = NULL;
  }
  else *success = FALSE;

  if(nuc_dsq   != NULL) free(nuc_dsq);
  if(amino_dsq != NULL) free(amino_dsq);
  
  return eslOK;

  ERROR:
    if(nuc_index != NULL) free(nuc_index);
    if(nuc_dsq   != NULL) free(nuc_dsq);
    if(amino_dsq != NULL) free(amino_dsq);
    return status;
   
}

int
align_spliced_path (SPLICE_PIPELINE *pli, P7_OPROFILE *om, P7_PROFILE *gm, ESL_SQ *target_seq, ESL_GENCODE *gcode, float fs_prob)
{

  int       i;
  int       splice_cnt;
  float     filtersc;
  float     envsc;
  float     seq_score;
  float     oasc;
  float     P;
  float     domcorrection;
  float    null2[p7_MAXCODE];
  P7_HIT   *hit;
  P7_TRACE *tr;
  int       status;

  hit          = p7_hit_Create_empty();
  hit->dcl     = p7_domain_Create_empty();
  hit->dcl->tr = NULL;
  hit->dcl->scores_per_pos = NULL;
  tr = p7_trace_CreateWithPP();

  p7_oprofile_ReconfigUnihit(om, pli->amino_sq->n);
  p7_omx_GrowTo(pli->fwd, om->M, pli->amino_sq->n, pli->amino_sq->n);
  p7_omx_GrowTo(pli->bwd, om->M, pli->amino_sq->n, pli->amino_sq->n);

  p7_bg_SetLength(pli->bg, pli->amino_sq->n);
  if (pli->do_biasfilter)
    p7_bg_FilterScore(pli->bg, pli->amino_sq->dsq, pli->amino_sq->n, &filtersc);
  else
    p7_bg_NullOne  (pli->bg, pli->amino_sq->dsq, pli->amino_sq->n, &filtersc);

  p7_Forward (pli->amino_sq->dsq, pli->amino_sq->n, om, pli->fwd, &envsc);

  seq_score = (envsc-filtersc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

  if (P > pli->F3) {
    p7_hit_Destroy(hit);
    p7_trace_Destroy(tr);
    return eslOK;
  }
  p7_Backward(pli->amino_sq->dsq, pli->amino_sq->n, om, pli->fwd, pli->bwd, NULL);

  if((status = p7_Decoding(om, pli->fwd, pli->bwd, pli->bwd)) == eslERANGE) goto ERROR;

  p7_OptimalAccuracy(om, pli->bwd, pli->fwd, &oasc);
  p7_OATrace        (om, pli->bwd, pli->fwd, tr);

  p7_trace_Index(tr);

  p7_splice_ComputeAliScores(hit->dcl, tr, pli->amino_sq->dsq, gm, fs_prob);

  hit->dcl->tr = p7_trace_splice_Convert(tr, pli->orig_nuc_idx, &splice_cnt);

  hit->dcl->ad = p7_alidisplay_splice_Create(hit->dcl->tr, 0, om, target_seq, pli->amino_sq, hit->dcl->scores_per_pos, tr->sqfrom[0], splice_cnt);

  p7_Null2_ByExpectation(om, pli->bwd, null2);
  domcorrection = 0.;
  for (i = 1; i <= pli->amino_sq->n; i++)
    domcorrection += logf(null2[pli->amino_sq->dsq[i]]);
  
  hit->dcl->domcorrection = domcorrection;

  hit->dcl->ihmm = hit->dcl->ad->hmmfrom;
  hit->dcl->jhmm = hit->dcl->ad->hmmto;

  /* Convert sqfrom, sqto to full sequence coords */
  if(target_seq->start < target_seq->end) {
    hit->dcl->ad->sqfrom  = hit->dcl->ad->sqfrom + target_seq->start - 1;
    hit->dcl->ad->sqto    = hit->dcl->ad->sqto   + target_seq->start - 1;
  } else {
    hit->dcl->ad->sqto    = target_seq->n - hit->dcl->ad->sqto   + target_seq->end;
    hit->dcl->ad->sqfrom  = target_seq->n - hit->dcl->ad->sqfrom + target_seq->end;
  }

  hit->dcl->iali = hit->dcl->ad->sqfrom;
  hit->dcl->jali = hit->dcl->ad->sqto;

  hit->dcl->envsc = envsc;
  hit->dcl->oasc  = oasc;
  hit->dcl->dombias       = 0.0;
  hit->dcl->bitscore      = 0.0;
  hit->dcl->lnP           = 0.0;
  hit->dcl->is_reported   = FALSE;
  hit->dcl->is_included   = FALSE;

  pli->hit = hit;

  p7_trace_Destroy(tr);
  return eslOK;

  ERROR:
    p7_trace_Destroy(tr);
    p7_hit_Destroy(hit);
    return status;
}

int
p7_splice_AlignFrameshiftPath(SPLICE_GRAPH *graph, SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_TOPHITS *tophits, P7_FS_PROFILE *gm_fs, ESL_GENCODE *gcode, ESL_SQ *path_seq, int64_t db_nuc_cnt, int *success)
{

  int          i;
  int          pos, seq_pos;
  int          exon;
  int          shift;
  int          seq_idx;
  int          ali_len;
  int          replace_node;
  int          remove_node;
  int          contains_orig;
  int          path_seq_len;
  float        dom_bias;
  float        nullsc;
  float        dom_score;
  double       dom_lnP;
  int         *nuc_index;
  ESL_DSQ     *nuc_dsq;
  ESL_SQ      *nuc_seq;
  P7_HIT      *replace_hit;
  P7_HIT      *remove_hit;
  int          status;
 
  nuc_index = NULL;
  nuc_dsq   = NULL;

  path_seq_len = 0;
  for (i = 0; i < path->path_len; i++)
    path_seq_len += abs(path->upstream_spliced_nuc_end[i+1] - path->downstream_spliced_nuc_start[i]) + 1;

  ESL_ALLOC(nuc_index, sizeof(int64_t) * (path_seq_len+2));
  ESL_ALLOC(nuc_dsq,   sizeof(ESL_DSQ) * (path_seq_len+2));

  /* Copy spliced nucleotides into single sequence and track their original indicies */
  nuc_index[0] = -1;
  nuc_dsq[0]   = eslDSQ_SENTINEL;
  seq_idx   = 1;

  for (i = 0; i < path->path_len; i++) {
    if (path->revcomp) {
      for (pos = path->downstream_spliced_nuc_start[i]; pos >= path->upstream_spliced_nuc_end[i+1]; pos--) {
        seq_pos = path_seq->n - pos + path_seq->end;

        nuc_index[seq_idx] = seq_pos;
        nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;
      }
    }
    else {
      for (pos = path->downstream_spliced_nuc_start[i]; pos <= path->upstream_spliced_nuc_end[i+1]; pos++) {
        seq_pos = pos - path_seq->start + 1;

        nuc_index[seq_idx] = seq_pos;
        nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;
      }
    }
  }

  nuc_index[seq_idx] = -1;
  nuc_dsq[seq_idx]   = eslDSQ_SENTINEL;

  nuc_seq = esl_sq_CreateDigitalFrom(gcode->aa_abc, path_seq->name, nuc_dsq, path_seq_len, NULL,NULL,NULL);

  pli->nuc_sq       = nuc_seq;
  pli->amino_sq     = NULL;
  pli->orig_nuc_idx = nuc_index;

  /* Algin the spliced Amino sequence */
  status = align_spliced_path_frameshift(pli, gm_fs, path_seq, gcode);

  /* Alignment failed */
  if(pli->hit == NULL ) {
    if(nuc_dsq   != NULL) free(nuc_dsq);
    *success = FALSE;
     
     return eslOK;
  }

  /* adjust all coords in hit and path */
  if(path->revcomp) {
    pli->hit->dcl->ad->sqfrom += 2;
    pli->hit->dcl->ienv = path_seq->n - pli->orig_nuc_idx[1]              + path_seq->end;
    pli->hit->dcl->jenv = path_seq->n - pli->orig_nuc_idx[pli->nuc_sq->n] + path_seq->end;
  }
  else {
    pli->hit->dcl->ienv = pli->orig_nuc_idx[1]              + path_seq->start -1;
    pli->hit->dcl->jenv = pli->orig_nuc_idx[pli->nuc_sq->n] + path_seq->start -1;
  }
 
  if(pli->hit->dcl->ad->exon_cnt > 1) {
    ali_len = 0;
    for( i = 0; i < pli->hit->dcl->ad->exon_cnt; i++) 
      ali_len += abs(pli->hit->dcl->ad->exon_seq_ends[i] - pli->hit->dcl->ad->exon_seq_starts[i]) + 1;
  }
  else 
    ali_len = abs(pli->hit->dcl->ad->sqto - pli->hit->dcl->ad->sqfrom) + 1;

  /* Adjust spliced hit score from nuc_len to gm_fs->max_length*3 */
  dom_score  = pli->hit->dcl->envsc;
  
  dom_score -= 2 * log(2. / ((pli->nuc_sq->n/3.)+2));
  dom_score += 2 * log(2. / (gm_fs->max_length+2));
  dom_score -= (pli->nuc_sq->n-ali_len)      * log((float) (pli->nuc_sq->n/3.) / (float) ((pli->nuc_sq->n/3.)+2));
  dom_score += (ESL_MAX(pli->nuc_sq->n, gm_fs->max_length*3)-ali_len) * log((float) gm_fs->max_length / (float) (gm_fs->max_length+2));

  /* Bias calculation and adjustments */
  if(pli->do_null2)
    dom_bias = p7_FLogsum(0.0, log(pli->bg->omega) + pli->hit->dcl->domcorrection);
  else
    dom_bias = 0.;

  p7_bg_SetLength(pli->bg, gm_fs->max_length*3);
  p7_bg_NullOne  (pli->bg, pli->nuc_sq->dsq, gm_fs->max_length*3, &nullsc);
  dom_score = (dom_score - (nullsc + dom_bias))  / eslCONST_LOG2;
  
  /* Add splice signal penalties */
  pli->hit->dcl->jali = pli->hit->dcl->ad->sqto;
  for(i = 0; i < path->path_len; i++)
    dom_score += path->signal_scores[i];
  
  dom_lnP   = esl_exp_logsurv(dom_score, gm_fs->evparam[p7_FTAUFS], gm_fs->evparam[p7_FLAMBDA]);
    
  /* E-value adjusment */
  dom_lnP += log((float)db_nuc_cnt / (float)gm_fs->max_length);

    
   if ((pli->by_E && exp(dom_lnP) <= pli->E) || ((!pli->by_E) && dom_score >= pli->T)) {

    *success = TRUE;

    if ( path->path_len > pli->hit->dcl->ad->exon_cnt) {
      /* Shift the path to start at the first hit that was inculded in the alignment
       * and end at the last hit that was included in the alignment */
      if(path->revcomp) {
        for(shift = 0; shift < path->path_len; shift++) {
          if(path->upstream_spliced_nuc_end[shift+1] <= pli->hit->dcl->iali) break;
        }
      }
      else {
        for(shift = 0; shift < path->path_len; shift++) {
           if(path->upstream_spliced_nuc_end[shift+1] >= pli->hit->dcl->iali) break;
        }
      }
      /* Ensure path will still contains an original hit after shifting */
      contains_orig = FALSE;
      for(exon = shift; exon < pli->hit->dcl->ad->exon_cnt ; exon ++ ) {
        if(path->node_id[exon] < graph->split_N) contains_orig = TRUE;
      }

      if(!contains_orig) {
        *success = FALSE;
        if(nuc_dsq   != NULL) free(nuc_dsq);
        
        return eslOK;
      }

      /* Shift path to start at frist hit that is in alignment */
      path->path_len = pli->hit->dcl->ad->exon_cnt;

      for(exon = 0; exon < path->path_len; exon++) {
        path->node_id[exon]                        = path->node_id[shift+exon];
        path->upstream_spliced_amino_end[exon]     = path->upstream_spliced_amino_end[shift+exon];
        path->downstream_spliced_amino_start[exon] = path->downstream_spliced_amino_start[shift+exon];
        path->upstream_spliced_nuc_end[exon]       = path->upstream_spliced_nuc_end[shift+exon];
        path->downstream_spliced_nuc_start[exon]   = path->downstream_spliced_nuc_start[shift+exon];
        path->hit_scores[exon]                     = path->hit_scores[shift+exon];
        path->edge_scores[exon]                    = path->edge_scores[shift+exon];
        path->signal_scores[exon]                  = path->signal_scores[shift+exon];
        path->hits[exon]                           = path->hits[shift+exon];
        path->split[exon]                          = path->split[shift+exon];
      }
      path->downstream_spliced_nuc_start[0]   = pli->hit->dcl->iali;
      path->downstream_spliced_amino_start[0] = pli->hit->dcl->ihmm;

      for(exon = 1; exon < path->path_len; exon++) {
        if (path->downstream_spliced_nuc_start[exon] > pli->hit->dcl->jali)
          break;
      }

      path->upstream_spliced_nuc_end[exon]   = pli->hit->dcl->jali;
      path->upstream_spliced_amino_end[exon] = pli->hit->dcl->jhmm;
    }
  
    /* Find first original hit to copy info */ 
    i = 0;
    while( path->node_id[i] >= graph->split_N) {

      pli->hit->dcl->ad->exon_orig[i] = FALSE;
      pli->hit->dcl->ad->exon_split[i] = path->split[i];
      i++;
    }

    pli->hit->dcl->ad->exon_orig[i] = TRUE;
    pli->hit->dcl->ad->exon_split[i] = path->split[i];

    replace_node = path->node_id[i];
    replace_hit  = tophits->hit[graph->orig_hit_idx[replace_node]];

    p7_domain_Destroy(replace_hit->dcl);

    replace_hit->dcl        = pli->hit->dcl;
    replace_hit->frameshift = TRUE;

    replace_hit->flags =  0;
    replace_hit->flags |= p7_IS_REPORTED;
    replace_hit->flags |= p7_IS_INCLUDED;
    replace_hit->nreported = 1;
    replace_hit->nincluded = 1;

    replace_hit->dcl->bitscore    = dom_score;
    replace_hit->dcl->lnP         = dom_lnP;
    replace_hit->dcl->dombias     = dom_bias;
    replace_hit->dcl->is_reported = TRUE;
    replace_hit->dcl->is_included = TRUE;

    replace_hit->pre_score = pli->hit->dcl->envsc  / eslCONST_LOG2;
    replace_hit->pre_lnP   = esl_exp_logsurv (replace_hit->pre_score, gm_fs->evparam[p7_FTAUFS], gm_fs->evparam[p7_FLAMBDA]);

    replace_hit->sum_score  = replace_hit->score  = dom_score;
    replace_hit->sum_lnP    = replace_hit->lnP    = dom_lnP;

    replace_hit->sortkey    = pli->inc_by_E ? -dom_lnP : dom_score;

    /* Set all original hits in alignment to unreportable */
    i++;
    for(  ; i < path->path_len; i++) {
      remove_node = path->node_id[i];
      if(remove_node >= graph->split_N) {
        pli->hit->dcl->ad->exon_orig[i] = FALSE;
        pli->hit->dcl->ad->exon_split[i] = path->split[i];
      }
      else {
        pli->hit->dcl->ad->exon_orig[i] = TRUE;
        pli->hit->dcl->ad->exon_split[i] = path->split[i];
        if(graph->orig_hit_idx[remove_node] == graph->orig_hit_idx[replace_node])
          continue;

        remove_hit = tophits->hit[graph->orig_hit_idx[remove_node]];

        if(remove_hit->flags & p7_IS_REPORTED ) {
          tophits->nreported--;
          remove_hit->flags &= ~p7_IS_REPORTED;
          remove_hit->dcl->is_reported = FALSE;
        }
        if((remove_hit->flags & p7_IS_INCLUDED)) {
          tophits->nincluded--;
          remove_hit->flags &= ~p7_IS_INCLUDED;
          remove_hit->dcl->is_included = FALSE;
        }
      }
    }
    pli->hit->dcl = NULL;
  }
  else *success = FALSE;
  
  if(nuc_dsq   != NULL) free(nuc_dsq);
  
  return eslOK;

  ERROR:
    if(nuc_index != NULL) free(nuc_index);
    if(nuc_dsq   != NULL) free(nuc_dsq);
    return status;

}

int
align_spliced_path_frameshift (SPLICE_PIPELINE *pli, P7_FS_PROFILE *gm_fs, ESL_SQ *target_seq, ESL_GENCODE *gcode)
{

  int       i, z;
  int       t, u, v, w, x;
  int       splice_cnt;
  int       codon_idx;
  float     nullsc;
  float     envsc;
  float     seq_score;
  float     oasc;
  float     P;
  float     domcorrection;
  float    null2[p7_MAXCODE];
  P7_GMX   *gxppfs = NULL;
  P7_HIT   *hit;
  P7_TRACE *tr;
  int       status;

  tr = p7_trace_fs_CreateWithPP();

  p7_fs_ReconfigUnihit(gm_fs, pli->nuc_sq->n);
  p7_gmx_fs_GrowTo(pli->gfwd, gm_fs->M, pli->nuc_sq->n, pli->nuc_sq->n, p7P_CODONS);
  p7_gmx_fs_GrowTo(pli->gbwd, gm_fs->M, pli->nuc_sq->n, pli->nuc_sq->n, 0);
 
  p7_bg_SetLength(pli->bg, pli->nuc_sq->n);
  p7_bg_NullOne  (pli->bg, pli->nuc_sq->dsq, pli->nuc_sq->n, &nullsc);

  p7_Forward_Frameshift(pli->nuc_sq->dsq, gcode, pli->nuc_sq->n, gm_fs, pli->gfwd, &envsc);

  seq_score = (envsc-nullsc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score, gm_fs->evparam[p7_FTAUFS],  gm_fs->evparam[p7_FLAMBDA]);

  if (P > pli->F3) {

    p7_trace_fs_Destroy(tr);
    return eslOK;
  }

  p7_Backward_Frameshift(pli->nuc_sq->dsq, gcode, pli->nuc_sq->n, gm_fs, pli->gbwd, NULL);

  if ((gxppfs = p7_gmx_fs_Create(gm_fs->M, pli->nuc_sq->n, pli->nuc_sq->n, p7P_CODONS)) == NULL) { status = eslFAIL; goto ERROR; }
  p7_Decoding_Frameshift(gm_fs, pli->gfwd, pli->gbwd, gxppfs);

  p7_OptimalAccuracy_Frameshift(gm_fs, gxppfs, pli->gbwd, &oasc);
  p7_OATrace_Frameshift(gm_fs, gxppfs, pli->gbwd, pli->gfwd, tr);   /* <tr>'s seq coords are offset by i-1, rel to orig dsq */

  p7_trace_Index(tr);

  hit          = p7_hit_Create_empty();
  hit->dcl     = p7_domain_Create_empty();
  hit->dcl->tr = NULL;

  p7_splice_ComputeAliScores_fs(hit->dcl, tr, pli->nuc_sq->dsq, gm_fs, target_seq->abc);

  if((hit->dcl->tr = p7_trace_splice_fs_Convert(tr, pli->orig_nuc_idx, &splice_cnt)) == NULL) { status = eslFAIL; goto ERROR; }

  if((hit->dcl->ad = p7_alidisplay_splice_fs_Create(hit->dcl->tr, 0, gm_fs, target_seq, pli->nuc_sq->dsq, gcode, hit->dcl->scores_per_pos, pli->orig_nuc_idx, tr->sqfrom[0], splice_cnt)) == NULL) { status = eslFAIL; goto ERROR; }

  p7_Null2_fs_ByExpectation(gm_fs, pli->gfwd, null2);
  domcorrection = 0.;
  t = u = v = w = x = -1;
  z = 0;
  i = 1;

  while(i <= pli->nuc_sq->n) {
    if(esl_abc_XIsCanonical(target_seq->abc, pli->nuc_sq->dsq[i])) x = pli->nuc_sq->dsq[i];
    else                                                           x = p7P_MAXCODONS;

    switch (tr->st[z]) {
      case p7T_N:
      case p7T_C:
      case p7T_J:  if(tr->i[z] == i && i > 2) i++;
                   z++;   break;
      case p7T_X:
      case p7T_S:
      case p7T_B:
      case p7T_E:
      case p7T_T:
      case p7T_D:  z++;   break;
      case p7T_M:  if(tr->i[z] == i)
                   {
                     if     (tr->c[z] == 1) { codon_idx = p7P_CODON1(x);             codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_QC2); }
                     else if(tr->c[z] == 2) { codon_idx = p7P_CODON2(w, x);          codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_QC1); }
                     else if(tr->c[z] == 3) { codon_idx = p7P_CODON3(v, w, x);       codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_C); }
                     else if(tr->c[z] == 4) { codon_idx = p7P_CODON4(u, v, w, x);    codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_QC1); }
                     else if(tr->c[z] == 5) { codon_idx = p7P_CODON5(t, u, v, w, x); codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_QC2); }
                     domcorrection += logf(null2[p7P_AMINO(gm_fs, tr->k[z], codon_idx)]);
                     z++;
                   }
                   i++;  break;
      case p7T_I:  if(tr->i[z] == i)
                   {
                     codon_idx = p7P_CODON3(v, w, x);
                     codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_C);
                     domcorrection += logf(null2[p7P_AMINO(gm_fs, tr->k[z], codon_idx)]);
                     z++;
                   }
                   i++;  break;
    }
    t = u;
    u = w;
    v = w;
    w = x;
  }

  hit->dcl->domcorrection = domcorrection;

  hit->dcl->ihmm = hit->dcl->ad->hmmfrom;
  hit->dcl->jhmm = hit->dcl->ad->hmmto;

  /* Convert sqfrom, sqto to full sequence coords */
  if(target_seq->start < target_seq->end) {
    hit->dcl->ad->sqfrom  = hit->dcl->ad->sqfrom + target_seq->start - 1;
    hit->dcl->ad->sqto    = hit->dcl->ad->sqto   + target_seq->start - 1;
  } else {
    hit->dcl->ad->sqto    = target_seq->n - hit->dcl->ad->sqto   + target_seq->end;
    hit->dcl->ad->sqfrom  = target_seq->n - hit->dcl->ad->sqfrom + target_seq->end;
  }
  hit->dcl->iali = hit->dcl->ad->sqfrom;
  hit->dcl->jali = hit->dcl->ad->sqto;

  hit->dcl->envsc = envsc;
  hit->dcl->oasc  = oasc;
  hit->dcl->dombias       = 0.0;
  hit->dcl->bitscore      = 0.0;
  hit->dcl->lnP           = 0.0;
  hit->dcl->is_reported   = FALSE;
  hit->dcl->is_included   = FALSE;

  pli->hit = hit;


  p7_trace_fs_Destroy(tr);
  p7_gmx_Destroy(gxppfs);
  return eslOK;

  ERROR:
    p7_trace_fs_Destroy(tr);
    p7_gmx_Destroy(gxppfs);
    p7_hit_Destroy(hit);
    return status;
}


int 
p7_splice_ReleaseHits(SPLICE_GRAPH *graph, int *hits_processed, int *num_hits_processed, int range_bound_min, int range_bound_max)
{

  int        i,j;
  int        hit_min, hit_max;
  int        overlap_min, overlap_max;
  int        overlap_len;
  int        hit_idx, split_idx;
  int        num_hits;
  P7_TOPHITS *th;
  P7_HIT     *hit;

  th = graph->th;
  num_hits = *num_hits_processed;

  for(i = 0; i < th->N; i++) {
    if(!graph->node_in_graph[i]) continue;

    hit = th->hit[i];

    hit_min = ESL_MIN(hit->dcl->iali, hit->dcl->jali);
    hit_max = ESL_MAX(hit->dcl->iali, hit->dcl->jali);
    overlap_min = ESL_MAX(hit_min, range_bound_min);
    overlap_max = ESL_MIN(hit_max, range_bound_max);
    overlap_len = overlap_max - overlap_min + 1;
    
    if(overlap_len > 0) { 
      /* If hit in path range set as not in garph */
      graph->node_in_graph[i] = FALSE;

      /* If hit in path range is original hit mark as processed */
      if(i < graph->orig_N) {
        num_hits++;
        hit_idx = graph->orig_hit_idx[i];
        hits_processed[hit_idx] = TRUE;
      }
      /* If hit in path range is split hit mark all split hits from 
       * the same original hit as not in graph and processed */
      else if(i < graph->split_N) {
        num_hits++;
        hit_idx = graph->orig_hit_idx[i];
        for(j = graph->orig_N; j < graph->split_N; j++) {
          if(i == j) continue;
          split_idx = graph->orig_hit_idx[j];
          if(hit_idx == split_idx && graph->node_in_graph[j]) {
            graph->node_in_graph[j] = FALSE;
            num_hits++;
          }
        }
        hits_processed[hit_idx] = TRUE;
      }
    }
  }

  *num_hits_processed = num_hits;
  return eslOK;

}

int
p7_splice_EnforceRangeBounds(SPLICE_GRAPH *graph, int64_t bound_min, int64_t bound_max) {

  int     up, down;
  int64_t up_hit_min, up_hit_max;
  int64_t down_hit_min, down_hit_max;
  int overlap_min, overlap_max, overlap_len;
  P7_HIT *up_hit;
  P7_HIT *down_hit;
  P7_TOPHITS  *th;
  SPLICE_EDGE *tmp_edge;

  th = graph->th;

  for(up = 0; up < th->N; up++) {
    for(down = 0; down < th->N; down++) {
      if( up == down) continue;

      tmp_edge = p7_splicegraph_GetEdge(graph, up, down);
      if(tmp_edge == NULL || tmp_edge->splice_score == -eslINFINITY) continue;

      up_hit   = th->hit[up];
      down_hit = th->hit[down];

      up_hit_min   = ESL_MIN(up_hit->dcl->iali, up_hit->dcl->jali);
      up_hit_max   = ESL_MAX(up_hit->dcl->iali, up_hit->dcl->jali);
      down_hit_min = ESL_MIN(down_hit->dcl->iali, down_hit->dcl->jali);
      down_hit_max = ESL_MAX(down_hit->dcl->iali, down_hit->dcl->jali);
      overlap_min = ESL_MAX(bound_min, ESL_MIN(up_hit_min, down_hit_min));
      overlap_max = ESL_MIN(bound_max, ESL_MAX(up_hit_max, down_hit_max));
      overlap_len = overlap_max - overlap_min + 1;
      
      if(overlap_len > 0) {
        tmp_edge->splice_score = -eslINFINITY;
        graph->num_edges--;
      } 
    }
  }

  return eslOK;

}

ESL_SQ* 
p7_splice_GetSubSequence(const ESL_SQFILE *seq_file, char* seqname, int64_t seq_min, int64_t seq_max, int revcomp)
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
p7_splice_ComputeAliScores(P7_DOMAIN *dom, P7_TRACE *tr, ESL_DSQ *amino_dsq, const P7_PROFILE *gm, float fs_prob)
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

  /* To keep scores consitstant when splcing frameshift and non frameshift hits all amino 
   * emissions scores are multiplied ny the stadard codon prbability from the fs model */
  while (z1 < z2) {
    if (tr->st[z1] == p7T_M) {
      dom->scores_per_pos[n] = p7P_MSC(gm, k, amino_dsq[i]) + log(1. - fs_prob*4);
      if (tr->st[z1-1] == p7T_I) 
        dom->scores_per_pos[n] += p7P_TSC(gm, k-1, p7P_IM);
      else if (tr->st[z1-1] == p7T_D)
        dom->scores_per_pos[n] += p7P_TSC(gm, k-1, p7P_DM); 
      i++; k++; z1++; n++;
 
      while(z1 < z2 && tr->st[z1] == p7T_M) {
        dom->scores_per_pos[n] =  p7P_MSC(gm, k,   amino_dsq[i]) + log(1. - fs_prob*4);
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
  int indel;
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
        tr->fs++;
      }
      else if(c == 2) {
        if(esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i]))
          codon_idx = p7P_CODON2(nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_QC1;
        tr->fs++;
      }
      else if(c == 3) {
        if(esl_abc_XIsCanonical(abc, nuc_dsq[i-2]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i]))
          codon_idx = p7P_CODON3(nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_C;
        /* record stop codon */
        indel = p7P_INDEL(gm_fs, k, codon_idx);
        if(indel == p7P_XXx || indel == p7P_XxX || indel == p7P_xXX) tr->fs++;
      }  
      else if(c == 4) {
        if(esl_abc_XIsCanonical(abc, nuc_dsq[i-3]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-2]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(abc, nuc_dsq[i]))
          codon_idx = p7P_CODON4(nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_QC1;
        tr->fs++;
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
        tr->fs++;
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



