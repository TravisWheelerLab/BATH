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

#ifdef HMMER_THREADS
#include <pthread.h>
#include "esl_threads.h"
#endif /*HMMER_THREADS*/

#include "hmmer.h"
#include "p7_splice.h"


static int serial_loop(SPLICE_WORKER_INFO *info, P7_TOPHITS *tophits, P7_TOPHITS *seed_hits);
#ifdef HMMER_THREADS
static int thread_loop(SPLICE_WORKER_INFO *info, P7_TOPHITS *tophits, P7_TOPHITS *seed_hits, int infocnt);
static void* splice_thread(void *arg);
#endif /*HMMER_THREADS*/
static int align_spliced_path (SPLICE_PIPELINE *pli, P7_OPROFILE *om, P7_PROFILE *gm, ESL_SQ *target_seq, ESL_GENCODE *gcode, float fs_prob);
static int align_spliced_path_frameshift (SPLICE_PIPELINE *pli, P7_FS_PROFILE *gm_fs, ESL_SQ *target_seq, ESL_GENCODE *gcode);
static int hits_spliceable(P7_DOMAIN *upstream, P7_DOMAIN *downstream);
static float hmm_overlap(P7_DOMAIN *dom_1, P7_DOMAIN *dom_2);
static float hmm_overlap2(SPLICE_PATH2 *path1, int step1, SPLICE_PATH2 *path2, int step2);
static float seq_overlap(P7_DOMAIN *dom_1, P7_DOMAIN *dom_2, int revcomp);
static float seq_overlap2(SPLICE_PATH2 *path1, int step1, SPLICE_PATH2 *path2, int step2); 
static int confirm_overlap(SPLICE_PATH2 *path1, int step1, SPLICE_PATH2 *path2, int step2, int confirm_side);
static int confirm_split(SPLICE_PATH2 *path1, SPLICE_PATH2 *path2);
static int remove_upstream(SPLICE_PATH2 *path1, SPLICE_PATH2 *path2, int *step2); 

static int assess_paths(SPLICE_WORKER_INFO *info, SPLICE_PATH2 *path1, SPLICE_PATH2 *path2, SPLICE_PATH2 *path3, P7_HMM *sub_hmm, ESL_SQ *path_seq, int full_intron, int *step);

/*  Function: p7_splice_SpliceHits
 *  Synopsis: Splcing pipeline
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
p7_splice_SpliceHits(P7_TOPHITS *tophits, P7_TOPHITS *seed_hits, P7_HMM *hmm, P7_OPROFILE *om, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, ESL_GETOPTS *go, ESL_GENCODE *gcode, ESL_SQFILE *seq_file, int64_t db_nuc_cnt)
{

  int                i;
  int                ncpus;
  int                infocnt;
  int                revcomp;
  int                seq_min, seq_max;
  ESL_SQ             *ali_seq;
  SPLICE_WORKER_INFO *info;
  int                 status;

  ali_seq = NULL;
  info    = NULL; 
 
  printf("\nQuery %s LENG %d\n",  gm->name, gm->M);
  fflush(stdout);


  /* Get the number of threads */
  ncpus = 0;
#ifdef HMMER_THREADS
  ncpus = ESL_MIN(esl_opt_GetInteger(go, "--cpu"), esl_threads_GetCPUCount());
#endif

//TODO
//ncpus = 4;
  /* Intialize data for threads */
  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(*info) * infocnt);
  for (i = 0; i < infocnt; ++i)
  {
	info[i].hmm            = p7_hmm_Clone(hmm);
    info[i].om             = p7_oprofile_Clone(om);
    info[i].gm             = p7_profile_Clone(gm);
    info[i].gm_fs          = p7_profile_fs_Clone(gm_fs);
	info[i].pli            = p7_splicepipeline_Create(go, gm->M, gm->M * 3);
	info[i].pli->bg        = p7_bg_Create(om->abc);
	if (info[i].pli->do_biasfilter)
	  p7_bg_SetFilter(info[i].pli->bg, om->M, om->compo);
	info[i].tophits        = tophits;
    info[i].gcode          = gcode;
	info[i].seq_file       = seq_file;
	info[i].db_nuc_cnt     = db_nuc_cnt;
  }

  for (i = 0; i < infocnt; ++i)
    info[i].seeds = seed_hits;

  /* Begin splicing */
#ifdef HMMER_THREADS
  if (ncpus > 0) thread_loop(info, tophits, seed_hits, infocnt); 
  else
#endif  /*HMMER_THREADS*/
    serial_loop (info, tophits, seed_hits); 
  
  /* Clean up */
  for (i = 0; i < infocnt; ++i)
  {
 	p7_hmm_Destroy(info[i].hmm);
	p7_oprofile_Destroy(info[i].om);
    p7_profile_Destroy(info[i].gm);
    p7_profile_fs_Destroy(info[i].gm_fs);
    p7_splicepipeline_Destroy(info[i].pli);
  }
  if(info           != NULL) free(info);

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
      ali_seq = p7_splice_GetSubSequence(seq_file, tophits->hit[i]->name, seq_min, seq_max, revcomp, NULL);
      tophits->hit[i]->dcl->ad = p7_alidisplay_fs_Create(tophits->hit[i]->dcl->tr, 0, gm_fs, ali_seq, gcode);
      tophits->hit[i]->dcl->ad->exon_cnt = 1;
      tophits->hit[i]->dcl->ad->sqfrom = tophits->hit[i]->dcl->iali;
      tophits->hit[i]->dcl->ad->sqto   = tophits->hit[i]->dcl->jali;  
      esl_sq_Destroy(ali_seq);
    }
  }

  return eslOK;

  ERROR:
	if(info != NULL) free(info);
    return status;
}



int  
serial_loop (SPLICE_WORKER_INFO *info, P7_TOPHITS *tophits, P7_TOPHITS *seed_hits) 
{

  int              g, h;
  int              num_graphs;
  int              revcomp, curr_revcomp;
  int64_t          seqidx, curr_seqidx;
  SPLICE_GRAPH    **graphs;
  SPLICE_GRAPH    *graph;
  int              status;

  graphs           = NULL;
  graph            = NULL;

  /*Find the number of unique sequences/strands with spliceable hits - this is the number of graphs to build*/
  num_graphs = 0;

  seqidx  = -1;
  revcomp = -1;

  for(h = 0; h < tophits->N; h++) {
	curr_seqidx = tophits->hit[h]->seqidx;
    curr_revcomp = 1;
	if (tophits->hit[h]->dcl->iali < tophits->hit[h]->dcl->jali)
      curr_revcomp = 0;

    if(curr_seqidx != seqidx || curr_revcomp != revcomp) {
      if(!(tophits->hit[h]->flags & p7_IS_DUPLICATE)) {
        if((tophits->hit[h]->flags & p7_IS_REPORTED) || exp(tophits->hit[h]->sum_lnP) < info->pli->S2) {

	      num_graphs++;
	      seqidx  = curr_seqidx; 
	      revcomp = curr_revcomp;
	    }
	  }
	}
  }

  if(num_graphs == 0) return eslOK;

  ESL_ALLOC(graphs, sizeof(SPLICE_GRAPH*) * num_graphs);

  /* Set the seqidx and revcomp of all graphs */
  seqidx  = -1;
  revcomp = -1;
  g = 0;
  for(h = 0; h < tophits->N; h++) {
    curr_seqidx = tophits->hit[h]->seqidx;
    curr_revcomp = 1;
    if (tophits->hit[h]->dcl->iali < tophits->hit[h]->dcl->jali)
      curr_revcomp = 0;

    if(curr_seqidx != seqidx || curr_revcomp != revcomp) {
      if(!(tophits->hit[h]->flags & p7_IS_DUPLICATE)) {
        if((tophits->hit[h]->flags & p7_IS_REPORTED) || exp(tophits->hit[h]->sum_lnP) < info->pli->S2) {
          graphs[g] = p7_splicegraph_Create();
		  graphs[g]->seqidx  = curr_seqidx;
		  graphs[g]->revcomp = curr_revcomp;
          g++;
          seqidx  = curr_seqidx;
          revcomp = curr_revcomp;
        }
      }
    }
  }

  /* Loop through all graphs and splice */
  for(g = 0; g < num_graphs; g++) {

	graph = graphs[g];

    /* Add all BATH hits from the correct seqidx and strand as nodes to graph */
    p7_splice_AddOriginals(info, graph, tophits);
    p7_splice_AddSeeds(graph, seed_hits);

    info->graph     = graph;
	info->thread_id = -1;

    p7_splice_SpliceGraph(info);

	p7_splicegraph_Destroy(graph);
    
  }
 
  if(graphs != NULL) free(graphs);

  return eslOK;

  ERROR:
    if(graphs != NULL)           free(graphs); 
    return status;
}

#ifdef HMMER_THREADS
int  
thread_loop (SPLICE_WORKER_INFO *info, P7_TOPHITS *tophits, P7_TOPHITS *seed_hits, int infocnt) 
{
  int              i, g, h;
  int              num_graphs;
  int              revcomp, curr_revcomp;
  int              graph_idx;
  int64_t          seqidx, curr_seqidx;
  SPLICE_GRAPH    **graphs;
  SPLICE_GRAPH    *graph;
  pthread_t       *threads;
  pthread_mutex_t  mutex;
  int              status;

  graphs           = NULL;
  graph            = NULL;
  threads          = NULL;

  /*Find the number of unique sequences/strands with spliceable hits - this is the number of graphs to build*/
  num_graphs = 0;

  seqidx  = -1;
  revcomp = -1;

  for(h = 0; h < tophits->N; h++) {
	curr_seqidx = tophits->hit[h]->seqidx;
    curr_revcomp = 1;
	if (tophits->hit[h]->dcl->iali < tophits->hit[h]->dcl->jali)
      curr_revcomp = 0;

    if(curr_seqidx != seqidx || curr_revcomp != revcomp) {
      if(!(tophits->hit[h]->flags & p7_IS_DUPLICATE)) {
        if((tophits->hit[h]->flags & p7_IS_REPORTED) || exp(tophits->hit[h]->sum_lnP) < info->pli->S2) {

	      num_graphs++;
	      seqidx  = curr_seqidx; 
	      revcomp = curr_revcomp;
	    }
	  }
	}
  }

  if(num_graphs == 0) return eslOK;

  ESL_ALLOC(graphs, sizeof(SPLICE_GRAPH*) * num_graphs);

  /* Set the seqidx and revcomp of all graphs */
  seqidx  = -1;
  revcomp = -1;
  g = 0;
  for(h = 0; h < tophits->N; h++) {
    curr_seqidx = tophits->hit[h]->seqidx;
    curr_revcomp = 1;
    if (tophits->hit[h]->dcl->iali < tophits->hit[h]->dcl->jali)
      curr_revcomp = 0;

    if(curr_seqidx != seqidx || curr_revcomp != revcomp) {
      if(!(tophits->hit[h]->flags & p7_IS_DUPLICATE)) {
        if((tophits->hit[h]->flags & p7_IS_REPORTED) || exp(tophits->hit[h]->sum_lnP) < info->pli->S2) {
          graphs[g] = p7_splicegraph_Create();
		  graphs[g]->seqidx  = curr_seqidx;
		  graphs[g]->revcomp = curr_revcomp;
          g++;
          seqidx  = curr_seqidx;
          revcomp = curr_revcomp;
        }
      }
    }
  }

  /* Loop though graphs and add BATH hits */
  for(g = 0; g < num_graphs; g++) {

	graph = graphs[g];

    p7_splice_AddOriginals(info, graph, tophits);
    p7_splice_AddSeeds(graph, seed_hits);

  }
 
  /* Initialize mutex */
  if (pthread_mutex_init(&mutex, NULL)) { status = eslFAIL; goto ERROR; }

  /* Allocate space for threads */
  ESL_ALLOC(threads, sizeof(pthread_t) * infocnt);

  /* Add data needed for multithreading and spin up threads to splice graphs */
  graph_idx = 0;
  for(i = 0; i < infocnt; i++) {
    info[i].mutex = &mutex;
    info[i].graph_idx  = &graph_idx;
	info[i].num_graphs = num_graphs; 
    info[i].graphs     = graphs;
	info[i].thread_id  = i;
	if (pthread_create(&threads[i], NULL, &splice_thread, &info[i])) { status = eslFAIL; goto ERROR; }
  }
  
  /* Wait for threads to complete */
  for (i = 0; i < infocnt; i++) {
    pthread_join(threads[i], NULL);
  }

  /* Clean up */
  pthread_mutex_destroy(&mutex);
  if (threads != NULL) free(threads); 
  for(g = 0; g < num_graphs; g++) p7_splicegraph_Destroy(graphs[g]);
  if(graphs != NULL) free(graphs);

  return eslOK;

  ERROR:
    if(threads != NULL) free(threads);
    if(graphs != NULL)  free(graphs); 
    pthread_mutex_destroy(&mutex);
    return status;

}

void* 
splice_thread(void *arg)
{

  SPLICE_WORKER_INFO *info = (SPLICE_WORKER_INFO *)arg;
  int graph_idx;

  impl_Init();
  
  while (1) {
    /* Get next graph to process */
    pthread_mutex_lock(info->mutex);
    graph_idx = (*info->graph_idx)++;
    pthread_mutex_unlock(info->mutex);
    
    /* Exit if no more graphs to process */
    if (graph_idx >= info->num_graphs) break;

    info->graph =  info->graphs[graph_idx];
    p7_splice_SpliceGraph(info); 
 
  }
  return NULL;
}
#endif  /*HMMER_THREADS*/


/*  Function: p7_splice_AddOriginals
 *  Synopsis: Add BATH hits to graph
 *
 *  Purpose : Find all hits in <tophits> that match the seqidx 
 *            and revcomp of the <graph> and add them as nodes.
 *
 * Returns:   <eslOK>.
 *
 */
int 
p7_splice_AddOriginals(SPLICE_WORKER_INFO *info, SPLICE_GRAPH *graph, const P7_TOPHITS *tophits)
{

  int     i;
  int     hit_cnt;
  P7_HIT *curr_hit;

  /* Get a count of the hits that will be added to the graph */
  hit_cnt = 0;
  for (i = 0; i < tophits->N; i++) {
    curr_hit = tophits->hit[i];
    if (curr_hit->seqidx != graph->seqidx) continue;
    if (graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
    if ((!graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;
    
    if ((tophits->hit[i]->flags & p7_IS_DUPLICATE)) continue;
    if(!(tophits->hit[i]->flags & p7_IS_REPORTED) && exp(tophits->hit[i]->sum_lnP) >= info->pli->S2) continue;

    hit_cnt++;
  }
  
  p7_splicegraph_CreateNodes(graph, hit_cnt);    
  
  /*Add all hits from current sequence and strand to grapgh*/
  for (i = 0; i < tophits->N; i++) {

    curr_hit = tophits->hit[i];

    if (curr_hit->seqidx != graph->seqidx) continue;

    if (graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
    if ((!graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;

    if ((curr_hit->flags & p7_IS_DUPLICATE)) continue;
	if(!(curr_hit->flags & p7_IS_REPORTED) && exp(curr_hit->sum_lnP) >= info->pli->S2) continue;

    if(curr_hit->flags & p7_IS_REPORTED) graph->reportable[graph->th->N] = TRUE;
    else                                 graph->reportable[graph->th->N] = FALSE;

    graph->node_in_graph[graph->th->N] = TRUE;
    graph->orig_hit_idx[graph->th->N]  = i;
 
    p7_splicegraph_AddNode(graph, curr_hit); 

    graph->split_orig_id[graph->th->N-1] = graph->th->N-1;
  }

  graph->orig_N = graph->recover_N = graph->th->N;
  graph->seqname = graph->th->hit[0]->name;

  return eslOK;
 
}

/* Only add seed that are upstream of one orignal hit and downstream fo another */
int 
p7_splice_AddSeeds(SPLICE_GRAPH *graph, const P7_TOPHITS *seed_hits)
{

  int     i;
  int     h1, h2;
  int     gap_len;
  P7_HIT *curr_hit;
  P7_TOPHITS *th;

  th = graph->th;

  if(graph->orig_N < 2) return eslOK;

  for (i = 0; i < seed_hits->N; i++) {

    curr_hit = &(seed_hits->unsrt[i]);

    /* Is the seed hit on the same sequence and strand as the graph */
    if (curr_hit->seqidx != graph->seqidx) continue;
    if (graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
    if ((!graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;

    for(h1 = 0; h1 < graph->orig_N; h1++) {
      /* Is the seed hit upstream of an original hit */
      if(p7_splice_HitUpstream(curr_hit->dcl, th->hit[h1]->dcl, graph->revcomp)) {
       //printf("seed %d upstream of hit %d\n", i+1, h1+1); 
        if(graph->revcomp) gap_len = curr_hit->dcl->jali - th->hit[h1]->dcl->iali - 1;
        else               gap_len = th->hit[h1]->dcl->iali - curr_hit->dcl->jali - 1;
        if(gap_len > MAX_INTRON_LENG) continue;
       
        for(h2 = 0; h2 < graph->orig_N; h2++) {
          if(h2 == h1) continue;

          /* Is the seed hit downstream of an original hit */
          if(p7_splice_HitUpstream(th->hit[h2]->dcl, curr_hit->dcl, graph->revcomp)) {
              
             if(graph->revcomp) gap_len = th->hit[h2]->dcl->jali - curr_hit->dcl->iali - 1;
             else               gap_len = curr_hit->dcl->iali - th->hit[h2]->dcl->jali - 1;
             if(gap_len > MAX_INTRON_LENG) continue;

             /* Repurpose domain "is_included" for seed hits added to graph */
             curr_hit->dcl->is_included = TRUE;
             p7_splicegraph_AddNode(graph, curr_hit);
             h1 = graph->orig_N;
             h2 = graph->orig_N;
          }
        }
      }
    }
  }
  graph->recover_N = graph->num_nodes; 
  return eslOK;
 
}


/*  Function: p7_splice_SpliceGraph
 *  Synopsis: Main splicing function
 *
 *  Purpose : Given a splice graph with BATH hits from a single target 
 *            sequence and strand, add recovered and missing exons, 
 *            find spice boundries and align the spliced sequence
 *
 * Returns:   <eslOK> on success. 
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int 
p7_splice_SpliceGraph(SPLICE_WORKER_INFO *info) 
{
  
  int                s;
  int                seq_min, seq_max;
  int64_t            path_min, path_max;
  SPLICE_PATH2       *final_path;
  SPLICE_PATH2       *orig_path;
  ESL_SQ            *path_seq;
  SPLICE_GRAPH      *graph;
  SPLICE_PIPELINE   *pli;
  P7_TOPHITS        *tophits;
  P7_OPROFILE       *om;
  P7_PROFILE        *gm;
  P7_FS_PROFILE     *gm_fs;
  ESL_GENCODE       *gcode;
  ESL_SQFILE        *seq_file;

  path_seq         = NULL;

  graph      = info->graph;
  pli        = info->pli;
  tophits    = info->tophits;
  om         = info->om;
  gm         = info->gm;
  gm_fs      = info->gm_fs;
  gcode      = info->gcode;
  seq_file   = info->seq_file;

  printf("\nQuery %s Target %s strand %c seqidx %ld\n", gm->name, graph->seqname, (graph->revcomp ? '-' : '+'), graph->seqidx);
  fflush(stdout);

//if(info->thread_id >= 0) pthread_mutex_lock(info->mutex);
//printf("RECOVER\n");
//p7_splicegraph_DumpHits(stdout, graph);
//fflush(stdout);
//if(info->thread_id >= 0) pthread_mutex_unlock(info->mutex);

  /* Create edges between original and recovered nodes */
  p7_splice_CreateUnsplicedEdges(graph);

  /* Build paths from orignal hit nodes and edge so that every node appears in one and only one path */
  orig_path = p7_splicepath_GetBestPath_Unspliced2(graph);


  while(orig_path != NULL) {

    for(s = 0; s < orig_path->path_len; s++) {
      graph->node_in_graph[orig_path->node_id[s]] = FALSE;
    }
   
    /* Break edges that overlap the path so that paths do not intertwine */ 
    path_min = ESL_MIN(orig_path->iali[0], orig_path->jali[orig_path->path_len-1]);
    path_max = ESL_MAX(orig_path->iali[0], orig_path->jali[orig_path->path_len-1]);
    p7_splice_EnforceRangeBounds(graph, path_min, path_max);

    p7_splice_ExtendPath2(info->seeds, orig_path, graph);

   
    seq_min = ESL_MIN(orig_path->iali[0], orig_path->jali[orig_path->path_len-1]) - ALIGNMENT_EXT*3;
    seq_max = ESL_MAX(orig_path->iali[0], orig_path->jali[orig_path->path_len-1]) + ALIGNMENT_EXT*3;
    path_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, orig_path->revcomp, info);
    
//if(info->thread_id >= 0) pthread_mutex_lock(info->mutex);
p7_splicepath_Dump2(stdout,orig_path);
//if(info->thread_id >= 0) pthread_mutex_unlock(info->mutex);

    final_path = p7_splice_FindExons2(info, orig_path, path_seq);
p7_splicepath_Dump2(stdout,final_path);
       //if(!path->frameshift)
      if(final_path->path_len > 1)
        p7_splice_AlignPath2(graph, final_path, pli, tophits, om, gm, gcode, path_seq, info->db_nuc_cnt, gm_fs->fs, info);

      //if(path->frameshift)
        //p7_splice_AlignFrameshiftPath (graph, path, pli, tophits, gm_fs, gcode, path_seq, info->db_nuc_cnt, info);

      esl_sq_Destroy(path_seq);
      p7_splicepath_Destroy2(orig_path);
      p7_splicepath_Destroy2(final_path);
      p7_splicepipeline_Reuse(pli);

      orig_path = p7_splicepath_GetBestPath_Unspliced2(graph);
  }
  
  return eslOK;

}




/*  Function: p7_splice_ExtendPath
 *  Synopsis: Add seed hits to beginning or end of path 
 *
 *  Purpose : Recover potential exons that were discarded by
 *            the BATH Viterbi filter but saved in <ifo->sh>
 *            and add them to beginning or end of splice path
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_splice_ExtendPath(P7_TOPHITS *seed_hits, SPLICE_PATH *path, SPLICE_GRAPH *graph)
{
  int i,s;
  int gap_len;
  P7_HIT *first_hit;
  P7_HIT *last_hit; 
  P7_HIT *curr_hit;
  SPLICE_GRAPH *tmp_graph;
  SPLICE_PATH  *tmp_path;
  
  first_hit = path->hits[0];
  last_hit  = path->hits[path->path_len-1];

  /* Only extend from begining if first hit is an original hit */
  if(path->node_id[0] < graph->orig_N) {
    tmp_graph = p7_splicegraph_Create();
    tmp_graph->seqidx  = graph->seqidx;
    tmp_graph->revcomp = graph->revcomp;
    tmp_graph->seqname = graph->seqname; 
   
    p7_splicegraph_CreateNodes(tmp_graph, 1);  
    p7_splicegraph_AddNode(tmp_graph, first_hit);
  
    tmp_graph->reportable[0]    = TRUE;
    tmp_graph->orig_hit_idx[0]  = path->node_id[0];
    tmp_graph->split_orig_id[0] = 0;
  
    tmp_graph->orig_N = 1;
    
    for(i = 0; i < seed_hits->N; i++) {
      
      curr_hit = &(seed_hits->unsrt[i]);
      
      /* skip seeds already added to a graph */
      if(curr_hit->dcl->is_included) continue;
  
      /* Is the seed hit on the same sequence and strand as the tmp_graph */
      if (curr_hit->seqidx != tmp_graph->seqidx) continue;
      if (tmp_graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
      if ((!tmp_graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;
  
      /*If the seed is upstream of the first hit, add it to the tmp graph */
      if(p7_splice_HitUpstream(curr_hit->dcl, first_hit->dcl, tmp_graph->revcomp)) { 
        if(graph->revcomp) gap_len = curr_hit->dcl->jali  - first_hit->dcl->iali - 1;
        else               gap_len = first_hit->dcl->iali - curr_hit->dcl->jali - 1;
        if(gap_len > MAX_INTRON_EXT || gap_len < ( first_hit->dcl->ihmm - curr_hit->dcl->jhmm - 1) * 3) continue;

        p7_splicegraph_AddNode(tmp_graph, curr_hit);
        tmp_graph->reportable[tmp_graph->num_nodes-1]    = FALSE;
        tmp_graph->orig_hit_idx[tmp_graph->num_nodes-1]  = -1;
        tmp_graph->split_orig_id[tmp_graph->num_nodes-1] = -1;  
      }
          
    }
    tmp_graph->recover_N = tmp_graph->num_nodes;
 
    p7_splice_CreateUnsplicedEdges(tmp_graph);
    tmp_path = p7_splicepath_GetBestPath_Extension(graph,tmp_graph); 
//printf("tmp graph %d\n", tmp_graph->num_nodes);
//p7_splicegraph_DumpHits(stdout, tmp_graph);    
//printf("TMP PATH\n");
//p7_splicepath_Dump(stdout,tmp_path);
//fflush(stdout);

    /* Set the most upstream hit in the tmp_path to start at the minimum 
     * hmm for all hits in tmp_path and add it to the orginal path  */
    for(s = tmp_path->path_len-2; s >= 0; s--) {
        
      tmp_path->hits[s]->dcl->is_included = TRUE;
      p7_splicepath_Insert(path, tmp_path->hits[s], tmp_path->edge_scores[s], 0); 
      p7_splicegraph_AddNode(graph, curr_hit);
      path->node_id[0] = graph->num_nodes-1; 
    }

    p7_splicepath_Destroy(tmp_path);
    p7_splicegraph_Destroy(tmp_graph);
  }

  /* Only extend from end if last hit is an original hit */
  if(path->node_id[path->path_len-1] < graph->orig_N) {
    tmp_graph = p7_splicegraph_Create();
    tmp_graph->seqidx  = graph->seqidx;
    tmp_graph->revcomp = graph->revcomp;
    tmp_graph->seqname = graph->seqname;

    p7_splicegraph_CreateNodes(tmp_graph, 1);
    p7_splicegraph_AddNode(tmp_graph, last_hit);
 
    tmp_graph->reportable[0]    = TRUE;
    tmp_graph->orig_hit_idx[0]  = path->node_id[path->path_len-1];
    tmp_graph->split_orig_id[0] = 0;

    tmp_graph->orig_N = 1;

    for(i = 0; i < seed_hits->N; i++) {

      curr_hit = &(seed_hits->unsrt[i]);

      /* skip seeds already added to a graph */
      if(curr_hit->dcl->is_included) continue;

      /* Is the seed hit on the same sequence and strand as the tmp_graph */
      if (curr_hit->seqidx != tmp_graph->seqidx) continue;
      if (tmp_graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
      if ((!tmp_graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;

      /*If the seed is upstream of the first hit, add it to the tmp graph */
      if(p7_splice_HitUpstream(last_hit->dcl, curr_hit->dcl, tmp_graph->revcomp)) {
        if(graph->revcomp) gap_len = last_hit->dcl->jali - curr_hit->dcl->iali - 1;
        else               gap_len = curr_hit->dcl->iali - last_hit->dcl->jali - 1;
        if(gap_len > MAX_INTRON_EXT || gap_len < ( first_hit->dcl->ihmm - curr_hit->dcl->jhmm - 1) * 3) continue;

        p7_splicegraph_AddNode(tmp_graph, curr_hit);
        tmp_graph->reportable[tmp_graph->num_nodes-1]    = FALSE;
        tmp_graph->orig_hit_idx[tmp_graph->num_nodes-1]  = -1;
        tmp_graph->split_orig_id[tmp_graph->num_nodes-1] = -1;
      }
    }   
    tmp_graph->recover_N = tmp_graph->num_nodes;

    p7_splice_CreateUnsplicedEdges(tmp_graph);
    tmp_path = p7_splicepath_GetBestPath_Extension(graph, tmp_graph);
//printf("TMP PATH\n");
//p7_splicepath_Dump(stdout,tmp_path);
    /* Set the most downstream hit in the tmp_path to end at the maximum  
     * hmm for all hits in tmp_path and add it to the orginal path  */
    /* Add any seed hits from tmp_path to real path */
    for(s = 1; s < tmp_path->path_len; s++) { 
      tmp_path->hits[s]->dcl->is_included = TRUE;
      p7_splicepath_Insert(path, tmp_path->hits[s], tmp_path->edge_scores[s], path->path_len);
      p7_splicegraph_AddNode(graph, curr_hit);
      path->node_id[path->path_len-1] = graph->num_nodes-1;
    }
    
    p7_splicepath_Destroy(tmp_path);
    p7_splicegraph_Destroy(tmp_graph);
  }
  return eslOK; 

}

int
p7_splice_ExtendPath2(P7_TOPHITS *seed_hits, SPLICE_PATH2 *path, SPLICE_GRAPH *graph)
{
  int i,s;
  int gap_len;
  P7_HIT *first_hit;
  P7_HIT *last_hit; 
  P7_HIT *curr_hit;
  SPLICE_GRAPH *tmp_graph;
  SPLICE_PATH  *tmp_path;
  
  first_hit = graph->th->hit[path->node_id[0]];
  last_hit  = graph->th->hit[path->node_id[path->path_len-1]];

  /* Only extend from begining if first hit is an original hit */
  if(path->node_id[0] < graph->orig_N) {
    tmp_graph = p7_splicegraph_Create();
    tmp_graph->seqidx  = graph->seqidx;
    tmp_graph->revcomp = graph->revcomp;
    tmp_graph->seqname = graph->seqname; 
   
    p7_splicegraph_CreateNodes(tmp_graph, 1);  
    p7_splicegraph_AddNode(tmp_graph, first_hit);
  
    tmp_graph->reportable[0]    = TRUE;
    tmp_graph->orig_hit_idx[0]  = path->node_id[0];
    tmp_graph->split_orig_id[0] = 0;
  
    tmp_graph->orig_N = 1;
    
    for(i = 0; i < seed_hits->N; i++) {
      
      curr_hit = &(seed_hits->unsrt[i]);
      
      /* skip seeds already added to a graph */
      if(curr_hit->dcl->is_included) continue;
  
      /* Is the seed hit on the same sequence and strand as the tmp_graph */
      if (curr_hit->seqidx != tmp_graph->seqidx) continue;
      if (tmp_graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
      if ((!tmp_graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;
  
      /*If the seed is upstream of the first hit, add it to the tmp graph */
      if(p7_splice_HitUpstream(curr_hit->dcl, first_hit->dcl, tmp_graph->revcomp)) { 
        if(graph->revcomp) gap_len = curr_hit->dcl->jali  - first_hit->dcl->iali - 1;
        else               gap_len = first_hit->dcl->iali - curr_hit->dcl->jali - 1;
        if(gap_len > MAX_INTRON_EXT || gap_len < ( first_hit->dcl->ihmm - curr_hit->dcl->jhmm - 1) * 3) continue;

        p7_splicegraph_AddNode(tmp_graph, curr_hit);
        tmp_graph->reportable[tmp_graph->num_nodes-1]    = FALSE;
        tmp_graph->orig_hit_idx[tmp_graph->num_nodes-1]  = -1;
        tmp_graph->split_orig_id[tmp_graph->num_nodes-1] = -1;  
      }
          
    }
    tmp_graph->recover_N = tmp_graph->num_nodes;
 
    p7_splice_CreateUnsplicedEdges(tmp_graph);
    tmp_path = p7_splicepath_GetBestPath_Extension(graph,tmp_graph); 

    /* Set the most upstream hit in the tmp_path to start at the minimum 
     * hmm for all hits in tmp_path and add it to the orginal path  */
    for(s = tmp_path->path_len-2; s >= 0; s--) {
        
      tmp_path->hits[s]->dcl->is_included = TRUE;
      p7_splicepath_Insert2(path, 0); 
      p7_splicegraph_AddNode(graph, curr_hit);
      path->node_id[0] = graph->num_nodes-1; 
      path->ihmm[0]    =  tmp_path->hits[s]->dcl->ihmm;
      path->jhmm[0]    =  tmp_path->hits[s]->dcl->jhmm;
      path->iali[0]    =  tmp_path->hits[s]->dcl->iali;
      path->jali[0]    =  tmp_path->hits[s]->dcl->jali;
    }

    p7_splicepath_Destroy(tmp_path);
    p7_splicegraph_Destroy(tmp_graph);
  }

  /* Only extend from end if last hit is an original hit */
  if(path->node_id[path->path_len-1] < graph->orig_N) {
    tmp_graph = p7_splicegraph_Create();
    tmp_graph->seqidx  = graph->seqidx;
    tmp_graph->revcomp = graph->revcomp;
    tmp_graph->seqname = graph->seqname;

    p7_splicegraph_CreateNodes(tmp_graph, 1);
    p7_splicegraph_AddNode(tmp_graph, last_hit);
 
    tmp_graph->reportable[0]    = TRUE;
    tmp_graph->orig_hit_idx[0]  = path->node_id[path->path_len-1];
    tmp_graph->split_orig_id[0] = 0;

    tmp_graph->orig_N = 1;

    for(i = 0; i < seed_hits->N; i++) {

      curr_hit = &(seed_hits->unsrt[i]);

      /* skip seeds already added to a graph */
      if(curr_hit->dcl->is_included) continue;

      /* Is the seed hit on the same sequence and strand as the tmp_graph */
      if (curr_hit->seqidx != tmp_graph->seqidx) continue;
      if (tmp_graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
      if ((!tmp_graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;

      /*If the seed is upstream of the first hit, add it to the tmp graph */
      if(p7_splice_HitUpstream(last_hit->dcl, curr_hit->dcl, tmp_graph->revcomp)) {
        if(graph->revcomp) gap_len = last_hit->dcl->jali - curr_hit->dcl->iali - 1;
        else               gap_len = curr_hit->dcl->iali - last_hit->dcl->jali - 1;
        if(gap_len > MAX_INTRON_EXT || gap_len < ( first_hit->dcl->ihmm - curr_hit->dcl->jhmm - 1) * 3) continue;

        p7_splicegraph_AddNode(tmp_graph, curr_hit);
        tmp_graph->reportable[tmp_graph->num_nodes-1]    = FALSE;
        tmp_graph->orig_hit_idx[tmp_graph->num_nodes-1]  = -1;
        tmp_graph->split_orig_id[tmp_graph->num_nodes-1] = -1;
      }
    }   
    tmp_graph->recover_N = tmp_graph->num_nodes;

    p7_splice_CreateUnsplicedEdges(tmp_graph);
    tmp_path = p7_splicepath_GetBestPath_Extension(graph, tmp_graph);

    /* Set the most downstream hit in the tmp_path to end at the maximum  
     * hmm for all hits in tmp_path and add it to the orginal path  */
    /* Add any seed hits from tmp_path to real path */
    for(s = 1; s < tmp_path->path_len; s++) { 
      tmp_path->hits[s]->dcl->is_included = TRUE;
      p7_splicepath_Insert2(path, path->path_len);
      p7_splicegraph_AddNode(graph, curr_hit);
      path->node_id[path->path_len-1] = graph->num_nodes-1;
      path->ihmm[path->path_len-1]    =  tmp_path->hits[s]->dcl->ihmm;
      path->jhmm[path->path_len-1]    =  tmp_path->hits[s]->dcl->jhmm;
      path->iali[path->path_len-1]    =  tmp_path->hits[s]->dcl->iali;
      path->jali[path->path_len-1]    =  tmp_path->hits[s]->dcl->jali;
    }
    
    p7_splicepath_Destroy(tmp_path);
    p7_splicegraph_Destroy(tmp_graph);
  }
  return eslOK; 

}


/*  Function: p7_splice_CreateUnsplicedEdges
 *  Synopsis: Add unspliced edges
 *
 *  Purpose : Given a splice graph with nodes, add unspliced
 *            edges between any hits that are up/down stream of
 *            each other. Edge scores are zero unless the hits 
 *            overlap in amino positions, when the edge score 
 *            the minimum lost score to eleimate the overlap. 
 *
 * Returns:   <eslOK>.
 *
 */
int
p7_splice_CreateUnsplicedEdges(SPLICE_GRAPH *graph) 
{
  int up, down;
  int seq_gap_len;
  int amino_gap_len;
  P7_TOPHITS  *th;
  SPLICE_EDGE *edge;

  th = graph->th;

  for(up = 0; up < graph->num_nodes; up++) {
    
    for(down = 0; down < graph->num_nodes; down++) {
      if(up == down) continue;
      
      /* Are hits upstream/downstream */
      if(!p7_splice_HitUpstream(th->hit[up]->dcl, th->hit[down]->dcl, graph->revcomp)) continue;
      
      if(graph->revcomp)
        seq_gap_len = th->hit[up]->dcl->jali - th->hit[down]->dcl->iali - 1;
      else
        seq_gap_len = th->hit[down]->dcl->iali - th->hit[up]->dcl->jali - 1;        

      if(seq_gap_len > MAX_INTRON_LENG) continue;

      amino_gap_len = th->hit[down]->dcl->ihmm - th->hit[up]->dcl->jhmm - 1;
      if(seq_gap_len < amino_gap_len * 3) continue;     

      edge = p7_splicegraph_AddEdge(graph, up, down);
      edge->splice_score       = 0.;

      /* If hits overlap, find the minimum lost socre to remove the overlap */
      p7_spliceedge_AliScoreEdge(edge, th->hit[up]->dcl, th->hit[down]->dcl); 

      edge->upstream_spliced_amino_end     = th->hit[up]->dcl->jhmm;
      edge->downstream_spliced_amino_start = th->hit[down]->dcl->ihmm; 
      edge->upstream_spliced_nuc_end       = th->hit[up]->dcl->jali;
      edge->downstream_spliced_nuc_start   = th->hit[down]->dcl->iali;
 
    }
  } 
  
  return eslOK;

}

/*  Function: p7_splice_FindExons
 *  Synopsis: Add unspliced edges
 *
 *  Purpose : Given an unspliced path of original and recovered 
 *            nodes, align the path target sub-sequence to the 
 *            model using spliced viterbi and missing exons and 
 *            locate likely splice boundries.
 *
 * Returns:   <eslOK>.
 *
 */
int
p7_splice_FindExons(SPLICE_WORKER_INFO *info, SPLICE_PATH *path, ESL_SQ *path_seq)
{

  int e, s;
  int num_hits;
  int seq_start;
  int seq_len;
  int intron_len;
  int amino_gap;
  int             *remove_idx;
  int64_t         *nuc_index;
  P7_HMM          *sub_hmm;
  P7_HMM          *hmm;
  P7_FS_PROFILE   *gm_fs;
  ESL_GENCODE     *gcode;
  ESL_SQ          *sub_seq;
  SPLICE_GRAPH    *graph;
  SPLICE_PIPELINE *pli;
  P7_HIT          **exons;

  graph = info->graph;
  pli   = info->pli;
  gcode = info->gcode;
  gm_fs = info->gm_fs;
  hmm   = info->hmm;

  remove_idx = NULL;
  nuc_index  = NULL;

  /* If there is only one hit in the path realign it with non-frameshift
   * spliced viterbi to see if it contains more than one exon. If the orginal
   * hit has frameshifts realign a second time with frameshift spliced viterbi */
  if(path->path_len == 1) {

    sub_hmm     = p7_splice_GetSubHMM(hmm, path->hits[0]->dcl->ihmm, path->hits[0]->dcl->jhmm);
    sub_hmm->fs = 0.;

    /* Get all hits (exons) in path */
     exons = p7_splice_AlignExons(pli, sub_hmm, gm_fs, pli->bg, path_seq, gcode, graph->revcomp, path->hits[0]->dcl->ihmm, 0, remove_idx, &num_hits);
    //printf("  %4s %10s %10s %5s %5s\n", "exon", "iali", "jali", "ihmm", "jhmm");
    for(e = 0; e < num_hits; e++) {
      if(graph->revcomp) {
        exons[e]->dcl->iali = path_seq->n - exons[e]->dcl->iali + path_seq->end;
        exons[e]->dcl->jali = path_seq->n - exons[e]->dcl->jali + path_seq->end;
      }
      else {
        exons[e]->dcl->iali = path_seq->start + exons[e]->dcl->iali - 1;
        exons[e]->dcl->jali = path_seq->start + exons[e]->dcl->jali - 1;
      }
      p7_splicegraph_AddNode(graph, exons[e]);
      graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[path->node_id[0]];
      graph->split_orig_id[graph->num_nodes-1] = path->node_id[0];
    }

    free(exons);

    p7_hmm_Destroy(sub_hmm);

    return eslOK;
  }


    /*For each neighboring pair of nodes in the path, create a sub sequence consisting of the 
   * nucs from the start of the upstream hit to the end of the downstream hit. If the sequence 
   * gap between the hits is large (over MAX_INTRON_INCL) include one the first and last 
   * MAX_INTRON_INCL/2 nucleotides.  If the sequence gap is trimmed and there is a gap in the 
   * aminos include enough degenerate nucleotides to cover the amino gap so they have a neutral 
   * place to align.  Run this sequence against the hmm positions from the start of the upstream 
   * hit to the end of the downstream hit.  After the first pair, use the coordinates of the 
   * returned exons from the pervious hit to determine the start coordinates of the upstream hit*/

  for(s = 1; s < path->path_len; s++) {

    sub_seq = p7_splice_GetSplicedSequence(path, path_seq, s-1, s, TRUE, &remove_idx, &nuc_index);
    
    intron_len = llabs(path->hits[s]->dcl->iali - path->hits[s-1]->dcl->jali) - 1;
    amino_gap =  ESL_MAX(0, path->hits[s]->dcl->ihmm - path->hits[s-1]->dcl->jhmm - 1);

    /* Crete sub hmm */
    sub_hmm     = p7_splice_GetSubHMM(hmm, path->hits[s-1]->dcl->ihmm, path->hits[s]->dcl->jhmm);
    sub_hmm->fs = 0.;

    /* Align sub seq to sub hmm and get exons */
    exons = p7_splice_AlignExons(pli, sub_hmm, gm_fs, pli->bg, sub_seq, gcode, graph->revcomp, path->hits[s-1]->dcl->ihmm, 1, remove_idx, &num_hits);  

    
    for(e = 0; e < num_hits; e++) {
      exons[e]->dcl->iali = nuc_index[exons[e]->dcl->iali];
      exons[e]->dcl->jali = nuc_index[exons[e]->dcl->jali];
     //printf("s-1 %d s %d e %d iali %d jali %d ihmm %d jhm %d\n", path->node_id[s-1]+1, path->node_id[s]+1, e, exons[e]->dcl->iali, exons[e]->dcl->jali, exons[e]->dcl->ihmm, exons[e]->dcl->jhmm); 
    }

    /* If there is only one exon it is likely that one of the nodes is a false posiitive and the 
     * resulting exon will be skewd by having to cover the hmm postions of both nodes. If the 
     * exon covers both hits we keep it. If the Exon covers only the first hit we remove the 
     * second node from the path, discard the exon and start the next search at s-1 again. If 
     * the Exon only covers the seecond hit we discard the exon and start the next seach at s */ 

    
   if(num_hits == 1) {
      
      if(seq_overlap(exons[0]->dcl, path->hits[s-1]->dcl, graph->revcomp) >= 0.0 &&
          hmm_overlap(exons[0]->dcl, path->hits[s-1]->dcl) >= 0.0) {
  
     
        if(seq_overlap(exons[0]->dcl, path->hits[s]->dcl, graph->revcomp) >= 0.0 &&
           hmm_overlap(exons[0]->dcl, path->hits[s]->dcl) >= 0.0) {
          
          p7_splicegraph_AddNode(graph, exons[0]);
          if(path->node_id[s-1] < graph->orig_N  && path->node_id[s-1] != -1) {
            graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[path->node_id[s-1]];
            graph->split_orig_id[graph->num_nodes-1] = path->node_id[s-1];        
          }
  
          if(path->node_id[s] < graph->orig_N && path->node_id[s] != -1) {
            graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[path->node_id[s]];
            graph->split_orig_id[graph->num_nodes-1] = path->node_id[s];
          }
  
    
          path->hits[s] = exons[0]; 
        }
   
        else {
          
          p7_splicepath_Remove(path, s);
          s--; 

          p7_trace_splice_Destroy(exons[0]->dcl->tr);
          free(exons[0]->dcl->scores_per_pos);
          p7_hit_Destroy(exons[0]);
        }
      } 
  
      else {
        p7_splicepath_Remove(path, s-1);
        s = ESL_MAX(0, s-2);
        p7_trace_splice_Destroy(exons[0]->dcl->tr);
        free(exons[0]->dcl->scores_per_pos);
        p7_hit_Destroy(exons[0]);
      }
    }   
    else if (num_hits == 2){

      for(e = 0; e < num_hits; e++) {
            p7_splicegraph_AddNode(graph, exons[e]);
            if(path->node_id[s-1] < graph->orig_N  && path->node_id[s-1] != -1) {
              if(seq_overlap(exons[e]->dcl, path->hits[s-1]->dcl, graph->revcomp) >= 0.0 &&
                 hmm_overlap(exons[e]->dcl, path->hits[s-1]->dcl) >= 0.0) {
                graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[path->node_id[s-1]];
                graph->split_orig_id[graph->num_nodes-1] = path->node_id[s-1];
              }
            }
            if(path->node_id[s] < graph->orig_N  && path->node_id[s] != -1) {
              if(seq_overlap(exons[e]->dcl, path->hits[s]->dcl, graph->revcomp) >= 0.0 &&
                 hmm_overlap(exons[e]->dcl, path->hits[s]->dcl) >= 0.0) {
                graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[path->node_id[s]];
                graph->split_orig_id[graph->num_nodes-1] = path->node_id[s];
              }
            }
          }
          path->hits[s] = exons[num_hits-1];
/*
      if(seq_overlap(exons[0]->dcl, path->hits[s-1]->dcl, graph->revcomp) >= 0.0 &&
         hmm_overlap(exons[0]->dcl, path->hits[s-1]->dcl) >= 0.0) {
        if(seq_overlap(exons[1]->dcl, path->hits[s]->dcl, graph->revcomp) >= 0.0 &&
           hmm_overlap(exons[1]->dcl, path->hits[s]->dcl) >= 0.0) {  

          for(e = 0; e < num_hits; e++) {
            p7_splicegraph_AddNode(graph, exons[e]);
            if(path->node_id[s-1] < graph->orig_N  && path->node_id[s-1] != -1) {
              if(seq_overlap(exons[e]->dcl, path->hits[s-1]->dcl, graph->revcomp) >= 0.0 &&
                 hmm_overlap(exons[e]->dcl, path->hits[s-1]->dcl) >= 0.0) {
                graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[path->node_id[s-1]];
                graph->split_orig_id[graph->num_nodes-1] = path->node_id[s-1];
              }
            }
            if(path->node_id[s] < graph->orig_N  && path->node_id[s] != -1) {
              if(seq_overlap(exons[e]->dcl, path->hits[s]->dcl, graph->revcomp) >= 0.0 &&
                 hmm_overlap(exons[e]->dcl, path->hits[s]->dcl) >= 0.0) {
                graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[path->node_id[s]];
                graph->split_orig_id[graph->num_nodes-1] = path->node_id[s];
              }
            }
          }
          path->hits[s] = exons[num_hits-1];
        }
        else {
          p7_splicepath_Remove(path, s);
          s--;
          for(e = 0; e < num_hits; e++) {
            p7_trace_splice_Destroy(exons[e]->dcl->tr);
            free(exons[e]->dcl->scores_per_pos);
            p7_hit_Destroy(exons[e]);
          }
        }
      }
      else {
        p7_splicepath_Remove(path, s-1);
        s = ESL_MAX(0, s-2);
        for(e = 0; e < num_hits; e++) {
          p7_trace_splice_Destroy(exons[e]->dcl->tr);
          free(exons[e]->dcl->scores_per_pos);
          p7_hit_Destroy(exons[e]);
        }
      }
*/
    }
    else if (num_hits > 2) {

      if(intron_len <= MAX_INTRON_INCL + (amino_gap*3)) {
        for(e = 0; e < num_hits; e++) {
          p7_splicegraph_AddNode(graph, exons[e]);
          if(path->node_id[s-1] < graph->orig_N  && path->node_id[s-1] != -1) {
            if(seq_overlap(exons[e]->dcl, path->hits[s-1]->dcl, graph->revcomp) >= 0.0 &&
               hmm_overlap(exons[e]->dcl, path->hits[s-1]->dcl) >= 0.0) {
              graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[path->node_id[s-1]];
              graph->split_orig_id[graph->num_nodes-1] = path->node_id[s-1];
            }
          }
          if(path->node_id[s] < graph->orig_N  && path->node_id[s] != -1) {
            if(seq_overlap(exons[e]->dcl, path->hits[s]->dcl, graph->revcomp) >= 0.0 &&
               hmm_overlap(exons[e]->dcl, path->hits[s]->dcl) >= 0.0) {
              graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[path->node_id[s]];
              graph->split_orig_id[graph->num_nodes-1] = path->node_id[s];
            }
          }
        }
        path->hits[s] = exons[num_hits-1];
  
      }
      else {
        for(e = 0; e < num_hits; e++) {
          p7_trace_splice_Destroy(exons[e]->dcl->tr);
          free(exons[e]->dcl->scores_per_pos);
          p7_hit_Destroy(exons[e]);
        }
        free(exons);
  
        esl_sq_Destroy(sub_seq);
        seq_start = llabs(path->hits[s-1]->dcl->iali - path_seq->start);
        seq_len   = llabs(path->hits[s]->dcl->jali - path->hits[s-1]->dcl->iali) + 1;
        printf("seq_len %d\n", seq_len);
        sub_seq   = esl_sq_CreateDigitalFrom(path_seq->abc, NULL, path_seq->dsq+seq_start, seq_len, NULL,NULL,NULL);
        sub_seq->start = path->hits[s-1]->dcl->iali;
        sub_seq->end   = path->hits[s]->dcl->jali; 
       
        if(remove_idx != NULL) free(remove_idx);
        remove_idx = NULL;
       
        exons = p7_splice_AlignExons(pli, sub_hmm, gm_fs, pli->bg, sub_seq, gcode, graph->revcomp, path->hits[s-1]->dcl->ihmm, 0, remove_idx, &num_hits);      
  
        for(e = 0; e < num_hits; e++) {
        
          if(graph->revcomp) {
            exons[e]->dcl->iali = sub_seq->n - exons[e]->dcl->iali + sub_seq->end;
            exons[e]->dcl->jali = sub_seq->n - exons[e]->dcl->jali + sub_seq->end;
          }
          else {
            exons[e]->dcl->iali = sub_seq->start + exons[e]->dcl->iali - 1;
            exons[e]->dcl->jali = sub_seq->start + exons[e]->dcl->jali - 1;
          }    
      
       
          p7_splicegraph_AddNode(graph, exons[e]);
          if(path->node_id[s-1] < graph->orig_N  && path->node_id[s-1] != -1) {
            if(seq_overlap(exons[e]->dcl, path->hits[s-1]->dcl, graph->revcomp) >= 0.0 &&
               hmm_overlap(exons[e]->dcl, path->hits[s-1]->dcl) >= 0.0) {
              graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[path->node_id[s-1]];
              graph->split_orig_id[graph->num_nodes-1] = path->node_id[s-1];
            }
          }
          if(path->node_id[s] < graph->orig_N  && path->node_id[s] != -1) {
            if(seq_overlap(exons[e]->dcl, path->hits[s]->dcl, graph->revcomp) >= 0.0 &&
               hmm_overlap(exons[e]->dcl, path->hits[s]->dcl) >= 0.0) {
              graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[path->node_id[s]];
              graph->split_orig_id[graph->num_nodes-1] = path->node_id[s];
            }
          }
        }
        path->hits[s] = exons[num_hits-1];
      }     
    }

    free(exons);
    esl_sq_Destroy(sub_seq);
    p7_hmm_Destroy(sub_hmm);
    if(remove_idx != NULL) free(remove_idx);
    if(nuc_index  != NULL) free(nuc_index);
    remove_idx = NULL;
    nuc_index  = NULL;
    
  } 

  return eslOK;

}




SPLICE_PATH2*
p7_splice_FindExons2(SPLICE_WORKER_INFO *info, SPLICE_PATH2 *path, ESL_SQ *path_seq)
{

  int i, s;
  int s_end;
  int intron_len;
  int amino_gap;
  int full_intron;
  int contains_anchor;
  int             *remove_idx;
  int64_t         *nuc_index;
  P7_HMM          *sub_hmm;
  P7_HMM          *hmm;
  P7_FS_PROFILE   *gm_fs;
  ESL_GENCODE     *gcode;
  ESL_SQ          *sub_seq;
  SPLICE_PIPELINE *pli;
  SPLICE_PATH2 *tmp_path;
  SPLICE_PATH2 *ret_path;
  SPLICE_GRAPH *graph;

  graph = info->graph;
  pli   = info->pli;
  gcode = info->gcode;
  gm_fs = info->gm_fs;
  hmm   = info->hmm;

  remove_idx = NULL;
  nuc_index  = NULL;

  /* If there is only one hit in th path we must still relaign it in case that hit aligned through an intron */
  if(path->path_len == 1) {

    sub_hmm     = p7_splice_GetSubHMM(hmm, path->ihmm[0], path->jhmm[0]);
    sub_hmm->fs = 0.;

    tmp_path = p7_splice_AlignExons2(pli, sub_hmm, gm_fs, pli->bg, path_seq, gcode, 0, 0);

    for(s = 0; s < tmp_path->path_len; s++) {
      
      tmp_path->node_id[s] = path->node_id[0];

      tmp_path->ihmm[s] = tmp_path->ihmm[s] + path->ihmm[0] - 1;
      tmp_path->jhmm[s] = tmp_path->jhmm[s] + path->ihmm[0] - 1;
      if(path->revcomp) {
         tmp_path->iali[s] = path_seq->n - tmp_path->iali[s] + path_seq->end;
         tmp_path->jali[s] = path_seq->n - tmp_path->jali[s] + path_seq->end;
      }
      else {
         tmp_path->iali[s] = path_seq->start + tmp_path->iali[s] - 1;
         tmp_path->jali[s] = path_seq->start + tmp_path->jali[s] - 1;
      }
    }

    p7_hmm_Destroy(sub_hmm);
    tmp_path->revcomp = path->revcomp;
    tmp_path->frameshift = path->frameshift; 
    return tmp_path;
  }

  /* Start by splicing from the first to the last anchor node */
  s = 0;
  while (path->node_id[s] >= graph->orig_N) s++;
  
  /* Create the return path and set its most upstream coords. */
  ret_path = p7_splicepath_Create2(1);

  ret_path->node_id[0] = path->node_id[s];
  ret_path->iali[0]    = path->iali[s];
  ret_path->jali[0]    = path->jali[s];
  ret_path->ihmm[0]    = path->ihmm[s];
  ret_path->jhmm[0]    = path->jhmm[s];

  s_end = path->path_len-1;
  while(path->node_id[s_end] >= graph->orig_N) s_end--;

  /* For each neighboring pair of nodes in the path, align the sub hmm to the sub seqqunce (with or without part of the intron removed). */
  s++;
  while(s <= s_end) {
    
printf("RET PATH\n");
p7_splicepath_Dump2(stdout, ret_path);
printf("s-1 %d s %d\n", s, s+1);
printf("ORIG PATH\n");   
//p7_splicepath_Dump2(stdout, path);

    s_end = path->path_len-1;
    while(path->node_id[s_end] >= graph->orig_N) s_end--;   
    

    //TODO If the amino gap is sufficently large start by searching the whole intron
    sub_seq = p7_splice_GetSplicedSequence2(path, path_seq, s-1, s, TRUE, &remove_idx, &nuc_index);

    intron_len = llabs(path->iali[s] - path->jali[s-1]) - 1;
    amino_gap =  ESL_MAX(0, path->ihmm[s] - path->jhmm[s-1] - 1);

    if(intron_len <= MAX_INTRON_INCL + (amino_gap*3)) full_intron = TRUE;
    else                                              full_intron = FALSE;

    /* Crete sub hmm */
    sub_hmm     = p7_splice_GetSubHMM(hmm, path->ihmm[s-1], path->jhmm[s]);
    sub_hmm->fs = 0.;

    /* Align sub seq to sub hmm and get exons */
    tmp_path = p7_splice_AlignExons2(pli, sub_hmm, gm_fs, pli->bg, sub_seq, gcode, remove_idx[0], remove_idx[1]);
   
    /* Adjust coords */ 
    for(i = 0; i < tmp_path->path_len; i++) {
      if(tmp_path->iali[i] == -1) continue; 

      tmp_path->ihmm[i] = tmp_path->ihmm[i] + path->ihmm[s-1] - 1;
      tmp_path->jhmm[i] = tmp_path->jhmm[i] + path->ihmm[s-1] - 1;
      tmp_path->iali[i] = nuc_index[tmp_path->iali[i]];
      tmp_path->jali[i] = nuc_index[tmp_path->jali[i]];
    }
//printf("TMP PATH\n");
//p7_splicepath_Dump2(stdout, tmp_path);    
    assess_paths(info, tmp_path, ret_path, path, sub_hmm, path_seq, full_intron, &s);  
    
    p7_splicepath_Destroy2(tmp_path);
    esl_sq_Destroy(sub_seq);
    p7_hmm_Destroy(sub_hmm);
    if(remove_idx != NULL) free(remove_idx);
    if(nuc_index  != NULL) free(nuc_index);
    remove_idx = NULL;
    nuc_index  = NULL; 
  }       

  ret_path->revcomp = path->revcomp;
  ret_path->frameshift = path->frameshift;

  //TODO double check node assignments

  /* Double check node assisgnments */
  for(s = 0; s < path->path_len; s++) {
    if(path->node_id[s] >= graph->orig_N) continue;
    for(i = 0; i < ret_path->path_len; i++) {
      if(ret_path->node_id[i] < graph->orig_N) continue;
      if(seq_overlap2(ret_path, i, path, s) > 0.0 && hmm_overlap2(ret_path, i, path, s) > 0.0)
        ret_path->node_id[i] = path->node_id[s];
    }
  }
   
  contains_anchor = FALSE;
  for(i = 0; i < ret_path->path_len; i++) {
    if(ret_path->node_id[i] < graph->orig_N)  
      contains_anchor = TRUE;
  }
  if(!contains_anchor) ret_path->path_len = 0;

  return ret_path;

}

/*
 * path1 = tmp_path
 * path2 = ret_path
 * path3 = path
 */
int 
assess_paths(SPLICE_WORKER_INFO *info, SPLICE_PATH2 *path1, SPLICE_PATH2 *path2, SPLICE_PATH2 *path3, P7_HMM *sub_hmm, ESL_SQ *path_seq, int full_intron, int *step)
{
  int i, s; 
  int seq_start;
  int seq_len;
  int upstream_confirmed;
  int downstream_confirmed;
  int split_confirmed;
  ESL_SQ *sub_seq;
  SPLICE_PATH2 *tmp_path;
    
  s = *step;

  upstream_confirmed  = confirm_overlap(path1, 0, path3, s-1, p7_CONFIRM_START);
  downstream_confirmed = confirm_overlap(path1, path1->path_len-1, path3, s, p7_CONFIRM_END); 
  

//printf("upstream_confirmed %d downstream_confirmed %d\n", upstream_confirmed, downstream_confirmed);

  /* Remove any unconfirmed nodes */
  if(!downstream_confirmed) 
    p7_splicepath_Remove2(path3, s); 
  if(!upstream_confirmed)
    remove_upstream(path2, path3, &s);

  if(upstream_confirmed && downstream_confirmed) {
/*    
    split_confirmed = confirm_split(path1, path2);
    printf("split_confirmed %d\n", split_confirmed);
    if(!split_confirmed) {
      if(path2->node_id[path2->path_len-1] >= info->graph->orig_N) {
         remove_upstream(path2, path3, &s);        
      }
      else {
         i = path2->path_len-1;
         while (i >= 0 && path2->node_id[i] == path2->node_id[path2->path_len-1]) i--;
         if(i == 0 && path2->node_id[i] == path2->node_id[path2->path_len-1]) {
            path2->path_len = 1;
            split_confirmed = TRUE;
         }
          else if(ath2->node_id[1] >= info->graph->orig_N
         else if(
          
         }
      }
    }
  }


  if(upstream_confirmed && downstream_confirmed && split_confirmed) {
*/
    if(path1->path_len == 1) {

      path2->jali[path2->path_len-1] = path1->jali[0];
      path2->jhmm[path2->path_len-1] = path1->jhmm[0];

      s++; 
    }
    else if(path1->path_len == 2) {

      path2->jali[path2->path_len-1] = path1->jali[0];
      path2->jhmm[path2->path_len-1] = path1->jhmm[0];       

      p7_splicepath_Grow2(path2);

      path2->node_id[path2->path_len] = path3->node_id[s];
      path2->iali[path2->path_len]    = path1->iali[1]; 
      path2->jali[path2->path_len]    = path1->jali[1];
      path2->ihmm[path2->path_len]    = path1->ihmm[1];
      path2->jhmm[path2->path_len]    = path1->jhmm[1];
      path2->path_len++;

      s++;
      
    }   
    else if(full_intron) {
//printf("TMP PATH\n");
//p7_splicepath_Dump2(stdout, path1);

      path2->jali[path2->path_len-1] = path1->jali[0];
      path2->jhmm[path2->path_len-1] = path1->jhmm[0];  

      for(i = 1; i < path1->path_len; i++) {
        p7_splicepath_Grow2(path2);
        path2->node_id[path2->path_len] = -1;  

        if(seq_overlap2(path1, i, path3, s-1) > 0.0 && hmm_overlap2(path1, i, path3, s-1) > 0.0)  
          path2->node_id[path2->path_len] = path3->node_id[s-1];
        else if(seq_overlap2(path1, i, path3, s) > 0.0 && hmm_overlap2(path1, i, path3, s) > 0.0)
          path2->node_id[path2->path_len] = path3->node_id[s];

        path2->iali[path2->path_len] = path1->iali[i];
        path2->jali[path2->path_len] = path1->jali[i];
        path2->ihmm[path2->path_len] = path1->ihmm[i];
        path2->jhmm[path2->path_len] = path1->jhmm[i];
        path2->path_len++;     
      }

      s++;
    } 
    else {
      //TODO Figure a way to cut very long introns up in to smalled parts 

      seq_start = llabs(path3->iali[s-1] - path_seq->start);
      seq_len   = llabs(path3->jali[s]   - path3->iali[s-1]) + 1;

      sub_seq   = esl_sq_CreateDigitalFrom(path_seq->abc, NULL, path_seq->dsq+seq_start, seq_len, NULL,NULL,NULL);
      sub_seq->start = path3->iali[s-1];
      sub_seq->end   = path3->jali[s];

      tmp_path = p7_splice_AlignExons2(info->pli, sub_hmm, info->gm_fs, info->pli->bg, sub_seq, info->gcode, 0, 0);
      
      for(i = 0; i < tmp_path->path_len; i++) {

        tmp_path->ihmm[i] = tmp_path->ihmm[i] + path3->ihmm[s-1] - 1;
        tmp_path->jhmm[i] = tmp_path->jhmm[i] + path3->ihmm[s-1] - 1;
        if(path3->revcomp) {
          tmp_path->iali[i] = sub_seq->n - tmp_path->iali[i] + sub_seq->end;
          tmp_path->jali[i] = sub_seq->n - tmp_path->jali[i] + sub_seq->end;
        }
        else {
          tmp_path->iali[i] = sub_seq->start + tmp_path->iali[i] - 1;
          tmp_path->jali[i] = sub_seq->start + tmp_path->jali[i] - 1;
        }
      }       
      assess_paths(info, tmp_path, path2, path3, sub_hmm, path_seq, TRUE, &s);      
      p7_splicepath_Destroy2(tmp_path);
      esl_sq_Destroy(sub_seq);
    }
  }
 
  *step = s;

  return eslOK;
}
P7_HIT**
p7_splice_AlignExons(SPLICE_PIPELINE *pli, P7_HMM *sub_hmm, const P7_FS_PROFILE *gm_fs, P7_BG *bg, ESL_SQ *ali_seq, const ESL_GENCODE *gcode, int revcomp, int hmm_start, int num_introns, int *removed_idx, int *num_exons)
{
  int         i, r, y, z;
  int         z1, z2;
  int         intron_cnt;
  int         exon_cnt;
  int         start_new;
  P7_HIT       *new_hit;
  P7_HIT      **ret_hits;
  P7_TRACE     *tr;
  P7_FS_PROFILE *sub_fs_model;
  int status;

  
  sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
  p7_ProfileConfig_fs(sub_hmm, bg, gcode, sub_fs_model, ali_seq->n, p7_UNIGLOBAL);
   
  p7_gmx_fs_GrowTo(pli->vit, sub_fs_model->M, ali_seq->n, ali_seq->n, p7P_SPLICE);
  tr = p7_trace_fs_Create();

  if(sub_hmm->fs > 0.0) {
    p7_splicepipline_GrowIndex(pli->sig_idx, sub_fs_model->M, ali_seq->n);
    p7_spliceviterbi_translated_semiglobal(pli, ali_seq->dsq, gcode, ali_seq->n, sub_fs_model, pli->vit);
    p7_splicevitebi_translated_semiglobal_trace(pli, ali_seq->dsq, ali_seq->n, gcode, sub_fs_model, pli->vit, tr);   
  }
  else {
    p7_splicepipline_GrowIndex(pli->sig_idx, sub_fs_model->M, ali_seq->n);
    p7_spliceviterbi_translated_semiglobal(pli, ali_seq->dsq, gcode, ali_seq->n, sub_fs_model, pli->vit);
    p7_splicevitebi_translated_semiglobal_trace(pli, ali_seq->dsq, ali_seq->n, gcode, sub_fs_model, pli->vit, tr);
  }
  //p7_trace_fs_Dump(stdout, tr, NULL, NULL, NULL);
  
  /* Find number of introns in trace */
  intron_cnt = 0;
  for(z = 0; z < tr->N; z++)
    if(tr->st[z] == p7T_P) intron_cnt++;

  ret_hits = NULL;
  ESL_ALLOC(ret_hits, sizeof(P7_HIT*) * (intron_cnt+1));

  /* Find first M state - start of first hit */
  for(z1 = 0; z1 < tr->N; z1++) if(tr->st[z1] == p7T_M) break;

  /* Find last M state state - end of last hit */
  for(z2 = tr->N-1; z1 >= 0; z2--) if(tr->st[z2] == p7T_M) break;

  exon_cnt = 0;
  start_new = TRUE;
  
  z = z1;
  while(z < z2) {

    if(start_new) {

      /* Save z value - currently set to fist M state in exon */
      y = z;

      /*Find end of exon */
      while(tr->st[z] != p7T_P && tr->st[z] != p7T_E) z++;
      z--;
     /*Trace back to last M state of exon*/
      //while(tr->st[z] != p7T_M) z--;


      /* If an exon crosses the bounrdy of a removed intron it is 
       * almost certainly a false positive and can be discarded */
      for(r = 0; r < num_introns; r++) {
       
        if(removed_idx[r*2] >= 0 ) {
          if      (tr->i[y] - tr->c[y] + 1 <= removed_idx[r*2]   && tr->i[z] >= removed_idx[r*2])   { start_new = FALSE; break; }
          else if (tr->i[y] - tr->c[y] + 1 <= removed_idx[r*2+1] && tr->i[z] >= removed_idx[r*2+1]) { start_new = FALSE; break; }
          else if (tr->i[y] - tr->c[y] + 1 >  removed_idx[r*2]   && tr->i[z] <  removed_idx[r*2+1]) { start_new = FALSE; break; }
        }
      }
      if( !start_new) continue;

      /* Create new hit and  set ihmm and iali coords*/
      new_hit          = p7_hit_Create_empty();
      new_hit->dcl     = p7_domain_Create_empty();
      new_hit->dcl->tr = p7_trace_fs_Create();

      new_hit->dcl->ihmm = tr->k[y] + hmm_start - 1;
      new_hit->dcl->jhmm = tr->k[z] + hmm_start - 1;
      new_hit->dcl->iali = tr->i[y] - tr->c[y] + 1;
      new_hit->dcl->jali = tr->i[z];

      
      /* Append starting special states */
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_S , 0, 0, 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_N , 0, tr->i[y]-3, 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_B , 0, tr->i[y]-3, 0);

      /* Append all states between first and last M state */
      for(i = y; i <= z; i++)  
        p7_trace_fs_Append(new_hit->dcl->tr, tr->st[i], tr->k[i]+hmm_start-1, tr->i[i], tr->c[i]);

      /* Append ending special states */
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_E, tr->k[z]+hmm_start-1, tr->i[z], 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_C, 0, tr->i[z], 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_T, 0, 0, 0);

      new_hit->dcl->ad             = NULL;
      new_hit->dcl->scores_per_pos = NULL;
      
      p7_splice_ComputeAliScores_fs(new_hit->dcl, new_hit->dcl->tr, ali_seq, gm_fs, bg, TRUE);

      ret_hits[exon_cnt] = new_hit;
      exon_cnt++;

      start_new = FALSE;
    }

    z++;
    if(tr->st[z] == p7T_M) start_new = TRUE;

  }

  *num_exons = exon_cnt;

  p7_profile_fs_Destroy(sub_fs_model);
  p7_trace_fs_Destroy(tr);

  return ret_hits;


  ERROR:
   if(ret_hits != NULL) free(ret_hits);
   return NULL;

}


SPLICE_PATH2*
p7_splice_AlignExons2(SPLICE_PIPELINE *pli, P7_HMM *sub_hmm, const P7_FS_PROFILE *gm_fs, P7_BG *bg, ESL_SQ *ali_seq, const ESL_GENCODE *gcode, int removed_start, int removed_end)
{
  int         y, z;
  int         z1, z2;
  int         intron_cnt;
  int         step_cnt;
  int         start_new;
  P7_TRACE     *tr;
  P7_FS_PROFILE *sub_fs_model;
  SPLICE_PATH2   *ret_path;


  sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
  p7_ProfileConfig_fs(sub_hmm, bg, gcode, sub_fs_model, ali_seq->n, p7_UNIGLOBAL);

  p7_gmx_fs_GrowTo(pli->vit, sub_fs_model->M, ali_seq->n, ali_seq->n, p7P_SPLICE);
  tr = p7_trace_fs_Create();

  if(sub_hmm->fs > 0.0) {
    p7_splicepipline_GrowIndex(pli->sig_idx, sub_fs_model->M, ali_seq->n);
    p7_spliceviterbi_translated_semiglobal(pli, ali_seq->dsq, gcode, ali_seq->n, sub_fs_model, pli->vit);
    p7_splicevitebi_translated_semiglobal_trace(pli, ali_seq->dsq, ali_seq->n, gcode, sub_fs_model, pli->vit, tr);
  }
  else {
    p7_splicepipline_GrowIndex(pli->sig_idx, sub_fs_model->M, ali_seq->n);
    p7_spliceviterbi_translated_semiglobal(pli, ali_seq->dsq, gcode, ali_seq->n, sub_fs_model, pli->vit);
    p7_splicevitebi_translated_semiglobal_trace(pli, ali_seq->dsq, ali_seq->n, gcode, sub_fs_model, pli->vit, tr);
  }
  //p7_trace_fs_Dump(stdout, tr, NULL, NULL, NULL);

  /* Find number of introns in trace */
  intron_cnt = 0;
  for(z = 0; z < tr->N; z++)
    if(tr->st[z] == p7T_P) intron_cnt++;

  ret_path = p7_splicepath_Create2(intron_cnt+1);

  /* Find first M state - start of first hit */
  for(z1 = 0; z1 < tr->N; z1++) if(tr->st[z1] == p7T_M) break;

  /* Find last M state state - end of last hit */
  for(z2 = tr->N-1; z1 >= 0; z2--) if(tr->st[z2] == p7T_M) break;

  step_cnt = 0;
  start_new = TRUE;

  z = z1;
  while(z <= z2) {

    if(start_new) {
      
      /* Save z value - currently set to fist M state in exon */
      y = z;

      /*Find end of exon */
      while(tr->st[z] != p7T_P && tr->st[z] != p7T_E) z++;
      if(tr->st[z] == p7T_E) while(tr->st[z] != p7T_M) z--;
      else z--;
  
      /* If an exon crosses the bounrdy of a removed intron we set all indicies to -1 */
      if(removed_start > 0 ) {
        if( (tr->i[y] - tr->c[y] + 1 <= removed_start && tr->i[z] >= removed_start)   ||
            (tr->i[y] - tr->c[y] + 1 <= removed_end   && tr->i[z] >= removed_end) ||
            (tr->i[y] - tr->c[y] + 1 >  removed_start && tr->i[z] <  removed_end)) {

          ret_path->node_id[step_cnt] = -1;
          ret_path->iali[step_cnt] = -1; 
          ret_path->jali[step_cnt] = -1;
          ret_path->ihmm[step_cnt] = -1;
          ret_path->jhmm[step_cnt] = -1;
          step_cnt++;
          start_new = FALSE;
          continue;
        }
      }

      /* If this is the first step in path the i coords are set my the first M state at
       * trace postion y, otherwise they are set by the P state at trace postion y-1 */
      ret_path->node_id[step_cnt] = -1;
      if(step_cnt == 0) {
        ret_path->iali[step_cnt] = tr->i[y] - tr->c[y] + 1; 
        ret_path->ihmm[step_cnt] = tr->k[y];        
      }
      else {
        /* Determine which splice codon type we have */
        ret_path->ihmm[step_cnt] = tr->k[y-1];
        if( tr->c[y-1] == 0)
          ret_path->iali[step_cnt] = tr->i[y-1] - 2; 
        else if( tr->c[y-1] == 1)
          ret_path->iali[step_cnt] = tr->i[y-1] - 1;
        else
          ret_path->iali[step_cnt] = tr->i[y-1];
      }

      ret_path->jhmm[step_cnt] = tr->k[z];

      if(step_cnt == ret_path->path_len-1) {
        ret_path->jali[step_cnt]   = tr->i[z];
      }
      else {
        /* Determine which splice codon type we have */
        if( tr->c[z+1] == 0)
          ret_path->jali[step_cnt]   = tr->i[z];
        else if( tr->c[z+1] == 1)
          ret_path->jali[step_cnt]   = tr->i[z]+1;
        else if( tr->c[z+1] == 2)
          ret_path->jali[step_cnt]   = tr->i[z]+2;
      }

      step_cnt++;

      start_new = FALSE;
    }
    
    z++;
    if(tr->st[z] == p7T_M) start_new = TRUE;

  }

  p7_profile_fs_Destroy(sub_fs_model);
  p7_trace_fs_Destroy(tr);

  return ret_path;

}




float
hmm_overlap(P7_DOMAIN *dom_1, P7_DOMAIN *dom_2) 
{

  int dom_1_hmm_len;
  int overlap_hmm_start; 
  int overlap_hmm_end;
  int overlap_hmm_len;


  /* Find the hmm coords overlap betweeen dom_1 and dom_2. */
  overlap_hmm_start = ESL_MAX(dom_1->ihmm, dom_2->ihmm);
  overlap_hmm_end   = ESL_MIN(dom_1->jhmm, dom_2->jhmm);

  overlap_hmm_len   = overlap_hmm_end - overlap_hmm_start + 1;

  /* If the length of the overlap is more than <pct_overlap> of the length of dom_1 return TRUE, otherwise return FALSE */
  dom_1_hmm_len = dom_1->jhmm - dom_1->ihmm + 1;

  return (float) overlap_hmm_len / (float) dom_1_hmm_len; 

}

float
hmm_overlap2(SPLICE_PATH2 *path1, int step1, SPLICE_PATH2 *path2, int step2)
{

  int hmm_len_1;
  int overlap_hmm_start;
  int overlap_hmm_end;
  int overlap_hmm_len;

  overlap_hmm_start = ESL_MAX(path1->ihmm[step1], path2->ihmm[step2]);
  overlap_hmm_end   = ESL_MIN(path1->jhmm[step1], path2->jhmm[step2]);
  overlap_hmm_len   = overlap_hmm_end - overlap_hmm_start + 1;
  hmm_len_1         = path1->jhmm[step1] - path1->ihmm[step1] + 1;

  return (float) overlap_hmm_len / (float) hmm_len_1;

}


float 
seq_overlap(P7_DOMAIN *dom_1, P7_DOMAIN *dom_2, int revcomp) 
{ 

  int dom_1_seq_len;
  int overlap_seq_start;
  int overlap_seq_end;
  int overlap_seq_len;

  /* Find the seq coords overlap betweeen dom_1 and dom_2. */
  if(revcomp) {
    overlap_seq_start = ESL_MAX(dom_1->jali, dom_2->jali);
    overlap_seq_end   = ESL_MIN(dom_1->iali, dom_2->iali);
  } 
  else {
    overlap_seq_start = ESL_MAX(dom_1->iali, dom_2->iali);
    overlap_seq_end   = ESL_MIN(dom_1->jali, dom_2->jali);
  }
  overlap_seq_len = overlap_seq_end - overlap_seq_start + 1;
  
  /* If the length of the overlap is more than <pct_overlap> of the length of dom_1 return TRUE, otherwise return FALSE */
  dom_1_seq_len = llabs(dom_1->jali - dom_1->iali) + 1;
  return (float) overlap_seq_len / (float) dom_1_seq_len;
  
}

float
seq_overlap2(SPLICE_PATH2 *path1, int step1, SPLICE_PATH2 *path2, int step2) 
{

  int seq_len_1;
  int overlap_seq_start;
  int overlap_seq_end;
  int overlap_seq_len;

   if(path2->revcomp) {
    overlap_seq_start = ESL_MAX(path1->jali[step1], path2->jali[step2]);
    overlap_seq_end   = ESL_MIN(path1->iali[step1], path2->iali[step2]);
  }
  else {
    overlap_seq_start = ESL_MAX(path1->iali[step1], path2->iali[step2]);
    overlap_seq_end   = ESL_MIN(path1->jali[step1], path2->jali[step2]);
  }
  overlap_seq_len = overlap_seq_end - overlap_seq_start + 1;
  seq_len_1 = llabs(path1->jali[step1] - path1->iali[step1]) + 1;

  return (float) overlap_seq_len / (float) seq_len_1;

}

/*
 * path1 = tmp_path
 * path2 = path
 */
int
confirm_overlap(SPLICE_PATH2 *path1, int step1, SPLICE_PATH2 *path2, int step2, int confirm_side)
{    
  int hmm_dist;
  int seq_dist;
    
  if(path1->iali[step1] == -1) return FALSE; //path 1 is invalid

  if(seq_overlap2(path1, step1, path2, step2) <= 0.0) return FALSE; // no sequence overlap 
  if(hmm_overlap2(path1, step1, path2, step2) <= 0.0) return FALSE; // no model overlap 
 
  if(confirm_side == p7_CONFIRM_START) {
    if(path1->iali[step1] % 3 != path2->iali[step2] % 3) return FALSE; // start on differnent frames

    hmm_dist = abs(path1->ihmm[step1] - path2->ihmm[step2]);
    seq_dist = llabs(path1->iali[step1] - path2->iali[step2]);
    if(hmm_dist > 3 ||  seq_dist > 9) return FALSE; // starts do not match
    
  }
  else if(confirm_side == p7_CONFIRM_END) {
    if(path1->jali[step1] % 3 != path2->jali[step2] % 3); // end on differnent frames
  
    hmm_dist = abs(path1->jhmm[step1] - path2->jhmm[step2]);    
    seq_dist = llabs(path1->jali[step1] - path2->jali[step2]); 
    if(hmm_dist > 3 ||  seq_dist > 9) return FALSE; // ends do not match
    
  }
 
  return TRUE; 
}

/*
 * path1 = tmp_path
 * path2 = ret_path
 */
int
confirm_split(SPLICE_PATH2 *path1, SPLICE_PATH2 *path2)
{
  int i, s;
  int split_cnt;

  /*If the last hit in the added to path2 was a split hit, 
   * does the first hit in path 1 agree with that split */
  if(path2->path_len < 2) return TRUE;
  if(path2->node_id[path2->path_len-1] == -1) return TRUE;
  if(path2->node_id[path2->path_len-1] != path2->node_id[path2->path_len-2]) return TRUE;
  
  s = path2->path_len-1;
  while(path2->node_id[s] == path2->node_id[path2->path_len-1]) s--;
  s++;
 
  split_cnt = path2->path_len - s;
  printf("split_cnt %d\n", split_cnt);
  if(path1->path_len < split_cnt) return FALSE; // not enough hits to cover split
  if(!confirm_overlap(path1, 0, path2, s, p7_CONFIRM_END)) return FALSE;
  s++;
  i = 1;
  while(s < path2->path_len-1) {
    if(!confirm_overlap(path1, i, path2, s, p7_CONFIRM_START)) return FALSE;
    if(!confirm_overlap(path1, i, path2, s, p7_CONFIRM_END)) return FALSE;
    i++;
    s++;
  } 
  if(!confirm_overlap(path1, i, path2, s, p7_CONFIRM_START)) return FALSE;

  return TRUE;

}


/*
 * path1 = ret_path
 * path2 = path
 */
int
remove_upstream(SPLICE_PATH2 *path1, SPLICE_PATH2 *path2, int *step2) 
{

  int i, s;
  int removed_node;
  int split_node;
  s = *step2;
  
  removed_node = path2->node_id[s-1];
  p7_splicepath_Remove2(path2, s-1);
//printf("removed s-1 %d\n", s);
  /* We removed the fist node in path2 and must restart path1 */ 
  if(s == 1 || path2->path_len == 1) {
    path1->node_id[0] = path2->node_id[0];
    path1->iali[0]    = path2->iali[0];
    path1->jali[0]    = path2->jali[0];
    path1->ihmm[0]    = path2->ihmm[0];
    path1->jhmm[0]    = path2->jhmm[0];
    path1->path_len   = 1;
  }
  /* Move path end to the first hit from the previous upstream node */
  else {
    /* Move end of path1 to before removed node */
//printf("removed_node %d\n", removed_node+1);
    while(path1->node_id[path1->path_len-1] == removed_node) path1->path_len--;
//p7_splicepath_Dump2(stdout, path1);
    /* Move end of path1 to befoer any -1 nodes */
    while(path1->node_id[path1->path_len-1] == -1) path1->path_len--;
    /* Move end of path1 to begining of split nodes */
    split_node = path1->node_id[path1->path_len-1];
    while(path1->node_id[path1->path_len-1] == split_node) {
      path1->path_len--;
      if(path1->path_len == 0) break;
    }
    path1->path_len++;

    /*Move s to after node in path 2 with same node id as last node in path 1*/
    for(i = 0; i < path2->path_len; i++) 
      if(path2->node_id[i] == path1->node_id[path1->path_len-1]) break;
    s = i+1;
  }
  
  *step2 = s;

  return eslOK;
}


int
p7_splice_ConnectGraph(SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, const P7_FS_PROFILE *gm_fs, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQFILE *seq_file, SPLICE_WORKER_INFO *info, ESL_SQ **path_seq_accumulator, int num_paths) 
{
  
  int p;
  int up, down;
  int up_min, up_max;
  int down_min, down_max;
  int seq_min, seq_max;
  int seq_len;
  P7_TOPHITS  *th;
  ESL_DSQ     *splice_dsq;
  ESL_SQ      *splice_seq;
  SPLICE_EDGE *edge;
  SPLICE_EDGE *tmp_edge;

  th = graph->th;

  for(up = 0; up < graph->num_nodes; up++) {
    if(!graph->node_in_graph[up]) continue;

    for(down = 0; down < graph->num_nodes; down++) {
      if(!graph->node_in_graph[down]) continue;

      /*Nodes that are split from the same original hit will already have edges */
      //if(graph->split_orig_id[up] != -1 && graph->split_orig_id[up] == graph->split_orig_id[down]) continue;
     
      if(!p7_splice_HitUpstream(th->hit[up]->dcl, th->hit[down]->dcl, graph->revcomp)) continue;
   
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
      seq_len  = seq_max - seq_min + 1;

      splice_seq = NULL;
      /* Check the existing path_seqs to find thw splice_seq. If this edge 
       * is not in an existing path_seq, fetch a neew subseqeucne */
      for(p = 0; p < num_paths; p++) {
 
        if(graph->revcomp) {
          if(seq_min >= path_seq_accumulator[p]->end && seq_max <= path_seq_accumulator[p]->start) {
            splice_dsq = path_seq_accumulator[p]->dsq+(path_seq_accumulator[p]->start-seq_max);
            splice_seq = esl_sq_CreateDigitalFrom(path_seq_accumulator[p]->abc, NULL, splice_dsq, seq_len, NULL,NULL,NULL);

            splice_seq->start = path_seq_accumulator[p]->start - (path_seq_accumulator[p]->start-seq_max);
            splice_seq->end   = path_seq_accumulator[p]->end + (seq_min-path_seq_accumulator[p]->end);

            p = num_paths;
          }

        }
        else { 
          if(seq_min >= path_seq_accumulator[p]->start && seq_max <= path_seq_accumulator[p]->end) {
            splice_dsq = path_seq_accumulator[p]->dsq+(seq_min-path_seq_accumulator[p]->start);
            splice_seq = esl_sq_CreateDigitalFrom(path_seq_accumulator[p]->abc, NULL, splice_dsq, seq_len, NULL,NULL,NULL);        

            splice_seq->start = path_seq_accumulator[p]->start + (seq_min-path_seq_accumulator[p]->start);
            splice_seq->end   = path_seq_accumulator[p]->end - (path_seq_accumulator[p]->end-seq_max);

            p = num_paths;
          }
        }
      }

      if(splice_seq == NULL)
        splice_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, graph->revcomp, info);
  
     // printf("up %d down %d \n", up+1, down+1);
      edge = p7_spliceedge_ConnectNodes(pli, th->hit[up]->dcl, th->hit[down]->dcl, gm_fs, hmm, gcode, splice_seq, graph->revcomp);
       
      if(edge != NULL) {
        tmp_edge = p7_splicegraph_AddEdge(graph, up, down);
        tmp_edge->frameshift = edge->frameshift;
        tmp_edge->bypass_checked = edge->bypass_checked;
        tmp_edge->upstream_nuc_start = edge->upstream_nuc_start;
        tmp_edge->upstream_nuc_end = edge->upstream_nuc_end;
        tmp_edge->downstream_nuc_start = edge->downstream_nuc_start;
        tmp_edge->downstream_nuc_end = edge->downstream_nuc_end;
        tmp_edge->upstream_spliced_amino_end = edge->upstream_spliced_amino_end;
        tmp_edge->downstream_spliced_amino_start = edge->downstream_spliced_amino_start;
        tmp_edge->upstream_spliced_nuc_end = edge->upstream_spliced_nuc_end;
        tmp_edge->downstream_spliced_nuc_start = edge->downstream_spliced_nuc_start;
        tmp_edge->splice_score = edge->splice_score;
        tmp_edge->signal_score = edge->signal_score; 
        free(edge);  
      }
      esl_sq_Destroy(splice_seq);   
    }

  }

  return eslOK;

}


int
p7_splice_HitUpstream(P7_DOMAIN *upstream, P7_DOMAIN *downstream, int revcomp) {

  if(upstream->ihmm > downstream->ihmm || upstream->jhmm > downstream->jhmm)
    return FALSE;
    
 
  if (( revcomp  && upstream->jali <= downstream->iali) ||
     ((!revcomp) && upstream->jali >= downstream->iali))
    return FALSE;

  return TRUE;
  
}


int
hits_spliceable(P7_DOMAIN *upstream, P7_DOMAIN *downstream) {
  

  /* Are nodes close enough on hmm coords */ 
  if(upstream->jhmm + MAX_SP_AMINO_GAP < downstream->ihmm) return FALSE;

  /* Are nodes close enough on seq coords */
  if(llabs(downstream->iali - upstream->jali) > MAX_INTRON_LENG) return FALSE;
 
  if(hmm_overlap(upstream, downstream) > 0.5) return FALSE;
  if(hmm_overlap(downstream, upstream) > 0.5) return FALSE;
  
  return TRUE;
  
}










int
p7_splice_AlignPath(SPLICE_GRAPH *graph, SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_TOPHITS *tophits, P7_OPROFILE *om, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *path_seq, int64_t db_nuc_cnt, float fs_prob, SPLICE_WORKER_INFO *info)
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
  int          path_start_pos, path_end_pos;
  int          ext_start_pos, ext_end_pos;
  float        dom_bias;
  float        nullsc;
  float        dom_score;
  double       dom_lnP;
  int64_t     *nuc_index;
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
  for (i = 0; i < path->path_len; i++) {
    path_seq_len += abs(path->upstream_spliced_nuc_end[i+1] - path->downstream_spliced_nuc_start[i]) + 1;
    //printf("i %d path->end[i+1] %d path->start[i] %d abs(end-start) +1 %d\n", i+1, path->upstream_spliced_nuc_end[i+1], path->downstream_spliced_nuc_start[i], abs(path->upstream_spliced_nuc_end[i+1] - path->downstream_spliced_nuc_start[i]) + 1);
  }
//  printf("1 path_seq_len %d\n", path_seq_len);
  /* If the spliced seqeunce length is non-mod 3 send to farmshift alignment */
  if(path_seq_len % 3 != 0) {
    path->frameshift = TRUE;
    return eslOK;
  }
// printf("2 path_seq_len %d\n", path_seq_len);
  
  /* Working backward from the start of the path find the first instance of a stop codon in the extended upstream region of path_seq */
  ext_start_pos = 1;
  if (path->revcomp) {
    path_start_pos = path_seq->n - path->downstream_spliced_nuc_start[0] + path_seq->end;
    for(pos = path->downstream_spliced_nuc_start[0]+3; pos <= path->downstream_spliced_nuc_start[0] + ALIGNMENT_EXT*3; pos+=3) {
      seq_pos = path_seq->n - pos + path_seq->end;   
      if(seq_pos < 1) {
        ext_start_pos = seq_pos+3;
        break;
      }
      amino = esl_gencode_GetTranslation(gcode,&path_seq->dsq[seq_pos]);
      
      if(esl_abc_XIsNonresidue(gm->abc, amino)) { 
        ext_start_pos = seq_pos+3;
        break;
      }
    }
  }
  else {
    path_start_pos = path->downstream_spliced_nuc_start[0] - path_seq->start + 1;
    for(pos = path->downstream_spliced_nuc_start[0]-3; pos >= path->downstream_spliced_nuc_start[0] - ALIGNMENT_EXT*3; pos-=3) {
      
      seq_pos = pos - path_seq->start + 1;
      
      if(seq_pos < 1) {
        ext_start_pos = seq_pos+3;
        break;
      }
      amino = esl_gencode_GetTranslation(gcode,&path_seq->dsq[seq_pos]);
 
      if(esl_abc_XIsNonresidue(gm->abc, amino)) {
        ext_start_pos = seq_pos+3;
        break;
      }
    }
    
  }


  /* Working forward from the end of the path find the first instance of a stop codon in the extended downstream region of path_seq */ 
  ext_end_pos = path_seq->n;
  if (path->revcomp) {
    path_end_pos =  path_seq->n - path->upstream_spliced_nuc_end[path->path_len] + path_seq->end;
    
    for(pos = path->upstream_spliced_nuc_end[path->path_len]-1; pos >= path->upstream_spliced_nuc_end[path->path_len] - ALIGNMENT_EXT*3; pos-=3) {
      seq_pos = path_seq->n - pos + path_seq->end;
      if(seq_pos > path_seq->n-2) {
        ext_end_pos = seq_pos-1;
        break;
      }
      
      amino = esl_gencode_GetTranslation(gcode,&path_seq->dsq[seq_pos]);
 
      if(esl_abc_XIsNonresidue(gm->abc, amino)) {
        ext_end_pos = seq_pos-1;
        break;
      } 
 
    }
  }
  else {
    path_end_pos = path->upstream_spliced_nuc_end[path->path_len] - path_seq->start + 1;

    for(pos = path->upstream_spliced_nuc_end[path->path_len]+1; pos <= path->upstream_spliced_nuc_end[path->path_len] +  ALIGNMENT_EXT*3; pos+=3) {

      seq_pos = pos - path_seq->start + 1;
      if(seq_pos > path_seq->n-2 ) {
        ext_end_pos = seq_pos-1;
        break;
      }

      amino = esl_gencode_GetTranslation(gcode,&path_seq->dsq[seq_pos]);

      if(esl_abc_XIsNonresidue(gm->abc, amino)) {
        ext_end_pos = seq_pos-1;
        break;
      }
    }
  }

  path_seq_len += (path_start_pos - ext_start_pos) + (ext_end_pos - path_end_pos);
 
  ESL_ALLOC(nuc_index, sizeof(int64_t) * (path_seq_len+2));
  ESL_ALLOC(nuc_dsq,   sizeof(ESL_DSQ) * (path_seq_len+2));

  nuc_index[0] = -1;
  nuc_dsq[0]   = eslDSQ_SENTINEL;
  seq_idx   = 1;

  /*Add upstream extension nucleotides */
  for(seq_pos = ext_start_pos; seq_pos < path_start_pos; seq_pos++) {
    nuc_index[seq_idx] = seq_pos;
    nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
    seq_idx++;
  }

  /* Copy spliced nucleotides into single sequence and track their original indicies */
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

  /*Add downstream extension nucleotides */
  for(seq_pos = path_end_pos+1; seq_pos <= ext_end_pos; seq_pos++) {
    nuc_index[seq_idx] = seq_pos;
    nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
    seq_idx++;
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
    //printf("seq_idx %d\n", seq_idx);  
    if(esl_abc_XIsNonresidue(gm->abc, amino)) {
      path->frameshift = TRUE;
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
  amino_seq   = esl_sq_CreateDigitalFrom(gcode->aa_abc, path_seq->name, amino_dsq, amino_len,    NULL,NULL,NULL);
  nuc_seq     = esl_sq_CreateDigitalFrom(gcode->nt_abc, path_seq->name, nuc_dsq,   path_seq_len, NULL,NULL,NULL);

  pli->nuc_sq       = nuc_seq;
  pli->amino_sq     = amino_seq;
  pli->orig_nuc_idx = nuc_index;

  //p7_omx_SetDumpMode(stdout, pli->bwd, TRUE);
  /* Algin the splices Amino sequence */
  status = align_spliced_path(pli, om, gm, path_seq, gcode, fs_prob);

   /* Alignment failed */
  if(pli->hit == NULL || pli->hit->dcl->ad->exon_cnt == 1) {

    if(nuc_dsq   != NULL) free(nuc_dsq);
    if(amino_dsq != NULL) free(amino_dsq);
     
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
        if(graph->split_orig_id[path->node_id[exon]] >= 0) contains_orig = TRUE;
      }

      if(!contains_orig) {
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
    while( graph->split_orig_id[path->node_id[i]] < 0 ) {
      pli->hit->dcl->ad->exon_orig[i] = FALSE;
      pli->hit->dcl->ad->exon_split[i] = path->split[i];
      i++;
    }
    pli->hit->dcl->ad->exon_orig[i] = TRUE;
    pli->hit->dcl->ad->exon_split[i] = path->split[i];

#ifdef HMMER_THREADS
  if(info->thread_id >= 0) pthread_mutex_lock(info->mutex);
#endif /*HMMER_THREADS*/

    replace_node = graph->split_orig_id[path->node_id[i]];
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
      remove_node = graph->split_orig_id[path->node_id[i]];
      if(remove_node < 0) {
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

#ifdef HMMER_THREADS
  if(info->thread_id >= 0) pthread_mutex_unlock(info->mutex);
#endif /*HMMER_THREADS*/

  }

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
p7_splice_AlignPath2(SPLICE_GRAPH *graph, SPLICE_PATH2 *path, SPLICE_PIPELINE *pli, P7_TOPHITS *tophits, P7_OPROFILE *om, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *path_seq, int64_t db_nuc_cnt, float fs_prob, SPLICE_WORKER_INFO *info)
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
  int          path_start_pos, path_end_pos;
  int          ext_start_pos, ext_end_pos;
  float        dom_bias;
  float        nullsc;
  float        dom_score;
  double       dom_lnP;
  int64_t     *nuc_index;
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
  for (i = 0; i < path->path_len; i++) {
    path_seq_len += abs(path->jali[i] - path->iali[i]) + 1;
  }
  
  /* If the spliced seqeunce length is non-mod 3 send to farmshift alignment */
  if(path_seq_len % 3 != 0) {
    path->frameshift = TRUE;
    ESL_XEXCEPTION(eslFAIL, "NOT mod 3");
    return eslOK;
  }
//printf("start %d end %d n %d\n", path_seq->start, path_seq->end, path_seq->n);
  /* Working backward from the start of the path find the first instance of a stop codon in the extended upstream region of path_seq */
  //printf("path->revcomp %d\n", path->revcomp);
  if (path->revcomp) {
    path_start_pos = path_seq->n - path->iali[0] + path_seq->end;
    ext_start_pos  = path_seq->n - (path->iali[0] + ALIGNMENT_EXT*3) + path_seq->end;
    for(pos = path->iali[0]+3; pos <= path->iali[0] + ALIGNMENT_EXT*3; pos+=3) {

      seq_pos = path_seq->n - pos + path_seq->end;   
      if(seq_pos < 1) {
        ext_start_pos = seq_pos+3;
        break;
      }
       
      amino = esl_gencode_GetTranslation(gcode,&path_seq->dsq[seq_pos]);
      
      if(esl_abc_XIsNonresidue(gm->abc, amino)) { 
        ext_start_pos = seq_pos+3;
        break;
      }
    }
  }
  else {
    path_start_pos = path->iali[0] - path_seq->start + 1;
    ext_start_pos  = (path->iali[0] - ALIGNMENT_EXT*3) - path_seq->start + 1;
    //printf("ext_start_pos %d\n", ext_start_pos);
    for(pos = path->iali[0]-3; pos >= path->iali[0] - ALIGNMENT_EXT*3; pos-=3) {
       
      seq_pos = pos - path_seq->start + 1;
      
      if(seq_pos < 1) {
        ext_start_pos = seq_pos+3;
        break;
      }
      amino = esl_gencode_GetTranslation(gcode,&path_seq->dsq[seq_pos]);
 
      if(esl_abc_XIsNonresidue(gm->abc, amino)) {
        ext_start_pos = seq_pos+3;
        break;
      }
    }
   //printf("ext_start_pos %d\n", ext_start_pos); 
  }


  /* Working forward from the end of the path find the first instance of a stop codon in the extended downstream region of path_seq */ 
  if (path->revcomp) {
    path_end_pos =  path_seq->n - path->jali[path->path_len-1] + path_seq->end;
    ext_end_pos  =  path_seq->n - (path->jali[path->path_len-1] - ALIGNMENT_EXT*3) + path_seq->end; 
    for(pos = path->jali[path->path_len-1]-1; pos >= path->jali[path->path_len-1] - ALIGNMENT_EXT*3; pos-=3) {
      seq_pos = path_seq->n - pos + path_seq->end;
      if(seq_pos > path_seq->n-2) {
        ext_end_pos = seq_pos-1;
        break;
      }
      
      amino = esl_gencode_GetTranslation(gcode,&path_seq->dsq[seq_pos]);
 
      if(esl_abc_XIsNonresidue(gm->abc, amino)) {
        ext_end_pos = seq_pos-1;
        break;
      } 
 
    }
  }
  else {
    path_end_pos = path->jali[path->path_len-1] - path_seq->start + 1;
    ext_end_pos  = (path->jali[path->path_len-1] +  ALIGNMENT_EXT*3) - path_seq->start + 1;
    //printf("ext_end_pos %d\n", ext_end_pos);
    for(pos = path->jali[path->path_len-1]+1; pos <= path->jali[path->path_len-1] +  ALIGNMENT_EXT*3; pos+=3) {

      seq_pos = pos - path_seq->start + 1;
      if(seq_pos > path_seq->n-2 ) {
        ext_end_pos = seq_pos-1;
        break;
      }

      amino = esl_gencode_GetTranslation(gcode,&path_seq->dsq[seq_pos]);

      if(esl_abc_XIsNonresidue(gm->abc, amino)) {
        ext_end_pos = seq_pos-1;
        break;
      }
    }
    //printf("ext_end_pos %d\n", ext_end_pos);
  }
 
  path_seq_len += (path_start_pos - ext_start_pos) + (ext_end_pos - path_end_pos);
 
  ESL_ALLOC(nuc_index, sizeof(int64_t) * (path_seq_len+2));
  ESL_ALLOC(nuc_dsq,   sizeof(ESL_DSQ) * (path_seq_len+2));

  nuc_index[0] = -1;
  nuc_dsq[0]   = eslDSQ_SENTINEL;
  seq_idx   = 1;
   
  /*Add upstream extension nucleotides */
  for(seq_pos = ext_start_pos; seq_pos < path_start_pos; seq_pos++) {
    nuc_index[seq_idx] = seq_pos;
    nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
    seq_idx++;
  }
  
  /* Copy spliced nucleotides into single sequence and track their original indicies */
  for (i = 0; i < path->path_len; i++) {
    if (path->revcomp) {
      for (pos = path->iali[i]; pos >= path->jali[i]; pos--) {
        seq_pos = path_seq->n - pos + path_seq->end;
        nuc_index[seq_idx] = seq_pos;
        nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;
      }
    }
    else {
      for (pos = path->iali[i]; pos <= path->jali[i]; pos++) {

        seq_pos = pos - path_seq->start + 1;
        nuc_index[seq_idx] = seq_pos;
        nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;
      }

        
    }
  }

  /*Add downstream extension nucleotides */
  for(seq_pos = path_end_pos+1; seq_pos <= ext_end_pos; seq_pos++) {
    nuc_index[seq_idx] = seq_pos;
    nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
    seq_idx++;
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
      path->frameshift = TRUE;
    //printf("seq_idx %d\n", seq_idx); 
      ESL_XEXCEPTION(eslFAIL, "STOP codon");
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
  amino_seq   = esl_sq_CreateDigitalFrom(gcode->aa_abc, path_seq->name, amino_dsq, amino_len,    NULL,NULL,NULL);
  nuc_seq     = esl_sq_CreateDigitalFrom(gcode->nt_abc, path_seq->name, nuc_dsq,   path_seq_len, NULL,NULL,NULL);

  pli->nuc_sq       = nuc_seq;
  pli->amino_sq     = amino_seq;
  pli->orig_nuc_idx = nuc_index;

  //p7_omx_SetDumpMode(stdout, pli->bwd, TRUE);
  /* Algin the splices Amino sequence */
  status = align_spliced_path(pli, om, gm, path_seq, gcode, fs_prob);

   /* Alignment failed */
  if(pli->hit == NULL || pli->hit->dcl->ad->exon_cnt == 1) {

    if(nuc_dsq   != NULL) free(nuc_dsq);
    if(amino_dsq != NULL) free(amino_dsq);
     
     return eslOK;
  }

  /* adjust the hit coords */
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
  //for(i = 0; i < path->path_len; i++)
   //  dom_score += path->signal_scores[i];

  dom_lnP   = esl_exp_logsurv(dom_score, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

  /* E-value adjusment */
  dom_lnP += log((float)db_nuc_cnt / (float)om->max_length);

  if ((pli->by_E && exp(dom_lnP) <= pli->E) || ((!pli->by_E) && dom_score >= pli->T)) {

    if ( path->path_len >  pli->hit->dcl->ad->exon_cnt) {
      /* Shift the path to start at the first hit that was inculded in the alignment
       * and end at the last hit that was included in the alignment */
      if(path->revcomp) {
        for(shift = 0; shift < path->path_len; shift++) {
          if(path->jali[shift] <= pli->hit->dcl->iali) break;
        }
      }
      else {
        for(shift = 0; shift < path->path_len; shift++) {
           if(path->jali[shift] >= pli->hit->dcl->iali) break;
        }
      }

      /* Ensure path will still contains an original hit after shifting */
      contains_orig = FALSE;
      for(exon = shift; exon < pli->hit->dcl->ad->exon_cnt ; exon ++ ) {
        if(path->node_id[exon] < graph->orig_N && path->node_id[exon] > -1) contains_orig = TRUE;
      }

      if(!contains_orig) {
        if(nuc_dsq   != NULL) free(nuc_dsq);
        if(amino_dsq != NULL) free(amino_dsq);
        
        return eslOK;
      }

      /* Shift path to start at frist hits that is in alignment */
      path->path_len = pli->hit->dcl->ad->exon_cnt;

      for(exon = 0; exon < path->path_len; exon++) {
        path->node_id[exon] = path->node_id[shift+exon];
        path->iali[exon]    = path->iali[shift+exon];
        path->jali[exon]    = path->jali[shift+exon];
        path->ihmm[exon]    = path->ihmm[shift+exon];
        path->jhmm[exon]    = path->jhmm[shift+exon];
      }
      path->iali[0] = pli->hit->dcl->iali;
      path->ihmm[0] = pli->hit->dcl->ihmm;

      for(exon = 1; exon < path->path_len; exon++) {
        if (path->iali[exon] > pli->hit->dcl->jali)
          break;
      }
 
      path->path_len =  pli->hit->dcl->ad->exon_cnt;

      path->jali[exon-1] = pli->hit->dcl->jali;
      path->jhmm[exon-1] = pli->hit->dcl->jhmm;
    }

    /* Find the first original hit in path to copy info*/
    i = 0;
    while( path->node_id[i] >= graph->orig_N || path->node_id[i] == -1) {
      pli->hit->dcl->ad->exon_orig[i] = FALSE;
      pli->hit->dcl->ad->exon_split[i] = FALSE;
      i++;
    }
    pli->hit->dcl->ad->exon_orig[i] = TRUE;
    if(i < path->path_len-2 && path->node_id[i] == path->node_id[i+1])
      pli->hit->dcl->ad->exon_split[i] = TRUE;
    else 
      pli->hit->dcl->ad->exon_split[i] = FALSE;


#ifdef HMMER_THREADS
  if(info->thread_id >= 0) pthread_mutex_lock(info->mutex);
#endif /*HMMER_THREADS*/

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
      if(remove_node < 0 || remove_node >= graph->orig_N) {
        pli->hit->dcl->ad->exon_orig[i] = FALSE;
        pli->hit->dcl->ad->exon_split[i] = FALSE; 
      }
      else {
        pli->hit->dcl->ad->exon_orig[i] = TRUE;
        if((i > 0                && path->node_id[i] == path->node_id[i-1]) ||
           (i < path->path_len-2 && path->node_id[i] == path->node_id[i+1]))
          pli->hit->dcl->ad->exon_split[i] = TRUE;
        else
          pli->hit->dcl->ad->exon_split[i] = FALSE;

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

#ifdef HMMER_THREADS
  if(info->thread_id >= 0) pthread_mutex_unlock(info->mutex);
#endif /*HMMER_THREADS*/

  }

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
  float bwdsc;
  p7_Backward(pli->amino_sq->dsq, pli->amino_sq->n, om, pli->fwd, pli->bwd, &bwdsc);

  if((status = p7_Decoding(om, pli->fwd, pli->bwd, pli->bwd)) == eslERANGE) goto ERROR;

  p7_OptimalAccuracy(om, pli->bwd, pli->fwd, &oasc);
  p7_OATrace        (om, pli->bwd, pli->fwd, tr);

  p7_trace_Index(tr);

  p7_splice_ComputeAliScores(hit->dcl, tr, pli->amino_sq->dsq, gm, pli->bg, fs_prob, FALSE);

  hit->dcl->tr = p7_trace_splice_Convert(tr, pli->orig_nuc_idx, &splice_cnt);

  hit->dcl->ad = p7_alidisplay_splice_Create(hit->dcl->tr, 0, om, target_seq, pli->amino_sq, hit->dcl->scores_per_pos, tr->sqfrom[0], splice_cnt);

  p7_Null2_ByExpectation(om, pli->bwd, null2);
  
  domcorrection = 0.;
  for (i = 1; i <= pli->amino_sq->n; i++) {
    domcorrection += logf(null2[pli->amino_sq->dsq[i]]);
  }

  /* Work around for underflow in the backward matrix */
  if( isnan(domcorrection) ) {
     
    p7_Null2_ByTrace(om, tr, tr->tfrom[0], tr->tto[0], pli->bwd , null2);

	domcorrection = 0.;
    for (i = 1; i <= pli->amino_sq->n; i++) {
      domcorrection += logf(null2[pli->amino_sq->dsq[i]]);
    }

  }
  hit->dcl->domcorrection = ESL_MAX(0.0, domcorrection);

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
p7_splice_AlignFrameshiftPath(SPLICE_GRAPH *graph, SPLICE_PATH *path, SPLICE_PIPELINE *pli, P7_TOPHITS *tophits, P7_FS_PROFILE *gm_fs, ESL_GENCODE *gcode, ESL_SQ *path_seq, int64_t db_nuc_cnt, SPLICE_WORKER_INFO *info)
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
  int          path_start_pos, path_end_pos;
  float        dom_bias;
  float        nullsc;
  float        dom_score;
  double       dom_lnP;
  int64_t     *nuc_index;
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

  /* Find out how may exteded nucletides are present in path_seq */
  if (path->revcomp) {
    path_start_pos = path_seq->n - path->downstream_spliced_nuc_start[0] + path_seq->end;
    path_end_pos   = path_seq->n - path->upstream_spliced_nuc_end[path->path_len] + path_seq->end;    
  }
  else {
    path_start_pos = path->downstream_spliced_nuc_start[0] - path_seq->start + 1;
    path_end_pos   = path->upstream_spliced_nuc_end[path->path_len] - path_seq->start + 1;  
  }

  path_seq_len += path_start_pos - 1;
  path_seq_len += path_seq->n - path_end_pos;

  ESL_ALLOC(nuc_index, sizeof(int64_t) * (path_seq_len+2));
  ESL_ALLOC(nuc_dsq,   sizeof(ESL_DSQ) * (path_seq_len+2));

  nuc_index[0] = -1;
  nuc_dsq[0]   = eslDSQ_SENTINEL;
  seq_idx   = 1;

  /*Add upstream extension nucleotides */
  for(seq_pos = 1; seq_pos < path_start_pos; seq_pos++) {
    nuc_index[seq_idx] = seq_pos;
    nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
    seq_idx++;
  } 

  /* Copy spliced nucleotides into single sequence and track their original indicies */
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

  /*Add downstream extension nucleotides */
  for(seq_pos = path_end_pos+1; seq_pos <= path_seq->n; seq_pos++) {
    nuc_index[seq_idx] = seq_pos;
    nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
    seq_idx++;
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
      ali_len += llabs(pli->hit->dcl->ad->exon_seq_ends[i] - pli->hit->dcl->ad->exon_seq_starts[i]) + 1;
  }
  else 
    ali_len = llabs(pli->hit->dcl->ad->sqto - pli->hit->dcl->ad->sqfrom) + 1;

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
        if(graph->split_orig_id[path->node_id[exon]] >= 0) contains_orig = TRUE;
      }

      if(!contains_orig) {
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
    while( graph->split_orig_id[path->node_id[i]] < 0) {

      pli->hit->dcl->ad->exon_orig[i] = FALSE;
      pli->hit->dcl->ad->exon_split[i] = path->split[i];
      i++;
    }

    pli->hit->dcl->ad->exon_orig[i] = TRUE;
    pli->hit->dcl->ad->exon_split[i] = path->split[i];

#ifdef HMMER_THREADS
  if(info->thread_id >= 0) pthread_mutex_lock(info->mutex);
#endif /*HMMER_THREADS*/

    replace_node = graph->split_orig_id[path->node_id[i]];
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
      remove_node = graph->split_orig_id[path->node_id[i]];
      if( remove_node < 0) {
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
#ifdef HMMER_THREADS
  if(info->thread_id >= 0) pthread_mutex_unlock(info->mutex);
#endif /*HMMER_THREADS*/
  }
  
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
 
  gxppfs = p7_gmx_fs_Create(gm_fs->M, pli->nuc_sq->n, pli->nuc_sq->n, p7P_CODONS);
  p7_Decoding_Frameshift(gm_fs, pli->gfwd, pli->gbwd, gxppfs);

  p7_OptimalAccuracy_Frameshift(gm_fs, gxppfs, pli->gbwd, &oasc);
  p7_OATrace_Frameshift(gm_fs, gxppfs, pli->gbwd, pli->gfwd, tr);   /* <tr>'s seq coords are offset by i-1, rel to orig dsq */

  p7_trace_fs_Index(tr);

  hit          = p7_hit_Create_empty();
  hit->dcl     = p7_domain_Create_empty();
  hit->dcl->tr = NULL;

  p7_splice_ComputeAliScores_fs(hit->dcl, tr, pli->nuc_sq, gm_fs, pli->bg, FALSE);

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

  hit->dcl->domcorrection = ESL_MAX(0.0, domcorrection);

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
p7_splice_RemoveHits(SPLICE_GRAPH *graph, int range_bound_min, int range_bound_max)
{

  int        i;
  int        hit_min, hit_max;
  int        overlap_min, overlap_max;
  int        overlap_len;
  P7_TOPHITS *th;
  P7_HIT     *hit;

  th = graph->th;

  for(i = 0; i < th->N; i++) {
    if(!graph->node_in_graph[i]) continue;

    hit = th->hit[i];

    hit_min = ESL_MIN(hit->dcl->iali, hit->dcl->jali);
    hit_max = ESL_MAX(hit->dcl->iali, hit->dcl->jali);
    overlap_min = ESL_MAX(hit_min, range_bound_min);
    overlap_max = ESL_MIN(hit_max, range_bound_max);
    overlap_len = overlap_max - overlap_min + 1;
    
    if(overlap_len > 0)  graph->node_in_graph[i] = FALSE;
  }

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
        tmp_edge->downstream_node_id = -1;
        graph->tot_edges--;
      } 
    }
  }

  return eslOK;

}

ESL_SQ* 
p7_splice_GetSubSequence(const ESL_SQFILE *seq_file, char* seqname, int64_t seq_min, int64_t seq_max, int revcomp, SPLICE_WORKER_INFO *info)
{

  
  ESL_SQ     *target_seq;
  ESL_SQFILE *tmp_file;

#ifdef HMMER_THREADS
  if(info && info->thread_id >= 0) pthread_mutex_lock(info->mutex);
#endif /*HMMER_THREADS*/
	
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

#ifdef HMMER_THREADS
  if(info && info->thread_id >= 0) pthread_mutex_unlock(info->mutex);
#endif /*HMMER_THREADS*/

  esl_sq_SetName(target_seq, seqname);
 
  if(revcomp)
   esl_sq_ReverseComplement(target_seq);

  esl_sq_Digitize(target_seq->abc, target_seq);

  return target_seq;
}

ESL_SQ* 
p7_splice_GetSplicedSequence(SPLICE_PATH *path, ESL_SQ *path_seq, int up_node, int down_node, int intron_include, int **remove_idx, int64_t **nuc_index)
{

  int a, r, s, x;
  int sub_path_len;
  int intron_len;
  int amino_gap;
  int seq_idx, seq_pos, true_idx;
  int     *r_idx;
  int64_t *n_idx;
  ESL_DSQ *splice_dsq;
  ESL_SQ  *splice_seq;
  int status;

  splice_dsq = NULL;
  splice_seq = NULL;
 
  r_idx = *remove_idx;
  n_idx = *nuc_index;

  /* Get the path_len of all nodes and introns (if included)*/
  sub_path_len = 0; 
  for(s = up_node; s < down_node; s++) { 

    sub_path_len += llabs(path->hits[s]->dcl->jali - path->hits[s]->dcl->iali)+1;

    if(intron_include) {
      intron_len = llabs(path->hits[s+1]->dcl->iali - path->hits[s]->dcl->jali) - 1;
      amino_gap =  ESL_MAX(0, path->hits[s+1]->dcl->ihmm - path->hits[s]->dcl->jhmm - 1);

      /* inron length */
      if(intron_len <= MAX_INTRON_INCL + (amino_gap*3)) //short intron - add full length
        sub_path_len += intron_len;
      else
        sub_path_len += MAX_INTRON_INCL + (amino_gap*3);
    }
    
  } 
  sub_path_len += llabs(path->hits[down_node]->dcl->jali - path->hits[down_node]->dcl->iali)+1;

  ESL_ALLOC(r_idx, sizeof(int)     * (down_node-up_node+1));
  ESL_ALLOC(n_idx,  sizeof(int64_t) * (sub_path_len+2));
  ESL_ALLOC(splice_dsq, sizeof(ESL_DSQ) * (sub_path_len+2));
  
  n_idx[0]   = -1;
  splice_dsq[0]  = eslDSQ_SENTINEL;
  seq_idx   = 1;

  seq_pos  = llabs(path->hits[up_node]->dcl->iali - path_seq->start) + 1;
  true_idx = path->hits[up_node]->dcl->iali;
  
  r = 0;
  for(s = up_node; s < down_node; s++) {
 
    /* Add upstream hit nucleotides */ 
    if(path->revcomp) {
      while(true_idx >= path->hits[s]->dcl->jali) {
        n_idx[seq_idx] = true_idx;
        splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;
        seq_pos++;
        true_idx--;
      }
    }
    else {
      while(true_idx <= path->hits[s]->dcl->jali) {
        n_idx[seq_idx] = true_idx;
        splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;
        seq_pos++;
        true_idx++;
      }
    } 

    intron_len = llabs(path->hits[s+1]->dcl->iali - path->hits[s]->dcl->jali) - 1;
    amino_gap =  ESL_MAX(0, path->hits[s+1]->dcl->ihmm - path->hits[s]->dcl->jhmm - 1); 
    
    /* Add full seqeunce gap if short */
    if(intron_len <= MAX_INTRON_INCL + (amino_gap*3)) {
      r_idx[r*2]   = -1;
      r_idx[r*2+1] = -1;

      if(path->revcomp) {
        while(true_idx > path->hits[s+1]->dcl->iali) {

          n_idx[seq_idx] = true_idx;
          splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
          seq_idx++;
          seq_pos++;
          true_idx--;
        }
      }
      else {
        while(true_idx < path->hits[s+1]->dcl->iali) {
        n_idx[seq_idx] = true_idx;
          splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
          seq_idx++;
          seq_pos++;
          true_idx++;
        }
      }
    } 
    /* Add partial sequence gap if long */
    else {
      /* Add first MAX_INTRON_INCL/2 nucleotides of seqeunce gap */
      for(x = 0; x < MAX_INTRON_INCL/2; x++) {
        if(path->revcomp) {
          n_idx[seq_idx] = true_idx;
          splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
          seq_idx++;
          seq_pos++;
          true_idx--;
        }
        else {
          n_idx[seq_idx] = true_idx;
          splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
          seq_idx++;
          seq_pos++;
          true_idx++;
        }
      }
      r_idx[r*2] = seq_idx-1;    

      /* Add degenerate nucleotides to cover amino gap */
      for(a = 0; a < amino_gap; a++) {
        n_idx[seq_idx] = -1;
        splice_dsq[seq_idx]   =  path_seq->abc->Kp-3;
        seq_idx++;

        n_idx[seq_idx] = -1;
        splice_dsq[seq_idx]   =  path_seq->abc->Kp-3;
        seq_idx++;

        n_idx[seq_idx] = -1;
        splice_dsq[seq_idx]   =  path_seq->abc->Kp-3;
        seq_idx++;
      }

      /* Skip middle portion of sequence gap  */
      if(path->revcomp) {
        while(true_idx > path->hits[s+1]->dcl->iali+MAX_INTRON_INCL/2) {
          seq_pos++;
          true_idx--;
        }
      }
      else {
        while(true_idx < path->hits[s+1]->dcl->iali-MAX_INTRON_INCL/2) {
          seq_pos++;
          true_idx++;
        }
      }
      r_idx[r*2+1] = seq_idx;

      /* Add last MAX_INTRON_INCL/2 nucleotides of seqeunce gap */
      if(path->revcomp) {
        while(true_idx > path->hits[s+1]->dcl->iali) {

          n_idx[seq_idx] = true_idx;
          splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
          seq_idx++;
          seq_pos++;
          true_idx--;
        }
      }
      else {
        while(true_idx < path->hits[s+1]->dcl->iali) {
          n_idx[seq_idx] = true_idx;
          splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
          seq_idx++;
          seq_pos++;
          true_idx++;
        }
      }  
    }
    r++;
  }
  /* Add final hit nucleotides */
  if(path->revcomp) {
    while(true_idx >= path->hits[down_node]->dcl->jali) {
      n_idx[seq_idx] = true_idx;
      splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
      seq_idx++;
      seq_pos++;
      true_idx--;
    }
  }
  else {
    while(true_idx <= path->hits[down_node]->dcl->jali) {
      n_idx[seq_idx] = true_idx;
      splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
      seq_idx++;
      seq_pos++;
      true_idx++;
    }
  }

  n_idx[seq_idx] = -1;
  splice_dsq[seq_idx]   = eslDSQ_SENTINEL;

  splice_seq   = esl_sq_CreateDigitalFrom(path_seq->abc, NULL, splice_dsq, seq_idx-1, NULL,NULL,NULL);
  splice_seq->start = path->hits[up_node]->dcl->iali;
  splice_seq->end   = path->hits[down_node]->dcl->jali;

  *nuc_index  = n_idx;
  *remove_idx = r_idx;
  
  if(splice_dsq != NULL) free(splice_dsq);

  return splice_seq;

  ERROR:
    if(r_idx != NULL) free(r_idx);
    if(n_idx  != NULL) free(n_idx);
    if(splice_dsq != NULL) free(splice_dsq);
    esl_sq_Destroy(splice_seq);   
    return NULL;
}

ESL_SQ* 
p7_splice_GetSplicedSequence2(SPLICE_PATH2 *path, ESL_SQ *path_seq, int up_node, int down_node, int intron_include, int **remove_idx, int64_t **nuc_index)
{

  int a, r, s, x;
  int sub_path_len;
  int intron_len;
  int amino_gap;
  int seq_idx, seq_pos, true_idx;
  int     *r_idx;
  int64_t *n_idx;
  ESL_DSQ *splice_dsq;
  ESL_SQ  *splice_seq;
  int status;

  splice_dsq = NULL;
  splice_seq = NULL;
 
  r_idx = *remove_idx;
  n_idx = *nuc_index;

  /* Get the path_len of all nodes and introns (if included)*/
  sub_path_len = 0; 
  for(s = up_node; s < down_node; s++) { 

    sub_path_len += llabs(path->jali[s] - path->iali[s])+1;

    if(intron_include) {
      intron_len = llabs(path->iali[s+1] - path->jali[s]) - 1;
      amino_gap =  ESL_MAX(0, path->ihmm[s+1] - path->jhmm[s] - 1);

      /* inron length */
      if(intron_len <= MAX_INTRON_INCL + (amino_gap*3)) //short intron - add full length
        sub_path_len += intron_len;
      else
        sub_path_len += MAX_INTRON_INCL + (amino_gap*3);
    }
    
  } 
  sub_path_len += llabs(path->jali[down_node] - path->iali[down_node])+1;

  ESL_ALLOC(r_idx, sizeof(int)     * (down_node-up_node+1));
  ESL_ALLOC(n_idx,  sizeof(int64_t) * (sub_path_len+2));
  ESL_ALLOC(splice_dsq, sizeof(ESL_DSQ) * (sub_path_len+2));
  
  n_idx[0]   = -1;
  splice_dsq[0]  = eslDSQ_SENTINEL;
  seq_idx   = 1;

  seq_pos  = llabs(path->iali[up_node] - path_seq->start) + 1;
  true_idx = path->iali[up_node];
  
  r = 0;
  for(s = up_node; s < down_node; s++) {
 
    /* Add upstream hit nucleotides */ 
    if(path->revcomp) {
      while(true_idx >= path->jali[s]) {
        n_idx[seq_idx] = true_idx;
        splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;
        seq_pos++;
        true_idx--;
      }
    }
    else {
      while(true_idx <= path->jali[s]) {
        n_idx[seq_idx] = true_idx;
        splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;
        seq_pos++;
        true_idx++;
      }
    } 

    intron_len = llabs(path->iali[s+1] - path->jali[s]) - 1;
    amino_gap =  ESL_MAX(0, path->ihmm[s+1] - path->jhmm[s] - 1); 
    
    /* Add full seqeunce gap if short */
    if(intron_len <= MAX_INTRON_INCL + (amino_gap*3)) {
      r_idx[r*2]   = 0;
      r_idx[r*2+1] = 0;

      if(path->revcomp) {
        while(true_idx > path->iali[s+1]) {

          n_idx[seq_idx] = true_idx;
          splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
          seq_idx++;
          seq_pos++;
          true_idx--;
        }
      }
      else {
        while(true_idx < path->iali[s+1]) {
        n_idx[seq_idx] = true_idx;
          splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
          seq_idx++;
          seq_pos++;
          true_idx++;
        }
      }
    } 
    /* Add partial sequence gap if long */
    else {
      /* Add first MAX_INTRON_INCL/2 nucleotides of seqeunce gap */
      for(x = 0; x < MAX_INTRON_INCL/2; x++) {
        if(path->revcomp) {
          n_idx[seq_idx] = true_idx;
          splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
          seq_idx++;
          seq_pos++;
          true_idx--;
        }
        else {
          n_idx[seq_idx] = true_idx;
          splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
          seq_idx++;
          seq_pos++;
          true_idx++;
        }
      }
      r_idx[r*2] = seq_idx-1;    

      /* Add degenerate nucleotides to cover amino gap */
      for(a = 0; a < amino_gap; a++) {
        n_idx[seq_idx] = -1;
        splice_dsq[seq_idx]   =  path_seq->abc->Kp-3;
        seq_idx++;

        n_idx[seq_idx] = -1;
        splice_dsq[seq_idx]   =  path_seq->abc->Kp-3;
        seq_idx++;

        n_idx[seq_idx] = -1;
        splice_dsq[seq_idx]   =  path_seq->abc->Kp-3;
        seq_idx++;
      }

      /* Skip middle portion of sequence gap  */
      if(path->revcomp) {
        while(true_idx > path->iali[s+1] + MAX_INTRON_INCL/2) {
          seq_pos++;
          true_idx--;
        }
      }
      else {
        while(true_idx < path->iali[s+1] - MAX_INTRON_INCL/2) {
          seq_pos++;
          true_idx++;
        }
      }
      r_idx[r*2+1] = seq_idx;

      /* Add last MAX_INTRON_INCL/2 nucleotides of seqeunce gap */
      if(path->revcomp) {
        while(true_idx > path->iali[s+1]) {

          n_idx[seq_idx] = true_idx;
          splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
          seq_idx++;
          seq_pos++;
          true_idx--;
        }
      }
      else {
        while(true_idx < path->iali[s+1]) {
          n_idx[seq_idx] = true_idx;
          splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
          seq_idx++;
          seq_pos++;
          true_idx++;
        }
      }  
    }
    r++;
  }
  /* Add final hit nucleotides */
  if(path->revcomp) {
    while(true_idx >= path->jali[down_node]) {
      n_idx[seq_idx] = true_idx;
      splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
      seq_idx++;
      seq_pos++;
      true_idx--;
    }
  }
  else {
    while(true_idx <= path->jali[down_node]) {
      n_idx[seq_idx] = true_idx;
      splice_dsq[seq_idx]   = path_seq->dsq[seq_pos];
      seq_idx++;
      seq_pos++;
      true_idx++;
    }
  }

  n_idx[seq_idx] = -1;
  splice_dsq[seq_idx]   = eslDSQ_SENTINEL;

  splice_seq   = esl_sq_CreateDigitalFrom(path_seq->abc, NULL, splice_dsq, seq_idx-1, NULL,NULL,NULL);
  splice_seq->start = path->iali[up_node];
  splice_seq->end   = path->jali[down_node];

  *nuc_index  = n_idx;
  *remove_idx = r_idx;
  
  if(splice_dsq != NULL) free(splice_dsq);

  return splice_seq;

  ERROR:
    if(r_idx != NULL) free(r_idx);
    if(n_idx  != NULL) free(n_idx);
    if(splice_dsq != NULL) free(splice_dsq);
    esl_sq_Destroy(splice_seq);   
    return NULL;
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
p7_splice_ComputeAliScores(P7_DOMAIN *dom, P7_TRACE *tr, ESL_DSQ *amino_dsq, const P7_PROFILE *gm, P7_BG *bg, float fs_prob, int do_bias)
{

  int status;
  int i, k, n;
  int z1, z2;
  int N;
  int first_i;
  float bias;

  fs_prob = 0.;
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
  first_i = i;

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
     else ESL_XEXCEPTION(eslFAIL, "Impossible state from p7_splice_ComputeAliScores()");
  }

  dom->aliscore = 0.0;
  for (n=0; n<N; n++) dom->aliscore += dom->scores_per_pos[n];

  if(do_bias) {
    p7_bg_SetLength(bg, i-first_i);
    p7_bg_FilterScore(bg, amino_dsq+first_i-1, i-first_i, &bias);
    dom->aliscore -= bias;
  }

  return eslOK;

  ERROR:
    return status;
}



int
p7_splice_ComputeAliScores_fs(P7_DOMAIN *dom, P7_TRACE *tr, ESL_SQ *nuc_sq, const P7_FS_PROFILE *gm_fs, P7_BG *bg, int do_bias)
{

  int a, i, k, c, n;
  int codon_idx;
  int indel;
  int amino;
  int z1, z2;
  int N;
  float bias;
  ESL_DSQ *nuc_dsq; 
  ESL_DSQ *amino_dsq;
  int status;

  nuc_dsq   = nuc_sq->dsq;
  amino_dsq = NULL;

  if(do_bias) {
    ESL_ALLOC(amino_dsq, sizeof(ESL_DSQ) * tr->N);
    amino_dsq[0] = eslDSQ_SENTINEL;
  }

  if(tr->ndom == 0) p7_trace_fs_Index(tr);

  for(z1 = 0;       z1 < tr->N; z1++) if(tr->st[z1] == p7T_M) break;
  for(z2 = tr->N-1; z2 >= 0;    z2--) if(tr->st[z2] == p7T_M || tr->st[z2] == p7T_D) break;
    
  N = z2 - z1 + 1;
  if(dom->scores_per_pos == NULL)
    ESL_ALLOC( dom->scores_per_pos, sizeof(float) * N);
  else
    ESL_REALLOC( dom->scores_per_pos, sizeof(float) * N);
  for (n=0; n<N; n++)  dom->scores_per_pos[n] = 0.0; 

  k  = tr->hmmfrom[0];
  n  = 0;
  a  = 1;
  while (z1<=z2) {
    
    if (tr->st[z1] == p7T_M) {
      i = tr->i[z1];
      c = tr->c[z1];
      if(c == 1) {
        if(esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i]))
          codon_idx = p7P_CODON1(nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_QC2;
        tr->fs++;
      }
      else if(c == 2) {
        if(esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i]))
          codon_idx = p7P_CODON2(nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_QC1;
        tr->fs++;
      }
      else if(c == 3) {
        if(esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-2]) &&
         esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i]))
          codon_idx = p7P_CODON3(nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_C;
        /* record stop codon */
        indel = p7P_INDEL(gm_fs, k, codon_idx);
        if(indel == p7P_XXx || indel == p7P_XxX || indel == p7P_xXX) tr->fs++; 
      }  
      else if(c == 4) {
        if(esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-3]) &&
         esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-2]) &&
         esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i]))
          codon_idx = p7P_CODON4(nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_QC1;
        tr->fs++;
      }
      else if(c == 5) {
        if(esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-4]) &&
         esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-3]) &&
         esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-2]) &&
         esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-1]) &&
         esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i]))
          codon_idx = p7P_CODON5(nuc_dsq[i-4], nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_QC2;
        tr->fs++;
      }
      
      dom->scores_per_pos[n] = p7P_MSC_CODON(gm_fs, k, codon_idx);
      amino = p7P_AMINO(gm_fs, k, codon_idx);
     
      if(do_bias) {
        amino_dsq[a] = amino;
        a++;
      }
      
      if (tr->st[z1-1] == p7T_I) 
        dom->scores_per_pos[n] += p7P_TSC(gm_fs, k-1, p7P_IM);
      else if (tr->st[z1-1] == p7T_D)
        dom->scores_per_pos[n] += p7P_TSC(gm_fs, k-1, p7P_DM);
      k++; z1++; n++;

      while(z1 < z2 && tr->st[z1] == p7T_M) {
        c = tr->c[z1];
        i = tr->i[z1];       
           
        if(c == 1) {
          if(esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i]))
            codon_idx = p7P_CODON1(nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_QC2;
          tr->fs++;
        }
        else if(c == 2) {
          if(esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-1]) &&
           esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i]))
            codon_idx = p7P_CODON2(nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_QC1;
          tr->fs++;
        }
        else if(c == 3) {
          if(esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-2]) &&
           esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-1]) &&
           esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i]))
            codon_idx = p7P_CODON3(nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_C;
          indel = p7P_INDEL(gm_fs, k, codon_idx);
          if(indel == p7P_XXx || indel == p7P_XxX || indel == p7P_xXX) tr->fs++;
        }
        else if(c == 4) {
          if(esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-3]) &&
           esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-2]) &&
           esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-1]) &&
           esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i]))
            codon_idx = p7P_CODON4(nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_QC1;
          tr->fs++;
        }
        else if(c == 5) {
          if(esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-4]) &&
           esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-3]) &&
           esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-2]) &&
           esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-1]) &&
           esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i]))
            codon_idx = p7P_CODON5(nuc_dsq[i-4], nuc_dsq[i-3], nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
          else
            codon_idx    = p7P_DEGEN_QC2;
          tr->fs++;
        }
             
        dom->scores_per_pos[n] = p7P_MSC_CODON(gm_fs, k, codon_idx);

        if(do_bias) {
          amino = p7P_AMINO(gm_fs, k, codon_idx);
          amino_dsq[a] = amino;
          a++;
        }
        
        dom->scores_per_pos[n] += p7P_TSC(gm_fs, k-1, p7P_MM);
        k++; z1++; n++;
      }
    }
    else if (tr->st[z1] == p7T_I) {
      if(do_bias) {
        if(esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-2]) &&
           esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i-1]) &&
           esl_abc_XIsCanonical(nuc_sq->abc, nuc_dsq[i]))
          codon_idx = p7P_CODON3(nuc_dsq[i-2], nuc_dsq[i-1], nuc_dsq[i]);
        else
          codon_idx    = p7P_DEGEN_C;   

        amino = p7P_AMINO(gm_fs, k, codon_idx);
        amino_dsq[a] = amino;
        a++;
      }

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
    else if (tr->st[z1] == p7T_P) {
    
    }
    else ESL_XEXCEPTION(eslFAIL, "Impossible state from p7_splice_ComputeAliScores_fs()");
  }

  dom->aliscore = 0.0;
  for (n=0; n<N; n++)  dom->aliscore += dom->scores_per_pos[n];  
  
  if(do_bias) {
    amino_dsq[a] = eslDSQ_SENTINEL;   
    p7_bg_SetLength(bg, a-1); 
    p7_bg_FilterScore(bg, amino_dsq, a-1, &bias); 

    dom->aliscore -= bias;
    free(amino_dsq);
  }
  
  return eslOK;

  ERROR:
    if(amino_dsq != NULL) free(amino_dsq);
    return status;
}


