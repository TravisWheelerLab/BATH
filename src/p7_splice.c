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


static int serial_loop(SPLICE_WORKER_INFO *info, P7_TOPHITS *tophits, SPLICE_SAVED_HITS *saved_hits);
#ifdef HMMER_THREADS
static int thread_loop(SPLICE_WORKER_INFO *info, P7_TOPHITS *tophits, SPLICE_SAVED_HITS *saved_hits, int infocnt);
static void* splice_thread(void *arg);
#endif /*HMMER_THREADS*/
static int add_split_exons(SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, P7_HIT **exons, const P7_HMM *hmm, const P7_FS_PROFILE *gm_fs, const P7_PROFILE *gm, const ESL_GENCODE *gcode, ESL_SQ *ali_seq, float orig_aliscore, int orig_node, int num_hits, int frameshift);
static SPLICE_EDGE** connect_split(SPLICE_PIPELINE *pli, P7_HIT **split_hits, const P7_HMM *hmm, const P7_FS_PROFILE *gm_fs, const P7_PROFILE *gm, const ESL_GENCODE *gcode, ESL_SQ *ali_seq, float orig_aliscore, int revcomp, int frameshift, int *num_hits); 
static int align_spliced_path (SPLICE_PIPELINE *pli, P7_OPROFILE *om, P7_PROFILE *gm, ESL_SQ *target_seq, ESL_GENCODE *gcode, float fs_prob);
static int align_spliced_path_frameshift (SPLICE_PIPELINE *pli, P7_FS_PROFILE *gm_fs, ESL_SQ *target_seq, ESL_GENCODE *gcode);
static int hit_upstream(P7_DOMAIN *upstream, P7_DOMAIN *downstream, int revcomp);
static int hits_spliceable(P7_DOMAIN *upstream, P7_DOMAIN *downstream);
static float hmm_overlap(P7_DOMAIN *dom_1, P7_DOMAIN *dom_2);
static float seq_overlap(P7_DOMAIN *dom_1, P7_DOMAIN *dom_2, int revcomp);


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
p7_splice_SpliceHits(P7_TOPHITS *tophits, SPLICE_SAVED_HITS *saved_hits, P7_HMM *hmm, P7_OPROFILE *om, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, ESL_GETOPTS *go, ESL_GENCODE *gcode, ESL_SQFILE *seq_file, int64_t db_nuc_cnt)
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
  
  /* Remove Duplicates from saved hits list */
  p7_splicehits_RemoveDuplicates(saved_hits);

  printf("\nQuery %s LENG %d\n",  gm->name, gm->M);
  fflush(stdout);

  /* Get the number of threads */
  ncpus = 0;
#ifdef HMMER_THREADS
  ncpus = ESL_MIN(esl_opt_GetInteger(go, "--cpu"), esl_threads_GetCPUCount());
#endif
  
  /* Intialize data for threads */
  infocnt = (ncpus == 0) ? 1 : ncpus-1;
  ESL_ALLOC(info, sizeof(*info) * infocnt);
  for (i = 0; i < infocnt; ++i)
  {
	info[i].hmm            = p7_hmm_Clone(hmm);
    info[i].om             = p7_oprofile_Clone(om);
    info[i].gm             = p7_profile_Clone(gm);
    info[i].gm_fs          = p7_profile_fs_Clone(gm_fs);
	info[i].pli            = p7_splicepipeline_Create(go, om->M, om->M * 3);
	info[i].pli->bg        = p7_bg_Create(om->abc);
	if (info[i].pli->do_biasfilter)
	  p7_bg_SetFilter(info[i].pli->bg, om->M, om->compo);
    info[i].sh             = saved_hits;
	info[i].tophits        = tophits;
    info[i].gcode          = gcode;
	info[i].seq_file       = seq_file;
	info[i].db_nuc_cnt     = db_nuc_cnt;
  }

  /* Begin splicing */
#ifdef HMMER_THREADS
  if (ncpus > 0)  thread_loop(info, tophits, saved_hits, infocnt); 
  else
#endif  /*HMMER_THREADS*/
    serial_loop (info, tophits, saved_hits); 
  
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
serial_loop (SPLICE_WORKER_INFO *info, P7_TOPHITS *tophits, SPLICE_SAVED_HITS *saved_hits) 
{

  int              g, h;
  int              num_graphs;
  int              revcomp, curr_revcomp;
  int              first, last;
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
        if((tophits->hit[h]->flags & p7_IS_REPORTED) || exp(tophits->hit[h]->sum_lnP) < 0.1) {

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
        if((tophits->hit[h]->flags & p7_IS_REPORTED) || exp(tophits->hit[h]->sum_lnP) < 0.1) {
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
    p7_splice_AddOriginals(graph, tophits);

//printf("ORIGINAL\n");
//p7_splicegraph_DumpHits(stdout, graph); 
//fflush(stdout);
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
thread_loop (SPLICE_WORKER_INFO *info, P7_TOPHITS *tophits, SPLICE_SAVED_HITS *saved_hits, int infocnt) 
{
  int              i, g, h;
  int              num_graphs;
  int              revcomp, curr_revcomp;
  int              first, last;
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
        if((tophits->hit[h]->flags & p7_IS_REPORTED) || exp(tophits->hit[h]->sum_lnP) < 0.1) {

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
        if((tophits->hit[h]->flags & p7_IS_REPORTED) || exp(tophits->hit[h]->sum_lnP) < 0.1) {
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

    p7_splice_AddOriginals(graph, tophits);

//printf("ORIGINAL\n");
//p7_splicegraph_DumpHits(stdout, graph);
//fflush(stdout);
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
p7_splice_AddOriginals(SPLICE_GRAPH *graph, const P7_TOPHITS *tophits)
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
    if(!(tophits->hit[i]->flags & p7_IS_REPORTED) && exp(tophits->hit[i]->sum_lnP) >= 0.1) continue;

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
	if(!(curr_hit->flags & p7_IS_REPORTED) && exp(curr_hit->sum_lnP) >= 0.1) continue;

    if(curr_hit->flags & p7_IS_REPORTED) graph->reportable[graph->th->N] = TRUE;
    else                                 graph->reportable[graph->th->N] = FALSE;

    graph->node_in_graph[graph->th->N]     = TRUE;
    graph->orig_hit_idx[graph->th->N] = i;
    
    p7_splicegraph_AddNode(graph, curr_hit); 
  }

  graph->orig_N  = graph->th->N;
  graph->seqname = graph->th->hit[0]->name;

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
  
  int                h, p, s;
  int                first, last;
  int                seq_min, seq_max;
  int                num_paths;
  int                orig_idx;
  int64_t            path_min, path_max;
  SPLICE_PATH       *path;
  SPLICE_PATH       **path_accumulator;
  ESL_SQ            *path_seq;
  SPLICE_GRAPH      *graph;
  SPLICE_SAVED_HITS *saved_hits;
  SPLICE_PIPELINE   *pli;
  P7_TOPHITS        *tophits;
  P7_HMM            *hmm;
  P7_OPROFILE       *om;
  P7_PROFILE        *gm;
  P7_FS_PROFILE     *gm_fs;
  ESL_GENCODE       *gcode;
  ESL_SQFILE        *seq_file;
  int                status;

  path             = NULL;
  path_accumulator = NULL;
  path_seq         = NULL;

  graph      = info->graph;
  saved_hits = info->sh;
  pli        = info->pli;
  tophits    = info->tophits;
  hmm        = info->hmm;
  om         = info->om;
  gm         = info->gm;
  gm_fs      = info->gm_fs;
  gcode      = info->gcode;
  seq_file   = info->seq_file;

  /* Find the first and last index of the current seqidx and strand in the saved hits list */
  first = 0;
  while(first < saved_hits->N && (saved_hits->srt[first]->seqidx != graph->seqidx || saved_hits->srt[first]->strand != graph->revcomp)) first++;
  last = first;
  while(last < saved_hits->N && saved_hits->srt[last]->seqidx == graph->seqidx && saved_hits->srt[last]->strand == graph->revcomp) last++;
  last--; 
 
  printf("\nQuery %s Target %s strand %c first %d last %d seqidx %d\n", gm->name, graph->seqname, (graph->revcomp ? '-' : '+'), first, last, graph->seqidx);
  fflush(stdout);
  
  p7_splice_RecoverHits(info, first, last);

// if(info->thread_id >= 0) pthread_mutex_lock(info->mutex);
//printf("RECOVER\n");
//p7_splicegraph_DumpHits(stdout, graph);
//fflush(stdout);
// if(info->thread_id >= 0) pthread_mutex_unlock(info->mutex);

  /* Create edges between original and recovered nodes */
  p7_splice_CreateEdges(graph);

  /* Allocate space for paths */
  ESL_ALLOC(path_accumulator, sizeof(SPLICE_PATH*) * graph->num_nodes);
  num_paths = 0;
   

  /* Build paths from orignal hit nodes and edge so that every node appears in one and only one path */
  num_paths = 0;
  path = p7_splicepath_GetBestPath_Unspliced(graph);
 
  while(path != NULL) {
//if(info->thread_id >= 0) pthread_mutex_lock(info->mutex);
//p7_splicepath_Dump(stdout,path);
//if(info->thread_id >= 0) pthread_mutex_unlock(info->mutex);
    path_accumulator[num_paths] = path;
    num_paths++;

    for(s = 0; s < path->path_len; s++) {
      graph->node_in_graph[path->node_id[s]] = FALSE;
    }
   
    /* Break edges that overlap the path so that paths do not intertwine */ 
    path_min = ESL_MIN(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]);
    path_max = ESL_MAX(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]);
    p7_splice_EnforceRangeBounds(graph, path_min, path_max);

    path = p7_splicepath_GetBestPath_Unspliced(graph);
  }

  /* Mark any remaining revocered hits as no in graph */
  for(h = graph->orig_N; h < graph->num_nodes; h++)    
    graph->node_in_graph[h] = FALSE;

  /* Align paths with spliced Viterbi to find missing exon and internal introns */
  for(p = 0; p < num_paths; p++) {

    p7_splice_FindExons(info, path_accumulator[p]); 
    p7_splicepath_Destroy(path_accumulator[p]);
  }

//if(info->thread_id >= 0) pthread_mutex_lock(info->mutex);
//p7_splicegraph_DumpHits(stdout, graph);
//fflush(stdout);
//if(info->thread_id >= 0) pthread_mutex_unlock(info->mutex);

  /* Create splice edges */    
  p7_splice_ConnectGraph(graph, pli, gm_fs, hmm, gcode, seq_file, info);
  /* Create spliced paths */
  num_paths = 0; 
  path = p7_splicepath_GetBestPath(graph); 

  while(path != NULL) {
//if(info->thread_id >= 0) pthread_mutex_lock(info->mutex);
//p7_splicepath_Dump(stdout,path);
//if(info->thread_id >= 0) pthread_mutex_unlock(info->mutex);  
    path_accumulator[num_paths] = path;
    num_paths++;

    for(s = 0; s < path->path_len; s++) {
      graph->node_in_graph[path->node_id[s]] = FALSE;
      orig_idx = graph->split_orig_id[path->node_id[s]];
      if(orig_idx < 0 || orig_idx >= graph->orig_N) continue;
		for(h = 0; h < graph->th->N; h++) {
        if(graph->split_orig_id[h] == orig_idx) 
          graph->node_in_graph[h] = FALSE;
      }
    }
    /*Break edge the overlap the path so that paths do not intertwine */
    path_min = ESL_MIN(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]);
    path_max = ESL_MAX(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]);
    p7_splice_RemoveHits(graph, path_min, path_max);
    p7_splice_EnforceRangeBounds(graph, path_min, path_max);

    path = p7_splicepath_GetBestPath(graph);
  }

  /* Align final paths */
  for(p = 0; p < num_paths; p++) {
 
    path = path_accumulator[p];

    if(path->path_len > 1) {
      seq_min = ESL_MIN(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]) - ALIGNMENT_EXT*3;
      seq_max = ESL_MAX(path->downstream_spliced_nuc_start[0], path->upstream_spliced_nuc_end[path->path_len]) + ALIGNMENT_EXT*3;
      
      path_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, path->revcomp, info);

      if(!path->frameshift)
        p7_splice_AlignPath(graph, path, pli, tophits, om, gm, gcode, path_seq, info->db_nuc_cnt, gm_fs->fs, info);
       
      if(path->frameshift)
        p7_splice_AlignFrameshiftPath (graph, path, pli, tophits, gm_fs, gcode, path_seq, info->db_nuc_cnt, info);
         
      esl_sq_Destroy(path_seq);
    }
   
    p7_splicepipeline_Reuse(pli);
    p7_splicepath_Destroy(path);

  }

  if(path_accumulator != NULL) free(path_accumulator);
  return eslOK;

  ERROR:
    if(path_accumulator != NULL) free(path_accumulator);
	return status;

}

/*  Function: p7_splice_RecoverHits
 *  Synopsis: Add revcoverd hits to graph
 *
 *  Purpose : Recover potential exons that were discarded by 
 *            the BATH forward filter but saved in <ifo->sh> 
 *            and add them to splice graph 
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_splice_RecoverHits(SPLICE_WORKER_INFO *info, int first, int last) 
{

  int i, j, s;
  int z1, z2;
  int gap_len;
  int seq_min, seq_max;
  int amino_len;
  int nuc_seq_idx;
  int amino;
  P7_DOMAIN         *tmp_dom;
  P7_HIT            *new_hit;
  P7_TOPHITS        *th;
  P7_HMM            *sub_hmm;
  P7_PROFILE        *sub_model;
  P7_GMX            *vit_mx;
  ESL_SQ            *hit_nuc_seq; 
  ESL_DSQ           *hit_aa_dsq;
  SPLICE_GRAPH      *graph;
  SPLICE_SAVED_HITS *sh;
  int  status;

  graph    = info->graph;
  sh       = info->sh;

  th = graph->th;

  p7_splicehits_AssignNodes(info->graph, info->sh, first, last);
   
  tmp_dom = p7_domain_Create_empty();

  for(i = first; i <= last; i++) {
    
    if(sh->srt[i]->duplicate) continue;
    if(sh->srt[i]->node_id >= 0) continue;
    
    tmp_dom->ihmm = sh->srt[i]->hmm_start;
    tmp_dom->jhmm = sh->srt[i]->hmm_end;    
    tmp_dom->iali = sh->srt[i]->seq_start;    
    tmp_dom->jali = sh->srt[i]->seq_end;

    /* Check if hit info is for a hit that is spliceable upstream or downstream of any hits int the graph*/
    for(j = 0; j < graph->orig_N; j++) {
      if(!graph->node_in_graph[j]) continue;
      if(hit_upstream(tmp_dom, th->hit[j]->dcl, graph->revcomp)) {
        if(graph->revcomp) gap_len = tmp_dom->jali - th->hit[j]->dcl->iali - 1;
        else               gap_len = th->hit[j]->dcl->iali - tmp_dom->jali - 1;
        if(gap_len > MAX_INTRON_LENG) continue;
      }
      else if(hit_upstream(th->hit[j]->dcl, tmp_dom, graph->revcomp)) {
        if(graph->revcomp) gap_len = th->hit[j]->dcl->jali - tmp_dom->iali - 1;
        else               gap_len = tmp_dom->iali - th->hit[j]->dcl->jali - 1;

        if(gap_len > MAX_INTRON_LENG) continue;
      }
      else continue;
            
      new_hit          = p7_hit_Create_empty();
      new_hit->dcl     = p7_domain_Create_empty();
      new_hit->dcl->tr = p7_trace_fs_Create();

      seq_min = ESL_MIN(tmp_dom->iali, tmp_dom->jali);
      seq_max = ESL_MAX(tmp_dom->iali, tmp_dom->jali); 
      hit_nuc_seq = p7_splice_GetSubSequence(info->seq_file, graph->seqname, seq_min, seq_max, graph->revcomp, info);

      amino_len = hit_nuc_seq->n/3;
      ESL_ALLOC(hit_aa_dsq, sizeof(ESL_DSQ) * (amino_len+2));
      hit_aa_dsq[0] = eslDSQ_SENTINEL;

      nuc_seq_idx = 1;

      for(s = 1; s <= amino_len; s++) {
        amino = esl_gencode_GetTranslation(info->gcode, &hit_nuc_seq->dsq[nuc_seq_idx]);
        hit_aa_dsq[s] = amino;
        nuc_seq_idx+=3;
      }
      hit_aa_dsq[amino_len+1] = eslDSQ_SENTINEL;    
      
      sub_hmm = p7_splice_GetSubHMM(info->hmm, tmp_dom->ihmm, tmp_dom->jhmm);
      sub_hmm->fs = 0.;
      
      sub_model = p7_profile_Create(sub_hmm->M, sub_hmm->abc);
      p7_ProfileConfig(sub_hmm, info->pli->bg, sub_model, amino_len, p7_UNILOCAL);
   
      vit_mx = p7_gmx_Create(sub_model->M, amino_len);
    
      p7_GViterbi(hit_aa_dsq, amino_len, sub_model, vit_mx, NULL);
      p7_GTrace(hit_aa_dsq, amino_len, sub_model, vit_mx, new_hit->dcl->tr);
      p7_splice_ComputeAliScores(new_hit->dcl, new_hit->dcl->tr, hit_aa_dsq, sub_model, info->pli->bg, info->hmm->fs, TRUE);      
      
      if(graph->revcomp)
        p7_trace_fs_Convert(new_hit->dcl->tr, hit_nuc_seq->start, 1);
      else
        p7_trace_fs_Convert(new_hit->dcl->tr, hit_nuc_seq->start, 1);

      for(z1 = 0; z1 < new_hit->dcl->tr->N; z1++)    if(new_hit->dcl->tr->st[z1] == p7T_M) break;
      for(z2 = new_hit->dcl->tr->N-1; z2 >= 0; z2--) if(new_hit->dcl->tr->st[z2] == p7T_M) break;
       
      new_hit->dcl->ihmm = tmp_dom->ihmm + new_hit->dcl->tr->k[z1] - 1;
      new_hit->dcl->jhmm = tmp_dom->ihmm + new_hit->dcl->tr->k[z2] - 1; 
      if(graph->revcomp) {
        new_hit->dcl->iali = new_hit->dcl->tr->i[z1]-2;
        new_hit->dcl->jali = new_hit->dcl->iali - (new_hit->dcl->tr->i[z2] - new_hit->dcl->iali);
      }
      else {
        new_hit->dcl->iali = new_hit->dcl->tr->i[z1] - 2;
        new_hit->dcl->jali = new_hit->dcl->tr->i[z2];
      }
      p7_splicegraph_AddNode(graph, new_hit);
      
      sh->srt[i]->node_id = graph->num_nodes-1;

      p7_hmm_Destroy(sub_hmm);
      p7_profile_Destroy(sub_model);
      p7_gmx_Destroy(vit_mx);
      esl_sq_Destroy(hit_nuc_seq);
      free(hit_aa_dsq);
      break;
    } 
  
  }     

  free(tmp_dom);
  return eslOK;

  ERROR:
	free(tmp_dom);
	free(hit_aa_dsq);
	return status;
}



/*  Function: p7_splice_CreateEdges
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
p7_splice_CreateEdges(SPLICE_GRAPH *graph) 
{
  int up, down;
  int seq_gap_len;
  P7_TOPHITS  *th;
  SPLICE_EDGE *edge;

  th = graph->th;

  for(up = 0; up < graph->num_nodes; up++) {
    
    for(down = 0; down < graph->num_nodes; down++) {
      if(up == down) continue;

      /* Are hits upstream/downstream */
      if(!hit_upstream(th->hit[up]->dcl, th->hit[down]->dcl, graph->revcomp)) continue;

      if(graph->revcomp)
        seq_gap_len = th->hit[up]->dcl->jali - th->hit[down]->dcl->iali - 1;
      else
        seq_gap_len = th->hit[down]->dcl->iali - th->hit[up]->dcl->jali - 1;        

      if(seq_gap_len > MAX_INTRON_LENG) continue;

      edge = p7_spliceedge_Create();
      edge->splice_score       = 0.;
      edge->upstream_node_id   = up;
      edge->downstream_node_id = down;

      /* If hits overlap, find the minimum lost socre to remove the overlap */
      p7_spliceedge_AliScoreEdge(edge, th->hit[up]->dcl, th->hit[down]->dcl); 

      edge->upstream_spliced_amino_end     = th->hit[up]->dcl->jhmm;
      edge->downstream_spliced_amino_start = th->hit[down]->dcl->ihmm; 
      edge->upstream_spliced_nuc_end       = th->hit[up]->dcl->jali;
      edge->downstream_spliced_nuc_start   = th->hit[down]->dcl->iali;
 
      p7_splicegraph_AddEdge(graph, edge);      
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
p7_splice_FindExons(SPLICE_WORKER_INFO *info, SPLICE_PATH *path)
{

  int e, s;
  int s1, s2;
  int num_hits;
  int path_min, path_max;
  int path_len;
  int seq_min, seq_max;
  int first_split, last_split;
  int split_exons;
  P7_HMM          *sub_hmm;  
  P7_HMM          *hmm;
  P7_PROFILE      *gm;
  P7_FS_PROFILE   *gm_fs;
  ESL_GENCODE     *gcode;
  ESL_SQ          *ali_seq;
  SPLICE_GRAPH    *graph;
  SPLICE_PIPELINE *pli;
  P7_HIT          **exons;

  graph = info->graph;
  pli   = info->pli;
  gcode = info->gcode;
  gm_fs = info->gm_fs;
  gm    = info->gm;
  hmm   = info->hmm;

  /* If there is only one hit in the path realign it with non-frameshift 
   * spliced viterbi to see if it contains more than one exon. If the orginal 
   * hit has frameshifts realign a second time with frameshift spliced viterbi */ 
  if(path->path_len == 1) {

    seq_min = ESL_MIN(path->hits[0]->dcl->iali, path->hits[0]->dcl->jali);
    seq_max = ESL_MAX(path->hits[0]->dcl->iali, path->hits[0]->dcl->jali);

    ali_seq = p7_splice_GetSubSequence(info->seq_file, graph->seqname, seq_min, seq_max, graph->revcomp, info);    
    
    sub_hmm     = p7_splice_GetSubHMM(hmm, path->hits[0]->dcl->ihmm, path->hits[0]->dcl->jhmm);
    sub_hmm->fs = 0.;

    /* Get all hits (exons) in path */    
    exons = p7_splice_AlignExons(sub_hmm, gm_fs, pli->bg, ali_seq, gcode, graph->revcomp, path->hits[0]->dcl->ihmm, &num_hits); 

    add_split_exons(graph, pli, exons, hmm, gm_fs, gm, gcode, ali_seq, path->hits[0]->dcl->aliscore, path->node_id[0], num_hits, FALSE);

    free(exons);

    if(path->hits[0]->dcl->tr->fs) {
      sub_hmm->fs = hmm->fs;
      
      exons = p7_splice_AlignExons(sub_hmm, gm_fs, pli->bg, ali_seq, gcode, graph->revcomp, path->hits[0]->dcl->ihmm, &num_hits);

      add_split_exons(graph, pli, exons, hmm, gm_fs, gm, gcode, ali_seq, path->hits[0]->dcl->aliscore, path->node_id[0], num_hits, TRUE);
      
      free(exons); 
    }

    esl_sq_Destroy(ali_seq);
    p7_hmm_Destroy(sub_hmm);    

    return eslOK;
  }

  path_len = 0; 
  s1 = 0;
  s2 = 0;
  while(s1 < path->path_len) {
    while(path_len <  MAX_INTRON_LENG) {
      s2++;
      if(s2 == path->path_len) break;
      
      path_min = ESL_MIN(path->hits[s1]->dcl->iali, path->hits[s2]->dcl->jali);
      path_max = ESL_MAX(path->hits[s1]->dcl->iali, path->hits[s2]->dcl->jali);
      path_len = path_max - path_min + 1;
    }

    path_min = ESL_MIN(path->hits[s1]->dcl->iali, path->hits[s2-1]->dcl->jali);
    path_max = ESL_MAX(path->hits[s1]->dcl->iali, path->hits[s2-1]->dcl->jali);
    path_len = path_max - path_min + 1;

    ali_seq = p7_splice_GetSubSequence(info->seq_file, graph->seqname, path_min, path_max, graph->revcomp, info);
    sub_hmm     = p7_splice_GetSubHMM(hmm, path->hits[s1]->dcl->ihmm, path->hits[s2-1]->dcl->jhmm);
    sub_hmm->fs = 0.;
    exons = p7_splice_AlignExons(sub_hmm, gm_fs, pli->bg, ali_seq, gcode, graph->revcomp, path->hits[s1]->dcl->ihmm, &num_hits);

    first_split = 0;
    last_split  = 0;
    e = 0;
    while ( e < num_hits) {
      split_exons = -1;

      /* Check if the current exon overlap suffcenitly with any of the original hits */
      s = s1;
      while( s < s2) {
        if(hmm_overlap(exons[e]->dcl, path->hits[s]->dcl) > 0.5 &&
           seq_overlap(exons[e]->dcl, path->hits[s]->dcl, graph->revcomp) > 0.5) {
          first_split = e;

          /*Find any additional hits that ovelap with the same orignal hit */
          while(e+1 < num_hits && hmm_overlap(exons[e+1]->dcl, path->hits[s]->dcl) > 0.5 &&
                seq_overlap(exons[e+1]->dcl, path->hits[s]->dcl, graph->revcomp) > 0.5) {
            e++;
          }
          last_split = e;
     
          /* Check that ant splits of original hits are better scoring than the 
           * unsplit alternaitve and add them to the graph */
          split_exons = last_split - first_split + 1;      
          
          add_split_exons(graph, pli, exons+first_split, hmm, gm_fs, gm, gcode, ali_seq, path->hits[s]->dcl->aliscore, path->node_id[s], split_exons, FALSE);
      
          s = s2;
        }
        s++;    
      }
    
      /* If the current exon does not overlap with any orignal hits add it to graph */
      if(split_exons == -1) 
        p7_splicegraph_AddNode(graph, exons[e]);    
     
      e++;
    }

    free(exons);

    if(path->frameshift) {
      sub_hmm->fs = hmm->fs;

      exons = p7_splice_AlignExons(sub_hmm, gm_fs, pli->bg, ali_seq, gcode, graph->revcomp, path->hits[s1]->dcl->ihmm, &num_hits);

      first_split = 0;
      last_split  = 0;
      e = 0;
      while ( e < num_hits) {
  
        split_exons = -1;
  
        /* Check if the current exon overlap suffcenitly with any of the original hits */
        s = s1;
        while( s < s2) {
          if(hmm_overlap(exons[e]->dcl, path->hits[s]->dcl) > 0.5 &&
             seq_overlap(exons[e]->dcl, path->hits[s]->dcl, graph->revcomp) > 0.5) {
            first_split = e;
  
            /*Find any additional hits that ovelap with the same orignal hit */
            while(e+1 < num_hits && hmm_overlap(exons[e+1]->dcl, path->hits[s]->dcl) > 0.5 &&
                  seq_overlap(exons[e+1]->dcl, path->hits[s]->dcl, graph->revcomp) > 0.5) {
              e++;
            }
            last_split = e;
  
            /* Check that ant splits of original hits are better scoring than the unsplit alternaitve and add them to the graph */
            split_exons = last_split - first_split + 1;
            add_split_exons(graph, pli, exons+first_split, hmm, gm_fs, gm, gcode, ali_seq, path->hits[s]->dcl->aliscore, path->node_id[s], split_exons, FALSE);
  
            s = s2;
          }
          s++;
        }
  
        /* If the current exon does not overlap with any orignal hits add it to graph */
        if(split_exons == -1)
          p7_splicegraph_AddNode(graph, exons[e]);
  
        e++;
      }
  
      free(exons);

    }

    s1 = s2;
    esl_sq_Destroy(ali_seq);
    p7_hmm_Destroy(sub_hmm);
  }
 
  return eslOK;

}


P7_HIT**
p7_splice_AlignExons(P7_HMM *sub_hmm, const P7_FS_PROFILE *gm_fs, P7_BG *bg, ESL_SQ *ali_seq, const ESL_GENCODE *gcode, int revcomp, int hmm_start, int *num_exons)
{
  int         i, y, z;
  int         z1, z2;
  int         intron_cnt;
  int         exon_cnt;
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
  p7_ProfileConfig_fs(sub_hmm, bg, gcode, sub_fs_model, ali_seq->n, p7_UNIGLOBAL);

  vit_mx = p7_gmx_fs_Create(sub_fs_model->M, ali_seq->n, ali_seq->n, p7P_SPLICE);
  tr = p7_trace_fs_Create();

  if(sub_hmm->fs > 0.0) {
    p7_sp_fs_semiglobal_Viterbi(ali_seq->dsq, gcode, ali_seq->n, sub_fs_model, vit_mx);
    p7_sp_fs_semiglobal_VTrace(ali_seq->dsq, ali_seq->n, gcode, sub_fs_model, vit_mx, tr);
  }
  else {
    p7_sp_trans_semiglobal_Viterbi(ali_seq->dsq, gcode, ali_seq->n, sub_fs_model, vit_mx);
    p7_sp_trans_semiglobal_VTrace(ali_seq->dsq, ali_seq->n, gcode, sub_fs_model, vit_mx, tr);
  }

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

  exon_cnt = 0;
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
      
      if(revcomp) {
        iali = ali_seq->n - tr->i[y] + ali_seq->end + tr->c[y] - 1;
        jali = ali_seq->n - tr->i[z] + ali_seq->end;
      }
      else {
        iali = ali_seq->start + tr->i[y] - tr->c[y];
        jali = ali_seq->start + tr->i[z] - 1;
      }
      
      /* Create new hit and  set ihmm and iali coords*/
      new_hit          = p7_hit_Create_empty();
      new_hit->dcl     = p7_domain_Create_empty();
      new_hit->dcl->tr = p7_trace_fs_Create();

      new_hit->dcl->ihmm = ihmm;
      new_hit->dcl->jhmm = jhmm;
      new_hit->dcl->iali = iali;
      new_hit->dcl->jali = jali;
       //printf("ihmm %d jhmm %d iali %d jali %d\n", ihmm, jhmm, iali, jali); 
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
  p7_gmx_Destroy(vit_mx);
  p7_trace_fs_Destroy(tr);

  return ret_hits;


  ERROR:
   if(ret_hits != NULL) free(ret_hits);
   return NULL;

}


int
add_split_exons(SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, P7_HIT **exons, const P7_HMM *hmm, const P7_FS_PROFILE *gm_fs, const P7_PROFILE *gm, const ESL_GENCODE *gcode, ESL_SQ *ali_seq, float orig_aliscore, int orig_node, int num_hits, int frameshift)
{

  int h;
  SPLICE_EDGE **edges;

  if(num_hits == 1) { 
   
    p7_splicegraph_AddNode(graph, exons[0]);
    graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[orig_node];
    if(orig_node < graph->orig_N) graph->split_orig_id[graph->num_nodes-1] = orig_node;
  }
  else if(num_hits > 1) {
    
    edges = connect_split(pli, exons, hmm, gm_fs, gm, gcode, ali_seq, orig_aliscore, graph->revcomp, frameshift, &num_hits);
 
    if(num_hits == 1) {
      
      p7_splicegraph_AddNode(graph, exons[0]);
      graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[orig_node];
      if(orig_node < graph->orig_N) graph->split_orig_id[graph->num_nodes-1] = orig_node;
       
      free(edges);
    }
    else if (num_hits > 1) {
      
      p7_splicegraph_AddNode(graph, exons[0]);
      graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[orig_node];
      if(orig_node < graph->orig_N) graph->split_orig_id[graph->num_nodes-1] = orig_node; 
     
      for(h = 1; h < num_hits; h++) {
     
        p7_splicegraph_AddNode(graph, exons[h]); 
        graph->orig_hit_idx[graph->num_nodes-1] = graph->orig_hit_idx[orig_node];
        if(orig_node < graph->orig_N) graph->split_orig_id[graph->num_nodes-1] = orig_node;

        edges[h-1]->upstream_node_id   = graph->num_nodes-2;
        edges[h-1]->downstream_node_id = graph->num_nodes-1; 
        p7_splicegraph_AddEdge(graph, edges[h-1]);
      }

      free(edges);

    }
  }

  return eslOK;
}

SPLICE_EDGE**
connect_split(SPLICE_PIPELINE *pli, P7_HIT **split_hits, const P7_HMM *hmm, const P7_FS_PROFILE *gm_fs, const P7_PROFILE *gm, const ESL_GENCODE *gcode, ESL_SQ *ali_seq, float orig_aliscore, int revcomp, int frameshift, int *num_hits) 
{

  int k, s, x, z;
  int hit_cnt;
  int edge_found;
  int seq_i, seq_j;
  int seq_len;
  int tmp_iali;
  P7_HMM        *sub_hmm; 
  P7_FS_PROFILE *sub_fs_model;
  P7_GMX        *vit_mx; 
  P7_TRACE      *tr;    
  SPLICE_EDGE   **split_edges;
  P7_DOMAIN      *merged_dom;
  int status;

  sub_hmm      = NULL;
  sub_fs_model = NULL;
  vit_mx       = NULL;
  split_edges  = NULL;

  hit_cnt = *num_hits;
  
  ESL_ALLOC(split_edges, sizeof(SPLICE_EDGE*) * (hit_cnt-1));

  /* Splice all split hits */
  for(s = 1; s < hit_cnt; s++) { 
    split_edges[s-1] = p7_spliceedge_ConnectNodes(pli, split_hits[s-1]->dcl, split_hits[s]->dcl, gm_fs, hmm, gcode, ali_seq, revcomp);        
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
 
	/* If the edge does not exist or the sum of the split hits and edge score is less than the orignal hit score, merge the split hits */
    if(split_edges[s] == NULL || split_hits[s]->dcl->aliscore + split_hits[s+1]->dcl->aliscore + split_edges[s]->splice_score < orig_aliscore) {
      
      sub_hmm = p7_splice_GetSubHMM(hmm, split_hits[s]->dcl->ihmm, split_hits[s+1]->dcl->jhmm);  
      if(!frameshift) sub_hmm->fs = 0.;
      
      seq_i = split_hits[s]->dcl->tr->sqfrom[0];
      seq_j = split_hits[s+1]->dcl->tr->sqto[0]; 
      seq_len = seq_j - seq_i + 1;  
//printf("seq_i %d\n", seq_i );     
      sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
      ESL_REALLOC(sub_fs_model->tsc,  sizeof(float)   * (sub_hmm->M+1) * p7P_NTRANS);
      p7_ProfileConfig_fs(sub_hmm, pli->bg, gcode, sub_fs_model, seq_len, p7_UNIGLOBAL);

      vit_mx = p7_gmx_fs_Create(sub_fs_model->M, seq_len, seq_len, p7P_CODONS);
      tr     = p7_trace_fs_Create();

      p7_fs_semiglobal_Viterbi(ali_seq->dsq+seq_i-1, gcode, seq_len, sub_fs_model, vit_mx);  
      p7_fs_semiglobal_VTrace(ali_seq->dsq+seq_i-1, seq_len, gcode, sub_fs_model, vit_mx, tr);

      /* shift tr->i to full seq coords */ 
      for (z = 0; z < tr->N; z++)
        if(tr->i[z] >= 0 ) tr->i[z] += seq_i-1;

      /* adjust trace k coords to full model */
      for(k = 0; k < tr->N; k++)
        if(tr->k[k] > 0) tr->k[k] += split_hits[s]->dcl->ihmm-1;

	  merged_dom = p7_domain_Create_empty();
      merged_dom->scores_per_pos = NULL;

      p7_splice_ComputeAliScores_fs(merged_dom, tr, ali_seq, gm_fs, pli->bg, TRUE);
      /* If the edge does not exist or the sum of the split hits and edge score is less than the merged hit score, replace the split hit with the merged hit */   
	  if(split_edges[s] == NULL || split_hits[s]->dcl->aliscore + split_hits[s+1]->dcl->aliscore + split_edges[s]->splice_score < merged_dom->aliscore) {
        p7_trace_fs_Destroy(split_hits[s]->dcl->tr);
        split_hits[s]->dcl->tr = tr;

        free(split_hits[s]->dcl->scores_per_pos);
        split_hits[s]->dcl->scores_per_pos = merged_dom->scores_per_pos; 
        merged_dom->scores_per_pos = NULL; 

        /* Set new hmm and ali coords */
        split_hits[s]->dcl->ihmm = tr->hmmfrom[0];
        split_hits[s]->dcl->jhmm = tr->hmmto[0];
        
        tmp_iali = split_hits[s]->dcl->iali;
        if(revcomp) {
          split_hits[s]->dcl->iali = tmp_iali - tr->sqfrom[0]  + seq_i;
          split_hits[s]->dcl->jali = tmp_iali - tr->sqto[0]    + seq_i;
        }
        else {
          split_hits[s]->dcl->iali = tmp_iali + tr->sqfrom[0] - seq_i;
          split_hits[s]->dcl->jali = tmp_iali + tr->sqto[0]   - seq_i;
        } 

        /*Destroy downstream hit */
        p7_trace_fs_Destroy(split_hits[s+1]->dcl->tr);
        free(split_hits[s+1]->dcl->scores_per_pos);
        p7_hit_Destroy(split_hits[s+1]);

        /* shift hits */
        for(x = s+1; x < hit_cnt-1; x++)
          split_hits[x] = split_hits[x+1];
        hit_cnt--;

        /*If there were pervious edges, they can no longer be trusted and we need to start over */
        if(hit_cnt > 1) { 
          for(x = 0; x < hit_cnt; x++ ) {
            if(split_edges[x] != NULL) free(split_edges[x]);
          }
          free(split_edges);
          split_edges = connect_split(pli, split_hits, hmm, gm_fs, gm, gcode, ali_seq, orig_aliscore, revcomp, frameshift, &hit_cnt); 
          s = 0;
        }        
        else {
          if(split_edges[s] != NULL) free(split_edges[s]);
          for(x = s; x < hit_cnt-1; x++)
            split_edges[x] = split_edges[x+1];
        }
      }
      else {
        p7_trace_fs_Destroy(tr);
        s++;
      }
      p7_domain_Destroy(merged_dom);
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
  dom_1_hmm_len = dom_1->jhmm - dom_1->ihmm - 1;

  return (float) overlap_hmm_len / (float) dom_1_hmm_len; 

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



int
p7_splice_ConnectGraph(SPLICE_GRAPH *graph, SPLICE_PIPELINE *pli, const P7_FS_PROFILE *gm_fs, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQFILE *seq_file, SPLICE_WORKER_INFO *info) 
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
       
      /*Nodes that are split from the same original hit will already have edges */
      if(graph->split_orig_id[up] != -1 && graph->split_orig_id[up] == graph->split_orig_id[down]) continue;
     
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
      splice_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, graph->revcomp, info);
  
     // printf("up %d down %d \n", up+1, down+1);
      edge = p7_spliceedge_ConnectNodes(pli, th->hit[up]->dcl, th->hit[down]->dcl, gm_fs, hmm, gcode, splice_seq, graph->revcomp);
       
      if(edge != NULL) {
        //printf("edge->splice_score %f\n", edge->splice_score);     
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
hits_spliceable(P7_DOMAIN *upstream, P7_DOMAIN *downstream) {
  

  /* Are nodes close enough on hmm coords */ 
  if(upstream->jhmm + MAX_AMINO_EXT < downstream->ihmm) return FALSE;

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

  /*Add upstream extention nucleotides */
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
  //      if (seq_idx == 772) printf("step %d pos %d seq_idx %d\n", i+1, pos, seq_idx);
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

  /*Add downstream extention nucleotides */
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

  /*Add upstream extention nucleotides */
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

  /*Add downstream extention nucleotides */
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
        graph->num_edges--;
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

#ifdef HMMER_THREADS
  if(info && info->thread_id >= 0) pthread_mutex_unlock(info->mutex);
#endif /*HMMER_THREADS*/

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
  
  N = tr->tto[0] - tr->tfrom[0] - 1;
  if(dom->scores_per_pos == NULL)
    ESL_ALLOC( dom->scores_per_pos, sizeof(float) * N);
  else
    ESL_REALLOC( dom->scores_per_pos, sizeof(float) * N);
  for (n=0; n<N; n++)  dom->scores_per_pos[n] = 0.0; 

  k  = tr->hmmfrom[0];
  z1 = tr->tfrom[0]+1;
  while(tr->st[z1] != p7T_M) z1++;
  z2 = tr->tto[0];
  while(tr->st[z2] != p7T_M) z2--;
  n  = 0;
  a  = 1;
  while (z1<z2) {
    
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


