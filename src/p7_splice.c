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
static float hmm_overlap(SPLICE_PATH *path1, int step1, SPLICE_PATH *path2, int step2);
static float seq_overlap(SPLICE_PATH *path1, int step1, SPLICE_PATH *path2, int step2); 


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
ncpus = ESL_MIN(ncpus,4);
  /* Intialize data for threads */
  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(*info) * infocnt);

  for (i = 0; i < infocnt; ++i)
  {
	info[i].hmm            = p7_hmm_Clone(hmm);
    info[i].om             = p7_oprofile_Clone(om);
    info[i].gm             = p7_profile_Clone(gm);
    info[i].gm_fs          = p7_profile_fs_Clone(gm_fs);
	info[i].pli            = p7_splicepipeline_Create(go, 100, 300);
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
        if((tophits->hit[h]->flags & p7_IS_REPORTED) || exp(tophits->hit[h]->sum_lnP) < info->pli->F3) {

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
        if((tophits->hit[h]->flags & p7_IS_REPORTED) || exp(tophits->hit[h]->sum_lnP) < info->pli->F3) {
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
    p7_splice_AddAnchors(info, graph, tophits);
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
        if((tophits->hit[h]->flags & p7_IS_REPORTED) || exp(tophits->hit[h]->sum_lnP) < info->pli->F3) {

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
        if((tophits->hit[h]->flags & p7_IS_REPORTED) || exp(tophits->hit[h]->sum_lnP) < info->pli->F3) {
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

    p7_splice_AddAnchors(info, graph, tophits);
    p7_splice_AddSeeds(graph, seed_hits);

  }
 
  /* Initialize mutex */
  if (pthread_mutex_init(&mutex, NULL)) { status = eslFAIL; goto ERROR; }

  /* Allocate space for threads */
  ESL_ALLOC(threads, sizeof(pthread_t) * infocnt);

  /* Add data needed for multi-threading and spin up threads to splice graphs */
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


/*  Function: p7_splice_AddAnchors
 *  Synopsis: Add BATH hits to graph
 *
 *  Purpose : Find all hits in <tophits> that match the seqidx 
 *            and revcomp of the <graph> and add them as nodes.
 *
 * Returns:   <eslOK>.
 *
 */
int 
p7_splice_AddAnchors(SPLICE_WORKER_INFO *info, SPLICE_GRAPH *graph, const P7_TOPHITS *tophits)
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
     
    if ((curr_hit->flags & p7_IS_DUPLICATE)) continue;
    if(!(curr_hit->flags & p7_IS_REPORTED) && exp(curr_hit->sum_lnP) >= info->pli->F3) continue;
     
    hit_cnt++;
  }
  
  p7_splicegraph_CreateNodes(graph, hit_cnt);    
  
  /*Add all hits from current sequence and strand to graph*/
  for (i = 0; i < tophits->N; i++) {

    curr_hit = tophits->hit[i];

    if (curr_hit->seqidx != graph->seqidx) continue;

    if (graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
    if ((!graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;

    if ((curr_hit->flags & p7_IS_DUPLICATE)) continue;
	if(!(curr_hit->flags & p7_IS_REPORTED) && exp(curr_hit->sum_lnP) >= info->pli->F3) continue;

    if(curr_hit->flags & p7_IS_REPORTED) graph->reportable[graph->th->N] = TRUE;
    else                                 graph->reportable[graph->th->N] = FALSE;

    graph->node_in_graph[graph->th->N] = TRUE;
    graph->orig_hit_idx[graph->th->N]  = i;
 
    p7_splicegraph_AddNode(graph, curr_hit); 

    graph->split_orig_id[graph->th->N-1] = graph->th->N-1;
  }

  graph->anchor_N = graph->recover_N = graph->th->N;
  graph->seqname = graph->th->hit[0]->name;

  return eslOK;
 
}

/* Only add seed that are upstream of one original hit and downstream from another */
int 
p7_splice_AddSeeds(SPLICE_GRAPH *graph, const P7_TOPHITS *seed_hits)
{

  int     i;
  int     h1, h2;
  int     gap_len;
  P7_HIT *curr_hit;
  P7_TOPHITS *th;

  th = graph->th;

  if(graph->anchor_N < 2) return eslOK;

  for (i = 0; i < seed_hits->N; i++) {

    curr_hit = &(seed_hits->unsrt[i]);

    /* Dont add seeds that didn't pass forward */
    if(!curr_hit->dcl->is_reported) continue;

    /* Is the seed hit on the same sequence and strand as the graph */
    if (curr_hit->seqidx != graph->seqidx) continue;
    if (graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
    if ((!graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;

    for(h1 = 0; h1 < graph->anchor_N; h1++) {
      /* Is the seed hit upstream of an original hit */
      if(p7_splice_HitUpstream(curr_hit->dcl, th->hit[h1]->dcl, graph->revcomp)) {
       
        if(graph->revcomp) gap_len = curr_hit->dcl->jali - th->hit[h1]->dcl->iali - 1;
        else               gap_len = th->hit[h1]->dcl->iali - curr_hit->dcl->jali - 1;
        if(gap_len > MAX_INTRON_LENG) continue;
       
        for(h2 = 0; h2 < graph->anchor_N; h2++) {
          if(h2 == h1) continue;

          /* Is the seed hit downstream of an original hit */
          if(p7_splice_HitUpstream(th->hit[h2]->dcl, curr_hit->dcl, graph->revcomp)) {
              
             if(graph->revcomp) gap_len = th->hit[h2]->dcl->jali - curr_hit->dcl->iali - 1;
             else               gap_len = curr_hit->dcl->iali - th->hit[h2]->dcl->jali - 1;
             if(gap_len > MAX_INTRON_LENG) continue;

             /* Re-purpose domain "is_included" for seed hits added to graph */
             curr_hit->dcl->is_included = TRUE;
             p7_splicegraph_AddNode(graph, curr_hit);
             h1 = graph->anchor_N;
             h2 = graph->anchor_N;
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
 *            find spice boundaries and align the spliced sequence
 *
 * Returns:   <eslOK> on success. 
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int 
p7_splice_SpliceGraph(SPLICE_WORKER_INFO *info) 
{
  
  int                h, s;
  int                success;
  int64_t            seq_min, seq_max;
  int64_t            hit_min, hit_max;
  int64_t            node_min, node_max;
  int64_t            path_min, path_max;
  SPLICE_PATH       *final_path;
  SPLICE_PATH       *orig_path;
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
  p7_splice_CreateUnsplicedEdges(graph, gm);
 
  /* Build paths from orignal hit nodes and edge so that every node appears in one and only one path */
  orig_path = p7_splicepath_GetBestPath_Unspliced(graph);
//p7_splicegraph_DumpGraph(stdout, graph);

  while(orig_path != NULL) {

    success = FALSE;

//TODO
    p7_splice_ExtendPath(info->seeds, gm, orig_path, graph);
    seq_min = ESL_MIN(orig_path->iali[0], orig_path->jali[orig_path->path_len-1]) - ALIGNMENT_EXT;
    seq_max = ESL_MAX(orig_path->iali[0], orig_path->jali[orig_path->path_len-1]) + ALIGNMENT_EXT;
    path_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, orig_path->revcomp, info);
    
//if(info->thread_id >= 0) pthread_mutex_lock(info->mutex);
//printf("FIRST PATH \n");
//p7_splicepath_Dump(stdout,orig_path);
//p7_splicepath_DumpScores(stdout,orig_path,graph);
//fflush(stdout);
//if(info->thread_id >= 0) pthread_mutex_unlock(info->mutex);

    
    final_path = p7_splice_FindExons(info, orig_path, path_seq);
 

//printf("FINAL PATH\n");
//p7_splicepath_Dump(stdout,final_path);

//fflush(stdout);
       //if(!path->frameshift)
    
    if(final_path != NULL) {
      
      p7_splice_AlignPath(graph, final_path, orig_path, pli, tophits, om, gm, gcode, path_seq, info->db_nuc_cnt, gm_fs->fs, info, &success);
     
      if(success) {
        /* Break edges that overlap the hit so that paths do not intertwine */
        hit_min = ESL_MIN(pli->hit->dcl->iali, pli->hit->dcl->jali); 
        hit_max = ESL_MAX(pli->hit->dcl->iali, pli->hit->dcl->jali); 
        p7_splice_EnforceRangeBounds(graph, hit_min, hit_max); 
        pli->hit->dcl = NULL;
      
        for(h = 0; h < graph->num_nodes; h++) {
          node_min = ESL_MIN(graph->th->hit[h]->dcl->iali, graph->th->hit[h]->dcl->jali);
          node_max = ESL_MAX(graph->th->hit[h]->dcl->iali, graph->th->hit[h]->dcl->jali);

          if(node_min <= hit_max && node_max >= hit_min)  
            graph->node_in_graph[h] = FALSE;
        }
      }
      else {
        /* Break edges that overlap the path so that paths do not intertwine */ 
        path_min = ESL_MIN(orig_path->iali[0], orig_path->jali[orig_path->path_len-1]);
        path_max = ESL_MAX(orig_path->iali[0], orig_path->jali[orig_path->path_len-1]);
        p7_splice_EnforceRangeBounds(graph, path_min, path_max);

        for(s = 0; s < orig_path->path_len; s++) 
          graph->node_in_graph[orig_path->node_id[s]] = FALSE;
         
      }
    }
   
    esl_sq_Destroy(path_seq);
    p7_splicepath_Destroy(orig_path);
    p7_splicepath_Destroy(final_path);
    p7_splicepipeline_Reuse(pli);
      orig_path = p7_splicepath_GetBestPath_Unspliced(graph);
//p7_splicegraph_DumpHits(stdout, graph); 
//p7_splicegraph_DumpGraph(stdout, graph);

    
  }

printf("\nComplete Query %s Target %s strand %c seqidx %ld\n", gm->name, graph->seqname, (graph->revcomp ? '-' : '+'), graph->seqidx);
  fflush(stdout);  
  return eslOK;

}




/*  Function: p7_splice_ExtendPath
 *  Synopsis: Add seed hits to beginning or end of path 
 *
 *  Purpose : Recover potential exons that were discarded by
 *            the BATH Forward filter but saved in <ifo->sh>
 *            and add them to beginning or end of splice path
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_splice_ExtendPath(P7_TOPHITS *seed_hits, P7_PROFILE *gm, SPLICE_PATH *path, SPLICE_GRAPH *graph)
{
  int i,s,n;
  int gap_len;
  int skip;
  P7_HIT *first_hit;
  P7_HIT *last_hit; 
  P7_HIT *curr_hit;
  SPLICE_GRAPH *tmp_graph;
  SPLICE_PATH  *tmp_path;
  SPLICE_EDGE  *tmp_edge;
  SPLICE_EDGE  *edge;
  
  first_hit = graph->th->hit[path->node_id[0]];
  last_hit  = graph->th->hit[path->node_id[path->path_len-1]];

  /* Only extend from beginning if first hit is an original hit */
  if(path->node_id[0] < graph->anchor_N) {
    tmp_graph = p7_splicegraph_Create();
    tmp_graph->seqidx  = graph->seqidx;
    tmp_graph->revcomp = graph->revcomp;
    tmp_graph->seqname = graph->seqname; 
  
    /*Create a new graph for the first node and the possible upstream extensions */ 
    p7_splicegraph_CreateNodes(tmp_graph, 1);  
    p7_splicegraph_AddNode(tmp_graph, first_hit);
  
    tmp_graph->reportable[0]    = TRUE;
    tmp_graph->orig_hit_idx[0]  = path->node_id[0];
    tmp_graph->split_orig_id[0] = 0;
  
    tmp_graph->anchor_N = 1;
    
    /* Look for seed his to add to the tmp graph */
    for(i = 0; i < seed_hits->N; i++) {
      
      curr_hit = &(seed_hits->unsrt[i]);
      
      /* skip seeds already added to a graph */
      if(curr_hit->dcl->is_included) continue;

      /* Is the seed hit on the same sequence and strand as the tmp_graph */
      if (curr_hit->seqidx != tmp_graph->seqidx) continue;
      if (tmp_graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
      if ((!tmp_graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;
  
      /*If the seed is upstream of the first hit add it to the tmp graph unless there are any anchor nodes between */
      if(p7_splice_HitUpstream(curr_hit->dcl, first_hit->dcl, tmp_graph->revcomp)) { 
        skip = FALSE;
        for(n = 0; n < graph->anchor_N; n++) {
          if(!graph->node_in_graph[n]) continue;
          if(p7_splice_HitUpstream(curr_hit->dcl, graph->th->hit[n]->dcl, tmp_graph->revcomp)) {
            if(p7_splice_HitBetween(curr_hit->dcl, graph->th->hit[n]->dcl, first_hit->dcl, tmp_graph->revcomp)) skip = TRUE;
          }
        }
        if(skip) continue;
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

    p7_splice_CreateExtensionEdges(tmp_graph, gm); 
    tmp_path = p7_splicepath_GetBestPath_Extension(graph,tmp_graph); 

    /* Set the most upstream hit in the tmp_path to start at the minimum 
     * hmm for all hits in tmp_path and add it to the original path  */
    for(s = tmp_path->path_len-2; s >= 0; s--) {

      tmp_edge = p7_splicegraph_GetEdge(tmp_graph, tmp_path->node_id[s], tmp_path->node_id[s+1]);
    
      curr_hit = tmp_graph->th->hit[tmp_path->node_id[s]];     
      curr_hit->dcl->is_included = TRUE;

      p7_splicegraph_AddNode(graph, curr_hit);

      edge = p7_splicegraph_AddEdge(graph, graph->num_nodes-1, path->node_id[0]);
      edge->upstream_amino_end     = tmp_edge->upstream_amino_end;
      edge->downstream_amino_start = tmp_edge->downstream_amino_start;
      edge->upstream_nuc_end       = tmp_edge->upstream_nuc_end;
      edge->downstream_nuc_start   = tmp_edge->downstream_nuc_start;
      edge->splice_score           = tmp_edge->splice_score;
      edge->jump_edge              = FALSE;

      p7_splicepath_Insert(path, 0); 
      path->node_id[0]  = graph->num_nodes-1; 
      path->ihmm[0]     =  tmp_path->ihmm[s];
      path->jhmm[0]     =  tmp_path->jhmm[s];
      path->iali[0]     =  tmp_path->iali[s];
      path->jali[0]     =  tmp_path->jali[s];
      path->aliscore[0] = tmp_path->aliscore[s];
    }

    p7_splicepath_Destroy(tmp_path);
    p7_splicegraph_Destroy(tmp_graph);
  }

  /* Only extend from end if last hit is an original hit */
  if(path->node_id[path->path_len-1] < graph->anchor_N) {
    tmp_graph = p7_splicegraph_Create();
    tmp_graph->seqidx  = graph->seqidx;
    tmp_graph->revcomp = graph->revcomp;
    tmp_graph->seqname = graph->seqname;

    p7_splicegraph_CreateNodes(tmp_graph, 1);
    p7_splicegraph_AddNode(tmp_graph, last_hit);
 
    tmp_graph->reportable[0]    = TRUE;
    tmp_graph->orig_hit_idx[0]  = path->node_id[path->path_len-1];
    tmp_graph->split_orig_id[0] = 0;

    tmp_graph->anchor_N = 1;

    for(i = 0; i < seed_hits->N; i++) {

      curr_hit = &(seed_hits->unsrt[i]);

      /* skip seeds already added to a graph */
      if(curr_hit->dcl->is_included) continue;

      /* Is the seed hit on the same sequence and strand as the tmp_graph */
      if (curr_hit->seqidx != tmp_graph->seqidx) continue;
      if (tmp_graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
      if ((!tmp_graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;

      /*If the seed is upstream of the first hit, add it to the tmp graph unless there are any anchor nodes between*/
      if(p7_splice_HitUpstream(last_hit->dcl, curr_hit->dcl, tmp_graph->revcomp)) {
        skip = FALSE;
        for(n = 0; n < graph->anchor_N; n++) {
          if(!graph->node_in_graph[n]) continue;
          if(p7_splice_HitUpstream(graph->th->hit[n]->dcl, curr_hit->dcl, tmp_graph->revcomp)) {
            if(p7_splice_HitBetween(last_hit->dcl, graph->th->hit[n]->dcl, curr_hit->dcl, tmp_graph->revcomp)) skip = TRUE;
          }
        }
        if(skip) continue;

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

    p7_splice_CreateExtensionEdges(tmp_graph, gm); 
    tmp_path = p7_splicepath_GetBestPath_Extension(graph, tmp_graph);

    /* Set the most downstream hit in the tmp_path to end at the maximum  
     * hmm for all hits in tmp_path and add it to the original path  */
    /* Add any seed hits from tmp_path to real path */
    for(s = 1; s < tmp_path->path_len; s++) { 
    
      tmp_edge = p7_splicegraph_GetEdge(tmp_graph, tmp_path->node_id[s-1], tmp_path->node_id[s]);
 
      curr_hit = tmp_graph->th->hit[tmp_path->node_id[s]];
      curr_hit->dcl->is_included = TRUE;
      p7_splicegraph_AddNode(graph, curr_hit);

      edge = p7_splicegraph_AddEdge(graph, path->node_id[path->path_len-1], graph->num_nodes-1);
      edge->upstream_amino_end     = tmp_edge->upstream_amino_end;
      edge->downstream_amino_start = tmp_edge->downstream_amino_start;
      edge->upstream_nuc_end       = tmp_edge->upstream_nuc_end;
      edge->downstream_nuc_start   = tmp_edge->downstream_nuc_start;
      edge->splice_score           = tmp_edge->splice_score;
      edge->jump_edge              = FALSE;

      p7_splicepath_Insert(path, path->path_len);
      path->node_id[path->path_len-1]  = graph->num_nodes-1;
      path->ihmm[path->path_len-1]     =  tmp_path->ihmm[s];
      path->jhmm[path->path_len-1]     =  tmp_path->jhmm[s];
      path->iali[path->path_len-1]     =  tmp_path->iali[s];
      path->jali[path->path_len-1]     =  tmp_path->jali[s];
      path->aliscore[path->path_len-1] = tmp_path->aliscore[s];
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
p7_splice_CreateUnsplicedEdges(SPLICE_GRAPH *graph, P7_PROFILE *gm) 
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
      if (( graph->revcomp  && th->hit[up]->dcl->iali <= th->hit[down]->dcl->iali) ||
         ((!graph->revcomp) && th->hit[up]->dcl->iali >= th->hit[down]->dcl->iali)) continue;

      if (( graph->revcomp  && th->hit[up]->dcl->jali <= th->hit[down]->dcl->jali) ||
         ((!graph->revcomp) && th->hit[up]->dcl->jali >= th->hit[down]->dcl->jali)) continue;

      if(graph->revcomp)
        seq_gap_len = th->hit[up]->dcl->jali - th->hit[down]->dcl->iali - 1;
      else
        seq_gap_len = th->hit[down]->dcl->iali - th->hit[up]->dcl->jali - 1;        

      if(seq_gap_len > MAX_INTRON_LENG) continue;

      amino_gap_len = th->hit[down]->dcl->ihmm - th->hit[up]->dcl->jhmm - 1;
      if(amino_gap_len > 0 && seq_gap_len < amino_gap_len) continue;

      if(th->hit[up]->dcl->ihmm > th->hit[down]->dcl->jhmm ) { 
        /* If hits are upstream/ downstream on the sequence but are the opposite on th modle create a jump edge to prevent paths from skipping over other nodes */
        if( up < graph->anchor_N && down < graph->anchor_N) {
          edge = p7_splicegraph_AddEdge(graph, up, down);
          edge->splice_score   = (th->hit[up]->dcl->aliscore + th->hit[down]->dcl->aliscore) * -1; 
          // TODO 
          //edge->splice_score   = -eslCONST_LOG2 + p7P_TSC(gm, th->hit[down]->dcl->ihmm-1, p7P_BM);
          edge->jump_edge      = TRUE;
          edge->upstream_amino_end     = th->hit[up]->dcl->jhmm;
          edge->downstream_amino_start = th->hit[down]->dcl->ihmm;
          edge->upstream_nuc_end       = th->hit[up]->dcl->jali;
          edge->downstream_nuc_start   = th->hit[down]->dcl->iali;
        }

      }
      else if(th->hit[up]->dcl->ihmm < th->hit[down]->dcl->ihmm || 
              th->hit[up]->dcl->jhmm < th->hit[down]->dcl->jhmm) { 
              
        edge = p7_splicegraph_AddEdge(graph, up, down);
        edge->jump_edge      = FALSE; 
        /* If hits overlap, find the minimum lost score to remove the overlap */
        
        p7_spliceedge_AliScoreEdge(edge, gm, th->hit[up]->dcl, th->hit[down]->dcl); 
       
        edge->upstream_amino_end     = th->hit[up]->dcl->jhmm;
        edge->downstream_amino_start = th->hit[down]->dcl->ihmm; 
        edge->upstream_nuc_end       = th->hit[up]->dcl->jali;
        edge->downstream_nuc_start   = th->hit[down]->dcl->iali;

        /* If the edge has an hmm overlap and cost of eliminating that overlap is greater than the B->M entry 
         * for the downstream hit then these hits are better off seperate and we will remove the edge */
        if(edge->splice_score < -eslCONST_LOG2 + p7P_TSC(gm, th->hit[down]->dcl->ihmm-1, p7P_BM)) {
          graph->num_edges[up]--;
          graph->tot_edges--;
        }
      }   
    }
  } 

  
  return eslOK;

}

int
p7_splice_CreateExtensionEdges(SPLICE_GRAPH *graph, P7_PROFILE *gm) 
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
      if(seq_gap_len < amino_gap_len) continue;
      edge = p7_splicegraph_AddEdge(graph, up, down);

      /* If hits overlap, find the minimum lost socre to remove the overlap */
      p7_spliceedge_AliScoreEdge(edge, gm, th->hit[up]->dcl, th->hit[down]->dcl); 

      edge->upstream_amino_end     = th->hit[up]->dcl->jhmm;
      edge->downstream_amino_start = th->hit[down]->dcl->ihmm; 
      edge->upstream_nuc_end       = th->hit[up]->dcl->jali;
      edge->downstream_nuc_start   = th->hit[down]->dcl->iali;
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
 *            locate splice boundries.
 *
 * Returns:   <SPLICE_PATH>.
 *
 */
SPLICE_PATH*
p7_splice_FindExons(SPLICE_WORKER_INFO *info, SPLICE_PATH *path, ESL_SQ *path_seq)
{

  int i, s;
  int i_start, i_end;
  int k_start, k_end;
  int next_i_start, next_k_start;
  int contains_anchor;
  int s_start, s_end;
  float           ali_score;
  P7_FS_PROFILE   *gm_fs;
  ESL_GENCODE     *gcode;
  SPLICE_EDGE     *edge;
  SPLICE_PIPELINE *pli;
  SPLICE_PATH *tmp_path;
  SPLICE_PATH *ret_path;
  SPLICE_GRAPH *graph;
  int status;
  
  graph = info->graph;
  pli   = info->pli;
  gcode = info->gcode;
  gm_fs = info->gm_fs;

  ret_path   = NULL;
 
  /* If there is only one hit in the path we must still realign it in case that hit aligned through an intron */
  if(path->path_len == 1) {
   
     i_start = path->iali[0];
     i_end   = path->jali[0];
     k_start = path->ihmm[0];
     k_end   = path->jhmm[0];
     
    /* Covert sequence positions to sub sequence positions */
    if(path->revcomp) {
      i_start = path_seq->n + path_seq->end - i_start;
      i_end   = path_seq->n + path_seq->end - i_end;
    }
    else  {
      i_start = i_start - path_seq->start + 1;
      i_end   = i_end   - path_seq->start + 1;
    }
     
     ret_path = p7_splice_AlignExons(pli, gm_fs, pli->bg, path_seq, gcode, i_start, i_end, k_start, k_end, FALSE, FALSE, &next_i_start, &next_k_start, &ali_score);
    if(ret_path == NULL)     ESL_XEXCEPTION(eslFAIL, "Splice Viterbi Fail");
  
    for(s = 0; s < ret_path->path_len; s++) {
      
      ret_path->node_id[s] = path->node_id[0];

      ret_path->ihmm[s] = ret_path->ihmm[s] + path->ihmm[0] - 1;
      ret_path->jhmm[s] = ret_path->jhmm[s] + path->ihmm[0] - 1;
      if(path->revcomp) {
        ret_path->iali[s] = path_seq->n - ret_path->iali[s] + path_seq->end;
        ret_path->jali[s] = path_seq->n - ret_path->jali[s] + path_seq->end;
      }
      else {
        ret_path->iali[s] = path_seq->start + ret_path->iali[s] - 1;
        ret_path->jali[s] = path_seq->start + ret_path->jali[s] - 1;
      }
    }
    ret_path->revcomp = path->revcomp;
    ret_path->frameshift = path->frameshift; 
    return ret_path;
  }

  for(s_start = 0; s_start < path->path_len; s_start++) if(path->node_id[s_start] < graph->anchor_N) break;
  for(s_end = path->path_len-1; s_end >= 0; s_end--) if(path->node_id[s_end] < graph->anchor_N) break;

  next_i_start = next_k_start = 0; 

  //TODO for extensions swithc to doing pairs and impement aliscore checks and save splice sites
  /* If their are upstream extention nodes see if they are recoved by spliced viterbi */ 
  if(s_start != 0) {
    k_start = path->ihmm[0]; 
    i_start = path->iali[0]; 
    k_end   = path->jhmm[s_start];
    i_end   = path->jali[s_start];
    
    edge =  p7_splicegraph_GetEdge(graph, path->node_id[0], path->node_id[1]);

    /* Covert sequence positions to sub sequence positions */
    if(path->revcomp) {
      i_start = path_seq->n + path_seq->end - i_start;
      i_end   = path_seq->n + path_seq->end - i_end;
    }
    else  {
      i_start = i_start - path_seq->start + 1; 
      i_end   = i_end   - path_seq->start + 1;
    }

    if(k_end <= k_start || i_end <= i_start + 2) {
      edge->splice_score = -eslINFINITY;
      return NULL;
    } 

    ali_score = -eslINFINITY;
    tmp_path = p7_splice_AlignExons(pli, gm_fs, pli->bg, path_seq, gcode, i_start, i_end, k_start, k_end, TRUE, FALSE, &next_i_start, &next_k_start, &ali_score);

    if(tmp_path == NULL) {
      edge->splice_score = -eslINFINITY;
      p7_splicepath_Destroy(tmp_path);
      p7_splicepath_Destroy(ret_path);
      return NULL;
    }

    for(s = 0; s < tmp_path->path_len; s++) {

      tmp_path->iali[s] = path->revcomp ? path_seq->n - tmp_path->iali[s] + path_seq->end : 
                                          path_seq->start + tmp_path->iali[s] - 1;
      tmp_path->jali[s] = path->revcomp ? path_seq->n - tmp_path->jali[s] + path_seq->end : 
                                          path_seq->start + tmp_path->jali[s] - 1;

    }

    ret_path = p7_splicepath_Clone(tmp_path);

    next_i_start = path->revcomp ? path_seq->n - next_i_start + path_seq->end :
                                   path_seq->start + next_i_start - 1;    

    p7_gmx_Reuse(pli->vit);
    p7_splicepath_Destroy(tmp_path);
  }


  for(s = s_start+1; s <= s_end; s++) {

    k_start = next_k_start == 0 ? path->ihmm[s-1] : next_k_start; 
    i_start = next_i_start == 0 ? path->iali[s-1] : next_i_start;
    k_end   = path->jhmm[s];
    i_end   = path->jali[s];
    
    /* Check if we have spliced this pair of nodes before */
    edge =  p7_splicegraph_GetEdge(graph, path->node_id[s-1], path->node_id[s]);
    if(i_start == edge->i_start && k_start == edge->k_start) {

      if(ret_path == NULL) {
        ret_path = p7_splicepath_Create(2);
        ret_path->iali[0] = i_start;
        ret_path->ihmm[0] = k_start;
      } 
      else p7_splicepath_Insert(ret_path, ret_path->path_len);

      ret_path->jali[ret_path->path_len-2] = edge->upstream_nuc_end;
      ret_path->jhmm[ret_path->path_len-2] = edge->upstream_amino_end;
      ret_path->iali[ret_path->path_len-1] = edge->downstream_nuc_start;
      ret_path->ihmm[ret_path->path_len-1] = edge->downstream_amino_start;
      ret_path->jali[ret_path->path_len-1] = i_end;
      ret_path->jhmm[ret_path->path_len-1] = k_end;

      ret_path->node_id[ret_path->path_len-2] = 0;
      ret_path->node_id[ret_path->path_len-1] = 0;
    
      next_k_start = edge->next_k_start; 
      next_i_start = edge->next_i_start;

//printf("Redo \n");
//p7_splicepath_Dump(stdout,ret_path);

      continue;
    }

    /* Covert sequence positions to sub sequence positions */
    if(path->revcomp) {
      i_start = path_seq->n + path_seq->end - i_start;
      i_end   = path_seq->n + path_seq->end - i_end;
    }
    else  {
      i_start = i_start - path_seq->start + 1;
      i_end   = i_end   - path_seq->start + 1;
    }
    
    if(k_end <= k_start || i_end <= i_start + 2) {
      //printf("gm_fs->name %s path_seq->name %s revcomp %d\n", gm_fs->name, path_seq->name, path->revcomp);
      edge->splice_score = -eslINFINITY;
      p7_splicepath_Destroy(ret_path);
      return NULL;
    }
 
    ali_score = -eslINFINITY; 
    tmp_path = p7_splice_AlignExons(pli, gm_fs, pli->bg, path_seq, gcode, i_start, i_end, k_start, k_end, FALSE, FALSE, &next_i_start, &next_k_start, &ali_score);
     
    /* If the spliced aligment score is less than the score of the two hits plus the B->M penalty,       
     * these exons are better off as seperate hits so we erase the edge and return NULL */ 

    if(tmp_path == NULL || ali_score < path->aliscore[s-1] + path->aliscore[s] + -eslCONST_LOG2 + p7P_TSC(gm_fs, path->ihmm[s]-1, p7P_BM)) {
      
      edge->splice_score = -eslINFINITY;

      p7_splicepath_Destroy(tmp_path);
      p7_splicepath_Destroy(ret_path); 
      return NULL;
    }

    /* Converpt path coords back to full sequence */
    for(i = 0; i < tmp_path->path_len; i++) {
      tmp_path->iali[i] = path->revcomp ? path_seq->n - tmp_path->iali[i] + path_seq->end :
                                          path_seq->start + tmp_path->iali[i] - 1;
      tmp_path->jali[i] = path->revcomp ? path_seq->n - tmp_path->jali[i] + path_seq->end :
                                          path_seq->start + tmp_path->jali[i] - 1;
    }
 
    next_i_start = path->revcomp ? path_seq->n - next_i_start + path_seq->end :
                                   path_seq->start + next_i_start - 1;

 
    /* If exactly two exons wer found record the splice site on the edge */ 
    if(tmp_path->path_len == 2) {
      edge->upstream_nuc_end       = tmp_path->jali[0];
      edge->upstream_amino_end     = tmp_path->jhmm[0];
      edge->downstream_nuc_start   = tmp_path->iali[1];   
      edge->downstream_amino_start = tmp_path->ihmm[1]; 
      edge->i_start = i_start;
      edge->k_start = k_start;
      edge->next_i_start = next_i_start;
      edge->next_k_start = next_k_start;
    }
    else {
      edge->i_start = -1;
      edge->k_start = -1;
      edge->next_i_start = -1;
      edge->next_k_start = -1;
    }
    //TODO if more than 2 were found add the new exons and edges
    
    if(ret_path == NULL) ret_path = p7_splicepath_Clone(tmp_path);
    else {
      ret_path->jali[ret_path->path_len-1] = tmp_path->jali[0];
      ret_path->jhmm[ret_path->path_len-1] = tmp_path->jhmm[0];

      for(i = 1; i < tmp_path->path_len; i++) {   
        p7_splicepath_Insert(ret_path, ret_path->path_len);
        ret_path->iali[ret_path->path_len-1] = tmp_path->iali[i];
        ret_path->jali[ret_path->path_len-1] = tmp_path->jali[i];
        ret_path->ihmm[ret_path->path_len-1] = tmp_path->ihmm[i];
        ret_path->jhmm[ret_path->path_len-1] = tmp_path->jhmm[i]; 
        
      }
    } 

//printf("New \n");
//p7_splicepath_Dump(stdout,ret_path);
  
    //p7_splicepath_Dump(stdout,tmp_path);
    p7_gmx_Reuse(pli->vit);
    p7_splicepath_Destroy(tmp_path); 
  }

  if(s_end != path->path_len-1) {
    k_start = next_k_start == 0 ? path->ihmm[s_end] : next_k_start;
    i_start = next_i_start == 0 ? path->iali[s_end] : next_i_start;
    k_end   = path->jhmm[path->path_len-1];
    i_end   = path->jali[path->path_len-1];

    edge =  p7_splicegraph_GetEdge(graph, path->node_id[path->path_len-2], path->node_id[path->path_len-1]);

    /* Covert sequence positions to sub sequence positions */
    if(path->revcomp) {
      i_start = path_seq->n + path_seq->end - i_start;
      i_end   = path_seq->n + path_seq->end - i_end;
    }
    else  {
      i_start = i_start - path_seq->start + 1;
      i_end   = i_end   - path_seq->start + 1;
    }

    if(k_end <= k_start || i_end <= i_start + 2) {
      edge->splice_score = -eslINFINITY;
      p7_splicepath_Destroy(ret_path);
      return NULL;
    }

    ali_score = -eslINFINITY;
    tmp_path = p7_splice_AlignExons(pli, gm_fs, pli->bg, path_seq, gcode, i_start, i_end, k_start, k_end, FALSE, TRUE, &next_i_start, &next_k_start, &ali_score);
    

    if(tmp_path == NULL) {
      edge->splice_score = -eslINFINITY;
      p7_splicepath_Destroy(tmp_path);
      p7_splicepath_Destroy(ret_path);
      return NULL;
    }

    for(s = 0; s < tmp_path->path_len; s++) {
      
      tmp_path->iali[s] = path->revcomp ? path_seq->n - tmp_path->iali[s] + path_seq->end :
                                          path_seq->start + tmp_path->iali[s] - 1;
      tmp_path->jali[s] = path->revcomp ? path_seq->n - tmp_path->jali[s] + path_seq->end :
                                          path_seq->start + tmp_path->jali[s] - 1;
    }

    if(ret_path == NULL) ret_path = p7_splicepath_Clone(tmp_path);
    else {
      ret_path->jali[ret_path->path_len-1] = tmp_path->jali[0];
      ret_path->jhmm[ret_path->path_len-1] = tmp_path->jhmm[0];

      for(i = 1; i < tmp_path->path_len; i++) {
        p7_splicepath_Insert(ret_path, ret_path->path_len);
        ret_path->iali[ret_path->path_len-1] = tmp_path->iali[i];
        ret_path->jali[ret_path->path_len-1] = tmp_path->jali[i];
        ret_path->ihmm[ret_path->path_len-1] = tmp_path->ihmm[i];
        ret_path->jhmm[ret_path->path_len-1] = tmp_path->jhmm[i];

      }
    }

  
    p7_gmx_Reuse(pli->vit);
    p7_splicepath_Destroy(tmp_path);
    
  }
 
  if(ret_path != NULL) {
    ret_path->revcomp = path->revcomp;

    for(s = 0; s < path->path_len; s++) {
      if(path->node_id[s] >= graph->anchor_N) continue;
      for(i = 0; i < ret_path->path_len; i++) {
        if(ret_path->node_id[i] >= 0) continue;
        if(seq_overlap(ret_path, i, path, s) > 0.0 && hmm_overlap(ret_path, i, path, s) > 0.0)
          ret_path->node_id[i] = path->node_id[s];
      }
    }

    contains_anchor = FALSE;
    for(i = 0; i < ret_path->path_len; i++) {
      if(ret_path->node_id[i] >= 0) contains_anchor = TRUE;
    }
    if(!contains_anchor) {
      p7_splicepath_Destroy(ret_path);
      ret_path = NULL;
    }
  }

  return ret_path;

  ERROR:
    return NULL;
}


SPLICE_PATH*
p7_splice_AlignExons(SPLICE_PIPELINE *pli, P7_FS_PROFILE *gm_fs, P7_BG *bg, ESL_SQ *path_seq, const ESL_GENCODE *gcode, int i_start, int i_end, int k_start, int k_end, int extend_up, int extend_down, int *next_i_start, int *next_k_start, float *ali_score)
{
 
  int         L = i_end - i_start + 1;
  int         M = k_end - k_start + 1;
  int         y, z;
  int         z1, z2;
  int         intron_cnt;
  int         step_cnt;
  int         start_new;
  P7_GMX       *gx = pli->vit;
  P7_TRACE     *tr;
  SPLICE_PATH   *ret_path;

//TODO dont need Lx for p7_gmx_sp_GrowTo
//  char strand;
//  strand = path_seq->start < path_seq->end ? '+' : '-';
  //printf("seq %s strand %c M %d L %d\n", path_seq->name, strand, M, L);
  p7_gmx_sp_GrowTo(pli->vit, M, L, L);
  p7_splicepipline_GrowIndex(pli->sig_idx, M, L, ALIGNMENT_MODE);
  p7_fs_ReconfigLength(gm_fs, L);
  
  if(extend_up) 
    p7_spliceviterbi_translated_semiglobal_extendup(pli, path_seq->dsq, gcode, gm_fs, pli->vit, i_start, i_end, k_start, k_end);
  else if(extend_down) 
    p7_spliceviterbi_translated_semiglobal_extenddown(pli, path_seq->dsq, gcode, gm_fs, pli->vit, i_start, i_end, k_start, k_end);
  else {
    p7_spliceviterbi_translated_semiglobal(pli, path_seq->dsq, gcode, gm_fs, pli->vit, i_start, i_end, k_start, k_end);
    /* If the hits were in different frames and no splice site was able to pull score 
     * from the upstream frame to the downstream frame the spliceing is a failure */
    if(gx->xmx[L*p7G_NXCELLS+p7G_C] == -eslINFINITY) return NULL; 
  }

  tr = p7_trace_fs_Create();
  p7_splicevitebi_translated_semiglobal_trace(pli, path_seq->dsq, gcode, gm_fs, pli->vit, tr, i_start, i_end, k_start, k_end, ali_score);
  //p7_trace_fs_Dump(stdout, tr, NULL, NULL, NULL);

  /* Find number of introns in trace */
  intron_cnt = 0;
  for(z = 0; z < tr->N; z++)
    if(tr->st[z] == p7T_P) intron_cnt++;

  /* Find first M state - start of first hit */
  for(z1 = 0; z1 < tr->N; z1++) if(tr->st[z1] == p7T_M) break;
  if (z1 == tr->N) {
    p7_trace_Destroy(tr);
    return NULL;
  }

  /* Find last M state state - end of last hit */
  for(z2 = tr->N-1; z2 >= 0; z2--) if(tr->st[z2] == p7T_M) break;
  if (z2 == -1) {
    p7_trace_Destroy(tr);
    return NULL;
  }

  ret_path = p7_splicepath_Create(intron_cnt+1);

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
  
      
      /* If this is the first step in path the i coords are set my the first M state at
       * trace postion y, otherwise they are set by the P state at trace postion y-1 */
      
      ret_path->node_id[step_cnt] = -1;
      if(step_cnt == 0) {
        ret_path->iali[step_cnt] = tr->i[y] - tr->c[y] + 1; 
        ret_path->ihmm[step_cnt] = tr->k[y];        
        *next_i_start = ret_path->iali[step_cnt]; 
        *next_k_start = ret_path->ihmm[step_cnt];
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
        *next_i_start = tr->i[y] - tr->c[y] + 1;
        *next_k_start = tr->k[y];
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
          ret_path->jali[step_cnt]   = tr->i[z] + 1;
        else if( tr->c[z+1] == 2)
          ret_path->jali[step_cnt]   = tr->i[z] + 2;
      }

      step_cnt++;

      start_new = FALSE;
    }
    
    z++;
    if(tr->st[z] == p7T_M) start_new = TRUE;

  }

  p7_trace_fs_Destroy(tr);

  return ret_path;

}



float
hmm_overlap(SPLICE_PATH *path1, int step1, SPLICE_PATH *path2, int step2)
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
seq_overlap(SPLICE_PATH *path1, int step1, SPLICE_PATH *path2, int step2) 
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


int
p7_splice_HitUpstream(P7_DOMAIN *upstream, P7_DOMAIN *downstream, int revcomp) {

  if(upstream->ihmm > downstream->ihmm || upstream->jhmm > downstream->jhmm)
    return FALSE;
    
  if (( revcomp  && upstream->iali <= downstream->iali) ||
     ((!revcomp) && upstream->iali >= downstream->iali))
    return FALSE;

  if (( revcomp  && upstream->jali <= downstream->jali) ||
     ((!revcomp) && upstream->jali >= downstream->jali))
    return FALSE;

  return TRUE;
  
}

int
p7_splice_HitBetween(P7_DOMAIN *up, P7_DOMAIN *mid, P7_DOMAIN *down, int revcomp) {

  if (( revcomp  && up->iali <= mid->iali) ||
     ((!revcomp) && up->iali >= mid->iali))
    return FALSE;

  if (( revcomp  && mid->jali <= down->jali) ||
     ((!revcomp) && mid->jali >= down->jali))
    return FALSE;

  return TRUE;
  
}




int
p7_splice_AlignPath(SPLICE_GRAPH *graph, SPLICE_PATH *path, SPLICE_PATH *orig_path, SPLICE_PIPELINE *pli, P7_TOPHITS *tophits, P7_OPROFILE *om, P7_PROFILE *gm, ESL_GENCODE *gcode, ESL_SQ *path_seq, int64_t db_nuc_cnt, float fs_prob, SPLICE_WORKER_INFO *info, int *success)
{

  int          i;
  int          s1, s2;
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
if(info->thread_id >= 0) pthread_mutex_lock(info->mutex);
p7_splicepath_Dump(stdout,path);
    ESL_XEXCEPTION(eslFAIL, "NOT mod 3");
if(info->thread_id >= 0) pthread_mutex_unlock(info->mutex); 
    return eslOK;
  }
//printf("start %d end %d n %d\n", path_seq->start, path_seq->end, path_seq->n);
  /* Working backward from the start of the path find the first instance of a stop codon in the extended upstream region of path_seq */
  //printf("path->revcomp %d\n", path->revcomp);
  if (path->revcomp) {
    path_start_pos = path_seq->n - path->iali[0] + path_seq->end;
    ext_start_pos  = path_seq->n - (path->iali[0] + ALIGNMENT_EXT) + path_seq->end;
    for(pos = path->iali[0]+3; pos <= path->iali[0] + ALIGNMENT_EXT; pos+=3) {

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
    ext_start_pos  = (path->iali[0] - ALIGNMENT_EXT) - path_seq->start + 1;
    //printf("ext_start_pos %d\n", ext_start_pos);
    for(pos = path->iali[0]-3; pos >= path->iali[0] - ALIGNMENT_EXT; pos-=3) {
       
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
    ext_end_pos  =  path_seq->n - (path->jali[path->path_len-1] - ALIGNMENT_EXT) + path_seq->end; 
    for(pos = path->jali[path->path_len-1]-1; pos >= path->jali[path->path_len-1] - ALIGNMENT_EXT; pos-=3) {
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
    ext_end_pos  = (path->jali[path->path_len-1] +  ALIGNMENT_EXT) - path_seq->start + 1;
    //printf("ext_end_pos %d %d\n", ext_end_pos, path->jali[path->path_len-1] +  ALIGNMENT_EXT);
    for(pos = path->jali[path->path_len-1]+1; pos <= path->jali[path->path_len-1] +  ALIGNMENT_EXT; pos+=3) {

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
    //printf("ext_end_pos %d %d\n", ext_end_pos, path_seq->start + ext_end_pos - 1 );
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
  //printf("seq_pos %d true %d\n", seq_pos, path_seq->start + seq_pos-1);
  /*Add downstream extension nucleotides */
  for(seq_pos = path_end_pos+1; seq_pos <= ext_end_pos; seq_pos++) {
    
    nuc_index[seq_idx] = seq_pos;
    nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
    seq_idx++;
  } 
//printf("seq_pos %d true %d\n", seq_pos, path_seq->start + seq_pos-1);
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
    printf("seq %s\n", path_seq->name); 
    fflush(stdout);
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
//printf("%10f \n", exp(dom_lnP));
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
        if(path->node_id[exon] < graph->anchor_N && path->node_id[exon] > -1) contains_orig = TRUE;
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
    while( path->node_id[i] >= graph->anchor_N || path->node_id[i] == -1) {
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
      if(remove_node < 0 || remove_node >= graph->anchor_N) {
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
    *success = TRUE; 

#ifdef HMMER_THREADS
  if(info->thread_id >= 0) pthread_mutex_unlock(info->mutex);
#endif /*HMMER_THREADS*/

  }

  /* Check of any nodes that got merged and therefore still need their original hits set to unreported */
  if(success) {
    for(s1 = 0; s1 < orig_path->path_len; s1++) {
      if(orig_path->node_id[s1] >= graph->anchor_N) continue;
      if(graph->th->hit[orig_path->node_id[s1]]->dcl->ad != NULL) continue;
      if(graph->th->hit[orig_path->node_id[s1]]->flags & p7_IS_REPORTED ) {
        for(s2 = 0; s2 < path->path_len; s2++) {
          if(path->node_id[s2] >= graph->anchor_N) continue;
          if(seq_overlap(orig_path, s1, path, s2) > 0.0 && 
             hmm_overlap(orig_path, s1, path, s2) > 0.0) {
#ifdef HMMER_THREADS
if(info->thread_id >= 0) pthread_mutex_lock(info->mutex);
#endif /*HMMER_THREADS*/
            tophits->nreported--;
            graph->th->hit[orig_path->node_id[s1]]->flags &= ~p7_IS_REPORTED;   
            graph->th->hit[orig_path->node_id[s1]]->dcl->is_reported = FALSE;
            if(graph->th->hit[orig_path->node_id[s1]]->flags & p7_IS_INCLUDED) {
              tophits->nincluded--;
              graph->th->hit[orig_path->node_id[s1]]->flags &= ~p7_IS_INCLUDED;
              graph->th->hit[orig_path->node_id[s1]]->dcl->is_included = FALSE;
           }
#ifdef HMMER_THREADS
if(info->thread_id >= 0) pthread_mutex_unlock(info->mutex);
#endif /*HMMER_THREADS*/
          }      
        }
      }
    }
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
  float     bwdsc;
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
  p7_Backward(pli->amino_sq->dsq, pli->amino_sq->n, om, pli->fwd, pli->bwd, &bwdsc);

  if((status = p7_Decoding(om, pli->fwd, pli->bwd, pli->bwd)) == eslERANGE) { 

    p7_gmx_fs_GrowTo(pli->gfwd, gm->M, pli->amino_sq->n, pli->amino_sq->n, p7P_CODONS);
    p7_ReconfigUnihit(gm, pli->amino_sq->n); 

    p7_GForward(pli->amino_sq->dsq, pli->amino_sq->n, gm, pli->gfwd, &envsc);

    p7_gmx_Reuse(pli->gfwd); 
    p7_GViterbi(pli->amino_sq->dsq, pli->amino_sq->n, gm, pli->gfwd, NULL);
    p7_GTrace(pli->amino_sq->dsq, pli->amino_sq->n, gm, pli->gfwd, tr);
        
  }
  else {
    p7_OptimalAccuracy(om, pli->bwd, pli->fwd, &oasc);
    p7_OATrace        (om, pli->bwd, pli->fwd, tr);
    
  }

  seq_score = (envsc-filtersc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);
  if (P > pli->F3) {
    p7_hit_Destroy(hit);
    p7_trace_Destroy(tr);
    return eslOK;
  }

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

}
/*
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

  // Find out how may exteded nucletides are present in path_seq 
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

  // Add upstream extension nucleotides 
  for(seq_pos = 1; seq_pos < path_start_pos; seq_pos++) {
    nuc_index[seq_idx] = seq_pos;
    nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
    seq_idx++;
  } 

  // Copy spliced nucleotides into single sequence and track their original indicies 
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

  // Add downstream extension nucleotides 
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

  // Algin the spliced Amino sequence 
  status = align_spliced_path_frameshift(pli, gm_fs, path_seq, gcode);

  // Alignment failed 
  if(pli->hit == NULL ) {
    if(nuc_dsq   != NULL) free(nuc_dsq);
     
     return eslOK;
  }

  // adjust all coords in hit and path 
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

  // Adjust spliced hit score from nuc_len to gm_fs->max_length*3 
  dom_score  = pli->hit->dcl->envsc;
  
  dom_score -= 2 * log(2. / ((pli->nuc_sq->n/3.)+2));
  dom_score += 2 * log(2. / (gm_fs->max_length+2));
  dom_score -= (pli->nuc_sq->n-ali_len)      * log((float) (pli->nuc_sq->n/3.) / (float) ((pli->nuc_sq->n/3.)+2));
  dom_score += (ESL_MAX(pli->nuc_sq->n, gm_fs->max_length*3)-ali_len) * log((float) gm_fs->max_length / (float) (gm_fs->max_length+2));

  // Bias calculation and adjustments 
  if(pli->do_null2)
    dom_bias = p7_FLogsum(0.0, log(pli->bg->omega) + pli->hit->dcl->domcorrection);
  else
    dom_bias = 0.;

  p7_bg_SetLength(pli->bg, gm_fs->max_length*3);
  p7_bg_NullOne  (pli->bg, pli->nuc_sq->dsq, gm_fs->max_length*3, &nullsc);
  dom_score = (dom_score - (nullsc + dom_bias))  / eslCONST_LOG2;
  
  // Add splice signal penalties 
  pli->hit->dcl->jali = pli->hit->dcl->ad->sqto;
  for(i = 0; i < path->path_len; i++)
    dom_score += path->signal_scores[i];
  
  dom_lnP   = esl_exp_logsurv(dom_score, gm_fs->evparam[p7_FTAUFS], gm_fs->evparam[p7_FLAMBDA]);

  // E-value adjusment 
  dom_lnP += log((float)db_nuc_cnt / (float)gm_fs->max_length);
    
   if ((pli->by_E && exp(dom_lnP) <= pli->E) || ((!pli->by_E) && dom_score >= pli->T)) {


    if ( path->path_len > pli->hit->dcl->ad->exon_cnt) {
      // Shift the path to start at the first hit that was inculded in the alignment and end at the last hit that was included in the alignment 
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
      // Ensure path will still contains an original hit after shifting 
      contains_orig = FALSE;
      for(exon = shift; exon < pli->hit->dcl->ad->exon_cnt ; exon ++ ) {
        if(graph->split_orig_id[path->node_id[exon]] >= 0) contains_orig = TRUE;
      }

      if(!contains_orig) {
        if(nuc_dsq   != NULL) free(nuc_dsq);
        
        return eslOK;
      }

      // Shift path to start at frist hit that is in alignment 
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
  
    // Find first original hit to copy info  
    i = 0;
    while( graph->split_orig_id[path->node_id[i]] < 0) {

      pli->hit->dcl->ad->exon_orig[i] = FALSE;
      pli->hit->dcl->ad->exon_split[i] = path->split[i];
      i++;
    }

    pli->hit->dcl->ad->exon_orig[i] = TRUE;
    pli->hit->dcl->ad->exon_split[i] = path->split[i];

//TODO add back commenting 
#ifdef HMMER_THREADS
  if(info->thread_id >= 0) pthread_mutex_lock(info->mutex);
#endif HMMER_THREADS

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

    // Set all original hits in alignment to unreportable 
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

/TODO add back commenting
#ifdef HMMER_THREADS
  if(info->thread_id >= 0) pthread_mutex_unlock(info->mutex);
#endif HMMER_THREADS
  }
  
  if(nuc_dsq   != NULL) free(nuc_dsq);
  
  return eslOK;

  ERROR:
    if(nuc_index != NULL) free(nuc_index);
    if(nuc_dsq   != NULL) free(nuc_dsq);
    return status;

}

*/
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

    sub_path_len += llabs(path->jali[s] - path->iali[s])+1;

    if(intron_include) {
      intron_len = llabs(path->iali[s+1] - path->jali[s]) - 1;
      amino_gap =  ESL_MAX(0, path->ihmm[s+1] - path->jhmm[s] - 1);

      /* intron length */
      if(intron_len <= MIN_INTRON_RMV + (amino_gap*3)) //short intron - add full length
        sub_path_len += intron_len;
      else {
        sub_path_len += MAX_INTRON_INC + (amino_gap*3);
      }
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
    if(intron_len <= MIN_INTRON_RMV + (amino_gap*3)) {
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
      /* Add first MAX_INTRON_INC/2 nucleotides of seqeunce gap */
      for(x = 0; x < MAX_INTRON_INC/2; x++) {
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

      /* Add degenerate nucleotides to cover amino*/ 
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
        while(true_idx > path->iali[s+1] + MAX_INTRON_INC/2) {
          seq_pos++;
          true_idx--;
        }
      }
      else {
        while(true_idx < path->iali[s+1] - MAX_INTRON_INC/2) {
          seq_pos++;
          true_idx++;
        }
      }
      r_idx[r*2+1] = seq_idx;

      /* Add last MAX_INTRON_INC/2 nucleotides of seqeunce gap */
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

ESL_SQ* 
p7_splice_GetSplicedSequence2(SPLICE_PATH *path, ESL_SQ *path_seq, int *exon_idx, int exon_cnt, int *remove_idx, int *remove_cnt, int64_t **nuc_index)
{

  int i;
  int max_seq_len;
  int seq_idx, seq_pos, true_idx;
  int r_cnt;
  int64_t *n_idx;
  ESL_DSQ *splice_dsq;
  ESL_SQ  *splice_seq;
  int status;

  splice_dsq = NULL;
  splice_seq = NULL;
 
  n_idx = *nuc_index;

  /* Get the max lenght of the cub sequence - all exon regions plus MIN_INTRON_RMV for each inrtron region */
  max_seq_len = exon_idx[(exon_cnt-1)*2+1] - exon_idx[(exon_cnt-1)*2] + 1; 
  
  for(i = exon_cnt-2; i >= 0; i--) {
    max_seq_len += MIN_INTRON_RMV;
    max_seq_len += exon_idx[i*2+1] - exon_idx[i*2] + 1;
  }

  ESL_ALLOC(n_idx,  sizeof(int64_t) * (max_seq_len+2));
  ESL_ALLOC(splice_dsq, sizeof(ESL_DSQ) * (max_seq_len+2));

  //printf("max_seq_len %d\n", max_seq_len);  
  n_idx[0]   = -1;
  splice_dsq[0]  = eslDSQ_SENTINEL;
  seq_idx   = 1;

  seq_pos  = exon_idx[(exon_cnt-1)*2];
  if(path->revcomp) true_idx = path_seq->n - exon_idx[(exon_cnt-1)*2] + path_seq->end;
  else true_idx = path_seq->start + exon_idx[(exon_cnt-1)*2] - 1;
  
  r_cnt = 0;
  for(i = exon_cnt-1; i > 0; i--) {

    /* add the exon region */
    while(seq_pos <= exon_idx[i*2+1]) {
      //printf("exon %d seq_idx %d seq_pos %d true_idx %d nuc %d\n", exon_cnt - i, seq_idx, seq_pos, true_idx, path_seq->dsq[seq_pos]);
      n_idx[seq_idx]      = true_idx;
      splice_dsq[seq_idx] = path_seq->dsq[seq_pos];
      seq_idx++;
      seq_pos++;
      if(path->revcomp) true_idx--;
      else              true_idx++;
    }   

    //printf("exon %d seq_idx %d\n", i+1, seq_idx);
   // printf("exon_idx[i*2+1] %d exon_idx[(i-1)*2] %d\n", exon_idx[i*2+1], exon_idx[(i-1)*2]);
    /* If the intron region is less than of equal in length MIN_INTRON_RMV, include full intron region */
    if(exon_idx[(i-1)*2] - exon_idx[i*2+1] - 1 <= MIN_INTRON_RMV)
    {
      while(seq_pos < exon_idx[(i-1)*2]) {
        //printf("exon %d seq_idx %d seq_pos %d true_idx %d nuc %d\n", exon_cnt - i, seq_idx, seq_pos, true_idx, path_seq->dsq[seq_pos]);
        n_idx[seq_idx]      = true_idx;
        splice_dsq[seq_idx] = path_seq->dsq[seq_pos];
        seq_idx++;
        seq_pos++;
        if(path->revcomp) true_idx--;
        else              true_idx++;
      }
    }
    /* If the intron region is longer than MIN_INTRON_RMV, include (MAX_INTRON_INC/2) from the start 
     * and (MAX_INTRON_INC/2) from the end aand record the location of the removed intron region */
    else {
      /* include start */
      while(seq_pos <= exon_idx[i*2+1] + (MAX_INTRON_INC/2)) {
//        printf("exon %d seq_idx %d seq_pos %d true_idx %d nuc %d\n", exon_cnt - i, seq_idx, seq_pos, true_idx, path_seq->dsq[seq_pos]);
        n_idx[seq_idx]      = true_idx;
        splice_dsq[seq_idx] = path_seq->dsq[seq_pos];
        seq_idx++;
        seq_pos++;
        if(path->revcomp) true_idx--;
        else              true_idx++;   
      }
      remove_idx[r_cnt] = seq_idx;  
      r_cnt++;

      /* exclude middle */
      while(seq_pos < exon_idx[(i-1)*2] - (MAX_INTRON_INC/2)) {
        seq_pos++;
        if(path->revcomp) true_idx--;        
        else              true_idx++;
      }
      
      /* include end */
      while(seq_pos < exon_idx[(i-1)*2]) {
    //    printf("exon %d seq_idx %d seq_pos %d true_idx %d nuc %d\n", exon_cnt - (i-1), seq_idx, seq_pos, true_idx, path_seq->dsq[seq_pos]);
        n_idx[seq_idx]      = true_idx;
        splice_dsq[seq_idx] = path_seq->dsq[seq_pos];
        seq_idx++;
        seq_pos++;
        if(path->revcomp) true_idx--;
        else              true_idx++;
      }
    }
  //  printf("intron %d seq_idx %d\n", i+1, seq_idx);
  } 

  /* add the exon region */
  while(seq_pos <= exon_idx[1]) {
   // printf("exon %d seq_idx %d seq_pos %d true_idx %d nuc %d\n", exon_cnt, seq_idx, seq_pos, true_idx, path_seq->dsq[seq_pos]);
    n_idx[seq_idx]      = true_idx;
    splice_dsq[seq_idx] = path_seq->dsq[seq_pos];
    seq_idx++;
    seq_pos++;
    if(path->revcomp) true_idx--;
    else              true_idx++;
  }

  n_idx[seq_idx] = -1;
  splice_dsq[seq_idx]   = eslDSQ_SENTINEL;

  splice_seq   = esl_sq_CreateDigitalFrom(path_seq->abc, NULL, splice_dsq, seq_idx-1, NULL,NULL,NULL);
  if(path->revcomp) {
    splice_seq->start = path_seq->n - exon_idx[(exon_cnt-1)*2] + path_seq->end;
    splice_seq->end   = path_seq->n - exon_idx[1]              + path_seq->end;
  }
  else {              
    splice_seq->start = path_seq->start + exon_idx[(exon_cnt-1)*2] - 1;
    splice_seq->end   = path_seq->start + exon_idx[1]              - 1; 
  }

  *nuc_index  = n_idx;
  *remove_cnt = r_cnt;
 
  if(splice_dsq != NULL) free(splice_dsq);

  return splice_seq;

  ERROR:
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
  int i, k, n, z;
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
  z = z1;
  while (z < z2) {
    if (tr->st[z] == p7T_M) {
     
      dom->scores_per_pos[n] = p7P_MSC(gm, k, amino_dsq[i]) + log(1. - fs_prob*4);
      if (tr->st[z-1] == p7T_I) 
        dom->scores_per_pos[n] += p7P_TSC(gm, k-1, p7P_IM);
      else if (tr->st[z-1] == p7T_D)
        dom->scores_per_pos[n] += p7P_TSC(gm, k-1, p7P_DM); 

      i++; k++; z++; n++;
 
      while(z < z2 && tr->st[z] == p7T_M) {
        dom->scores_per_pos[n] =  p7P_MSC(gm, k,   amino_dsq[i]) + log(1. - fs_prob*4);
        dom->scores_per_pos[n] += p7P_TSC(gm, k-1, p7P_MM);

        i++; k++; z++; n++;
      }

    }
    else if (tr->st[z] == p7T_I) {
      dom->scores_per_pos[n] = p7P_TSC(gm, k, p7P_MI);
      i++; z++; n++;
      while (z < z2 && tr->st[z] == p7T_I) {
        dom->scores_per_pos[n] = p7P_TSC(gm, k, p7P_II);
        i++; z++; n++;
      }
     }
     else if (tr->st[z] == p7T_D) {
       dom->scores_per_pos[n] = p7P_TSC(gm, k-1, p7P_MD);
       k++; z++; n++;
       while (z < z2 && tr->st[z] == p7T_D)  {
         dom->scores_per_pos[n] = p7P_TSC(gm, k-1, p7P_DD);
         k++; z++; n++;
       } 
     }
     else ESL_XEXCEPTION(eslFAIL, "Impossible state from p7_splice_ComputeAliScores()");
  }

  dom->aliscore = dom->scores_per_pos[0];
  for (n=1; n<N; n++)
    dom->aliscore += dom->scores_per_pos[n];
 
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


