/* Create spliced alignments
 *   1. Split hits from tophits by target sequence and strand
 *   2. Add each set of hits to a graph as anchor nodes
 *   3. Add additional hits from saved_hits as seed nodes
 *   4. Add edges between nodes that are upstream/downstream compatible
 *   5. Find the best path that contains at least one anchor node
 *   6. Find splice sites with spliced Viterbi
 *   7. Realign spliced exons to model with Fwd/Bwd
 *   8. Return to step 5 until no more paths with anchor nodes exist.
 *
 * Contents:
 *    1. Graph Creation
 *    2. Path Finding and Splicing
 *    3. Path Alignenment 
 *    3. Helper Functions
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


/*****************************************************************
 * 1. Graph Creation
 *****************************************************************/

static int serial_loop(SPLICE_WORKER_INFO *info, P7_TOPHITS *tophits, P7_TOPHITS *seed_hits);
#ifdef HMMER_THREADS
static int thread_loop(SPLICE_WORKER_INFO *info, P7_TOPHITS *tophits, P7_TOPHITS *seed_hits, int infocnt);
static void* splice_thread(void *arg);
#endif /*HMMER_THREADS*/

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

  /* Intialize data for threads */
  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(*info) * infocnt);

  for (i = 0; i < infocnt; ++i)
  {
	info[i].hmm            = p7_hmm_Clone(hmm);
    info[i].om             = p7_oprofile_Clone(om);
    info[i].gm             = p7_profile_Clone(gm);
    info[i].gm_fs          = p7_profile_fs_Clone(gm_fs);
	info[i].pli            = p7_splicepipeline_Create(go, 100, 100);
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

  graph->anchor_N = graph->th->N;
  graph->seqname = graph->th->hit[0]->name;

  return eslOK;
 
}

/*  Function: p7_splice_AddSeeds
 *  Synopsis: Add additional hits that passed the Forward filter to graph
 *
 *  Purpose : Find all hits in <seed_hits> that match the seqidx
 *            and revcomp of the <graph>, that are not already 
 *            added to the graph by p7_splice_AddAnchors(), and 
 *            are upsteam of one anchor node and downstream of 
 *            another and add them to the graph
 *
 * Returns:   <eslOK>.
 *
 */
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
  return eslOK;
 
}

/*  Function: p7_splice_ComputeAliScores
 *  Synopsis: Compute the per-postion scoring used for node and edge scores 
 *
 *  Purpose : For a protein to protein alignment - use the trace to 
 *            compute the emission and transtions scores for each 
 *            postions in the alignment. Sum these scores to produce 
 *            the aliscore and add bias adgustment if requested.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_splice_ComputeAliScores(P7_DOMAIN *dom, P7_TRACE *tr, ESL_DSQ *amino_dsq, const P7_PROFILE *gm, P7_BG *bg, float fs_prob, int do_bias)
{

  int status;
  int i, k, n, z;
  int z1, z2;
  int N;
  int first_i;
  float bias;

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
   * emissions scores are multiplied ny the stadard codon probability from the fs model */
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


/*  Function: p7_splice_ComputeAliScores_fs
 *  Synopsis: Compute the per-postion scoring used for node and edge scores
 *
 *  Purpose : For a translated and frameshift aware alignment - use the 
 *            trace to compute the emission and transtions scores for 
 *            each postions in the alignment. Sum these scores to produce 
 *            the aliscore and add bias adgustment if requested.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
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
    i = tr->i[z1];
    c = tr->c[z1];

    if (tr->st[z1] == p7T_M) {
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


/*  Function: p7_splice_HitUpstream
 *  Synopsis: Determine if one hit is upstream of another
 *
 *  Purpose : Check the sequence and hmm coordinates of two hits 
 *            to determine if one is upstream of another 
 *
 * Returns:   <TRUE> if <upstream> is indeed upstream of <downstream> and <FALSE> otherwise.
 *
 */
int
p7_splice_HitUpstream(P7_DOMAIN *upstream, P7_DOMAIN *downstream, int revcomp) 
{

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


/*  Function: p7_splice_HitiBetween
 *
 *  Synopsis: Determine if one hit is between to other his on the sequence
 *
 * Returns:   <TRUE> if <mid> is indeed between <up> and <down> and <FALSE> otherwise 
 *
 */
int
p7_splice_HitBetween(P7_DOMAIN *up, P7_DOMAIN *mid, P7_DOMAIN *down, int revcomp) 
{

  if (( revcomp  && up->iali <= mid->iali) ||
     ((!revcomp) && up->iali >= mid->iali))
    return FALSE;

  if (( revcomp  && mid->jali <= down->jali) ||
     ((!revcomp) && mid->jali >= down->jali))
    return FALSE;

  return TRUE;
  
}



/*****************************************************************
 * 2. Path Finding and Splicing
 *****************************************************************/

/*  Function: p7_splice_SpliceGraph
 *  Synopsis: Main splicing function
 *
 *  Purpose : Given a splice graph with BATH hits from a single target 
 *            sequence and strand, find and splice the best path(s) and 
 *            align the spliced sequnce to the model
 *
 * Returns:   <eslOK>. 
 *
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
  SPLICE_PATH       *orig_path;
  SPLICE_PATH       *copy_path;
  SPLICE_PATH       *spliced_path;
  ESL_SQ            *path_seq;
  SPLICE_GRAPH      *graph;
  SPLICE_PIPELINE   *pli;
  P7_PROFILE        *gm;
  ESL_SQFILE        *seq_file;

  path_seq         = NULL;

  graph      = info->graph;
  pli        = info->pli;
  gm         = info->gm;
  seq_file   = info->seq_file;

  printf("\nQuery %s Target %s strand %c seqidx %lld\n", gm->name, graph->seqname, (graph->revcomp ? '-' : '+'), graph->seqidx);
  fflush(stdout);

//printf("RECOVER\n");
//p7_splicegraph_DumpHits(stdout, graph);
//fflush(stdout);

  /* Create edges between original and recovered nodes */
  p7_splice_CreateUnsplicedEdges(graph, gm);
 
  /* Build paths from orignal hit nodes and edge so that every node appears in one and only one path */
  orig_path = p7_splicepath_GetBestPath(graph);

//p7_splicegraph_DumpGraph(stdout, graph);

  /* Find paths until no more exist with at least one anchor node */
  while(orig_path != NULL) {

    seq_min = ESL_MIN(orig_path->iali[0], orig_path->jali[orig_path->path_len-1]) - ALIGNMENT_EXT;
    seq_max = ESL_MAX(orig_path->iali[0], orig_path->jali[orig_path->path_len-1]) + ALIGNMENT_EXT;
    path_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, orig_path->revcomp, info);
    
//printf("FIRST PATH \n");
//p7_splicepath_Dump(stdout,orig_path);
//p7_splicepath_DumpScores(stdout,orig_path,graph);
//fflush(stdout);

    /* Create dopy of orig_path so orig)path does not get altered by p7_splice_SpliceExons() */
    copy_path = p7_splicepath_Clone(orig_path);
    spliced_path = p7_splice_SpliceExons(info, copy_path, path_seq);
    
//printf("SPLICED PATH\n");
//p7_splicepath_Dump(stdout,spliced_path);
//fflush(stdout);

    if(spliced_path != NULL) {

      /* Add additional nodes to the begining and end of spliced_path */
      p7_splice_ExtendPath(info->seeds, gm, orig_path, spliced_path, graph);

//printf("EXTEND PATH\n");
//p7_splicepath_Dump(stdout,spliced_path);
//fflush(stdout);

      /* If extension nodes were added, splice them */
      if(spliced_path->extension[0] || spliced_path->extension[spliced_path->path_len-1]) {
  
        esl_sq_Destroy(path_seq);
  
        seq_min = ESL_MIN(spliced_path->iali[0], spliced_path->jali[spliced_path->path_len-1]) - ALIGNMENT_EXT;
        seq_max = ESL_MAX(spliced_path->iali[0], spliced_path->jali[spliced_path->path_len-1]) + ALIGNMENT_EXT;
        path_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, spliced_path->revcomp, info);
       
        p7_splice_SpliceExtensions(info, spliced_path, path_seq);
        
      }
      else if (spliced_path->path_len == 1) p7_splice_SpliceSingle(info, spliced_path, path_seq);        
      
//printf("FINAL PATH\n");
//p7_splicepath_Dump(stdout,spliced_path);
//fflush(stdout);

      success = FALSE;

      if(spliced_path->path_len > 1) {

        /* Recheck node assignments to maxiimize the number of anchor nodes */
        for(h = 0; h < graph->anchor_N; h++) {
          if(!graph->node_in_graph[h]) continue;
          for(s = 0; s < spliced_path->path_len; s++) {
            if(spliced_path->node_id[s] >= graph->anchor_N) {
              if(p7_splicegraph_NodeOverlap(graph, h, spliced_path, s))
                spliced_path->node_id[s] = h;
            }
          }
        }
        /* Create the final alginement */  
        p7_splice_AlignSplicedPath(info, orig_path, spliced_path, path_seq, &success);

      }
 
      if(success) {
        /* Break edges that overlap the hit so that paths do not intertwine */
        hit_min = ESL_MIN(pli->hit->dcl->iali, pli->hit->dcl->jali); 
        hit_max = ESL_MAX(pli->hit->dcl->iali, pli->hit->dcl->jali); 
        p7_splice_EnforceBounds(graph, hit_min, hit_max); 
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
        if(spliced_path->path_len > 1) {
          path_min = ESL_MIN(orig_path->iali[0], orig_path->jali[orig_path->path_len-1]);
          path_max = ESL_MAX(orig_path->iali[0], orig_path->jali[orig_path->path_len-1]);
          p7_splice_EnforceBounds(graph, path_min, path_max);
        }
        for(s = 0; s < orig_path->path_len; s++) 
          graph->node_in_graph[orig_path->node_id[s]] = FALSE;
      }
    } 

    esl_sq_Destroy(path_seq);
    p7_splicepath_Destroy(orig_path);
    p7_splicepath_Destroy(copy_path);
    p7_splicepath_Destroy(spliced_path);
    p7_splicepipeline_Reuse(pli);

//p7_splicegraph_DumpHits(stdout, graph);
    orig_path = p7_splicepath_GetBestPath(graph);

//p7_splicegraph_DumpHits(stdout, graph); 
//p7_splicegraph_DumpGraph(stdout, graph);
    
  }
printf("\nComplete Query %s Target %s strand %c seqidx %ld\n", gm->name, graph->seqname, (graph->revcomp ? '-' : '+'), graph->seqidx);
  fflush(stdout);  
  return eslOK;

}


/*  Function: p7_splice_CreateUnsplicedEdges
 *  Synopsis: Add unspliced edges to splice graph
 *
 *  Purpose : Given a splice graph with nodes, add unspliced
 *            edges between any hits that are up/down stream of
 *            each other. Edge scores are zero unless the hits 
 *            overlap in amino positions, when the edge score is
 *            the minimum lost score to eliminate the overlap. 
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

      if(amino_gap_len > MAX_AMINO_GAP) continue;

      if(amino_gap_len > 0 && seq_gap_len < amino_gap_len) continue;

      if(th->hit[up]->dcl->ihmm >= th->hit[down]->dcl->jhmm ) { 

        /* If hits are upstream/ downstream on the sequence but are the opposite on th modle create a jump edge to prevent paths from skipping over other nodes */
        if( up < graph->anchor_N && down < graph->anchor_N) {
          edge = p7_splicegraph_AddEdge(graph, up, down);
          edge->splice_score   = (th->hit[up]->dcl->aliscore + th->hit[down]->dcl->aliscore) * -1; 
       
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

        /* If hits overlap, find the minimum lost score to remove the overlap */
        p7_splicegraph_AliScoreEdge(edge, gm, th->hit[up]->dcl, th->hit[down]->dcl); 
       
        edge->upstream_amino_end     = th->hit[up]->dcl->jhmm;
        edge->downstream_amino_start = th->hit[down]->dcl->ihmm; 
        edge->upstream_nuc_end       = th->hit[up]->dcl->jali;
        edge->downstream_nuc_start   = th->hit[down]->dcl->iali;

        /* If the edge has an hmm overlap and cost of eliminating that overlap is greater than the B->M entry 
         * for the downstream hit then these hits are better off seperate and we will remove the edge */
        if(edge->splice_score < -eslCONST_LOG2 + p7P_TSC(gm, th->hit[down]->dcl->ihmm-1, p7P_BM)) {
          graph->num_edges[up]--;
        }
      }   
    }
  } 
  
  return eslOK;

}



/*  Function: p7_splice_ExtendPath
 *  Synopsis: Add seed hits to beginning or end of path 
 *
 *  Purpose : Recover potential exons that were discarded by
 *            the BATH Forward filter from either the original 
 *            path of the seed_hits, and add them to the 
 *            beginning or end of the spliced path
 *
 * Returns:   <eslOK> on success.
 *
 */
int
p7_splice_ExtendPath(P7_TOPHITS *seed_hits, P7_PROFILE *gm, SPLICE_PATH *path, SPLICE_PATH *spliced_path, SPLICE_GRAPH *graph)
{
  int i,s,n;
  int skip;
  int seeds_in_graph;
  P7_HIT *first_hit;
  P7_HIT *last_hit; 
  P7_HIT *curr_hit;
  SPLICE_GRAPH *tmp_graph;
  SPLICE_PATH  *tmp_path;
  SPLICE_EDGE  *tmp_edge; 
  SPLICE_EDGE  *edge;
 
   
  /************************
  * EXTEND UP
  *************************/

  first_hit = graph->th->hit[spliced_path->node_id[0]];
 
  tmp_graph = p7_splicegraph_Create();
  tmp_graph->seqidx  = graph->seqidx;
  tmp_graph->revcomp = graph->revcomp;
  tmp_graph->seqname = graph->seqname; 

  /* Create a new graph for the first node and the possible upstream extensions */ 
  p7_splicegraph_CreateNodes(tmp_graph, 1);  
  p7_splicegraph_AddNode(tmp_graph, first_hit);

  tmp_graph->reportable[0]    = TRUE;
  tmp_graph->orig_hit_idx[0]  = spliced_path->node_id[0];
  tmp_graph->split_orig_id[0] = 0;

  tmp_graph->anchor_N = 1;
 
  /* Add any existing upstream extension nodes to the tmp graph */
  s = 0;
  while(path->extension[s]) {

    curr_hit = graph->th->hit[path->node_id[s]];
   
    p7_splicegraph_AddNode(tmp_graph, curr_hit);
    tmp_graph->reportable[tmp_graph->num_nodes-1]    = FALSE;
    tmp_graph->orig_hit_idx[tmp_graph->num_nodes-1]  = path->node_id[s];
    tmp_graph->split_orig_id[tmp_graph->num_nodes-1] = -1;
 
    /* Set tmp_nodes so we don't make new edges for them */
    if(graph->tmp_node[path->node_id[s]]) tmp_graph->tmp_node[tmp_graph->num_nodes-1] = TRUE;
    s++;
  }

  seeds_in_graph = tmp_graph->num_nodes;

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

      p7_splicegraph_AddNode(tmp_graph, curr_hit);
      tmp_graph->reportable[tmp_graph->num_nodes-1]    = FALSE;
      tmp_graph->orig_hit_idx[tmp_graph->num_nodes-1]  = -1;
      tmp_graph->split_orig_id[tmp_graph->num_nodes-1] = -1;  
    }
  }

  p7_splice_CreateExtensionEdges(graph, tmp_graph, gm); 
  tmp_path = p7_splicepath_GetBestPath_Extension(tmp_graph); 

  /* Set the most upstream hit in the tmp_path to start at the minimum 
   * hmm for all hits in tmp_path and add it to the original path  */
  for(s = tmp_path->path_len-2; s >= 0; s--) {

    curr_hit = tmp_graph->th->hit[tmp_path->node_id[s]];     

    /* If this is a new seed add it to the main graph */
    if(tmp_path->node_id[s] >= seeds_in_graph) {
      curr_hit->dcl->is_included = TRUE;
  
      p7_splicegraph_AddNode(graph, curr_hit);
   
      tmp_edge = p7_splicegraph_GetEdge(tmp_graph, tmp_path->node_id[s], tmp_path->node_id[s+1]);
  
      edge = p7_splicegraph_AddEdge(graph, graph->num_nodes-1, spliced_path->node_id[0]);
      edge->upstream_amino_end     = tmp_edge->upstream_amino_end;
      edge->downstream_amino_start = tmp_edge->downstream_amino_start;
      edge->upstream_nuc_end       = tmp_edge->upstream_nuc_end;
      edge->downstream_nuc_start   = tmp_edge->downstream_nuc_start;
      edge->splice_score           = tmp_edge->splice_score;
  
      p7_splicepath_Insert(spliced_path, 0); 
      spliced_path->node_id[0]   = graph->num_nodes-1; 
      spliced_path->extension[0] = TRUE;
      spliced_path->ihmm[0]      = tmp_path->ihmm[s];
      spliced_path->jhmm[0]      = tmp_path->jhmm[s];
      spliced_path->iali[0]      = tmp_path->iali[s];
      spliced_path->jali[0]      = tmp_path->jali[s];
      spliced_path->aliscore[0]  = tmp_path->aliscore[s];

    }
    else {

      if(!p7_splicegraph_EdgeExists(graph, tmp_graph->orig_hit_idx[tmp_path->node_id[s]], spliced_path->node_id[0])) { 
        tmp_edge = p7_splicegraph_GetEdge(tmp_graph, tmp_path->node_id[s], tmp_path->node_id[s+1]);
        edge = p7_splicegraph_AddEdge(graph, tmp_graph->orig_hit_idx[tmp_path->node_id[s]], spliced_path->node_id[0]);
        edge->upstream_amino_end     = tmp_edge->upstream_amino_end;
        edge->downstream_amino_start = tmp_edge->downstream_amino_start;
        edge->upstream_nuc_end       = tmp_edge->upstream_nuc_end;
        edge->downstream_nuc_start   = tmp_edge->downstream_nuc_start;
        edge->splice_score           = tmp_edge->splice_score;
      }

      p7_splicepath_Insert(spliced_path, 0); 
      spliced_path->node_id[0]   = tmp_graph->orig_hit_idx[tmp_path->node_id[s]]; 
      spliced_path->extension[0] = TRUE;
      spliced_path->ihmm[0]      = tmp_path->ihmm[s];
      spliced_path->jhmm[0]      = tmp_path->jhmm[s];
      spliced_path->iali[0]      = tmp_path->iali[s];
      spliced_path->jali[0]      = tmp_path->jali[s];
      spliced_path->aliscore[0]  = tmp_path->aliscore[s];
   }
  }

  /*Unset tmp_nodes so they don't get destroyed */
  for(i = 0; i < tmp_graph->num_nodes; i++)
    tmp_graph->tmp_node[i] = FALSE;

  p7_splicepath_Destroy(tmp_path);
  p7_splicegraph_Destroy(tmp_graph);


  /************************
  * EXTEND DOWN
  ************************/

  last_hit  = graph->th->hit[spliced_path->node_id[spliced_path->path_len-1]];

  tmp_graph = p7_splicegraph_Create();
  tmp_graph->seqidx  = graph->seqidx;
  tmp_graph->revcomp = graph->revcomp;
  tmp_graph->seqname = graph->seqname;

  p7_splicegraph_CreateNodes(tmp_graph, 1);
  p7_splicegraph_AddNode(tmp_graph, last_hit);
 
  tmp_graph->reportable[0]    = TRUE;
  tmp_graph->orig_hit_idx[0]  = spliced_path->node_id[spliced_path->path_len-1];
  tmp_graph->split_orig_id[0] = 0;

  tmp_graph->anchor_N = 1;

   /* Add any existing downstream extension nodes to the tmp graph then remove them from the path */
   s = path->path_len-1;
  while(path->extension[s]) {

    curr_hit = graph->th->hit[path->node_id[s]];

    p7_splicegraph_AddNode(tmp_graph, curr_hit);
    tmp_graph->reportable[tmp_graph->num_nodes-1]    = FALSE;
    tmp_graph->orig_hit_idx[tmp_graph->num_nodes-1]  = path->node_id[s];
    tmp_graph->split_orig_id[tmp_graph->num_nodes-1] = -1;

    /* Set tmp_nodes so we don't make new edges for them */
    if(graph->tmp_node[path->node_id[s]]) tmp_graph->tmp_node[tmp_graph->num_nodes-1] = TRUE;
    s--;
  } 

  seeds_in_graph = tmp_graph->num_nodes;

  /* Add new node from seed_hits */
  for(i = 0; i < seed_hits->N; i++) {

    curr_hit = &(seed_hits->unsrt[i]);

    /* skip seeds already added to a graph */
    if(curr_hit->dcl->is_included) continue;

    /* Is the seed hit on the same sequence and strand as the tmp_graph */
    if (curr_hit->seqidx != tmp_graph->seqidx) continue;
    if (tmp_graph->revcomp    && curr_hit->dcl->iali < curr_hit->dcl->jali) continue;
    if ((!tmp_graph->revcomp) && curr_hit->dcl->iali > curr_hit->dcl->jali) continue;

    /*If the seed is downstream of the last hit, add it to the tmp graph unless there are any anchor nodes between*/
    if(p7_splice_HitUpstream(last_hit->dcl, curr_hit->dcl, tmp_graph->revcomp)) {
      skip = FALSE;
      for(n = 0; n < graph->anchor_N; n++) {
        if(!graph->node_in_graph[n]) continue;
        if(p7_splice_HitUpstream(graph->th->hit[n]->dcl, curr_hit->dcl, tmp_graph->revcomp)) {
          if(p7_splice_HitBetween(last_hit->dcl, graph->th->hit[n]->dcl, curr_hit->dcl, tmp_graph->revcomp)) skip = TRUE;
        }
      }
      if(skip) continue;

      p7_splicegraph_AddNode(tmp_graph, curr_hit);
      tmp_graph->reportable[tmp_graph->num_nodes-1]    = FALSE;
      tmp_graph->orig_hit_idx[tmp_graph->num_nodes-1]  = -1;
      tmp_graph->split_orig_id[tmp_graph->num_nodes-1] = -1;
    }
  }   

  p7_splice_CreateExtensionEdges(graph, tmp_graph, gm); 


  tmp_path = p7_splicepath_GetBestPath_Extension(tmp_graph);

  /* Set the most downstream hit in the tmp_path to end at the maximum  
   * hmm for all hits in tmp_path and add it to the original path  */
  /* Add any seed hits from tmp_path to real path */
  for(s = 1; s < tmp_path->path_len; s++) { 
  
    curr_hit = tmp_graph->th->hit[tmp_path->node_id[s]];

    if(tmp_path->node_id[s] >= seeds_in_graph) {
      curr_hit->dcl->is_included = TRUE;
      p7_splicegraph_AddNode(graph, curr_hit);
  
 
      tmp_edge = p7_splicegraph_GetEdge(tmp_graph, tmp_path->node_id[s-1], tmp_path->node_id[s]);
       
      edge = p7_splicegraph_AddEdge(graph, spliced_path->node_id[spliced_path->path_len-1], graph->num_nodes-1);
      
      edge->upstream_amino_end     = tmp_edge->upstream_amino_end;
      edge->downstream_amino_start = tmp_edge->downstream_amino_start;
      edge->upstream_nuc_end       = tmp_edge->upstream_nuc_end;
      edge->downstream_nuc_start   = tmp_edge->downstream_nuc_start;
      edge->splice_score           = tmp_edge->splice_score;

      p7_splicepath_Insert(path, path->path_len);
      path->node_id[path->path_len-1]   = graph->num_nodes-1;
      path->extension[path->path_len-1] = TRUE;
      path->ihmm[path->path_len-1]      =  tmp_path->ihmm[s];
      path->jhmm[path->path_len-1]      =  tmp_path->jhmm[s];
      path->iali[path->path_len-1]      =  tmp_path->iali[s];
      path->jali[path->path_len-1]      =  tmp_path->jali[s];
      path->aliscore[path->path_len-1]  = tmp_path->aliscore[s];
 
    }
    else {

      if(!p7_splicegraph_EdgeExists(graph, spliced_path->node_id[spliced_path->path_len-1], tmp_graph->orig_hit_idx[tmp_path->node_id[s]])) {
        tmp_edge = p7_splicegraph_GetEdge(tmp_graph, tmp_path->node_id[s-1], tmp_path->node_id[s]);
        edge = p7_splicegraph_AddEdge(graph, spliced_path->node_id[spliced_path->path_len-1], tmp_graph->orig_hit_idx[tmp_path->node_id[s]]);
        edge->upstream_amino_end     = tmp_edge->upstream_amino_end;
        edge->downstream_amino_start = tmp_edge->downstream_amino_start;
        edge->upstream_nuc_end       = tmp_edge->upstream_nuc_end;
        edge->downstream_nuc_start   = tmp_edge->downstream_nuc_start;
        edge->splice_score           = tmp_edge->splice_score;
      }

      p7_splicepath_Insert(spliced_path, spliced_path->path_len);
      spliced_path->node_id[spliced_path->path_len-1]   = tmp_graph->orig_hit_idx[tmp_path->node_id[s]];
      spliced_path->extension[spliced_path->path_len-1] = TRUE;
      spliced_path->ihmm[spliced_path->path_len-1]      = tmp_path->ihmm[s];
      spliced_path->jhmm[spliced_path->path_len-1]      = tmp_path->jhmm[s];
      spliced_path->iali[spliced_path->path_len-1]      = tmp_path->iali[s];
      spliced_path->jali[spliced_path->path_len-1]      = tmp_path->jali[s];
      spliced_path->aliscore[spliced_path->path_len-1]  = tmp_path->aliscore[s];
   }
  }

    /*Unset tmp_nodes so they don't get destroyed */
  for(i = 0; i < tmp_graph->num_nodes; i++)
    tmp_graph->tmp_node[i] = FALSE;

  p7_splicepath_Destroy(tmp_path);
  p7_splicegraph_Destroy(tmp_graph);

  return eslOK; 

}


/*  Function: p7_splice_CreateExtensionEdges
 *  Synopsis: Add unspliced extension edges to a graph of potential extension nodes
 *
 *  Purpose : Given a splice graph with the first or last anchor 
 *            node from a path and a set of potential extention 
 *            nodes, add unspliced edges between any hits that 
 *            are up/down stream of each other. Edge scores are 
 *            zero unless the hits overlap in amino positions, 
 *            when the edge score is the minimum lost score to 
 *            eliminate the overlap.
 *
 * Returns:   <eslOK>.
 *
 */
int
p7_splice_CreateExtensionEdges(SPLICE_GRAPH *orig_graph, SPLICE_GRAPH *extension_graph, P7_PROFILE *gm) 
{
  int up, down;
  int seq_gap_len;
  int amino_gap_len;
  P7_TOPHITS  *th;
  SPLICE_EDGE *edge;
  SPLICE_EDGE *orig_edge;

  th = extension_graph->th;

  for(up = 0; up < extension_graph->num_nodes; up++) {
    
    for(down = 0; down < extension_graph->num_nodes; down++) {
      if(up == down) continue;

      /* Are hits upstream/downstream */
      if(!p7_splice_HitUpstream(th->hit[up]->dcl, th->hit[down]->dcl, extension_graph->revcomp)) continue;
    
      if(extension_graph->revcomp)
        seq_gap_len = th->hit[up]->dcl->jali - th->hit[down]->dcl->iali - 1;
      else
        seq_gap_len = th->hit[down]->dcl->iali - th->hit[up]->dcl->jali - 1;        

      if(seq_gap_len > MAX_INTRON_EXT) continue;
     
      amino_gap_len = th->hit[down]->dcl->ihmm - th->hit[up]->dcl->jhmm - 1;

      if(amino_gap_len > MAX_AMINO_GAP) continue;

      if(amino_gap_len > 0 && seq_gap_len < amino_gap_len) continue;


      /* If these nodes alreagy have an edge in the original graph, copy it */
      if(extension_graph->orig_hit_idx[up] >= 0 && extension_graph->orig_hit_idx[down] >= 0) {
 
        orig_edge =  p7_splicegraph_GetEdge(orig_graph, extension_graph->orig_hit_idx[up],  extension_graph->orig_hit_idx[down]);

        if(orig_edge != NULL) {

          edge = p7_splicegraph_AddEdge(extension_graph, up, down);
          edge->upstream_amino_end     = orig_edge->upstream_amino_end;
          edge->downstream_amino_start = orig_edge->downstream_amino_start;
          edge->upstream_nuc_end       = orig_edge->upstream_nuc_end;
          edge->downstream_nuc_start   = orig_edge->downstream_nuc_start;
  
          edge->i_start      = orig_edge->i_start;
          edge->k_start      = orig_edge->k_start;
          edge->i_end        = orig_edge->i_end;
          edge->k_end        = orig_edge->k_end;
          edge->next_i_start = orig_edge->next_i_start;
          edge->next_k_start = orig_edge->next_k_start;
        
          edge->splice_score = orig_edge->splice_score;  
        }
      } 
      else if(!extension_graph->tmp_node[up] && !extension_graph->tmp_node[down])
      {
        /* If hits overlap, find the minimum lost socre to remove the overlap */

        edge = p7_splicegraph_AddEdge(extension_graph, up, down);
        p7_splicegraph_AliScoreEdge(edge, gm, th->hit[up]->dcl, th->hit[down]->dcl); 

        edge->upstream_amino_end     = th->hit[up]->dcl->jhmm;
        edge->downstream_amino_start = th->hit[down]->dcl->ihmm; 
        edge->upstream_nuc_end       = th->hit[up]->dcl->jali;
        edge->downstream_nuc_start   = th->hit[down]->dcl->iali;
      }
    }
  } 

  
  return eslOK;

}


/*  Function: p7_splice_SpliceExons
 *  Synopsis: Find the best splice site (if any) between each pair of node in the path
 *
 *  Purpose : From the first to the last anchor node in the path, recover 
 *            the splice site by either performing the alignment or by
 *            retreveing the site from a pervious alignment. If no site 
 *            is found or the spliced alignment score is less that the 
 *            combined unspliced alignment sores, break the edge and 
 *            return null. Save all valid splice sites coodinates to 
 *            prevent repeated alignment.
 *
 *
 * Returns:   <SPLICE_PATH> on success and NULL on failure. 
 *
 */
SPLICE_PATH*
p7_splice_SpliceExons(SPLICE_WORKER_INFO *info, SPLICE_PATH *orig_path, ESL_SQ *path_seq)
{

  int i, s;
  int s_start, s_end;
  int i_start, i_end;
  int k_start, k_end;
  int next_i_start, next_k_start;
  float           ali_score;
  P7_PROFILE      *gm;
  SPLICE_EDGE     *edge;
  SPLICE_EDGE     *next_edge;
  SPLICE_PIPELINE *pli;
  SPLICE_PATH *tmp_path;
  SPLICE_PATH *ret_path;
  SPLICE_GRAPH *graph;
  
  graph = info->graph;
  pli   = info->pli;
  gm    = info->gm;

  ret_path   = NULL;
 
  for(s_start = 0; s_start < orig_path->path_len; s_start++) if(!orig_path->extension[s_start]) break;
  for(s_end = orig_path->path_len-1; s_end >= 0; s_end--)    if(!orig_path->extension[s_end])   break;

  if (s_start == s_end) {
    ret_path = p7_splicepath_Create(1);

    ret_path->node_id[0] = orig_path->node_id[s_start];

    ret_path->iali[0] = orig_path->iali[s_start];
    ret_path->jali[0] = orig_path->jali[s_start];
    ret_path->ihmm[0] = orig_path->ihmm[s_start];
    ret_path->jhmm[0] = orig_path->jhmm[s_start];

    ret_path->extension[0] = FALSE;

    ret_path->revcomp    = orig_path->revcomp;
    ret_path->frameshift = orig_path->frameshift;
    
    return ret_path;
  }

  next_i_start = next_k_start = 0; 
  /* Splice each pair or nodes from the first to the last anchor */  
  for(s = s_start+1; s < s_end+1; s++) {

    /* After the first pair of nodes the start coordinates are set by the previous search */
    k_start = next_k_start == 0 ? orig_path->ihmm[s-1] : next_k_start; 
    i_start = next_i_start == 0 ? orig_path->iali[s-1] : next_i_start;
    k_end   = orig_path->jhmm[s];
    i_end   = orig_path->jali[s];
  
    /* Check if we have spliced this pair of nodes before and retrieve data from edge */
    edge =  p7_splicegraph_GetEdge(graph, orig_path->node_id[s-1], orig_path->node_id[s]);

    if(i_start == edge->i_start && k_start == edge->k_start) {

      /* Create path and/or insert step */
      if(ret_path == NULL) {
        ret_path = p7_splicepath_Create(2);
        ret_path->iali[0] = i_start;
        ret_path->ihmm[0] = k_start;

        ret_path->extension[0] = FALSE;
      } 
      else p7_splicepath_Insert(ret_path, ret_path->path_len);

      /* Recover splice coordinates and add to path */
      ret_path->jali[ret_path->path_len-2] = edge->upstream_nuc_end;
      ret_path->jhmm[ret_path->path_len-2] = edge->upstream_amino_end;
      ret_path->iali[ret_path->path_len-1] = edge->downstream_nuc_start;
      ret_path->ihmm[ret_path->path_len-1] = edge->downstream_amino_start;
      ret_path->jali[ret_path->path_len-1] = i_end;
      ret_path->jhmm[ret_path->path_len-1] = k_end;

      ret_path->node_id[ret_path->path_len-2] = orig_path->node_id[s-1];
      ret_path->node_id[ret_path->path_len-1] = orig_path->node_id[s];

      ret_path->extension[ret_path->path_len-2] = FALSE;
      ret_path->extension[ret_path->path_len-1] = FALSE;   

      next_k_start = edge->next_k_start; 
      next_i_start = edge->next_i_start;
    
      continue;
    }
    else {
      edge->i_start = i_start;
      edge->k_start = k_start;
    }
    
    /* Covert sequence positions to sub sequence positions */
    if(orig_path->revcomp) {
      i_start = path_seq->n + path_seq->end - i_start;
      i_end   = path_seq->n + path_seq->end - i_end;
    }
    else  {
      i_start = i_start - path_seq->start + 1;
      i_end   = i_end   - path_seq->start + 1;
    }

    /* If the previous search has moved the start positions to after 
       the next steps end postions break the edge and return NULL */
    if(k_end <= k_start || i_end <= i_start) {

      edge->splice_score = -eslINFINITY;
      p7_splicepath_Destroy(ret_path);
      return NULL;
    }

    /* Align the pair of nodes to find the splice site */
    ali_score = -eslINFINITY; 

    tmp_path = p7_splice_AlignExons(info, orig_path, path_seq, s, i_start, i_end, k_start, k_end, &next_i_start, &next_k_start, &ali_score);

    /* If the alignment failed or the spliced aligment score is less than 
     * the score of the two hits plus the B->M penalty these exons are 
     * better off as seperate hits so we erase the edge and return NULL */ 
    if(tmp_path == NULL) {

      /* refetch edge in in case of realloc */
      edge =  p7_splicegraph_GetEdge(graph, orig_path->node_id[s-1], orig_path->node_id[s]); 
      edge->splice_score = -eslINFINITY;
      p7_splicepath_Destroy(ret_path);
      return NULL;
    } 

    /* Add new splice coordinates to path */
    if(ret_path == NULL) ret_path = p7_splicepath_Clone(tmp_path);
    else {
      ret_path->jali[ret_path->path_len-1] = tmp_path->jali[0];
      ret_path->jhmm[ret_path->path_len-1] = tmp_path->jhmm[0];
     
      ret_path->extension[ret_path->path_len-1] = FALSE; 
      for(i = 1; i < tmp_path->path_len; i++) {   
        p7_splicepath_Insert(ret_path, ret_path->path_len);
        ret_path->iali[ret_path->path_len-1] = tmp_path->iali[i];
        ret_path->jali[ret_path->path_len-1] = tmp_path->jali[i];
        ret_path->ihmm[ret_path->path_len-1] = tmp_path->ihmm[i];
        ret_path->jhmm[ret_path->path_len-1] = tmp_path->jhmm[i]; 

        ret_path->node_id[ret_path->path_len-1]   = tmp_path->node_id[i]; 
        ret_path->extension[ret_path->path_len-1] = FALSE;
      }
    } 

    /* When hits are merged we need to remove one from the path */
    if(tmp_path->path_len == 1) {
      if(s != s_end) {
   
        if(!p7_splicegraph_EdgeExists(graph, orig_path->node_id[s-1], orig_path->node_id[s+1])) {
          next_edge = p7_splicegraph_AddEdge(graph, orig_path->node_id[s-1], orig_path->node_id[s+1]);
          next_edge->upstream_amino_end     = orig_path->jhmm[s-1];
          next_edge->downstream_amino_start = orig_path->ihmm[s+1]; 
          next_edge->upstream_nuc_end       = orig_path->jali[s-1];
          next_edge->downstream_nuc_start   = orig_path->iali[s+1];
          if(!graph->tmp_node[orig_path->node_id[s-1]] && !graph->tmp_node[orig_path->node_id[s+1]])
            p7_splicegraph_AliScoreEdge(next_edge, gm, graph->th->hit[orig_path->node_id[s-1]]->dcl, graph->th->hit[orig_path->node_id[s+1]]->dcl);
        }
        p7_splicepath_Remove(orig_path, s);
        s_end--;
        s--;
      }
    }
    p7_gmx_Reuse(pli->vit);
    p7_splicepath_Destroy(tmp_path); 
  }
 

  ret_path->revcomp    = orig_path->revcomp;
  ret_path->frameshift = orig_path->frameshift;
  return ret_path;

}


/*  Function: p7_splice_SpliceExtensions
 *  Synopsis: Find the splice sites between the first or last anchor node in a path any extstion nodes
 *
 *  Purpose : If any extension nodes were found by p7_splice_ExtendPath()
 *            test whaether they can be successfull spliced to the first or
 *            last anchor node in the path. 
 *
 * Returns:   eslOK 
 *
 */
int
p7_splice_SpliceExtensions(SPLICE_WORKER_INFO *info, SPLICE_PATH *spliced_path, ESL_SQ *path_seq)
{

  int i;
  int i_start, i_end;
  int k_start, k_end;
  int next_i_end;
  int next_k_end;
  int s_start, s_end;
  int tmp_path_len;
  SPLICE_PIPELINE *pli;
  SPLICE_PATH  *tmp_path;
  SPLICE_GRAPH *graph;
  SPLICE_EDGE  *edge; 

  pli   = info->pli;
  graph = info->graph;

  for(s_start = 0; s_start < spliced_path->path_len; s_start++) if(!spliced_path->extension[s_start]) break;
  for(s_end = spliced_path->path_len-1; s_end >= 0; s_end--)    if(!spliced_path->extension[s_end])   break;

  /**************************
  * Downstream Extension
  *************************/
  next_i_end = next_k_end = 0;
  /* If there are downstream extension nodes see if they are recoved by spliced viterbi */
  if(s_end != spliced_path->path_len-1) {

    /* if last anchor node has an upstream splce site - get teh start coordinates from the edge */
    if(s_end == s_start) {
       k_start = spliced_path->ihmm[s_end];
       i_start = spliced_path->iali[s_end];
    }
    else {
      edge = p7_splicegraph_GetEdge(graph, spliced_path->node_id[s_end-1], spliced_path->node_id[s_end]);
      k_start = edge->next_k_start;
      i_start = edge->next_i_start;
    }

    k_end   = spliced_path->jhmm[spliced_path->path_len-1];
    i_end   = spliced_path->jali[spliced_path->path_len-1];
    
    /* Covert sequence positions to sub sequence positions */
    if(spliced_path->revcomp) {
      i_start = path_seq->n + path_seq->end - i_start;
      i_end   = path_seq->n + path_seq->end - i_end;
    }
    else  {
      i_start = i_start - path_seq->start + 1;
      i_end   = i_end   - path_seq->start + 1;
    } 

    /* Algin downstream extension region to see if splice site(s) are found */
    tmp_path = p7_splice_AlignExtendDown(info, spliced_path, path_seq, s_end, i_start, i_end, k_start, k_end, &next_i_end, &next_k_end);

    /* Remove unspliced downstream extensions */
    tmp_path_len = spliced_path->path_len;
    for(i = s_end+1; i < tmp_path_len; i++) 
      p7_splicepath_Remove(spliced_path, spliced_path->path_len-1);

    /* Add new spliced downstream extensions */
    if(tmp_path != NULL) {
      /*Insert new exons and recorded splice boundaries */   
      spliced_path->jali[spliced_path->path_len-1] = tmp_path->jali[0];
      spliced_path->jhmm[spliced_path->path_len-1] = tmp_path->jhmm[0];
  
      for(i = 1; i < tmp_path->path_len; i++) {
        p7_splicepath_Insert(spliced_path, spliced_path->path_len);
        spliced_path->iali[spliced_path->path_len-1] = tmp_path->iali[i];
        spliced_path->jali[spliced_path->path_len-1] = tmp_path->jali[i];
        spliced_path->ihmm[spliced_path->path_len-1] = tmp_path->ihmm[i];
        spliced_path->jhmm[spliced_path->path_len-1] = tmp_path->jhmm[i];
      }
      p7_splicepath_Destroy(tmp_path);
    }
    p7_gmx_Reuse(pli->vit);
  }

  /**************************
   * Upstream Extension 
   *************************/
  /* If there are upstream extension nodes see if they are recoved by spliced viterbi */ 
  if(s_start != 0) {

    k_start = spliced_path->ihmm[0]; 
    i_start = spliced_path->iali[0]; 
   
    /* if first anchor node has a down stream splce site - get the end coordinates from the edge */
    if(s_start == spliced_path->path_len-1) {  // only one anchor with no downstream extension

      k_end   = spliced_path->jhmm[s_start];
      i_end   = spliced_path->jali[s_start];
    }
    else if( s_end == s_start ) { // only one anchor with downstream extension
      k_end = next_k_end;
      i_end = next_i_end;
    } 
    else { // more than one anchor
      edge = p7_splicegraph_GetEdge(graph, spliced_path->node_id[s_start], spliced_path->node_id[s_start+1]); 

      k_end = edge->k_end;
      i_end = edge->i_end;
    }

    /* Covert sequence positions to sub sequence positions */
    if(spliced_path->revcomp) {
      i_start = path_seq->n + path_seq->end - i_start;
      i_end   = path_seq->n + path_seq->end - i_end;
    }
    else  {
      i_start = i_start - path_seq->start + 1; 
      i_end   = i_end   - path_seq->start + 1;
    }
  
    /* Algin upstream extension region to see if splice site(s) are found */
    tmp_path = p7_splice_AlignExtendUp(info, spliced_path, path_seq, s_start, i_start, i_end, k_start, k_end);

    /* Remove unspliced upstream extensions */
    for(i = 0; i < s_start; i++)
      p7_splicepath_Remove(spliced_path, 0);

    /* Add new spliced upstream exenstions */
    if(tmp_path != NULL) {
      spliced_path->iali[0] = tmp_path->iali[tmp_path->path_len-1];
      spliced_path->ihmm[0] = tmp_path->ihmm[tmp_path->path_len-1];   
   
      for(i = tmp_path->path_len-2; i >= 0; i--) {

        p7_splicepath_Insert(spliced_path, 0);
        spliced_path->iali[0] = tmp_path->iali[i];
        spliced_path->jali[0] = tmp_path->jali[i];
        spliced_path->ihmm[0] = tmp_path->ihmm[i];
        spliced_path->jhmm[0] = tmp_path->jhmm[i];

      }
      p7_splicepath_Destroy(tmp_path);
    }
    p7_gmx_Reuse(pli->vit);
  }

  return eslOK;

}


/*  Function: p7_splice_SpliceSingle
 *  Synopsis: Find any internal spilce site that might exist in a single node path
 *
 *  Purpose : Since two or more exons with short introns between them can sometimes 
 *            be aligned as a single hit we check single node paths for splice sites
 *
 *
 * Returns:   eslOK
 *
 */
int
p7_splice_SpliceSingle(SPLICE_WORKER_INFO *info, SPLICE_PATH *spliced_path, ESL_SQ *path_seq)
{

  int i;
  int i_start, i_end;
  int k_start, k_end;
  SPLICE_PIPELINE *pli;
  SPLICE_PATH *tmp_path;
  
  pli   = info->pli;

  i_start = spliced_path->iali[0];
  i_end   = spliced_path->jali[0];
  k_start = spliced_path->ihmm[0];
  k_end   = spliced_path->jhmm[0];
   
  
  /* Covert sequence positions to sub sequence positions */
  if(spliced_path->revcomp) {
    i_start = path_seq->n + path_seq->end - i_start;
    i_end   = path_seq->n + path_seq->end - i_end;
  }
  else  {
    i_start = i_start - path_seq->start + 1;
    i_end   = i_end   - path_seq->start + 1;
  }
  /* Find any splice sites */
  tmp_path = p7_splice_AlignSingle(info, spliced_path, path_seq, i_start, i_end, k_start, k_end);

  /* If more than one exon was found addd them and their splice sites to the path */
  if(tmp_path != NULL) {
    spliced_path->jali[0] = tmp_path->jali[0];
    spliced_path->jhmm[0] = tmp_path->jhmm[0];

    for(i = 1; i < tmp_path->path_len; i++) {
      p7_splicepath_Insert(spliced_path, spliced_path->path_len);
      spliced_path->iali[spliced_path->path_len-1] = tmp_path->iali[i];
      spliced_path->jali[spliced_path->path_len-1] = tmp_path->jali[i];
      spliced_path->ihmm[spliced_path->path_len-1] = tmp_path->ihmm[i];
      spliced_path->jhmm[spliced_path->path_len-1] = tmp_path->jhmm[i];
      
      spliced_path->node_id[spliced_path->path_len-1] = spliced_path->node_id[0];
    }

     p7_splicepath_Destroy(tmp_path);
  }
  
  p7_gmx_Reuse(pli->vit);

  return eslOK;
  
}


/*  Function: p7_splice_AlignExons
 *  Synopsis: Use global spliced Viterbi to find the splice site between two exons 
 *
 *  Purpose : Align two neighboring exons from a splice path to find their 
 *            splice sites (if any) and any new exons not in the path. 
 *            Add new exons to the graph, and record spliced and unspliced 
 *            coordinates.
 *
 * Returns:   <SPLICE_PATH> on success and NULL on failure.
 *
 */
SPLICE_PATH*
p7_splice_AlignExons(SPLICE_WORKER_INFO *info, SPLICE_PATH *orig_path, ESL_SQ *path_seq, int down, int i_start, int i_end, int k_start, int k_end, int *next_i_start, int *next_k_start,float *ali_score)
{
 
  int         L = i_end - i_start + 1;
  int         M = k_end - k_start + 1;
  int         up = down - 1;
  int         s, y, z;
  int         z1, z2;
  int         intron_cnt;
  int         step_cnt;
  int         start_new;
  P7_GMX      *gx; 
  P7_TRACE    *tr;
  SPLICE_PATH *ret_path;
  SPLICE_PATH *tmp_path;
  SPLICE_EDGE *edge;
  SPLICE_GRAPH *graph;
  SPLICE_PIPELINE *pli;
  P7_FS_PROFILE *gm_fs;
  ESL_GENCODE *gcode;
  P7_HIT      *hit;

  graph = info->graph;
  pli   = info->pli;
  gm_fs = info->gm_fs;
  gcode = info->gcode;

  gx = pli->vit;

  ret_path = NULL;
  tmp_path = NULL;

  p7_gmx_sp_GrowTo(pli->vit, M, L, L);
  p7_splicepipline_GrowIndex(pli->sig_idx, M, L, ALIGNMENT_MODE);
  p7_fs_ReconfigLength(gm_fs, L);
   
  p7_spliceviterbi_translated_global(pli, path_seq->dsq, gcode, gm_fs, pli->vit, i_start, i_end, k_start, k_end);

  /* If the hits were in different frames and no splice site was able to pull score 
   * from the upstream frame to the downstream frame the spliceing is a failure */
  if(gx->xmx[L*p7G_NXCELLS+p7G_C] == -eslINFINITY) return NULL; 

  tr = p7_trace_fs_Create();
  p7_splicevitebi_translated_trace(pli, path_seq->dsq, gcode, gm_fs, pli->vit, tr, i_start, i_end, k_start, k_end, ali_score);

  /* Find number of introns in trace */
  intron_cnt = 0;
  for(z = 0; z < tr->N; z++)
    if(tr->st[z] == p7T_P) intron_cnt++;

  /* Find first M state - start of first hit */
  for(z1 = 0; z1 < tr->N; z1++) if(tr->st[z1] == p7T_M) break;
  if (z1 == tr->N) {
    p7_trace_fs_Destroy(tr);
    return NULL;
  }

  /* Find last M state state - end of last hit */
  for(z2 = tr->N-1; z2 >= 0; z2--) if(tr->st[z2] == p7T_M) break;
  if (z2 == -1) {
    p7_trace_fs_Destroy(tr);
    return NULL;
  }

  tmp_path = p7_splicepath_Create(intron_cnt+1);
  ret_path = p7_splicepath_Create(intron_cnt+1);

  step_cnt = 0;
  start_new = TRUE;
  
  z = z1;
  while(z <= z2) {

    if(start_new) {
      
      /* Save z value - currently set to first M state in exon */
      y = z;

      /*Find end of exon */
      while(tr->st[z] != p7T_P && tr->st[z] != p7T_E) z++;
      if(tr->st[z] == p7T_E) while(tr->st[z] != p7T_M) z--;
      else z--;
  
      
      /* If this is the first step in path the i coords are set my the first M state at
       * trace postion y, otherwise they are set by the P state at trace postion y-1 */
      tmp_path->node_id[step_cnt] = -1;
      ret_path->node_id[step_cnt] = -1;
      if(step_cnt == 0) {
        tmp_path->iali[step_cnt] = tr->i[y] - tr->c[y] + 1;
        tmp_path->ihmm[step_cnt] = tr->k[y];

        ret_path->iali[step_cnt] = tr->i[y] - tr->c[y] + 1; 
        ret_path->ihmm[step_cnt] = tr->k[y];        

      }
      else {
        /* Determine which splice codon type we have and assign the P state (y-1) amino to the correct exon */
        
        if( tr->c[y-1] == 0) {
          ret_path->iali[step_cnt]   = tr->i[y-1] - 2; 
          ret_path->ihmm[step_cnt]   = tr->k[y-1];
        }
        else if( tr->c[y-1] == 1) { 
          ret_path->iali[step_cnt]   = tr->i[y-1] - 1;
          ret_path->ihmm[step_cnt]   = tr->k[y-1];
        }
        else {
          ret_path->iali[step_cnt]   = tr->i[y-1];
          ret_path->ihmm[step_cnt]   = tr->k[y];
          ret_path->jhmm[step_cnt-1] = tr->k[y-1];
        }
        tmp_path->iali[step_cnt] = tr->i[y] - tr->c[y] + 1;
        tmp_path->ihmm[step_cnt] = tr->k[y];
      }

      tmp_path->jhmm[step_cnt] = tr->k[z];
      ret_path->jhmm[step_cnt] = tr->k[z];

      if(step_cnt == ret_path->path_len-1) {
        tmp_path->jali[step_cnt]   = tr->i[z];
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

        tmp_path->jali[step_cnt] = tr->i[z];
      }

      ret_path->extension[step_cnt] = FALSE;

      step_cnt++;

      start_new = FALSE;
    }
    
    z++;
    if(tr->st[z] == p7T_M) start_new = TRUE;

  }

  tmp_path->revcomp  = orig_path->revcomp;
  ret_path->revcomp  = orig_path->revcomp;

  /* Convert to true coordinates */
  for(s = 0; s < ret_path->path_len; s++) {
    if(orig_path->revcomp) {
      tmp_path->iali[s] = path_seq->n - tmp_path->iali[s] + path_seq->end;
      tmp_path->jali[s] = path_seq->n - tmp_path->jali[s] + path_seq->end;

      ret_path->iali[s] = path_seq->n - ret_path->iali[s] + path_seq->end;
      ret_path->jali[s] = path_seq->n - ret_path->jali[s] + path_seq->end;
    }
    else {
      tmp_path->iali[s] = path_seq->start + tmp_path->iali[s] - 1;
      tmp_path->jali[s] = path_seq->start + tmp_path->jali[s] - 1;

      ret_path->iali[s] = path_seq->start + ret_path->iali[s] - 1;
      ret_path->jali[s] = path_seq->start + ret_path->jali[s] - 1;
    }
  }


  /* Assign node ids */
  if(tmp_path->path_len == 1) {
    tmp_path->node_id[0] = orig_path->node_id[up];
    ret_path->node_id[0] = orig_path->node_id[up];

    /*If up and down merged into one hit beark the edge betwen up and down */
    edge = p7_splicegraph_GetEdge(graph, orig_path->node_id[up], orig_path->node_id[down]);
    edge->splice_score = -eslINFINITY;
  }
  else {
    tmp_path->node_id[0]                    = orig_path->node_id[up];
    tmp_path->node_id[tmp_path->path_len-1] = orig_path->node_id[down];

    ret_path->node_id[0]                    = orig_path->node_id[up];
    ret_path->node_id[ret_path->path_len-1] = orig_path->node_id[down];
  }

  /*If we have new nodes we need to break the old edges so that the new edges are used in future paths */
  if(tmp_path->path_len > 2) {
    edge = p7_splicegraph_GetEdge(graph, orig_path->node_id[up], orig_path->node_id[down]); 
    edge->splice_score = -eslINFINITY;   
  }

  /* Now that we have node assignments we
   * 1. can add nodes with no node assignment
   * 2. add edges that don't currently exist
   * 3. record start and next start coordinates from tmp_path onto edges
   * 4. record splice sites from ret_path onto edges
   */
  for(s = 0; s < tmp_path->path_len; s++) {
    
    if(tmp_path->node_id[s] == -1) {
  
      /* 1. add nodes */
      hit          = p7_hit_Create_empty();
      hit->dcl     = p7_domain_Create_empty();
      hit->dcl->tr = p7_trace_fs_Create();

      hit->dcl->iali = tmp_path->iali[s];
      hit->dcl->jali = tmp_path->jali[s];
      hit->dcl->ihmm = tmp_path->ihmm[s];
      hit->dcl->jhmm = tmp_path->jhmm[s];

      hit->dcl->aliscore = 1.; // any positive score ensures it is included in next path

      p7_splicegraph_AddNode(graph, hit);
      
      graph->tmp_node[graph->num_nodes-1] = TRUE;
 
      tmp_path->node_id[s] = graph->num_nodes-1;
      ret_path->node_id[s] = graph->num_nodes-1;

      /* Connect to upstream node (if any) */
      if(s != 0) {
        /* 2. add edges */
        edge = p7_splicegraph_AddEdge(graph, tmp_path->node_id[s-1], tmp_path->node_id[s]); 
        
        /* 3. add start coords */
        edge->i_start = tmp_path->iali[s-1];
        edge->k_start = tmp_path->ihmm[s-1];

        edge->i_end = tmp_path->jali[s-1];
        edge->k_end = tmp_path->jhmm[s-1];

        edge->next_i_start = tmp_path->iali[s];
        edge->next_k_start = tmp_path->ihmm[s];
        
        /* 4. Add splice coords */
        edge->upstream_nuc_end       = ret_path->jali[s-1];
        edge->upstream_amino_end     = ret_path->jhmm[s-1];
        edge->downstream_nuc_start   = ret_path->iali[s];
        edge->downstream_amino_start = ret_path->ihmm[s];

      }
    }
    else if(s != 0) {
      edge = p7_splicegraph_GetEdge(graph, tmp_path->node_id[s-1], tmp_path->node_id[s]);

      if(edge == NULL) {
        /* 2. add edges */
        edge = p7_splicegraph_AddEdge(graph, tmp_path->node_id[s-1], tmp_path->node_id[s]);

        /* 3. add start coords */
        edge->i_start = tmp_path->iali[s-1];
        edge->k_start = tmp_path->ihmm[s-1];

        edge->i_end = tmp_path->jali[s-1];
        edge->k_end = tmp_path->jhmm[s-1];

        edge->next_i_start = tmp_path->iali[s];
        edge->next_k_start = tmp_path->ihmm[s];

        /* 4. Add splice coords */
        edge->upstream_nuc_end       = ret_path->jali[s-1];
        edge->upstream_amino_end     = ret_path->jhmm[s-1];
        edge->downstream_nuc_start   = ret_path->iali[s];
        edge->downstream_amino_start = ret_path->ihmm[s];

      }
      else {
        
        /* 3. add start coords */
        edge->i_start = tmp_path->iali[s-1];
        edge->k_start = tmp_path->ihmm[s-1];

        edge->i_end = tmp_path->jali[s-1];
        edge->k_end = tmp_path->jhmm[s-1];
       
        edge->next_i_start = tmp_path->iali[s];
        edge->next_k_start = tmp_path->ihmm[s];

        /* 4. Add splice coords */
        edge->upstream_nuc_end       = ret_path->jali[s-1];
        edge->upstream_amino_end     = ret_path->jhmm[s-1];
        edge->downstream_nuc_start   = ret_path->iali[s];
        edge->downstream_amino_start = ret_path->ihmm[s];
      }
    } 
  }

  *next_k_start = tmp_path->ihmm[tmp_path->path_len-1];
  *next_i_start = tmp_path->iali[tmp_path->path_len-1];

  p7_trace_fs_Destroy(tr);
  p7_splicepath_Destroy(tmp_path);

  return ret_path;

}


/*  Function: p7_splice_AlignExtendDown
 *  Synopsis: Use semi-global spliced Viterbi to find exons downstream of an anchor node
 *
 *  Purpose : Align the last anchor node in a path and any downstream estension 
 *            nodes to find splice sites (if any) and any new exons not in the path.
 *            Add new exons to the graph, and record spliced and unspliced
 *            coordinates.
 *
 * Returns:   <SPLICE_PATH> on success and NULL on failure.
 *
 */
SPLICE_PATH*
p7_splice_AlignExtendDown(SPLICE_WORKER_INFO *info, SPLICE_PATH *spliced_path, ESL_SQ *path_seq, int s_end, int i_start, int i_end, int k_start, int k_end, int *next_i_end, int *next_k_end)
{
 
  int         L = i_end - i_start + 1;
  int         M = k_end - k_start + 1;
  int         s, y, z;
  int         s1, s2;
  int         z1, z2;
  int         intron_cnt;
  int         step_cnt;
  int         start_new;
  float       ali_score;
  P7_TRACE    *tr;
  SPLICE_PATH *ret_path; 
  SPLICE_PATH *tmp_path;
  SPLICE_EDGE *edge;
  P7_HIT      *hit;
  SPLICE_GRAPH *graph;
  SPLICE_PIPELINE *pli;
  P7_FS_PROFILE *gm_fs;
  ESL_GENCODE *gcode;

  graph = info->graph;
  pli   = info->pli;
  gm_fs = info->gm_fs;
  gcode = info->gcode;  

  p7_gmx_sp_GrowTo(pli->vit, M, L, L);
  p7_splicepipline_GrowIndex(pli->sig_idx, M, L, ALIGNMENT_MODE);
  p7_fs_ReconfigLength(gm_fs, L);
  
   p7_spliceviterbi_translated_semiglobal_extenddown(pli, path_seq->dsq, gcode, gm_fs, pli->vit, i_start, i_end, k_start, k_end);

  tr = p7_trace_fs_Create();
  p7_splicevitebi_translated_trace(pli, path_seq->dsq, gcode, gm_fs, pli->vit, tr, i_start, i_end, k_start, k_end, &ali_score);
  //p7_trace_fs_Dump(stdout, tr, NULL, NULL, NULL);

  /* Find number of introns in trace */
  intron_cnt = 0;
  for(z = 0; z < tr->N; z++)
    if(tr->st[z] == p7T_P) intron_cnt++;

  if(intron_cnt == 0) {
    p7_trace_fs_Destroy(tr);
    return NULL;
  }

  /* Find first M state - start of first hit */
  for(z1 = 0; z1 < tr->N; z1++) if(tr->st[z1] == p7T_M) break;
  if (z1 == tr->N) {
    p7_trace_fs_Destroy(tr);
    return NULL;
  }

  /* Find last M state state - end of last hit */
  for(z2 = tr->N-1; z2 >= 0; z2--) if(tr->st[z2] == p7T_M) break;
  if (z2 == -1) {
    p7_trace_fs_Destroy(tr);
    return NULL;
  }

  ret_path = p7_splicepath_Create(intron_cnt+1);
  tmp_path = p7_splicepath_Create(intron_cnt+1);
  step_cnt = 0;
  start_new = TRUE;
  
  z = z1;
  while(z <= z2) {

    if(start_new) {
      
      /* Save z value - currently set to first M state in exon */
      y = z;

      /*Find end of exon */
      while(tr->st[z] != p7T_P && tr->st[z]  != p7T_E) z++;
      if(tr->st[z] == p7T_E) while(tr->st[z] != p7T_M) z--;
      else z--;
  
      /* If this is the first step in path the i coords are set my the first M state at
       * trace postion y, otherwise they are set by the P state at trace postion y-1 */
      ret_path->node_id[step_cnt] = -1;
      tmp_path->node_id[step_cnt] = -1;
      if(step_cnt == 0) {
        ret_path->iali[step_cnt] = tr->i[y] - tr->c[y] + 1;
        ret_path->ihmm[step_cnt] = tr->k[y];        

        tmp_path->iali[step_cnt] = tr->i[y] - tr->c[y] + 1;
        tmp_path->ihmm[step_cnt] = tr->k[y];
      }
      else {
        /* Determine which splice codon type we have and assign the P state (y-1) amino to the correct exon */

        if( tr->c[y-1] == 0) {
          ret_path->iali[step_cnt]   = tr->i[y-1] - 2;
          ret_path->ihmm[step_cnt]   = tr->k[y-1];
        }
        else if( tr->c[y-1] == 1) {
          ret_path->iali[step_cnt]   = tr->i[y-1] - 1;
          ret_path->ihmm[step_cnt]   = tr->k[y-1];
        }
        else {
          ret_path->iali[step_cnt]   = tr->i[y-1];
          ret_path->ihmm[step_cnt]   = tr->k[y];
          ret_path->jhmm[step_cnt-1] = tr->k[y-1];
        }
        tmp_path->iali[step_cnt] = tr->i[y] - tr->c[y] + 1;
        tmp_path->ihmm[step_cnt] = tr->k[y];
      }
     
      ret_path->jhmm[step_cnt] = tr->k[z];

      if(step_cnt == ret_path->path_len-1) {
        ret_path->jali[step_cnt]   = tr->i[z];
        tmp_path->jali[step_cnt]   = tr->i[z];
      }
      else {
        /* Determine which splice codon type we have */
        if( tr->c[z+1] == 0)
          ret_path->jali[step_cnt] = tr->i[z];
        else if( tr->c[z+1] == 1)
          ret_path->jali[step_cnt] = tr->i[z] + 1;
        else if( tr->c[z+1] == 2)
          ret_path->jali[step_cnt] = tr->i[z] + 2;
   
        tmp_path->jali[step_cnt] = tr->i[z];
      }

      if(step_cnt == 0) {
        *next_i_end = tr->i[z];
        *next_k_end = tr->k[z];
      }

      step_cnt++;

      start_new = FALSE;
    }
    
    z++;
    if(tr->st[z] == p7T_M) start_new = TRUE;

  }
  
  ret_path->revcomp = spliced_path->revcomp;
  tmp_path->revcomp = spliced_path->revcomp;

  /* Convert to true coordinates */
  for(s = 0; s < ret_path->path_len; s++) {
    if(spliced_path->revcomp) {
      ret_path->iali[s] = path_seq->n - ret_path->iali[s] + path_seq->end;
      ret_path->jali[s] = path_seq->n - ret_path->jali[s] + path_seq->end;

      tmp_path->iali[s] = path_seq->n - tmp_path->iali[s] + path_seq->end;
      tmp_path->jali[s] = path_seq->n - tmp_path->jali[s] + path_seq->end;
     
    }
    else {
      ret_path->iali[s] = path_seq->start + ret_path->iali[s] - 1;
      ret_path->jali[s] = path_seq->start + ret_path->jali[s] - 1;

      tmp_path->iali[s] = path_seq->start + tmp_path->iali[s] - 1;
      tmp_path->jali[s] = path_seq->start + tmp_path->jali[s] - 1;
 
    }
  }

  /* Assign node ids */
  ret_path->node_id[0] = s_end;
  tmp_path->node_id[0] = s_end;

  for(s1 = 1; s1 < tmp_path->path_len; s1++) {
    for(s2 = s_end+1; s2 < spliced_path->path_len; s2++) {
      if(p7_splicegraph_NodeOverlap(graph, spliced_path->node_id[s2], tmp_path, s1)) {
        ret_path->node_id[s1] = spliced_path->node_id[s2];
        tmp_path->node_id[s1] = spliced_path->node_id[s2]; 
      }
    } 
  }

  /* Now that we have node assignments we
   * 1. can add nodes with no node assignment
   * 2. add edges that don't currently exist
   * 3. record start and next start coordinates from tmp_path onto edges
   * 4. record splice sites from ret_path onto edges
   */
  for(s = 1; s < tmp_path->path_len; s++) {   
    if(tmp_path->node_id[s] == -1) {  
      hit          = p7_hit_Create_empty();
      hit->dcl     = p7_domain_Create_empty();
      hit->dcl->tr = p7_trace_fs_Create();

      hit->dcl->iali = tmp_path->iali[s];
      hit->dcl->jali = tmp_path->jali[s];
      hit->dcl->ihmm = tmp_path->ihmm[s];
      hit->dcl->jhmm = tmp_path->jhmm[s];

      hit->dcl->aliscore = 1.; // any positive score ensures it is included in next path

      p7_splicegraph_AddNode(graph, hit);

      graph->tmp_node[graph->num_nodes-1] = TRUE;

      tmp_path->node_id[s] = graph->num_nodes-1;
      ret_path->node_id[s] = graph->num_nodes-1;

      edge = p7_splicegraph_AddEdge(graph, tmp_path->node_id[s-1], tmp_path->node_id[s]);

        /* 3. add start coords */
        edge->i_start = tmp_path->iali[s-1];
        edge->k_start = tmp_path->ihmm[s-1];

        edge->i_end = tmp_path->jali[s-1];
        edge->k_end = tmp_path->jhmm[s-1];

        edge->next_i_start = tmp_path->iali[s];
        edge->next_k_start = tmp_path->ihmm[s];

        /* 4. Add splice coords */
        edge->upstream_nuc_end       = ret_path->jali[s-1];
        edge->upstream_amino_end     = ret_path->jhmm[s-1];
        edge->downstream_nuc_start   = ret_path->iali[s];
        edge->downstream_amino_start = ret_path->ihmm[s];
    }
    else {
      edge = p7_splicegraph_GetEdge(graph, tmp_path->node_id[s-1], tmp_path->node_id[s]);

      /* 3. add start coords */
      edge->i_start = tmp_path->iali[s-1];
      edge->k_start = tmp_path->ihmm[s-1];

      edge->i_end = tmp_path->jali[s-1];
      edge->k_end = tmp_path->jhmm[s-1];

      edge->next_i_start = tmp_path->iali[s];
      edge->next_k_start = tmp_path->ihmm[s];

      /* 4. Add splice coords */
      edge->upstream_nuc_end       = ret_path->jali[s-1];
      edge->upstream_amino_end     = ret_path->jhmm[s-1];
      edge->downstream_nuc_start   = ret_path->iali[s];
      edge->downstream_amino_start = ret_path->ihmm[s];

    }
  }


  if(spliced_path->revcomp) *next_i_end = path_seq->n - *next_i_end + path_seq->end;
  else                      *next_i_end = path_seq->start + *next_i_end - 1;

  p7_trace_fs_Destroy(tr);
  p7_splicepath_Destroy(tmp_path);

  return ret_path;

}


/*  Function: p7_splice_AlignExtendUp
 *  Synopsis: Use semi-global spliced Viterbi to find exons upstream of an anchor node
 *
 *  Purpose : Align the first anchor node in a path and any upstream extension
 *            nodes to find splice sites (if any) and any new exons not in the path.
 *            Add new exons to the graph, and record spliced and unspliced
 *            coordinates.
 *
 * Returns:   <SPLICE_PATH> on success and NULL on failure.
 *
 */
SPLICE_PATH*
p7_splice_AlignExtendUp(SPLICE_WORKER_INFO *info, SPLICE_PATH *spliced_path, ESL_SQ *path_seq, int s_start, int i_start, int i_end, int k_start, int k_end)
{
 
  int         L = i_end - i_start + 1;
  int         M = k_end - k_start + 1;
  int         s, y, z;
  int         z1, z2;
  int         s1, s2;
  int         intron_cnt;
  int         step_cnt;
  int         start_new;
  float       ali_score;
  P7_TRACE    *tr;
  SPLICE_PATH *ret_path; 
  SPLICE_PATH *tmp_path;
  SPLICE_EDGE *edge;
  P7_HIT      *hit;
  SPLICE_GRAPH *graph;
  SPLICE_PIPELINE *pli;
  P7_FS_PROFILE *gm_fs;
  ESL_GENCODE *gcode;

  graph = info->graph;
  pli   = info->pli;
  gm_fs = info->gm_fs;
  gcode = info->gcode;

  p7_gmx_sp_GrowTo(pli->vit, M, L, L);
  p7_splicepipline_GrowIndex(pli->sig_idx, M, L, ALIGNMENT_MODE);
  p7_fs_ReconfigLength(gm_fs, L);
  
  p7_spliceviterbi_translated_semiglobal_extendup(pli, path_seq->dsq, gcode, gm_fs, pli->vit, i_start, i_end, k_start, k_end);

  tr = p7_trace_fs_Create();
  p7_splicevitebi_translated_trace(pli, path_seq->dsq, gcode, gm_fs, pli->vit, tr, i_start, i_end, k_start, k_end, &ali_score);
  //p7_trace_fs_Dump(stdout, tr, NULL, NULL, NULL);

  /* Find number of introns in trace */
  intron_cnt = 0;
  for(z = 0; z < tr->N; z++)
    if(tr->st[z] == p7T_P) intron_cnt++;

  if(intron_cnt == 0) {
    p7_trace_fs_Destroy(tr);
    return NULL;
  }

  /* Find first M state - start of first hit */
  for(z1 = 0; z1 < tr->N; z1++) if(tr->st[z1] == p7T_M) break;
  if (z1 == tr->N) {
    p7_trace_fs_Destroy(tr);
    return NULL;
  }

  /* Find last M state state - end of last hit */
  for(z2 = tr->N-1; z2 >= 0; z2--) if(tr->st[z2] == p7T_M) break;
  if (z2 == -1) {
    p7_trace_fs_Destroy(tr);
    return NULL;
  }

  ret_path = p7_splicepath_Create(intron_cnt+1);
  tmp_path = p7_splicepath_Create(intron_cnt+1);

  step_cnt = 0;
  start_new = TRUE;
  
  z = z1;
  while(z <= z2) {

    if(start_new) {
      
      /* Save z value - currently set to first M state in exon */
      y = z;

      /*Find end of exon */
      while(tr->st[z] != p7T_P && tr->st[z]  != p7T_E) z++;
      if(tr->st[z] == p7T_E) while(tr->st[z] != p7T_M) z--;
      else z--;
  
      /* If this is the first step in path the i coords are set my the first M state at
       * trace postion y, otherwise they are set by the P state at trace postion y-1 */
      ret_path->node_id[step_cnt] = -1;
      tmp_path->node_id[step_cnt] = -1;
      if(step_cnt == 0) {
        ret_path->iali[step_cnt] = tr->i[y] - tr->c[y] + 1;
        ret_path->ihmm[step_cnt] = tr->k[y];        

        tmp_path->iali[step_cnt] = tr->i[y] - tr->c[y] + 1;
        tmp_path->ihmm[step_cnt] = tr->k[y];
      }
      else {
        /* Determine which splice codon type we have and assign the P state (y-1) amino to the correct exon */

        if( tr->c[y-1] == 0) {
          ret_path->iali[step_cnt]   = tr->i[y-1] - 2;
          ret_path->ihmm[step_cnt]   = tr->k[y-1];
        }
        else if( tr->c[y-1] == 1) {
          ret_path->iali[step_cnt]   = tr->i[y-1] - 1;
          ret_path->ihmm[step_cnt]   = tr->k[y-1];
        }
        else {
          ret_path->iali[step_cnt]   = tr->i[y-1];
          ret_path->ihmm[step_cnt]   = tr->k[y];
          ret_path->jhmm[step_cnt-1] = tr->k[y-1];
        }

        tmp_path->iali[step_cnt] = tr->i[y] - tr->c[y] + 1;
        tmp_path->ihmm[step_cnt] = tr->k[y];
      }
     
      ret_path->jhmm[step_cnt] = tr->k[z];

      if(step_cnt == ret_path->path_len-1) {
        ret_path->jali[step_cnt]   = tr->i[z];
        tmp_path->jali[step_cnt]   = tr->i[z];
      }
      else {
        /* Determine which splice codon type we have */
        if( tr->c[z+1] == 0)
          ret_path->jali[step_cnt] = tr->i[z];
        else if( tr->c[z+1] == 1)
          ret_path->jali[step_cnt] = tr->i[z] + 1;
        else if( tr->c[z+1] == 2)
          ret_path->jali[step_cnt] = tr->i[z] + 2;

        tmp_path->jali[step_cnt] = tr->i[z];
   
      }

      step_cnt++;

      start_new = FALSE;
    }
    
    z++;
    if(tr->st[z] == p7T_M) start_new = TRUE;

  }
  
  ret_path->revcomp = spliced_path->revcomp;
  tmp_path->revcomp = spliced_path->revcomp;  

  /* Convert to true coordinates */
  for(s = 0; s < ret_path->path_len; s++) {
    if(spliced_path->revcomp) {
      ret_path->iali[s] = path_seq->n - ret_path->iali[s] + path_seq->end;
      ret_path->jali[s] = path_seq->n - ret_path->jali[s] + path_seq->end;

      tmp_path->iali[s] = path_seq->n - tmp_path->iali[s] + path_seq->end;
      tmp_path->jali[s] = path_seq->n - tmp_path->jali[s] + path_seq->end;
    }
    else {
      ret_path->iali[s] = path_seq->start + ret_path->iali[s] - 1;
      ret_path->jali[s] = path_seq->start + ret_path->jali[s] - 1;

      tmp_path->iali[s] = path_seq->start + tmp_path->iali[s] - 1;
      tmp_path->jali[s] = path_seq->start + tmp_path->jali[s] - 1;
    }
  }

  /* Assign node ids */
  ret_path->node_id[ret_path->path_len-1] = s_start;
  tmp_path->node_id[tmp_path->path_len-1] = s_start;

  for(s1 = 0; s1 < tmp_path->path_len-1; s1++) {  
    for(s2 = s_start-1; s2 >= 0; s2--) {
      if(p7_splicegraph_NodeOverlap(graph, spliced_path->node_id[s2], tmp_path, s1)) {
        ret_path->node_id[s1] = spliced_path->node_id[s2];
        tmp_path->node_id[s1] = spliced_path->node_id[s2];
      } 
    }
  }

  /* Now that we have node assignments we
   * 1. can add nodes with no node assignment
   * 2. add edges that don't currently exist
   * 3. record start and next start coordinates from tmp_path onto edges
   * 4. record splice sites from ret_path onto edges
   */
  for(s = tmp_path->path_len-2; s >= 0; s--) {
  
    if(tmp_path->node_id[s] == -1) {
      hit          = p7_hit_Create_empty();
      hit->dcl     = p7_domain_Create_empty();
      hit->dcl->tr = p7_trace_fs_Create();

      hit->dcl->iali = tmp_path->iali[s];
      hit->dcl->jali = tmp_path->jali[s];
      hit->dcl->ihmm = tmp_path->ihmm[s];
      hit->dcl->jhmm = tmp_path->jhmm[s];

      hit->dcl->aliscore = 1.; // any positive score ensures it is included in next path

      p7_splicegraph_AddNode(graph, hit);

      graph->tmp_node[graph->num_nodes-1] = TRUE;

      tmp_path->node_id[s] = graph->num_nodes-1;
      ret_path->node_id[s] = graph->num_nodes-1;

      edge = p7_splicegraph_AddEdge(graph, tmp_path->node_id[s], tmp_path->node_id[s+1]);

      /* 3. add start coords */
      edge->i_start = tmp_path->iali[s];
      edge->k_start = tmp_path->ihmm[s];

      edge->i_end = tmp_path->jali[s];
      edge->k_end = tmp_path->jhmm[s];

      edge->next_i_start = tmp_path->iali[s+1];
      edge->next_k_start = tmp_path->ihmm[s+1];
    
      /* 4. Add splice coords */
      edge->upstream_nuc_end       = ret_path->jali[s];
      edge->upstream_amino_end     = ret_path->jhmm[s];
      edge->downstream_nuc_start   = ret_path->iali[s+1];
      edge->downstream_amino_start = ret_path->ihmm[s+1];
    }
    else {

      edge = p7_splicegraph_GetEdge(graph, tmp_path->node_id[s], tmp_path->node_id[s+1]);

      /* 3. add start coords */
      edge->i_start = tmp_path->iali[s];
      edge->k_start = tmp_path->ihmm[s];

      edge->i_end = tmp_path->jali[s];
      edge->k_end = tmp_path->jhmm[s];

      edge->next_i_start = tmp_path->iali[s+1];
      edge->next_k_start = tmp_path->ihmm[s+1];

      /* 4. Add splice coords */
      edge->upstream_nuc_end       = ret_path->jali[s];
      edge->upstream_amino_end     = ret_path->jhmm[s];
      edge->downstream_nuc_start   = ret_path->iali[s+1];
      edge->downstream_amino_start = ret_path->ihmm[s+1];

    }
  
  }

  p7_trace_fs_Destroy(tr);
  p7_splicepath_Destroy(tmp_path);

  return ret_path;

}



/*  Function: p7_splice_AlignSingle
 *  Synopsis: Use global spliced Viterbi to any internal introns in a single hit
 *
 *  Purpose : Align the single anchor node in a path to find splice sites 
 *            (if any) and any new exons not in the path. Add new exons to 
 *            the graph, and record spliced and unspliced coordinates.
 *
 * Returns:   <SPLICE_PATH> on success and NULL on failure.
 *
 */
SPLICE_PATH*
p7_splice_AlignSingle(SPLICE_WORKER_INFO *info, SPLICE_PATH *spliced_path, ESL_SQ *path_seq, int i_start, int i_end, int k_start, int k_end)
{
 
  int         L = i_end - i_start + 1;
  int         M = k_end - k_start + 1;
  int         s, y, z;
  int         z1, z2;
  int         intron_cnt;
  int         step_cnt;
  int         start_new;
  float       ali_score;
  P7_TRACE    *tr;
  SPLICE_PATH *ret_path;
  SPLICE_PIPELINE *pli;
  P7_FS_PROFILE *gm_fs;
  ESL_GENCODE *gcode;

  pli   = info->pli;
  gm_fs = info->gm_fs;
  gcode = info->gcode;

  p7_gmx_sp_GrowTo(pli->vit, M, L, L);
  p7_splicepipline_GrowIndex(pli->sig_idx, M, L, ALIGNMENT_MODE);
  p7_fs_ReconfigLength(gm_fs, L);
  
  p7_spliceviterbi_translated_global(pli, path_seq->dsq, gcode, gm_fs, pli->vit, i_start, i_end, k_start, k_end);

  tr = p7_trace_fs_Create();
  p7_splicevitebi_translated_trace(pli, path_seq->dsq, gcode, gm_fs, pli->vit, tr, i_start, i_end, k_start, k_end, &ali_score);
  //p7_trace_fs_Dump(stdout, tr, NULL, NULL, NULL);

  /* Find number of introns in trace */
  intron_cnt = 0;
  for(z = 0; z < tr->N; z++)
    if(tr->st[z] == p7T_P) intron_cnt++;

  if(intron_cnt == 0) {
    p7_trace_fs_Destroy(tr); 
    return NULL;
  }

  /* Find first M state - start of first hit */
  for(z1 = 0; z1 < tr->N; z1++) if(tr->st[z1] == p7T_M) break;
  if (z1 == tr->N) {
    p7_trace_fs_Destroy(tr);
    return NULL;
  }

  /* Find last M state state - end of last hit */
  for(z2 = tr->N-1; z2 >= 0; z2--) if(tr->st[z2] == p7T_M) break;
  if (z2 == -1) {
    p7_trace_fs_Destroy(tr);
    return NULL;
  }

  ret_path = p7_splicepath_Create(intron_cnt+1);

  step_cnt = 0;
  start_new = TRUE;
  
  z = z1;
  while(z <= z2) {

    if(start_new) {
      
      /* Save z value - currently set to first M state in exon */
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
      }
      else {
        /* Determine which splice codon type we have and assign the P state (y-1) amino to the correct exon */

        if( tr->c[y-1] == 0) {
          ret_path->iali[step_cnt]   = tr->i[y-1] - 2;
          ret_path->ihmm[step_cnt]   = tr->k[y-1];
        }
        else if( tr->c[y-1] == 1) {
          ret_path->iali[step_cnt]   = tr->i[y-1] - 1;
          ret_path->ihmm[step_cnt]   = tr->k[y-1];
        }
        else {
          ret_path->iali[step_cnt]   = tr->i[y-1];
          ret_path->ihmm[step_cnt]   = tr->k[y];
          ret_path->jhmm[step_cnt-1] = tr->k[y-1];
        }
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

  ret_path->revcomp = spliced_path->revcomp;

  /* Convert to true coordinates */
  for(s = 0; s < ret_path->path_len; s++) {

    if(spliced_path->revcomp) {
      ret_path->iali[s] = path_seq->n - ret_path->iali[s] + path_seq->end;
      ret_path->jali[s] = path_seq->n - ret_path->jali[s] + path_seq->end;
    }
    else {
      ret_path->iali[s] = path_seq->start + ret_path->iali[s] - 1;
      ret_path->jali[s] = path_seq->start + ret_path->jali[s] - 1;
    }
  } 

  ret_path->node_id[0] = spliced_path->node_id[0];

  p7_trace_fs_Destroy(tr);

  return ret_path;

}


/*  Function: p7_splice_EnforceBounds
 *  Synopsis: Remove edges that cross a sequence boundry 
 *
 *  Purpose : After a path has been ailgned, break all edges that cross the
 *            path's sequnce boudries to prevent interwoven alignments. 
 *
 * Returns:   <eslOK> 
 *
 */
int
p7_splice_EnforceBounds(SPLICE_GRAPH *graph, int64_t bound_min, int64_t bound_max) 
{

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
      } 
    }
  }

  return eslOK;

}

/*****************************************************************
 * 3. Path Alginment
 *****************************************************************/

/*  Function: p7_splice_AlignSplicedPath
 *  Synopsis: Create a spliced alignment from a spliced path
 *
 *  Purpose : Given a fully spliced SPLICE_PATH create the protein sequence produced by 
 *            that splicing, algin it to the model using protein to protein alignment.  
 *            Convert that protein to protein alginment to a spliced alignment and 
 *            replace one of path's anchor nodes' original hits with the new spliced hit. 
 *            Mark the original hits of any ther anchor nodes in the path as unreportable.
 *
 * Returns:   <eslOK>.
 *
 */
int
p7_splice_AlignSplicedPath(SPLICE_WORKER_INFO *info, SPLICE_PATH *orig_path, SPLICE_PATH *spliced_path, ESL_SQ *path_seq, int *success)
{
  int           e, i, s;
  int           s1, s2;
  int           shift;
  int           orf_len;
  int           remove_node;
  int           replace_node;
  int           contains_anchor;
  int           seq_min, seq_max;
  float         dom_bias;
  float         nullsc;
  float         dom_score;
  double        dom_lnP;
  ESL_SQ       *new_seq;
  P7_HIT       *replace_hit;
  P7_HIT       *remove_hit;
  P7_TOPHITS *tophits;
  SPLICE_GRAPH *graph;
  SPLICE_PIPELINE *pli;
  P7_OPROFILE *om;
  ESL_SQFILE *seq_file;
  int          status;


  tophits = info->tophits;
  graph = info->graph;
  pli = info->pli;
  om = info->om;
  seq_file = info->seq_file;
 
  /* Create an amino sequence from the spliced nuc sequece whose coords are in the spliced_path */ 
  p7_splice_CreateSplicedSequnce(info, spliced_path, path_seq);
 
  /* Algin the spliced amino sequence */
  status = p7_splice_AlignSplicedSequence(info, spliced_path, path_seq);

  if(status == eslEINACCURATE) {
    /* path has been altered, redo alignmnt with new path */
    p7_splicepipeline_Reuse(pli);

    /* Check if we need to refetch the path_seq */
    if(( spliced_path->revcomp && (spliced_path->iali[0]                        > path_seq->start - ALIGNMENT_EXT   || 
                                   spliced_path->jali[spliced_path->path_len-1]         < path_seq->end   + ALIGNMENT_EXT)) ||
       (!spliced_path->revcomp && (spliced_path->iali[0]                        < path_seq->start + ALIGNMENT_EXT   || 
                                   spliced_path->jali[spliced_path->path_len-1] > path_seq->end   - ALIGNMENT_EXT))) {
  
      seq_min = ESL_MIN(spliced_path->iali[0], spliced_path->jali[spliced_path->path_len-1]) - ALIGNMENT_EXT;
      seq_max = ESL_MAX(spliced_path->iali[0], spliced_path->jali[spliced_path->path_len-1]) + ALIGNMENT_EXT;
      new_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, seq_min, seq_max, spliced_path->revcomp, info);
    
      p7_splice_AlignSplicedPath(info, orig_path, spliced_path, new_seq, success);
      esl_sq_Destroy(new_seq);
    }
    else 
      p7_splice_AlignSplicedPath(info, orig_path, spliced_path, path_seq, success);

    return eslOK;
    
  } 
  else if(pli->hit == NULL || pli->hit->dcl->ad->exon_cnt == 1) return eslOK;

  /* adjust the hit coords */
  if(spliced_path->revcomp) {
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
  dom_score -= 2 * log(2. / (pli->amino_sq->n+2));
  dom_score += 2 * log(2. / (om->max_length+2));
  dom_score -= (pli->amino_sq->n-orf_len)      * log((float) (pli->amino_sq->n) / (float) (pli->amino_sq->n+2));
  dom_score += (om->max_length-orf_len) * log((float) om->max_length / (float) (om->max_length+2));
 
  /* Bias calculation and adjustments */
  if(pli->do_null2) 
    dom_bias = p7_FLogsum(0.0, log(pli->bg->omega) + pli->hit->dcl->domcorrection);
  else 
    dom_bias = 0.;

  p7_bg_SetLength(pli->bg, om->max_length);
  p7_bg_NullOne  (pli->bg, pli->amino_sq->dsq, om->max_length, &nullsc);
  dom_score = (dom_score - (nullsc + dom_bias))  / eslCONST_LOG2;

  dom_lnP   = esl_exp_logsurv(dom_score, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

  /* E-value adjusment */
  dom_lnP += log((float)info->db_nuc_cnt / (float)om->max_length);

  if ((pli->by_E && exp(dom_lnP) <= pli->E) || ((!pli->by_E) && dom_score >= pli->T)) {

    /*If the first or last exon does not pass the reporting threshold, 
     * break the edge to the corresponding node and realign */
    for(e = 0; e < pli->hit->dcl->ad->exon_cnt; e++) 
      pli->hit->dcl->ad->exon_lnP[e] += log((float)info->db_nuc_cnt / (float)om->max_length);

    
    if ( spliced_path->path_len >  pli->hit->dcl->ad->exon_cnt) {
      /* Shift the path to start at the first hit that was inculded in the alignment
       * and end at the last hit that was included in the alignment */
      if(spliced_path->revcomp) {
        for(shift = 0; shift < spliced_path->path_len; shift++) {
          if(spliced_path->jali[shift] <= pli->hit->dcl->iali) break;
        }
      }
      else {
        for(shift = 0; shift < spliced_path->path_len; shift++) {
           if(spliced_path->jali[shift] >= pli->hit->dcl->iali) break;
        }
      }

      /* Shift path to start at frist hits that is in alignment */
      for(s = 0; s < shift; s++) 
        p7_splicepath_Remove(spliced_path, 0);

      spliced_path->iali[0] = pli->hit->dcl->iali;
      spliced_path->ihmm[0] = pli->hit->dcl->ihmm;

      spliced_path->path_len =  pli->hit->dcl->ad->exon_cnt;

      spliced_path->jali[spliced_path->path_len-1] = pli->hit->dcl->jali;
      spliced_path->jhmm[spliced_path->path_len-1] = pli->hit->dcl->jhmm;

      /*Make sure path still contains an anchor */
      contains_anchor = FALSE;
      for(s = 0; s < spliced_path->path_len; s++) 
        if(spliced_path->node_id[s] < graph->anchor_N && spliced_path->node_id[s] >= 0) contains_anchor = TRUE;

      if(!contains_anchor) return eslOK;       
  
    }

    /* Find the first original hit in path to copy info*/
    i = 0;
    while( spliced_path->node_id[i] >= graph->anchor_N || spliced_path->node_id[i] == -1) {
      pli->hit->dcl->ad->exon_orig[i] = FALSE;
      pli->hit->dcl->ad->exon_split[i] = FALSE;
      i++;
    }
    pli->hit->dcl->ad->exon_orig[i] = TRUE;
    if(i < spliced_path->path_len-2 && spliced_path->node_id[i] == spliced_path->node_id[i+1])
      pli->hit->dcl->ad->exon_split[i] = TRUE;
    else 
      pli->hit->dcl->ad->exon_split[i] = FALSE;


#ifdef HMMER_THREADS
  if(info->thread_id >= 0) pthread_mutex_lock(info->mutex);
#endif /*HMMER_THREADS*/

    replace_node = spliced_path->node_id[i];
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
    for(  ; i < spliced_path->path_len; i++) {
      remove_node = spliced_path->node_id[i];
      if(remove_node < 0 || remove_node >= graph->anchor_N) {
        pli->hit->dcl->ad->exon_orig[i] = FALSE;
        pli->hit->dcl->ad->exon_split[i] = FALSE; 
      }
      else {
        pli->hit->dcl->ad->exon_orig[i] = TRUE;
        if((i > 0                && spliced_path->node_id[i] == spliced_path->node_id[i-1]) ||
           (i < spliced_path->path_len-2 && spliced_path->node_id[i] == spliced_path->node_id[i+1]))
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
      if(graph->th->hit[orig_path->node_id[s1]]->dcl->ad != NULL) continue; // skip hits that have been replaced with spliced hits 
      if(graph->th->hit[orig_path->node_id[s1]]->flags & p7_IS_REPORTED ) {
        for(s2 = 0; s2 < spliced_path->path_len; s2++) {
          if(spliced_path->node_id[s2] >= graph->anchor_N) continue;
          if(p7_splicegraph_NodeOverlap(graph, orig_path->node_id[s1], spliced_path, s2)) {
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

  return eslOK;
   
}


/*  Function: p7_splice_CreateSplicedSequnce
 *  Synopsis: Create a protein sequnce from a spliced path
 *
 *  Purpose : Given a fully spliced SPLICE_PATH create the nucleotide sequence produced by
 *            that splicing, keeping track of the coodinates on the path_seq, and then
 *            translate that into a protein sequence.
 *
 * Returns:   <eslOK> on success.
 *
 *  Throws:    <eslEMEM> on allocation failure.
 */

int
p7_splice_CreateSplicedSequnce(SPLICE_WORKER_INFO *info, SPLICE_PATH *spliced_path, ESL_SQ *path_seq) 
{

  int      i, s;
  int      path_seq_len;
  int      pos, seq_pos;
  int      seq_idx;
  int      path_start_pos;
  int      path_end_pos;
  int      ext_start_pos;
  int      ext_end_pos;
  int      amino_len;
  ESL_DSQ  amino;
  int64_t  *nuc_index;
  ESL_DSQ  *nuc_dsq;
  ESL_DSQ  *amino_dsq;
  ESL_SQ   *nuc_seq;
  ESL_SQ   *amino_seq;
  ESL_GENCODE *gcode;
  SPLICE_PIPELINE *pli;
  int status;  

  gcode = info->gcode;
  pli   = info->pli;

  nuc_index = NULL;
  nuc_dsq   = NULL;
  amino_dsq = NULL;

  path_seq_len = 0;
  for (s = 0; s < spliced_path->path_len; s++) 
    path_seq_len += llabs(spliced_path->jali[s] - spliced_path->iali[s]) + 1;
  
  /* If the spliced seqeunce length is non-mod 3 send to farmshift alignment */
  if(path_seq_len % 3 != 0) {
    spliced_path->frameshift = TRUE;
    p7_splicepath_Dump(stdout, spliced_path);
    ESL_XEXCEPTION(eslFAIL, "NOT mod 3");
    return eslOK;
  }


//TODO - Test if adding extantions improves results
  /* Working backward from the start of the path find the first instance of a stop codon in the extended upstream region of path_seq */
  if (spliced_path->revcomp) {
    path_start_pos = path_seq->n - spliced_path->iali[0] + path_seq->end;
    ext_start_pos  = path_seq->n - (spliced_path->iali[0] + ALIGNMENT_EXT) + path_seq->end;
    for(pos = spliced_path->iali[0]+3; pos <= spliced_path->iali[0] + ALIGNMENT_EXT; pos+=3) {

      seq_pos = path_seq->n - pos + path_seq->end;   
      if(seq_pos < 1) {
        ext_start_pos = seq_pos+3;
        break;
      }
      
      amino = esl_gencode_GetTranslation(gcode,&path_seq->dsq[seq_pos]);
      
      if(esl_abc_XIsNonresidue(gcode->aa_abc, amino)) { 
        ext_start_pos = seq_pos+3;
        break;
      }
    }
  }
  else {
    path_start_pos = spliced_path->iali[0] - path_seq->start + 1;
    ext_start_pos  = (spliced_path->iali[0] - ALIGNMENT_EXT) - path_seq->start + 1;
  
    for(pos = spliced_path->iali[0]-3; pos >= spliced_path->iali[0] - ALIGNMENT_EXT; pos-=3) {
       
      seq_pos = pos - path_seq->start + 1;
      
      if(seq_pos < 1) {
        ext_start_pos = seq_pos+3;
        break;
      }
      amino = esl_gencode_GetTranslation(gcode,&path_seq->dsq[seq_pos]);
 
      if(esl_abc_XIsNonresidue(gcode->aa_abc, amino)) {
        ext_start_pos = seq_pos+3;
        break;
      }
    }
  }

  /* Working forward from the end of the path find the first instance of a stop codon in the extended downstream region of path_seq */ 
  if (spliced_path->revcomp) {
    path_end_pos =  path_seq->n - spliced_path->jali[spliced_path->path_len-1] + path_seq->end;
    ext_end_pos  =  path_seq->n - (spliced_path->jali[spliced_path->path_len-1] - ALIGNMENT_EXT) + path_seq->end; 
    for(pos = spliced_path->jali[spliced_path->path_len-1]-1; pos >= spliced_path->jali[spliced_path->path_len-1] - ALIGNMENT_EXT; pos-=3) {
      seq_pos = path_seq->n - pos + path_seq->end;
      if(seq_pos > path_seq->n-2) {
        ext_end_pos = seq_pos-1;
        break;
      }
      
      amino = esl_gencode_GetTranslation(gcode,&path_seq->dsq[seq_pos]);
 
      if(esl_abc_XIsNonresidue(gcode->aa_abc, amino)) {
        ext_end_pos = seq_pos-1;
        break;
      } 
 
    }
  }
  else {
    path_end_pos = spliced_path->jali[spliced_path->path_len-1] - path_seq->start + 1;
    ext_end_pos  = (spliced_path->jali[spliced_path->path_len-1] +  ALIGNMENT_EXT) - path_seq->start + 1;
    
    for(pos = spliced_path->jali[spliced_path->path_len-1]+1; pos <= spliced_path->jali[spliced_path->path_len-1] +  ALIGNMENT_EXT; pos+=3) {

      seq_pos = pos - path_seq->start + 1;
      if(seq_pos > path_seq->n-2 ) {
        ext_end_pos = seq_pos-1;
        break;
      }

      amino = esl_gencode_GetTranslation(gcode,&path_seq->dsq[seq_pos]);

      if(esl_abc_XIsNonresidue(gcode->aa_abc, amino)) {
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
  for (s = 0; s < spliced_path->path_len; s++) {
    if (spliced_path->revcomp) {
      for (pos = spliced_path->iali[s]; pos >= spliced_path->jali[s]; pos--) {
        seq_pos = path_seq->n - pos + path_seq->end;
        nuc_index[seq_idx] = seq_pos;
        nuc_dsq[seq_idx]   = path_seq->dsq[seq_pos];
        seq_idx++;
      }
    }
    else {
      for (pos = spliced_path->iali[s]; pos <= spliced_path->jali[s]; pos++) {

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
    
    if(esl_abc_XIsNonresidue(gcode->aa_abc, amino)) {
      spliced_path->frameshift = TRUE;
  
      p7_splicepath_Dump(stdout, spliced_path); 
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

  if(nuc_dsq   != NULL) free(nuc_dsq);
  if(amino_dsq != NULL) free(amino_dsq);

  return eslOK;

  ERROR:
    if(nuc_index != NULL) free(nuc_index);
    if(nuc_dsq   != NULL) free(nuc_dsq);
    if(amino_dsq != NULL) free(amino_dsq);
    return status;

}


/*  Function: p7_splice_AlignSplicedSequnce
 *  Synopsis: Align a spliced protein sequence to a protein model
 *
 *  Purpose : Given a spliced protein sequence produced and a model, align them using Fwd/Bwd.
 *            Address underflow errors caused by low scoreing exons by removing them from 
 *            the path. Produce a spliced alignment and provide scores, p-values, and average 
 *            posterior probablities for all exons in the alignment. 
 *
 * Returns:   <eslOK> on success and <eslEINACCURATE> if an underflow error required altering the path.
 *
 */
int
p7_splice_AlignSplicedSequence(SPLICE_WORKER_INFO *info, SPLICE_PATH *spliced_path, ESL_SQ *path_seq) 
{

  int       i, e;
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
  SPLICE_GRAPH *graph;
  SPLICE_PIPELINE *pli;
  P7_OPROFILE *om;
  P7_PROFILE *gm;
  int       status;

  graph = info->graph;
  pli   = info->pli;
  om    = info->om;
  gm    = info->gm;

  tr           = p7_trace_CreateWithPP();
  hit          = p7_hit_Create_empty();
  hit->dcl     = p7_domain_Create_empty();
  hit->dcl->tr = NULL;
  hit->dcl->scores_per_pos = NULL;

  p7_oprofile_ReconfigUnihit(om, pli->amino_sq->n);

  p7_omx_GrowTo(pli->fwd, om->M, pli->amino_sq->n, pli->amino_sq->n);
  p7_omx_GrowTo(pli->bwd, om->M, pli->amino_sq->n, pli->amino_sq->n);
  p7_omx_GrowTo(pli->pp,  om->M, pli->amino_sq->n, pli->amino_sq->n);

  p7_bg_SetLength(pli->bg, pli->amino_sq->n);
  if (pli->do_biasfilter)
    p7_bg_FilterScore(pli->bg, pli->amino_sq->dsq, pli->amino_sq->n, &filtersc);
  else
    p7_bg_NullOne  (pli->bg, pli->amino_sq->dsq, pli->amino_sq->n, &filtersc);

  p7_Forward (pli->amino_sq->dsq, pli->amino_sq->n, om, pli->fwd, &envsc);
  p7_Backward(pli->amino_sq->dsq, pli->amino_sq->n, om, pli->fwd, pli->bwd, NULL);

  if((status = p7_Decoding(om, pli->fwd, pli->bwd, pli->pp)) == eslERANGE) { 
    /* This is a rare event usually caused by a low probability exon somewhere in the path. 
     * If we can find the offending exon and cut the path in two at that point then we can 
     * save the good exons, but to do that we need an alignment so we createone with Viterbi */
    
    p7_gmx_fs_GrowTo(pli->gfwd, gm->M, pli->amino_sq->n, pli->amino_sq->n, p7P_CODONS);
    p7_ReconfigUnihit(gm, pli->amino_sq->n); 

    p7_GViterbi(pli->amino_sq->dsq, pli->amino_sq->n, gm, pli->gfwd, NULL);
    p7_GTrace(pli->amino_sq->dsq, pli->amino_sq->n, gm, pli->gfwd, tr);
   
    p7_trace_Index(tr);
    hit->dcl->tr = p7_trace_splice_Convert(tr, pli->orig_nuc_idx, &splice_cnt);

    if(splice_cnt == 0) {
      p7_hit_Destroy(hit);
      p7_trace_Destroy(tr);
      return eslOK;
    }
    
    hit->dcl->ad = p7_alidisplay_splice_Create(hit->dcl->tr, 0, om, path_seq, pli->amino_sq, hit->dcl->scores_per_pos, tr->sqfrom[0], splice_cnt);

    p7_splice_ScoreExons(pli, tr, hit->dcl->ad, om, NULL, FALSE, FALSE); 
    status = p7_splice_FixDecodingErrors(graph, spliced_path, hit->dcl->ad, path_seq);

    p7_trace_splice_Destroy(hit->dcl->tr);
    p7_alidisplay_Destroy(hit->dcl->ad);
  
    hit->dcl->tr = NULL;
    hit->dcl->ad = NULL;

    p7_hit_Destroy(hit);
    p7_trace_Destroy(tr);

    return status;
  }
  
  p7_OptimalAccuracy(om, pli->pp, pli->bwd, &oasc); /* <bwd> is now overwritten with OA scores */
  p7_OATrace        (om, pli->pp, pli->bwd, tr);

  p7_trace_Index(tr);
  hit->dcl->tr = p7_trace_splice_Convert(tr, pli->orig_nuc_idx, &splice_cnt);
   
  if(splice_cnt == 0) {
    p7_trace_splice_Destroy(hit->dcl->tr);
    hit->dcl->tr = NULL;
    p7_hit_Destroy(hit);
    p7_trace_Destroy(tr);
    return eslOK;
  } 

  
  seq_score = (envsc-filtersc) / eslCONST_LOG2;
  P = esl_exp_surv(seq_score, om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

  if (P > pli->F3) {
    p7_trace_splice_Destroy(hit->dcl->tr);
    hit->dcl->tr = NULL;
    p7_hit_Destroy(hit);
    p7_trace_Destroy(tr);
    return eslOK;
  }

  p7_Null2_ByExpectation(om, pli->pp, null2);
  domcorrection = 0.;
  for (i = 1; i <= pli->amino_sq->n; i++) {
    domcorrection += logf(null2[pli->amino_sq->dsq[i]]);
  }

  hit->dcl->domcorrection = ESL_MAX(0.0, domcorrection);


  hit->dcl->ad = p7_alidisplay_splice_Create(hit->dcl->tr, 0, om, path_seq, pli->amino_sq, hit->dcl->scores_per_pos, tr->sqfrom[0], splice_cnt);

  p7_splice_ScoreExons(pli, tr, hit->dcl->ad, om, null2, pli->do_null2, TRUE);

  /*Check for 0 posterior probablity caused by low quailty exons */
  for(e = 0; e < hit->dcl->ad->exon_cnt; e++) {

    if(hit->dcl->ad->exon_pp[e] == 0.) {
      status = p7_splice_FixDecodingErrors(graph, spliced_path, hit->dcl->ad, path_seq);
     
      p7_trace_splice_Destroy(hit->dcl->tr);
      p7_alidisplay_Destroy(hit->dcl->ad);

      hit->dcl->tr = NULL;
      hit->dcl->ad = NULL;

      p7_hit_Destroy(hit);
      p7_trace_Destroy(tr);

      return status; 

    }
  }
  
  hit->dcl->ihmm = hit->dcl->ad->hmmfrom;
  hit->dcl->jhmm = hit->dcl->ad->hmmto;

  /* Convert sqfrom, sqto to full sequence coords */
  if(spliced_path->revcomp) {
    hit->dcl->ad->sqto    = path_seq->n - hit->dcl->ad->sqto   + path_seq->end;
    hit->dcl->ad->sqfrom  = path_seq->n - hit->dcl->ad->sqfrom + path_seq->end;
  }
  else {
    hit->dcl->ad->sqfrom  = hit->dcl->ad->sqfrom + path_seq->start - 1;
    hit->dcl->ad->sqto    = hit->dcl->ad->sqto   + path_seq->start - 1;
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

/*  Function: p7_splice_AlignSplicedSequnce
 *  Synopsis: If a plsiced alignment experienced an underflow error, find and remove the offending exon 
 *
 *  Purpose : Low scoring exon can cause underflow errors in the posterior probaility matrix. 
 *            Here we find the most likely culprit and split the path at that point.  
 *            The remaining path can then be realigned. If the remaining path has only one 
 *            exon or contains no anchor nodes, the path is rejected. 
 *
 * Returns:   <eslEINACCURATE> if the remaining path needs to be realgined and <eslOK> if not. 
 *
 */
int
p7_splice_FixDecodingErrors(SPLICE_GRAPH *graph, SPLICE_PATH *spliced_path, P7_ALIDISPLAY *ad, ESL_SQ *path_seq)
{

  int e, s;
  int shift;
  int contains_anchor;
  int min_idx;
  float min_score;
 
  /* If the alignmnt has rejected some exons we first remove those from the path */ 
  if(spliced_path->path_len > ad->exon_cnt) {

    /* Convert sqfrom, sqto to full sequence coords */
    if(spliced_path->revcomp) {
      ad->sqto    = path_seq->n - ad->sqto   + path_seq->end;
      ad->sqfrom  = path_seq->n - ad->sqfrom + path_seq->end;
    }
    else {
      ad->sqfrom  = ad->sqfrom + path_seq->start - 1;
      ad->sqto    = ad->sqto   + path_seq->start - 1;
    }

    if(spliced_path->revcomp) {
      for(shift = 0; shift < spliced_path->path_len; shift++)
        if(spliced_path->jali[shift] <= ad->sqfrom) break;
    }
    else {
      for(shift = 0; shift < spliced_path->path_len; shift++)
         if(spliced_path->jali[shift] >= ad->sqfrom) break;
    }

    /* Shift path to start at frist hits that is in alignment */
    for(s = 0; s < shift; s++)
      p7_splicepath_Remove(spliced_path, 0);

    spliced_path->iali[0] = ad->sqfrom;
    spliced_path->ihmm[0] = ad->hmmfrom;

    spliced_path->path_len =  ad->exon_cnt;

    spliced_path->jali[spliced_path->path_len-1] = ad->sqto;
    spliced_path->jhmm[spliced_path->path_len-1] = ad->hmmto;

    /*Remove any starting or ending non-anchor nodes */
    if(spliced_path->path_len == 1) return eslOK;

    while(spliced_path->node_id[0] < 0 || spliced_path->node_id[0] >= graph->anchor_N) {
      p7_splicepath_Remove(spliced_path, 0);
      if(spliced_path->path_len == 1) return eslOK;
    }
    spliced_path->iali[0] = graph->th->hit[spliced_path->node_id[0]]->dcl->iali;
    spliced_path->ihmm[0] = graph->th->hit[spliced_path->node_id[0]]->dcl->ihmm;

    while(spliced_path->node_id[spliced_path->path_len-1] < 0 || spliced_path->node_id[spliced_path->path_len-1] >= graph->anchor_N) { 
      spliced_path->path_len--;
      spliced_path->jali[spliced_path->path_len-1] = graph->th->hit[spliced_path->node_id[spliced_path->path_len-1]]->dcl->jali;
      spliced_path->jhmm[spliced_path->path_len-1] = graph->th->hit[spliced_path->node_id[spliced_path->path_len-1]]->dcl->jhmm;
      if(spliced_path->path_len == 1) return eslOK;
    }

    
  }
  else {
    /* Use the Exon scores to find the weakest place in the path */

    min_idx   = 0;
    min_score = ad->exon_score[0];
    for(e = 0; e < ad->exon_cnt; e++) {  
      if(isnan(ad->exon_score[e]) || ad->exon_score[e] == -eslINFINITY) {
        min_idx = e;
        break;
      }
      if( ad->exon_score[e] < min_score) {
          min_score = ad->exon_score[e];
          min_idx = e;
      }  
    }

    if(min_idx == 0) {
      p7_splicepath_Remove(spliced_path, 0);
      if(spliced_path->path_len == 1) return eslOK;

      while(spliced_path->node_id[0] < 0 || graph->tmp_node[spliced_path->node_id[0]]) {
        p7_splicepath_Remove(spliced_path, 0);
        if(spliced_path->path_len == 1) return eslOK;
      }
      spliced_path->iali[0] = graph->th->hit[spliced_path->node_id[0]]->dcl->iali;
      spliced_path->ihmm[0] = graph->th->hit[spliced_path->node_id[0]]->dcl->ihmm;
    }
    else if (min_idx > 0) {
      spliced_path->path_len = min_idx;
      if(spliced_path->path_len == 1) return eslOK; 
 
      while(spliced_path->node_id[spliced_path->path_len-1] < 0 ||  graph->tmp_node[spliced_path->node_id[spliced_path->path_len-1]]) {
        spliced_path->path_len--; 
        if(spliced_path->path_len == 1) return eslOK;
      }
      spliced_path->jali[spliced_path->path_len-1] = graph->th->hit[spliced_path->node_id[spliced_path->path_len-1]]->dcl->jali;
      spliced_path->jhmm[spliced_path->path_len-1] = graph->th->hit[spliced_path->node_id[spliced_path->path_len-1]]->dcl->jhmm;
    }
  
  }
 
  /*Make sure path still contains an anchor */
  contains_anchor = FALSE;
  for(s = 0; s < spliced_path->path_len; s++)
    if(spliced_path->node_id[s] < graph->anchor_N && spliced_path->node_id[s] >= 0) contains_anchor = TRUE;

  if(!contains_anchor) return eslOK;
 
  return eslEINACCURATE;
   
}


/*  Function: p7_splice_ScoreExons
 *  Synopsis: Produce Forward scores, p-values and average posterior probailities for each exon, 
 *
 *  Purpose : Extract approximate Forward scores for each exon from the fwd matrix. 
 *            Use these scores to produce p-values. Also find the averapge posertior 
 *            probablity of all emited aminos in each exon's alignment.
 *
 * Returns:   <eslOK>.
 *
 */
int
p7_splice_ScoreExons(SPLICE_PIPELINE *pli, P7_TRACE *tr, P7_ALIDISPLAY *ad, P7_OPROFILE *om, float *null2, int do_null2, int do_pp)
{

  int i, e, z;
  int start_i, end_i;
  int remainder;
  int exon_nuc_len;
  int exon_amino_len; 
  float scale;
  float start_score;
  float end_score;
  float exon_score;
  float exon_pp;
  float nullsc;
  P7_OMX* fwd;
  P7_BG *bg;
  ESL_SQ *amino_sq;

  fwd      = pli->fwd;
  bg       = pli->bg;
  amino_sq = pli->amino_sq;

//p7_trace_Dump(stdout, tr, NULL, NULL);

  p7_bg_SetLength(bg, om->max_length);
  p7_bg_NullOne  (bg, amino_sq->dsq, om->max_length, &nullsc);

  for(z = 0; z < tr->N; z++) if(tr->st[z] == p7T_M) break;
  start_i = tr->i[z]-1; // start_i is the -1 from the first amino position in the exon 
  
  scale = 0.;
  for(i = 0; i <= start_i; i++)
    scale += log(fwd->xmx[i*p7X_NXCELLS+p7X_SCALE]);

  if(start_i == 0) start_score = 0.;
  else             start_score = log(fwd->xmx[start_i*p7X_NXCELLS+p7X_C]) + scale;
  exon_nuc_len = llabs(ad->exon_seq_ends[0] - ad->exon_seq_starts[0]) + 1;

  remainder = exon_nuc_len % 3;

  /* If a codon is split at tke splice site assign the full codon to the exon that contributed two nucleotides */
  if(remainder == 1) exon_nuc_len--;
  if(remainder == 2) exon_nuc_len++;

  exon_amino_len = exon_nuc_len / 3;
  
  end_i = start_i + exon_amino_len;
  
  ad->exon_bias[0] = 0.;
  for(i = start_i+1; i <= end_i; i++) 
    scale += log(fwd->xmx[i*p7X_NXCELLS+p7X_SCALE]); 
  
  if(do_null2) {
    ad->exon_bias[0] = 0.;
    for(i = start_i+1; i <= end_i; i++) 
      ad->exon_bias[0] += logf(null2[amino_sq->dsq[i]]);
    
    ad->exon_bias[0] = ESL_MAX(0., p7_FLogsum(0.0, log(bg->omega) + ad->exon_bias[0]));
  }
  else ad->exon_bias[0] = 0.;

  end_score = log(fwd->xmx[end_i*p7X_NXCELLS+p7X_C]) + scale;

  exon_score = (end_score - start_score);

  /* Subtract out old N->B perobablity and add in new N->B and C->T probabilities */
  exon_score -=     log(2.0 / ((float) amino_sq->n + 2.0));
  exon_score += 2 * log(2.0 / ((float) om->max_length + 2.0));
  ad->exon_score[0] = (exon_score - (nullsc + ad->exon_bias[0]))  / eslCONST_LOG2;

//  printf("start_score %f end_score %f exon_score %f nullsc %f ad->exon_bias[0] %f ad->exon_score[0] %f\n", start_score, end_score, exon_score, nullsc,  ad->exon_bias[0], ad->exon_score[0]); 
  ad->exon_lnP[0] = esl_exp_logsurv (ad->exon_score[0],  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

  if(do_pp) {
    exon_pp = tr->pp[z];
    z++;
    while(tr->i[z] <= end_i)  { exon_pp += tr->pp[z]; z++; }
    
    ad->exon_pp[0] = exon_pp / (float) exon_amino_len;

  }
  else ad->exon_pp[0] = -eslINFINITY;
 //printf("ad->exon_score[e] %f ad->exon_pp[e] %f exp(ad->exon_lnP[e]) %f\n", ad->exon_score[0],  ad->exon_pp[0], exp(ad->exon_lnP[0]));

  for(e = 1; e < ad->exon_cnt; e++) {

    start_i = end_i;
    start_score = end_score;

    exon_nuc_len = llabs(ad->exon_seq_ends[e] - ad->exon_seq_starts[e]) + 1;
    /* Handle posssible upstream codon split */
    if(remainder == 1) exon_nuc_len++;
    if(remainder == 2) exon_nuc_len--;

    /* Handle posssible downstream codon split */
    remainder = exon_nuc_len % 3;
    if(remainder == 1) exon_nuc_len--;
    if(remainder == 2) exon_nuc_len++;
   
    exon_amino_len = exon_nuc_len / 3;
 
    end_i = start_i + exon_amino_len;
   
    ad->exon_bias[e] = 0.;   
   
    for(i = start_i+1; i <= end_i; i++) { 
      scale += log(fwd->xmx[i*p7X_NXCELLS+p7X_SCALE]);

    }
    if(do_null2) {
      for(i = start_i+1; i <= end_i; i++)
        ad->exon_bias[e] += logf(null2[amino_sq->dsq[i]]);
      ad->exon_bias[e] = ESL_MAX(0., p7_FLogsum(0.0, log(bg->omega) + ad->exon_bias[e]));
    }

    ad->exon_bias[e] = ESL_MAX(ad->exon_bias[e], 0.);

    end_score = log(fwd->xmx[end_i*p7X_NXCELLS+p7X_C]) + scale;

    exon_score = (end_score - start_score);

    /* Subtract out old N->B perobablity and add in new N->B and C->T probabilities */
    exon_score -=     log(2.0 / ((float) amino_sq->n + 2.0));
    exon_score += 2 * log(2.0 / ((float) om->max_length + 2.0));

    ad->exon_score[e] = (exon_score - (nullsc + ad->exon_bias[e]))  / eslCONST_LOG2;  

    ad->exon_lnP[e] = esl_exp_logsurv (ad->exon_score[e],  om->evparam[p7_FTAU], om->evparam[p7_FLAMBDA]);

    if(do_pp) {
      exon_pp = tr->pp[z];
      z++;
     
      while(tr->i[z] <= end_i && tr->st[z] != p7T_E) { exon_pp += tr->pp[z]; z++; }

      ad->exon_pp[e] = exon_pp / (float) exon_amino_len;

    }
    else ad->exon_pp[e] = -eslINFINITY;
 //printf("start_score %f end_score %f exon_score %f nullsc %f ad->exon_bias[0] %f ad->exon_score[0] %f\n", start_score, end_score, exon_score, nullsc,  ad->exon_bias[e], ad->exon_score[e]);
   // printf("ad->exon_score[e] %f ad->exon_pp[e] %f exp(ad->exon_lnP[e]) %f\n", ad->exon_score[e],  ad->exon_pp[e], exp(ad->exon_lnP[e]));
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





