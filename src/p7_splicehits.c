/* SPLICE_SAVED_HITS: list of SPLICE_HIT_INFO for a set of hits
 *
 * Contents:
 *    1. The SPLICE_SAVED_HITS object.
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

static int hit_sorter(const void *vh1, const void *vh2);

/*****************************************************************
 * 1. The SPLICE_SAVED_HITS structure.
 *****************************************************************/

/* Function:  p7_splichits_CreateSavedHits()
 * Synopsis:  Allocates a saved hits list.
 *
 * Purpose:   Allocates a new <SPLICE_SAVED_HITS> 
 *
 * Returns:   a pointer to the new <SPLICE_SAVED_HITS> structure 
 *            on success.
 *
 * Throws:    <NULL> on allocation error.
 */
SPLICE_SAVED_HITS*
p7_splicehits_CreateSavedHits(void)
{
  int default_nalloc = 256;
  int status;

  SPLICE_SAVED_HITS *saved_hits = NULL;
  ESL_ALLOC(saved_hits, sizeof(SPLICE_SAVED_HITS));

  ESL_ALLOC(saved_hits->unsrt, sizeof(SPLICE_HIT_INFO)  * default_nalloc);
  ESL_ALLOC(saved_hits->srt,   sizeof(SPLICE_HIT_INFO*) * default_nalloc);

  saved_hits->Nalloc = default_nalloc;
  saved_hits->N      = 0;

  saved_hits->is_sorted = FALSE;

  return saved_hits;

ERROR:
    p7_splicehits_DestroySavedHits(saved_hits);
    return NULL;

}


/* Function:  p7_splicehits_GrowSavedList()
 * Synopsis:  Reallocates a larger saved hits list, if needed.
 *
 * Purpose:   If <SPLICE_SAVED_HITS> cannot hold another
 *            <SPLICE_HIT_INFO>, doubles the internal 
 *            allocation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *   
 */
int
p7_splicehits_GrowSavedHits(SPLICE_SAVED_HITS *saved_hits)
{
  int status;

  if(saved_hits->N == saved_hits->Nalloc) {

    saved_hits->Nalloc *= 2;

    ESL_REALLOC(saved_hits->unsrt, sizeof(SPLICE_HIT_INFO)  * saved_hits->Nalloc);
    ESL_REALLOC(saved_hits->srt,   sizeof(SPLICE_HIT_INFO*) * saved_hits->Nalloc);
  }


  return eslOK;

  ERROR:
    p7_splicehits_DestroySavedHits(saved_hits); 
    return status;

}


/* Function: p7_splicehits_DestroySavedHits(); 
 *
 * Purpose:  Frees a <SPLICE_SAVED_HITS>
 */
void 
p7_splicehits_DestroySavedHits(SPLICE_SAVED_HITS *saved_hits)
{

  if(saved_hits == NULL) return;
  if(saved_hits->srt != NULL) free(saved_hits->srt);
  if(saved_hits->unsrt != NULL) free(saved_hits->unsrt);
  free(saved_hits);

  return;
}

/* Function:  p7_splicehits_CreateNext()
 * Synopsis:  Get pointer to new structure for recording hit info.
 *
 * Purpose:   Ask the saved hits object <saved_hits> to do any 
 *            necessary internal allocation to add a new,
 *            empty hit info to its list; return a pointer to
 *            this new <SPLICE_HIT_INFO> structure for data to be 
 *            filled in by the caller.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_splicehits_CreateNext(SPLICE_SAVED_HITS *saved_hits, SPLICE_HIT_INFO **ret_info)
{
  int status;
  
  if ((status = p7_splicehits_GrowSavedHits(saved_hits)) != eslOK) goto ERROR;
  *ret_info = &(saved_hits->unsrt[saved_hits->N]);
  saved_hits->N++;

  saved_hits->is_sorted = FALSE;

  return eslOK;

  ERROR:
  *ret_info = NULL;
   return status;

}

/* Function:  p7_splicehits_SortSavedHits()
 * Synopsis:  Sorts a saved hits list 
 *
 * Purpose:   Sorts a saved hits list by sequence index, strand,  
 *            (forward strand first), and sequence start postion
 *
 * Returns:   <eslOK> on success.
 */
int
p7_splicehits_SortSavedHits(SPLICE_SAVED_HITS *sh) 
{
  int i;
  if (sh->is_sorted) return eslOK;
  for (i = 0; i < sh->N; i++) sh->srt[i] = sh->unsrt + i;
  if (sh->N > 1)  qsort(sh->srt, sh->N, sizeof(SPLICE_HIT_INFO*), hit_sorter);
  sh->is_sorted = TRUE;
  return eslOK;
}

int
hit_sorter(const void *vh1, const void *vh2) {

 SPLICE_HIT_INFO *h1 = *((SPLICE_HIT_INFO **) vh1);
 SPLICE_HIT_INFO *h2 = *((SPLICE_HIT_INFO **) vh2);

 /* first key, seqidx (unique id for sequences), low to high */
 if      (h1->seqidx > h2->seqidx) return  1;
 else if (h1->seqidx < h2->seqidx) return -1;

 /* strand key, strand (0 for forward, 1 for reverse), low to high */
 if      (h1->strand > h2->strand) return  1;
 else if (h1->strand < h2->strand) return -1;

 if      (h1->seq_start > h2->seq_start) return  1;   // sort primarily from smallest to largest start pos
 else if (h1->seq_start < h2->seq_start) return -1;
 else if (h1->seq_end   < h2->seq_end)   return  1;   // secondarily, larger to smallest end position (i.e. longer hit first)
 else if (h1->seq_end   > h2->seq_end)   return -1;
 else                                    return  0;
 
}

/* Function:  p7_splicehits_MergeSavedHits()
 * Synopsis:  Merge two saved hits lists.
 *
 * Purpose:   Merge <sh2> into <sh1>. Upon return, <sh1>
 *            contains the sorted, merged list. <sh2>
 *            is effectively destroyed; caller should
 *            not access it further, and may as well free
 *            it immediately.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure, and
 *            both <sh1> and <sh2> remain valid.
 */
int
p7_splicehits_MergeSavedHits(SPLICE_SAVED_HITS *sh1, SPLICE_SAVED_HITS *sh2)
{

  int     i, j, k;
  SPLICE_HIT_INFO **new_srt = NULL;
  SPLICE_HIT_INFO *ori1 = sh1->unsrt;
  SPLICE_HIT_INFO *new2;
  uint64_t Nalloc;  
  int      status;

  if(sh2->N <= 0) return eslOK;

  p7_splicehits_SortSavedHits(sh1);
  p7_splicehits_SortSavedHits(sh2);

  Nalloc = sh1->N + sh2->N;
  ESL_REALLOC(sh1->unsrt, sizeof(SPLICE_HIT_INFO)  * Nalloc);
  ESL_ALLOC  (new_srt,    sizeof(SPLICE_HIT_INFO*) * Nalloc); 

 for (i = 0; i < sh1->N; i++)
   sh1->srt[i] = sh1->unsrt + (sh1->srt[i] - ori1);

  new2 = sh1->unsrt + sh1->N;
  memcpy(new2, sh2->unsrt, sizeof(SPLICE_HIT_INFO) * sh2->N);
  
  for (i=0,j=0,k=0; i < sh1->N && j < sh2->N ; k++)
    new_srt[k] = (hit_sorter(&sh1->srt[i], &sh2->srt[j]) > 0) ? new2 + (sh2->srt[j++] - sh2->unsrt) : sh1->srt[i++];  

  while (i < sh1->N) new_srt[k++] = sh1->srt[i++];
  while (j < sh2->N) new_srt[k++] = new2 + (sh2->srt[j++] - sh2->unsrt);

  free(sh1->srt);
  sh1->srt    = new_srt;
  sh1->Nalloc = Nalloc;
  sh1->N     += sh2->N;

  
  return eslOK;

  ERROR:
    if (new_srt != NULL) free(new_srt);
    return status;
}



int
p7_splicehits_AssignNodes(SPLICE_GRAPH *graph, SPLICE_SAVED_HITS *sh, int first, int last) 
{

  int     i,h;
  int     hit_hmm_min, hit_hmm_max;
  int     info_hmm_min, info_hmm_max;
  int     hmm_overlap_min, hmm_overlap_max;
  int     hmm_overlap_len;
  int     max_hmm_overlap;
  int     max_seq_overlap;
  int64_t hit_seq_min, hit_seq_max;
  int64_t info_seq_min, info_seq_max;
  int64_t seq_overlap_min, seq_overlap_max;
  int64_t seq_overlap_len;
  P7_TOPHITS *th;

  th = graph->th;
  
  for(i = first; i <= last; i++) {
    info_hmm_min = sh->srt[i]->hmm_start;
    info_hmm_max = sh->srt[i]->hmm_end;
    info_seq_min = ESL_MIN(sh->srt[i]->seq_start, sh->srt[i]->seq_end);
    info_seq_max = ESL_MAX(sh->srt[i]->seq_start, sh->srt[i]->seq_end);

    max_hmm_overlap = 0;
    max_seq_overlap = 0;
  
    for(h = 0; h < th->N; h++) {
      hit_hmm_min = th->hit[h]->dcl->ihmm;
      hit_hmm_max = th->hit[h]->dcl->jhmm;
      hit_seq_min = ESL_MIN(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
      hit_seq_max = ESL_MAX(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
 
      hmm_overlap_min = ESL_MAX(hit_hmm_min, info_hmm_min);
      hmm_overlap_max = ESL_MIN(hit_hmm_max, info_hmm_max);
      hmm_overlap_len = hmm_overlap_max - hmm_overlap_min + 1;
       
      seq_overlap_min = ESL_MAX(hit_seq_min,info_seq_min);
      seq_overlap_max = ESL_MIN(hit_seq_max, info_seq_max);
      seq_overlap_len = seq_overlap_max - seq_overlap_min + 1;

      if(hmm_overlap_len > max_hmm_overlap && seq_overlap_len > max_seq_overlap) {
        max_hmm_overlap = hmm_overlap_len;
        max_seq_overlap = seq_overlap_len;
             
        sh->srt[i]->node_id = h;
      }
    }
  }
 
  return eslOK;
}

int 
p7_splicehits_RemoveDuplicates(SPLICE_SAVED_HITS *sh) 
{

  int     i, j;
  int     s_i, s_j, e_i, e_j, len_i, len_j;
  int     intersect_alistart, intersect_aliend, intersect_alilen;
  int     intersect_hmmstart, intersect_hmmend, intersect_hmmlen;
  int     remove;

  if(sh->N < 2) return eslOK;

  if(!sh->is_sorted) p7_splicehits_SortSavedHits(sh);
  
  j = 0;
  for (i = 1; i < sh->N; i++)  
  {

    if(sh->srt[j]->seqidx != sh->srt[i]->seqidx || 
       sh->srt[j]->strand != sh->srt[i]->strand) {
       j = i; 
       continue;
    }

    s_j   = sh->srt[j]->seq_start;
    e_j   = sh->srt[j]->seq_end;      
   
    if(sh->srt[j]->strand) ESL_SWAP(s_j, e_j, int);
    len_j = e_j - s_j + 1;

    s_i   = sh->srt[i]->seq_start;
    e_i   = sh->srt[i]->seq_end;

    if (sh->srt[j]->strand) ESL_SWAP(s_i, e_i, int);
    len_i = e_i - s_i + 1 ; 

    intersect_alistart  = ESL_MAX(s_i, s_j);
    intersect_aliend    = ESL_MIN(e_i, e_j);
    intersect_alilen    = intersect_aliend - intersect_alistart + 1;

    intersect_hmmstart = ESL_MAX(sh->srt[i]->hmm_start, sh->srt[j]->hmm_start);   
    intersect_hmmend   = ESL_MIN(sh->srt[i]->hmm_end,   sh->srt[j]->hmm_end);
    intersect_hmmlen   = intersect_hmmend - intersect_hmmstart + 1;
 
    if(  intersect_hmmlen > 0              && // hmm corrds overlap and  
      (( s_i >= s_j-3 && s_i <= s_j+3)     || // at least one side is essentially flush
       ( e_i >= e_j-3 && e_i <= e_j+3)     ||
       ( intersect_alilen >= len_i * 0.95) || // or one of the hits covers >90% of the other
       ( intersect_alilen >= len_j * 0.95))) {

      remove = len_i > len_j ? j : i;

      sh->srt[remove]->duplicate = TRUE; 
    }
    else j = i;
  }


  return eslOK;
}

/*****************************************************************
 * 2. Debugging tools.
 *****************************************************************/


/* Function:  p7_splicephits_Dump()
 *
 * Purpose: Dumps hit info for each hit in <sh>
 *
 */
void
p7_splicehits_Dump(FILE *fp, SPLICE_SAVED_HITS *sh)
{
  int i;
  SPLICE_HIT_INFO *hi;

  fprintf(fp, "  Saved Hit Info\n");
  fprintf(fp, "  %9s %5s %9s %9s %10s %10s %4s %5s\n", "SeqIdx", "Strand", "hmm_start", "hmm_end", "seq_start", "seq_end", "node", "dup");
  for(i = 0; i < sh->N ; i++) {
    if(sh->is_sorted) hi = sh->srt[i];
    else              hi = &sh->unsrt[i];
 
    fprintf(fp, "  %9" PRId64 " %5s %9d %9d %10" PRId64 " %10" PRId64 " %4d %5s\n", 
    hi->seqidx, (hi->strand ? "-" : "+"), hi->hmm_start, hi->hmm_end, hi->seq_start, hi->seq_end, hi->node_id+1, (hi->duplicate? "TRUE" : "FALSE"));
  
  }

  return;
}

