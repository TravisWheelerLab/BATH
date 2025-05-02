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
  void    *p;
  SPLICE_HIT_INFO **new_srt = NULL;
  SPLICE_HIT_INFO *ori1;
  SPLICE_HIT_INFO *new2;
  uint64_t Nalloc;  
  int      status;

  if(sh2->N <= 0) return eslOK;

  p7_splicehits_SortSavedHits(sh1);
  p7_splicehits_SortSavedHits(sh2);

  Nalloc = sh1->N + sh2->N;
  ESL_RALLOC(sh1->unsrt, p, sizeof(SPLICE_HIT_INFO)  * Nalloc);
  ESL_ALLOC (new_srt,       sizeof(SPLICE_HIT_INFO*) * Nalloc); 

  ori1 = sh1->unsrt; 
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
  fprintf(fp, "  %9s %5s %9s %9s %10s %10s %9s \n", "SeqIdx", "Strand", "hmm_start", "hmm_end", "seq_start", "seq_end", "ali_score");
  for(i = 0; i < sh->N ; i++) {
    if(sh->is_sorted) hi = sh->srt[i];
    else              hi = &sh->unsrt[i];
 
    fprintf(fp, "  %9d %5s %9d %9d %10" PRId64 " %10" PRId64 " %9.2f %10.2f\n", 
    hi->seqidx, (hi->strand ? "-" : "+"), hi->hmm_start, hi->hmm_end, hi->seq_start, hi->seq_end, hi->aliscore);
  
  }

  return;
}

