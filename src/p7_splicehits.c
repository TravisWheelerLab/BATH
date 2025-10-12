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
p7_splicehits_RemoveDuplicates(SPLICE_SAVED_HITS *sh, P7_TOPHITS *th, double F3) 
{

  int     i, j;
  int     j_start;
  int     strand;
  int     s_i, s_j, e_i, e_j, len_i, len_j;
  int     intersect_alistart, intersect_aliend, intersect_alilen;
  int     intersect_hmmstart, intersect_hmmend, intersect_hmmlen;
  int     remove; 

  if(sh->N < 2) return eslOK;

  if(!sh->is_sorted) p7_splicehits_SortSavedHits(sh);
  
  /* If two saved hits overlap set the shorter lone as duplicate */ 
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

      remove     = len_i > len_j ? j : i;
 
      sh->srt[remove]->duplicate = TRUE; 
    }
    else j = i;
  }

  j_start = 0;
  /* If a saved hit overlaps with a top hit, set the saved hit as duplicate */
  for (i = 0; i < th->N; i++)  
  {

   if ((th->hit[i]->flags & p7_IS_DUPLICATE)) continue;
   if(!(th->hit[i]->flags & p7_IS_REPORTED) && exp(th->hit[i]->sum_lnP) >= F3) continue;
    
    s_i   = th->hit[i]->dcl[0].iali;
    e_i   = th->hit[i]->dcl[0].jali;
    strand = (s_i < e_i ? p7_NOCOMPLEMENT : p7_COMPLEMENT);
    if (strand) ESL_SWAP(s_i, e_i, int);
    len_i = e_i - s_i + 1 ; 

    while(sh->srt[j_start]->seqidx != th->hit[i]->seqidx || 
          sh->srt[j_start]->strand != strand) j_start++; 

    j = j_start;
    while(sh->srt[j]->seqidx == th->hit[i]->seqidx &&
          sh->srt[j]->strand == strand) {
      s_j   = sh->srt[j]->seq_start;
      e_j   = sh->srt[j]->seq_end;      

      if(sh->srt[j]->strand) ESL_SWAP(s_j, e_j, int);
      len_j = e_j - s_j + 1;

      intersect_alistart  = ESL_MAX(s_i, s_j);
      intersect_aliend    = ESL_MIN(e_i, e_j);
      intersect_alilen    = intersect_aliend - intersect_alistart + 1;

      intersect_hmmstart = ESL_MAX(th->hit[i]->dcl[0].ihmm, sh->srt[j]->hmm_start);   
      intersect_hmmend   = ESL_MIN(th->hit[i]->dcl[0].jhmm, sh->srt[j]->hmm_end);
      intersect_hmmlen   = intersect_hmmend - intersect_hmmstart + 1;
 
      if(  intersect_hmmlen > 0              && // hmm corrds overlap and  
        (( s_i >= s_j-3 && s_i <= s_j+3)     || // at least one side is essentially flush
         ( e_i >= e_j-3 && e_i <= e_j+3)     ||
         ( intersect_alilen >= len_i * 0.95) || // or one of the hits covers >90% of the other
         ( intersect_alilen >= len_j * 0.95))) {

        sh->srt[j]->duplicate = TRUE; 
      }
      j++;
      if(j == sh->N) break;
    }
  }
 

  return eslOK;
}

P7_TOPHITS*
p7_splicehits_GetSeedHits(SPLICE_SAVED_HITS *sh, const P7_TOPHITS *th, P7_HMM *hmm, P7_FS_PROFILE *gm_fs, ESL_SQFILE *seq_file, ESL_GENCODE *gcode, double F3) 
{

  int i, h, y, z;
  int i_start;
  int strand;
  int hit_min, hit_max;
  int seed_min, seed_max;
  int seq_max;
  int last_seqidx, last_strand;
  int window_len;
  char         *seqname;
  P7_HIT       *hit;
  P7_TOPHITS   *seed_hits;
  P7_BG        *bg;
  ESL_ALPHABET *abcDNA;
  ESL_SQFILE   *dbfp;
  ESL_SQ       *dbsq_dna;
  int status;

  window_len = (1024 * 256);

  
  bg = p7_bg_Create(gm_fs->abc);
  p7_bg_SetFilter(bg, gm_fs->M, gm_fs->compo);

  seed_hits = p7_tophits_Create();
  hit = NULL;

  dbfp  = NULL;
  
  /* Open sequence file */
  status = esl_sqfile_Open(seq_file->filename,seq_file->format,NULL,&dbfp); 
  if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",          seq_file);
  else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",            seq_file);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, seq_file);

  abcDNA = esl_alphabet_Create(eslDNA);
  dbsq_dna    = esl_sq_CreateDigital(abcDNA);

  esl_sqfile_SetDigital(dbfp, abcDNA);
  esl_sqfile_OpenSSI(dbfp,NULL);
 
  last_seqidx = -1;
  last_strand = -1; 
  i_start = 0;
  for(h = 0; h < th->N; h++) {

    if ((th->hit[h]->flags & p7_IS_DUPLICATE)) continue;
    if(!(th->hit[h]->flags & p7_IS_REPORTED) && exp(th->hit[h]->sum_lnP) >= F3) continue;

    strand = (th->hit[h]->dcl->iali < th->hit[h]->dcl->jali ? p7_NOCOMPLEMENT : p7_COMPLEMENT);   
    hit_min  = ESL_MIN(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
    hit_max  = ESL_MAX(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
   
    i = i_start;
    while(i < sh->N) {

      if(sh->srt[i]->seqidx < th->hit[h]->seqidx) { i++; continue; }
      if(sh->srt[i]->strand < strand)       { i++; continue; }
        
      /*If this is the first hit on a new sequence or strand, save 
       * the first postion in the saved hits that coresponds to 
       * that sequence and strand for future use */     

      if(last_seqidx != th->hit[h]->seqidx || last_strand != strand) i_start = i;

      if(sh->srt[i]->seqidx > th->hit[h]->seqidx) break;
      if(sh->srt[i]->strand > strand)       break;

      if(sh->srt[i]->duplicate) { i++; continue; }
      if(sh->srt[i]->is_seed)   { i++; continue; }

      seed_max = ESL_MAX(sh->srt[i]->seq_start, sh->srt[i]->seq_end);
      if(hit_min - seed_max > MAX_INTRON_LENG) { i++; continue; }

      seed_min = ESL_MIN(sh->srt[i]->seq_start, sh->srt[i]->seq_end);
      if(seed_min - hit_max > MAX_INTRON_LENG) break;

      // Is saved hit upstearm and within of top hit  
      if(sh->srt[i]->hmm_start <= th->hit[h]->dcl->ihmm ||
         sh->srt[i]->hmm_end   <= th->hit[h]->dcl->jhmm) {

        if (( strand  && sh->srt[i]->seq_end > th->hit[h]->dcl->iali) ||
           ((!strand) && sh->srt[i]->seq_end < th->hit[h]->dcl->iali))  {
          sh->srt[i]->is_seed = TRUE;
          sh->srt[i]->hit_id = h;
          i++; 
          continue;
        }
       
      }

      // Is saved hit downstearm and within of top hit 
	  if(th->hit[h]->dcl->ihmm <= sh->srt[i]->hmm_start ||
         th->hit[h]->dcl->jhmm <= sh->srt[i]->hmm_end) {

        if (( strand  && th->hit[h]->dcl->iali > sh->srt[i]->seq_end ) ||
           ((!strand) && th->hit[h]->dcl->iali < sh->srt[i]->seq_end)) { 
          sh->srt[i]->is_seed = TRUE;
          sh->srt[i]->hit_id = h;
          i++;
          continue;
        }
      }
      i++;  
    } 
    last_seqidx = th->hit[h]->seqidx;
    last_strand = strand;  
  }

  /* Add all saved hits where is_seed = TRUE to seed_hits */
  last_seqidx = -1;
  last_strand = -1;
  for(i = 0 ; i < sh->N; i++) {
   
    if(!sh->srt[i]->is_seed) continue;
    //printf("sh->srt[i]->seq_start %d sh->srt[i]->seq_end %d\n", sh->srt[i]->seq_start, sh->srt[i]->seq_end); 
    //fflush(stdout);
    /* If the saved hit is on a new sequence or strand sove to beginging of current sequence */
    if(sh->srt[i]->seqidx != last_seqidx ||
       sh->srt[i]->strand != last_strand) {

      if(last_strand != -1 ) esl_sq_Reuse(dbsq_dna);

      seqname = th->hit[sh->srt[i]->hit_id]->name;
      
      esl_sqfile_PositionByKey(dbfp, seqname);
      
      status = esl_sqio_ReadWindow(dbfp, 0, window_len, dbsq_dna);
      if(sh->srt[i]->strand == p7_COMPLEMENT) 
        esl_sq_ReverseComplement(dbsq_dna);
      seq_max = ESL_MAX(dbsq_dna->start, dbsq_dna->end);
    }
    
    seed_min = ESL_MIN(sh->srt[i]->seq_start, sh->srt[i]->seq_end);
    seed_max = ESL_MAX(sh->srt[i]->seq_start, sh->srt[i]->seq_end);

    while(status == eslOK && seed_max > seq_max) {
      status = esl_sqio_ReadWindow(dbfp, gm_fs->max_length*3, window_len, dbsq_dna);
      seq_max = ESL_MAX(dbsq_dna->start, dbsq_dna->end);
 
      if(seed_max <= seq_max && sh->srt[i]->strand == p7_COMPLEMENT)    
        esl_sq_ReverseComplement(dbsq_dna);
    }
    //printf("seqname %s seqidx %d dbsq_dna->start %d dbsq_dna->end %d\n", seqname, sh->srt[i]->seqidx, dbsq_dna->start, dbsq_dna->end);
    switch(status) {
      case eslEFORMAT:
        esl_fatal("Parse failed (sequence file %s):\n%s\n",
                   dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
        break;
      case eslEOD:
        if(seed_max > seq_max) esl_fatal("Failed to find subsequence in sequence file %s", dbfp->filename);
        else esl_sq_Reuse(dbsq_dna);
        break;
      case eslEOF:
      case eslOK:
        /* do nothing */
         break;
      default:
        esl_fatal("Unexpected error %d reading sequence file %s", status, dbfp->filename);
    }
    
    //printf("dbsq_dna->start %d dbsq_dna->end %d sh->srt[i]->seq_start %d sh->srt[i]->seq_end %d\n", dbsq_dna->start, dbsq_dna->end, sh->srt[i]->seq_start, sh->srt[i]->seq_end);
    p7_tophits_CreateNextHit(seed_hits, &hit);
    hit->seqidx  = sh->srt[i]->seqidx;
    hit->dcl     = p7_domain_Create_empty();
    hit->dcl->tr = p7_trace_fs_Create();
    
    hit->dcl->ihmm = sh->srt[i]->hmm_start;
    hit->dcl->jhmm = sh->srt[i]->hmm_end;
    hit->dcl->iali = sh->srt[i]->seq_start;
    hit->dcl->jali = sh->srt[i]->seq_end;
   // printf("ihmm %d jhmm %d iali %d jali %d\n", hit->dcl->ihmm, hit->dcl->jhmm, hit->dcl->iali, hit->dcl->jali); 
    /* Create trace for seed hit */
    p7_trace_fs_Append(hit->dcl->tr, p7T_S , 0, 0, 0);
    p7_trace_fs_Append(hit->dcl->tr, p7T_N , 0, 0, 0);
    p7_trace_fs_Append(hit->dcl->tr, p7T_B , 0, 0, 0);

    y = llabs(hit->dcl->iali - dbsq_dna->start) + 3;
    for(z = hit->dcl->ihmm; z <= hit->dcl->jhmm; z++) {
      p7_trace_fs_Append(hit->dcl->tr, p7T_M, z, y, 3);
      y+=3;
    }
    
    p7_trace_fs_Append(hit->dcl->tr, p7T_E, z-1, y-=3, 0);
    p7_trace_fs_Append(hit->dcl->tr, p7T_C, 0, y-=3, 0);
    p7_trace_fs_Append(hit->dcl->tr, p7T_T, 0, 0, 0);
   
    hit->dcl->scores_per_pos = NULL;
    p7_splice_ComputeAliScores_fs(hit->dcl, hit->dcl->tr, dbsq_dna, gm_fs, bg, TRUE);

    last_seqidx = sh->srt[i]->seqidx;
    last_strand  = sh->srt[i]->strand;     
  }
 
  p7_bg_Destroy(bg); 
  esl_alphabet_Destroy(abcDNA);
  esl_sq_Destroy(dbsq_dna); 
  esl_sqfile_Close(dbfp);
  
  //p7_splicehits_Dump(stdout, sh); 
  return seed_hits;
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
  fprintf(fp, "  %10s %9s %5s %9s %9s %10s %10s %5s %5s\n", "hit", "SeqIdx", "Strand", "hmm_start", "hmm_end", "seq_start", "seq_end", "dup", "seed");
  for(i = 0; i < sh->N ; i++) {
    if(sh->is_sorted) hi = sh->srt[i];
    else              hi = &sh->unsrt[i];
 
    fprintf(fp, "  %10d %9" PRId64 " %5s %9d %9d %10" PRId64 " %10" PRId64 " %5s %5s\n", 
    i+1, hi->seqidx, (hi->strand ? "-" : "+"), hi->hmm_start, hi->hmm_end, hi->seq_start, hi->seq_end, (hi->duplicate? "TRUE" : "FALSE"), (hi->is_seed? "TRUE" : "FALSE"));
  
  }

  return;
}

