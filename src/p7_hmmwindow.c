/* The Plan7 HMMWINDOW data structure, which holds a compact representation
 * of substitution scores and maximal extensions, used by nhmmer.
 *
 * Repurposed to store and manage ungapped exon seeds for spliced alignments
 *
 * Contents:
 *   1. The P7_HMMWINDOW object: allocation, initialization, destruction.
 *   2. Splcing Specific Functions
 *   3. Debugging tools 
 *
 */
#include "p7_config.h"
#include <string.h>
#include "hmmer.h"
#include "p7_splice.h"

/*********************************************************************
 *# 1. The P7_MSVDATA object: allocation, initialization, destruction.
 *********************************************************************/

/* Function:  p7_splicehits_CreateLists()
 * Synopsis:  Allocates an P7_HMM_WINDOWLIST with windows.
 *
 * Returns:   a pointer to the new <P7_HMM_WINDOWLIST> 
 *            structure on success.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_HMM_WINDOWLIST*
p7_hmmwindow_CreateList(void)
{

  P7_HMM_WINDOWLIST *hw = NULL;
  int status;
  
  ESL_ALLOC(hw, sizeof(P7_HMM_WINDOWLIST));

  hw->windows = NULL;
  p7_hmmwindow_init(hw);

  return hw;

ERROR:
  p7_hmmwindow_DestroyList(hw);
  return NULL;

}

/* Function:  p7_hmmwindow_init()
 *
 * Synopsis:  initialize the object used to store a list of sequence windows
 *
 * Returns:   eslEMEM in event of allocation failure, otherwise eslOK
 */
int
p7_hmmwindow_init (P7_HMM_WINDOWLIST *list) {
  int status;
  list->size   = 10000;
  list->count  = 0;
  list->sorted = 0;
  ESL_ALLOC(list->windows, list->size * sizeof(P7_HMM_WINDOW));

  return eslOK;

ERROR:
  return eslEMEM;

}

/* Function:  p7_hmmwindow_new()
 *
 * Synopsis:  Return a pointer to the next window element on the list
 *
 * Purpose:   Accepts <id>, <pos>, <k>, <length>, <score>,
 *            and <complementarity>, assigns those to the next window
 *            element, then returns it, increasing the size of the
 *            list, if necessary.
 *
 * Returns:   NULL in event of allocation failure, otherwise pointer to
 *            the next seed diagonal
 */
P7_HMM_WINDOW *
p7_hmmwindow_new (P7_HMM_WINDOWLIST *list, uint32_t id, uint32_t pos, uint16_t k, uint32_t length, float score, uint8_t complementarity, uint32_t target_len) {
  int status;
  P7_HMM_WINDOW *window;

  if (list->count == list->size) {
    list->size *= 4;
    ESL_REALLOC(list->windows, list->size * sizeof(P7_HMM_WINDOW));
  }
  window = list->windows + list->count;

  window->id               = id;
  window->n                = pos;
  window->k                = k;
  window->length           = length;
  window->score            = score;
  window->complementarity  = complementarity;
  window->target_len       = target_len;
  window->duplicate        = FALSE;
  window->pass_forward     = FALSE;
  window->is_seed          = FALSE;
  list->count++;
 
  list->sorted = FALSE;
  
  return window;

ERROR:
  return NULL;
}


/* Function: p7_hmmwindow_DestroyList();
 *
 * Purpose:  Frees a <P7_HMM_WINDOWLIST>
 */
void
p7_hmmwindow_DestroyList(P7_HMM_WINDOWLIST *hw)
{

  if(hw == NULL) return;

  if(hw->windows != NULL) free(hw->windows);
  free(hw);

  return;
}


/* window_sorter(): qsort's pawn, below */
static int
window_pos_sorter(const void *vw1, const void *vw2)
{
  P7_HMM_WINDOW w1 = *((P7_HMM_WINDOW *) vw1);  /* don't ask. don't change. Don't Panic. */
  P7_HMM_WINDOW w2 = *((P7_HMM_WINDOW *) vw2);

  if      (w1.n > w2.n) return  1;
  else if (w1.n < w2.n) return -1;
  else                  return  0;
}



/* Function:  p7_hmmwindow_SortByStart()
 * Synopsis:  Sorts an hmm window list by sequence position.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_hmmwindow_SortByStart(P7_HMM_WINDOWLIST *w)
{
  if (w->count > 1)  qsort(w->windows, w->count, sizeof(P7_HMM_WINDOW), window_pos_sorter);
  return eslOK;
}


int
window_seq_sorter(const void *vw1, const void *vw2) {

 P7_HMM_WINDOW w1 = *((P7_HMM_WINDOW *) vw1);  /* don't ask. don't change. Don't Panic. */
 P7_HMM_WINDOW w2 = *((P7_HMM_WINDOW *) vw2);

 /* first key, seqidx (unique id for sequences), low to high */
 if      (w1.id > w2.id) return  1;
 else if (w1.id < w2.id) return -1;

 /* strand key, strand (0 for forward, 1 for reverse), low to high */
 if      (w1.complementarity > w2.complementarity) return  1;
 else if (w1.complementarity < w2.complementarity) return -1;

 if      (w1.n      > w2.n)      return  1;   // sort primarily from smallest to largest start pos
 else if (w1.n      < w2.n)      return -1;
 else if (w1.length < w2.length) return  1;   // secondarily, larger to smallest end position (i.e. longer hit first)
 else if (w1.length > w2.length) return -1;
 else                            return  0;

}



/* Function:  p7_hmmwindow_SortBySeq()
 * Synopsis:  Sorts an  hmm window list by Sequence id, strand and sequence position.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_hmmwindow_SortBySeq(P7_HMM_WINDOWLIST *w)
{
  if(w->sorted) return eslOK;

  if (w->count > 1) qsort(w->windows, w->count, sizeof(P7_HMM_WINDOW), window_seq_sorter); 

  w->sorted = TRUE;

  return eslOK;
}

/*********************************************************************
 *# 1. Splicing Functions
 *********************************************************************/


/* Function:  p7_hmmwindow_Merge()
 * Synopsis:  Merge two hmm window lists.
 *
 * Purpose:   Merge <hw2> into <hw1>. Upon return, <hw1>
 *            contains the sorted, merged list. 
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure, and
 *            both <hw1> and <hw2> remain valid.
 */
int
p7_hmmwindow_Merge(P7_HMM_WINDOWLIST *hw1, P7_HMM_WINDOWLIST *hw2)
{

  P7_HMM_WINDOW *end;
  int status;

  if(hw2->count == 0) return eslOK;
 
   if(hw1->count + hw2->count > hw1->size) {
     hw1->size += hw2->size;
     ESL_REALLOC(hw1->windows, sizeof(P7_HMM_WINDOW) * hw1->size);
   }

   end = hw1->windows + hw1->count;
   memcpy(end, hw2->windows, sizeof(P7_HMM_WINDOW) * hw2->count);
   hw1->count += hw2->count;

   hw1->sorted = FALSE;
   p7_hmmwindow_SortBySeq(hw1);

   return eslOK;

   ERROR:
     return status;
}



/* Function:  p7_hmmwindow_RemoveDuplicates()
 * Synopsis:  Remove duplicate entires in P7_HMM_WINDOWLIST
 *
 * Purpose:   Find any duplicate (ie overlapping) hmm windwow
 *            and "remove" the shorter one by setting duplicate
 *            to TRUE. Also check for overlaps with hits in the 
 *            P7_TOPHITS, and set those as windwow as duplicate.
 *
 * Returns:   <eslOK> on success.
 *
 */
int
p7_hmmwindow_RemoveDuplicates(P7_HMM_WINDOWLIST *hw, P7_TOPHITS *th, double F3) 
{

  int     i, j;
  int     j_start;
  int64_t s_i, s_j, e_i, e_j;  
  int     len_i, len_j;
  int64_t intersect_alistart, intersect_aliend;
  int     intersect_alilen;
  int     intersect_hmmstart, intersect_hmmend, intersect_hmmlen;
  int     strand;
  int     remove;

  j = 0;
  for (i = 1; i < hw->count; i++) {
    
    if(hw->windows[j].id != hw->windows[i].id ||
       hw->windows[j].complementarity != hw->windows[i].complementarity) {
      j = i;
      continue;
    }
 
    s_j = hw->windows[j].n;
    e_j = hw->windows[j].n + hw->windows[j].length - 1;

    if(hw->windows[j].complementarity) ESL_SWAP(s_j, e_j, int);
    len_j = e_j - s_j + 1;

    s_i = hw->windows[i].n;
    e_i = hw->windows[i].n + hw->windows[i].length - 1;   
    len_i = e_i - s_i + 1;

    intersect_alistart  = ESL_MAX(s_i, s_j);
    intersect_aliend    = ESL_MIN(e_i, e_j);
    intersect_alilen    = intersect_aliend - intersect_alistart + 1;

    intersect_hmmstart = ESL_MIN(hw->windows[j].k - (hw->windows[j].length/3) + 1,  hw->windows[i].k - (hw->windows[i].length/3) + 1);
    intersect_hmmend   = ESL_MAX(hw->windows[j].k, hw->windows[i].k);
    intersect_hmmlen   = intersect_hmmend - intersect_hmmstart + 1;

    if(  intersect_hmmlen > 0              && // hmm corrds overlap and
      (( s_i >= s_j-3 && s_i <= s_j+3)     || // at least one side is essentially flush
       ( e_i >= e_j-3 && e_i <= e_j+3)     ||
       ( intersect_alilen >= len_i * 0.95) || // or one of the hits covers >90% of the other
       ( intersect_alilen >= len_j * 0.95))) {

      remove     = len_i > len_j ? j : i;
      hw->windows[remove].duplicate = TRUE;
    }
    else j = i;
  }

  j_start = 0; 
  for (i = 0; i < th->N; i++)
  {
    if ((th->hit[i]->flags & p7_IS_DUPLICATE)) continue;
    if(!(th->hit[i]->flags & p7_IS_REPORTED) && exp(th->hit[i]->sum_lnP) >= F3) continue;

    s_i   = th->hit[i]->dcl[0].iali;
    e_i   = th->hit[i]->dcl[0].jali;

    strand = (s_i < e_i ? p7_NOCOMPLEMENT : p7_COMPLEMENT);
    if (strand) ESL_SWAP(s_i, e_i, int);
    len_i = e_i - s_i + 1 ;
    
     while(j_start < hw->count &&
          (hw->windows[j_start].id != th->hit[i]->seqidx ||
           hw->windows[j_start].complementarity != strand)) j_start++;

    j = j_start;
    while(j < hw->count &&
          hw->windows[j].id == th->hit[i]->seqidx &&
          hw->windows[j].complementarity == strand) {

      if( hw->windows[j].duplicate)    { j++; continue; }
      if(!hw->windows[j].pass_forward) { j++; continue; }

      s_j = hw->windows[j].n;
      e_j = hw->windows[j].n + hw->windows[j].length - 1;
      len_j = e_j - s_j + 1; 

      intersect_alistart = ESL_MAX(s_i, s_j);
      intersect_aliend   = ESL_MIN(e_i, e_j);
      intersect_alilen   = intersect_aliend - intersect_alistart + 1;

      intersect_hmmstart = ESL_MAX(th->hit[i]->dcl[0].ihmm, hw->windows[i].k - (hw->windows[i].length/3) + 1);
      intersect_hmmend   = ESL_MIN(th->hit[i]->dcl[0].jhmm, hw->windows[i].k);
      intersect_hmmlen   = intersect_hmmend - intersect_hmmstart + 1;

      if(  intersect_hmmlen > 0              && // hmm corrds overlap and
        (( s_i >= s_j-3 && s_i <= s_j+3)     || // at least one side is essentially flush
         ( e_i >= e_j-3 && e_i <= e_j+3)     ||
         ( intersect_alilen >= len_i * 0.9) || // or one of the hits covers >90% of the other
         ( intersect_alilen >= len_j * 0.9))) {
         hw->windows[j].duplicate = TRUE;
      }
      j++;
    }
  }

  return eslOK;

}


/* Function:  p7_hmmwindow_GetSeedHits()
 * Synopsis:  Find seed hits in P7_HMM_WINDOWLIST and transfer to P7_TOPHITS
 *
 * Purpose:   Find any hmm winodws in <hw> that are upstream or 
 *            downstream of a hit in the <th> and sdd them to 
 *            the <seed_hits>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    p7_Fail on sequence file read error
 */
P7_TOPHITS*
p7_hmmwindow_GetSeedHits(P7_HMM_WINDOWLIST *hw, const P7_TOPHITS *th, P7_HMM *hmm, P7_FS_PROFILE *gm_fs, ESL_SQFILE *seq_file, ESL_GENCODE *gcode, double F3)
{

  int i, h, y, z;
  int i_start;
  int strand;
  int hit_min, hit_max;
  int window_min, window_max;
  int hmm_start, hmm_end;
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
 
   /* Find all windows that within MAX_INTRON_LENG upstream or downstram of a top hit */
   for(h = 0; h < th->N; h++) {

    if ((th->hit[h]->flags & p7_IS_DUPLICATE)) continue;
    if(!(th->hit[h]->flags & p7_IS_REPORTED) && exp(th->hit[h]->sum_lnP) >= F3) continue;

    strand = (th->hit[h]->dcl->iali < th->hit[h]->dcl->jali ? p7_NOCOMPLEMENT : p7_COMPLEMENT);
    hit_min  = ESL_MIN(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
    hit_max  = ESL_MAX(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);

    i = i_start;
    while(i < hw->count) {

      if(hw->windows[i].id < th->hit[h]->seqidx)  { i++; continue; }
      if(hw->windows[i].complementarity < strand) { i++; continue; }

      /*If this is the first hit on a new sequence or strand, save
       * the first postion in the saved hits that coresponds to
       * that sequence and strand for future use */

      if(last_seqidx != th->hit[h]->seqidx || last_strand != strand) i_start = i;

      if(hw->windows[i].id > th->hit[h]->seqidx) break;
      if(hw->windows[i].complementarity > strand)       break;

      if(hw->windows[i].duplicate) { i++; continue; }
      if(hw->windows[i].is_seed)   { i++; continue; }

      window_min = hw->windows[i].n;
      window_max = hw->windows[i].n + hw->windows[i].length - 1;

      if(hit_min - window_max > MAX_INTRON_LENG) { i++; continue; }
      if(window_min - hit_max > MAX_INTRON_LENG) break;

      hmm_start = hw->windows[i].k - (hw->windows[i].length/3) + 1;
      hmm_end   = hw->windows[i].k;
      /* Is saved hit upstream of top hit  */
      if(hmm_start <= th->hit[h]->dcl->ihmm ||
         hmm_end   <= th->hit[h]->dcl->jhmm) {

        if (( strand  && window_min > th->hit[h]->dcl->iali) ||
           ((!strand) && window_max < th->hit[h]->dcl->iali))  {
          hw->windows[i].is_seed = TRUE;
          i++;
          continue;
        }
      }

      // Is saved hit downstearm of top hit
      if(th->hit[h]->dcl->ihmm <= hmm_start ||
         th->hit[h]->dcl->jhmm <= hmm_end) {

        if (( strand  && th->hit[h]->dcl->iali > window_min ) ||
           ((!strand) && th->hit[h]->dcl->iali < window_max)) {
          hw->windows[i].is_seed = TRUE;
          i++;
          continue;
        }
      }
      i++;
    }
    last_seqidx = th->hit[h]->seqidx;
    last_strand = strand;
  }

  /* Add all windows where is_seed = TRUE to seed_hits */
  last_seqidx = -1;
  last_strand = -1;
  i_start = 0;
  for(i = 0 ; i < hw->count; i++) {

    if(!hw->windows[i].is_seed) continue;

    /* If we have a new sequence or strand move to the start of the sequence */
    if(hw->windows[i].id              != last_seqidx ||
       hw->windows[i].complementarity != last_strand) {

      if(last_strand != -1 ) esl_sq_Reuse(dbsq_dna);

      for(h = 0; h < th->N; h++) {
        if(th->hit[h]->seqidx == hw->windows[i].id) {
          seqname = th->hit[h]->name;
          break;
        }
      }

      esl_sqfile_PositionByKey(dbfp, seqname);

      status = esl_sqio_ReadWindow(dbfp, 0, window_len, dbsq_dna);
      if(hw->windows[i].complementarity == p7_COMPLEMENT)
        esl_sq_ReverseComplement(dbsq_dna);
      seq_max = ESL_MAX(dbsq_dna->start, dbsq_dna->end);
    }

    window_min = hw->windows[i].n;
    window_max = hw->windows[i].n + hw->windows[i].length - 1;    
    /* Fetch the next winow until it contains the current window */
    while(status == eslOK && window_max > seq_max) {
      status = esl_sqio_ReadWindow(dbfp, gm_fs->max_length*3, window_len, dbsq_dna);
      seq_max = ESL_MAX(dbsq_dna->start, dbsq_dna->end);

      if(window_max <= seq_max && hw->windows[i].complementarity)
        esl_sq_ReverseComplement(dbsq_dna);
    }
   
    /* handle errors */
    switch(status) {
      case eslEFORMAT:
        esl_fatal("Parse failed (sequence file %s):\n%s\n",
                   dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
        break;
      case eslEOD:
        if(window_max > seq_max) esl_fatal("Failed to find subsequence in sequence file %s", dbfp->filename);
        else esl_sq_Reuse(dbsq_dna);
        break;
      case eslEOF:
      case eslOK:
        /* do nothing */
         break;
      default:
        esl_fatal("Unexpected error %d reading sequence file %s", status, dbfp->filename);
    }

    p7_tophits_CreateNextHit(seed_hits, &hit);
    hit->seqidx  = hw->windows[i].id;
    hit->dcl     = p7_domain_Create_empty();
    hit->dcl->tr = p7_trace_fs_Create();

    /* repurpose is_reported for hits that passed the forward filter */
    if(hw->windows[i].pass_forward) hit->dcl->is_reported = TRUE;

    hit->dcl->ihmm = hw->windows[i].k - (hw->windows[i].length/3) + 1;
    hit->dcl->jhmm = hw->windows[i].k;

    if(hw->windows[i].complementarity) {
      hit->dcl->iali = hw->windows[i].n + hw->windows[i].length - 1;
      hit->dcl->jali = hw->windows[i].n;
    }
    else {
      hit->dcl->iali = hw->windows[i].n;
      hit->dcl->iali = hw->windows[i].n + hw->windows[i].length - 1;
    }
 
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
 
    last_seqidx = hw->windows[i].id;
    last_strand = hw->windows[i].complementarity;    
  }

  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abcDNA);
  esl_sq_Destroy(dbsq_dna);
  esl_sqfile_Close(dbfp);

  return seed_hits;
}


/*****************************************************************
 * 3. Debugging tools.
 *****************************************************************/

/* Function:  p7_hmmwindow_Dump()
 *
 * Purpose: Dumps info for each window in <hw>
 *
 */
void
p7_hmwwindow_Dump(FILE *fp, P7_HMM_WINDOWLIST *hw)
{
  int i;
  P7_HMM_WINDOW *w;

  fprintf(fp, "  Window Info\n");
  fprintf(fp, "  %10s %9s %5s %9s %9s %10s %10s %5s %5s %5s\n", "hit", "SeqIdx", "Strand", "hmm_start", "hmm_end", "seq_start", "seq_end", "fwd", "dup", "seed");
  for(i = 0; i < hw->count ; i++) {
    w = &(hw->windows[i]);
    if(w->complementarity)
      fprintf(fp, "  %10d %9" PRId64 " %5s %9d %9d %10" PRId64 " %10" PRId64 " %5s %5s %5s\n",
      i+1, w->id, "-", w->k-(w->length/3)+1, w->k, w->n+w->length-1, w->n, (w->pass_forward? "TRUE" : "FALSE"), (w->duplicate? "TRUE" : "FALSE"), (w->is_seed? "TRUE" : "FALSE"));
    else
      fprintf(fp, "  %10d %9" PRId64 " %5s %9d %9d %10" PRId64 " %10" PRId64 " %5s %5s %5s\n",
      i+1, w->id, "+", w->k-(w->length/3)+1, w->k, w->n, w->n+w->length-1, (w->pass_forward? "TRUE" : "FALSE"), (w->duplicate? "TRUE" : "FALSE"), (w->is_seed? "TRUE" : "FALSE"));

  }

  return;
}
