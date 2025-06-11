/* SPLICE_EDGE: edge structure with all relavatnd coodiate and score data on the splicing of two hits
 *
 * Contents:
 *    1. The SPLICE_EDGE object.
 *    2. Access routines.
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

static int get_overlap_nuc_coords (SPLICE_EDGE *edge, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const ESL_SQ *splice_seq, int revcomp);
static int find_optimal_splice_site (SPLICE_EDGE *edge, SPLICE_PIPELINE *pli, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq);
static int select_splice_option (SPLICE_EDGE *edge, P7_PROFILE *sub_model, P7_FS_PROFILE *sub_fs_model, const ESL_GENCODE *gcode, const ESL_SQ *target_seq, float signal_score, int up_nuc_pos, int down_nuc_pos);

/*****************************************************************
 * 1. The SPLICE_EDGE structure.
 *****************************************************************/


/* Function:  p7_splicegraph_CreateEdge)
 * Synopsis:  Allocates a splice edge.
 *
 * Purpose:   Allocates a new <SPLICE_EDGE>.
 *
 * Returns:   a pointer to the new <SPLICE_EDGE> structure
 *            on success.
 *
 * Throws:    <NULL> on allocation error.
 */
SPLICE_EDGE*
p7_spliceedge_Create(void)
{

  int status;
  SPLICE_EDGE *edge;

  edge = NULL;
  ESL_ALLOC(edge, sizeof(SPLICE_EDGE));

  edge->frameshift = FALSE;
  edge->bypass_checked = FALSE;

  edge->splice_score = -eslINFINITY;
  edge->signal_score = -eslINFINITY;

  edge->prev = NULL;
  edge->next = NULL;

  return edge;

  ERROR:
    if (edge != NULL) free(edge);
    return NULL;
}

SPLICE_EDGE*
p7_spliceedge_ConnectHits(SPLICE_PIPELINE *pli, const P7_DOMAIN *upstream_domain, const P7_DOMAIN *downstream_domain, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq, int revcomp, int full_intron) 
{

  int      i, j, z;
  int      z1, z2;
  int      upstream_len;
  int      downstream_len;
  int      intron_len;
  int      intron_chars;
  int      amino_gap;
  int      nuc_extention;
  int      splice_len;
  int      intron_cnt;
  int      introns_in_gap;
  int      overlap_min, overlap_max;
  int      overlap_len;
  int      total_distance;
  int      distance_from_upstream;
  int      distance_from_downstream;
  int      intron_start, intron_end;
  P7_HMM        *sub_hmm;
  P7_FS_PROFILE *sub_fs_model;
  P7_GMX        *vit_mx;
  P7_TRACE      *tr;
  ESL_DSQ       *splice_dsq;
  SPLICE_EDGE   *edge;
  int           status;

  splice_dsq = NULL;
  edge       = NULL; 

  upstream_len   = abs(upstream_domain->jali - upstream_domain->iali) + 1;
  downstream_len = abs(downstream_domain->jali - downstream_domain->iali) + 1;
  intron_len     = splice_seq->n - (upstream_len + downstream_len);

  amino_gap = (downstream_domain->ihmm - upstream_domain->jhmm - 1) / 2;
  nuc_extention = ESL_MAX(amino_gap * 3, MAX_AMINO_EXT_SHORT * 3);

  /*Should we align through the full intron length or only enought nucleotides to cove the amino gap */
  if(full_intron) intron_chars = intron_len;
  else            intron_chars   = ESL_MIN(intron_len, nuc_extention*2);
 
  splice_len = upstream_len + downstream_len + intron_chars;
  splice_dsq = NULL;

  if(intron_chars == intron_len) {
    splice_dsq = splice_seq->dsq;
  }
  else {
    ESL_ALLOC(splice_dsq, sizeof(ESL_DSQ) * (splice_len+2));
    splice_dsq[0] = eslDSQ_SENTINEL;
    for(i = 1; i <= upstream_len + nuc_extention; i++) 
      splice_dsq[i] = splice_seq->dsq[i];
    j = downstream_len + nuc_extention - 1;
    for(; i <= splice_len; i++) {
      splice_dsq[i] = splice_seq->dsq[splice_seq->n-j];
      j--;
    }
    splice_dsq[i] = eslDSQ_SENTINEL;
  }
  
  sub_hmm     = p7_splice_GetSubHMM(hmm, upstream_domain->ihmm, downstream_domain->jhmm);
  if(!upstream_domain->tr->fs && !downstream_domain->tr->fs)
    sub_hmm->fs = 0.;
 
  sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
  p7_ProfileConfig_fs(sub_hmm, pli->bg, gcode, sub_fs_model, splice_len, p7_UNIGLOBAL);
  vit_mx = p7_gmx_fs_Create(sub_fs_model->M, splice_len, splice_len, p7P_SPLICE); 

  tr = p7_trace_fs_Create();

  if(!upstream_domain->tr->fs && !downstream_domain->tr->fs) {
    p7_sp_trans_semiglobal_Viterbi(splice_dsq, gcode, splice_len, sub_fs_model, vit_mx);
    p7_sp_trans_semiglobal_VTrace(splice_dsq, splice_len, gcode, sub_fs_model, vit_mx, tr);
  }
  else {
    p7_sp_fs_semiglobal_Viterbi(splice_dsq, gcode, splice_len, sub_fs_model, vit_mx);
    p7_sp_fs_semiglobal_VTrace(splice_dsq, splice_len, gcode, sub_fs_model, vit_mx, tr);
  }


  p7_trace_fs_Index(tr);
  /* If the traces starts after the end of the upstream hit or ends 
   * before the start of the downstream hit this is a failed splice */
  if(tr->sqfrom[0] > upstream_len || tr->sqto[0] < upstream_len+intron_chars+1) {
    p7_hmm_Destroy(sub_hmm);
    p7_profile_fs_Destroy(sub_fs_model);
    p7_gmx_Destroy(vit_mx);
    p7_trace_fs_Destroy(tr);
    if(intron_chars != intron_len) free(splice_dsq);
    return NULL;
  }

  /* Find number of introns in trace */
  intron_cnt = 0;
  for(z = 0; z < tr->N; z++)
    if(tr->st[z] == p7T_R) intron_cnt++;
  // printf("intron_cnt %d\n", intron_cnt);
  /* If no intron states were used this is a failed splice */
  if(intron_cnt == 0) {
    p7_hmm_Destroy(sub_hmm);
    p7_profile_fs_Destroy(sub_fs_model);
    p7_gmx_Destroy(vit_mx);
    p7_trace_fs_Destroy(tr);
    if(intron_chars != intron_len) free(splice_dsq);
    return NULL;
  }
  /* If multiple introns overlap the gap between hits this suggests 
   * there is an intervening exon and no edge should be created */
  if(intron_cnt > 1) 
  {  
    introns_in_gap = 0;
    z2 = 0;
    for(i = 0 ; i < intron_cnt; i++) {
      for(z1 = z2;   z1 < tr->N; z1++) if(tr->st[z1] == p7T_R) break;
      for(z2 = z1+1; z2 < tr->N; z2++) if(tr->st[z2] != p7T_P) break;
      overlap_min = ESL_MAX(tr->i[z1-1], upstream_len+1);
      overlap_max = ESL_MIN(tr->i[z2],   upstream_len+intron_chars);
      overlap_len = overlap_max - overlap_min + 1;
    //  printf("overlap_len %d\n", overlap_len);
      if(overlap_len > 1) 
        introns_in_gap++;
    }
    if(introns_in_gap > 1) {
      p7_hmm_Destroy(sub_hmm);
      p7_profile_fs_Destroy(sub_fs_model);
      p7_gmx_Destroy(vit_mx);
      p7_trace_fs_Destroy(tr);
      if(intron_chars != intron_len) free(splice_dsq);
      return NULL;
    }
  }
  
  /* If more than on set of splice states were used then we need 
   * to find which one is closest to the existing hit bounds */
  total_distance = splice_len*2;
  z2 = 0;
  for(i = 0; i < intron_cnt; i++) {
    for(z1 = z2;   z1 < tr->N; z1++) if(tr->st[z1] == p7T_R) break;
    for(z2 = z1+1; z2 < tr->N; z2++) if(tr->st[z2] == p7T_M) break;
    for(         ; z1 > 0    ; z1--) if(tr->st[z1] == p7T_M) break;
    if(i == 0)
    distance_from_upstream   = abs(tr->i[z1] - upstream_len);
    distance_from_downstream = abs(tr->i[z2] - (upstream_len+intron_chars));

    if(distance_from_upstream+distance_from_downstream < total_distance) {
       total_distance = distance_from_upstream+distance_from_downstream;
       intron_start = z1; 
       intron_end   = z2;
    } 
  }
  
  edge = p7_spliceedge_Create();
  edge->overlap_amino_start = ESL_MAX(upstream_domain->ihmm,   (upstream_domain->ihmm + tr->k[intron_start] - 1) - MIN_AMINO_OVERLAP/2);

  edge->overlap_amino_end   = ESL_MIN(downstream_domain->jhmm, (upstream_domain->ihmm + tr->k[intron_end])     + MIN_AMINO_OVERLAP/2);
   //printf("edge->overlap_amino_start %d edge->overlap_amino_end %d\n", edge->overlap_amino_start, edge->overlap_amino_end); 
  get_overlap_nuc_coords(edge, upstream_domain, downstream_domain, splice_seq, revcomp);
  
  if(upstream_domain->tr->fs || downstream_domain->tr->fs)
    edge->frameshift = TRUE;

  edge->upstream_nuc_end     = ESL_MIN(edge->upstream_nuc_end + 2,     splice_seq->n);
  edge->downstream_nuc_start = ESL_MAX(edge->downstream_nuc_start - 2, 1);

  if ((status = find_optimal_splice_site (edge, pli, upstream_domain, downstream_domain, hmm, gcode, splice_seq)) != eslOK) goto ERROR;

  /* If no splice was found or the splice score plus the minimum hit score is less than zero, try again with a larger overlap */
  if(ESL_MIN(upstream_domain->aliscore, downstream_domain->aliscore) + edge->splice_score < 0) {
	
    /* If the edge failed using the MIN_AMINO_OVERLAP extntion tray again with the MAX_AMINO_OVERLAP */
    edge->overlap_amino_start = ESL_MAX(upstream_domain->ihmm,   (upstream_domain->ihmm + tr->k[intron_start] - 1) - MAX_AMINO_OVERLAP/2);
    edge->overlap_amino_end   = ESL_MIN(downstream_domain->jhmm, (upstream_domain->ihmm + tr->k[intron_end])     + MAX_AMINO_OVERLAP/2);
  
    get_overlap_nuc_coords(edge, upstream_domain, downstream_domain, splice_seq, revcomp);

    edge->upstream_nuc_end     = ESL_MIN(edge->upstream_nuc_end + 2,     splice_seq->n);
    edge->downstream_nuc_start = ESL_MAX(edge->downstream_nuc_start - 2, 1);
  
    if ((status = find_optimal_splice_site (edge, pli, upstream_domain, downstream_domain, hmm, gcode, splice_seq)) != eslOK) goto ERROR;

    if(edge->splice_score == -eslINFINITY) {    
      p7_hmm_Destroy(sub_hmm);
      p7_profile_fs_Destroy(sub_fs_model);
      p7_gmx_Destroy(vit_mx);
      p7_trace_fs_Destroy(tr);
      if(intron_chars != intron_len) free(splice_dsq);
      free(edge);
      return NULL;
    }
  }

  if(revcomp) {
    edge->upstream_spliced_nuc_end     = splice_seq->n - edge->upstream_spliced_nuc_end     + splice_seq->end;
    edge->downstream_spliced_nuc_start = splice_seq->n - edge->downstream_spliced_nuc_start + splice_seq->end;
  }
  else {
    edge->upstream_spliced_nuc_end     = edge->upstream_spliced_nuc_end     + splice_seq->start - 1;
    edge->downstream_spliced_nuc_start = edge->downstream_spliced_nuc_start + splice_seq->start - 1;
  }

  p7_hmm_Destroy(sub_hmm);
  p7_profile_fs_Destroy(sub_fs_model);
  p7_gmx_Destroy(vit_mx);
  p7_trace_fs_Destroy(tr);
  if(intron_chars != intron_len) free(splice_dsq);

  return edge;
   
  
  ERROR:
    p7_hmm_Destroy(sub_hmm);
    p7_profile_fs_Destroy(sub_fs_model);
    p7_gmx_Destroy(vit_mx);
    p7_trace_fs_Destroy(tr);
    if(splice_dsq != NULL) free(splice_dsq);
    if(edge != NULL) free(edge);
    return NULL;

}

SPLICE_EDGE*
p7_spliceedge_ConnectSplits(SPLICE_PIPELINE *pli, const P7_DOMAIN *upstream_domain, const P7_DOMAIN *downstream_domain, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq, int frameshift, int revcomp) 
{

  SPLICE_EDGE   *edge;

  edge = p7_spliceedge_Create();
  edge->overlap_amino_start = ESL_MAX(upstream_domain->ihmm,   upstream_domain->jhmm   - MIN_AMINO_OVERLAP/2);
  edge->overlap_amino_end   = ESL_MIN(downstream_domain->jhmm, downstream_domain->ihmm + MIN_AMINO_OVERLAP/2);
   //printf("up ihmm %d jhmm %d down ihmm %d jhmm %d\n", upstream_domain->ihmm, upstream_domain->jhmm,  downstream_domain->ihmm, downstream_domain->jhmm);
   //printf("up iali %d jali %d down iali %d jali %d\n", upstream_domain->iali, upstream_domain->jali,  downstream_domain->iali, downstream_domain->jali);
   //printf("edge->overlap_amino_start %d edge->overlap_amino_end %d \n", edge->overlap_amino_start, edge->overlap_amino_end);
  get_overlap_nuc_coords(edge, upstream_domain, downstream_domain, splice_seq, revcomp);
  
  edge->upstream_nuc_end     = ESL_MIN(edge->upstream_nuc_end + 2,     splice_seq->n);
  edge->downstream_nuc_start = ESL_MAX(edge->downstream_nuc_start - 2, 1);

  edge->frameshift = frameshift;

  find_optimal_splice_site (edge, pli, upstream_domain, downstream_domain, hmm, gcode, splice_seq);

  if(ESL_MIN(upstream_domain->aliscore, downstream_domain->aliscore) + edge->splice_score < 0) { 

    edge->overlap_amino_start = ESL_MAX(upstream_domain->ihmm,   upstream_domain->jhmm   - MAX_AMINO_OVERLAP/2);
    edge->overlap_amino_end   = ESL_MIN(downstream_domain->jhmm, downstream_domain->ihmm + MAX_AMINO_OVERLAP/2);

    get_overlap_nuc_coords(edge, upstream_domain, downstream_domain, splice_seq, revcomp);

    edge->upstream_nuc_end     = ESL_MIN(edge->upstream_nuc_end + 2,     splice_seq->n);
    edge->downstream_nuc_start = ESL_MAX(edge->downstream_nuc_start - 2, 1);

    find_optimal_splice_site (edge, pli, upstream_domain, downstream_domain, hmm, gcode, splice_seq);

    if(edge->splice_score == -eslINFINITY) {
      free(edge);
      return NULL;
    }
  }

  if(revcomp) {
    edge->upstream_spliced_nuc_end     = splice_seq->n - edge->upstream_spliced_nuc_end     + splice_seq->end;
    edge->downstream_spliced_nuc_start = splice_seq->n - edge->downstream_spliced_nuc_start + splice_seq->end;
  }
  else {
    edge->upstream_spliced_nuc_end     = edge->upstream_spliced_nuc_end     + splice_seq->start - 1;
    edge->downstream_spliced_nuc_start = edge->downstream_spliced_nuc_start + splice_seq->start - 1;
  }

  
  return edge;

}

int
get_overlap_nuc_coords (SPLICE_EDGE *edge, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const ESL_SQ *splice_seq, int revcomp)
{

  int       z1,z2;
  int       strand;
  int       curr_hmm_pos;
  P7_TRACE *up_trace;
  P7_TRACE *down_trace;
  strand = (revcomp? -1 : 1);

  up_trace = upstream->tr;
  down_trace = downstream->tr;

 /***************UPSTREAM*******************/

  for (z2 = up_trace->N-1 ; z2 >= 0; z2--) if (up_trace->st[z2] == p7T_M) break;
  edge->upstream_trace_end = z2;
  curr_hmm_pos = upstream->jhmm;
  edge->upstream_nuc_end   = upstream->jali;
  edge->upstream_nuc_start = upstream->jali;
  
  /* If the amino overlap end is outside the upstream hit find the upstream seq 
   * overlap end by adding the appropriate number of nucs to the hit jali coords */
  if(edge->overlap_amino_end > upstream->jhmm) { 
    edge->upstream_nuc_end += (edge->overlap_amino_end - upstream->jhmm) * strand * 3;
  }
  /* If the amino overlap end is inside the upstream hit 
   * use the trace to find the upstream seq overlap end */
  else {
    while(curr_hmm_pos > edge->overlap_amino_end) {
      if      (up_trace->st[z2] == p7T_M) {
        edge->upstream_nuc_end -= up_trace->c[z2] * strand;
        curr_hmm_pos--;
      }    
      else if (up_trace->st[z2] == p7T_I)
        edge->upstream_nuc_end -= 3 * strand;
      else if (up_trace->st[z2] == p7T_D)
        curr_hmm_pos--;
      else
        ESL_EXCEPTION(eslFAIL, "splice boundry not found in trace");
     
       z2--;
    }    
  }
   
  /* If the amino overlap start is outside the upstream hit find the upstream seq
   * overlap start by adding the appropriate number of nucs to the hit jali coords */
  if(edge->overlap_amino_start > upstream->jhmm) {
    edge->upstream_nuc_start += (edge->overlap_amino_start - upstream->jhmm) * strand * 3 - 2; 
    edge->upstream_trace_start = edge->upstream_trace_end+1;
  }
  /* If the amino overlap start is inside the upstream hit
   * use the trace to find the upstream seq overlap start */
  else {
    edge->upstream_nuc_start = (revcomp ? 
                               ESL_MAX(edge->upstream_nuc_end, edge->upstream_nuc_start) : 
                               ESL_MIN(edge->upstream_nuc_end, edge->upstream_nuc_start));
    
    while(curr_hmm_pos > edge->overlap_amino_start) {
  
      if      (up_trace->st[z2] == p7T_M) {
        edge->upstream_nuc_start -= up_trace->c[z2] * strand;
        //if(up_trace->c[z2] != 3) edge->frameshift = TRUE;
        curr_hmm_pos--;
      }
      else if (up_trace->st[z2] == p7T_I)
        edge->upstream_nuc_start -= 3 * strand;
      else if (up_trace->st[z2] == p7T_D)
        curr_hmm_pos--;
      else
        ESL_EXCEPTION(eslFAIL, "splice boundry not found in trace");
  
      z2--;
    }
    edge->upstream_nuc_start -= (up_trace->c[z2]-1) * strand;
    edge->upstream_trace_start = z2;
  }
//printf("edge->upstream_nuc_start %d edge->upstream_nuc_end %d\n", edge->upstream_nuc_start, edge->upstream_nuc_end); 
  /* Adjust upstream nuc overlap coordrs to be on the splice_seq */
  if(revcomp) {
    edge->upstream_nuc_start = splice_seq->start - edge->upstream_nuc_start + 1;
    edge->upstream_nuc_end   = splice_seq->start - edge->upstream_nuc_end   + 1;
  }
  else {
    edge->upstream_nuc_start = edge->upstream_nuc_start - splice_seq->start + 1;
    edge->upstream_nuc_end   = edge->upstream_nuc_end   - splice_seq->start + 1;
  }


  /***************DOWNSTREAM*******************/

  for (z1 = 0; z1 < down_trace->N; z1++) if (down_trace->st[z1] == p7T_M) break;
  edge->downstream_trace_start = z1;
  curr_hmm_pos = downstream->ihmm;
  edge->downstream_nuc_start = downstream->iali;
  edge->downstream_nuc_end   = downstream->iali; //- ((down_trace->c[z1]-1) * strand);
 
  /* If the amino overlap start is outside the downstream hit find the downstream seq
   * overlap end by subtracting the appropriate number of nucs from the hit iali coords */
  if(edge->overlap_amino_start < downstream->ihmm) {
    edge->downstream_nuc_start -= (downstream->ihmm - edge->overlap_amino_start) * strand * 3; 
  }
  /* If the amino overlap start is inside the downstream hit
   * use the trace to find the downstream seq overlap start */
  else {
    while(curr_hmm_pos < edge->overlap_amino_start) {
      if (down_trace->st[z1] == p7T_M) {
        edge->downstream_nuc_start += down_trace->c[z1] * strand;
        curr_hmm_pos++;
      }
      else if (down_trace->st[z1] == p7T_I)
        edge->downstream_nuc_start += 3 * strand;
      else if (down_trace->st[z1] == p7T_D)
        curr_hmm_pos++;
      else 
        ESL_EXCEPTION(eslFAIL, "splice boundry not found in trace"); 
      
      z1++; 
    }
  }  
  
  /* If the amino overlap end is outside the downstream hit find the downstream seq
   * overlap end by subreacting the appropriate number of nucs to the hit iali coords */
  if(edge->overlap_amino_end < downstream->ihmm) { 
    edge->downstream_nuc_end -= (downstream->ihmm - edge->overlap_amino_end) * strand * 3; 
    edge->downstream_trace_end = edge->downstream_trace_start-1;   
  
  }
  else { 
    edge->downstream_nuc_end = (revcomp ?
                               ESL_MIN(edge->downstream_nuc_start, downstream->iali) :
                               ESL_MAX(edge->downstream_nuc_start, downstream->iali));
   
    while(curr_hmm_pos <= edge->overlap_amino_end) {
      if (down_trace->st[z1] == p7T_M) {
        edge->downstream_nuc_end += down_trace->c[z1] * strand;
        
        //if(down_trace->c[z1] != 3) edge->frameshift = TRUE;
        curr_hmm_pos++;
      }
      else if (down_trace->st[z1] == p7T_I)
        edge->downstream_nuc_end += 3 * strand;
      else if (down_trace->st[z1] == p7T_D)
        curr_hmm_pos++;
      else
        ESL_EXCEPTION(eslFAIL, "splice boundry not found in trace");
      
      z1++;
    }
    edge->downstream_nuc_end += (down_trace->c[z1]-1) * strand; 
   
    edge->downstream_trace_end = z1-1;
  }
  //printf("edge->downstream_nuc_start %d edge->downstream_nuc_end %d\n", edge->downstream_nuc_start, edge->downstream_nuc_end);
  /* Adjust downstream nuc overlap coordrs to be on the splice_seq */ 
  if(revcomp) {
    edge->downstream_nuc_start = splice_seq->start - edge->downstream_nuc_start + 1;
    edge->downstream_nuc_end   = splice_seq->start - edge->downstream_nuc_end   + 1;
  }
  else {
    edge->downstream_nuc_start = edge->downstream_nuc_start - splice_seq->start + 1;
    edge->downstream_nuc_end   = edge->downstream_nuc_end   - splice_seq->start + 1;
  }

 return eslOK;
}


int
find_optimal_splice_site (SPLICE_EDGE *edge, SPLICE_PIPELINE *pli, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq)
{

  int            n, z;
  int            unp,dnp;
  float          lost_sc;
  float         *tp_0;
  float         *tp_M;
  P7_HMM        *sub_hmm;
  P7_PROFILE    *sub_model;
  P7_FS_PROFILE *sub_fs_model;
  int            status;

  
  /* Get the summed ali score for the overlap region covered by the existing alignments */
  lost_sc = 0.;
  n = upstream->tr->tto[0] - upstream->tr->tfrom[0] - 2; 
  //p7_trace_fs_Dump(stdout, upstream->tr,NULL, NULL);
  //printf("edge->upstream_trace_end %d edge->upstream_trace_start %d\n", edge->upstream_trace_end, edge->upstream_trace_start);
  for(z = edge->upstream_trace_end; z >= edge->upstream_trace_start; z--) {
    lost_sc += upstream->scores_per_pos[n];
    n--;
  }
  
  //p7_trace_fs_Dump(stdout, downstream->tr,NULL, NULL);
  //printf("edge->downstream_trace_end %d edge->downstream_trace_start %d\n", edge->downstream_trace_end, edge->downstream_trace_start);
  n = 0;
  for(z = edge->downstream_trace_start; z <= edge->downstream_trace_end; z++) {
    lost_sc += downstream->scores_per_pos[n];
    n++;
  }
   
  /*Get a submodel that covers the overlap region */
  sub_hmm    = p7_splice_GetSubHMM(hmm, edge->overlap_amino_start, edge->overlap_amino_end);
  sub_model = NULL;
  sub_fs_model = NULL;

  if(edge->frameshift) {
    sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
    ESL_REALLOC(sub_fs_model->tsc,  sizeof(float)   * (sub_hmm->M+1) * p7P_NTRANS);
    p7_ProfileConfig_fs(sub_hmm, pli->bg, gcode, sub_fs_model, 100, p7_UNIGLOBAL);
    tp_0 = sub_fs_model->tsc;
    tp_M = sub_fs_model->tsc + sub_fs_model->M * p7P_NTRANS;

    /* Add extra transtion position for insetions at begining and end of alignments*/
    tp_0[p7P_MM] = log(hmm->t[edge->overlap_amino_start-1][p7H_MM]);
    tp_0[p7P_MI] = log(hmm->t[edge->overlap_amino_start-1][p7H_MI]);
    tp_0[p7P_MD] = log(hmm->t[edge->overlap_amino_start-1][p7H_MD]);
    tp_0[p7P_IM] = log(hmm->t[edge->overlap_amino_start-1][p7H_IM]);
    tp_0[p7P_II] = log(hmm->t[edge->overlap_amino_start-1][p7H_II]);
    tp_0[p7P_DM] = log(hmm->t[edge->overlap_amino_start-1][p7H_DM]);
    tp_0[p7P_DD] = log(hmm->t[edge->overlap_amino_start-1][p7H_DD]);


    tp_M[p7P_MM] = log(hmm->t[edge->overlap_amino_end][p7H_MM]);
    tp_M[p7P_MI] = log(hmm->t[edge->overlap_amino_end][p7H_MI]);
    tp_M[p7P_MD] = log(hmm->t[edge->overlap_amino_end][p7H_MD]);
    tp_M[p7P_IM] = log(hmm->t[edge->overlap_amino_end][p7H_IM]);
    tp_M[p7P_II] = log(hmm->t[edge->overlap_amino_end][p7H_II]);
    tp_M[p7P_DM] = log(hmm->t[edge->overlap_amino_end][p7H_DM]);
    tp_M[p7P_DD] = log(hmm->t[edge->overlap_amino_end][p7H_DD]);
  }

  sub_model  = p7_profile_Create (sub_hmm->M, sub_hmm->abc);
  ESL_REALLOC(sub_model->tsc,  sizeof(float)   * (sub_hmm->M+1) * p7P_NTRANS);
  p7_ProfileConfig(sub_hmm, pli->bg, sub_model, 100, p7_UNIGLOBAL);
  tp_0 = sub_model->tsc;
  tp_M = sub_model->tsc + sub_model->M * p7P_NTRANS;

  /* Add extra transtion position for insetions at begining and end of alignments*/
  tp_0[p7P_MM] = log(hmm->t[edge->overlap_amino_start-1][p7H_MM]);
  tp_0[p7P_MI] = log(hmm->t[edge->overlap_amino_start-1][p7H_MI]);
  tp_0[p7P_MD] = log(hmm->t[edge->overlap_amino_start-1][p7H_MD]);
  tp_0[p7P_IM] = log(hmm->t[edge->overlap_amino_start-1][p7H_IM]);
  tp_0[p7P_II] = log(hmm->t[edge->overlap_amino_start-1][p7H_II]);
  tp_0[p7P_DM] = log(hmm->t[edge->overlap_amino_start-1][p7H_DM]);
  tp_0[p7P_DD] = log(hmm->t[edge->overlap_amino_start-1][p7H_DD]);


  tp_M[p7P_MM] = log(hmm->t[edge->overlap_amino_end][p7H_MM]);
  tp_M[p7P_MI] = log(hmm->t[edge->overlap_amino_end][p7H_MI]);
  tp_M[p7P_MD] = log(hmm->t[edge->overlap_amino_end][p7H_MD]);
  tp_M[p7P_IM] = log(hmm->t[edge->overlap_amino_end][p7H_IM]);
  tp_M[p7P_II] = log(hmm->t[edge->overlap_amino_end][p7H_II]);
  tp_M[p7P_DM] = log(hmm->t[edge->overlap_amino_end][p7H_DM]);
  tp_M[p7P_DD] = log(hmm->t[edge->overlap_amino_end][p7H_DD]);

  /* Scan the overlap nucleotides for splice signals */
  for(unp = edge->upstream_nuc_start; unp < edge->upstream_nuc_end; unp++) {
    /*GT-AG*/
    if(splice_seq->dsq[unp] == 2 && splice_seq->dsq[unp+1] == 3) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {
        if(splice_seq->dsq[dnp-1] == 0 && splice_seq->dsq[dnp] == 2) {
          //printf("unp-1 %d dnp+1 %d\n", unp-1+splice_seq->start-1, dnp+1+splice_seq->start-1);
          if (dnp - unp > MIN_INTRON_LEN)
            select_splice_option(edge, sub_model, sub_fs_model, gcode, splice_seq, pli->signal_scores[p7S_GTAG], unp-1, dnp+1);
        }
      }
    }
    /*GC-AG*/
    if(splice_seq->dsq[unp] == 2 && splice_seq->dsq[unp+1] == 1) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {
        if(splice_seq->dsq[dnp-1] == 0 && splice_seq->dsq[dnp] == 2) {
          if (dnp - unp > MIN_INTRON_LEN)
            select_splice_option(edge, sub_model, sub_fs_model, gcode, splice_seq, pli->signal_scores[p7S_GCAG], unp-1, dnp+1);
        }
      }
    }
    /*AT-AC*/
    if(splice_seq->dsq[unp] == 0 && splice_seq->dsq[unp+1] == 3) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {
        if(splice_seq->dsq[dnp-1] == 0 && splice_seq->dsq[dnp] == 1) {
          if (dnp - unp > MIN_INTRON_LEN)
            select_splice_option(edge, sub_model, sub_fs_model, gcode, splice_seq, pli->signal_scores[p7S_ATAC], unp-1, dnp+1);
        }
      }
    }
  }
  //printf("edge->splice_score %f edge->signal_score  %f lost_sc %f\n", edge->splice_score, edge->signal_score, lost_sc);
  if(edge->splice_score != -eslINFINITY)
    edge->splice_score -= lost_sc;
   
  if(sub_hmm       != NULL) p7_hmm_Destroy(sub_hmm);
  if(sub_model     != NULL) p7_profile_Destroy(sub_model);
  if(sub_fs_model  != NULL) p7_profile_fs_Destroy(sub_fs_model);

  return eslOK;

  ERROR:
    if(sub_hmm       != NULL) p7_hmm_Destroy(sub_hmm);
    if(sub_model     != NULL) p7_profile_Destroy(sub_model);
    if(sub_fs_model  != NULL) p7_profile_fs_Destroy(sub_fs_model);
    return status;

}

int
select_splice_option (SPLICE_EDGE *edge, P7_PROFILE *sub_model, P7_FS_PROFILE *sub_fs_model, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq, float signal_score, int up_nuc_pos, int down_nuc_pos)
{

  int         i,k,z;
  int         nuc_seq_idx;
  int         nuc_seq_len;
  int         amino_len;
  int         amino;
  int         upstream_nuc_cnt;
  int         upstream_amino_cnt;
  int         upstream_amino_end;
  float       sum_ali_sc;
  float       vitsc;
  float       overlap_sc;
  ESL_DSQ    *nuc_dsq;
  ESL_DSQ    *amino_dsq;
  P7_GMX     *vit_mx;
  P7_TRACE   *tr;
  int         status;

  nuc_dsq   = NULL;
  amino_dsq = NULL;
  vit_mx    = NULL;
  tr        = NULL;

  /*Get the overlap nucleotides that correspond to the current splice signal*/
  nuc_seq_len  = up_nuc_pos - edge->upstream_nuc_start + 1;
  //printf("nuc_seq_len %d\n", nuc_seq_len);
  nuc_seq_len += edge->downstream_nuc_end - down_nuc_pos + 1;
 //printf("nuc_seq_len %d\n", nuc_seq_len); 

  /* If the splice signal is right at the overlap boundries all amino acid positions are deletions */
  if(nuc_seq_len == 0) {
    sum_ali_sc = 0;
    for(k = 0; k < sub_model->M; k++)
      sum_ali_sc += sub_model->tsc[k * p7P_NTRANS + p7H_DD];

    overlap_sc = sum_ali_sc + signal_score;
    
    if (overlap_sc > edge->splice_score) {
      edge->signal_score = signal_score;
      edge->splice_score = overlap_sc;
      edge->upstream_spliced_nuc_end = up_nuc_pos;
      edge->downstream_spliced_nuc_start = down_nuc_pos;
      edge->upstream_spliced_amino_end = edge->overlap_amino_start - 1;
      edge->downstream_spliced_amino_start = edge->overlap_amino_start;
    }
    return eslOK;
  }


  /* When not using the frameshift algorithm we can only
   * allow mod three length nucelotide sequences  */
  //printf("frameshift %d nuc_seq_len mod 3 %d\n", edge->frameshift, nuc_seq_len % 3); 
  if((!edge->frameshift) && nuc_seq_len % 3) return eslOK;

  ESL_ALLOC(nuc_dsq,   sizeof(ESL_DSQ) * (nuc_seq_len+2));
  nuc_seq_idx = 0;
  nuc_dsq[nuc_seq_idx] = eslDSQ_SENTINEL;
  nuc_seq_idx++;
  for(i = edge->upstream_nuc_start; i <= up_nuc_pos; i++) {
    nuc_dsq[nuc_seq_idx] = splice_seq->dsq[i];
    nuc_seq_idx++;
  }
  upstream_nuc_cnt = nuc_seq_idx-1;
  upstream_amino_cnt = (nuc_seq_idx-1) / 3;
  if((nuc_seq_idx-1) % 3)  upstream_amino_cnt++;

  for (i = down_nuc_pos; i <= edge->downstream_nuc_end; i++) {
    nuc_dsq[nuc_seq_idx] = splice_seq->dsq[i];
    nuc_seq_idx++;
  }
  nuc_dsq[nuc_seq_idx] = eslDSQ_SENTINEL;

  if(nuc_seq_len % 3) {
     vit_mx = p7_gmx_fs_Create(sub_fs_model->M, nuc_seq_len, nuc_seq_len, p7P_CODONS);
     p7_fs_ReconfigLength(sub_fs_model, nuc_seq_len);
     p7_fs_global_Viterbi(nuc_dsq, gcode, nuc_seq_len, sub_fs_model, vit_mx, &vitsc);
  }
  else {
    /* Translate overalp nucleotides to amino sequence */
    amino_len = nuc_seq_len / 3;
    ESL_ALLOC(amino_dsq, sizeof(ESL_DSQ) * (amino_len+2));

    amino_dsq[0] = eslDSQ_SENTINEL;
    nuc_seq_idx = 1;
    for(i = 1; i <= amino_len; i++) {
      amino = esl_gencode_GetTranslation(gcode,&nuc_dsq[nuc_seq_idx]);
      amino_dsq[i] = amino;
      nuc_seq_idx+=3;
    }
    amino_dsq[amino_len+1] = eslDSQ_SENTINEL;

    /* Align translated overlap amino acids to submodel */
    vit_mx = p7_gmx_Create(sub_model->M, amino_len);
    p7_ReconfigLength(sub_model, amino_len);
    p7_global_Viterbi(amino_dsq, amino_len, sub_model, vit_mx, &vitsc);
  }

    /* vitsc from global alignment is the same as ali_score */
  if (vitsc != -eslINFINITY) {
    overlap_sc = vitsc + signal_score;

    if (overlap_sc > edge->splice_score) {
      /*use trace to find splice point*/

      if(nuc_seq_len % 3) {
        tr = p7_trace_fs_Create();
        p7_fs_global_Trace(nuc_dsq, nuc_seq_len, sub_fs_model, vit_mx, tr);

        z = 0;
        while(z < tr->N) {
          if(tr->i[z] >= upstream_nuc_cnt) {
            upstream_amino_end = tr->k[z] + edge->overlap_amino_start - 1;
            z = tr->N;
          }
          z++;
        }
      }
      else {
        tr = p7_trace_Create();
        p7_global_Trace(amino_dsq, amino_len, sub_model, vit_mx, tr);

        z = 0;
        while(z < tr->N) {
          if(tr->i[z] == upstream_amino_cnt) {
            upstream_amino_end = tr->k[z] + edge->overlap_amino_start - 1;
            z = tr->N;
          }
          z++;
        }
      }
      edge->signal_score = signal_score;
      edge->splice_score = overlap_sc;
      edge->upstream_spliced_nuc_end = up_nuc_pos;
      edge->downstream_spliced_nuc_start = down_nuc_pos;
      edge->upstream_spliced_amino_end = upstream_amino_end;
      edge->downstream_spliced_amino_start = upstream_amino_end + 1;
      
    }
  }
  /* Case where all translated aminos are stop codons - treat them as insertions and all hmm positions as deletions */
  else {
    for(i = 0; i < amino_len; i++)
      sum_ali_sc += sub_model->tsc[p7H_II];
    for(k = 0; k < sub_model->M; k++)
      sum_ali_sc += sub_model->tsc[k * p7P_NTRANS + p7H_DD];

    overlap_sc = sum_ali_sc + signal_score;

    if (overlap_sc > edge->splice_score) {
      edge->signal_score = signal_score;
      edge->splice_score = overlap_sc;
      edge->upstream_spliced_nuc_end = up_nuc_pos;
      edge->downstream_spliced_nuc_start = down_nuc_pos;
      edge->upstream_spliced_amino_end = edge->overlap_amino_start - 1;
      edge->downstream_spliced_amino_start = edge->overlap_amino_start;
    }
  }

  if(nuc_dsq   != NULL) free(nuc_dsq);
  if(amino_dsq != NULL) free(amino_dsq);
  if(vit_mx    != NULL) p7_gmx_Destroy(vit_mx);
  if(nuc_seq_len % 3) { if(tr != NULL) p7_trace_fs_Destroy(tr); }
  else                 { if(tr != NULL) p7_trace_Destroy(tr); }
  return eslOK;

  ERROR:
    if(nuc_dsq   != NULL) free(nuc_dsq);
    if(amino_dsq != NULL) free(amino_dsq);
    if(vit_mx    != NULL) p7_gmx_Destroy(vit_mx);
    if(nuc_seq_len % 3) { if(tr != NULL) p7_trace_fs_Destroy(tr); }
    else                 { if(tr != NULL) p7_trace_Destroy(tr); }
    return status;
}

