/* SPLICE_GAP: gap structure defines region in graph where nodes (hits) may be missing 
 *
 * Contents:
 *    1. The SPLICE_GAP object.
 *    2. Gap filling algorithim 
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

/*****************************************************************
 * 1. The SPLICE_GAP structure.
 *****************************************************************/


/* Function:  p7_splicegap_Create()
 * Synopsis:  Allocates a splice gap.
 *
 * Purpose:   Allocates a new <SPLICE_GAP>.
 *
 * Returns:   a pointer to the new <SPLICE_GAP> structure
 *            on success.
 *
 * Throws:    <NULL> on allocation error.
 */
SPLICE_GAP*
p7_splicegap_Create(void)
{

  int status;
  SPLICE_GAP *gap;

  gap = NULL;
  ESL_ALLOC(gap, sizeof(SPLICE_GAP));

  return gap;

  ERROR:
    if (gap != NULL) free(gap);
    return NULL;
}

P7_HIT**
p7_splicegap_AlignGap(SPLICE_GRAPH *graph, SPLICE_GAP *gap, const P7_HMM *hmm, const P7_BG *bg, const ESL_GENCODE *gcode, const ESL_SQFILE *seq_file, int *num_hits) 
{

  int           i, h, y, z;
  int           z1, z2;
  int           hit_cnt;
  int           intron_cnt;
  int           start_new;
  int           duplicate;
  int           ihmm, jhmm;
  int           iali, jali;
  int           new_hit_min, new_hit_max;
  int           new_hit_len;
  int           old_hit_min, old_hit_max;
  int           old_hit_len;
  int           hmm_overlap_min, hmm_overlap_max;
  int           hmm_overlap_len;
  int           seq_overlap_min, seq_overlap_max;
  int           seq_overlap_len;          
  int           codon;
  P7_TOPHITS    *th;
  P7_HMM        *sub_hmm;
  P7_FS_PROFILE *sub_fs_model;
  P7_GMX        *vit_mx;
  P7_TRACE      *tr;
  P7_HIT        *new_hit;
  P7_HIT        **gap_hits;
  ESL_SQ        *gap_seq;
  int           status;

  gap_hits = NULL;
  hit_cnt  = 0;
  th = graph->th;

  /* Fetch the gap subsequence */ 
  gap_seq = p7_splice_GetSubSequence(seq_file, graph->seqname, gap->seq_min, gap->seq_max, graph->revcomp);
   
  /* Fetch the gap sub-hmm and convert to codon model */
  sub_hmm = p7_splice_GetSubHMM(hmm, gap->hmm_min, gap->hmm_max);
  sub_hmm->fs = 0.; //set frameshift probability to 0 to align standard codons only

  sub_fs_model = p7_profile_fs_Create(sub_hmm->M, sub_hmm->abc);
  p7_ProfileConfig_fs(sub_hmm, bg, gcode, sub_fs_model, gap_seq->n*3, p7_GLOBAL); //len*3 to get single nucleotide loops

  vit_mx = p7_gmx_fs_Create(sub_fs_model->M, gap_seq->n, gap_seq->n, p7P_SPLICE);
  tr = p7_trace_fs_Create();

  /* Align the gap using translated viterbi - global with respect to sub model */
  p7_sp_trans_semiglobal_Viterbi(gap_seq->dsq, gcode, gap_seq->n, sub_fs_model, vit_mx);
  p7_sp_trans_semiglobal_VTrace(gap_seq->dsq, gap_seq->n, gcode, sub_fs_model, vit_mx, tr);

   /* Find number of introns in trace */
  intron_cnt = 0;
  for(z = 0; z < tr->N; z++)
    if(tr->st[z] == p7T_R) intron_cnt++;
  
  ESL_ALLOC(gap_hits, sizeof(P7_HIT*) * (intron_cnt+1));
 
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
      ihmm = tr->k[y] + gap->hmm_min - 1;
      jhmm = tr->k[z] + gap->hmm_min - 1;

      if(graph->revcomp) {
        iali = gap_seq->n - tr->i[y] + gap_seq->end + 2;
        jali = gap_seq->n - tr->i[z] + gap_seq->end;
      }
      else {
        iali = gap_seq->start + tr->i[y] - 3;
        jali = gap_seq->start + tr->i[z] - 1;
      }

      /* check if this new hit is a duplicate */
      duplicate = FALSE;

      new_hit_min = ESL_MIN(iali, jali);
      new_hit_max = ESL_MAX(iali, jali);
      new_hit_len = new_hit_max - new_hit_min + 1;

      for (h = 0 ; h < th->N; h++) {
        /* If hit is not in gap skip it */
        if(th->hit[h]->dcl->jhmm < gap->hmm_min || th->hit[h]->dcl->ihmm > gap->hmm_max) continue;

        hmm_overlap_min = ESL_MAX(ihmm, th->hit[h]->dcl->ihmm);
        hmm_overlap_max = ESL_MIN(jhmm, th->hit[h]->dcl->jhmm);
        hmm_overlap_len = hmm_overlap_max - hmm_overlap_min + 1;

        if(hmm_overlap_len < 1) continue;

        old_hit_min = ESL_MIN(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
        old_hit_max = ESL_MAX(th->hit[h]->dcl->iali, th->hit[h]->dcl->jali);
        old_hit_len = old_hit_max - old_hit_min + 1;

        seq_overlap_min = ESL_MAX(new_hit_min, old_hit_min);
        seq_overlap_max = ESL_MIN(new_hit_max, old_hit_max);
        seq_overlap_len = seq_overlap_max - seq_overlap_min + 1;

        if(seq_overlap_len >= new_hit_len * 0.95 ||
           seq_overlap_len >= old_hit_len * 0.95) {
          duplicate = TRUE;
          break;
        }
      }
      /* If hit is duplicate - skip it */
      if(duplicate) { z++; start_new = FALSE; continue; }

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
      for(i = y; i <= z; i++) {
        codon = (tr->st[i] == p7T_M ? 3 : 0);
        p7_trace_fs_Append(new_hit->dcl->tr, tr->st[i], tr->k[i], tr->i[i], codon);
      }

      /* Append ending special states */
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_E, tr->k[z], tr->i[z], 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_C, 0, tr->i[z], 0);
      p7_trace_fs_Append(new_hit->dcl->tr, p7T_T, 0, 0, 0);

      p7_splice_ComputeAliScores_fs(new_hit->dcl, new_hit->dcl->tr, gap_seq->dsq, sub_fs_model, gap_seq->abc);

      /* Adjust k postions to full HMM */
      for(i = 0; i < new_hit->dcl->tr->N; i++) {
        if(new_hit->dcl->tr->k[i] != 0)
          new_hit->dcl->tr->k[i] += gap->hmm_min - 1;
      }
      gap_hits[hit_cnt] = new_hit;
      hit_cnt++;

      start_new = FALSE;
    }
     
    z++;
    if(tr->st[z] == p7T_M) start_new = TRUE;
  } 

  esl_sq_Destroy(gap_seq);
  p7_hmm_Destroy(sub_hmm);
  p7_profile_fs_Destroy(sub_fs_model);
  p7_gmx_Destroy(vit_mx);
  p7_trace_fs_Destroy(tr);

  *num_hits = hit_cnt;
  return gap_hits;

  ERROR:
    esl_sq_Destroy(gap_seq);
    p7_hmm_Destroy(sub_hmm);
    p7_profile_fs_Destroy(sub_fs_model);
    p7_gmx_Destroy(vit_mx);
    p7_trace_fs_Destroy(tr);
    for(h = 0; h < hit_cnt; h++) {
      p7_trace_fs_Destroy(gap_hits[h]->dcl->tr);
      free(gap_hits[h]->dcl->scores_per_pos);
      p7_hit_Destroy(gap_hits[h]);
    }  
    if(gap_hits != NULL) free(gap_hits);
    return NULL;

}
