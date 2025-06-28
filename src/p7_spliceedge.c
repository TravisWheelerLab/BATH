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

static int get_overlap_nuc_coords (SPLICE_EDGE *edge, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const ESL_SQ *splice_seq, const ESL_GENCODE *gcode, int revcomp);
static int find_optimal_splice_site (SPLICE_EDGE *edge, SPLICE_PIPELINE *pli, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream,  const P7_FS_PROFILE *gm_fs, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq);
static int select_splice_option (SPLICE_EDGE *edge, P7_FS_PROFILE *sub_fs_model, const P7_FS_PROFILE *gm_fs, const ESL_GENCODE *gcode, const ESL_SQ *target_seq, float signal_score, int up_nuc_pos, int down_nuc_pos);

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
p7_spliceedge_ConnectNodes(SPLICE_PIPELINE *pli, const P7_DOMAIN *upstream_domain, const P7_DOMAIN *downstream_domain, const P7_FS_PROFILE *gm_fs, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq, int revcomp) 
{

  SPLICE_EDGE   *edge;

  edge = p7_spliceedge_Create();
  p7_spliceedge_AliScoreEdge(edge, upstream_domain, downstream_domain);

  if(upstream_domain->tr->fs || downstream_domain->tr->fs) edge->frameshift = TRUE;

//printf("edge->overlap_amino_start %d edge->overlap_amino_end %d\n", edge->overlap_amino_start, edge->overlap_amino_end);  
  get_overlap_nuc_coords(edge, upstream_domain, downstream_domain, splice_seq, gcode, revcomp);
  
  edge->upstream_nuc_end     = ESL_MIN(edge->upstream_nuc_end + 2,     splice_seq->n);
  edge->downstream_nuc_start = ESL_MAX(edge->downstream_nuc_start - 2, 1);

 
  find_optimal_splice_site (edge, pli, upstream_domain, downstream_domain, gm_fs, hmm, gcode, splice_seq);

  if(ESL_MIN(upstream_domain->aliscore, downstream_domain->aliscore) + edge->splice_score < 0) { 
    //printf("edge->overlap_amino_start %d edge->overlap_amino_end %d\n", edge->overlap_amino_start, edge->overlap_amino_end);
// printf("edge->upstream_nuc_start %d edge->upstream_nuc_end %d\n", edge->upstream_nuc_start + splice_seq->start - 1, edge->upstream_nuc_end + splice_seq->start - 1);
// printf("edge->downstream_nuc_start %d edge->downstream_nuc_end %d\n", edge->downstream_nuc_start + splice_seq->start - 1, edge->downstream_nuc_end + splice_seq->start - 1);
//    printf("edge->splice_score %f \n", edge->splice_score);
    edge->overlap_amino_start = ESL_MAX(upstream_domain->ihmm+1,   edge->overlap_amino_start - MIN_AMINO_OVERLAP/2);
    edge->overlap_amino_end   = ESL_MIN(downstream_domain->jhmm-1, edge->overlap_amino_end   + MIN_AMINO_OVERLAP/2);

    get_overlap_nuc_coords(edge, upstream_domain, downstream_domain, splice_seq, gcode, revcomp);

    edge->upstream_nuc_end     = ESL_MIN(edge->upstream_nuc_end + 2,     splice_seq->n);
    edge->downstream_nuc_start = ESL_MAX(edge->downstream_nuc_start - 2, 1);

    find_optimal_splice_site (edge, pli, upstream_domain, downstream_domain, gm_fs, hmm, gcode, splice_seq);

    if(edge->splice_score == -eslINFINITY) {
      free(edge);
      return NULL;
    }
  }

//printf("edge->splice_score %f \n", edge->splice_score); 
//printf("edge->overlap_amino_start %d edge->overlap_amino_end %d\n", edge->overlap_amino_start, edge->overlap_amino_end);
//printf("edge->upstream_nuc_start %d edge->upstream_nuc_end %d\n", edge->upstream_nuc_start + splice_seq->start - 1, edge->upstream_nuc_end + splice_seq->start - 1);
//printf("edge->downstream_nuc_start %d edge->downstream_nuc_end %d\n", edge->downstream_nuc_start + splice_seq->start - 1, edge->downstream_nuc_end + splice_seq->start - 1);
  /* Convert splice coords to full sequence */
  if(revcomp) {
    edge->upstream_spliced_nuc_end     = splice_seq->n - edge->upstream_spliced_nuc_end     + splice_seq->end;
    edge->downstream_spliced_nuc_start = splice_seq->n - edge->downstream_spliced_nuc_start + splice_seq->end;
  }
  else {
    edge->upstream_spliced_nuc_end     = edge->upstream_spliced_nuc_end     + splice_seq->start - 1;
    edge->downstream_spliced_nuc_start = edge->downstream_spliced_nuc_start + splice_seq->start - 1;
  }
//  printf("edge->upstream_spliced_nuc_end %d edge->downstream_spliced_nuc_start %d\n", edge->upstream_spliced_nuc_end, edge->downstream_spliced_nuc_start);
  return edge;

}


int
p7_spliceedge_AliScoreEdge(SPLICE_EDGE *edge, const P7_DOMAIN *upstream_dom, const P7_DOMAIN *downstream_dom)
{

  int k, p, s, z;
  int us, ds;
  int k_start, k_end;
  int up_N, down_N;
  float lost_sc;
  float min_lost_sc;
  float *spp;
  float *upstream_suffix_sum;
  float *downstream_prefix_sum;
  P7_TRACE    *tr;
  int status;

  upstream_suffix_sum   = NULL;
  downstream_prefix_sum = NULL;

  edge->overlap_amino_start = ESL_MAX(upstream_dom->ihmm+1,   upstream_dom->jhmm   - MIN_AMINO_OVERLAP/2);
  edge->overlap_amino_end   = ESL_MIN(downstream_dom->jhmm-1, downstream_dom->ihmm + MIN_AMINO_OVERLAP/2);

  if(downstream_dom->ihmm <= upstream_dom->jhmm) {

    up_N = upstream_dom->jhmm - upstream_dom->ihmm + 1;
    ESL_ALLOC(upstream_suffix_sum, sizeof(float) * up_N);
    esl_vec_FSet(upstream_suffix_sum, up_N, 0.0);

    tr  = upstream_dom->tr;
    spp = upstream_dom->scores_per_pos;

    z = tr->tfrom[0]+1;
    while(tr->st[z] != p7T_M) z++;
    k = tr->k[z];
    s = 0;
    p = 0;
    while (k <= upstream_dom->jhmm) {
      if(tr->k[z] == k) {
        upstream_suffix_sum[s] += spp[p];
        z++; p++;
      }
      else {
        k++; s++;
      }
    }

    for(s = up_N-2; s >= 0; s--)
      upstream_suffix_sum[s] += upstream_suffix_sum[s+1];

    down_N  = downstream_dom->jhmm - downstream_dom->ihmm + 1;
    ESL_ALLOC(downstream_prefix_sum, sizeof(float) * down_N);
    
    esl_vec_FSet(downstream_prefix_sum, down_N, 0.0);

    tr  = downstream_dom->tr;
    spp = downstream_dom->scores_per_pos;

    z = tr->tfrom[0]+1;
    while(tr->st[z] != p7T_M) z++;
    k = tr->k[z];
    s = 0;
    p = 0;
    while (k <= downstream_dom->jhmm) {
       if(tr->k[z] == k) {
        downstream_prefix_sum[s] += spp[p];
        z++; p++;
      }
      else {
        k++; s++;
      }
    }

    for(s = 1; s < down_N; s++)
      downstream_prefix_sum[s] += downstream_prefix_sum[s-1];

    k_start = ESL_MAX(downstream_dom->ihmm, upstream_dom->ihmm+1);
    k_end   = ESL_MIN(upstream_dom->jhmm,   downstream_dom->jhmm-1);

    /* Find the minimum lost_sc be checking each k from k_start to e_end */
    us = k_start - upstream_dom->ihmm;
    ds = k_start - downstream_dom->ihmm - 1;
    min_lost_sc = upstream_suffix_sum[us];  // lost_sc if all k belong to downstream
    if(ds >= 0)
      min_lost_sc += downstream_prefix_sum[ds];  // include any lost_sc from downstream overextention
    edge->overlap_amino_start = k_start;
    edge->overlap_amino_end   = k_start + 1;

    us++;
    ds++;
    for(k = k_start; k < k_end; k++) {
      /* lost_sc if k_start to k belong to upstream and k+1 to k_end belong to downstream */
      lost_sc = upstream_suffix_sum[us] +  downstream_prefix_sum[ds];
      if(lost_sc < min_lost_sc) {
        min_lost_sc = lost_sc;
        edge->overlap_amino_start = k;
        edge->overlap_amino_end   = k+1;
      }
      us++;
      ds++;
    }

    lost_sc = downstream_prefix_sum[ds]; // lost_sc if all k belong to upstream
    if(k_end < upstream_dom->jhmm)
      lost_sc += upstream_suffix_sum[us]; // include any lost_sc from upstream overextention
    if(lost_sc < min_lost_sc) {
      min_lost_sc = lost_sc;
      edge->overlap_amino_start = k_end;
      edge->overlap_amino_end   = k_end+1;
    }

    edge->splice_score -= min_lost_sc;

    edge->overlap_amino_start = ESL_MAX(upstream_dom->ihmm+1,   edge->overlap_amino_start - MIN_AMINO_OVERLAP/2);
    edge->overlap_amino_end   = ESL_MIN(downstream_dom->jhmm-1, edge->overlap_amino_end   + MIN_AMINO_OVERLAP/2);
    
    if(upstream_suffix_sum   != NULL) free(upstream_suffix_sum);
    if(downstream_prefix_sum != NULL) free(downstream_prefix_sum);

  }

  return eslOK;

  ERROR:
    if(upstream_suffix_sum   != NULL) free(upstream_suffix_sum);
    if(downstream_prefix_sum != NULL) free(downstream_prefix_sum);
    return status;

}


int
get_overlap_nuc_coords (SPLICE_EDGE *edge, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const ESL_SQ *splice_seq, const ESL_GENCODE *gcode, int revcomp)
{

  int       i, x;
  int       z1,z2;
  int       codon;
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
  
  /* If the amino overlap end is outside the upstream hit find the upstream seq 
   * overlap end by adding the appropriate number of nucs to the hit jali coords */
  if(edge->overlap_amino_end > upstream->jhmm) { 
    edge->upstream_nuc_end += (edge->overlap_amino_end - upstream->jhmm) * strand * 3;
    /* If this is a non-frameshift edge - make sure we don't add stop codons */
    if(!edge->frameshift) {
      if(revcomp) {
        for(i = upstream->jali; i >= edge->upstream_nuc_end; i-=3) {
          /* get sub-seq coords */
          x = splice_seq->start - i + 1;
          codon = 16 * splice_seq->dsq[x-2] + 4 * splice_seq->dsq[x-1] + splice_seq->dsq[x]; 
          /*If we have found a stop codon end the upstream sequence at the previous codon */
		  if(gcode->basic[codon] == gcode->aa_abc->Kp-2) {
            edge->upstream_nuc_end = i+3;
            break;     
          }
        }
      }
      else {

        for(i = upstream->jali; i <= edge->upstream_nuc_end; i+=3) {
          /* get sub-seq coords */
          x = i - splice_seq->start + 1;
          codon = 16 * splice_seq->dsq[x-2] + 4 * splice_seq->dsq[x-1] + splice_seq->dsq[x]; 
          /*If we have found a stop codon end the upstream sequence at the previous codon */
		  if(gcode->basic[codon] == gcode->aa_abc->Kp-2) {
            edge->upstream_nuc_end = i-3;
            break;     
          }
        }
      }
    }
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
   
 
  /* Use the trace to find the upstream seq overlap start */
  edge->upstream_nuc_start = (revcomp ? 
                             ESL_MAX(edge->upstream_nuc_end, upstream->jali) : 
                             ESL_MIN(edge->upstream_nuc_end, upstream->jali));
  
  while(curr_hmm_pos > edge->overlap_amino_start) {

    if      (up_trace->st[z2] == p7T_M) {
      edge->upstream_nuc_start -= up_trace->c[z2] * strand;
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
  edge->upstream_nuc_start += strand;
  edge->upstream_trace_start = z2;
  
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
 
  /* If the amino overlap start is outside the downstream hit find the downstream seq
   * overlap end by subtracting the appropriate number of nucs from the hit iali coords */
  if(edge->overlap_amino_start < downstream->ihmm) {
    edge->downstream_nuc_start -= (downstream->ihmm - edge->overlap_amino_start) * strand * 3; 
    /* If this is a non-frameshift edge - make sure we don't add stop codons */
    if(!edge->frameshift) {
      if(revcomp) {
        for(i = downstream->iali; i <= edge->downstream_nuc_start; i+=3) {
          /* get sub-seq coords */
          x = splice_seq->start - i + 1;
          codon = 16 * splice_seq->dsq[x] + 4 * splice_seq->dsq[x+1] + splice_seq->dsq[x+2]; 
          /*If we have found a stop codon end the upstream sequence at the previous codon */
		  if(gcode->basic[codon] == gcode->aa_abc->Kp-2) {
            edge->downstream_nuc_start = i-3;
            break;     
          }
        }
      }
      else {

        for(i = downstream->iali; i >= edge->downstream_nuc_start; i-=3) {
          /* get sub-seq coords */
          x = i - splice_seq->start + 1;
          codon = 16 * splice_seq->dsq[x] + 4 * splice_seq->dsq[x+1] + splice_seq->dsq[x+2]; 
          //printf("i %d codon %d%d%d amino %d\n", i, splice_seq->dsq[x], splice_seq->dsq[x+1], splice_seq->dsq[x+2], gcode->basic[codon]);
          /*If we have found a stop codon start the downstream sequence at the next codon */
		  if(gcode->basic[codon] == gcode->aa_abc->Kp-2) {
            edge->downstream_nuc_start = i+3;
            break;     
          }
        }
      }
    }

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
 
  /* Use the trace to find the downstream seq overlap end */ 
  edge->downstream_nuc_end = (revcomp ?
                             ESL_MIN(edge->downstream_nuc_start, downstream->iali) :
                             ESL_MAX(edge->downstream_nuc_start, downstream->iali));

  
  while(curr_hmm_pos <= edge->overlap_amino_end) {
    if (down_trace->st[z1] == p7T_M) {
      edge->downstream_nuc_end += down_trace->c[z1] * strand;
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
  edge->downstream_nuc_end -= strand;  
  edge->downstream_trace_end = z1-1;
 
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
find_optimal_splice_site (SPLICE_EDGE *edge, SPLICE_PIPELINE *pli, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const P7_FS_PROFILE *gm_fs, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq)
{

  int            n, z;
  int            unp,dnp;
  float          lost_sc;
  float         *tp_0;
  float         *tp_M;
  P7_HMM        *sub_hmm;
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
  if(!edge->frameshift) sub_hmm->fs = 0.;

  sub_fs_model = NULL;
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


  /* Scan the overlap nucleotides for splice signals */
  for(unp = edge->upstream_nuc_start; unp < edge->upstream_nuc_end; unp++) {
    /*GT-AG*/
    if(splice_seq->dsq[unp] == 2 && splice_seq->dsq[unp+1] == 3) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {
        if(splice_seq->dsq[dnp-1] == 0 && splice_seq->dsq[dnp] == 2) {
          
          if (dnp - unp > MIN_INTRON_LENG)
            select_splice_option(edge, sub_fs_model, gm_fs, gcode, splice_seq, pli->signal_scores[p7S_GTAG], unp-1, dnp+1);
//           printf("unp-1 %d dnp+1 %d\n", splice_seq->start-unp, splice_seq->start-dnp+2);
        }
      }
    }
    /*GC-AG*/
    if(splice_seq->dsq[unp] == 2 && splice_seq->dsq[unp+1] == 1) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {
        if(splice_seq->dsq[dnp-1] == 0 && splice_seq->dsq[dnp] == 2) {
          if (dnp - unp > MIN_INTRON_LENG)
            select_splice_option(edge, sub_fs_model, gm_fs, gcode, splice_seq, pli->signal_scores[p7S_GCAG], unp-1, dnp+1);
        }
      }
    }
    /*AT-AC*/
    if(splice_seq->dsq[unp] == 0 && splice_seq->dsq[unp+1] == 3) {
      for(dnp = edge->downstream_nuc_end; dnp > edge->downstream_nuc_start; dnp--) {
        if(splice_seq->dsq[dnp-1] == 0 && splice_seq->dsq[dnp] == 1) {
          if (dnp - unp > MIN_INTRON_LENG)
            select_splice_option(edge, sub_fs_model, gm_fs, gcode, splice_seq, pli->signal_scores[p7S_ATAC], unp-1, dnp+1);
        }
      }
    }
  }

  //printf("edge->splice_score %f edge->signal_score  %f lost_sc %f\n", edge->splice_score, edge->signal_score, lost_sc);
  if(edge->splice_score != -eslINFINITY)
    edge->splice_score -= lost_sc;
   
  if(sub_hmm       != NULL) p7_hmm_Destroy(sub_hmm);
  if(sub_fs_model  != NULL) p7_profile_fs_Destroy(sub_fs_model);

  return eslOK;

  ERROR:
    if(sub_hmm       != NULL) p7_hmm_Destroy(sub_hmm);
    if(sub_fs_model  != NULL) p7_profile_fs_Destroy(sub_fs_model);
    return status;

}

int
select_splice_option (SPLICE_EDGE *edge, P7_FS_PROFILE *sub_fs_model, const P7_FS_PROFILE *gm_fs, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq, float signal_score, int up_nuc_pos, int down_nuc_pos)
{

  int         i,k,z;
  int         nuc_seq_idx;
  int         nuc_seq_len;
  int         amino;
  int         upstream_nuc_cnt;
  int         upstream_amino_cnt;
  int         upstream_amino_end;
  float       sum_ali_sc;
  float       vitsc;
  float       overlap_sc;
  ESL_DSQ    *nuc_dsq;
  ESL_SQ     *nuc_sq;
  P7_GMX     *vit_mx;
  P7_TRACE   *tr;
  P7_DOMAIN  *tmp_dom;
  int         status;

  nuc_dsq   = NULL;
  nuc_sq    = NULL;
  vit_mx    = NULL;
  tr        = NULL;
  tmp_dom   = NULL;

  /*Get the overlap nucleotides that correspond to the current splice signal*/
  nuc_seq_len  = up_nuc_pos - edge->upstream_nuc_start + 1;
  //printf("nuc_seq_len %d\n", nuc_seq_len);
  nuc_seq_len += edge->downstream_nuc_end - down_nuc_pos + 1;
 //printf("nuc_seq_len %d\n", nuc_seq_len); 

  /* If the splice signal is right at the overlap boundries all amino acid positions are deletions */
  if(nuc_seq_len == 0) {
    sum_ali_sc = 0;
    for(k = 0; k < sub_fs_model->M; k++)
      sum_ali_sc += sub_fs_model->tsc[k * p7P_NTRANS + p7H_DD];

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

  if(!edge->frameshift) {
     /* Check for stop codons */
    for(i = 1; i <= nuc_seq_len; i+=3) {
      amino = esl_gencode_GetTranslation(gcode,&nuc_dsq[i]);
      if(amino ==  gcode->aa_abc->Kp-2) {
        if(nuc_dsq   != NULL) free(nuc_dsq);
        return eslOK;
      }
    }
  }

  nuc_sq = esl_sq_CreateDigitalFrom(gcode->nt_abc, NULL, nuc_dsq, nuc_seq_len, NULL,NULL,NULL);

  vit_mx = p7_gmx_fs_Create(sub_fs_model->M, nuc_seq_len, nuc_seq_len, p7P_CODONS);
  p7_fs_ReconfigLength(sub_fs_model, nuc_seq_len);
  p7_fs_global_Viterbi(nuc_dsq, gcode, nuc_seq_len, sub_fs_model, vit_mx, &vitsc);

  tr = p7_trace_fs_Create();
  p7_fs_global_Trace(nuc_dsq, nuc_seq_len, sub_fs_model, vit_mx, tr);

  for(k = 0; k < tr->N; k++)
    if(tr->k[k] > 0) tr->k[k] += edge->overlap_amino_start-1;

  tmp_dom = p7_domain_Create_empty();
  tmp_dom->scores_per_pos = NULL;
  p7_splice_ComputeAliScores_fs(tmp_dom, tr, nuc_sq, gm_fs, NULL, FALSE);

  overlap_sc = tmp_dom->aliscore + signal_score;

  if (overlap_sc > edge->splice_score) {
    /*use trace to find splice point*/

    if(upstream_nuc_cnt == 0) {
      upstream_amino_end = tr->hmmfrom[0]-1;
    }
    else {
      z = 0;
      while(z < tr->N) {
        if(tr->i[z] >= upstream_nuc_cnt) {
          upstream_amino_end = tr->k[z];
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
  
  if(nuc_sq  != NULL) esl_sq_Destroy(nuc_sq);
  if(nuc_dsq != NULL) free(nuc_dsq);
  if(vit_mx  != NULL) p7_gmx_Destroy(vit_mx);
  if(tr      != NULL) p7_trace_fs_Destroy(tr); 
  if(tmp_dom != NULL) p7_domain_Destroy(tmp_dom);
  return eslOK;

  ERROR:
    if(nuc_sq  != NULL) esl_sq_Destroy(nuc_sq);
    if(nuc_dsq != NULL) free(nuc_dsq);
    if(nuc_dsq != NULL) free(nuc_dsq);
    if(vit_mx  != NULL) p7_gmx_Destroy(vit_mx);
    if(tr      != NULL) p7_trace_fs_Destroy(tr); 
    if(tmp_dom != NULL) p7_domain_Destroy(tmp_dom);
    return status;
}

