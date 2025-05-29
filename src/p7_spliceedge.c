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
p7_spliceedge_ConnectHits(SPLICE_PIPELINE *pli, const P7_DOMAIN *upstream_domain, const P7_DOMAIN *downstream_domain, const P7_HMM *hmm, const ESL_GENCODE *gcode, const ESL_SQ *splice_seq, int revcomp) {

  int          up_amino_start, up_amino_end;
  int          down_amino_start, down_amino_end;
  int          overlap;
  int          num_ext_aminos;
  SPLICE_EDGE *edge1;
  SPLICE_EDGE *edge2;
  int          status;
  
  edge1 = NULL;
  edge2 = NULL;

  up_amino_start = upstream_domain->ihmm;
  up_amino_end   = upstream_domain->jhmm;

  down_amino_start = downstream_domain->ihmm;
  down_amino_end   = downstream_domain->jhmm;
  
  /* Ensure overlap region is at least MIN_AMINO_OVERLAP long plus any gap extention positions */
  overlap = up_amino_end - down_amino_start + 1;
  if(overlap >= MIN_AMINO_OVERLAP)  // hits naturally have sufficent overlap
    num_ext_aminos = 0;
  else if( overlap > 0)             //hits overlap nuterally bu still need to be extended
    num_ext_aminos = MIN_AMINO_OVERLAP - overlap;
  else                              // hits do not overlap - gap extention + MIN_AMINO_OVERLAP needed
    num_ext_aminos = (down_amino_start - upstream_domain->jhmm)*2 + MIN_AMINO_OVERLAP;

  
  /* If there is a natural overlap of <= MAX_MAINO_OVERLAP we will search the enrire overlap*/  
  if(overlap <= MAX_AMINO_OVERLAP) {

    edge1 = p7_spliceedge_Create();
    edge1->overlap_amino_start = down_amino_start;
    edge1->overlap_amino_end   = up_amino_end;
  
    /* If the hits do not overlap by at least MIN_AMINO_OVERLAP hmm positions, extend them */
    if (num_ext_aminos > 0) {
      num_ext_aminos = (num_ext_aminos+1)/2;
      edge1->overlap_amino_start -= num_ext_aminos;
      edge1->overlap_amino_end   += num_ext_aminos;
    }
  
    if(edge1->overlap_amino_end > down_amino_end)   edge1->overlap_amino_end = down_amino_end;
    if(edge1->overlap_amino_start < up_amino_start) edge1->overlap_amino_start = up_amino_start;

    get_overlap_nuc_coords(edge1, upstream_domain, downstream_domain, splice_seq, revcomp);
  
    /* Add extra nucleotides for splice sites */
    edge1->upstream_nuc_end     = ESL_MIN(edge1->upstream_nuc_end + 2, splice_seq->n);
    edge1->downstream_nuc_start = ESL_MAX(edge1->downstream_nuc_start - 2, 1);

    if ((status = find_optimal_splice_site (edge1, pli, upstream_domain, downstream_domain, hmm, gcode, splice_seq)) != eslOK) goto ERROR;

    if(edge1->splice_score == -eslINFINITY) {
      free(edge1);
      return NULL;
    }

    if(revcomp) {
      edge1->upstream_spliced_nuc_end     = splice_seq->n - edge1->upstream_spliced_nuc_end     + splice_seq->end;
      edge1->downstream_spliced_nuc_start = splice_seq->n - edge1->downstream_spliced_nuc_start + splice_seq->end;
    }
    else {
      edge1->upstream_spliced_nuc_end     = edge1->upstream_spliced_nuc_end     + splice_seq->start - 1;
      edge1->downstream_spliced_nuc_start = edge1->downstream_spliced_nuc_start + splice_seq->start - 1;
    }

    return edge1;

  }
  /* If the natural overalp is > MAX_MAINO_OVERLAP we will serach the first 
   * and last MAX_MAINO_OVERLAP/2 sections of the total overlap region */
  else { 
   
    edge1 = p7_spliceedge_Create();
    edge2 = p7_spliceedge_Create();

    edge1->overlap_amino_start = down_amino_start;
    edge1->overlap_amino_end   = down_amino_start + MAX_AMINO_OVERLAP/2;
    get_overlap_nuc_coords(edge1, upstream_domain, downstream_domain, splice_seq, revcomp);
    /* Add extra nucleotides for splice sites */
    edge1->upstream_nuc_end     = ESL_MIN(edge1->upstream_nuc_end + 2, splice_seq->n);
    edge1->downstream_nuc_start = ESL_MAX(edge1->downstream_nuc_start - 2, 1);

    if ((status = find_optimal_splice_site (edge1, pli, upstream_domain, downstream_domain, hmm, gcode, splice_seq)) != eslOK) goto ERROR;

    edge2->overlap_amino_start = up_amino_end - MAX_AMINO_OVERLAP/2;
    edge2->overlap_amino_end   = up_amino_end; 
    
    get_overlap_nuc_coords(edge2, upstream_domain, downstream_domain, splice_seq, revcomp);

    /* Add extra nucleotides for splice sites */
    edge2->upstream_nuc_end     = ESL_MIN(edge2->upstream_nuc_end + 2, splice_seq->n);
    edge2->downstream_nuc_start = ESL_MAX(edge2->downstream_nuc_start - 2, 1);

    if ((status = find_optimal_splice_site (edge2, pli, upstream_domain, downstream_domain, hmm, gcode, splice_seq)) != eslOK) goto ERROR;    
 
    if(edge1->splice_score == -eslINFINITY && edge2->splice_score == -eslINFINITY) {
      free(edge1);
      free(edge2);
      return NULL;
    }
    else if (edge1->splice_score > edge2->splice_score) {
      free(edge2);

      if(revcomp) {
        edge1->upstream_spliced_nuc_end     = splice_seq->n - edge1->upstream_spliced_nuc_end     + splice_seq->end;
        edge1->downstream_spliced_nuc_start = splice_seq->n - edge1->downstream_spliced_nuc_start + splice_seq->end;
      }
      else {
        edge1->upstream_spliced_nuc_end     = edge1->upstream_spliced_nuc_end     + splice_seq->start - 1;
        edge1->downstream_spliced_nuc_start = edge1->downstream_spliced_nuc_start + splice_seq->start - 1;
      }

      return edge1;
    }
    else {
      free(edge1);

      if(revcomp) {
        edge2->upstream_spliced_nuc_end     = splice_seq->n - edge2->upstream_spliced_nuc_end     + splice_seq->end;
        edge2->downstream_spliced_nuc_start = splice_seq->n - edge2->downstream_spliced_nuc_start + splice_seq->end;
      }
      else {
        edge2->upstream_spliced_nuc_end     = edge2->upstream_spliced_nuc_end     + splice_seq->start - 1;
        edge2->downstream_spliced_nuc_start = edge2->downstream_spliced_nuc_start + splice_seq->start - 1;
      }
      return edge2;
    }
     
  }
   
  
  ERROR:
    if(edge1 != NULL) free(edge1);
    if(edge2 != NULL) free(edge2);
    return NULL;

}

int
get_overlap_nuc_coords (SPLICE_EDGE *edge, const P7_DOMAIN *upstream, const P7_DOMAIN *downstream, const ESL_SQ *splice_seq, int revcomp)
{

  int       z1,z2;
  int       strand;
  int       extention_nuc_count;
  int       curr_hmm_pos;
  P7_TRACE *up_trace;
  P7_TRACE *down_trace;
  strand = (revcomp? -1 : 1);

  up_trace = upstream->tr;
  down_trace = downstream->tr;

 /***************UPSTREAM*******************/

 /* Get number of nucleotides that need to be added to the end of the
  * upstream edge nucleotide range if the end of the amino range
  * had to be extended to meet the mimimum overlap */
  if(edge->overlap_amino_end > upstream->jhmm) 
    extention_nuc_count = (edge->overlap_amino_end - upstream->jhmm) * strand * 3;
  else
    extention_nuc_count = 0;

  edge->upstream_nuc_end = upstream->jali + extention_nuc_count;
  
  for (z2 = up_trace->N-1 ; z2 >= 0; z2--) if (up_trace->st[z2] == p7T_M) break;
  edge->upstream_trace_end = z2;

  curr_hmm_pos = upstream->jhmm;
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
     
     if(curr_hmm_pos == edge->overlap_amino_end)
       break;
     z2--;
  }

  /* Set edge->upstream_nucl_start to the begining of the upstream
   * nulceotides that correpond to overlapping aminos */
  if(extention_nuc_count)
    edge->upstream_nuc_start = upstream->jali + strand;
  else
    edge->upstream_nuc_start =  edge->upstream_nuc_end + strand; 

  while(curr_hmm_pos >= edge->overlap_amino_start) {

    if      (up_trace->st[z2] == p7T_M) {
      edge->upstream_nuc_start -= up_trace->c[z2] * strand;
      if(up_trace->c[z2] != 3) edge->frameshift = TRUE;
      curr_hmm_pos--;
    }
    else if (up_trace->st[z2] == p7T_I)
      edge->upstream_nuc_start -= 3 * strand;
    else if (up_trace->st[z2] == p7T_D)
      curr_hmm_pos--;
    else
      ESL_EXCEPTION(eslFAIL, "splice boundry not found in trace");

    if(curr_hmm_pos == edge->overlap_amino_start)
      break;
    z2--;

  }
  edge->upstream_trace_start = z2;

  if(revcomp) {
    edge->upstream_nuc_start = splice_seq->start - edge->upstream_nuc_start + 1;
    edge->upstream_nuc_end   = splice_seq->start - edge->upstream_nuc_end   + 1;
  }
  else {
    edge->upstream_nuc_start = edge->upstream_nuc_start - splice_seq->start + 1;
    edge->upstream_nuc_end   = edge->upstream_nuc_end   - splice_seq->start + 1;
  }

  /***************DOWNSTREAM*******************/

 /* Get number of nucleotides that need to be added to the start of the
  * downstream edge nucleotide range if the end of the amino range
  * had to be extended to meet the mimimum overlap */
  if(edge->overlap_amino_start < downstream->ihmm)
    extention_nuc_count = (downstream->ihmm - edge->overlap_amino_start) * strand * 3;
  else
    extention_nuc_count = 0;

  edge->downstream_nuc_start = downstream->iali - extention_nuc_count;

  for (z1 = 0; z1 < down_trace->N; z1++) if (down_trace->st[z1] == p7T_M) break;
  edge->downstream_trace_start = z1;

  curr_hmm_pos = downstream->ihmm;
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
    
    if(curr_hmm_pos == edge->overlap_amino_start)
      break;
    z1++; 
  }

 /* Set edge->downstream_nucl_end to the end of the downstream
  * nulceotides that correpond to overlapping aminos */
  if(extention_nuc_count)
    edge->downstream_nuc_end = downstream->iali - strand;
  else 
    edge->downstream_nuc_end = edge->downstream_nuc_start - strand;
  
  while(curr_hmm_pos <= edge->overlap_amino_end) {
    if (down_trace->st[z1] == p7T_M) {
      edge->downstream_nuc_end += down_trace->c[z1] * strand;
      if(down_trace->c[z1] != 3) edge->frameshift = TRUE;
      curr_hmm_pos++;
    }
    else if (down_trace->st[z1] == p7T_I)
      edge->downstream_nuc_end += 3 * strand;
    else if (down_trace->st[z1] == p7T_D)
      curr_hmm_pos++;
    else
      ESL_EXCEPTION(eslFAIL, "splice boundry not found in trace");
    
    if(curr_hmm_pos > edge->overlap_amino_end)
      break;
    z1++;
  }
  edge->downstream_trace_end = z1;

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

  int            n,z;
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
  z = edge->upstream_trace_end;
  
  n = z - (upstream->tr->tfrom[0]+1);
  while(upstream->tr->k[z] >= edge->overlap_amino_start) {
    lost_sc += upstream->scores_per_pos[n];
    n--;
    z--;
  }
  
  z = edge->downstream_trace_start;
  n = 0;

  while(downstream->tr->k[z] && downstream->tr->k[z] <= edge->overlap_amino_end) {
    lost_sc += downstream->scores_per_pos[n];
    n++;
    z++;
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
  nuc_seq_len += edge->downstream_nuc_end - down_nuc_pos + 1;

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

