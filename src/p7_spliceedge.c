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

  return edge;

  ERROR:
    if (edge != NULL) free(edge);
    return NULL;
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
    if(us == up_N) us--;
    lost_sc = downstream_prefix_sum[ds]; // lost_sc if all k belong to upstream
    if(k_end < upstream_dom->jhmm) {
 
      lost_sc += upstream_suffix_sum[us]; // include any lost_sc from upstream overextention
    }
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



