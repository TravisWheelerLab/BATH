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
  edge->jump_edge = FALSE;

  edge->splice_score = -eslINFINITY;
  edge->signal_score = -eslINFINITY;

  return edge;

  ERROR:
    if (edge != NULL) free(edge);
    return NULL;
}


int
p7_spliceedge_AliScoreEdge(SPLICE_EDGE *edge, const P7_PROFILE *gm, const P7_DOMAIN *upstream_dom, const P7_DOMAIN *downstream_dom)
{

  int k, p, s, z;
  int us, ds;
  int k_start, k_end;
  int k_trans;
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
  
  /* If the hits overlap on the hmm coords use the the per postition 
   * scores to find the optimal place to switch from one hit to the 
   * next and subtrct the lost score from the edge score */
  if(downstream_dom->ihmm <= upstream_dom->jhmm) {

    /* Create an array of running sums of the per position sore suffixes for the upstream hit */
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
    /* Sum all postions (M & I) from each k */
    while (k <= upstream_dom->jhmm) {
      if(tr->k[z] == k) {
        upstream_suffix_sum[s] += spp[p];
        z++; p++;
      }
      else {
        k++; s++;
      }
    }
    /* Sum suffixes */
    for(s = up_N-2; s >= 0; s--)
      upstream_suffix_sum[s] += upstream_suffix_sum[s+1];

    /* Create an array of running sums of the per position sore prefixes for the downstream hit */    
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
     /* Sum all postions (M & I) from each k */
    while (k <= downstream_dom->jhmm) {
       if(tr->k[z] == k) {
        downstream_prefix_sum[s] += spp[p];
        z++; p++;
      }
      else {
        k++; s++;
      }
    }
    /* Sum prefixes */
    for(s = 1; s < down_N; s++)
      downstream_prefix_sum[s] += downstream_prefix_sum[s-1];

    /*Step through the overlap postions to find the best (minimiun lost score) transition point */
    k_start = ESL_MAX(downstream_dom->ihmm, upstream_dom->ihmm+1);
    k_end   = ESL_MIN(upstream_dom->jhmm,   downstream_dom->jhmm-1);
    
    us = k_start - upstream_dom->ihmm;
    ds = k_start - downstream_dom->ihmm - 1;
    min_lost_sc = upstream_suffix_sum[us];  // lost_sc if all k belong to downstream
    if(ds >= 0)
      min_lost_sc += downstream_prefix_sum[ds];  // include any lost_sc from downstream overextention
    k_trans = k_start;

    us++;
    ds++;
    for(k = k_start; k < k_end; k++) {
      /* lost_sc if k_start to k belong to upstream and k+1 to k_end belong to downstream */
      lost_sc = upstream_suffix_sum[us] +  downstream_prefix_sum[ds];
      if(lost_sc < min_lost_sc) {
        min_lost_sc = lost_sc;
        k_trans = k;
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
      k_trans = k_end;
    }

    /* edge score equals lost score plus M->M trastion for the postion where we transtion from up to down*/
    edge->splice_score -= min_lost_sc + gm->tsc[k_trans * p7P_NTRANS + p7P_MM];
    
    if(upstream_suffix_sum   != NULL) free(upstream_suffix_sum);
    if(downstream_prefix_sum != NULL) free(downstream_prefix_sum);

  }
  else {

    /* Edge score is sum of M->M trasntions for all missing positions */
    for(k = upstream_dom->jhmm; k < downstream_dom->ihmm; k++)
      edge->splice_score +=  gm->tsc[k * p7P_NTRANS + p7P_MM];
  }

  return eslOK;

  ERROR:
    if(upstream_suffix_sum   != NULL) free(upstream_suffix_sum);
    if(downstream_prefix_sum != NULL) free(downstream_prefix_sum);
    return status;

}



