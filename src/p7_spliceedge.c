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
  edge->jump_edge = FALSE;

  edge->i_start = -1;
  edge->k_start = -1;

  edge->next_i_start = -1;
  edge->next_k_start = -1;
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

  int p, s, z;
  int overlap_start;
  int overlap_end;
  int overlap_len;
  int z1, z2;
  float min_lost_sc;
  float curr_lost_sc;
  float upstream_lost;
  float downstream_lost;
  float *spp;
  float *upstream_suffix_sum;
  float *downstream_prefix_sum;
  P7_TRACE    *tr;
  int status;

  upstream_suffix_sum   = NULL;
  downstream_prefix_sum = NULL;

  /* return if no there is no hmm overlap */
  if(downstream_dom->ihmm > upstream_dom->jhmm)  return eslOK; 

  overlap_start = ESL_MAX(upstream_dom->ihmm, downstream_dom->ihmm);
  overlap_end   = ESL_MIN(upstream_dom->jhmm, downstream_dom->jhmm);
  if(overlap_start == upstream_dom->ihmm)   overlap_start++;
  if(overlap_end   == downstream_dom->jhmm) overlap_end--;
   
  overlap_len   = overlap_end - overlap_start + 1;

  if(overlap_len < 1) {
    edge->splice_score = -eslINFINITY;
    return eslOK;
  }

  ESL_ALLOC(upstream_suffix_sum,   sizeof(float) * overlap_len);
  ESL_ALLOC(downstream_prefix_sum, sizeof(float) * overlap_len);
  
  esl_vec_FSet(upstream_suffix_sum,   overlap_len, 0.0);
  esl_vec_FSet(downstream_prefix_sum, overlap_len, 0.0);

  /* Fill the upstream array */
  tr  = upstream_dom->tr;
  spp = upstream_dom->scores_per_pos;
  
  z1 = tr->tfrom[0];
  z2 = tr->tto[0];
  while(tr->st[z2] != p7T_M) z2--;

  /* Move to the end of the overlap in the trace and the scores_per_pos array */
  p = z2 - z1 - 2;
  z = z2;
  while(tr->k[z] > overlap_end) { z--; p--;}

  s = overlap_len - 1;
  while(tr->k[z]  >= overlap_start) {
    upstream_suffix_sum[s] += spp[p];
    z--;
    p--;
    if(tr->k[z] < tr->k[z+1]) s--;
  }
 
  for(s = overlap_len - 2; s >= 0; s--)
    upstream_suffix_sum[s] += upstream_suffix_sum[s+1];

  upstream_lost = 0.;
  /*Calculate any lost score from the end of the upstream hit */
  if(upstream_dom->jhmm > overlap_end) {
    z = z1+1;
    p = 0;
    while(tr->k[z] <= overlap_end) { z++; p++; }
    while(z < z2) {
      upstream_lost += spp[p];
      z++;
      p++;
    }
  }
 
  /* Fill the downstream array */  
  tr  = downstream_dom->tr;
  spp = downstream_dom->scores_per_pos;

  z1 = tr->tfrom[0];
  z2 = tr->tto[0];
 
 /* Move to the start of the overlap in the trace and the scores_per_pos array */
  p = 0;
  z = z1+1;
  while(tr->k[z] < overlap_start) { z++; p++;}
  
  s = 0;
  while(tr->k[z] <= overlap_end && z < z2) {
    downstream_prefix_sum[s] += spp[p];
    z++;
    p++;
    if(tr->k[z] > tr->k[z-1]) s++;
  }

  for(s = 1; s < overlap_len; s++)
    downstream_prefix_sum[s] +=  downstream_prefix_sum[s-1];
  
  downstream_lost = 0.;
  /*Calculate any lost score from the start of the downstream hit */
  if(downstream_dom->ihmm < overlap_start) {
    z = z2-1;
    p = z2 - z1 - 2;
    while(tr->k[z] >= overlap_start) { z--; p--; }
    while(z > z1) {
      downstream_lost += spp[p];
      z--;
      p--;
    }
  } 


  /* Find the minimiom score loss to eliminate the overlap*/ 
  min_lost_sc = upstream_suffix_sum[0];
  for(s = 1; s < overlap_len; s++) {
    curr_lost_sc = upstream_suffix_sum[s] + downstream_prefix_sum[s-1];
    min_lost_sc = ESL_MIN(min_lost_sc, curr_lost_sc);
  }
  min_lost_sc = ESL_MIN(min_lost_sc, downstream_prefix_sum[overlap_len-1]);

  edge->splice_score -= (min_lost_sc + upstream_lost + downstream_lost);
    
  if(upstream_suffix_sum   != NULL) free(upstream_suffix_sum);
  if(downstream_prefix_sum != NULL) free(downstream_prefix_sum);
    
  return eslOK;

  ERROR:
    if(upstream_suffix_sum   != NULL) free(upstream_suffix_sum);
    if(downstream_prefix_sum != NULL) free(downstream_prefix_sum);
    return status;

}



