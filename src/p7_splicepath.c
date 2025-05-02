/* SPLICE_PATH: a set of spliced hits 
 *
 * Contents:
 *    1. The SPLICE_PATH object.
 *    2. Debugging tools.
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
 * 1. The SPLICE_PATH structure.
 *****************************************************************/

/* Function:  p7_splicepath_Create()
 *
 * Purpose:   Allocates a splice path with room for 2x <path_len> 
 *            hits.
 *
 * Returns:   a pointer to the new <SPLICE_PATH> structure
 *            on success.
 *
 * Throws:    <NULL> on allocation error.
 */
SPLICE_PATH*
p7_splicepath_Create(int path_len)
{

  SPLICE_PATH *path;
  int          status;

  path = NULL;
  ESL_ALLOC(path, sizeof(SPLICE_PATH));

  path->alloc_len = path_len*2;
  path->path_len  = path_len;

  path->node_id = NULL;
  ESL_ALLOC(path->node_id,                        sizeof(int)*(path_len*2));
  path->split = NULL;
  ESL_ALLOC(path->split,                          sizeof(int)*(path_len*2));

  path->upstream_spliced_amino_end     = NULL;
  path->downstream_spliced_amino_start = NULL;
  ESL_ALLOC(path->upstream_spliced_amino_end,     sizeof(int)*(path_len*2+1));
  ESL_ALLOC(path->downstream_spliced_amino_start, sizeof(int)*(path_len*2+1));

  path->upstream_spliced_nuc_end     = NULL;
  path->downstream_spliced_nuc_start = NULL;
  ESL_ALLOC(path->upstream_spliced_nuc_end,       sizeof(int)*(path_len*2+1));
  ESL_ALLOC(path->downstream_spliced_nuc_start,   sizeof(int)*(path_len*2+1));

  path->hit_scores = NULL;
  ESL_ALLOC(path->hit_scores, sizeof(float)*path_len*2);

  path->edge_scores = NULL;
  ESL_ALLOC(path->edge_scores, sizeof(float)*path_len*2);

  path->signal_scores = NULL;
  ESL_ALLOC(path->signal_scores, sizeof(float)*path_len*2);

  path->hits = NULL;
  ESL_ALLOC(path->hits,          sizeof(P7_HIT*)*path_len*2);

  return path;

  ERROR:
    p7_splicepath_Destroy(path);
    return NULL;
}

/* Function:  p7_splicepath_Grow()
 * Synopsis:  Reallocates a larger splice path, if needed.
 *
 * Purpose:   If <SPLICE_PATH> cannot hold another hit,
 *            doubles the internal allocation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *
 */
int
p7_splicepath_Grow(SPLICE_PATH *path)
{

  int status;

  if(path->path_len < path->alloc_len) return eslOK;

  path->alloc_len *= 2;

  ESL_REALLOC(path->node_id,                        sizeof(int)     * path->alloc_len);
  ESL_REALLOC(path->split,                          sizeof(int)     * path->alloc_len);
  ESL_REALLOC(path->upstream_spliced_amino_end,     sizeof(int)     * (path->alloc_len+1));
  ESL_REALLOC(path->downstream_spliced_amino_start, sizeof(int)     * (path->alloc_len+1));
  ESL_REALLOC(path->upstream_spliced_nuc_end,       sizeof(int)     * (path->alloc_len+1));
  ESL_REALLOC(path->downstream_spliced_nuc_start,   sizeof(int)     * (path->alloc_len+1));
  ESL_REALLOC(path->hit_scores,                     sizeof(float)   * path->alloc_len);
  ESL_REALLOC(path->edge_scores,                    sizeof(float)   * path->alloc_len);
  ESL_REALLOC(path->signal_scores,                  sizeof(float)   * path->alloc_len);
  ESL_REALLOC(path->hits,                           sizeof(P7_HIT*) * path->alloc_len);

  return eslOK;

  ERROR:
    p7_splicepath_Destroy(path);
    return status;

}

/* Function: p7_splicepath_Destroy()
 *
 * Purpose:  Frees a <SPLICE_PATH>
 */
void
p7_splicepath_Destroy(SPLICE_PATH *path)
{

   int i;

   if(path == NULL) return;

   /* Destroy split hits */
   for(i = 0; i < path->path_len; i++) {
     if(path->split[i]) {

       p7_alidisplay_Destroy(path->hits[i]->dcl->ad);
       p7_trace_fs_Destroy(path->hits[i]->dcl->tr);
       free(path->hits[i]->dcl->scores_per_pos);
       p7_hit_Destroy(path->hits[i]);
     }
   }

   if(path->node_id != NULL) free(path->node_id);
   if(path->split != NULL)   free(path->split);

   if(path->upstream_spliced_amino_end     != NULL)
     free(path->upstream_spliced_amino_end);
   if(path->downstream_spliced_amino_start != NULL)
     free(path->downstream_spliced_amino_start);

   if(path->upstream_spliced_nuc_end       != NULL)
     free(path->upstream_spliced_nuc_end);
   if(path->downstream_spliced_nuc_start   != NULL)
     free(path->downstream_spliced_nuc_start);

   if(path->hit_scores    != NULL) free(path->hit_scores);
   if(path->edge_scores   != NULL) free(path->edge_scores);
   if(path->signal_scores != NULL) free(path->signal_scores);
   if(path->hits          != NULL) free(path->hits);
   if(path                != NULL) free(path);

   return;
}


/*****************************************************************
 * 2. Debugging tools.
 *****************************************************************/


/* Function:  p7_splicepath_Dump()
 *
 * Purpose: Dumps splice coords and score data for 
 *          each hit in <path> 
 *
 */
void
p7_splicepath_Dump(FILE *fp, SPLICE_PATH *path)
{

  int i;

  fprintf(fp, "  Path Length  %d\n", path->path_len);
  fprintf(fp, "  %4s %4s %9s %9s %10s %10s %9s %10s \n", "Step", "Node", "hmm_start", "hmm_end", "seq_start", "seq_end", "hit_score", "edge_score");
  for(i = 0; i < path->path_len; i++) {
    fprintf(fp, "  %4d %4d %9d %9d %10d %10d %9.2f %10.2f\n", i+1, path->node_id[i]+1,
      path->downstream_spliced_amino_start[i], path->upstream_spliced_amino_end[i+1], path->downstream_spliced_nuc_start[i], path->upstream_spliced_nuc_end[i+1], path->hit_scores[i], path->edge_scores[i]);
  }

  fprintf(fp, "\n");

  return;
}
