/* SPLICE_BOUNDS: list of hmm and sequence coordinate boundries to constrain paths
 *
 * Contents:
 *    1. The SPLICE_BOUNDS object.
 *
 */


#include "p7_config.h"

#include <string.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_splice.h"

/*****************************************************************
 * 1. The SPLICE_BOUNDS structure.
 *****************************************************************/

/* Function:  p7_splicebounds_Create()
 *
 * Purpose:   Allocates a splice bounds object
 *
 * Returns:   a pointer to the new <SPLICE_BOUNDS> structure
 *            on success.
 *
 * Throws:    <NULL> on allocation error.
 */
SPLICE_BOUNDS* 
p7_splicebounds_Create(int allocN)
{
  SPLICE_BOUNDS *bounds;
  int            status;
 
  bounds = NULL;
  ESL_ALLOC(bounds, sizeof(SPLICE_BOUNDS));

  bounds->bound_hmm_mins = NULL;
  bounds->bound_hmm_maxs = NULL;
  bounds->bound_seq_mins = NULL;
  bounds->bound_seq_maxs = NULL;

  ESL_ALLOC(bounds->bound_hmm_mins, sizeof(int)     * allocN);
  ESL_ALLOC(bounds->bound_hmm_maxs, sizeof(int)     * allocN);
  ESL_ALLOC(bounds->bound_seq_mins, sizeof(int64_t) * allocN);
  ESL_ALLOC(bounds->bound_seq_maxs, sizeof(int64_t) * allocN);

  bounds->N = 0;
  bounds->allocN = allocN;

  return bounds;

  ERROR:
    p7_splicebounds_Destroy(bounds);
    return NULL;

}

/* Function:  p7_splicebounds_GrowTo()
 * Synopsis:  Grow a splice bounds object.
 *
 * Purpose:   Allocates a larger <SPLICE_BOUNDS>
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_splicebounds_GorwTo(SPLICE_BOUNDS *bounds, int allocN)
{

  int status;

  if(allocN <= bounds->allocN) return eslOK;

  ESL_REALLOC(bounds->bound_hmm_mins, sizeof(int)     * allocN);
  ESL_REALLOC(bounds->bound_hmm_maxs, sizeof(int)     * allocN);
  ESL_REALLOC(bounds->bound_seq_mins, sizeof(int64_t) * allocN);
  ESL_REALLOC(bounds->bound_seq_maxs, sizeof(int64_t) * allocN); 

  bounds->allocN = allocN;

  return eslOK;

  ERROR:
    p7_splicebounds_Destroy(bounds);
    return status;

}


/* Function: p7_splicebounds_Destroy()
 *
 * Purpose:  Frees a <SPLICE_BOUNDS>
 */
void 
p7_splicebounds_Destroy(SPLICE_BOUNDS *bounds)
{

  if(bounds == NULL) return;

  if(bounds->bound_hmm_mins != NULL) free(bounds->bound_hmm_mins);
  if(bounds->bound_hmm_maxs != NULL) free(bounds->bound_hmm_maxs);
  if(bounds->bound_seq_mins != NULL) free(bounds->bound_seq_mins);
  if(bounds->bound_seq_maxs != NULL) free(bounds->bound_seq_maxs);
 
  free(bounds);

 return;
  
}

/* Function:  p7_splicebounds_ADD()
 * Synopsis:  add new bound to <SPLICE_BOUNDS>.
 *
 * Returns:   <eslOK> on success.
 *
 */
int 
p7_splicebounds_Add(SPLICE_BOUNDS *bounds, int64_t seq_min, int64_t seq_max, int hmm_min, int hmm_max)
{

  if(bounds->N == bounds->allocN) 
    p7_splicebounds_GorwTo(bounds, bounds->allocN *2);

  bounds->bound_hmm_mins[bounds->N] = hmm_min;
  bounds->bound_hmm_maxs[bounds->N] = hmm_max;
  bounds->bound_seq_mins[bounds->N] = seq_min;
  bounds->bound_seq_maxs[bounds->N] = seq_max;

  bounds->N++;
 
  return eslOK;
}
