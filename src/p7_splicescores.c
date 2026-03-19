/* SPLICE_SCORES: 
 * hardcoded splice signal scores,
 * hardcoded donor/acceptor validation arrays,
 * space for P state score storage
 *
 * Contents:
 *    1. The SPLICE_SCORES object.
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

/* Splice singal probabilities taken from
 * "Comprehensive splice-site analysis using comparative genomics",
 * Nihar Sheth et al., 2006
 */
void
p7_splicescores_SignalScores(float *f)
{
  f[0] = log(0.9921);     /* GT-AG */
  f[1] = log(0.0073);     /* GC-AG */
  f[2] = log(0.0006);     /* AT-AC */
  return;
}



/*****************************************************************
 * 1. The SPLICE_SCORES structure.
 *****************************************************************/


/* Function:  p7_splicescores_Create()
 *
 * Purpose:   Allocates a <SPLICE_SCORES>, initializes
 *            hardcoded data, and allocates space P-sate
 *            score staorage for an <M_hint> length modle 
 *
 * Returns:   a pointer to the new <SPLICE_SCORES> structure
 *            on success.
 *
 * Throws:    <NULL> on allocation error.
 */
SPLICE_SCORES* 
p7_splicescores_Create(int M_hint)
{
  int i;
  SPLICE_SCORES *splice_scores;
  int status;
 
  splice_scores = NULL;
  ESL_ALLOC(splice_scores, sizeof(SPLICE_SCORES));

  splice_scores->allocM = M_hint;

  splice_scores->score_mem = NULL;
  splice_scores->score = NULL;
  ESL_ALLOC(splice_scores->score_mem, sizeof(float)  * M_hint * SIGNAL_MEM_SIZE);
  ESL_ALLOC(splice_scores->score,     sizeof(float*) * M_hint);

  for(i = 0; i < M_hint; i++)
    splice_scores->score[i] = splice_scores->score_mem + (i * SIGNAL_MEM_SIZE);
  
  splice_scores->signal_scores = NULL; 
  ESL_ALLOC(splice_scores->signal_scores, sizeof(float) * p7S_SPLICE_SIGNALS);
  p7_splicescores_SignalScores(splice_scores->signal_scores);  

  return splice_scores;

  ERROR:
    p7_splicescores_Destroy(splice_scores);
    return NULL;

}

/* Function:  p7_splicescores_GrowTo() 
 * Synopsis:  Grow P state score stroage
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_splicescores_GrowTo(SPLICE_SCORES *splice_scores, int M)
{

  int i;
  int status;

  if(M <= splice_scores->allocM) return eslOK;

  ESL_REALLOC(splice_scores->score_mem, sizeof(float) * M * SIGNAL_MEM_SIZE);
  ESL_REALLOC(splice_scores->score,     sizeof(float*) * M);

  for(i = 0; i < M; i++)
    splice_scores->score[i] = splice_scores->score_mem + (i * SIGNAL_MEM_SIZE);

  splice_scores->allocM = M; 

  return eslOK;

  ERROR:
  p7_splicescores_Destroy(splice_scores);
  return status;
}


/* Function: p7_splicescores_Destroy()
 *
 * Purpose:  Frees a <SPLICE_SCORES>
 */
void 
p7_splicescores_Destroy(SPLICE_SCORES *splice_scores)
{

  if(splice_scores == NULL) return;

  if(splice_scores->score         != NULL) free(splice_scores->score);
  if(splice_scores->score_mem     != NULL) free(splice_scores->score_mem);

  if(splice_scores->signal_scores != NULL) free(splice_scores->signal_scores);

  free(splice_scores);

  return;
  
}


