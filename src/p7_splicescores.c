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
p7_splicepipeline_SignalScores(float *f)
{
  f[0] = log(0.9921);     /* GT-AG */
  f[1] = log(0.0073);     /* GC-AG */
  f[2] = log(0.0006);     /* AT-AC */
  return;
}

void
p7_splicepipeline_DonorGT(float *f)
{
  /* G = 2; T = 3; */
  /* x*4 + y */
  f[0]  = -eslINFINITY;
  f[1]  = -eslINFINITY;
  f[2]  = -eslINFINITY;
  f[3]  = -eslINFINITY;
  f[4]  = -eslINFINITY;
  f[5]  = -eslINFINITY;
  f[6]  = -eslINFINITY;
  f[7]  = -eslINFINITY;
  f[8]  = -eslINFINITY;
  f[9]  = -eslINFINITY;
  f[10] = -eslINFINITY;
  f[11] = 0.0f;
  f[12] = -eslINFINITY;
  f[13] = -eslINFINITY;
  f[14] = -eslINFINITY;
  f[15] = -eslINFINITY;
  return;
}

void
p7_splicepipeline_DonorGC(float *f)
{
  /* G = 2; C = 1; */
  /* x*4 + y */
  f[0]  = -eslINFINITY;
  f[1]  = -eslINFINITY;
  f[2]  = -eslINFINITY;
  f[3]  = -eslINFINITY;
  f[4]  = -eslINFINITY;
  f[5]  = -eslINFINITY;
  f[6]  = -eslINFINITY;
  f[7]  = -eslINFINITY;
  f[8]  = -eslINFINITY;
  f[9]  = 0.0f;
  f[10] = -eslINFINITY;
  f[11] = -eslINFINITY;
  f[12] = -eslINFINITY;
  f[13] = -eslINFINITY;
  f[14] = -eslINFINITY;
  f[15] = -eslINFINITY;
  return;
}

void
p7_splicepipeline_DonorAT(float *f)
{
  /* A = 0; T = 3; */
  /* x*4 + y */
  f[0]  = -eslINFINITY;
  f[1]  = -eslINFINITY;
  f[2]  = -eslINFINITY;
  f[3]  = 0.0f;
  f[4]  = -eslINFINITY;
  f[5]  = -eslINFINITY;
  f[6]  = -eslINFINITY;
  f[7]  = -eslINFINITY;
  f[8]  = -eslINFINITY;
  f[9]  = -eslINFINITY;
  f[10] = -eslINFINITY;
  f[11] = -eslINFINITY;
  f[12] = -eslINFINITY;
  f[13] = -eslINFINITY;
  f[14] = -eslINFINITY;
  f[15] = -eslINFINITY;
  return;
}

void
p7_splicepipeline_AcceptorAG(float *f)
{
  /* A = 0; G = 2; */
  /* x*4 + y */	
  f[0]  = -eslINFINITY;      
  f[1]  = -eslINFINITY;
  f[2]  = 0.0f;
  f[3]  = -eslINFINITY;
  f[4]  = -eslINFINITY;
  f[5]  = -eslINFINITY;
  f[6]  = -eslINFINITY;
  f[7]  = -eslINFINITY;
  f[8]  = -eslINFINITY;
  f[9]  = -eslINFINITY;
  f[10] = -eslINFINITY;
  f[11] = -eslINFINITY;
  f[12] = -eslINFINITY;
  f[13] = -eslINFINITY;
  f[14] = -eslINFINITY;
  f[15] = -eslINFINITY;
  return;
}

void
p7_splicepipeline_AcceptorAC(float *f)
{
  /* A = 0; C = 1 */
  /* x*4 + y */
  f[0]  = -eslINFINITY;
  f[1]  = 0.0f;
  f[2]  = -eslINFINITY;
  f[3]  = -eslINFINITY;
  f[4]  = -eslINFINITY;
  f[5]  = -eslINFINITY;
  f[6]  = -eslINFINITY;
  f[7]  = -eslINFINITY;
  f[8]  = -eslINFINITY;
  f[9]  = -eslINFINITY;
  f[10] = -eslINFINITY;
  f[11] = -eslINFINITY;
  f[12] = -eslINFINITY;
  f[13] = -eslINFINITY;
  f[14] = -eslINFINITY;
  f[15] = -eslINFINITY;
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
  ESL_ALLOC(splice_scores->score,     sizeof(float*) * SIGNAL_MEM_SIZE);

  for(i = 0; i < SIGNAL_MEM_SIZE; i++) 
    splice_scores->score[i] = splice_scores->score_mem + (i * M_hint);
  
  splice_scores->signal_scores = NULL; 
  ESL_ALLOC(splice_scores->signal_scores, sizeof(float) * p7S_SPLICE_SIGNALS);
  p7_splicepipeline_SignalScores(splice_scores->signal_scores);  

  splice_scores->donor_GT = splice_scores->donor_GC = splice_scores->donor_AT = NULL;
  splice_scores->acceptor_AG  = splice_scores->acceptor_AC = NULL;

  ESL_ALLOC(splice_scores->donor_GT, sizeof(float) * 16);
  ESL_ALLOC(splice_scores->donor_GC, sizeof(float) * 16);
  ESL_ALLOC(splice_scores->donor_AT, sizeof(float) * 16);

  ESL_ALLOC(splice_scores->acceptor_AG, sizeof(float) * 16);
  ESL_ALLOC(splice_scores->acceptor_AC, sizeof(float) * 16);

  p7_splicepipeline_DonorGT(splice_scores->donor_GT);
  p7_splicepipeline_DonorGC(splice_scores->donor_GC);
  p7_splicepipeline_DonorAT(splice_scores->donor_AT);

  p7_splicepipeline_AcceptorAG(splice_scores->acceptor_AG);
  p7_splicepipeline_AcceptorAC(splice_scores->acceptor_AC);

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

  for(i = 0; i < SIGNAL_MEM_SIZE; i++) 
    splice_scores->score[i] = splice_scores->score_mem + (i * M);

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

  if(splice_scores->donor_GT      != NULL) free(splice_scores->donor_GT);
  if(splice_scores->donor_GC      != NULL) free(splice_scores->donor_GC);
  if(splice_scores->donor_AT      != NULL) free(splice_scores->donor_AT);

  if(splice_scores->acceptor_AG   != NULL) free(splice_scores->acceptor_AG);
  if(splice_scores->acceptor_AC   != NULL) free(splice_scores->acceptor_AC);

  free(splice_scores);

  return;
  
}


