/* SPLICE_SITE_IDX: list of indices in a sequence where donor splice signals can be found
 *
 * Contents:
 *    1.  allocation and destruction of the SPLICE_SITE_IDX
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
 * 1. The SPLICE_SITE_IDX structure.
 *****************************************************************/

/* Function:  p7_splicesiteidx_Create()
 * Synopsis:  Allocates a splice site idx.
 *
 * Purpose:   Allocates a new <SPLICE_SITE_IDX> with memory for 
 *            the expected occurance of 2-mers in a nucleotide 
 *            sequence of length <L>
 *
 * Returns:   a pointer to the new <SPLICE_SITE_IDX> structure 
 *            on success.
 *
 * Throws:    <NULL> on allocation error.
 */
SPLICE_SITE_IDX*
p7_splicesiteidx_Create(int L);
{
  int i;
  int status;


  SPLICE_SITE_IDX *signal_sites = NULL;
  ESL_ALLOC(signal_sites, sizeof(SPLICE_SITE_IDX));
 
  ESL_ALLOC(signal_sites->index,  sizeof(int*) * p7S_SPLICE_SIGNALS);
  ESL_ALLOC(signal_sites->N,      sizeof(int)  * p7S_SPLICE_SIGNALS);
  ESL_ALLOC(signal_sites->allocN, sizeof(int)  * p7S_SPLICE_SIGNALS);
  
  for(i = 0; i < p7S_SPLICE_SIGNALS; i++) {
    signal_sites->N[i]       = 0;
    signal_sites->allocN[i]  = (L-1) / 16;
    ESL_ALLOC(signal_sites->index[i], sizeof(int) * signal_sites->allocN[i]);
  }

  return signal_sites;

ERROR:
    p7_splicesiteidx_Destroy(signal_sites);
    return NULL;

}


/* Function:  p7_splicesiteidx_Grow() 
 * Synopsis:  Reallocates a larger SPLICE_SITE_IDX, if needed.
 *
 * Purpose:   If <signal_sites> cannot hold another index for signal 
 *            <signal> increase the allocation by the expected 
 *            occurance of 2-mers in a nucleotide seq of length <L>
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 *   
 */
int
p7_splicesiteidx_Grow(SPLICE_SITE_IDX *signal_sites, int L, int signal)
{
  int i,n;
  int new_alloc;
  int status;

  if(signal_sites->N[signal] == signal_sites->allocN) {

    signal_sites->allocN[signal] += (L-1) / 16;
    ESL_REALLOC(signal_sites->index[signal], sizeof(int) * signal_sites->allocN[signal]); 

  }

  return eslOK;
  
  ERROR:
    p7_splicesiteidx_Destroy(signal_sites);  
    return status;

}


/* Function:  p7_splicesiteidx_Destroy()
 *
 * Purpose:  Frees a <SPLICE_SITE_IDX>
 */
void 
p7_splicesiteidx_Destroy(SPLICE_SITE_IDX *signal_sites)
{

  int i;

  if(signal_sites == NULL) return;

  free(signal_sites->N);
  free(signal_sites->allocN);
    
  for(i = 0; i < p7S_SPLICE_SIGNALS; i++) {
    free(signal_sites->index[i]);
  }
  free(signal_sites->index);
  free(signal_sites); 

  signal_sites = NULL;

  return;
}


