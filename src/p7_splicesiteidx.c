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
 * Purpose:   Allocates a new <SPLICE_SITE_IDX>  
 *
 * Returns:   a pointer to the new <SPLICE_SITE_IDX> structure 
 *            on success.
 *
 * Throws:    <NULL> on allocation error.
 */
SPLICE_SITE_IDX*
p7_splicesiteidx_Create()
{
  int i;
  int status;

  SPLICE_SITE_IDX *signal_sites = NULL;
  ESL_ALLOC(signal_sites, sizeof(SPLICE_SITE_IDX));
 
  ESL_ALLOC(signal_sites->index_mem, sizeof(int) * p7S_SPLICE_SIGNALS * SIGNAL_MEM_SIZE);
  ESL_ALLOC(signal_sites->index,     sizeof(int*) * p7S_SPLICE_SIGNALS);
  
  for(i = 0; i < p7S_SPLICE_SIGNALS; i++) {
    signal_sites->index[i] = signal_sites->index_mem + (p7S_SPLICE_SIGNALS * SIGNAL_MEM_SIZE); 
  }

  return signal_sites;

  ERROR:
    p7_splicesiteidx_Destroy(signal_sites);
    return NULL;

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

  if(signal_sites->index     != NULL) free(signal_sites->index);
  if(signal_sites->index_mem != NULL) free(signal_sites->index_mem);

  signal_sites = NULL;

  return;
}


