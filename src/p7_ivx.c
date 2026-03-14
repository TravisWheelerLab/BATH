/* P7_IVX implementation: a generic intermediate value matrix for a frameshift aware model
 *
 * Contents:
 *   1. The <P7_IVX> frameshift aware object
 */
#include "p7_config.h"
#include "hmmer.h"

/*****************************************************************
 * 1. The <P7_IVX> object.
 *****************************************************************/

/* Function:  p7_ivx_Create()
 * Synopsis:  Allocate a new <P7_IMX>.
 *
 * Purpose:   Allocate a reusable, resizeable <P7_IVX> for models up to
 *            size <allocM>
 *
 * Returns:   a pointer to the new <P7_IVX>.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_IVX *
p7_ivx_Create(int allocM, int allocC)
{
  int     status;
  P7_IVX *iv = NULL;

  ESL_ALLOC(iv, sizeof(P7_IVX));
  iv->ivx = NULL;

  ESL_ALLOC(iv->ivx, sizeof(float) * (allocM+1) * allocC);

  iv->allocM = allocM;
  iv->allocC = allocC;

  return iv;

 ERROR:
  if (iv != NULL) p7_ivx_Destroy(iv);
  return NULL;
}


/* Function:  p7_ivx_GrowTo()
 * Synopsis:  Assure that intermadiate values matrix is big enough.
 *
 * Returns:   <eslOK> on success, and <iv> may be reallocated upon
 *            return; any data that may have been in <iv> must be
 *            assumed to be invalidated.
 *
 * Throws:    <eslEMEM> on allocation failure, and any data that may
 *            have been in <gx> must be assumed to be invalidated.
 */
int
p7_ivx_GrowTo(P7_IVX *iv, int M, int C)
{
  int      status;
  void    *p;

  if((M+1) * C > (iv->allocM+1) * iv->allocC) {
    ESL_RALLOC(iv->ivx, p, sizeof(float) * (M+1) * C);
    iv->allocM = M;
    iv->allocC = C;
  }

  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_ivx_Destroy()
 * Synopsis:  Frees an intermediate values matrix.
 *
 * Purpose:   Frees a <P7_IVX>.
 *
 * Returns:   (void)
 */
void
p7_ivx_Destroy(P7_IVX *iv)
{
  if (iv == NULL) return;

  if (iv->ivx != NULL) free(iv->ivx);
  free(iv);
  return;
}


