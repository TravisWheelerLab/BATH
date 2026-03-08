/* Routines for the P7_FS_OPROFILE structure:  
 * a search profile in an optimized implementation of framwshift aware search.
 * 
 * Contents:
 *   1. The P7_FS_OPROFILE object: allocation, initialization, destruction.
 *   2. Conversion from generic P7_FS_PROFILE to optimized P7_FS_OPROFILE
 *   3. Conversion from optimized P7_FS_OPROFILE to compact score arrays
 *   4. Debugging and development utilities.
 *   5. Benchmark driver.
 *   6. Unit tests.
 *   7. Test driver.
 *   8. Example.
 */
#include "p7_config.h"

#include <stdio.h>
#include <string.h>
#include <math.h>		/* roundf() */

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_random.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_sse.h"

/*********************************************************************
 * 1. The P7_FS_OPROFILE structure: a frameshift aware score profile.
 *********************************************************************/

/* Function:  p7_fs_oprofile_Create()
 * Synopsis:  Allocate an optimized frameshift profile structure.
 *
 * Purpose:   Allocate a P7_FS_OPROFILE for profiles of up to <allocM>
 *            nodes, digital alphabet <abc>, and <codon_lengths> codon
 *            length variants (1, 3, or 5).
 *
 *            The emission table <rfv> has one row per codon/quasicodon
 *            type plus one row per alphabet symbol:
 *              codon_lengths == 1: p7P_MAXCODONS1 + abc->Kp rows
 *              codon_lengths == 3: p7P_MAXCODONS3 + abc->Kp rows
 *              codon_lengths == 5: p7P_MAXCODONS5 + abc->Kp rows
 *
 * Throws:    <NULL> on allocation error or invalid <codon_lengths>.
 */
P7_FS_OPROFILE *
p7_fs_oprofile_Create(int allocM, const ESL_ALPHABET *abc, int codon_lengths)
{
  int              status;
  P7_FS_OPROFILE  *om_fs = NULL;
  int              nqf;          /* # of float vectors per stripe      */
  int              ncodon_rows;  /* # of rows in rfv emission table    */
  int              c;

  if      (codon_lengths == 1) ncodon_rows = p7P_MAXCODONS1 + abc->Kp;
  else if (codon_lengths == 3) ncodon_rows = p7P_MAXCODONS3 + abc->Kp;
  else if (codon_lengths == 5) ncodon_rows = p7P_MAXCODONS5 + abc->Kp;
  else    ESL_XEXCEPTION(eslEINVAL, "codon_lengths must be 1, 3 or 5");

  nqf = p7O_NQF(allocM);

  /* level 0 */
  ESL_ALLOC(om_fs, sizeof(P7_FS_OPROFILE));
  om_fs->rfv_mem = NULL;
  om_fs->tfv_mem = NULL;
  om_fs->rfv     = NULL;
  om_fs->tfv     = NULL;
  om_fs->clone   = 0;

  /* level 1: raw memory (+15 bytes for manual 16-byte alignment) */
  ESL_ALLOC(om_fs->rfv_mem, sizeof(__m128) * nqf * ncodon_rows + 15);
  ESL_ALLOC(om_fs->tfv_mem, sizeof(__m128) * nqf * p7O_NTRANS  + 15);

  /* pointer array: one __m128 * per emission row */
  ESL_ALLOC(om_fs->rfv, sizeof(__m128 *) * ncodon_rows);

  /* align vector memory on 16-byte boundaries */
  om_fs->rfv[0] = (__m128 *) (((unsigned long int) om_fs->rfv_mem + 15) & (~0xf));
  om_fs->tfv    = (__m128 *) (((unsigned long int) om_fs->tfv_mem + 15) & (~0xf));

  /* set remaining row pointers */
  for (c = 1; c < ncodon_rows; c++)
    om_fs->rfv[c] = om_fs->rfv[0] + (c * nqf);

  om_fs->allocQ4 = nqf;

  /* Frameshift-specific initializations */
  om_fs->codon_lengths = codon_lengths;
  om_fs->fs            = 0.0f;

  for (c = 0; c < p7_NOFFSETS; c++) om_fs->offs[c]    = -1;
  for (c = 0; c < p7_NEVPARAM; c++) om_fs->evparam[c] = p7_EVPARAM_UNSET;
  for (c = 0; c < p7_NCUTOFFS; c++) om_fs->cutoff[c]  = p7_CUTOFF_UNSET;
  for (c = 0; c < p7_MAXABET;  c++) om_fs->compo[c]   = p7_COMPO_UNSET;

  om_fs->name = NULL;
  om_fs->acc  = NULL;
  om_fs->desc = NULL;

  /* Always allocate RF, MM, CS, consensus annotation; rely on leading \0 to
   * signal unused. Zeroed to suppress valgrind warnings on fwrite() in io. */
  ESL_ALLOC(om_fs->rf,        sizeof(char) * (allocM+2));
  ESL_ALLOC(om_fs->mm,        sizeof(char) * (allocM+2));
  ESL_ALLOC(om_fs->cs,        sizeof(char) * (allocM+2));
  ESL_ALLOC(om_fs->consensus, sizeof(char) * (allocM+2));
  memset(om_fs->rf,        '\0', sizeof(char) * (allocM+2));
  memset(om_fs->mm,        '\0', sizeof(char) * (allocM+2));
  memset(om_fs->cs,        '\0', sizeof(char) * (allocM+2));
  memset(om_fs->consensus, '\0', sizeof(char) * (allocM+2));

  om_fs->abc        = abc;
  om_fs->L          = 0;
  om_fs->M          = 0;
  om_fs->max_length = -1;
  om_fs->allocM     = allocM;
  om_fs->mode       = p7_NO_MODE;
  om_fs->nj         = 0.0f;
  return om_fs;

 ERROR:
  p7_fs_oprofile_Destroy(om_fs);
  return NULL;
}

/* Function:  p7_fs_oprofile_IsLocal()
 * Synopsis:  Returns TRUE if profile is in local alignment mode.
 */
int
p7_fs_oprofile_IsLocal(const P7_FS_OPROFILE *om_fs)
{
  if (om_fs->mode == p7_LOCAL || om_fs->mode == p7_UNILOCAL) return TRUE;
  return FALSE;
}


/* Function:  p7_fs_oprofile_Destroy()
 * Synopsis:  Frees an optimized frameshift profile structure.
 */
void
p7_fs_oprofile_Destroy(P7_FS_OPROFILE *om_fs)
{
  if (om_fs == NULL) return;

  if (om_fs->clone == 0)
    {
      if (om_fs->rfv_mem   != NULL) free(om_fs->rfv_mem);
      if (om_fs->tfv_mem   != NULL) free(om_fs->tfv_mem);
      if (om_fs->rfv       != NULL) free(om_fs->rfv);
      if (om_fs->name      != NULL) free(om_fs->name);
      if (om_fs->acc       != NULL) free(om_fs->acc);
      if (om_fs->desc      != NULL) free(om_fs->desc);
      if (om_fs->rf        != NULL) free(om_fs->rf);
      if (om_fs->mm        != NULL) free(om_fs->mm);
      if (om_fs->cs        != NULL) free(om_fs->cs);
      if (om_fs->consensus != NULL) free(om_fs->consensus);
    }

  free(om_fs);
}

/* Function:  p7_fs_oprofile_Clone()
 * Synopsis:  Allocate a cloned copy of an optimized frameshift profile.
 *
 * Purpose:   Shallow copy of <om_fs> for use in multiple threads. The
 *            clone shares all allocated memory with the original; no
 *            vector data is reallocated. <p7_fs_oprofile_Destroy()> on
 *            the clone frees only the shell struct, not the shared memory.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_FS_OPROFILE *
p7_fs_oprofile_Clone(const P7_FS_OPROFILE *om_fs)
{
  int              status;
  P7_FS_OPROFILE  *om2 = NULL;

  ESL_ALLOC(om2, sizeof(P7_FS_OPROFILE));
  memcpy(om2, om_fs, sizeof(P7_FS_OPROFILE));

  om2->clone = 1;

  return om2;

 ERROR:
  p7_fs_oprofile_Destroy(om2);
  return NULL;
}

/*----------------- end, P7_FS_OPROFILE structure ------------------*/



/*********************************************************************
 * 2. Conversion from generic P7_FS_PROFILE to optimized P7_FS_OPROFILE
 *********************************************************************/

/*------------ end, conversions to P7_FS_OPROFILE ------------------*/

/***********************************************************************
*   3. Conversion from optimized P7_FS_OPROFILE to compact score arrays
 ***********************************************************************/


/*------------ end, conversions from P7_FS_OPROFILE ------------------*/


/*********************************************************************
 * 4. Debugging and development utilities.
 *********************************************************************/

/*------------ end, P7_FS_OPROFILE debugging tools  ----------------*/



/*****************************************************************
 * 5. Benchmark driver.
 *****************************************************************/

/*---------------- end, benchmark driver ------------------------*/




  
/*****************************************************************
 * 6. Unit tests
 *****************************************************************/

/*------------------- end, unit tests ---------------------------*/




/*****************************************************************
 * 7. Test driver
 *****************************************************************/

/*------------------- end, test driver --------------------------*/


/*****************************************************************
 * 8. Example
 *****************************************************************/
/*----------------------- end, example --------------------------*/


