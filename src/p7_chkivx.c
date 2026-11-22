/* P7_CHKIVX implementation: values saved at each checkpoint of a P7_CHKGMX
 *
 * Contents:
 *   1. The <P7_CHKIVX> object
 *   2. Saving and loading checkpoints
 *   3. Unit tests
 *   4. Test driver
 */
#include "p7_config.h"

#include <string.h>

#include "easel.h"

#include "hmmer.h"

/*****************************************************************
 * 1. The <P7_CHKIVX> object.
 *****************************************************************/

/* Function:  p7_chkivx_Create()
 * Synopsis:  Allocate a new <P7_CHKIVX>.
 *
 * Purpose:   Allocate a reusable, resizeable <P7_CHKIVX> for models up
 *            to size <allocM>: <CI> working int channels, and <K>
 *            checkpoints that each hold <CF> float and <CI> int
 *            channels.
 *
 * Returns:   a pointer to the new <P7_CHKIVX>.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_CHKIVX *
p7_chkivx_Create(int allocM, int CF, int CI, int K)
{
  P7_CHKIVX *cv = NULL;
  int        status;

  ESL_ALLOC(cv, sizeof(P7_CHKIVX));
  cv->ivx = NULL;
  cv->ckf = NULL;
  cv->cki = NULL;

  ESL_ALLOC(cv->ivx, sizeof(int)   * ESL_MAX(1, (allocM+1) * CI));
  ESL_ALLOC(cv->ckf, sizeof(float) * ESL_MAX(1, (size_t) K * (allocM+1) * CF));
  ESL_ALLOC(cv->cki, sizeof(int)   * ESL_MAX(1, (size_t) K * (allocM+1) * CI));

  cv->allocM  = allocM;
  cv->allocCF = CF;
  cv->allocCI = CI;
  cv->allocK  = K;

  return cv;

 ERROR:
  if (cv != NULL) p7_chkivx_Destroy(cv);
  return NULL;
}


/* Function:  p7_chkivx_GrowTo()
 * Synopsis:  Assure that checkpoint storage is big enough.
 *
 * Purpose:   Assures that <cv> has room for models up to size <M>, <CF>
 *            float and <CI> int channels, and <K> checkpoints;
 *            reallocates if necessary.
 *
 * Returns:   <eslOK> on success, and <cv> may be reallocated upon
 *            return; any data that may have been in <cv> must be
 *            assumed to be invalidated.
 *
 * Throws:    <eslEMEM> on allocation failure, and any data that may
 *            have been in <cv> must be assumed to be invalidated.
 */
int
p7_chkivx_GrowTo(P7_CHKIVX *cv, int M, int CF, int CI, int K)
{
  void *p;
  int   status;

  if (M <= cv->allocM && CF <= cv->allocCF && CI <= cv->allocCI && K <= cv->allocK) return eslOK;

  M  = ESL_MAX(M,  cv->allocM);
  CF = ESL_MAX(CF, cv->allocCF);
  CI = ESL_MAX(CI, cv->allocCI);
  K  = ESL_MAX(K,  cv->allocK);

  ESL_RALLOC(cv->ivx, p, sizeof(int)   * ESL_MAX(1, (M+1) * CI));
  ESL_RALLOC(cv->ckf, p, sizeof(float) * ESL_MAX(1, (size_t) K * (M+1) * CF));
  ESL_RALLOC(cv->cki, p, sizeof(int)   * ESL_MAX(1, (size_t) K * (M+1) * CI));

  cv->allocM  = M;
  cv->allocCF = CF;
  cv->allocCI = CI;
  cv->allocK  = K;

  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_chkivx_Destroy()
 * Synopsis:  Frees checkpoint storage.
 *
 * Purpose:   Frees a <P7_CHKIVX>. Safe to call on NULL.
 *
 * Returns:   (void)
 */
void
p7_chkivx_Destroy(P7_CHKIVX *cv)
{
  if (cv == NULL) return;

  if (cv->ivx != NULL) free(cv->ivx);
  if (cv->ckf != NULL) free(cv->ckf);
  if (cv->cki != NULL) free(cv->cki);
  free(cv);
  return;
}


/*****************************************************************
 * 2. Saving and loading checkpoints.
 *****************************************************************/

/* Function:  p7_chkivx_Save()
 * Synopsis:  Save working values into a checkpoint.
 *
 * Purpose:   Copy the first (<M>+1)*<allocCF> floats of <iv> (model
 *            positions 0..M of <allocCF> channels, in P7_IVX layout) and
 *            the working ints of <cv> for model positions 0..M into
 *            checkpoint <k>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <k>, <M> or <iv> don't fit <cv>.
 */
int
p7_chkivx_Save(P7_CHKIVX *cv, const P7_IVX *iv, int M, int k)
{
  if (k < 0 || k >= cv->allocK || M > cv->allocM || (M+1) * cv->allocCF > (iv->allocM+1) * iv->allocC)
    ESL_EXCEPTION(eslEINVAL, "checkpoint %d for M=%d doesn't fit", k, M);

  memcpy(cv->ckf + (size_t) k * (cv->allocM+1) * cv->allocCF, iv->ivx, sizeof(float) * (M+1) * cv->allocCF);
  memcpy(cv->cki + (size_t) k * (cv->allocM+1) * cv->allocCI, cv->ivx, sizeof(int)   * (M+1) * cv->allocCI);
  return eslOK;
}


/* Function:  p7_chkivx_Load()
 * Synopsis:  Restore working values from a checkpoint.
 *
 * Purpose:   Copy checkpoint <k> back into the first (<M>+1)*<allocCF>
 *            floats of <iv> and the working ints of <cv>, for model
 *            positions 0..M.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <k>, <M> or <iv> don't fit <cv>.
 */
int
p7_chkivx_Load(P7_CHKIVX *cv, P7_IVX *iv, int M, int k)
{
  if (k < 0 || k >= cv->allocK || M > cv->allocM || (M+1) * cv->allocCF > (iv->allocM+1) * iv->allocC)
    ESL_EXCEPTION(eslEINVAL, "checkpoint %d for M=%d doesn't fit", k, M);

  memcpy(iv->ivx, cv->ckf + (size_t) k * (cv->allocM+1) * cv->allocCF, sizeof(float) * (M+1) * cv->allocCF);
  memcpy(cv->ivx, cv->cki + (size_t) k * (cv->allocM+1) * cv->allocCI, sizeof(int)   * (M+1) * cv->allocCI);
  return eslOK;
}


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7CHKIVX_TESTDRIVE
#include "p7_splice.h"

/* Fill the working values with values tagged by <tag>: floats
 * tag*100000 + n, ints -(tag*100000 + n).
 */
static void
write_values(P7_IVX *iv, P7_CHKIVX *cv, int M, int tag)
{
  int n;
  for (n = 0; n < (M+1) * cv->allocCF; n++) iv->ivx[n] =  (float) (tag * 100000 + n);
  for (n = 0; n < (M+1) * cv->allocCI; n++) cv->ivx[n] = -(tag * 100000 + n);
}

static void
check_values(const P7_IVX *iv, const P7_CHKIVX *cv, int M, int tag)
{
  int n;
  for (n = 0; n < (M+1) * cv->allocCF; n++) if (iv->ivx[n] !=  (float) (tag * 100000 + n)) esl_fatal("float value %d wrong for checkpoint %d", n, tag);
  for (n = 0; n < (M+1) * cv->allocCI; n++) if (cv->ivx[n] != -(tag * 100000 + n))         esl_fatal("int value %d wrong for checkpoint %d",   n, tag);
}

/* utest_saveload():
 * Save distinct values into each of <K> checkpoints, overwrite the
 * working values, then load every checkpoint back (last to first)
 * and check its values.
 */
static void
utest_saveload(P7_IVX *iv, P7_CHKIVX *cv, int M, int K)
{
  int k;

  for (k = 0; k < K; k++)
    {
      write_values(iv, cv, M, k + 1);
      if (p7_chkivx_Save(cv, iv, M, k) != eslOK) esl_fatal("Save(%d) failed", k);
    }
  write_values(iv, cv, M, 99);
  for (k = K - 1; k >= 0; k--)
    {
      if (p7_chkivx_Load(cv, iv, M, k) != eslOK) esl_fatal("Load(%d) failed", k);
      check_values(iv, cv, M, k + 1);
    }

  if (p7_chkivx_Save(cv, iv, M, K) != eslEINVAL) esl_fatal("Save past the last checkpoint accepted");
}

static void
utest_chkivx(void)
{
  P7_IVX    *iv = NULL;
  P7_CHKIVX *cv = NULL;

  if ((iv = p7_ivx_Create(13, SIGNAL_MEM_SIZE))                         == NULL) esl_fatal("P7_IVX Create failed");
  if ((cv = p7_chkivx_Create(13, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, 4))  == NULL) esl_fatal("Create failed");
  utest_saveload(iv, cv, 13, 4);
  utest_saveload(iv, cv,  5, 4);                                 /* smaller M, same allocation */

  if (p7_ivx_GrowTo(iv, 70, SIGNAL_MEM_SIZE)                          != eslOK) esl_fatal("P7_IVX GrowTo failed");
  if (p7_chkivx_GrowTo(cv, 70, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, 9)   != eslOK) esl_fatal("GrowTo failed");
  utest_saveload(iv, cv, 70, 9);

  p7_ivx_Destroy(iv);
  p7_chkivx_Destroy(cv);
}

#endif /*p7CHKIVX_TESTDRIVE*/


/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7CHKIVX_TESTDRIVE
/*
   gcc -g -Wall -std=gnu99 -o p7_chkivx_utest -I. -L. -I../easel -L../easel -Dp7CHKIVX_TESTDRIVE p7_chkivx.c -lhmmer -leasel -lm
   ./p7_chkivx_utest
 */
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default   env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for P7_CHKIVX checkpoint values";

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);

  esl_exception_SetHandler(&esl_nonfatal_handler);   /* the tests check that bad arguments are refused */
  utest_chkivx();

  esl_getopts_Destroy(go);
  printf("All tests passed.\n");
  return eslOK;
}
#endif /*p7CHKIVX_TESTDRIVE*/
