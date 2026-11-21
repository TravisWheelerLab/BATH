/* p7_chkoivx.c — NEON (128-bit) implementations of P7_CHKOIVX management.
 */
#include "p7_config.h"

#include <string.h>

#include <arm_neon.h>

#include "easel.h"
#include "hmmer.h"
#include "impl_neon.h"


/* ck_vec(): first vector of channel <c> of checkpoint <k>. Float
 * channels are c = 0..allocCF-1, int channels follow them.
 */
static inline float32x4_t *
ck_vec(const P7_CHKOIVX *cv, int k, int c)
{
  float32x4_t *base = (float32x4_t *) (((unsigned long int) ((char *) cv->ck_mem + 15)) & (~0xf));
  return base + ((size_t) k * (cv->allocCF + cv->allocCI) + c) * cv->allocQ4;
}

/* set_ivx_ptrs(): point the working int channels into ivx_mem. */
static void
set_ivx_ptrs(P7_CHKOIVX *cv, int CI, int nq)
{
  int c;

  if (CI == 0) return;
  cv->ivx[0] = (int32x4_t *) (((unsigned long int) cv->ivx_mem + 15) & (~0xf));
  for (c = 1; c < CI; c++)
    cv->ivx[c] = cv->ivx[0] + (c * nq);
}


/*****************************************************************
 * 1. The P7_CHKOIVX lifecycle functions.
 *****************************************************************/

/* Function:  p7_chkoivx_Create()
 *
 * Purpose:   Allocate a <P7_CHKOIVX> for a model of up to <M_hint>
 *            nodes using NEON (128-bit) vectors: <CI> working int32
 *            channels, and <K> checkpoints that each hold <CF> float
 *            and <CI> int channels. Channels are p7O_NQF(M_hint)
 *            vectors.
 *
 * Returns:   Pointer to the new <P7_CHKOIVX> on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_CHKOIVX *
p7_chkoivx_Create(int M_hint, int CF, int CI, int K)
{
  P7_CHKOIVX *cv;
  int         nq = p7O_NQF(M_hint);
  int         status;

  cv = NULL;
  ESL_ALLOC(cv, sizeof(P7_CHKOIVX));

  cv->ivx     = NULL;
  cv->ivx_mem = NULL;
  cv->ck_mem  = NULL;
  cv->allocQ4 = 0;

  cv->allocM  = M_hint;
  cv->allocCF = CF;
  cv->allocCI = CI;
  cv->allocK  = K;

  ESL_ALLOC(cv->ivx_mem, sizeof(int32x4_t) * CI * nq + 15);
  ESL_ALLOC(cv->ivx,     sizeof(int32x4_t *) * ESL_MAX(CI, 1));
  ESL_ALLOC(cv->ck_mem,  sizeof(float32x4_t) * (size_t) K * (CF + CI) * nq + 15);

  set_ivx_ptrs(cv, CI, nq);
  cv->allocQ4 = nq;

  return cv;

 ERROR:
  p7_chkoivx_Destroy(cv);
  return NULL;
}


/* Function:  p7_chkoivx_GrowTo()
 *
 * Purpose:   Ensure that <cv> has
 *            storage for at least <M> model nodes, <CF> float and <CI>
 *            int channels, and <K> checkpoints. Returns immediately if
 *            the existing allocation is sufficient. Otherwise frees and
 *            re-allocates the NEON backing memory and re-establishes the
 *            pointer array.
 *
 *            Any data in <cv> must be considered invalid after a
 *            successful reallocation.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; <cv> is destroyed.
 */
int
p7_chkoivx_GrowTo(P7_CHKOIVX *cv, int M, int CF, int CI, int K)
{
  int   nq = p7O_NQF(M);
  void *p;
  int   status;

  if (nq <= cv->allocQ4 && CF <= cv->allocCF && CI <= cv->allocCI && K <= cv->allocK) return eslOK;

  nq = ESL_MAX(nq, cv->allocQ4);
  CF = ESL_MAX(CF, cv->allocCF);
  CI = ESL_MAX(CI, cv->allocCI);
  K  = ESL_MAX(K,  cv->allocK);

  if (cv->ivx_mem != NULL) { free(cv->ivx_mem); cv->ivx_mem = NULL; }
  if (cv->ck_mem  != NULL) { free(cv->ck_mem);  cv->ck_mem  = NULL; }

  if (CI > cv->allocCI) {
    ESL_RALLOC(cv->ivx, p, sizeof(int32x4_t *) * CI);
  }

  ESL_ALLOC(cv->ivx_mem, sizeof(int32x4_t) * CI * nq + 15);
  ESL_ALLOC(cv->ck_mem,  sizeof(float32x4_t) * (size_t) K * (CF + CI) * nq + 15);

  set_ivx_ptrs(cv, CI, nq);

  cv->allocM      = ESL_MAX(M, cv->allocM);
  cv->allocCF     = CF;
  cv->allocCI     = CI;
  cv->allocK      = K;
  cv->allocQ4 = nq;

  return eslOK;

 ERROR:
  p7_chkoivx_Destroy(cv);
  return status;
}


/* Function:  p7_chkoivx_Destroy()
 *
 * Purpose:   Free a <P7_CHKOIVX> that was allocated by
 *            p7_chkoivx_Create(). Safe to call on NULL.
 */
void
p7_chkoivx_Destroy(P7_CHKOIVX *cv)
{
  if (cv == NULL) return;

  if (cv->ivx_mem != NULL) free(cv->ivx_mem);
  if (cv->ivx     != NULL) free(cv->ivx);
  if (cv->ck_mem  != NULL) free(cv->ck_mem);

  free(cv);
}


/*****************************************************************
 * 2. Saving and loading checkpoints.
 *****************************************************************/

/* Function:  p7_chkoivx_Save()
 *
 * Purpose:   Copy the first <allocCF> float channels of <ov> and the
 *            working int channels of <cv> into checkpoint <k>, for a
 *            model of length <M>.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <k>, <M> or <ov> don't fit <cv>.
 */
int
p7_chkoivx_Save(P7_CHKOIVX *cv, const P7_OIVX *ov, int M, int k)
{
  size_t nbytes = sizeof(float32x4_t) * p7O_NQF(M);
  int    c;

  if (k < 0 || k >= cv->allocK || p7O_NQF(M) > cv->allocQ4 || ov->allocC < cv->allocCF || p7O_NQF(M) > ov->allocQ4)
    ESL_EXCEPTION(eslEINVAL, "checkpoint %d for M=%d doesn't fit", k, M);

  for (c = 0; c < cv->allocCF; c++) memcpy(ck_vec(cv, k, c),               ov->ivx[c], nbytes);
  for (c = 0; c < cv->allocCI; c++) memcpy(ck_vec(cv, k, cv->allocCF + c), cv->ivx[c], nbytes);
  return eslOK;
}


/* Function:  p7_chkoivx_Load()
 *
 * Purpose:   Copy checkpoint <k> back into the first <allocCF> float
 *            channels of <ov> and the working int channels of <cv>, for
 *            a model of length <M>.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <k>, <M> or <ov> don't fit <cv>.
 */
int
p7_chkoivx_Load(P7_CHKOIVX *cv, P7_OIVX *ov, int M, int k)
{
  size_t nbytes = sizeof(float32x4_t) * p7O_NQF(M);
  int    c;

  if (k < 0 || k >= cv->allocK || p7O_NQF(M) > cv->allocQ4 || ov->allocC < cv->allocCF || p7O_NQF(M) > ov->allocQ4)
    ESL_EXCEPTION(eslEINVAL, "checkpoint %d for M=%d doesn't fit", k, M);

  for (c = 0; c < cv->allocCF; c++) memcpy(ov->ivx[c], ck_vec(cv, k, c),               nbytes);
  for (c = 0; c < cv->allocCI; c++) memcpy(cv->ivx[c], ck_vec(cv, k, cv->allocCF + c), nbytes);
  return eslOK;
}


/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef p7CHKOIVX_TESTDRIVE
#include "p7_splice.h"

typedef struct {
  P7_CHKOIVX *(*Create) (int, int, int, int);
  int         (*GrowTo) (P7_CHKOIVX *, int, int, int, int);
  int         (*Save)   (P7_CHKOIVX *, const P7_OIVX *, int, int);
  int         (*Load)   (P7_CHKOIVX *, P7_OIVX *, int, int);
  void        (*Destroy)(P7_CHKOIVX *);
  P7_OIVX    *(*OCreate)(int, int);
  void        (*ODestroy)(P7_OIVX *);
  int           V;     /* floats (or int32s) per vector */
  const char   *name;
} CHKOIVX_ISA;

static float *
float_channel(const P7_OIVX *ov, const CHKOIVX_ISA *isa, int c)
{
  return (float *) ov->ivx[c];
}

static int32_t *
int_channel(const P7_CHKOIVX *cv, const CHKOIVX_ISA *isa, int c)
{
  return (int32_t *) cv->ivx[c];
}

/* Fill the working channels with values tagged by <tag>: floats
 * tag*100000 + c*1000 + n, ints -(tag*100000 + c*1000 + n).
 */
static void
write_channels(P7_OIVX *ov, P7_CHKOIVX *cv, const CHKOIVX_ISA *isa, int CF, int CI, int n, int tag)
{
  int c, i;
  for (c = 0; c < CF; c++) { float   *p = float_channel(ov, isa, c); for (i = 0; i < n; i++) p[i] =  (float) (tag * 100000 + c * 1000 + i); }
  for (c = 0; c < CI; c++) { int32_t *p = int_channel  (cv, isa, c); for (i = 0; i < n; i++) p[i] = -(tag * 100000 + c * 1000 + i); }
}

static void
check_channels(const P7_OIVX *ov, const P7_CHKOIVX *cv, const CHKOIVX_ISA *isa, int CF, int CI, int n, int tag)
{
  int c, i;
  for (c = 0; c < CF; c++) { float   *p = float_channel(ov, isa, c); for (i = 0; i < n; i++) if (p[i] !=  (float) (tag * 100000 + c * 1000 + i)) esl_fatal("%s: float channel %d wrong for checkpoint %d", isa->name, c, tag); }
  for (c = 0; c < CI; c++) { int32_t *p = int_channel  (cv, isa, c); for (i = 0; i < n; i++) if (p[i] != -(tag * 100000 + c * 1000 + i)) esl_fatal("%s: int channel %d wrong for checkpoint %d",   isa->name, c, tag); }
}

/* utest_saveload():
 * Save distinct values into each of <K> checkpoints, overwrite the
 * working channels, then load every checkpoint back (last to first)
 * and check its values.
 */
static void
utest_saveload(P7_OIVX *ov, P7_CHKOIVX *cv, const CHKOIVX_ISA *isa, int M, int CF, int CI, int K)
{
  int n = ESL_MAX(2, (M - 1) / isa->V + 1) * isa->V;
  int k;

  for (k = 0; k < K; k++)
    {
      write_channels(ov, cv, isa, CF, CI, n, k + 1);
      if (isa->Save(cv, ov, M, k) != eslOK) esl_fatal("%s: Save(%d) failed", isa->name, k);
    }
  write_channels(ov, cv, isa, CF, CI, n, 99);
  for (k = K - 1; k >= 0; k--)
    {
      if (isa->Load(cv, ov, M, k) != eslOK) esl_fatal("%s: Load(%d) failed", isa->name, k);
      check_channels(ov, cv, isa, CF, CI, n, k + 1);
    }

  if (isa->Save(cv, ov, M, K) != eslEINVAL) esl_fatal("%s: Save past the last checkpoint accepted", isa->name);
}

static void
utest_isa(const CHKOIVX_ISA *isa)
{
  P7_OIVX    *ov = NULL;
  P7_CHKOIVX *cv = NULL;

  if ((ov = isa->OCreate(13, SIGNAL_MEM_SIZE))                       == NULL) esl_fatal("%s: P7_OIVX Create failed", isa->name);
  if ((cv = isa->Create(13, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, 4))    == NULL) esl_fatal("%s: Create failed", isa->name);
  utest_saveload(ov, cv, isa, 13, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, 4);
  utest_saveload(ov, cv, isa,  5, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, 4);   /* smaller M, same allocation */
  isa->ODestroy(ov);

  if ((ov = isa->OCreate(70, SIGNAL_MEM_SIZE))                                    == NULL) esl_fatal("%s: P7_OIVX Create failed", isa->name);
  if (isa->GrowTo(cv, 70, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, 9)                    != eslOK) esl_fatal("%s: GrowTo failed", isa->name);
  utest_saveload(ov, cv, isa, 70, SIGNAL_MEM_SIZE, SIGNAL_MEM_SIZE, 9);
  isa->ODestroy(ov);
  isa->Destroy(cv);
}

#endif /*p7CHKOIVX_TESTDRIVE*/


/*****************************************************************
 * 4. Test driver.
 *****************************************************************/
#ifdef p7CHKOIVX_TESTDRIVE
/*
   make chkoivx_utest
   ./chkoivx_utest
 */
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default   env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE,  NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for P7_CHKOIVX checkpoint vectors";

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);

  esl_exception_SetHandler(&esl_nonfatal_handler);   /* the tests check that bad arguments are refused */

  { CHKOIVX_ISA neon = { p7_chkoivx_Create, p7_chkoivx_GrowTo, p7_chkoivx_Save, p7_chkoivx_Load, p7_chkoivx_Destroy, p7_oivx_Create, p7_oivx_Destroy, 4, "NEON" };
    utest_isa(&neon); }

  esl_getopts_Destroy(go);
  printf("All tests passed.\n");
  return eslOK;
}
#endif /*p7CHKOIVX_TESTDRIVE*/
