/* p7_fs_oprofile_avx.c — AVX2 (256-bit) implementations of P7_FS_OPROFILE management.
 *
 * All public functions carry the _avx suffix and are selected at runtime
 * by impl_Init() in p7_oprofile.c when AVX2 is the fastest available ISA.
 *
 * Ported from p7_fs_oprofile_sse.c with:
 *   - Public function names suffixed with _avx
 *   - __m128 / p7O_NQF / allocQ4 / rfv/tfv/rfv_mem/tfv_mem
 *       → __m256 / p7O_NQF_AVX / allocQ4_avx / rfv_avx/tfv_avx/rfv_mem_avx/tfv_mem_avx
 *   - 4-lane union → 8-lane union
 *   - esl_sse_expf/logf → local avx_expf/avx_logf helpers
 *   - Memory alignment padding +15/~0xf → +31/~0x1f (32-byte AVX alignment)
 *   - sizeof(__m128) → sizeof(__m256)
 */
#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <stdio.h>
#include <string.h>
#include <math.h>

#include <immintrin.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_sse.h"       /* esl_sse_expf / esl_sse_logf used by avx helpers */
#include "esl_avx.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"

/* Forward declarations of internal static helpers */
static int fs_fb_conversion    (const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs);
static int fs_fb_conversion_log(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs);

/* exp/log on __m256: apply esl_sse_expf/logf to the two 128-bit halves */
static inline __m256
avx_expf(__m256 v)
{
  __m128 lo = esl_sse_expf(_mm256_castps256_ps128(v));
  __m128 hi = esl_sse_expf(_mm256_extractf128_ps(v, 1));
  return _mm256_set_m128(hi, lo);
}

static inline __m256
avx_logf(__m256 v)
{
  __m128 lo = esl_sse_logf(_mm256_castps256_ps128(v));
  __m128 hi = esl_sse_logf(_mm256_extractf128_ps(v, 1));
  return _mm256_set_m128(hi, lo);
}


/*****************************************************************
 * 1. P7_FS_OPROFILE: allocation, initialization, destruction.
 *****************************************************************/

/* Function:  p7_fs_oprofile_Create_avx()
 * Synopsis:  Allocate an optimized frameshift profile structure (AVX2).
 */
P7_FS_OPROFILE *
p7_fs_oprofile_Create_avx(int allocM, const ESL_ALPHABET *abc, int codon_lengths)
{
  int             status;
  P7_FS_OPROFILE *om_fs = NULL;
  int             nqf;
  int             ncodon_rows;
  int             c;

  if      (codon_lengths == 1) ncodon_rows = p7P_MAXCODONS1 + abc->Kp;
  else if (codon_lengths == 3) ncodon_rows = p7P_MAXCODONS3 + abc->Kp;
  else if (codon_lengths == 5) ncodon_rows = p7P_MAXCODONS5 + abc->Kp;
  else    ESL_XEXCEPTION(eslEINVAL, "codon_lengths must be 1, 3 or 5");

  nqf = p7O_NQF_AVX(allocM);

  ESL_ALLOC(om_fs, sizeof(P7_FS_OPROFILE));

  /* SSE fields — not used on the AVX path */
  om_fs->rfv_mem = NULL;
  om_fs->tfv_mem = NULL;
  om_fs->rfv     = NULL;
  om_fs->tfv     = NULL;
  om_fs->allocQ4 = 0;

  /* AVX fields */
  om_fs->rfv_mem_avx = NULL;
  om_fs->tfv_mem_avx = NULL;
  om_fs->rfv_avx     = NULL;
  om_fs->tfv_avx     = NULL;
  om_fs->allocQ4_avx = 0;

#ifdef eslENABLE_AVX512
  om_fs->rfv_mem_avx512 = NULL;
  om_fs->tfv_mem_avx512 = NULL;
  om_fs->rfv_avx512     = NULL;
  om_fs->tfv_avx512     = NULL;
  om_fs->allocQ4_avx512 = 0;
#endif

  om_fs->clone = 0;

  /* Raw memory (+31 bytes for manual 32-byte alignment) */
  ESL_ALLOC(om_fs->rfv_mem_avx, sizeof(__m256) * nqf * ncodon_rows + 31);
  ESL_ALLOC(om_fs->tfv_mem_avx, sizeof(__m256) * nqf * p7O_NTRANS  + 31);

  /* Pointer array: one __m256 * per emission row */
  ESL_ALLOC(om_fs->rfv_avx, sizeof(__m256 *) * ncodon_rows);

  /* Align to 32-byte boundaries */
  om_fs->rfv_avx[0] = (__m256 *) (((unsigned long int) om_fs->rfv_mem_avx + 31) & (~0x1f));
  om_fs->tfv_avx    = (__m256 *) (((unsigned long int) om_fs->tfv_mem_avx + 31) & (~0x1f));

  /* Set remaining row pointers */
  for (c = 1; c < ncodon_rows; c++)
    om_fs->rfv_avx[c] = om_fs->rfv_avx[0] + (c * nqf);

  om_fs->allocQ4_avx = nqf;

  /* Frameshift-specific initializations */
  om_fs->codon_lengths = codon_lengths;
  om_fs->fsprob        = 0.0f;

  for (c = 0; c < p7_NOFFSETS; c++) om_fs->offs[c]    = -1;
  for (c = 0; c < p7_NEVPARAM; c++) om_fs->evparam[c] = p7_EVPARAM_UNSET;
  for (c = 0; c < p7_NCUTOFFS; c++) om_fs->cutoff[c]  = p7_CUTOFF_UNSET;
  for (c = 0; c < p7_MAXABET;  c++) om_fs->compo[c]   = p7_COMPO_UNSET;

  om_fs->name = NULL;
  om_fs->acc  = NULL;
  om_fs->desc = NULL;

  ESL_ALLOC(om_fs->rf,        sizeof(char) * (allocM + 2));
  ESL_ALLOC(om_fs->mm,        sizeof(char) * (allocM + 2));
  ESL_ALLOC(om_fs->cs,        sizeof(char) * (allocM + 2));
  ESL_ALLOC(om_fs->consensus, sizeof(char) * (allocM + 2));
  memset(om_fs->rf,        '\0', sizeof(char) * (allocM + 2));
  memset(om_fs->mm,        '\0', sizeof(char) * (allocM + 2));
  memset(om_fs->cs,        '\0', sizeof(char) * (allocM + 2));
  memset(om_fs->consensus, '\0', sizeof(char) * (allocM + 2));

  om_fs->abc        = abc;
  om_fs->L          = 0;
  om_fs->M          = 0;
  om_fs->max_length = -1;
  om_fs->allocM     = allocM;
  om_fs->mode       = p7_NO_MODE;
  om_fs->nj         = 0.0f;
  return om_fs;

 ERROR:
  p7_fs_oprofile_Destroy_avx(om_fs);
  return NULL;
}


/* Function:  p7_fs_oprofile_Destroy_avx()
 * Synopsis:  Free an optimized frameshift profile (AVX fields only).
 */
void
p7_fs_oprofile_Destroy_avx(P7_FS_OPROFILE *om_fs)
{
  if (om_fs == NULL) return;
  if (om_fs->clone == 0)
    {
      if (om_fs->rfv_mem_avx != NULL) free(om_fs->rfv_mem_avx);
      if (om_fs->tfv_mem_avx != NULL) free(om_fs->tfv_mem_avx);
      if (om_fs->rfv_avx     != NULL) free(om_fs->rfv_avx);
      if (om_fs->name        != NULL) free(om_fs->name);
      if (om_fs->acc         != NULL) free(om_fs->acc);
      if (om_fs->desc        != NULL) free(om_fs->desc);
      if (om_fs->rf          != NULL) free(om_fs->rf);
      if (om_fs->mm          != NULL) free(om_fs->mm);
      if (om_fs->cs          != NULL) free(om_fs->cs);
      if (om_fs->consensus   != NULL) free(om_fs->consensus);
    }
  free(om_fs);
}


/* Function:  p7_fs_oprofile_Clone_avx()
 * Synopsis:  Shallow clone for multithreading (AVX2).
 */
P7_FS_OPROFILE *
p7_fs_oprofile_Clone_avx(const P7_FS_OPROFILE *om_fs)
{
  int             status;
  P7_FS_OPROFILE *om2 = NULL;

  ESL_ALLOC(om2, sizeof(P7_FS_OPROFILE));
  memcpy(om2, om_fs, sizeof(P7_FS_OPROFILE));
  om2->clone = 1;
  return om2;

 ERROR:
  p7_fs_oprofile_Destroy_avx(om2);
  return NULL;
}


/*****************************************************************
 * 2. Conversion from P7_FS_PROFILE to P7_FS_OPROFILE (AVX2).
 *****************************************************************/

static int
fs_fb_conversion(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs)
{
  int   M  = gm_fs->M;
  int   nq = p7O_NQF_AVX(M);
  int   ncodon_rows;
  int   c, q, k, kb, z, t, tg, j;
  union { __m256 v; float x[8]; } tmp;

  if (nq > om_fs->allocQ4_avx) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  if      (om_fs->codon_lengths == 1) ncodon_rows = p7P_MAXCODONS1 + gm_fs->abc->Kp;
  else if (om_fs->codon_lengths == 3) ncodon_rows = p7P_MAXCODONS3 + gm_fs->abc->Kp;
  else if (om_fs->codon_lengths == 5) ncodon_rows = p7P_MAXCODONS5 + gm_fs->abc->Kp;
  else ESL_EXCEPTION(eslEINVAL, "codon_lengths must be 1, 3, or 5");

  /* Striped match emission scores: odds ratios */
  for (c = 0; c < ncodon_rows; c++)
    for (k = 1, q = 0; q < nq; q++, k++) {
      for (z = 0; z < 8; z++) tmp.x[z] = (k + z*nq <= M) ? gm_fs->rsc[c][k + z*nq] : -eslINFINITY;
      om_fs->rfv_avx[c][q] = avx_expf(tmp.v);
    }

  /* Transition scores (all but DD) */
  for (j = 0, k = 1, q = 0; q < nq; q++, k++) {
    for (t = p7O_BM; t <= p7O_II; t++) {
      switch (t) {
      case p7O_BM: tg = p7P_BM; kb = k-1; break;
      case p7O_MM: tg = p7P_MM; kb = k-1; break;
      case p7O_IM: tg = p7P_IM; kb = k-1; break;
      case p7O_DM: tg = p7P_DM; kb = k-1; break;
      case p7O_MD: tg = p7P_MD; kb = k;   break;
      case p7O_MI: tg = p7P_MI; kb = k;   break;
      case p7O_II: tg = p7P_II; kb = k;   break;
      default:     tg = 0;      kb = k;   break;
      }
      for (z = 0; z < 8; z++) tmp.x[z] = (kb + z*nq < M) ? p7P_TSC(gm_fs, kb + z*nq, tg) : -eslINFINITY;
      om_fs->tfv_avx[j++] = avx_expf(tmp.v);
    }
  }

  /* DD transitions */
  for (k = 1, q = 0; q < nq; q++, k++) {
    for (z = 0; z < 8; z++) tmp.x[z] = (k + z*nq < M) ? p7P_TSC(gm_fs, k + z*nq, p7P_DD) : -eslINFINITY;
    om_fs->tfv_avx[j++] = avx_expf(tmp.v);
  }

  /* Special state (ENJC) transition costs: pspace odds ratios */
  om_fs->xf[p7O_E][p7O_LOOP] = expf(gm_fs->xsc[p7P_E][p7P_LOOP]);
  om_fs->xf[p7O_E][p7O_MOVE] = expf(gm_fs->xsc[p7P_E][p7P_MOVE]);
  om_fs->xf[p7O_N][p7O_LOOP] = expf(gm_fs->xsc[p7P_N][p7P_LOOP]);
  om_fs->xf[p7O_N][p7O_MOVE] = expf(gm_fs->xsc[p7P_N][p7P_MOVE]);
  om_fs->xf[p7O_C][p7O_LOOP] = expf(gm_fs->xsc[p7P_C][p7P_LOOP]);
  om_fs->xf[p7O_C][p7O_MOVE] = expf(gm_fs->xsc[p7P_C][p7P_MOVE]);
  om_fs->xf[p7O_J][p7O_LOOP] = expf(gm_fs->xsc[p7P_J][p7P_LOOP]);
  om_fs->xf[p7O_J][p7O_MOVE] = expf(gm_fs->xsc[p7P_J][p7P_MOVE]);
  return eslOK;
}

static int
fs_fb_conversion_log(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs)
{
  int   M  = gm_fs->M;
  int   nq = p7O_NQF_AVX(M);
  int   ncodon_rows;
  int   c, q, k, kb, z, t, tg, j;
  union { __m256 v; float x[8]; } tmp;

  if (nq > om_fs->allocQ4_avx) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  if      (om_fs->codon_lengths == 1) ncodon_rows = p7P_MAXCODONS1 + gm_fs->abc->Kp;
  else if (om_fs->codon_lengths == 3) ncodon_rows = p7P_MAXCODONS3 + gm_fs->abc->Kp;
  else if (om_fs->codon_lengths == 5) ncodon_rows = p7P_MAXCODONS5 + gm_fs->abc->Kp;
  else ESL_EXCEPTION(eslEINVAL, "codon_lengths must be 1, 3, or 5");

  /* Striped match emission scores in log-space: copy directly */
  for (c = 0; c < ncodon_rows; c++)
    for (k = 1, q = 0; q < nq; q++, k++) {
      for (z = 0; z < 8; z++) tmp.x[z] = (k + z*nq <= M) ? gm_fs->rsc[c][k + z*nq] : -eslINFINITY;
      om_fs->rfv_avx[c][q] = tmp.v;
    }

  /* Transition scores in log-space */
  for (j = 0, k = 1, q = 0; q < nq; q++, k++) {
    for (t = p7O_BM; t <= p7O_II; t++) {
      switch (t) {
      case p7O_BM: tg = p7P_BM; kb = k-1; break;
      case p7O_MM: tg = p7P_MM; kb = k-1; break;
      case p7O_IM: tg = p7P_IM; kb = k-1; break;
      case p7O_DM: tg = p7P_DM; kb = k-1; break;
      case p7O_MD: tg = p7P_MD; kb = k;   break;
      case p7O_MI: tg = p7P_MI; kb = k;   break;
      case p7O_II: tg = p7P_II; kb = k;   break;
      default:     tg = 0;      kb = k;   break;
      }
      for (z = 0; z < 8; z++) tmp.x[z] = (kb + z*nq < M) ? p7P_TSC(gm_fs, kb + z*nq, tg) : -eslINFINITY;
      om_fs->tfv_avx[j++] = tmp.v;
    }
  }

  /* DD transitions in log-space */
  for (k = 1, q = 0; q < nq; q++, k++) {
    for (z = 0; z < 8; z++) tmp.x[z] = (k + z*nq < M) ? p7P_TSC(gm_fs, k + z*nq, p7P_DD) : -eslINFINITY;
    om_fs->tfv_avx[j++] = tmp.v;
  }

  /* Special state transitions in log-space */
  om_fs->xf[p7O_E][p7O_LOOP] = gm_fs->xsc[p7P_E][p7P_LOOP];
  om_fs->xf[p7O_E][p7O_MOVE] = gm_fs->xsc[p7P_E][p7P_MOVE];
  om_fs->xf[p7O_N][p7O_LOOP] = gm_fs->xsc[p7P_N][p7P_LOOP];
  om_fs->xf[p7O_N][p7O_MOVE] = gm_fs->xsc[p7P_N][p7P_MOVE];
  om_fs->xf[p7O_C][p7O_LOOP] = gm_fs->xsc[p7P_C][p7P_LOOP];
  om_fs->xf[p7O_C][p7O_MOVE] = gm_fs->xsc[p7P_C][p7P_MOVE];
  om_fs->xf[p7O_J][p7O_LOOP] = gm_fs->xsc[p7P_J][p7P_LOOP];
  om_fs->xf[p7O_J][p7O_MOVE] = gm_fs->xsc[p7P_J][p7P_MOVE];
  return eslOK;
}


/* Function:  p7_fs_oprofile_Convert_avx()
 */
int
p7_fs_oprofile_Convert_avx(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs)
{
  int status, z;

  if (gm_fs->abc->type != om_fs->abc->type) ESL_EXCEPTION(eslEINVAL, "alphabets don't match");
  if (gm_fs->M         >  om_fs->allocM)    ESL_EXCEPTION(eslEINVAL, "optimized profile is too small");

  om_fs->mode          = gm_fs->mode;
  om_fs->L             = gm_fs->L;
  om_fs->M             = gm_fs->M;
  om_fs->nj            = gm_fs->nj;
  om_fs->max_length    = gm_fs->max_length;
  om_fs->codon_lengths = gm_fs->codon_lengths;
  om_fs->fsprob        = gm_fs->fsprob;

  if ((status = fs_fb_conversion(gm_fs, om_fs)) != eslOK) return status;

  if (om_fs->name != NULL) free(om_fs->name);
  if (om_fs->acc  != NULL) free(om_fs->acc);
  if (om_fs->desc != NULL) free(om_fs->desc);
  if ((status = esl_strdup(gm_fs->name, -1, &om_fs->name)) != eslOK) goto ERROR;
  if ((status = esl_strdup(gm_fs->acc,  -1, &om_fs->acc))  != eslOK) goto ERROR;
  if ((status = esl_strdup(gm_fs->desc, -1, &om_fs->desc)) != eslOK) goto ERROR;
  strcpy(om_fs->rf,        gm_fs->rf);
  strcpy(om_fs->mm,        gm_fs->mm);
  strcpy(om_fs->cs,        gm_fs->cs);
  strcpy(om_fs->consensus, gm_fs->consensus);
  for (z = 0; z < p7_NEVPARAM; z++) om_fs->evparam[z] = gm_fs->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) om_fs->cutoff[z]  = gm_fs->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) om_fs->compo[z]   = gm_fs->compo[z];
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_fs_oprofile_Convert_Log_avx()
 */
int
p7_fs_oprofile_Convert_Log_avx(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs)
{
  int status, z;

  if (gm_fs->abc->type != om_fs->abc->type) ESL_EXCEPTION(eslEINVAL, "alphabets don't match");
  if (gm_fs->M         >  om_fs->allocM)    ESL_EXCEPTION(eslEINVAL, "optimized profile is too small");

  om_fs->mode          = gm_fs->mode;
  om_fs->L             = gm_fs->L;
  om_fs->M             = gm_fs->M;
  om_fs->nj            = gm_fs->nj;
  om_fs->max_length    = gm_fs->max_length;
  om_fs->codon_lengths = gm_fs->codon_lengths;
  om_fs->fsprob        = gm_fs->fsprob;

  if ((status = fs_fb_conversion_log(gm_fs, om_fs)) != eslOK) return status;

  if (om_fs->name != NULL) free(om_fs->name);
  if (om_fs->acc  != NULL) free(om_fs->acc);
  if (om_fs->desc != NULL) free(om_fs->desc);
  if ((status = esl_strdup(gm_fs->name, -1, &om_fs->name)) != eslOK) goto ERROR;
  if ((status = esl_strdup(gm_fs->acc,  -1, &om_fs->acc))  != eslOK) goto ERROR;
  if ((status = esl_strdup(gm_fs->desc, -1, &om_fs->desc)) != eslOK) goto ERROR;
  strcpy(om_fs->rf,        gm_fs->rf);
  strcpy(om_fs->mm,        gm_fs->mm);
  strcpy(om_fs->cs,        gm_fs->cs);
  strcpy(om_fs->consensus, gm_fs->consensus);
  for (z = 0; z < p7_NEVPARAM; z++) om_fs->evparam[z] = gm_fs->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) om_fs->cutoff[z]  = gm_fs->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) om_fs->compo[z]   = gm_fs->compo[z];
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_fs_oprofile_SubConvert_Log_avx()
 */
int
p7_fs_oprofile_SubConvert_Log_avx(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs, int k_start, int k_end)
{
  int   Mp = k_end - k_start + 1;
  int   nq = p7O_NQF_AVX(Mp);
  int   ncodon_rows;
  int   c, q, k, kb, z, t, tg, j;
  union { __m256 v; float x[8]; } tmp;

  if (nq > om_fs->allocQ4_avx) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold sub-region");

  if      (om_fs->codon_lengths == 1) ncodon_rows = p7P_MAXCODONS1 + gm_fs->abc->Kp;
  else if (om_fs->codon_lengths == 3) ncodon_rows = p7P_MAXCODONS3 + gm_fs->abc->Kp;
  else if (om_fs->codon_lengths == 5) ncodon_rows = p7P_MAXCODONS5 + gm_fs->abc->Kp;
  else ESL_EXCEPTION(eslEINVAL, "codon_lengths must be 1, 3, or 5");

  /* Re-stripe emission scores for sub-region */
  for (c = 0; c < ncodon_rows; c++)
    for (k = 1, q = 0; q < nq; q++, k++) {
      for (z = 0; z < 8; z++)
        tmp.x[z] = (k + z*nq <= Mp) ? gm_fs->rsc[c][k_start - 1 + k + z*nq] : -eslINFINITY;
      om_fs->rfv_avx[c][q] = tmp.v;
    }

  /* Re-stripe transition scores */
  for (j = 0, k = 1, q = 0; q < nq; q++, k++) {
    for (t = p7O_BM; t <= p7O_II; t++) {
      switch (t) {
      case p7O_BM: tg = p7P_BM; kb = k-1; break;
      case p7O_MM: tg = p7P_MM; kb = k-1; break;
      case p7O_IM: tg = p7P_IM; kb = k-1; break;
      case p7O_DM: tg = p7P_DM; kb = k-1; break;
      case p7O_MD: tg = p7P_MD; kb = k;   break;
      case p7O_MI: tg = p7P_MI; kb = k;   break;
      case p7O_II: tg = p7P_II; kb = k;   break;
      default:     tg = 0;      kb = k;   break;
      }
      for (z = 0; z < 8; z++)
        tmp.x[z] = (kb + z*nq < Mp) ? p7P_TSC(gm_fs, k_start - 1 + kb + z*nq, tg) : -eslINFINITY;
      om_fs->tfv_avx[j++] = tmp.v;
    }
  }

  /* DD transitions */
  for (k = 1, q = 0; q < nq; q++, k++) {
    for (z = 0; z < 8; z++)
      tmp.x[z] = (k + z*nq < Mp) ? p7P_TSC(gm_fs, k_start - 1 + k + z*nq, p7P_DD) : -eslINFINITY;
    om_fs->tfv_avx[j++] = tmp.v;
  }

  om_fs->M = Mp;

  om_fs->xf[p7O_E][p7O_LOOP] = gm_fs->xsc[p7P_E][p7P_LOOP];
  om_fs->xf[p7O_E][p7O_MOVE] = gm_fs->xsc[p7P_E][p7P_MOVE];

  return eslOK;
}


/*****************************************************************
 * 3. Reconfig functions (AVX2).
 *****************************************************************/

int
p7_fs_oprofile_ReconfigLength_avx(P7_FS_OPROFILE *om_fs, int L)
{
  float pmove = (2.0f + om_fs->nj) / ((float)L + 2.0f + om_fs->nj);
  float ploop = 1.0f - pmove;

  om_fs->xf[p7O_N][p7O_LOOP] = om_fs->xf[p7O_C][p7O_LOOP] = om_fs->xf[p7O_J][p7O_LOOP] = ploop;
  om_fs->xf[p7O_N][p7O_MOVE] = om_fs->xf[p7O_C][p7O_MOVE] = om_fs->xf[p7O_J][p7O_MOVE] = pmove;

  om_fs->L = L;
  return eslOK;
}

int
p7_fs_oprofile_ReconfigLength_Log_avx(P7_FS_OPROFILE *om_fs, int L)
{
  float pmove = (2.0f + om_fs->nj) / ((float)L + 2.0f + om_fs->nj);
  float ploop = 1.0f - pmove;

  om_fs->xf[p7O_N][p7O_LOOP] = om_fs->xf[p7O_C][p7O_LOOP] = om_fs->xf[p7O_J][p7O_LOOP] = logf(ploop);
  om_fs->xf[p7O_N][p7O_MOVE] = om_fs->xf[p7O_C][p7O_MOVE] = om_fs->xf[p7O_J][p7O_MOVE] = logf(pmove);

  om_fs->L = L;
  return eslOK;
}

int
p7_fs_oprofile_ReconfigMultihit_avx(P7_FS_OPROFILE *om_fs, int L)
{
  om_fs->xf[p7O_E][p7O_MOVE] = 0.5f;
  om_fs->xf[p7O_E][p7O_LOOP] = 0.5f;
  om_fs->nj = 1.0f;
  return p7_fs_oprofile_ReconfigLength_avx(om_fs, L);
}

int
p7_fs_oprofile_ReconfigUnihit_avx(P7_FS_OPROFILE *om_fs, int L)
{
  om_fs->xf[p7O_E][p7O_MOVE] = 1.0f;
  om_fs->xf[p7O_E][p7O_LOOP] = 0.0f;
  om_fs->nj = 0.0f;
  return p7_fs_oprofile_ReconfigLength_avx(om_fs, L);
}


/*****************************************************************
 * 4. Logify (AVX2).
 *****************************************************************/

int
p7_fs_oprofile_Logify_avx(P7_FS_OPROFILE *om_fs)
{
  int Q = p7O_NQF_AVX(om_fs->M);
  int ncodon_rows;
  int c, j, s, t;

  if      (om_fs->codon_lengths == 1) ncodon_rows = p7P_MAXCODONS1 + om_fs->abc->Kp;
  else if (om_fs->codon_lengths == 3) ncodon_rows = p7P_MAXCODONS3 + om_fs->abc->Kp;
  else if (om_fs->codon_lengths == 5) ncodon_rows = p7P_MAXCODONS5 + om_fs->abc->Kp;
  else ESL_EXCEPTION(eslEINVAL, "codon_lengths must be 1, 3, or 5");

  for (c = 0; c < ncodon_rows; c++)
    for (j = 0; j < Q; j++)
      om_fs->rfv_avx[c][j] = avx_logf(om_fs->rfv_avx[c][j]);

  for (j = 0; j < p7O_NTRANS * Q; j++)
    om_fs->tfv_avx[j] = avx_logf(om_fs->tfv_avx[j]);

  for (s = 0; s < p7O_NXSTATES; s++)
    for (t = 0; t < p7O_NXTRANS; t++)
      om_fs->xf[s][t] = logf(om_fs->xf[s][t]);

  /* om_fs->fsprob is already log-odds; leave unchanged. */
  return eslOK;
}

#endif /* eslENABLE_AVX */
