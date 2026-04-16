/* p7_oprofile_avx.c — AVX2 (256-bit) implementations of P7_OPROFILE management.
 *
 * All public functions carry the _avx suffix and are selected at runtime
 * by impl_Init() in p7_oprofile.c when AVX2 is the fastest available ISA.
 *
 * Systematic differences from p7_oprofile_sse.c:
 *   - __m128i / __m128  →  __m256i / __m256
 *   - p7O_NQB/NQW/NQF   →  p7O_NQB_AVX / p7O_NQW_AVX / p7O_NQF_AVX
 *   - rbv/sbv/rwv/twv/rfv/tfv  →  *_avx counterparts
 *   - rbv_mem/… → rbv_mem_avx/…
 *   - allocQ16/8/4 → allocQ16_avx/allocQ8_avx/allocQ4_avx
 *   - alignment padding +15/~0xf → +31/~0x1f  (32-byte AVX alignment)
 *   - lane counts: 16/8/4 → 32/16/8
 *   - _mm_*  →  _mm256_*  intrinsics
 *   - esl_sse_expf/logf replaced by local avx_expf/avx_logf (two SSE calls)
 */
#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <stdio.h>
#include <string.h>
#include <math.h>

#include <immintrin.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_sse.h"       /* esl_sse_expf / esl_sse_logf for avx helpers */
#include "esl_avx.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"

/* Forward declarations of internal static helpers */
static uint8_t unbiased_byteify(P7_OPROFILE *om, float sc);
static uint8_t biased_byteify  (P7_OPROFILE *om, float sc);
static int16_t wordify         (P7_OPROFILE *om, float sc);
static int     sf_conversion   (P7_OPROFILE *om);
static int     mf_conversion   (const P7_PROFILE *gm, P7_OPROFILE *om);
static int     vf_conversion   (const P7_PROFILE *gm, P7_OPROFILE *om);
static int     fb_conversion   (const P7_PROFILE *gm, P7_OPROFILE *om);
static int     fb_conversion_log(const P7_PROFILE *gm, P7_OPROFILE *om);

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
 * 1. P7_OPROFILE: allocation, initialization, destruction.
 *****************************************************************/

/* Function:  p7_oprofile_Create_avx()
 * Synopsis:  Allocate an optimized profile structure (AVX2).
 */
P7_OPROFILE *
p7_oprofile_Create_avx(int allocM, const ESL_ALPHABET *abc)
{
  int          status;
  P7_OPROFILE *om     = NULL;
  int          nqb    = p7O_NQB_AVX(allocM);   /* # uint8  vectors: ceil(M/32) */
  int          nqw    = p7O_NQW_AVX(allocM);   /* # int16  vectors: ceil(M/16) */
  int          nqf    = p7O_NQF_AVX(allocM);   /* # float  vectors: ceil(M/8)  */
  int          nqs    = nqb + p7O_EXTRA_SB;
  int          x;

  ESL_ALLOC(om, sizeof(P7_OPROFILE));

  /* SSE fields — not used on the AVX path */
  om->rbv_mem = NULL;  om->sbv_mem = NULL;
  om->rwv_mem = NULL;  om->twv_mem = NULL;
  om->rfv_mem = NULL;  om->tfv_mem = NULL;
  om->rbv     = NULL;  om->sbv     = NULL;
  om->rwv     = NULL;  om->twv     = NULL;
  om->rfv     = NULL;  om->tfv     = NULL;
  om->allocQ16 = 0;    om->allocQ8 = 0;    om->allocQ4 = 0;

  /* AVX fields */
  om->rbv_mem_avx = NULL;  om->sbv_mem_avx = NULL;
  om->rwv_mem_avx = NULL;  om->twv_mem_avx = NULL;
  om->rfv_mem_avx = NULL;  om->tfv_mem_avx = NULL;
  om->rbv_avx     = NULL;  om->sbv_avx     = NULL;
  om->rwv_avx     = NULL;  om->twv_avx     = NULL;
  om->rfv_avx     = NULL;  om->tfv_avx     = NULL;

#ifdef eslENABLE_AVX512
  om->rbv_mem_avx512 = NULL;  om->sbv_mem_avx512 = NULL;
  om->rwv_mem_avx512 = NULL;  om->twv_mem_avx512 = NULL;
  om->rfv_mem_avx512 = NULL;  om->tfv_mem_avx512 = NULL;
  om->rbv_avx512 = NULL;  om->sbv_avx512 = NULL;
  om->rwv_avx512 = NULL;  om->twv_avx512 = NULL;
  om->rfv_avx512 = NULL;  om->tfv_avx512 = NULL;
  om->allocQ4_avx512 = 0;  om->allocQ8_avx512 = 0;  om->allocQ16_avx512 = 0;
#endif
  om->clone = 0;

  /* Allocate AVX vector memory: +31 for 32-byte alignment */
  ESL_ALLOC(om->rbv_mem_avx, sizeof(__m256i) * nqb * abc->Kp          + 31);
  ESL_ALLOC(om->sbv_mem_avx, sizeof(__m256i) * nqs * abc->Kp          + 31);
  ESL_ALLOC(om->rwv_mem_avx, sizeof(__m256i) * nqw * abc->Kp          + 31);
  ESL_ALLOC(om->twv_mem_avx, sizeof(__m256i) * nqw * p7O_NTRANS       + 31);
  ESL_ALLOC(om->rfv_mem_avx, sizeof(__m256)  * nqf * abc->Kp          + 31);
  ESL_ALLOC(om->tfv_mem_avx, sizeof(__m256)  * nqf * p7O_NTRANS       + 31);

  /* Allocate row-pointer arrays */
  ESL_ALLOC(om->rbv_avx, sizeof(__m256i *) * abc->Kp);
  ESL_ALLOC(om->sbv_avx, sizeof(__m256i *) * abc->Kp);
  ESL_ALLOC(om->rwv_avx, sizeof(__m256i *) * abc->Kp);
  ESL_ALLOC(om->rfv_avx, sizeof(__m256  *) * abc->Kp);

  /* Align to 32-byte boundaries */
  om->rbv_avx[0] = (__m256i *) (((unsigned long int) om->rbv_mem_avx + 31) & (~0x1f));
  om->sbv_avx[0] = (__m256i *) (((unsigned long int) om->sbv_mem_avx + 31) & (~0x1f));
  om->rwv_avx[0] = (__m256i *) (((unsigned long int) om->rwv_mem_avx + 31) & (~0x1f));
  om->twv_avx    = (__m256i *) (((unsigned long int) om->twv_mem_avx + 31) & (~0x1f));
  om->rfv_avx[0] = (__m256  *) (((unsigned long int) om->rfv_mem_avx + 31) & (~0x1f));
  om->tfv_avx    = (__m256  *) (((unsigned long int) om->tfv_mem_avx + 31) & (~0x1f));

  for (x = 1; x < abc->Kp; x++) {
    om->rbv_avx[x] = om->rbv_avx[0] + (x * nqb);
    om->sbv_avx[x] = om->sbv_avx[0] + (x * nqs);
    om->rwv_avx[x] = om->rwv_avx[0] + (x * nqw);
    om->rfv_avx[x] = om->rfv_avx[0] + (x * nqf);
  }
  om->allocQ16_avx = nqb;
  om->allocQ8_avx  = nqw;
  om->allocQ4_avx  = nqf;

  om->tbm_b = 0;  om->tec_b = 0;  om->tjb_b = 0;
  om->scale_b = 0.0f;  om->base_b = 0;  om->bias_b = 0;
  om->scale_w = 0.0f;  om->base_w = 0;  om->ddbound_w = 0;
  om->ncj_roundoff = 0.0f;

  for (x = 0; x < p7_NOFFSETS; x++) om->offs[x]    = -1;
  for (x = 0; x < p7_NEVPARAM; x++) om->evparam[x] = p7_EVPARAM_UNSET;
  for (x = 0; x < p7_NCUTOFFS; x++) om->cutoff[x]  = p7_CUTOFF_UNSET;
  for (x = 0; x < p7_MAXABET;  x++) om->compo[x]   = p7_COMPO_UNSET;

  om->name = NULL;  om->acc = NULL;  om->desc = NULL;

  ESL_ALLOC(om->rf,        sizeof(char) * (allocM + 2));
  ESL_ALLOC(om->mm,        sizeof(char) * (allocM + 2));
  ESL_ALLOC(om->cs,        sizeof(char) * (allocM + 2));
  ESL_ALLOC(om->consensus, sizeof(char) * (allocM + 2));
  memset(om->rf,        '\0', sizeof(char) * (allocM + 2));
  memset(om->mm,        '\0', sizeof(char) * (allocM + 2));
  memset(om->cs,        '\0', sizeof(char) * (allocM + 2));
  memset(om->consensus, '\0', sizeof(char) * (allocM + 2));

  om->abc        = abc;
  om->L          = 0;
  om->M          = 0;
  om->max_length = -1;
  om->allocM     = allocM;
  om->mode       = p7_NO_MODE;
  om->nj         = 0.0f;
  return om;

 ERROR:
  p7_oprofile_Destroy_avx(om);
  return NULL;
}


/* Function:  p7_oprofile_Destroy_avx()
 * Synopsis:  Free an optimized profile (AVX2 fields).
 */
void
p7_oprofile_Destroy_avx(P7_OPROFILE *om)
{
  if (om == NULL) return;
  if (om->clone == 0)
    {
      if (om->rbv_mem_avx != NULL) free(om->rbv_mem_avx);
      if (om->sbv_mem_avx != NULL) free(om->sbv_mem_avx);
      if (om->rwv_mem_avx != NULL) free(om->rwv_mem_avx);
      if (om->twv_mem_avx != NULL) free(om->twv_mem_avx);
      if (om->rfv_mem_avx != NULL) free(om->rfv_mem_avx);
      if (om->tfv_mem_avx != NULL) free(om->tfv_mem_avx);
      if (om->rbv_avx     != NULL) free(om->rbv_avx);
      if (om->sbv_avx     != NULL) free(om->sbv_avx);
      if (om->rwv_avx     != NULL) free(om->rwv_avx);
      if (om->rfv_avx     != NULL) free(om->rfv_avx);
      if (om->name        != NULL) free(om->name);
      if (om->acc         != NULL) free(om->acc);
      if (om->desc        != NULL) free(om->desc);
      if (om->rf          != NULL) free(om->rf);
      if (om->mm          != NULL) free(om->mm);
      if (om->cs          != NULL) free(om->cs);
      if (om->consensus   != NULL) free(om->consensus);
    }
  free(om);
}


/* Function:  p7_oprofile_Clone_avx()
 * Synopsis:  Shallow clone (borrow all pointers) for multithreading.
 */
P7_OPROFILE *
p7_oprofile_Clone_avx(const P7_OPROFILE *om1)
{
  int          status;
  P7_OPROFILE *om2 = NULL;

  ESL_ALLOC(om2, sizeof(P7_OPROFILE));
  memcpy(om2, om1, sizeof(P7_OPROFILE));
  om2->clone = 1;
  return om2;

 ERROR:
  p7_oprofile_Destroy_avx(om2);
  return NULL;
}


/* Function:  p7_oprofile_UpdateFwdEmissionScores_avx()
 */
int
p7_oprofile_UpdateFwdEmissionScores_avx(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr)
{
  int   M   = om->M;
  int   k, q, x, z;
  int   nq  = p7O_NQF_AVX(M);
  int   K   = om->abc->K;
  int   Kp  = om->abc->Kp;
  union { __m256 v; float x[8]; } tmp;

  for (k = 1, q = 0; q < nq; q++, k++) {
    for (x = 0; x < K; x++) {
      for (z = 0; z < 8; z++) {
        if (k + z*nq <= M) sc_arr[z*Kp + x] = (om->mm && om->mm[(k+z*nq)]=='m') ? 0 :
                                                log((double)(fwd_emissions[Kp*(k+z*nq)+x])/bg->f[x]);
        else               sc_arr[z*Kp + x] = -eslINFINITY;
        tmp.x[z] = sc_arr[z*Kp + x];
      }
      om->rfv_avx[x][q] = avx_expf(tmp.v);
    }
    for (z = 0; z < 8; z++) {
      sc_arr[z*Kp + K]        = -eslINFINITY;
      sc_arr[z*Kp + (Kp - 2)] = -eslINFINITY;
      sc_arr[z*Kp + (Kp - 1)] = -eslINFINITY;
    }
    for (z = 0; z < 8; z++)
      esl_abc_FExpectScVec(om->abc, sc_arr + (z*Kp), bg->f);
    for (x = K; x < Kp; x++) {
      for (z = 0; z < 8; z++) tmp.x[z] = sc_arr[z*Kp + x];
      om->rfv_avx[x][q] = avx_expf(tmp.v);
    }
  }
  return eslOK;
}


/* Function:  p7_oprofile_UpdateVitEmissionScores_avx()
 */
int
p7_oprofile_UpdateVitEmissionScores_avx(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr)
{
  int   M   = om->M;
  int   k, q, x, z, idx;
  int   nq  = p7O_NQW_AVX(M);
  int   K   = om->abc->K;
  int   Kp  = om->abc->Kp;
  union { __m256i v; int16_t i[16]; } tmp;

  for (k = 1, q = 0; q < nq; q++, k++) {
    for (x = 0; x < K; x++) {
      for (z = 0; z < 16; z++) {
        idx = z*Kp + x;
        if (k + z*nq <= M) {
          sc_arr[idx] = (om->mm && om->mm[(k+z*nq)]=='m') ? 0 :
                         log((double)(fwd_emissions[Kp*(k+z*nq)+x])/bg->f[x]);
          tmp.i[z] = wordify(om, sc_arr[idx]);
        } else {
          sc_arr[idx] = -eslINFINITY;
          tmp.i[z]    = -32768;
        }
      }
      om->rwv_avx[x][q] = tmp.v;
    }
    for (z = 0; z < 16; z++)
      esl_abc_FExpectScVec(om->abc, sc_arr + (z*Kp), bg->f);
    for (x = K; x < Kp; x++) {
      for (z = 0; z < 16; z++) {
        idx = z*Kp + x;
        if (x==K || x>Kp-3 || sc_arr[idx] == -eslINFINITY) tmp.i[z] = -32768;
        else                                                  tmp.i[z] = wordify(om, sc_arr[idx]);
      }
      om->rwv_avx[x][q] = tmp.v;
    }
  }
  return eslOK;
}


/* Function:  p7_oprofile_UpdateMSVEmissionScores_avx()
 */
int
p7_oprofile_UpdateMSVEmissionScores_avx(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr)
{
  int   M   = om->M;
  int   k, q, x, z, idx;
  int   nq  = p7O_NQB_AVX(M);
  int   K   = om->abc->K;
  int   Kp  = om->abc->Kp;
  float max = 0.0f;
  union { __m256i v; uint8_t i[32]; } tmp;

  for (x = 0; x < K; x++) max = ESL_MAX(max, esl_vec_FMax(fwd_emissions + Kp, M*Kp));
  om->scale_b = 3.0f / eslCONST_LOG2;
  om->base_b  = 190;
  om->bias_b  = unbiased_byteify(om, -1.0f * max);

  for (k = 1, q = 0; q < nq; q++, k++) {
    for (x = 0; x < K; x++) {
      for (z = 0; z < 32; z++) {
        idx = z*Kp + x;
        if (k + z*nq <= M) {
          sc_arr[idx] = (om->mm && om->mm[(k+z*nq)]=='m') ? 0 :
                         log((double)(fwd_emissions[Kp*(k+z*nq)+x])/bg->f[x]);
          tmp.i[z] = biased_byteify(om, sc_arr[idx]);
        } else {
          sc_arr[idx] = -eslINFINITY;
          tmp.i[z]    = 255;
        }
      }
      om->rbv_avx[x][q] = tmp.v;
    }
    for (z = 0; z < 32; z++) {
      sc_arr[z*Kp + K]        = -eslINFINITY;
      sc_arr[z*Kp + (Kp - 2)] = -eslINFINITY;
      sc_arr[z*Kp + (Kp - 1)] = -eslINFINITY;
    }
    for (z = 0; z < 32; z++)
      esl_abc_FExpectScVec(om->abc, sc_arr + (z*Kp), bg->f);
    for (x = K; x < Kp; x++) {
      for (z = 0; z < 32; z++) {
        idx = z*Kp + x;
        tmp.i[z] = (x==K || x>Kp-3 || sc_arr[idx]==-eslINFINITY) ? 255 : biased_byteify(om, sc_arr[idx]);
      }
      om->rbv_avx[x][q] = tmp.v;
    }
  }
  sf_conversion(om);
  return eslOK;
}


/*****************************************************************
 * 2. Conversion from P7_PROFILE to P7_OPROFILE (AVX2).
 *****************************************************************/

static uint8_t
unbiased_byteify(P7_OPROFILE *om, float sc)
{
  uint8_t b;
  sc = -1.0f * roundf(om->scale_b * sc);
  b  = (sc > 255.) ? 255 : (uint8_t) sc;
  return b;
}

static uint8_t
biased_byteify(P7_OPROFILE *om, float sc)
{
  uint8_t b;
  sc = -1.0f * roundf(om->scale_b * sc);
  b  = (sc > 255 - om->bias_b) ? 255 : (uint8_t) sc + om->bias_b;
  return b;
}

static int16_t
wordify(P7_OPROFILE *om, float sc)
{
  float rv = roundf(om->scale_w * sc);
  if      (rv >=  32767.0f) return  32767;
  else if (rv <= -32768.0f) return -32768;
  else return (int16_t) rv;
}

/* sf_conversion(): build the SSV scoring table (sbv_avx) from rbv_avx. */
static int
sf_conversion(P7_OPROFILE *om)
{
  int     M   = om->M;
  int     nq  = p7O_NQB_AVX(M);
  int     x, q;
  __m256i tmp  = _mm256_set1_epi8((int8_t)(om->bias_b + 127));
  __m256i tmp2 = _mm256_set1_epi8(127);

  for (x = 0; x < om->abc->Kp; x++) {
    for (q = 0;  q < nq;               q++) om->sbv_avx[x][q] = _mm256_xor_si256(_mm256_subs_epu8(tmp, om->rbv_avx[x][q]), tmp2);
    for (q = nq; q < nq+p7O_EXTRA_SB;  q++) om->sbv_avx[x][q] = om->sbv_avx[x][q % nq];
  }
  return eslOK;
}

/* mf_conversion(): fill MSV byte scoring table (rbv_avx) from generic profile. */
static int
mf_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int     M   = gm->M;
  int     nq  = p7O_NQB_AVX(M);
  float   max = 0.0f;
  int     x, q, k, z;
  union { __m256i v; uint8_t i[32]; } tmp;

  if (nq > om->allocQ16_avx) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  for (x = 0; x < gm->abc->K; x++) max = ESL_MAX(max, esl_vec_FMax(gm->rsc[x], (M+1)*2));
  om->scale_b = 3.0f / eslCONST_LOG2;
  om->base_b  = 190;
  om->bias_b  = unbiased_byteify(om, -1.0f * max);

  for (x = 0; x < gm->abc->Kp; x++)
    for (q = 0, k = 1; q < nq; q++, k++) {
      for (z = 0; z < 32; z++) tmp.i[z] = ((k+z*nq <= M) ? biased_byteify(om, p7P_MSC(gm, k+z*nq, x)) : 255);
      om->rbv_avx[x][q] = tmp.v;
    }

  om->tbm_b = unbiased_byteify(om, logf(2.0f / ((float)gm->M * (float)(gm->M+1))));
  om->tec_b = unbiased_byteify(om, logf(0.5f));
  om->tjb_b = unbiased_byteify(om, logf(3.0f / (float)(gm->L+3)));

  sf_conversion(om);
  return eslOK;
}

/* vf_conversion(): fill Viterbi filter int16 scoring tables (rwv_avx, twv_avx). */
static int
vf_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int     M   = gm->M;
  int     nq  = p7O_NQW_AVX(M);
  int     x, q, k, kb, z, t, tg, j, ddtmp;
  int16_t maxval, val;
  union { __m256i v; int16_t i[16]; } tmp;

  if (nq > om->allocQ8_avx) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  om->scale_w = 500.0f / eslCONST_LOG2;
  om->base_w  = 12000;

  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 1, q = 0; q < nq; q++, k++) {
      for (z = 0; z < 16; z++) tmp.i[z] = ((k+z*nq <= M) ? wordify(om, p7P_MSC(gm, k+z*nq, x)) : -32768);
      om->rwv_avx[x][q] = tmp.v;
    }

  for (j = 0, k = 1, q = 0; q < nq; q++, k++) {
    for (t = p7O_BM; t <= p7O_II; t++) {
      switch (t) {
      case p7O_BM: tg = p7P_BM; kb = k-1; maxval =  0; break;
      case p7O_MM: tg = p7P_MM; kb = k-1; maxval =  0; break;
      case p7O_IM: tg = p7P_IM; kb = k-1; maxval =  0; break;
      case p7O_DM: tg = p7P_DM; kb = k-1; maxval =  0; break;
      case p7O_MD: tg = p7P_MD; kb = k;   maxval =  0; break;
      case p7O_MI: tg = p7P_MI; kb = k;   maxval =  0; break;
      case p7O_II: tg = p7P_II; kb = k;   maxval = -1; break;
      default: tg = 0; kb = 0; maxval = 0; break;
      }
      for (z = 0; z < 16; z++) {
        val      = ((kb+z*nq < M) ? wordify(om, p7P_TSC(gm, kb+z*nq, tg)) : -32768);
        tmp.i[z] = (val <= maxval) ? val : maxval;
      }
      om->twv_avx[j++] = tmp.v;
    }
  }

  for (k = 1, q = 0; q < nq; q++, k++) {
    for (z = 0; z < 16; z++) tmp.i[z] = ((k+z*nq < M) ? wordify(om, p7P_TSC(gm, k+z*nq, p7P_DD)) : -32768);
    om->twv_avx[j++] = tmp.v;
  }

  om->xw[p7O_E][p7O_LOOP] = wordify(om, gm->xsc[p7P_E][p7P_LOOP]);
  om->xw[p7O_E][p7O_MOVE] = wordify(om, gm->xsc[p7P_E][p7P_MOVE]);
  om->xw[p7O_N][p7O_MOVE] = wordify(om, gm->xsc[p7P_N][p7P_MOVE]);
  om->xw[p7O_N][p7O_LOOP] = 0;
  om->xw[p7O_C][p7O_MOVE] = wordify(om, gm->xsc[p7P_C][p7P_MOVE]);
  om->xw[p7O_C][p7O_LOOP] = 0;
  om->xw[p7O_J][p7O_MOVE] = wordify(om, gm->xsc[p7P_J][p7P_MOVE]);
  om->xw[p7O_J][p7O_LOOP] = 0;
  om->ncj_roundoff = 0.0f;

  om->ddbound_w = -32768;
  for (k = 2; k < M-1; k++) {
    ddtmp          = (int) wordify(om, p7P_TSC(gm, k,   p7P_DD));
    ddtmp         += (int) wordify(om, p7P_TSC(gm, k+1, p7P_DM));
    ddtmp         -= (int) wordify(om, p7P_TSC(gm, k+1, p7P_BM));
    om->ddbound_w  = ESL_MAX(om->ddbound_w, ddtmp);
  }
  return eslOK;
}

/* fb_conversion(): fill Forward/Backward float scoring tables (rfv_avx, tfv_avx). */
static int
fb_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int   M  = gm->M;
  int   nq = p7O_NQF_AVX(M);
  int   x, q, k, kb, z, t, tg, j;
  union { __m256 v; float x[8]; } tmp;

  if (nq > om->allocQ4_avx) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 1, q = 0; q < nq; q++, k++) {
      for (z = 0; z < 8; z++) tmp.x[z] = (k+z*nq <= M) ? p7P_MSC(gm, k+z*nq, x) : -eslINFINITY;
      om->rfv_avx[x][q] = avx_expf(tmp.v);
    }

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
      default: tg = 0; kb = 0; break;
      }
      for (z = 0; z < 8; z++) tmp.x[z] = (kb+z*nq < M) ? p7P_TSC(gm, kb+z*nq, tg) : -eslINFINITY;
      om->tfv_avx[j++] = avx_expf(tmp.v);
    }
  }
  for (k = 1, q = 0; q < nq; q++, k++) {
    for (z = 0; z < 8; z++) tmp.x[z] = (k+z*nq < M) ? p7P_TSC(gm, k+z*nq, p7P_DD) : -eslINFINITY;
    om->tfv_avx[j++] = avx_expf(tmp.v);
  }

  om->xf[p7O_E][p7O_LOOP] = expf(gm->xsc[p7P_E][p7P_LOOP]);
  om->xf[p7O_E][p7O_MOVE] = expf(gm->xsc[p7P_E][p7P_MOVE]);
  om->xf[p7O_N][p7O_LOOP] = expf(gm->xsc[p7P_N][p7P_LOOP]);
  om->xf[p7O_N][p7O_MOVE] = expf(gm->xsc[p7P_N][p7P_MOVE]);
  om->xf[p7O_C][p7O_LOOP] = expf(gm->xsc[p7P_C][p7P_LOOP]);
  om->xf[p7O_C][p7O_MOVE] = expf(gm->xsc[p7P_C][p7P_MOVE]);
  om->xf[p7O_J][p7O_LOOP] = expf(gm->xsc[p7P_J][p7P_LOOP]);
  om->xf[p7O_J][p7O_MOVE] = expf(gm->xsc[p7P_J][p7P_MOVE]);
  return eslOK;
}

/* fb_conversion_log(): like fb_conversion() but stores log probs directly. */
static int
fb_conversion_log(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int   M  = gm->M;
  int   nq = p7O_NQF_AVX(M);
  int   x, q, k, kb, z, t, tg, j;
  union { __m256 v; float x[8]; } tmp;

  if (nq > om->allocQ4_avx) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 1, q = 0; q < nq; q++, k++) {
      for (z = 0; z < 8; z++) tmp.x[z] = (k+z*nq <= M) ? p7P_MSC(gm, k+z*nq, x) : -eslINFINITY;
      om->rfv_avx[x][q] = tmp.v;
    }

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
      default: tg = 0; kb = 0; break;
      }
      for (z = 0; z < 8; z++) tmp.x[z] = (kb+z*nq < M) ? p7P_TSC(gm, kb+z*nq, tg) : -eslINFINITY;
      om->tfv_avx[j++] = tmp.v;
    }
  }
  for (k = 1, q = 0; q < nq; q++, k++) {
    for (z = 0; z < 8; z++) tmp.x[z] = (k+z*nq < M) ? p7P_TSC(gm, k+z*nq, p7P_DD) : -eslINFINITY;
    om->tfv_avx[j++] = tmp.v;
  }

  om->xf[p7O_E][p7O_LOOP] = gm->xsc[p7P_E][p7P_LOOP];
  om->xf[p7O_E][p7O_MOVE] = gm->xsc[p7P_E][p7P_MOVE];
  om->xf[p7O_N][p7O_LOOP] = gm->xsc[p7P_N][p7P_LOOP];
  om->xf[p7O_N][p7O_MOVE] = gm->xsc[p7P_N][p7P_MOVE];
  om->xf[p7O_C][p7O_LOOP] = gm->xsc[p7P_C][p7P_LOOP];
  om->xf[p7O_C][p7O_MOVE] = gm->xsc[p7P_C][p7P_MOVE];
  om->xf[p7O_J][p7O_LOOP] = gm->xsc[p7P_J][p7P_LOOP];
  om->xf[p7O_J][p7O_MOVE] = gm->xsc[p7P_J][p7P_MOVE];
  return eslOK;
}


/* Function:  p7_oprofile_Convert_avx()
 */
int
p7_oprofile_Convert_avx(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int status, z;

  om->mode       = gm->mode;
  om->L          = gm->L;
  om->M          = gm->M;
  om->nj         = gm->nj;
  om->max_length = gm->max_length;

  if (gm->abc->type != om->abc->type) ESL_EXCEPTION(eslEINVAL, "alphabets of the two profiles don't match");
  if (gm->M         >  om->allocM)    ESL_EXCEPTION(eslEINVAL, "oprofile is too small");

  if ((status = mf_conversion    (gm, om)) != eslOK) return status;
  if ((status = vf_conversion    (gm, om)) != eslOK) return status;
  if ((status = fb_conversion    (gm, om)) != eslOK) return status;

  if (om->name != NULL) free(om->name);
  if (om->acc  != NULL) free(om->acc);
  if (om->desc != NULL) free(om->desc);
  if ((status = esl_strdup(gm->name, -1, &om->name)) != eslOK) goto ERROR;
  if ((status = esl_strdup(gm->acc,  -1, &om->acc))  != eslOK) goto ERROR;
  if ((status = esl_strdup(gm->desc, -1, &om->desc)) != eslOK) goto ERROR;
  strcpy(om->rf,        gm->rf);
  strcpy(om->mm,        gm->mm);
  strcpy(om->cs,        gm->cs);
  strcpy(om->consensus, gm->consensus);
  for (z = 0; z < p7_NEVPARAM; z++) om->evparam[z] = gm->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) om->cutoff[z]  = gm->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) om->compo[z]   = gm->compo[z];
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_oprofile_Convert_Log_avx()
 */
int
p7_oprofile_Convert_Log_avx(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int status, z;

  om->mode       = gm->mode;
  om->L          = gm->L;
  om->M          = gm->M;
  om->nj         = gm->nj;
  om->max_length = gm->max_length;

  if (gm->abc->type != om->abc->type) ESL_EXCEPTION(eslEINVAL, "alphabets of the two profiles don't match");
  if (gm->M         >  om->allocM)    ESL_EXCEPTION(eslEINVAL, "oprofile is too small");

  if ((status = mf_conversion     (gm, om)) != eslOK) return status;
  if ((status = vf_conversion     (gm, om)) != eslOK) return status;
  if ((status = fb_conversion_log (gm, om)) != eslOK) return status;

  if (om->name != NULL) free(om->name);
  if (om->acc  != NULL) free(om->acc);
  if (om->desc != NULL) free(om->desc);
  if ((status = esl_strdup(gm->name, -1, &om->name)) != eslOK) goto ERROR;
  if ((status = esl_strdup(gm->acc,  -1, &om->acc))  != eslOK) goto ERROR;
  if ((status = esl_strdup(gm->desc, -1, &om->desc)) != eslOK) goto ERROR;
  strcpy(om->rf,        gm->rf);
  strcpy(om->mm,        gm->mm);
  strcpy(om->cs,        gm->cs);
  strcpy(om->consensus, gm->consensus);
  for (z = 0; z < p7_NEVPARAM; z++) om->evparam[z] = gm->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) om->cutoff[z]  = gm->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) om->compo[z]   = gm->compo[z];
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_oprofile_Logify_avx()
 * Synopsis:  Convert exp-space rfv_avx/tfv_avx tables to log space in place.
 */
int
p7_oprofile_Logify_avx(P7_OPROFILE *om)
{
  int Q = p7O_NQF_AVX(om->M);
  int j, x, s, t;

  for (x = 0; x < om->abc->Kp; x++)
    for (j = 0; j < Q; j++)
      om->rfv_avx[x][j] = avx_logf(om->rfv_avx[x][j]);

  for (j = 0; j < 8 * Q; j++)
    om->tfv_avx[j] = avx_logf(om->tfv_avx[j]);

  for (s = 0; s < p7O_NXSTATES; s++)
    for (t = 0; t < p7O_NXTRANS; t++)
      om->xf[s][t] = logf(om->xf[s][t]);

  return eslOK;
}


/*****************************************************************
 * 3. ReconfigLength and friends (AVX2).
 *****************************************************************/

int
p7_oprofile_ReconfigLength_avx(P7_OPROFILE *om, int L)
{
  int status;
  if ((status = p7_oprofile_ReconfigMSVLength_avx (om, L)) != eslOK) return status;
  if ((status = p7_oprofile_ReconfigRestLength_avx(om, L)) != eslOK) return status;
  return eslOK;
}

int
p7_oprofile_ReconfigMSVLength_avx(P7_OPROFILE *om, int L)
{
  om->tjb_b = unbiased_byteify(om, logf(3.0f / (float)(L + 3)));
  return eslOK;
}

int
p7_oprofile_ReconfigRestLength_avx(P7_OPROFILE *om, int L)
{
  float pmove = (2.0f + om->nj) / ((float)L + 2.0f + om->nj);
  float ploop = 1.0f - pmove;

  om->xf[p7O_N][p7O_LOOP] = om->xf[p7O_C][p7O_LOOP] = om->xf[p7O_J][p7O_LOOP] = ploop;
  om->xf[p7O_N][p7O_MOVE] = om->xf[p7O_C][p7O_MOVE] = om->xf[p7O_J][p7O_MOVE] = pmove;

  om->xw[p7O_N][p7O_MOVE] = om->xw[p7O_C][p7O_MOVE] = om->xw[p7O_J][p7O_MOVE] = wordify(om, logf(pmove));

  om->L = L;
  return eslOK;
}

int
p7_oprofile_ReconfigLength_Log_avx(P7_OPROFILE *om, int L)
{
  float pmove = (2.0f + om->nj) / ((float)L + 2.0f + om->nj);
  float ploop = 1.0f - pmove;

  om->xf[p7O_N][p7O_LOOP] = om->xf[p7O_C][p7O_LOOP] = om->xf[p7O_J][p7O_LOOP] = logf(ploop);
  om->xf[p7O_N][p7O_MOVE] = om->xf[p7O_C][p7O_MOVE] = om->xf[p7O_J][p7O_MOVE] = logf(pmove);

  om->L = L;
  return eslOK;
}

int
p7_oprofile_ReconfigMultihit_avx(P7_OPROFILE *om, int L)
{
  om->xf[p7O_E][p7O_MOVE] = 0.5f;
  om->xf[p7O_E][p7O_LOOP] = 0.5f;
  om->nj = 1.0f;
  om->xw[p7O_E][p7O_MOVE] = wordify(om, -eslCONST_LOG2);
  om->xw[p7O_E][p7O_LOOP] = wordify(om, -eslCONST_LOG2);
  return p7_oprofile_ReconfigLength_avx(om, L);
}

int
p7_oprofile_ReconfigMultihit_Log_avx(P7_OPROFILE *om, int L)
{
  om->xf[p7O_E][p7O_MOVE] = -eslCONST_LOG2;
  om->xf[p7O_E][p7O_LOOP] = -eslCONST_LOG2;
  om->nj = 1.0f;
  return p7_oprofile_ReconfigLength_Log_avx(om, L);
}

int
p7_oprofile_ReconfigUnihit_avx(P7_OPROFILE *om, int L)
{
  om->xf[p7O_E][p7O_MOVE] = 1.0f;
  om->xf[p7O_E][p7O_LOOP] = 0.0f;
  om->nj = 0.0f;
  om->xw[p7O_E][p7O_MOVE] = 0;
  om->xw[p7O_E][p7O_LOOP] = -32768;
  return p7_oprofile_ReconfigLength_avx(om, L);
}

int
p7_oprofile_ReconfigUnihit_Log_avx(P7_OPROFILE *om, int L)
{
  om->xf[p7O_E][p7O_MOVE] = 0.0f;
  om->xf[p7O_E][p7O_LOOP] = -eslINFINITY;
  om->nj = 0.0f;
  return p7_oprofile_ReconfigLength_Log_avx(om, L);
}

/*****************************************************************
 * Score array extraction (AVX2).
 *****************************************************************/

int
p7_oprofile_GetFwdTransitionArray_avx(const P7_OPROFILE *om, int type, float *arr)
{
  int nq = p7O_NQF_AVX(om->M);
  int i, j;
  union { __m256 v; float x[8]; } tmp;

  for (i = 0; i < nq; i++) {
    tmp.v = om->tfv_avx[(type == p7O_DD ? nq*7+i : type+7*i)];
    for (j = 0; j < 8; j++)
      if (i+1+j*nq < om->M+1) arr[i+1+j*nq] = tmp.x[j];
  }
  return eslOK;
}

int
p7_oprofile_GetSSVEmissionScoreArray_avx(const P7_OPROFILE *om, uint8_t *arr)
{
  int x, q, z, k;
  union { __m256i v; uint8_t i[32]; } tmp;
  int M        = om->M;
  int K        = om->abc->Kp;
  int nq       = p7O_NQB_AVX(M);
  int cell_cnt = (om->M + 1) * K;

  for (x = 0; x < K; x++)
    for (q = 0, k = 1; q < nq; q++, k++) {
      tmp.v = om->rbv_avx[x][q];
      for (z = 0; z < 32; z++)
        if ((K*(k+z*nq)+x) < cell_cnt) arr[K*(k+z*nq)+x] = tmp.i[z];
    }
  return eslOK;
}

int
p7_oprofile_GetFwdEmissionScoreArray_avx(const P7_OPROFILE *om, float *arr)
{
  int x, q, z, k;
  union { __m256 v; float f[8]; } tmp;
  int M        = om->M;
  int K        = om->abc->Kp;
  int nq       = p7O_NQF_AVX(M);
  int cell_cnt = (om->M + 1) * K;

  for (x = 0; x < K; x++)
    for (q = 0, k = 1; q < nq; q++, k++) {
      tmp.v = avx_logf(om->rfv_avx[x][q]);
      for (z = 0; z < 8; z++)
        if ((K*(k+z*nq)+x) < cell_cnt) arr[K*(k+z*nq)+x] = tmp.f[z];
    }
  return eslOK;
}

int
p7_oprofile_GetFwdEmissionArray_avx(const P7_OPROFILE *om, P7_BG *bg, float *arr)
{
  int x, q, z, k;
  union { __m256 v; float f[8]; } tmp;
  int M        = om->M;
  int Kp       = om->abc->Kp;
  int K        = om->abc->K;
  int nq       = p7O_NQF_AVX(M);
  int cell_cnt = (om->M + 1) * Kp;

  for (x = 0; x < K; x++)
    for (q = 0, k = 1; q < nq; q++, k++) {
      tmp.v = om->rfv_avx[x][q];
      for (z = 0; z < 8; z++)
        if ((Kp*(k+z*nq)+x) < cell_cnt) arr[Kp*(k+z*nq)+x] = tmp.f[z] * bg->f[x];
    }

  for (x = 0; x <= M; x++)
    esl_abc_FExpectScVec(om->abc, arr + Kp*x, bg->f);
  return eslOK;
}

/* Function:  p7_oprofile_FGetEmission_avx()
 * Synopsis:  Retrieve match odds ratio [k][x] from AVX-allocated profile.
 */
float
p7_oprofile_FGetEmission_avx(const P7_OPROFILE *om, int k, int x)
{
  union { __m256 v; float p[8]; } u;
  int Q = p7O_NQF_AVX(om->M);
  u.v = om->rfv_avx[x][(k-1) % Q];
  return u.p[(k-1)/Q];
}


/* Function:  p7_fs_oprofile_FGetEmission_avx()
 * Synopsis:  Retrieve float match emission score [k][c] from AVX-allocated FS profile.
 */
float
p7_fs_oprofile_FGetEmission_avx(const P7_FS_OPROFILE *om_fs, int k, int c)
{
  union { __m256 v; float p[8]; } u;
  int Q = p7O_NQF_AVX(om_fs->M);
  u.v = om_fs->rfv_avx[c][(k-1) % Q];
  return u.p[(k-1)/Q];
}

#endif /* eslENABLE_AVX */
