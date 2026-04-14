/* p7_oprofile_sse.c — SSE (128-bit) implementations of P7_OPROFILE management.
 *
 * All public functions carry the _sse suffix and are selected at runtime
 * by impl_Init() in p7_oprofile.c when SSE is the fastest available ISA.
 *
 * This file is a direct port of impl_sse/p7_oprofile.c with:
 *   - Public function names suffixed with _sse
 *   - Internal cross-calls updated to use _sse variants
 *   - Include changed from impl_sse.h to impl_avx.h
 *   - Benchmark/unit-test/example drivers omitted (live in p7_oprofile.c)
 */
#include "p7_config.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#include <xmmintrin.h>
#include <emmintrin.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_sse.h"
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


/*****************************************************************
 * 1. P7_OPROFILE: allocation, initialization, destruction.
 *****************************************************************/

/* Function:  p7_oprofile_Create_sse()
 * Synopsis:  Allocate an optimized profile structure (SSE).
 */
P7_OPROFILE *
p7_oprofile_Create_sse(int allocM, const ESL_ALPHABET *abc)
{
  int          status;
  P7_OPROFILE *om  = NULL;
  int          nqb = p7O_NQB(allocM);
  int          nqw = p7O_NQW(allocM);
  int          nqf = p7O_NQF(allocM);
  int          nqs = nqb + p7O_EXTRA_SB;
  int          x;

  ESL_ALLOC(om, sizeof(P7_OPROFILE));
  om->rbv_mem = NULL;  om->sbv_mem = NULL;
  om->rwv_mem = NULL;  om->twv_mem = NULL;
  om->rfv_mem = NULL;  om->tfv_mem = NULL;
  om->rbv     = NULL;  om->sbv     = NULL;
  om->rwv     = NULL;  om->twv     = NULL;
  om->rfv     = NULL;  om->tfv     = NULL;
#ifdef eslENABLE_AVX
  om->rbv_mem_avx = NULL;  om->sbv_mem_avx = NULL;
  om->rwv_mem_avx = NULL;  om->twv_mem_avx = NULL;
  om->rfv_mem_avx = NULL;  om->tfv_mem_avx = NULL;
  om->rbv_avx = NULL;  om->sbv_avx = NULL;
  om->rwv_avx = NULL;  om->twv_avx = NULL;
  om->rfv_avx = NULL;  om->tfv_avx = NULL;
  om->allocQ4_avx = 0;  om->allocQ8_avx = 0;  om->allocQ16_avx = 0;
#endif
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

  ESL_ALLOC(om->rbv_mem, sizeof(__m128i) * nqb * abc->Kp          + 15);
  ESL_ALLOC(om->sbv_mem, sizeof(__m128i) * nqs * abc->Kp          + 15);
  ESL_ALLOC(om->rwv_mem, sizeof(__m128i) * nqw * abc->Kp          + 15);
  ESL_ALLOC(om->twv_mem, sizeof(__m128i) * nqw * p7O_NTRANS       + 15);
  ESL_ALLOC(om->rfv_mem, sizeof(__m128)  * nqf * abc->Kp          + 15);
  ESL_ALLOC(om->tfv_mem, sizeof(__m128)  * nqf * p7O_NTRANS       + 15);

  ESL_ALLOC(om->rbv, sizeof(__m128i *) * abc->Kp);
  ESL_ALLOC(om->sbv, sizeof(__m128i *) * abc->Kp);
  ESL_ALLOC(om->rwv, sizeof(__m128i *) * abc->Kp);
  ESL_ALLOC(om->rfv, sizeof(__m128  *) * abc->Kp);

  om->rbv[0] = (__m128i *) (((unsigned long int) om->rbv_mem + 15) & (~0xf));
  om->sbv[0] = (__m128i *) (((unsigned long int) om->sbv_mem + 15) & (~0xf));
  om->rwv[0] = (__m128i *) (((unsigned long int) om->rwv_mem + 15) & (~0xf));
  om->twv    = (__m128i *) (((unsigned long int) om->twv_mem + 15) & (~0xf));
  om->rfv[0] = (__m128  *) (((unsigned long int) om->rfv_mem + 15) & (~0xf));
  om->tfv    = (__m128  *) (((unsigned long int) om->tfv_mem + 15) & (~0xf));

  for (x = 1; x < abc->Kp; x++) {
    om->rbv[x] = om->rbv[0] + (x * nqb);
    om->sbv[x] = om->sbv[0] + (x * nqs);
    om->rwv[x] = om->rwv[0] + (x * nqw);
    om->rfv[x] = om->rfv[0] + (x * nqf);
  }
  om->allocQ16 = nqb;
  om->allocQ8  = nqw;
  om->allocQ4  = nqf;

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
  p7_oprofile_Destroy_sse(om);
  return NULL;
}


/* Function:  p7_oprofile_Destroy_sse()
 * Synopsis:  Free an optimized profile (SSE fields only).
 */
void
p7_oprofile_Destroy_sse(P7_OPROFILE *om)
{
  if (om == NULL) return;
  if (om->clone == 0)
    {
      if (om->rbv_mem   != NULL) free(om->rbv_mem);
      if (om->sbv_mem   != NULL) free(om->sbv_mem);
      if (om->rwv_mem   != NULL) free(om->rwv_mem);
      if (om->twv_mem   != NULL) free(om->twv_mem);
      if (om->rfv_mem   != NULL) free(om->rfv_mem);
      if (om->tfv_mem   != NULL) free(om->tfv_mem);
      if (om->rbv       != NULL) free(om->rbv);
      if (om->sbv       != NULL) free(om->sbv);
      if (om->rwv       != NULL) free(om->rwv);
      if (om->rfv       != NULL) free(om->rfv);
      if (om->name      != NULL) free(om->name);
      if (om->acc       != NULL) free(om->acc);
      if (om->desc      != NULL) free(om->desc);
      if (om->rf        != NULL) free(om->rf);
      if (om->mm        != NULL) free(om->mm);
      if (om->cs        != NULL) free(om->cs);
      if (om->consensus != NULL) free(om->consensus);
    }
  free(om);
}


/* Function:  p7_oprofile_Clone_sse()
 * Synopsis:  Shallow clone (borrow all pointers) for multithreading.
 */
P7_OPROFILE *
p7_oprofile_Clone_sse(const P7_OPROFILE *om1)
{
  int          status;
  P7_OPROFILE *om2 = NULL;

  ESL_ALLOC(om2, sizeof(P7_OPROFILE));
  memcpy(om2, om1, sizeof(P7_OPROFILE));
  om2->clone = 1;
  return om2;

 ERROR:
  p7_oprofile_Destroy_sse(om2);
  return NULL;
}


/* Function:  p7_oprofile_UpdateFwdEmissionScores_sse()
 */
int
p7_oprofile_UpdateFwdEmissionScores_sse(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr)
{
  int   M   = om->M;
  int   k, q, x, z;
  int   nq  = p7O_NQF(M);
  int   K   = om->abc->K;
  int   Kp  = om->abc->Kp;
  union { __m128 v; float x[4]; } tmp;

  for (k = 1, q = 0; q < nq; q++, k++) {
    for (x = 0; x < K; x++) {
      for (z = 0; z < 4; z++) {
        if (k + z*nq <= M) sc_arr[z*Kp + x] = (om->mm && om->mm[(k+z*nq)]=='m') ? 0 :
                                                log((double)(fwd_emissions[Kp*(k+z*nq)+x])/bg->f[x]);
        else               sc_arr[z*Kp + x] = -eslINFINITY;
        tmp.x[z] = sc_arr[z*Kp + x];
      }
      om->rfv[x][q] = esl_sse_expf(tmp.v);
    }
    for (z = 0; z < 4; z++) {
      sc_arr[z*Kp + K]        = -eslINFINITY;
      sc_arr[z*Kp + (Kp - 2)] = -eslINFINITY;
      sc_arr[z*Kp + (Kp - 1)] = -eslINFINITY;
    }
    for (z = 0; z < 4; z++)
      esl_abc_FExpectScVec(om->abc, sc_arr + (z*Kp), bg->f);
    for (x = K; x < Kp; x++) {
      for (z = 0; z < 4; z++) tmp.x[z] = sc_arr[z*Kp + x];
      om->rfv[x][q] = esl_sse_expf(tmp.v);
    }
  }
  return eslOK;
}


/* Function:  p7_oprofile_UpdateVitEmissionScores_sse()
 */
int
p7_oprofile_UpdateVitEmissionScores_sse(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr)
{
  int   M   = om->M;
  int   k, q, x, z, idx;
  int   nq  = p7O_NQW(M);
  int   K   = om->abc->K;
  int   Kp  = om->abc->Kp;
  union { __m128i v; int16_t i[8]; } tmp;

  for (k = 1, q = 0; q < nq; q++, k++) {
    for (x = 0; x < K; x++) {
      for (z = 0; z < 8; z++) {
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
      om->rwv[x][q] = tmp.v;
    }
    for (z = 0; z < 8; z++)
      esl_abc_FExpectScVec(om->abc, sc_arr + (z*Kp), bg->f);
    for (x = K; x < Kp; x++) {
      for (z = 0; z < 8; z++) {
        idx = z*Kp + x;
        if (x==K || x>Kp-3 || sc_arr[idx] == -eslINFINITY) tmp.i[z] = -32768;
        else                                                  tmp.i[z] = wordify(om, sc_arr[idx]);
      }
      om->rwv[x][q] = tmp.v;
    }
  }
  return eslOK;
}


/* Function:  p7_oprofile_UpdateMSVEmissionScores_sse()
 */
int
p7_oprofile_UpdateMSVEmissionScores_sse(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr)
{
  int   M   = om->M;
  int   k, q, x, z, idx;
  int   nq  = p7O_NQB(M);
  int   K   = om->abc->K;
  int   Kp  = om->abc->Kp;
  float max = 0.0f;
  union { __m128i v; uint8_t i[16]; } tmp;

  for (x = 0; x < K; x++) max = ESL_MAX(max, esl_vec_FMax(fwd_emissions + Kp, M*Kp));
  om->scale_b = 3.0f / eslCONST_LOG2;
  om->base_b  = 190;
  om->bias_b  = unbiased_byteify(om, -1.0f * max);

  for (k = 1, q = 0; q < nq; q++, k++) {
    for (x = 0; x < K; x++) {
      for (z = 0; z < 16; z++) {
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
      om->rbv[x][q] = tmp.v;
    }
    for (z = 0; z < 16; z++) {
      sc_arr[z*Kp + K]        = -eslINFINITY;
      sc_arr[z*Kp + (Kp - 2)] = -eslINFINITY;
      sc_arr[z*Kp + (Kp - 1)] = -eslINFINITY;
    }
    for (z = 0; z < 16; z++)
      esl_abc_FExpectScVec(om->abc, sc_arr + (z*Kp), bg->f);
    for (x = K; x < Kp; x++) {
      for (z = 0; z < 16; z++) {
        idx = z*Kp + x;
        tmp.i[z] = (x==K || x>Kp-3 || sc_arr[idx]==-eslINFINITY) ? 255 : biased_byteify(om, sc_arr[idx]);
      }
      om->rbv[x][q] = tmp.v;
    }
  }
  sf_conversion(om);
  return eslOK;
}


/*****************************************************************
 * 2. Conversion from P7_PROFILE to P7_OPROFILE (SSE).
 *****************************************************************/

static uint8_t
unbiased_byteify(P7_OPROFILE *om, float sc)
{
  uint8_t b;
  sc = -1.0f * roundf(om->scale_b * sc);        /* convert to integer cost in a float */
  b  = (sc > 255.) ? 255 : (uint8_t) sc;        /* cast, saturate to unsigned char cost */
  return b;
}

static uint8_t
biased_byteify(P7_OPROFILE *om, float sc)
{
  uint8_t b;
  sc = -1.0f * roundf(om->scale_b * sc);        /* convert to integer cost in a float */
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

static int
sf_conversion(P7_OPROFILE *om)
{
  int     M   = om->M;
  int     nq  = p7O_NQB(M);
  int     x, q;
  __m128i tmp  = _mm_set1_epi8((int8_t)(om->bias_b + 127));
  __m128i tmp2 = _mm_set1_epi8(127);

  for (x = 0; x < om->abc->Kp; x++) {
    for (q = 0;  q < nq;              q++) om->sbv[x][q] = _mm_xor_si128(_mm_subs_epu8(tmp, om->rbv[x][q]), tmp2);
    for (q = nq; q < nq+p7O_EXTRA_SB; q++) om->sbv[x][q] = om->sbv[x][q % nq];
  }
  return eslOK;
}

static int
mf_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int     M   = gm->M;
  int     nq  = p7O_NQB(M);
  float   max = 0.0f;
  int     x, q, k, z;
  union { __m128i v; uint8_t i[16]; } tmp;

  if (nq > om->allocQ16) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  for (x = 0; x < gm->abc->K; x++) max = ESL_MAX(max, esl_vec_FMax(gm->rsc[x], (M+1)*2));
  om->scale_b = 3.0f / eslCONST_LOG2;
  om->base_b  = 190;
  om->bias_b  = unbiased_byteify(om, -1.0f * max);

  for (x = 0; x < gm->abc->Kp; x++)
    for (q = 0, k = 1; q < nq; q++, k++) {
      for (z = 0; z < 16; z++) tmp.i[z] = ((k+z*nq <= M) ? biased_byteify(om, p7P_MSC(gm, k+z*nq, x)) : 255);
      om->rbv[x][q] = tmp.v;
    }

  om->tbm_b = unbiased_byteify(om, logf(2.0f / ((float)gm->M * (float)(gm->M+1))));
  om->tec_b = unbiased_byteify(om, logf(0.5f));
  om->tjb_b = unbiased_byteify(om, logf(3.0f / (float)(gm->L+3)));

  sf_conversion(om);
  return eslOK;
}

static int
vf_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int     M   = gm->M;
  int     nq  = p7O_NQW(M);
  int     x, q, k, kb, z, t, tg, j, ddtmp;
  int16_t maxval, val;
  union { __m128i v; int16_t i[8]; } tmp;

  if (nq > om->allocQ8) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  om->scale_w = 500.0f / eslCONST_LOG2;
  om->base_w  = 12000;

  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 1, q = 0; q < nq; q++, k++) {
      for (z = 0; z < 8; z++) tmp.i[z] = ((k+z*nq <= M) ? wordify(om, p7P_MSC(gm, k+z*nq, x)) : -32768);
      om->rwv[x][q] = tmp.v;
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
      for (z = 0; z < 8; z++) {
        val      = ((kb+z*nq < M) ? wordify(om, p7P_TSC(gm, kb+z*nq, tg)) : -32768);
        tmp.i[z] = (val <= maxval) ? val : maxval;
      }
      om->twv[j++] = tmp.v;
    }
  }

  for (k = 1, q = 0; q < nq; q++, k++) {
    for (z = 0; z < 8; z++) tmp.i[z] = ((k+z*nq < M) ? wordify(om, p7P_TSC(gm, k+z*nq, p7P_DD)) : -32768);
    om->twv[j++] = tmp.v;
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

static int
fb_conversion(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int   M  = gm->M;
  int   nq = p7O_NQF(M);
  int   x, q, k, kb, z, t, tg, j;
  union { __m128 v; float x[4]; } tmp;

  if (nq > om->allocQ4) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 1, q = 0; q < nq; q++, k++) {
      for (z = 0; z < 4; z++) tmp.x[z] = (k+z*nq <= M) ? p7P_MSC(gm, k+z*nq, x) : -eslINFINITY;
      om->rfv[x][q] = esl_sse_expf(tmp.v);
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
      for (z = 0; z < 4; z++) tmp.x[z] = (kb+z*nq < M) ? p7P_TSC(gm, kb+z*nq, tg) : -eslINFINITY;
      om->tfv[j++] = esl_sse_expf(tmp.v);
    }
  }
  for (k = 1, q = 0; q < nq; q++, k++) {
    for (z = 0; z < 4; z++) tmp.x[z] = (k+z*nq < M) ? p7P_TSC(gm, k+z*nq, p7P_DD) : -eslINFINITY;
    om->tfv[j++] = esl_sse_expf(tmp.v);
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

static int
fb_conversion_log(const P7_PROFILE *gm, P7_OPROFILE *om)
{
  int   M  = gm->M;
  int   nq = p7O_NQF(M);
  int   x, q, k, kb, z, t, tg, j;
  union { __m128 v; float x[4]; } tmp;

  if (nq > om->allocQ4) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 1, q = 0; q < nq; q++, k++) {
      for (z = 0; z < 4; z++) tmp.x[z] = (k+z*nq <= M) ? p7P_MSC(gm, k+z*nq, x) : -eslINFINITY;
      om->rfv[x][q] = tmp.v;
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
      for (z = 0; z < 4; z++) tmp.x[z] = (kb+z*nq < M) ? p7P_TSC(gm, kb+z*nq, tg) : -eslINFINITY;
      om->tfv[j++] = tmp.v;
    }
  }
  for (k = 1, q = 0; q < nq; q++, k++) {
    for (z = 0; z < 4; z++) tmp.x[z] = (k+z*nq < M) ? p7P_TSC(gm, k+z*nq, p7P_DD) : -eslINFINITY;
    om->tfv[j++] = tmp.v;
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


/* Function:  p7_oprofile_Convert_sse()
 */
int
p7_oprofile_Convert_sse(const P7_PROFILE *gm, P7_OPROFILE *om)
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


/* Function:  p7_oprofile_Convert_Log_sse()
 */
int
p7_oprofile_Convert_Log_sse(const P7_PROFILE *gm, P7_OPROFILE *om)
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


/* Function:  p7_oprofile_Logify_sse()
 */
int
p7_oprofile_Logify_sse(P7_OPROFILE *om)
{
  int Q = p7O_NQF(om->M);
  int j, x, s, t;

  for (x = 0; x < om->abc->Kp; x++)
    for (j = 0; j < Q; j++)
      om->rfv[x][j] = esl_sse_logf(om->rfv[x][j]);

  for (j = 0; j < 8 * Q; j++)
    om->tfv[j] = esl_sse_logf(om->tfv[j]);

  for (s = 0; s < p7O_NXSTATES; s++)
    for (t = 0; t < p7O_NXTRANS; t++)
      om->xf[s][t] = logf(om->xf[s][t]);

  return eslOK;
}


/*****************************************************************
 * 3. ReconfigLength and friends (SSE).
 *****************************************************************/

int
p7_oprofile_ReconfigLength_sse(P7_OPROFILE *om, int L)
{
  int status;
  if ((status = p7_oprofile_ReconfigMSVLength_sse (om, L)) != eslOK) return status;
  if ((status = p7_oprofile_ReconfigRestLength_sse(om, L)) != eslOK) return status;
  return eslOK;
}

int
p7_oprofile_ReconfigMSVLength_sse(P7_OPROFILE *om, int L)
{
  om->tjb_b = unbiased_byteify(om, logf(3.0f / (float)(L + 3)));
  return eslOK;
}

int
p7_oprofile_ReconfigRestLength_sse(P7_OPROFILE *om, int L)
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
p7_oprofile_ReconfigLength_Log_sse(P7_OPROFILE *om, int L)
{
  float pmove = (2.0f + om->nj) / ((float)L + 2.0f + om->nj);
  float ploop = 1.0f - pmove;

  om->xf[p7O_N][p7O_LOOP] = om->xf[p7O_C][p7O_LOOP] = om->xf[p7O_J][p7O_LOOP] = logf(ploop);
  om->xf[p7O_N][p7O_MOVE] = om->xf[p7O_C][p7O_MOVE] = om->xf[p7O_J][p7O_MOVE] = logf(pmove);

  om->L = L;
  return eslOK;
}

int
p7_oprofile_ReconfigMultihit_sse(P7_OPROFILE *om, int L)
{
  om->xf[p7O_E][p7O_MOVE] = 0.5f;
  om->xf[p7O_E][p7O_LOOP] = 0.5f;
  om->nj = 1.0f;
  om->xw[p7O_E][p7O_MOVE] = wordify(om, -eslCONST_LOG2);
  om->xw[p7O_E][p7O_LOOP] = wordify(om, -eslCONST_LOG2);
  return p7_oprofile_ReconfigLength_sse(om, L);
}

int
p7_oprofile_ReconfigMultihit_Log_sse(P7_OPROFILE *om, int L)
{
  om->xf[p7O_E][p7O_MOVE] = -eslCONST_LOG2;
  om->xf[p7O_E][p7O_LOOP] = -eslCONST_LOG2;
  om->nj = 1.0f;
  return p7_oprofile_ReconfigLength_Log_sse(om, L);
}

int
p7_oprofile_ReconfigUnihit_sse(P7_OPROFILE *om, int L)
{
  om->xf[p7O_E][p7O_MOVE] = 1.0f;
  om->xf[p7O_E][p7O_LOOP] = 0.0f;
  om->nj = 0.0f;
  om->xw[p7O_E][p7O_MOVE] = 0;
  om->xw[p7O_E][p7O_LOOP] = -32768;
  return p7_oprofile_ReconfigLength_sse(om, L);
}

int
p7_oprofile_ReconfigUnihit_Log_sse(P7_OPROFILE *om, int L)
{
  om->xf[p7O_E][p7O_MOVE] = 0.0f;
  om->xf[p7O_E][p7O_LOOP] = -eslINFINITY;
  om->nj = 0.0f;
  return p7_oprofile_ReconfigLength_Log_sse(om, L);
}


/*****************************************************************
 * 4. Compact score array extraction (SSE).
 *****************************************************************/

int
p7_oprofile_GetFwdTransitionArray(const P7_OPROFILE *om, int type, float *arr)
{
  int nq = p7O_NQF(om->M);
  int i, j;
  union { __m128 v; float x[4]; } tmp;

  for (i = 0; i < nq; i++) {
    tmp.v = om->tfv[(type == p7O_DD ? nq*7+i : type+7*i)];
    for (j = 0; j < 4; j++)
      if (i+1+j*nq < om->M+1) arr[i+1+j*nq] = tmp.x[j];
  }
  return eslOK;
}

int
p7_oprofile_GetSSVEmissionScoreArray(const P7_OPROFILE *om, uint8_t *arr)
{
  int x, q, z, k;
  union { __m128i v; uint8_t i[16]; } tmp;
  int M        = om->M;
  int K        = om->abc->Kp;
  int nq       = p7O_NQB(M);
  int cell_cnt = (om->M + 1) * K;

  for (x = 0; x < K; x++)
    for (q = 0, k = 1; q < nq; q++, k++) {
      tmp.v = om->rbv[x][q];
      for (z = 0; z < 16; z++)
        if ((K*(k+z*nq)+x) < cell_cnt) arr[K*(k+z*nq)+x] = tmp.i[z];
    }
  return eslOK;
}

int
p7_oprofile_GetFwdEmissionScoreArray(const P7_OPROFILE *om, float *arr)
{
  int x, q, z, k;
  union { __m128 v; float f[4]; } tmp;
  int M        = om->M;
  int K        = om->abc->Kp;
  int nq       = p7O_NQF(M);
  int cell_cnt = (om->M + 1) * K;

  for (x = 0; x < K; x++)
    for (q = 0, k = 1; q < nq; q++, k++) {
      tmp.v = esl_sse_logf(om->rfv[x][q]);
      for (z = 0; z < 4; z++)
        if ((K*(k+z*nq)+x) < cell_cnt) arr[K*(k+z*nq)+x] = tmp.f[z];
    }
  return eslOK;
}

int
p7_oprofile_GetFwdEmissionArray(const P7_OPROFILE *om, P7_BG *bg, float *arr)
{
  int x, q, z, k;
  union { __m128 v; float f[4]; } tmp;
  int M        = om->M;
  int Kp       = om->abc->Kp;
  int K        = om->abc->K;
  int nq       = p7O_NQF(M);
  int cell_cnt = (om->M + 1) * Kp;

  for (x = 0; x < K; x++)
    for (q = 0, k = 1; q < nq; q++, k++) {
      tmp.v = om->rfv[x][q];
      for (z = 0; z < 4; z++)
        if ((Kp*(k+z*nq)+x) < cell_cnt) arr[Kp*(k+z*nq)+x] = tmp.f[z] * bg->f[x];
    }

  for (x = 0; x <= M; x++)
    esl_abc_FExpectScVec(om->abc, arr + Kp*x, bg->f);
  return eslOK;
}


/*****************************************************************
 * 5. Debug comparison utilities (SSE).
 *****************************************************************/

int
p7_profile_SameAsMF_sse(const P7_OPROFILE *om, P7_PROFILE *gm)
{
  int    k, x;
  float  tbm = roundf(om->scale_b * (log(2.0f / ((float) gm->M * (float) (gm->M+1)))));

  /* Transitions */
  esl_vec_FSet(gm->tsc, p7P_NTRANS * gm->M, -eslINFINITY);
  for (k = 1; k <  gm->M; k++) p7P_TSC(gm, k, p7P_MM) = 0.0f;
  for (k = 0; k <  gm->M; k++) p7P_TSC(gm, k, p7P_BM) = tbm;

  /* Emissions */
  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 0; k <= gm->M; k++)
      {
        gm->rsc[x][k*2]   = (gm->rsc[x][k*2] <= -eslINFINITY) ? -eslINFINITY : roundf(om->scale_b * gm->rsc[x][k*2]);
        gm->rsc[x][k*2+1] = 0;  /* insert score: MSV makes it zero no matter what. */
      }

  /* Specials */
  for (k = 0; k < p7P_NXSTATES; k++)
    for (x = 0; x < p7P_NXTRANS; x++)
      gm->xsc[k][x] = (gm->xsc[k][x] <= -eslINFINITY) ? -eslINFINITY : roundf(om->scale_b * gm->xsc[k][x]);

  /* NN, CC, JJ hardcoded 0 in limited precision */
  gm->xsc[p7P_N][p7P_LOOP] = gm->xsc[p7P_J][p7P_LOOP] = gm->xsc[p7P_C][p7P_LOOP] = 0;

  return eslOK;
}

int
p7_profile_SameAsVF_sse(const P7_OPROFILE *om, P7_PROFILE *gm)
{
  int k, x;

  /* Transitions */
  for (x = 0; x < gm->M*p7P_NTRANS; x++)
    gm->tsc[x] = (gm->tsc[x] <= -eslINFINITY) ? -eslINFINITY : roundf(om->scale_w * gm->tsc[x]);

  /* Enforce the rule that no II can be 0; max of -1 */
  for (x = p7P_II; x < gm->M*p7P_NTRANS; x += p7P_NTRANS)
    if (gm->tsc[x] == 0.0) gm->tsc[x] = -1.0;

  /* Emissions */
  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 0; k <= gm->M; k++)
      {
        gm->rsc[x][k*2]   = (gm->rsc[x][k*2]   <= -eslINFINITY) ? -eslINFINITY : roundf(om->scale_w * gm->rsc[x][k*2]);
        gm->rsc[x][k*2+1] = 0.0;  /* insert score: VF makes it zero no matter what. */
      }

  /* Specials */
  for (k = 0; k < p7P_NXSTATES; k++)
    for (x = 0; x < p7P_NXTRANS; x++)
      gm->xsc[k][x] = (gm->xsc[k][x] <= -eslINFINITY) ? -eslINFINITY : roundf(om->scale_w * gm->xsc[k][x]);

  /* 3nat approximation: NN, CC, JJ hardcoded 0 in limited precision */
  gm->xsc[p7P_N][p7P_LOOP] = gm->xsc[p7P_J][p7P_LOOP] = gm->xsc[p7P_C][p7P_LOOP] = 0.0;

  return eslOK;
}
