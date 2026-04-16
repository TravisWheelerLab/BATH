/* p7_omx_avx.c — AVX2 (256-bit) P7_OMX matrix lifecycle.
 *
 * All public functions carry the _avx suffix and are selected at
 * runtime by impl_Init() when AVX2 is the fastest available ISA.
 *
 * Key differences from p7_omx_sse.c:
 *   - dp_mem is sized for __m256 / __m256i (32 bytes each)
 *   - dpf_avx / dpw_avx / dpb_avx used instead of dpf / dpw / dpb
 *   - allocQ4_avx / allocQ8_avx / allocQ16_avx
 *   - p7O_NQF_AVX / p7O_NQW_AVX / p7O_NQB_AVX macros
 *   - 32-byte alignment (+31 / & ~0x1f)
 */
#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <stdio.h>
#include <math.h>
#include <float.h>

#include <immintrin.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"


/*****************************************************************
 * 1. P7_OMX AVX lifecycle.
 *****************************************************************/

/* Function:  p7_omx_Create_avx()
 *
 * Purpose:   Allocate a <P7_OMX> for models up to size <allocM>
 *            and sequences up to length <allocL>/<allocXL> using
 *            AVX2 256-bit vectors.  SSE row-pointer fields are set
 *            to NULL and their allocQ fields to 0.
 *
 *            Pass <allocL=0, allocXL=0>  for a one-row scorer.
 *            Pass <allocL=0, allocXL=L>  for a parser (X states only).
 *            Pass <allocL=allocXL=L>     for a full traceback matrix.
 *
 * Returns:   pointer to new <P7_OMX> on success.
 * Throws:    NULL on allocation failure.
 */
P7_OMX *
p7_omx_Create_avx(int allocM, int allocL, int allocXL)
{
  P7_OMX *ox    = NULL;
  int     i;
  int     status;

  ESL_ALLOC(ox, sizeof(P7_OMX));
  ox->dp_mem = NULL;
  ox->xmx    = NULL;
  ox->x_mem  = NULL;

  /* SSE fields — unused on AVX path */
#ifdef eslENABLE_SSE
  ox->dpb    = NULL;
  ox->dpw    = NULL;
  ox->dpf    = NULL;
#endif
  ox->allocQ4  = 0;
  ox->allocQ8  = 0;
  ox->allocQ16 = 0;

  /* AVX fields */
  ox->dpf_avx = NULL;
  ox->dpw_avx = NULL;
  ox->dpb_avx = NULL;
  ox->allocQ4_avx  = p7O_NQF_AVX(allocM);
  ox->allocQ8_avx  = p7O_NQW_AVX(allocM);
  ox->allocQ16_avx = p7O_NQB_AVX(allocM);

#ifdef eslENABLE_AVX512
  ox->dpf_avx512 = NULL;
  ox->dpw_avx512 = NULL;
  ox->dpb_avx512 = NULL;
  ox->allocQ4_avx512  = 0;
  ox->allocQ8_avx512  = 0;
  ox->allocQ16_avx512 = 0;
#endif

  ox->nscells  = p7X_NSCELLS;
  ox->allocR   = allocL + 1;
  ox->validR   = ox->allocR;
  ox->ncells   = (size_t) ox->allocR * ox->allocQ4_avx * 4;

  /* Single block for all DP rows; +31 for 32-byte alignment */
  ESL_ALLOC(ox->dp_mem, sizeof(__m256) * ox->allocR * ox->allocQ4_avx * p7X_NSCELLS + 31);
  ESL_ALLOC(ox->dpb_avx, sizeof(__m256i *) * ox->allocR);
  ESL_ALLOC(ox->dpw_avx, sizeof(__m256i *) * ox->allocR);
  ESL_ALLOC(ox->dpf_avx, sizeof(__m256  *) * ox->allocR);

  ox->dpb_avx[0] = (__m256i *) (((unsigned long int) ((char *) ox->dp_mem + 31)) & (~0x1f));
  ox->dpw_avx[0] = (__m256i *) (((unsigned long int) ((char *) ox->dp_mem + 31)) & (~0x1f));
  ox->dpf_avx[0] = (__m256  *) (((unsigned long int) ((char *) ox->dp_mem + 31)) & (~0x1f));

  for (i = 1; i <= allocL; i++) {
    ox->dpf_avx[i] = ox->dpf_avx[0] + i * ox->allocQ4_avx  * p7X_NSCELLS;
    ox->dpw_avx[i] = ox->dpw_avx[0] + i * ox->allocQ8_avx  * p7X_NSCELLS;
    ox->dpb_avx[i] = ox->dpb_avx[0] + i * ox->allocQ16_avx;
  }

  ox->allocXR = allocXL + 1;
  ESL_ALLOC(ox->x_mem, sizeof(float) * ox->allocXR * p7X_NXCELLS + 15);
  ox->xmx = (float *) (((unsigned long int) ((char *) ox->x_mem + 15)) & (~0xf));

  ox->M              = 0;
  ox->L              = 0;
  ox->totscale       = 0.0;
  ox->has_own_scales = TRUE;
#if eslDEBUGLEVEL > 0
  ox->debugging = FALSE;
  ox->dfp       = NULL;
#endif
  return ox;

 ERROR:
  p7_omx_Destroy_avx(ox);
  return NULL;
}


/* Function:  p7_omx_GrowTo_avx()
 *
 * Purpose:   Assure that AVX-allocated <P7_OMX> <ox> can hold a model
 *            of up to <allocM> nodes and a sequence of up to <allocL>
 *            residues.  Re-allocates if necessary.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_omx_GrowTo_avx(P7_OMX *ox, int allocM, int allocL, int allocXL)
{
  void  *p;
  int    nqf    = p7O_NQF_AVX(allocM);
  int    nqw    = p7O_NQW_AVX(allocM);
  int    nqb    = p7O_NQB_AVX(allocM);
  size_t ncells = (size_t)(allocL + 1) * nqf * 4;
  int    reset_row_pointers = FALSE;
  int    i;
  int    status;

  if (ox->allocQ4_avx * 4 >= allocM && ox->validR > allocL && ox->allocXR >= allocXL + 1)
    return eslOK;

  if (ncells > ox->ncells) {
    ESL_RALLOC(ox->dp_mem, p, sizeof(__m256) * (allocL + 1) * nqf * p7X_NSCELLS + 31);
    ox->ncells = ncells;
    reset_row_pointers = TRUE;
  }

  if (allocXL + 1 >= ox->allocXR) {
    ESL_RALLOC(ox->x_mem, p, sizeof(float) * (allocXL + 1) * p7X_NXCELLS + 15);
    ox->allocXR = allocXL + 1;
    ox->xmx     = (float *) (((unsigned long int) ((char *) ox->x_mem + 15)) & (~0xf));
  }

  if (allocL >= ox->allocR) {
    ESL_RALLOC(ox->dpb_avx, p, sizeof(__m256i *) * (allocL + 1));
    ESL_RALLOC(ox->dpw_avx, p, sizeof(__m256i *) * (allocL + 1));
    ESL_RALLOC(ox->dpf_avx, p, sizeof(__m256  *) * (allocL + 1));
    ox->allocR          = allocL + 1;
    reset_row_pointers  = TRUE;
  }

  if (allocM > ox->allocQ4_avx * 4) reset_row_pointers = TRUE;
  if (allocL >= ox->validR)          reset_row_pointers = TRUE;

  if (reset_row_pointers) {
    ox->dpb_avx[0] = (__m256i *) (((unsigned long int) ((char *) ox->dp_mem + 31)) & (~0x1f));
    ox->dpw_avx[0] = (__m256i *) (((unsigned long int) ((char *) ox->dp_mem + 31)) & (~0x1f));
    ox->dpf_avx[0] = (__m256  *) (((unsigned long int) ((char *) ox->dp_mem + 31)) & (~0x1f));

    ox->validR = ESL_MIN((int)(ox->ncells / (nqf * 4)), ox->allocR);
    for (i = 1; i < ox->validR; i++) {
      ox->dpb_avx[i] = ox->dpb_avx[0] + i * nqb;
      ox->dpw_avx[i] = ox->dpw_avx[0] + i * nqw * p7X_NSCELLS;
      ox->dpf_avx[i] = ox->dpf_avx[0] + i * nqf * p7X_NSCELLS;
    }

    ox->allocQ4_avx  = nqf;
    ox->allocQ8_avx  = nqw;
    ox->allocQ16_avx = nqb;
  }

  ox->M = 0;
  ox->L = 0;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_omx_Create_dpf_avx()
 *
 * Purpose:   Like p7_omx_Create_avx(), but allocates only float (dpf_avx)
 *            rows with <nscells> vectors per stripe position.
 *            dpw_avx and dpb_avx are not allocated.
 *
 * Returns:   pointer to new <P7_OMX> on success, NULL on failure.
 */
P7_OMX *
p7_omx_Create_dpf_avx(int allocM, int allocL, int allocXL, int nscells)
{
  P7_OMX *ox    = NULL;
  int     i;
  int     status;

  ESL_ALLOC(ox, sizeof(P7_OMX));
  ox->dp_mem = NULL;
  ox->xmx    = NULL;
  ox->x_mem  = NULL;

#ifdef eslENABLE_SSE
  ox->dpb    = NULL;
  ox->dpw    = NULL;
  ox->dpf    = NULL;
#endif
  ox->allocQ4  = 0;
  ox->allocQ8  = 0;
  ox->allocQ16 = 0;

  ox->dpf_avx = NULL;
  ox->dpw_avx = NULL;
  ox->dpb_avx = NULL;
  ox->allocQ4_avx  = p7O_NQF_AVX(allocM);
  ox->allocQ8_avx  = 0;
  ox->allocQ16_avx = 0;

#ifdef eslENABLE_AVX512
  ox->dpf_avx512 = NULL;
  ox->dpw_avx512 = NULL;
  ox->dpb_avx512 = NULL;
  ox->allocQ4_avx512  = 0;
  ox->allocQ8_avx512  = 0;
  ox->allocQ16_avx512 = 0;
#endif

  ox->nscells  = nscells;
  ox->allocR   = allocL + 1;
  ox->validR   = ox->allocR;
  ox->ncells   = (size_t) ox->allocR * ox->allocQ4_avx * 4;

  ESL_ALLOC(ox->dp_mem, sizeof(__m256) * ox->allocR * ox->allocQ4_avx * nscells + 31);
  ESL_ALLOC(ox->dpf_avx, sizeof(__m256 *) * ox->allocR);

  ox->dpf_avx[0] = (__m256 *) (((unsigned long int) ((char *) ox->dp_mem + 31)) & (~0x1f));
  for (i = 1; i <= allocL; i++)
    ox->dpf_avx[i] = ox->dpf_avx[0] + i * ox->allocQ4_avx * nscells;

  ox->allocXR = allocXL + 1;
  ESL_ALLOC(ox->x_mem, sizeof(float) * ox->allocXR * p7X_NXCELLS + 15);
  ox->xmx = (float *) (((unsigned long int) ((char *) ox->x_mem + 15)) & (~0xf));

  ox->M              = 0;
  ox->L              = 0;
  ox->totscale       = 0.0;
  ox->has_own_scales = TRUE;
#if eslDEBUGLEVEL > 0
  ox->debugging = FALSE;
  ox->dfp       = NULL;
#endif
  return ox;

 ERROR:
  p7_omx_Destroy_avx(ox);
  return NULL;
}


/* Function:  p7_omx_GrowTo_dpf_avx()
 *
 * Purpose:   Like p7_omx_GrowTo_avx(), but for dpf_avx-only matrices.
 *            Uses ox->nscells as the per-stripe cell count.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on failure.
 */
int
p7_omx_GrowTo_dpf_avx(P7_OMX *ox, int allocM, int allocL, int allocXL)
{
  void  *p;
  int    nqf    = p7O_NQF_AVX(allocM);
  size_t ncells = (size_t)(allocL + 1) * nqf * 4;
  int    reset_row_pointers = FALSE;
  int    i;
  int    status;

  if (ox->allocQ4_avx * 4 >= allocM && ox->validR > allocL && ox->allocXR >= allocXL + 1)
    return eslOK;

  if (ncells > ox->ncells) {
    ESL_RALLOC(ox->dp_mem, p, sizeof(__m256) * (allocL + 1) * nqf * ox->nscells + 31);
    ox->ncells = ncells;
    reset_row_pointers = TRUE;
  }

  if (allocXL + 1 >= ox->allocXR) {
    ESL_RALLOC(ox->x_mem, p, sizeof(float) * (allocXL + 1) * p7X_NXCELLS + 15);
    ox->allocXR = allocXL + 1;
    ox->xmx     = (float *) (((unsigned long int) ((char *) ox->x_mem + 15)) & (~0xf));
  }

  if (allocL >= ox->allocR) {
    ESL_RALLOC(ox->dpf_avx, p, sizeof(__m256 *) * (allocL + 1));
    ox->allocR          = allocL + 1;
    reset_row_pointers  = TRUE;
  }

  if (allocM > ox->allocQ4_avx * 4) reset_row_pointers = TRUE;
  if (allocL >= ox->validR)          reset_row_pointers = TRUE;

  if (reset_row_pointers) {
    ox->dpf_avx[0] = (__m256 *) (((unsigned long int) ((char *) ox->dp_mem + 31)) & (~0x1f));
    ox->validR  = ESL_MIN((int)(ox->ncells / (nqf * 4)), ox->allocR);
    for (i = 1; i < ox->validR; i++)
      ox->dpf_avx[i] = ox->dpf_avx[0] + i * nqf * ox->nscells;
    ox->allocQ4_avx = nqf;
  }

  ox->M = 0;
  ox->L = 0;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_omx_FDeconvert_avx()
 *
 * Purpose:   Convert 32-bit float values in AVX2 optimized DP matrix <ox>
 *            to generic matrix <gx>.  Caller provides <gx> with sufficient
 *            space for the <ox->M> by <ox->L> matrix.
 *
 *            If ox->nscells == p7X_NSCELLS_FS, the full 8-cell frameshift
 *            layout is extracted; otherwise the standard 3-cell path is used.
 */
int
p7_omx_FDeconvert_avx(P7_OMX *ox, P7_GMX *gx)
{
  int Q = p7O_NQF_AVX(ox->M);
  int i, q, r, k, c;
  union { __m256 v; float p[8]; } u;
  float  **dp  = gx->dp;
  float   *xmx = gx->xmx;

  if (ox->nscells == p7X_NSCELLS_FS)
    {
      for (i = 0; i <= ox->L; i++)
        {
          for (c = 0; c < 6; c++) MMX_FS(i,0,c) = -eslINFINITY;
          IMX_FS(i,0) = DMX_FS(i,0) = -eslINFINITY;
          for (q = 0; q < Q; q++)
            {
              u.v = DMO_FS(ox->dpf_avx[i],q);
              for (r = 0; r < 8; r++) { k = Q*r+q+1; if (k <= ox->M) DMX_FS(i,k)   = u.p[r]; }
              u.v = IMO_FS(ox->dpf_avx[i],q);
              for (r = 0; r < 8; r++) { k = Q*r+q+1; if (k <= ox->M) IMX_FS(i,k)   = u.p[r]; }
              for (c = 0; c < 6; c++) {
                u.v = MMO_FS(ox->dpf_avx[i],q,c);
                for (r = 0; r < 8; r++) { k = Q*r+q+1; if (k <= ox->M) MMX_FS(i,k,c) = u.p[r]; }
              }
            }
          XMX(i,p7G_E) = ox->xmx[i*p7X_NXCELLS+p7X_E];
          XMX(i,p7G_N) = ox->xmx[i*p7X_NXCELLS+p7X_N];
          XMX(i,p7G_J) = ox->xmx[i*p7X_NXCELLS+p7X_J];
          XMX(i,p7G_B) = ox->xmx[i*p7X_NXCELLS+p7X_B];
          XMX(i,p7G_C) = ox->xmx[i*p7X_NXCELLS+p7X_C];
        }
    }
  else
    {
      for (i = 0; i <= ox->L; i++)
        {
          MMX(i,0) = DMX(i,0) = IMX(i,0) = -eslINFINITY;
          for (q = 0; q < Q; q++)
            {
              u.v = MMO(ox->dpf_avx[i],q);
              for (r = 0; r < 8; r++) { k = Q*r+q+1; if (k <= ox->M) MMX(i,k) = u.p[r]; }
              u.v = DMO(ox->dpf_avx[i],q);
              for (r = 0; r < 8; r++) { k = Q*r+q+1; if (k <= ox->M) DMX(i,k) = u.p[r]; }
              u.v = IMO(ox->dpf_avx[i],q);
              for (r = 0; r < 8; r++) { k = Q*r+q+1; if (k <= ox->M) IMX(i,k) = u.p[r]; }
            }
          XMX(i,p7G_E) = ox->xmx[i*p7X_NXCELLS+p7X_E];
          XMX(i,p7G_N) = ox->xmx[i*p7X_NXCELLS+p7X_N];
          XMX(i,p7G_J) = ox->xmx[i*p7X_NXCELLS+p7X_J];
          XMX(i,p7G_B) = ox->xmx[i*p7X_NXCELLS+p7X_B];
          XMX(i,p7G_C) = ox->xmx[i*p7X_NXCELLS+p7X_C];
        }
    }

  gx->L = ox->L;
  gx->M = ox->M;
  return eslOK;
}


/* Function:  p7_omx_Destroy_avx()
 *
 * Purpose:   Free a <P7_OMX> allocated by p7_omx_Create_avx() or
 *            p7_omx_Create_dpf_avx().  Safe to call on NULL.
 */
void
p7_omx_Destroy_avx(P7_OMX *ox)
{
  if (ox == NULL) return;

  if (ox->x_mem   != NULL) free(ox->x_mem);
  if (ox->dp_mem  != NULL) free(ox->dp_mem);
  if (ox->dpf_avx != NULL) free(ox->dpf_avx);
  if (ox->dpw_avx != NULL) free(ox->dpw_avx);
  if (ox->dpb_avx != NULL) free(ox->dpb_avx);
  free(ox);
}

#endif /* eslENABLE_AVX */
