/* P7_OMX SSE implementations.
 * Optimised DP matrix lifecycle and debug routines for the SSE (128-bit) path.
 * Ported from impl_sse/p7_omx.c with _sse suffix for runtime dispatch.
 *
 * Contents:
 *    1. P7_OMX SSE lifecycle (Create / GrowTo / Destroy).
 *    2. Debug / dump routines (SSE-specific unpack).
 */

#include "p7_config.h"

#include <stdio.h>
#include <math.h>
#include <float.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"

/*****************************************************************
 * 1. P7_OMX SSE lifecycle.
 *****************************************************************/

/* Function:  p7_omx_Create_sse()
 *
 * Purpose:   Allocate a <P7_OMX> for models up to size <allocM>
 *            and sequences up to length <allocL>/<allocXL> using
 *            SSE 128-bit vectors.  AVX/AVX-512 row-pointer fields
 *            are set to NULL and their allocQ fields to 0.
 *
 *            Pass <allocL=0> and <allocXL=0> for a one-row scorer.
 *            Pass <allocL=0>, <allocXL=L> for a parser (X states only).
 *            Pass <allocL=allocXL=L> for a full traceback matrix.
 *
 * Returns:   pointer to new <P7_OMX> on success.
 * Throws:    NULL on allocation failure.
 */
P7_OMX *
p7_omx_Create_sse(int allocM, int allocL, int allocXL)
{
#ifdef eslENABLE_SSE
  P7_OMX *ox    = NULL;
  int     i;
  int     status;

  ESL_ALLOC(ox, sizeof(P7_OMX));
  ox->dp_mem = NULL;
  ox->dpb    = NULL;
  ox->dpw    = NULL;
  ox->dpf    = NULL;
  ox->xmx    = NULL;
  ox->x_mem  = NULL;

#ifdef eslENABLE_AVX
  ox->dpf_avx = NULL;
  ox->dpw_avx = NULL;
  ox->dpb_avx = NULL;
  ox->allocQ4_avx  = 0;
  ox->allocQ8_avx  = 0;
  ox->allocQ16_avx = 0;
#endif
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
  ox->allocQ4  = p7O_NQF(allocM);
  ox->allocQ8  = p7O_NQW(allocM);
  ox->allocQ16 = p7O_NQB(allocM);
  ox->ncells   = (size_t) ox->allocR * ox->allocQ4 * 4;

  ESL_ALLOC(ox->dp_mem, sizeof(__m128) * ox->allocR * ox->allocQ4 * p7X_NSCELLS + 15);
  ESL_ALLOC(ox->dpb,    sizeof(__m128i *) * ox->allocR);
  ESL_ALLOC(ox->dpw,    sizeof(__m128i *) * ox->allocR);
  ESL_ALLOC(ox->dpf,    sizeof(__m128  *) * ox->allocR);

  ox->dpb[0] = (__m128i *) (((unsigned long int) ((char *) ox->dp_mem + 15)) & (~0xf));
  ox->dpw[0] = (__m128i *) (((unsigned long int) ((char *) ox->dp_mem + 15)) & (~0xf));
  ox->dpf[0] = (__m128  *) (((unsigned long int) ((char *) ox->dp_mem + 15)) & (~0xf));

  for (i = 1; i <= allocL; i++) {
    ox->dpf[i] = ox->dpf[0] + i * ox->allocQ4  * p7X_NSCELLS;
    ox->dpw[i] = ox->dpw[0] + i * ox->allocQ8  * p7X_NSCELLS;
    ox->dpb[i] = ox->dpb[0] + i * ox->allocQ16;
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
  p7_omx_Destroy_sse(ox);
  return NULL;
#else
  return NULL;
#endif
}


/* Function:  p7_omx_GrowTo_sse()
 *
 * Purpose:   Assure that SSE-allocated <P7_OMX> <ox> can hold a model
 *            of up to <allocM> nodes and a sequence of up to <allocL>
 *            residues.  Re-allocates if necessary.  Any data previously
 *            in <ox> must be considered invalid after a reallocation.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_omx_GrowTo_sse(P7_OMX *ox, int allocM, int allocL, int allocXL)
{
#ifdef eslENABLE_SSE
  void  *p;
  int    nqf    = p7O_NQF(allocM);
  int    nqw    = p7O_NQW(allocM);
  int    nqb    = p7O_NQB(allocM);
  size_t ncells = (size_t)(allocL + 1) * nqf * 4;
  int    reset_row_pointers = FALSE;
  int    i;
  int    status;

  if (ox->allocQ4 * 4 >= allocM && ox->validR > allocL && ox->allocXR >= allocXL + 1)
    return eslOK;

  if (ncells > ox->ncells) {
    ESL_RALLOC(ox->dp_mem, p, sizeof(__m128) * (allocL + 1) * nqf * p7X_NSCELLS + 15);
    ox->ncells = ncells;
    reset_row_pointers = TRUE;
  }

  if (allocXL + 1 >= ox->allocXR) {
    ESL_RALLOC(ox->x_mem, p, sizeof(float) * (allocXL + 1) * p7X_NXCELLS + 15);
    ox->allocXR = allocXL + 1;
    ox->xmx     = (float *) (((unsigned long int) ((char *) ox->x_mem + 15)) & (~0xf));
  }

  if (allocL >= ox->allocR) {
    ESL_RALLOC(ox->dpb, p, sizeof(__m128i *) * (allocL + 1));
    ESL_RALLOC(ox->dpw, p, sizeof(__m128i *) * (allocL + 1));
    ESL_RALLOC(ox->dpf, p, sizeof(__m128  *) * (allocL + 1));
    ox->allocR         = allocL + 1;
    reset_row_pointers = TRUE;
  }

  if (allocM > ox->allocQ4 * 4) reset_row_pointers = TRUE;
  if (allocL >= ox->validR)      reset_row_pointers = TRUE;

  if (reset_row_pointers) {
    ox->dpb[0] = (__m128i *) (((unsigned long int) ((char *) ox->dp_mem + 15)) & (~0xf));
    ox->dpw[0] = (__m128i *) (((unsigned long int) ((char *) ox->dp_mem + 15)) & (~0xf));
    ox->dpf[0] = (__m128  *) (((unsigned long int) ((char *) ox->dp_mem + 15)) & (~0xf));

    ox->validR = ESL_MIN((int)(ox->ncells / (nqf * 4)), ox->allocR);
    for (i = 1; i < ox->validR; i++) {
      ox->dpb[i] = ox->dpb[0] + i * nqb;
      ox->dpw[i] = ox->dpw[0] + i * nqw * p7X_NSCELLS;
      ox->dpf[i] = ox->dpf[0] + i * nqf * p7X_NSCELLS;
    }

    ox->allocQ4  = nqf;
    ox->allocQ8  = nqw;
    ox->allocQ16 = nqb;
  }

  ox->M = 0;
  ox->L = 0;
  return eslOK;

 ERROR:
  return status;
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_omx_Create_dpf_sse()
 *
 * Purpose:   Like p7_omx_Create_sse(), but allocates only float (dpf) rows
 *            with <nscells> vectors per stripe position.  dpw and dpb are
 *            not allocated.  Pass p7X_NSCELLS (3) for standard layout or
 *            p7X_NSCELLS_FS (8) for frameshift layout.
 *
 * Returns:   pointer to new <P7_OMX> on success, NULL on failure.
 */
P7_OMX *
p7_omx_Create_dpf_sse(int allocM, int allocL, int allocXL, int nscells)
{
#ifdef eslENABLE_SSE
  P7_OMX *ox    = NULL;
  int     i;
  int     status;

  ESL_ALLOC(ox, sizeof(P7_OMX));
  ox->dp_mem = NULL;
  ox->dpb    = NULL;
  ox->dpw    = NULL;
  ox->dpf    = NULL;
  ox->xmx    = NULL;
  ox->x_mem  = NULL;

#ifdef eslENABLE_AVX
  ox->dpf_avx = NULL;
  ox->dpw_avx = NULL;
  ox->dpb_avx = NULL;
  ox->allocQ4_avx  = 0;
  ox->allocQ8_avx  = 0;
  ox->allocQ16_avx = 0;
#endif
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
  ox->allocQ4  = p7O_NQF(allocM);
  ox->allocQ8  = 0;
  ox->allocQ16 = 0;
  ox->ncells   = (size_t) ox->allocR * ox->allocQ4 * 4;

  ESL_ALLOC(ox->dp_mem, sizeof(__m128) * ox->allocR * ox->allocQ4 * nscells + 15);
  ESL_ALLOC(ox->dpf,    sizeof(__m128 *) * ox->allocR);

  ox->dpf[0] = (__m128 *) (((unsigned long int) ((char *) ox->dp_mem + 15)) & (~0xf));
  for (i = 1; i <= allocL; i++)
    ox->dpf[i] = ox->dpf[0] + i * ox->allocQ4 * nscells;

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
  p7_omx_Destroy_sse(ox);
  return NULL;
#else
  return NULL;
#endif
}


/* Function:  p7_omx_GrowTo_dpf_sse()
 *
 * Purpose:   Like p7_omx_GrowTo_sse(), but for dpf-only matrices.
 *            Uses ox->nscells as the per-stripe cell count.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on failure.
 */
int
p7_omx_GrowTo_dpf_sse(P7_OMX *ox, int allocM, int allocL, int allocXL)
{
#ifdef eslENABLE_SSE
  void  *p;
  int    nqf    = p7O_NQF(allocM);
  size_t ncells = (size_t)(allocL + 1) * nqf * 4;
  int    reset_row_pointers = FALSE;
  int    i;
  int    status;

  if (ox->allocQ4 * 4 >= allocM && ox->validR > allocL && ox->allocXR >= allocXL + 1)
    return eslOK;

  if (ncells > ox->ncells) {
    ESL_RALLOC(ox->dp_mem, p, sizeof(__m128) * (allocL + 1) * nqf * ox->nscells + 15);
    ox->ncells = ncells;
    reset_row_pointers = TRUE;
  }

  if (allocXL + 1 >= ox->allocXR) {
    ESL_RALLOC(ox->x_mem, p, sizeof(float) * (allocXL + 1) * p7X_NXCELLS + 15);
    ox->allocXR = allocXL + 1;
    ox->xmx     = (float *) (((unsigned long int) ((char *) ox->x_mem + 15)) & (~0xf));
  }

  if (allocL >= ox->allocR) {
    ESL_RALLOC(ox->dpf, p, sizeof(__m128 *) * (allocL + 1));
    ox->allocR         = allocL + 1;
    reset_row_pointers = TRUE;
  }

  if (allocM > ox->allocQ4 * 4) reset_row_pointers = TRUE;
  if (allocL >= ox->validR)      reset_row_pointers = TRUE;

  if (reset_row_pointers) {
    ox->dpf[0] = (__m128 *) (((unsigned long int) ((char *) ox->dp_mem + 15)) & (~0xf));
    ox->validR  = ESL_MIN((int)(ox->ncells / (nqf * 4)), ox->allocR);
    for (i = 1; i < ox->validR; i++)
      ox->dpf[i] = ox->dpf[0] + i * nqf * ox->nscells;
    ox->allocQ4 = nqf;
  }

  ox->M = 0;
  ox->L = 0;
  return eslOK;

 ERROR:
  return status;
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_omx_Destroy_sse()
 *
 * Purpose:   Free a <P7_OMX> allocated by p7_omx_Create_sse() or
 *            p7_omx_Create_dpf_sse().  Safe to call on NULL.
 */
void
p7_omx_Destroy_sse(P7_OMX *ox)
{
  if (ox == NULL) return;

  if (ox->x_mem  != NULL) free(ox->x_mem);
  if (ox->dp_mem != NULL) free(ox->dp_mem);
#ifdef eslENABLE_SSE
  if (ox->dpf    != NULL) free(ox->dpf);
  if (ox->dpw    != NULL) free(ox->dpw);
  if (ox->dpb    != NULL) free(ox->dpb);
#endif
  free(ox);
}


/*****************************************************************
 * 2. Debug / dump routines (SSE-specific vector unpack).
 *****************************************************************/

/* Function:  p7_omx_FDeconvert_sse()
 *
 * Purpose:   Convert 32-bit float values in SSE optimized DP matrix <ox>
 *            to generic matrix <gx>.  Caller provides <gx> with sufficient
 *            space for the <ox->M> by <ox->L> matrix.
 *
 *            If ox->nscells == p7X_NSCELLS_FS, the full 8-cell frameshift
 *            layout is extracted; otherwise the standard 3-cell path is used.
 */
int
p7_omx_FDeconvert_sse(P7_OMX *ox, P7_GMX *gx)
{
#ifdef eslENABLE_SSE
  int Q = p7O_NQF(ox->M);
  int i, q, r, k, c;
  union { __m128 v; float p[4]; } u;
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
              u.v = DMO_FS(ox->dpf[i],q);
              for (r = 0; r < 4; r++) { k = Q*r+q+1; if (k <= ox->M) DMX_FS(i,k)   = u.p[r]; }
              u.v = IMO_FS(ox->dpf[i],q);
              for (r = 0; r < 4; r++) { k = Q*r+q+1; if (k <= ox->M) IMX_FS(i,k)   = u.p[r]; }
              for (c = 0; c < 6; c++) {
                u.v = MMO_FS(ox->dpf[i],q,c);
                for (r = 0; r < 4; r++) { k = Q*r+q+1; if (k <= ox->M) MMX_FS(i,k,c) = u.p[r]; }
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
              u.v = MMO(ox->dpf[i],q);
              for (r = 0; r < 4; r++) { k = Q*r+q+1; if (k <= ox->M) MMX(i,k) = u.p[r]; }
              u.v = DMO(ox->dpf[i],q);
              for (r = 0; r < 4; r++) { k = Q*r+q+1; if (k <= ox->M) DMX(i,k) = u.p[r]; }
              u.v = IMO(ox->dpf[i],q);
              for (r = 0; r < 4; r++) { k = Q*r+q+1; if (k <= ox->M) IMX(i,k) = u.p[r]; }
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
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_omx_FGetMDI_sse()
 * Synopsis:  Read one float DP cell from an SSE-allocated P7_OMX.
 */
float
p7_omx_FGetMDI_sse(const P7_OMX *ox, int s, int i, int k)
{
#ifdef eslENABLE_SSE
  union { __m128 v; float p[4]; } u;
  int Q = p7O_NQF(ox->M);
  u.v = ox->dpf[i][p7X_NSCELLS * ((k-1) % Q) + s];
  return u.p[(k-1)/Q];
#else
  return 0.0f;
#endif
}


/* Function:  p7_omx_FSetMDI_sse()
 * Synopsis:  Write one float DP cell in an SSE-allocated P7_OMX.
 */
void
p7_omx_FSetMDI_sse(const P7_OMX *ox, int s, int i, int k, float val)
{
#ifdef eslENABLE_SSE
  union { __m128 v; float p[4]; } u;
  int Q = p7O_NQF(ox->M);
  int q = p7X_NSCELLS * ((k-1) % Q) + s;
  u.v           = ox->dpf[i][q];
  u.p[(k-1)/Q]  = val;
  ox->dpf[i][q] = u.v;
#endif
}
