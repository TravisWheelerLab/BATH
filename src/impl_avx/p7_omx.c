/* P7_OMX dispatcher and ISA-independent functions.
 * Provides the non-suffixed extern functions declared in impl_avx.h.
 * ISA-specific work is delegated to p7_omx_sse.c (and future _avx.c/_avx512.c).
 * Debug dump routines are implemented directly here using SSE ifdefs.
 *
 * Contents:
 *    1. Dispatcher wrappers (Create / GrowTo / Destroy / FDeconvert).
 *    2. ISA-independent functions (Reuse, SetDumpMode).
 *    3. Debug dump routines (SSE vector unpack).
 */

#include "p7_config.h"

#include <stdio.h>
#include <math.h>
#include <float.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"

/*****************************************************************
 * 1. Dispatcher wrappers.
 *****************************************************************/

/* Function:  p7_omx_Create()
 *
 * Purpose:   Allocate a <P7_OMX> for models up to size <allocM> and
 *            sequences up to <allocL>/<allocXL>.  Dispatches to the
 *            fastest available ISA implementation detected at runtime.
 *
 * Returns:   pointer to new <P7_OMX> on success, NULL on failure.
 */
P7_OMX *
p7_omx_Create(int allocM, int allocL, int allocXL)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_omx_Create_avx512(allocM, allocL, allocXL);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_omx_Create_sse(allocM, allocL, allocXL);
#endif
#ifdef eslENABLE_SSE
  return p7_omx_Create_sse(allocM, allocL, allocXL);
#else
  return NULL;
#endif
}


/* Function:  p7_omx_GrowTo()
 *
 * Purpose:   Assure that <P7_OMX> <ox> can accommodate a model of up
 *            to <allocM> nodes and a sequence up to <allocL>/<allocXL>.
 *            Dispatches to the ISA-specific implementation that allocated <ox>.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_omx_GrowTo(P7_OMX *ox, int allocM, int allocL, int allocXL)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_omx_GrowTo_avx512(ox, allocM, allocL, allocXL);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_omx_GrowTo_sse(ox, allocM, allocL, allocXL);
#endif
#ifdef eslENABLE_SSE
  return p7_omx_GrowTo_sse(ox, allocM, allocL, allocXL);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_omx_Create_dpf()
 *
 * Purpose:   Like p7_omx_Create(), but allocates only float (dpf) rows.
 *            Pass p7X_NSCELLS (3) or p7X_NSCELLS_FS (8) as <nscells>.
 *
 * Returns:   pointer to new <P7_OMX> on success, NULL on failure.
 */
P7_OMX *
p7_omx_Create_dpf(int allocM, int allocL, int allocXL, int nscells)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_omx_Create_dpf_avx512(allocM, allocL, allocXL, nscells);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_omx_Create_dpf_sse(allocM, allocL, allocXL, nscells);
#endif
#ifdef eslENABLE_SSE
  return p7_omx_Create_dpf_sse(allocM, allocL, allocXL, nscells);
#else
  return NULL;
#endif
}


/* Function:  p7_omx_GrowTo_dpf()
 *
 * Purpose:   Like p7_omx_GrowTo(), but for dpf-only matrices.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_omx_GrowTo_dpf(P7_OMX *ox, int allocM, int allocL, int allocXL)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_omx_GrowTo_dpf_avx512(ox, allocM, allocL, allocXL);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_omx_GrowTo_dpf_sse(ox, allocM, allocL, allocXL);
#endif
#ifdef eslENABLE_SSE
  return p7_omx_GrowTo_dpf_sse(ox, allocM, allocL, allocXL);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_omx_FDeconvert()
 *
 * Purpose:   Convert float values in optimized DP matrix <ox> to
 *            generic matrix <gx>.  Dispatches to the ISA-specific
 *            implementation.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_omx_FDeconvert(P7_OMX *ox, P7_GMX *gx)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_omx_FDeconvert_avx512(ox, gx);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_omx_FDeconvert_sse(ox, gx);
#endif
#ifdef eslENABLE_SSE
  return p7_omx_FDeconvert_sse(ox, gx);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_omx_Destroy()
 *
 * Purpose:   Free <P7_OMX> <ox>.  Dispatches to the ISA-specific
 *            implementation.  Safe to call on NULL.
 */
void
p7_omx_Destroy(P7_OMX *ox)
{
  if (ox == NULL) return;

  /* Destroy via the same ISA that allocated ox.
   * We detect which ISA's fields are populated by checking the allocQ fields. */
#ifdef eslENABLE_AVX512
  if (ox->allocQ4_avx512 > 0 || ox->allocQ8_avx512 > 0 || ox->allocQ16_avx512 > 0) {
    p7_omx_Destroy_avx512(ox);
    return;
  }
#endif
#ifdef eslENABLE_AVX
  if (ox->allocQ4_avx > 0 || ox->allocQ8_avx > 0 || ox->allocQ16_avx > 0) {
    p7_omx_Destroy_sse(ox);
    return;
  }
#endif
#ifdef eslENABLE_SSE
  p7_omx_Destroy_sse(ox);
#else
  free(ox);
#endif
}


/*****************************************************************
 * 2. ISA-independent functions.
 *****************************************************************/

/* Function:  p7_omx_Reuse()
 *
 * Purpose:   Recycle <ox> for re-use without reallocating.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_omx_Reuse(P7_OMX *ox)
{
  ox->M              = 0;
  ox->L              = 0;
  ox->totscale       = 0.0;
  ox->has_own_scales = TRUE;
#if eslDEBUGLEVEL > 0
  ox->debugging = FALSE;
  ox->dfp       = NULL;
#endif
  return eslOK;
}


/* Function:  p7_omx_SetDumpMode()
 *
 * Purpose:   Enable or disable row-by-row DP matrix dumping for debugging.
 *            When enabled, each DP routine dumps its row immediately after
 *            computing it.  Has no effect unless compiled with eslDEBUGLEVEL > 0.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_omx_SetDumpMode(FILE *fp, P7_OMX *ox, int truefalse)
{
#if eslDEBUGLEVEL > 0
  ox->debugging = truefalse;
  ox->dfp       = fp;
#endif
  return eslOK;
}


/*****************************************************************
 * 3. Debug dump routines.
 *    These use SSE vector types for unpack; they run on any ISA
 *    path because the data is always accessible via the SSE dpf/dpw/dpb
 *    pointers on SSE-built matrices.
 *****************************************************************/

/* Function:  p7_omx_DumpMFRow()
 *
 * Purpose:   Dump current MSV (uint8) row of <ox> to its debug stream.
 *            Prints a header on row 0.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_omx_DumpMFRow(P7_OMX *ox, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC)
{
#ifdef eslENABLE_SSE
  __m128i *dp = ox->dpb[0];
  int      M  = ox->M;
  int      Q  = p7O_NQB(M);
  uint8_t *v  = NULL;
  int      q, z, k;
  union { __m128i v; uint8_t i[16]; } tmp;
  int      status;

  ESL_ALLOC(v, sizeof(uint8_t) * (Q * 16 + 1));
  v[0] = 0;

  if (rowi == 0)
    {
      fprintf(ox->dfp, "       ");
      for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", k);
      fprintf(ox->dfp, "%3s %3s %3s %3s %3s\n", "E", "N", "J", "B", "C");
      fprintf(ox->dfp, "       ");
      for (k = 0; k <= M + 5; k++) fprintf(ox->dfp, "%3s ", "---");
      fprintf(ox->dfp, "\n");
    }

  for (q = 0; q < Q; q++) {
    tmp.v = dp[q];
    for (z = 0; z < 16; z++) v[q + Q*z + 1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d M ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", v[k]);
  fprintf(ox->dfp, "%3d %3d %3d %3d %3d\n", xE, xN, xJ, xB, xC);

  fprintf(ox->dfp, "%4d I ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", 0);
  fprintf(ox->dfp, "\n");

  fprintf(ox->dfp, "%4d D ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%3d ", 0);
  fprintf(ox->dfp, "\n\n");

  free(v);
  return eslOK;

 ERROR:
  free(v);
  return status;
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_omx_DumpVFRow()
 *
 * Purpose:   Dump current ViterbiFilter (int16) row of <ox>.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_omx_DumpVFRow(P7_OMX *ox, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC)
{
#ifdef eslENABLE_SSE
  __m128i *dp = ox->dpw[0];
  int      M  = ox->M;
  int      Q  = p7O_NQW(M);
  int16_t *v  = NULL;
  int      q, z, k;
  union { __m128i v; int16_t i[8]; } tmp;
  int      status;

  ESL_ALLOC(v, sizeof(int16_t) * (Q * 8 + 1));
  v[0] = 0;

  if (rowi == 0)
    {
      fprintf(ox->dfp, "       ");
      for (k = 0; k <= M; k++) fprintf(ox->dfp, "%6d ", k);
      fprintf(ox->dfp, "%6s %6s %6s %6s %6s\n", "E", "N", "J", "B", "C");
      fprintf(ox->dfp, "       ");
      for (k = 0; k <= M + 5; k++) fprintf(ox->dfp, "%6s ", "------");
      fprintf(ox->dfp, "\n");
    }

  for (q = 0; q < Q; q++) {
    tmp.v = MMXo(q);
    for (z = 0; z < 8; z++) v[q + Q*z + 1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d M ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%6d ", v[k]);
  fprintf(ox->dfp, "%6d %6d %6d %6d %6d\n", xE, xN, xJ, xB, xC);

  for (q = 0; q < Q; q++) {
    tmp.v = IMXo(q);
    for (z = 0; z < 8; z++) v[q + Q*z + 1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d I ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%6d ", v[k]);
  fprintf(ox->dfp, "\n");

  for (q = 0; q < Q; q++) {
    tmp.v = DMXo(q);
    for (z = 0; z < 8; z++) v[q + Q*z + 1] = tmp.i[z];
  }
  fprintf(ox->dfp, "%4d D ", rowi);
  for (k = 0; k <= M; k++) fprintf(ox->dfp, "%6d ", v[k]);
  fprintf(ox->dfp, "\n\n");

  free(v);
  return eslOK;

 ERROR:
  free(v);
  return status;
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_omx_DumpFBRow()
 *
 * Purpose:   Dump current Forward/Backward (float) row of <ox>.
 *            <logify>: TRUE to print log(score), FALSE for raw.
 *            <width>, <precision>: format (8,5 for pspace; 5,2 for lspace).
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_omx_DumpFBRow(P7_OMX *ox, int logify, int rowi, int width, int precision,
                 float xE, float xN, float xJ, float xB, float xC)
{
#ifdef eslENABLE_SSE
  __m128 *dp;
  int     M  = ox->M;
  int     Q  = p7O_NQF(M);
  float  *v  = NULL;
  int     q, z, k;
  union { __m128 v; float x[4]; } tmp;
  int     status;

  dp = (ox->allocR == 1) ? ox->dpf[0] : ox->dpf[rowi];

  ESL_ALLOC(v, sizeof(float) * (Q * 4 + 1));
  v[0] = 0.0f;

  if (rowi == 0)
    {
      fprintf(ox->dfp, "      ");
      for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*d ", width, k);
      fprintf(ox->dfp, "%*s %*s %*s %*s %*s\n",
              width, "E", width, "N", width, "J", width, "B", width, "C");
      fprintf(ox->dfp, "      ");
      for (k = 0; k <= M + 5; k++) fprintf(ox->dfp, "%*s ", width, "--------");
      fprintf(ox->dfp, "\n");
    }

  for (q = 0; q < Q; q++) {
    tmp.v = MMXo(q);
    for (z = 0; z < 4; z++) v[q + Q*z + 1] = tmp.x[z];
  }
  fprintf(ox->dfp, "%3d M ", rowi);
  if (logify) for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k] == 0.0f ? -eslINFINITY : log(v[k]));
  else        for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k]);
  if (logify)
    fprintf(ox->dfp, "%*.*f %*.*f %*.*f %*.*f %*.*f\n",
            width, precision, xE == 0.0f ? -eslINFINITY : log(xE),
            width, precision, xN == 0.0f ? -eslINFINITY : log(xN),
            width, precision, xJ == 0.0f ? -eslINFINITY : log(xJ),
            width, precision, xB == 0.0f ? -eslINFINITY : log(xB),
            width, precision, xC == 0.0f ? -eslINFINITY : log(xC));
  else
    fprintf(ox->dfp, "%*.*f %*.*f %*.*f %*.*f %*.*f\n",
            width, precision, xE, width, precision, xN,
            width, precision, xJ, width, precision, xB,
            width, precision, xC);

  for (q = 0; q < Q; q++) {
    tmp.v = IMXo(q);
    for (z = 0; z < 4; z++) v[q + Q*z + 1] = tmp.x[z];
  }
  fprintf(ox->dfp, "%3d I ", rowi);
  if (logify) for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k] == 0.0f ? -eslINFINITY : log(v[k]));
  else        for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k]);
  fprintf(ox->dfp, "\n");

  for (q = 0; q < Q; q++) {
    tmp.v = DMXo(q);
    for (z = 0; z < 4; z++) v[q + Q*z + 1] = tmp.x[z];
  }
  fprintf(ox->dfp, "%3d D ", rowi);
  if (logify) for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k] == 0.0f ? -eslINFINITY : log(v[k]));
  else        for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k]);
  fprintf(ox->dfp, "\n\n");

  free(v);
  return eslOK;

 ERROR:
  free(v);
  return status;
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_omx_DumpFBRow_FS()
 *
 * Purpose:   Like p7_omx_DumpFBRow() but for frameshift (8-cell) matrices.
 *            <i> is used for row label; header printed when i == 0.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_omx_DumpFBRow_FS(P7_OMX *ox, int logify, int i, int rowi, int width, int precision,
                    float xE, float xN, float xJ, float xB, float xC)
{
#ifdef eslENABLE_SSE
  __m128 *dp;
  int     M  = ox->M;
  int     Q  = p7O_NQF(M);
  float  *v  = NULL;
  int     q, z, k;
  union { __m128 v; float x[4]; } tmp;
  int     status;

  dp = (ox->allocR == 1) ? ox->dpf[0] : ox->dpf[rowi];

  ESL_ALLOC(v, sizeof(float) * (Q * 4 + 1));
  v[0] = 0.0f;

  if (i == 0)
    {
      fprintf(ox->dfp, "      ");
      for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*d ", width, k);
      fprintf(ox->dfp, "%*s %*s %*s %*s %*s\n",
              width, "E", width, "N", width, "J", width, "B", width, "C");
      fprintf(ox->dfp, "      ");
      for (k = 0; k <= M + 5; k++) fprintf(ox->dfp, "%*s ", width, "--------");
      fprintf(ox->dfp, "\n");
    }

  for (q = 0; q < Q; q++) {
    tmp.v = MMXo(q);
    for (z = 0; z < 4; z++) v[q + Q*z + 1] = tmp.x[z];
  }
  fprintf(ox->dfp, "%3d M ", i);
  if (logify) for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k] == 0.0f ? -eslINFINITY : log(v[k]));
  else        for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k]);
  if (logify)
    fprintf(ox->dfp, "%*.*f %*.*f %*.*f %*.*f %*.*f\n",
            width, precision, xE == 0.0f ? -eslINFINITY : log(xE),
            width, precision, xN == 0.0f ? -eslINFINITY : log(xN),
            width, precision, xJ == 0.0f ? -eslINFINITY : log(xJ),
            width, precision, xB == 0.0f ? -eslINFINITY : log(xB),
            width, precision, xC == 0.0f ? -eslINFINITY : log(xC));
  else
    fprintf(ox->dfp, "%*.*f %*.*f %*.*f %*.*f %*.*f\n",
            width, precision, xE, width, precision, xN,
            width, precision, xJ, width, precision, xB,
            width, precision, xC);

  for (q = 0; q < Q; q++) {
    tmp.v = IMXo(q);
    for (z = 0; z < 4; z++) v[q + Q*z + 1] = tmp.x[z];
  }
  fprintf(ox->dfp, "%3d I ", i);
  if (logify) for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k] == 0.0f ? -eslINFINITY : log(v[k]));
  else        for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k]);
  fprintf(ox->dfp, "\n");

  for (q = 0; q < Q; q++) {
    tmp.v = DMXo(q);
    for (z = 0; z < 4; z++) v[q + Q*z + 1] = tmp.x[z];
  }
  fprintf(ox->dfp, "%3d D ", i);
  if (logify) for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k] == 0.0f ? -eslINFINITY : log(v[k]));
  else        for (k = 0; k <= M; k++) fprintf(ox->dfp, "%*.*f ", width, precision, v[k]);
  fprintf(ox->dfp, "\n\n");

  free(v);
  return eslOK;

 ERROR:
  free(v);
  return status;
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_omx_Dump()
 *
 * Purpose:   Dump the complete float (MDI) DP matrix stored in <ox>
 *            to stream <fp>.  <ox> must be a full-matrix allocation.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_omx_Dump(FILE *fp, P7_OMX *ox)
{
#ifdef eslENABLE_SSE
  int    M  = ox->M;
  int    L  = ox->L;
  int    Q  = p7O_NQF(M);
  int    width     = 9;
  int    precision = 4;
  float *v  = NULL;
  int    i, q, z, k;
  union { __m128 v; float x[4]; } tmp;
  int    status;

  ESL_ALLOC(v, sizeof(float) * (Q * 4 + 1));

  fprintf(fp, "      ");
  for (k = 0; k <= M; k++) fprintf(fp, "%*d ", width, k);
  fprintf(fp, "%*s %*s %*s %*s %*s\n",
          width, "E", width, "N", width, "J", width, "B", width, "C");
  fprintf(fp, "      ");
  for (k = 0; k <= M + 5; k++) fprintf(fp, "%*.*s ", width, width, "----------");
  fprintf(fp, "\n");

  for (i = 0; i <= L; i++)
    {
      __m128 *dp = ox->dpf[i];

      v[0] = 0.0f;
      for (q = 0; q < Q; q++) {
        tmp.v = MMO(dp, q);
        for (z = 0; z < 4; z++) { k = q + z*Q + 1; if (k <= M) v[k] = tmp.x[z]; }
      }
      fprintf(fp, "%3d M ", i);
      for (k = 0; k <= M; k++) fprintf(fp, "%*.*f ", width, precision, v[k]);
      if (ox->xmx != NULL)
        fprintf(fp, "%*.*f %*.*f %*.*f %*.*f %*.*f",
                width, precision, ox->xmx[i * p7X_NXCELLS + p7X_E],
                width, precision, ox->xmx[i * p7X_NXCELLS + p7X_N],
                width, precision, ox->xmx[i * p7X_NXCELLS + p7X_J],
                width, precision, ox->xmx[i * p7X_NXCELLS + p7X_B],
                width, precision, ox->xmx[i * p7X_NXCELLS + p7X_C]);
      fprintf(fp, "\n");

      v[0] = 0.0f;
      for (q = 0; q < Q; q++) {
        tmp.v = IMO(dp, q);
        for (z = 0; z < 4; z++) { k = q + z*Q + 1; if (k <= M) v[k] = tmp.x[z]; }
      }
      fprintf(fp, "%3d I ", i);
      for (k = 0; k <= M; k++) fprintf(fp, "%*.*f ", width, precision, v[k]);
      fprintf(fp, "\n");

      v[0] = 0.0f;
      for (q = 0; q < Q; q++) {
        tmp.v = DMO(dp, q);
        for (z = 0; z < 4; z++) { k = q + z*Q + 1; if (k <= M) v[k] = tmp.x[z]; }
      }
      fprintf(fp, "%3d D ", i);
      for (k = 0; k <= M; k++) fprintf(fp, "%*.*f ", width, precision, v[k]);
      fprintf(fp, "\n\n");
    }

  free(v);
  return eslOK;

 ERROR:
  free(v);
  return status;
#else
  return eslENORESULT;
#endif
}
