/* Frameshift posterior decoding algorithms; SSE implementations.
 * Ported from impl_sse/decoding_fs.c with _sse suffix for runtime dispatch.
 *
 * Contents:
 *   1. p7_Decoding_Frameshift_sse()
 *   2. p7_DomainDecoding_Frameshift_sse()
 */
#include "p7_config.h"

#include <string.h>

#ifdef eslENABLE_SSE

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>  /* SSE  */
#include <emmintrin.h>  /* SSE2 */

#include "easel.h"
#include "esl_sse.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_DecodingRow_Frameshift_sse()
 * Synopsis:  Posterior decoding of one frameshift row.
 *
 * Purpose:   Posterior probabilities of row <i> (1..L) into <ppv> and the
 *            special states <ppx> (E, N, J, B, C of the row), from Forward
 *            row <i> in <fv>, Backward row <i> in <bv>, the Backward
 *            special states <bxmx>, the Forward special states of row
 *            i-3 in <fx3> (read for i > 2), the cumulative log scale
 *            products <log_sfwd> and <log_sbck>, and <log_inv_Z>, as
 *            <p7_Decoding_Frameshift_sse()> computes them. <ppv> may be
 *            <fv>. Same calculation as each row of
 *            <p7_Decoding_Frameshift_sse()>, which calls it.
 *
 * Returns:   <eslOK> on success; <eslERANGE> if a scale factor or the
 *            special states' posteriors overflow.
 */
int
p7_DecodingRow_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const float *bxmx, const float *fx3, const float *log_sfwd, const float *log_sbck,
                         float log_inv_Z, const __m128 *fv, const __m128 *bv, __m128 *ppv, float *ppx, int i)
{
  int    Q = p7O_NQF(om_fs->M);
  int    q, c;
  float  factor_mdi, factor_njc;
  __m128 scv;
  __m128 zerov    = _mm_setzero_ps();
  __m128  fI, fC[6], bM, bI;
  float  N_odds = om_fs->xf[p7O_N][p7O_LOOP];
  float  J_odds = om_fs->xf[p7O_J][p7O_LOOP];
  float  C_odds = om_fs->xf[p7O_C][p7O_LOOP];
  float  bck_N, bck_J, bck_C;
  float  N_pp, J_pp, C_pp;

  /* Backward special states at row i */
  bck_N = bxmx[i*p7X_NXCELLS + p7X_N];
  bck_J = bxmx[i*p7X_NXCELLS + p7X_J];
  bck_C = bxmx[i*p7X_NXCELLS + p7X_C];

  /* Scale factor for same-row M/I products */
  factor_mdi = expf(log_sfwd[i] + log_sbck[i] + log_inv_Z);

  if (isinf(factor_mdi)) return eslERANGE;

  /* M/I products f * b; scaled to posteriors below */
  for (q = 0; q < Q; q++)
    {
      bM = bv[q * p7X_NSCELLS + p7X_M];
      bI = bv[q * p7X_NSCELLS + p7X_I];

      fI = fv[q * p7X_NSCELLS_FS + p7X_FS_I];
      for (c = 0; c < 6; c++)
        fC[c] = fv[q * p7X_NSCELLS_FS + p7X_FS_M + c];

      ppv[q * p7X_NSCELLS_FS + p7X_FS_D] = zerov;
      ppv[q * p7X_NSCELLS_FS + p7X_FS_I] = _mm_mul_ps(fI, bI);
      for (c = 0; c < 6; c++)
        ppv[q * p7X_NSCELLS_FS + p7X_FS_M + c] = _mm_mul_ps(fC[c], bM);
    }

  /* N, J, C posteriors */
  if (i > 2) {
    factor_njc = expf(log_sfwd[i-3] + log_sbck[i] + log_inv_Z);
    N_pp = fx3[p7X_N] * bck_N * N_odds * factor_njc;
    J_pp = fx3[p7X_J] * bck_J * J_odds * factor_njc;
    C_pp = fx3[p7X_C] * bck_C * C_odds * factor_njc;
  } else {
    float factor_nsmall = expf(log_sbck[i] + log_inv_Z);
    N_pp = bck_N * factor_nsmall;
    J_pp = 0.0f;
    C_pp = 0.0f;
  }
  if (isinf(N_pp + J_pp + C_pp)) return eslERANGE;

  /* Scale the M/I products by the row factor to get posteriors */
  scv = _mm_set1_ps(factor_mdi);
  for (q = 0; q < Q; q++)
    for (c = 0; c < p7X_NSCELLS_FS; c++)
      ppv[q * p7X_NSCELLS_FS + c] = _mm_mul_ps(ppv[q * p7X_NSCELLS_FS + c], scv);

  /* Special states */
  ppx[p7X_E] = 0.0f;
  ppx[p7X_B] = 0.0f;
  ppx[p7X_N] = N_pp;
  ppx[p7X_J] = J_pp;
  ppx[p7X_C] = C_pp;
  return eslOK;
}


/* Function:  p7_Decoding_Frameshift_sse()
 *
 * Purpose:   SSE implementation of frameshift posterior decoding.
 *            Given filled forward matrix <fwd> and backward matrix <bck>,
 *            overwrites <fwd> in-place with the posterior decoding matrix.
 *
 *            Each value at row i is the posterior probability that the
 *            state's emission ends at nucleotide i (M_C0..M_C5, I, N, J,
 *            C; D is zero). Emissions span 1-5 nucleotides, so the values
 *            at one row do not sum to one; over all rows they sum to the
 *            expected number of emissions.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if numeric range is exceeded.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_Decoding_Frameshift_sse(const P7_FS_OPROFILE *om_fs, P7_OMX *fwd, const P7_OMX *bck)
{
  int    L = fwd->L;
  int    M = om_fs->M;
  int    Q = p7O_NQF(M);
  int    i, q;
  float *log_sfwd = NULL;
  float *log_sbck = NULL;
  float  log_inv_Z;
  __m128 zerov    = _mm_setzero_ps();
  __m128 *dpc;
  float  lagx[4][p7X_NXCELLS];          /* fwd special states of the last 4 rows, saved before overwrite */
  int    status;

  ESL_ALLOC(log_sfwd, sizeof(float) * (L + 2));
  ESL_ALLOC(log_sbck, sizeof(float) * (L + 2));

  /* Cumulative log scale products */
  log_sfwd[0] = logf(fwd->xmx[p7X_SCALE]);
  for (i = 1; i <= L; i++)
    log_sfwd[i] = log_sfwd[i-1] + logf(fwd->xmx[i*p7X_NXCELLS + p7X_SCALE]);

  log_sbck[L+1] = 0.0f;
  for (i = L; i >= 0; i--)
    log_sbck[i] = log_sbck[i+1] + logf(bck->xmx[i*p7X_NXCELLS + p7X_SCALE]);

  log_inv_Z = -p7_FLogsum(
                 logf(bck->xmx[p7X_N]) + log_sbck[0],
                 p7_FLogsum(
                   logf(bck->xmx[1*p7X_NXCELLS + p7X_N]) + log_sbck[1],
                   logf(bck->xmx[2*p7X_NXCELLS + p7X_N]) + log_sbck[2]));

  /* Save fwd row 0 special states before zeroing row 0 */
  memcpy(lagx[0], fwd->xmx, sizeof(float) * p7X_NXCELLS);

  /* Zero row 0 of the output PP matrix */
  dpc = fwd->dpf[0];
  for (q = 0; q < Q * p7X_NSCELLS_FS; q++) dpc[q] = zerov;
  fwd->xmx[p7X_E] = fwd->xmx[p7X_N] = fwd->xmx[p7X_J] = fwd->xmx[p7X_B] = fwd->xmx[p7X_C] = 0.0f;

  for (i = 1; i <= L; i++)
    {
      /* Save current fwd special states into lag buffer before overwriting */
      memcpy(lagx[i % 4], fwd->xmx + i*p7X_NXCELLS, sizeof(float) * p7X_NXCELLS);

      /* fwd special states 3 rows back: (i-3+4)%4 == (i+1)%4 */
      if ((status = p7_DecodingRow_Frameshift_sse(om_fs, bck->xmx, lagx[(i + 1) % 4], log_sfwd, log_sbck, log_inv_Z,
                                              fwd->dpf[i], bck->dpf[i], fwd->dpf[i], fwd->xmx + i*p7X_NXCELLS, i)) != eslOK)
        { free(log_sfwd); free(log_sbck); return status; }
    }

  free(log_sfwd);
  free(log_sbck);
  return eslOK;

 ERROR:
  if (log_sfwd) free(log_sfwd);
  if (log_sbck) free(log_sbck);
  return status;
}


/* Function:  p7_DomainDecoding_Frameshift_sse()
 *
 * Purpose:   SSE implementation of frameshift domain decoding.
 *            Fills <ddef->btot[i]>, <ddef->etot[i]>, <ddef->mocc[i]>
 *            for i = 0..L, then sets <ddef->L>.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_DomainDecoding_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *oxf, const P7_OMX *oxb,
                                  P7_DOMAINDEF *ddef)
{
  int    L   = oxf->L;
  float *log_sfwd = NULL;
  float *log_sbck = NULL;
  float  log_inv_Z;
  float  njcp;
  int    i;
  int    status;

  ESL_ALLOC(log_sfwd, sizeof(float) * (L+2));
  ESL_ALLOC(log_sbck, sizeof(float) * (L+2));

  /* Cumulative log forward scale products */
  log_sfwd[0] = logf(oxf->xmx[0*p7X_NXCELLS + p7X_SCALE]);
  for (i = 1; i <= L; i++)
    log_sfwd[i] = log_sfwd[i-1] + logf(oxf->xmx[i*p7X_NXCELLS + p7X_SCALE]);

  /* Cumulative log backward scale products */
  log_sbck[L+1] = 0.0f;
  for (i = L; i >= 0; i--)
    log_sbck[i] = log_sbck[i+1] + logf(oxb->xmx[i*p7X_NXCELLS + p7X_SCALE]);

  log_inv_Z = -p7_FLogsum(
                 logf(oxb->xmx[0*p7X_NXCELLS + p7X_N]) + log_sbck[0],
                 p7_FLogsum(
                   logf(oxb->xmx[1*p7X_NXCELLS + p7X_N]) + log_sbck[1],
                   logf(oxb->xmx[2*p7X_NXCELLS + p7X_N]) + log_sbck[2]));

  /* Positions 0, 1, 2: no domains can have started or ended yet */
  ddef->btot[0] = 0.;
  ddef->btot[1] = 0.;
  ddef->btot[2] = 0.;
  ddef->etot[0] = 0.;
  ddef->etot[1] = 0.;
  ddef->etot[2] = 0.;
  ddef->mocc[0] = 0.;
  ddef->mocc[1] = 0.;
  ddef->mocc[2] = 0.;

  for (i = 3; i <= L; i++)
    {
      ddef->btot[i] = ddef->btot[i-3]
        + oxf->xmx[(i-3)*p7X_NXCELLS + p7X_B] * oxb->xmx[(i-3)*p7X_NXCELLS + p7X_B]
          * expf(log_sfwd[i-3] + log_sbck[i-3] + log_inv_Z);

      ddef->etot[i] = ddef->etot[i-3]
        + oxf->xmx[i*p7X_NXCELLS + p7X_E] * oxb->xmx[i*p7X_NXCELLS + p7X_E]
          * expf(log_sfwd[i] + log_sbck[i] + log_inv_Z);

      njcp = 0.;

      /* N state */
      njcp += oxf->xmx[(i-3)*p7X_NXCELLS + p7X_N] * oxb->xmx[i*p7X_NXCELLS + p7X_N]
              * om_fs->xf[p7O_N][p7O_LOOP] * expf(log_sfwd[i-3] + log_sbck[i] + log_inv_Z);
      if (i < L)
        njcp += oxf->xmx[(i-2)*p7X_NXCELLS + p7X_N] * oxb->xmx[(i+1)*p7X_NXCELLS + p7X_N]
                * om_fs->xf[p7O_N][p7O_LOOP] * expf(log_sfwd[i-2] + log_sbck[i+1] + log_inv_Z);
      if (i < L-1)
        njcp += oxf->xmx[(i-1)*p7X_NXCELLS + p7X_N] * oxb->xmx[(i+2)*p7X_NXCELLS + p7X_N]
                * om_fs->xf[p7O_N][p7O_LOOP] * expf(log_sfwd[i-1] + log_sbck[i+2] + log_inv_Z);

      /* J state */
      njcp += oxf->xmx[(i-3)*p7X_NXCELLS + p7X_J] * oxb->xmx[i*p7X_NXCELLS + p7X_J]
              * om_fs->xf[p7O_J][p7O_LOOP] * expf(log_sfwd[i-3] + log_sbck[i] + log_inv_Z);
      if (i < L)
        njcp += oxf->xmx[(i-2)*p7X_NXCELLS + p7X_J] * oxb->xmx[(i+1)*p7X_NXCELLS + p7X_J]
                * om_fs->xf[p7O_J][p7O_LOOP] * expf(log_sfwd[i-2] + log_sbck[i+1] + log_inv_Z);
      if (i < L-1)
        njcp += oxf->xmx[(i-1)*p7X_NXCELLS + p7X_J] * oxb->xmx[(i+2)*p7X_NXCELLS + p7X_J]
                * om_fs->xf[p7O_J][p7O_LOOP] * expf(log_sfwd[i-1] + log_sbck[i+2] + log_inv_Z);

      /* C state */
      njcp += oxf->xmx[(i-3)*p7X_NXCELLS + p7X_C] * oxb->xmx[i*p7X_NXCELLS + p7X_C]
              * om_fs->xf[p7O_C][p7O_LOOP] * expf(log_sfwd[i-3] + log_sbck[i] + log_inv_Z);
      if (i < L)
        njcp += oxf->xmx[(i-2)*p7X_NXCELLS + p7X_C] * oxb->xmx[(i+1)*p7X_NXCELLS + p7X_C]
                * om_fs->xf[p7O_C][p7O_LOOP] * expf(log_sfwd[i-2] + log_sbck[i+1] + log_inv_Z);
      if (i < L-1)
        njcp += oxf->xmx[(i-1)*p7X_NXCELLS + p7X_C] * oxb->xmx[(i+2)*p7X_NXCELLS + p7X_C]
                * om_fs->xf[p7O_C][p7O_LOOP] * expf(log_sfwd[i-1] + log_sbck[i+2] + log_inv_Z);

      ddef->mocc[i] = 1. - njcp;
    }

  ddef->L = L;
  free(log_sfwd);
  free(log_sbck);
  return eslOK;

 ERROR:
  if (log_sfwd) free(log_sfwd);
  if (log_sbck) free(log_sbck);
  return status;
}

#endif /* eslENABLE_SSE */
