/* Posterior decoding algorithms; SSE implementations for impl_avx dispatch.
 * These are the _sse-suffixed versions called by the runtime dispatchers
 * in decoding.c and decoding_fs.c.
 */

#ifdef eslENABLE_SSE
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>

#include "easel.h"
#include "esl_sse.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Decoding_sse()
 *
 * Purpose:   SSE implementation of posterior decoding of residue assignment.
 *            See p7_Decoding() documentation in decoding.c for full details.
 *
 * Args:      om   - profile
 *            oxf  - filled Forward matrix
 *            oxb  - filled Backward matrix (overwritten with PP on success)
 *            pp   - RESULT: posterior decoding matrix
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if numeric range is exceeded.
 */
int
p7_Decoding_sse(const P7_OPROFILE *om, const P7_OMX *oxf, P7_OMX *oxb, P7_OMX *pp)
{
  __m128 *ppv;
  __m128 *fv;
  __m128 *bv;
  __m128  totrv;
  int    L  = oxf->L;
  int    M  = om->M;
  int    Q  = p7O_NQF(M);
  int    i,q;
  float  scaleproduct = 1.0 / oxb->xmx[p7X_N];

  pp->M = M;
  pp->L = L;

  ppv = pp->dpf[0];
  for (q = 0; q < Q; q++) {
    *ppv = _mm_setzero_ps(); ppv++;
    *ppv = _mm_setzero_ps(); ppv++;
    *ppv = _mm_setzero_ps(); ppv++;
  }
  pp->xmx[p7X_E] = 0.0;
  pp->xmx[p7X_N] = 0.0;
  pp->xmx[p7X_J] = 0.0;
  pp->xmx[p7X_C] = 0.0;
  pp->xmx[p7X_B] = 0.0;

  for (i = 1; i <= L; i++)
    {
      ppv   =  pp->dpf[i];
      fv    = oxf->dpf[i];
      bv    = oxb->dpf[i];
      totrv = _mm_set1_ps(scaleproduct * oxf->xmx[i*p7X_NXCELLS+p7X_SCALE]);

      for (q = 0; q < Q; q++)
        {
          /* M */
          *ppv = _mm_mul_ps(*fv,  *bv);
          *ppv = _mm_mul_ps(*ppv,  totrv);
          ppv++;  fv++;  bv++;

          /* D */
          *ppv = _mm_setzero_ps();
          ppv++;  fv++;  bv++;

          /* I */
          *ppv = _mm_mul_ps(*fv,  *bv);
          *ppv = _mm_mul_ps(*ppv,  totrv);
          ppv++;  fv++;  bv++;
        }
      pp->xmx[i*p7X_NXCELLS+p7X_E] = 0.0;
      pp->xmx[i*p7X_NXCELLS+p7X_N] = oxf->xmx[(i-1)*p7X_NXCELLS+p7X_N] * oxb->xmx[i*p7X_NXCELLS+p7X_N] * om->xf[p7O_N][p7O_LOOP] * scaleproduct;
      pp->xmx[i*p7X_NXCELLS+p7X_J] = oxf->xmx[(i-1)*p7X_NXCELLS+p7X_J] * oxb->xmx[i*p7X_NXCELLS+p7X_J] * om->xf[p7O_J][p7O_LOOP] * scaleproduct;
      pp->xmx[i*p7X_NXCELLS+p7X_C] = oxf->xmx[(i-1)*p7X_NXCELLS+p7X_C] * oxb->xmx[i*p7X_NXCELLS+p7X_C] * om->xf[p7O_C][p7O_LOOP] * scaleproduct;
      pp->xmx[i*p7X_NXCELLS+p7X_B] = 0.0;

      if (oxb->has_own_scales) scaleproduct *= oxf->xmx[i*p7X_NXCELLS+p7X_SCALE] /  oxb->xmx[i*p7X_NXCELLS+p7X_SCALE];
    }

  if (isinf(scaleproduct)) return eslERANGE;
  else                     return eslOK;
}


/* Function:  p7_DomainDecoding_sse()
 *
 * Purpose:   SSE implementation of posterior decoding of domain location.
 *            See p7_DomainDecoding() documentation in decoding.c for full details.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> on numeric overflow.
 */
int
p7_DomainDecoding_sse(const P7_OPROFILE *om, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef)
{
  int   L             = oxf->L;
  float scaleproduct  = 1.0 / oxb->xmx[p7X_N];
  float njcp;
  int   i;

  ddef->btot[0] = 0.0;
  ddef->etot[0] = 0.0;
  ddef->mocc[0] = 0.0;
  for (i = 1; i <= L; i++)
    {
      /* scaleproduct is prod_j=0^i-2 now */
      ddef->btot[i] = ddef->btot[i-1] +
        (oxf->xmx[(i-1)*p7X_NXCELLS+p7X_B] * oxb->xmx[(i-1)*p7X_NXCELLS+p7X_B] * oxf->xmx[(i-1)*p7X_NXCELLS+p7X_SCALE] * scaleproduct);

      if (oxb->has_own_scales) scaleproduct *= oxf->xmx[(i-1)*p7X_NXCELLS+p7X_SCALE] /  oxb->xmx[(i-1)*p7X_NXCELLS+p7X_SCALE];
      /* scaleproduct is prod_j=0^i-1 now */

      ddef->etot[i] = ddef->etot[i-1] +
        (oxf->xmx[i*p7X_NXCELLS+p7X_E] * oxb->xmx[i*p7X_NXCELLS+p7X_E] * oxf->xmx[i*p7X_NXCELLS+p7X_SCALE] * scaleproduct);

      njcp  = oxf->xmx[(i-1)*p7X_NXCELLS+p7X_N] * oxb->xmx[i*p7X_NXCELLS+p7X_N] * om->xf[p7O_N][p7O_LOOP] * scaleproduct;
      njcp += oxf->xmx[(i-1)*p7X_NXCELLS+p7X_J] * oxb->xmx[i*p7X_NXCELLS+p7X_J] * om->xf[p7O_J][p7O_LOOP] * scaleproduct;
      njcp += oxf->xmx[(i-1)*p7X_NXCELLS+p7X_C] * oxb->xmx[i*p7X_NXCELLS+p7X_C] * om->xf[p7O_C][p7O_LOOP] * scaleproduct;
      ddef->mocc[i] = 1. - njcp;
    }
  ddef->L = oxf->L;

  if (isinf(scaleproduct)) return eslERANGE;
  else                     return eslOK;
}


/* Function:  p7_Decoding_Frameshift_sse()
 *
 * Purpose:   SSE implementation of frameshift-aware posterior decoding.
 *            Overwrites <fwd> with posterior probability matrix in-place.
 *            See p7_Decoding_Frameshift() documentation for full details.
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
  int    i, q, c;
  float *log_sfwd = NULL;
  float *log_sbck = NULL;
  float  log_inv_Z;
  float  factor_mdi, factor_njc;
  __m128 scv, dv;
  __m128 zerov    = _mm_setzero_ps();
  float  raw_denom;
  __m128 *dpc, *dbv;
  __m128  fI, fC[6], bM, bI;
  /* Circular lag buffer depth 4: nlag[i%4] = fwd_N(i) saved before overwrite */
  float  nlag[4], jlag[4], clag[4];
  float  fN3, fJ3, fC3;
  float  bck_N, bck_J, bck_C;
  float  N_odds = om_fs->xf[p7O_N][p7O_LOOP];
  float  J_odds = om_fs->xf[p7O_J][p7O_LOOP];
  float  C_odds = om_fs->xf[p7O_C][p7O_LOOP];
  float  N_pp, J_pp, C_pp, inv_denom;
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
  nlag[0] = fwd->xmx[p7X_N];
  jlag[0] = fwd->xmx[p7X_J];
  clag[0] = fwd->xmx[p7X_C];
  nlag[1] = nlag[2] = nlag[3] = 0.0f;
  jlag[1] = jlag[2] = jlag[3] = 0.0f;
  clag[1] = clag[2] = clag[3] = 0.0f;

  /* Zero row 0 of the output PP matrix */
  dpc = fwd->dpf[0];
  for (q = 0; q < Q * p7X_NSCELLS_FS; q++) dpc[q] = zerov;
  fwd->xmx[p7X_E] = fwd->xmx[p7X_N] = fwd->xmx[p7X_J] = fwd->xmx[p7X_B] = fwd->xmx[p7X_C] = 0.0f;

  for (i = 1; i <= L; i++)
    {
      /* Save current fwd special states into lag buffer before overwriting */
      nlag[i % 4] = fwd->xmx[i*p7X_NXCELLS + p7X_N];
      jlag[i % 4] = fwd->xmx[i*p7X_NXCELLS + p7X_J];
      clag[i % 4] = fwd->xmx[i*p7X_NXCELLS + p7X_C];

      /* fwd special states 3 rows back: (i-3+4)%4 == (i+1)%4 */
      fN3 = nlag[(i + 1) % 4];
      fJ3 = jlag[(i + 1) % 4];
      fC3 = clag[(i + 1) % 4];

      /* Backward special states at row i */
      bck_N = bck->xmx[i*p7X_NXCELLS + p7X_N];
      bck_J = bck->xmx[i*p7X_NXCELLS + p7X_J];
      bck_C = bck->xmx[i*p7X_NXCELLS + p7X_C];

      /* Scale factor for same-row M/I products */
      factor_mdi = expf(log_sfwd[i] + log_sbck[i] + log_inv_Z);

      if (isinf(factor_mdi)) { free(log_sfwd); free(log_sbck); return eslERANGE; }

      /* ------- Pass 1: compute raw M/I products, accumulate denom ------- */
      dv  = zerov;
      dpc = fwd->dpf[i];
      dbv = bck->dpf[i];

      for (q = 0; q < Q; q++)
        {
          bM = dbv[q * p7X_NSCELLS + p7X_M];
          bI = dbv[q * p7X_NSCELLS + p7X_I];

          fI = dpc[q * p7X_NSCELLS_FS + p7X_FS_I];
          for (c = 0; c < 6; c++)
            fC[c] = dpc[q * p7X_NSCELLS_FS + p7X_FS_M + c];

          dpc[q * p7X_NSCELLS_FS + p7X_FS_D] = zerov;
          dpc[q * p7X_NSCELLS_FS + p7X_FS_I] = _mm_mul_ps(fI, bI);
          for (c = 0; c < 6; c++)
            dpc[q * p7X_NSCELLS_FS + p7X_FS_M + c] = _mm_mul_ps(fC[c], bM);

          /* Only C0 and I contribute to the per-row normalization denominator */
          dv = _mm_add_ps(dv, _mm_add_ps(dpc[q * p7X_NSCELLS_FS + p7X_FS_M],
                                          dpc[q * p7X_NSCELLS_FS + p7X_FS_I]));
        }

      /* Horizontal sum of dv to get scalar raw MDI denom */
      dv = _mm_add_ps(dv, _mm_movehl_ps(dv, dv));
      dv = _mm_add_ss(dv, _mm_shuffle_ps(dv, dv, _MM_SHUFFLE(0,0,0,1)));
      raw_denom = _mm_cvtss_f32(dv);

      /* N, J, C contributions to denominator */
      if (i > 2) {
        factor_njc = expf(log_sfwd[i-3] + log_sbck[i] + log_inv_Z);
        N_pp = fN3 * bck_N * N_odds * factor_njc;
        J_pp = fJ3 * bck_J * J_odds * factor_njc;
        C_pp = fC3 * bck_C * C_odds * factor_njc;
      } else {
        float factor_nsmall = expf(log_sbck[i] + log_inv_Z);
        N_pp = bck_N * factor_nsmall;
        J_pp = 0.0f;
        C_pp = 0.0f;
      }
      inv_denom = 1.0f / (raw_denom * factor_mdi + N_pp + J_pp + C_pp);

      if (isinf(inv_denom)) { free(log_sfwd); free(log_sbck); return eslERANGE; }

      /* ------- Pass 2: normalize all M/I vectors ------- */
      scv = _mm_set1_ps(factor_mdi * inv_denom);
      for (q = 0; q < Q; q++)
        for (c = 0; c < p7X_NSCELLS_FS; c++)
          fwd->dpf[i][q * p7X_NSCELLS_FS + c] = _mm_mul_ps(fwd->dpf[i][q * p7X_NSCELLS_FS + c], scv);

      /* Normalized special states */
      fwd->xmx[i*p7X_NXCELLS + p7X_E] = 0.0f;
      fwd->xmx[i*p7X_NXCELLS + p7X_B] = 0.0f;
      fwd->xmx[i*p7X_NXCELLS + p7X_N] = N_pp * inv_denom;
      fwd->xmx[i*p7X_NXCELLS + p7X_J] = J_pp * inv_denom;
      fwd->xmx[i*p7X_NXCELLS + p7X_C] = C_pp * inv_denom;
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
 * Purpose:   SSE implementation of frameshift-aware posterior decoding of
 *            domain location. Fills ddef->btot, etot, mocc arrays.
 *            See p7_DomainDecoding_Frameshift() documentation for full details.
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
      /* btot[i]: cumulative expected B-state entries up to position i */
      ddef->btot[i] = ddef->btot[i-3]
        + oxf->xmx[(i-3)*p7X_NXCELLS + p7X_B] * oxb->xmx[(i-3)*p7X_NXCELLS + p7X_B]
          * expf(log_sfwd[i-3] + log_sbck[i-3] + log_inv_Z);

      /* etot[i]: cumulative expected E-state exits up to position i */
      ddef->etot[i] = ddef->etot[i-3]
        + oxf->xmx[i*p7X_NXCELLS + p7X_E] * oxb->xmx[i*p7X_NXCELLS + p7X_E]
          * expf(log_sfwd[i] + log_sbck[i] + log_inv_Z);

      /* mocc[i]: 1 - P(residue i in N/J/C flanking states) */
      njcp = 0.;

      /* N state — three codon-offset contributions */
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
