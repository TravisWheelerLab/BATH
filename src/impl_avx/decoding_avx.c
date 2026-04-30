/* Posterior decoding algorithms; AVX2 implementations for impl_avx dispatch.
 * These are the _avx-suffixed versions called by the runtime dispatchers
 * in decoding.c and decoding_fs.c.
 */

#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <stdio.h>
#include <math.h>

#include <immintrin.h>

#include "easel.h"
#include "esl_avx.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Decoding_avx()
 *
 * Purpose:   AVX2 implementation of posterior decoding of residue assignment.
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
p7_Decoding_avx(const P7_OPROFILE *om, const P7_OMX *oxf, P7_OMX *oxb, P7_OMX *pp)
{
  __m256 *ppv;
  __m256 *fv;
  __m256 *bv;
  __m256  totrv;
  int    L  = oxf->L;
  int    M  = om->M;
  int    Q  = p7O_NQF_AVX(M);
  int    i,q;
  float  scaleproduct = 1.0 / oxb->xmx[p7X_N];

  pp->M = M;
  pp->L = L;

  ppv = pp->dpf_avx[0];
  for (q = 0; q < Q; q++) {
    *ppv = _mm256_setzero_ps(); ppv++;
    *ppv = _mm256_setzero_ps(); ppv++;
    *ppv = _mm256_setzero_ps(); ppv++;
  }
  pp->xmx[p7X_E] = 0.0;
  pp->xmx[p7X_N] = 0.0;
  pp->xmx[p7X_J] = 0.0;
  pp->xmx[p7X_C] = 0.0;
  pp->xmx[p7X_B] = 0.0;

  for (i = 1; i <= L; i++)
    {
      ppv   =  pp->dpf_avx[i];
      fv    = oxf->dpf_avx[i];
      bv    = oxb->dpf_avx[i];
      totrv = _mm256_set1_ps(scaleproduct * oxf->xmx[i*p7X_NXCELLS+p7X_SCALE]);

      for (q = 0; q < Q; q++)
        {
          /* M */
          *ppv = _mm256_mul_ps(*fv,  *bv);
          *ppv = _mm256_mul_ps(*ppv,  totrv);
          ppv++;  fv++;  bv++;

          /* D */
          *ppv = _mm256_setzero_ps();
          ppv++;  fv++;  bv++;

          /* I */
          *ppv = _mm256_mul_ps(*fv,  *bv);
          *ppv = _mm256_mul_ps(*ppv,  totrv);
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


/* Function:  p7_DomainDecoding_avx()
 *
 * Purpose:   AVX2 implementation of posterior decoding of domain location.
 *            See p7_DomainDecoding() documentation in decoding.c for full details.
 *            This implementation is scalar (no SIMD) — it only accesses xmx arrays.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> on numeric overflow.
 */
int
p7_DomainDecoding_avx(const P7_OPROFILE *om, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef)
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


#endif /* eslENABLE_AVX */
