/* Null2 biased composition correction; SSE implementations for impl_avx dispatch.
 * These are the _sse-suffixed versions called by the runtime dispatchers
 * in null2.c and null2_fs.c.
 */

#include "p7_config.h"

#ifdef eslENABLE_SSE

#include <stdlib.h>
#include <string.h>

#include <xmmintrin.h>
#include <emmintrin.h>

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Null2_ByExpectation_sse()
 *
 * Purpose:   SSE implementation of null2 estimation from posterior probabilities.
 *            See p7_Null2_ByExpectation() documentation for details.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_Null2_ByExpectation_sse(const P7_OPROFILE *om, const P7_OMX *pp, float *null2)
{
  int      Ld   = pp->L;
  int      Q    = p7O_NQF(om->M);
  float   *xmx  = pp->xmx;
  int      i,q;

  memcpy(pp->dpf[0], pp->dpf[1], sizeof(__m128) * 3 * Q);
  XMXo(0,p7X_N) = XMXo(1,p7X_N);
  XMXo(0,p7X_C) = XMXo(1,p7X_C);
  XMXo(0,p7X_J) = XMXo(1,p7X_J);

  for (i = 2; i <= Ld; i++)
    {
      for (q = 0; q < Q; q++)
        {
          pp->dpf[0][q*3 + p7X_M] = _mm_add_ps(pp->dpf[i][q*3 + p7X_M], pp->dpf[0][q*3 + p7X_M]);
          pp->dpf[0][q*3 + p7X_I] = _mm_add_ps(pp->dpf[i][q*3 + p7X_I], pp->dpf[0][q*3 + p7X_I]);
        }
      XMXo(0,p7X_N) += XMXo(i,p7X_N);
      XMXo(0,p7X_C) += XMXo(i,p7X_C);
      XMXo(0,p7X_J) += XMXo(i,p7X_J);
    }

  return p7_Null2_FromCounts_sse(om, pp->dpf[0], pp->xmx, Ld, null2);
}


/* Function:  p7_Null2_FromCounts_sse()
 * Synopsis:  Null2 odds ratios from expected state usage, SSE.
 *
 * Purpose:   The second half of <p7_Null2_ByExpectation_sse()>: from the
 *            posterior probabilities of M and I summed over rows 1..Ld
 *            of an envelope of <Ld> residues (<cnt>, one row of striped
 *            cells), and those of N, C and J (<cntx>, indexed by
 *            p7X_N etc.), calculate the null2 odds ratios into <null2>
 *            (<0..Kp-1>, caller allocated). <cnt> and <cntx> are
 *            overwritten.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_Null2_FromCounts_sse(const P7_OPROFILE *om, __m128 *cnt, float *cntx, int Ld, float *null2)
{
  int      Q    = p7O_NQF(om->M);
  float    norm;
  __m128  *rp;
  __m128   sv;
  float    xfactor;
  int      q,x;

  norm = 1.0 / (float) Ld;
  sv   = _mm_set1_ps(norm);
  for (q = 0; q < Q; q++)
    {
      cnt[q*3 + p7X_M] = _mm_mul_ps(cnt[q*3 + p7X_M], sv);
      cnt[q*3 + p7X_I] = _mm_mul_ps(cnt[q*3 + p7X_I], sv);
    }
  cntx[p7X_N] *= norm;
  cntx[p7X_C] *= norm;
  cntx[p7X_J] *= norm;

  xfactor = cntx[p7X_N] + cntx[p7X_C] + cntx[p7X_J];
  for (x = 0; x < om->abc->K; x++)
    {
      sv = _mm_setzero_ps();
      rp = om->rfv[x];
      for (q = 0; q < Q; q++)
        {
          sv = _mm_add_ps(sv, _mm_mul_ps(cnt[q*3 + p7X_M], *rp)); rp++;
          sv = _mm_add_ps(sv,            cnt[q*3 + p7X_I]);  /* insert odds implicitly 1.0 */
        }
      esl_sse_hsum_ps(sv, &(null2[x]));
      null2[x] += xfactor;
    }

  esl_abc_FAvgScVec(om->abc, null2);
  null2[om->abc->K]    = 1.0;
  null2[om->abc->Kp-2] = 1.0;
  null2[om->abc->Kp-1] = 1.0;

  return eslOK;
}


#endif /* eslENABLE_SSE */
