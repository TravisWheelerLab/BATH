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
  int      M    = om->M;
  int      Ld   = pp->L;
  int      Q    = p7O_NQF(M);
  float   *xmx  = pp->xmx;
  float    norm;
  __m128  *rp;
  __m128   sv;
  float    xfactor;
  int      i,q,x;

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

  norm = 1.0 / (float) Ld;
  sv   = _mm_set1_ps(norm);
  for (q = 0; q < Q; q++)
    {
      pp->dpf[0][q*3 + p7X_M] = _mm_mul_ps(pp->dpf[0][q*3 + p7X_M], sv);
      pp->dpf[0][q*3 + p7X_I] = _mm_mul_ps(pp->dpf[0][q*3 + p7X_I], sv);
    }
  XMXo(0,p7X_N) *= norm;
  XMXo(0,p7X_C) *= norm;
  XMXo(0,p7X_J) *= norm;

  xfactor = XMXo(0, p7X_N) + XMXo(0, p7X_C) + XMXo(0, p7X_J);
  for (x = 0; x < om->abc->K; x++)
    {
      sv = _mm_setzero_ps();
      rp = om->rfv[x];
      for (q = 0; q < Q; q++)
        {
          sv = _mm_add_ps(sv, _mm_mul_ps(pp->dpf[0][q*3 + p7X_M], *rp)); rp++;
          sv = _mm_add_ps(sv,            pp->dpf[0][q*3 + p7X_I]);  /* insert odds implicitly 1.0 */
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


/* Function:  p7_Null2_ByTrace_sse()
 *
 * Purpose:   SSE implementation of null2 estimation from a trace (sampling method).
 *            See p7_Null2_ByTrace() documentation for details.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_Null2_ByTrace_sse(const P7_OPROFILE *om, const P7_TRACE *tr, int zstart, int zend,
                     P7_OMX *wrk, float *null2)
{
  union { __m128 v; float p[4]; } u;
  int    Q  = p7O_NQF(om->M);
  int    Ld = 0;
  float *xmx = wrk->xmx;
  float  norm;
  float  xfactor;
  __m128 sv;
  __m128 *rp;
  int    q, r;
  int    x;
  int    z;

  for (q = 0; q < Q; q++)
    {
      wrk->dpf[0][q*3 + p7X_M] = _mm_setzero_ps();
      wrk->dpf[0][q*3 + p7X_I] = _mm_setzero_ps();
    }
  XMXo(0,p7X_N) =  0.0;
  XMXo(0,p7X_C) =  0.0;
  XMXo(0,p7X_J) =  0.0;

  for (z = zstart; z <= zend; z++)
    {
      if (tr->i[z] == 0) continue;
      Ld++;
      if (tr->k[z] > 0)
        {
          q = p7X_NSCELLS * ( (tr->k[z] - 1) % Q) + p7X_M;
          r = (tr->k[z] - 1) / Q;
          u.v            = wrk->dpf[0][q];
          u.p[r]        += 1.0;
          wrk->dpf[0][q] = u.v;
        }
      else
        {
          switch (tr->st[z]) {
          case p7T_N: XMXo(0,p7X_N) += 1.0; break;
          case p7T_C: XMXo(0,p7X_C) += 1.0; break;
          case p7T_J: XMXo(0,p7X_J) += 1.0; break;
          }
        }
    }
  norm = 1.0 / (float) Ld;
  sv = _mm_set1_ps(norm);
  for (q = 0; q < Q; q++)
    {
      wrk->dpf[0][q*3 + p7X_M] = _mm_mul_ps(wrk->dpf[0][q*3 + p7X_M], sv);
      wrk->dpf[0][q*3 + p7X_I] = _mm_mul_ps(wrk->dpf[0][q*3 + p7X_I], sv);
    }
  XMXo(0,p7X_N) *= norm;
  XMXo(0,p7X_C) *= norm;
  XMXo(0,p7X_J) *= norm;

  xfactor =  XMXo(0,p7X_N) + XMXo(0,p7X_C) + XMXo(0,p7X_J);
  for (x = 0; x < om->abc->K; x++)
    {
      sv = _mm_setzero_ps();
      rp = om->rfv[x];
      for (q = 0; q < Q; q++)
        {
          sv = _mm_add_ps(sv, _mm_mul_ps(wrk->dpf[0][q*3 + p7X_M], *rp)); rp++;
          sv = _mm_add_ps(sv,            wrk->dpf[0][q*3 + p7X_I]);  /* insert odds implicitly 1.0 */
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
