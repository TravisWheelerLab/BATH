/* vitscore_avx.c — Viterbi score (one-row, no traceback); AVX2 implementation.
 * Ported from vitscore_sse.c using 256-bit vectors.
 *
 * Contents:
 *    1. p7_ViterbiScore_avx()
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


/*****************************************************************
 * Local helpers (float, AVX2)
 *****************************************************************/

/* Shift __m256 right by 1 float element, inserting infv[0] at position 0. */
static inline __m256
avx_rightshift_ps(__m256 v, __m256 infv)
{
  return _mm256_blend_ps(esl_avx_rightshiftz_float(v), infv, 0x01);
}

/* Horizontal max of 8 floats in a __m256, returned by pointer. */
static inline void
avx_hmax_ps(__m256 v, float *ret_max)
{
  __m128 hi   = _mm256_extractf128_ps(v, 1);
  __m128 lo   = _mm256_castps256_ps128(v);
  __m128 max4 = _mm_max_ps(lo, hi);
  max4 = _mm_max_ps(max4, _mm_shuffle_ps(max4, max4, _MM_SHUFFLE(2,3,0,1)));
  max4 = _mm_max_ps(max4, _mm_shuffle_ps(max4, max4, _MM_SHUFFLE(1,0,3,2)));
  _mm_store_ss(ret_max, max4);
}

/* Return nonzero if any a[i] > b[i] (float, 8 lanes). */
static inline int
avx_any_gt_ps(__m256 a, __m256 b)
{
  return (_mm256_movemask_ps(_mm256_cmp_ps(a, b, _CMP_GT_OS)) != 0);
}


/*****************************************************************
 * 1. p7_ViterbiScore_avx()
 *****************************************************************/

/* Function:  p7_ViterbiScore_avx()
 *
 * Purpose:   AVX2 implementation of one-row Viterbi score (no traceback).
 *            Uses 256-bit vectors: Q = ceil(M/8) stripes of 8 float lanes.
 *            See p7_ViterbiScore() documentation for details.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_ViterbiScore_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register __m256 mpv, dpv, ipv;
  register __m256 sv;
  register __m256 dcv;
  register __m256 xEv;
  register __m256 xBv;
  register __m256 Dmaxv;
  __m256   infv;
  float    xN, xE, xB, xC, xJ;
  float    Dmax;
  int i;
  int q;
  int Q       = p7O_NQF_AVX(om->M);
  __m256 *dp  = ox->dpf_avx[0];   /* one-row: always dpf_avx[0] */
  __m256 *rsc;
  __m256 *tsc;

  if (Q > ox->allocQ4_avx) ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M = om->M;

  infv = _mm256_set1_ps(-eslINFINITY);
  for (q = 0; q < Q; q++)
    MMO(dp,q) = IMO(dp,q) = DMO(dp,q) = infv;
  xN = 0.;
  xB = om->xf[p7O_N][p7O_MOVE];
  xE = -eslINFINITY;
  xJ = -eslINFINITY;
  xC = -eslINFINITY;

  for (i = 1; i <= L; i++)
    {
      rsc   = om->rfv_avx[dsq[i]];
      tsc   = om->tfv_avx;
      dcv   = infv;
      xEv   = infv;
      Dmaxv = infv;
      xBv   = _mm256_set1_ps(xB);

      mpv = avx_rightshift_ps(MMO(dp,Q-1), infv);
      dpv = avx_rightshift_ps(DMO(dp,Q-1), infv);
      ipv = avx_rightshift_ps(IMO(dp,Q-1), infv);

      for (q = 0; q < Q; q++)
        {
          sv   =                   _mm256_add_ps(xBv, *tsc);  tsc++;   /* B->Mk     */
          sv   = _mm256_max_ps(sv, _mm256_add_ps(mpv, *tsc)); tsc++;   /* Mk-1->Mk  */
          sv   = _mm256_max_ps(sv, _mm256_add_ps(ipv, *tsc)); tsc++;   /* Ik-1->Mk  */
          sv   = _mm256_max_ps(sv, _mm256_add_ps(dpv, *tsc)); tsc++;   /* Dk-1->Mk  */
          sv   = _mm256_add_ps(sv, *rsc);                     rsc++;   /* + emission */
          xEv  = _mm256_max_ps(xEv, sv);

          mpv = MMO(dp,q);
          dpv = DMO(dp,q);
          ipv = IMO(dp,q);

          MMO(dp,q) = sv;
          DMO(dp,q) = dcv;

          dcv   = _mm256_add_ps(sv, *tsc); tsc++;                      /* Mk->Dk+1  */
          Dmaxv = _mm256_max_ps(dcv, Dmaxv);

          sv        =                   _mm256_add_ps(mpv, *tsc);  tsc++;  /* Mk->Ik  */
          sv        = _mm256_max_ps(sv, _mm256_add_ps(ipv, *tsc)); tsc++;  /* Ik->Ik  */
          IMO(dp,q) = _mm256_add_ps(sv, *rsc);                     rsc++;  /* + emission */
        }

      avx_hmax_ps(xEv, &xE);
      xN = xN +  om->xf[p7O_N][p7O_LOOP];
      xC = ESL_MAX(xC + om->xf[p7O_C][p7O_LOOP],  xE + om->xf[p7O_E][p7O_MOVE]);
      xJ = ESL_MAX(xJ + om->xf[p7O_J][p7O_LOOP],  xE + om->xf[p7O_E][p7O_LOOP]);
      xB = ESL_MAX(xJ + om->xf[p7O_J][p7O_MOVE],  xN + om->xf[p7O_N][p7O_MOVE]);

      /* Lazy F loop. */
      avx_hmax_ps(Dmaxv, &Dmax);
      if (Dmax + om->ddbound_w > xB)
        {
          dcv = avx_rightshift_ps(dcv, infv);
          tsc = om->tfv_avx + 7*Q;
          for (q = 0; q < Q; q++)
            {
              DMO(dp,q) = _mm256_max_ps(dcv, DMO(dp,q));
              dcv       = _mm256_add_ps(DMO(dp,q), *tsc); tsc++;
            }

          do {
            dcv = avx_rightshift_ps(dcv, infv);
            tsc = om->tfv_avx + 7*Q;
            for (q = 0; q < Q; q++)
              {
                if (! avx_any_gt_ps(dcv, DMO(dp,q))) break;
                DMO(dp,q) = _mm256_max_ps(dcv, DMO(dp,q));
                dcv       = _mm256_add_ps(DMO(dp,q), *tsc); tsc++;
              }
          } while (q == Q);
        }
      else
        {
          dcv       = avx_rightshift_ps(dcv, infv);
          DMO(dp,0) = dcv;
        }

    } /* end loop over sequence residues 1..L */

  *ret_sc = xC + om->xf[p7O_C][p7O_MOVE];
  return eslOK;
}

#endif /* eslENABLE_AVX */
