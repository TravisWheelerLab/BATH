/* Viterbi score (one-row, no traceback); SSE implementation for impl_avx dispatch.
 * This is the _sse-suffixed version called by the runtime dispatcher in vitscore.c.
 */

#ifdef eslENABLE_SSE
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>
#include <emmintrin.h>

#include "easel.h"
#include "esl_sse.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_ViterbiScore_sse()
 *
 * Purpose:   SSE implementation of one-row Viterbi score (no traceback).
 *            See p7_ViterbiScore() documentation for details.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_ViterbiScore_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register __m128 mpv, dpv, ipv;
  register __m128 sv;
  register __m128 dcv;
  register __m128 xEv;
  register __m128 xBv;
  register __m128 Dmaxv;
  __m128  infv;
  float    xN, xE, xB, xC, xJ;
  float    Dmax;
  int i;
  int q;
  int Q      = p7O_NQF(om->M);
  __m128 *dp = ox->dpf[0];   /* one-row: always dpf[0] */
  __m128 *rsc;
  __m128 *tsc;

  if (Q > ox->allocQ4) ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  ox->M = om->M;

  infv = _mm_set1_ps(-eslINFINITY);
  for (q = 0; q < Q; q++)
    MMO(dp,q) = IMO(dp,q) = DMO(dp,q) = infv;
  xN = 0.;
  xB = om->xf[p7O_N][p7O_MOVE];
  xE = -eslINFINITY;
  xJ = -eslINFINITY;
  xC = -eslINFINITY;

  for (i = 1; i <= L; i++)
    {
      rsc   = om->rf[dsq[i]];
      tsc   = om->tf;
      dcv   = infv;
      xEv   = infv;
      Dmaxv = infv;
      xBv   = _mm_set1_ps(xB);

      mpv = esl_sse_rightshift_ps(MMO(dp,Q-1), infv);
      dpv = esl_sse_rightshift_ps(DMO(dp,Q-1), infv);
      ipv = esl_sse_rightshift_ps(IMO(dp,Q-1), infv);

      for (q = 0; q < Q; q++)
        {
          sv   =                _mm_add_ps(xBv, *tsc);  tsc++;
          sv   = _mm_max_ps(sv, _mm_add_ps(mpv, *tsc)); tsc++;
          sv   = _mm_max_ps(sv, _mm_add_ps(ipv, *tsc)); tsc++;
          sv   = _mm_max_ps(sv, _mm_add_ps(dpv, *tsc)); tsc++;
          sv   = _mm_add_ps(sv, *rsc);                  rsc++;
          xEv  = _mm_max_ps(xEv, sv);

          mpv = MMO(dp,q);
          dpv = DMO(dp,q);
          ipv = IMO(dp,q);

          MMO(dp,q) = sv;
          DMO(dp,q) = dcv;

          dcv   = _mm_add_ps(sv, *tsc); tsc++;
          Dmaxv = _mm_max_ps(dcv, Dmaxv);

          sv        =                _mm_add_ps(mpv, *tsc);  tsc++;
          sv        = _mm_max_ps(sv, _mm_add_ps(ipv, *tsc)); tsc++;
          IMO(dp,q) = _mm_add_ps(sv, *rsc);                  rsc++;
        }

      esl_sse_hmax_ps(xEv, &xE);
      xN = xN +  om->xf[p7O_N][p7O_LOOP];
      xC = ESL_MAX(xC + om->xf[p7O_C][p7O_LOOP],  xE + om->xf[p7O_E][p7O_MOVE]);
      xJ = ESL_MAX(xJ + om->xf[p7O_J][p7O_LOOP],  xE + om->xf[p7O_E][p7O_LOOP]);
      xB = ESL_MAX(xJ + om->xf[p7O_J][p7O_MOVE],  xN + om->xf[p7O_N][p7O_MOVE]);

      /* Lazy F loop. */
      esl_sse_hmax_ps(Dmaxv, &Dmax);
      if (Dmax + om->ddbound_f > xB)
        {
          dcv = esl_sse_rightshift_ps(dcv, infv);
          tsc = om->tf + 7*Q;
          for (q = 0; q < Q; q++)
            {
              DMO(dp,q) = _mm_max_ps(dcv, DMO(dp,q));
              dcv       = _mm_add_ps(DMO(dp,q), *tsc); tsc++;
            }

          do {
            dcv = esl_sse_rightshift_ps(dcv, infv);
            tsc = om->tf + 7*Q;
            for (q = 0; q < Q; q++)
              {
                if (! esl_sse_any_gt_ps(dcv, DMO(dp,q))) break;
                DMO(dp,q) = _mm_max_ps(dcv, DMO(dp,q));
                dcv       = _mm_add_ps(DMO(dp,q), *tsc); tsc++;
              }
          } while (q == Q);
        }
      else
        {
          dcv       = esl_sse_rightshift_ps(dcv, infv);
          DMO(dp,0) = dcv;
        }

    } /* end loop over sequence residues 1..L */

  *ret_sc = xC + om->xf[p7O_C][p7O_MOVE];
  return eslOK;
}

#endif /* eslENABLE_SSE */
