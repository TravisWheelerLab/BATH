/* Viterbi full-matrix and traceback; SSE implementations for impl_avx dispatch.
 * These are the _sse-suffixed versions called by the runtime dispatcher in viterbi.c.
 */

#ifdef eslENABLE_SSE
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>
#include <emmintrin.h>

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Viterbi_sse()
 *
 * Purpose:   SSE implementation of full-matrix Viterbi, enabling traceback.
 *            See p7_Viterbi() documentation for details.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_Viterbi_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register __m128 mpv, dpv, ipv;
  register __m128 sv;
  register __m128 dcv;
  register __m128 xEv;
  register __m128 xBv;
  __m128  infv;
  float    xN, xE, xB, xC, xJ;
  int i;
  int q;
  int Q       = p7O_NQF(om->M);
  __m128 *dpp;
  __m128 *dp;
  __m128 *rsc;
  __m128 *tsc;
  float  *xmx = ox->xmx;

  if (Q > ox->allocQ4)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (Q)");
  if (L >= ox->validR)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (rows)");
  if (ox->xmx == NULL)  ESL_EXCEPTION(eslEINVAL, "DP matrix has no xmx special state storage");
  ox->M = om->M;
  ox->L = L;

  infv = _mm_set1_ps(-eslINFINITY);
  dp   = ox->dpf[0];
  for (q = 0; q < Q; q++)
    MMO(dp,q) = IMO(dp,q) = DMO(dp,q) = infv;

  xN = 0.0f;
  xB = om->xf[p7O_N][p7O_MOVE];
  xE = -eslINFINITY;
  xJ = -eslINFINITY;
  xC = -eslINFINITY;

  XMXo(0, p7X_N) = xN;
  XMXo(0, p7X_B) = xB;
  XMXo(0, p7X_E) = xE;
  XMXo(0, p7X_J) = xJ;
  XMXo(0, p7X_C) = xC;

  for (i = 1; i <= L; i++)
    {
      rsc  = om->rfv[dsq[i]];
      tsc  = om->tfv;
      dcv  = infv;
      xEv  = infv;
      xBv  = _mm_set1_ps(xB);

      dpp  = ox->dpf[i-1];
      dp   = ox->dpf[i];

      mpv = esl_sse_rightshift_ps(MMO(dpp, Q-1), infv);
      dpv = esl_sse_rightshift_ps(DMO(dpp, Q-1), infv);
      ipv = esl_sse_rightshift_ps(IMO(dpp, Q-1), infv);

      for (q = 0; q < Q; q++)
        {
          sv   =                _mm_add_ps(xBv, *tsc);  tsc++;  /* B->Mk  */
          sv   = _mm_max_ps(sv, _mm_add_ps(mpv, *tsc)); tsc++;  /* Mk-1->Mk */
          sv   = _mm_max_ps(sv, _mm_add_ps(ipv, *tsc)); tsc++;  /* Ik-1->Mk */
          sv   = _mm_max_ps(sv, _mm_add_ps(dpv, *tsc)); tsc++;  /* Dk-1->Mk */
          sv   = _mm_add_ps(sv, *rsc);                  rsc++;  /* + emission */
          xEv  = _mm_max_ps(xEv, sv);

          mpv = MMO(dpp, q);
          dpv = DMO(dpp, q);
          ipv = IMO(dpp, q);

          MMO(dp, q) = sv;
          DMO(dp, q) = dcv;

          dcv = _mm_add_ps(sv, *tsc); tsc++;  /* Mk->Dk+1 */

          sv        =                _mm_add_ps(mpv, *tsc);  tsc++;  /* Mk->Ik   */
          sv        = _mm_max_ps(sv, _mm_add_ps(ipv, *tsc)); tsc++;  /* Ik->Ik   */
          IMO(dp,q) = sv;
        }

      esl_sse_hmax_ps(xEv, &xE);
      xN = xN +  om->xf[p7O_N][p7O_LOOP];
      xC = ESL_MAX(xC + om->xf[p7O_C][p7O_LOOP],  xE + om->xf[p7O_E][p7O_MOVE]);
      xJ = ESL_MAX(xJ + om->xf[p7O_J][p7O_LOOP],  xE + om->xf[p7O_E][p7O_LOOP]);
      xB = ESL_MAX(xJ + om->xf[p7O_J][p7O_MOVE],  xN + om->xf[p7O_N][p7O_MOVE]);

      XMXo(i, p7X_E) = xE;
      XMXo(i, p7X_N) = xN;
      XMXo(i, p7X_J) = xJ;
      XMXo(i, p7X_B) = xB;
      XMXo(i, p7X_C) = xC;

      /* D->D propagation. */
      dcv = esl_sse_rightshift_ps(dcv, infv);
      tsc = om->tfv + 7*Q;
      for (q = 0; q < Q; q++)
        {
          DMO(dp,q) = _mm_max_ps(dcv, DMO(dp,q));
          dcv       = _mm_add_ps(DMO(dp,q), *tsc); tsc++;
        }

      do {
        dcv = esl_sse_rightshift_ps(dcv, infv);
        tsc = om->tfv + 7*Q;
        for (q = 0; q < Q; q++)
          {
            if (! esl_sse_any_gt_ps(dcv, DMO(dp,q))) break;
            DMO(dp,q) = _mm_max_ps(dcv, DMO(dp,q));
            dcv       = _mm_add_ps(DMO(dp,q), *tsc); tsc++;
          }
      } while (q == Q);

    } /* end loop over sequence residues 1..L */

  if (ret_sc != NULL) *ret_sc = xC + om->xf[p7O_C][p7O_MOVE];
  return eslOK;
}


/* Function:  p7_Viterbi_Trace_sse()
 *
 * Purpose:   SSE implementation of Viterbi traceback from an lspace full matrix.
 *            See p7_Viterbi_Trace() documentation for details.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if an impossible state is reached.
 */
int
p7_Viterbi_Trace_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *ox,
                     P7_TRACE *tr)
{
  int   Q     = p7O_NQF(om->M);
  int   M     = om->M;
  int   i     = L;
  int   k     = 0;
  float r_tol = 1e-5;
  float a_tol = 1e-4;
  int   sprv, scur;
  int   status;

#define OMMo(i,k)  ((k)<1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_M, (i), (k)))
#define ODMo(i,k)  ((k)<1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_D, (i), (k)))
#define OIMo(i,k)  ((k)<1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_I, (i), (k)))
#define OXMXo(i,s)  (ox->xmx[(i)*p7X_NXCELLS+(s)])

  if ((status = p7_trace_Append(tr, p7T_T, k, i)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_C, k, i)) != eslOK) return status;
  sprv = p7T_C;

  while (sprv != p7T_S) {
    switch (sprv) {

    case p7T_C:
      if (OXMXo(i, p7X_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
      if      (esl_FCompare(OXMXo(i,p7X_C), OXMXo(i-1,p7X_C) + om->xf[p7O_C][p7O_LOOP], r_tol, a_tol) == eslOK) scur = p7T_C;
      else if (esl_FCompare(OXMXo(i,p7X_C), OXMXo(i,  p7X_E) + om->xf[p7O_E][p7O_MOVE], r_tol, a_tol) == eslOK) scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7T_E:
      if (OXMXo(i, p7X_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
      if (p7_oprofile_IsLocal(om)) {
        scur = p7T_M;
        for (k = M; k >= 1; k--) if (esl_FCompare(OXMXo(i,p7X_E), OMMo(i,k), r_tol, a_tol) == eslOK) break;
        if (k == 0) ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      } else {
        if      (esl_FCompare(OXMXo(i,p7X_E), OMMo(i,M), r_tol, a_tol) == eslOK) { scur = p7T_M; k = M; }
        else if (esl_FCompare(OXMXo(i,p7X_E), ODMo(i,M), r_tol, a_tol) == eslOK) { scur = p7T_D; k = M; }
        else ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      }
      break;

    case p7T_M: {
      union { __m128 v; float p[4]; } u;
      int   q   = (k-1) % Q,  r   = (k-1) / Q;
      float rsc = (u.v = om->rfv[dsq[i]][q],        u.p[r]);
      float tBM = (u.v = om->tfv[7*q + p7O_BM],     u.p[r]);
      float tMM = (u.v = om->tfv[7*q + p7O_MM],     u.p[r]);
      float tIM = (u.v = om->tfv[7*q + p7O_IM],     u.p[r]);
      float tDM = (u.v = om->tfv[7*q + p7O_DM],     u.p[r]);
      if (OMMo(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k, i);
      if      (esl_FCompare(OMMo(i,k), OXMXo(i-1,p7X_B) + tBM + rsc, r_tol, a_tol) == eslOK) scur = p7T_B;
      else if (esl_FCompare(OMMo(i,k), OMMo(i-1,k-1)   + tMM + rsc, r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(OMMo(i,k), OIMo(i-1,k-1)   + tIM + rsc, r_tol, a_tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare(OMMo(i,k), ODMo(i-1,k-1)   + tDM + rsc, r_tol, a_tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k, i);
      k--; i--;
      break; }

    case p7T_D: {
      union { __m128 v; float p[4]; } u;
      int   qp  = (k-2) % Q, rp  = (k-2) / Q;
      float tMD = (k > 1) ? (u.v = om->tfv[7*qp + p7O_MD], u.p[rp]) : -eslINFINITY;
      float tDD = (k > 1) ? (u.v = om->tfv[7*Q  + qp],     u.p[rp]) : -eslINFINITY;
      if (ODMo(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k, i);
      if      (esl_FCompare(ODMo(i,k), OMMo(i,k-1) + tMD, r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(ODMo(i,k), ODMo(i,k-1) + tDD, r_tol, a_tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k, i);
      k--;
      break; }

    case p7T_I: {
      union { __m128 v; float p[4]; } u;
      int   q   = (k-1) % Q, r   = (k-1) / Q;
      float tMI = (u.v = om->tfv[7*q + p7O_MI], u.p[r]);
      float tII = (u.v = om->tfv[7*q + p7O_II], u.p[r]);
      if (OIMo(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k, i);
      if      (esl_FCompare(OIMo(i,k), OMMo(i-1,k) + tMI, r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(OIMo(i,k), OIMo(i-1,k) + tII, r_tol, a_tol) == eslOK) scur = p7T_I;
      else ESL_EXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k, i);
      i--;
      break; }

    case p7T_N:
      if (OXMXo(i, p7X_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = (i == 0) ? p7T_S : p7T_N;
      break;

    case p7T_B:
      if (OXMXo(i, p7X_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);
      if      (esl_FCompare(OXMXo(i,p7X_B), OXMXo(i,p7X_N) + om->xf[p7O_N][p7O_MOVE], r_tol, a_tol) == eslOK) scur = p7T_N;
      else if (esl_FCompare(OXMXo(i,p7X_B), OXMXo(i,p7X_J) + om->xf[p7O_J][p7O_MOVE], r_tol, a_tol) == eslOK) scur = p7T_J;
      else ESL_EXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    case p7T_J:
      if (OXMXo(i, p7X_J) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible J reached at i=%d", i);
      if      (esl_FCompare(OXMXo(i,p7X_J), OXMXo(i-1,p7X_J) + om->xf[p7O_J][p7O_LOOP], r_tol, a_tol) == eslOK) scur = p7T_J;
      else if (esl_FCompare(OXMXo(i,p7X_J), OXMXo(i,  p7X_E) + om->xf[p7O_E][p7O_LOOP], r_tol, a_tol) == eslOK) scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "J at i=%d couldn't be traced", i);
      break;

    default: ESL_EXCEPTION(eslEINVAL, "bogus state in Viterbi traceback");
    }

    if ((status = p7_trace_Append(tr, scur, k, i)) != eslOK) return status;
    if ((scur == p7T_N || scur == p7T_J || scur == p7T_C) && scur == sprv) i--;
    sprv = scur;
  }

  tr->M = M;
  tr->L = L;

#undef OMMo
#undef ODMo
#undef OIMo
#undef OXMXo

  return p7_trace_Reverse(tr);
}

#endif /* eslENABLE_SSE */
