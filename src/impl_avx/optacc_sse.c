/* Optimal accuracy alignment; SSE implementations for impl_avx dispatch.
 * These are the _sse-suffixed versions called by the runtime dispatchers
 * in optacc.c.
 */

#include "p7_config.h"

#ifdef eslENABLE_SSE

#include <float.h>

#include <xmmintrin.h>
#include <emmintrin.h>

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"


/*****************************************************************
 * 1. Standard optimal accuracy DP fill and traceback
 *****************************************************************/

/* Function:  p7_OptimalAccuracy_sse()
 *
 * Purpose:   SSE implementation of optimal accuracy DP fill.
 *            See p7_OptimalAccuracy() documentation for details.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_OptimalAccuracy_sse(const P7_OPROFILE *om, const P7_OMX *pp, P7_OMX *ox, float *ret_e)
{
  register __m128 mpv, dpv, ipv;
  register __m128 sv;
  register __m128 xEv;
  register __m128 xBv;
  register __m128 dcv;
  float  *xmx = ox->xmx;
  __m128 *dpc = ox->dpf[0];
  __m128 *dpp;
  __m128 *ppp;
  __m128 *tp;
  __m128 zerov = _mm_setzero_ps();
  __m128 infv  = _mm_set1_ps(-eslINFINITY);
  int M = om->M;
  int Q = p7O_NQF(M);
  int q;
  int j;
  int i;
  float t1, t2;

  ox->M = om->M;
  ox->L = pp->L;
  for (q = 0; q < Q; q++) MMO(dpc, q) = IMO(dpc,q) = DMO(dpc,q) = infv;
  XMXo(0, p7X_E)    = -eslINFINITY;
  XMXo(0, p7X_N)    = 0.;
  XMXo(0, p7X_J)    = -eslINFINITY;
  XMXo(0, p7X_B)    = 0.;
  XMXo(0, p7X_C)    = -eslINFINITY;

  for (i = 1; i <= pp->L; i++)
    {
      dpp = dpc;
      dpc = ox->dpf[i];
      ppp = pp->dpf[i];
      tp  = om->tfv;
      dcv = infv;
      xEv = infv;
      xBv = _mm_set1_ps(XMXo(i-1, p7X_B));

      mpv = esl_sse_rightshift_ps(MMO(dpp,Q-1), infv);
      dpv = esl_sse_rightshift_ps(DMO(dpp,Q-1), infv);
      ipv = esl_sse_rightshift_ps(IMO(dpp,Q-1), infv);
      for (q = 0; q < Q; q++)
        {
          sv  =                _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), xBv);  tp++;
          sv  = _mm_max_ps(sv, _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), mpv)); tp++;
          sv  = _mm_max_ps(sv, _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), ipv)); tp++;
          sv  = _mm_max_ps(sv, _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), dpv)); tp++;
          sv  = _mm_add_ps(sv, *ppp);                                      ppp += 2;
          xEv = _mm_max_ps(xEv, sv);

          mpv = MMO(dpp,q);
          dpv = DMO(dpp,q);
          ipv = IMO(dpp,q);

          MMO(dpc,q) = sv;
          DMO(dpc,q) = dcv;

          dcv = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), sv); tp++;

          sv         =                _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), mpv);   tp++;
          sv         = _mm_max_ps(sv, _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), ipv));  tp++;
          IMO(dpc,q) = _mm_add_ps(sv, *ppp);                                       ppp++;
        }

      dcv = esl_sse_rightshift_ps(dcv, infv);
      tp  = om->tfv + 7*Q;
      for (q = 0; q < Q; q++)
        {
          DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
          dcv         = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), DMO(dpc,q));   tp++;
        }

      for (j = 1; j < 4; j++)
        {
          dcv = esl_sse_rightshift_ps(dcv, infv);
          tp  = om->tfv + 7*Q;
          for (q = 0; q < Q; q++)
            {
              DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
              dcv         = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), dcv);   tp++;
            }
        }

      for (q = 0; q < Q; q++) xEv = _mm_max_ps(xEv, DMO(dpc,q));

      esl_sse_hmax_ps(xEv, &(XMXo(i,p7X_E)));

      t1 = ( (om->xf[p7O_J][p7O_LOOP] == 0.0) ? 0.0 : ox->xmx[(i-1)*p7X_NXCELLS+p7X_J] + pp->xmx[i*p7X_NXCELLS+p7X_J]);
      t2 = ( (om->xf[p7O_E][p7O_LOOP] == 0.0) ? 0.0 : ox->xmx[   i *p7X_NXCELLS+p7X_E]);
      ox->xmx[i*p7X_NXCELLS+p7X_J] = ESL_MAX(t1, t2);

      t1 = ( (om->xf[p7O_C][p7O_LOOP] == 0.0) ? 0.0 : ox->xmx[(i-1)*p7X_NXCELLS+p7X_C] + pp->xmx[i*p7X_NXCELLS+p7X_C]);
      t2 = ( (om->xf[p7O_E][p7O_MOVE] == 0.0) ? 0.0 : ox->xmx[   i *p7X_NXCELLS+p7X_E]);
      ox->xmx[i*p7X_NXCELLS+p7X_C] = ESL_MAX(t1, t2);

      ox->xmx[i*p7X_NXCELLS+p7X_N] = ((om->xf[p7O_N][p7O_LOOP] == 0.0) ? 0.0 : ox->xmx[(i-1)*p7X_NXCELLS+p7X_N] + pp->xmx[i*p7X_NXCELLS+p7X_N]);

      t1 = ( (om->xf[p7O_N][p7O_MOVE] == 0.0) ? 0.0 : ox->xmx[i*p7X_NXCELLS+p7X_N]);
      t2 = ( (om->xf[p7O_J][p7O_MOVE] == 0.0) ? 0.0 : ox->xmx[i*p7X_NXCELLS+p7X_J]);
      ox->xmx[i*p7X_NXCELLS+p7X_B] = ESL_MAX(t1, t2);
    }

  *ret_e = ox->xmx[pp->L*p7X_NXCELLS+p7X_C];
  return eslOK;
}


/* OA traceback helpers (standard, non-FS) */

static inline float oa_get_postprob(const P7_OMX *pp, int scur, int sprv, int k, int i);
static inline int oa_select_m(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int k);
static inline int oa_select_d(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int k);
static inline int oa_select_i(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int k);
static inline int oa_select_n(int i);
static inline int oa_select_c(const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, int i);
static inline int oa_select_j(const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, int i);
static inline int oa_select_e(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int *ret_k);
static inline int oa_select_b(const P7_OPROFILE *om,                   const P7_OMX *ox, int i);

static inline float
oa_get_postprob(const P7_OMX *pp, int scur, int sprv, int k, int i)
{
  int     Q     = p7O_NQF(pp->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  union { __m128 v; float p[4]; } u;

  switch (scur) {
  case p7T_M: u.v = MMO(pp->dpf[i], q); return u.p[r];
  case p7T_I: u.v = IMO(pp->dpf[i], q); return u.p[r];
  case p7T_N: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_N];
  case p7T_C: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_C];
  case p7T_J: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_J];
  default:    return 0.0;
  }
}

static inline int
oa_select_m(const P7_OPROFILE *om, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m128 *tp    = om->tfv + 7*q;
  __m128  xBv   = _mm_set1_ps(ox->xmx[(i-1)*p7X_NXCELLS+p7X_B]);
  __m128  mpv, dpv, ipv;
  union { __m128 v; float p[4]; } u, tv;
  float   path[4];
  int     state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };

  if (q > 0) {
    mpv = ox->dpf[i-1][(q-1)*3 + p7X_M];
    dpv = ox->dpf[i-1][(q-1)*3 + p7X_D];
    ipv = ox->dpf[i-1][(q-1)*3 + p7X_I];
  } else {
    mpv = esl_sse_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_M]);
    dpv = esl_sse_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_D]);
    ipv = esl_sse_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_I]);
  }

  u.v = xBv;  tv.v = *tp;  path[3] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);  tp++;
  u.v = mpv;  tv.v = *tp;  path[0] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);  tp++;
  u.v = ipv;  tv.v = *tp;  path[1] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);  tp++;
  u.v = dpv;  tv.v = *tp;  path[2] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
oa_select_d(const P7_OPROFILE *om, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  union { __m128 v; float p[4]; } mpv, dpv, tmdv, tddv;
  float   path[2];

  if (q > 0) {
    mpv.v  = ox->dpf[i][(q-1)*3 + p7X_M];
    dpv.v  = ox->dpf[i][(q-1)*3 + p7X_D];
    tmdv.v = om->tfv[7*(q-1) + p7O_MD];
    tddv.v = om->tfv[7*Q + (q-1)];
  } else {
    mpv.v  = esl_sse_rightshiftz_float(ox->dpf[i][(Q-1)*3 + p7X_M]);
    dpv.v  = esl_sse_rightshiftz_float(ox->dpf[i][(Q-1)*3 + p7X_D]);
    tmdv.v = esl_sse_rightshiftz_float(om->tfv[7*(Q-1) + p7O_MD]);
    tddv.v = esl_sse_rightshiftz_float(om->tfv[8*Q-1]);
  }

  path[0] = ((tmdv.p[r] == 0.0) ? -eslINFINITY : mpv.p[r]);
  path[1] = ((tddv.p[r] == 0.0) ? -eslINFINITY : dpv.p[r]);
  return  ((path[0] >= path[1]) ? p7T_M : p7T_D);
}

static inline int
oa_select_i(const P7_OPROFILE *om, const P7_OMX *ox, int i, int k)
{
  int     Q    = p7O_NQF(ox->M);
  int     q    = (k-1) % Q;
  int     r    = (k-1) / Q;
  __m128 *tp   = om->tfv + 7*q + p7O_MI;
  union { __m128 v; float p[4]; } tv, mpv, ipv;
  float   path[2];

  mpv.v = ox->dpf[i-1][q*3 + p7X_M]; tv.v = *tp;  path[0] = ((tv.p[r] == 0.0) ? -eslINFINITY : mpv.p[r]);  tp++;
  ipv.v = ox->dpf[i-1][q*3 + p7X_I]; tv.v = *tp;  path[1] = ((tv.p[r] == 0.0) ? -eslINFINITY : ipv.p[r]);
  return  ((path[0] >= path[1]) ? p7T_M : p7T_I);
}

static inline int
oa_select_n(int i)
{
  return ((i==0) ? p7T_S : p7T_N);
}

static inline int
oa_select_c(const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, int i)
{
  float path[2];
  path[0] = ( (om->xf[p7O_C][p7O_LOOP] == 0.0) ? -eslINFINITY : ox->xmx[(i-1)*p7X_NXCELLS+p7X_C] + pp->xmx[i*p7X_NXCELLS+p7X_C]);
  path[1] = ( (om->xf[p7O_E][p7O_MOVE] == 0.0) ? -eslINFINITY : ox->xmx[   i *p7X_NXCELLS+p7X_E]);
  return  ((path[0] > path[1]) ? p7T_C : p7T_E);
}

static inline int
oa_select_j(const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, int i)
{
  float path[2];
  path[0] = ( (om->xf[p7O_J][p7O_LOOP] == 0.0) ? -eslINFINITY : ox->xmx[(i-1)*p7X_NXCELLS+p7X_J] + pp->xmx[i*p7X_NXCELLS+p7X_J]);
  path[1] = ( (om->xf[p7O_E][p7O_LOOP] == 0.0) ? -eslINFINITY : ox->xmx[   i *p7X_NXCELLS+p7X_E]);
  return  ((path[0] > path[1]) ? p7T_J : p7T_E);
}

static inline int
oa_select_e(const P7_OPROFILE *om, const P7_OMX *ox, int i, int *ret_k)
{
  int     Q     = p7O_NQF(ox->M);
  __m128 *dp    = ox->dpf[i];
  union { __m128 v; float p[4]; } u;
  float  max   = -eslINFINITY;
  int    smax, kmax;
  int    q,r;

  for (q = 0; q < Q; q++)
    {
      u.v   = *dp; dp++;  for (r = 0; r < 4; r++) if (u.p[r] >= max) { max = u.p[r]; smax = p7T_M; kmax = r*Q + q + 1; }
      u.v   = *dp; dp+=2; for (r = 0; r < 4; r++) if (u.p[r] > max)  { max = u.p[r]; smax = p7T_D; kmax = r*Q + q + 1; }
    }
  *ret_k = kmax;
  return smax;
}

static inline int
oa_select_b(const P7_OPROFILE *om, const P7_OMX *ox, int i)
{
  float path[2];
  path[0] = ( (om->xf[p7O_N][p7O_MOVE] == 0.0) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_N]);
  path[1] = ( (om->xf[p7O_J][p7O_MOVE] == 0.0) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_J]);
  return  ((path[0] > path[1]) ? p7T_N : p7T_J);
}


/* Function:  p7_OATrace_sse()
 *
 * Purpose:   SSE implementation of optimal accuracy traceback.
 *            See p7_OATrace() documentation for details.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> if trace is not empty.
 */
int
p7_OATrace_sse(const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, P7_TRACE *tr)
{
  int   i   = ox->L;
  int   k   = 0;
  int   s0, s1;
  float postprob;
  int   status;

  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace not empty; needs to be Reuse()'d?");

  if ((status = p7_trace_AppendWithPP(tr, p7T_T, k, i, 0.0)) != eslOK) return status;
  if ((status = p7_trace_AppendWithPP(tr, p7T_C, k, i, 0.0)) != eslOK) return status;

  s0 = tr->st[tr->N-1];
  while (s0 != p7T_S)
    {
      switch (s0) {
      case p7T_M: s1 = oa_select_m(om,     ox, i, k);  k--; i--; break;
      case p7T_D: s1 = oa_select_d(om,     ox, i, k);  k--;      break;
      case p7T_I: s1 = oa_select_i(om,     ox, i, k);       i--; break;
      case p7T_N: s1 = oa_select_n(i);                           break;
      case p7T_C: s1 = oa_select_c(om, pp, ox, i);               break;
      case p7T_J: s1 = oa_select_j(om, pp, ox, i);               break;
      case p7T_E: s1 = oa_select_e(om,     ox, i, &k);           break;
      case p7T_B: s1 = oa_select_b(om,     ox, i);               break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }
      if (s1 == -1) ESL_EXCEPTION(eslEINVAL, "OA traceback choice failed");

      postprob = oa_get_postprob(pp, s1, s0, k, i);
      if ((status = p7_trace_AppendWithPP(tr, s1, k, i, postprob)) != eslOK) return status;

      if ( (s1 == p7T_N || s1 == p7T_J || s1 == p7T_C) && s1 == s0) i--;
      s0 = s1;
    }
  tr->M = om->M;
  tr->L = ox->L;
  return p7_trace_Reverse(tr);
}
#endif /* eslENABLE_SSE */
