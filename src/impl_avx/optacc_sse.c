/* Optimal accuracy alignment; SSE implementations for impl_avx dispatch.
 * These are the _sse-suffixed versions called by the runtime dispatchers
 * in optacc.c and optacc_fs.c.
 */

#ifdef eslENABLE_SSE
#include "p7_config.h"

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

/*****************************************************************
 * 2. Frameshift optimal accuracy DP fill and traceback
 *****************************************************************/

/* Function:  p7_OptimalAccuracy_Frameshift_sse()
 *
 * Purpose:   SSE implementation of frameshift optimal accuracy DP fill.
 *            See p7_OptimalAccuracy_Frameshift() documentation for details.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_OptimalAccuracy_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                                   P7_OMX *ox, float *ret_e)
{
  int    M   = om_fs->M;
  int    L   = pp->L;
  int    Q   = p7O_NQF(M);
  int    i, q, j;
  int    qM  = (M-1) % Q;
  int    rM  = (M-1) / Q;
  float  xN, xE, xB, xJ, xC;
  float  t1, t2;
  register __m128 mpv1, dpv1, ipv1;
  register __m128 mpv2, dpv2, ipv2;
  register __m128 mpv3, dpv3, ipv3;
  register __m128 mpv4, dpv4, ipv4;
  register __m128 mpv5, dpv5, ipv5;
  register __m128 sv1, sv2, sv3, sv4, sv5, sv;
  register __m128 xEv, xBv1, xBv2, xBv3, xBv4, xBv5;
  register __m128 dcv;
  __m128 bm, mm, im, dm, md, mi, ii;
  __m128 *tp;
  __m128 *dpc;
  __m128 *dp1, *dp2, *dp3, *dp4, *dp5;
  __m128 *ppc;
  __m128  zerov = _mm_setzero_ps();
  __m128  infv  = _mm_set1_ps(-eslINFINITY);
  union { __m128 v; float p[4]; } imfix;

  ox->M = M;
  ox->L = L;

  dpc = ox->dpf[0];
  for (q = 0; q < Q; q++) { MMO(dpc,q) = DMO(dpc,q) = IMO(dpc,q) = infv; }
  ox->xmx[0*p7X_NXCELLS+p7X_E] = -eslINFINITY;
  ox->xmx[0*p7X_NXCELLS+p7X_N] =  0.0f;
  ox->xmx[0*p7X_NXCELLS+p7X_J] = -eslINFINITY;
  ox->xmx[0*p7X_NXCELLS+p7X_B] =  0.0f;
  ox->xmx[0*p7X_NXCELLS+p7X_C] = -eslINFINITY;

  for (i = 1; i <= L; i++)
    {
      dpc = ox->dpf[i];
      ppc = pp->dpf[i];

      dp1 = ox->dpf[i-1];
      dp2 = (i >= 2) ? ox->dpf[i-2] : ox->dpf[0];
      dp3 = (i >= 3) ? ox->dpf[i-3] : ox->dpf[0];
      dp4 = (i >= 4) ? ox->dpf[i-4] : ox->dpf[0];
      dp5 = (i >= 5) ? ox->dpf[i-5] : ox->dpf[0];

      xBv1 = _mm_set1_ps(ox->xmx[(i-1)*p7X_NXCELLS + p7X_B]);
      xBv2 = (i >= 2) ? _mm_set1_ps(ox->xmx[(i-2)*p7X_NXCELLS + p7X_B]) : infv;
      xBv3 = (i >= 3) ? _mm_set1_ps(ox->xmx[(i-3)*p7X_NXCELLS + p7X_B]) : infv;
      xBv4 = (i >= 4) ? _mm_set1_ps(ox->xmx[(i-4)*p7X_NXCELLS + p7X_B]) : infv;
      xBv5 = (i >= 5) ? _mm_set1_ps(ox->xmx[(i-5)*p7X_NXCELLS + p7X_B]) : infv;

      mpv1 = esl_sse_rightshift_ps(MMO(dp1, Q-1), infv);
      dpv1 = esl_sse_rightshift_ps(DMO(dp1, Q-1), infv);
      ipv1 = esl_sse_rightshift_ps(IMO(dp1, Q-1), infv);
      mpv2 = esl_sse_rightshift_ps(MMO(dp2, Q-1), infv);
      dpv2 = esl_sse_rightshift_ps(DMO(dp2, Q-1), infv);
      ipv2 = esl_sse_rightshift_ps(IMO(dp2, Q-1), infv);
      mpv3 = esl_sse_rightshift_ps(MMO(dp3, Q-1), infv);
      dpv3 = esl_sse_rightshift_ps(DMO(dp3, Q-1), infv);
      ipv3 = esl_sse_rightshift_ps(IMO(dp3, Q-1), infv);
      mpv4 = esl_sse_rightshift_ps(MMO(dp4, Q-1), infv);
      dpv4 = esl_sse_rightshift_ps(DMO(dp4, Q-1), infv);
      ipv4 = esl_sse_rightshift_ps(IMO(dp4, Q-1), infv);
      mpv5 = esl_sse_rightshift_ps(MMO(dp5, Q-1), infv);
      dpv5 = esl_sse_rightshift_ps(DMO(dp5, Q-1), infv);
      ipv5 = esl_sse_rightshift_ps(IMO(dp5, Q-1), infv);

      dcv = infv;
      xEv = infv;
      tp  = om_fs->tfv;

      for (q = 0; q < Q; q++)
        {
          bm = _mm_cmpgt_ps(*tp, zerov);  tp++;
          mm = _mm_cmpgt_ps(*tp, zerov);  tp++;
          im = _mm_cmpgt_ps(*tp, zerov);  tp++;
          dm = _mm_cmpgt_ps(*tp, zerov);  tp++;
          md = _mm_cmpgt_ps(*tp, zerov);  tp++;
          mi = _mm_cmpgt_ps(*tp, zerov);  tp++;
          ii = _mm_cmpgt_ps(*tp, zerov);  tp++;

          sv1 =                _mm_and_ps(bm, xBv1);
          sv1 = _mm_max_ps(sv1, _mm_and_ps(mm, mpv1));
          sv1 = _mm_max_ps(sv1, _mm_and_ps(im, ipv1));
          sv1 = _mm_max_ps(sv1, _mm_and_ps(dm, dpv1));
          sv1 = _mm_add_ps(sv1, MMO_FS(ppc, q, p7X_FS_C1));

          sv2 =                _mm_and_ps(bm, xBv2);
          sv2 = _mm_max_ps(sv2, _mm_and_ps(mm, mpv2));
          sv2 = _mm_max_ps(sv2, _mm_and_ps(im, ipv2));
          sv2 = _mm_max_ps(sv2, _mm_and_ps(dm, dpv2));
          sv2 = _mm_add_ps(sv2, MMO_FS(ppc, q, p7X_FS_C2));

          sv3 =                _mm_and_ps(bm, xBv3);
          sv3 = _mm_max_ps(sv3, _mm_and_ps(mm, mpv3));
          sv3 = _mm_max_ps(sv3, _mm_and_ps(im, ipv3));
          sv3 = _mm_max_ps(sv3, _mm_and_ps(dm, dpv3));
          sv3 = _mm_add_ps(sv3, MMO_FS(ppc, q, p7X_FS_C3));

          sv4 =                _mm_and_ps(bm, xBv4);
          sv4 = _mm_max_ps(sv4, _mm_and_ps(mm, mpv4));
          sv4 = _mm_max_ps(sv4, _mm_and_ps(im, ipv4));
          sv4 = _mm_max_ps(sv4, _mm_and_ps(dm, dpv4));
          sv4 = _mm_add_ps(sv4, MMO_FS(ppc, q, p7X_FS_C4));

          sv5 =                _mm_and_ps(bm, xBv5);
          sv5 = _mm_max_ps(sv5, _mm_and_ps(mm, mpv5));
          sv5 = _mm_max_ps(sv5, _mm_and_ps(im, ipv5));
          sv5 = _mm_max_ps(sv5, _mm_and_ps(dm, dpv5));
          sv5 = _mm_add_ps(sv5, MMO_FS(ppc, q, p7X_FS_C5));

          sv = _mm_max_ps(_mm_max_ps(sv1, sv2),
                          _mm_max_ps(_mm_max_ps(sv3, sv4), sv5));
          xEv = _mm_max_ps(xEv, sv);

          mpv1 = MMO(dp1, q);  dpv1 = DMO(dp1, q);  ipv1 = IMO(dp1, q);
          mpv2 = MMO(dp2, q);  dpv2 = DMO(dp2, q);  ipv2 = IMO(dp2, q);
          mpv3 = MMO(dp3, q);  dpv3 = DMO(dp3, q);  ipv3 = IMO(dp3, q);
          mpv4 = MMO(dp4, q);  dpv4 = DMO(dp4, q);  ipv4 = IMO(dp4, q);
          mpv5 = MMO(dp5, q);  dpv5 = DMO(dp5, q);  ipv5 = IMO(dp5, q);

          MMO(dpc, q) = sv;
          DMO(dpc, q) = dcv;
          dcv = _mm_and_ps(md, sv);

          sv  =                _mm_and_ps(mi, mpv3);
          sv  = _mm_max_ps(sv, _mm_and_ps(ii, ipv3));
          IMO(dpc, q) = _mm_add_ps(sv, IMO_FS(ppc, q));
        }

      imfix.v      = IMO(dpc, qM);
      imfix.p[rM]  = -eslINFINITY;
      IMO(dpc, qM) = imfix.v;

      dcv = esl_sse_rightshift_ps(dcv, infv);
      tp  = om_fs->tfv + 7*Q;
      for (q = 0; q < Q; q++)
        {
          DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
          dcv         = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), DMO(dpc, q));  tp++;
        }

      for (j = 1; j < 4; j++)
        {
          dcv = esl_sse_rightshift_ps(dcv, infv);
          tp  = om_fs->tfv + 7*Q;
          for (q = 0; q < Q; q++)
            {
              DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
              dcv         = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), dcv);  tp++;
            }
        }

      for (q = 0; q < Q; q++) xEv = _mm_max_ps(xEv, DMO(dpc, q));

      esl_sse_hmax_ps(xEv, &xE);
      ox->xmx[i*p7X_NXCELLS+p7X_E] = xE;

      if (i > 2)
        xN = (om_fs->xf[p7O_N][p7O_LOOP] == 0.0f) ? 0.0f :
             ox->xmx[(i-3)*p7X_NXCELLS+p7X_N] + pp->xmx[i*p7X_NXCELLS+p7X_N];
      else
        xN = (om_fs->xf[p7O_N][p7O_LOOP] == 0.0f) ? 0.0f :
             pp->xmx[i*p7X_NXCELLS+p7X_N];
      ox->xmx[i*p7X_NXCELLS+p7X_N] = xN;

      if (i > 2) {
        t1 = (om_fs->xf[p7O_J][p7O_LOOP] == 0.0f) ? 0.0f :
             ox->xmx[(i-3)*p7X_NXCELLS+p7X_J] + pp->xmx[i*p7X_NXCELLS+p7X_J];
        t2 = (om_fs->xf[p7O_E][p7O_LOOP] == 0.0f) ? 0.0f : xE;
        xJ = ESL_MAX(t1, t2);
      } else {
        xJ = (om_fs->xf[p7O_E][p7O_LOOP] == 0.0f) ? 0.0f : xE;
      }
      ox->xmx[i*p7X_NXCELLS+p7X_J] = xJ;

      if (i > 2) {
        t1 = (om_fs->xf[p7O_C][p7O_LOOP] == 0.0f) ? 0.0f :
             ox->xmx[(i-3)*p7X_NXCELLS+p7X_C] + pp->xmx[i*p7X_NXCELLS+p7X_C];
        t2 = (om_fs->xf[p7O_E][p7O_MOVE] == 0.0f) ? 0.0f : xE;
        xC = ESL_MAX(t1, t2);
      } else {
        xC = (om_fs->xf[p7O_E][p7O_MOVE] == 0.0f) ? 0.0f : xE;
      }
      ox->xmx[i*p7X_NXCELLS+p7X_C] = xC;

      t1 = (om_fs->xf[p7O_N][p7O_MOVE] == 0.0f) ? 0.0f : xN;
      t2 = (om_fs->xf[p7O_J][p7O_MOVE] == 0.0f) ? 0.0f : xJ;
      ox->xmx[i*p7X_NXCELLS+p7X_B] = ESL_MAX(t1, t2);
    }

  *ret_e = ox->xmx[L    *p7X_NXCELLS+p7X_C] +
           ox->xmx[(L-1)*p7X_NXCELLS+p7X_C] +
           ox->xmx[(L-2)*p7X_NXCELLS+p7X_C];
  return eslOK;
}


/* OA traceback helpers (FS variants) */

static inline float oa_get_postprob_fs(const P7_OMX *pp, int scur, int sprv, int k, int i);
static inline int oa_select_m_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k);
static inline int oa_select_d_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k);
static inline int oa_select_i_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k);
static inline int oa_select_n_fs(int i);
static inline int oa_select_c_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp, const P7_OMX *ox, int i);
static inline int oa_select_j_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp, const P7_OMX *ox, int i);
static inline int oa_select_e_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int *ret_k);
static inline int oa_select_b_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i);
static inline int oa_select_codon_fs(const P7_OMX *pp, int i, int k);

static inline float
oa_get_postprob_fs(const P7_OMX *pp, int scur, int sprv, int k, int i)
{
  int   Q = p7O_NQF(pp->M);
  int   q = (k-1) % Q;
  int   r = (k-1) / Q;
  union { __m128 v; float p[4]; } u;

  switch (scur) {
  case p7T_M:
    u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C0);  return u.p[r];
  case p7T_I:
    u.v = IMO_FS(pp->dpf[i], q);              return u.p[r];
  case p7T_N: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_N]; break;
  case p7T_C: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_C]; break;
  case p7T_J: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_J]; break;
  default: break;
  }
  return 0.0f;
}

static inline int
oa_select_m_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m128 *tp    = om_fs->tfv + 7*q;
  union { __m128 v; float p[4]; } u, tv;
  float   path[4];
  int     state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };
  __m128  mpv, dpv, ipv;

  if (q > 0) {
    mpv = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_M];
    dpv = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_D];
    ipv = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_I];
  } else {
    __m128 neg_inf = _mm_set1_ps(-eslINFINITY);
    mpv = esl_sse_rightshift_ps(ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_M], neg_inf);
    dpv = esl_sse_rightshift_ps(ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_D], neg_inf);
    ipv = esl_sse_rightshift_ps(ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_I], neg_inf);
  }

  tv.v = *tp;  path[3] = ((tv.p[r] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_B]);  tp++;
  u.v  = mpv;  tv.v = *tp;  path[0] = ((tv.p[r] == 0.0f) ? -eslINFINITY : u.p[r]);  tp++;
  u.v  = ipv;  tv.v = *tp;  path[1] = ((tv.p[r] == 0.0f) ? -eslINFINITY : u.p[r]);  tp++;
  u.v  = dpv;  tv.v = *tp;  path[2] = ((tv.p[r] == 0.0f) ? -eslINFINITY : u.p[r]);
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
oa_select_d_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  union { __m128 v; float p[4]; } mpv, dpv, tmdv, tddv;
  float   path[2];

  if (q > 0) {
    mpv.v  = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_M];
    dpv.v  = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_D];
    tmdv.v = om_fs->tfv[7*(q-1) + p7O_MD];
    tddv.v = om_fs->tfv[7*Q     + (q-1)];
  } else {
    __m128 neg_inf = _mm_set1_ps(-eslINFINITY);
    mpv.v  = esl_sse_rightshift_ps(ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_M], neg_inf);
    dpv.v  = esl_sse_rightshift_ps(ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_D], neg_inf);
    tmdv.v = esl_sse_rightshiftz_float(om_fs->tfv[7*(Q-1) + p7O_MD]);
    tddv.v = esl_sse_rightshiftz_float(om_fs->tfv[8*Q - 1]);
  }

  path[0] = ((tmdv.p[r] == 0.0f) ? -eslINFINITY : mpv.p[r]);
  path[1] = ((tddv.p[r] == 0.0f) ? -eslINFINITY : dpv.p[r]);
  return ((path[0] >= path[1]) ? p7T_M : p7T_D);
}

static inline int
oa_select_i_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q      = p7O_NQF(ox->M);
  int     q      = (k-1) % Q;
  int     r      = (k-1) / Q;
  int     prev_i = (i >= 3) ? i-3 : 0;
  __m128 *tp     = om_fs->tfv + 7*q + p7O_MI;
  union { __m128 v; float p[4]; } tv, mpv, ipv;
  float   path[2];

  mpv.v = ox->dpf[prev_i][q*p7X_NSCELLS + p7X_M];
  tv.v  = *tp;  path[0] = ((tv.p[r] == 0.0f) ? -eslINFINITY : mpv.p[r]);  tp++;
  ipv.v = ox->dpf[prev_i][q*p7X_NSCELLS + p7X_I];
  tv.v  = *tp;  path[1] = ((tv.p[r] == 0.0f) ? -eslINFINITY : ipv.p[r]);
  return ((path[0] >= path[1]) ? p7T_M : p7T_I);
}

static inline int
oa_select_n_fs(int i)
{
  return ((i == 0) ? p7T_S : p7T_N);
}

static inline int
oa_select_c_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp, const P7_OMX *ox, int i)
{
  int   L     = ox->L;
  float t1    = (om_fs->xf[p7O_C][p7O_LOOP] == 0.0f) ? 0.0f : 1.0f;
  float t2    = (om_fs->xf[p7O_E][p7O_MOVE] == 0.0f) ? 0.0f : 1.0f;
  float path[4];
  int   state[4] = { p7T_C, p7T_C, p7T_C, p7T_E };

  if (i < 4) return p7T_E;

  path[0] = (t1 == 0.0f) ? -eslINFINITY :
            ox->xmx[(i-3)*p7X_NXCELLS+p7X_C] + pp->xmx[ i   *p7X_NXCELLS+p7X_C];
  path[1] = (i < L && t1 != 0.0f) ?
            ox->xmx[(i-2)*p7X_NXCELLS+p7X_C] + pp->xmx[(i+1)*p7X_NXCELLS+p7X_C] : -eslINFINITY;
  path[2] = (i < L-1 && t1 != 0.0f) ?
            ox->xmx[(i-1)*p7X_NXCELLS+p7X_C] + pp->xmx[(i+2)*p7X_NXCELLS+p7X_C] : -eslINFINITY;
  path[3] = (t2 == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_E];
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
oa_select_j_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp, const P7_OMX *ox, int i)
{
  float path[2];
  int   state[2] = { p7T_J, p7T_E };

  if (i <= 5) return p7T_E;

  path[0] = (om_fs->xf[p7O_J][p7O_LOOP] == 0.0f) ? -eslINFINITY :
            ox->xmx[i*p7X_NXCELLS+p7X_J] + pp->xmx[i*p7X_NXCELLS+p7X_J];
  path[1] = (om_fs->xf[p7O_E][p7O_LOOP] == 0.0f) ? -eslINFINITY :
            ox->xmx[i*p7X_NXCELLS+p7X_E];
  return state[esl_vec_FArgMax(path, 2)];
}

static inline int
oa_select_e_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int *ret_k)
{
  int    Q    = p7O_NQF(ox->M);
  union { __m128 v; float p[4]; } u;
  float  max  = -eslINFINITY;
  int    smax = p7T_M, kmax = 1;
  int    k, q, r;

  for (k = 1; k <= ox->M; k++) {
    q = (k-1) % Q;
    r = (k-1) / Q;
    u.v = ox->dpf[i][q * p7X_NSCELLS + p7X_M];
    if (u.p[r] > max) { max = u.p[r]; smax = p7T_M; kmax = k; }
    u.v = ox->dpf[i][q * p7X_NSCELLS + p7X_D];
    if (u.p[r] > max) { max = u.p[r]; smax = p7T_D; kmax = k; }
  }
  *ret_k = kmax;
  return smax;
}

static inline int
oa_select_b_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i)
{
  float path[2];
  path[0] = (om_fs->xf[p7O_N][p7O_MOVE] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_N];
  path[1] = (om_fs->xf[p7O_J][p7O_MOVE] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_J];
  return ((path[0] > path[1]) ? p7T_N : p7T_J);
}

static inline int
oa_select_codon_fs(const P7_OMX *pp, int i, int k)
{
  int   Q = p7O_NQF(pp->M);
  int   q = (k-1) % Q;
  int   r = (k-1) / Q;
  union { __m128 v; float p[4]; } u;
  float codon[5];

  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C1);  codon[0] = u.p[r];
  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C2);  codon[1] = u.p[r];
  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C3);  codon[2] = u.p[r];
  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C4);  codon[3] = u.p[r];
  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C5);  codon[4] = u.p[r];
  return esl_vec_FArgMax(codon, 5) + 1;
}


/* Function:  p7_OATrace_Frameshift_sse()
 *
 * Purpose:   SSE implementation of frameshift optimal accuracy traceback.
 *            See p7_OATrace_Frameshift() documentation for details.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if trace is not empty.
 *            <eslEMEM> on allocation error.
 */
int
p7_OATrace_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                           const P7_OMX *ox, P7_TRACE *tr)
{
  int   i    = ox->L;
  int   k    = 0;
  int   c    = 0;
  float postprob;
  int   sprv, scur;
  int   status;

  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace not empty; needs to be Reuse()'d?");

  if ((status = p7_trace_fs_AppendWithPP(tr, p7T_T, k, i, c, 0.0f))     != eslOK) return status;
  if ((status = p7_trace_fs_AppendWithPP(tr, p7T_C, k, i, c, 0.0f))     != eslOK) return status;

  sprv = p7T_C;
  while (sprv != p7T_S)
    {
      switch (sprv) {
      case p7T_M: scur = oa_select_m_fs(om_fs,     ox, i, k);        k--;    break;
      case p7T_D: scur = oa_select_d_fs(om_fs,     ox, i, k);        k--;    break;
      case p7T_I: scur = oa_select_i_fs(om_fs,     ox, i, k);        i -= 3; break;
      case p7T_N: scur = oa_select_n_fs(i);                                   break;
      case p7T_C: scur = oa_select_c_fs(om_fs, pp, ox, i);                    break;
      case p7T_J: scur = oa_select_j_fs(om_fs, pp, ox, i);                    break;
      case p7T_E: scur = oa_select_e_fs(om_fs,     ox, i, &k);                break;
      case p7T_B: scur = oa_select_b_fs(om_fs,     ox, i);                    break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in OA FS traceback");
      }
      if (scur == -1) ESL_EXCEPTION(eslEINVAL, "OA FS SSE traceback choice failed");

      postprob = oa_get_postprob_fs(pp, scur, sprv, k, i);

      if (scur == p7T_M) c = oa_select_codon_fs(pp, i, k);
      else               c = 0;

      if ((status = p7_trace_fs_AppendWithPP(tr, scur, k, i, c, postprob)) != eslOK) return status;

      if ( (scur == p7T_N || scur == p7T_C || scur == p7T_J) && scur == sprv) i--;
      sprv = scur;
      i   -= c;
    }

  tr->M = om_fs->M;
  tr->L = ox->L;
  return p7_trace_fs_Reverse(tr);
}

#endif /* eslENABLE_SSE */
