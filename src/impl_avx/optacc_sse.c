/* Optimal accuracy alignment; SSE implementations for impl_avx dispatch.
 * These are the _sse-suffixed versions called by the runtime dispatchers
 * in optacc.c.
 */

#include "p7_config.h"

#include <string.h>

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

/* oa_row_sse()
 * OA row <i> into <dpc> and row i of the OA special states <xmx>, from
 * OA row i-1 in <dpp> and the posterior row <ppp> with its special
 * states <ppx> (E, N, J, B, C of row i). Row 0 is the initialization
 * (<ppp>, <ppx>, <dpp> unused). <dpc> may be <ppp>: each posterior
 * vector is read before the OA vector at its place is written. Same
 * calculation as each row of p7_OptimalAccuracy_sse(), which calls it.
 */
static void
oa_row_sse(const P7_OPROFILE *om, const __m128 *ppp, const float *ppx, float *xmx, __m128 *dpc, const __m128 *dpp, int i)
{
  register __m128 mpv, dpv, ipv;
  register __m128 sv;
  register __m128 xEv;
  register __m128 xBv;
  register __m128 dcv;
  __m128 *tp;
  __m128 zerov = _mm_setzero_ps();
  __m128 infv  = _mm_set1_ps(-eslINFINITY);
  int M = om->M;
  int Q = p7O_NQF(M);
  int q;
  int j;
  float t1, t2;


  if (i == 0)
    {
      for (q = 0; q < Q; q++) MMO(dpc, q) = IMO(dpc,q) = DMO(dpc,q) = infv;
      XMXo(0, p7X_E)    = -eslINFINITY;
      XMXo(0, p7X_N)    = 0.;
      XMXo(0, p7X_J)    = -eslINFINITY;
      XMXo(0, p7X_B)    = 0.;
      XMXo(0, p7X_C)    = -eslINFINITY;
      return;
    }

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

  t1 = ( (om->xf[p7O_J][p7O_LOOP] == 0.0) ? 0.0 : xmx[(i-1)*p7X_NXCELLS+p7X_J] + ppx[p7X_J]);
  t2 = ( (om->xf[p7O_E][p7O_LOOP] == 0.0) ? 0.0 : xmx[   i *p7X_NXCELLS+p7X_E]);
  xmx[i*p7X_NXCELLS+p7X_J] = ESL_MAX(t1, t2);

  t1 = ( (om->xf[p7O_C][p7O_LOOP] == 0.0) ? 0.0 : xmx[(i-1)*p7X_NXCELLS+p7X_C] + ppx[p7X_C]);
  t2 = ( (om->xf[p7O_E][p7O_MOVE] == 0.0) ? 0.0 : xmx[   i *p7X_NXCELLS+p7X_E]);
  xmx[i*p7X_NXCELLS+p7X_C] = ESL_MAX(t1, t2);

  xmx[i*p7X_NXCELLS+p7X_N] = ((om->xf[p7O_N][p7O_LOOP] == 0.0) ? 0.0 : xmx[(i-1)*p7X_NXCELLS+p7X_N] + ppx[p7X_N]);

  t1 = ( (om->xf[p7O_N][p7O_MOVE] == 0.0) ? 0.0 : xmx[i*p7X_NXCELLS+p7X_N]);
  t2 = ( (om->xf[p7O_J][p7O_MOVE] == 0.0) ? 0.0 : xmx[i*p7X_NXCELLS+p7X_J]);
  xmx[i*p7X_NXCELLS+p7X_B] = ESL_MAX(t1, t2);
}


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
  int i;

  ox->M = om->M;
  ox->L = pp->L;
  oa_row_sse(om, NULL, NULL, ox->xmx, ox->dpf[0], NULL, 0);
  for (i = 1; i <= pp->L; i++)
    oa_row_sse(om, pp->dpf[i], pp->xmx + i*p7X_NXCELLS, ox->xmx, ox->dpf[i], ox->dpf[i-1], i);

  *ret_e = ox->xmx[pp->L*p7X_NXCELLS+p7X_C];
  return eslOK;
}


/* OA traceback helpers (standard, non-FS) */

static inline float oa_get_postprob(const P7_OPROFILE *om, const P7_OMX *ox, int scur, int sprv, int k, int i);
static inline int oa_select_m(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int k);
static inline int oa_select_d(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int k);
static inline int oa_select_i(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int k);
static inline int oa_select_n(int i);
static inline int oa_select_c(const P7_OPROFILE *om,                   const P7_OMX *ox, int i);
static inline int oa_select_j(const P7_OPROFILE *om,                   const P7_OMX *ox, int i);
static inline int oa_select_e(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int *ret_k);
static inline int oa_select_b(const P7_OPROFILE *om,                   const P7_OMX *ox, int i);

/* oa_get_postprob()
 * Posterior probability of state <scur> at (<i>,<k>), recovered from
 * the OA matrix <ox>. The fill set an M or I cell to its best
 * predecessor score plus the cell's posterior, with each predecessor
 * term set to 0 when its transition is not allowed and the M, I, D
 * cells of column 0 taken as -inf; the posterior is the cell's score
 * minus that same predecessor score. The N, C and J loops add the
 * posterior to the previous row's score.
 */
static inline float
oa_get_postprob(const P7_OPROFILE *om, const P7_OMX *ox, int scur, int sprv, int k, int i)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;		/* (q,r) is the position of cell (i,k)    */
  int     r     = (k-1) / Q;
  int     qp    = (k-2) % Q;		/* (qp,rp) is the position of cell (i,k-1) */
  int     rp    = (k-2) / Q;
  __m128  *tp    = om->tfv + 7*q;
  union { __m128 v; float p[4]; } u, tv;
  float   mv, iv, dv, prv;

  switch (scur) {
  case p7T_M:
    if (k > 1) {
      u.v = ox->dpf[i-1][qp*3 + p7X_M];  mv = u.p[rp];
      u.v = ox->dpf[i-1][qp*3 + p7X_I];  iv = u.p[rp];
      u.v = ox->dpf[i-1][qp*3 + p7X_D];  dv = u.p[rp];
    } else {
      mv = iv = dv = -eslINFINITY;
    }
    tv.v = tp[p7O_BM];  prv =              ((tv.p[r] > 0.0) ? ox->xmx[(i-1)*p7X_NXCELLS+p7X_B] : 0.0);
    tv.v = tp[p7O_MM];  prv = ESL_MAX(prv, ((tv.p[r] > 0.0) ? mv : 0.0));
    tv.v = tp[p7O_IM];  prv = ESL_MAX(prv, ((tv.p[r] > 0.0) ? iv : 0.0));
    tv.v = tp[p7O_DM];  prv = ESL_MAX(prv, ((tv.p[r] > 0.0) ? dv : 0.0));
    u.v  = ox->dpf[i][q*3 + p7X_M];
    return u.p[r] - prv;
  case p7T_I:
    u.v = ox->dpf[i-1][q*3 + p7X_M];  tv.v = tp[p7O_MI];  prv =              ((tv.p[r] > 0.0) ? u.p[r] : 0.0);
    u.v = ox->dpf[i-1][q*3 + p7X_I];  tv.v = tp[p7O_II];  prv = ESL_MAX(prv, ((tv.p[r] > 0.0) ? u.p[r] : 0.0));
    u.v = ox->dpf[i][q*3 + p7X_I];
    return u.p[r] - prv;
  case p7T_N: if (sprv == scur) return ox->xmx[i*p7X_NXCELLS+p7X_N] - ox->xmx[(i-1)*p7X_NXCELLS+p7X_N];
  case p7T_C: if (sprv == scur) return ox->xmx[i*p7X_NXCELLS+p7X_C] - ox->xmx[(i-1)*p7X_NXCELLS+p7X_C];
  case p7T_J: if (sprv == scur) return ox->xmx[i*p7X_NXCELLS+p7X_J] - ox->xmx[(i-1)*p7X_NXCELLS+p7X_J];
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
oa_select_c(const P7_OPROFILE *om, const P7_OMX *ox, int i)
{
  /* The fill set C(i) to the larger of C(i-1) + pp_C(i) and E(i), each
   * term 0 when its transition is not allowed. The C->C loop is
   * taken when it is allowed and either E->C is not allowed (and
   * C(i-1) is reachable) or the loop term is larger, which is
   * C(i) > E(i); on a tie E wins, as in the fill.
   */
  if (om->xf[p7O_C][p7O_LOOP]  == 0.0) return p7T_E;
  if (om->xf[p7O_E][p7O_MOVE] == 0.0) return ((ox->xmx[(i-1)*p7X_NXCELLS+p7X_C] > -eslINFINITY) ? p7T_C : p7T_E);
  return ((ox->xmx[i*p7X_NXCELLS+p7X_C] > ox->xmx[i*p7X_NXCELLS+p7X_E]) ? p7T_C : p7T_E);
}

static inline int
oa_select_j(const P7_OPROFILE *om, const P7_OMX *ox, int i)
{
  /* The fill set J(i) to the larger of J(i-1) + pp_J(i) and E(i), each
   * term 0 when its transition is not allowed. The J->J loop is
   * taken when it is allowed and either E->J is not allowed (and
   * J(i-1) is reachable) or the loop term is larger, which is
   * J(i) > E(i); on a tie E wins, as in the fill.
   */
  if (om->xf[p7O_J][p7O_LOOP]  == 0.0) return p7T_E;
  if (om->xf[p7O_E][p7O_LOOP] == 0.0) return ((ox->xmx[(i-1)*p7X_NXCELLS+p7X_J] > -eslINFINITY) ? p7T_J : p7T_E);
  return ((ox->xmx[i*p7X_NXCELLS+p7X_J] > ox->xmx[i*p7X_NXCELLS+p7X_E]) ? p7T_J : p7T_E);
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
 *            Reads only the OA matrix <ox>; posterior probabilities are
 *            recovered from it by subtraction, so <ox> may have been
 *            written over the posterior decoding matrix.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> if trace is not empty.
 */
int
p7_OATrace_sse(const P7_OPROFILE *om, const P7_OMX *ox, P7_TRACE *tr)
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
      case p7T_C: s1 = oa_select_c(om,     ox, i);               break;
      case p7T_J: s1 = oa_select_j(om,     ox, i);               break;
      case p7T_E: s1 = oa_select_e(om,     ox, i, &k);           break;
      case p7T_B: s1 = oa_select_b(om,     ox, i);               break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }
      if (s1 == -1) ESL_EXCEPTION(eslEINVAL, "OA traceback choice failed");

      postprob = oa_get_postprob(om, ox, s1, s0, k, i);
      if ((status = p7_trace_AppendWithPP(tr, s1, k, i, postprob)) != eslOK) return status;

      if ( (s1 == p7T_N || s1 == p7T_J || s1 == p7T_C) && s1 == s0) i--;
      s0 = s1;
    }
  tr->M = om->M;
  tr->L = ox->L;
  return p7_trace_Reverse(tr);
}

/*****************************************************************
 * OA sweep over a checkpointed Backward, and its trace
 *****************************************************************/

/* oat_blend_sse()
 * Lanes of <b> where <mask> is set, of <a> elsewhere (SSE2).
 */
static inline __m128
oat_blend_sse(__m128 a, __m128 b, __m128 mask)
{
  return _mm_or_ps(_mm_and_ps(mask, b), _mm_andnot_ps(mask, a));
}

/* oat_ppcodes4()
 * Posterior codes (p7OAT_PPCODE(), 0..10) of the 4 lanes of <v> as 4
 * lane masks, bit b of each lane's code in <m[b]>.
 */
static inline void
oat_ppcodes4(__m128 v, int m[4])
{
  union { __m128 v; float p[4]; } u;
  int r, c;

  u.v  = v;
  m[0] = m[1] = m[2] = m[3] = 0;
  for (r = 0; r < 4; r++)
    {
      c     = p7OAT_PPCODE(u.p[r]);
      m[0] |= ( c       & 1) << r;
      m[1] |= ((c >> 1) & 1) << r;
      m[2] |= ((c >> 2) & 1) << r;
      m[3] |= ((c >> 3) & 1) << r;
    }
}

/* oat_store_row_sse()
 * Store into <ox> what the OA trace reads from OA row <i> (<dpc>, with
 * the OA special states in <xmx>): for each k, the predecessors of an
 * M at (i+1,k), an I at (i+1,k) and a D at (i,k), chosen with the same
 * comparisons, shifts and ties as p7_OATrace_sse() makes on the OA
 * matrix; the posterior codes of M and I at (i,k) from the posterior
 * row <ppp> (all 0 for row 0, <ppp> NULL); the row's OA special states
 * and the posteriors of N, J, C (<ppx>, NULL for row 0); and the cell
 * an E at row i comes from, in p7_OATrace_sse()'s scan order.
 */
static void
oat_store_row_sse(const P7_OPROFILE *om, P7_OATMX *ox, const __m128 *dpc, const float *xmx, const __m128 *ppp, const float *ppx, int i)
{
  int      M     = om->M;
  int      Q     = p7O_NQF(M);
  __m128   zerov = _mm_setzero_ps();
  __m128   infv  = _mm_set1_ps(-eslINFINITY);
  __m128   xBv   = _mm_set1_ps(xmx[i*p7X_NXCELLS+p7X_B]);
  __m128   mpv, ipv, dpv, tmdv, tddv;
  __m128   pM, pI, pD, pB, best, gt, b0, b1;
  __m128  *tp;
  union { __m128 v; float p[4]; } u;
  float    max;
  int      kmax, smax;
  int      mm[4], mi[4];
  int      q, r, c;

  for (q = 0; q < Q; q++)
    {
      tp = om->tfv + 7*q;

      /* row i at k-1 */
      if (q > 0) {
        mpv  = MMO(dpc, q-1);
        ipv  = IMO(dpc, q-1);
        dpv  = DMO(dpc, q-1);
        tmdv = om->tfv[7*(q-1) + p7O_MD];
        tddv = om->tfv[7*Q + (q-1)];
      } else {
        mpv  = esl_sse_rightshiftz_float(MMO(dpc, Q-1));
        ipv  = esl_sse_rightshiftz_float(IMO(dpc, Q-1));
        dpv  = esl_sse_rightshiftz_float(DMO(dpc, Q-1));
        tmdv = esl_sse_rightshiftz_float(om->tfv[7*(Q-1) + p7O_MD]);
        tddv = esl_sse_rightshiftz_float(om->tfv[8*Q-1]);
      }

      /* M at (i+1,k): first best of M, I, D, B */
      pB   = oat_blend_sse(xBv, infv, _mm_cmpeq_ps(tp[p7O_BM], zerov));
      pM   = oat_blend_sse(mpv, infv, _mm_cmpeq_ps(tp[p7O_MM], zerov));
      pI   = oat_blend_sse(ipv, infv, _mm_cmpeq_ps(tp[p7O_IM], zerov));
      pD   = oat_blend_sse(dpv, infv, _mm_cmpeq_ps(tp[p7O_DM], zerov));
      best = pM;
      gt   = _mm_cmpgt_ps(pI, best);  best = oat_blend_sse(best, pI, gt);  b0 = gt;                    b1 = zerov;
      gt   = _mm_cmpgt_ps(pD, best);  best = oat_blend_sse(best, pD, gt);  b0 = _mm_andnot_ps(gt, b0); b1 = _mm_or_ps(b1, gt);
      gt   = _mm_cmpgt_ps(pB, best);                                           b0 = _mm_or_ps(b0, gt);    b1 = _mm_or_ps(b1, gt);
      OATO_SET_STRIPE4(ox, i, p7OAT_MP, q, _mm_movemask_ps(b0));
      OATO_SET_STRIPE4(ox, i, p7OAT_MP+1, q, _mm_movemask_ps(b1));

      /* I at (i+1,k): M on ties */
      pM = oat_blend_sse(MMO(dpc, q), infv, _mm_cmpeq_ps(tp[p7O_MI], zerov));
      pI = oat_blend_sse(IMO(dpc, q), infv, _mm_cmpeq_ps(tp[p7O_II], zerov));
      OATO_SET_STRIPE4(ox, i, p7OAT_IP, q, _mm_movemask_ps(_mm_cmpgt_ps(pI, pM)));

      /* D at (i,k): M on ties */
      pM = oat_blend_sse(mpv, infv, _mm_cmpeq_ps(tmdv, zerov));
      pD = oat_blend_sse(dpv, infv, _mm_cmpeq_ps(tddv, zerov));
      OATO_SET_STRIPE4(ox, i, p7OAT_DP, q, _mm_movemask_ps(_mm_cmpgt_ps(pD, pM)));

      /* posterior codes of M and I */
      if (ppp != NULL) { oat_ppcodes4(ppp[q*3 + p7X_M], mm); oat_ppcodes4(ppp[q*3 + p7X_I], mi); }
      else             { mm[0] = mm[1] = mm[2] = mm[3] = mi[0] = mi[1] = mi[2] = mi[3] = 0; }
      for (c = 0; c < 4; c++)
        {
          OATO_SET_STRIPE4(ox, i, p7OAT_PPM + c, q, mm[c]);
          OATO_SET_STRIPE4(ox, i, p7OAT_PPI + c, q, mi[c]);
        }
    }

  ox->xmx[i*p7OAT_NXCELLS + p7OAT_N]   = xmx[i*p7X_NXCELLS + p7X_N];
  ox->xmx[i*p7OAT_NXCELLS + p7OAT_J]   = xmx[i*p7X_NXCELLS + p7X_J];
  ox->xmx[i*p7OAT_NXCELLS + p7OAT_C]   = xmx[i*p7X_NXCELLS + p7X_C];
  ox->xmx[i*p7OAT_NXCELLS + p7OAT_B]   = xmx[i*p7X_NXCELLS + p7X_B];
  ox->xmx[i*p7OAT_NXCELLS + p7OAT_E]   = xmx[i*p7X_NXCELLS + p7X_E];
  ox->xmx[i*p7OAT_NXCELLS + p7OAT_PPN] = ppx ? ppx[p7X_N] : 0.0f;
  ox->xmx[i*p7OAT_NXCELLS + p7OAT_PPJ] = ppx ? ppx[p7X_J] : 0.0f;
  ox->xmx[i*p7OAT_NXCELLS + p7OAT_PPC] = ppx ? ppx[p7X_C] : 0.0f;

  /* the cell an E at row i comes from, as oa_select_e() scans the row */
  max  = -eslINFINITY;
  smax = p7T_M;
  kmax = 1;
  for (q = 0; q < Q; q++)
    {
      u.v = MMO(dpc, q); for (r = 0; r < 4; r++) if (u.p[r] >= max) { max = u.p[r]; smax = p7T_M; kmax = r*Q + q + 1; }
      u.v = DMO(dpc, q); for (r = 0; r < 4; r++) if (u.p[r] >  max) { max = u.p[r]; smax = p7T_D; kmax = r*Q + q + 1; }
    }
  ox->ek[i] = (smax == p7T_M) ? kmax : -kmax;
}


/* Function:  p7_OptimalAccuracySweep_sse()
 * Synopsis:  OA alignment in one forward sweep over a checkpointed Backward, SSE.
 *
 * Purpose:   For profile <om> and sequence <dsq> of length <L>, with the
 *            checkpointed Backward matrix <cx> filled by
 *            <p7_BackwardChk_sse()>, compute in one sweep through rows
 *            1..L the Forward, posterior and optimal accuracy (OA) rows,
 *            computing each block of Backward rows again from its
 *            checkpoint as the sweep reaches it, and store in <ox> what
 *            <p7_OATraceSweep_sse()> needs for the OA trace. Rows of
 *            Forward, posterior and OA values are kept only as long as
 *            the next row needs them.
 *
 *            Each row is calculated as <p7_Forward_sse()>,
 *            <p7_Decoding_sse()> and <p7_OptimalAccuracy_sse()> calculate
 *            it. The checkpointed Backward has its own scale factors, so
 *            the posteriors are those <p7_Decoding_sse()> gives for a
 *            Backward matrix with its own scale factors: normalized by
 *            the Backward score, and equal to those from a Backward
 *            scaled like Forward up to rounding.
 *
 *            If <opt_null2> is non-NULL, the null2 odds ratios of the
 *            sequence are also calculated into it, as
 *            <p7_Null2_ByExpectation_sse()> calculates them from the
 *            posteriors: the posterior rows are summed as the sweep
 *            computes them, in the same order.
 *
 *            <ox> must have been created with <p7OAT_NPLANES> planes;
 *            its layout is set to 4 lanes and it is grown for <om->M> by
 *            <L>. On return <cx> has its last block loaded.
 *
 * Args:      dsq       - sequence in digitized form, 1..L
 *            L         - length of dsq (>= 1)
 *            om        - profile
 *            cx        - checkpointed Backward matrix from p7_BackwardChk_sse()
 *            ox        - RESULT: what the OA trace reads
 *            opt_null2 - optRETURN: null2 odds ratios, <0..Kp-1>; caller allocated
 *            opt_fsc   - optRETURN: Forward score in nats
 *            ret_oasc  - RETURN: OA score, the expected number of
 *                        correctly aligned residues
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; <eslEINVAL> if <ox> or <cx>
 *            isn't set up for this computation; <eslERANGE> if the
 *            posterior scale correction overflows (as <p7_Decoding_sse()>).
 */
int
p7_OptimalAccuracySweep_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, P7_OATMX *ox,
                            float *opt_null2, float *opt_fsc, float *ret_oasc)
{
  int      M      = om->M;
  int      Q      = p7O_NQF(M);
  P7_OMX  *wk     = NULL;   /* rows 0,1: Forward; 2,3: OA; 4: posteriors; 5: posterior sums; xmx: Forward specials */
  float   *oxmx   = NULL;
  float    ppx[p7X_NXCELLS];
  float    cntx[p7X_NXCELLS];
  float    totscale = 0.0f;
  float    scaleproduct;
  __m128  *pp, *cnt;
  int      i, q, b;
  int      status;

  if (ox->nplanes != p7OAT_NPLANES) ESL_EXCEPTION(eslEINVAL, "P7_OATMX not set up for the standard OA sweep");
  if (cx->L != L || cx->M != M || cx->dir != p7_CHK_BACKWARD || cx->nscells != p7X_NSCELLS) ESL_EXCEPTION(eslEINVAL, "checkpointed Backward matrix doesn't match");
  ox->lanes = 4;
  if ((status = p7_oatmx_GrowTo(ox, M, L)) != eslOK) return status;

  if ((wk = p7_omx_Create_sse(M, 5, L)) == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(oxmx, sizeof(float) * (L+1) * p7X_NXCELLS);
  pp  = wk->dpf[4];
  cnt = wk->dpf[5];

  scaleproduct = 1.0 / cx->xmx[p7X_N];
  for (q = 0; q < p7X_NXCELLS; q++) cntx[q] = 0.0f;

  p7_ForwardRow_sse(dsq, om, wk->xmx, &totscale, wk->dpf[0], NULL, 0);
  oa_row_sse(om, NULL, NULL, oxmx, wk->dpf[2], NULL, 0);
  oat_store_row_sse(om, ox, wk->dpf[2], oxmx, NULL, NULL, 0);

  for (i = 1; i <= L; i++)
    {
      b = (i / cx->blockR) * cx->blockR;
      if (b != cx->b && (status = p7_BackwardChkLoad_sse(dsq, om, cx, b)) != eslOK) goto ERROR;

      p7_ForwardRow_sse(dsq, om, wk->xmx, &totscale, wk->dpf[i%2], wk->dpf[(i-1)%2], i);
      p7_DecodingRow_sse(om, wk->xmx, cx->xmx, wk->dpf[i%2], cx->dpf[i], pp, ppx, i, scaleproduct);
      scaleproduct *= wk->xmx[i*p7X_NXCELLS+p7X_SCALE] / cx->xmx[i*p7X_NXCELLS+p7X_SCALE];

      oa_row_sse(om, pp, ppx, oxmx, wk->dpf[2 + i%2], wk->dpf[2 + (i-1)%2], i);
      oat_store_row_sse(om, ox, wk->dpf[2 + i%2], oxmx, pp, ppx, i);

      /* posterior sums for null2, added up as p7_Null2_ByExpectation_sse() does */
      if (i == 1)
        {
          memcpy(cnt, pp, sizeof(__m128) * 3 * Q);
          cntx[p7X_N] = ppx[p7X_N];
          cntx[p7X_C] = ppx[p7X_C];
          cntx[p7X_J] = ppx[p7X_J];
        }
      else
        {
          for (q = 0; q < Q; q++)
            {
              cnt[q*3 + p7X_M] = _mm_add_ps(pp[q*3 + p7X_M], cnt[q*3 + p7X_M]);
              cnt[q*3 + p7X_I] = _mm_add_ps(pp[q*3 + p7X_I], cnt[q*3 + p7X_I]);
            }
          cntx[p7X_N] += ppx[p7X_N];
          cntx[p7X_C] += ppx[p7X_C];
          cntx[p7X_J] += ppx[p7X_J];
        }
    }
  if (isinf(scaleproduct)) { status = eslERANGE; goto ERROR; }

  if (opt_null2 != NULL) p7_Null2_FromCounts_sse(om, cnt, cntx, L, opt_null2);
  if (opt_fsc   != NULL) *opt_fsc = totscale + log(wk->xmx[L*p7X_NXCELLS+p7X_C] * om->xf[p7O_C][p7O_MOVE]);
  *ret_oasc = oxmx[L*p7X_NXCELLS+p7X_C];

  p7_omx_Destroy_sse(wk);
  free(oxmx);
  return eslOK;

 ERROR:
  if (wk) p7_omx_Destroy_sse(wk);
  free(oxmx);
  return status;
}


/* Function:  p7_OATraceSweep_sse()
 * Synopsis:  OA trace from what the SSE OA sweep stored.
 *
 * Purpose:   The OA trace for profile <om> from <ox>, filled by
 *            <p7_OptimalAccuracySweep_sse()>, into the empty trace <tr>
 *            (created with <p7_trace_CreateWithPP()>). The choices are
 *            those <p7_OATrace_sse()> makes on the same OA scores. The
 *            posterior probability of each residue is the value its PP
 *            character stands for (<p7OAT_PPVALUE()>), from the
 *            posterior row, as <p7_alidisplay_EncodePostProb()> prints
 *            it.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if <tr> isn't
 *            empty or a bogus state is reached.
 */
int
p7_OATraceSweep_sse(const P7_OPROFILE *om, const P7_OATMX *ox, P7_TRACE *tr)
{
  static const int state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };   /* indexed by p7OAT_FROM_* */
  int   i = ox->L;
  int   k = 0;
  int   s0, s1;
  float path[2];
  float postprob;
  int   status;

  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace not empty; needs to be Reuse()'d?");

  if ((status = p7_trace_AppendWithPP(tr, p7T_T, k, i, 0.0)) != eslOK) return status;
  if ((status = p7_trace_AppendWithPP(tr, p7T_C, k, i, 0.0)) != eslOK) return status;

  s0 = tr->st[tr->N-1];
  while (s0 != p7T_S)
    {
      switch (s0) {
      case p7T_M: s1 = state[OATO_MPRED(ox, i-1, k)];                              k--; i--; break;
      case p7T_D: s1 = OATO_DPRED(ox, i, k) ? p7T_D : p7T_M;                        k--;      break;
      case p7T_I: s1 = (OATO_IPRED(ox, i-1, k) == p7OAT_FROM_I) ? p7T_I : p7T_M;         i--; break;
      case p7T_N: s1 = (i == 0) ? p7T_S : p7T_N;                                               break;
      case p7T_C:
        if      (om->xf[p7O_C][p7O_LOOP] == 0.0) s1 = p7T_E;
        else if (om->xf[p7O_E][p7O_MOVE] == 0.0) s1 = (ox->xmx[(i-1)*p7OAT_NXCELLS+p7OAT_C] > -eslINFINITY) ? p7T_C : p7T_E;
        else    s1 = (ox->xmx[i*p7OAT_NXCELLS+p7OAT_C] > ox->xmx[i*p7OAT_NXCELLS+p7OAT_E]) ? p7T_C : p7T_E;
        break;
      case p7T_J:
        if      (om->xf[p7O_J][p7O_LOOP] == 0.0) s1 = p7T_E;
        else if (om->xf[p7O_E][p7O_LOOP] == 0.0) s1 = (ox->xmx[(i-1)*p7OAT_NXCELLS+p7OAT_J] > -eslINFINITY) ? p7T_J : p7T_E;
        else    s1 = (ox->xmx[i*p7OAT_NXCELLS+p7OAT_J] > ox->xmx[i*p7OAT_NXCELLS+p7OAT_E]) ? p7T_J : p7T_E;
        break;
      case p7T_E:
        s1 = (ox->ek[i] > 0) ? p7T_M : p7T_D;
        k  = abs(ox->ek[i]);
        break;
      case p7T_B:
        path[0] = (om->xf[p7O_N][p7O_MOVE] == 0.0) ? -eslINFINITY : ox->xmx[i*p7OAT_NXCELLS+p7OAT_N];
        path[1] = (om->xf[p7O_J][p7O_MOVE] == 0.0) ? -eslINFINITY : ox->xmx[i*p7OAT_NXCELLS+p7OAT_J];
        s1 = (path[0] > path[1]) ? p7T_N : p7T_J;
        break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }

      switch (s1) {
      case p7T_M: postprob = p7OAT_PPVALUE(OATO_PPM(ox, i, k));                              break;
      case p7T_I: postprob = p7OAT_PPVALUE(OATO_PPI(ox, i, k));                              break;
      case p7T_N: postprob = (s1 == s0) ? ox->xmx[i*p7OAT_NXCELLS+p7OAT_PPN] : 0.0f;         break;
      case p7T_C: postprob = (s1 == s0) ? ox->xmx[i*p7OAT_NXCELLS+p7OAT_PPC] : 0.0f;         break;
      case p7T_J: postprob = (s1 == s0) ? ox->xmx[i*p7OAT_NXCELLS+p7OAT_PPJ] : 0.0f;         break;
      default:    postprob = 0.0f;                                                           break;
      }
      if ((status = p7_trace_AppendWithPP(tr, s1, k, i, postprob)) != eslOK) return status;

      if ( (s1 == p7T_N || s1 == p7T_J || s1 == p7T_C) && s1 == s0) i--;
      s0 = s1;
    }
  tr->M = om->M;
  tr->L = ox->L;
  return p7_trace_Reverse(tr);
}

#endif /* eslENABLE_SSE */
