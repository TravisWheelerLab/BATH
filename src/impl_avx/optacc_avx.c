/* Optimal accuracy alignment; AVX2 implementations for impl_avx dispatch.
 * These are the _avx-suffixed versions called by the runtime dispatchers
 * in optacc.c.
 */

#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <float.h>
#include <string.h>

#include <immintrin.h>

#include "easel.h"
#include "esl_avx.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"


/* avx_rightshift_ps()
 * Right-shift a 256-bit float vector by one position and insert infv at position 0.
 * Mirrors esl_sse_rightshift_ps() for AVX2.
 */
static inline __m256
avx_rightshift_ps(__m256 v, __m256 infv)
{
  return _mm256_blend_ps(esl_avx_rightshiftz_float(v), infv, 0x01);
}

/* avx_hmax_ps()
 * Compute horizontal max of all 8 floats in a 256-bit register.
 * Stores result in *ret. No equivalent in esl_avx.h.
 */
static inline void
avx_hmax_ps(__m256 v, float *ret)
{
  __m128 hi   = _mm256_extractf128_ps(v, 1);
  __m128 lo   = _mm256_castps256_ps128(v);
  __m128 max4 = _mm_max_ps(lo, hi);
  max4 = _mm_max_ps(max4, _mm_shuffle_ps(max4, max4, _MM_SHUFFLE(2,3,0,1)));
  max4 = _mm_max_ps(max4, _mm_shuffle_ps(max4, max4, _MM_SHUFFLE(1,0,3,2)));
  _mm_store_ss(ret, max4);
}


/*****************************************************************
 * 1. Standard optimal accuracy DP fill and traceback
 *****************************************************************/

/* oa_row_avx()
 * OA row <i> into <dpc> and row i of the OA special states <xmx>, from
 * OA row i-1 in <dpp> and the posterior row <ppp> with its special
 * states <ppx> (E, N, J, B, C of row i). Row 0 is the initialization
 * (<ppp>, <ppx>, <dpp> unused). <dpc> may be <ppp>: each posterior
 * vector is read before the OA vector at its place is written. Same
 * calculation as each row of p7_OptimalAccuracy_avx(), which calls it.
 */
static void
oa_row_avx(const P7_OPROFILE *om, const __m256 *ppp, const float *ppx, float *xmx, __m256 *dpc, const __m256 *dpp, int i)
{
  register __m256 mpv, dpv, ipv;
  register __m256 sv;
  register __m256 xEv;
  register __m256 xBv;
  register __m256 dcv;
  __m256 *tp;
  __m256 zerov = _mm256_setzero_ps();
  __m256 infv  = _mm256_set1_ps(-eslINFINITY);
  int M = om->M;
  int Q = p7O_NQF_AVX(M);
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

  tp  = om->tfv_avx;
  dcv = infv;
  xEv = infv;
  xBv = _mm256_set1_ps(XMXo(i-1, p7X_B));

  mpv = avx_rightshift_ps(MMO(dpp,Q-1), infv);
  dpv = avx_rightshift_ps(DMO(dpp,Q-1), infv);
  ipv = avx_rightshift_ps(IMO(dpp,Q-1), infv);
  for (q = 0; q < Q; q++)
    {
      sv  =                _mm256_and_ps(_mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ), xBv);  tp++;
      sv  = _mm256_max_ps(sv, _mm256_and_ps(_mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ), mpv)); tp++;
      sv  = _mm256_max_ps(sv, _mm256_and_ps(_mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ), ipv)); tp++;
      sv  = _mm256_max_ps(sv, _mm256_and_ps(_mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ), dpv)); tp++;
      sv  = _mm256_add_ps(sv, *ppp);                                                       ppp += 2;
      xEv = _mm256_max_ps(xEv, sv);

      mpv = MMO(dpp,q);
      dpv = DMO(dpp,q);
      ipv = IMO(dpp,q);

      MMO(dpc,q) = sv;
      DMO(dpc,q) = dcv;

      dcv = _mm256_and_ps(_mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ), sv); tp++;

      sv         =                _mm256_and_ps(_mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ), mpv);   tp++;
      sv         = _mm256_max_ps(sv, _mm256_and_ps(_mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ), ipv));  tp++;
      IMO(dpc,q) = _mm256_add_ps(sv, *ppp);                                                       ppp++;
    }

  dcv = avx_rightshift_ps(dcv, infv);
  tp  = om->tfv_avx + 7*Q;
  for (q = 0; q < Q; q++)
    {
      DMO(dpc, q) = _mm256_max_ps(dcv, DMO(dpc, q));
      dcv         = _mm256_and_ps(_mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ), DMO(dpc,q));   tp++;
    }

  for (j = 1; j < 8; j++)
    {
      dcv = avx_rightshift_ps(dcv, infv);
      tp  = om->tfv_avx + 7*Q;
      for (q = 0; q < Q; q++)
        {
          DMO(dpc, q) = _mm256_max_ps(dcv, DMO(dpc, q));
          dcv         = _mm256_and_ps(_mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ), dcv);   tp++;
        }
    }

  for (q = 0; q < Q; q++) xEv = _mm256_max_ps(xEv, DMO(dpc,q));

  avx_hmax_ps(xEv, &(XMXo(i,p7X_E)));

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


/* Function:  p7_OptimalAccuracy_avx()
 *
 * Purpose:   AVX2 implementation of optimal accuracy DP fill.
 *            See p7_OptimalAccuracy() documentation for details.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_OptimalAccuracy_avx(const P7_OPROFILE *om, const P7_OMX *pp, P7_OMX *ox, float *ret_e)
{
  int i;

  ox->M = om->M;
  ox->L = pp->L;
  oa_row_avx(om, NULL, NULL, ox->xmx, ox->dpf_avx[0], NULL, 0);
  for (i = 1; i <= pp->L; i++)
    oa_row_avx(om, pp->dpf_avx[i], pp->xmx + i*p7X_NXCELLS, ox->xmx, ox->dpf_avx[i], ox->dpf_avx[i-1], i);

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
  int     Q     = p7O_NQF_AVX(ox->M);
  int     q     = (k-1) % Q;		/* (q,r) is the position of cell (i,k)    */
  int     r     = (k-1) / Q;
  int     qp    = (k-2) % Q;		/* (qp,rp) is the position of cell (i,k-1) */
  int     rp    = (k-2) / Q;
  __m256  *tp    = om->tfv_avx + 7*q;
  union { __m256 v; float p[8]; } u, tv;
  float   mv, iv, dv, prv;

  switch (scur) {
  case p7T_M:
    if (k > 1) {
      u.v = ox->dpf_avx[i-1][qp*3 + p7X_M];  mv = u.p[rp];
      u.v = ox->dpf_avx[i-1][qp*3 + p7X_I];  iv = u.p[rp];
      u.v = ox->dpf_avx[i-1][qp*3 + p7X_D];  dv = u.p[rp];
    } else {
      mv = iv = dv = -eslINFINITY;
    }
    tv.v = tp[p7O_BM];  prv =              ((tv.p[r] > 0.0) ? ox->xmx[(i-1)*p7X_NXCELLS+p7X_B] : 0.0);
    tv.v = tp[p7O_MM];  prv = ESL_MAX(prv, ((tv.p[r] > 0.0) ? mv : 0.0));
    tv.v = tp[p7O_IM];  prv = ESL_MAX(prv, ((tv.p[r] > 0.0) ? iv : 0.0));
    tv.v = tp[p7O_DM];  prv = ESL_MAX(prv, ((tv.p[r] > 0.0) ? dv : 0.0));
    u.v  = ox->dpf_avx[i][q*3 + p7X_M];
    return u.p[r] - prv;
  case p7T_I:
    u.v = ox->dpf_avx[i-1][q*3 + p7X_M];  tv.v = tp[p7O_MI];  prv =              ((tv.p[r] > 0.0) ? u.p[r] : 0.0);
    u.v = ox->dpf_avx[i-1][q*3 + p7X_I];  tv.v = tp[p7O_II];  prv = ESL_MAX(prv, ((tv.p[r] > 0.0) ? u.p[r] : 0.0));
    u.v = ox->dpf_avx[i][q*3 + p7X_I];
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
  int     Q     = p7O_NQF_AVX(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m256 *tp    = om->tfv_avx + 7*q;
  __m256  xBv   = _mm256_set1_ps(ox->xmx[(i-1)*p7X_NXCELLS+p7X_B]);
  __m256  mpv, dpv, ipv;
  union { __m256 v; float p[8]; } u, tv;
  float   path[4];
  int     state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };

  if (q > 0) {
    mpv = ox->dpf_avx[i-1][(q-1)*3 + p7X_M];
    dpv = ox->dpf_avx[i-1][(q-1)*3 + p7X_D];
    ipv = ox->dpf_avx[i-1][(q-1)*3 + p7X_I];
  } else {
    mpv = esl_avx_rightshiftz_float(ox->dpf_avx[i-1][(Q-1)*3 + p7X_M]);
    dpv = esl_avx_rightshiftz_float(ox->dpf_avx[i-1][(Q-1)*3 + p7X_D]);
    ipv = esl_avx_rightshiftz_float(ox->dpf_avx[i-1][(Q-1)*3 + p7X_I]);
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
  int     Q     = p7O_NQF_AVX(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  union { __m256 v; float p[8]; } mpv, dpv, tmdv, tddv;
  float   path[2];

  if (q > 0) {
    mpv.v  = ox->dpf_avx[i][(q-1)*3 + p7X_M];
    dpv.v  = ox->dpf_avx[i][(q-1)*3 + p7X_D];
    tmdv.v = om->tfv_avx[7*(q-1) + p7O_MD];
    tddv.v = om->tfv_avx[7*Q + (q-1)];
  } else {
    mpv.v  = esl_avx_rightshiftz_float(ox->dpf_avx[i][(Q-1)*3 + p7X_M]);
    dpv.v  = esl_avx_rightshiftz_float(ox->dpf_avx[i][(Q-1)*3 + p7X_D]);
    tmdv.v = esl_avx_rightshiftz_float(om->tfv_avx[7*(Q-1) + p7O_MD]);
    tddv.v = esl_avx_rightshiftz_float(om->tfv_avx[8*Q-1]);
  }

  path[0] = ((tmdv.p[r] == 0.0) ? -eslINFINITY : mpv.p[r]);
  path[1] = ((tddv.p[r] == 0.0) ? -eslINFINITY : dpv.p[r]);
  return  ((path[0] >= path[1]) ? p7T_M : p7T_D);
}

static inline int
oa_select_i(const P7_OPROFILE *om, const P7_OMX *ox, int i, int k)
{
  int     Q    = p7O_NQF_AVX(ox->M);
  int     q    = (k-1) % Q;
  int     r    = (k-1) / Q;
  __m256 *tp   = om->tfv_avx + 7*q + p7O_MI;
  union { __m256 v; float p[8]; } tv, mpv, ipv;
  float   path[2];

  mpv.v = ox->dpf_avx[i-1][q*3 + p7X_M]; tv.v = *tp;  path[0] = ((tv.p[r] == 0.0) ? -eslINFINITY : mpv.p[r]);  tp++;
  ipv.v = ox->dpf_avx[i-1][q*3 + p7X_I]; tv.v = *tp;  path[1] = ((tv.p[r] == 0.0) ? -eslINFINITY : ipv.p[r]);
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
  int     Q     = p7O_NQF_AVX(ox->M);
  __m256 *dp    = ox->dpf_avx[i];
  union { __m256 v; float p[8]; } u;
  float  max   = -eslINFINITY;
  int    smax  = p7T_M;
  int    kmax  = 1;
  int    q,r;

  for (q = 0; q < Q; q++)
    {
      u.v   = *dp; dp++;  for (r = 0; r < 8; r++) if (u.p[r] >= max) { max = u.p[r]; smax = p7T_M; kmax = r*Q + q + 1; }
      u.v   = *dp; dp+=2; for (r = 0; r < 8; r++) if (u.p[r] > max)  { max = u.p[r]; smax = p7T_D; kmax = r*Q + q + 1; }
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


/* Function:  p7_OATrace_avx()
 *
 * Purpose:   AVX2 implementation of optimal accuracy traceback.
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
p7_OATrace_avx(const P7_OPROFILE *om, const P7_OMX *ox, P7_TRACE *tr)
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

/* oat_ppcodes8()
 * Posterior codes (p7OAT_PPCODE(), 0..10) of the 8 lanes of <v> as 4
 * lane masks, bit b of each lane's code in <m[b]>.
 */
static inline void
oat_ppcodes8(__m256 v, int m[4])
{
  union { __m256 v; float p[8]; } u;
  int r, c;

  u.v  = v;
  m[0] = m[1] = m[2] = m[3] = 0;
  for (r = 0; r < 8; r++)
    {
      c     = p7OAT_PPCODE(u.p[r]);
      m[0] |= ( c       & 1) << r;
      m[1] |= ((c >> 1) & 1) << r;
      m[2] |= ((c >> 2) & 1) << r;
      m[3] |= ((c >> 3) & 1) << r;
    }
}

/* oat_store_row_avx()
 * Store into <ox> what the OA trace reads from OA row <i> (<dpc>, with
 * the OA special states in <xmx>): for each k, the predecessors of an
 * M at (i+1,k), an I at (i+1,k) and a D at (i,k), chosen with the same
 * comparisons, shifts and ties as p7_OATrace_avx() makes on the OA
 * matrix; the posterior codes of M and I at (i,k) from the posterior
 * row <ppp> (all 0 for row 0, <ppp> NULL); the row's OA special states
 * and the posteriors of N, J, C (<ppx>, NULL for row 0); and the cell
 * an E at row i comes from, in p7_OATrace_avx()'s scan order.
 */
static void
oat_store_row_avx(const P7_OPROFILE *om, P7_OATMX *ox, const __m256 *dpc, const float *xmx, const __m256 *ppp, const float *ppx, int i)
{
  int      M     = om->M;
  int      Q     = p7O_NQF_AVX(M);
  __m256   zerov = _mm256_setzero_ps();
  __m256   infv  = _mm256_set1_ps(-eslINFINITY);
  __m256   xBv   = _mm256_set1_ps(xmx[i*p7X_NXCELLS+p7X_B]);
  __m256   mpv, ipv, dpv, tmdv, tddv;
  __m256   pM, pI, pD, pB, best, gt, b0, b1;
  __m256  *tp;
  union { __m256 v; float p[8]; } u;
  float    max;
  int      kmax, smax;
  int      mm[4], mi[4];
  int      q, r, c;

  for (q = 0; q < Q; q++)
    {
      tp = om->tfv_avx + 7*q;

      /* row i at k-1 */
      if (q > 0) {
        mpv  = MMO(dpc, q-1);
        ipv  = IMO(dpc, q-1);
        dpv  = DMO(dpc, q-1);
        tmdv = om->tfv_avx[7*(q-1) + p7O_MD];
        tddv = om->tfv_avx[7*Q + (q-1)];
      } else {
        mpv  = esl_avx_rightshiftz_float(MMO(dpc, Q-1));
        ipv  = esl_avx_rightshiftz_float(IMO(dpc, Q-1));
        dpv  = esl_avx_rightshiftz_float(DMO(dpc, Q-1));
        tmdv = esl_avx_rightshiftz_float(om->tfv_avx[7*(Q-1) + p7O_MD]);
        tddv = esl_avx_rightshiftz_float(om->tfv_avx[8*Q-1]);
      }

      /* M at (i+1,k): first best of M, I, D, B */
      pB   = _mm256_blendv_ps(xBv, infv, _mm256_cmp_ps(tp[p7O_BM], zerov, _CMP_EQ_OQ));
      pM   = _mm256_blendv_ps(mpv, infv, _mm256_cmp_ps(tp[p7O_MM], zerov, _CMP_EQ_OQ));
      pI   = _mm256_blendv_ps(ipv, infv, _mm256_cmp_ps(tp[p7O_IM], zerov, _CMP_EQ_OQ));
      pD   = _mm256_blendv_ps(dpv, infv, _mm256_cmp_ps(tp[p7O_DM], zerov, _CMP_EQ_OQ));
      best = pM;
      gt   = _mm256_cmp_ps(pI, best, _CMP_GT_OQ);  best = _mm256_blendv_ps(best, pI, gt);  b0 = gt;                    b1 = zerov;
      gt   = _mm256_cmp_ps(pD, best, _CMP_GT_OQ);  best = _mm256_blendv_ps(best, pD, gt);  b0 = _mm256_andnot_ps(gt, b0); b1 = _mm256_or_ps(b1, gt);
      gt   = _mm256_cmp_ps(pB, best, _CMP_GT_OQ);                                           b0 = _mm256_or_ps(b0, gt);    b1 = _mm256_or_ps(b1, gt);
      OATO_STRIPE8(ox, i, p7OAT_MP,   q) = (uint8_t) _mm256_movemask_ps(b0);
      OATO_STRIPE8(ox, i, p7OAT_MP+1, q) = (uint8_t) _mm256_movemask_ps(b1);

      /* I at (i+1,k): M on ties */
      pM = _mm256_blendv_ps(MMO(dpc, q), infv, _mm256_cmp_ps(tp[p7O_MI], zerov, _CMP_EQ_OQ));
      pI = _mm256_blendv_ps(IMO(dpc, q), infv, _mm256_cmp_ps(tp[p7O_II], zerov, _CMP_EQ_OQ));
      OATO_STRIPE8(ox, i, p7OAT_IP, q) = (uint8_t) _mm256_movemask_ps(_mm256_cmp_ps(pI, pM, _CMP_GT_OQ));

      /* D at (i,k): M on ties */
      pM = _mm256_blendv_ps(mpv, infv, _mm256_cmp_ps(tmdv, zerov, _CMP_EQ_OQ));
      pD = _mm256_blendv_ps(dpv, infv, _mm256_cmp_ps(tddv, zerov, _CMP_EQ_OQ));
      OATO_STRIPE8(ox, i, p7OAT_DP, q) = (uint8_t) _mm256_movemask_ps(_mm256_cmp_ps(pD, pM, _CMP_GT_OQ));

      /* posterior codes of M and I */
      if (ppp != NULL) { oat_ppcodes8(ppp[q*3 + p7X_M], mm); oat_ppcodes8(ppp[q*3 + p7X_I], mi); }
      else             { mm[0] = mm[1] = mm[2] = mm[3] = mi[0] = mi[1] = mi[2] = mi[3] = 0; }
      for (c = 0; c < 4; c++)
        {
          OATO_STRIPE8(ox, i, p7OAT_PPM + c, q) = (uint8_t) mm[c];
          OATO_STRIPE8(ox, i, p7OAT_PPI + c, q) = (uint8_t) mi[c];
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
      u.v = MMO(dpc, q); for (r = 0; r < 8; r++) if (u.p[r] >= max) { max = u.p[r]; smax = p7T_M; kmax = r*Q + q + 1; }
      u.v = DMO(dpc, q); for (r = 0; r < 8; r++) if (u.p[r] >  max) { max = u.p[r]; smax = p7T_D; kmax = r*Q + q + 1; }
    }
  ox->ek[i] = (smax == p7T_M) ? kmax : -kmax;
}


/* Function:  p7_OptimalAccuracySweep_avx()
 * Synopsis:  OA alignment in one forward sweep over a checkpointed Backward, AVX.
 *
 * Purpose:   For profile <om> and sequence <dsq> of length <L>, with the
 *            checkpointed Backward matrix <cx> filled by
 *            <p7_BackwardChk_avx()>, compute in one sweep through rows
 *            1..L the Forward, posterior and optimal accuracy (OA) rows,
 *            computing each block of Backward rows again from its
 *            checkpoint as the sweep reaches it, and store in <ox> what
 *            <p7_OATraceSweep_avx()> needs for the OA trace. Rows of
 *            Forward, posterior and OA values are kept only as long as
 *            the next row needs them.
 *
 *            Each row is calculated as <p7_Forward_avx()>,
 *            <p7_Decoding_avx()> and <p7_OptimalAccuracy_avx()> calculate
 *            it. The checkpointed Backward has its own scale factors, so
 *            the posteriors are those <p7_Decoding_avx()> gives for a
 *            Backward matrix with its own scale factors: normalized by
 *            the Backward score, and equal to those from a Backward
 *            scaled like Forward up to rounding.
 *
 *            If <opt_null2> is non-NULL, the null2 odds ratios of the
 *            sequence are also calculated into it, as
 *            <p7_Null2_ByExpectation_avx()> calculates them from the
 *            posteriors: the posterior rows are summed as the sweep
 *            computes them, in the same order.
 *
 *            <ox> must have been created with <p7OAT_NPLANES> planes;
 *            its layout is set to 8 lanes and it is grown for <om->M> by
 *            <L>. On return <cx> has its last block loaded.
 *
 * Args:      dsq       - sequence in digitized form, 1..L
 *            L         - length of dsq (>= 1)
 *            om        - profile
 *            cx        - checkpointed Backward matrix from p7_BackwardChk_avx()
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
 *            posterior scale correction overflows (as <p7_Decoding_avx()>).
 */
int
p7_OptimalAccuracySweep_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, P7_OATMX *ox,
                            float *opt_null2, float *opt_fsc, float *ret_oasc)
{
  int      M      = om->M;
  int      Q      = p7O_NQF_AVX(M);
  P7_OMX  *wk     = NULL;   /* rows 0,1: Forward; 2,3: OA; 4: posteriors; 5: posterior sums; xmx: Forward specials */
  float   *oxmx   = NULL;
  float    ppx[p7X_NXCELLS];
  float    cntx[p7X_NXCELLS];
  float    totscale = 0.0f;
  float    scaleproduct;
  __m256  *pp, *cnt;
  int      i, q, b;
  int      status;

  if (ox->nplanes != p7OAT_NPLANES) ESL_EXCEPTION(eslEINVAL, "P7_OATMX not set up for the standard OA sweep");
  if (cx->L != L || cx->M != M || cx->dir != p7_CHK_BACKWARD || cx->nscells != p7X_NSCELLS) ESL_EXCEPTION(eslEINVAL, "checkpointed Backward matrix doesn't match");
  ox->lanes = 8;
  if ((status = p7_oatmx_GrowTo(ox, M, L)) != eslOK) return status;

  if ((wk = p7_omx_Create_avx(M, 5, L)) == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(oxmx, sizeof(float) * (L+1) * p7X_NXCELLS);
  pp  = wk->dpf_avx[4];
  cnt = wk->dpf_avx[5];

  scaleproduct = 1.0 / cx->xmx[p7X_N];
  for (q = 0; q < p7X_NXCELLS; q++) cntx[q] = 0.0f;

  p7_ForwardRow_avx(dsq, om, wk->xmx, &totscale, wk->dpf_avx[0], NULL, 0);
  oa_row_avx(om, NULL, NULL, oxmx, wk->dpf_avx[2], NULL, 0);
  oat_store_row_avx(om, ox, wk->dpf_avx[2], oxmx, NULL, NULL, 0);

  for (i = 1; i <= L; i++)
    {
      b = (i / cx->blockR) * cx->blockR;
      if (b != cx->b && (status = p7_BackwardChkLoad_avx(dsq, om, cx, b)) != eslOK) goto ERROR;

      p7_ForwardRow_avx(dsq, om, wk->xmx, &totscale, wk->dpf_avx[i%2], wk->dpf_avx[(i-1)%2], i);
      p7_DecodingRow_avx(om, wk->xmx, cx->xmx, wk->dpf_avx[i%2], cx->dpf_avx[i], pp, ppx, i, scaleproduct);
      scaleproduct *= wk->xmx[i*p7X_NXCELLS+p7X_SCALE] / cx->xmx[i*p7X_NXCELLS+p7X_SCALE];

      oa_row_avx(om, pp, ppx, oxmx, wk->dpf_avx[2 + i%2], wk->dpf_avx[2 + (i-1)%2], i);
      oat_store_row_avx(om, ox, wk->dpf_avx[2 + i%2], oxmx, pp, ppx, i);

      /* posterior sums for null2, added up as p7_Null2_ByExpectation_avx() does */
      if (i == 1)
        {
          memcpy(cnt, pp, sizeof(__m256) * 3 * Q);
          cntx[p7X_N] = ppx[p7X_N];
          cntx[p7X_C] = ppx[p7X_C];
          cntx[p7X_J] = ppx[p7X_J];
        }
      else
        {
          for (q = 0; q < Q; q++)
            {
              cnt[q*3 + p7X_M] = _mm256_add_ps(pp[q*3 + p7X_M], cnt[q*3 + p7X_M]);
              cnt[q*3 + p7X_I] = _mm256_add_ps(pp[q*3 + p7X_I], cnt[q*3 + p7X_I]);
            }
          cntx[p7X_N] += ppx[p7X_N];
          cntx[p7X_C] += ppx[p7X_C];
          cntx[p7X_J] += ppx[p7X_J];
        }
    }
  if (isinf(scaleproduct)) { status = eslERANGE; goto ERROR; }

  if (opt_null2 != NULL) p7_Null2_FromCounts_avx(om, cnt, cntx, L, opt_null2);
  if (opt_fsc   != NULL) *opt_fsc = totscale + log(wk->xmx[L*p7X_NXCELLS+p7X_C] * om->xf[p7O_C][p7O_MOVE]);
  *ret_oasc = oxmx[L*p7X_NXCELLS+p7X_C];

  p7_omx_Destroy_avx(wk);
  free(oxmx);
  return eslOK;

 ERROR:
  if (wk) p7_omx_Destroy_avx(wk);
  free(oxmx);
  return status;
}


/* Function:  p7_OATraceSweep_avx()
 * Synopsis:  OA trace from what the AVX OA sweep stored.
 *
 * Purpose:   The OA trace for profile <om> from <ox>, filled by
 *            <p7_OptimalAccuracySweep_avx()>, into the empty trace <tr>
 *            (created with <p7_trace_CreateWithPP()>). The choices are
 *            those <p7_OATrace_avx()> makes on the same OA scores. The
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
p7_OATraceSweep_avx(const P7_OPROFILE *om, const P7_OATMX *ox, P7_TRACE *tr)
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

#endif /* eslENABLE_AVX */
