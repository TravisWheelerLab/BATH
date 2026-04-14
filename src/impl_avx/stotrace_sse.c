/* Stochastic traceback; SSE implementations for impl_avx dispatch.
 * These are the _sse-suffixed versions called by the runtime dispatchers
 * in stotrace.c and stotrace_fs.c.
 */

#ifdef eslENABLE_SSE
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>
#include <emmintrin.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"


/*****************************************************************
 * 1. Standard stochastic traceback (p7_StochasticTrace_sse)
 *****************************************************************/

static inline int select_m(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const P7_OMX *ox, int i, int k);
static inline int select_d(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const P7_OMX *ox, int i, int k);
static inline int select_i(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const P7_OMX *ox, int i, int k);
static inline int select_n(int i);
static inline int select_c(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const P7_OMX *ox, int i);
static inline int select_j(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const P7_OMX *ox, int i);
static inline int select_e(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const P7_OMX *ox, int i, int *ret_k);
static inline int select_b(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const P7_OMX *ox, int i);


/* M(i,k) is reached from B(i-1), M(i-1,k-1), D(i-1,k-1), or I(i-1,k-1). */
static inline int
select_m(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m128 *tp    = om->tfv + 7*q;
  __m128  xBv   = _mm_set1_ps(ox->xmx[(i-1)*p7X_NXCELLS+p7X_B]);
  __m128  mpv, dpv, ipv;
  union { __m128 v; float p[4]; } u;
  float   path[4];
  int     state[4] = { p7T_B, p7T_M, p7T_I, p7T_D };

  if (q > 0) {
    mpv = ox->dpf[i-1][(q-1)*3 + p7X_M];
    dpv = ox->dpf[i-1][(q-1)*3 + p7X_D];
    ipv = ox->dpf[i-1][(q-1)*3 + p7X_I];
  } else {
    mpv = esl_sse_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_M]);
    dpv = esl_sse_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_D]);
    ipv = esl_sse_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_I]);
  }

  u.v = _mm_mul_ps(xBv, *tp); tp++;  path[0] = u.p[r];
  u.v = _mm_mul_ps(mpv, *tp); tp++;  path[1] = u.p[r];
  u.v = _mm_mul_ps(ipv, *tp); tp++;  path[2] = u.p[r];
  u.v = _mm_mul_ps(dpv, *tp);        path[3] = u.p[r];
  esl_vec_FNorm(path, 4);
  return state[esl_rnd_FChoose(rng, path, 4)];
}

/* D(i,k) is reached from M(i, k-1) or D(i,k-1). */
static inline int
select_d(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m128  mpv, dpv;
  __m128  tmdv, tddv;
  union { __m128 v; float p[4]; } u;
  float   path[2];
  int     state[2] = { p7T_M, p7T_D };

  if (q > 0) {
    mpv  = ox->dpf[i][(q-1)*3 + p7X_M];
    dpv  = ox->dpf[i][(q-1)*3 + p7X_D];
    tmdv = om->tfv[7*(q-1) + p7O_MD];
    tddv = om->tfv[7*Q + (q-1)];
  } else {
    mpv  = esl_sse_rightshiftz_float(ox->dpf[i][(Q-1)*3 + p7X_M]);
    dpv  = esl_sse_rightshiftz_float(ox->dpf[i][(Q-1)*3 + p7X_D]);
    tmdv = esl_sse_rightshiftz_float(om->tfv[7*(Q-1) + p7O_MD]);
    tddv = esl_sse_rightshiftz_float(om->tfv[8*Q-1]);
  }

  u.v = _mm_mul_ps(mpv, tmdv); path[0] = u.p[r];
  u.v = _mm_mul_ps(dpv, tddv); path[1] = u.p[r];
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}

/* I(i,k) is reached from M(i-1, k) or I(i-1,k). */
static inline int
select_i(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const P7_OMX *ox, int i, int k)
{
  int     Q    = p7O_NQF(ox->M);
  int     q    = (k-1) % Q;
  int     r    = (k-1) / Q;
  __m128  mpv  = ox->dpf[i-1][q*3 + p7X_M];
  __m128  ipv  = ox->dpf[i-1][q*3 + p7X_I];
  __m128 *tp   = om->tfv + 7*q + p7O_MI;
  union { __m128 v; float p[4]; } u;
  float   path[2];
  int     state[2] = { p7T_M, p7T_I };

  u.v = _mm_mul_ps(mpv, *tp); tp++;  path[0] = u.p[r];
  u.v = _mm_mul_ps(ipv, *tp);        path[1] = u.p[r];
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}

/* N(i) must come from N(i-1) for i>0; else it comes from S */
static inline int
select_n(int i)
{
  if (i == 0) return p7T_S;
  else        return p7T_N;
}

/* C(i) is reached from E(i) or C(i-1). */
static inline int
select_c(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const P7_OMX *ox, int i)
{
  float path[2];
  int   state[2] = { p7T_C, p7T_E };

  path[0] = ox->xmx[(i-1)*p7X_NXCELLS+p7X_C] * om->xf[p7O_C][p7O_LOOP];
  path[1] = ox->xmx[    i*p7X_NXCELLS+p7X_E] * om->xf[p7O_E][p7O_MOVE] * ox->xmx[i*p7X_NXCELLS+p7X_SCALE];
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}

/* J(i) is reached from E(i) or J(i-1). */
static inline int
select_j(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const P7_OMX *ox, int i)
{
  float path[2];
  int   state[2] = { p7T_J, p7T_E };

  path[0] = ox->xmx[(i-1)*p7X_NXCELLS+p7X_J] * om->xf[p7O_J][p7O_LOOP];
  path[1] = ox->xmx[    i*p7X_NXCELLS+p7X_E] * om->xf[p7O_E][p7O_LOOP] * ox->xmx[i*p7X_NXCELLS+p7X_SCALE];
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}

/* E(i) is reached from any M(i, k=1..M) or D(i, k=2..M).
 * Uses on-the-fly FChoose with E(i) as normalizer.
 */
static inline int
select_e(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const P7_OMX *ox, int i, int *ret_k)
{
  int    Q     = p7O_NQF(ox->M);
  double sum   = 0.0;
  double roll  = esl_random(rng);
  double norm  = 1.0 / ox->xmx[i*p7X_NXCELLS+p7X_E];
  __m128 xEv   = _mm_set1_ps(norm);
  union { __m128 v; float p[4]; } u;
  int    q,r;

  while (1) {
    for (q = 0; q < Q; q++)
      {
        u.v = _mm_mul_ps(ox->dpf[i][q*3 + p7X_M], xEv);
        for (r = 0; r < 4; r++) {
          sum += u.p[r];
          if (roll < sum) { *ret_k = r*Q + q + 1; return p7T_M;}
        }

        u.v = _mm_mul_ps(ox->dpf[i][q*3 + p7X_D], xEv);
        for (r = 0; r < 4; r++) {
          sum += u.p[r];
          if (roll < sum) { *ret_k = r*Q + q + 1; return p7T_D;}
        }
      }
    ESL_DASSERT1((sum > 0.99));
  }
  /*UNREACHED*/
  ESL_EXCEPTION(-1, "unreached code was reached. universe collapses.");
}

/* B(i) is reached from N(i) or J(i). */
static inline int
select_b(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const P7_OMX *ox, int i)
{
  float path[2];
  int   state[2] = { p7T_N, p7T_J };

  path[0] = ox->xmx[i*p7X_NXCELLS+p7X_N] * om->xf[p7O_N][p7O_MOVE];
  path[1] = ox->xmx[i*p7X_NXCELLS+p7X_J] * om->xf[p7O_J][p7O_MOVE];
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}


/* Function:  p7_StochasticTrace_sse()
 *
 * Purpose:   SSE implementation of stochastic backtrace from a Forward matrix.
 *            See p7_StochasticTrace() documentation for details.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> on various internal failures.
 */
int
p7_StochasticTrace_sse(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,
                       const P7_OMX *ox, P7_TRACE *tr)
{
  int   i;
  int   k;
  int   s0, s1;
  int   status;

  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace not empty; needs to be Reuse()'d?");

  i = L;
  k = 0;
  if ((status = p7_trace_Append(tr, p7T_T, k, i)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_C, k, i)) != eslOK) return status;
  s0 = tr->st[tr->N-1];
  while (s0 != p7T_S)
    {
      switch (s0) {
      case p7T_M: s1 = select_m(rng, om, ox, i, k);  k--; i--; break;
      case p7T_D: s1 = select_d(rng, om, ox, i, k);  k--;      break;
      case p7T_I: s1 = select_i(rng, om, ox, i, k);       i--; break;
      case p7T_N: s1 = select_n(i);                            break;
      case p7T_C: s1 = select_c(rng, om, ox, i);               break;
      case p7T_J: s1 = select_j(rng, om, ox, i);               break;
      case p7T_E: s1 = select_e(rng, om, ox, i, &k);           break;
      case p7T_B: s1 = select_b(rng, om, ox, i);               break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }
      if (s1 == -1) ESL_EXCEPTION(eslEINVAL, "Stochastic traceback choice failed");

      if ((status = p7_trace_Append(tr, s1, k, i)) != eslOK) return status;

      if ( (s1 == p7T_N || s1 == p7T_J || s1 == p7T_C) && s1 == s0) i--;
      s0 = s1;
    }

  tr->M = om->M;
  tr->L = L;
  return p7_trace_Reverse(tr);
}


/*****************************************************************
 * 2. Frameshift stochastic traceback (p7_StochasticTrace_Frameshift_sse)
 *****************************************************************/

static inline int select_m_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k);
static inline int select_d_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k);
static inline int select_i_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k);
static inline int select_n_fs(int i);
static inline int select_c_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i);
static inline int select_j_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i);
static inline int select_e_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int *ret_k);
static inline int select_b_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i);
static inline int select_codon_len_fs(ESL_RANDOMNESS *rng, const P7_OMX *ox, int i, int k, int Q);


/* M(i,k) is reached from B(i), M(i,k-1), I(i,k-1), or D(i,k-1). Same row i. */
static inline int
select_m_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m128 *tp    = om_fs->tfv + 7*q;
  __m128  xBv   = _mm_set1_ps(ox->xmx[i*p7X_NXCELLS + p7X_B]);
  __m128  mpv, dpv, ipv;
  union { __m128 v; float p[4]; } u;
  float   path[4];
  int     state[4] = { p7T_B, p7T_M, p7T_I, p7T_D };

  if (q > 0) {
    mpv = MMO_FS(ox->dpf[i], q-1, p7X_FS_C0);
    dpv = DMO_FS(ox->dpf[i], q-1);
    ipv = IMO_FS(ox->dpf[i], q-1);
  } else {
    mpv = esl_sse_rightshiftz_float(MMO_FS(ox->dpf[i], Q-1, p7X_FS_C0));
    dpv = esl_sse_rightshiftz_float(DMO_FS(ox->dpf[i], Q-1));
    ipv = esl_sse_rightshiftz_float(IMO_FS(ox->dpf[i], Q-1));
  }

  u.v = _mm_mul_ps(xBv, *tp); tp++;  path[0] = u.p[r];
  u.v = _mm_mul_ps(mpv, *tp); tp++;  path[1] = u.p[r];
  u.v = _mm_mul_ps(ipv, *tp); tp++;  path[2] = u.p[r];
  u.v = _mm_mul_ps(dpv, *tp);        path[3] = u.p[r];
  esl_vec_FNorm(path, 4);
  return state[esl_rnd_FChoose(rng, path, 4)];
}

/* D(i,k) is reached from M(i,k-1) or D(i,k-1). Same row i. */
static inline int
select_d_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m128  mpv, dpv;
  __m128  tmdv, tddv;
  union { __m128 v; float p[4]; } u;
  float   path[2];
  int     state[2] = { p7T_M, p7T_D };

  if (q > 0) {
    mpv  = MMO_FS(ox->dpf[i], q-1, p7X_FS_C0);
    dpv  = DMO_FS(ox->dpf[i], q-1);
    tmdv = om_fs->tfv[7*(q-1) + p7O_MD];
    tddv = om_fs->tfv[7*Q + (q-1)];
  } else {
    mpv  = esl_sse_rightshiftz_float(MMO_FS(ox->dpf[i], Q-1, p7X_FS_C0));
    dpv  = esl_sse_rightshiftz_float(DMO_FS(ox->dpf[i], Q-1));
    tmdv = esl_sse_rightshiftz_float(om_fs->tfv[7*(Q-1) + p7O_MD]);
    tddv = esl_sse_rightshiftz_float(om_fs->tfv[8*Q-1]);
  }

  u.v = _mm_mul_ps(mpv, tmdv); path[0] = u.p[r];
  u.v = _mm_mul_ps(dpv, tddv); path[1] = u.p[r];
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}

/* I(i,k) is reached from M(i-3,k) or I(i-3,k). Inserts emit 3-nucleotide codons. */
static inline int
select_i_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q    = p7O_NQF(ox->M);
  int     q    = (k-1) % Q;
  int     r    = (k-1) / Q;
  __m128  mpv  = MMO_FS(ox->dpf[i-3], q, p7X_FS_C0);
  __m128  ipv  = IMO_FS(ox->dpf[i-3], q);
  __m128 *tp   = om_fs->tfv + 7*q + p7O_MI;
  union { __m128 v; float p[4]; } u;
  float   path[2];
  int     state[2] = { p7T_M, p7T_I };

  u.v = _mm_mul_ps(mpv, *tp); tp++;  path[0] = u.p[r];
  u.v = _mm_mul_ps(ipv, *tp);        path[1] = u.p[r];
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}

/* N(i) must come from N(i-1) for i>0, or S for i==0. */
static inline int
select_n_fs(int i)
{
  return (i == 0) ? p7T_S : p7T_N;
}

/* C(i) is reached from E(i) or C at one of three codon offsets (i-3, i-2, i-1).
 * Scale factor adjustment uses cumscale(i-3) as common denominator.
 */
static inline int
select_c_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i)
{
  float loop = om_fs->xf[p7O_C][p7O_LOOP];
  float s_im2, s_im1, s_i;
  float path[4];

  if (i < 4) return p7T_E;

  s_im2 = ox->xmx[(i-2)*p7X_NXCELLS + p7X_SCALE];
  s_im1 = ox->xmx[(i-1)*p7X_NXCELLS + p7X_SCALE];
  s_i   = ox->xmx[ i   *p7X_NXCELLS + p7X_SCALE];

  path[0] = ox->xmx[(i-3)*p7X_NXCELLS + p7X_C] * loop;
  path[1] = ox->xmx[(i-2)*p7X_NXCELLS + p7X_C] * loop  * s_im2;
  path[2] = ox->xmx[(i-1)*p7X_NXCELLS + p7X_C] * loop  * s_im2 * s_im1;
  path[3] = ox->xmx[ i   *p7X_NXCELLS + p7X_E] * om_fs->xf[p7O_E][p7O_MOVE] * s_im2 * s_im1 * s_i;
  esl_vec_FNorm(path, 4);
  return (esl_rnd_FChoose(rng, path, 4) < 3) ? p7T_C : p7T_E;
}

/* J(i) is reached from E(i) or J at one of three codon offsets (i-3, i-2, i-1).
 * Same scale-factor adjustment logic as select_c_fs().
 */
static inline int
select_j_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i)
{
  float loop = om_fs->xf[p7O_J][p7O_LOOP];
  float s_im2, s_im1, s_i;
  float path[4];

  if (i < 4) return p7T_E;

  s_im2 = ox->xmx[(i-2)*p7X_NXCELLS + p7X_SCALE];
  s_im1 = ox->xmx[(i-1)*p7X_NXCELLS + p7X_SCALE];
  s_i   = ox->xmx[ i   *p7X_NXCELLS + p7X_SCALE];

  path[0] = ox->xmx[(i-3)*p7X_NXCELLS + p7X_J] * loop;
  path[1] = ox->xmx[(i-2)*p7X_NXCELLS + p7X_J] * loop  * s_im2;
  path[2] = ox->xmx[(i-1)*p7X_NXCELLS + p7X_J] * loop  * s_im2 * s_im1;
  path[3] = ox->xmx[ i   *p7X_NXCELLS + p7X_E] * om_fs->xf[p7O_E][p7O_LOOP] * s_im2 * s_im1 * s_i;
  esl_vec_FNorm(path, 4);
  return (esl_rnd_FChoose(rng, path, 4) < 3) ? p7T_J : p7T_E;
}

/* E(i) is reached from any M(i,k=1..M) or D(i,k=2..M).
 * Uses on-the-fly FChoose with E(i) as normalizer.
 */
static inline int
select_e_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int *ret_k)
{
  int    Q    = p7O_NQF(ox->M);
  double sum  = 0.0;
  double roll = esl_random(rng);
  double norm = 1.0 / ox->xmx[i*p7X_NXCELLS + p7X_E];
  __m128 xEv  = _mm_set1_ps((float) norm);
  union { __m128 v; float p[4]; } u;
  int    q, r;

  while (1) {
    for (q = 0; q < Q; q++)
      {
        u.v = _mm_mul_ps(MMO_FS(ox->dpf[i], q, p7X_FS_C0), xEv);
        for (r = 0; r < 4; r++) {
          sum += u.p[r];
          if (roll < sum) { *ret_k = r*Q + q + 1; return p7T_M; }
        }

        u.v = _mm_mul_ps(DMO_FS(ox->dpf[i], q), xEv);
        for (r = 0; r < 4; r++) {
          sum += u.p[r];
          if (roll < sum) { *ret_k = r*Q + q + 1; return p7T_D; }
        }
      }
    ESL_DASSERT1((sum > 0.99));
  }
  /*UNREACHED*/
  ESL_EXCEPTION(-1, "unreached code was reached. universe collapses.");
}

/* B(i) is reached from N(i) or J(i). Same row, no scale adjustment. */
static inline int
select_b_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i)
{
  float path[2];
  int   state[2] = { p7T_N, p7T_J };

  path[0] = ox->xmx[i*p7X_NXCELLS + p7X_N] * om_fs->xf[p7O_N][p7O_MOVE];
  path[1] = ox->xmx[i*p7X_NXCELLS + p7X_J] * om_fs->xf[p7O_J][p7O_MOVE];
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}

/* select_codon_len_fs()
 * For M state at (i,k): stochastically select codon length c in 1..5.
 * All cells are at the same row i, so no scale-factor adjustment needed.
 */
static inline int
select_codon_len_fs(ESL_RANDOMNESS *rng, const P7_OMX *ox, int i, int k, int Q)
{
  int   q = (k-1) % Q;
  int   r = (k-1) / Q;
  union { __m128 v; float p[4]; } u;
  float path[5];

  u.v = MMO_FS(ox->dpf[i], q, p7X_FS_C1); path[0] = u.p[r];
  u.v = MMO_FS(ox->dpf[i], q, p7X_FS_C2); path[1] = u.p[r];
  u.v = MMO_FS(ox->dpf[i], q, p7X_FS_C3); path[2] = u.p[r];
  u.v = MMO_FS(ox->dpf[i], q, p7X_FS_C4); path[3] = u.p[r];
  u.v = MMO_FS(ox->dpf[i], q, p7X_FS_C5); path[4] = u.p[r];
  esl_vec_FNorm(path, 5);
  return esl_rnd_FChoose(rng, path, 5) + 1;  /* returns 1..5 */
}


/* Function:  p7_StochasticTrace_Frameshift_sse()
 *
 * Purpose:   SSE implementation of frameshift-aware stochastic backtrace.
 *            See p7_StochasticTrace_Frameshift() documentation for details.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on various internal failures.
 */
int
p7_StochasticTrace_Frameshift_sse(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L,
                                   const P7_FS_OPROFILE *om_fs, const P7_OMX *ox,
                                   P7_TRACE *tr)
{
  int Q = p7O_NQF(ox->M);
  int i = L;
  int k = 0;
  int c = 0;
  int s0, s1;
  int status;

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i, c)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i, c)) != eslOK) return status;
  s0 = p7T_C;

  while (s0 != p7T_S)
    {
      switch (s0) {
      case p7T_M: s1 = select_m_fs(rng, om_fs, ox, i, k);  k--;    break;
      case p7T_D: s1 = select_d_fs(rng, om_fs, ox, i, k);  k--;    break;
      case p7T_I: s1 = select_i_fs(rng, om_fs, ox, i, k);  i -= 3; break;
      case p7T_N: s1 = select_n_fs(i);                              break;
      case p7T_C: s1 = select_c_fs(rng, om_fs, ox, i);             break;
      case p7T_J: s1 = select_j_fs(rng, om_fs, ox, i);             break;
      case p7T_E: s1 = select_e_fs(rng, om_fs, ox, i, &k);         break;
      case p7T_B: s1 = select_b_fs(rng, om_fs, ox, i);             break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }
      if (s1 == -1) ESL_EXCEPTION(eslEINVAL, "Stochastic traceback choice failed");

      /* For M states: stochastically select codon length c (1..5). */
      if (s1 == p7T_M) {
        c = select_codon_len_fs(rng, ox, i, k, Q);
        if (i - c < 0) s1 = p7T_B;
      } else {
        c = 0;
      }

      if ((status = p7_trace_fs_Append(tr, s1, k, i, c)) != eslOK) return status;

      /* For NCJ states: deferred i decrement. */
      if ((s1 == p7T_N || s1 == p7T_C || s1 == p7T_J) && s1 == s0) i--;

      s0 = s1;
      i -= c;
    }

  tr->M = om_fs->M;
  tr->L = L;
  return p7_trace_fs_Reverse(tr);
}

#endif /* eslENABLE_SSE */
