/* Frameshift stochastic traceback; AVX2 implementation.
 * Ported from stotrace_fs_sse.c with _avx suffix for runtime dispatch.
 *
 * Contents:
 *   1. p7_StochasticTrace_Frameshift_avx() — main traceback
 *   2. select_*_fs() helpers
 */
#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <stdio.h>
#include <math.h>

#include <immintrin.h>  /* AVX2 */

#include "easel.h"
#include "esl_random.h"
#include "esl_avx.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"

/* Forward declarations of static helpers */
static inline int select_m_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k);
static inline int select_d_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k);
static inline int select_i_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k);
static inline int select_n_fs(int i);
static inline int select_c_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i);
static inline int select_j_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i);
static inline int select_e_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int *ret_k);
static inline int select_b_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i);
static inline int select_codon_len_fs(ESL_RANDOMNESS *rng, const P7_OMX *ox, int i, int k, int Q);


/*****************************************************************
 * 1. Stochastic trace implementation
 *****************************************************************/

/* Function:  p7_StochasticTrace_Frameshift_avx()
 *
 * Purpose:   AVX2 implementation of stochastic backtrace from a frameshift-aware
 *            Forward matrix <ox> (8-cell FS layout).  Samples one alignment path
 *            proportional to its probability and writes it into <tr>.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on various internal failures.
 */
int
p7_StochasticTrace_Frameshift_avx(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L,
                                    const P7_FS_OPROFILE *om_fs, const P7_OMX *ox,
                                    P7_TRACE *tr)
{
  int Q = p7O_NQF_AVX(ox->M);
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

      if (s1 == p7T_M) {
        c = select_codon_len_fs(rng, ox, i, k, Q);
        if (i - c < 0) s1 = p7T_B;
      } else {
        c = 0;
      }

      if ((status = p7_trace_fs_Append(tr, s1, k, i, c)) != eslOK) return status;

      if ((s1 == p7T_N || s1 == p7T_C || s1 == p7T_J) && s1 == s0) i--;

      s0 = s1;
      i -= c;
    }

  tr->M = om_fs->M;
  tr->L = L;
  return p7_trace_fs_Reverse(tr);
}


/*****************************************************************
 * 2. Selection helpers
 *****************************************************************/

/* M(i,k) is reached from B(i), M(i,k-1), I(i,k-1), or D(i,k-1). Same row i. */
static inline int
select_m_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF_AVX(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m256 *tp    = om_fs->tfv_avx + 7*q;
  __m256  xBv   = _mm256_set1_ps(ox->xmx[i*p7X_NXCELLS + p7X_B]);
  __m256  mpv, dpv, ipv;
  union { __m256 v; float p[8]; } u;
  float   path[4];
  int     state[4] = { p7T_B, p7T_M, p7T_I, p7T_D };

  if (q > 0) {
    mpv = MMO_FS(ox->dpf_avx[i], q-1, p7X_FS_C0);
    dpv = DMO_FS(ox->dpf_avx[i], q-1);
    ipv = IMO_FS(ox->dpf_avx[i], q-1);
  } else {
    mpv = esl_avx_rightshiftz_float(MMO_FS(ox->dpf_avx[i], Q-1, p7X_FS_C0));
    dpv = esl_avx_rightshiftz_float(DMO_FS(ox->dpf_avx[i], Q-1));
    ipv = esl_avx_rightshiftz_float(IMO_FS(ox->dpf_avx[i], Q-1));
  }

  u.v = _mm256_mul_ps(xBv, *tp); tp++;  path[0] = u.p[r];
  u.v = _mm256_mul_ps(mpv, *tp); tp++;  path[1] = u.p[r];
  u.v = _mm256_mul_ps(ipv, *tp); tp++;  path[2] = u.p[r];
  u.v = _mm256_mul_ps(dpv, *tp);        path[3] = u.p[r];
  esl_vec_FNorm(path, 4);
  return state[esl_rnd_FChoose(rng, path, 4)];
}

/* D(i,k) is reached from M(i,k-1) or D(i,k-1). Same row i. */
static inline int
select_d_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF_AVX(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m256  mpv, dpv;
  __m256  tmdv, tddv;
  union { __m256 v; float p[8]; } u;
  float   path[2];
  int     state[2] = { p7T_M, p7T_D };

  if (q > 0) {
    mpv  = MMO_FS(ox->dpf_avx[i], q-1, p7X_FS_C0);
    dpv  = DMO_FS(ox->dpf_avx[i], q-1);
    tmdv = om_fs->tfv_avx[7*(q-1) + p7O_MD];
    tddv = om_fs->tfv_avx[7*Q + (q-1)];
  } else {
    mpv  = esl_avx_rightshiftz_float(MMO_FS(ox->dpf_avx[i], Q-1, p7X_FS_C0));
    dpv  = esl_avx_rightshiftz_float(DMO_FS(ox->dpf_avx[i], Q-1));
    tmdv = esl_avx_rightshiftz_float(om_fs->tfv_avx[7*(Q-1) + p7O_MD]);
    tddv = esl_avx_rightshiftz_float(om_fs->tfv_avx[8*Q-1]);
  }

  u.v = _mm256_mul_ps(mpv, tmdv); path[0] = u.p[r];
  u.v = _mm256_mul_ps(dpv, tddv); path[1] = u.p[r];
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}

/* I(i,k) is reached from M(i-3,k) or I(i-3,k). */
static inline int
select_i_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q    = p7O_NQF_AVX(ox->M);
  int     q    = (k-1) % Q;
  int     r    = (k-1) / Q;
  __m256  mpv  = MMO_FS(ox->dpf_avx[i-3], q, p7X_FS_C0);
  __m256  ipv  = IMO_FS(ox->dpf_avx[i-3], q);
  __m256 *tp   = om_fs->tfv_avx + 7*q + p7O_MI;
  union { __m256 v; float p[8]; } u;
  float   path[2];
  int     state[2] = { p7T_M, p7T_I };

  u.v = _mm256_mul_ps(mpv, *tp); tp++;  path[0] = u.p[r];
  u.v = _mm256_mul_ps(ipv, *tp);        path[1] = u.p[r];
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}

/* N(i) comes from N(i-1) for i>0, or S for i==0. */
static inline int
select_n_fs(int i)
{
  return (i == 0) ? p7T_S : p7T_N;
}

/* C(i) is reached from E(i) or C at codon offsets i-3, i-2, i-1.
 * Scale-factor adjustment uses cumscale(i-3) as common denominator.
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

/* J(i) is reached from E(i) or J at codon offsets i-3, i-2, i-1. */
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
 * Uses on-the-fly FChoose algorithm with E(i) as normalizer.
 */
static inline int
select_e_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int *ret_k)
{
  int    Q    = p7O_NQF_AVX(ox->M);
  double sum  = 0.0;
  double roll = esl_random(rng);
  double norm = 1.0 / ox->xmx[i*p7X_NXCELLS + p7X_E];
  __m256 xEv  = _mm256_set1_ps((float) norm);
  union { __m256 v; float p[8]; } u;
  int    q, r;

  while (1) {
    for (q = 0; q < Q; q++)
      {
        u.v = _mm256_mul_ps(MMO_FS(ox->dpf_avx[i], q, p7X_FS_C0), xEv);
        for (r = 0; r < 8; r++) {
          sum += u.p[r];
          if (roll < sum) { *ret_k = r*Q + q + 1; return p7T_M; }
        }

        u.v = _mm256_mul_ps(DMO_FS(ox->dpf_avx[i], q), xEv);
        for (r = 0; r < 8; r++) {
          sum += u.p[r];
          if (roll < sum) { *ret_k = r*Q + q + 1; return p7T_D; }
        }
      }
    ESL_DASSERT1((sum > 0.99));
  }
  /*UNREACHED*/
  ESL_EXCEPTION(-1, "unreached code was reached");
}

/* B(i) is reached from N(i) or J(i). Same row, no scale adjustment needed. */
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

/* select_codon_len_fs():
 * Stochastically select codon length c in 1..5 from C1..C5 cells at (i,k).
 */
static inline int
select_codon_len_fs(ESL_RANDOMNESS *rng, const P7_OMX *ox, int i, int k, int Q)
{
  int   q = (k-1) % Q;
  int   r = (k-1) / Q;
  union { __m256 v; float p[8]; } u;
  float path[5];

  u.v = MMO_FS(ox->dpf_avx[i], q, p7X_FS_C1); path[0] = u.p[r];
  u.v = MMO_FS(ox->dpf_avx[i], q, p7X_FS_C2); path[1] = u.p[r];
  u.v = MMO_FS(ox->dpf_avx[i], q, p7X_FS_C3); path[2] = u.p[r];
  u.v = MMO_FS(ox->dpf_avx[i], q, p7X_FS_C4); path[3] = u.p[r];
  u.v = MMO_FS(ox->dpf_avx[i], q, p7X_FS_C5); path[4] = u.p[r];
  esl_vec_FNorm(path, 5);
  return esl_rnd_FChoose(rng, path, 5) + 1;  /* returns 1..5 */
}

#endif /* eslENABLE_AVX */
