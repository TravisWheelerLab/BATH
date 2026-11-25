/* Frameshift stochastic traceback; SSE implementation.
 * Ported from impl_sse/stotrace_fs.c with _sse suffix for runtime dispatch.
 *
 * Contents:
 *   1. p7_StochasticTrace_Frameshift_sse() — main traceback
 *   2. select_*_fs() helpers
 */
#include "p7_config.h"

#ifdef eslENABLE_SSE

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>  /* SSE  */
#include <emmintrin.h>  /* SSE2 */

#include "easel.h"
#include "esl_random.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"

/* What the choice functions read of a frameshift Forward matrix: the
 * model length, the row pointers and the special states. Filled from a
 * full Forward matrix (P7_OMX) or from the loaded rows of a checkpointed
 * one (P7_CHKOMX).
 */
typedef struct {
  int            M;
  __m128 * const *dpf;
  const float   *xmx;
} STO_FS_ROWS;

/* Forward declarations of static helpers */
static inline int select_m_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, int i, int k);
static inline int select_d_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, int i, int k);
static inline int select_i_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, int i, int k);
static inline int select_n_fs(int i);
static inline int select_c_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, int i);
static inline int select_j_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, int i);
static inline int select_e_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, int i, int *ret_k);
static inline int select_b_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, int i);
static inline int select_codon_len_fs(ESL_RANDOMNESS *rng, const STO_FS_ROWS *ox, int i, int k, int Q);


/*****************************************************************
 * 1. Stochastic trace implementation
 *****************************************************************/

/* sto_fs_step()
 * One step of a stochastic trace back through the frameshift Forward
 * rows in <ox>: from state <*s0> at row <*i> and model position <*k>,
 * choose the state before it (and for M, its codon length), append it
 * to <tr>, and update <*s0>, <*i> and <*k>. Reads rows <*i> and <*i>-3
 * of <ox>.
 */
static int
sto_fs_step(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, P7_TRACE *tr, int *s0, int *i, int *k)
{
  int Q = p7O_NQF(ox->M);
  int c, s1;
  int status;

  switch (*s0) {
  case p7T_M: s1 = select_m_fs(rng, om_fs, ox, *i, *k);  (*k)--;    break;
  case p7T_D: s1 = select_d_fs(rng, om_fs, ox, *i, *k);  (*k)--;    break;
  case p7T_I: s1 = select_i_fs(rng, om_fs, ox, *i, *k);  (*i) -= 3; break;
  case p7T_N: s1 = select_n_fs(*i);                                 break;
  case p7T_C: s1 = select_c_fs(rng, om_fs, ox, *i);                 break;
  case p7T_J: s1 = select_j_fs(rng, om_fs, ox, *i);                 break;
  case p7T_E: s1 = select_e_fs(rng, om_fs, ox, *i, k);              break;
  case p7T_B: s1 = select_b_fs(rng, om_fs, ox, *i);                 break;
  default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
  }
  if (s1 == -1) ESL_EXCEPTION(eslEINVAL, "Stochastic traceback choice failed");

  /* For M states: stochastically select codon length c (1..5) from C1..C5 cells.
   * If the codon would extend before the sequence start, revert to B.
   */
  if (s1 == p7T_M) {
    c = select_codon_len_fs(rng, ox, *i, *k, Q);
    if (*i - c < 0) s1 = p7T_B;
  } else {
    c = 0;
  }

  if ((status = p7_trace_fs_Append(tr, s1, *k, *i, c)) != eslOK) return status;

  /* For NCJ states: deferred i decrement (one nucleotide per step). */
  if ((s1 == p7T_N || s1 == p7T_C || s1 == p7T_J) && s1 == *s0) (*i)--;
  *s0 = s1;
  *i -= c;
  return eslOK;
}

/* Function:  p7_StochasticTrace_Frameshift_sse()
 *
 * Purpose:   SSE implementation of stochastic backtrace from a frameshift-aware
 *            Forward matrix <ox> (8-cell FS layout).  Samples one alignment path
 *            proportional to its probability and writes it into <tr>.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> on various internal failures.
 */
int
p7_StochasticTrace_Frameshift_sse(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L,
                                   const P7_FS_OPROFILE *om_fs, const P7_OMX *ox,
                                   P7_TRACE *tr)
{
  STO_FS_ROWS fx;
  int         i = L;
  int         k = 0;
  int         s0;
  int         status;

  fx.M   = ox->M;
  fx.dpf = ox->dpf;
  fx.xmx = ox->xmx;

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i, 0)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i, 0)) != eslOK) return status;
  s0 = p7T_C;
  while (s0 != p7T_S)
    if ((status = sto_fs_step(rng, om_fs, &fx, tr, &s0, &i, &k)) != eslOK) return status;

  tr->M = om_fs->M;
  tr->L = L;
  return p7_trace_fs_Reverse(tr);
}


/* Function:  p7_StochasticTraceChk_Frameshift_sse()
 * Synopsis:  Sample tracebacks from a checkpointed frameshift Forward matrix.
 *
 * Purpose:   Sample <ntr> stochastic tracebacks of the 5-codon profile
 *            <om_fs> against <dsq> of length <L> from the checkpointed
 *            Forward matrix <cx> that <p7_ForwardChk_Frameshift_sse()>
 *            filled: trace <t> into <tr[t]>, with its own random number
 *            generator <rng[t]>. Each trace makes the same choices, in
 *            the same order, as <p7_StochasticTrace_Frameshift_sse()> on
 *            the full Forward matrix with the same generator, so the
 *            traces are identical to it.
 *
 *            All traces are sampled together. The blocks of <cx> are
 *            loaded from the last to the first, each computed again from
 *            its checkpoint by <p7_ForwardChkLoad_Frameshift_sse()> with
 *            <ov> as workspace, and every trace is taken back through the
 *            rows of the loaded block before the next block is loaded, so
 *            the Forward matrix is computed again only once for all the
 *            traces. On return block 0 of <cx> is loaded.
 *
 * Args:      rng   - <ntr> random number generators, one per trace
 *            ntr   - number of traces
 *            dsq   - digital sequence, 1..L
 *            L     - length of dsq
 *            om_fs - 5-codon frameshift profile
 *            cx    - checkpointed Forward matrix from p7_ForwardChk_Frameshift_sse()
 *            ov    - workspace for the profile's M
 *            tr    - <ntr> empty traces; RETURN: the sampled traces
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if a trace isn't
 *            empty, <cx> isn't a checkpointed frameshift Forward matrix for
 *            this comparison, or a choice fails.
 */
int
p7_StochasticTraceChk_Frameshift_sse(ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                     P7_CHKOMX *cx, P7_OIVX *ov, P7_TRACE **tr)
{
  STO_FS_ROWS fx;
  int        *s0 = NULL;
  int        *i  = NULL;
  int        *k  = NULL;
  int         b, t;
  int         status;

  if (cx->dir != p7_CHK_FORWARD || cx->nscells != p7X_NSCELLS_FS || cx->L != L || cx->M != om_fs->M)
    ESL_EXCEPTION(eslEINVAL, "checkpointed frameshift Forward matrix doesn't match");
  for (t = 0; t < ntr; t++)
    if (tr[t]->N != 0) ESL_EXCEPTION(eslEINVAL, "trace not empty; needs to be Reuse()'d?");

  ESL_ALLOC(s0, sizeof(int) * ESL_MAX(1, ntr));
  ESL_ALLOC(i,  sizeof(int) * ESL_MAX(1, ntr));
  ESL_ALLOC(k,  sizeof(int) * ESL_MAX(1, ntr));
  for (t = 0; t < ntr; t++)
    {
      i[t] = L;
      k[t] = 0;
      if ((status = p7_trace_fs_Append(tr[t], p7T_T, 0, L, 0)) != eslOK) goto ERROR;
      if ((status = p7_trace_fs_Append(tr[t], p7T_C, 0, L, 0)) != eslOK) goto ERROR;
      s0[t] = p7T_C;
    }

  fx.M   = om_fs->M;
  fx.dpf = cx->dpf;
  fx.xmx = cx->xmx;

  /* a step from row i reads rows i and i-3, which are loaded with the
   * block holding row i and the 5 rows saved before it */
  for (b = (L / cx->blockR) * cx->blockR; b >= 0; b -= cx->blockR)
    {
      if (b != cx->b && (status = p7_ForwardChkLoad_Frameshift_sse(dsq, om_fs, cx, ov, b)) != eslOK) goto ERROR;
      for (t = 0; t < ntr; t++)
        while (s0[t] != p7T_S && i[t] >= b)
          if ((status = sto_fs_step(rng[t], om_fs, &fx, tr[t], &s0[t], &i[t], &k[t])) != eslOK) goto ERROR;
    }

  for (t = 0; t < ntr; t++)
    {
      tr[t]->M = om_fs->M;
      tr[t]->L = L;
      if ((status = p7_trace_fs_Reverse(tr[t])) != eslOK) goto ERROR;
    }
  free(s0);
  free(i);
  free(k);
  return eslOK;

 ERROR:
  free(s0);
  free(i);
  free(k);
  return status;
}


/*****************************************************************
 * 2. Selection helpers
 *****************************************************************/

/* M(i,k) is reached from B(i), M(i,k-1), I(i,k-1), or D(i,k-1). Same row i. */
static inline int
select_m_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, int i, int k)
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
select_d_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, int i, int k)
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

/* I(i,k) is reached from M(i-3,k) or I(i-3,k). */
static inline int
select_i_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, int i, int k)
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
select_c_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, int i)
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
select_j_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, int i)
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
select_e_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, int i, int *ret_k)
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
  ESL_EXCEPTION(-1, "unreached code was reached");
}

/* B(i) is reached from N(i) or J(i). Same row, no scale adjustment needed. */
static inline int
select_b_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const STO_FS_ROWS *ox, int i)
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
select_codon_len_fs(ESL_RANDOMNESS *rng, const STO_FS_ROWS *ox, int i, int k, int Q)
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

#endif /* eslENABLE_SSE */
