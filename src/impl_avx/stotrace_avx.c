/* Stochastic traceback; AVX2 implementation for impl_avx dispatch.
 * This is the _avx-suffixed version called by the runtime dispatcher
 * in stotrace.c.
 */

#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <stdio.h>
#include <math.h>

#include <immintrin.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_avx.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"


/* What the choice functions read of a Forward matrix: the model length,
 * the row pointers and the special states. Filled from a full Forward
 * matrix (P7_OMX) or from the loaded rows of a checkpointed one
 * (P7_CHKOMX).
 */
typedef struct {
  int            M;
  __m256 * const *dpf;
  const float   *xmx;
} STO_ROWS;

static inline int select_m(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *ox, int i, int k);
static inline int select_d(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *ox, int i, int k);
static inline int select_i(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *ox, int i, int k);
static inline int select_n(int i);
static inline int select_c(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *ox, int i);
static inline int select_j(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *ox, int i);
static inline int select_e(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *ox, int i, int *ret_k);
static inline int select_b(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *ox, int i);


/* M(i,k) is reached from B(i-1), M(i-1,k-1), D(i-1,k-1), or I(i-1,k-1). */
static inline int
select_m(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *ox, int i, int k)
{
  int     Q     = p7O_NQF_AVX(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m256 *tp    = om->tfv_avx + 7*q;
  __m256  xBv   = _mm256_set1_ps(ox->xmx[(i-1)*p7X_NXCELLS+p7X_B]);
  __m256  mpv, dpv, ipv;
  union { __m256 v; float p[8]; } u;
  float   path[4];
  int     state[4] = { p7T_B, p7T_M, p7T_I, p7T_D };

  if (q > 0) {
    mpv = ox->dpf[i-1][(q-1)*3 + p7X_M];
    dpv = ox->dpf[i-1][(q-1)*3 + p7X_D];
    ipv = ox->dpf[i-1][(q-1)*3 + p7X_I];
  } else {
    mpv = esl_avx_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_M]);
    dpv = esl_avx_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_D]);
    ipv = esl_avx_rightshiftz_float(ox->dpf[i-1][(Q-1)*3 + p7X_I]);
  }

  u.v = _mm256_mul_ps(xBv, *tp); tp++;  path[0] = u.p[r];
  u.v = _mm256_mul_ps(mpv, *tp); tp++;  path[1] = u.p[r];
  u.v = _mm256_mul_ps(ipv, *tp); tp++;  path[2] = u.p[r];
  u.v = _mm256_mul_ps(dpv, *tp);        path[3] = u.p[r];
  esl_vec_FNorm(path, 4);
  return state[esl_rnd_FChoose(rng, path, 4)];
}

/* D(i,k) is reached from M(i, k-1) or D(i,k-1). */
static inline int
select_d(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *ox, int i, int k)
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
    mpv  = ox->dpf[i][(q-1)*3 + p7X_M];
    dpv  = ox->dpf[i][(q-1)*3 + p7X_D];
    tmdv = om->tfv_avx[7*(q-1) + p7O_MD];
    tddv = om->tfv_avx[7*Q + (q-1)];
  } else {
    mpv  = esl_avx_rightshiftz_float(ox->dpf[i][(Q-1)*3 + p7X_M]);
    dpv  = esl_avx_rightshiftz_float(ox->dpf[i][(Q-1)*3 + p7X_D]);
    tmdv = esl_avx_rightshiftz_float(om->tfv_avx[7*(Q-1) + p7O_MD]);
    tddv = esl_avx_rightshiftz_float(om->tfv_avx[8*Q-1]);
  }

  u.v = _mm256_mul_ps(mpv, tmdv); path[0] = u.p[r];
  u.v = _mm256_mul_ps(dpv, tddv); path[1] = u.p[r];
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}

/* I(i,k) is reached from M(i-1, k) or I(i-1,k). */
static inline int
select_i(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *ox, int i, int k)
{
  int     Q    = p7O_NQF_AVX(ox->M);
  int     q    = (k-1) % Q;
  int     r    = (k-1) / Q;
  __m256  mpv  = ox->dpf[i-1][q*3 + p7X_M];
  __m256  ipv  = ox->dpf[i-1][q*3 + p7X_I];
  __m256 *tp   = om->tfv_avx + 7*q + p7O_MI;
  union { __m256 v; float p[8]; } u;
  float   path[2];
  int     state[2] = { p7T_M, p7T_I };

  u.v = _mm256_mul_ps(mpv, *tp); tp++;  path[0] = u.p[r];
  u.v = _mm256_mul_ps(ipv, *tp);        path[1] = u.p[r];
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
select_c(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *ox, int i)
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
select_j(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *ox, int i)
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
select_e(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *ox, int i, int *ret_k)
{
  int    Q     = p7O_NQF_AVX(ox->M);
  double sum   = 0.0;
  double roll  = esl_random(rng);
  double norm  = 1.0 / ox->xmx[i*p7X_NXCELLS+p7X_E];
  __m256 xEv   = _mm256_set1_ps(norm);
  union { __m256 v; float p[8]; } u;
  int    q,r;

  while (1) {
    for (q = 0; q < Q; q++)
      {
        u.v = _mm256_mul_ps(ox->dpf[i][q*3 + p7X_M], xEv);
        for (r = 0; r < 8; r++) {
          sum += u.p[r];
          if (roll < sum) { *ret_k = r*Q + q + 1; return p7T_M;}
        }

        u.v = _mm256_mul_ps(ox->dpf[i][q*3 + p7X_D], xEv);
        for (r = 0; r < 8; r++) {
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
select_b(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *ox, int i)
{
  float path[2];
  int   state[2] = { p7T_N, p7T_J };

  path[0] = ox->xmx[i*p7X_NXCELLS+p7X_N] * om->xf[p7O_N][p7O_MOVE];
  path[1] = ox->xmx[i*p7X_NXCELLS+p7X_J] * om->xf[p7O_J][p7O_MOVE];
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}


/* sto_step()
 * One step of a stochastic trace back through the Forward rows in <fx>:
 * from state <*s0> at row <*i> and model position <*k>, choose the state
 * before it, append that state to <tr>, and update <*s0>, <*i> and <*k>.
 * Reads rows <*i> and <*i>-1 of <fx>.
 */
static int
sto_step(ESL_RANDOMNESS *rng, const P7_OPROFILE *om, const STO_ROWS *fx, P7_TRACE *tr, int *s0, int *i, int *k)
{
  int s1;
  int status;

  switch (*s0) {
  case p7T_M: s1 = select_m(rng, om, fx, *i, *k);  (*k)--; (*i)--; break;
  case p7T_D: s1 = select_d(rng, om, fx, *i, *k);  (*k)--;         break;
  case p7T_I: s1 = select_i(rng, om, fx, *i, *k);          (*i)--; break;
  case p7T_N: s1 = select_n(*i);                                   break;
  case p7T_C: s1 = select_c(rng, om, fx, *i);                      break;
  case p7T_J: s1 = select_j(rng, om, fx, *i);                      break;
  case p7T_E: s1 = select_e(rng, om, fx, *i, k);                   break;
  case p7T_B: s1 = select_b(rng, om, fx, *i);                      break;
  default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
  }
  if (s1 == -1) ESL_EXCEPTION(eslEINVAL, "Stochastic traceback choice failed");

  if ((status = p7_trace_Append(tr, s1, *k, *i)) != eslOK) return status;

  if ( (s1 == p7T_N || s1 == p7T_J || s1 == p7T_C) && s1 == *s0) (*i)--;
  *s0 = s1;
  return eslOK;
}

/* Function:  p7_StochasticTrace_avx()
 *
 * Purpose:   AVX2 implementation of stochastic backtrace from a Forward matrix.
 *            See p7_StochasticTrace() documentation for details.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> on various internal failures.
 */
int
p7_StochasticTrace_avx(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,
                        const P7_OMX *ox, P7_TRACE *tr)
{
  int   i;
  int   k;
  int   s0;
  STO_ROWS fx;
  int   status;

  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace not empty; needs to be Reuse()'d?");

  fx.M   = ox->M;
  fx.dpf = ox->dpf_avx;
  fx.xmx = ox->xmx;

  i = L;
  k = 0;
  if ((status = p7_trace_Append(tr, p7T_T, k, i)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_C, k, i)) != eslOK) return status;
  s0 = tr->st[tr->N-1];
  while (s0 != p7T_S)
    if ((status = sto_step(rng, om, &fx, tr, &s0, &i, &k)) != eslOK) return status;

  tr->M = om->M;
  tr->L = L;
  return p7_trace_Reverse(tr);
}

/* Function:  p7_StochasticTraceChk_avx()
 * Synopsis:  Sample tracebacks from a checkpointed Forward matrix.
 *
 * Purpose:   Sample <ntr> stochastic tracebacks of profile <om> against
 *            <dsq> of length <L> from the checkpointed Forward matrix
 *            <cx> that <p7_ForwardChk_avx()> filled: trace <t> into
 *            <tr[t]>, with its own random number generator <rng[t]>.
 *            Each trace makes the same choices, in the same order, as
 *            <p7_StochasticTrace_avx()> on the full Forward matrix with
 *            the same generator, so the traces are identical to it.
 *
 *            All traces are sampled together. The blocks of <cx> are
 *            loaded from the last to the first, each computed again from
 *            its checkpoint by <p7_ForwardChkLoad_avx()>, and every trace
 *            is taken back through the rows of the loaded block before
 *            the next block is loaded, so the Forward matrix is computed
 *            again only once for all the traces. On return block 0 of
 *            <cx> is loaded.
 *
 * Args:      rng - <ntr> random number generators, one per trace
 *            ntr - number of traces
 *            dsq - digital sequence, 1..L
 *            L   - length of dsq
 *            om  - profile
 *            cx  - checkpointed Forward matrix from p7_ForwardChk_avx()
 *            tr  - <ntr> empty traces; RETURN: the sampled traces
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if a trace isn't
 *            empty, <cx> isn't a checkpointed Forward matrix for this
 *            comparison, or a choice fails.
 */
int
p7_StochasticTraceChk_avx(ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, P7_TRACE **tr)
{
  STO_ROWS fx;
  int     *s0 = NULL;
  int     *i  = NULL;
  int     *k  = NULL;
  int      b, t;
  int      status;

  if (cx->dir != p7_CHK_FORWARD || cx->L != L || cx->M != om->M) ESL_EXCEPTION(eslEINVAL, "checkpointed Forward matrix doesn't match");
  for (t = 0; t < ntr; t++)
    if (tr[t]->N != 0) ESL_EXCEPTION(eslEINVAL, "trace not empty; needs to be Reuse()'d?");

  ESL_ALLOC(s0, sizeof(int) * ESL_MAX(1, ntr));
  ESL_ALLOC(i,  sizeof(int) * ESL_MAX(1, ntr));
  ESL_ALLOC(k,  sizeof(int) * ESL_MAX(1, ntr));
  for (t = 0; t < ntr; t++)
    {
      i[t] = L;
      k[t] = 0;
      if ((status = p7_trace_Append(tr[t], p7T_T, 0, L)) != eslOK) goto ERROR;
      if ((status = p7_trace_Append(tr[t], p7T_C, 0, L)) != eslOK) goto ERROR;
      s0[t] = p7T_C;
    }

  fx.M   = om->M;
  fx.dpf = cx->dpf_avx;
  fx.xmx = cx->xmx;

  /* a step from row i reads rows i and i-1, which are loaded with the
   * block holding row i and the row saved before it */
  for (b = (L / cx->blockR) * cx->blockR; b >= 0; b -= cx->blockR)
    {
      if (b != cx->b && (status = p7_ForwardChkLoad_avx(dsq, om, cx, b)) != eslOK) goto ERROR;
      for (t = 0; t < ntr; t++)
        while (s0[t] != p7T_S && i[t] >= b)
          if ((status = sto_step(rng[t], om, &fx, tr[t], &s0[t], &i[t], &k[t])) != eslOK) goto ERROR;
    }

  for (t = 0; t < ntr; t++)
    {
      tr[t]->M = om->M;
      tr[t]->L = L;
      if ((status = p7_trace_Reverse(tr[t])) != eslOK) goto ERROR;
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


#endif /* eslENABLE_AVX */
