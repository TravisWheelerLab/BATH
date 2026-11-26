/* Frameshift stochastic traceback; AVX2 implementation.
 * Ported from stotrace_fs_sse.c with _avx suffix for runtime dispatch.
 *
 * Contents:
 *   1. p7_StochasticTraceChk_Frameshift_avx() — traceback from a checkpointed Forward matrix
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
static inline void select_m_joint_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, const ESL_DSQ *dsq, int i, int k, int *ret_c, int *ret_prv);
static inline int select_d_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int i, int k);
static inline int select_i_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int i, int k);
static inline int select_n_fs(int i);
static inline int select_c_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int i);
static inline int select_c_end_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int L);
static inline int select_j_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int i);
static inline int select_e_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int i, int *ret_k);
static inline int select_b_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int i);


/*****************************************************************
 * 1. Stochastic trace implementation
 *****************************************************************/

/* sto_fs_step()
 * One step of a stochastic trace back through the loaded 3-cell rows
 * (M, D, I) of the checkpointed frameshift Forward matrix <cx> for the
 * sequence <dsq>: from state <*s0> at row <*i> and
 * model position <*k>, choose the state before it, append it to <tr>,
 * and update <*s0>, <*i> and <*k>. An M state is appended with codon
 * length 0; the step from it chooses its codon length and the state
 * before it together (select_m_joint_fs()) and sets the codon length in
 * <tr>. Reads rows <*i>-5 .. <*i> of <cx>.
 */
static int
sto_fs_step(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, const ESL_DSQ *dsq, P7_TRACE *tr, int *s0, int *i, int *k)
{
  int c, s1, z;
  int status;

  switch (*s0) {
  case p7T_M: select_m_joint_fs(rng, om_fs, cx, dsq, *i, *k, &c, &s1);  tr->c[tr->N-1] = c;  (*i) -= c;  (*k)--;  break;
  case p7T_D: s1 = select_d_fs(rng, om_fs, cx, *i, *k);  (*k)--;    break;
  case p7T_I: s1 = select_i_fs(rng, om_fs, cx, *i, *k);  (*i) -= 3; break;
  case p7T_N: s1 = select_n_fs(*i);                                 break;
  case p7T_C:   /* first C step: the row the alignment ends at; a C loop appends its three nucleotides */
    if (tr->N == 2)
      for (z = select_c_end_fs(rng, om_fs, cx, *i); z > 0; z--)
        { if ((status = p7_trace_fs_Append(tr, p7T_C, *k, *i, 0)) != eslOK) return status; (*i)--; }
    s1 = select_c_fs(rng, om_fs, cx, *i);
    if (s1 == p7T_C)
      for (z = 0; z < 2; z++)
        { if ((status = p7_trace_fs_Append(tr, p7T_C, *k, *i, 0)) != eslOK) return status; (*i)--; }
    break;
  case p7T_J:   /* a J loop appends its three nucleotides */
    s1 = select_j_fs(rng, om_fs, cx, *i);
    if (s1 == p7T_J)
      for (z = 0; z < 2; z++)
        { if ((status = p7_trace_fs_Append(tr, p7T_J, *k, *i, 0)) != eslOK) return status; (*i)--; }
    break;
  case p7T_E: s1 = select_e_fs(rng, om_fs, cx, *i, k);              break;
  case p7T_B: s1 = select_b_fs(rng, om_fs, cx, *i);                 break;
  default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
  }
  if (s1 == -1) ESL_EXCEPTION(eslEINVAL, "Stochastic traceback choice failed");

  if ((status = p7_trace_fs_Append(tr, s1, *k, *i, 0)) != eslOK) return status;

  /* For NCJ states: deferred i decrement (one nucleotide per step). */
  if ((s1 == p7T_N || s1 == p7T_C || s1 == p7T_J) && s1 == *s0) (*i)--;
  *s0 = s1;
  return eslOK;
}

/* Function:  p7_StochasticTraceChk_Frameshift_avx()
 * Synopsis:  Sample tracebacks from a checkpointed frameshift Forward matrix.
 *
 * Purpose:   Sample <ntr> stochastic tracebacks of the 5-codon profile
 *            <om_fs> against <dsq> of length <L> from the checkpointed
 *            Forward matrix <cx> that <p7_ForwardChk_Frameshift_avx()>
 *            filled: trace <t> into <tr[t]>, with its own random number
 *            generator <rng[t]>.
 *
 *            All traces are sampled together. The blocks of <cx> are
 *            loaded from the last to the first, each computed again from
 *            its checkpoint by <p7_ForwardChkLoad_Frameshift_avx()> with
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
 *            cx    - checkpointed Forward matrix from p7_ForwardChk_Frameshift_avx()
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
p7_StochasticTraceChk_Frameshift_avx(ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
                                     P7_CHKOMX *cx, P7_OIVX *ov, P7_TRACE **tr)
{
  int        *s0 = NULL;
  int        *i  = NULL;
  int        *k  = NULL;
  int         b, t;
  int         status;

  if (cx->dir != p7_CHK_FORWARD || cx->nscells != p7X_NSCELLS || cx->L != L || cx->M != om_fs->M)
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

  /* a step from row i reads rows i-5..i, which are loaded with the
   * block holding row i and the 5 rows saved before it */
  for (b = (L / cx->blockR) * cx->blockR; b >= 0; b -= cx->blockR)
    {
      if (b != cx->b && (status = p7_ForwardChkLoad_Frameshift_avx(dsq, om_fs, cx, ov, b)) != eslOK) goto ERROR;
      for (t = 0; t < ntr; t++)
        while (s0[t] != p7T_S && i[t] >= b)
          if ((status = sto_fs_step(rng[t], om_fs, cx, dsq, tr[t], &s0[t], &i[t], &k[t])) != eslOK) goto ERROR;
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

/* select_m_joint_fs()
 * For M at row <i>, model position <k>: choose its codon length c and the
 * state before it together, from the 20 terms the Forward sums for this
 * cell: X(i-c, k-1) * t(X->M_k) * e_c(k), for X in B, M, I, D and c = 1..5
 * (c <= i), where e_c is the emission of the codon of length c that ends
 * at <i> in <dsq>. Each term is taken from row i-c of <cx> and put on a
 * common scale with the scale factors of the rows between. Sets <*ret_c>
 * and <*ret_prv>.
 */
static inline void
select_m_joint_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, const ESL_DSQ *dsq, int i, int k, int *ret_c, int *ret_prv)
{
  int     Q     = p7O_NQF_AVX(cx->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  int     st[4] = { p7T_B, p7T_M, p7T_I, p7T_D };
  __m256  *tp    = om_fs->tfv_avx + 7*q;   /* transitions into M_k: BM, MM, IM, DM */
  __m256   mpv, dpv, ipv, xBv;
  union { __m256 v; float p[8]; } u;
  double  w[20], lsc[6];
  double  e, mx;
  float   tw[4];
  int     cidx[6];
  int     t, x1, x2, x3, x4, x5;
  int     c, p, z;

  x1 = (            dsq[i]   < p7P_MAXNUC) ? dsq[i]   : p7P_MAXCODONS5;
  x2 = (i >= 2 && dsq[i-1] < p7P_MAXNUC) ? dsq[i-1] : p7P_MAXCODONS5;
  x3 = (i >= 3 && dsq[i-2] < p7P_MAXNUC) ? dsq[i-2] : p7P_MAXCODONS5;
  x4 = (i >= 4 && dsq[i-3] < p7P_MAXNUC) ? dsq[i-3] : p7P_MAXCODONS5;
  x5 = (i >= 5 && dsq[i-4] < p7P_MAXNUC) ? dsq[i-4] : p7P_MAXCODONS5;
  cidx[1] = p7P_MINIDX(p7P_CODON1_FS5(x1),                 p7P_DEGEN5_QC2);
  cidx[2] = p7P_MINIDX(p7P_CODON2_FS5(x2, x1),             p7P_DEGEN5_QC1);
  cidx[3] = p7P_MINIDX(p7P_CODON3_FS5(x3, x2, x1),         p7P_DEGEN5_C);
  cidx[4] = p7P_MINIDX(p7P_CODON4_FS5(x4, x3, x2, x1),     p7P_DEGEN5_QC1);
  cidx[5] = p7P_MINIDX(p7P_CODON5_FS5(x5, x4, x3, x2, x1), p7P_DEGEN5_QC2);

  /* lsc[c]: sum of the log scale factors of rows i-4 .. i-c, which puts a
   * value of row i-c on the scale of row i-5 */
  lsc[5] = 0.;
  for (c = 4; c >= 1; c--)
    lsc[c] = lsc[c+1] + ((i-c >= 0) ? log(cx->xmx[(i-c)*p7X_NXCELLS + p7X_SCALE]) : 0.);

  for (c = 1; c <= 5; c++)
    {
      p = i - c;
      for (z = 0; z < 4; z++) w[(c-1)*4 + z] = -eslINFINITY;
      if (p < 0) continue;
      xBv = _mm256_set1_ps(cx->xmx[p*p7X_NXCELLS + p7X_B]);
      if (q > 0) {
        mpv = MMO(cx->dpf_avx[p], q-1);
        dpv = DMO(cx->dpf_avx[p], q-1);
        ipv = IMO(cx->dpf_avx[p], q-1);
      } else {
        mpv = esl_avx_rightshiftz_float(MMO(cx->dpf_avx[p], Q-1));
        dpv = esl_avx_rightshiftz_float(DMO(cx->dpf_avx[p], Q-1));
        ipv = esl_avx_rightshiftz_float(IMO(cx->dpf_avx[p], Q-1));
      }
      u.v = _mm256_mul_ps(xBv, tp[0]); tw[0] = u.p[r];
      u.v = _mm256_mul_ps(mpv, tp[1]); tw[1] = u.p[r];
      u.v = _mm256_mul_ps(ipv, tp[2]); tw[2] = u.p[r];
      u.v = _mm256_mul_ps(dpv, tp[3]); tw[3] = u.p[r];
      u.v = om_fs->rfv_avx[cidx[c]][q];
      e   = u.p[r];
      for (z = 0; z < 4; z++)
        if (tw[z] > 0. && e > 0.) w[(c-1)*4 + z] = log(tw[z]) + log(e) + lsc[c];
    }

  mx = -eslINFINITY;
  for (z = 0; z < 20; z++) if (w[z] > mx) mx = w[z];
  for (z = 0; z < 20; z++) w[z] = (w[z] == -eslINFINITY) ? 0. : exp(w[z] - mx);
  esl_vec_DNorm(w, 20);
  z = esl_rnd_DChoose(rng, w, 20);
  *ret_c   = z / 4 + 1;
  *ret_prv = st[z % 4];
}


/* D(i,k) is reached from M(i,k-1) or D(i,k-1). Same row i. */
static inline int
select_d_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int i, int k)
{
  int     Q     = p7O_NQF_AVX(cx->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m256  mpv, dpv;
  __m256  tmdv, tddv;
  union { __m256 v; float p[8]; } u;
  float   path[2];
  int     state[2] = { p7T_M, p7T_D };

  if (q > 0) {
    mpv  = MMO(cx->dpf_avx[i], q-1);
    dpv  = DMO(cx->dpf_avx[i], q-1);
    tmdv = om_fs->tfv_avx[7*(q-1) + p7O_MD];
    tddv = om_fs->tfv_avx[7*Q + (q-1)];
  } else {
    mpv  = esl_avx_rightshiftz_float(MMO(cx->dpf_avx[i], Q-1));
    dpv  = esl_avx_rightshiftz_float(DMO(cx->dpf_avx[i], Q-1));
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
select_i_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int i, int k)
{
  int     Q    = p7O_NQF_AVX(cx->M);
  int     q    = (k-1) % Q;
  int     r    = (k-1) / Q;
  __m256  mpv  = MMO(cx->dpf_avx[i-3], q);
  __m256  ipv  = IMO(cx->dpf_avx[i-3], q);
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

/* C(i) is reached from C(i-3) or E(i). Values of rows i-2..i are put on
 * the scale of row i-3 with their scale factors. */
static inline int
select_c_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int i)
{
  float path[2];

  if (i < 3) return p7T_E;
  path[0] = cx->xmx[(i-3)*p7X_NXCELLS + p7X_C] * om_fs->xf[p7O_C][p7O_LOOP];
  path[1] = cx->xmx[ i   *p7X_NXCELLS + p7X_E] * om_fs->xf[p7O_E][p7O_MOVE]
          * cx->xmx[(i-2)*p7X_NXCELLS + p7X_SCALE] * cx->xmx[(i-1)*p7X_NXCELLS + p7X_SCALE] * cx->xmx[i*p7X_NXCELLS + p7X_SCALE];
  esl_vec_FNorm(path, 2);
  return (esl_rnd_FChoose(rng, path, 2) == 0) ? p7T_C : p7T_E;
}

/* The row the alignment ends at, as an offset 0, 1 or 2 back from row <L>:
 * the Forward score sums C(L) + C(L-1) tCC + C(L-2) tCC. */
static inline int
select_c_end_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int L)
{
  float loop = om_fs->xf[p7O_C][p7O_LOOP];
  float path[3];

  path[2] = (L >= 2) ? cx->xmx[(L-2)*p7X_NXCELLS + p7X_C] * loop : 0.;
  path[1] = (L >= 1) ? cx->xmx[(L-1)*p7X_NXCELLS + p7X_C] * loop * cx->xmx[(L-1)*p7X_NXCELLS + p7X_SCALE] : 0.;
  path[0] = cx->xmx[L*p7X_NXCELLS + p7X_C] * ((L >= 1) ? cx->xmx[(L-1)*p7X_NXCELLS + p7X_SCALE] : 1.) * cx->xmx[L*p7X_NXCELLS + p7X_SCALE];
  esl_vec_FNorm(path, 3);
  return esl_rnd_FChoose(rng, path, 3);
}



/* J(i) is reached from J(i-3) or E(i), with the scale handling of select_c_fs(). */
static inline int
select_j_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int i)
{
  float path[2];

  if (i < 3) return p7T_E;
  path[0] = cx->xmx[(i-3)*p7X_NXCELLS + p7X_J] * om_fs->xf[p7O_J][p7O_LOOP];
  path[1] = cx->xmx[ i   *p7X_NXCELLS + p7X_E] * om_fs->xf[p7O_E][p7O_LOOP]
          * cx->xmx[(i-2)*p7X_NXCELLS + p7X_SCALE] * cx->xmx[(i-1)*p7X_NXCELLS + p7X_SCALE] * cx->xmx[i*p7X_NXCELLS + p7X_SCALE];
  esl_vec_FNorm(path, 2);
  return (esl_rnd_FChoose(rng, path, 2) == 0) ? p7T_J : p7T_E;
}


/* E(i) is reached from any M(i,k=1..M) or D(i,k=2..M).
 * Uses on-the-fly FChoose algorithm with E(i) as normalizer.
 */
static inline int
select_e_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int i, int *ret_k)
{
  int    Q    = p7O_NQF_AVX(cx->M);
  double sum  = 0.0;
  double roll = esl_random(rng);
  double norm = 1.0 / cx->xmx[i*p7X_NXCELLS + p7X_E];
  __m256 xEv  = _mm256_set1_ps((float) norm);
  union { __m256 v; float p[8]; } u;
  int    q, r;

  while (1) {
    for (q = 0; q < Q; q++)
      {
        u.v = _mm256_mul_ps(MMO(cx->dpf_avx[i], q), xEv);
        for (r = 0; r < 8; r++) {
          sum += u.p[r];
          if (roll < sum) { *ret_k = r*Q + q + 1; return p7T_M; }
        }

        u.v = _mm256_mul_ps(DMO(cx->dpf_avx[i], q), xEv);
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
select_b_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int i)
{
  float path[2];
  int   state[2] = { p7T_N, p7T_J };

  path[0] = cx->xmx[i*p7X_NXCELLS + p7X_N] * om_fs->xf[p7O_N][p7O_MOVE];
  path[1] = cx->xmx[i*p7X_NXCELLS + p7X_J] * om_fs->xf[p7O_J][p7O_MOVE];
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}


#endif /* eslENABLE_AVX */
