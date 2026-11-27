/* NEON implementation of frameshift stochastic backtrace of a checkpointed Forward matrix.
 * (Compare generic version, p7_GStochasticTraceChk_Frameshift().)
 *
 * Contents:
 *    1. Stochastic trace implementation.
 *    2. Selection of steps in the traceback.
 *    3. Benchmark driver.
 *    4. Unit tests.
 *    5. Test driver.
 *
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <arm_neon.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_neon.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_neon.h"

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
 * 1. Stochastic trace implementation.
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

/* Function:  p7_StochasticTraceChk_Frameshift()
 * Synopsis:  Sample tracebacks from a checkpointed frameshift Forward matrix.
 *
 * Purpose:   Sample <ntr> stochastic tracebacks of the 5-codon profile
 *            <om_fs> against <dsq> of length <L> from the checkpointed
 *            Forward matrix <cx> that <p7_ForwardChk_Frameshift_5Codons()>
 *            filled: trace <t> into <tr[t]>, with its own random number
 *            generator <rng[t]>.
 *
 *            All traces are sampled together. The blocks of <cx> are
 *            loaded from the last to the first, each computed again from
 *            its checkpoint by <p7_ForwardChkLoad_Frameshift_5Codons()> with
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
 *            cx    - checkpointed Forward matrix from p7_ForwardChk_Frameshift_5Codons()
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
p7_StochasticTraceChk_Frameshift(ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,
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
      if (b != cx->b && (status = p7_ForwardChkLoad_Frameshift_5Codons(dsq, om_fs, cx, ov, b)) != eslOK) goto ERROR;
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
/*------------------ end, stochastic traceback ------------------*/


/*****************************************************************
 * 2. Selection of steps in the traceback
 *****************************************************************/
/* The guts of the stochastic backtrace is broken into select_?_fs() helpers.
 * Each randomly selects one of the possible predecessors and returns the
 * state type.
 *
 * Key differences from the standard NEON select_*() functions:
 *  - The rows read are the loaded rows of a checkpointed Forward matrix,
 *    3 cells per stripe (p7X_NSCELLS): M (the C0 sum), D and I.
 *  - M(i,k) comes from B, M, I or D at row i-c, model position k-1, for
 *    codon lengths c = 1..5; the codon length and the state are chosen
 *    together.
 *  - I(i,k) connects from M(i-3,k) and I(i-3,k) (3-nucleotide codon).
 *  - C(i) and J(i) connect from C(i-3) / J(i-3) and E(i); the first C step
 *    also chooses whether the alignment ends at row L, L-1 or L-2. Scale
 *    factors from cx->xmx[j*p7X_NXCELLS+p7X_SCALE] bring values at
 *    different rows to a common scale before normalization.
 *  - All dp values are in pspace (linear).
 */

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
  int     Q     = p7O_NQF(cx->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  int     st[4] = { p7T_B, p7T_M, p7T_I, p7T_D };
  float32x4_t  *tp    = om_fs->tfv + 7*q;   /* transitions into M_k: BM, MM, IM, DM */
  float32x4_t  zerov = vdupq_n_f32(0.0f);
  float32x4_t   mpv, dpv, ipv, xBv;
  union { float32x4_t v; float p[4]; } u;
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
      xBv = vdupq_n_f32(cx->xmx[p*p7X_NXCELLS + p7X_B]);
      if (q > 0) {
        mpv = MMO(cx->dpf[p], q-1);
        dpv = DMO(cx->dpf[p], q-1);
        ipv = IMO(cx->dpf[p], q-1);
      } else {
        mpv = vextq_f32(zerov, MMO(cx->dpf[p], Q-1), 3);
        dpv = vextq_f32(zerov, DMO(cx->dpf[p], Q-1), 3);
        ipv = vextq_f32(zerov, IMO(cx->dpf[p], Q-1), 3);
      }
      u.v = vmulq_f32(xBv, tp[0]); tw[0] = u.p[r];
      u.v = vmulq_f32(mpv, tp[1]); tw[1] = u.p[r];
      u.v = vmulq_f32(ipv, tp[2]); tw[2] = u.p[r];
      u.v = vmulq_f32(dpv, tp[3]); tw[3] = u.p[r];
      u.v = om_fs->rfv[cidx[c]][q];
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
  int          Q     = p7O_NQF(cx->M);
  int          q     = (k-1) % Q;
  int          r     = (k-1) / Q;
  float32x4_t  zerov = vdupq_n_f32(0.0f);
  float32x4_t  mpv, dpv;
  float32x4_t  tmdv, tddv;
  union { float32x4_t v; float p[4]; } u;
  float   path[2];
  int     state[2] = { p7T_M, p7T_D };

  if (q > 0) {
    mpv  = MMO(cx->dpf[i], q-1);
    dpv  = DMO(cx->dpf[i], q-1);
    tmdv = om_fs->tfv[7*(q-1) + p7O_MD];
    tddv = om_fs->tfv[7*Q + (q-1)];       /* DD transitions follow 7*Q block */
  } else {
    mpv  = vextq_f32(zerov, MMO(cx->dpf[i], Q-1), 3);
    dpv  = vextq_f32(zerov, DMO(cx->dpf[i], Q-1), 3);
    tmdv = vextq_f32(zerov, om_fs->tfv[7*(Q-1) + p7O_MD], 3);
    tddv = vextq_f32(zerov, om_fs->tfv[8*Q-1], 3);
  }

  u.v = vmulq_f32(mpv, tmdv); path[0] = u.p[r];
  u.v = vmulq_f32(dpv, tddv); path[1] = u.p[r];
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}

/* I(i,k) is reached from M(i-3,k) or I(i-3,k). Inserts emit 3-nucleotide codons. */
static inline int
select_i_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int i, int k)
{
  int          Q    = p7O_NQF(cx->M);
  int          q    = (k-1) % Q;
  int          r    = (k-1) / Q;
  float32x4_t  mpv  = MMO(cx->dpf[i-3], q);
  float32x4_t  ipv  = IMO(cx->dpf[i-3], q);
  float32x4_t *tp   = om_fs->tfv + 7*q + p7O_MI;
  union { float32x4_t v; float p[4]; } u;
  float   path[2];
  int     state[2] = { p7T_M, p7T_I };

  u.v = vmulq_f32(mpv, *tp); tp++;  path[0] = u.p[r];  /* M * T_MI */
  u.v = vmulq_f32(ipv, *tp);        path[1] = u.p[r];  /* I * T_II */
  esl_vec_FNorm(path, 2);
  return state[esl_rnd_FChoose(rng, path, 2)];
}

/* N(i) must come from N(i-1) for i>0, or S for i==0. */
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
 * Uses the on-the-fly FChoose algorithm (avoids 2M-1 allocation) with E(i) as normalizer.
 */
static inline int
select_e_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_CHKOMX *cx, int i, int *ret_k)
{
  int          Q    = p7O_NQF(cx->M);
  double       sum  = 0.0;
  double       roll = esl_random(rng);
  double       norm = 1.0 / cx->xmx[i*p7X_NXCELLS + p7X_E];
  float32x4_t  xEv  = vdupq_n_f32((float) norm);
  union { float32x4_t v; float p[4]; } u;
  int    q, r;

  while (1) {
    for (q = 0; q < Q; q++)
      {
        u.v = vmulq_f32(MMO(cx->dpf[i], q), xEv);
        for (r = 0; r < 4; r++) {
          sum += u.p[r];
          if (roll < sum) { *ret_k = r*Q + q + 1; return p7T_M; }
        }

        u.v = vmulq_f32(DMO(cx->dpf[i], q), xEv);
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

/*---------------------- end, step selection --------------------*/

/*****************************************************************
 * 3. Benchmark
 *****************************************************************/
/*----------------- end, benchmark ------------------------------*/


/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7STOTRACE_FS_TESTDRIVE
#include "esl_getopts.h"
#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_stotrace_fs()
 * Tests:
 *   1. Generic and NEON Forward scores agree to within <tolerance> nats.
 *   2. Every trace sampled by p7_StochasticTraceChk_Frameshift() from the
 *      checkpointed Forward matrix (blocks of 5..12 rows) passes
 *      p7_trace_fs_Validate().
 */
static void
utest_stotrace_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
                  ESL_GENCODE *gcode, P7_BG *bgAA,
                  P7_CODONTABLE *codon_table, int M, int N)
{
  char            *msg    = "utest_stotrace_fs failed";
  int              ntr    = 10;
  P7_HMM          *hmm    = NULL;
  P7_PROFILE      *gm     = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE   *gm_fs5 = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE  *om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
  ESL_SQ          *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ         *dsq    = NULL;
  ESL_RANDOMNESS **rng    = malloc(sizeof(ESL_RANDOMNESS *) * ntr);
  P7_TRACE       **tr     = malloc(sizeof(P7_TRACE *) * ntr);
  P7_IVX          *iv5    = p7_ivx_Create(M, p7P_5CODONS);
  P7_OIVX         *ov5    = p7_oivx_Create(M, p7P_5CODONS);
  P7_GMX          *gxf    = p7_gmx_Create(M, M, M, p7G_NSCELLS_FS);
  P7_OMX          *fwd    = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS_FS);
  P7_CHKOMX       *cx     = p7_chkomx_Create(M, M, 8, 5, p7X_NSCELLS, p7_CHK_FORWARD);
  char             errbuf[eslERRBUFSIZE];
  float            gfsc, ofsc;
  float            tolerance;
  int              i, j, t, curr_L;

  tolerance = (p7_FLogsumError(-0.4, -0.5) > 0.0001) ? 0.1f : 0.001f;   /* nats; the generic Forward's logsum table rounds */

  if (!gm || !gm_fs5 || !om_fs5 || !sq || !rng || !tr || !iv5 || !gxf || !fwd || !cx) esl_fatal(msg);
  for (t = 0; t < ntr; t++) { rng[t] = esl_randomness_Create(1); tr[t] = p7_trace_fs_Create(); }

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig(hmm, bgAA, gm, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, M, p7_LOCAL);
  p7_fs_oprofile_Convert(gm_fs5, om_fs5);

  while (N--)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, NULL);
      curr_L = sq->n * 3;

      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (curr_L + 2))) == NULL) esl_fatal("malloc failed");
      j = 1;
      for (i = 1; i <= (int)sq->n; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }

      p7_fs_oprofile_ReconfigLength(om_fs5, sq->n);
      p7_fs_ReconfigLength(gm_fs5, sq->n);

      /* Generic forward */
      p7_gmx_GrowTo(gxf, M, curr_L, curr_L);
      p7_ivx_GrowTo(iv5, M, p7P_5CODONS);
      if (p7_GForward_Frameshift_5Codons(dsq, curr_L, gm_fs5, gxf, iv5, &gfsc) != eslOK) esl_fatal(msg);

      /* NEON forward */
      p7_omx_GrowTo_dpf(fwd, M, curr_L, curr_L);
      { int s = p7_Forward_Frameshift_5Codons(dsq, curr_L, om_fs5, fwd, ov5, &ofsc);
        if (s == eslERANGE) continue;
        if (s != eslOK)     esl_fatal(msg); }

      /* Forward scores must agree */
      if (fabsf(gfsc - ofsc) > tolerance)
        esl_fatal("%s: generic fwd %.4f != NEON fwd %.4f (tol %.4f)", msg, gfsc, ofsc, tolerance);

      /* traces from the checkpointed Forward matrix must be valid */
      if (p7_chkomx_GrowTo(cx, M, curr_L, 5 + esl_rnd_Roll(r, 8))      != eslOK) esl_fatal(msg);
      if (p7_ForwardChk_Frameshift_5Codons(dsq, curr_L, om_fs5, cx, ov5, NULL) != eslOK) esl_fatal(msg);
      for (t = 0; t < ntr; t++) esl_randomness_Init(rng[t], 1 + esl_rnd_Roll(r, 1000000));
      if (p7_StochasticTraceChk_Frameshift(rng, ntr, dsq, curr_L, om_fs5, cx, ov5, tr) != eslOK)
        esl_fatal("%s: NEON stochastic trace failed", msg);
      for (t = 0; t < ntr; t++)
        {
          if (p7_trace_fs_Validate(tr[t], abcDNA, dsq, errbuf) != eslOK)
            esl_fatal("%s: NEON trace invalid: %s", msg, errbuf);
          p7_trace_Reuse(tr[t]);
        }
    }

  for (t = 0; t < ntr; t++) { esl_randomness_Destroy(rng[t]); p7_trace_fs_Destroy(tr[t]); }
  free(rng); free(tr);
  free(dsq);
  p7_chkomx_Destroy(cx);
  p7_omx_Destroy(fwd);
  p7_ivx_Destroy(iv5);
  p7_gmx_Destroy(gxf);
  p7_oivx_Destroy(ov5);
  p7_hmm_Destroy(hmm);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_fs5);
  p7_fs_oprofile_Destroy(om_fs5);
}


/* utest_stotrace_posterior_fs()
 * NEON stochastic traces must follow the posterior distribution. As the
 * generic utest_stotrace_posterior_frameshift(): for the back-translated
 * consensus of a sharpened model between random flanks with two nucleotides
 * deleted and two inserted, and frameshift probability 0.05, the mean number
 * per trace of M states of each codon length and of I states, over <N>
 * traces sampled together by p7_StochasticTraceChk_Frameshift() (blocks of
 * 5..12 rows), must be within 5 standard errors (plus 0.01) of the
 * posterior sums from the generic p7_GDecoding_Frameshift(), and every
 * trace must be valid.
 */
static void
utest_stotrace_posterior_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA,
                            P7_CODONTABLE *codon_table, int M, int N)
{
  char            *msg    = "frameshift stochastic trace posterior unit test failed";
  P7_HMM          *hmm    = NULL;
  P7_BG           *bgDNA  = p7_bg_Create(abcDNA);
  P7_FS_PROFILE   *gm_fs5 = p7_profile_fs_Create(M, abcAA, p7P_5CODONS);
  P7_FS_OPROFILE  *om_fs5 = NULL;
  int              L      = 3*M + 120;
  P7_GMX          *fx     = p7_gmx_Create(M, L, L, p7G_NSCELLS_FS);
  P7_GMX          *bx     = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_IVX          *iv     = p7_ivx_Create(M, p7P_5CODONS);
  P7_CHKOMX       *cx     = NULL;
  P7_OIVX         *ov     = NULL;
  ESL_SQ          *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ         *dsq    = malloc(sizeof(ESL_DSQ) * (L+2));
  ESL_DSQ         *tmp    = malloc(sizeof(ESL_DSQ) * (L+2));
  ESL_RANDOMNESS **rng    = malloc(sizeof(ESL_RANDOMNESS *) * N);
  P7_TRACE       **tr     = malloc(sizeof(P7_TRACE *) * N);
  double           expct[6], sum[6], sumsq[6], mean, se;
  char             errbuf[eslERRBUFSIZE];
  int              cnt[6];
  int              i, j, k, x, n, z, c, Ln, start, del1, del2, ins1, ins2;

  for (n = 0; n < N; n++) { rng[n] = esl_randomness_Create(1 + esl_rnd_Roll(r, 1000000)); tr[n] = p7_trace_fs_Create(); }

  if (p7_hmm_Sample(r, M, abcAA, &hmm) != eslOK) esl_fatal(msg);
  for (k = 1; k <= M; k++)
    {
      x = esl_vec_FArgMax(hmm->mat[k], abcAA->K);
      esl_vec_FSet(hmm->mat[k], abcAA->K, 0.15 / (abcAA->K - 1));
      hmm->mat[k][x] = 0.85;
      if (k < M) { hmm->t[k][p7H_MM] = 0.9; hmm->t[k][p7H_MI] = 0.05; hmm->t[k][p7H_MD] = 0.05; }
    }
  hmm->fsprob = 0.05;
  p7_emit_SimpleConsensus(hmm, sq);

  Ln = 3*M + 100;
  esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, tmp);
  start = 1 + esl_rnd_Roll(r, Ln - 3*M + 1);
  for (i = 1, j = start; i <= M; i++, j += 3) p7_codontable_GetCodon(codon_table, r, sq->dsq[i], tmp+j);
  del1 = start + M/2;
  ins1 = start + M   + esl_rnd_Roll(r, M/4);
  del2 = start + 3*M/2;
  ins2 = start + 2*M + esl_rnd_Roll(r, M/4);
  for (i = 1, j = 1; i <= Ln; i++)
    {
      if (i == del1 || i == del2) continue;
      if (i == ins1 || i == ins2) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = tmp[i];
    }
  Ln = j - 1;

  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, Ln/3, p7_LOCAL) != eslOK) esl_fatal(msg);
  om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
  p7_fs_oprofile_Convert(gm_fs5, om_fs5);
  cx     = p7_chkomx_Create(M, L, 5 + esl_rnd_Roll(r, 8), 5, p7X_NSCELLS, p7_CHK_FORWARD);
  ov     = p7_oivx_Create(M, p7P_5CODONS);
  if (p7_chkomx_GrowTo(cx, M, Ln, cx->blockR)                  != eslOK) esl_fatal(msg);
  if (p7_ForwardChk_Frameshift_5Codons(dsq, Ln, om_fs5, cx, ov, NULL)    != eslOK) esl_fatal(msg);
  if (p7_StochasticTraceChk_Frameshift(rng, N, dsq, Ln, om_fs5, cx, ov, tr) != eslOK) esl_fatal(msg);

  for (n = 0; n < 6; n++) { sum[n] = sumsq[n] = 0.; }
  for (n = 0; n < N; n++)
    {
      if (p7_trace_fs_Validate(tr[n], abcDNA, dsq, errbuf) != eslOK) esl_fatal("%s: invalid trace: %s", msg, errbuf);
      for (c = 0; c < 6; c++) cnt[c] = 0;
      for (z = 0; z < tr[n]->N; z++)
        {
          if      (tr[n]->st[z] == p7T_M) { if (tr[n]->c[z] < 1 || tr[n]->c[z] > 5) esl_fatal("%s: M with codon length %d", msg, tr[n]->c[z]); cnt[tr[n]->c[z]]++; }
          else if (tr[n]->st[z] == p7T_I) cnt[0]++;
        }
      for (c = 0; c < 6; c++) { sum[c] += cnt[c]; sumsq[c] += (double) cnt[c] * cnt[c]; }
    }

  if (p7_GForward_Frameshift_5Codons (dsq, Ln, gm_fs5, fx, iv, NULL) != eslOK) esl_fatal(msg);
  if (p7_GBackward_Frameshift_5Codons(dsq, Ln, gm_fs5, bx, iv, NULL) != eslOK) esl_fatal(msg);
  if (p7_GDecoding_Frameshift(gm_fs5, fx, bx)               != eslOK) esl_fatal(msg);
  for (c = 0; c < 6; c++) expct[c] = 0.;
  for (i = 1; i <= Ln; i++)
    for (k = 1; k <= M; k++)
      {
        for (c = 1; c <= 5; c++) expct[c] += fx->dp[i][k*p7G_NSCELLS_FS + p7G_M + c];
        expct[0] += fx->dp[i][k*p7G_NSCELLS_FS + p7G_I];
      }
  for (c = 0; c < 6; c++)
    {
      mean = sum[c] / N;
      se   = sqrt(ESL_MAX(0., sumsq[c] / N - mean * mean) / N);
      if (fabs(mean - expct[c]) > 5. * se + 0.01)
        esl_fatal("%s: %s %d: mean count %.4f per trace, posterior sum %.4f (s.e. %.4f)", msg, c ? "M codon length" : "I", c, mean, expct[c], se);
    }

  for (n = 0; n < N; n++) { esl_randomness_Destroy(rng[n]); p7_trace_fs_Destroy(tr[n]); }
  free(rng); free(tr);
  free(dsq); free(tmp);
  esl_sq_Destroy(sq);
  p7_ivx_Destroy(iv);
  p7_gmx_Destroy(fx);
  p7_gmx_Destroy(bx);
  p7_profile_fs_Destroy(gm_fs5);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bgDNA);
  p7_chkomx_Destroy(cx);
  p7_oivx_Destroy(ov);
  p7_fs_oprofile_Destroy(om_fs5);
}

/* utest_stotrace_chk_fs()
 * Stochastic traces from a checkpointed frameshift Forward matrix in one
 * block against the same matrix in blocks of 5..12 rows. For random DNA
 * and the back-translated consensus of a model with sharpened match
 * emissions between random flanks, <ntr> traces sampled together by
 * p7_StochasticTraceChk_Frameshift() from each must be identical (states,
 * k, i and codon lengths) with the same per-trace random number
 * generators, and valid.
 */
static void
utest_stotrace_chk_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA,
                      P7_CODONTABLE *codon_table, int M, int L, int N)
{
  char            *msg    = "checkpointed frameshift stochastic trace unit test failed";
  int              ntr    = 20;
  P7_HMM          *hmm    = NULL;
  P7_BG           *bgDNA  = p7_bg_Create(abcDNA);
  P7_FS_PROFILE   *gm_fs5 = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE  *om_fs5 = NULL;
  P7_OIVX         *ov5    = NULL;
  P7_CHKOMX       *cx1    = NULL;
  P7_CHKOMX       *cx2    = NULL;
  ESL_SQ          *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ         *dsq    = malloc(sizeof(ESL_DSQ) * (L+2));
  ESL_RANDOMNESS **rng    = malloc(sizeof(ESL_RANDOMNESS *) * ntr);
  P7_TRACE       **tr1    = malloc(sizeof(P7_TRACE *) * ntr);
  P7_TRACE       **tr2    = malloc(sizeof(P7_TRACE *) * ntr);
  uint32_t        *seed   = malloc(sizeof(uint32_t) * ntr);
  char             errbuf[eslERRBUFSIZE];
  int              ntested = 0;
  int              n, Ln, blockR, start, i, j, t, z, k, x;

  for (t = 0; t < ntr; t++) { rng[t] = esl_randomness_Create(1); tr1[t] = p7_trace_fs_Create(); tr2[t] = p7_trace_fs_Create(); }

  if (p7_hmm_Sample(r, M, abcAA, &hmm) != eslOK) esl_fatal(msg);
  for (k = 1; k <= M; k++)
    {
      x = esl_vec_FArgMax(hmm->mat[k], abcAA->K);
      esl_vec_FSet(hmm->mat[k], abcAA->K, 0.15 / (abcAA->K - 1));
      hmm->mat[k][x] = 0.85;
      if (k < M) { hmm->t[k][p7H_MM] = 0.9; hmm->t[k][p7H_MI] = 0.05; hmm->t[k][p7H_MD] = 0.05; }
    }
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_LOCAL) != eslOK) esl_fatal(msg);
  p7_emit_SimpleConsensus(hmm, sq);
  om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
  p7_fs_oprofile_Convert(gm_fs5, om_fs5);
  ov5    = p7_oivx_Create(M, p7P_5CODONS);
  cx1    = p7_chkomx_Create(M, L, L+1, 5, p7X_NSCELLS, p7_CHK_FORWARD);
  cx2    = p7_chkomx_Create(M, L, 8,   5, p7X_NSCELLS, p7_CHK_FORWARD);

  for (n = 0; n < N; n++)
    {
      Ln = 15 + esl_rnd_Roll(r, L - 14);
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, dsq);
      if (n % 2 == 1 && 3*M <= L)
        {
          Ln    = ESL_MAX(Ln, 3*M);
          esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, dsq);
          start = 1 + esl_rnd_Roll(r, Ln - 3*M + 1);
          for (i = 1, j = start; i <= M; i++, j += 3)
            p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq+j);
        }
      blockR = 5 + esl_rnd_Roll(r, 8);

      p7_fs_oprofile_ReconfigLength(om_fs5, Ln/3);
      if (p7_chkomx_GrowTo(cx1, M, Ln, Ln+1)                    != eslOK) esl_fatal(msg);
      if (p7_chkomx_GrowTo(cx2, M, Ln, blockR)                  != eslOK) esl_fatal(msg);
      if (p7_ForwardChk_Frameshift_5Codons(dsq, Ln, om_fs5, cx1, ov5, NULL) != eslOK) continue;
      if (p7_ForwardChk_Frameshift_5Codons(dsq, Ln, om_fs5, cx2, ov5, NULL) != eslOK) esl_fatal(msg);

      for (t = 0; t < ntr; t++) { seed[t] = 1 + esl_rnd_Roll(r, 1000000); esl_randomness_Init(rng[t], seed[t]); }
      if (p7_StochasticTraceChk_Frameshift(rng, ntr, dsq, Ln, om_fs5, cx1, ov5, tr1) != eslOK) esl_fatal(msg);
      for (t = 0; t < ntr; t++) esl_randomness_Init(rng[t], seed[t]);
      if (p7_StochasticTraceChk_Frameshift(rng, ntr, dsq, Ln, om_fs5, cx2, ov5, tr2) != eslOK) esl_fatal(msg);
      for (t = 0; t < ntr; t++)
        {
          if (p7_trace_fs_Validate(tr1[t], abcDNA, dsq, errbuf) != eslOK) esl_fatal("%s: invalid trace: %s", msg, errbuf);
          if (tr1[t]->N != tr2[t]->N) esl_fatal("%s: trace %d lengths differ (%d %d)", msg, t, tr1[t]->N, tr2[t]->N);
          for (z = 0; z < tr1[t]->N; z++)
            if (tr1[t]->st[z] != tr2[t]->st[z] || tr1[t]->k[z] != tr2[t]->k[z] || tr1[t]->i[z] != tr2[t]->i[z] || tr1[t]->c[z] != tr2[t]->c[z])
              esl_fatal("%s: trace %d differs at %d (block size %d)", msg, t, z, blockR);
          p7_trace_Reuse(tr1[t]);
          p7_trace_Reuse(tr2[t]);
        }
      ntested++;
    }
  if (ntested == 0) esl_fatal("%s: no sequence was tested", msg);

  for (t = 0; t < ntr; t++) { esl_randomness_Destroy(rng[t]); p7_trace_fs_Destroy(tr1[t]); p7_trace_fs_Destroy(tr2[t]); }
  free(rng); free(tr1); free(tr2); free(seed); free(dsq);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bgDNA);
  p7_profile_fs_Destroy(gm_fs5);
  p7_oivx_Destroy(ov5);
  p7_chkomx_Destroy(cx1);
  p7_chkomx_Destroy(cx2);
  p7_fs_oprofile_Destroy(om_fs5);
}
#endif /*p7STOTRACE_FS_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/



/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7STOTRACE_FS_TESTDRIVE
/*
  gcc -g -Wall -march=armv8-a -std=gnu99 -o stotrace_fs_utest -I.. -L.. -I../../easel -L../../easel -Dp7STOTRACE_FS_TESTDRIVE stotrace_fs.c -lhmmer -leasel -lm
  ./stotrace_fs_utest
*/
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_neon.h"

static ESL_OPTIONS options[] = {
  /* name    type           default env range toggles reqs incomp help                                       docgroup */
  { "-h",  eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",  eslARG_INT,     "42", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M",  eslARG_INT,     "72", NULL, NULL, NULL, NULL, NULL, "size of random models to sample",                0 },
  { "-N",  eslARG_INT,     "10", NULL, NULL, NULL, NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for NEON frameshift stochastic traceback";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA  = esl_alphabet_Create(eslAMINO);
  ESL_ALPHABET   *abcDNA = esl_alphabet_Create(eslDNA);
  ESL_GENCODE    *gcode  = esl_gencode_Create(abcDNA, abcAA);
  P7_CODONTABLE  *ct     = p7_codontable_Create(gcode);
  P7_BG          *bgAA   = p7_bg_Create(abcAA);
  int             M      = esl_opt_GetInteger(go, "-M");
  int             N      = esl_opt_GetInteger(go, "-N");

  p7_FLogsumInit();

  utest_stotrace_fs(r, abcAA, abcDNA, gcode, bgAA, ct, M, N);
  utest_stotrace_chk_fs(r, abcAA, abcDNA, gcode, bgAA, ct, M, 3*M+150, N);
  utest_stotrace_posterior_fs(r, abcAA, abcDNA, gcode, bgAA, ct, 45, 2000);

  p7_bg_Destroy(bgAA);
  p7_codontable_Destroy(ct);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcAA);
  esl_alphabet_Destroy(abcDNA);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7STOTRACE_FS_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/


/*****************************************************************
 * 6. Example.
 *****************************************************************/

/*------------------------ end, example -------------------------*/
