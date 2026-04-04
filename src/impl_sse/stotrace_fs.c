/* SSE implementation of frameshift stochastic backtrace of a Forward matrix.
 * (Compare generic version, p7_GStochasticTrace_Frameshift().)
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

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_random.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_sse.h"

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
 * 1. Stochastic trace implementation.
 *****************************************************************/

/* Function:  p7_StochasticTrace_Frameshift()
 * Synopsis:  Sample a traceback from a frameshift-aware Forward matrix; SSE version.
 *
 * Purpose:   Identical to <p7_GStochasticTrace_Frameshift()> except that
 *            <om_fs>, <ox> are SSE optimized versions of the frameshift
 *            profile and Forward matrix. See <p7_GStochasticTrace_Frameshift()>
 *            documentation.
 *
 *            The Forward matrix <ox> must have been computed by
 *            <p7_Forward_Frameshift()> and uses the 8-cell FS layout
 *            (p7X_NSCELLS_FS=8: D, I, M_C0..M_C5 per stripe).
 *
 *            Scale factors stored in ox->xmx[i*p7X_NXCELLS + p7X_SCALE] are
 *            used to normalize values across rows when comparing states at
 *            different sequence positions (as in C/J codon-length choices).
 *
 * Args:      rng   - source of random numbers
 *            dsq   - digital sequence, 1..L
 *            L     - length of dsq
 *            om_fs - optimized frameshift profile
 *            ox    - Forward matrix (FS 8-cell layout)
 *            tr    - RETURN: sampled traceback; caller provides initial alloc
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> on various internal failures (impossible state reached,
 *            bogus state in traceback).
 */
int
p7_StochasticTrace_Frameshift(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L,
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

      /* For M states: stochastically select codon length c (1..5) from C1..C5 cells.
       * If the codon would extend before the sequence start, revert to B.
       */
      if (s1 == p7T_M) {
        c = select_codon_len_fs(rng, ox, i, k, Q);
        if (i - c < 0) s1 = p7T_B;  /* codon would start before position 1 */
      } else {
        c = 0;
      }

      if ((status = p7_trace_fs_Append(tr, s1, k, i, c)) != eslOK) return status;

      /* For NCJ states: deferred i decrement (one nucleotide per loop iteration). */
      if ((s1 == p7T_N || s1 == p7T_C || s1 == p7T_J) && s1 == s0) i--;

      s0 = s1;
      i -= c;
    }

  tr->M = om_fs->M;
  tr->L = L;
  return p7_trace_fs_Reverse(tr);
}
/*------------------ end, stochastic traceback ------------------*/


/*****************************************************************
 * 2. Selection of steps in the traceback
 *****************************************************************/
/* The guts of the stochastic backtrace is broken into select_?_fs() helpers.
 * Each randomly selects one of the possible predecessors and returns the
 * state type.
 *
 * Key differences from the standard SSE select_*() functions:
 *  - The FS dp matrix has 8 cells per stripe (p7X_NSCELLS_FS) vs 3.
 *    Match cells are accessed via MMO_FS(dp,q,c), delete via DMO_FS, insert via IMO_FS.
 *  - In the FS model, M(i,k) connects from B(i), M(i,k-1), D(i,k-1), I(i,k-1)
 *    all at the same row i (codon offsets are handled separately via codon length c).
 *  - I(i,k) connects from M(i-3,k) and I(i-3,k) (3-nucleotide codon).
 *  - C(i) and J(i) each connect from three codon-length offsets (i-3, i-2, i-1)
 *    as well as E(i).  Scale factors from ox->xmx[j*p7X_NXCELLS+p7X_SCALE] are
 *    used to bring values at different rows to a common scale before normalization.
 *  - All dp values are in pspace (linear), so esl_vec_FNorm is used (not FLogNorm).
 */

/* M(i,k) is reached from B(i), M(i,k-1), I(i,k-1), or D(i,k-1).
 * All predecessors are at the same row i (codon position is resolved later via c).
 */
static inline int
select_m_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;     /* stripe of current k */
  int     r     = (k-1) / Q;     /* element within stripe */
  __m128 *tp    = om_fs->tfv + 7*q;  /* transitions: BM, MM, IM, DM, MD, MI, II */
  __m128  xBv   = _mm_set1_ps(ox->xmx[i*p7X_NXCELLS + p7X_B]);
  __m128  mpv, dpv, ipv;
  union { __m128 v; float p[4]; } u;
  float   path[4];
  int     state[4] = { p7T_B, p7T_M, p7T_I, p7T_D };

  /* Predecessor k-1: stripe q-1, or wrap to Q-1 with rightshiftz. */
  if (q > 0) {
    mpv = MMO_FS(ox->dpf[i], q-1, p7X_FS_C0);
    dpv = DMO_FS(ox->dpf[i], q-1);
    ipv = IMO_FS(ox->dpf[i], q-1);
  } else {
    mpv = esl_sse_rightshiftz_float(MMO_FS(ox->dpf[i], Q-1, p7X_FS_C0));
    dpv = esl_sse_rightshiftz_float(DMO_FS(ox->dpf[i], Q-1));
    ipv = esl_sse_rightshiftz_float(IMO_FS(ox->dpf[i], Q-1));
  }

  u.v = _mm_mul_ps(xBv, *tp); tp++;  path[0] = u.p[r];  /* B  * T_BM */
  u.v = _mm_mul_ps(mpv, *tp); tp++;  path[1] = u.p[r];  /* M' * T_MM */
  u.v = _mm_mul_ps(ipv, *tp); tp++;  path[2] = u.p[r];  /* I' * T_IM */
  u.v = _mm_mul_ps(dpv, *tp);        path[3] = u.p[r];  /* D' * T_DM */
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
    tddv = om_fs->tfv[7*Q + (q-1)];       /* DD transitions follow 7*Q block */
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

  u.v = _mm_mul_ps(mpv, *tp); tp++;  path[0] = u.p[r];  /* M * T_MI */
  u.v = _mm_mul_ps(ipv, *tp);        path[1] = u.p[r];  /* I * T_II */
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
 *
 * Scale factor adjustment: stored values at different rows carry different cumulative
 * scale factors.  Using cumscale(i-3) as the common denominator:
 *   path[0] = C(i-3) * T_CL                                (no adjustment)
 *   path[1] = C(i-2) * T_CL * S(i-2)
 *   path[2] = C(i-1) * T_CL * S(i-2) * S(i-1)
 *   path[3] = E(i)   * T_EM * S(i-2) * S(i-1) * S(i)
 * where S(j) = ox->xmx[j*p7X_NXCELLS + p7X_SCALE].
 */
static inline int
select_c_fs(ESL_RANDOMNESS *rng, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i)
{
  float loop = om_fs->xf[p7O_C][p7O_LOOP];
  float s_im2, s_im1, s_i;
  float path[4];

  if (i < 4) return p7T_E;  /* too close to start; can only come from E */

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
 * Uses the on-the-fly FChoose algorithm (avoids 2M-1 allocation) with E(i) as normalizer.
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

/* select_codon_len_fs()
 * For M state at (i,k): stochastically select codon length c in 1..5 from
 * the C1..C5 cells in the FS dp matrix.  All cells are at the same row i,
 * so no scale-factor adjustment is needed.
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
 *   1. Generic and SSE forward scores agree to tolerance.
 *   2. Each SSE stochastic trace is structurally valid: N >= 3, starts with
 *      S, ends with T.
 *   3. Both the generic and SSE stochastic trace functions complete without
 *      error on every sequence.
 */
static void
utest_stotrace_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
                  ESL_GENCODE *gcode, P7_BG *bgAA,
                  P7_CODONTABLE *codon_table, int M, int N)
{
  char           *msg    = "utest_stotrace_fs failed";
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_fs5 = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE *om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq    = NULL;
  P7_TRACE       *tr     = p7_trace_fs_Create();
  P7_IVX         *iv5    = p7_ivx_Create(M, p7P_5CODONS);
  P7_OIVX        *ov5    = p7_oivx_Create(M, p7P_5CODONS);
  P7_GMX         *gxf    = p7_gmx_Create(M, M, M, p7G_NSCELLS_FS);
  P7_OMX         *fwd    = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS_FS);
  float           gfsc, ofsc;
  float           tolerance;
  int             i, j, curr_L, idx;

  tolerance = (p7_FLogsumError(-0.4, -0.5) > 0.0001) ? 1.0f : 0.001f;

  if (!gm || !gm_fs5 || !om_fs5 || !sq || !tr || !iv5 || !gxf || !fwd) esl_fatal(msg);

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
      if (p7_GForward_Frameshift(dsq, curr_L, gm_fs5, gxf, iv5, &gfsc) != eslOK) esl_fatal(msg);

      /* SSE forward */
      p7_omx_GrowTo_dpf(fwd, M, curr_L, curr_L);
      { int s = p7_Forward_Frameshift(dsq, curr_L, om_fs5, fwd, ov5, &ofsc);
        if (s == eslERANGE) continue;
        if (s != eslOK)     esl_fatal(msg); }

      /* Forward scores must agree */
      if (esl_FCompare_old(gfsc, ofsc, tolerance) != eslOK)
        esl_fatal("%s: generic fwd %.4f != SSE fwd %.4f (tol %.4f)", msg, gfsc, ofsc, tolerance);

      /* Sample SSE stochastic traces and validate each with p7_trace_fs_Validate */
      for (idx = 0; idx < 10; idx++)
        {
          char errbuf[eslERRBUFSIZE];
          if (p7_StochasticTrace_Frameshift(r, dsq, curr_L, om_fs5, fwd, tr) != eslOK)
            esl_fatal("%s: SSE stochastic trace failed", msg);
          if (p7_trace_fs_Validate(tr, abcDNA, dsq, errbuf) != eslOK)
            esl_fatal("%s: SSE trace invalid: %s", msg, errbuf);
          p7_trace_Reuse(tr);
        }

      /* Verify generic stochastic trace also runs cleanly */
      if (p7_GStochasticTrace_Frameshift(r, dsq, curr_L, gm_fs5, gxf, tr) != eslOK)
        esl_fatal("%s: generic stochastic trace failed", msg);
      p7_trace_Reuse(tr);
    }

  free(dsq);
  p7_trace_fs_Destroy(tr);
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
#endif /*p7STOTRACE_FS_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/



/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7STOTRACE_FS_TESTDRIVE
/*
  gcc -g -Wall -msse2 -std=gnu99 -o stotrace_fs_utest -I.. -L.. -I../../easel -L../../easel -Dp7STOTRACE_FS_TESTDRIVE stotrace_fs.c -lhmmer -leasel -lm
  ./stotrace_fs_utest
*/
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name    type           default env range toggles reqs incomp help                                       docgroup */
  { "-h",  eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",  eslARG_INT,     "42", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M",  eslARG_INT,     "72", NULL, NULL, NULL, NULL, NULL, "size of random models to sample",                0 },
  { "-N",  eslARG_INT,     "10", NULL, NULL, NULL, NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE frameshift stochastic traceback";

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


