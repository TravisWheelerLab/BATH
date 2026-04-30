/* Viterbi full matrix implementation; SSE version.
 *
 * This is a SIMD vectorized, striped, interleaved, full-matrix O(M*L)
 * memory implementation of the Viterbi algorithm, for calculating an
 * accurate Viterbi score with traceback.
 *
 * This implementation retains all DP rows in ox->dpf[0..L] and all
 * special-state scores in ox->xmx, allowing subsequent traceback via
 * p7_OTrace().
 *
 * The optimized profile must be configured to contain lspace float
 * scores (i.e., p7_oprofile_Logify() must have been called).
 *
 * Contents:
 *   1. Viterbi implementation.
 *   2. p7_Viterbi_Trace() implementation.
 *   3. Benchmark driver.
 *   4. Unit tests.
 *   5. Test driver.
 *
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_sse.h"

/*****************************************************************
 * 1. Viterbi implementation
 *****************************************************************/

/* Function:  p7_Viterbi()
 * Synopsis:  Full-matrix SSE Viterbi, enabling traceback.
 *
 * Purpose:   Calculates the Viterbi score for sequence <dsq> of length <L>
 *            residues, using optimized profile <om>, and a preallocated
 *            full-matrix DP matrix <ox>. Return the Viterbi score (in nats)
 *            in <ret_sc>.
 *
 *            Unlike <p7_ViterbiScore()>, the complete DP matrix is retained
 *            in <ox->dpf[0..L]> and special-state scores are stored in
 *            <ox->xmx>, enabling subsequent traceback by <p7_OTrace()>.
 *
 *            The model <om> must be configured with lspace float scores,
 *            not pspace float scores (call <p7_oprofile_Logify()> first).
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues
 *            om      - optimized profile (lspace float scores)
 *            ox      - DP matrix, allocated for >= L+1 rows and >= om->M cols
 *            ret_sc  - RETURN: Viterbi score (in nats)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small.
 */
int
p7_Viterbi(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register __m128 mpv, dpv, ipv;   /* previous row values                                       */
  register __m128 sv;		   /* temp storage of 1 curr row value in progress              */
  register __m128 dcv;		   /* delayed storage of D(i,q+1)                               */
  register __m128 xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register __m128 xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  __m128  infv;			   /* -eslINFINITY in a vector                                  */
  float    xN, xE, xB, xC, xJ;	   /* special states' scores                                    */
  int i;			   /* counter over sequence positions 1..L                      */
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q       = p7O_NQF(om->M);   /* segment length: # of vectors                              */
  __m128 *dpp;			   /* pointer to previous DP row ox->dpf[i-1]                  */
  __m128 *dp;			   /* pointer to current  DP row ox->dpf[i]                    */
  __m128 *rsc;			   /* will point at om->rfv[x] for residue x[i]                */
  __m128 *tsc;			   /* will point into (and step thru) om->tfv                  */
  float  *xmx = ox->xmx;	   /* special state scores [0..L][p7X_NXCELLS]                 */

  /* Check that the DP matrix is large enough. */
  if (Q > ox->allocQ4)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (Q)");
  if (L >= ox->validR)  ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small (rows)");
  if (ox->xmx == NULL)  ESL_EXCEPTION(eslEINVAL, "DP matrix has no xmx special state storage");
  ox->M = om->M;
  ox->L = L;

  /* Initialization of row 0. */
  infv = _mm_set1_ps(-eslINFINITY);
  dp   = ox->dpf[0];
  for (q = 0; q < Q; q++)
    MMO(dp,q) = IMO(dp,q) = DMO(dp,q) = infv;

  xN = 0.0f;
  xB = om->xf[p7O_N][p7O_MOVE];
  xE = -eslINFINITY;
  xJ = -eslINFINITY;
  xC = -eslINFINITY;

  XMXo(0, p7X_N) = xN;
  XMXo(0, p7X_B) = xB;
  XMXo(0, p7X_E) = xE;
  XMXo(0, p7X_J) = xJ;
  XMXo(0, p7X_C) = xC;

  for (i = 1; i <= L; i++)
    {
      rsc  = om->rfv[dsq[i]];
      tsc  = om->tfv;
      dcv  = infv;
      xEv  = infv;
      xBv  = _mm_set1_ps(xB);

      dpp  = ox->dpf[i-1];	/* previous row (read-only) */
      dp   = ox->dpf[i];	/* current row  (write)     */

      /* Shift the last vector of the previous row into the initial mpv/dpv/ipv.
       * Right-shift by one float position (4 bytes), padding with -inf at left. */
      mpv = esl_sse_rightshift_ps(MMO(dpp, Q-1), infv);
      dpv = esl_sse_rightshift_ps(DMO(dpp, Q-1), infv);
      ipv = esl_sse_rightshift_ps(IMO(dpp, Q-1), infv);

      for (q = 0; q < Q; q++)
	{
	  /* Calculate M(i,q); hold in sv (not yet stored). */
	  sv   =                _mm_add_ps(xBv, *tsc);  tsc++;  /* B->Mk  */
	  sv   = _mm_max_ps(sv, _mm_add_ps(mpv, *tsc)); tsc++;  /* Mk-1->Mk */
	  sv   = _mm_max_ps(sv, _mm_add_ps(ipv, *tsc)); tsc++;  /* Ik-1->Mk */
	  sv   = _mm_max_ps(sv, _mm_add_ps(dpv, *tsc)); tsc++;  /* Dk-1->Mk */
	  sv   = _mm_add_ps(sv, *rsc);                  rsc++;  /* + emission */
	  xEv  = _mm_max_ps(xEv, sv);		                 /* Mk->E    */

	  /* Load {MDI}(i-1,q) before we overwrite anything at q. */
	  mpv = MMO(dpp, q);
	  dpv = DMO(dpp, q);
	  ipv = IMO(dpp, q);

	  /* Store M(i,q) and D(i,q) [D value came from dcv of prior q]. */
	  MMO(dp, q) = sv;
	  DMO(dp, q) = dcv;

	  /* Partially compute D(i,q+1): M->D transition only; D->D added below. */
	  dcv = _mm_add_ps(sv, *tsc); tsc++;  /* Mk->Dk+1 */

	  /* Calculate and store I(i,q). Insert emission = log(1) = 0 (background freq). */
	  sv        =                _mm_add_ps(mpv, *tsc);  tsc++;  /* Mk->Ik   */
	  sv        = _mm_max_ps(sv, _mm_add_ps(ipv, *tsc)); tsc++;  /* Ik->Ik   */
	  IMO(dp,q) = sv;
	}

      /* Special states (E must be done; B must follow N,J). */
      esl_sse_hmax_ps(xEv, &xE);
      xN = xN +  om->xf[p7O_N][p7O_LOOP];
      xC = ESL_MAX(xC + om->xf[p7O_C][p7O_LOOP],  xE + om->xf[p7O_E][p7O_MOVE]);
      xJ = ESL_MAX(xJ + om->xf[p7O_J][p7O_LOOP],  xE + om->xf[p7O_E][p7O_LOOP]);
      xB = ESL_MAX(xJ + om->xf[p7O_J][p7O_MOVE],  xN + om->xf[p7O_N][p7O_MOVE]);

      XMXo(i, p7X_E) = xE;
      XMXo(i, p7X_N) = xN;
      XMXo(i, p7X_J) = xJ;
      XMXo(i, p7X_B) = xB;
      XMXo(i, p7X_C) = xC;

      /* D->D propagation.
       *
       * After the main q-loop, dcv holds M(i,Q-1)->D partial result.
       * Rightshift wraps it to position 0 (with -inf padding), then sweep
       * all Q stripes to propagate D->D transitions.  Repeat until no D
       * value improves (handles paths that cross a segment boundary).
       */
      dcv = esl_sse_rightshift_ps(dcv, infv);
      tsc = om->tfv + 7*Q;	/* DD transition block starts after 7*Q main transitions */
      for (q = 0; q < Q; q++)
	{
	  DMO(dp,q) = _mm_max_ps(dcv, DMO(dp,q));
	  dcv       = _mm_add_ps(DMO(dp,q), *tsc); tsc++;
	}

      /* Additional passes: D->D paths may cross one segment boundary. */
      do {
	dcv = esl_sse_rightshift_ps(dcv, infv);
	tsc = om->tfv + 7*Q;
	for (q = 0; q < Q; q++)
	  {
	    if (! esl_sse_any_gt_ps(dcv, DMO(dp,q))) break;
	    DMO(dp,q) = _mm_max_ps(dcv, DMO(dp,q));
	    dcv       = _mm_add_ps(DMO(dp,q), *tsc); tsc++;
	  }
      } while (q == Q);

    } /* end loop over sequence residues 1..L */

  /* C->T */
  if (ret_sc != NULL) *ret_sc = xC + om->xf[p7O_C][p7O_MOVE];
  return eslOK;
}
/*------------------ end, p7_Viterbi() -------------------------*/



/*****************************************************************
 * 2. p7_Viterbi_Trace() implementation
 *****************************************************************/


/* Function:  p7_Viterbi_Trace()
 * Synopsis:  Traceback from an lspace Viterbi DP matrix; SSE version.
 *
 * Purpose:   Recover the optimal alignment path from DP matrix <ox>
 *            filled by <p7_Viterbi()>.  The trace is written into <tr>
 *            (caller provides an initial allocation).
 *
 *            The profile <om> must be in lspace (p7_oprofile_Logify()
 *            called), matching the state in which <p7_Viterbi()> was run.
 *
 * Args:      dsq   - digital target sequence, 1..L
 *            L     - length of dsq
 *            om    - optimized profile (lspace float scores)
 *            ox    - Viterbi DP matrix from p7_Viterbi()
 *            tr    - RETURN: optimal traceback; caller provides alloc
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if an impossible state is reached.
 */
int
p7_Viterbi_Trace(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *ox,
                 P7_TRACE *tr)
{
  int   Q     = p7O_NQF(om->M);
  int   M     = om->M;
  int   i     = L;
  int   k     = 0;
  float r_tol = 1e-5;
  float a_tol = 1e-4;
  int   sprv, scur;
  int   status;

#define OMMo(i,k)  ((k)<1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_M, (i), (k)))
#define ODMo(i,k)  ((k)<1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_D, (i), (k)))
#define OIMo(i,k)  ((k)<1 ? -eslINFINITY : p7_omx_FGetMDI(ox, p7X_I, (i), (k)))
#define OXMXo(i,s)  (ox->xmx[(i)*p7X_NXCELLS+(s)])

#if eslDEBUGLEVEL > 0
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif

  if ((status = p7_trace_Append(tr, p7T_T, k, i)) != eslOK) return status;
  if ((status = p7_trace_Append(tr, p7T_C, k, i)) != eslOK) return status;
  sprv = p7T_C;

  while (sprv != p7T_S) {
    switch (sprv) {

    case p7T_C:
      if (OXMXo(i, p7X_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
      if      (esl_FCompare(OXMXo(i,p7X_C), OXMXo(i-1,p7X_C) + om->xf[p7O_C][p7O_LOOP], r_tol, a_tol) == eslOK) scur = p7T_C;
      else if (esl_FCompare(OXMXo(i,p7X_C), OXMXo(i,  p7X_E) + om->xf[p7O_E][p7O_MOVE], r_tol, a_tol) == eslOK) scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
      break;

    case p7T_E:
      if (OXMXo(i, p7X_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
      if (p7_oprofile_IsLocal(om)) {
        scur = p7T_M;
        for (k = M; k >= 1; k--) if (esl_FCompare(OXMXo(i,p7X_E), OMMo(i,k), r_tol, a_tol) == eslOK) break;
        if (k == 0) ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      } else {
        if      (esl_FCompare(OXMXo(i,p7X_E), OMMo(i,M), r_tol, a_tol) == eslOK) { scur = p7T_M; k = M; }
        else if (esl_FCompare(OXMXo(i,p7X_E), ODMo(i,M), r_tol, a_tol) == eslOK) { scur = p7T_D; k = M; }
        else ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      }
      break;

    case p7T_M: {
      union { __m128 v; float p[4]; } u;
      int   q   = (k-1) % Q,  r   = (k-1) / Q;
      float rsc = (u.v = om->rfv[dsq[i]][q],        u.p[r]);
      float tBM = (u.v = om->tfv[7*q + p7O_BM],     u.p[r]);
      float tMM = (u.v = om->tfv[7*q + p7O_MM],     u.p[r]);
      float tIM = (u.v = om->tfv[7*q + p7O_IM],     u.p[r]);
      float tDM = (u.v = om->tfv[7*q + p7O_DM],     u.p[r]);
      if (OMMo(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k, i);
      if      (esl_FCompare(OMMo(i,k), OXMXo(i-1,p7X_B) + tBM + rsc, r_tol, a_tol) == eslOK) scur = p7T_B;
      else if (esl_FCompare(OMMo(i,k), OMMo(i-1,k-1)   + tMM + rsc, r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(OMMo(i,k), OIMo(i-1,k-1)   + tIM + rsc, r_tol, a_tol) == eslOK) scur = p7T_I;
      else if (esl_FCompare(OMMo(i,k), ODMo(i-1,k-1)   + tDM + rsc, r_tol, a_tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "M at k=%d,i=%d couldn't be traced", k, i);
      k--; i--;
      break; }

    case p7T_D: {
      union { __m128 v; float p[4]; } u;
      int   qp  = (k-2) % Q, rp  = (k-2) / Q;   /* stripe coords of k-1 */
      float tMD = (k > 1) ? (u.v = om->tfv[7*qp + p7O_MD], u.p[rp]) : -eslINFINITY;
      float tDD = (k > 1) ? (u.v = om->tfv[7*Q  + qp],     u.p[rp]) : -eslINFINITY;
      if (ODMo(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k, i);
      if      (esl_FCompare(ODMo(i,k), OMMo(i,k-1) + tMD, r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(ODMo(i,k), ODMo(i,k-1) + tDD, r_tol, a_tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k, i);
      k--;
      break; }

    case p7T_I: {
      union { __m128 v; float p[4]; } u;
      int   q   = (k-1) % Q, r   = (k-1) / Q;
      float tMI = (u.v = om->tfv[7*q + p7O_MI], u.p[r]);
      float tII = (u.v = om->tfv[7*q + p7O_II], u.p[r]);
      if (OIMo(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k, i);
      if      (esl_FCompare(OIMo(i,k), OMMo(i-1,k) + tMI, r_tol, a_tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare(OIMo(i,k), OIMo(i-1,k) + tII, r_tol, a_tol) == eslOK) scur = p7T_I;
      else ESL_EXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k, i);
      i--;
      break; }

    case p7T_N:
      if (OXMXo(i, p7X_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = (i == 0) ? p7T_S : p7T_N;
      break;

    case p7T_B:
      if (OXMXo(i, p7X_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);
      if      (esl_FCompare(OXMXo(i,p7X_B), OXMXo(i,p7X_N) + om->xf[p7O_N][p7O_MOVE], r_tol, a_tol) == eslOK) scur = p7T_N;
      else if (esl_FCompare(OXMXo(i,p7X_B), OXMXo(i,p7X_J) + om->xf[p7O_J][p7O_MOVE], r_tol, a_tol) == eslOK) scur = p7T_J;
      else ESL_EXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    case p7T_J:
      if (OXMXo(i, p7X_J) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible J reached at i=%d", i);
      if      (esl_FCompare(OXMXo(i,p7X_J), OXMXo(i-1,p7X_J) + om->xf[p7O_J][p7O_LOOP], r_tol, a_tol) == eslOK) scur = p7T_J;
      else if (esl_FCompare(OXMXo(i,p7X_J), OXMXo(i,  p7X_E) + om->xf[p7O_E][p7O_LOOP], r_tol, a_tol) == eslOK) scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "J at i=%d couldn't be traced", i);
      break;

    default: ESL_EXCEPTION(eslEINVAL, "bogus state in Viterbi traceback");
    } /* end switch */

    if ((status = p7_trace_Append(tr, scur, k, i)) != eslOK) return status;
    if ((scur == p7T_N || scur == p7T_J || scur == p7T_C) && scur == sprv) i--;
    sprv = scur;
  } /* end traceback, at S state */

  tr->M = M;
  tr->L = L;

#undef OMMo
#undef ODMo
#undef OIMo
#undef OXMXo

  return p7_trace_Reverse(tr);
}
/*------------------ end, p7_Viterbi_Trace() --------------------*/


/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/
#ifdef p7VITERBI_BENCHMARK
/*
   gcc -g -O3 -msse2 -std=gnu99 -o viterbi_benchmark -I.. -L.. -I../../easel -L../../easel -Dp7VITERBI_BENCHMARK viterbi.c -lhmmer -leasel -lm
   icc  -O3 -static -o viterbi_benchmark -I.. -L.. -I../../easel -L../../easel -Dp7VITERBI_BENCHMARK viterbi.c -lhmmer -leasel -lm

   ./viterbi_benchmark <hmmfile>           runs benchmark
   ./viterbi_benchmark -c -N100 <hmmfile>  compare scores to generic Viterbi
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "20000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  { "--notrace", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "Do not benchmark the trace",                       0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for SSE Viterbi";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_OMX         *ox      = NULL;
  P7_GMX         *gx      = NULL;
  P7_TRACE       *tr      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc, sc2;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_UNILOCAL);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert_Log(gm, om);
  p7_oprofile_ReconfigLength_Log(om, L);

  ox = p7_omx_Create(gm->M, L, L);
  gx = p7_gmx_Create(gm->M, L, L, p7G_NSCELLS);
  tr = p7_trace_fs_Create();

  /* Baseline time: just sequence generation */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark time */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_Viterbi(dsq, L, om, ox, &sc);

	  if (! esl_opt_GetBoolean(go, "--notrace"))
	  {
        p7_Viterbi_Trace(dsq, L, om, ox, tr);  
	    p7_trace_Reuse(tr);
	  }
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n", gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_trace_fs_Destroy(tr);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7VITERBI_BENCHMARK*/
/*------------------- end, benchmark driver ---------------------*/



/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef p7VITERBI_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_viterbi()
 *
 * Compare p7_Viterbi() scores to p7_GViterbi() scores.  Differences
 * should be negligible floating-point roundoff only (< 0.001 nats).
 * The optimized profile must be in lspace (p7_oprofile_Logify()).
 */
static void
utest_viterbi(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char        *msg = "SSE p7_Viterbi() unit test failed";
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *ox  = p7_omx_Create(M, L, L);
  P7_GMX      *gx  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_TRACE    *tr  = p7_trace_fs_Create();
  P7_TRACE    *trg = p7_trace_fs_Create();
  char         errbuf[eslERRBUFSIZE];
  float        sc1, sc2, tsc1, tsc2;

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  p7_oprofile_Logify(om);
  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_Viterbi (dsq, L, om, ox, &sc1);
      p7_GViterbi(dsq, L, gm, gx, &sc2);
      if (fabs(sc1 - sc2) > 0.001) esl_fatal("%s: Viterbi scores disagree: %.4f vs %.4f", msg, sc1, sc2);

      p7_Viterbi_Trace(dsq, L, om, ox, tr);
      p7_GTrace(dsq, L, gm, gx, trg);

      if (p7_trace_Validate(tr,  abc, dsq, errbuf) != eslOK) esl_fatal("%s: SSE trace invalid: %s",     msg, errbuf);
      if (p7_trace_Validate(trg, abc, dsq, errbuf) != eslOK) esl_fatal("%s: generic trace invalid: %s", msg, errbuf);

      p7_trace_Score(tr,  dsq, gm, &tsc1);
      p7_trace_Score(trg, dsq, gm, &tsc2);
      if (fabs(tsc1 - tsc2) > 0.001) esl_fatal("%s: trace scores disagree: %.4f vs %.4f", msg, tsc1, tsc2);

      p7_trace_Reuse(tr);
      p7_trace_Reuse(trg);
    }

  free(dsq);
  p7_trace_fs_Destroy(tr);
  p7_trace_fs_Destroy(trg);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7VITERBI_TESTDRIVE*/
/*-------------------- end, unit tests --------------------------*/

/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7VITERBI_TESTDRIVE
/*
   gcc -g -Wall -msse2 -std=gnu99 -I.. -L.. -I../../easel -L../../easel -o viterbi_utest -Dp7VITERBI_TESTDRIVE viterbi.c -lhmmer -leasel -lm
   ./viterbi_utest
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE p7_Viterbi() full-matrix implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r   = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc = NULL;
  P7_BG          *bg  = NULL;
  int             M   = esl_opt_GetInteger(go, "-M");
  int             L   = esl_opt_GetInteger(go, "-L");
  int             N   = esl_opt_GetInteger(go, "-N");

  /* First round of tests for DNA alphabets. */
  if ((abc = esl_alphabet_Create(eslDNA))   == NULL) esl_fatal("failed to create alphabet");
  if ((bg  = p7_bg_Create(abc))             == NULL) esl_fatal("failed to create null model");

  utest_viterbi(r, abc, bg, M, L, N);   /* normal sized models  */
  utest_viterbi(r, abc, bg, 1, L, 10);  /* size 1 models        */
  utest_viterbi(r, abc, bg, M, 1, 10);  /* size 1 sequences     */

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets. */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL) esl_fatal("failed to create alphabet");
  if ((bg  = p7_bg_Create(abc))             == NULL) esl_fatal("failed to create null model");

  utest_viterbi(r, abc, bg, M, L, N);
  utest_viterbi(r, abc, bg, 1, L, 10);
  utest_viterbi(r, abc, bg, M, 1, 10);

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*p7VITERBI_TESTDRIVE*/
/*--------------------- end, test driver ------------------------*/


