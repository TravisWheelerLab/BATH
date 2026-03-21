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
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"

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
  *ret_sc = xC + om->xf[p7O_C][p7O_MOVE];
  return eslOK;
}
/*------------------ end, p7_Viterbi() -------------------------*/



/*****************************************************************
 * 2. Benchmark driver.
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
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "compare scores to generic implementation (debug)", 0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "20000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  { "--notrace", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "Do not benchmark the trace",                       0 }
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
  p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);
  p7_oprofile_Logify(om);

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

      if (esl_opt_GetBoolean(go, "-c"))
	  {
	    p7_GViterbi(dsq, L, gm, gx, &sc2);
	    printf("%.4f %.4f\n", sc, sc2);
	  }

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
  float        sc1, sc2;

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  p7_oprofile_Logify(om);
  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_Viterbi (dsq, L, om, ox, &sc1);
      p7_GViterbi(dsq, L, gm, gx, &sc2);
      if (fabs(sc1 - sc2) > 0.001) esl_fatal(msg);

      p7_Viterbi_Trace(dsq, L, om, ox, tr);
	  p7_GViterbi_Trace(dsq, L, gm, gx, trg);

	  if (p7_trace_fs_Validate(tr, abc, dsq, NULL)   != eslOK) esl_fatal(msg);
	  if (p7_trace_fs_Validate(trg, abc, dsq, NULL)  != eslOK) esl_fatal(msg);
	  if (p7_trace_Compare(tr, trg, 0.)              != eslOK) esl_fatal(msg);
	  
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


