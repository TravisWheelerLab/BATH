/* Viterbi filter implementation; NEON version.
 *
 * This is a SIMD vectorized, striped, interleaved, one-row, reduced
 * precision (epi16) implementation of the Viterbi algorithm.
 *
 * It calculates a close approximation of the Viterbi score, in
 * limited precision (signed words: 16 bits) and range. It may overflow on
 * high scoring sequences, but this indicates that the sequence is a
 * high-scoring hit worth examining more closely anyway.  It will not
 * underflow, in local alignment mode.
 *
 * Contents:
 *   1. Viterbi filter implementation.
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 *   5. Example.
 */
#include <p7_config.h>

#include <stdio.h>
#include <math.h>

#include <arm_neon.h>		/* NEON  */

#include "easel.h"
#include "esl_neon.h"
#include "esl_gumbel.h"

#include "hmmer.h"
#include "impl_neon.h"


/*****************************************************************
 * 1. Viterbi filter implementation.
 *****************************************************************/

/* Function:  p7_ViterbiFilter()
 * Synopsis:  Calculates Viterbi score, vewy vewy fast, in limited precision.
 * Incept:    SRE, Tue Nov 27 09:15:24 2007 [Janelia]
 *
 * Purpose:   Calculates an approximation of the Viterbi score for sequence
 *            <dsq> of length <L> residues, using optimized profile <om>,
 *            and a preallocated one-row DP matrix <ox>. Return the
 *            estimated Viterbi score (in nats) in <ret_sc>.
 *
 *            Score may overflow (and will, on high-scoring
 *            sequences), but will not underflow.
 *
 *            The model must be in a local alignment mode; other modes
 *            cannot provide the necessary guarantee of no underflow.
 *
 *            This is a striped SIMD Viterbi implementation using AMD NEON
 *            integer intrinsics \citep{Farrar07}, in reduced
 *            precision (signed words, 16 bits).
 *
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues
 *            om      - optimized profile
 *            ox      - DP matrix
 *            ret_sc  - RETURN: Viterbi score (in nats)
 *
 * Returns:   <eslOK> on success;
 *            <eslERANGE> if the score overflows; in this case
 *            <*ret_sc> is <eslINFINITY>, and the sequence can
 *            be treated as a high-scoring hit.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small, or if
 *            profile isn't in a local alignment mode. (Must be in local
 *            alignment mode because that's what helps us guarantee
 *            limited dynamic range.)
 *
 * Xref:      [Farrar07] for ideas behind striped SIMD DP.
 *            J2/46-47 for layout of HMMER's striped SIMD DP.
 *            J2/50 for single row DP.
 *            J2/60 for reduced precision (epu8)
 *            J2/65 for initial benchmarking
 *            J2/66 for precision maximization
 *            J4/138-140 for reimplementation in 16-bit precision
 */
int
p7_ViterbiFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc)
{
  register int16x8_t mpv, dpv, ipv;  /* previous row values                                       */
  register int16x8_t sv;		   /* temp storage of 1 curr row value in progress              */
  register int16x8_t dcv;		   /* delayed storage of D(i,q+1)                               */
  register int16x8_t xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register int16x8_t xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register int16x8_t Dmaxv;          /* keeps track of maximum D cell on row                      */
  int16_t  xE, xB, xC, xJ, xN;	   /* special states' scores                                    */
  int16_t  Dmax;		   /* maximum D cell score on row                               */
  int i;			   /* counter over sequence positions 1..L                      */
  int q;			   /* counter over vectors 0..nq-1                              */
  int Q        = p7O_NQW(om->M);   /* segment length: # of vectors                              */
  int16x8_t *dp  = ox->dpw[0];	   /* using {MDI}MX(q) macro requires initialization of <dp>    */
  int16x8_t *rsc;			   /* will point at om->ru[x] for residue x[i]                  */
  int16x8_t *tsc;			   /* will point into (and step thru) om->tu                    */

  int16x8_t negInfv;

  /* Check that the DP matrix is ok for us. */
  if (Q > ox->allocQ8)                                 ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  if (om->mode != p7_LOCAL && om->mode != p7_UNILOCAL) ESL_EXCEPTION(eslEINVAL, "Fast filter only works for local alignment");
  ox->M   = om->M;

  /* -infinity is -32768 */
  negInfv = vmovq_n_s16(-32768);  /* negInfv = 16-byte vector, with -32768 in all lanes, for a VEXT operation. */

  /* Initialization. In unsigned arithmetic, -infinity is -32768
   */
  for (q = 0; q < Q; q++)
    MMXo(q) = IMXo(q) = DMXo(q) = vmovq_n_s16(-32768);
  xN   = om->base_w;
  xB   = xN + om->xw[p7O_N][p7O_MOVE];
  xJ   = -32768;
  xC   = -32768;
  xE   = -32768;

#if eslDEBUGLEVEL > 0
  if (ox->debugging) p7_omx_DumpVFRow(ox, 0, xE, 0, xJ, xB, xC); /* first 0 is <rowi>: do header. second 0 is xN: always 0 here. */
#endif

  for (i = 1; i <= L; i++)
    {
      rsc   = om->rwv[dsq[i]];
      tsc   = om->twv;
      dcv   = vmovq_n_s16(-32768);      /* "-infinity" */
      xEv   = vmovq_n_s16(-32768);
      Dmaxv = vmovq_n_s16(-32768);
      xBv   = vmovq_n_s16(xB);

      /* Right shifts by 1 value (2 bytes). 4,8,12,x becomes x,4,8,12.
       * Because ia32 is littlendian, this means a left bit shift.
       * Zeros shift on automatically; replace it with -32768.
       */
      mpv = MMXo(Q-1);  mpv = vextq_s16(negInfv, mpv, 7);
      dpv = DMXo(Q-1);  dpv = vextq_s16(negInfv, dpv, 7);
      ipv = IMXo(Q-1);  ipv = vextq_s16(negInfv, ipv, 7);

      for (q = 0; q < Q; q++)
      {
        /* Calculate new MMXo(i,q); don't store it yet, hold it in sv. */
        sv   =                vqaddq_s16(xBv, *tsc);  tsc++;
        sv   = vmaxq_s16 (sv, vqaddq_s16(mpv, *tsc)); tsc++;
        sv   = vmaxq_s16 (sv, vqaddq_s16(ipv, *tsc)); tsc++;
        sv   = vmaxq_s16 (sv, vqaddq_s16(dpv, *tsc)); tsc++;
        sv   = vqaddq_s16(sv, *rsc);                  rsc++;
        xEv  = vmaxq_s16 (xEv, sv);

        /* Load {MDI}(i-1,q) into mpv, dpv, ipv;
         * {MDI}MX(q) is then the current, not the prev row
         */
        mpv = MMXo(q);
        dpv = DMXo(q);
        ipv = IMXo(q);

        /* Do the delayed stores of {MD}(i,q) now that memory is usable */
        MMXo(q) = sv;
        DMXo(q) = dcv;

        /* Calculate the next D(i,q+1) partially: M->D only;
               * delay storage, holding it in dcv
         */
        dcv   = vqaddq_s16(sv, *tsc);  tsc++;
        Dmaxv = vmaxq_s16(dcv, Dmaxv);

        /* Calculate and store I(i,q) */
        sv     =                vqaddq_s16(mpv, *tsc);  tsc++;
        IMXo(q)= vmaxq_s16 (sv, vqaddq_s16(ipv, *tsc)); tsc++;
      }

      /* Now the "special" states, which start from Mk->E (->C, ->J->B) */
      xE = esl_neon_hmax_s16((esl_neon_128i_t) xEv);
      if (xE >= 32767) { *ret_sc = eslINFINITY; return eslERANGE; }	/* immediately detect overflow */
      xN = xN + om->xw[p7O_N][p7O_LOOP];
      xC = ESL_MAX(xC + om->xw[p7O_C][p7O_LOOP], xE + om->xw[p7O_E][p7O_MOVE]);
      xJ = ESL_MAX(xJ + om->xw[p7O_J][p7O_LOOP], xE + om->xw[p7O_E][p7O_LOOP]);
      xB = ESL_MAX(xJ + om->xw[p7O_J][p7O_MOVE], xN + om->xw[p7O_N][p7O_MOVE]);
      /* and now xB will carry over into next i, and xC carries over after i=L */

      /* Finally the "lazy F" loop (sensu [Farrar07]). We can often
       * prove that we don't need to evaluate any D->D paths at all.
       *
       * The observation is that if we can show that on the next row,
       * B->M(i+1,k) paths always dominate M->D->...->D->M(i+1,k) paths
       * for all k, then we don't need any D->D calculations.
       *
       * The test condition is:
       *      max_k D(i,k) + max_k ( TDD(k-2) + TDM(k-1) - TBM(k) ) < xB(i)
       * So:
       *   max_k (TDD(k-2) + TDM(k-1) - TBM(k)) is precalc'ed in om->dd_bound;
       *   max_k D(i,k) is why we tracked Dmaxv;
       *   xB(i) was just calculated above.
       */
      Dmax = esl_neon_hmax_s16((esl_neon_128i_t) Dmaxv);
      if (Dmax + om->ddbound_w > xB)
	{
	  /* Now we're obligated to do at least one complete DD path to be sure. */
	  /* dcv has carried through from end of q loop above */
    dcv = vextq_s16(negInfv, dcv, 7);
	  tsc = om->twv + 7*Q;	/* set tsc to start of the DD's */
	  for (q = 0; q < Q; q++)
	    {
	      DMXo(q) = vmaxq_s16(dcv, DMXo(q));
	      dcv     = vqaddq_s16(DMXo(q), *tsc); tsc++;
	    }

	  /* We may have to do up to three more passes; the check
	   * is for whether crossing a segment boundary can improve
	   * our score.
	   */
	  do {
      dcv = vextq_s16(negInfv, dcv, 7);
	    tsc = om->twv + 7*Q;	/* set tsc to start of the DD's */
	    for (q = 0; q < Q; q++)
	      {
		if (! esl_neon_any_gt_s16((esl_neon_128i_t) dcv, (esl_neon_128i_t) DMXo(q))) break;
		DMXo(q) = vmaxq_s16(dcv, DMXo(q));
		dcv     = vqaddq_s16(DMXo(q), *tsc);   tsc++;
	      }
	  } while (q == Q);
	}
      else  /* not calculating DD? then just store the last M->D vector calc'ed.*/
	{
	  DMXo(0) = vextq_s16(negInfv, dcv, 7);
	}

#if eslDEBUGLEVEL > 0
      if (ox->debugging) p7_omx_DumpVFRow(ox, i, xE, 0, xJ, xB, xC);
#endif
    } /* end loop over sequence residues 1..L */

  /* finally C->T */
  if (xC > -32768)
    {
      *ret_sc = (float) xC + (float) om->xw[p7O_C][p7O_MOVE] - (float) om->base_w;
      /* *ret_sc += L * om->ncj_roundoff;  see J4/150 for rationale: superceded by -3.0nat approximation*/
      *ret_sc /= om->scale_w;
      *ret_sc -= 3.0; /* the NN/CC/JJ=0,-3nat approximation: see J5/36. That's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ contrib */
    }
  else  *ret_sc = -eslINFINITY;
  return eslOK;
}
/*---------------- end, p7_ViterbiFilter() ----------------------*/



/* Function:  p7_ViterbiFilter_BATH()
 * Synopsis:  Viterbi filter returning both a score and above-threshold diagonal windows.
 *
 * Purpose:   Runs the standard Viterbi DP (identical to p7_ViterbiFilter) without
 *            resetting the dp matrix, so the global Viterbi score accumulates normally.
 *            Additionally, at each row i where the row-maximum xE >= <sc_thresh>
 *            (derived from <filtersc> and <P> using Viterbi Gumbel parameters), the
 *            model position k_start responsible for that maximum is identified.  A
 *            forward diagonal extension in SSV score space then traces the diagonal
 *            from (i, k_start) until the match score drops, giving k_end and i_end.
 *            One window (n=i, k=k_end, k_length=k_end-k_start+1) is appended to
 *            <windowlist> and threshold checks are suppressed until i > i_end to
 *            avoid re-triggering on the same alignment.
 *
 *            The global Viterbi score is returned in <ret_sc> exactly as
 *            p7_ViterbiFilter() would return it, so the caller can apply a
 *            bias-corrected F2 filter using the windows for local composition.
 *
 * Args:      dsq        - digital target sequence, 1..L
 *            L          - length of dsq in residues
 *            om         - optimized profile
 *            ox         - DP matrix (one-row)
 *            ssvdata    - precomputed SSV score data (ssv_scores, for diagonal extension)
 *            filtersc   - bias-corrected null score (nats); used for sc_thresh
 *            P          - p-value threshold (F2); used for sc_thresh
 *            windowlist - RETURN: above-threshold diagonal windows appended here
 *            ret_sc     - RETURN: Viterbi score (nats)
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if score overflows; <*ret_sc> is <eslINFINITY>.
 *
 * Throws:    <eslEINVAL> if <ox> is too small or profile is not in local mode.
 */
int
p7_ViterbiFilter_BATH(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox,
                      const P7_SCOREDATA *ssvdata, float filtersc, double P,
                      P7_HMM_WINDOWLIST *windowlist, float *ret_sc)
{
  register int16x8_t mpv, dpv, ipv;
  register int16x8_t sv;
  register int16x8_t dcv;
  register int16x8_t xEv;
  register int16x8_t xBv;
  register int16x8_t Dmaxv;
  int16x8_t negInfv;
  int16_t   xE, xB, xC, xJ, xN;
  int16_t   Dmax;
  int       i;
  int       q;
  int       Q        = p7O_NQW(om->M);
  int16x8_t *dp      = ox->dpw[0];
  int16x8_t *rsc;
  int16x8_t *tsc;

  int16_t   sc_thresh;       /* Viterbi score threshold, int16 space              */
  int       sc_ext_thresh;   /* SSV score threshold for diagonal extension        */
  float     invP;
  int       z, k;
  int       skip_until = 0;  /* suppress xE check for i <= skip_until             */
  union { int16x8_t v; int16_t i[8]; } tmp;

  /* Viterbi threshold: invert P using Viterbi Gumbel params + filtersc */
  invP      = esl_gumbel_invsurv(P, om->evparam[p7_VMU], om->evparam[p7_VLAMBDA]);
  sc_thresh = (int16_t) ceil( ( (filtersc + eslCONST_LOG2 * invP + 3.0) * om->scale_w )
              - (float)om->xw[p7O_E][p7O_MOVE] - (float)om->xw[p7O_C][p7O_MOVE] + (float)om->base_w );

  /* SSV threshold for forward diagonal extension: invert P using MSV Gumbel params */
  invP          = esl_gumbel_invsurv(P, om->evparam[p7_MMU], om->evparam[p7_MLAMBDA]);
  sc_ext_thresh = (int) ceil( ( (filtersc + eslCONST_LOG2 * invP + 3.0) * om->scale_b )
                  + om->base_b + om->tec_b + om->tjb_b );

  if (Q > ox->allocQ8)                                 ESL_EXCEPTION(eslEINVAL, "DP matrix allocated too small");
  if (om->mode != p7_LOCAL && om->mode != p7_UNILOCAL) ESL_EXCEPTION(eslEINVAL, "Fast filter only works for local alignment");
  ox->M = om->M;

  negInfv = vmovq_n_s16(-32768);

  for (q = 0; q < Q; q++)
    MMXo(q) = IMXo(q) = DMXo(q) = vmovq_n_s16(-32768);
  xN = om->base_w;
  xB = xN + om->xw[p7O_N][p7O_MOVE];
  xJ = -32768;
  xC = -32768;
  xE = -32768;

  for (i = 1; i <= L; i++)
  {
    rsc   = om->rwv[dsq[i]];
    tsc   = om->twv;
    dcv   = vmovq_n_s16(-32768);
    xEv   = vmovq_n_s16(-32768);
    Dmaxv = vmovq_n_s16(-32768);
    xBv   = vmovq_n_s16(xB);

    mpv = MMXo(Q-1);  mpv = vextq_s16(negInfv, mpv, 7);
    dpv = DMXo(Q-1);  dpv = vextq_s16(negInfv, dpv, 7);
    ipv = IMXo(Q-1);  ipv = vextq_s16(negInfv, ipv, 7);

    for (q = 0; q < Q; q++)
    {
      sv   =                vqaddq_s16(xBv, *tsc);  tsc++;
      sv   = vmaxq_s16(sv,  vqaddq_s16(mpv, *tsc)); tsc++;
      sv   = vmaxq_s16(sv,  vqaddq_s16(ipv, *tsc)); tsc++;
      sv   = vmaxq_s16(sv,  vqaddq_s16(dpv, *tsc)); tsc++;
      sv   = vqaddq_s16(sv, *rsc);                  rsc++;
      xEv  = vmaxq_s16(xEv, sv);

      mpv     = MMXo(q);
      dpv     = DMXo(q);
      ipv     = IMXo(q);

      MMXo(q) = sv;
      DMXo(q) = dcv;

      dcv   = vqaddq_s16(sv, *tsc);  tsc++;
      Dmaxv = vmaxq_s16(dcv, Dmaxv);

      sv      =                vqaddq_s16(mpv, *tsc); tsc++;
      IMXo(q) = vmaxq_s16(sv,  vqaddq_s16(ipv, *tsc)); tsc++;
    }

    xE = esl_neon_hmax_s16((esl_neon_128i_t) xEv);
    if (xE >= 32767) { *ret_sc = eslINFINITY; return eslERANGE; }

    /* Special state updates — always performed, never skipped (no dp reset) */
    xN = xN +  om->xw[p7O_N][p7O_LOOP];
    xC = ESL_MAX(xC + om->xw[p7O_C][p7O_LOOP], xE + om->xw[p7O_E][p7O_MOVE]);
    xJ = ESL_MAX(xJ + om->xw[p7O_J][p7O_LOOP], xE + om->xw[p7O_E][p7O_LOOP]);
    xB = ESL_MAX(xJ + om->xw[p7O_J][p7O_MOVE], xN + om->xw[p7O_N][p7O_MOVE]);

    if (i > skip_until && xE >= sc_thresh)
    {
      /* Find k_start: first k in scan order where MMXo score == xE */
      int k_start = 0;
      for (q = 0; q < Q && k_start == 0; q++) {
        tmp.v = MMXo(q);
        for (z = 0; z < 8; z++) {
          k = q + Q*z + 1;
          if (tmp.i[z] == xE && k <= om->M) { k_start = k; break; }
        }
      }

      /* Forward diagonal extension in SSV score space from (i, k_start).
       * Starts at sc_ext_thresh (the SSV threshold) and extends M->M until
       * score stops improving for 5 consecutive steps. */
      int max_k_end     = k_start;
      int max_i_end     = i;
      int sc_ext        = sc_ext_thresh;
      int max_sc_ext    = sc_ext;
      int pos_since_max = 0;
      int kk = k_start + 1;
      int nn = i + 1;
      while (kk <= om->M && nn <= L) {
        sc_ext += om->bias_b - ssvdata->ssv_scores[kk * om->abc->Kp + dsq[nn]];
        if (sc_ext >= max_sc_ext) {
          max_sc_ext    = sc_ext;
          max_k_end     = kk;
          max_i_end     = nn;
          pos_since_max = 0;
        } else {
          if (++pos_since_max == 5) break;
        }
        kk++;
        nn++;
      }

      p7_hmmwindow_new(windowlist, 0, i, max_k_end, max_k_end - k_start + 1, 0.0, p7_NOCOMPLEMENT, L);
      skip_until = max_i_end;
    }

    /* Lazy F loop */
    Dmax = esl_neon_hmax_s16((esl_neon_128i_t) Dmaxv);
    if (Dmax + om->ddbound_w > xB)
    {
      dcv = vextq_s16(negInfv, dcv, 7);
      tsc = om->twv + 7*Q;
      for (q = 0; q < Q; q++)
      {
        DMXo(q) = vmaxq_s16(dcv, DMXo(q));
        dcv     = vqaddq_s16(DMXo(q), *tsc); tsc++;
      }
      do {
        dcv = vextq_s16(negInfv, dcv, 7);
        tsc = om->twv + 7*Q;
        for (q = 0; q < Q; q++)
        {
          if (!esl_neon_any_gt_s16((esl_neon_128i_t) dcv, (esl_neon_128i_t) DMXo(q))) break;
          DMXo(q) = vmaxq_s16(dcv, DMXo(q));
          dcv     = vqaddq_s16(DMXo(q), *tsc); tsc++;
        }
      } while (q == Q);
    }
    else
    {
      DMXo(0) = vextq_s16(negInfv, dcv, 7);
    }
  }

  /* Global Viterbi score via C->T */
  if (xC > -32768) {
    *ret_sc  = (float) xC + (float) om->xw[p7O_C][p7O_MOVE] - (float) om->base_w;
    *ret_sc /= om->scale_w;
    *ret_sc -= 3.0;
  } else *ret_sc = -eslINFINITY;

  return eslOK;
}
/*---------------- end, p7_ViterbiFilter_BATH() -----------------*/



/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
#ifdef p7VITFILTER_BENCHMARK
/* -c, -x are used for debugging, testing; see msvfilter.c for explanation */

/*
   gcc -g -O3 -march=armv8-a -std=gnu99 -o vitfilter_benchmark -I.. -L.. -I../../easel -L../../easel -Dp7VITFILTER_BENCHMARK vitfilter.c -lhmmer -leasel -lm
   ./benchmark-vitfilter <hmmfile>          runs benchmark
   ./benchmark-vitfilter -N100 -c <hmmfile> compare scores to generic impl
   ./benchmark-vitfilter -N100 -x <hmmfile> compare scores to exact emulation
 */
#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_neon.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-x", "compare scores to generic implementation (debug)", 0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-c", "equate scores to trusted implementation (debug)",  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for Viterbi filter";


int
main(int argc, char **argv)
{
  ESL_GETOPTS       *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char              *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH     *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS    *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET      *abc     = NULL;
  P7_HMMFILE        *hfp     = NULL;
  P7_HMM            *hmm     = NULL;
  P7_BG             *bg      = NULL;
  P7_PROFILE        *gm      = NULL;
  P7_OPROFILE       *om      = NULL;
  P7_SCOREDATA      *data    = NULL;
  P7_OMX            *ox      = NULL;
  P7_GMX            *gx      = NULL;
  P7_HMM_WINDOWLIST  wlist;
  int                L       = esl_opt_GetInteger(go, "-L");
  int                N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ           *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int                i;
  float              sc1, sc2;
  float              nullsc;
  double             base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)           != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);
  data = p7_hmm_ScoreDataCreate(om, NULL);

  if (esl_opt_GetBoolean(go, "-x")) p7_profile_SameAsVF(om, gm);

  ox = p7_omx_Create(gm->M, 0, 0);
  gx = p7_gmx_Create(gm->M, L, L, p7G_NSCELLS);
  p7_hmmwindow_init(&wlist);

  /* Get a baseline time: how long it takes just to generate the sequences */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark p7_ViterbiFilter() */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_ViterbiFilter(dsq, L, om, ox, &sc1);

      if (esl_opt_GetBoolean(go, "-c"))
	{
	  p7_GViterbi(dsq, L, gm, gx, &sc2);
	  printf("%.4f %.4f\n", sc1, sc2);
	}

      if (esl_opt_GetBoolean(go, "-x"))
	{
	  p7_GViterbi(dsq, L, gm, gx, &sc2);
	  sc2 /= om->scale_w;
	  if (om->mode == p7_UNILOCAL)   sc2 -= 2.0; /* that's ~ L \log \frac{L}{L+2}, for our NN,CC,JJ */
	  else if (om->mode == p7_LOCAL) sc2 -= 3.0; /* that's ~ L \log \frac{L}{L+3}, for our NN,CC,JJ */
	  printf("%.4f %.4f\n", sc1, sc2);
	}
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# ViterbiFilter       CPU time: ");
  printf("# M    = %d\n", gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  /* Benchmark p7_ViterbiFilter_BATH() */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
      p7_bg_NullOne(bg, dsq, L, &nullsc);
      wlist.count = 0;
      p7_ViterbiFilter_BATH(dsq, L, om, ox, data, nullsc, 0.2, &wlist, &sc1);
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# ViterbiFilter_BATH  CPU time: ");
  printf("# M    = %d\n", gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(wlist.windows);
  free(dsq);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_hmm_ScoreDataDestroy(data);
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
#endif /*p7VITFILTER_BENCHMARK*/
/*---------------- end, benchmark driver ------------------------*/




/*****************************************************************
 * 3. Unit tests.
 *****************************************************************/
#ifdef p7VITFILTER_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* ViterbiFilter() unit test
 *
 * We can check that scores are identical (within machine error) to
 * scores of generic DP with scores rounded the same way.  Do this for
 * a random model of length <M>, for <N> test sequences of length <L>.
 *
 * We assume that we don't accidentally generate a high-scoring random
 * sequence that overflows ViterbiFilter()'s limited range.
 *
 */
static void
utest_viterbi_filter(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_DSQ     *dsq = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *ox  = p7_omx_Create(M, 0, 0);
  P7_GMX      *gx  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  float sc1, sc2;

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  p7_profile_SameAsVF(om, gm);	/* round and scale the scores in <gm> the same as in <om> */

#if 0
  p7_oprofile_Dump(stdout, om);              // dumps the optimized profile
  p7_omx_SetDumpMode(stdout, ox, TRUE);      // makes the fast DP algorithms dump their matrices
#endif

  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      p7_ViterbiFilter(dsq, L, om, ox, &sc1);
      p7_GViterbi     (dsq, L, gm, gx, &sc2);

#if 0
      p7_gmx_Dump(stdout, gx, p7_DEFAULT);   // dumps a generic DP matrix
      break;
#endif

      sc2 /= om->scale_w;
      sc2 -= 3.0;

      if (fabs(sc1-sc2) > 0.001) esl_fatal("viterbi filter unit test failed: scores differ (%.2f, %.2f)", sc1, sc2);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}

/* utest_viterbi_filter_bath()
 *
 * Checks two properties of p7_ViterbiFilter_BATH():
 *
 * 1. Score agreement: the returned Viterbi score must match
 *    p7_ViterbiFilter() to within floating-point rounding (0.001 nats),
 *    because both run the same DP without resetting the matrix.
 *
 * 2. Window validity: every window appended to the windowlist must have
 *    n in [1,L], k in [1,M], and k_length in [1,M].
 */
static void
utest_viterbi_filter_bath(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char              msg[]  = "ViterbiFilter_BATH unit test failed";
  P7_HMM           *hmm   = NULL;
  P7_PROFILE       *gm    = NULL;
  P7_OPROFILE      *om    = NULL;
  P7_SCOREDATA     *data  = NULL;
  ESL_DSQ          *dsq   = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX           *ox    = p7_omx_Create(M, 0, 0);
  P7_HMM_WINDOWLIST wlist;
  float             sc1, sc2;
  float             nullsc;
  int               w;

  p7_hmmwindow_init(&wlist);
  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  data = p7_hmm_ScoreDataCreate(om, NULL);

  while (N--)
    {
      esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);

      /* Reference score from standard ViterbiFilter */
      p7_ViterbiFilter(dsq, L, om, ox, &sc1);

      /* Score + windows from ViterbiFilter_BATH; use nullsc as filtersc */
      p7_bg_NullOne(bg, dsq, L, &nullsc);
      wlist.count = 0;
      p7_ViterbiFilter_BATH(dsq, L, om, ox, data, nullsc, 0.001, &wlist, &sc2);

      /* 1. Scores must agree */
      if (sc1 != eslINFINITY && sc2 != eslINFINITY)
        if (fabs(sc1 - sc2) > 0.001)
          esl_fatal("%s: scores differ: ViterbiFilter=%.4f  ViterbiFilter_BATH=%.4f", msg, sc1, sc2);

      /* 2. Every window must have valid boundaries */
      for (w = 0; w < wlist.count; w++)
        {
          if (wlist.windows[w].n        <  1) esl_fatal("%s: window %d n=%d < 1",             msg, w, wlist.windows[w].n);
          if (wlist.windows[w].n        >  L) esl_fatal("%s: window %d n=%d > L=%d",          msg, w, wlist.windows[w].n, L);
          if (wlist.windows[w].k        <  1) esl_fatal("%s: window %d k=%d < 1",             msg, w, wlist.windows[w].k);
          if (wlist.windows[w].k        >  M) esl_fatal("%s: window %d k=%d > M=%d",          msg, w, wlist.windows[w].k, M);
          if (wlist.windows[w].k_length <  1) esl_fatal("%s: window %d k_length=%d < 1",      msg, w, wlist.windows[w].k_length);
          if (wlist.windows[w].k_length >  M) esl_fatal("%s: window %d k_length=%d > M=%d",   msg, w, wlist.windows[w].k_length, M);
        }
    }

  free(wlist.windows);
  free(dsq);
  p7_hmm_ScoreDataDestroy(data);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(ox);
  p7_profile_Destroy(gm);
  p7_oprofile_Destroy(om);
}
#endif /*p7VITFILTER_TESTDRIVE*/


/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7VITFILTER_TESTDRIVE
/*
   gcc -g -Wall -march=armv8-a -std=gnu99 -o vitfilter_utest -I.. -L.. -I../../easel -L../../easel -Dp7VITFILTER_TESTDRIVE vitfilter.c -lhmmer -leasel -lm
*/
#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "hmmer.h"
#include "impl_neon.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "-L",        eslARG_INT,    "200", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the NEON implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = NULL;
  P7_BG          *bg   = NULL;
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");

  /* First round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))            == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiFilter() tests, DNA\n");
  utest_viterbi_filter(r, abc, bg, M, L, N);
  utest_viterbi_filter(r, abc, bg, 1, L, 10);
  utest_viterbi_filter(r, abc, bg, M, 1, 10);

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiFilter_BATH() tests, DNA\n");
  utest_viterbi_filter_bath(r, abc, bg, M, L, N);
  utest_viterbi_filter_bath(r, abc, bg, 1, L, 10);
  utest_viterbi_filter_bath(r, abc, bg, M, 1, 10);

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiFilter() tests, protein\n");
  utest_viterbi_filter(r, abc, bg, M, L, N);
  utest_viterbi_filter(r, abc, bg, 1, L, 10);
  utest_viterbi_filter(r, abc, bg, M, 1, 10);

  if (esl_opt_GetBoolean(go, "-v")) printf("ViterbiFilter_BATH() tests, protein\n");
  utest_viterbi_filter_bath(r, abc, bg, M, L, N);
  utest_viterbi_filter_bath(r, abc, bg, 1, L, 10);
  utest_viterbi_filter_bath(r, abc, bg, M, 1, 10);

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*VITFILTER_TESTDRIVE*/
/*---------------- end, test driver -----------------------------*/



/*****************************************************************
 * 5. Example
 *****************************************************************/
#ifdef p7VITFILTER_EXAMPLE
/* A minimal example.
   Also useful for debugging on small HMMs and sequences.
   gcc -g -Wall -march=armv8-a -std=gnu99 -o vitfilter_example -I.. -L.. -I../../easel -L../../easel -Dp7VITFILTER_EXAMPLE vitfilter.c -lhmmer -leasel -lm
   ./vitfilter_example <hmmfile> <seqfile>
 */
#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_gumbel.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "impl_neon.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-1",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "output in one line awkable format",                0 },
  { "-P",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "output in profmark format",                        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of Viterbi filter algorithm";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_OMX         *ox      = NULL;
  P7_GMX         *gx      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  float           vfraw, nullsc, vfscore;
  float           graw, gscore;
  double          P, gP;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)           != eslOK) p7_Fail("Failed to read HMM");

  /* Read in one sequence */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_Open(seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);

  /* create default null model, then create and optimize profile */
  bg = p7_bg_Create(abc);
  p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);
  p7_oprofile_Convert(gm, om);

  /* allocate DP matrices, both a generic and an optimized one */
  ox = p7_omx_Create(gm->M, 0, sq->n);
  gx = p7_gmx_Create(gm->M, sq->n, sq->n, p7G_NSCELLS);

  /* Useful to place and compile in for debugging:
     p7_oprofile_Dump(stdout, om);         dumps the optimized profile
     p7_omx_SetDumpMode(ox, TRUE);         makes the fast DP algorithms dump their matrices
     p7_gmx_Dump(stdout, gx, p7_DEFAULT);  dumps a generic DP matrix
  */

  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK)
    {
      p7_oprofile_ReconfigLength(om, sq->n);
      p7_ReconfigLength(gm,          sq->n);
      p7_bg_SetLength(bg,            sq->n);
      p7_omx_GrowTo(ox, om->M, 0,    sq->n);
      p7_gmx_GrowTo(gx, gm->M,       sq->n, sq->n);

      p7_ViterbiFilter  (sq->dsq, sq->n, om, ox, &vfraw);
      p7_bg_NullOne (bg, sq->dsq, sq->n, &nullsc);
      vfscore = (vfraw - nullsc) / eslCONST_LOG2;
      P        = esl_gumbel_surv(vfscore,  om->evparam[p7_VMU],  om->evparam[p7_VLAMBDA]);

      p7_GViterbi       (sq->dsq, sq->n, gm, gx, &graw);
      gscore   = (graw - nullsc) / eslCONST_LOG2;
      gP       = esl_gumbel_surv(gscore,  gm->evparam[p7_VMU],  gm->evparam[p7_VLAMBDA]);

      if (esl_opt_GetBoolean(go, "-1"))
	{
	  printf("%-30s\t%-20s\t%9.2g\t%7.2f\t%9.2g\t%7.2f\n", sq->name, hmm->name, P, vfscore, gP, gscore);
	}
      else if (esl_opt_GetBoolean(go, "-P"))
	{ /* output suitable for direct use in profmark benchmark postprocessors: */
	  printf("%g\t%.2f\t%s\t%s\n", P, vfscore, sq->name, hmm->name);
	}
      else
	{
	  printf("target sequence:      %s\n",        sq->name);
	  printf("vit filter raw score: %.2f nats\n", vfraw);
	  printf("null score:           %.2f nats\n", nullsc);
	  printf("per-seq score:        %.2f bits\n", vfscore);
	  printf("P-value:              %g\n",        P);
	  printf("GViterbi raw score:   %.2f nats\n", graw);
	  printf("GViterbi seq score:   %.2f bits\n", gscore);
	  printf("GViterbi P-value:     %g\n",        gP);
	}

      esl_sq_Reuse(sq);
    }

  /* cleanup */
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  p7_omx_Destroy(ox);
  p7_gmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7VITFILTER_EXAMPLE*/
/*-------------------- end, example -----------------------------*/
