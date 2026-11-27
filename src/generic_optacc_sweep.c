/* Optimal accuracy alignment in a forward sweep, generic (non-SIMD) version.
 *
 * The Backward matrix is first computed checkpointed (<p7_GBackwardChk()>).
 * One sweep through the rows 1..L then computes, for each row, the
 * Forward row, the Backward row (its block computed again from its
 * checkpoint), the posterior probabilities and the optimal accuracy
 * (OA) scores, keeping only the rows the next row needs. What the OA
 * trace reads is stored in a <P7_OATMX> as each row is finished, and
 * <p7_GOATraceSweep()> traces from it alone, so no full-size matrix is
 * kept. Traces are identical to those of <p7_GForward()>,
 * <p7_GBackward()>, <p7_GDecoding()>, <p7_GOptimalAccuracy()> and
 * <p7_GOATrace()>, with each residue's posterior probability given as
 * the value its PP character stands for. The frameshift version is in
 * generic_optacc_sweep_frameshift.c.
 *
 * Computing posterior probabilities from Backward rows as they are
 * finished, without a full Backward matrix, follows HMMER4's Backward
 * filter; a checkpointed Backward matrix that a forward sweep reads
 * row by row, computing Forward, posteriors and OA together, is an
 * elaboration of that idea.
 *
 * Contents:
 *   1. OA sweep and trace
 *   2. Unit tests
 *   3. Test driver
 */
#include "p7_config.h"

#include <float.h>
#include <string.h>

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#define TSCDELTA(s,k) ( (tsc[(k) * p7P_NTRANS + (s)] == -eslINFINITY) ? FLT_MIN : 1.0)


/*****************************************************************
 * 1. OA sweep and trace
 *****************************************************************/

/* oa_row()
 * OA row <i> (1..L) into dp[i] and xmx row i, from row i-1 and the
 * posterior row <pp> / <ppx> of row i. Same calculation as
 * p7_GOptimalAccuracy().
 */
static void
oa_row(const P7_PROFILE *gm, const float *pp, const float *ppx, float **dp, float *xmx, int i)
{
  float const *tsc = gm->tsc;
  int          M   = gm->M;
  float        esc = p7_profile_IsLocal(gm) ? 1.0 : 0.0;
  float        t1, t2;
  int          k;

  MMX(i,0) = IMX(i,0) = DMX(i,0) = XMX(i,p7G_E) = -eslINFINITY;

  for (k = 1; k < M; k++)
    {
      MMX(i,k)     = ESL_MAX(ESL_MAX(TSCDELTA(p7P_MM, k-1) * (MMX(i-1,k-1)  + pp[k*p7G_NSCELLS + p7G_M]),
                                     TSCDELTA(p7P_IM, k-1) * (IMX(i-1,k-1)  + pp[k*p7G_NSCELLS + p7G_M])),
                             ESL_MAX(TSCDELTA(p7P_DM, k-1) * (DMX(i-1,k-1)  + pp[k*p7G_NSCELLS + p7G_M]),
                                     TSCDELTA(p7P_BM, k-1) * (XMX(i-1,p7G_B)+ pp[k*p7G_NSCELLS + p7G_M])));

      XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), esc * MMX(i,k));

      IMX(i,k)     = ESL_MAX(TSCDELTA(p7P_MI, k) * (MMX(i-1,k) + pp[k*p7G_NSCELLS + p7G_I]),
                             TSCDELTA(p7P_II, k) * (IMX(i-1,k) + pp[k*p7G_NSCELLS + p7G_I]));

      DMX(i,k)     = ESL_MAX(TSCDELTA(p7P_MD, k-1) * MMX(i,k-1),
                             TSCDELTA(p7P_DD, k-1) * DMX(i,k-1));
    }

  MMX(i,M)     = ESL_MAX(ESL_MAX(TSCDELTA(p7P_MM, M-1) * (MMX(i-1,M-1)  + pp[M*p7G_NSCELLS + p7G_M]),
                                 TSCDELTA(p7P_IM, M-1) * (IMX(i-1,M-1)  + pp[M*p7G_NSCELLS + p7G_M])),
                         ESL_MAX(TSCDELTA(p7P_DM, M-1) * (DMX(i-1,M-1)  + pp[M*p7G_NSCELLS + p7G_M]),
                                 TSCDELTA(p7P_BM, M-1) * (XMX(i-1,p7G_B)+ pp[M*p7G_NSCELLS + p7G_M])));
  IMX(i,M)     = -eslINFINITY;
  DMX(i,M)     = ESL_MAX(TSCDELTA(p7P_MD, M-1) * MMX(i,M-1),
                         TSCDELTA(p7P_DD, M-1) * DMX(i,M-1));
  XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), ESL_MAX(MMX(i,M), DMX(i, M)));

  t1 = ( (gm->xsc[p7P_J][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  t2 = ( (gm->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  XMX(i, p7G_J) = ESL_MAX( t1 * (XMX(i-1,p7G_J) + ppx[p7G_J]),
                           t2 * XMX(i,  p7G_E));

  t1 = ( (gm->xsc[p7P_C][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  t2 = ( (gm->xsc[p7P_E][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  XMX(i,p7G_C) = ESL_MAX( t1 * (XMX(i-1,p7G_C) + ppx[p7G_C]),
                          t2 * XMX(i,  p7G_E));

  t1 = ( (gm->xsc[p7P_N][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  XMX(i,p7G_N) = t1 *  (XMX(i-1,p7G_N) + ppx[p7G_N]);

  t1 = ( (gm->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  t2 = ( (gm->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  XMX(i,p7G_B) = ESL_MAX( t1 * XMX(i,  p7G_N),
                          t2 * XMX(i,  p7G_J));
}

/* oa_row_0()
 * OA row 0 into dp[0] and xmx row 0, as p7_GOptimalAccuracy() sets it.
 */
static void
oa_row_0(const P7_PROFILE *gm, float **dp, float *xmx)
{
  int k;

  XMX(0,p7G_N) = 0.;
  XMX(0,p7G_B) = 0.;
  XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) = -eslINFINITY;
  for (k = 0; k <= gm->M; k++)
    MMX(0,k) = IMX(0,k) = DMX(0,k) = -eslINFINITY;
}

/* oat_store_row()
 * Store into <ox> what the OA trace reads from OA row <i>: for each k,
 * the predecessors of an M at (i+1,k), an I at (i+1,k) and a D at
 * (i,k), chosen as p7_GOATrace() chooses them from the OA matrix; the
 * posterior codes of M and I at (i,k) from the posterior row <pp>
 * (all 0 for row 0); the row's OA special states and posteriors of N,
 * J, C (<ppx>, NULL for row 0); and the cell an E at row i comes from.
 */
static void
oat_store_row(const P7_PROFILE *gm, float **dp, const float *xmx, const float *pp, const float *ppx, P7_OATMX *ox, int i)
{
  static const int from[4] = { p7OAT_FROM_M, p7OAT_FROM_I, p7OAT_FROM_D, p7OAT_FROM_B };
  float const *tsc  = gm->tsc;
  int          M    = gm->M;
  float        path[4];
  float        max;
  int          smax, kmax;
  int          k;

  for (k = 1; k <= M; k++)
    {
      path[0] = TSCDELTA(p7P_MM, k-1) * MMX(i,k-1);
      path[1] = TSCDELTA(p7P_IM, k-1) * IMX(i,k-1);
      path[2] = TSCDELTA(p7P_DM, k-1) * DMX(i,k-1);
      path[3] = TSCDELTA(p7P_BM, k-1) * XMX(i,p7G_B);
      OATX_SET_MPRED(ox, i, k, from[esl_vec_FArgMax(path, 4)]);

      if (k < M)
        {
          path[0] = TSCDELTA(p7P_MI, k) * MMX(i,k);
          path[1] = TSCDELTA(p7P_II, k) * IMX(i,k);
          OATX_SET_IPRED(ox, i, k, (path[0] >= path[1]) ? p7OAT_FROM_M : p7OAT_FROM_I);
        }
      else OATX_SET_IPRED(ox, i, k, p7OAT_FROM_M);

      path[0] = TSCDELTA(p7P_MD, k-1) * MMX(i,k-1);
      path[1] = TSCDELTA(p7P_DD, k-1) * DMX(i,k-1);
      OATX_SET_DPRED(ox, i, k, (path[0] >= path[1]) ? 0 : 1);

      OATX_SET_PPM(ox, i, k, pp ? p7OAT_PPCODE(pp[k*p7G_NSCELLS + p7G_M]) : 0);
      OATX_SET_PPI(ox, i, k, pp ? p7OAT_PPCODE(pp[k*p7G_NSCELLS + p7G_I]) : 0);
    }

  OATX_XMX(ox, i, p7OAT_N)   = XMX(i, p7G_N);
  OATX_XMX(ox, i, p7OAT_J)   = XMX(i, p7G_J);
  OATX_XMX(ox, i, p7OAT_C)   = XMX(i, p7G_C);
  OATX_XMX(ox, i, p7OAT_B)   = XMX(i, p7G_B);
  OATX_XMX(ox, i, p7OAT_E)   = XMX(i, p7G_E);
  OATX_XMX(ox, i, p7OAT_PPN) = ppx ? ppx[p7G_N] : 0.0;
  OATX_XMX(ox, i, p7OAT_PPJ) = ppx ? ppx[p7G_J] : 0.0;
  OATX_XMX(ox, i, p7OAT_PPC) = ppx ? ppx[p7G_C] : 0.0;

  /* the cell an E at row i comes from, as p7_GOATrace()'s select_e() picks it */
  if (! p7_profile_IsLocal(gm))
    ox->ek[i] = (MMX(i,M) >= DMX(i,M)) ? M : -M;
  else
    {
      max  = -eslINFINITY;
      smax = -1;
      kmax = 0;
      for (k = 1; k <= M; k++)
        {
          if (MMX(i,k) >= max) { max = MMX(i,k); smax = p7T_M; kmax = k; }
          if (DMX(i,k) >  max) { max = DMX(i,k); smax = p7T_D; kmax = k; }
        }
      ox->ek[i] = (smax == -1) ? 0 : ((smax == p7T_M) ? kmax : -kmax);
    }
}


/* Function:  p7_GOptimalAccuracySweep()
 * Synopsis:  OA alignment in one forward sweep over a checkpointed Backward.
 *
 * Purpose:   For profile <gm> and sequence <dsq> of length <L>, with the
 *            checkpointed Backward matrix <cx> filled by
 *            <p7_GBackwardChk()>, compute in one sweep through rows 1..L
 *            the Forward, posterior and optimal accuracy (OA) rows,
 *            computing each block of Backward rows again from its
 *            checkpoint as the sweep reaches it, and store in <ox> what
 *            <p7_GOATraceSweep()> needs for the OA trace.
 *
 *            The posteriors are normalized by the Backward score, as
 *            <p7_GDecoding()> does. Rows of Forward, posterior and OA
 *            values are kept only as long as the next row needs them.
 *
 *            If <opt_null2> is non-NULL, the null2 odds ratios of the
 *            sequence are also calculated into it, as
 *            <p7_GNull2_ByExpectation()> calculates them from the
 *            posteriors: the posterior rows are summed as the sweep
 *            computes them, in the same order.
 *
 *            <ox> must have been created with <p7OAT_NPLANES> planes and
 *            1 lane (the generic layout); it is grown for <gm->M> by <L>.
 *            On return <cx> has its last block loaded.
 *
 * Args:      dsq      - sequence in digitized form, 1..L
 *            L        - length of dsq (>= 1)
 *            gm       - profile
 *            cx       - checkpointed Backward matrix from p7_GBackwardChk()
 *            ox       - RESULT: what the OA trace reads
 *            opt_null2- optRETURN: null2 odds ratios, <0..Kp-1>; caller allocated
 *            opt_fsc  - optRETURN: Forward score in nats
 *            ret_oasc - RETURN: OA score, the expected number of
 *                       correctly aligned residues
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; <eslEINVAL> if <ox> or <cx>
 *            isn't set up for this computation.
 */
int
p7_GOptimalAccuracySweep(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_CHKGMX *cx, P7_OATMX *ox, float *opt_null2, float *opt_fsc, float *ret_oasc)
{
  int     M     = gm->M;
  int     W     = (M+1) * p7G_NSCELLS;
  float  *mem   = NULL;
  float **fdp   = NULL;
  float **odp   = NULL;
  float  *fxmx  = NULL;
  float  *oxmx  = NULL;
  float  *pp;
  float  *cnt;
  float   ppx[p7G_NXCELLS];
  float   cntx[p7G_NXCELLS];
  float   overall_sc;
  int     i, b;
  int     status;

  if (ox->lanes != 1 || ox->nplanes != p7OAT_NPLANES) ESL_EXCEPTION(eslEINVAL, "P7_OATMX not set up for the generic standard OA sweep");
  if (cx->L != L || cx->M != M || cx->dir != p7_CHK_BACKWARD) ESL_EXCEPTION(eslEINVAL, "checkpointed Backward matrix doesn't match");
  if ((status = p7_oatmx_GrowTo(ox, M, L)) != eslOK) return status;

  /* two Forward rows, two OA rows, one posterior row and the posterior sums for null2; row pointers and special states for all rows */
  ESL_ALLOC(mem,  sizeof(float)   * W * 6);
  ESL_ALLOC(fdp,  sizeof(float *) * (L+1));
  ESL_ALLOC(odp,  sizeof(float *) * (L+1));
  ESL_ALLOC(fxmx, sizeof(float)   * (L+1) * p7G_NXCELLS);
  ESL_ALLOC(oxmx, sizeof(float)   * (L+1) * p7G_NXCELLS);
  pp  = mem + 4*W;
  cnt = mem + 5*W;

  p7_FLogsumInit();

  overall_sc = cx->xmx[p7G_N];   /* Backward score */

  /* the sweep */
  fdp[0] = mem;
  odp[0] = mem + 2*W;
  p7_GForwardRow(dsq, gm, fdp, fxmx, 0);
  oa_row_0(gm, odp, oxmx);
  oat_store_row(gm, odp, oxmx, NULL, NULL, ox, 0);

  for (i = 1; i <= L; i++)
    {
      fdp[i] = mem + (i%2)*W;
      odp[i] = mem + (2 + i%2)*W;

      b = (i / cx->blockR) * cx->blockR;
      if (b != cx->b && (status = p7_GBackwardChkLoad(dsq, gm, cx, b)) != eslOK) goto ERROR;

      p7_GForwardRow(dsq, gm, fdp, fxmx, i);
      p7_GDecodingRow(gm, fdp[i], fxmx + (i-1)*p7G_NXCELLS, cx->dp[i], cx->xmx + i*p7G_NXCELLS, overall_sc, pp, ppx);
      oa_row(gm, pp, ppx, odp, oxmx, i);
      oat_store_row(gm, odp, oxmx, pp, ppx, ox, i);

      /* posterior sums for null2, added up as p7_GNull2_ByExpectation() does */
      if (i == 1) { esl_vec_FCopy(pp, W, cnt); esl_vec_FCopy(ppx, p7G_NXCELLS, cntx); }
      else        { esl_vec_FAdd (cnt, pp, W); esl_vec_FAdd (cntx, ppx, p7G_NXCELLS); }
    }

  if (opt_null2 != NULL) p7_GNull2_FromCounts(gm, cnt, cntx, L, opt_null2);
  if (opt_fsc   != NULL) *opt_fsc = fxmx[L*p7G_NXCELLS + p7G_C] + gm->xsc[p7P_C][p7P_MOVE];
  *ret_oasc = oxmx[L*p7G_NXCELLS + p7G_C];

  free(mem); free(fdp); free(odp); free(fxmx); free(oxmx);
  return eslOK;

 ERROR:
  free(mem); free(fdp); free(odp); free(fxmx); free(oxmx);
  return status;
}


/* Function:  p7_GOATraceSweep()
 * Synopsis:  OA trace from what the OA sweep stored.
 *
 * Purpose:   The OA trace for profile <gm> from <ox>, filled by
 *            <p7_GOptimalAccuracySweep()>, into the empty trace <tr>
 *            (created with <p7_trace_CreateWithPP()>). The trace is
 *            identical to <p7_GOATrace()>'s from the full OA matrix.
 *            The posterior probability of each residue is the value its
 *            PP character stands for (<p7OAT_PPVALUE()>), from the
 *            posterior matrix, as <p7_alidisplay_EncodePostProb()>
 *            prints it.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if a choice
 *            fails or <tr> isn't empty.
 */
int
p7_GOATraceSweep(const P7_PROFILE *gm, const P7_OATMX *ox, P7_TRACE *tr)
{
  static const int state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };   /* indexed by p7OAT_FROM_* */
  int   i = ox->L;
  int   k = 0;
  float t1, t2;
  float postprob;
  int   sprv, scur;
  int   status;

  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");

  if ((status = p7_trace_AppendWithPP(tr, p7T_T, k, i, 0.0)) != eslOK) return status;
  if ((status = p7_trace_AppendWithPP(tr, p7T_C, k, i, 0.0)) != eslOK) return status;

  sprv = p7T_C;
  while (sprv != p7T_S)
    {
      switch (sprv) {
      case p7T_M: scur = state[OATX_MPRED(ox, i-1, k)];               k--; i--; break;
      case p7T_D: scur = OATX_DPRED(ox, i, k) ? p7T_D : p7T_M;         k--;      break;
      case p7T_I: scur = (OATX_IPRED(ox, i-1, k) == p7OAT_FROM_I) ? p7T_I : p7T_M; i--; break;
      case p7T_N: scur = (i == 0) ? p7T_S : p7T_N;                               break;
      case p7T_C:
        t2   = ( (gm->xsc[p7P_E][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
        scur = (OATX_XMX(ox, i, p7OAT_C) > t2 * OATX_XMX(ox, i, p7OAT_E)) ? p7T_C : p7T_E;
        break;
      case p7T_J:
        t2   = ( (gm->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
        scur = (OATX_XMX(ox, i, p7OAT_J) > t2 * OATX_XMX(ox, i, p7OAT_E)) ? p7T_J : p7T_E;
        break;
      case p7T_E:
        if (ox->ek[i] == 0) ESL_EXCEPTION(eslEINVAL, "OA traceback choice failed");
        scur = (ox->ek[i] > 0) ? p7T_M : p7T_D;
        k    = abs(ox->ek[i]);
        break;
      case p7T_B:
        t1   = ( (gm->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
        t2   = ( (gm->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
        scur = (t1 * OATX_XMX(ox, i, p7OAT_N) > t2 * OATX_XMX(ox, i, p7OAT_J)) ? p7T_N : p7T_J;
        break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }

      switch (scur) {
      case p7T_M: postprob = p7OAT_PPVALUE(OATX_PPM(ox, i, k));                          break;
      case p7T_I: postprob = p7OAT_PPVALUE(OATX_PPI(ox, i, k));                          break;
      case p7T_N: postprob = (sprv == scur) ? OATX_XMX(ox, i, p7OAT_PPN) : 0.0;          break;
      case p7T_C: postprob = (sprv == scur) ? OATX_XMX(ox, i, p7OAT_PPC) : 0.0;          break;
      case p7T_J: postprob = (sprv == scur) ? OATX_XMX(ox, i, p7OAT_PPJ) : 0.0;          break;
      default:    postprob = 0.0;                                                        break;
      }
      if ((status = p7_trace_AppendWithPP(tr, scur, k, i, postprob)) != eslOK) return status;

      if ( (scur == p7T_N || scur == p7T_J || scur == p7T_C) && scur == sprv) i--;
      sprv = scur;
    }
  tr->M = gm->M;
  tr->L = ox->L;
  return p7_trace_Reverse(tr);
}


/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef p7GENERIC_OPTACC_SWEEP_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_sweep()
 * The OA sweep and its trace against the full-matrix path
 * (p7_GForward, p7_GBackward, p7_GDecoding, p7_GOptimalAccuracy,
 * p7_GOATrace) for random sequences and sequences emitted from the
 * model, with Backward blocks of 1..8 rows and of the whole sequence:
 * same OA score and null2, same states, k and i on the trace, and each
 * residue's PP character equal to the one the posterior matrix gives it.
 */
static void
utest_sweep(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N, int mode)
{
  char       *msg  = "generic OA sweep unit test failed";
  P7_HMM     *hmm  = NULL;
  P7_PROFILE *gm   = p7_profile_Create(M, abc);
  ESL_SQ     *sq   = esl_sq_CreateDigital(abc);
  ESL_DSQ    *dsq  = NULL;
  P7_GMX     *fwd  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX     *bck  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX     *pp   = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX     *oa   = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_CHKGMX  *cx   = p7_chkgmx_Create(M, L, 1, 1, p7G_NSCELLS, p7_CHK_BACKWARD);
  P7_OATMX   *ox   = p7_oatmx_Create(M, L, p7OAT_NPLANES, 1);
  P7_TRACE   *tr1  = p7_trace_CreateWithPP();
  P7_TRACE   *tr2  = p7_trace_CreateWithPP();
  float      *null2a = malloc(sizeof(float) * abc->Kp);
  float      *null2b = malloc(sizeof(float) * abc->Kp);
  float       fsc, fsc2, oasc1, oasc2, p;
  int         n, Ln, z, blockR, x;

  if (p7_hmm_Sample(r, M, abc, &hmm)           != eslOK) esl_fatal(msg);
  if (p7_ProfileConfig(hmm, bg, gm, L, mode)   != eslOK) esl_fatal(msg);

  for (n = 0; n < N; n++)
    {
      if (n % 2 == 0)
        {
          Ln  = 1 + esl_rnd_Roll(r, L);
          dsq = realloc(dsq, sizeof(ESL_DSQ) * (Ln+2));
          esl_rsq_xfIID(r, bg->f, abc->K, Ln, dsq);
        }
      else
        {
          do { esl_sq_Reuse(sq); p7_ProfileEmit(r, hmm, gm, bg, sq, NULL); } while (sq->n < 1 || sq->n > L);
          Ln  = sq->n;
          dsq = realloc(dsq, sizeof(ESL_DSQ) * (Ln+2));
          memcpy(dsq, sq->dsq, sizeof(ESL_DSQ) * (Ln+2));
        }
      p7_ReconfigLength(gm, Ln);

      /* full-matrix path */
      p7_GForward (dsq, Ln, gm, fwd, &fsc);
      p7_GBackward(dsq, Ln, gm, bck, NULL);
      p7_GDecoding(gm, fwd, bck, pp);
      p7_GOptimalAccuracy(gm, pp, oa, &oasc1);
      p7_trace_Reuse(tr1);
      if (p7_GOATrace(gm, oa, tr1) != eslOK) esl_fatal(msg);
      p7_GNull2_ByExpectation(gm, pp, null2a);                   /* uses row 0 of <pp> as workspace */

      /* sweep */
      blockR = (n % 3 == 2) ? Ln+1 : 1 + esl_rnd_Roll(r, 8);
      if (p7_chkgmx_GrowTo(cx, M, Ln, blockR)             != eslOK) esl_fatal(msg);
      if (p7_GBackwardChk(dsq, Ln, gm, cx, NULL)           != eslOK) esl_fatal(msg);
      if (p7_GOptimalAccuracySweep(dsq, Ln, gm, cx, ox, null2b, &fsc2, &oasc2) != eslOK) esl_fatal(msg);
      p7_trace_Reuse(tr2);
      if (p7_GOATraceSweep(gm, ox, tr2) != eslOK) esl_fatal(msg);

      if (fsc2  != fsc)   esl_fatal("%s: Forward scores differ (%f %f)", msg, fsc2, fsc);
      if (oasc2 != oasc1) esl_fatal("%s: OA scores differ (%f %f)", msg, oasc2, oasc1);
      for (x = 0; x < abc->Kp; x++)
        if (null2b[x] != null2a[x]) esl_fatal("%s: null2 differs for residue %d (%g %g)", msg, x, null2b[x], null2a[x]);
      if (tr1->N != tr2->N) esl_fatal("%s: trace lengths differ (%d %d)", msg, tr1->N, tr2->N);
      for (z = 0; z < tr1->N; z++)
        {
          if (tr1->st[z] != tr2->st[z] || tr1->k[z] != tr2->k[z] || tr1->i[z] != tr2->i[z])
            esl_fatal("%s: traces differ at %d", msg, z);
          switch (tr1->st[z]) {
          case p7T_M: p = pp->dp[tr1->i[z]][tr1->k[z]*p7G_NSCELLS + p7G_M]; break;
          case p7T_I: p = pp->dp[tr1->i[z]][tr1->k[z]*p7G_NSCELLS + p7G_I]; break;
          case p7T_N: case p7T_C: case p7T_J:
            /* an N, C or J emits residue i when the state before it is the same one */
            p = (z > 0 && tr1->st[z-1] == tr1->st[z]) ? pp->xmx[tr1->i[z]*p7G_NXCELLS + (tr1->st[z] == p7T_N ? p7G_N : tr1->st[z] == p7T_C ? p7G_C : p7G_J)] : 0.0;
            break;
          default: p = 0.0; break;
          }
          if (p7_alidisplay_EncodePostProb(tr2->pp[z]) != p7_alidisplay_EncodePostProb(p))
            esl_fatal("%s: PP character differs at %d (state %d, i %d, k %d)", msg, z, tr1->st[z], tr1->i[z], tr1->k[z]);
        }
    }

  free(dsq);
  free(null2a);
  free(null2b);
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr1);
  p7_trace_Destroy(tr2);
  p7_oatmx_Destroy(ox);
  p7_chkgmx_Destroy(cx);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(bck);
  p7_gmx_Destroy(pp);
  p7_gmx_Destroy(oa);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}

#endif /*p7GENERIC_OPTACC_SWEEP_TESTDRIVE*/


/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef p7GENERIC_OPTACC_SWEEP_TESTDRIVE
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name  type         default  env   range togs  reqs  incomp  help                                   docgroup */
  { "-h",  eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",        0 },
  { "-s",  eslARG_INT,     "42", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",               0 },
  { "-L",  eslARG_INT,    "100", NULL, NULL, NULL, NULL, NULL, "maximum length of sequences",                 0 },
  { "-M",  eslARG_INT,     "50", NULL, NULL, NULL, NULL, NULL, "length of sampled models",                    0 },
  { "-N",  eslARG_INT,     "40", NULL, NULL, NULL, NULL, NULL, "number of sequences per test",                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for the generic OA sweep";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r   = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg  = p7_bg_Create(abc);
  int             M   = esl_opt_GetInteger(go, "-M");
  int             L   = esl_opt_GetInteger(go, "-L");
  int             N   = esl_opt_GetInteger(go, "-N");

  p7_FLogsumInit();

  utest_sweep(r, abc, bg, M, L, N,  p7_LOCAL);
  utest_sweep(r, abc, bg, M, L, N,  p7_UNILOCAL);
  utest_sweep(r, abc, bg, M, L, N,  p7_GLOCAL);
  utest_sweep(r, abc, bg, 1, L, 10, p7_LOCAL);
  utest_sweep(r, abc, bg, M, 1, 10, p7_LOCAL);


  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_OPTACC_SWEEP_TESTDRIVE*/
