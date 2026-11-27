/* Optimal accuracy alignment in one sweep over a checkpointed Backward
 * matrix, and its trace; standard model, NEON implementation.
 *
 * Contents:
 *    1. p7_OptimalAccuracySweep(), p7_OATraceSweep()
 *
 * Forward rows are computed by p7_ForwardRow(), posteriors by
 * p7_DecodingRow() and OA rows by p7_OptimalAccuracyRow() (optacc.c).
 */
#include <p7_config.h>

#include <float.h>
#include <string.h>

#include <arm_neon.h>

#include "easel.h"
#include "esl_neon.h"
#include "esl_vectorops.h"

#include "hmmer.h"


/*****************************************************************
 * 1. OA sweep over a checkpointed Backward, and its trace
 *****************************************************************/

/* oat_movemask()
 * The 4-bit lane mask of comparison result <m>: bit r set when lane r is.
 */
static inline int
oat_movemask(uint32x4_t m)
{
  return (int) ((vgetq_lane_u32(m, 0) & 1) | (vgetq_lane_u32(m, 1) & 2) | (vgetq_lane_u32(m, 2) & 4) | (vgetq_lane_u32(m, 3) & 8));
}

/* oat_ppcodes4()
 * Posterior codes (p7OAT_PPCODE(), 0..10) of the 4 lanes of <v> as 4
 * lane masks, bit b of each lane's code in <m[b]>.
 */
static inline void
oat_ppcodes4(float32x4_t v, int m[4])
{
  union { float32x4_t v; float p[4]; } u;
  int r, c;

  u.v  = v;
  m[0] = m[1] = m[2] = m[3] = 0;
  for (r = 0; r < 4; r++)
    {
      c     = p7OAT_PPCODE(u.p[r]);
      m[0] |= ( c       & 1) << r;
      m[1] |= ((c >> 1) & 1) << r;
      m[2] |= ((c >> 2) & 1) << r;
      m[3] |= ((c >> 3) & 1) << r;
    }
}

/* oat_store_row()
 * Store into <ox> what the OA trace reads from OA row <i> (<dpc>, with
 * the OA special states in <xmx>): for each k, the predecessors of an
 * M at (i+1,k), an I at (i+1,k) and a D at (i,k), chosen with the same
 * comparisons, shifts and ties as p7_OATrace() makes on the OA
 * matrix; the posterior codes of M and I at (i,k) from the posterior
 * row <ppp> (all 0 for row 0, <ppp> NULL); the row's OA special states
 * and the posteriors of N, J, C (<ppx>, NULL for row 0); and the cell
 * an E at row i comes from, in p7_OATrace()'s scan order.
 */
static void
oat_store_row(const P7_OPROFILE *om, P7_OATMX *ox, const float32x4_t *dpc, const float *xmx, const float32x4_t *ppp, const float *ppx, int i)
{
  int      M     = om->M;
  int      Q     = p7O_NQF(M);
  float32x4_t   zerov = vmovq_n_f32(0.0f);
  float32x4_t   infv  = vmovq_n_f32(-eslINFINITY);
  float32x4_t   xBv   = vmovq_n_f32(xmx[i*p7X_NXCELLS+p7X_B]);
  float32x4_t   mpv, ipv, dpv, tmdv, tddv;
  float32x4_t pM, pI, pD, pB, best;
  uint32x4_t  gt, b0, b1;
  float32x4_t  *tp;
  union { float32x4_t v; float p[4]; } u;
  float    max;
  int      kmax, smax;
  int      mm[4], mi[4];
  int      q, r, c;

  for (q = 0; q < Q; q++)
    {
      tp = om->tfv + 7*q;

      /* row i at k-1 */
      if (q > 0) {
        mpv  = MMO(dpc, q-1);
        ipv  = IMO(dpc, q-1);
        dpv  = DMO(dpc, q-1);
        tmdv = om->tfv[7*(q-1) + p7O_MD];
        tddv = om->tfv[7*Q + (q-1)];
      } else {
        mpv  = esl_neon_rightshift_float((esl_neon_128f_t) MMO(dpc, Q-1), (esl_neon_128f_t) zerov).f32x4;
        ipv  = esl_neon_rightshift_float((esl_neon_128f_t) IMO(dpc, Q-1), (esl_neon_128f_t) zerov).f32x4;
        dpv  = esl_neon_rightshift_float((esl_neon_128f_t) DMO(dpc, Q-1), (esl_neon_128f_t) zerov).f32x4;
        tmdv = esl_neon_rightshift_float((esl_neon_128f_t) om->tfv[7*(Q-1) + p7O_MD], (esl_neon_128f_t) zerov).f32x4;
        tddv = esl_neon_rightshift_float((esl_neon_128f_t) om->tfv[8*Q-1], (esl_neon_128f_t) zerov).f32x4;
      }

      /* M at (i+1,k): first best of M, I, D, B */
      pB   = vbslq_f32(vceqq_f32(tp[p7O_BM], zerov), infv, xBv);
      pM   = vbslq_f32(vceqq_f32(tp[p7O_MM], zerov), infv, mpv);
      pI   = vbslq_f32(vceqq_f32(tp[p7O_IM], zerov), infv, ipv);
      pD   = vbslq_f32(vceqq_f32(tp[p7O_DM], zerov), infv, dpv);
      best = pM;
      gt   = vcgtq_f32(pI, best);  best = vbslq_f32(gt, pI, best);  b0 = gt;                  b1 = vdupq_n_u32(0);
      gt   = vcgtq_f32(pD, best);  best = vbslq_f32(gt, pD, best);  b0 = vbicq_u32(b0, gt);    b1 = vorrq_u32(b1, gt);
      gt   = vcgtq_f32(pB, best);                                           b0 = vorrq_u32(b0, gt);    b1 = vorrq_u32(b1, gt);
      OATO_SET_STRIPE4(ox, i, p7OAT_MP, q, oat_movemask(b0));
      OATO_SET_STRIPE4(ox, i, p7OAT_MP+1, q, oat_movemask(b1));

      /* I at (i+1,k): M on ties */
      pM = vbslq_f32(vceqq_f32(tp[p7O_MI], zerov), infv, MMO(dpc, q));
      pI = vbslq_f32(vceqq_f32(tp[p7O_II], zerov), infv, IMO(dpc, q));
      OATO_SET_STRIPE4(ox, i, p7OAT_IP, q, oat_movemask(vcgtq_f32(pI, pM)));

      /* D at (i,k): M on ties */
      pM = vbslq_f32(vceqq_f32(tmdv, zerov), infv, mpv);
      pD = vbslq_f32(vceqq_f32(tddv, zerov), infv, dpv);
      OATO_SET_STRIPE4(ox, i, p7OAT_DP, q, oat_movemask(vcgtq_f32(pD, pM)));

      /* posterior codes of M and I */
      if (ppp != NULL) { oat_ppcodes4(ppp[q*3 + p7X_M], mm); oat_ppcodes4(ppp[q*3 + p7X_I], mi); }
      else             { mm[0] = mm[1] = mm[2] = mm[3] = mi[0] = mi[1] = mi[2] = mi[3] = 0; }
      for (c = 0; c < 4; c++)
        {
          OATO_SET_STRIPE4(ox, i, p7OAT_PPM + c, q, mm[c]);
          OATO_SET_STRIPE4(ox, i, p7OAT_PPI + c, q, mi[c]);
        }
    }

  ox->xmx[i*p7OAT_NXCELLS + p7OAT_N]   = xmx[i*p7X_NXCELLS + p7X_N];
  ox->xmx[i*p7OAT_NXCELLS + p7OAT_J]   = xmx[i*p7X_NXCELLS + p7X_J];
  ox->xmx[i*p7OAT_NXCELLS + p7OAT_C]   = xmx[i*p7X_NXCELLS + p7X_C];
  ox->xmx[i*p7OAT_NXCELLS + p7OAT_B]   = xmx[i*p7X_NXCELLS + p7X_B];
  ox->xmx[i*p7OAT_NXCELLS + p7OAT_E]   = xmx[i*p7X_NXCELLS + p7X_E];
  ox->xmx[i*p7OAT_NXCELLS + p7OAT_PPN] = ppx ? ppx[p7X_N] : 0.0f;
  ox->xmx[i*p7OAT_NXCELLS + p7OAT_PPJ] = ppx ? ppx[p7X_J] : 0.0f;
  ox->xmx[i*p7OAT_NXCELLS + p7OAT_PPC] = ppx ? ppx[p7X_C] : 0.0f;

  /* the cell an E at row i comes from, as oa_select_e() scans the row */
  max  = -eslINFINITY;
  smax = p7T_M;
  kmax = 1;
  for (q = 0; q < Q; q++)
    {
      u.v = MMO(dpc, q); for (r = 0; r < 4; r++) if (u.p[r] >= max) { max = u.p[r]; smax = p7T_M; kmax = r*Q + q + 1; }
      u.v = DMO(dpc, q); for (r = 0; r < 4; r++) if (u.p[r] >  max) { max = u.p[r]; smax = p7T_D; kmax = r*Q + q + 1; }
    }
  ox->ek[i] = (smax == p7T_M) ? kmax : -kmax;
}


/* Function:  p7_OptimalAccuracySweep()
 * Synopsis:  OA alignment in one forward sweep over a checkpointed Backward, NEON.
 *
 * Purpose:   For profile <om> and sequence <dsq> of length <L>, with the
 *            checkpointed Backward matrix <cx> filled by
 *            <p7_BackwardChk()>, compute in one sweep through rows
 *            1..L the Forward, posterior and optimal accuracy (OA) rows,
 *            computing each block of Backward rows again from its
 *            checkpoint as the sweep reaches it, and store in <ox> what
 *            <p7_OATraceSweep()> needs for the OA trace. Rows of
 *            Forward, posterior and OA values are kept only as long as
 *            the next row needs them.
 *
 *            Each row is calculated as <p7_Forward()>,
 *            <p7_Decoding()> and <p7_OptimalAccuracy()> calculate
 *            it. The checkpointed Backward has its own scale factors, so
 *            the posteriors are those <p7_Decoding()> gives for a
 *            Backward matrix with its own scale factors: normalized by
 *            the Backward score, and equal to those from a Backward
 *            scaled like Forward up to rounding.
 *
 *            If <opt_null2> is non-NULL, the null2 odds ratios of the
 *            sequence are also calculated into it, as
 *            <p7_Null2_ByExpectation()> calculates them from the
 *            posteriors: the posterior rows are summed as the sweep
 *            computes them, in the same order.
 *
 *            <ox> must have been created with <p7OAT_NPLANES> planes;
 *            its layout is set to 4 lanes and it is grown for <om->M> by
 *            <L>. On return <cx> has its last block loaded.
 *
 * Args:      dsq       - sequence in digitized form, 1..L
 *            L         - length of dsq (>= 1)
 *            om        - profile
 *            cx        - checkpointed Backward matrix from p7_BackwardChk()
 *            ox        - RESULT: what the OA trace reads
 *            opt_null2 - optRETURN: null2 odds ratios, <0..Kp-1>; caller allocated
 *            opt_fsc   - optRETURN: Forward score in nats
 *            ret_oasc  - RETURN: OA score, the expected number of
 *                        correctly aligned residues
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; <eslEINVAL> if <ox> or <cx>
 *            isn't set up for this computation; <eslERANGE> if the
 *            posterior scale correction overflows (as <p7_Decoding()>).
 */
int
p7_OptimalAccuracySweep(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, P7_OATMX *ox,
                            float *opt_null2, float *opt_fsc, float *ret_oasc)
{
  int      M      = om->M;
  int      Q      = p7O_NQF(M);
  P7_OMX  *wk     = NULL;   /* rows 0,1: Forward; 2,3: OA; 4: posteriors; 5: posterior sums; xmx: Forward specials */
  float   *oxmx   = NULL;
  float    ppx[p7X_NXCELLS];
  float    cntx[p7X_NXCELLS];
  float    totscale = 0.0f;
  float    scaleproduct;
  float32x4_t  *pp, *cnt;
  int      i, q, b;
  int      status;

  if (ox->nplanes != p7OAT_NPLANES) ESL_EXCEPTION(eslEINVAL, "P7_OATMX not set up for the standard OA sweep");
  if (cx->L != L || cx->M != M || cx->dir != p7_CHK_BACKWARD || cx->nscells != p7X_NSCELLS) ESL_EXCEPTION(eslEINVAL, "checkpointed Backward matrix doesn't match");
  ox->lanes = 4;
  if ((status = p7_oatmx_GrowTo(ox, M, L)) != eslOK) return status;

  if ((wk = p7_omx_Create(M, 5, L)) == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(oxmx, sizeof(float) * (L+1) * p7X_NXCELLS);
  pp  = wk->dpf[4];
  cnt = wk->dpf[5];

  scaleproduct = 1.0 / cx->xmx[p7X_N];
  for (q = 0; q < p7X_NXCELLS; q++) cntx[q] = 0.0f;

  p7_ForwardRow(dsq, om, wk->xmx, &totscale, wk->dpf[0], NULL, 0);
  p7_OptimalAccuracyRow(om, NULL, NULL, oxmx, wk->dpf[2], NULL, 0);
  oat_store_row(om, ox, wk->dpf[2], oxmx, NULL, NULL, 0);

  for (i = 1; i <= L; i++)
    {
      b = (i / cx->blockR) * cx->blockR;
      if (b != cx->b && (status = p7_BackwardChkLoad(dsq, om, cx, b)) != eslOK) goto ERROR;

      p7_ForwardRow(dsq, om, wk->xmx, &totscale, wk->dpf[i%2], wk->dpf[(i-1)%2], i);
      p7_DecodingRow(om, wk->xmx, cx->xmx, wk->dpf[i%2], cx->dpf[i], pp, ppx, i, scaleproduct);
      scaleproduct *= wk->xmx[i*p7X_NXCELLS+p7X_SCALE] / cx->xmx[i*p7X_NXCELLS+p7X_SCALE];

      p7_OptimalAccuracyRow(om, pp, ppx, oxmx, wk->dpf[2 + i%2], wk->dpf[2 + (i-1)%2], i);
      oat_store_row(om, ox, wk->dpf[2 + i%2], oxmx, pp, ppx, i);

      /* posterior sums for null2, added up as p7_Null2_ByExpectation() does */
      if (i == 1)
        {
          memcpy(cnt, pp, sizeof(float32x4_t) * 3 * Q);
          cntx[p7X_N] = ppx[p7X_N];
          cntx[p7X_C] = ppx[p7X_C];
          cntx[p7X_J] = ppx[p7X_J];
        }
      else
        {
          for (q = 0; q < Q; q++)
            {
              cnt[q*3 + p7X_M] = vaddq_f32(pp[q*3 + p7X_M], cnt[q*3 + p7X_M]);
              cnt[q*3 + p7X_I] = vaddq_f32(pp[q*3 + p7X_I], cnt[q*3 + p7X_I]);
            }
          cntx[p7X_N] += ppx[p7X_N];
          cntx[p7X_C] += ppx[p7X_C];
          cntx[p7X_J] += ppx[p7X_J];
        }
    }
  if (isinf(scaleproduct)) { status = eslERANGE; goto ERROR; }

  if (opt_null2 != NULL) p7_Null2_FromCounts(om, cnt, cntx, L, opt_null2);
  if (opt_fsc   != NULL) *opt_fsc = totscale + log(wk->xmx[L*p7X_NXCELLS+p7X_C] * om->xf[p7O_C][p7O_MOVE]);
  *ret_oasc = oxmx[L*p7X_NXCELLS+p7X_C];

  p7_omx_Destroy(wk);
  free(oxmx);
  return eslOK;

 ERROR:
  if (wk) p7_omx_Destroy(wk);
  free(oxmx);
  return status;
}


/* Function:  p7_OATraceSweep()
 * Synopsis:  OA trace from what the NEON OA sweep stored.
 *
 * Purpose:   The OA trace for profile <om> from <ox>, filled by
 *            <p7_OptimalAccuracySweep()>, into the empty trace <tr>
 *            (created with <p7_trace_CreateWithPP()>). The choices are
 *            those <p7_OATrace()> makes on the same OA scores. The
 *            posterior probability of each residue is the value its PP
 *            character stands for (<p7OAT_PPVALUE()>), from the
 *            posterior row, as <p7_alidisplay_EncodePostProb()> prints
 *            it.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if <tr> isn't
 *            empty or a bogus state is reached.
 */
int
p7_OATraceSweep(const P7_OPROFILE *om, const P7_OATMX *ox, P7_TRACE *tr)
{
  static const int state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };   /* indexed by p7OAT_FROM_* */
  int   i = ox->L;
  int   k = 0;
  int   s0, s1;
  float path[2];
  float postprob;
  int   status;

  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace not empty; needs to be Reuse()'d?");

  if ((status = p7_trace_AppendWithPP(tr, p7T_T, k, i, 0.0)) != eslOK) return status;
  if ((status = p7_trace_AppendWithPP(tr, p7T_C, k, i, 0.0)) != eslOK) return status;

  s0 = tr->st[tr->N-1];
  while (s0 != p7T_S)
    {
      switch (s0) {
      case p7T_M: s1 = state[OATO_MPRED(ox, i-1, k)];                              k--; i--; break;
      case p7T_D: s1 = OATO_DPRED(ox, i, k) ? p7T_D : p7T_M;                        k--;      break;
      case p7T_I: s1 = (OATO_IPRED(ox, i-1, k) == p7OAT_FROM_I) ? p7T_I : p7T_M;         i--; break;
      case p7T_N: s1 = (i == 0) ? p7T_S : p7T_N;                                               break;
      case p7T_C:
        if      (om->xf[p7O_C][p7O_LOOP] == 0.0) s1 = p7T_E;
        else if (om->xf[p7O_E][p7O_MOVE] == 0.0) s1 = (ox->xmx[(i-1)*p7OAT_NXCELLS+p7OAT_C] > -eslINFINITY) ? p7T_C : p7T_E;
        else    s1 = (ox->xmx[i*p7OAT_NXCELLS+p7OAT_C] > ox->xmx[i*p7OAT_NXCELLS+p7OAT_E]) ? p7T_C : p7T_E;
        break;
      case p7T_J:
        if      (om->xf[p7O_J][p7O_LOOP] == 0.0) s1 = p7T_E;
        else if (om->xf[p7O_E][p7O_LOOP] == 0.0) s1 = (ox->xmx[(i-1)*p7OAT_NXCELLS+p7OAT_J] > -eslINFINITY) ? p7T_J : p7T_E;
        else    s1 = (ox->xmx[i*p7OAT_NXCELLS+p7OAT_J] > ox->xmx[i*p7OAT_NXCELLS+p7OAT_E]) ? p7T_J : p7T_E;
        break;
      case p7T_E:
        s1 = (ox->ek[i] > 0) ? p7T_M : p7T_D;
        k  = abs(ox->ek[i]);
        break;
      case p7T_B:
        path[0] = (om->xf[p7O_N][p7O_MOVE] == 0.0) ? -eslINFINITY : ox->xmx[i*p7OAT_NXCELLS+p7OAT_N];
        path[1] = (om->xf[p7O_J][p7O_MOVE] == 0.0) ? -eslINFINITY : ox->xmx[i*p7OAT_NXCELLS+p7OAT_J];
        s1 = (path[0] > path[1]) ? p7T_N : p7T_J;
        break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }

      switch (s1) {
      case p7T_M: postprob = p7OAT_PPVALUE(OATO_PPM(ox, i, k));                              break;
      case p7T_I: postprob = p7OAT_PPVALUE(OATO_PPI(ox, i, k));                              break;
      case p7T_N: postprob = (s1 == s0) ? ox->xmx[i*p7OAT_NXCELLS+p7OAT_PPN] : 0.0f;         break;
      case p7T_C: postprob = (s1 == s0) ? ox->xmx[i*p7OAT_NXCELLS+p7OAT_PPC] : 0.0f;         break;
      case p7T_J: postprob = (s1 == s0) ? ox->xmx[i*p7OAT_NXCELLS+p7OAT_PPJ] : 0.0f;         break;
      default:    postprob = 0.0f;                                                           break;
      }
      if ((status = p7_trace_AppendWithPP(tr, s1, k, i, postprob)) != eslOK) return status;

      if ( (s1 == p7T_N || s1 == p7T_J || s1 == p7T_C) && s1 == s0) i--;
      s0 = s1;
    }
  tr->M = om->M;
  tr->L = ox->L;
  return p7_trace_Reverse(tr);
}
