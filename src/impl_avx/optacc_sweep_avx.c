/* Optimal accuracy alignment in one sweep over a checkpointed Backward
 * matrix, and its trace; standard model, AVX2 implementation.
 *
 * Contents:
 *    1. p7_OptimalAccuracySweep_avx(), p7_OATraceSweep_avx()
 *
 * Forward rows are computed by p7_ForwardRow_avx(), posteriors by
 * p7_DecodingRow_avx() and OA rows by p7_OptimalAccuracyRow_avx() (optacc_avx.c).
 */
#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <float.h>
#include <string.h>

#include <immintrin.h>

#include "easel.h"
#include "esl_avx.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"


/*****************************************************************
 * 1. OA sweep over a checkpointed Backward, and its trace
 *****************************************************************/

/* oat_ppcodes8()
 * Posterior codes (p7OAT_PPCODE(), 0..10) of the 8 lanes of <v> as 4
 * lane masks, bit b of each lane's code in <m[b]>.
 */
static inline void
oat_ppcodes8(__m256 v, int m[4])
{
  union { __m256 v; float p[8]; } u;
  int r, c;

  u.v  = v;
  m[0] = m[1] = m[2] = m[3] = 0;
  for (r = 0; r < 8; r++)
    {
      c     = p7OAT_PPCODE(u.p[r]);
      m[0] |= ( c       & 1) << r;
      m[1] |= ((c >> 1) & 1) << r;
      m[2] |= ((c >> 2) & 1) << r;
      m[3] |= ((c >> 3) & 1) << r;
    }
}

/* oat_store_row_avx()
 * Store into <ox> what the OA trace reads from OA row <i> (<dpc>, with
 * the OA special states in <xmx>): for each k, the predecessors of an
 * M at (i+1,k), an I at (i+1,k) and a D at (i,k), chosen with the same
 * comparisons, shifts and ties as p7_OATrace_avx() makes on the OA
 * matrix; the posterior codes of M and I at (i,k) from the posterior
 * row <ppp> (all 0 for row 0, <ppp> NULL); the row's OA special states
 * and the posteriors of N, J, C (<ppx>, NULL for row 0); and the cell
 * an E at row i comes from, in p7_OATrace_avx()'s scan order.
 */
static void
oat_store_row_avx(const P7_OPROFILE *om, P7_OATMX *ox, const __m256 *dpc, const float *xmx, const __m256 *ppp, const float *ppx, int i)
{
  int      M     = om->M;
  int      Q     = p7O_NQF_AVX(M);
  __m256   zerov = _mm256_setzero_ps();
  __m256   infv  = _mm256_set1_ps(-eslINFINITY);
  __m256   xBv   = _mm256_set1_ps(xmx[i*p7X_NXCELLS+p7X_B]);
  __m256   mpv, ipv, dpv, tmdv, tddv;
  __m256   pM, pI, pD, pB, best, gt, b0, b1;
  __m256  *tp;
  union { __m256 v; float p[8]; } u;
  float    max;
  int      kmax, smax;
  int      mm[4], mi[4];
  int      q, r, c;

  for (q = 0; q < Q; q++)
    {
      tp = om->tfv_avx + 7*q;

      /* row i at k-1 */
      if (q > 0) {
        mpv  = MMO(dpc, q-1);
        ipv  = IMO(dpc, q-1);
        dpv  = DMO(dpc, q-1);
        tmdv = om->tfv_avx[7*(q-1) + p7O_MD];
        tddv = om->tfv_avx[7*Q + (q-1)];
      } else {
        mpv  = esl_avx_rightshiftz_float(MMO(dpc, Q-1));
        ipv  = esl_avx_rightshiftz_float(IMO(dpc, Q-1));
        dpv  = esl_avx_rightshiftz_float(DMO(dpc, Q-1));
        tmdv = esl_avx_rightshiftz_float(om->tfv_avx[7*(Q-1) + p7O_MD]);
        tddv = esl_avx_rightshiftz_float(om->tfv_avx[8*Q-1]);
      }

      /* M at (i+1,k): first best of M, I, D, B */
      pB   = _mm256_blendv_ps(xBv, infv, _mm256_cmp_ps(tp[p7O_BM], zerov, _CMP_EQ_OQ));
      pM   = _mm256_blendv_ps(mpv, infv, _mm256_cmp_ps(tp[p7O_MM], zerov, _CMP_EQ_OQ));
      pI   = _mm256_blendv_ps(ipv, infv, _mm256_cmp_ps(tp[p7O_IM], zerov, _CMP_EQ_OQ));
      pD   = _mm256_blendv_ps(dpv, infv, _mm256_cmp_ps(tp[p7O_DM], zerov, _CMP_EQ_OQ));
      best = pM;
      gt   = _mm256_cmp_ps(pI, best, _CMP_GT_OQ);  best = _mm256_blendv_ps(best, pI, gt);  b0 = gt;                    b1 = zerov;
      gt   = _mm256_cmp_ps(pD, best, _CMP_GT_OQ);  best = _mm256_blendv_ps(best, pD, gt);  b0 = _mm256_andnot_ps(gt, b0); b1 = _mm256_or_ps(b1, gt);
      gt   = _mm256_cmp_ps(pB, best, _CMP_GT_OQ);                                           b0 = _mm256_or_ps(b0, gt);    b1 = _mm256_or_ps(b1, gt);
      OATO_STRIPE8(ox, i, p7OAT_MP,   q) = (uint8_t) _mm256_movemask_ps(b0);
      OATO_STRIPE8(ox, i, p7OAT_MP+1, q) = (uint8_t) _mm256_movemask_ps(b1);

      /* I at (i+1,k): M on ties */
      pM = _mm256_blendv_ps(MMO(dpc, q), infv, _mm256_cmp_ps(tp[p7O_MI], zerov, _CMP_EQ_OQ));
      pI = _mm256_blendv_ps(IMO(dpc, q), infv, _mm256_cmp_ps(tp[p7O_II], zerov, _CMP_EQ_OQ));
      OATO_STRIPE8(ox, i, p7OAT_IP, q) = (uint8_t) _mm256_movemask_ps(_mm256_cmp_ps(pI, pM, _CMP_GT_OQ));

      /* D at (i,k): M on ties */
      pM = _mm256_blendv_ps(mpv, infv, _mm256_cmp_ps(tmdv, zerov, _CMP_EQ_OQ));
      pD = _mm256_blendv_ps(dpv, infv, _mm256_cmp_ps(tddv, zerov, _CMP_EQ_OQ));
      OATO_STRIPE8(ox, i, p7OAT_DP, q) = (uint8_t) _mm256_movemask_ps(_mm256_cmp_ps(pD, pM, _CMP_GT_OQ));

      /* posterior codes of M and I */
      if (ppp != NULL) { oat_ppcodes8(ppp[q*3 + p7X_M], mm); oat_ppcodes8(ppp[q*3 + p7X_I], mi); }
      else             { mm[0] = mm[1] = mm[2] = mm[3] = mi[0] = mi[1] = mi[2] = mi[3] = 0; }
      for (c = 0; c < 4; c++)
        {
          OATO_STRIPE8(ox, i, p7OAT_PPM + c, q) = (uint8_t) mm[c];
          OATO_STRIPE8(ox, i, p7OAT_PPI + c, q) = (uint8_t) mi[c];
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
      u.v = MMO(dpc, q); for (r = 0; r < 8; r++) if (u.p[r] >= max) { max = u.p[r]; smax = p7T_M; kmax = r*Q + q + 1; }
      u.v = DMO(dpc, q); for (r = 0; r < 8; r++) if (u.p[r] >  max) { max = u.p[r]; smax = p7T_D; kmax = r*Q + q + 1; }
    }
  ox->ek[i] = (smax == p7T_M) ? kmax : -kmax;
}


/* Function:  p7_OptimalAccuracySweep_avx()
 * Synopsis:  OA alignment in one forward sweep over a checkpointed Backward, AVX.
 *
 * Purpose:   For profile <om> and sequence <dsq> of length <L>, with the
 *            checkpointed Backward matrix <cx> filled by
 *            <p7_BackwardChk_avx()>, compute in one sweep through rows
 *            1..L the Forward, posterior and optimal accuracy (OA) rows,
 *            computing each block of Backward rows again from its
 *            checkpoint as the sweep reaches it, and store in <ox> what
 *            <p7_OATraceSweep_avx()> needs for the OA trace. Rows of
 *            Forward, posterior and OA values are kept only as long as
 *            the next row needs them.
 *
 *            Each row is calculated as <p7_Forward_avx()>,
 *            <p7_Decoding_avx()> and <p7_OptimalAccuracy_avx()> calculate
 *            it. The checkpointed Backward has its own scale factors, so
 *            the posteriors are those <p7_Decoding_avx()> gives for a
 *            Backward matrix with its own scale factors: normalized by
 *            the Backward score, and equal to those from a Backward
 *            scaled like Forward up to rounding.
 *
 *            If <opt_null2> is non-NULL, the null2 odds ratios of the
 *            sequence are also calculated into it, as
 *            <p7_Null2_ByExpectation_avx()> calculates them from the
 *            posteriors: the posterior rows are summed as the sweep
 *            computes them, in the same order.
 *
 *            <ox> must have been created with <p7OAT_NPLANES> planes;
 *            its layout is set to 8 lanes and it is grown for <om->M> by
 *            <L>. On return <cx> has its last block loaded.
 *
 * Args:      dsq       - sequence in digitized form, 1..L
 *            L         - length of dsq (>= 1)
 *            om        - profile
 *            cx        - checkpointed Backward matrix from p7_BackwardChk_avx()
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
 *            posterior scale correction overflows (as <p7_Decoding_avx()>).
 */
int
p7_OptimalAccuracySweep_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, P7_OATMX *ox,
                            float *opt_null2, float *opt_fsc, float *ret_oasc)
{
  int      M      = om->M;
  int      Q      = p7O_NQF_AVX(M);
  P7_OMX  *wk     = NULL;   /* rows 0,1: Forward; 2,3: OA; 4: posteriors; 5: posterior sums; xmx: Forward specials */
  float   *oxmx   = NULL;
  float    ppx[p7X_NXCELLS];
  float    cntx[p7X_NXCELLS];
  float    totscale = 0.0f;
  float    scaleproduct;
  __m256  *pp, *cnt;
  int      i, q, b;
  int      status;

  if (ox->nplanes != p7OAT_NPLANES) ESL_EXCEPTION(eslEINVAL, "P7_OATMX not set up for the standard OA sweep");
  if (cx->L != L || cx->M != M || cx->dir != p7_CHK_BACKWARD || cx->nscells != p7X_NSCELLS) ESL_EXCEPTION(eslEINVAL, "checkpointed Backward matrix doesn't match");
  ox->lanes = 8;
  if ((status = p7_oatmx_GrowTo(ox, M, L)) != eslOK) return status;

  if ((wk = p7_omx_Create_avx(M, 5, L)) == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(oxmx, sizeof(float) * (L+1) * p7X_NXCELLS);
  pp  = wk->dpf_avx[4];
  cnt = wk->dpf_avx[5];

  scaleproduct = 1.0 / cx->xmx[p7X_N];
  for (q = 0; q < p7X_NXCELLS; q++) cntx[q] = 0.0f;

  p7_ForwardRow_avx(dsq, om, wk->xmx, &totscale, wk->dpf_avx[0], NULL, 0);
  p7_OptimalAccuracyRow_avx(om, NULL, NULL, oxmx, wk->dpf_avx[2], NULL, 0);
  oat_store_row_avx(om, ox, wk->dpf_avx[2], oxmx, NULL, NULL, 0);

  for (i = 1; i <= L; i++)
    {
      b = (i / cx->blockR) * cx->blockR;
      if (b != cx->b && (status = p7_BackwardChkLoad_avx(dsq, om, cx, b)) != eslOK) goto ERROR;

      p7_ForwardRow_avx(dsq, om, wk->xmx, &totscale, wk->dpf_avx[i%2], wk->dpf_avx[(i-1)%2], i);
      p7_DecodingRow_avx(om, wk->xmx, cx->xmx, wk->dpf_avx[i%2], cx->dpf_avx[i], pp, ppx, i, scaleproduct);
      scaleproduct *= wk->xmx[i*p7X_NXCELLS+p7X_SCALE] / cx->xmx[i*p7X_NXCELLS+p7X_SCALE];

      p7_OptimalAccuracyRow_avx(om, pp, ppx, oxmx, wk->dpf_avx[2 + i%2], wk->dpf_avx[2 + (i-1)%2], i);
      oat_store_row_avx(om, ox, wk->dpf_avx[2 + i%2], oxmx, pp, ppx, i);

      /* posterior sums for null2, added up as p7_Null2_ByExpectation_avx() does */
      if (i == 1)
        {
          memcpy(cnt, pp, sizeof(__m256) * 3 * Q);
          cntx[p7X_N] = ppx[p7X_N];
          cntx[p7X_C] = ppx[p7X_C];
          cntx[p7X_J] = ppx[p7X_J];
        }
      else
        {
          for (q = 0; q < Q; q++)
            {
              cnt[q*3 + p7X_M] = _mm256_add_ps(pp[q*3 + p7X_M], cnt[q*3 + p7X_M]);
              cnt[q*3 + p7X_I] = _mm256_add_ps(pp[q*3 + p7X_I], cnt[q*3 + p7X_I]);
            }
          cntx[p7X_N] += ppx[p7X_N];
          cntx[p7X_C] += ppx[p7X_C];
          cntx[p7X_J] += ppx[p7X_J];
        }
    }
  if (isinf(scaleproduct)) { status = eslERANGE; goto ERROR; }

  if (opt_null2 != NULL) p7_Null2_FromCounts_avx(om, cnt, cntx, L, opt_null2);
  if (opt_fsc   != NULL) *opt_fsc = totscale + log(wk->xmx[L*p7X_NXCELLS+p7X_C] * om->xf[p7O_C][p7O_MOVE]);
  *ret_oasc = oxmx[L*p7X_NXCELLS+p7X_C];

  p7_omx_Destroy_avx(wk);
  free(oxmx);
  return eslOK;

 ERROR:
  if (wk) p7_omx_Destroy_avx(wk);
  free(oxmx);
  return status;
}


/* Function:  p7_OATraceSweep_avx()
 * Synopsis:  OA trace from what the AVX OA sweep stored.
 *
 * Purpose:   The OA trace for profile <om> from <ox>, filled by
 *            <p7_OptimalAccuracySweep_avx()>, into the empty trace <tr>
 *            (created with <p7_trace_CreateWithPP()>). The choices are
 *            those <p7_OATrace_avx()> makes on the same OA scores. The
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
p7_OATraceSweep_avx(const P7_OPROFILE *om, const P7_OATMX *ox, P7_TRACE *tr)
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


#endif /* eslENABLE_AVX */
