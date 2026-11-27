/* Frameshift optimal accuracy alignment in one sweep over a
 * checkpointed Backward matrix, and its trace; AVX2 implementation.
 *
 * Contents:
 *    1. p7_OptimalAccuracySweep_Frameshift_avx(), p7_OATraceSweep_Frameshift_avx()
 *
 * Forward rows are computed by p7_ForwardRow_Frameshift_avx(), posteriors
 * by p7_DecodingRow_Frameshift_avx() and OA rows by
 * p7_OptimalAccuracyRow_Frameshift_avx() (optacc_fs_avx.c).
 */
#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <float.h>
#include <string.h>

#include <immintrin.h>  /* AVX2 */

#include "easel.h"
#include "esl_avx.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"


/* avx_rightshift_ps: right-shift a __m256 by one float lane, filling lane 0 with infv. */
static inline __m256
avx_rightshift_ps(__m256 v, __m256 infv)
{
  return _mm256_blend_ps(esl_avx_rightshiftz_float(v), infv, 0x01);
}


/*****************************************************************
 * 1. OA sweep over a checkpointed Backward, and its trace
 *****************************************************************/

/* oat_ppcodes8_fs()
 * Posterior codes (p7OAT_PPCODE(), 0..10) of the 8 lanes of <v> as 4
 * lane masks, bit b of each lane's code in <m[b]>.
 */
static inline void
oat_ppcodes8_fs(__m256 v, int m[4])
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

/* oat_store_row_fs_avx()
 * Store into <ox> what the frameshift OA trace reads from OA row <i>
 * (<dpc>, in the OA cell layout above, with the OA special states in
 * <xmx>): for each k, the predecessors of an M at (i+c,k), an I at
 * (i+3,k) and a D at (i,k), chosen with the same comparisons, shifts
 * and ties as p7_OATrace_Frameshift_avx() makes on the OA matrix; the
 * codon length of an M at (i,k), from the fill; the posterior codes of
 * M (pp_C0) and I at (i,k) from the posterior row <ppp> (all 0 for row
 * 0, <ppp> NULL); the row's OA special states and the posteriors of N,
 * J, C (<ppx>, NULL for row 0); and the cell an E at row i comes from,
 * in p7_OATrace_Frameshift_avx()'s scan order.
 */
static void
oat_store_row_fs_avx(const P7_FS_OPROFILE *om_fs, P7_OATMX *ox, const __m256 *dpc, const float *xmx, const __m256 *ppp, const float *ppx, int i)
{
  int      M     = om_fs->M;
  int      Q     = p7O_NQF_AVX(M);
  __m256   zerov = _mm256_setzero_ps();
  __m256   infv  = _mm256_set1_ps(-eslINFINITY);
  __m256   xBv   = _mm256_set1_ps(xmx[i*p7X_NXCELLS+p7X_B]);
  __m256   mpv, ipv, dpv, tmdv, tddv;
  __m256   pM, pI, pD, pB, best, gt, b0, b1;
  __m256   cv;
  __m256  *tp;
  union { __m256 v; float p[8]; } u;
  float    max;
  int      kmax, smax;
  int      mm[4], mi[4], mc;
  int      q, r, c, k;

  for (q = 0; q < Q; q++)
    {
      tp = om_fs->tfv_avx + 7*q;

      /* row i at k-1 */
      if (q > 0) {
        mpv  = MMO_OA(dpc, q-1);
        ipv  = IMO_OA(dpc, q-1);
        dpv  = DMO_OA(dpc, q-1);
        tmdv = om_fs->tfv_avx[7*(q-1) + p7O_MD];
        tddv = om_fs->tfv_avx[7*Q + (q-1)];
      } else {
        mpv  = avx_rightshift_ps(MMO_OA(dpc, Q-1), infv);
        ipv  = avx_rightshift_ps(IMO_OA(dpc, Q-1), infv);
        dpv  = avx_rightshift_ps(DMO_OA(dpc, Q-1), infv);
        tmdv = esl_avx_rightshiftz_float(om_fs->tfv_avx[7*(Q-1) + p7O_MD]);
        tddv = esl_avx_rightshiftz_float(om_fs->tfv_avx[8*Q - 1]);
      }

      /* M at (i+c,k): first best of M, I, D, B */
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

      /* I at (i+3,k): M on ties */
      pM = _mm256_blendv_ps(MMO_OA(dpc, q), infv, _mm256_cmp_ps(tp[p7O_MI], zerov, _CMP_EQ_OQ));
      pI = _mm256_blendv_ps(IMO_OA(dpc, q), infv, _mm256_cmp_ps(tp[p7O_II], zerov, _CMP_EQ_OQ));
      OATO_STRIPE8(ox, i, p7OAT_IP, q) = (uint8_t) _mm256_movemask_ps(_mm256_cmp_ps(pI, pM, _CMP_GT_OQ));

      /* D at (i,k): M on ties */
      pM = _mm256_blendv_ps(mpv, infv, _mm256_cmp_ps(tmdv, zerov, _CMP_EQ_OQ));
      pD = _mm256_blendv_ps(dpv, infv, _mm256_cmp_ps(tddv, zerov, _CMP_EQ_OQ));
      OATO_STRIPE8(ox, i, p7OAT_DP, q) = (uint8_t) _mm256_movemask_ps(_mm256_cmp_ps(pD, pM, _CMP_GT_OQ));

      /* codon length of an M here (1..5), as the fill stored it */
      cv = CODONO_OA(dpc, q);
      for (c = 0; c < 3; c++)
        {
          u.v = cv; mc = 0;
          for (r = 0; r < 8; r++) mc |= ((((int) u.p[r]) >> c) & 1) << r;
          OATO_STRIPE8(ox, i, p7OAT_COD + c, q) = (uint8_t) mc;
        }

      /* posterior codes of M and I */
      if (ppp != NULL) { oat_ppcodes8_fs(MMO_FS(ppp, q, p7X_FS_C0), mm); oat_ppcodes8_fs(IMO_FS(ppp, q), mi); }
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

  /* the cell an E at row i comes from, as select_e_fs() scans the row */
  max  = -eslINFINITY;
  smax = p7T_M;
  kmax = 1;
  for (k = 1; k <= M; k++)
    {
      q = (k-1) % Q;
      r = (k-1) / Q;
      u.v = MMO_OA(dpc, q); if (u.p[r] > max) { max = u.p[r]; smax = p7T_M; kmax = k; }
      u.v = DMO_OA(dpc, q); if (u.p[r] > max) { max = u.p[r]; smax = p7T_D; kmax = k; }
    }
  ox->ek[i] = (smax == p7T_M) ? kmax : -kmax;
}


/* Function:  p7_OptimalAccuracySweep_Frameshift_avx()
 * Synopsis:  Frameshift OA alignment in one forward sweep over a checkpointed Backward, AVX.
 *
 * Purpose:   For the 5-codon frameshift profile <om_fs> and DNA sequence
 *            <dsq> of length <L>, with the checkpointed Backward matrix
 *            <cx> filled by <p7_BackwardChk_Frameshift_avx()>, compute in
 *            one sweep through rows 1..L the Forward, posterior and
 *            optimal accuracy (OA) rows, computing each block of
 *            Backward rows again from its checkpoint as the sweep reaches
 *            it, and store in <ox> what <p7_OATraceSweep_Frameshift_avx()>
 *            needs for the OA trace. <ov> is workspace for 5 codon
 *            lengths. Rows of Forward, posterior and OA values are kept
 *            only as long as later rows need them (up to 5 rows back).
 *
 *            Each row is calculated as <p7_Forward_Frameshift_avx()>,
 *            <p7_Decoding_Frameshift_avx()> and
 *            <p7_OptimalAccuracy_Frameshift_avx()> calculate it; the
 *            posteriors are those <p7_Decoding_Frameshift_avx()> gives for
 *            the checkpointed Backward, which has its own scale factors.
 *
 *            If <opt_null2> is non-NULL, the null2 odds ratios of the
 *            sequence are also calculated into it, as
 *            <p7_Null2_fs_ByExpectation_avx()> calculates them from the
 *            posteriors: the posterior rows are summed as the sweep
 *            computes them, in the same order.
 *
 *            <ox> must have been created with <p7OAT_NPLANES_FS> planes;
 *            its layout is set to 8 lanes and it is grown for <om_fs->M>
 *            by <L>. On return <cx> has its last block loaded.
 *
 * Args:      dsq       - DNA sequence in digitized form, 1..L
 *            L         - length of dsq
 *            om_fs     - 5-codon frameshift profile
 *            cx        - checkpointed Backward matrix from p7_BackwardChk_Frameshift_avx()
 *            ov        - workspace for the profile's M and 5 codon lengths
 *            ox        - RESULT: what the OA trace reads
 *            opt_null2 - optRETURN: null2 odds ratios, <0..Kp-1>; caller allocated
 *            opt_fsc   - optRETURN: Forward score in nats
 *            ret_oasc  - RETURN: OA score, the expected number of
 *                        correctly aligned residues
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; <eslEINVAL> if the arguments
 *            aren't set up for this computation; <eslERANGE> if the
 *            Forward score or a posterior scale factor is out of range (as
 *            <p7_Forward_Frameshift_avx()> and <p7_Decoding_Frameshift_avx()>).
 */
int
p7_OptimalAccuracySweep_Frameshift_avx(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, P7_OATMX *ox,
                                       float *opt_null2, float *opt_fsc, float *ret_oasc)
{
  int            M        = om_fs->M;
  int            Q        = p7O_NQF_AVX(M);
  P7_OMX        *wk       = NULL;   /* rows 0..3: Forward ring; 4..9: OA ring; 10: posteriors; 11: posterior sums; xmx: Forward specials */
  P7_OIVX       *bv       = NULL;   /* workspace for the Backward block reloads, apart from Forward's codon sums in <ov> */
  float         *oxmx     = NULL;
  float         *log_sfwd = NULL;
  float         *log_sbck = NULL;
  float          log_inv_Z;
  float          ppx[p7X_NXCELLS];
  float          cntx[p7X_NXCELLS];
  float          totscale = 0.0f;
  float          xCtot;
  P7_FS_FWDRING  rg;
  __m256        *pp, *cnt;
  __m256        *oa[6];
  int            i, q, b;
  int            status;

  if (om_fs->codon_lengths != 5)        ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");
  if (ox->nplanes != p7OAT_NPLANES_FS)  ESL_EXCEPTION(eslEINVAL, "P7_OATMX not set up for the frameshift OA sweep");
  if (cx->L != L || cx->M != M || cx->dir != p7_CHK_BACKWARD || cx->nscells != p7X_NSCELLS) ESL_EXCEPTION(eslEINVAL, "checkpointed Backward matrix doesn't match");
  ox->lanes = 8;
  if ((status = p7_oatmx_GrowTo(ox, M, L)) != eslOK) return status;

  if ((wk = p7_omx_Create_dpf_avx(M, 11, L, p7X_NSCELLS_FS)) == NULL) { status = eslEMEM; goto ERROR; }
  if ((bv = p7_oivx_Create_avx(M, p7P_5CODONS))              == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(oxmx,     sizeof(float) * (L+1) * p7X_NXCELLS);
  ESL_ALLOC(log_sfwd, sizeof(float) * (L+2));
  ESL_ALLOC(log_sbck, sizeof(float) * (L+2));
  for (q = 0; q < 6; q++) oa[q] = wk->dpf_avx[4+q];
  pp  = wk->dpf_avx[10];
  cnt = wk->dpf_avx[11];
  for (q = 0; q < p7X_NXCELLS; q++) cntx[q] = 0.0f;

  /* Backward's cumulative log scale products and normalization, as p7_Decoding_Frameshift_avx() computes them */
  log_sbck[L+1] = 0.0f;
  for (i = L; i >= 0; i--)
    log_sbck[i] = log_sbck[i+1] + logf(cx->xmx[i*p7X_NXCELLS + p7X_SCALE]);
  log_inv_Z = -p7_FLogsum(
                 logf(cx->xmx[p7X_N]) + log_sbck[0],
                 p7_FLogsum(
                   logf(cx->xmx[1*p7X_NXCELLS + p7X_N]) + log_sbck[1],
                   logf(cx->xmx[2*p7X_NXCELLS + p7X_N]) + log_sbck[2]));

  p7_ForwardRow_Frameshift_avx(dsq, om_fs, wk->xmx, &totscale, &rg, ov->ivx_avx, wk->dpf_avx[0], NULL, NULL, 0);
  log_sfwd[0] = logf(wk->xmx[p7X_SCALE]);
  p7_OptimalAccuracyRow_Frameshift_avx(om_fs, NULL, NULL, oxmx, oa[0], NULL, NULL, NULL, NULL, NULL, 0);
  oat_store_row_fs_avx(om_fs, ox, oa[0], oxmx, NULL, NULL, 0);

  for (i = 1; i <= L; i++)
    {
      b = (i / cx->blockR) * cx->blockR;
      if (b != cx->b && (status = p7_BackwardChkLoad_Frameshift_avx(dsq, om_fs, cx, bv, b)) != eslOK) goto ERROR;

      p7_ForwardRow_Frameshift_avx(dsq, om_fs, wk->xmx, &totscale, &rg, ov->ivx_avx, wk->dpf_avx[i%4], wk->dpf_avx[(i-1)%4],
                                   (i >= 3) ? wk->dpf_avx[(i-3)%4] : NULL, i);
      log_sfwd[i] = log_sfwd[i-1] + logf(wk->xmx[i*p7X_NXCELLS + p7X_SCALE]);

      if ((status = p7_DecodingRow_Frameshift_avx(om_fs, cx->xmx, (i > 2) ? wk->xmx + (i-3)*p7X_NXCELLS : NULL, log_sfwd, log_sbck, log_inv_Z,
                                                  wk->dpf_avx[i%4], cx->dpf_avx[i], pp, ppx, i)) != eslOK) goto ERROR;

      /* row 0 (all -inf M, I, D) stands in for rows before it */
      p7_OptimalAccuracyRow_Frameshift_avx(om_fs, pp, ppx, oxmx, oa[i%6], oa[(i-1)%6],
                    (i >= 2) ? oa[(i-2)%6] : oa[0],
                    (i >= 3) ? oa[(i-3)%6] : oa[0],
                    (i >= 4) ? oa[(i-4)%6] : oa[0],
                    (i >= 5) ? oa[(i-5)%6] : oa[0], i);
      oat_store_row_fs_avx(om_fs, ox, oa[i%6], oxmx, pp, ppx, i);

      /* posterior sums for null2, added up as p7_Null2_fs_ByExpectation_avx() does */
      if (i == 1)
        {
          memcpy(cnt, pp, sizeof(__m256) * p7X_NSCELLS_FS * Q);
          cntx[p7X_N] = ppx[p7X_N];
          cntx[p7X_C] = ppx[p7X_C];
          cntx[p7X_J] = ppx[p7X_J];
        }
      else
        {
          for (q = 0; q < Q; q++)
            {
              MMO_FS(cnt, q, p7X_FS_C0) = _mm256_add_ps(MMO_FS(pp, q, p7X_FS_C0), MMO_FS(cnt, q, p7X_FS_C0));
              IMO_FS(cnt, q)             = _mm256_add_ps(IMO_FS(pp, q),             IMO_FS(cnt, q));
            }
          cntx[p7X_N] += ppx[p7X_N];
          cntx[p7X_C] += ppx[p7X_C];
          cntx[p7X_J] += ppx[p7X_J];
        }
    }

  /* Forward score, as p7_Forward_Frameshift_avx() computes it */
  xCtot = rg.xC[   L    %PARSER_ROWS_FWD]
        + rg.xC[((L-1)%PARSER_ROWS_FWD+PARSER_ROWS_FWD)%PARSER_ROWS_FWD]*om_fs->xf[p7O_C][p7O_LOOP]
        + rg.xC[((L-2)%PARSER_ROWS_FWD+PARSER_ROWS_FWD)%PARSER_ROWS_FWD]*om_fs->xf[p7O_C][p7O_LOOP];
  if (isnan(xCtot) || isinf(xCtot) == 1 || (L > 1 && xCtot == 0.0f)) { status = eslERANGE; goto ERROR; }

  if (opt_null2 != NULL) p7_Null2_fs_FromCounts_avx(om_fs, cnt, cntx, opt_null2);
  if (opt_fsc   != NULL) *opt_fsc = totscale + logf(xCtot*om_fs->xf[p7O_C][p7O_MOVE]);
  *ret_oasc = oxmx[L*p7X_NXCELLS+p7X_C] + oxmx[(L-1)*p7X_NXCELLS+p7X_C] + oxmx[(L-2)*p7X_NXCELLS+p7X_C];

  p7_omx_Destroy_avx(wk);
  p7_oivx_Destroy_avx(bv);
  free(oxmx); free(log_sfwd); free(log_sbck);
  return eslOK;

 ERROR:
  if (wk) p7_omx_Destroy_avx(wk);
  if (bv) p7_oivx_Destroy_avx(bv);
  free(oxmx); free(log_sfwd); free(log_sbck);
  return status;
}


/* Function:  p7_OATraceSweep_Frameshift_avx()
 * Synopsis:  Frameshift OA trace from what the AVX OA sweep stored.
 *
 * Purpose:   The frameshift OA trace for profile <om_fs> from <ox>, filled
 *            by <p7_OptimalAccuracySweep_Frameshift_avx()>, into the empty
 *            trace <tr> (created with <p7_trace_fs_CreateWithPP()>). The
 *            choices are those <p7_OATrace_Frameshift_avx()> makes on the
 *            same OA scores. The posterior probability of each residue is
 *            the value its PP character stands for (<p7OAT_PPVALUE()>),
 *            from the posterior row, as <p7_alidisplay_EncodePostProb()>
 *            prints it.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if <tr> isn't
 *            empty or a bogus state is reached.
 */
int
p7_OATraceSweep_Frameshift_avx(const P7_FS_OPROFILE *om_fs, const P7_OATMX *ox, P7_TRACE *tr)
{
  static const int state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };   /* indexed by p7OAT_FROM_* */
  int   L = ox->L;
  int   i = L;
  int   k = 0;
  int   c = 0;
  float t1, t2;
  float path[4];
  float postprob;
  int   sprv, scur;
  int   status;

  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace not empty; needs to be Reuse()'d?");

  if ((status = p7_trace_fs_AppendWithPP(tr, p7T_T, k, i, c, 0.0f)) != eslOK) return status;
  if ((status = p7_trace_fs_AppendWithPP(tr, p7T_C, k, i, c, 0.0f)) != eslOK) return status;

  sprv = p7T_C;
  while (sprv != p7T_S)
    {
      switch (sprv) {
      case p7T_M: scur = state[OATO_MPRED(ox, i, k)];                                                   k--;    break;
      case p7T_D: scur = OATO_DPRED(ox, i, k) ? p7T_D : p7T_M;                                           k--;    break;
      case p7T_I: scur = (OATO_IPRED(ox, (i >= 3) ? i-3 : 0, k) == p7OAT_FROM_I) ? p7T_I : p7T_M;        i -= 3; break;
      case p7T_N: scur = (i == 0) ? p7T_S : p7T_N;                                                                break;
      case p7T_C:
        t1 = (om_fs->xf[p7O_C][p7O_LOOP] == 0.0f) ? 0.0f : 1.0f;
        t2 = (om_fs->xf[p7O_E][p7O_MOVE] == 0.0f) ? 0.0f : 1.0f;
        if (i < 4) { scur = p7T_E; break; }
        path[0] = (t1 == 0.0f) ? -eslINFINITY :
                  ox->xmx[(i-3)*p7OAT_NXCELLS+p7OAT_C] + ox->xmx[ i   *p7OAT_NXCELLS+p7OAT_PPC];
        path[1] = (i < L && t1 != 0.0f) ?
                  ox->xmx[(i-2)*p7OAT_NXCELLS+p7OAT_C] + ox->xmx[(i+1)*p7OAT_NXCELLS+p7OAT_PPC] : -eslINFINITY;
        path[2] = (i < L-1 && t1 != 0.0f) ?
                  ox->xmx[(i-1)*p7OAT_NXCELLS+p7OAT_C] + ox->xmx[(i+2)*p7OAT_NXCELLS+p7OAT_PPC] : -eslINFINITY;
        path[3] = (t2 == 0.0f) ? -eslINFINITY : ox->xmx[i*p7OAT_NXCELLS+p7OAT_E];
        scur = (esl_vec_FArgMax(path, 4) == 3) ? p7T_E : p7T_C;
        break;
      case p7T_J:
        if (i <= 5) { scur = p7T_E; break; }
        path[0] = (om_fs->xf[p7O_J][p7O_LOOP] == 0.0f) ? -eslINFINITY :
                  ox->xmx[i*p7OAT_NXCELLS+p7OAT_J] + ox->xmx[i*p7OAT_NXCELLS+p7OAT_PPJ];
        path[1] = (om_fs->xf[p7O_E][p7O_LOOP] == 0.0f) ? -eslINFINITY :
                  ox->xmx[i*p7OAT_NXCELLS+p7OAT_E];
        scur = (esl_vec_FArgMax(path, 2) == 0) ? p7T_J : p7T_E;
        break;
      case p7T_E:
        scur = (ox->ek[i] > 0) ? p7T_M : p7T_D;
        k    = abs(ox->ek[i]);
        break;
      case p7T_B:
        path[0] = (om_fs->xf[p7O_N][p7O_MOVE] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7OAT_NXCELLS+p7OAT_N];
        path[1] = (om_fs->xf[p7O_J][p7O_MOVE] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7OAT_NXCELLS+p7OAT_J];
        scur = (path[0] > path[1]) ? p7T_N : p7T_J;
        break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in OA FS traceback");
      }

      switch (scur) {
      case p7T_M: postprob = p7OAT_PPVALUE(OATO_PPM(ox, i, k));                                break;
      case p7T_I: postprob = p7OAT_PPVALUE(OATO_PPI(ox, i, k));                                break;
      case p7T_N: postprob = (sprv == scur) ? ox->xmx[i*p7OAT_NXCELLS+p7OAT_PPN] : 0.0f;       break;
      case p7T_C: postprob = (sprv == scur) ? ox->xmx[i*p7OAT_NXCELLS+p7OAT_PPC] : 0.0f;       break;
      case p7T_J: postprob = (sprv == scur) ? ox->xmx[i*p7OAT_NXCELLS+p7OAT_PPJ] : 0.0f;       break;
      default:    postprob = 0.0f;                                                             break;
      }
      c = (scur == p7T_M) ? OATO_CODON(ox, i, k) : 0;

      if ((status = p7_trace_fs_AppendWithPP(tr, scur, k, i, c, postprob)) != eslOK) return status;

      if ( (scur == p7T_N || scur == p7T_C || scur == p7T_J) && scur == sprv) i--;
      sprv = scur;
      i   -= c;
    }

  tr->M = om_fs->M;
  tr->L = L;
  return p7_trace_fs_Reverse(tr);
}

#endif /* eslENABLE_AVX */
