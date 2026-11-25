/* Frameshift optimal accuracy alignment in a forward sweep, generic
 * (non-SIMD) version.
 *
 * The Backward matrix is first computed checkpointed
 * (<p7_GBackwardChk_Frameshift()>). One sweep through the rows 1..L
 * then computes, for each row, the Forward row, the Backward row (its
 * block computed again from its checkpoint), the posterior
 * probabilities and the optimal accuracy (OA) scores, keeping only the
 * rows later rows read (up to 5 back). What the OA trace reads is
 * stored in a <P7_OATMX> as each row is finished, and
 * <p7_GOATraceSweep_Frameshift()> traces from it alone, so no
 * full-size matrix is kept. Traces are identical to those of
 * <p7_GForward_Frameshift()>, <p7_GBackward_Frameshift()>,
 * <p7_GDecoding_Frameshift()>, <p7_GOptimalAccuracy_Frameshift()> and
 * <p7_GOATrace_Frameshift()>, with each residue's posterior probability
 * given as the value its PP character stands for. The standard-model
 * version is in generic_optacc_sweep.c.
 *
 * Computing posterior probabilities from Backward rows as they are
 * finished, without a full Backward matrix, follows HMMER4's Backward
 * filter; a checkpointed Backward matrix that a forward sweep reads
 * row by row, computing Forward, posteriors and OA together, is an
 * elaboration of that idea.
 *
 * Contents:
 *   1. OA sweep and trace, frameshift model (5 codon lengths)
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
 * 1. OA sweep and trace, frameshift model (5 codon lengths)
 *****************************************************************/

#define IVX5(i,k)     (ivx[((k)*p7P_5CODONS) + (i)])
#define MSC_FS1(k)    (rsc1[(k)])
#define MSC_FS2(k)    (rsc2[(k)])
#define MSC_FS3(k)    (rsc3[(k)])
#define MSC_FS4(k)    (rsc4[(k)])
#define MSC_FS5(k)    (rsc5[(k)])
#define MMX_OA(i,k)   (dp[(i)][(k) * p7G_NSCELLS_FS + p7G_M + p7G_C1])
#define IMX_OA(i,k)   (dp[(i)][(k) * p7G_NSCELLS_FS + p7G_M + p7G_C2])
#define DMX_OA(i,k)   (dp[(i)][(k) * p7G_NSCELLS_FS + p7G_D])
#define PPM_FS(k,c)   (pp[(k) * p7G_NSCELLS_FS + p7G_M + (c)])
#define PPI_FS(k)     (pp[(k) * p7G_NSCELLS_FS + p7G_I])

/* decode_fs_row()
 * Posterior probabilities of row <i> (1..L) of the frameshift model
 * into <pp> (one row, [k*p7G_NSCELLS_FS + s]) and <ppx> (N, J, C of
 * the row), from Forward row <i> in <fdp>, the Forward special states
 * of row i-3 in <fx3> (used for i > 2), the Backward row <i> in <bdp>
 * and its special states in <bxp>, with overall score <overall_sc>.
 * Same calculation as p7_GDecoding_Frameshift().
 */
static void
decode_fs_row(const P7_FS_PROFILE *gm_fs5, const float *fdp, const float *fx3, const float *bdp, const float *bxp,
              float overall_sc, float *pp, float *ppx, int i)
{
  int M = gm_fs5->M;
  int k, c;

  for (c = p7G_C0; c <= p7G_C5; c++) PPM_FS(0,c) = 0.0;
  PPI_FS(0) = pp[0*p7G_NSCELLS_FS + p7G_D] = 0.0;

  for (k = 1; k < M; k++)
    {
      for (c = p7G_C0; c <= p7G_C5; c++)
        PPM_FS(k,c) = expf(fdp[k*p7G_NSCELLS_FS + p7G_M + c] + bdp[k*p7G_NSCELLS + p7G_M] - overall_sc);
      PPI_FS(k) = expf(fdp[k*p7G_NSCELLS_FS + p7G_I] + bdp[k*p7G_NSCELLS + p7G_I] - overall_sc);
      pp[k*p7G_NSCELLS_FS + p7G_D] = 0.0;
    }
  for (c = p7G_C0; c <= p7G_C5; c++)
    PPM_FS(M,c) = expf(fdp[M*p7G_NSCELLS_FS + p7G_M + c] + bdp[M*p7G_NSCELLS + p7G_M] - overall_sc);
  PPI_FS(M) = 0.0;
  pp[M*p7G_NSCELLS_FS + p7G_D] = 0.0;

  ppx[p7G_E] = 0.0;
  ppx[p7G_B] = 0.0;
  if (i > 2)
    {
      ppx[p7G_N] = expf(fx3[p7G_N] + bxp[p7G_N] + gm_fs5->xsc[p7P_N][p7P_LOOP] - overall_sc);
      ppx[p7G_C] = expf(fx3[p7G_C] + bxp[p7G_C] + gm_fs5->xsc[p7P_C][p7P_LOOP] - overall_sc);
      ppx[p7G_J] = expf(fx3[p7G_J] + bxp[p7G_J] + gm_fs5->xsc[p7P_J][p7P_LOOP] - overall_sc);
    }
  else
    {
      ppx[p7G_N] = expf(bxp[p7G_N] - overall_sc);
      ppx[p7G_C] = 0.0;
      ppx[p7G_J] = 0.0;
    }
}

/* oa_fs_max()
 * For oa_fs_row(): the best OA path into an M at (i,k) through codon
 * length <c>, from row i-c (cells at k-1 and B), as
 * p7_GOptimalAccuracy_Frameshift() computes its max1..max5.
 */
static inline float
oa_fs_max(const float *tsc, float **dp, const float *xmx, const float *pp, int i, int k, int c)
{
  return ESL_MAX( TSCDELTA(p7P_MM, k-1) * (MMX_OA(i-c,k-1)   + PPM_FS(k, p7G_C0 + c)),
         ESL_MAX( TSCDELTA(p7P_IM, k-1) * (IMX_OA(i-c,k-1)   + PPM_FS(k, p7G_C0 + c)),
         ESL_MAX( TSCDELTA(p7P_DM, k-1) * (DMX_OA(i-c,k-1)   + PPM_FS(k, p7G_C0 + c)),
                  TSCDELTA(p7P_BM, k-1) * (XMX(i-c,p7G_B)    + PPM_FS(k, p7G_C0 + c)))));
}

/* oa_fs_row()
 * OA row <i> (1..L) of the frameshift model into dp[i] (OA cells) and
 * xmx row i, from OA rows i-1..i-5 and the posterior row <pp> / <ppx>
 * of row i. Same calculation as p7_GOptimalAccuracy_Frameshift().
 */
static void
oa_fs_row(const P7_FS_PROFILE *gm_fs5, const float *pp, const float *ppx, float **dp, float *xmx, int i)
{
  float const *tsc  = gm_fs5->tsc;
  int          M    = gm_fs5->M;
  float        esc  = p7_fs_profile_IsLocal(gm_fs5) ? 1.0 : 0.0;
  float        nn   = ((gm_fs5->xsc[p7P_N][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  float        jj   = ((gm_fs5->xsc[p7P_J][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  float        cc   = ((gm_fs5->xsc[p7P_C][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  float        nb   = ((gm_fs5->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  float        jb   = ((gm_fs5->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  float        ej   = ((gm_fs5->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
  float        ec   = ((gm_fs5->xsc[p7P_E][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
  float        max1, max2, max3, max4, max5;
  int          k;

  MMX_OA(i,0) = IMX_OA(i,0) = DMX_OA(i,0) = XMX(i,p7G_E) = -eslINFINITY;

  if (i == 1)
    {
      for (k = 1; k < M; k++) {
        MMX_OA(1,k)  = TSCDELTA(p7P_BM,k-1) * PPM_FS(k,p7G_C1);
        IMX_OA(1,k)  = -eslINFINITY;
        DMX_OA(1,k)  = ESL_MAX( TSCDELTA(p7P_MD, k-1) * MMX_OA(1,k-1),
                                TSCDELTA(p7P_DD, k-1) * DMX_OA(1,k-1));
        XMX(1,p7G_E) = ESL_MAX(XMX(1,p7G_E), esc * MMX_OA(1,k));
      }
      MMX_OA(1,M) = TSCDELTA(p7P_BM, M-1) * PPM_FS(M,p7G_C1);
      IMX_OA(1,M)  = -eslINFINITY;
      DMX_OA(1,M)  = ESL_MAX( TSCDELTA(p7P_MD, M-1) * MMX_OA(1,M-1),
                              TSCDELTA(p7P_DD, M-1) * DMX_OA(1,M-1));
      XMX(1,p7G_E) = ESL_MAX(XMX(1,p7G_E), ESL_MAX(MMX_OA(1,M), DMX_OA(1, M)));

      XMX(1,p7G_J) = ej * XMX(1,p7G_E);
      XMX(1,p7G_C) = ec * XMX(1,p7G_E);
      XMX(1,p7G_N) = nn * ppx[p7G_N];
      XMX(1,p7G_B) = ESL_MAX(nb * XMX(1,p7G_N), jb * XMX(1,p7G_J));
      return;
    }

  if (i == 2)
    {
      for (k = 1; k <= M; k++) {
        max1 = oa_fs_max(tsc, dp, xmx, pp, 2, k, 1);
        max2 = TSCDELTA(p7P_BM, k-1) * (XMX(0,p7G_B) + PPM_FS(k,p7G_C2));
        MMX_OA(2,k)  = ESL_MAX( max1, max2);
        IMX_OA(2,k)  = -eslINFINITY;
        DMX_OA(2,k)  = ESL_MAX( TSCDELTA(p7P_MD, k-1) * MMX_OA(2,k-1),
                                TSCDELTA(p7P_DD, k-1) * DMX_OA(2,k-1));
        if (k < M) XMX(2,p7G_E) = ESL_MAX(XMX(2,p7G_E), esc * MMX_OA(2,k));
      }
      XMX(2,p7G_E) = ESL_MAX(XMX(2,p7G_E), ESL_MAX(MMX_OA(2,M), DMX_OA(2, M)));

      XMX(2,p7G_J) = ej * XMX(2,p7G_E);
      XMX(2,p7G_C) = ec * XMX(2,p7G_E);
      XMX(2,p7G_N) = nn * ppx[p7G_N];
      XMX(2,p7G_B) = ESL_MAX(nb * XMX(2,p7G_N), jb * XMX(2,p7G_J));
      return;
    }

  for (k = 1; k <= M; k++)
    {
      max1 = oa_fs_max(tsc, dp, xmx, pp, i, k, 1);
      max2 = oa_fs_max(tsc, dp, xmx, pp, i, k, 2);
      max3 = oa_fs_max(tsc, dp, xmx, pp, i, k, 3);
      if (i >= 5)
        {
          max4 = oa_fs_max(tsc, dp, xmx, pp, i, k, 4);
          max5 = oa_fs_max(tsc, dp, xmx, pp, i, k, 5);
          MMX_OA(i,k) = ESL_MAX(ESL_MAX(max1, max2),
                                ESL_MAX(ESL_MAX(max3, max4), max5));
        }
      else
        {
          max4 = (i == 4) ? TSCDELTA(p7P_BM, k-1) * (XMX(0,p7G_B) + PPM_FS(k,p7G_C4)) : -eslINFINITY;
          MMX_OA(i,k) = ESL_MAX( ESL_MAX( max1, max2),
                                 ESL_MAX( max3, max4));
        }

      if (k < M)
        IMX_OA(i,k) = ESL_MAX( TSCDELTA(p7P_MI, k) * (MMX_OA(i-3,k) + PPI_FS(k)),
                               TSCDELTA(p7P_II, k) * (IMX_OA(i-3,k) + PPI_FS(k)));
      else
        IMX_OA(i,k) = -eslINFINITY;

      DMX_OA(i,k) = ESL_MAX( TSCDELTA(p7P_MD, k-1) * MMX_OA(i,k-1),
                             TSCDELTA(p7P_DD, k-1) * DMX_OA(i,k-1));

      if (k < M) XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), esc * MMX_OA(i,k));
    }
  XMX(i,p7G_E) = ESL_MAX(XMX(i,p7G_E), ESL_MAX(MMX_OA(i,M), DMX_OA(i, M)));

  XMX(i,p7G_J) = ESL_MAX( jj * (XMX(i-3,p7G_J) + ppx[p7G_J]),
                          ej *  XMX(i,  p7G_E));
  XMX(i,p7G_C) = ESL_MAX( cc * (XMX(i-3,p7G_C) + ppx[p7G_C]),
                          ec *  XMX(i,  p7G_E));
  XMX(i,p7G_N) = nn * (XMX(i-3,p7G_N) + ppx[p7G_N]);
  XMX(i,p7G_B) = ESL_MAX( nb * XMX(i,  p7G_N), jb * XMX(i,  p7G_J));
}

/* oat_store_row_fs()
 * Store into <ox> what the frameshift OA trace reads from OA row <i>:
 * for each k, the predecessors of an M at (i+c,k), an I at (i+3,k) and
 * a D at (i,k), chosen as p7_GOATrace_Frameshift() chooses them from
 * the OA matrix; the codon length an M at (i,k) takes, the argmax of
 * pp_C1..C5 (the first on ties); the posterior codes of M (pp_C0) and
 * I at (i,k) from the posterior row <pp> (all 0 for row 0); the row's
 * OA special states and posteriors of N, J, C (<ppx>, NULL for row 0);
 * and the cell an E at row i comes from.
 */
static void
oat_store_row_fs(const P7_FS_PROFILE *gm_fs5, float **dp, const float *xmx, const float *pp, const float *ppx, P7_OATMX *ox, int i)
{
  static const int from[4] = { p7OAT_FROM_M, p7OAT_FROM_I, p7OAT_FROM_D, p7OAT_FROM_B };
  float const *tsc  = gm_fs5->tsc;
  int          M    = gm_fs5->M;
  float        path[5];
  float        max;
  int          smax, kmax;
  int          k, c;

  for (k = 1; k <= M; k++)
    {
      path[0] = TSCDELTA(p7P_MM, k-1) * MMX_OA(i,k-1);
      path[1] = TSCDELTA(p7P_IM, k-1) * IMX_OA(i,k-1);
      path[2] = TSCDELTA(p7P_DM, k-1) * DMX_OA(i,k-1);
      path[3] = TSCDELTA(p7P_BM, k-1) * XMX(i,p7G_B);
      OATX_SET_MPRED(ox, i, k, from[esl_vec_FArgMax(path, 4)]);

      if (k < M)
        {
          path[0] = TSCDELTA(p7P_MI, k) * MMX_OA(i,k);
          path[1] = TSCDELTA(p7P_II, k) * IMX_OA(i,k);
          OATX_SET_IPRED(ox, i, k, (path[0] >= path[1]) ? p7OAT_FROM_M : p7OAT_FROM_I);
        }
      else OATX_SET_IPRED(ox, i, k, p7OAT_FROM_M);

      path[0] = TSCDELTA(p7P_MD, k-1) * MMX_OA(i,k-1);
      path[1] = TSCDELTA(p7P_DD, k-1) * DMX_OA(i,k-1);
      OATX_SET_DPRED(ox, i, k, (path[0] >= path[1]) ? 0 : 1);

      if (pp)
        {
          for (c = 0; c < 5; c++) path[c] = PPM_FS(k, p7G_C1 + c);
          OATX_SET_CODON(ox, i, k, esl_vec_FArgMax(path, 5) + 1);
          OATX_SET_PPM  (ox, i, k, p7OAT_PPCODE(PPM_FS(k, p7G_C0)));
          OATX_SET_PPI  (ox, i, k, p7OAT_PPCODE(PPI_FS(k)));
        }
      else
        {
          OATX_SET_CODON(ox, i, k, 1);
          OATX_SET_PPM  (ox, i, k, 0);
          OATX_SET_PPI  (ox, i, k, 0);
        }
    }

  OATX_XMX(ox, i, p7OAT_N)   = XMX(i, p7G_N);
  OATX_XMX(ox, i, p7OAT_J)   = XMX(i, p7G_J);
  OATX_XMX(ox, i, p7OAT_C)   = XMX(i, p7G_C);
  OATX_XMX(ox, i, p7OAT_B)   = XMX(i, p7G_B);
  OATX_XMX(ox, i, p7OAT_E)   = XMX(i, p7G_E);
  OATX_XMX(ox, i, p7OAT_PPN) = ppx ? ppx[p7G_N] : 0.0;
  OATX_XMX(ox, i, p7OAT_PPJ) = ppx ? ppx[p7G_J] : 0.0;
  OATX_XMX(ox, i, p7OAT_PPC) = ppx ? ppx[p7G_C] : 0.0;

  /* the cell an E at row i comes from, as p7_GOATrace_Frameshift()'s select_e() picks it */
  if (! p7_fs_profile_IsLocal(gm_fs5))
    ox->ek[i] = (MMX_OA(i,M) >= DMX_OA(i,M)) ? M : -M;
  else
    {
      max  = -eslINFINITY;
      smax = -1;
      kmax = 0;
      for (k = 1; k <= M; k++)
        {
          if (MMX_OA(i,k) > max) { max = MMX_OA(i,k); smax = p7T_M; kmax = k; }
          if (DMX_OA(i,k) > max) { max = DMX_OA(i,k); smax = p7T_D; kmax = k; }
        }
      ox->ek[i] = (smax == -1) ? 0 : ((smax == p7T_M) ? kmax : -kmax);
    }
}


/* Function:  p7_GOptimalAccuracySweep_Frameshift()
 * Synopsis:  Frameshift OA alignment in one forward sweep over a checkpointed Backward.
 *
 * Purpose:   For the 5-codon frameshift profile <gm_fs5> and DNA
 *            sequence <dsq> of length <L>, with the checkpointed
 *            Backward matrix <cx> filled by <p7_GBackwardChk_Frameshift()>,
 *            compute in one sweep through rows 1..L the Forward,
 *            posterior and optimal accuracy (OA) rows, computing each
 *            block of Backward rows again from its checkpoint as the
 *            sweep reaches it, and store in <ox> what
 *            <p7_GOATraceSweep_Frameshift()> needs for the OA trace.
 *            <iv> is workspace for 5 codon lengths. Rows of Forward,
 *            posterior and OA values are kept only as long as later
 *            rows need them (up to 5 rows back).
 *
 *            If <opt_null2> is non-NULL, the null2 odds ratios of the
 *            sequence are also calculated into it, as
 *            <p7_GNull2_fs_ByExpectation()> calculates them from the
 *            posteriors: the posterior rows are summed as the sweep
 *            computes them, in the same order.
 *
 *            <ox> must have been created with <p7OAT_NPLANES_FS> planes
 *            and 1 lane (the generic layout); it is grown for
 *            <gm_fs5->M> by <L>. On return <cx> has its last block
 *            loaded.
 *
 * Args:      dsq      - DNA sequence in digitized form, 1..L
 *            L        - length of dsq (>= 3)
 *            gm_fs5   - 5-codon frameshift profile
 *            cx       - checkpointed Backward matrix from p7_GBackwardChk_Frameshift()
 *            iv       - workspace, allocated for 5 codon lengths
 *            ox       - RESULT: what the OA trace reads
 *            opt_null2- optRETURN: null2 odds ratios, <0..Kp-1>; caller allocated
 *            opt_fsc  - optRETURN: Forward score in nats
 *            ret_oasc - RETURN: OA score, the expected number of
 *                       correctly aligned residues
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; <eslEINVAL> if the
 *            arguments aren't set up for this computation.
 */
int
p7_GOptimalAccuracySweep_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs5, P7_CHKGMX *cx, P7_IVX *iv,
                                    P7_OATMX *ox, float *opt_null2, float *opt_fsc, float *ret_oasc)
{
  int     M     = gm_fs5->M;
  int     W     = (M+1) * p7G_NSCELLS_FS;
  float  *mem   = NULL;
  float **fdp   = NULL;
  float **odp   = NULL;
  float  *fxmx  = NULL;
  float  *oxmx  = NULL;
  P7_IVX *bv    = NULL;
  float  *pp;
  float  *cnt;
  float   ppx[p7G_NXCELLS];
  float   cntx[p7G_NXCELLS];
  float   overall_sc;
  int     i, k, b;
  int     status;

  if (gm_fs5->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");
  if (L < 3)                      ESL_EXCEPTION(eslEINVAL, "sequence length must be >= 3");
  if (ox->lanes != 1 || ox->nplanes != p7OAT_NPLANES_FS) ESL_EXCEPTION(eslEINVAL, "P7_OATMX not set up for the generic frameshift OA sweep");
  if (cx->L != L || cx->M != M || cx->dir != p7_CHK_BACKWARD) ESL_EXCEPTION(eslEINVAL, "checkpointed Backward matrix doesn't match");
  if ((status = p7_oatmx_GrowTo(ox, M, L)) != eslOK) return status;

  /* four Forward rows (i-1 and i-3 are read), six OA rows (i-1..i-5), one posterior row, the posterior sums for null2 */
  ESL_ALLOC(mem,  sizeof(float)   * W * 12);
  ESL_ALLOC(fdp,  sizeof(float *) * (L+1));
  ESL_ALLOC(odp,  sizeof(float *) * (L+1));
  ESL_ALLOC(fxmx, sizeof(float)   * (L+1) * p7G_NXCELLS);
  ESL_ALLOC(oxmx, sizeof(float)   * (L+1) * p7G_NXCELLS);
  if ((bv = p7_ivx_Create(M, p7P_5CODONS)) == NULL) { status = eslEMEM; goto ERROR; }
  pp  = mem + 10*W;
  cnt = mem + 11*W;

  p7_FLogsumInit();

  /* posteriors are normalized by the Backward score, as p7_GDecoding_Frameshift() does */
  overall_sc = p7_FLogsum( cx->xmx[0*p7G_NXCELLS+p7G_N],
               p7_FLogsum( cx->xmx[1*p7G_NXCELLS+p7G_N],
                           cx->xmx[2*p7G_NXCELLS+p7G_N]));

  fdp[0] = mem;
  odp[0] = mem + 4*W;
  p7_GForwardRow_Frameshift(dsq, gm_fs5, fdp, fxmx, iv->ivx, 0);
  oxmx[0*p7G_NXCELLS + p7G_N] = 0.;
  oxmx[0*p7G_NXCELLS + p7G_B] = 0.;
  oxmx[0*p7G_NXCELLS + p7G_E] = oxmx[0*p7G_NXCELLS + p7G_C] = oxmx[0*p7G_NXCELLS + p7G_J] = -eslINFINITY;
  for (k = 0; k <= M; k++)
    odp[0][k*p7G_NSCELLS_FS + p7G_M + p7G_C1] = odp[0][k*p7G_NSCELLS_FS + p7G_M + p7G_C2] = odp[0][k*p7G_NSCELLS_FS + p7G_D] = -eslINFINITY;
  oat_store_row_fs(gm_fs5, odp, oxmx, NULL, NULL, ox, 0);

  for (i = 1; i <= L; i++)
    {
      fdp[i] = mem + (i%4)*W;
      odp[i] = mem + (4 + i%6)*W;

      b = (i / cx->blockR) * cx->blockR;
      if (b != cx->b && (status = p7_GBackwardChkLoad_Frameshift(dsq, gm_fs5, cx, bv, b)) != eslOK) goto ERROR;

      p7_GForwardRow_Frameshift(dsq, gm_fs5, fdp, fxmx, iv->ivx, i);
      decode_fs_row(gm_fs5, fdp[i], (i > 2) ? fxmx + (i-3)*p7G_NXCELLS : NULL, cx->dp[i], cx->xmx + i*p7G_NXCELLS, overall_sc, pp, ppx, i);
      oa_fs_row(gm_fs5, pp, ppx, odp, oxmx, i);
      oat_store_row_fs(gm_fs5, odp, oxmx, pp, ppx, ox, i);

      /* posterior sums for null2, added up as p7_GNull2_fs_ByExpectation() does */
      if (i == 1) { esl_vec_FCopy(pp, W, cnt); esl_vec_FCopy(ppx, p7G_NXCELLS, cntx); }
      else        { esl_vec_FAdd (cnt, pp, W); esl_vec_FAdd (cntx, ppx, p7G_NXCELLS); }
    }

  if (opt_null2 != NULL) p7_GNull2_fs_FromCounts(gm_fs5, cnt, cntx, opt_null2);

  if (opt_fsc != NULL)
    *opt_fsc = p7_FLogsum( fxmx[L*p7G_NXCELLS + p7G_C],
               p7_FLogsum( fxmx[(L-1)*p7G_NXCELLS + p7G_C] + gm_fs5->xsc[p7P_C][p7P_LOOP],
                           fxmx[(L-2)*p7G_NXCELLS + p7G_C] + gm_fs5->xsc[p7P_C][p7P_LOOP])) +
               gm_fs5->xsc[p7P_C][p7P_MOVE];
  *ret_oasc = oxmx[L*p7G_NXCELLS + p7G_C] + oxmx[(L-1)*p7G_NXCELLS + p7G_C] + oxmx[(L-2)*p7G_NXCELLS + p7G_C];

  free(mem); free(fdp); free(odp); free(fxmx); free(oxmx); p7_ivx_Destroy(bv);
  return eslOK;

 ERROR:
  free(mem); free(fdp); free(odp); free(fxmx); free(oxmx); p7_ivx_Destroy(bv);
  return status;
}


/* Function:  p7_GOATraceSweep_Frameshift()
 * Synopsis:  Frameshift OA trace from what the OA sweep stored.
 *
 * Purpose:   The frameshift OA trace for profile <gm_fs5> from <ox>,
 *            filled by <p7_GOptimalAccuracySweep_Frameshift()>, into the
 *            empty trace <tr> (created with <p7_trace_fs_CreateWithPP()>).
 *            The trace is identical to <p7_GOATrace_Frameshift()>'s from
 *            the full OA matrix. The posterior probability of each
 *            residue is the value its PP character stands for
 *            (<p7OAT_PPVALUE()>), as <p7_alidisplay_EncodePostProb()>
 *            prints it.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if a choice
 *            fails or <tr> isn't empty.
 */
int
p7_GOATraceSweep_Frameshift(const P7_FS_PROFILE *gm_fs5, const P7_OATMX *ox, P7_TRACE *tr)
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

  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");

  if ((status = p7_trace_fs_AppendWithPP(tr, p7T_T, k, i, c, 0.0f)) != eslOK) return status;
  if ((status = p7_trace_fs_AppendWithPP(tr, p7T_C, k, i, c, 0.0f)) != eslOK) return status;

  sprv = p7T_C;
  while (sprv != p7T_S)
    {
      switch (sprv) {
      case p7T_M: scur = state[OATX_MPRED(ox, i, k)];                                k--;  break;
      case p7T_D: scur = OATX_DPRED(ox, i, k) ? p7T_D : p7T_M;                        k--;  break;
      case p7T_I: scur = (OATX_IPRED(ox, i-3, k) == p7OAT_FROM_I) ? p7T_I : p7T_M;    i-=3; break;
      case p7T_N: scur = (i == 0) ? p7T_S : p7T_N;                                          break;
      case p7T_C:
        t1 = ( (gm_fs5->xsc[p7P_C][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
        t2 = ( (gm_fs5->xsc[p7P_E][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
        if (i < 4) { scur = p7T_E; break; }
        path[0] = t1 * (OATX_XMX(ox, i-3, p7OAT_C) + OATX_XMX(ox, i, p7OAT_PPC));
        path[1] = (i < L)   ? t1 * (OATX_XMX(ox, i-2, p7OAT_C) + OATX_XMX(ox, i+1, p7OAT_PPC)) : FLT_MIN;
        path[2] = (i < L-1) ? t1 * (OATX_XMX(ox, i-1, p7OAT_C) + OATX_XMX(ox, i+2, p7OAT_PPC)) : FLT_MIN;
        path[3] = t2 * OATX_XMX(ox, i, p7OAT_E);
        scur = (esl_vec_FArgMax(path, 4) == 3) ? p7T_E : p7T_C;
        break;
      case p7T_J:
        t1 = ( (gm_fs5->xsc[p7P_J][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
        t2 = ( (gm_fs5->xsc[p7P_E][p7P_LOOP] == -eslINFINITY) ? FLT_MIN : 1.0);
        if (i <= 5) { scur = p7T_E; break; }
        path[0] = t1 * (OATX_XMX(ox, i, p7OAT_J) + OATX_XMX(ox, i, p7OAT_PPJ));
        path[1] = t2 *  OATX_XMX(ox, i, p7OAT_E);
        scur = (esl_vec_FArgMax(path, 2) == 0) ? p7T_J : p7T_E;
        break;
      case p7T_E:
        if (ox->ek[i] == 0) ESL_EXCEPTION(eslEINVAL, "OA traceback choice failed");
        scur = (ox->ek[i] > 0) ? p7T_M : p7T_D;
        k    = abs(ox->ek[i]);
        break;
      case p7T_B:
        t1   = ( (gm_fs5->xsc[p7P_N][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
        t2   = ( (gm_fs5->xsc[p7P_J][p7P_MOVE] == -eslINFINITY) ? FLT_MIN : 1.0);
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
      c = (scur == p7T_M) ? OATX_CODON(ox, i, k) : 0;

      if ((status = p7_trace_fs_AppendWithPP(tr, scur, k, i, c, postprob)) != eslOK) return status;

      if ( (scur == p7T_N || scur == p7T_C || scur == p7T_J) && scur == sprv) i--;
      sprv = scur;
      i   -= c;
    }
  tr->M = gm_fs5->M;
  tr->L = L;
  return p7_trace_fs_Reverse(tr);
}

/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef p7GENERIC_OPTACC_SWEEP_FRAMESHIFT_TESTDRIVE
#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_sweep_fs()
 * The frameshift OA sweep and its trace against the full-matrix path
 * (p7_GForward_Frameshift, p7_GBackward_Frameshift,
 * p7_GDecoding_Frameshift, p7_GOptimalAccuracy_Frameshift,
 * p7_GOATrace_Frameshift) on random DNA and back-translated sequences
 * emitted from the model, with Backward blocks of 5..12 rows and of the
 * whole sequence: same OA score and null2, same states, k, i and codon
 * lengths on the trace, and each residue's PP character equal to the one
 * the posterior matrix gives it.
 */
static void
utest_sweep_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA,
               P7_CODONTABLE *codon_table, int M, int L, int N, int mode)
{
  char          *msg  = "generic frameshift OA sweep unit test failed";
  P7_HMM        *hmm  = NULL;
  P7_PROFILE    *gm   = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE *gm5  = p7_profile_fs_Create(M, abcAA, p7P_5CODONS);
  ESL_SQ        *sq   = esl_sq_CreateDigital(abcAA);
  ESL_DSQ       *dsq  = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_GMX        *fwd  = p7_gmx_Create(M, L, L, p7G_NSCELLS_FS);
  P7_GMX        *bck  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX        *oa   = p7_gmx_Create(M, L, L, p7G_NSCELLS_FS);
  P7_CHKGMX     *cx   = p7_chkgmx_Create(M, L, 5, 5, p7G_NSCELLS, p7_CHK_BACKWARD);
  P7_OATMX      *ox   = p7_oatmx_Create(M, L, p7OAT_NPLANES_FS, 1);
  P7_IVX        *iv   = p7_ivx_Create(M, p7P_5CODONS);
  P7_TRACE      *tr   = p7_trace_Create();
  P7_TRACE      *tr1  = p7_trace_fs_CreateWithPP();
  P7_TRACE      *tr2  = p7_trace_fs_CreateWithPP();
  float         *null2a = malloc(sizeof(float) * abcAA->Kp);
  float         *null2b = malloc(sizeof(float) * abcAA->Kp);
  float          fsc, fsc2, oasc1, oasc2, p;
  int            n, Ln, z, i, j, blockR, x;

  if (p7_hmm_Sample(r, M, abcAA, &hmm) != eslOK) esl_fatal(msg);
  if (p7_ProfileConfig(hmm, bgAA, gm, L/3, mode) != eslOK) esl_fatal(msg);
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm5, L/3, mode) != eslOK) esl_fatal(msg);

  for (n = 0; n < N; n++)
    {
      if (n % 2 == 0)
        {
          Ln = 15 + esl_rnd_Roll(r, L - 14);
          esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, dsq);
        }
      else
        {
          do { esl_sq_Reuse(sq); p7_trace_Reuse(tr); p7_ProfileEmit(r, hmm, gm, bgAA, sq, tr); } while (sq->n < 5 || 3*sq->n > L);
          Ln = 3 * sq->n;
          for (i = 1, j = 1; i <= sq->n; i++, j += 3)
            p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq+j);
          /* now and then a frameshift: delete or insert a nucleotide */
          if (n % 4 == 1 && Ln > 20)
            {
              x = 5 + esl_rnd_Roll(r, Ln - 10);
              memmove(dsq + x, dsq + x + 1, Ln - x);
              Ln--;
            }
        }
      dsq[0] = dsq[Ln+1] = eslDSQ_SENTINEL;
      p7_fs_ReconfigLength(gm5, Ln/3);

      /* full-matrix path */
      p7_GForward_Frameshift (dsq, Ln, gm5, fwd, iv, &fsc);
      p7_GBackward_Frameshift(dsq, Ln, gm5, bck, iv, NULL);
      p7_GDecoding_Frameshift(gm5, fwd, bck);                      /* <fwd> now holds the posteriors */
      p7_GOptimalAccuracy_Frameshift(gm5, fwd, oa, &oasc1);
      p7_trace_Reuse(tr1);
      if (p7_GOATrace_Frameshift(gm5, fwd->xmx, oa, tr1) != eslOK) esl_fatal(msg);
      p7_GNull2_fs_ByExpectation(gm5, fwd, null2a);              /* uses row 0 of <fwd> as workspace */

      /* sweep */
      blockR = (n % 3 == 2) ? Ln+1 : 5 + esl_rnd_Roll(r, 8);
      if (p7_chkgmx_GrowTo(cx, M, Ln, blockR)                       != eslOK) esl_fatal(msg);
      if (p7_GBackwardChk_Frameshift(dsq, Ln, gm5, cx, iv, NULL)     != eslOK) esl_fatal(msg);
      if (p7_GOptimalAccuracySweep_Frameshift(dsq, Ln, gm5, cx, iv, ox, null2b, &fsc2, &oasc2) != eslOK) esl_fatal(msg);
      p7_trace_Reuse(tr2);
      if (p7_GOATraceSweep_Frameshift(gm5, ox, tr2) != eslOK) esl_fatal(msg);

      if (fsc2  != fsc)   esl_fatal("%s: Forward scores differ (%f %f)", msg, fsc2, fsc);
      if (oasc2 != oasc1) esl_fatal("%s: OA scores differ (%f %f)", msg, oasc2, oasc1);
      for (x = 0; x < abcAA->Kp; x++)
        if (null2b[x] != null2a[x]) esl_fatal("%s: null2 differs for residue %d (%g %g)", msg, x, null2b[x], null2a[x]);
      if (tr1->N != tr2->N) esl_fatal("%s: trace lengths differ (%d %d)", msg, tr1->N, tr2->N);
      for (z = 0; z < tr1->N; z++)
        {
          if (tr1->st[z] != tr2->st[z] || tr1->k[z] != tr2->k[z] || tr1->i[z] != tr2->i[z] || tr1->c[z] != tr2->c[z])
            esl_fatal("%s: traces differ at %d", msg, z);
          switch (tr1->st[z]) {
          case p7T_M: p = fwd->dp[tr1->i[z]][tr1->k[z]*p7G_NSCELLS_FS + p7G_M + p7G_C0]; break;
          case p7T_I: p = fwd->dp[tr1->i[z]][tr1->k[z]*p7G_NSCELLS_FS + p7G_I];          break;
          case p7T_N: case p7T_C: case p7T_J:
            /* an N, C or J emits residue i when the state before it is the same one */
            p = (z > 0 && tr1->st[z-1] == tr1->st[z]) ? fwd->xmx[tr1->i[z]*p7G_NXCELLS + (tr1->st[z] == p7T_N ? p7G_N : tr1->st[z] == p7T_C ? p7G_C : p7G_J)] : 0.0;
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
  p7_trace_Destroy(tr);
  p7_trace_fs_Destroy(tr1);
  p7_trace_fs_Destroy(tr2);
  p7_ivx_Destroy(iv);
  p7_oatmx_Destroy(ox);
  p7_chkgmx_Destroy(cx);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(bck);
  p7_gmx_Destroy(oa);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm5);
  p7_hmm_Destroy(hmm);
}
#endif /*p7GENERIC_OPTACC_SWEEP_FRAMESHIFT_TESTDRIVE*/


/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef p7GENERIC_OPTACC_SWEEP_FRAMESHIFT_TESTDRIVE
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
static char banner[] = "unit test driver for the generic frameshift OA sweep";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go  = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r   = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc = esl_alphabet_Create(eslAMINO);
  ESL_ALPHABET   *dna = esl_alphabet_Create(eslDNA);
  ESL_GENCODE    *gc  = esl_gencode_Create(dna, abc);
  P7_CODONTABLE  *ct  = p7_codontable_Create(gc);
  P7_BG          *bg  = p7_bg_Create(abc);
  P7_BG          *bgd = p7_bg_Create(dna);
  int             M   = esl_opt_GetInteger(go, "-M");
  int             L   = esl_opt_GetInteger(go, "-L");
  int             N   = esl_opt_GetInteger(go, "-N");

  p7_FLogsumInit();

  utest_sweep_fs(r, abc, dna, gc, bg, bgd, ct, M, 3*L, N,  p7_LOCAL);
  utest_sweep_fs(r, abc, dna, gc, bg, bgd, ct, M, 3*L, N,  p7_UNILOCAL);
  utest_sweep_fs(r, abc, dna, gc, bg, bgd, ct, M, 3*L, N,  p7_GLOCAL);
  utest_sweep_fs(r, abc, dna, gc, bg, bgd, ct, 1, 60,  10, p7_LOCAL);

  p7_codontable_Destroy(ct);
  esl_gencode_Destroy(gc);
  p7_bg_Destroy(bgd);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(dna);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_OPTACC_SWEEP_FRAMESHIFT_TESTDRIVE*/
