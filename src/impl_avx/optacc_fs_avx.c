/* Frameshift optimal accuracy alignment; AVX2 implementation.
 * Ported from optacc_fs_sse.c with _avx suffix for runtime dispatch.
 *
 * Contents:
 *   1. p7_OptimalAccuracy_Frameshift_avx() — DP fill
 *   2. p7_OATrace_Frameshift_avx()         — traceback
 */
#include "p7_config.h"

#ifdef eslENABLE_AVX

#include <float.h>

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

/* avx_hmax_ps: horizontal max across all 8 lanes of a __m256, return scalar. */
static inline void
avx_hmax_ps(__m256 v, float *ret_max)
{
  __m128 a = _mm256_extractf128_ps(v, 0);
  __m128 b = _mm256_extractf128_ps(v, 1);
  a = _mm_max_ps(a, b);
  a = _mm_max_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(2,3,0,1)));
  a = _mm_max_ps(a, _mm_shuffle_ps(a, a, _MM_SHUFFLE(1,0,3,2)));
  _mm_store_ss(ret_max, a);
}


/*****************************************************************
 * 1. Optimal accuracy alignment, DP fill
 *****************************************************************/

/* Function:  p7_OptimalAccuracy_Frameshift_avx()
 *
 * Purpose:   AVX2 implementation of the DP fill step of optimal accuracy
 *            alignment for the frameshift model.  Five codon lengths (1-5 nt)
 *            are supported.  The posterior decoding matrix <pp> uses the
 *            8-cell FS layout; the OA result matrix <ox> uses the standard
 *            3-cell layout.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_OptimalAccuracy_Frameshift_avx(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                                   P7_OMX *ox, float *ret_e)
{
  int    M   = om_fs->M;
  int    L   = pp->L;
  int    Q   = p7O_NQF_AVX(M);
  int    i, q, j;
  int    qM  = (M-1) % Q;   /* stripe containing k=M */
  int    rM  = (M-1) / Q;   /* lane  containing k=M */
  float  xN, xE, xJ, xC;
  float  t1, t2;
  register __m256 mpv1, dpv1, ipv1;   /* right-shifted prev row i-1, then updated to current q */
  register __m256 mpv2, dpv2, ipv2;   /* right-shifted prev row i-2 */
  register __m256 mpv3, dpv3, ipv3;   /* right-shifted prev row i-3 */
  register __m256 mpv4, dpv4, ipv4;   /* right-shifted prev row i-4 */
  register __m256 mpv5, dpv5, ipv5;   /* right-shifted prev row i-5 */
  register __m256 sv1, sv2, sv3, sv4, sv5, sv;
  register __m256 xEv, xBv1, xBv2, xBv3, xBv4, xBv5;
  register __m256 dcv;
  __m256 bm, mm, im, dm, md, mi, ii;
  __m256 *tp;
  __m256 *dpc;                          /* current OA row (3-cell layout)     */
  __m256 *dp1, *dp2, *dp3, *dp4, *dp5; /* previous OA rows                   */
  __m256 *ppc;                          /* current pp row (8-cell FS layout)  */
  __m256  zerov = _mm256_setzero_ps();
  __m256  infv  = _mm256_set1_ps(-eslINFINITY);
  union { __m256 v; float p[8]; } imfix;   /* for I(i,M) = -inf fixup */

  ox->M = M;
  ox->L = L;

  /* Initialize row 0: MDI = -inf; N=0, B=0, E/J/C = -inf */
  dpc = ox->dpf_avx[0];
  for (q = 0; q < Q; q++) { MMO(dpc,q) = DMO(dpc,q) = IMO(dpc,q) = infv; }
  ox->xmx[0*p7X_NXCELLS+p7X_E] = -eslINFINITY;
  ox->xmx[0*p7X_NXCELLS+p7X_N] =  0.0f;
  ox->xmx[0*p7X_NXCELLS+p7X_J] = -eslINFINITY;
  ox->xmx[0*p7X_NXCELLS+p7X_B] =  0.0f;
  ox->xmx[0*p7X_NXCELLS+p7X_C] = -eslINFINITY;

  for (i = 1; i <= L; i++)
    {
      dpc = ox->dpf_avx[i];
      ppc = pp->dpf_avx[i];

      /* Previous OA rows; row 0 (all -inf MDI) serves as sentinel for out-of-bounds lags */
      dp1 = ox->dpf_avx[i-1];
      dp2 = (i >= 2) ? ox->dpf_avx[i-2] : ox->dpf_avx[0];
      dp3 = (i >= 3) ? ox->dpf_avx[i-3] : ox->dpf_avx[0];
      dp4 = (i >= 4) ? ox->dpf_avx[i-4] : ox->dpf_avx[0];
      dp5 = (i >= 5) ? ox->dpf_avx[i-5] : ox->dpf_avx[0];

      /* B values from previous rows; -inf sentinel when lag exceeds i */
      xBv1 = _mm256_set1_ps(ox->xmx[(i-1)*p7X_NXCELLS + p7X_B]);
      xBv2 = (i >= 2) ? _mm256_set1_ps(ox->xmx[(i-2)*p7X_NXCELLS + p7X_B]) : infv;
      xBv3 = (i >= 3) ? _mm256_set1_ps(ox->xmx[(i-3)*p7X_NXCELLS + p7X_B]) : infv;
      xBv4 = (i >= 4) ? _mm256_set1_ps(ox->xmx[(i-4)*p7X_NXCELLS + p7X_B]) : infv;
      xBv5 = (i >= 5) ? _mm256_set1_ps(ox->xmx[(i-5)*p7X_NXCELLS + p7X_B]) : infv;

      /* Right-shift last stripe of each prev row for k=0 (k-1 access at q=0) */
      mpv1 = avx_rightshift_ps(MMO(dp1, Q-1), infv);
      dpv1 = avx_rightshift_ps(DMO(dp1, Q-1), infv);
      ipv1 = avx_rightshift_ps(IMO(dp1, Q-1), infv);

      mpv2 = avx_rightshift_ps(MMO(dp2, Q-1), infv);
      dpv2 = avx_rightshift_ps(DMO(dp2, Q-1), infv);
      ipv2 = avx_rightshift_ps(IMO(dp2, Q-1), infv);

      mpv3 = avx_rightshift_ps(MMO(dp3, Q-1), infv);
      dpv3 = avx_rightshift_ps(DMO(dp3, Q-1), infv);
      ipv3 = avx_rightshift_ps(IMO(dp3, Q-1), infv);

      mpv4 = avx_rightshift_ps(MMO(dp4, Q-1), infv);
      dpv4 = avx_rightshift_ps(DMO(dp4, Q-1), infv);
      ipv4 = avx_rightshift_ps(IMO(dp4, Q-1), infv);

      mpv5 = avx_rightshift_ps(MMO(dp5, Q-1), infv);
      dpv5 = avx_rightshift_ps(DMO(dp5, Q-1), infv);
      ipv5 = avx_rightshift_ps(IMO(dp5, Q-1), infv);

      dcv = infv;
      xEv = infv;
      tp  = om_fs->tfv_avx;

      for (q = 0; q < Q; q++)
        {
          /* Read 7 transition masks per stripe: BM, MM, IM, DM, MD, MI, II */
          bm = _mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ);  tp++;
          mm = _mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ);  tp++;
          im = _mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ);  tp++;
          dm = _mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ);  tp++;
          md = _mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ);  tp++;
          mi = _mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ);  tp++;
          ii = _mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ);  tp++;

          /* M contributions for codon lag 1 (from row i-1, k-1 access) */
          sv1 =                  _mm256_and_ps(bm, xBv1);
          sv1 = _mm256_max_ps(sv1, _mm256_and_ps(mm, mpv1));
          sv1 = _mm256_max_ps(sv1, _mm256_and_ps(im, ipv1));
          sv1 = _mm256_max_ps(sv1, _mm256_and_ps(dm, dpv1));
          sv1 = _mm256_add_ps(sv1, MMO_FS(ppc, q, p7X_FS_C1));

          /* M contributions for codon lag 2 (from row i-2, k-1 access) */
          sv2 =                  _mm256_and_ps(bm, xBv2);
          sv2 = _mm256_max_ps(sv2, _mm256_and_ps(mm, mpv2));
          sv2 = _mm256_max_ps(sv2, _mm256_and_ps(im, ipv2));
          sv2 = _mm256_max_ps(sv2, _mm256_and_ps(dm, dpv2));
          sv2 = _mm256_add_ps(sv2, MMO_FS(ppc, q, p7X_FS_C2));

          /* M contributions for codon lag 3 (from row i-3, k-1 access) */
          sv3 =                  _mm256_and_ps(bm, xBv3);
          sv3 = _mm256_max_ps(sv3, _mm256_and_ps(mm, mpv3));
          sv3 = _mm256_max_ps(sv3, _mm256_and_ps(im, ipv3));
          sv3 = _mm256_max_ps(sv3, _mm256_and_ps(dm, dpv3));
          sv3 = _mm256_add_ps(sv3, MMO_FS(ppc, q, p7X_FS_C3));

          /* M contributions for codon lag 4 (from row i-4, k-1 access) */
          sv4 =                  _mm256_and_ps(bm, xBv4);
          sv4 = _mm256_max_ps(sv4, _mm256_and_ps(mm, mpv4));
          sv4 = _mm256_max_ps(sv4, _mm256_and_ps(im, ipv4));
          sv4 = _mm256_max_ps(sv4, _mm256_and_ps(dm, dpv4));
          sv4 = _mm256_add_ps(sv4, MMO_FS(ppc, q, p7X_FS_C4));

          /* M contributions for codon lag 5 (from row i-5, k-1 access) */
          sv5 =                  _mm256_and_ps(bm, xBv5);
          sv5 = _mm256_max_ps(sv5, _mm256_and_ps(mm, mpv5));
          sv5 = _mm256_max_ps(sv5, _mm256_and_ps(im, ipv5));
          sv5 = _mm256_max_ps(sv5, _mm256_and_ps(dm, dpv5));
          sv5 = _mm256_add_ps(sv5, MMO_FS(ppc, q, p7X_FS_C5));

          /* Best M(i,k) over all 5 codon lengths */
          sv = _mm256_max_ps(_mm256_max_ps(sv1, sv2),
                             _mm256_max_ps(_mm256_max_ps(sv3, sv4), sv5));

          xEv = _mm256_max_ps(xEv, sv);

          /* Update prev-row ptrs from right-shifted (k-1) to current q (k access) */
          mpv1 = MMO(dp1, q);  dpv1 = DMO(dp1, q);  ipv1 = IMO(dp1, q);
          mpv2 = MMO(dp2, q);  dpv2 = DMO(dp2, q);  ipv2 = IMO(dp2, q);
          mpv3 = MMO(dp3, q);  dpv3 = DMO(dp3, q);  ipv3 = IMO(dp3, q);
          mpv4 = MMO(dp4, q);  dpv4 = DMO(dp4, q);  ipv4 = IMO(dp4, q);
          mpv5 = MMO(dp5, q);  dpv5 = DMO(dp5, q);  ipv5 = IMO(dp5, q);

          MMO(dpc, q) = sv;
          DMO(dpc, q) = dcv;
          dcv = _mm256_and_ps(md, sv);

          /* I(i,k) = max(MI*M(i-3,k), II*I(i-3,k)) + pp_I(i,k); lag-3, same column */
          sv  =                  _mm256_and_ps(mi, mpv3);
          sv  = _mm256_max_ps(sv, _mm256_and_ps(ii, ipv3));
          IMO(dpc, q) = _mm256_add_ps(sv, IMO_FS(ppc, q));
        }

      /* I(i,M) has no I state; explicitly set to -inf */
      imfix.v      = IMO(dpc, qM);
      imfix.p[rM]  = -eslINFINITY;
      IMO(dpc, qM) = imfix.v;

      /* DD propagation: first pass */
      dcv = avx_rightshift_ps(dcv, infv);
      tp  = om_fs->tfv_avx + 7*Q;
      for (q = 0; q < Q; q++)
        {
          DMO(dpc, q) = _mm256_max_ps(dcv, DMO(dpc, q));
          dcv         = _mm256_and_ps(_mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ), DMO(dpc, q));  tp++;
        }

      /* DD propagation: 7 more lazy passes (AVX: 8 lanes, need 8 total passes) */
      for (j = 1; j < 8; j++)
        {
          dcv = avx_rightshift_ps(dcv, infv);
          tp  = om_fs->tfv_avx + 7*Q;
          for (q = 0; q < Q; q++)
            {
              DMO(dpc, q) = _mm256_max_ps(dcv, DMO(dpc, q));
              dcv         = _mm256_and_ps(_mm256_cmp_ps(*tp, zerov, _CMP_GT_OQ), dcv);  tp++;
            }
        }

      /* D->E */
      for (q = 0; q < Q; q++) xEv = _mm256_max_ps(xEv, DMO(dpc, q));

      xE = -eslINFINITY;
      avx_hmax_ps(xEv, &xE);
      ox->xmx[i*p7X_NXCELLS+p7X_E] = xE;

      /* N state: lag-3 */
      if (i > 2)
        xN = (om_fs->xf[p7O_N][p7O_LOOP] == 0.0f) ? 0.0f :
             ox->xmx[(i-3)*p7X_NXCELLS+p7X_N] + pp->xmx[i*p7X_NXCELLS+p7X_N];
      else
        xN = (om_fs->xf[p7O_N][p7O_LOOP] == 0.0f) ? 0.0f :
             pp->xmx[i*p7X_NXCELLS+p7X_N];
      ox->xmx[i*p7X_NXCELLS+p7X_N] = xN;

      /* J state: lag-3 */
      if (i > 2) {
        t1 = (om_fs->xf[p7O_J][p7O_LOOP] == 0.0f) ? 0.0f :
             ox->xmx[(i-3)*p7X_NXCELLS+p7X_J] + pp->xmx[i*p7X_NXCELLS+p7X_J];
        t2 = (om_fs->xf[p7O_E][p7O_LOOP] == 0.0f) ? 0.0f : xE;
        xJ = ESL_MAX(t1, t2);
      } else {
        xJ = (om_fs->xf[p7O_E][p7O_LOOP] == 0.0f) ? 0.0f : xE;
      }
      ox->xmx[i*p7X_NXCELLS+p7X_J] = xJ;

      /* C state: lag-3 */
      if (i > 2) {
        t1 = (om_fs->xf[p7O_C][p7O_LOOP] == 0.0f) ? 0.0f :
             ox->xmx[(i-3)*p7X_NXCELLS+p7X_C] + pp->xmx[i*p7X_NXCELLS+p7X_C];
        t2 = (om_fs->xf[p7O_E][p7O_MOVE] == 0.0f) ? 0.0f : xE;
        xC = ESL_MAX(t1, t2);
      } else {
        xC = (om_fs->xf[p7O_E][p7O_MOVE] == 0.0f) ? 0.0f : xE;
      }
      ox->xmx[i*p7X_NXCELLS+p7X_C] = xC;

      /* B state */
      t1 = (om_fs->xf[p7O_N][p7O_MOVE] == 0.0f) ? 0.0f : xN;
      t2 = (om_fs->xf[p7O_J][p7O_MOVE] == 0.0f) ? 0.0f : xJ;
      ox->xmx[i*p7X_NXCELLS+p7X_B] = ESL_MAX(t1, t2);
    }

  *ret_e = ox->xmx[L    *p7X_NXCELLS+p7X_C] +
           ox->xmx[(L-1)*p7X_NXCELLS+p7X_C] +
           ox->xmx[(L-2)*p7X_NXCELLS+p7X_C];
  return eslOK;
}


/*****************************************************************
 * 2. OA traceback
 *****************************************************************/

static inline float
get_postprob_fs(const P7_OMX *pp, int scur, int sprv, int k, int i)
{
  int   Q = p7O_NQF_AVX(pp->M);
  int   q = (k-1) % Q;
  int   r = (k-1) / Q;
  union { __m256 v; float p[8]; } u;

  switch (scur) {
  case p7T_M:
    u.v = MMO_FS(pp->dpf_avx[i], q, p7X_FS_C0);  return u.p[r];
  case p7T_I:
    u.v = IMO_FS(pp->dpf_avx[i], q);              return u.p[r];
  case p7T_N: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_N]; break;
  case p7T_C: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_C]; break;
  case p7T_J: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_J]; break;
  default: break;
  }
  return 0.0f;
}

static inline int
select_m_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF_AVX(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m256 *tp    = om_fs->tfv_avx + 7*q;
  union { __m256 v; float p[8]; } u, tv;
  float   path[4];
  int     state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };
  __m256  mpv, dpv, ipv;

  if (q > 0) {
    mpv = ox->dpf_avx[i][(q-1)*p7X_NSCELLS + p7X_M];
    dpv = ox->dpf_avx[i][(q-1)*p7X_NSCELLS + p7X_D];
    ipv = ox->dpf_avx[i][(q-1)*p7X_NSCELLS + p7X_I];
  } else {
    __m256 neg_inf = _mm256_set1_ps(-eslINFINITY);
    mpv = avx_rightshift_ps(ox->dpf_avx[i][(Q-1)*p7X_NSCELLS + p7X_M], neg_inf);
    dpv = avx_rightshift_ps(ox->dpf_avx[i][(Q-1)*p7X_NSCELLS + p7X_D], neg_inf);
    ipv = avx_rightshift_ps(ox->dpf_avx[i][(Q-1)*p7X_NSCELLS + p7X_I], neg_inf);
  }

  tv.v = *tp;  path[3] = ((tv.p[r] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_B]);  tp++;
  u.v  = mpv;  tv.v = *tp;  path[0] = ((tv.p[r] == 0.0f) ? -eslINFINITY : u.p[r]);  tp++;
  u.v  = ipv;  tv.v = *tp;  path[1] = ((tv.p[r] == 0.0f) ? -eslINFINITY : u.p[r]);  tp++;
  u.v  = dpv;  tv.v = *tp;  path[2] = ((tv.p[r] == 0.0f) ? -eslINFINITY : u.p[r]);
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
select_d_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF_AVX(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  union { __m256 v; float p[8]; } mpv, dpv, tmdv, tddv;
  float   path[2];

  if (q > 0) {
    mpv.v  = ox->dpf_avx[i][(q-1)*p7X_NSCELLS + p7X_M];
    dpv.v  = ox->dpf_avx[i][(q-1)*p7X_NSCELLS + p7X_D];
    tmdv.v = om_fs->tfv_avx[7*(q-1) + p7O_MD];
    tddv.v = om_fs->tfv_avx[7*Q     + (q-1)];
  } else {
    __m256 neg_inf = _mm256_set1_ps(-eslINFINITY);
    mpv.v  = avx_rightshift_ps(ox->dpf_avx[i][(Q-1)*p7X_NSCELLS + p7X_M], neg_inf);
    dpv.v  = avx_rightshift_ps(ox->dpf_avx[i][(Q-1)*p7X_NSCELLS + p7X_D], neg_inf);
    tmdv.v = esl_avx_rightshiftz_float(om_fs->tfv_avx[7*(Q-1) + p7O_MD]);
    tddv.v = esl_avx_rightshiftz_float(om_fs->tfv_avx[8*Q - 1]);
  }

  path[0] = ((tmdv.p[r] == 0.0f) ? -eslINFINITY : mpv.p[r]);
  path[1] = ((tddv.p[r] == 0.0f) ? -eslINFINITY : dpv.p[r]);
  return ((path[0] >= path[1]) ? p7T_M : p7T_D);
}

static inline int
select_i_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q      = p7O_NQF_AVX(ox->M);
  int     q      = (k-1) % Q;
  int     r      = (k-1) / Q;
  int     prev_i = (i >= 3) ? i-3 : 0;
  __m256 *tp     = om_fs->tfv_avx + 7*q + p7O_MI;
  union { __m256 v; float p[8]; } tv, mpv, ipv;
  float   path[2];

  mpv.v = ox->dpf_avx[prev_i][q*p7X_NSCELLS + p7X_M];
  tv.v  = *tp;  path[0] = ((tv.p[r] == 0.0f) ? -eslINFINITY : mpv.p[r]);  tp++;
  ipv.v = ox->dpf_avx[prev_i][q*p7X_NSCELLS + p7X_I];
  tv.v  = *tp;  path[1] = ((tv.p[r] == 0.0f) ? -eslINFINITY : ipv.p[r]);
  return ((path[0] >= path[1]) ? p7T_M : p7T_I);
}

static inline int
select_n_fs(int i)
{
  return ((i == 0) ? p7T_S : p7T_N);
}

static inline int
select_c_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp, const P7_OMX *ox, int i)
{
  int   L     = ox->L;
  float t1    = (om_fs->xf[p7O_C][p7O_LOOP] == 0.0f) ? 0.0f : 1.0f;
  float t2    = (om_fs->xf[p7O_E][p7O_MOVE] == 0.0f) ? 0.0f : 1.0f;
  float path[4];
  int   state[4] = { p7T_C, p7T_C, p7T_C, p7T_E };

  if (i < 4) return p7T_E;

  path[0] = (t1 == 0.0f) ? -eslINFINITY :
            ox->xmx[(i-3)*p7X_NXCELLS+p7X_C] + pp->xmx[ i   *p7X_NXCELLS+p7X_C];
  path[1] = (i < L && t1 != 0.0f) ?
            ox->xmx[(i-2)*p7X_NXCELLS+p7X_C] + pp->xmx[(i+1)*p7X_NXCELLS+p7X_C] : -eslINFINITY;
  path[2] = (i < L-1 && t1 != 0.0f) ?
            ox->xmx[(i-1)*p7X_NXCELLS+p7X_C] + pp->xmx[(i+2)*p7X_NXCELLS+p7X_C] : -eslINFINITY;
  path[3] = (t2 == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_E];
  return state[esl_vec_FArgMax(path, 4)];
}

static inline int
select_j_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp, const P7_OMX *ox, int i)
{
  float path[2];
  int   state[2] = { p7T_J, p7T_E };

  if (i <= 5) return p7T_E;

  path[0] = (om_fs->xf[p7O_J][p7O_LOOP] == 0.0f) ? -eslINFINITY :
            ox->xmx[i*p7X_NXCELLS+p7X_J] + pp->xmx[i*p7X_NXCELLS+p7X_J];
  path[1] = (om_fs->xf[p7O_E][p7O_LOOP] == 0.0f) ? -eslINFINITY :
            ox->xmx[i*p7X_NXCELLS+p7X_E];
  return state[esl_vec_FArgMax(path, 2)];
}

static inline int
select_e_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int *ret_k)
{
  int    Q    = p7O_NQF_AVX(ox->M);
  union { __m256 v; float p[8]; } u;
  float  max  = -eslINFINITY;
  int    smax = p7T_M, kmax = 1;
  int    k, q, r;

  for (k = 1; k <= ox->M; k++) {
    q = (k-1) % Q;
    r = (k-1) / Q;
    u.v = ox->dpf_avx[i][q * p7X_NSCELLS + p7X_M];
    if (u.p[r] > max) { max = u.p[r]; smax = p7T_M; kmax = k; }
    u.v = ox->dpf_avx[i][q * p7X_NSCELLS + p7X_D];
    if (u.p[r] > max) { max = u.p[r]; smax = p7T_D; kmax = k; }
  }
  *ret_k = kmax;
  return smax;
}

static inline int
select_b_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i)
{
  float path[2];
  path[0] = (om_fs->xf[p7O_N][p7O_MOVE] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_N];
  path[1] = (om_fs->xf[p7O_J][p7O_MOVE] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_J];
  return ((path[0] > path[1]) ? p7T_N : p7T_J);
}

static inline int
select_codon_fs(const P7_OMX *pp, int i, int k)
{
  int   Q = p7O_NQF_AVX(pp->M);
  int   q = (k-1) % Q;
  int   r = (k-1) / Q;
  union { __m256 v; float p[8]; } u;
  float codon[5];

  u.v = MMO_FS(pp->dpf_avx[i], q, p7X_FS_C1);  codon[0] = u.p[r];
  u.v = MMO_FS(pp->dpf_avx[i], q, p7X_FS_C2);  codon[1] = u.p[r];
  u.v = MMO_FS(pp->dpf_avx[i], q, p7X_FS_C3);  codon[2] = u.p[r];
  u.v = MMO_FS(pp->dpf_avx[i], q, p7X_FS_C4);  codon[3] = u.p[r];
  u.v = MMO_FS(pp->dpf_avx[i], q, p7X_FS_C5);  codon[4] = u.p[r];
  return esl_vec_FArgMax(codon, 5) + 1;
}


/* Function:  p7_OATrace_Frameshift_avx()
 *
 * Purpose:   AVX2 implementation of optimal accuracy traceback for the
 *            frameshift model.  Traces back through OA matrix <ox> annotated
 *            with posterior probabilities from FS decoding matrix <pp>.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if trace is not empty or a bogus state is reached.
 *            <eslEMEM> on allocation error.
 */
int
p7_OATrace_Frameshift_avx(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                           const P7_OMX *ox, P7_TRACE *tr)
{
  int   i    = ox->L;
  int   k    = 0;
  int   c    = 0;
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
      case p7T_M: scur = select_m_fs(om_fs,     ox, i, k);       k--;    break;
      case p7T_D: scur = select_d_fs(om_fs,     ox, i, k);       k--;    break;
      case p7T_I: scur = select_i_fs(om_fs,     ox, i, k);       i -= 3; break;
      case p7T_N: scur = select_n_fs(i);                                  break;
      case p7T_C: scur = select_c_fs(om_fs, pp, ox, i);                   break;
      case p7T_J: scur = select_j_fs(om_fs, pp, ox, i);                   break;
      case p7T_E: scur = select_e_fs(om_fs,     ox, i, &k);               break;
      case p7T_B: scur = select_b_fs(om_fs,     ox, i);                   break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in OA FS traceback");
      }
      if (scur == -1) ESL_EXCEPTION(eslEINVAL, "OA FS AVX traceback choice failed");

      postprob = get_postprob_fs(pp, scur, sprv, k, i);

      if (scur == p7T_M) c = select_codon_fs(pp, i, k);
      else               c = 0;

      if ((status = p7_trace_fs_AppendWithPP(tr, scur, k, i, c, postprob)) != eslOK) return status;

      if ( (scur == p7T_N || scur == p7T_C || scur == p7T_J) && scur == sprv) i--;
      sprv = scur;
      i   -= c;
    }

  tr->M = om_fs->M;
  tr->L = ox->L;
  return p7_trace_fs_Reverse(tr);
}

#endif /* eslENABLE_AVX */
