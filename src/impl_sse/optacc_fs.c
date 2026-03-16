/* Optimal accuracy framshift alignment; SSE version.
 * 
 * Contents:
 *   1. Optimal accuracy alignment, DP fill
 *   2. OA traceback
 *   4. Unit tests
 *   5. Test driver
 * 
 */
#include "p7_config.h"

#include <float.h>

#include <xmmintrin.h>
#include <emmintrin.h>

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_sse.h"


/*****************************************************************
 * 1. Optimal accuracy alignment, DP fill
 *****************************************************************/

/* Function:  p7_OptimalAccuracy_Frameshift_SSE()
 * Synopsis:  DP fill of optimal accuracy alignment, frameshift SSE version.
 *
 * Purpose:   Calculates the fill step of the optimal accuracy decoding
 *            algorithm for a frameshift-aware profile, using SSE SIMD
 *            parallelism over model positions.
 *
 *            The posterior decoding matrix <pp> must have the 8-cell FS
 *            layout (p7X_NSCELLS_FS=8; from p7_Decoding_Frameshift_SSE).
 *            The OA result matrix <ox> uses the standard 3-cell layout
 *            (p7X_NSCELLS=3; M, D, I per stripe position).
 *
 *            Five codon lengths (1-5 nt) are supported. M(i,k) takes
 *            contributions from rows i-1..i-5; I(i,k) from row i-3;
 *            special states N/J/C are computed with a lag of 3.
 *
 * Args:      om_fs  - optimized frameshift profile
 *            pp     - posterior decoding matrix (8-cell FS layout)
 *            ox     - RESULT: OA DP matrix (3-cell standard layout)
 *            ret_e  - RETURN: expected number of correctly decoded positions
 *
 * Returns:   eslOK on success.
 */
int
p7_OptimalAccuracy_Frameshift_SSE(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                                   P7_OMX *ox, float *ret_e)
{
  int    M   = om_fs->M;
  int    L   = pp->L;
  int    Q   = p7O_NQF(M);
  int    i, q, j;
  int    qM  = (M-1) % Q;   /* stripe containing k=M */
  int    rM  = (M-1) / Q;   /* lane  containing k=M */
  float  xN, xE, xB, xJ, xC;
  float  t1, t2;
  register __m128 mpv1, dpv1, ipv1;   /* right-shifted prev row i-1, then updated to current q */
  register __m128 mpv2, dpv2, ipv2;   /* right-shifted prev row i-2 */
  register __m128 mpv3, dpv3, ipv3;   /* right-shifted prev row i-3 */
  register __m128 mpv4, dpv4, ipv4;   /* right-shifted prev row i-4 */
  register __m128 mpv5, dpv5, ipv5;   /* right-shifted prev row i-5 */
  register __m128 sv1, sv2, sv3, sv4, sv5, sv;
  register __m128 xEv, xBv1, xBv2, xBv3, xBv4, xBv5;
  register __m128 dcv;
  __m128 bm, mm, im, dm, md, mi, ii;
  __m128 *tp;
  __m128 *dpc;                         /* current OA row (3-cell layout)     */
  __m128 *dp1, *dp2, *dp3, *dp4, *dp5; /* previous OA rows                   */
  __m128 *ppc;                         /* current pp row (8-cell FS layout)  */
  __m128  zerov = _mm_setzero_ps();
  __m128  infv  = _mm_set1_ps(-eslINFINITY);
  union { __m128 v; float p[4]; } imfix;   /* for I(i,M) = -inf fixup */

  ox->M = M;
  ox->L = L;

  /* Initialize row 0: MDI = -inf; N=0, B=0, E/J/C = -inf */
  dpc = ox->dpf[0];
  for (q = 0; q < Q; q++) { MMO(dpc,q) = DMO(dpc,q) = IMO(dpc,q) = infv; }
  ox->xmx[0*p7X_NXCELLS+p7X_E] = -eslINFINITY;
  ox->xmx[0*p7X_NXCELLS+p7X_N] =  0.0f;
  ox->xmx[0*p7X_NXCELLS+p7X_J] = -eslINFINITY;
  ox->xmx[0*p7X_NXCELLS+p7X_B] =  0.0f;
  ox->xmx[0*p7X_NXCELLS+p7X_C] = -eslINFINITY;

  for (i = 1; i <= L; i++)
    {
      dpc = ox->dpf[i];
      ppc = pp->dpf[i];

      /* Previous OA rows; row 0 (all -inf MDI) serves as sentinel for out-of-bounds lags */
      dp1 = ox->dpf[i-1];
      dp2 = (i >= 2) ? ox->dpf[i-2] : ox->dpf[0];
      dp3 = (i >= 3) ? ox->dpf[i-3] : ox->dpf[0];
      dp4 = (i >= 4) ? ox->dpf[i-4] : ox->dpf[0];
      dp5 = (i >= 5) ? ox->dpf[i-5] : ox->dpf[0];

      /* B values from previous rows; -inf sentinel when lag exceeds i */
      xBv1 = _mm_set1_ps(ox->xmx[(i-1)*p7X_NXCELLS + p7X_B]);
      xBv2 = (i >= 2) ? _mm_set1_ps(ox->xmx[(i-2)*p7X_NXCELLS + p7X_B]) : infv;
      xBv3 = (i >= 3) ? _mm_set1_ps(ox->xmx[(i-3)*p7X_NXCELLS + p7X_B]) : infv;
      xBv4 = (i >= 4) ? _mm_set1_ps(ox->xmx[(i-4)*p7X_NXCELLS + p7X_B]) : infv;
      xBv5 = (i >= 5) ? _mm_set1_ps(ox->xmx[(i-5)*p7X_NXCELLS + p7X_B]) : infv;

      /* Right-shift last stripe of each prev row for k=0 (k-1 access at q=0) */
      mpv1 = esl_sse_rightshift_ps(MMO(dp1, Q-1), infv);
      dpv1 = esl_sse_rightshift_ps(DMO(dp1, Q-1), infv);
      ipv1 = esl_sse_rightshift_ps(IMO(dp1, Q-1), infv);

      mpv2 = esl_sse_rightshift_ps(MMO(dp2, Q-1), infv);
      dpv2 = esl_sse_rightshift_ps(DMO(dp2, Q-1), infv);
      ipv2 = esl_sse_rightshift_ps(IMO(dp2, Q-1), infv);

      mpv3 = esl_sse_rightshift_ps(MMO(dp3, Q-1), infv);
      dpv3 = esl_sse_rightshift_ps(DMO(dp3, Q-1), infv);
      ipv3 = esl_sse_rightshift_ps(IMO(dp3, Q-1), infv);

      mpv4 = esl_sse_rightshift_ps(MMO(dp4, Q-1), infv);
      dpv4 = esl_sse_rightshift_ps(DMO(dp4, Q-1), infv);
      ipv4 = esl_sse_rightshift_ps(IMO(dp4, Q-1), infv);

      mpv5 = esl_sse_rightshift_ps(MMO(dp5, Q-1), infv);
      dpv5 = esl_sse_rightshift_ps(DMO(dp5, Q-1), infv);
      ipv5 = esl_sse_rightshift_ps(IMO(dp5, Q-1), infv);

      dcv = infv;
      xEv = infv;
      tp  = om_fs->tfv;

      for (q = 0; q < Q; q++)
        {
          /* Read 7 transition masks per stripe: BM, MM, IM, DM, MD, MI, II */
          bm = _mm_cmpgt_ps(*tp, zerov);  tp++;
          mm = _mm_cmpgt_ps(*tp, zerov);  tp++;
          im = _mm_cmpgt_ps(*tp, zerov);  tp++;
          dm = _mm_cmpgt_ps(*tp, zerov);  tp++;
          md = _mm_cmpgt_ps(*tp, zerov);  tp++;
          mi = _mm_cmpgt_ps(*tp, zerov);  tp++;
          ii = _mm_cmpgt_ps(*tp, zerov);  tp++;

          /* M contributions for codon lag 1 (from row i-1, k-1 access) */
          sv1 =                _mm_and_ps(bm, xBv1);
          sv1 = _mm_max_ps(sv1, _mm_and_ps(mm, mpv1));
          sv1 = _mm_max_ps(sv1, _mm_and_ps(im, ipv1));
          sv1 = _mm_max_ps(sv1, _mm_and_ps(dm, dpv1));
          sv1 = _mm_add_ps(sv1, MMO_FS(ppc, q, p7X_FS_C1));

          /* M contributions for codon lag 2 (from row i-2, k-1 access) */
          sv2 =                _mm_and_ps(bm, xBv2);
          sv2 = _mm_max_ps(sv2, _mm_and_ps(mm, mpv2));
          sv2 = _mm_max_ps(sv2, _mm_and_ps(im, ipv2));
          sv2 = _mm_max_ps(sv2, _mm_and_ps(dm, dpv2));
          sv2 = _mm_add_ps(sv2, MMO_FS(ppc, q, p7X_FS_C2));

          /* M contributions for codon lag 3 (from row i-3, k-1 access) */
          sv3 =                _mm_and_ps(bm, xBv3);
          sv3 = _mm_max_ps(sv3, _mm_and_ps(mm, mpv3));
          sv3 = _mm_max_ps(sv3, _mm_and_ps(im, ipv3));
          sv3 = _mm_max_ps(sv3, _mm_and_ps(dm, dpv3));
          sv3 = _mm_add_ps(sv3, MMO_FS(ppc, q, p7X_FS_C3));

          /* M contributions for codon lag 4 (from row i-4, k-1 access) */
          sv4 =                _mm_and_ps(bm, xBv4);
          sv4 = _mm_max_ps(sv4, _mm_and_ps(mm, mpv4));
          sv4 = _mm_max_ps(sv4, _mm_and_ps(im, ipv4));
          sv4 = _mm_max_ps(sv4, _mm_and_ps(dm, dpv4));
          sv4 = _mm_add_ps(sv4, MMO_FS(ppc, q, p7X_FS_C4));

          /* M contributions for codon lag 5 (from row i-5, k-1 access) */
          sv5 =                _mm_and_ps(bm, xBv5);
          sv5 = _mm_max_ps(sv5, _mm_and_ps(mm, mpv5));
          sv5 = _mm_max_ps(sv5, _mm_and_ps(im, ipv5));
          sv5 = _mm_max_ps(sv5, _mm_and_ps(dm, dpv5));
          sv5 = _mm_add_ps(sv5, MMO_FS(ppc, q, p7X_FS_C5));

          /* Best M(i,k) over all 5 codon lengths */
          sv = _mm_max_ps(_mm_max_ps(sv1, sv2),
                          _mm_max_ps(_mm_max_ps(sv3, sv4), sv5));

          xEv = _mm_max_ps(xEv, sv);

          /* Update prev-row ptrs from right-shifted (k-1) to current q (k access),
           * for the I-state calculation and for the next stripe's right-shift init */
          mpv1 = MMO(dp1, q);  dpv1 = DMO(dp1, q);  ipv1 = IMO(dp1, q);
          mpv2 = MMO(dp2, q);  dpv2 = DMO(dp2, q);  ipv2 = IMO(dp2, q);
          mpv3 = MMO(dp3, q);  dpv3 = DMO(dp3, q);  ipv3 = IMO(dp3, q);
          mpv4 = MMO(dp4, q);  dpv4 = DMO(dp4, q);  ipv4 = IMO(dp4, q);
          mpv5 = MMO(dp5, q);  dpv5 = DMO(dp5, q);  ipv5 = IMO(dp5, q);

          MMO(dpc, q) = sv;
          DMO(dpc, q) = dcv;
          dcv = _mm_and_ps(md, sv);

          /* I(i,k) = max(MI*M(i-3,k), II*I(i-3,k)) + pp_I(i,k); lag-3, same column */
          sv  =                _mm_and_ps(mi, mpv3);
          sv  = _mm_max_ps(sv, _mm_and_ps(ii, ipv3));
          IMO(dpc, q) = _mm_add_ps(sv, IMO_FS(ppc, q));
        }

      /* I(i,M) has no I state; MI/II transitions stored as 0 cause AND(-inf)=0.
       * The scalar explicitly sets IMX(i,M)=-inf; replicate that here.
       */
      imfix.v      = IMO(dpc, qM);
      imfix.p[rM]  = -eslINFINITY;
      IMO(dpc, qM) = imfix.v;

      /* DD propagation: first pass (includes M->D carry from main loop) */
      dcv = esl_sse_rightshift_ps(dcv, infv);
      tp  = om_fs->tfv + 7*Q;   /* DD transitions follow the 7*Q regular transitions */
      for (q = 0; q < Q; q++)
        {
          DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
          dcv         = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), DMO(dpc, q));  tp++;
        }

      /* DD propagation: 3 more lazy passes */
      for (j = 1; j < 4; j++)
        {
          dcv = esl_sse_rightshift_ps(dcv, infv);
          tp  = om_fs->tfv + 7*Q;
          for (q = 0; q < Q; q++)
            {
              DMO(dpc, q) = _mm_max_ps(dcv, DMO(dpc, q));
              dcv         = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), dcv);  tp++;
            }
        }

      /* D->E */
      for (q = 0; q < Q; q++) xEv = _mm_max_ps(xEv, DMO(dpc, q));

      /* E state */
      esl_sse_hmax_ps(xEv, &xE);
      ox->xmx[i*p7X_NXCELLS+p7X_E] = xE;

      /* N state: lag-3 for i>=3; for i<=2 only pp contributes (no prior N term) */
      if (i > 2)
        xN = (om_fs->xf[p7O_N][p7O_LOOP] == 0.0f) ? 0.0f :
             ox->xmx[(i-3)*p7X_NXCELLS+p7X_N] + pp->xmx[i*p7X_NXCELLS+p7X_N];
      else
        xN = (om_fs->xf[p7O_N][p7O_LOOP] == 0.0f) ? 0.0f :
             pp->xmx[i*p7X_NXCELLS+p7X_N];
      ox->xmx[i*p7X_NXCELLS+p7X_N] = xN;

      /* J state: lag-3 for i>=3 */
      if (i > 2) {
        t1 = (om_fs->xf[p7O_J][p7O_LOOP] == 0.0f) ? 0.0f :
             ox->xmx[(i-3)*p7X_NXCELLS+p7X_J] + pp->xmx[i*p7X_NXCELLS+p7X_J];
        t2 = (om_fs->xf[p7O_E][p7O_LOOP] == 0.0f) ? 0.0f : xE;
        xJ = ESL_MAX(t1, t2);
      } else {
        xJ = (om_fs->xf[p7O_E][p7O_LOOP] == 0.0f) ? 0.0f : xE;
      }
      ox->xmx[i*p7X_NXCELLS+p7X_J] = xJ;

      /* C state: lag-3 for i>=3 */
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

/*------------------- end, OA DP fill ---------------------------*/





/*****************************************************************
 * 2. OA traceback
 *****************************************************************/

/* get_postprob_fs:
 *   Return the posterior probability annotation for state <scur> at
 *   position (i,k), given previous state <sprv>.  The pp matrix uses
 *   the 8-cell FS layout.
 */
static inline float
get_postprob_fs(const P7_OMX *pp, int scur, int sprv, int k, int i)
{
  int   Q = p7O_NQF(pp->M);
  int   q = (k-1) % Q;
  int   r = (k-1) / Q;
  union { __m128 v; float p[4]; } u;

  switch (scur) {
  case p7T_M:
    u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C0);  return u.p[r];
  case p7T_I:
    u.v = IMO_FS(pp->dpf[i], q);              return u.p[r];
  case p7T_N: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_N]; break;
  case p7T_C: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_C]; break;
  case p7T_J: if (sprv == scur) return pp->xmx[i*p7X_NXCELLS+p7X_J]; break;
  default: break;
  }
  return 0.0f;
}

/* select_m_fs_sse:
 *   M(i,k) is reached from B(i-c), M(i-c,k-1), I(i-c,k-1), or D(i-c,k-1)
 *   for some codon length c.  To determine the STATE TYPE (independent of c),
 *   we compare accumulated OA scores at the current row i, column k-1,
 *   following the same convention as the scalar p7_OATrace_Frameshift.
 */
static inline int
select_m_fs_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m128 *tp    = om_fs->tfv + 7*q;   /* transitions to dest stripe q: BM,MM,IM,DM,... */
  union { __m128 v; float p[4]; } u, tv;
  float   path[4];
  int     state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };
  __m128  mpv, dpv, ipv;

  /* Source M/I/D at row i, column k-1 (right-shift if k-1 crosses stripe boundary) */
  if (q > 0) {
    mpv = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_M];
    dpv = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_D];
    ipv = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_I];
  } else {
    mpv = esl_sse_rightshiftz_float(ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_M]);
    dpv = esl_sse_rightshiftz_float(ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_D]);
    ipv = esl_sse_rightshiftz_float(ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_I]);
  }

  /* tp[0]=BM, tp[1]=MM, tp[2]=IM, tp[3]=DM; ordered so M beats I beats D beats B in ties */
  tv.v = *tp;  path[3] = ((tv.p[r] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_B]);  tp++;
  u.v  = mpv;  tv.v = *tp;  path[0] = ((tv.p[r] == 0.0f) ? -eslINFINITY : u.p[r]);  tp++;
  u.v  = ipv;  tv.v = *tp;  path[1] = ((tv.p[r] == 0.0f) ? -eslINFINITY : u.p[r]);  tp++;
  u.v  = dpv;  tv.v = *tp;  path[2] = ((tv.p[r] == 0.0f) ? -eslINFINITY : u.p[r]);
  return state[esl_vec_FArgMax(path, 4)];
}

/* select_d_fs_sse:
 *   D(i,k) is reached from M(i,k-1) or D(i,k-1) within the same row.
 */
static inline int
select_d_fs_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  union { __m128 v; float p[4]; } mpv, dpv, tmdv, tddv;
  float   path[2];

  if (q > 0) {
    mpv.v  = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_M];
    dpv.v  = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_D];
    tmdv.v = om_fs->tfv[7*(q-1) + p7O_MD];
    tddv.v = om_fs->tfv[7*Q     + (q-1)];   /* DD block starts at 7*Q */
  } else {
    mpv.v  = esl_sse_rightshiftz_float(ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_M]);
    dpv.v  = esl_sse_rightshiftz_float(ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_D]);
    tmdv.v = esl_sse_rightshiftz_float(om_fs->tfv[7*(Q-1) + p7O_MD]);
    tddv.v = esl_sse_rightshiftz_float(om_fs->tfv[8*Q - 1]);
  }

  path[0] = ((tmdv.p[r] == 0.0f) ? -eslINFINITY : mpv.p[r]);
  path[1] = ((tddv.p[r] == 0.0f) ? -eslINFINITY : dpv.p[r]);
  return ((path[0] >= path[1]) ? p7T_M : p7T_D);
}

/* select_i_fs_sse:
 *   I(i,k) is reached from M(i-3,k) or I(i-3,k) (lag-3, same model column).
 */
static inline int
select_i_fs_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q      = p7O_NQF(ox->M);
  int     q      = (k-1) % Q;
  int     r      = (k-1) / Q;
  int     prev_i = (i >= 3) ? i-3 : 0;
  __m128 *tp     = om_fs->tfv + 7*q + p7O_MI;
  union { __m128 v; float p[4]; } tv, mpv, ipv;
  float   path[2];

  mpv.v = ox->dpf[prev_i][q*p7X_NSCELLS + p7X_M];
  tv.v  = *tp;  path[0] = ((tv.p[r] == 0.0f) ? -eslINFINITY : mpv.p[r]);  tp++;
  ipv.v = ox->dpf[prev_i][q*p7X_NSCELLS + p7X_I];
  tv.v  = *tp;  path[1] = ((tv.p[r] == 0.0f) ? -eslINFINITY : ipv.p[r]);
  return ((path[0] >= path[1]) ? p7T_M : p7T_I);
}

/* select_n_fs_sse: N(i) comes from N(i-1) for i>0, else from S. */
static inline int
select_n_fs_sse(int i)
{
  return ((i == 0) ? p7T_S : p7T_N);
}

/* select_c_fs_sse:
 *   C(i) comes from C(i-3)+pp(i), C(i-2)+pp(i+1), C(i-1)+pp(i+2), or E(i).
 *   Returns p7T_C or p7T_E.
 */
static inline int
select_c_fs_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp, const P7_OMX *ox, int i)
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

/* select_j_fs_sse:
 *   J(i) comes from J(i)+pp(i) (loop self-transition) or E(i).
 *   Returns p7T_J or p7T_E.
 */
static inline int
select_j_fs_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp, const P7_OMX *ox, int i)
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

/* select_e_fs_sse:
 *   E(i) is reached from M(i,k) or D(i,k) for any k=1..M.
 *   Scan all stripes to find the argmax; M beats D in ties.
 */
static inline int
select_e_fs_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int *ret_k)
{
  int    Q    = p7O_NQF(ox->M);
  __m128 *dp  = ox->dpf[i];
  union { __m128 v; float p[4]; } u;
  float  max  = -eslINFINITY;
  int    smax = p7T_M, kmax = 1;
  int    q, r;

  for (q = 0; q < Q; q++) {
    u.v = *dp; dp++;    /* M */
    for (r = 0; r < 4; r++) if (u.p[r] >= max && r*Q+q+1 <= ox->M) { max = u.p[r]; smax = p7T_M; kmax = r*Q+q+1; }
    u.v = *dp; dp += 2; /* D (skip I) */
    for (r = 0; r < 4; r++) if (u.p[r] >  max && r*Q+q+1 <= ox->M) { max = u.p[r]; smax = p7T_D; kmax = r*Q+q+1; }
  }
  *ret_k = kmax;
  return smax;
}

/* select_b_fs_sse: B(i) comes from N(i) or J(i). */
static inline int
select_b_fs_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i)
{
  float path[2];
  path[0] = (om_fs->xf[p7O_N][p7O_MOVE] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_N];
  path[1] = (om_fs->xf[p7O_J][p7O_MOVE] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_J];
  return ((path[0] > path[1]) ? p7T_N : p7T_J);
}

/* select_codon_fs_sse:
 *   At M(i,k), choose the codon length c=1..5 with the highest posterior
 *   probability from the FS pp matrix.
 */
static inline int
select_codon_fs_sse(const P7_OMX *pp, int i, int k)
{
  int   Q = p7O_NQF(pp->M);
  int   q = (k-1) % Q;
  int   r = (k-1) / Q;
  union { __m128 v; float p[4]; } u;
  float codon[5];

  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C1);  codon[0] = u.p[r];
  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C2);  codon[1] = u.p[r];
  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C3);  codon[2] = u.p[r];
  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C4);  codon[3] = u.p[r];
  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C5);  codon[4] = u.p[r];
  return esl_vec_FArgMax(codon, 5) + 1;
}


/* Function:  p7_OATrace_Frameshift_SSE()
 * Synopsis:  Optimal accuracy decoding traceback, frameshift SSE version.
 *
 * Purpose:   Traceback through the OA DP matrix <ox> filled by
 *            <p7_OptimalAccuracy_Frameshift_SSE>, annotated with
 *            posterior probabilities from the FS decoding matrix <pp>.
 *
 *            The OA matrix <ox> uses standard 3-cell layout; the pp matrix
 *            uses 8-cell FS layout.  At each M state, <select_codon_fs_sse>
 *            picks the codon length c=1..5 with the highest posterior, and
 *            i is decremented by c after the step.  I states use lag-3.
 *            N/J/C special states use lag-3 nucleotide stepping.
 *
 * Args:      om_fs - optimized frameshift profile
 *            pp    - posterior decoding matrix (8-cell FS layout)
 *            ox    - OA DP matrix (3-cell standard layout, from p7_OptimalAccuracy_Frameshift_SSE)
 *            tr    - RESULT: OA traceback (must be empty; allocated with p7_trace_fs_CreateWithPP)
 *
 * Returns:   eslOK on success.
 *
 * Throws:    eslEINVAL if <tr> is not empty.
 *            eslEMEM on allocation error.
 */
int
p7_OATrace_Frameshift_SSE(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                           const P7_OMX *ox, P7_TRACE *tr)
{
  int   i    = ox->L;
  int   k    = 0;
  int   c    = 0;
  float postprob;
  int   sprv, scur;
  int   status;

  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace not empty; needs to be Reuse()'d?");

  if ((status = p7_trace_fs_Append      (tr, p7T_T, k, i, c))           != eslOK) return status;
  if ((status = p7_trace_fs_Append      (tr, p7T_C, k, i, c))           != eslOK) return status;

  sprv = p7T_C;
  while (sprv != p7T_S)
    {
      switch (sprv) {
      case p7T_M: scur = select_m_fs_sse(om_fs,     ox, i, k);       k--;   break;
      case p7T_D: scur = select_d_fs_sse(om_fs,     ox, i, k);       k--;   break;
      case p7T_I: scur = select_i_fs_sse(om_fs,     ox, i, k);       i -= 3; break;
      case p7T_N: scur = select_n_fs_sse(i);                                  break;
      case p7T_C: scur = select_c_fs_sse(om_fs, pp, ox, i);                   break;
      case p7T_J: scur = select_j_fs_sse(om_fs, pp, ox, i);                   break;
      case p7T_E: scur = select_e_fs_sse(om_fs,     ox, i, &k);               break;
      case p7T_B: scur = select_b_fs_sse(om_fs,     ox, i);                   break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in OA FS traceback");
      }
      if (scur == -1) ESL_EXCEPTION(eslEINVAL, "OA FS SSE traceback choice failed");

      postprob = get_postprob_fs(pp, scur, sprv, k, i);

      if (scur == p7T_M) c = select_codon_fs_sse(pp, i, k);
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

/*---------------------- end, OA traceback ----------------------*/




/*****************************************************************
 * 3. Benchmark driver
 *****************************************************************/



/*---------------- end, benchmark driver ------------------------*/




/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7OPTACC_FS_TESTDRIVE
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_optacc_fs:
 *   Compare p7_OptimalAccuracy_Frameshift_SSE and p7_OATrace_Frameshift_SSE
 *   against their scalar counterparts.
 *
 *   Strategy: run the full SSE pipeline (Forward->Backward->Decoding->OA->Trace)
 *   and the scalar pipeline (via deconverted pp matrix).  Compare OA scores,
 *   OA DP matrices, and tracebacks for equality within tolerance.
 *
 *   Note: Like utest_optacc in optacc.c, this test can fail for innocent
 *   reasons when numerical ties in the OA matrix lead to different traceback
 *   choices between SSE and scalar implementations.
 */
static void
utest_optacc_fs(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA,
                ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode,
                P7_BG *bgAA, P7_BG *bgDNA, int M, int L, int N)
{
  char            *msg   = "optacc_fs SSE unit test failed";
  P7_HMM          *hmm   = NULL;
  P7_FS_PROFILE   *gm_fs = NULL;
  P7_FS_OPROFILE  *om_fs = NULL;
  P7_GMX          *gx_pp  = NULL;   /* scalar pp (8-cell FS) */
  P7_GMX          *gx_oa  = NULL;   /* scalar OA (3-cell)    */
  P7_GMX          *gx_oa2 = NULL;   /* SSE OA deconverted    */
  P7_OMX          *ox_fwd = NULL;   /* SSE forward/pp (8-cell FS) */
  P7_OMX          *ox_bck = NULL;   /* SSE backward (8-cell FS)   */
  P7_TRACE        *tr_sse = NULL;   /* SSE OA traceback           */
  P7_TRACE        *tr_ref = NULL;   /* scalar OA traceback        */
  ESL_DSQ         *dsq    = NULL;
  float            fsc, bsc, accscore_ref, accscore_sse;
  float            tol    = 0.01f;
  float            pptol  = 0.01f;
  int              seq_L;

  p7_FLogsumInit();

  if (p7_hmm_Sample(r, M, abcAA, &hmm)                              != eslOK) esl_fatal(msg);

  gm_fs = p7_profile_fs_Create(hmm->M, abcAA, p7P_5CODONS);
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs, L, p7_UNILOCAL) != eslOK) esl_fatal(msg);

  om_fs = p7_fs_oprofile_Create(hmm->M, abcAA, p7P_5CODONS);
  if (p7_fs_oprofile_Convert(gm_fs, om_fs)                          != eslOK) esl_fatal(msg);

  seq_L = L;
  dsq   = malloc(sizeof(ESL_DSQ) * (seq_L + 2));

  /* SSE: forward/backward use 8-cell FS layout; OA uses standard 3-cell */
  ox_fwd = p7_omx_Create_dpf(hmm->M, seq_L, seq_L, p7X_NSCELLS_FS);
  ox_bck = p7_omx_Create_dpf(hmm->M, seq_L, seq_L, p7X_NSCELLS);

  /* Scalar: pp uses 8-cell FS layout; OA matrices use standard 3-cell */
  gx_pp  = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS_FS);
  gx_oa  = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS);
  gx_oa2 = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS);

  tr_sse = p7_trace_fs_CreateWithPP();
  tr_ref = p7_trace_fs_CreateWithPP();

  while (N--)
    {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, seq_L, dsq);
      printf("N %d\n", N);
      /* SSE pipeline: Forward -> Backward -> Decoding (in-place into ox_fwd) -> OA -> Trace */
      if (p7_Forward_Frameshift_SSE (dsq, seq_L, om_fs,         ox_fwd, &fsc) != eslOK) esl_fatal(msg);
      if (p7_Backward_Frameshift_SSE(dsq, seq_L, om_fs, ox_fwd, ox_bck, &bsc) != eslOK) esl_fatal(msg);
      if (esl_FCompare_old(fsc, bsc, tol)                                      != eslOK) esl_fatal(msg);
      if (p7_Decoding_Frameshift_SSE(om_fs, ox_fwd, ox_bck)                    != eslOK) esl_fatal(msg);
      /* ox_fwd now holds the SSE posterior probability matrix */
      if (p7_OptimalAccuracy_Frameshift_SSE(om_fs, ox_fwd, ox_bck, &accscore_sse) != eslOK) esl_fatal(msg);
      if (p7_OATrace_Frameshift_SSE(om_fs, ox_fwd, ox_bck, tr_sse)              != eslOK) esl_fatal(msg);
//       p7_trace_fs_Dump(stdout, tr_sse, gm_fs, dsq, abcDNA);
      /* Reference: convert SSE pp to scalar GMX, run scalar OA and trace */
      if (p7_omx_FDeconvert(ox_fwd, gx_pp)                                    != eslOK) esl_fatal(msg);
      if (p7_OptimalAccuracy_Frameshift(gm_fs, gx_pp, gx_oa, &accscore_ref)   != eslOK) esl_fatal(msg);
      if (p7_OATrace_Frameshift(gm_fs, gx_pp, gx_oa, tr_ref)                  != eslOK) esl_fatal(msg);
//      p7_trace_fs_Dump(stdout, tr_ref, gm_fs, dsq, abcDNA);
      /* Convert SSE OA to scalar GMX for matrix-level comparison */
      if (p7_omx_FDeconvert(ox_bck, gx_oa2)                                    != eslOK) esl_fatal(msg);
      
      printf("111111\n");
      if (esl_FCompare_old(accscore_sse, accscore_ref, tol) != eslOK) esl_fatal(msg);
      printf("2222222\n");
      p7_gmx_Dump(stdout, gx_oa, p7_DEFAULT);
      p7_gmx_Dump(stdout, gx_oa2, p7_DEFAULT);

      if (p7_gmx_Compare(gx_oa, gx_oa2, tol)               != eslOK) esl_fatal(msg);
      printf("3333333\n");
      if (p7_trace_Compare(tr_sse, tr_ref, pptol)           != eslOK) esl_fatal(msg);

      if (esl_opt_GetBoolean(go, "--traces"))
        {
          p7_trace_fs_Dump(stdout, tr_sse, gm_fs, dsq, abcDNA);
          p7_trace_fs_Dump(stdout, tr_ref, gm_fs, dsq, abcDNA);
        }

      p7_trace_Reuse(tr_sse);
      p7_trace_Reuse(tr_ref);
    }

  free(dsq);
  p7_trace_fs_Destroy(tr_sse);
  p7_trace_fs_Destroy(tr_ref);
  p7_omx_Destroy(ox_fwd);
  p7_omx_Destroy(ox_bck);
  p7_gmx_Destroy(gx_pp);
  p7_gmx_Destroy(gx_oa);
  p7_gmx_Destroy(gx_oa2);
  p7_fs_oprofile_Destroy(om_fs);
  p7_profile_fs_Destroy(gm_fs);
  p7_hmm_Destroy(hmm);
}
#endif /*p7OPTACC_FS_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/




/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7OPTACC_FS_TESTDRIVE
/*
   gcc -g -Wall -msse2 -std=gnu99 -o optacc_fs_utest \
       -I.. -L.. -I../../easel -L../../easel \
       -Dp7OPTACC_FS_TESTDRIVE optacc_fs.c -lhmmer -leasel -lm
   ./optacc_fs_utest
*/
#include "p7_config.h"
#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_getopts.h"
#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  { "-h",       eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help",              0 },
  { "-s",       eslARG_INT,    "42", NULL, NULL, NULL, NULL, NULL, "set random number seed",       0 },
  { "-L",       eslARG_INT,   "150", NULL, NULL, NULL, NULL, NULL, "length of random target seqs", 0 },
  { "-M",       eslARG_INT,    "45", NULL, NULL, NULL, NULL, NULL, "size of random models",        0 },
  { "-N",       eslARG_INT,    "20", NULL, NULL, NULL, NULL, NULL, "number of test sequences",     0 },
  { "--traces", eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "dump all tracebacks",          0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE frameshift optimal accuracy";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA  = esl_alphabet_Create(eslAMINO);
  ESL_ALPHABET   *abcDNA = esl_alphabet_Create(eslDNA);
  ESL_GENCODE    *gcode  = esl_gencode_Create(abcDNA, abcAA);
  P7_BG          *bgAA   = p7_bg_Create(abcAA);
  P7_BG          *bgDNA  = p7_bg_Create(abcDNA);
  int             M      = esl_opt_GetInteger(go, "-M");
  int             L      = esl_opt_GetInteger(go, "-L");
  int             N      = esl_opt_GetInteger(go, "-N");

  utest_optacc_fs(go, r, abcAA, abcDNA, gcode, bgAA, bgDNA, M, L, N);

  p7_bg_Destroy(bgAA);
  p7_bg_Destroy(bgDNA);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcAA);
  esl_alphabet_Destroy(abcDNA);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return eslOK;
}
#endif /*p7OPTACC_FS_TESTDRIVE*/
/*------------------ end, test driver ---------------------------*/




