/* Optimal accuracy frameshift alignment; NEON version.
 *
 * Contents:
 *   1. Optimal accuracy alignment, DP fill
 *   2. OA traceback
 *   3. Benchmark driver
 *   4. Unit tests
 *   5. Test driver
 */
#include "p7_config.h"

#include <float.h>

#include <arm_neon.h>

#include "easel.h"
#include "esl_neon.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_neon.h"


/*****************************************************************
 * 1. Optimal accuracy alignment, DP fill
 *****************************************************************/

/* Function:  p7_OptimalAccuracy_Frameshift()
 * Synopsis:  DP fill of optimal accuracy alignment, frameshift NEON version.
 *
 * Purpose:   Calculates the fill step of the optimal accuracy decoding
 *            algorithm for a frameshift-aware profile, using NEON SIMD
 *            parallelism over model positions.
 *
 *            The posterior decoding matrix <pp> must have the 8-cell FS
 *            layout (p7X_NSCELLS_FS=8; from p7_Decoding_Frameshift).
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
p7_OptimalAccuracy_Frameshift(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                                   P7_OMX *ox, float *ret_e)
{
  int    M   = om_fs->M;
  int    L   = pp->L;
  int    Q   = p7O_NQF(M);
  int    i, q, j;
  int    qM  = (M-1) % Q;   /* stripe containing k=M */
  int    rM  = (M-1) / Q;   /* lane  containing k=M */
  float  xN, xE, xJ, xC;
  float  t1, t2;
  register float32x4_t mpv1, dpv1, ipv1;   /* right-shifted prev row i-1, then updated to current q */
  register float32x4_t mpv2, dpv2, ipv2;   /* right-shifted prev row i-2 */
  register float32x4_t mpv3, dpv3, ipv3;   /* right-shifted prev row i-3 */
  register float32x4_t mpv4, dpv4, ipv4;   /* right-shifted prev row i-4 */
  register float32x4_t mpv5, dpv5, ipv5;   /* right-shifted prev row i-5 */
  register float32x4_t sv1, sv2, sv3, sv4, sv5, sv;
  register float32x4_t xEv, xBv1, xBv2, xBv3, xBv4, xBv5;
  register float32x4_t dcv;
  float32x4_t bm, mm, im, dm, md, mi, ii;
  float32x4_t *tp;
  float32x4_t *dpc;                           /* current OA row (3-cell layout)     */
  float32x4_t *dp1, *dp2, *dp3, *dp4, *dp5;  /* previous OA rows                   */
  float32x4_t *ppc;                           /* current pp row (8-cell FS layout)  */
  float32x4_t  zerov = vdupq_n_f32(0.0f);
  float32x4_t  infv  = vdupq_n_f32(-eslINFINITY);
  union { float32x4_t v; float p[4]; } imfix;   /* for I(i,M) = -inf fixup */

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
      xBv1 = vdupq_n_f32(ox->xmx[(i-1)*p7X_NXCELLS + p7X_B]);
      xBv2 = (i >= 2) ? vdupq_n_f32(ox->xmx[(i-2)*p7X_NXCELLS + p7X_B]) : infv;
      xBv3 = (i >= 3) ? vdupq_n_f32(ox->xmx[(i-3)*p7X_NXCELLS + p7X_B]) : infv;
      xBv4 = (i >= 4) ? vdupq_n_f32(ox->xmx[(i-4)*p7X_NXCELLS + p7X_B]) : infv;
      xBv5 = (i >= 5) ? vdupq_n_f32(ox->xmx[(i-5)*p7X_NXCELLS + p7X_B]) : infv;

      /* Right-shift last stripe of each prev row for k=0 (k-1 access at q=0) */
      mpv1 = vextq_f32(infv, MMO(dp1, Q-1), 3);
      dpv1 = vextq_f32(infv, DMO(dp1, Q-1), 3);
      ipv1 = vextq_f32(infv, IMO(dp1, Q-1), 3);

      mpv2 = vextq_f32(infv, MMO(dp2, Q-1), 3);
      dpv2 = vextq_f32(infv, DMO(dp2, Q-1), 3);
      ipv2 = vextq_f32(infv, IMO(dp2, Q-1), 3);

      mpv3 = vextq_f32(infv, MMO(dp3, Q-1), 3);
      dpv3 = vextq_f32(infv, DMO(dp3, Q-1), 3);
      ipv3 = vextq_f32(infv, IMO(dp3, Q-1), 3);

      mpv4 = vextq_f32(infv, MMO(dp4, Q-1), 3);
      dpv4 = vextq_f32(infv, DMO(dp4, Q-1), 3);
      ipv4 = vextq_f32(infv, IMO(dp4, Q-1), 3);

      mpv5 = vextq_f32(infv, MMO(dp5, Q-1), 3);
      dpv5 = vextq_f32(infv, DMO(dp5, Q-1), 3);
      ipv5 = vextq_f32(infv, IMO(dp5, Q-1), 3);

      dcv = infv;
      xEv = infv;
      tp  = om_fs->tfv;

      for (q = 0; q < Q; q++)
        {
          /* Read 7 transition masks per stripe: BM, MM, IM, DM, MD, MI, II
           * Mask: lane is all-1 if transition > 0, else all-0 */
          bm = vreinterpretq_f32_u32(vcgtq_f32(*tp, zerov));  tp++;
          mm = vreinterpretq_f32_u32(vcgtq_f32(*tp, zerov));  tp++;
          im = vreinterpretq_f32_u32(vcgtq_f32(*tp, zerov));  tp++;
          dm = vreinterpretq_f32_u32(vcgtq_f32(*tp, zerov));  tp++;
          md = vreinterpretq_f32_u32(vcgtq_f32(*tp, zerov));  tp++;
          mi = vreinterpretq_f32_u32(vcgtq_f32(*tp, zerov));  tp++;
          ii = vreinterpretq_f32_u32(vcgtq_f32(*tp, zerov));  tp++;

          /* M contributions for codon lag 1 (from row i-1, k-1 access) */
          sv1 =                vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(bm), vreinterpretq_u32_f32(xBv1)));
          sv1 = vmaxq_f32(sv1, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(mm), vreinterpretq_u32_f32(mpv1))));
          sv1 = vmaxq_f32(sv1, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(im), vreinterpretq_u32_f32(ipv1))));
          sv1 = vmaxq_f32(sv1, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(dm), vreinterpretq_u32_f32(dpv1))));
          sv1 = vaddq_f32(sv1, MMO_FS(ppc, q, p7X_FS_C1));

          /* M contributions for codon lag 2 (from row i-2, k-1 access) */
          sv2 =                vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(bm), vreinterpretq_u32_f32(xBv2)));
          sv2 = vmaxq_f32(sv2, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(mm), vreinterpretq_u32_f32(mpv2))));
          sv2 = vmaxq_f32(sv2, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(im), vreinterpretq_u32_f32(ipv2))));
          sv2 = vmaxq_f32(sv2, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(dm), vreinterpretq_u32_f32(dpv2))));
          sv2 = vaddq_f32(sv2, MMO_FS(ppc, q, p7X_FS_C2));

          /* M contributions for codon lag 3 (from row i-3, k-1 access) */
          sv3 =                vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(bm), vreinterpretq_u32_f32(xBv3)));
          sv3 = vmaxq_f32(sv3, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(mm), vreinterpretq_u32_f32(mpv3))));
          sv3 = vmaxq_f32(sv3, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(im), vreinterpretq_u32_f32(ipv3))));
          sv3 = vmaxq_f32(sv3, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(dm), vreinterpretq_u32_f32(dpv3))));
          sv3 = vaddq_f32(sv3, MMO_FS(ppc, q, p7X_FS_C3));

          /* M contributions for codon lag 4 (from row i-4, k-1 access) */
          sv4 =                vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(bm), vreinterpretq_u32_f32(xBv4)));
          sv4 = vmaxq_f32(sv4, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(mm), vreinterpretq_u32_f32(mpv4))));
          sv4 = vmaxq_f32(sv4, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(im), vreinterpretq_u32_f32(ipv4))));
          sv4 = vmaxq_f32(sv4, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(dm), vreinterpretq_u32_f32(dpv4))));
          sv4 = vaddq_f32(sv4, MMO_FS(ppc, q, p7X_FS_C4));

          /* M contributions for codon lag 5 (from row i-5, k-1 access) */
          sv5 =                vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(bm), vreinterpretq_u32_f32(xBv5)));
          sv5 = vmaxq_f32(sv5, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(mm), vreinterpretq_u32_f32(mpv5))));
          sv5 = vmaxq_f32(sv5, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(im), vreinterpretq_u32_f32(ipv5))));
          sv5 = vmaxq_f32(sv5, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(dm), vreinterpretq_u32_f32(dpv5))));
          sv5 = vaddq_f32(sv5, MMO_FS(ppc, q, p7X_FS_C5));

          /* Best M(i,k) over all 5 codon lengths */
          sv = vmaxq_f32(vmaxq_f32(sv1, sv2),
                         vmaxq_f32(vmaxq_f32(sv3, sv4), sv5));

          xEv = vmaxq_f32(xEv, sv);

          /* Update prev-row ptrs from right-shifted (k-1) to current q (k access),
           * for the I-state calculation and for the next stripe's right-shift init */
          mpv1 = MMO(dp1, q);  dpv1 = DMO(dp1, q);  ipv1 = IMO(dp1, q);
          mpv2 = MMO(dp2, q);  dpv2 = DMO(dp2, q);  ipv2 = IMO(dp2, q);
          mpv3 = MMO(dp3, q);  dpv3 = DMO(dp3, q);  ipv3 = IMO(dp3, q);
          mpv4 = MMO(dp4, q);  dpv4 = DMO(dp4, q);  ipv4 = IMO(dp4, q);
          mpv5 = MMO(dp5, q);  dpv5 = DMO(dp5, q);  ipv5 = IMO(dp5, q);

          MMO(dpc, q) = sv;
          DMO(dpc, q) = dcv;
          dcv = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(md), vreinterpretq_u32_f32(sv)));

          /* I(i,k) = max(MI*M(i-3,k), II*I(i-3,k)) + pp_I(i,k); lag-3, same column */
          sv  =                vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(mi), vreinterpretq_u32_f32(mpv3)));
          sv  = vmaxq_f32(sv, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(ii), vreinterpretq_u32_f32(ipv3))));
          IMO(dpc, q) = vaddq_f32(sv, IMO_FS(ppc, q));
        }

      /* I(i,M) has no I state; MI/II transitions stored as 0 cause AND(0)=0.
       * The scalar explicitly sets IMX(i,M)=-inf; replicate that here.
       */
      imfix.v      = IMO(dpc, qM);
      imfix.p[rM]  = -eslINFINITY;
      IMO(dpc, qM) = imfix.v;

      /* DD propagation: first pass (includes M->D carry from main loop) */
      dcv = vextq_f32(infv, dcv, 3);
      tp  = om_fs->tfv + 7*Q;   /* DD transitions follow the 7*Q regular transitions */
      for (q = 0; q < Q; q++)
        {
          DMO(dpc, q) = vmaxq_f32(dcv, DMO(dpc, q));
          dcv         = vreinterpretq_f32_u32(vandq_u32(
                          vreinterpretq_u32_f32(vreinterpretq_f32_u32(vcgtq_f32(*tp, zerov))),
                          vreinterpretq_u32_f32(DMO(dpc, q))));  tp++;
        }

      /* DD propagation: 3 more lazy passes */
      for (j = 1; j < 4; j++)
        {
          dcv = vextq_f32(infv, dcv, 3);
          tp  = om_fs->tfv + 7*Q;
          for (q = 0; q < Q; q++)
            {
              DMO(dpc, q) = vmaxq_f32(dcv, DMO(dpc, q));
              dcv         = vreinterpretq_f32_u32(vandq_u32(
                              vreinterpretq_u32_f32(vreinterpretq_f32_u32(vcgtq_f32(*tp, zerov))),
                              vreinterpretq_u32_f32(dcv)));  tp++;
            }
        }

      /* D->E */
      for (q = 0; q < Q; q++) xEv = vmaxq_f32(xEv, DMO(dpc, q));

      /* E state */
      xE = esl_neon_hmax_f32((esl_neon_128f_t){ .f32x4 = xEv });
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
  union { float32x4_t v; float p[4]; } u;

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

/* select_m_fs:
 *   M(i,k) is reached from B(i-c), M(i-c,k-1), I(i-c,k-1), or D(i-c,k-1)
 *   for some codon length c.  To determine the STATE TYPE (independent of c),
 *   we compare accumulated OA scores at the current row i, column k-1,
 *   following the same convention as the scalar p7_OATrace_Frameshift.
 */
static inline int
select_m_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  float32x4_t *tp = om_fs->tfv + 7*q;   /* transitions to dest stripe q: BM,MM,IM,DM,... */
  union { float32x4_t v; float p[4]; } u, tv;
  float   path[4];
  int     state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };
  float32x4_t mpv, dpv, ipv;

  /* Source M/I/D at row i, column k-1 (right-shift if k-1 crosses stripe boundary).
   * Use -inf fill (not zero) so k=0 boundary values behave like MMX(i,0)=-inf (scalar). */
  if (q > 0) {
    mpv = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_M];
    dpv = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_D];
    ipv = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_I];
  } else {
    float32x4_t neg_inf = vdupq_n_f32(-eslINFINITY);
    mpv = vextq_f32(neg_inf, ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_M], 3);
    dpv = vextq_f32(neg_inf, ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_D], 3);
    ipv = vextq_f32(neg_inf, ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_I], 3);
  }

  /* tp[0]=BM, tp[1]=MM, tp[2]=IM, tp[3]=DM; ordered so M beats I beats D beats B in ties */
  tv.v = *tp;  path[3] = ((tv.p[r] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_B]);  tp++;
  u.v  = mpv;  tv.v = *tp;  path[0] = ((tv.p[r] == 0.0f) ? -eslINFINITY : u.p[r]);  tp++;
  u.v  = ipv;  tv.v = *tp;  path[1] = ((tv.p[r] == 0.0f) ? -eslINFINITY : u.p[r]);  tp++;
  u.v  = dpv;  tv.v = *tp;  path[2] = ((tv.p[r] == 0.0f) ? -eslINFINITY : u.p[r]);
  return state[esl_vec_FArgMax(path, 4)];
}

/* select_d_fs:
 *   D(i,k) is reached from M(i,k-1) or D(i,k-1) within the same row.
 */
static inline int
select_d_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  union { float32x4_t v; float p[4]; } mpv, dpv, tmdv, tddv;
  float   path[2];

  if (q > 0) {
    mpv.v  = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_M];
    dpv.v  = ox->dpf[i][(q-1)*p7X_NSCELLS + p7X_D];
    tmdv.v = om_fs->tfv[7*(q-1) + p7O_MD];
    tddv.v = om_fs->tfv[7*Q     + (q-1)];   /* DD block starts at 7*Q */
  } else {
    /* Use -inf fill so k=0 boundary gives -inf (matching scalar DMX(i,0)=-inf). */
    float32x4_t neg_inf = vdupq_n_f32(-eslINFINITY);
    float32x4_t zerov   = vdupq_n_f32(0.0f);
    mpv.v  = vextq_f32(neg_inf, ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_M], 3);
    dpv.v  = vextq_f32(neg_inf, ox->dpf[i][(Q-1)*p7X_NSCELLS + p7X_D], 3);
    tmdv.v = vextq_f32(zerov,   om_fs->tfv[7*(Q-1) + p7O_MD], 3);
    tddv.v = vextq_f32(zerov,   om_fs->tfv[8*Q - 1], 3);
  }

  path[0] = ((tmdv.p[r] == 0.0f) ? -eslINFINITY : mpv.p[r]);
  path[1] = ((tddv.p[r] == 0.0f) ? -eslINFINITY : dpv.p[r]);
  return ((path[0] >= path[1]) ? p7T_M : p7T_D);
}

/* select_i_fs:
 *   I(i,k) is reached from M(i-3,k) or I(i-3,k) (lag-3, same model column).
 */
static inline int
select_i_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q      = p7O_NQF(ox->M);
  int     q      = (k-1) % Q;
  int     r      = (k-1) / Q;
  int     prev_i = (i >= 3) ? i-3 : 0;
  float32x4_t *tp = om_fs->tfv + 7*q + p7O_MI;
  union { float32x4_t v; float p[4]; } tv, mpv, ipv;
  float   path[2];

  mpv.v = ox->dpf[prev_i][q*p7X_NSCELLS + p7X_M];
  tv.v  = *tp;  path[0] = ((tv.p[r] == 0.0f) ? -eslINFINITY : mpv.p[r]);  tp++;
  ipv.v = ox->dpf[prev_i][q*p7X_NSCELLS + p7X_I];
  tv.v  = *tp;  path[1] = ((tv.p[r] == 0.0f) ? -eslINFINITY : ipv.p[r]);
  return ((path[0] >= path[1]) ? p7T_M : p7T_I);
}

/* select_n_fs: N(i) comes from N(i-1) for i>0, else from S. */
static inline int
select_n_fs(int i)
{
  return ((i == 0) ? p7T_S : p7T_N);
}

/* select_c_fs:
 *   C(i) comes from C(i-3)+pp(i), C(i-2)+pp(i+1), C(i-1)+pp(i+2), or E(i).
 *   Returns p7T_C or p7T_E.
 */
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

/* select_j_fs:
 *   J(i) comes from J(i)+pp(i) (loop self-transition) or E(i).
 *   Returns p7T_J or p7T_E.
 */
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

/* select_e_fs:
 *   E(i) is reached from M(i,k) or D(i,k) for any k=1..M.
 *   Iterate in sequential k=1..M order with strict >, matching the scalar
 *   select_e() tie-breaking (first/lowest k with the max value wins).
 */
static inline int
select_e_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int *ret_k)
{
  int    Q    = p7O_NQF(ox->M);
  union { float32x4_t v; float p[4]; } u;
  float  max  = -eslINFINITY;
  int    smax = p7T_M, kmax = 1;
  int    k, q, r;

  for (k = 1; k <= ox->M; k++) {
    q = (k-1) % Q;
    r = (k-1) / Q;
    u.v = ox->dpf[i][q * p7X_NSCELLS + p7X_M];
    if (u.p[r] > max) { max = u.p[r]; smax = p7T_M; kmax = k; }
    u.v = ox->dpf[i][q * p7X_NSCELLS + p7X_D];
    if (u.p[r] > max) { max = u.p[r]; smax = p7T_D; kmax = k; }
  }
  *ret_k = kmax;
  return smax;
}

/* select_b_fs: B(i) comes from N(i) or J(i). */
static inline int
select_b_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i)
{
  float path[2];
  path[0] = (om_fs->xf[p7O_N][p7O_MOVE] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_N];
  path[1] = (om_fs->xf[p7O_J][p7O_MOVE] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_J];
  return ((path[0] > path[1]) ? p7T_N : p7T_J);
}

/* select_codon_fs:
 *   At M(i,k), choose the codon length c=1..5 with the highest posterior
 *   probability from the FS pp matrix.
 */
static inline int
select_codon_fs(const P7_OMX *pp, int i, int k)
{
  int   Q = p7O_NQF(pp->M);
  int   q = (k-1) % Q;
  int   r = (k-1) / Q;
  union { float32x4_t v; float p[4]; } u;
  float codon[5];

  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C1);  codon[0] = u.p[r];
  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C2);  codon[1] = u.p[r];
  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C3);  codon[2] = u.p[r];
  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C4);  codon[3] = u.p[r];
  u.v = MMO_FS(pp->dpf[i], q, p7X_FS_C5);  codon[4] = u.p[r];
  return esl_vec_FArgMax(codon, 5) + 1;
}


/* Function:  p7_OATrace_Frameshift()
 * Synopsis:  Optimal accuracy decoding traceback, frameshift NEON version.
 *
 * Purpose:   Traceback through the OA DP matrix <ox> filled by
 *            <p7_OptimalAccuracy_Frameshift>, annotated with
 *            posterior probabilities from the FS decoding matrix <pp>.
 *
 *            The OA matrix <ox> uses standard 3-cell layout; the pp matrix
 *            uses 8-cell FS layout.  At each M state, <select_codon_fs>
 *            picks the codon length c=1..5 that maximizes the joint OA
 *            predecessor score plus posterior pp_Cc(i,k), and i is
 *            decremented by c after the step.  I states use lag-3.
 *            N/J/C special states use lag-3 nucleotide stepping.
 *
 * Args:      om_fs - optimized frameshift profile
 *            pp    - posterior decoding matrix (8-cell FS layout)
 *            ox    - OA DP matrix (3-cell standard layout, from p7_OptimalAccuracy_Frameshift)
 *            tr    - RESULT: OA traceback (must be empty; allocated with p7_trace_fs_CreateWithPP)
 *
 * Returns:   eslOK on success.
 *
 * Throws:    eslEINVAL if <tr> is not empty.
 *            eslEMEM on allocation error.
 */
int
p7_OATrace_Frameshift(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                           const P7_OMX *ox, P7_TRACE *tr)
{
  int   i    = ox->L;
  int   k    = 0;
  int   c    = 0;
  float postprob;
  int   sprv, scur;
  int   status;

  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace not empty; needs to be Reuse()'d?");

  if ((status = p7_trace_fs_AppendWithPP(tr, p7T_T, k, i, c, 0.0f))     != eslOK) return status;
  if ((status = p7_trace_fs_AppendWithPP(tr, p7T_C, k, i, c, 0.0f))     != eslOK) return status;

  sprv = p7T_C;
  while (sprv != p7T_S)
    {
      switch (sprv) {
      case p7T_M: scur = select_m_fs(om_fs,     ox, i, k);       k--;   break;
      case p7T_D: scur = select_d_fs(om_fs,     ox, i, k);       k--;   break;
      case p7T_I: scur = select_i_fs(om_fs,     ox, i, k);       i -= 3; break;
      case p7T_N: scur = select_n_fs(i);                                  break;
      case p7T_C: scur = select_c_fs(om_fs, pp, ox, i);                   break;
      case p7T_J: scur = select_j_fs(om_fs, pp, ox, i);                   break;
      case p7T_E: scur = select_e_fs(om_fs,     ox, i, &k);               break;
      case p7T_B: scur = select_b_fs(om_fs,     ox, i);                   break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in OA FS traceback");
      }
      if (scur == -1) ESL_EXCEPTION(eslEINVAL, "OA FS NEON traceback choice failed");

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

/*---------------------- end, OA traceback ----------------------*/


/*****************************************************************
 * 3. Benchmark driver
 *****************************************************************/
#ifdef p7OPTACC_FS_BENCHMARK
/*
   gcc -g -O3 -march=armv8-a -std=gnu99 -o optacc_fs_benchmark \
       -I.. -L.. -I../../easel -L../../easel \
       -Dp7OPTACC_FS_BENCHMARK optacc_fs.c -lhmmer -leasel -lm
   ./optacc_fs_benchmark <hmmfile>
   ./optacc_fs_benchmark --notrace <hmmfile>   (benchmark DP fill only)

   Fill pipeline (Forward->Backward->Decoding) is run once to produce a
   posterior probability matrix before the timed loop.  p7_Decoding_Frameshift()
   overwrites the forward matrix in-place with the pp matrix; because pp is
   read-only in p7_OptimalAccuracy_Frameshift() and p7_OATrace_Frameshift(),
   the timed loop measures only OA fill and trace cost.
 */
#include "p7_config.h"

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
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-L",        eslARG_INT,   "1200", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs (nucleotides)",       0 },
  { "-N",        eslARG_INT,   "5000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  { "--notrace", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark the DP fill stage",                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for optimal accuracy alignment, NEON frameshift version";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA   = NULL;
  ESL_ALPHABET   *abcDNA  = NULL;
  ESL_GENCODE    *gcode   = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bgAA    = NULL;
  P7_BG          *bgDNA   = NULL;
  P7_FS_PROFILE  *gm_fs5  = NULL;
  P7_FS_OPROFILE *om_fs5  = NULL;
  P7_OMX         *ox_fwd  = NULL;   /* forward matrix; overwritten in-place by Decoding -> pp */
  P7_OMX         *ox_bck  = NULL;   /* backward matrix; overwritten by OptimalAccuracy -> OA  */
  P7_TRACE       *tr      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           fsc, bsc, accscore;
  double          Mcs;

  p7_FLogsumInit();

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)          != eslOK) p7_Fail("Failed to read HMM");

  abcDNA = esl_alphabet_Create(eslDNA);
  gcode  = esl_gencode_Create(abcDNA, abcAA);
  bgAA   = p7_bg_Create(abcAA);                 p7_bg_SetLength(bgAA, L/3);
  bgDNA  = p7_bg_Create(abcDNA);                p7_bg_SetLength(bgDNA, L);

  gm_fs5 = p7_profile_fs_Create(hmm->M, abcAA, p7P_5CODONS);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_UNILOCAL);

  om_fs5 = p7_fs_oprofile_Create(hmm->M, abcAA, p7P_5CODONS);
  p7_fs_oprofile_Convert(gm_fs5, om_fs5);
  p7_fs_oprofile_ReconfigLength(om_fs5, L/3);

  /* ox_fwd: 8-cell FS layout (Forward input; becomes pp after Decoding)  */
  /* ox_bck: standard 3-cell layout (Backward input; becomes OA after OA) */
  ox_fwd = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS_FS);
  ox_bck = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS);
  P7_OIVX *ov5 = p7_oivx_Create(hmm->M, p7P_5CODONS);
  tr     = p7_trace_fs_CreateWithPP();

  /* Fill pipeline once to produce pp matrix in ox_fwd */
  esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
  p7_Forward_Frameshift (dsq, L, om_fs5, ox_fwd, ov5, &fsc);
  p7_Backward_Frameshift(dsq, L, om_fs5, ox_fwd, ox_bck, ov5, &bsc);
  p7_Decoding_Frameshift(om_fs5, ox_fwd, ox_bck);  /* ox_fwd now holds pp */

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      p7_OptimalAccuracy_Frameshift(om_fs5, ox_fwd, ox_bck, &accscore);

      if (! esl_opt_GetBoolean(go, "--notrace"))
        {
          p7_OATrace_Frameshift(om_fs5, ox_fwd, ox_bck, tr);
          p7_trace_Reuse(tr);
        }
    }
  esl_stopwatch_Stop(w);

  Mcs = (double) N * (double) L * (double) hmm->M * 1e-6 / w->user;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   hmm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_trace_fs_Destroy(tr);
  p7_omx_Destroy(ox_bck);
  p7_omx_Destroy(ox_fwd);
  p7_oivx_Destroy(ov5);
  p7_fs_oprofile_Destroy(om_fs5);
  p7_profile_fs_Destroy(gm_fs5);
  p7_bg_Destroy(bgDNA);
  p7_bg_Destroy(bgAA);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7OPTACC_FS_BENCHMARK*/
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
 *   Compare p7_OptimalAccuracy_Frameshift and p7_OATrace_Frameshift
 *   against their scalar counterparts.
 *
 *   Strategy: run the full NEON pipeline (Forward->Backward->Decoding->OA->Trace)
 *   and the scalar pipeline (via deconverted pp matrix).  Compare OA scores,
 *   OA DP matrices, and tracebacks for equality within tolerance.
 *
 *   Note: Like utest_optacc in optacc.c, this test can fail for innocent
 *   reasons when numerical ties in the OA matrix lead to different traceback
 *   choices between NEON and scalar implementations.
 */
static void
utest_optacc_fs(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA,
                ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode,
                P7_BG *bgAA, P7_BG *bgDNA, int M, int L, int N)
{
  char            *msg   = "optacc_fs NEON unit test failed";
  P7_HMM          *hmm   = NULL;
  P7_FS_PROFILE   *gm_fs = NULL;
  P7_FS_OPROFILE  *om_fs = NULL;
  P7_GMX          *gx_pp  = NULL;   /* scalar pp (8-cell FS) */
  P7_GMX          *gx_oa  = NULL;   /* scalar OA (3-cell)    */
  P7_GMX          *gx_oa2 = NULL;   /* NEON OA deconverted   */
  P7_OMX          *ox_fwd = NULL;   /* NEON forward/pp (8-cell FS) */
  P7_OMX          *ox_bck = NULL;   /* NEON backward (3-cell)      */
  P7_TRACE        *tr_neon = NULL;  /* NEON OA traceback           */
  P7_TRACE        *tr_ref  = NULL;  /* scalar OA traceback         */
  ESL_DSQ         *dsq    = NULL;
  float            fsc, bsc, accscore_ref, accscore_neon;
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

  /* NEON: forward/backward use 8-cell FS layout; OA uses standard 3-cell */
  ox_fwd = p7_omx_Create_dpf(hmm->M, seq_L, seq_L, p7X_NSCELLS_FS);
  ox_bck = p7_omx_Create_dpf(hmm->M, seq_L, seq_L, p7X_NSCELLS);
  P7_OIVX *ov5 = p7_oivx_Create(hmm->M, p7P_5CODONS);

  /* Scalar: pp uses 8-cell FS layout; OA matrices use standard 3-cell */
  gx_pp  = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS_FS);
  gx_oa  = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS);
  gx_oa2 = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS);

  tr_neon = p7_trace_fs_CreateWithPP();
  tr_ref  = p7_trace_fs_CreateWithPP();

  while (N--)
    {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, seq_L, dsq);

      /* NEON pipeline: Forward -> Backward -> Decoding (in-place into ox_fwd) -> OA -> Trace */
      if (p7_Forward_Frameshift (dsq, seq_L, om_fs,         ox_fwd,        ov5, &fsc) == eslERANGE) continue;
      if (p7_Backward_Frameshift(dsq, seq_L, om_fs, ox_fwd, ox_bck, ov5, &bsc) == eslERANGE) continue;
      if (esl_FCompare_old(fsc, bsc, tol)                                        != eslOK) esl_fatal(msg);
      if (p7_Decoding_Frameshift(om_fs, ox_fwd, ox_bck)                      != eslOK) esl_fatal(msg);
      /* ox_fwd now holds the NEON posterior probability matrix */
      if (p7_OptimalAccuracy_Frameshift(om_fs, ox_fwd, ox_bck, &accscore_neon) != eslOK) esl_fatal(msg);
      if (p7_OATrace_Frameshift(om_fs, ox_fwd, ox_bck, tr_neon)               != eslOK) esl_fatal(msg);

      /* Reference: convert NEON pp to scalar GMX, run scalar OA and trace */
      if (p7_omx_FDeconvert(ox_fwd, gx_pp)                                      != eslOK) esl_fatal(msg);
      if (p7_GOptimalAccuracy_Frameshift(gm_fs, gx_pp, gx_oa, &accscore_ref)    != eslOK) esl_fatal(msg);
      if (p7_GOATrace_Frameshift(gm_fs, gx_pp, gx_oa, tr_ref)                   != eslOK) esl_fatal(msg);

      /* Convert NEON OA to scalar GMX for matrix-level comparison */
      if (p7_omx_FDeconvert(ox_bck, gx_oa2)                                     != eslOK) esl_fatal(msg);

      if (esl_FCompare_old(accscore_neon, accscore_ref, tol) != eslOK) esl_fatal(msg);
      if (p7_gmx_Compare(gx_oa, gx_oa2, tol)                != eslOK) esl_fatal(msg);
      if (p7_trace_Compare(tr_neon, tr_ref, pptol)           != eslOK) esl_fatal(msg);

      if (esl_opt_GetBoolean(go, "--traces"))
        {
          p7_trace_fs_Dump(stdout, tr_neon, gm_fs, dsq, abcDNA);
          p7_trace_fs_Dump(stdout, tr_ref,  gm_fs, dsq, abcDNA);
        }

      p7_trace_Reuse(tr_neon);
      p7_trace_Reuse(tr_ref);
    }

  free(dsq);
  p7_trace_fs_Destroy(tr_neon);
  p7_trace_fs_Destroy(tr_ref);
  p7_omx_Destroy(ox_fwd);
  p7_omx_Destroy(ox_bck);
  p7_oivx_Destroy(ov5);
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
   gcc -g -Wall -march=armv8-a -std=gnu99 -o optacc_fs_utest \
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
#include "impl_neon.h"

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
static char banner[] = "test driver for NEON frameshift optimal accuracy";

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
