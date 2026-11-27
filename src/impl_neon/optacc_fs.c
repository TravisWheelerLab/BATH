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
#include <string.h>

#include <arm_neon.h>

#include "easel.h"
#include "esl_neon.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_neon.h"


/* codon_argmax_fs()
 * For the fill: the codon length the trace takes for an M in stripe <q>
 * of posterior row <ppc>, the argmax of pp_C1..C5 (the first one on
 * ties), as a float in each lane.
 */
static inline float32x4_t
codon_argmax_fs(const float32x4_t *ppc, int q)
{
  float32x4_t best = MMO_FS(ppc, q, p7X_FS_C1);
  float32x4_t cv   = vdupq_n_f32(1.0f);
  uint32x4_t  gt;
  int         c;

  for (c = p7X_FS_C2; c <= p7X_FS_C5; c++)
    {
      gt   = vcgtq_f32(MMO_FS(ppc, q, c), best);
      best = vmaxq_f32(best, MMO_FS(ppc, q, c));
      cv   = vbslq_f32(gt, vdupq_n_f32((float) c), cv);
    }
  return cv;
}

/*****************************************************************
 * 1. Optimal accuracy alignment, DP fill
 *****************************************************************/

/* p7_OptimalAccuracyRow_Frameshift()
 * OA row <i> into <dpc> and row i of the OA special states <xmx>, from
 * OA rows i-1..i-5 in <dp1>..<dp5> (row 0 for rows before it) and the
 * posterior row <ppc> with its special states <ppx> (E, N, J, B, C of
 * row i); see the OA cell layout above. Row 0 is the initialization
 * (<ppc>, <ppx>, <dp1>..<dp5> unused). <dpc> may be <ppc>. Same
 * calculation as each row of p7_OptimalAccuracy_Frameshift(), which calls it.
 */
void
p7_OptimalAccuracyRow_Frameshift(const P7_FS_OPROFILE *om_fs, const float32x4_t *ppc, const float *ppx, float *xmx, float32x4_t *dpc,
           const float32x4_t *dp1, const float32x4_t *dp2, const float32x4_t *dp3, const float32x4_t *dp4, const float32x4_t *dp5, int i)
{
  int    M   = om_fs->M;
  int    Q   = p7O_NQF(M);
  int    q, j;
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
  float32x4_t cdv;
  float32x4_t bm, mm, im, dm, md, mi, ii;
  float32x4_t *tp;
  float32x4_t  zerov = vdupq_n_f32(0.0f);
  float32x4_t  infv  = vdupq_n_f32(-eslINFINITY);
  union { float32x4_t v; float p[4]; } imfix;   /* for I(i,M) = -inf fixup */

  if (i == 0)
    {
      for (q = 0; q < Q; q++) { MMO_OA(dpc,q) = DMO_OA(dpc,q) = IMO_OA(dpc,q) = infv; }
      xmx[0*p7X_NXCELLS+p7X_E] = -eslINFINITY;
      xmx[0*p7X_NXCELLS+p7X_N] =  0.0f;
      xmx[0*p7X_NXCELLS+p7X_J] = -eslINFINITY;
      xmx[0*p7X_NXCELLS+p7X_B] =  0.0f;
      xmx[0*p7X_NXCELLS+p7X_C] = -eslINFINITY;
      return;
    }

  /* B values from previous rows; -inf sentinel when lag exceeds i */
  xBv1 = vdupq_n_f32(xmx[(i-1)*p7X_NXCELLS + p7X_B]);
  xBv2 = (i >= 2) ? vdupq_n_f32(xmx[(i-2)*p7X_NXCELLS + p7X_B]) : infv;
  xBv3 = (i >= 3) ? vdupq_n_f32(xmx[(i-3)*p7X_NXCELLS + p7X_B]) : infv;
  xBv4 = (i >= 4) ? vdupq_n_f32(xmx[(i-4)*p7X_NXCELLS + p7X_B]) : infv;
  xBv5 = (i >= 5) ? vdupq_n_f32(xmx[(i-5)*p7X_NXCELLS + p7X_B]) : infv;

  /* Right-shift last stripe of each prev row for k=0 (k-1 access at q=0) */
  mpv1 = vextq_f32(infv, MMO_OA(dp1, Q-1), 3);
  dpv1 = vextq_f32(infv, DMO_OA(dp1, Q-1), 3);
  ipv1 = vextq_f32(infv, IMO_OA(dp1, Q-1), 3);

  mpv2 = vextq_f32(infv, MMO_OA(dp2, Q-1), 3);
  dpv2 = vextq_f32(infv, DMO_OA(dp2, Q-1), 3);
  ipv2 = vextq_f32(infv, IMO_OA(dp2, Q-1), 3);

  mpv3 = vextq_f32(infv, MMO_OA(dp3, Q-1), 3);
  dpv3 = vextq_f32(infv, DMO_OA(dp3, Q-1), 3);
  ipv3 = vextq_f32(infv, IMO_OA(dp3, Q-1), 3);

  mpv4 = vextq_f32(infv, MMO_OA(dp4, Q-1), 3);
  dpv4 = vextq_f32(infv, DMO_OA(dp4, Q-1), 3);
  ipv4 = vextq_f32(infv, IMO_OA(dp4, Q-1), 3);

  mpv5 = vextq_f32(infv, MMO_OA(dp5, Q-1), 3);
  dpv5 = vextq_f32(infv, DMO_OA(dp5, Q-1), 3);
  ipv5 = vextq_f32(infv, IMO_OA(dp5, Q-1), 3);

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
      mpv1 = MMO_OA(dp1, q);  dpv1 = DMO_OA(dp1, q);  ipv1 = IMO_OA(dp1, q);
      mpv2 = MMO_OA(dp2, q);  dpv2 = DMO_OA(dp2, q);  ipv2 = IMO_OA(dp2, q);
      mpv3 = MMO_OA(dp3, q);  dpv3 = DMO_OA(dp3, q);  ipv3 = IMO_OA(dp3, q);
      mpv4 = MMO_OA(dp4, q);  dpv4 = DMO_OA(dp4, q);  ipv4 = IMO_OA(dp4, q);
      mpv5 = MMO_OA(dp5, q);  dpv5 = DMO_OA(dp5, q);  ipv5 = IMO_OA(dp5, q);

      /* C3: the codon length for an M here, from pp_C1..C5 (read above).
       * C0, I, C4 and C5 keep the posteriors, copied when <ox> isn't <pp>. */
      cdv = codon_argmax_fs(ppc, q);
      if (dpc != ppc)
        {
          MMO_FS(dpc, q, p7X_FS_C0) = MMO_FS(ppc, q, p7X_FS_C0);
          MMO_FS(dpc, q, p7X_FS_C4) = MMO_FS(ppc, q, p7X_FS_C4);
          MMO_FS(dpc, q, p7X_FS_C5) = MMO_FS(ppc, q, p7X_FS_C5);
          IMO_FS(dpc, q)            = IMO_FS(ppc, q);
        }
      CODONO_OA(dpc, q) = cdv;

      MMO_OA(dpc, q) = sv;
      DMO_OA(dpc, q) = dcv;
      dcv = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(md), vreinterpretq_u32_f32(sv)));

      /* I(i,k) = max(MI*M(i-3,k), II*I(i-3,k)) + pp_I(i,k); lag-3, same column */
      sv  =                vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(mi), vreinterpretq_u32_f32(mpv3)));
      sv  = vmaxq_f32(sv, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(ii), vreinterpretq_u32_f32(ipv3))));
      IMO_OA(dpc, q) = vaddq_f32(sv, IMO_FS(ppc, q));
    }

  /* I(i,M) has no I state; MI/II transitions stored as 0 cause AND(0)=0.
   * The scalar explicitly sets IMX(i,M)=-inf; replicate that here.
   */
  imfix.v      = IMO_OA(dpc, qM);
  imfix.p[rM]  = -eslINFINITY;
  IMO_OA(dpc, qM) = imfix.v;

  /* DD propagation: first pass (includes M->D carry from main loop) */
  dcv = vextq_f32(infv, dcv, 3);
  tp  = om_fs->tfv + 7*Q;   /* DD transitions follow the 7*Q regular transitions */
  for (q = 0; q < Q; q++)
    {
      DMO_OA(dpc, q) = vmaxq_f32(dcv, DMO_OA(dpc, q));
      dcv         = vreinterpretq_f32_u32(vandq_u32(
                      vreinterpretq_u32_f32(vreinterpretq_f32_u32(vcgtq_f32(*tp, zerov))),
                      vreinterpretq_u32_f32(DMO_OA(dpc, q))));  tp++;
    }

  /* DD propagation: 3 more lazy passes */
  for (j = 1; j < 4; j++)
    {
      dcv = vextq_f32(infv, dcv, 3);
      tp  = om_fs->tfv + 7*Q;
      for (q = 0; q < Q; q++)
        {
          DMO_OA(dpc, q) = vmaxq_f32(dcv, DMO_OA(dpc, q));
          dcv         = vreinterpretq_f32_u32(vandq_u32(
                          vreinterpretq_u32_f32(vreinterpretq_f32_u32(vcgtq_f32(*tp, zerov))),
                          vreinterpretq_u32_f32(dcv)));  tp++;
        }
    }

  /* D->E */
  for (q = 0; q < Q; q++) xEv = vmaxq_f32(xEv, DMO_OA(dpc, q));

  /* E state */
  xE = esl_neon_hmax_f32((esl_neon_128f_t){ .f32x4 = xEv });
  xmx[i*p7X_NXCELLS+p7X_E] = xE;

  /* N state: lag-3 for i>=3; for i<=2 only pp contributes (no prior N term) */
  if (i > 2)
    xN = (om_fs->xf[p7O_N][p7O_LOOP] == 0.0f) ? 0.0f :
         xmx[(i-3)*p7X_NXCELLS+p7X_N] + ppx[p7X_N];
  else
    xN = (om_fs->xf[p7O_N][p7O_LOOP] == 0.0f) ? 0.0f :
         ppx[p7X_N];
  xmx[i*p7X_NXCELLS+p7X_N] = xN;

  /* J state: lag-3 for i>=3 */
  if (i > 2) {
    t1 = (om_fs->xf[p7O_J][p7O_LOOP] == 0.0f) ? 0.0f :
         xmx[(i-3)*p7X_NXCELLS+p7X_J] + ppx[p7X_J];
    t2 = (om_fs->xf[p7O_E][p7O_LOOP] == 0.0f) ? 0.0f : xE;
    xJ = ESL_MAX(t1, t2);
  } else {
    xJ = (om_fs->xf[p7O_E][p7O_LOOP] == 0.0f) ? 0.0f : xE;
  }
  xmx[i*p7X_NXCELLS+p7X_J] = xJ;

  /* C state: lag-3 for i>=3 */
  if (i > 2) {
    t1 = (om_fs->xf[p7O_C][p7O_LOOP] == 0.0f) ? 0.0f :
         xmx[(i-3)*p7X_NXCELLS+p7X_C] + ppx[p7X_C];
    t2 = (om_fs->xf[p7O_E][p7O_MOVE] == 0.0f) ? 0.0f : xE;
    xC = ESL_MAX(t1, t2);
  } else {
    xC = (om_fs->xf[p7O_E][p7O_MOVE] == 0.0f) ? 0.0f : xE;
  }
  xmx[i*p7X_NXCELLS+p7X_C] = xC;

  /* B state */
  t1 = (om_fs->xf[p7O_N][p7O_MOVE] == 0.0f) ? 0.0f : xN;
  t2 = (om_fs->xf[p7O_J][p7O_MOVE] == 0.0f) ? 0.0f : xJ;
  xmx[i*p7X_NXCELLS+p7X_B] = ESL_MAX(t1, t2);
}


/* Function:  p7_OptimalAccuracy_Frameshift()
 * Synopsis:  DP fill of optimal accuracy alignment, frameshift NEON version.
 *
 * Purpose:   Calculates the fill step of the optimal accuracy decoding
 *            algorithm for a frameshift-aware profile, using NEON SIMD
 *            parallelism over model positions.
 *
 *            The posterior decoding matrix <pp> must have the 8-cell FS
 *            layout (p7X_NSCELLS_FS=8; from p7_Decoding_Frameshift).
 *            The OA result matrix <ox> uses the same layout and may be
 *            <pp> itself; see the OA cell layout above. The caller copies
 *            <pp->xmx> for the trace before the fill.
 *
 *            Five codon lengths (1-5 nt) are supported. M(i,k) takes
 *            contributions from rows i-1..i-5; I(i,k) from row i-3;
 *            special states N/J/C are computed with a lag of 3.
 *
 * Args:      om_fs  - optimized frameshift profile
 *            pp     - posterior decoding matrix (8-cell FS layout)
 *            ox     - RESULT: OA DP matrix (8-cell FS layout); may be <pp>
 *            ret_e  - RETURN: expected number of correctly decoded positions
 *
 * Returns:   eslOK on success.
 */
int
p7_OptimalAccuracy_Frameshift(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                              P7_OMX *ox, float *ret_e)
{
  int    L   = pp->L;
  int    i;

  ox->M = om_fs->M;
  ox->L = L;

  p7_OptimalAccuracyRow_Frameshift(om_fs, NULL, NULL, ox->xmx, ox->dpf[0], NULL, NULL, NULL, NULL, NULL, 0);
  /* row 0 (all -inf M, I, D) stands in for rows before it */
  for (i = 1; i <= L; i++)
    p7_OptimalAccuracyRow_Frameshift(om_fs, pp->dpf[i], pp->xmx + i*p7X_NXCELLS, ox->xmx, ox->dpf[i],
              ox->dpf[i-1],
              (i >= 2) ? ox->dpf[i-2] : ox->dpf[0],
              (i >= 3) ? ox->dpf[i-3] : ox->dpf[0],
              (i >= 4) ? ox->dpf[i-4] : ox->dpf[0],
              (i >= 5) ? ox->dpf[i-5] : ox->dpf[0], i);

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
 *   position (i,k), given previous state <sprv>: pp_C0 or pp_I from the
 *   OA matrix <ox>, N/C/J from the copy <ppx> of the special states.
 */
static inline float
get_postprob_fs(const P7_OMX *ox, const float *ppx, int scur, int sprv, int k, int i)
{
  int   Q = p7O_NQF(ox->M);
  int   q = (k-1) % Q;
  int   r = (k-1) / Q;
  union { float32x4_t v; float p[4]; } u;

  switch (scur) {
  case p7T_M:
    u.v = MMO_FS(ox->dpf[i], q, p7X_FS_C0);  return u.p[r];
  case p7T_I:
    u.v = IMO_FS(ox->dpf[i], q);              return u.p[r];
  case p7T_N: if (sprv == scur) return ppx[i*p7X_NXCELLS+p7X_N]; break;
  case p7T_C: if (sprv == scur) return ppx[i*p7X_NXCELLS+p7X_C]; break;
  case p7T_J: if (sprv == scur) return ppx[i*p7X_NXCELLS+p7X_J]; break;
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
    mpv = MMO_OA(ox->dpf[i], q-1);
    dpv = DMO_OA(ox->dpf[i], q-1);
    ipv = IMO_OA(ox->dpf[i], q-1);
  } else {
    float32x4_t neg_inf = vdupq_n_f32(-eslINFINITY);
    mpv = vextq_f32(neg_inf, MMO_OA(ox->dpf[i], Q-1), 3);
    dpv = vextq_f32(neg_inf, DMO_OA(ox->dpf[i], Q-1), 3);
    ipv = vextq_f32(neg_inf, IMO_OA(ox->dpf[i], Q-1), 3);
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
    mpv.v  = MMO_OA(ox->dpf[i], q-1);
    dpv.v  = DMO_OA(ox->dpf[i], q-1);
    tmdv.v = om_fs->tfv[7*(q-1) + p7O_MD];
    tddv.v = om_fs->tfv[7*Q     + (q-1)];   /* DD block starts at 7*Q */
  } else {
    /* Use -inf fill so k=0 boundary gives -inf (matching scalar DMX(i,0)=-inf). */
    float32x4_t neg_inf = vdupq_n_f32(-eslINFINITY);
    float32x4_t zerov   = vdupq_n_f32(0.0f);
    mpv.v  = vextq_f32(neg_inf, MMO_OA(ox->dpf[i], Q-1), 3);
    dpv.v  = vextq_f32(neg_inf, DMO_OA(ox->dpf[i], Q-1), 3);
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

  mpv.v = MMO_OA(ox->dpf[prev_i], q);
  tv.v  = *tp;  path[0] = ((tv.p[r] == 0.0f) ? -eslINFINITY : mpv.p[r]);  tp++;
  ipv.v = IMO_OA(ox->dpf[prev_i], q);
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
select_c_fs(const P7_FS_OPROFILE *om_fs, const float *ppx, const P7_OMX *ox, int i)
{
  int   L     = ox->L;
  float t1    = (om_fs->xf[p7O_C][p7O_LOOP] == 0.0f) ? 0.0f : 1.0f;
  float t2    = (om_fs->xf[p7O_E][p7O_MOVE] == 0.0f) ? 0.0f : 1.0f;
  float path[4];
  int   state[4] = { p7T_C, p7T_C, p7T_C, p7T_E };

  if (i < 4) return p7T_E;

  path[0] = (t1 == 0.0f) ? -eslINFINITY :
            ox->xmx[(i-3)*p7X_NXCELLS+p7X_C] + ppx[ i   *p7X_NXCELLS+p7X_C];
  path[1] = (i < L && t1 != 0.0f) ?
            ox->xmx[(i-2)*p7X_NXCELLS+p7X_C] + ppx[(i+1)*p7X_NXCELLS+p7X_C] : -eslINFINITY;
  path[2] = (i < L-1 && t1 != 0.0f) ?
            ox->xmx[(i-1)*p7X_NXCELLS+p7X_C] + ppx[(i+2)*p7X_NXCELLS+p7X_C] : -eslINFINITY;
  path[3] = (t2 == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_E];
  return state[esl_vec_FArgMax(path, 4)];
}

/* select_j_fs:
 *   J(i) comes from J(i)+pp(i) (loop self-transition) or E(i).
 *   Returns p7T_J or p7T_E.
 */
static inline int
select_j_fs(const P7_FS_OPROFILE *om_fs, const float *ppx, const P7_OMX *ox, int i)
{
  float path[2];
  int   state[2] = { p7T_J, p7T_E };

  if (i <= 5) return p7T_E;

  path[0] = (om_fs->xf[p7O_J][p7O_LOOP] == 0.0f) ? -eslINFINITY :
            ox->xmx[i*p7X_NXCELLS+p7X_J] + ppx[i*p7X_NXCELLS+p7X_J];
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
    u.v = MMO_OA(ox->dpf[i], q);
    if (u.p[r] > max) { max = u.p[r]; smax = p7T_M; kmax = k; }
    u.v = DMO_OA(ox->dpf[i], q);
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
 *   At M(i,k), the codon length c=1..5 with the highest posterior, which
 *   the fill stored in the cell's C3 slot.
 */
static inline int
select_codon_fs(const P7_OMX *ox, int i, int k)
{
  int   Q = p7O_NQF(ox->M);
  int   q = (k-1) % Q;
  int   r = (k-1) / Q;
  union { float32x4_t v; float p[4]; } u;

  u.v = CODONO_OA(ox->dpf[i], q);
  return (int) u.p[r];
}


/* Function:  p7_OATrace_Frameshift()
 * Synopsis:  Optimal accuracy decoding traceback, frameshift NEON version.
 *
 * Purpose:   Traceback through the OA DP matrix <ox> filled by
 *            <p7_OptimalAccuracy_Frameshift>, annotated with posterior
 *            probabilities: those of M and I from <ox>, those of the
 *            special states from <ppx>, the caller's copy of the
 *            decoding matrix's <xmx>.
 *
 *            At each M state, <select_codon_fs> takes the codon length
 *            c=1..5 that the fill stored (the argmax of pp_C1..C5), and
 *            i is decremented by c after the step.  I states use lag-3.
 *            N/J/C special states use lag-3 nucleotide stepping.
 *
 * Args:      om_fs - optimized frameshift profile
 *            ppx   - copy of the posterior special states (decoding matrix <xmx>)
 *            ox    - OA DP matrix (8-cell FS layout, from p7_OptimalAccuracy_Frameshift)
 *            tr    - RESULT: OA traceback (must be empty; allocated with p7_trace_fs_CreateWithPP)
 *
 * Returns:   eslOK on success.
 *
 * Throws:    eslEINVAL if <tr> is not empty.
 *            eslEMEM on allocation error.
 */
int
p7_OATrace_Frameshift(const P7_FS_OPROFILE *om_fs, const float *ppx,
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
      case p7T_C: scur = select_c_fs(om_fs, ppx, ox, i);                   break;
      case p7T_J: scur = select_j_fs(om_fs, ppx, ox, i);                   break;
      case p7T_E: scur = select_e_fs(om_fs,     ox, i, &k);               break;
      case p7T_B: scur = select_b_fs(om_fs,     ox, i);                   break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in OA FS traceback");
      }
      if (scur == -1) ESL_EXCEPTION(eslEINVAL, "OA FS NEON traceback choice failed");

      postprob = get_postprob_fs(ox, ppx, scur, sprv, k, i);

      if (scur == p7T_M) c = select_codon_fs(ox, i, k);
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
  P7_OMX         *ox_bck  = NULL;   /* backward matrix                                  */
  P7_OMX         *ox_oa   = NULL;   /* OA matrix (8-cell FS layout)                     */
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
  ox_fwd = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS_FS);
  ox_bck = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS);
  ox_oa  = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS_FS);
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
      p7_OptimalAccuracy_Frameshift(om_fs5, ox_fwd, ox_oa, &accscore);

      if (! esl_opt_GetBoolean(go, "--notrace"))
        {
          p7_OATrace_Frameshift(om_fs5, ox_fwd->xmx, ox_oa, tr);
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
  p7_omx_Destroy(ox_oa);
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
#include <string.h>

#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"

/* oa_fs_compare()
 * Compare two frameshift OA matrices (scalar, or SIMD deconverted): every
 * slot of rows 1..L, the OA slots (C1, C2, D) of row 0, which is all the
 * fill writes there, and the special states. With <tol> 0 the values must
 * be identical.
 */
static int
oa_fs_compare(const P7_GMX *gx1, const P7_GMX *gx2, float tol)
{
  int i, k, s, x;

  if (gx1->M != gx2->M || gx1->L != gx2->L) return eslFAIL;

  for (i = 0; i <= gx1->L; i++)
    {
      for (k = 1; k <= gx1->M; k++)
        for (s = 0; s < p7G_NSCELLS_FS; s++)
          if (i > 0 || s == p7G_D || s == p7G_M + p7G_C1 || s == p7G_M + p7G_C2)
            if (esl_FCompare_old(gx1->dp[i][k*p7G_NSCELLS_FS + s], gx2->dp[i][k*p7G_NSCELLS_FS + s], tol) != eslOK) return eslFAIL;
      for (x = 0; x < p7G_NXCELLS; x++)
        if (esl_FCompare_old(gx1->xmx[i*p7G_NXCELLS + x], gx2->xmx[i*p7G_NXCELLS + x], tol) != eslOK) return eslFAIL;
    }
  return eslOK;
}

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
  P7_GMX          *gx_oa  = NULL;   /* scalar OA (8-cell FS) */
  P7_GMX          *gx_oa2 = NULL;   /* NEON OA deconverted (8-cell) */
  P7_OMX          *ox_fwd = NULL;   /* NEON forward/pp (8-cell FS) */
  P7_OMX          *ox_bck = NULL;   /* NEON backward (3-cell)      */
  P7_OMX          *ox_oa  = NULL;   /* SIMD OA (8-cell FS)         */
  P7_TRACE        *tr_neon = NULL;  /* NEON OA traceback           */
  P7_TRACE        *tr_ref  = NULL;  /* scalar OA traceback         */
  P7_TRACE        *tr_ip   = NULL;  /* scalar OA traceback, in place */
  float           *ppx     = NULL;  /* copy of the pp special states */
  ESL_DSQ         *dsq    = NULL;
  float            fsc, bsc, accscore_ref, accscore_neon, accscore_ip, accscore_simd_ip;
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

  /* SIMD: Forward and pp use the 8-cell FS layout, Backward 3-cell, OA 8-cell */
  ox_fwd = p7_omx_Create_dpf(hmm->M, seq_L, seq_L, p7X_NSCELLS_FS);
  ox_bck = p7_omx_Create_dpf(hmm->M, seq_L, seq_L, p7X_NSCELLS);
  ox_oa  = p7_omx_Create_dpf(hmm->M, seq_L, seq_L, p7X_NSCELLS_FS);
  P7_OIVX *ov5 = p7_oivx_Create(hmm->M, p7P_5CODONS);

  /* Scalar: pp, OA and the deconverted SIMD OA use the 8-cell FS layout */
  gx_pp  = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS_FS);
  gx_oa  = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS_FS);
  gx_oa2 = p7_gmx_Create(hmm->M, seq_L, seq_L, p7G_NSCELLS_FS);

  tr_neon = p7_trace_fs_CreateWithPP();
  tr_ref  = p7_trace_fs_CreateWithPP();
  tr_ip   = p7_trace_fs_CreateWithPP();
  ppx     = malloc(sizeof(float) * (seq_L+1) * p7X_NXCELLS);   /* p7X_NXCELLS >= p7G_NXCELLS */

  while (N--)
    {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, seq_L, dsq);

      /* NEON pipeline: Forward -> Backward -> Decoding (in-place into ox_fwd) -> OA -> Trace */
      if (p7_Forward_Frameshift (dsq, seq_L, om_fs,         ox_fwd,        ov5, &fsc) == eslERANGE) continue;
      if (p7_Backward_Frameshift(dsq, seq_L, om_fs, ox_fwd, ox_bck, ov5, &bsc) == eslERANGE) continue;
      if (esl_FCompare_old(fsc, bsc, tol)                                        != eslOK) esl_fatal(msg);
      if (p7_Decoding_Frameshift(om_fs, ox_fwd, ox_bck)                      != eslOK) esl_fatal(msg);
      /* ox_fwd now holds the NEON posterior probability matrix */
      if (p7_OptimalAccuracy_Frameshift(om_fs, ox_fwd, ox_oa, &accscore_neon)  != eslOK) esl_fatal(msg);
      if (p7_OATrace_Frameshift(om_fs, ox_fwd->xmx, ox_oa, tr_neon)          != eslOK) esl_fatal(msg);

      /* Reference: convert NEON pp to scalar GMX, run scalar OA and trace */
      if (p7_omx_FDeconvert(ox_fwd, gx_pp)                                      != eslOK) esl_fatal(msg);
      if (p7_GOptimalAccuracy_Frameshift(gm_fs, gx_pp, gx_oa, &accscore_ref)    != eslOK) esl_fatal(msg);
      if (p7_GOATrace_Frameshift(gm_fs, gx_pp->xmx, gx_oa, tr_ref)              != eslOK) esl_fatal(msg);

      /* In place: the scalar OA written over the posteriors, with a copy of
       * their special states for the trace, matches the separate OA exactly */
      esl_vec_FCopy(gx_pp->xmx, (seq_L+1) * p7G_NXCELLS, ppx);
      if (p7_GOptimalAccuracy_Frameshift(gm_fs, gx_pp, gx_pp, &accscore_ip)  != eslOK) esl_fatal(msg);
      if (p7_GOATrace_Frameshift(gm_fs, ppx, gx_pp, tr_ip)                   != eslOK) esl_fatal(msg);
      if (accscore_ip != accscore_ref)                                                 esl_fatal(msg);
      if (oa_fs_compare(gx_pp, gx_oa, 0.0)                                    != eslOK) esl_fatal(msg);
      if (p7_trace_Compare(tr_ip, tr_ref, 0.0)                               != eslOK) esl_fatal(msg);

      /* Convert NEON OA to scalar GMX for matrix-level comparison */
      if (p7_omx_FDeconvert(ox_oa, gx_oa2)                                      != eslOK) esl_fatal(msg);

      if (esl_FCompare_old(accscore_neon, accscore_ref, tol) != eslOK) esl_fatal(msg);
      if (oa_fs_compare(gx_oa, gx_oa2, tol)                 != eslOK) esl_fatal(msg);
      if (p7_trace_Compare(tr_neon, tr_ref, pptol)           != eslOK) esl_fatal(msg);

      /* In place: the SIMD OA written over the posteriors, with a copy of
       * their special states for the trace, matches the separate OA exactly */
      p7_trace_Reuse(tr_ip);
      esl_vec_FCopy(ox_fwd->xmx, (seq_L+1) * p7X_NXCELLS, ppx);
      if (p7_OptimalAccuracy_Frameshift(om_fs, ox_fwd, ox_fwd, &accscore_simd_ip) != eslOK) esl_fatal(msg);
      if (p7_OATrace_Frameshift(om_fs, ppx, ox_fwd, tr_ip)                        != eslOK) esl_fatal(msg);
      if (accscore_simd_ip != accscore_neon)                                                esl_fatal(msg);
      if (p7_omx_FDeconvert(ox_fwd, gx_pp)                                        != eslOK) esl_fatal(msg);
      if (oa_fs_compare(gx_pp, gx_oa2, 0.0)                                       != eslOK) esl_fatal(msg);
      if (p7_trace_Compare(tr_ip, tr_neon, 0.0)                                   != eslOK) esl_fatal(msg);

      if (esl_opt_GetBoolean(go, "--traces"))
        {
          p7_trace_fs_Dump(stdout, tr_neon, gm_fs, dsq, abcDNA);
          p7_trace_fs_Dump(stdout, tr_ref,  gm_fs, dsq, abcDNA);
        }

      p7_trace_Reuse(tr_neon);
      p7_trace_Reuse(tr_ref);
      p7_trace_Reuse(tr_ip);
    }

  free(dsq);
  p7_trace_fs_Destroy(tr_neon);
  p7_trace_fs_Destroy(tr_ref);
  p7_trace_fs_Destroy(tr_ip);
  free(ppx);
  p7_omx_Destroy(ox_fwd);
  p7_omx_Destroy(ox_bck);
  p7_omx_Destroy(ox_oa);
  p7_oivx_Destroy(ov5);
  p7_gmx_Destroy(gx_pp);
  p7_gmx_Destroy(gx_oa);
  p7_gmx_Destroy(gx_oa2);
  p7_fs_oprofile_Destroy(om_fs);
  p7_profile_fs_Destroy(gm_fs);
  p7_hmm_Destroy(hmm);
}
/* pp_cell_fs()
 * Value of cell <s> (frameshift layout, e.g. p7X_FS_M + p7X_FS_C0) of
 * model position <k> on row <i> of <ox>.
 */
static float
pp_cell_fs(const P7_OMX *ox, int i, int k, int s)
{
  int Q = p7O_NQF(ox->M);
  int q = (k-1) % Q;
  int r = (k-1) / Q;
  const float *v = (const float *) (ox->dpf[i] + q * p7X_NSCELLS_FS + s);
  return v[r];
}

/* utest_sweep_fs()
 * The frameshift OA sweep over a checkpointed Backward against the
 * full-matrix path on the same Backward: p7_Forward_Frameshift(), the
 * checkpointed Backward filled in one block and copied into a Backward
 * matrix, p7_Decoding_Frameshift(), p7_OptimalAccuracy_Frameshift(),
 * p7_OATrace_Frameshift() and p7_Null2_fs_ByExpectation(). The sweep,
 * with Backward blocks of 5..12 rows, must give the same Forward and OA
 * scores, null2, states, k, i and codon lengths on the trace, and PP
 * characters equal to those the posterior matrix gives. Sequences
 * alternate random DNA and the back-translated consensus of a model
 * with sharpened match emissions and M->M transitions between random
 * flanks, so that rows are rescaled.
 */
static void
utest_sweep_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA,
               int M, int L, int N)
{
  char           *msg  = "frameshift OA sweep unit test failed";
  P7_CODONTABLE  *ct   = p7_codontable_Create(gcode);
  P7_HMM         *hmm  = NULL;
  P7_FS_PROFILE  *gm   = p7_profile_fs_Create(M, abcAA, p7P_5CODONS);
  P7_FS_OPROFILE *om   = NULL;
  P7_OIVX        *ov   = NULL;
  ESL_SQ         *sq   = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq  = malloc(sizeof(ESL_DSQ) * (L+2));
  float          *n2a  = malloc(sizeof(float) * abcAA->Kp);
  float          *n2b  = malloc(sizeof(float) * abcAA->Kp);
  float          *ppx  = malloc(sizeof(float) * (L+1) * p7X_NXCELLS);
  P7_OMX         *fwd  = NULL, *bck = NULL, *oa = NULL;
  P7_CHKOMX      *cx1  = NULL, *cx2 = NULL;
  P7_OATMX       *ox   = p7_oatmx_Create(M, L, p7OAT_NPLANES_FS, 1);
  P7_TRACE       *tr1  = p7_trace_fs_CreateWithPP();
  P7_TRACE       *tr2  = p7_trace_fs_CreateWithPP();
  int             Q    = p7O_NQF(M);
  size_t          rowb = sizeof(float32x4_t) * p7X_NSCELLS * (size_t) Q;
  float           fsc, fsc2, oasc1, oasc2, p;
  int             n, Ln, blockR, i, j, k, x, z, ntested = 0;

  if (p7_hmm_Sample(r, M, abcAA, &hmm) != eslOK) esl_fatal(msg);
  for (k = 1; k <= M; k++)
    {
      x = esl_vec_FArgMax(hmm->mat[k], abcAA->K);
      esl_vec_FSet(hmm->mat[k], abcAA->K, 0.15 / (abcAA->K - 1));
      hmm->mat[k][x] = 0.85;
      if (k < M) { hmm->t[k][p7H_MM] = 0.9; hmm->t[k][p7H_MI] = 0.05; hmm->t[k][p7H_MD] = 0.05; }
    }
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm, L/3, p7_LOCAL) != eslOK) esl_fatal(msg);
  p7_emit_SimpleConsensus(hmm, sq);

  om  = p7_fs_oprofile_Create(M, abcAA, p7P_5CODONS);  p7_fs_oprofile_Convert(gm, om);
  ov  = p7_oivx_Create(M, p7P_5CODONS);
  fwd = p7_omx_Create_dpf(M, L, L, p7X_NSCELLS_FS); bck = p7_omx_Create_dpf(M, L, L, p7X_NSCELLS); oa = p7_omx_Create_dpf(M, L, L, p7X_NSCELLS_FS);
  cx1 = p7_chkomx_Create(M, L, L+1, 5, p7X_NSCELLS, p7_CHK_BACKWARD);
  cx2 = p7_chkomx_Create(M, L, 5,   5, p7X_NSCELLS, p7_CHK_BACKWARD);

  for (n = 0; n < N; n++)
    {
      Ln = 15 + esl_rnd_Roll(r, L - 14);
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, dsq);
      if (n % 2 == 1 && 3*M <= L)
        {
          Ln = ESL_MAX(Ln, 3*M);
          esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, dsq);
          for (i = 1, j = 1 + esl_rnd_Roll(r, Ln - 3*M + 1); i <= M; i++, j += 3)
            p7_codontable_GetCodon(ct, r, sq->dsq[i], dsq+j);
        }
      dsq[0] = dsq[Ln+1] = eslDSQ_SENTINEL;
      blockR = 5 + esl_rnd_Roll(r, 8);
      p7_trace_Reuse(tr1);
      p7_trace_Reuse(tr2);

      p7_fs_oprofile_ReconfigLength(om, Ln/3);
      if (p7_Forward_Frameshift(dsq, Ln, om, fwd, ov, &fsc)       != eslOK) continue;
      if (p7_chkomx_GrowTo(cx1, M, Ln, Ln+1)                       != eslOK) esl_fatal(msg);
      if (p7_BackwardChk_Frameshift(dsq, Ln, om, cx1, ov, NULL)    != eslOK) continue;
      for (i = 0; i <= Ln; i++) memcpy(bck->dpf[i], cx1->dpf[i], rowb);
      memcpy(bck->xmx, cx1->xmx, sizeof(float) * (Ln+1) * p7X_NXCELLS);
      bck->M = M; bck->L = Ln;

      if (p7_Decoding_Frameshift(om, fwd, bck)                     != eslOK) esl_fatal(msg);
      memcpy(ppx, fwd->xmx, sizeof(float) * (Ln+1) * p7X_NXCELLS);
      p7_OptimalAccuracy_Frameshift(om, fwd, oa, &oasc1);
      if (p7_OATrace_Frameshift(om, ppx, oa, tr1)                  != eslOK) esl_fatal(msg);
      if (p7_chkomx_GrowTo(cx2, M, Ln, blockR)                     != eslOK) esl_fatal(msg);
      if (p7_BackwardChk_Frameshift(dsq, Ln, om, cx2, ov, NULL)    != eslOK) esl_fatal(msg);
      if (p7_OptimalAccuracySweep_Frameshift(dsq, Ln, om, cx2, ov, ox, n2b, &fsc2, &oasc2) != eslOK) esl_fatal(msg);
      if (p7_OATraceSweep_Frameshift(om, ox, tr2)                  != eslOK) esl_fatal(msg);

      if (fsc2  != fsc)    esl_fatal("%s: Forward scores differ (%f %f)", msg, fsc2, fsc);
      if (oasc2 != oasc1)  esl_fatal("%s: OA scores differ (%f %f)", msg, oasc2, oasc1);
      if (tr1->N != tr2->N) esl_fatal("%s: trace lengths differ (%d %d)", msg, tr1->N, tr2->N);
      for (z = 0; z < tr1->N; z++)
        {
          if (tr1->st[z] != tr2->st[z] || tr1->k[z] != tr2->k[z] || tr1->i[z] != tr2->i[z] || tr1->c[z] != tr2->c[z])
            esl_fatal("%s: traces differ at %d", msg, z);
          switch (tr1->st[z]) {
          case p7T_M: p = pp_cell_fs(fwd, tr1->i[z], tr1->k[z], p7X_FS_M + p7X_FS_C0); break;
          case p7T_I: p = pp_cell_fs(fwd, tr1->i[z], tr1->k[z], p7X_FS_I); break;
          case p7T_N: case p7T_C: case p7T_J:
            /* an N, C or J emits residue i when the state before it is the same one */
            p = (z > 0 && tr1->st[z-1] == tr1->st[z]) ? ppx[tr1->i[z]*p7X_NXCELLS + (tr1->st[z] == p7T_N ? p7X_N : tr1->st[z] == p7T_C ? p7X_C : p7X_J)] : 0.0;
            break;
          default: p = 0.0; break;
          }
          if (p7_alidisplay_EncodePostProb(tr2->pp[z]) != p7_alidisplay_EncodePostProb(p))
            esl_fatal("%s: PP character differs at %d (state %d, i %d, k %d)", msg, z, tr1->st[z], tr1->i[z], tr1->k[z]);
        }

      /* null2 from the full posterior matrix, which it uses as workspace, so last */
      p7_Null2_fs_ByExpectation(om, fwd, n2a);
      for (x = 0; x < abcAA->Kp; x++)
        if (n2b[x] != n2a[x]) esl_fatal("%s: null2 differs for residue %d (%g %g)", msg, x, n2b[x], n2a[x]);
      ntested++;
    }
  if (ntested == 0) esl_fatal("%s: no sequence was tested", msg);

  free(dsq); free(n2a); free(n2b); free(ppx);
  esl_sq_Destroy(sq);
  p7_trace_fs_Destroy(tr1);
  p7_trace_fs_Destroy(tr2);
  p7_oatmx_Destroy(ox);
  p7_omx_Destroy(fwd); p7_omx_Destroy(bck); p7_omx_Destroy(oa);
  p7_chkomx_Destroy(cx1); p7_chkomx_Destroy(cx2);
  p7_oivx_Destroy(ov);
  p7_fs_oprofile_Destroy(om);
  p7_profile_fs_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_codontable_Destroy(ct);
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
  utest_sweep_fs(r, abcAA, abcDNA, gcode, bgAA, bgDNA, M, 3*M+150, N);
  utest_sweep_fs(r, abcAA, abcDNA, gcode, bgAA, bgDNA, 1, 60,      10);

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
