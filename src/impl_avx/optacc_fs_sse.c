/* Frameshift optimal accuracy alignment; SSE implementation.
 * Ported from impl_sse/optacc_fs.c with _sse suffix for runtime dispatch.
 *
 * Contents:
 *   1. p7_OptimalAccuracy_Frameshift_sse() — DP fill
 *   2. p7_OATrace_Frameshift_sse()         — traceback
 */
#include "p7_config.h"

#ifdef eslENABLE_SSE

#include <float.h>
#include <string.h>

#include <xmmintrin.h>  /* SSE  */
#include <emmintrin.h>  /* SSE2 */

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"


/* codon_argmax_fs()
 * For the fill: the codon length the trace takes for an M in stripe <q>
 * of posterior row <ppc>, the argmax of pp_C1..C5 (the first one on
 * ties), as a float in each lane.
 */
static inline __m128
codon_argmax_fs(const __m128 *ppc, int q)
{
  __m128 best = MMO_FS(ppc, q, p7X_FS_C1);
  __m128 cv   = _mm_set1_ps(1.0f);
  __m128 gt;
  int    c;

  for (c = p7X_FS_C2; c <= p7X_FS_C5; c++)
    {
      gt   = _mm_cmpgt_ps(MMO_FS(ppc, q, c), best);
      best = _mm_max_ps(best, MMO_FS(ppc, q, c));
      cv   = esl_sse_select_ps(cv, _mm_set1_ps((float) c), gt);
    }
  return cv;
}

/*****************************************************************
 * 1. Optimal accuracy alignment, DP fill
 *****************************************************************/

/* p7_OptimalAccuracyRow_Frameshift_sse()
 * OA row <i> into <dpc> and row i of the OA special states <xmx>, from
 * OA rows i-1..i-5 in <dp1>..<dp5> (row 0 for rows before it) and the
 * posterior row <ppc> with its special states <ppx> (E, N, J, B, C of
 * row i); see the OA cell layout above. Row 0 is the initialization
 * (<ppc>, <ppx>, <dp1>..<dp5> unused). <dpc> may be <ppc>. Same
 * calculation as each row of p7_OptimalAccuracy_Frameshift_sse(), which calls it.
 */
void
p7_OptimalAccuracyRow_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const __m128 *ppc, const float *ppx, float *xmx, __m128 *dpc,
           const __m128 *dp1, const __m128 *dp2, const __m128 *dp3, const __m128 *dp4, const __m128 *dp5, int i)
{
  int    M   = om_fs->M;
  int    Q   = p7O_NQF(M);
  int    q, j;
  int    qM  = (M-1) % Q;   /* stripe containing k=M */
  int    rM  = (M-1) / Q;   /* lane  containing k=M */
  float  xN, xE, xJ, xC;
  float  t1, t2;
  register __m128 mpv1, dpv1, ipv1;   /* right-shifted prev row i-1, then updated to current q */
  register __m128 mpv2, dpv2, ipv2;   /* right-shifted prev row i-2 */
  register __m128 mpv3, dpv3, ipv3;   /* right-shifted prev row i-3 */
  register __m128 mpv4, dpv4, ipv4;   /* right-shifted prev row i-4 */
  register __m128 mpv5, dpv5, ipv5;   /* right-shifted prev row i-5 */
  register __m128 sv1, sv2, sv3, sv4, sv5, sv;
  register __m128 xEv, xBv1, xBv2, xBv3, xBv4, xBv5;
  register __m128 dcv;
  __m128 cdv;
  __m128 bm, mm, im, dm, md, mi, ii;
  __m128 *tp;
  __m128  zerov = _mm_setzero_ps();
  __m128  infv  = _mm_set1_ps(-eslINFINITY);
  union { __m128 v; float p[4]; } imfix;   /* for I(i,M) = -inf fixup */

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
  xBv1 = _mm_set1_ps(xmx[(i-1)*p7X_NXCELLS + p7X_B]);
  xBv2 = (i >= 2) ? _mm_set1_ps(xmx[(i-2)*p7X_NXCELLS + p7X_B]) : infv;
  xBv3 = (i >= 3) ? _mm_set1_ps(xmx[(i-3)*p7X_NXCELLS + p7X_B]) : infv;
  xBv4 = (i >= 4) ? _mm_set1_ps(xmx[(i-4)*p7X_NXCELLS + p7X_B]) : infv;
  xBv5 = (i >= 5) ? _mm_set1_ps(xmx[(i-5)*p7X_NXCELLS + p7X_B]) : infv;

  /* Right-shift last stripe of each prev row for k=0 (k-1 access at q=0) */
  mpv1 = esl_sse_rightshift_ps(MMO_OA(dp1, Q-1), infv);
  dpv1 = esl_sse_rightshift_ps(DMO_OA(dp1, Q-1), infv);
  ipv1 = esl_sse_rightshift_ps(IMO_OA(dp1, Q-1), infv);

  mpv2 = esl_sse_rightshift_ps(MMO_OA(dp2, Q-1), infv);
  dpv2 = esl_sse_rightshift_ps(DMO_OA(dp2, Q-1), infv);
  ipv2 = esl_sse_rightshift_ps(IMO_OA(dp2, Q-1), infv);

  mpv3 = esl_sse_rightshift_ps(MMO_OA(dp3, Q-1), infv);
  dpv3 = esl_sse_rightshift_ps(DMO_OA(dp3, Q-1), infv);
  ipv3 = esl_sse_rightshift_ps(IMO_OA(dp3, Q-1), infv);

  mpv4 = esl_sse_rightshift_ps(MMO_OA(dp4, Q-1), infv);
  dpv4 = esl_sse_rightshift_ps(DMO_OA(dp4, Q-1), infv);
  ipv4 = esl_sse_rightshift_ps(IMO_OA(dp4, Q-1), infv);

  mpv5 = esl_sse_rightshift_ps(MMO_OA(dp5, Q-1), infv);
  dpv5 = esl_sse_rightshift_ps(DMO_OA(dp5, Q-1), infv);
  ipv5 = esl_sse_rightshift_ps(IMO_OA(dp5, Q-1), infv);

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

      /* Update prev-row ptrs from right-shifted (k-1) to current q (k access) */
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
      dcv = _mm_and_ps(md, sv);

      /* I(i,k) = max(MI*M(i-3,k), II*I(i-3,k)) + pp_I(i,k); lag-3, same column */
      sv  =                _mm_and_ps(mi, mpv3);
      sv  = _mm_max_ps(sv, _mm_and_ps(ii, ipv3));
      IMO_OA(dpc, q) = _mm_add_ps(sv, IMO_FS(ppc, q));
    }

  /* I(i,M) has no I state; explicitly set to -inf */
  imfix.v      = IMO_OA(dpc, qM);
  imfix.p[rM]  = -eslINFINITY;
  IMO_OA(dpc, qM) = imfix.v;

  /* DD propagation: first pass */
  dcv = esl_sse_rightshift_ps(dcv, infv);
  tp  = om_fs->tfv + 7*Q;
  for (q = 0; q < Q; q++)
    {
      DMO_OA(dpc, q) = _mm_max_ps(dcv, DMO_OA(dpc, q));
      dcv         = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), DMO_OA(dpc, q));  tp++;
    }

  /* DD propagation: 3 more lazy passes */
  for (j = 1; j < 4; j++)
    {
      dcv = esl_sse_rightshift_ps(dcv, infv);
      tp  = om_fs->tfv + 7*Q;
      for (q = 0; q < Q; q++)
        {
          DMO_OA(dpc, q) = _mm_max_ps(dcv, DMO_OA(dpc, q));
          dcv         = _mm_and_ps(_mm_cmpgt_ps(*tp, zerov), dcv);  tp++;
        }
    }

  /* D->E */
  for (q = 0; q < Q; q++) xEv = _mm_max_ps(xEv, DMO_OA(dpc, q));

  esl_sse_hmax_ps(xEv, &xE);
  xmx[i*p7X_NXCELLS+p7X_E] = xE;

  /* N state: lag-3 */
  if (i > 2)
    xN = (om_fs->xf[p7O_N][p7O_LOOP] == 0.0f) ? 0.0f :
         xmx[(i-3)*p7X_NXCELLS+p7X_N] + ppx[p7X_N];
  else
    xN = (om_fs->xf[p7O_N][p7O_LOOP] == 0.0f) ? 0.0f :
         ppx[p7X_N];
  xmx[i*p7X_NXCELLS+p7X_N] = xN;

  /* J state: lag-3 */
  if (i > 2) {
    t1 = (om_fs->xf[p7O_J][p7O_LOOP] == 0.0f) ? 0.0f :
         xmx[(i-3)*p7X_NXCELLS+p7X_J] + ppx[p7X_J];
    t2 = (om_fs->xf[p7O_E][p7O_LOOP] == 0.0f) ? 0.0f : xE;
    xJ = ESL_MAX(t1, t2);
  } else {
    xJ = (om_fs->xf[p7O_E][p7O_LOOP] == 0.0f) ? 0.0f : xE;
  }
  xmx[i*p7X_NXCELLS+p7X_J] = xJ;

  /* C state: lag-3 */
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


/* Function:  p7_OptimalAccuracy_Frameshift_sse()
 *
 * Purpose:   SSE implementation of the DP fill step of optimal accuracy
 *            alignment for the frameshift model.  Five codon lengths (1-5 nt)
 *            are supported.  The posterior decoding matrix <pp> and the
 *            OA matrix <ox> use the 8-cell FS layout, and <ox> may be
 *            <pp> itself; see the OA cell layout above. The caller copies
 *            <pp->xmx> for the trace before the fill.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_OptimalAccuracy_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,
                                  P7_OMX *ox, float *ret_e)
{
  int    L   = pp->L;
  int    i;

  ox->M = om_fs->M;
  ox->L = L;

  p7_OptimalAccuracyRow_Frameshift_sse(om_fs, NULL, NULL, ox->xmx, ox->dpf[0], NULL, NULL, NULL, NULL, NULL, 0);
  /* row 0 (all -inf M, I, D) stands in for rows before it */
  for (i = 1; i <= L; i++)
    p7_OptimalAccuracyRow_Frameshift_sse(om_fs, pp->dpf[i], pp->xmx + i*p7X_NXCELLS, ox->xmx, ox->dpf[i],
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


/*****************************************************************
 * 2. OA traceback
 *****************************************************************/

static inline float
get_postprob_fs(const P7_OMX *ox, const float *ppx, int scur, int sprv, int k, int i)
{
  int   Q = p7O_NQF(ox->M);
  int   q = (k-1) % Q;
  int   r = (k-1) / Q;
  union { __m128 v; float p[4]; } u;

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

static inline int
select_m_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  __m128 *tp    = om_fs->tfv + 7*q;
  union { __m128 v; float p[4]; } u, tv;
  float   path[4];
  int     state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };
  __m128  mpv, dpv, ipv;

  if (q > 0) {
    mpv = MMO_OA(ox->dpf[i], q-1);
    dpv = DMO_OA(ox->dpf[i], q-1);
    ipv = IMO_OA(ox->dpf[i], q-1);
  } else {
    __m128 neg_inf = _mm_set1_ps(-eslINFINITY);
    mpv = esl_sse_rightshift_ps(MMO_OA(ox->dpf[i], Q-1), neg_inf);
    dpv = esl_sse_rightshift_ps(DMO_OA(ox->dpf[i], Q-1), neg_inf);
    ipv = esl_sse_rightshift_ps(IMO_OA(ox->dpf[i], Q-1), neg_inf);
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
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;
  int     r     = (k-1) / Q;
  union { __m128 v; float p[4]; } mpv, dpv, tmdv, tddv;
  float   path[2];

  if (q > 0) {
    mpv.v  = MMO_OA(ox->dpf[i], q-1);
    dpv.v  = DMO_OA(ox->dpf[i], q-1);
    tmdv.v = om_fs->tfv[7*(q-1) + p7O_MD];
    tddv.v = om_fs->tfv[7*Q     + (q-1)];
  } else {
    __m128 neg_inf = _mm_set1_ps(-eslINFINITY);
    mpv.v  = esl_sse_rightshift_ps(MMO_OA(ox->dpf[i], Q-1), neg_inf);
    dpv.v  = esl_sse_rightshift_ps(DMO_OA(ox->dpf[i], Q-1), neg_inf);
    tmdv.v = esl_sse_rightshiftz_float(om_fs->tfv[7*(Q-1) + p7O_MD]);
    tddv.v = esl_sse_rightshiftz_float(om_fs->tfv[8*Q - 1]);
  }

  path[0] = ((tmdv.p[r] == 0.0f) ? -eslINFINITY : mpv.p[r]);
  path[1] = ((tddv.p[r] == 0.0f) ? -eslINFINITY : dpv.p[r]);
  return ((path[0] >= path[1]) ? p7T_M : p7T_D);
}

static inline int
select_i_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int k)
{
  int     Q      = p7O_NQF(ox->M);
  int     q      = (k-1) % Q;
  int     r      = (k-1) / Q;
  int     prev_i = (i >= 3) ? i-3 : 0;
  __m128 *tp     = om_fs->tfv + 7*q + p7O_MI;
  union { __m128 v; float p[4]; } tv, mpv, ipv;
  float   path[2];

  mpv.v = MMO_OA(ox->dpf[prev_i], q);
  tv.v  = *tp;  path[0] = ((tv.p[r] == 0.0f) ? -eslINFINITY : mpv.p[r]);  tp++;
  ipv.v = IMO_OA(ox->dpf[prev_i], q);
  tv.v  = *tp;  path[1] = ((tv.p[r] == 0.0f) ? -eslINFINITY : ipv.p[r]);
  return ((path[0] >= path[1]) ? p7T_M : p7T_I);
}

static inline int
select_n_fs(int i)
{
  return ((i == 0) ? p7T_S : p7T_N);
}

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

static inline int
select_e_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i, int *ret_k)
{
  int    Q    = p7O_NQF(ox->M);
  union { __m128 v; float p[4]; } u;
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

static inline int
select_b_fs(const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, int i)
{
  float path[2];
  path[0] = (om_fs->xf[p7O_N][p7O_MOVE] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_N];
  path[1] = (om_fs->xf[p7O_J][p7O_MOVE] == 0.0f) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_J];
  return ((path[0] > path[1]) ? p7T_N : p7T_J);
}

static inline int
select_codon_fs(const P7_OMX *ox, int i, int k)
{
  int   Q = p7O_NQF(ox->M);
  int   q = (k-1) % Q;
  int   r = (k-1) / Q;
  union { __m128 v; float p[4]; } u;

  u.v = CODONO_OA(ox->dpf[i], q);
  return (int) u.p[r];
}


/* Function:  p7_OATrace_Frameshift_sse()
 *
 * Purpose:   SSE implementation of optimal accuracy traceback for the
 *            frameshift model.  Traces back through OA matrix <ox>,
 *            reading the posteriors of M and I from <ox> and those of the
 *            special states from <ppx>, the caller's copy of the
 *            decoding matrix's <xmx>.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if trace is not empty or a bogus state is reached.
 *            <eslEMEM> on allocation error.
 */
int
p7_OATrace_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const float *ppx,
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
      case p7T_C: scur = select_c_fs(om_fs, ppx, ox, i);                   break;
      case p7T_J: scur = select_j_fs(om_fs, ppx, ox, i);                   break;
      case p7T_E: scur = select_e_fs(om_fs,     ox, i, &k);               break;
      case p7T_B: scur = select_b_fs(om_fs,     ox, i);                   break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in OA FS traceback");
      }
      if (scur == -1) ESL_EXCEPTION(eslEINVAL, "OA FS SSE traceback choice failed");

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


#endif /* eslENABLE_SSE */
