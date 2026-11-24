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



/* The OA matrix uses the frameshift (p7X_NSCELLS_FS) cell layout, so the
 * fill can write it over the posterior matrix. Each filled cell holds
 * the OA scores of M, I and D in the C1, C2 and D slots, the codon length
 * the trace takes for an M there in C3, and the posteriors pp_C0 and
 * pp_I, which the trace reports, in C0 and I. C4 and C5 keep pp_C4 and
 * pp_C5, which nothing reads after the fill.
 */
#define MMO_OA(dp,q)   MMO_FS(dp,q,p7X_FS_C1)
#define IMO_OA(dp,q)   MMO_FS(dp,q,p7X_FS_C2)
#define DMO_OA(dp,q)   DMO_FS(dp,q)
#define CODON_OA(dp,q) MMO_FS(dp,q,p7X_FS_C3)

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

/* oa_fs_row_sse()
 * OA row <i> into <dpc> and row i of the OA special states <xmx>, from
 * OA rows i-1..i-5 in <dp1>..<dp5> (row 0 for rows before it) and the
 * posterior row <ppc> with its special states <ppx> (E, N, J, B, C of
 * row i); see the OA cell layout above. Row 0 is the initialization
 * (<ppc>, <ppx>, <dp1>..<dp5> unused). <dpc> may be <ppc>. Same
 * calculation as each row of p7_OptimalAccuracy_Frameshift_sse(), which calls it.
 */
static void
oa_fs_row_sse(const P7_FS_OPROFILE *om_fs, const __m128 *ppc, const float *ppx, float *xmx, __m128 *dpc,
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
      CODON_OA(dpc, q) = cdv;

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

  oa_fs_row_sse(om_fs, NULL, NULL, ox->xmx, ox->dpf[0], NULL, NULL, NULL, NULL, NULL, 0);
  /* row 0 (all -inf M, I, D) stands in for rows before it */
  for (i = 1; i <= L; i++)
    oa_fs_row_sse(om_fs, pp->dpf[i], pp->xmx + i*p7X_NXCELLS, ox->xmx, ox->dpf[i],
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

  u.v = CODON_OA(ox->dpf[i], q);
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

/*****************************************************************
 * OA sweep over a checkpointed Backward, and its trace
 *****************************************************************/

/* oat_ppcodes4_fs()
 * Posterior codes (p7OAT_PPCODE(), 0..10) of the 4 lanes of <v> as 4
 * lane masks, bit b of each lane's code in <m[b]>.
 */
static inline void
oat_ppcodes4_fs(__m128 v, int m[4])
{
  union { __m128 v; float p[4]; } u;
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

/* oat_store_row_fs_sse()
 * Store into <ox> what the frameshift OA trace reads from OA row <i>
 * (<dpc>, in the OA cell layout above, with the OA special states in
 * <xmx>): for each k, the predecessors of an M at (i+c,k), an I at
 * (i+3,k) and a D at (i,k), chosen with the same comparisons, shifts
 * and ties as p7_OATrace_Frameshift_sse() makes on the OA matrix; the
 * codon length of an M at (i,k), from the fill; the posterior codes of
 * M (pp_C0) and I at (i,k) from the posterior row <ppp> (all 0 for row
 * 0, <ppp> NULL); the row's OA special states and the posteriors of N,
 * J, C (<ppx>, NULL for row 0); and the cell an E at row i comes from,
 * in p7_OATrace_Frameshift_sse()'s scan order.
 */
/* oat_blend_fs_sse()
 * Lanes of <b> where <mask> is set, of <a> elsewhere (SSE2).
 */
static inline __m128
oat_blend_fs_sse(__m128 a, __m128 b, __m128 mask)
{
  return _mm_or_ps(_mm_and_ps(mask, b), _mm_andnot_ps(mask, a));
}

static void
oat_store_row_fs_sse(const P7_FS_OPROFILE *om_fs, P7_OATMX *ox, const __m128 *dpc, const float *xmx, const __m128 *ppp, const float *ppx, int i)
{
  int      M     = om_fs->M;
  int      Q     = p7O_NQF(M);
  __m128   zerov = _mm_setzero_ps();
  __m128   infv  = _mm_set1_ps(-eslINFINITY);
  __m128   xBv   = _mm_set1_ps(xmx[i*p7X_NXCELLS+p7X_B]);
  __m128   mpv, ipv, dpv, tmdv, tddv;
  __m128   pM, pI, pD, pB, best, gt, b0, b1;
  __m128   cv;
  __m128  *tp;
  union { __m128 v; float p[4]; } u;
  float    max;
  int      kmax, smax;
  int      mm[4], mi[4], mc;
  int      q, r, c, k;

  for (q = 0; q < Q; q++)
    {
      tp = om_fs->tfv + 7*q;

      /* row i at k-1 */
      if (q > 0) {
        mpv  = MMO_OA(dpc, q-1);
        ipv  = IMO_OA(dpc, q-1);
        dpv  = DMO_OA(dpc, q-1);
        tmdv = om_fs->tfv[7*(q-1) + p7O_MD];
        tddv = om_fs->tfv[7*Q + (q-1)];
      } else {
        mpv  = esl_sse_rightshift_ps(MMO_OA(dpc, Q-1), infv);
        ipv  = esl_sse_rightshift_ps(IMO_OA(dpc, Q-1), infv);
        dpv  = esl_sse_rightshift_ps(DMO_OA(dpc, Q-1), infv);
        tmdv = esl_sse_rightshiftz_float(om_fs->tfv[7*(Q-1) + p7O_MD]);
        tddv = esl_sse_rightshiftz_float(om_fs->tfv[8*Q - 1]);
      }

      /* M at (i+c,k): first best of M, I, D, B */
      pB   = oat_blend_fs_sse(xBv, infv, _mm_cmpeq_ps(tp[p7O_BM], zerov));
      pM   = oat_blend_fs_sse(mpv, infv, _mm_cmpeq_ps(tp[p7O_MM], zerov));
      pI   = oat_blend_fs_sse(ipv, infv, _mm_cmpeq_ps(tp[p7O_IM], zerov));
      pD   = oat_blend_fs_sse(dpv, infv, _mm_cmpeq_ps(tp[p7O_DM], zerov));
      best = pM;
      gt   = _mm_cmpgt_ps(pI, best);  best = oat_blend_fs_sse(best, pI, gt);  b0 = gt;                    b1 = zerov;
      gt   = _mm_cmpgt_ps(pD, best);  best = oat_blend_fs_sse(best, pD, gt);  b0 = _mm_andnot_ps(gt, b0); b1 = _mm_or_ps(b1, gt);
      gt   = _mm_cmpgt_ps(pB, best);                                           b0 = _mm_or_ps(b0, gt);    b1 = _mm_or_ps(b1, gt);
      OATO_SET_STRIPE4(ox, i, p7OAT_MP, q, _mm_movemask_ps(b0));
      OATO_SET_STRIPE4(ox, i, p7OAT_MP+1, q, _mm_movemask_ps(b1));

      /* I at (i+3,k): M on ties */
      pM = oat_blend_fs_sse(MMO_OA(dpc, q), infv, _mm_cmpeq_ps(tp[p7O_MI], zerov));
      pI = oat_blend_fs_sse(IMO_OA(dpc, q), infv, _mm_cmpeq_ps(tp[p7O_II], zerov));
      OATO_SET_STRIPE4(ox, i, p7OAT_IP, q, _mm_movemask_ps(_mm_cmpgt_ps(pI, pM)));

      /* D at (i,k): M on ties */
      pM = oat_blend_fs_sse(mpv, infv, _mm_cmpeq_ps(tmdv, zerov));
      pD = oat_blend_fs_sse(dpv, infv, _mm_cmpeq_ps(tddv, zerov));
      OATO_SET_STRIPE4(ox, i, p7OAT_DP, q, _mm_movemask_ps(_mm_cmpgt_ps(pD, pM)));

      /* codon length of an M here (1..5), as the fill stored it */
      cv = CODON_OA(dpc, q);
      for (c = 0; c < 3; c++)
        {
          u.v = cv; mc = 0;
          for (r = 0; r < 4; r++) mc |= ((((int) u.p[r]) >> c) & 1) << r;
          OATO_SET_STRIPE4(ox, i, p7OAT_COD + c, q, mc);
        }

      /* posterior codes of M and I */
      if (ppp != NULL) { oat_ppcodes4_fs(MMO_FS(ppp, q, p7X_FS_C0), mm); oat_ppcodes4_fs(IMO_FS(ppp, q), mi); }
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


/* Function:  p7_OptimalAccuracySweep_Frameshift_sse()
 * Synopsis:  Frameshift OA alignment in one forward sweep over a checkpointed Backward, SSE.
 *
 * Purpose:   For the 5-codon frameshift profile <om_fs> and DNA sequence
 *            <dsq> of length <L>, with the checkpointed Backward matrix
 *            <cx> filled by <p7_BackwardChk_Frameshift_sse()>, compute in
 *            one sweep through rows 1..L the Forward, posterior and
 *            optimal accuracy (OA) rows, computing each block of
 *            Backward rows again from its checkpoint as the sweep reaches
 *            it, and store in <ox> what <p7_OATraceSweep_Frameshift_sse()>
 *            needs for the OA trace. <ov> is workspace for 5 codon
 *            lengths. Rows of Forward, posterior and OA values are kept
 *            only as long as later rows need them (up to 5 rows back).
 *
 *            Each row is calculated as <p7_Forward_Frameshift_sse()>,
 *            <p7_Decoding_Frameshift_sse()> and
 *            <p7_OptimalAccuracy_Frameshift_sse()> calculate it; the
 *            posteriors are those <p7_Decoding_Frameshift_sse()> gives for
 *            the checkpointed Backward, which has its own scale factors.
 *
 *            If <opt_null2> is non-NULL, the null2 odds ratios of the
 *            sequence are also calculated into it, as
 *            <p7_Null2_fs_ByExpectation_sse()> calculates them from the
 *            posteriors: the posterior rows are summed as the sweep
 *            computes them, in the same order.
 *
 *            <ox> must have been created with <p7OAT_NPLANES_FS> planes;
 *            its layout is set to 4 lanes and it is grown for <om_fs->M>
 *            by <L>. On return <cx> has its last block loaded.
 *
 * Args:      dsq       - DNA sequence in digitized form, 1..L
 *            L         - length of dsq
 *            om_fs     - 5-codon frameshift profile
 *            cx        - checkpointed Backward matrix from p7_BackwardChk_Frameshift_sse()
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
 *            <p7_Forward_Frameshift_sse()> and <p7_Decoding_Frameshift_sse()>).
 */
int
p7_OptimalAccuracySweep_Frameshift_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, P7_OATMX *ox,
                                       float *opt_null2, float *opt_fsc, float *ret_oasc)
{
  int            M        = om_fs->M;
  int            Q        = p7O_NQF(M);
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
  __m128        *pp, *cnt;
  __m128        *oa[6];
  int            i, q, b;
  int            status;

  if (om_fs->codon_lengths != 5)        ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");
  if (ox->nplanes != p7OAT_NPLANES_FS)  ESL_EXCEPTION(eslEINVAL, "P7_OATMX not set up for the frameshift OA sweep");
  if (cx->L != L || cx->M != M || cx->dir != p7_CHK_BACKWARD || cx->nscells != p7X_NSCELLS) ESL_EXCEPTION(eslEINVAL, "checkpointed Backward matrix doesn't match");
  ox->lanes = 4;
  if ((status = p7_oatmx_GrowTo(ox, M, L)) != eslOK) return status;

  if ((wk = p7_omx_Create_dpf_sse(M, 11, L, p7X_NSCELLS_FS)) == NULL) { status = eslEMEM; goto ERROR; }
  if ((bv = p7_oivx_Create_sse(M, p7P_5CODONS))              == NULL) { status = eslEMEM; goto ERROR; }
  ESL_ALLOC(oxmx,     sizeof(float) * (L+1) * p7X_NXCELLS);
  ESL_ALLOC(log_sfwd, sizeof(float) * (L+2));
  ESL_ALLOC(log_sbck, sizeof(float) * (L+2));
  for (q = 0; q < 6; q++) oa[q] = wk->dpf[4+q];
  pp  = wk->dpf[10];
  cnt = wk->dpf[11];
  for (q = 0; q < p7X_NXCELLS; q++) cntx[q] = 0.0f;

  /* Backward's cumulative log scale products and normalization, as p7_Decoding_Frameshift_sse() computes them */
  log_sbck[L+1] = 0.0f;
  for (i = L; i >= 0; i--)
    log_sbck[i] = log_sbck[i+1] + logf(cx->xmx[i*p7X_NXCELLS + p7X_SCALE]);
  log_inv_Z = -p7_FLogsum(
                 logf(cx->xmx[p7X_N]) + log_sbck[0],
                 p7_FLogsum(
                   logf(cx->xmx[1*p7X_NXCELLS + p7X_N]) + log_sbck[1],
                   logf(cx->xmx[2*p7X_NXCELLS + p7X_N]) + log_sbck[2]));

  p7_ForwardRow_Frameshift_sse(dsq, om_fs, wk->xmx, &totscale, &rg, ov->ivx, wk->dpf[0], NULL, NULL, 0);
  log_sfwd[0] = logf(wk->xmx[p7X_SCALE]);
  oa_fs_row_sse(om_fs, NULL, NULL, oxmx, oa[0], NULL, NULL, NULL, NULL, NULL, 0);
  oat_store_row_fs_sse(om_fs, ox, oa[0], oxmx, NULL, NULL, 0);

  for (i = 1; i <= L; i++)
    {
      b = (i / cx->blockR) * cx->blockR;
      if (b != cx->b && (status = p7_BackwardChkLoad_Frameshift_sse(dsq, om_fs, cx, bv, b)) != eslOK) goto ERROR;

      p7_ForwardRow_Frameshift_sse(dsq, om_fs, wk->xmx, &totscale, &rg, ov->ivx, wk->dpf[i%4], wk->dpf[(i-1)%4],
                                   (i >= 3) ? wk->dpf[(i-3)%4] : NULL, i);
      log_sfwd[i] = log_sfwd[i-1] + logf(wk->xmx[i*p7X_NXCELLS + p7X_SCALE]);

      if ((status = p7_DecodingRow_Frameshift_sse(om_fs, cx->xmx, (i > 2) ? wk->xmx + (i-3)*p7X_NXCELLS : NULL, log_sfwd, log_sbck, log_inv_Z,
                                                  wk->dpf[i%4], cx->dpf[i], pp, ppx, i)) != eslOK) goto ERROR;

      /* row 0 (all -inf M, I, D) stands in for rows before it */
      oa_fs_row_sse(om_fs, pp, ppx, oxmx, oa[i%6], oa[(i-1)%6],
                    (i >= 2) ? oa[(i-2)%6] : oa[0],
                    (i >= 3) ? oa[(i-3)%6] : oa[0],
                    (i >= 4) ? oa[(i-4)%6] : oa[0],
                    (i >= 5) ? oa[(i-5)%6] : oa[0], i);
      oat_store_row_fs_sse(om_fs, ox, oa[i%6], oxmx, pp, ppx, i);

      /* posterior sums for null2, added up as p7_Null2_fs_ByExpectation_sse() does */
      if (i == 1)
        {
          memcpy(cnt, pp, sizeof(__m128) * p7X_NSCELLS_FS * Q);
          cntx[p7X_N] = ppx[p7X_N];
          cntx[p7X_C] = ppx[p7X_C];
          cntx[p7X_J] = ppx[p7X_J];
        }
      else
        {
          for (q = 0; q < Q; q++)
            {
              MMO_FS(cnt, q, p7X_FS_C0) = _mm_add_ps(MMO_FS(pp, q, p7X_FS_C0), MMO_FS(cnt, q, p7X_FS_C0));
              IMO_FS(cnt, q)             = _mm_add_ps(IMO_FS(pp, q),             IMO_FS(cnt, q));
            }
          cntx[p7X_N] += ppx[p7X_N];
          cntx[p7X_C] += ppx[p7X_C];
          cntx[p7X_J] += ppx[p7X_J];
        }
    }

  /* Forward score, as p7_Forward_Frameshift_sse() computes it */
  xCtot = rg.xC[   L    %PARSER_ROWS_FWD]
        + rg.xC[((L-1)%PARSER_ROWS_FWD+PARSER_ROWS_FWD)%PARSER_ROWS_FWD]*om_fs->xf[p7O_C][p7O_LOOP]
        + rg.xC[((L-2)%PARSER_ROWS_FWD+PARSER_ROWS_FWD)%PARSER_ROWS_FWD]*om_fs->xf[p7O_C][p7O_LOOP];
  if (isnan(xCtot) || isinf(xCtot) == 1 || (L > 1 && xCtot == 0.0f)) { status = eslERANGE; goto ERROR; }

  if (opt_null2 != NULL) p7_Null2_fs_FromCounts_sse(om_fs, cnt, cntx, opt_null2);
  if (opt_fsc   != NULL) *opt_fsc = totscale + logf(xCtot*om_fs->xf[p7O_C][p7O_MOVE]);
  *ret_oasc = oxmx[L*p7X_NXCELLS+p7X_C] + oxmx[(L-1)*p7X_NXCELLS+p7X_C] + oxmx[(L-2)*p7X_NXCELLS+p7X_C];

  p7_omx_Destroy_sse(wk);
  p7_oivx_Destroy_sse(bv);
  free(oxmx); free(log_sfwd); free(log_sbck);
  return eslOK;

 ERROR:
  if (wk) p7_omx_Destroy_sse(wk);
  if (bv) p7_oivx_Destroy_sse(bv);
  free(oxmx); free(log_sfwd); free(log_sbck);
  return status;
}


/* Function:  p7_OATraceSweep_Frameshift_sse()
 * Synopsis:  Frameshift OA trace from what the SSE OA sweep stored.
 *
 * Purpose:   The frameshift OA trace for profile <om_fs> from <ox>, filled
 *            by <p7_OptimalAccuracySweep_Frameshift_sse()>, into the empty
 *            trace <tr> (created with <p7_trace_fs_CreateWithPP()>). The
 *            choices are those <p7_OATrace_Frameshift_sse()> makes on the
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
p7_OATraceSweep_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const P7_OATMX *ox, P7_TRACE *tr)
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

#endif /* eslENABLE_SSE */
