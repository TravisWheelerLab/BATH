/* Optimal accuracy alignment; NEON version.
 *
 * Contents:
 *   1. Optimal accuracy alignment, DP fill
 *   2. OA traceback
 *   3. Benchmark driver
 *   4. Unit tests
 *   5. Test driver
 *   6. Example
 *
 * SRE, Mon Aug 18 20:01:01 2008 [Casa de Gatos]
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
 * 1. Optimal accuracy alignment, DP fill
 *****************************************************************/

/* oa_row()
 * OA row <i> into <dpc> and row i of the OA special states <xmx>, from
 * OA row i-1 in <dpp> and the posterior row <ppp> with its special
 * states <ppx> (E, N, J, B, C of row i). Row 0 is the initialization
 * (<ppp>, <ppx>, <dpp> unused). <dpc> may be <ppp>: each posterior
 * vector is read before the OA vector at its place is written. Same
 * calculation as each row of p7_OptimalAccuracy(), which calls it.
 */
static void
oa_row(const P7_OPROFILE *om, const float32x4_t *ppp, const float *ppx, float *xmx, float32x4_t *dpc, const float32x4_t *dpp, int i)
{
  register float32x4_t mpv, dpv, ipv;   /* previous row values                                       */
  register float32x4_t sv;                 /* temp storage of 1 curr row value in progress              */
  register float32x4_t xEv;                /* E state: keeps max for Mk->E as we go                     */
  register float32x4_t xBv;                /* B state: splatted vector of B[i-1] for B->Mk calculations */
  register float32x4_t dcv;
  float32x4_t *tp;                         /* quads in the <om->tfv> transition scores                  */
  float32x4_t zerov = vmovq_n_f32(0.0);
  float32x4_t infv  = vmovq_n_f32(-eslINFINITY);
  int M = om->M;
  int Q = p7O_NQF(M);
  int q;
  int j;
  float t1, t2;


  if (i == 0)
    {
      for (q = 0; q < Q; q++) MMO(dpc, q) = IMO(dpc,q) = DMO(dpc,q) = infv;
      XMXo(0, p7X_E)    = -eslINFINITY;
      XMXo(0, p7X_N)    = 0.;
      XMXo(0, p7X_J)    = -eslINFINITY;
      XMXo(0, p7X_B)    = 0.;
      XMXo(0, p7X_C)    = -eslINFINITY;
      return;
    }

    tp  = om->tfv;            /* transition probabilities */
    dcv = infv;
    xEv = infv;
    xBv = vmovq_n_f32(XMXo(i-1, p7X_B));

    mpv = esl_neon_rightshift_float((esl_neon_128f_t) MMO(dpp,Q-1), (esl_neon_128f_t) infv).f32x4;  /* Right shifts by 4 bytes. 4,8,12,x becomes x,4,8,12. */
    dpv = esl_neon_rightshift_float((esl_neon_128f_t) DMO(dpp,Q-1), (esl_neon_128f_t) infv).f32x4;
    ipv = esl_neon_rightshift_float((esl_neon_128f_t) IMO(dpp,Q-1), (esl_neon_128f_t) infv).f32x4;
    for (q = 0; q < Q; q++)
      {
        sv  =               vreinterpretq_f32_u32(vandq_u32(vcgtq_f32(*tp, zerov), vreinterpretq_u32_f32(xBv)));  tp++;
  sv  = vmaxq_f32(sv, vreinterpretq_f32_u32(vandq_u32(vcgtq_f32(*tp, zerov), vreinterpretq_u32_f32(mpv)))); tp++;
  sv  = vmaxq_f32(sv, vreinterpretq_f32_u32(vandq_u32(vcgtq_f32(*tp, zerov), vreinterpretq_u32_f32(ipv)))); tp++;
        sv  = vmaxq_f32(sv, vreinterpretq_f32_u32(vandq_u32(vcgtq_f32(*tp, zerov), vreinterpretq_u32_f32(dpv)))); tp++;
        sv  = vaddq_f32(sv, *ppp);                                                                                ppp += 2;
        xEv = vmaxq_f32(xEv, sv);

        mpv = MMO(dpp,q);
        dpv = DMO(dpp,q);
        ipv = IMO(dpp,q);

        MMO(dpc,q) = sv;
        DMO(dpc,q) = dcv;

        dcv = vreinterpretq_f32_u32(vandq_u32(vcgtq_f32(*tp, zerov), vreinterpretq_u32_f32(sv))); tp++;

        sv         =               vreinterpretq_f32_u32(vandq_u32(vcgtq_f32(*tp, zerov), vreinterpretq_u32_f32(mpv)));   tp++;
        sv         = vmaxq_f32(sv, vreinterpretq_f32_u32(vandq_u32(vcgtq_f32(*tp, zerov), vreinterpretq_u32_f32(ipv))));  tp++;
        IMO(dpc,q) = vaddq_f32(sv, *ppp);                                   ppp++;
      }

    /* dcv has carried through from end of q loop above; store it
     * in first pass, we add M->D and D->D path into DMX
     */
    dcv = esl_neon_rightshift_float((esl_neon_128f_t) dcv, (esl_neon_128f_t) infv).f32x4;
    tp  = om->tfv + 7*Q;      /* set tp to start of the DD's */
    for (q = 0; q < Q; q++)
      {
        DMO(dpc, q) = vmaxq_f32(dcv, DMO(dpc, q));
        dcv         = vreinterpretq_f32_u32(vandq_u32(vcgtq_f32(*tp, zerov), vreinterpretq_u32_f32(DMO(dpc,q))));   tp++;
      }

    /* fully serialized D->D; can optimize later */
    for (j = 1; j < 4; j++)
      {
        dcv = esl_neon_rightshift_float((esl_neon_128f_t) dcv, (esl_neon_128f_t) infv).f32x4;
        tp  = om->tfv + 7*Q;
        for (q = 0; q < Q; q++)
          {
            DMO(dpc, q) = vmaxq_f32(dcv, DMO(dpc, q));
            dcv         = vreinterpretq_f32_u32(vandq_u32(vcgtq_f32(*tp, zerov), vreinterpretq_u32_f32(dcv)));   tp++;
          }
      }

    /* D->E paths */
    for (q = 0; q < Q; q++) xEv = vmaxq_f32(xEv, DMO(dpc,q));

    /* Specials */
    XMXo(i,p7X_E) = esl_neon_hmax_f32((esl_neon_128f_t) xEv);

    t1 = ( (om->xf[p7O_J][p7O_LOOP] == 0.0) ? 0.0 : xmx[(i-1)*p7X_NXCELLS+p7X_J] + ppx[p7X_J]);
    t2 = ( (om->xf[p7O_E][p7O_LOOP] == 0.0) ? 0.0 : xmx[   i *p7X_NXCELLS+p7X_E]);
    xmx[i*p7X_NXCELLS+p7X_J] = ESL_MAX(t1, t2);

    t1 = ( (om->xf[p7O_C][p7O_LOOP] == 0.0) ? 0.0 : xmx[(i-1)*p7X_NXCELLS+p7X_C] + ppx[p7X_C]);
    t2 = ( (om->xf[p7O_E][p7O_MOVE] == 0.0) ? 0.0 : xmx[   i *p7X_NXCELLS+p7X_E]);
    xmx[i*p7X_NXCELLS+p7X_C] = ESL_MAX(t1, t2);

    xmx[i*p7X_NXCELLS+p7X_N] = ((om->xf[p7O_N][p7O_LOOP] == 0.0) ? 0.0 : xmx[(i-1)*p7X_NXCELLS+p7X_N] + ppx[p7X_N]);

    t1 = ( (om->xf[p7O_N][p7O_MOVE] == 0.0) ? 0.0 : xmx[i*p7X_NXCELLS+p7X_N]);
    t2 = ( (om->xf[p7O_J][p7O_MOVE] == 0.0) ? 0.0 : xmx[i*p7X_NXCELLS+p7X_J]);
    xmx[i*p7X_NXCELLS+p7X_B] = ESL_MAX(t1, t2);
}


/* Function:  p7_OptimalAccuracy()
 * Synopsis:  DP fill of an optimal accuracy alignment calculation.
 * Incept:    SRE, Mon Aug 18 11:04:48 2008 [Janelia]
 *
 * Purpose:   Calculates the fill step of the optimal accuracy decoding
 *            algorithm \citep{Kall05}.
 *
 *            Caller provides the posterior decoding matrix <pp>,
 *            which was calculated by Forward/Backward on a target sequence
 *            of length <pp->L> using the query model <om>.
 *
 *            Caller also provides a DP matrix <ox>, allocated for a full
 *            <om->M> by <L> comparison. The routine fills this in
 *            with OA scores. <ox> may be <pp> itself: each cell's
 *            posterior is read before the cell's OA score is written,
 *            and no other posterior is read, so the OA matrix can
 *            overwrite the posterior matrix.
 *
 * Args:      om    - query profile
 *            pp    - posterior decoding matrix created by <p7_Decoding()>
 *            ox    - RESULT: caller provided DP matrix for <om->M> by <L>; may be <pp>
 *            ret_e - RETURN: expected number of correctly decoded positions
 *
 * Returns:   <eslOK> on success, and <*ret_e> contains the final OA
 *            score, which is the expected number of correctly decoded
 *            positions in the target sequence (up to <L>).
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_OptimalAccuracy(const P7_OPROFILE *om, const P7_OMX *pp, P7_OMX *ox, float *ret_e)
{
  int i;

  ox->M = om->M;
  ox->L = pp->L;
  oa_row(om, NULL, NULL, ox->xmx, ox->dpf[0], NULL, 0);
  for (i = 1; i <= pp->L; i++)
    oa_row(om, pp->dpf[i], pp->xmx + i*p7X_NXCELLS, ox->xmx, ox->dpf[i], ox->dpf[i-1], i);

  *ret_e = ox->xmx[pp->L*p7X_NXCELLS+p7X_C];
  return eslOK;
}
/*------------------- end, OA DP fill ---------------------------*/





/*****************************************************************
 * 2. OA traceback
 *****************************************************************/

static inline float get_postprob(const P7_OPROFILE *om, const P7_OMX *ox, int scur, int sprv, int k, int i);

static inline int select_m(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int k);
static inline int select_d(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int k);
static inline int select_i(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int k);
static inline int select_n(int i);
static inline int select_c(const P7_OPROFILE *om,                   const P7_OMX *ox, int i);
static inline int select_j(const P7_OPROFILE *om,                   const P7_OMX *ox, int i);
static inline int select_e(const P7_OPROFILE *om,                   const P7_OMX *ox, int i, int *ret_k);
static inline int select_b(const P7_OPROFILE *om,                   const P7_OMX *ox, int i);


/* Function:  p7_OATrace()
 * Synopsis:  Optimal accuracy decoding: traceback.
 * Incept:    SRE, Mon Aug 18 13:53:33 2008 [Janelia]
 *
 * Purpose:   The traceback stage of the optimal accuracy decoding algorithm
 *            \citep{Kall05}.
 *
 *            Caller provides the OA DP matrix <ox> that was just
 *            calculated by <p7_OptimalAccuracy()> for a target
 *            sequence and the query model <om>. The traceback needs
 *            only <ox>, so the OA matrix may have been written over the
 *            posterior decoding matrix, and the target sequence itself
 *            need not be provided. The posterior probability of each
 *            residue's state is recovered from <ox> as the difference
 *            between the OA score of its cell and the OA score of the
 *            predecessor that the fill chose, which equals the
 *            posterior up to float rounding.
 *
 *            The resulting optimal accuracy decoding traceback is put
 *            in a caller-provided traceback structure <tr>, which the
 *            caller has allocated for optional posterior probability
 *            annotation on residues (with <p7_trace_CreateWithPP()>,
 *            generally). This structure will be reallocated
 *            internally if necessary.
 *
 * Args:      om  - profile
 *            ox  - OA matrix to trace, LxM
 *            tr  - storage for the recovered traceback
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 *            <eslEINVAL> if the trace <tr> isn't empty (needs to be Reuse()'d).
 */
int
p7_OATrace(const P7_OPROFILE *om, const P7_OMX *ox, P7_TRACE *tr)
{
  int   i   = ox->L;		/* position in sequence 1..L */
  int   k   = 0;		/* position in model 1..M */
  int   s0, s1;			/* choice of a state */
  float postprob;
  int   status;

  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace not empty; needs to be Reuse()'d?");

  if ((status = p7_trace_AppendWithPP(tr, p7T_T, k, i, 0.0)) != eslOK) return status;
  if ((status = p7_trace_AppendWithPP(tr, p7T_C, k, i, 0.0)) != eslOK) return status;

  s0 = tr->st[tr->N-1];
  while (s0 != p7T_S)
    {
      switch (s0) {
      case p7T_M: s1 = select_m(om,     ox, i, k);  k--; i--; break;
      case p7T_D: s1 = select_d(om,     ox, i, k);  k--;      break;
      case p7T_I: s1 = select_i(om,     ox, i, k);       i--; break;
      case p7T_N: s1 = select_n(i);                           break;
      case p7T_C: s1 = select_c(om,     ox, i);               break;
      case p7T_J: s1 = select_j(om,     ox, i);               break;
      case p7T_E: s1 = select_e(om,     ox, i, &k);           break;
      case p7T_B: s1 = select_b(om,     ox, i);               break;
      default: ESL_EXCEPTION(eslEINVAL, "bogus state in traceback");
      }
      if (s1 == -1) ESL_EXCEPTION(eslEINVAL, "OA traceback choice failed");

      postprob = get_postprob(om, ox, s1, s0, k, i);
      if ((status = p7_trace_AppendWithPP(tr, s1, k, i, postprob)) != eslOK) return status;

      if ( (s1 == p7T_N || s1 == p7T_J || s1 == p7T_C) && s1 == s0) i--;
      s0 = s1;
    } /* end traceback, at S state */
  tr->M = om->M;
  tr->L = ox->L;
  return p7_trace_Reverse(tr);
}

/* get_postprob()
 * Posterior probability of state <scur> at (<i>,<k>), recovered from
 * the OA matrix <ox>. The fill set an M or I cell to its best
 * predecessor score plus the cell's posterior, with each predecessor
 * term set to 0 when its transition is not allowed and the M, I, D
 * cells of column 0 taken as -inf; the posterior is the cell's score
 * minus that same predecessor score. The N, C and J loops add the
 * posterior to the previous row's score.
 */
static inline float
get_postprob(const P7_OPROFILE *om, const P7_OMX *ox, int scur, int sprv, int k, int i)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;		/* (q,r) is the position of cell (i,k)    */
  int     r     = (k-1) / Q;
  int     qp    = (k-2) % Q;		/* (qp,rp) is the position of cell (i,k-1) */
  int     rp    = (k-2) / Q;
  float32x4_t  *tp    = om->tfv + 7*q;
  union { float32x4_t v; float p[4]; } u, tv;
  float   mv, iv, dv, prv;

  switch (scur) {
  case p7T_M:
    if (k > 1) {
      u.v = ox->dpf[i-1][qp*3 + p7X_M];  mv = u.p[rp];
      u.v = ox->dpf[i-1][qp*3 + p7X_I];  iv = u.p[rp];
      u.v = ox->dpf[i-1][qp*3 + p7X_D];  dv = u.p[rp];
    } else {
      mv = iv = dv = -eslINFINITY;
    }
    tv.v = tp[p7O_BM];  prv =              ((tv.p[r] > 0.0) ? ox->xmx[(i-1)*p7X_NXCELLS+p7X_B] : 0.0);
    tv.v = tp[p7O_MM];  prv = ESL_MAX(prv, ((tv.p[r] > 0.0) ? mv : 0.0));
    tv.v = tp[p7O_IM];  prv = ESL_MAX(prv, ((tv.p[r] > 0.0) ? iv : 0.0));
    tv.v = tp[p7O_DM];  prv = ESL_MAX(prv, ((tv.p[r] > 0.0) ? dv : 0.0));
    u.v  = ox->dpf[i][q*3 + p7X_M];
    return u.p[r] - prv;
  case p7T_I:
    u.v = ox->dpf[i-1][q*3 + p7X_M];  tv.v = tp[p7O_MI];  prv =              ((tv.p[r] > 0.0) ? u.p[r] : 0.0);
    u.v = ox->dpf[i-1][q*3 + p7X_I];  tv.v = tp[p7O_II];  prv = ESL_MAX(prv, ((tv.p[r] > 0.0) ? u.p[r] : 0.0));
    u.v = ox->dpf[i][q*3 + p7X_I];
    return u.p[r] - prv;
  case p7T_N: if (sprv == scur) return ox->xmx[i*p7X_NXCELLS+p7X_N] - ox->xmx[(i-1)*p7X_NXCELLS+p7X_N];
  case p7T_C: if (sprv == scur) return ox->xmx[i*p7X_NXCELLS+p7X_C] - ox->xmx[(i-1)*p7X_NXCELLS+p7X_C];
  case p7T_J: if (sprv == scur) return ox->xmx[i*p7X_NXCELLS+p7X_J] - ox->xmx[(i-1)*p7X_NXCELLS+p7X_J];
  default:    return 0.0;
  }
}

/* M(i,k) is reached from B(i-1), M(i-1,k-1), D(i-1,k-1), or I(i-1,k-1). */
static inline int
select_m(const P7_OPROFILE *om, const P7_OMX *ox, int i, int k)
{
  int          Q     = p7O_NQF(ox->M);
  int          q     = (k-1) % Q;		/* (q,r) is position of the current DP cell M(i,k) */
  int          r     = (k-1) / Q;
  float32x4_t *tp    = om->tfv + 7*q;       	/* *tp now at start of transitions to cur cell M(i,k) */
  float32x4_t  xBv   = vmovq_n_f32(ox->xmx[(i-1)*p7X_NXCELLS+p7X_B]);
  float32x4_t  zerov = vmovq_n_f32(0.0);
  float32x4_t  mpv, dpv, ipv;
  union { float32x4_t v; float p[4]; } u, tv;
  float   path[4];
  int     state[4] = { p7T_M, p7T_I, p7T_D, p7T_B };

  if (q > 0) {
    mpv = ox->dpf[i-1][(q-1)*3 + p7X_M];
    dpv = ox->dpf[i-1][(q-1)*3 + p7X_D];
    ipv = ox->dpf[i-1][(q-1)*3 + p7X_I];
  } else {
    mpv = esl_neon_rightshift_float((esl_neon_128f_t) ox->dpf[i-1][(Q-1)*3 + p7X_M], (esl_neon_128f_t) zerov).f32x4;
    dpv = esl_neon_rightshift_float((esl_neon_128f_t) ox->dpf[i-1][(Q-1)*3 + p7X_D], (esl_neon_128f_t) zerov).f32x4;
    ipv = esl_neon_rightshift_float((esl_neon_128f_t) ox->dpf[i-1][(Q-1)*3 + p7X_I], (esl_neon_128f_t) zerov).f32x4;
  }

  /* paths are numbered so that most desirable choice in case of tie is first. */
  u.v = xBv;  tv.v = *tp;  path[3] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);  tp++;
  u.v = mpv;  tv.v = *tp;  path[0] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);  tp++;
  u.v = ipv;  tv.v = *tp;  path[1] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);  tp++;
  u.v = dpv;  tv.v = *tp;  path[2] = ((tv.p[r] == 0.0) ?  -eslINFINITY : u.p[r]);
  return state[esl_vec_FArgMax(path, 4)];
}


/* D(i,k) is reached from M(i, k-1) or D(i,k-1). */
static inline int
select_d(const P7_OPROFILE *om, const P7_OMX *ox, int i, int k)
{
  int     Q     = p7O_NQF(ox->M);
  int     q     = (k-1) % Q;		/* (q,r) is position of the current DP cell D(i,k) */
  int     r     = (k-1) / Q;
  float32x4_t  zerov = vmovq_n_f32(0.0);
  union { float32x4_t v; float p[4]; } mpv, dpv, tmdv, tddv;
  float   path[2];

  if (q > 0) {
    mpv.v  = ox->dpf[i][(q-1)*3 + p7X_M];
    dpv.v  = ox->dpf[i][(q-1)*3 + p7X_D];
    tmdv.v = om->tfv[7*(q-1) + p7O_MD];
    tddv.v = om->tfv[7*Q + (q-1)];
  } else {
    mpv.v  = esl_neon_rightshift_float((esl_neon_128f_t) ox->dpf[i][(Q-1)*3 + p7X_M], (esl_neon_128f_t) zerov).f32x4;
    dpv.v  = esl_neon_rightshift_float((esl_neon_128f_t) ox->dpf[i][(Q-1)*3 + p7X_D], (esl_neon_128f_t) zerov).f32x4;
    tmdv.v = esl_neon_rightshift_float((esl_neon_128f_t) om->tfv[7*(Q-1) + p7O_MD],   (esl_neon_128f_t) zerov).f32x4;
    tddv.v = esl_neon_rightshift_float((esl_neon_128f_t) om->tfv[8*Q-1],              (esl_neon_128f_t) zerov).f32x4;
  }

  path[0] = ((tmdv.p[r] == 0.0) ? -eslINFINITY : mpv.p[r]);
  path[1] = ((tddv.p[r] == 0.0) ? -eslINFINITY : dpv.p[r]);
  return  ((path[0] >= path[1]) ? p7T_M : p7T_D);
}

/* I(i,k) is reached from M(i-1, k) or I(i-1,k). */
static inline int
select_i(const P7_OPROFILE *om, const P7_OMX *ox, int i, int k)
{
  int          Q   = p7O_NQF(ox->M);
  int          q   = (k-1) % Q;		/* (q,r) is position of the current DP cell D(i,k) */
  int          r   = (k-1) / Q;
  float32x4_t *tp  = om->tfv + 7*q + p7O_MI;
  union { float32x4_t v; float p[4]; } tv, mpv, ipv;
  float   path[2];

  mpv.v = ox->dpf[i-1][q*3 + p7X_M]; tv.v = *tp;  path[0] = ((tv.p[r] == 0.0) ? -eslINFINITY : mpv.p[r]);  tp++;
  ipv.v = ox->dpf[i-1][q*3 + p7X_I]; tv.v = *tp;  path[1] = ((tv.p[r] == 0.0) ? -eslINFINITY : ipv.p[r]);
  return  ((path[0] >= path[1]) ? p7T_M : p7T_I);
}

/* N(i) must come from N(i-1) for i>0; else it comes from S */
static inline int
select_n(int i)
{
  return ((i==0) ? p7T_S : p7T_N);
}

/* C(i) is reached from E(i) or C(i-1). */
static inline int
select_c(const P7_OPROFILE *om, const P7_OMX *ox, int i)
{
  /* The fill set C(i) to the larger of C(i-1) + pp_C(i) and E(i), each
   * term 0 when its transition is not allowed. The C->C loop is
   * taken when it is allowed and either E->C is not allowed (and
   * C(i-1) is reachable) or the loop term is larger, which is
   * C(i) > E(i); on a tie E wins, as in the fill.
   */
  if (om->xf[p7O_C][p7O_LOOP]  == 0.0) return p7T_E;
  if (om->xf[p7O_E][p7O_MOVE] == 0.0) return ((ox->xmx[(i-1)*p7X_NXCELLS+p7X_C] > -eslINFINITY) ? p7T_C : p7T_E);
  return ((ox->xmx[i*p7X_NXCELLS+p7X_C] > ox->xmx[i*p7X_NXCELLS+p7X_E]) ? p7T_C : p7T_E);
}

/* J(i) is reached from E(i) or J(i-1). */
static inline int
select_j(const P7_OPROFILE *om, const P7_OMX *ox, int i)
{
  /* The fill set J(i) to the larger of J(i-1) + pp_J(i) and E(i), each
   * term 0 when its transition is not allowed. The J->J loop is
   * taken when it is allowed and either E->J is not allowed (and
   * J(i-1) is reachable) or the loop term is larger, which is
   * J(i) > E(i); on a tie E wins, as in the fill.
   */
  if (om->xf[p7O_J][p7O_LOOP]  == 0.0) return p7T_E;
  if (om->xf[p7O_E][p7O_LOOP] == 0.0) return ((ox->xmx[(i-1)*p7X_NXCELLS+p7X_J] > -eslINFINITY) ? p7T_J : p7T_E);
  return ((ox->xmx[i*p7X_NXCELLS+p7X_J] > ox->xmx[i*p7X_NXCELLS+p7X_E]) ? p7T_J : p7T_E);
}

/* E(i) is reached from any M(i, k=1..M) or D(i, k=2..M). */
/* This assumes all M_k->E, D_k->E are 1.0 */
static inline int
select_e(const P7_OPROFILE *om, const P7_OMX *ox, int i, int *ret_k)
{
  int     Q          = p7O_NQF(ox->M);
  float32x4_t *dp    = ox->dpf[i];
  union { float32x4_t v; float p[4]; } u;
  float  max   = -eslINFINITY;
  int    smax, kmax;
  int    q,r;

  /* precedence rules in case of ties here are a little tricky: M beats D: note the >= max!  */
  for (q = 0; q < Q; q++)
    {
      u.v   = *dp; dp++;  for (r = 0; r < 4; r++) if (u.p[r] >= max) { max = u.p[r]; smax = p7T_M; kmax = r*Q + q + 1; }
      u.v   = *dp; dp+=2; for (r = 0; r < 4; r++) if (u.p[r] > max)  { max = u.p[r]; smax = p7T_D; kmax = r*Q + q + 1; }
    }
  *ret_k = kmax;
  return smax;
}


/* B(i) is reached from N(i) or J(i). */
static inline int
select_b(const P7_OPROFILE *om, const P7_OMX *ox, int i)
{
  float path[2];

  path[0] = ( (om->xf[p7O_N][p7O_MOVE] == 0.0) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_N]);
  path[1] = ( (om->xf[p7O_J][p7O_MOVE] == 0.0) ? -eslINFINITY : ox->xmx[i*p7X_NXCELLS+p7X_J]);
  return  ((path[0] > path[1]) ? p7T_N : p7T_J);
}
/*---------------------- end, OA traceback ----------------------*/





/*****************************************************************
 * OA sweep over a checkpointed Backward, and its trace
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
  oa_row(om, NULL, NULL, oxmx, wk->dpf[2], NULL, 0);
  oat_store_row(om, ox, wk->dpf[2], oxmx, NULL, NULL, 0);

  for (i = 1; i <= L; i++)
    {
      b = (i / cx->blockR) * cx->blockR;
      if (b != cx->b && (status = p7_BackwardChkLoad(dsq, om, cx, b)) != eslOK) goto ERROR;

      p7_ForwardRow(dsq, om, wk->xmx, &totscale, wk->dpf[i%2], wk->dpf[(i-1)%2], i);
      p7_DecodingRow(om, wk->xmx, cx->xmx, wk->dpf[i%2], cx->dpf[i], pp, ppx, i, scaleproduct);
      scaleproduct *= wk->xmx[i*p7X_NXCELLS+p7X_SCALE] / cx->xmx[i*p7X_NXCELLS+p7X_SCALE];

      oa_row(om, pp, ppx, oxmx, wk->dpf[2 + i%2], wk->dpf[2 + (i-1)%2], i);
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

/*****************************************************************
 * 3. Benchmark driver
 *****************************************************************/
#ifdef p7OPTACC_BENCHMARK
/*
   icc  -O3 -static -o optacc_benchmark -I.. -L.. -I../../easel -L../../easel -Dp7OPTACC_BENCHMARK optacc.c -lhmmer -leasel -lm

   ./optacc_benchmark <hmmfile>         runs benchmark on optimal accuracy fill and trace
   ./optacc_benchmark -c -N1 <hmmfile>  compare scores of NEON version to generic impl
   ./optacc_benchmark -x -N1 <hmmfile>  test that scores match trusted implementation.

                    RRM_1 (M=72)       Caudal_act (M=136)     SMC_N (M=1151)
                 -----------------    ------------------     ---------------
   20 Aug 08:     13.11u (110 Mc/s)     23.39u (116 Mc/s)    332.62u (69 Mc/s)

 */
#include <p7_config.h>

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
  { "-c",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-x", "compare scores to generic implementation (debug)", 0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-x",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-c", "equate scores to trusted implementation (debug)",  0 },
  { "-L",        eslARG_INT,    "400", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                     0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  { "--notrace", eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark the DP fill stage",                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for optimal accuracy alignment, NEON version";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_GMX         *gx1     = NULL;
  P7_GMX         *gx2     = NULL;
  P7_OMX         *ox1     = NULL;
  P7_OMX         *ox2     = NULL;
  P7_TRACE       *tr      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           fsc, bsc, accscore;
  float           fsc_g, bsc_g, accscore_g;
  double          Mcs;

  p7_FLogsumInit();

  if (p7_hmmfile_Open(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)           != eslOK) p7_Fail("Failed to read HMM");

  bg = p7_bg_Create(abc);                 p7_bg_SetLength(bg, L);
  gm = p7_profile_Create(hmm->M, abc);    p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL);
  om = p7_oprofile_Create(gm->M, abc);    p7_oprofile_Convert(gm, om);
  p7_oprofile_ReconfigLength(om, L);

  if (esl_opt_GetBoolean(go, "-x") && p7_FLogsumError(-0.4, -0.5) > 0.0001)
    p7_Fail("-x here requires p7_Logsum() recompiled in slow exact mode");

  ox1 = p7_omx_Create(gm->M, L, L);
  ox2 = p7_omx_Create(gm->M, L, L);
  tr  = p7_trace_CreateWithPP();

  esl_rsq_xfIID(r, bg->f, abc->K, L, dsq);
  p7_Forward (dsq, L, om, ox1,      &fsc);
  p7_Backward(dsq, L, om, ox1, ox2, &bsc);
  p7_Decoding(om, ox1, ox2, ox2);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      p7_OptimalAccuracy(om, ox2, ox1, &accscore);

      if (! esl_opt_GetBoolean(go, "--notrace"))
	{
	  p7_OATrace(om, ox1, tr);
	  p7_trace_Reuse(tr);
	}
    }
  esl_stopwatch_Stop(w);

  Mcs        = (double) N * (double) L * (double) gm->M * 1e-6 / (double) w->user;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm->M);
  printf("# %.1f Mc/s\n", Mcs);

  if (esl_opt_GetBoolean(go, "-c") || esl_opt_GetBoolean(go, "-x") )
    {
      gx1 = p7_gmx_Create(gm->M, L, L, p7G_NSCELLS);
      gx2 = p7_gmx_Create(gm->M, L, L, p7G_NSCELLS);

      p7_GForward (dsq, L, gm, gx1, &fsc_g);
      p7_GBackward(dsq, L, gm, gx2, &bsc_g);
      p7_GDecoding(gm, gx1, gx2, gx2);
      p7_GOptimalAccuracy(gm, gx2, gx1, &accscore_g);

      printf("generic:  fwd=%8.4f  bck=%8.4f  acc=%8.4f\n", fsc_g, bsc_g, accscore_g);
      printf("NEON:      fwd=%8.4f  bck=%8.4f  acc=%8.4f\n", fsc,   bsc,   accscore);

      p7_gmx_Destroy(gx1);
      p7_gmx_Destroy(gx2);
    }

  free(dsq);
  p7_omx_Destroy(ox1);
  p7_omx_Destroy(ox2);
  p7_trace_Destroy(tr);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_alphabet_Destroy(abc);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7OPTACC_BENCHMARK*/
/*---------------- end, benchmark driver ------------------------*/




/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7OPTACC_TESTDRIVE
#include <math.h>

#include "esl_randomseq.h"
#include "esl_vectorops.h"

#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
/*
 * 1. Compare accscore to GOptimalAccuracy().
 * 2. Compare trace to GOATrace().
 *
 * Note: This test is subject to some expected noise and can fail
 * for entirely innocent reasons. Generic Forward/Backward calculations with
 * p7_GForward(), p7_GBackward() use coarse-grain table lookups to sum
 * log probabilities, and sufficient roundoff error can accumulate to
 * change the optimal accuracy traceback, causing this test to fail.
 * So, if optacc_utest fails, before you go looking for bugs, first
 * go to ../logsum.c, change the #ifdef to activate the slow/accurate
 * version, recompile and rerun optacc_utest. If the failure goes away,
 * you can ignore it.   - SRE, Wed Dec 17 09:45:31 2008
 */
static void
utest_optacc(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char        *msg = "optimal accuracy unit test failed";
  P7_HMM      *hmm = NULL;
  P7_PROFILE  *gm  = NULL;
  P7_OPROFILE *om  = NULL;
  ESL_SQ      *sq  = esl_sq_CreateDigital(abc);
  P7_OMX      *ox1 = p7_omx_Create(M, L, L);
  P7_OMX      *ox2 = p7_omx_Create(M, L, L);
  P7_GMX      *gx1 = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX      *gx2 = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_TRACE    *tr  = p7_trace_CreateWithPP();
  P7_TRACE    *trg = p7_trace_CreateWithPP();
  P7_TRACE    *tro = p7_trace_CreateWithPP();
  float        accscore_o;
  float        fsc, bsc, accscore;
  float        fsc_g, bsc_g, accscore_g, accscore_g2;
  float        pptol = 0.01;
  float        sctol = 0.001;
  float        gtol;

  p7_FLogsumInit();
  gtol = ( (p7_FLogsumError(-0.4, -0.5) > 0.0001) ?  0.1 : 0.001);

  if (p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om)!= eslOK) esl_fatal(msg);
  while (N--)
    {
      if (p7_ProfileEmit(r, hmm, gm, bg, sq, tro)         != eslOK) esl_fatal(msg);

      if (p7_omx_GrowTo(ox1, M, sq->n, sq->n)             != eslOK) esl_fatal(msg);
      if (p7_omx_GrowTo(ox2, M, sq->n, sq->n)             != eslOK) esl_fatal(msg);
      if (p7_gmx_GrowTo(gx1, M, sq->n, sq->n)             != eslOK) esl_fatal(msg);
      if (p7_gmx_GrowTo(gx2, M, sq->n, sq->n)             != eslOK) esl_fatal(msg);

      if (p7_Forward (sq->dsq, sq->n, om, ox1,      &fsc) != eslOK) esl_fatal(msg);
      if (p7_Backward(sq->dsq, sq->n, om, ox1, ox2, &bsc) != eslOK) esl_fatal(msg);
      if (p7_Decoding(om, ox1, ox2, ox2)                  != eslOK) esl_fatal(msg);
      if (p7_OptimalAccuracy(om, ox2, ox1, &accscore)     != eslOK) esl_fatal(msg);

#if 0
      p7_omx_FDeconvert(ox1, gx1);
      p7_gmx_Dump(stdout, gx1, p7_DEFAULT);
      p7_omx_FDeconvert(ox2, gx1);
      p7_gmx_Dump(stdout, gx1, p7_DEFAULT);
#endif
      if (p7_OATrace(om, ox1, tr)                    != eslOK) esl_fatal(msg);

      if (p7_GForward (sq->dsq, sq->n, gm, gx1, &fsc_g)   != eslOK) esl_fatal(msg);
      if (p7_GBackward(sq->dsq, sq->n, gm, gx2, &bsc_g)   != eslOK) esl_fatal(msg);

#if 0
      p7_gmx_Dump(stdout, gx1, p7_DEFAULT); /* fwd */
      p7_gmx_Dump(stdout, gx2, p7_DEFAULT); /* bck */
#endif

      if (p7_GDecoding(gm, gx1, gx2, gx2)                 != eslOK) esl_fatal(msg);
      if (p7_GOptimalAccuracy(gm, gx2, gx1, &accscore_g)  != eslOK) esl_fatal(msg);

#if 0
      p7_gmx_Dump(stdout, gx1, p7_DEFAULT); /* oa */
      p7_gmx_Dump(stdout, gx2, p7_DEFAULT); /* pp */
#endif
      if (p7_GOATrace(gm, gx1, trg)                       != eslOK) esl_fatal(msg);

      if (p7_trace_SetPP(tro, gx2)                        != eslOK) esl_fatal(msg);

      if (esl_opt_GetBoolean(go, "--traces"))
	{
	  p7_trace_Dump(stdout, tro, gm, sq->dsq);
	  p7_trace_Dump(stdout, tr,  gm, sq->dsq);
	  p7_trace_Dump(stdout, trg, gm, sq->dsq);
	}

      if (p7_trace_Validate(tr,  abc, sq->dsq, NULL)      != eslOK) esl_fatal(msg);
      if (p7_trace_Validate(trg, abc, sq->dsq, NULL)      != eslOK) esl_fatal(msg);
      if (p7_trace_Compare(tr, trg, pptol)                != eslOK) esl_fatal(msg);

      accscore_o  = p7_trace_GetExpectedAccuracy(tro); /* according to gx2; see p7_trace_SetPP() call above */
      accscore_g2 = p7_trace_GetExpectedAccuracy(trg);

#if 0
      printf("%f %f %f %f\n", accscore, accscore_g, accscore_g2, accscore_o);
#endif

      if (esl_FCompare_old(fsc,        bsc,         sctol)    != eslOK) esl_fatal(msg);
      if (esl_FCompare_old(fsc_g,      bsc_g,       gtol)     != eslOK) esl_fatal(msg);
      if (esl_FCompare_old(fsc,        fsc_g,       gtol)     != eslOK) esl_fatal(msg);
      if (esl_FCompare_old(accscore,   accscore_g,  gtol)     != eslOK) esl_fatal(msg);
      if (esl_FCompare_old(accscore_g, accscore_g2, gtol)     != eslOK) esl_fatal(msg);
      if (accscore_g2 < accscore_o)                                 esl_fatal(msg);
      /* the above deserves explanation:
       *  - accscore_o is the accuracy of the originally emitted trace, according
       *      to the generic posterior decoding matrix <gx2>. This is a lower bound
       *      on the expected # of accurately aligned residues found by a DP
       *      optimization.
       *  - accscore is the accuracy found by the fast (vector) code DP implementation.
       *  - accscore_g is the accuracy found by the generic DP implementation.
       *      accscore and accscore_g should be nearly identical,
       *      within tolerance of roundoff error accumulation and
       *      the imprecision of Logsum() tables.
       *  - accscore_g2 is the accuracy of the traceback identified by the generic
       *      DP implementation. It should be identical (within order-of-evaluation
       *      roundoff error) to accscore_g.
       *
       * the "accscore_g2 < accscore_o" test is carefully contrived.
       * accscore_o is a theoretical lower bound but because of fp error,
       * accscore and (much more rarely) even accscore_g can exceed accscore_o.
       * accscore_g2, however, is calculated with identical order of evaluation
       * as accscore_o if the optimal trace does turn out to be identical to
       * the originally emitted trace. It should be extremely unlikely (though
       * not impossible) for accscore_o to exceed accscore_g2. (The DP algorithm
       * would have to identify a trace that was different than the original trace,
       * which the DP algorithm, by order-of-evaluation, assigned higher accuracy,
       * but order-of-evaluation in traceback dependent code assigned lower accuracy.
       * [xref J5/29]
       */

      esl_sq_Reuse(sq);
      p7_trace_Reuse(tr);
      p7_trace_Reuse(trg);
      p7_trace_Reuse(tro);
    }

  p7_trace_Destroy(tro);
  p7_trace_Destroy(trg);
  p7_trace_Destroy(tr);
  p7_gmx_Destroy(gx2);
  p7_gmx_Destroy(gx1);
  p7_omx_Destroy(ox2);
  p7_omx_Destroy(ox1);
  esl_sq_Destroy(sq);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}


/* utest_inplace()
 * On sequences emitted from a sampled profile, the OA fill written over
 * the posterior matrix (<ox> = <pp>) must give the same OA matrix and
 * score as a fill into a separate matrix, and the traceback from either
 * matrix must be the same. The traceback's posterior probabilities,
 * recovered from the OA matrix by subtraction, must agree with the
 * posterior matrix within <pptol>.
 */
static void
utest_inplace(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char        *msg   = "optimal accuracy in-place unit test failed";
  P7_HMM      *hmm   = NULL;
  P7_PROFILE  *gm    = NULL;
  P7_OPROFILE *om    = NULL;
  ESL_SQ      *sq    = esl_sq_CreateDigital(abc);
  P7_OMX      *ox1   = p7_omx_Create(M, L, L);
  P7_OMX      *ox2   = p7_omx_Create(M, L, L);
  P7_GMX      *gpp   = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX      *gox1  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX      *gox2  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_TRACE    *tro   = p7_trace_Create();
  P7_TRACE    *tr1   = p7_trace_CreateWithPP();
  P7_TRACE    *tr2   = p7_trace_CreateWithPP();
  float        pptol = 0.001;
  float        oasc1, oasc2;
  int          z;

  p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om);
  while (N--)
    {
      if (p7_ProfileEmit(r, hmm, gm, bg, sq, tro)                  != eslOK) esl_fatal(msg);
      if (p7_gmx_GrowTo(gpp,  M, sq->n, sq->n)                     != eslOK) esl_fatal(msg);
      if (p7_gmx_GrowTo(gox1, M, sq->n, sq->n)                     != eslOK) esl_fatal(msg);
      if (p7_gmx_GrowTo(gox2, M, sq->n, sq->n)                     != eslOK) esl_fatal(msg);
      if (p7_omx_GrowTo(ox1, M, sq->n, sq->n)                  != eslOK) esl_fatal(msg);
      if (p7_omx_GrowTo(ox2, M, sq->n, sq->n)                  != eslOK) esl_fatal(msg);
      if (p7_Forward (sq->dsq, sq->n, om, ox1,      NULL)      != eslOK) esl_fatal(msg);
      if (p7_Backward(sq->dsq, sq->n, om, ox1, ox2, NULL)      != eslOK) esl_fatal(msg);
      if (p7_Decoding(om, ox1, ox2, ox2)                       != eslOK) esl_fatal(msg);
      if (p7_omx_FDeconvert(ox2, gpp)                          != eslOK) esl_fatal(msg);
      if (p7_OptimalAccuracy(om, ox2, ox1, &oasc1)             != eslOK) esl_fatal(msg);
      if (p7_OATrace(om, ox1, tr1)                             != eslOK) esl_fatal(msg);
      if (p7_OptimalAccuracy(om, ox2, ox2, &oasc2)             != eslOK) esl_fatal(msg);
      if (p7_OATrace(om, ox2, tr2)                             != eslOK) esl_fatal(msg);
      if (p7_omx_FDeconvert(ox1, gox1)                         != eslOK) esl_fatal(msg);
      if (p7_omx_FDeconvert(ox2, gox2)                         != eslOK) esl_fatal(msg);

      if (oasc1 != oasc2)                                        esl_fatal("%s: OA scores differ", msg);
      if (p7_gmx_Compare(gox2, gox1, 0.0)              != eslOK) esl_fatal("%s: OA matrices differ", msg);
      if (tr1->N != tr2->N)                                      esl_fatal("%s: traces differ", msg);
      for (z = 0; z < tr1->N; z++)
        if (tr1->st[z] != tr2->st[z] || tr1->k[z] != tr2->k[z] || tr1->i[z] != tr2->i[z] || tr1->pp[z] != tr2->pp[z])
          esl_fatal("%s: traces differ", msg);
      for (z = 0; z < tr1->N; z++)
        {
          if      (tr1->st[z] == p7T_M) { if (fabs(tr1->pp[z] - gpp->dp[tr1->i[z]][tr1->k[z]*p7G_NSCELLS + p7G_M]) > pptol) esl_fatal("%s: M posterior", msg); }
          else if (tr1->st[z] == p7T_I) { if (fabs(tr1->pp[z] - gpp->dp[tr1->i[z]][tr1->k[z]*p7G_NSCELLS + p7G_I]) > pptol) esl_fatal("%s: I posterior", msg); }
          else if ((tr1->st[z] == p7T_N || tr1->st[z] == p7T_C || tr1->st[z] == p7T_J) && tr1->i[z] > 0 && tr1->st[z-1] == tr1->st[z])
            {
              if (fabs(tr1->pp[z] - gpp->xmx[tr1->i[z]*p7G_NXCELLS + (tr1->st[z] == p7T_N ? p7G_N : tr1->st[z] == p7T_C ? p7G_C : p7G_J)]) > pptol) esl_fatal("%s: NCJ posterior", msg);
            }
        }

      esl_sq_Reuse(sq);
      p7_trace_Reuse(tro);
      p7_trace_Reuse(tr1);
      p7_trace_Reuse(tr2);
    }

  p7_trace_Destroy(tro);
  p7_trace_Destroy(tr1);
  p7_trace_Destroy(tr2);
  p7_gmx_Destroy(gpp);
  p7_gmx_Destroy(gox1);
  p7_gmx_Destroy(gox2);
  p7_omx_Destroy(ox1);
  p7_omx_Destroy(ox2);
  p7_oprofile_Destroy(om);
  esl_sq_Destroy(sq);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}
/* utest_sweep()
 * The OA sweep over a checkpointed Backward against the full-matrix
 * path on the same Backward: p7_Forward(), the checkpointed Backward
 * filled in one block and copied into a Backward matrix with its own
 * scale factors, p7_Decoding(), p7_OptimalAccuracy(), p7_OATrace() and
 * p7_Null2_ByExpectation(). The sweep, with Backward blocks of 1..8
 * rows, must give the same Forward and OA scores, null2, states, k and
 * i on the trace, and PP characters equal to those the posterior
 * matrix gives. Sequences alternate random ones and the consensus of a
 * model with sharpened match emissions and M->M transitions between
 * random flanks, so that rows are rescaled.
 */
static void
utest_sweep(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N)
{
  char        *msg  = "OA sweep unit test failed";
  P7_HMM      *hmm  = NULL;
  P7_PROFILE  *gm   = NULL;
  P7_OPROFILE *om   = NULL;
  ESL_SQ      *sq   = esl_sq_CreateDigital(abc);
  ESL_DSQ     *dsq  = malloc(sizeof(ESL_DSQ) * (L+2));
  float       *n2a  = malloc(sizeof(float) * abc->Kp);
  float       *n2b  = malloc(sizeof(float) * abc->Kp);
  P7_OMX      *fwd  = NULL, *bck = NULL, *pp = NULL, *oa = NULL;
  P7_CHKOMX   *cx1  = NULL, *cx2 = NULL;
  P7_OATMX    *ox   = p7_oatmx_Create(M, L, p7OAT_NPLANES, 1);
  P7_TRACE    *tr1  = p7_trace_CreateWithPP();
  P7_TRACE    *tr2  = p7_trace_CreateWithPP();
  int          Q    = p7O_NQF(M);
  size_t       rowb = sizeof(float32x4_t) * 3 * (size_t) Q;
  float        fsc, fsc2, oasc1, oasc2, p;
  int          n, Ln, blockR, i, k, x, z;

  if (p7_hmm_Sample(r, M, abc, &hmm) != eslOK) esl_fatal(msg);
  for (k = 1; k <= M; k++)
    {
      x = esl_vec_FArgMax(hmm->mat[k], abc->K);
      esl_vec_FSet(hmm->mat[k], abc->K, 0.15 / (abc->K - 1));
      hmm->mat[k][x] = 0.85;
      if (k < M) { hmm->t[k][p7H_MM] = 0.9; hmm->t[k][p7H_MI] = 0.05; hmm->t[k][p7H_MD] = 0.05; }
    }
  gm = p7_profile_Create(M, abc);
  if (p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL) != eslOK) esl_fatal(msg);
  p7_emit_SimpleConsensus(hmm, sq);

  om  = p7_oprofile_Create(M, abc);  p7_oprofile_Convert(gm, om);
  fwd = p7_omx_Create(M, L, L); bck = p7_omx_Create(M, L, L); pp = p7_omx_Create(M, L, L); oa = p7_omx_Create(M, L, L);
  cx1 = p7_chkomx_Create(M, L, L+1, 1, p7X_NSCELLS, p7_CHK_BACKWARD);
  cx2 = p7_chkomx_Create(M, L, 1,   1, p7X_NSCELLS, p7_CHK_BACKWARD);

  for (n = 0; n < N; n++)
    {
      Ln = 1 + esl_rnd_Roll(r, L);
      esl_rsq_xfIID(r, bg->f, abc->K, Ln, dsq);
      if (n % 2 == 1 && M <= L)
        {
          Ln = ESL_MAX(Ln, M);
          esl_rsq_xfIID(r, bg->f, abc->K, Ln, dsq);
          memcpy(dsq + 1 + esl_rnd_Roll(r, Ln - M + 1), sq->dsq + 1, sizeof(ESL_DSQ) * M);
        }
      blockR = 1 + esl_rnd_Roll(r, 8);
      p7_trace_Reuse(tr1);
      p7_trace_Reuse(tr2);

      p7_oprofile_ReconfigLength(om, Ln);
      if (p7_Forward(dsq, Ln, om, fwd, &fsc)                    != eslOK) esl_fatal(msg);
      if (p7_chkomx_GrowTo(cx1, M, Ln, Ln+1)                    != eslOK) esl_fatal(msg);
      if (p7_BackwardChk(dsq, Ln, om, cx1, NULL)                != eslOK) esl_fatal(msg);
      for (i = 0; i <= Ln; i++) memcpy(bck->dpf[i], cx1->dpf[i], rowb);
      memcpy(bck->xmx, cx1->xmx, sizeof(float) * (Ln+1) * p7X_NXCELLS);
      bck->M = M; bck->L = Ln; bck->has_own_scales = TRUE;

      if (p7_Decoding(om, fwd, bck, pp)                         != eslOK) esl_fatal(msg);
      p7_OptimalAccuracy(om, pp, oa, &oasc1);
      if (p7_OATrace(om, oa, tr1)                               != eslOK) esl_fatal(msg);
      if (p7_chkomx_GrowTo(cx2, M, Ln, blockR)                  != eslOK) esl_fatal(msg);
      if (p7_BackwardChk(dsq, Ln, om, cx2, NULL)                != eslOK) esl_fatal(msg);
      if (p7_OptimalAccuracySweep(dsq, Ln, om, cx2, ox, n2b, &fsc2, &oasc2) != eslOK) esl_fatal(msg);
      if (p7_OATraceSweep(om, ox, tr2)                          != eslOK) esl_fatal(msg);

      if (fsc2  != fsc)    esl_fatal("%s: Forward scores differ (%f %f)", msg, fsc2, fsc);
      if (oasc2 != oasc1)  esl_fatal("%s: OA scores differ (%f %f)", msg, oasc2, oasc1);
      if (tr1->N != tr2->N) esl_fatal("%s: trace lengths differ (%d %d)", msg, tr1->N, tr2->N);
      for (z = 0; z < tr1->N; z++)
        {
          if (tr1->st[z] != tr2->st[z] || tr1->k[z] != tr2->k[z] || tr1->i[z] != tr2->i[z])
            esl_fatal("%s: traces differ at %d", msg, z);
          switch (tr1->st[z]) {
          case p7T_M: p = p7_omx_FGetMDI(pp, p7X_M, tr1->i[z], tr1->k[z]); break;
          case p7T_I: p = p7_omx_FGetMDI(pp, p7X_I, tr1->i[z], tr1->k[z]); break;
          case p7T_N: case p7T_C: case p7T_J:
            /* an N, C or J emits residue i when the state before it is the same one */
            p = (z > 0 && tr1->st[z-1] == tr1->st[z]) ? pp->xmx[tr1->i[z]*p7X_NXCELLS + (tr1->st[z] == p7T_N ? p7X_N : tr1->st[z] == p7T_C ? p7X_C : p7X_J)] : 0.0;
            break;
          default: p = 0.0; break;
          }
          if (p7_alidisplay_EncodePostProb(tr2->pp[z]) != p7_alidisplay_EncodePostProb(p))
            esl_fatal("%s: PP character differs at %d (state %d, i %d, k %d)", msg, z, tr1->st[z], tr1->i[z], tr1->k[z]);
        }

      /* null2 from the full posterior matrix, which it uses as workspace, so last */
      p7_Null2_ByExpectation(om, pp, n2a);
      for (x = 0; x < abc->Kp; x++)
        if (n2b[x] != n2a[x]) esl_fatal("%s: null2 differs for residue %d (%g %g)", msg, x, n2b[x], n2a[x]);
    }

  free(dsq); free(n2a); free(n2b);
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr1);
  p7_trace_Destroy(tr2);
  p7_oatmx_Destroy(ox);
  p7_omx_Destroy(fwd); p7_omx_Destroy(bck); p7_omx_Destroy(pp); p7_omx_Destroy(oa);
  p7_chkomx_Destroy(cx1); p7_chkomx_Destroy(cx2);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}


#endif /*p7OPTACC_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/




/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7OPTACC_TESTDRIVE

/* Failures in this test are to be expected, if you change the defaults.
 * The default RNG seed of 41 is very carefully chosen! Most seeds will
 * cause this test to fail. (Only 13 and 41 work for seeds 1..50.)
 *
 * The generic fwd/bck algorithms work in log space and suffer from a
 * small amount of imprecision, not only from the use of FLogsum()'s
 * table-driven approximation, but even (apparently) inherent in log()
 * and exp(). To minimize this, the generic decoding algorithm burns
 * time renormalizing each row. Still, it's a problem. See notes at
 * the header of utest_optacc() for more info.
 *
 * Another expected failure mode is when a fwd, bck nat score are close to
 * 0.0; FCompare() can evaluate two close-to-zero numbers as "different"
 * even if their absolute diff is negligible. (I suppose I could fix this.)
 */

/*
     ./optacc_utest
 */
#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"
#include "impl_neon.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "41", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,     "50", NULL, NULL,  NULL,  NULL, NULL, "size of random sequences to sample",             0 },
  { "-M",        eslARG_INT,     "45", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,     "20", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  { "--traces",  eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "dump all tracebacks",                            0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for NEON Forward, Backward implementations";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = NULL;
  P7_BG          *bg   = NULL;
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");

  /* first round of tests for DNA alphabets.  */
  if ((abc = esl_alphabet_Create(eslDNA)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))            == NULL)  esl_fatal("failed to create null model");

  utest_optacc(go, r, abc, bg, M, L, N);   /* normal sized models */
  utest_optacc(go, r, abc, bg, 1, L, 10);  /* size 1 models       */
  utest_optacc(go, r, abc, bg, M, 1, 10);  /* size 1 sequences    */

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  /* Second round of tests for amino alphabets.  */
  if ((abc = esl_alphabet_Create(eslAMINO)) == NULL)  esl_fatal("failed to create alphabet");
  if ((bg = p7_bg_Create(abc))              == NULL)  esl_fatal("failed to create null model");

  utest_optacc(go, r, abc, bg, M, L, N);
  utest_optacc(go, r, abc, bg, 1, L, 10);
  utest_optacc(go, r, abc, bg, M, 1, 10);
  /* after every utest_optacc() call, which passes only for the default seed's draws */
  utest_inplace(r, abc, bg, M, L, N);
  utest_inplace(r, abc, bg, 1, L, 10);
  utest_inplace(r, abc, bg, M, 1, 10);
  utest_sweep(r, abc, bg, M, 200, N);
  utest_sweep(r, abc, bg, 1, 200, 10);

  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}
#endif /*p7OPTACC_TESTDRIVE*/
/*------------------ end, test driver ---------------------------*/




/*****************************************************************
 * 6. Example
 *****************************************************************/
#ifdef p7OPTACC_EXAMPLE
/*
   gcc -g -Wall -o optacc_example -Dp7OPTACC_EXAMPLE -I.. -I../../easel -L.. -L../../easel optacc.c -lhmmer -leasel -lm
   ./optacc_example <hmmfile> <seqfile>
*/
#include <p7_config.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"
#include "impl_neon.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-d",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "dump posterior residue decoding matrix",           0 },
  { "-m",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "dump OA matrix",                                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile> <seqfile>";
static char banner[] = "example of optimal accuracy alignment, NEON implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  char           *seqfile = esl_opt_GetArg(go, 2);
  ESL_ALPHABET   *abc     = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bg      = NULL;
  P7_PROFILE     *gm      = NULL;
  P7_OPROFILE    *om      = NULL;
  P7_OMX         *ox1     = NULL;
  P7_OMX         *ox2     = NULL;
  P7_GMX         *gx      = NULL;
  ESL_SQ         *sq      = NULL;
  ESL_SQFILE     *sqfp    = NULL;
  P7_TRACE       *tr      = NULL;
  int             format  = eslSQFILE_UNKNOWN;
  char            errbuf[eslERRBUFSIZE];
  float           fsc, bsc;
  float           accscore;
  int             status;

  /* Read in one HMM */
  if (p7_hmmfile_Open(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abc, &hmm)           != eslOK) p7_Fail("Failed to read HMM");
  p7_hmmfile_Close(hfp);

  /* Read in one sequence */
  sq     = esl_sq_CreateDigital(abc);
  status = esl_sqfile_OpenDigital(abc, seqfile, format, NULL, &sqfp);
  if      (status == eslENOTFOUND) p7_Fail("No such file.");
  else if (status == eslEFORMAT)   p7_Fail("Format unrecognized.");
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect stdin or .gz.");
  else if (status != eslOK)        p7_Fail("Open failed, code %d.", status);
  if  (esl_sqio_Read(sqfp, sq) != eslOK) p7_Fail("Failed to read sequence");
  esl_sqfile_Close(sqfp);

  /* Configure a profile from the HMM */
  bg = p7_bg_Create(abc);                 p7_bg_SetLength(bg, sq->n);
  gm = p7_profile_Create(hmm->M, abc);    p7_ProfileConfig(hmm, bg, gm, sq->n, p7_LOCAL); /* multihit local: H3 default */
  om = p7_oprofile_Create(gm->M, abc);    p7_oprofile_Convert(gm, om);

  /* Allocations */
  ox1 = p7_omx_Create(gm->M, sq->n, sq->n);
  ox2 = p7_omx_Create(gm->M, sq->n, sq->n);
  gx  = p7_gmx_Create(gm->M, sq->n, sq->n, p7G_NSCELLS);
  tr  = p7_trace_CreateWithPP();
  p7_FLogsumInit();

  /* Run Forward, Backward; do OA fill and trace */
  p7_Forward (sq->dsq, sq->n, om, ox1,      &fsc);
  p7_Backward(sq->dsq, sq->n, om, ox1, ox2, &bsc);
  p7_Decoding(om, ox1, ox2, ox2);                   /* <gx2> is now the posterior decoding matrix */
  p7_OptimalAccuracy(om, ox2, ox1, &accscore);	    /* <gx1> is now the OA matrix */
  p7_OATrace(om, ox1, tr);

  if (esl_opt_GetBoolean(go, "-d")) { p7_omx_FDeconvert(ox2, gx);  p7_gmx_Dump(stdout, gx, p7_DEFAULT); }
  if (esl_opt_GetBoolean(go, "-m")) { p7_omx_FDeconvert(ox1, gx);  p7_gmx_Dump(stdout, gx, p7_DEFAULT); }

  p7_trace_Dump(stdout, tr, gm, sq->dsq);

  if (p7_trace_Validate(tr, abc, sq->dsq, errbuf) != eslOK) p7_Die("trace fails validation:\n%s\n", errbuf);

  printf("fwd = %.4f nats\n", fsc);
  printf("bck = %.4f nats\n", bsc);
  printf("acc = %.4f (%.2f%%)\n", accscore, accscore * 100. / (float) sq->n);

  /* Cleanup */
  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  p7_omx_Destroy(ox1);
  p7_omx_Destroy(ox2);
  p7_gmx_Destroy(gx);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7OPTACC_EXAMPLE*/
/*-------------------- end, example -----------------------------*/
