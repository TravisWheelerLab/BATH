/* Checkpointed frameshift Forward and Backward algorithms, 5 codon
 * lengths; NEON implementation.
 *
 * Contents:
 *    1. p7_ForwardChk_Frameshift(),  p7_ForwardChkLoad_Frameshift():  checkpointed Forward
 *       p7_BackwardChk_Frameshift(), p7_BackwardChkLoad_Frameshift(): checkpointed Backward
 *
 * Rows are computed by p7_ForwardRow3_Frameshift(),
 * p7_BackwardRowL_Frameshift(), p7_BackwardRow_Frameshift() and
 * p7_BackwardRow0_Frameshift() (fwdback_fs.c).
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <arm_neon.h>		/* NEON */

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"
#include "esl_neon.h"

#include "hmmer.h"
#include "impl_neon.h"


/*****************************************************************
 * 1. Checkpointed frameshift Forward and Backward
 *****************************************************************/

/* fwd_fs_ring_load()
 * Ring entries before row <b> (> 0) is computed again: rows b-4 .. b-1
 * (from 0), from the special states stored in <xmx>, each multiplied by
 * the inverse scale factors of the rows computed after it up to b-1, in
 * the same order as when the rows were first computed, so the entries
 * are identical.
 */
static void
fwd_fs_ring_load(P7_FS_FWDRING *rg, const float *xmx, int b)
{
  float sf;
  int   j, r, s;

  for (r = 0; r < PARSER_ROWS_FWD; r++)
    rg->xN[r] = rg->xB[r] = rg->xJ[r] = rg->xC[r] = 0.0f;
  for (j = ESL_MAX(0, b-PARSER_ROWS_FWD); j < b; j++)
    {
      s = j % PARSER_ROWS_FWD;
      rg->xN[s] = xmx[j*p7X_NXCELLS+p7X_N];
      rg->xB[s] = xmx[j*p7X_NXCELLS+p7X_B];
      rg->xJ[s] = xmx[j*p7X_NXCELLS+p7X_J];
      rg->xC[s] = xmx[j*p7X_NXCELLS+p7X_C];
      for (r = j+1; r < b; r++)
        if (xmx[r*p7X_NXCELLS+p7X_SCALE] > 1.0f)
          {
            sf = 1.0f / xmx[r*p7X_NXCELLS+p7X_SCALE];
            rg->xN[s] *= sf; rg->xB[s] *= sf; rg->xJ[s] *= sf; rg->xC[s] *= sf;
          }
    }
}


/* forward_ivx_load()
 * Codon partial sums before row <b> (> 0) is computed again: rows
 * b-4 .. b-1 (from 1) into <ivx>, each from the row before it (core
 * cells in the 3-cell rows <dp>, xB in <xmx>) as p7_ForwardRow3_Frameshift() computes
 * it, then multiplied by the inverse scale factors of its own row and
 * the rows after it up to b-1, in the same order as when the rows were
 * first computed, so the values are identical. The other slots are set
 * to 0.
 */
static void
forward_ivx_load(const P7_FS_OPROFILE *om_fs, float32x4_t * const *dp, const float *xmx, float32x4_t **ivx, int b)
{
  float32x4_t   zerov = vmovq_n_f32(0.0f);
  float32x4_t   mpv, dpv, ipv, xBv, sv, sfv;
  float32x4_t  *tp;
  float32x4_t  *dpp;
  int  Q = p7O_NQF(om_fs->M);
  int  q, r, j, slot;

  for (r = 0; r < p7P_5CODONS; r++)
    for (q = 0; q < Q; q++)
      ivx[r][q] = zerov;

  for (r = ESL_MAX(1, b-4); r < b; r++)
    {
      slot = r % p7P_5CODONS;
      dpp  = dp[r-1];
      xBv  = vmovq_n_f32(xmx[(r-1)*p7X_NXCELLS+p7X_B]);
      mpv  = vextq_f32(zerov, MMO(dpp,Q-1), 3);
      dpv  = vextq_f32(zerov, DMO(dpp,Q-1), 3);
      ipv  = vextq_f32(zerov, IMO(dpp,Q-1), 3);
      tp   = om_fs->tfv;
      for (q = 0; q < Q; q++)
        {
          sv  =            vmulq_f32(xBv, *tp);  tp++;
          sv  = vaddq_f32(sv,   vmulq_f32(mpv, *tp)); tp++;
          sv  = vaddq_f32(sv,   vmulq_f32(ipv, *tp)); tp++;
          sv  = vaddq_f32(sv,   vmulq_f32(dpv, *tp)); tp++;
          tp += 3;
          ivx[slot][q] = sv;
          mpv = MMO(dpp,q); dpv = DMO(dpp,q); ipv = IMO(dpp,q);
        }
      for (j = r; j < b; j++)
        if (xmx[j*p7X_NXCELLS+p7X_SCALE] > 1.0f)
          {
            sfv = vmovq_n_f32(1.0f / xmx[j*p7X_NXCELLS+p7X_SCALE]);
            for (q = 0; q < Q; q++) ivx[slot][q] = vmulq_f32(ivx[slot][q], sfv);
          }
    }
}


/* Function:  p7_ForwardChk_Frameshift()
 * Synopsis:  Frameshift Forward in a checkpointed matrix.
 *
 * Purpose:   The Forward algorithm of <p7_Forward_Frameshift()> for the
 *            5-codon profile <om_fs> against <dsq> of length <L>, computed
 *            in the checkpointed matrix <cx>, with <ov> as workspace. Each
 *            row is computed by <p7_ForwardRow3_Frameshift()> and
 *            holds the M (C0 sum), D and I cells, so these cells,
 *            special states, scale factors and score are identical to
 *            <p7_Forward_Frameshift()>.
 *
 *            Rows are computed from 0 up to <L>, one block at a time;
 *            checkpoints of the 5 rows before each block are saved, so
 *            <p7_ForwardChkLoad_Frameshift()> can compute any block
 *            again, with identical values. On return the block holding
 *            row <L> is loaded and the special states and scale factors
 *            of all rows are in <cx->xmx>.
 *
 *            The caller creates <cx> with <p7_CHK_FORWARD> (or sets that
 *            direction with <p7_chkomx_SetDirection()>), <p7X_NSCELLS>
 *            vectors per stripe and at least 5 saved rows, and grows it
 *            for <om_fs->M> by <L> with the block size it wants
 *            (<p7_chkomx_GrowTo()>).
 *
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq (>= 1)
 *            om_fs  - 5-codon frameshift profile
 *            cx     - checkpointed DP matrix
 *            ov     - workspace for the profile's M
 *            opt_sc - optRETURN: Forward lod score in nats
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if the score underflows to 0 (L > 1); <*opt_sc>
 *            is then -infinity.
 *
 * Throws:    <eslEINVAL> if <cx> isn't set up for this computation or
 *            <om_fs> isn't a 5-codon profile.
 *            <eslERANGE> if the score is NaN or overflows.
 */
int
p7_ForwardChk_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, float *opt_sc)
{
  P7_FS_FWDRING rg;
  float         totscale = 0.0f;
  float         xCL, xCLm1, xCLm2, xCtot;
  int           b, i;
  int           status;

  if (om_fs->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");
  if (cx->dir != p7_CHK_FORWARD || cx->saveR < 5 || cx->nscells != p7X_NSCELLS)
    ESL_EXCEPTION(eslEINVAL, "matrix not set up for the checkpointed frameshift Forward");
  if (L < 1) ESL_EXCEPTION(eslEINVAL, "sequence length must be >= 1");

  cx->M = om_fs->M;
  cx->L = L;

  for (b = 0; b <= L; b += cx->blockR)
    {
      if ((status = p7_chkomx_StartBlock(cx, b)) != eslOK) return status;
      for (i = b; i <= ESL_MIN(L, b + cx->blockR - 1); i++)
        p7_ForwardRow3_Frameshift(dsq, om_fs, cx->xmx, &totscale, &rg, ov->ivx, cx->dpf[i],
                                     (i >= 1) ? cx->dpf[i-1] : NULL, (i >= 3) ? cx->dpf[i-3] : NULL, i);
    }

  xCL   = rg.xC[   L    %PARSER_ROWS_FWD];
  xCLm1 = rg.xC[((L-1)%PARSER_ROWS_FWD+PARSER_ROWS_FWD)%PARSER_ROWS_FWD];
  xCLm2 = rg.xC[((L-2)%PARSER_ROWS_FWD+PARSER_ROWS_FWD)%PARSER_ROWS_FWD];
  xCtot = xCL + xCLm1*om_fs->xf[p7O_C][p7O_LOOP] + xCLm2*om_fs->xf[p7O_C][p7O_LOOP];

  if      (isnan(xCtot))         ESL_EXCEPTION(eslERANGE, "forward score is NaN");
  else if (isinf(xCtot) == 1)    ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");
  else if (L > 1 && xCtot == 0.0f) { if (opt_sc) *opt_sc = -eslINFINITY; return eslERANGE; }
  if (opt_sc != NULL) *opt_sc = totscale + logf(xCtot * om_fs->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}


/* Function:  p7_ForwardChkLoad_Frameshift()
 * Synopsis:  Compute one block of a checkpointed frameshift Forward again.
 *
 * Purpose:   After <p7_ForwardChk_Frameshift()> has filled <cx> for
 *            <om_fs> and <dsq>, compute block <b> (a multiple of
 *            <cx->blockR>) again from the checkpoint before it, with <ov>
 *            as workspace: the codon partial sums of the 4 rows before
 *            the block are computed from the saved rows and the special
 *            states of the rows before it are taken from <cx->xmx>
 *            first. Values are identical to those of the fill. Blocks are
 *            loaded in decreasing order, as a computation that runs
 *            backward through the rows reads them
 *            (<p7_chkomx_LoadBlock()>).
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if block <b> can't be loaded.
 */
int
p7_ForwardChkLoad_Frameshift(const ESL_DSQ *dsq, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, int b)
{
  P7_FS_FWDRING rg;
  float         totscale = 0.0f;   /* unused: the scale factors are already in <cx->xmx> */
  int           i;
  int           status;

  if ((status = p7_chkomx_LoadBlock(cx, b)) != eslOK) return status;
  if (b > 0)
    {
      fwd_fs_ring_load(&rg, cx->xmx, b);
      forward_ivx_load(om_fs, cx->dpf, cx->xmx, ov->ivx, b);
    }
  for (i = b; i <= ESL_MIN(cx->L, b + cx->blockR - 1); i++)
    p7_ForwardRow3_Frameshift(dsq, om_fs, cx->xmx, &totscale, &rg, ov->ivx, cx->dpf[i],
                                 (i >= 1) ? cx->dpf[i-1] : NULL, (i >= 3) ? cx->dpf[i-3] : NULL, i);
  return eslOK;
}


/* bwd_fs_ring_init()
 * Ring entries before row L is computed: 0, except C for the rows L+1
 * and L+2 past the end of the sequence, which rows L-1 and L-2 read.
 */
static void
bwd_fs_ring_init(P7_FS_BWDRING *rg, const P7_FS_OPROFILE *om_fs, int L)
{
  int r;
  for (r=0;r<PARSER_ROWS_BWD;r++)
    rg->xN[r]=rg->xB[r]=rg->xJ[r]=rg->xC[r]=0.0f;

  rg->xC[(L+1)%PARSER_ROWS_BWD]=om_fs->xf[p7O_C][p7O_MOVE];
  rg->xC[(L+2)%PARSER_ROWS_BWD]=om_fs->xf[p7O_C][p7O_MOVE];
}


/* bwd_fs_ring_load()
 * Ring entries before row <top> (< L) is computed again: rows <top>+1
 * .. <top>+3, the ones rows <top> .. <top>-2 read, from the special
 * states stored in <xmx> (the starting values of bwd_fs_ring_init()
 * for rows past <L>), each multiplied by the inverse scale factors of
 * the rows computed after it down to <top>+1, in the same order as when
 * the rows were first computed, so the entries are identical.
 */
static void
bwd_fs_ring_load(P7_FS_BWDRING *rg, const P7_FS_OPROFILE *om_fs, const float *xmx, int L, int top)
{
  float sf;
  int   j, r, s;

  bwd_fs_ring_init(rg, om_fs, L);
  for (j = top+1; j <= top+3; j++)
    {
      s = j % PARSER_ROWS_BWD;
      if (j <= L)
        {
          rg->xN[s] = xmx[j*p7X_NXCELLS+p7X_N];
          rg->xB[s] = xmx[j*p7X_NXCELLS+p7X_B];
          rg->xJ[s] = xmx[j*p7X_NXCELLS+p7X_J];
          rg->xC[s] = xmx[j*p7X_NXCELLS+p7X_C];
        }
      for (r = ESL_MIN(j-1, L); r > top; r--)
        if (xmx[r*p7X_NXCELLS+p7X_SCALE] > 1.0f)
          {
            sf = 1.0f / xmx[r*p7X_NXCELLS+p7X_SCALE];
            rg->xN[s] *= sf; rg->xB[s] *= sf; rg->xJ[s] *= sf; rg->xC[s] *= sf;
          }
    }
}


/* backward_frameshift_chk_block()
 * Compute rows <b> .. <b>+<blockR>-1 (up to <cx->L>) of the loaded
 * block <b> of the checkpointed frameshift Backward matrix <cx>, top
 * row first, with scale factors chosen by <scale_mode>. <rg> must hold
 * the ring entries for the block's top row (bwd_fs_ring_init() for the
 * block holding row L, else the state left by the block above or
 * bwd_fs_ring_load()). <zero_row> is a row of zeros for rows past L;
 * <ivxf> is workspace. For block 0 this includes row 0, whose core
 * cells are set to 0. If <totscale> is non-NULL the logs of the scale
 * factors are added to it. Returns the sum of the N values of rows 0,
 * 1 and 2 at row 0's scale (the Backward score before scaling) for
 * block 0, else 0.
 */
static float
backward_frameshift_chk_block(const ESL_DSQ *dsq, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_FS_BWDRING *rg,
                                  float32x4_t *zero_row, float32x4_t *ivxf, int b, float *totscale, int scale_mode)
{
  float32x4_t *dp[6];
  int     L    = cx->L;
  int     top  = ESL_MIN(L, b + cx->blockR - 1);
  int     Q    = p7O_NQF(om_fs->M);
  float   xN   = 0.0f;
  int     i, d, q;

  for (i = top; i >= ESL_MAX(b, 1); i--)
    {
      if (i == L)
        p7_BackwardRowL_Frameshift(om_fs, NULL, cx->xmx, totscale, rg, cx->dpf[L], L, scale_mode);
      else
        {
          for (d = 0; d < 6; d++) dp[d] = (i+d <= L) ? cx->dpf[i+d] : zero_row;
          p7_BackwardRow_Frameshift(dsq, L, om_fs, NULL, cx->xmx, totscale, NULL, rg, dp, ivxf, i, scale_mode);
        }
    }

  if (b == 0)
    {
      dp[0] = NULL;
      for (d = 1; d < 6; d++) dp[d] = (d <= L) ? cx->dpf[d] : zero_row;
      xN  = p7_BackwardRow0_Frameshift(dsq, L, om_fs, cx->xmx, rg, dp, ivxf);
      xN += rg->xN[1%PARSER_ROWS_BWD] + rg->xN[2%PARSER_ROWS_BWD];
      for (q = 0; q < Q * p7X_NSCELLS; q++) cx->dpf[0][q] = vmovq_n_f32(0.0f);
    }
  return xN;
}


/* Function:  p7_BackwardChk_Frameshift()
 * Synopsis:  Frameshift Backward in a checkpointed matrix, with its own scale factors.
 *
 * Purpose:   The Backward algorithm of <p7_Backward_Frameshift()> for
 *            the 5-codon profile <om_fs> and sequence <dsq> of length
 *            <L>, computed in the checkpointed matrix <cx>, with <ov>
 *            as workspace.
 *
 *            This version is meant to be run before Forward, and it uses
 *            its own scale factors: each row is scaled by xB when xB >
 *            1e4 (else not scaled), not by the Forward matrix's factor for
 *            the row, so no Forward matrix is needed. The factors are
 *            stored in <cx->xmx> with the special states. Values and
 *            score equal those of <p7_Backward_Frameshift()> up to
 *            rounding.
 *
 *            Rows are computed from <L> down to 0, one block at a time;
 *            checkpoints of the 5 rows after each block are saved, so
 *            <p7_BackwardChkLoad_Frameshift()> can compute any block
 *            again, with identical values. On return block 0 is loaded
 *            and the special states and scale factors of all rows are in
 *            <cx->xmx>. Row 0 holds only special states; its core cells
 *            are set to 0.
 *
 *            The caller creates <cx> with <p7_CHK_BACKWARD>, <p7X_NSCELLS>
 *            vectors per stripe and at least 5 saved rows, and grows it
 *            for <om_fs->M> by <L> with the block size it wants
 *            (<p7_chkomx_GrowTo()>).
 *
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq (>= 1)
 *            om_fs  - 5-codon frameshift profile
 *            cx     - checkpointed DP matrix
 *            ov     - workspace for the profile's M
 *            opt_sc - optRETURN: Backward lod score in nats
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if the score underflows to 0; <*opt_sc> is then
 *            -infinity.
 *
 * Throws:    <eslEINVAL> if <cx> isn't set up for this computation or
 *            <om_fs> isn't a 5-codon profile.
 *            <eslERANGE> if the score is NaN or overflows.
 *            <eslEMEM> on allocation failure.
 */
int
p7_BackwardChk_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, float *opt_sc)
{
  P7_FS_BWDRING rg;
  float32x4_t     *zero_mem = NULL;
  float32x4_t     *zero_row;
  int         Q        = p7O_NQF(om_fs->M);
  float       totscale = 0.0f;
  float       xNtot    = 0.0f;
  int         b, q;
  int         status;

  if (om_fs->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");
  if (cx->dir != p7_CHK_BACKWARD || cx->saveR < 5 || cx->nscells != p7X_NSCELLS)
    ESL_EXCEPTION(eslEINVAL, "matrix not set up for the checkpointed frameshift Backward");
  if (L < 1) ESL_EXCEPTION(eslEINVAL, "sequence length must be >= 1");

  ESL_ALLOC(zero_mem, sizeof(float32x4_t) * Q * p7X_NSCELLS + 15);
  zero_row = (float32x4_t *) (((unsigned long int) zero_mem + 15) & (~0xf));
  for (q = 0; q < Q * p7X_NSCELLS; q++) zero_row[q] = vmovq_n_f32(0.0f);

  cx->M = om_fs->M;
  cx->L = L;

  bwd_fs_ring_init(&rg, om_fs, L);
  for (b = (L / cx->blockR) * cx->blockR; b >= 0; b -= cx->blockR)
    {
      if ((status = p7_chkomx_StartBlock(cx, b)) != eslOK) goto ERROR;
      xNtot = backward_frameshift_chk_block(dsq, om_fs, cx, &rg, zero_row, ov->ivx[0], b, &totscale, p7_BWD_SCALE_OWN);
    }
  free(zero_mem);

  if      (isnan(xNtot))      ESL_EXCEPTION(eslERANGE, "backward score is NaN");
  else if (xNtot == 0.0f)     { if (opt_sc != NULL) *opt_sc = -eslINFINITY; return eslERANGE; }
  else if (isinf(xNtot) == 1) ESL_EXCEPTION(eslERANGE, "backward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = totscale + logf(xNtot);
  return eslOK;

 ERROR:
  free(zero_mem);
  return status;
}


/* Function:  p7_BackwardChkLoad_Frameshift()
 * Synopsis:  Compute one block of a checkpointed frameshift Backward matrix again.
 *
 * Purpose:   Load block <b> (a multiple of <cx->blockR>) of the Backward
 *            matrix that <p7_BackwardChk_Frameshift()> filled in <cx>
 *            for profile <om_fs> and sequence <dsq>, computing its rows
 *            again from the checkpoint after it with the scale factors
 *            stored in <cx->xmx>, with <ov> as workspace. Afterwards rows
 *            <b> .. <b>+<blockR>-1 (up to <cx->L>) hold the same values as
 *            when <p7_BackwardChk_Frameshift()> computed them. Blocks
 *            can be loaded in any order.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <b> isn't the start of a block.
 *            <eslEMEM> on allocation failure.
 */
int
p7_BackwardChkLoad_Frameshift(const ESL_DSQ *dsq, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, int b)
{
  P7_FS_BWDRING rg;
  float32x4_t     *zero_mem = NULL;
  float32x4_t     *zero_row;
  int         Q        = p7O_NQF(om_fs->M);
  int         top;
  int         q;
  int         status;

  if ((status = p7_chkomx_LoadBlock(cx, b)) != eslOK) return status;

  ESL_ALLOC(zero_mem, sizeof(float32x4_t) * Q * p7X_NSCELLS + 15);
  zero_row = (float32x4_t *) (((unsigned long int) zero_mem + 15) & (~0xf));
  for (q = 0; q < Q * p7X_NSCELLS; q++) zero_row[q] = vmovq_n_f32(0.0f);

  top = ESL_MIN(cx->L, b + cx->blockR - 1);
  if (top == cx->L) bwd_fs_ring_init(&rg, om_fs, cx->L);
  else              bwd_fs_ring_load(&rg, om_fs, cx->xmx, cx->L, top);

  backward_frameshift_chk_block(dsq, om_fs, cx, &rg, zero_row, ov->ivx[0], b, NULL, p7_BWD_SCALE_STORED);
  free(zero_mem);
  return eslOK;

 ERROR:
  free(zero_mem);
  return status;
}


/*-------------- end, forward/backward  -----------------*/

