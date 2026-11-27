/* Checkpointed Backward and Forward algorithms, standard model, NEON
 * implementation.
 *
 * Contents:
 *    1. p7_BackwardChk(), p7_BackwardChkLoad(): checkpointed Backward
 *       p7_ForwardChk(),  p7_ForwardChkLoad():  checkpointed Forward
 *
 * Rows are computed by p7_BackwardRow() and p7_ForwardRow()
 * (fwdback.c).
 */
#include <p7_config.h>

#include <math.h>

#include <arm_neon.h>		/* NEON */

#include "easel.h"
#include "esl_neon.h"

#include "hmmer.h"
#include "impl_neon.h"


/*****************************************************************
 * 1. Checkpointed Backward and Forward
 *****************************************************************/

/* backward_chk_block()
 * Compute rows <b> .. <b>+<blockR>-1 (up to <cx->L>) of the loaded
 * block <b> of the checkpointed Backward matrix <cx>, top row first,
 * with scale factors chosen by <scale_mode>. For block 0 this includes
 * row 0, whose core cells are set to 0. If <totscale> is non-NULL the
 * logs of the scale factors are added to it. Returns row 0's N value
 * for block 0, else 0.
 */
static float
backward_chk_block(const ESL_DSQ *dsq, const P7_OPROFILE *om, P7_CHKOMX *cx, int b, float *totscale, int scale_mode)
{
  int   L   = cx->L;
  int   top = ESL_MIN(L, b + cx->blockR - 1);
  float xN  = 0.0f;
  int   i;

  for (i = top; i >= b; i--)
    xN = p7_BackwardRow(dsq, L, om, NULL, cx->xmx, totscale, NULL, cx->dpf[i], (i < L) ? cx->dpf[i+1] : NULL, i, scale_mode);
  return xN;
}


/* Function:  p7_BackwardChk()
 * Synopsis:  Backward in a checkpointed matrix, with its own scale factors.
 *
 * Purpose:   The Backward algorithm of <p7_Backward()> for profile
 *            <om> and sequence <dsq> of length <L>, computed in the
 *            checkpointed matrix <cx>.
 *
 *            This version is meant to be run before Forward, and it uses
 *            its own scale factors: each row is scaled by xB when xB >
 *            1e4 (else not scaled), not by the Forward matrix's factor for
 *            the row, so no Forward matrix is needed. The factors are
 *            stored in <cx->xmx> with the special states. Values and
 *            score equal those of <p7_Backward()> up to rounding.
 *
 *            Rows are computed from <L> down to 0, one block at a time;
 *            checkpoints of the row after each block are saved, so
 *            <p7_BackwardChkLoad()> can compute any block again, with
 *            identical values. On return block 0 is loaded and the special
 *            states and scale factors of all rows are in <cx->xmx>. Row 0
 *            holds only special states; its core cells are set to 0.
 *
 *            The caller creates <cx> with <p7_CHK_BACKWARD>, <p7X_NSCELLS>
 *            vectors per stripe and at least 1 saved row, and grows it for
 *            <om->M> by <L> with the block size it wants
 *            (<p7_chkomx_GrowTo()>).
 *
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq (>= 1)
 *            om     - profile
 *            cx     - checkpointed DP matrix
 *            opt_sc - optRETURN: Backward lod score in nats
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <cx> isn't set up for this computation.
 *            <eslERANGE> if the score is NaN, underflows to 0 or overflows.
 */
int
p7_BackwardChk(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc)
{
  float totscale = 0.0f;
  float xN       = 0.0f;
  int   b;
  int   status;

  if (cx->dir != p7_CHK_BACKWARD || cx->saveR < 1 || cx->nscells != p7X_NSCELLS)
    ESL_EXCEPTION(eslEINVAL, "matrix not set up for the checkpointed Backward");
  if (L < 1) ESL_EXCEPTION(eslEINVAL, "sequence length must be >= 1");

  cx->M = om->M;
  cx->L = L;

  for (b = (L / cx->blockR) * cx->blockR; b >= 0; b -= cx->blockR)
    {
      if ((status = p7_chkomx_StartBlock(cx, b)) != eslOK) return status;
      xN = backward_chk_block(dsq, om, cx, b, &totscale, p7_BWD_SCALE_OWN);
    }

  if      (isnan(xN))      ESL_EXCEPTION(eslERANGE, "backward score is NaN");
  else if (xN == 0.0f)     ESL_EXCEPTION(eslERANGE, "backward score underflow (is 0.0)");
  else if (isinf(xN) == 1) ESL_EXCEPTION(eslERANGE, "backward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = totscale + log(xN);
  return eslOK;
}


/* Function:  p7_BackwardChkLoad()
 * Synopsis:  Compute one block of a checkpointed Backward matrix again.
 *
 * Purpose:   Load block <b> (a multiple of <cx->blockR>) of the Backward
 *            matrix that <p7_BackwardChk()> filled in <cx> for profile
 *            <om> and sequence <dsq>, computing its rows again from the
 *            checkpoint after it with the scale factors stored in
 *            <cx->xmx>. Afterwards rows <b> .. <b>+<blockR>-1 (up to
 *            <cx->L>) hold the same values as when <p7_BackwardChk()>
 *            computed them. Blocks can be loaded in any order.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <b> isn't the start of a block.
 */
int
p7_BackwardChkLoad(const ESL_DSQ *dsq, const P7_OPROFILE *om, P7_CHKOMX *cx, int b)
{
  int status;

  if ((status = p7_chkomx_LoadBlock(cx, b)) != eslOK) return status;
  backward_chk_block(dsq, om, cx, b, NULL, p7_BWD_SCALE_STORED);
  return eslOK;
}

/* Function:  p7_ForwardChk()
 * Synopsis:  Checkpointed Forward, NEON.
 *
 * Purpose:   Forward of profile <om> against <dsq> of length <L> into
 *            the checkpointed matrix <cx>. Each row is computed as in
 *            <p7_Forward()>, by <p7_ForwardRow()>, so cells,
 *            special states, scale factors and score are identical to
 *            it.
 *
 *            Rows are computed from 0 up to <L>, one block at a time;
 *            checkpoints of the rows before each block are saved, so
 *            <p7_ForwardChkLoad()> can compute any block again, with
 *            identical values. On return the block holding row <L> is
 *            loaded and the special states and scale factors of all rows
 *            are in <cx->xmx>.
 *
 *            The caller creates <cx> with <p7_CHK_FORWARD> (or sets that
 *            direction with <p7_chkomx_SetDirection()>), <p7X_NSCELLS>
 *            vectors per stripe and at least 1 saved row, and grows it
 *            for <om->M> by <L> with the block size it wants
 *            (<p7_chkomx_GrowTo()>).
 *
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq (>= 1)
 *            om     - profile
 *            cx     - checkpointed DP matrix
 *            opt_sc - optRETURN: Forward lod score in nats
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <cx> isn't set up for this computation.
 *            <eslERANGE> if the score is NaN, underflows to 0 or overflows.
 */
int
p7_ForwardChk(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc)
{
  float totscale = 0.0f;
  float xC;
  int   b, i;
  int   status;

  if (cx->dir != p7_CHK_FORWARD || cx->saveR < 1 || cx->nscells != p7X_NSCELLS)
    ESL_EXCEPTION(eslEINVAL, "matrix not set up for the checkpointed Forward");
  if (L < 1) ESL_EXCEPTION(eslEINVAL, "sequence length must be >= 1");

  cx->M = om->M;
  cx->L = L;

  for (b = 0; b <= L; b += cx->blockR)
    {
      if ((status = p7_chkomx_StartBlock(cx, b)) != eslOK) return status;
      for (i = b; i <= ESL_MIN(L, b + cx->blockR - 1); i++)
        p7_ForwardRow(dsq, om, cx->xmx, &totscale, cx->dpf[i], (i > 0 ? cx->dpf[i-1] : NULL), i);
    }
  xC = cx->xmx[L*p7X_NXCELLS+p7X_C];

  if      (isnan(xC))      ESL_EXCEPTION(eslERANGE, "forward score is NaN");
  else if (xC == 0.0f)     ESL_EXCEPTION(eslERANGE, "forward score underflow (is 0.0)");
  else if (isinf(xC) == 1) ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = totscale + log(xC * om->xf[p7O_C][p7O_MOVE]);
  return eslOK;
}


/* Function:  p7_ForwardChkLoad()
 * Synopsis:  Compute one block of a checkpointed Forward again.
 *
 * Purpose:   After <p7_ForwardChk()> has filled <cx>, compute block
 *            <b> (a multiple of <cx->blockR>) again from the checkpoint
 *            before it, loading the block and its saved rows into <cx>.
 *            Values are identical to those of the fill. Blocks are loaded
 *            in decreasing order, as a computation that runs backward
 *            through the rows reads them (<p7_chkomx_LoadBlock()>).
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if block <b> can't be loaded.
 */
int
p7_ForwardChkLoad(const ESL_DSQ *dsq, const P7_OPROFILE *om, P7_CHKOMX *cx, int b)
{
  float totscale = 0.0f;   /* unused: the scale factors are already in <cx->xmx> */
  int   i;
  int   status;

  if ((status = p7_chkomx_LoadBlock(cx, b)) != eslOK) return status;
  for (i = b; i <= ESL_MIN(cx->L, b + cx->blockR - 1); i++)
    p7_ForwardRow(dsq, om, cx->xmx, &totscale, cx->dpf[i], (i > 0 ? cx->dpf[i-1] : NULL), i);
  return eslOK;
}
