/* Checkpointed Forward and Backward algorithms, standard model,
 * generic (non-SIMD) versions.
 *
 * The Forward matrix is computed from row 0 up to row L in blocks of
 * rows held in a P7_CHKGMX created with <p7_CHK_FORWARD>. Only one
 * block is kept, together with checkpoints of the rows before each
 * block, so any block can be computed again later; a computation that
 * reads the Forward rows in decreasing order loads the blocks again
 * from last to first. The values are identical to <p7_GForward()>.
 *
 * The Backward matrix is the mirror image, in a P7_CHKGMX created with
 * <p7_CHK_BACKWARD>: rows L down to 0, checkpoints of the rows after
 * each block, blocks loaded again from first to last for a computation
 * that reads them in increasing order (a forward sweep). The values
 * are identical to <p7_GBackward()>. The special states are kept for
 * every row in both.
 *
 * Computing posterior probabilities from Backward rows as they are
 * finished, without a full Backward matrix, follows HMMER4's Backward
 * filter; checkpointing the Backward matrix for a forward sweep is an
 * elaboration of that idea.
 *
 * The frameshift checkpointed Backward is in
 * generic_fwdback_chk_frameshift.c.
 *
 * Contents:
 *   1. Checkpointed Forward and Backward
 *   2. Unit tests
 *   3. Test driver
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"

#include "hmmer.h"


/*****************************************************************
 * 1. Checkpointed Forward and Backward
 *****************************************************************/

/* Function:  p7_GForwardChk()
 * Synopsis:  Checkpointed Forward, standard model.
 *
 * Purpose:   Forward of profile <gm> against <dsq> of length <L> into
 *            the checkpointed matrix <cx>. Each row is computed by
 *            <p7_GForwardRow()>, as in <p7_GForward()>, so values and
 *            score are identical to it.
 *
 *            Rows are computed from 0 up to <L>, one block at a time;
 *            checkpoints of the rows before each block are saved, so
 *            <p7_GForwardChkLoad()> can compute any block again. On
 *            return the block holding row <L> is loaded and the special
 *            states of all rows are in <cx->xmx>.
 *
 *            The caller creates <cx> with <p7_CHK_FORWARD> (or sets that
 *            direction with <p7_chkgmx_SetDirection()>), <p7G_NSCELLS>
 *            cells per model position and at least 1 saved row, and grows
 *            it for <gm->M> by <L> with the block size it wants
 *            (<p7_chkgmx_GrowTo()>).
 *
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gm     - profile
 *            cx     - checkpointed DP matrix
 *            opt_sc - optRETURN: Forward lod score in nats
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <cx> isn't set up for this computation.
 */
int
p7_GForwardChk(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_CHKGMX *cx, float *opt_sc)
{
  float  *xmx = cx->xmx;
  int     b, i;
  int     status;

  if (cx->dir != p7_CHK_FORWARD || cx->saveR < 1 || cx->nscells != p7G_NSCELLS)
    ESL_EXCEPTION(eslEINVAL, "matrix not set up for the checkpointed Forward");

  p7_FLogsumInit();
  cx->M = gm->M;
  cx->L = L;

  for (b = 0; b <= L; b += cx->blockR)
    {
      if ((status = p7_chkgmx_StartBlock(cx, b)) != eslOK) return status;
      for (i = b; i <= ESL_MIN(L, b + cx->blockR - 1); i++)
        p7_GForwardRow(dsq, gm, cx->dp, cx->xmx, i);
    }

  if (opt_sc != NULL) *opt_sc = XMX(L,p7G_C) + gm->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}


/* Function:  p7_GForwardChkLoad()
 * Synopsis:  Compute one block of a checkpointed Forward matrix again.
 *
 * Purpose:   Load block <b> (a multiple of <cx->blockR>) of the Forward
 *            matrix that <p7_GForwardChk()> filled in <cx> for profile
 *            <gm> and sequence <dsq>, computing its rows again from the
 *            checkpoint before it. Afterwards rows <b> .. <b>+<blockR>-1
 *            (up to <cx->L>) hold the same values as rows of the full
 *            matrix from <p7_GForward()>. Blocks are loaded in decreasing
 *            order, as a computation that runs backward through the rows
 *            reads them (<p7_chkgmx_LoadBlock()>).
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if block <b> can't be loaded.
 */
int
p7_GForwardChkLoad(const ESL_DSQ *dsq, const P7_PROFILE *gm, P7_CHKGMX *cx, int b)
{
  int i;
  int status;

  if ((status = p7_chkgmx_LoadBlock(cx, b)) != eslOK) return status;
  for (i = b; i <= ESL_MIN(cx->L, b + cx->blockR - 1); i++)
    p7_GForwardRow(dsq, gm, cx->dp, cx->xmx, i);
  return eslOK;
}


/* Function:  p7_GBackwardChk()
 * Synopsis:  The Backward algorithm in a checkpointed matrix.
 *
 * Purpose:   The Backward algorithm of <p7_GBackward()>, computed in the
 *            checkpointed matrix <cx>. Rows are computed from <L> down
 *            to 0, one block at a time; checkpoints of the rows after
 *            each block are saved, so <p7_GBackwardChkLoad()> can compute
 *            any block again. On return block 0 is loaded and the special
 *            states of all rows are in <cx->xmx>.
 *
 *            The caller creates <cx> with <p7_CHK_BACKWARD>, <p7G_NSCELLS>
 *            cells per model position and at least 1 saved row, and grows
 *            it for <gm->M> by <L> with the block size it wants
 *            (<p7_chkgmx_GrowTo()>).
 *
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq (>= 1)
 *            gm     - profile
 *            cx     - checkpointed DP matrix
 *            opt_sc - optRETURN: Backward lod score in nats
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <cx> isn't set up for a backward computation
 *            with at least one saved row.
 */
int
p7_GBackwardChk(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_CHKGMX *cx, float *opt_sc)
{
  float **dp  = cx->dp;
  float  *xmx = cx->xmx;
  int     b, i;
  int     status;

  if (cx->dir != p7_CHK_BACKWARD || cx->saveR < 1 || cx->nscells != p7G_NSCELLS)
    ESL_EXCEPTION(eslEINVAL, "matrix not set up for the checkpointed Backward");

  p7_FLogsumInit();
  cx->M = gm->M;
  cx->L = L;

  for (b = (L / cx->blockR) * cx->blockR; b >= 0; b -= cx->blockR)
    {
      if ((status = p7_chkgmx_StartBlock(cx, b)) != eslOK) return status;
      for (i = ESL_MIN(L, b + cx->blockR - 1); i >= b; i--)
        p7_GBackwardRow(dsq, L, gm, dp, xmx, i);
    }

  if (opt_sc != NULL) *opt_sc = XMX(0,p7G_N);
  return eslOK;
}


/* Function:  p7_GBackwardChkLoad()
 * Synopsis:  Compute one block of a checkpointed Backward matrix again.
 *
 * Purpose:   Load block <b> (a multiple of <cx->blockR>) of the Backward
 *            matrix that <p7_GBackwardChk()> filled in <cx> for profile
 *            <gm> and sequence <dsq>, computing its rows again from the
 *            checkpoint after it. Afterwards rows <b> .. <b>+<blockR>-1
 *            (up to <cx->L>) hold the same values as rows of the full
 *            matrix from <p7_GBackward()>. Blocks can be loaded in any
 *            order.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <b> isn't the start of a block.
 */
int
p7_GBackwardChkLoad(const ESL_DSQ *dsq, const P7_PROFILE *gm, P7_CHKGMX *cx, int b)
{
  int i;
  int status;

  if ((status = p7_chkgmx_LoadBlock(cx, b)) != eslOK) return status;
  for (i = ESL_MIN(cx->L, b + cx->blockR - 1); i >= b; i--)
    p7_GBackwardRow(dsq, cx->L, gm, cx->dp, cx->xmx, i);
  return eslOK;
}


/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef p7GENERIC_FWDBACK_CHK_TESTDRIVE
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"

/* compare_block()
 * Rows of the loaded block <b> of <cx>, model positions 1..M, must
 * equal the same rows of the full matrix <gx> exactly. Column 0 and
 * the I cell at k = M (there is no I_M state) are left out: the full
 * Backward doesn't set them in every row, and nothing reads them.
 * Used for both the Backward and the Forward matrices.
 */
static void
compare_block(const char *msg, const P7_CHKGMX *cx, const P7_GMX *gx, int b)
{
  int i, k, s;

  for (i = b; i <= ESL_MIN(cx->L, b + cx->blockR - 1); i++)
    for (k = 1; k <= cx->M; k++)
      for (s = 0; s < p7G_NSCELLS; s++)
        if (!(k == cx->M && s == p7G_I) && cx->dp[i][k*p7G_NSCELLS + s] != gx->dp[i][k*p7G_NSCELLS + s])
          esl_fatal("%s: row %d, k %d, state %d differs (block size %d)", msg, i, k, s, cx->blockR);
}

/* compare_specials(): special states of all rows 0..L must be identical. */
static void
compare_specials(const char *msg, const P7_CHKGMX *cx, const P7_GMX *gx)
{
  int i, x;

  for (i = 0; i <= cx->L; i++)
    for (x = 0; x < p7G_NXCELLS; x++)
      if (cx->xmx[i*p7G_NXCELLS + x] != gx->xmx[i*p7G_NXCELLS + x])
        esl_fatal("%s: special state %d of row %d differs (block size %d)", msg, x, i, cx->blockR);
}

/* utest_fwdback_chk():
 * For random sequences of several lengths, compute the full Backward
 * matrix with p7_GBackward() and the checkpointed one with block sizes
 * from 1 row to more than L. The scores and special states must be
 * identical, and every block loaded again with p7_GBackwardChkLoad(),
 * first to last, must hold the same rows as the full matrix. The same
 * matrix, turned around by p7_chkgmx_SetDirection(), is also filled by
 * p7_GForwardChk(): score, special states and every block loaded again
 * with p7_GForwardChkLoad(), last to first, must equal p7_GForward().
 */
static void
utest_fwdback_chk(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, P7_PROFILE *gm, int nseq)
{
  char      *msg      = "checkpointed Backward unit test failed";
  int        lengths[] = { 1, 2, 3, 7, 50, 201 };
  int        blocks[]  = { 1, 2, 3, 8, 32, 300 };
  int        nL       = sizeof(lengths) / sizeof(int);
  int        nB       = sizeof(blocks)  / sizeof(int);
  P7_GMX    *gx       = p7_gmx_Create(gm->M, 201, 201, p7G_NSCELLS);
  P7_GMX    *gxf      = p7_gmx_Create(gm->M, 201, 201, p7G_NSCELLS);
  P7_CHKGMX *cx       = p7_chkgmx_Create(gm->M, 201, 8, 1, p7G_NSCELLS, p7_CHK_BACKWARD);
  ESL_DSQ   *dsq      = malloc(sizeof(ESL_DSQ) * (201 + 2));
  float      sc, fsc, chk_sc;
  int        n, j, q, b, L;

  if (gx == NULL || gxf == NULL || cx == NULL || dsq == NULL) esl_fatal("%s: allocation failed", msg);

  for (n = 0; n < nseq; n++)
    for (j = 0; j < nL; j++)
      {
        L = lengths[j];
        if (esl_rsq_xfIID(r, bg->f, abc->K, L, dsq)         != eslOK) esl_fatal("%s: seq generation failed", msg);
        if (p7_ReconfigLength(gm, L)                         != eslOK) esl_fatal(msg);
        if (p7_GBackward(dsq, L, gm, gx, &sc)                != eslOK) esl_fatal(msg);
        if (p7_GForward (dsq, L, gm, gxf, &fsc)              != eslOK) esl_fatal(msg);

        for (q = 0; q < nB; q++)
          {
            if (p7_chkgmx_GrowTo(cx, gm->M, L, blocks[q])    != eslOK) esl_fatal(msg);
            if (p7_GBackwardChk(dsq, L, gm, cx, &chk_sc)     != eslOK) esl_fatal(msg);
            if (chk_sc != sc) esl_fatal("%s: scores differ (%f %f), L %d, block size %d", msg, sc, chk_sc, L, blocks[q]);
            compare_specials(msg, cx, gx);
            compare_block(msg, cx, gx, 0);                    /* block 0 is loaded when the fill ends */

            for (b = 0; b <= L; b += cx->blockR)
              {
                if (p7_GBackwardChkLoad(dsq, gm, cx, b)      != eslOK) esl_fatal(msg);
                compare_block(msg, cx, gx, b);
              }

            /* the same matrix turned around for the checkpointed Forward:
             * the fill and every block loaded again, last to first, are
             * identical to p7_GForward() */
            if (p7_chkgmx_SetDirection(cx, p7_CHK_FORWARD, 1) != eslOK) esl_fatal(msg);
            if (p7_chkgmx_GrowTo(cx, gm->M, L, blocks[q])      != eslOK) esl_fatal(msg);
            if (p7_GForwardChk(dsq, L, gm, cx, &chk_sc)        != eslOK) esl_fatal(msg);
            if (chk_sc != fsc) esl_fatal("%s: Forward scores differ (%f %f), L %d, block size %d", msg, fsc, chk_sc, L, blocks[q]);
            compare_specials(msg, cx, gxf);
            for (b = (L / cx->blockR) * cx->blockR; b >= 0; b -= cx->blockR)
              {
                if (p7_GForwardChkLoad(dsq, gm, cx, b)        != eslOK) esl_fatal(msg);
                compare_block(msg, cx, gxf, b);
              }
            if (p7_chkgmx_SetDirection(cx, p7_CHK_BACKWARD, 1) != eslOK) esl_fatal(msg);
          }
      }

  free(dsq);
  p7_chkgmx_Destroy(cx);
  p7_gmx_Destroy(gx);
  p7_gmx_Destroy(gxf);
}
#endif /*p7GENERIC_FWDBACK_CHK_TESTDRIVE*/


/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef p7GENERIC_FWDBACK_CHK_TESTDRIVE
/*
   gcc -g -Wall -std=gnu99 -o generic_fwdback_chk_utest -I. -L. -I../easel -L../easel -Dp7GENERIC_FWDBACK_CHK_TESTDRIVE generic_fwdback_chk.c -lhmmer -leasel -lm
   ./generic_fwdback_chk_utest
 */
static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M",        eslARG_INT,     "45", NULL, NULL,  NULL,  NULL, NULL, "size of random models",                          0 },
  { "-N",        eslARG_INT,      "3", NULL, NULL,  NULL,  NULL, NULL, "number of sequences of each length",             0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for the generic checkpointed Forward and Backward";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  int             M      = esl_opt_GetInteger(go, "-M");
  int             nseq   = esl_opt_GetInteger(go, "-N");
  ESL_ALPHABET   *abc    = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg     = p7_bg_Create(abc);
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = NULL;

  p7_FLogsumInit();

  if (p7_hmm_Sample(r, M, abc, &hmm)             != eslOK) esl_fatal("failed to sample an HMM");
  if ((gm = p7_profile_Create(hmm->M, abc))      == NULL)  esl_fatal("failed to create profile");
  if (p7_ProfileConfig(hmm, bg, gm, 100, p7_LOCAL) != eslOK) esl_fatal("failed to config profile");

  utest_fwdback_chk(r, abc, bg, gm, nseq);

  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bg);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  printf("All tests passed.\n");
  return eslOK;
}
#endif /*p7GENERIC_FWDBACK_CHK_TESTDRIVE*/
