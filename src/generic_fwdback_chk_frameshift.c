/* Checkpointed Forward and Backward algorithms, frameshift model (5
 * codon lengths), generic (non-SIMD) versions.
 *
 * The frameshift Forward matrix is computed from row 0 up to row L in
 * blocks of rows held in a P7_CHKGMX created with <p7_CHK_FORWARD>, as
 * the standard checkpointed Forward in generic_fwdback_chk.c is. Each
 * row reads rows i-1 and i-3 and the codon partial sums of rows
 * i-1..i-4, which are computed from rows i-2..i-5, so checkpoints save
 * 5 rows and the partial sums are computed again from them when a
 * block is loaded. The values are identical to
 * <p7_GForward_Frameshift()>.
 *
 * The frameshift Backward matrix is the mirror image, computed from row
 * L down to row 0 in a P7_CHKGMX created with <p7_CHK_BACKWARD>. Each
 * row reads the 5 rows after it, so checkpoints save 5 rows. The values
 * are identical to <p7_GBackward_Frameshift()>.
 *
 * Computing posterior probabilities from Backward rows as they are
 * finished, without a full Backward matrix, follows HMMER4's Backward
 * filter; checkpointing the Backward matrix for a forward sweep is an
 * elaboration of that idea.
 *
 * Contents:
 *   1. Checkpointed Forward, frameshift model
 *   2. Checkpointed Backward, frameshift model
 *   3. Unit tests
 *   4. Test driver
 */
#include "p7_config.h"

#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"

#include "hmmer.h"

#define IVX5(i,k)   (ivx[((k)*p7P_5CODONS) + (i)])


/*****************************************************************
 * 1. Checkpointed Forward, frameshift model
 *****************************************************************/

/* forward_ivx_row()
 * The codon partial sums of Forward row <r> (>= 1) into <ivx>, from row
 * r-1 (core cells and special states) of the Forward matrix <dp>,
 * <xmx>. Same calculation as p7_GForwardRow_Frameshift() does for row
 * <r>.
 */
static void
forward_ivx_row(const P7_FS_PROFILE *gm_fs5, float **dp, const float *xmx, float *ivx, int r)
{
  float const *tsc  = gm_fs5->tsc;
  int          M    = gm_fs5->M;
  int          slot = r % p7P_5CODONS;
  int          k;

  if (r <= 2)
    for (k = 1; k <= M; k++)
      IVX5(slot,k) = XMX_FS(r-1,p7G_B) + TSC(p7P_BM,k-1);
  else
    for (k = 1; k <= M; k++)
      IVX5(slot,k) = p7_FLogsum(MMX_FS(r-1,k-1,p7G_C0) + TSC(p7P_MM,k-1),
                     p7_FLogsum(IMX_FS(r-1,k-1)        + TSC(p7P_IM,k-1),
                     p7_FLogsum(DMX_FS(r-1,k-1)        + TSC(p7P_DM,k-1),
                                XMX_FS(r-1,p7G_B)      + TSC(p7P_BM,k-1))));
}


/* Function:  p7_GForwardChk_Frameshift()
 * Synopsis:  The frameshift Forward algorithm in a checkpointed matrix.
 *
 * Purpose:   The Forward algorithm of <p7_GForward_Frameshift()> (5
 *            codon lengths) for <gm_fs5> against <dsq> of length <L>,
 *            computed in the checkpointed matrix <cx>, with <iv> as
 *            workspace. Each row is computed by
 *            <p7_GForwardRow_Frameshift()>, so values and score are
 *            identical to <p7_GForward_Frameshift()>.
 *
 *            Rows are computed from 0 up to <L>, one block at a time;
 *            checkpoints of the 5 rows before each block are saved, so
 *            <p7_GForwardChkLoad_Frameshift()> can compute any block
 *            again. On return the block holding row <L> is loaded and
 *            the special states of all rows are in <cx->xmx>.
 *
 *            The caller creates <cx> with <p7_CHK_FORWARD> (or sets that
 *            direction with <p7_chkgmx_SetDirection()>),
 *            <p7G_NSCELLS_FS> cells per model position and at least 5
 *            saved rows, and grows it for <gm_fs5->M> by <L> with the
 *            block size it wants (<p7_chkgmx_GrowTo()>).
 *
 *            Caller must have initialized the log-sum calculation with
 *            a call to <p7_FLogsumInit()>.
 *
 * Args:      dsq    - nucleotide sequence in digitized form, 1..L
 *            L      - length of dsq (>= 1)
 *            gm_fs5 - frameshift aware codon profile, 5 codon lengths
 *            cx     - checkpointed DP matrix
 *            iv     - intermediate value matrix, workspace
 *            opt_sc - optRETURN: Forward lod score in nats
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <gm_fs5> isn't for 5 codon lengths, or if
 *            <cx> isn't set up for a forward computation with at least
 *            5 saved rows.
 */
int
p7_GForwardChk_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs5, P7_CHKGMX *cx, P7_IVX *iv, float *opt_sc)
{
  float  *xmx = cx->xmx;
  int     b, i;
  int     status;

  if (gm_fs5->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");
  if (cx->dir != p7_CHK_FORWARD || cx->saveR < 5 || cx->nscells != p7G_NSCELLS_FS)
    ESL_EXCEPTION(eslEINVAL, "matrix not set up for the checkpointed frameshift Forward");

  cx->M = gm_fs5->M;
  cx->L = L;

  for (b = 0; b <= L; b += cx->blockR)
    {
      if ((status = p7_chkgmx_StartBlock(cx, b)) != eslOK) return status;
      for (i = b; i <= ESL_MIN(L, b + cx->blockR - 1); i++)
        p7_GForwardRow_Frameshift(dsq, gm_fs5, cx->dp, cx->xmx, iv->ivx, i);
    }

  if (opt_sc != NULL) *opt_sc = p7_FLogsum( XMX_FS(L,p7G_C),
                                p7_FLogsum( XMX_FS(L-1,p7G_C) + gm_fs5->xsc[p7P_C][p7P_LOOP],
                                            (L >= 2 ? XMX_FS(L-2,p7G_C) + gm_fs5->xsc[p7P_C][p7P_LOOP] : -eslINFINITY))) +
                                            gm_fs5->xsc[p7P_C][p7P_MOVE];
  return eslOK;
}


/* Function:  p7_GForwardChkLoad_Frameshift()
 * Synopsis:  Compute one block of a checkpointed frameshift Forward matrix again.
 *
 * Purpose:   Load block <b> (a multiple of <cx->blockR>) of the frameshift
 *            Forward matrix that <p7_GForwardChk_Frameshift()> filled in
 *            <cx> for <gm_fs5> and <dsq>, computing its rows again from
 *            the checkpoint before it, with <iv> as workspace: the codon
 *            partial sums of the 4 rows before the block are computed
 *            from the saved rows first. Afterwards rows <b> ..
 *            <b>+<blockR>-1 (up to <cx->L>) hold the same values as rows
 *            of the full matrix from <p7_GForward_Frameshift()>. Blocks
 *            are loaded in decreasing order, as a computation that runs
 *            backward through the rows reads them
 *            (<p7_chkgmx_LoadBlock()>).
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if block <b> can't be loaded.
 */
int
p7_GForwardChkLoad_Frameshift(const ESL_DSQ *dsq, const P7_FS_PROFILE *gm_fs5, P7_CHKGMX *cx, P7_IVX *iv, int b)
{
  float *ivx = iv->ivx;
  int    i, k, r;
  int    status;

  if ((status = p7_chkgmx_LoadBlock(cx, b)) != eslOK) return status;

  if (b > 0)
    {
      for (k = 0; k <= gm_fs5->M; k++)
        for (r = 0; r < p7P_5CODONS; r++)
          IVX5(r,k) = -eslINFINITY;
      for (r = ESL_MAX(1, b-4); r < b; r++)
        forward_ivx_row(gm_fs5, cx->dp, cx->xmx, ivx, r);
    }

  for (i = b; i <= ESL_MIN(cx->L, b + cx->blockR - 1); i++)
    p7_GForwardRow_Frameshift(dsq, gm_fs5, cx->dp, cx->xmx, ivx, i);
  return eslOK;
}


/*****************************************************************
 * 2. Checkpointed Backward, frameshift model
 *****************************************************************/

/* Function:  p7_GBackwardChk_Frameshift()
 * Synopsis:  The frameshift Backward algorithm in a checkpointed matrix.
 *
 * Purpose:   The Backward algorithm of <p7_GBackward_Frameshift()> (5
 *            codon lengths), computed in the checkpointed matrix <cx>.
 *            Rows are computed from <L> down to 0, one block at a time;
 *            checkpoints of the 5 rows after each block are saved, so
 *            <p7_GBackwardChkLoad_Frameshift()> can compute any block
 *            again. On return block 0 is loaded and the special states
 *            of all rows are in <cx->xmx>.
 *
 *            The caller creates <cx> with <p7_CHK_BACKWARD>, <p7G_NSCELLS>
 *            cells per model position and at least 5 saved rows, and
 *            grows it for <gm_fs5->M> by <L> with the block size it wants
 *            (<p7_chkgmx_GrowTo()>).
 *
 * Args:      dsq    - nucleotide sequence in digitized form, 1..L
 *            L      - length of dsq (>= 1)
 *            gm_fs5 - frameshift aware codon profile, 5 codon lengths
 *            cx     - checkpointed DP matrix
 *            iv     - intermediate value matrix, workspace
 *            opt_sc - optRETURN: Backward lod score in nats
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <gm_fs5> isn't for 5 codon lengths, or if
 *            <cx> isn't set up for a backward computation with at least
 *            5 saved rows.
 */
int
p7_GBackwardChk_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs5, P7_CHKGMX *cx, P7_IVX *iv, float *opt_sc)
{
  float **dp  = cx->dp;
  float  *xmx = cx->xmx;
  int     b, i;
  int     status;

  if (gm_fs5->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");
  if (cx->dir != p7_CHK_BACKWARD || cx->saveR < 5 || cx->nscells != p7G_NSCELLS)
    ESL_EXCEPTION(eslEINVAL, "matrix not set up for the checkpointed frameshift Backward");

  cx->M = gm_fs5->M;
  cx->L = L;

  for (b = (L / cx->blockR) * cx->blockR; b >= 0; b -= cx->blockR)
    {
      if ((status = p7_chkgmx_StartBlock(cx, b)) != eslOK) return status;
      for (i = ESL_MIN(L, b + cx->blockR - 1); i >= b; i--)
        p7_GBackwardRow_Frameshift(dsq, L, gm_fs5, dp, xmx, iv->ivx, i);
    }

  if (opt_sc != NULL) *opt_sc =  p7_FLogsum( XMX(0,p7G_N),
                                 p7_FLogsum( XMX(1,p7G_N),
                                             (L >= 2 ? XMX(2,p7G_N) : -eslINFINITY)));
  return eslOK;
}


/* Function:  p7_GBackwardChkLoad_Frameshift()
 * Synopsis:  Compute one block of a checkpointed frameshift Backward matrix again.
 *
 * Purpose:   Load block <b> (a multiple of <cx->blockR>) of the frameshift
 *            Backward matrix that <p7_GBackwardChk_Frameshift()> filled
 *            in <cx> for <gm_fs5> and <dsq>, computing its rows again
 *            from the checkpoint after it, with <iv> as workspace.
 *            Afterwards rows <b> .. <b>+<blockR>-1 (up to <cx->L>) hold
 *            the same values as rows of the full matrix from
 *            <p7_GBackward_Frameshift()>. Blocks can be loaded in any
 *            order.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <b> isn't the start of a block.
 */
int
p7_GBackwardChkLoad_Frameshift(const ESL_DSQ *dsq, const P7_FS_PROFILE *gm_fs5, P7_CHKGMX *cx, P7_IVX *iv, int b)
{
  int i;
  int status;

  if ((status = p7_chkgmx_LoadBlock(cx, b)) != eslOK) return status;
  for (i = ESL_MIN(cx->L, b + cx->blockR - 1); i >= b; i--)
    p7_GBackwardRow_Frameshift(dsq, cx->L, gm_fs5, cx->dp, cx->xmx, iv->ivx, i);
  return eslOK;
}


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7GENERIC_FWDBACK_CHK_FRAMESHIFT_TESTDRIVE
#include "esl_getopts.h"
#include "esl_gencode.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"

/* compare_block()
 * Every cell, k = 0..M, of the rows of the loaded block <b> of <cx>
 * must equal the same rows of the full matrix <gx> exactly. <what>
 * names the matrix in the error message.
 */
static void
compare_block(const char *msg, const char *what, const P7_CHKGMX *cx, const P7_GMX *gx, int b)
{
  int i, j;

  for (i = b; i <= ESL_MIN(cx->L, b + cx->blockR - 1); i++)
    for (j = 0; j < (cx->M+1) * cx->nscells; j++)
      if (cx->dp[i][j] != gx->dp[i][j])
        esl_fatal("%s: %s row %d, k %d, cell %d differs (block size %d)", msg, what, i, j / cx->nscells, j % cx->nscells, cx->blockR);
}

/* compare_specials(): special states of all rows 0..L must be identical. */
static void
compare_specials(const char *msg, const char *what, const P7_CHKGMX *cx, const P7_GMX *gx)
{
  int i, x;

  for (i = 0; i <= cx->L; i++)
    for (x = 0; x < p7G_NXCELLS; x++)
      if (cx->xmx[i*p7G_NXCELLS + x] != gx->xmx[i*p7G_NXCELLS + x])
        esl_fatal("%s: %s special state %d of row %d differs (block size %d)", msg, what, x, i, cx->blockR);
}

/* utest_fwdback_chk_frameshift():
 * For random DNA sequences of several lengths, p7_GForward_Frameshift()
 * against p7_GForwardChk_Frameshift() and p7_GBackward_Frameshift()
 * against p7_GBackwardChk_Frameshift(), with 5 saved rows, including
 * lengths below 6 where rows near 0 and L have their own
 * initialization, and block sizes from 5 rows to more than L. Scores
 * and special states must be identical, and every block loaded again,
 * Forward last to first and Backward first to last, must hold the same
 * rows as the full matrix. Sequences contain a few N's so the codon
 * lookups for non-ACGT nucleotides are exercised.
 */
static void
utest_fwdback_chk_frameshift(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA, P7_HMM *hmm, int nseq)
{
  char          *msg      = "checkpointed frameshift Forward/Backward unit test failed";
  int            lengths[] = { 1, 2, 3, 4, 5, 6, 7, 11, 60, 301 };
  int            blocks[]  = { 5, 6, 7, 13, 64, 400 };
  int            nL       = sizeof(lengths) / sizeof(int);
  int            nB       = sizeof(blocks)  / sizeof(int);
  P7_FS_PROFILE *gm_fs5   = p7_profile_fs_Create(hmm->M, abcAA, p7P_5CODONS);
  P7_GMX        *gx       = p7_gmx_Create(hmm->M, 301, 301, p7G_NSCELLS);
  P7_GMX        *fx       = p7_gmx_Create(hmm->M, 301, 301, p7G_NSCELLS_FS);
  P7_CHKGMX     *cx       = p7_chkgmx_Create(hmm->M, 301, 8, 5, p7G_NSCELLS, p7_CHK_BACKWARD);
  P7_CHKGMX     *fcx      = p7_chkgmx_Create(hmm->M, 301, 8, 5, p7G_NSCELLS_FS, p7_CHK_FORWARD);
  P7_IVX        *iv       = p7_ivx_Create(hmm->M, p7P_5CODONS);
  ESL_DSQ       *dsq      = malloc(sizeof(ESL_DSQ) * (301 + 2));
  float          sc, fsc, chk_sc;
  int            n, j, q, b, i, L;

  if (gm_fs5 == NULL || gx == NULL || fx == NULL || cx == NULL || fcx == NULL || iv == NULL || dsq == NULL) esl_fatal("%s: allocation failed", msg);

  for (n = 0; n < nseq; n++)
    for (j = 0; j < nL; j++)
      {
        L = lengths[j];
        if (esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq)   != eslOK) esl_fatal("%s: seq generation failed", msg);
        for (i = 1; i <= L; i++)
          if (esl_rnd_Roll(r, 40) == 0) dsq[i] = esl_abc_DigitizeSymbol(abcDNA, 'N');
        if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, ESL_MAX(1, L/3), p7_UNILOCAL) != eslOK) esl_fatal(msg);
        if (p7_GForward_Frameshift (dsq, L, gm_fs5, fx, iv, &fsc)  != eslOK) esl_fatal(msg);
        if (p7_GBackward_Frameshift(dsq, L, gm_fs5, gx, iv, &sc)   != eslOK) esl_fatal(msg);

        for (q = 0; q < nB; q++)
          {
            if (p7_chkgmx_GrowTo(cx, hmm->M, L, blocks[q])                 != eslOK) esl_fatal(msg);
            if (p7_GBackwardChk_Frameshift(dsq, L, gm_fs5, cx, iv, &chk_sc) != eslOK) esl_fatal(msg);
            if (chk_sc != sc) esl_fatal("%s: scores differ (%f %f), L %d, block size %d", msg, sc, chk_sc, L, blocks[q]);
            compare_specials(msg, "Backward", cx, gx);
            compare_block(msg, "Backward", cx, gx, 0);

            for (b = 0; b <= L; b += cx->blockR)
              {
                if (p7_GBackwardChkLoad_Frameshift(dsq, gm_fs5, cx, iv, b) != eslOK) esl_fatal(msg);
                compare_block(msg, "Backward", cx, gx, b);
              }

            if (p7_chkgmx_GrowTo(fcx, hmm->M, L, blocks[q])                != eslOK) esl_fatal(msg);
            if (p7_GForwardChk_Frameshift(dsq, L, gm_fs5, fcx, iv, &chk_sc) != eslOK) esl_fatal(msg);
            if (chk_sc != fsc) esl_fatal("%s: Forward scores differ (%f %f), L %d, block size %d", msg, fsc, chk_sc, L, blocks[q]);
            compare_specials(msg, "Forward", fcx, fx);
            compare_block(msg, "Forward", fcx, fx, (L / fcx->blockR) * fcx->blockR);

            for (b = (L / fcx->blockR) * fcx->blockR; b >= 0; b -= fcx->blockR)
              {
                if (p7_GForwardChkLoad_Frameshift(dsq, gm_fs5, fcx, iv, b) != eslOK) esl_fatal(msg);
                compare_block(msg, "Forward", fcx, fx, b);
              }
          }
      }

  free(dsq);
  p7_ivx_Destroy(iv);
  p7_chkgmx_Destroy(cx);
  p7_chkgmx_Destroy(fcx);
  p7_gmx_Destroy(gx);
  p7_gmx_Destroy(fx);
  p7_profile_fs_Destroy(gm_fs5);
}

/* utest_stotrace_chk_frameshift():
 * Stochastic traces from a checkpointed frameshift Forward matrix against
 * p7_GStochasticTrace_Frameshift() on the full Forward matrix. For random
 * DNA and the back-translated consensus of a model with sharpened match
 * emissions between random flanks, <ntr> traces sampled together by
 * p7_GStochasticTraceChk_Frameshift() must be identical (states, k, i and
 * codon lengths) to those sampled one at a time with the same per-trace
 * random number generators, with the checkpointed matrix in one block and
 * in blocks of 5..12 rows.
 */
static void
utest_stotrace_chk_frameshift(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA, int M, int L, int N)
{
  char            *msg    = "checkpointed frameshift stochastic trace unit test failed";
  int              ntr    = 20;
  P7_HMM          *hmm    = NULL;
  P7_CODONTABLE   *ct     = p7_codontable_Create(gcode);
  P7_FS_PROFILE   *gm_fs5 = p7_profile_fs_Create(M, abcAA, p7P_5CODONS);
  P7_GMX          *fx     = p7_gmx_Create(M, L, L, p7G_NSCELLS_FS);
  P7_CHKGMX       *cx     = p7_chkgmx_Create(M, L, L+1, 5, p7G_NSCELLS_FS, p7_CHK_FORWARD);
  P7_IVX          *iv     = p7_ivx_Create(M, p7P_5CODONS);
  ESL_SQ          *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ         *dsq    = malloc(sizeof(ESL_DSQ) * (L+2));
  ESL_RANDOMNESS  *rng1   = esl_randomness_Create(1);
  ESL_RANDOMNESS **rng    = malloc(sizeof(ESL_RANDOMNESS *) * ntr);
  P7_TRACE        *tr1    = p7_trace_fs_Create();
  P7_TRACE       **tr2    = malloc(sizeof(P7_TRACE *) * ntr);
  uint32_t        *seed   = malloc(sizeof(uint32_t) * ntr);
  int              n, Ln, blockR, start, i, j, t, z, k, x;

  for (t = 0; t < ntr; t++) { rng[t] = esl_randomness_Create(1); tr2[t] = p7_trace_fs_Create(); }

  if (p7_hmm_Sample(r, M, abcAA, &hmm) != eslOK) esl_fatal(msg);
  for (k = 1; k <= M; k++)
    {
      x = esl_vec_FArgMax(hmm->mat[k], abcAA->K);
      esl_vec_FSet(hmm->mat[k], abcAA->K, 0.15 / (abcAA->K - 1));
      hmm->mat[k][x] = 0.85;
      if (k < M) { hmm->t[k][p7H_MM] = 0.9; hmm->t[k][p7H_MI] = 0.05; hmm->t[k][p7H_MD] = 0.05; }
    }
  p7_emit_SimpleConsensus(hmm, sq);

  for (n = 0; n < N; n++)
    {
      Ln = 15 + esl_rnd_Roll(r, L - 14);
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, dsq);
      if (n % 2 == 1 && 3*M <= L)
        {
          Ln    = ESL_MAX(Ln, 3*M);
          esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, dsq);
          start = 1 + esl_rnd_Roll(r, Ln - 3*M + 1);
          for (i = 1, j = start; i <= M; i++, j += 3)
            p7_codontable_GetCodon(ct, r, sq->dsq[i], dsq+j);
        }
      blockR = (n % 3 == 2) ? Ln+1 : 5 + esl_rnd_Roll(r, 8);

      if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, Ln/3, p7_LOCAL) != eslOK) esl_fatal(msg);
      if (p7_GForward_Frameshift(dsq, Ln, gm_fs5, fx, iv, NULL)         != eslOK) esl_fatal(msg);
      if (p7_chkgmx_GrowTo(cx, M, Ln, blockR)                           != eslOK) esl_fatal(msg);
      if (p7_GForwardChk_Frameshift(dsq, Ln, gm_fs5, cx, iv, NULL)      != eslOK) esl_fatal(msg);

      for (t = 0; t < ntr; t++) { seed[t] = 1 + esl_rnd_Roll(r, 1000000); esl_randomness_Init(rng[t], seed[t]); }
      if (p7_GStochasticTraceChk_Frameshift(rng, ntr, dsq, Ln, gm_fs5, cx, iv, tr2) != eslOK) esl_fatal(msg);
      for (t = 0; t < ntr; t++)
        {
          esl_randomness_Init(rng1, seed[t]);
          if (p7_GStochasticTrace_Frameshift(rng1, dsq, Ln, gm_fs5, fx, tr1) != eslOK) esl_fatal(msg);
          if (tr1->N != tr2[t]->N) esl_fatal("%s: trace %d lengths differ (%d %d)", msg, t, tr1->N, tr2[t]->N);
          for (z = 0; z < tr1->N; z++)
            if (tr1->st[z] != tr2[t]->st[z] || tr1->k[z] != tr2[t]->k[z] || tr1->i[z] != tr2[t]->i[z] || tr1->c[z] != tr2[t]->c[z])
              esl_fatal("%s: trace %d differs at %d (block size %d)", msg, t, z, blockR);
          p7_trace_Reuse(tr1);
          p7_trace_Reuse(tr2[t]);
        }
    }

  for (t = 0; t < ntr; t++) { esl_randomness_Destroy(rng[t]); p7_trace_fs_Destroy(tr2[t]); }
  free(rng); free(tr2); free(seed); free(dsq);
  esl_randomness_Destroy(rng1);
  p7_trace_fs_Destroy(tr1);
  esl_sq_Destroy(sq);
  p7_ivx_Destroy(iv);
  p7_chkgmx_Destroy(cx);
  p7_gmx_Destroy(fx);
  p7_profile_fs_Destroy(gm_fs5);
  p7_codontable_Destroy(ct);
  p7_hmm_Destroy(hmm);
}

/* utest_stotrace_posterior_frameshift():
 * Stochastic traces must follow the posterior distribution. For the
 * back-translated consensus of a model with sharpened match emissions
 * between random flanks, with two nucleotides deleted and two inserted,
 * and the model's frameshift probability raised to 0.05 so that codon
 * lengths other than 3 occur, <N> traces are sampled with
 * p7_GStochasticTrace_Frameshift(). For each codon length c = 1..5 the
 * mean number of M states with that codon length per trace, and the mean
 * number of I states, must be within 5 standard errors (plus 0.01) of
 * the sum of the matching posterior probabilities from
 * p7_GDecoding_Frameshift(). Every trace must pass p7_trace_fs_Validate().
 */
static void
utest_stotrace_posterior_frameshift(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA, int M, int N)
{
  char           *msg    = "frameshift stochastic trace posterior unit test failed";
  P7_HMM         *hmm    = NULL;
  P7_CODONTABLE  *ct     = p7_codontable_Create(gcode);
  P7_FS_PROFILE  *gm_fs5 = p7_profile_fs_Create(M, abcAA, p7P_5CODONS);
  int             L      = 3*M + 120;
  P7_GMX         *fx     = p7_gmx_Create(M, L, L, p7G_NSCELLS_FS);
  P7_GMX         *bx     = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_IVX         *iv     = p7_ivx_Create(M, p7P_5CODONS);
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq    = malloc(sizeof(ESL_DSQ) * (L+2));
  ESL_DSQ        *tmp    = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_TRACE       *tr     = p7_trace_fs_Create();
  double          expct[6], sum[6], sumsq[6], mean, se;
  char            errbuf[eslERRBUFSIZE];
  int             cnt[6];
  int             i, j, k, x, n, z, c, Ln, start, del1, del2, ins1, ins2;

  if (p7_hmm_Sample(r, M, abcAA, &hmm) != eslOK) esl_fatal(msg);
  for (k = 1; k <= M; k++)
    {
      x = esl_vec_FArgMax(hmm->mat[k], abcAA->K);
      esl_vec_FSet(hmm->mat[k], abcAA->K, 0.15 / (abcAA->K - 1));
      hmm->mat[k][x] = 0.85;
      if (k < M) { hmm->t[k][p7H_MM] = 0.9; hmm->t[k][p7H_MI] = 0.05; hmm->t[k][p7H_MD] = 0.05; }
    }
  hmm->fsprob = 0.05;
  p7_emit_SimpleConsensus(hmm, sq);

  /* consensus between random flanks, then two deletions and two insertions inside it */
  Ln = 3*M + 100;
  esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, tmp);
  start = 1 + esl_rnd_Roll(r, Ln - 3*M + 1);
  for (i = 1, j = start; i <= M; i++, j += 3) p7_codontable_GetCodon(ct, r, sq->dsq[i], tmp+j);
  del1 = start + M/2;
  ins1 = start + M   + esl_rnd_Roll(r, M/4);
  del2 = start + 3*M/2;
  ins2 = start + 2*M + esl_rnd_Roll(r, M/4);
  for (i = 1, j = 1; i <= Ln; i++)
    {
      if (i == del1 || i == del2) continue;
      if (i == ins1 || i == ins2) dsq[j++] = esl_rnd_Roll(r, 4);
      dsq[j++] = tmp[i];
    }
  Ln = j - 1;

  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, Ln/3, p7_LOCAL) != eslOK) esl_fatal(msg);
  if (p7_GForward_Frameshift (dsq, Ln, gm_fs5, fx, iv, NULL) != eslOK) esl_fatal(msg);

  for (n = 0; n < 6; n++) { sum[n] = sumsq[n] = 0.; }
  for (n = 0; n < N; n++)
    {
      if (p7_GStochasticTrace_Frameshift(r, dsq, Ln, gm_fs5, fx, tr) != eslOK) esl_fatal(msg);
      if (p7_trace_fs_Validate(tr, abcDNA, dsq, errbuf)                != eslOK) esl_fatal("%s: invalid trace: %s", msg, errbuf);
      for (c = 0; c < 6; c++) cnt[c] = 0;
      for (z = 0; z < tr->N; z++)
        {
          if      (tr->st[z] == p7T_M) { if (tr->c[z] < 1 || tr->c[z] > 5) esl_fatal("%s: M with codon length %d", msg, tr->c[z]); cnt[tr->c[z]]++; }
          else if (tr->st[z] == p7T_I) cnt[0]++;
        }
      for (c = 0; c < 6; c++) { sum[c] += cnt[c]; sumsq[c] += (double) cnt[c] * cnt[c]; }
      p7_trace_Reuse(tr);
    }

  if (p7_GBackward_Frameshift(dsq, Ln, gm_fs5, bx, iv, NULL) != eslOK) esl_fatal(msg);
  if (p7_GDecoding_Frameshift(gm_fs5, fx, bx)               != eslOK) esl_fatal(msg);
  for (c = 0; c < 6; c++) expct[c] = 0.;
  for (i = 1; i <= Ln; i++)
    for (k = 1; k <= M; k++)
      {
        for (c = 1; c <= 5; c++) expct[c] += fx->dp[i][k*p7G_NSCELLS_FS + p7G_M + c];
        expct[0] += fx->dp[i][k*p7G_NSCELLS_FS + p7G_I];
      }

  for (c = 0; c < 6; c++)
    {
      mean = sum[c] / N;
      se   = sqrt(ESL_MAX(0., sumsq[c] / N - mean * mean) / N);
      if (fabs(mean - expct[c]) > 5. * se + 0.01)
        esl_fatal("%s: %s: mean count %.4f per trace, posterior sum %.4f (s.e. %.4f)", msg, c ? "M codon length" : "I", mean, expct[c], se);
    }

  free(dsq); free(tmp);
  p7_trace_fs_Destroy(tr);
  esl_sq_Destroy(sq);
  p7_ivx_Destroy(iv);
  p7_gmx_Destroy(fx);
  p7_gmx_Destroy(bx);
  p7_profile_fs_Destroy(gm_fs5);
  p7_codontable_Destroy(ct);
  p7_hmm_Destroy(hmm);
}
#endif /*p7GENERIC_FWDBACK_CHK_FRAMESHIFT_TESTDRIVE*/


/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7GENERIC_FWDBACK_CHK_FRAMESHIFT_TESTDRIVE
/*
   gcc -g -Wall -std=gnu99 -o generic_fwdback_chk_frameshift_utest -I. -L. -I../easel -L../easel -Dp7GENERIC_FWDBACK_CHK_FRAMESHIFT_TESTDRIVE generic_fwdback_chk_frameshift.c -lhmmer -leasel -lm
   ./generic_fwdback_chk_frameshift_utest
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
static char banner[] = "unit test driver for the generic checkpointed frameshift Forward and Backward, and stochastic traces from the Forward";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  int             M      = esl_opt_GetInteger(go, "-M");
  int             nseq   = esl_opt_GetInteger(go, "-N");
  ESL_ALPHABET   *abcAA  = esl_alphabet_Create(eslAMINO);
  ESL_ALPHABET   *abcDNA = esl_alphabet_Create(eslDNA);
  ESL_GENCODE    *gcode  = esl_gencode_Create(abcDNA, abcAA);
  P7_BG          *bgAA   = p7_bg_Create(abcAA);
  P7_BG          *bgDNA  = p7_bg_Create(abcDNA);
  P7_HMM         *hmm    = NULL;

  p7_FLogsumInit();

  if (p7_hmm_Sample(r, M, abcAA, &hmm) != eslOK) esl_fatal("failed to sample an HMM");

  utest_fwdback_chk_frameshift(r, abcAA, abcDNA, gcode, bgAA, bgDNA, hmm, nseq);
  utest_stotrace_chk_frameshift(r, abcAA, abcDNA, gcode, bgAA, bgDNA, M, 3*M+150, 10);
  utest_stotrace_posterior_frameshift(r, abcAA, abcDNA, gcode, bgAA, bgDNA, M, 2000);

  p7_hmm_Destroy(hmm);
  p7_bg_Destroy(bgAA);
  p7_bg_Destroy(bgDNA);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcAA);
  esl_alphabet_Destroy(abcDNA);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  printf("All tests passed.\n");
  return eslOK;
}
#endif /*p7GENERIC_FWDBACK_CHK_FRAMESHIFT_TESTDRIVE*/
