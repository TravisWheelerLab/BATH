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

#include "easel.h"
#include "esl_alphabet.h"

#include "hmmer.h"

#define MSC_FS1(k)  (rsc1[(k)])
#define MSC_FS2(k)  (rsc2[(k)])
#define MSC_FS3(k)  (rsc3[(k)])
#define MSC_FS4(k)  (rsc4[(k)])
#define MSC_FS5(k)  (rsc5[(k)])
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

/* nuc_at()
 * The nucleotide at position <j> of <dsq> for codon lookup: the
 * placeholder <p7P_MAXCODONS5> past <L> or for anything but A,C,G,T.
 */
static inline int
nuc_at(const ESL_DSQ *dsq, int L, int j)
{
  return (j <= L && dsq[j] < p7P_MAXNUC) ? dsq[j] : p7P_MAXCODONS5;
}

/* backward_row_frameshift()
 * Frameshift Backward row <i> of profile <gm_fs5> against <dsq> of
 * length <L>: core cells through row pointers <dp>, special states in
 * <xmx>, with <ivx> (M+1 floats) as workspace. Reads rows i+1..i+5
 * (core) and the special states of row i+3. Same calculation as
 * p7_GBackward_Frameshift(), with the codon context of row <i> read
 * from dsq[i+1..i+5] instead of carried from row to row.
 */
static void
backward_row_frameshift(const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs5, float **dp, float *xmx, float *ivx, int i)
{
  float const *tsc = gm_fs5->tsc;
  float const *rsc1, *rsc2, *rsc3, *rsc4, *rsc5;
  int          M   = gm_fs5->M;
  float        esc = p7_fs_profile_IsLocal(gm_fs5) ? 0 : -eslINFINITY;
  int          k;
  int          c1, c2, c3, c4, c5;
  int          t, u, v, w, x;

  x = nuc_at(dsq, L, i+1);
  w = nuc_at(dsq, L, i+2);
  v = nuc_at(dsq, L, i+3);
  u = nuc_at(dsq, L, i+4);
  t = nuc_at(dsq, L, i+5);

  if (i == L)
    {
      XMX(L,p7G_J) = XMX(L,p7G_B) = XMX(L,p7G_N) = -eslINFINITY; /* need to enter and exit model */
      XMX(L,p7G_C) = gm_fs5->xsc[p7P_C][p7P_MOVE];                   /* C<-T         */
      XMX(L,p7G_E) = XMX(L,p7G_C) + gm_fs5->xsc[p7P_E][p7P_MOVE];    /* E<-C, no tail */
      MMX(L,M)     = DMX(L,M) = XMX(L,p7G_E);                    /* {MD}_M <- E (prob 1.0) */
      IMX(L,M)     = -eslINFINITY;                               /* no I_M state        */

      for (k = M-1; k >= 1; k--)
        {
          MMX(L,k) = p7_FLogsum( XMX(L,p7G_E) + esc,
                                 DMX(L, k+1)  + TSC(p7P_MD,k));
          DMX(L,k) = p7_FLogsum( XMX(L,p7G_E) + esc,
                                 DMX(L, k+1)  + TSC(p7P_DD,k));
          IMX(L,k) = -eslINFINITY;
        }
      MMX(L,0) = IMX(L,0) = DMX(L,0)  = -eslINFINITY;
    }
  else if (i >= 1 && i == L-1)
    {
      c1 = p7P_CODON1_FS5(x);
      c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
      rsc1 = gm_fs5->rsc[c1];

      ivx[1]       = MMX(L,1) + MSC_FS1(1);
      XMX(L-1,p7G_B) = ivx[1]   + TSC(p7P_BM,0);
      for (k = 2; k <= M; k++)
        {
          ivx[k] =  MMX(L,k) + MSC_FS1(k);
          XMX(L-1,p7G_B) = p7_FLogsum( XMX(L-1,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
        }

      XMX(L-1,p7G_J) = XMX(L-1,p7G_B) + gm_fs5->xsc[p7P_J][p7P_MOVE];
      XMX(L-1,p7G_N) = XMX(L-1,p7G_B) + gm_fs5->xsc[p7P_N][p7P_MOVE];
      XMX(L-1,p7G_C) = gm_fs5->xsc[p7P_C][p7P_LOOP] + gm_fs5->xsc[p7P_C][p7P_MOVE];
      XMX(L-1,p7G_E) = p7_FLogsum( XMX(L-1,p7G_J) + gm_fs5->xsc[p7P_E][p7P_LOOP],
                                   XMX(L-1,p7G_C) + gm_fs5->xsc[p7P_E][p7P_MOVE]);

      MMX(L-1,M)  = DMX(L-1,M) = XMX(L-1,p7G_E);
      for (k = M-1; k >= 1; k--)
        {
          MMX(L-1,k) = p7_FLogsum( DMX(L-1,k+1)   + TSC(p7P_MD,k),
                       p7_FLogsum( ivx[k+1]       + TSC(p7P_MM,k),
                                   XMX(L-1,p7G_E)  + esc));

          DMX(L-1,k) = p7_FLogsum( p7_FLogsum( XMX(L-1,p7G_E) + esc,
                                               DMX(L-1, k+1)  + TSC(p7P_DD,k)),
                                               ivx[k+1]       + TSC(p7P_DM,k));

          IMX(L-1,k) = ivx[k+1]            + TSC(p7P_IM,k);
        }
      MMX(L-1,0) = IMX(L-1,0) = DMX(L-1,0)  = -eslINFINITY;
    }
  else if (i >= 1 && i == L-2)
    {
      c1 = p7P_CODON1_FS5(x);
      c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
      c2 = p7P_CODON2_FS5(x, w);
      c2 = p7P_MINIDX(c2, p7P_DEGEN5_QC1);
      rsc1 = gm_fs5->rsc[c1];
      rsc2 = gm_fs5->rsc[c2];

      ivx[1] =                     MMX(L-1,1) + MSC_FS1(1);
      ivx[1] = p7_FLogsum( ivx[1], MMX(L,1)   + MSC_FS2(1));

      XMX(L-2,p7G_B)   =  ivx[1] + TSC(p7P_BM,0);
      for (k = 2; k <= M; k++)
        {
          ivx[k] =                     MMX(L-1,k) + MSC_FS1(k);
          ivx[k] = p7_FLogsum( ivx[k], MMX(L,k)   + MSC_FS2(k));
          XMX(L-2,p7G_B) = p7_FLogsum( XMX(L-2,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
        }

      XMX(L-2,p7G_J) = XMX(L-2,p7G_B) + gm_fs5->xsc[p7P_J][p7P_MOVE];
      XMX(L-2,p7G_N) = XMX(L-2,p7G_B) + gm_fs5->xsc[p7P_N][p7P_MOVE];
      XMX(L-2,p7G_C) = gm_fs5->xsc[p7P_C][p7P_LOOP] + gm_fs5->xsc[p7P_C][p7P_MOVE];
      XMX(L-2,p7G_E) = p7_FLogsum( XMX(L-2,p7G_J) + gm_fs5->xsc[p7P_E][p7P_LOOP],
                                   XMX(L-2,p7G_C) + gm_fs5->xsc[p7P_E][p7P_MOVE]);

      MMX(L-2,M)  = DMX(L-2,M) = XMX(L-2,p7G_E);
      for (k = M-1; k >= 1; k--)
        {
          MMX(L-2,k) = p7_FLogsum( DMX(L-2,k+1)   + TSC(p7P_MD,k),
                       p7_FLogsum( ivx[k+1]       + TSC(p7P_MM,k),
                                   XMX(L-2,p7G_E)  + esc));

          DMX(L-2,k) = p7_FLogsum( p7_FLogsum( XMX(L-2,p7G_E) + esc,
                                               DMX(L-2, k+1)  + TSC(p7P_DD,k)),
                                               ivx[k+1]       + TSC(p7P_DM,k));

          IMX(L-2,k) = ivx[k+1]            + TSC(p7P_IM,k);
        }
      MMX(L-2,0) = IMX(L-2,0) = DMX(L-2,0)  = -eslINFINITY;
    }
  else if (i >= 1 && (i == L-3 || i == L-4))
    {
      c1 = p7P_CODON1_FS5(x);
      c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
      c2 = p7P_CODON2_FS5(x, w);
      c2 = p7P_MINIDX(c2, p7P_DEGEN5_QC1);
      c3 = p7P_CODON3_FS5(x, w, v);
      c3 = p7P_MINIDX(c3, p7P_DEGEN5_C);
      c4 = p7P_CODON4_FS5(x, w, v, u);
      c4 = p7P_MINIDX(c4, p7P_DEGEN5_QC1);
      rsc1 = gm_fs5->rsc[c1];
      rsc2 = gm_fs5->rsc[c2];
      rsc3 = gm_fs5->rsc[c3];
      rsc4 = gm_fs5->rsc[c4];

      ivx[1] =                    MMX(i+1,1) + MSC_FS1(1);
      ivx[1] = p7_FLogsum(ivx[1], MMX(i+2,1) + MSC_FS2(1));
      ivx[1] = p7_FLogsum(ivx[1], MMX(i+3,1) + MSC_FS3(1));
      if( i == L-4 )
        ivx[1] = p7_FLogsum(ivx[1], MMX(i+4,1) + MSC_FS4(1));

      XMX(i,p7G_B)   =  ivx[1] + TSC(p7P_BM,0);

      for (k = 2; k <= M; k++)
        {
          ivx[k] =                    MMX(i+1,k) + MSC_FS1(k);
          ivx[k] = p7_FLogsum(ivx[k], MMX(i+2,k) + MSC_FS2(k));
          ivx[k] = p7_FLogsum(ivx[k], MMX(i+3,k) + MSC_FS3(k));
          if( i == L-4 )
            ivx[k] = p7_FLogsum(ivx[k], MMX(i+4,k) + MSC_FS4(k));

          XMX(i,p7G_B) = p7_FLogsum( XMX(i,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
        }

      XMX(i,p7G_J) = p7_FLogsum( XMX(i+3,p7G_J) + gm_fs5->xsc[p7P_J][p7P_LOOP],
                                 XMX(i,  p7G_B) + gm_fs5->xsc[p7P_J][p7P_MOVE]);
      XMX(i,p7G_C) =             XMX(i+3,p7G_C) + gm_fs5->xsc[p7P_C][p7P_LOOP];
      XMX(i,p7G_N) = p7_FLogsum( XMX(i+3,p7G_N) + gm_fs5->xsc[p7P_N][p7P_LOOP],
                                 XMX(i,  p7G_B) + gm_fs5->xsc[p7P_N][p7P_MOVE]);
      XMX(i,p7G_E) = p7_FLogsum(XMX(i,p7G_J) + gm_fs5->xsc[p7P_E][p7P_LOOP],
                                XMX(i,p7G_C) + gm_fs5->xsc[p7P_E][p7P_MOVE]);

      MMX(i,M)     =            DMX(i,M) = XMX(i,p7G_E); /* {MD}_M <- E (prob 1.0) */
      IMX(i,M)     =            -eslINFINITY;            /* no I_M state        */

      for (k = M-1; k >= 1; k--)
        {
          MMX(i,k) = p7_FLogsum( p7_FLogsum( DMX(i,k+1)   + TSC(p7P_MD,k),
                                 p7_FLogsum( IMX(i+3,k)   + TSC(p7P_MI,k),
                                             ivx[k+1]      + TSC(p7P_MM,k))),
                                             XMX(i,p7G_E) + esc);

          DMX(i,k) = p7_FLogsum( p7_FLogsum( XMX(i,p7G_E) + esc,
                                             DMX(i, k+1)  + TSC(p7P_DD,k)),
                                             ivx[k+1]      + TSC(p7P_DM,k));

          IMX(i,k) = p7_FLogsum(           IMX(i+3,k  )   + TSC(p7P_II,k),
                                           ivx[k+1]          + TSC(p7P_IM,k));
        }
      MMX(i,0) = IMX(i,0) = DMX(i,0)  = -eslINFINITY;
    }
  else if (i >= 1) /* 1 <= i <= L-5 */
    {
      c1 = p7P_CODON1_FS5(x);
      c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
      c2 = p7P_CODON2_FS5(x, w);
      c2 = p7P_MINIDX(c2, p7P_DEGEN5_QC1);
      c3 = p7P_CODON3_FS5(x, w, v);
      c3 = p7P_MINIDX(c3, p7P_DEGEN5_C);
      c4 = p7P_CODON4_FS5(x, w, v, u);
      c4 = p7P_MINIDX(c4, p7P_DEGEN5_QC1);
      c5 = p7P_CODON5_FS5(x, w, v, u, t);
      c5 = p7P_MINIDX(c5, p7P_DEGEN5_QC2);
      rsc1 = gm_fs5->rsc[c1];
      rsc2 = gm_fs5->rsc[c2];
      rsc3 = gm_fs5->rsc[c3];
      rsc4 = gm_fs5->rsc[c4];
      rsc5 = gm_fs5->rsc[c5];

      ivx[1] = p7_FLogsum( MMX(i+1,1) + MSC_FS1(1),
               p7_FLogsum( MMX(i+2,1) + MSC_FS2(1),
               p7_FLogsum( MMX(i+3,1) + MSC_FS3(1),
               p7_FLogsum( MMX(i+4,1) + MSC_FS4(1),
                           MMX(i+5,1) + MSC_FS5(1)))));

      XMX(i,p7G_B) = ivx[1] + TSC(p7P_BM,0);

      for (k = 2; k <= M; k++) {
        ivx[k] = p7_FLogsum( MMX(i+1,k) + MSC_FS1(k),
                 p7_FLogsum( MMX(i+2,k) + MSC_FS2(k),
                 p7_FLogsum( MMX(i+3,k) + MSC_FS3(k),
                 p7_FLogsum( MMX(i+4,k) + MSC_FS4(k),
                             MMX(i+5,k) + MSC_FS5(k)))));

        XMX(i,p7G_B) = p7_FLogsum( XMX(i, p7G_B), ivx[k] + TSC(p7P_BM,k-1));
      }

      XMX(i,p7G_J) = p7_FLogsum( XMX(i+3,p7G_J) + gm_fs5->xsc[p7P_J][p7P_LOOP],
                                 XMX(i,  p7G_B) + gm_fs5->xsc[p7P_J][p7P_MOVE]);
      XMX(i,p7G_C) =             XMX(i+3,p7G_C) + gm_fs5->xsc[p7P_C][p7P_LOOP];
      XMX(i,p7G_N) = p7_FLogsum( XMX(i+3,p7G_N) + gm_fs5->xsc[p7P_N][p7P_LOOP],
                                 XMX(i,  p7G_B) + gm_fs5->xsc[p7P_N][p7P_MOVE]);
      XMX(i,p7G_E) = p7_FLogsum(XMX(i,p7G_J) + gm_fs5->xsc[p7P_E][p7P_LOOP],
                                XMX(i,p7G_C) + gm_fs5->xsc[p7P_E][p7P_MOVE]);

      MMX(i,M)     = DMX(i,M) = XMX(i,p7G_E); /* {MD}_M <- E (prob 1.0) */
      IMX(i,M)     = -eslINFINITY;            /* no I_M state        */

      for (k = M-1; k >= 1; k--)
        {
          MMX(i,k) = p7_FLogsum( p7_FLogsum( DMX(i,k+1)   + TSC(p7P_MD,k),
                                 p7_FLogsum( IMX(i+3,k)   + TSC(p7P_MI,k),
                                             ivx[k+1]      + TSC(p7P_MM,k))),
                                             XMX(i,p7G_E) + esc);

          DMX(i,k) = p7_FLogsum( p7_FLogsum( XMX(i,p7G_E) + esc,
                                             DMX(i, k+1)  + TSC(p7P_DD,k)),
                                             ivx[k+1]      + TSC(p7P_DM,k));

          IMX(i,k) = p7_FLogsum(             IMX(i+3,k  ) + TSC(p7P_II,k),
                                             ivx[k+1]      + TSC(p7P_IM,k));
        }
      MMX(i,0) = IMX(i,0) = DMX(i,0)  = -eslINFINITY;
    }
  else /* i == 0: only N,B states are reachable. */
    {
      c1 = p7P_CODON1_FS5(x);
      c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
      c2 = p7P_CODON2_FS5(x, w);
      c2 = p7P_MINIDX(c2, p7P_DEGEN5_QC1);
      c3 = p7P_CODON3_FS5(x, w, v);
      c3 = p7P_MINIDX(c3, p7P_DEGEN5_C);
      c4 = p7P_CODON4_FS5(x, w, v, u);
      c4 = p7P_MINIDX(c4, p7P_DEGEN5_QC1);
      c5 = p7P_CODON5_FS5(x, w, v, u, t);
      c5 = p7P_MINIDX(c5, p7P_DEGEN5_QC2);
      rsc1 = gm_fs5->rsc[c1];
      rsc2 = gm_fs5->rsc[c2];
      rsc3 = gm_fs5->rsc[c3];
      rsc4 = gm_fs5->rsc[c4];
      rsc5 = gm_fs5->rsc[c5];

      /* Rows past L do not exist; when L < 5 they contribute nothing. */
      ivx[1] = p7_FLogsum(               MMX(1,1) + MSC_FS1(1),
               p7_FLogsum( (L >= 2 ? MMX(2,1) + MSC_FS2(1) : -eslINFINITY),
               p7_FLogsum( (L >= 3 ? MMX(3,1) + MSC_FS3(1) : -eslINFINITY),
               p7_FLogsum( (L >= 4 ? MMX(4,1) + MSC_FS4(1) : -eslINFINITY),
                           (L >= 5 ? MMX(5,1) + MSC_FS5(1) : -eslINFINITY)))));

      XMX(0,p7G_B) = ivx[1] + TSC(p7P_BM,0);

      for (k = 2; k <= M; k++)
        {
          ivx[k] =  p7_FLogsum(               MMX(1,k) + MSC_FS1(k),
                    p7_FLogsum( (L >= 2 ? MMX(2,k) + MSC_FS2(k) : -eslINFINITY),
                    p7_FLogsum( (L >= 3 ? MMX(3,k) + MSC_FS3(k) : -eslINFINITY),
                    p7_FLogsum( (L >= 4 ? MMX(4,k) + MSC_FS4(k) : -eslINFINITY),
                                (L >= 5 ? MMX(5,k) + MSC_FS5(k) : -eslINFINITY)))));

          XMX(0,p7G_B) = p7_FLogsum(XMX(0, p7G_B), ivx[k] + TSC(p7P_BM,k-1));
        }

      XMX(0,p7G_J) = XMX(0,p7G_C) = XMX(0,p7G_E) = -eslINFINITY;

      XMX(0,p7G_N) = p7_FLogsum( (L >= 3 ? XMX(3,p7G_N) + gm_fs5->xsc[p7P_N][p7P_LOOP] : -eslINFINITY),
                                 XMX(0,  p7G_B) + gm_fs5->xsc[p7P_N][p7P_MOVE]);

      for (k = M; k >= 0; k--)
        MMX(0,k) = DMX(0,k) =  IMX(0,k) = -eslINFINITY;
    }
}


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
        backward_row_frameshift(dsq, L, gm_fs5, dp, xmx, iv->ivx, i);
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
    backward_row_frameshift(dsq, cx->L, gm_fs5, cx->dp, cx->xmx, iv->ivx, i);
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

/* compare_block_fs(), compare_specials_fs()
 * The same for the Forward matrices, which set every cell, so all
 * cells k = 0..M of the loaded block's rows are compared.
 */
static void
compare_block_fs(const char *msg, const P7_CHKGMX *cx, const P7_GMX *gx, int b)
{
  int i, j;

  for (i = b; i <= ESL_MIN(cx->L, b + cx->blockR - 1); i++)
    for (j = 0; j < (cx->M+1) * p7G_NSCELLS_FS; j++)
      if (cx->dp[i][j] != gx->dp[i][j])
        esl_fatal("%s: Forward row %d, k %d, cell %d differs (block size %d)", msg, i, j / p7G_NSCELLS_FS, j % p7G_NSCELLS_FS, cx->blockR);
}

static void
compare_specials_fs(const char *msg, const P7_CHKGMX *cx, const P7_GMX *gx)
{
  int i, x;

  for (i = 0; i <= cx->L; i++)
    for (x = 0; x < p7G_NXCELLS; x++)
      if (cx->xmx[i*p7G_NXCELLS + x] != gx->xmx[i*p7G_NXCELLS + x])
        esl_fatal("%s: Forward special state %d of row %d differs (block size %d)", msg, x, i, cx->blockR);
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
            compare_specials(msg, cx, gx);
            compare_block(msg, cx, gx, 0);

            for (b = 0; b <= L; b += cx->blockR)
              {
                if (p7_GBackwardChkLoad_Frameshift(dsq, gm_fs5, cx, iv, b) != eslOK) esl_fatal(msg);
                compare_block(msg, cx, gx, b);
              }

            if (p7_chkgmx_GrowTo(fcx, hmm->M, L, blocks[q])                != eslOK) esl_fatal(msg);
            if (p7_GForwardChk_Frameshift(dsq, L, gm_fs5, fcx, iv, &chk_sc) != eslOK) esl_fatal(msg);
            if (chk_sc != fsc) esl_fatal("%s: Forward scores differ (%f %f), L %d, block size %d", msg, fsc, chk_sc, L, blocks[q]);
            compare_specials_fs(msg, fcx, fx);
            compare_block_fs(msg, fcx, fx, (L / fcx->blockR) * fcx->blockR);

            for (b = (L / fcx->blockR) * fcx->blockR; b >= 0; b -= fcx->blockR)
              {
                if (p7_GForwardChkLoad_Frameshift(dsq, gm_fs5, fcx, iv, b) != eslOK) esl_fatal(msg);
                compare_block_fs(msg, fcx, fx, b);
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
static char banner[] = "unit test driver for the generic checkpointed frameshift Forward and Backward";

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
