/* Forward/Backward Frameshift algorithms, full matrices.
 *
 * Contents:
 *   1. Forward, 5 and 3 codon lengths.
 *   2. Backward, 5 and 3 codon lengths.
 *   3. Benchmark driver.
 *   4. Unit tests.
 *   5. Test driver.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"

#include "hmmer.h"

#define IVX3(i,k) (ivx[((k)*p7P_3CODONS) + (i)])
#define IVX5(i,k) (ivx[((k)*p7P_5CODONS) + (i)])

/* Emission score accessors: rscN must be set to gm->rsc[cN] before the k-loop */
#define MSC_FS1(k)  (rsc1[(k)])
#define MSC_FS2(k)  (rsc2[(k)])
#define MSC_FS3(k)  (rsc3[(k)])
#define MSC_FS4(k)  (rsc4[(k)])
#define MSC_FS5(k)  (rsc5[(k)])

/*****************************************************************
 * 1. Forward, 5 and 3 codon lengths.
 *****************************************************************/

/* nt_at_fs5()
 * The nucleotide at position <j> of <dsq> for codon lookup: the
 * placeholder <p7P_MAXCODONS5> before position 1 or for anything but
 * A,C,G,T.
 */
static inline int
nt_at_fs5(const ESL_DSQ *dsq, int j)
{
  return (j >= 1 && dsq[j] < p7P_MAXNUC) ? dsq[j] : p7P_MAXCODONS5;
}

/* Function:  p7_GForwardRow_Frameshift_5Codons()
 * Synopsis:  One row of the frameshift Forward algorithm.
 *
 * Purpose:   Forward row <i> (0..L) of the 5-codon profile <gm_fs5>
 *            against the DNA sequence <dsq> into <dp>[i] and row <i> of
 *            <xmx>, from rows i-1 and i-3 (core cells and special
 *            states) and the codon partial sums of rows i-1..i-4 in
 *            <ivx> (<p7P_5CODONS> rows of a P7_IVX), where row <i>'s own
 *            partial sums are written. The codon context of row <i> is
 *            read from dsq[i-4..i]. Row 0 is the initialization, of
 *            <dp>[0], row 0 of <xmx> and all of <ivx>. Same calculation
 *            as each row of <p7_GForward_Frameshift_5Codons()>, which calls it.
 *
 *            Caller must have initialized the log-sum calculation with
 *            a call to <p7_FLogsumInit()>.
 */
void
p7_GForwardRow_Frameshift_5Codons(const ESL_DSQ *dsq, const P7_FS_PROFILE *gm_fs5, float **dp, float *xmx, float *ivx, int i)
{
  float const *tsc  = gm_fs5->tsc;
  float const *rsc1, *rsc2, *rsc3, *rsc4, *rsc5;
  int          M    = gm_fs5->M;
  int          k;
  int          c1, c2, c3, c4, c5;
  int          t, u, v, w, x;
  int          ivx_1, ivx_2, ivx_3, ivx_4, ivx_5;
  float        esc  = p7_fs_profile_IsLocal(gm_fs5) ? 0 : -eslINFINITY;

  if (i == 0)
    {
      for (k = 0; k <= M; k++)
        for (c1 = 0; c1 < p7P_5CODONS; c1++)
          IVX5(c1,k) = -eslINFINITY;
      XMX_FS(0,p7G_N) = 0.;
      XMX_FS(0,p7G_B) = gm_fs5->xsc[p7P_N][p7P_MOVE];
      XMX_FS(0,p7G_E) = XMX_FS(0,p7G_J) = XMX_FS(0,p7G_C) = -eslINFINITY;
      for (k = 0; k <= M; k++) {
        MMX_FS(0,k,p7G_C1) = MMX_FS(0,k,p7G_C2) = MMX_FS(0,k,p7G_C3) = MMX_FS(0,k,p7G_C4) = MMX_FS(0,k,p7G_C5) = -eslINFINITY;
        MMX_FS(0,k,p7G_C0) = IMX_FS(0,k) = DMX_FS(0,k) = -eslINFINITY;
      }
      return;
    }

  x = nt_at_fs5(dsq, i);
  w = nt_at_fs5(dsq, i-1);
  v = nt_at_fs5(dsq, i-2);
  u = nt_at_fs5(dsq, i-3);
  t = nt_at_fs5(dsq, i-4);

  if (i == 1)
    {
      XMX_FS(1,p7G_N) = 0.;
      XMX_FS(1,p7G_B) = gm_fs5->xsc[p7P_N][p7P_MOVE];
      XMX_FS(1,p7G_E) = -eslINFINITY;
      MMX_FS(1,0,p7G_C1) = MMX_FS(1,0,p7G_C2) = MMX_FS(1,0,p7G_C3) = MMX_FS(1,0,p7G_C4) = MMX_FS(1,0,p7G_C5) = -eslINFINITY;
      MMX_FS(1,0,p7G_C0) = IMX_FS(1,0) = DMX_FS(1,0) = -eslINFINITY;

      c1 = p7P_CODON1_FS5(x);
      c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
      rsc1 = gm_fs5->rsc[c1];
      for (k = 1; k <= M; k++) {
        IVX5(1,k) = XMX_FS(0,p7G_B) + TSC(p7P_BM,k-1);

        MMX_FS(1,k,p7G_C1) = IVX5(1,k) + MSC_FS1(k);
        MMX_FS(1,k,p7G_C2) = -eslINFINITY;
        MMX_FS(1,k,p7G_C3) = -eslINFINITY;
        MMX_FS(1,k,p7G_C4) = -eslINFINITY;
        MMX_FS(1,k,p7G_C5) = -eslINFINITY;
        MMX_FS(1,k,p7G_C0) = MMX_FS(1,k,p7G_C1);
        IMX_FS(1,k)  = -eslINFINITY;
        DMX_FS(1,k)  = p7_FLogsum(MMX_FS(1,k-1,p7G_C0) + TSC(p7P_MD,k-1),
                                  DMX_FS(1,k-1)        + TSC(p7P_DD,k-1));

        XMX_FS(1,p7G_E) = p7_FLogsum(MMX_FS(1,k,p7G_C0) + esc,
                          p7_FLogsum(DMX_FS(1,k) + esc,
                                     XMX_FS(1,p7G_E)));
      }
      XMX_FS(1,p7G_J) = XMX_FS(1,p7G_E) + gm_fs5->xsc[p7P_E][p7P_LOOP];
      XMX_FS(1,p7G_C) = XMX_FS(1,p7G_E) + gm_fs5->xsc[p7P_E][p7P_MOVE];
      return;
    }

  if (i == 2)
    {
      XMX_FS(2,p7G_N) = 0.;
      XMX_FS(2,p7G_B) = gm_fs5->xsc[p7P_N][p7P_MOVE];
      XMX_FS(2,p7G_E) = -eslINFINITY;
      MMX_FS(2,0,p7G_C1) = MMX_FS(2,0,p7G_C2) = MMX_FS(2,0,p7G_C3) = MMX_FS(2,0,p7G_C4) = MMX_FS(2,0,p7G_C5) = -eslINFINITY;
      MMX_FS(2,0,p7G_C0) = IMX_FS(2,0) = DMX_FS(2,0) = -eslINFINITY;

      c1 = p7P_CODON1_FS5(x);
      c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
      c2 = p7P_CODON2_FS5(w, x);
      c2 = p7P_MINIDX(c2, p7P_DEGEN5_QC1);
      rsc1 = gm_fs5->rsc[c1];
      rsc2 = gm_fs5->rsc[c2];

      for (k = 1; k <= M; k++) {
        IVX5(2,k) = p7_FLogsum(MMX_FS(1,k-1,p7G_C0) + TSC(p7P_MM,k-1),
                    p7_FLogsum(IMX_FS(1,k-1)        + TSC(p7P_IM,k-1),
                    p7_FLogsum(DMX_FS(1,k-1)        + TSC(p7P_DM,k-1),
                               XMX_FS(1,p7G_B)      + TSC(p7P_BM,k-1))));
        MMX_FS(2,k,p7G_C1) = IVX5(2,k) + MSC_FS1(k);
        MMX_FS(2,k,p7G_C2) = IVX5(1,k) + MSC_FS2(k);
        MMX_FS(2,k,p7G_C3) = -eslINFINITY;
        MMX_FS(2,k,p7G_C4) = -eslINFINITY;
        MMX_FS(2,k,p7G_C5) = -eslINFINITY;
        MMX_FS(2,k,p7G_C0) = p7_FLogsum( MMX_FS(2,k,p7G_C1), MMX_FS(2,k,p7G_C2));
        IMX_FS(2,k)  = -eslINFINITY;
        DMX_FS(2,k)  = p7_FLogsum(MMX_FS(2,k-1,p7G_C0) + TSC(p7P_MD,k-1),
                                  DMX_FS(2,k-1)        + TSC(p7P_DD,k-1));

        XMX_FS(2,p7G_E) = p7_FLogsum(MMX_FS(2,k,p7G_C0) + esc,
                          p7_FLogsum(DMX_FS(2,k)        + esc,
                                     XMX_FS(2,p7G_E)));
      }
      XMX_FS(2,p7G_J) = XMX_FS(2,p7G_E) + gm_fs5->xsc[p7P_E][p7P_LOOP];
      XMX_FS(2,p7G_C) = XMX_FS(2,p7G_E) + gm_fs5->xsc[p7P_E][p7P_MOVE];
      return;
    }

  c1 = p7P_CODON1_FS5(x);
  c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
  c2 = p7P_CODON2_FS5(w, x);
  c2 = p7P_MINIDX(c2, p7P_DEGEN5_QC1);
  c3 = p7P_CODON3_FS5(v, w, x);
  c3 = p7P_MINIDX(c3, p7P_DEGEN5_C);
  c4 = p7P_CODON4_FS5(u, v, w, x);
  c4 = p7P_MINIDX(c4, p7P_DEGEN5_QC1);
  rsc1 = gm_fs5->rsc[c1];
  rsc2 = gm_fs5->rsc[c2];
  rsc3 = gm_fs5->rsc[c3];
  rsc4 = gm_fs5->rsc[c4];

  ivx_1 = i     % p7P_5CODONS;
  ivx_2 = (i-1) % p7P_5CODONS;
  ivx_3 = (i-2) % p7P_5CODONS;
  ivx_4 = (i-3) % p7P_5CODONS;

  if (i < 5)   /* rows 3 and 4 */
    {
      MMX_FS(i,0,p7G_C1) = MMX_FS(i,0,p7G_C2) = MMX_FS(i,0,p7G_C3) = MMX_FS(i,0,p7G_C4) = MMX_FS(i,0,p7G_C5) = -eslINFINITY;
      MMX_FS(i,0,p7G_C0) = IMX_FS(i,0) = DMX_FS(i,0) = -eslINFINITY;

      XMX_FS(i,p7G_E) = -eslINFINITY;

      for (k = 1; k < M; k++)
        {
          IVX5(ivx_1,k) = p7_FLogsum(MMX_FS(i-1,k-1,p7G_C0) + TSC(p7P_MM,k-1),
                          p7_FLogsum(IMX_FS(i-1,k-1)        + TSC(p7P_IM,k-1),
                          p7_FLogsum(DMX_FS(i-1,k-1)        + TSC(p7P_DM,k-1),
                                     XMX_FS(i-1,p7G_B)      + TSC(p7P_BM,k-1))));

          MMX_FS(i,k,p7G_C1) = IVX5(ivx_1,k) + MSC_FS1(k);
          MMX_FS(i,k,p7G_C2) = IVX5(ivx_2,k) + MSC_FS2(k);
          MMX_FS(i,k,p7G_C3) = IVX5(ivx_3,k) + MSC_FS3(k);
          MMX_FS(i,k,p7G_C4) = -eslINFINITY;
          MMX_FS(i,k,p7G_C5) = -eslINFINITY;
          if( i == 4 )
            MMX_FS(i,k,p7G_C4) = IVX5(ivx_4,k) + MSC_FS4(k);

          MMX_FS(i,k,p7G_C0) = p7_FLogsum( MMX_FS(i,k,p7G_C1),
                               p7_FLogsum( MMX_FS(i,k,p7G_C2),
                               p7_FLogsum( MMX_FS(i,k,p7G_C3),
                                           MMX_FS(i,k,p7G_C4))));

          IMX_FS(i,k) = p7_FLogsum(MMX_FS(i-3,k,p7G_C0) + TSC(p7P_MI,k),
                                   IMX_FS(i-3,k)        + TSC(p7P_II,k));

          DMX_FS(i,k) = p7_FLogsum(MMX_FS(i,k-1,p7G_C0) + TSC(p7P_MD,k-1),
                                   DMX_FS(i,k-1)        + TSC(p7P_DD,k-1));

          XMX_FS(i,p7G_E) = p7_FLogsum(MMX_FS(i,k,p7G_C0) + esc,
                            p7_FLogsum(DMX_FS(i,k)        + esc,
                                       XMX_FS(i,p7G_E)));
        }

      IVX5(ivx_1,M) = p7_FLogsum(MMX_FS(i-1,M-1,p7G_C0) + TSC(p7P_MM,M-1),
                      p7_FLogsum(IMX_FS(i-1,M-1)        + TSC(p7P_IM,M-1),
                      p7_FLogsum(DMX_FS(i-1,M-1)        + TSC(p7P_DM,M-1),
                                 XMX_FS(i-1,p7G_B)      + TSC(p7P_BM,M-1))));

      MMX_FS(i,M,p7G_C1) = IVX5(ivx_1,M) + MSC_FS1(M);
      MMX_FS(i,M,p7G_C2) = IVX5(ivx_2,M) + MSC_FS2(M);
      MMX_FS(i,M,p7G_C3) = IVX5(ivx_3,M) + MSC_FS3(M);
      MMX_FS(i,M,p7G_C4) = -eslINFINITY;
      MMX_FS(i,M,p7G_C5) = -eslINFINITY;
      if( i == 4 )
        MMX_FS(i,M,p7G_C4) = IVX5(ivx_4,M) + MSC_FS4(M);

      MMX_FS(i,M,p7G_C0) = p7_FLogsum( MMX_FS(i,M,p7G_C1),
                           p7_FLogsum( MMX_FS(i,M,p7G_C2),
                           p7_FLogsum( MMX_FS(i,M,p7G_C3),
                                       MMX_FS(i,M,p7G_C4))));

      IMX_FS(i,M) = -eslINFINITY;

      DMX_FS(i,M) = p7_FLogsum(MMX_FS(i,M-1,p7G_C0) + TSC(p7P_MD,M-1),
                               DMX_FS(i,M-1)        + TSC(p7P_DD,M-1));

      XMX_FS(i,p7G_E) = p7_FLogsum(MMX_FS(i,M,p7G_C0),
                        p7_FLogsum(DMX_FS(i,M),
                                   XMX_FS(i,p7G_E)));

      XMX_FS(i,p7G_J) = p7_FLogsum(XMX_FS(i-3,p7G_J) + gm_fs5->xsc[p7P_J][p7P_LOOP],
                                   XMX_FS(i,p7G_E)   + gm_fs5->xsc[p7P_E][p7P_LOOP]);
      XMX_FS(i,p7G_C) = p7_FLogsum(XMX_FS(i-3,p7G_C) + gm_fs5->xsc[p7P_C][p7P_LOOP],
                                   XMX_FS(i,p7G_E)   + gm_fs5->xsc[p7P_E][p7P_MOVE]);
      XMX_FS(i,p7G_N) =            XMX_FS(i-3,p7G_N) + gm_fs5->xsc[p7P_N][p7P_LOOP];
      XMX_FS(i,p7G_B) = p7_FLogsum(XMX_FS(i,p7G_N) + gm_fs5->xsc[p7P_N][p7P_MOVE],
                                   XMX_FS(i,p7G_J) + gm_fs5->xsc[p7P_J][p7P_MOVE]);
      return;
    }

  /* rows 5..L */
  c5 = p7P_CODON5_FS5(t, u, v, w, x);
  c5 = p7P_MINIDX(c5, p7P_DEGEN5_QC2);
  rsc5 = gm_fs5->rsc[c5];
  ivx_5 = (i-4) % p7P_5CODONS;

  MMX_FS(i,0,p7G_C0) = MMX_FS(i,0,p7G_C1) = MMX_FS(i,0,p7G_C2) = MMX_FS(i,0,p7G_C3)
  = MMX_FS(i,0,p7G_C4) = MMX_FS(i,0,p7G_C5) = IMX_FS(i,0) = DMX_FS(i,0) = -eslINFINITY;

  XMX_FS(i, p7G_E) = -eslINFINITY;

  for (k = 1; k < M; k++)
    {
      IVX5(ivx_1,k) = p7_FLogsum(MMX_FS(i-1,k-1,p7G_C0)   + TSC(p7P_MM,k-1),
                      p7_FLogsum(IMX_FS(i-1,k-1)          + TSC(p7P_IM,k-1),
                      p7_FLogsum(DMX_FS(i-1,k-1)          + TSC(p7P_DM,k-1),
                                 XMX_FS(i-1,p7G_B)        + TSC(p7P_BM,k-1))));

      MMX_FS(i,k,p7G_C1) = IVX5(ivx_1,k) + MSC_FS1(k);
      MMX_FS(i,k,p7G_C2) = IVX5(ivx_2,k) + MSC_FS2(k);
      MMX_FS(i,k,p7G_C3) = IVX5(ivx_3,k) + MSC_FS3(k);
      MMX_FS(i,k,p7G_C4) = IVX5(ivx_4,k) + MSC_FS4(k);
      MMX_FS(i,k,p7G_C5) = IVX5(ivx_5,k) + MSC_FS5(k);

      MMX_FS(i,k,p7G_C0) =  p7_FLogsum(p7_FLogsum(MMX_FS(i,k,p7G_C1),
                            p7_FLogsum(MMX_FS(i,k,p7G_C2), MMX_FS(i,k,p7G_C3))),
                            p7_FLogsum(MMX_FS(i,k,p7G_C4), MMX_FS(i,k,p7G_C5)));

      IMX_FS(i,k) = p7_FLogsum(MMX_FS(i-3,k,p7G_C0) + TSC(p7P_MI,k),
                               IMX_FS(i-3,k)        + TSC(p7P_II,k));

      DMX_FS(i,k) = p7_FLogsum(MMX_FS(i,k-1,p7G_C0) + TSC(p7P_MD,k-1),
                               DMX_FS(i,k-1)        + TSC(p7P_DD,k-1));

      XMX_FS(i,p7G_E) = p7_FLogsum(MMX_FS(i,k,p7G_C0) + esc,
                        p7_FLogsum(DMX_FS(i,k)        + esc,
                                   XMX_FS(i,p7G_E)));
    }

  IVX5(ivx_1,M) = p7_FLogsum(MMX_FS(i-1,M-1,p7G_C0)   + TSC(p7P_MM,M-1),
                  p7_FLogsum(IMX_FS(i-1,M-1)          + TSC(p7P_IM,M-1),
                  p7_FLogsum(DMX_FS(i-1,M-1)          + TSC(p7P_DM,M-1),
                             XMX_FS(i-1,p7G_B)        + TSC(p7P_BM,M-1))));

  MMX_FS(i,M,p7G_C1) = IVX5(ivx_1,M) + MSC_FS1(M);
  MMX_FS(i,M,p7G_C2) = IVX5(ivx_2,M) + MSC_FS2(M);
  MMX_FS(i,M,p7G_C3) = IVX5(ivx_3,M) + MSC_FS3(M);
  MMX_FS(i,M,p7G_C4) = IVX5(ivx_4,M) + MSC_FS4(M);
  MMX_FS(i,M,p7G_C5) = IVX5(ivx_5,M) + MSC_FS5(M);

  MMX_FS(i,M,p7G_C0) =  p7_FLogsum(p7_FLogsum(MMX_FS(i,M,p7G_C1),
                        p7_FLogsum(MMX_FS(i,M,p7G_C2), MMX_FS(i,M,p7G_C3))),
                        p7_FLogsum(MMX_FS(i,M,p7G_C4), MMX_FS(i,M,p7G_C5)));

  IMX_FS(i,M) = -eslINFINITY;

  DMX_FS(i,M) = p7_FLogsum(MMX_FS(i,M-1,p7G_C0) + TSC(p7P_MD,M-1),
                           DMX_FS(i,M-1) + TSC(p7P_DD,M-1));

  XMX_FS(i,p7G_E) = p7_FLogsum(p7_FLogsum(MMX_FS(i,M,p7G_C0),
                                          DMX_FS(i,M)),
                                          XMX_FS(i,p7G_E));

  XMX_FS(i,p7G_J) = p7_FLogsum(XMX_FS(i-3,p7G_J) + gm_fs5->xsc[p7P_J][p7P_LOOP],
                               XMX_FS(i,p7G_E)   + gm_fs5->xsc[p7P_E][p7P_LOOP]);
  XMX_FS(i,p7G_C) = p7_FLogsum(XMX_FS(i-3,p7G_C) + gm_fs5->xsc[p7P_C][p7P_LOOP],
                               XMX_FS(i,p7G_E)   + gm_fs5->xsc[p7P_E][p7P_MOVE]);
  XMX_FS(i,p7G_N) =            XMX_FS(i-3,p7G_N) + gm_fs5->xsc[p7P_N][p7P_LOOP];
  XMX_FS(i,p7G_B) = p7_FLogsum(XMX_FS(i,p7G_N) + gm_fs5->xsc[p7P_N][p7P_MOVE],
                               XMX_FS(i,p7G_J) + gm_fs5->xsc[p7P_J][p7P_MOVE]);
}


/* Function:  p7_GForward_Frameshift_5Codons()
 * Synopsis:  The frameshift Forward algorithm, 5 codon lengths.
 *
 * Purpose:   The Forward dynamic programming algorithm for frameshift
 *            aware translated comparison of the DNA sequence <dsq> of
 *            length <L> to the profile <gm_fs5>, with five codon
 *            lengths (1 to 5 nucleotides).
 *
 *            Caller provides the DP matrix <gx>, allocated for at least
 *            <gm_fs5->M> by <L> with <p7G_NSCELLS_FS> cells per position
 *            (M_C0..M_C5, I, D; M_C1..M_C5 are the M state by codon
 *            length, M_C0 their sum), and the intermediate value matrix
 *            <iv> with <p7P_5CODONS> rows. Each row is computed by
 *            <p7_GForwardRow_Frameshift_5Codons()>. The Forward matrix
 *            is returned in <gx> and the Forward score in <opt_sc>.
 *
 *            The Forward score is in lod score form. To convert to a
 *            bitscore, the caller needs to subtract a null model lod
 *            score, then convert to bits.
 *
 *            Caller must have initialized the log-sum calculation
 *            with a call to <p7_FLogsumInit()>.
 *
 * Args:      dsq    - nucleotide sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gm_fs5 - frameshift aware codon profile, 5 codon lengths
 *            gx     - DP matrix with room for an MxL alignment
 *            iv     - intermediate value matrix
 *            opt_sc - optRETURN: Forward lod score in nats
 *
 * Return:    <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <gm_fs5> isn't for 5 codon lengths or <gx>
 *            doesn't have <p7G_NSCELLS_FS> cells per position.
 */
int
p7_GForward_Frameshift_5Codons(const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs5, P7_GMX *gx, P7_IVX *iv, float *opt_sc)
{
  float **dp  = gx->dp;
  float  *xmx = gx->xmx;
  int     i;

  if (gm_fs5->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");
  if (gx->nscells != p7G_NSCELLS_FS) ESL_EXCEPTION(eslEINVAL, "DP matrix not allocated with p7G_NSCELLS_FS cells");

  for (i = 0; i <= L; i++)
    p7_GForwardRow_Frameshift_5Codons(dsq, gm_fs5, dp, xmx, iv->ivx, i);

  if (opt_sc != NULL) *opt_sc = p7_FLogsum( XMX_FS(L,p7G_C),
                                p7_FLogsum( XMX_FS(L-1,p7G_C) + gm_fs5->xsc[p7P_C][p7P_LOOP],
                                            (L >= 2 ? XMX_FS(L-2,p7G_C) + gm_fs5->xsc[p7P_C][p7P_LOOP] : -eslINFINITY))) +
                                            gm_fs5->xsc[p7P_C][p7P_MOVE];

  gx->M = gm_fs5->M;
  gx->L = L;
  return eslOK;
}


/* nt_at_fs3()
 * The nucleotide at position <j> of <dsq> for 3-codon lookup: the
 * placeholder <p7P_MAXCODONS3> before position 1 or for anything but
 * A,C,G,T.
 */
static inline int
nt_at_fs3(const ESL_DSQ *dsq, int j)
{
  return (j >= 1 && dsq[j] < p7P_MAXNUC) ? dsq[j] : p7P_MAXCODONS3;
}

/* Function:  p7_GForwardRow_Frameshift_3Codons()
 * Synopsis:  One row of the frameshift Forward algorithm, 3 codon lengths.
 *
 * Purpose:   Forward row <i> (0..L) of the 3-codon profile <gm_fs3>
 *            against the DNA sequence <dsq> into <dp>[i] (<p7G_NSCELLS>
 *            cells per position) and row <i> of <xmx>, from rows i-2
 *            and i-3 (core cells and special states) and the codon
 *            partial sums of rows i-1 and i-2 in <ivx> (<p7P_3CODONS>
 *            rows of a P7_IVX), where row <i>'s own partial sums are
 *            written. The codon context of row <i> is read from
 *            dsq[i-3..i]. Row 0 is the initialization, of <dp>[0], row
 *            0 of <xmx> and all of <ivx>. Same calculation as each row
 *            of <p7_GForward_Frameshift_3Codons()>, which calls it.
 *
 *            Caller must have initialized the log-sum calculation with
 *            a call to <p7_FLogsumInit()>.
 */
void
p7_GForwardRow_Frameshift_3Codons(const ESL_DSQ *dsq, const P7_FS_PROFILE *gm_fs3, float **dp, float *xmx, float *ivx, int i)
{
  float const *tsc  = gm_fs3->tsc;
  float const *rsc2, *rsc3, *rsc4;
  int          M    = gm_fs3->M;
  int          k;
  int          c2, c3, c4;
  int          u, v, w, x;
  int          ivx_2, ivx_3, ivx_4;
  float        esc  = p7_fs_profile_IsLocal(gm_fs3) ? 0 : -eslINFINITY;

  if (i == 0)
    for (k = 0; k <= M; k++)
      for (c2 = 0; c2 < p7P_3CODONS; c2++)
        IVX3(c2,k) = -eslINFINITY;

  if (i <= 1)
    {
      XMX(i,p7G_N) = 0.;                           /* S->N, p=1          */
      XMX(i,p7G_B) = gm_fs3->xsc[p7P_N][p7P_MOVE]; /* S->N->B, no N-tail */
      XMX(i,p7G_E) = XMX(i,p7G_J) = XMX(i,p7G_C) = -eslINFINITY;
      for (k = 0; k <= M; k++)
        MMX(i,k) = IMX(i,k) = DMX(i,k) = -eslINFINITY;
      return;
    }

  x = nt_at_fs3(dsq, i);
  w = nt_at_fs3(dsq, i-1);
  v = nt_at_fs3(dsq, i-2);
  u = nt_at_fs3(dsq, i-3);

  /* Length 2 quasi-codon */
  c2   = p7P_CODON2_FS3(w, x);
  c2   = p7P_MINIDX(c2, p7P_DEGEN3_QC1);
  rsc2 = gm_fs3->rsc[c2];

  if (i == 2)
    {
      XMX(2,p7G_E) = -eslINFINITY;
      MMX(2,0) = IMX(2,0) = DMX(2,0) = -eslINFINITY;

      for (k = 1; k <= M; k++)
        {
          IVX3(2,k) = XMX(0,p7G_B) + TSC(p7P_BM,k-1);
          MMX(2,k)  = IVX3(2,k) + MSC_FS2(k);
          IMX(2,k)  = -eslINFINITY;
          DMX(2,k)  = p7_FLogsum(MMX(2,k-1) + TSC(p7P_MD,k-1),
                                 DMX(2,k-1) + TSC(p7P_DD,k-1));
          XMX(2,p7G_E) = p7_FLogsum(MMX(2,k) + esc,
                         p7_FLogsum(DMX(2,k) + esc,
                                    XMX(2,p7G_E)));
        }

      XMX(2,p7G_N) = 0.; /* S->N, p=1 */
      XMX(2,p7G_J) = XMX(2,p7G_E) + gm_fs3->xsc[p7P_E][p7P_LOOP];
      XMX(2,p7G_C) = XMX(2,p7G_E) + gm_fs3->xsc[p7P_E][p7P_MOVE];
      XMX(2,p7G_B) = p7_FLogsum(XMX(2,p7G_N) + gm_fs3->xsc[p7P_N][p7P_MOVE],
                                XMX(2,p7G_J) + gm_fs3->xsc[p7P_J][p7P_MOVE]);
      return;
    }

  /* Length 3 codon */
  c3 = p7P_CODON3_FS3(v, w, x);
  c3 = p7P_MINIDX(c3, p7P_DEGEN3_C);

  /* Length 4 quasi-codon */
  c4 = p7P_CODON4_FS3(u, v, w, x);
  c4 = p7P_MINIDX(c4, p7P_DEGEN3_QC1);
  rsc3 = gm_fs3->rsc[c3];
  rsc4 = gm_fs3->rsc[c4];

  ivx_2 = i     % p7P_3CODONS;
  ivx_3 = (i-1) % p7P_3CODONS;
  ivx_4 = (i-2) % p7P_3CODONS;

  MMX(i,0) = IMX(i,0) = DMX(i,0) = -eslINFINITY;
  XMX(i,p7G_E) = -eslINFINITY;

  for (k = 1; k < M; k++)
    {
      /* Only the partial sum for the length 2 codon is new in each row;
       * those of rows i-1 and i-2 serve the length 3 and 4 codons. */
      IVX3(ivx_2,k) = p7_FLogsum(MMX(i-2,k-1)   + TSC(p7P_MM,k-1),
                      p7_FLogsum(IMX(i-2,k-1)   + TSC(p7P_IM,k-1),
                      p7_FLogsum(DMX(i-2,k-1)   + TSC(p7P_DM,k-1),
                                 XMX(i-2,p7G_B) + TSC(p7P_BM,k-1))));

      MMX(i,k) =                      IVX3(ivx_2,k) + MSC_FS2(k);
      MMX(i,k) = p7_FLogsum(MMX(i,k), IVX3(ivx_3,k) + MSC_FS3(k));
      MMX(i,k) = p7_FLogsum(MMX(i,k), IVX3(ivx_4,k) + MSC_FS4(k));

      IMX(i,k) = p7_FLogsum(MMX(i-3,k) + TSC(p7P_MI,k),
                            IMX(i-3,k) + TSC(p7P_II,k));

      DMX(i,k) = p7_FLogsum(MMX(i,k-1) + TSC(p7P_MD,k-1),
                            DMX(i,k-1) + TSC(p7P_DD,k-1));

      XMX(i,p7G_E) = p7_FLogsum(MMX(i,k) + esc,
                     p7_FLogsum(DMX(i,k) + esc,
                                XMX(i,p7G_E)));
    }

  IVX3(ivx_2,M) = p7_FLogsum(MMX(i-2,M-1)   + TSC(p7P_MM,M-1),
                  p7_FLogsum(IMX(i-2,M-1)   + TSC(p7P_IM,M-1),
                  p7_FLogsum(DMX(i-2,M-1)   + TSC(p7P_DM,M-1),
                             XMX(i-2,p7G_B) + TSC(p7P_BM,M-1))));

  MMX(i,M) =                      IVX3(ivx_2,M) + MSC_FS2(M);
  MMX(i,M) = p7_FLogsum(MMX(i,M), IVX3(ivx_3,M) + MSC_FS3(M));
  MMX(i,M) = p7_FLogsum(MMX(i,M), IVX3(ivx_4,M) + MSC_FS4(M));

  IMX(i,M) = -eslINFINITY;

  DMX(i,M) = p7_FLogsum(MMX(i,M-1) + TSC(p7P_MD,M-1),
                        DMX(i,M-1) + TSC(p7P_DD,M-1));

  XMX(i,p7G_E) = p7_FLogsum(MMX(i,M),
                 p7_FLogsum(DMX(i,M),
                            XMX(i,p7G_E)));

  XMX(i,p7G_N) =            XMX(i-3,p7G_N) + gm_fs3->xsc[p7P_N][p7P_LOOP];

  XMX(i,p7G_J) = p7_FLogsum(XMX(i-3,p7G_J) + gm_fs3->xsc[p7P_J][p7P_LOOP],
                            XMX(i,p7G_E)   + gm_fs3->xsc[p7P_E][p7P_LOOP]);

  XMX(i,p7G_C) = p7_FLogsum(XMX(i-3,p7G_C) + gm_fs3->xsc[p7P_C][p7P_LOOP],
                            XMX(i,p7G_E)   + gm_fs3->xsc[p7P_E][p7P_MOVE]);

  XMX(i,p7G_B) = p7_FLogsum(XMX(i,p7G_N) + gm_fs3->xsc[p7P_N][p7P_MOVE],
                            XMX(i,p7G_J) + gm_fs3->xsc[p7P_J][p7P_MOVE]);
}


/* Function:  p7_GForward_Frameshift_3Codons()
 * Synopsis:  The frameshift Forward algorithm, 3 codon lengths.
 *
 * Purpose:   The Forward dynamic programming algorithm for frameshift
 *            aware translated comparison of the DNA sequence <dsq> of
 *            length <L> to the profile <gm_fs3>, with three codon
 *            lengths (2, 3 and 4 nucleotides).
 *
 *            Caller provides the DP matrix <gx>, allocated for at least
 *            <gm_fs3->M> by <L> with <p7G_NSCELLS> cells per position
 *            (M, I, D), and the intermediate value matrix <iv> with
 *            <p7P_3CODONS> rows. Each row is computed by
 *            <p7_GForwardRow_Frameshift_3Codons()>. The Forward matrix
 *            is returned in <gx> and the Forward score in <opt_sc>.
 *
 *            The Forward score is in lod score form. To convert to a
 *            bitscore, the caller needs to subtract a null model lod
 *            score, then convert to bits.
 *
 *            Caller must have initialized the log-sum calculation
 *            with a call to <p7_FLogsumInit()>.
 *
 * Args:      dsq    - nucleotide sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gm_fs3 - frameshift aware codon profile, 3 codon lengths
 *            gx     - DP matrix with room for an MxL alignment
 *            iv     - intermediate value matrix
 *            opt_sc - optRETURN: Forward lod score in nats
 *
 * Return:    <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <gm_fs3> isn't for 3 codon lengths or <gx>
 *            doesn't have <p7G_NSCELLS> cells per position.
 */
int
p7_GForward_Frameshift_3Codons(const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs3, P7_GMX *gx, P7_IVX *iv, float *opt_sc)
{
  float *xmx = gx->xmx;
  int    i;

  if (gm_fs3->codon_lengths != 3) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 3 codon lengths");
  if (gx->nscells != p7G_NSCELLS)    ESL_EXCEPTION(eslEINVAL, "DP matrix not allocated with p7G_NSCELLS cells");

  for (i = 0; i <= L; i++)
    p7_GForwardRow_Frameshift_3Codons(dsq, gm_fs3, gx->dp, xmx, iv->ivx, i);

  if (opt_sc != NULL) *opt_sc = p7_FLogsum( XMX(L,p7G_C),
                                p7_FLogsum( XMX(L-1,p7G_C) + gm_fs3->xsc[p7P_C][p7P_LOOP],
                                            (L >= 2 ? XMX(L-2,p7G_C) + gm_fs3->xsc[p7P_C][p7P_LOOP] : -eslINFINITY))) +
                                            gm_fs3->xsc[p7P_C][p7P_MOVE];
  gx->M = gm_fs3->M;
  gx->L = L;
  return eslOK;
}


/*****************************************************************
 * 2. Backward, 5 and 3 codon lengths.
 *****************************************************************/

/* nuc_at_fs5()
 * The nucleotide at position <j> of <dsq> for codon lookup: the
 * placeholder <p7P_MAXCODONS5> past <L> or for anything but A,C,G,T.
 */
static inline int
nuc_at_fs5(const ESL_DSQ *dsq, int L, int j)
{
  return (j <= L && dsq[j] < p7P_MAXNUC) ? dsq[j] : p7P_MAXCODONS5;
}

/* Function:  p7_GBackwardRow_Frameshift_5Codons()
 * Synopsis:  One row of the frameshift Backward algorithm.
 *
 * Purpose:   Frameshift Backward row <i> (0..L) of the 5-codon profile
 *            <gm_fs5> against <dsq> of length <L>: core cells through
 *            row pointers <dp>, special states in <xmx>, with <ivx>
 *            (M+1 floats) as workspace. Reads rows i+1..i+5 (core) and
 *            the special states of row i+3; row <L> is the
 *            initialization. The codon context of row <i> is read from
 *            dsq[i+1..i+5]. Same calculation as each row of
 *            <p7_GBackward_Frameshift_5Codons()>, which calls it.
 *
 *            Caller must have initialized the log-sum calculation with
 *            a call to <p7_FLogsumInit()>.
 */
void
p7_GBackwardRow_Frameshift_5Codons(const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs5, float **dp, float *xmx, float *ivx, int i)
{
  float const *tsc = gm_fs5->tsc;
  float const *rsc1, *rsc2, *rsc3, *rsc4, *rsc5;
  int          M   = gm_fs5->M;
  float        esc = p7_fs_profile_IsLocal(gm_fs5) ? 0 : -eslINFINITY;
  int          k;
  int          c1, c2, c3, c4, c5;
  int          t, u, v, w, x;

  x = nuc_at_fs5(dsq, L, i+1);
  w = nuc_at_fs5(dsq, L, i+2);
  v = nuc_at_fs5(dsq, L, i+3);
  u = nuc_at_fs5(dsq, L, i+4);
  t = nuc_at_fs5(dsq, L, i+5);

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
      IMX(L-1,M)  = -eslINFINITY;                       /* no I_M state        */
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
      IMX(L-2,M)  = -eslINFINITY;                       /* no I_M state        */
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


/* Function:  p7_GBackward_Frameshift_5Codons()
 * Synopsis:  The frameshift Backward algorithm, 5 codon lengths.
 *
 * Purpose:   The Backward dynamic programming algorithm for frameshift
 *            aware translated comparison of the DNA sequence <dsq> of
 *            length <L> to the profile <gm_fs5>, with five codon
 *            lengths (1 to 5 nucleotides).
 *
 *            Caller provides the DP matrix <gx>, allocated for at least
 *            <gm_fs5->M> by <L> with <p7G_NSCELLS> cells per position,
 *            and the intermediate value matrix <iv> (at least M+1
 *            floats). Each row is computed by
 *            <p7_GBackwardRow_Frameshift_5Codons()>. The Backward
 *            matrix is returned in <gx> and the Backward score in
 *            <opt_sc>.
 *
 *            The Backward score is in lod score form. To convert to a
 *            bitscore, the caller needs to subtract a null model lod
 *            score, then convert to bits.
 *
 *            Caller must have initialized the log-sum calculation
 *            with a call to <p7_FLogsumInit()>.
 *
 * Args:      dsq    - nucleotide sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gm_fs5 - frameshift aware codon profile, 5 codon lengths
 *            gx     - DP matrix with room for an MxL alignment
 *            iv     - intermediate value matrix
 *            opt_sc - optRETURN: Backward lod score in nats
 *
 * Return:    <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <gm_fs5> isn't for 5 codon lengths or <gx>
 *            doesn't have <p7G_NSCELLS> cells per position.
 */
int
p7_GBackward_Frameshift_5Codons(const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs5, P7_GMX *gx, P7_IVX *iv, float *opt_sc)
{
  float *xmx = gx->xmx;
  int    i;

  if (gm_fs5->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 5 codon lengths");
  if (gx->nscells != p7G_NSCELLS)    ESL_EXCEPTION(eslEINVAL, "DP matrix not allocated with p7G_NSCELLS cells");

  for (i = L; i >= 0; i--)
    p7_GBackwardRow_Frameshift_5Codons(dsq, L, gm_fs5, gx->dp, gx->xmx, iv->ivx, i);

  if (opt_sc != NULL) *opt_sc =  p7_FLogsum( XMX(0,p7G_N),
                                 p7_FLogsum( XMX(1,p7G_N),
                                             (L >= 2 ? XMX(2,p7G_N) : -eslINFINITY)));
  gx->M = gm_fs5->M;
  gx->L = L;
  return eslOK;
}


/* nuc_at_fs3()
 * The nucleotide at position <j> of <dsq> for 3-codon lookup: the
 * placeholder <p7P_MAXCODONS3> past <L> or for anything but A,C,G,T.
 */
static inline int
nuc_at_fs3(const ESL_DSQ *dsq, int L, int j)
{
  return (j <= L && dsq[j] < p7P_MAXNUC) ? dsq[j] : p7P_MAXCODONS3;
}

/* Function:  p7_GBackwardRow_Frameshift_3Codons()
 * Synopsis:  One row of the frameshift Backward algorithm, 3 codon lengths.
 *
 * Purpose:   Frameshift Backward row <i> (0..L) of the 3-codon profile
 *            <gm_fs3> against <dsq> of length <L>: core cells
 *            (<p7G_NSCELLS> per position) through row pointers <dp>,
 *            special states in <xmx>, with <ivx> (M+1 floats) as
 *            workspace. Reads rows i+2..i+4 (core) and the special
 *            states of row i+3; row <L> is the initialization. The
 *            codon context of row <i> is read from dsq[i+1..i+4]. Same
 *            calculation as each row of
 *            <p7_GBackward_Frameshift_3Codons()>, which calls it.
 *
 *            Caller must have initialized the log-sum calculation with
 *            a call to <p7_FLogsumInit()>.
 */
void
p7_GBackwardRow_Frameshift_3Codons(const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs3, float **dp, float *xmx, float *ivx, int i)
{
  float const *tsc = gm_fs3->tsc;
  float const *rsc2, *rsc3, *rsc4;
  int          M   = gm_fs3->M;
  float        esc = p7_fs_profile_IsLocal(gm_fs3) ? 0 : -eslINFINITY;
  int          k;
  int          c2, c3, c4;
  int          u, v, w, x;

  if (i >= L-1 && i >= 1)
    {
      XMX(i,p7G_C) = (i == L ? gm_fs3->xsc[p7P_C][p7P_MOVE] : gm_fs3->xsc[p7P_C][p7P_LOOP] + gm_fs3->xsc[p7P_C][p7P_MOVE]);
      XMX(i,p7G_J) = XMX(i,p7G_B) = XMX(i,p7G_N) = -eslINFINITY;  /* need to enter and exit model */
      XMX(i,p7G_E) = XMX(i,p7G_C) + gm_fs3->xsc[p7P_E][p7P_MOVE]; /* E<-C, no tail */
      MMX(i,M)     = DMX(i,M) = XMX(i,p7G_E);                     /* {MD}_M <- E (prob 1.0) */
      IMX(i,M)     = -eslINFINITY;                                /* no I_M state        */

      for (k = M-1; k >= 1; k--)
        {
          MMX(i,k) = p7_FLogsum( XMX(i,p7G_E) + esc,
                                 DMX(i,k+1)   + TSC(p7P_MD,k));
          DMX(i,k) = p7_FLogsum( XMX(i,p7G_E) + esc,
                                 DMX(i,k+1)   + TSC(p7P_DD,k));
          IMX(i,k) = -eslINFINITY;
        }
      MMX(i,0) = IMX(i,0) = DMX(i,0) = -eslINFINITY;
      return;
    }

  x = nuc_at_fs3(dsq, L, i+1);
  w = nuc_at_fs3(dsq, L, i+2);
  v = nuc_at_fs3(dsq, L, i+3);
  u = nuc_at_fs3(dsq, L, i+4);

  c2 = p7P_CODON2_FS3(x, w);
  c2 = p7P_MINIDX(c2, p7P_DEGEN3_QC1);
  rsc2 = gm_fs3->rsc[c2];

  if (i == 0)
    {
      /* At i=0, only N,B states are reachable. Rows past L do not
       * exist; when L < 4 they contribute nothing. */
      c3 = p7P_CODON3_FS3(x, w, v);
      c3 = p7P_MINIDX(c3, p7P_DEGEN3_C);
      c4 = p7P_CODON4_FS3(x, w, v, u);
      c4 = p7P_MINIDX(c4, p7P_DEGEN3_QC1);
      rsc3 = gm_fs3->rsc[c3];
      rsc4 = gm_fs3->rsc[c4];

      ivx[1] = p7_FLogsum( (L >= 2 ? MMX(2,1) + MSC_FS2(1) : -eslINFINITY),
               p7_FLogsum( (L >= 3 ? MMX(3,1) + MSC_FS3(1) : -eslINFINITY),
                           (L >= 4 ? MMX(4,1) + MSC_FS4(1) : -eslINFINITY)));

      XMX(0,p7G_B) = ivx[1] + TSC(p7P_BM,0);

      for (k = 2; k <= M; k++)
        {
          ivx[k] = p7_FLogsum( (L >= 2 ? MMX(2,k) + MSC_FS2(k) : -eslINFINITY),
                   p7_FLogsum( (L >= 3 ? MMX(3,k) + MSC_FS3(k) : -eslINFINITY),
                               (L >= 4 ? MMX(4,k) + MSC_FS4(k) : -eslINFINITY)));

          XMX(0,p7G_B) = p7_FLogsum(XMX(0,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
        }

      XMX(0,p7G_J) = XMX(0,p7G_C) = XMX(0,p7G_E) = -eslINFINITY;

      XMX(0,p7G_N) = p7_FLogsum( (L >= 3 ? XMX(3,p7G_N) + gm_fs3->xsc[p7P_N][p7P_LOOP] : -eslINFINITY),
                                 XMX(0,p7G_B) + gm_fs3->xsc[p7P_N][p7P_MOVE]);

      for (k = M; k >= 0; k--)
        MMX(0,k) = DMX(0,k) = IMX(0,k) = -eslINFINITY;
      return;
    }

  if (i == L-2)
    {
      ivx[1]       = MMX(L,1) + MSC_FS2(1);
      XMX(i,p7G_B) = ivx[1]   + TSC(p7P_BM,0);
      for (k = 2; k <= M; k++)
        {
          ivx[k] = MMX(L,k) + MSC_FS2(k);
          XMX(i,p7G_B) = p7_FLogsum( XMX(i,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
        }

      XMX(i,p7G_J) = XMX(i,p7G_B) + gm_fs3->xsc[p7P_J][p7P_MOVE];
      XMX(i,p7G_N) = XMX(i,p7G_B) + gm_fs3->xsc[p7P_N][p7P_MOVE];
      XMX(i,p7G_C) = gm_fs3->xsc[p7P_C][p7P_LOOP] + gm_fs3->xsc[p7P_C][p7P_MOVE];

      XMX(i,p7G_E) = p7_FLogsum(XMX(i,p7G_J) + gm_fs3->xsc[p7P_E][p7P_LOOP],
                                XMX(i,p7G_C) + gm_fs3->xsc[p7P_E][p7P_MOVE]);

      MMX(i,M) = DMX(i,M) = XMX(i,p7G_E);
      IMX(i,M) = -eslINFINITY;

      for (k = M-1; k >= 1; k--)
        {
          MMX(i,k) = p7_FLogsum( DMX(i,k+1)   + TSC(p7P_MD,k),
                     p7_FLogsum( ivx[k+1]     + TSC(p7P_MM,k),
                                 XMX(i,p7G_E) + esc));

          DMX(i,k) = p7_FLogsum( p7_FLogsum( XMX(i,p7G_E) + esc,
                                             DMX(i,k+1)   + TSC(p7P_DD,k)),
                                             ivx[k+1]     + TSC(p7P_DM,k));

          IMX(i,k) = ivx[k+1] + TSC(p7P_IM,k);
        }

      MMX(i,0) = IMX(i,0) = DMX(i,0) = -eslINFINITY;
      return;
    }

  /* find correct index for looking up scores of codons and quasicodons */
  c3 = p7P_CODON3_FS3(x, w, v);
  c3 = p7P_MINIDX(c3, p7P_DEGEN3_C);
  c4 = p7P_CODON4_FS3(x, w, v, u);
  c4 = p7P_MINIDX(c4, p7P_DEGEN3_QC1);
  rsc3 = gm_fs3->rsc[c3];
  rsc4 = gm_fs3->rsc[c4];

  if (i > L-5)   /* rows L-3 and L-4: row i+4 exists only for L-4 */
    {
      ivx[1] =                     MMX(i+2,1) + MSC_FS2(1);
      ivx[1] = p7_FLogsum( ivx[1], MMX(i+3,1) + MSC_FS3(1));
      if (i == L-4)
        ivx[1] = p7_FLogsum( ivx[1], MMX(i+4,1) + MSC_FS4(1));

      XMX(i,p7G_B) = ivx[1] + TSC(p7P_BM,0);

      for (k = 2; k <= M; k++)
        {
          ivx[k] =                     MMX(i+2,k) + MSC_FS2(k);
          ivx[k] = p7_FLogsum( ivx[k], MMX(i+3,k) + MSC_FS3(k));
          if (i == L-4)
            ivx[k] = p7_FLogsum( ivx[k], MMX(i+4,k) + MSC_FS4(k));

          XMX(i,p7G_B) = p7_FLogsum( XMX(i,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
        }
    }
  else
    {
      ivx[1] = p7_FLogsum( MMX(i+2,1) + MSC_FS2(1),
               p7_FLogsum( MMX(i+3,1) + MSC_FS3(1),
                           MMX(i+4,1) + MSC_FS4(1)));

      XMX(i,p7G_B) = ivx[1] + TSC(p7P_BM,0);

      for (k = 2; k <= M; k++)
        {
          ivx[k] = p7_FLogsum( MMX(i+2,k) + MSC_FS2(k),
                   p7_FLogsum( MMX(i+3,k) + MSC_FS3(k),
                               MMX(i+4,k) + MSC_FS4(k)));

          XMX(i,p7G_B) = p7_FLogsum( XMX(i,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
        }
    }

  XMX(i,p7G_J) = p7_FLogsum( XMX(i+3,p7G_J) + gm_fs3->xsc[p7P_J][p7P_LOOP],
                             XMX(i,  p7G_B) + gm_fs3->xsc[p7P_J][p7P_MOVE]);

  XMX(i,p7G_C) =             XMX(i+3,p7G_C) + gm_fs3->xsc[p7P_C][p7P_LOOP];

  XMX(i,p7G_N) = p7_FLogsum( XMX(i+3,p7G_N) + gm_fs3->xsc[p7P_N][p7P_LOOP],
                             XMX(i,  p7G_B) + gm_fs3->xsc[p7P_N][p7P_MOVE]);

  XMX(i,p7G_E) = p7_FLogsum( XMX(i,p7G_J) + gm_fs3->xsc[p7P_E][p7P_LOOP],
                             XMX(i,p7G_C) + gm_fs3->xsc[p7P_E][p7P_MOVE]);

  MMX(i,M) = DMX(i,M) = XMX(i,p7G_E); /* {MD}_M <- E (prob 1.0) */
  IMX(i,M) = -eslINFINITY;            /* no I_M state        */

  if (i > L-5)
    for (k = M-1; k >= 1; k--)
      {
        MMX(i,k) = p7_FLogsum( DMX(i,k+1)   + TSC(p7P_MD,k),
                   p7_FLogsum( IMX(i+3,k)   + TSC(p7P_MI,k),
                   p7_FLogsum( ivx[k+1]     + TSC(p7P_MM,k),
                               XMX(i,p7G_E) + esc)));

        DMX(i,k) = p7_FLogsum( DMX(i,k+1)   + TSC(p7P_DD,k),
                   p7_FLogsum( XMX(i,p7G_E) + esc,
                               ivx[k+1]     + TSC(p7P_DM,k)));

        IMX(i,k) = p7_FLogsum( IMX(i+3,k)   + TSC(p7P_II,k),
                               ivx[k+1]     + TSC(p7P_IM,k));
      }
  else
    for (k = M-1; k >= 1; k--)
      {
        MMX(i,k) = p7_FLogsum( p7_FLogsum( DMX(i,k+1)   + TSC(p7P_MD,k),
                               p7_FLogsum( IMX(i+3,k)   + TSC(p7P_MI,k),
                                           ivx[k+1]     + TSC(p7P_MM,k))),
                                           XMX(i,p7G_E) + esc);

        DMX(i,k) = p7_FLogsum( p7_FLogsum( XMX(i,p7G_E) + esc,
                                           DMX(i,k+1)   + TSC(p7P_DD,k)),
                                           ivx[k+1]     + TSC(p7P_DM,k));

        IMX(i,k) = p7_FLogsum(             IMX(i+3,k)   + TSC(p7P_II,k),
                                           ivx[k+1]     + TSC(p7P_IM,k));
      }

  MMX(i,0) = IMX(i,0) = DMX(i,0) = -eslINFINITY;
}


/* Function:  p7_GBackward_Frameshift_3Codons()
 * Synopsis:  The frameshift Backward algorithm, 3 codon lengths.
 *
 * Purpose:   The Backward dynamic programming algorithm for frameshift
 *            aware translated comparison of the DNA sequence <dsq> of
 *            length <L> to the profile <gm_fs3>, with three codon
 *            lengths (2, 3 and 4 nucleotides).
 *
 *            Caller provides the DP matrix <gx>, allocated for at least
 *            <gm_fs3->M> by <L> with <p7G_NSCELLS> cells per position,
 *            and the intermediate value matrix <iv> (at least M+1
 *            floats). Each row is computed by
 *            <p7_GBackwardRow_Frameshift_3Codons()>. The Backward
 *            matrix is returned in <gx> and the Backward score in
 *            <opt_sc>.
 *
 *            The Backward score is in lod score form. To convert to a
 *            bitscore, the caller needs to subtract a null model lod
 *            score, then convert to bits.
 *
 *            Caller must have initialized the log-sum calculation
 *            with a call to <p7_FLogsumInit()>.
 *
 * Args:      dsq    - nucleotide sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gm_fs3 - frameshift aware codon profile, 3 codon lengths
 *            gx     - DP matrix with room for an MxL alignment
 *            iv     - intermediate value matrix
 *            opt_sc - optRETURN: Backward lod score in nats
 *
 * Return:    <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <gm_fs3> isn't for 3 codon lengths or <gx>
 *            doesn't have <p7G_NSCELLS> cells per position.
 */
int
p7_GBackward_Frameshift_3Codons(const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs3, P7_GMX *gx, P7_IVX *iv, float *opt_sc)
{
  float *xmx = gx->xmx;
  int    i;

  if (gm_fs3->codon_lengths != 3) ESL_EXCEPTION(eslEINVAL, "profile not allocated for 3 codon lengths");
  if (gx->nscells != p7G_NSCELLS)    ESL_EXCEPTION(eslEINVAL, "DP matrix not allocated with p7G_NSCELLS cells");

  for (i = L; i >= 0; i--)
    p7_GBackwardRow_Frameshift_3Codons(dsq, L, gm_fs3, gx->dp, xmx, iv->ivx, i);

  if (opt_sc != NULL) *opt_sc =  p7_FLogsum( XMX(0,p7G_N),
                                 p7_FLogsum( XMX(1,p7G_N),
                                             (L >= 2 ? XMX(2,p7G_N) : -eslINFINITY)));
  gx->M = gm_fs3->M;
  gx->L = L;
  return eslOK;
}



/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/
#ifdef p7GENERIC_FWDBACK_FRAMESHIFT_BENCHMARK
/*
   gcc -g -O3      -o generic_fwdback_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_FWDBACK_FRAMESHIFT_BENCHMARK generic_fwdback_frameshift.c -lhmmer -leasel -lm
   icc -O3 -static -o generic_fwdback_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_FWDBACK_FRAMESHIFT_BENCHMARK generic_fwdback_frameshift.c -lhmmer -leasel -lm
   ./generic_fwdback_frameshift_benchmark <hmmfile>
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-L",        eslARG_INT,   "1200", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs",                   0 },
  { "-N",        eslARG_INT,   "2000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                   0 },
  { "-B",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark Backward_Frameshift()",           0 },
  { "-F",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark Forward_Frameshift()",            0 },
  { "-T",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-V", "only benchmark 3 codon lengths",                 0 },
  { "-V",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, "-T", "only benchmark 5 codon lengths",                 0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for generic Forward/Backward";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA   = NULL;
  ESL_ALPHABET   *abcDNA  = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bgDNA   = NULL;
  P7_BG          *bgAA    = NULL;
  P7_FS_PROFILE  *gm_fs5  = NULL;
  P7_FS_PROFILE  *gm_fs3  = NULL;
  P7_GMX         *fwd     = NULL;
  P7_GMX         *fwd3    = NULL;
  P7_GMX         *bck     = NULL;
  P7_IVX         *iv      = NULL;
  ESL_GENCODE    *gcode   = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           sc;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)            != eslOK) p7_Fail("Failed to read HMM");

  abcDNA = esl_alphabet_Create(eslDNA);
  gcode = esl_gencode_Create(abcDNA, abcAA);
  bgDNA = p7_bg_Create(abcDNA);
  bgAA  = p7_bg_Create(abcAA);
  gm_fs5 = p7_profile_fs_Create(hmm->M, abcAA, p7P_5CODONS);
  gm_fs3 = p7_profile_fs_Create(hmm->M, abcAA, p7P_3CODONS);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_UNILOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs3, L/3, p7_UNILOCAL);

  fwd   = p7_gmx_Create(gm_fs5->M, L, L, p7G_NSCELLS_FS);
  fwd3  = p7_gmx_Create(gm_fs3->M, L, L, p7G_NSCELLS);
  bck   = p7_gmx_Create(gm_fs5->M, L, L, p7G_NSCELLS);
  iv    = p7_ivx_Create(gm_fs5->M, p7P_5CODONS);

  /* Baseline time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++) esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark time. */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
  {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
      if (! esl_opt_GetBoolean(go, "-B"))  {
        if (! esl_opt_GetBoolean(go, "-T")) p7_GForward_Frameshift_5Codons(dsq, L, gm_fs5, fwd,  iv, &sc);
        if (! esl_opt_GetBoolean(go, "-V")) p7_GForward_Frameshift_3Codons(dsq, L, gm_fs3, fwd3, iv, &sc);
      }
      if (! esl_opt_GetBoolean(go, "-F"))  {
        if (! esl_opt_GetBoolean(go, "-T")) p7_GBackward_Frameshift_5Codons(dsq, L, gm_fs5, bck, iv, NULL);
        if (! esl_opt_GetBoolean(go, "-V")) p7_GBackward_Frameshift_3Codons(dsq, L, gm_fs3, bck, iv, NULL);
     }

     p7_gmx_Reuse(fwd);
     p7_gmx_Reuse(fwd3);
     p7_gmx_Reuse(bck);
    
  }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm_fs5->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm_fs5->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_gmx_Destroy(bck);
  p7_gmx_Destroy(fwd3);
  p7_gmx_Destroy(fwd);
  p7_ivx_Destroy(iv);
  p7_profile_fs_Destroy(gm_fs5);
  p7_profile_fs_Destroy(gm_fs3);
  p7_bg_Destroy(bgDNA);
  p7_bg_Destroy(bgAA);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcAA);
  esl_alphabet_Destroy(abcDNA);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_FWDBACK_FRAMESHIFT_BENCHMARK*/
/*----------------- end, benchmark ------------------------------*/




/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7GENERIC_FWDBACK_FRAMESHIFT_TESTDRIVE
#include <string.h>
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_vectorops.h"

static void
utest_forward_fs(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, P7_CODONTABLE *codon_table, P7_BG *bgAA, P7_HMM *hmm, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs3, P7_FS_PROFILE *gm_fs5, int nseq, int L)
{
  int         i,j;
  float       avg_sc_rnd;
  float       avg_sc_gen;
  ESL_SQ      *sq;
  P7_TRACE    *tr;
  ESL_DSQ     *dsqAA  = NULL;
  ESL_DSQ     *dsqDNA = NULL;
  P7_GMX      *fwd    = NULL;
  P7_GMX      *fwd3   = NULL;
  P7_GMX      *bck    = NULL;
  P7_GMX      *vit    = NULL;
  P7_IVX      *iv     = NULL;
  int         idx;
  float       vsc, fsc, bsc;
  float       fsc3, bsc3;
  float       nullsc;

  p7_FLogsumInit();
 
  if ((sq     = esl_sq_CreateDigital(abcAA))                                  == NULL)  esl_fatal("sequence creation failed");
  if ((dsqAA  = malloc(sizeof(ESL_DSQ) *((L/3)+2)))                           == NULL)  esl_fatal("malloc failed");
  if ((dsqDNA = malloc(sizeof(ESL_DSQ) *(L+2)))                               == NULL)  esl_fatal("malloc failed");
  if ((fwd    = p7_gmx_Create(gm_fs5->M, L,               L, p7G_NSCELLS_FS)) == NULL)  esl_fatal("matrix creation failed");
  if ((fwd3   = p7_gmx_Create(gm_fs3->M, L,               L, p7G_NSCELLS))    == NULL)  esl_fatal("matrix creation failed");
  if ((bck    = p7_gmx_Create(gm_fs5->M, L,               L, p7G_NSCELLS))    == NULL)  esl_fatal("matrix creation failed");
  if ((vit    = p7_gmx_Create(gm_fs5->M, L,               L, p7G_NSCELLS))    == NULL)  esl_fatal("matrix creation failed");
  if ((tr     = p7_trace_Create())                                            == NULL)  esl_fatal("trace creation failed");
  if ((iv     = p7_ivx_Create(gm_fs5->M, p7P_5CODONS))                        == NULL)  esl_fatal("ivx creation failed");

  /* Compare Viterbi, Forward, and Backward scores when aligning to DNA sequences reverse translated 
     from randomly generated Amino Acid sequences. Keep track of the average Forward score */
  avg_sc_rnd = 0.;
  for (idx = 0; idx < nseq; idx++)
    {
      if (esl_rsq_xfIID(r, bgAA->f, abcAA->K, (L/3), dsqAA) != eslOK) esl_fatal("seq generation failed");
      j = 1; 
	  for(i = 1; i <= (L/3); i++) {
        p7_codontable_GetCodon(codon_table, r, dsqAA[i], dsqDNA+j);
        j+=3;
	  }
      if (p7_GViterbi_Frameshift(dsqDNA, L, gm_fs5, vit, iv, &vsc)      != eslOK) esl_fatal("viterbi failed");
      if (p7_GForward_Frameshift_5Codons(dsqDNA, L, gm_fs5, fwd, iv, &fsc)      != eslOK) esl_fatal("forward failed");
      if (p7_GBackward_Frameshift_5Codons(dsqDNA, L, gm_fs5, bck, iv, &bsc)     != eslOK) esl_fatal("backward failed");
     
      if (fsc < vsc)             esl_fatal("Forward score can't be less than Viterbi score"); 
      if (fabs(fsc-bsc) > 0.001) esl_fatal("Forward/Backward failed: %f %f\n", fsc, bsc);

      if (p7_bg_fs_NullOne(bgAA, dsqAA, (L/3), &nullsc)      != eslOK) esl_fatal("null score failed");
   
      avg_sc_rnd += fsc - nullsc;

      p7_gmx_Reuse(bck);
      if (p7_GForward_Frameshift_3Codons (dsqDNA, L, gm_fs3, fwd3, iv, &fsc3) != eslOK) esl_fatal("forward 3 failed");
      if (p7_GBackward_Frameshift_3Codons(dsqDNA, L, gm_fs3, bck,  iv, &bsc3) != eslOK) esl_fatal("backward 3 failed");

      if (fabs(fsc3-bsc3) > 0.001) esl_fatal("Forward 3 /Backward 3 failed: %f %f\n", fsc3, bsc3);

      if (esl_opt_GetBoolean(go, "--vv")) 
        printf("utest_forward_fs: Forward score on random sequence len %d: %.4f (total so far: %.4f)\n", L, fsc, avg_sc_rnd);

      p7_gmx_Reuse(fwd);
      p7_gmx_Reuse(fwd3);
      p7_gmx_Reuse(bck);
      p7_gmx_Reuse(vit);
    }

  avg_sc_rnd /= (float) nseq;


  /* Get the average forward score on DNA sequence reverse translated from Amino Acid 
   * sequence generated by the model and compare to the average on random sequence */
  avg_sc_gen = 0.;
  for (idx = 0; idx < nseq; idx++)
    {
      if (p7_ProfileEmit(r, hmm, gm, bgAA, sq, tr)     != eslOK) esl_fatal("profile emission failed");
          
	  if(dsqDNA != NULL) free(dsqDNA);
      if ((dsqDNA = malloc(sizeof(ESL_DSQ) *(sq->n*3+2))) == NULL)  esl_fatal("malloc failed"); 	  

	  j = 1;
      for(i = 1; i <= sq->n; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsqDNA+j);
        j+=3;
      }
      p7_ReconfigLength(gm, sq->n);

	  p7_gmx_GrowTo(fwd, gm->M, (sq->n*3), (sq->n*3));

	  if (p7_GForward_Frameshift_5Codons(dsqDNA, (sq->n*3), gm_fs5, fwd, iv, &fsc)      != eslOK) esl_fatal("forward failed");

	  p7_bg_SetLength(bgAA,  sq->n);

      p7_bg_fs_NullOne(bgAA, dsqAA, sq->n, &nullsc);	  

	  avg_sc_gen += fsc - nullsc;

	  if (esl_opt_GetBoolean(go, "--vv"))
        printf("utest_forward_fs: Forward score on generated sequence len %d: %.4f (total so far: %.4f)\n", (int)(sq->n*3), fsc, avg_sc_gen);
    }

  avg_sc_gen /= (float) nseq;

  if (avg_sc_rnd > avg_sc_gen) esl_fatal("Average Forwrd scores on random sequence (%f) and better than on generated seqeucen (%f)", avg_sc_rnd, avg_sc_gen);

  esl_sq_Destroy(sq);
  p7_gmx_Destroy(fwd3);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(bck);
  p7_gmx_Destroy(vit);
  p7_ivx_Destroy(iv);
  p7_trace_Destroy(tr);
  free(dsqAA);
  free(dsqDNA);
  return;
}


/* utest_fwdback_agree_fs()
 * The frameshift Forward and Backward (5 codon lengths) must sum the
 * same paths: on random DNA and on the back-translated consensus of
 * <hmm> with two nucleotides deleted and two inserted, with the model's
 * frameshift probability raised to 0.05 so that codon lengths other
 * than 3 carry weight, the two scores must agree to within 0.01 nats.
 */
static void
utest_fwdback_agree_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_CODONTABLE *codon_table, P7_BG *bgAA, P7_BG *bgDNA, P7_HMM *hmm, int nseq)
{
  char          *msg     = "frameshift Forward/Backward agreement unit test failed";
  int            M       = hmm->M;
  int            L       = 3*M + 100;
  float          fsprob  = hmm->fsprob;
  P7_FS_PROFILE *gm_fs5  = p7_profile_fs_Create(M, abcAA, 5);
  P7_GMX        *fx      = p7_gmx_Create(M, L, L, p7G_NSCELLS_FS);
  P7_GMX        *bx      = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_IVX        *iv      = p7_ivx_Create(M, p7P_5CODONS);
  ESL_SQ        *sq      = esl_sq_CreateDigital(abcAA);
  ESL_DSQ       *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  ESL_DSQ       *tmp     = malloc(sizeof(ESL_DSQ) * (L+2));
  float          fsc, bsc;
  int            n, i, j, Ln, start, del1, del2, ins1, ins2;

  hmm->fsprob = 0.05;
  p7_emit_SimpleConsensus(hmm, sq);
  for (n = 0; n < nseq; n++)
    {
      Ln = 3*M + 60;
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, Ln, tmp);
      if (n % 2 == 1)
        {
          start = 1 + esl_rnd_Roll(r, Ln - 3*M + 1);
          for (i = 1, j = start; i <= M; i++, j += 3) p7_codontable_GetCodon(codon_table, r, sq->dsq[i], tmp+j);
        }
      del1 = 1 + esl_rnd_Roll(r, Ln);  del2 = 1 + esl_rnd_Roll(r, Ln);
      ins1 = 1 + esl_rnd_Roll(r, Ln);  ins2 = 1 + esl_rnd_Roll(r, Ln);
      for (i = 1, j = 1; i <= Ln; i++)
        {
          if (i == del1 || i == del2) continue;
          if (i == ins1 || i == ins2) dsq[j++] = esl_rnd_Roll(r, 4);
          dsq[j++] = tmp[i];
        }
      Ln = j - 1;

      if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, Ln/3, p7_LOCAL)          != eslOK) esl_fatal(msg);
      if (p7_GForward_Frameshift_5Codons (dsq, Ln, gm_fs5, fx, iv, &fsc) != eslOK) esl_fatal(msg);
      if (p7_GBackward_Frameshift_5Codons(dsq, Ln, gm_fs5, bx, iv, &bsc) != eslOK) esl_fatal(msg);
      if (fabs(fsc - bsc) > 0.01)
        esl_fatal("%s: Forward %.4f, Backward %.4f nats", msg, fsc, bsc);
    }

  hmm->fsprob = fsprob;
  free(dsq);
  free(tmp);
  esl_sq_Destroy(sq);
  p7_ivx_Destroy(iv);
  p7_gmx_Destroy(fx);
  p7_gmx_Destroy(bx);
  p7_profile_fs_Destroy(gm_fs5);
}
#endif /*p7GENERIC_FWDBACK_FRAMESHIFT_TESTDRIVE*/
/*------------------------- end, unit tests ---------------------*/

/*****************************************************************
 * 5. Test driver.
 *****************************************************************/

/* gcc -g -Wall -Dp7GENERIC_FWDBACK_FRAMESHIFT_TESTDRIVE -I. -I../easel -L. -L../easel -o generic_fwdback_frameshift_utest generic_fwdback_frameshift.c -lhmmer -leasel -lm
 */
#ifdef p7GENERIC_FWDBACK_FRAMESHIFT_TESTDRIVE
#include "easel.h"
#include "esl_getopts.h"
#include "esl_msa.h"

#include "p7_config.h"
#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-v",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be verbose",                                     0 },
  { "--vv",      eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "be very verbose",                                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "unit test driver for the generic Forward/Backward implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA  = NULL;
  ESL_ALPHABET   *abcDNA = NULL;
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = NULL;
  P7_FS_PROFILE  *gm_fs5 = NULL;
  P7_FS_PROFILE  *gm_fs3 = NULL;
  P7_BG          *bgAA   = NULL;
  P7_BG          *bgDNA  = NULL;
  ESL_GENCODE    *gcode  = NULL;
  P7_CODONTABLE  *ct     = NULL;
  int             M      = 100;
  int             L      = 600; // must be a multiple of 3
  int             nseq   = 20;
  char            errbuf[eslERRBUFSIZE];

  p7_FLogsumInit();

  if ((abcAA  = esl_alphabet_Create(eslAMINO))                  == NULL)  esl_fatal("failed to create alphabet");
  if ((abcDNA = esl_alphabet_Create(eslDNA))                    == NULL)  esl_fatal("failed to create alphabet");
  if ((gcode  = esl_gencode_Create(abcDNA,abcAA))               == NULL)  esl_fatal("failed to create gencode");
  if ((ct     = p7_codontable_Create(gcode))                    == NULL)  esl_fatal("failed to create codon table");
  if (p7_hmm_Sample(r, M, abcAA, &hmm)                          != eslOK) esl_fatal("failed to sample an HMM");
  if ((bgAA = p7_bg_Create(abcAA))                              == NULL)  esl_fatal("failed to create null model");
  if ((bgDNA = p7_bg_Create(abcDNA))                            == NULL)  esl_fatal("failed to create null model");
  if ((p7_bg_SetLength(bgAA,  L/3))                             != eslOK) esl_fatal("failed to config background");
  if ((p7_bg_SetLength(bgDNA,  L))                              != eslOK) esl_fatal("failed to config background");
  if ((gm = p7_profile_Create(hmm->M, abcAA))                   == NULL)  esl_fatal("failed to create profile");
  if ((p7_ProfileConfig(hmm, bgAA, gm, L, p7_LOCAL))            != eslOK) esl_fatal("failed to config profile");
  if ((gm_fs5 = p7_profile_fs_Create(hmm->M, abcAA, 5))           == NULL)  esl_fatal("failed to create profile");
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_LOCAL) != eslOK) esl_fatal("failed to config profile");
  if ((gm_fs3 = p7_profile_fs_Create(hmm->M, abcAA, 3))           == NULL)  esl_fatal("failed to create profile");
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs3, L/3, p7_LOCAL) != eslOK) esl_fatal("failed to config profile");
  if (p7_hmm_Validate    (hmm, errbuf, 0.0001)      != eslOK) esl_fatal("whoops, HMM is bad!: %s", errbuf);

  utest_forward_fs    (go, r, abcAA, ct, bgAA, hmm, gm, gm_fs3, gm_fs5, nseq, L);
  utest_fwdback_agree_fs(r, abcAA, abcDNA, gcode, ct, bgAA, bgDNA, hmm, nseq);

  p7_profile_fs_Destroy(gm_fs5);
  p7_profile_fs_Destroy(gm_fs3);
  p7_profile_Destroy(gm);
  p7_bg_Destroy(bgAA);
  p7_bg_Destroy(bgDNA);
  p7_hmm_Destroy(hmm);
  p7_codontable_Destroy(ct);
  esl_alphabet_Destroy(abcAA);
  esl_alphabet_Destroy(abcDNA);
  esl_gencode_Destroy(gcode);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_FWDBACK_FRAMESHIFT_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/



