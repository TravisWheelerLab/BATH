/* Frameshift aware Viterbi algorithm; generic (non-SIMD) version.
 *
 * Full mattix, log space, Frameshift aware translated Viterbi 
 * algorithm with 5 codon lengths, plus trace back. 
 *
 * Contents:
 *    1. Frameshift Viterbi and trace back.
 *    2. Benchmark driver.
 *    3. Unit tests.
 *    4. Test driver.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#define IVX(i,k) (ivx[((k)*p7P_5CODONS) + (i)])

/*****************************************************************
 * 1. Frameshift Viterbi and trace back.
 *****************************************************************/

/* Function:  p7_GViterbi_Frameshift()
 * Synopsis:  The Viterbi algorithm, with frameshift awareness.
 * 
 * Purpose:   The Viterbi dynamic programming algorithm for aligning 
 *            DNA to proteins with frameshift awareness. 
 *
 *            Given a digital nucleotide sequence <dsq> of length <L>, 
 *            a frameshift aware codon profile <gm_fs5>, and DP matrix 
 *            <gx> allocated for at least <L> by <gm_fs5->M> cells; 
 *            calculate the maximum scoring path by fremashift-aware 
 *            Viterbi; return the Viterbi score in <ret_sc>, and the
 *            Viterbi matrix is in <gx>.
 *            
 *            The caller may then retrieve the Viterbi path by calling
 *            <p7_fs_VTrace()>.
 *           
 *            The Viterbi lod score is returned in nats. The caller
 *            needs to subtract a null model lod score, then convert
 *            to bits.
 *           
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gm_fs5  - profile. 
 *            gx     - DP matrix with room for an MxL alignment
 *            iv     - intermediate value matrix
 *            opt_sc - optRETURN: Viterbi lod score in nats
 *           
 * Return:   <eslOK> on success.
 */
int
p7_GViterbi_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs5, P7_GMX *gx, P7_IVX *iv, float *opt_sc)
{
  float const *tsc  = gm_fs5->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float       *ivx  = iv->ivx;
  int          M    = gm_fs5->M;
  int          i,k;
  int          c1, c2, c3, c4, c5;
  int          t, u, v, w, x;
  int          ivx_1, ivx_2, ivx_3, ivx_4, ivx_5;
  float        esc  = p7_fs_profile_IsLocal(gm_fs5) ? 0 : -eslINFINITY;
  int status;

  if(gm_fs5->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "proflie not allocated for 5 codon lengths");
  
  /* Initialization of the zero row.  */
  XMX_FS(0,p7G_N) = 0.;                                  /* S->N, p=1            */
  XMX_FS(0,p7G_B) = gm_fs5->xsc[p7P_N][p7P_MOVE];         /* S->N->B, no N-tail   */
  XMX_FS(0,p7G_E) = XMX_FS(0,p7G_J) = XMX_FS(0,p7G_C) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_FS(0,k,p7G_C0) = MMX_FS(0,k,p7G_C1) = MMX_FS(0,k,p7G_C2) = MMX_FS(0,k,p7G_C3) =
    MMX_FS(0,k,p7G_C4) = MMX_FS(0,k,p7G_C5) = IMX_FS(0,k)        = DMX_FS(0,k)        = -eslINFINITY;

   /* Initialization for row 1 */
  XMX_FS(1,p7G_N) = 0.;
  XMX_FS(1,p7G_B) = gm_fs5->xsc[p7P_N][p7P_MOVE];
  XMX_FS(1,p7G_E) = -eslINFINITY;
  MMX_FS(1,0,p7G_C1) = MMX_FS(1,0,p7G_C2) = MMX_FS(1,0,p7G_C3) = MMX_FS(1,0,p7G_C4) = MMX_FS(1,0,p7G_C5) = -eslINFINITY;
  MMX_FS(1,0,p7G_C0) = IMX_FS(1,0) = DMX_FS(1,0) = -eslINFINITY;

  if(dsq[1] < p7P_MAXNUC) x = dsq[1];
  else                                            x = p7P_MAXCODONS5;

  c1 = p7P_CODON1_FS5(x);
  c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);
  for (k = 1; k <= M; k++) {
    IVX(1,k) = XMX_FS(0,p7G_B) + TSC(p7P_BM,k-1);

    MMX_FS(1,k,p7G_C1) = IVX(1,k) + p7P_MSC_CODON(gm_fs5, k, c1);
    MMX_FS(1,k,p7G_C2) = -eslINFINITY;
    MMX_FS(1,k,p7G_C3) = -eslINFINITY;
    MMX_FS(1,k,p7G_C4) = -eslINFINITY;
    MMX_FS(1,k,p7G_C5) = -eslINFINITY;
    MMX_FS(1,k,p7G_C0) = MMX_FS(1,k,p7G_C1);
    IMX_FS(1,k)  = -eslINFINITY;
    DMX_FS(1,k)  = ESL_MAX(MMX_FS(1,k-1,p7G_C0) + TSC(p7P_MD,k-1),
                           DMX_FS(1,k-1)        + TSC(p7P_DD,k-1));

    XMX_FS(1,p7G_E) = ESL_MAX(MMX_FS(1,k,p7G_C0) + esc,
                      ESL_MAX(DMX_FS(1,k)        + esc,
                              XMX_FS(1,p7G_E)));
  }

  XMX_FS(1,p7G_J) = XMX_FS(1,p7G_E) + gm_fs5->xsc[p7P_E][p7P_LOOP];
  XMX_FS(1,p7G_C) = XMX_FS(1,p7G_E) + gm_fs5->xsc[p7P_E][p7P_MOVE];

   /* Initialization for row 2 */
  XMX_FS(2,p7G_N) = 0.;
  XMX_FS(2,p7G_B) = gm_fs5->xsc[p7P_N][p7P_MOVE];
  XMX_FS(2,p7G_E) = -eslINFINITY;
  MMX_FS(2,0,p7G_C1) = MMX_FS(2,0,p7G_C2) = MMX_FS(2,0,p7G_C3) = MMX_FS(2,0,p7G_C4) = MMX_FS(2,0,p7G_C5) = -eslINFINITY;
  MMX_FS(2,0,p7G_C0) = IMX_FS(2,0) = DMX_FS(2,0) = -eslINFINITY;

  w = x;
  if(dsq[2] < p7P_MAXNUC) x = dsq[2];
  else                                            x = p7P_MAXCODONS5;

  c1 = p7P_CODON1_FS5(x);
  c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);

  c2 = p7P_CODON2_FS5(w, x);
  c2 = p7P_MINIDX(c2, p7P_DEGEN5_QC1);

  for (k = 1; k <= M; k++) {
    IVX(2,k) = XMX_FS(1,p7G_B) + TSC(p7P_BM,k-1);
    MMX_FS(2,k,p7G_C1) = IVX(2,k) + p7P_MSC_CODON(gm_fs5, k, c1);
    MMX_FS(2,k,p7G_C2) = IVX(1,k) + p7P_MSC_CODON(gm_fs5, k, c2);
    MMX_FS(2,k,p7G_C3) = -eslINFINITY;
    MMX_FS(2,k,p7G_C4) = -eslINFINITY;
    MMX_FS(2,k,p7G_C5) = -eslINFINITY;
    MMX_FS(2,k,p7G_C0) = ESL_MAX( MMX_FS(2,k,p7G_C1), MMX_FS(2,k,p7G_C2));
    IMX_FS(2,k)  = -eslINFINITY;
    DMX_FS(2,k)  = ESL_MAX(MMX_FS(2,k-1,p7G_C0) + TSC(p7P_MD,k-1),
                           DMX_FS(2,k-1)        + TSC(p7P_DD,k-1));

    XMX_FS(2,p7G_E) = ESL_MAX(MMX_FS(2,k,p7G_C0) + esc,
                      ESL_MAX(DMX_FS(2,k)        + esc,
                              XMX_FS(2,p7G_E)));
  }

  XMX_FS(2,p7G_J) = XMX_FS(2,p7G_E) + gm_fs5->xsc[p7P_E][p7P_LOOP];
  XMX_FS(2,p7G_C) = XMX_FS(2,p7G_E) + gm_fs5->xsc[p7P_E][p7P_MOVE];

  t = u = v = p7P_MAXCODONS5;
  /* Initialization for rows 3 and 4 */
  for(i = 3; i < 5; i++)
  {
    u = v;
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder value */
    if(dsq[i] < p7P_MAXNUC) x = dsq[i];
    else                                            x = p7P_MAXCODONS5;

    /* find correct index for looking up scores of codons and quasicodons */
    c1 = p7P_CODON1_FS5(x);
    c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);

    c2 = p7P_CODON2_FS5(w, x);
    c2 = p7P_MINIDX(c2, p7P_DEGEN5_QC1);

    c3 = p7P_CODON3_FS5(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN5_C);

    c4 = p7P_CODON4_FS5(u, v, w, x);
    c4 = p7P_MINIDX(c4, p7P_DEGEN5_QC1);

    ivx_1 = i     % p7P_5CODONS;
    ivx_2 = (i-1) % p7P_5CODONS;
    ivx_3 = (i-2) % p7P_5CODONS;
    ivx_4 = (i-3) % p7P_5CODONS;

    MMX_FS(i,0,p7G_C1) = MMX_FS(i,0,p7G_C2) = MMX_FS(i,0,p7G_C3) = MMX_FS(i,0,p7G_C4) = MMX_FS(i,0,p7G_C5) = -eslINFINITY;
    MMX_FS(i,0,p7G_C0) = IMX_FS(i,0) = DMX_FS(i,0) = -eslINFINITY;

    XMX_FS(i,p7G_E) = -eslINFINITY;

    for (k = 1; k < M; k++)
    {
      IVX(ivx_1,k) = ESL_MAX(MMX_FS(i-1,k-1,p7G_C0) + TSC(p7P_MM,k-1),
                     ESL_MAX(IMX_FS(i-1,k-1)        + TSC(p7P_IM,k-1),
                     ESL_MAX(DMX_FS(i-1,k-1)        + TSC(p7P_DM,k-1),
                             XMX_FS(i-1,p7G_B)      + TSC(p7P_BM,k-1))));

      MMX_FS(i,k,p7G_C1) = IVX(ivx_1,k) + p7P_MSC_CODON(gm_fs5, k, c1);
      MMX_FS(i,k,p7G_C2) = IVX(ivx_2,k) + p7P_MSC_CODON(gm_fs5, k, c2);
      MMX_FS(i,k,p7G_C3) = IVX(ivx_3,k) + p7P_MSC_CODON(gm_fs5, k, c3);
      MMX_FS(i,k,p7G_C4) = -eslINFINITY;
      MMX_FS(i,k,p7G_C5) = -eslINFINITY;
      if( i == 4 )
        MMX_FS(i,k,p7G_C4) = IVX(ivx_4,k) + p7P_MSC_CODON(gm_fs5, k, c4);

      MMX_FS(i,k,p7G_C0) = ESL_MAX( MMX_FS(i,k,p7G_C1),
                           ESL_MAX( MMX_FS(i,k,p7G_C2),
                           ESL_MAX( MMX_FS(i,k,p7G_C3),
                                    MMX_FS(i,k,p7G_C4))));

      IMX_FS(i,k) = ESL_MAX(MMX_FS(i-3,k,p7G_C0) + TSC(p7P_MI,k),
                            IMX_FS(i-3,k)        + TSC(p7P_II,k));

      DMX_FS(i,k) = ESL_MAX(MMX_FS(i,k-1,p7G_C0) + TSC(p7P_MD,k-1),
                            DMX_FS(i,k-1)        + TSC(p7P_DD,k-1));

      XMX_FS(i,p7G_E) = ESL_MAX(MMX_FS(i,k,p7G_C0) + esc,
                        ESL_MAX(DMX_FS(i,k)        + esc,
                                XMX_FS(i,p7G_E)));

    }

    IVX(ivx_1,M) = ESL_MAX(MMX_FS(i-1,M-1,p7G_C0) + TSC(p7P_MM,M-1),
                   ESL_MAX(IMX_FS(i-1,M-1)        + TSC(p7P_IM,M-1),
                   ESL_MAX(DMX_FS(i-1,M-1)        + TSC(p7P_DM,M-1),
                           XMX_FS(i-1,p7G_B)      + TSC(p7P_BM,M-1))));

    MMX_FS(i,M,p7G_C1) = IVX(ivx_1,M) + p7P_MSC_CODON(gm_fs5, M, c1);
    MMX_FS(i,M,p7G_C2) = IVX(ivx_2,M) + p7P_MSC_CODON(gm_fs5, M, c2);
    MMX_FS(i,M,p7G_C3) = IVX(ivx_3,M) + p7P_MSC_CODON(gm_fs5, M, c3);
    MMX_FS(i,M,p7G_C4) = -eslINFINITY;
    MMX_FS(i,M,p7G_C5) = -eslINFINITY;
    if( i == 4 )
      MMX_FS(i,M,p7G_C4) = IVX(ivx_4,M) + p7P_MSC_CODON(gm_fs5, M, c4);

    MMX_FS(i,M,p7G_C0) = ESL_MAX( MMX_FS(i,M,p7G_C1),
                         ESL_MAX( MMX_FS(i,M,p7G_C2),
                         ESL_MAX( MMX_FS(i,M,p7G_C3),
                                  MMX_FS(i,M,p7G_C4))));

    IMX_FS(i,M) = -eslINFINITY;

    DMX_FS(i,M) = ESL_MAX(MMX_FS(i,M-1,p7G_C0) + TSC(p7P_MD,M-1),
                          DMX_FS(i,M-1)        + TSC(p7P_DD,M-1));

    XMX_FS(i,p7G_E) = ESL_MAX(MMX_FS(i,M,p7G_C0),
                      ESL_MAX(DMX_FS(i,M),
                              XMX_FS(i,p7G_E)));

    XMX_FS(i,p7G_J) = ESL_MAX(XMX_FS(i-3,p7G_J) + gm_fs5->xsc[p7P_J][p7P_LOOP],
                              XMX_FS(i,p7G_E)   + gm_fs5->xsc[p7P_E][p7P_LOOP]);

    XMX_FS(i,p7G_C) = ESL_MAX(XMX_FS(i-3,p7G_C) + gm_fs5->xsc[p7P_C][p7P_LOOP],
                              XMX_FS(i,p7G_E)   + gm_fs5->xsc[p7P_E][p7P_MOVE]);

    XMX_FS(i,p7G_N) =         XMX_FS(i-3,p7G_N) + gm_fs5->xsc[p7P_N][p7P_LOOP];

    XMX_FS(i,p7G_B) = ESL_MAX(XMX_FS(i,p7G_N) + gm_fs5->xsc[p7P_N][p7P_MOVE],
                              XMX_FS(i,p7G_J) + gm_fs5->xsc[p7P_J][p7P_MOVE]);

  }


   /* Main DP recursion */
  for (i = 5; i <= L; i++) 
  {
    MMX_FS(i,0,p7G_C0) = MMX_FS(i,0,p7G_C1) = MMX_FS(i,0,p7G_C2) = MMX_FS(i,0,p7G_C3)
    = MMX_FS(i,0,p7G_C4) = MMX_FS(i,0,p7G_C5) = IMX_FS(i,0) = DMX_FS(i,0) = -eslINFINITY;

    XMX_FS(i, p7G_E) = -eslINFINITY;

    /* Reasign nucleotide to correct temporary holders for use in emissions array */
    t = u;
    u = v;
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder value */
    if(dsq[i] < p7P_MAXNUC) x = dsq[i];
    else                                            x = p7P_MAXCODONS5;

    /* find correct index for looking up scores of codons and quasicodons */
    c1 = p7P_CODON1_FS5(x);
    c1 = p7P_MINIDX(c1, p7P_DEGEN5_QC2);

    c2 = p7P_CODON2_FS5(w, x);
    c2 = p7P_MINIDX(c2, p7P_DEGEN5_QC1);

    c3 = p7P_CODON3_FS5(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN5_C);

    c4 = p7P_CODON4_FS5(u, v, w, x);
    c4 = p7P_MINIDX(c4, p7P_DEGEN5_QC1);

    c5 = p7P_CODON5_FS5(t, u, v, w, x);
    c5 = p7P_MINIDX(c5, p7P_DEGEN5_QC2);

    ivx_1 = i     % p7P_5CODONS;
    ivx_2 = (i-1) % p7P_5CODONS;
    ivx_3 = (i-2) % p7P_5CODONS;
    ivx_4 = (i-3) % p7P_5CODONS;
    ivx_5 = (i-4) % p7P_5CODONS;

    for (k = 1; k < M; k++)
    {

      IVX(ivx_1,k) = ESL_MAX(MMX_FS(i-1,k-1,p7G_C0)   + TSC(p7P_MM,k-1),
                     ESL_MAX(IMX_FS(i-1,k-1)          + TSC(p7P_IM,k-1),
                     ESL_MAX(DMX_FS(i-1,k-1)          + TSC(p7P_DM,k-1),
                             XMX_FS(i-1,p7G_B)        + TSC(p7P_BM,k-1))));

      MMX_FS(i,k,p7G_C1) = IVX(ivx_1,k) + p7P_MSC_CODON(gm_fs5, k, c1);

      MMX_FS(i,k,p7G_C2) = IVX(ivx_2,k) + p7P_MSC_CODON(gm_fs5, k, c2);

      MMX_FS(i,k,p7G_C3) = IVX(ivx_3,k) + p7P_MSC_CODON(gm_fs5, k, c3);

      MMX_FS(i,k,p7G_C4) = IVX(ivx_4,k) + p7P_MSC_CODON(gm_fs5, k, c4);

      MMX_FS(i,k,p7G_C5) = IVX(ivx_5,k) + p7P_MSC_CODON(gm_fs5, k, c5);

      MMX_FS(i,k,p7G_C0) =  ESL_MAX(ESL_MAX(MMX_FS(i,k,p7G_C1),
                            ESL_MAX(MMX_FS(i,k,p7G_C2), MMX_FS(i,k,p7G_C3))),
                            ESL_MAX(MMX_FS(i,k,p7G_C4), MMX_FS(i,k,p7G_C5)));

      /* insert state */
      IMX_FS(i,k) = ESL_MAX(MMX_FS(i-3,k,p7G_C0) + TSC(p7P_MI,k),
                            IMX_FS(i-3,k)        + TSC(p7P_II,k));

      /* delete state */
      DMX_FS(i,k) = ESL_MAX(MMX_FS(i,k-1,p7G_C0) + TSC(p7P_MD,k-1),
                            DMX_FS(i,k-1)        + TSC(p7P_DD,k-1));

      /* E state update */
      XMX_FS(i,p7G_E) = ESL_MAX(MMX_FS(i,k,p7G_C0) + esc,
                        ESL_MAX(DMX_FS(i,k)        + esc,
                                XMX_FS(i,p7G_E)));
    }

    /* unrolled match state M_M */
    IVX(ivx_1,M) = ESL_MAX(MMX_FS(i-1,M-1,p7G_C0)   + TSC(p7P_MM,M-1),
                   ESL_MAX(IMX_FS(i-1,M-1)          + TSC(p7P_IM,M-1),
                   ESL_MAX(DMX_FS(i-1,M-1)          + TSC(p7P_DM,M-1),
                           XMX_FS(i-1,p7G_B)        + TSC(p7P_BM,M-1))));

    MMX_FS(i,M,p7G_C1) = IVX(ivx_1,M) + p7P_MSC_CODON(gm_fs5, M, c1);

    MMX_FS(i,M,p7G_C2) = IVX(ivx_2,M) + p7P_MSC_CODON(gm_fs5, M, c2);

    MMX_FS(i,M,p7G_C3) = IVX(ivx_3,M) + p7P_MSC_CODON(gm_fs5, M, c3);

    MMX_FS(i,M,p7G_C4) = IVX(ivx_4,M) + p7P_MSC_CODON(gm_fs5, M, c4);

    MMX_FS(i,M,p7G_C5) = IVX(ivx_5,M) + p7P_MSC_CODON(gm_fs5, M, c5);

    MMX_FS(i,M,p7G_C0) =  ESL_MAX(ESL_MAX(MMX_FS(i,M,p7G_C1),
                          ESL_MAX(MMX_FS(i,M,p7G_C2), MMX_FS(i,M,p7G_C3))),
                          ESL_MAX(MMX_FS(i,M,p7G_C4), MMX_FS(i,M,p7G_C5)));

    IMX_FS(i,M) = -eslINFINITY;

    /* unrolled delete state D_M */
    DMX_FS(i,M) = ESL_MAX(MMX_FS(i,M-1,p7G_C0) + TSC(p7P_MD,M-1),
                          DMX_FS(i,M-1) + TSC(p7P_DD,M-1));

    /* unrolled E state update */
    XMX_FS(i,p7G_E) = ESL_MAX(ESL_MAX(MMX_FS(i,M,p7G_C0),
                                      DMX_FS(i,M)),
                                      XMX_FS(i,p7G_E));

    /* J, C and N states */
    XMX_FS(i,p7G_J) = ESL_MAX(XMX_FS(i-3,p7G_J) + gm_fs5->xsc[p7P_J][p7P_LOOP],
                              XMX_FS(i,p7G_E)   + gm_fs5->xsc[p7P_E][p7P_LOOP]);

    XMX_FS(i,p7G_C) = ESL_MAX(XMX_FS(i-3,p7G_C) + gm_fs5->xsc[p7P_C][p7P_LOOP],
                              XMX_FS(i,p7G_E)   + gm_fs5->xsc[p7P_E][p7P_MOVE]);
    XMX_FS(i,p7G_N) =         XMX_FS(i-3,p7G_N) + gm_fs5->xsc[p7P_N][p7P_LOOP];
    XMX_FS(i,p7G_B) = ESL_MAX(XMX_FS(i,p7G_N) + gm_fs5->xsc[p7P_N][p7P_MOVE],
                              XMX_FS(i,p7G_J) + gm_fs5->xsc[p7P_J][p7P_MOVE]);
  }
  
  /* T state (not stored) */
  if (opt_sc != NULL) *opt_sc = ESL_MAX( XMX_FS(L,p7G_C),
                                ESL_MAX( XMX_FS(L-1,p7G_C) + gm_fs5->xsc[p7P_C][p7P_LOOP],
                                         XMX_FS(L-2,p7G_C) + gm_fs5->xsc[p7P_C][p7P_LOOP])) + gm_fs5->xsc[p7P_C][p7P_MOVE];

  gx->M = gm_fs5->M;
  gx->L = L;

  return eslOK;

}



/* Function: p7_GVTrace_Frameshift()
 *
 * Purpose:  Traceback of a Frameshift aware Viterbi matrix: retrieval
 *           of optimum alignment.
 *
 * Args:     dsq    - digital sequence aligned to, 1..L
 *           L      - length of <dsq>
 *           gm     - profile
 *           mx     - Viterbi matrix to trace, L x M
 *           tr     - storage for the recovered traceback.
 *
 * Return:   <eslOK> on success.
 *           <eslFAIL> if even the optimal path has zero probability;
 *           in this case, the trace is set blank (<tr->N = 0>).
 *
 */
int
p7_GVTrace_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs5, const P7_GMX *gx, P7_TRACE *tr)
{
  int          i   = L;     /* position in seq (1..L)         */
  int          k   = 0;     /* position in model (1..M)       */
  int          M   = gm_fs5->M;
  int          c   = 0;
  int          prev_c;
  float      **dp  = gx->dp;    /* so {MDI}MX() macros work       */
  float       *xmx = gx->xmx;   /* so XMX() macro works           */
  float        tol = 1e-5;  /* floating point "equality" test */
  float const *tsc = gm_fs5->tsc;
  int     sprv, scur;       /* previous, current state in trace */
  float   path[4];
  int     state[4] = {p7T_M, p7T_I, p7T_D, p7T_B};
  float   match_codon[5];
  int     status;

  if(gm_fs5->codon_lengths != 5) ESL_EXCEPTION(eslEINVAL, "proflie not allocated for 5 codon lengths");

#if eslDEBUGLEVEL > 0
  if (tr->N != 0) ESL_EXCEPTION(eslEINVAL, "trace isn't empty: forgot to Reuse()?");
#endif

  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i, c)) != eslOK) return status;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i, c)) != eslOK) return status;

  sprv = p7T_C;
 
  while (sprv != p7T_S) {

    switch (sprv) {
    case p7T_C:     /* C(i) comes from C(i-1) or E(i) */
      if   (XMX_FS(i,p7G_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);

      if      (XMX_FS(i, p7G_C) < XMX_FS(i-2, p7G_C) || XMX_FS(i, p7G_C) < XMX_FS(i-1, p7G_C))                      scur = p7T_C;
      else if (esl_FCompare_old(XMX_FS(i, p7G_C), XMX_FS(i-3, p7G_C) + gm_fs5->xsc[p7P_C][p7P_LOOP], tol) == eslOK)  scur = p7T_C;
      else if (esl_FCompare_old(XMX_FS(i, p7G_C), XMX_FS(i,   p7G_E) + gm_fs5->xsc[p7P_E][p7P_MOVE], tol) == eslOK)  scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
 
      break;

    case p7T_E:     /* E connects from any M state. k set here */
      if (XMX_FS(i, p7G_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);

      if (p7_fs_profile_IsLocal(gm_fs5))
      {
        scur = p7T_M;     /* can't come from D, in a *local* Viterbi trace. */
        for (k = M; k >= 1; k--) if (esl_FCompare_old(XMX_FS(i, p7G_E), MMX_FS(i,k,p7G_C0), tol) == eslOK) break;
        if (k == 0) ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      }
      else          /* glocal mode: we either come from D_M or M_M */
      {
        if      (esl_FCompare_old(XMX_FS(i, p7G_E), MMX_FS(i,M,p7G_C0), tol) == eslOK) { scur = p7T_M; k = M; }
        else if (esl_FCompare_old(XMX_FS(i, p7G_E), DMX_FS(i,M),        tol) == eslOK) { scur = p7T_D; k = M; }
        else    ESL_EXCEPTION(eslFAIL, "E at i=%d couldn't be traced", i);
      }
    
      
      break;

    case p7T_M:         /* M connects from i-1,k-1, or B */
      if (MMX_FS(i,k, p7G_C0) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);

      path[0] = MMX_FS(i-prev_c,k-1,p7G_C0) + TSC(p7P_MM, k-1);
      path[1] = IMX_FS(i-prev_c,k-1)        + TSC(p7P_IM, k-1);
      path[2] = DMX_FS(i-prev_c,k-1)        + TSC(p7P_DM, k-1);
      path[3] = XMX_FS(i-prev_c,p7G_B)      + TSC(p7P_BM, k-1);
      
      scur = state[esl_vec_FArgMax(path, 4)]; 
      k--; i-=prev_c;
      break;

    case p7T_D:         /* D connects from M,D at i,k-1 */
      if (DMX_FS(i, k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);

      if      (esl_FCompare_old(DMX_FS(i,k), MMX_FS(i, k-1, p7G_C0) + TSC(p7P_MD, k-1), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(DMX_FS(i,k), DMX_FS(i, k-1)         + TSC(p7P_DD, k-1), tol) == eslOK) scur = p7T_D;
      else ESL_EXCEPTION(eslFAIL, "D at k=%d,i=%d couldn't be traced", k,i);
      k--;
      break;

    case p7T_I:         /* I connects from M,I at i-1,k*/
      if (IMX_FS(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);

      if      (esl_FCompare_old(IMX_FS(i,k), MMX_FS(i-3,k,p7G_C0) + TSC(p7P_MI, k), tol) == eslOK) scur = p7T_M;
      else if (esl_FCompare_old(IMX_FS(i,k), IMX_FS(i-3,k)        + TSC(p7P_II, k), tol) == eslOK) scur = p7T_I;
      else ESL_EXCEPTION(eslFAIL, "I at k=%d,i=%d couldn't be traced", k,i);
      i-=3;
      break;

    case p7T_N:         /* N connects from S, N */
      if (XMX_FS(i, p7G_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = ( (i == 0) ? p7T_S : p7T_N);
      break;

    case p7T_B:         /* B connects from N, J */
      if (XMX_FS(i,p7G_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

      if      (esl_FCompare_old(XMX_FS(i,p7G_B), XMX_FS(i, p7G_N) + gm_fs5->xsc[p7P_N][p7P_MOVE], tol) == eslOK) scur = p7T_N;
      else if (esl_FCompare_old(XMX_FS(i,p7G_B), XMX_FS(i, p7G_J) + gm_fs5->xsc[p7P_J][p7P_MOVE], tol) == eslOK) scur = p7T_J;
      else  ESL_EXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    case p7T_J:         /* J connects from E(i) or J(i-1) */
      if (XMX_FS(i,p7G_J) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible J reached at i=%d", i);
   
      if      (esl_FCompare_old(XMX_FS(i,p7G_J), XMX_FS(i-3,p7G_J) + gm_fs5->xsc[p7P_J][p7P_LOOP], tol) == eslOK) scur = p7T_J;
      else if (esl_FCompare_old(XMX_FS(i,p7G_J), XMX_FS(i,  p7G_E) + gm_fs5->xsc[p7P_E][p7P_LOOP], tol) == eslOK) scur = p7T_E;
      else  ESL_EXCEPTION(eslFAIL, "J at i=%d couldn't be traced", i);
      break;

    default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over statetype[tpos-1] */

    if(scur == p7T_M) {
       match_codon[0] = MMX_FS(i,k,p7G_C1);
       match_codon[1] = MMX_FS(i,k,p7G_C2);
       match_codon[2] = MMX_FS(i,k,p7G_C3);
       match_codon[3] = MMX_FS(i,k,p7G_C4);
       match_codon[4] = MMX_FS(i,k,p7G_C5);

       c = esl_vec_FArgMax(match_codon, 5) + 1;
    }
   
    if ((status = p7_trace_fs_Append(tr, scur, k, i, c)) != eslOK) return status;
  
    /* For NCJ, we had to defer i decrement. */
    if ( (scur == p7T_N || scur == p7T_C) && scur == sprv) i--; 
    if (  scur == p7T_J && scur == sprv) i-=3;
    prev_c = c;
    c = 0;
    sprv = scur;
  } /* end traceback, at S state */

  tr->M = gm_fs5->M;
  tr->L = L;

  return p7_trace_fs_Reverse(tr);

}

/*-------------------- end, Viterbi -----------------------------*/

/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/


/*----------------- end, benchmark ------------------------------*/



/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7GENERIC_VITERBI_FRAMESHIFT_TESTDRIVE
#include <string.h>
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_basic: build a tiny protein HMM from a Stockholm alignment,
 * reverse-translate one codon per residue into a DNA target sequence,
 * run frameshift Viterbi, run the traceback, and validate the trace.
 * Also confirms that Forward score >= Viterbi score on the same input.
 */
static void
utest_basic(ESL_GETOPTS *go)
{
  /* Minimal Stockholm alignment of two identical 4-residue seqs */
  char           *query = "# STOCKHOLM 1.0\n\nseq1 MAFY\nseq2 MAFY\n//\n";
  int             fmt   = eslMSAFILE_STOCKHOLM;
  ESL_ALPHABET   *abcAA  = NULL;
  ESL_ALPHABET   *abcDNA = NULL;
  ESL_MSA        *msa    = NULL;
  P7_HMM         *hmm    = NULL;
  P7_BG          *bgAA   = NULL;
  P7_FS_PROFILE  *gm_fs5 = NULL;
  P7_PRIOR       *pri    = NULL;
  ESL_GENCODE    *gcode  = NULL;
  P7_CODONTABLE  *ct     = NULL;
  ESL_RANDOMNESS *r      = NULL;
  /* DNA target: ATG(M) GCT(A) TTT(F) TAT(Y) = MAFY */
  char           *targ   = "ATGGCTTTTTAT";
  int             L      = strlen(targ);  /* 12 nt */
  ESL_DSQ        *dsq    = NULL;
  P7_GMX         *vit    = NULL;
  P7_GMX         *fwd    = NULL;
  P7_IVX         *iv     = NULL;
  P7_TRACE       *tr     = NULL;
  char            errbuf[eslERRBUFSIZE];
  float           vsc, fsc;

  if ((abcAA  = esl_alphabet_Create(eslAMINO))          == NULL) esl_fatal("failed to create AA alphabet");
  if ((abcDNA = esl_alphabet_Create(eslDNA))             == NULL) esl_fatal("failed to create DNA alphabet");
  if ((r      = esl_randomness_CreateFast(42))           == NULL) esl_fatal("failed to create RNG");
  if ((pri    = p7_prior_CreateAmino())                  == NULL) esl_fatal("failed to create prior");
  if ((msa    = esl_msa_CreateFromString(query, fmt))    == NULL) esl_fatal("failed to create MSA");
  if (esl_msa_Digitize(abcAA, msa, NULL)                != eslOK) esl_fatal("failed to digitize MSA");
  if (p7_Fastmodelmaker(msa, 0.5, NULL, &hmm, NULL)     != eslOK) esl_fatal("failed to build HMM");
  if (p7_ParameterEstimation(hmm, pri)                  != eslOK) esl_fatal("failed to parameterize HMM");
  if (p7_hmm_SetConsensus(hmm, NULL)                    != eslOK) esl_fatal("failed to set consensus");
  if ((bgAA   = p7_bg_Create(abcAA))                    == NULL) esl_fatal("failed to create background");
  if ((gcode  = esl_gencode_Create(abcDNA, abcAA))      == NULL) esl_fatal("failed to create gencode");
  if ((ct     = p7_codontable_Create(gcode))             == NULL) esl_fatal("failed to create codon table");
  if ((gm_fs5 = p7_profile_fs_Create(hmm->M, abcAA, 5)) == NULL) esl_fatal("failed to create fs profile");
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_UNILOCAL) != eslOK) esl_fatal("failed to configure fs profile");

  if (esl_abc_CreateDsq(abcDNA, targ, &dsq)             != eslOK) esl_fatal("failed to digitize DNA sequence");

  if ((vit = p7_gmx_Create(gm_fs5->M, L, L, p7P_5CODONS))    == NULL) esl_fatal("failed to create Viterbi matrix");
  if ((fwd = p7_gmx_Create(gm_fs5->M, L, L, p7G_NSCELLS_FS)) == NULL) esl_fatal("failed to create Forward matrix");
  if ((iv  = p7_ivx_Create(gm_fs5->M, p7P_5CODONS))           == NULL) esl_fatal("failed to create IVX");
  if ((tr  = p7_trace_Create())                                == NULL) esl_fatal("failed to create trace");

  if (p7_GViterbi_Frameshift(dsq, L, gm_fs5, vit, iv, &vsc) != eslOK) esl_fatal("Viterbi failed");
  if (esl_opt_GetBoolean(go, "-v")) printf("utest_basic: Viterbi score: %.4f\n", vsc);

  if (p7_GVTrace_Frameshift(dsq, L, gm_fs5, vit, tr)        != eslOK) esl_fatal("traceback failed");
  if (p7_trace_fs_Validate(tr, abcDNA, dsq, errbuf)         != eslOK) esl_fatal("trace invalid: %s", errbuf);
  if (esl_opt_GetBoolean(go, "-v")) p7_trace_Dump(stdout, tr, NULL, dsq);

  p7_FLogsumInit();
  if (p7_GForward_Frameshift(dsq, L, gm_fs5, fwd, iv, &fsc) != eslOK) esl_fatal("Forward failed");
  if (esl_opt_GetBoolean(go, "-v")) printf("utest_basic: Forward score: %.4f\n", fsc);
  if (fsc < vsc - 0.001) esl_fatal("utest_basic: Forward score (%.4f) < Viterbi score (%.4f)", fsc, vsc);

  p7_trace_Destroy(tr);
  p7_ivx_Destroy(iv);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(vit);
  free(dsq);
  p7_profile_fs_Destroy(gm_fs5);
  p7_codontable_Destroy(ct);
  esl_gencode_Destroy(gcode);
  p7_bg_Destroy(bgAA);
  p7_hmm_Destroy(hmm);
  esl_msa_Destroy(msa);
  p7_prior_Destroy(pri);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
}


/* utest_viterbi_fs: Viterbi validation on random DNA sequences.
 *
 * For each sequence:
 *   1. Reverse-translate a random AA sequence to get in-frame DNA.
 *   2. Run Viterbi and traceback; validate the trace structure.
 *   3. Confirm Forward score >= Viterbi score.
 *   4. Accumulate (Viterbi - null) scores to check the average is <= 0
 *      (random sequences should not score better than the null model).
 */
static void
utest_viterbi_fs(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA,
                 P7_CODONTABLE *ct, P7_BG *bgAA, P7_FS_PROFILE *gm_fs5, int nseq, int L)
{
  int       i, j, idx;
  float     avg_sc = 0.;
  ESL_DSQ  *dsqAA  = NULL;
  ESL_DSQ  *dsqDNA = NULL;
  P7_GMX   *vit    = NULL;
  P7_GMX   *fwd    = NULL;
  P7_IVX   *iv     = NULL;
  P7_TRACE *tr     = NULL;
  char      errbuf[eslERRBUFSIZE];
  float     vsc, fsc, nullsc;

  if ((dsqAA  = malloc(sizeof(ESL_DSQ) * ((L/3)+2))) == NULL) esl_fatal("malloc failed");
  if ((dsqDNA = malloc(sizeof(ESL_DSQ) * (L+2)))     == NULL) esl_fatal("malloc failed");
  if ((vit    = p7_gmx_Create(gm_fs5->M, L, L, p7P_5CODONS))    == NULL) esl_fatal("matrix creation failed");
  if ((fwd    = p7_gmx_Create(gm_fs5->M, L, L, p7G_NSCELLS_FS)) == NULL) esl_fatal("matrix creation failed");
  if ((iv     = p7_ivx_Create(gm_fs5->M, p7P_5CODONS))           == NULL) esl_fatal("ivx creation failed");
  if ((tr     = p7_trace_Create())                                == NULL) esl_fatal("trace creation failed");

  for (idx = 0; idx < nseq; idx++)
    {
      /* Build a random in-frame DNA sequence from a random AA sequence */
      if (esl_rsq_xfIID(r, bgAA->f, abcAA->K, L/3, dsqAA) != eslOK) esl_fatal("seq generation failed");
      j = 1;
      for (i = 1; i <= L/3; i++) {
        p7_codontable_GetCodon(ct, r, dsqAA[i], dsqDNA+j);
        j += 3;
      }
      dsqDNA[0] = dsqDNA[L+1] = eslDSQ_SENTINEL;

      if (p7_GViterbi_Frameshift(dsqDNA, L, gm_fs5, vit, iv, &vsc) != eslOK) esl_fatal("Viterbi failed");
      if (p7_GVTrace_Frameshift(dsqDNA, L, gm_fs5, vit, tr)        != eslOK) esl_fatal("traceback failed");
      if (p7_trace_fs_Validate(tr, abcDNA, dsqDNA, errbuf)         != eslOK) esl_fatal("trace invalid: %s", errbuf);

      if (p7_GForward_Frameshift(dsqDNA, L, gm_fs5, fwd, iv, &fsc) != eslOK) esl_fatal("Forward failed");
      if (fsc < vsc - 0.001) esl_fatal("Forward score (%.4f) < Viterbi score (%.4f)", fsc, vsc);

      if (p7_bg_fs_NullOne(bgAA, dsqAA, L/3, &nullsc) != eslOK) esl_fatal("null score failed");
      avg_sc += vsc - nullsc;

      if (esl_opt_GetBoolean(go, "--vv"))
        printf("utest_viterbi_fs: vsc %.4f fsc %.4f nullsc %.4f\n", vsc, fsc, nullsc);

      p7_trace_Reuse(tr);
      p7_gmx_Reuse(vit);
      p7_gmx_Reuse(fwd);
    }

  avg_sc /= (float) nseq;
  if (avg_sc > 0.) esl_fatal("Viterbi scores have positive expectation (%.4f nats)", avg_sc);

  p7_trace_Destroy(tr);
  p7_ivx_Destroy(iv);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(vit);
  free(dsqDNA);
  free(dsqAA);
}

#endif /*p7GENERIC_VITERBI_FRAMESHIFT_TESTDRIVE*/
/*----------------- end, unit tests -----------------------------*/


/*****************************************************************
 * 4. Test driver.
 *****************************************************************/
/* gcc -g -Wall -Dp7GENERIC_VITERBI_FRAMESHIFT_TESTDRIVE -I. -I../easel -L. -L../easel -o generic_viterbi_frameshift_utest generic_viterbi_frameshift.c -lhmmer -leasel -lm
 */
#ifdef p7GENERIC_VITERBI_FRAMESHIFT_TESTDRIVE
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
static char banner[] = "unit test driver for the generic frameshift Viterbi implementation";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA  = NULL;
  ESL_ALPHABET   *abcDNA = NULL;
  P7_HMM         *hmm    = NULL;
  P7_FS_PROFILE  *gm_fs5 = NULL;
  P7_BG          *bgAA   = NULL;
  ESL_GENCODE    *gcode  = NULL;
  P7_CODONTABLE  *ct     = NULL;
  int             M      = 100;
  int             L      = 600;  /* must be a multiple of 3 */
  int             nseq   = 20;
  char            errbuf[eslERRBUFSIZE];

  p7_FLogsumInit();

  if ((abcAA  = esl_alphabet_Create(eslAMINO))                         == NULL) esl_fatal("failed to create AA alphabet");
  if ((abcDNA = esl_alphabet_Create(eslDNA))                           == NULL) esl_fatal("failed to create DNA alphabet");
  if ((gcode  = esl_gencode_Create(abcDNA, abcAA))                     == NULL) esl_fatal("failed to create gencode");
  if ((ct     = p7_codontable_Create(gcode))                           == NULL) esl_fatal("failed to create codon table");
  if (p7_hmm_Sample(r, M, abcAA, &hmm)                                != eslOK) esl_fatal("failed to sample HMM");
  if ((bgAA   = p7_bg_Create(abcAA))                                   == NULL) esl_fatal("failed to create background");
  if (p7_bg_SetLength(bgAA, L/3)                                       != eslOK) esl_fatal("failed to set bg length");
  if ((gm_fs5 = p7_profile_fs_Create(hmm->M, abcAA, 5))                == NULL) esl_fatal("failed to create fs profile");
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_LOCAL)    != eslOK) esl_fatal("failed to configure fs profile");
  if (p7_hmm_Validate(hmm, errbuf, 0.0001)                             != eslOK) esl_fatal("HMM invalid: %s", errbuf);

  utest_basic      (go);
  utest_viterbi_fs (go, r, abcAA, abcDNA, ct, bgAA, gm_fs5, nseq, L);

  fprintf(stderr, "All tests passed.\n");

  p7_profile_fs_Destroy(gm_fs5);
  p7_bg_Destroy(bgAA);
  p7_hmm_Destroy(hmm);
  p7_codontable_Destroy(ct);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7GENERIC_VITERBI_FRAMESHIFT_TESTDRIVE*/
/*----------------- end, test driver -----------------------------*/
