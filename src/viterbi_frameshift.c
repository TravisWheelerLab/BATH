/* Frameshift aware Viterbi algorithm and traceback.*/
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#define IVX(i,k) (ivx[((k)*p7P_5CODONS) + (i)])

/* Function:  p7_Viterbi_Frameshift()
 * Synopsis:  The Viterbi algorithm, with frameshift awareness.
 * 
 * Purpose:   The Viterbi dynamic programming algorithm for aligning 
 *            DNA to proteins with frameshift awareness. 
 *
 *            Given a digital nucleotide sequence <dsq> of length <L>, 
 *            a frameshift aware codon profile <gm_fs>, and DP matrix 
 *            <gx> allocated for at least <L> by <gm_fs->M> cells; 
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
 *            gcode  - genetic code for translation
 *            L      - length of dsq
 *            gm_fs  - profile. 
 *            gx     - DP matrix with room for an MxL alignment
 *            iv     - intermediate value matrix
 *            opt_sc - optRETURN: Viterbi lod score in nats
 *           
 * Return:   <eslOK> on success.
 */
int
p7_Viterbi_Frameshift(const ESL_DSQ *dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx, P7_IVX *iv, float *opt_sc)
{
  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float       *ivx  = iv->ivx;
  int          M    = gm_fs->M;
  int          i,k,c;
  int          c1, c2, c3, c4, c5;
  int          t, u, v, w, x;
  int          ivx_1, ivx_2, ivx_3, ivx_4, ivx_5;
  float        esc  = p7_fs_profile_IsLocal(gm_fs) ? 0 : -eslINFINITY;
  
  /* Initialization of the zero row.  */
  XMX_FS(0,p7G_N) = 0.;                                  /* S->N, p=1            */
  XMX_FS(0,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE];         /* S->N->B, no N-tail   */
  XMX_FS(0,p7G_E) = XMX_FS(0,p7G_J) = XMX_FS(0,p7G_C) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_FS(0,k,p7G_C0) = MMX_FS(0,k,p7G_C1) = MMX_FS(0,k,p7G_C2) = MMX_FS(0,k,p7G_C3) =
    MMX_FS(0,k,p7G_C4) = MMX_FS(0,k,p7G_C5) = IMX_FS(0,k)        = DMX_FS(0,k)        = -eslINFINITY;

   /* Initialization for row 1 */
  XMX_FS(1,p7G_N) = 0.;
  XMX_FS(1,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE];
  XMX_FS(1,p7G_E) = -eslINFINITY;
  MMX_FS(1,0,p7G_C1) = MMX_FS(1,0,p7G_C2) = MMX_FS(1,0,p7G_C3) = MMX_FS(1,0,p7G_C4) = MMX_FS(1,0,p7G_C5) = -eslINFINITY;
  MMX_FS(1,0,p7G_C0) = IMX_FS(1,0) = DMX_FS(1,0) = -eslINFINITY;

  if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[1])) x = dsq[1];
  else                                            x = p7P_MAXCODONS;

  c1 = p7P_CODON1(x);
  c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);
  for (k = 1; k <= M; k++) {
    IVX(1,k) = XMX_FS(0,p7G_B) + TSC(p7P_BM,k-1);

    MMX_FS(1,k,p7G_C1) = IVX(1,k) + p7P_MSC_CODON(gm_fs, k, c1);
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

  XMX_FS(1,p7G_J) = XMX_FS(1,p7G_E) + gm_fs->xsc[p7P_E][p7P_LOOP];
  XMX_FS(1,p7G_C) = XMX_FS(1,p7G_E) + gm_fs->xsc[p7P_E][p7P_MOVE];

   /* Initialization for row 2 */
  XMX_FS(2,p7G_N) = 0.;
  XMX_FS(2,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE];
  XMX_FS(2,p7G_E) = -eslINFINITY;
  MMX_FS(2,0,p7G_C1) = MMX_FS(2,0,p7G_C2) = MMX_FS(2,0,p7G_C3) = MMX_FS(2,0,p7G_C4) = MMX_FS(2,0,p7G_C5) = -eslINFINITY;
  MMX_FS(2,0,p7G_C0) = IMX_FS(2,0) = DMX_FS(2,0) = -eslINFINITY;

  w = x;
  if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[2])) x = dsq[2];
  else                                            x = p7P_MAXCODONS;

  c1 = p7P_CODON1(x);
  c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

  c2 = p7P_CODON2(w, x);
  c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

  for (k = 1; k <= M; k++) {
    IVX(2,k) = XMX_FS(1,p7G_B) + TSC(p7P_BM,k-1);
    MMX_FS(2,k,p7G_C1) = IVX(2,k) + p7P_MSC_CODON(gm_fs, k, c1);
    MMX_FS(2,k,p7G_C2) = IVX(1,k) + p7P_MSC_CODON(gm_fs, k, c2);
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

  XMX_FS(2,p7G_J) = XMX_FS(2,p7G_E) + gm_fs->xsc[p7P_E][p7P_LOOP];
  XMX_FS(2,p7G_C) = XMX_FS(2,p7G_E) + gm_fs->xsc[p7P_E][p7P_MOVE];

  t = u = v = p7P_MAXCODONS;
  /* Initialization for rows 3 and 4 */
  for(i = 3; i < 5; i++)
  {
    u = v;
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder value */
    if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[i])) x = dsq[i];
    else                                            x = p7P_MAXCODONS;

    /* find correct index for looking up scores of codons and quasicodons */
    c1 = p7P_CODON1(x);
    c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

    c2 = p7P_CODON2(w, x);
    c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    c4 = p7P_CODON4(u, v, w, x);
    c4 = p7P_MINIDX(c4, p7P_DEGEN_QC1);

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

      MMX_FS(i,k,p7G_C1) = IVX(ivx_1,k) + p7P_MSC_CODON(gm_fs, k, c1);
      MMX_FS(i,k,p7G_C2) = IVX(ivx_2,k) + p7P_MSC_CODON(gm_fs, k, c2);
      MMX_FS(i,k,p7G_C3) = IVX(ivx_3,k) + p7P_MSC_CODON(gm_fs, k, c3);
      MMX_FS(i,k,p7G_C4) = -eslINFINITY;
      MMX_FS(i,k,p7G_C5) = -eslINFINITY;
      if( i == 4 )
        MMX_FS(i,k,p7G_C4) = IVX(ivx_4,k) + p7P_MSC_CODON(gm_fs, k, c4);

      MMX_FS(i,k,p7G_C0) = ESL_MAX( MMX_FS(i,k,p7G_C1),
                           ESL_MAX( MMX_FS(i,k,p7G_C2),
                           ESL_MAX( MMX_FS(i,k,p7G_C3),
                                    MMX_FS(i,k,p7G_C4))));

      IMX_FS(i,k) = ESL_MAX(MMX_FS(i-3,k,p7G_C0) + TSC(p7P_MI,k),
                            IMX_FS(i-3,k)        + TSC(p7P_II,k));

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

    MMX_FS(i,M,p7G_C1) = IVX(ivx_1,M) + p7P_MSC_CODON(gm_fs, M, c1);
    MMX_FS(i,M,p7G_C2) = IVX(ivx_2,M) + p7P_MSC_CODON(gm_fs, M, c2);
    MMX_FS(i,M,p7G_C3) = IVX(ivx_3,M) + p7P_MSC_CODON(gm_fs, M, c3);
    MMX_FS(i,M,p7G_C4) = -eslINFINITY;
    MMX_FS(i,M,p7G_C5) = -eslINFINITY;
    if( i == 4 )
      MMX_FS(i,M,p7G_C4) = IVX(ivx_4,M) + p7P_MSC_CODON(gm_fs, M, c4);

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

    XMX_FS(i,p7G_J) = ESL_MAX(XMX_FS(i-3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP],
                              XMX_FS(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_LOOP]);

    XMX_FS(i,p7G_C) = ESL_MAX(XMX_FS(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                              XMX_FS(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);

    XMX_FS(i,p7G_N) =         XMX_FS(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];

    XMX_FS(i,p7G_B) = ESL_MAX(XMX_FS(i,p7G_N) + gm_fs->xsc[p7P_N][p7P_MOVE],
                              XMX_FS(i,p7G_J) + gm_fs->xsc[p7P_J][p7P_MOVE]);

  }


   /* Main DP recursion */
  for (i = 5; i <= L; i++) 
  {
    MMX_FS(i,0,p7G_C0) = MMX_FS(i,0,p7G_C1) = MMX_FS(i,0,p7G_C2) = MMX_FS(i,0,p7G_C3)
    = MMX_FS(i,0,p7G_C4) = MMX_FS(i,0,p7G_C5) = IMX_FS(i,0) = DMX_FS(i,0) = -eslINFINITY;

    XMX_FS(i, p7G_E) = -eslINFINITY;

    /* Reasign nucluotide to correct temporary holders for use in emissions array */
    t = u;
    u = v;
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder vlaue */
    if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[i])) x = dsq[i];
    else                                            x = p7P_MAXCODONS;

    /* find correct index for looking up scores of codons and quasicodons */
    c1 = p7P_CODON1(x);
    c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

    c2 = p7P_CODON2(w, x);
    c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

    c3 = p7P_CODON3(v, w, x);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    c4 = p7P_CODON4(u, v, w, x);
    c4 = p7P_MINIDX(c4, p7P_DEGEN_QC1);

    c5 = p7P_CODON5(t, u, v, w, x);
    c5 = p7P_MINIDX(c5, p7P_DEGEN_QC2);

    ivx_1 = i     % p7P_5CODONS;
    ivx_2 = (i-1) % p7P_5CODONS;
    ivx_3 = (i-2) % p7P_5CODONS;
    ivx_4 = (i-3) % p7P_5CODONS;
    ivx_5 = (i-5) % p7P_5CODONS;

    for (k = 1; k < M; k++)
    {

      IVX(ivx_1,k) = ESL_MAX(MMX_FS(i-1,k-1,p7G_C0)   + TSC(p7P_MM,k-1),
                     ESL_MAX(IMX_FS(i-1,k-1)          + TSC(p7P_IM,k-1),
                     ESL_MAX(DMX_FS(i-1,k-1)          + TSC(p7P_DM,k-1),
                             XMX_FS(i-1,p7G_B)        + TSC(p7P_BM,k-1))));

      MMX_FS(i,k,p7G_C1) = IVX(ivx_1,k) + p7P_MSC_CODON(gm_fs, k, c1);

      MMX_FS(i,k,p7G_C2) = IVX(ivx_2,k) + p7P_MSC_CODON(gm_fs, k, c2);

      MMX_FS(i,k,p7G_C3) = IVX(ivx_3,k) + p7P_MSC_CODON(gm_fs, k, c3);

      MMX_FS(i,k,p7G_C4) = IVX(ivx_4,k) + p7P_MSC_CODON(gm_fs, k, c4);

      MMX_FS(i,k,p7G_C5) = IVX(ivx_5,k) + p7P_MSC_CODON(gm_fs, k, c5);

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

    MMX_FS(i,M,p7G_C1) = IVX(ivx_1,M) + p7P_MSC_CODON(gm_fs, M, c1);

    MMX_FS(i,M,p7G_C2) = IVX(ivx_2,M) + p7P_MSC_CODON(gm_fs, M, c2);

    MMX_FS(i,M,p7G_C3) = IVX(ivx_3,M) + p7P_MSC_CODON(gm_fs, M, c3);

    MMX_FS(i,M,p7G_C4) = IVX(ivx_4,M) + p7P_MSC_CODON(gm_fs, M, c4);

    MMX_FS(i,M,p7G_C5) = IVX(ivx_5,M) + p7P_MSC_CODON(gm_fs, M, c5);

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
    XMX_FS(i,p7G_J) = ESL_MAX(XMX_FS(i-3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP],
                              XMX_FS(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_LOOP]);

    XMX_FS(i,p7G_C) = ESL_MAX(XMX_FS(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                              XMX_FS(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);
    XMX_FS(i,p7G_N) =         XMX_FS(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];
    XMX_FS(i,p7G_B) = ESL_MAX(XMX_FS(i,p7G_N) + gm_fs->xsc[p7P_N][p7P_MOVE],
                              XMX_FS(i,p7G_J) + gm_fs->xsc[p7P_J][p7P_MOVE]);
  }
  
  /* T state (not stored) */
  if (opt_sc != NULL) *opt_sc = ESL_MAX( XMX_FS(L,p7G_C),
                                ESL_MAX( XMX_FS(L-1,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                                         XMX_FS(L-2,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP])) + gm_fs->xsc[p7P_C][p7P_MOVE];

  gx->M = gm_fs->M;
  gx->L = L;

  return eslOK;

}



/* Function: p7_VTrace_Frameshift()
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
p7_VTrace_Frameshift(const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs, const P7_GMX *gx, P7_TRACE *tr)
{
  int          i   = L;     /* position in seq (1..L)         */
  int          k   = 0;     /* position in model (1..M)       */
  int          M   = gm_fs->M;
  int          c   = 0;
  int          prev_c;
  float      **dp  = gx->dp;    /* so {MDI}MX() macros work       */
  float       *xmx = gx->xmx;   /* so XMX() macro works           */
  float        tol = 1e-5;  /* floating point "equality" test */
  float const *tsc = gm_fs->tsc;
  int     sprv, scur;       /* previous, current state in trace */
  float   path[4];
  int     state[4] = {p7T_M, p7T_I, p7T_D, p7T_B};
  float   match_codon[5];
  int     status;

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
      else if (esl_FCompare_old(XMX_FS(i, p7G_C), XMX_FS(i-3, p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP], tol) == eslOK)  scur = p7T_C;
      else if (esl_FCompare_old(XMX_FS(i, p7G_C), XMX_FS(i,   p7G_E) + gm_fs->xsc[p7P_E][p7P_MOVE], tol) == eslOK)  scur = p7T_E;
      else ESL_EXCEPTION(eslFAIL, "C at i=%d couldn't be traced", i);
 
      break;

    case p7T_E:     /* E connects from any M state. k set here */
      if (XMX_FS(i, p7G_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);

      if (p7_fs_profile_IsLocal(gm_fs))
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

      if      (esl_FCompare_old(XMX_FS(i,p7G_B), XMX_FS(i, p7G_N) + gm_fs->xsc[p7P_N][p7P_MOVE], tol) == eslOK) scur = p7T_N;
      else if (esl_FCompare_old(XMX_FS(i,p7G_B), XMX_FS(i, p7G_J) + gm_fs->xsc[p7P_J][p7P_MOVE], tol) == eslOK) scur = p7T_J;
      else  ESL_EXCEPTION(eslFAIL, "B at i=%d couldn't be traced", i);
      break;

    case p7T_J:         /* J connects from E(i) or J(i-1) */
      if (XMX_FS(i,p7G_J) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible J reached at i=%d", i);
   
      if      (esl_FCompare_old(XMX_FS(i,p7G_J), XMX_FS(i-3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP], tol) == eslOK) scur = p7T_J;
      else if (esl_FCompare_old(XMX_FS(i,p7G_J), XMX_FS(i,  p7G_E) + gm_fs->xsc[p7P_E][p7P_LOOP], tol) == eslOK) scur = p7T_E;
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

  tr->M = gm_fs->M;
  tr->L = L;

  return p7_trace_fs_Reverse(tr);
}





