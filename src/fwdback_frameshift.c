/* Forward/Backward Frameshift algorithms;
 * 
 * Contents:
 *   1. Forward, Backward, implementations.  
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gencode.h"

#include "hmmer.h"

#define IVX3(i,k) (ivx[((k)*p7P_3CODONS) + (i)])
#define IVX5(i,k) (ivx[((k)*p7P_5CODONS) + (i)])
/*****************************************************************
 * 1. Forward, Backward, implementations.
 *****************************************************************/

/* Function:  p7_Forward_Frameshift()
 * Synopsis:  The Frameshift Aware Forward algorithm.
 *
 * Purpose:   The Forward dynamic programming algorithm for frameshift
 *            aware translated comparison between a dna sequence and a
 *            framshift-aware codon HMM.  
 *
 *            Given a digital sequence <dsq> of length <L>, a profile
 *            <gm_fs>, and DP matrix <gx> allocated with <L> cells for 
 *            all 5 special states <N,B,E,J,C>, and <L> by <gm_fs-M> 
 *            cells for all core model states <M,I,D> and all 5 codon
 *            lengths; calculate the probability of the sequence
 *            given the model using the Forward algorithm; use the 
 *            intermediate value matrix to store partial calculations; 
 *            return the Forward matrix in <gx>, and the Forward score 
 *            in <ret_sc>.
 *            
 *            The Forward score is in lod score form.  To convert to a
 *            bitscore, the caller needs to subtract a null model lod
 *            score, then convert to bits.
 *           
 *            Caller must have initialized the log-sum calculation
 *            with a call to <p7_FLogsumInit()>.
 *
 * Args:      dsq    - nucleotide sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gcode  - genetic code table 
 *            gm_fs  - a frameshift-aware codon profile. 
 *            gx     - DP matrix with room for an MxL alignment
 *            iv     - intermediate value matrix 
 *            opt_sc - optRETURN: Forward lod score in nats
 *           
 * Return:    <eslOK> on success.
 */
int
p7_Forward_Frameshift(const ESL_DSQ *dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx, P7_IVX *iv, float *opt_sc)
{ 

  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;           
  float       *ivx  = iv->ivx;
  int          M    = gm_fs->M;
  int          i, k;  
  int          c1, c2, c3, c4, c5;
  int          t, u, v, w, x;  
  int          ivx_1, ivx_2, ivx_3, ivx_4, ivx_5;
  float        esc  = p7_fs_profile_IsLocal(gm_fs) ? 0 : -eslINFINITY;

  for (i = 0; i < p7P_5CODONS; i++) {
    for(k = 0; k <= M; k++)
      IVX5(i,k) = -eslINFINITY;
  }

  /* Initialization for row 0 */
  XMX_FS(0,p7G_N) = 0.; //* S->N, p=1            */
  XMX_FS(0,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE];                   /* S->N->B, no N-tail   */
  XMX_FS(0,p7G_E) = XMX_FS(0,p7G_J) = XMX_FS(0,p7G_C) = -eslINFINITY;
  for (k = 0; k <= M; k++) {
    MMX_FS(0,k,p7G_C1) = MMX_FS(0,k,p7G_C2) = MMX_FS(0,k,p7G_C3) = MMX_FS(0,k,p7G_C4) = MMX_FS(0,k,p7G_C5) = -eslINFINITY;
    MMX_FS(0,k,p7G_C0) = IMX_FS(0,k) = DMX_FS(0,k) = -eslINFINITY;
  }

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
    IVX5(1,k) = XMX_FS(0,p7G_B) + TSC(p7P_BM,k-1);  

    MMX_FS(1,k,p7G_C1) = IVX5(1,k) + p7P_MSC_CODON(gm_fs, k, c1);
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
    IVX5(2,k) = XMX_FS(1,p7G_B) + TSC(p7P_BM,k-1);  
	MMX_FS(2,k,p7G_C1) = IVX5(2,k) + p7P_MSC_CODON(gm_fs, k, c1);
    MMX_FS(2,k,p7G_C2) = IVX5(1,k) + p7P_MSC_CODON(gm_fs, k, c2);
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
      IVX5(ivx_1,k) = p7_FLogsum(MMX_FS(i-1,k-1,p7G_C0) + TSC(p7P_MM,k-1),
                      p7_FLogsum(IMX_FS(i-1,k-1)        + TSC(p7P_IM,k-1),
                      p7_FLogsum(DMX_FS(i-1,k-1)        + TSC(p7P_DM,k-1),
                                 XMX_FS(i-1,p7G_B)      + TSC(p7P_BM,k-1))));

      MMX_FS(i,k,p7G_C1) = IVX5(ivx_1,k) + p7P_MSC_CODON(gm_fs, k, c1);
      MMX_FS(i,k,p7G_C2) = IVX5(ivx_2,k) + p7P_MSC_CODON(gm_fs, k, c2);
      MMX_FS(i,k,p7G_C3) = IVX5(ivx_3,k) + p7P_MSC_CODON(gm_fs, k, c3);
      MMX_FS(i,k,p7G_C4) = -eslINFINITY;
      MMX_FS(i,k,p7G_C5) = -eslINFINITY;
      if( i == 4 )
        MMX_FS(i,k,p7G_C4) = IVX5(ivx_4,k) + p7P_MSC_CODON(gm_fs, k, c4);

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

    MMX_FS(i,M,p7G_C1) = IVX5(ivx_1,M) + p7P_MSC_CODON(gm_fs, M, c1);
    MMX_FS(i,M,p7G_C2) = IVX5(ivx_2,M) + p7P_MSC_CODON(gm_fs, M, c2);
    MMX_FS(i,M,p7G_C3) = IVX5(ivx_3,M) + p7P_MSC_CODON(gm_fs, M, c3);
    MMX_FS(i,M,p7G_C4) = -eslINFINITY;
    MMX_FS(i,M,p7G_C5) = -eslINFINITY;
    if( i == 4 )
      MMX_FS(i,M,p7G_C4) = IVX5(ivx_4,M) + p7P_MSC_CODON(gm_fs, M, c4);

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

    XMX_FS(i,p7G_J) = p7_FLogsum(XMX_FS(i-3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP],
                              XMX_FS(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_LOOP]);

    XMX_FS(i,p7G_C) = p7_FLogsum(XMX_FS(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                                 XMX_FS(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);

    XMX_FS(i,p7G_N) =            XMX_FS(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP]; 

    XMX_FS(i,p7G_B) = p7_FLogsum(XMX_FS(i,p7G_N) + gm_fs->xsc[p7P_N][p7P_MOVE],
                                 XMX_FS(i,p7G_J) + gm_fs->xsc[p7P_J][p7P_MOVE]);

  }

  /* Main Recusion. Done as a pull */
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
      
      IVX5(ivx_1,k) = p7_FLogsum(MMX_FS(i-1,k-1,p7G_C0)   + TSC(p7P_MM,k-1), 
                      p7_FLogsum(IMX_FS(i-1,k-1)          + TSC(p7P_IM,k-1),
                      p7_FLogsum(DMX_FS(i-1,k-1)          + TSC(p7P_DM,k-1),
                                 XMX_FS(i-1,p7G_B)        + TSC(p7P_BM,k-1))));

      MMX_FS(i,k,p7G_C1) = IVX5(ivx_1,k) + p7P_MSC_CODON(gm_fs, k, c1); 

      MMX_FS(i,k,p7G_C2) = IVX5(ivx_2,k) + p7P_MSC_CODON(gm_fs, k, c2);
    
      MMX_FS(i,k,p7G_C3) = IVX5(ivx_3,k) + p7P_MSC_CODON(gm_fs, k, c3);
      
      MMX_FS(i,k,p7G_C4) = IVX5(ivx_4,k) + p7P_MSC_CODON(gm_fs, k, c4);

      MMX_FS(i,k,p7G_C5) = IVX5(ivx_5,k) + p7P_MSC_CODON(gm_fs, k, c5); 

      MMX_FS(i,k,p7G_C0) =  p7_FLogsum(p7_FLogsum(MMX_FS(i,k,p7G_C1), 
                            p7_FLogsum(MMX_FS(i,k,p7G_C2), MMX_FS(i,k,p7G_C3))),
                            p7_FLogsum(MMX_FS(i,k,p7G_C4), MMX_FS(i,k,p7G_C5)));

      /* insert state */
      IMX_FS(i,k) = p7_FLogsum(MMX_FS(i-3,k,p7G_C0) + TSC(p7P_MI,k),
                               IMX_FS(i-3,k)        + TSC(p7P_II,k));
    
      /* delete state */
      DMX_FS(i,k) = p7_FLogsum(MMX_FS(i,k-1,p7G_C0) + TSC(p7P_MD,k-1),
                               DMX_FS(i,k-1)        + TSC(p7P_DD,k-1));

      /* E state update */
      XMX_FS(i,p7G_E) = p7_FLogsum(MMX_FS(i,k,p7G_C0) + esc,
                        p7_FLogsum(DMX_FS(i,k)        + esc,  
                                   XMX_FS(i,p7G_E)));
    }

    /* unrolled match state M_M */
    IVX5(ivx_1,M) = p7_FLogsum(MMX_FS(i-1,M-1,p7G_C0)   + TSC(p7P_MM,M-1), 
                    p7_FLogsum(IMX_FS(i-1,M-1)          + TSC(p7P_IM,M-1),
                    p7_FLogsum(DMX_FS(i-1,M-1)          + TSC(p7P_DM,M-1),
                               XMX_FS(i-1,p7G_B)        + TSC(p7P_BM,M-1))));

    MMX_FS(i,M,p7G_C1) = IVX5(ivx_1,M) + p7P_MSC_CODON(gm_fs, M, c1);

    MMX_FS(i,M,p7G_C2) = IVX5(ivx_2,M) + p7P_MSC_CODON(gm_fs, M, c2);
    
    MMX_FS(i,M,p7G_C3) = IVX5(ivx_3,M) + p7P_MSC_CODON(gm_fs, M, c3); 

    MMX_FS(i,M,p7G_C4) = IVX5(ivx_4,M) + p7P_MSC_CODON(gm_fs, M, c4); 

    MMX_FS(i,M,p7G_C5) = IVX5(ivx_5,M) + p7P_MSC_CODON(gm_fs, M, c5);

    MMX_FS(i,M,p7G_C0) =  p7_FLogsum(p7_FLogsum(MMX_FS(i,M,p7G_C1), 
                          p7_FLogsum(MMX_FS(i,M,p7G_C2), MMX_FS(i,M,p7G_C3))),
                          p7_FLogsum(MMX_FS(i,M,p7G_C4), MMX_FS(i,M,p7G_C5)));

    IMX_FS(i,M) = -eslINFINITY;

    /* unrolled delete state D_M */
    DMX_FS(i,M) = p7_FLogsum(MMX_FS(i,M-1,p7G_C0) + TSC(p7P_MD,M-1),
                             DMX_FS(i,M-1) + TSC(p7P_DD,M-1));

    /* unrolled E state update */
    XMX_FS(i,p7G_E) = p7_FLogsum(p7_FLogsum(MMX_FS(i,M,p7G_C0),
                                            DMX_FS(i,M)),
                                            XMX_FS(i,p7G_E));

    /* J, C and N states */
    XMX_FS(i,p7G_J) = p7_FLogsum(XMX_FS(i-3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP],
                                 XMX_FS(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_LOOP]);
    XMX_FS(i,p7G_C) = p7_FLogsum(XMX_FS(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                                 XMX_FS(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);
    XMX_FS(i,p7G_N) =            XMX_FS(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];
    XMX_FS(i,p7G_B) = p7_FLogsum(XMX_FS(i,p7G_N) + gm_fs->xsc[p7P_N][p7P_MOVE],
                                 XMX_FS(i,p7G_J) + gm_fs->xsc[p7P_J][p7P_MOVE]);
	
  }

  if (opt_sc != NULL) *opt_sc = p7_FLogsum( XMX_FS(L,p7G_C),
                                p7_FLogsum( XMX_FS(L-1,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                                            XMX_FS(L-2,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP])) + 
                                            gm_fs->xsc[p7P_C][p7P_MOVE];
 
  gx->M = M;
  gx->L = L;
 
  return eslOK;

}


/* Function:  p7_ForwardParser_Frameshift_3Codons()
 * Synopsis:  The frameshift-aware Forward algorithm using 3 codon lengths - low memory.
 *
 * Purpose:   The Forward dynamic programming algorithm for frameshift
 *            aware translated comparison between a dna sequence and an
 *            frameshift aware codon HMM, using three codon lengths,
 *            with a minimal sized DP matrix.
 *
 *            Given a digital sequence <dsq> of length <L>, a profile
 *            <gm_fs>, and DP matrix <gx> allocated with <L> cells for
 *            all 5 special states <N,B,E,J,C>, and at least 2 by
 *            <gm_fs->M> cells for all core model states <M,I,D>;
 *            calculate the probability of the sequence given the model
 *            and 3 permissible codon lengths using the Forward algorithm -
 *            reusing the rows in the core matrix; use the intermdiate
 *            value matrix <iv> to store partial calculations; return the
 *            Forward matrix in <gx>, and the Forward score in <ret_sc>.
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
 *            gcode  - genetic code table
 *            gm_fs  - framshift-aware codon profile.
 *            gx     - DP matrix
 *            iv     - intermedite value matrix
 *            opt_sc - optRETURN: Forward lod score in nats
 *
 * Return:    <eslOK> on success.
 */



/* Function:  p7_ForwardParser_Frameshift_3Codons()
 * Synopsis:  The frameshift-aware Forward algorithm using 3 codon lengths - low memory.
 *
 * Purpose:   The Forward dynamic programming algorithm for frameshift
 *            aware translated comparison between a dna sequence and an
 *            frameshift aware codon HMM, using three codon lengths, 
 *            with a minimal sized DP matrix. 
 *
 *            Given a digital sequence <dsq> of length <L>, a profile
 *            <gm_fs>, and DP matrix <gx> allocated with <L> cells for
 *            the special states <N,B,E,J,C>, and <PARSER_ROWS_FWD> by 
 *            <gm_fs->M> cells for all core model states <M,I,D>; 
 *            calculate the probability of the sequence given the model 
 *            and 3 permissible codon lengths using the Forward algorithm - 
 *            reusing the rows in the core matrix; use the intermdiate 
 *            value matrix <iv> to store partial calculations; return the 
 *            Forward matrix in <gx>, and the Forward score in <ret_sc>.
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
 *            gcode  - genetic code table
 *            gm_fs  - framshift-aware codon profile. 
 *            gx     - DP matrix 
 *            iv     - intermedite value matrix 
 *            opt_sc - optRETURN: Forward lod score in nats
 *           
 * Return:    <eslOK> on success.
 */
int
p7_ForwardParser_Frameshift_3Codons(const ESL_DSQ *dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx, P7_IVX *iv, float *opt_sc)
{ 

  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;           
  float       *ivx  = iv->ivx;
  int          M    = gm_fs->M;
  int          i, k;
  int          c2, c3, c4;  
  int          u, v, w, x;
  float        esc  = p7_fs_profile_IsLocal(gm_fs) ? 0 : -eslINFINITY;
  int curr, prev2, prev3;
  int ivx_2, ivx_3, ivx_4;

  for (i = 0; i < p7P_3CODONS; i++) {
    for(k = 0; k <= M; k++) { 
      IVX3(i,k) = -eslINFINITY;
    }
  }

  /* Initialization of rows 0-1 */
  for(i = 0; i < 2; i++ ) {
    XMX(i,p7G_N) = 0.;                          /* S->N, p=1          */
    XMX(i,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE]; /* S->N->B, no N-tail */
    XMX(i,p7G_E) = XMX(i,p7G_J) = XMX(i,p7G_C) = -eslINFINITY;
    for (k = 0; k <= M; k++) 
      MMX(i,k) = IMX(i,k) = DMX(i,k) = -eslINFINITY;
  }

  /* Initialization of row 2 */
  u = v = p7P_MAXCODONS;
  
  if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[1])) w = dsq[1];
  else                                            w = p7P_MAXCODONS;

  if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[2])) x = dsq[2];
  else                                            x = p7P_MAXCODONS;

  /* Length 2 quasi-codon */
  c2 = p7P_CODON2(w, x);
  c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);
 
  XMX(2,p7G_E) = -eslINFINITY;
  MMX(2,0) = IMX(2,0) = DMX(2,0) = -eslINFINITY;
  
  for (k = 1; k <= M; k++)
  {
    IVX3(2,k) = XMX(0,p7G_B) + TSC(p7P_BM,k-1);
    MMX(2,k)  = IVX3(2,k) + p7P_MSC_CODON(gm_fs, k, c2); 
    IMX(2,k) = -eslINFINITY;
    DMX(2,k) = p7_FLogsum(MMX(2,k-1)        + TSC(p7P_MD,k-1),
                          DMX(2,k-1)        + TSC(p7P_DD,k-1));
    XMX(2,p7G_E) = p7_FLogsum(MMX(2,k) + esc,
                   p7_FLogsum(DMX(2,k) + esc,
                              XMX(2,p7G_E)));
  }

  XMX(2,p7G_N) = 0.; /* S->N, p=1            */
  XMX(2,p7G_J) = XMX(2,p7G_E) + gm_fs->xsc[p7P_E][p7P_LOOP];
  XMX(2,p7G_C) = XMX(2,p7G_E) + gm_fs->xsc[p7P_E][p7P_MOVE]; 
  XMX(2,p7G_B) = p7_FLogsum(XMX(2,p7G_N) + gm_fs->xsc[p7P_N][p7P_MOVE],
                            XMX(2,p7G_J) + gm_fs->xsc[p7P_J][p7P_MOVE]);

  /* Main Recusion */
  for(i = 3; i <= L; i++) {

    u = v;
    v = w;
    w = x;
	
    /* if new nucleotide is not A,C,G, or T set it to placeholder value */
    if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[i])) x = dsq[i];
    else                                            x = p7P_MAXCODONS;

	/* Length 2 quasi-codon */
	c2 = p7P_CODON2(w, x);
	c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

	/* Length 3 codon */
    c3 = p7P_CODON3(v, w, x);
	c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

	/* Length 4 quasi-codon */
	c4 = p7P_CODON4(u, v, w, x);
	c4 = p7P_MINIDX(c4, p7P_DEGEN_QC1);

	curr  = i     % PARSER_ROWS_FWD;
    prev2 = (i-2) % PARSER_ROWS_FWD;
    prev3 = (i-3) % PARSER_ROWS_FWD;

	ivx_2 = i     % p7P_3CODONS;
	ivx_3 = (i-1) % p7P_3CODONS;
	ivx_4 = (i-2) % p7P_3CODONS;

	MMX(curr,0) = IMX(curr,0) = DMX(curr,0) = -eslINFINITY;
    
	XMX(i,p7G_E) = -eslINFINITY;

    for (k = 1; k < M; k++)
    {
    /* For every row the 1 nuc codon translations must be reclacuated (i-1).
     * The other codon transtions can be reused with i-1 becoming 1-2,
     * i-2 becoming i-3 and so on */
      IVX3(ivx_2,k) = p7_FLogsum(MMX(prev2,k-1) + TSC(p7P_MM,k-1),  
                      p7_FLogsum(IMX(prev2,k-1) + TSC(p7P_IM,k-1),
                      p7_FLogsum(DMX(prev2,k-1) + TSC(p7P_DM,k-1),
                                 XMX(i-2,p7G_B) + TSC(p7P_BM,k-1))));	
           
      MMX(curr,k) =                         IVX3(ivx_2,k) + p7P_MSC_CODON(gm_fs, k, c2);
      MMX(curr,k) = p7_FLogsum(MMX(curr,k), IVX3(ivx_3,k) + p7P_MSC_CODON(gm_fs, k, c3));
      MMX(curr,k) = p7_FLogsum(MMX(curr,k), IVX3(ivx_4,k) + p7P_MSC_CODON(gm_fs, k, c4));

	  IMX(curr,k) = p7_FLogsum(MMX(prev3,k)  + TSC(p7P_MI,k),
                               IMX(prev3,k)  + TSC(p7P_II,k));

      DMX(curr,k) = p7_FLogsum(MMX(curr,k-1) + TSC(p7P_MD,k-1),
                               DMX(curr,k-1) + TSC(p7P_DD,k-1));

      XMX(i,p7G_E) = p7_FLogsum(MMX(curr,k) + esc,
                     p7_FLogsum(DMX(curr,k) + esc,
                                XMX(i,p7G_E)));
	}

	IVX3(ivx_2,M) = p7_FLogsum(MMX(prev2,M-1) + TSC(p7P_MM,M-1),  
                    p7_FLogsum(IMX(prev2,M-1) + TSC(p7P_IM,M-1),
                    p7_FLogsum(DMX(prev2,M-1) + TSC(p7P_DM,M-1),
                               XMX(i-2,p7G_B) + TSC(p7P_BM,M-1))));
    
	MMX(curr,M) =                         IVX3(ivx_2,M) + p7P_MSC_CODON(gm_fs, M, c2);
    MMX(curr,M) = p7_FLogsum(MMX(curr,M), IVX3(ivx_3,M) + p7P_MSC_CODON(gm_fs, M, c3));
    MMX(curr,M) = p7_FLogsum(MMX(curr,M), IVX3(ivx_4,M) + p7P_MSC_CODON(gm_fs, M, c4));

	IMX(curr,M) = -eslINFINITY;

    DMX(curr,M) = p7_FLogsum(MMX(curr,M-1) + TSC(p7P_MD,M-1),
                             DMX(curr,M-1)        + TSC(p7P_DD,M-1));

    XMX(i,p7G_E) = p7_FLogsum(MMX(curr,M),
                   p7_FLogsum(DMX(curr,M),
                              XMX(i,p7G_E)));

    XMX(i,p7G_N) =            XMX(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];
    
    XMX(i,p7G_J) = p7_FLogsum(XMX(i-3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP],
                              XMX(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_LOOP]);

    XMX(i,p7G_C) = p7_FLogsum(XMX(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                              XMX(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);

    XMX(i,p7G_B) = p7_FLogsum(XMX(i,p7G_N) + gm_fs->xsc[p7P_N][p7P_MOVE],
                              XMX(i,p7G_J) + gm_fs->xsc[p7P_J][p7P_MOVE]);
  
  }
   
  if (opt_sc != NULL) *opt_sc = p7_FLogsum( XMX(L,p7G_C),
                                p7_FLogsum( XMX(L-1,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                                            XMX(L-2,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP])) +
		                                    gm_fs->xsc[p7P_C][p7P_MOVE];
  gx->M = M;
  gx->L = L;
 
  return eslOK;

}


/* Function:  p7_ForwardParser_Frameshift_5Codons()
 * Synopsis:  The frameshift-aware Forward algorithm using 5 codon lengths - low memory.
 *
 * Purpose:   The Forward dynamic programming algorithm for frameshift
 *            aware translated comparison between a dna sequence and an
 *            amino acid HMM. 
 *
 *            Given a digital sequence <dsq> of length <L>, a profile
 *            <gmi_fs>, and DP matrix <gx> allocated with <L> cells for
 *            the special states <N,B,E,J,C>, and <PARSER_ROWS_FWD> by 
 *            <gm_fs->M> cells for all core model states <M,I,D>; 
 *            calculate th probability of the sequence given the model 
 *            and 5 permissible codon lengths using the Forward algorithm 
 *            - reusing the rows in the core matrix; use the intermdiate
 *            value matrix <iv> to store partial calculations; return the 
 *            Forward  matrix in <gx>, and the Forward score in <ret_sc>.
 *           
 *            The Forward score is in lod score form.  To convert to a
 *            bitscore, the caller needs to subtract a null model lod
 *            score, then convert to bits.
 *           
 *            Caller must have initialized the log-sum calculation
 *            with a call to <p7_FLogsumInit()>.
 *
 * Args:      dsq    - nucleotide sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gcode  - genetic code table
 *            gm_fs  - framshift aware codon profile. 
 *            gx     - DP matrix with room for an MxL alignment
 *            iv     - intermediate value martix
 *            opt_sc - optRETURN: Forward lod score in nats
 *           
 * Return:    <eslOK> on success.
 */
int
p7_ForwardParser_Frameshift_5Codons(const ESL_DSQ *dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx, P7_IVX *iv, float *opt_sc)
{ 

  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;           
  float       *ivx  = iv->ivx;
  int          M    = gm_fs->M;
  int          i, k;
  int          c1, c2, c3, c4, c5;  
  int          t, u, v, w, x;
  float        esc  = p7_fs_profile_IsLocal(gm_fs) ? 0 : -eslINFINITY;
  int curr, prev1, prev3;
  int ivx_1, ivx_2, ivx_3, ivx_4, ivx_5;

  for (i = 0; i < p7P_5CODONS; i++) {
    for(k = 0; k <= M; k++) {
      IVX5(i,k) = -eslINFINITY;
    }
  }

  /* Initialization for row 0 */
  XMX(0,p7G_N) = 0.; /* S->N, p=1            */
  XMX(0,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE];                   /* S->N->B, no N-tail   */
  XMX(0,p7G_E) = XMX(0,p7G_J) = XMX(0,p7G_C) = -eslINFINITY;
  for (k = 0; k <= M; k++) 
    MMX(0,k) = IMX(0,k) = DMX(0,k)        = -eslINFINITY;

  /* Initialization for row 1 */
  XMX(1,p7G_N) = 0.; 
  XMX(1,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE];                   
  XMX(1,p7G_E) = -eslINFINITY;
  MMX(1,0) = IMX(1,0) = DMX(1,0) = -eslINFINITY;
  
  if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[1])) x = dsq[1];
  else                                            x = p7P_MAXCODONS;
  
  c1 = p7P_CODON1(x);
  c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);
  for (k = 1; k < M; k++) {
    IVX5(1,k) = XMX(0,p7G_B) + TSC(p7P_BM,k-1);  
	MMX(1,k)  = IVX5(1,k) + p7P_MSC_CODON(gm_fs, k, c1);
	IMX(1,k)  = -eslINFINITY;
    DMX(1,k)  = p7_FLogsum(MMX(1,k-1) + TSC(p7P_MD,k-1),
                           DMX(1,k-1) + TSC(p7P_DD,k-1));

	XMX(1,p7G_E) = p7_FLogsum(MMX(1,k) + esc,
                   p7_FLogsum(DMX(1,k) + esc,
                              XMX(1,p7G_E)));
  }
 
  XMX(1,p7G_J) = XMX(1,p7G_E) + gm_fs->xsc[p7P_E][p7P_LOOP];
  XMX(1,p7G_C) = XMX(1,p7G_E) + gm_fs->xsc[p7P_E][p7P_MOVE];
   

  /* Initialization for row 2 */
  XMX(2,p7G_N) = 0.; 
  XMX(2,p7G_B) = gm_fs->xsc[p7P_N][p7P_MOVE];                   
  XMX(2,p7G_E) = -eslINFINITY;
  MMX(2,0) = IMX(2,0) = DMX(2,0) = -eslINFINITY;
  
  w = x;
  if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[2])) x = dsq[2];
  else                                            x = p7P_MAXCODONS;
  
  c1 = p7P_CODON1(x);
  c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

  c2 = p7P_CODON2(w, x);
  c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

  for (k = 1; k < M; k++) {
    IVX5(2,k) = XMX(1,p7G_B) + TSC(p7P_BM,k-1);  
	MMX(2,k)  = IVX5(2,k) + p7P_MSC_CODON(gm_fs, k, c1);
	MMX(2,k)  = p7_FLogsum(MMX(2,k), IVX5(1,k) + p7P_MSC_CODON(gm_fs, k, c2)); //IVX5(1,k) now holds the i-2 transtion
	IMX(2,k)  = -eslINFINITY;
    DMX(2,k)  = p7_FLogsum(MMX(2,k-1) + TSC(p7P_MD,k-1),
                           DMX(2,k-1) + TSC(p7P_DD,k-1));

	XMX(2,p7G_E) = p7_FLogsum(MMX(2,k) + esc,
                   p7_FLogsum(DMX(2,k) + esc,
                              XMX(2,p7G_E)));
  }
 
  XMX(2,p7G_J) = XMX(2,p7G_E) + gm_fs->xsc[p7P_E][p7P_LOOP];
  XMX(2,p7G_C) = XMX(2,p7G_E) + gm_fs->xsc[p7P_E][p7P_MOVE];
  
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

    curr  = i     % PARSER_ROWS_FWD;
    prev1 = (i-1) % PARSER_ROWS_FWD;
    prev3 = (i-3) % PARSER_ROWS_FWD;

    ivx_1 = i     % p7P_5CODONS;
    ivx_2 = (i-1) % p7P_5CODONS;
    ivx_3 = (i-2) % p7P_5CODONS;
    ivx_4 = (i-3) % p7P_5CODONS;

    MMX(curr,0) = IMX(curr,0) = DMX(curr,0) = -eslINFINITY;

    XMX(i,p7G_E) = -eslINFINITY;
   
    for (k = 1; k < M; k++)
    {
      IVX5(ivx_1,k) = p7_FLogsum(MMX(prev1,k-1)   + TSC(p7P_MM,k-1),
                      p7_FLogsum(IMX(prev1,k-1)   + TSC(p7P_IM,k-1),
                      p7_FLogsum(DMX(prev1,k-1)   + TSC(p7P_DM,k-1),
                                 XMX(i-1,p7G_B)   + TSC(p7P_BM,k-1))));

      MMX(curr,k) =                         IVX5(ivx_1,k) + p7P_MSC_CODON(gm_fs, k, c1);
      MMX(curr,k) = p7_FLogsum(MMX(curr,k), IVX5(ivx_2,k) + p7P_MSC_CODON(gm_fs, k, c2));
      MMX(curr,k) = p7_FLogsum(MMX(curr,k), IVX5(ivx_3,k) + p7P_MSC_CODON(gm_fs, k, c3));
      if( i == 4 )
        MMX(curr,k) = p7_FLogsum(MMX(curr,k), IVX5(ivx_4,k) + p7P_MSC_CODON(gm_fs, k, c4));
      
      IMX(curr,k) = p7_FLogsum(MMX(prev3,k)        + TSC(p7P_MI,k),
                               IMX(prev3,k)        + TSC(p7P_II,k));

      DMX(curr,k) = p7_FLogsum(MMX(curr,k-1)        + TSC(p7P_MD,k-1),
                               DMX(curr,k-1)        + TSC(p7P_DD,k-1));

      XMX(i,p7G_E) = p7_FLogsum(MMX(curr,k) + esc,
                     p7_FLogsum(DMX(curr,k) + esc, 
                                XMX(i,p7G_E)));

    }

    IVX5(ivx_1,M) = p7_FLogsum(MMX(prev1,M-1) + TSC(p7P_MM,M-1),
                    p7_FLogsum(IMX(prev1,M-1) + TSC(p7P_IM,M-1),
                    p7_FLogsum(DMX(prev1,M-1) + TSC(p7P_DM,M-1),
                               XMX(i-1,p7G_B) + TSC(p7P_BM,M-1))));

    MMX(curr,M) =                         IVX5(ivx_1,M) + p7P_MSC_CODON(gm_fs, M, c1);
    MMX(curr,M) = p7_FLogsum(MMX(curr,M), IVX5(ivx_2,M) + p7P_MSC_CODON(gm_fs, M, c2));
    MMX(curr,M) = p7_FLogsum(MMX(curr,M), IVX5(ivx_3,M) + p7P_MSC_CODON(gm_fs, M, c3));
    if( i == 4 )
      MMX(curr,M) = p7_FLogsum(MMX(curr,M), IVX5(ivx_4,M) + p7P_MSC_CODON(gm_fs, M, c4));

    IMX(curr,M) = -eslINFINITY;
    
    DMX(curr,M) = p7_FLogsum(MMX(curr,M-1) + TSC(p7P_MD,M-1),
                             DMX(curr,M-1)        + TSC(p7P_DD,M-1));

    XMX(i,p7G_E) = p7_FLogsum(MMX(curr,M),
                    p7_FLogsum(DMX(curr,M), 
                               XMX(i,p7G_E)));

    XMX(i,p7G_J) = p7_FLogsum(XMX(i-3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP],
                              XMX(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_LOOP]);

    XMX(i,p7G_C) = p7_FLogsum(XMX(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                              XMX(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);

    XMX(i,p7G_N) =            XMX(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP]; 

    XMX(i,p7G_B) = p7_FLogsum(XMX(i,p7G_N) + gm_fs->xsc[p7P_N][p7P_MOVE],
                              XMX(i,p7G_J) + gm_fs->xsc[p7P_J][p7P_MOVE]);

  }

  /* Main Recusion. Done as a pull. */
 
  for (i = 5; i <= L; i++) 
  {
    curr  = i     % PARSER_ROWS_FWD;
    prev1 = (i-1) % PARSER_ROWS_FWD;
    prev3 = (i-3) % PARSER_ROWS_FWD;

	ivx_1 = i     % p7P_5CODONS;
    ivx_2 = (i-1) % p7P_5CODONS;
    ivx_3 = (i-2) % p7P_5CODONS;
    ivx_4 = (i-3) % p7P_5CODONS;
    ivx_5 = (i-5) % p7P_5CODONS;

    MMX(curr,0) =  IMX(curr,0)        = DMX(curr,0)        = -eslINFINITY;
 
    XMX(i, p7G_E) = -eslINFINITY;

    /* Reasign nucluotide to correct temporary holders for use in emissions array */ 
    t = u;
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

    c5 = p7P_CODON5(t, u, v, w, x);
    c5 = p7P_MINIDX(c5, p7P_DEGEN_QC2);

    for (k = 1; k < M; k++)
    {  

      IVX5(ivx_1,k) = p7_FLogsum(MMX(prev1,k-1) + TSC(p7P_MM,k-1), 
                      p7_FLogsum(IMX(prev1,k-1) + TSC(p7P_IM,k-1),
                      p7_FLogsum(DMX(prev1,k-1) + TSC(p7P_DM,k-1),
                                 XMX(i-1,p7G_B) + TSC(p7P_BM,k-1))));

      MMX(curr,k) =                          IVX5(ivx_1,k) + p7P_MSC_CODON(gm_fs, k, c1); 

      MMX(curr,k) =  p7_FLogsum(MMX(curr,k), IVX5(ivx_2,k) + p7P_MSC_CODON(gm_fs, k, c2));
    
      MMX(curr,k) =  p7_FLogsum(MMX(curr,k), IVX5(ivx_3,k) + p7P_MSC_CODON(gm_fs, k, c3));

      MMX(curr,k) =  p7_FLogsum(MMX(curr,k), IVX5(ivx_4,k) + p7P_MSC_CODON(gm_fs, k, c4));

      MMX(curr,k) =  p7_FLogsum(MMX(curr,k), IVX5(ivx_5,k) + p7P_MSC_CODON(gm_fs, k, c5)); 

      /* insert state */
      IMX(curr,k) = p7_FLogsum(MMX(prev3,k) + TSC(p7P_MI,k),
                              IMX(prev3,k)  + TSC(p7P_II,k));
    
      /* delete state */
      DMX(curr,k) = p7_FLogsum(MMX(curr,k-1) + TSC(p7P_MD,k-1),
                               DMX(curr,k-1) + TSC(p7P_DD,k-1));

      /* E state update */
      XMX(i,p7G_E) = p7_FLogsum(MMX(curr,k) + esc,
                     p7_FLogsum(DMX(curr,k) + esc,  
                                XMX(i,p7G_E)));
    }

    /* unrolled match state M_M */
    IVX5(ivx_1,M) = p7_FLogsum(MMX(prev1,M-1) + TSC(p7P_MM,M-1), 
                    p7_FLogsum(IMX(prev1,M-1) + TSC(p7P_IM,M-1),
                    p7_FLogsum(DMX(prev1,M-1) + TSC(p7P_DM,M-1),
                               XMX(i-1,p7G_B) + TSC(p7P_BM,M-1))));

    MMX(curr,M) =                         IVX5(ivx_1,M) + p7P_MSC_CODON(gm_fs, M, c1);

    MMX(curr,M) = p7_FLogsum(MMX(curr,M), IVX5(ivx_2,M) + p7P_MSC_CODON(gm_fs, M, c2));
    
    MMX(curr,M) = p7_FLogsum(MMX(curr,M), IVX5(ivx_3,M) + p7P_MSC_CODON(gm_fs, M, c3)); 

    MMX(curr,M) = p7_FLogsum(MMX(curr,M), IVX5(ivx_4,M) + p7P_MSC_CODON(gm_fs, M, c4)); 

    MMX(curr,M) = p7_FLogsum(MMX(curr,M), IVX5(ivx_5,M) + p7P_MSC_CODON(gm_fs, M, c5));

    IMX(curr,M) = -eslINFINITY;

    /* unrolled delete state D_M */
    DMX(curr,M) = p7_FLogsum(MMX(curr,M-1) + TSC(p7P_MD,M-1),
                             DMX(curr,M-1) + TSC(p7P_DD,M-1));

    /* unrolled E state update */
    XMX(i,p7G_E) = p7_FLogsum(MMX(curr,M),
                   p7_FLogsum(DMX(curr,M),
                              XMX(i,p7G_E)));

    /* J, C and N states */
    XMX(i,p7G_J) = p7_FLogsum(XMX(i-3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP],
                              XMX(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_LOOP]);

    XMX(i,p7G_C) = p7_FLogsum(XMX(i-3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                              XMX(i,p7G_E)   + gm_fs->xsc[p7P_E][p7P_MOVE]);

    XMX(i,p7G_N) =            XMX(i-3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP];

    XMX(i,p7G_B) = p7_FLogsum(XMX(i,p7G_N) + gm_fs->xsc[p7P_N][p7P_MOVE],
                              XMX(i,p7G_J) + gm_fs->xsc[p7P_J][p7P_MOVE]);
  }


  if (opt_sc != NULL) *opt_sc = p7_FLogsum( XMX(L,p7G_C),
                                p7_FLogsum( XMX(L-1,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP],
                                            XMX(L-2,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP])) 
		                                   +  gm_fs->xsc[p7P_C][p7P_MOVE];
  gx->M = M;
  gx->L = L;
 
  return eslOK;
}

/* Function:  p7_Backward_Frameshift()
 * Synopsis:  The Backward algorithm.
 *
 * Purpose:   The Frameshift Aware Backward dynamic programming algorithm.
 * 
 *            Given a digital sequence <dsq> of length <L>, a profile
 *            <gm_fs>, and DP matrix <gx> allocated for at least <gm->M>
 *            by <L> cells; calculate the probability of the sequence
 *            given the model using the Backward algorithm;  use the
 *            intermediate value matrix to store partial calculations;
 *            return the Backward matrix in <gx>, and the Backward score 
 *            in <ret_sc>.
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
 *            gcode  - genetic code table
 *            gm_fs  - framshift aware codon profile.
 *            gx     - DP matrix with room for an MxL alignment
 *            iv     - intermediate value martix
 *            opt_sc - optRETURN: Backward lod score in nats
 *           
 * Return:    <eslOK> on success.
 */
int
p7_Backward_Frameshift(const ESL_DSQ *dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx, P7_IVX *iv, float *opt_sc)
{

  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;           
  float       *ivx   = iv->ivx;
  int          M    = gm_fs->M;
  int          i, k;;  
  int          c1, c2, c3, c4, c5;
  float        esc  = p7_fs_profile_IsLocal(gm_fs) ? 0 : -eslINFINITY;
  int          t, u, v, w, x;

  for(k = 0; k <= M; k++)
    ivx[k] = -eslINFINITY;

  /* Initialization of row L  */
  XMX(L,p7G_J) = XMX(L,p7G_B) = XMX(L,p7G_N) = -eslINFINITY; /* need to enter and exit model */
  XMX(L,p7G_C) = gm_fs->xsc[p7P_C][p7P_MOVE];                   /* C<-T         */
  XMX(L,p7G_E) = XMX(L,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE];    /* E<-C, no tail */
  MMX(L,M)     = DMX(L,M) = XMX(L,p7G_E);                    /* {MD}_M <- E (prob 1.0) */
  IMX(L,M)     = -eslINFINITY;                               /* no I_M state        */

  /* Initialization of core model for row L */
  for (k = M-1; k >= 1; k--) 
  {
    /*L comes form E & D state only */
    MMX(L,k) = p7_FLogsum( XMX(L,p7G_E) + esc,
                           DMX(L, k+1)  + TSC(p7P_MD,k));

    DMX(L,k) = p7_FLogsum( XMX(L,p7G_E) + esc,
                           DMX(L, k+1)  + TSC(p7P_DD,k));

    /* no I state at L to L-2 */
    IMX(L,k) = -eslINFINITY;
  }
  MMX(L,0) = IMX(L,0) = DMX(L,0)  = -eslINFINITY;

  /* Initialization of row L - 1  */
  if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[L])) x = dsq[L];
  else                                            x = p7P_MAXCODONS;

  c1 = p7P_CODON1(x);
  c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

  ivx[1]       = MMX(L,1) + p7P_MSC_CODON(gm_fs, 1, c1);
  XMX(L-1,p7G_B) = ivx[1]   + TSC(p7P_BM,0);
  for (k = 2; k <= M; k++)
  {
    ivx[k] =  MMX(L,k) + p7P_MSC_CODON(gm_fs, k, c1);

    XMX(L-1,p7G_B) = p7_FLogsum( XMX(L-1,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
  }

  XMX(L-1,p7G_J) = XMX(L-1,p7G_B) + gm_fs->xsc[p7P_J][p7P_MOVE];
  XMX(L-1,p7G_N) = XMX(L-1,p7G_B) + gm_fs->xsc[p7P_N][p7P_MOVE];
  XMX(L-1,p7G_C) = gm_fs->xsc[p7P_C][p7P_LOOP] + gm_fs->xsc[p7P_C][p7P_MOVE];
  XMX(L-1,p7G_E) = p7_FLogsum( XMX(L-1,p7G_J) + gm_fs->xsc[p7P_E][p7P_LOOP],
				               XMX(L-1,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE]);

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

  /* Initialization of row L - 2  */
  w = x;

  if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[L-1])) x = dsq[L-1];
  else                                            x = p7P_MAXCODONS;

  c1 = p7P_CODON1(x);
  c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

  c2 = p7P_CODON2(x, w);
  c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

  ivx[1] =                     MMX(L-1,1) + p7P_MSC_CODON(gm_fs, 1, c1);
  ivx[1] = p7_FLogsum( ivx[1], MMX(L,1)   + p7P_MSC_CODON(gm_fs, 1, c2));

  XMX(L-2,p7G_B)   =  ivx[1] + TSC(p7P_BM,0);
  for (k = 2; k <= M; k++)
  {
    ivx[k] =                     MMX(L-1,k) + p7P_MSC_CODON(gm_fs, k, c1);
    ivx[k] = p7_FLogsum( ivx[k], MMX(L,k)   + p7P_MSC_CODON(gm_fs, k, c2));
   
    XMX(L-2,p7G_B) = p7_FLogsum( XMX(L-2,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
  }
  
  XMX(L-2,p7G_J) = XMX(L-2,p7G_B) + gm_fs->xsc[p7P_J][p7P_MOVE];
  XMX(L-2,p7G_N) = XMX(L-2,p7G_B) + gm_fs->xsc[p7P_N][p7P_MOVE];
  XMX(L-2,p7G_C) = gm_fs->xsc[p7P_C][p7P_LOOP] + gm_fs->xsc[p7P_C][p7P_MOVE];
  XMX(L-2,p7G_E) = p7_FLogsum( XMX(L-2,p7G_J) + gm_fs->xsc[p7P_E][p7P_LOOP],
                               XMX(L-2,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE]);

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
  
  /* Initialization of rows L-3 and L-4  */
  t = u = v = p7P_MAXCODONS;
  for (i = L-3; i > L-5; i--)
  {
    u = v;
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder value */  
    if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[i+1])) x = dsq[i+1];
    else                                              x = p7P_MAXCODONS; 
  
    /* find correct index for looking up scores of codons and quasicodons */ 
    c1 = p7P_CODON1(x);
    c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

    c2 = p7P_CODON2(x, w);
    c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

    c3 = p7P_CODON3(x, w, v);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    c4 = p7P_CODON4(x, w, v, u);
    c4 = p7P_MINIDX(c4, p7P_DEGEN_QC1);

    ivx[1] =                    MMX(i+1,1) + p7P_MSC_CODON(gm_fs, 1, c1);
    ivx[1] = p7_FLogsum(ivx[1], MMX(i+2,1) + p7P_MSC_CODON(gm_fs, 1, c2));
    ivx[1] = p7_FLogsum(ivx[1], MMX(i+3,1) + p7P_MSC_CODON(gm_fs, 1, c3));
    if( i == L-4 )
      ivx[1] = p7_FLogsum(ivx[1], MMX(i+4,1) + p7P_MSC_CODON(gm_fs, 1, c4));
  
    XMX(i,p7G_B)   =  ivx[1] + TSC(p7P_BM,0);
 
    for (k = 2; k <= M; k++) 
    {
      ivx[k] =                    MMX(i+1,k) + p7P_MSC_CODON(gm_fs, k, c1);
      ivx[k] = p7_FLogsum(ivx[k], MMX(i+2,k) + p7P_MSC_CODON(gm_fs, k, c2));
      ivx[k] = p7_FLogsum(ivx[k], MMX(i+3,k) + p7P_MSC_CODON(gm_fs, k, c3));
      if( i == L-4 )
        ivx[k] = p7_FLogsum(ivx[k], MMX(i+4,k) + p7P_MSC_CODON(gm_fs, k, c4));

      XMX(i,p7G_B) = p7_FLogsum( XMX(i,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
    }

    XMX(i,p7G_J) = p7_FLogsum( XMX(i+3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP],
                               XMX(i,  p7G_B) + gm_fs->xsc[p7P_J][p7P_MOVE]);

    XMX(i,p7G_C) =             XMX(i+3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP];

    XMX(i,p7G_N) = p7_FLogsum( XMX(i+3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP],
                               XMX(i,  p7G_B) + gm_fs->xsc[p7P_N][p7P_MOVE]);

    XMX(i,p7G_E) = p7_FLogsum(XMX(i,p7G_J) + gm_fs->xsc[p7P_E][p7P_LOOP],
                              XMX(i,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE]);

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

  /* Main recursion */
  for (i = L-5; i > 0; i--)
  {
    t = u;
    u = v;
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder value */  
    if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[i+1])) x = dsq[i+1];
    else                                              x = p7P_MAXCODONS; 
  
    /* find correct index for looking up scores of codons and quasicodons */ 
    c1 = p7P_CODON1(x);
    c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

    c2 = p7P_CODON2(x, w);
    c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

    c3 = p7P_CODON3(x, w, v);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    c4 = p7P_CODON4(x, w, v, u);
    c4 = p7P_MINIDX(c4, p7P_DEGEN_QC1);

    c5 = p7P_CODON5(x, w, v, u, t);
    c5 = p7P_MINIDX(c5, p7P_DEGEN_QC2);

    ivx[1] = p7_FLogsum( MMX(i+1,1) + p7P_MSC_CODON(gm_fs, 1, c1), 
            p7_FLogsum( MMX(i+2,1) + p7P_MSC_CODON(gm_fs, 1, c2), 
            p7_FLogsum( MMX(i+3,1) + p7P_MSC_CODON(gm_fs, 1, c3),
            p7_FLogsum( MMX(i+4,1) + p7P_MSC_CODON(gm_fs, 1, c4),
                        MMX(i+5,1) + p7P_MSC_CODON(gm_fs, 1, c5)))));
 
    XMX(i,p7G_B) = ivx[1] + TSC(p7P_BM,0);

    for (k = 2; k <= M; k++) {
      ivx[k] = p7_FLogsum( MMX(i+1,k) + p7P_MSC_CODON(gm_fs, k, c1),
              p7_FLogsum( MMX(i+2,k) + p7P_MSC_CODON(gm_fs, k, c2),
              p7_FLogsum( MMX(i+3,k) + p7P_MSC_CODON(gm_fs, k, c3),
              p7_FLogsum( MMX(i+4,k) + p7P_MSC_CODON(gm_fs, k, c4),
                          MMX(i+5,k) + p7P_MSC_CODON(gm_fs, k, c5)))));

      XMX(i,p7G_B) = p7_FLogsum( XMX(i, p7G_B), ivx[k] + TSC(p7P_BM,k-1));  
    }
  
    XMX(i,p7G_J) = p7_FLogsum( XMX(i+3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP],
                               XMX(i,  p7G_B) + gm_fs->xsc[p7P_J][p7P_MOVE]);
    XMX(i,p7G_C) =             XMX(i+3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP];
    XMX(i,p7G_N) = p7_FLogsum( XMX(i+3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP],
                               XMX(i,  p7G_B) + gm_fs->xsc[p7P_N][p7P_MOVE]);
    XMX(i,p7G_E) = p7_FLogsum(XMX(i,p7G_J) + gm_fs->xsc[p7P_E][p7P_LOOP],
                              XMX(i,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE]);
    
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

  /* At i=0, only N,B states are reachable. */
    t = u;
    u = v;
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder value */  
    if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[1])) x = dsq[1];
    else                                              x = p7P_MAXCODONS; 
  
    /* find correct index for looking up scores of codons and quasicodons */ 
    c1 = p7P_CODON1(x);
    c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

    c2 = p7P_CODON2(x, w);
    c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

    c3 = p7P_CODON3(x, w, v);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    c4 = p7P_CODON4(x, w, v, u);
    c4 = p7P_MINIDX(c4, p7P_DEGEN_QC1);

    c5 = p7P_CODON5(x, w, v, u, t);
    c5 = p7P_MINIDX(c5, p7P_DEGEN_QC2);

  ivx[1] = p7_FLogsum( MMX(1,1) + p7P_MSC_CODON(gm_fs, 1, c1), 
          p7_FLogsum( MMX(2,1) + p7P_MSC_CODON(gm_fs, 1, c2),
          p7_FLogsum( MMX(3,1) + p7P_MSC_CODON(gm_fs, 1, c3),
          p7_FLogsum( MMX(4,1) + p7P_MSC_CODON(gm_fs, 1, c4),
                      MMX(5,1) + p7P_MSC_CODON(gm_fs, 1, c5)))));
                      

  XMX(0,p7G_B) = ivx[1] + TSC(p7P_BM,0); 
    
  for (k = 2; k <= M; k++) 
  {
   ivx[k] =  p7_FLogsum( MMX(1,k) + p7P_MSC_CODON(gm_fs, k, c1),
            p7_FLogsum( MMX(2,k) + p7P_MSC_CODON(gm_fs, k, c2),
            p7_FLogsum( MMX(3,k) + p7P_MSC_CODON(gm_fs, k, c3),
            p7_FLogsum( MMX(4,k) + p7P_MSC_CODON(gm_fs, k, c4),
                        MMX(5,k) + p7P_MSC_CODON(gm_fs, k, c5)))));
 
    XMX(0,p7G_B) = p7_FLogsum(XMX(0, p7G_B), ivx[k] + TSC(p7P_BM,k-1)); 
  }

  XMX(0,p7G_J) = XMX(0,p7G_C) = XMX(0,p7G_E) = -eslINFINITY;
 
  XMX(0,p7G_N) = p7_FLogsum( XMX(3,p7G_N)   + gm_fs->xsc[p7P_N][p7P_LOOP],
                             XMX(0,  p7G_B) + gm_fs->xsc[p7P_N][p7P_MOVE]); 

  for (k = M; k >= 0; k--)
    MMX(0,k) = DMX(0,k) =  IMX(0,k) = -eslINFINITY;           

  if (opt_sc != NULL) *opt_sc =  p7_FLogsum( XMX(0,p7G_N),
                                 p7_FLogsum( XMX(1,p7G_N),
                                             XMX(2,p7G_N)));
  
  gx->M = M;
  gx->L = L;
  
  return eslOK;

}

/* Function:  p7_BackwardPraser_Frameshift_3Codons()
 * Synopsis:  The Backward algorithm.
 *
 * Purpose:   The Backward dynamic programming algorithm - low memeory, 3 codon lengths.
 * 
 *            Given a digital sequence <dsq> of length <L>, a profile <gm_fs>, 
 *            and DP matrix <gx> allocated for at least <L>  cells for the 
 *            special states, and <PARSER_ROWS_BWD> by <gm_fs->M> cells for 
 *            the core model states <M,I,D>; calculate the probability of the 
 *            sequence given the model using the Backward algorithm; use the 
 *            intermediate value matrix to store partial calculations; return 
 *            the Backward matrix in <gx>, and the Backward score in <ret_sc>.
 *           
 *            The Backward score is in lod score form. To convert to a
 *            bitscore, the caller needs to subtract a null model lod
 *            score, then convert to bits.
 *
 * Args:      dsq    - nucleotide sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gcode  - genetic code table
 *            gm_fs  - framshift aware codon profile.
 *            gx     - DP matrix with room for an MxL alignment
 *            iv     - intermediate value martix
 *            opt_sc - optRETURN: Backward lod score in nats
 *           
 * Return:    <eslOK> on success.
 */
int
p7_BackwardParser_Frameshift_3Codons(const ESL_DSQ *dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx, P7_IVX *iv, float *opt_sc)
{

  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  float       *ivx  = iv->ivx;
  int          M    = gm_fs->M;
  int          i, k;
  int          c2, c3, c4;
  float        esc  = p7_fs_profile_IsLocal(gm_fs) ? 0 : -eslINFINITY;
  int          u, v, w, x;
  int          curr, prev2, prev3, prev4;

  for(k = 0; k <= M; k++)
    ivx[k] = -eslINFINITY;

  /* Initialization of row L and L-1  */
  XMX(L,p7G_C)   =                               gm_fs->xsc[p7P_C][p7P_MOVE];
  XMX(L-1,p7G_C) = gm_fs->xsc[p7P_C][p7P_LOOP] + gm_fs->xsc[p7P_C][p7P_MOVE];
  for(i = L; i >= L-1; i--) {
    curr = i % PARSER_ROWS_BWD;
    XMX(i,p7G_J) = XMX(i,p7G_B) = XMX(i,p7G_N) = -eslINFINITY; /* need to enter and exit model */
    XMX(i,p7G_E) = XMX(i,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE];    /* E<-C, no tail */
    MMX(curr,M)  = DMX(curr,M) = XMX(i,p7G_E);                    /* {MD}_M <- E (prob 1.0) */
    IMX(curr,M)  = -eslINFINITY;                               /* no I_M state        */

    for (k = M-1; k >= 1; k--) 
    {
   
      MMX(curr,k) = p7_FLogsum( XMX(i,p7G_E)    + esc,
                                DMX(curr, k+1)  + TSC(p7P_MD,k));

      DMX(curr,k) = p7_FLogsum( XMX(i,p7G_E)    + esc,
                                DMX(curr, k+1)  + TSC(p7P_DD,k));

      IMX(curr,k) = -eslINFINITY;
    }
    MMX(curr,0) = IMX(curr,0) = DMX(curr,0)  = -eslINFINITY;
  }
  
  for (k = M-1; k >= 1; k--) 
  {
    MMX(curr,k) = p7_FLogsum( XMX(L-1,p7G_E) + esc,
                              DMX(curr, k+1) + TSC(p7P_MD,k));

    DMX(curr,k) = p7_FLogsum( XMX(L-1,p7G_E)    + esc,
                              DMX(curr, k+1)  + TSC(p7P_DD,k));

    IMX(curr,k) = -eslINFINITY;
  }
  /* Check for degenerate nucleotides */

  /* Initialization of row L-2  */
  if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[L])) w = dsq[L];
  else                                              w = p7P_MAXCODONS;

  if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[L-1])) x = dsq[L-1];
  else                                              x = p7P_MAXCODONS;
  
  
  c2 = p7P_CODON2(x, w);
  c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

  curr  = (L-2) % PARSER_ROWS_BWD;
  prev2 = L     % PARSER_ROWS_BWD;

  ivx[1]       = MMX(prev2,1) + p7P_MSC_CODON(gm_fs, 1, c2); 
  XMX(L-2,p7G_B) = ivx[1]      + TSC(p7P_BM,0);
  for (k = 2; k <= M; k++)
  {
    ivx[k] =  MMX(prev2,k) + p7P_MSC_CODON(gm_fs, k, c2);

    XMX(L-2,p7G_B) = p7_FLogsum( XMX(L-2,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
  }

  XMX(L-2,p7G_J) = XMX(L-2,p7G_B) + gm_fs->xsc[p7P_J][p7P_MOVE];
  XMX(L-2,p7G_N) = XMX(L-2,p7G_B) + gm_fs->xsc[p7P_N][p7P_MOVE];
  XMX(L-2,p7G_C) = gm_fs->xsc[p7P_C][p7P_LOOP] + gm_fs->xsc[p7P_C][p7P_MOVE];
  
  XMX(L-2,p7G_E) = p7_FLogsum(XMX(L-2,p7G_J) + gm_fs->xsc[p7P_E][p7P_LOOP],
                              XMX(L-2,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE]);

  MMX(curr,M)  = DMX(curr,M) = XMX(L-2,p7G_E); 
  IMX(curr,M)  = -eslINFINITY;       

  for (k = M-1; k >= 1; k--)
  {
    MMX(curr,k) = p7_FLogsum( DMX(curr,k+1)   + TSC(p7P_MD,k),
                  p7_FLogsum( ivx[k+1]        + TSC(p7P_MM,k),
                              XMX(L-2,p7G_E)    + esc));

    DMX(curr,k) = p7_FLogsum( p7_FLogsum( XMX(L-2,p7G_E) + esc,
                                          DMX(curr, k+1)  + TSC(p7P_DD,k)),
                                          ivx[k+1]        + TSC(p7P_DM,k));

    IMX(curr,k) = ivx[k+1]            + TSC(p7P_IM,k);
  }

  MMX(curr,0) = IMX(curr,0) = DMX(curr,0)  = -eslINFINITY;

  /* Initialization of rows L-3 and L-4  */
  u = v = p7P_MAXCODONS; 
  for (i = L-3; i > L-5; i--)
  {
    u = v;
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder value */  
    if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[i+1])) x = dsq[i+1];
    else                                              x = p7P_MAXCODONS; 
  
    /* find correct index for looking up scores of codons and quasicodons */ 
    c2 = p7P_CODON2(x, w);
    c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

    c3 = p7P_CODON3(x, w, v);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    c4 = p7P_CODON4(x, w, v, u);
    c4 = p7P_MINIDX(c4, p7P_DEGEN_QC1);

    curr  =  i    % PARSER_ROWS_BWD;
    prev2 = (i+2) % PARSER_ROWS_BWD;
    prev3 = (i+3) % PARSER_ROWS_BWD;
    prev4 = (i+4) % PARSER_ROWS_BWD;

    ivx[1] =                     MMX(prev2,1) + p7P_MSC_CODON(gm_fs, 1, c2);
 
    ivx[1] = p7_FLogsum( ivx[1], MMX(prev3,1) + p7P_MSC_CODON(gm_fs, 1, c3));

    if( i == L-4 )
      ivx[1] = p7_FLogsum( ivx[1], MMX(prev4,1) + p7P_MSC_CODON(gm_fs, 1, c4));
  
    XMX(i,p7G_B)   =  ivx[1] + TSC(p7P_BM,0);

     for (k = 2; k <= M; k++) 
    {
      ivx[k] =                     MMX(prev2,k) + p7P_MSC_CODON(gm_fs, k, c2);
      
      ivx[k] = p7_FLogsum( ivx[k], MMX(prev3,k) + p7P_MSC_CODON(gm_fs, k, c3));

      if( i == L-4 )
        ivx[k]  = p7_FLogsum( ivx[k], MMX(prev4,k) + p7P_MSC_CODON(gm_fs, k, c4));

      XMX(i,p7G_B) = p7_FLogsum( XMX(i,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
    }

    XMX(i,p7G_J) = p7_FLogsum( XMX(i+3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP],
                               XMX(i,  p7G_B) + gm_fs->xsc[p7P_J][p7P_MOVE]);

    XMX(i,p7G_C) =             XMX(i+3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP];

    XMX(i,p7G_N) = p7_FLogsum( XMX(i+3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP],
                               XMX(i,  p7G_B) + gm_fs->xsc[p7P_N][p7P_MOVE]);

    XMX(i,p7G_E) = p7_FLogsum(XMX(i,p7G_J) + gm_fs->xsc[p7P_E][p7P_LOOP],
                              XMX(i,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE]);

    MMX(curr,M)     = DMX(curr,M) = XMX(i,p7G_E); /* {MD}_M <- E (prob 1.0) */
    IMX(curr,M)     = -eslINFINITY;            /* no I_M state        */

    for (k = M-1; k >= 1; k--)
    {
      MMX(curr,k) = p7_FLogsum( DMX(curr,k+1) + TSC(p7P_MD,k),
                    p7_FLogsum( IMX(prev3,k)  + TSC(p7P_MI,k),
                    p7_FLogsum( ivx[k+1]      + TSC(p7P_MM,k),
                                XMX(i,p7G_E)  + esc)));
        
      DMX(curr,k) = p7_FLogsum( DMX(curr,k+1) + TSC(p7P_DD,k),
                    p7_FLogsum( XMX(i,p7G_E)  + esc,
                                ivx[k+1]      + TSC(p7P_DM,k)));

      IMX(curr,k) = p7_FLogsum( IMX(prev3,k)  + TSC(p7P_II,k),
                                ivx[k+1]       + TSC(p7P_IM,k));
    }

    MMX(curr,0) = IMX(curr,0) = DMX(curr,0)  = -eslINFINITY;

  }

  /* Main recursion */
  for (i = L-5; i > 0; i--)
  {
    u = v;
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder value */  
    if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[i+1])) x = dsq[i+1];
    else                                              x = p7P_MAXCODONS; 
  
    /* find correct index for looking up scores of codons and quasicodons */ 
    c2 = p7P_CODON2(x, w);
    c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

    c3 = p7P_CODON3(x, w, v);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    c4 = p7P_CODON4(x, w, v, u);
    c4 = p7P_MINIDX(c4, p7P_DEGEN_QC1);

    curr  =  i    % PARSER_ROWS_BWD;
    prev2 = (i+2) % PARSER_ROWS_BWD;
    prev3 = (i+3) % PARSER_ROWS_BWD;
    prev4 = (i+4) % PARSER_ROWS_BWD;

    ivx[1] = p7_FLogsum( MMX(prev2,1) + p7P_MSC_CODON(gm_fs, 1, c2), 
             p7_FLogsum( MMX(prev3,1) + p7P_MSC_CODON(gm_fs, 1, c3),
                         MMX(prev4,1) + p7P_MSC_CODON(gm_fs, 1, c4)));
            
    XMX(i,p7G_B) = ivx[1] + TSC(p7P_BM,0);

    for (k = 2; k <= M; k++) {
      ivx[k] = p7_FLogsum( MMX(prev2,k) + p7P_MSC_CODON(gm_fs, k, c2),
               p7_FLogsum( MMX(prev3,k) + p7P_MSC_CODON(gm_fs, k, c3),
                           MMX(prev4,k) + p7P_MSC_CODON(gm_fs, k, c4)));

      XMX(i,p7G_B) = p7_FLogsum( XMX(i, p7G_B), ivx[k] + TSC(p7P_BM,k-1));  
    }
  
    XMX(i,p7G_J) = p7_FLogsum( XMX(i+3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP],
                               XMX(i,  p7G_B) + gm_fs->xsc[p7P_J][p7P_MOVE]);

    XMX(i,p7G_C) =             XMX(i+3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP];

    XMX(i,p7G_N) = p7_FLogsum( XMX(i+3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP],
                               XMX(i,  p7G_B) + gm_fs->xsc[p7P_N][p7P_MOVE]);

    XMX(i,p7G_E) = p7_FLogsum( XMX(i,p7G_J) + gm_fs->xsc[p7P_E][p7P_LOOP],
                               XMX(i,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE]);
    
    MMX(curr,M)     = DMX(curr,M) = XMX(i,p7G_E); /* {MD}_M <- E (prob 1.0) */
    IMX(curr,M)     = -eslINFINITY;            /* no I_M state        */

    for (k = M-1; k >= 1; k--)
    {
      MMX(curr,k) = p7_FLogsum( p7_FLogsum( DMX(curr,k+1)   + TSC(p7P_MD,k),
                                p7_FLogsum( IMX(prev3,k)   + TSC(p7P_MI,k),
                                            ivx[k+1]      + TSC(p7P_MM,k))),
                                            XMX(i,p7G_E) + esc);

      DMX(curr,k) = p7_FLogsum( p7_FLogsum( XMX(i,p7G_E) + esc,
                                            DMX(curr, k+1)  + TSC(p7P_DD,k)),
                                            ivx[k+1]      + TSC(p7P_DM,k));

      IMX(curr,k) = p7_FLogsum(             IMX(prev3,k  ) + TSC(p7P_II,k),
                                            ivx[k+1]      + TSC(p7P_IM,k));
    }

    MMX(curr,0) = IMX(curr,0) = DMX(curr,0)  = -eslINFINITY;

  }

  /* At i=0, only N,B states are reachable. */
  u = v;
  v = w;
  w = x;

  /* if new nucleotide is not A,C,G, or T set it to placeholder value */  
  if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[1])) x = dsq[1];
  else                                              x = p7P_MAXCODONS; 
  
  /* find correct index for looking up scores of codons and quasicodons */ 
  c2 = p7P_CODON2(x, w);
  c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

  c3 = p7P_CODON3(x, w, v);
  c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

  c4 = p7P_CODON4(x, w, v, u);
  c4 = p7P_MINIDX(c4, p7P_DEGEN_QC1);

  prev2 = 2 % PARSER_ROWS_BWD;
  prev3 = 3 % PARSER_ROWS_BWD;
  prev4 = 4 % PARSER_ROWS_BWD;

  ivx[1] = p7_FLogsum( MMX(prev2,1) + p7P_MSC_CODON(gm_fs, 1, c2), 
           p7_FLogsum( MMX(prev3,1) + p7P_MSC_CODON(gm_fs, 1, c3),
                       MMX(prev4,1) + p7P_MSC_CODON(gm_fs, 1, c4)));
                      
  XMX(0,p7G_B) = ivx[1] + TSC(p7P_BM,0); 

  for (k = 2; k <= M; k++) 
  {
   ivx[k] =  p7_FLogsum( MMX(prev2,k) + p7P_MSC_CODON(gm_fs, k, c2),
             p7_FLogsum( MMX(prev3,k) + p7P_MSC_CODON(gm_fs, k, c3),
                         MMX(prev4,k) + p7P_MSC_CODON(gm_fs, k, c4)));
 
    XMX(0,p7G_B) = p7_FLogsum(XMX(0, p7G_B), ivx[k] + TSC(p7P_BM,k-1)); 
  }

  XMX(0,p7G_J) = XMX(0,p7G_C) = XMX(0,p7G_E) = -eslINFINITY;
 
  XMX(0,p7G_N) = p7_FLogsum( XMX(3,p7G_N)   + gm_fs->xsc[p7P_N][p7P_LOOP],
                             XMX(0,  p7G_B) + gm_fs->xsc[p7P_N][p7P_MOVE]); 

  for (k = M; k >= 0; k--)
    MMX(0,k) = DMX(0,k) =  IMX(0,k) = -eslINFINITY;           

  if (opt_sc != NULL) *opt_sc =  p7_FLogsum( XMX(0,p7G_N),
                                 p7_FLogsum( XMX(1,p7G_N),
                                             XMX(2,p7G_N)));
  gx->M = M;
  gx->L = L;
  
  return eslOK;

}

/* Function:  p7_BackwardPraser_Frameshift_5Codons()
 * Synopsis:  The Backward algorithm.
 *
 * Purpose:   The Backward dynamic programming algorithm - low memeory, 5 codon lengths.
 *
 *            Given a digital sequence <dsq> of length <L>, a profile <gm_fs>,
 *            and DP matrix <gx> allocated for at least <L>  cells for the
 *            special states, and <PARSER_ROWS_BWD> by <gm_fs->M> cells for
 *            the core model states <M,I,D>; calculate the probability of the
 *            sequence given the model using the Backward algorithm; use the
 *            intermediate value matrix to store partial calculations; return
 *            the Backward matrix in <gx>, and the Backward score in <ret_sc>.
 *
 *            The Backward score is in lod score form. To convert to a
 *            bitscore, the caller needs to subtract a null model lod
 *            score, then convert to bits.
 *
 * Args:      dsq    - nucleotide sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gcode  - genetic code table
 *            gm_fs  - framshift aware codon profile.
 *            gx     - DP matrix with room for an MxL alignment
 *            iv     - intermediate value martix
 *            opt_sc - optRETURN: Backward lod score in nats
 *
 * Return:    <eslOK> on success.
 */
int
p7_BackwardParser_Frameshift_5Codons(const ESL_DSQ *dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx, P7_IVX *iv, float *opt_sc)
{

  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;           
  float       *ivx   = iv->ivx;
  int          M    = gm_fs->M;
  int          i, k;;  
  int          c1, c2, c3, c4, c5;
  float        esc  = p7_fs_profile_IsLocal(gm_fs) ? 0 : -eslINFINITY;
  int          t, u, v, w, x;
  int          curr, prev1, prev2, prev3, prev4, prev5;

  for(k = 0; k <= M; k++)
    ivx[k] = -eslINFINITY;

  /* Initialization of row L  */
  curr  =  L    % PARSER_ROWS_BWD;
  XMX(L,p7G_J) = XMX(L,p7G_B) = XMX(L,p7G_N) = -eslINFINITY; /* need to enter and exit model */
  XMX(L,p7G_C) = gm_fs->xsc[p7P_C][p7P_MOVE];                   /* C<-T         */
  XMX(L,p7G_E) = XMX(L,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE];    /* E<-C, no tail */
  MMX(curr,M)     = DMX(curr,M) = XMX(L,p7G_E);                    /* {MD}_M <- E (prob 1.0) */
  IMX(curr,M)     = -eslINFINITY;                               /* no I_M state        */

  /* Initialization of core model for row L */
  for (k = M-1; k >= 1; k--) 
  {
    /*L comes form E & D state only */
    MMX(curr,k) = p7_FLogsum( XMX(L,p7G_E) + esc,
                              DMX(curr, k+1)  + TSC(p7P_MD,k));
 
    DMX(curr,k) = p7_FLogsum( XMX(L,p7G_E) + esc,
                              DMX(curr, k+1)  + TSC(p7P_DD,k));

    /* no I state at L to L-2 */
    IMX(curr,k) = -eslINFINITY;
  }
  MMX(curr,0) = IMX(curr,0) = DMX(curr,0)  = -eslINFINITY;

  /* Initialization of row L - 1  */
  if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[L])) x = dsq[L];
  else                                            x = p7P_MAXCODONS;

  c1 = p7P_CODON1(x);
  c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);
  
  curr  = (L-1)    % PARSER_ROWS_BWD;
  prev1 =  L       % PARSER_ROWS_BWD;
  
  ivx[1]       = MMX(prev1,1) + p7P_MSC_CODON(gm_fs, 1, c1);
  XMX(L-1,p7G_B) = ivx[1]   + TSC(p7P_BM,0);
  for (k = 2; k <= M; k++)
  {
    ivx[k] =  MMX(prev1,k) + p7P_MSC_CODON(gm_fs, k, c1);

    XMX(L-1,p7G_B) = p7_FLogsum( XMX(L-1,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
  }

  XMX(L-1,p7G_J) = XMX(L-1,p7G_B) + gm_fs->xsc[p7P_J][p7P_MOVE];
  XMX(L-1,p7G_N) = XMX(L-1,p7G_B) + gm_fs->xsc[p7P_N][p7P_MOVE];
  XMX(L-1,p7G_C) = gm_fs->xsc[p7P_C][p7P_LOOP] + gm_fs->xsc[p7P_C][p7P_MOVE];
  XMX(L-1,p7G_E) = p7_FLogsum( XMX(L-1,p7G_J) + gm_fs->xsc[p7P_E][p7P_LOOP],
				               XMX(L-1,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE]);

  MMX(curr,M)  = DMX(curr,M) = XMX(L-1,p7G_E); 
  for (k = M-1; k >= 1; k--)
  {
    MMX(curr,k) = p7_FLogsum( DMX(curr,k+1)   + TSC(p7P_MD,k),
                 p7_FLogsum( ivx[k+1]       + TSC(p7P_MM,k),
                             XMX(L-1,p7G_E)  + esc));

    DMX(curr,k) = p7_FLogsum( p7_FLogsum( XMX(L-1,p7G_E) + esc,
                                          DMX(curr, k+1)  + TSC(p7P_DD,k)),
                                          ivx[k+1]       + TSC(p7P_DM,k));

    IMX(curr,k) = ivx[k+1]            + TSC(p7P_IM,k);
  }

  MMX(curr,0) = IMX(curr,0) = DMX(curr,0)  = -eslINFINITY;

  /* Initialization of row L - 2  */
  w = x;

  if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[L-1])) x = dsq[L-1];
  else                                            x = p7P_MAXCODONS;

  c1 = p7P_CODON1(x);
  c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

  c2 = p7P_CODON2(x, w);
  c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

  curr  = (L-2)    % PARSER_ROWS_BWD;
  prev1 = (L-1)      % PARSER_ROWS_BWD;
  prev2 =  L       % PARSER_ROWS_BWD;
  ivx[1] =                     MMX(prev1,1) + p7P_MSC_CODON(gm_fs, 1, c1);
  ivx[1] = p7_FLogsum( ivx[1], MMX(prev2,1) + p7P_MSC_CODON(gm_fs, 1, c2));

  XMX(L-2,p7G_B)   =  ivx[1] + TSC(p7P_BM,0);
  for (k = 2; k <= M; k++)
  {
    ivx[k] =                     MMX(prev1,k) + p7P_MSC_CODON(gm_fs, k, c1);
    ivx[k] = p7_FLogsum( ivx[k], MMX(prev2,k) + p7P_MSC_CODON(gm_fs, k, c2));
   
    XMX(L-2,p7G_B) = p7_FLogsum( XMX(L-2,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
  }
  
  XMX(L-2,p7G_J) = XMX(L-2,p7G_B) + gm_fs->xsc[p7P_J][p7P_MOVE];
  XMX(L-2,p7G_N) = XMX(L-2,p7G_B) + gm_fs->xsc[p7P_N][p7P_MOVE];
  XMX(L-2,p7G_C) = gm_fs->xsc[p7P_C][p7P_LOOP] + gm_fs->xsc[p7P_C][p7P_MOVE];
  XMX(L-2,p7G_E) = p7_FLogsum( XMX(L-2,p7G_J) + gm_fs->xsc[p7P_E][p7P_LOOP],
                               XMX(L-2,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE]);

  MMX(curr,M)  = DMX(curr,M) = XMX(L-2,p7G_E);
  for (k = M-1; k >= 1; k--)
  {
    MMX(curr,k) = p7_FLogsum( DMX(curr,k+1)   + TSC(p7P_MD,k),
                 p7_FLogsum( ivx[k+1]       + TSC(p7P_MM,k),
                             XMX(L-2,p7G_E)  + esc));

    DMX(curr,k) = p7_FLogsum( p7_FLogsum( XMX(L-2,p7G_E) + esc,
                                          DMX(curr, k+1)  + TSC(p7P_DD,k)),
                                          ivx[k+1]       + TSC(p7P_DM,k));

    IMX(curr,k) = ivx[k+1]            + TSC(p7P_IM,k);
  }
  MMX(curr,0) = IMX(curr,0) = DMX(curr,0)  = -eslINFINITY;
  
  /* Initialization of rows L-3 and L-4  */
  t = u = v = p7P_MAXCODONS;
  for (i = L-3; i > L-5; i--)
  {
    u = v;
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder value */  
    if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[i+1])) x = dsq[i+1];
    else                                              x = p7P_MAXCODONS; 
  
    /* find correct index for looking up scores of codons and quasicodons */ 
    c1 = p7P_CODON1(x);
    c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

    c2 = p7P_CODON2(x, w);
    c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

    c3 = p7P_CODON3(x, w, v);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    c4 = p7P_CODON4(x, w, v, u);
    c4 = p7P_MINIDX(c4, p7P_DEGEN_QC1);

    curr  =  i    % PARSER_ROWS_BWD;
    prev1 = (i+1) % PARSER_ROWS_BWD;
    prev2 = (i+2) % PARSER_ROWS_BWD;
    prev3 = (i+3) % PARSER_ROWS_BWD;
    prev4 = (i+4) % PARSER_ROWS_BWD; 

    ivx[1] =                    MMX(prev1,1) + p7P_MSC_CODON(gm_fs, 1, c1);
    ivx[1] = p7_FLogsum(ivx[1], MMX(prev2,1) + p7P_MSC_CODON(gm_fs, 1, c2));
    ivx[1] = p7_FLogsum(ivx[1], MMX(prev3,1) + p7P_MSC_CODON(gm_fs, 1, c3));
    if( i == L-4 )
      ivx[1] = p7_FLogsum(ivx[1], MMX(prev4,1) + p7P_MSC_CODON(gm_fs, 1, c4));
  
    XMX(i,p7G_B)   =  ivx[1] + TSC(p7P_BM,0);
 
    for (k = 2; k <= M; k++) 
    {
      ivx[k] =                    MMX(prev1,k) + p7P_MSC_CODON(gm_fs, k, c1);
      ivx[k] = p7_FLogsum(ivx[k], MMX(prev2,k) + p7P_MSC_CODON(gm_fs, k, c2));
      ivx[k] = p7_FLogsum(ivx[k], MMX(prev3,k) + p7P_MSC_CODON(gm_fs, k, c3));
      if( i == L-4 )
        ivx[k] = p7_FLogsum(ivx[k], MMX(prev4,k) + p7P_MSC_CODON(gm_fs, k, c4));

      XMX(i,p7G_B) = p7_FLogsum( XMX(i,p7G_B), ivx[k] + TSC(p7P_BM,k-1));
    }

    XMX(i,p7G_J) = p7_FLogsum( XMX(i+3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP],
                               XMX(i,  p7G_B) + gm_fs->xsc[p7P_J][p7P_MOVE]);

    XMX(i,p7G_C) =             XMX(i+3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP];

    XMX(i,p7G_N) = p7_FLogsum( XMX(i+3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP],
                               XMX(i,  p7G_B) + gm_fs->xsc[p7P_N][p7P_MOVE]);

    XMX(i,p7G_E) = p7_FLogsum(XMX(i,p7G_J) + gm_fs->xsc[p7P_E][p7P_LOOP],
                              XMX(i,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE]);

    MMX(curr,M)     =            DMX(curr,M) = XMX(i,p7G_E); /* {MD}_M <- E (prob 1.0) */
 
    IMX(curr,M)     =            -eslINFINITY;            /* no I_M state        */

    for (k = M-1; k >= 1; k--)
    {
      MMX(curr,k) = p7_FLogsum( p7_FLogsum( DMX(curr,k+1)   + TSC(p7P_MD,k),
                                p7_FLogsum( IMX(prev3,k)   + TSC(p7P_MI,k),
                                            ivx[k+1]      + TSC(p7P_MM,k))),
                                            XMX(i,p7G_E) + esc);

      DMX(curr,k) = p7_FLogsum( p7_FLogsum( XMX(i,p7G_E) + esc,
                                            DMX(curr, k+1)  + TSC(p7P_DD,k)),
                                            ivx[k+1]     + TSC(p7P_DM,k));

      IMX(curr,k) = p7_FLogsum(           IMX(prev3,k  )   + TSC(p7P_II,k),
                                          ivx[k+1]         + TSC(p7P_IM,k));
    }

    MMX(curr,0) = IMX(curr,0) = DMX(curr,0)  = -eslINFINITY;

  }

  /* Main recursion */
  for (i = L-5; i > 0; i--)
  {
    t = u;
    u = v;
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder value */  
    if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[i+1])) x = dsq[i+1];
    else                                              x = p7P_MAXCODONS; 
  
    /* find correct index for looking up scores of codons and quasicodons */ 
    c1 = p7P_CODON1(x);
    c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

    c2 = p7P_CODON2(x, w);
    c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

    c3 = p7P_CODON3(x, w, v);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    c4 = p7P_CODON4(x, w, v, u);
    c4 = p7P_MINIDX(c4, p7P_DEGEN_QC1);

    c5 = p7P_CODON5(x, w, v, u, t);
    c5 = p7P_MINIDX(c5, p7P_DEGEN_QC2);

    curr  =  i    % PARSER_ROWS_BWD;
    prev1 = (i+1) % PARSER_ROWS_BWD;
    prev2 = (i+2) % PARSER_ROWS_BWD;
    prev3 = (i+3) % PARSER_ROWS_BWD;
    prev4 = (i+4) % PARSER_ROWS_BWD;
    prev5 = (i+5) % PARSER_ROWS_BWD;

    ivx[1] = p7_FLogsum( MMX(prev1,1) + p7P_MSC_CODON(gm_fs, 1, c1), 
            p7_FLogsum( MMX(prev2,1) + p7P_MSC_CODON(gm_fs, 1, c2), 
            p7_FLogsum( MMX(prev3,1) + p7P_MSC_CODON(gm_fs, 1, c3),
            p7_FLogsum( MMX(prev4,1) + p7P_MSC_CODON(gm_fs, 1, c4),
                        MMX(prev5,1) + p7P_MSC_CODON(gm_fs, 1, c5)))));
 
    XMX(i,p7G_B) = ivx[1] + TSC(p7P_BM,0);

    for (k = 2; k <= M; k++) {
      ivx[k] = p7_FLogsum( MMX(prev1,k) + p7P_MSC_CODON(gm_fs, k, c1),
              p7_FLogsum( MMX(prev2,k) + p7P_MSC_CODON(gm_fs, k, c2),
              p7_FLogsum( MMX(prev3,k) + p7P_MSC_CODON(gm_fs, k, c3),
              p7_FLogsum( MMX(prev4,k) + p7P_MSC_CODON(gm_fs, k, c4),
                          MMX(prev5,k) + p7P_MSC_CODON(gm_fs, k, c5)))));

      XMX(i,p7G_B) = p7_FLogsum( XMX(i, p7G_B), ivx[k] + TSC(p7P_BM,k-1));  
    }
  
    XMX(i,p7G_J) = p7_FLogsum( XMX(i+3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP],
                               XMX(i,  p7G_B) + gm_fs->xsc[p7P_J][p7P_MOVE]);
    XMX(i,p7G_C) =             XMX(i+3,p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP];
    XMX(i,p7G_N) = p7_FLogsum( XMX(i+3,p7G_N) + gm_fs->xsc[p7P_N][p7P_LOOP],
                               XMX(i,  p7G_B) + gm_fs->xsc[p7P_N][p7P_MOVE]);
    XMX(i,p7G_E) = p7_FLogsum(XMX(i,p7G_J) + gm_fs->xsc[p7P_E][p7P_LOOP],
                              XMX(i,p7G_C) + gm_fs->xsc[p7P_E][p7P_MOVE]);
    
    MMX(curr,M)     = DMX(curr,M) = XMX(i,p7G_E); /* {MD}_M <- E (prob 1.0) */
    IMX(curr,M)     = -eslINFINITY;            /* no I_M state        */

    for (k = M-1; k >= 1; k--)
    {
      MMX(curr,k) = p7_FLogsum( p7_FLogsum( DMX(curr,k+1)   + TSC(p7P_MD,k),
                             p7_FLogsum( IMX(prev3,k)   + TSC(p7P_MI,k),
                                         ivx[k+1]      + TSC(p7P_MM,k))),
                                         XMX(i,p7G_E) + esc);

      DMX(curr,k) = p7_FLogsum( p7_FLogsum( XMX(i,p7G_E) + esc,
                                         DMX(curr, k+1)  + TSC(p7P_DD,k)),
                                         ivx[k+1]      + TSC(p7P_DM,k));

      IMX(curr,k) = p7_FLogsum(             IMX(prev3,k  ) + TSC(p7P_II,k),
                                         ivx[k+1]      + TSC(p7P_IM,k));
    }

    MMX(curr,0) = IMX(curr,0) = DMX(curr,0)  = -eslINFINITY;

  }

  /* At i=0, only N,B states are reachable. */
    t = u;
    u = v;
    v = w;
    w = x;

    /* if new nucleotide is not A,C,G, or T set it to placeholder value */  
    if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[1])) x = dsq[1];
    else                                              x = p7P_MAXCODONS; 
  
    /* find correct index for looking up scores of codons and quasicodons */ 
    c1 = p7P_CODON1(x);
    c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

    c2 = p7P_CODON2(x, w);
    c2 = p7P_MINIDX(c2, p7P_DEGEN_QC1);

    c3 = p7P_CODON3(x, w, v);
    c3 = p7P_MINIDX(c3, p7P_DEGEN_C);

    c4 = p7P_CODON4(x, w, v, u);
    c4 = p7P_MINIDX(c4, p7P_DEGEN_QC1);

    c5 = p7P_CODON5(x, w, v, u, t);
    c5 = p7P_MINIDX(c5, p7P_DEGEN_QC2);

    prev1 = 1 % PARSER_ROWS_BWD;
    prev2 = 2 % PARSER_ROWS_BWD;
    prev3 = 3 % PARSER_ROWS_BWD;
    prev4 = 4 % PARSER_ROWS_BWD;
    prev5 = 5 % PARSER_ROWS_BWD;

  ivx[1] = p7_FLogsum( MMX(prev1,1) + p7P_MSC_CODON(gm_fs, 1, c1), 
          p7_FLogsum( MMX(prev2,1) + p7P_MSC_CODON(gm_fs, 1, c2),
          p7_FLogsum( MMX(prev3,1) + p7P_MSC_CODON(gm_fs, 1, c3),
          p7_FLogsum( MMX(prev4,1) + p7P_MSC_CODON(gm_fs, 1, c4),
                      MMX(prev5,1) + p7P_MSC_CODON(gm_fs, 1, c5)))));
                      

  XMX(0,p7G_B) = ivx[1] + TSC(p7P_BM,0); 
    
  for (k = 2; k <= M; k++) 
  {
   ivx[k] =  p7_FLogsum( MMX(prev1,k) + p7P_MSC_CODON(gm_fs, k, c1),
            p7_FLogsum( MMX(prev2,k) + p7P_MSC_CODON(gm_fs, k, c2),
            p7_FLogsum( MMX(prev3,k) + p7P_MSC_CODON(gm_fs, k, c3),
            p7_FLogsum( MMX(prev4,k) + p7P_MSC_CODON(gm_fs, k, c4),
                        MMX(prev5,k) + p7P_MSC_CODON(gm_fs, k, c5)))));
 
    XMX(0,p7G_B) = p7_FLogsum(XMX(0, p7G_B), ivx[k] + TSC(p7P_BM,k-1)); 
  }

  XMX(0,p7G_J) = XMX(0,p7G_C) = XMX(0,p7G_E) = -eslINFINITY;
 
  XMX(0,p7G_N) = p7_FLogsum( XMX(3,p7G_N)   + gm_fs->xsc[p7P_N][p7P_LOOP],
                             XMX(0,  p7G_B) + gm_fs->xsc[p7P_N][p7P_MOVE]); 

  if (opt_sc != NULL) *opt_sc =  p7_FLogsum( XMX(0,p7G_N),
                                 p7_FLogsum( XMX(1,p7G_N),
                                             XMX(2,p7G_N)));
  
  gx->M = M;
  gx->L = L;
  
  return eslOK;

}



/*------------- end: forward, backward, hybrid ------------------*/

/*****************************************************************
 * 2. Benchmark driver.
 *****************************************************************/
#ifdef p7FWDBACK_FRAMESHIFT_BENCHMARK
/*
   gcc -g -O2      -o fwdback_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7FWDBACK_FRAMESHIFT_BENCHMARK fwdback_frameshift.c -lhmmer -leasel -lm
   icc -O3 -static -o fwdback_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7FWDBACK_FRAMESHIFT_BENCHMARK fwdback_frameshift.c -lhmmer -leasel -lm
   ./fwdback_frameshift_benchmark <hmmfile>
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
  { "-P",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark Parsers",                         0 },
  { "-U",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "only benchmark full matrix version",             0 },
  { "-T",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  "-P", NULL, "only benchmark 3 codon length Forward Parser",   0 },
  { "-V",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  "-P", NULL, "only benchmark 5 codon length Forward Parser",   0 },
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
  P7_FS_PROFILE  *gm_fs   = NULL;
  P7_GMX         *fwd_p   = NULL;
  P7_GMX         *bck_p   = NULL;
  P7_GMX         *fwd     = NULL;
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
  bgDNA = p7_bg_fs_Create(abcDNA);
  bgAA  = p7_bg_fs_Create(abcAA);
  gm_fs = p7_profile_fs_Create(hmm->M, abcAA);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs, L/3, p7_UNILOCAL);
  fwd_p = p7_gmx_fs_Create(gm_fs->M, PARSER_ROWS_FWD, L, 0);
  bck_p = p7_gmx_fs_Create(gm_fs->M, PARSER_ROWS_BWD, L, 0);
  fwd   = p7_gmx_fs_Create(gm_fs->M, L, L, p7P_5CODONS);
  bck   = p7_gmx_fs_Create(gm_fs->M, L, L, 0);
  iv    = p7_ivx_Create(gm_fs->M, p7P_5CODONS);

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
        if (! esl_opt_GetBoolean(go, "-P"))                                   p7_Forward_Frameshift (dsq, gcode, L, gm_fs, fwd, iv, &sc);
        if (! esl_opt_GetBoolean(go, "-U") && ! esl_opt_GetBoolean(go, "-T")) p7_ForwardParser_Frameshift_5Codons(dsq, gcode, L, gm_fs, fwd_p, iv, &sc);
        p7_gmx_Reuse(fwd_p);
        if (! esl_opt_GetBoolean(go, "-U") && ! esl_opt_GetBoolean(go, "-V")) p7_ForwardParser_Frameshift_3Codons(dsq, gcode, L, gm_fs, fwd_p, iv, &sc);
      } 
      if (! esl_opt_GetBoolean(go, "-F"))  {
        if (! esl_opt_GetBoolean(go, "-P"))                                   p7_Backward_Frameshift(dsq, gcode, L, gm_fs, bck, iv, NULL); 
        if (! esl_opt_GetBoolean(go, "-U"))                                   p7_BackwardParser_Frameshift_3Codons(dsq, gcode, L, gm_fs, bck_p, iv, NULL);
     }

     p7_gmx_Reuse(fwd_p);
     p7_gmx_Reuse(bck_p);
     p7_gmx_Reuse(fwd);
     p7_gmx_Reuse(bck);
    
  }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;
  Mcs        = (double) N * (double) L * (double) gm_fs->M * 1e-6 / (double) bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   gm_fs->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_gmx_Destroy(bck_p);
  p7_gmx_Destroy(fwd_p);
  p7_gmx_Destroy(bck);
  p7_gmx_Destroy(fwd);
  p7_ivx_Destroy(iv);
  p7_profile_fs_Destroy(gm_fs);
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
#endif /*p7FWDBACK_FRAMESHIFT_BENCHMARK*/
/*----------------- end, benchmark ------------------------------*/




/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7FWDBACK_FRAMESHIFT_TESTDRIVE
#include <string.h>
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"

static void
utest_forward_fs(ESL_GETOPTS *go, ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_GENCODE *gcode, P7_CODONTABLE *codon_table, P7_BG *bgAA, P7_HMM *hmm, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, int nseq, int L)
{
  int         i,j;
  float       avg_sc_rnd;
  float       avg_sc_gen;
  ESL_SQ      *sq;
  P7_TRACE    *tr;
  ESL_DSQ     *dsqAA  = NULL;
  ESL_DSQ     *dsqDNA = NULL;
  P7_GMX      *fwd_p  = NULL;
  P7_GMX      *bck_p  = NULL;
  P7_GMX      *fwd    = NULL;
  P7_GMX      *bck    = NULL;
  P7_IVX      *iv     = NULL;
  int         idx;
  float       fsc, bsc;
  float       fsc_p, bsc_p;
  float       nullsc;

  p7_FLogsumInit();
 
  if ((sq     = esl_sq_CreateDigital(abcAA))                                 == NULL)  esl_fatal("sequence creation failed");
  if ((dsqAA  = malloc(sizeof(ESL_DSQ) *((L/3)+2)))                          == NULL)  esl_fatal("malloc failed");
  if ((dsqDNA = malloc(sizeof(ESL_DSQ) *(L+2)))                              == NULL)  esl_fatal("malloc failed");
  if ((fwd_p  = p7_gmx_fs_Create(gm_fs->M, PARSER_ROWS_FWD, L, 0))           == NULL)  esl_fatal("matrix creation failed");
  if ((bck_p  = p7_gmx_fs_Create(gm_fs->M, PARSER_ROWS_BWD, L, 0))           == NULL)  esl_fatal("matrix creation failed");
  if ((fwd    = p7_gmx_fs_Create(gm_fs->M, L,               L, p7P_5CODONS)) == NULL)  esl_fatal("matrix creation failed");
  if ((bck    = p7_gmx_fs_Create(gm_fs->M, L,               L, 0))           == NULL)  esl_fatal("matrix creation failed");
  if ((tr     = p7_trace_Create())                                           == NULL)  esl_fatal("trace creation failed");
  if ((iv     = p7_ivx_Create(gm_fs->M, p7P_5CODONS))                        == NULL)  esl_fatal("ivx creation failed");


  /* Compare Forward and Backward scores when aligneing to DAN sequences reverse translated 
     from randomly generate Amni Acid Sequences. Keep track of the average Forward score */
  avg_sc_rnd = 0.;
  for (idx = 0; idx < nseq; idx++)
    {
      if (esl_rsq_xfIID(r, bgAA->f, abcAA->K, (L/3), dsqAA) != eslOK) esl_fatal("seq generation failed");
      j = 1; 
	  for(i = 1; i <= (L/3); i++) {
        p7_codontable_GetCodon(codon_table, r, dsqAA[i], dsqDNA+j);
        j+=3;
	  }

      if (p7_Forward_Frameshift(dsqDNA, gcode, L, gm_fs, fwd, iv, &fsc)      != eslOK) esl_fatal("forward failed");
      if (p7_Backward_Frameshift(dsqDNA, gcode, L, gm_fs, bck, iv, &bsc)     != eslOK) esl_fatal("backward failed");
      
      if (fabs(fsc-bsc) > 0.001) esl_fatal("Forward/Backward failed: %f %f\n", fsc, bsc);

      if (p7_bg_fs_NullOne(bgAA, dsqAA, (L/3), &nullsc)      != eslOK) esl_fatal("null score failed");
   
      avg_sc_rnd += fsc - nullsc;

      if (p7_ForwardParser_Frameshift_5Codons(dsqDNA, gcode, L, gm_fs, fwd_p, iv, &fsc_p) != eslOK) esl_fatal("forward parser 5 failed");

      if (fabs(fsc-fsc_p) > 0.001) esl_fatal("Forward/Forward Parser failed: %f %f\n", fsc, fsc_p);

      if (p7_BackwardParser_Frameshift_5Codons(dsqDNA, gcode, L, gm_fs, bck_p, iv, &bsc_p) != eslOK) esl_fatal("backward parser 5 failed");

      if (fabs(bsc-bsc_p) > 0.001) esl_fatal("Backward/Backward Parser failed: %f %f\n", bsc, bsc_p);
      if (fabs(fsc_p-bsc_p) > 0.001) esl_fatal("Forward Parser 5 /Backward Parser 5 failed: %f %f\n", fsc_p, bsc_p);

	  p7_gmx_Reuse(fwd_p);
	  p7_gmx_Reuse(bck_p);
      if (p7_ForwardParser_Frameshift_3Codons(dsqDNA, gcode, L, gm_fs, fwd_p, iv, &fsc_p) != eslOK) esl_fatal("forward parser 3 failed");
      if (p7_BackwardParser_Frameshift_3Codons(dsqDNA, gcode, L, gm_fs, bck_p, iv, &bsc_p) != eslOK) esl_fatal("backward parser failed");
       
      if (fabs(fsc_p-bsc_p) > 0.001) esl_fatal("Forward Parser 3 /Backward Parser 3 failed: %f %f\n", fsc_p, bsc_p);

      if (esl_opt_GetBoolean(go, "--vv")) 
        printf("utest_forward_fs: Forward score on random sequence len %d: %.4f (total so far: %.4f)\n", L, fsc, avg_sc_rnd);
    }

  avg_sc_rnd /= (float) nseq;


  /* Get the average forward score on DNA sequence reverse tranlated from Amino Acid 
   * sequence generated by the model and compare to the averge on random sequence */
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

	  p7_gmx_fs_GrowTo(fwd, gm->M, (sq->n*3), (sq->n*3), p7P_5CODONS);

	  if (p7_Forward_Frameshift(dsqDNA, gcode, (sq->n*3), gm_fs, fwd, iv, &fsc)      != eslOK) esl_fatal("forward failed");

	  p7_bg_SetLength(bgAA,  sq->n);

      p7_bg_fs_NullOne(bgAA, dsqAA, sq->n, &nullsc);	  

	  avg_sc_gen += fsc - nullsc;

	  if (esl_opt_GetBoolean(go, "--vv"))
        printf("utest_forward_fs: Forward score on generated sequence len %d: %.4f (total so far: %.4f)\n", (int)(sq->n*3), fsc, avg_sc_gen);
    }

  avg_sc_gen /= (float) nseq;

  if (avg_sc_rnd > avg_sc_gen) esl_fatal("Average Forwrd scores on random sequence (%f) and better than on generated seqeucen (%f)", avg_sc_rnd, avg_sc_gen);

  esl_sq_Destroy(sq);
  p7_gmx_Destroy(fwd_p);
  p7_gmx_Destroy(bck_p);
  p7_gmx_Destroy(fwd);
  p7_gmx_Destroy(bck);
  p7_ivx_Destroy(iv);
  p7_trace_Destroy(tr);
  free(dsqAA);
  free(dsqDNA);
  return;
}



#endif /*p7FWDBACK_FRAMESHIFT_TESTDRIVE*/
/*------------------------- end, unit tests ---------------------*/

/*****************************************************************
 * 4. Test driver.
 *****************************************************************/

/* gcc -g -Wall -Dp7FWDBACK_FRAMESHIFT_TESTDRIVE -I. -I../easel -L. -L../easel -o fwdback_frameshift_utest fwdback_frameshift.c -lhmmer -leasel -lm
 */
#ifdef p7FWDBACK_FRAMESHIFT_TESTDRIVE
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
  P7_FS_PROFILE  *gm_fs  = NULL;
  P7_BG          *bgAA   = NULL;
  P7_BG          *bgDNA  = NULL;
  ESL_GENCODE    *gcode  = NULL;
  P7_CODONTABLE  *ct     = NULL;
  int             M      = 100;
  int             L      = 600; // must by a multiple of 3
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
  if ((gm_fs = p7_profile_fs_Create(hmm->M, abcAA))             == NULL)  esl_fatal("failed to create profile");
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs, L/3, p7_LOCAL) != eslOK) esl_fatal("failed to config profile");
  if (p7_hmm_Validate    (hmm, errbuf, 0.0001)      != eslOK) esl_fatal("whoops, HMM is bad!: %s", errbuf);

  utest_forward_fs    (go, r, abcAA, gcode, ct, bgAA, hmm, gm, gm_fs, nseq, L);

  p7_profile_fs_Destroy(gm_fs);
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
#endif /*p7FWDBACK_FRAMESHIFT_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/



