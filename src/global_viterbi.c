/* Global Viterbi algorithm; both protien to protien and translated 
 * frameshift-aware versions.
 *
 * These implementations are modified to enforce a global alignment
 * necessary for finding optimal splice sites by the 
 * select_splice_option() function of p7_splice.c 
 *
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_gumbel.h"

#include "hmmer.h"

#define IVX(i,k,c) (iv[((k)*p7P_CODONS)+L+3-(i)+(c)])

/* Function:  p7_global_Viterbi()
 * Synopsis:  Global Viterbi algorithm.
 *
 * Purpose:   The Viterbi algorithm for global alignment.
 *
 *            Given a digital sequence <dsq> of length <L>, a profile
 *            <gm>, and DP matrix <gx> allocated for at least <L>
 *            by <gm->M> cells; calculate the maximum scoring path 
 *            that algines all positions in gm and all residues in gm; 
 *            return the Viterbi score in <ret_sc>, and the Viterbi 
 *             matrix is in <gx>.
 *
 *            The caller may then retrieve the Viterbi path by calling
 *            <p7_global_Trace()>.
 *
 *            The Viterbi lod score is returned in nats. The caller
 *            needs to subtract a null model lod score, then convert
 *            to bits.
 *
 * Args:      dsq    - sequence in digitized form, 1..L
 *            L      - length of dsq
 *            gm     - profile.
 *            gx     - DP matrix with room for an MxL alignment
 *            opt_sc - optRETURN: Viterbi lod score in nats
 *
 * Return:   <eslOK> on success.
 */
int
p7_global_Viterbi(const ESL_DSQ *dsq, int L, const P7_PROFILE *gm, P7_GMX *gx, float *opt_sc)
{
  float const *tsc  = gm->tsc;
  float       *rsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  int          M    = gm->M;
  int          i,k;
  float        sc;
  

  /* Initialization of the zero row.  */
  XMX(0,p7G_N) = 0;                                           /* S->N, p=1            */
  XMX(0,p7G_B) = 0;                                           /* S->N->B[0], p=1      */
  XMX(0,p7G_E) = XMX(0,p7G_C) = XMX(0,p7G_J) = -eslINFINITY;  /* need seq to get here */
  for (k = 0; k <= M; k++)
    MMX(0,k) = IMX(0,k) = -eslINFINITY;                       /* need seq to get here */

  i = 0;
  /* global allows zeroth row D state */
  DMX(0,0) = -eslINFINITY;
  

  /* B->D */
  DMX(0,1) = 0.;
  for (k = 2; k <= M; k++)
    DMX(0,k) =  DMX(0,k-1) + TSC(p7P_DD,k-1); 

  /* Initialization of the first row.  */
  rsc = gm->rsc[dsq[1]];

  MMX(1,0) = DMX(1,0) = -eslINFINITY;
 
  /* B->I */
  IMX(1,0) = 0.;
  
  /* B->M */
  MMX(1,1) = MSC(1);

  IMX(1,1) = -eslINFINITY;
  DMX(1,1) = -eslINFINITY;
 
  for (k = 2; k <= M; k++)
  {
    MMX(1,k) = DMX(0,k-1) + TSC(p7P_DM,k-1) + MSC(k);
    IMX(1,k) = -eslINFINITY;
    DMX(1,k) = MMX(1,k-1) + TSC(p7P_MD,k-1);
                       
  }

  XMX(1,p7G_E) = -eslINFINITY;
  XMX(1,p7G_J) = -eslINFINITY;
  XMX(1,p7G_C) = -eslINFINITY;
  XMX(1,p7G_N) = -eslINFINITY;
  XMX(1,p7G_B) = -eslINFINITY;   
  
  /* DP recursion */
  for (i = 2; i <= L; i++)
  {
    rsc = gm->rsc[dsq[i]];

    MMX(i,0) = DMX(i,0) = -eslINFINITY;

    IMX(i,0) = IMX(i-1,0) + TSC(p7P_II,0);

    for (k = 1; k <= M; k++)
    {
      /* match state */
      sc       = ESL_MAX( MMX(i-1,k-1)   + TSC(p7P_MM,k-1),
                          IMX(i-1,k-1)   + TSC(p7P_IM,k-1));
      sc       = ESL_MAX(sc, DMX(i-1,k-1)   + TSC(p7P_DM,k-1));
      MMX(i,k) = sc + MSC(k);

      /* insert state */
      IMX(i,k) = ESL_MAX(MMX(i-1,k) + TSC(p7P_MI,k),
                         IMX(i-1,k) + TSC(p7P_II,k));

      /* delete state */
      DMX(i,k) = ESL_MAX(MMX(i,k-1) + TSC(p7P_MD,k-1),
                         DMX(i,k-1) + TSC(p7P_DD,k-1));
    }

    XMX(i,p7G_E) = -eslINFINITY;
    XMX(i,p7G_J) = -eslINFINITY;
    XMX(i,p7G_C) = -eslINFINITY;
    XMX(i,p7G_N) = -eslINFINITY;
    XMX(i,p7G_B) = -eslINFINITY;
  }

  /* terminal special states for row L */
  XMX(L,p7G_E) = ESL_MAX( MMX(L,M) , 
                 ESL_MAX( IMX(L,M),  DMX(L,M)));

  XMX(L,p7G_C) = XMX(L, p7G_E);
   
  /* T state (not stored) */
  if (opt_sc != NULL) *opt_sc = XMX(L,p7G_C);
  gx->M = M;
  gx->L = L;

 // p7_gmx_Dump(stdout, gx, p7_DEFAULT);
  return eslOK;
}




/* Function:  p7_fs_global_Viterbi()
 * Synopsis:  The global Viterbi algorithm, with frameshift awareness.
 * 
 * Purpose:   The global Viterbi dynamic programming algorithm for aligning 
 *            DNA to proteins with frameshift awareness. 
 *
 *            Given a digital nucleotide sequence <dsq> of length <L>, 
 *            a frameshift aware codon profile <gm_fs>, and DP matrix 
 *            <gx> allocated for at least <L> by <gm_fs->M> cells; 
 *            calculate the maximum scoring global path by 
 *            fremashift-aware Viterbi; return the Viterbi score in 
 *            <ret_sc>, and the Viterbi matrix is in <gx>.
 *            
 *            The caller may then retrieve the Viterbi path by calling
 *            <p7_fs_global_Trace()>.
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
 *            opt_sc - optRETURN: Viterbi lod score in nats
 *           
 * Return:   <eslOK> on success.
 */
int
p7_fs_global_Viterbi(const ESL_DSQ *dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *gm_fs, P7_GMX *gx, float *opt_sc)
{
  float const *tsc  = gm_fs->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  int          M    = gm_fs->M;
  int          i,k,c;
  int          c1, c2, c3, c4, c5;
  int          t, u, v, w, x;
  float       *iv   = NULL;  
  int          status;
 
   /* Allocation and initalization of invermediate value array */
  ESL_ALLOC(iv,  sizeof(float)   * p7P_CODONS * (M+1 + L+1) );

  for (c = 0; c < p7P_CODONS; c++) {
    for(k = 0; k <= M; k++) {
      for(i = 0; i < ESL_MIN(L, 5); i++)
        IVX(i,k,c) = -eslINFINITY;
    }
  }

  /* Initialization of the zero row.  */
  XMX_FS(0,p7G_N) = 0.;                                  /* S->N, p=1            */
  XMX_FS(0,p7G_B) = 0.;                                  /* S->N->B[0], p=1   */
  XMX_FS(0,p7G_E) = XMX_FS(0,p7G_J) = XMX_FS(0,p7G_C) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_FS(0,k,p7G_C0) = MMX_FS(0,k,p7G_C1) = MMX_FS(0,k,p7G_C2) = MMX_FS(0,k,p7G_C3) =
    MMX_FS(0,k,p7G_C4) = MMX_FS(0,k,p7G_C5) = IMX_FS(0,k)        = -eslINFINITY;

  /* global allows zeroth row D state */
  DMX_FS(0,0) = -eslINFINITY;
  DMX_FS(0,1) = 0.; //B->D
  for (k = 2; k <= M; k++)
    DMX_FS(0,k) = DMX_FS(0,k-1) + TSC(p7P_DD,k-1);  
 
  if(esl_abc_XIsCanonical(gcode->nt_abc, dsq[1])) x = dsq[1];
  else                                            x = p7P_MAXCODONS;

  c1 = p7P_CODON1(x);
  
  c1 = p7P_MINIDX(c1, p7P_DEGEN_QC2);

  MMX_FS(1,0,p7G_C0) = MMX_FS(1,0,p7G_C1) = MMX_FS(1,0,p7G_C2) = MMX_FS(1,0,p7G_C3) =
  MMX_FS(1,0,p7G_C4) = MMX_FS(1,0,p7G_C5) =  IMX_FS(1,0) = DMX_FS(1,0)        = -eslINFINITY;

  /* Initialization of the first row.  */
  IVX(1,1,p7P_C1) = 0.; // B->M
 
  MMX_FS(1,1,p7G_C1) = IVX(1,1,p7P_C1) + p7P_MSC_CODON(gm_fs, 1, c1);
  MMX_FS(1,1,p7G_C0) = MMX_FS(1,1,p7G_C1);
  MMX_FS(1,1,p7G_C2) = MMX_FS(1,1,p7G_C3) = MMX_FS(1,1,p7G_C4) = MMX_FS(1,1,p7G_C5) = -eslINFINITY; 

  IMX_FS(1,1) = -eslINFINITY;
  DMX_FS(1,1) = -eslINFINITY;
 
  for (k = 2; k <= M; k++)
  {
    IVX(1,k,p7P_C1) = DMX_FS(0,k-1) + TSC(p7P_DM,k-1);
    MMX_FS(1,k,p7G_C1) = IVX(1,k,p7P_C1) + p7P_MSC_CODON(gm_fs, k, c1);

    MMX_FS(1,k,p7G_C2) = MMX_FS(1,k,p7G_C3) = MMX_FS(1,k,p7G_C4) = MMX_FS(1,k,p7G_C5) = -eslINFINITY;
    MMX_FS(1,k,p7G_C0) = MMX_FS(1,k,p7G_C1);

    IMX_FS(1,k) = -eslINFINITY;
    DMX_FS(1,k) = MMX_FS(1,k-1,p7G_C0) + TSC(p7P_MD,k-1);
  }

  XMX_FS(1,p7G_E) = -eslINFINITY;
  XMX_FS(1,p7G_J) = -eslINFINITY;
  XMX_FS(1,p7G_C) = -eslINFINITY;
  XMX_FS(1,p7G_N) = -eslINFINITY;
  XMX_FS(1,p7G_B) = -eslINFINITY;

  
  /*Special cases for the rows 2-5 */
  t = u = v = w = -1;
  for(i = 2; i <= ESL_MIN(5,L); i++)
  {
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
    
    MMX_FS(i,0,p7G_C0) = MMX_FS(i,0,p7G_C1) = MMX_FS(i,0,p7G_C2) = MMX_FS(i,0,p7G_C3) =
    MMX_FS(i,0,p7G_C4) = MMX_FS(i,0,p7G_C5) = DMX_FS(i,0)        = -eslINFINITY;

    
    /* global allows zeroth column I state */
    if(i == 3) 
       IMX_FS(i,0) = 0.; //B->I
    else
       IMX_FS(i,0) = -eslINFINITY;

    IVX(i,1,p7P_C1) = IMX_FS(i-1,0)          + TSC(p7P_IM,0);

    MMX_FS(i,1,p7G_C1) = IVX(i,1,p7P_C1) + p7P_MSC_CODON(gm_fs, 1, c1);
   
    if (i==2)
      MMX_FS(i,1,p7G_C2) =                   p7P_MSC_CODON(gm_fs, 1, c2); //B->M
    else
      MMX_FS(i,1,p7G_C2) = IVX(i,1,p7P_C2) + p7P_MSC_CODON(gm_fs, 1, c2);
    
    if( i == 3 )
      MMX_FS(i,1,p7G_C3) =                   p7P_MSC_CODON(gm_fs, 1, c3); //B->M     
    else if (i > 3) 
      MMX_FS(i,1,p7G_C3) = IVX(i,1,p7P_C3) + p7P_MSC_CODON(gm_fs, 1, c3);
    else
      MMX_FS(i,1,p7G_C3) = -eslINFINITY;
   
    if(i == 4)
      MMX_FS(i,1,p7G_C4) =                    p7P_MSC_CODON(gm_fs, 1, c4); //B->M
    else if (i > 4)
       MMX_FS(i,1,p7G_C4) = IVX(i,1,p7P_C4) + p7P_MSC_CODON(gm_fs, 1, c4);
    else
      MMX_FS(i,1,p7G_C4) = -eslINFINITY;

    if(i == 5)
      MMX_FS(i,1,p7G_C5) = p7P_MSC_CODON(gm_fs, 1, c5);  //B->M 
    else
      MMX_FS(i,1,p7G_C5) = -eslINFINITY;    

     MMX_FS(i,1,p7G_C0) =  ESL_MAX(ESL_MAX(MMX_FS(i,1,p7G_C1),
                           ESL_MAX(MMX_FS(i,1,p7G_C2), MMX_FS(i,1,p7G_C3))),
                           ESL_MAX(MMX_FS(i,1,p7G_C4), MMX_FS(i,1,p7G_C5)));

    /* insert state */
    if ( i > 2)
      IMX_FS(i,1) = ESL_MAX(MMX_FS(i-3,1,p7G_C0) + TSC(p7P_MI,1),
                            IMX_FS(i-3,1)        + TSC(p7P_II,1));
    else
      IMX_FS(i,1) = -eslINFINITY;

    /* delete state */
    DMX_FS(i,1) =  -eslINFINITY; 
                     

    /* Initialization of the states reacheable at row i */
    for (k = 2; k <= M; k++)
    {
      IVX(i,k,p7P_C1) = ESL_MAX(MMX_FS(i-1,k-1,p7G_C0)   + TSC(p7P_MM,k-1),
                        ESL_MAX(IMX_FS(i-1,k-1)          + TSC(p7P_IM,k-1),
                                DMX_FS(i-1,k-1)          + TSC(p7P_DM,k-1)));

      MMX_FS(i,k,p7G_C1) = IVX(i,k,p7P_C1) + p7P_MSC_CODON(gm_fs, k, c1);
      MMX_FS(i,k,p7G_C2) = IVX(i,k,p7P_C2) + p7P_MSC_CODON(gm_fs, k, c2);
        
      if( i > 2 )
        MMX_FS(i,k,p7G_C3) = IVX(i,k,p7P_C3) + p7P_MSC_CODON(gm_fs, k, c3);
      else
        MMX_FS(i,k,p7G_C3) = -eslINFINITY;
      if( i > 3 )
        MMX_FS(i,k,p7G_C4) = IVX(i,k,p7P_C4) + p7P_MSC_CODON(gm_fs, k, c4);
      else
        MMX_FS(i,k,p7G_C4) = -eslINFINITY;

      if(i == 5)
         MMX_FS(i,k,p7G_C5) = IVX(i,k,p7P_C5) + p7P_MSC_CODON(gm_fs, k, c5);
      else
        MMX_FS(i,k,p7G_C5) = -eslINFINITY;

      MMX_FS(i,k,p7G_C0) =  ESL_MAX(ESL_MAX(MMX_FS(i,k,p7G_C1),
                            ESL_MAX(MMX_FS(i,k,p7G_C2), MMX_FS(i,k,p7G_C3))),
                            ESL_MAX(MMX_FS(i,k,p7G_C4), MMX_FS(i,k,p7G_C5)));

      /* insert state */
      if ( i > 2)
        IMX_FS(i,k) = ESL_MAX(MMX_FS(i-3,k,p7G_C0) + TSC(p7P_MI,k),
                              IMX_FS(i-3,k)        + TSC(p7P_II,k));
      else
        IMX_FS(i,k) = -eslINFINITY;

      /* delete state */
      DMX_FS(i,k) = ESL_MAX(MMX_FS(i,k-1,p7G_C0) + TSC(p7P_MD,k-1),
                            DMX_FS(i,k-1)        + TSC(p7P_DD,k-1));
    }

    XMX_FS(i,p7G_E) = -eslINFINITY;
    XMX_FS(i,p7G_J) = -eslINFINITY; 
    XMX_FS(i,p7G_C) = -eslINFINITY;
    XMX_FS(i,p7G_N) = -eslINFINITY;
    XMX_FS(i,p7G_B) = -eslINFINITY;
  } 

  /* Main DP recursion */
  for (i = 6; i <= L; i++) 
  {
    MMX_FS(i,0,p7G_C0) = MMX_FS(i,0,p7G_C1) = MMX_FS(i,0,p7G_C2) = MMX_FS(i,0,p7G_C3)
    = MMX_FS(i,0,p7G_C4) = MMX_FS(i,0,p7G_C5) = DMX_FS(i,0) = -eslINFINITY;

    IMX_FS(i,0) = IMX_FS(i-3,0) + TSC(p7P_II,0);
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

    for (k = 1; k <= M; k++)
    {

      IVX(i,k,p7P_C1) = ESL_MAX(MMX_FS(i-1,k-1,p7G_C0)   + TSC(p7P_MM,k-1),
                        ESL_MAX(IMX_FS(i-1,k-1)          + TSC(p7P_IM,k-1),
                                DMX_FS(i-1,k-1)          + TSC(p7P_DM,k-1)));

      MMX_FS(i,k,p7G_C1) = IVX(i,k,p7P_C1) + p7P_MSC_CODON(gm_fs, k, c1);

      MMX_FS(i,k,p7G_C2) = IVX(i,k,p7P_C2) + p7P_MSC_CODON(gm_fs, k, c2);

      MMX_FS(i,k,p7G_C3) = IVX(i,k,p7P_C3) + p7P_MSC_CODON(gm_fs, k, c3);

      MMX_FS(i,k,p7G_C4) = IVX(i,k,p7P_C4) + p7P_MSC_CODON(gm_fs, k, c4);

      MMX_FS(i,k,p7G_C5) = IVX(i,k,p7P_C5) + p7P_MSC_CODON(gm_fs, k, c5);

      MMX_FS(i,k,p7G_C0) =  ESL_MAX(ESL_MAX(MMX_FS(i,k,p7G_C1),
                            ESL_MAX(MMX_FS(i,k,p7G_C2), MMX_FS(i,k,p7G_C3))),
                            ESL_MAX(MMX_FS(i,k,p7G_C4), MMX_FS(i,k,p7G_C5)));

      /* insert state */
      IMX_FS(i,k) = ESL_MAX(MMX_FS(i-3,k,p7G_C0) + TSC(p7P_MI,k),
                            IMX_FS(i-3,k)        + TSC(p7P_II,k));

      /* delete state */
      DMX_FS(i,k) = ESL_MAX(MMX_FS(i,k-1,p7G_C0) + TSC(p7P_MD,k-1),
                            DMX_FS(i,k-1)        + TSC(p7P_DD,k-1));

    }
  
    XMX_FS(i,p7G_E) = -eslINFINITY;
    XMX_FS(i,p7G_J) = -eslINFINITY;
    XMX_FS(i,p7G_C) = -eslINFINITY;
    XMX_FS(i,p7G_N) = -eslINFINITY;
    XMX_FS(i,p7G_B) = -eslINFINITY;
  }

  /* change terminal special states for row L */
  XMX_FS(L,p7G_E) = ESL_MAX( MMX_FS(L,M,p7G_C0),
                    ESL_MAX( IMX_FS(L,M), DMX_FS(L,M)));

  XMX(L,p7G_C) = XMX_FS(L,p7G_E);

  /* T state (not stored) */
  if (opt_sc != NULL) *opt_sc = XMX_FS(L,p7G_C);
                                
  gx->M = M;
  gx->L = L;

  if (iv != NULL) free(iv);
  //p7_gmx_fs_Dump(stdout, gx, p7_DEFAULT, FALSE);
  return eslOK;

  ERROR:
  if (iv != NULL) free(iv);
  return status;
}


int
p7_fs_semiglobal_Viterbi(const ESL_DSQ *sub_dsq, const ESL_GENCODE *gcode, int L, const P7_FS_PROFILE *sub_gm, P7_GMX *gx)
{
  float const *tsc  = sub_gm->tsc;
  float      **dp   = gx->dp;
  float       *xmx  = gx->xmx;
  int          M    = sub_gm->M;
  int          i,k,c;
  int          c1, c2, c3, c4, c5;
  int          t, u, v, w, x;
  float       *iv   = NULL;
  int          status;

  /* Allocation and initalization of intermediate value array */
  ESL_ALLOC(iv,  sizeof(float)   * p7P_CODONS * (M+1 + L+1) );
  for (c = 0; c < p7P_CODONS; c++) {
    for(k = 0; k <= M; k++)
      IVX(5,k,c) = -eslINFINITY;
  }

  /* Initialization of the zero row.  */
  XMX(0,p7G_N) = 0.;                                  /* S->N, p=1            */
  XMX(0,p7G_B) = sub_gm->xsc[p7P_N][p7P_MOVE];            /* S->N->B, no N-tail   */ 
  XMX(0,p7G_E) = XMX(0,p7G_J) = XMX(0,p7G_C) = -eslINFINITY;
  for (k = 0; k <= M; k++)
    MMX_SP(0,k) = IMX_SP(0,k) = -eslINFINITY;

  DMX_SP(0,0) = -eslINFINITY;
  DMX_SP(0,1) = XMX(0,p7G_B) + TSC(p7P_BM,0); 
  for (k = 2; k <= M; k++)
    DMX_SP(0,k) = DMX_SP(0,k-1) + TSC(p7P_DD,k-1);

  /*Special cases for the first 4 rows */
  t = u = v = w = x = -1;
  for(i = 1; i <= 4; i++)
  {
    u = v;
    v = w;
    w = x;

    if( i < 3 ) 
      XMX(i,p7G_N) = 0.;
    else 
      XMX(i,p7G_N) =  XMX(i-3,p7G_N) + sub_gm->xsc[p7P_N][p7P_LOOP];
    XMX(i,p7G_B) =  XMX(i,p7G_N)   + sub_gm->xsc[p7P_N][p7P_MOVE]; 
 
    /* if new nucleotide is not A,C,G, or T set it to placeholder vlaue */
    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
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

    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = -eslINFINITY;
 
    IVX(i,1,p7P_C1) = XMX(i-1,p7G_B) + TSC(p7P_BM,0);
  
    MMX_SP(i,1) = IVX(i,1,p7P_C1) + p7P_MSC_CODON(sub_gm, 1, c1);
    if( i > 1 ) MMX_SP(i,1) = ESL_MAX(MMX_SP(i,1), IVX(i,1,p7P_C2) + p7P_MSC_CODON(sub_gm, 1, c2));
    if( i > 2 ) MMX_SP(i,1) = ESL_MAX(MMX_SP(i,1), IVX(i,1,p7P_C3) + p7P_MSC_CODON(sub_gm, 1, c3));
    if( i > 3 ) MMX_SP(i,1) = ESL_MAX(MMX_SP(i,1), IVX(i,1,p7P_C4) + p7P_MSC_CODON(sub_gm, 1, c4));                    
    
    if ( i > 2) {
     if(p7P_MSC_CODON(sub_gm, 1, c3) == -eslINFINITY) IMX_SP(i,1) = -eslINFINITY;
     else                                             IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,1),
                                                                            IMX_SP(i-3,1) + TSC(p7P_II,1));
    }
    else        IMX_SP(i,1) = -eslINFINITY;                  

    DMX_SP(i,1) = XMX(i,p7G_B) + TSC(p7P_BM,0);
   
    /* Initialization of the states reacheable at row i */
    for (k = 2; k < M; k++)
    {
      IVX(i,k,p7P_C1) = ESL_MAX(MMX_SP(i-1,k-1) + TSC(p7P_MM,k-1),
                        ESL_MAX(IMX_SP(i-1,k-1) + TSC(p7P_IM,k-1),
                                DMX_SP(i-1,k-1) + TSC(p7P_DM,k-1)));

      MMX_SP(i,k) = IVX(i,k,p7P_C1) + p7P_MSC_CODON(sub_gm, k, c1);   
      if( i > 1 ) MMX_SP(i,k) = ESL_MAX(MMX_SP(i,k), IVX(i,k,p7P_C2) + p7P_MSC_CODON(sub_gm, k, c2));
      if( i > 2 ) MMX_SP(i,k) = ESL_MAX(MMX_SP(i,k), IVX(i,k,p7P_C3) + p7P_MSC_CODON(sub_gm, k, c3));
      if( i > 3 ) MMX_SP(i,k) = ESL_MAX(MMX_SP(i,k), IVX(i,k,p7P_C4) + p7P_MSC_CODON(sub_gm, k, c4));

      if ( i > 2) { 
       if(p7P_MSC_CODON(sub_gm, k, c3) == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;
       else                                             IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,k),
                                                                              IMX_SP(i-3,k) + TSC(p7P_II,k));
      }
      else        IMX_SP(i,k) = -eslINFINITY;

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,k-1));         
    } 
  
    IVX(i,M,p7P_C1) = ESL_MAX(MMX_SP(i-1,M-1) + TSC(p7P_MM,M-1),
                      ESL_MAX(IMX_SP(i-1,M-1) + TSC(p7P_IM,M-1),
                              DMX_SP(i-1,M-1) + TSC(p7P_DM,k=M-1))); 

    MMX_SP(i,M) = IVX(i,M,p7P_C1) + p7P_MSC_CODON(sub_gm, M, c1);
    if( i > 1 ) MMX_SP(i,M) = ESL_MAX(MMX_SP(i,M), IVX(i,M,p7P_C2) + p7P_MSC_CODON(sub_gm, M, c2));
    if( i > 2 ) MMX_SP(i,M) = ESL_MAX(MMX_SP(i,M), IVX(i,M,p7P_C3) + p7P_MSC_CODON(sub_gm, M, c3));
    if( i > 3 ) MMX_SP(i,M) = ESL_MAX(MMX_SP(i,M), IVX(i,M,p7P_C4) + p7P_MSC_CODON(sub_gm, M, c4));

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,M-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,M-1));

    XMX(i,p7G_E) = ESL_MAX(MMX_SP(i,M), DMX_SP(i,M));
    
    if (i < 3)
       XMX(i,p7G_C) = XMX(i,p7G_E)   + sub_gm->xsc[p7P_E][p7P_MOVE];
    else
      XMX(i,p7G_C) = ESL_MAX(XMX(i-3,p7G_C) + sub_gm->xsc[p7P_C][p7P_LOOP],
                             XMX(i,p7G_E)   + sub_gm->xsc[p7P_E][p7P_MOVE]);

    XMX(i,p7G_J) = -eslINFINITY;

  }

    /* Main DP recursion */
  for (i = 5; i <= L; i++) 
  {
    t = u;
    u = v;
    v = w;
    w = x;

    if(esl_abc_XIsCanonical(gcode->nt_abc, sub_dsq[i])) x = sub_dsq[i];
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

    XMX(i,p7G_N) = XMX(i-3,p7G_N) + sub_gm->xsc[p7P_N][p7P_LOOP];
    XMX(i,p7G_B) = XMX(i,p7G_N)   + sub_gm->xsc[p7P_N][p7P_MOVE];
    
    MMX_SP(i,0) = IMX_SP(i,0) = DMX_SP(i,0) = -eslINFINITY;

    IVX(i,1,p7P_C1) = XMX(i-1,p7G_B) + TSC(p7P_BM,0);

    MMX_SP(i,1) = IVX(i,1,p7P_C1) + p7P_MSC_CODON(sub_gm, 1, c1);
    MMX_SP(i,1) = ESL_MAX(MMX_SP(i,1), IVX(i,1,p7P_C2) + p7P_MSC_CODON(sub_gm, 1, c2));
    MMX_SP(i,1) = ESL_MAX(MMX_SP(i,1), IVX(i,1,p7P_C3) + p7P_MSC_CODON(sub_gm, 1, c3));
    MMX_SP(i,1) = ESL_MAX(MMX_SP(i,1), IVX(i,1,p7P_C4) + p7P_MSC_CODON(sub_gm, 1, c4));
    MMX_SP(i,1) = ESL_MAX(MMX_SP(i,1), IVX(i,1,p7P_C5) + p7P_MSC_CODON(sub_gm, 1, c5));

    if(p7P_MSC_CODON(sub_gm, 1, c3) == -eslINFINITY) IMX_SP(i,1) = -eslINFINITY;
    else                                             IMX_SP(i,1) = ESL_MAX(MMX_SP(i-3,1) + TSC(p7P_MI,1),
                                                                           IMX_SP(i-3,1) + TSC(p7P_II,1));

    DMX_SP(i,1) = XMX(i,p7G_B) + TSC(p7P_BM,0); 
    
    for (k = 2; k < M; k++) {

      IVX(i,k,p7P_C1) = ESL_MAX(MMX_SP(i-1,k-1) + TSC(p7P_MM,k-1),
                        ESL_MAX(IMX_SP(i-1,k-1) + TSC(p7P_IM,k-1),
                                DMX_SP(i-1,k-1) + TSC(p7P_DM,k-1)));

      MMX_SP(i,k) = IVX(i,k,p7P_C1) + p7P_MSC_CODON(sub_gm, k, c1);
      MMX_SP(i,k) = ESL_MAX(MMX_SP(i,k), IVX(i,k,p7P_C2) + p7P_MSC_CODON(sub_gm, k, c2));
      MMX_SP(i,k) = ESL_MAX(MMX_SP(i,k), IVX(i,k,p7P_C3) + p7P_MSC_CODON(sub_gm, k, c3));
      MMX_SP(i,k) = ESL_MAX(MMX_SP(i,k), IVX(i,k,p7P_C4) + p7P_MSC_CODON(sub_gm, k, c4));
      MMX_SP(i,k) = ESL_MAX(MMX_SP(i,k), IVX(i,k,p7P_C5) + p7P_MSC_CODON(sub_gm, k, c5));

      if(p7P_MSC_CODON(sub_gm, k, c3) == -eslINFINITY) IMX_SP(i,k) = -eslINFINITY;
      else                                             IMX_SP(i,k) = ESL_MAX(MMX_SP(i-3,k) + TSC(p7P_MI,k),
                                                       IMX_SP(i-3,k) + TSC(p7P_II,k));

      DMX_SP(i,k) = ESL_MAX(MMX_SP(i,k-1) + TSC(p7P_MD,k-1),
                            DMX_SP(i,k-1) + TSC(p7P_DD,k-1));
   
    }

    IVX(i,M,p7P_C1) = ESL_MAX(MMX_SP(i-1,M-1) + TSC(p7P_MM,M-1),
                      ESL_MAX(IMX_SP(i-1,M-1) + TSC(p7P_IM,M-1),
                              DMX_SP(i-1,M-1) + TSC(p7P_DM,M-1)));

    MMX_SP(i,M) = IVX(i,M,p7P_C1) + p7P_MSC_CODON(sub_gm, M, c1);
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i,M), IVX(i,M,p7P_C2) + p7P_MSC_CODON(sub_gm, M, c2));
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i,M), IVX(i,M,p7P_C3) + p7P_MSC_CODON(sub_gm, M, c3));
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i,M), IVX(i,M,p7P_C4) + p7P_MSC_CODON(sub_gm, M, c4));
    MMX_SP(i,M) = ESL_MAX(MMX_SP(i,M), IVX(i,M,p7P_C5) + p7P_MSC_CODON(sub_gm, M, c5));

    IMX_SP(i,M) = -eslINFINITY;

    DMX_SP(i,M) = ESL_MAX(MMX_SP(i,M-1) + TSC(p7P_MD,M-1),
                          DMX_SP(i,M-1) + TSC(p7P_DD,M-1));

    XMX(i,p7G_E) = ESL_MAX(MMX_SP(i,M), DMX_SP(i,M));
    XMX(i,p7G_C) = ESL_MAX(XMX(i-3,p7G_C) + sub_gm->xsc[p7P_C][p7P_LOOP],
                           XMX(i,p7G_E)   + sub_gm->xsc[p7P_E][p7P_MOVE]);

    XMX(i,p7G_J) = -eslINFINITY;    

  }
  
  gx->M = sub_gm->M;
  gx->L = L;

  if (iv != NULL) free(iv);
  return eslOK;

  ERROR:
    if (iv != NULL) free(iv);
    return status;    
}

/*-------------------- end, viterbi -----------------------------*/




