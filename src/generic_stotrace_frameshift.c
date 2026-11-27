/* Stochastic traceback of a checkpointed Forward matrix; frameshift aware version.
 * 
 * Contents:
 *   1. Stochastic traceback implementation.
 *   2. Benchmark driver.
 *
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "hmmer.h"

/*****************************************************************
 * 1. Stochastic traceback implementation.
 *****************************************************************/

/* nt_at()
 * The nucleotide at position <j> of <dsq> for codon lookup: the
 * placeholder <p7P_MAXCODONS5> before position 1 or for anything but
 * A,C,G,T.
 */
static inline int
nt_at(const ESL_DSQ *dsq, int j)
{
  return (j >= 1 && dsq[j] < p7P_MAXNUC) ? dsq[j] : p7P_MAXCODONS5;
}

/* select_m_joint_fs()
 * For M at row <i>, model position <k>: choose its codon length c and the
 * state before it together, from the 20 terms the Forward sums for this
 * cell: X(i-c, k-1) + t(X->M_k) + e_c(k), for X in B, M, I, D and c = 1..5
 * (c <= i), where e_c is the emission score of the codon of length c that
 * ends at <i>, as p7_GForwardRow_Frameshift_5Codons() sums them. Sets <*ret_c>
 * and <*ret_prv>. <sc> is workspace for 20 scores.
 */
static void
select_m_joint_fs(ESL_RANDOMNESS *r, const P7_FS_PROFILE *gm_fs, float * const *dp, const float *xmx, const ESL_DSQ *dsq,
                  float *sc, int i, int k, int *ret_c, int *ret_prv)
{
  float const *tsc = gm_fs->tsc;
  int          st[4] = { p7T_B, p7T_M, p7T_I, p7T_D };
  int          cidx[6];
  int          t, u, v, w, x;
  int          c, p, z;
  float        e;

  x = nt_at(dsq, i);   w = nt_at(dsq, i-1); v = nt_at(dsq, i-2);
  u = nt_at(dsq, i-3); t = nt_at(dsq, i-4);
  cidx[1] = p7P_MINIDX(p7P_CODON1_FS5(x),             p7P_DEGEN5_QC2);
  cidx[2] = p7P_MINIDX(p7P_CODON2_FS5(w, x),          p7P_DEGEN5_QC1);
  cidx[3] = p7P_MINIDX(p7P_CODON3_FS5(v, w, x),       p7P_DEGEN5_C);
  cidx[4] = p7P_MINIDX(p7P_CODON4_FS5(u, v, w, x),    p7P_DEGEN5_QC1);
  cidx[5] = p7P_MINIDX(p7P_CODON5_FS5(t, u, v, w, x), p7P_DEGEN5_QC2);

  for (c = 1; c <= 5; c++)
    {
      p = i - c;
      for (z = 0; z < 4; z++) sc[(c-1)*4 + z] = -eslINFINITY;
      if (p < 0) continue;
      e = gm_fs->rsc[cidx[c]][k];
      sc[(c-1)*4 + 0] = XMX_FS(p,p7G_B)      + TSC(p7P_BM,k-1) + e;
      sc[(c-1)*4 + 1] = MMX(p,k-1) + TSC(p7P_MM,k-1) + e;
      sc[(c-1)*4 + 2] = IMX(p,k-1)        + TSC(p7P_IM,k-1) + e;
      sc[(c-1)*4 + 3] = DMX(p,k-1)        + TSC(p7P_DM,k-1) + e;
    }
  esl_vec_FLogNorm(sc, 20);
  z = esl_rnd_FChoose(r, sc, 20);
  *ret_c   = z / 4 + 1;
  *ret_prv = st[z % 4];
}


/* gsto_fs_step()
 * One step of a stochastic trace back through the frameshift Forward
 * 3-cell rows <dp> (D, I and M, the C0 sum, of a checkpointed Forward
 * matrix) and special states <xmx> of the sequence <dsq>: from state
 * <*ret_sprv> at row <*ret_i> and model position <*ret_k>, choose the
 * state before it, append it to <tr>, and update <*ret_sprv>, <*ret_i>
 * and <*ret_k>. An M state is appended with codon length 0; the step
 * from it chooses its codon length and the state before it together
 * (select_m_joint_fs()) and sets the codon length in <tr>. <sc> is
 * workspace for 2M+1 (and at least 20) scores. Reads the core cells of
 * rows <*ret_i>-5..<*ret_i>.
 */
static int
gsto_fs_step(ESL_RANDOMNESS *r, const P7_FS_PROFILE *gm_fs, float * const *dp, const float *xmx, const ESL_DSQ *dsq, float *sc, P7_TRACE *tr,
             int *ret_sprv, int *ret_i, int *ret_k)
{
  float const *tsc  = gm_fs->tsc;
  int          M    = gm_fs->M;
  int          i    = *ret_i;
  int          k    = *ret_k;
  int          sprv = *ret_sprv;
  int          scur = p7T_BOGUS;
  int          c, z;
  int          status;

  switch (sprv) {
    /* C(i) comes from C(i-3) or E(i). On the first C step, from row L,
     * the alignment ends at row L, L-1 or L-2 (as the Forward score sums
     * C(L) + C(L-1) tCC + C(L-2) tCC); the nucleotides after it are
     * appended as C. A C loop appends its three nucleotides as C. */
    case p7T_C:
      if (tr->N == 2)
        {
          sc[0] = XMX_FS(i, p7G_C);
          sc[1] = (i >= 1) ? XMX_FS(i-1, p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP] : -eslINFINITY;
          sc[2] = (i >= 2) ? XMX_FS(i-2, p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP] : -eslINFINITY;
          esl_vec_FLogNorm(sc, 3);
          for (z = esl_rnd_FChoose(r, sc, 3); z > 0; z--)
            { if ((status = p7_trace_fs_Append(tr, p7T_C, k, i, 0)) != eslOK) return status; i--; }
        }
      if (XMX_FS(i,p7G_C) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible C reached at i=%d", i);

      sc[0] = (i >= 3) ? XMX_FS(i-3, p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP] : -eslINFINITY;
      sc[1] = XMX_FS(i, p7G_E) + gm_fs->xsc[p7P_E][p7P_MOVE];
      esl_vec_FLogNorm(sc, 2);
      if (esl_rnd_FChoose(r, sc, 2) == 0)
        {
          for (z = 0; z < 2; z++)
            { if ((status = p7_trace_fs_Append(tr, p7T_C, k, i, 0)) != eslOK) return status; i--; }
          scur = p7T_C;     /* the third nucleotide is appended below */
        }
      else scur = p7T_E;
      break;

    /* E connects from any M or D state. k set here */
    case p7T_E:  
      if (XMX_FS(i, p7G_E) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
    if (p7_fs_profile_IsLocal(gm_fs)) 
    { /* local models come from any M, D */
      sc[0] = sc[M+1] = -eslINFINITY;
      for (k = 1; k <= M; k++) sc[k]   = MMX(i,k);
      for (k = 2; k <= M; k++) sc[k+M] = DMX(i,k);
      esl_vec_FLogNorm(sc, 2*M+1); /* now sc is a prob vector */
      k = esl_rnd_FChoose(r, sc, 2*M+1);
      if (k <= M)    scur = p7T_M;
      else { k -= M; scur = p7T_D; }
    } 
    else 
    {     /* glocal models come from M_M or D_M  */
      k     = M;
      sc[0] = MMX(i,M);
      sc[1] = DMX(i,M);
      esl_vec_FLogNorm(sc, 2); /* now sc is a prob vector */
      scur = (esl_rnd_FChoose(r, sc, 2) == 0) ? p7T_M : p7T_D;
    }
    break;

    /* M at row i: its codon length c and the state before it, at row
     * i-c and k-1, are chosen together */
    case p7T_M:
     select_m_joint_fs(r, gm_fs, dp, xmx, dsq, sc, i, k, &c, &scur);
     tr->c[tr->N-1] = c;
     i -= c;
     k--; 
     break;

    /* D connects from M,D at i,k-1 */
    case p7T_D:
      if (DMX(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);
      sc[0] = MMX(i,k-1) + TSC(p7P_MD, k-1);
      sc[1] = DMX(i,k-1) + TSC(p7P_DD, k-1);
      esl_vec_FLogNorm(sc, 2); 
      scur = (esl_rnd_FChoose(r, sc, 2) == 0) ? p7T_M : p7T_D;
      k--;
      break;

    /* I connects from M,I at i-1,k */
    case p7T_I:
      if (IMX(i,k) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);
      sc[0] = MMX(i-3,k) + TSC(p7P_MI, k);
      sc[1] = IMX(i-3,k) + TSC(p7P_II, k);

      esl_vec_FLogNorm(sc, 2); 
      scur = (esl_rnd_FChoose(r, sc, 2) == 0) ? p7T_M : p7T_I;
      i-=3;
      break;

    /* N connects from S, N */
    case p7T_N:
      if (XMX_FS(i, p7G_N) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
      scur = (i == 0) ? p7T_S : p7T_N;
      break;

    /* B connects from N, J */
    case p7T_B:    
      if (XMX_FS(i,p7G_B) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

      sc[0] = XMX_FS(i, p7G_N) + gm_fs->xsc[p7P_N][p7P_MOVE];
      sc[1] = XMX_FS(i, p7G_J) + gm_fs->xsc[p7P_J][p7P_MOVE];
      esl_vec_FLogNorm(sc, 2); 
      scur = (esl_rnd_FChoose(r, sc, 2) == 0) ? p7T_N : p7T_J;
      break;

    /* J(i) comes from J(i-3) or E(i); a J loop appends its three nucleotides as J */
    case p7T_J:  
      if (XMX_FS(i,p7G_J) == -eslINFINITY) ESL_EXCEPTION(eslFAIL, "impossible J reached at i=%d", i);

      sc[0] = (i >= 3) ? XMX_FS(i-3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP] : -eslINFINITY;
      sc[1] = XMX_FS(i,  p7G_E) + gm_fs->xsc[p7P_E][p7P_LOOP];
      esl_vec_FLogNorm(sc, 2);
      if (esl_rnd_FChoose(r, sc, 2) == 0)
        {
          for (z = 0; z < 2; z++)
            { if ((status = p7_trace_fs_Append(tr, p7T_J, k, i, 0)) != eslOK) return status; i--; }
          scur = p7T_J;     /* the third nucleotide is appended below */
        }
      else scur = p7T_E;
      break;
    default: ESL_EXCEPTION(eslFAIL, "bogus state in traceback");
    } /* end switch over statetype[tpos-1] */

    /* Append this state and the current i,k to be explained to the growing trace */
    if ((status = p7_trace_fs_Append(tr, scur, k, i, 0)) != eslOK) return status;

    /* For NCJ, we had to defer i decrement. */
    if ( (scur == p7T_N || scur == p7T_C || scur == p7T_J) && scur == sprv) i--;

  sprv = scur;

  *ret_sprv = sprv;
  *ret_i    = i;
  *ret_k    = k;
  return eslOK;
}


/* Function:  p7_GStochasticTraceChk_Frameshift()
 * Synopsis:  Sample tracebacks from a checkpointed frameshift Forward matrix.
 *
 * Purpose:   Sample <ntr> stochastic tracebacks of the 5-codon profile
 *            <gm_fs5> against <dsq> of length <L> from the checkpointed
 *            Forward matrix <cx> that <p7_GForwardChk_Frameshift_5Codons()>
 *            filled: trace <t> into <tr[t]>, with its own random number
 *            generator <r[t]>.
 *
 *            All traces are sampled together. The blocks of <cx> are
 *            loaded from the last to the first, each computed again from
 *            its checkpoint by <p7_GForwardChkLoad_Frameshift_5Codons()> with
 *            <iv> as workspace, and every trace is taken back through the
 *            rows of the loaded block before the next block is loaded, so
 *            the Forward matrix is computed again only once for all the
 *            traces. On return block 0 of <cx> is loaded.
 *
 * Args:      r      - <ntr> random number generators, one per trace
 *            ntr    - number of traces
 *            dsq    - digital sequence, 1..L
 *            L      - length of dsq
 *            gm_fs5 - 5-codon frameshift profile
 *            cx     - checkpointed Forward matrix from p7_GForwardChk_Frameshift_5Codons()
 *            iv     - intermediate value matrix, workspace
 *            tr     - <ntr> traces; RETURN: the sampled traces
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if <cx> isn't a
 *            checkpointed frameshift Forward matrix for this comparison;
 *            <eslFAIL> if a trace reaches an impossible state.
 */
int
p7_GStochasticTraceChk_Frameshift(ESL_RANDOMNESS **r, int ntr, const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs5,
                                  P7_CHKGMX *cx, P7_IVX *iv, P7_TRACE **tr)
{
  float *sc   = NULL;
  int   *sprv = NULL;
  int   *i    = NULL;
  int   *k    = NULL;
  int    b, t;
  int    status;

  if (cx->dir != p7_CHK_FORWARD || cx->nscells != p7G_NSCELLS || cx->L != L || cx->M != gm_fs5->M)
    ESL_EXCEPTION(eslEINVAL, "checkpointed frameshift Forward matrix doesn't match");

  ESL_ALLOC(sc,   sizeof(float) * ESL_MAX(2*gm_fs5->M+1, 20));
  ESL_ALLOC(sprv, sizeof(int)   * ESL_MAX(1, ntr));
  ESL_ALLOC(i,    sizeof(int)   * ESL_MAX(1, ntr));
  ESL_ALLOC(k,    sizeof(int)   * ESL_MAX(1, ntr));
  for (t = 0; t < ntr; t++)
    {
      i[t] = L;
      k[t] = 0;
      if ((status = p7_trace_fs_Append(tr[t], p7T_T, 0, L, 0)) != eslOK) goto ERROR;
      if ((status = p7_trace_fs_Append(tr[t], p7T_C, 0, L, 0)) != eslOK) goto ERROR;
      sprv[t] = p7T_C;
    }

  /* a step from row i reads rows i-5..i, which are loaded with the
   * block holding row i and the 5 rows saved before it */
  for (b = (L / cx->blockR) * cx->blockR; b >= 0; b -= cx->blockR)
    {
      if (b != cx->b && (status = p7_GForwardChkLoad_Frameshift_5Codons(dsq, gm_fs5, cx, iv, b)) != eslOK) goto ERROR;
      for (t = 0; t < ntr; t++)
        while (sprv[t] != p7T_S && i[t] >= b)
          if ((status = gsto_fs_step(r[t], gm_fs5, cx->dp, cx->xmx, dsq, sc, tr[t], &sprv[t], &i[t], &k[t])) != eslOK) goto ERROR;
    }

  for (t = 0; t < ntr; t++)
    {
      if ((status = p7_trace_fs_Reverse(tr[t])) != eslOK) goto ERROR;
      tr[t]->M = gm_fs5->M;
      tr[t]->L = L;
    }
  free(sc);
  free(sprv);
  free(i);
  free(k);
  return eslOK;

 ERROR:
  free(sc);
  free(sprv);
  free(i);
  free(k);
  return status;
}



/*------------------- end, stochastic trace ---------------------*/





/*****************************************************************
 * 2. Benchmark driver
 *****************************************************************/
#ifdef p7GENERIC_STOTRACE_FRAMESHIFT_BENCHMARK
/*
   gcc -g -O2      -o generic_stotrace_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_STOTRACE_FRAMESHIFT_BENCHMARK generic_stotrace_frameshift.c -lhmmer -leasel -lm
   icc -O3 -static -o generic_stotrace_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7GENERIC_STOTRACE_FRAMESHIFT_BENCHMARK generic_stotrace_frameshift.c -lhmmer -leasel -lm
   ./generic_stotrace_frameshift_benchmark <hmmfile>
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
  { "-L",        eslARG_INT,   "1200", NULL, "n>0", NULL,  NULL, NULL, "length of random target seq" ,                   0 },
  { "-N",        eslARG_INT,  "50000", NULL, "n>0", NULL,  NULL, NULL, "number of sampled tracebacks",                   0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for generic stochastic trace";

int 
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = p7_CreateDefaultApp(options, 1, argc, argv, banner, usage);
  char           *hmmfile = esl_opt_GetArg(go, 1);
  ESL_STOPWATCH  *w       = esl_stopwatch_Create();
  ESL_RANDOMNESS *r       = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA   = NULL;
  ESL_ALPHABET   *abcDNA  = NULL;
  ESL_GENCODE    *gcode   = NULL;
  P7_HMMFILE     *hfp     = NULL;
  P7_HMM         *hmm     = NULL;
  P7_BG          *bgAA    = NULL;
  P7_BG          *bgDNA   = NULL;
  P7_FS_PROFILE  *gm_fs5   = NULL;
  P7_CHKGMX      *cx      = NULL;
  P7_IVX         *iv      = NULL;
  ESL_RANDOMNESS **rng    = NULL;
  P7_TRACE      **tr      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             t;
  float           sc, fsc;
  float           bestsc  = -eslINFINITY;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)          != eslOK) p7_Fail("Failed to read HMM");

  abcDNA = esl_alphabet_Create(eslDNA);
  gcode  = esl_gencode_Create(abcDNA, abcAA);
  bgAA   = p7_bg_Create(abcAA);
  bgDNA  = p7_bg_Create(abcDNA);
  p7_bg_SetLength(bgAA, L/3);
  p7_bg_SetLength(bgDNA, L/3);
  gm_fs5 = p7_profile_fs_Create(hmm->M, abcAA, 5);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_UNILOCAL);
  cx  = p7_chkgmx_Create(gm_fs5->M, L, L+1, 5, p7G_NSCELLS, p7_CHK_FORWARD);
  iv  = p7_ivx_Create(gm_fs5->M, p7P_5CODONS); 
  rng = malloc(sizeof(ESL_RANDOMNESS *) * N);
  tr  = malloc(sizeof(P7_TRACE *) * N);
  for (t = 0; t < N; t++) { rng[t] = esl_randomness_Create(1 + esl_rnd_Roll(r, 1000000)); tr[t] = p7_trace_fs_Create(); }
  esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);

  p7_chkgmx_GrowTo(cx, gm_fs5->M, L, L+1);
  p7_GForwardChk_Frameshift_5Codons(dsq, L, gm_fs5, cx, iv, &fsc);

  esl_stopwatch_Start(w);
  p7_GStochasticTraceChk_Frameshift(rng, N, dsq, L, gm_fs5, cx, iv, tr);
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");

  free(dsq);
  for (t = 0; t < N; t++) { esl_randomness_Destroy(rng[t]); p7_trace_fs_Destroy(tr[t]); }
  free(rng); free(tr);
  p7_chkgmx_Destroy(cx);
  p7_ivx_Destroy(iv);
  p7_profile_fs_Destroy(gm_fs5);
  p7_bg_Destroy(bgAA);
  p7_bg_Destroy(bgDNA);
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
#endif /*p7GENERIC_STOTRACE_FRAMESHIFT_BENCHMARK*/
/*---------------------- end, benchmark -------------------------*/


