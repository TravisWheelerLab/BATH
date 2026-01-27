/* Stochastic traceback of a Forward matrix; frameshift aware version.
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

/* Function:  p7_StochasticTrace_Frameshift()
 * Synopsis:  Stochastic traceback of a frameshift-aware Forward matrix.
 *
 * Purpose:   Stochastic traceback of Forward matrix <gx> to
 *            sample an alignment of digital sequence <dsq>
 *            (of length <L>) to the profile <gm_fs>. 
 *            
 *            The sampled traceback is returned in <tr>, which the
 *            caller must have at least made an initial allocation of
 *            (the <tr> will be grown as needed here).
 *
 * Args:      r      - source of random numbers
 *            dsq    - digital sequence aligned to, 1..L 
 *            L      - length of dsq
 *            gm_fs    - profile
 *            mx     - Forward matrix to trace, L x M
 *            tr     - storage for the recovered traceback.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_StochasticTrace_Frameshift(ESL_RANDOMNESS *r, const ESL_DSQ *dsq, int L, const P7_FS_PROFILE *gm_fs, const P7_GMX *gx, P7_TRACE *tr)
{
  int     status;
  int     i;      /* position in seq (1..L) */
  int     k;      /* position in model (1..M) */ 
  int     c;
  int     M   = gm_fs->M;
  float **dp  = gx->dp;
  float  *xmx = gx->xmx;
  float const *tsc  = gm_fs->tsc;
  float  *sc;      /* scores of possible choices: up to 2M-1, in the case of exits to E  */
  int     scur, sprv;

  /* we'll index M states as 1..M, and D states as 2..M = M+2..2M: M0, D1 are impossibles. */
  ESL_ALLOC(sc, sizeof(float) * (2*M+1)); 

  k = 0;
  i = L;    
  c = 0;  
  if ((status = p7_trace_fs_Append(tr, p7T_T, k, i, c)) != eslOK) goto ERROR;
  if ((status = p7_trace_fs_Append(tr, p7T_C, k, i, c)) != eslOK) goto ERROR;
  sprv = p7T_C;

  while (sprv != p7T_S) 
  { 
    switch (sprv) {
      /* check all three frames of C as well as E(i) */
      case p7T_C:
        if   (XMX_FS(i,p7G_C) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible C reached at i=%d", i);
        if   (i < 4) { scur = p7T_E; break; }
        
        sc[0] = XMX_FS(i-3, p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP];
        sc[1] = XMX_FS(i-2, p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP];
        sc[2] = XMX_FS(i-1, p7G_C) + gm_fs->xsc[p7P_C][p7P_LOOP];
        sc[3] = XMX_FS(i,   p7G_E) + gm_fs->xsc[p7P_E][p7P_MOVE];
        
        esl_vec_FLogNorm(sc, 4);
        scur = (esl_rnd_FChoose(r, sc, 4) < 3) ? p7T_C : p7T_E;
        break;
  
      /* E connects from any M or D state. k set here */
      case p7T_E:  
        if (XMX_FS(i, p7G_E) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible E reached at i=%d", i);
      if (p7_fs_profile_IsLocal(gm_fs)) 
      { /* local models come from any M, D */
        sc[0] = sc[M+1] = -eslINFINITY;
        for (k = 1; k <= M; k++) sc[k]   = MMX_FS(i,k,p7G_C0);
        for (k = 2; k <= M; k++) sc[k+M] = DMX_FS(i,k);
        esl_vec_FLogNorm(sc, 2*M+1); /* now sc is a prob vector */
        k = esl_rnd_FChoose(r, sc, 2*M+1);
        if (k <= M)    scur = p7T_M;
        else { k -= M; scur = p7T_D; }
      } 
      else 
      {     /* glocal models come from M_M or D_M  */
        k     = M;
        sc[0] = MMX_FS(i,M,p7G_C0);
        sc[1] = DMX_FS(i,M);
        esl_vec_FLogNorm(sc, 2); /* now sc is a prob vector */
        scur = (esl_rnd_FChoose(r, sc, 2) == 0) ? p7T_M : p7T_D;
      }
      break;

      /* M connects from {MDI} i-1,k-1, or B */
      case p7T_M:
       if (MMX_FS(i,k,p7G_C0) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible M reached at k=%d,i=%d", k,i);
       
       sc[0] = XMX_FS(i,p7G_B)          + TSC(p7P_BM, k-1);
       sc[1] = MMX_FS(i,k-1,p7G_C0)  + TSC(p7P_MM, k-1);
       sc[2] = IMX_FS(i,k-1)         + TSC(p7P_IM, k-1);
       sc[3] = DMX_FS(i,k-1)         + TSC(p7P_DM, k-1);
       esl_vec_FLogNorm(sc, 4);
       switch (esl_rnd_FChoose(r, sc, 4)) {
          case 0: scur = p7T_B;  break;
          case 1: scur = p7T_M;  break;
          case 2: scur = p7T_I;  break;
          case 3: scur = p7T_D;  break;
          default: ESL_XEXCEPTION(eslFAIL, "bogus state in traceback");
         }
       k--; 
       break;

      /* D connects from M,D at i,k-1 */
      case p7T_D:
        if (DMX_FS(i, k) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible D reached at k=%d,i=%d", k,i);
        sc[0] = MMX_FS(i, k-1,p7G_C0) + TSC(p7P_MD, k-1);
        sc[1] = DMX_FS(i, k-1) + TSC(p7P_DD, k-1);
        esl_vec_FLogNorm(sc, 2); 
        scur = (esl_rnd_FChoose(r, sc, 2) == 0) ? p7T_M : p7T_D;
        k--;
        break;

      /* I connects from M,I at i-1,k */
      case p7T_I:
        if (IMX_FS(i,k) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible I reached at k=%d,i=%d", k,i);
        sc[0] = MMX_FS(i-3,k, p7G_C0) + TSC(p7P_MI, k);
        sc[1] = IMX_FS(i-3,k) + TSC(p7P_II, k);

        esl_vec_FLogNorm(sc, 2); 
        scur = (esl_rnd_FChoose(r, sc, 2) == 0) ? p7T_M : p7T_I;
        i-=3;
        break;

      /* N connects from S, N */
      case p7T_N:
        if (XMX_FS(i, p7G_N) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible N reached at i=%d", i);
        scur = (i == 0) ? p7T_S : p7T_N;
        break;

      /* B connects from N, J */
      case p7T_B:    
        if (XMX_FS(i,p7G_B) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible B reached at i=%d", i);

        sc[0] = XMX_FS(i, p7G_N) + gm_fs->xsc[p7P_N][p7P_MOVE];
        sc[1] = XMX_FS(i, p7G_J) + gm_fs->xsc[p7P_J][p7P_MOVE];
        esl_vec_FLogNorm(sc, 2); 
        scur = (esl_rnd_FChoose(r, sc, 2) == 0) ? p7T_N : p7T_J;
        break;

      /* J connects from E(i) or J(i-1) */
      case p7T_J:  
        if (XMX_FS(i,p7G_J) == -eslINFINITY) ESL_XEXCEPTION(eslFAIL, "impossible J reached at i=%d", i);
        
        if   (i < 4) { scur = p7T_E; break; }
        sc[0] = XMX_FS(i-3,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP];
        sc[1] = XMX_FS(i-2,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP];
        sc[2] = XMX_FS(i-1,p7G_J) + gm_fs->xsc[p7P_J][p7P_LOOP]; 
        sc[3] = XMX_FS(i,  p7G_E) + gm_fs->xsc[p7P_E][p7P_LOOP];
        esl_vec_FLogNorm(sc, 4); 
        scur = (esl_rnd_FChoose(r, sc, 4) < 3) ? p7T_J : p7T_E;  
        
        break;
      default: ESL_XEXCEPTION(eslFAIL, "bogus state in traceback");
      } /* end switch over statetype[tpos-1] */

      if(scur == p7T_M)
      {
       sc[0] = MMX_FS(i,k,p7G_C1);
       sc[1] = MMX_FS(i,k,p7G_C2);
       sc[2] = MMX_FS(i,k,p7G_C3);
       sc[3] = MMX_FS(i,k,p7G_C4);
       sc[4] = MMX_FS(i,k,p7G_C5);
       esl_vec_FLogNorm(sc, 5);
       c = esl_rnd_FChoose(r, sc, 5) + 1;
       if(i - c < 1) scur = p7T_B;
      }
      else c = 0; 
      /* Append this state and the current i,k to be explained to the growing trace */
      if ((status = p7_trace_fs_Append(tr, scur, k, i, c)) != eslOK) goto ERROR;

      /* For NCJ, we had to defer i decrement. */
      if ( (scur == p7T_N || scur == p7T_C || scur == p7T_J) && scur == sprv) i--;

      sprv = scur;
      i-=c;
    } /* end traceback, at S state */

  if ((status = p7_trace_fs_Reverse(tr)) != eslOK) goto ERROR;
  tr->M = gm_fs->M;
  tr->L = L;
  free(sc);
  return eslOK;

 ERROR:
  if (sc != NULL) free(sc);
  return status;
}



/*------------------- end, stochastic trace ---------------------*/





/*****************************************************************
 * 2. Benchmark driver
 *****************************************************************/
#ifdef p7STOTRACE_FRAMESHIFT_BENCHMARK
/*
   gcc -g -O2      -o stotrace_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7STOTRACE_FRAMESHIFT_BENCHMARK stotrace_frameshift.c -lhmmer -leasel -lm
   icc -O3 -static -o stotrace_frameshift_benchmark -I. -L. -I../easel -L../easel -Dp7STOTRACE_FRAMESHIFT_BENCHMARK stotrace_frameshift.c -lhmmer -leasel -lm
   ./stotrace_frameshift_benchmark <hmmfile>
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
  P7_FS_PROFILE  *gm_fs   = NULL;
  P7_GMX         *fwd     = NULL;
  P7_IVX         *iv      = NULL;
  P7_TRACE       *tr      = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
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
  gm_fs = p7_profile_fs_Create(hmm->M, abcAA);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs, L/3, p7_UNILOCAL);
  fwd = p7_gmx_fs_Create(gm_fs->M, L, L, p7P_5CODONS);
  iv  = p7_ivx_Create(gm_fs->M, p7P_5CODONS); 
  tr  = p7_trace_fs_Create();
  esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);

  p7_Forward_Frameshift(dsq, gcode, L, gm_fs, fwd, iv, &fsc);

  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      p7_StochasticTrace_Frameshift(r, dsq, L, gm_fs, fwd, tr);
      p7_trace_Reuse(tr);
    }
  esl_stopwatch_Stop(w);
  esl_stopwatch_Display(stdout, w, "# CPU time: ");

  free(dsq);
  p7_trace_fs_Destroy(tr);
  p7_gmx_Destroy(fwd);
  p7_ivx_Destroy(iv);
  p7_profile_fs_Destroy(gm_fs);
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
#endif /*p7STOTRACE_FRAMESHIFT_BENCHMARK*/
/*---------------------- end, benchmark -------------------------*/


