/* Frameshift aware posterior decoding algorithms; SSE versions.
 *
 * Contents:
 *   1. Posterior decoding algorithms.
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 */

#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <x86intrin.h>

#include "easel.h"
#include "esl_sse.h"

#include "hmmer.h"
#include "impl_sse.h"

/*****************************************************************
 * 1. Posterior decoding algorithms.
 *****************************************************************/

/* Function:  p7_DomainDecoding_Frameshift_SSE()
 * Synopsis:  SSE posterior decoding of domain location for frameshift-aware alignments.
 *
 * Purpose:   Identical to <p7_DomainDecoding_Frameshift()> except that <om_fs>,
 *            <oxf>, <oxb> are SSE optimized versions. See
 *            <p7_DomainDecoding_Frameshift()> documentation for more info.
 *
 *            Fills <ddef->btot[i]>, <ddef->etot[i]>, <ddef->mocc[i]> for
 *            i = 0..L, then sets <ddef->L>.
 *
 *            <btot[i]> accumulates expected B-state entries: the contribution
 *            from position i-3 (a B-entry that emits a codon ending near i)
 *            is added to <btot[i-3]>, maintaining three independent running
 *            sums (one per residue-modulo-3 class).
 *
 *            <etot[i]> accumulates expected E-state exits similarly.
 *
 *            <mocc[i]> is the probability that residue i is in the core
 *            model (not in N/J/C flanking states), computed as 1 minus the
 *            sum of N, J, C occupancy probabilities.  Each of N, J, C
 *            contributes three terms because a 3-nucleotide codon spanning
 *            residues i-2..i, i-1..i+1, or i..i+2 all involve residue i.
 *
 *            Scale factors are handled using precomputed cumulative products
 *            sfwd[i] (product of all fwd scale factors at rows 0..i) and
 *            sbck[i] (product of all bck scale factors at rows i..L).  The
 *            posterior probability of any (fwd_j, bck_k) product is then
 *            fwd_stored(j) * bck_stored(k) * sfwd[j] * sbck[k] * inv_Z,
 *            where inv_Z = 1 / ((N(0)+N(1)+N(2)) * sbck[0]).  This handles
 *            both has_own_scales=FALSE (bck reuses fwd scale factors) and
 *            has_own_scales=TRUE (bck has its own scale factors).
 *
 * Args:      om_fs - optimized frameshift profile (for N/J/C transition odds ratios)
 *            oxf   - filled Forward DP matrix (parser mode, all xmx rows stored)
 *            oxb   - filled Backward DP matrix (parser mode, all xmx rows stored)
 *            ddef  - container for the results
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_DomainDecoding_Frameshift_SSE(const P7_FS_OPROFILE *om_fs, const P7_OMX *oxf, const P7_OMX *oxb,
                                  P7_DOMAINDEF *ddef)
{
  int    L   = oxf->L;
  float *log_sfwd = NULL;  /* log_sfwd[i] = sum of log(fwd_scale[0..i])        */
  float *log_sbck = NULL;  /* log_sbck[i] = sum of log(bck_scale[i..L])        */
  float  log_inv_Z;        /* -log(Z)                                          */
  float  njcp;
  int    i;
  int    status;

  ESL_ALLOC(log_sfwd, sizeof(float) * (L+2));
  ESL_ALLOC(log_sbck, sizeof(float) * (L+2));

  /* Cumulative log forward scale products: log_sfwd[i] = sum(log(fwd_scale_0..i)) */
  log_sfwd[0] = logf(oxf->xmx[0*p7X_NXCELLS + p7X_SCALE]);
  for (i = 1; i <= L; i++)
    log_sfwd[i] = log_sfwd[i-1] + logf(oxf->xmx[i*p7X_NXCELLS + p7X_SCALE]);

  /* Cumulative log backward scale products: log_sbck[i] = sum(log(bck_scale_i..L))
   * When has_own_scales=FALSE, bck reuses fwd scale factors; when TRUE, bck has
   * its own.  Either way, oxb->xmx[i*p7X_NXCELLS+p7X_SCALE] holds the right value. */
  log_sbck[L+1] = 0.0f;
  for (i = L; i >= 0; i--)
    log_sbck[i] = log_sbck[i+1] + logf(oxb->xmx[i*p7X_NXCELLS + p7X_SCALE]);

  /* Normalization in log space.
   * The three stored backward N values at rows 0, 1, 2 each carry different
   * cumulative scale factors (log_sbck[0], log_sbck[1], log_sbck[2]).  We
   * must restore each to its true magnitude before summing (LogSumExp). */
  log_inv_Z = -p7_FLogsum(
                 logf(oxb->xmx[0*p7X_NXCELLS + p7X_N]) + log_sbck[0],
                 p7_FLogsum(
                   logf(oxb->xmx[1*p7X_NXCELLS + p7X_N]) + log_sbck[1],
                   logf(oxb->xmx[2*p7X_NXCELLS + p7X_N]) + log_sbck[2]));

  /* Positions 0, 1, 2: no domains can have started or ended yet */
  ddef->btot[0] = 0.;
  ddef->btot[1] = 0.;
  ddef->btot[2] = 0.;
  ddef->etot[0] = 0.;
  ddef->etot[1] = 0.;
  ddef->etot[2] = 0.;
  ddef->mocc[0] = 0.;
  ddef->mocc[1] = 0.;
  ddef->mocc[2] = 0.;

  for (i = 3; i <= L; i++)
    {
      /* btot[i]: cumulative expected B-state entries up to position i.
       * B(i-3) contributes: fwd_B(i-3) * bck_B(i-3) (same-row product). */
      ddef->btot[i] = ddef->btot[i-3]
        + oxf->xmx[(i-3)*p7X_NXCELLS + p7X_B] * oxb->xmx[(i-3)*p7X_NXCELLS + p7X_B]
          * expf(log_sfwd[i-3] + log_sbck[i-3] + log_inv_Z);

      /* etot[i]: cumulative expected E-state exits up to position i.
       * E(i) contributes: fwd_E(i) * bck_E(i) (same-row product). */
      ddef->etot[i] = ddef->etot[i-3]
        + oxf->xmx[i*p7X_NXCELLS + p7X_E] * oxb->xmx[i*p7X_NXCELLS + p7X_E]
          * expf(log_sfwd[i] + log_sbck[i] + log_inv_Z);

      /* mocc[i]: probability residue i is in the core model = 1 - P(i in N/J/C).
       *
       * For each of N, J, C: three codon-length offsets contribute to residue i.
       * Codon ending at i:     fwd(i-3) * T_XL * bck(i),   rows (i-3, i)
       * Codon ending at i+1:   fwd(i-2) * T_XL * bck(i+1), rows (i-2, i+1)
       * Codon ending at i+2:   fwd(i-1) * T_XL * bck(i+2), rows (i-1, i+2)
       *
       * Cross-row posterior: fwd_stored(j) * bck_stored(k) * exp(log_sfwd[j] + log_sbck[k] + log_inv_Z)
       */
      njcp = 0.;

      /* N state */
      njcp += oxf->xmx[(i-3)*p7X_NXCELLS + p7X_N] * oxb->xmx[i*p7X_NXCELLS + p7X_N]
              * om_fs->xf[p7O_N][p7O_LOOP] * expf(log_sfwd[i-3] + log_sbck[i] + log_inv_Z);
      if (i < L)
        njcp += oxf->xmx[(i-2)*p7X_NXCELLS + p7X_N] * oxb->xmx[(i+1)*p7X_NXCELLS + p7X_N]
                * om_fs->xf[p7O_N][p7O_LOOP] * expf(log_sfwd[i-2] + log_sbck[i+1] + log_inv_Z);
      if (i < L-1)
        njcp += oxf->xmx[(i-1)*p7X_NXCELLS + p7X_N] * oxb->xmx[(i+2)*p7X_NXCELLS + p7X_N]
                * om_fs->xf[p7O_N][p7O_LOOP] * expf(log_sfwd[i-1] + log_sbck[i+2] + log_inv_Z);

      /* J state */
      njcp += oxf->xmx[(i-3)*p7X_NXCELLS + p7X_J] * oxb->xmx[i*p7X_NXCELLS + p7X_J]
              * om_fs->xf[p7O_J][p7O_LOOP] * expf(log_sfwd[i-3] + log_sbck[i] + log_inv_Z);
      if (i < L)
        njcp += oxf->xmx[(i-2)*p7X_NXCELLS + p7X_J] * oxb->xmx[(i+1)*p7X_NXCELLS + p7X_J]
                * om_fs->xf[p7O_J][p7O_LOOP] * expf(log_sfwd[i-2] + log_sbck[i+1] + log_inv_Z);
      if (i < L-1)
        njcp += oxf->xmx[(i-1)*p7X_NXCELLS + p7X_J] * oxb->xmx[(i+2)*p7X_NXCELLS + p7X_J]
                * om_fs->xf[p7O_J][p7O_LOOP] * expf(log_sfwd[i-1] + log_sbck[i+2] + log_inv_Z);

      /* C state */
      njcp += oxf->xmx[(i-3)*p7X_NXCELLS + p7X_C] * oxb->xmx[i*p7X_NXCELLS + p7X_C]
              * om_fs->xf[p7O_C][p7O_LOOP] * expf(log_sfwd[i-3] + log_sbck[i] + log_inv_Z);
      if (i < L)
        njcp += oxf->xmx[(i-2)*p7X_NXCELLS + p7X_C] * oxb->xmx[(i+1)*p7X_NXCELLS + p7X_C]
                * om_fs->xf[p7O_C][p7O_LOOP] * expf(log_sfwd[i-2] + log_sbck[i+1] + log_inv_Z);
      if (i < L-1)
        njcp += oxf->xmx[(i-1)*p7X_NXCELLS + p7X_C] * oxb->xmx[(i+2)*p7X_NXCELLS + p7X_C]
                * om_fs->xf[p7O_C][p7O_LOOP] * expf(log_sfwd[i-1] + log_sbck[i+2] + log_inv_Z);

      ddef->mocc[i] = 1. - njcp;
    }

  ddef->L = L;
  free(log_sfwd);
  free(log_sbck);
  return eslOK;

 ERROR:
  if (log_sfwd) free(log_sfwd);
  if (log_sbck) free(log_sbck);
  return status;
}

/*------------------ end, posterior decoding --------------------*/

/*****************************************************************
 * 2. Benchmark driver
 *****************************************************************/

/*------------------ end, benchmark driver ----------------------*/




/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7DOMDEF_FS_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"
static void
utest_domdef(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA, P7_CODONTABLE *codon_table, int M, int N)
{
  int  i,j;
  int  curr_L;
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_fs3 = p7_profile_fs_Create(M, abcAA, 3);
  P7_FS_OPROFILE *om_fs3 = p7_fs_oprofile_Create(M, abcAA, 3);
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq    = NULL;
  P7_TRACE       *tr     = p7_trace_Create();
  P7_OMX         *fwd    = p7_omx_Create(M, PARSER_ROWS_FWD, M);
  P7_OMX         *bwd    = p7_omx_Create(M, PARSER_ROWS_BWD, M);
  P7_GMX         *fgx    = p7_gmx_fs_Create(M, PARSER_ROWS_FWD, M, 0);
  P7_GMX         *bgx    = p7_gmx_fs_Create(M, PARSER_ROWS_BWD, M, 0);
  P7_IVX         *iv     = p7_ivx_Create(M, p7P_3CODONS);
  P7_DOMAINDEF   *oddef  = NULL; 
  P7_DOMAINDEF   *gddef  = NULL;
  float tolerance;
  int status;

    float fsc3, bsc3;
  float generic_fsc3, generic_bsc3;
  p7_FLogsumInit();
  if (p7_FLogsumError(-0.4, -0.5) > 0.0001) tolerance = 1.0;  /* weaker test against generic   */
  else tolerance = 0.001;   /* stronger test: FLogsum() is in slow exact mode. */  

  ESL_ALLOC(oddef, sizeof(P7_DOMAINDEF));
  ESL_ALLOC(gddef, sizeof(P7_DOMAINDEF));

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig(hmm, bgAA, gm, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs3, M, p7_LOCAL);
  p7_fs_oprofile_Convert(gm_fs3, om_fs3);
  p7_fs_oprofile_ReconfigLength(om_fs3, M);

    while (N--)
    {

      p7_ProfileEmit(r, hmm, gm, bgAA, sq, tr);
      curr_L = sq->n*3;

      if(dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) *(curr_L+2))) == NULL)  esl_fatal("malloc failed");

      j = 1;
      for(i = 1; i <= sq->n; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq+j);
        j+=3;
      }

      p7_fs_oprofile_ReconfigLength(om_fs3, sq->n);
      p7_fs_ReconfigLength(gm_fs3, sq->n);

      p7_omx_GrowTo(fwd, M, PARSER_ROWS_FWD, curr_L);
      p7_omx_GrowTo(bwd, M, PARSER_ROWS_BWD, curr_L);

      p7_ForwardParser_Frameshift_3Codons_SSE(dsq, gcode, curr_L, om_fs3, fwd, &fsc3);
      p7_BackwardParser_Frameshift_3Codons_SSE(dsq, gcode, curr_L, om_fs3, fwd, bwd, &generic_fsc3);

      oddef->mocc = oddef->btot = oddef->etot = NULL;
      ESL_ALLOC(oddef->mocc, sizeof(float) * (curr_L+1));
      ESL_ALLOC(oddef->btot, sizeof(float) * (curr_L+1));
      ESL_ALLOC(oddef->etot, sizeof(float) * (curr_L+1));
      
      p7_DomainDecoding_Frameshift_SSE(om_fs3, fwd, bwd, oddef);

      p7_gmx_fs_GrowTo(fgx, M, PARSER_ROWS_FWD, curr_L, 0);
      p7_gmx_fs_GrowTo(bgx, M, PARSER_ROWS_BWD, curr_L, 0);

      p7_ForwardParser_Frameshift_3Codons(dsq, gcode, curr_L, gm_fs3, fgx, iv, &bsc3);
      p7_BackwardParser_Frameshift_3Codons(dsq, gcode, curr_L, gm_fs3, bgx, iv, &generic_bsc3);

      gddef->mocc = gddef->btot = gddef->etot = NULL;
      ESL_ALLOC(gddef->mocc, sizeof(float) * (curr_L+1));
      ESL_ALLOC(gddef->btot, sizeof(float) * (curr_L+1));
      ESL_ALLOC(gddef->etot, sizeof(float) * (curr_L+1));
      p7_DomainDecoding_Frameshift(gm_fs3, fgx, bgx, gddef);
     
      for(i = 0; i <= curr_L; i++) {
        if(fabs(oddef->mocc[i] - gddef->mocc[i]) > tolerance) esl_fatal("domain def fs unit test failed at mocc");
        if(fabs(oddef->btot[i] - gddef->btot[i]) > tolerance) esl_fatal("domain def fs unit test failed at btot");
        if(fabs(oddef->etot[i] - gddef->etot[i]) > tolerance) esl_fatal("domain def fs unit test failed at etot");
      }

      free(oddef->mocc);
      free(oddef->btot);
      free(oddef->etot);
      free(gddef->mocc);
      free(gddef->btot);
      free(gddef->etot);

    }

  free(dsq);
  free(oddef);
  free(gddef);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bwd);
  p7_gmx_Destroy(fgx);
  p7_gmx_Destroy(bgx);
  p7_ivx_Destroy(iv);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_fs3);
  p7_fs_oprofile_Destroy(om_fs3);

ERROR:
  return;
}
#endif /*p7DOMDEF_FS_TESTDRIVE*/
/*--------------------- end, unit tests -------------------------*/




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7DOMDEF_FS_TESTDRIVE
/* 
  gcc -g -Wall -msse2 -std=gnu99 -o domdef_fs_utest -I.. -L.. -I../../easel -L../../easel -Dp7DOMDEF_FS_TESTDRIVE decoding_fs.c -lhmmer -leasel -lm
 ./domdef_fs_utest
*/
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",           0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                  0 },
  { "-M",        eslARG_INT,    "145", NULL, NULL,  NULL,  NULL, NULL, "size of random models to sample",                0 },
  { "-N",        eslARG_INT,    "100", NULL, NULL,  NULL,  NULL, NULL, "number of random sequences to sample",           0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE frameshift domain definition ";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_GENCODE    *gcode  = NULL;
  ESL_ALPHABET   *abcAA  = NULL;
  ESL_ALPHABET   *abcDNA = NULL;
  P7_BG          *bgAA   = NULL;
  P7_BG          *bgDNA  = NULL;
  P7_CODONTABLE  *ct     = NULL;
  int             M      = esl_opt_GetInteger(go, "-M");
  int             N      = esl_opt_GetInteger(go, "-N");

  if ((abcDNA = esl_alphabet_Create(eslDNA))      == NULL)  esl_fatal("failed to create alphabet");
  if ((bgDNA = p7_bg_Create(abcDNA))              == NULL)  esl_fatal("failed to create null model");
  if ((abcAA = esl_alphabet_Create(eslAMINO))     == NULL)  esl_fatal("failed to create alphabet");
  if ((bgAA = p7_bg_Create(abcAA))                == NULL)  esl_fatal("failed to create null model");
  if ((gcode  = esl_gencode_Create(abcDNA,abcAA)) == NULL)  esl_fatal("failed to create gencode");
  if ((ct     = p7_codontable_Create(gcode))      == NULL)  esl_fatal("failed to create codon table");

  utest_domdef(r, abcAA, abcDNA, gcode, bgAA, bgDNA, ct, M, N);

  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  p7_bg_Destroy(bgDNA);
  p7_bg_Destroy(bgAA);
  esl_gencode_Destroy(gcode);
  p7_codontable_Destroy(ct);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}

#endif /*p7DOMDEF_FS_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/


