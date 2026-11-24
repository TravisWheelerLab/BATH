/* Frameshift aware posterior decoding algorithms; NEON versions.
 *
 * Contents:
 *   1. Posterior decoding algorithms.
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 */

#include "p7_config.h"

#include <string.h>

#include <stdio.h>
#include <math.h>

#include <arm_neon.h>

#include "easel.h"
#include "esl_neon.h"

#include "hmmer.h"
#include "impl_neon.h"

/*****************************************************************
 * 1. Posterior decoding algorithms.
 *****************************************************************/

/* Function:  p7_DecodingRow_Frameshift()
 * Synopsis:  Posterior decoding of one frameshift row.
 *
 * Purpose:   Posterior probabilities of row <i> (1..L) into <ppv> and the
 *            special states <ppx> (E, N, J, B, C of the row), from Forward
 *            row <i> in <fv>, Backward row <i> in <bv>, the Backward
 *            special states <bxmx>, the Forward special states of row
 *            i-3 in <fx3> (read for i > 2), the cumulative log scale
 *            products <log_sfwd> and <log_sbck>, and <log_inv_Z>, as
 *            <p7_Decoding_Frameshift()> computes them. <ppv> may be
 *            <fv>. Same calculation as each row of
 *            <p7_Decoding_Frameshift()>, which calls it.
 *
 * Returns:   <eslOK> on success; <eslERANGE> if a scale factor or the
 *            special states' posteriors overflow.
 */
int
p7_DecodingRow_Frameshift(const P7_FS_OPROFILE *om_fs, const float *bxmx, const float *fx3, const float *log_sfwd, const float *log_sbck,
                         float log_inv_Z, const float32x4_t *fv, const float32x4_t *bv, float32x4_t *ppv, float *ppx, int i)
{
  int    Q = p7O_NQF(om_fs->M);
  int    q, c;
  float  factor_mdi, factor_njc;
  float32x4_t scv;
  float32x4_t zerov    = vdupq_n_f32(0.0f);
  float32x4_t  fI, fC[6], bM, bI;
  float  N_odds = om_fs->xf[p7O_N][p7O_LOOP];
  float  J_odds = om_fs->xf[p7O_J][p7O_LOOP];
  float  C_odds = om_fs->xf[p7O_C][p7O_LOOP];
  float  bck_N, bck_J, bck_C;
  float  N_pp, J_pp, C_pp;

  /* Backward special states at row i */
  bck_N = bxmx[i*p7X_NXCELLS + p7X_N];
  bck_J = bxmx[i*p7X_NXCELLS + p7X_J];
  bck_C = bxmx[i*p7X_NXCELLS + p7X_C];

  /* Scale factor for same-row M/I products */
  factor_mdi = expf(log_sfwd[i] + log_sbck[i] + log_inv_Z);

  if (isinf(factor_mdi)) return eslERANGE;

  /* M/I products f * b; scaled to posteriors below */
  for (q = 0; q < Q; q++)
    {
      bM = bv[q * p7X_NSCELLS + p7X_M];
      bI = bv[q * p7X_NSCELLS + p7X_I];

      fI = fv[q * p7X_NSCELLS_FS + p7X_FS_I];
      for (c = 0; c < 6; c++)
        fC[c] = fv[q * p7X_NSCELLS_FS + p7X_FS_M + c];

      ppv[q * p7X_NSCELLS_FS + p7X_FS_D] = zerov;
      ppv[q * p7X_NSCELLS_FS + p7X_FS_I] = vmulq_f32(fI, bI);
      for (c = 0; c < 6; c++)
        ppv[q * p7X_NSCELLS_FS + p7X_FS_M + c] = vmulq_f32(fC[c], bM);
    }

  /* N, J, C posteriors */
  if (i > 2) {
    factor_njc = expf(log_sfwd[i-3] + log_sbck[i] + log_inv_Z);
    N_pp = fx3[p7X_N] * bck_N * N_odds * factor_njc;
    J_pp = fx3[p7X_J] * bck_J * J_odds * factor_njc;
    C_pp = fx3[p7X_C] * bck_C * C_odds * factor_njc;
  } else {
    float factor_nsmall = expf(log_sbck[i] + log_inv_Z);
    N_pp = bck_N * factor_nsmall;
    J_pp = 0.0f;
    C_pp = 0.0f;
  }
  if (isinf(N_pp + J_pp + C_pp)) return eslERANGE;

  /* Scale the M/I products by the row factor to get posteriors */
  scv = vmovq_n_f32(factor_mdi);
  for (q = 0; q < Q; q++)
    for (c = 0; c < p7X_NSCELLS_FS; c++)
      ppv[q * p7X_NSCELLS_FS + c] = vmulq_f32(ppv[q * p7X_NSCELLS_FS + c], scv);

  /* Special states */
  ppx[p7X_E] = 0.0f;
  ppx[p7X_B] = 0.0f;
  ppx[p7X_N] = N_pp;
  ppx[p7X_J] = J_pp;
  ppx[p7X_C] = C_pp;
  return eslOK;
}


/* Function:  p7_Decoding_Frameshift()
 * Synopsis:  NEON posterior decoding of residue assignments; frameshift full-matrix version.
 *
 * Purpose:   NEON equivalent of <p7_GDecoding_Frameshift()>.  Given filled forward
 *            matrix <fwd> (8-cell FS layout, from <p7_Forward_Frameshift()>)
 *            and backward matrix <bck> (3-cell layout, from
 *            <p7_Backward_Frameshift()>), overwrites <fwd> with the posterior
 *            decoding matrix.
 *
 *            Each value at row i is the posterior probability that the
 *            state's emission ends at nucleotide i (M_C0..M_C5, I, N, J,
 *            C; D is zero). Emissions span 1-5 nucleotides, so the values
 *            at one row do not sum to one; over all rows they sum to the
 *            expected number of emissions.
 *
 *            Scale factors are handled via precomputed cumulative log products
 *            log_sfwd[i] and log_sbck[i].  N/J/C cross-row contributions
 *            (fwd[i-3] * bck[i]) are adjusted by factor_njc to bring them to the
 *            same scale as the same-row M/I products.
 *
 * Args:      om_fs - optimized frameshift profile (N/J/C transition odds ratios)
 *            fwd   - filled Forward matrix (OVERWRITTEN with PP on return)
 *            bck   - filled Backward matrix (read-only)
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEMEM> on allocation failure.
 */
int
p7_Decoding_Frameshift(const P7_FS_OPROFILE *om_fs, P7_OMX *fwd, const P7_OMX *bck)
{
  int    L = fwd->L;
  int    M = om_fs->M;
  int    Q = p7O_NQF(M);
  int    i, q;
  float *log_sfwd = NULL;
  float *log_sbck = NULL;
  float  log_inv_Z;
  float32x4_t zerov    = vdupq_n_f32(0.0f);
  float32x4_t *dpc;
  float  lagx[4][p7X_NXCELLS];          /* fwd special states of the last 4 rows, saved before overwrite */
  int    status;

  ESL_ALLOC(log_sfwd, sizeof(float) * (L + 2));
  ESL_ALLOC(log_sbck, sizeof(float) * (L + 2));

  /* Cumulative log scale products */
  log_sfwd[0] = logf(fwd->xmx[p7X_SCALE]);
  for (i = 1; i <= L; i++)
    log_sfwd[i] = log_sfwd[i-1] + logf(fwd->xmx[i*p7X_NXCELLS + p7X_SCALE]);

  log_sbck[L+1] = 0.0f;
  for (i = L; i >= 0; i--)
    log_sbck[i] = log_sbck[i+1] + logf(bck->xmx[i*p7X_NXCELLS + p7X_SCALE]);

  log_inv_Z = -p7_FLogsum(
                 logf(bck->xmx[p7X_N]) + log_sbck[0],
                 p7_FLogsum(
                   logf(bck->xmx[1*p7X_NXCELLS + p7X_N]) + log_sbck[1],
                   logf(bck->xmx[2*p7X_NXCELLS + p7X_N]) + log_sbck[2]));

  /* Save fwd row 0 special states before zeroing row 0 */
  memcpy(lagx[0], fwd->xmx, sizeof(float) * p7X_NXCELLS);

  /* Zero row 0 of the output PP matrix */
  dpc = fwd->dpf[0];
  for (q = 0; q < Q * p7X_NSCELLS_FS; q++) dpc[q] = zerov;
  fwd->xmx[p7X_E] = fwd->xmx[p7X_N] = fwd->xmx[p7X_J] = fwd->xmx[p7X_B] = fwd->xmx[p7X_C] = 0.0f;

  for (i = 1; i <= L; i++)
    {
      /* Save current fwd special states into lag buffer before overwriting */
      memcpy(lagx[i % 4], fwd->xmx + i*p7X_NXCELLS, sizeof(float) * p7X_NXCELLS);

      /* fwd special states 3 rows back: (i-3+4)%4 == (i+1)%4 */
      if ((status = p7_DecodingRow_Frameshift(om_fs, bck->xmx, lagx[(i + 1) % 4], log_sfwd, log_sbck, log_inv_Z,
                                              fwd->dpf[i], bck->dpf[i], fwd->dpf[i], fwd->xmx + i*p7X_NXCELLS, i)) != eslOK)
        { free(log_sfwd); free(log_sbck); return status; }
    }

  free(log_sfwd);
  free(log_sbck);
  return eslOK;

 ERROR:
  if (log_sfwd) free(log_sfwd);
  if (log_sbck) free(log_sbck);
  return status;
}


/* Function:  p7_DomainDecoding_Frameshift()
 * Synopsis:  NEON posterior decoding of domain location for frameshift-aware alignments.
 *
 * Purpose:   Identical to <p7_GDomainDecoding_Frameshift()> except that <om_fs>,
 *            <oxf>, <oxb> are NEON optimized versions. See
 *            <p7_GDomainDecoding_Frameshift()> documentation for more info.
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
p7_DomainDecoding_Frameshift(const P7_FS_OPROFILE *om_fs, const P7_OMX *oxf, const P7_OMX *oxb,
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
#ifdef p7DECODING_FS_BENCHMARK
/*
   gcc -g -O3 -march=armv8-a -std=gnu99 -o decoding_fs_benchmark -I.. -L.. -I../../easel -L../../easel -Dp7DECODING_FS_BENCHMARK decoding_fs.c -lhmmer -leasel -lm
   ./decoding_fs_benchmark <hmmfile>

   Because p7_Decoding_Frameshift() overwrites the forward matrix in-place with
   posterior probabilities, we cannot pre-fill fwd/bck once and loop over
   Decoding alone.  Instead we time a baseline loop (Forward+Backward only) and
   a benchmark loop (Forward+Backward+Decoding), then subtract to isolate the
   decoding cost.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"
#include "esl_stopwatch.h"

#include "hmmer.h"
#include "impl_neon.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-L",        eslARG_INT,   "1200", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs (nucleotides)",       0 },
  { "-N",        eslARG_INT,   "5000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for posterior residue decoding, NEON frameshift version";

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
  P7_FS_PROFILE  *gm_fs5  = NULL;
  P7_FS_OPROFILE *om_fs5  = NULL;
  P7_OMX         *fwd     = NULL;
  P7_OMX         *bck     = NULL;
  int             L       = esl_opt_GetInteger(go, "-L");
  int             N       = esl_opt_GetInteger(go, "-N");
  ESL_DSQ        *dsq     = malloc(sizeof(ESL_DSQ) * (L+2));
  int             i;
  float           fsc, bsc;
  double          base_time, bench_time, Mcs;

  if (p7_hmmfile_OpenE(hmmfile, NULL, &hfp, NULL) != eslOK) p7_Fail("Failed to open HMM file %s", hmmfile);
  if (p7_hmmfile_Read(hfp, &abcAA, &hmm)          != eslOK) p7_Fail("Failed to read HMM");

  abcDNA = esl_alphabet_Create(eslDNA);
  gcode  = esl_gencode_Create(abcDNA, abcAA);
  bgAA   = p7_bg_Create(abcAA);                 p7_bg_SetLength(bgAA, L/3);
  bgDNA  = p7_bg_Create(abcDNA);                p7_bg_SetLength(bgDNA, L);

  gm_fs5 = p7_profile_fs_Create(hmm->M, abcAA, p7P_5CODONS);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_UNILOCAL);

  om_fs5 = p7_fs_oprofile_Create(hmm->M, abcAA, p7P_5CODONS);
  p7_fs_oprofile_Convert(gm_fs5, om_fs5);
  p7_fs_oprofile_ReconfigLength(om_fs5, L/3);

  fwd = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS_FS);
  bck = p7_omx_Create_dpf(hmm->M, L, L, p7X_NSCELLS);
  P7_OIVX *ov5 = p7_oivx_Create(hmm->M, p7P_5CODONS);

  /* Baseline: time Forward+Backward alone */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
      p7_Forward_Frameshift (dsq, L, om_fs5, fwd,      ov5, &fsc);
      p7_Backward_Frameshift(dsq, L, om_fs5, fwd, bck, ov5, &bsc);
    }
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark: Forward+Backward+Decoding */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
      p7_Forward_Frameshift (dsq, L, om_fs5, fwd,      ov5, &fsc);
      p7_Backward_Frameshift(dsq, L, om_fs5, fwd, bck, ov5, &bsc);
      p7_Decoding_Frameshift(om_fs5, fwd, bck);  /* overwrites fwd in-place */
    }
  esl_stopwatch_Stop(w);
  bench_time = w->user - base_time;

  Mcs = (double) N * (double) L * (double) hmm->M * 1e-6 / bench_time;
  esl_stopwatch_Display(stdout, w, "# CPU time: ");
  printf("# M    = %d\n",   hmm->M);
  printf("# %.1f Mc/s\n", Mcs);

  free(dsq);
  p7_omx_Destroy(bck);
  p7_omx_Destroy(fwd);
  p7_oivx_Destroy(ov5);
  p7_fs_oprofile_Destroy(om_fs5);
  p7_profile_fs_Destroy(gm_fs5);
  p7_bg_Destroy(bgDNA);
  p7_bg_Destroy(bgAA);
  p7_hmm_Destroy(hmm);
  p7_hmmfile_Close(hfp);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  esl_stopwatch_Destroy(w);
  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7DECODING_FS_BENCHMARK*/
/*------------------ end, benchmark driver ----------------------*/


/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7DECODING_FS_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* utest_decoding_fs()
 * Compare p7_Decoding_Frameshift() against the generic
 * p7_GDecoding_Frameshift() on <N> random DNA sequences of length <L>,
 * as utest_decoding() does for the standard model. Each sequence is
 * decoded with the generic and the NEON full-matrix Forward/Backward,
 * and the two posterior matrices are compared within the relative
 * <tolerance>.
 */
static void
utest_decoding_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA, int M, int L, int N, float tolerance)
{
  char           *msg    = "decoding fs unit test failed";
  P7_HMM         *hmm    = NULL;
  P7_FS_PROFILE  *gm_fs5 = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE *om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
  ESL_DSQ        *dsq    = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_GMX         *gxf    = p7_gmx_Create(M, L, L, p7G_NSCELLS_FS);
  P7_GMX         *gxb    = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX         *gxpp   = p7_gmx_Create(M, L, L, p7G_NSCELLS_FS);
  P7_IVX         *iv5    = p7_ivx_Create(M, p7P_5CODONS);
  P7_OIVX        *ov5    = p7_oivx_Create(M, p7P_5CODONS);
  P7_OMX         *fwd    = p7_omx_Create_dpf(M, L, L, p7X_NSCELLS_FS);
  P7_OMX         *bck    = p7_omx_Create_dpf(M, L, L, p7X_NSCELLS);
  float           fsc, bsc;

  if (p7_hmm_Sample(r, M, abcAA, &hmm)                                  != eslOK) esl_fatal(msg);
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, L/3, p7_LOCAL)      != eslOK) esl_fatal(msg);
  if (p7_fs_oprofile_Convert(gm_fs5, om_fs5)                            != eslOK) esl_fatal(msg);
  if (p7_fs_oprofile_ReconfigLength(om_fs5, L/3)                        != eslOK) esl_fatal(msg);

  while (N--)
    {
      if (esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq)                 != eslOK) esl_fatal(msg);

      /* generic: gxf becomes the posterior matrix */
      if (p7_GForward_Frameshift (dsq, L, gm_fs5, gxf, iv5, &fsc)       != eslOK) esl_fatal(msg);
      if (p7_GBackward_Frameshift(dsq, L, gm_fs5, gxb, iv5, &bsc)       != eslOK) esl_fatal(msg);
      if (p7_GDecoding_Frameshift(gm_fs5, gxf, gxb)                     != eslOK) esl_fatal(msg);

      /* NEON: fwd becomes the posterior matrix */
      if (p7_Forward_Frameshift (dsq, L, om_fs5, fwd, ov5, &fsc)        != eslOK) esl_fatal(msg);
      if (p7_Backward_Frameshift(dsq, L, om_fs5, fwd, bck, ov5, &bsc)   != eslOK) esl_fatal(msg);
      if (p7_Decoding_Frameshift(om_fs5, fwd, bck)                      != eslOK) esl_fatal(msg);
      if (p7_omx_FDeconvert(fwd, gxpp)                                  != eslOK) esl_fatal(msg);

      if (p7_gmx_Compare(gxpp, gxf, tolerance)                          != eslOK) esl_fatal(msg);
    }

  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_gmx_Destroy(gxf);
  p7_gmx_Destroy(gxb);
  p7_gmx_Destroy(gxpp);
  p7_ivx_Destroy(iv5);
  p7_oivx_Destroy(ov5);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bck);
  p7_profile_fs_Destroy(gm_fs5);
  p7_fs_oprofile_Destroy(om_fs5);
}

static void
utest_domdef(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA, int M, int L, int N)
{
  char           *msg    = "domain def fs unit test failed";
  P7_HMM         *hmm    = NULL;
  P7_FS_PROFILE  *gm_fs3 = p7_profile_fs_Create(M, abcAA, 3);
  P7_FS_OPROFILE *om_fs3 = p7_fs_oprofile_Create(M, abcAA, 3);
  ESL_DSQ        *dsq    = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX         *fwd    = p7_omx_Create(M, PARSER_ROWS_FWD, L);
  P7_OMX         *bwd    = p7_omx_Create(M, PARSER_ROWS_BWD, L);
  P7_GMX         *fgx    = p7_gmx_Create(M, PARSER_ROWS_FWD, L, p7G_NSCELLS);
  P7_GMX         *bgx    = p7_gmx_Create(M, PARSER_ROWS_BWD, L, p7G_NSCELLS);
  P7_IVX         *iv     = p7_ivx_Create(M, p7P_3CODONS);
  P7_OIVX        *ov3    = p7_oivx_Create(M, p7P_3CODONS);
  P7_DOMAINDEF    oddef, gddef;
  float           tolerance;
  int             i;

  if (p7_FLogsumError(-0.4, -0.5) > 0.0001) tolerance = 1.0;    /* weaker test against generic   */
  else                                      tolerance = 0.001;  /* stronger test: FLogsum() is in slow exact mode. */

  oddef.mocc = malloc(sizeof(float) * (L+1));
  oddef.btot = malloc(sizeof(float) * (L+1));
  oddef.etot = malloc(sizeof(float) * (L+1));
  gddef.mocc = malloc(sizeof(float) * (L+1));
  gddef.btot = malloc(sizeof(float) * (L+1));
  gddef.etot = malloc(sizeof(float) * (L+1));

  if (p7_hmm_Sample(r, M, abcAA, &hmm)                             != eslOK) esl_fatal(msg);
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs3, L/3, p7_LOCAL) != eslOK) esl_fatal(msg);
  if (p7_fs_oprofile_Convert(gm_fs3, om_fs3)                       != eslOK) esl_fatal(msg);
  if (p7_fs_oprofile_ReconfigLength(om_fs3, L/3)                   != eslOK) esl_fatal(msg);

  while (N--)
    {
      if (esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq)                      != eslOK) esl_fatal(msg);

      if (p7_ForwardParser_Frameshift_3Codons (dsq, L, om_fs3, fwd,      ov3, NULL) != eslOK) esl_fatal(msg);
      if (p7_BackwardParser_Frameshift_3Codons(dsq, L, om_fs3, fwd, bwd, ov3, NULL) != eslOK) esl_fatal(msg);
      if (p7_DomainDecoding_Frameshift(om_fs3, fwd, bwd, &oddef)             != eslOK) esl_fatal(msg);

      if (p7_GForwardParser_Frameshift_3Codons (dsq, L, gm_fs3, fgx, iv, NULL) != eslOK) esl_fatal(msg);
      if (p7_GBackwardParser_Frameshift_3Codons(dsq, L, gm_fs3, bgx, iv, NULL) != eslOK) esl_fatal(msg);
      if (p7_GDomainDecoding_Frameshift(gm_fs3, fgx, bgx, &gddef)            != eslOK) esl_fatal(msg);

      for (i = 0; i <= L; i++)
        {
          if (fabs(oddef.mocc[i] - gddef.mocc[i]) > tolerance) esl_fatal("domain def fs unit test failed at mocc");
          if (fabs(oddef.btot[i] - gddef.btot[i]) > tolerance) esl_fatal("domain def fs unit test failed at btot");
          if (fabs(oddef.etot[i] - gddef.etot[i]) > tolerance) esl_fatal("domain def fs unit test failed at etot");
        }
    }

  free(oddef.mocc);  free(oddef.btot);  free(oddef.etot);
  free(gddef.mocc);  free(gddef.btot);  free(gddef.etot);
  free(dsq);
  p7_hmm_Destroy(hmm);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bwd);
  p7_gmx_Destroy(fgx);
  p7_gmx_Destroy(bgx);
  p7_ivx_Destroy(iv);
  p7_oivx_Destroy(ov3);
  p7_profile_fs_Destroy(gm_fs3);
  p7_fs_oprofile_Destroy(om_fs3);
}
#endif /*p7DECODING_FS_TESTDRIVE*/

/*--------------------- end, unit tests -------------------------*/


/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7DECODING_FS_TESTDRIVE
/*
  gcc -g -Wall -march=armv8-a -std=gnu99 -o decoding_fs_utest -I.. -L.. -I../../easel -L../../easel -Dp7DECODING_FS_TESTDRIVE decoding_fs.c -lhmmer -leasel -lm
  ./decoding_fs_utest
*/
#include "p7_config.h"

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_random.h"

#include "hmmer.h"
#include "impl_neon.h"

static ESL_OPTIONS options[] = {
  /* name  type         default  env   range togs  reqs  incomp  help                                         docgrp */
  { "-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",                            0 },
  { "-s",  eslARG_INT,     "42",  NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",                  0 },
  { "-t",  eslARG_REAL,  "0.01",  NULL, NULL, NULL, NULL, NULL, "relative tolerance for the decoding test",       0 },
  { "-L",  eslARG_INT,    "120",  NULL, NULL, NULL, NULL, NULL, "length of random DNA seqs (nucleotides)",        0 },
  { "-M",  eslARG_INT,     "40",  NULL, NULL, NULL, NULL, NULL, "length of sampled test profile",                 0 },
  { "-N",  eslARG_INT,     "10",  NULL, NULL, NULL, NULL, NULL, "number of random test sequences",                0 },
  { 0,0,0,0,0,0,0,0,0,0},
};

static char usage[]  = "[-options]";
static char banner[] = "test driver for NEON frameshift domain definition";

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
  int             M      = esl_opt_GetInteger(go, "-M");
  int             N      = esl_opt_GetInteger(go, "-N");
  int             L      = esl_opt_GetInteger(go, "-L");
  float           tol    = esl_opt_GetReal   (go, "-t");

  if ((abcDNA = esl_alphabet_Create(eslDNA))      == NULL)  esl_fatal("failed to create alphabet");
  if ((bgDNA = p7_bg_Create(abcDNA))              == NULL)  esl_fatal("failed to create null model");
  if ((abcAA = esl_alphabet_Create(eslAMINO))     == NULL)  esl_fatal("failed to create alphabet");
  if ((bgAA = p7_bg_Create(abcAA))                == NULL)  esl_fatal("failed to create null model");
  if ((gcode  = esl_gencode_Create(abcDNA,abcAA)) == NULL)  esl_fatal("failed to create gencode");

  p7_FLogsumInit();

  utest_domdef      (r, abcAA, abcDNA, gcode, bgAA, bgDNA, M, L, N);
  utest_decoding_fs (r, abcAA, abcDNA, gcode, bgAA, bgDNA, M, L, N, tol);

  esl_alphabet_Destroy(abcDNA);
  esl_alphabet_Destroy(abcAA);
  p7_bg_Destroy(bgDNA);
  p7_bg_Destroy(bgAA);
  esl_gencode_Destroy(gcode);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  return eslOK;
}

#endif /*p7DECODING_FS_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/
