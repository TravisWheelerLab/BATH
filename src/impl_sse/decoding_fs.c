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

/* Function:  p7_Decoding_Frameshift()
 * Synopsis:  SSE posterior decoding of residue assignments; frameshift full-matrix version.
 *
 * Purpose:   SSE equivalent of <p7_Decoding_Frameshift()>.  Given filled forward
 *            matrix <fwd> (8-cell FS layout, from <p7_Forward_Frameshift()>)
 *            and backward matrix <bck> (3-cell layout, from
 *            <p7_Backward_Frameshift()>), overwrites <fwd> with the posterior
 *            decoding matrix.
 *
 *            At each row i, all posterior probabilities (M_C0..M_C5, I, N, J, C)
 *            are normalized to sum to 1 using a per-row denominator, exactly as in
 *            the generic <p7_Decoding_Frameshift()>.  Only M_C0 and I contribute
 *            to the denominator (D states are zero; M_C1..C5 are normalized but do
 *            not add to the denominator).
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
  int    i, q, c;
  float *log_sfwd = NULL;
  float *log_sbck = NULL;
  float  log_inv_Z;
  float  factor_mdi, factor_njc;
  __m128 scv, dv;
  __m128 zerov    = _mm_setzero_ps();
  float  raw_denom;
  __m128 *dpc, *dbv;
  __m128  fI, fC[6], bM, bI;
  /* Circular lag buffer depth 4: nlag[i%4] = fwd_N(i) saved before overwrite */
  float  nlag[4], jlag[4], clag[4];
  float  fN3, fJ3, fC3;
  float  bck_N, bck_J, bck_C;
  float  N_odds = om_fs->xf[p7O_N][p7O_LOOP];
  float  J_odds = om_fs->xf[p7O_J][p7O_LOOP];
  float  C_odds = om_fs->xf[p7O_C][p7O_LOOP];
  float  N_pp, J_pp, C_pp, inv_denom;
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
  nlag[0] = fwd->xmx[p7X_N];
  jlag[0] = fwd->xmx[p7X_J];
  clag[0] = fwd->xmx[p7X_C];
  nlag[1] = nlag[2] = nlag[3] = 0.0f;
  jlag[1] = jlag[2] = jlag[3] = 0.0f;
  clag[1] = clag[2] = clag[3] = 0.0f;

  /* Zero row 0 of the output PP matrix */
  dpc = fwd->dpf[0];
  for (q = 0; q < Q * p7X_NSCELLS_FS; q++) dpc[q] = zerov;
  fwd->xmx[p7X_E] = fwd->xmx[p7X_N] = fwd->xmx[p7X_J] = fwd->xmx[p7X_B] = fwd->xmx[p7X_C] = 0.0f;

  for (i = 1; i <= L; i++)
    {
      /* Save current fwd special states into lag buffer before overwriting */
      nlag[i % 4] = fwd->xmx[i*p7X_NXCELLS + p7X_N];
      jlag[i % 4] = fwd->xmx[i*p7X_NXCELLS + p7X_J];
      clag[i % 4] = fwd->xmx[i*p7X_NXCELLS + p7X_C];

      /* fwd special states 3 rows back: (i-3+4)%4 == (i+1)%4 */
      fN3 = nlag[(i + 1) % 4];
      fJ3 = jlag[(i + 1) % 4];
      fC3 = clag[(i + 1) % 4];

      /* Backward special states at row i */
      bck_N = bck->xmx[i*p7X_NXCELLS + p7X_N];
      bck_J = bck->xmx[i*p7X_NXCELLS + p7X_J];
      bck_C = bck->xmx[i*p7X_NXCELLS + p7X_C];

      /* Scale factor for same-row M/I products */
      factor_mdi = expf(log_sfwd[i] + log_sbck[i] + log_inv_Z);

      if (isinf(factor_mdi)) { free(log_sfwd); free(log_sbck); return eslERANGE; }

      /* ------- Pass 1: compute raw M/I products, accumulate denom ------- */
      dv  = zerov;
      dpc = fwd->dpf[i];
      dbv = bck->dpf[i];

      for (q = 0; q < Q; q++)
        {
          bM = dbv[q * p7X_NSCELLS + p7X_M];
          bI = dbv[q * p7X_NSCELLS + p7X_I];

          fI = dpc[q * p7X_NSCELLS_FS + p7X_FS_I];
          for (c = 0; c < 6; c++)
            fC[c] = dpc[q * p7X_NSCELLS_FS + p7X_FS_M + c];

          dpc[q * p7X_NSCELLS_FS + p7X_FS_D] = zerov;
          dpc[q * p7X_NSCELLS_FS + p7X_FS_I] = _mm_mul_ps(fI, bI);
          for (c = 0; c < 6; c++)
            dpc[q * p7X_NSCELLS_FS + p7X_FS_M + c] = _mm_mul_ps(fC[c], bM);

          /* Only C0 and I contribute to the per-row normalization denominator */
          dv = _mm_add_ps(dv, _mm_add_ps(dpc[q * p7X_NSCELLS_FS + p7X_FS_M],
                                          dpc[q * p7X_NSCELLS_FS + p7X_FS_I]));
        }

      /* Horizontal sum of dv to get scalar raw MDI denom */
      dv = _mm_add_ps(dv, _mm_movehl_ps(dv, dv));
      dv = _mm_add_ss(dv, _mm_shuffle_ps(dv, dv, _MM_SHUFFLE(0,0,0,1)));
      raw_denom = _mm_cvtss_f32(dv);

      /* N, J, C contributions to denominator */
      if (i > 2) {
        factor_njc = expf(log_sfwd[i-3] + log_sbck[i] + log_inv_Z);
        N_pp = fN3 * bck_N * N_odds * factor_njc;
        J_pp = fJ3 * bck_J * J_odds * factor_njc;
        C_pp = fC3 * bck_C * C_odds * factor_njc;
      } else {
        float factor_nsmall = expf(log_sbck[i] + log_inv_Z);
        N_pp = bck_N * factor_nsmall;
        J_pp = 0.0f;
        C_pp = 0.0f;
      }
      inv_denom = 1.0f / (raw_denom * factor_mdi + N_pp + J_pp + C_pp); 

      if(isinf(inv_denom)) { free(log_sfwd); free(log_sbck); return eslERANGE; }

      /* ------- Pass 2: normalize all M/I vectors ------- */
      scv = _mm_set1_ps(factor_mdi * inv_denom);
      for (q = 0; q < Q; q++)
        for (c = 0; c < p7X_NSCELLS_FS; c++)
          fwd->dpf[i][q * p7X_NSCELLS_FS + c] = _mm_mul_ps(fwd->dpf[i][q * p7X_NSCELLS_FS + c], scv);

      /* Normalized special states */
      fwd->xmx[i*p7X_NXCELLS + p7X_E] = 0.0f;
      fwd->xmx[i*p7X_NXCELLS + p7X_B] = 0.0f;
      fwd->xmx[i*p7X_NXCELLS + p7X_N] = N_pp * inv_denom;
      fwd->xmx[i*p7X_NXCELLS + p7X_J] = J_pp * inv_denom;
      fwd->xmx[i*p7X_NXCELLS + p7X_C] = C_pp * inv_denom;
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
   gcc -g -O2      -o decoding_fs_benchmark -I.. -L.. -I../../easel -L../../easel -Dp7DECODING_FS_BENCHMARK decoding_fs.c -lhmmer -leasel -lm
   icc -O3 -static -o decoding_fs_benchmark -I.. -L.. -I../../easel -L../../easel -Dp7DECODING_FS_BENCHMARK decoding_fs.c -lhmmer -leasel -lm
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
#include "impl_sse.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,  NULL,  NULL, NULL, "show brief help on version and usage",             0 },
  { "-s",        eslARG_INT,     "42", NULL, NULL,  NULL,  NULL, NULL, "set random number seed to <n>",                    0 },
  { "-L",        eslARG_INT,   "1200", NULL, "n>0", NULL,  NULL, NULL, "length of random target seqs (nucleotides)",       0 },
  { "-N",        eslARG_INT,   "5000", NULL, "n>0", NULL,  NULL, NULL, "number of random target seqs",                     0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "benchmark driver for posterior residue decoding, SSE frameshift version";

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

  /* Baseline: time Forward+Backward alone */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
      p7_Forward_Frameshift (dsq, L, om_fs5, fwd, &fsc);
      p7_Backward_Frameshift(dsq, L, om_fs5, fwd, bck, &bsc);
    }
  esl_stopwatch_Stop(w);
  base_time = w->user;

  /* Benchmark: Forward+Backward+Decoding */
  esl_stopwatch_Start(w);
  for (i = 0; i < N; i++)
    {
      esl_rsq_xfIID(r, bgDNA->f, abcDNA->K, L, dsq);
      p7_Forward_Frameshift (dsq, L, om_fs5, fwd, &fsc);
      p7_Backward_Frameshift(dsq, L, om_fs5, fwd, bck, &bsc);
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
 * Compare p7_Decoding_Frameshift() against the generic p7_Decoding_Frameshift().
 * Both are given the same sequence and profile.  We run the generic full-matrix
 * forward/backward to produce a generic PP matrix (gx1), and the SSE full-matrix
 * forward/backward to produce an SSE PP matrix (fwd).  We compare the special-state
 * posteriors (N, J, C per row) between the two.
 */
static void
utest_decoding_fs(ESL_RANDOMNESS *r, ESL_ALPHABET *abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode, P7_BG *bgAA, P7_BG *bgDNA, P7_CODONTABLE *codon_table, int M, int N)
{
  int  i, j;
  int  curr_L;
  char           *msg    = "decoding fs unit test failed";
  P7_HMM         *hmm    = NULL;
  P7_PROFILE     *gm     = p7_profile_Create(M, abcAA);
  P7_FS_PROFILE  *gm_fs5 = p7_profile_fs_Create(M, abcAA, 5);
  P7_FS_OPROFILE *om_fs5 = p7_fs_oprofile_Create(M, abcAA, 5);
  ESL_SQ         *sq     = esl_sq_CreateDigital(abcAA);
  ESL_DSQ        *dsq    = NULL;
  P7_TRACE       *tr     = p7_trace_Create();
  P7_GMX         *gx1    = p7_gmx_Create(M, M, M, p7G_NSCELLS_FS);
  P7_GMX         *gx2    = p7_gmx_Create(M, M, M, p7G_NSCELLS);
  P7_IVX         *iv5    = p7_ivx_Create(M, p7P_5CODONS);
  P7_OMX         *fwd    = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS_FS);
  P7_OMX         *bck    = p7_omx_Create_dpf(M, M, M, p7X_NSCELLS);
  P7_GMX         *gxpp   = p7_gmx_Create(M, M, M, p7G_NSCELLS_FS);
  float           fsc, bsc;
  float           tolerance;

  tolerance = (p7_FLogsumError(-0.4, -0.5) > 0.0001) ? 0.2f : 0.001f;

  p7_hmm_Sample(r, M, abcAA, &hmm);
  p7_ProfileConfig(hmm, bgAA, gm, M, p7_LOCAL);
  p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs5, M, p7_LOCAL);
  p7_fs_oprofile_Convert(gm_fs5, om_fs5);
  p7_fs_oprofile_ReconfigLength(om_fs5, M);

  while (N--)
    {
      p7_ProfileEmit(r, hmm, gm, bgAA, sq, tr);
      curr_L = sq->n * 3;

      if (dsq != NULL) free(dsq);
      if ((dsq = malloc(sizeof(ESL_DSQ) * (curr_L + 2))) == NULL) esl_fatal("malloc failed");

      j = 1;
      for (i = 1; i <= sq->n; i++) {
        p7_codontable_GetCodon(codon_table, r, sq->dsq[i], dsq + j);
        j += 3;
      }

      p7_fs_oprofile_ReconfigLength(om_fs5, sq->n);
      p7_fs_ReconfigLength(gm_fs5, sq->n);

      /* Generic reference: p7_Decoding_Frameshift overwrites gx1 with PP  */
      p7_gmx_GrowTo(gx1, M, curr_L, curr_L);
      p7_gmx_GrowTo(gx2, M, curr_L, curr_L);
      p7_gmx_GrowTo(gxpp, M, curr_L, curr_L);
      p7_ivx_GrowTo(iv5, M, p7P_5CODONS);

      p7_GForward_Frameshift (dsq, curr_L, gm_fs5, gx1, iv5, &fsc);
      p7_GBackward_Frameshift(dsq, curr_L, gm_fs5, gx2, iv5, &bsc);
      p7_GDecoding_Frameshift(gm_fs5, gx1, gx2);   /* gx1 is now the PP matrix */

      /*  simd under test: p7_Decoding_Frameshift overwrites fwd with PP  */
      p7_omx_GrowTo_dpf(fwd, M, curr_L, curr_L);
      p7_omx_GrowTo_dpf(bck, M, curr_L, curr_L);

      if (p7_Forward_Frameshift (dsq, curr_L, om_fs5, fwd, &fsc) == eslERANGE) continue;
      if (p7_Backward_Frameshift(dsq, curr_L, om_fs5, fwd, bck, &bsc) == eslERANGE) continue;
      p7_Decoding_Frameshift(om_fs5, fwd, bck);  /* fwd is now the PP matrix */

      /* --- Deconvert SSE PP to generic layout and compare --- */
      p7_omx_FDeconvert(fwd, gxpp);
      if (p7_gmx_Compare(gxpp, gx1, tolerance) != eslOK) esl_fatal(msg);
    }

  free(dsq);
  esl_sq_Destroy(sq);
  p7_hmm_Destroy(hmm);
  p7_gmx_Destroy(gx1);
  p7_gmx_Destroy(gx2);
  p7_ivx_Destroy(iv5);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bck);
  p7_trace_Destroy(tr);
  p7_profile_Destroy(gm);
  p7_profile_fs_Destroy(gm_fs5);
  p7_fs_oprofile_Destroy(om_fs5);
}

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
  P7_GMX         *fgx    = p7_gmx_Create(M, PARSER_ROWS_FWD, M, p7G_NSCELLS);
  P7_GMX         *bgx    = p7_gmx_Create(M, PARSER_ROWS_BWD, M, p7G_NSCELLS);
  P7_IVX         *iv     = p7_ivx_Create(M, p7P_3CODONS);
  P7_DOMAINDEF   *oddef  = NULL; 
  P7_DOMAINDEF   *gddef  = NULL;
  float tolerance;
  int status;

  float fsc3, bsc3;
  float generic_fsc3, generic_bsc3;
  
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

      if (p7_ForwardParser_Frameshift_3Codons(dsq, curr_L, om_fs3, fwd, &fsc3)          == eslERANGE) continue;
      if (p7_BackwardParser_Frameshift_3Codons(dsq, curr_L, om_fs3, fwd, bwd, &generic_fsc3) == eslERANGE) continue;

      oddef->mocc = oddef->btot = oddef->etot = NULL;
      ESL_ALLOC(oddef->mocc, sizeof(float) * (curr_L+1));
      ESL_ALLOC(oddef->btot, sizeof(float) * (curr_L+1));
      ESL_ALLOC(oddef->etot, sizeof(float) * (curr_L+1));
      
      p7_DomainDecoding_Frameshift(om_fs3, fwd, bwd, oddef);

      p7_gmx_GrowTo(fgx, M, PARSER_ROWS_FWD, curr_L);
      p7_gmx_GrowTo(bgx, M, PARSER_ROWS_BWD, curr_L);

      p7_GForwardParser_Frameshift_3Codons(dsq, curr_L, gm_fs3, fgx, iv, &bsc3);
      p7_GBackwardParser_Frameshift_3Codons(dsq, curr_L, gm_fs3, bgx, iv, &generic_bsc3);

      gddef->mocc = gddef->btot = gddef->etot = NULL;
      ESL_ALLOC(gddef->mocc, sizeof(float) * (curr_L+1));
      ESL_ALLOC(gddef->btot, sizeof(float) * (curr_L+1));
      ESL_ALLOC(gddef->etot, sizeof(float) * (curr_L+1));
      p7_GDomainDecoding_Frameshift(gm_fs3, fgx, bgx, gddef);
     
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
#endif /*p7DECODING_FS_TESTDRIVE*/


/*--------------------- end, unit tests -------------------------*/




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7DECODING_FS_TESTDRIVE
/* 
  gcc -g -Wall -msse2 -std=gnu99 -o decoding_fs_utest -I.. -L.. -I../../easel -L../../easel -Dp7DECODING_FS_TESTDRIVE decoding_fs.c -lhmmer -leasel -lm
 ./decoding_fs_utest
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

  p7_FLogsumInit();

  utest_domdef      (r, abcAA, abcDNA, gcode, bgAA, bgDNA, ct, M, N);
  utest_decoding_fs (r, abcAA, abcDNA, gcode, bgAA, bgDNA, ct, M, N);

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

#endif /*p7DECODING_FS_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/


