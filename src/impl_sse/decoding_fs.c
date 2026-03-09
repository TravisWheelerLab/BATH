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
  float *sfwd = NULL;   /* sfwd[i] = product of fwd scale factors at rows 0..i  */
  float *sbck = NULL;   /* sbck[i] = product of bck scale factors at rows i..L  */
  float  bck_total;     /* N(0)+N(1)+N(2) from backward = total bck probability */
  float  inv_Z;         /* 1 / (bck_total * sbck[0])                            */
  float  njcp;
  int    i;
  int    status;

  ESL_ALLOC(sfwd, sizeof(float) * (L+2));
  ESL_ALLOC(sbck, sizeof(float) * (L+2));

  /* Cumulative forward scale products: sfwd[i] = prod(fwd_scale_0 .. fwd_scale_i) */
  sfwd[0] = oxf->xmx[0*p7X_NXCELLS + p7X_SCALE];
  for (i = 1; i <= L; i++)
    sfwd[i] = sfwd[i-1] * oxf->xmx[i*p7X_NXCELLS + p7X_SCALE];

  /* Cumulative backward scale products: sbck[i] = prod(bck_scale_i .. bck_scale_L)
   * When has_own_scales=FALSE, bck uses the fwd scale factors; when TRUE, bck has
   * its own.  Either way, oxb->xmx[i*p7X_NXCELLS+p7X_SCALE] holds the right value. */
  sbck[L+1] = 1.0f;
  for (i = L; i >= 0; i--)
    sbck[i] = sbck[i+1] * oxb->xmx[i*p7X_NXCELLS + p7X_SCALE];

  /* Normalization: bck_total = N(0)+N(1)+N(2) from backward, each divided by sbck[0].
   * inv_Z = 1/Z where Z = bck_total * sbck[0] = true total probability. */
  bck_total = oxb->xmx[0*p7X_NXCELLS + p7X_N]
            + oxb->xmx[1*p7X_NXCELLS + p7X_N]
            + oxb->xmx[2*p7X_NXCELLS + p7X_N];
  inv_Z = 1.0f / (bck_total * sbck[0]);

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
          * sfwd[i-3] * sbck[i-3] * inv_Z;

      /* etot[i]: cumulative expected E-state exits up to position i.
       * E(i) contributes: fwd_E(i) * bck_E(i) (same-row product). */
      ddef->etot[i] = ddef->etot[i-3]
        + oxf->xmx[i*p7X_NXCELLS + p7X_E] * oxb->xmx[i*p7X_NXCELLS + p7X_E]
          * sfwd[i] * sbck[i] * inv_Z;

      /* mocc[i]: probability residue i is in the core model = 1 - P(i in N/J/C).
       *
       * For each of N, J, C: three codon-length offsets contribute to residue i.
       * Codon ending at i:     fwd(i-3) * T_XL * bck(i),   rows (i-3, i)
       * Codon ending at i+1:   fwd(i-2) * T_XL * bck(i+1), rows (i-2, i+1)
       * Codon ending at i+2:   fwd(i-1) * T_XL * bck(i+2), rows (i-1, i+2)
       *
       * General cross-row posterior: fwd_stored(j) * bck_stored(k) * sfwd[j] * sbck[k] * inv_Z
       */
      njcp = 0.;

      /* N state */
      njcp += oxf->xmx[(i-3)*p7X_NXCELLS + p7X_N] * oxb->xmx[i*p7X_NXCELLS + p7X_N]
              * om_fs->xf[p7O_N][p7O_LOOP] * sfwd[i-3] * sbck[i] * inv_Z;
      if (i < L)
        njcp += oxf->xmx[(i-2)*p7X_NXCELLS + p7X_N] * oxb->xmx[(i+1)*p7X_NXCELLS + p7X_N]
                * om_fs->xf[p7O_N][p7O_LOOP] * sfwd[i-2] * sbck[i+1] * inv_Z;
      if (i < L-1)
        njcp += oxf->xmx[(i-1)*p7X_NXCELLS + p7X_N] * oxb->xmx[(i+2)*p7X_NXCELLS + p7X_N]
                * om_fs->xf[p7O_N][p7O_LOOP] * sfwd[i-1] * sbck[i+2] * inv_Z;

      /* J state */
      njcp += oxf->xmx[(i-3)*p7X_NXCELLS + p7X_J] * oxb->xmx[i*p7X_NXCELLS + p7X_J]
              * om_fs->xf[p7O_J][p7O_LOOP] * sfwd[i-3] * sbck[i] * inv_Z;
      if (i < L)
        njcp += oxf->xmx[(i-2)*p7X_NXCELLS + p7X_J] * oxb->xmx[(i+1)*p7X_NXCELLS + p7X_J]
                * om_fs->xf[p7O_J][p7O_LOOP] * sfwd[i-2] * sbck[i+1] * inv_Z;
      if (i < L-1)
        njcp += oxf->xmx[(i-1)*p7X_NXCELLS + p7X_J] * oxb->xmx[(i+2)*p7X_NXCELLS + p7X_J]
                * om_fs->xf[p7O_J][p7O_LOOP] * sfwd[i-1] * sbck[i+2] * inv_Z;

      /* C state */
      njcp += oxf->xmx[(i-3)*p7X_NXCELLS + p7X_C] * oxb->xmx[i*p7X_NXCELLS + p7X_C]
              * om_fs->xf[p7O_C][p7O_LOOP] * sfwd[i-3] * sbck[i] * inv_Z;
      if (i < L)
        njcp += oxf->xmx[(i-2)*p7X_NXCELLS + p7X_C] * oxb->xmx[(i+1)*p7X_NXCELLS + p7X_C]
                * om_fs->xf[p7O_C][p7O_LOOP] * sfwd[i-2] * sbck[i+1] * inv_Z;
      if (i < L-1)
        njcp += oxf->xmx[(i-1)*p7X_NXCELLS + p7X_C] * oxb->xmx[(i+2)*p7X_NXCELLS + p7X_C]
                * om_fs->xf[p7O_C][p7O_LOOP] * sfwd[i-1] * sbck[i+2] * inv_Z;

      ddef->mocc[i] = 1. - njcp;
    }

  ddef->L = L;
  free(sfwd);
  free(sbck);
  return eslOK;

 ERROR:
  if (sfwd) free(sfwd);
  if (sbck) free(sbck);
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
/*--------------------- end, unit tests -------------------------*/




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
/*-------------------- end, test driver -------------------------*/


