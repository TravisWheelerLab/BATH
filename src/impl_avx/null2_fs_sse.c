/* Frameshift null2 model; SSE implementation.
 * Ported from impl_sse/null2_fs.c with _sse suffix for runtime dispatch.
 *
 * Contents:
 *   1. p7_Null2_fs_ByExpectation_sse()
 */
#include "p7_config.h"

#ifdef eslENABLE_SSE

#include <stdlib.h>
#include <string.h>

#include <xmmintrin.h>  /* SSE  */
#include <emmintrin.h>  /* SSE2 */

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Null2_fs_ByExpectation_sse()
 *
 * Purpose:   SSE implementation of null2 estimation by expectation for
 *            the frameshift model. Uses the 8-cell FS layout posterior
 *            matrix <pp> (from p7_Decoding_Frameshift_sse()) and the
 *            frameshift profile <om_fs> to fill <null2[0..Kp-1]>.
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <om_fs->codon_lengths> is not 1, 3, or 5.
 */
int
p7_Null2_fs_ByExpectation_sse(const P7_FS_OPROFILE *om_fs, P7_OMX *pp, float *null2)
{
  int      M    = om_fs->M;
  int      Ld   = pp->L;
  int      Q    = p7O_NQF(M);
  float   *xmx  = pp->xmx;   /* enables use of XMXo(i,s) macro */
  float    norm;
  __m128  *rp;
  __m128   sv;
  float    xfactor;
  int      amino_offset;
  int      i, q, x;

  /* Determine the amino acid emission row offset in rfv[] */
  if      (om_fs->codon_lengths == 5) amino_offset = p7P_MAXCODONS5;
  else if (om_fs->codon_lengths == 3) amino_offset = p7P_MAXCODONS3;
  else if (om_fs->codon_lengths == 1) amino_offset = p7P_MAXCODONS1;
  else    ESL_EXCEPTION(eslEINVAL, "codon_lengths must be 1, 3, or 5");

  /* Calculate expected # of times that each emitting state was used
   * in generating the Ld residues in this domain.
   * The 0 row in <pp> is used to hold these numbers.
   * We only need M_C0 (total match) and I; copy full row for simplicity.
   */
  memcpy(pp->dpf[0], pp->dpf[1], sizeof(__m128) * p7X_NSCELLS_FS * Q);
  XMXo(0, p7X_N) = XMXo(1, p7X_N);
  XMXo(0, p7X_C) = XMXo(1, p7X_C);  /* 0.0 */
  XMXo(0, p7X_J) = XMXo(1, p7X_J);  /* 0.0 */

  for (i = 2; i <= Ld; i++)
    {
      for (q = 0; q < Q; q++)
        {
          MMO_FS(pp->dpf[0], q, p7X_FS_C0) = _mm_add_ps(MMO_FS(pp->dpf[i], q, p7X_FS_C0),
                                                          MMO_FS(pp->dpf[0], q, p7X_FS_C0));
          IMO_FS(pp->dpf[0], q)             = _mm_add_ps(IMO_FS(pp->dpf[i], q),
                                                          IMO_FS(pp->dpf[0], q));
        }
      XMXo(0, p7X_N) += XMXo(i, p7X_N);
      XMXo(0, p7X_C) += XMXo(i, p7X_C);
      XMXo(0, p7X_J) += XMXo(i, p7X_J);
    }

  /* Convert those expected #'s to frequencies, to use as posterior weights. */
  norm = 1.0 / (float) Ld;
  sv   = _mm_set1_ps(norm);
  for (q = 0; q < Q; q++)
    {
      MMO_FS(pp->dpf[0], q, p7X_FS_C0) = _mm_mul_ps(MMO_FS(pp->dpf[0], q, p7X_FS_C0), sv);
      IMO_FS(pp->dpf[0], q)             = _mm_mul_ps(IMO_FS(pp->dpf[0], q),             sv);
    }
  XMXo(0, p7X_N) *= norm;
  XMXo(0, p7X_C) *= norm;
  XMXo(0, p7X_J) *= norm;

  /* Calculate null2's emission odds, by taking posterior weighted sum
   * over all emission vectors used in paths explaining the domain.
   * Amino acid emission odds are in rfv[amino_offset + x].
   * Insert emission odds are implicitly 1.0.
   */
  xfactor = XMXo(0, p7X_N) + XMXo(0, p7X_C) + XMXo(0, p7X_J);
  for (x = 0; x < om_fs->abc->K; x++)
    {
      sv = _mm_setzero_ps();
      rp = om_fs->rfv[amino_offset + x];
      for (q = 0; q < Q; q++)
        {
          sv = _mm_add_ps(sv, _mm_mul_ps(MMO_FS(pp->dpf[0], q, p7X_FS_C0), *rp)); rp++;
          sv = _mm_add_ps(sv, IMO_FS(pp->dpf[0], q));  /* insert odds implicitly 1.0 */
        }
      esl_sse_hsum_ps(sv, &(null2[x]));
      null2[x] += xfactor;
    }

  /* make valid scores for all degeneracies, by averaging the odds ratios. */
  esl_abc_FAvgScVec(om_fs->abc, null2);
  null2[om_fs->abc->K]    = 1.0;        /* gap character    */
  null2[om_fs->abc->Kp-2] = 1.0;        /* nonresidue "*"   */
  null2[om_fs->abc->Kp-1] = 1.0;        /* missing data "~" */

  return eslOK;
}

#endif /* eslENABLE_SSE */
