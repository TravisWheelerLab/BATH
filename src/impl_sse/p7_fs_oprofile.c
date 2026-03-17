/* Routines for the P7_FS_OPROFILE structure:  
 * a search profile in an optimized implementation of framwshift aware search.
 * 
 * Contents:
 *   1. The P7_FS_OPROFILE object: allocation, initialization, destruction.
 *   2. Conversion from generic P7_FS_PROFILE to optimized P7_FS_OPROFILE
 *   3. Conversion from optimized P7_FS_OPROFILE to compact score arrays
 *   4. Debugging and development utilities.
 *   5. Benchmark driver.
 *   6. Unit tests.
 *   7. Test driver.
 *   8. Example.
 */
#include "p7_config.h"

#include <stdio.h>
#include <string.h>
#include <math.h>		/* roundf() */

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_random.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_sse.h"

/*********************************************************************
 * 1. The P7_FS_OPROFILE structure: a frameshift aware score profile.
 *********************************************************************/

/* Function:  p7_fs_oprofile_Create()
 * Synopsis:  Allocate an optimized frameshift profile structure.
 *
 * Purpose:   Allocate a P7_FS_OPROFILE for profiles of up to <allocM>
 *            nodes, digital alphabet <abc>, and <codon_lengths> codon
 *            length variants (1, 3, or 5).
 *
 *            The emission table <rfv> has one row per codon/quasicodon
 *            type plus one row per alphabet symbol:
 *              codon_lengths == 1: p7P_MAXCODONS1 + abc->Kp rows
 *              codon_lengths == 3: p7P_MAXCODONS3 + abc->Kp rows
 *              codon_lengths == 5: p7P_MAXCODONS5 + abc->Kp rows
 *
 * Throws:    <NULL> on allocation error or invalid <codon_lengths>.
 */
P7_FS_OPROFILE *
p7_fs_oprofile_Create(int allocM, const ESL_ALPHABET *abc, int codon_lengths)
{
  int              status;
  P7_FS_OPROFILE  *om_fs = NULL;
  int              nqf;          /* # of float vectors per stripe      */
  int              ncodon_rows;  /* # of rows in rfv emission table    */
  int              c;

  if      (codon_lengths == 1) ncodon_rows = p7P_MAXCODONS1 + abc->Kp;
  else if (codon_lengths == 3) ncodon_rows = p7P_MAXCODONS3 + abc->Kp;
  else if (codon_lengths == 5) ncodon_rows = p7P_MAXCODONS5 + abc->Kp;
  else    ESL_XEXCEPTION(eslEINVAL, "codon_lengths must be 1, 3 or 5");

  nqf = p7O_NQF(allocM);

  /* level 0 */
  ESL_ALLOC(om_fs, sizeof(P7_FS_OPROFILE));
  om_fs->rfv_mem = NULL;
  om_fs->tfv_mem = NULL;
  om_fs->rfv     = NULL;
  om_fs->tfv     = NULL;
  om_fs->clone   = 0;

  /* level 1: raw memory (+15 bytes for manual 16-byte alignment) */
  ESL_ALLOC(om_fs->rfv_mem, sizeof(__m128) * nqf * ncodon_rows + 15);
  ESL_ALLOC(om_fs->tfv_mem, sizeof(__m128) * nqf * p7O_NTRANS  + 15);

  /* pointer array: one __m128 * per emission row */
  ESL_ALLOC(om_fs->rfv, sizeof(__m128 *) * ncodon_rows);

  /* align vector memory on 16-byte boundaries */
  om_fs->rfv[0] = (__m128 *) (((unsigned long int) om_fs->rfv_mem + 15) & (~0xf));
  om_fs->tfv    = (__m128 *) (((unsigned long int) om_fs->tfv_mem + 15) & (~0xf));

  /* set remaining row pointers */
  for (c = 1; c < ncodon_rows; c++)
    om_fs->rfv[c] = om_fs->rfv[0] + (c * nqf);

  om_fs->allocQ4 = nqf;

  /* Frameshift-specific initializations */
  om_fs->codon_lengths = codon_lengths;
  om_fs->fs            = 0.0f;

  for (c = 0; c < p7_NOFFSETS; c++) om_fs->offs[c]    = -1;
  for (c = 0; c < p7_NEVPARAM; c++) om_fs->evparam[c] = p7_EVPARAM_UNSET;
  for (c = 0; c < p7_NCUTOFFS; c++) om_fs->cutoff[c]  = p7_CUTOFF_UNSET;
  for (c = 0; c < p7_MAXABET;  c++) om_fs->compo[c]   = p7_COMPO_UNSET;

  om_fs->name = NULL;
  om_fs->acc  = NULL;
  om_fs->desc = NULL;

  /* Always allocate RF, MM, CS, consensus annotation; rely on leading \0 to
   * signal unused. Zeroed to suppress valgrind warnings on fwrite() in io. */
  ESL_ALLOC(om_fs->rf,        sizeof(char) * (allocM+2));
  ESL_ALLOC(om_fs->mm,        sizeof(char) * (allocM+2));
  ESL_ALLOC(om_fs->cs,        sizeof(char) * (allocM+2));
  ESL_ALLOC(om_fs->consensus, sizeof(char) * (allocM+2));
  memset(om_fs->rf,        '\0', sizeof(char) * (allocM+2));
  memset(om_fs->mm,        '\0', sizeof(char) * (allocM+2));
  memset(om_fs->cs,        '\0', sizeof(char) * (allocM+2));
  memset(om_fs->consensus, '\0', sizeof(char) * (allocM+2));

  om_fs->abc        = abc;
  om_fs->L          = 0;
  om_fs->M          = 0;
  om_fs->max_length = -1;
  om_fs->allocM     = allocM;
  om_fs->mode       = p7_NO_MODE;
  om_fs->nj         = 0.0f;
  return om_fs;

 ERROR:
  p7_fs_oprofile_Destroy(om_fs);
  return NULL;
}

/* Function:  p7_fs_oprofile_IsLocal()
 * Synopsis:  Returns TRUE if profile is in local alignment mode.
 */
int
p7_fs_oprofile_IsLocal(const P7_FS_OPROFILE *om_fs)
{
  if (om_fs->mode == p7_LOCAL || om_fs->mode == p7_UNILOCAL) return TRUE;
  return FALSE;
}


/* Function:  p7_fs_oprofile_Destroy()
 * Synopsis:  Frees an optimized frameshift profile structure.
 */
void
p7_fs_oprofile_Destroy(P7_FS_OPROFILE *om_fs)
{
  if (om_fs == NULL) return;

  if (om_fs->clone == 0)
    {
      if (om_fs->rfv_mem   != NULL) free(om_fs->rfv_mem);
      if (om_fs->tfv_mem   != NULL) free(om_fs->tfv_mem);
      if (om_fs->rfv       != NULL) free(om_fs->rfv);
      if (om_fs->name      != NULL) free(om_fs->name);
      if (om_fs->acc       != NULL) free(om_fs->acc);
      if (om_fs->desc      != NULL) free(om_fs->desc);
      if (om_fs->rf        != NULL) free(om_fs->rf);
      if (om_fs->mm        != NULL) free(om_fs->mm);
      if (om_fs->cs        != NULL) free(om_fs->cs);
      if (om_fs->consensus != NULL) free(om_fs->consensus);
    }

  free(om_fs);
}

/* Function:  p7_fs_oprofile_Clone()
 * Synopsis:  Allocate a cloned copy of an optimized frameshift profile.
 *
 * Purpose:   Shallow copy of <om_fs> for use in multiple threads. The
 *            clone shares all allocated memory with the original; no
 *            vector data is reallocated. <p7_fs_oprofile_Destroy()> on
 *            the clone frees only the shell struct, not the shared memory.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_FS_OPROFILE *
p7_fs_oprofile_Clone(const P7_FS_OPROFILE *om_fs)
{
  int              status;
  P7_FS_OPROFILE  *om2 = NULL;

  ESL_ALLOC(om2, sizeof(P7_FS_OPROFILE));
  memcpy(om2, om_fs, sizeof(P7_FS_OPROFILE));

  om2->clone = 1;

  return om2;

 ERROR:
  p7_fs_oprofile_Destroy(om2);
  return NULL;
}

/*----------------- end, P7_FS_OPROFILE structure ------------------*/



/*********************************************************************
 * 2. Conversion from generic P7_FS_PROFILE to optimized P7_FS_OPROFILE
 *********************************************************************/

/* fs_fb_conversion()
 *
 * Fills the Forward/Backward float vector parts of <om_fs> from the
 * generic frameshift profile <gm_fs>. Scores are converted from
 * log-odds (nats) to odds ratios (pspace floats) for use in the
 * Forward/Backward algorithm.
 *
 * Emission scores: all codon/quasicodon rows (0..ncodon_rows-1) are
 * striped into rfv[c][q], with positions k=1..M distributed across
 * stripes of 4 in the SSE float vectors.
 *
 * Transition scores: striped into tfv[] identically to P7_OPROFILE's
 * fb_conversion(). BM, MM, IM, DM vectors are rotated by -1 (start
 * from k=0); MD, MI, II are straight; DD follows at the end.
 *
 * Special state (ENJC) transitions: stored as scalar odds ratios in xf[][].
 *
 * Returns:   <eslOK> on success.
 * Throws:    <eslEINVAL> if <om_fs> is too small to hold the conversion.
 */
static int
fs_fb_conversion(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs)
{
  int     M   = gm_fs->M;          /* number of model nodes                                      */
  int     nq  = p7O_NQF(M);        /* number of striped float vectors per row                    */
  int     ncodon_rows;              /* total rows in rfv (codon + amino acid rows)                */
  int     c;                        /* counter over codon/aa emission rows                        */
  int     q;                        /* counter over stripes, 0..nq-1                              */
  int     k;                        /* counter over model nodes, 1..M                             */
  int     kb;                       /* possibly offset base k for rotated transition vectors      */
  int     z;                        /* counter within one SIMD float vector (0..3)                */
  int     t;                        /* counter over transitions p7O_{BM..II}                      */
  int     tg;                       /* transition index into gm_fs->tsc                           */
  int     j;                        /* running index into om_fs->tfv                              */
  union { __m128 v; float x[4]; } tmp; /* for building SIMD vectors element-by-element           */

  if (nq > om_fs->allocQ4) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  if      (om_fs->codon_lengths == 1) ncodon_rows = p7P_MAXCODONS1 + gm_fs->abc->Kp;
  else if (om_fs->codon_lengths == 3) ncodon_rows = p7P_MAXCODONS3 + gm_fs->abc->Kp;
  else if (om_fs->codon_lengths == 5) ncodon_rows = p7P_MAXCODONS5 + gm_fs->abc->Kp;
  else ESL_EXCEPTION(eslEINVAL, "codon_lengths must be 1, 3, or 5");

  /* Striped match emission scores: odds ratios from log-odds scores.
   * gm_fs->rsc[c][k] is the log-odds emission score for codon row c at position k.
   * Positions beyond M are set to -inf before exponentiation (-> 0.0).
   */
  for (c = 0; c < ncodon_rows; c++)
    for (k = 1, q = 0; q < nq; q++, k++)
      {
        for (z = 0; z < 4; z++) tmp.x[z] = (k + z*nq <= M) ? gm_fs->rsc[c][k + z*nq] : -eslINFINITY;
        om_fs->rfv[c][q] = esl_sse_expf(tmp.v);
      }

  /* Transition scores (all but DD), striped and interleaved.
   * BM, MM, IM, DM are rotated by -1: base index kb = k-1, starting from k=0.
   * MD, MI, II are straight:          base index kb = k.
   */
  for (j = 0, k = 1, q = 0; q < nq; q++, k++)
    {
      for (t = p7O_BM; t <= p7O_II; t++)
        {
          switch (t) {
          case p7O_BM: tg = p7P_BM; kb = k-1; break;
          case p7O_MM: tg = p7P_MM; kb = k-1; break;
          case p7O_IM: tg = p7P_IM; kb = k-1; break;
          case p7O_DM: tg = p7P_DM; kb = k-1; break;
          case p7O_MD: tg = p7P_MD; kb = k;   break;
          case p7O_MI: tg = p7P_MI; kb = k;   break;
          case p7O_II: tg = p7P_II; kb = k;   break;
          default:     tg = 0;      kb = k;   break; /* unreachable; suppresses compiler warning */
          }
          for (z = 0; z < 4; z++) tmp.x[z] = (kb + z*nq < M) ? p7P_TSC(gm_fs, kb + z*nq, tg) : -eslINFINITY;
          om_fs->tfv[j++] = esl_sse_expf(tmp.v);
        }
    }

  /* DD transitions follow at the end of tfv (j is already positioned there). */
  for (k = 1, q = 0; q < nq; q++, k++)
    {
      for (z = 0; z < 4; z++) tmp.x[z] = (k + z*nq < M) ? p7P_TSC(gm_fs, k + z*nq, p7P_DD) : -eslINFINITY;
      om_fs->tfv[j++] = esl_sse_expf(tmp.v);
    }

  /* Special state (ENJC) transition costs: pspace float odds ratios. */
  om_fs->xf[p7O_E][p7O_LOOP] = expf(gm_fs->xsc[p7P_E][p7P_LOOP]);
  om_fs->xf[p7O_E][p7O_MOVE] = expf(gm_fs->xsc[p7P_E][p7P_MOVE]);
  om_fs->xf[p7O_N][p7O_LOOP] = expf(gm_fs->xsc[p7P_N][p7P_LOOP]);
  om_fs->xf[p7O_N][p7O_MOVE] = expf(gm_fs->xsc[p7P_N][p7P_MOVE]);
  om_fs->xf[p7O_C][p7O_LOOP] = expf(gm_fs->xsc[p7P_C][p7P_LOOP]);
  om_fs->xf[p7O_C][p7O_MOVE] = expf(gm_fs->xsc[p7P_C][p7P_MOVE]);
  om_fs->xf[p7O_J][p7O_LOOP] = expf(gm_fs->xsc[p7P_J][p7P_LOOP]);
  om_fs->xf[p7O_J][p7O_MOVE] = expf(gm_fs->xsc[p7P_J][p7P_MOVE]);

  return eslOK;
}


/* Function:  p7_fs_oprofile_Convert()
 * Synopsis:  Convert a generic frameshift profile to an optimized one.
 *
 * Purpose:   Convert a generic frameshift profile <gm_fs> to an optimized
 *            profile <om_fs>, where <om_fs> has already been allocated for
 *            at least <gm_fs->M> nodes with the same alphabet and
 *            <codon_lengths> setting.
 *
 *            Sets all emission and transition score vectors, copies all
 *            metadata annotation fields, and records the model's current
 *            length and mode configuration.
 *
 * Args:      gm_fs - generic frameshift profile to convert
 *            om_fs - allocated optimized profile to receive the result
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if the two profiles are incompatible (different
 *            alphabets, or <om_fs> too small for <gm_fs->M> nodes).
 *            <eslEMEM> on allocation failure (string duplication).
 */
int
p7_fs_oprofile_Convert(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs)
{
  int status, z;

  if (gm_fs->abc->type != om_fs->abc->type) ESL_EXCEPTION(eslEINVAL, "alphabets of the two profiles don't match");
  if (gm_fs->M         >  om_fs->allocM)    ESL_EXCEPTION(eslEINVAL, "optimized profile is too small");

  /* Set configuration fields first; fs_fb_conversion() may use them. */
  om_fs->mode          = gm_fs->mode;
  om_fs->L             = gm_fs->L;
  om_fs->M             = gm_fs->M;
  om_fs->nj            = gm_fs->nj;
  om_fs->max_length    = gm_fs->max_length;
  om_fs->codon_lengths = gm_fs->codon_lengths;
  om_fs->fs            = gm_fs->fs;

  if ((status = fs_fb_conversion(gm_fs, om_fs)) != eslOK) return status;

  if (om_fs->name != NULL) free(om_fs->name);
  if (om_fs->acc  != NULL) free(om_fs->acc);
  if (om_fs->desc != NULL) free(om_fs->desc);
  if ((status = esl_strdup(gm_fs->name, -1, &(om_fs->name))) != eslOK) goto ERROR;
  if ((status = esl_strdup(gm_fs->acc,  -1, &(om_fs->acc)))  != eslOK) goto ERROR;
  if ((status = esl_strdup(gm_fs->desc, -1, &(om_fs->desc))) != eslOK) goto ERROR;
  strcpy(om_fs->rf,        gm_fs->rf);
  strcpy(om_fs->mm,        gm_fs->mm);
  strcpy(om_fs->cs,        gm_fs->cs);
  strcpy(om_fs->consensus, gm_fs->consensus);
  for (z = 0; z < p7_NEVPARAM; z++) om_fs->evparam[z] = gm_fs->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) om_fs->cutoff[z]  = gm_fs->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) om_fs->compo[z]   = gm_fs->compo[z];

  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_fs_oprofile_ReconfigLength()
 * Synopsis:  Set the target sequence length of an optimized frameshift profile.
 *
 * Purpose:   Given an already configured profile <om_fs>, quickly reset its
 *            expected length distribution for a new mean target sequence
 *            length of <L>.
 *
 *            Updates only the N, C, J special-state loop/move transition
 *            probabilities in <om_fs->xf>. The E-state transitions are
 *            length-independent and are left unchanged.
 *
 *            This does not affect the null model; call <p7_bg_SetLength()>
 *            separately to keep them synchronized.
 *
 * Args:      om_fs - optimized frameshift profile to reconfigure
 *            L     - new target sequence length
 *
 * Returns:   <eslOK> on success.
 */
int
p7_fs_oprofile_ReconfigLength(P7_FS_OPROFILE *om_fs, int L)
{
  float pmove, ploop;

  pmove = (2.0f + om_fs->nj) / ((float) L + 2.0f + om_fs->nj); /* 2/(L+2) for sw; 3/(L+3) for fs */
  ploop = 1.0f - pmove;

  om_fs->xf[p7O_N][p7O_LOOP] = om_fs->xf[p7O_C][p7O_LOOP] = om_fs->xf[p7O_J][p7O_LOOP] = ploop;
  om_fs->xf[p7O_N][p7O_MOVE] = om_fs->xf[p7O_C][p7O_MOVE] = om_fs->xf[p7O_J][p7O_MOVE] = pmove;

  om_fs->L = L;
  return eslOK;
}


/* Function:  p7_fs_oprofile_ReconfigMultihit()
 * Synopsis:  Quickly reconfig model into multihit mode for target length <L>.
 *
 * Purpose:   Given a profile <om_fs> that's already been configured once,
 *            quickly reconfigure it into a multihit mode for target 
 *            length <L>. 
 *            
 *            This gets called in domain definition, when we need to
 *            flip the model in and out of unihit mode to
 *            process individual domains.
 *            
 * Note:      You can't just flip uni/multi mode alone, because that
 *            parameterization also affects target length
 *            modeling. You need to make sure uni vs. multi choice is
 *            made before the length model is set, and you need to
 *            make sure the length model is recalculated if you change
 *            the uni/multi mode. Hence, these functions call
 *            <p7_fs_oprofile_ReconfigLength()>.
 */
int
p7_fs_oprofile_ReconfigMultihit(P7_FS_OPROFILE *om_fs, int L)
{
  om_fs->xf[p7O_E][p7O_MOVE] = 0.5;
  om_fs->xf[p7O_E][p7O_LOOP] = 0.5;
  om_fs->nj = 1.0f;

  return p7_fs_oprofile_ReconfigLength(om_fs, L);
}

/* Function:  p7_fs_oprofile_ReconfigUnihit()
 * Synopsis:  Quickly reconfig model into unihit mode for target length <L>.
 *
 * Purpose:   Given a profile <om_fs> that's already been configured once,
 *            quickly reconfigure it into a unihit mode for target 
 *            length <L>. 
 *            
 *            This gets called in domain definition, when we need to
 *            flip the model in and out of unihit <L=0> mode to
 *            process individual domains.
 */
int
p7_fs_oprofile_ReconfigUnihit(P7_FS_OPROFILE *om_fs, int L)
{
  om_fs->xf[p7O_E][p7O_MOVE] = 1.0f;
  om_fs->xf[p7O_E][p7O_LOOP] = 0.0f;
  om_fs->nj = 0.0f;

  return p7_fs_oprofile_ReconfigLength(om_fs, L);
}


/*------------ end, conversions to P7_FS_OPROFILE ------------------*/

/***********************************************************************
*   3. Conversion from optimized P7_FS_OPROFILE to compact score arrays
 ***********************************************************************/


/*------------ end, conversions from P7_FS_OPROFILE ------------------*/


/*********************************************************************
 * 4. Debugging and development utilities.
 *********************************************************************/

/*------------ end, P7_FS_OPROFILE debugging tools  ----------------*/



/*****************************************************************
 * 5. Benchmark driver.
 *****************************************************************/

/*---------------- end, benchmark driver ------------------------*/




  
/*****************************************************************
 * 6. Unit tests
 *****************************************************************/

/*------------------- end, unit tests ---------------------------*/




/*****************************************************************
 * 7. Test driver
 *****************************************************************/

/*------------------- end, test driver --------------------------*/


/*****************************************************************
 * 8. Example
 *****************************************************************/
/*----------------------- end, example --------------------------*/


