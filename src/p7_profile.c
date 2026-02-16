/* Routines for the P7_PROFILE structure - Plan 7's search profile
 *                                         
 *    1. The P7_PROFILE object: allocation, initialization, destruction.
 *    2. Access methods.
 *    3. Debugging and development code.
 *    4. Unit tests.
 *    5. Test driver.
 *
 * See also: 
 *   modelconfig.c : routines that configure a profile given an HMM
 */

#include "p7_config.h"

#include <string.h>
#ifdef HMMER_MPI
#include <mpi.h>
#endif

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"


/*****************************************************************
 * 1. The P7_PROFILE object: allocation, initialization, destruction.
 *****************************************************************/

/* Function:  p7_profile_Create()
 * Synopsis:  Allocates a profile.
 *
 * Purpose:   Allocates for a profile of up to <M> nodes, for digital
 *            alphabet <abc>.
 *            
 *            Because this function might be in the critical path (in
 *            hmmscan, for example), we leave much of the model
 *            unintialized, including scores and length model
 *            probabilities. The <p7_ProfileConfig()> call is what
 *            sets these. 
 *            
 *            The alignment mode is set to <p7_NO_MODE>.  The
 *            reference pointer <gm->abc> is set to <abc>.
 *
 * Returns:   a pointer to the new profile.
 *
 * Throws:    <NULL> on allocation error.
 *
 * Xref:      STL11/125.
 */
P7_PROFILE *
p7_profile_Create(int allocM, const ESL_ALPHABET *abc)
{
  P7_PROFILE *gm = NULL;
  int         x;
  int         status;

  /* level 0 */
  ESL_ALLOC(gm, sizeof(P7_PROFILE));
  gm->tsc       = NULL;
  gm->rsc       = NULL;
  gm->rf        = NULL;
  gm->mm        = NULL;
  gm->cs        = NULL;
  gm->consensus = NULL;
  
  /* level 1 */
  ESL_ALLOC(gm->tsc,       sizeof(float)   * allocM * p7P_NTRANS); 
  ESL_ALLOC(gm->rsc,       sizeof(float *) * abc->Kp);
  ESL_ALLOC(gm->rf,        sizeof(char)    * (allocM+2)); /* yes, +2: each is (0)1..M, +trailing \0  */
  ESL_ALLOC(gm->mm,        sizeof(char)    * (allocM+2));
  ESL_ALLOC(gm->cs,        sizeof(char)    * (allocM+2));
  ESL_ALLOC(gm->consensus, sizeof(char)    * (allocM+2));
  gm->rsc[0] = NULL;
  
  /* level 2 */
  ESL_ALLOC(gm->rsc[0], sizeof(float) * abc->Kp * (allocM+1) * p7P_NR);
  for (x = 1; x < abc->Kp; x++) 
    gm->rsc[x] = gm->rsc[0] + x * (allocM+1) * p7P_NR;

  /* Initialize some edge pieces of memory that are never used,
   * and are only present for indexing convenience.
   */
  esl_vec_FSet(gm->tsc, p7P_NTRANS, -eslINFINITY);     /* node 0 nonexistent, has no transitions  */
  if (allocM > 1) {
    p7P_TSC(gm, 1, p7P_DM) = -eslINFINITY;             /* delete state D_1 is wing-retracted      */
    p7P_TSC(gm, 1, p7P_DD) = -eslINFINITY;
  }
  for (x = 0; x < abc->Kp; x++) {        
    p7P_MSC(gm, 0,      x) = -eslINFINITY;             /* no emissions from nonexistent M_0... */
    p7P_ISC(gm, 0,      x) = -eslINFINITY;             /* or I_0... */
    /* I_M is initialized in profile config, when we know actual M, not just allocated max M   */
  }
  x = esl_abc_XGetGap(abc);                         /* no emission can emit/score gap characters */
  esl_vec_FSet(gm->rsc[x], (allocM+1)*p7P_NR, -eslINFINITY);
  x = esl_abc_XGetMissing(abc);                        /* no emission can emit/score missing data characters */
  esl_vec_FSet(gm->rsc[x], (allocM+1)*p7P_NR, -eslINFINITY);

  /* Set remaining info  */
  gm->mode             = p7_NO_MODE;
  gm->L                = 0;
  gm->allocM           = allocM;
  gm->M                = 0;
  gm->max_length       = -1;
  gm->nj               = 0.0f;

  gm->roff             = -1;
  gm->eoff             = -1;
  gm->offs[p7_MOFFSET] = -1;
  gm->offs[p7_FOFFSET] = -1;
  gm->offs[p7_POFFSET] = -1;

  gm->name             = NULL;
  gm->acc              = NULL;
  gm->desc             = NULL;
  gm->rf[0]            = 0;     /* RF line is optional annotation; this flags that it's not set yet */
  gm->mm[0]            = 0;     /* likewise for MM annotation line */
  gm->cs[0]            = 0;     /* likewise for CS annotation line */
  gm->consensus[0]     = 0;
  
  for (x = 0; x < p7_NEVPARAM; x++) gm->evparam[x] = p7_EVPARAM_UNSET;
  for (x = 0; x < p7_NCUTOFFS; x++) gm->cutoff[x]  = p7_CUTOFF_UNSET;
  for (x = 0; x < p7_MAXABET;  x++) gm->compo[x]   = p7_COMPO_UNSET;

  gm->abc         = abc;
  return gm;

 ERROR:
  p7_profile_Destroy(gm);
  return NULL;
}

/* Function:  p7_profile_fs5_Create()
 * Synopsis:  Allocates a 5 codon length frameshift aware profile.
 *
 * Purpose:   Allocates for a 5 codon length frameshift profile of up 
 *            to <M> nodes, for digital alphabet <abc>.
 *
 *            We leave much of the model unintialized, including 
 *            scores and length model probabilities. 
 *            The <p7_ProfileConfig_fs5()> call is what sets these.
 *
 *            The alignment mode is set to <p7_NO_MODE>.  The
 *            reference pointer <gm->abc> is set to <abc>.
 *
 * Returns:   a pointer to the new profile.
 *
 * Throws:    <NULL> on allocation error.
 *
 * p7_profile_fs5_CreateXref:      STL11/125.
 */
P7_FS_PROFILE *
p7_profile_fs5_Create(int allocM, const ESL_ALPHABET *abc)
{
  P7_FS_PROFILE *gm_fs5 = NULL;
  int         x;
  int         status;

  /* level 0 */
  ESL_ALLOC(gm_fs5, sizeof(P7_FS_PROFILE));
  gm_fs5->tsc       = NULL;
  gm_fs5->rsc       = NULL;
  gm_fs5->codons    = NULL;
  gm_fs5->indel_pos = NULL;
  gm_fs5->rf        = NULL;
  gm_fs5->mm        = NULL;
  gm_fs5->cs        = NULL;
  gm_fs5->consensus = NULL;

  /* level 1 */
  ESL_ALLOC(gm_fs5->tsc,       sizeof(float)     * allocM * p7P_NTRANS);
  ESL_ALLOC(gm_fs5->rsc,       sizeof(float *)   * (allocM+1));
  ESL_ALLOC(gm_fs5->codons,    sizeof(ESL_DSQ *) * (allocM+1));
  ESL_ALLOC(gm_fs5->indel_pos, sizeof(ESL_DSQ *) * (allocM+1));
  ESL_ALLOC(gm_fs5->rf,        sizeof(char)      * (allocM+2)); /* yes, +2: each is (0)1..M, +trailing \0  */
  ESL_ALLOC(gm_fs5->mm,        sizeof(char)      * (allocM+2));
  ESL_ALLOC(gm_fs5->cs,        sizeof(char)      * (allocM+2));
  ESL_ALLOC(gm_fs5->consensus, sizeof(char)      * (allocM+2));
  gm_fs5->rsc[0]       = NULL;
  gm_fs5->codons[0]    = NULL;
  gm_fs5->indel_pos[0] = NULL;


  /* level 2 */
  ESL_ALLOC(gm_fs5->rsc[0], sizeof(float) * (allocM+1) * (p7P_MAXCODONS5 + abc->Kp));

  for (x = 1; x <= allocM; x++)   
    gm_fs5->rsc[x] = gm_fs5->rsc[0] + x * (p7P_MAXCODONS5 + abc->Kp);

  ESL_ALLOC(gm_fs5->codons[0], sizeof(ESL_DSQ) * (allocM+1) * (p7P_MAXCODONS5+1)); /* +1 for trailing \0 */

  for (x = 1; x <= allocM; x++)
    gm_fs5->codons[x] = gm_fs5->codons[0] + x * p7P_MAXCODONS5;

  ESL_ALLOC(gm_fs5->indel_pos[0], sizeof(ESL_DSQ) * (allocM+1) * (p7P_MAXCODONS5+1)); /* +1 for trailing \0 */

  for (x = 1; x <= allocM; x++)
    gm_fs5->indel_pos[x] = gm_fs5->indel_pos[0] + x * p7P_MAXCODONS5;

  /* Initialize some edge pieces of memory that are never used,
   * and are only present for indexing convenience. */
  esl_vec_FSet(gm_fs5->tsc, p7P_NTRANS, -eslINFINITY);     /* node 0 nonexistent, has no transitions  */
  if (allocM > 1) {
    p7P_TSC(gm_fs5, 1, p7P_DM) = -eslINFINITY;             /* delete state D_1 is wing-retracted      */
    p7P_TSC(gm_fs5, 1, p7P_DD) = -eslINFINITY;
  }

  for (x = 0; x < (p7P_MAXCODONS5 + abc->Kp); x++) 
    p7P_MSC_CODON(gm_fs5, 0,      x) = -eslINFINITY;             /* no emissions from nonexistent M_0... */
  
  /* Set remaining info  */
  gm_fs5->mode             = p7_NO_MODE;
  gm_fs5->L                = 0;
  gm_fs5->allocM           = allocM;
  gm_fs5->M                = 0;
  gm_fs5->max_length       = -1;
  gm_fs5->nj               = 0.0f;
  gm_fs5->abc              = abc;

  gm_fs5->roff             = -1;
  gm_fs5->eoff             = -1;
  gm_fs5->offs[p7_MOFFSET] = -1;
  gm_fs5->offs[p7_FOFFSET] = -1;
  gm_fs5->offs[p7_POFFSET] = -1;

  gm_fs5->codon_lengths    = 5;
  gm_fs5->fs               = 0.;

  gm_fs5->name             = NULL;
  gm_fs5->acc              = NULL;
  gm_fs5->desc             = NULL;
  gm_fs5->rf[0]            = 0;     /* RF line is optional annotation; this flags that it's not set yet */
  gm_fs5->mm[0]            = 0;     /* likewise for MM annotation line */
  gm_fs5->cs[0]            = 0;     /* likewise for CS annotation line */
  gm_fs5->consensus[0]     = 0;

  for (x = 0; x < p7_NEVPARAM; x++) gm_fs5->evparam[x] = p7_EVPARAM_UNSET;
  for (x = 0; x < p7_NCUTOFFS; x++) gm_fs5->cutoff[x]  = p7_CUTOFF_UNSET;
  for (x = 0; x < p7_MAXABET;  x++) gm_fs5->compo[x]   = p7_COMPO_UNSET;

  return gm_fs5;

 ERROR:
  p7_profile_fs_Destroy(gm_fs5);
  return NULL;
}


/* Function:  p7_profile_fs3_Create()
 * Synopsis:  Allocates a 3 codon length frameshift aware profile.
 *
 * Purpose:   Allocates for a 3 codon length frameshift profile of up
 *            to <M> nodes, for digital alphabet <abc>.
 *
 *            We leave much of the model unintialized, including
 *            scores and length model probabilities.
 *            The <p7_ProfileConfig_fs5()> call is what sets these.
 * p7_profile_fs5_CreateXref:      STL11/125.
 */
P7_FS_PROFILE *
p7_profile_fs3_Create(int allocM, const ESL_ALPHABET *abc)
{
  P7_FS_PROFILE *gm_fs3 = NULL;
  int         x;
  int         status;

  /* level 0 */
  ESL_ALLOC(gm_fs3, sizeof(P7_FS_PROFILE));
  gm_fs3->tsc       = NULL;
  gm_fs3->rsc       = NULL;
  gm_fs3->codons    = NULL;
  gm_fs3->indel_pos = NULL;
  gm_fs3->rf        = NULL;
  gm_fs3->mm        = NULL;
  gm_fs3->cs        = NULL;
  gm_fs3->consensus = NULL;

  /* level 1 */
  ESL_ALLOC(gm_fs3->tsc,       sizeof(float)     * allocM * p7P_NTRANS);
  ESL_ALLOC(gm_fs3->rsc,       sizeof(float *)   * (allocM+1));
  ESL_ALLOC(gm_fs3->codons,    sizeof(ESL_DSQ *) * (allocM+1));
  ESL_ALLOC(gm_fs3->indel_pos, sizeof(ESL_DSQ *) * (allocM+1));
  ESL_ALLOC(gm_fs3->rf,        sizeof(char)      * (allocM+2)); /* yes, +2: each is (0)1..M, +trailing \0  */
  ESL_ALLOC(gm_fs3->mm,        sizeof(char)      * (allocM+2));
  ESL_ALLOC(gm_fs3->cs,        sizeof(char)      * (allocM+2));
  ESL_ALLOC(gm_fs3->consensus, sizeof(char)      * (allocM+2));
  gm_fs3->rsc[0]       = NULL;
  gm_fs3->codons[0]    = NULL;
  gm_fs3->indel_pos[0] = NULL;


  /* level 2 */
  ESL_ALLOC(gm_fs3->rsc[0], sizeof(float) * (allocM+1) * (p7P_MAXCODONS3 + abc->Kp));

  for (x = 1; x <= allocM; x++)   
    gm_fs3->rsc[x] = gm_fs3->rsc[0] + x * (p7P_MAXCODONS3 + abc->Kp);

  ESL_ALLOC(gm_fs3->codons[0], sizeof(ESL_DSQ) * (allocM+1) * (p7P_MAXCODONS3+1)); /* +1 for trailing \0 */

  for (x = 1; x <= allocM; x++)
    gm_fs3->codons[x] = gm_fs3->codons[0] + x * p7P_MAXCODONS3;

  ESL_ALLOC(gm_fs3->indel_pos[0], sizeof(ESL_DSQ) * (allocM+1) * (p7P_MAXCODONS3+1)); /* +1 for trailing \0 */

  for (x = 1; x <= allocM; x++)
    gm_fs3->indel_pos[x] = gm_fs3->indel_pos[0] + x * p7P_MAXCODONS3;

  /* Initialize some edge pieces of memory that are never used,
   * and are only present for indexing convenience. */
  esl_vec_FSet(gm_fs3->tsc, p7P_NTRANS, -eslINFINITY);     /* node 0 nonexistent, has no transitions  */
  if (allocM > 1) {
    p7P_TSC(gm_fs3, 1, p7P_DM) = -eslINFINITY;             /* delete state D_1 is wing-retracted      */
    p7P_TSC(gm_fs3, 1, p7P_DD) = -eslINFINITY;
  }

  for (x = 0; x < (p7P_MAXCODONS3 + abc->Kp); x++) 
    p7P_MSC_CODON(gm_fs3, 0,      x) = -eslINFINITY;             /* no emissions from nonexistent M_0... */
  
  /* Set remaining info  */
  gm_fs3->mode             = p7_NO_MODE;
  gm_fs3->L                = 0;
  gm_fs3->allocM           = allocM;
  gm_fs3->M                = 0;
  gm_fs3->max_length       = -1;
  gm_fs3->nj               = 0.0f;
  gm_fs3->abc              = abc;

  gm_fs3->roff             = -1;
  gm_fs3->eoff             = -1;
  gm_fs3->offs[p7_MOFFSET] = -1;
  gm_fs3->offs[p7_FOFFSET] = -1;
  gm_fs3->offs[p7_POFFSET] = -1;
  
  gm_fs3->codon_lengths    = 3;
  gm_fs3->fs               = 0.; 

  gm_fs3->name             = NULL;
  gm_fs3->acc              = NULL;
  gm_fs3->desc             = NULL;
  gm_fs3->rf[0]            = 0;     /* RF line is optional annotation; this flags that it's not set yet */
  gm_fs3->mm[0]            = 0;     /* likewise for MM annotation line */
  gm_fs3->cs[0]            = 0;     /* likewise for CS annotation line */
  gm_fs3->consensus[0]     = 0;

  for (x = 0; x < p7_NEVPARAM; x++) gm_fs3->evparam[x] = p7_EVPARAM_UNSET;
  for (x = 0; x < p7_NCUTOFFS; x++) gm_fs3->cutoff[x]  = p7_CUTOFF_UNSET;
  for (x = 0; x < p7_MAXABET;  x++) gm_fs3->compo[x]   = p7_COMPO_UNSET;

  return gm_fs3;

 ERROR:
  p7_profile_fs_Destroy(gm_fs3);
  return NULL;
}
  

/* Function:  p7_profile_tr_Create()
 * Synopsis:  Allocates a translated profile.
 *
 * Purpose:   Allocates a translated profile of up to <M> nodes, 
 *            for digital alphabet <abc>.
 *
 *            We leave much of the model unintialized, including 
 *            scores and length model probabilities. The 
 *            <p7_ProfileConfig_tr()> call is what sets these.
 *
 *            The alignment mode is set to <p7_NO_MODE>.  The
 *            reference pointer <gm->abc> is set to <abc>.
 *
 * Returns:   a pointer to the new profile.
 *
 * Throws:    <NULL> on allocation error.
 *
 * Xref:      STL11/125.
 */
P7_FS_PROFILE *
p7_profile_tr_Create(int allocM, const ESL_ALPHABET *abc)
{

  P7_FS_PROFILE *gm_tr = NULL;
  int         x;
  int         status;

  /* level 0 */
  ESL_ALLOC(gm_tr, sizeof(P7_FS_PROFILE));
  gm_tr->tsc       = NULL;
  gm_tr->rsc       = NULL;
  gm_tr->codons    = NULL;
  gm_tr->indel_pos = NULL;
  gm_tr->rf        = NULL;
  gm_tr->mm        = NULL;
  gm_tr->cs        = NULL;
  gm_tr->consensus = NULL;

  /* level 1 */
  ESL_ALLOC(gm_tr->tsc,       sizeof(float)     * allocM * p7P_NTRANS);
  ESL_ALLOC(gm_tr->rsc,       sizeof(float *)   * (allocM+1));
  ESL_ALLOC(gm_tr->codons,    sizeof(ESL_DSQ *) * (allocM+1));
  ESL_ALLOC(gm_tr->rf,        sizeof(char)      * (allocM+2)); /* yes, +2: each is (0)1..M, +trailing \0  */
  ESL_ALLOC(gm_tr->mm,        sizeof(char)      * (allocM+2));
  ESL_ALLOC(gm_tr->cs,        sizeof(char)      * (allocM+2));
  ESL_ALLOC(gm_tr->consensus, sizeof(char)      * (allocM+2));
  gm_tr->rsc[0]       = NULL;
  gm_tr->codons[0]    = NULL;

  /* level 2 */
  ESL_ALLOC(gm_tr->rsc[0], sizeof(float) * (allocM+1) * (p7P_MAXCODONS1 + abc->Kp));

  for (x = 1; x <= allocM; x++)   
    gm_tr->rsc[x] = gm_tr->rsc[0] + x * (p7P_MAXCODONS1 + abc->Kp);

  ESL_ALLOC(gm_tr->codons[0], sizeof(ESL_DSQ) * (allocM+1) * (p7P_MAXCODONS1+1)); /* +1 for trailing \0 */

  for (x = 1; x <= allocM; x++)
    gm_tr->codons[x] = gm_tr->codons[0] + x * p7P_MAXCODONS1;

  for (x = 1; x <= allocM; x++)

  /* Initialize some edge pieces of memory that are never used,
   * and are only present for indexing convenience. */
  esl_vec_FSet(gm_tr->tsc, p7P_NTRANS, -eslINFINITY);     /* node 0 nonexistent, has no transitions  */
  if (allocM > 1) {
    p7P_TSC(gm_tr, 1, p7P_DM) = -eslINFINITY;             /* delete state D_1 is wing-retracted      */
    p7P_TSC(gm_tr, 1, p7P_DD) = -eslINFINITY;
  }

  for (x = 0; x < (p7P_MAXCODONS1 + abc->Kp); x++) 
    p7P_MSC_CODON(gm_tr, 0,      x) = -eslINFINITY;             /* no emissions from nonexistent M_0... */
  
  /* Set remaining info  */
  gm_tr->mode             = p7_NO_MODE;
  gm_tr->L                = 0;
  gm_tr->allocM           = allocM;
  gm_tr->M                = 0;
  gm_tr->max_length       = -1;
  gm_tr->nj               = 0.0f;
  gm_tr->abc              = abc;

  gm_tr->roff             = -1;
  gm_tr->eoff             = -1;
  gm_tr->offs[p7_MOFFSET] = -1;
  gm_tr->offs[p7_FOFFSET] = -1;
  gm_tr->offs[p7_POFFSET] = -1;
  
  gm_tr->codon_lengths    = 1;
  gm_tr->fs               = 0.;

  gm_tr->name            = NULL;
  gm_tr->acc             = NULL;
  gm_tr->desc            = NULL;
  gm_tr->rf[0]           = 0;     /* RF line is optional annotation; this flags that it's not set yet */
  gm_tr->mm[0]           = 0;     /* likewise for MM annotation line */
  gm_tr->cs[0]           = 0;     /* likewise for CS annotation line */
  gm_tr->consensus[0]    = 0;

  for (x = 0; x < p7_NEVPARAM; x++) gm_tr->evparam[x] = p7_EVPARAM_UNSET;
  for (x = 0; x < p7_NCUTOFFS; x++) gm_tr->cutoff[x]  = p7_CUTOFF_UNSET;
  for (x = 0; x < p7_MAXABET;  x++) gm_tr->compo[x]   = p7_COMPO_UNSET;

  return gm_tr;

 ERROR:
  p7_profile_fs_Destroy(gm_tr);
  return NULL;


}



/* Function:  p7_profile_Copy()
 * Synopsis:  Copy a profile.
 *
 * Purpose:   Copies profile <src> to profile <dst>, where <dst>
 *            has already been allocated to be of sufficient size.
 *
 * Returns:   <eslOK> on success.
 * 
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if <dst> is too small 
 *            to fit <src>.
 */
int
p7_profile_Copy(const P7_PROFILE *src, P7_PROFILE *dst)
{
  int x,z;
  int status;

  if (src->M > dst->allocM) ESL_EXCEPTION(eslEINVAL, "destination profile is too small to hold a copy of source profile");

  esl_vec_FCopy(src->tsc, src->M*p7P_NTRANS, dst->tsc);
  for (x = 0; x < src->abc->Kp;   x++) esl_vec_FCopy(src->rsc[x], (src->M+1)*p7P_NR, dst->rsc[x]);
  for (x = 0; x < p7P_NXSTATES;   x++) esl_vec_FCopy(src->xsc[x], p7P_NXTRANS,       dst->xsc[x]);

  dst->mode        = src->mode;
  dst->L           = src->L;
  dst->allocM      = src->allocM;
  dst->M           = src->M;
  dst->max_length  = src->max_length;
  dst->nj          = src->nj;

  dst->roff        = src->roff;
  dst->eoff        = src->eoff;
  for (x = 0; x < p7_NOFFSETS; ++x) dst->offs[x] = src->offs[x];

  if (dst->name != NULL) free(dst->name);
  if (dst->acc  != NULL) free(dst->acc);
  if (dst->desc != NULL) free(dst->desc);

  if ((status = esl_strdup(src->name,      -1, &(dst->name)))      != eslOK) return status;
  if ((status = esl_strdup(src->acc,       -1, &(dst->acc)))       != eslOK) return status;
  if ((status = esl_strdup(src->desc,      -1, &(dst->desc)))      != eslOK) return status;

  strcpy(dst->rf,        src->rf);         /* RF is optional: if it's not set, *rf=0, and strcpy still works fine */
  strcpy(dst->mm,        src->mm);         /* MM is also optional annotation */
  strcpy(dst->cs,        src->cs);         /* CS is also optional annotation */
  strcpy(dst->consensus, src->consensus);  /* consensus though is always present on a valid profile */

  for (z = 0; z < p7_NEVPARAM; z++) dst->evparam[z] = src->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) dst->cutoff[z]  = src->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) dst->compo[z]   = src->compo[z];
  return eslOK;
}

/* Function:  p7_profile_fs5_Copy()
 * Synopsis:  Copy a 5 cocon length frameshift profile.
 *
 * Purpose:   Copies profile <src> to profile <dst>, where <dst>
 *            has already been allocated to be of sufficient size.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if <dst> is too small
 *            to fit <src>.
 */
int
p7_profile_fs5_Copy(const P7_FS_PROFILE *src, P7_FS_PROFILE *dst)
{
  int x,z;
  int status;

  if (src->M > dst->allocM)       ESL_EXCEPTION(eslEINVAL, "destination profile is too small to hold a copy of source profile");
  if (dst->codon_lengths != 5)    ESL_EXCEPTION(eslEINVAL, "destination proflie not allocated for 5 codon lengths");

  esl_vec_FCopy(src->tsc, src->M*p7P_NTRANS, dst->tsc);
  for (x = 0; x <= src->M;      x++) { esl_vec_FCopy( src->rsc[x],       (p7P_MAXCODONS5 + src->abc->Kp), dst->rsc[x]);       }
  for (x = 0; x < p7P_NXSTATES; x++) { esl_vec_FCopy( src->xsc[x],       p7P_NXTRANS,                    dst->xsc[x]);       }
  for (x = 0; x <= src->M;      x++) { esl_abc_dsqcpy(src->codons[x],    p7P_MAXCODONS5,                  dst->codons[x]);    }
  for (x = 0; x <= src->M;      x++) { esl_abc_dsqcpy(src->indel_pos[x], p7P_MAXCODONS5,                  dst->indel_pos[x]); }

  dst->mode        = src->mode;
  dst->L           = src->L;
  dst->allocM      = src->allocM;
  dst->M           = src->M;
  dst->max_length  = src->max_length;
  dst->nj          = src->nj;

  dst->roff        = src->roff;
  dst->eoff        = src->eoff;
  for (x = 0; x < p7_NOFFSETS; ++x) dst->offs[x] = src->offs[x];

  if (dst->name != NULL) free(dst->name);
  if (dst->acc  != NULL) free(dst->acc);
  if (dst->desc != NULL) free(dst->desc);

  if ((status = esl_strdup(src->name,      -1, &(dst->name)))      != eslOK) return status;
  if ((status = esl_strdup(src->acc,       -1, &(dst->acc)))       != eslOK) return status;
  if ((status = esl_strdup(src->desc,      -1, &(dst->desc)))      != eslOK) return status;

  strcpy(dst->rf,        src->rf);         /* RF is optional: if it's not set, *rf=0, and strcpy still works fine */
  strcpy(dst->mm,        src->mm);         /* MM is also optional annotation */
  strcpy(dst->cs,        src->cs);         /* CS is also optional annotation */
  strcpy(dst->consensus, src->consensus);  /* consensus though is always present on a valid profile */

  for (z = 0; z < p7_NEVPARAM; z++) dst->evparam[z] = src->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) dst->cutoff[z]  = src->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) dst->compo[z]   = src->compo[z];
  return eslOK;
}



/* Function:  p7_profile_fs3_Copy()
 * Synopsis:  Copy a 3 cocon length frameshift profile.
 *
 * Purpose:   Copies profile <src> to profile <dst>, where <dst>
 *            has already been allocated to be of sufficient size.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if <dst> is too small
 *            to fit <src>.
 */
int
p7_profile_fs3_Copy(const P7_FS_PROFILE *src, P7_FS_PROFILE *dst)
{
  int x,z;
  int status;

  if (src->M > dst->allocM)       ESL_EXCEPTION(eslEINVAL, "destination profile is too small to hold a copy of source profile");
  if (dst->codon_lengths != 3)    ESL_EXCEPTION(eslEINVAL, "destination proflie not allocated for 3 codon lengths");

  esl_vec_FCopy(src->tsc, src->M*p7P_NTRANS, dst->tsc);
  for (x = 0; x <= src->M;      x++) { esl_vec_FCopy( src->rsc[x],       (p7P_MAXCODONS3 + src->abc->Kp), dst->rsc[x]);       }
  for (x = 0; x < p7P_NXSTATES; x++) { esl_vec_FCopy( src->xsc[x],       p7P_NXTRANS,                    dst->xsc[x]);       }
  for (x = 0; x <= src->M;      x++) { esl_abc_dsqcpy(src->codons[x],    p7P_MAXCODONS3,                  dst->codons[x]);    }
  for (x = 0; x <= src->M;      x++) { esl_abc_dsqcpy(src->indel_pos[x], p7P_MAXCODONS3,                  dst->indel_pos[x]); }

  dst->mode        = src->mode;
  dst->L           = src->L;
  dst->allocM      = src->allocM;
  dst->M           = src->M;
  dst->max_length  = src->max_length;
  dst->nj          = src->nj;
  dst->fs          = src->fs;

  dst->roff        = src->roff;
  dst->eoff        = src->eoff;
  for (x = 0; x < p7_NOFFSETS; ++x) dst->offs[x] = src->offs[x];

  if (dst->name != NULL) free(dst->name);
  if (dst->acc  != NULL) free(dst->acc);
  if (dst->desc != NULL) free(dst->desc);

  if ((status = esl_strdup(src->name,      -1, &(dst->name)))      != eslOK) return status;
  if ((status = esl_strdup(src->acc,       -1, &(dst->acc)))       != eslOK) return status;
  if ((status = esl_strdup(src->desc,      -1, &(dst->desc)))      != eslOK) return status;

  strcpy(dst->rf,        src->rf);         /* RF is optional: if it's not set, *rf=0, and strcpy still works fine */
  strcpy(dst->mm,        src->mm);         /* MM is also optional annotation */
  strcpy(dst->cs,        src->cs);         /* CS is also optional annotation */
  strcpy(dst->consensus, src->consensus);  /* consensus though is always present on a valid profile */

  for (z = 0; z < p7_NEVPARAM; z++) dst->evparam[z] = src->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) dst->cutoff[z]  = src->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) dst->compo[z]   = src->compo[z];
  return eslOK;
}





/* Function:  p7_profile_tr_Copy()
 * Synopsis:  Copy a translated profile.
 *
 * Purpose:   Copies profile <src> to profile <dst>, where <dst>
 *            has already been allocated to be of sufficient size.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error; <eslEINVAL> if <dst> is too small
 *            to fit <src>.
 */
int
p7_profile_tr_Copy(const P7_FS_PROFILE *src, P7_FS_PROFILE *dst)
{
  int x,z;
  int status;

  if (src->M > dst->allocM)       ESL_EXCEPTION(eslEINVAL, "destination profile is too small to hold a copy of source profile");
  if (dst->codon_lengths != 1)    ESL_EXCEPTION(eslEINVAL, "destination proflie not allocated for 1 codon length");

  esl_vec_FCopy(src->tsc, src->M*p7P_NTRANS, dst->tsc);
  for (x = 0; x <= src->M;      x++) { esl_vec_FCopy( src->rsc[x],       (p7P_MAXCODONS1 + src->abc->Kp), dst->rsc[x]);       }
  for (x = 0; x < p7P_NXSTATES; x++) { esl_vec_FCopy( src->xsc[x],       p7P_NXTRANS,                    dst->xsc[x]);       }
  for (x = 0; x <= src->M;      x++) { esl_abc_dsqcpy(src->codons[x],    p7P_MAXCODONS1,                  dst->codons[x]);    }

  dst->mode        = src->mode;
  dst->L           = src->L;
  dst->allocM      = src->allocM;
  dst->M           = src->M;
  dst->max_length  = src->max_length;
  dst->nj          = src->nj;
  dst->fs          = src->fs;

  dst->roff        = src->roff;
  dst->eoff        = src->eoff;
  for (x = 0; x < p7_NOFFSETS; ++x) dst->offs[x] = src->offs[x];

  if (dst->name != NULL) free(dst->name);
  if (dst->acc  != NULL) free(dst->acc);
  if (dst->desc != NULL) free(dst->desc);

  if ((status = esl_strdup(src->name,      -1, &(dst->name)))      != eslOK) return status;
  if ((status = esl_strdup(src->acc,       -1, &(dst->acc)))       != eslOK) return status;
  if ((status = esl_strdup(src->desc,      -1, &(dst->desc)))      != eslOK) return status;

  strcpy(dst->rf,        src->rf);         /* RF is optional: if it's not set, *rf=0, and strcpy still works fine */
  strcpy(dst->mm,        src->mm);         /* MM is also optional annotation */
  strcpy(dst->cs,        src->cs);         /* CS is also optional annotation */
  strcpy(dst->consensus, src->consensus);  /* consensus though is always present on a valid profile */

  for (z = 0; z < p7_NEVPARAM; z++) dst->evparam[z] = src->evparam[z];
  for (z = 0; z < p7_NCUTOFFS; z++) dst->cutoff[z]  = src->cutoff[z];
  for (z = 0; z < p7_MAXABET;  z++) dst->compo[z]   = src->compo[z];
  return eslOK;
}





/* Function:  p7_profile_Clone()
 * Synopsis:  Duplicates a profile.
 *
 * Purpose:   Duplicate profile <gm>; return a pointer
 *            to the newly allocated copy.
 */
P7_PROFILE *
p7_profile_Clone(const P7_PROFILE *gm)
{
  P7_PROFILE *g2 = NULL;
  int         status;

  if ((g2 = p7_profile_Create(gm->allocM, gm->abc)) == NULL) return NULL;
  if ((status = p7_profile_Copy(gm, g2)) != eslOK) goto ERROR;
  return g2;
  
 ERROR:
  p7_profile_Destroy(g2);
  return NULL;
}

/* Function:  p7_profile_fs_Clone()
 * Synopsis:  Duplicates a profile.
 *
 * Purpose:   Duplicate profile <gm>; return a pointer
 *            to the newly allocated copy.
 */
P7_FS_PROFILE *
p7_profile_fs_Clone(const P7_FS_PROFILE *gm_fs)
{
  P7_FS_PROFILE *g2 = NULL;
  int         status;

  if(gm_fs->codon_lengths == 5) {
    if ((g2 = p7_profile_fs5_Create(gm_fs->allocM, gm_fs->abc)) == NULL) return NULL;
    if ((status = p7_profile_fs5_Copy(gm_fs, g2)) != eslOK) goto ERROR;
  }
  else if(gm_fs->codon_lengths == 3) {
    if ((g2 = p7_profile_fs3_Create(gm_fs->allocM, gm_fs->abc)) == NULL) return NULL;
    if ((status = p7_profile_fs3_Copy(gm_fs, g2)) != eslOK) goto ERROR;
  }
  else if(gm_fs->codon_lengths == 1) {
    if ((g2 = p7_profile_tr_Create(gm_fs->allocM, gm_fs->abc)) == NULL) return NULL;
    if ((status = p7_profile_tr_Copy(gm_fs, g2)) != eslOK) goto ERROR;
  }
  else ESL_XEXCEPTION(eslEINVAL, "invalid number of codon lengths %d (1,3 and 5 acceppted)", gm_fs->codon_lengths);
 
  return g2;

 ERROR:
  p7_profile_fs_Destroy(g2);
  return NULL;
}

/* Function:  p7_profile_GetFwdEmissionArray()
 * Synopsis:  Retrieve Fwd (float) residue emission values from an optimized
 *            profile into an array
 *
 * Purpose:   Extract an implicitly 2D array of 32-bit float Fwd residue
 *            emission values from an optimized profile <om>, converting
 *            back to emission values based on the background. <arr> must
 *            be allocated by the calling function to be of size
 *            ( om->abc->Kp * ( om->M  + 1 )), and indexing into the array
 *            is done as  [om->abc->Kp * i +  c ] for character c at
 *            position i.
 *
 * Args:      <om>   - optimized profile, containing transition information
 *            <bg>   - background frequencies
 *            <arr>  - preallocated array into which scores will be placed
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
int
p7_profile_GetFwdEmissionArray(const P7_PROFILE *gm, P7_BG *bg, float *arr )
{
  int i, j;

  for (i = 1; i <= gm->M; i++) {
    for (j=0; j<gm->abc->Kp; j++) {
      arr[i*gm->abc->Kp + j] =  bg->f[j] * exp( gm->rsc[j][(i) * p7P_NR     + p7P_MSC]);
    }
  }

  return eslOK;
}


/* Function:  p7_profile_SetNullEmissions()
 * Synopsis:  Set all emission scores to zero (experimental).
 *
 * Purpose:   Set all emission scores in profile <gm> to zero.
 *            This makes the profile a null model, with all the same
 *            length distributions as the original model, but
 *            the emission probabilities of the background.
 *            
 *            Written to test the idea that score statistics will be
 *            even better behaved when using a null model with the
 *            same length distribution as the search model.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_profile_SetNullEmissions(P7_PROFILE *gm)
{
  int x;
  for (x = 0; x <= gm->abc->K; x++)                esl_vec_FSet(gm->rsc[x], (gm->M+1)*p7P_NR, 0.0);   /* canonicals    */
  for (x = gm->abc->K+1; x <= gm->abc->Kp-3; x++)  esl_vec_FSet(gm->rsc[x], (gm->M+1)*p7P_NR, 0.0);   /* noncanonicals */
  return eslOK;
}


/* Function:  p7_profile_Reuse()
 * Synopsis:  Prepare profile to be re-used for a new HMM.
 *
 * Purpose:   Prepare profile <gm>'s memory to be re-used
 *            for a new HMM.
 */
int
p7_profile_Reuse(P7_PROFILE *gm)
{
  /* name, acc, desc annotation is dynamically allocated for each HMM */
  if (gm->name != NULL) { free(gm->name); gm->name = NULL; }
  if (gm->acc  != NULL) { free(gm->acc);  gm->acc  = NULL; }
  if (gm->desc != NULL) { free(gm->desc); gm->desc = NULL; }

  /* set annotations to empty strings */
  gm->rf[0]        = 0;
  gm->mm[0]        = 0;
  gm->cs[0]        = 0;
  gm->consensus[0] = 0;
      
  /* reset some other things, but leave the rest alone. */
  gm->mode = p7_NO_MODE;
  gm->L    = 0;
  gm->M    = 0;
  gm->nj   = 0.0f;

  gm->roff             = -1;
  gm->eoff             = -1;
  gm->offs[p7_MOFFSET] = -1;
  gm->offs[p7_FOFFSET] = -1;
  gm->offs[p7_POFFSET] = -1;

  return eslOK;
}


/* Function:  p7_profile_Sizeof()
 * Synopsis:  Return the allocated size of a P7_PROFILE.
 *
 * Purpose:   Return the allocated size of a <P7_PROFILE>, in bytes.
 */
size_t
p7_profile_Sizeof(P7_PROFILE *gm)
{
  size_t n = 0;

  /* these mirror malloc()'s in p7_profile_Create(); maintain one:one correspondence for maintainability */
  n += sizeof(P7_PROFILE);
  n += sizeof(float)   * gm->allocM * p7P_NTRANS;             /* gm->tsc       */
  n += sizeof(float *) * gm->abc->Kp;                        /* gm->rsc       */
  n += sizeof(char)    * (gm->allocM+2);                /* gm->rf        */
  n += sizeof(char)    * (gm->allocM+2);                /* gm->mm        */
  n += sizeof(char)    * (gm->allocM+2);                /* gm->cs        */
  n += sizeof(char)    * (gm->allocM+2);                /* gm->consensus */

  n += sizeof(float) * gm->abc->Kp * (gm->allocM+1) * p7P_NR; /* gm->rsc[0]    */

  return n;
}


/* Function:  p7_profile_Destroy()
 * Synopsis:  Frees a profile.
 *
 * Purpose:   Frees a profile <gm>.
 *
 * Returns:   (void).
 *
 * Xref:      STL11/125.
 */
void
p7_profile_Destroy(P7_PROFILE *gm)
{
  if (gm != NULL) {
    if (gm->rsc   != NULL && gm->rsc[0] != NULL) free(gm->rsc[0]);
    if (gm->tsc       != NULL) free(gm->tsc);
    if (gm->rsc       != NULL) free(gm->rsc);
    if (gm->name      != NULL) free(gm->name);
    if (gm->acc       != NULL) free(gm->acc);
    if (gm->desc      != NULL) free(gm->desc);
    if (gm->rf        != NULL) free(gm->rf);
    if (gm->mm        != NULL) free(gm->mm);
    if (gm->cs        != NULL) free(gm->cs);
    if (gm->consensus != NULL) free(gm->consensus);
    free(gm);
  }
  return;
}

/* Function:  p7_profile_fs_Destroy()
 * Synopsis:  Frees a profile.
 *
 * Purpose:   Frees a profile <gm>.
 *
 * Returns:   (void).
 *
 * Xref:      STL11/125.
 */
void
p7_profile_fs_Destroy(P7_FS_PROFILE *gm)
{
  if (gm != NULL) {
    if (gm->rsc       != NULL && gm->rsc[0] != NULL) free(gm->rsc[0]);
    if (gm->codons    != NULL && gm->codons[0] != NULL) free(gm->codons[0]);
    if (gm->indel_pos != NULL && gm->indel_pos[0] != NULL) free(gm->indel_pos[0]);
    if (gm->tsc       != NULL) free(gm->tsc);
    if (gm->rsc       != NULL) free(gm->rsc);
    if (gm->codons    != NULL) free(gm->codons);
    if (gm->indel_pos != NULL) free(gm->indel_pos);
    if (gm->name      != NULL) free(gm->name);
    if (gm->acc       != NULL) free(gm->acc);
    if (gm->desc      != NULL) free(gm->desc);
    if (gm->rf        != NULL) free(gm->rf);
    if (gm->mm        != NULL) free(gm->mm);
    if (gm->cs        != NULL) free(gm->cs);
    if (gm->consensus != NULL) free(gm->consensus);
    free(gm);
  }
  return;
}


/*****************************************************************
 * 2. Access methods.
 *****************************************************************/

/* Function:  p7_profile_IsLocal()
 * Synopsis:  Return TRUE if profile is in a local alignment mode.
 *
 * Purpose:   Return <TRUE> if profile is in a local alignment mode.
 */
int
p7_profile_IsLocal(const P7_PROFILE *gm)
{
  if (gm->mode == p7_UNILOCAL || gm->mode == p7_LOCAL) return TRUE;
  return FALSE;
}

/* Function:  p7_fs_profile_IsLocal()
 * Synopsis:  Return TRUE if profile is in a local alignment mode.
 *
 * Purpose:   Return <TRUE> if profile is in a local alignment mode.
 */
int
p7_fs_profile_IsLocal(const P7_FS_PROFILE *gm)
{
  if (gm->mode == p7_UNILOCAL || gm->mode == p7_LOCAL) return TRUE;
  return FALSE;
}

/* Function:  p7_profile_IsMultihit()
 * Synopsis:  Return TRUE if profile is in a multihit alignment mode.
 *
 * Purpose:   Return <TRUE> if profile is in a multihit alignment mode.
 */
int
p7_profile_IsMultihit(const P7_PROFILE *gm)
{
  if (gm->mode == p7_LOCAL || gm->mode == p7_GLOCAL) return TRUE;
  return FALSE;
}

/* Function:  p7_fs_profile_IsMultihit()
 * Synopsis:  Return TRUE if profile is in a multihit alignment mode.
 *
 * Purpose:   Return <TRUE> if profile is in a multihit alignment mode.
 */
int
p7_fs_profile_IsMultihit(const P7_FS_PROFILE *gm)
{
  if (gm->mode == p7_LOCAL || gm->mode == p7_GLOCAL ) return TRUE;
   return FALSE;
}


/* Function:  p7_profile_GetT()
 *
 * Purpose:   Convenience function that looks up a transition score in
 *            profile <gm> for a transition from state type <st1> in
 *            node <k1> to state type <st2> in node <k2>. For unique
 *            state types that aren't in nodes (<p7T_S>, for example), the
 *            <k> value is ignored, though it would be customarily passed as 0.
 *            Return the transition score in <ret_tsc>.
 *            
 *            This function would almost always be called on profile
 *            traces, of course, but it's possible to call it
 *            on core traces (for example, if you were to try to 
 *            trace_Dump() during HMM construction, and you wanted
 *            to see detailed profile scores for that trace). Core traces
 *            can contain <p7T_X> "states" used solely to signal
 *            a sequence fragment, treated as missing data. Transitions
 *            involving <p7T_X> states are assigned zero score here.
 *            Other transitions that occur only in core traces
 *            (B->I0, B->D1, I_M->E) also silently get a zero score.
 *            This is safe, because we would only ever use this number
 *            for display, not as a log probability somewhere.
 *
 * Returns:   <eslOK> on success, and <*ret_tsc> contains the requested
 *            transition score.            
 * 
 * Throws:    <eslEINVAL> if a nonexistent transition is requested. Now
 *            <*ret_tsc> is set to $-\infty$.
 *            
 */
int
p7_profile_GetT(const P7_PROFILE *gm, char st1, int k1, char st2, int k2, float *ret_tsc)
{
  float tsc = 0.0f;
  int   status;

  /* Detect transitions that can only come from core traces;
   * return 0.0 as a special case (this is only done for displaying
   * "scores" in trace dumps, during debugging.)
   */
  if (st1 == p7T_X || st2 == p7T_X) return eslOK;
  if (st1 == p7T_B && st2 == p7T_I) return eslOK;
  if (st1 == p7T_B && st2 == p7T_D) return eslOK;
  if (st1 == p7T_I && st2 == p7T_E) return eslOK;

  /* Now we're sure this is a profile trace, as it should usually be. */
  switch (st1) {
  case p7T_S:  break;
  case p7T_T:  break;
  case p7T_N:
    switch (st2) {
    case p7T_B: tsc =  gm->xsc[p7P_N][p7P_MOVE]; break;
    case p7T_N: tsc =  gm->xsc[p7P_N][p7P_LOOP]; break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DecodeStatetype(st1), p7_hmm_DecodeStatetype(st2));
    }
    break;

  case p7T_B:
    switch (st2) {
    case p7T_M: tsc = p7P_TSC(gm, k2-1, p7P_BM); break; /* remember, B->Mk is stored in [k-1][p7P_BM] */
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DecodeStatetype(st1), p7_hmm_DecodeStatetype(st2));
    }
    break;

  case p7T_M:
    switch (st2) {
    case p7T_M: tsc = p7P_TSC(gm, k1, p7P_MM); break;
    case p7T_I: tsc = p7P_TSC(gm, k1, p7P_MI); break;
    case p7T_D: tsc = p7P_TSC(gm, k1, p7P_MD); break;
    case p7T_E: 
      if (k1 != gm->M && ! p7_profile_IsLocal(gm)) ESL_EXCEPTION(eslEINVAL, "local end transition (M%d of %d) in non-local model", k1, gm->M);
      tsc = 0.0f;    /* by def'n in H3 local alignment */
      break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", p7_hmm_DecodeStatetype(st1), k1, p7_hmm_DecodeStatetype(st2));
    }
    break;

  case p7T_D:
    switch (st2) {
    case p7T_M: tsc = p7P_TSC(gm, k1, p7P_DM); break;
    case p7T_D: tsc = p7P_TSC(gm, k1, p7P_DD); break;
    case p7T_E: 
      if (k1 != gm->M && ! p7_profile_IsLocal(gm)) ESL_EXCEPTION(eslEINVAL, "local end transition (D%d of %d) in non-local model", k1, gm->M);
      tsc = 0.0f;    /* by def'n in H3 local alignment */
      break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", p7_hmm_DecodeStatetype(st1), k1, p7_hmm_DecodeStatetype(st2));
    }
    break;

  case p7T_I:
    switch (st2) {
    case p7T_M: tsc = p7P_TSC(gm, k1, p7P_IM); break;
    case p7T_I: tsc = p7P_TSC(gm, k1, p7P_II); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", p7_hmm_DecodeStatetype(st1), k1, p7_hmm_DecodeStatetype(st2));
    }
    break;

  case p7T_E:
    switch (st2) {
    case p7T_C: tsc = gm->xsc[p7P_E][p7P_MOVE]; break;
    case p7T_J: tsc = gm->xsc[p7P_E][p7P_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DecodeStatetype(st1), p7_hmm_DecodeStatetype(st2));
    }
    break;

  case p7T_J:
    switch (st2) {
    case p7T_B: tsc = gm->xsc[p7P_J][p7P_MOVE]; break;
    case p7T_J: tsc = gm->xsc[p7P_J][p7P_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DecodeStatetype(st1), p7_hmm_DecodeStatetype(st2));
    }
    break;

  case p7T_C:
    switch (st2) {
    case p7T_T:  tsc = gm->xsc[p7P_C][p7P_MOVE]; break;
    case p7T_C:  tsc = gm->xsc[p7P_C][p7P_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DecodeStatetype(st1), p7_hmm_DecodeStatetype(st2));
    }
    break;

  default: ESL_XEXCEPTION(eslEINVAL, "bad state type %d in traceback", st1);
  }

  *ret_tsc = tsc;
  return eslOK;

 ERROR:
  *ret_tsc = -eslINFINITY;
  return status;
}

int
p7_profile_fs_GetT(const P7_FS_PROFILE *gm_fs, char st1, int k1, char st2, int k2, float *ret_tsc)
{
  float tsc = 0.0f;
  int   status;

  /* Detect transitions that can only come from core traces;
   * return 0.0 as a special case (this is only done for displaying
   * "scores" in trace dumps, during debugging.)
   */
  if (st1 == p7T_X || st2 == p7T_X) return eslOK;
  if (st1 == p7T_B && st2 == p7T_I) return eslOK;
  if (st1 == p7T_B && st2 == p7T_D) return eslOK;
  if (st1 == p7T_I && st2 == p7T_E) return eslOK;
  if (st1 == p7T_M && st2 == p7T_R) return eslOK;
  if (st1 == p7T_D && st2 == p7T_R) return eslOK;
  if (st1 == p7T_R && st2 == p7T_P) return eslOK;
  if (st1 == p7T_P && st2 == p7T_P) return eslOK;
  if (st1 == p7T_P && st2 == p7T_M) return eslOK;
  
  /* Now we're sure this is a profile trace, as it should usually be. */
  switch (st1) {
  case p7T_S:  break;
  case p7T_T:  break;
  case p7T_N:
    switch (st2) {
    case p7T_B: tsc =  gm_fs->xsc[p7P_N][p7P_MOVE]; break;
    case p7T_N: tsc =  gm_fs->xsc[p7P_N][p7P_LOOP]; break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DecodeStatetype(st1), p7_hmm_DecodeStatetype(st2));
    }
    break;

  case p7T_B:
    switch (st2) {
    case p7T_M: tsc = p7P_TSC(gm_fs, k2-1, p7P_BM); break; /* remember, B->Mk is stored in [k-1][p7P_BM] */
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DecodeStatetype(st1), p7_hmm_DecodeStatetype(st2));
    }
    break;

  case p7T_M:
    switch (st2) {
    case p7T_M: tsc = p7P_TSC(gm_fs, k1, p7P_MM); break;
    case p7T_I: tsc = p7P_TSC(gm_fs, k1, p7P_MI); break;
    case p7T_D: tsc = p7P_TSC(gm_fs, k1, p7P_MD); break;
    case p7T_E: 
      if (k1 != gm_fs->M && ! p7_fs_profile_IsLocal(gm_fs)) ESL_EXCEPTION(eslEINVAL, "local end transition (M%d of %d) in non-local model", k1, gm_fs->M);
      tsc = 0.0f;    /* by def'n in H3 local alignment */
      break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", p7_hmm_DecodeStatetype(st1), k1, p7_hmm_DecodeStatetype(st2));
    }
    break;

  case p7T_D:
    switch (st2) {
    case p7T_M: tsc = p7P_TSC(gm_fs, k1, p7P_DM); break;
    case p7T_D: tsc = p7P_TSC(gm_fs, k1, p7P_DD); break;
    case p7T_E: 
      if (k1 != gm_fs->M && ! p7_fs_profile_IsLocal(gm_fs)) ESL_EXCEPTION(eslEINVAL, "local end transition (D%d of %d) in non-local model", k1, gm_fs->M);
      tsc = 0.0f;    /* by def'n in H3 local alignment */
      break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", p7_hmm_DecodeStatetype(st1), k1, p7_hmm_DecodeStatetype(st2));
    }
    break;

  case p7T_I:
    switch (st2) {
    case p7T_M: tsc = p7P_TSC(gm_fs, k1, p7P_IM); break;
    case p7T_I: tsc = p7P_TSC(gm_fs, k1, p7P_II); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", p7_hmm_DecodeStatetype(st1), k1, p7_hmm_DecodeStatetype(st2));
    }
    break;

  case p7T_E:
    switch (st2) {
    case p7T_C: tsc = gm_fs->xsc[p7P_E][p7P_MOVE]; break;
    case p7T_J: tsc = gm_fs->xsc[p7P_E][p7P_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DecodeStatetype(st1), p7_hmm_DecodeStatetype(st2));
    }
    break;

  case p7T_J:
    switch (st2) {
    case p7T_B: tsc = gm_fs->xsc[p7P_J][p7P_MOVE]; break;
    case p7T_J: tsc = gm_fs->xsc[p7P_J][p7P_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DecodeStatetype(st1), p7_hmm_DecodeStatetype(st2));
    }
    break;

  case p7T_C:
    switch (st2) {
    case p7T_T:  tsc = gm_fs->xsc[p7P_C][p7P_MOVE]; break;
    case p7T_C:  tsc = gm_fs->xsc[p7P_C][p7P_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", p7_hmm_DecodeStatetype(st1), p7_hmm_DecodeStatetype(st2));
    }
    break;

  default: ESL_XEXCEPTION(eslEINVAL, "bad state type %d in traceback", st1);
  }

  *ret_tsc = tsc;
  return eslOK;

 ERROR:
  *ret_tsc = -eslINFINITY;
  return status;
}



/*****************************************************************
 * 3. Debugging and development code.
 *****************************************************************/

/* Function:  p7_profile_Validate()
 *
 * Purpose:   Validates the internals of the generic profile structure
 *            <gm>.
 *            
 *            TODO: currently this function is incomplete, and only
 *            validates the entry distribution.
 *            
 * Returns:   <eslOK> if <gm> internals look fine. Returns <eslFAIL>
 *            if something is wrong, and leaves an error message in
 *            <errbuf> if caller passed it non-<NULL>.
 */
int
p7_profile_Validate(const P7_PROFILE *gm, char *errbuf, float tol)
{
  int     status;
  int     k;
  double *pstart = NULL;

  ESL_ALLOC(pstart, sizeof(double) * (gm->M+1));
  pstart[0] = 0.0;

  /* Validate the entry distribution.
   * In a glocal model, this is an explicit probability distribution,
   * corresponding to left wing retraction.
   * In a local model, this is an implicit probability distribution,
   * corresponding to the implicit local alignment model, and we have
   * to calculate the M(M+1)/2 fragment probabilities accordingly.
   */
  if (p7_profile_IsLocal(gm))
    {        /* the code block below is also in emit.c:sample_endpoints */
      for (k = 1; k <= gm->M; k++)
  pstart[k] = exp(p7P_TSC(gm, k-1, p7P_BM)) * (gm->M - k + 1); /* multiply p_ij by the number of exits j */
    }
  else
    {
      for (k = 1; k <= gm->M; k++)
  pstart[k] = exp(p7P_TSC(gm, k-1, p7P_BM));
    }

  if (esl_vec_DValidate(pstart, gm->M+1, tol, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "profile entry distribution is not normalized properly");
  free(pstart);
  return eslOK;

 ERROR:
  if (pstart != NULL) free(pstart);
  return eslFAIL;
}

int
p7_profile_fs_Validate(const P7_FS_PROFILE *gm_fs, char *errbuf, float tol)
{
  int     status;
  int     k;
  double *pstart = NULL;

  ESL_ALLOC(pstart, sizeof(double) * (gm_fs->M+1));
  pstart[0] = 0.0;

  /* Validate the entry distribution.
   * In a glocal model, this is an explicit probability distribution,
   * corresponding to left wing retraction.
   * In a local model, this is an implicit probability distribution,
   * corresponding to the implicit local alignment model, and we have
   * to calculate the M(M+1)/2 fragment probabilities accordingly.
   */
  if (p7_fs_profile_IsLocal(gm_fs))
    {        /* the code block below is also in emit.c:sample_endpoints */
      for (k = 1; k <= gm_fs->M; k++)
  pstart[k] = exp(p7P_TSC(gm_fs, k-1, p7P_BM)) * (gm_fs->M - k + 1); /* multiply p_ij by the number of exits j */
    }
  else
    {
      for (k = 1; k <= gm_fs->M; k++)
  pstart[k] = exp(p7P_TSC(gm_fs, k-1, p7P_BM));
    }

  if (esl_vec_DValidate(pstart, gm_fs->M+1, tol, NULL) != eslOK) ESL_XFAIL(eslFAIL, errbuf, "profile entry distribution is not normalized properly");
  free(pstart);
  return eslOK;

 ERROR:
  if (pstart != NULL) free(pstart);
  return eslFAIL;
}




/* Function:  p7_profile_Compare()
 * Synopsis:  Compare two profiles for equality.
 *
 * Purpose:   Compare two profiles <gm1> and <gm2> to each other.
 *            Return <eslOK> if they're identical, and <eslFAIL> if
 *            they differ. Floating-point probabilities are 
 *            compared for equality within a fractional tolerance
 *            <tol>.  Only compares the scores, not any annotation
 *            on the profiles.
 */
int
p7_profile_Compare(P7_PROFILE *gm1, P7_PROFILE *gm2, float tol)
{
  int x;

  if (gm1->mode != gm2->mode) return eslFAIL;
  if (gm1->M    != gm2->M)    return eslFAIL;

  if (esl_vec_FCompare(gm1->tsc, gm2->tsc, gm1->M*p7P_NTRANS, tol)         != eslOK) return eslFAIL;
  for (x = 0; x < gm1->abc->Kp; x++) 
    if (esl_vec_FCompare(gm1->rsc[x], gm2->rsc[x], (gm1->M+1)*p7P_NR, tol) != eslOK) return eslFAIL;

  for (x = 0; x < p7P_NXSTATES; x++)
    if (esl_vec_FCompare(gm1->xsc[x], gm2->xsc[x], p7P_NXTRANS, tol)       != eslOK) return eslFAIL;

  return eslOK;
}



/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7PROFILE_TESTDRIVE
#include "esl_alphabet.h"
#include "esl_random.h"

static void
utest_Compare(void)
{
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(42);
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_HMM         *hmm  = NULL;
  P7_BG          *bg   = NULL;
  P7_PROFILE     *gm   = NULL;
  P7_PROFILE     *gm2  = NULL;
  int             M    = 200;
  int             L    = 400;

  p7_hmm_Sample(r, M, abc, &hmm); /* master and worker's sampled profiles are identical */
  bg  = p7_bg_Create(abc);
  gm  = p7_profile_Create(hmm->M, abc);
  gm2 = p7_profile_Create(hmm->M, abc);
  p7_ProfileConfig(hmm, bg, gm,  400, p7_LOCAL);
  p7_ProfileConfig(hmm, bg, gm2, 400, p7_LOCAL);
  p7_ReconfigLength(gm,  L);
  p7_ReconfigLength(gm2, L);

  if (p7_profile_Compare(gm, gm2, 0.001) != eslOK) p7_Die("identical profile comparison failed");
  
  p7_profile_Destroy(gm);
  p7_profile_Destroy(gm2);
  p7_bg_Destroy(bg);
  p7_hmm_Destroy(hmm);
  esl_alphabet_Destroy(abc);
  esl_randomness_Destroy(r);
  return;
}


#endif /*p7PROFILE_TESTDRIVE*/

/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7PROFILE_TESTDRIVE

/* gcc -o profile_utest -g -Wall -I../easel -L../easel -I. -L. -Dp7PROFILE_TESTDRIVE p7_profile.c -lhmmer -leasel -lm
 * ./profile_utest
 */
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range toggles reqs incomp  help                                       docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage",              0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for p7_profile.c";

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);

  utest_Compare();

  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7PROFILE_TESTDRIVE*/


