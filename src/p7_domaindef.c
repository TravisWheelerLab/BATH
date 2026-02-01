/* Definition of multidomain structure of a target sequence, and
 * Contents:
 *    1. The P7_DOMAINDEF object: allocation, reuse, destruction
 *    2. Routines inferring domain structure of a target sequence
 *    3. Internal routines 
 *    
 */
#include "p7_config.h"

#include <math.h>
#include <string.h>

#include "easel.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_vectorops.h"
#include "esl_exponential.h"

#include "hmmer.h"

static int is_multidomain_region  (P7_DOMAINDEF *ddef, int i, int j);
static int is_multidomain_region_fs  (P7_DOMAINDEF *ddef, int i, int j);
static int region_trace_ensemble  (P7_DOMAINDEF *ddef, const P7_OPROFILE *om, const ESL_DSQ *dsq, int ireg, int jreg, const P7_OMX *fwd, P7_OMX *wrk, int *ret_nc);
static int region_trace_ensemble_frameshift  (P7_DOMAINDEF *ddef, const P7_FS_PROFILE *gm, const ESL_DSQ *dsq, const ESL_ALPHABET *abc, int ireg, int jreg, const P7_GMX *fwd, P7_GMX *wrk, int *ret_nc);
static int rescore_isolated_domain_frameshift(P7_DOMAINDEF *ddef, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, ESL_SQ *windowsq,  
           P7_GMX *gx1, P7_GMX *gx2, P7_IVX *iv, int i, int j, int null2_is_done, P7_BG *bg, 
           ESL_GENCODE *gcode, int do_biasfilter);
static int rescore_isolated_domain_bath(P7_DOMAINDEF *ddef, P7_OPROFILE *om, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs,
	   const ESL_SQ *orfsq, const ESL_SQ *windowsq, const int64_t ntsqlen, const ESL_GENCODE *gcode, 
           P7_OMX *ox1, P7_OMX *ox2, int i, int j, int null2_is_done, P7_BG *bg);



/*****************************************************************
 * 1. The P7_DOMAINDEF object: allocation, reuse, destruction
 *****************************************************************/

/* Function:  p7_domaindef_Create_BATH()
 * Synopsis:  Creates a new <P7_DOMAINDEF> object.
 *
 * Purpose:   Creates a new <P7_DOMAINDEF> object, with <r> registered
 *            as its random number generator, using default settings
 *            for all thresholds.
 *
 * Returns:   a pointer to the new <P7_DOMAINDEF> object.
 *
 * Throws:    <NULL> on allocation failure.
 */
P7_DOMAINDEF *
p7_domaindef_Create_BATH(ESL_RANDOMNESS *r, ESL_GETOPTS *go)
{
  P7_DOMAINDEF *ddef   = NULL;
  int           Lalloc = 512;  /* this initial alloc doesn't matter much; space is realloced as needed */
  int           nalloc = 32;
  int           status;

  /* level 1 alloc */
  ESL_ALLOC(ddef, sizeof(P7_DOMAINDEF));
  ddef->mocc = ddef->btot = ddef->etot = NULL;
  ddef->n2sc = NULL;
  ddef->sp   = NULL;
  ddef->tr   = NULL;
  ddef->dcl  = NULL;
  
  /* level 2 alloc: posterior prob arrays */
  ESL_ALLOC(ddef->mocc, sizeof(float) * (Lalloc+1));
  ESL_ALLOC(ddef->btot, sizeof(float) * (Lalloc+1));
  ESL_ALLOC(ddef->etot, sizeof(float) * (Lalloc+1));
  ESL_ALLOC(ddef->n2sc, sizeof(float) * (Lalloc+1));
  ddef->mocc[0] = ddef->etot[0] = ddef->btot[0] = 0.;
  ddef->n2sc[0] = 0.;
  ddef->Lalloc  = Lalloc;
  ddef->L       = 0;
  /* level 2 alloc: results storage */
  ESL_ALLOC(ddef->dcl, sizeof(P7_DOMAIN) * nalloc);
  ddef->nalloc = nalloc;
  ddef->ndom   = 0;

  ddef->nexpected  = 0.0;
  ddef->nregions   = 0;
  ddef->nclustered = 0;
  ddef->noverlaps  = 0;
  ddef->nenvelopes = 0;

  /* default thresholds */
  ddef->rt1           = 0.25;
  ddef->rt2           = 0.10;
  ddef->rt3           = 0.20;
  ddef->nsamples      = 200;
  ddef->min_overlap   = 0.8;
  ddef->of_smaller    = TRUE;
  ddef->max_diagdiff  = 4;
  ddef->min_posterior = 0.25;
  ddef->min_endpointp = 0.02;

  /* allocate reusable, growable objects that domain def reuses for each seq */
  ddef->sp  = p7_spensemble_Create(1024, 64, 32); /* init allocs = # sampled pairs; max endpoint range; # of domains */
  ddef->tr  = p7_trace_fs_CreateWithPP();
  ddef->gtr = p7_trace_fs_Create();

  /* keep a copy of ptr to the RNG */
  ddef->r            = r;  
  ddef->do_reseeding = TRUE;

  ddef->fstbl = esl_opt_IsUsed(go, "--fstblout"); /* TRUE to produce tabular frameshift location output */

  return ddef;
  
 ERROR:
  p7_domaindef_Destroy_BATH(ddef);
  return NULL;
}

/* p7_domaindef_GrowTo()
 * Synopsis:  Reallocates a <P7_DOMAINDEF> for new seq length <L>
 * Incept:    SRE, Fri Jan 25 13:27:24 2008 [Janelia]
 *
 * Purpose:   Reallocates a <P7_DOMAINDEF> object <ddef> so that
 *            it can hold a sequence of up to <L> residues. 
 *
 *            (This might be a no-op, if <ddef> is already large
 *            enough.)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure. In this case, the
 *            data in <ddef> are unaffected.
 */
int
p7_domaindef_GrowTo(P7_DOMAINDEF *ddef, int L)
{
  void *p;
  int   status;

  if (L <= ddef->Lalloc) return eslOK;

  ESL_RALLOC(ddef->mocc, p, sizeof(float) * (L+1));
  ESL_RALLOC(ddef->btot, p, sizeof(float) * (L+1));
  ESL_RALLOC(ddef->etot, p, sizeof(float) * (L+1));
  ESL_RALLOC(ddef->n2sc, p, sizeof(float) * (L+1));
  ddef->Lalloc = L;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_domaindef_Reuse()
 * Synopsis:  Prepare to reuse a <P7_DOMAINDEF> on a new sequence.
 * Incept:    SRE, Fri Jan 25 13:48:36 2008 [Janelia]
 *
 * Purpose:   Prepare a <P7_DOMAINDEF> object <ddef> to be reused on
 *            a new sequence, reusing as much memory as possible.
 *            
 * Note:      Because of the way we handle alidisplays, handing them off to
 *            the caller, we don't reuse their memory; any unused
 *            alidisplays are destroyed. It's not really possible to
 *            reuse alidisplay memory. We need alidisplays to persist
 *            until all sequences have been processed and we're
 *            writing our final output to the user.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_domaindef_Reuse(P7_DOMAINDEF *ddef)
{
  int status;
  int d;

  /* If ddef->dcl is NULL, we turned the domain list over to a P7_HIT
   * for permanent storage, and we need to allocate a new one;
   * else, reuse the one we've got.
   */
  if (ddef->dcl == NULL)  
    ESL_ALLOC(ddef->dcl, sizeof(P7_DOMAIN) * ddef->nalloc);
  else
    {
      for (d = 0; d < ddef->ndom; d++) {
  p7_alidisplay_Destroy(ddef->dcl[d].ad); ddef->dcl[d].ad             = NULL;
  free(ddef->dcl[d].scores_per_pos);      ddef->dcl[d].scores_per_pos = NULL;
      }
      
    }
  ddef->ndom = 0;
  ddef->L    = 0;

  ddef->nexpected  = 0.0;
  ddef->nregions   = 0;
  ddef->nclustered = 0;
  ddef->noverlaps  = 0;
  ddef->nenvelopes = 0;

  p7_spensemble_Reuse(ddef->sp);
  p7_trace_Reuse(ddef->tr);  /* probable overkill; should already have been called */
  p7_trace_Reuse(ddef->gtr);  /* likewise */
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_domaindef_DumpPosteriors()
 * Synopsis:  Output posteriors that define domain structure to a stream.
 * Incept:    SRE, Fri Feb 29 08:32:14 2008 [Janelia]
 *
 * Purpose:   Output the vectors from <ddef> that are used to
 *            define domain structure to a stream <ofp>, in xmgrace format.
 *            
 *            There are four vectors. The first set is
 *            <mocc[1..i..L]>, the probability that residue <i> is
 *            emitted by the core model (is in a domain). The second
 *            set is <btot[1..i..L]>, the cumulative expected number
 *            of times that a domain uses a B state (starts) at or
 *            before position <i>. The third set is <etot[1..i..L]>,
 *            the cumulative expected number of times that a domain
 *            uses an E state (ends) at or before position <i>. The
 *            fourth set is <n2sc[1..i..L]>, the score of residue i
 *            under the ad hoc null2 model; this is a measure of local
 *            biased composition.
 *            
 *            These three fields will only be available after a call
 *            to domain definition by
 *            <p7_domaindef_ByPosteriorHeuristics()>.
 *
 * Returns:   <eslOK> on success
 *            
 * Xref:      J2/126
 */
int
p7_domaindef_DumpPosteriors(FILE *ofp, P7_DOMAINDEF *ddef)
{
  int i;
  fprintf(ofp, "# mocc btot etot n2sc\n");
  for (i = 0; i <= ddef->L; i++) {
    fprintf(ofp, "%d %f ", i, ddef->mocc[i]);
    fprintf(ofp, "%f ", ddef->btot[i]);
    fprintf(ofp, "%f ", ddef->etot[i]);
    fprintf(ofp, "%f ", ddef->n2sc[i]);
    fprintf(ofp, "\n");
  }
  return eslOK;
}


/* Function:  p7_domaindef_Destroy_BATH() 
 * Synopsis:  Destroys a <P7_DOMAINDEF>.
 * Incept:    SRE, Fri Jan 25 13:52:46 2008 [Janelia]
 *
 * Purpose:   Destroys a <P7_DOMAINDEF>.
 */
void
p7_domaindef_Destroy_BATH(P7_DOMAINDEF *ddef)
{
  int d;
  if (ddef == NULL) return;

  if (ddef->mocc != NULL) free(ddef->mocc);
  if (ddef->btot != NULL) free(ddef->btot);
  if (ddef->etot != NULL) free(ddef->etot);
  if (ddef->n2sc != NULL) free(ddef->n2sc);

  if (ddef->dcl  != NULL) {
    for (d = 0; d < ddef->ndom; d++) {
      if (ddef->dcl[d].scores_per_pos) free(ddef->dcl[d].scores_per_pos);
      p7_alidisplay_Destroy(ddef->dcl[d].ad);
    }
    free(ddef->dcl);
  }

  p7_spensemble_Destroy(ddef->sp);
  p7_trace_fs_Destroy(ddef->tr);
  p7_trace_fs_Destroy(ddef->gtr);
  free(ddef);
  return;
}

/*****************************************************************
 * 2. Routines inferring domain structure of a target sequence
 *****************************************************************/

/* Function:  p7_domaindef_ByPosteriorHeuristics_Frameshift_BATH() 
 * Synopsis:  Define "domains" in a DNA window using posterior probs
 *            with frameshift awareness.
 *
 * Purpose:   Given a DNA sequence <sq> and protien model (<gm>, <gm_fs>) 
 *            for which we have already calculated a Forward and 
 *            Backward parsing matrices <gxf> and <gxb>; use posterior 
 *            probability heuristics to determine an annotated "domain". 
 *            In this context domains are simply subsequences of the DNA 
 *            that have demonstrated a high prbability of homology to the 
 *            core modle. For each domain found, score it (with null2
 *            calculations) and obtain an optimal accuracy alignment,
 *            using <fwd> and <bck> matrices as workspace for the
 *            necessary full-matrix DP calculations. Caller provides a
 *            new or reused <ddef> object to hold these results. 
 *            Upon return, <ddef> contains the definitions of all the
 *            domains: their bounds, their null-corrected Forward
 *            scores, and their optimal posterior accuracy alignments.
 *            
 * Returns:   <eslOK> on success.           
 *           
 * Throws:    <eslEMEM> on allocation failure. 
 */
int
p7_domaindef_ByPosteriorHeuristics_Frameshift_BATH(ESL_SQ *windowsq, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, 
           P7_GMX *gxf, P7_GMX *gxb, P7_GMX *fwd, P7_GMX *bck, P7_IVX *iv, P7_DOMAINDEF *ddef, P7_BG *bg, 
	  ESL_GENCODE *gcode, int64_t window_start, int do_biasfilter
)
{

  int i, j;
  int triggered;
  int start;
  int end;
  int d;
  int i2,j2;
  int last_j2;
  int nc;
  int saveL     = gm_fs->L;     /* Save the length config of <gm_fs>; will restore upon return */
  int save_mode = gm_fs->mode;  /* Likewise for the mode. */
  int status;
 
  if ((status = p7_domaindef_GrowTo(ddef, windowsq->n))      != eslOK) return status;          /* ddef's btot,etot,mocc now ready for seq of length n            */
  if ((status = p7_DomainDecoding_Frameshift(gm_fs, gxf, gxb, ddef)) != eslOK) return status;  /* ddef->{btot,etot,mocc} now made.                               */

  esl_vec_FSet(ddef->n2sc, windowsq->n+1, 0.0);                                                /* ddef->n2sc null2 scores are initialized                        */
  ddef->nexpected = ddef->btot[windowsq->n];                                                   /* posterior expectation for # of domains (same as etot[sq->n])   */
  p7_fs_ReconfigUnihit(gm_fs, saveL/3);                                                          /* process each domain in unihit mode, regardless of om->mode     */

  i         = -1;
  triggered = FALSE;
  start     = FALSE;
  end       = FALSE;

  for (j = 1; j < gxf->L; j++)
  { 
    if (! triggered)  // no domain start found
    {
      if       (ddef->mocc[j]  >= ddef->rt1 ) triggered = TRUE; // possible domain start found
      d = j;
    }
    else 
    {
      /* Frameshift aware domain starts must be evident in all three frames */
      while(d > 1 && !start ) 
      { 
        d--;
        if(d > 3 && ddef->mocc[d] - (ddef->btot[d] - ddef->btot[d-3]) < ddef->rt2) 
        { 
          d--;
          if(d > 3 && ddef->mocc[d] - (ddef->btot[d] - ddef->btot[d-3]) < ddef->rt2)
          { 
            d--;
            if(d > 3 && ddef->mocc[d] - (ddef->btot[d] - ddef->btot[d-3]) < ddef->rt2) 
            { 
              d--;
              start = TRUE;
            }
          } 
        }
      }
      
      i = ESL_MAX(1, d-3); // one codon wiggle room
      d = j+1;

      /* Frameshift aware domain ends must be evident in all three frames */
      while(d < gxf->L && !end) 
      {
        d++;
        if(d < gxf->L && ddef->mocc[d] - (ddef->etot[d] - ddef->etot[d-3]) < ddef->rt2)
        {
          d++;
          if(d < gxf->L && ddef->mocc[d] - (ddef->etot[d] - ddef->etot[d-3]) < ddef->rt2) 
          {
            d++;
            if(d < gxf->L && ddef->mocc[d] - (ddef->etot[d] - ddef->etot[d-3]) < ddef->rt2) 
            {
              d++;
              end = TRUE;  
            }
          }
        }
      }
      
      j = ESL_MIN(gxf->L, d+3); // one codon wiggle room

      if(j - i + 1 < 12) {
        i     = -1;
        triggered = FALSE;
        start     = FALSE;
        end       = FALSE; 
        continue; 
      }

      /* We have a region i..j to evaluate. */
      p7_gmx_fs_GrowTo(fwd, gm_fs->M, j-i+1, j-i+1, p7P_5CODONS);
      p7_gmx_fs_GrowTo(bck, gm_fs->M, j-i+1, j-i+1, 0);
      ddef->nregions++;

      if (is_multidomain_region_fs(ddef, i, j))
      {
	
        /* This region appears to contain more than one domain, so we have to
        * resolve it by cluster analysis of posterior trace samples, to define
        * one or more domain envelopes.
        */
        ddef->nclustered++;
       
       /* Resolve the region into domains by stochastic trace
        * clustering; assign position-specific null2 model by
        * stochastic trace clustering; there is redundancy
        * here; we will consolidate later if null2 strategy
        * works
        */
        p7_ivx_GrowTo(iv, gm_fs->M, p7P_5CODONS); 
        p7_fs_ReconfigMultihit(gm_fs, saveL/3);
        p7_Forward_Frameshift(windowsq->dsq+i-1, gcode, j-i+1, gm_fs, fwd, iv, NULL);

        region_trace_ensemble_frameshift(ddef, gm_fs, windowsq->dsq, windowsq->abc, i, j, fwd, bck, &nc);

        p7_fs_ReconfigUnihit(gm_fs, saveL/3);
       
        /* ddef->n2sc is now set on i..j by the traceback-dependent method */
        last_j2 = 0;
	
        for (d = 0; d < nc; d++) {
         
          p7_spensemble_GetClusterCoords(ddef->sp, d, &i2, &j2, NULL, NULL, NULL);
         if (i2 <= last_j2) ddef->noverlaps++;

         /* Note that k..m coords on model are available, but
          * we're currently ignoring them.  This leads to a
          * rare clustering bug that we eventually need to fix
          * properly [xref J3/32]: two different regions in one
          * profile HMM might have hit same seq domain, and
          * when we now go to calculate an OA trace, nothing
          * constrains us to find the two different alignments
          * to the HMM; in fact, because OA is optimal, we'll
          * find one and the *same* alignment, leading to an
          * apparent duplicate alignment in the output.
          *
          * Registered as #h74, Dec 2009, after EBI finds and
          * reports it.  #h74 is worked around in p7_tophits.c
          * by hiding all but one envelope with an identical
          * alignment, in the rare event that this
          * happens. [xref J5/130].
          */
          ddef->nenvelopes++;         
         
         i2 = ESL_MAX(1,i2); // Hacky bug fix to prevent 0 index - real fix requires changes to region_trace_ensemble_frameshift() 

         if (rescore_isolated_domain_frameshift(ddef, gm, gm_fs, windowsq, fwd, bck, iv, i2, j2, FALSE, bg, gcode, do_biasfilter) == eslOK) last_j2 = j2;
        }

        p7_spensemble_Reuse(ddef->sp);
        p7_trace_Reuse(ddef->tr);

     } else {
	
      ddef->nenvelopes++;
       
      rescore_isolated_domain_frameshift(ddef, gm, gm_fs, windowsq, fwd, bck, iv, i, j, FALSE, bg, gcode, do_biasfilter);
    }
    i     = -1;
    triggered = FALSE;
    start     = FALSE;
    end       = FALSE;

   } 
  }
  /* Restore model to uni/multihit mode, and to its original length model */
  if (p7_IsMulti(save_mode)) p7_fs_ReconfigMultihit(gm_fs, saveL/3); 
  else                       p7_fs_ReconfigUnihit(gm_fs, saveL/3); 

  return eslOK;
}


/* Function:  p7_domaindef_ByPosteriorHeuristics_BATH() 
 * Synopsis:  Define domains in a sequence using posterior probs.
 * Incept:    SRE, Sat Feb 23 08:17:44 2008 [Janelia]
 *
 * Purpose:   Given an ORF <orfsq> and model <om> for which we have
 *            already calculated a Forward and Backward parsing
 *            matrices <oxf> and <oxb>; use posterior probability
 *            heuristics to determine an annotated domain structure;
 *            and for each domain found, score it (with null2
 *            calculations) and obtain an optimal accuracy alignment,
 *            using <fwd> and <bck> matrices as workspace for the
 *            necessary full-matrix DP calculations. Caller provides a
 *            new or reused <ddef> object to hold these results.
 *            Upon return, <ddef> contains the definitions of all the
 *            domains: their bounds, their null-corrected Forward
 *            scores, and their optimal posterior accuracy alignments.
 *             
 * Returns:   <eslOK> on success.           
 *            
 *            <eslERANGE> on numeric overflow in posterior
 *            decoding. This should not be possible for multihit
 *            models.
 */
int
p7_domaindef_ByPosteriorHeuristics_BATH(const ESL_SQ *orfsq, const ESL_SQ *windowsq, const int64_t ntsqlen, const ESL_GENCODE *gcode, 
	   P7_OPROFILE *om, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, P7_OMX *tmp_fwd, P7_OMX *fwd, P7_OMX *bck, P7_DOMAINDEF *ddef, P7_BG *bg)
{
  int i, j;
  int triggered;
  int d;
  int i2,j2;
  int last_j2;
  int nc;
  int saveL     = om->L;  /* Save the length config of <om>; will restore upon return */
  int save_mode = om->mode;  /* Likewise for the mode. */
  int status;

  if ((status = p7_domaindef_GrowTo(ddef, orfsq->n))      != eslOK) return status;  /* ddef's btot,etot,mocc now ready for seq of length n */
  if ((status = p7_DomainDecoding(om, tmp_fwd, bck, ddef)) != eslOK) return status;  /* ddef->{btot,etot,mocc} now made.                    */

  esl_vec_FSet(ddef->n2sc, orfsq->n+1, 0.0);          /* ddef->n2sc null2 scores are initialized                        */
  ddef->nexpected = ddef->btot[orfsq->n];             /* posterior expectation for # of domains (same as etot[orfsq->n])   */

  p7_oprofile_ReconfigUnihit(om, saveL);     /* process each domain in unihit mode, regardless of om->mode     */
  i     = -1;
  triggered = FALSE;

  for (j = 1; j <= orfsq->n; j++)
  {

    if (! triggered)
    {      /* xref J2/101 for what the logic below is: */

      if       (ddef->mocc[j] - (ddef->btot[j] - ddef->btot[j-1]) <  ddef->rt2) i = j;
      else if  (i == -1)                                                        i = j;
      if       (ddef->mocc[j]                                     >= ddef->rt1) triggered = TRUE;
    }
    else if (ddef->mocc[j] - (ddef->etot[j] - ddef->etot[j-1])  <  ddef->rt2)
    {
        /* We have a region i..j to evaluate. */
        p7_omx_GrowTo(fwd, om->M, j-i+1, j-i+1);
        p7_omx_GrowTo(bck, om->M, j-i+1, j-i+1);
        ddef->nregions++;
        
        if (is_multidomain_region(ddef, i, j))
        {  
        /* This region appears to contain more than one domain, so we have to
             * resolve it by cluster analysis of posterior trace samples, to define
             * one or more domain envelopes.
             */
            ddef->nclustered++;

            /* Resolve the region into domains by stochastic trace
             * clustering; assign position-specific null2 model by
             * stochastic trace clustering; there is redundancy
             * here; we will consolidate later if null2 strategy
             * works
             */
            p7_oprofile_ReconfigMultihit(om, saveL);
            p7_Forward(orfsq->dsq+i-1, j-i+1, om, fwd, NULL);
		
            region_trace_ensemble(ddef, om, orfsq->dsq, i, j, fwd, bck, &nc);
            p7_oprofile_ReconfigUnihit(om, saveL);
            /* ddef->n2sc is now set on i..j by the traceback-dependent method */
	    
            last_j2 = 0;
            if(nc == 0) ddef->nenvelopes++;
            for (d = 0; d < nc; d++) {
                  p7_spensemble_GetClusterCoords(ddef->sp, d, &i2, &j2, NULL, NULL, NULL);
                  if (i2 <= last_j2) ddef->noverlaps++;

                  /* Note that k..m coords on model are available, but
                     * we're currently ignoring them.  This leads to a
                     * rare clustering bug that we eventually need to fix
                     * properly [xref J3/32]: two different regions in one
                     * profile HMM might have hit same seq domain, and
                     * when we now go to calculate an OA trace, nothing
                     * constrains us to find the two different alignments
                     * to the HMM; in fact, because OA is optimal, we'll
                     * find one and the *same* alignment, leading to an
                     * apparent duplicate alignment in the output.
                     *
                     * Registered as #h74, Dec 2009, after EBI finds and
                     * reports it.  #h74 is worked around in p7_tophits.c
                     * by hiding all but one envelope with an identical
                     * alignment, in the rare event that this
                     * happens. [xref J5/130].
                  */
                  ddef->nenvelopes++;

                  /*the !long_target argument will cause the function to recompute null2
                   * scores if this is part of a long_target (nhmmer) pipeline */
		  
                  if (rescore_isolated_domain_bath(ddef, om, gm, gm_fs, orfsq, windowsq, ntsqlen, gcode, fwd, bck, i2, j2, TRUE, bg) == eslOK)
                       last_j2 = j2;

          }
            p7_spensemble_Reuse(ddef->sp);
            p7_trace_Reuse(ddef->tr);
        }
        else
        {
		
            /* The region looks simple, single domain; convert the region to an envelope. */
            ddef->nenvelopes++;
	    
            rescore_isolated_domain_bath(ddef, om, gm, gm_fs, orfsq, windowsq, ntsqlen, gcode, fwd, bck, i, j, FALSE, bg);
        }
        i     = -1;
        triggered = FALSE;
    }

  }

  /* Restore model to uni/multihit mode, and to its original length model */
  if (p7_IsMulti(save_mode)) p7_oprofile_ReconfigMultihit(om, saveL); 
  else                       p7_oprofile_ReconfigUnihit  (om, saveL); 
  return eslOK;
}

/*****************************************************************
 * 3. Internal routines 
 *****************************************************************/

/* is_multidomain_region()
 * SRE, Fri Feb  8 11:35:04 2008 [Janelia]
 *
 * This defines the trigger for when we need to hand a "region" off to
 * a deeper analysis (using stochastic tracebacks and clustering)
 * because there's reason to suspect it may encompass two or more
 * domains. 
 * 
 * The criterion is to find the split point z at which the expected
 * number of E occurrences preceding B occurrences is maximized, and
 * if that number is greater than the heuristic threshold <ddef->rt3>,
 * then return TRUE. In other words, we're checking to see if there's
 * any point in the region at which it looks like an E was followed by
 * a B, as expected for a multidomain interpretation of the region.
 * 
 * More precisely: return TRUE if  \max_z [ \min (B(z), E(z)) ]  >= rt3
 * where
 *   E(z) = expected number of E states occurring in region before z is emitted
 *        = \sum_{y=i}^{z} eocc[i]  =  etot[z] - etot[i-1]
 *   B(z) = expected number of B states occurring in region after z is emitted
 *        = \sum_{y=z}^{j} bocc[i]  =  btot[j] - btot[z-1]               
 *        
 *        
 * Because this relies on the <ddef->etot> and <ddef->btot> arrays,
 * <p7_DomainDecoding()> needs to have been called first.
 *
 * Xref:    J2/101.  
 */
static int
is_multidomain_region(P7_DOMAINDEF *ddef, int i, int j)
{
  int   z;
  float max;
  float expected_n;
  max = -1.0;
  for (z = i; z <= j; z++)
    { 
      expected_n = ESL_MIN( (ddef->etot[z] - ddef->etot[i-1]), (ddef->btot[j] - ddef->btot[z-1]) );
      max        = ESL_MAX(max, expected_n);
    }
  return ( (max >= ddef->rt3) ? TRUE : FALSE);
}

/* is_multidomain_region_fs() - BATH
 *
 * This function is supposed to define the trigger for when we need 
 * to hand a "region" of a DNA window off to a deeper analysis (using 
 * stochastic tracebacks and clustering) because there's reason to 
 * suspect it may encompass two or more seperate homologous regions. 
 * 
 * The criterion is to find the split point z at which the expected
 * number of E occurrences preceding B occurrences is maximized, and
 * if that number is greater than the heuristic threshold <ddef->rt3>,
 * then return TRUE. In other words, we're checking to see if there's
 * any point in the region at which it looks like an E was followed by
 * a B, as expected for a "multidomain" interpretation of the region.
 * 
 * More precisely: return TRUE if  \max_z [ \min (B(z), E(z)) ]  >= rt3
 * where
 *   E(z) = expected number of E states occurring in region before z is emitted
 *        = \sum_{y=i}^{z} eocc[i]  =  etot[z] - etot[i-1]
 *   B(z) = expected number of B states occurring in region after z is emitted
 *        = \sum_{y=z}^{j} bocc[i]  =  btot[j] - btot[z-1]               
 *        
 *        
 * Because this relies on the <ddef->etot> and <ddef->btot> arrays,
 * <p7_DomainDecoding_Frameshift()> needs to have been called first.
 *
 * Xref:    J2/101.  
 */
static int
is_multidomain_region_fs(P7_DOMAINDEF *ddef, int i, int j)
{
  int   z,f;
  float max;
  float expected_n;
  max = -1.0;
  
  f = (j-i+1) % 3;
  
  for (z = i+2; z <= j-f; z+=3)
    { 
      expected_n = ESL_MIN( (ddef->etot[z] - ddef->etot[i-1]), (ddef->btot[j-f] - ddef->btot[z-3]) );
      max        = ESL_MAX(max, expected_n);
  }

  f = (j-i) % 3;
  for (z = i+3; z <= j-f; z+=3)
    {
      expected_n = ESL_MIN( (ddef->etot[z] - ddef->etot[i]), (ddef->btot[j-f] - ddef->btot[z-3]) );
      max        = ESL_MAX(max, expected_n);
   }
  f = (j-i-1) % 3;

  for (z = i+4; z <= j-f; z+=3)
    {
      expected_n = ESL_MIN( (ddef->etot[z] - ddef->etot[i+1]), (ddef->btot[j-f] - ddef->btot[z-3]) );
      max        = ESL_MAX(max, expected_n);
  }

 return ( (max >= ddef->rt3) ? TRUE : FALSE);
}


/* region_trace_ensemble()
 * SRE, Fri Feb  8 11:49:44 2008 [Janelia]
 *
 * Here, we've decided that region <ireg>..<jreg> in sequence <dsq> might be
 * composed of more than one domain, and we're going to use clustering
 * of a posterior ensemble of stochastic tracebacks to sort it out.
 * 
 * Caller provides a filled Forward matrix in <fwd> for the sequence
 * region <dsq+ireg-1>, length <jreg-ireg+1>, for the model <om>
 * configured in multihit mode with its target length distribution
 * set to the total length of <dsq>: i.e., the same model
 * configuration used to score the complete sequence (if it weren't
 * multihit, we wouldn't be worried about multiple domains).
 * 
 * Caller also provides a DP matrix in <wrk> containing at least one
 * row, for use as temporary workspace. (This will typically be the
 * caller's Backwards matrix, which we haven't yet used at this point
 * in the processing pipeline.)
 * 
 * Caller provides <ddef>, which defines heuristic parameters that
 * control the clustering, and provides working space for the
 * calculation and the answers. The <ddef->sp> object must have been
 * reused (i.e., it needs to be fresh; we're going to use it here);
 * the caller needs to Reuse() it specifically, because it can't just
 * Reuse() the whole <ddef>, when it's in the process of analyzing
 * regions.
 * 
 * Upon return, <*ret_nc> contains the number of clusters that were
 * defined.
 * 
 * The caller can retrieve info on each cluster by calling
 * <p7_spensemble_GetClusterCoords(ddef->sp...)> on the
 * <P7_SPENSEMBLE> object in <ddef>.
 * 
 * Other information on what's happened in working memory:
 * 
 * <ddef->n2sc[ireg..jreg]> now contains log f'(x_i) / f(x_i) null2 scores
 *    for each residue.
 *
 * <ddef->sp> gets filled in, and upon return, it's holding the answers 
 *    (the cluster definitions). When the caller is done retrieving those
 *    answers, it needs to <esl_spensemble_Reuse()> it before calling
 *    <region_trace_ensemble()> again.
 *    
 * <ddef->tr> is used as working memory for sampled traces.
 *    
 * <wrk> has had its zero row clobbered as working space for a null2 calculation.
 */
static int
region_trace_ensemble(P7_DOMAINDEF *ddef, const P7_OPROFILE *om, const ESL_DSQ *dsq, int ireg, int jreg, 
          const P7_OMX *fwd, P7_OMX *wrk, int *ret_nc)
{
  int    Lr  = jreg-ireg+1;
  int    t, d, d2;
  int    nov, n;
  int    nc;
  int    pos;
  float  null2[p7_MAXCODE];

  esl_vec_FSet(ddef->n2sc+ireg, Lr, 0.0); /* zero the null2 scores in region */

  /* By default, we make results reproducible by forcing a reset of
   * the RNG to its originally seeded state.
   */
  if (ddef->do_reseeding) 
    esl_randomness_Init(ddef->r, esl_randomness_GetSeed(ddef->r));

  /* Collect an ensemble of sampled traces; calculate null2 odds ratios from these */
  for (t = 0; t < ddef->nsamples; t++)
    {
      p7_StochasticTrace(ddef->r, dsq+ireg-1, Lr, om, fwd, ddef->tr);
      p7_trace_Index(ddef->tr);

      pos = 1;
      for (d = 0; d < ddef->tr->ndom; d++)
  {
    p7_spensemble_Add(ddef->sp, t, ddef->tr->sqfrom[d]+ireg-1, ddef->tr->sqto[d]+ireg-1, ddef->tr->hmmfrom[d], ddef->tr->hmmto[d]);

    p7_Null2_ByTrace(om, ddef->tr, ddef->tr->tfrom[d], ddef->tr->tto[d], wrk, null2);
    /* residues outside domains get bumped +1: because f'(x) = f(x), so f'(x)/f(x) = 1 in these segments */
    for (; pos <= ddef->tr->sqfrom[d]; pos++) ddef->n2sc[ireg+pos-1] += 1.0;

    /* Residues inside domains get bumped by their null2 ratio */
    for (; pos <= ddef->tr->sqto[d];   pos++) ddef->n2sc[ireg+pos-1] += null2[dsq[ireg+pos-1]];
  }
      /* the remaining residues in the region outside any domains get +1 */
      for (; pos <= Lr; pos++)  ddef->n2sc[ireg+pos-1] += 1.0;

      p7_trace_Reuse(ddef->tr);        
    }

  /* Convert the accumulated n2sc[] ratios in this region to log odds null2 scores on each residue. */
  for (pos = ireg; pos <= jreg; pos++)
    ddef->n2sc[pos] = logf(ddef->n2sc[pos] / (float) ddef->nsamples);

  /* Cluster the ensemble of traces to break region into envelopes. */
  p7_spensemble_Cluster(ddef->sp, ddef->min_overlap, ddef->of_smaller, ddef->max_diagdiff, ddef->min_posterior, ddef->min_endpointp, &nc);

  /* A little hacky now. Remove "dominated" domains relative to seq coords. */
  for (d = 0; d < nc; d++) 
    ddef->sp->assignment[d] = 0; /* overload <assignment> to flag that a domain is dominated */

  /* who dominates who? (by post prob) */
  for (d = 0; d < nc; d++)
    {
      for (d2 = d+1; d2 < nc; d2++)
  {
    nov = ESL_MIN(ddef->sp->sigc[d].j, ddef->sp->sigc[d2].j) - ESL_MAX(ddef->sp->sigc[d].i, ddef->sp->sigc[d2].i) + 1;
    if (nov == 0) break;
    n   = ESL_MIN(ddef->sp->sigc[d].j - ddef->sp->sigc[d].i + 1,  ddef->sp->sigc[d2].j - ddef->sp->sigc[d2].i + 1);
    if ((float) nov / (float) n >= 0.8) /* overlap */
      {
        if (ddef->sp->sigc[d].prob > ddef->sp->sigc[d2].prob) ddef->sp->assignment[d2] = 1;
        else                                                  ddef->sp->assignment[d]  = 1;
      }
  }
    }
      
  /* shrink the sigc list, removing dominated domains */
  d = 0;
  for (d2 = 0; d2 < nc; d2++)
    {
      if (ddef->sp->assignment[d2]) continue; /* skip domain d2, it's dominated. */
      if (d != d2) memcpy(ddef->sp->sigc + d, ddef->sp->sigc + d2, sizeof(struct p7_spcoord_s));
      d++;
    }
  ddef->sp->nc = d;
  *ret_nc = d;
  return eslOK;
}

/* region_trace_ensemble_frameshift() - BATH
 *
 * Here, we've decided that region <ireg>..<jreg> in dna window <dsq> 
 * might be composed of more than one homologous region or domain, and 
 * we're going to use clustering of a posterior ensemble of stochastic 
 * tracebacks to sort it out.
 * 
 * Caller provides a filled Forward matrix in <fwd> for the sequence
 * region <dsq+ireg-1>, length <jreg-ireg+1>, for the model <gm_fs>
 * configured in multihit mode with its target length distribution
 * set to the total length of <dsq>: i.e., the same model
 * configuration used to score the complete sequence (if it weren't
 * multihit, we wouldn't be worried about multiple domains).
 * 
 * Caller also provides a DP matrix in <wrk> containing at least one
 * row, for use as temporary workspace. (This will typically be the
 * caller's Backwards matrix, which we haven't yet used at this point
 * in the processing pipeline.)
 * 
 * Caller provides <ddef>, which defines heuristic parameters that
 * control the clustering, and provides working space for the
 * calculation and the answers. The <ddef->sp> object must have been
 * reused (i.e., it needs to be fresh; we're going to use it here);
 * the caller needs to Reuse() it specifically, because it can't just
 * Reuse() the whole <ddef>, when it's in the process of analyzing
 * regions.
 * 
 * Upon return, <*ret_nc> contains the number of clusters that were
 * defined.
 * 
 * The caller can retrieve info on each cluster by calling
 * <p7_spensemble_GetClusterCoords(ddef->sp...)> on the
 * <P7_SPENSEMBLE> object in <ddef>.
 * 
 * Other information on what's happened in working memory:
 * 
 * <ddef->n2sc[ireg..jreg]> now contains log f'(x_i) / f(x_i) null2 scores
 *    for each residue.
 *
 * <ddef->sp> gets filled in, and upon return, it's holding the answers 
 *    (the cluster definitions). When the caller is done retrieving those
 *    answers, it needs to <esl_spensemble_Reuse()> it before calling
 *    <region_trace_ensemble()> again.
 *    
 * <ddef->tr> is used as working memory for sampled traces.
 *    
 * <wrk> has had its zero row clobbered as working space for a null2 calculation.
 */
static int
region_trace_ensemble_frameshift(P7_DOMAINDEF *ddef, const P7_FS_PROFILE *gm_fs, const ESL_DSQ *dsq, const ESL_ALPHABET *abc, int ireg, int jreg, const P7_GMX *fwd, P7_GMX *wrk, int *ret_nc)
{
  int    Lr  = jreg-ireg+1;
  int    t, d, d2;
  int    nov, n;
  int    nc;
  float   n2sc[Lr];

  esl_vec_FSet(n2sc, Lr, 0.0); /* zero the null2 scores in region */

  /* By default, we make results reproducible by forcing a reset of
   * the RNG to its originally seeded state.
   */
  if (ddef->do_reseeding) 
    esl_randomness_Init(ddef->r, esl_randomness_GetSeed(ddef->r));
  /* Collect an ensemble of sampled traces; calculate null2 odds ratios from these */
  for (t = 0; t < ddef->nsamples; t++)
    {

      p7_StochasticTrace_Frameshift(ddef->r, dsq+ireg-1, Lr, gm_fs, fwd, ddef->tr);

      p7_trace_fs_Index(ddef->tr);
 
    for (d = 0; d < ddef->tr->ndom; d++)
      p7_spensemble_Add(ddef->sp, t, ddef->tr->sqfrom[d]+ireg-1, ddef->tr->sqto[d]+ireg-1, ddef->tr->hmmfrom[d], ddef->tr->hmmto[d]);
     
     p7_trace_Reuse(ddef->tr);   
  }

  /* Cluster the ensemble of traces to break region into envelopes. */
  p7_spensemble_fs_Cluster(ddef->sp, ddef->min_overlap, ddef->of_smaller, ddef->max_diagdiff, ddef->min_posterior, ddef->min_endpointp, &nc);
 	
  /* A little hacky now. Remove "dominated" domains relative to seq coords. */
  for (d = 0; d < nc; d++) 
    ddef->sp->assignment[d] = 0; /* overload <assignment> to flag that a domain is dominated */

  /* who dominates who? (by post prob) */
  for (d = 0; d < nc; d++)
  {
    for (d2 = d+1; d2 < nc; d2++)
    {
    nov = ESL_MIN(ddef->sp->sigc[d].j, ddef->sp->sigc[d2].j) - ESL_MAX(ddef->sp->sigc[d].i, ddef->sp->sigc[d2].i) + 1;
    if (nov == 0) break;
    n   = ESL_MIN(ddef->sp->sigc[d].j - ddef->sp->sigc[d].i + 1,  ddef->sp->sigc[d2].j - ddef->sp->sigc[d2].i + 1);
    if ((float) nov / (float) n >= 0.8) /* overlap */
      {
        if (ddef->sp->sigc[d].prob > ddef->sp->sigc[d2].prob) ddef->sp->assignment[d2] = 1;
        else                                                  ddef->sp->assignment[d]  = 1;
      }
    }
  }

  /* shrink the sigc list, removing dominated domains */
  d = 0;

  for (d2 = 0; d2 < nc; d2++)
    {
      if (ddef->sp->assignment[d2]) continue; /* skip domain d2, it's dominated. */
      if (d != d2) memcpy(ddef->sp->sigc + d, ddef->sp->sigc + d2, sizeof(struct p7_spcoord_s));
      d++;
    }
  ddef->sp->nc = d;
 
  *ret_nc = d;
  return eslOK;
}



 /* rescore_isolated_domain_frameshift() - BATH
 *
 * We have isolated a single domain's envelope from <i>..<j> in DNA 
 * window <windowsq>, and now we want to score it in isolation and 
 * obtain an alignment display for it.
 * 
 * The caller provides model <gm_fs> configured in unilocal mode; by
 * using unilocal (as opposed to multilocal), we're going to force the
 * identification of a single domain in this envelope now.
 * 
 * The alignment is an optimal accuracy alignment (sensu IH Holmes),
 * also obtained in unilocal mode.
 * 
 * The caller provides DP matrices <gx1> and <gx2> with sufficient
 * space to hold Forward and Backward calculations for this domain
 * against the model. (The caller will typically already have matrices
 * sufficient for the complete sequence lying around, and can just use
 * those.) A third matrix <gxppfs> will need to be created because the 
 * frameshift aware posterior probability algorithim does not allow 
 * gx2 to be overwriten. It will be destroyed again before exit. 
 *
 * The caller also provides a <P7_DOMAINDEF> object (ddef)
 * which is (efficiently, we trust) managing any necessary temporary
 * working space and heuristic thresholds.
 *
 * Returns <eslOK> if a domain was successfully identified, scored,
 * and aligned in the envelope; if so, the relavant information is
 * registered in <ddef>, in <ddef->dcl>.
 *
 * Throws:    <eslEMEM> on allocation failure. 
 * 
 */
static int
rescore_isolated_domain_frameshift(P7_DOMAINDEF *ddef, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, ESL_SQ *windowsq, 
		                   P7_GMX *gx1, P7_GMX *gx2, P7_IVX *iv, int i, int j, int null2_is_done, P7_BG *bg, 
				   ESL_GENCODE *gcode, int do_biasfilter)
{

  P7_DOMAIN     *dom           = NULL;
  P7_GMX        *gxppfs;
  int            Ld            = j-i+1;
  int            n_holder;
  float          domcorrection = 0.0;
  float          envsc, oasc;
  int            codon_idx;  
  int            z;
  int            pos;
  float          null2[p7_MAXCODE];
  int            status;
  ESL_DSQ        t, u, v, w, x;
  ESL_DSQ       *dsq_holder;

  if (Ld < 15) return eslOK;
  
  p7_fs_ReconfigLength(gm_fs, Ld/3);
  p7_ivx_GrowTo(iv, gm_fs->M, p7P_5CODONS); 
 
  dsq_holder = windowsq->dsq;
  windowsq->dsq = windowsq->dsq+i-1;
  n_holder = windowsq->n;
  windowsq->n = Ld;
  windowsq->L = Ld;
  
  windowsq->dsq = dsq_holder;
  windowsq->n = n_holder; 
  windowsq->L = n_holder;  
   
  /* Forward */ 
  p7_Forward_Frameshift(windowsq->dsq+i-1, gcode, Ld, gm_fs, gx1, iv, &envsc);
  
  /* Backward */
  p7_Backward_Frameshift(windowsq->dsq+i-1, gcode, Ld, gm_fs, gx2, iv, NULL);

  /* Posterior Probabilities */
  if ((gxppfs = p7_gmx_fs_Create(gm_fs->M, Ld, Ld, p7P_5CODONS)) == NULL) goto ERROR;
  p7_Decoding_Frameshift(gm_fs, gx1, gx2, gxppfs);      

  /* Find an optimal accuracy alignment */
  p7_OptimalAccuracy_Frameshift(gm_fs, gxppfs, gx2, &oasc);      
  p7_OATrace_Frameshift(gm_fs, gxppfs, gx2, gx1, ddef->tr);   /* <tr>'s seq coords are offset by i-1, rel to orig dsq */

  /* hack the trace's sq coords to be correct w.r.t. original dsq */
  for (z = 0; z < ddef->tr->N; z++)    
    if (ddef->tr->i[z] >= 0) ddef->tr->i[z] += i-1;
  
  /* get ptr to next empty domain structure in domaindef's results */
  if (ddef->ndom == ddef->nalloc) {
    ESL_REALLOC(ddef->dcl, sizeof(P7_DOMAIN) * (ddef->nalloc*2));
    ddef->nalloc *= 2;
  }

  dom = &(ddef->dcl[ddef->ndom]);
  dom->ad             = p7_alidisplay_fs_Create(ddef->tr, 0, gm, gm_fs, windowsq, gcode);
  dom->scores_per_pos = NULL; 
  
   /* Compute bias correction
   * Is null2 set already for this i..j? (It is, if we're in a domain that
   * was defined by stochastic traceback clustering in a multidomain region;
   * it isn't yet, if we're in a simple one-domain region).  -  this possibility
   * is currently out of order.
   *
   * If null2 is not set, do it now, by the expectation (posterior decoding) method.
   * Framshift aware bias score are usually low and quite often zero scince they 
   * take all three frames into account - just as the Forward score does
   */
  
  if (!null2_is_done)
  { 
    p7_Null2_fs_ByExpectation(gm_fs, gxppfs, null2);

    t = u = v = w = x = -1;
    z = 0;
    pos = i;
  
    while(pos <= j)
    {
      if(esl_abc_XIsCanonical(windowsq->abc, windowsq->dsq[pos])) x = windowsq->dsq[pos];
      else if(esl_abc_XIsDegenerate(windowsq->abc, windowsq->dsq[pos]))
      {
        for(x = 0; x < windowsq->abc->K; x++)
          if(windowsq->abc->degen[windowsq->dsq[pos]][x]) break;
      }      
      
      switch (ddef->tr->st[z]) {
        case p7T_N:
        case p7T_C:
        case p7T_J:  ddef->n2sc[pos]  = 0.0;
                     if(ddef->tr->i[z] == pos && pos > i+1) pos++;
                     z++;   break;
        case p7T_X:
        case p7T_S:
        case p7T_B:
        case p7T_E:
        case p7T_T:
        case p7T_D:  z++;   break;
        case p7T_M:  if(ddef->tr->i[z] == pos)
                     {  
                       if(ddef->tr->c[z] == 1)      { codon_idx = p7P_CODON1(x);             codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_QC2); }      
                       else if(ddef->tr->c[z] == 2) { codon_idx = p7P_CODON2(w, x);          codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_QC1); }
                       else if(ddef->tr->c[z] == 3) { codon_idx = p7P_CODON3(v, w, x);       codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_C); } 
                       else if(ddef->tr->c[z] == 4) { codon_idx = p7P_CODON4(u, v, w, x);    codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_QC1); }
                       else if(ddef->tr->c[z] == 5) { codon_idx = p7P_CODON5(t, u, v, w, x); codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_QC2); }
                       ddef->n2sc[pos]  = logf(null2[p7P_AMINO(gm_fs, ddef->tr->k[z], codon_idx)]);
                       z++; 
                     }
                     else 
                       ddef->n2sc[pos]  = 0.0;
                     pos++;  break;
        case p7T_I:  if(ddef->tr->i[z] == pos)
                     {
                       codon_idx = p7P_CODON3(v, w, x);       
                       codon_idx = p7P_MINIDX(codon_idx, p7P_DEGEN_C);
                       ddef->n2sc[pos]  = logf(null2[p7P_AMINO(gm_fs, ddef->tr->k[z], codon_idx)]);
                       z++;
                     }
                     else 
                        ddef->n2sc[pos]  = 0.0;
                     pos++;  break;
      }
      t = u;
      u = w;
      v = w;
      w = x;
    } 
  } 

  for (pos = i; pos <= j; pos++) 
    domcorrection   += ddef->n2sc[pos];         /* domcorrection is in units of NATS */
  
  dom->domcorrection = ESL_MAX(0., domcorrection); /* in units of NATS */
   
  if(windowsq->start < windowsq->end)
  {
    dom->iali          = dom->ad->sqfrom;
    dom->jali          = dom->ad->sqto;
    dom->ienv          = i;
    dom->jenv          = j;
 }
  else
  {
    dom->iali          = dom->ad->sqto - 2;
    dom->jali          = dom->ad->sqfrom;
    dom->ienv          = j;
    dom->jenv          = i;
 }
  
  dom->envsc         = envsc;         /* in units of NATS */
  dom->oasc          = oasc;        /* in units of expected # of correctly aligned residues */
  dom->dombias       = 0.0; /* gets set later, using bg->omega and dombias */
  dom->bitscore      = 0.0; /* gets set later by caller, using envsc, null score, and dombias */
  dom->lnP           = 0.0; /* gets set later by caller, using bitscore */
  dom->is_reported   = FALSE; /* gets set later by caller */
  dom->is_included   = FALSE; /* gets set later by caller */
  dom->tr            = p7_trace_fs_Clone(ddef->tr); 
  
  ddef->ndom++;
  p7_trace_Reuse(ddef->tr);
  p7_gmx_Destroy(gxppfs);
  
  return eslOK;

 ERROR:
  p7_trace_Reuse(ddef->tr);
  p7_gmx_Destroy(gxppfs);
  return eslEMEM;
}

/* rescore_isolated_domain_bath() 
 *
 * We have isolated a single domain's envelope from <i>..<j> in
 * translated ORF sequence <orfsq>, and now we want to score it 
 * in isolation and map it back to the DNA window to obtain an 
 * alignment display for it.
 * 
 * The caller provides model <om> configured in unilocal mode; by
 * using unilocal (as opposed to multilocal), we're going to force the
 * identification of a single domain in this envelope now. Models <gm> 
 * and <gm_fs> are also provided for use in creasting the alignment
 * display.
 * 
 * The alignment is an optimal accuracy alignment (sensu IH Holmes),
 * also obtained in unilocal mode.
 * 
 * The caller provides DP matrices <ox1> and <ox2> with sufficient
 * space to hold Forward and Backward calculations for this domain
 * against the model. (The caller will typically already have matrices
 * sufficient for the complete sequence lying around, and can just use
 * those.) The caller also provides a <P7_DOMAINDEF> object (ddef)
 * which is (efficiently, we trust) managing any necessary temporary
 * working space and heuristic thresholds.
 *
 * Returns <eslOK> if a domain was successfully identified, scored,
 * and aligned in the envelope; if so, the per-domain information is
 * registered in <ddef>, in <ddef->dcl>.
 * 
 * Returns <eslFAIL> if domain is not successfully identified.  This
 * is rare; one way it can happen is if posterior decoding calculation
 * overflows, which can occur on highly repetitive sequence
 * {J3/119-121}. Beware: as a result, it is possible to have
 * <ddef->ndom = 0>, for nonzero region(s)/envelope(s). See {iss131}.
 * 
 */
static int
rescore_isolated_domain_bath(P7_DOMAINDEF *ddef, P7_OPROFILE *om, P7_PROFILE *gm, P7_FS_PROFILE *gm_fs, 
		                      const ESL_SQ *orfsq, const ESL_SQ *windowsq, const int64_t ntsqlen, const ESL_GENCODE *gcode, 
		                      P7_OMX *ox1, P7_OMX *ox2, int i, int j, int null2_is_done, P7_BG *bg)
{

  P7_DOMAIN     *dom           = NULL;
  int            Ld            = j-i+1;
  float          domcorrection = 0.0;
  float          envsc, oasc, bcksc;
  int            z;
  int            pos;
  float          null2[p7_MAXCODE];
  int            status;
 
  p7_oprofile_ReconfigLength(om, orfsq->n);
  p7_omx_GrowTo(ox1, om->M, orfsq->n, orfsq->n); 
  p7_omx_GrowTo(ox2, om->M, orfsq->n, orfsq->n); 
	
  p7_Forward (orfsq->dsq + i-1, Ld, om,      ox1, &envsc);
  p7_Backward(orfsq->dsq + i-1, Ld, om, ox1, ox2, &bcksc);
   
  status = p7_Decoding(om, ox1, ox2, ox2);      /* <ox2> is now overwritten with post probabilities     */
  if (status == eslERANGE) return eslFAIL;      /* rare: numeric overflow; domain is assumed to be repetitive garbage [J3/119-121] */
 
  /* Find an optimal accuracy alignment */
  p7_OptimalAccuracy(om, ox2, ox1, &oasc);      /* <ox1> is now overwritten with OA scores              */

  p7_OATrace        (om, ox2, ox1, ddef->tr);   /* <tr>'s seq coords are offset by i-1, rel to orig dsq */
  
   
  /* get ptr to next empty domain structure in domaindef's results */
  if (ddef->ndom == ddef->nalloc) {
    ESL_REALLOC(ddef->dcl, sizeof(P7_DOMAIN) * (ddef->nalloc*2));
    ddef->nalloc *= 2;
  }
    
  /* hack the trace's sq coords to be correct w.r.t. original dsq */
  for (z = 0; z < ddef->tr->N; z++)
    if (ddef->tr->i[z] > 0) ddef->tr->i[z] += i-1;

  if(orfsq->start < orfsq->end)
    p7_trace_fs_Convert(ddef->tr, orfsq->start, windowsq->start);
  else
    p7_trace_fs_Convert(ddef->tr, ntsqlen - orfsq->start + 1, windowsq->start);

  dom = &(ddef->dcl[ddef->ndom]);
  dom->ad             = p7_alidisplay_fs_Create(ddef->tr, 0, gm, gm_fs, windowsq, gcode);
  
  dom->scores_per_pos = NULL;  

  if (!null2_is_done) {   
    p7_Null2_ByExpectation(om, ox2, null2);
    for (pos = i; pos <= j; pos++) 
      ddef->n2sc[pos]  = logf(null2[orfsq->dsq[pos]]);
  }
  
  for (pos = i; pos <= j; pos++)
    domcorrection   += ddef->n2sc[pos];

  dom->domcorrection = domcorrection; /* in units of NATS */
	
  dom->iali          = dom->ad->sqfrom;
  dom->jali          = dom->ad->sqto;
  dom->ienv          = i; 
  dom->jenv          = j; 

  dom->envsc         = envsc;         /* in units of NATS */
  dom->oasc          = oasc;          /* in units of expected # of correctly aligned residues */
  dom->dombias       = 0.0;           /* gets set later, using bg->omega and dombias */
  dom->bitscore      = 0.0;           /* gets set later by caller, using envsc, null score, and dombias */
  dom->lnP           = 0.0;           /* gets set later by caller, using bitscore */
  dom->is_reported   = FALSE;         /* gets set later by caller */
  dom->is_included   = FALSE;         /* gets set later by caller */
  dom->tr            = p7_trace_fs_Clone(ddef->tr); 

  ddef->ndom++;
  p7_trace_Reuse(ddef->tr);

  return eslOK;

 ERROR:
  p7_trace_Reuse(ddef->tr);
  return status;
}



   

