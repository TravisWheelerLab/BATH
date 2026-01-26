/* P7_GMX implementation: a generic dynamic programming matrix for a frameshift aware model
 *
 * Contents:
 *   1. The <P7_GMX> frameshift aware object
 *   2. The <P7_IVX> frameshift aware object
 *   3. Debugging aids
 *   4. Unit tests
 *   5. Test driver
 */
#include "p7_config.h"
#include "hmmer.h"

/*****************************************************************
 *= 1. The <P7_GMX> object.
 *****************************************************************/

/* Function:  p7_gmx_fs_Create()
 * Synopsis:  Allocate a new <P7_GMX>.
 *
 * Purpose:   Allocate a reusable, resizeable <P7_GMX> for models up to
 *            size <allocM> and sequences up to length <allocL> with 
 *            cells for all frameshift codons.
 *            
 *            We've set this up so it should be easy to allocate
 *            aligned memory, though we're not doing this yet.
 *
 * Returns:   a pointer to the new <P7_GMX>.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_GMX *
p7_gmx_fs_Create(int allocM, int allocL, int allocLx, int allocC)
{
  int     status;
  P7_GMX *gx = NULL;
  int     i;

  /* don't try to make large allocs on 32-bit systems */
  if ( (uint64_t) (allocM+1) * (uint64_t) (allocLx+1) * sizeof(float) * (p7G_NSCELLS + allocC) > SIZE_MAX / 2)
    return NULL;

  /* level 1: the structure itself */
  ESL_ALLOC(gx, sizeof(P7_GMX));
  gx->dp     = NULL;
  gx->xmx    = NULL;
  gx->dp_mem = NULL;

  /* level 2: row pointers, 0.1..L; and dp cell memory  */
  ESL_ALLOC(gx->dp,      sizeof(float *) * (allocL+1));
  ESL_ALLOC(gx->xmx,     sizeof(float)   * (allocLx+1) * p7G_NXCELLS);
  ESL_ALLOC(gx->dp_mem,  sizeof(float)   * (allocL+1) * (allocM+1) *  (p7G_NSCELLS + allocC));
  
  /* Set the row pointers. */
  for (i = 0; i <= allocL; i++)  
    gx->dp[i] = gx->dp_mem + i * (allocM+1) * (p7G_NSCELLS + allocC);
   
  /* Initialize memory that's allocated but unused, only to keep
   * valgrind and friends happy.
   */
  for (i = 0; i <= allocL; i++) 
    { 
      gx->dp[i][0      * (p7G_NSCELLS + allocC) + p7G_M + p7G_C0] = -eslINFINITY; /* M_0 Codon 0*/
      if(allocC) {// allocC should equal either 0 or 5
        gx->dp[i][0      * (p7G_NSCELLS + allocC) + p7G_M + p7G_C1] = -eslINFINITY; /* M_0 Codon 1*/
        gx->dp[i][0      * (p7G_NSCELLS + allocC) + p7G_M + p7G_C2] = -eslINFINITY; /* M_0 Codon 2*/
        gx->dp[i][0      * (p7G_NSCELLS + allocC) + p7G_M + p7G_C3] = -eslINFINITY; /* M_0 Codon 3*/
        gx->dp[i][0      * (p7G_NSCELLS + allocC) + p7G_M + p7G_C4] = -eslINFINITY; /* M_0 Codon 4*/
        gx->dp[i][0      * (p7G_NSCELLS + allocC) + p7G_M + p7G_C5] = -eslINFINITY; /* M_0 Codon 5*/
      }
      gx->dp[i][0      * (p7G_NSCELLS + allocC) + p7G_I] = -eslINFINITY; /* I_0 */      
      gx->dp[i][0      * (p7G_NSCELLS + allocC) + p7G_D] = -eslINFINITY; /* D_0 */
      gx->dp[i][1      * (p7G_NSCELLS + allocC) + p7G_D] = -eslINFINITY; /* D_1 */
      gx->dp[i][allocM * (p7G_NSCELLS + allocC) + p7G_I] = -eslINFINITY; /* I_M */
    }
 
  gx->M      = 0;
  gx->L      = 0;
  gx->allocW = allocM+1;
  gx->allocR = allocLx+1;
  gx->validR = allocL+1;
  gx->allocC = allocC;
  gx->ncells = (uint64_t) (allocM+1)* (uint64_t) (allocLx+1);

  return gx;

 ERROR:
  if (gx != NULL) p7_gmx_Destroy(gx);
  return NULL;
}

/* Function:  p7_gmx_fs_GrowTo()
 * Synopsis:  Assure that DP matrix is big enough.
 *
 * Purpose:   Assures that a DP matrix <gx> is allocated
 *            for a model of size up to <M> and a sequence of
 *            length up to <L> with all frameshift cells 
 *            included; reallocates if necessary.
 *            
 *            This function does not respect the configured
 *            <RAMLIMIT>; it will allocate what it's told to
 *            allocate. 
 *
 * Returns:   <eslOK> on success, and <gx> may be reallocated upon
 *            return; any data that may have been in <gx> must be 
 *            assumed to be invalidated.
 *
 * Throws:    <eslEMEM> on allocation failure, and any data that may
 *            have been in <gx> must be assumed to be invalidated.
 */
int
p7_gmx_fs_GrowTo(P7_GMX *gx, int M, int L, int Lx, int C)
{
  int      status;
  void    *p;
  int      i;
  uint64_t ncells;
  int      do_reset = FALSE;

 if (M < gx->allocW && L < gx->validR && Lx < gx->allocR && C <= gx->allocC) return eslOK;

  /* must we realloc the 2D matrices? (or can we get away with just
   * jiggering the row pointers, if we are growing in one dimension
   * while shrinking in another?)
   */
  ncells = (uint64_t) (M+1) * (uint64_t) (Lx+1);

  /* don't try to make large allocs on 32-bit systems */
  if ( ncells * sizeof(float) * (p7G_NSCELLS + ESL_MAX(C, gx->allocC)) > SIZE_MAX / 2) return eslEMEM;

  if (ncells > gx->ncells || C > gx->allocC) 
    {
      ESL_RALLOC(gx->dp_mem, p, sizeof(float) * ncells * (p7G_NSCELLS + C));
      gx->ncells = ncells;
      gx->allocC = C;
      do_reset   = TRUE;
    }
  /* must we reallocate the row pointers? */
  if (Lx >= gx->allocR)
    {
      ESL_RALLOC(gx->xmx, p, sizeof(float)   * (Lx+1) * p7G_NXCELLS);
      gx->allocR = Lx+1;		/* allocW will also get set, in the do_reset block */
    }

   if(L >= gx->validR)
   {
     ESL_RALLOC(gx->dp,  p, sizeof(float *) * (L+1));
	 gx->validR = L+1; 
     do_reset   = TRUE;
   } 

  /* must we widen the rows? */
  if (M >= gx->allocW) do_reset = TRUE;

  /* resize the rows and reset all the valid row pointers.*/
  if (do_reset)
    {
      gx->allocW = M+1;
      for (i = 0; i < gx->validR; i++) 
	    gx->dp[i] = gx->dp_mem + i * (gx->allocW) * (p7G_NSCELLS + gx->allocC);
    }

  gx->M      = 0;
  gx->L      = 0;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_gmx_fs_Sizeof()
 * Synopsis:  Returns the allocation size of DP matrix, in bytes.
 */
size_t
p7_gmx_fs_Sizeof(P7_GMX *gx)
{
  size_t n = 0;

  n += sizeof(P7_GMX);
  n += gx->ncells * (p7G_NSCELLS + gx->allocC) * sizeof(float); /* main dp cells: gx->dp_mem */
  n += gx->validR * sizeof(float *);             /* row ptrs:      gx->dp[]   */
  n += gx->allocR * p7G_NXCELLS * sizeof(float); /* specials:      gx->xmx    */
  return n;
}

/*****************************************************************
 * 2. The <P7_IVX> object.
 *****************************************************************/

/* Function:  p7_ivx_Create()
 * Synopsis:  Allocate a new <P7_IMX>.
 *
 * Purpose:   Allocate a reusable, resizeable <P7_IVX> for models up to
 *            size <allocM>
 *
 * Returns:   a pointer to the new <P7_IVX>.
 *
 * Throws:    <NULL> on allocation error.
 */
P7_IVX *
p7_ivx_Create(int allocM, int allocC)
{
  int     status;
  P7_IVX *iv = NULL;

  ESL_ALLOC(iv, sizeof(P7_IVX));
  iv->ivx = NULL;

  ESL_ALLOC(iv->ivx, sizeof(float) * (allocM+1) * allocC);

  iv->allocM = allocM;
  iv->allocC = allocC;

  return iv;

 ERROR:
  if (iv != NULL) p7_ivx_Destroy(iv);
  return NULL;
}


/* Function:  p7_ivx_GrowTo()
 * Synopsis:  Assure that intermadiate values matrix is big enough.
 *
 * Returns:   <eslOK> on success, and <iv> may be reallocated upon
 *            return; any data that may have been in <iv> must be
 *            assumed to be invalidated.
 *
 * Throws:    <eslEMEM> on allocation failure, and any data that may
 *            have been in <gx> must be assumed to be invalidated.
 */
int
p7_ivx_GrowTo(P7_IVX *iv, int M, int C)
{
  int      status;
  void    *p;

  if((M+1) * C > (iv->allocM+1) * iv->allocC) {
    ESL_RALLOC(iv->ivx, p, sizeof(float) * (M+1) * C);
    iv->allocM = M;
    iv->allocC = C;
  }

  return eslOK;

 ERROR:
  return status;
}

/* Function:  p7_ivx_Destroy()
 * Synopsis:  Frees an intermediate values matrix.
 *
 * Purpose:   Frees a <P7_IVX>.
 *
 * Returns:   (void)
 */
void
p7_ivx_Destroy(P7_IVX *iv)
{
  if (iv == NULL) return;

  if (iv->ivx != NULL) free(iv->ivx);
  free(iv);
  return;
}



/*****************************************************************
 * 3. Debugging aids
 *****************************************************************/


/* Function:  p7_gmx_fs_Compare()
 * Synopsis:  Compare two Frameshift DP matrices for equality within given tolerance.
 *
 * Purpose:   Compare all the values in DP matrices <gx1> and <gx2> using
 *            <esl_FCompare_old()> and relative epsilon <tolerance>. If any
 *            value pairs differ by more than the acceptable <tolerance>
 *            return <eslFAIL>.  If all value pairs are identical within
 *            tolerance, return <eslOK>.
 */
int
p7_gmx_fs_Compare(P7_GMX *gx1, P7_GMX *gx2, int C, float tolerance)
{
  int i,k,c,x;
  if (gx1->M != gx2->M) return eslFAIL;
  if (gx1->L != gx2->L) return eslFAIL;
  
  for (i = 0; i <= gx1->L; i++)
  {
      for (k = 1; k <= gx1->M; k++) /* k=0 is a boundary; doesn't need to be checked */
      {
        if (esl_FCompare_old(gx1->dp[i][k * (p7G_NSCELLS+C) + p7G_M],  gx2->dp[i][k * (p7G_NSCELLS+C) + p7G_M], tolerance) != eslOK) return eslFAIL;
        if (esl_FCompare_old(gx1->dp[i][k * (p7G_NSCELLS+C) + p7G_I],  gx2->dp[i][k * (p7G_NSCELLS+C) + p7G_I], tolerance) != eslOK) return eslFAIL;
        if (esl_FCompare_old(gx1->dp[i][k * (p7G_NSCELLS+C) + p7G_D],  gx2->dp[i][k * (p7G_NSCELLS+C) + p7G_D], tolerance) != eslOK) return eslFAIL;
        for(c = p7G_NSCELLS; c < (p7G_NSCELLS+C); c++)
           if (esl_FCompare_old(gx1->dp[i][k * (p7G_NSCELLS+C) + p7G_M + c],  gx2->dp[i][k * (p7G_NSCELLS+C) + c], tolerance) != eslOK) return eslFAIL;
      }
      for (x = 0; x < p7G_NXCELLS; x++)
        if (esl_FCompare_old(gx1->xmx[i * p7G_NXCELLS + x], gx2->xmx[i * p7G_NXCELLS + x], tolerance) != eslOK) return eslFAIL;
  }
  return eslOK;
}


/* Function:  p7_gmx_fs_Dump()
 * Synopsis:  Dump a frameshift DP matrix to a stream, for diagnostics.
 *
 * Purpose:   Dump matrix <gx> to stream <fp> for diagnostics.
 *
 *            <flags> control some optional output behaviors, as follows:
 *              | <p7_HIDE_SPECIALS> | don't show scores for <ENJBC> states  |
 *              | <p7_SHOW_LOG>      | <gx> is in probs; show as log probs   |
 *
 *            int sci_note = TRUE to print values in scientific notation
 */
int
p7_gmx_fs_Dump(FILE *ofp, P7_GMX *gx, int flags, int scientific)
{
  if (scientific) return p7_gmx_fs_DumpWindow_Scientific(ofp, gx, 0, gx->L, 0, gx->M, flags);
  else            return p7_gmx_fs_DumpWindow(ofp, gx, 0, gx->L, 0, gx->M, flags);
}


/* Function:  p7_gmx_fs_DumpWindow()
 * Synopsis:  Dump a window of a frameshift DP matrix to a stream for diagnostics.
 *
 * Purpose:   Dump a window of matrix <gx> to stream <fp> for diagnostics,
 *            from row <istart> to <iend>, from column <kstart> to <kend>.
 *            
 *            Asking for <0..L,0..M> with <flags=p7_SHOW_SPECIALS> is the
 *            same as <p7_gmx_Dump()>.
 *            
 *            <flags> control some optional output behaviors, as follows:
 *              | <p7_HIDE_SPECIALS> | don't show scores for <ENJBC> states  |
 *              | <p7_SHOW_LOG>      | <gx> is in probs; show as log probs   |
 *  
 * Returns:   <eslOK> on success.
 */
int
p7_gmx_fs_DumpWindow(FILE *ofp, P7_GMX *gx, int istart, int iend, int kstart, int kend, int flags)
{
  int   width     = 9;
  int   precision = 4;
  int   i, k, c, x;
  float val;

  /* Header */
  fprintf(ofp, "     ");
  for (k = kstart; k <= kend;  k++) {
    fprintf(ofp, "%*d_0", width, k);
    for(c = 1; c <= gx->allocC; c++) 
      fprintf(ofp, "%*d_%d", width, k, c);
  }

  if (! (flags & p7_HIDE_SPECIALS)) fprintf(ofp, "%*s %*s %*s %*s %*s\n", width, "E", width, "N", width, "J", width, "B", width, "C");
  fprintf(ofp, "      ");
  for (k = kstart; k <= kend; k++) { 
    fprintf(ofp, "%*.*s ", width, width, "----------");
    for (c = 1; c <= gx->allocC; c++)  	
      fprintf(ofp, "%*.*s ", width, width, "----------");
  }

  if (! (flags & p7_HIDE_SPECIALS)) 
    for (x = 0; x < 5; x++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");
  
  /* DP matrix data */
  for (i = istart; i <= iend; i++)
  {
      fprintf(ofp, "%3d M ", i);
      for (k = kstart; k <= kend;        k++)  
	{
	  val = gx->dp[i][k * (p7G_NSCELLS + gx->allocC) + p7G_M + p7G_C0];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*.*f ", width, precision, val);
         
	  val = gx->dp[i][k * (p7G_NSCELLS + gx->allocC)  + p7G_M + p7G_C1];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*.*f ", width, precision, val);

	  val = gx->dp[i][k * (p7G_NSCELLS + gx->allocC) + p7G_M + p7G_C2];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*.*f ", width, precision, val);

       	  val = gx->dp[i][k * (p7G_NSCELLS + gx->allocC)  + p7G_M + p7G_C3];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*.*f ", width, precision, val);

	  val = gx->dp[i][k * (p7G_NSCELLS + gx->allocC)  + p7G_M + p7G_C4];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*.*f ", width, precision, val);

	  val = gx->dp[i][k * (p7G_NSCELLS + gx->allocC) + p7G_M + p7G_C5];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*.*f ", width, precision, val);
	}
      if (! (flags & p7_HIDE_SPECIALS))
	{
    	  for (x = 0;  x <  p7G_NXCELLS; x++) 
	    {
	      val = gx->xmx[i * p7G_NXCELLS + x];
	      if (flags & p7_SHOW_LOG) val = log(val);
	      fprintf(ofp, "%*.*f ", width, precision, val);
	    }
	}
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d I ", i);
      for (k = kstart; k <= kend;        k++) 
	{
	  val = gx->dp[i][k * (p7G_NSCELLS + gx->allocC) + p7G_I];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*.*f ", width, precision, val);
    	  fprintf(ofp, "%*.*f ", width, precision, 0.0);
	  fprintf(ofp, "%*.*f ", width, precision, 0.0);
	  fprintf(ofp, "%*.*f ", width, precision, 0.0);
	  fprintf(ofp, "%*.*f ", width, precision, 0.0);
	  fprintf(ofp, "%*.*f ", width, precision, 0.0);
	}
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d D ", i);
      for (k = kstart; k <= kend;        k++) 
	{
	  val =  gx->dp[i][k * (p7G_NSCELLS + gx->allocC)  + p7G_D];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*.*f ", width, precision, val);
	  fprintf(ofp, "%*.*f ", width, precision, 0.0);
	  fprintf(ofp, "%*.*f ", width, precision, 0.0);
	  fprintf(ofp, "%*.*f ", width, precision, 0.0);
	  fprintf(ofp, "%*.*f ", width, precision, 0.0);
	  fprintf(ofp, "%*.*f ", width, precision, 0.0);
        }
      fprintf(ofp, "\n\n");
  }
  return eslOK;
}

/* Function:  p7_gmx_fs_DumpWindow_Scientific()
 * Synopsis:  Dump a window of a frameshift DP matrix to a stream for diagnostics
 *            with values printed in scientific notaion.
 *
 * Purpose:   Dump a window of matrix <gx> to stream <fp> for diagnostics,
 *            from row <istart> to <iend>, from column <kstart> to <kend>.
 *            
 *            Asking for <0..L,0..M> with <flags=p7_SHOW_SPECIALS> is the
 *            same as <p7_gmx_Dump()>.
 *            
 *            <flags> control some optional output behaviors, as follows:
 *              | <p7_HIDE_SPECIALS> | don't show scores for <ENJBC> states  |
 *              | <p7_SHOW_LOG>      | <gx> is in probs; show as log probs   |
 *  
 * Returns:   <eslOK> on success.
 */
int
p7_gmx_fs_DumpWindow_Scientific(FILE *ofp, P7_GMX *gx, int istart, int iend, int kstart, int kend, int flags)
{
  int   width     = 9;
  int   i, k, c, x;
  float val;

  /* Header */
  fprintf(ofp, "     ");
  for (k = kstart; k <= kend;  k++) {
    fprintf(ofp, "%*d_0", width, k);
    for(c = 1; c <= gx->allocC; c++) 
      fprintf(ofp, "%*d_%d", width, k, c);
  }

  if (! (flags & p7_HIDE_SPECIALS)) fprintf(ofp, "%*s %*s %*s %*s %*s\n", width, "E", width, "N", width, "J", width, "B", width, "C");
  fprintf(ofp, "      ");
  for (k = kstart; k <= kend; k++) { 
    fprintf(ofp, "%*.*s ", width, width, "----------");
    for (c = 1; c <= gx->allocC; c++)  	
      fprintf(ofp, "%*.*s ", width, width, "----------");
  }

  if (! (flags & p7_HIDE_SPECIALS)) 
    for (x = 0; x < 5; x++) fprintf(ofp, "%*.*s ", width, width, "----------");
  fprintf(ofp, "\n");
  
  /* DP matrix data */
  for (i = istart; i <= iend; i++)
  {
      fprintf(ofp, "%3d M ", i);
      for (k = kstart; k <= kend;        k++)  
	{
	  val = gx->dp[i][k * (p7G_NSCELLS + gx->allocC) + p7G_M + p7G_C0];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*g ", width, val);
         
	  val = gx->dp[i][k * (p7G_NSCELLS + gx->allocC)  + p7G_M + p7G_C1];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*g ", width, val);

	  val = gx->dp[i][k * (p7G_NSCELLS + gx->allocC) + p7G_M + p7G_C2];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*g ", width, val);

       	  val = gx->dp[i][k * (p7G_NSCELLS + gx->allocC)  + p7G_M + p7G_C3];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*g ", width, val);

	  val = gx->dp[i][k * (p7G_NSCELLS + gx->allocC)  + p7G_M + p7G_C4];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*g ", width, val);

	  val = gx->dp[i][k * (p7G_NSCELLS + gx->allocC) + p7G_M + p7G_C5];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*g ", width, val);
	}
      if (! (flags & p7_HIDE_SPECIALS))
	{
    	  for (x = 0;  x <  p7G_NXCELLS; x++) 
	    {
	      val = gx->xmx[i * p7G_NXCELLS + x];
	      if (flags & p7_SHOW_LOG) val = log(val);
	      fprintf(ofp, "%*g ", width, val);
	    }
	}
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d I ", i);
      for (k = kstart; k <= kend;        k++) 
	{
	  val = gx->dp[i][k * (p7G_NSCELLS + gx->allocC) + p7G_I];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*g ", width, val);
    	  fprintf(ofp, "%*g ", width, 0.0);
	  fprintf(ofp, "%*g ", width, 0.0);
	  fprintf(ofp, "%*g ", width, 0.0);
	  fprintf(ofp, "%*g ", width, 0.0);
	  fprintf(ofp, "%*g ", width, 0.0);
	}
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d D ", i);
      for (k = kstart; k <= kend;        k++) 
	{
	  val =  gx->dp[i][k * (p7G_NSCELLS + gx->allocC)  + p7G_D];
	  if (flags & p7_SHOW_LOG) val = log(val);
	  fprintf(ofp, "%*g ", width, val);
	  fprintf(ofp, "%*g ", width, 0.0);
	  fprintf(ofp, "%*g ", width, 0.0);
	  fprintf(ofp, "%*g ", width, 0.0);
	  fprintf(ofp, "%*g ", width, 0.0);
	  fprintf(ofp, "%*g ", width, 0.0);
        }
      fprintf(ofp, "\n\n");
  }
  return eslOK;
}


/* Function:  p7_gmx_fs_ParserDumpw()
 * Synopsis:  Dump a single row of a frameshift DP matrix from a parser function.
 *
 * Purpose:   Dump a single row of matrix <gx> to stream <fp> for diagnostics,
 *            from column <kstart> to <kend>.
 *            
 *            Asking for <i, curr,0..M> with <flags=p7_SHOW_SPECIALS>
 *            
 *            <flags> control some optional output behaviors, as follows:
 *              | <p7_HIDE_SPECIALS> | don't show scores for <ENJBC> states  |
 *              | <p7_SHOW_LOG>      | <gx> is in probs; show as log probs   |
 *  
 * Returns:   <eslOK> on success.
 */
int
p7_gmx_fs_ParserDump(FILE *ofp, P7_GMX *gx, int i, int curr, int kstart, int kend, int flags)
{
  int   width     = 9;
  int   precision = 4;
  int   k, c, x;
  float val;

  /* Header */
  if(i == 0) {
    fprintf(ofp, "     ");
    for (k = kstart; k <= kend;  k++) {
      fprintf(ofp, "%*d,", width, k);
    }

    if (! (flags & p7_HIDE_SPECIALS)) fprintf(ofp, "%*s, %*s, %*s, %*s, %*s,\n", width, "E", width, "N", width, "J", width, "B", width, "C");
    fprintf(ofp, "      ");
    for (k = kstart; k <= kend; k++) { 
      fprintf(ofp, "%*.*s, ", width, width, "----------");
      for (c = 1; c <= gx->allocC; c++)  	
        fprintf(ofp, "%*.*s, ", width, width, "----------");
    }

    if (! (flags & p7_HIDE_SPECIALS)) 
      for (x = 0; x < 5; x++) fprintf(ofp, "%*.*s, ", width, width, "----------");
    fprintf(ofp, "\n");
  }

  /* DP matrix data */
  
  fprintf(ofp, "%3d, M, ", i);
      for (k = kstart; k <= kend;        k++)  
      {
	    val = gx->dp[curr][k * (p7G_NSCELLS + gx->allocC) + p7G_M + p7G_C0];
	    if (flags & p7_SHOW_LOG) val = log(val);
	    fprintf(ofp, "%*.*f, ", width, precision, val);
	  }
      if (! (flags & p7_HIDE_SPECIALS))
	  {
    	  for (x = 0;  x <  p7G_NXCELLS; x++) 
	    {
	      val = gx->xmx[i * p7G_NXCELLS + x];
	      if (flags & p7_SHOW_LOG) val = log(val);
	      fprintf(ofp, "%*.*f, ", width, precision, val);
	    }
	  }
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d, I, ", i);
      for (k = kstart; k <= kend;        k++) 
	  {
	    val = gx->dp[curr][k * (p7G_NSCELLS + gx->allocC) + p7G_I];
	    if (flags & p7_SHOW_LOG) val = log(val);
	    fprintf(ofp, "%*.*f, ", width, precision, val);
	  }
      fprintf(ofp, "\n");

      fprintf(ofp, "%3d, D, ", i);
      for (k = kstart; k <= kend;        k++) 
	  {
	    val =  gx->dp[curr][k * (p7G_NSCELLS + gx->allocC)  + p7G_D];
	    if (flags & p7_SHOW_LOG) val = log(val);
	    fprintf(ofp, "%*.*f, ", width, precision, val);
      }
      fprintf(ofp, "\n\n");
  
  return eslOK;
}



/*****************************************************************
 * 4. Unit tests
 *****************************************************************/
#ifdef p7GMX_FS_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

static void
gmx_testpattern(P7_GMX *gx, int M, int Lx)
{
  int i,k,s,n, n2;

  /* Write a test pattern, via the dp[i] pointers */
  n = 0;
  for (i = 0; i <= Lx; i++)
    for (k = 0; k <= M; k++)
      for (s = 0; s < (p7G_NSCELLS+gx->allocC); s++)
	    gx->dp[i][k*(p7G_NSCELLS+gx->allocC)+s] = n++;

  /* Read it back, via the dp[i] pointers */
  n = 0;
  for (i = 0; i <= Lx; i++)
    for (k = 0; k <= M; k++)
      for (s = 0; s < (p7G_NSCELLS+gx->allocC); s++)
		{
		  if (gx->dp[i][k*(p7G_NSCELLS+gx->allocC)+s] != n) esl_fatal("gmx_fs unit test failed: test pattern corrupted");
		  n++;
		}
   
  /* Reading it back via the dp_mem vector itself ought to be the same */
  if (gx->allocR == gx->validR && gx->ncells == (int64_t) gx->validR * (int64_t) gx->allocW) {
    n2 = 0;
    for (i = 0; i < gx->ncells; i++)
      for (s = 0; s < (p7G_NSCELLS+gx->allocC); s++)
	  {
	    if (gx->dp_mem[i*(p7G_NSCELLS+gx->allocC)+s] != n2) esl_fatal("gmx_fs unit test failed: test pattern corrupted (2nd test)");
		  n2++;
	  }  
  
    /* and the number of cells ought to match too */
    if (n != n2) esl_fatal("gmx_fs unit test failed: unexpected # of cells");

  }
}


static void
utest_GrowTo(void)
{
  int     M, L, Lx, C;
  P7_GMX *gx = NULL;
  M = 20;  L = 20;  Lx = 4;   C = 0; gx= p7_gmx_fs_Create(M, L, Lx, C);  gmx_testpattern(gx, M, Lx);
  M = 40;  L = 20;  Lx = 4;   C = 0; p7_gmx_fs_GrowTo(gx, M, L, Lx, C);  gmx_testpattern(gx, M, Lx);  /* grow in M only */
  M = 40;  L = 40;  Lx = 4;   C = 0; p7_gmx_fs_GrowTo(gx, M, L, Lx, C);  gmx_testpattern(gx, M, Lx);  /* grow in L only */
  M = 40;  L = 40;  Lx = 40;  C = 0; p7_gmx_fs_GrowTo(gx, M, L, Lx, C);  gmx_testpattern(gx, M, Lx);  /* grow in Lx only */
  M = 40;  L = 40;  Lx = 40;  C = 3; p7_gmx_fs_GrowTo(gx, M, L, Lx, C);  gmx_testpattern(gx, M, Lx);  /* grow in C only */
  M = 80;  L = 10;  Lx = 10;  C = 3; p7_gmx_fs_GrowTo(gx, M, L, Lx, C);  gmx_testpattern(gx, M, Lx);  /* grow in M, but with enough ncells */
  M = 10;  L = 80;  Lx = 80;  C = 3; p7_gmx_fs_GrowTo(gx, M, L, Lx, C);  gmx_testpattern(gx, M, Lx);  /* grow in Lx, but with enough ncells */
  M = 10;  L = 40;  Lx = 40;  C = 5; p7_gmx_fs_GrowTo(gx, M, L, Lx, C);  gmx_testpattern(gx, M, Lx);  /* grow in C, but with enough ncells */
  M = 100; L = 100; Lx = 100; C = 6; p7_gmx_fs_GrowTo(gx, M, L, Lx, C);  gmx_testpattern(gx, M, Lx);  /* grow in all */
  M = 100; L = 100; Lx = 100; C = 3; p7_gmx_fs_GrowTo(gx, M, L, Lx, C);  gmx_testpattern(gx, M, Lx);  /* shrink in C*/
  
  p7_gmx_Destroy(gx);
}

static void
utest_Compare(ESL_RANDOMNESS *r, P7_FS_PROFILE *gm_fs, P7_BG *bgDNA, ESL_GENCODE *gcode, int L, float tolerance)
{
  char         *msg = "gmx_Compare unit test failure";
  ESL_DSQ      *dsq = malloc(sizeof(ESL_DSQ) *(L+2));
  P7_GMX       *gx1 = p7_gmx_fs_Create(gm_fs->M, L, L, p7P_5CODONS);
  P7_GMX       *gx2 = p7_gmx_fs_Create(5, 4, 4, 0);
  P7_IVX       *iv1 = p7_ivx_Create(gm_fs->M, p7P_5CODONS);
  P7_IVX       *iv2 = p7_ivx_Create(5, p7P_3CODONS);
  float         fsc;

  if (!r || !gm_fs || !dsq || !gx1 || !gx2 )                                esl_fatal(msg);
  if (esl_rsq_xfIID(r, bgDNA->f, gcode->nt_abc->K, L, dsq)        != eslOK) esl_fatal(msg);
  if (p7_gmx_fs_GrowTo(gx2, gm_fs->M, L, L, p7P_5CODONS)          != eslOK) esl_fatal(msg);
  if (p7_Forward_Frameshift(dsq, gcode, L, gm_fs, gx1, iv1, &fsc) != eslOK) esl_fatal(msg);
  if (p7_Forward_Frameshift(dsq, gcode, L, gm_fs, gx2, iv1, &fsc) != eslOK) esl_fatal(msg);
  if (p7_gmx_fs_Compare(gx1, gx2, p7P_5CODONS, tolerance)         != eslOK) esl_fatal(msg);   

  p7_gmx_Reuse(gx2);
  if (p7_ivx_GrowTo(iv2, gm_fs->M, p7P_5CODONS)                   != eslOK) esl_fatal(msg);
  if (p7_Forward_Frameshift(dsq, gcode, L, gm_fs, gx2, iv2, &fsc) != eslOK) esl_fatal(msg);  
  if (p7_gmx_fs_Compare(gx1, gx2, p7P_5CODONS, tolerance)         != eslOK) esl_fatal(msg);

  p7_gmx_Destroy(gx1);
  p7_gmx_Destroy(gx2);
  p7_ivx_Destroy(iv1);
  p7_ivx_Destroy(iv2);
  free(dsq);
}

#endif /*p7GMX_FS_TESTDRIVE*/
/*------------------- end, unit tests ---------------------------*/


/*****************************************************************
 * 5. Test driver
 *****************************************************************/
#ifdef p7GMX_FS_TESTDRIVE
/*
  gcc -o p7_gmx_fs_utest -msse2 -g -Wall -I. -L. -I../easel -L../easel -Dp7GMX_FS_TESTDRIVE p7_gmx_fs.c -lhmmer -leasel -lm
  ./p7_gmx_fs_utest
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  { "-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",                  0},
  { "-s",  eslARG_INT,     "42",  NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",        0 },
  { "-t",  eslARG_REAL,  "1e-5",  NULL, NULL, NULL, NULL, NULL, "floating point comparison tolerance",  0 },
  { "-L",  eslARG_INT,    "120",  NULL, NULL, NULL, NULL, NULL, "length of sampled sequences",          0 },
  { "-M",  eslARG_INT,     "40",  NULL, NULL, NULL, NULL, NULL, "length of sampled test profile",       0 },
  { 0,0,0,0,0,0,0,0,0,0},
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for p7_gmx.c";

int 
main(int argc, char **argv)
{
  char           *msg    = "p7_gmx unit test driver failed";
  ESL_GETOPTS    *go     = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r      = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abcAA  = esl_alphabet_Create(eslAMINO);
  ESL_ALPHABET   *abcDNA = esl_alphabet_Create(eslDNA);
  ESL_GENCODE    *gcode  = esl_gencode_Create(abcDNA, abcAA);
  P7_BG          *bgAA   = p7_bg_Create(abcAA);
  P7_BG          *bgDNA  = p7_bg_Create(abcDNA);
  P7_HMM         *hmm   = NULL;
  P7_FS_PROFILE  *gm_fs = NULL;
  int             M     = esl_opt_GetInteger(go, "-M");
  int             L     = esl_opt_GetInteger(go, "-L");
  float           tol   = esl_opt_GetReal   (go, "-t");

  p7_FLogsumInit();

  if (p7_hmm_Sample(r, M, abcAA, &hmm)                              != eslOK) esl_fatal(msg);
  if ((gm_fs = p7_profile_fs_Create(hmm->M, abcAA))                 == NULL)  esl_fatal(msg);
  if (p7_bg_SetLength(bgAA, L/3)                                    != eslOK) esl_fatal(msg);
  if (p7_bg_SetLength(bgDNA, L)                                     != eslOK) esl_fatal(msg);
  if (p7_ProfileConfig_fs(hmm, bgAA, gcode, gm_fs, L/3, p7_UNILOCAL) != eslOK) esl_fatal(msg);

  utest_GrowTo();
  utest_Compare(r, gm_fs, bgDNA, gcode, L, tol);

  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  esl_gencode_Destroy(gcode);
  esl_alphabet_Destroy(abcAA);
  esl_alphabet_Destroy(abcDNA);
  p7_bg_Destroy(bgAA);
  p7_bg_Destroy(bgDNA);
  p7_hmm_Destroy(hmm);
  p7_profile_fs_Destroy(gm_fs);
  return eslOK;
}
#endif /*p7GMX_FS_TESTDRIVE*/
/*------------------ end, test driver ---------------------------*/




