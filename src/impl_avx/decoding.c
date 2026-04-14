/* Posterior decoding dispatcher.
 * Provides the non-suffixed p7_Decoding() and p7_DomainDecoding() declared
 * in impl_avx.h.  Delegates to the fastest available ISA implementation
 * at runtime.
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/* Function:  p7_Decoding()
 *
 * Purpose:   Dispatch posterior decoding of residue assignment to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> if numeric range is exceeded.
 */
int
p7_Decoding(const P7_OPROFILE *om, const P7_OMX *oxf, P7_OMX *oxb, P7_OMX *pp)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_Decoding_avx512(om, oxf, oxb, pp);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_Decoding_sse(om, oxf, oxb, pp);
#endif
#ifdef eslENABLE_SSE
  return p7_Decoding_sse(om, oxf, oxb, pp);
#else
  return eslENORESULT;
#endif
}


/* Function:  p7_DomainDecoding()
 *
 * Purpose:   Dispatch posterior decoding of domain location to the fastest ISA.
 *
 * Returns:   <eslOK> on success.
 *            <eslERANGE> on numeric overflow.
 */
int
p7_DomainDecoding(const P7_OPROFILE *om, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) return p7_DomainDecoding_avx512(om, oxf, oxb, ddef);
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    return p7_DomainDecoding_sse(om, oxf, oxb, ddef);
#endif
#ifdef eslENABLE_SSE
  return p7_DomainDecoding_sse(om, oxf, oxb, ddef);
#else
  return eslENORESULT;
#endif
}
/*****************************************************************
 * 3. Unit tests
 *****************************************************************/
#ifdef p7DECODING_TESTDRIVE
#include "esl_random.h"
#include "esl_randomseq.h"

/* compare results to GDecoding(). */
static void
utest_decoding(ESL_RANDOMNESS *r, ESL_ALPHABET *abc, P7_BG *bg, int M, int L, int N, float tolerance)
{
  char        *msg  = "decoding unit test failed";
  P7_HMM      *hmm  = NULL;
  P7_PROFILE  *gm   = NULL;
  P7_OPROFILE *om   = NULL;
  ESL_DSQ     *dsq  = malloc(sizeof(ESL_DSQ) * (L+2));
  P7_OMX      *fwd  = p7_omx_Create(M, L, L);
  P7_OMX      *bck  = p7_omx_Create(M, L, L);
  P7_OMX      *pp   = p7_omx_Create(M, L, L);
  P7_GMX      *gxf  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX      *gxb  = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX      *gxp1 = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  P7_GMX      *gxp2 = p7_gmx_Create(M, L, L, p7G_NSCELLS);
  float fsc1, fsc2;
  float bsc1, bsc2;

  if (p7_oprofile_Sample(r, abc, bg, M, L, &hmm, &gm, &om) != eslOK) esl_fatal(msg);
  while (N--)
    {
      if (esl_rsq_xfIID(r, bg->f, abc->K, L, dsq) != eslOK) esl_fatal(msg);
      if (p7_Forward       (dsq, L, om, fwd,      &fsc1) != eslOK) esl_fatal(msg);
      if (p7_Backward      (dsq, L, om, fwd, bck, &bsc1) != eslOK) esl_fatal(msg);
      if (p7_Decoding(om, fwd, bck, pp)                  != eslOK) esl_fatal(msg);
      if (p7_omx_FDeconvert(pp, gxp1)                    != eslOK) esl_fatal(msg);
      
      if (p7_GForward (dsq, L, gm, gxf, &fsc2)           != eslOK) esl_fatal(msg);
      if (p7_GBackward(dsq, L, gm, gxb, &bsc2)           != eslOK) esl_fatal(msg);
      if (p7_GDecoding(gm, gxf, gxb, gxp2)               != eslOK) esl_fatal(msg);

      // p7_gmx_Dump(stdout, gxp1, p7_DEFAULT);
      // p7_gmx_Dump(stdout, gxp2, p7_DEFAULT);

      if (p7_gmx_Compare(gxp1, gxp2, tolerance)          != eslOK) esl_fatal(msg);
    }

  p7_gmx_Destroy(gxp1);
  p7_gmx_Destroy(gxp2);
  p7_gmx_Destroy(gxf);
  p7_gmx_Destroy(gxb);
  p7_omx_Destroy(fwd);
  p7_omx_Destroy(bck);
  p7_omx_Destroy(pp);
  free(dsq);
  p7_oprofile_Destroy(om);
  p7_profile_Destroy(gm);
  p7_hmm_Destroy(hmm);
}
#endif /*p7DECODING_TESTDRIVE*/
/*--------------------- end, unit tests -------------------------*/




/*****************************************************************
 * 4. Test driver
 *****************************************************************/
#ifdef p7DECODING_TESTDRIVE

/*
  gcc -o decoding_utest -msse2 -g -Wall -I.. -L.. -I../../easel -L../../easel -Dp7DECODING_TESTDRIVE decoding.c -lhmmer -leasel -lm
  ./decoding_utest
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include "easel.h"
#include "esl_getopts.h"
#include "esl_random.h"
#include "esl_randomseq.h"

#include "hmmer.h"
#include "impl_avx.h"

static ESL_OPTIONS options[] = {
  /* name  type         default  env   range togs  reqs  incomp  help                docgrp */
  { "-h",  eslARG_NONE,    FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",                  0},
  { "-s",  eslARG_INT,     "42",  NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",        0 },
  { "-t",  eslARG_REAL,  "0.01",  NULL, NULL, NULL, NULL, NULL, "floating point comparison tolerance",  0 },
  { "-L",  eslARG_INT,     "40",  NULL, NULL, NULL, NULL, NULL, "length of sampled sequences",          0 },
  { "-M",  eslARG_INT,     "40",  NULL, NULL, NULL, NULL, NULL, "length of sampled test profile",       0 },
  { "-N",  eslARG_INT,     "10",  NULL, NULL, NULL, NULL, NULL, "number of sampled test sequences",     0 },
  { 0,0,0,0,0,0,0,0,0,0},
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for SSE posterior decoding";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go   = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r    = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));
  ESL_ALPHABET   *abc  = esl_alphabet_Create(eslAMINO);
  P7_BG          *bg   = p7_bg_Create(abc);
  int             M    = esl_opt_GetInteger(go, "-M");
  int             L    = esl_opt_GetInteger(go, "-L");
  int             N    = esl_opt_GetInteger(go, "-N");
  float           tol  = esl_opt_GetReal   (go, "-t");
  
  p7_FLogsumInit();

  utest_decoding(r, abc, bg, M, L, N, tol);
  
  esl_getopts_Destroy(go);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  p7_bg_Destroy(bg);
  return eslOK;
}
#endif /*p7DECODING_TESTDRIVE*/
/*-------------------- end, test driver -------------------------*/


