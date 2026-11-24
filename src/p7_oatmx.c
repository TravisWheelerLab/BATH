/* P7_OATMX : what the optimal accuracy trace reads
 *
 * Contents:
 *   1. The <P7_OATMX> object
 *   2. Unit tests
 *   3. Test driver
 */
#include "p7_config.h"

#include <string.h>

#include "easel.h"

#include "hmmer.h"


/*****************************************************************
 * 1. The <P7_OATMX> object
 *****************************************************************/

/* oat_rowbytes(): bytes per plane row for a model of <M> positions. */
static int
oat_rowbytes(int M, int lanes)
{
  return (p7OAT_NQ(M, lanes) * lanes + 7) / 8;
}


/* Function:  p7_oatmx_Create()
 * Synopsis:  Allocate a <P7_OATMX>.
 *
 * Purpose:   Allocate a <P7_OATMX> for models up to <allocM> positions
 *            and sequences up to length <allocL>, with <nplanes> bit
 *            planes per row (<p7OAT_NPLANES>, or <p7OAT_NPLANES_FS> for
 *            the frameshift model) laid out for <lanes>: 1 for the
 *            generic code (cell k is bit k-1), else the number of floats
 *            per SIMD vector (striped layout).
 *
 * Returns:   a pointer to the new <P7_OATMX>.
 *
 * Throws:    <NULL> on allocation failure or bad arguments.
 */
P7_OATMX *
p7_oatmx_Create(int allocM, int allocL, int nplanes, int lanes)
{
  P7_OATMX *ox = NULL;
  int       status;

  if (allocM < 1 || allocL < 0 || nplanes < 1 || nplanes > p7OAT_NPLANES_FS || lanes < 1)
    ESL_XEXCEPTION(eslEINVAL, "bad arguments to p7_oatmx_Create()");

  ESL_ALLOC(ox, sizeof(P7_OATMX));
  ox->mem     = NULL;
  ox->xmx     = NULL;
  ox->ek      = NULL;

  ox->M        = 0;
  ox->L        = 0;
  ox->nplanes  = nplanes;
  ox->lanes    = lanes;
  ox->rowbytes = oat_rowbytes(allocM, lanes);

  ox->allocB   = (size_t) (allocL+1) * nplanes * ox->rowbytes;
  ox->allocXR  = allocL+1;
  ESL_ALLOC(ox->mem, sizeof(uint8_t) * ox->allocB);
  ESL_ALLOC(ox->xmx, sizeof(float)   * ox->allocXR * p7OAT_NXCELLS);
  ESL_ALLOC(ox->ek,  sizeof(int)     * ox->allocXR);
  return ox;

 ERROR:
  p7_oatmx_Destroy(ox);
  return NULL;
}


/* Function:  p7_oatmx_GrowTo()
 * Synopsis:  Make a <P7_OATMX> big enough for <M> by <L>.
 *
 * Purpose:   Reallocate <ox> if needed for a model of <M> positions and
 *            a sequence of length <L>, and set <ox->M>, <ox->L> and the
 *            plane row size for them. Contents are not kept.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; <eslEINVAL> on bad
 *            arguments.
 */
int
p7_oatmx_GrowTo(P7_OATMX *ox, int M, int L)
{
  int    rowbytes = oat_rowbytes(M, ox->lanes);
  size_t nbytes   = (size_t) (L+1) * ox->nplanes * rowbytes;
  int    status;

  if (M < 1 || L < 0) ESL_EXCEPTION(eslEINVAL, "bad arguments to p7_oatmx_GrowTo()");

  if (nbytes > ox->allocB)
    {
      ESL_REALLOC(ox->mem, sizeof(uint8_t) * nbytes);
      ox->allocB = nbytes;
    }
  if (L+1 > ox->allocXR)
    {
      ESL_REALLOC(ox->xmx, sizeof(float) * (L+1) * p7OAT_NXCELLS);
      ESL_REALLOC(ox->ek,  sizeof(int)   * (L+1));
      ox->allocXR = L+1;
    }
  ox->M        = M;
  ox->L        = L;
  ox->rowbytes = rowbytes;
  return eslOK;

 ERROR:
  return status;
}


/* Function:  p7_oatmx_Sizeof()
 * Synopsis:  Bytes allocated for a <P7_OATMX>.
 */
size_t
p7_oatmx_Sizeof(const P7_OATMX *ox)
{
  return sizeof(P7_OATMX) + ox->allocB + (size_t) ox->allocXR * (p7OAT_NXCELLS * sizeof(float) + sizeof(int));
}


/* Function:  p7_oatmx_Destroy()
 * Synopsis:  Free a <P7_OATMX>.
 */
void
p7_oatmx_Destroy(P7_OATMX *ox)
{
  if (ox == NULL) return;
  free(ox->mem);
  free(ox->xmx);
  free(ox->ek);
  free(ox);
}


/*****************************************************************
 * 2. Unit tests
 *****************************************************************/
#ifdef p7OATMX_TESTDRIVE
#include "esl_random.h"

/* oat_testbit(): bit index of cell <k> for the layout of <ox>, written
 * out here independently of the OATX_* and OATO_* macros. */
static int
oat_testbit(const P7_OATMX *ox, int k)
{
  int Q;
  if (ox->lanes == 1) return k-1;
  Q = ESL_MAX(2, (ox->M - 1) / ox->lanes + 1);
  return ((k-1) % Q) * ox->lanes + (k-1) / Q;
}

/* utest_fields()
 * Write random values into every field of every cell, then read them
 * back, for each layout; the generic layout is also read through the
 * OATX_* macros. Setting a field must not change other fields or
 * other cells.
 */
static void
utest_fields(ESL_RANDOMNESS *r)
{
  char     *msg      = "P7_OATMX fields unit test failed";
  int       lanes[3] = { 1, 4, 8 };
  int       nbits[7] = { 2, 1, 1, 4, 4, 3 };
  int       p0[7]    = { p7OAT_MP, p7OAT_IP, p7OAT_DP, p7OAT_PPM, p7OAT_PPI, p7OAT_COD };
  int       Ms[4]    = { 1, 7, 33, 150 };
  int       Ls[3]    = { 0, 1, 40 };
  P7_OATMX *ox       = NULL;
  int      *val      = NULL;
  int       a, b, c, fs, nf, f, i, k, M, L;

  for (fs = 0; fs <= 1; fs++)
    for (a = 0; a < 3; a++)
      {
        nf = fs ? 6 : 5;
        if ((ox = p7_oatmx_Create(10, 2, fs ? p7OAT_NPLANES_FS : p7OAT_NPLANES, lanes[a])) == NULL) esl_fatal(msg);
        for (b = 0; b < 4; b++)
          for (c = 0; c < 3; c++)
            {
              M = Ms[b];
              L = Ls[c];
              if (p7_oatmx_GrowTo(ox, M, L) != eslOK) esl_fatal(msg);
              if (ox->rowbytes * 8 < p7OAT_NQ(M, ox->lanes) * ox->lanes) esl_fatal("%s: plane rows too short", msg);
              memset(ox->mem, 0xff, (size_t) (L+1) * ox->nplanes * ox->rowbytes);

              if ((val = realloc(val, sizeof(int) * (L+1) * (M+1) * nf)) == NULL) esl_fatal(msg);
              for (i = 0; i <= L; i++)
                for (k = 1; k <= M; k++)
                  for (f = 0; f < nf; f++)
                    {
                      val[(i*(M+1) + k)*nf + f] = esl_rnd_Roll(r, 1 << nbits[f]);
                      p7_oatmx_SetField(ox, i, oat_testbit(ox, k), p0[f], nbits[f], val[(i*(M+1) + k)*nf + f]);
                    }
              for (i = 0; i <= L; i++)
                for (k = 1; k <= M; k++)
                  for (f = 0; f < nf; f++)
                    if (p7_oatmx_GetField(ox, i, oat_testbit(ox, k), p0[f], nbits[f]) != val[(i*(M+1) + k)*nf + f])
                      esl_fatal("%s: lanes %d, M %d, row %d, k %d, field %d", msg, ox->lanes, M, i, k, f);

              if (ox->lanes == 1)
                for (i = 0; i <= L; i++)
                  for (k = 1; k <= M; k++)
                    {
                      int *v = val + (i*(M+1) + k)*nf;
                      if (OATX_MPRED(ox,i,k) != v[0] || OATX_IPRED(ox,i,k) != v[1] || OATX_DPRED(ox,i,k) != v[2] ||
                          OATX_PPM(ox,i,k)   != v[3] || OATX_PPI(ox,i,k)   != v[4] || (fs && OATX_CODON(ox,i,k) != v[5]))
                        esl_fatal("%s: OATX macros, M %d, row %d, k %d", msg, M, i, k);
                    }
            }

        /* the OATX setters, read back through p7_oatmx_GetField() */
        if (ox->lanes == 1)
          {
            if (p7_oatmx_GrowTo(ox, 37, 5) != eslOK) esl_fatal(msg);
            for (i = 0; i <= 5; i++)
              for (k = 1; k <= 37; k++)
                {
                  OATX_SET_MPRED(ox,i,k, (i+k) % 4);
                  OATX_SET_IPRED(ox,i,k, (i+k) % 2);
                  OATX_SET_DPRED(ox,i,k, (i*k) % 2);
                  OATX_SET_PPM  (ox,i,k, (i+2*k) % 11);
                  OATX_SET_PPI  (ox,i,k, (2*i+k) % 11);
                  if (fs) OATX_SET_CODON(ox,i,k, 1 + (i+k) % 5);
                }
            for (i = 0; i <= 5; i++)
              for (k = 1; k <= 37; k++)
                if (p7_oatmx_GetField(ox, i, k-1, p7OAT_MP,  2) != (i+k) % 4    ||
                    p7_oatmx_GetField(ox, i, k-1, p7OAT_IP,  1) != (i+k) % 2    ||
                    p7_oatmx_GetField(ox, i, k-1, p7OAT_DP,  1) != (i*k) % 2    ||
                    p7_oatmx_GetField(ox, i, k-1, p7OAT_PPM, 4) != (i+2*k) % 11 ||
                    p7_oatmx_GetField(ox, i, k-1, p7OAT_PPI, 4) != (2*i+k) % 11 ||
                    (fs && p7_oatmx_GetField(ox, i, k-1, p7OAT_COD, 3) != 1 + (i+k) % 5))
                  esl_fatal("%s: OATX setters, row %d, k %d", msg, i, k);
          }
        p7_oatmx_Destroy(ox);
      }
  free(val);
}

/* utest_ppcodes()
 * Posterior codes print as p7_alidisplay_EncodePostProb() prints the
 * probability, and the probability for each code prints as that code.
 */
static void
utest_ppcodes(ESL_RANDOMNESS *r)
{
  char  *msg = "P7_OATMX posterior code unit test failed";
  char   sym[11] = { '0','1','2','3','4','5','6','7','8','9','*' };
  float  p;
  int    c, n;

  for (c = 0; c <= 10; c++)
    if (p7_alidisplay_EncodePostProb(p7OAT_PPVALUE(c)) != sym[c]) esl_fatal("%s: value of code %d", msg, c);

  for (n = 0; n < 100000; n++)
    {
      p = (n < 1001) ? (float) n / 1000.0f : esl_random(r);
      if (p7_alidisplay_EncodePostProb(p) != sym[p7OAT_PPCODE(p)]) esl_fatal("%s: code of %f", msg, p);
    }
  if (p7OAT_PPCODE(-1e-7f) != 0) esl_fatal("%s: small negative probability", msg);
}

/* utest_grow()
 * Growing keeps the plane row size right for each model size, and
 * Sizeof() covers what is allocated.
 */
static void
utest_grow(void)
{
  char     *msg = "P7_OATMX grow unit test failed";
  P7_OATMX *ox  = p7_oatmx_Create(5, 3, p7OAT_NPLANES, 8);

  if (ox == NULL) esl_fatal(msg);
  if (p7_oatmx_GrowTo(ox, 300, 1000) != eslOK) esl_fatal(msg);
  if (ox->rowbytes != p7OAT_NQ(300, 8))                                    esl_fatal("%s: rowbytes", msg);
  if (ox->allocB   <  (size_t) 1001 * p7OAT_NPLANES * ox->rowbytes)       esl_fatal("%s: allocB", msg);
  if (ox->allocXR  <  1001)                                                esl_fatal("%s: allocXR", msg);
  if (p7_oatmx_Sizeof(ox) < ox->allocB)                                    esl_fatal("%s: Sizeof", msg);
  if (p7_oatmx_GrowTo(ox, 3, 2) != eslOK)                                  esl_fatal(msg);
  if (ox->rowbytes != 2 || ox->M != 3 || ox->L != 2)                      esl_fatal("%s: shrink", msg);
  p7_oatmx_Destroy(ox);
}
#endif /*p7OATMX_TESTDRIVE*/


/*****************************************************************
 * 3. Test driver
 *****************************************************************/
#ifdef p7OATMX_TESTDRIVE
#include "esl_getopts.h"

static ESL_OPTIONS options[] = {
  /* name  type         default  env   range togs  reqs  incomp  help                           docgroup */
  { "-h",  eslARG_NONE,   FALSE, NULL, NULL, NULL, NULL, NULL, "show brief help on version and usage", 0 },
  { "-s",  eslARG_INT,     "42", NULL, NULL, NULL, NULL, NULL, "set random number seed to <n>",        0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options]";
static char banner[] = "test driver for the P7_OATMX optimal accuracy trace matrix";

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go = p7_CreateDefaultApp(options, 0, argc, argv, banner, usage);
  ESL_RANDOMNESS *r  = esl_randomness_CreateFast(esl_opt_GetInteger(go, "-s"));

  utest_fields(r);
  utest_ppcodes(r);
  utest_grow();

  esl_randomness_Destroy(r);
  esl_getopts_Destroy(go);
  return 0;
}
#endif /*p7OATMX_TESTDRIVE*/
