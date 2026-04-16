/* p7_oprofile.c — runtime dispatch for P7_OPROFILE management.
 *
 * Contains:
 *   1. Function pointer variable definitions (one definition each; extern in impl_avx.h).
 *   2. impl_Init(): CPU detection, flush-zero mode, function pointer assignment.
 *   3. ISA-independent functions (IsLocal, Dump, Sample, Compare, etc.).
 *
 * ISA-specific implementations live in:
 *   p7_oprofile_sse.c    — SSE (128-bit) versions, suffix _sse
 *   p7_oprofile_avx.c    — AVX-256 versions, suffix _avx      [TODO]
 *   p7_oprofile_avx512.c — AVX-512 versions, suffix _avx512   [TODO]
 */
#include "p7_config.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#ifdef eslENABLE_SSE
#include <xmmintrin.h>
#include <emmintrin.h>
#ifdef __SSE3__
#include <pmmintrin.h>
#endif
#endif

#include "easel.h"
#include "esl_cpu.h"
#include "esl_random.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "impl_avx.h"


/*****************************************************************
 * 1. Function pointer variable definitions.
 *
 * Each is initialized to NULL; impl_Init() assigns them to the
 * fastest available ISA implementation.  All are declared extern
 * in impl_avx.h.
 *****************************************************************/

/* P7_OPROFILE lifecycle */
P7_OPROFILE *(*p7_oprofile_Create)            (int M, const ESL_ALPHABET *abc)                                                    = NULL;
void         (*p7_oprofile_Destroy)           (P7_OPROFILE *om)                                                                   = NULL;
P7_OPROFILE *(*p7_oprofile_Clone)             (const P7_OPROFILE *om)                                                             = NULL;
int          (*p7_oprofile_Convert)           (const P7_PROFILE *gm, P7_OPROFILE *om)                                            = NULL;
int          (*p7_oprofile_Convert_Log)       (const P7_PROFILE *gm, P7_OPROFILE *om)                                            = NULL;
int          (*p7_oprofile_ReconfigLength)    (P7_OPROFILE *om, int L)                                                            = NULL;
int          (*p7_oprofile_ReconfigLength_Log)(P7_OPROFILE *om, int L)                                                            = NULL;
int          (*p7_oprofile_ReconfigMSVLength) (P7_OPROFILE *om, int L)                                                            = NULL;
int          (*p7_oprofile_ReconfigRestLength)(P7_OPROFILE *om, int L)                                                            = NULL;
int          (*p7_oprofile_ReconfigMultihit)  (P7_OPROFILE *om, int L)                                                            = NULL;
int          (*p7_oprofile_ReconfigMultihit_Log)(P7_OPROFILE *om, int L)                                                          = NULL;
int          (*p7_oprofile_ReconfigUnihit)    (P7_OPROFILE *om, int L)                                                            = NULL;
int          (*p7_oprofile_ReconfigUnihit_Log)(P7_OPROFILE *om, int L)                                                            = NULL;
int          (*p7_oprofile_UpdateFwdEmissionScores)(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr)              = NULL;
int          (*p7_oprofile_UpdateVitEmissionScores)(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr)              = NULL;
int          (*p7_oprofile_UpdateMSVEmissionScores)(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr)              = NULL;

/*****************************************************************
 * 2. impl_Init(): CPU detection and function pointer assignment.
 *
 * Call once at program startup before any profile operations.
 * Selects the fastest available ISA (AVX-512 > AVX-256 > SSE)
 * and wires all function pointers accordingly.
 *****************************************************************/
void
impl_Init(void)
{
  /* Set flush-zero mode on x86 to avoid sub-normal fp slowdowns */
#ifdef HAVE_FLUSH_ZERO_MODE
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif
#ifdef _PMMINTRIN_H_INCLUDED
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif

#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512())
    {
      /* P7_OPROFILE */
      p7_oprofile_Create                = p7_oprofile_Create_avx512;
      p7_oprofile_Destroy               = p7_oprofile_Destroy_avx512;
      p7_oprofile_Clone                 = p7_oprofile_Clone_avx512;
      p7_oprofile_Convert               = p7_oprofile_Convert_avx512;
      p7_oprofile_Convert_Log           = p7_oprofile_Convert_Log_avx512;
      p7_oprofile_ReconfigLength        = p7_oprofile_ReconfigLength_avx512;
      p7_oprofile_ReconfigLength_Log    = p7_oprofile_ReconfigLength_Log_avx512;
      p7_oprofile_ReconfigMSVLength     = p7_oprofile_ReconfigMSVLength_avx512;
      p7_oprofile_ReconfigRestLength    = p7_oprofile_ReconfigRestLength_avx512;
      p7_oprofile_ReconfigMultihit      = p7_oprofile_ReconfigMultihit_avx512;
      p7_oprofile_ReconfigMultihit_Log  = p7_oprofile_ReconfigMultihit_Log_avx512;
      p7_oprofile_ReconfigUnihit        = p7_oprofile_ReconfigUnihit_avx512;
      p7_oprofile_ReconfigUnihit_Log    = p7_oprofile_ReconfigUnihit_Log_avx512;
      p7_oprofile_UpdateFwdEmissionScores = p7_oprofile_UpdateFwdEmissionScores_avx512;
      p7_oprofile_UpdateVitEmissionScores = p7_oprofile_UpdateVitEmissionScores_avx512;
      p7_oprofile_UpdateMSVEmissionScores = p7_oprofile_UpdateMSVEmissionScores_avx512;
      /* P7_FS_OPROFILE */
      p7_fs_oprofile_Create             = p7_fs_oprofile_Create_avx512;
      p7_fs_oprofile_Destroy            = p7_fs_oprofile_Destroy_avx512;
      p7_fs_oprofile_Clone              = p7_fs_oprofile_Clone_avx512;
      p7_fs_oprofile_Convert            = p7_fs_oprofile_Convert_avx512;
      p7_fs_oprofile_Convert_Log        = p7_fs_oprofile_Convert_Log_avx512;
      p7_fs_oprofile_SubConvert_Log     = p7_fs_oprofile_SubConvert_Log_avx512;
      p7_fs_oprofile_ReconfigLength     = p7_fs_oprofile_ReconfigLength_avx512;
      p7_fs_oprofile_ReconfigLength_Log = p7_fs_oprofile_ReconfigLength_Log_avx512;
      p7_fs_oprofile_ReconfigMultihit   = p7_fs_oprofile_ReconfigMultihit_avx512;
      p7_fs_oprofile_ReconfigUnihit     = p7_fs_oprofile_ReconfigUnihit_avx512;
      p7_fs_oprofile_Logify             = p7_fs_oprofile_Logify_avx512;
      /* P7_OIVX */
      p7_oivx_Create                    = p7_oivx_Create_avx512;
      p7_oivx_GrowTo                    = p7_oivx_GrowTo_avx512;
      p7_oivx_Destroy                   = p7_oivx_Destroy_avx512;
      return;
    }
#endif

#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())
    {
      /* P7_OPROFILE — AVX2 implementations */
      p7_oprofile_Create                = p7_oprofile_Create_avx;
      p7_oprofile_Destroy               = p7_oprofile_Destroy_avx;
      p7_oprofile_Clone                 = p7_oprofile_Clone_avx;
      p7_oprofile_Convert               = p7_oprofile_Convert_avx;
      p7_oprofile_Convert_Log           = p7_oprofile_Convert_Log_avx;
      p7_oprofile_ReconfigLength        = p7_oprofile_ReconfigLength_avx;
      p7_oprofile_ReconfigLength_Log    = p7_oprofile_ReconfigLength_Log_avx;
      p7_oprofile_ReconfigMSVLength     = p7_oprofile_ReconfigMSVLength_avx;
      p7_oprofile_ReconfigRestLength    = p7_oprofile_ReconfigRestLength_avx;
      p7_oprofile_ReconfigMultihit      = p7_oprofile_ReconfigMultihit_avx;
      p7_oprofile_ReconfigMultihit_Log  = p7_oprofile_ReconfigMultihit_Log_avx;
      p7_oprofile_ReconfigUnihit        = p7_oprofile_ReconfigUnihit_avx;
      p7_oprofile_ReconfigUnihit_Log    = p7_oprofile_ReconfigUnihit_Log_avx;
      p7_oprofile_UpdateFwdEmissionScores = p7_oprofile_UpdateFwdEmissionScores_avx;
      p7_oprofile_UpdateVitEmissionScores = p7_oprofile_UpdateVitEmissionScores_avx;
      p7_oprofile_UpdateMSVEmissionScores = p7_oprofile_UpdateMSVEmissionScores_avx;
      /* P7_FS_OPROFILE — AVX2 implementations */
      p7_fs_oprofile_Create             = p7_fs_oprofile_Create_avx;
      p7_fs_oprofile_Destroy            = p7_fs_oprofile_Destroy_avx;
      p7_fs_oprofile_Clone              = p7_fs_oprofile_Clone_avx;
      p7_fs_oprofile_Convert            = p7_fs_oprofile_Convert_avx;
      p7_fs_oprofile_Convert_Log        = p7_fs_oprofile_Convert_Log_avx;
      p7_fs_oprofile_SubConvert_Log     = p7_fs_oprofile_SubConvert_Log_avx;
      p7_fs_oprofile_ReconfigLength     = p7_fs_oprofile_ReconfigLength_avx;
      p7_fs_oprofile_ReconfigLength_Log = p7_fs_oprofile_ReconfigLength_Log_avx;
      p7_fs_oprofile_ReconfigMultihit   = p7_fs_oprofile_ReconfigMultihit_avx;
      p7_fs_oprofile_ReconfigUnihit     = p7_fs_oprofile_ReconfigUnihit_avx;
      p7_fs_oprofile_Logify             = p7_fs_oprofile_Logify_avx;
      /* P7_OIVX — AVX2 implementations */
      p7_oivx_Create                    = p7_oivx_Create_avx;
      p7_oivx_GrowTo                    = p7_oivx_GrowTo_avx;
      p7_oivx_Destroy                   = p7_oivx_Destroy_avx;
      return;
    }
#endif

#ifdef eslENABLE_SSE
  /* Fallback: SSE (always available on x86 targets we support) */
  p7_oprofile_Create                = p7_oprofile_Create_sse;
  p7_oprofile_Destroy               = p7_oprofile_Destroy_sse;
  p7_oprofile_Clone                 = p7_oprofile_Clone_sse;
  p7_oprofile_Convert               = p7_oprofile_Convert_sse;
  p7_oprofile_Convert_Log           = p7_oprofile_Convert_Log_sse;
  p7_oprofile_ReconfigLength        = p7_oprofile_ReconfigLength_sse;
  p7_oprofile_ReconfigLength_Log    = p7_oprofile_ReconfigLength_Log_sse;
  p7_oprofile_ReconfigMSVLength     = p7_oprofile_ReconfigMSVLength_sse;
  p7_oprofile_ReconfigRestLength    = p7_oprofile_ReconfigRestLength_sse;
  p7_oprofile_ReconfigMultihit      = p7_oprofile_ReconfigMultihit_sse;
  p7_oprofile_ReconfigMultihit_Log  = p7_oprofile_ReconfigMultihit_Log_sse;
  p7_oprofile_ReconfigUnihit        = p7_oprofile_ReconfigUnihit_sse;
  p7_oprofile_ReconfigUnihit_Log    = p7_oprofile_ReconfigUnihit_Log_sse;
  p7_oprofile_UpdateFwdEmissionScores = p7_oprofile_UpdateFwdEmissionScores_sse;
  p7_oprofile_UpdateVitEmissionScores = p7_oprofile_UpdateVitEmissionScores_sse;
  p7_oprofile_UpdateMSVEmissionScores = p7_oprofile_UpdateMSVEmissionScores_sse;
  /* P7_FS_OPROFILE */
  p7_fs_oprofile_Create             = p7_fs_oprofile_Create_sse;
  p7_fs_oprofile_Destroy            = p7_fs_oprofile_Destroy_sse;
  p7_fs_oprofile_Clone              = p7_fs_oprofile_Clone_sse;
  p7_fs_oprofile_Convert            = p7_fs_oprofile_Convert_sse;
  p7_fs_oprofile_Convert_Log        = p7_fs_oprofile_Convert_Log_sse;
  p7_fs_oprofile_SubConvert_Log     = p7_fs_oprofile_SubConvert_Log_sse;
  p7_fs_oprofile_ReconfigLength     = p7_fs_oprofile_ReconfigLength_sse;
  p7_fs_oprofile_ReconfigLength_Log = p7_fs_oprofile_ReconfigLength_Log_sse;
  p7_fs_oprofile_ReconfigMultihit   = p7_fs_oprofile_ReconfigMultihit_sse;
  p7_fs_oprofile_ReconfigUnihit     = p7_fs_oprofile_ReconfigUnihit_sse;
  p7_fs_oprofile_Logify             = p7_fs_oprofile_Logify_sse;
  /* P7_OIVX */
  p7_oivx_Create                    = p7_oivx_Create_sse;
  p7_oivx_GrowTo                    = p7_oivx_GrowTo_sse;
  p7_oivx_Destroy                   = p7_oivx_Destroy_sse;
#endif
}


/*****************************************************************
 * 3. ISA-independent functions.
 *
 * These do not access any SIMD vector fields, so they work
 * regardless of which ISA is active.
 *****************************************************************/

/* Function:  p7_oprofile_IsLocal()
 * Synopsis:  Returns TRUE if profile is in local alignment mode.
 */
int
p7_oprofile_IsLocal(const P7_OPROFILE *om)
{
  if (om->mode == p7_LOCAL || om->mode == p7_UNILOCAL) return TRUE;
  return FALSE;
}

/* Function:  p7_oprofile_Dump()
 * Synopsis:  Dump a P7_OPROFILE to a stream for debugging.
 *
 * Note: Accesses scalar fields only; vector fields would need
 *       ISA-specific handling.  Add per-ISA vector dump later.
 */
int
p7_oprofile_Dump(FILE *fp, const P7_OPROFILE *om)
{
  fprintf(fp, "P7_OPROFILE: %s  M=%d  mode=%d\n", om->name ? om->name : "(unnamed)", om->M, om->mode);
  fprintf(fp, "  scale_b=%.4f  base_b=%u  bias_b=%u\n", om->scale_b, om->base_b, om->bias_b);
  fprintf(fp, "  scale_w=%.4f  base_w=%d  ddbound_w=%d\n", om->scale_w, om->base_w, om->ddbound_w);
  fprintf(fp, "  allocM=%d  allocQ4=%d  allocQ8=%d  allocQ16=%d\n",
          om->allocM, om->allocQ4, om->allocQ8, om->allocQ16);
#ifdef eslENABLE_AVX
  fprintf(fp, "  allocQ4_avx=%d  allocQ8_avx=%d  allocQ16_avx=%d\n",
          om->allocQ4_avx, om->allocQ8_avx, om->allocQ16_avx);
#endif
#ifdef eslENABLE_AVX512
  fprintf(fp, "  allocQ4_avx512=%d  allocQ8_avx512=%d  allocQ16_avx512=%d\n",
          om->allocQ4_avx512, om->allocQ8_avx512, om->allocQ16_avx512);
#endif
  return eslOK;
}

/* Function:  p7_oprofile_Sample()
 * Synopsis:  Sample a random optimized profile (for unit testing).
 */
int
p7_oprofile_Sample(ESL_RANDOMNESS *rng, const ESL_ALPHABET *abc, const P7_BG *bg,
                   int M, int L, P7_HMM **opt_hmm, P7_PROFILE **opt_gm, P7_OPROFILE **ret_om)
{
  P7_HMM     *hmm = NULL;
  P7_PROFILE *gm  = NULL;
  P7_OPROFILE *om = NULL;
  int          status;

  if ((status = p7_hmm_Sample(rng, M, abc, &hmm))        != eslOK) goto ERROR;
  if ((gm = p7_profile_Create(hmm->M, abc))              == NULL)  { status = eslEMEM; goto ERROR; }
  if ((status = p7_ProfileConfig(hmm, bg, gm, L, p7_LOCAL)) != eslOK) goto ERROR;
  if ((om = (*p7_oprofile_Create)(M, abc))               == NULL)  { status = eslEMEM; goto ERROR; }
  if ((status = (*p7_oprofile_Convert)(gm, om))          != eslOK) goto ERROR;

  if (opt_hmm != NULL) *opt_hmm = hmm; else p7_hmm_Destroy(hmm);
  if (opt_gm  != NULL) *opt_gm  = gm;  else p7_profile_Destroy(gm);
  *ret_om = om;
  return eslOK;

 ERROR:
  if (hmm != NULL) p7_hmm_Destroy(hmm);
  if (gm  != NULL) p7_profile_Destroy(gm);
  if (om  != NULL) (*p7_oprofile_Destroy)(om);
  *ret_om = NULL;
  return status;
}

/* Function:  p7_oprofile_Compare()
 * Synopsis:  Compare two profiles for equality (scalar fields only).
 *
 * Note: Vector field comparison requires ISA-specific code; add later.
 */
int
p7_oprofile_Compare(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg)
{
  if (om1->M    != om2->M)    ESL_FAIL(eslFAIL, errmsg, "M mismatch");
  if (om1->mode != om2->mode) ESL_FAIL(eslFAIL, errmsg, "mode mismatch");
  if (esl_FCompare(om1->scale_b, om2->scale_b, tol, 0.0) != eslOK) ESL_FAIL(eslFAIL, errmsg, "scale_b mismatch");
  if (esl_FCompare(om1->scale_w, om2->scale_w, tol, 0.0) != eslOK) ESL_FAIL(eslFAIL, errmsg, "scale_w mismatch");
  return eslOK;
}

/* p7_profile_SameAsMF(), p7_profile_SameAsVF() — defined in impl_sse/p7_oprofile.c
 * and called only from unit tests.  When building impl_avx, the _sse versions
 * of these are used via the existing SSE compilation unit (p7_oprofile_sse.c).
 * These stubs forward to the SSE implementations.
 */
int
p7_profile_SameAsMF(const P7_OPROFILE *om, P7_PROFILE *gm)
{
#ifdef eslENABLE_SSE
  /* Delegates to the SSE implementation since it accesses rbv/twv */
  extern int p7_profile_SameAsMF_sse(const P7_OPROFILE *om, P7_PROFILE *gm);
  return p7_profile_SameAsMF_sse(om, gm);
#else
  ESL_EXCEPTION(eslENORESULT, "p7_profile_SameAsMF requires SSE");
#endif
}

int
p7_profile_SameAsVF(const P7_OPROFILE *om, P7_PROFILE *gm)
{
#ifdef eslENABLE_SSE
  extern int p7_profile_SameAsVF_sse(const P7_OPROFILE *om, P7_PROFILE *gm);
  return p7_profile_SameAsVF_sse(om, gm);
#else
  ESL_EXCEPTION(eslENORESULT, "p7_profile_SameAsVF requires SSE");
#endif
}

/* p7_oprofile_Logify() — ISA-independent dispatch */
int
p7_oprofile_Logify(P7_OPROFILE *om)
{
#ifdef eslENABLE_AVX512
  if (om->allocQ4_avx512 > 0) { extern int p7_oprofile_Logify_avx512(P7_OPROFILE *om); return p7_oprofile_Logify_avx512(om); }
#endif
#ifdef eslENABLE_AVX
  if (om->allocQ4_avx    > 0) { extern int p7_oprofile_Logify_avx(P7_OPROFILE *om);    return p7_oprofile_Logify_avx(om); }
#endif
#ifdef eslENABLE_SSE
  { extern int p7_oprofile_Logify_sse(P7_OPROFILE *om); return p7_oprofile_Logify_sse(om); }
#endif
  ESL_EXCEPTION(eslENORESULT, "p7_oprofile_Logify: no SIMD implementation available");
}

/* Block management — ISA-independent */
P7_OM_BLOCK *
p7_oprofile_CreateBlock(int size)
{
  int          status;
  P7_OM_BLOCK *block = NULL;
  ESL_ALLOC(block, sizeof(P7_OM_BLOCK));
  block->count    = 0;
  block->listSize = 0;
  block->list     = NULL;
  ESL_ALLOC(block->list, sizeof(P7_OPROFILE *) * size);
  block->listSize = size;
  return block;
 ERROR:
  if (block != NULL) free(block);
  return NULL;
}

void
p7_oprofile_DestroyBlock(P7_OM_BLOCK *block)
{
  int i;
  if (block == NULL) return;
  if (block->list != NULL) {
    for (i = 0; i < block->count; i++)
      if (block->list[i] != NULL)
        (*p7_oprofile_Destroy)(block->list[i]);
    free(block->list);
  }
  free(block);
}
