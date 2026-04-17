/* p7_fs_oprofile.c — runtime dispatch for P7_FS_OPROFILE management.
 *
 * Contains:
 *   1. Self-patching dispatcher function pointer definitions (extern in impl_avx.h).
 *   2. ISA-independent functions (IsLocal).
 *
 * ISA-specific implementations live in:
 *   p7_fs_oprofile_sse.c — SSE (128-bit) versions, suffix _sse
 *   p7_fs_oprofile_avx.c — AVX-256 versions, suffix _avx
 *
 * Each pointer is initialized to a static dispatcher.  On the first call,
 * the dispatcher detects the ISA, patches the pointer, and calls the
 * implementation.  All subsequent calls skip dispatch entirely.
 */
#include "p7_config.h"

#include "easel.h"
#include "esl_cpu.h"

#include "hmmer.h"
#include "impl_avx.h"


/*****************************************************************
 * 1. Self-patching dispatcher function pointer definitions.
 *
 * Forward declarations, then global pointer = dispatcher,
 * then dispatcher body.  All are declared extern in impl_avx.h.
 *****************************************************************/

static P7_FS_OPROFILE *p7_fs_oprofile_Create_Dispatcher         (int M, const ESL_ALPHABET *abc, int codon_lengths);
static void            p7_fs_oprofile_Destroy_Dispatcher        (P7_FS_OPROFILE *om_fs);
static P7_FS_OPROFILE *p7_fs_oprofile_Clone_Dispatcher          (const P7_FS_OPROFILE *om_fs);
static int             p7_fs_oprofile_Convert_Dispatcher        (const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs);
static int             p7_fs_oprofile_Convert_Log_Dispatcher    (const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs);
static int             p7_fs_oprofile_SubConvert_Log_Dispatcher (const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs, int k_start, int k_end);
static int             p7_fs_oprofile_ReconfigLength_Dispatcher (P7_FS_OPROFILE *om_fs, int L);
static int             p7_fs_oprofile_ReconfigLength_Log_Dispatcher(P7_FS_OPROFILE *om_fs, int L);
static int             p7_fs_oprofile_ReconfigMultihit_Dispatcher(P7_FS_OPROFILE *om_fs, int L);
static int             p7_fs_oprofile_ReconfigUnihit_Dispatcher  (P7_FS_OPROFILE *om_fs, int L);
static int             p7_fs_oprofile_Logify_Dispatcher          (P7_FS_OPROFILE *om_fs);

P7_FS_OPROFILE *(*p7_fs_oprofile_Create)         (int M, const ESL_ALPHABET *abc, int codon_lengths)                             = p7_fs_oprofile_Create_Dispatcher;
void            (*p7_fs_oprofile_Destroy)        (P7_FS_OPROFILE *om_fs)                                                         = p7_fs_oprofile_Destroy_Dispatcher;
P7_FS_OPROFILE *(*p7_fs_oprofile_Clone)          (const P7_FS_OPROFILE *om_fs)                                                   = p7_fs_oprofile_Clone_Dispatcher;
int             (*p7_fs_oprofile_Convert)        (const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs)                             = p7_fs_oprofile_Convert_Dispatcher;
int             (*p7_fs_oprofile_Convert_Log)    (const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs)                             = p7_fs_oprofile_Convert_Log_Dispatcher;
int             (*p7_fs_oprofile_SubConvert_Log) (const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs, int k_start, int k_end)     = p7_fs_oprofile_SubConvert_Log_Dispatcher;
int             (*p7_fs_oprofile_ReconfigLength) (P7_FS_OPROFILE *om_fs, int L)                                                  = p7_fs_oprofile_ReconfigLength_Dispatcher;
int             (*p7_fs_oprofile_ReconfigLength_Log)(P7_FS_OPROFILE *om_fs, int L)                                               = p7_fs_oprofile_ReconfigLength_Log_Dispatcher;
int             (*p7_fs_oprofile_ReconfigMultihit)(P7_FS_OPROFILE *om_fs, int L)                                                 = p7_fs_oprofile_ReconfigMultihit_Dispatcher;
int             (*p7_fs_oprofile_ReconfigUnihit)  (P7_FS_OPROFILE *om_fs, int L)                                                 = p7_fs_oprofile_ReconfigUnihit_Dispatcher;
int             (*p7_fs_oprofile_Logify)          (P7_FS_OPROFILE *om_fs)                                                        = p7_fs_oprofile_Logify_Dispatcher;

static P7_FS_OPROFILE *
p7_fs_oprofile_Create_Dispatcher(int M, const ESL_ALPHABET *abc, int codon_lengths)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_fs_oprofile_Create = p7_fs_oprofile_Create_avx512; return p7_fs_oprofile_Create_avx512(M, abc, codon_lengths); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_fs_oprofile_Create = p7_fs_oprofile_Create_avx;    return p7_fs_oprofile_Create_avx(M, abc, codon_lengths); }
#endif
#ifdef eslENABLE_SSE
  p7_fs_oprofile_Create = p7_fs_oprofile_Create_sse;
  return p7_fs_oprofile_Create_sse(M, abc, codon_lengths);
#else
  p7_Die("p7_fs_oprofile_Create: no SIMD implementation available");
  return NULL;
#endif
}

static void
p7_fs_oprofile_Destroy_Dispatcher(P7_FS_OPROFILE *om_fs)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_fs_oprofile_Destroy = p7_fs_oprofile_Destroy_avx512; p7_fs_oprofile_Destroy_avx512(om_fs); return; }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_fs_oprofile_Destroy = p7_fs_oprofile_Destroy_avx;    p7_fs_oprofile_Destroy_avx(om_fs); return; }
#endif
#ifdef eslENABLE_SSE
  p7_fs_oprofile_Destroy = p7_fs_oprofile_Destroy_sse;
  p7_fs_oprofile_Destroy_sse(om_fs);
#else
  p7_Die("p7_fs_oprofile_Destroy: no SIMD implementation available");
#endif
}

static P7_FS_OPROFILE *
p7_fs_oprofile_Clone_Dispatcher(const P7_FS_OPROFILE *om_fs)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_fs_oprofile_Clone = p7_fs_oprofile_Clone_avx512; return p7_fs_oprofile_Clone_avx512(om_fs); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_fs_oprofile_Clone = p7_fs_oprofile_Clone_avx;    return p7_fs_oprofile_Clone_avx(om_fs); }
#endif
#ifdef eslENABLE_SSE
  p7_fs_oprofile_Clone = p7_fs_oprofile_Clone_sse;
  return p7_fs_oprofile_Clone_sse(om_fs);
#else
  p7_Die("p7_fs_oprofile_Clone: no SIMD implementation available");
  return NULL;
#endif
}

static int
p7_fs_oprofile_Convert_Dispatcher(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_fs_oprofile_Convert = p7_fs_oprofile_Convert_avx512; return p7_fs_oprofile_Convert_avx512(gm_fs, om_fs); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_fs_oprofile_Convert = p7_fs_oprofile_Convert_avx;    return p7_fs_oprofile_Convert_avx(gm_fs, om_fs); }
#endif
#ifdef eslENABLE_SSE
  p7_fs_oprofile_Convert = p7_fs_oprofile_Convert_sse;
  return p7_fs_oprofile_Convert_sse(gm_fs, om_fs);
#else
  p7_Die("p7_fs_oprofile_Convert: no SIMD implementation available");
  return eslENORESULT;
#endif
}

static int
p7_fs_oprofile_Convert_Log_Dispatcher(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_fs_oprofile_Convert_Log = p7_fs_oprofile_Convert_Log_avx512; return p7_fs_oprofile_Convert_Log_avx512(gm_fs, om_fs); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_fs_oprofile_Convert_Log = p7_fs_oprofile_Convert_Log_avx;    return p7_fs_oprofile_Convert_Log_avx(gm_fs, om_fs); }
#endif
#ifdef eslENABLE_SSE
  p7_fs_oprofile_Convert_Log = p7_fs_oprofile_Convert_Log_sse;
  return p7_fs_oprofile_Convert_Log_sse(gm_fs, om_fs);
#else
  p7_Die("p7_fs_oprofile_Convert_Log: no SIMD implementation available");
  return eslENORESULT;
#endif
}

static int
p7_fs_oprofile_SubConvert_Log_Dispatcher(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs, int k_start, int k_end)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_fs_oprofile_SubConvert_Log = p7_fs_oprofile_SubConvert_Log_avx512; return p7_fs_oprofile_SubConvert_Log_avx512(gm_fs, om_fs, k_start, k_end); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_fs_oprofile_SubConvert_Log = p7_fs_oprofile_SubConvert_Log_avx;    return p7_fs_oprofile_SubConvert_Log_avx(gm_fs, om_fs, k_start, k_end); }
#endif
#ifdef eslENABLE_SSE
  p7_fs_oprofile_SubConvert_Log = p7_fs_oprofile_SubConvert_Log_sse;
  return p7_fs_oprofile_SubConvert_Log_sse(gm_fs, om_fs, k_start, k_end);
#else
  p7_Die("p7_fs_oprofile_SubConvert_Log: no SIMD implementation available");
  return eslENORESULT;
#endif
}

static int
p7_fs_oprofile_ReconfigLength_Dispatcher(P7_FS_OPROFILE *om_fs, int L)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_fs_oprofile_ReconfigLength = p7_fs_oprofile_ReconfigLength_avx512; return p7_fs_oprofile_ReconfigLength_avx512(om_fs, L); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_fs_oprofile_ReconfigLength = p7_fs_oprofile_ReconfigLength_avx;    return p7_fs_oprofile_ReconfigLength_avx(om_fs, L); }
#endif
#ifdef eslENABLE_SSE
  p7_fs_oprofile_ReconfigLength = p7_fs_oprofile_ReconfigLength_sse;
  return p7_fs_oprofile_ReconfigLength_sse(om_fs, L);
#else
  p7_Die("p7_fs_oprofile_ReconfigLength: no SIMD implementation available");
  return eslENORESULT;
#endif
}

static int
p7_fs_oprofile_ReconfigLength_Log_Dispatcher(P7_FS_OPROFILE *om_fs, int L)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_fs_oprofile_ReconfigLength_Log = p7_fs_oprofile_ReconfigLength_Log_avx512; return p7_fs_oprofile_ReconfigLength_Log_avx512(om_fs, L); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_fs_oprofile_ReconfigLength_Log = p7_fs_oprofile_ReconfigLength_Log_avx;    return p7_fs_oprofile_ReconfigLength_Log_avx(om_fs, L); }
#endif
#ifdef eslENABLE_SSE
  p7_fs_oprofile_ReconfigLength_Log = p7_fs_oprofile_ReconfigLength_Log_sse;
  return p7_fs_oprofile_ReconfigLength_Log_sse(om_fs, L);
#else
  p7_Die("p7_fs_oprofile_ReconfigLength_Log: no SIMD implementation available");
  return eslENORESULT;
#endif
}

static int
p7_fs_oprofile_ReconfigMultihit_Dispatcher(P7_FS_OPROFILE *om_fs, int L)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_fs_oprofile_ReconfigMultihit = p7_fs_oprofile_ReconfigMultihit_avx512; return p7_fs_oprofile_ReconfigMultihit_avx512(om_fs, L); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_fs_oprofile_ReconfigMultihit = p7_fs_oprofile_ReconfigMultihit_avx;    return p7_fs_oprofile_ReconfigMultihit_avx(om_fs, L); }
#endif
#ifdef eslENABLE_SSE
  p7_fs_oprofile_ReconfigMultihit = p7_fs_oprofile_ReconfigMultihit_sse;
  return p7_fs_oprofile_ReconfigMultihit_sse(om_fs, L);
#else
  p7_Die("p7_fs_oprofile_ReconfigMultihit: no SIMD implementation available");
  return eslENORESULT;
#endif
}

static int
p7_fs_oprofile_ReconfigUnihit_Dispatcher(P7_FS_OPROFILE *om_fs, int L)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_fs_oprofile_ReconfigUnihit = p7_fs_oprofile_ReconfigUnihit_avx512; return p7_fs_oprofile_ReconfigUnihit_avx512(om_fs, L); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_fs_oprofile_ReconfigUnihit = p7_fs_oprofile_ReconfigUnihit_avx;    return p7_fs_oprofile_ReconfigUnihit_avx(om_fs, L); }
#endif
#ifdef eslENABLE_SSE
  p7_fs_oprofile_ReconfigUnihit = p7_fs_oprofile_ReconfigUnihit_sse;
  return p7_fs_oprofile_ReconfigUnihit_sse(om_fs, L);
#else
  p7_Die("p7_fs_oprofile_ReconfigUnihit: no SIMD implementation available");
  return eslENORESULT;
#endif
}

static int
p7_fs_oprofile_Logify_Dispatcher(P7_FS_OPROFILE *om_fs)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_fs_oprofile_Logify = p7_fs_oprofile_Logify_avx512; return p7_fs_oprofile_Logify_avx512(om_fs); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_fs_oprofile_Logify = p7_fs_oprofile_Logify_avx;    return p7_fs_oprofile_Logify_avx(om_fs); }
#endif
#ifdef eslENABLE_SSE
  p7_fs_oprofile_Logify = p7_fs_oprofile_Logify_sse;
  return p7_fs_oprofile_Logify_sse(om_fs);
#else
  p7_Die("p7_fs_oprofile_Logify: no SIMD implementation available");
  return eslENORESULT;
#endif
}


/*****************************************************************
 * 2. ISA-independent functions.
 *****************************************************************/

/* Function:  p7_fs_oprofile_IsLocal()
 * Synopsis:  Returns TRUE if FS profile is in local alignment mode.
 */
int
p7_fs_oprofile_IsLocal(const P7_FS_OPROFILE *om_fs)
{
  if (om_fs->mode == p7_LOCAL || om_fs->mode == p7_UNILOCAL) return TRUE;
  return FALSE;
}
