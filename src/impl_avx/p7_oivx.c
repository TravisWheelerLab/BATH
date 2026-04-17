/* p7_oivx.c — runtime dispatch for P7_OIVX management.
 *
 * Contains:
 *   1. Self-patching dispatcher function pointer definitions (extern in impl_avx.h).
 *
 * ISA-specific implementations live in:
 *   p7_oivx_sse.c — SSE (128-bit) versions, suffix _sse
 *   p7_oivx_avx.c — AVX-256 versions, suffix _avx
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
 *****************************************************************/

static P7_OIVX *p7_oivx_Create_Dispatcher (int M_hint, int C);
static int      p7_oivx_GrowTo_Dispatcher (P7_OIVX *ov, int M, int C);
static void     p7_oivx_Destroy_Dispatcher(P7_OIVX *ov);

P7_OIVX *(*p7_oivx_Create) (int M_hint, int C)          = p7_oivx_Create_Dispatcher;
int       (*p7_oivx_GrowTo) (P7_OIVX *ov, int M, int C) = p7_oivx_GrowTo_Dispatcher;
void      (*p7_oivx_Destroy)(P7_OIVX *ov)               = p7_oivx_Destroy_Dispatcher;

static P7_OIVX *
p7_oivx_Create_Dispatcher(int M_hint, int C)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_oivx_Create = p7_oivx_Create_avx512; return p7_oivx_Create_avx512(M_hint, C); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_oivx_Create = p7_oivx_Create_avx;    return p7_oivx_Create_avx(M_hint, C); }
#endif
#ifdef eslENABLE_SSE
  p7_oivx_Create = p7_oivx_Create_sse;
  return p7_oivx_Create_sse(M_hint, C);
#else
  p7_Die("p7_oivx_Create: no SIMD implementation available");
  return NULL;
#endif
}

static int
p7_oivx_GrowTo_Dispatcher(P7_OIVX *ov, int M, int C)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_oivx_GrowTo = p7_oivx_GrowTo_avx512; return p7_oivx_GrowTo_avx512(ov, M, C); }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_oivx_GrowTo = p7_oivx_GrowTo_avx;    return p7_oivx_GrowTo_avx(ov, M, C); }
#endif
#ifdef eslENABLE_SSE
  p7_oivx_GrowTo = p7_oivx_GrowTo_sse;
  return p7_oivx_GrowTo_sse(ov, M, C);
#else
  p7_Die("p7_oivx_GrowTo: no SIMD implementation available");
  return eslENORESULT;
#endif
}

static void
p7_oivx_Destroy_Dispatcher(P7_OIVX *ov)
{
#ifdef eslENABLE_AVX512
  if (esl_cpu_has_avx512()) { p7_oivx_Destroy = p7_oivx_Destroy_avx512; p7_oivx_Destroy_avx512(ov); return; }
#endif
#ifdef eslENABLE_AVX
  if (esl_cpu_has_avx())    { p7_oivx_Destroy = p7_oivx_Destroy_avx;    p7_oivx_Destroy_avx(ov); return; }
#endif
#ifdef eslENABLE_SSE
  p7_oivx_Destroy = p7_oivx_Destroy_sse;
  p7_oivx_Destroy_sse(ov);
#else
  p7_Die("p7_oivx_Destroy: no SIMD implementation available");
#endif
}
