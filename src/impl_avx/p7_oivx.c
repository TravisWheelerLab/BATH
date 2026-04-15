/* p7_oivx.c — runtime dispatch for P7_OIVX management.
 *
 * Contains:
 *   1. Function pointer variable definitions (one definition each; extern in impl_avx.h).
 *
 * ISA-specific implementations live in:
 *   p7_oivx_sse.c — SSE (128-bit) versions, suffix _sse
 *   p7_oivx_avx.c — AVX-256 versions, suffix _avx
 *
 * Function pointers are assigned by impl_Init() in p7_oprofile.c at startup.
 */
#include "p7_config.h"

#include "easel.h"

#include "hmmer.h"
#include "impl_avx.h"


/*****************************************************************
 * 1. Function pointer variable definitions.
 *
 * Each is initialized to NULL; impl_Init() assigns them to the
 * fastest available ISA implementation.  All are declared extern
 * in impl_avx.h.
 *****************************************************************/

P7_OIVX *(*p7_oivx_Create) (int M_hint, int C)          = NULL;
int       (*p7_oivx_GrowTo) (P7_OIVX *ov, int M, int C) = NULL;
void      (*p7_oivx_Destroy)(P7_OIVX *ov)               = NULL;
