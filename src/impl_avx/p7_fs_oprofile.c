/* p7_fs_oprofile.c — runtime dispatch for P7_FS_OPROFILE management.
 *
 * Contains:
 *   1. Function pointer variable definitions (one definition each; extern in impl_avx.h).
 *   2. ISA-independent functions (IsLocal).
 *
 * ISA-specific implementations live in:
 *   p7_fs_oprofile_sse.c — SSE (128-bit) versions, suffix _sse
 *   p7_fs_oprofile_avx.c — AVX-256 versions, suffix _avx
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

P7_FS_OPROFILE *(*p7_fs_oprofile_Create)         (int M, const ESL_ALPHABET *abc, int codon_lengths)                             = NULL;
void            (*p7_fs_oprofile_Destroy)        (P7_FS_OPROFILE *om_fs)                                                         = NULL;
P7_FS_OPROFILE *(*p7_fs_oprofile_Clone)          (const P7_FS_OPROFILE *om_fs)                                                   = NULL;
int             (*p7_fs_oprofile_Convert)        (const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs)                             = NULL;
int             (*p7_fs_oprofile_Convert_Log)    (const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs)                             = NULL;
int             (*p7_fs_oprofile_SubConvert_Log) (const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs, int k_start, int k_end)     = NULL;
int             (*p7_fs_oprofile_ReconfigLength) (P7_FS_OPROFILE *om_fs, int L)                                                  = NULL;
int             (*p7_fs_oprofile_ReconfigLength_Log)(P7_FS_OPROFILE *om_fs, int L)                                               = NULL;
int             (*p7_fs_oprofile_ReconfigMultihit)(P7_FS_OPROFILE *om_fs, int L)                                                 = NULL;
int             (*p7_fs_oprofile_ReconfigUnihit)  (P7_FS_OPROFILE *om_fs, int L)                                                 = NULL;
int             (*p7_fs_oprofile_Logify)          (P7_FS_OPROFILE *om_fs)                                                        = NULL;


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
