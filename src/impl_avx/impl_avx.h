/* x86 SIMD optimized implementation: structures, declarations, and macros.
 * Supports SSE, AVX-256, and AVX-512 selected at runtime via CPU detection.
 *
 * Based on impl_sse/impl_sse.h; extended for runtime dispatch across
 * SSE (128-bit), AVX-256 (256-bit), and AVX-512 (512-bit) on x86.
 */
#ifndef P7_IMPL_AVX_INCLUDED
#define P7_IMPL_AVX_INCLUDED

#include "p7_config.h"

#include <stdint.h>      /* uint8_t, int16_t, etc.  */
#include <sys/types.h>   /* off_t                   */

#include "esl_alphabet.h"
#include "esl_random.h"
#include "esl_cpu.h"

#ifdef eslENABLE_SSE
#include <xmmintrin.h>   /* SSE  */
#include <emmintrin.h>   /* SSE2 */
#ifdef __SSE3__
#include <pmmintrin.h>   /* DENORMAL_MODE */
#endif
#include "esl_sse.h"
#endif

#ifdef eslENABLE_AVX
#include <immintrin.h>   /* AVX, AVX2 */
#include "esl_avx.h"
#endif

#ifdef eslENABLE_AVX512
#include <immintrin.h>   /* AVX-512 */
#include "esl_avx512.h"
#endif

#include "hmmer.h"

/*****************************************************************
 * Q macros: number of SIMD vectors needed per DP row
 *
 * SSE  (128-bit): 16 uint8, 8 int16, 4 float per vector
 * AVX  (256-bit): 32 uint8, 16 int16, 8 float per vector
 * AVX512 (512-bit): 64 uint8, 32 int16, 16 float per vector
 *
 * Always at least 2 vectors; striped implementation requires it.
 *****************************************************************/

/* SSE */
#define p7O_NQB(M)        ( ESL_MAX(2, ((((M)-1) / 16) + 1)))  /* 16 uchars/vector  */
#define p7O_NQW(M)        ( ESL_MAX(2, ((((M)-1) / 8)  + 1)))  /*  8 int16s/vector  */
#define p7O_NQF(M)        ( ESL_MAX(2, ((((M)-1) / 4)  + 1)))  /*  4 floats/vector  */

/* AVX-256 */
#define p7O_NQB_AVX(M)    ( ESL_MAX(2, ((((M)-1) / 32) + 1)))  /* 32 uchars/vector  */
#define p7O_NQW_AVX(M)    ( ESL_MAX(2, ((((M)-1) / 16) + 1)))  /* 16 int16s/vector  */
#define p7O_NQF_AVX(M)    ( ESL_MAX(2, ((((M)-1) / 8)  + 1)))  /*  8 floats/vector  */

/* AVX-512 */
#define p7O_NQB_AVX512(M) ( ESL_MAX(2, ((((M)-1) / 64) + 1)))  /* 64 uchars/vector  */
#define p7O_NQW_AVX512(M) ( ESL_MAX(2, ((((M)-1) / 32) + 1)))  /* 32 int16s/vector  */
#define p7O_NQF_AVX512(M) ( ESL_MAX(2, ((((M)-1) / 16) + 1)))  /* 16 floats/vector  */

#define p7O_EXTRA_SB 17    /* see ssvfilter.c for explanation */


/*****************************************************************
 * 1. P7_OPROFILE: an optimized score profile
 *
 * Each ISA section is compiled in if the corresponding eslENABLE_*
 * flag is set.  All three can be active simultaneously in a single
 * binary (runtime dispatch model).
 *****************************************************************/

#define p7O_NXSTATES  4    /* special states stored: ENJC                       */
#define p7O_NXTRANS   2    /* special states all have 2 transitions: move, loop */
#define p7O_NTRANS    8    /* 7 core transitions + BMk entry                    */
enum p7o_xstates_e      { p7O_E    = 0, p7O_N    = 1,  p7O_J  = 2,  p7O_C  = 3 };
enum p7o_xtransitions_e { p7O_MOVE = 0, p7O_LOOP = 1 };
enum p7o_tsc_e          { p7O_BM   = 0, p7O_MM   = 1,  p7O_IM = 2,  p7O_DM = 3,
                          p7O_MD   = 4, p7O_MI   = 5,  p7O_II = 6,  p7O_DD = 7 };

typedef struct p7_oprofile_s {

  /* MSVFilter: scaled, biased uchars */
#ifdef eslENABLE_SSE
  __m128i **rbv;           /* match scores [x][q]   [Kp][Q16]          */
  __m128i **sbv;           /* match scores for ssvfilter                */
#endif
#ifdef eslENABLE_AVX
  __m256i **rbv_avx;
  __m256i **sbv_avx;
#endif
#ifdef eslENABLE_AVX512
  __m512i **rbv_avx512;
  __m512i **sbv_avx512;
#endif
  uint8_t   tbm_b;         /* constant B->Mk cost:    scaled log 2/M(M+1)  */
  uint8_t   tec_b;         /* constant E->C  cost:    scaled log 0.5       */
  uint8_t   tjb_b;         /* constant NCJ move cost: scaled log 3/(L+3)   */
  float     scale_b;       /* typically 3 / log2: scores scale to 1/3 bits */
  uint8_t   base_b;        /* typically +190: offset of uchar scores       */
  uint8_t   bias_b;        /* positive bias to emission scores, >=0        */

  /* ViterbiFilter: scaled signed 16-bit integer words */
#ifdef eslENABLE_SSE
  __m128i **rwv;           /* [x][q]  [Kp][Q8]                            */
  __m128i  *twv;           /* transition score blocks [8*Q8]              */
#endif
#ifdef eslENABLE_AVX
  __m256i **rwv_avx;
  __m256i  *twv_avx;
#endif
#ifdef eslENABLE_AVX512
  __m512i **rwv_avx512;
  __m512i  *twv_avx512;
#endif
  int16_t   xw[p7O_NXSTATES][p7O_NXTRANS]; /* NECJ state transition costs */
  float     scale_w;       /* score units: typically 500/log(2), 1/500 bits */
  int16_t   base_w;        /* offset of sword scores: typically +12000      */
  int16_t   ddbound_w;     /* threshold for lazy DD evaluation              */
  float     ncj_roundoff;  /* missing precision on NN,CC,JJ after rounding  */

  /* Forward/Backward: IEEE754 single-precision floats */
#ifdef eslENABLE_SSE
  __m128 **rfv;            /* [x][q]  [Kp][Q4]                            */
  __m128  *tfv;            /* transition probability blocks [8*Q4]        */
#endif
#ifdef eslENABLE_AVX
  __m256 **rfv_avx;
  __m256  *tfv_avx;
#endif
#ifdef eslENABLE_AVX512
  __m512 **rfv_avx512;
  __m512  *tfv_avx512;
#endif
  float    xf[p7O_NXSTATES][p7O_NXTRANS]; /* NECJ transition costs        */

  /* Raw malloc pointers (before alignment) */
#ifdef eslENABLE_SSE
  __m128i  *rbv_mem;
  __m128i  *sbv_mem;
  __m128i  *rwv_mem;
  __m128i  *twv_mem;
  __m128   *tfv_mem;
  __m128   *rfv_mem;
#endif
#ifdef eslENABLE_AVX
  __m256i  *rbv_mem_avx;
  __m256i  *sbv_mem_avx;
  __m256i  *rwv_mem_avx;
  __m256i  *twv_mem_avx;
  __m256   *tfv_mem_avx;
  __m256   *rfv_mem_avx;
#endif
#ifdef eslENABLE_AVX512
  __m512i  *rbv_mem_avx512;
  __m512i  *sbv_mem_avx512;
  __m512i  *rwv_mem_avx512;
  __m512i  *twv_mem_avx512;
  __m512   *tfv_mem_avx512;
  __m512   *rfv_mem_avx512;
#endif

  /* Disk offset information for hmmpfam's fast model retrieval */
  off_t  offs[p7_NOFFSETS];   /* p7_{MFP}OFFSET, or -1                    */
  off_t  roff;                /* record offset (start of record); -1 if none */
  off_t  eoff;                /* offset to last byte of record; -1 if unknown */

  /* Annotation copied from parent profile */
  char  *name;
  char  *acc;
  char  *desc;
  char  *rf;
  char  *mm;
  char  *cs;
  char  *consensus;
  float  evparam[p7_NEVPARAM];
  float  cutoff[p7_NCUTOFFS];
  float  compo[p7_MAXABET];
  const ESL_ALPHABET *abc;

  /* Configuration, size, allocation */
  int    L;
  int    M;
  int    max_length;
  int    allocM;
  int    allocQ4;            /* p7O_NQF(allocM):    SSE float stripe count  */
  int    allocQ8;            /* p7O_NQW(allocM):    SSE int16 stripe count  */
  int    allocQ16;           /* p7O_NQB(allocM):    SSE uint8 stripe count  */
  int    allocQ4_avx;        /* p7O_NQF_AVX(allocM)                        */
  int    allocQ8_avx;        /* p7O_NQW_AVX(allocM)                        */
  int    allocQ16_avx;       /* p7O_NQB_AVX(allocM)                        */
  int    allocQ4_avx512;     /* p7O_NQF_AVX512(allocM)                     */
  int    allocQ8_avx512;     /* p7O_NQW_AVX512(allocM)                     */
  int    allocQ16_avx512;    /* p7O_NQB_AVX512(allocM)                     */
  int    mode;               /* currently must be p7_LOCAL                 */
  float  nj;                 /* expected # of J's: 0 or 1                  */
  int    clone;              /* if set, pointers borrowed; do not free      */

} P7_OPROFILE;

typedef struct {
  int            count;
  int            listSize;
  P7_OPROFILE  **list;
} P7_OM_BLOCK;

/* retrieve match odds ratio [k][x] — SSE version */
#ifdef eslENABLE_SSE
static inline float
p7_oprofile_FGetEmission(const P7_OPROFILE *om, int k, int x)
{
  union { __m128 v; float p[4]; } u;
  int Q = p7O_NQF(om->M);
  u.v = om->rfv[x][(k-1) % Q];
  return u.p[(k-1)/Q];
}
#endif


/*****************************************************************
 * 2. P7_FS_OPROFILE: optimized frameshift score profile
 *
 * Stores only IEEE754 single-precision float vectors (no byte or
 * word variants).  The frameshift algorithm has no MSV/Viterbi
 * filter stage, so only float vectors are needed.
 *
 * Emission layout: rfv[c][q] for codon index c and stripe q.
 * Transition layout: tfv[p7O_NTRANS * allocQ4] in stripe-major order.
 * Special state (ENJC) costs: scalar xf[][].
 *****************************************************************/
typedef struct p7_fs_oprofile_s {

  /* IEEE754 single-precision float vectors */
#ifdef eslENABLE_SSE
  __m128  **rfv;           /* match emission scores [c][q]                */
  __m128   *tfv;           /* transition score blocks [p7O_NTRANS*allocQ4] */
  __m128   *rfv_mem;       /* raw malloc before 16-byte alignment          */
  __m128   *tfv_mem;
  int       allocQ4;       /* p7O_NQF(allocM): SSE float stripe count      */
#endif
#ifdef eslENABLE_AVX
  __m256  **rfv_avx;
  __m256   *tfv_avx;
  __m256   *rfv_mem_avx;
  __m256   *tfv_mem_avx;
  int       allocQ4_avx;   /* p7O_NQF_AVX(allocM)                         */
#endif
#ifdef eslENABLE_AVX512
  __m512  **rfv_avx512;
  __m512   *tfv_avx512;
  __m512   *rfv_mem_avx512;
  __m512   *tfv_mem_avx512;
  int       allocQ4_avx512; /* p7O_NQF_AVX512(allocM)                     */
#endif

  float     xf[p7O_NXSTATES][p7O_NXTRANS]; /* ENJC special state transition costs */

  /* Frameshift-specific parameters */
  int       codon_lengths; /* number of codon lengths                      */
  float     fsprob;        /* frameshift penalty (log-odds score)          */

  /* Disk offset information */
  off_t  offs[p7_NOFFSETS];
  off_t  roff;
  off_t  eoff;

  /* Annotation copied from parent profile */
  char  *name;
  char  *acc;
  char  *desc;
  char  *rf;
  char  *mm;
  char  *cs;
  char  *consensus;
  float  evparam[p7_NEVPARAM];
  float  cutoff[p7_NCUTOFFS];
  float  compo[p7_MAXABET];
  const ESL_ALPHABET *abc;

  /* Configuration, size, allocation */
  int    L;
  int    M;
  int    max_length;
  int    allocM;
  int    mode;
  float  nj;
  int    clone;

} P7_FS_OPROFILE;

/* Retrieve float match emission score for model position k, codon index c.
 * Used for display/debugging; inner loops use the striped vectors directly.
 * ISA-specific versions select the right vector width.
 */
#ifdef eslENABLE_SSE
static inline float
p7_fs_oprofile_FGetEmission(const P7_FS_OPROFILE *om_fs, int k, int c)
{
  union { __m128 v; float p[4]; } u;
  int Q = p7O_NQF(om_fs->M);
  u.v = om_fs->rfv[c][(k-1) % Q];
  return u.p[(k-1)/Q];
}
#endif

/* OSSX macros: index donor IVX arrays by stripe q.
 * Slot offsets match SPLICE_OFFSET_1=3, SPLICE_OFFSET_2=18,
 * p7S_SPLICE_SIGNALS=3 from p7_splice.h.
 */
#define OSSX0(q, signal)          (don_ovx[(signal)][q])
#define OSSX1(q, signal, nuc1)    (don_ovx[3  + (nuc1) * 3 + (signal)][q])
#define OSSX2(q, signal, nuc3)    (don_ovx[18 + (nuc3) * 3 + (signal)][q])


/*****************************************************************
 * 3. P7_OIVX: vectorized intermediate-value matrix
 *
 * Stores vectors in [channel][stripe] layout for the splice
 * algorithm.  Each ISA has its own pointer and stripe count.
 *****************************************************************/
typedef struct p7_oivx_s {
  int      allocM;         /* max model length currently allocated for    */
  int      allocC;         /* number of channels                          */
#ifdef eslENABLE_SSE
  __m128 **ivx;            /* [allocC][allocQ4] striped float vectors     */
  __m128  *ivx_mem;        /* flat backing memory (+15 bytes for alignment)*/
  int      allocQ4;        /* p7O_NQF(allocM)                             */
#endif
#ifdef eslENABLE_AVX
  __m256 **ivx_avx;
  __m256  *ivx_mem_avx;
  int      allocQ4_avx;    /* p7O_NQF_AVX(allocM)                        */
#endif
#ifdef eslENABLE_AVX512
  __m512 **ivx_avx512;
  __m512  *ivx_mem_avx512;
  int      allocQ4_avx512; /* p7O_NQF_AVX512(allocM)                     */
#endif
} P7_OIVX;

/* OIVXo(c,q): access channel c, stripe q.
 * Caller must declare a local pointer to the appropriate ivx array:
 *   __m128 **ivx = ov->ivx;          (SSE)
 *   __m256 **ivx = ov->ivx_avx;      (AVX)
 *   __m512 **ivx = ov->ivx_avx512;   (AVX-512)
 */
#define OIVXo(c, q)  (ivx[(c)][(q)])


/*****************************************************************
 * 4. P7_OMX: one-row dynamic programming matrix
 *
 * dp_mem is a single allocation shared by dpf/dpw/dpb.  Each ISA
 * has typed pointer arrays pointing into that memory; only the
 * pointers for the active ISA are populated at runtime.
 *****************************************************************/

enum p7x_scells_e { p7X_M = 0, p7X_D = 1, p7X_I = 2 };
#define p7X_NSCELLS 3

/* Frameshift full-matrix cell layout */
enum p7x_fscells_e  { p7X_FS_D = 0, p7X_FS_I = 1, p7X_FS_M = 2 };
enum p7x_fscodons_e {
  p7X_FS_C0 = 0,
  p7X_FS_C1 = 1,
  p7X_FS_C2 = 2,
  p7X_FS_C3 = 3,
  p7X_FS_C4 = 4,
  p7X_FS_C5 = 5,
};
#define p7X_NSCELLS_FS 8   /* D + I + M_C0..M_C5 */

enum p7x_xcells_e { p7X_E = 0, p7X_N = 1, p7X_J = 2, p7X_B = 3, p7X_C = 4, p7X_SCALE = 5 };
#define p7X_NXCELLS 6

typedef struct p7_omx_s {
  int       M;
  int       L;
  int       nscells;   /* p7X_NSCELLS (3) or p7X_NSCELLS_FS (8) */

  /* Main DP matrix pointers — only the active ISA's pointers are set */
#ifdef eslENABLE_SSE
  __m128  **dpf;       /* float vectors  [0..L][0..Q4-1][MDI]   */
  __m128i **dpw;       /* int16 vectors  [0..L][0..Q8-1][MDI]   */
  __m128i **dpb;       /* uint8 vectors  [0..L][0..Q16-1]       */
#endif
#ifdef eslENABLE_AVX
  __m256  **dpf_avx;
  __m256i **dpw_avx;
  __m256i **dpb_avx;
#endif
#ifdef eslENABLE_AVX512
  __m512  **dpf_avx512;
  __m512i **dpw_avx512;
  __m512i **dpb_avx512;
#endif

  void     *dp_mem;    /* single allocation shared by all dp arrays      */
  int       allocR;    /* allocated row count; allocR >= L+1             */
  int       validR;    /* rows actually pointing at dp_mem               */
  int       allocQ4;   /* SSE float stripe count   (allocQ4*4  >= M)     */
  int       allocQ8;   /* SSE int16 stripe count   (allocQ8*8  >= M)     */
  int       allocQ16;  /* SSE uint8 stripe count   (allocQ16*16 >= M)    */
  int       allocQ4_avx;
  int       allocQ8_avx;
  int       allocQ16_avx;
  int       allocQ4_avx512;
  int       allocQ8_avx512;
  int       allocQ16_avx512;
  size_t    ncells;    /* current allocation size of dp_mem, in cells    */

  /* X states (full/parser) or NULL (scorer) */
  float    *xmx;       /* [0..L][ENJBCS]; indexed [i*p7X_NXCELLS+s]    */
  void     *x_mem;
  int       allocXR;
  float     totscale;
  int       has_own_scales;

  int       debugging;
  FILE     *dfp;
} P7_OMX;

/* DP cell access macros.
 * Each algorithm file declares a local typed dp pointer, e.g.:
 *   __m128  **dp = ox->dpf;       (SSE float)
 *   __m256  **dp = ox->dpf_avx;   (AVX float)
 * Then uses these ISA-independent macros:
 */
#define MMXo(q)        (dp[(q) * p7X_NSCELLS    + p7X_M])
#define DMXo(q)        (dp[(q) * p7X_NSCELLS    + p7X_D])
#define IMXo(q)        (dp[(q) * p7X_NSCELLS    + p7X_I])
#define XMXo(i,s)      (xmx[(i) * p7X_NXCELLS   + s])

#define MMO(dp,q)      ((dp)[(q) * p7X_NSCELLS    + p7X_M])
#define DMO(dp,q)      ((dp)[(q) * p7X_NSCELLS    + p7X_D])
#define IMO(dp,q)      ((dp)[(q) * p7X_NSCELLS    + p7X_I])

#define MMO_FS(dp,q,c) ((dp)[(q) * p7X_NSCELLS_FS + p7X_FS_M + (c)])
#define DMO_FS(dp,q)   ((dp)[(q) * p7X_NSCELLS_FS + p7X_FS_D])
#define IMO_FS(dp,q)   ((dp)[(q) * p7X_NSCELLS_FS + p7X_FS_I])

/* Debug accessors — SSE version; AVX algorithm files use local dp pointer instead */
#ifdef eslENABLE_SSE
static inline float
p7_omx_FGetMDI(const P7_OMX *ox, int s, int i, int k)
{
  union { __m128 v; float p[4]; } u;
  int Q = p7O_NQF(ox->M);
  u.v = ox->dpf[i][p7X_NSCELLS * ((k-1) % Q) + s];
  return u.p[(k-1)/Q];
}

static inline void
p7_omx_FSetMDI(const P7_OMX *ox, int s, int i, int k, float val)
{
  union { __m128 v; float p[4]; } u;
  int Q = p7O_NQF(ox->M);
  int q = p7X_NSCELLS * ((k-1) % Q) + s;
  u.v        = ox->dpf[i][q];
  u.p[(k-1)/Q] = val;
  ox->dpf[i][q] = u.v;
}
#endif


/*****************************************************************
 * 5. External API declarations
 *****************************************************************/

/* p7_omx.c — dispatch wrappers and ISA-independent functions */
extern P7_OMX      *p7_omx_Create    (int allocM, int allocL, int allocXL);
extern int          p7_omx_GrowTo    (P7_OMX *ox, int allocM, int allocL, int allocXL);
extern P7_OMX      *p7_omx_Create_dpf(int allocM, int allocL, int allocXL, int nscells);
extern int          p7_omx_GrowTo_dpf(P7_OMX *ox, int allocM, int allocL, int allocXL);
extern int          p7_omx_FDeconvert(P7_OMX *ox, P7_GMX *gx);
extern int          p7_omx_Reuse     (P7_OMX *ox);
extern void         p7_omx_Destroy   (P7_OMX *ox);
extern int          p7_omx_SetDumpMode(FILE *fp, P7_OMX *ox, int truefalse);
extern int          p7_omx_Dump      (FILE *fp, P7_OMX *ox);
extern int          p7_omx_DumpMFRow (P7_OMX *ox, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);
extern int          p7_omx_DumpVFRow (P7_OMX *ox, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC);
extern int          p7_omx_DumpFBRow (P7_OMX *ox, int logify, int rowi, int width, int precision, float xE, float xN, float xJ, float xB, float xC);
extern int          p7_omx_DumpFBRow_FS(P7_OMX *ox, int logify, int i, int rowi, int width, int precision, float xE, float xN, float xJ, float xB, float xC);

/* p7_omx_sse.c — SSE implementations */
extern P7_OMX      *p7_omx_Create_sse    (int allocM, int allocL, int allocXL);
extern int          p7_omx_GrowTo_sse    (P7_OMX *ox, int allocM, int allocL, int allocXL);
extern P7_OMX      *p7_omx_Create_dpf_sse(int allocM, int allocL, int allocXL, int nscells);
extern int          p7_omx_GrowTo_dpf_sse(P7_OMX *ox, int allocM, int allocL, int allocXL);
extern int          p7_omx_FDeconvert_sse(P7_OMX *ox, P7_GMX *gx);
extern void         p7_omx_Destroy_sse   (P7_OMX *ox);

/* p7_oprofile.c — function pointers (dispatched at runtime) */
extern P7_OPROFILE *(*p7_oprofile_Create)(int M, const ESL_ALPHABET *abc);
extern int          p7_oprofile_IsLocal(const P7_OPROFILE *om);
extern void       (*p7_oprofile_Destroy)(P7_OPROFILE *om);
extern P7_OPROFILE *(*p7_oprofile_Clone)(const P7_OPROFILE *om);
extern int        (*p7_oprofile_Convert)(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int        (*p7_oprofile_Convert_Log)(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int        (*p7_oprofile_ReconfigLength)(P7_OPROFILE *om, int L);
extern int        (*p7_oprofile_ReconfigLength_Log)(P7_OPROFILE *om, int L);
extern int        (*p7_oprofile_ReconfigMSVLength)(P7_OPROFILE *om, int L);
extern int        (*p7_oprofile_ReconfigRestLength)(P7_OPROFILE *om, int L);
extern int        (*p7_oprofile_ReconfigMultihit)(P7_OPROFILE *om, int L);
extern int        (*p7_oprofile_ReconfigMultihit_Log)(P7_OPROFILE *om, int L);
extern int        (*p7_oprofile_ReconfigUnihit)(P7_OPROFILE *om, int L);
extern int        (*p7_oprofile_ReconfigUnihit_Log)(P7_OPROFILE *om, int L);
extern int        (*p7_oprofile_UpdateFwdEmissionScores)(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);
extern int        (*p7_oprofile_UpdateVitEmissionScores)(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);
extern int        (*p7_oprofile_UpdateMSVEmissionScores)(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);
extern int          p7_oprofile_Dump(FILE *fp, const P7_OPROFILE *om);
extern int          p7_oprofile_Sample(ESL_RANDOMNESS *r, const ESL_ALPHABET *abc, const P7_BG *bg, int M, int L,
                                       P7_HMM **opt_hmm, P7_PROFILE **opt_gm, P7_OPROFILE **ret_om);
extern int          p7_oprofile_Compare(const P7_OPROFILE *om1, const P7_OPROFILE *om2, float tol, char *errmsg);
extern int          p7_profile_SameAsMF(const P7_OPROFILE *om, P7_PROFILE *gm);
extern int          p7_profile_SameAsVF(const P7_OPROFILE *om, P7_PROFILE *gm);
extern int          p7_oprofile_GetFwdTransitionArray(const P7_OPROFILE *om, int type, float *arr);
extern int          p7_oprofile_GetSSVEmissionScoreArray(const P7_OPROFILE *om, uint8_t *arr);
extern int          p7_oprofile_GetFwdEmissionScoreArray(const P7_OPROFILE *om, float *arr);
extern int          p7_oprofile_GetFwdEmissionArray(const P7_OPROFILE *om, P7_BG *bg, float *arr);
extern int          p7_oprofile_Logify(P7_OPROFILE *om);

/* ISA-specific p7_oprofile implementations (called via function pointers above) */
#ifdef eslENABLE_SSE
extern P7_OPROFILE *p7_oprofile_Create_sse(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy_sse(P7_OPROFILE *om);
extern P7_OPROFILE *p7_oprofile_Clone_sse(const P7_OPROFILE *om);
extern int          p7_oprofile_Convert_sse(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_Convert_Log_sse(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_ReconfigLength_sse(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigLength_Log_sse(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigMSVLength_sse(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigRestLength_sse(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigMultihit_sse(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigMultihit_Log_sse(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigUnihit_sse(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigUnihit_Log_sse(P7_OPROFILE *om, int L);
extern int          p7_oprofile_UpdateFwdEmissionScores_sse(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);
extern int          p7_oprofile_UpdateVitEmissionScores_sse(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);
extern int          p7_oprofile_UpdateMSVEmissionScores_sse(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);
#endif
#ifdef eslENABLE_AVX
extern P7_OPROFILE *p7_oprofile_Create_avx(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy_avx(P7_OPROFILE *om);
extern P7_OPROFILE *p7_oprofile_Clone_avx(const P7_OPROFILE *om);
extern int          p7_oprofile_Convert_avx(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_Convert_Log_avx(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_ReconfigLength_avx(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigLength_Log_avx(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigMSVLength_avx(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigRestLength_avx(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigMultihit_avx(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigMultihit_Log_avx(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigUnihit_avx(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigUnihit_Log_avx(P7_OPROFILE *om, int L);
extern int          p7_oprofile_UpdateFwdEmissionScores_avx(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);
extern int          p7_oprofile_UpdateVitEmissionScores_avx(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);
extern int          p7_oprofile_UpdateMSVEmissionScores_avx(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);
#endif
#ifdef eslENABLE_AVX512
extern P7_OPROFILE *p7_oprofile_Create_avx512(int M, const ESL_ALPHABET *abc);
extern void         p7_oprofile_Destroy_avx512(P7_OPROFILE *om);
extern P7_OPROFILE *p7_oprofile_Clone_avx512(const P7_OPROFILE *om);
extern int          p7_oprofile_Convert_avx512(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_Convert_Log_avx512(const P7_PROFILE *gm, P7_OPROFILE *om);
extern int          p7_oprofile_ReconfigLength_avx512(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigLength_Log_avx512(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigMSVLength_avx512(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigRestLength_avx512(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigMultihit_avx512(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigMultihit_Log_avx512(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigUnihit_avx512(P7_OPROFILE *om, int L);
extern int          p7_oprofile_ReconfigUnihit_Log_avx512(P7_OPROFILE *om, int L);
extern int          p7_oprofile_UpdateFwdEmissionScores_avx512(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);
extern int          p7_oprofile_UpdateVitEmissionScores_avx512(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);
extern int          p7_oprofile_UpdateMSVEmissionScores_avx512(P7_OPROFILE *om, P7_BG *bg, float *fwd_emissions, float *sc_arr);
#endif

/* p7_fs_oprofile.c — function pointers (dispatched at runtime) */
extern P7_FS_OPROFILE *(*p7_fs_oprofile_Create)(int M, const ESL_ALPHABET *abc, int codon_lengths);
extern int              p7_fs_oprofile_IsLocal(const P7_FS_OPROFILE *om_fs);
extern void           (*p7_fs_oprofile_Destroy)(P7_FS_OPROFILE *om_fs);
extern P7_FS_OPROFILE *(*p7_fs_oprofile_Clone)(const P7_FS_OPROFILE *om_fs);
extern int            (*p7_fs_oprofile_Convert)(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs);
extern int            (*p7_fs_oprofile_Convert_Log)(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs);
extern int            (*p7_fs_oprofile_SubConvert_Log)(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs, int k_start, int k_end);
extern int            (*p7_fs_oprofile_ReconfigLength)(P7_FS_OPROFILE *om_fs, int L);
extern int            (*p7_fs_oprofile_ReconfigLength_Log)(P7_FS_OPROFILE *om_fs, int L);
extern int            (*p7_fs_oprofile_ReconfigMultihit)(P7_FS_OPROFILE *om_fs, int L);
extern int            (*p7_fs_oprofile_ReconfigUnihit)(P7_FS_OPROFILE *om_fs, int L);
extern int            (*p7_fs_oprofile_Logify)(P7_FS_OPROFILE *om_fs);

/* p7_oivx.c — function pointers (dispatched at runtime) */
extern P7_OIVX *(*p7_oivx_Create)(int M_hint, int C);
extern int      (*p7_oivx_GrowTo)(P7_OIVX *ov, int M, int C);
extern void     (*p7_oivx_Destroy)(P7_OIVX *ov);

/* decoding.c */
extern int p7_Decoding      (const P7_OPROFILE *om, const P7_OMX *oxf,       P7_OMX *oxb, P7_OMX *pp);
extern int p7_DomainDecoding(const P7_OPROFILE *om, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef);

/* decoding_fs.c */
extern int p7_Decoding_Frameshift   (const P7_FS_OPROFILE *om_fs, P7_OMX *fwd, const P7_OMX *bck);
extern int p7_DomainDecoding_Frameshift(const P7_FS_OPROFILE *om_fs, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef);

/* decoding_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Decoding_sse               (const P7_OPROFILE *om,    const P7_OMX *oxf,       P7_OMX *oxb,       P7_OMX *pp);
extern int p7_DomainDecoding_sse         (const P7_OPROFILE *om,    const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef);
extern int p7_Decoding_Frameshift_sse    (const P7_FS_OPROFILE *om_fs,                   P7_OMX *fwd, const P7_OMX *bck);
extern int p7_DomainDecoding_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef);
#endif

/* fwdback.c — dispatch wrappers */
extern int p7_Forward       (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
extern int p7_ForwardParser (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
extern int p7_Backward      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);
extern int p7_BackwardParser(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);
/* fwdback_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Forward_sse       (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
extern int p7_ForwardParser_sse (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
extern int p7_Backward_sse      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);
extern int p7_BackwardParser_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);
#endif

/* fwdback_fs.c */
extern int p7_ForwardParser_Frameshift_3Codons (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int p7_BackwardParser_Frameshift_3Codons(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);
extern int p7_ForwardParser_Frameshift_5Codons (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int p7_BackwardParser_Frameshift_5Codons(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);
extern int p7_Forward_Frameshift               (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int p7_Backward_Frameshift              (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);

/* fwdback_fs_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_ForwardParser_Frameshift_3Codons_sse (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int p7_BackwardParser_Frameshift_3Codons_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);
extern int p7_ForwardParser_Frameshift_5Codons_sse (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int p7_BackwardParser_Frameshift_5Codons_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);
extern int p7_Forward_Frameshift_sse               (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int p7_Backward_Frameshift_sse              (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);
#endif

/* p7_oprofile.c — block management (ISA-independent) */
extern P7_OM_BLOCK *p7_oprofile_CreateBlock(int size);
extern void         p7_oprofile_DestroyBlock(P7_OM_BLOCK *block);

/* msvfilter.c — dispatch wrappers */
extern int p7_MSVFilter      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_SSVFilter_BATH (const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *msvdata, P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist);
/* msvfilter_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_MSVFilter_sse      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_SSVFilter_BATH_sse (const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *ssvdata, P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist);
#endif

/* null2.c */
extern int p7_Null2_ByExpectation(const P7_OPROFILE *om, const P7_OMX *pp, float *null2);
extern int p7_Null2_ByTrace      (const P7_OPROFILE *om, const P7_TRACE *tr, int zstart, int zend, P7_OMX *wrk, float *null2);

/* null2_fs.c */
extern int p7_Null2_fs_ByExpectation(const P7_FS_OPROFILE *om_fs, P7_OMX *pp, float *null2);

/* null2_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Null2_ByExpectation_sse    (const P7_OPROFILE    *om,    const P7_OMX *pp, float *null2);
extern int p7_Null2_ByTrace_sse          (const P7_OPROFILE    *om,    const P7_TRACE *tr, int zstart, int zend, P7_OMX *wrk, float *null2);
extern int p7_Null2_fs_ByExpectation_sse (const P7_FS_OPROFILE *om_fs, P7_OMX *pp, float *null2);
#endif

/* optacc.c */
extern int p7_OptimalAccuracy(const P7_OPROFILE *om, const P7_OMX *pp,       P7_OMX *ox, float *ret_e);
extern int p7_OATrace        (const P7_OPROFILE *om, const P7_OMX *pp, const P7_OMX *ox, P7_TRACE *tr);

/* optacc_fs.c */
extern int p7_OptimalAccuracy_Frameshift(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp, P7_OMX *ox, float *ret_e);
extern int p7_OATrace_Frameshift        (const P7_FS_OPROFILE *om_fs, const P7_OMX *pp, const P7_OMX *ox, P7_TRACE *tr);

/* optacc_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_OptimalAccuracy_sse           (const P7_OPROFILE    *om,    const P7_OMX *pp,       P7_OMX *ox, float *ret_e);
extern int p7_OATrace_sse                   (const P7_OPROFILE    *om,    const P7_OMX *pp, const P7_OMX *ox, P7_TRACE *tr);
extern int p7_OptimalAccuracy_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,       P7_OMX *ox, float *ret_e);
extern int p7_OATrace_Frameshift_sse        (const P7_FS_OPROFILE *om_fs, const P7_OMX *pp, const P7_OMX *ox, P7_TRACE *tr);
#endif

/* ssvfilter.c — dispatch wrapper */
extern int p7_SSVFilter(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc);
/* ssvfilter_sse.c — SSE implementation */
#ifdef eslENABLE_SSE
extern int p7_SSVFilter_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc);
#endif

/* stotrace.c */
extern int p7_StochasticTrace(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *ox, P7_TRACE *tr);

/* stotrace_fs.c */
extern int p7_StochasticTrace_Frameshift(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, P7_TRACE *tr);

/* stotrace_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_StochasticTrace_sse          (ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_OPROFILE    *om,    const P7_OMX *ox, P7_TRACE *tr);
extern int p7_StochasticTrace_Frameshift_sse(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, P7_TRACE *tr);
#endif

/* vitfilter.c — dispatch wrappers */
extern int p7_ViterbiFilter     (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_ViterbiFilter_BATH(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *ssvdata, float filtersc, double P, P7_HMM_WINDOWLIST *windowlist, float *ret_sc);
/* vitfilter_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_ViterbiFilter_sse     (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_ViterbiFilter_BATH_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *ssvdata, float filtersc, double P, P7_HMM_WINDOWLIST *windowlist, float *ret_sc);
#endif

/* vitfilter_fs.c */
extern int p7_Viterbi_Frameshift      (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                  P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int p7_Viterbi_Frameshift_Trace(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, P7_TRACE *tr);
/* vitfilter_fs_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Viterbi_Frameshift_sse      (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                  P7_OMX *ox, P7_OIVX *ov, float *opt_sc);
extern int p7_Viterbi_Frameshift_Trace_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, P7_TRACE *tr);
#endif

/* viterbi.c */
extern int p7_Viterbi      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_Viterbi_Trace(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *ox, P7_TRACE *tr);
/* viterbi_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Viterbi_sse      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_Viterbi_Trace_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *ox, P7_TRACE *tr);
#endif

/* viterbi_sp.c */
extern int p7_Viterbi_Spliced    (const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_OMX *ox, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov, int i_start, int i_end, int min_intron, int global_start, int global_end);
extern int p7_Viterbi_SplicedTrace(const ESL_DSQ *sub_dsq, const P7_OMX *ox, const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end, int min_intron, float *vitsc);
/* viterbi_sp_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Viterbi_Spliced_sse    (const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_OMX *ox, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov, int i_start, int i_end, int min_intron, int global_start, int global_end);
extern int p7_Viterbi_SplicedTrace_sse(const ESL_DSQ *sub_dsq, const P7_OMX *ox, const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end, int min_intron, float *vitsc);
#endif

/* vitscore.c */
extern int p7_ViterbiScore(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
/* vitscore_sse.c — SSE implementation */
#ifdef eslENABLE_SSE
extern int p7_ViterbiScore_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
#endif


/*****************************************************************
 * 6. Runtime dispatch initializer
 *
 * Call impl_Init() once at startup (before any profile operations).
 * It detects available ISA features, sets SSE flush-zero mode, and
 * wires the function pointers to the fastest available implementations.
 *****************************************************************/
extern void impl_Init(void);


#endif /* P7_IMPL_AVX_INCLUDED */
