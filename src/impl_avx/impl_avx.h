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
#ifdef __SSE2__
#include "esl_sse.h"
#endif
#endif

#ifdef eslENABLE_AVX
#include <immintrin.h>   /* AVX, AVX2 */
#ifdef __AVX__
#include "esl_avx.h"
#endif
#endif

#ifdef eslENABLE_AVX512
#include <immintrin.h>   /* AVX-512 */
#ifdef __AVX512F__
#include "esl_avx512.h"
#endif
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

/* retrieve match odds ratio [k][x]: self-patching dispatcher in p7_oprofile.c */
extern float p7_oprofile_FGetEmission(const P7_OPROFILE *om, int k, int x);


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
 * Self-patching dispatcher in p7_oprofile.c.
 */
extern float p7_fs_oprofile_FGetEmission(const P7_FS_OPROFILE *om_fs, int k, int c);

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


/* P7_CHKOIVX: vectors saved at each checkpoint of a P7_CHKOMX
 *
 * Holds <allocCI> working int32 vector channels, and for each of
 * <allocK> checkpoints a copy of <allocCF> float channels (taken
 * from a P7_OIVX) and of the <allocCI> int channels. Vectors use the
 * same [channel][stripe] layout as P7_OIVX, so int32 lane j lines up
 * with float lane j. Each ISA has its own pointers and stripe count.
 */
typedef struct p7_chkoivx_s {
  int        allocM;       /* max model length currently allocated for          */
  int        allocCF;      /* float channels saved per checkpoint               */
  int        allocCI;      /* int channels, working and saved per checkpoint    */
  int        allocK;       /* number of checkpoints                             */
#ifdef eslENABLE_SSE
  __m128i  **ivx;          /* [allocCI][allocQ4] working int32 vectors          */
  void      *ivx_mem;      /* backing memory for ivx (+15 bytes for alignment)  */
  void      *ck_mem;       /* [allocK][allocCF+allocCI][allocQ4] saved vectors  */
  int        allocQ4;      /* p7O_NQF(allocM)                                   */
#endif
#ifdef eslENABLE_AVX
  __m256i  **ivx_avx;
  void      *ivx_mem_avx;
  void      *ck_mem_avx;
  int        allocQ4_avx;  /* p7O_NQF_AVX(allocM)                               */
#endif
} P7_CHKOIVX;


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

/* The OA matrix uses the frameshift (p7X_NSCELLS_FS) cell layout, so the
 * fill can write it over the posterior matrix. Each filled cell holds
 * the OA scores of M, I and D in the C1, C2 and D slots, the codon length
 * the trace takes for an M there in C3, and the posteriors pp_C0 and
 * pp_I, which the trace reports, in C0 and I. C4 and C5 keep pp_C4 and
 * pp_C5, which nothing reads after the fill.
 */
#define MMO_OA(dp,q)    MMO_FS(dp,q,p7X_FS_C1)
#define IMO_OA(dp,q)    MMO_FS(dp,q,p7X_FS_C2)
#define DMO_OA(dp,q)    DMO_FS(dp,q)
#define CODONO_OA(dp,q) MMO_FS(dp,q,p7X_FS_C3)


/* DP cell accessors: self-patching dispatchers in p7_omx.c. */
extern float (*p7_omx_FGetMDI)(const P7_OMX *ox, int s, int i, int k);
extern void  (*p7_omx_FSetMDI)(const P7_OMX *ox, int s, int i, int k, float val);


/* OATO_*: cells of a P7_OATMX (hmmer.h) in the striped layout, where
 * <ox->lanes> is the number of floats per vector. Cell k is bit
 * q*lanes + r of each plane row, q = (k-1) % Q and r = (k-1) / Q with
 * Q = p7OAT_NQ(M, lanes), the stripe count of the SIMD profile.
 */
#define OATO_BIT(ox,k)            ( (((k)-1) % p7OAT_NQ((ox)->M, (ox)->lanes)) * (ox)->lanes + ((k)-1) / p7OAT_NQ((ox)->M, (ox)->lanes) )
#define OATO_MPRED(ox,i,k)        p7_oatmx_GetField((ox), (i), OATO_BIT(ox,k), p7OAT_MP,  2)
#define OATO_IPRED(ox,i,k)        p7_oatmx_GetField((ox), (i), OATO_BIT(ox,k), p7OAT_IP,  1)
#define OATO_DPRED(ox,i,k)        p7_oatmx_GetField((ox), (i), OATO_BIT(ox,k), p7OAT_DP,  1)
#define OATO_PPM(ox,i,k)          p7_oatmx_GetField((ox), (i), OATO_BIT(ox,k), p7OAT_PPM, 4)
#define OATO_PPI(ox,i,k)          p7_oatmx_GetField((ox), (i), OATO_BIT(ox,k), p7OAT_PPI, 4)
#define OATO_CODON(ox,i,k)        p7_oatmx_GetField((ox), (i), OATO_BIT(ox,k), p7OAT_COD, 3)
#define OATO_SET_MPRED(ox,i,k,v)  p7_oatmx_SetField((ox), (i), OATO_BIT(ox,k), p7OAT_MP,  2, (v))
#define OATO_SET_IPRED(ox,i,k,v)  p7_oatmx_SetField((ox), (i), OATO_BIT(ox,k), p7OAT_IP,  1, (v))
#define OATO_SET_DPRED(ox,i,k,v)  p7_oatmx_SetField((ox), (i), OATO_BIT(ox,k), p7OAT_DP,  1, (v))
#define OATO_SET_PPM(ox,i,k,v)    p7_oatmx_SetField((ox), (i), OATO_BIT(ox,k), p7OAT_PPM, 4, (v))
#define OATO_SET_PPI(ox,i,k,v)    p7_oatmx_SetField((ox), (i), OATO_BIT(ox,k), p7OAT_PPI, 4, (v))
#define OATO_SET_CODON(ox,i,k,v)  p7_oatmx_SetField((ox), (i), OATO_BIT(ox,k), p7OAT_COD, 3, (v))

/* Writing one plane's bits for all lanes of stripe q at once, from a
 * lane mask <m> (bit r = lane r):
 */
/* AVX, 8 lanes: the stripe is byte q of the plane row. */
#define OATO_STRIPE8(ox,i,p,q)    ( p7OAT_ROW(ox,i,p)[q] )
/* 4 lanes (SSE): the stripe is half of byte q/2 of the plane row. */
#define OATO_SET_STRIPE4(ox,i,p,q,m) do {                                              \
    uint8_t *oato_c = p7OAT_ROW(ox,i,p) + ((q) >> 1);                                  \
    int      oato_s = ((q) & 1) * 4;                                                   \
    *oato_c = (uint8_t) ((*oato_c & ~(0xf << oato_s)) | (((m) & 0xf) << oato_s));      \
  } while (0)


/* P7_CHKOMX: checkpointed float dynamic programming matrix
 *
 * Rows are computed in blocks of <blockR> rows; the block starting at
 * row b covers rows b .. b+blockR-1, for b a multiple of blockR. Only
 * one block is stored at a time, together with <saveR> rows next to
 * it: the rows before it for a computation that runs forward through
 * the rows (dir = p7_CHK_FORWARD), the rows after it for one that runs
 * backward (dir = p7_CHK_BACKWARD). The row pointer arrays have an
 * entry for every row 0..L: loaded rows point into dp_mem, all other
 * rows are NULL. When the computation moves into the next block, the
 * <saveR> rows it needs from the block it leaves are copied into a
 * checkpoint in ck_mem, so that the block can be loaded again later.
 * X states are kept for every row.
 *
 * dp_mem holds <allocR> >= <saveR> + <blockR> rows: forward, slots
 * 0..saveR-1 for the rows before the loaded block, then the block;
 * backward, the block in slots 0..blockR-1, then the rows after it. Rows
 * are <allocQ4_avx> (or <allocQ4>) stripes of <nscells> vectors.
 */
typedef struct p7_chkomx_s {
  int       M;
  int       L;
  int       nscells;   /* vectors per stripe position, e.g. p7X_NSCELLS        */
  int       dir;       /* p7_CHK_FORWARD or p7_CHK_BACKWARD                    */
  int       blockR;    /* rows per block                                       */
  int       saveR;     /* rows saved at each checkpoint, next to the block       */
  int       b;         /* first row of the loaded block; -1 if none            */

  /* Row pointers [0..L] — only the active ISA's pointers are set */
#ifdef eslENABLE_SSE
  __m128  **dpf;
#endif
#ifdef eslENABLE_AVX
  __m256  **dpf_avx;
#endif

  void     *dp_mem;    /* storage for allocR rows                              */
  int       allocR;    /* rows of the current width that fit in dp_mem: >= saveR + blockR */
  int       allocP;    /* entries in the row pointer array; allocP >= L+1      */
  int       allocQ4;   /* SSE float stripes per row     (allocQ4*4 >= M)       */
  int       allocQ4_avx;
  size_t    nvec;      /* current allocation size of dp_mem, in vectors        */
  size_t    nckvec;    /* current allocation size of ck_mem, in vectors        */

  /* Checkpoints: saveR rows each */
  void     *ck_mem;
  int       allocK;    /* checkpoints of the current width that fit in ck_mem  */
  int       nK;        /* checkpoints saved                                    */

  /* X states for every row */
  float    *xmx;       /* [0..L][ENJBCS]; indexed [i*p7X_NXCELLS+s]            */
  void     *x_mem;
  int       allocXR;
} P7_CHKOMX;

/* Scale factor a SIMD Backward row takes (p7_BackwardRow() and the
 * frameshift Backward rows):
 *   p7_BWD_SCALE_FORWARD: the Forward matrix's factor for the row, or
 *     the row's own factor once its xB has exceeded 1e16 (full Backward);
 *   p7_BWD_SCALE_OWN:     the row's own factor, xB when xB > 1e4, else 1
 *     (checkpointed Backward);
 *   p7_BWD_SCALE_STORED:  the factor already stored for the row, when a
 *     row is computed again (checkpointed Backward Load).
 */
enum p7_bwdscale_e { p7_BWD_SCALE_FORWARD = 0, p7_BWD_SCALE_OWN = 1, p7_BWD_SCALE_STORED = 2 };


/*****************************************************************
 * 5. External API declarations
 *****************************************************************/

/* p7_omx.c — self-patching dispatchers and ISA-independent functions */
extern P7_OMX     *(*p7_omx_Create)    (int allocM, int allocL, int allocXL);
extern int         (*p7_omx_GrowTo)    (P7_OMX *ox, int allocM, int allocL, int allocXL);
extern P7_OMX     *(*p7_omx_Create_dpf)(int allocM, int allocL, int allocXL, int nscells);
extern int         (*p7_omx_GrowTo_dpf)(P7_OMX *ox, int allocM, int allocL, int allocXL);
extern int         (*p7_omx_FDeconvert)(P7_OMX *ox, P7_GMX *gx);
extern int          p7_omx_Reuse     (P7_OMX *ox);
extern void       (*p7_omx_Destroy)  (P7_OMX *ox);
extern int          p7_omx_SetDumpMode(FILE *fp, P7_OMX *ox, int truefalse);
extern int          p7_omx_Dump      (FILE *fp, P7_OMX *ox);
extern int          p7_omx_DumpMFRow (P7_OMX *ox, int rowi, uint8_t xE, uint8_t xN, uint8_t xJ, uint8_t xB, uint8_t xC);
extern int          p7_omx_DumpVFRow (P7_OMX *ox, int rowi, int16_t xE, int16_t xN, int16_t xJ, int16_t xB, int16_t xC);
extern int          p7_omx_DumpFBRow (P7_OMX *ox, int logify, int rowi, int width, int precision, float xE, float xN, float xJ, float xB, float xC);
extern int          p7_omx_DumpFBRow_FS(P7_OMX *ox, int logify, int i, int rowi, int width, int precision, float xE, float xN, float xJ, float xB, float xC);

/* p7_omx_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern P7_OMX      *p7_omx_Create_sse    (int allocM, int allocL, int allocXL);
extern int          p7_omx_GrowTo_sse    (P7_OMX *ox, int allocM, int allocL, int allocXL);
extern P7_OMX      *p7_omx_Create_dpf_sse(int allocM, int allocL, int allocXL, int nscells);
extern int          p7_omx_GrowTo_dpf_sse(P7_OMX *ox, int allocM, int allocL, int allocXL);
extern int          p7_omx_FDeconvert_sse(P7_OMX *ox, P7_GMX *gx);
extern float        p7_omx_FGetMDI_sse   (const P7_OMX *ox, int s, int i, int k);
extern void         p7_omx_Destroy_sse   (P7_OMX *ox);
#endif
/* p7_omx_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern P7_OMX      *p7_omx_Create_avx    (int allocM, int allocL, int allocXL);
extern int          p7_omx_GrowTo_avx    (P7_OMX *ox, int allocM, int allocL, int allocXL);
extern P7_OMX      *p7_omx_Create_dpf_avx(int allocM, int allocL, int allocXL, int nscells);
extern int          p7_omx_GrowTo_dpf_avx(P7_OMX *ox, int allocM, int allocL, int allocXL);
extern void         p7_omx_Destroy_avx   (P7_OMX *ox);
extern int          p7_omx_FDeconvert_avx(P7_OMX *ox, P7_GMX *gx);
#endif

/* p7_chkomx.c — self-patching dispatchers */
extern P7_CHKOMX *(*p7_chkomx_Create)    (int allocM, int allocL, int blockR, int saveR, int nscells, int dir);
extern int        (*p7_chkomx_GrowTo)    (P7_CHKOMX *cx, int allocM, int allocL, int blockR);
extern int        (*p7_chkomx_StartBlock)(P7_CHKOMX *cx, int b);
extern int        (*p7_chkomx_LoadBlock) (P7_CHKOMX *cx, int b);
extern float      (*p7_chkomx_FGetMDI)   (const P7_CHKOMX *cx, int s, int i, int k);
extern void       (*p7_chkomx_Destroy)   (P7_CHKOMX *cx);
extern int          p7_chkomx_SetDirection(P7_CHKOMX *cx, int dir, int saveR);

/* p7_chkomx_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern P7_CHKOMX   *p7_chkomx_Create_sse    (int allocM, int allocL, int blockR, int saveR, int nscells, int dir);
extern int          p7_chkomx_GrowTo_sse    (P7_CHKOMX *cx, int allocM, int allocL, int blockR);
extern int          p7_chkomx_StartBlock_sse(P7_CHKOMX *cx, int b);
extern int          p7_chkomx_LoadBlock_sse (P7_CHKOMX *cx, int b);
extern float        p7_chkomx_FGetMDI_sse   (const P7_CHKOMX *cx, int s, int i, int k);
extern void         p7_chkomx_Destroy_sse   (P7_CHKOMX *cx);
#endif

/* p7_chkomx_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern P7_CHKOMX   *p7_chkomx_Create_avx    (int allocM, int allocL, int blockR, int saveR, int nscells, int dir);
extern int          p7_chkomx_GrowTo_avx    (P7_CHKOMX *cx, int allocM, int allocL, int blockR);
extern int          p7_chkomx_StartBlock_avx(P7_CHKOMX *cx, int b);
extern int          p7_chkomx_LoadBlock_avx (P7_CHKOMX *cx, int b);
extern float        p7_chkomx_FGetMDI_avx   (const P7_CHKOMX *cx, int s, int i, int k);
extern void         p7_chkomx_Destroy_avx   (P7_CHKOMX *cx);
#endif

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
extern int          p7_oprofile_GetFwdTransitionArray_sse(const P7_OPROFILE *om, int type, float *arr);
extern int          p7_oprofile_GetSSVEmissionScoreArray_sse(const P7_OPROFILE *om, uint8_t *arr);
extern int          p7_oprofile_GetFwdEmissionScoreArray_sse(const P7_OPROFILE *om, float *arr);
extern int          p7_oprofile_GetFwdEmissionArray_sse(const P7_OPROFILE *om, P7_BG *bg, float *arr);
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
extern int          p7_oprofile_Logify_avx(P7_OPROFILE *om);
extern int          p7_oprofile_GetFwdTransitionArray_avx(const P7_OPROFILE *om, int type, float *arr);
extern int          p7_oprofile_GetSSVEmissionScoreArray_avx(const P7_OPROFILE *om, uint8_t *arr);
extern int          p7_oprofile_GetFwdEmissionScoreArray_avx(const P7_OPROFILE *om, float *arr);
extern int          p7_oprofile_GetFwdEmissionArray_avx(const P7_OPROFILE *om, P7_BG *bg, float *arr);
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

/* p7_fs_oprofile_avx.c — AVX2 implementations */
extern P7_FS_OPROFILE *p7_fs_oprofile_Create_avx(int allocM, const ESL_ALPHABET *abc, int codon_lengths);
extern void            p7_fs_oprofile_Destroy_avx(P7_FS_OPROFILE *om_fs);
extern P7_FS_OPROFILE *p7_fs_oprofile_Clone_avx(const P7_FS_OPROFILE *om_fs);
extern int             p7_fs_oprofile_Convert_avx(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs);
extern int             p7_fs_oprofile_Convert_Log_avx(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs);
extern int             p7_fs_oprofile_SubConvert_Log_avx(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs, int k_start, int k_end);
extern int             p7_fs_oprofile_ReconfigLength_avx(P7_FS_OPROFILE *om_fs, int L);
extern int             p7_fs_oprofile_ReconfigLength_Log_avx(P7_FS_OPROFILE *om_fs, int L);
extern int             p7_fs_oprofile_ReconfigMultihit_avx(P7_FS_OPROFILE *om_fs, int L);
extern int             p7_fs_oprofile_ReconfigUnihit_avx(P7_FS_OPROFILE *om_fs, int L);
extern int             p7_fs_oprofile_Logify_avx(P7_FS_OPROFILE *om_fs);

/* p7_fs_oprofile_sse.c — SSE implementations (used as SSE fallback in dispatch) */
extern P7_FS_OPROFILE *p7_fs_oprofile_Create_sse(int allocM, const ESL_ALPHABET *abc, int codon_lengths);
extern void            p7_fs_oprofile_Destroy_sse(P7_FS_OPROFILE *om_fs);
extern P7_FS_OPROFILE *p7_fs_oprofile_Clone_sse(const P7_FS_OPROFILE *om_fs);
extern int             p7_fs_oprofile_Convert_sse(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs);
extern int             p7_fs_oprofile_Convert_Log_sse(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs);
extern int             p7_fs_oprofile_SubConvert_Log_sse(const P7_FS_PROFILE *gm_fs, P7_FS_OPROFILE *om_fs, int k_start, int k_end);
extern int             p7_fs_oprofile_ReconfigLength_sse(P7_FS_OPROFILE *om_fs, int L);
extern int             p7_fs_oprofile_ReconfigLength_Log_sse(P7_FS_OPROFILE *om_fs, int L);
extern int             p7_fs_oprofile_ReconfigMultihit_sse(P7_FS_OPROFILE *om_fs, int L);
extern int             p7_fs_oprofile_ReconfigUnihit_sse(P7_FS_OPROFILE *om_fs, int L);
extern int             p7_fs_oprofile_Logify_sse(P7_FS_OPROFILE *om_fs);

/* p7_oivx.c — function pointers (dispatched at runtime) */
extern P7_OIVX *(*p7_oivx_Create)(int M_hint, int C);
extern int      (*p7_oivx_GrowTo)(P7_OIVX *ov, int M, int C);
extern void     (*p7_oivx_Destroy)(P7_OIVX *ov);

/* p7_oivx_avx.c — AVX2 implementations */
extern P7_OIVX *p7_oivx_Create_avx(int M_hint, int C);
extern int      p7_oivx_GrowTo_avx(P7_OIVX *ov, int M, int C);
extern void     p7_oivx_Destroy_avx(P7_OIVX *ov);

/* p7_oivx_sse.c — SSE implementations (used as SSE fallback in dispatch) */
extern P7_OIVX *p7_oivx_Create_sse(int M_hint, int C);
extern int      p7_oivx_GrowTo_sse(P7_OIVX *ov, int M, int C);
extern void     p7_oivx_Destroy_sse(P7_OIVX *ov);

/* p7_chkoivx.c — function pointers (dispatched at runtime) */
extern P7_CHKOIVX *(*p7_chkoivx_Create) (int M_hint, int CF, int CI, int K);
extern int         (*p7_chkoivx_GrowTo) (P7_CHKOIVX *cv, int M, int CF, int CI, int K);
extern int         (*p7_chkoivx_Save)   (P7_CHKOIVX *cv, const P7_OIVX *ov, int M, int k);
extern int         (*p7_chkoivx_Load)   (P7_CHKOIVX *cv, P7_OIVX *ov, int M, int k);
extern void        (*p7_chkoivx_Destroy)(P7_CHKOIVX *cv);

/* p7_chkoivx_avx.c — AVX2 implementations */
extern P7_CHKOIVX *p7_chkoivx_Create_avx (int M_hint, int CF, int CI, int K);
extern int         p7_chkoivx_GrowTo_avx (P7_CHKOIVX *cv, int M, int CF, int CI, int K);
extern int         p7_chkoivx_Save_avx   (P7_CHKOIVX *cv, const P7_OIVX *ov, int M, int k);
extern int         p7_chkoivx_Load_avx   (P7_CHKOIVX *cv, P7_OIVX *ov, int M, int k);
extern void        p7_chkoivx_Destroy_avx(P7_CHKOIVX *cv);

/* p7_chkoivx_sse.c — SSE implementations (used as SSE fallback in dispatch) */
extern P7_CHKOIVX *p7_chkoivx_Create_sse (int M_hint, int CF, int CI, int K);
extern int         p7_chkoivx_GrowTo_sse (P7_CHKOIVX *cv, int M, int CF, int CI, int K);
extern int         p7_chkoivx_Save_sse   (P7_CHKOIVX *cv, const P7_OIVX *ov, int M, int k);
extern int         p7_chkoivx_Load_sse   (P7_CHKOIVX *cv, P7_OIVX *ov, int M, int k);
extern void        p7_chkoivx_Destroy_sse(P7_CHKOIVX *cv);

/* decoding.c — self-patching dispatchers */
extern int (*p7_Decoding)      (const P7_OPROFILE *om, const P7_OMX *oxf,       P7_OMX *oxb, P7_OMX *pp);
extern int (*p7_DomainDecoding)(const P7_OPROFILE *om, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef);

/* decoding_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Decoding_sse               (const P7_OPROFILE *om,    const P7_OMX *oxf,       P7_OMX *oxb,       P7_OMX *pp);
extern void p7_DecodingRow_sse           (const P7_OPROFILE *om, const float *fxmx, const float *bxmx, const __m128 *fv, const __m128 *bv, __m128 *ppv, float *ppx, int i, float scaleproduct);
extern int p7_DomainDecoding_sse         (const P7_OPROFILE *om,    const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef);
#endif
/* decoding_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_Decoding_avx               (const P7_OPROFILE *om,    const P7_OMX *oxf,       P7_OMX *oxb,       P7_OMX *pp);
extern void p7_DecodingRow_avx           (const P7_OPROFILE *om, const float *fxmx, const float *bxmx, const __m256 *fv, const __m256 *bv, __m256 *ppv, float *ppx, int i, float scaleproduct);
extern int p7_DomainDecoding_avx         (const P7_OPROFILE *om,    const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef);
#endif

/* decoding_fs.c — self-patching dispatchers */
extern int (*p7_Decoding_Frameshift)   (const P7_FS_OPROFILE *om_fs, P7_OMX *fwd, const P7_OMX *bck);
extern int (*p7_DomainDecoding_Frameshift)(const P7_FS_OPROFILE *om_fs, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef);

/* decoding_fs_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Decoding_Frameshift_sse    (const P7_FS_OPROFILE *om_fs,                   P7_OMX *fwd, const P7_OMX *bck);
extern int p7_DecodingRow_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const float *bxmx, const float *fx3, const float *log_sfwd, const float *log_sbck, float log_inv_Z, const __m128 *fv, const __m128 *bv, __m128 *ppv, float *ppx, int i);
extern int p7_DomainDecoding_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef);
#endif
/* decoding_fs_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_Decoding_Frameshift_avx    (const P7_FS_OPROFILE *om_fs,                   P7_OMX *fwd, const P7_OMX *bck);
extern int p7_DecodingRow_Frameshift_avx(const P7_FS_OPROFILE *om_fs, const float *bxmx, const float *fx3, const float *log_sfwd, const float *log_sbck, float log_inv_Z, const __m256 *fv, const __m256 *bv, __m256 *ppv, float *ppx, int i);
extern int p7_DomainDecoding_Frameshift_avx(const P7_FS_OPROFILE *om_fs, const P7_OMX *oxf, const P7_OMX *oxb, P7_DOMAINDEF *ddef);
#endif



/* fwdback.c — self-patching dispatchers */
extern int (*p7_Forward)         (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
extern int (*p7_ForwardParser)   (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
extern int (*p7_Backward)        (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);
extern int (*p7_BackwardParser)  (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);
extern int (*p7_BackwardChk)     (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc);
extern int (*p7_BackwardChkLoad) (const ESL_DSQ *dsq,        const P7_OPROFILE *om, P7_CHKOMX *cx, int b);
extern int (*p7_ForwardChk)      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc);
extern int (*p7_ForwardChkLoad)  (const ESL_DSQ *dsq,        const P7_OPROFILE *om, P7_CHKOMX *cx, int b);
/* fwdback_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Forward_sse         (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
extern int p7_ForwardParser_sse   (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
extern int p7_Backward_sse        (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);
extern int p7_BackwardParser_sse  (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);
extern void p7_ForwardRow_sse      (const ESL_DSQ *dsq, const P7_OPROFILE *om, float *xmx, float *totscale, __m128 *dpc, const __m128 *dpp, int i);
extern float p7_BackwardRow_sse    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const float *fwd_xmx, float *xmx, float *totscale, int *has_own_scales, __m128 *dpc, const __m128 *dpp, int i, int scale_mode);
extern int p7_BackwardChk_sse     (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc);
extern int p7_BackwardChkLoad_sse (const ESL_DSQ *dsq,        const P7_OPROFILE *om, P7_CHKOMX *cx, int b);
extern int p7_ForwardChk_sse      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc);
extern int p7_ForwardChkLoad_sse  (const ESL_DSQ *dsq,        const P7_OPROFILE *om, P7_CHKOMX *cx, int b);
#endif

/* fwdback_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_Forward_avx         (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
extern int p7_ForwardParser_avx   (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om,                    P7_OMX *fwd, float *opt_sc);
extern int p7_Backward_avx        (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);
extern int p7_BackwardParser_avx  (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *fwd, P7_OMX *bck, float *opt_sc);
extern void p7_ForwardRow_avx      (const ESL_DSQ *dsq, const P7_OPROFILE *om, float *xmx, float *totscale, __m256 *dpc, const __m256 *dpp, int i);
extern float p7_BackwardRow_avx    (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const float *fwd_xmx, float *xmx, float *totscale, int *has_own_scales, __m256 *dpc, const __m256 *dpp, int i, int scale_mode);
extern int p7_BackwardChk_avx     (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc);
extern int p7_BackwardChkLoad_avx (const ESL_DSQ *dsq,        const P7_OPROFILE *om, P7_CHKOMX *cx, int b);
extern int p7_ForwardChk_avx      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, float *opt_sc);
extern int p7_ForwardChkLoad_avx  (const ESL_DSQ *dsq,        const P7_OPROFILE *om, P7_CHKOMX *cx, int b);
#endif

/* fwdback_fs.c — self-patching dispatchers */
extern int (*p7_ForwardParser_Frameshift_3Codons) (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int (*p7_BackwardParser_Frameshift_3Codons)(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);
extern int (*p7_ForwardParser_Frameshift_5Codons) (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int (*p7_BackwardParser_Frameshift_5Codons)(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);
extern int (*p7_Forward_Frameshift)               (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int (*p7_Backward_Frameshift)              (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);
extern int (*p7_ForwardChk_Frameshift)            (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, float *opt_sc);
extern int (*p7_ForwardChkLoad_Frameshift)        (const ESL_DSQ *dsq,        const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, int b);
extern int (*p7_BackwardChk_Frameshift)           (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, float *opt_sc);
extern int (*p7_BackwardChkLoad_Frameshift)       (const ESL_DSQ *dsq,        const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, int b);

/* fwdback_fs_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_ForwardParser_Frameshift_3Codons_sse (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int p7_BackwardParser_Frameshift_3Codons_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);
extern int p7_ForwardParser_Frameshift_5Codons_sse (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int p7_BackwardParser_Frameshift_5Codons_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);
extern int p7_Forward_Frameshift_sse               (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern void p7_ForwardRow_Frameshift_sse         (const ESL_DSQ *dsq, const P7_FS_OPROFILE *om_fs, float *xmx, float *totscale, P7_FS_FWDRING *rg, __m128 **ivx, __m128 *dpc, const __m128 *dpp1, const __m128 *dpp3, int i);
extern void p7_ForwardRow3_Frameshift_sse        (const ESL_DSQ *dsq, const P7_FS_OPROFILE *om_fs, float *xmx, float *totscale, P7_FS_FWDRING *rg, __m128 **ivx, __m128 *dpc, const __m128 *dpp1, const __m128 *dpp3, int i);
extern void p7_BackwardRowL_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const float *fwd_xmx, float *xmx, float *totscale, P7_FS_BWDRING *rg, __m128 *dpc, int L, int scale_mode);
extern void p7_BackwardRow_Frameshift_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const float *fwd_xmx, float *xmx, float *totscale, int *has_own_scales, P7_FS_BWDRING *rg, __m128 *const *dp, __m128 *ivxf, int i, int scale_mode);
extern float p7_BackwardRow0_Frameshift_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, float *xmx, const P7_FS_BWDRING *rg, __m128 *const *dp, __m128 *ivxf);
extern int p7_Backward_Frameshift_sse              (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);
#endif
/* fwdback_chk_fs_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_ForwardChk_Frameshift_sse            (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, float *opt_sc);
extern int p7_ForwardChkLoad_Frameshift_sse        (const ESL_DSQ *dsq,        const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, int b);
extern int p7_BackwardChk_Frameshift_sse           (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, float *opt_sc);
extern int p7_BackwardChkLoad_Frameshift_sse       (const ESL_DSQ *dsq,        const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, int b);
#endif

/* fwdback_fs_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_ForwardParser_Frameshift_3Codons_avx (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int p7_BackwardParser_Frameshift_3Codons_avx(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);
extern int p7_ForwardParser_Frameshift_5Codons_avx (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int p7_BackwardParser_Frameshift_5Codons_avx(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);
extern int p7_Forward_Frameshift_avx               (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                    P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern void p7_ForwardRow_Frameshift_avx         (const ESL_DSQ *dsq, const P7_FS_OPROFILE *om_fs, float *xmx, float *totscale, P7_FS_FWDRING *rg, __m256 **ivx, __m256 *dpc, const __m256 *dpp1, const __m256 *dpp3, int i);
extern void p7_ForwardRow3_Frameshift_avx        (const ESL_DSQ *dsq, const P7_FS_OPROFILE *om_fs, float *xmx, float *totscale, P7_FS_FWDRING *rg, __m256 **ivx, __m256 *dpc, const __m256 *dpp1, const __m256 *dpp3, int i);
extern void p7_BackwardRowL_Frameshift_avx(const P7_FS_OPROFILE *om_fs, const float *fwd_xmx, float *xmx, float *totscale, P7_FS_BWDRING *rg, __m256 *dpc, int L, int scale_mode);
extern void p7_BackwardRow_Frameshift_avx(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const float *fwd_xmx, float *xmx, float *totscale, int *has_own_scales, P7_FS_BWDRING *rg, __m256 *const *dp, __m256 *ivxf, int i, int scale_mode);
extern float p7_BackwardRow0_Frameshift_avx(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, float *xmx, const P7_FS_BWDRING *rg, __m256 *const *dp, __m256 *ivxf);
extern int p7_Backward_Frameshift_avx              (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *fwd, P7_OMX *bck, P7_OIVX *ov, float *opt_sc);
#endif
/* fwdback_chk_fs_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_ForwardChk_Frameshift_avx            (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, float *opt_sc);
extern int p7_ForwardChkLoad_Frameshift_avx        (const ESL_DSQ *dsq,        const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, int b);
extern int p7_BackwardChk_Frameshift_avx           (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, float *opt_sc);
extern int p7_BackwardChkLoad_Frameshift_avx       (const ESL_DSQ *dsq,        const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, int b);
#endif

/* p7_oprofile.c — block management (ISA-independent) */
extern P7_OM_BLOCK *p7_oprofile_CreateBlock(int size);
extern void         p7_oprofile_DestroyBlock(P7_OM_BLOCK *block);

/* msvfilter.c — self-patching dispatchers */
extern int (*p7_MSVFilter)      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int (*p7_SSVFilter_BATH) (const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *msvdata, P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist);
/* msvfilter_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_MSVFilter_sse      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_SSVFilter_BATH_sse (const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *ssvdata, P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist);
#endif
/* msvfilter_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_MSVFilter_avx      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_SSVFilter_BATH_avx (const ESL_DSQ *dsq, int L, P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *ssvdata, P7_BG *bg, double P, P7_HMM_WINDOWLIST *windowlist);
#endif

/* null2.c — self-patching dispatchers */
extern int (*p7_Null2_ByExpectation)(const P7_OPROFILE *om, const P7_OMX *pp, float *null2);

/* null2_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Null2_ByExpectation_sse    (const P7_OPROFILE    *om,    const P7_OMX *pp, float *null2);
extern int p7_Null2_FromCounts_sse      (const P7_OPROFILE    *om,    __m128 *cnt, float *cntx, int Ld, float *null2);
#endif
/* null2_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_Null2_ByExpectation_avx   (const P7_OPROFILE    *om,    const P7_OMX *pp, float *null2);
extern int p7_Null2_FromCounts_avx      (const P7_OPROFILE    *om,    __m256 *cnt, float *cntx, int Ld, float *null2);
#endif

/* null2_fs.c — self-patching dispatcher */
extern int (*p7_Null2_fs_ByExpectation)(const P7_FS_OPROFILE *om_fs, P7_OMX *pp, float *null2);

/* null2_fs_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Null2_fs_ByExpectation_sse (const P7_FS_OPROFILE *om_fs, P7_OMX *pp, float *null2);
extern int p7_Null2_fs_FromCounts_sse(const P7_FS_OPROFILE *om_fs, __m128 *cnt, float *cntx, float *null2);
#endif
/* null2_fs_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_Null2_fs_ByExpectation_avx(const P7_FS_OPROFILE *om_fs, P7_OMX *pp, float *null2);
extern int p7_Null2_fs_FromCounts_avx(const P7_FS_OPROFILE *om_fs, __m256 *cnt, float *cntx, float *null2);
#endif

/* optacc.c — self-patching dispatchers */
extern int (*p7_OptimalAccuracy)(const P7_OPROFILE *om, const P7_OMX *pp,       P7_OMX *ox, float *ret_e);
extern int (*p7_OATrace)        (const P7_OPROFILE *om,                   const P7_OMX *ox, P7_TRACE *tr);
extern int (*p7_OptimalAccuracySweep)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, P7_OATMX *ox, float *opt_null2, float *opt_fsc, float *ret_oasc);
extern int (*p7_OATraceSweep)        (const P7_OPROFILE *om, const P7_OATMX *ox, P7_TRACE *tr);

/* optacc_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern void p7_OptimalAccuracyRow_sse       (const P7_OPROFILE *om, const __m128 *ppp, const float *ppx, float *xmx, __m128 *dpc, const __m128 *dpp, int i);
extern int p7_OptimalAccuracy_sse           (const P7_OPROFILE    *om,    const P7_OMX *pp,       P7_OMX *ox, float *ret_e);
extern int p7_OptimalAccuracySweep_sse      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, P7_OATMX *ox, float *opt_null2, float *opt_fsc, float *ret_oasc);
extern int p7_OATraceSweep_sse              (const P7_OPROFILE    *om,    const P7_OATMX *ox, P7_TRACE *tr);
extern int p7_OATrace_sse                   (const P7_OPROFILE    *om,                      const P7_OMX *ox, P7_TRACE *tr);
#endif
/* optacc_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern void p7_OptimalAccuracyRow_avx       (const P7_OPROFILE *om, const __m256 *ppp, const float *ppx, float *xmx, __m256 *dpc, const __m256 *dpp, int i);
extern int p7_OptimalAccuracy_avx           (const P7_OPROFILE    *om,    const P7_OMX *pp,       P7_OMX *ox, float *ret_e);
extern int p7_OptimalAccuracySweep_avx      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, P7_OATMX *ox, float *opt_null2, float *opt_fsc, float *ret_oasc);
extern int p7_OATraceSweep_avx              (const P7_OPROFILE    *om,    const P7_OATMX *ox, P7_TRACE *tr);
extern int p7_OATrace_avx                   (const P7_OPROFILE    *om,                      const P7_OMX *ox, P7_TRACE *tr);
#endif

/* optacc_fs.c — self-patching dispatchers */
extern int (*p7_OptimalAccuracy_Frameshift)(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp, P7_OMX *ox, float *ret_e);
extern int (*p7_OATrace_Frameshift)        (const P7_FS_OPROFILE *om_fs, const float *ppx, const P7_OMX *ox, P7_TRACE *tr);
extern int (*p7_OptimalAccuracySweep_Frameshift)(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, P7_OATMX *ox, float *opt_null2, float *opt_fsc, float *ret_oasc);
extern int (*p7_OATraceSweep_Frameshift)        (const P7_FS_OPROFILE *om_fs, const P7_OATMX *ox, P7_TRACE *tr);

/* optacc_fs_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_OptimalAccuracy_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,       P7_OMX *ox, float *ret_e);
extern int p7_OATrace_Frameshift_sse        (const P7_FS_OPROFILE *om_fs, const float *ppx, const P7_OMX *ox, P7_TRACE *tr);
extern void p7_OptimalAccuracyRow_Frameshift_sse(const P7_FS_OPROFILE *om_fs, const __m128 *ppc, const float *ppx, float *xmx, __m128 *dpc, const __m128 *dp1, const __m128 *dp2, const __m128 *dp3, const __m128 *dp4, const __m128 *dp5, int i);
#endif
/* optacc_sweep_fs_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_OptimalAccuracySweep_Frameshift_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, P7_OATMX *ox, float *opt_null2, float *opt_fsc, float *ret_oasc);
extern int p7_OATraceSweep_Frameshift_sse        (const P7_FS_OPROFILE *om_fs, const P7_OATMX *ox, P7_TRACE *tr);
#endif
/* optacc_fs_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_OptimalAccuracy_Frameshift_avx(const P7_FS_OPROFILE *om_fs, const P7_OMX *pp,       P7_OMX *ox, float *ret_e);
extern int p7_OATrace_Frameshift_avx        (const P7_FS_OPROFILE *om_fs, const float *ppx, const P7_OMX *ox, P7_TRACE *tr);
extern void p7_OptimalAccuracyRow_Frameshift_avx(const P7_FS_OPROFILE *om_fs, const __m256 *ppc, const float *ppx, float *xmx, __m256 *dpc, const __m256 *dp1, const __m256 *dp2, const __m256 *dp3, const __m256 *dp4, const __m256 *dp5, int i);
#endif
/* optacc_sweep_fs_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_OptimalAccuracySweep_Frameshift_avx(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, P7_OATMX *ox, float *opt_null2, float *opt_fsc, float *ret_oasc);
extern int p7_OATraceSweep_Frameshift_avx        (const P7_FS_OPROFILE *om_fs, const P7_OATMX *ox, P7_TRACE *tr);
#endif


/* ssvfilter.c — self-patching dispatcher */
extern int (*p7_SSVFilter)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc);
/* ssvfilter_sse.c — SSE implementation */
#ifdef eslENABLE_SSE
extern int p7_SSVFilter_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc);
#endif
/* ssvfilter_avx.c — AVX2 implementation */
#ifdef eslENABLE_AVX
extern int p7_SSVFilter_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, float *ret_sc);
#endif

/* stotrace.c — self-patching dispatcher */
extern int (*p7_StochasticTrace)(ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *ox, P7_TRACE *tr);
extern int (*p7_StochasticTraceChk)(ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, P7_TRACE **tr);

/* stotrace_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_StochasticTrace_sse          (ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_OPROFILE    *om,    const P7_OMX *ox, P7_TRACE *tr);
extern int p7_StochasticTraceChk_sse       (ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, P7_TRACE **tr);
#endif
/* stotrace_avx.c — AVX2 implementation */
#ifdef eslENABLE_AVX
extern int p7_StochasticTrace_avx           (ESL_RANDOMNESS *rng, const ESL_DSQ *dsq, int L, const P7_OPROFILE    *om,    const P7_OMX *ox, P7_TRACE *tr);
extern int p7_StochasticTraceChk_avx        (ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_CHKOMX *cx, P7_TRACE **tr);
#endif

/* stotrace_fs.c — self-patching dispatcher */
extern int (*p7_StochasticTraceChk_Frameshift)(ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, P7_TRACE **tr);

/* stotrace_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_StochasticTraceChk_Frameshift_sse(ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, P7_TRACE **tr);
#endif
/* stotrace_avx.c — AVX2 implementation */
#ifdef eslENABLE_AVX
extern int p7_StochasticTraceChk_Frameshift_avx(ESL_RANDOMNESS **rng, int ntr, const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, P7_CHKOMX *cx, P7_OIVX *ov, P7_TRACE **tr);
#endif

/* vitfilter.c — self-patching dispatchers */
extern int (*p7_ViterbiFilter)     (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int (*p7_ViterbiFilter_BATH)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *ssvdata, float filtersc, double P, P7_HMM_WINDOWLIST *windowlist, float *ret_sc);
/* vitfilter_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_ViterbiFilter_sse     (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_ViterbiFilter_BATH_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *ssvdata, float filtersc, double P, P7_HMM_WINDOWLIST *windowlist, float *ret_sc);
#endif
/* vitfilter_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_ViterbiFilter_avx     (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_ViterbiFilter_BATH_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, const P7_SCOREDATA *ssvdata, float filtersc, double P, P7_HMM_WINDOWLIST *windowlist, float *ret_sc);
#endif

/* viterbi_fs.c — self-patching dispatchers */
extern int (*p7_Viterbi_Frameshift)      (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                  P7_OMX *ox,  P7_OIVX *ov, float *opt_sc);
extern int (*p7_Viterbi_Frameshift_Trace)(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, P7_TRACE *tr);
/* viterbi_fs_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Viterbi_Frameshift_sse      (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                  P7_OMX *ox, P7_OIVX *ov, float *opt_sc);
extern int p7_Viterbi_Frameshift_Trace_sse(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, P7_TRACE *tr);
#endif
/* viterbi_fs_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_Viterbi_Frameshift_avx      (const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs,                  P7_OMX *ox, P7_OIVX *ov, float *opt_sc);
extern int p7_Viterbi_Frameshift_Trace_avx(const ESL_DSQ *dsq, int L, const P7_FS_OPROFILE *om_fs, const P7_OMX *ox, P7_TRACE *tr);
#endif

/* viterbi.c — self-patching dispatchers */
extern int (*p7_Viterbi)      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int (*p7_Viterbi_Trace)(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *ox, P7_TRACE *tr);
/* viterbi_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Viterbi_sse      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_Viterbi_Trace_sse(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *ox, P7_TRACE *tr);
#endif
/* viterbi_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_Viterbi_avx      (const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, P7_OMX *ox, float *ret_sc);
extern int p7_Viterbi_Trace_avx(const ESL_DSQ *dsq, int L, const P7_OPROFILE *om, const P7_OMX *ox, P7_TRACE *tr);
#endif

/* viterbi_sp.c — self-patching dispatchers */
extern int (*p7_Viterbi_Spliced)    (const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_OMX *ox, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov, int i_start, int i_end, int min_intron, int global_start, int global_end);
extern int (*p7_Viterbi_SplicedTrace)(const ESL_DSQ *sub_dsq, const P7_OMX *ox, const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end, int min_intron, float *vitsc, float *max_drop);
/* viterbi_sp_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Viterbi_Spliced_sse    (const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_OMX *ox, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov, int i_start, int i_end, int min_intron, int global_start, int global_end);
extern int p7_Viterbi_SplicedTrace_sse(const ESL_DSQ *sub_dsq, const P7_OMX *ox, const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end, int min_intron, float *vitsc, float *max_drop);
#endif
/* viterbi_sp_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_Viterbi_Spliced_avx    (const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_OMX *ox, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov, int i_start, int i_end, int min_intron, int global_start, int global_end);
extern int p7_Viterbi_SplicedTrace_avx(const ESL_DSQ *sub_dsq, const P7_OMX *ox, const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end, int min_intron, float *vitsc, float *max_drop);
#endif

/* viterbi_sp_chk.c — self-patching dispatchers */
extern int (*p7_Viterbi_SplicedChk)     (const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov, int i_start, int i_end, int min_intron, int global_start, int global_end);
extern int (*p7_Viterbi_SplicedTraceChk)(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv, const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov, P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end, int min_intron, int global_start, int global_end, float *vitsc, float *max_drop);
/* viterbi_sp_chk_sse.c — SSE implementations */
#ifdef eslENABLE_SSE
extern int p7_Viterbi_SplicedChk_sse     (const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov, int i_start, int i_end, int min_intron, int global_start, int global_end);
extern int p7_Viterbi_SplicedTraceChk_sse(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv, const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov, P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end, int min_intron, int global_start, int global_end, float *vitsc, float *max_drop);
#endif
/* viterbi_sp_chk_avx.c — AVX2 implementations */
#ifdef eslENABLE_AVX
extern int p7_Viterbi_SplicedChk_avx     (const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov, int i_start, int i_end, int min_intron, int global_start, int global_end);
extern int p7_Viterbi_SplicedTraceChk_avx(const ESL_DSQ *sub_dsq, const P7_FS_OPROFILE *om_tr, P7_CHKOMX *cx, P7_CHKOIVX *cv, const P7_FS_PROFILE *gm_tr, const float *signal_scores, P7_OIVX *acc_ov, P7_OIVX *don_ov, P7_TRACE *tr, int i_start, int i_end, int k_start, int k_end, int min_intron, int global_start, int global_end, float *vitsc, float *max_drop);
#endif


/*****************************************************************
 * 6. Runtime dispatch initializer
 *
 * Call impl_Init() once at startup (before any profile operations).
 * It detects available ISA features, sets SSE flush-zero mode, and
 * wires the function pointers to the fastest available implementations.
 *****************************************************************/
static inline void
impl_Init(void)
{
#ifdef HAVE_FLUSH_ZERO_MODE
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif
#ifdef _PMMINTRIN_H_INCLUDED
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
}


#endif /* P7_IMPL_AVX_INCLUDED */
