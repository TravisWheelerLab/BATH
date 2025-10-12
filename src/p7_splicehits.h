/* Structs and MACROs for bilding a splice graph */

#include "p7_config.h"


typedef struct _splice_hit_info
{

  int64_t  seqidx;
  int      strand; /* 0 = forward strand, 1 = reverse strand */
  int      duplicate;
  int      is_seed;
  int      viterbi;

  int hmm_start;
  int hmm_end;
  int64_t seq_start;
  int64_t seq_end;

  uint64_t hit_id;  /* anchor hit for which splice hit is seed */ 

} SPLICE_HIT_INFO;

typedef struct _splice_saved_hits
{
  SPLICE_HIT_INFO **srt;
  SPLICE_HIT_INFO *unsrt;

  uint64_t Nalloc;
  uint64_t N;

  int is_sorted;

} SPLICE_SAVED_HITS;



