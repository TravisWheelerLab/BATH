/* Structs and MACROs for bilding a splice graph */

#include "p7_config.h"


typedef struct _splice_hit_info
{

  int64_t  seqidx;
  int      strand; /* 0 = forward strand, 1 = reverse strand */
  int      duplicate;

  int hmm_start;
  int hmm_end;
  int64_t seq_start;
  int64_t seq_end;

  float aliscore;

  int node_id;  /* node in splice graph that corresponds to the hit */

} SPLICE_HIT_INFO;

typedef struct _splice_saved_hits
{
  SPLICE_HIT_INFO **srt;
  SPLICE_HIT_INFO *unsrt;

  uint64_t Nalloc;
  uint64_t N;

  int is_sorted;

} SPLICE_SAVED_HITS;



