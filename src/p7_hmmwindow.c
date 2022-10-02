/* The Plan7 HMMWINDOW data structure, which holds a compact representation
 * of substitution scores and maximal extensions, used by nhmmer.
 *
 * Contents:
 *   1. The P7_HMMWINDOW object: allocation, initialization, destruction.
 */
#include "p7_config.h"
#include "hmmer.h"

/*********************************************************************
 *# 1. The P7_MSVDATA object: allocation, initialization, destruction.
 *********************************************************************/


/* Function:  p7_hmmwindow_init()
 *
 * Synopsis:  initialize the object used to store a list of sequence windows
 *
 * Returns:   eslEMEM in event of allocation failure, otherwise eslOK
 */
int
p7_hmmwindow_init (P7_HMM_WINDOWLIST *list) {
  int status;
  list->size = 10000;
  list->count = 0;
  ESL_ALLOC(list->windows, list->size * sizeof(P7_HMM_WINDOW));

  return eslOK;

ERROR:
  return eslEMEM;

}

/* Function:  p7_hmmwindow_new()
 *
 * Synopsis:  Return a pointer to the next window element on the list
 *
 * Purpose:   Accepts <id>, <pos>, <fm_pos>, <k>, <length>, <score>,
 *            and <complementarity>, assigns those to the next window
 *            element, then returns it, increasing the size of the
 *            list, if necessary.
 *
 * Returns:   NULL in event of allocation failure, otherwise pointer to
 *            the next seed diagonal
 */
P7_HMM_WINDOW *
p7_hmmwindow_new (P7_HMM_WINDOWLIST *list, uint32_t id, uint32_t pos, uint32_t fm_pos, uint16_t k, uint32_t length, float score, uint8_t complementarity, uint32_t target_len) {
  int status;
  P7_HMM_WINDOW *window;

  if (list->count == list->size) {
    list->size *= 4;
    ESL_REALLOC(list->windows, list->size * sizeof(P7_HMM_WINDOW));
  }
  window = list->windows + list->count;

  window->id               = id;
  window->n                = pos;
  window->fm_n             = fm_pos;
  window->k                = k;
  window->length           = length;
  window->score            = score;
  window->complementarity  = complementarity;
  window->target_len       = target_len;
  list->count++;

  return window;

ERROR:
  return NULL;
}

/* window_sorter(): qsort's pawn, below */
static int
window_sorter(const void *vw1, const void *vw2)
{
  P7_HMM_WINDOW w1 = *((P7_HMM_WINDOW *) vw1);  /* don't ask. don't change. Don't Panic. */
  P7_HMM_WINDOW w2 = *((P7_HMM_WINDOW *) vw2);

  if      (w1.n > w2.n) return  1;
  else if (w1.n < w2.n) return -1;
  else                    return  0;
}



/* Function:  p7_hmmwindow_SortByStart()
 * Synopsis:  Sorts a hit list.
 *
 * Purpose:   Sorts a top hit list. After this call,
 *            <h->hit[i]> points to the i'th ranked
 *            <P7_HIT> for all <h->N> hits.
 *
 * Returns:   <eslOK> on success.
 */
int
p7_hmmwindow_SortByStart(P7_HMM_WINDOWLIST *w)
{
  if (w->count > 1)  qsort(w->windows, w->count, sizeof(P7_HMM_WINDOW), window_sorter);
  return eslOK;
}


