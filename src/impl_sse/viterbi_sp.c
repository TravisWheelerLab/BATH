/* SSE-accelerated Spliced Viterbi algorithms; full matrix.
 *
 * Contents:
 *   1. p7_Viterbi_SplicedGlobal()
 *   2. Benchmark driver.
 *   3. Unit tests.
 *   4. Test driver.
 */
#include "p7_config.h"

#include <stdio.h>
#include <math.h>

#include <xmmintrin.h>		/* SSE  */
#include <emmintrin.h>		/* SSE2 */

#include "easel.h"
#include "esl_sse.h"
#include "esl_vectorops.h"

#include "hmmer.h"
#include "p7_splice.h"
#include "impl_sse.h"

/* Probability-space constant for the P->M transition (= exp(log(4.58e-5))). */
#define TSC_P_PROB  4.58e-5f

/*****************************************************************
 * 1. p7_Viterbi_SplicedGlobal()
 *****************************************************************/

/*------------------ end, p7_Viterbi_SplicedGlobal() ---------------*/



/*****************************************************************
 * 3. Benchmark driver.
 *****************************************************************/



/*----------------- end, benchmark driver -----------------------*/



/*****************************************************************
 * 4. Unit tests.
 *****************************************************************/


/*----------------- end, unit tests -----------------------------*/



/*****************************************************************
 * 5. Test driver.
 *****************************************************************/


/*----------------- end, test driver ----------------------------*/
