#ifndef NCSU_PENALTIES_H
#define NCSU_PENALTIES_H

/* written by vbabin-at-ncsu-dot-edu in 05/2008 */

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

struct _ncsu_penalty_pair_t {
    size_t set_id;
    size_t element_id;
};

typedef struct _ncsu_penalty_pair_t ncsu_penalty_pair_t;

/*
*                           1  2  3
*  example: N = 3, sizes = (2, 3, 1) and goal = 2
*
*  (i.e., there are 2 elements with penalty = 1; 3 elements
*   with penalty equal 2 and 1 element of penalty 3)
*
*  result: (s = set_id, e = element_id)
*    1 : (s = 0, e = 0), (s = 0, e = 1)
*    2 : (s = 1, e = 0)
*    3 : (s = 1, e = 1)
*    4 : (s = 1, e = 2)
*
*/

typedef void (* ncsu_penalty_func_t)
    (size_t npairs, const ncsu_penalty_pair_t* /* npairs <= N */, void*);

void ncsu_penalties(size_t N /* > 0 */, const size_t *sizes /* !NULL */,
                    size_t goal /* <= N */,
                    ncsu_penalty_func_t /* !NULL */, void*);

size_t ncsu_penalties_size(size_t N, const size_t *sizes, size_t goal);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* NCSU_PENALTIES_H */
