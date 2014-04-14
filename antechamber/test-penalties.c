#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include "ncsu-penalties.h"

#define NR 3
#define NC 5

static const size_t maxaps = 7;

/* penalties matrix */
static const size_t M[NR][NC] = {
    {0, 1, 2, 0, 8},
    {1, 2, 3, 4, 5},
    {1, 1, 0, 2, 3}
};

typedef struct {
    size_t addr;
    size_t valence;
} pair_t;

#define TO_LU(x) ((unsigned long) (x))

static void* do_malloc(size_t n)
{
    void *ptr = malloc(n);
    if (ptr == NULL) {
        fprintf(stdout,
            "** Error ** : failed to malloc() %lu bytes\n", TO_LU(n));
        exit(EXIT_FAILURE);
    }
    return ptr;
}

static void penalties_cb(size_t np, const ncsu_penalty_pair_t* p, void* data)
{
    size_t i;

    pair_t **pairs = (pair_t **) data;

    for (i = 0; i < np; ++i) {
        size_t penalty = p[i].set_id + 1;
        const pair_t *pair = &(pairs[penalty][p[i].element_id]);
        printf(" (%lu|%lu, %lu)", TO_LU(penalty),
               TO_LU(pair->addr), TO_LU(pair->valence));
    }

    putchar('\n');
}

int main(int argc, char** argv)
{
    size_t i, j;

    static size_t *idx = 0;
    static size_t *npairs = 0;
    static pair_t **pairs = NULL;

    puts("\n penalties matrix :\n");

    for (i = 0; i < NR; ++i) {
        for (j = 0; j < NC; ++j) {
            printf(" %2lu", TO_LU(M[i][j]));
        }
        putchar('\n');
    }

    printf("\n maxaps = %lu\n", TO_LU(maxaps));

    npairs = (size_t *) do_malloc((maxaps + 1)*sizeof(size_t));
    for (i = 0; i <= maxaps; ++i) {
        npairs[i] = 0;
    }

    for (i = 0; i < NR; ++i) {
        for (j = 0; j < NC; ++j) {
            if (M[i][j] <= maxaps)
                ++npairs[M[i][j]];
        }
    }

    pairs = (pair_t **) do_malloc((maxaps + 1)*sizeof(pair_t*));
    for (i = 0; i <= maxaps; ++i) {
        pairs[i] = (pair_t *) do_malloc(npairs[i]*sizeof(pair_t));
    }

    idx = (size_t *) do_malloc((maxaps + 1)*sizeof(size_t));
    for (i = 0; i <= maxaps; ++i) {
        idx[i] = 0;
    }

    for (i = 0; i < NR; ++i) {
        for (j = 0; j < NC; ++j) {
            if (M[i][j] <= maxaps) {
                assert(idx[M[i][j]] < npairs[M[i][j]]);
                pairs[M[i][j]][idx[M[i][j]]].addr = i;
                pairs[M[i][j]][idx[M[i][j]]].valence = j;
                ++idx[M[i][j]];
            }
        }
    }

    /* done with initialization */

    for (i = 0; i <= maxaps; ++i) {
        printf("\nm = %lu :", TO_LU(i));
        for (j = 0; j < npairs[i]; ++j) {
            printf(" (%lu, %lu)",
                TO_LU(pairs[i][j].addr), TO_LU(pairs[i][j].valence));
        }
    }

    putchar('\n');
    putchar('\n');

    for (i = 1; i <= maxaps; ++i) {
        printf(">>  calling  ncsu_penalties() for goal = %lu >>\n", TO_LU(i));
        printf("    ||  expected size is %lu ||\n",
            TO_LU(ncsu_penalties_size(maxaps, npairs + 1, i)));
        ncsu_penalties(maxaps, npairs + 1, i, penalties_cb, (void *) pairs);
        printf("<< back from ncsu_penalties() for goal = %lu <<\n\n", TO_LU(i));
    }

    return 0;
}
