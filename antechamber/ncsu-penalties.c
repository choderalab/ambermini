/* written by vbabin-at-ncsu-dot-edu in 05/2008 */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include "ncsu-penalties.h"

/******************************************************************************/

static void *xmalloc(size_t n)
{
    void *ptr;

    ptr = malloc(n);
    if (ptr == NULL) {
        fputs("** Error ** : out of memory\n", stdout);
        exit(EXIT_FAILURE);
    }

    return ptr;
}

/******************************************************************************/

static size_t binomial_coefficient(size_t n, size_t k)
{
    size_t i, t, b;

    assert(n > 0);
    assert(k > 0);
    assert(n >= k);

    for (i = k, t = 1; i != n; ++i)
        t *= i + 1;

    for (i = 0, b = 1; i + k != n; ++i)
        b *= i + 1;

    return (t/b);
}

/******************************************************************************/

/*
*  Iterates through all k-combinations of (0, ..., n - 1), i.e.,
*  for n = 4 and k = 2 the output is (0, 1), (0, 2), (0, 3),
*  (1, 2), (1, 3), (2, 3)
*/

typedef struct {
    size_t  n, k; /* elements go from 0 to n - 1 ; 1 <= k <= n */
    size_t  nc, *c; /* size == nc >= k */
} combination_t;

static void combination_init(combination_t *self)
{
    assert(self != NULL);

    self->n = 0;
    self->k = 0;

    self->nc = 0;
    self->c = NULL;
}

static void combination_fini(combination_t *self)
{
    assert(self != NULL);

    if (self->c != NULL) {
        assert(self->nc > 0);
        free(self->c);
    }
}

static void combination_setup(combination_t *self, size_t n, size_t k)
{
    assert(n != 0);
    assert(k != 0);

    assert(k <= n);
    assert(self != NULL);

    self->n = n;
    self->k = k;

    if (k > self->nc) {
        if (self->c != NULL)
            free(self->c);
        self->c = (size_t *) xmalloc(k*sizeof(size_t));
        self->nc = k;
    }

    self->c[0] = self->n; /* impossible value */
}

static size_t* combination_next(combination_t *self)
{
    size_t r, n, i, mark;

    assert(self != NULL);

    assert(self->n > 0);
    assert(self->k > 0);

    assert(self->k <= self->n);
    assert(self->c != NULL);
    assert(self->nc >= self->k);

    if (self->c[0] == self->n) {
        for (i = 0; i < self->k; ++i)
            self->c[i] = i;

        return self->c;
    }

    r = self->k - 1;
    n = self->n - 1;

    for (mark = self->k;; --r, --n) {
        if (self->c[r] == n) {
            if (r != 0) {
                mark = --r;
                ++r;
                continue;
            } else {
                self->c[0] = self->n;
                return NULL;
            }
        } else {
            if (mark < self->k) {
                size_t nm = self->c[mark] + 1;
                for (i = mark; i != self->k; ++i, ++nm)
                    self->c[i] = nm;

                return self->c;
            }
            for (i = 0; i < self->n; ++i) {
                if (self->c[r] == i) {
                    self->c[r] = ++i;
                    return self->c;
                }
            }
        }
    }

    return NULL; /* should not be reached */
}

#define combination_is_null(self) ((self)->c[0] == (self)->n)

#ifdef ENABLE_COMBINATION_MAIN
int main(int argc, char** argv)
{
    const size_t *c;
    combination_t cmb;

    combination_init(&cmb);

    combination_setup(&cmb, 7, 3);
    for (c = combination_next(&cmb); c != NULL; c = combination_next(&cmb)) {
        printf("%lu %lu %lu\n",
            (unsigned long) c[0],
            (unsigned long) c[1],
            (unsigned long) c[2]);
    }

    combination_fini(&cmb);
    exit(EXIT_SUCCESS);
}
#endif /* ENABLE_COMBINATION_MAIN */

/******************************************************************************/

/*
* Iterates through non-negative integer solutions (c[0], ..., c[K - 1]) of
* c[0]*1 + ... + c[K - 1]*K = N, (both K and N are positive integers). The
* algorithm is from:
*
* Frank Stockmal
* "Algorithm 95: Generation of partitions in part-count form",
* Communications of the ACM, Volume 5, Issue 6 (June 1962),
* Page 344, 1962
*/

typedef struct {
    size_t K, N;
    size_t nc, *c; /* size == nc >= K */
} int_part_t;

static void int_part_init(int_part_t *self)
{
    assert(self != NULL);

    self->K = 0;
    self->N = 0;

    self->nc = 0;
    self->c = NULL;
}

static void int_part_fini(int_part_t *self)
{
    assert(self != NULL);

    if (self->c != NULL) {
        assert(self->nc > 0);
        free(self->c);
    }
}

static void int_part_setup(int_part_t *self, size_t K, size_t N)
{
    assert(K != 0);
    assert(N != 0);
    assert(self != NULL);

    self->N = N;
    self->K = K;

    if (K > self->nc) {
        if (self->c != NULL)
            free(self->c);
        self->c = (size_t *) xmalloc(K*sizeof(size_t));
        self->nc = K;
    }

    self->c[0] = N + 1; /* impossible value */
}

static const size_t* int_part_next(int_part_t *self)
{
    size_t i, j, a;

    assert(self != NULL);
    assert(self->K != 0);
    assert(self->N != 0);
    assert(self->c != NULL);
    assert(self->nc >= self->K);

    /* K == 1 case first */
    if (self->K == 1 && self->c[0] == self->N) {
        self->c[0] = self->N + 1;
        return NULL;
    }

    if (self->c[0] == self->N + 1)
        goto L3;

    j = 1;
    a = self->c[0];

L0: if (a <= j)
        goto L2;

    ++(self->c[j]);
    self->c[0] = a - j - 1;

L1: for (i = 1; i < j; ++i)
        self->c[i] = 0;
    goto L4;

L2: if (j + 1 == self->K) {
        self->c[0] = self->N + 1; /* impossible value */
        return NULL;
    }

    a += (j + 1)*self->c[j];
    ++j; goto L0;

L3: /* the first partition : N*1 + 0*2 + ... + 0*K */
    self->c[0] = self->N;
    j = self->K;
    goto L1;

L4:
    return self->c;
}

/******************************************************************************/

#ifdef ENABLE_INT_PART_MAIN
int main(int argc, char** argv)
{
    size_t N = 0, K = 0;

    const size_t *t;
    int_part_t iter;

    int_part_init(&iter); /* constructor */

    if (argc > 1) {
        char* endptr;
        N = (size_t) strtoul(argv[1], &endptr, 10);
        if (endptr == argv[1])
            N = 0;
    }

    if (argc > 2) {
        char* endptr;
        K = (size_t) strtoul(argv[2], &endptr, 10);
        if (endptr == argv[2])
            K = 0;
    }

    if (N == 0)
        N = 7;

    if (K == 0)
        K = 3;

    int_part_setup(&iter, K, N);
    for (t = int_part_next(&iter); t != NULL; t = int_part_next(&iter)) {
        printf("%lu = ", (unsigned long) iter.N);
        for (K = 0; K + 1 < iter.K; ++K)
            printf("%lu*%lu + ", (unsigned long) t[K], (unsigned long) K + 1);
        printf("%lu*%lu\n", (unsigned long) t[K], (unsigned long) K + 1);
    }

    int_part_fini(&iter); /* destructor */

    exit(EXIT_SUCCESS);
}
#endif /* ENABLE_INT_PART_MAIN */

/******************************************************************************/

size_t ncsu_penalties_size(size_t N, const size_t *sizes, size_t g)
{
    int_part_t prt;

    size_t i, result;
    const size_t *tpl;

    assert(g > 0);
    assert(g <= N);
    assert(sizes != NULL);

    result = 0;

    int_part_init(&prt);
    int_part_setup(&prt, N, g);

    /* go through partitions */
    for (tpl = int_part_next(&prt); tpl != NULL; tpl = int_part_next(&prt)) {

        /* check if 'tpl' is applicable */
        for (i = 0; i < g; ++i) {
            if (tpl[i] > sizes[i]) {
                i = g + 1; /* impossible value */
                break;
            }
        }

        /* too few elements in some set, skip this partition */
        if (i == g + 1)
            continue;

        /* count the size for this partition (keep it in 'N') */
        for (N = 1, i = 0; i < g; ++i) {
            if (tpl[i] != 0)
                N *= binomial_coefficient(sizes[i], tpl[i]);
        }

        /* accumulate it */
        result += N;
    }

    int_part_fini(&prt);

    return result;
}

/******************************************************************************/

void ncsu_penalties(size_t N, const size_t *sizes, size_t g,
                    ncsu_penalty_func_t cb, void *cb_data)
{
    size_t i;
    int_part_t prt;
    combination_t *cmb;
    ncsu_penalty_pair_t *pairs;

    const size_t *tpl;

    assert(g > 0);
    assert(g <= N);
    assert(sizes != NULL);

    assert(cb != NULL);

    pairs = (ncsu_penalty_pair_t *) xmalloc(N*sizeof(ncsu_penalty_pair_t));

    cmb = (combination_t *) xmalloc(N*sizeof(combination_t));
    for (i = 0; i < N; ++i)
        combination_init(cmb + i);

    int_part_init(&prt);
    int_part_setup(&prt, N, g);

    /* go through partitions */
    for (tpl = int_part_next(&prt); tpl != NULL; tpl = int_part_next(&prt)) {
        size_t last = 0;

        /* check if 'tpl' is applicable */
        for (i = 0; i < g; ++i) {
            if (tpl[i] > sizes[i]) {
                i = g + 1; /* impossible value */
                break;
            }
        }

        /* too few elements in some set, skip this partition */
        if (i == g + 1)
            continue;

        /* setup the combinations */
        for (i = 0; i < g; ++i) {
            if (tpl[i] != 0) {
                assert(tpl[i] <= N);
                combination_setup(cmb + i, sizes[i], tpl[i]);
                (void) combination_next(cmb + i);
                last = i;
            }
        }
        /* go through Cartesian product of the applicable combinations */
        while (1) {
            size_t npairs = 0;

            /* report current element */
            for (i = 0; i < g; ++i) {            
                if (tpl[i] != 0) { /* use these combination */
                    size_t l;
                    for (l = 0; l < cmb[i].k; ++l) {
                        assert(npairs < N);
                        pairs[npairs].set_id = i;
                        pairs[npairs].element_id = cmb[i].c[l];
                        ++npairs;
                    }
                } /* tpl[i] != 0 */
            }

            assert(npairs > 0);
            (*cb)(npairs, pairs, cb_data);

            /* advance to the next one */
            for (i = 0; i < g; ++i) {
                if (tpl[i] != 0) {
                    if (combination_next(cmb + i) != NULL) {
                        break;
                    } else {
                        if (i == last)
                            break;
                        (void) combination_next(cmb + i);
                    }
                } /* tpl[i] != 0 */
            }

            /* advance to a next partition */
            if (i == last && combination_is_null(cmb + i))
                break;
        } /* while (1) */
    }

    /* call destructors, free allocated memory */

    int_part_fini(&prt);

    for (i = 0; i < N; ++i)
        combination_fini(cmb + i);

    free(cmb);
    free(pairs);
}

/******************************************************************************/
