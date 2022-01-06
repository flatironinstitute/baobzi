#ifndef BAOBZI_H
#define BAOBZI_H

#include "baobzi/macros.h"

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    void *obj;
    int DIM;
    int ORDER;
    double (*f_)(const double *);
    double (*eval)(const void *, const double *);
    void (*free)(void *);
} baobzi_t;

double baobzi_eval(const baobzi_t *func, const double *x);
void baobzi_free(baobzi_t *func);
baobzi_t baobzi_init(double (*)(const double *), uint16_t dim, uint16_t order, const double *center,
                     const double *half_length, const double tol);

BAOBZI_DECLS(1, 6, 0)
BAOBZI_DECLS(1, 8, 0)
BAOBZI_DECLS(1, 10, 0)
BAOBZI_DECLS(1, 12, 0)
BAOBZI_DECLS(1, 14, 0)
BAOBZI_DECLS(1, 16, 0)
BAOBZI_DECLS(2, 6, 0)
BAOBZI_DECLS(2, 8, 0)
BAOBZI_DECLS(2, 10, 0)
BAOBZI_DECLS(2, 12, 0)
BAOBZI_DECLS(2, 14, 0)
BAOBZI_DECLS(2, 16, 0)
BAOBZI_DECLS(3, 6, 0)
BAOBZI_DECLS(3, 8, 0)
BAOBZI_DECLS(3, 10, 0)
BAOBZI_DECLS(3, 12, 0)
BAOBZI_DECLS(3, 14, 0)
BAOBZI_DECLS(3, 16, 0)

BAOBZI_DECLS(1, 6, 1)
BAOBZI_DECLS(1, 8, 1)
BAOBZI_DECLS(1, 10, 1)
BAOBZI_DECLS(1, 12, 1)
BAOBZI_DECLS(1, 14, 1)
BAOBZI_DECLS(1, 16, 1)
BAOBZI_DECLS(2, 6, 1)
BAOBZI_DECLS(2, 8, 1)
BAOBZI_DECLS(2, 10, 1)
BAOBZI_DECLS(2, 12, 1)
BAOBZI_DECLS(2, 14, 1)
BAOBZI_DECLS(2, 16, 1)
BAOBZI_DECLS(3, 6, 1)
BAOBZI_DECLS(3, 8, 1)
BAOBZI_DECLS(3, 10, 1)
BAOBZI_DECLS(3, 12, 1)
BAOBZI_DECLS(3, 14, 1)
BAOBZI_DECLS(3, 16, 1)

BAOBZI_DECLS(1, 6, 2)
BAOBZI_DECLS(1, 8, 2)
BAOBZI_DECLS(1, 10, 2)
BAOBZI_DECLS(1, 12, 2)
BAOBZI_DECLS(1, 14, 2)
BAOBZI_DECLS(1, 16, 2)
BAOBZI_DECLS(2, 6, 2)
BAOBZI_DECLS(2, 8, 2)
BAOBZI_DECLS(2, 10, 2)
BAOBZI_DECLS(2, 12, 2)
BAOBZI_DECLS(2, 14, 2)
BAOBZI_DECLS(2, 16, 2)
BAOBZI_DECLS(3, 6, 2)
BAOBZI_DECLS(3, 8, 2)
BAOBZI_DECLS(3, 10, 2)
BAOBZI_DECLS(3, 12, 2)
BAOBZI_DECLS(3, 14, 2)
BAOBZI_DECLS(3, 16, 2)

BAOBZI_DECLS(1, 6, 3)
BAOBZI_DECLS(1, 8, 3)
BAOBZI_DECLS(1, 10, 3)
BAOBZI_DECLS(1, 12, 3)
BAOBZI_DECLS(1, 14, 3)
BAOBZI_DECLS(1, 16, 3)
BAOBZI_DECLS(2, 6, 3)
BAOBZI_DECLS(2, 8, 3)
BAOBZI_DECLS(2, 10, 3)
BAOBZI_DECLS(2, 12, 3)
BAOBZI_DECLS(2, 14, 3)
BAOBZI_DECLS(2, 16, 3)
BAOBZI_DECLS(3, 6, 3)
BAOBZI_DECLS(3, 8, 3)
BAOBZI_DECLS(3, 10, 3)
BAOBZI_DECLS(3, 12, 3)
BAOBZI_DECLS(3, 14, 3)
BAOBZI_DECLS(3, 16, 3)


#ifdef __cplusplus
}
#endif

#endif
