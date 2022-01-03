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

BAOBZI_DECLS(1, 6)
BAOBZI_DECLS(1, 8)
BAOBZI_DECLS(1, 10)
BAOBZI_DECLS(1, 12)
BAOBZI_DECLS(1, 14)
BAOBZI_DECLS(1, 16)
BAOBZI_DECLS(2, 6)
BAOBZI_DECLS(2, 8)
BAOBZI_DECLS(2, 10)
BAOBZI_DECLS(2, 12)
BAOBZI_DECLS(2, 14)
BAOBZI_DECLS(2, 16)

#ifdef __cplusplus
}
#endif

#endif
