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
    void (*save)(const void *, const char *);
    void (*free)(void *);
} baobzi_struct;

typedef baobzi_struct* baobzi_t;

double baobzi_eval(const baobzi_t func, const double *x);
void baobzi_save(const baobzi_t func, const char *filename);
baobzi_t baobzi_restore(double (*fin)(const double *), const char *filename);
baobzi_t baobzi_free(baobzi_t func);
baobzi_t baobzi_init(double (*)(const double *), uint16_t dim, uint16_t order, const double *center,
                     const double *half_length, const double tol);

#include "baobzi/baobzi_decls.h"

#ifdef __cplusplus
}
#endif

#endif
