#include "baobzi.h"
#include "baobzi/macros.h"

#include <cstdint>
#include <iostream>

enum ISET { GENERIC, AVX, AVX2, AVX512 };

extern "C" {
double baobzi_eval(const baobzi_t *func, const double *x) { return func->eval(func->obj, x); };

void baobzi_free(baobzi_t *func) {
    func->free(func->obj);
    func->obj = nullptr;
    func->eval = nullptr;
    func->free = nullptr;
};

baobzi_t baobzi_init(double (*fin)(const double *), uint16_t dim, uint16_t order, const double *center,
                     const double *half_length, const double tol) {
    baobzi_t res;
    res.f_ = fin;
    res.DIM = dim;
    res.ORDER = order;

    int iset = ISET::GENERIC;
    if (__builtin_cpu_supports("avx2"))
        iset = ISET::AVX;
    if (__builtin_cpu_supports("avx2"))
        iset = ISET::AVX2;
    if (__builtin_cpu_supports("avx512f"))
        iset = ISET::AVX512;

    switch (BAOBZI_JOIN(dim, order, iset)) {
        BAOBZI_CASE(1, 6, 0)
        BAOBZI_CASE(1, 8, 0)
        BAOBZI_CASE(1, 10, 0)
        BAOBZI_CASE(1, 12, 0)
        BAOBZI_CASE(1, 14, 0)
        BAOBZI_CASE(1, 16, 0)
        BAOBZI_CASE(2, 6, 0)
        BAOBZI_CASE(2, 8, 0)
        BAOBZI_CASE(2, 10, 0)
        BAOBZI_CASE(2, 12, 0)
        BAOBZI_CASE(2, 14, 0)
        BAOBZI_CASE(2, 16, 0)
        BAOBZI_CASE(3, 6, 0)
        BAOBZI_CASE(3, 8, 0)
        BAOBZI_CASE(3, 10, 0)
        BAOBZI_CASE(3, 12, 0)
        BAOBZI_CASE(3, 14, 0)
        BAOBZI_CASE(3, 16, 0)
        BAOBZI_CASE(1, 6, 1)
        BAOBZI_CASE(1, 8, 1)
        BAOBZI_CASE(1, 10, 1)
        BAOBZI_CASE(1, 12, 1)
        BAOBZI_CASE(1, 14, 1)
        BAOBZI_CASE(1, 16, 1)
        BAOBZI_CASE(2, 6, 1)
        BAOBZI_CASE(2, 8, 1)
        BAOBZI_CASE(2, 10, 1)
        BAOBZI_CASE(2, 12, 1)
        BAOBZI_CASE(2, 14, 1)
        BAOBZI_CASE(2, 16, 1)
        BAOBZI_CASE(3, 6, 1)
        BAOBZI_CASE(3, 8, 1)
        BAOBZI_CASE(3, 10, 1)
        BAOBZI_CASE(3, 12, 1)
        BAOBZI_CASE(3, 14, 1)
        BAOBZI_CASE(3, 16, 1)
        BAOBZI_CASE(1, 6, 2)
        BAOBZI_CASE(1, 8, 2)
        BAOBZI_CASE(1, 10, 2)
        BAOBZI_CASE(1, 12, 2)
        BAOBZI_CASE(1, 14, 2)
        BAOBZI_CASE(1, 16, 2)
        BAOBZI_CASE(2, 6, 2)
        BAOBZI_CASE(2, 8, 2)
        BAOBZI_CASE(2, 10, 2)
        BAOBZI_CASE(2, 12, 2)
        BAOBZI_CASE(2, 14, 2)
        BAOBZI_CASE(2, 16, 2)
        BAOBZI_CASE(3, 6, 2)
        BAOBZI_CASE(3, 8, 2)
        BAOBZI_CASE(3, 10, 2)
        BAOBZI_CASE(3, 12, 2)
        BAOBZI_CASE(3, 14, 2)
        BAOBZI_CASE(3, 16, 2)
        BAOBZI_CASE(1, 6, 3)
        BAOBZI_CASE(1, 8, 3)
        BAOBZI_CASE(1, 10, 3)
        BAOBZI_CASE(1, 12, 3)
        BAOBZI_CASE(1, 14, 3)
        BAOBZI_CASE(1, 16, 3)
        BAOBZI_CASE(2, 6, 3)
        BAOBZI_CASE(2, 8, 3)
        BAOBZI_CASE(2, 10, 3)
        BAOBZI_CASE(2, 12, 3)
        BAOBZI_CASE(2, 14, 3)
        BAOBZI_CASE(2, 16, 3)
        BAOBZI_CASE(3, 6, 3)
        BAOBZI_CASE(3, 8, 3)
        BAOBZI_CASE(3, 10, 3)
        BAOBZI_CASE(3, 12, 3)
        BAOBZI_CASE(3, 14, 3)
        BAOBZI_CASE(3, 16, 3)
    default: {
        std::cerr << "BAOBZI ERROR: Unable to initialize Baobzi function with variables (DIM, ORDER): (" << dim << ", "
                  << order << ")\n";
        break;
    }
    }

    return res;
}
}
