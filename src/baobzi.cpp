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
    if (__builtin_cpu_supports("avx"))
        iset = ISET::AVX;
    if (__builtin_cpu_supports("avx2"))
        iset = ISET::AVX2;
    if (__builtin_cpu_supports("avx512f"))
        iset = ISET::AVX512;

    switch (BAOBZI_JOIN(dim, order, iset)) {
        #include "baobzi/baobzi_cases.h"
    default: {
        std::cerr << "BAOBZI ERROR: Unable to initialize Baobzi function with variables (DIM, ORDER): (" << dim << ", "
                  << order << ")\n";
        break;
    }
    }

    return res;
}
}
