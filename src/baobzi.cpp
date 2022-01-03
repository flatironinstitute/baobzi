#include "baobzi.h"
#include "baobzi.hpp"
#include "baobzi/macros.h"

#include <cstdint>
#include <iostream>
#include <new>

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

    switch (BAOBZI_JOIN(dim, order)) {
        BAOBZI_CASE(1, 6)
        BAOBZI_CASE(1, 8)
        BAOBZI_CASE(1, 10)
        BAOBZI_CASE(1, 12)
        BAOBZI_CASE(1, 14)
        BAOBZI_CASE(1, 16)
        BAOBZI_CASE(2, 6)
        BAOBZI_CASE(2, 8)
        BAOBZI_CASE(2, 10)
        BAOBZI_CASE(2, 12)
        BAOBZI_CASE(2, 14)
        BAOBZI_CASE(2, 16)
    default: {
        std::cerr << "BAOBZI ERROR: Unable to initialize Baobzi function with variables (DIM, ORDER): (" << dim << ", "
                  << order << ")\n";
        break;
    }
    }

    return res;
}
}
