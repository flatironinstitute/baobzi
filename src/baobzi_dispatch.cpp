#include <baobzi/baobzi_binding_impl.hpp>

#ifndef BAOBZI_ISA
#error "BAOBZI_ISA not defined"
#endif

constexpr baobzi::baobzi_isa_t ISA = baobzi::baobzi_isa_t(BAOBZI_ISA);

extern "C" {
baobzi_t baobzi_init(const baobzi_input_t *input, const double *center, const double *half_width_in) {
    return baobzi::baobzi_init<ISA>(input, center, half_width_in);
}

void baobzi_eval_multi(const baobzi_t f, const double *x, double *res, int ntrg) {
    return baobzi::baobzi_eval_multi<ISA>(f, x, res, ntrg);
}

void baobzi_stats(const baobzi_t f) { return baobzi::baobzi_stats<ISA>(f); }

baobzi_t baobzi_free(baobzi_t f) { return baobzi::baobzi_free<ISA>(f); }
}
