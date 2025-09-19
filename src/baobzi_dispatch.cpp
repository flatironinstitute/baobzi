#include <baobzi/baobzi_binding_impl.hpp>

#ifndef BAOBZI_ISA
#error "BAOBZI_ISA not defined"
#endif

constexpr baobzi::baobzi_isa_t ISA = baobzi::baobzi_isa_t(BAOBZI_ISA);

template baobzi_t baobzi::baobzi_init<ISA>(const baobzi_input_t *input, const double *center,
                                           const double *half_width_in);

template void baobzi::baobzi_eval_multi<ISA>(const baobzi_t f, const double *input_point, double *output_point,
                                             int ntrg);

template void baobzi::baobzi_stats<ISA>(const baobzi_t f);

template baobzi_t baobzi::baobzi_free<ISA>(const baobzi_t f);
