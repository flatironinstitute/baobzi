#ifndef BAOBZI_BINDING_HPP
#define BAOBZI_BINDING_HPP

#include <baobzi.hpp>

namespace baobzi {
typedef enum {
    GENERIC = 0,
    AVX = 1,
    AVX2 = 2,
    AVX512 = 3,
} baobzi_isa_t;

template <baobzi_isa_t ISA>
baobzi_t baobzi_init(const baobzi_input_t *input, const double *center, const double *half_width_in);

template <baobzi_isa_t ISA>
void baobzi_eval_multi(const baobzi_t f, const double *input_point, double *output_point, int ntrg);

template <baobzi_isa_t ISA>
void baobzi_stats(const baobzi_t f);

template <baobzi_isa_t ISA>
baobzi_t baobzi_free(const baobzi_t f);

} // namespace baobzi
#endif
