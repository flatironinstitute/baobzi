#include "baobzi_template.hpp"
#include "baobzi.h"
#include "baobzi/macros.h"

namespace baobzi {
template
typename Function<2, 14, 1>::CoeffVec Function<2, 14, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 14, 1>::VanderMat> Function<2, 14, 1>::VLU_;
}

extern "C" {
BAOBZI_DEFS(2, 14, 1)
}
