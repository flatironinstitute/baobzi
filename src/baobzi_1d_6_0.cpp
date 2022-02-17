#include "baobzi_template.hpp"
#include "baobzi.h"
#include "baobzi/macros.h"

namespace baobzi {
template
typename Function<1, 6, 0>::CoeffVec Function<1, 6, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 6, 0>::VanderMat> Function<1, 6, 0>::VLU_;
}

extern "C" {
BAOBZI_DEFS(1, 6, 0)
}
