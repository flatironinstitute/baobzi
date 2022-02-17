#include "baobzi_template.hpp"
#include "baobzi.h"
#include "baobzi/macros.h"

namespace baobzi {
template
typename Function<1, 16, 3>::CoeffVec Function<1, 16, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 16, 3>::VanderMat> Function<1, 16, 3>::VLU_;
}

extern "C" {
BAOBZI_DEFS(1, 16, 3)
}
