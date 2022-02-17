#include "baobzi_template.hpp"
#include "baobzi.h"
#include "baobzi/macros.h"

namespace baobzi {
template
typename Function<2, 16, 2>::CoeffVec Function<2, 16, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 16, 2>::VanderMat> Function<2, 16, 2>::VLU_;
}

extern "C" {
BAOBZI_DEFS(2, 16, 2)
}
