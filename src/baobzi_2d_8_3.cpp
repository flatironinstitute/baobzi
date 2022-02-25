#include "baobzi_template.hpp"
#include "baobzi.h"
#include "baobzi/macros.h"

namespace baobzi {
template
typename Function<2, 8, 3>::VecOrderD Function<2, 8, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 8, 3>::VanderMat> Function<2, 8, 3>::VLU_;
}

extern "C" {
BAOBZI_DEFS(2, 8, 3)
}
