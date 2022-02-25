#include "baobzi_template.hpp"
#include "baobzi.h"
#include "baobzi/macros.h"

namespace baobzi {
template
typename Function<3, 10, 3>::VecOrderD Function<3, 10, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 10, 3>::VanderMat> Function<3, 10, 3>::VLU_;
}

extern "C" {
BAOBZI_DEFS(3, 10, 3)
}
