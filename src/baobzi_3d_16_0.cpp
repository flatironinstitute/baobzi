#include "baobzi_template.hpp"
#include "baobzi.h"
#include "baobzi/macros.h"

namespace baobzi {
template
typename Function<3, 16, 0>::VecOrderD Function<3, 16, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 16, 0>::VanderMat> Function<3, 16, 0>::VLU_;
}

extern "C" {
BAOBZI_DEFS(3, 16, 0)
}
