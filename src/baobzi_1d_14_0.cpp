#include "baobzi_template.hpp"
#include "baobzi.h"
#include "baobzi/macros.h"

namespace baobzi {
template
typename Function<1, 14, 0>::VecOrderD Function<1, 14, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 14, 0>::VanderMat> Function<1, 14, 0>::VLU_;
}

extern "C" {
BAOBZI_DEFS(1, 14, 0)
}
