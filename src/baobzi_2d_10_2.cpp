#include "baobzi_template.hpp"
#include "baobzi.h"
#include "baobzi/macros.h"

namespace baobzi {
template
typename Function<2, 10, 2>::VecOrderD Function<2, 10, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 10, 2>::VanderMat> Function<2, 10, 2>::VLU_;
}

extern "C" {
BAOBZI_DEFS(2, 10, 2)
}
