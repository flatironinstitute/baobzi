#include "baobzi_template.hpp"
#include "baobzi.h"
#include "baobzi/macros.h"

namespace baobzi {
template
typename Function<1, 6, 0>::VecOrderD Function<1, 6, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 6, 0>::VanderMat> Function<1, 6, 0>::VLU_;

template
typename Function<1, 8, 0>::VecOrderD Function<1, 8, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 8, 0>::VanderMat> Function<1, 8, 0>::VLU_;

template
typename Function<1, 10, 0>::VecOrderD Function<1, 10, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 10, 0>::VanderMat> Function<1, 10, 0>::VLU_;

template
typename Function<1, 12, 0>::VecOrderD Function<1, 12, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 12, 0>::VanderMat> Function<1, 12, 0>::VLU_;

template
typename Function<1, 14, 0>::VecOrderD Function<1, 14, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 14, 0>::VanderMat> Function<1, 14, 0>::VLU_;

template
typename Function<1, 16, 0>::VecOrderD Function<1, 16, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 16, 0>::VanderMat> Function<1, 16, 0>::VLU_;

template
typename Function<2, 6, 0>::VecOrderD Function<2, 6, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 6, 0>::VanderMat> Function<2, 6, 0>::VLU_;

template
typename Function<2, 8, 0>::VecOrderD Function<2, 8, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 8, 0>::VanderMat> Function<2, 8, 0>::VLU_;

template
typename Function<2, 10, 0>::VecOrderD Function<2, 10, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 10, 0>::VanderMat> Function<2, 10, 0>::VLU_;

template
typename Function<2, 12, 0>::VecOrderD Function<2, 12, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 12, 0>::VanderMat> Function<2, 12, 0>::VLU_;

template
typename Function<2, 14, 0>::VecOrderD Function<2, 14, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 14, 0>::VanderMat> Function<2, 14, 0>::VLU_;

template
typename Function<2, 16, 0>::VecOrderD Function<2, 16, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 16, 0>::VanderMat> Function<2, 16, 0>::VLU_;

template
typename Function<3, 6, 0>::VecOrderD Function<3, 6, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 6, 0>::VanderMat> Function<3, 6, 0>::VLU_;

template
typename Function<3, 8, 0>::VecOrderD Function<3, 8, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 8, 0>::VanderMat> Function<3, 8, 0>::VLU_;

template
typename Function<3, 10, 0>::VecOrderD Function<3, 10, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 10, 0>::VanderMat> Function<3, 10, 0>::VLU_;

template
typename Function<3, 12, 0>::VecOrderD Function<3, 12, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 12, 0>::VanderMat> Function<3, 12, 0>::VLU_;

template
typename Function<3, 14, 0>::VecOrderD Function<3, 14, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 14, 0>::VanderMat> Function<3, 14, 0>::VLU_;

template
typename Function<3, 16, 0>::VecOrderD Function<3, 16, 0>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 16, 0>::VanderMat> Function<3, 16, 0>::VLU_;

}

extern "C" {
BAOBZI_DEFS(1, 6, 0)
BAOBZI_DEFS(1, 8, 0)
BAOBZI_DEFS(1, 10, 0)
BAOBZI_DEFS(1, 12, 0)
BAOBZI_DEFS(1, 14, 0)
BAOBZI_DEFS(1, 16, 0)
BAOBZI_DEFS(2, 6, 0)
BAOBZI_DEFS(2, 8, 0)
BAOBZI_DEFS(2, 10, 0)
BAOBZI_DEFS(2, 12, 0)
BAOBZI_DEFS(2, 14, 0)
BAOBZI_DEFS(2, 16, 0)
BAOBZI_DEFS(3, 6, 0)
BAOBZI_DEFS(3, 8, 0)
BAOBZI_DEFS(3, 10, 0)
BAOBZI_DEFS(3, 12, 0)
BAOBZI_DEFS(3, 14, 0)
BAOBZI_DEFS(3, 16, 0)
}
