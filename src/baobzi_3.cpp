#include "baobzi_template.hpp"
#include "baobzi.h"
#include "baobzi/macros.h"

namespace baobzi {
template
typename Function<1, 6, 3>::VecOrderD Function<1, 6, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 6, 3>::VanderMat> Function<1, 6, 3>::VLU_;

template
typename Function<1, 8, 3>::VecOrderD Function<1, 8, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 8, 3>::VanderMat> Function<1, 8, 3>::VLU_;

template
typename Function<1, 10, 3>::VecOrderD Function<1, 10, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 10, 3>::VanderMat> Function<1, 10, 3>::VLU_;

template
typename Function<1, 12, 3>::VecOrderD Function<1, 12, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 12, 3>::VanderMat> Function<1, 12, 3>::VLU_;

template
typename Function<1, 14, 3>::VecOrderD Function<1, 14, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 14, 3>::VanderMat> Function<1, 14, 3>::VLU_;

template
typename Function<1, 16, 3>::VecOrderD Function<1, 16, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 16, 3>::VanderMat> Function<1, 16, 3>::VLU_;

template
typename Function<2, 6, 3>::VecOrderD Function<2, 6, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 6, 3>::VanderMat> Function<2, 6, 3>::VLU_;

template
typename Function<2, 8, 3>::VecOrderD Function<2, 8, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 8, 3>::VanderMat> Function<2, 8, 3>::VLU_;

template
typename Function<2, 10, 3>::VecOrderD Function<2, 10, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 10, 3>::VanderMat> Function<2, 10, 3>::VLU_;

template
typename Function<2, 12, 3>::VecOrderD Function<2, 12, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 12, 3>::VanderMat> Function<2, 12, 3>::VLU_;

template
typename Function<2, 14, 3>::VecOrderD Function<2, 14, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 14, 3>::VanderMat> Function<2, 14, 3>::VLU_;

template
typename Function<2, 16, 3>::VecOrderD Function<2, 16, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 16, 3>::VanderMat> Function<2, 16, 3>::VLU_;

template
typename Function<3, 6, 3>::VecOrderD Function<3, 6, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 6, 3>::VanderMat> Function<3, 6, 3>::VLU_;

template
typename Function<3, 8, 3>::VecOrderD Function<3, 8, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 8, 3>::VanderMat> Function<3, 8, 3>::VLU_;

template
typename Function<3, 10, 3>::VecOrderD Function<3, 10, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 10, 3>::VanderMat> Function<3, 10, 3>::VLU_;

template
typename Function<3, 12, 3>::VecOrderD Function<3, 12, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 12, 3>::VanderMat> Function<3, 12, 3>::VLU_;

template
typename Function<3, 14, 3>::VecOrderD Function<3, 14, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 14, 3>::VanderMat> Function<3, 14, 3>::VLU_;

template
typename Function<3, 16, 3>::VecOrderD Function<3, 16, 3>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 16, 3>::VanderMat> Function<3, 16, 3>::VLU_;

}

extern "C" {
BAOBZI_DEFS(1, 6, 3)
BAOBZI_DEFS(1, 8, 3)
BAOBZI_DEFS(1, 10, 3)
BAOBZI_DEFS(1, 12, 3)
BAOBZI_DEFS(1, 14, 3)
BAOBZI_DEFS(1, 16, 3)
BAOBZI_DEFS(2, 6, 3)
BAOBZI_DEFS(2, 8, 3)
BAOBZI_DEFS(2, 10, 3)
BAOBZI_DEFS(2, 12, 3)
BAOBZI_DEFS(2, 14, 3)
BAOBZI_DEFS(2, 16, 3)
BAOBZI_DEFS(3, 6, 3)
BAOBZI_DEFS(3, 8, 3)
BAOBZI_DEFS(3, 10, 3)
BAOBZI_DEFS(3, 12, 3)
BAOBZI_DEFS(3, 14, 3)
BAOBZI_DEFS(3, 16, 3)
}
