#include "baobzi_template.hpp"
#include "baobzi.h"
#include "baobzi/macros.h"

namespace baobzi {
template
typename Function<1, 6, 1>::VecOrderD Function<1, 6, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 6, 1>::VanderMat> Function<1, 6, 1>::VLU_;

template
typename Function<1, 8, 1>::VecOrderD Function<1, 8, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 8, 1>::VanderMat> Function<1, 8, 1>::VLU_;

template
typename Function<1, 10, 1>::VecOrderD Function<1, 10, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 10, 1>::VanderMat> Function<1, 10, 1>::VLU_;

template
typename Function<1, 12, 1>::VecOrderD Function<1, 12, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 12, 1>::VanderMat> Function<1, 12, 1>::VLU_;

template
typename Function<1, 14, 1>::VecOrderD Function<1, 14, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 14, 1>::VanderMat> Function<1, 14, 1>::VLU_;

template
typename Function<1, 16, 1>::VecOrderD Function<1, 16, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 16, 1>::VanderMat> Function<1, 16, 1>::VLU_;

template
typename Function<2, 6, 1>::VecOrderD Function<2, 6, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 6, 1>::VanderMat> Function<2, 6, 1>::VLU_;

template
typename Function<2, 8, 1>::VecOrderD Function<2, 8, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 8, 1>::VanderMat> Function<2, 8, 1>::VLU_;

template
typename Function<2, 10, 1>::VecOrderD Function<2, 10, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 10, 1>::VanderMat> Function<2, 10, 1>::VLU_;

template
typename Function<2, 12, 1>::VecOrderD Function<2, 12, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 12, 1>::VanderMat> Function<2, 12, 1>::VLU_;

template
typename Function<2, 14, 1>::VecOrderD Function<2, 14, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 14, 1>::VanderMat> Function<2, 14, 1>::VLU_;

template
typename Function<2, 16, 1>::VecOrderD Function<2, 16, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 16, 1>::VanderMat> Function<2, 16, 1>::VLU_;

template
typename Function<3, 6, 1>::VecOrderD Function<3, 6, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 6, 1>::VanderMat> Function<3, 6, 1>::VLU_;

template
typename Function<3, 8, 1>::VecOrderD Function<3, 8, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 8, 1>::VanderMat> Function<3, 8, 1>::VLU_;

template
typename Function<3, 10, 1>::VecOrderD Function<3, 10, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 10, 1>::VanderMat> Function<3, 10, 1>::VLU_;

template
typename Function<3, 12, 1>::VecOrderD Function<3, 12, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 12, 1>::VanderMat> Function<3, 12, 1>::VLU_;

template
typename Function<3, 14, 1>::VecOrderD Function<3, 14, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 14, 1>::VanderMat> Function<3, 14, 1>::VLU_;

template
typename Function<3, 16, 1>::VecOrderD Function<3, 16, 1>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 16, 1>::VanderMat> Function<3, 16, 1>::VLU_;

}

extern "C" {
BAOBZI_DEFS(1, 6, 1)
BAOBZI_DEFS(1, 8, 1)
BAOBZI_DEFS(1, 10, 1)
BAOBZI_DEFS(1, 12, 1)
BAOBZI_DEFS(1, 14, 1)
BAOBZI_DEFS(1, 16, 1)
BAOBZI_DEFS(2, 6, 1)
BAOBZI_DEFS(2, 8, 1)
BAOBZI_DEFS(2, 10, 1)
BAOBZI_DEFS(2, 12, 1)
BAOBZI_DEFS(2, 14, 1)
BAOBZI_DEFS(2, 16, 1)
BAOBZI_DEFS(3, 6, 1)
BAOBZI_DEFS(3, 8, 1)
BAOBZI_DEFS(3, 10, 1)
BAOBZI_DEFS(3, 12, 1)
BAOBZI_DEFS(3, 14, 1)
BAOBZI_DEFS(3, 16, 1)
}
