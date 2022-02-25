#include "baobzi_template.hpp"
#include "baobzi.h"
#include "baobzi/macros.h"

namespace baobzi {
template
typename Function<1, 6, 2>::VecOrderD Function<1, 6, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 6, 2>::VanderMat> Function<1, 6, 2>::VLU_;

template
typename Function<1, 8, 2>::VecOrderD Function<1, 8, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 8, 2>::VanderMat> Function<1, 8, 2>::VLU_;

template
typename Function<1, 10, 2>::VecOrderD Function<1, 10, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 10, 2>::VanderMat> Function<1, 10, 2>::VLU_;

template
typename Function<1, 12, 2>::VecOrderD Function<1, 12, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 12, 2>::VanderMat> Function<1, 12, 2>::VLU_;

template
typename Function<1, 14, 2>::VecOrderD Function<1, 14, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 14, 2>::VanderMat> Function<1, 14, 2>::VLU_;

template
typename Function<1, 16, 2>::VecOrderD Function<1, 16, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<1, 16, 2>::VanderMat> Function<1, 16, 2>::VLU_;

template
typename Function<2, 6, 2>::VecOrderD Function<2, 6, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 6, 2>::VanderMat> Function<2, 6, 2>::VLU_;

template
typename Function<2, 8, 2>::VecOrderD Function<2, 8, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 8, 2>::VanderMat> Function<2, 8, 2>::VLU_;

template
typename Function<2, 10, 2>::VecOrderD Function<2, 10, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 10, 2>::VanderMat> Function<2, 10, 2>::VLU_;

template
typename Function<2, 12, 2>::VecOrderD Function<2, 12, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 12, 2>::VanderMat> Function<2, 12, 2>::VLU_;

template
typename Function<2, 14, 2>::VecOrderD Function<2, 14, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 14, 2>::VanderMat> Function<2, 14, 2>::VLU_;

template
typename Function<2, 16, 2>::VecOrderD Function<2, 16, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<2, 16, 2>::VanderMat> Function<2, 16, 2>::VLU_;

template
typename Function<3, 6, 2>::VecOrderD Function<3, 6, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 6, 2>::VanderMat> Function<3, 6, 2>::VLU_;

template
typename Function<3, 8, 2>::VecOrderD Function<3, 8, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 8, 2>::VanderMat> Function<3, 8, 2>::VLU_;

template
typename Function<3, 10, 2>::VecOrderD Function<3, 10, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 10, 2>::VanderMat> Function<3, 10, 2>::VLU_;

template
typename Function<3, 12, 2>::VecOrderD Function<3, 12, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 12, 2>::VanderMat> Function<3, 12, 2>::VLU_;

template
typename Function<3, 14, 2>::VecOrderD Function<3, 14, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 14, 2>::VanderMat> Function<3, 14, 2>::VLU_;

template
typename Function<3, 16, 2>::VecOrderD Function<3, 16, 2>::cosarray_;
template
Eigen::PartialPivLU<typename Function<3, 16, 2>::VanderMat> Function<3, 16, 2>::VLU_;

}

extern "C" {
BAOBZI_DEFS(1, 6, 2)
BAOBZI_DEFS(1, 8, 2)
BAOBZI_DEFS(1, 10, 2)
BAOBZI_DEFS(1, 12, 2)
BAOBZI_DEFS(1, 14, 2)
BAOBZI_DEFS(1, 16, 2)
BAOBZI_DEFS(2, 6, 2)
BAOBZI_DEFS(2, 8, 2)
BAOBZI_DEFS(2, 10, 2)
BAOBZI_DEFS(2, 12, 2)
BAOBZI_DEFS(2, 14, 2)
BAOBZI_DEFS(2, 16, 2)
BAOBZI_DEFS(3, 6, 2)
BAOBZI_DEFS(3, 8, 2)
BAOBZI_DEFS(3, 10, 2)
BAOBZI_DEFS(3, 12, 2)
BAOBZI_DEFS(3, 14, 2)
BAOBZI_DEFS(3, 16, 2)
}
