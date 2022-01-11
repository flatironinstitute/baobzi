#ifndef BAOBZI_HPP
#define BAOBZI_HPP

#include <fstream>
#include <iostream>
#include <queue>
#include <vector>

#include <msgpack.hpp>
#define EIGEN_MATRIX_PLUGIN "eigen_matrix_plugin.h"

#define EIGEN_MAX_ALIGN_BYTES 64
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <unsupported/Eigen/CXX11/Tensor>

namespace baobzi {

template <int D, int ISET>
struct Box {
    using VEC = Eigen::Vector<double, D>;
    VEC center;
    VEC half_length;

    Box<D, ISET>() = default;
    Box<D, ISET>(const VEC &x, const VEC &hl) : center(x), half_length(hl) {}

    bool contains(const VEC &x) const {
        VEC dx = (x - center).array().abs();
        for (int i = 0; i < D; ++i)
            if (dx[i] > half_length[i])
                return false;
        return true;
    }
    MSGPACK_DEFINE(center, half_length);
};

inline double standard_error(const Eigen::Ref<Eigen::MatrixXd> &coeffs) {
    double maxcoeff = 0.0;
    double scaling_factor = 1.0;
    if (coeffs.cols() == 1) {
        int n = coeffs.size();
        for (auto i = n - 2; i < n; ++i)
            maxcoeff = std::max(std::abs(coeffs(i, 0)), maxcoeff);
        scaling_factor = std::max(scaling_factor, std::abs(coeffs(0, 0)));
    } else {
        int n = coeffs.rows();
        for (auto i = 0; i < n; ++i)
            maxcoeff = std::max(std::abs(coeffs(i, n - i - 1)), maxcoeff);

        scaling_factor = std::max(scaling_factor, std::abs(coeffs(n - 1, 0)));
        scaling_factor = std::max(scaling_factor, std::abs(coeffs(0, n - 1)));
    }

    return maxcoeff / scaling_factor;
}

template <int ORDER, int ISET>
inline double cheb_eval(const Eigen::Vector<double, 1> &x, const Box<1, ISET> &box,
                        const std::vector<double, Eigen::aligned_allocator<double>> &coeffs_raw) {
    double xd = x[0] - box.center[0] / box.half_length[0];

    Eigen::Vector<double, ORDER> Tn;
    Tn[0] = 1.0;
    Tn[1] = xd;
    xd *= 2.0;
    for (int i = 2; i < ORDER; ++i)
        Tn[i] = xd * Tn[i - 1] - Tn[i - 2];

    Eigen::Map<const Eigen::Vector<double, ORDER>> coeffs(coeffs_raw.data());

    return coeffs.dot(Tn);
}

template <int ORDER, int ISET>
inline double cheb_eval(const Eigen::Vector2d &x, const Box<2, ISET> &box,
                        const std::vector<double, Eigen::aligned_allocator<double>> &coeffs_raw) {
    Eigen::Vector2d xinterp = (x - box.center).array() / box.half_length.array();
    Eigen::Matrix<double, 2, ORDER> Tns;
    Tns.col(0).setOnes();
    Tns.col(1) = xinterp;
    xinterp *= 2.0;
    for (int i = 2; i < ORDER; ++i)
        Tns.col(i) = xinterp.array() * Tns.col(i - 1).array() - Tns.col(i - 2).array();

    Eigen::Map<const Eigen::Matrix<double, ORDER, ORDER>> coeffs(coeffs_raw.data());

    return Tns.row(0).transpose().dot(coeffs * Tns.row(1).transpose());
}

template <int ORDER, int ISET>
inline double cheb_eval(const Eigen::Vector3d &x, const Box<3, ISET> &box,
                        const std::vector<double, Eigen::aligned_allocator<double>> &coeffs_raw) {

    Eigen::Vector3d xinterp = (x - box.center).array() / box.half_length.array();

    Eigen::Vector<double, ORDER> Tn[3];
    Tn[0][0] = Tn[1][0] = Tn[2][0] = 1.0;
    for (int i = 0; i < 3; ++i) {
        Tn[i][1] = xinterp[i];
        xinterp[i] *= 2.0;
        for (int j = 2; j < ORDER; ++j)
            Tn[i][j] = xinterp[i] * Tn[i][j - 1] - Tn[i][j - 2];
    }

    double res = 0.0;
    using map_t = Eigen::Map<const Eigen::Matrix<double, ORDER, ORDER>>;
    for (int i = 0; i < ORDER; ++i)
        res += Tn[0][i] * Tn[1].dot(map_t(coeffs_raw.data() + i * ORDER * ORDER) * Tn[2]);

    return res;
}

template <int ORDER, int ISET>
inline double cheb_eval(const Eigen::Vector4d &x, const Box<4, ISET> &box,
                        const std::vector<double, Eigen::aligned_allocator<double>> &coeffs_raw) {

    Eigen::Vector4d xinterp = (x - box.center).array() / box.half_length.array();

    Eigen::Vector<double, ORDER> Tn[4];
    Tn[0][0] = Tn[1][0] = Tn[2][0] = Tn[3][0] = 1.0;
    for (int i = 0; i < 4; ++i) {
        Tn[i][1] = xinterp[i];
        xinterp[i] *= 2.0;
        for (int j = 2; j < ORDER; ++j)
            Tn[i][j] = xinterp[i] * Tn[i][j - 1] - Tn[i][j - 2];
    }

    double res = 0.0;
    for (int i = 0; i < ORDER; ++i)
        for (int j = 0; j < ORDER; ++j)
            for (int k = 0; k < ORDER; ++k)
                for (int l = 0; l < ORDER; ++l)
                    res += Tn[0][i] * Tn[1][j] * Tn[2][k] * Tn[3][l];

    return res;
}

template <int DIM, int ORDER, int ISET>
class Function;

template <int D, int ORDER, int ISET>
class Node {
  public:
    using VEC = Eigen::Vector<double, D>;
    using CoeffVec = Eigen::Vector<double, ORDER>;
    std::vector<double, Eigen::aligned_allocator<double>> coeffs_;
    using Func = Function<D, ORDER, ISET>;
    Box<D, ISET> box_;
    uint64_t first_child_idx = -1;
    bool leaf_ = false;

    Node<D, ORDER, ISET>() = default;

    Node<D, ORDER, ISET>(const Box<D, ISET> &box) : box_(box) {}

    inline bool is_leaf() const { return leaf_; }

    bool fit(double (*f)(const double *), double tol) {
        if constexpr (D == 1) {
            Eigen::Vector<double, ORDER> F;
            CoeffVec xvec =
                Func::get_cheb_nodes(box_.center[0] - box_.half_length[0], box_.center[0] + box_.half_length[0]);

            for (int i = 0; i < ORDER; ++i)
                F(i) = f(&xvec[i]);

            Eigen::Vector<double, ORDER> coeffs = Func::VLU_.solve(F);

            if (standard_error(coeffs) > tol)
                return false;

            coeffs_.resize(coeffs.size());
            for (int i = 0; i < coeffs.size(); ++i)
                coeffs_[i] = coeffs(i);

            leaf_ = true;
            return true;
        }
        if constexpr (D == 2) {
            Eigen::Matrix<double, ORDER, ORDER> F;
            CoeffVec xvec =
                Func::get_cheb_nodes(box_.center[0] - box_.half_length[0], box_.center[0] + box_.half_length[0]);
            CoeffVec yvec =
                Func::get_cheb_nodes(box_.center[1] - box_.half_length[1], box_.center[1] + box_.half_length[1]);

            for (int i = 0; i < ORDER; ++i) {
                for (int j = 0; j < ORDER; ++j) {
                    double x[2] = {xvec[i], yvec[j]};
                    F(i, j) = f(x);
                }
            }

            Eigen::Matrix<double, ORDER, ORDER> coeffs = Func::VLU_.solve(F);
            coeffs = Func::VLU_.solve(coeffs.transpose()).transpose();

            if (standard_error(coeffs) > tol)
                return false;

            coeffs_.resize(coeffs.size());
            for (int i = 0; i < coeffs.size(); ++i)
                coeffs_[i] = coeffs(i);

            leaf_ = true;
            return true;
        }
        if constexpr (D == 3) {
            Eigen::Tensor<double, 3> F(ORDER, ORDER, ORDER);

            CoeffVec xvec =
                Func::get_cheb_nodes(box_.center[0] - box_.half_length[0], box_.center[0] + box_.half_length[0]);
            CoeffVec yvec =
                Func::get_cheb_nodes(box_.center[1] - box_.half_length[1], box_.center[1] + box_.half_length[1]);
            CoeffVec zvec =
                Func::get_cheb_nodes(box_.center[2] - box_.half_length[2], box_.center[2] + box_.half_length[2]);

            for (int i = 0; i < ORDER; ++i) {
                for (int j = 0; j < ORDER; ++j) {
                    for (int k = 0; k < ORDER; ++k) {
                        double x[3] = {xvec[i], yvec[j], zvec[k]};
                        F(i, j, k) = f(x);
                    }
                }
            }

            coeffs_.resize(ORDER * ORDER * ORDER);
            Eigen::Tensor<double, 3> coeffs_tensor(ORDER, ORDER, ORDER);
            using matrix_t = Eigen::Matrix<double, ORDER, ORDER>;
            using matrix_map_t = Eigen::Map<matrix_t>;
            using tensor_t = Eigen::Tensor<double, 2>;
            for (int block = 0; block < ORDER; ++block) {
                tensor_t F_block_tensor = F.chip(block, 2);
                matrix_map_t F_block(F_block_tensor.data());

                matrix_t coeffs_tmp = Func::VLU_.solve(F_block);
                coeffs_tmp = Func::VLU_.solve(coeffs_tmp.transpose()).transpose();
                coeffs_tensor.chip(block, 2) = Eigen::TensorMap<tensor_t>(coeffs_tmp.data(), ORDER, ORDER);
            }
            for (int block = 0; block < ORDER; ++block) {
                Eigen::Tensor<double, 2> coeffs_tmp = coeffs_tensor.chip(block, 0);
                matrix_map_t coeffs_ysolve(coeffs_tmp.data());
                matrix_map_t(coeffs_.data() + block * ORDER * ORDER) = Func::VLU_.solve(coeffs_ysolve.transpose()).transpose();
            }

            for (int i = 0; i < ORDER; ++i) {
                for (int j = 0; j < ORDER; ++j) {
                    for (int k = 0; k < ORDER; ++k) {
                        VEC point =
                            (box_.center - box_.half_length).array() +
                            2.0 * VEC{(double)i, (double)j, (double)k}.array() * box_.half_length.array() / ORDER;

                        const double test_val = eval(point);
                        const double actual_val = f(point.data());
                        const double rel_error = std::abs((actual_val - test_val) / actual_val);

                        if (fabs(actual_val) > 1E-16 && rel_error > tol) {
                            coeffs_.clear();
                            coeffs_.shrink_to_fit();
                            return false;
                        }
                    }
                }
            }

            leaf_ = true;
            return true;
        }
        if constexpr (D == 4) {
            using tensor4_t = Eigen::Tensor<double, 4>;
            tensor4_t F(ORDER, ORDER, ORDER, ORDER);

            Eigen::Matrix<double, ORDER, 4> rvec;
            for (int i = 0; i < 4; ++i)
                rvec.col(i) =
                    Func::get_cheb_nodes(box_.center[i] - box_.half_length[i], box_.center[i] + box_.half_length[i]);

            for (int i = 0; i < ORDER; ++i)
                for (int j = 0; j < ORDER; ++j)
                    for (int k = 0; k < ORDER; ++k)
                        for (int l = 0; l < ORDER; ++l)
                            F(i, j, k, l) = f(rvec.col(i).data());

            coeffs_.resize(ORDER * ORDER * ORDER * ORDER);
            tensor4_t coeffs_tensor(ORDER, ORDER, ORDER, ORDER);
            using matrix_t = Eigen::Matrix<double, ORDER, ORDER>;
            using tensor3_t = Eigen::Tensor<double, 3>;
            using tensor2_t = Eigen::Tensor<double, 2>;

            using matrix_map_t = Eigen::Map<matrix_t>;
            for (int i_block = 0; i_block < ORDER; ++i_block) {
                tensor3_t F_cube_tensor = F.chip(i_block, 3);
                for (int j_block = 0; j_block < ORDER; ++j_block) {
                    tensor2_t F_block_tensor = F_cube_tensor.chip(j_block, 2);
                    matrix_map_t F_block(F_block_tensor.data());

                    matrix_t coeffs_tmp = Func::VLU_.solve(F_block);
                    coeffs_tmp = Func::VLU_.solve(coeffs_tmp.transpose()).transpose();
                    coeffs_tensor.chip(i_block, 3).chip(j_block, 2) =
                        Eigen::TensorMap<tensor2_t>(coeffs_tmp.data(), ORDER, ORDER);
                }
                for (int j_block = 0; j_block < ORDER; ++j_block) {
                    tensor2_t coeffs_tmp = coeffs_tensor.chip(i_block, 3).chip(j_block, 1);
                    matrix_map_t coeffs_ysolve(coeffs_tmp.data());
                    matrix_map_t(coeffs_.data() + i_block * ORDER * ORDER * ORDER) =
                        Func::VLU_.solve(coeffs_ysolve.transpose()).transpose();
                }
            }

            for (int i = 0; i < ORDER; ++i) {
                for (int j = 0; j < ORDER; ++j) {
                    for (int k = 0; k < ORDER; ++k) {
                        for (int l = 0; l < ORDER; ++l) {
                            VEC point =
                                (box_.center - box_.half_length).array() +
                                2.0 * VEC{(double)i, (double)j, (double)k, (double)l}.array() * box_.half_length.array() / ORDER;

                            const double test_val = eval(point);
                            const double actual_val = f(point.data());
                            const double rel_error = std::abs((actual_val - test_val) / actual_val);

                            if (fabs(actual_val) > 1E-16 && rel_error > tol) {
                                coeffs_.clear();
                                coeffs_.shrink_to_fit();
                                return false;
                            }
                        }
                    }
                }
            }

            leaf_ = true;
            return true;
        }        
    }

    inline double eval(const VEC &x) const { return cheb_eval<ORDER, ISET>(x, box_, coeffs_); }
    MSGPACK_DEFINE(box_, first_child_idx, leaf_, coeffs_);
};

template <int DIM, int ORDER, int ISET>
struct FunctionTree {
    static constexpr int NChild = 1 << DIM;
    static constexpr int Dim = DIM;
    static constexpr int Order = ORDER;

    using VEC = Eigen::Vector<double, DIM>;
    std::vector<Node<DIM, ORDER, ISET>> nodes_;

    FunctionTree<DIM, ORDER, ISET>(double (*f)(const double *), const Box<DIM, ISET> &box, double tol) {
        std::queue<Box<DIM, ISET>> q;
        VEC half_width = box.half_length * 0.5;
        q.push(box);

        uint64_t curr_child_idx = 1;
        while (!q.empty()) {
            int n_next = q.size();
            int node_index = nodes_.size();
            for (int i = 0; i < n_next; ++i) {
                Box<DIM, ISET> box = q.front();
                q.pop();

                nodes_.push_back(Node<DIM, ORDER, ISET>(box));
            }

#pragma omp parallel for
            for (size_t i = 0; i < n_next; ++i) {
                auto &node = nodes_[i + node_index];
                node.fit(f, tol);
            }

            for (int i = 0; i < n_next; ++i) {
                auto &node = nodes_[i + node_index];
                if (!node.is_leaf()) {
                    node.first_child_idx = curr_child_idx;
                    curr_child_idx += NChild;

                    VEC &center = node.box_.center;
                    for (uint64_t child = 0; child < NChild; ++child) {
                        VEC offset;

                        // Extract sign of each offset component from the bits of child
                        // Basically: permute all possible offsets
                        for (int j = 0; j < DIM; ++j) {
                            double signed_hw[2] = {-half_width[j], half_width[j]};
                            offset[j] = signed_hw[(child >> j) & 1];
                        }

                        q.push(Box<DIM, ISET>(center + offset, half_width));
                    }
                }
            }

            half_width *= 0.5;
        }
    }

    FunctionTree<DIM, ORDER, ISET>() = default;

    inline const Node<DIM, ORDER, ISET> &find_node_traverse(const VEC &x) const {
        auto *node = &nodes_[0];
        while (!node->is_leaf()) {
            uint64_t child_idx = 0;
            for (int i = 0; i < DIM; ++i)
                child_idx = child_idx | ((x[i] > node->box_.center[i]) << i);

            node = &nodes_[node->first_child_idx + child_idx];
        }

        return *node;
    }

    inline double eval(const VEC &x) const { return find_node_traverse(x).eval(x); }

    MSGPACK_DEFINE(nodes_);
};

template <int DIM, int ORDER, int ISET = 0>
class Function {
  public:
    static constexpr int NChild = 1 << DIM;
    static constexpr int Dim = DIM;
    static constexpr int Order = ORDER;
    static constexpr int ISet = ISET;

    using VEC = Eigen::Vector<double, DIM>;
    using CoeffVec = Eigen::Vector<double, ORDER>;
    using VanderMat = Eigen::Matrix<double, ORDER, ORDER>;

    using DBox = Box<DIM, ISET>;

    static CoeffVec cosarray_;
    static Eigen::PartialPivLU<VanderMat> VLU_;

    double (*f_)(const double *);
    DBox box_;
    double tol_;
    VEC lower_left_;

    std::vector<FunctionTree<DIM, ORDER, ISET>> subtrees_;
    Eigen::Vector<int, DIM> n_subtrees_;
    VEC bin_size_;

    static VanderMat calc_vandermonde() {
        VanderMat V;

        for (int j = 0; j < ORDER; ++j) {
            V(0, j) = 1;
            V(1, j) = cosarray_(j);
        }

        for (int i = 2; i < ORDER; ++i) {
            for (int j = 0; j < ORDER; ++j) {
                V(i, j) = double(2) * V(i - 1, j) * cosarray_(j) - V(i - 2, j);
            }
        }

        return V.transpose();
    }

    static inline CoeffVec get_cheb_nodes(double lb, double ub) {
        return 0.5 * ((lb + ub) + (ub - lb) * cosarray_.array());
    }

    static void init_statics() {
        for (int i = 0; i < ORDER; ++i)
            cosarray_[ORDER - i - 1] = cos(M_PI * (i + 0.5) / ORDER);
        VLU_ = Eigen::PartialPivLU<VanderMat>(calc_vandermonde());
    }

    Function<DIM, ORDER, ISET>(double (*f)(const double *), const double *xp, const double *lp, double tol)
        : f_(f), box_(VEC(xp), VEC(lp)), tol_(tol) {
        init_statics();

        VEC l(lp);
        VEC x(xp);
        std::queue<DBox> q;

        for (int i = 0; i < DIM; ++i)
            n_subtrees_[i] = l[i] / l.minCoeff();

        uint8_t max_depth_ = 0;
        q.push(DBox(x, l));

        // Half-width of next children
        VEC half_width = l * 0.5;
        while (!q.empty()) {
            int n_next = q.size();

            std::vector<Node<DIM, ORDER, ISET>> nodes;
            for (int i = 0; i < n_next; ++i) {
                DBox box = q.front();
                q.pop();

                nodes.emplace_back(Node<DIM, ORDER, ISET>(box));
            }

#pragma omp parallel for
            for (auto &node : nodes)
                node.fit(f, tol);

            for (auto &node : nodes) {
                if (!node.is_leaf()) {
                    VEC &center = node.box_.center;
                    for (unsigned child = 0; child < NChild; ++child) {
                        VEC offset;

                        // Extract sign of each offset component from the bits of child
                        // Basically: permute all possible offsets
                        for (int j = 0; j < DIM; ++j) {
                            double signed_hw[2] = {-half_width[j], half_width[j]};
                            offset[j] = signed_hw[(child >> j) & 1];
                        }

                        q.push(DBox(center + offset, half_width));
                    }
                }
            }

            if (!q.empty())
                max_depth_++;

            half_width *= 0.5;
            if ((1 << (DIM * max_depth_)) == q.size())
                n_subtrees_ *= 2;
            else
                break;
        }

        for (int j = 0; j < DIM; ++j)
            bin_size_[j] = 2.0 * box_.half_length[j] / n_subtrees_[j];
        lower_left_ = box_.center - box_.half_length;

        subtrees_.reserve(n_subtrees_.prod());
        for (int i_bin = 0; i_bin < n_subtrees_.prod(); ++i_bin) {
            Eigen::Vector<int, DIM> bins = get_bins(i_bin);

            VEC parent_center = (bins.template cast<double>().array() + 0.5) * bin_size_.array() + lower_left_.array();

            Box<DIM, ISET> root_box = {parent_center, 0.5 * bin_size_};
            subtrees_.push_back(FunctionTree<DIM, ORDER, ISET>(f_, root_box, tol_));
        }
    }

    Function<DIM, ORDER, ISET>() { init_statics(); };

    inline Eigen::Vector<int, DIM> get_bins(const int i_bin) const {
        if constexpr (DIM == 1)
            return Eigen::Vector<int, DIM>{i_bin};
        else if constexpr (DIM == 2)
            return Eigen::Vector<int, DIM>{i_bin % n_subtrees_[0], i_bin / n_subtrees_[0]};
        else if constexpr (DIM == 3)
            return Eigen::Vector<int, DIM>{i_bin % n_subtrees_[0], (i_bin / n_subtrees_[0]) % n_subtrees_[1],
                                           i_bin / (n_subtrees_[0] * n_subtrees_[1])};
        else if constexpr (DIM == 4)
            return Eigen::Vector<int, DIM>{i_bin % n_subtrees_[0], (i_bin / n_subtrees_[0]) % n_subtrees_[1],
                                           i_bin / (n_subtrees_[0] * n_subtrees_[1]) % n_subtrees_[2],
                                           i_bin / (n_subtrees_[0] * n_subtrees_[1] * n_subtrees_[2])};
    }

    inline int get_linear_bin(const Eigen::Vector<double, 1> &x) const {
        const double x_bin = x[0] - lower_left_[0];
        return x_bin / bin_size_[0];
    }

    inline int get_linear_bin(const Eigen::Vector2d &x) const {
        const VEC x_bin = x - lower_left_;
        const Eigen::Vector<int, DIM> bin = (x_bin.array() / bin_size_.array()).template cast<int>();
        return bin[0] + n_subtrees_[0] * bin[1];
    }

    inline int get_linear_bin(const Eigen::Vector3d &x) const {
        const VEC x_bin = x - lower_left_;
        const Eigen::Vector<int, DIM> bin = (x_bin.array() / bin_size_.array()).template cast<int>();
        return bin[0] + n_subtrees_[0] * bin[1] + n_subtrees_[0] * n_subtrees_[1] * bin[2];
    }

    inline int get_linear_bin(const Eigen::Vector4d &x) const {
        const VEC x_bin = x - lower_left_;
        const Eigen::Vector<int, DIM> bin = (x_bin.array() / bin_size_.array()).template cast<int>();
        return bin[0] + n_subtrees_[0] * bin[1] + n_subtrees_[0] * n_subtrees_[1] * bin[2] +
               n_subtrees_[0] * n_subtrees_[1] * n_subtrees_[2] * bin[3];
    }

    inline const Node<DIM, ORDER, ISET> &find_node(const VEC &x) const {
        return subtrees_[get_linear_bin(x)].find_node_traverse(x);
    }

    inline double eval(const VEC &x) const { return find_node(x).eval(x); }
    inline double eval(const double *xp) const { return eval(VEC(xp)); }

    inline double operator()(const VEC &x) const { return eval(x); }
    inline double operator()(const double *x) const { return eval(x); }

    void save(const char *filename) {
        std::ofstream ofs(filename, std::ofstream::binary | std::ofstream::out);
        std::array<int, 2> params{Dim, Order};
        msgpack::pack(ofs, params);
        msgpack::pack(ofs, *this);
    }

    MSGPACK_DEFINE_MAP(box_, subtrees_, n_subtrees_, tol_, lower_left_, bin_size_);
};

template <int DIM, int ORDER, int ISET>
typename Function<DIM, ORDER, ISET>::CoeffVec Function<DIM, ORDER, ISET>::cosarray_;

template <int DIM, int ORDER, int ISET>
Eigen::PartialPivLU<typename Function<DIM, ORDER, ISET>::VanderMat> Function<DIM, ORDER, ISET>::VLU_;
} // namespace baobzi

#endif
