#ifndef BAOBZI_HPP
#define BAOBZI_HPP

#include <queue>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

namespace baobzi {

template <int ORDER>
struct chebyshev_nodes {
    using arr_type = Eigen::Vector<double, ORDER>;
    static const arr_type cosarray;

    static arr_type get_cos_array() {
        arr_type cosarray(ORDER);
        for (int i = 0; i < ORDER; ++i)
            cosarray[ORDER - i - 1] = cos(M_PI * (i + 0.5) / ORDER);
        return cosarray;
    }

    static inline arr_type get(double lb, double ub) { return 0.5 * ((lb + ub) + (ub - lb) * cosarray.array()); }
};

template <int ORDER>
const typename chebyshev_nodes<ORDER>::arr_type
    chebyshev_nodes<ORDER>::cosarray = chebyshev_nodes<ORDER>::get_cos_array();

template <int D>
struct Box {
    using VEC = Eigen::Vector<double, D>;
    VEC center;
    VEC half_length;

    Box<D>() = default;
    Box<D>(const VEC &x, const VEC &hl) : center(x), half_length(hl) {}

    bool contains(const VEC &x) const {
        VEC dx = (x - center).array().abs();
        for (int i = 0; i < D; ++i)
            if (dx[i] > half_length[i])
                return false;
        return true;
    }
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

template <int ORDER>
inline double cheb_eval(const Eigen::Vector<double, 1> &x, const Box<1> &box,
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

template <int ORDER>
inline double cheb_eval(const Eigen::Vector2d &x, const Box<2> &box,
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

template <int ORDER>
inline double cheb_eval(const Eigen::Vector3d &x, const Box<3> &box,
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

template <int D, int ORDER>
class Node {
  public:
    using VEC = Eigen::Vector<double, D>;
    using VanderMat = Eigen::Matrix<double, ORDER, ORDER>;
    using CoeffVec = Eigen::Vector<double, ORDER>;
    static const Eigen::PartialPivLU<VanderMat> VLU_;
    std::vector<double, Eigen::aligned_allocator<double>> coeffs_;
    Box<D> box_;
    uint64_t first_child_idx;
    bool leaf_ = false;

    static VanderMat calc_vandermonde() {
        VanderMat V;

        auto x = chebyshev_nodes<ORDER>::get_cos_array();
        for (int j = 0; j < ORDER; ++j) {
            V(0, j) = 1;
            V(1, j) = x(j);
        }

        for (int i = 2; i < ORDER; ++i) {
            for (int j = 0; j < ORDER; ++j) {
                V(i, j) = double(2) * V(i - 1, j) * x(j) - V(i - 2, j);
            }
        }
        V = V.transpose().eval();
        return V;
    }

    Node<D, ORDER>() = default;

    Node<D, ORDER>(const Box<D> &box) : box_(box) {}

    inline bool is_leaf() const { return leaf_; }

    bool fit(double (*f)(const double *), double tol) {
        if constexpr (D == 1) {
            Eigen::Vector<double, ORDER> F;
            CoeffVec xvec =
                chebyshev_nodes<ORDER>::get(box_.center[0] - box_.half_length[0], box_.center[0] + box_.half_length[0]);

            for (int i = 0; i < ORDER; ++i)
                F(i) = f(&xvec[i]);

            Eigen::Vector<double, ORDER> coeffs = VLU_.solve(F);

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
                chebyshev_nodes<ORDER>::get(box_.center[0] - box_.half_length[0], box_.center[0] + box_.half_length[0]);
            CoeffVec yvec =
                chebyshev_nodes<ORDER>::get(box_.center[1] - box_.half_length[1], box_.center[1] + box_.half_length[1]);

            for (int i = 0; i < ORDER; ++i) {
                for (int j = 0; j < ORDER; ++j) {
                    double x[2] = {xvec[i], yvec[j]};
                    F(i, j) = f(x);
                }
            }

            Eigen::Matrix<double, ORDER, ORDER> coeffs = VLU_.solve(F);
            coeffs = VLU_.solve(coeffs.transpose()).transpose();

            if (standard_error(coeffs) > tol)
                return false;

            coeffs_.resize(coeffs.size());
            for (int i = 0; i < coeffs.size(); ++i)
                coeffs_[i] = coeffs(i);

            leaf_ = true;
            return true;
        }
        if constexpr (D == 3) {
            std::vector<double, Eigen::aligned_allocator<double>> F(ORDER * ORDER * ORDER);
            CoeffVec xvec =
                chebyshev_nodes<ORDER>::get(box_.center[0] - box_.half_length[0], box_.center[0] + box_.half_length[0]);
            CoeffVec yvec =
                chebyshev_nodes<ORDER>::get(box_.center[1] - box_.half_length[1], box_.center[1] + box_.half_length[1]);
            CoeffVec zvec =
                chebyshev_nodes<ORDER>::get(box_.center[2] - box_.half_length[2], box_.center[2] + box_.half_length[2]);

            for (int i = 0; i < ORDER; ++i) {
                for (int j = 0; j < ORDER; ++j) {
                    for (int k = 0; k < ORDER; ++k) {
                        double x[3] = {xvec[i], yvec[j], zvec[k]};
                        F[i * ORDER * ORDER + j * ORDER + k] = f(x);
                    }
                }
            }

            std::vector<double, Eigen::aligned_allocator<double>> coeffs_z(ORDER * ORDER * ORDER);
            std::vector<double, Eigen::aligned_allocator<double>> coeffs_y(ORDER * ORDER * ORDER);
            std::vector<double, Eigen::aligned_allocator<double>> coeffs_x(ORDER * ORDER * ORDER);
            for (int block = 0; block < ORDER; ++block) {
                using matrix_t = Eigen::Matrix<double, ORDER, ORDER>;
                using map_t = Eigen::Map<matrix_t>;
                int block_offset = block * ORDER * ORDER;
                map_t F_block(F.data() + block_offset);
                map_t coeffs_zsolve(coeffs_z.data() + block_offset);
                map_t coeffs_ysolve(coeffs_y.data() + block_offset);
                coeffs_zsolve = VLU_.solve(F_block);
                coeffs_ysolve = VLU_.solve(coeffs_zsolve.transpose()).transpose();
            }
            for (int block = 0; block < ORDER; ++block) {
                using matrix_t = Eigen::Matrix<double, ORDER, ORDER>;
                using map_t = Eigen::Map<matrix_t, 0, Eigen::OuterStride<ORDER * ORDER>>;
                int block_offset = block * ORDER;
                matrix_t coeffs_ysolve = map_t(coeffs_y.data() + block_offset);

                Eigen::Map<matrix_t> coeffs_xsolve(coeffs_x.data() + block * ORDER * ORDER);
                coeffs_xsolve = VLU_.solve(coeffs_ysolve.transpose()).transpose();
            }

            coeffs_ = coeffs_x;

            for (int i = 0; i < ORDER; ++i) {
                for (int j = 0; j < ORDER; ++j) {
                    for (int k = 0; k < ORDER; ++k) {
                        VEC point =
                            (box_.center - box_.half_length).array() +
                            2.0 * VEC{(double)i, (double)j, (double)k}.array() * box_.half_length.array() / ORDER;

                        double test_val = this->eval(point);
                        double actual_val = f(point.data());
                        double rel_error = std::abs((actual_val - test_val) / test_val);

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
    }

    inline double eval(const VEC &x) const { return cheb_eval<ORDER>(x, box_, coeffs_); }
};

template <int DIM, int ORDER>
struct FunctionTree {
    static constexpr int NChild = 1 << DIM;
    static constexpr int Dim = DIM;
    static constexpr int Order = ORDER;

    using VEC = Eigen::Vector<double, DIM>;
    std::vector<Node<DIM, ORDER>> nodes_;

    FunctionTree<DIM, ORDER>(double (*f)(const double *), const Box<DIM> &box, double tol) {
        std::queue<Box<DIM>> q;
        VEC half_width = box.half_length * 0.5;
        q.push(box);

        uint64_t curr_child_idx = 1;
        while (!q.empty()) {
            int n_next = q.size();

            for (int i = 0; i < n_next; ++i) {
                Box<DIM> box = q.front();
                q.pop();

                nodes_.push_back(Node<DIM, ORDER>(box));
                auto &node = nodes_.back();
                if (!node.fit(f, tol)) {
                    node.first_child_idx = curr_child_idx;
                    curr_child_idx += NChild;

                    VEC &center = box.center;
                    for (uint64_t child = 0; child < NChild; ++child) {
                        VEC offset;

                        // Extract sign of each offset component from the bits of child
                        // Basically: permute all possible offsets
                        for (int j = 0; j < DIM; ++j) {
                            double signed_hw[2] = {-half_width[j], half_width[j]};
                            offset[j] = signed_hw[(child >> j) & 1];
                        }

                        q.push(Box<DIM>(center + offset, half_width));
                    }
                }
            }

            half_width *= 0.5;
        }
    }

    inline const Node<DIM, ORDER> &find_node_traverse(const VEC &x) const {
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
};

template <int DIM, int ORDER>
class Function {
  public:
    static constexpr int NChild = 1 << DIM;
    static constexpr int Dim = DIM;
    static constexpr int Order = ORDER;

    using VEC = Eigen::Vector<double, DIM>;
    using CoeffVec = Eigen::Vector<double, ORDER>;
    using DBox = Box<DIM>;

    const double tol_;
    double (*f_)(const double *);
    DBox box_;
    VEC lower_left_;

    std::vector<FunctionTree<DIM, ORDER>> subtrees_;
    Eigen::Vector<int, DIM> n_subtrees_;
    VEC bin_size_;

    Function<DIM, ORDER>(double (*f)(const double *), const double *xp, const double *lp, double tol)
        : f_(f), box_(VEC(xp), VEC(lp)), tol_(tol) {
        VEC l(lp);
        VEC x(xp);
        std::queue<DBox> q;

        for (int i = 0; i < DIM; ++i)
            n_subtrees_[i] = l[i] / l.minCoeff();

        uint8_t max_depth_ = 0;
        q.push(DBox(x, l));

        bool fill_subtrees = false;

        // Half-width of next children
        VEC half_width = l * 0.5;
        while (!q.empty()) {
            int n_next = q.size();

            for (int i = 0; i < n_next; ++i) {
                DBox box = q.front();
                q.pop();

                Node<DIM, ORDER> node(box);
                if (!node.fit(f, tol_)) {
                    VEC &center = box.center;
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

            Box<DIM> root_box = {parent_center, 0.5 * bin_size_};
            subtrees_.push_back(FunctionTree<DIM, ORDER>(f_, root_box, tol_));
        }
    }

    inline Eigen::Vector<int, DIM> get_bins(const int i_bin) const {
        if constexpr (DIM == 1)
            return Eigen::Vector<int, DIM>{i_bin};
        else if constexpr (DIM == 2)
            return Eigen::Vector<int, DIM>{i_bin % n_subtrees_[0], i_bin / n_subtrees_[0]};
        else if constexpr (DIM == 3)
            return Eigen::Vector<int, DIM>{i_bin % n_subtrees_[0], (i_bin / n_subtrees_[0]) % n_subtrees_[1],
                                           i_bin / (n_subtrees_[0] * n_subtrees_[1])};
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

    inline const Node<DIM, ORDER> &find_node(const VEC &x) const {
        return subtrees_[get_linear_bin(x)].find_node_traverse(x);
    }

    inline double eval(const VEC &x) const { return find_node(x).eval(x); }
    inline double eval(const double *xp) const { return eval(VEC(xp)); }

    inline double operator()(const VEC &x) const { return eval(x); }
    inline double operator()(const double *x) const { return eval(x); }
};

template <int D, int ORDER>
const Eigen::PartialPivLU<typename Node<D, ORDER>::VanderMat>
    Node<D, ORDER>::VLU_ = Eigen::PartialPivLU<typename Node<D, ORDER>::VanderMat>(Node<D, ORDER>::calc_vandermonde());
} // namespace baobzi

#endif
