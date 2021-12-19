#ifndef BAOBZI_HPP
#define BAOBZI_HPP

#include <queue>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/LU>

namespace baobzi {

static std::array<uint64_t, 5> dilate_masks[4] = {
    {},
    {},
    {
        0x0000FFFF0000FFFF,
        0x00FF00FF00FF00FF,
        0x0F0F0F0F0F0F0F0F,
        0x3333333333333333,
        0x5555555555555555,
    },
    {
        0xFFFF00000000FFFF,
        0x00FF0000FF0000FF,
        0xF00F00F00F00F00F,
        0x30C30C30C30C30C3,
        0x9249249249249249,
    },
};

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

template <int D, int ORDER>
class Node {
  public:
    using VEC = Eigen::Vector<double, D>;
    using VanderMat = Eigen::Matrix<double, ORDER, ORDER>;
    using CoeffVec = Eigen::Vector<double, ORDER>;
    static const Eigen::PartialPivLU<VanderMat> VLU_;
    std::vector<double> coeffs_;
    Box<D> box_;
    uint64_t first_child_idx;

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

    bool is_leaf() const { return coeffs_.size(); }

    bool fit(double (*f)(VEC), double tol) {
        if constexpr (D == 1) {
            Eigen::Vector<double, ORDER> F;
            CoeffVec xvec =
                chebyshev_nodes<ORDER>::get(box_.center[0] - box_.half_length[0], box_.center[0] + box_.half_length[0]);

            for (int i = 0; i < ORDER; ++i)
                F(i) = f(Eigen::Vector<double, 1>{xvec[i]});

            Eigen::Vector<double, ORDER> coeffs = VLU_.solve(F);

            if (standard_error(coeffs) > tol)
                return false;

            coeffs_.resize(coeffs.size());
            for (int i = 0; i < coeffs.size(); ++i)
                coeffs_[i] = coeffs(i);

            return true;
        }
        if constexpr (D == 2) {
            Eigen::Matrix<double, ORDER, ORDER> F;
            CoeffVec xvec =
                chebyshev_nodes<ORDER>::get(box_.center[0] - box_.half_length[0], box_.center[0] + box_.half_length[0]);
            CoeffVec yvec =
                chebyshev_nodes<ORDER>::get(box_.center[1] - box_.half_length[1], box_.center[1] + box_.half_length[1]);

            for (int i = 0; i < ORDER; ++i)
                for (int j = 0; j < ORDER; ++j)
                    F(i, j) = f({xvec[i], yvec[j]});

            Eigen::Matrix<double, ORDER, ORDER> coeffs = VLU_.solve(F);
            coeffs = VLU_.solve(coeffs.transpose()).transpose();

            if (standard_error(coeffs) > tol)
                return false;

            coeffs_.resize(coeffs.size());
            for (int i = 0; i < coeffs.size(); ++i)
                coeffs_[i] = coeffs(i);

            return true;
        }
        if constexpr (D == 3) {
            std::vector<double> F(ORDER * ORDER * ORDER);
            CoeffVec xvec =
                chebyshev_nodes<ORDER>::get(box_.center[0] - box_.half_length[0], box_.center[0] + box_.half_length[0]);
            CoeffVec yvec =
                chebyshev_nodes<ORDER>::get(box_.center[1] - box_.half_length[1], box_.center[1] + box_.half_length[1]);
            CoeffVec zvec =
                chebyshev_nodes<ORDER>::get(box_.center[2] - box_.half_length[2], box_.center[2] + box_.half_length[2]);

            for (int i = 0; i < ORDER; ++i)
                for (int j = 0; j < ORDER; ++j)
                    for (int k = 0; k < ORDER; ++k)
                        F[i * ORDER * ORDER + j * ORDER + k] = f({xvec[i], yvec[j], zvec[k]});

            std::vector<double> coeffs_z(ORDER * ORDER * ORDER);
            std::vector<double> coeffs_y(ORDER * ORDER * ORDER);
            std::vector<double> coeffs_x(ORDER * ORDER * ORDER);
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
                        VEC point = box_.center - box_.half_length + i * box_.half_length / ORDER;

                        double test_val = this->eval(point);
                        double actual_val = f(point);
                        double rel_error = std::abs(this->eval(point) - f(point));

                        if (fabs(actual_val) > 1E-16 && rel_error > tol) {
                            coeffs_.clear();
                            coeffs_.shrink_to_fit();
                            return false;
                        }
                    }
                }
            }

            return true;
        }
    }

    double eval(const VEC &x) const {
        VEC xinterp = (x - box_.center).array() / box_.half_length.array();

        if constexpr (D == 1) {
            double xd = xinterp[0];
            CoeffVec Tn;
            Tn[0] = 1.0;
            Tn[1] = xd;
            xd *= 2.0;
            for (int i = 2; i < ORDER; ++i)
                Tn[i] = xd * Tn[i - 1] - Tn[i - 2];

            Eigen::Map<const Eigen::Vector<double, ORDER>> coeffs(coeffs_.data());

            return coeffs.dot(Tn);
        }
        if constexpr (D == 2) {
            CoeffVec Tnx, Tny;
            Eigen::Matrix<double, 2, ORDER> Tns;
            Tns.col(0).setOnes();
            Tns.col(1) = xinterp;
            xinterp *= 2.0;
            for (int i = 2; i < ORDER; ++i)
                Tns.col(i) = xinterp.array() * Tns.col(i - 1).array() - Tns.col(i - 2).array();

            Eigen::Map<const Eigen::Matrix<double, ORDER, ORDER>> coeffs(coeffs_.data());

            return Tns.row(0).transpose().dot(coeffs * Tns.row(1).transpose());
        }
        if constexpr (D == 3) {
            CoeffVec Tn[3];
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
                res += Tn[0][i] * Tn[1].dot(map_t(coeffs_.data() + i * ORDER * ORDER) * Tn[2]);

            return res;
        }
    }
};

template <int b>
struct Int64Hash {
    size_t operator()(const uint64_t &t) const {
        constexpr uint64_t mask = (1 << b) - 1;
        return t & mask;
    }
};

template <int DIM>
inline uint64_t calculate_key(const Eigen::Vector<double, DIM> &x_scaled, int depth) {
    const std::array<uint64_t, 5> &dilate_mask = dilate_masks[DIM];

    // root node always one
    uint64_t key = 1 << (depth * DIM);
    for (int dim = 0; dim < DIM; ++dim) {
        uint32_t dim_key = x_scaled[dim] * (1 << depth);
        constexpr unsigned int base_offset = sizeof(dim_key) * 8 * (DIM - 1);
        for (int i = 0; i < 5; ++i)
            dim_key = (dim_key | (dim_key << (base_offset >> (i + 1)))) & dilate_mask[i];
        key = key | (dim_key << dim);
    }

    return key;
}

template <int DIM, int ORDER>
struct FunctionTree {
    static constexpr int NChild = 1 << DIM;
    static constexpr int Dim = DIM;
    static constexpr int Order = ORDER;

    using VEC = Eigen::Vector<double, DIM>;
    Box<DIM> box_;
    std::vector<Node<DIM, ORDER>> nodes_;

    FunctionTree<DIM, ORDER>() : box_(Box<DIM>()){};

    FunctionTree<DIM, ORDER>(const Box<DIM> &box, double (*f)(VEC), double tol) : box_(box) {
        using key_t = uint64_t;
        std::queue<std::pair<Box<DIM>, key_t>> q;

        VEC half_width = box.half_length * 0.5;

        q.push(std::make_pair(box, 1));

        uint64_t parent_idx = 0;
        while (!q.empty()) {
            int n_next = q.size();

            uint64_t curr_child_idx = parent_idx + q.size();

            for (int i = 0; i < n_next; ++i) {
                auto [box, node_key] = q.front();
                q.pop();
                Node<DIM, ORDER> test_node(box);
                nodes_.push_back(Node<DIM, ORDER>(box));

                auto &node = nodes_.back();
                if (!node.fit(f, tol)) {
                    node.first_child_idx = curr_child_idx;
                    curr_child_idx += NChild;

                    key_t parent = node_key << DIM;
                    VEC &center = box.center;
                    for (key_t child = 0; child < NChild; ++child) {
                        VEC offset;

                        // Extract sign of each offset component from the bits of child
                        // Basically: permute all possible offsets
                        for (int j = 0; j < DIM; ++j) {
                            double signed_hw[2] = {-half_width[j], half_width[j]};
                            offset[j] = signed_hw[(child >> j) & 1];
                        }

                        q.push(std::make_pair(Box<DIM>(center + offset, half_width), parent + child));
                    }
                }
            }

            half_width *= 0.5;
        }
    }

    const Node<DIM, ORDER> &find_node_traverse(const VEC &x) const {
        auto *node = &nodes_[0];
        while (!node->is_leaf()) {
            uint64_t child_idx = 0;
            for (int i = 0; i < DIM; ++i)
                child_idx = child_idx | ((x[i] > node->box_.center[i]) << i);

            node = &nodes_[node->first_child_idx + child_idx];
        }

        return *node;
    }

    double eval(const VEC &x) { return find_node_traverse(x).eval(x); }

    // uint64_t find_node_key(const VEC &x) const {
    //     const VEC &center = box_.center;
    //     const VEC &half_width = box_.half_length;

    //     const VEC x_scaled = 0.5 * (x - center).array() / box_.half_length.array() + VEC::Constant(0.5).array();
    //     uint64_t m_max = calculate_key(x_scaled, max_depth_);

    //     uint64_t rel_depth = (max_depth_ - max_full_) * DIM;
    //     // Our guess is below the actual depth of the target node, move up tree until we hit it
    //     while (flat_map_[m_max >> (rel_depth - DIM)] != std::numeric_limits<uint64_t>::max())
    //         rel_depth -= DIM;

    //     return m_max >> rel_depth;
    // }
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
    double (*f_)(VEC);
    DBox box_;
    VEC lower_left_;

    std::vector<FunctionTree<DIM, ORDER>> subtrees_;
    Eigen::Vector<int, DIM> n_subtrees_;
    VEC bin_size_;

    Function<DIM, ORDER>(const VEC &x, const VEC &l, double (*f)(VEC), double tol) : f_(f), box_(x, l), tol_(tol) {
        using key_t = uint64_t;
        std::queue<std::pair<DBox, key_t>> q;

        for (int i = 0; i < DIM; ++i)
            n_subtrees_[i] = l[i] / l.minCoeff();

        uint8_t max_depth_ = 0;
        q.push(std::make_pair(DBox(x, l), 1));

        bool fill_subtrees = false;

        // Half-width of next children
        VEC half_width = l * 0.5;
        uint64_t parent_idx = 0;
        std::vector<Node<DIM, ORDER>> nodes;
        while (!q.empty()) {
            int n_next = q.size();

            uint64_t curr_child_idx = parent_idx + q.size();

            for (int i = 0; i < n_next; ++i) {
                auto [box, node_key] = q.front();
                q.pop();
                Node<DIM, ORDER> test_node(box);
                nodes.push_back(Node<DIM, ORDER>(box));

                auto &node = nodes.back();
                if (!node.fit(f, tol_)) {
                    node.first_child_idx = curr_child_idx;
                    curr_child_idx += NChild;

                    key_t parent = node_key << DIM;
                    VEC &center = box.center;
                    for (key_t child = 0; child < NChild; ++child) {
                        VEC offset;

                        // Extract sign of each offset component from the bits of child
                        // Basically: permute all possible offsets
                        for (int j = 0; j < DIM; ++j) {
                            double signed_hw[2] = {-half_width[j], half_width[j]};
                            offset[j] = signed_hw[(child >> j) & 1];
                        }

                        q.push(std::make_pair(DBox(center + offset, half_width), parent + child));
                    }
                }
            }

            if (!q.empty())
                max_depth_++;

            half_width *= 0.5;
            if ((1 << (DIM * max_depth_)) == q.size())
                n_subtrees_ *= 2;
            else {
                for (int j = 0; j < DIM; ++j)
                    bin_size_[j] = 2.0 * box_.half_length[j] / n_subtrees_[j];
                lower_left_ = box_.center - box_.half_length;

                subtrees_.reserve(n_subtrees_.prod());
                for (int i_bin = 0; i_bin < n_subtrees_.prod(); ++i_bin) {
                    Eigen::Vector<int, DIM> bins = get_bins(i_bin);

                    VEC parent_center =
                        (bins.template cast<double>().array() + 0.5) * bin_size_.array() + lower_left_.array();

                    Box<DIM> root_box = {parent_center, 0.5 * bin_size_};
                    subtrees_.push_back(FunctionTree<DIM, ORDER>(root_box, f_, tol_));
                }
            }
        }
    }

    Eigen::Vector<int, DIM> get_bins(const int i_bin) {
        if constexpr (DIM == 1)
            return Eigen::Vector<int, DIM>{i_bin};
        else if constexpr (DIM == 2)
            return Eigen::Vector<int, DIM>{i_bin % n_subtrees_[0], i_bin / n_subtrees_[0]};
        else if constexpr (DIM == 3)
            return Eigen::Vector<int, DIM>{i_bin % n_subtrees_[0], (i_bin / n_subtrees_[0]) % n_subtrees_[1],
                                           i_bin / (n_subtrees_[0] * n_subtrees_[1])};
    }

    inline int get_linear_bin(const VEC &x) const {
        const VEC x_bin = x - lower_left_;
        const Eigen::Vector<int, DIM> bins = (x_bin.array() / bin_size_.array()).template cast<int>();
        if constexpr (DIM == 1)
            return bins[0];
        else if constexpr (DIM == 2)
            return bins[0] + n_subtrees_[0] * bins[1];
        else if constexpr (DIM == 3)
            return bins[0] + n_subtrees_[0] * bins[1] + n_subtrees_[0] * n_subtrees_[1] * bins[2];
    }

    const Node<DIM, ORDER> &find_node(const VEC &x) const { return subtrees_[get_linear_bin(x)].find_node_traverse(x); }

    const double eval(const VEC &x) const { return find_node(x).eval(x); }

    double operator()(const VEC &x) const { return eval(x); }
};

template <int D, int ORDER>
const Eigen::PartialPivLU<typename Node<D, ORDER>::VanderMat>
    Node<D, ORDER>::VLU_ = Eigen::PartialPivLU<typename Node<D, ORDER>::VanderMat>(Node<D, ORDER>::calc_vandermonde());
} // namespace baobzi

#endif
