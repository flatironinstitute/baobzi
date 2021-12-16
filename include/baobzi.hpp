#ifndef BAOBZI_HPP
#define BAOBZI_HPP

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
Eigen::VectorXd get_chebyshev_nodes(double lb, double ub) {
    static bool first_call = true;
    static std::array<double, ORDER> cosarray;
    if (first_call) {
        for (int i = 0; i < ORDER; ++i)
            cosarray[i] = cos(M_PI * (i + 0.5) / ORDER);
        first_call = false;
    }

    Eigen::VectorXd res(ORDER);
    for (int i = 0; i < ORDER; ++i)
        res[ORDER - i - 1] = 0.5 * ((lb + ub) + (ub - lb) * cosarray[i]);

    return res;
}

template <int D>
class Box {
  public:
    using VEC = Eigen::Vector<double, D>;
    VEC center;
    double half_length;

    Box<D>() = default;
    Box<D>(const VEC &x, double hl) : center(x), half_length(hl) {}

    bool contains(const VEC &x) const {
        VEC dx = (x - center).array().abs();
        for (int i = 0; i < D; ++i)
            if (dx[i] > half_length)
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
        for (auto i = n - 2; i < n; ++i)
            maxcoeff = std::max(std::abs(coeffs(n - 1, i)), maxcoeff);
        for (auto i = n - 2; i < n; ++i)
            maxcoeff = std::max(std::abs(coeffs(i, n - 1)), maxcoeff);

        scaling_factor = std::max(scaling_factor, std::abs(coeffs(n - 1, 0)));
        scaling_factor = std::max(scaling_factor, std::abs(coeffs(0, n - 1)));
    }

    return maxcoeff / scaling_factor;
}

template <int D, int ORDER>
class LinearNode {
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

        auto x = get_chebyshev_nodes<ORDER>(-1, 1);
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

    LinearNode<D, ORDER>() = default;

    LinearNode<D, ORDER>(const Box<D> &box) : box_(box) {}

    bool is_leaf() const { return coeffs_.size(); }

    bool fit(double (*f)(VEC)) {
        if constexpr (D == 2) {
            Eigen::Matrix<double, ORDER, ORDER> F;
            CoeffVec xvec =
                get_chebyshev_nodes<ORDER>(box_.center[0] - box_.half_length, box_.center[0] + box_.half_length);
            CoeffVec yvec =
                get_chebyshev_nodes<ORDER>(box_.center[1] - box_.half_length, box_.center[1] + box_.half_length);

            for (int i = 0; i < ORDER; ++i)
                for (int j = 0; j < ORDER; ++j)
                    F(i, j) = f({xvec[i], yvec[j]});

            Eigen::Matrix<double, ORDER, ORDER> coeffs = VLU_.solve(F);
            coeffs = VLU_.solve(coeffs.transpose()).transpose();

            double tol_ = 1E-12;
            if (standard_error(coeffs) > tol_)
                return false;

            coeffs_.resize(coeffs.size());
            for (int i = 0; i < coeffs.size(); ++i)
                coeffs_[i] = coeffs(i);

            return true;
        }
        if constexpr (D == 3) {
            std::vector<double> F(ORDER * ORDER * ORDER);
            CoeffVec xvec =
                get_chebyshev_nodes<ORDER>(box_.center[0] - box_.half_length, box_.center[0] + box_.half_length);
            CoeffVec yvec =
                get_chebyshev_nodes<ORDER>(box_.center[1] - box_.half_length, box_.center[1] + box_.half_length);
            CoeffVec zvec =
                get_chebyshev_nodes<ORDER>(box_.center[2] - box_.half_length, box_.center[2] + box_.half_length);

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

            double tol_ = 1E-10;
            for (int i = 0; i < 8; ++i) {
                for (int j = 0; j < 8; ++j) {
                    for (int k = 0; k < 8; ++k) {
                        VEC point;
                        for (int l = 0; l < 3; ++l)
                            point[l] = box_.center[l] - box_.half_length + i * box_.half_length / 4.0;

                        double test_val = this->eval(point);
                        double actual_val = f(point);
                        double rel_error = std::abs(this->eval(point) - f(point));

                        if (fabs(actual_val) > 1E-16 && rel_error > tol_) {
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
        if constexpr (D == 1) {
            VEC Tns[ORDER];
            Tns[0] = VEC::Ones();
            VEC xinterp = (x - box_.center) / box_.half_length;
            Tns[1] = xinterp;
            for (int i = 2; i < ORDER; ++i) {
                Tns[i] = 2. * xinterp.array() * Tns[i - 1].array() - Tns[i - 2].array();
            }

            double res = 0.0;
            for (int i = 0; i < ORDER; ++i) {
                res += coeffs_[i] * Tns[i][0];
            }

            return res;
        }
        if constexpr (D == 2) {
            VEC xinterp = (x - box_.center) / box_.half_length;
            CoeffVec Tnx, Tny;
            Tnx[0] = Tny[0] = 1.0;
            Tnx[1] = xinterp[0];
            Tny[1] = xinterp[1];
            for (int i = 2; i < ORDER; ++i) {
                Tnx[i] = 2. * xinterp[0] * Tnx[i - 1] - Tnx[i - 2];
                Tny[i] = 2. * xinterp[1] * Tny[i - 1] - Tny[i - 2];
            }

            Eigen::Map<const Eigen::Matrix<double, ORDER, ORDER>> coeffs(coeffs_.data());

            return Tnx.dot(coeffs * Tny);
        }
        if constexpr (D == 3) {
            VEC xinterp = (x - box_.center) / box_.half_length;
            CoeffVec Tn[3];
            Tn[0][0] = Tn[1][0] = Tn[2][0] = 1.0;
            for (int i = 0; i < 3; ++i) {
                Tn[i][1] = xinterp[i];
                for (int j = 2; j < ORDER; ++j)
                    Tn[i][j] = 2. * xinterp[i] * Tn[i][j - 1] - Tn[i][j - 2];
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
class LinearTree {
  public:
    static constexpr int NChild = 1 << DIM;
    static constexpr int Dim = DIM;
    static constexpr int Order = ORDER;

    using VEC = Eigen::Vector<double, DIM>;
    using CoeffVec = Eigen::Vector<double, ORDER>;
    using DBox = Box<DIM>;

    std::vector<LinearNode<DIM, ORDER>> nodes_;
    std::vector<uint64_t> flat_map_;
    double (*f_)(VEC);
    DBox box_;
    uint16_t max_depth_;
    uint64_t parent_idx = 0;
    uint64_t child_idx_base = 0;
    uint16_t max_full_ = 0;

    LinearTree<DIM, ORDER>(const VEC &x, double l, double (*f)(VEC)) : f_(f), box_(x, l) {
        using key_t = uint64_t;
        std::queue<std::pair<DBox, key_t>> q;

        max_depth_ = 0;
        q.push(std::make_pair(DBox(x, l), 1));

        // Half-width of next children
        double half_width = l * 0.5;
        while (!q.empty()) {
            int n_next = q.size();

            uint64_t curr_child_idx = parent_idx + q.size();

            int old_size = flat_map_.size();
            flat_map_.resize(1 << ((max_depth_ + 1) * DIM));
            for (int i = old_size; i < flat_map_.size(); i++)
                flat_map_[i] = std::numeric_limits<uint64_t>::max();

            for (int i = 0; i < n_next; ++i) {
                auto [box, node_key] = q.front();
                q.pop();
                LinearNode<DIM, ORDER> test_node(box);
                flat_map_[node_key] = parent_idx++;
                // node_map_[node_key] = parent_idx++;
                nodes_.push_back(LinearNode<DIM, ORDER>(box));

                auto &node = nodes_.back();
                if (!node.fit(f)) {
                    node.first_child_idx = curr_child_idx;
                    curr_child_idx += NChild;

                    key_t parent = node_key << DIM;
                    VEC &center = box.center;
                    double signed_hw[2] = {-half_width, half_width};
                    for (key_t child = 0; child < NChild; ++child) {
                        VEC offset;

                        // Extract sign of each offset component from the bits of child
                        // Basically: permute all possible offsets
                        for (int j = 0; j < DIM; ++j)
                            offset[j] = signed_hw[(child >> j) & 1];

                        q.push(std::make_pair(DBox(center + offset, half_width), parent + child));
                    }
                }
            }

            if (!q.empty())
                max_depth_++;
            half_width *= 0.5;
            if ((1 << (DIM * max_depth_)) == q.size())
                max_full_ = max_depth_;
            std::cout << max_depth_ << " " << q.size() << " " << max_full_ << std::endl;
        }
    }

    uint64_t find_node_key(const VEC &x) const {
        VEC center = box_.center;
        double half_width = box_.half_length;

        const VEC x_scaled = 0.5 * (x - center) / box_.half_length + VEC::Ones() * 0.5;
        uint64_t m_max = calculate_key(x_scaled, max_depth_);

        uint64_t rel_depth = (max_depth_ - max_full_) * DIM;
        // Our guess is below the actual depth of the target node, move up tree until we hit it
        while (flat_map_[m_max >> (rel_depth - DIM)] != std::numeric_limits<uint64_t>::max())
            rel_depth -= DIM;

        return m_max >> rel_depth;
    }

    const LinearNode<DIM, ORDER> &find_node(const VEC &x) const { return nodes_[flat_map_[find_node_key(x)]]; }

    const LinearNode<DIM, ORDER> &find_node_traverse(const VEC &x) const {
        auto *node = &nodes_[0];
        while (!node->is_leaf()) {
            uint64_t child_idx = 0;
            for (int i = 0; i < DIM; ++i)
                child_idx = child_idx | ((x[i] > node->box_.center[i]) << i);

            node = &nodes_[node->first_child_idx + child_idx];
        }

        return *node;
    }

    const double eval(const VEC &x) const { return find_node(x).eval(x); }

    double operator()(const VEC &x) const { return eval(x); }
};
}

#endif
