#ifndef BAOBZI_TEMPLATE_HPP
#define BAOBZI_TEMPLATE_HPP

#include <fstream>
#include <iostream>
#include <mutex>
#include <queue>
#include <vector>

#include <msgpack.hpp>
#define EIGEN_MATRIX_PLUGIN "baobzi/eigen_matrix_plugin.h"

#define EIGEN_MAX_ALIGN_BYTES 64
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <unsupported/Eigen/CXX11/Tensor>

#include <baobzi/header.h>

/// Namespace for baobzi
namespace baobzi {

template <int DIM, int ORDER, int ISET>
class Function;

/// @brief Structure to represent geometric portion of Baobzi nodes
/// @tparam D number of dimensions of box
/// @tparam ISET Instruction set index (dummy variable to force alignment for different instruction sets)
template <int D, int ISET>
struct Box {
    using VEC = Eigen::Vector<double, D>;
    VEC center;          ///< Center of box
    VEC half_length;     ///< Half the dimension of the box
    VEC inv_half_length; ///< 1.0 / half the dimension of the box

    Box<D, ISET>() = default; ///< default constructor for msgpack happiness
    /// @brief Constructor, just copies x, hl over
    Box<D, ISET>(const VEC &x, const VEC &hl)
        : center(x), half_length(hl), inv_half_length(VEC::Ones().array() / hl.array()) {}

    /// @brief Check if point lies inside box
    /// @param[in] x point to check
    /// @returns true if point in box, false otherwise
    bool contains(const VEC &x) const {
        VEC dx = (x - center).array().abs();
        return !((dx > half_length).any());
    }

    /// @brief MSGPACK serialization magic
    MSGPACK_DEFINE(center, half_length, inv_half_length);
};

/// @brief Return an estimate of the error for a given set of coefficientsj
/// @param[in] coeffs one or two dimensional Vector/Matrix of coefficients
/// @returns estimation of error given those coefficients
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

/// @brief Evaluate chebyshev polynomial given a box and a point inside that box
/// @tparam DIM dim of chebyshev polynomial to evaluate
/// @tparam ORDER order of chebyshev polynomial to evaluate
/// @tparam ISET Instruction set index (dummy variable to force alignment for different instruction sets)
/// @param[in] x position of point to evaluate
/// @param[in] box box that x lives in
/// @param[in] coeffs_raw flat column-major vector of coefficients
/// @returns value of interpolating function at x
template <int DIM, int ORDER, int ISET>
inline double cheb_eval(const Eigen::Vector<double, DIM> &x, const Box<DIM, ISET> &box,
                        const std::vector<double, Eigen::aligned_allocator<double>> &coeffs_raw);

template <int ORDER, int ISET>
inline double cheb_eval(const Eigen::Vector<double, 1> &x, const Box<1, ISET> &box,
                        const std::vector<double, Eigen::aligned_allocator<double>> &coeffs_raw) {
    double xd = (x[0] - box.center[0]) * box.inv_half_length[0];

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
    Eigen::Vector2d xinterp = (x - box.center).array() * box.inv_half_length.array();
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

    Eigen::Vector3d xinterp = (x - box.center).array() * box.inv_half_length.array();

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

/// @brief Node in baobzi::FunctionTree. If leaf, contains evaluation data, otherwise children
/// @tparam D dimension of function
/// @tparam ORDER order of evaluation polynomial
/// @tparam ISET instruction set index (dummy variable to force alignment for different instruction sets)
template <int D, int ORDER, int ISET>
class Node {
  public:
    using VEC = Eigen::Vector<double, D>;                          ///< D dimensional vector type
    using CoeffVec = Eigen::Vector<double, ORDER>;                 ///< ORDER dimensional vector type
    std::vector<double, Eigen::aligned_allocator<double>> coeffs_; ///< Flattened chebyshev coeffs
    using Func = Function<D, ORDER, ISET>;                         ///< Type of boabzi function this belongs to
    Box<D, ISET> box_;                                             ///< Geometric position/size of this node
    uint64_t first_child_idx = -1; ///< First child's index in a flattened list of all nodes
    bool leaf_ = false;            ///< Helper variable to determine if node is a leaf

    Node<D, ORDER, ISET>() = default; ///< Default constructor for msgpack happiness

    /// @brief Construct node from box (without fitting)
    /// @param [in] box box this node represents
    Node<D, ORDER, ISET>(const Box<D, ISET> &box) : box_(box) {}

    /// @brief check if node is leaf
    /// @return true if leaf, false otherwise
    inline bool is_leaf() const { return leaf_; }

    /// @brief Fit node to a given tolerance. If fit succeeds, set leaf and coeffs, otherwise ... don't
    ///
    /// Modifies: Node::leaf_, Node::coeffs_
    /// @param[in] input parameters for fit (function, tol, etc)
    /// @returns true if fit successful, false if not good enough
    bool fit(const baobzi_input_t *input) {
        if constexpr (D == 1) {
            Eigen::Vector<double, ORDER> F;
            CoeffVec xvec =
                Func::get_cheb_nodes(box_.center[0] - box_.half_length[0], box_.center[0] + box_.half_length[0]);

            for (int i = 0; i < ORDER; ++i)
                F(i) = input->func(&xvec[i], input->data);

            Eigen::Vector<double, ORDER> coeffs = Func::VLU_.solve(F);

            if (standard_error(coeffs) > input->tol)
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
                    F(i, j) = input->func(x, input->data);
                }
            }

            Eigen::Matrix<double, ORDER, ORDER> coeffs = Func::VLU_.solve(F);
            coeffs = Func::VLU_.solve(coeffs.transpose()).transpose();

            if (standard_error(coeffs) > input->tol)
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
                        F(i, j, k) = input->func(x, input->data);
                    }
                }
            }

            coeffs_.resize(ORDER * ORDER * ORDER);
            Eigen::Tensor<double, 3> coeffs_tensor(ORDER, ORDER, ORDER);
            using matrix_t = Eigen::Matrix<double, ORDER, ORDER>;
            using map_t = Eigen::Map<matrix_t>;
            using tensor_t = Eigen::Tensor<double, 2>;
            for (int block = 0; block < ORDER; ++block) {
                tensor_t F_block_tensor = F.chip(block, 2);
                map_t F_block(F_block_tensor.data());

                matrix_t coeffs_tmp = Func::VLU_.solve(F_block);
                coeffs_tmp = Func::VLU_.solve(coeffs_tmp.transpose()).transpose();
                coeffs_tensor.chip(block, 2) = Eigen::TensorMap<tensor_t>(coeffs_tmp.data(), ORDER, ORDER);
            }
            for (int block = 0; block < ORDER; ++block) {
                Eigen::Tensor<double, 2> coeffs_tmp = coeffs_tensor.chip(block, 0);
                map_t coeffs_ysolve(coeffs_tmp.data());
                map_t(coeffs_.data() + block * ORDER * ORDER) = Func::VLU_.solve(coeffs_ysolve.transpose()).transpose();
            }

            for (int i = 0; i < ORDER; ++i) {
                for (int j = 0; j < ORDER; ++j) {
                    for (int k = 0; k < ORDER; ++k) {
                        VEC point =
                            (box_.center - box_.half_length).array() +
                            2.0 * VEC{(double)i, (double)j, (double)k}.array() * box_.half_length.array() / ORDER;

                        const double test_val = eval(point);
                        const double actual_val = input->func(point.data(), input->data);
                        const double rel_error = std::abs((actual_val - test_val) / actual_val);

                        if (fabs(actual_val) > 1E-16 && rel_error > input->tol) {
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

    /// @brief eval node at point x
    /// @param[in] x point to evaluate at
    /// @returns function approximation at x
    inline double eval(const VEC &x) const { return cheb_eval<ORDER, ISET>(x, box_, coeffs_); }

    /// @brief MSGPACK serialization magic
    MSGPACK_DEFINE(box_, first_child_idx, leaf_, coeffs_);
};

/// @brief Represent a function in some domain as a tree of chebyshev nodes
/// @tparam DIM dimension of function
/// @tparam ORDER order of evaluation polynomial
/// @tparam ISET instruction set index (dummy variable to force alignment for different instruction sets)
template <int DIM, int ORDER, int ISET>
struct FunctionTree {
    static constexpr int NChild = 1 << DIM; ///< Number of children each node potentially has (2^D)
    static constexpr int Dim = DIM;         ///< Dimension of tree
    static constexpr int Order = ORDER;     ///< Order of tree

    using VEC = Eigen::Vector<double, DIM>;     ///< D dimensional vector type
    std::vector<Node<DIM, ORDER, ISET>> nodes_; ///< Flat list of all nodes in Tree (leaf or otherwise)

    /// @brief Construct tree
    /// @param[in] input parameters for fit (function, tol, etc)
    /// @param[in] box box that this tree lives in
    FunctionTree<DIM, ORDER, ISET>(const baobzi_input_t *input, const Box<DIM, ISET> &box) {
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

            for (size_t i = 0; i < n_next; ++i) {
                auto &node = nodes_[i + node_index];
                node.fit(input);
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

    FunctionTree<DIM, ORDER, ISET>() = default; ///< Default constructor for msgpack happiness

    /// @brief Find leaf node containing a point via standard pointer traversal
    /// @param[in] x point that the node will contain
    /// @return leaf node containing point x
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

    /// @brief eval function approximation at point
    /// @param[in] x point to evaluate function at
    /// @returns function approximation at point x
    inline double eval(const VEC &x) const { return find_node_traverse(x).eval(x); }

    /// @brief msgpack serialization magic
    MSGPACK_DEFINE(nodes_);
};

/// @brief Represents a function in some domain as a grid of baobzi::FunctionTree objects
/// @tparam DIM dimension of function
/// @tparam ORDER order of evaluation polynomial
/// @tparam ISET instruction set index (dummy variable to force alignment for different instruction sets)
template <int DIM, int ORDER, int ISET = 0>
class Function {
  public:
    static constexpr int NChild = 1 << DIM; ///< Number of children each node potentially has (2^D)
    static constexpr int Dim = DIM;         ///< Input dimension of function
    static constexpr int Order = ORDER;     ///< Order of polynomial representation
    static constexpr int ISet = ISET;       ///< Instruction set (dummy param)
    static std::mutex statics_mutex;        ///< mutex for locking vandermonde/chebyshev initialization

    using VEC = Eigen::Vector<double, DIM>;                ///< D dimensional vector type
    using CoeffVec = Eigen::Vector<double, ORDER>;         ///< Order dimensional vector type
    using VanderMat = Eigen::Matrix<double, ORDER, ORDER>; ///< VanderMonde Matrix type

    using DBox = Box<DIM, ISET>; ///< D dimensional box type

    static CoeffVec cosarray_;                  ///< Cached array of cosine values at chebyshev nodes
    static Eigen::PartialPivLU<VanderMat> VLU_; ///< Cached LU decomposition of Vandermonde matrix

    DBox box_;       ///< box representing the domain of our function
    double tol_;     ///< Desired relative tolerance of our approximation
    VEC lower_left_; ///< Bottom 'corner' of our domain

    std::vector<FunctionTree<DIM, ORDER, ISET>> subtrees_; ///< Grid of FunctionTree objects that do the work
    Eigen::Vector<int, DIM> n_subtrees_;                   ///< Number of subtrees in each linear dimension of our space
    VEC bin_size_;                                         ///< Linear dimensions of the bins that our subtrees live

    /// @brief calculate vandermonde matrix
    /// @return Vandermonde matrix for chebyshev polynomials with order=ORDER
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

    /// @brief calculate chebyshev nodes on bounds [lb, ub]
    /// @param[in] lb lower bound
    /// @param[in] ub upper bound
    /// @returns vector of chebyshev nodes scaled within [lb, ub]
    static inline CoeffVec get_cheb_nodes(double lb, double ub) {
        return 0.5 * ((lb + ub) + (ub - lb) * cosarray_.array());
    }

    /// @brief initialize static class variables
    ///
    /// Modifies baobzi::Function::cosarray_, baobzi::Function::VLU_
    static void init_statics() {
        static bool is_initialized = false;
        std::lock_guard<std::mutex> lock(statics_mutex);
        if (is_initialized)
            return;

        for (int i = 0; i < ORDER; ++i)
            cosarray_[ORDER - i - 1] = cos(M_PI * (i + 0.5) / ORDER);
        VLU_ = Eigen::PartialPivLU<VanderMat>(calc_vandermonde());
        is_initialized = true;
    }

    /// @brief Construct our Function object (fits recursively, can be slow)
    /// @param[in] input parameters for fit (function, tol, etc)
    /// @param[in] xp [dim] center of function domain
    /// @param[in] lp [dim] half length of function domain
    Function<DIM, ORDER, ISET>(const baobzi_input_t *input, const double *xp, const double *lp)
        : box_(VEC(xp), VEC(lp)), tol_(input->tol) {
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

        // Breadth first search. Step through each level of the tree and test fit all of the nodes
        // We exit when a level isn't completely filled with parent nodes (rather than leaves)
        // This way we can always avoid redundant traversals by jumping straight to a root node of a subtree
        while (!q.empty()) {
            int n_next = q.size();

            std::vector<Node<DIM, ORDER, ISET>> nodes;
            for (int i = 0; i < n_next; ++i) {
                DBox box = q.front();
                q.pop();

                nodes.emplace_back(Node<DIM, ORDER, ISET>(box));
            }

            for (int i = 0; i < nodes.size(); ++i)
                nodes[i].fit(input);

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
            subtrees_.push_back(FunctionTree<DIM, ORDER, ISET>(input, root_box));
        }
    }

    /// @brief default constructor for msgpack magic
    Function<DIM, ORDER, ISET>() { init_statics(); };

    /// @brief convert linear bin index to [dim] bin vector
    /// @param[in] i_bin linear index
    /// @returns [dim] bin vector
    inline Eigen::Vector<int, DIM> get_bins(const int i_bin) const {
        if constexpr (DIM == 1)
            return Eigen::Vector<int, DIM>{i_bin};
        else if constexpr (DIM == 2)
            return Eigen::Vector<int, DIM>{i_bin % n_subtrees_[0], i_bin / n_subtrees_[0]};
        else if constexpr (DIM == 3)
            return Eigen::Vector<int, DIM>{i_bin % n_subtrees_[0], (i_bin / n_subtrees_[0]) % n_subtrees_[1],
                                           i_bin / (n_subtrees_[0] * n_subtrees_[1])};
    }

    /// @brief find linear index of bin at a point
    /// @param[in] x [1] position to find bin
    /// @returns linear index of bin that x lives in
    inline int get_linear_bin(const Eigen::Vector<double, 1> &x) const {
        const double x_bin = x[0] - lower_left_[0];
        return x_bin / bin_size_[0];
    }

    /// @brief find linear index of bin at a point
    /// @param[in] x [2] position to find bin
    /// @returns linear index of bin that x lives in
    inline int get_linear_bin(const Eigen::Vector2d &x) const {
        const VEC x_bin = x - lower_left_;
        const Eigen::Vector<int, DIM> bin = (x_bin.array() / bin_size_.array()).template cast<int>();
        return bin[0] + n_subtrees_[0] * bin[1];
    }

    /// @brief find linear index of bin at a point
    /// @param[in] x [3] position to find bin
    /// @returns linear index of bin that x lives in
    inline int get_linear_bin(const Eigen::Vector3d &x) const {
        const VEC x_bin = x - lower_left_;
        const Eigen::Vector<int, DIM> bin = (x_bin.array() / bin_size_.array()).template cast<int>();
        return bin[0] + n_subtrees_[0] * bin[1] + n_subtrees_[0] * n_subtrees_[1] * bin[2];
    }

    /// @brief get constant reference to leaf node that contains a point
    /// @param[in] x point of interest
    /// @returns constant reference to leaf node that contains x
    inline const Node<DIM, ORDER, ISET> &find_node(const VEC &x) const {
        return subtrees_[get_linear_bin(x)].find_node_traverse(x);
    }

    /// @brief eval function approximation at point
    /// @param[in] x point to evaluate function at
    /// @returns function approximation at point x
    inline double eval(const VEC &x) const { return find_node(x).eval(x); }

    /// @brief eval function approximation at point
    /// @param[in] xp [DIM] point to evaluate function at
    /// @returns function approximation at point xp
    inline double eval(const double *xp) const { return eval(VEC(xp)); }

    /// @brief eval function approximation at ntrg points
    /// @param[in] xp [DIM * ntrg] array of points to evaluate function at
    /// @param[out] res [ntrg] array of results
    inline void eval(const double *xp, double *res, int ntrg) const {
        for (int i = 0; i < ntrg; i++)
            res[i] = eval(VEC(xp + DIM * i));
    }

    /// @brief eval function approximation at point
    /// @param[in] x [DIM] point to evaluate function at
    /// @returns function approximation at point x
    inline double operator()(const VEC &x) const { return eval(x); }

    /// @brief eval function approximation at point
    /// @param[in] x point to evaluate function at
    /// @returns function approximation at point x
    inline double operator()(const double *x) const { return eval(x); }

    /// @brief eval function approximation at ntrg points
    /// @param[in] xp [DIM * ntrg] array of points to evaluate function at
    /// @param[out] res [DIM * ntrg] array of results
    inline void operator()(const double *xp, double *res, int ntrg) const { eval(xp, res, ntrg); }

    /// @brief save function approximation to file
    /// @param[in] filename path to save file at
    void save(const char *filename) {
        std::ofstream ofs(filename, std::ofstream::binary | std::ofstream::out);
        baobzi_header_t params{Dim, Order, BAOBZI_HEADER_VERSION};
        msgpack::pack(ofs, params);
        msgpack::pack(ofs, *this);
    }

    /// @brief msgpack serialization magic
    MSGPACK_DEFINE_MAP(box_, subtrees_, n_subtrees_, tol_, lower_left_, bin_size_);
};

template <int DIM, int ORDER, int ISET>
std::mutex Function<DIM, ORDER, ISET>::statics_mutex;

template <int DIM, int ORDER, int ISET>
typename Function<DIM, ORDER, ISET>::CoeffVec Function<DIM, ORDER, ISET>::cosarray_;

template <int DIM, int ORDER, int ISET>
Eigen::PartialPivLU<typename Function<DIM, ORDER, ISET>::VanderMat> Function<DIM, ORDER, ISET>::VLU_;
} // namespace baobzi

#endif
