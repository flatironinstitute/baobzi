#ifndef BAOBZI_TEMPLATE_HPP
#define BAOBZI_TEMPLATE_HPP
#include <functional>
#define _USE_MATH_DEFINES

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
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

namespace msgpack {
MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS) {
    namespace adaptor {

    // Place class template specialization here
    template <>
    struct convert<baobzi_header_t> {
        msgpack::object const &operator()(msgpack::object const &o, baobzi_header_t &v) const {
            if (o.type != msgpack::type::ARRAY)
                throw msgpack::type_error();
            if (o.via.array.size != 3)
                throw msgpack::type_error();
            v = baobzi_header_t{.dim = o.via.array.ptr[0].as<int>(),
                                .order = o.via.array.ptr[1].as<int>(),
                                .version = o.via.array.ptr[2].as<int>()};
            return o;
        }
    };

    template <>
    struct pack<baobzi_header_t> {
        template <typename Stream>
        packer<Stream> &operator()(msgpack::packer<Stream> &o, baobzi_header_t const &v) const {
            o.pack_array(3);
            o.pack(v.dim);
            o.pack(v.order);
            o.pack(v.version);
            return o;
        }
    };

    } // namespace adaptor
} // MSGPACK_API_VERSION_NAMESPACE(MSGPACK_DEFAULT_API_NS)
} // namespace msgpack

/// Namespace for baobzi
namespace baobzi {
static baobzi_input_t input_default = {.func = NULL,
                                       .data = NULL,
                                       .dim = 1,
                                       .output_dim = 1,
                                       .order = 8,
                                       .tol = 1E-10,
                                       .minimum_leaf_fraction = 0.0,
                                       .split_multi_eval = 0,
                                       .min_depth = 0,
                                       .max_depth = 50};

using raw_leaf_node = struct {
    double a;
    double L;
    const double *coeffs;
};

struct leaf_compare {
    bool operator()(baobzi::raw_leaf_node a, baobzi::raw_leaf_node b) { return a.a < b.a; };
    bool operator()(baobzi::raw_leaf_node a, double b) { return a.a < b; };
};

using index_t = uint32_t; ///< Type specifying indexing into flattened tree

class MaxDepthExceeded : public std::exception {
    virtual const char *what() const throw() { return "Baobzi fit error: tree depth exceeded max allowed input depth"; }
};

template <int DIM, int ORDER, int ISET, typename T>
class Function;

/// @brief Structure to represent geometric portion of Baobzi nodes
/// @tparam DIM number of dimensions of box
/// @tparam ISET Instruction set index (dummy variable to force alignment for different instruction sets)
template <int DIM, int ISET, typename T = double>
struct Box {
    using VecDimD = Eigen::Vector<T, DIM>; ///< DIM dimensional vector type

    VecDimD center;          ///< Center of box
    VecDimD inv_half_length; ///< 1.0 / half the dimension of the box

    Box<DIM, ISET, T>() = default; ///< default constructor for msgpack happiness
    /// @brief Constructor, just copies x, hl over
    Box<DIM, ISET, T>(const VecDimD &x, const VecDimD &hl)
        : center(x), inv_half_length(VecDimD::Ones().array() / hl.array()) {}

    /// @brief return vector of box half lengths along each dimension
    inline VecDimD half_length() const { return VecDimD::Ones().array() / inv_half_length.array(); }

    /// @brief MSGPACK serialization magic
    MSGPACK_DEFINE(center, inv_half_length);
};

/// @brief Return an estimate of the error for a given set of coefficients
/// @param[in] coeffs one or two dimensional Vector/Matrix of coefficients
/// @returns estimation of error given those coefficients
template <typename T>
inline T standard_error(const Eigen::Ref<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> &coeffs) {
    T maxcoeff = 0.0;
    T scaling_factor = 1.0;
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
/// @param[in] x position of point to evaluate (pre-normalized on interval from -1:1)
/// @param[in] coeffs_raw flat vector of coefficients
/// @returns value of interpolating function at x
template <int DIM, int ORDER, int ISET, typename T = double>
inline T cheb_eval(const Eigen::Vector<T, DIM> &x, const double *coeffs_raw);

template <int ORDER, int ISET, typename T = double>
inline T cheb_eval(const Eigen::Vector<T, 1> &x, const T *c) {
    // note (RB): uses clenshaw's method to avoid direct calculation of recurrence relation of
    // T_i, where res = \Sum_i T_i c_i
    const T x2 = 2 * x[0];

    T c0 = c[0];
    T c1 = c[1];
    for (int i = 2; i < ORDER; ++i) {
        T tmp = c1;
        c1 = c[i] - c0;
        c0 = tmp + c0 * x2;
    }

    return c1 + c0 * x[0];
}

template <int ORDER, int ISET, typename T = double>
inline T cheb_eval(const Eigen::Vector<T, 2> &x, const T *coeffs_raw) {
    // note (RB): There is code to do this with clenshaw's method (twice), but it doesn't seem
    // faster (isolated tests shows it's 3x faster, but that doesn't bear fruit in production
    // and this is, imho, clearer)
    Eigen::Matrix<T, 2, ORDER> Tns;
    Tns.col(0).setOnes();
    Tns.col(1) = x;
    for (int i = 2; i < ORDER; ++i)
        Tns.col(i) = 2 * x.array() * Tns.col(i - 1).array() - Tns.col(i - 2).array();

    Eigen::Map<const Eigen::Matrix<T, ORDER, ORDER>> coeffs(coeffs_raw);

    return Tns.row(0).transpose().dot(coeffs * Tns.row(1).transpose());
}

template <int ORDER, int ISET, typename T = double>
inline T cheb_eval(const Eigen::Vector<T, 3> &x, const T *coeffs_raw) {
    Eigen::Vector<T, ORDER> Tn[3];
    Tn[0][0] = Tn[1][0] = Tn[2][0] = 1.0;
    for (int i = 0; i < 3; ++i) {
        Tn[i][1] = x[i];
        for (int j = 2; j < ORDER; ++j)
            Tn[i][j] = 2 * x[i] * Tn[i][j - 1] - Tn[i][j - 2];
    }

    T res = 0.0;
    using map_t = Eigen::Map<const Eigen::Matrix<T, ORDER, ORDER>>;
    for (int i = 0; i < ORDER; ++i)
        res += Tn[0][i] * Tn[1].dot(map_t(coeffs_raw + i * ORDER * ORDER) * Tn[2]);
    return res;
}

/// @brief Node in baobzi::FunctionTree. If leaf, contains evaluation data, otherwise children
/// @tparam DIM dimension of function
/// @tparam ORDER order of evaluation polynomial
/// @tparam ISET instruction set index (dummy variable to force alignment for different instruction sets)
template <int DIM, int ORDER, int ISET, typename T = double>
class Node {
  public:
    using VecDimD = Eigen::Vector<T, DIM>;      ///< D dimensional vector type
    using VecOrderD = Eigen::Vector<T, ORDER>;  ///< ORDER dimensional vector type
    using Func = Function<DIM, ORDER, ISET, T>; ///< Type of boabzi function this belongs to

    Box<DIM, ISET, T> box_;                                       ///< Geometric position/size of this node
    uint64_t coeff_offset = std::numeric_limits<uint64_t>::max(); ///< Flattened chebyshev coeffs
    uint32_t first_child_idx = -1; ///< First child's index in a flattened list of all nodes
    uint32_t output_dim = 1;

    Node<DIM, ORDER, ISET, T>() = default; ///< Default constructor for msgpack happiness

    /// @brief Construct node from box (without fitting)
    /// @param [in] box box this node represents
    Node<DIM, ORDER, ISET, T>(const Box<DIM, ISET, T> &box) : box_(box) {}

    /// @brief check if node is leaf
    /// @return true if leaf, false otherwise
    inline bool is_leaf() const { return coeff_offset != std::numeric_limits<uint64_t>::max(); }

    /// @brief Fit node to a given tolerance. If fit succeeds, set leaf and coeffs, otherwise ... don't
    ///
    /// @param[in] input parameters for fit (function, tol, etc)
    /// @returns coefficient vector list if fit successful, empty list if not good enough
    std::vector<T> fit(const baobzi_input_t *input, std::function<void(const double *, double *, const void *)> func,
                       const std::vector<T> &samples) {
        output_dim = input->output_dim;

        const VecDimD half_length = box_.half_length();
        const auto lb = (box_.center - half_length);
        const auto ub = (box_.center + half_length);

        if constexpr (DIM == 1) {
            Eigen::Vector<double, ORDER> xvec = Func::get_cheb_nodes(lb[0], ub[0]).template cast<double>();

            std::vector<T> coeffs_stl(ORDER * output_dim);
            Eigen::MatrixXd actual_vals(output_dim, ORDER);
            for (int i = 0; i < ORDER; ++i)
                func(&xvec[i], actual_vals.data() + i * output_dim, input->data);

            for (int i_dim = 0; i_dim < output_dim; ++i_dim) {
                Eigen::Vector<T, ORDER> F = actual_vals.row(i_dim);
                Eigen::Vector<T, ORDER> coeffs = Func::VLU_.solve(F);

                if (standard_error<T>(coeffs) > input->tol)
                    return std::vector<T>();

                for (int i = 0; i < coeffs.size(); ++i)
                    coeffs_stl[i + ORDER * i_dim] = coeffs(ORDER - i - 1);
            }

            for (const auto &sample : samples) {
                if (sample < lb[0] || sample >= ub[0])
                    continue;

                const VecDimD x(sample);
                const VecDimD xinterp = (x - box_.center).array() * box_.inv_half_length.array();

                T actual_val[output_dim];
                func(x.data(), actual_val, input->data);

                for (int i = 0; i < output_dim; ++i) {
                    if (actual_val[i] < 1E-200)
                        continue;

                    T test_val = cheb_eval<ORDER, ISET, T>(xinterp, coeffs_stl.data() + i * ORDER);
                    T rel_error = std::fabs((actual_val[i] - test_val) / actual_val[i]);
                    if (rel_error > input->tol)
                        return std::vector<T>();
                }
            }

            coeff_offset = 0;
            return coeffs_stl;
        }
        if constexpr (DIM == 2) {
            Eigen::Matrix<T, ORDER, ORDER> F;
            VecOrderD xvec = Func::get_cheb_nodes(box_.center[0] - half_length[0], box_.center[0] + half_length[0]);
            VecOrderD yvec = Func::get_cheb_nodes(box_.center[1] - half_length[1], box_.center[1] + half_length[1]);

            for (int i = 0; i < ORDER; ++i) {
                for (int j = 0; j < ORDER; ++j) {
                    double x[2] = {xvec[i], yvec[j]};
                    func(x, &F(i, j), input->data);
                }
            }

            Eigen::Matrix<T, ORDER, ORDER> coeffs = Func::VLU_.solve(F);
            coeffs = Func::VLU_.solve(coeffs.transpose()).transpose();

            if (standard_error<T>(coeffs) > input->tol)
                return std::vector<T>();

            std::vector<T> coeffs_stl(coeffs.size());
            for (int i = 0; i < coeffs.size(); ++i)
                coeffs_stl[i] = coeffs(i);

            coeff_offset = 0;
            return coeffs_stl;
        }
        if constexpr (DIM == 3) {
            Eigen::Tensor<T, 3> F(ORDER, ORDER, ORDER);

            VecOrderD xvec = Func::get_cheb_nodes(box_.center[0] - half_length[0], box_.center[0] + half_length[0]);
            VecOrderD yvec = Func::get_cheb_nodes(box_.center[1] - half_length[1], box_.center[1] + half_length[1]);
            VecOrderD zvec = Func::get_cheb_nodes(box_.center[2] - half_length[2], box_.center[2] + half_length[2]);

            for (int i = 0; i < ORDER; ++i) {
                for (int j = 0; j < ORDER; ++j) {
                    for (int k = 0; k < ORDER; ++k) {
                        T x[3] = {xvec[i], yvec[j], zvec[k]};
                        func(x, &F(i, j, k), input->data);
                    }
                }
            }

            std::vector<T> coeffs(ORDER * ORDER * ORDER);
            Eigen::Tensor<T, 3> coeffs_tensor(ORDER, ORDER, ORDER);
            using matrix_t = Eigen::Matrix<T, ORDER, ORDER>;
            using map_t = Eigen::Map<matrix_t>;
            using tensor_t = Eigen::Tensor<T, 2>;
            for (int block = 0; block < ORDER; ++block) {
                tensor_t F_block_tensor = F.chip(block, 2);
                map_t F_block(F_block_tensor.data());

                matrix_t coeffs_tmp = Func::VLU_.solve(F_block);
                coeffs_tmp = Func::VLU_.solve(coeffs_tmp.transpose()).transpose();
                coeffs_tensor.chip(block, 2) = Eigen::TensorMap<tensor_t>(coeffs_tmp.data(), ORDER, ORDER);
            }
            for (int block = 0; block < ORDER; ++block) {
                Eigen::Tensor<T, 2> coeffs_tmp = coeffs_tensor.chip(block, 0);
                map_t coeffs_ysolve(coeffs_tmp.data());
                map_t(coeffs.data() + block * ORDER * ORDER) = Func::VLU_.solve(coeffs_ysolve.transpose()).transpose();
            }

            // Hack to use local coefficient array rather than global one
            coeff_offset = 0;
            for (int i = 0; i < ORDER; ++i) {
                for (int j = 0; j < ORDER; ++j) {
                    for (int k = 0; k < ORDER; ++k) {
                        VecDimD point = (box_.center - half_length).array() +
                                        2.0 * VecDimD{(T)i, (T)j, (T)k}.array() * half_length.array() / ORDER;

                        const T test_val = eval(point, coeffs.data());
                        T actual_val; // FIXME will break with vector-valued funcs
                        func(point.data(), &actual_val, input->data);
                        const T rel_error = std::abs((actual_val - test_val) / actual_val);

                        if (fabs(actual_val) > 1E-16 && rel_error > input->tol) {
                            coeff_offset = std::numeric_limits<uint64_t>::max();
                            return std::vector<T>();
                        }
                    }
                }
            }

            return coeffs;
        }
    }

    /// @brief eval node at point x
    /// @param[in] x point to evaluate at
    /// @param[in] coeffs flat/global coefficient array
    /// @returns function approximation at x
    inline T eval(const VecDimD &x, const T *coeffs) const {
        const VecDimD xinterp = (x - box_.center).array() * box_.inv_half_length.array();
        return cheb_eval<ORDER, ISET, T>(xinterp, coeffs + coeff_offset);
    }

    void eval(const VecDimD &x, T *res, const T *coeffs) const {
        const VecDimD xinterp = (x - box_.center).array() * box_.inv_half_length.array();
        for (int i = 0; i < output_dim; ++i)
            res[i] = cheb_eval<ORDER, ISET, T>(xinterp, coeffs + coeff_offset + i * ORDER);
    }

    /// @brief Calculate memory usage of self (including unused space from vector allocation)
    /// @returns size in bytes of object instance
    inline std::size_t memory_usage() const { return sizeof(*this); }

    /// @brief MSGPACK serialization magic
    MSGPACK_DEFINE(box_, first_child_idx, coeff_offset, output_dim);
};

/// @brief Represent a function in some domain as a tree of chebyshev nodes
/// @tparam DIM dimension of function
/// @tparam ORDER order of evaluation polynomial
/// @tparam ISET instruction set index (dummy variable to force alignment for different instruction sets)
template <int DIM, int ORDER, int ISET, typename T = double>
struct FunctionTree {
    static constexpr int NChild = 1 << DIM; ///< Number of children each node potentially has (2^D)
    static constexpr int Dim = DIM;         ///< Dimension of tree
    static constexpr int Order = ORDER;     ///< Order of tree

    using node_t = Node<DIM, ORDER, ISET, T>; ///< DIM,ORDER node type
    using box_t = Box<DIM, ISET, T>;          ///< DIM box type
    using VecDimD = Eigen::Vector<T, DIM>;    ///< D dimensional vector type

    std::vector<node_t> nodes_; ///< Flat list of all nodes in Tree (leaf or otherwise)
    int max_depth_;             ///< Maximum depth of tree

    /// @brief Construct tree
    /// @param[in] input parameters for fit (function, tol, etc)
    /// @param[in] coeffs flat/global coefficient vector
    /// @param[in] box box that this tree lives in
    FunctionTree<DIM, ORDER, ISET, T>(const baobzi_input_t *input, const Box<DIM, ISET, T> &box, std::vector<T> &coeffs,
                                      std::function<void(const double *, double *, const void *)> func,
                                      const std::vector<T> &samples) {
        std::queue<Box<DIM, ISET, T>> q;
        VecDimD half_width = box.half_length() * 0.5;
        q.push(box);

        index_t curr_child_idx = 1;
        max_depth_ = 0;
        while (!q.empty()) {
            int n_next = q.size();
            int node_index = nodes_.size();
            for (int i = 0; i < n_next; ++i) {
                box_t box = q.front();
                q.pop();

                nodes_.push_back(node_t(box));

                auto &node = nodes_[i + node_index];
                std::vector new_coeffs = node.fit(input, func, samples);

                if (node.is_leaf()) {
                    node.coeff_offset = coeffs.size();
                    coeffs.insert(std::end(coeffs), std::begin(new_coeffs), std::end(new_coeffs));
                } else if (!node.is_leaf()) {
                    node.first_child_idx = curr_child_idx;
                    curr_child_idx += NChild;

                    VecDimD &center = node.box_.center;
                    for (index_t child = 0; child < NChild; ++child) {
                        VecDimD offset;

                        // Extract sign of each offset component from the bits of child
                        // Basically: permute all possible offsets
                        for (int j = 0; j < DIM; ++j) {
                            T signed_hw[2] = {-half_width[j], half_width[j]};
                            offset[j] = signed_hw[(child >> j) & 1];
                        }

                        q.push(Box<DIM, ISET, T>(center + offset, half_width));
                    }
                }
            }

            if (!q.empty())
                max_depth_++;
            if (max_depth_ > input->max_depth)
                throw MaxDepthExceeded();

            half_width *= 0.5;
        }
    }

    FunctionTree<DIM, ORDER, ISET, T>() = default; ///< Default constructor for msgpack happiness

    /// @brief Find leaf node containing a point via standard pointer traversal
    /// @param[in] x point that the node will contain
    /// @return leaf node containing point x
    inline const node_t &find_node_traverse(const VecDimD &x) const {
        auto *node = &nodes_[0];
        auto *next_node = &nodes_[node->first_child_idx]; // attempt to force preload of potential next node
        while (!node->is_leaf()) {
            index_t child_idx = 0;
            for (int i = 0; i < DIM; ++i)
                child_idx = child_idx | ((x[i] > node->box_.center[i]) << i);

            node = next_node + child_idx;
            next_node = &nodes_[node->first_child_idx];
        }

        return *node;
    }

    /// @brief Get index of node at point x (relative to local nodes_ array)
    /// @param[in] x [DIM] point to lookup
    /// @returns index of node in nodes_ array containing x
    inline std::size_t get_node_index(const VecDimD &x) const {
        index_t curr_index = 0;
        while (!nodes_[curr_index].is_leaf()) {
            index_t child_idx = 0;
            for (int i = 0; i < DIM; ++i)
                child_idx = child_idx | ((x[i] > nodes_[curr_index].box_.center[i]) << i);

            curr_index = nodes_[curr_index].first_child_idx + child_idx;
        }

        return curr_index;
    }

    /// @brief Calculate total number of nodes in instance
    /// @return number of nodes in instance
    inline std::size_t size() const { return nodes_.size(); }

    /// @brief Calculate lowest depth of any node in instance (relative subtree node)
    /// @return lowest depth of all contained nodes
    inline int max_depth() const { return max_depth_; }

    /// @brief Calculate memory usage of self (including all contained nodes)
    /// @returns size in bytes of object instance
    inline std::size_t memory_usage() const {
        std::size_t memory_usage = sizeof(*this);
        for (const auto &node : nodes_)
            memory_usage += node.memory_usage();
        return memory_usage;
    }

    /// @brief eval function approximation at point
    /// @param[in] x point to evaluate function at
    /// @param[in] coeffs flat/global coefficient array
    /// @returns function approximation at point x
    inline T eval(const VecDimD &x, const T *coeffs) const { return find_node_traverse(x).eval(x, coeffs); }

    /// @brief msgpack serialization magic
    MSGPACK_DEFINE(nodes_);
};

/// @brief Represents a function in some domain as a grid of baobzi::FunctionTree objects
/// @tparam DIM dimension of function
/// @tparam ORDER order of evaluation polynomial
/// @tparam ISET instruction set index (dummy variable to force alignment for different instruction sets)
template <int DIM, int ORDER, int ISET = 0, typename T = double>
class Function {
  public:
    using VecDimD = Eigen::Vector<T, DIM>;            ///< DIM dimensional vector type
    using VecOrderD = Eigen::Vector<T, ORDER>;        ///< Order dimensional vector type
    using VanderMat = Eigen::Matrix<T, ORDER, ORDER>; ///< VanderMonde Matrix type
    using node_t = Node<DIM, ORDER, ISET, T>;         ///< DIM,ORDER Node type (duh)
    using box_t = Box<DIM, ISET, T>;                  ///< DIM dimensional box type

    static constexpr int NChild = 1 << DIM; ///< Number of children each node potentially has (2^D)
    static constexpr int Dim = DIM;         ///< Input dimension of function
    static constexpr int Order = ORDER;     ///< Order of polynomial representation
    static constexpr int ISet = ISET;       ///< Instruction set (dummy param)
    uint32_t output_dim_ = 1;

    static std::mutex statics_mutex;            ///< mutex for locking vandermonde/chebyshev initialization
    static VecOrderD cosarray_;                 ///< Cached array of cosine values at chebyshev nodes
    static Eigen::PartialPivLU<VanderMat> VLU_; ///< Cached LU decomposition of Vandermonde matrix

    baobzi_input_t input_;
    box_t box_;           ///< box representing the domain of our function
    T tol_;               ///< Desired relative tolerance of our approximation
    VecDimD lower_left_;  ///< Bottom 'corner' of our domain
    VecDimD upper_right_; ///< Upper 'corner' of our domain

    std::vector<FunctionTree<DIM, ORDER, ISET, T>> subtrees_; ///< Grid of FunctionTree objects that do the work
    Eigen::Vector<int, DIM> n_subtrees_;    ///< Number of subtrees in each linear dimension of our space
    std::vector<int> subtree_node_offsets_; ///< n_subtrees array of offsets for where in the global array of node
                                            ///< pointers the global node pointer array starts
    std::vector<node_t *> node_pointers_;   ///< Vector of pointers to every node from every subtree
    VecDimD inv_bin_size_;                  ///< Inverse linear dimensions of the bins that our subtrees live

    std::vector<T> coeffs_; ///< Flat vector of all chebyshev coefficients from all leaf nodes

    bool split_multi_eval_ = true; ///< Split node-search and evaluation when evaluating multiple points

    /// Structure containing info about self creation :D
    struct {
        uint16_t base_depth = 0;   ///< depth of subtrees
        uint64_t n_evals_root = 0; ///< number of function evals before subtree calls
        uint32_t t_elapsed = 0;    ///< time in milliseconds to create object
    } stats_;

    /// @brief Calculate memory_usage of this object in bytes
    /// @returns Memory usage of baobzi object in bytes
    std::size_t memory_usage() const {
        std::size_t mem = sizeof(*this);
        mem += subtree_node_offsets_.capacity() * sizeof(subtree_node_offsets_[0]);
        mem += node_pointers_.capacity() * sizeof(node_pointers_[0]);
        mem += coeffs_.capacity() * sizeof(T);
        for (const auto &subtree : subtrees_)
            mem += subtree.memory_usage();
        return mem;
    }

    /// @brief Calculate and print various information about object instance to stdout
    void print_stats() const {
        std::size_t n_nodes = 0;
        std::size_t n_leaves = 0;
        std::size_t n_subtrees = subtrees_.size();
        int max_depth = 0;
        std::size_t mem = memory_usage();
        for (const auto &subtree : subtrees_) {
            n_nodes += subtree.size();
            max_depth = std::max(max_depth, subtree.max_depth());
            for (const auto &node : subtree.nodes_)
                n_leaves += node.is_leaf();
        }

        std::cout << "Baobzi function mapping " << DIM << " to " << output_dim_ << std::endl;
        std::cout << "Tree represented by " << n_nodes << " nodes, of which " << n_leaves << " are leaves\n";
        std::cout << "Nodes are distributed across " << n_subtrees << " subtrees at an initial depth of "
                  << stats_.base_depth << " with a maximum subtree depth of " << max_depth << "\n";
        std::cout << "Total function evaluations required for fit: "
                  << n_nodes * (int)std::pow(ORDER, DIM) + stats_.n_evals_root << std::endl;
        std::cout << "Total time to create tree: " << stats_.t_elapsed << " milliseconds\n";
        std::cout << "Approximate memory usage of tree: " << (T)mem / (1024 * 1024) << " MiB" << std::endl;
    }

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
                V(i, j) = T(2) * V(i - 1, j) * cosarray_(j) - V(i - 2, j);
            }
        }

        return V.transpose();
    }

    /// @brief calculate chebyshev nodes on bounds [lb, ub]
    /// @param[in] lb lower bound
    /// @param[in] ub upper bound
    /// @returns vector of chebyshev nodes scaled within [lb, ub]
    static inline VecOrderD get_cheb_nodes(T lb, T ub) { return 0.5 * ((lb + ub) + (ub - lb) * cosarray_.array()); }

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
    /// @param[in] samples list of points to force fit check
    Function<DIM, ORDER, ISET, T>(const baobzi_input_t *input, const T *xp, const T *lp,
                                  std::function<void(const double *, double *, const void *)> func = nullptr,
                                  const std::vector<T> &samples = {})
        : box_(VecDimD(xp), VecDimD(lp)), tol_(input->tol), split_multi_eval_(input->split_multi_eval),
          output_dim_(input->output_dim), input_(*input) {
        auto t_start = std::chrono::steady_clock::now();
        init_statics();

        if (!func)
            func = input->func;

        VecDimD l(lp);
        VecDimD x(xp);
        std::queue<box_t> q;
        std::queue<box_t> maybe_q;

        for (int i = 0; i < DIM; ++i)
            n_subtrees_[i] = l[i] / l.minCoeff();

        q.push(box_t(x, l));

        // Half-width of next children
        VecDimD half_width = l * 0.5;

        // Breadth first search. Step through each level of the tree and test fit all of the nodes
        // We exit when a level isn't completely filled with parent nodes (rather than leaves)
        // This way we can always avoid redundant traversals by jumping straight to a root node of a subtree
        while (!q.empty()) {
            int n_next = q.size();

            auto add_node_children_to_queue = [](std::queue<box_t> &theq, const VecDimD &center,
                                                 const VecDimD &half_width) {
                for (unsigned child = 0; child < NChild; ++child) {
                    VecDimD offset;

                    // Extract sign of each offset component from the bits of child
                    // Basically: permute all possible offsets
                    for (int j = 0; j < DIM; ++j) {
                        T signed_hw[2] = {-half_width[j], half_width[j]};
                        offset[j] = signed_hw[(child >> j) & 1];
                    }

                    theq.push(box_t(center + offset, half_width));
                }
            };

            std::vector<node_t> nodes;
            T leaf_fraction = 0.0;
            for (int i = 0; i < n_next; ++i) {
                box_t box = q.front();
                q.pop();

                nodes.emplace_back(node_t(box));
                auto &node = nodes.back();
                node.fit(input, func, samples);

                if (!node.is_leaf() || stats_.base_depth < input->min_depth) {
                    add_node_children_to_queue(q, node.box_.center, half_width);
                } else {
                    leaf_fraction += 1.0;
                    add_node_children_to_queue(maybe_q, node.box_.center, half_width);
                }
            }
            stats_.n_evals_root += nodes.size() * std::pow(ORDER, DIM);

            leaf_fraction /= nodes.size();
            if (leaf_fraction < input->minimum_leaf_fraction) {
                while (!maybe_q.empty()) {
                    box_t box = maybe_q.front();
                    maybe_q.pop();
                    q.push(box);
                }
            }

            half_width *= 0.5;
            if ((1 << (DIM * (stats_.base_depth + 1))) == q.size()) {
                n_subtrees_ *= 2;
                stats_.base_depth++;
                if (stats_.base_depth > input->max_depth)
                    throw MaxDepthExceeded();
            } else
                break;
        }

        VecDimD bin_size;
        VecDimD half_length = box_.half_length();
        for (int j = 0; j < DIM; ++j) {
            bin_size[j] = 2.0 * half_length[j] / n_subtrees_[j];
            inv_bin_size_[j] = 0.5 * n_subtrees_[j] / half_length[j];
        }
        lower_left_ = box_.center - half_length;
        upper_right_ = box_.center + half_length;

        subtrees_.reserve(n_subtrees_.prod());

        auto input_local = *input;
        input_local.max_depth -= stats_.base_depth;
        for (int i_bin = 0; i_bin < n_subtrees_.prod(); ++i_bin) {
            Eigen::Vector<int, DIM> bins = get_bins(i_bin);

            VecDimD parent_center = (bins.template cast<T>().array() + 0.5) * bin_size.array() + lower_left_.array();

            Box<DIM, ISET, T> root_box = {parent_center, 0.5 * bin_size};
            subtrees_.push_back(FunctionTree<DIM, ORDER, ISET, T>(&input_local, root_box, coeffs_, func, samples));
        }

        auto t_end = std::chrono::steady_clock::now();
        auto t_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
        stats_.t_elapsed = t_elapsed.count();
        build_cache();
    }

    /// @brief Build any intermediate state necessary for computation
    void build_cache() {
        subtree_node_offsets_.resize(n_subtrees_.prod());
        subtree_node_offsets_[0] = 0;
        for (int i = 1; i < subtree_node_offsets_.size(); ++i)
            subtree_node_offsets_[i] = subtree_node_offsets_[i - 1] + subtrees_[i - 1].size();

        auto n_nodes_tot = std::accumulate(subtrees_.begin(), subtrees_.end(), (std::size_t)0,
                                           [](size_t prior, auto &subtree) { return prior + subtree.size(); });

        node_pointers_.resize(n_nodes_tot);

        int i = 0;
        for (auto &subtree : subtrees_)
            for (node_t &node : subtree.nodes_)
                node_pointers_[i++] = &node;
    }

    /// @brief default constructor for msgpack magic
    Function<DIM, ORDER, ISET, T>() { init_statics(); };

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
    inline int get_linear_bin(const Eigen::Vector<T, 1> &x) const {
        const T x_bin = x[0] - lower_left_[0];
        return x_bin * inv_bin_size_[0];
    }

    /// @brief find linear index of bin at a point
    /// @param[in] x [2] position to find bin
    /// @returns linear index of bin that x lives in
    inline int get_linear_bin(const Eigen::Vector<T, 2> &x) const {
        const VecDimD x_bin = x - lower_left_;
        const Eigen::Vector<int, DIM> bin = (x_bin.array() * inv_bin_size_.array()).template cast<int>();
        return bin[0] + n_subtrees_[0] * bin[1];
    }

    /// @brief find linear index of bin at a point
    /// @param[in] x [3] position to find bin
    /// @returns linear index of bin that x lives in
    inline int get_linear_bin(const Eigen::Vector<T, 3> &x) const {
        const VecDimD x_bin = x - lower_left_;
        const Eigen::Vector<int, DIM> bin = (x_bin.array() * inv_bin_size_.array()).template cast<int>();
        return bin[0] + n_subtrees_[0] * bin[1] + n_subtrees_[0] * n_subtrees_[1] * bin[2];
    }

    /// @brief get constant reference to leaf node that contains a point
    /// @param[in] x point of interest
    /// @returns constant reference to leaf node that contains x
    inline const node_t &find_node(const VecDimD &x) const {
        return subtrees_[get_linear_bin(x)].find_node_traverse(x);
    }

    /// @brief eval function approximation at point
    /// @param[in] x point to evaluate function at
    /// @returns function approximation at point x
    inline T eval(const VecDimD &x) const {
        return (x.array() < lower_left_.array()).any() || (x.array() >= upper_right_.array()).any()
                   ? NAN
                   : find_node(x).eval(x, coeffs_.data());
    }

    /// @brief eval function approximation at point
    /// @param[in] xp [DIM] point to evaluate function at
    /// @returns function approximation at point xp
    inline T eval(const T *xp) const { return eval(VecDimD(xp)); }

    /// @brief get index of node (across all subnodes)
    /// @param[in] x [DIM] point to find the node of
    /// @returns index in global node array
    inline std::size_t get_global_node_index(const VecDimD &x) const {
        int i_sub = get_linear_bin(x);
        return subtree_node_offsets_[i_sub] + subtrees_[i_sub].get_node_index(x);
    }

    inline void eval(const VecDimD &x, T *res) const {
        if ((x.array() < lower_left_.array()).any() || (x.array() >= upper_right_.array()).any()) {
            for (int i = 0; i < output_dim_; ++i)
                res[i] = NAN;
        }

        find_node(x).eval(x, res, coeffs_.data());
    }

    inline void eval(const T *xp, T *res) const { eval(VecDimD(xp), res); }

    /// @brief eval function approximation at n_trg points
    /// @param[in] xp [DIM * n_trg] array of points to evaluate function at
    /// @param[out] res [n_trg] array of results
    /// @param[in] n_trg number of points to evaluate
    inline void eval(const T *xp, T *res, int n_trg) const {
        if (split_multi_eval_) {
            std::vector<std::pair<node_t *, VecDimD>> node_map(n_trg);
            for (int i = 0; i < n_trg; ++i) {
                VecDimD xi = VecDimD(xp + DIM * i);
                node_t *node_ptr =
                    (xi.array() < lower_left_.array()).any() || (xi.array() >= upper_right_.array()).any()
                        ? nullptr
                        : node_pointers_[get_global_node_index(xi)];

                node_map[i] = std::make_pair(node_ptr, xi);
            }

            for (int i_trg = 0; i_trg < n_trg; i_trg++) {
                res[i_trg] = node_map[i_trg].first == nullptr
                                 ? NAN
                                 : node_map[i_trg].first->eval(node_map[i_trg].second, coeffs_.data());
            }
        } else {
            for (int i_trg = 0; i_trg < n_trg; i_trg++) {
                res[i_trg] = eval(VecDimD(xp + DIM * i_trg));
            }
        }
    }

    /// @brief eval function approximation at point
    /// @param[in] x [DIM] point to evaluate function at
    /// @returns function approximation at point x
    inline T operator()(const VecDimD &x) const { return output_dim_ > 1 ? NAN : eval(x); }

    /// @brief eval function approximation at point
    /// @param[in] x point to evaluate function at
    /// @returns function approximation at point x
    inline T operator()(const T *x) const { return output_dim_ > 1 ? NAN : eval(x); }

    /// @brief eval function approximation at n_trg points
    /// @param[in] xp [DIM * n_trg] array of points to evaluate function at
    /// @param[out] res [DIM * n_trg] array of results
    /// @param[in] n_trg number of points to evaluate
    inline void operator()(const T *xp, T *res, int n_trg) const { eval(xp, res, n_trg); }

    /// @brief eval function approximation at 1 point by reference
    /// @param[in] xp [DIM] array of points to evaluate function at
    /// @param[out] res [DIM] array of results
    inline void operator()(const T *xp, T *res) const { eval(xp, res); }

    /// @brief save function approximation to file
    /// @param[in] filename path to save file at
    void save(const char *filename) const {
        std::ofstream ofs(filename, std::ofstream::binary | std::ofstream::out);
        baobzi_header_t params{Dim, Order, BAOBZI_HEADER_VERSION};
        msgpack::pack(ofs, params);
        msgpack::pack(ofs, *this);
    }

    std::vector<raw_leaf_node> get_leaves() const {
        std::vector<raw_leaf_node> leaves;

        for (const auto &subtree : subtrees_) {
            for (const auto &node : subtree.nodes_) {
                if (!node.is_leaf())
                    continue;

                double L = 2.0 * node.box_.half_length()[0];
                double a = node.box_.center[0] - 0.5 * L;
                const double *coeffs = node.coeff_offset + coeffs_.data();
                leaves.emplace_back(raw_leaf_node{a, L, coeffs});
            }
        }

        std::sort(leaves.begin(), leaves.end(), leaf_compare());
        return leaves;
    }

    std::pair<VecDimD, VecDimD> get_bounds() const { return std::make_pair(lower_left_, upper_right_); }

    Function<DIM, ORDER, ISET> shallow_copy() const {
        Function<DIM, ORDER, ISET> other;
        other.n_subtrees_ = n_subtrees_;
        other.lower_left_ = lower_left_;
        other.inv_bin_size_ = inv_bin_size_;
        other.box_ = box_;
        other.inv_bin_size_ = inv_bin_size_;
        other.split_multi_eval_ = split_multi_eval_;
        other.input_ = input_;
        other.input_.func = nullptr;
        other.input_.data = nullptr;

        return other;
    }

    std::pair<baobzi_input_t, box_t> operator_helper(const Function<DIM, ORDER, ISET> &B) const {
        const auto &A = *this;
        double lbound =
            std::max(A.box_.center[0] - A.box_.half_length()[0], B.box_.center[0] - B.box_.half_length()[0]);
        double rbound =
            std::min(A.box_.center[0] + A.box_.half_length()[0], B.box_.center[0] + B.box_.half_length()[0]);
        box_t new_box{VecDimD(0.5 * (rbound + lbound)), VecDimD(0.5 * (rbound - lbound))};
        baobzi_input_t new_input = input_;
        new_input.func = nullptr;
        new_input.data = nullptr;
        return {new_input, new_box};
    }

    Function<DIM, ORDER, ISET> operator+(const Function<DIM, ORDER, ISET> &B) const {
        static_assert(DIM == 1, "Baobzi: Function addition only defined for 1D functions");
        const auto &A = *this;
        auto [new_input, new_box] = operator_helper(B);
        const int output_dim = input_.output_dim;
        const auto func = [&A, &B, output_dim](const double *x, double *y, const void *data) {
            T res2[output_dim];
            A(x, y);
            B(x, res2);
            for (int i = 0; i < output_dim; ++i)
                y[i] += res2[i];
        };

        return Function<DIM, ORDER, ISET>(&new_input, new_box.center.data(), new_box.half_length().data(), func, {});
    }

    Function<DIM, ORDER, ISET> operator*(const Function<DIM, ORDER, ISET> &B) const {
        static_assert(DIM == 1, "Baobzi: Function multiplication only defined for 1D functions");
        const auto &A = *this;
        auto [new_input, new_box] = operator_helper(B);
        const int output_dim = input_.output_dim;
        const auto func = [&A, &B, output_dim](const double *x, double *y, const void *data) {
            T res2[output_dim];
            A(x, y);
            B(x, res2);
            for (int i = 0; i < output_dim; ++i)
                y[i] *= res2[i];
        };

        return Function<DIM, ORDER, ISET>(&new_input, new_box.center.data(), new_box.half_length().data(), func, {});
    }

    Function<DIM, ORDER, ISET> operator/(const Function<DIM, ORDER, ISET> &B) const {
        static_assert(DIM == 1, "Baobzi: Function division only defined for 1D functions");
        const auto &A = *this;
        auto [new_input, new_box] = operator_helper(B);
        const int output_dim = input_.output_dim;
        const auto func = [&A, &B, output_dim](const double *x, double *y, const void *data) {
            T res2[output_dim];
            A(x, y);
            B(x, res2);
            for (int i = 0; i < output_dim; ++i)
                y[i] /= res2[i];
        };

        return Function<DIM, ORDER, ISET>(&new_input, new_box.center.data(), new_box.half_length().data(), func, {});
    }

    Function<DIM, ORDER, ISET> operator-(const Function<DIM, ORDER, ISET> &B) const {
        static_assert(DIM == 1, "Baobzi: Function subtraction only defined for 1D functions");
        const auto &A = *this;
        auto [new_input, new_box] = operator_helper(B);
        const int output_dim = input_.output_dim;
        const auto func = [&A, &B, output_dim](const double *x, double *y, const void *data) {
            T res2[output_dim];
            A(x, y);
            B(x, res2);
            for (int i = 0; i < output_dim; ++i)
                y[i] -= res2[i];
        };

        return Function<DIM, ORDER, ISET>(&new_input, new_box.center.data(), new_box.half_length().data(), func, {});
    }

    template <typename U>
    Function<DIM, ORDER, ISET> operator*(const U &scale_factor) const {
        Function<DIM, ORDER, ISET> copy = *this;

        Eigen::Map<Eigen::VectorXd> coeffs(copy.coeffs_.data(), copy.coeffs_.size());
        coeffs *= scale_factor;

        return copy;
    }

    template <typename U>
    friend Function<DIM, ORDER, ISET> operator*(const U &scale_factor, const Function<DIM, ORDER, ISET> &func) {
        return func * scale_factor;
    }

    template <typename U>
    Function<DIM, ORDER, ISET> operator/(const U &divisor) const {
        Function<DIM, ORDER, ISET> copy = *this;

        Eigen::Map<Eigen::VectorXd> coeffs(copy.coeffs_.data(), copy.coeffs_.size());
        coeffs /= divisor;

        return copy;
    }

    template <typename U>
    Function<DIM, ORDER, ISET> operator+(const U &shift) const {
        static_assert(DIM == 1, "Baobzi: Scalar addition only defined for 1D functions");
        Function<DIM, ORDER, ISET> copy = *this;

        for (std::size_t i = ORDER - 1; i < copy.coeffs_.size(); i += ORDER)
            copy.coeffs_[i] += shift;

        return copy;
    }

    template <typename U>
    friend Function<DIM, ORDER, ISET> operator+(const U &shift, const Function<DIM, ORDER, ISET> &func) {
        return func + shift;
    }

    template <typename U>
    Function<1, ORDER, ISET> operator-(const U &shift) const {
        return *this + (-shift);
    }

    /// @brief msgpack serialization magic
    MSGPACK_DEFINE_MAP(box_, subtrees_, n_subtrees_, tol_, lower_left_, upper_right_, inv_bin_size_, coeffs_,
                       split_multi_eval_, output_dim_);
};

template <int DIM, int ORDER, int ISET, typename T>
std::mutex Function<DIM, ORDER, ISET, T>::statics_mutex;

template <int DIM, int ORDER, int ISET, typename T>
typename Function<DIM, ORDER, ISET, T>::VecOrderD Function<DIM, ORDER, ISET, T>::cosarray_;

template <int DIM, int ORDER, int ISET, typename T>
Eigen::PartialPivLU<typename Function<DIM, ORDER, ISET, T>::VanderMat> Function<DIM, ORDER, ISET, T>::VLU_;
} // namespace baobzi

#endif
