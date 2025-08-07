#ifndef BAOBZI_TEMPLATE_HPP
#define BAOBZI_TEMPLATE_HPP

#define _USE_MATH_DEFINES

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <tuple>
#include <type_traits>
#include <vector>

#include <polyfit/fast_eval.hpp>

#include <baobzi.h>
#include <msgpack.hpp>

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

template <int ORDER, class Func>
class Function;

inline auto inverse_array(const auto &arr) {
    std::remove_cvref_t<decltype(arr)> res;
    for (int i = 0; i < arr.size(); ++i)
        res[i] = 1.0 / arr[i];
    return res;
}

inline auto prod(const auto &arr) {
    typename std::remove_cvref_t<decltype(arr)>::value_type res{1};
    for (const auto &el : arr)
        res *= el;
    return res;
}

inline auto scale(const auto &arr, auto factor) {
    auto res = arr;
    for (auto &el : res)
        el *= factor;
    return res;
}

template <typename T>
constexpr int get_tuple_size() {
    if constexpr (has_tuple_size_v<T>)
        return std::tuple_size_v<T>;
    else
        return 1;
};

/// @brief Structure to represent geometric portion of Baobzi nodes
/// @tparam DIM number of dimensions of box
/// @tparam ISET Instruction set index (dummy variable to force alignment for different instruction sets)
template <int DIM, int ISET, typename T = double>
struct Box {
    using VecDimD = std::array<T, DIM>; ///< DIM dimensional vector type

    VecDimD center;          ///< Center of box
    VecDimD inv_half_length; ///< 1.0 / half the dimension of the box

    Box<DIM, ISET, T>() = default; ///< default constructor for msgpack happiness
    /// @brief Constructor, just copies x, hl over
    Box<DIM, ISET, T>(const VecDimD &x, const VecDimD &hl) : center(x), inv_half_length(inverse_array(hl)) {}

    /// @brief return vector of box half lengths along each dimension
    inline VecDimD half_length() const { return inverse_array(inv_half_length); }

    /// @brief MSGPACK serialization magic
    MSGPACK_DEFINE(center, inv_half_length);
};

/// @brief Return an estimate of the error for a given set of coefficients
/// @param[in] coeffs one or two dimensional Vector/Matrix of coefficients
/// @returns estimation of error given those coefficients
inline double standard_error(const auto &polyfit, baobzi_tol_t tol_type) {
    using T = std::remove_cvref_t<decltype(polyfit)>::InputType;
    constexpr int input_dim = get_tuple_size<T>();

    T maxcoeff{0.0};
    T scaling_factor{1.0};
    const auto &coeffs = polyfit.coeffs();

    if constexpr (input_dim == 1) {
        int n = coeffs.size();
        for (auto i = 0; i < 2; ++i)
            maxcoeff = std::max(std::abs(coeffs[i]), maxcoeff);
        scaling_factor = std::max(scaling_factor, std::abs(coeffs[n - 1]));
    } else {
        throw std::runtime_error("Baobzi standard_error error: only scalar functions are currently supported");
        // int n = coeffs.size() / output_dim;
        // for (auto i = 0; i < n; ++i)
        //     maxcoeff = std::max(std::abs(coeffs(i, n - i - 1)), maxcoeff);

        // scaling_factor = std::max(scaling_factor, std::abs(coeffs(n - 1, 0)));
        // scaling_factor = std::max(scaling_factor, std::abs(coeffs(0, n - 1)));
    }

    if (tol_type == BAOBZI_TOL_RELATIVE)
        return maxcoeff / scaling_factor;
    else
        return maxcoeff;
}

/// @brief Node in baobzi::FunctionTree. If leaf, contains evaluation data, otherwise children
/// @tparam DIM dimension of function
/// @tparam ORDER order of evaluation polynomial
/// @tparam ISET instruction set index (dummy variable to force alignment for different instruction sets)
template <class Func, int ORDER, int ISET = 0>
class Node {
  public:
    using input_type = poly_eval::function_traits<Func>::arg0_type;
    using output_type = poly_eval::function_traits<Func>::result_type;
    using value_type = value_type_or_identity<input_type>::type;
    using PolyEvalType = poly_eval::FuncEval<Func, ORDER>;
    static constexpr int DIM = get_tuple_size<input_type>();

    using VecDimD = std::array<value_type, DIM>;     ///< D dimensional vector type
    using VecOrderD = std::array<value_type, ORDER>; ///< ORDER dimensional vector type

    Box<DIM, ISET, value_type> box_;                              ///< Geometric position/size of this node
    uint64_t poly_eval_id = std::numeric_limits<uint64_t>::max(); ///< Position of poly_eval object in global array
    uint32_t first_child_idx = -1; ///< First child's index in a flattened list of all nodes
    uint32_t output_dim = 1;

    Node<Func, ORDER, ISET>() = default; ///< Default constructor for msgpack happiness

    /// @brief Construct node from box (without fitting)
    /// @param [in] box box this node represents
    Node<Func, ORDER, ISET>(const Box<DIM, ISET, value_type> &box) : box_(box) {}

    /// @brief check if node is leaf
    /// @return true if leaf, false otherwise
    inline bool is_leaf() const { return poly_eval_id != std::numeric_limits<uint64_t>::max(); }

    /// @brief Fit node to a given tolerance. If fit succeeds, set leaf and coeffs, otherwise ... don't
    ///
    /// @param[in] input parameters for fit (function, tol, etc)
    /// @returns coefficient vector list if fit successful, empty list if not good enough
    std::vector<PolyEvalType> fit(const baobzi_input_t &input, const Func &func,
                                  const std::vector<value_type> &samples) {
        if (samples.size())
            throw std::runtime_error("Baobzi fit error: sample points not yet supported");

        output_dim = input.output_dim;
        if constexpr (DIM == 1) {
            const auto half_length = box_.half_length()[0];
            const auto lb = (box_.center[0] - half_length);
            const auto ub = (box_.center[0] + half_length);

            std::vector<PolyEvalType> poly_evals;
            for (int i_dim = 0; i_dim < output_dim; ++i_dim) {
                auto poly = poly_eval::make_func_eval<ORDER>(func, lb, ub);

                if (standard_error(poly, input.tol_type) > input.tol)
                    return std::vector<PolyEvalType>();
                else
                    poly_evals.push_back(poly);
            }

            poly_eval_id = 0;
            return poly_evals;
        }

        throw std::runtime_error("Baobzi fit error: only 1D functions are currently supported");
        // if constexpr (DIM == 2) {
        //     Eigen::Matrix<T, ORDER, ORDER> F;
        //     VecOrderD xvec = Func::get_cheb_nodes(box_.center[0] - half_length[0], box_.center[0] + half_length[0]);
        //     VecOrderD yvec = Func::get_cheb_nodes(box_.center[1] - half_length[1], box_.center[1] + half_length[1]);

        //     for (int i = 0; i < ORDER; ++i) {
        //         for (int j = 0; j < ORDER; ++j) {
        //             double x[2] = {xvec[i], yvec[j]};
        //             func(x, &F(i, j), input->data);
        //         }
        //     }

        //     Eigen::Matrix<T, ORDER, ORDER> coeffs = Func::VLU_.solve(F);
        //     coeffs = Func::VLU_.solve(coeffs.transpose()).transpose();

        //     if (standard_error<T>(coeffs, input->tol_type) > input->tol)
        //         return std::vector<T>();

        //     std::vector<T> coeffs_stl(coeffs.size());
        //     for (int i = 0; i < coeffs.size(); ++i)
        //         coeffs_stl[i] = coeffs(i);

        //     coeff_offset = 0;
        //     return coeffs_stl;
        // }
        // if constexpr (DIM == 3) {
        //     Eigen::Tensor<T, 3> F(ORDER, ORDER, ORDER);

        //     VecOrderD xvec = Func::get_cheb_nodes(box_.center[0] - half_length[0], box_.center[0] + half_length[0]);
        //     VecOrderD yvec = Func::get_cheb_nodes(box_.center[1] - half_length[1], box_.center[1] + half_length[1]);
        //     VecOrderD zvec = Func::get_cheb_nodes(box_.center[2] - half_length[2], box_.center[2] + half_length[2]);

        //     for (int i = 0; i < ORDER; ++i) {
        //         for (int j = 0; j < ORDER; ++j) {
        //             for (int k = 0; k < ORDER; ++k) {
        //                 T x[3] = {xvec[i], yvec[j], zvec[k]};
        //                 func(x, &F(i, j, k), input->data);
        //             }
        //         }
        //     }

        //     std::vector<T> coeffs(ORDER * ORDER * ORDER);
        //     Eigen::Tensor<T, 3> coeffs_tensor(ORDER, ORDER, ORDER);
        //     using matrix_t = Eigen::Matrix<T, ORDER, ORDER>;
        //     using map_t = Eigen::Map<matrix_t>;
        //     using tensor_t = Eigen::Tensor<T, 2>;
        //     for (int block = 0; block < ORDER; ++block) {
        //         tensor_t F_block_tensor = F.chip(block, 2);
        //         map_t F_block(F_block_tensor.data());

        //         matrix_t coeffs_tmp = Func::VLU_.solve(F_block);
        //         coeffs_tmp = Func::VLU_.solve(coeffs_tmp.transpose()).transpose();
        //         coeffs_tensor.chip(block, 2) = Eigen::TensorMap<tensor_t>(coeffs_tmp.data(), ORDER, ORDER);
        //     }
        //     for (int block = 0; block < ORDER; ++block) {
        //         Eigen::Tensor<T, 2> coeffs_tmp = coeffs_tensor.chip(block, 0);
        //         map_t coeffs_ysolve(coeffs_tmp.data());
        //         map_t(coeffs.data() + block * ORDER * ORDER) =
        //         Func::VLU_.solve(coeffs_ysolve.transpose()).transpose();
        //     }

        //     // Hack to use local coefficient array rather than global one
        //     coeff_offset = 0;
        //     for (int i = 0; i < ORDER; ++i) {
        //         for (int j = 0; j < ORDER; ++j) {
        //             for (int k = 0; k < ORDER; ++k) {
        //                 VecDimD point = (box_.center - half_length).array() +
        //                                 2.0 * VecDimD{(T)i, (T)j, (T)k}.array() * half_length.array() / ORDER;

        //                 const T test_val = eval(point, coeffs.data());
        //                 T actual_val; // FIXME will break with vector-valued funcs
        //                 func(point.data(), &actual_val, input->data);
        //                 const T rel_error = std::abs((actual_val - test_val) / actual_val);

        //                 if (fabs(actual_val) > 1E-16 && rel_error > input->tol) {
        //                     coeff_offset = std::numeric_limits<uint64_t>::max();
        //                     return std::vector<T>();
        //                 }
        //             }
        //         }
        //     }

        //     return coeffs;
        // }
    }

    // /// @brief eval node at point x
    // /// @param[in] x point to evaluate at
    // /// @param[in] coeffs flat/global coefficient array
    // /// @returns function approximation at x
    // inline output_type eval(const input_type &x, const PolyEvalType *coeffs) const {
    //     const VecDimD xinterp = (x - box_.center).array() * box_.inv_half_length.array();
    //     return cheb_eval<ORDER, ISET, T>(xinterp, coeffs + coeff_offset);
    // }

    // // void eval(const VecDimD &x, T *res, const T *coeffs) const {
    // //     const VecDimD xinterp = (x - box_.center).array() * box_.inv_half_length.array();
    // //     for (int i = 0; i < output_dim; ++i)
    // //         res[i] = cheb_eval<ORDER, ISET, T>(xinterp, coeffs + coeff_offset + i * ORDER);
    // // }

    /// @brief Calculate memory usage of self (including unused space from vector allocation)
    /// @returns size in bytes of object instance
    inline std::size_t memory_usage() const { return sizeof(*this); }

    /// @brief MSGPACK serialization magic
    MSGPACK_DEFINE(box_, first_child_idx, poly_eval_id, output_dim);
};

/// @brief Represent a function in some domain as a tree of chebyshev nodes
/// @tparam DIM dimension of function
/// @tparam ORDER order of evaluation polynomial
/// @tparam ISET instruction set index (dummy variable to force alignment for different instruction sets)
template <int ORDER, class Func>
struct FunctionTree {
    using input_type = poly_eval::function_traits<Func>::arg0_type;
    using value_type = value_type_or_identity<input_type>::type;
    using PolyEvalType = poly_eval::FuncEval<Func, ORDER>;
    static constexpr int ISET = 0; // fixme
    static constexpr int DIM = get_tuple_size<input_type>();
    static constexpr int NChild = 1 << DIM; ///< Number of children each node potentially has (2^D)
    static constexpr int Dim = DIM;         ///< Dimension of tree
    static constexpr int Order = ORDER;     ///< Order of tree

    using node_t = Node<Func, ORDER, ISET>;      ///< DIM,ORDER node type
    using box_t = Box<DIM, ISET, value_type>;    ///< DIM box type
    using VecDimD = std::array<value_type, DIM>; ///< D dimensional vector type

    std::vector<node_t> nodes_; ///< Flat list of all nodes in Tree (leaf or otherwise)
    int max_depth_;             ///< Maximum depth of tree

    /// @brief Construct tree
    /// @param[in] input parameters for fit (function, tol, etc)
    /// @param[in] coeffs flat/global coefficient vector
    /// @param[in] box box that this tree lives in
    FunctionTree<ORDER, Func>(const baobzi_input_t &input, const Box<DIM, ISET, value_type> &box,
                              std::vector<PolyEvalType> &polyfits, const Func &func) {
        std::queue<Box<DIM, ISET, input_type>> q;
        VecDimD half_width = scale(box.half_length(), 0.5);
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
                std::vector new_polyfits = node.fit(input, func, {});

                if (node.is_leaf()) {
                    node.poly_eval_id = polyfits.size();
                    for (auto &pf : new_polyfits)
                        polyfits.emplace_back(std::move(pf));
                } else if (!node.is_leaf()) {
                    node.first_child_idx = curr_child_idx;
                    curr_child_idx += NChild;

                    const VecDimD &center = node.box_.center;
                    for (index_t child = 0; child < NChild; ++child) {
                        VecDimD center_offset;

                        // Extract sign of each offset component from the bits of child
                        // Basically: permute all possible offsets
                        for (int j = 0; j < DIM; ++j) {
                            value_type signed_hw[2] = {-half_width[j], half_width[j]};
                            center_offset[j] = center[i] + signed_hw[(child >> j) & 1];
                        }

                        q.push(Box<DIM, ISET, value_type>(center_offset, half_width));
                    }
                }
            }

            if (!q.empty())
                max_depth_++;
            if (max_depth_ > input.max_depth)
                throw MaxDepthExceeded();

            half_width = scale(half_width, 0.5);
        }
    }

    FunctionTree<ORDER, Func>() = default; ///< Default constructor for msgpack happiness

    /// @brief Find leaf node containing a point via standard pointer traversal
    /// @param[in] x point that the node will contain
    /// @return leaf node containing point x
    inline const node_t &find_node_traverse(const input_type &x) const { return nodes_[get_node_index(x)]; }

    /// @brief Get index of node at point x (relative to local nodes_ array)
    /// @param[in] x [DIM] point to lookup
    /// @returns index of node in nodes_ array containing x
    inline std::size_t get_node_index(const input_type &x) const {
        index_t curr_index = 0;
        while (!nodes_[curr_index].is_leaf()) {
            index_t child_idx = 0;

            if constexpr (has_tuple_size_v<input_type>)
                for (int i = 0; i < DIM; ++i)
                    child_idx = child_idx | ((x[i] > nodes_[curr_index].box_.center[i]) << i);
            else
                child_idx = (x > nodes_[curr_index].box_.center[0]) ? 1 : 0;

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
    inline value_type eval(const VecDimD &x, const std::vector<PolyEvalType> &polyfits) const {
        return polyfits[get_node_index(x)](x);
    }

    /// @brief msgpack serialization magic
    MSGPACK_DEFINE(nodes_);
};

/// @brief Represents a function in some domain as a grid of baobzi::FunctionTree objects
/// @tparam DIM dimension of function
/// @tparam ORDER order of evaluation polynomial
/// @tparam ISET instruction set index (dummy variable to force alignment for different instruction sets)
template <int ORDER, class Func>
class Function {
  public:
    using input_type = poly_eval::function_traits<Func>::arg0_type;
    using output_type = poly_eval::function_traits<Func>::result_type;
    using value_type = value_type_or_identity<input_type>::type;
    using PolyEvalType = poly_eval::FuncEval<Func, ORDER>;
    static constexpr int DIM = get_tuple_size<input_type>();
    static constexpr int NChild = 1 << DIM; ///< Number of children each node potentially has (2^D)
    static constexpr int Dim = DIM;         ///< Dimension of tree
    static constexpr int Order = ORDER;     ///< Order of tree
    static constexpr int ISET = 0;
    static constexpr int ISet = ISET; ///< Instruction set (dummy param)

    using node_t = Node<Func, ORDER, ISET>;      ///< DIM,ORDER node type
    using box_t = Box<DIM, ISET, value_type>;    ///< DIM box type
    using VecDimD = std::array<value_type, DIM>; ///< D dimensional vector type

    uint32_t output_dim_ = 1;

    baobzi_input_t input_;
    box_t box_;           ///< box representing the domain of our function
    value_type tol_;      ///< Desired relative tolerance of our approximation
    VecDimD lower_left_;  ///< Bottom 'corner' of our domain
    VecDimD upper_right_; ///< Upper 'corner' of our domain

    std::vector<FunctionTree<ORDER, Func>> subtrees_; ///< Grid of FunctionTree objects that do the work
    std::array<int, DIM> n_subtrees_;                 ///< Number of subtrees in each linear dimension of our space
    std::vector<int> subtree_node_offsets_; ///< n_subtrees array of offsets for where in the global array of node
                                            ///< pointers the global node pointer array starts
    std::vector<node_t *> node_pointers_;   ///< Vector of pointers to every node from every subtree
    VecDimD inv_bin_size_;                  ///< Inverse linear dimensions of the bins that our subtrees live

    std::vector<PolyEvalType> polyfits_; ///< Flat vector of all chebyshev coefficients from all leaf nodes

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
        mem += subtree_node_offsets_.capacity() * sizeof(typename decltype(subtree_node_offsets_)::value_type);
        mem += node_pointers_.capacity() * sizeof(node_t *);
        mem += polyfits_.capacity() * sizeof(PolyEvalType);
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
        std::cout << "Approximate memory usage of tree: " << (value_type)mem / (1024 * 1024) << " MiB" << std::endl;
    }

    /// @brief Construct our Function object (fits recursively, can be slow)
    /// @param[in] input parameters for fit (function, tol, etc)
    /// @param[in] xp [dim] center of function domain
    /// @param[in] lp [dim] half length of function domain
    /// @param[in] samples list of points to force fit check
    Function<Func, ORDER, ISET>(const baobzi_input_t &input, const input_type &xp, const input_type &lp,
                                const Func &func)
        : box_(VecDimD{xp}, VecDimD{lp}), tol_(input.tol), split_multi_eval_(input.split_multi_eval),
          output_dim_(input.output_dim), input_(input) {
        auto t_start = std::chrono::steady_clock::now();

        VecDimD l, x;
        if constexpr (has_tuple_size_v<input_type>) {
            l = lp;
            x = xp;
        } else {
            l[0] = lp;
            x[0] = xp;
        }

        std::queue<box_t> q;
        std::queue<box_t> maybe_q;

        auto lmin = *std::min_element(l.begin(), l.end());
        for (int i = 0; i < DIM; ++i)
            n_subtrees_[i] = l[i] / lmin;

        q.push(box_t(x, l));

        // Half-width of next children
        VecDimD half_width = scale(l, 0.5);

        // Breadth first search. Step through each level of the tree and test fit all of the nodes
        // We exit when a level isn't completely filled with parent nodes (rather than leaves)
        // This way we can always avoid redundant traversals by jumping straight to a root node of a subtree
        while (!q.empty()) {
            int n_next = q.size();

            auto add_node_children_to_queue = [](std::queue<box_t> &theq, const VecDimD &center,
                                                 const VecDimD &half_width) {
                for (unsigned child = 0; child < NChild; ++child) {
                    VecDimD offset_center;

                    // Extract sign of each offset component from the bits of child
                    // Basically: permute all possible offsets
                    for (int j = 0; j < DIM; ++j) {
                        value_type signed_hw[2] = {-half_width[j], half_width[j]};
                        offset_center[j] = center[j] + signed_hw[(child >> j) & 1];
                    }

                    theq.push(box_t(offset_center, half_width));
                }
            };

            std::vector<node_t> nodes;
            value_type leaf_fraction = 0.0;
            for (int i = 0; i < n_next; ++i) {
                box_t box = q.front();
                q.pop();

                nodes.emplace_back(node_t(box));
                auto &node = nodes.back();
                node.fit(input, func, {});

                if (!node.is_leaf() || stats_.base_depth < input.min_depth) {
                    add_node_children_to_queue(q, node.box_.center, half_width);
                } else {
                    leaf_fraction += 1.0;
                    add_node_children_to_queue(maybe_q, node.box_.center, half_width);
                }
            }
            stats_.n_evals_root += nodes.size() * std::pow(ORDER, DIM);

            leaf_fraction /= nodes.size();
            if (leaf_fraction < input.minimum_leaf_fraction) {
                while (!maybe_q.empty()) {
                    box_t box = maybe_q.front();
                    maybe_q.pop();
                    q.push(box);
                }
            }

            half_width = scale(half_width, 0.5);
            if ((1 << (DIM * (stats_.base_depth + 1))) == q.size()) {
                n_subtrees_ = scale(n_subtrees_, 2);
                stats_.base_depth++;
                if (stats_.base_depth > input.max_depth)
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
        for (int i = 0; i < DIM; ++i) {
            lower_left_[i] = box_.center[i] - half_length[i];
            upper_right_[i] = box_.center[i] + half_length[i];
        }

        subtrees_.reserve(prod(n_subtrees_));

        auto input_local = input;
        input_local.max_depth -= stats_.base_depth;
        for (int i_bin = 0; i_bin < prod(n_subtrees_); ++i_bin) {
            std::array<int, DIM> bins = get_bins(i_bin);

            VecDimD parent_center;
            for (int i = 0; i < DIM; ++i)
                parent_center[i] = (bins[i] + value_type{0.5}) * bin_size[i] + lower_left_[i];

            Box<DIM, ISET, value_type> root_box = {parent_center, scale(bin_size, 0.5)};
            subtrees_.emplace_back(input_local, root_box, polyfits_, func);
        }

        auto t_end = std::chrono::steady_clock::now();
        auto t_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
        stats_.t_elapsed = t_elapsed.count();
        build_cache();
    }

    /// @brief Build any intermediate state necessary for computation
    void build_cache() {
        subtree_node_offsets_.resize(prod(n_subtrees_));
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
    Function<Func, DIM, ORDER, ISET>() = default;

    /// @brief convert linear bin index to [dim] bin vector
    /// @param[in] i_bin linear index
    /// @returns [dim] bin vector
    inline std::array<int, DIM> get_bins(const int i_bin) const {
        if constexpr (DIM == 1)
            return std::array<int, DIM>{i_bin};
        else if constexpr (DIM == 2)
            return std::array<int, DIM>{i_bin % n_subtrees_[0], i_bin / n_subtrees_[0]};
        else if constexpr (DIM == 3)
            return std::array<int, DIM>{i_bin % n_subtrees_[0], (i_bin / n_subtrees_[0]) % n_subtrees_[1],
                                        i_bin / (n_subtrees_[0] * n_subtrees_[1])};
    }

    /// @brief find linear index of bin at a point
    /// @param[in] x [1] position to find bin
    /// @returns linear index of bin that x lives in
    inline int get_linear_bin(const input_type &x) const {
        if constexpr (DIM == 1) {
            const value_type x_bin = [this, &x]() {
                if constexpr (has_tuple_size_v<input_type>)
                    return x[0] - lower_left_[0];
                else
                    return x - lower_left_[0];
            }();
            return x_bin * inv_bin_size_[0];
        } else {
            std::array<int, DIM> bin;
            for (int i = 0; i < DIM; ++i)
                bin[i] = (x[i] - lower_left_[i]) * inv_bin_size_[i];

            if constexpr (DIM == 2)
                return bin[0] + n_subtrees_[0] * bin[1];
            else if constexpr (DIM == 3)
                return bin[0] + n_subtrees_[0] * bin[1] + n_subtrees_[0] * n_subtrees_[1] * bin[2];
        }
    }

    /// @brief get constant reference to leaf node that contains a point
    /// @param[in] x point of interest
    /// @returns constant reference to leaf node that contains x
    inline const node_t &find_node(const input_type &x) const {
        return subtrees_[get_linear_bin(x)].find_node_traverse(x);
    }

    inline output_type eval(input_type x) const {
        for (int i = 0; i < DIM; ++i) {
            if (x < lower_left_[i] || x >= upper_right_[i])
                return NAN;
        }
        return polyfits_[find_node(x).poly_eval_id](x);
    }

    /// @brief get index of node (across all subnodes)
    /// @param[in] x [DIM] point to find the node of
    /// @returns index in global node array
    inline std::size_t get_global_node_index(const input_type &x) const {
        const int i_sub = get_linear_bin(x);
        return subtree_node_offsets_[i_sub] + subtrees_[i_sub].get_node_index(x);
    }

    inline void eval(const input_type &x, output_type *res) const {
        if ((x.array() < lower_left_.array()).any() || (x.array() >= upper_right_.array()).any()) {
            for (int i = 0; i < output_dim_; ++i)
                res[i] = NAN;
        }

        *res = find_node(x).eval(x);
    }

    /// @brief eval function approximation at n_trg points
    /// @param[in] xp [DIM * n_trg] array of points to evaluate function at
    /// @param[out] res [n_trg] array of results
    /// @param[in] n_trg number of points to evaluate
    inline void eval(const value_type *xp, value_type *res, int n_trg) const {
        if (split_multi_eval_) {
            std::vector<std::pair<node_t *, input_type>> node_map(n_trg);
            for (int i = 0; i < n_trg; ++i) {
                value_type xi = *(xp + DIM * i);
                node_t *node_ptr = [xi, this]() {
                    if (xi < lower_left_[0] || xi >= upper_right_[0])
                        return (node_t *)nullptr;
                    return node_pointers_[get_global_node_index(xi)];
                }();

                node_map[i] = std::make_pair(node_ptr, xi);
            }

            for (int i_trg = 0; i_trg < n_trg; i_trg++) {
                res[i_trg] = node_map[i_trg].first == nullptr
                                 ? NAN
                                 : polyfits_[node_map[i_trg].first->poly_eval_id](node_map[i_trg].second);
            }
        } else {
            for (int i_trg = 0; i_trg < n_trg; i_trg++)
                res[i_trg] = eval(*(xp + DIM * i_trg));
        }
    }

    inline void operator()(const value_type *xp, value_type *res, int n_trg) const { eval(xp, res, n_trg); }

    /// @brief eval function approximation at point
    /// @param[in] x [DIM] point to evaluate function at
    /// @returns function approximation at point x
    inline output_type operator()(const input_type &x) const { return output_dim_ > 1 ? NAN : eval(x); }

    // /// @brief save function approximation to file
    // /// @param[in] filename path to save file at
    // void save(const char *filename) const {
    //     std::ofstream ofs(filename, std::ofstream::binary | std::ofstream::out);
    //     baobzi_header_t params{Dim, Order, BAOBZI_HEADER_VERSION};
    //     msgpack::pack(ofs, params);
    //     msgpack::pack(ofs, *this);
    // }

    // std::vector<raw_leaf_node> get_leaves() const {
    //     std::vector<raw_leaf_node> leaves;

    //     for (const auto &subtree : subtrees_) {
    //         for (const auto &node : subtree.nodes_) {
    //             if (!node.is_leaf())
    //                 continue;

    //             double L = 2.0 * node.box_.half_length()[0];
    //             double a = node.box_.center[0] - 0.5 * L;
    //             const double *coeffs = node.coeff_offset + coeffs_.data();
    //             leaves.emplace_back(raw_leaf_node{a, L, coeffs});
    //         }
    //     }

    //     std::sort(leaves.begin(), leaves.end(), leaf_compare());
    //     return leaves;
    // }

    std::pair<VecDimD, VecDimD> get_bounds() const { return std::make_pair(lower_left_, upper_right_); }

    // Function<DIM, ORDER, ISET> shallow_copy() const {
    //     Function<DIM, ORDER, ISET> other;
    //     other.n_subtrees_ = n_subtrees_;
    //     other.lower_left_ = lower_left_;
    //     other.inv_bin_size_ = inv_bin_size_;
    //     other.box_ = box_;
    //     other.inv_bin_size_ = inv_bin_size_;
    //     other.split_multi_eval_ = split_multi_eval_;
    //     other.input_ = input_;
    //     other.input_.func = nullptr;
    //     other.input_.data = nullptr;

    //     return other;
    // }

    /// @brief msgpack serialization magic
    MSGPACK_DEFINE_MAP(box_, subtrees_, n_subtrees_, tol_, lower_left_, upper_right_, inv_bin_size_, polyfits_,
                       split_multi_eval_, output_dim_);
};
} // namespace baobzi

#endif
