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
#include <queue>
#include <tuple>
#include <type_traits>
#include <vector>

#include <polyfit/fast_eval.hpp>

#include <baobzi.h>

/// Namespace for baobzi
namespace baobzi {

class MaxDepthExceeded : public std::exception {
    virtual const char *what() const throw() { return "Baobzi fit error: tree depth exceeded max allowed input depth"; }
};

namespace detail {
using index_t = uint32_t; ///< Type specifying indexing into flattened tree

template <std::size_t Order, class Func>
class Function;

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
/// @tparam Dim number of dimensions of box
/// @tparam T type of coordinates (e.g., double, float)
template <int Dim, typename T>
struct Box {
    const std::array<T, Dim> center;      ///< Center of box
    const std::array<T, Dim> half_length; ///< half the dimension of the box

    /// @brief Constructor, just copies x, hl over
    Box(const auto &x, const auto &hl) : center{x}, half_length{hl} {}
};

/// @brief Return an estimate of the error for a given set of coefficients
/// @param[in] coeffs one or two dimensional Vector/Matrix of coefficients
/// @returns estimation of error given those coefficients
inline double tail_error_estimate(int i_dim, const auto &polyfit, baobzi_tol_t tol_type) {
    using input_type = std::remove_cvref_t<decltype(polyfit)>::InputType;
    constexpr int input_dim = get_tuple_size<input_type>();
    static_assert(input_dim == 1 || input_dim == 2, "tail_error_estimate only implemented for 1D and 2D input");
    using T = value_type_or_identity<input_type>::type;

    T maxcoeff{0.0};
    T scaling_factor{1.0};

    if constexpr (input_dim == 1) {
        const auto &coeffs = polyfit.coeffs();

        int n = coeffs.size();
        for (auto i = 0; i < 2; ++i)
            maxcoeff = std::max(std::abs(coeffs[i]), maxcoeff);
        scaling_factor = std::max(scaling_factor, std::abs(coeffs[n - 1]));
    } else if constexpr (input_dim == 2) {
        const int n = polyfit.degree();
        for (auto i = 0; i < n; ++i)
            maxcoeff = std::max(std::abs(polyfit.coeff_at(i_dim, i, n - i - 1)), maxcoeff);

        scaling_factor = std::max(scaling_factor, std::abs(polyfit.coeff_at(i_dim, n - 1, 0)));
        scaling_factor = std::max(scaling_factor, std::abs(polyfit.coeff_at(i_dim, 0, n - 1)));
    }

    if (tol_type == BAOBZI_TOL_RELATIVE)
        return maxcoeff / scaling_factor;
    else
        return maxcoeff;
}

/// @brief Node in baobzi::FunctionTree. If leaf, contains evaluation data, otherwise children
/// @tparam Func function type to evaluate at this node
/// @tparam Order order of evaluation polynomial
template <class Func, std::size_t Order>
class Node {
  public:
    using input_type_cv = typename poly_eval::function_traits<Func>::arg0_type;
    using input_type = typename std::remove_cvref_t<typename poly_eval::function_traits<Func>::arg0_type>;
    using output_type = poly_eval::function_traits<Func>::result_type;
    using value_type = value_type_or_identity<input_type>::type;
    using poly_eval_type = std::conditional<has_tuple_size_v<input_type>, poly_eval::FuncEvalND<Func, Order>,
                                            poly_eval::FuncEval<Func, Order>>::type;

    static constexpr int input_dim = get_tuple_size<input_type>();
    static constexpr int output_dim = get_tuple_size<output_type>();

    using dim_array_t = std::array<value_type, input_dim>; ///< input_dim dimensional vector type
    using order_array_t = std::array<value_type, Order>;   ///< Order dimensional vector type

    std::array<value_type, input_dim> center;                     ///< Center of the node
    uint64_t poly_eval_id = std::numeric_limits<uint64_t>::max(); ///< Position of poly_eval object in global array
    uint32_t first_child_idx = -1; ///< First child's index in a flattened list of all nodes

    /// @brief Construct node from box (without fitting)
    /// @param [in] box box this node represents
    Node(const Box<input_dim, value_type> &box) : center{box.center} {}

    /// @brief check if node is leaf
    /// @return true if leaf, false otherwise
    inline bool is_leaf() const { return poly_eval_id != std::numeric_limits<uint64_t>::max(); }

    /// @brief Fit node to a given tolerance. If fit succeeds, set leaf and coeffs, otherwise ... don't
    ///
    /// @param[in] input parameters for fit (function, tol, etc)
    /// @returns coefficient vector list if fit successful, empty list if not good enough
    bool fit(const baobzi_input_t &input, const Func &func, const std::array<value_type, input_dim> &half_length,
             const std::vector<value_type> &samples, std::vector<poly_eval_type> &polyfits) {
        if (samples.size())
            throw std::runtime_error("Baobzi fit error: sample points not yet supported");

        const auto n_polyfit_before = polyfits.size();
        input_type lb, ub;
        if constexpr (input_dim == 1) {
            lb = center[0] - half_length[0];
            ub = center[0] + half_length[0];
        } else {
            for (int i = 0; i < input_dim; ++i) {
                lb[i] = center[i] - half_length[i];
                ub[i] = center[i] + half_length[i];
            }
        }

        auto rollback_and_fail = [&polyfits, n_polyfit_before]() {
            while (polyfits.size() != n_polyfit_before)
                polyfits.pop_back();
            return false;
        };

        for (int i_dim = 0; i_dim < output_dim; ++i_dim) {
            if constexpr (input_dim == 1)
                polyfits.emplace_back(func, lb, ub, nullptr);
            else
                polyfits.emplace_back(func, lb, ub);

            if constexpr (input_dim <= 2) {
                if (tail_error_estimate(i_dim, polyfits.back(), input.tol_type) > input.tol)
                    return rollback_and_fail();
            } else {
                // For higher dimensions, we need to sample for error, as the tail estimate is not
                // implemented. Here we sample uniformly in each dimension.
                constexpr int n_sample_1d = Order;
                for (int linear_index = 0; linear_index < poly_eval::detail::constexpr_power<n_sample_1d, input_dim>();
                     ++linear_index) {
                    std::array<double, input_dim> sample_point;
                    int curr_index = linear_index;
                    for (int dim = 0; dim < input_dim; ++dim) {
                        const double dx = 2.0 * half_length[dim] / 5;
                        sample_point[dim] = center[dim] - half_length[dim] + dx / 2.0 + dx * (curr_index % n_sample_1d);
                        curr_index /= n_sample_1d;
                    }

                    const std::array<double, output_dim> actual = func(sample_point);
                    const std::array<double, output_dim> approx = polyfits.back()(sample_point);
                    for (int j = 0; j < output_dim; ++j)
                        if (std::abs(1.0 - approx[j] / actual[j]) > input.tol)
                            return rollback_and_fail();
                }
            }
        }

        poly_eval_id = n_polyfit_before;
        return true;
    }

    /// @brief Calculate memory usage of self (including unused space from vector allocation)
    /// @returns size in bytes of object instance
    inline std::size_t memory_usage() const { return sizeof(*this); }
};

/// @brief Represent a function in some domain as a tree of chebyshev nodes
/// @tparam Order order of evaluation polynomial
/// @tparam Func input function type to fit
template <std::size_t Order, class Func>
struct FunctionTree {
    using input_type_cv = typename poly_eval::function_traits<Func>::arg0_type;
    using input_type = typename std::remove_cvref_t<typename poly_eval::function_traits<Func>::arg0_type>;
    using value_type = value_type_or_identity<input_type>::type;
    using poly_eval_type = std::conditional<has_tuple_size_v<input_type>, poly_eval::FuncEvalND<Func, Order>,
                                            poly_eval::FuncEval<Func, Order>>::type;
    using output_type = poly_eval::function_traits<Func>::result_type;

    static constexpr int output_dim = get_tuple_size<output_type>();
    static constexpr int input_dim = get_tuple_size<input_type>();
    static constexpr int n_child = 1 << input_dim; ///< Number of children each node potentially has (2^D)

    using node_t = Node<Func, Order>;                      ///< Func,Order node type
    using box_t = Box<input_dim, value_type>;              ///< input_dim box type
    using dim_array_t = std::array<value_type, input_dim>; ///< input_dim dimensional vector type

    /// @brief Construct tree
    /// @param[in] input parameters for fit (function, tol, etc)
    /// @param[in] coeffs flat/global coefficient vector
    /// @param[in] box box that this tree lives in
    FunctionTree(const baobzi_input_t &input, const Box<input_dim, value_type> &box,
                 std::vector<poly_eval_type> &polyfits, const Func &func) {
        std::queue<box_t> q;
        dim_array_t half_width = scale(box.half_length, 0.5);
        q.push(box);

        index_t curr_child_idx = 1;
        max_depth_ = 0;
        while (!q.empty()) {
            int n_next = q.size();
            int node_index = nodes_.size();
            for (int i = 0; i < n_next; ++i) {
                box_t box = q.front();
                q.pop();

                nodes_.emplace_back(box);

                auto &node = nodes_[i + node_index];
                const auto poly_id = polyfits.size();
                bool successful_fit = node.fit(input, func, box.half_length, {}, polyfits);

                if (successful_fit) {
                    assert(node.poly_eval_id == poly_id);
                    assert(polyfits.size() == poly_id + output_dim);
                    if constexpr (input_dim == 2)
                        polyfits.back()(box.center);
                } else {
                    node.first_child_idx = curr_child_idx;
                    curr_child_idx += n_child;

                    const dim_array_t &center = node.center;
                    for (index_t child = 0; child < n_child; ++child) {
                        dim_array_t center_offset;

                        // Extract sign of each offset component from the bits of child
                        // Basically: permute all possible offsets
                        for (int j = 0; j < input_dim; ++j) {
                            value_type signed_hw[2] = {-half_width[j], half_width[j]};
                            center_offset[j] = center[j] + signed_hw[(child >> j) & 1];
                        }

                        q.push(Box<input_dim, value_type>(center_offset, half_width));
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

    /// @brief Find leaf node containing a point via standard pointer traversal
    /// @param[in] x point that the node will contain
    /// @return leaf node containing point x
    inline const node_t &find_node(const input_type &x) const { return nodes_[get_node_index(x)]; }

    /// @brief Get index of node at point x (relative to local nodes_ array)
    /// @param[in] x [input_dim] point to lookup
    /// @returns index of node in nodes_ array containing x
    inline std::size_t get_node_index(const input_type &x) const {
        index_t curr_index = 0;
        while (!nodes_[curr_index].is_leaf()) {
            index_t child_idx = 0;

            if constexpr (has_tuple_size_v<input_type>)
                for (int i = 0; i < input_dim; ++i)
                    child_idx = child_idx | ((x[i] > nodes_[curr_index].center[i]) << i);
            else
                child_idx = x > nodes_[curr_index].center[0];

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

    inline auto &get_nodes() { return nodes_; }
    inline auto &get_nodes() const { return nodes_; }

  private:
    std::vector<node_t> nodes_; ///< Flat list of all nodes in Tree (leaf or otherwise)
    int max_depth_;             ///< Maximum depth of tree
};
} // namespace detail

/// @brief Represents a function in some domain as a grid of baobzi::FunctionTree objects
/// @tparam Order order of evaluation polynomial
/// @tparam Func function type to evaluate at this node
template <std::size_t Order, class Func>
class Function {
  public:
    using input_type_cv = typename poly_eval::function_traits<Func>::arg0_type;
    using input_type = typename std::remove_cvref_t<typename poly_eval::function_traits<Func>::arg0_type>;
    using output_type = poly_eval::function_traits<Func>::result_type;
    using value_type = value_type_or_identity<input_type>::type;
    using poly_eval_type = std::conditional<has_tuple_size_v<input_type>, poly_eval::FuncEvalND<Func, Order>,
                                            poly_eval::FuncEval<Func, Order>>::type;

    static constexpr int input_dim = detail::get_tuple_size<input_type>();
    static constexpr int output_dim = detail::get_tuple_size<output_type>();
    static constexpr int n_child = 1 << input_dim; ///< Number of children each node potentially has (2^D)

    using node_t = detail::Node<Func, Order>;              ///< Func,Order node type
    using box_t = detail::Box<input_dim, value_type>;      ///< input_dim box type
    using dim_array_t = std::array<value_type, input_dim>; ///< input_dim dimensional vector type

    /// @brief Calculate memory_usage of this object in bytes
    /// @returns Memory usage of baobzi object in bytes
    std::size_t memory_usage() const {
        std::size_t mem = sizeof(*this);
        mem += subtree_node_offsets_.capacity() * sizeof(typename decltype(subtree_node_offsets_)::value_type);
        mem += node_pointers_.capacity() * sizeof(node_t *);
        mem += polyfits_.capacity() * sizeof(poly_eval_type);
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
            for (const auto &node : subtree.get_nodes())
                n_leaves += node.is_leaf();
        }

        std::cout << "Baobzi function mapping " << input_dim << " to " << output_dim_ << std::endl;
        std::cout << "Tree represented by " << n_nodes << " nodes, of which " << n_leaves << " are leaves\n";
        std::cout << "Nodes are distributed across " << n_subtrees << " subtrees at an initial depth of "
                  << stats_.base_depth << " with a maximum subtree depth of " << max_depth << "\n";
        std::cout << "Total function evaluations required for fit: "
                  << n_nodes * (int)std::pow(Order, input_dim) + stats_.n_evals_root << std::endl;
        std::cout << "Total time to create tree: " << stats_.t_elapsed << " milliseconds\n";
        std::cout << "Approximate memory usage of tree: " << (value_type)mem / (1024 * 1024) << " MiB" << std::endl;
    }

    /// @brief Construct our Function object (fits recursively, can be slow)
    /// @param[in] input parameters for fit (function, tol, etc)
    /// @param[in] xp [dim] center of function domain
    /// @param[in] lp [dim] half length of function domain
    /// @param[in] samples list of points to force fit check
    Function(const baobzi_input_t &input, const input_type center, const input_type half_width_in, const Func &func)
        : box_(dim_array_t{center}, dim_array_t{half_width_in}), tol_(input.tol),
          split_multi_eval_(input.split_multi_eval), output_dim_(input.output_dim), input_(input) {
        auto t_start = std::chrono::steady_clock::now();

        dim_array_t lvec{half_width_in}, xvec{center};
        std::queue<box_t> q;
        std::queue<box_t> maybe_q;

        auto hlmin = *std::min_element(lvec.begin(), lvec.end());
        for (int i = 0; i < input_dim; ++i)
            n_subtrees_[i] = lvec[i] / hlmin;

        q.push(box_t(xvec, lvec));

        // Half-width of next children
        dim_array_t half_width = detail::scale(lvec, 0.5);

        // Breadth first search. Step through each level of the tree and test fit all of the nodes
        // We exit when a level isn't completely filled with parent nodes (rather than leaves)
        // This way we can always avoid redundant traversals by jumping straight to a root node of a subtree
        while (!q.empty()) {
            int n_next = q.size();

            auto add_node_children_to_queue = [](std::queue<box_t> &theq, const dim_array_t &center,
                                                 const dim_array_t &half_width) {
                for (unsigned child = 0; child < n_child; ++child) {
                    dim_array_t offset_center;

                    // Extract sign of each offset component from the bits of child
                    // Basically: permute all possible offsets
                    for (int j = 0; j < input_dim; ++j) {
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
                std::vector<poly_eval_type> dummy;
                node.fit(input, func, box.half_length, {}, dummy);
                if (node.poly_eval_id)
                    node.poly_eval_id = 0;

                if (!node.is_leaf() || stats_.base_depth < input.min_depth) {
                    add_node_children_to_queue(q, node.center, half_width);
                } else {
                    leaf_fraction += 1.0;
                    add_node_children_to_queue(maybe_q, node.center, half_width);
                }
            }
            stats_.n_evals_root += nodes.size() * std::pow(Order, input_dim);

            leaf_fraction /= nodes.size();
            if (leaf_fraction < input.minimum_leaf_fraction) {
                while (!maybe_q.empty()) {
                    box_t box = maybe_q.front();
                    maybe_q.pop();
                    q.push(box);
                }
            }

            half_width = detail::scale(half_width, 0.5);
            if ((1 << (input_dim * (stats_.base_depth + 1))) == q.size()) {
                n_subtrees_ = detail::scale(n_subtrees_, 2);
                stats_.base_depth++;
                if (stats_.base_depth > input.max_depth)
                    throw MaxDepthExceeded();
            } else
                break;
        }

        dim_array_t bin_size;
        for (int j = 0; j < input_dim; ++j) {
            bin_size[j] = 2.0 * box_.half_length[j] / n_subtrees_[j];
            inv_bin_size_[j] = 0.5 * n_subtrees_[j] / box_.half_length[j];
        }
        for (int i = 0; i < input_dim; ++i) {
            lower_left_[i] = box_.center[i] - box_.half_length[i];
            upper_right_[i] = box_.center[i] + box_.half_length[i];
        }

        subtrees_.reserve(detail::prod(n_subtrees_));

        auto input_local = input;
        input_local.max_depth -= stats_.base_depth;
        for (int i_bin = 0; i_bin < detail::prod(n_subtrees_); ++i_bin) {
            std::array<int, input_dim> bins = get_bins(i_bin);

            dim_array_t parent_center;
            for (int i = 0; i < input_dim; ++i)
                parent_center[i] = (bins[i] + value_type{0.5}) * bin_size[i] + lower_left_[i];

            detail::Box<input_dim, value_type> root_box = {parent_center, detail::scale(bin_size, 0.5)};
            subtrees_.emplace_back(input_local, root_box, polyfits_, func);
        }

        auto t_end = std::chrono::steady_clock::now();
        auto t_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start);
        stats_.t_elapsed = t_elapsed.count();
        build_cache();
    }

    /// @brief Build any intermediate state necessary for computation
    void build_cache() {
        subtree_node_offsets_.resize(detail::prod(n_subtrees_));
        subtree_node_offsets_[0] = 0;
        for (int i = 1; i < subtree_node_offsets_.size(); ++i)
            subtree_node_offsets_[i] = subtree_node_offsets_[i - 1] + subtrees_[i - 1].size();

        auto n_nodes_tot = std::accumulate(subtrees_.begin(), subtrees_.end(), (std::size_t)0,
                                           [](size_t prior, auto &subtree) { return prior + subtree.size(); });

        node_pointers_.resize(n_nodes_tot);

        int i = 0;
        for (auto &subtree : subtrees_)
            for (auto &node : subtree.get_nodes())
                node_pointers_[i++] = &node;
    }

    /// @brief convert linear bin index to [dim] bin vector
    /// @param[in] i_bin linear index
    /// @returns [dim] bin vector
    inline std::array<int, input_dim> get_bins(const int i_bin) const {
        if constexpr (input_dim == 1)
            return std::array<int, input_dim>{i_bin};
        else if constexpr (input_dim == 2)
            return std::array<int, input_dim>{i_bin % n_subtrees_[0], i_bin / n_subtrees_[0]};
        else if constexpr (input_dim == 3)
            return std::array<int, input_dim>{i_bin % n_subtrees_[0], (i_bin / n_subtrees_[0]) % n_subtrees_[1],
                                              i_bin / (n_subtrees_[0] * n_subtrees_[1])};
        else if constexpr (input_dim == 4)
            return std::array<int, input_dim>{i_bin % n_subtrees_[0], (i_bin / n_subtrees_[0]) % n_subtrees_[1],
                                              (i_bin / (n_subtrees_[0] * n_subtrees_[1])) % n_subtrees_[2],
                                              i_bin / (n_subtrees_[0] * n_subtrees_[1] * n_subtrees_[2])};
        else if constexpr (input_dim == 5)
            return std::array<int, input_dim>{
                i_bin % n_subtrees_[0], (i_bin / n_subtrees_[0]) % n_subtrees_[1],
                (i_bin / (n_subtrees_[0] * n_subtrees_[1])) % n_subtrees_[2],
                (i_bin / (n_subtrees_[0] * n_subtrees_[1] * n_subtrees_[2])) % n_subtrees_[3],
                i_bin / (n_subtrees_[0] * n_subtrees_[1] * n_subtrees_[2] * n_subtrees_[3])};
    }

    /// @brief find linear index of bin at a point
    /// @param[in] x [1] position to find bin
    /// @returns linear index of bin that x lives in
    inline int get_linear_bin(const input_type &x) const {
        if constexpr (input_dim == 1) {
            const value_type x_bin = [this, &x]() {
                if constexpr (has_tuple_size_v<input_type>)
                    return x[0] - lower_left_[0];
                else
                    return x - lower_left_[0];
            }();
            return x_bin * inv_bin_size_[0];
        } else {
            std::array<int, input_dim> bin;
            for (int i = 0; i < input_dim; ++i)
                bin[i] = (x[i] - lower_left_[i]) * inv_bin_size_[i];

            if constexpr (input_dim == 2)
                return bin[0] + n_subtrees_[0] * bin[1];
            else if constexpr (input_dim == 3)
                return bin[0] + n_subtrees_[0] * bin[1] + n_subtrees_[0] * n_subtrees_[1] * bin[2];
            else if constexpr (input_dim == 4)
                return bin[0] + n_subtrees_[0] * bin[1] + n_subtrees_[0] * n_subtrees_[1] * bin[2] +
                       n_subtrees_[0] * n_subtrees_[1] * n_subtrees_[2] * bin[3];
            else if constexpr (input_dim == 5)
                return bin[0] + n_subtrees_[0] * bin[1] + n_subtrees_[0] * n_subtrees_[1] * bin[2] +
                       n_subtrees_[0] * n_subtrees_[1] * n_subtrees_[2] * bin[3] +
                       n_subtrees_[0] * n_subtrees_[1] * n_subtrees_[2] * n_subtrees_[3] * bin[4];
        }
    }

    /// @brief get constant reference to leaf node that contains a point
    /// @param[in] x point of interest
    /// @returns constant reference to leaf node that contains x
    inline const node_t &find_node(const input_type &x) const { return subtrees_[get_linear_bin(x)].find_node(x); }

    /// @brief get index of node (across all subnodes)
    /// @param[in] x [input_dim] point to find the node of
    /// @returns index in global node array
    inline std::size_t get_global_node_index(const input_type &x) const {
        const int i_sub = get_linear_bin(x);
        return subtree_node_offsets_[i_sub] + subtrees_[i_sub].get_node_index(x);
    }

    /// @brief eval function approximation at n_trg points
    /// @param[in] xp [input_dim * n_trg] array of points to evaluate function at
    /// @param[out] res [n_trg] array of results
    /// @param[in] n_trg number of points to evaluate
    inline void operator()(const value_type *xp, value_type *res, int n_trg) const {
        if (split_multi_eval_) {
            std::vector<std::pair<node_t *, input_type>> node_map(n_trg);
            for (int i = 0; i < n_trg; ++i) {
                value_type xi = *(xp + input_dim * i);
                node_t *node_ptr = [xi, this]() -> node_t * {
                    if (xi < lower_left_[0] || xi >= upper_right_[0])
                        return nullptr;
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
                res[i_trg] = (*this)(*(xp + input_dim * i_trg));
        }
    }

    /// @brief eval function approximation at point
    /// @param[in] x [input_dim] point to evaluate function at
    /// @returns function approximation at point x
    inline output_type operator()(const input_type &x) const {
        if constexpr (input_dim == 1) {
            if (x < lower_left_[0] || x >= upper_right_[0])
                return NAN;
        } else {
            for (int i = 0; i < input_dim; ++i)
                if (x[i] < lower_left_[i] || x[i] >= upper_right_[i])
                    return output_type{NAN};
        }

        return polyfits_[find_node(x).poly_eval_id](x);
    }

    std::pair<dim_array_t, dim_array_t> get_bounds() const { return std::make_pair(lower_left_, upper_right_); }

  private:
    uint32_t output_dim_ = 1;

    baobzi_input_t input_;
    box_t box_;               ///< box representing the domain of our function
    value_type tol_;          ///< Desired relative tolerance of our approximation
    dim_array_t lower_left_;  ///< Bottom 'corner' of our domain
    dim_array_t upper_right_; ///< Upper 'corner' of our domain

    std::vector<detail::FunctionTree<Order, Func>> subtrees_; ///< Grid of FunctionTree objects that do the work
    std::array<int, input_dim> n_subtrees_; ///< Number of subtrees in each linear dimension of our space
    std::vector<int> subtree_node_offsets_; ///< n_subtrees array of offsets for where in the global array of node
                                            ///< pointers the global node pointer array starts
    std::vector<node_t *> node_pointers_;   ///< Vector of pointers to every node from every subtree
    dim_array_t inv_bin_size_;              ///< Inverse linear dimensions of the bins that our subtrees live

    std::vector<poly_eval_type> polyfits_; ///< Flat vector of all chebyshev coefficients from all leaf nodes

    bool split_multi_eval_ = true; ///< Split node-search and evaluation when evaluating multiple points

    /// Structure containing info about self creation :D
    struct {
        uint16_t base_depth = 0;   ///< depth of subtrees
        uint64_t n_evals_root = 0; ///< number of function evals before subtree calls
        uint32_t t_elapsed = 0;    ///< time in milliseconds to create object
    } stats_;
};

template <std::size_t Order, class Func>
Function<Order, Func>
make_function(const baobzi_input_t &input,
              const std::remove_cvref_t<typename poly_eval::function_traits<Func>::arg0_type> center,
              const std::remove_cvref_t<typename poly_eval::function_traits<Func>::arg0_type> half_width_in,
              const Func &func) {
    return Function<Order, Func>(input, center, half_width_in, func);
}
} // namespace baobzi

#endif
