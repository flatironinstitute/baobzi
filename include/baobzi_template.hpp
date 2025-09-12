#ifndef BAOBZI_TEMPLATE_HPP
#define BAOBZI_TEMPLATE_HPP

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

#include <polyfit/fast_eval.hpp>

#include <baobzi.h>

namespace baobzi {

class MaxDepthExceeded : public std::exception {
    virtual const char *what() const throw() { return "Baobzi fit error: tree depth exceeded max allowed input depth"; }
};

namespace detail {
/// @brief Wrapper class treat arrays and scalars as a common type with basic arithmetic
/// operators
///
/// @tparam T type to extract value_type from
/// @tparam N number of elements. If 1, T is scalar, otherwise tuple-like
template <typename T, std::size_t N>
class Value {
    using storage_t = std::conditional_t<N == 1, T, std::array<T, N>>;
    storage_t data_;

  public:
    // Scalar constructor
    template <std::size_t M = N, typename = std::enable_if_t<M == 1>>
    Value(const T &val) : data_(val) {}
    Value(const std::array<T, 1> &arr) { data_ = arr[0]; }

    // Array constructor
    template <std::size_t M = N, typename = std::enable_if_t<M != 1>>
    Value(const std::array<T, N> &arr) : data_(arr) {}
    Value() = default;

    // Assignment
    inline Value &operator=(const Value &other) {
        data_ = other.data_;
        return *this;
    }

    // Arithmetic operators
    inline Value operator+(const Value &rhs) const {
        if constexpr (N == 1) {
            return Value(data_ + rhs.data_);
        } else {
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i)
                result[i] = data_[i] + rhs.data_[i];
            return Value(result);
        }
    }

    inline Value operator+(const T &rhs) const {
        if constexpr (N == 1) {
            return Value(data_ + rhs);
        } else {
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i)
                result[i] = data_[i] + rhs;
            return Value(result);
        }
    }

    inline Value operator-(const Value &rhs) const {
        if constexpr (N == 1) {
            return Value(data_ - rhs.data_);
        } else {
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i)
                result[i] = data_[i] - rhs.data_[i];
            return Value(result);
        }
    }

    inline Value operator-(const T &rhs) const {
        if constexpr (N == 1) {
            return Value(data_ - rhs);
        } else {
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i)
                result[i] = data_[i] - rhs;
            return Value(result);
        }
    }

    inline Value operator*(const Value &rhs) const {
        if constexpr (N == 1) {
            return Value(data_ * rhs.data_);
        } else {
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i)
                result[i] = data_[i] * rhs.data_[i];
            return Value(result);
        }
    }

    inline Value operator*(const T &rhs) const {
        if constexpr (N == 1) {
            return Value(data_ * rhs);
        } else {
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i)
                result[i] = data_[i] * rhs;
            return Value(result);
        }
    }

    inline Value operator/(const Value &rhs) const {
        if constexpr (N == 1) {
            return Value(data_ / rhs.data_);
        } else {
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i)
                result[i] = data_[i] / rhs.data_[i];
            return Value(result);
        }
    }

    inline Value operator/(const T &rhs) const {
        if constexpr (N == 1) {
            return Value(data_ / rhs);
        } else {
            std::array<T, N> result;
            for (std::size_t i = 0; i < N; ++i)
                result[i] = data_[i] / rhs;
            return Value(result);
        }
    }

    inline T &operator[](size_t idx) {
        if constexpr (N == 1) {
            return data_;
        } else {
            return data_[idx];
        }
    }

    inline const T &operator[](size_t idx) const {
        if constexpr (N == 1) {
            return data_;
        } else {
            return data_[idx];
        }
    }

    inline T *begin() {
        if constexpr (N == 1) {
            return &data_;
        } else {
            return data_.data();
        }
    }

    inline T *end() {
        if constexpr (N == 1) {
            return &data_ + 1;
        } else {
            return data_.data() + N;
        }
    }

    inline T prod() const {
        if constexpr (N == 1) {
            return data_;
        } else {
            T result = 1;
            for (const auto &val : data_)
                result *= val;
            return result;
        }
    }

    // Automatic casting to scalar or array
    inline operator T() const {
        static_assert(N == 1, "Can only cast to scalar if N == 1");
        return data_;
    }

    inline operator std::array<T, N>() const {
        static_assert(N != 1, "Can only cast to array if N != 1");
        return data_;
    }

    // Access underlying data
    inline const T &scalar() const {
        static_assert(N == 1, "Not a scalar");
        return data_;
    }
    inline const std::array<T, N> &array() const {
        static_assert(N != 1, "Not an array");
        return data_;
    }

    inline const std::array<T, N> as_array() const {
        if constexpr (N == 1) {
            return std::array<T, N>{data_};
        } else {
            return data_;
        }
    }
};

template <typename T>
Value(const T &) -> Value<T, 1>;

template <typename T, std::size_t N>
Value(const std::array<T, N> &) -> Value<T, N>;

using index_t = uint32_t; ///< Type specifying indexing into flattened tree

template <std::size_t Degree, class Func>
class Function;

template <int EXP, typename T>
constexpr T powi(T base) {
    if constexpr (EXP == 0) {
        return T{1};
    } else if constexpr (EXP % 2 == 0) {
        const auto half = powi<EXP / 2>(base);
        return half * half;
    } else {
        return base * powi<EXP - 1>(base);
    }
}

/// @brief Trait to extract number of elements from input/output types
/// @tparam T type to extract size from
/// @returns number of elements in tuple-like type, or 1 for scalars
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
template <typename T, int Dim>
struct Box {
    const Value<T, Dim> center;      ///< Center of box
    const Value<T, Dim> half_length; ///< half the dimension of the box

    /// @brief Constructor, just copies x, hl over
    Box(const auto &x, const auto &hl) : center{x}, half_length{hl} {}
};

/// @brief Check if the 'tail estimate' of the error for a given set of coefficients is greater
/// than tol.
///
/// @tparam Polyfit polynomial fit of function as this node
/// @param[in] tol_type type of tolerance to use for error check. See baobzi_tol_t
/// @param[in] tol tolerance value to use for error check
/// @param[in] polyfit polynomial fit to evaluate at sample points
/// @returns true if estimated error greater than tol, false otherwise
template <class Polyfit>
inline bool tail_error_check(baobzi_tol_t tol_type, double tol, const Polyfit &polyfit) {
    constexpr int input_dim = get_tuple_size<typename Polyfit::InputType>();
    constexpr int output_dim = get_tuple_size<typename Polyfit::OutputType>();
    using T = value_type_or_identity<typename Polyfit::InputType>::type;
    if (input_dim > 2)
        throw std::runtime_error("Baobzi fit error: tail_error only implemented for 1D and 2D input");

    T maxcoeff{0.0};
    T scaling_factor{1.0};

    if constexpr (input_dim == 1) {
        const auto &coeffs = polyfit.coeffs();
        constexpr int N = Polyfit::kDegreeCompileTime;
        static_assert(output_dim == 1, "tail_error only implemented for single output in 1D");

        for (auto i = 0; i < 2; ++i)
            maxcoeff = std::max(std::abs(coeffs[i]), maxcoeff);
        scaling_factor = std::max(scaling_factor, std::abs(coeffs[N - 1]));
    } else if constexpr (input_dim == 2) {
        const int n = polyfit.degree();
        for (int i_dim = 0; i_dim < output_dim; ++i_dim) {
            for (auto i = 0; i < n; ++i)
                maxcoeff = std::max(std::abs(polyfit.coeff_at(i_dim, i, n - i - 1)), maxcoeff);

            scaling_factor = std::max(scaling_factor, std::abs(polyfit.coeff_at(i_dim, n - 1, 0)));
            scaling_factor = std::max(scaling_factor, std::abs(polyfit.coeff_at(i_dim, 0, n - 1)));
        }
    }

    if (tol_type == BAOBZI_TOL_RELATIVE_L2 || tol_type == BAOBZI_TOL_RELATIVE_MAX)
        return maxcoeff / scaling_factor > tol;
    else
        return maxcoeff > tol;
}

/// @brief Sample a polynomial fit on a uniform grid of points and compare to actual
/// function. Poor fit indicated by true, i.e. returns if measured error greater than tol.
///
/// @tparam Func (actual) input function type to evaluate at this node
/// @tparam Polyfit polynomial fit of function as this node
/// @param[in] n_sample_1d number of samples per dimension (total samples = n_sample_1d^input_dim)
/// @param[in] tol_type type of tolerance to use for error check. See baobzi_tol_t
/// @param[in] tol tolerance value to use for error check
/// @param[in] center_in [input_dim] center of box to sample in
/// @param[in] half_length_in [input_dim] half length of box to sample in
/// @param[in] func (actual) function to evaluate at sample points
/// @param[in] polyfit polynomial fit to evaluate at sample points
/// @returns true if sampled error greater than tol, false otherwise
template <class Func, class Polyfit>
inline bool
sample_error_check(int n_sample_1d, baobzi_tol_t tol_type, double tol, const typename Polyfit::InputType &center_in,
                   const typename Polyfit::InputType &half_length_in, const Func &func, const Polyfit &polyfit) {
    constexpr auto input_dim = get_tuple_size<typename Polyfit::InputType>();
    constexpr auto output_dim = get_tuple_size<typename Polyfit::OutputType>();
    const int n_samples = powi<input_dim>(n_sample_1d);
    const Value half_length = half_length_in;
    const Value center = center_in;

    double max_abs_err{0.0}, max_rel_err{0.0}, abs_err_l2{0.0}, direct_sum{0.0};
    for (int linear_index = 0; linear_index < n_samples; ++linear_index) {
        Value<double, input_dim> sample_point;
        int curr_index = linear_index;

        for (int dim = 0; dim < input_dim; ++dim) {
            const double dx = 2.0 * half_length[dim] / n_sample_1d;
            sample_point[dim] = center[dim] - half_length[dim] + dx / 2.0 + dx * (curr_index % n_sample_1d);
            curr_index /= n_sample_1d;
        }

        Value<double, output_dim> actual = func(sample_point);
        Value<double, output_dim> approx = polyfit(sample_point);

        for (int i = 0; i < output_dim; ++i) {
            const double abs_err = std::abs(approx[i] - actual[i]);
            max_abs_err = std::max((double)max_abs_err, abs_err);
            if (actual[i] != 0.0)
                max_rel_err = std::max(max_rel_err, std::abs(abs_err / actual[i]));
            abs_err_l2 += powi<2>(abs_err);
            direct_sum += powi<2>(actual[i]);
        }
    }

    switch (tol_type) {
    case BAOBZI_TOL_RELATIVE_L2:
        return std::sqrt(abs_err_l2 / direct_sum) > tol;
    case BAOBZI_TOL_ABSOLUTE_L2:
        return std::sqrt(abs_err_l2) / (n_samples * output_dim) > tol;
    case BAOBZI_TOL_RELATIVE_MAX:
        return max_rel_err > tol;
    case BAOBZI_TOL_ABSOLUTE_MAX:
        return max_abs_err > tol;
    default:
        throw std::runtime_error("Baobzi fit error: unknown tolerance type for sampling");
    }
}

/// @brief Node in baobzi::PolyTree. If leaf, contains evaluation data, otherwise children
/// @tparam Func function type to evaluate at this node
/// @tparam Degree of evaluation polynomial
template <class Func, std::size_t Degree>
class Node {
  public:
    using input_type = typename std::remove_cvref_t<typename poly_eval::function_traits<Func>::arg0_type>;
    using output_type = poly_eval::function_traits<Func>::result_type;
    using value_type = value_type_or_identity<input_type>::type;
    using poly_eval_type = std::conditional<has_tuple_size_v<input_type>, poly_eval::FuncEvalND<Func, Degree>,
                                            poly_eval::FuncEval<Func, Degree>>::type;

    static constexpr int input_dim = get_tuple_size<input_type>();
    static constexpr int output_dim = get_tuple_size<output_type>();

    Value<value_type, input_dim> center;                          ///< Center of the node
    uint64_t poly_eval_id = std::numeric_limits<uint64_t>::max(); ///< Position of poly_eval object in global array
    uint32_t first_child_idx = -1; ///< First child's index in a flattened list of all nodes

    /// @brief Construct node from box (without fitting)
    /// @param [in] box box this node represents
    Node(const Box<value_type, input_dim> &box) : center{box.center} {}

    /// @brief check if node is leaf
    /// @return true if leaf, false otherwise
    inline bool is_leaf() const { return poly_eval_id != std::numeric_limits<uint64_t>::max(); }

    /// @brief Fit node to a given tolerance. If fit succeeds, set leaf and coeffs, otherwise ... don't
    ///
    /// @param[in] input parameters for fit (function, tol, etc)
    /// @returns coefficient vector list if fit successful, empty list if not good enough
    bool fit(const baobzi_input_t &input, const Func &func, const Value<value_type, input_dim> &half_length,
             const std::vector<value_type> &samples, std::vector<poly_eval_type> &polyfits) {
        if (samples.size())
            throw std::runtime_error("Baobzi fit error: sample points not yet supported");

        const auto n_polyfit_before = polyfits.size();
        const input_type lb = center - half_length;
        const input_type ub = center + half_length;

        auto rollback_and_fail = [&polyfits, n_polyfit_before]() {
            while (polyfits.size() != n_polyfit_before)
                polyfits.pop_back();
            return false;
        };

        auto polyfit = polyfits.emplace_back(func, lb, ub);

        if (input.tol_type == BAOBZI_TOL_RELATIVE_TAIL || input.tol_type == BAOBZI_TOL_ABSOLUTE_TAIL) {
            if (tail_error_check(input.tol_type, input.tol, polyfit))
                return rollback_and_fail();
        } else {
            if (sample_error_check(input.n_samples_per_dim, input.tol_type, input.tol, center, half_length, func,
                                   polyfit))
                return rollback_and_fail();
        }

        poly_eval_id = n_polyfit_before;
        return true;
    }

    /// @brief Calculate memory usage of self (including unused space from vector allocation)
    /// @returns size in bytes of object instance
    inline std::size_t memory_usage() const { return sizeof(*this); }
};

/// @brief Represent a function in some domain as a tree of chebyshev nodes
/// @tparam Degree of evaluation polynomial
/// @tparam Func input function type to fit
template <std::size_t Degree, class Func>
struct PolyTree {
    using input_type = typename std::remove_cvref_t<typename poly_eval::function_traits<Func>::arg0_type>;
    using value_type = value_type_or_identity<input_type>::type;
    using poly_eval_type = std::conditional<has_tuple_size_v<input_type>, poly_eval::FuncEvalND<Func, Degree>,
                                            poly_eval::FuncEval<Func, Degree>>::type;
    using output_type = poly_eval::function_traits<Func>::result_type;

    static constexpr int output_dim = get_tuple_size<output_type>();
    static constexpr int input_dim = get_tuple_size<input_type>();
    static constexpr int n_child = 1 << input_dim; ///< Number of children each node potentially has (2^D)

    using node_t = Node<Func, Degree>;                        ///< Func,Degree node type
    using box_t = Box<value_type, input_dim>;                 ///< input_dim box type
    using dim_array_t = detail::Value<value_type, input_dim>; ///< input_dim dimensional vector type

    /// @brief Construct tree
    /// @param[in] input parameters for fit (function, tol, etc)
    /// @param[in] coeffs flat/global coefficient vector
    /// @param[in] box box that this tree lives in
    PolyTree(const baobzi_input_t &input, const Box<value_type, input_dim> &box, std::vector<poly_eval_type> &polyfits,
             const Func &func) {
        std::queue<box_t> q;
        dim_array_t half_width = box.half_length * 0.5;
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
                const bool successful_fit = node.fit(input, func, box.half_length, {}, polyfits);

                if (successful_fit) {
                    assert(node.poly_eval_id == poly_id);
                    assert(polyfits.size() == poly_id + 1);
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

                        q.push(box_t(center_offset, half_width));
                    }
                }
            }

            if (!q.empty())
                max_depth_++;
            if (max_depth_ > input.max_depth)
                throw MaxDepthExceeded();

            half_width = half_width * 0.5;
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

/// @brief Represents a function in some domain as a grid of baobzi::PolyTree objects
///
/// @tparam Degree of evaluation polynomial
/// @tparam Func function type to evaluate at this node
template <std::size_t Degree, class Func>
class Function {
  public:
    using input_type = typename std::remove_cvref_t<typename poly_eval::function_traits<Func>::arg0_type>;
    using output_type = poly_eval::function_traits<Func>::result_type;
    using value_type = value_type_or_identity<input_type>::type;
    using poly_eval_type = std::conditional<has_tuple_size_v<input_type>, poly_eval::FuncEvalND<Func, Degree>,
                                            poly_eval::FuncEval<Func, Degree>>::type;

    static constexpr int input_dim = detail::get_tuple_size<input_type>();   ///< Function input dimensions
    static constexpr int output_dim = detail::get_tuple_size<output_type>(); ///< Function output dimension
    static constexpr int n_child = 1 << input_dim; ///< Number of children each node potentially has (2^D)

    using node_t = detail::Node<Func, Degree>;                ///< Func,Degree node type
    using box_t = detail::Box<value_type, input_dim>;         ///< input_dim box type
    using dim_array_t = detail::Value<value_type, input_dim>; ///< input_dim dimensional vector type

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

        std::cout << "Baobzi function mapping " << input_dim << " to " << output_dim << std::endl;
        std::cout << "Tree represented by " << n_nodes << " nodes, of which " << n_leaves << " are leaves\n";
        std::cout << "Nodes are distributed across " << n_subtrees << " subtrees at an initial depth of "
                  << stats_.base_depth << " with a maximum subtree depth of " << max_depth << "\n";
        std::cout << "Total function evaluations required for fit: "
                  << n_nodes * (int)std::pow(Degree, input_dim) + stats_.n_evals_root << std::endl;
        std::cout << "Total time to create tree: " << stats_.t_elapsed << " milliseconds\n";
        std::cout << "Approximate memory usage of tree: " << (value_type)mem / (1024 * 1024) << " MiB" << std::endl;
    }

    /// @brief Construct our Function object (fits recursively, can be slow)
    /// @param[in] input parameters for fit (function, tol, etc)
    /// @param[in] center [input_dim] center of function domain
    /// @param[in] half_width_in [input_dim] half length of function domain
    /// @param[in] func function to fit
    Function(const baobzi_input_t &input, const input_type center, const input_type half_width_in, const Func &func)
        : box_(dim_array_t{center}, dim_array_t{half_width_in}), tol_(input.tol),
          split_multi_eval_(input.split_multi_eval), input_(input) {
        auto t_start = std::chrono::steady_clock::now();

        dim_array_t lvec{half_width_in};
        std::queue<box_t> q;
        std::queue<box_t> maybe_q;

        const auto hlmin = *std::min_element(lvec.begin(), lvec.end());
        for (int i = 0; i < input_dim; ++i)
            n_subtrees_[i] = lvec[i] / hlmin;

        q.push(box_t(center, lvec));

        // Half-width of next children
        dim_array_t half_width = lvec * 0.5;

        // Breadth first search. Step through each level of the tree and test fit all of the nodes
        // We exit when a level isn't completely filled with parent nodes (rather than leaves)
        // This way we can always avoid redundant traversals by jumping straight to a root node of a subtree
        while (!q.empty()) {
            int n_next = q.size();

            auto add_node_children_to_queue = [](std::queue<box_t> &theq, const dim_array_t &center,
                                                 const dim_array_t &half_width) {
                for (unsigned child = 0; child < n_child; ++child) {
                    detail::Value<double, input_dim> offset_center;

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
            stats_.n_evals_root += nodes.size() * std::pow(Degree, input_dim);

            leaf_fraction /= nodes.size();
            if (leaf_fraction < input.minimum_leaf_fraction) {
                while (!maybe_q.empty()) {
                    box_t box = maybe_q.front();
                    maybe_q.pop();
                    q.push(box);
                }
            }

            half_width = half_width * 0.5;
            if ((1 << (input_dim * (stats_.base_depth + 1))) == q.size()) {
                n_subtrees_ = n_subtrees_ * 2;
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
        lower_left_ = box_.center - box_.half_length;
        upper_right_ = box_.center + box_.half_length;

        subtrees_.reserve(n_subtrees_.prod());

        auto input_local = input;
        input_local.max_depth -= stats_.base_depth;
        for (int i_bin = 0; i_bin < n_subtrees_.prod(); ++i_bin) {
            std::array<int, input_dim> bins = get_bins(i_bin);

            dim_array_t parent_center;
            for (int i = 0; i < input_dim; ++i)
                parent_center[i] = (bins[i] + value_type{0.5}) * bin_size[i] + lower_left_[i];

            box_t root_box = {parent_center, bin_size * 0.5};
            subtrees_.emplace_back(input_local, root_box, polyfits_, func);
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
    /// @param[out] res [output_dim * n_trg] array of results
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
    baobzi_input_t input_;    ///< copy of input parameters
    box_t box_;               ///< box representing the domain of our function
    value_type tol_;          ///< Desired relative tolerance of our approximation
    dim_array_t lower_left_;  ///< Bottom 'corner' of our domain
    dim_array_t upper_right_; ///< Upper 'corner' of our domain

    std::vector<detail::PolyTree<Degree, Func>> subtrees_; ///< Grid of PolyTree objects that do the work
    detail::Value<int, input_dim> n_subtrees_;             ///< Number of subtrees in each linear dimension of our space
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

/// @brief Factory function to create baobzi::Function object with type deduction
///
/// @tparam Degree Degree of evaluation polynomial
/// @tparam Func Function to fit
/// @param[in] input parameters for fit (function, tol, etc)
/// @param[in] center [M] center of function domain
/// @param[in] half_width_in [M] half length of function domain
/// @param[in] func function to fit. Should map an input from R^M to R^N for some M,N
/// @returns baobzi::Function object representing func in the given domain
template <std::size_t Degree, class Func>
inline Function<Degree, Func>
make_function(const baobzi_input_t &input,
              const std::remove_cvref_t<typename poly_eval::function_traits<Func>::arg0_type> &center,
              const std::remove_cvref_t<typename poly_eval::function_traits<Func>::arg0_type> &half_width_in,
              const Func &func) {
    return Function<Degree, Func>(input, center, half_width_in, func);
}

} // namespace baobzi

#endif
