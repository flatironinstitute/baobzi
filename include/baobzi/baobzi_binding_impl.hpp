#include <baobzi/baobzi_binding.hpp>

#include <variant>

namespace baobzi {
using Degrees = std::integer_sequence<int, 6, 8, 10, 12, 14, 16>;
using InDims = std::integer_sequence<int, 2, 3, 4, 5>;
using OutDims = std::integer_sequence<int, 1, 2, 3, 4, 5>;

template <std::size_t Degree, std::size_t input_dim, std::size_t output_dim>
struct FunctionTypeHelper {
    using type =
        baobzi::Function<Degree, std::function<std::array<double, output_dim>(const std::array<double, input_dim> &)>>;
};

template <std::size_t Degree, std::size_t output_dim>
struct FunctionTypeHelper<Degree, 1, output_dim> {
    static_assert(output_dim == 1, "Only output_dim == 1 is supported for input_dim == 1");
    using type = baobzi::Function<Degree, std::function<double(double)>>;
};

template <std::size_t Degree, std::size_t input_dim, std::size_t output_dim>
using FunctionType = typename FunctionTypeHelper<Degree, input_dim, output_dim>::type;

// Helper to generate all combinations
template <typename... Ts>
struct type_list {};

template <typename... Lists>
struct concat;

template <typename... Ts>
struct concat<type_list<Ts...>> {
    using type = type_list<Ts...>;
};

template <typename... Ts1, typename... Ts2, typename... Rest>
struct concat<type_list<Ts1...>, type_list<Ts2...>, Rest...> {
    using type = typename concat<type_list<Ts1..., Ts2...>, Rest...>::type;
};

// Helper to generate all FunctionType<D, I, O> for a given D
template <int D, int... Is>
struct for_degree {
    template <int I>
    struct for_indim {
        template <int... Os>
        struct for_outdim {
            using type = type_list<FunctionType<D, I, Os>...>;
        };
        template <typename OutDims>
        struct apply;
        template <int... Os>
        struct apply<std::integer_sequence<int, Os...>> {
            using type = typename for_outdim<Os...>::type;
        };
    };

    template <typename InDims, typename OutDims>
    struct apply;

    template <int... I, typename OutDims>
    struct apply<std::integer_sequence<int, I...>, OutDims> {
        using type = typename concat<typename for_indim<I>::template apply<OutDims>::type...>::type;
    };
};

// Main combination generator
template <typename Degrees, typename InDims, typename OutDims>
struct make_combinations;

template <int... Ds, typename InDims, typename OutDims>
struct make_combinations<std::integer_sequence<int, Ds...>, InDims, OutDims> {
    using type = typename concat<typename for_degree<Ds>::template apply<InDims, OutDims>::type...>::type;
};

// All combinations as a type_list
using GeneratedCombinations = make_combinations<Degrees, InDims, OutDims>::type;
using ManualCombinations = type_list<FunctionType<6, 1, 1>, FunctionType<8, 1, 1>, FunctionType<10, 1, 1>,
                                     FunctionType<12, 1, 1>, FunctionType<14, 1, 1>, FunctionType<16, 1, 1>>;
using AllCombinations = concat<ManualCombinations, GeneratedCombinations>::type;

// Convert type_list to std::variant
template <typename>
struct to_variant;
template <typename... Ts>
struct to_variant<type_list<Ts...>> {
    using type = std::variant<Ts...>;
};

// List all supported combinations
using FunctionVariant = to_variant<AllCombinations>::type;
} // namespace baobzi

struct baobzi_function {
    baobzi::FunctionVariant fn;
};

namespace baobzi {
template <std::size_t input_dim, std::size_t output_dim>
struct WrapCFuncType {
    using type = std::function<std::array<double, output_dim>(const std::array<double, input_dim> &)>;
};

template <std::size_t output_dim>
struct WrapCFuncType<1, output_dim> {
    using type = std::function<double(double)>;
};

// For input_dim == 1
template <std::size_t input_dim, std::size_t output_dim>
typename std::enable_if<input_dim == 1, typename WrapCFuncType<input_dim, output_dim>::type>::type
wrap_c_func(baobzi_input_func_t c_func, const void *user_data) {
    return [=](double x) -> double {
        double y;
        c_func(&x, &y, user_data);
        return y;
    };
}

// For input_dim > 1
template <std::size_t input_dim, std::size_t output_dim>
typename std::enable_if<input_dim != 1, typename WrapCFuncType<input_dim, output_dim>::type>::type
wrap_c_func(baobzi_input_func_t c_func, const void *user_data) {
    return [=](const std::array<double, input_dim> &x) -> std::array<double, output_dim> {
        std::array<double, output_dim> y;
        c_func(x.data(), y.data(), user_data);
        return y;
    };
}

template <baobzi_isa_t ISA>
baobzi_t baobzi_free(baobzi_t f) {
    if (f)
        delete f;
    return nullptr;
}

template <typename Variant>
struct variant_types;

template <typename... Ts>
struct variant_types<std::variant<Ts...>> {
    template <typename F>
    static void for_each(F &&f) {
        (f.template operator()<Ts>(), ...);
    }
};

template <typename Variant, typename F>
void for_each_variant_type(F &&f) {
    variant_types<Variant>::for_each(std::forward<F>(f));
}

template <typename Variant, typename F>
bool dispatch(const baobzi_input_t *input, F &&f) {
    bool matched = false;
    for_each_variant_type<Variant>([&]<typename T>() {
        if (input->degree == T::degree && input->input_dim == T::input_dim && input->output_dim == T::output_dim) {
            matched = true;
            f.template operator()<T>();
        }
    });
    return matched;
}

template <baobzi_isa_t ISA>
baobzi_t baobzi_init(const baobzi_input_t *input, const double *center, const double *half_width_in) {
    baobzi_t result = nullptr;

    auto functor = [&]<typename T>() {
        constexpr int D = T::degree;
        constexpr int I = T::input_dim;
        constexpr int O = T::output_dim;
        auto lambda = wrap_c_func<I, O>(input->func, input->data);
        baobzi::detail::Value<double, I> c, h;
        std::copy(center, center + I, c.begin());
        std::copy(half_width_in, half_width_in + I, h.begin());
        result = new baobzi_function{
            FunctionVariant(baobzi::Function<D, decltype(lambda)>(*input, c.get(), h.get(), lambda))};
    };

    const bool found = dispatch<FunctionVariant>(input, functor);
    if (!found)
        std::cerr << "Baobzi: Can't fit function with (degree, input_dim, output_dim) = (" << input->degree << ", "
                  << input->input_dim << ", " << input->output_dim << ")\n";
    return result;
}

template <baobzi_isa_t ISA>
void baobzi_eval_multi(const baobzi_t f, const double *input_points, double *output_points, int ntrg) {
    std::visit(
        [&](auto &&fn) {
            using fn_t = std::decay_t<decltype(fn)>;
            fn(input_points, output_points, ntrg);
        },
        f->fn);
}

template <baobzi_isa_t ISA>
void baobzi_stats(const baobzi_t f) {
    std::visit(
        [&](auto &&fn) {
            using fn_t = std::decay_t<decltype(fn)>;
            fn.print_stats();
        },
        f->fn);
}

} // namespace baobzi
