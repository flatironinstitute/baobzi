#include <baobzi.hpp>

#include <memory>
#include <variant>

namespace baobzi {
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

// List all supported combinations
using FunctionVariant = std::variant<FunctionType<6, 1, 1>, FunctionType<8, 1, 1>, FunctionType<6, 2, 1>,
                                     FunctionType<6, 2, 2>, FunctionType<8, 2, 1>, FunctionType<8, 2, 2>>;

struct baobzi_function {
    FunctionVariant fn;
};

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

typedef enum {
    GENERIC = 0,
    SSE42 = 1,
    AVX2 = 2,
    AVX512 = 3,
} baobzi_isa_t;

template <baobzi_isa_t ISA>
std::unique_ptr<baobzi_function> baobzi_init(const baobzi_input_t *input, const double *center,
                                             const double *half_width_in) {
    // Helper macro for dispatch
#define DISPATCH(D, I, O)                                                                                              \
    if (input->degree == D && input->input_dim == I && input->output_dim == O) {                                       \
        auto lambda = wrap_c_func<I, O>(input->func, input->data);                                                     \
        baobzi::detail::Value<double, I> c, h;                                                                         \
        std::copy(center, center + I, c.begin());                                                                      \
        std::copy(half_width_in, half_width_in + I, h.begin());                                                        \
        return std::make_unique<baobzi_function>(baobzi_function{                                                      \
            FunctionVariant(baobzi::Function<D, decltype(lambda)>(*input, c.get(), h.get(), lambda))});                \
    }

    // clang-format off
    DISPATCH(6, 1, 1) DISPATCH(8, 1, 1)
    DISPATCH(6, 2, 1) DISPATCH(6, 2, 2) DISPATCH(8, 2, 1) DISPATCH(8, 2, 2)
    // clang-format on

#undef DISPATCH
        return nullptr;
}

template <baobzi_isa_t ISA>
void baobzi_free(baobzi_function *f) {
    delete f;
}

template <baobzi_isa_t ISA>
void baobzi_eval(const baobzi_function *f, const double *input_point, double *output_point) {
    std::visit(
        [&](auto &&fn) {
            using fn_t = std::decay_t<decltype(fn)>;
            if constexpr (fn_t::input_dim == 1) {
                double x = *input_point;
                *output_point = fn(x);
            } else {
                typename fn_t::input_type x;
                std::copy(input_point, input_point + x.size(), x.begin());
                auto y = fn(x);
                std::copy(y.begin(), y.end(), output_point);
            }
        },
        f->fn);
}
} // namespace baobzi
