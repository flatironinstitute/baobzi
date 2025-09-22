#ifdef BAOBZI_CPU_DISPATCH
#include <baobzi/baobzi_binding.hpp>
#else
#include <baobzi/baobzi_binding_impl.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include <dlfcn.h>

const struct baobzi_input_t baobzi_input_default;

#ifdef BAOBZI_CPU_DISPATCH
baobzi::baobzi_isa_t get_baobzi_isa() {
    using namespace baobzi;
    auto iset = baobzi_isa_t::GENERIC;
#ifdef BAOBZI_CPU_DISPATCH
    if (__builtin_cpu_supports("avx"))
        iset = baobzi_isa_t::AVX;
    if (__builtin_cpu_supports("avx2"))
        iset = baobzi_isa_t::AVX2;
    if (__builtin_cpu_supports("avx512f"))
        iset = baobzi_isa_t::AVX512;

    const char *iset_str_const = getenv("BAOBZI_ARCH");
    if (iset_str_const) {
        std::string iset_str(iset_str_const);
        std::transform(iset_str.begin(), iset_str.end(), iset_str.begin(),
                       [](unsigned char c) { return std::tolower(c); });

        if (iset_str == "generic")
            iset = baobzi_isa_t::GENERIC;
        else if (iset_str == "sse4.2")
            iset = baobzi_isa_t::AVX;
        else if (iset_str == "avx2")
            iset = baobzi_isa_t::AVX2;
        else if (iset_str == "avx512")
            iset = baobzi_isa_t::AVX512;
        else
            std::cerr << "Error: unable to parse BAOBZI_ARCH. Valid options are: GENERIC, AVX, AVX2, AVX512\n";
    }
#endif

    return iset;
}

std::string isa_to_string(baobzi::baobzi_isa_t isa) {
    switch (isa) {
    case baobzi::baobzi_isa_t::GENERIC:
        return "generic";
    case baobzi::baobzi_isa_t::AVX:
        return "avx";
    case baobzi::baobzi_isa_t::AVX2:
        return "avx2";
    case baobzi::baobzi_isa_t::AVX512:
        return "avx512";
    default:
        throw std::runtime_error("Unknown baobzi ISA");
    }
}

static const auto [baobzi_init_impl, baobzi_eval_impl, baobzi_eval_multi_impl, baobzi_stats_impl,
                   baobzi_free_impl] = []() {
    const auto libstr = "libbaobzi_" + isa_to_string(get_baobzi_isa()) + ".so";
    void *handle = dlopen(libstr.c_str(), RTLD_NOW);
    if (!handle)
        std::cerr << "Error: unable to open " << libstr << "with dlopen: " << dlerror() << "\n";

    auto init_ptr = (decltype(baobzi_init) *)dlsym(handle, "baobzi_init");
    auto eval_ptr = (decltype(baobzi_eval) *)dlsym(handle, "baobzi_eval");
    auto eval_multi_ptr = (decltype(baobzi_eval_multi) *)dlsym(handle, "baobzi_eval_multi");
    auto stats_ptr = (decltype(baobzi_stats) *)dlsym(handle, "baobzi_stats");
    auto free_ptr = (decltype(baobzi_free) *)dlsym(handle, "baobzi_free");

    return std::tuple{init_ptr, eval_ptr, eval_multi_ptr, stats_ptr, free_ptr};
}();
#endif

extern "C" {

baobzi_t baobzi_init(const baobzi_input_t *input, const double *center, const double *half_length) {
    const auto is_valid_func = [](const baobzi_input_t *input, const double *point) {
        if (!input->func)
            return false;

        std::vector<double> res(input->output_dim);
        try {
            input->func(point, res.data(), input->data);
        } catch (std::exception(e)) {
            return false;
        }

        return true;
    };

    if (input->tol <= 0.0) {
        std::cerr << "Baobzi error: Unable to initialize Baobzi due to invalid 'tol' parameter. Please supply "
                     "something greater than zero.\n";
        return nullptr;
    } else if (!is_valid_func(input, center)) {
        std::cerr
            << "BAOBZI ERROR: Unable to initialize Baobzi due to empty or invalid 'func' parameter. Please supply "
               "a valid function to fit.\n";
        return nullptr;
    }

#ifdef BAOBZI_CPU_DISPATCH
    return baobzi_init_impl(input, center, half_length);
#else
    return baobzi::baobzi_init<baobzi::baobzi_isa_t::GENERIC>(input, center, half_length);
#endif
}

void baobzi_eval(const baobzi_t func, const double *x, double *y) {
#ifdef BAOBZI_CPU_DISPATCH
    return baobzi_eval_impl(func, x, y);
#else
    return baobzi::baobzi_eval<baobzi::baobzi_isa_t::GENERIC>(func, x, y);
#endif
}

void baobzi_eval_multi(const baobzi_t func, const double *x, double *res, int ntrg) {
#ifdef BAOBZI_CPU_DISPATCH
    return baobzi_eval_multi_impl(func, x, res, ntrg);
#else
    return baobzi::baobzi_eval_multi<baobzi::baobzi_isa_t::GENERIC>(func, x, res, ntrg);
#endif
}

void baobzi_stats(baobzi_t func) {
#ifdef BAOBZI_CPU_DISPATCH
    return baobzi_stats_impl(func);
#else
    return baobzi::baobzi_stats<baobzi::baobzi_isa_t::GENERIC>(func);
#endif
}

baobzi_t baobzi_free(baobzi_t func) {
    if (!func)
        return nullptr;
#ifdef BAOBZI_CPU_DISPATCH
    return baobzi_free_impl(func);
#else
    return baobzi::baobzi_free<baobzi::baobzi_isa_t::GENERIC>(func);
#endif
}
}
