#include <baobzi/baobzi_binding.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

const struct baobzi_input_t baobzi_input_default;

static const auto baobzi_isa = []() {
    using namespace baobzi;
    int iset = baobzi_isa_t::GENERIC;
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
}();

inline std::string file_to_string(const std::string &path) {
    std::ostringstream buf;
    std::ifstream input(path.c_str());
    buf << input.rdbuf();
    return buf.str();
}

extern "C" {

void baobzi_eval(const baobzi_t func, const double *x, double *y) {
#ifdef BAOBZI_CPU_DISPATCH
    switch (baobzi_isa) {
    case baobzi::baobzi_isa_t::GENERIC:
        return baobzi::baobzi_eval_multi<baobzi::baobzi_isa_t::GENERIC>(func, x, y, 1);
    case baobzi::baobzi_isa_t::AVX:
        return baobzi::baobzi_eval_multi<baobzi::baobzi_isa_t::AVX>(func, x, y, 1);
    case baobzi::baobzi_isa_t::AVX2:
        return baobzi::baobzi_eval_multi<baobzi::baobzi_isa_t::AVX2>(func, x, y, 1);
    case baobzi::baobzi_isa_t::AVX512:
        return baobzi::baobzi_eval_multi<baobzi::baobzi_isa_t::AVX512>(func, x, y, 1);
    }
#else
    return baobzi::baobzi_eval<baobzi::baobzi_isa_t::GENERIC>(func, x, y);
#endif
}

void baobzi_eval_multi(const baobzi_t func, const double *x, double *res, int ntrg) {
#ifdef BAOBZI_CPU_DISPATCH
    switch (baobzi_isa) {
    case baobzi::baobzi_isa_t::GENERIC:
        return baobzi::baobzi_eval_multi<baobzi::baobzi_isa_t::GENERIC>(func, x, res, ntrg);
    case baobzi::baobzi_isa_t::AVX:
        return baobzi::baobzi_eval_multi<baobzi::baobzi_isa_t::AVX>(func, x, res, ntrg);
    case baobzi::baobzi_isa_t::AVX2:
        return baobzi::baobzi_eval_multi<baobzi::baobzi_isa_t::AVX2>(func, x, res, ntrg);
    case baobzi::baobzi_isa_t::AVX512:
        return baobzi::baobzi_eval_multi<baobzi::baobzi_isa_t::AVX512>(func, x, res, ntrg);
    }
#else
    return baobzi::baobzi_eval_multi<baobzi::baobzi_isa_t::GENERIC>(func, x, res, ntrg);
#endif
}

void baobzi_stats(baobzi_t func) {
#ifdef BAOBZI_CPU_DISPATCH
    switch (baobzi_isa) {
    case baobzi::baobzi_isa_t::GENERIC:
        return baobzi::baobzi_stats<baobzi::baobzi_isa_t::GENERIC>(func);
    case baobzi::baobzi_isa_t::AVX:
        return baobzi::baobzi_stats<baobzi::baobzi_isa_t::AVX>(func);
    case baobzi::baobzi_isa_t::AVX2:
        return baobzi::baobzi_stats<baobzi::baobzi_isa_t::AVX2>(func);
    case baobzi::baobzi_isa_t::AVX512:
        return baobzi::baobzi_stats<baobzi::baobzi_isa_t::AVX512>(func);
    }
#else
    return baobzi::baobzi_stats<baobzi::baobzi_isa_t::GENERIC>(func);
#endif
}

baobzi_t baobzi_free(baobzi_t func) {
    if (!func)
        return nullptr;
#ifdef BAOBZI_CPU_DISPATCH
    switch (baobzi_isa) {
    case baobzi::baobzi_isa_t::GENERIC:
        return baobzi::baobzi_free<baobzi::baobzi_isa_t::GENERIC>(func);
    case baobzi::baobzi_isa_t::AVX:
        return baobzi::baobzi_free<baobzi::baobzi_isa_t::AVX>(func);
    case baobzi::baobzi_isa_t::AVX2:
        return baobzi::baobzi_free<baobzi::baobzi_isa_t::AVX2>(func);
    case baobzi::baobzi_isa_t::AVX512:
        return baobzi::baobzi_free<baobzi::baobzi_isa_t::AVX512>(func);
    }
    return nullptr;
#else
    return baobzi::baobzi_free<baobzi::baobzi_isa_t::GENERIC>(func);
#endif
}

bool is_valid_func(const baobzi_input_t *input, const double *point) {
    if (!input->func)
        return false;

    double res[input->output_dim];
    try {
        input->func(point, res, input->data);
    } catch (std::exception(e)) {
        return false;
    }

    return true;
}

baobzi_t baobzi_init(const baobzi_input_t *input, const double *center, const double *half_length) {
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
    switch (baobzi_isa) {
    case baobzi::baobzi_isa_t::GENERIC:
        return baobzi::baobzi_init<baobzi::baobzi_isa_t::GENERIC>(input, center, half_length);
    case baobzi::baobzi_isa_t::AVX:
        return baobzi::baobzi_init<baobzi::baobzi_isa_t::AVX>(input, center, half_length);
    case baobzi::baobzi_isa_t::AVX2:
        return baobzi::baobzi_init<baobzi::baobzi_isa_t::AVX2>(input, center, half_length);
    case baobzi::baobzi_isa_t::AVX512:
        return baobzi::baobzi_init<baobzi::baobzi_isa_t::AVX512>(input, center, half_length);
    default:
        std::cerr << "Baobzi error: Unknown CPU instruction set. Using generic\n";
        return baobzi::baobzi_init<baobzi::baobzi_isa_t::GENERIC>(input, center, half_length);
    }
#else
    return baobzi::baobzi_init<baobzi::baobzi_isa_t::GENERIC>(input, center, half_length);
#endif
}
}
