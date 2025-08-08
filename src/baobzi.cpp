#include "baobzi.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

const struct baobzi_input_t baobzi_input_default;

inline std::string file_to_string(const std::string &path) {
    std::ostringstream buf;
    std::ifstream input(path.c_str());
    buf << input.rdbuf();
    return buf.str();
}

extern "C" {

int get_iset() {
    enum ISET { GENERIC, AVX, AVX2, AVX512 };

    int iset = ISET::GENERIC;
#ifdef BAOBZI_CPU_DISPATCH
    if (__builtin_cpu_supports("avx"))
        iset = ISET::AVX;
    if (__builtin_cpu_supports("avx2"))
        iset = ISET::AVX2;
    if (__builtin_cpu_supports("avx512f"))
        iset = ISET::AVX512;

    const char *iset_str_const = getenv("BAOBZI_ARCH");
    if (iset_str_const) {
        std::string iset_str(iset_str_const);
        std::transform(iset_str.begin(), iset_str.end(), iset_str.begin(),
                       [](unsigned char c) { return std::tolower(c); });

        if (iset_str == "generic")
            iset = ISET::GENERIC;
        else if (iset_str == "avx")
            iset = ISET::AVX;
        else if (iset_str == "avx2")
            iset = ISET::AVX2;
        else if (iset_str == "avx512")
            iset = ISET::AVX512;
        else
            std::cerr << "Error: unable to parse BAOBZI_ARCH. Valid options are: GENERIC, AVX, AVX2, AVX512\n";
    }
#endif

    return iset;
}

void baobzi_eval(const baobzi_t func, const double *x, double *y) { func->eval(func->obj, x, y); }

void baobzi_eval_multi(const baobzi_t func, const double *x, double *res, int ntrg) {
    func->eval_multi(func->obj, x, res, ntrg);
}

void baobzi_stats(baobzi_t func) {
    if (!func)
        return;
    func->stats(func->obj);
}

baobzi_t baobzi_free(baobzi_t func) {
    if (!func)
        return nullptr;
    if (func->obj)
        func->free(func->obj);
    free(func);
    return nullptr;
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

    baobzi_t res;
    try {
        res = (baobzi_t)malloc(sizeof(baobzi_struct));
        res->DIM = input->dim;
        res->ORDER = input->order;
        res->OUTPUT_DIM = input->output_dim;

        int iset = get_iset();

        switch (BAOBZI_JOIN(res->DIM, res->ORDER, iset)) {
#include "baobzi/baobzi_cases.h"
        default: {
            std::cerr << "Baobzi error: Unable to initialize Baobzi function with variables (DIM, ORDER): (" << res->DIM
                      << ", " << res->ORDER << ")\n";

            free(res);
            return nullptr;
        }
        }
    } catch (std::exception &e) {
        std::cerr << e.what() << std::endl;
        free(res);
        return nullptr;
    }
    return res;
}
}
