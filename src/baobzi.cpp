#include <msgpack.hpp>

#include "baobzi.h"
#include "baobzi.hpp"

#include <tuple>

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>

const baobzi_input_t baobzi_input_default = {
    .func = NULL, .data = NULL, .dim = 0, .order = 0, .tol = 0.0, .minimum_leaf_fraction = 0.0, .split_multi_eval = 1};

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
#ifdef __x86_64__
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
            std::cout << "Error: unable to parse BAOBZI_ARCH. Valid options are: GENERIC, AVX, AVX2, AVX512\n";
    }
#endif

    return iset;
}

double baobzi_eval(const baobzi_t func, const double *x) { return func->eval(func->obj, x); }

void baobzi_eval_multi(const baobzi_t func, const double *x, double *res, int ntrg) {
    func->eval_multi(func->obj, x, res, ntrg);
}

void baobzi_save(const baobzi_t func, const char *filename) { func->save(func->obj, filename); }

baobzi_header_t read_header(const char *addr, const std::size_t buflen, std::size_t *offset) {
    msgpack::object_handle oh;
    msgpack::unpack(oh, addr, buflen, *offset); // actually increments offset
    return oh.get().as<baobzi_header_t>();
}

baobzi_t baobzi_restore(const char *filename_cstr) {
    std::string filename(filename_cstr);
    baobzi_t res = (baobzi_t)malloc(sizeof(baobzi_struct));

    std::size_t offset = 0;
    std::string filedata_str = file_to_string(filename);
    baobzi_header_t header = read_header(filedata_str.data(), filedata_str.size(), &offset);

    msgpack::object_handle oh;
    msgpack::unpack(oh, filedata_str.data(), filedata_str.size(), offset);
    msgpack::object obj = oh.get();

    res->DIM = header.dim;
    res->ORDER = header.order;

    auto [dim, order, version] = std::make_tuple(header.dim, header.order, header.version);

    if (version != BAOBZI_HEADER_VERSION) {
        free(res);
        return nullptr;
    }

    int iset = get_iset();
    switch (BAOBZI_JOIN(header.dim, header.order, iset)) {
#include "baobzi/baobzi_cases_restore.h"
    default: {
        std::cerr << "BAOBZI ERROR: Unable to initialize Baobzi function with variables (DIM, ORDER): (" << dim << ", "
                  << order << ")\n";
        break;
    }
    }

    return res;
}

void baobzi_stats(baobzi_t func) {
    if (!func)
        return;
    func->stats(func->obj);
}

baobzi_t baobzi_free(baobzi_t func) {
    if (!func)
        return nullptr;
    func->free(func->obj);
    free(func);
    return nullptr;
}

bool is_valid_func(const baobzi_input_t *input, const double *point) {
    try {
        input->func(point, input->data);
    } catch (std::exception(e)) {
        return false;
    }

    return true;
}

baobzi_t baobzi_init(const baobzi_input_t *input, const double *center, const double *half_length) {
    if (input->tol <= 0.0) {
        std::cerr << "BAOBZI ERROR: Unable to initialize Baobzi due to invalid 'tol' parameter. Please supply "
                     "something greater than zero.\n";
        return nullptr;
    } else if (!input->func || !is_valid_func(input, center)) {
        std::cerr
            << "BAOBZI ERROR: Unable to initialize Baobzi due to empty or invalid 'func' parameter. Please supply "
               "a valid function to fit.\n";
        return nullptr;
    }

    baobzi_t res = (baobzi_t)malloc(sizeof(baobzi_struct));
    res->DIM = input->dim;
    res->ORDER = input->order;

    int iset = get_iset();

    switch (BAOBZI_JOIN(res->DIM, res->ORDER, iset)) {
#include "baobzi/baobzi_cases.h"
    default: {
        std::cerr << "BAOBZI ERROR: Unable to initialize Baobzi function with variables (DIM, ORDER): (" << res->DIM
                  << ", " << res->ORDER << ")\n";
        free(res);
        return nullptr;
    }
    }

    return res;
}
}
