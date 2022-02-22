#include <msgpack.hpp>

#include "baobzi.h"
#include "baobzi.hpp"

#include <exception>
#include <stdexcept>
#include <tuple>
#include <unistd.h>

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

const baobzi_input_t baobzi_input_default = {
    .func = NULL, .data = NULL, .dim = 0, .order = 0, .tol = 0.0, .minimum_leaf_fraction = 0.0};

enum ISET { GENERIC, AVX, AVX2, AVX512 };

/// @brief mmap a file
struct mmap_wrapper {
    char *addr;         ///< Underlying buffer
    int fd;             ///< Underlying file descriptor
    std::size_t buflen; ///< Size of file/buffer

    /// @brief Construct mmap_wrapper object
    /// @param[in] fname path to file to map
    mmap_wrapper(const std::string &fname) {
        fd = open(fname.c_str(), O_RDONLY);
        if (fd == -1)
            throw std::runtime_error("Unable to open baobzi file " + fname + ".");

        struct stat sb;
        if (fstat(fd, &sb) == -1)
            throw std::runtime_error("Error statting " + fname + ".");

        buflen = sb.st_size;

        addr = static_cast<char *>(mmap(NULL, buflen, PROT_READ, MAP_PRIVATE, fd, 0u));
        if (addr == MAP_FAILED)
            throw std::runtime_error("Error mapping " + fname + " for restore.");
    }

    /// @brief close file descriptor on destruction
    ~mmap_wrapper() { close(fd); }
};

extern "C" {
int get_iset() {
    int iset = ISET::GENERIC;
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

    return iset;
}

double baobzi_eval(const baobzi_t func, const double *x) { return func->eval(func->obj, x); }

void baobzi_eval_multi(const baobzi_t func, const double *x, double *res, int ntrg) {
    func->eval_multi(func->obj, x, res, ntrg);
}

void baobzi_save(const baobzi_t func, const char *filename) { func->save(func->obj, filename); }

baobzi_header_t read_header(const char *addr, std::size_t buflen, std::size_t *offset) {
    msgpack::object_handle oh;
    msgpack::unpack(oh, addr, buflen, *offset);
    return oh.get().as<baobzi_header_t>();
}

baobzi_header_t baobzi_read_header_from_file(const char *fname) {
    std::size_t offset = 0;
    std::string filename(fname);
    mmap_wrapper infile(filename);
    return read_header(infile.addr, infile.buflen, &offset);
}

baobzi_t baobzi_restore(const char *filename_cstr) {
    std::string filename(filename_cstr);
    baobzi_t res = (baobzi_t)malloc(sizeof(baobzi_struct));

    std::size_t offset = 0;
    mmap_wrapper infile(filename);

    baobzi_header_t header = read_header(infile.addr, infile.buflen, &offset);

    msgpack::object_handle oh;
    msgpack::unpack(oh, infile.addr, infile.buflen, offset);
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
