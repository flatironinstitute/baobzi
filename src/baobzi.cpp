#include "baobzi.hpp"
#include "baobzi.h"

#include <msgpack.hpp>
#include <tuple>
#include <unistd.h>
#define EIGEN_MATRIX_PLUGIN "eigen_matrix_plugin.h"

#include <cstdint>
#include <fstream>
#include <iostream>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>

enum ISET { GENERIC, AVX, AVX2, AVX512 };

struct mmap_wrapper {
    char *addr;
    int fd;
    std::size_t buflen;

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
    return iset;
}

double baobzi_eval(const baobzi_t func, const double *x) { return func->eval(func->obj, x); }

void baobzi_save(const baobzi_t func, const char *filename) { func->save(func->obj, filename); }

baobzi_header_t read_header(const char *addr, std::size_t buflen, std::size_t *offset) {
    msgpack::object_handle oh;
    msgpack::unpack(oh, addr, buflen, *offset);
    return oh.get().as<baobzi_header_t>();
}

baobzi_header_t read_header_from_file(const char *fname) {
    std::string filename(fname);
    mmap_wrapper infile(filename);
    std::size_t buflen = 0, offset = 0;
    msgpack::object_handle oh;
    msgpack::unpack(oh, infile.addr, buflen, offset);
    return oh.get().as<baobzi_header_t>();
}

baobzi_t baobzi_restore(double (*fin)(const double *), const char *filename_cstr) {
    std::string filename(filename_cstr);
    baobzi_t res = (baobzi_t)malloc(sizeof(baobzi_struct));

    std::size_t offset = 0;
    mmap_wrapper infile(filename);

    baobzi_header_t header = read_header(infile.addr, infile.buflen, &offset);

    msgpack::object_handle oh;
    msgpack::unpack(oh, infile.addr, infile.buflen, offset);
    msgpack::object obj = oh.get();

    res->f_ = fin;
    res->DIM = header.dim;
    res->ORDER = header.order;

    auto [dim, order, version] = std::make_tuple(header.dim, header.order, header.version);

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

baobzi_t baobzi_free(baobzi_t func) {
    func->free(func->obj);
    free(func);
    return nullptr;
}

baobzi_t baobzi_init(double (*fin)(const double *), uint16_t dim, uint16_t order, const double *center,
                     const double *half_length, const double tol) {
    baobzi_t res = (baobzi_t)malloc(sizeof(baobzi_struct));
    res->f_ = fin;
    res->DIM = dim;
    res->ORDER = order;

    int iset = get_iset();

    switch (BAOBZI_JOIN(dim, order, iset)) {
#include "baobzi/baobzi_cases.h"
    default: {
        std::cerr << "BAOBZI ERROR: Unable to initialize Baobzi function with variables (DIM, ORDER): (" << dim << ", "
                  << order << ")\n";
        break;
    }
    }

    return res;
}
}
