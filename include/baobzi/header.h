#ifndef BAOBZI_HEADER_H
#define BAOBZI_HEADER_H

#define BAOBZI_HEADER_VERSION 5

typedef void (*baobzi_input_func_t)(const double *, double *, const void *);

typedef struct {
    baobzi_input_func_t func;
    void *data;
    int dim;
    int output_dim;
    int order;
    double tol;
    double minimum_leaf_fraction;
    int split_multi_eval;
    int min_depth;
    int max_depth;
} baobzi_input_t;

#ifdef __cplusplus
#include <msgpack.hpp>
struct baobzi_header_t {
    int dim;
    int order;
    int version;
    MSGPACK_DEFINE(dim, order, version);
};
#else
/// @brief header for serialization
typedef struct {
    int dim;     ///< Dimension of function
    int order;   ///< Order of polynomial
    int version; ///< Version of output format (BAOBZI_HEADER_VERSION)
} baobzi_header_t;
#endif

#endif
