#ifndef BAOBZI_HEADER_H
#define BAOBZI_HEADER_H

#define BAOBZI_HEADER_VERSION 5

#ifdef __cplusplus
#include <msgpack.hpp>
#endif

typedef struct {
    int dim;     ///< Dimension of function
    int order;   ///< Order of polynomial
    int version; ///< Version of output format (BAOBZI_HEADER_VERSION)
#ifdef __cplusplus
    MSGPACK_DEFINE(dim, order, version);
#endif

} baobzi_header_t;

#endif
