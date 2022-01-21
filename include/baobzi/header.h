#ifndef BAOBZI_HEADER_H
#define BAOBZI_HEADER_H

#define BAOBZI_HEADER_VERSION 0

#ifdef BAOBZI_TEMPLATE_HPP
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
    int dim; ///< Dimension of function
    int order; ///< Order of polynomial
    int version; ///< Version of output format (BAOBZI_HEADER_VERSION)
} baobzi_header_t;
#endif


#endif
