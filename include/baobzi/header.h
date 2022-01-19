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
typedef struct {
    int dim;
    int order;
    int version;
} baobzi_header_t;
#endif


#endif
