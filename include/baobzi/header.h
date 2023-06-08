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

typedef struct {
    int dim;
    int order;
    int version;
} baobzi_header_t;

#endif
