#ifndef BAOBZI_H
#define BAOBZI_H

#include "baobzi/macros.h"

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*baobzi_input_func_t)(const double *, double *, const void *);

typedef enum {
    BAOBZI_TOL_RELATIVE_TAIL = 0,
    BAOBZI_TOL_ABSOLUTE_TAIL = 1,
    BAOBZI_TOL_RELATIVE_MAX = 2,
    BAOBZI_TOL_ABSOLUTE_MAX = 3,
    BAOBZI_TOL_RELATIVE_L2 = 4,
    BAOBZI_TOL_ABSOLUTE_L2 = 5,
} baobzi_tol_t;

/// @brief Input data type to define baobzi function
struct baobzi_input_t {
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
    baobzi_tol_t tol_type;
    int n_samples_per_dim;
#ifdef __cplusplus
    baobzi_input_t()
        : func(nullptr), data(nullptr), dim(0), output_dim(1), order(8), tol(0.0), minimum_leaf_fraction(0.0),
          split_multi_eval(1), min_depth(0), max_depth(50), tol_type(BAOBZI_TOL_RELATIVE_MAX),
          n_samples_per_dim(order) {}
#endif
};

typedef struct baobzi_input_t baobzi_input_t;

/// @brief Baobzi C structure for a common API through C bindings. All work is done through the
/// pointer baobzi_t though.
///
/// Contains pointers to the wrappers of the relevant template C++ functions for a
/// dim+order+instruction set
typedef struct {
    void *obj;                                                            ///< Actual baobzi::Function object
    int DIM;                                                              ///< Input dimension of our function
    int OUTPUT_DIM;                                                       ///< Output dimension of our function
    int ORDER;                                                            ///< Order of the polynomial
    void (*eval)(const void *, const double *, double *);                 ///< Pointer to evaluation function
    void (*eval_multi)(const void *, const double *, double *, int ntrg); ///< Pointer to multi-evaluation function
    void (*stats)(void *);                                                ///< pointer to stats function
    void (*free)(void *);                                                 ///< pointer to free function
} baobzi_struct;

/// Our type for the C API
typedef baobzi_struct *baobzi_t;

extern const baobzi_input_t baobzi_input_default;

/// @brief eval approximator at point x
/// @param[in] func initialized C baobzi object
/// @param[in] x point [DIM] to evaluate at
/// @param[out] y point [OUTPUT_DIM] to store result
/// @returns void
void baobzi_eval(const baobzi_t func, const double *x, double *y);

/// @brief eval function approximation at ntrg points
/// @param[in] func initialized C baobzi object
/// @param[in] x [DIM * ntrg] array of points to evaluate function at
/// @param[out] res [DIM * ntrg] array of results
/// @param[in] ntrg number of points to evaluate
/// @returns void
void baobzi_eval_multi(const baobzi_t func, const double *x, double *res, int ntrg);

/// @brief Print stats about baobzi object creation
/// @param[in] func initialized C baobzi object
/// @returns void
void baobzi_stats(baobzi_t func);

/// @brief free all memory associated with C baobzi object
/// @returns nullptr
baobzi_t baobzi_free(baobzi_t func);

/// @brief Construct C baobzi object from input function
/// @param[in] input pointer to baobzi_input_t object
/// @param[in] center [DIM] center of the domain
/// @param[in] half_length [DIM] half the size of the domain in each dimension
/// @returns initialized baobzi C object
baobzi_t baobzi_init(const baobzi_input_t *input, const double *center, const double *half_length);

#include "baobzi/baobzi_decls.h"

#ifdef __cplusplus
}
#endif

#endif
