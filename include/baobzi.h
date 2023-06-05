#ifndef BAOBZI_H
#define BAOBZI_H

#include "baobzi/header.h"
#include "baobzi/macros.h"

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

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
    void (*save)(const void *, const char *);                             ///< Pointer to save function
    void (*stats)(void *);                                                ///< pointer to stats function
    void (*free)(void *);                                                 ///< pointer to free function
} baobzi_struct;

/// Our type for the C API
typedef baobzi_struct *baobzi_t;

extern const baobzi_input_t baobzi_input_default;

/// @brief eval approximator at point x
/// @param[in] func initialized C baobzi object
/// @param[in] x point to evaluate at
/// @returns approximation of function at x

void baobzi_eval(const baobzi_t func, const double *x, double *y);

/// @brief eval function approximation at ntrg points
/// @param[in] func initialized C baobzi object
/// /// @param[in] xp [DIM * ntrg] array of points to evaluate function at
/// @param[out] res [DIM * ntrg] array of results
void baobzi_eval_multi(const baobzi_t func, const double *x, double *res, int ntrg);

/// @brief save approximator to file
/// @param[in] func initialized C baobzi object
/// @param[in] filename path to output file
void baobzi_save(const baobzi_t func, const char *filename);

/// @brief restore approximator from file
/// @param[in] filename path to serialized baobzi file
/// @returns initialized baobzi C object
baobzi_t baobzi_restore(const char *filename);

/// @brief Print stats about baobzi object creation
void baobzi_stats(baobzi_t func);

/// @brief free all memory associated with C baobzi object
/// @returns nullptr
baobzi_t baobzi_free(baobzi_t func);

/// @brief Construct C baobzi object from input function
/// @param[in] input pointer to baobzi_input_t object
/// @param[in] center [dim] center of the domain
/// @param[in] half_length [dim] half the size of the domain in each dimension
/// @returns initialized baobzi C object
baobzi_t baobzi_init(const baobzi_input_t *input, const double *center, const double *half_length);

#include "baobzi/baobzi_decls.h"

#ifdef __cplusplus
}
#endif

#endif
