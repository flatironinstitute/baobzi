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
    void *obj;                                    ///< Actual baobzi::Function object
    int DIM;                                      ///< Dimension of our function
    int ORDER;                                    ///< Order of the polynomial
    double (*eval)(const void *, const double *); ///< Pointer to evaluation function
    void (*save)(const void *, const char *);     ///< Pointer to save function
    void (*free)(void *);                         ///< pointer to save function
} baobzi_struct;

/// Our type for the C API
typedef baobzi_struct *baobzi_t;

extern const baobzi_input_t baobzi_input_default;

/// @brief eval approximator at point x
/// @param[in] func initialized C baobzi object
/// @param[in] x point to evaluate at
/// @returns approximation of function at x
double baobzi_eval(const baobzi_t func, const double *x);

/// @brief save approximator to file
/// @param[in] func initialized C baobzi object
/// @param[in] filename path to output file
void baobzi_save(const baobzi_t func, const char *filename);

/// @brief restore approximator from file
/// @param[in] filename path to serialized baobzi file
/// @returns initialized baobzi C object
baobzi_t baobzi_restore(const char *filename);

/// @brief free all memory associated with C baobzi object
/// @returns nullptr
baobzi_t baobzi_free(baobzi_t func);

/// @brief Construct C baobzi object from input function
/// @param[in] fin pointer to function to fit
/// @param[in] dim input dimension of the function
/// @param[in] order order of polynomial to fit
/// @param[in] center [dim] center of the domain
/// @param[in] half_length [dim] half the size of the domain in each dimension
/// @param[in] tol desired relative tolerance
/// @returns initialized baobzi C object
baobzi_t baobzi_init(const baobzi_input_t *input, const double *center, const double *half_length);

/// @brief read dim/order/version info from file
/// @param[in] fname path to input file
/// @returns file header
baobzi_header_t baobzi_read_header_from_file(const char *fname);

#include "baobzi/baobzi_decls.h"

#ifdef __cplusplus
}
#endif

#endif
