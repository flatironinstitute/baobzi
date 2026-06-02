#ifndef BAOBZI_H
#define BAOBZI_H

#define BAOBZI_DEFAULT(x)
#ifdef __cplusplus
#undef BAOBZI_DEFAULT
#define BAOBZI_DEFAULT(x) = x
#endif

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
    baobzi_input_func_t func BAOBZI_DEFAULT(nullptr);
    void *data BAOBZI_DEFAULT(nullptr);
    int input_dim BAOBZI_DEFAULT(0);
    int output_dim BAOBZI_DEFAULT(1);
    int degree BAOBZI_DEFAULT(8);
    double tol BAOBZI_DEFAULT(0.0);
    double minimum_leaf_fraction BAOBZI_DEFAULT(0.0);
    int split_multi_eval BAOBZI_DEFAULT(1);
    int min_depth BAOBZI_DEFAULT(0);
    int max_depth BAOBZI_DEFAULT(50);
    int tol_type BAOBZI_DEFAULT(BAOBZI_TOL_RELATIVE_MAX);
    int n_samples_per_dim BAOBZI_DEFAULT(8);
};

typedef struct baobzi_input_t baobzi_input_t;
typedef struct baobzi_function baobzi_function;
typedef baobzi_function *baobzi_t;

extern const baobzi_input_t baobzi_input_default;

/// @brief Construct C baobzi object from input function
/// @param[in] input pointer to baobzi_input_t object
/// @param[in] center [DIM] center of the domain
/// @param[in] half_length [DIM] half the size of the domain in each dimension
/// @returns initialized baobzi C object
baobzi_t baobzi_init(const baobzi_input_t *input, const double *center, const double *half_length);

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
void baobzi_stats(const baobzi_t func);

/// @brief free all memory associated with C baobzi object
/// @returns nullptr
baobzi_t baobzi_free(baobzi_t func);

#ifdef __cplusplus
}
#endif

#undef BAOBZI_DEFAULT

#endif
