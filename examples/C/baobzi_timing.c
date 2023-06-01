#include <baobzi.h>

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

struct timespec get_wtime() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts;
}

double get_wtime_diff(const struct timespec *ts, const struct timespec *tf) {
    return (tf->tv_sec - ts->tv_sec) + (tf->tv_nsec - ts->tv_nsec) * 1E-9;
}

void testfun_1d(const double *x, double *res, const void *data) {
    const double scale_factor = *(double *)data;
    *res = scale_factor * log(x[0]);
}

void testfun_1d2(const double *x, double *res, const void *data) {
    const double scale_factor = *(double *)data;
    res[0] = scale_factor * log(x[0]);
    res[1] = sin(x[0]);
}

void testfun_2d(const double *x, double *res, const void *data) {
    const double scale_factor = *(double *)data;
    *res = scale_factor * exp(cos(5.0 * x[0]) * sin(5.0 * x[1]));
}

void testfun_2d_2(const double *x, double *res, const void *data) {
    *res = exp(x[0] + 2 * sin(x[1])) * (x[0] * x[0] + log(2 + x[1]));
}

void testfun_3d(const double *x, double *res, const void *data) {
    *res = exp(x[0] + 2 * sin(x[1])) * (x[0] * x[0] + log(2 + x[1] * x[2]));
}

void time_function(const baobzi_t function, const double *x, int size, int n_runs) {
    const int ntrg = size / function->DIM;
    double *res = (double *)malloc(sizeof(double) * ntrg * function->OUTPUT_DIM);

    const struct timespec st = get_wtime();
    for (int i_run = 0; i_run < n_runs; ++i_run) {
        baobzi_eval_multi(function, x, res, ntrg);
    }
    const struct timespec ft = get_wtime();
    const double dt = get_wtime_diff(&st, &ft);
    const long n_eval = n_runs * ntrg;
    printf("Elapsed time: %g\nMevals/s: %g\n", dt, n_eval / (dt * 1E6));
    free(res);
}

void print_error(const baobzi_t function, baobzi_input_t *input, const double *x, int size) {
    double max_error = 0.0;
    double max_rel_error = 0.0;
    double mean_error = 0.0;
    double mean_rel_error = 0.0;

    size_t n_meas = 0;
    for (int i = 0; i < size; i += function->DIM) {
        const double *point = &x[i];

        double actual[input->output_dim];
        input->func(point, actual, input->data);
        double interp[input->output_dim];
        baobzi_eval(function, point, interp);
        for (int i = 0; i < input->output_dim; ++i) {
            double delta = actual[i] - interp[i];

            max_error = fmax(max_error, fabs(delta));

            if (fabs(actual[i]) > 1E-15) {
                double rel_error = fabs(interp[i] / actual[i] - 1.0);
                max_rel_error = fmax(max_rel_error, rel_error);
                mean_rel_error += fabs(rel_error);
                n_meas++;
            }

            mean_error += fabs(delta);
        }
    }
    mean_error = mean_error / n_meas;
    mean_rel_error = mean_rel_error / n_meas;

    printf("rel error max, mean: %g %g\n", max_rel_error, mean_rel_error);
    printf("abs error max, mean: %g %g\n", max_error, mean_error);
}

void test_func(baobzi_input_t *input, const double *xin, const double *hl, const double *center, int n_points,
               int n_runs) {
    // Scale test points to our domain
    double *x_transformed = (double *)malloc(n_points * input->dim * sizeof(double));
    for (int i = 0; i < input->dim * n_points; i += input->dim)
        for (int j = 0; j < input->dim; ++j)
            x_transformed[i + j] = hl[j] * (2.0 * xin[i + j] - 1.0) + center[j];

    // Create baobzi function approximator. Has pointers to relevant structures inside
    // This may take a while, since it fits the function on init
    baobzi_t func_approx = baobzi_init(input, center, hl);
    baobzi_stats(func_approx);

    char filename[256];
    sprintf(filename, "func_approx_%dd", input->dim);

    time_function(func_approx, x_transformed, n_points * input->dim, n_runs);
    print_error(func_approx, input, x_transformed, n_points * input->dim);
    baobzi_save(func_approx, filename);

    free(x_transformed);
    // DON'T FORGET TO FREE THE OBJECT WHEN YOU ARE TOTALLY DEFINITELY DONE WITH IT.
    // They can be HUGE. Also memory leaks :(
    baobzi_free(func_approx);
}

int main(int argc, char *argv[]) {
    srand(1);
    size_t n_points = (size_t)1E6;
    size_t n_runs = 10;
    int order = 10;

    if (argc >= 2)
        n_runs = atoi(argv[1]);
    if (argc >= 3)
        order = atoi(argv[2]);

    // Generate enough points for up to 5 dimensions (for later!)
    double *x = (double *)malloc(n_points * 5 * sizeof(double));
    for (size_t i = 0; i < n_points * 5; ++i)
        x[i] = ((double)rand()) / RAND_MAX;

    {
        printf("Testing on 1D function...\n");
        baobzi_input_t input = baobzi_input_default;
        double scale_factor = 1.5;
        input.dim = 1;
        input.order = order;
        input.tol = 1E-10;
        input.func = testfun_1d;
        input.data = &scale_factor;
        input.minimum_leaf_fraction = 1.0;
        input.split_multi_eval = 0;
        input.output_dim = 1;

        double hl[] = {1.0};
        double center[] = {2.0};
        test_func(&input, x, hl, center, n_points, n_runs);
        printf("\n\n");
    }

    {
        printf("Testing on 1D2 function...\n");
        baobzi_input_t input = baobzi_input_default;
        double scale_factor = 1.5;
        input.dim = 1;
        input.output_dim = 2;
        input.order = order;
        input.tol = 1E-10;
        input.func = testfun_1d2;
        input.data = &scale_factor;
        input.minimum_leaf_fraction = 1.0;
        input.split_multi_eval = 0;

        const double half_l[] = {1.0};
        const double center[] = {3.0};
        baobzi_init(&input, center, half_l);

        test_func(&input, x, half_l, center, n_points, n_runs);
        printf("\n\n");
    }

    {
        printf("Testing on 2D function...\n");
        baobzi_input_t input = baobzi_input_default;
        double scale_factor = 1.5;
        input.dim = 2;
        input.order = order;
        input.func = testfun_2d;
        input.tol = 1E-10; // Maximum relative error target
        input.data = &scale_factor;
        input.minimum_leaf_fraction = 0.0;
        input.split_multi_eval = 1;
        input.output_dim = 1;

        const double hl[2] = {1.0, 1.0};                     // half the length of the domain in each dimension
        const double center[2] = {hl[0] + 0.5, hl[1] + 2.0}; // center of the domain
        test_func(&input, x, hl, center, n_points, n_runs);
        printf("\n\n");
    }

    {
        printf("Testing on 3D function...\n");
        baobzi_input_t input = baobzi_input_default;
        input.dim = 3;
        input.order = 8;
        input.tol = 1E-12;
        input.func = testfun_3d;
        input.minimum_leaf_fraction = 0.0;
        input.split_multi_eval = 1;
        input.output_dim = 1;

        double hl[3] = {1.0, 1.0, 1.0};
        double center[3] = {hl[0] + 0.5, hl[1] + 2.0, hl[2] + 0.5};
        test_func(&input, x, hl, center, n_points, n_runs);
    }

    return 0;
}
