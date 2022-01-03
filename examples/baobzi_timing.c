#include <baobzi.h>

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

double testfun_1d(const double *x) { return log(x[0]); }
double testfun_2d(const double *x) { return exp(cos(5.0 * x[0]) * sin(5.0 * x[1])); }
double testfun_2d_2(const double *x) { return exp(x[0] + 2 * sin(x[1])) * (x[0] * x[0] + log(2 + x[1])); }
double testfun_3d(const double *x) { return exp(x[0] + 2 * sin(x[1])) * (x[0] * x[0] + log(2 + x[1] * x[2])); }

void time_function(const baobzi_t *function, const double *x, int size, int n_runs) {
    const double time = omp_get_wtime();
    double res = 0.0;
    for (int i_run = 0; i_run < n_runs; ++i_run) {
        for (int i = 0; i < size; i += function->DIM) {
            res += baobzi_eval(function, x);
        }
    }
    const double dt = omp_get_wtime() - time;
    const long n_eval = n_runs * (size / function->DIM);
    printf("%g %g %g\n", dt, n_eval / (dt * 1E6), res);
}

void print_error(const baobzi_t *function, const double *x, int size) {
    double max_error = 0.0;
    double max_rel_error = 0.0;
    double mean_error = 0.0;
    double mean_rel_error = 0.0;

    size_t n_meas = 0;
    for (int i = 0; i < size; i += function->DIM) {
        const double *point = &x[i];

        double actual = function->f_(point);
        double interp = baobzi_eval(function, point);
        double delta = actual - interp;

        max_error = fmax(max_error, fabs(delta));

        if (fabs(actual) > 1E-15) {
            double rel_error = fabs(interp / actual - 1.0);
            max_rel_error = fmax(max_rel_error, rel_error);
            mean_rel_error += fabs(rel_error);
            n_meas++;
        }

        mean_error += fabs(delta);
    }
    mean_error = mean_error / n_meas;
    mean_rel_error = mean_rel_error / n_meas;

    printf("rel error max, mean: %g %g\n", max_rel_error, mean_rel_error);
    printf("abs error max, mean: %g %g\n", max_error, mean_error);
}

#define JOIN_INTS(A, B) (((A) << 16) | (B))

int main(int argc, char *argv[]) {
    srand(1);
    size_t n_points = 1000000;
    size_t n_runs = 50;

    if (argc == 2)
        n_runs = atoi(argv[1]);

    double *x = (double *) aligned_alloc(32, n_points * 3 * sizeof(double));
    for (size_t i = 0; i < n_points * 3; ++i)
        x[i] = ((double)rand()) / RAND_MAX;

    {
        double hl[2] = {1.0, 1.0};
        double center2d[2] = {hl[0] + 0.5, hl[1] + 2.0};
        double *x_2d_transformed = (double *)aligned_alloc(32, n_points * 2 * sizeof(double));

        for (int i = 0; i < 2 * n_points; i += 2)
            for (int j = 0; j < 2; ++j)
                x_2d_transformed[i + j] = hl[j] * (2.0 * x[i + j] - 1.0) + center2d[j];

        baobzi_t func_approx_2d = baobzi_init(&testfun_2d, 2, 8, center2d, hl, 1E-8);

        time_function(&func_approx_2d, x_2d_transformed, n_points * 2, n_runs);
        print_error(&func_approx_2d, x_2d_transformed, n_points * 2);

        free(x_2d_transformed);
        baobzi_free(&func_approx_2d);
    }

    return 0;
}
