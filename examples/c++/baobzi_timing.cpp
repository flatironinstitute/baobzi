#include <baobzi_template.hpp>

#include <fstream>
#include <iostream>
#include <random>
#include <time.h>

using real_t = double;

struct timespec get_wtime() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts;
}

double get_wtime_diff(const struct timespec *ts, const struct timespec *tf) {
    return (tf->tv_sec - ts->tv_sec) + (tf->tv_nsec - ts->tv_nsec) * 1E-9;
}

void testfun_1d(const double *x, double *y, const void *data) {
    const double scale_factor = *(real_t *)data;
    *y = scale_factor * log(x[0]);
}
void testfun_2d3(const double *x, double *y, const void *data) {
    const double scale_factor = *(double *)data;
    *y = scale_factor * exp(cos(5.0 * x[0]) * sin(5.0 * x[1]));
}
void testfun_2d_2(const double *x, double *y, const void *data) {
    *y = exp(x[0] + 2 * sin(x[1])) * (x[0] * x[0] + log(2 + x[1]));
}
void testfun_3d1(const double *x, double *y, const void *data) {
    *y = exp(x[0] + 2 * sin(x[1])) * (x[0] * x[0] + log(2 + x[1] * x[2]));
}

template <int DIM, typename Function>
std::vector<real_t> time_function(const Function &function, const std::vector<real_t> &x, int n_runs) {
    const size_t n_points = x.size() / DIM;
    std::vector<real_t> res(n_points);
    const auto st = get_wtime();
    for (int i_run = 0; i_run < n_runs; ++i_run)
        function(x.data(), res.data(), n_points);
    const auto ft = get_wtime();

    const real_t dt = get_wtime_diff(&st, &ft);
    const long n_eval = n_runs * n_points;
    std::cout << "Elapsed time (s): " << dt << std::endl;
    std::cout << "Mevals/s: " << n_eval / (dt * 1E6) << std::endl;
    return res;
}

template <typename Function>
void print_error(const Function &function, baobzi_input_t &input, const std::vector<real_t> &x) {
    real_t max_error = 0.0;
    real_t max_rel_error = 0.0;
    real_t mean_error = 0.0;
    real_t mean_rel_error = 0.0;

    size_t n_meas = 0;
    for (int i = 0; i < x.size(); i += Function::Dim) {
        const real_t *point = &x[i];
        double pointd[Function::Dim];
        for (int j = 0; j < Function::Dim; ++j)
            pointd[j] = x[i + j];

        real_t actual;
        input.func(pointd, &actual, input.data);
        real_t interp = function.eval(point);
        real_t delta = actual - interp;

        max_error = std::max(max_error, std::fabs(delta));

        if (std::abs(actual) > 1E-15) {
            real_t rel_error = std::abs(interp / actual - 1.0);
            max_rel_error = std::max(max_rel_error, rel_error);
            mean_rel_error += std::abs(rel_error);
            n_meas++;
        }

        mean_error += std::abs(delta);
    }
    mean_error = mean_error / n_meas;
    mean_rel_error = mean_rel_error / n_meas;

    std::cout << "rel error max, mean: " << max_rel_error << " " << mean_rel_error << std::endl;
    std::cout << "abs error max, mean: " << max_error << " " << mean_error << std::endl;
}

int main(int argc, char *argv[]) {
    size_t n_points = 1000000;
    size_t n_runs = 50;

    if (argc == 2)
        n_runs = atoi(argv[1]);

    std::mt19937 gen(1);
    std::uniform_real_distribution<> dis(0, 1);
    std::vector<real_t> x(n_points * 3);
    for (size_t i = 0; i < n_points * 3; ++i)
        x[i] = dis(gen);

    {
        real_t hl = 1.0;
        real_t center = 2.0;
        std::vector<real_t> x_transformed(n_points);
        real_t scale_factor = 1.5;
        baobzi_input_t input;
        input.dim = 1;
        input.order = 8;
        input.data = &scale_factor;
        input.tol = 1E-10;
        input.func = testfun_1d;
        input.minimum_leaf_fraction = 1.0;
        input.split_multi_eval = 0;
        input.max_depth = 8;
        input.output_dim = 1;

        for (int i = 0; i < n_points; i++)
            x_transformed[i] = hl * (2.0 * x[i] - 1.0) + center;

        std::cout << "Testing on 1D function...\n";
        baobzi::Function<1, 6, 0, real_t> func_approx_1d(&input, &center, &hl);
        func_approx_1d.print_stats();

        time_function<1>(func_approx_1d, x_transformed, n_runs);
        print_error(func_approx_1d, input, x_transformed);
        std::cout << "\n";
    }

    {
        Eigen::Vector<real_t, 2> hl{1.0, 1.0};
        Eigen::Vector<real_t, 2> center2d = hl + Eigen::Vector<real_t, 2>{0.5, 2.0};
        std::vector<real_t> x_2d_transformed(n_points * 2);
        real_t scale_factor = 1.5;
        baobzi_input_t input;
        input.dim = 2;
        input.order = 10;
        input.data = &scale_factor;
        input.tol = 1E-10;
        input.func = testfun_2d3;
        input.minimum_leaf_fraction = 0.0;
        input.split_multi_eval = 1;
        input.max_depth = 50;
        input.output_dim = 1;

        for (int i = 0; i < 2 * n_points; i += 2)
            for (int j = 0; j < 2; ++j)
                x_2d_transformed[i + j] = hl[j] * (2.0 * x[i + j] - 1.0) + center2d[j];

        std::cout << "Testing on 2D function...\n";
        baobzi::Function<2, 10, 0, real_t> func_approx_2d(&input, center2d.data(), hl.data());
        func_approx_2d.print_stats();

        time_function<2>(func_approx_2d, x_2d_transformed, n_runs);
        print_error(func_approx_2d, input, x_2d_transformed);
    }

    return 0;
}
