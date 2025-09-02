#include <algorithm>
#include <baobzi_template.hpp>

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
void testfun_2d(const double *x, double *y, const void *data) {
    const double scale_factor = *(double *)data;
    *y = scale_factor * exp(cos(5.0 * x[0] * x[0]) * sin(5.0 * x[1]));
}
void testfun_2d_2(const double *x, double *y, const void *data) {
    *y = exp(x[0] + 2 * sin(x[1])) * (x[0] * x[0] + log(2 + x[1]));
}
void testfun_3d1(const double *x, double *y, const void *data) {
    *y = exp(x[0] + 2 * sin(x[1])) * (x[0] * x[0] + log(2 + x[1] * x[2]));
}

template <typename Function>
void time_function(const Function &function, const std::vector<real_t> &x, int n_runs) {
    constexpr int DIM = Function::input_dim;
    constexpr int OUT_DIM = Function::output_dim;
    const size_t n_points = x.size() / DIM;
    std::vector<real_t> res_arr(n_points * OUT_DIM);
    real_t *res = res_arr.data();

    const auto st = get_wtime();
    for (int i_run = 0; i_run < n_runs; ++i_run) {
        if constexpr (DIM == 1)
            function(x.data(), res, n_points);
        else {
            for (size_t i = 0; i < n_points; ++i) {
                std::array<double, DIM> point;
                for (int j = 0; j < DIM; ++j)
                    point[j] = x[i * DIM + j];
                auto out = function(point);
                for (int j = 0; j < OUT_DIM; ++j)
                    res[i * OUT_DIM + j] = out[j];
            }
        }
    }

    const auto ft = get_wtime();
    const real_t dt = get_wtime_diff(&st, &ft);
    const long n_eval = n_runs * n_points;
    volatile auto noopt = res;
    std::cout << "Elapsed time (s): " << dt << std::endl;
    std::cout << "Mevals/s: " << n_eval / (dt * 1E6) << std::endl;
}

template <typename Function>
void print_error(const Function &function, baobzi_input_t &input, const std::vector<real_t> &x) {
    real_t max_error = 0.0;
    real_t max_rel_error = 0.0;
    real_t mean_error = 0.0;
    real_t mean_rel_error = 0.0;
    constexpr int input_dim = Function::input_dim;
    constexpr int output_dim = Function::output_dim;
    using in_arr_t = std::array<double, input_dim>;
    using out_arr_t = std::array<double, output_dim>;

    size_t n_meas = 0;
    for (int i = 0; i < x.size(); i += input_dim) {
        const in_arr_t pointd = [&x, i]() {
            in_arr_t p;
            for (int j = 0; j < input_dim; ++j)
                p[j] = x[i + j];
            return p;
        }();

        out_arr_t actual;
        input.func(pointd.data(), actual.data(), input.data);
        out_arr_t interp{[&function, &pointd]() {
            if constexpr (input_dim == 1)
                return function(pointd[0]);
            else
                return function(pointd);
        }()};

        for (int j = 0; j < output_dim; ++j) {
            const double delta = actual[j] - interp[j];
            max_error = std::max(max_error, std::sqrt(delta));
            mean_error += std::fabs(delta);
        }

        for (int j = 0; j < output_dim; ++j) {
            if (std::abs(actual[j]) > 1E-100) {
                real_t rel_error = std::abs(interp[j] / actual[j] - 1.0);
                max_rel_error = std::max(max_rel_error, rel_error);
                mean_rel_error += std::abs(rel_error);
                n_meas++;
            }
        }
    }
    mean_error = mean_error / x.size();
    mean_rel_error = mean_rel_error / n_meas;

    std::cout << "rel error max, mean: " << max_rel_error << " " << mean_rel_error << std::endl;
    std::cout << "abs error max, mean: " << max_error << " " << mean_error << std::endl;
}

template <int DIM>
std::vector<double> transform(const std::vector<double> &x, int n_points, const std::array<double, DIM> &hl,
                              const std::array<double, DIM> &center) {
    std::vector<double> transformed(n_points * DIM);
    for (int i = 0; i < DIM * n_points; i += DIM)
        for (int j = 0; j < DIM; ++j)
            transformed[i + j] = hl[j] * (2.0 * x[i + j] - 1.0) + center[j];

    return transformed;
}

baobzi_input_t create_input(int dim, baobzi_input_func_t func) {
    static real_t scale_factor = 1.5;
    baobzi_input_t input;
    input.dim = dim;
    input.order = 8;
    input.data = &scale_factor;
    input.tol = 1E-10;
    input.func = func;
    input.minimum_leaf_fraction = 0.0;
    input.split_multi_eval = 0;
    input.max_depth = 50;
    input.output_dim = 1;
    input.tol_type = BAOBZI_TOL_RELATIVE;
    return input;
}

template <int DIM>
void test(int n_runs, int n_points, std::vector<double> &x);

template <>
void test<1>(int n_runs, int n_points, std::vector<double> &x) {
    auto input = create_input(1, testfun_1d);
    const real_t hl = 1.0;
    const real_t center = 2.0;
    std::vector<double> x_transformed = transform<1>(x, n_points, {hl}, {center});

    auto func = [input](double x) -> double {
        double y;
        input.func(&x, &y, input.data);
        return y;
    };

    std::cout << "Testing on 1D function...\n";
    baobzi::Function<6, decltype(func)> func_approx(input, center, hl, func);
    func_approx.print_stats();

    time_function(func_approx, x_transformed, n_runs);
    print_error(func_approx, input, x_transformed);
    std::cout << "\n";
}

template <> void test<2>(int n_runs, int n_points, std::vector<double> &x) {
    std::array<real_t, 2> hl{1.0, 1.0};
    // std::array<real_t, 2> center = {hl[0] + 0.5, hl[1] + 2.0};
    std::array<real_t, 2> center = {0.0, 0.0};
    auto input = create_input(2, testfun_2d);
    const auto x_transformed = transform<2>(x, n_points, hl, center);

    auto func = [input](const std::array<double, 2> &x) -> std::array<double, 1> {
        std::array<double, 1> y{0.0};
        input.func(x.data(), y.data(), input.data);
        return y;
    };

    std::cout << "Testing on 2D function...\n";
    baobzi::Function<10, decltype(func)> func_approx(input, center, hl, func);
    func_approx.print_stats();

    time_function(func_approx, x_transformed, n_runs);
    print_error(func_approx, input, x_transformed);
}

int main(int argc, char *argv[]) {
    size_t n_points = 1000000;
    size_t n_runs = 50;

    std::vector<int> run_dims{1, 2};

    if (argc >= 2)
        n_runs = atoi(argv[1]);
    if (argc > 2) {
        run_dims.clear();
        for (int i = 2; i < argc; ++i)
            run_dims.push_back(atoi(argv[i]));
    }

    std::mt19937 gen(1);
    std::uniform_real_distribution<> dis(0, 1);
    std::vector<real_t> x(n_points * 3);
    for (size_t i = 0; i < n_points * 3; ++i)
        x[i] = dis(gen);

    std::array<void (*)(int, int, std::vector<double> &), 2> runners{&test<1>, &test<2>};

    for (auto dim : run_dims) {
        if (dim < 1 || dim > 2) {
            std::cerr << "Only 1D and 2D tests are implemented\n";
            return 1;
        }

        runners[dim - 1](n_runs, n_points, x);
    }

    return 0;
}
