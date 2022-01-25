#include <baobzi_template.hpp>

#include <fstream>
#include <iostream>
#include <omp.h>
#include <random>

using aligned_vector = std::vector<double, Eigen::aligned_allocator<double>>;

double testfun_1d(const double *x, const void *data) { return log(x[0]); }
double testfun_2d(const double *x, const void *data) {
    const double scale_factor = *(double *)data;
    return scale_factor * exp(cos(5.0 * x[0]) * sin(5.0 * x[1]));
}
double testfun_2d_2(const double *x, const void *data) {
    return exp(x[0] + 2 * sin(x[1])) * (x[0] * x[0] + log(2 + x[1]));
}
double testfun_3d(const double *x, const void *data) {
    return exp(x[0] + 2 * sin(x[1])) * (x[0] * x[0] + log(2 + x[1] * x[2]));
}

template <int DIM, typename Function>
void time_function(const Function &function, const aligned_vector &x, int n_runs) {
    const double time = omp_get_wtime();
    double res = 0.0;
    for (int i_run = 0; i_run < n_runs; ++i_run) {
        for (int i = 0; i < x.size(); i += DIM) {
            res += function(&x[i]);
        }
    }
    const double dt = omp_get_wtime() - time;
    const long n_eval = n_runs * (x.size() / DIM);
    std::cout << dt << " " << n_eval / (dt * 1E6) << " " << res << std::endl;
}

template <typename Function>
void print_error(const Function &function, baobzi_input_func_t exact_function, const aligned_vector &x) {
    double max_error = 0.0;
    double max_rel_error = 0.0;
    double mean_error = 0.0;
    double mean_rel_error = 0.0;

    size_t n_meas = 0;
    for (int i = 0; i < x.size(); i += Function::Dim) {
        const double *point = &x[i];

        double actual = exact_function(point, nullptr);
        double interp = function.eval(point);
        double delta = actual - interp;

        max_error = std::max(max_error, std::fabs(delta));

        if (std::abs(actual) > 1E-15) {
            double rel_error = std::abs(interp / actual - 1.0);
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
    aligned_vector x(n_points * 3);
    for (size_t i = 0; i < n_points * 3; ++i)
        x[i] = dis(gen);

    // {
    //     aligned_vector x_1d_transformed(n_points);
    //     Eigen::Vector<double, 1> hl_1d{2.0};
    //     Eigen::Vector<double, 1> center_1d = hl_1d + Eigen::Vector<double, 1>{2.5};
    //     baobzi::Function<1, 8> func_approx_1d(center_1d, hl_1d, testfun_1d, 1E-8);
    //     for (int i = 0; i < n_points; i += 1)
    //         x_1d_transformed[i] = hl_1d[0] * (2.0 * x[i] - 1.0) + center_1d[0];

    //     time_function<1>(func_approx_1d, x_1d_transformed, n_runs);
    //     print_error(func_approx_1d, x_1d_transformed);
    // }

    {
        Eigen::Vector2d hl{1.0, 1.0};
        Eigen::Vector2d center2d = hl + Eigen::Vector2d{0.5, 2.0};
        aligned_vector x_2d_transformed(n_points * 2);
        double scale_factor = 1.5;
        baobzi_input_t input;
        input.dim = 2;
        input.order = 6;
        input.data = &scale_factor;
        input.tol = 1E-8;
        input.func = testfun_2d_2;

        for (int i = 0; i < 2 * n_points; i += 2)
            for (int j = 0; j < 2; ++j)
                x_2d_transformed[i + j] = hl[j] * (2.0 * x[i + j] - 1.0) + center2d[j];

        baobzi::Function<2, 6> func_approx_2d(&input, center2d.data(), hl.data());

        // time_function<2>(func_approx_2d.f_, x_2d_transformed, 1);
        time_function<2>(func_approx_2d, x_2d_transformed, n_runs);
        print_error(func_approx_2d, input.func, x_2d_transformed);
    }

    return 0;
}
