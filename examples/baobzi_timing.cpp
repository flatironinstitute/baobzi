#include <baobzi.hpp>

#include <fstream>
#include <iostream>
#include <omp.h>
#include <random>

double testfun_1d(Eigen::Vector<double, 1> x) { return cos(x[0]); }
double testfun_2d(Eigen::Vector<double, 2> x) { return exp(cos(5.0 * x[0]) * sin(5.0 * x[1])); }
double testfun_2d_2(Eigen::Vector<double, 2> x) { return exp(x[0] + 2 * sin(x[1])) * (x[0] * x[0] + log(2 + x[1])); }
double testfun_3d(Eigen::Vector<double, 3> x) {
    return exp(x[0] + 2 * sin(x[1])) * (x[0] * x[0] + log(2 + x[1] * x[2]));
}

template <int DIM, typename Function>
void time_function(const Function &function, const std::vector<double> &x, int n_runs) {
    const double time = omp_get_wtime();
    double res = 0.0;
    using VEC = Eigen::Vector<double, DIM>;
    for (int i_run = 0; i_run < n_runs; ++i_run) {
        for (int i = 0; i < x.size(); i += DIM) {
            const VEC point(&x[i]);
            res += function(point);
        }
    }
    const double dt = omp_get_wtime() - time;
    const long n_eval = n_runs * (x.size() / DIM);
    std::cout << dt << " " << n_eval / (dt * 1E6) << " " << res << std::endl;
}

template <typename Function>
void print_error(const Function &function, const std::vector<double> &x) {
    double max_error = 0.0;
    double max_rel_error = 0.0;
    double mean_error = 0.0;
    double mean_rel_error = 0.0;
    using VEC = typename Function::VEC;

    size_t n_meas = 0;
    for (int i = 0; i < x.size(); i += Function::Dim) {
        const VEC point(&x[i]);

        double actual = function.f_(point);
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

    std::random_device rd;
    std::mt19937 gen(1);
    std::uniform_real_distribution<> dis(0, 1);
    std::vector<double> x(n_points * 3);
    for (size_t i = 0; i < n_points * 3; ++i)
        x[i] = dis(gen);

    double time = omp_get_wtime();
    Eigen::Vector2d hl{1.0, 1.0};
    Eigen::Vector2d center2d = hl + Eigen::Vector2d{0.5, 2.0};
    std::vector<double> x_2d_transformed(n_points * 2);

    for (int i = 0; i < 2 * n_points; i += 2)
        for (int j = 0; j < 2; ++j)
            x_2d_transformed[i + j] = hl[j] * (2.0 * x[i + j] - 1.0) + center2d[j];

    baobzi::Function<2, 8> func_approx_2d(center2d, hl, testfun_2d_2, 1E-8);

    std::ofstream finterp("funinterp.2d");
    std::ofstream fun("fun.2d");

    finterp.precision(17);
    fun.precision(17);
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            Eigen::Vector<double, 2> point{center2d[0] + hl[0] * (i - 50) / 50.1,
                                           center2d[1] + hl[1] * (j - 50) / 50.1};
            auto node2d = func_approx_2d.find_node(point);
            assert(node2d.box_.contains(point));
            finterp << point[0] << " " << point[1] << " " << node2d.eval(point) << std::endl;
            fun << point[0] << " " << point[1] << " " << func_approx_2d.f_(point) << std::endl;
        }
    }

    time_function<2>(func_approx_2d.f_, x_2d_transformed, n_runs);
    time_function<2>(func_approx_2d, x_2d_transformed, n_runs);
    print_error(func_approx_2d, x_2d_transformed);

    return 0;
}
