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

template <typename Function>
void time_function(const Function &function, const std::vector<double> &x, int n_runs) {
    double time = omp_get_wtime();
    double res = 0.0;
    size_t n_el = x.size() / Function::Dim;
    using VEC = typename Function::VEC;
    const auto &hl = function.nodes_[0].box_.half_length;
    const auto &center = function.nodes_[0].box_.center;
    for (int i_run = 0; i_run < n_runs; ++i_run) {
        for (int i = 0; i < Function::Dim * n_el; i += Function::Dim) {
            const VEC xvec(&x[i]);
            const VEC point = hl.array() * (2.0 * xvec.array() - 1.0) + center.array();
            res += function(point);
        }
    }
    double dt = omp_get_wtime() - time;
    std::cout << dt << " " << n_runs * n_el / dt / 1E6 << " " << res << std::endl;
}

template <typename Function>
void print_error(const Function &function, const std::vector<double> &x) {
    double max_error = 0.0;
    double max_rel_error = 0.0;
    double rms_error = 0.0;
    double rms_rel_error = 0.0;
    using VEC = typename Function::VEC;
    const auto &hl = function.nodes_[0].box_.half_length;
    const auto &center = function.nodes_[0].box_.center;

    size_t n_meas = 0;
    for (int i = 0; i < x.size(); i += Function::Dim) {
        const VEC xvec(&x[i]);
        const VEC point = hl.array() * (2.0 * xvec.array() - 1.0) + center.array();

        double actual = function.f_(point);
        double interp = function.eval(point);
        double delta = actual - interp;

        max_error = std::max(max_error, std::fabs(delta));
        if (std::abs(actual) > 1E-15) {
            double rel_error = std::abs(interp / actual - 1.0);
            max_rel_error = std::max(max_rel_error, rel_error);
            rms_rel_error += rel_error * rel_error;
            n_meas++;
        }

        rms_error += delta * delta;
    }
    rms_error = sqrt(rms_error) / n_meas;
    rms_rel_error = sqrt(rms_rel_error) / n_meas;

    std::cout << max_rel_error << " " << rms_rel_error << std::endl;
    std::cout << max_error << " " << rms_error << std::endl;
}

int main(int argc, char *argv[]) {
    size_t n_el = 1000000;
    size_t n_runs = 50;

    std::random_device rd;
    std::mt19937 gen(1);
    std::uniform_real_distribution<> dis(0, 1);
    std::vector<double> x(n_el * 2);
    for (size_t i = 0; i < n_el * 2; ++i)
        x[i] = dis(gen);

    double time = omp_get_wtime();
    Eigen::Vector2d hl{1.0, 1.0};
    Eigen::Vector2d center2d = hl + Eigen::Vector2d{0.5, 2.0};

    baobzi::Function<2, 8> func_approx_2d(center2d, hl, testfun_2d_2);
    std::cout << omp_get_wtime() - time << std::endl;

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

    {
        time = omp_get_wtime();
        double res = 0.0;
        for (int i_run = 0; i_run < n_runs; ++i_run) {
            for (int i = 0; i < 2 * n_el; i += 2) {
                Eigen::Vector2d point{hl[0] * (2.0 * x[i] - 1.0) + center2d[0],
                                      hl[1] * (2.0 * x[i + 1] - 1.0) + center2d[1]};
                res += testfun_2d_2(point);
            }
        }
        double dt = omp_get_wtime() - time;
        std::cout << dt << " " << n_runs * n_el / dt / 1E6 << " " << res << std::endl;
    }

    time_function(func_approx_2d, x, n_runs);
    print_error(func_approx_2d, x);

    return 0;
}
