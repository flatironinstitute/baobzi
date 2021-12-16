#include <baobzi.hpp>

#include <iostream>
#include <fstream>
#include <random>
#include <omp.h>

using baobzi::LinearTree;

typedef double (*cfunc_double)(double, double);

double testfun_1d(Eigen::Vector<double, 1> x) { return cos(x[0]); }
double testfun_2d(Eigen::Vector<double, 2> x) { return exp(cos(5.0 * x[0]) * sin(5.0 * x[1])); }
double testfun_2d_2(Eigen::Vector<double, 2> x) { return exp(x[0] + 2 * sin(x[1])) * (x[0] * x[0] + log(2 + x[1])); }
double testfun_3d(Eigen::Vector<double, 3> x) {
    return exp(x[0] + 2 * sin(x[1])) * (x[0] * x[0] + log(2 + x[1] * x[2]));
}

using LinearTree1d = LinearTree<1, 8>;
using LinearTree2d = LinearTree<2, 8>;
using LinearTree3d = LinearTree<3, 6>;

void time_linear_tree_2d(const LinearTree2d &lineartree2d, const std::vector<double> &x, int n_runs) {
    double time = omp_get_wtime();
    double res = 0.0;
    size_t n_el = x.size() / 2;
    double hl = lineartree2d.nodes_[0].box_.half_length;
    auto center2d = lineartree2d.nodes_[0].box_.center;
    for (int i_run = 0; i_run < n_runs; ++i_run) {
        for (int i = 0; i < 2 * n_el; i += 2) {
            Eigen::Vector2d point{hl * (2.0 * x[i] - 1.0) + center2d[0], hl * (2.0 * x[i + 1] - 1.0) + center2d[1]};
            res += lineartree2d(point);
        }
    }
    double dt = omp_get_wtime() - time;
    std::cout << dt << " " << n_runs * n_el / dt / 1E6 << " " << res << std::endl;
}

int main(int argc, char *argv[]) {
    double time = omp_get_wtime();
    double hl = 1.0;
    //     Eigen::Vector3d center3d{0.5 + hl, 0.5 + 2 * hl, 0.5 + hl};
    //     LinearTree3d tree3d_test(center3d, hl, testfun_3d);

    //     std::cout << tree3d_test.eval(center3d) << " " << testfun_3d(center3d) << std::endl;

    //     size_t n_el = 1000000;
    //     size_t n_runs = 50;

    //     std::random_device rd;
    //     std::mt19937 gen(1);
    //     std::uniform_real_distribution<> dis(0, 1);
    //     std::vector<double> x(n_el * 3);
    //     for (size_t i = 0; i < n_el * 3; ++i)
    //         x[i] = dis(gen);

    //     {
    //         time = omp_get_wtime();
    //         double res = 0.0;
    //         for (int i_run = 0; i_run < n_runs; ++i_run) {
    // #pragma omp parallel for schedule(static) reduction(+:res)
    //             for (int i = 0; i < 3 * n_el; i += 3) {
    //                 Eigen::Vector3d point{
    //                     hl * (2.0 * x[i + 0] - 1.0) + center3d[0],
    //                     hl * (2.0 * x[i + 1] - 1.0) + center3d[1],
    //                     hl * (2.0 * x[i + 2] - 1.0) + center3d[2],
    //                 };
    //                 res += tree3d_test.eval(point);
    //             }
    //         }
    //         double dt = omp_get_wtime() - time;
    //         std::cout << dt << " " << n_runs * n_el / dt / 1E6 << " " << res << std::endl;
    //     }

    //     return 0;

    Eigen::Vector2d center2d{0.5 + hl, 0.5 + 2 * hl};
    LinearTree2d lineartree2d(center2d, hl, testfun_2d_2);
    std::cout << omp_get_wtime() - time << std::endl;

    std::ofstream finterp("funinterp.2d");
    std::ofstream fun("fun.2d");

    finterp.precision(17);
    fun.precision(17);
    for (int i = 0; i < 100; ++i) {
        for (int j = 0; j < 100; ++j) {
            Eigen::Vector<double, 2> point{center2d[0] + hl * (i - 50) / 50.1, center2d[1] + hl * (j - 50) / 50.1};
            auto node2d = lineartree2d.find_node(point);
            assert(node2d.box_.contains(point));
            finterp << point[0] << " " << point[1] << " " << node2d.eval(point) << std::endl;
            fun << point[0] << " " << point[1] << " " << lineartree2d.f_(point) << std::endl;
        }
    }

    size_t n_el = 1000000;
    size_t n_runs = 50;

    double start = omp_get_wtime();

    std::random_device rd;
    std::mt19937 gen(1);
    std::uniform_real_distribution<> dis(0, 1);
    std::vector<double> x(n_el * 2);
    for (size_t i = 0; i < n_el * 2; ++i)
        x[i] = dis(gen);

    {
        time = omp_get_wtime();
        double res = 0.0;
        for (int i_run = 0; i_run < n_runs; ++i_run) {
            for (int i = 0; i < 2 * n_el; i += 2) {
                Eigen::Vector2d point{hl * (2.0 * x[i] - 1.0) + center2d[0], hl * (2.0 * x[i + 1] - 1.0) + center2d[1]};
                res += testfun_2d_2(point);
            }
        }
        double dt = omp_get_wtime() - time;
        std::cout << dt << " " << n_runs * n_el / dt / 1E6 << " " << res << std::endl;
    }

    time_linear_tree_2d(lineartree2d, x, n_runs);

    double max_error = 0.0;
    double max_rel_error = 0.0;
    double rms_error = 0.0;
    double rms_rel_error = 0.0;
    for (int i = 0; i < 2 * n_el; i += 2) {
        Eigen::Vector2d point{hl * (2.0 * x[i] - 1.0) + center2d[0], hl * (2.0 * x[i + 1] - 1.0) + center2d[1]};
        double actual = testfun_2d_2(point);
        double interp = lineartree2d.eval(point);
        double delta = actual - interp;

        max_error = std::max(max_error, std::fabs(delta));
        if (std::abs(actual) > 1E-15) {
            double rel_error = std::abs(interp / actual - 1.0);
            max_rel_error = std::max(max_rel_error, rel_error);
            rms_rel_error += rel_error * rel_error;
        }

        rms_error += delta * delta;
    }
    rms_error = 0.5 * sqrt(rms_error) / n_el;
    rms_rel_error = 0.5 * sqrt(rms_rel_error) / n_el;

    std::cout << max_rel_error << " " << rms_rel_error << std::endl;
    std::cout << max_error << " " << rms_error << std::endl;

    // res = 0.0;
    // time = omp_get_wtime();
    // for (int i_run = 0; i_run < n_runs; ++i_run) {
    //     for (int i = 0; i < 2 * n_el; i += 2) {
    //         Eigen::Vector<double, 2> point{x[i], x[i + 1]};
    //         const auto &node2d = tree2d.find_node_index(point, 10);
    //         res += node2d;
    //     }
    // }
    // dt = omp_get_wtime() - time;
    // std::cout << dt << " " << n_runs * n_el / dt / 1E6 << " " << res << std::endl;

    // MyTree1d tree1d(Eigen::Vector<double, 1>{0.0}, 1.0, testfun_1d);

    // auto point1d = Eigen::Vector<double, 1>{0.4};
    // auto node1d = tree1d.find_node(point1d);

    // std::cout << node1d->eval(point1d) << std::endl;
    // std::cout << testfun_1d(point1d) << std::endl;
    // point1d[0] = 0.00001;
    // std::cout << node1d->eval(point1d) - testfun_1d(point1d) << std::endl;
    // std::cout << testfun_1d(point1d) << std::endl;
    return 0;
}
