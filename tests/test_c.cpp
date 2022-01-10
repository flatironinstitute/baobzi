#include <baobzi.h>

#include <catch2/catch_test_macros.hpp>
#include <cmath>

double testfun_2d(const double *x) { return exp(cos(5.0 * x[0]) * sin(5.0 * x[1])); }

TEST_CASE("2D evaluations", "[baobzi]") {
    const int dim = 2;
    const int order = 6;
    const double tol = 1E-10;
    const double half_l[2] = {M_PI / 5, M_PI / 5};
    const double center[2] = {-10.0, 3.0};

    baobzi_t test_fun = baobzi_init(testfun_2d, dim, order, center, half_l, tol);

    SECTION("evaluations at lower left") {
        double x[2] = {center[0] - half_l[0], center[1] - half_l[1]};
        double y_appx = baobzi_eval(&test_fun, x);
        double y_exact = test_fun.f_(x);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < tol);
    }

    SECTION("evaluations at center") {
        double y_appx = baobzi_eval(&test_fun, center);
        double y_exact = test_fun.f_(center);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < tol);
    }
}
