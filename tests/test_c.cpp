#include <baobzi.h>

#include <catch2/catch_test_macros.hpp>
#include <cmath>

double testfun_2d(const double *x) { return exp(cos(5.0 * x[0]) * sin(5.0 * x[1])); }
double testfun_4d(const double *x) { return x[0] * x[1] * x[2] * x[3]; }

TEST_CASE("2D evaluations", "[baobzi]") {
    const int dim = 2;
    const int order = 6;
    const double tol = 1E-10;
    const double half_l[2] = {M_PI / 5, M_PI / 5};
    const double center[2] = {-10.0, 3.0};

    baobzi_t baobzi_func = baobzi_init(testfun_2d, dim, order, center, half_l, tol);

    SECTION("evaluations at lower left") {
        double x[2] = {center[0] - half_l[0], center[1] - half_l[1]};
        double y_appx = baobzi_eval(&baobzi_func, x);
        double y_exact = baobzi_func.f_(x);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < tol);
    }

    SECTION("evaluations at center") {
        double y_appx = baobzi_eval(&baobzi_func, center);
        double y_exact = baobzi_func.f_(center);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < tol);
    }

    SECTION("save/restore") {
        const char *filename = "test_func_approx_2d.baobzi";
        baobzi_save(&baobzi_func, filename);
        baobzi_t baobzi_func_restored = baobzi_restore(baobzi_func.f_, filename);

        REQUIRE(baobzi_eval(&baobzi_func, center) == baobzi_eval(&baobzi_func_restored, center));

        baobzi_free(&baobzi_func_restored);
    }

    baobzi_free(&baobzi_func);
}

TEST_CASE("4D evaluations", "[baobzi]") {
    const int dim = 4;
    const int order = 6;
    const double tol = 1E-6;
    const double half_l[4] = {2.0, 2.0, 2.0, 2.0};
    const double center[4] = {0.0, 0.0, 0.0, 0.0};

    baobzi_t baobzi_func = baobzi_init(testfun_4d, dim, order, center, half_l, tol);

    baobzi_eval(&baobzi_func, center);

    baobzi_free(&baobzi_func);
}
