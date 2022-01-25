#include <baobzi.h>

#include <catch2/catch_test_macros.hpp>
#include <cmath>

double testfun_2d(const double *x, const void *data) { return exp(cos(5.0 * x[0]) * sin(5.0 * x[1])); }

TEST_CASE("2D evaluations", "[baobzi]") {
    baobzi_input_t input;
    input.dim = 2;
    input.order = 6;
    input.tol = 1E-10;
    input.func = testfun_2d;
    const double half_l[2] = {M_PI / 5, M_PI / 5};
    const double center[2] = {-10.0, 3.0};
    const void *data = nullptr;

    baobzi_t baobzi_func = baobzi_init(&input, center, half_l);

    SECTION("evaluations at lower left") {
        double x[2] = {center[0] - half_l[0], center[1] - half_l[1]};
        double y_appx = baobzi_eval(baobzi_func, x);
        double y_exact = testfun_2d(x, data);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < input.tol);
    }

    SECTION("evaluations at center") {
        double y_appx = baobzi_eval(baobzi_func, center);
        double y_exact = testfun_2d(center, data);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < input.tol);
    }

    SECTION("save/restore") {
        const char *filename = "test_func_approx_2d.baobzi";
        baobzi_save(baobzi_func, filename);
        baobzi_t baobzi_func_restored = baobzi_restore(filename);

        REQUIRE(baobzi_eval(baobzi_func, center) == baobzi_eval(baobzi_func_restored, center));

        baobzi_free(baobzi_func_restored);
    }

    baobzi_free(baobzi_func);
}
