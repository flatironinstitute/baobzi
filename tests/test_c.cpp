#include <baobzi.h>

#include <catch2/catch_test_macros.hpp>
#include <cmath>

void testfun_1d1(const double *x, double *res, const void *data) {
    const double scale_factor = *(double *)data;
    *res = scale_factor * log(x[0]);
}

void testfun_1d2(const double *x, double *res, const void *data) {
    const double scale_factor = *(double *)data;
    res[0] = scale_factor * log(x[0]);
    res[1] = sin(x[0]);
}

void testfun_2d3(const double *x, double *res, const void *data) {
    const double scale_factor = *(double *)data;
    *res = scale_factor * exp(cos(5.0 * x[0]) * sin(5.0 * x[1]));
}

void testfun_3d1(const double *x, double *res, const void *data) {
    const double scale_factor = *(double *)data;
    *res = scale_factor * exp(cos(5.0 * x[0]) * sin(5.0 * x[1]) * cos(4.0 * x[2]));
}

TEST_CASE("1D1 evaluations", "[baobzi]") {
    baobzi_input_t input = baobzi_input_default;
    const double scale_factor = 1.5;
    input.dim = 1;
    input.output_dim = 1;
    input.order = 8;
    input.tol = 1E-10;
    input.func = testfun_1d1;
    input.data = (void *)&scale_factor;
    const double half_l[] = {1.0};
    const double center[] = {3.0};

    baobzi_t baobzi_func = baobzi_init(&input, center, half_l);
    double y_appx, y_exact;

    SECTION("evaluations at lower left") {
        double x[] = {center[0] - half_l[0]};
        baobzi_eval(baobzi_func, x, &y_appx);
        testfun_1d1(x, &y_exact, input.data);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < input.tol);
    }

    SECTION("evaluations at center") {
        baobzi_eval(baobzi_func, center, &y_appx);
        testfun_1d1(center, &y_exact, input.data);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < input.tol);
    }

    SECTION("save/restore") {
        const char *filename = "test_func_approx_1d.baobzi";
        baobzi_save(baobzi_func, filename);
        baobzi_t baobzi_func_restored = baobzi_restore(filename);

        double y_old, y_restored;
        baobzi_eval(baobzi_func, center, &y_old);
        baobzi_eval(baobzi_func, center, &y_restored);
        REQUIRE(y_old == y_restored);

        baobzi_free(baobzi_func_restored);
    }

    baobzi_free(baobzi_func);
}

TEST_CASE("1D2 evaluations", "[baobzi]") {
    baobzi_input_t input = baobzi_input_default;
    const double scale_factor = 1.5;
    input.dim = 1;
    input.output_dim = 2;
    input.order = 8;
    input.tol = 1E-10;
    input.func = testfun_1d2;
    input.data = (void *)&scale_factor;
    const double half_l[] = {1.0};
    const double center[] = {3.0};

    baobzi_t baobzi_func = baobzi_init(&input, center, half_l);
    double y_appx[2], y_exact[2];

    SECTION("evaluations at lower left") {
        double x[] = {center[0] - half_l[0]};
        baobzi_eval(baobzi_func, x, y_appx);
        testfun_1d2(x, y_exact, input.data);

        REQUIRE(fabs((y_appx[0] - y_exact[0]) / y_exact[0]) < input.tol);
        REQUIRE(fabs((y_appx[1] - y_exact[1]) / y_exact[1]) < input.tol);
    }

    SECTION("evaluations at center") {
        baobzi_eval(baobzi_func, center, y_appx);
        testfun_1d2(center, y_exact, input.data);

        REQUIRE(fabs((y_appx[0] - y_exact[0]) / y_exact[0]) < input.tol);
        REQUIRE(fabs((y_appx[1] - y_exact[1]) / y_exact[1]) < input.tol);
    }

    SECTION("save/restore") {
        const char *filename = "test_func_approx_1d.baobzi";
        baobzi_save(baobzi_func, filename);
        baobzi_t baobzi_func_restored = baobzi_restore(filename);

        double y_old[2], y_restored[2];
        baobzi_eval(baobzi_func, center, y_old);
        baobzi_eval(baobzi_func, center, y_restored);
        REQUIRE(y_old[0] == y_restored[0]);
        REQUIRE(y_old[1] == y_restored[1]);

        baobzi_free(baobzi_func_restored);
    }

    baobzi_free(baobzi_func);
}

TEST_CASE("2D evaluations", "[baobzi]") {
    baobzi_input_t input = baobzi_input_default;
    const double scale_factor = 1.5;
    input.dim = 2;
    input.output_dim = 1;
    input.order = 6;
    input.tol = 1E-10;
    input.func = testfun_2d3;
    input.data = (void *)&scale_factor;
    const double half_l[2] = {M_PI / 5, M_PI / 5};
    const double center[2] = {-10.0, 3.0};

    baobzi_t baobzi_func = baobzi_init(&input, center, half_l);
    double y_appx, y_exact;

    SECTION("evaluations at lower left") {
        double x[2] = {center[0] - half_l[0], center[1] - half_l[1]};
        baobzi_eval(baobzi_func, x, &y_appx);
        testfun_2d3(x, &y_exact, input.data);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < input.tol);
    }

    SECTION("evaluations at center") {
        baobzi_eval(baobzi_func, center, &y_appx);
        testfun_2d3(center, &y_exact, input.data);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < input.tol);
    }

    SECTION("save/restore") {
        const char *filename = "test_func_approx_2d.baobzi";
        baobzi_save(baobzi_func, filename);
        baobzi_t baobzi_func_restored = baobzi_restore(filename);

        double y_old, y_restored;
        baobzi_eval(baobzi_func, center, &y_old);
        baobzi_eval(baobzi_func, center, &y_restored);
        REQUIRE(y_old == y_restored);

        baobzi_free(baobzi_func_restored);
    }

    baobzi_free(baobzi_func);
}

TEST_CASE("3D evaluations", "[baobzi]") {
    baobzi_input_t input = baobzi_input_default;
    const double scale_factor = 1.5;
    input.dim = 3;
    input.output_dim = 1;
    input.order = 12;
    input.tol = 1E-10;
    input.func = testfun_3d1;
    input.data = (void *)&scale_factor;
    const double half_l[] = {M_PI / 5, M_PI / 5, M_PI / 10.0};
    const double center[] = {-10.0, 3.0, 5.2};

    baobzi_t baobzi_func = baobzi_init(&input, center, half_l);
    double y_appx, y_exact;

    SECTION("evaluations at lower left") {
        double x[] = {center[0] - half_l[0], center[1] - half_l[1], center[2] - half_l[2]};
        baobzi_eval(baobzi_func, x, &y_appx);
        testfun_3d1(x, &y_exact, input.data);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < input.tol);
    }

    SECTION("evaluations at center") {
        baobzi_eval(baobzi_func, center, &y_appx);
        testfun_3d1(center, &y_exact, input.data);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < input.tol);
    }

    SECTION("save/restore") {
        const char *filename = "test_func_approx_3d.baobzi";
        baobzi_save(baobzi_func, filename);
        baobzi_t baobzi_func_restored = baobzi_restore(filename);

        double y_old, y_restored;
        baobzi_eval(baobzi_func, center, &y_old);
        baobzi_eval(baobzi_func, center, &y_restored);
        REQUIRE(y_old == y_restored);

        baobzi_free(baobzi_func_restored);
    }

    baobzi_free(baobzi_func);
}
