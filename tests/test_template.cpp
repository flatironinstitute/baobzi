#include <baobzi_template.hpp>

#include <catch2/catch_test_macros.hpp>
#include <cmath>

TEST_CASE("1D1 evaluations", "[baobzi_template]") {
    using BaobziFunc = baobzi::Function<1, 8, 0, double>;

    double scale_factor = 1.5;
    auto testfun_1d1 = [scale_factor](const double *x, double *res, const void *data) {
        *res = scale_factor * log(x[0]);
    };

    baobzi_input_t input;
    input.dim = 1;
    input.output_dim = 1;
    input.order = 8;
    input.tol = 1E-10;
    input.max_depth = 50;
    input.min_depth = 0;
    input.minimum_leaf_fraction = 0.0;
    input.split_multi_eval = false;
    input.func = nullptr;
    input.data = nullptr;
    const double half_l[] = {1.0};
    const double center[] = {3.0};

    auto baobzifunc = baobzi::Function<1, 8, 0, double>(&input, center, half_l, testfun_1d1, {});

    double y_appx, y_exact;
    SECTION("evaluations at lower left") {
        double x[] = {center[0] - half_l[0]};
        baobzifunc(x, &y_appx);
        testfun_1d1(x, &y_exact, nullptr);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < input.tol);
    }

    SECTION("evaluations at center") {
        baobzifunc(center, &y_appx);
        testfun_1d1(center, &y_exact, nullptr);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < input.tol);
    }

    SECTION("left scalar multiply") {
        BaobziFunc newfunc = 1.5 * baobzifunc;
        double postscaled, prescaled;
        baobzifunc(center, &postscaled);
        postscaled *= 1.5;

        newfunc(center, &prescaled);

        REQUIRE(fabs(prescaled - postscaled) < 1E-15);
    }

    SECTION("right scalar multiply") {
        BaobziFunc newfunc = baobzifunc * 1.5;
        double postscaled, prescaled;
        baobzifunc(center, &postscaled);
        postscaled *= 1.5;

        newfunc(center, &prescaled);

        REQUIRE(fabs(prescaled - postscaled) < 1E-15);
    }

    SECTION("scalar divide") {
        BaobziFunc newfunc = baobzifunc / 2.0;
        double postscaled, prescaled;
        baobzifunc(center, &postscaled);
        postscaled /= 2.0;

        newfunc(center, &prescaled);

        REQUIRE(fabs(prescaled - postscaled) < 1E-15);
    }

    SECTION("right scalar add") {
        BaobziFunc newfunc = baobzifunc + 2.0;
        double post, pre;
        baobzifunc(center, &post);
        post += 2.0;

        newfunc(center, &pre);

        REQUIRE(fabs(pre - post) < 1E-15);
    }

    SECTION("left scalar add") {
        BaobziFunc newfunc = 2.0 + baobzifunc;
        double post, pre;
        baobzifunc(center, &post);
        post += 2.0;

        newfunc(center, &pre);

        REQUIRE(fabs(pre - post) < 1E-15);
    }

    SECTION("scalar subtract right") {
        BaobziFunc newfunc = baobzifunc - 1.0;
        double post, pre;
        baobzifunc(center, &post);
        post -= 1.0;

        newfunc(center, &pre);

        REQUIRE(fabs(pre - post) < 1E-15);
    }

    SECTION("add two baobzi functions") {
        auto xfunc = [](const double *x, double *y, const void *) { *y = *x; };
        BaobziFunc xfit = BaobziFunc(&input, center, half_l, xfunc, {});

        BaobziFunc newfunc = xfit + baobzifunc;
        double post, pre;
        baobzifunc(center, &post);
        post += center[0];

        newfunc(center, &pre);

        REQUIRE(fabs(pre - post) < 1E-15);
    }

    SECTION("subtract two baobzi functions") {
        auto xfunc = [](const double *x, double *y, const void *) { *y = *x; };
        BaobziFunc xfit = BaobziFunc(&input, center, half_l, xfunc, {});

        BaobziFunc newfunc = baobzifunc - xfit;
        double post, pre;
        baobzifunc(center, &post);
        post -= center[0];

        newfunc(center, &pre);

        REQUIRE(fabs(pre - post) < 1E-15);
    }

    SECTION("multiply two baobzi functions") {
        auto xfunc = [](const double *x, double *y, const void *) { *y = *x; };
        BaobziFunc xfit = BaobziFunc(&input, center, half_l, xfunc, {});

        BaobziFunc newfunc = baobzifunc * xfit;
        double post, pre;
        baobzifunc(center, &post);
        post *= center[0];

        newfunc(center, &pre);

        REQUIRE(fabs((pre - post) / post) < 1E-12);
    }

    SECTION("divide two baobzi functions") {
        auto xfunc = [](const double *x, double *y, const void *) { *y = *x; };
        BaobziFunc xfit = BaobziFunc(&input, center, half_l, xfunc, {});

        BaobziFunc newfunc = baobzifunc / xfit;
        double post, pre;
        baobzifunc(center, &post);
        post /= center[0];

        newfunc(center, &pre);

        REQUIRE(fabs((pre - post) / post) < 1E-12);
    }
}
