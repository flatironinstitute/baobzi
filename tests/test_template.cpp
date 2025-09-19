#include <baobzi.hpp>

#include <catch2/catch_test_macros.hpp>
#include <cmath>

TEST_CASE("1D1 evaluations", "[baobzi_template]") {
    double scale_factor = 1.5;
    auto testfun_1d1 = [scale_factor](const double x) { return scale_factor * log(x); };

    baobzi_input_t input;
    input.input_dim = 1;
    input.tol = 1E-10;
    input.degree = 8;
    input.split_multi_eval = false;

    const double half_l = 1.0;
    const double center = 3.0;

    auto baobzifunc =
        baobzi::make_function<8>(baobzi_input_t{.tol = 1E-10, .split_multi_eval = false}, center, half_l, testfun_1d1);

    double y_appx, y_exact;
    SECTION("evaluations at lower left") {
        const double x = center - half_l;
        y_appx = baobzifunc(x);
        y_exact = testfun_1d1(x);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < input.tol);
    }

    SECTION("evaluations at center") {
        y_appx = baobzifunc(center);
        y_exact = testfun_1d1(center);

        REQUIRE(fabs((y_appx - y_exact) / y_exact) < input.tol);
    }

    // SECTION("left scalar multiply") {
    //     BaobziFunc newfunc = 1.5 * baobzifunc;
    //     double postscaled, prescaled;
    //     baobzifunc(center, &postscaled);
    //     postscaled *= 1.5;

    //     newfunc(center, &prescaled);

    //     REQUIRE(fabs(prescaled - postscaled) < 1E-15);
    // }

    // SECTION("right scalar multiply") {
    //     BaobziFunc newfunc = baobzifunc * 1.5;
    //     double postscaled, prescaled;
    //     baobzifunc(center, &postscaled);
    //     postscaled *= 1.5;

    //     newfunc(center, &prescaled);

    //     REQUIRE(fabs(prescaled - postscaled) < 1E-15);
    // }

    // SECTION("scalar divide") {
    //     BaobziFunc newfunc = baobzifunc / 2.0;
    //     double postscaled, prescaled;
    //     baobzifunc(center, &postscaled);
    //     postscaled /= 2.0;

    //     newfunc(center, &prescaled);

    //     REQUIRE(fabs(prescaled - postscaled) < 1E-15);
    // }

    // SECTION("right scalar add") {
    //     BaobziFunc newfunc = baobzifunc + 2.0;
    //     double post, pre;
    //     baobzifunc(center, &post);
    //     post += 2.0;

    //     newfunc(center, &pre);

    //     REQUIRE(fabs(pre - post) < 1E-15);
    // }

    // SECTION("left scalar add") {
    //     BaobziFunc newfunc = 2.0 + baobzifunc;
    //     double post, pre;
    //     baobzifunc(center, &post);
    //     post += 2.0;

    //     newfunc(center, &pre);

    //     REQUIRE(fabs(pre - post) < 1E-15);
    // }

    // SECTION("scalar subtract right") {
    //     BaobziFunc newfunc = baobzifunc - 1.0;
    //     double post, pre;
    //     baobzifunc(center, &post);
    //     post -= 1.0;

    //     newfunc(center, &pre);

    //     REQUIRE(fabs(pre - post) < 1E-15);
    // }

    // SECTION("add two baobzi functions") {
    //     auto xfunc = [](const double *x, double *y, const void *) { *y = *x; };
    //     BaobziFunc xfit = BaobziFunc(&input, center, half_l, xfunc, {});

    //     BaobziFunc newfunc = xfit + baobzifunc;
    //     double post, pre;
    //     baobzifunc(center, &post);
    //     post += center[0];

    //     newfunc(center, &pre);

    //     REQUIRE(fabs(pre - post) < 1E-15);
    // }

    // SECTION("subtract two baobzi functions") {
    //     auto xfunc = [](const double *x, double *y, const void *) { *y = *x; };
    //     BaobziFunc xfit = BaobziFunc(&input, center, half_l, xfunc, {});

    //     BaobziFunc newfunc = baobzifunc - xfit;
    //     double post, pre;
    //     baobzifunc(center, &post);
    //     post -= center[0];

    //     newfunc(center, &pre);

    //     REQUIRE(fabs(pre - post) < 1E-15);
    // }

    // SECTION("multiply two baobzi functions") {
    //     auto xfunc = [](const double *x, double *y, const void *) { *y = *x; };
    //     BaobziFunc xfit = BaobziFunc(&input, center, half_l, xfunc, {});

    //     BaobziFunc newfunc = baobzifunc * xfit;
    //     double post, pre;
    //     baobzifunc(center, &post);
    //     post *= center[0];

    //     newfunc(center, &pre);

    //     REQUIRE(fabs((pre - post) / post) < 1E-12);
    // }

    // SECTION("divide two baobzi functions") {
    //     auto xfunc = [](const double *x, double *y, const void *) { *y = *x; };
    //     BaobziFunc xfit = BaobziFunc(&input, center, half_l, xfunc, {});

    //     BaobziFunc newfunc = baobzifunc / xfit;
    //     double post, pre;
    //     baobzifunc(center, &post);
    //     post /= center[0];

    //     newfunc(center, &pre);

    //     REQUIRE(fabs((pre - post) / post) < 1E-12);
    // }

    // SECTION("sampling") {
    //     auto xfunc = [](const double *x, double *y, const void *) {
    //         constexpr double sigma2 = 1E-4;
    //         *y = exp(-0.5 * *x * *x / sigma2);
    //     };
    //     double center[] = {0.4};
    //     double half_l[] = {1.0};

    //     BaobziFunc gaussfit = BaobziFunc(&input, center, half_l, xfunc, {});
    //     BaobziFunc gaussfit_sample = BaobziFunc(&input, center, half_l, xfunc, {0.0});

    //     double miss, hit;
    //     double x = 0.0;
    //     gaussfit(&x, &miss);
    //     gaussfit_sample(&x, &hit);

    //     REQUIRE(miss < 1E-15);
    //     REQUIRE(std::fabs(hit - 1.0) <= 1E-15);
    // }
}
