///----------------------------------------
///     @file hyperbola_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for the focus-curve hyperbola fit
///           (src/core/hyperbola.cpp).
///  @details Synthesises a V-curve from known (p, a, b), feeds it to
///           find_best_hyperbola_fit(), and verifies recovered parameters.
///           Also exercises hfd_calc / steps_to_focus round-trips.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "core/hyperbola.h"

#include <cmath>
#include <vector>

using namespace astap::core;

///----------------------------------------
/// MARK: Helpers
///----------------------------------------

/// @brief Analytic hyperbola: y = a * sqrt(1 + ((x - p) / b)^2).
[[nodiscard]] static double analytic_hfd(double x, double p, double a, double b) {
	const double dx = (x - p) / b;
	return a * std::sqrt(1.0 + dx * dx);
}

/// @brief Build a symmetric focus curve sampled at 9 positions around `p`.
[[nodiscard]] static std::vector<FocusPoint> make_vcurve(double p, double a, double b) {
	const int offsets[] = {-300, -200, -100, -50, 0, 50, 100, 200, 300};
	std::vector<FocusPoint> data;
	data.reserve(sizeof(offsets) / sizeof(offsets[0]));
	for (int off : offsets) {
		const double x = p + off;
		data.push_back({.position = x, .hfd = analytic_hfd(x, p, a, b)});
	}
	return data;
}

///----------------------------------------
/// MARK: Curve-fit recovery
///----------------------------------------

TEST_CASE("find_best_hyperbola_fit recovers synthetic parameters") {
	SUBCASE("nominal p=5000, a=2.0, b=100") {
		constexpr double p = 5000.0;
		constexpr double a = 2.0;
		constexpr double b = 100.0;

		const auto data = make_vcurve(p, a, b);
		const auto fit  = find_best_hyperbola_fit(data);

		CHECK(fit.focus_position == doctest::Approx(p).epsilon(1e-3));
		CHECK(fit.a              == doctest::Approx(a).epsilon(1e-3));
		CHECK(fit.b              == doctest::Approx(b).epsilon(1e-2));
		// Mean relative error on a perfect hyperbola should be < 1%.
		CHECK(fit.mean_error      < 1e-2);
	}

	SUBCASE("displaced p=12345, a=3.5, b=250") {
		constexpr double p = 12345.0;
		constexpr double a = 3.5;
		constexpr double b = 250.0;

		const auto data = make_vcurve(p, a, b);
		const auto fit  = find_best_hyperbola_fit(data);

		CHECK(fit.focus_position == doctest::Approx(p).epsilon(1e-3));
		CHECK(fit.a              == doctest::Approx(a).epsilon(1e-2));
		CHECK(fit.b              == doctest::Approx(b).epsilon(1e-2));
	}
}

///----------------------------------------
/// MARK: Analytic helpers
///----------------------------------------

TEST_CASE("hfd_calc at perfect focus returns a") {
	CHECK(hfd_calc(5000.0, 5000.0, 2.0, 100.0) == doctest::Approx(2.0));
	CHECK(hfd_calc(  42.0,   42.0, 1.0,  50.0) == doctest::Approx(1.0));
}

TEST_CASE("hfd_calc is symmetric about perfect focus") {
	constexpr double p = 5000.0;
	constexpr double a = 2.0;
	constexpr double b = 100.0;
	for (int d : {50, 100, 150, 200, 300}) {
		const double y_minus = hfd_calc(p - d, p, a, b);
		const double y_plus  = hfd_calc(p + d, p, a, b);
		CAPTURE(d);
		CHECK(y_minus == doctest::Approx(y_plus));
	}
}

TEST_CASE("hfd_calc matches closed-form analytic hyperbola") {
	constexpr double p = 7000.0;
	constexpr double a = 2.5;
	constexpr double b = 80.0;
	for (int d : {-300, -200, -100, 0, 100, 200, 300}) {
		const double x = p + d;
		CAPTURE(x);
		CHECK(hfd_calc(x, p, a, b) == doctest::Approx(analytic_hfd(x, p, a, b)));
	}
}

TEST_CASE("steps_to_focus inverts hfd_calc") {
	constexpr double a = 2.0;
	constexpr double b = 100.0;
	// y(x=200) from analytic_hfd(200, 0, 2, 100) = 2*sqrt(5) = 4.472...
	constexpr double distance = 200.0;
	const double y = analytic_hfd(distance, 0.0, a, b);

	// steps_to_focus consumes hfd value and returns the positive branch
	// distance from the hyperbola minimum that produces that hfd.
	CHECK(steps_to_focus(y, a, b) == doctest::Approx(distance).epsilon(1e-6));
}

TEST_CASE("steps_to_focus at minimum HFD is 0") {
	constexpr double a = 3.0;
	constexpr double b = 150.0;
	CHECK(steps_to_focus(a, a, b) == doctest::Approx(0.0).epsilon(1e-10));
}
