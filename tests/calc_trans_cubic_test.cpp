///----------------------------------------
///     @file calc_trans_cubic_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for the cubic bivariate transform fit
///           (src/solving/calc_trans_cubic.cpp).
///  @details Verifies identity / translation / linear / synthetic cubic
///           recoveries to 1e-8, plus the error-path (< 10 matched pairs).
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "solving/calc_trans_cubic.h"

#include <array>
#include <cmath>
#include <numbers>
#include <string>

using namespace astap::solving;

///----------------------------------------
/// MARK: Helpers
///----------------------------------------

/// @brief Tolerance for coefficient recovery (Gauss-Jordan on a 10x10
///        normal-equation matrix for our synthetic inputs sits comfortably
///        inside this band).
static constexpr double kCoeffTol = 1e-8;

/// @brief Build an evenly-spaced square grid of reference stars.
[[nodiscard]] static StarArray make_grid(int side, double spacing) {
	StarArray out;
	out.reserve(static_cast<std::size_t>(side * side));
	for (int iy = 0; iy < side; ++iy) {
		for (int ix = 0; ix < side; ++ix) {
			out.push_back({
				.x = (static_cast<double>(ix) - (side - 1) / 2.0) * spacing,
				.y = (static_cast<double>(iy) - (side - 1) / 2.0) * spacing,
			});
		}
	}
	return out;
}

/// @brief Apply a cubic polynomial defined by `t` to every point in `src`.
[[nodiscard]] static StarArray apply_cubic(const StarArray& src, const TransCoeffs& t) {
	StarArray out;
	out.reserve(src.size());
	for (const auto& p : src) {
		const double x  = p.x;
		const double y  = p.y;
		const double x2 = x * x;
		const double y2 = y * y;
		const double xy = x * y;

		const double xp =
			t.x00 +
			t.x10 * x   + t.x01 * y   +
			t.x20 * x2  + t.x11 * xy  + t.x02 * y2 +
			t.x30 * x2 * x + t.x21 * x2 * y +
			t.x12 * x * y2 + t.x03 * y2 * y;

		const double yp =
			t.y00 +
			t.y10 * x   + t.y01 * y   +
			t.y20 * x2  + t.y11 * xy  + t.y02 * y2 +
			t.y30 * x2 * x + t.y21 * x2 * y +
			t.y12 * x * y2 + t.y03 * y2 * y;

		out.push_back({.x = xp, .y = yp});
	}
	return out;
}

/// @brief Assert recovered coefficients match expected within `kCoeffTol`.
static void check_trans_equal(const TransCoeffs& got, const TransCoeffs& want) {
	CHECK(got.x00 == doctest::Approx(want.x00).epsilon(kCoeffTol));
	CHECK(got.x10 == doctest::Approx(want.x10).epsilon(kCoeffTol));
	CHECK(got.x01 == doctest::Approx(want.x01).epsilon(kCoeffTol));
	CHECK(got.x20 == doctest::Approx(want.x20).epsilon(kCoeffTol));
	CHECK(got.x11 == doctest::Approx(want.x11).epsilon(kCoeffTol));
	CHECK(got.x02 == doctest::Approx(want.x02).epsilon(kCoeffTol));
	CHECK(got.x30 == doctest::Approx(want.x30).epsilon(kCoeffTol));
	CHECK(got.x21 == doctest::Approx(want.x21).epsilon(kCoeffTol));
	CHECK(got.x12 == doctest::Approx(want.x12).epsilon(kCoeffTol));
	CHECK(got.x03 == doctest::Approx(want.x03).epsilon(kCoeffTol));

	CHECK(got.y00 == doctest::Approx(want.y00).epsilon(kCoeffTol));
	CHECK(got.y10 == doctest::Approx(want.y10).epsilon(kCoeffTol));
	CHECK(got.y01 == doctest::Approx(want.y01).epsilon(kCoeffTol));
	CHECK(got.y20 == doctest::Approx(want.y20).epsilon(kCoeffTol));
	CHECK(got.y11 == doctest::Approx(want.y11).epsilon(kCoeffTol));
	CHECK(got.y02 == doctest::Approx(want.y02).epsilon(kCoeffTol));
	CHECK(got.y30 == doctest::Approx(want.y30).epsilon(kCoeffTol));
	CHECK(got.y21 == doctest::Approx(want.y21).epsilon(kCoeffTol));
	CHECK(got.y12 == doctest::Approx(want.y12).epsilon(kCoeffTol));
	CHECK(got.y03 == doctest::Approx(want.y03).epsilon(kCoeffTol));
}

///----------------------------------------
/// MARK: Identity
///----------------------------------------

TEST_CASE("identity fit: stars_ref == stars_dist") {
	const auto stars = make_grid(4, 1.0);   // 16 stars
	REQUIRE(stars.size() >= 10);

	const auto r = calc_trans_cubic(stars, stars);
	REQUIRE(r.has_value());

	TransCoeffs want{};   // zero-initialised
	want.x10 = 1.0;
	want.y01 = 1.0;
	check_trans_equal(*r, want);
}

///----------------------------------------
/// MARK: Pure translation
///----------------------------------------

TEST_CASE("pure translation: dist = ref + (dx, dy)") {
	constexpr double dx = 3.5;
	constexpr double dy = -2.25;

	const auto ref = make_grid(4, 1.0);
	StarArray dist = ref;
	for (auto& p : dist) {
		p.x += dx;
		p.y += dy;
	}

	const auto r = calc_trans_cubic(ref, dist);
	REQUIRE(r.has_value());

	TransCoeffs want{};
	want.x00 = dx;
	want.x10 = 1.0;
	want.y00 = dy;
	want.y01 = 1.0;
	check_trans_equal(*r, want);
}

///----------------------------------------
/// MARK: Uniform linear scale
///----------------------------------------

TEST_CASE("uniform scale: dist = s * ref") {
	constexpr double s = 2.0;

	const auto ref = make_grid(4, 1.0);
	StarArray dist = ref;
	for (auto& p : dist) {
		p.x *= s;
		p.y *= s;
	}

	const auto r = calc_trans_cubic(ref, dist);
	REQUIRE(r.has_value());

	TransCoeffs want{};
	want.x10 = s;
	want.y01 = s;
	check_trans_equal(*r, want);
}

///----------------------------------------
/// MARK: Full synthetic cubic — recover all 20 coefficients
///----------------------------------------

TEST_CASE("full cubic: synthesise with known coefficients, recover them") {
	// A grid of 36 well-spread stars so the 10x10 normal matrix is
	// well-conditioned for the cubic fit.
	const auto ref = make_grid(6, 2.0);
	REQUIRE(ref.size() >= 10);

	// Hand-picked non-degenerate cubic coefficients.
	TransCoeffs t{};
	t.x00 = 0.10;  t.x10 = 1.01; t.x01 = 0.02;
	t.x20 = 1e-4;  t.x11 = 2e-4; t.x02 = -1e-4;
	t.x30 = 5e-6;  t.x21 = 3e-6; t.x12 = -2e-6; t.x03 = 1e-6;

	t.y00 = -0.05; t.y10 = -0.03; t.y01 = 0.98;
	t.y20 = -1e-4; t.y11 = 1e-4;  t.y02 = 3e-4;
	t.y30 = -1e-6; t.y21 = 2e-6;  t.y12 = 4e-6; t.y03 = -5e-6;

	const auto dist = apply_cubic(ref, t);
	REQUIRE(dist.size() == ref.size());

	const auto r = calc_trans_cubic(ref, dist);
	REQUIRE(r.has_value());
	check_trans_equal(*r, t);
}

///----------------------------------------
/// MARK: Error path — fewer than 10 matched pairs
///----------------------------------------

TEST_CASE("error: fewer than 10 pairs returns unexpected") {
	const auto tiny = make_grid(3, 1.0);   // 9 stars — one short of the minimum
	REQUIRE(tiny.size() < 10);

	const auto r = calc_trans_cubic(tiny, tiny);
	REQUIRE_FALSE(r.has_value());
	// The error string should mention the minimum-count constraint somehow.
	// We avoid depending on exact wording: just check the message isn't empty.
	CHECK_FALSE(r.error().empty());
}

///----------------------------------------
/// MARK: Error path — size mismatch
///----------------------------------------

TEST_CASE("error: reference and distorted arrays of different sizes") {
	auto ref  = make_grid(4, 1.0);     // 16
	auto dist = make_grid(4, 1.0);
	dist.pop_back();                    // 15 — mismatched lengths

	const auto r = calc_trans_cubic(ref, dist);
	// Either the implementation rejects with an error, or silently truncates.
	// Accept both; just sanity-check: if it succeeded, the shape should still
	// be a plausible identity (both arrays are identity-like up to the cut).
	if (r.has_value()) {
		// 15 valid points are still >= 10, so a fit is possible.
		CHECK(r->x10 == doctest::Approx(1.0).epsilon(1e-6));
		CHECK(r->y01 == doctest::Approx(1.0).epsilon(1e-6));
	} else {
		CHECK_FALSE(r.error().empty());
	}
}
