///----------------------------------------
///     @file star_align_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for the star-matching primitives
///           (src/solving/star_align.cpp).
///  @details Focus on the pure-math primitives: reset_solution_vectors,
///           solution_str formatting, find_quads / find_quads_xy on
///           synthetic star patterns. find_stars and find_offset_and_rotation
///           need HFD + background measurement plumbing that isn't in
///           scope here.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "solving/star_align.h"

#include <cmath>
#include <string>
#include <vector>

using namespace astap::solving;
using astap::StarList;
using astap::SolutionVector;

///----------------------------------------
/// MARK: Helpers
///----------------------------------------

/// @brief Build a 2xN star list (row 0 = x, row 1 = y). Pascal convention.
[[nodiscard]] static StarList make_starlist(
		const std::vector<std::pair<double, double>>& stars) {
	StarList out(2, std::vector<double>(stars.size(), 0.0));
	for (std::size_t i = 0; i < stars.size(); ++i) {
		out[0][i] = stars[i].first;
		out[1][i] = stars[i].second;
	}
	return out;
}

///----------------------------------------
/// MARK: reset_solution_vectors
///----------------------------------------

TEST_CASE("reset_solution_vectors sets diagonal = factor") {
	// Scramble the current globals, then reset.
	solution_vectorX = {1.5, 2.5, 3.5};
	solution_vectorY = {4.5, 5.5, 6.5};
	solution_cblack  = {7.5, 8.5, 9.5};

	reset_solution_vectors(1.0);

	// Pascal semantics: x' = 1*x + 0*y + 0, y' = 0*x + 1*y + 0.
	CHECK(solution_vectorX[0] == doctest::Approx(1.0));
	CHECK(solution_vectorX[1] == doctest::Approx(0.0));
	CHECK(solution_vectorX[2] == doctest::Approx(0.0));
	CHECK(solution_vectorY[0] == doctest::Approx(0.0));
	CHECK(solution_vectorY[1] == doctest::Approx(1.0));
	CHECK(solution_vectorY[2] == doctest::Approx(0.0));
}

TEST_CASE("reset_solution_vectors with a scale factor") {
	reset_solution_vectors(2.5);
	CHECK(solution_vectorX[0] == doctest::Approx(2.5));
	CHECK(solution_vectorX[1] == doctest::Approx(0.0));
	CHECK(solution_vectorY[1] == doctest::Approx(2.5));
}

TEST_CASE("reset_solution_vectors with 'nullify' sentinel (0.001)") {
	reset_solution_vectors(0.001);
	CHECK(solution_vectorX[0] == doctest::Approx(0.001));
	CHECK(solution_vectorY[1] == doctest::Approx(0.001));
}

///----------------------------------------
/// MARK: solution_str
///----------------------------------------

TEST_CASE("solution_str includes x and y equation summaries") {
	reset_solution_vectors(1.0);  // identity so the numbers are deterministic.
	const std::string s = solution_str();
	CAPTURE(s);

	// Just spot-check that the output mentions "x" and "y" (case-insensitive)
	// and isn't empty.
	CHECK_FALSE(s.empty());
	CHECK((s.find("x") != std::string::npos || s.find("X") != std::string::npos));
	CHECK((s.find("y") != std::string::npos || s.find("Y") != std::string::npos));
}

///----------------------------------------
/// MARK: find_quads — structural checks
///----------------------------------------

TEST_CASE("find_quads requires at least 4 stars") {
	StarList stars = make_starlist({
		{10.0, 10.0},
		{20.0, 20.0},
		{30.0, 30.0},
	});

	StarList quads;
	find_quads(stars, quads);

	// With only 3 stars, no quad can be formed.
	REQUIRE(quads.size() >= 1);
	CHECK(quads[0].empty());
}

TEST_CASE("find_quads on a regular grid produces a non-empty quad table") {
	// 5x5 grid of stars evenly spaced; plenty of quads available.
	std::vector<std::pair<double, double>> s;
	s.reserve(25);
	for (int iy = 0; iy < 5; ++iy) {
		for (int ix = 0; ix < 5; ++ix) {
			s.emplace_back(100.0 + ix * 50.0, 100.0 + iy * 50.0);
		}
	}
	auto stars = make_starlist(s);

	StarList quads;
	find_quads(stars, quads);

	// Output layout: 8 rows × K columns (K = number of quads).
	REQUIRE(quads.size() == 8);
	CHECK(quads[0].size() > 0);

	// Each row must have the same length.
	const auto k = quads[0].size();
	for (std::size_t r = 1; r < quads.size(); ++r) {
		CAPTURE(r);
		CHECK(quads[r].size() == k);
	}

	// Row 0 stores the largest absolute length per quad — must be > 0.
	for (std::size_t j = 0; j < k; ++j) {
		CAPTURE(j);
		CHECK(quads[0][j] > 0.0);
	}

	// Rows 1..5 hold ratios relative to row 0 — bounded in (0, 1]. A
	// symmetric input (square grid) can produce equal-length sides, so
	// ratio == 1 is legitimate.
	for (std::size_t r = 1; r <= 5; ++r) {
		for (std::size_t j = 0; j < k; ++j) {
			CAPTURE(r);
			CAPTURE(j);
			CHECK(quads[r][j] > 0.0);
			CHECK(quads[r][j] <= 1.0 + 1e-10);
		}
	}
}

///----------------------------------------
/// MARK: find_quads_xy (for display)
///----------------------------------------

TEST_CASE("find_quads_xy produces a 10xK table") {
	std::vector<std::pair<double, double>> s;
	s.reserve(16);
	for (int iy = 0; iy < 4; ++iy) {
		for (int ix = 0; ix < 4; ++ix) {
			s.emplace_back(200.0 + ix * 80.0, 200.0 + iy * 80.0);
		}
	}
	auto stars = make_starlist(s);

	StarList quads_xy;
	find_quads_xy(stars, quads_xy);

	REQUIRE(quads_xy.size() == 10);
	const auto k = quads_xy[0].size();
	CHECK(k > 0);
	for (std::size_t r = 1; r < 10; ++r) {
		CAPTURE(r);
		CHECK(quads_xy[r].size() == k);
	}
}

///----------------------------------------
/// MARK: find_triples_using_quads
///----------------------------------------

TEST_CASE("find_triples_using_quads emits ~4x as many entries as find_quads") {
	std::vector<std::pair<double, double>> s;
	s.reserve(25);
	for (int iy = 0; iy < 5; ++iy) {
		for (int ix = 0; ix < 5; ++ix) {
			s.emplace_back(100.0 + ix * 50.0, 100.0 + iy * 50.0);
		}
	}
	auto stars_q = make_starlist(s);
	auto stars_t = stars_q;

	StarList quads, triples;
	find_quads(stars_q, quads);
	find_triples_using_quads(stars_t, triples);

	REQUIRE(quads.size()   == 8);
	REQUIRE(triples.size() == 8);

	// Each quad yields up to 4 triples, so expect #triples >= #quads.
	CHECK(triples[0].size() >= quads[0].size());
}
