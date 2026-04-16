///----------------------------------------
///     @file inspector_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for CCD inspector helpers (src/analysis/inspector.cpp).
///  @details Covers fnmodulo2, filter_hfd, and measure_star_aspect.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright (C) 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "analysis/inspector.h"

#include <cmath>
#include <vector>

using namespace astap::analysis;
using astap::ImageArray;
using astap::StarList;

///----------------------------------------
/// MARK: fnmodulo2 — basic wrapping
///----------------------------------------

TEST_CASE("fnmodulo2: basic wrapping") {
	CHECK(fnmodulo2(190.0, 360.0) == doctest::Approx(-170.0));
	CHECK(fnmodulo2(-200.0, 360.0) == doctest::Approx(160.0));
	CHECK(fnmodulo2(45.0, 360.0) == doctest::Approx(45.0));
}

///----------------------------------------
/// MARK: fnmodulo2 — exact boundary
///----------------------------------------

TEST_CASE("fnmodulo2: exact boundary") {
	// At exactly +/- range/2, the while-loop condition is > or <,
	// so 180 stays at 180 and -180 stays at -180.
	CHECK(fnmodulo2(180.0, 360.0) == doctest::Approx(180.0));
	CHECK(fnmodulo2(-180.0, 360.0) == doctest::Approx(-180.0));
}

///----------------------------------------
/// MARK: filter_hfd — 3 collinear stars
///----------------------------------------

TEST_CASE("filter_hfd: 3 collinear stars") {
	// Stars at (0,0), (10,0), (20,0) with HFD*100 values 100, 200, 300.
	StarList hfd(3);
	hfd[0] = {0.0, 10.0, 20.0};  // X
	hfd[1] = {0.0, 0.0, 0.0};    // Y
	hfd[2] = {100.0, 200.0, 300.0};  // HFD * 100

	float mean = 0, min_val = 0, max_val = 0;
	filter_hfd(hfd, 3, mean, min_val, max_val);

	// Each star's two nearest neighbours are the other two stars.
	// Star 0: median(100, 200, 300) = 200
	// Star 1: median(200, 100, 300) = 200
	// Star 2: median(300, 200, 100) = 200
	CHECK(hfd[2][0] == doctest::Approx(200.0));
	CHECK(hfd[2][1] == doctest::Approx(200.0));
	CHECK(hfd[2][2] == doctest::Approx(200.0));
}

///----------------------------------------
/// MARK: filter_hfd — output stats
///----------------------------------------

TEST_CASE("filter_hfd: output stats") {
	StarList hfd(3);
	hfd[0] = {0.0, 10.0, 20.0};
	hfd[1] = {0.0, 0.0, 0.0};
	hfd[2] = {100.0, 200.0, 300.0};

	float mean = 0, min_val = 0, max_val = 0;
	filter_hfd(hfd, 3, mean, min_val, max_val);

	CHECK(mean >= min_val);
	CHECK(mean <= max_val);
	CHECK(min_val <= max_val);
}

///----------------------------------------
/// MARK: measure_star_aspect — uniform image
///----------------------------------------

TEST_CASE("measure_star_aspect: uniform image returns failure") {
	ImageArray img(1);
	img[0].assign(64, std::vector<float>(64, 100.0f));

	// Star at centre with generous radius; star_bg=100, sd_bg=1.
	// No pixels exceed 7*sd_bg above background, so expect failure.
	auto result = measure_star_aspect(img, 32.0, 32.0, 10, 100.0, 1.0);
	CHECK(result.aspect == doctest::Approx(999.0));
}

///----------------------------------------
/// MARK: measure_star_aspect — out of bounds
///----------------------------------------

TEST_CASE("measure_star_aspect: star at edge returns failure") {
	ImageArray img(1);
	img[0].assign(64, std::vector<float>(64, 100.0f));

	// Star at edge with large radius — bounds check should fail.
	auto result = measure_star_aspect(img, 2.0, 2.0, 20, 100.0, 1.0);
	CHECK(result.aspect == doctest::Approx(999.0));
}
