///----------------------------------------
///     @file contour_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for astap::analysis contour helpers (src/analysis/contour.cpp).
///  @details Covers line_distance, trendline, and trendline_without_outliers.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright (C) 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "analysis/contour.h"

#include <cmath>

using namespace astap::analysis;
using astap::StarList;

///----------------------------------------
/// MARK: line_distance
///----------------------------------------

TEST_CASE("line_distance: horizontal line") {
	// y = 0*x + 5 (horizontal at y=5).
	CHECK(line_distance(3.0, 5.0, 0.0, 5.0) == doctest::Approx(0.0));
	CHECK(line_distance(3.0, 10.0, 0.0, 5.0) == doctest::Approx(5.0));
}

TEST_CASE("line_distance: y = x") {
	// slope=1, intercept=0.  Point (1,0): distance = |1*1 - 0 + 0| / sqrt(2)
	CHECK(line_distance(1.0, 0.0, 1.0, 0.0)
	      == doctest::Approx(1.0 / std::sqrt(2.0)).epsilon(1e-10));
}

///----------------------------------------
/// MARK: trendline — perfect fit
///----------------------------------------

TEST_CASE("trendline: perfect fit on y = 2x + 1") {
	StarList xy(2);
	for (int i = 0; i < 5; ++i) {
		double x = static_cast<double>(i);
		xy[0].push_back(x);
		xy[1].push_back(2.0 * x + 1.0);
	}

	double slope = 0.0, intercept = 0.0;
	trendline(xy, 5, slope, intercept);

	CHECK(slope == doctest::Approx(2.0).epsilon(1e-10));
	CHECK(intercept == doctest::Approx(1.0).epsilon(1e-10));
}

///----------------------------------------
/// MARK: trendline — noisy data
///----------------------------------------

TEST_CASE("trendline: approximate fit on y = 0.5x + 3 with noise") {
	StarList xy(2);
	// Small perturbations that average out.
	const double noise[] = {0.05, -0.03, 0.02, -0.04, 0.01, 0.03, -0.02, 0.01};
	for (int i = 0; i < 8; ++i) {
		double x = static_cast<double>(i);
		xy[0].push_back(x);
		xy[1].push_back(0.5 * x + 3.0 + noise[i]);
	}

	double slope = 0.0, intercept = 0.0;
	trendline(xy, 8, slope, intercept);

	CHECK(slope == doctest::Approx(0.5).epsilon(0.05));
	CHECK(intercept == doctest::Approx(3.0).epsilon(0.1));
}

///----------------------------------------
/// MARK: trendline_without_outliers
///----------------------------------------

TEST_CASE("trendline_without_outliers: rejects wild outlier on y = x") {
	StarList xy(2);
	// 10 points on y = x.
	for (int i = 0; i < 10; ++i) {
		double x = static_cast<double>(i);
		xy[0].push_back(x);
		xy[1].push_back(x);
	}
	// One wild outlier.
	xy[0].push_back(5.0);
	xy[1].push_back(100.0);

	double slope = 0.0, intercept = 0.0, sd = 0.0;
	trendline_without_outliers(xy, 11, slope, intercept, sd);

	CHECK(slope == doctest::Approx(1.0).epsilon(0.1));
	CHECK(intercept == doctest::Approx(0.0).epsilon(0.5));
	CHECK(sd > 0.0);
}
