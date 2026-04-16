///----------------------------------------
///     @file image_sharpness_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for astap::analysis::image_sharpness (src/analysis/image_sharpness.cpp).
///  @details Covers uniform images, high-contrast block patterns, minimum-size
///           images, and the zero-image epsilon guard.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright (C) 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "analysis/image_sharpness.h"

#include <cmath>
#include <limits>

using namespace astap::analysis;
using astap::ImageArray;

/// Helper: create a single-channel image filled with a constant value.
static ImageArray make_uniform(int rows, int cols, float value) {
	ImageArray img(1);
	img[0].assign(rows, std::vector<float>(cols, value));
	return img;
}

///----------------------------------------
/// MARK: Uniform image
///----------------------------------------

TEST_CASE("image_sharpness: uniform image returns large value") {
	auto img = make_uniform(16, 16, 1000.0f);
	double result = image_sharpness(img);
	// No contrast at all -- the inverted metric should be very large.
	CHECK(result > 1e10);
}

///----------------------------------------
/// MARK: High-contrast block pattern
///----------------------------------------

TEST_CASE("image_sharpness: high-contrast blocks return smaller value than uniform") {
	// 16x16 with alternating bright/dark 2x2 blocks so that each 4x4
	// processing region contains both bright and dark sub-quads.
	ImageArray img(1);
	img[0].resize(16, std::vector<float>(16, 0.0f));
	for (int r = 0; r < 16; ++r) {
		for (int c = 0; c < 16; ++c) {
			bool bright = ((r / 2) + (c / 2)) % 2 == 0;
			img[0][r][c] = bright ? 10000.0f : 100.0f;
		}
	}

	double contrast_result = image_sharpness(img);

	auto uniform = make_uniform(16, 16, 1000.0f);
	double uniform_result = image_sharpness(uniform);

	CHECK(contrast_result < uniform_result);
}

///----------------------------------------
/// MARK: Minimum 4x4 image
///----------------------------------------

TEST_CASE("image_sharpness: minimum 4x4 image does not crash") {
	auto img = make_uniform(4, 4, 500.0f);
	double result = image_sharpness(img);
	CHECK(std::isfinite(result));
}

///----------------------------------------
/// MARK: Zero image with epsilon guard
///----------------------------------------

TEST_CASE("image_sharpness: zero image returns finite value") {
	auto img = make_uniform(16, 16, 0.0f);
	double result = image_sharpness(img);
	// The 1e-18 guard in the denominator should prevent NaN/Inf.
	CHECK(std::isfinite(result));
}
