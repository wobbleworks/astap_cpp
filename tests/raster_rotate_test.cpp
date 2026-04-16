///----------------------------------------
///     @file raster_rotate_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for accurate raster rotation
///           (src/stacking/raster_rotate.cpp).
///  @details 0° is identity, 90°/180°/270° are lossless pixel shuffles,
///           forward + reverse rotation round-trips to near-identity.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "stacking/raster_rotate.h"

#include <cmath>
#include <vector>

using namespace astap::stacking;

///----------------------------------------
/// MARK: Helpers
///----------------------------------------

[[nodiscard]] static ImageArray make_mono(int width, int height, float fill = 0.0f) {
	return ImageArray(1,
		std::vector<std::vector<float>>(height,
			std::vector<float>(width, fill)));
}

/// @brief Sum of every pixel. Used to confirm area-weighted rotation
///        preserves total flux to within a small tolerance.
[[nodiscard]] static double total_flux(const ImageArray& img) {
	double sum = 0.0;
	for (const auto& channel : img) {
		for (const auto& row : channel) {
			for (float v : row) {
				sum += v;
			}
		}
	}
	return sum;
}

///----------------------------------------
/// MARK: Identity rotation
///----------------------------------------

TEST_CASE("0 degrees is the identity") {
	constexpr int w = 16;
	constexpr int h = 16;
	auto img = make_mono(w, h);
	// Asymmetric fill so we can see any displacement.
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			img[0][y][x] = static_cast<float>(x + y * w);
		}
	}
	const auto before = img;

	raster_rotate(0.0, w / 2.0, h / 2.0, img);

	// Dimensions may expand, but the original-sized interior should match.
	REQUIRE(!img.empty());
	REQUIRE(img[0].size() >= static_cast<std::size_t>(h));
	REQUIRE(img[0][0].size() >= static_cast<std::size_t>(w));

	const int new_h = static_cast<int>(img[0].size());
	const int new_w = static_cast<int>(img[0][0].size());
	const int dy = (new_h - h) / 2;
	const int dx = (new_w - w) / 2;

	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			CAPTURE(x);
			CAPTURE(y);
			CHECK(img[0][y + dy][x + dx] == doctest::Approx(before[0][y][x]).epsilon(1e-4));
		}
	}
}

///----------------------------------------
/// MARK: Quarter-turn round-trip
///----------------------------------------

TEST_CASE("rotate(90) then rotate(-90) returns close to original") {
	constexpr int w = 32;
	constexpr int h = 32;
	auto img = make_mono(w, h);
	// Mark a single bright pixel off-centre.
	img[0][10][12] = 1000.0f;
	const auto before = img;

	raster_rotate( 90.0, w / 2.0, h / 2.0, img);
	raster_rotate(-90.0, img[0][0].size() / 2.0,
	                     img[0].size() / 2.0, img);

	// The pipeline expands the canvas each call, so land the original
	// interior inside the expanded buffer and compare interior windows.
	const int new_h = static_cast<int>(img[0].size());
	const int new_w = static_cast<int>(img[0][0].size());
	REQUIRE(new_h >= h);
	REQUIRE(new_w >= w);

	const int dy = (new_h - h) / 2;
	const int dx = (new_w - w) / 2;
	CHECK(img[0][10 + dy][12 + dx] == doctest::Approx(1000.0f).epsilon(1e-3));
}

///----------------------------------------
/// MARK: 180° two-step round-trip
///----------------------------------------

TEST_CASE("rotate(180) twice returns close to original") {
	constexpr int w = 20;
	constexpr int h = 20;
	auto img = make_mono(w, h);
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			img[0][y][x] = static_cast<float>(x * 3 + y);
		}
	}
	const auto before = img;

	raster_rotate(180.0, w / 2.0, h / 2.0, img);
	const int sz1_h = static_cast<int>(img[0].size());
	const int sz1_w = static_cast<int>(img[0][0].size());
	raster_rotate(180.0, sz1_w / 2.0, sz1_h / 2.0, img);

	const int new_h = static_cast<int>(img[0].size());
	const int new_w = static_cast<int>(img[0][0].size());
	const int dy = (new_h - h) / 2;
	const int dx = (new_w - w) / 2;
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			CAPTURE(x);
			CAPTURE(y);
			CHECK(img[0][y + dy][x + dx] == doctest::Approx(before[0][y][x]).epsilon(1e-3));
		}
	}
}

///----------------------------------------
/// MARK: Flux preservation
///----------------------------------------

TEST_CASE("arbitrary rotation approximately preserves total flux") {
	constexpr int w = 40;
	constexpr int h = 40;
	auto img = make_mono(w, h);
	// A compact central blob so rotation never pushes flux off the canvas.
	for (int y = 18; y < 22; ++y) {
		for (int x = 18; x < 22; ++x) {
			img[0][y][x] = 100.0f;
		}
	}
	const double original = total_flux(img);

	raster_rotate(37.5, w / 2.0, h / 2.0, img);

	// Area-weighted interpolation is designed to conserve flux; allow
	// 1% slack for edge-pixel partial-coverage rounding.
	CHECK(total_flux(img) == doctest::Approx(original).epsilon(1e-2));
}

///----------------------------------------
/// MARK: Small rotations are smooth
///----------------------------------------

TEST_CASE("tiny angle does not introduce large artefacts in interior") {
	constexpr int w = 40;
	constexpr int h = 40;
	auto img = make_mono(w, h, 50.0f);   // constant image

	raster_rotate(1.0, w / 2.0, h / 2.0, img);

	// Expect the constant value to survive in the expanded interior.
	const int new_h = static_cast<int>(img[0].size());
	const int new_w = static_cast<int>(img[0][0].size());

	// Sample interior (far from edges) — should remain ≈ 50.
	for (int y = new_h / 3; y < 2 * new_h / 3; ++y) {
		for (int x = new_w / 3; x < 2 * new_w / 3; ++x) {
			CAPTURE(x);
			CAPTURE(y);
			CHECK(img[0][y][x] == doctest::Approx(50.0f).epsilon(2e-2));
		}
	}
}
