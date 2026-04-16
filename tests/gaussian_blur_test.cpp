///----------------------------------------
///     @file gaussian_blur_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for the separable Gaussian blur
///           (src/stacking/gaussian_blur.cpp).
///  @details Zero-radius passthrough, constant-image invariance, total-flux
///           conservation, symmetry preservation, and monotonic smoothing
///           of an impulse response.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "stacking/gaussian_blur.h"

#include <cmath>
#include <vector>

using namespace astap::stacking;

///----------------------------------------
/// MARK: Helpers
///----------------------------------------

/// @brief Make a single-channel (mono) image buffer of @p width x @p height.
[[nodiscard]] static ImageArray make_mono(int width, int height, float fill = 0.0f) {
	ImageArray img(1,
		std::vector<std::vector<float>>(height,
			std::vector<float>(width, fill)));
	return img;
}

/// @brief Sum every pixel across every channel.
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
/// MARK: Tests
///----------------------------------------

TEST_CASE("gaussian_blur2 with radius < 0.001 is a no-op") {
	constexpr int w = 8;
	constexpr int h = 8;
	auto img = make_mono(w, h);
	img[0][4][4] = 1000.0f;   // single bright pixel

	gaussian_blur2(img, 0.0);
	CHECK(img[0][4][4] == doctest::Approx(1000.0f));

	// Other pixels still zero.
	CHECK(img[0][0][0] == doctest::Approx(0.0f));
	CHECK(img[0][7][7] == doctest::Approx(0.0f));
}

TEST_CASE("gaussian_blur2 on a constant image leaves interior pixels unchanged") {
	// Edge pixels can shift slightly due to boundary handling; check only
	// the interior where the kernel fits entirely inside the image.
	constexpr int w = 20;
	constexpr int h = 20;
	auto img = make_mono(w, h, 42.0f);

	gaussian_blur2(img, 2.0);

	for (int y = 5; y < h - 5; ++y) {
		for (int x = 5; x < w - 5; ++x) {
			CAPTURE(x);
			CAPTURE(y);
			CHECK(img[0][y][x] == doctest::Approx(42.0f).epsilon(1e-4));
		}
	}
}

TEST_CASE("gaussian_blur2 approximately preserves total flux") {
	// A separable Gaussian with unit-sum weights conserves flux modulo
	// edge effects. Use a compact central blob so the kernel stays away
	// from the borders.
	constexpr int w = 64;
	constexpr int h = 64;
	auto img = make_mono(w, h);
	img[0][32][32] = 1000.0f;
	const double original = total_flux(img);

	gaussian_blur2(img, 3.0);
	CHECK(total_flux(img) == doctest::Approx(original).epsilon(1e-3));
}

TEST_CASE("gaussian_blur2 preserves 4-fold symmetry of a centred impulse") {
	constexpr int w = 33;  // odd → exact centre pixel
	constexpr int h = 33;
	auto img = make_mono(w, h);
	img[0][16][16] = 100.0f;

	gaussian_blur2(img, 2.5);

	// Pick a few offsets and confirm pixels equidistant from centre match.
	for (int d : {1, 2, 3, 5}) {
		CAPTURE(d);
		const float left  = img[0][16][16 - d];
		const float right = img[0][16][16 + d];
		const float up    = img[0][16 - d][16];
		const float down  = img[0][16 + d][16];
		CHECK(left  == doctest::Approx(right).epsilon(1e-4));
		CHECK(up    == doctest::Approx(down ).epsilon(1e-4));
		CHECK(left  == doctest::Approx(up   ).epsilon(1e-4));
	}
}

TEST_CASE("gaussian_blur2 produces monotonically decreasing values from a centre impulse") {
	constexpr int w = 33;
	constexpr int h = 33;
	auto img = make_mono(w, h);
	img[0][16][16] = 100.0f;

	gaussian_blur2(img, 3.0);

	// Moving outward along one axis, values should decrease monotonically.
	float prev = img[0][16][16];
	for (int d = 1; d <= 10; ++d) {
		const float v = img[0][16][16 + d];
		CAPTURE(d);
		CHECK(v <= prev);
		prev = v;
	}
}

TEST_CASE("gaussian_blur2 spreads a point source to non-adjacent pixels") {
	// A radius-3 Gaussian should meaningfully spread the impulse at least
	// a few pixels from the source.
	constexpr int w = 33;
	constexpr int h = 33;
	auto img = make_mono(w, h);
	img[0][16][16] = 1000.0f;

	gaussian_blur2(img, 3.0);

	// Pixel at distance 3 should be non-trivially non-zero.
	CHECK(img[0][16][16 + 3] > 0.1f);
	// Peak at centre should have dropped considerably (the 1000 has spread
	// over many neighbours).
	CHECK(img[0][16][16] < 1000.0f);
}
