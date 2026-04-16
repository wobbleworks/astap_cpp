///----------------------------------------
///     @file demosaic_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for Bayer demosaic helpers
///           (src/core/demosaic.cpp).
///  @details Focused on the pattern-resolver logic and the simplest
///           variants (demosaic_superpixel: 2x2 down-sample;
///           preserve_colour_saturated_bayer). The more complex
///           bilinear/astro-C/astro-M variants are exercised only for
///           "produces a valid 3-channel output" sanity.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "core/demosaic.h"

#include <vector>

using namespace astap::core;
using astap::ImageArray;
using astap::Header;

///----------------------------------------
/// MARK: Helpers
///----------------------------------------

/// @brief Build a 1-channel mosaic image with the RGGB pattern laid in
///        a repeating 2x2 tile. Each cell type gets a distinct value.
/// @details Pixel (x, y) is:
///            R at even x, even y  →  1000
///            G at odd  x, even y  →  2000
///            G at even x, odd  y  →  3000
///            B at odd  x, odd  y  →  4000
///          (Distinct G1 / G2 values so we can tell them apart.)
[[nodiscard]] static ImageArray make_rggb(int w, int h) {
	ImageArray img(1,
		std::vector<std::vector<float>>(h,
			std::vector<float>(w, 0.0f)));
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			const bool odd_x = (x & 1) != 0;
			const bool odd_y = (y & 1) != 0;
			if (!odd_x && !odd_y)      img[0][y][x] = 1000.0f;
			else if ( odd_x && !odd_y) img[0][y][x] = 2000.0f;
			else if (!odd_x &&  odd_y) img[0][y][x] = 3000.0f;
			else                       img[0][y][x] = 4000.0f;
		}
	}
	return img;
}

///----------------------------------------
/// MARK: get_demosaic_pattern
///----------------------------------------

TEST_CASE("get_demosaic_pattern defaults to RGGB for invalid hints") {
	CHECK(get_demosaic_pattern(-1,  0.0, 0.0, "TOP-DOWN") == 2);   // RGGB
	CHECK(get_demosaic_pattern( 99, 0.0, 0.0, "TOP-DOWN") == 2);
}

TEST_CASE("get_demosaic_pattern identity for clean inputs") {
	for (int p = 0; p <= 4; ++p) {
		CAPTURE(p);
		CHECK(get_demosaic_pattern(p, 0.0, 0.0, "TOP-DOWN") == p);
	}
}

TEST_CASE("X-Trans is never parity-flipped") {
	// Pattern 4 stays as 4 even with odd xbayroff / ybayroff.
	CHECK(get_demosaic_pattern(4, 1.0, 1.0, "BOT-UP") == 4);
	CHECK(get_demosaic_pattern(4, 0.0, 3.0, "TOP-DOWN") == 4);
}

TEST_CASE("get_demosaic_pattern applies XBAYROFF parity swap") {
	// XBAYROFF=1 swaps RGGB<->GRBG, BGGR<->GBRG.
	CHECK(get_demosaic_pattern(2, 1.0, 0.0, "TOP-DOWN") == 0);  // RGGB->GRBG
	CHECK(get_demosaic_pattern(0, 1.0, 0.0, "TOP-DOWN") == 2);  // GRBG->RGGB
	CHECK(get_demosaic_pattern(1, 1.0, 0.0, "TOP-DOWN") == 3);  // BGGR->GBRG
	CHECK(get_demosaic_pattern(3, 1.0, 0.0, "TOP-DOWN") == 1);  // GBRG->BGGR
}

TEST_CASE("get_demosaic_pattern applies YBAYROFF parity swap") {
	// YBAYROFF=1 swaps RGGB<->GBRG, BGGR<->GRBG.
	CHECK(get_demosaic_pattern(2, 0.0, 1.0, "TOP-DOWN") == 3);  // RGGB->GBRG
	CHECK(get_demosaic_pattern(3, 0.0, 1.0, "TOP-DOWN") == 2);  // GBRG->RGGB
	CHECK(get_demosaic_pattern(1, 0.0, 1.0, "TOP-DOWN") == 0);  // BGGR->GRBG
	CHECK(get_demosaic_pattern(0, 0.0, 1.0, "TOP-DOWN") == 1);  // GRBG->BGGR
}

TEST_CASE("BOT-UP row order adds +1 to ybayroff parity") {
	// Even ybayroff + BOT-UP → effective odd ybayroff → parity swap.
	CHECK(get_demosaic_pattern(2, 0.0, 0.0, "BOT-UP") == 3);     // RGGB -> GBRG
	// Odd ybayroff + BOT-UP → effective even → no swap.
	CHECK(get_demosaic_pattern(2, 0.0, 1.0, "BOT-UP") == 2);
}

///----------------------------------------
/// MARK: demosaic_superpixel (2x2 down-sample)
///----------------------------------------

TEST_CASE("demosaic_superpixel halves dimensions and separates channels (RGGB)") {
	constexpr int w = 8;
	constexpr int h = 8;
	auto img = make_rggb(w, h);

	Header head{};
	head.width   = w;
	head.height  = h;
	head.naxis   = 2;
	head.naxis3  = 1;

	demosaic_superpixel(img, head, /*pattern=*/2);   // RGGB

	CHECK(head.naxis  == 3);
	CHECK(head.naxis3 == 3);
	CHECK(head.width  == w / 2);
	CHECK(head.height == h / 2);

	REQUIRE(img.size() == 3);
	REQUIRE(img[0].size() == static_cast<std::size_t>(h / 2));
	REQUIRE(img[0][0].size() == static_cast<std::size_t>(w / 2));

	// Every super-pixel in the R channel should be 1000 (only R cells fed it).
	// Channel 0 = red, 1 = green, 2 = blue.
	for (int y = 0; y < h / 2; ++y) {
		for (int x = 0; x < w / 2; ++x) {
			CAPTURE(x);
			CAPTURE(y);
			CHECK(img[0][y][x] == doctest::Approx(1000.0f));   // red
			// Green is the average of the two green cells (2000 + 3000) / 2 = 2500.
			CHECK(img[1][y][x] == doctest::Approx(2500.0f));
			CHECK(img[2][y][x] == doctest::Approx(4000.0f));   // blue
		}
	}
}

///----------------------------------------
/// MARK: demosaic_bilinear_interpolation (spot check)
///----------------------------------------

TEST_CASE("demosaic_bilinear_interpolation keeps same size, produces 3 channels") {
	constexpr int w = 12;
	constexpr int h = 12;
	auto img = make_rggb(w, h);

	Header head{};
	head.width   = w;
	head.height  = h;
	head.naxis   = 2;
	head.naxis3  = 1;

	demosaic_bilinear_interpolation(img, head, /*pattern=*/2);

	CHECK(head.naxis  == 3);
	CHECK(head.naxis3 == 3);
	// Same dimensions.
	CHECK(head.width  == w);
	CHECK(head.height == h);
	REQUIRE(img.size() == 3);
	REQUIRE(img[0].size() == static_cast<std::size_t>(h));
	REQUIRE(img[0][0].size() == static_cast<std::size_t>(w));
}

///----------------------------------------
/// MARK: demosaic_astrosimple (interior pixel values preserved)
///----------------------------------------

TEST_CASE("demosaic_astrosimple keeps the cell's own colour at its pixel") {
	constexpr int w = 8;
	constexpr int h = 8;
	auto img = make_rggb(w, h);

	Header head{};
	head.width   = w;
	head.height  = h;
	head.naxis   = 2;
	head.naxis3  = 1;

	demosaic_astrosimple(img, head, /*pattern=*/2);

	CHECK(head.naxis3 == 3);

	// In RGGB, the red pixel at (0,0) should have kept its 1000 in channel 0.
	// (Other channels get approximations from neighbours; we don't constrain
	// them here.)
	CHECK(img[0][0][0] == doctest::Approx(1000.0f));
}
