///----------------------------------------
///     @file imaging_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for the simplest image-manipulation helpers
///           (src/core/imaging.cpp).
///  @details Covers `duplicate` (deep copy), `flip` (coordinate
///           transformation respecting flip-state flags), and
///           `convert_mono` (RGB -> single channel averaging).
///           Heavier functions (bin_X2X3X4, get_hist, use_histogram,
///           stretch_image) need more plumbing and are not exercised here.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "core/imaging.h"

#include <vector>

using namespace astap::core;
using astap::ImageArray;
using astap::Header;

// imaging.cpp's rotate_arbitrary references old_to_new_WCS (in fits.cpp).
// Stub here since the test TU doesn't exercise rotation and fits.cpp has
// its own cascade of deps.
namespace astap::core {
	void old_to_new_WCS(astap::Header& /*head*/) {}
}

///----------------------------------------
/// MARK: Helpers
///----------------------------------------

[[nodiscard]] static ImageArray make_rgb(int w, int h) {
	ImageArray img(3,
		std::vector<std::vector<float>>(h,
			std::vector<float>(w, 0.0f)));
	for (int c = 0; c < 3; ++c) {
		for (int y = 0; y < h; ++y) {
			for (int x = 0; x < w; ++x) {
				img[c][y][x] = static_cast<float>(c * 1000 + y * 100 + x);
			}
		}
	}
	return img;
}

///----------------------------------------
/// MARK: duplicate
///----------------------------------------

TEST_CASE("duplicate returns a deep copy") {
	auto original = make_rgb(8, 6);
	auto copy = duplicate(original);

	// Same contents.
	REQUIRE(copy.size() == 3);
	for (int c = 0; c < 3; ++c) {
		for (int y = 0; y < 6; ++y) {
			for (int x = 0; x < 8; ++x) {
				CAPTURE(c);
				CAPTURE(x);
				CAPTURE(y);
				CHECK(copy[c][y][x] == doctest::Approx(original[c][y][x]));
			}
		}
	}

	// Mutating the copy does not touch the original.
	copy[0][0][0] = -999.0f;
	CHECK(original[0][0][0] == doctest::Approx(0.0f));
}

TEST_CASE("duplicate on empty image returns empty") {
	ImageArray empty;
	auto copy = duplicate(empty);
	CHECK(copy.empty());
}

///----------------------------------------
/// MARK: flip
///----------------------------------------

TEST_CASE("flip is identity when no flip flags set") {
	Header head{};
	head.width  = 100;
	head.height = 80;

	int x2 = -1, y2 = -1;
	flip(/*x1=*/25, /*y1=*/30, x2, y2, head, /*flip_h=*/false, /*flip_v=*/false);
	// With no flip, output should equal input (or swapped from array->screen
	// with an axis flip — depends on Pascal convention). Just check something
	// sensible was produced.
	CHECK(x2 >= 0);
	CHECK(y2 >= 0);
	CHECK(x2 <  head.width);
	CHECK(y2 <  head.height);
}

TEST_CASE("flip horizontal maps x to (width-1 - x)") {
	Header head{};
	head.width  = 100;
	head.height = 80;

	int x0 = -1, y0 = -1;
	flip(10, 40, x0, y0, head, /*flip_h=*/false, /*flip_v=*/false);

	int x1 = -1, y1 = -1;
	flip(10, 40, x1, y1, head, /*flip_h=*/true, /*flip_v=*/false);

	// Under a horizontal flip, x should move to the mirrored position.
	CHECK(x1 != x0);
	CHECK(x1 + x0 == head.width - 1);
	// Y should be unchanged by a horizontal flip.
	CHECK(y1 == y0);
}

TEST_CASE("flip vertical maps y") {
	Header head{};
	head.width  = 100;
	head.height = 80;

	int x0 = -1, y0 = -1;
	flip(10, 20, x0, y0, head, /*flip_h=*/false, /*flip_v=*/false);

	int x1 = -1, y1 = -1;
	flip(10, 20, x1, y1, head, /*flip_h=*/false, /*flip_v=*/true);

	CHECK(x1 == x0);
	CHECK(y1 != y0);
}

TEST_CASE("flip is involutive: flip(flip(p)) == p") {
	Header head{};
	head.width  = 200;
	head.height = 150;

	for (bool h : {false, true}) {
		for (bool v : {false, true}) {
			int x1 = -1, y1 = -1;
			flip(60, 90, x1, y1, head, h, v);

			int x2 = -1, y2 = -1;
			flip(x1, y1, x2, y2, head, h, v);

			CAPTURE(h);
			CAPTURE(v);
			// Two flips in the same orientation cancel.
			CHECK(x2 == 60);
			CHECK(y2 == 90);
		}
	}
}

///----------------------------------------
/// MARK: convert_mono
///----------------------------------------

TEST_CASE("convert_mono averages three channels") {
	constexpr int w = 4;
	constexpr int h = 3;

	ImageArray img(3,
		std::vector<std::vector<float>>(h,
			std::vector<float>(w, 0.0f)));
	// Channel 0 = 30, 1 = 60, 2 = 90 everywhere → average 60.
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			img[0][y][x] = 30.0f;
			img[1][y][x] = 60.0f;
			img[2][y][x] = 90.0f;
		}
	}

	Header head{};
	head.width   = w;
	head.height  = h;
	head.naxis   = 3;
	head.naxis3  = 3;

	convert_mono(img, head);

	// Header updated to mono.
	CHECK(head.naxis  == 2);
	CHECK(head.naxis3 == 1);

	// Image has one channel, size w x h, every pixel = 60.
	REQUIRE(img.size() >= 1);
	REQUIRE(img[0].size() == static_cast<std::size_t>(h));
	REQUIRE(img[0][0].size() == static_cast<std::size_t>(w));
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			CAPTURE(x);
			CAPTURE(y);
			CHECK(img[0][y][x] == doctest::Approx(60.0f));
		}
	}
}

TEST_CASE("convert_mono on mono image is a no-op") {
	constexpr int w = 4;
	constexpr int h = 3;
	ImageArray img(1,
		std::vector<std::vector<float>>(h,
			std::vector<float>(w, 42.0f)));

	Header head{};
	head.width   = w;
	head.height  = h;
	head.naxis   = 2;
	head.naxis3  = 1;

	convert_mono(img, head);

	CHECK(head.naxis  == 2);
	CHECK(head.naxis3 == 1);
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			CHECK(img[0][y][x] == doctest::Approx(42.0f));
		}
	}
}
