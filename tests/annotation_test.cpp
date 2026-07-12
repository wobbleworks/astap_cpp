///----------------------------------------
///     @file annotation_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for annotation utilities (src/analysis/annotation.cpp).
///  @details Covers the 5x9 bitmap font tables, text-to-image rendering,
///           the rotate() helper, equatorial_standard, and find_object.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright (C) 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "analysis/annotation.h"

#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace astap::analysis;
using astap::ImageArray;

///----------------------------------------
/// MARK: Font table dimensions
///----------------------------------------

TEST_CASE("kFont5x9: dimensions and value range") {
	// 94 characters (ASCII 33..126), each 9 rows of 5 columns.
	static_assert(sizeof(kFont5x9) == 94 * 9 * 5);

	for (int ch = 0; ch < 94; ++ch) {
		for (int row = 0; row < 9; ++row) {
			for (int col = 0; col < 5; ++col) {
				CHECK_MESSAGE((kFont5x9[ch][row][col] == 0 ||
				               kFont5x9[ch][row][col] == 1),
				              "invalid pixel at [" << ch << "][" << row << "][" << col << "]");
			}
		}
	}
}

///----------------------------------------
/// MARK: Font '!' character
///----------------------------------------

TEST_CASE("kFont5x9: '!' character shape") {
	// kFont5x9[0] = '!' (ASCII 33). Column 2 is set in rows 0-6, blank row 7,
	// dot at row 8 col 2.
	for (int row = 0; row <= 6; ++row) {
		CHECK(kFont5x9[0][row][2] == 1);
	}
	CHECK(kFont5x9[0][7][2] == 0);
	CHECK(kFont5x9[0][8][2] == 1);
}

///----------------------------------------
/// MARK: Font 'A' character
///----------------------------------------

TEST_CASE("kFont5x9: 'A' character row 0") {
	// kFont5x9[32] = 'A' (ASCII 65). Row 0 = {0,0,1,0,0}.
	CHECK(kFont5x9[32][0][0] == 0);
	CHECK(kFont5x9[32][0][1] == 0);
	CHECK(kFont5x9[32][0][2] == 1);
	CHECK(kFont5x9[32][0][3] == 0);
	CHECK(kFont5x9[32][0][4] == 0);
}

///----------------------------------------
/// MARK: annotation_to_array writes pixels
///----------------------------------------

TEST_CASE("annotation_to_array: writes pixels for 'A'") {
	ImageArray img(1);
	img[0].assign(50, std::vector<float>(50, 0.0f));

	annotation_to_array("A", false, 1000, 1, 5, 20, img);

	// Some pixels at the expected location should be non-zero.
	bool found_nonzero = false;
	for (int r = 10; r <= 25; ++r) {
		for (int c = 5; c < 12; ++c) {
			if (img[0][r][c] != 0.0f) {
				found_nonzero = true;
				break;
			}
		}
		if (found_nonzero) break;
	}
	CHECK(found_nonzero);
}

///----------------------------------------
/// MARK: annotation_to_array transparent mode
///----------------------------------------

TEST_CASE("annotation_to_array: transparent mode preserves background") {
	ImageArray img(1);
	img[0].assign(50, std::vector<float>(50, 42.0f));

	annotation_to_array("A", true, 1000, 1, 5, 20, img);

	// In transparent mode, pixels where font=0 should retain 42.0.
	// Find at least one background pixel in the glyph bounding box.
	bool found_preserved = false;
	for (int r = 12; r <= 20; ++r) {
		for (int c = 5; c < 10; ++c) {
			if (img[0][r][c] == 42.0f) {
				found_preserved = true;
				break;
			}
		}
		if (found_preserved) break;
	}
	CHECK(found_preserved);
}

///----------------------------------------
/// MARK: rotate identity
///----------------------------------------

TEST_CASE("rotate: identity at rot=0") {
	auto [x2, y2] = rotate(0.0, 1.0, 0.0);
	// x2 = x*sin(0) + y*cos(0) = 0, y2 = -x*cos(0) + y*sin(0) = -1
	CHECK(x2 == doctest::Approx(0.0).epsilon(1e-12));
	CHECK(y2 == doctest::Approx(-1.0).epsilon(1e-12));
}

TEST_CASE("rotate: rot=0 with (0, 1)") {
	auto [x2, y2] = rotate(0.0, 0.0, 1.0);
	// x2 = 0*sin(0) + 1*cos(0) = 1, y2 = -0*cos(0) + 1*sin(0) = 0
	CHECK(x2 == doctest::Approx(1.0).epsilon(1e-12));
	CHECK(y2 == doctest::Approx(0.0).epsilon(1e-12));
}

///----------------------------------------
/// MARK: equatorial_standard at reference point
///----------------------------------------

TEST_CASE("equatorial_standard: at reference point returns (0, 0)") {
	double ra0 = 1.5;
	double dec0 = 0.8;
	auto [xx, yy] = equatorial_standard(ra0, dec0, ra0, dec0, 1.0);
	CHECK(xx == doctest::Approx(0.0).epsilon(1e-8));
	CHECK(yy == doctest::Approx(0.0).epsilon(1e-8));
}

///----------------------------------------
/// MARK: find_object with mock database
///----------------------------------------

TEST_CASE("find_object: finds object in mock database") {
	// Minimal CSV database matching the expected format.
	// Lines 0 and 1 are headers (skipped). Line 2+ are data.
	// Fields: RA_encoded, Dec_encoded, names, length, width, PA
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"432000,324000, M31/NGC224,178.0,63.0,35"
	};

	auto result = find_object("M31", db);
	REQUIRE(result.has_value());
	CHECK(result->name == "M31_NGC224");
	CHECK(result->length == doctest::Approx(178.0));
	CHECK(result->width == doctest::Approx(63.0));
	CHECK(result->pa == doctest::Approx(35.0));
}

TEST_CASE("find_object: returns nullopt for missing object") {
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"432000,324000, M31/NGC224,178.0,63.0,35"
	};

	auto result = find_object("M42", db);
	CHECK_FALSE(result.has_value());
}

///----------------------------------------
/// MARK: read_deepsky field parsing
///----------------------------------------

TEST_CASE("read_deepsky: parses a record inside the field") {
	// RA 432000 -> pi rad; Dec 0 -> 0 rad.
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"432000,0, M31/NGC224,178.0,63.0,35"
	};

	auto objs = read_deepsky(db, DeepSkySearch::WithinField,
	                         std::numbers::pi, 0.0, 0.1);
	REQUIRE(objs.size() == 1);
	CHECK(objs[0].name == "M31");
	CHECK(objs[0].name2 == "NGC224");
	CHECK(objs[0].name3.empty());
	CHECK(objs[0].ra == doctest::Approx(std::numbers::pi).epsilon(1e-9));
	CHECK(objs[0].dec == doctest::Approx(0.0).epsilon(1e-12));
	CHECK(objs[0].length == doctest::Approx(178.0));
	CHECK(objs[0].width == doctest::Approx(63.0));
	CHECK(objs[0].pa == doctest::Approx(35.0));
}

///----------------------------------------
/// MARK: read_deepsky FOV gate
///----------------------------------------

TEST_CASE("read_deepsky: FOV gate excludes distant objects") {
	// First object at the field centre (pi, 0); second far away (0, pi/2).
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"432000,0, M31,178.0,63.0,35",
		"0,324000, FAR,10.0,10.0,0"
	};

	auto objs = read_deepsky(db, DeepSkySearch::WithinField,
	                         std::numbers::pi, 0.0, 0.1);
	REQUIRE(objs.size() == 1);
	CHECK(objs[0].name == "M31");
}

TEST_CASE("read_deepsky: FullDatabase ignores the FOV gate") {
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"432000,0, M31,178.0,63.0,35",
		"0,324000, FAR,10.0,10.0,0"
	};

	auto objs = read_deepsky(db, DeepSkySearch::FullDatabase,
	                         std::numbers::pi, 0.0, 0.1);
	CHECK(objs.size() == 2);
}

///----------------------------------------
/// MARK: read_deepsky RA wraparound
///----------------------------------------

TEST_CASE("read_deepsky: RA difference wraps across 0/2pi") {
	// RA 863999 -> ~2pi, i.e. ~7.3e-6 rad short of a full turn.  A field
	// centred at RA 0 must include it via the wrapped delta_ra fold.
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"863999,0, WRAP,10.0,5.0,0"
	};

	auto objs = read_deepsky(db, DeepSkySearch::WithinField,
	                         0.0, 0.0, 0.001);
	REQUIRE(objs.size() == 1);
	CHECK(objs[0].name == "WRAP");
}

///----------------------------------------
/// MARK: read_deepsky unknown position angle
///----------------------------------------

TEST_CASE("read_deepsky: unparseable PA maps to 999") {
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"432000,0, NOPA,10.0,5.0,"
	};

	auto objs = read_deepsky(db, DeepSkySearch::WithinField,
	                         std::numbers::pi, 0.0, 0.1);
	REQUIRE(objs.size() == 1);
	CHECK(objs[0].pa == doctest::Approx(999.0));
}

///----------------------------------------
/// MARK: read_deepsky skips header lines
///----------------------------------------

TEST_CASE("read_deepsky: skips the two header lines") {
	// The first two lines look like valid records but must be treated as
	// headers and skipped; only the third record is returned.
	std::vector<std::string> db = {
		"432000,0,HDR1,1,1,1",
		"432000,0,HDR2,1,1,1",
		"432000,0, REAL,10.0,5.0,35"
	};

	auto objs = read_deepsky(db, DeepSkySearch::WithinField,
	                         std::numbers::pi, 0.0, 0.1);
	REQUIRE(objs.size() == 1);
	CHECK(objs[0].name == "REAL");
}
