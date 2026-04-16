///----------------------------------------
///     @file constellations_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for constellation data tables (src/analysis/constellations.cpp).
///  @details Validates array sizes, known entries, drawing-mode values,
///           coordinate ranges, and the Serpens Caput/Cauda duplication.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright (C) 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "analysis/constellations.h"

#include <algorithm>
#include <string_view>

using namespace astap::analysis;

///----------------------------------------
/// MARK: Array sizes
///----------------------------------------

TEST_CASE("constellation data: array sizes") {
	CHECK(constellation.size() == 603);
	CHECK(constPos.size() == 89);
	CHECK(constShortName.size() == 89);
}

///----------------------------------------
/// MARK: First entry is Andromeda alpha
///----------------------------------------

TEST_CASE("constellation data: first entry is Andromeda alpha") {
	CHECK(constellation[0].dm == -2);
	// Bayer designation starts with Greek alpha (UTF-8: 0xCE 0xB1).
	std::string_view bay = constellation[0].bay;
	REQUIRE(bay.size() >= 2);
	CHECK(static_cast<unsigned char>(bay[0]) == 0xCE);
	CHECK(static_cast<unsigned char>(bay[1]) == 0xB1);
}

///----------------------------------------
/// MARK: Known constellation short names
///----------------------------------------

TEST_CASE("constellation data: known short names") {
	CHECK(constShortName[0] == "And");
	CHECK(constShortName[59] == "Ori");
	CHECK(constShortName[88] == "Vul");
}

///----------------------------------------
/// MARK: All entries have valid dm
///----------------------------------------

TEST_CASE("constellation data: all entries have valid dm") {
	for (size_t i = 0; i < constellation.size(); ++i) {
		bool valid = (constellation[i].dm == kDrawModeStart ||
		              constellation[i].dm == kDrawModeDraw);
		CHECK_MESSAGE(valid, "invalid dm at index " << i);
	}
}

///----------------------------------------
/// MARK: RA/Dec ranges
///----------------------------------------

TEST_CASE("constellation data: RA and Dec ranges") {
	for (size_t i = 0; i < constellation.size(); ++i) {
		CHECK_MESSAGE(constellation[i].ra >= 0,
		              "RA < 0 at index " << i);
		CHECK_MESSAGE(constellation[i].ra <= 24000,
		              "RA > 24000 at index " << i);
		CHECK_MESSAGE(constellation[i].dec >= -9000,
		              "Dec < -9000 at index " << i);
		CHECK_MESSAGE(constellation[i].dec <= 9000,
		              "Dec > 9000 at index " << i);
	}
}

///----------------------------------------
/// MARK: Serpens appears twice
///----------------------------------------

TEST_CASE("constellation data: Serpens appears twice") {
	int count = 0;
	for (size_t i = 0; i < constShortName.size(); ++i) {
		if (constShortName[i] == "Ser") ++count;
	}
	CHECK(count == 2);
	// Serpens Caput and Cauda are at indices 75 and 76.
	CHECK(constShortName[75] == "Ser");
	CHECK(constShortName[76] == "Ser");
}
