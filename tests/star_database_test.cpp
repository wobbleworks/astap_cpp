///----------------------------------------
///     @file star_database_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for the HNSKY star-database area-resolution logic
///           (src/reference/star_database.cpp).
///  @details Exercises `find_areas` which maps (RA, Dec, FOV) to up to four
///           1-based area indices covering the field. I/O-bound reading of
///           actual .290 files is intentionally out of scope here.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "reference/star_database.h"

#include <cmath>
#include <numbers>
#include <string>

using namespace astap::reference;

///----------------------------------------
/// MARK: Constants
///----------------------------------------

static constexpr double kPi = std::numbers::pi;
static constexpr double kDeg2Rad = kPi / 180.0;

///----------------------------------------
/// MARK: find_areas — 290 database
///----------------------------------------

TEST_CASE("find_areas: south celestial pole resolves to area 1 (290-format)") {
	database_type = 290;
	int a1 = 0, a2 = 0, a3 = 0, a4 = 0;
	double f1 = 0, f2 = 0, f3 = 0, f4 = 0;

	// Small FOV at the south pole. Area 1 is the polar cap file 0101.290.
	find_areas(0.0, -kPi / 2, 1.0 * kDeg2Rad, a1, a2, a3, a4, f1, f2, f3, f4);
	CHECK(a1 == 1);
	CHECK(f1 > 0.0);
}

TEST_CASE("find_areas: north celestial pole resolves to area 290 (290-format)") {
	database_type = 290;
	int a1 = 0, a2 = 0, a3 = 0, a4 = 0;
	double f1 = 0, f2 = 0, f3 = 0, f4 = 0;

	find_areas(0.0, kPi / 2, 1.0 * kDeg2Rad, a1, a2, a3, a4, f1, f2, f3, f4);
	CHECK(a1 == 290);
	CHECK(f1 > 0.0);
}

TEST_CASE("find_areas: equator resolves to an area in range [1, 290]") {
	database_type = 290;
	int a1 = 0, a2 = 0, a3 = 0, a4 = 0;
	double f1 = 0, f2 = 0, f3 = 0, f4 = 0;

	// Small FOV near the equator.
	find_areas(1.0, 0.0, 1.0 * kDeg2Rad, a1, a2, a3, a4, f1, f2, f3, f4);
	CHECK(a1 >= 1);
	CHECK(a1 <= 290);
	CHECK(f1 > 0.0);
}

TEST_CASE("find_areas: large FOV near equator spans multiple areas") {
	database_type = 290;
	int a1 = 0, a2 = 0, a3 = 0, a4 = 0;
	double f1 = 0, f2 = 0, f3 = 0, f4 = 0;

	// 8° FOV — spans multiple ~15° cells near the equator.
	find_areas(1.0, 0.0, 8.0 * kDeg2Rad, a1, a2, a3, a4, f1, f2, f3, f4);
	CHECK(a1 >= 1);
	CHECK(f1 > 0.0);
	// Expect at least 2 non-zero areas covering the field.
	const int used = (a1 != 0) + (a2 != 0) + (a3 != 0) + (a4 != 0);
	CHECK(used >= 2);
}

TEST_CASE("find_areas: fractions of used areas sum close to 1.0") {
	database_type = 290;
	int a1 = 0, a2 = 0, a3 = 0, a4 = 0;
	double f1 = 0, f2 = 0, f3 = 0, f4 = 0;

	find_areas(2.0, 0.3, 4.0 * kDeg2Rad, a1, a2, a3, a4, f1, f2, f3, f4);

	// Only the fractions corresponding to non-zero area indices contribute.
	double total = 0.0;
	if (a1 != 0) total += f1;
	if (a2 != 0) total += f2;
	if (a3 != 0) total += f3;
	if (a4 != 0) total += f4;

	// Sum should be close to 1.0 (the full field). Allow 10% slack since
	// fractional coverage for the four quadrants is an approximation.
	CHECK(total == doctest::Approx(1.0).epsilon(0.1));
}

///----------------------------------------
/// MARK: find_areas — 1476 database
///----------------------------------------

TEST_CASE("find_areas: south pole resolves to area 1 (1476-format)") {
	database_type = 1476;
	int a1 = 0, a2 = 0, a3 = 0, a4 = 0;
	double f1 = 0, f2 = 0, f3 = 0, f4 = 0;

	find_areas(0.0, -kPi / 2, 1.0 * kDeg2Rad, a1, a2, a3, a4, f1, f2, f3, f4);
	CHECK(a1 == 1);
	CHECK(f1 > 0.0);
}

TEST_CASE("find_areas: north pole resolves to area 1476 (1476-format)") {
	database_type = 1476;
	int a1 = 0, a2 = 0, a3 = 0, a4 = 0;
	double f1 = 0, f2 = 0, f3 = 0, f4 = 0;

	find_areas(0.0, kPi / 2, 1.0 * kDeg2Rad, a1, a2, a3, a4, f1, f2, f3, f4);
	CHECK(a1 == 1476);
	CHECK(f1 > 0.0);
}

TEST_CASE("find_areas: equator area number in [1, 1476] for 1476-format") {
	database_type = 1476;
	int a1 = 0, a2 = 0, a3 = 0, a4 = 0;
	double f1 = 0, f2 = 0, f3 = 0, f4 = 0;

	find_areas(3.0, 0.0, 0.5 * kDeg2Rad, a1, a2, a3, a4, f1, f2, f3, f4);
	CHECK(a1 >= 1);
	CHECK(a1 <= 1476);
	CHECK(f1 > 0.0);
}

///----------------------------------------
/// MARK: RA wraparound stability
///----------------------------------------

TEST_CASE("find_areas: ra = 0 and ra = 2pi should pick the same primary area") {
	database_type = 290;
	int a1_l = 0, a2_l = 0, a3_l = 0, a4_l = 0;
	int a1_r = 0, a2_r = 0, a3_r = 0, a4_r = 0;
	double f1_l = 0, f2_l = 0, f3_l = 0, f4_l = 0;
	double f1_r = 0, f2_r = 0, f3_r = 0, f4_r = 0;

	find_areas(0.0,      0.0, 1.0 * kDeg2Rad, a1_l, a2_l, a3_l, a4_l, f1_l, f2_l, f3_l, f4_l);
	find_areas(2 * kPi,  0.0, 1.0 * kDeg2Rad, a1_r, a2_r, a3_r, a4_r, f1_r, f2_r, f3_r, f4_r);

	CHECK(a1_l == a1_r);
}

///----------------------------------------
/// MARK: get_database_passband — auto / Local mode (filter-driven)
///----------------------------------------

TEST_CASE("get_database_passband: auto mode infers passband from FITS filter") {
	SUBCASE("empty filter defaults to BP") {
		CHECK(get_database_passband("",       "auto") == "BP");
		CHECK(get_database_passband("",       "Local") == "BP");
	}

	SUBCASE("CV (clear-view) filter resolves to BP") {
		CHECK(get_database_passband("CV",     "auto") == "BP");
		CHECK(get_database_passband("cv",     "auto") == "BP"); // case-insensitive
	}

	SUBCASE("Sloan filters resolve to SG/SR/SI") {
		CHECK(get_database_passband("Sloan g", "auto") == "SG");
		CHECK(get_database_passband("Sloan r", "auto") == "SR");
		CHECK(get_database_passband("Sloan i", "auto") == "SI");
		// Plain 'S' with no colour letter falls back to BP.
		CHECK(get_database_passband("S",       "auto") == "BP");
	}

	SUBCASE("Johnson/Cousins G and V both map to Johnson-V") {
		CHECK(get_database_passband("G",       "auto") == "V");
		CHECK(get_database_passband("V",       "auto") == "V");
	}

	SUBCASE("Johnson-B and Cousins-R") {
		CHECK(get_database_passband("B",       "auto") == "B");
		CHECK(get_database_passband("R",       "auto") == "R");
	}

	SUBCASE("unknown filter name falls back to BP") {
		CHECK(get_database_passband("Xyzzy",   "auto") == "BP");
	}
}

///----------------------------------------
/// MARK: get_database_passband — manual mode (database-driven)
///----------------------------------------

TEST_CASE("get_database_passband: manual mode parses passband from database name") {
	SUBCASE("BP is checked before B") {
		// "Gaia BP" contains 'B' but we must pick "BP", not "B".
		CHECK(get_database_passband("V", "Gaia BP")    == "BP");
	}

	SUBCASE("V / B / R all recognised") {
		CHECK(get_database_passband("V", "Johnson V")   == "V");
		CHECK(get_database_passband("R", "Johnson B")   == "B");
		CHECK(get_database_passband("G", "Cousins R")   == "R");
	}

	SUBCASE("Sloan SG/SR/SI: Pascal priority matches B/R first") {
		// The original Pascal checks plain B and R before SG/SR/SI in manual
		// mode, so "Sloan SG" returns "B" (via the 'B' in "Sloan"), not "SG".
		// Preserved here intentionally; a database name like "SG_only" that
		// lacks B/V/R can still reach the Sloan branches.
		CHECK(get_database_passband("", "SG_only")      == "SG");
		CHECK(get_database_passband("", "SI_only")      == "SI");
	}

	SUBCASE("unrecognised database name returns double-question mark") {
		CHECK(get_database_passband("V", "unknown_db")  == "??");
	}
}
