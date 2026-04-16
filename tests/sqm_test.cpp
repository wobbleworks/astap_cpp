///----------------------------------------
///     @file sqm_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for the sky-quality classification helpers
///           (src/core/sqm.cpp).
///  @details Only the pure-math `bortle` classifier is exercised here.
///           `calculate_sqm` needs a full photometry/airmass pipeline
///           that is out of scope for a self-contained unit test.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "core/sqm.h"

#include <array>
#include <atomic>
#include <cstdint>
#include <string>
#include <vector>

using astap::core::bortle;

///----------------------------------------
/// MARK: External symbols sqm.cpp references at link time
///----------------------------------------

// link_stubs.cpp + globals.cpp + hjd.cpp + imaging.cpp satisfy most of
// what calculate_sqm pulls in. A few symbols from not-yet-linked modules
// (aberration, wcs, stacking) need direct stubs here.

namespace astap::core {
	void J2000_to_apparent(double /*jd*/, double& /*ra*/, double& /*dec*/) {}
	void dec_text_to_radians(std::string /*inp*/, double& /*dec*/, bool& err) {
		err = true;
	}
	void precession3(double, double, double&, double&) {}
	void old_to_new_WCS(astap::Header& /*head*/) {}

	// The signature of the get_background overload calculate_sqm wants is
	// the imaging.cpp histogram-aware version. Real impl isn't linked in.
	void get_background(int, const astap::ImageArray&, bool, bool,
	                    astap::Background& bck,
	                    int (&)[3][65536], int (&)[3],
	                    int, int,
	                    const std::string&,
	                    std::vector<std::string>&) {
		bck = astap::Background{};
	}
}

namespace astap::stacking {
	double julian_calc(int, int, double, double, double, double) {
		return 2451545.0;
	}
	void date_to_jd(std::string_view, std::string_view, double) {}
}

///----------------------------------------
/// MARK: bortle — classification boundaries
///
/// The Bortle scale classes, in ascending sky brightness (darker to
/// brighter):
///   Class 1: excellent dark-sky site        (~21.75-22.0 mag/arcsec^2)
///   Class 2: typical truly dark             (~21.5)
///   Class 3: rural                          (~21.0-21.3)
///   Class 4: rural/suburban transition      (~20.5-21.0)
///   Class 5: suburban                       (~19.5-20.5)
///   Class 6: bright suburban                (~19.0-19.5)
///   Class 7: suburban/urban transition      (~18.5-19.0)
///   Class 8-9: city                         (< 18.5)
///----------------------------------------

TEST_CASE("bortle returns a non-empty description for any plausible SQM") {
	for (double sqm = 16.0; sqm <= 23.0; sqm += 0.5) {
		CAPTURE(sqm);
		CHECK_FALSE(bortle(sqm).empty());
	}
}

TEST_CASE("bortle for a dark site (SQM >= 21.5) differs from city (SQM <= 18)") {
	const auto dark = bortle(21.8);
	const auto city = bortle(17.5);

	CAPTURE(dark);
	CAPTURE(city);
	CHECK(dark != city);
}

TEST_CASE("bortle is monotonic in SQM (darker sky -> lower class)") {
	// Bortle scale assigns a numeric class 1-9 (1 = best). Higher SQM
	// (darker) should map to a lower class than lower SQM (brighter).
	// We can't assume the exact string format the Pascal uses, so check
	// the numeric digit embedded in the description.
	const auto dark_desc   = bortle(22.0);
	const auto bright_desc = bortle(17.0);
	CAPTURE(dark_desc);
	CAPTURE(bright_desc);

	auto first_digit = [](std::string_view s) -> int {
		for (char c : s) {
			if (c >= '0' && c <= '9') return c - '0';
		}
		return -1;
	};

	const int dark_class   = first_digit(dark_desc);
	const int bright_class = first_digit(bright_desc);
	CHECK(dark_class   > 0);
	CHECK(bright_class > 0);
	CHECK(dark_class   < bright_class);
}

TEST_CASE("bortle boundary consistency: tiny SQM steps don't jump multiple classes") {
	auto first_digit = [](std::string_view s) -> int {
		for (char c : s) {
			if (c >= '0' && c <= '9') return c - '0';
		}
		return -1;
	};

	int prev_class = first_digit(bortle(16.0));
	CHECK(prev_class > 0);

	for (double sqm = 16.1; sqm <= 22.5; sqm += 0.1) {
		const int cls = first_digit(bortle(sqm));
		CAPTURE(sqm);
		CAPTURE(cls);
		CHECK(cls > 0);
		// Class number must either stay the same or drop by exactly 1
		// as sky gets darker. Never jump multiple classes per 0.1 step.
		CHECK(prev_class - cls <= 1);
		CHECK(prev_class - cls >= 0);
		prev_class = cls;
	}
}
