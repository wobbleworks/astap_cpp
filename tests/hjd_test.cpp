///----------------------------------------
///     @file hjd_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for pure-math helpers in src/core/hjd.cpp.
///  @details Covers airmass_calc (Pickering), atmospheric_absorption,
///           JD_to_HJD bound, equatorial -> galactic conversion, polar2
///           decomposition, and ra_az / az_ra round-trip. Functions that
///           depend on the `astap_main` runtime state (calc_jd_now,
///           calculate_az_alt, calculate_az_alt_basic) are intentionally
///           not exercised here.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "core/hjd.h"

#include <cmath>
#include <numbers>
#include <string>

using namespace astap::core;

///----------------------------------------
/// MARK: External-symbol stubs
///
/// hjd.cpp references several symbols that live in modules we don't want
/// to drag in for a pure-math unit test. We provide minimal stubs here so
/// the test binary links cleanly. The core-state globals (astap::head,
/// astap::jd_mid, etc.) come from the real globals.cpp, linked in via
/// CMakeLists.
///----------------------------------------

namespace astap::core {
	/// @brief No-op precession stub so sun() in JD_to_HJD links without
	///        pulling in core/fits.cpp's definition (and its cascade).
	void precession3(double /*jd0*/, double /*jd1*/,
	                 double& /*ra*/, double& /*dec*/) {}

	/// @brief No-op J2000->apparent stub (real impl in aberration.cpp).
	void J2000_to_apparent(double /*jd*/, double& /*ra*/, double& /*dec*/) {}

	/// @brief Stub for the RA/Dec text parser — not exercised here.
	void dec_text_to_radians(std::string /*inp*/,
	                         double& /*dec*/, bool& error_out) {
		error_out = true;
	}

	/// @brief Stub for the legacy memo2 log sink.
	void memo2_message(std::string_view /*msg*/) {}

	// Module-scope externs declared at the top of hjd.cpp. These live in
	// the astap::core namespace because hjd.cpp is wrapped in that
	// namespace when it declares them.
	double      site_lat_radians  = 999.0;
	double      site_long_radians = 999.0;
	double      wtime2actual      = 0.0;
	std::string sitelat;
	std::string sitelong;
	std::string lat_default;
	std::string long_default;
	std::string centaz;
	std::string centalt;
	double      focus_temp        = 100.0;
	double      pressure          = 1013.0;
}

namespace astap::stacking {
	/// @brief Stub returning J2000 epoch so calc_jd_now links.
	double julian_calc(int /*y*/, int /*m*/, double /*d*/,
	                   double /*h*/, double /*min*/, double /*s*/) {
		return 2451545.0;
	}
	/// @brief No-op date_to_jd stub.
	void date_to_jd(std::string_view /*date_obs*/,
	                std::string_view /*date_avg*/,
	                double /*exposure*/) {}
}

///----------------------------------------
/// MARK: airmass_calc (Pickering 2002)
///----------------------------------------

TEST_CASE("airmass_calc zenith is ~1") {
	CHECK(airmass_calc(90.0) == doctest::Approx(1.0).epsilon(1e-3));
}

TEST_CASE("airmass_calc 30 degrees is ~2") {
	CHECK(airmass_calc(30.0) == doctest::Approx(2.0).epsilon(1e-2));
}

TEST_CASE("airmass_calc just above horizon is large") {
	// Pickering at 1° altitude ≈ 26.3; at 0.5° ≈ 32; never more than ~38 at
	// non-trivial altitudes above the horizon threshold.
	CHECK(airmass_calc(1.0)   > 20.0);
	CHECK(airmass_calc(1.0)   < 50.0);
}

TEST_CASE("airmass_calc returns 999 sentinel at or below horizon") {
	// The Pascal implementation treats h < 1e-7 (including exact 0) as the
	// below-horizon sentinel.
	CHECK(airmass_calc(0.0)   == doctest::Approx(999.0));
	CHECK(airmass_calc(-5.0)  == doctest::Approx(999.0));
	CHECK(airmass_calc(-10.0) == doctest::Approx(999.0));
}

TEST_CASE("airmass_calc is monotonically decreasing with altitude") {
	double prev = airmass_calc(1.0);
	for (double alt = 5.0; alt <= 90.0; alt += 5.0) {
		const double am = airmass_calc(alt);
		CAPTURE(alt);
		CHECK(am < prev);   // strictly decreasing
		prev = am;
	}
}

///----------------------------------------
/// MARK: atmospheric_absorption
///----------------------------------------

TEST_CASE("atmospheric_absorption is 0.2811 mag per airmass") {
	// Schaefer 1992 / ICQ convention: 0.2811 mag/airmass.
	CHECK(atmospheric_absorption(1.0) == doctest::Approx(0.2811).epsilon(1e-3));
}

TEST_CASE("atmospheric_absorption scales linearly with airmass") {
	const double at_1 = atmospheric_absorption(1.0);
	CHECK(atmospheric_absorption(2.0) == doctest::Approx(2.0 * at_1));
	CHECK(atmospheric_absorption(3.5) == doctest::Approx(3.5 * at_1));
}

///----------------------------------------
/// MARK: JD_to_HJD
///----------------------------------------

TEST_CASE("JD_to_HJD correction is bounded by Earth-Sun light-time") {
	// Maximum correction is one Earth-Sun light-time (499 seconds) one-way,
	// i.e. 499 / 86400 ≈ 5.78e-3 days. Use 600 s bound for safety margin.
	constexpr double bound_days = 600.0 / 86400.0;

	constexpr double jd = 2460400.0;   // arbitrary modern date.
	for (double ra : {0.0, 1.0, 3.0, 5.0, 6.0}) {
		for (double dec : {-1.0, -0.5, 0.0, 0.5, 1.0}) {
			const double hjd = JD_to_HJD(jd, ra, dec);
			CAPTURE(ra);
			CAPTURE(dec);
			CHECK(std::abs(hjd - jd) < bound_days);
		}
	}
}

///----------------------------------------
/// MARK: equ_gal
///----------------------------------------

TEST_CASE("equ_gal North galactic pole has b ~= pi/2") {
	// NGP (J2000): RA 12h51m26.282s, Dec +27°07'42.01".
	//   RA  = 192.859508° = 3.366031 rad
	//   Dec = +27.128336° = +0.473478 rad
	constexpr double ra_ngp  = 3.366031;
	constexpr double dec_ngp = 0.473478;

	double l = 0, b = 0;
	equ_gal(ra_ngp, dec_ngp, l, b);
	// b should be π/2 within a few arcseconds.
	CHECK(b == doctest::Approx(std::numbers::pi / 2).epsilon(1e-3));
}

TEST_CASE("equ_gal Sagittarius A* is near galactic centre") {
	// Sgr A* (J2000): RA 17h45m40.04s, Dec -29°00'28.1".
	//   RA  ≈ 266.41670° = 4.649852 rad
	//   Dec ≈ -29.00781° = -0.506237 rad
	// Historical definition places Sgr A* a few arc-minutes off (l=0,b=0).
	constexpr double ra_sgra  = 4.649852;
	constexpr double dec_sgra = -0.506237;

	double l = 0, b = 0;
	equ_gal(ra_sgra, dec_sgra, l, b);
	// Expect l and b within ~0.01 rad (~34 arcmin) of (0, 0) — comfortable
	// margin over the known few-arcminute offset between Sgr A* and the
	// nominal galactic origin.
	CHECK(std::abs(l) < 0.01);
	CHECK(std::abs(b) < 0.01);
}

///----------------------------------------
/// MARK: polar2
///----------------------------------------

TEST_CASE("polar2 unit axis vectors") {
	double r = 0, theta = 0, phi = 0;

	SUBCASE("+X axis") {
		polar2(1, 0, 0, r, theta, phi);
		CHECK(r     == doctest::Approx(1.0));
		CHECK(theta == doctest::Approx(0.0));
		CHECK(phi   == doctest::Approx(0.0));
	}

	SUBCASE("+Y axis") {
		polar2(0, 1, 0, r, theta, phi);
		CHECK(r     == doctest::Approx(1.0));
		CHECK(theta == doctest::Approx(0.0));
		CHECK(phi   == doctest::Approx(std::numbers::pi / 2));
	}

	SUBCASE("+Z axis") {
		polar2(0, 0, 1, r, theta, phi);
		CHECK(r     == doctest::Approx(1.0));
		CHECK(theta == doctest::Approx(std::numbers::pi / 2));
	}
}

///----------------------------------------
/// MARK: ra_az <-> az_ra round-trip
///----------------------------------------

TEST_CASE("ra_az and az_ra are inverses") {
	// Observer at 37° N, sidereal time 2 rad, looking at various points.
	constexpr double lat  = 37.0 * std::numbers::pi / 180.0;
	constexpr double lon  = -122.0 * std::numbers::pi / 180.0;
	constexpr double t    = 2.0;

	struct Pos { double ra; double dec; };
	const Pos samples[] = {
		{0.5, 0.3},
		{2.1, -0.2},
		{4.7, 0.0},
		{1.0, 1.0},
	};

	for (const auto& p : samples) {
		double az = 0, alt = 0;
		ra_az(p.ra, p.dec, lat, lon, t, az, alt);

		double ra_back = 0, dec_back = 0;
		az_ra(az, alt, lat, lon, t, ra_back, dec_back);

		// Wrap the RA difference into [-π, π] to handle the 0 / 2π seam.
		const double two_pi = 2.0 * std::numbers::pi;
		double d_ra = ra_back - p.ra;
		while (d_ra >  std::numbers::pi) d_ra -= two_pi;
		while (d_ra < -std::numbers::pi) d_ra += two_pi;

		CAPTURE(p.ra);
		CAPTURE(p.dec);
		CHECK(std::abs(d_ra)            < 1e-9);
		CHECK(std::abs(dec_back - p.dec) < 1e-9);
	}
}
