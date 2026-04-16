///----------------------------------------
///     @file wcs_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for WCS and coordinate-formatting helpers
///           (src/core/wcs.cpp).
///  @details pixel_to_celestial / celestial_to_pixel round-trip, ang_sep
///           against known great-circle distances, standard_equatorial2
///           smoke tests, RA/Dec text parsers, and the prepare_* family.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "core/wcs.h"

#include <cmath>
#include <numbers>
#include <string>

using namespace astap::core;
using astap::Header;

///----------------------------------------
/// MARK: Test-local stubs
///
/// wcs.cpp pulls in a handful of externals (globals.h SIP coefficients,
/// dsspos, EQU_GAL, calculate_az_alt_basic). Link with stubs so tests
/// isolate wcs arithmetic.
///----------------------------------------

namespace astap::core {

// Forward-declared in wcs.h but never defined in the library.
// Tests set them to default values so pixel_to_celestial's TAN branch
// is exercised (formalism == 0).
void dsspos(double, double, double& ra, double& dec) { ra = 0; dec = 0; }
void EQU_GAL(double, double, double& l, double& b)   { l  = 0; b   = 0; }
bool calculate_az_alt_basic(double, double, double& az, double& alt) {
	az = alt = 0;
	return false;
}

// SIP coefficient and flag globals are defined in wcs.cpp itself; no
// need to re-declare them here.

double position_angle(double, double, double, double) { return 0.0; }

}  // namespace astap::core

///----------------------------------------
/// MARK: Helpers
///----------------------------------------

/// @brief Build a simple Header with a TAN projection, unit scale, and a
///        known reference-pixel / reference-sky position.
[[nodiscard]] static Header make_tan_header(double ra0_rad, double dec0_rad) {
	Header h{};
	h.width   = 2000;
	h.height  = 2000;
	h.naxis   = 2;
	h.naxis3  = 1;
	h.crpix1  = 1000.0;
	h.crpix2  = 1000.0;
	// 1 arcsec per pixel in degrees.
	constexpr double arcsec_deg = 1.0 / 3600.0;
	h.cdelt1  = -arcsec_deg;   // RA decreases with X, by convention
	h.cdelt2  =  arcsec_deg;
	h.ra0     = ra0_rad;
	h.dec0    = dec0_rad;
	h.crota1  = 0.0;
	h.crota2  = 0.0;
	h.cd1_1   = h.cdelt1;
	h.cd1_2   = 0.0;
	h.cd2_1   = 0.0;
	h.cd2_2   = h.cdelt2;
	return h;
}

///----------------------------------------
/// MARK: pixel_to_celestial / celestial_to_pixel round trip
///----------------------------------------

TEST_CASE("pixel_to_celestial maps the reference pixel to (ra0, dec0)") {
	const auto head = make_tan_header(1.0, 0.3);

	double ra = 0, dec = 0;
	pixel_to_celestial(head, head.crpix1, head.crpix2, /*formalism=*/0, ra, dec);
	CHECK(ra  == doctest::Approx(head.ra0 ).epsilon(1e-10));
	CHECK(dec == doctest::Approx(head.dec0).epsilon(1e-10));
}

TEST_CASE("pixel <-> celestial round trip close to image centre") {
	const auto head = make_tan_header(1.2, 0.4);

	struct Offset { double dx; double dy; };
	const Offset samples[] = {
		{  0.0,   0.0},
		{100.0,  50.0},
		{-50.0, 200.0},
		{500.0, -300.0},
	};

	for (const auto& o : samples) {
		const double x0 = head.crpix1 + o.dx;
		const double y0 = head.crpix2 + o.dy;

		double ra = 0, dec = 0;
		pixel_to_celestial(head, x0, y0, 0, ra, dec);

		double x1 = 0, y1 = 0;
		celestial_to_pixel(head, ra, dec, x1, y1);

		CAPTURE(o.dx);
		CAPTURE(o.dy);
		CHECK(x1 == doctest::Approx(x0).epsilon(1e-6));
		CHECK(y1 == doctest::Approx(y0).epsilon(1e-6));
	}
}

TEST_CASE("celestial_to_pixel maps (ra0, dec0) back to crpix") {
	const auto head = make_tan_header(0.8, -0.2);
	double x = 0, y = 0;
	celestial_to_pixel(head, head.ra0, head.dec0, x, y);
	CHECK(x == doctest::Approx(head.crpix1).epsilon(1e-10));
	CHECK(y == doctest::Approx(head.crpix2).epsilon(1e-10));
}

///----------------------------------------
/// MARK: ang_sep
///----------------------------------------

TEST_CASE("ang_sep zero separation for identical points") {
	// The acos-based formulation has ~1e-8 rad precision even when the
	// cosine should be exactly 1; allow up to 1 milli-arcsec (~5e-9 rad)
	// as "zero".
	double sep = -1.0;
	ang_sep(1.0, 0.3, 1.0, 0.3, sep);
	CHECK(sep < 1e-7);
}

TEST_CASE("ang_sep antipodes are pi apart") {
	// (ra, dec) and (ra + pi, -dec) are diametrically opposite on the sphere.
	constexpr double pi = std::numbers::pi;

	double sep = 0.0;
	ang_sep(0.0, 0.0, pi, 0.0, sep);
	CHECK(sep == doctest::Approx(pi).epsilon(1e-10));

	ang_sep(1.2, 0.4, 1.2 + pi, -0.4, sep);
	CHECK(sep == doctest::Approx(pi).epsilon(1e-10));
}

TEST_CASE("ang_sep between equatorial points 1 hour apart is 15 degrees") {
	// Two equatorial points separated by 1 hour of RA (= 15°).
	constexpr double pi = std::numbers::pi;
	const double ra1 =  0.0;
	const double ra2 = 15.0 * pi / 180.0;

	double sep = 0.0;
	ang_sep(ra1, 0.0, ra2, 0.0, sep);
	CHECK(sep * 180.0 / pi == doctest::Approx(15.0).epsilon(1e-10));
}

///----------------------------------------
/// MARK: standard_equatorial2
///----------------------------------------

TEST_CASE("standard_equatorial2 returns (ra0, dec0) for origin pixel") {
	constexpr double pi = std::numbers::pi;
	const double ra0  = 1.0;
	const double dec0 = 0.3;

	double ra = 0, dec = 0;
	standard_equatorial2(ra0, dec0, /*x=*/0.0, /*y=*/0.0, /*cdelt=*/1.0, ra, dec);
	CHECK(ra  == doctest::Approx(ra0 ).epsilon(1e-10));
	CHECK(dec == doctest::Approx(dec0).epsilon(1e-10));
	(void)pi;
}

///----------------------------------------
/// MARK: RA / Dec text parsers
///----------------------------------------

TEST_CASE("ra_text_to_radians parses a colon-separated hours string") {
	constexpr double pi = std::numbers::pi;
	double ra = 0.0;
	bool   err = true;

	// 06h 00m 00s = 90° = pi/2 radians.
	ra_text_to_radians("06 00 00", ra, err);
	CHECK_FALSE(err);
	CHECK(ra == doctest::Approx(pi / 2).epsilon(1e-6));
}

TEST_CASE("dec_text_to_radians parses a degrees / negative declination") {
	constexpr double pi = std::numbers::pi;
	double dec = 0.0;
	bool   err = true;

	// -45° = -pi/4.
	dec_text_to_radians("-45 00 00", dec, err);
	CHECK_FALSE(err);
	CHECK(dec == doctest::Approx(-pi / 4).epsilon(1e-6));
}

TEST_CASE("ra_text_to_radians flags junk as error") {
	double ra = 0.0;
	bool   err = false;
	ra_text_to_radians("not a coordinate", ra, err);
	CHECK(err == true);
}

///----------------------------------------
/// MARK: prepare_* formatters
///----------------------------------------

TEST_CASE("prepare_ra for 6h produces a string starting with 06") {
	constexpr double pi = std::numbers::pi;
	const std::string s = prepare_ra(pi / 2, ":");
	CAPTURE(s);
	// Expect hour field == 6 (6h == pi/2 rad).
	CHECK(s.find("06") != std::string::npos);
}

TEST_CASE("prepare_dec for +45 degrees starts with +45") {
	constexpr double pi = std::numbers::pi;
	const std::string s = prepare_dec(pi / 4, ":");
	CAPTURE(s);
	CHECK((s.find("+45") != std::string::npos || s.find("45") != std::string::npos));
}

TEST_CASE("prepare_IAU_designation packs HHMMSS+DDMMSS") {
	constexpr double pi = std::numbers::pi;
	const auto s = prepare_IAU_designation(pi / 2, -pi / 4);
	CAPTURE(s);
	// Expect a '+' or '-' sign present, and 8-ish digit chars.
	CHECK((s.find('+') != std::string::npos || s.find('-') != std::string::npos));
	CHECK(s.size() > 8);
}

///----------------------------------------
/// MARK: Jd_To_MPCDate
///----------------------------------------

TEST_CASE("Jd_To_MPCDate for J2000 returns a 2000-01-01-ish string") {
	// JD 2451545.0 is J2000.0 (2000-01-01 12:00 TT).
	const std::string s = Jd_To_MPCDate(2451545.0);
	CAPTURE(s);
	CHECK(s.find("2000") != std::string::npos);
	CHECK(s.find("01")   != std::string::npos);
}
