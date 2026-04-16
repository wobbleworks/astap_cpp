///----------------------------------------
///     @file asteroid_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for asteroid/comet helpers (src/analysis/asteroid.cpp).
///  @details Covers deltaT_calc, parallax_xyz, illum2, asteroid_magn_comp,
///           and convert_MPCORB_line. Link-isolation stubs are provided for
///           julian_calc and the ephem::{earth_state, propagate,
///           cartesian_to_spherical} symbols that asteroid.cpp depends on.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright (C) 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "analysis/asteroid.h"

#include <cmath>
#include <numbers>
#include <string>

// ---- Link-isolation stubs --------------------------------------------------
// These satisfy unresolved symbols from asteroid.cpp without pulling in the
// full stacking/core modules.

namespace astap::stacking {
[[nodiscard]] double julian_calc(int /*yyyy*/, int /*mm*/, double /*dd*/,
                                 double /*hours*/, double /*minutes*/,
                                 double /*seconds*/) {
	return 2451545.0;  // J2000 placeholder
}
} // namespace astap::stacking

namespace astap::core::ephem {
State earth_state(double /*jd_tt*/, ReferenceFrame /*frame*/) {
	return { .position = {1.0, 0.0, 0.0}, .velocity = {0.0, 0.0, 0.0} };
}
State propagate(const OrbitalElements& /*elements*/, double /*jd_tt*/) {
	return { .position = {2.0, 0.0, 0.0}, .velocity = {0.0, 0.0, 0.0} };
}
SphericalResult cartesian_to_spherical(const Vec3& xyz) noexcept {
	SphericalResult r;
	r.radius = std::sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
	r.dec    = std::asin(xyz[2] / (r.radius + 1e-30));
	r.ra     = std::atan2(xyz[1], xyz[0]);
	if (r.ra < 0.0) r.ra += 2.0 * std::numbers::pi;
	return r;
}
} // namespace astap::core::ephem

using namespace astap::analysis;

///----------------------------------------
/// MARK: deltaT_calc
///----------------------------------------

TEST_CASE("deltaT_calc: known epochs") {
	SUBCASE("J2020.0") {
		// JD 2458849.5 -> year ~2020, expect ~69.5 s as days.
		double dt = deltaT_calc(2458849.5);
		double dt_sec = dt * 86400.0;
		CHECK(dt_sec == doctest::Approx(69.5).epsilon(1.0));
	}

	SUBCASE("J2025.0") {
		// JD 2460676.5 -> year ~2025, expect 69..80 s.
		double dt = deltaT_calc(2460676.5);
		double dt_sec = dt * 86400.0;
		CHECK(dt_sec > 69.0);
		CHECK(dt_sec < 80.0);
	}

	SUBCASE("far outside range — year 1500") {
		// JD ~2268924 for year 1500 — returns 60s as days.
		double dt = deltaT_calc(2268924.0);
		double dt_sec = dt * 86400.0;
		CHECK(dt_sec == doctest::Approx(60.0).epsilon(0.01));
	}
}

///----------------------------------------
/// MARK: parallax_xyz
///----------------------------------------

TEST_CASE("parallax_xyz: equatorial observer correction magnitude") {
	double x = 1.0, y = 0.0, z = 0.0;
	double x0 = x, y0 = y, z0 = z;

	parallax_xyz(0.0, 0.0, x, y, z);  // wtime=0, lat=0 (equator)

	double dx = x - x0;
	double dy = y - y0;
	double dz = z - z0;
	double correction = std::sqrt(dx * dx + dy * dy + dz * dz);

	// Earth radius / AU ~ 4.3e-5 AU.
	CHECK(correction > 1e-5);
	CHECK(correction < 1e-4);
}

///----------------------------------------
/// MARK: illum2
///----------------------------------------

TEST_CASE("illum2: opposition geometry") {
	// Planet at (2,0,0), Earth at (1,0,0).
	double r_sp = 0, r_ep = 0, elong = 0, phi = 0, phase = 0;
	illum2(2.0, 0.0, 0.0,  1.0, 0.0, 0.0,
	       r_sp, r_ep, elong, phi, phase);

	CHECK(r_sp == doctest::Approx(2.0));
	CHECK(r_ep == doctest::Approx(1.0));
	CHECK(elong == doctest::Approx(180.0).epsilon(0.1));
	CHECK(phi == doctest::Approx(0.0).epsilon(0.1));
	CHECK(phase == doctest::Approx(100.0).epsilon(0.1));
}

TEST_CASE("illum2: quadrature geometry") {
	// Planet at (0,2,0), Earth at (1,0,0).
	double r_sp = 0, r_ep = 0, elong = 0, phi = 0, phase = 0;
	illum2(0.0, 2.0, 0.0,  1.0, 0.0, 0.0,
	       r_sp, r_ep, elong, phi, phase);

	CHECK(r_sp == doctest::Approx(2.0));
	CHECK(r_ep > 1.0);
	CHECK(elong > 0.0);
	CHECK(elong < 180.0);
	CHECK(phase >= 0.0);
	CHECK(phase <= 100.0);
}

///----------------------------------------
/// MARK: asteroid_magn_comp
///----------------------------------------

TEST_CASE("asteroid_magn_comp: zero phase angle") {
	// b=0 -> fully illuminated, correction should be ~0.
	double corr = asteroid_magn_comp(0.15, 0.0);
	CHECK(corr == doctest::Approx(0.0).epsilon(0.05));
}

TEST_CASE("asteroid_magn_comp: moderate phase angle") {
	// b=0.35 rad (~20 deg), correction should be positive (fainter).
	double corr = asteroid_magn_comp(0.15, 0.35);
	CHECK(corr > 0.0);
}

///----------------------------------------
/// MARK: convert_MPCORB_line
///----------------------------------------

TEST_CASE("convert_MPCORB_line: Ceres") {
	// Actual Ceres line from MPCORB format, padded to >= 195 chars.
	std::string line =
		"00001    3.4   0.15 K205V 162.68631   73.73161   80.28698   10.58862"
		"  0.0775571  0.21406009   2.7676569  0 MPO492748  6751 115 1801-2019"
		" 0.60 M-v 30h Williams   0000      (1) Ceres              20190915";
	// Pad to at least 200 chars.
	while (line.size() < 200) line += ' ';

	std::string desn, name;
	int yy = 0, mm = 0;
	double dd = 0, a_e = 0, a_a = 0, a_i = 0, a_ohm = 0, a_w = 0, a_M = 0;
	double h = 0, g = 0;

	convert_MPCORB_line(line, desn, name, yy, mm, dd, a_e, a_a, a_i, a_ohm, a_w, a_M, h, g);

	CHECK(desn == "00001");
	CHECK(h == doctest::Approx(3.4).epsilon(0.01));
	CHECK(g == doctest::Approx(0.15).epsilon(0.01));
	CHECK(a_a == doctest::Approx(2.7676569).epsilon(0.001));
	CHECK(a_i == doctest::Approx(10.58862).epsilon(0.01));
}

TEST_CASE("convert_MPCORB_line: line too short") {
	std::string short_line = "00001    3.4   0.15 K205V";

	std::string desn, name;
	int yy = 0, mm = 0;
	double dd = 0, a_e = 0, a_a = 0, a_i = 0, a_ohm = 0, a_w = 0, a_M = 0;
	double h = 0, g = 0;

	convert_MPCORB_line(short_line, desn, name, yy, mm, dd, a_e, a_a, a_i, a_ohm, a_w, a_M, h, g);

	CHECK(desn.empty());
}
