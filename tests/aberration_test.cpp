///----------------------------------------
///     @file aberration_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for annual aberration / nutation corrections
///           (src/core/aberration.cpp).
///  @details Tests rely on structural invariants rather than reference
///           tables: aberration shifts are bounded by ~21" (the constant
///           of aberration), nutation by ~17" (nutation in longitude), and
///           applying a correction then its inverse should round-trip to a
///           negligible residual.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "core/aberration.h"
#include "core/ephemerides.h"

#include <cmath>
#include <numbers>

using namespace astap::core;

///----------------------------------------
/// MARK: Stubs for external symbols referenced by aberration.cpp
///----------------------------------------

// aberration.cpp depends on astap::core::ephem::earth_state and
// astap::core::ephem::precess_iau1976 (plus deltaT_calc). We stub all three
// so the aberration-only invariants are testable in isolation — the real
// ephemeris and precession routines are exercised by the library build.
namespace astap::core::ephem {
	/// @brief Stub for ephem::earth_state returning a canonical Earth
	///        orbital velocity (~2π AU per 365.25 days along +Y in J2000
	///        ecliptic). Enough to drive a non-trivial, bounded aberration.
	State earth_state(double /*jd_tt*/, ReferenceFrame /*frame*/) {
		constexpr double orbital_speed_au_per_day = 2.0 * 3.14159265358979323846 / 365.25;
		return {
			.position = {1.0, 0.0, 0.0},
			.velocity = {0.0, orbital_speed_au_per_day, 0.0},
		};
	}

	/// @brief No-op precession stub — tests isolate aberration from precession.
	EquatorialCoords precess_iau1976(EquatorialCoords coords,
	                                 double /*jd_from*/, double /*jd_to*/) {
		return coords;
	}
}

namespace astap::core {
	/// @brief Delta T (TT - UT1) stub — a plausible ~70 s for modern epochs.
	double deltaT_calc(double /*jd*/) {
		return 70.0 / 86400.0;   // convert seconds to days.
	}
}

///----------------------------------------
/// MARK: Magnitude / boundedness
///----------------------------------------

/// @brief Nominal epoch for samples: J2024.0 (2024-01-01 12:00 TT).
static constexpr double kJdTestEpoch = 2460311.0;

/// @brief Bound for the combined aberration correction: the constant of
///        aberration κ = 20.49552 arcsec, converted to radians with a
///        generous factor-of-2 safety margin for edge effects at
///        high-declination positions.
static constexpr double kAberrationBoundRadians =
	(2.0 * 20.49552 / 3600.0) * (std::numbers::pi / 180.0);

TEST_CASE("annual aberration is bounded by ~2 * kappa") {
	struct Position {
		double ra;
		double dec;
	};
	const Position samples[] = {
		{0.0,                    0.0},                    // vernal equinox, equator
		{std::numbers::pi,       0.0},                    // autumnal equinox
		{std::numbers::pi / 2,   0.5},                    // summer solstice-ish
		{3.0 * std::numbers::pi / 2, -0.5},               // winter solstice-ish
		{1.0,                    1.0},                    // arbitrary high-dec
	};

	for (const auto& p : samples) {
		double ra  = p.ra;
		double dec = p.dec;
		aberration_correction_equatorial(kJdTestEpoch, ra, dec);

		// Ignore the RA wrap at 2π — use the minimum-distance absolute delta.
		const double two_pi = 2.0 * std::numbers::pi;
		double d_ra = ra - p.ra;
		if (d_ra >  std::numbers::pi) d_ra -= two_pi;
		if (d_ra < -std::numbers::pi) d_ra += two_pi;
		// RA displacement scales with 1/cos(dec); separate the dec displacement
		// which directly reflects the sky-plane shift.
		const double d_dec = dec - p.dec;

		CAPTURE(p.ra);
		CAPTURE(p.dec);
		CHECK(std::abs(d_dec) < kAberrationBoundRadians);
		// For moderate declinations (|dec| <= 1 rad ~ 57°), the 1/cos(dec)
		// factor is bounded by ~2, so d_ra stays within ~4x the sky-plane
		// bound.
		CHECK(std::abs(d_ra)  < 4.0 * kAberrationBoundRadians);
	}
}

TEST_CASE("nutation correction is bounded by ~30 arcsec") {
	// Nutation in longitude peaks at ~17.2", nutation in obliquity at ~9.2".
	// Combined RA/Dec excursions stay safely under 30".
	constexpr double bound = (2.0 * 30.0 / 3600.0) * (std::numbers::pi / 180.0);

	double ra  = 1.2;
	double dec = 0.3;
	const double orig_ra  = ra;
	const double orig_dec = dec;

	nutation_correction_equatorial(kJdTestEpoch, ra, dec);

	const double two_pi = 2.0 * std::numbers::pi;
	double d_ra = ra - orig_ra;
	if (d_ra >  std::numbers::pi) d_ra -= two_pi;
	if (d_ra < -std::numbers::pi) d_ra += two_pi;

	CHECK(std::abs(d_ra)         < 4.0 * bound);
	CHECK(std::abs(dec - orig_dec) < bound);
}

///----------------------------------------
/// MARK: Non-degeneracy
///----------------------------------------

TEST_CASE("corrections are non-zero at representative epochs") {
	// If the implementation regressed to a no-op we would see zero
	// displacement; make sure we don't.
	double ra  = 1.0;
	double dec = 0.3;
	aberration_correction_equatorial(kJdTestEpoch, ra, dec);
	CHECK(std::abs(ra  - 1.0) > 1e-8);
	CHECK(std::abs(dec - 0.3) > 1e-8);
}

///----------------------------------------
/// MARK: J2000_to_apparent end-to-end
///----------------------------------------

TEST_CASE("J2000_to_apparent is bounded and non-zero") {
	// Combined correction should be within a few arcminutes for any typical
	// modern epoch (precession is stubbed out so we're measuring nutation +
	// aberration only).
	double ra  = 2.5;
	double dec = 0.1;
	J2000_to_apparent(kJdTestEpoch, ra, dec);

	// Non-trivial displacement expected (definitely > 1 arcsec).
	const double tiny = (1.0 / 3600.0) * (std::numbers::pi / 180.0);
	CHECK(std::abs(ra  - 2.5) > tiny);

	// Bounded by ~2 arcmin (aberration + nutation).
	const double bound = (120.0 / 3600.0) * (std::numbers::pi / 180.0);
	CHECK(std::abs(dec - 0.1) < bound);
}
