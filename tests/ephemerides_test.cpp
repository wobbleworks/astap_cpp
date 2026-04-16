///----------------------------------------
///     @file ephemerides_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for the ephemerides module (src/core/ephemerides.cpp).
///  @details Covers the public API in astap::core::ephem:
///           cartesian_to_spherical / spherical_to_cartesian, precess_iau1976,
///           earth_state (Heliocentric and Barycentric), propagate for
///           elliptic/parabolic/hyperbolic orbits, and vsop::evaluate. Tests
///           rely on physical invariants (Earth-Sun distance bounds, Keplerian
///           period closure, frame invertibility) rather than reference tables.
///   @author Created by John Stephen on 4/16/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "core/ephemerides.h"

#include <cmath>
#include <numbers>

using namespace astap::core::ephem;
using doctest::Approx;

static constexpr double kPi    = std::numbers::pi;
static constexpr double kTwoPi = 2.0 * kPi;
static constexpr double kArcsec2Rad = kPi / (180.0 * 3600.0);

static constexpr double kJ2000 = 2451545.0;

///----------------------------------------
/// MARK: cartesian_to_spherical / spherical_to_cartesian
///----------------------------------------

TEST_CASE("cartesian_to_spherical: zero vector") {
	const auto s = cartesian_to_spherical({0.0, 0.0, 0.0});
	CHECK(s.radius == Approx(0.0));
	CHECK(s.ra     == Approx(0.0));
	CHECK(s.dec    == Approx(0.0));
}

TEST_CASE("cartesian_to_spherical: unit +X axis") {
	const auto s = cartesian_to_spherical({1.0, 0.0, 0.0});
	CHECK(s.radius == Approx(1.0));
	CHECK(s.ra     == Approx(0.0));
	CHECK(s.dec    == Approx(0.0));
}

TEST_CASE("cartesian_to_spherical: unit +Y axis (ra = π/2)") {
	const auto s = cartesian_to_spherical({0.0, 1.0, 0.0});
	CHECK(s.radius == Approx(1.0));
	CHECK(s.ra     == Approx(kPi / 2.0));
	CHECK(s.dec    == Approx(0.0));
}

TEST_CASE("cartesian_to_spherical: unit -X axis (ra = π)") {
	const auto s = cartesian_to_spherical({-1.0, 0.0, 0.0});
	CHECK(s.radius == Approx(1.0));
	CHECK(s.ra     == Approx(kPi));
	CHECK(s.dec    == Approx(0.0));
}

TEST_CASE("cartesian_to_spherical: ra wraps to [0, 2π)") {
	// atan2 returns (-π, π]; -Y direction should produce ra = 3π/2, not -π/2.
	const auto s = cartesian_to_spherical({0.0, -1.0, 0.0});
	CHECK(s.ra >= 0.0);
	CHECK(s.ra <  kTwoPi);
	CHECK(s.ra == Approx(3.0 * kPi / 2.0));
}

TEST_CASE("cartesian_to_spherical: +Z gives declination +π/2") {
	const auto s = cartesian_to_spherical({0.0, 0.0, 1.0});
	CHECK(s.radius == Approx(1.0));
	CHECK(s.dec    == Approx(kPi / 2.0));
}

TEST_CASE("cartesian_to_spherical: -Z gives declination -π/2") {
	const auto s = cartesian_to_spherical({0.0, 0.0, -1.0});
	CHECK(s.radius == Approx(1.0));
	CHECK(s.dec    == Approx(-kPi / 2.0));
}

TEST_CASE("spherical_to_cartesian: known unit directions") {
	SUBCASE("ra=0, dec=0 -> +X") {
		const auto v = spherical_to_cartesian(1.0, 0.0, 0.0);
		CHECK(v[0] == Approx(1.0));
		CHECK(v[1] == Approx(0.0));
		CHECK(v[2] == Approx(0.0));
	}
	SUBCASE("ra=π/2, dec=0 -> +Y") {
		const auto v = spherical_to_cartesian(1.0, 0.0, kPi / 2.0);
		CHECK(v[0] == Approx(0.0).epsilon(1e-12));
		CHECK(v[1] == Approx(1.0));
		CHECK(v[2] == Approx(0.0));
	}
	SUBCASE("dec=π/2 -> +Z") {
		const auto v = spherical_to_cartesian(2.0, kPi / 2.0, 0.0);
		CHECK(v[0] == Approx(0.0).epsilon(1e-12));
		CHECK(v[1] == Approx(0.0));
		CHECK(v[2] == Approx(2.0));
	}
}

TEST_CASE("coordinate conversion round-trips to identity") {
	// Sample ra/dec values covering the full sphere.
	struct Point { double r, ra, dec; };
	const Point samples[] = {
		{1.0,       0.3,        0.5},
		{2.5,       kPi,       -0.4},
		{0.5,       1.8,        kPi / 4},
		{10.0,      4.7,       -kPi / 3},
		{1.0,       0.0,        0.0},
	};
	for (const auto& p : samples) {
		const auto v = spherical_to_cartesian(p.r, p.dec, p.ra);
		const auto s = cartesian_to_spherical(v);
		CAPTURE(p.r); CAPTURE(p.ra); CAPTURE(p.dec);
		CHECK(s.radius == Approx(p.r).epsilon(1e-12));
		CHECK(s.ra     == Approx(p.ra).epsilon(1e-12));
		CHECK(s.dec    == Approx(p.dec).epsilon(1e-12));
	}
}

///----------------------------------------
/// MARK: precess_iau1976
///----------------------------------------

TEST_CASE("precess_iau1976: zero interval is identity") {
	const EquatorialCoords c{1.2, 0.3};
	const auto out = precess_iau1976(c, kJ2000, kJ2000);
	CHECK(out.ra  == Approx(c.ra).epsilon(1e-12));
	CHECK(out.dec == Approx(c.dec).epsilon(1e-12));
}

TEST_CASE("precess_iau1976: invertible (round-trip < 10 µas)") {
	const double jd_then = kJ2000 + 50.0 * 365.25;   // J2050
	const EquatorialCoords samples[] = {
		{0.5, 0.2},
		{3.1, -0.4},
		{6.0, 0.7},
		{0.0, 0.0},
	};
	constexpr double kTenMicroarcsecRad = 10e-6 * kArcsec2Rad;
	for (const auto& c : samples) {
		const auto fwd = precess_iau1976(c, kJ2000, jd_then);
		const auto back = precess_iau1976(fwd, jd_then, kJ2000);

		// Account for RA wrap: take the signed shortest-arc distance.
		double d_ra = back.ra - c.ra;
		while (d_ra >  kPi) d_ra -= kTwoPi;
		while (d_ra < -kPi) d_ra += kTwoPi;

		CAPTURE(c.ra); CAPTURE(c.dec);
		CHECK(std::abs(d_ra)             < kTenMicroarcsecRad);
		CHECK(std::abs(back.dec - c.dec) < kTenMicroarcsecRad);
	}
}

TEST_CASE("precess_iau1976: RA is normalised to [0, 2π)") {
	const auto out = precess_iau1976({0.001, 0.1}, kJ2000, kJ2000 + 365.25);
	CHECK(out.ra >= 0.0);
	CHECK(out.ra <  kTwoPi);
}

TEST_CASE("precess_iau1976: J2000->J2025 non-trivial and bounded") {
	// Over 25 years the annual precession of ~50"/yr accumulates to
	// ~21 arcmin in ra. The per-pair shift at the equator is ~21 arcmin;
	// at high declinations it's smaller in absolute terms but can be
	// amplified in RA by 1/cos(dec). Here we just check "non-zero" and
	// "well under 1 degree" for a low-dec equator point.
	const double jd_2025 = kJ2000 + 25.0 * 365.25;
	const auto out = precess_iau1976({2.0, 0.1}, kJ2000, jd_2025);

	double d_ra = out.ra - 2.0;
	while (d_ra >  kPi) d_ra -= kTwoPi;
	while (d_ra < -kPi) d_ra += kTwoPi;

	const double arcsec = 180.0 * 3600.0 / kPi;
	const double d_ra_arcsec = std::abs(d_ra) * arcsec;

	// Non-trivial: at least 10" shift.
	CHECK(d_ra_arcsec > 10.0);
	// Bounded: well under 1 degree.
	CHECK(d_ra_arcsec < 3600.0);
}

///----------------------------------------
/// MARK: earth_state
///----------------------------------------

TEST_CASE("earth_state: Sun-Earth distance is within [0.983, 1.017] AU") {
	// Sample a few epochs spanning a year; Earth's perihelion/aphelion are
	// 0.9833 / 1.0167 AU so any year-round sample should fall inside this
	// range.
	for (int k = 0; k < 12; ++k) {
		const double jd = kJ2000 + k * (365.25 / 12.0);
		const auto s = earth_state(jd, ReferenceFrame::Heliocentric);
		const double r = std::sqrt(s.position[0] * s.position[0]
		                          + s.position[1] * s.position[1]
		                          + s.position[2] * s.position[2]);
		CAPTURE(jd); CAPTURE(r);
		CHECK(r > 0.982);
		CHECK(r < 1.018);
	}
}

TEST_CASE("earth_state: speed ≈ Earth mean orbital velocity") {
	// Mean orbital speed = 2π AU/sidereal year = 2π / 365.2564 AU/day
	// ≈ 0.01721 AU/day. Instantaneous speed varies ~3% around this (ev).
	constexpr double mean_speed = 2.0 * std::numbers::pi / 365.256363004;
	for (int k = 0; k < 12; ++k) {
		const double jd = kJ2000 + k * (365.25 / 12.0);
		const auto s = earth_state(jd, ReferenceFrame::Heliocentric);
		const double v = std::sqrt(s.velocity[0] * s.velocity[0]
		                          + s.velocity[1] * s.velocity[1]
		                          + s.velocity[2] * s.velocity[2]);
		CAPTURE(jd); CAPTURE(v);
		CHECK(v > mean_speed * 0.95);
		CHECK(v < mean_speed * 1.05);
	}
}

TEST_CASE("earth_state: Z component non-trivial (ecliptic tilt)") {
	// Earth's position in equatorial coordinates has a Z component that
	// varies sinusoidally through the year (~sin(ε) × R ≈ 0.4 AU peak).
	double max_abs_z = 0.0;
	for (int k = 0; k < 36; ++k) {
		const double jd = kJ2000 + k * (365.25 / 36.0);
		const auto s = earth_state(jd, ReferenceFrame::Heliocentric);
		max_abs_z = std::max(max_abs_z, std::abs(s.position[2]));
	}
	// Peak should be roughly R·sin(23.44°) ≈ 0.40 AU. Demand at least 0.30.
	CHECK(max_abs_z > 0.30);
	CHECK(max_abs_z < 0.50);
}

TEST_CASE("earth_state: Barycentric offset within [0.001, 0.010] AU") {
	// Jupiter + Saturn combined mass-weighted reflex of the Sun peaks at
	// ~0.008 AU depending on their relative longitudes.
	for (int k = 0; k < 4; ++k) {
		const double jd = kJ2000 + k * 3000.0;   // sample across ~33 yr
		const auto helio = earth_state(jd, ReferenceFrame::Heliocentric);
		const auto bary  = earth_state(jd, ReferenceFrame::Barycentric);
		const double dx = bary.position[0] - helio.position[0];
		const double dy = bary.position[1] - helio.position[1];
		const double dz = bary.position[2] - helio.position[2];
		const double dmag = std::sqrt(dx * dx + dy * dy + dz * dz);
		CAPTURE(jd); CAPTURE(dmag);
		CHECK(dmag > 0.001);
		CHECK(dmag < 0.010);
	}
}

TEST_CASE("earth_state: Barycentric velocity offset is tiny but non-zero") {
	// Sun's reflex velocity around the SSB is dominated by Jupiter's
	// orbital motion: v_reflex ≈ (m_Jup / M_Sun) × v_Jup ≈ 10⁻⁵ AU/day.
	const auto helio = earth_state(kJ2000, ReferenceFrame::Heliocentric);
	const auto bary  = earth_state(kJ2000, ReferenceFrame::Barycentric);
	const double dvx = bary.velocity[0] - helio.velocity[0];
	const double dvy = bary.velocity[1] - helio.velocity[1];
	const double dvz = bary.velocity[2] - helio.velocity[2];
	const double dvmag = std::sqrt(dvx * dvx + dvy * dvy + dvz * dvz);
	CHECK(dvmag > 1e-7);   // non-zero
	CHECK(dvmag < 1e-4);   // but tiny (< 1% of Earth's own speed)
}

///----------------------------------------
/// MARK: propagate — elliptic orbits
///----------------------------------------

namespace {

// Ceres-like asteroid: a ≈ 2.77 AU, e ≈ 0.076, reasonable J2000 elements.
OrbitalElements ceres_like() {
	OrbitalElements e{};
	e.epoch_jd               = kJ2000;
	e.perihelion_distance    = 2.77 * (1.0 - 0.076);    // q = a(1-e)
	e.eccentricity           = 0.076;
	e.inclination            = 10.6 * kPi / 180.0;
	e.ascending_node         = 80.3 * kPi / 180.0;
	e.argument_of_perihelion = 73.6 * kPi / 180.0;
	e.mean_anomaly_at_epoch  = 77.3 * kPi / 180.0;
	return e;
}

}  // namespace

TEST_CASE("propagate: Keplerian period closure for elliptic orbit") {
	// After one full orbital period, the body should return to its
	// starting position to within floating-point tolerance.
	const auto elements = ceres_like();
	const double a = elements.perihelion_distance / (1.0 - elements.eccentricity);
	constexpr double kGaussK = 0.01720209895;
	const double period_days = kTwoPi / kGaussK * std::sqrt(a * a * a);

	const auto start = propagate(elements, elements.epoch_jd);
	const auto end   = propagate(elements, elements.epoch_jd + period_days);

	const double dx = end.position[0] - start.position[0];
	const double dy = end.position[1] - start.position[1];
	const double dz = end.position[2] - start.position[2];
	const double r_drift = std::sqrt(dx * dx + dy * dy + dz * dz);

	// Tolerance: ~1 part in 1e10 of the semimajor axis.
	CHECK(r_drift < 1e-8);
}

TEST_CASE("propagate: specific energy conserved for elliptic orbit") {
	// For a bound orbit around the Sun with GM = k², the specific
	// orbital energy ε = v²/2 − GM/r must equal −GM/(2a), invariant
	// under time propagation.
	const auto elements = ceres_like();
	const double a = elements.perihelion_distance / (1.0 - elements.eccentricity);
	constexpr double kGM = 0.01720209895 * 0.01720209895;
	const double expected_energy = -kGM / (2.0 * a);

	for (double dt : {0.0, 37.0, 333.0, 1500.0}) {
		const auto s = propagate(elements, elements.epoch_jd + dt);
		const double r  = std::sqrt(s.position[0] * s.position[0]
		                          + s.position[1] * s.position[1]
		                          + s.position[2] * s.position[2]);
		const double v2 = s.velocity[0] * s.velocity[0]
		                + s.velocity[1] * s.velocity[1]
		                + s.velocity[2] * s.velocity[2];
		const double energy = 0.5 * v2 - kGM / r;
		CAPTURE(dt);
		CHECK(energy == Approx(expected_energy).epsilon(1e-10));
	}
}

TEST_CASE("propagate: radius oscillates between q and apoapsis for elliptic") {
	// r must always satisfy q ≤ r ≤ a(1+e) = 2a − q for an elliptic orbit.
	const auto elements = ceres_like();
	const double q = elements.perihelion_distance;
	const double a = q / (1.0 - elements.eccentricity);
	const double apoapsis = a * (1.0 + elements.eccentricity);

	double r_min = 1e99, r_max = 0.0;
	for (int k = 0; k < 100; ++k) {
		const double jd = elements.epoch_jd + k * 40.0;  // sweep ~11 yr
		const auto s = propagate(elements, jd);
		const double r = std::sqrt(s.position[0] * s.position[0]
		                         + s.position[1] * s.position[1]
		                         + s.position[2] * s.position[2]);
		r_min = std::min(r_min, r);
		r_max = std::max(r_max, r);
	}
	// Every sample must respect the bounds.
	CHECK(r_min > q    - 1e-6);
	CHECK(r_max < apoapsis + 1e-6);
	// Sampling across a full period should actually visit near-perihelion
	// and near-apoapsis, so the observed range should span most of the
	// full oscillation.
	CHECK(r_max - r_min > 0.5 * (apoapsis - q));
}

///----------------------------------------
/// MARK: propagate — parabolic and hyperbolic orbits
///----------------------------------------

TEST_CASE("propagate: parabolic orbit returns finite state") {
	OrbitalElements parabolic{};
	parabolic.epoch_jd               = kJ2000;   // perihelion passage
	parabolic.perihelion_distance    = 1.5;      // AU
	parabolic.eccentricity           = 1.0;
	parabolic.inclination            = 0.3;
	parabolic.ascending_node         = 1.1;
	parabolic.argument_of_perihelion = 0.7;
	parabolic.mean_anomaly_at_epoch  = 0.0;      // at perihelion

	const auto at_peri = propagate(parabolic, kJ2000);
	const double r0 = std::sqrt(at_peri.position[0] * at_peri.position[0]
	                          + at_peri.position[1] * at_peri.position[1]
	                          + at_peri.position[2] * at_peri.position[2]);
	CHECK(r0 == Approx(parabolic.perihelion_distance).epsilon(1e-9));

	// After some days r must have increased (parabola leaves perihelion).
	const auto later = propagate(parabolic, kJ2000 + 200.0);
	const double r1 = std::sqrt(later.position[0] * later.position[0]
	                          + later.position[1] * later.position[1]
	                          + later.position[2] * later.position[2]);
	CHECK(r1 > r0);
	CHECK(std::isfinite(r1));
}

TEST_CASE("propagate: hyperbolic orbit, r strictly increasing past perihelion") {
	OrbitalElements hyperbolic{};
	hyperbolic.epoch_jd               = kJ2000;   // at perihelion
	hyperbolic.perihelion_distance    = 2.0;
	hyperbolic.eccentricity           = 1.4;
	hyperbolic.inclination            = 0.2;
	hyperbolic.ascending_node         = 0.0;
	hyperbolic.argument_of_perihelion = 0.0;
	hyperbolic.mean_anomaly_at_epoch  = 0.0;

	double r_prev = hyperbolic.perihelion_distance;
	for (int k = 1; k <= 5; ++k) {
		const double jd = kJ2000 + k * 500.0;
		const auto s = propagate(hyperbolic, jd);
		const double r = std::sqrt(s.position[0] * s.position[0]
		                         + s.position[1] * s.position[1]
		                         + s.position[2] * s.position[2]);
		CAPTURE(k); CAPTURE(r);
		CHECK(r > r_prev);
		r_prev = r;
	}
}

///----------------------------------------
/// MARK: vsop::evaluate
///----------------------------------------

TEST_CASE("vsop::evaluate: empty body returns zero") {
	vsop::Body empty{};   // all counts = 0
	const auto r = vsop::evaluate(empty, kJ2000);
	CHECK(r.longitude == Approx(0.0));
	CHECK(r.latitude  == Approx(0.0));
	CHECK(r.radius    == Approx(0.0));
}

TEST_CASE("vsop::evaluate: single constant term returns amplitude") {
	// One term with zero frequency and zero phase: result = A·cos(0) = A.
	vsop::Term t{ .amplitude = 1.5, .phase = 0.0, .frequency = 0.0 };
	vsop::Series series{ .terms = &t, .count = 1 };
	vsop::Component comp{ .series = &series, .count = 1 };

	vsop::Body body{};
	body.R = comp;   // populate radius component

	const auto r = vsop::evaluate(body, kJ2000);
	CHECK(r.radius == Approx(1.5));
	// L and B were not populated, so they stay zero.
	CHECK(r.longitude == Approx(0.0));
	CHECK(r.latitude  == Approx(0.0));
}

TEST_CASE("vsop::evaluate: τ⁰ and τ¹ time-powers accumulate correctly") {
	// τ=0 at J2000 → only τ⁰ contributes. At J2000 + 3652.5 d, τ = 0.01.
	vsop::Term t0{ .amplitude = 1.0, .phase = 0.0, .frequency = 0.0 };
	vsop::Term t1{ .amplitude = 5.0, .phase = 0.0, .frequency = 0.0 };

	vsop::Series s0{ .terms = &t0, .count = 1 };
	vsop::Series s1{ .terms = &t1, .count = 1 };
	std::array<vsop::Series, 2> series = { s0, s1 };

	vsop::Component comp{ .series = series.data(), .count = 2 };
	vsop::Body body{};
	body.R = comp;

	// At J2000 (τ=0): result = 1.0 + 5.0·0 = 1.0.
	const auto at_j2000 = vsop::evaluate(body, kJ2000);
	CHECK(at_j2000.radius == Approx(1.0));

	// τ = 0.01 → result = 1.0 + 5.0·0.01 = 1.05.
	const double jd_later = kJ2000 + 3652.5;   // 3652.5 d ≈ 0.01 millennia
	const auto at_later = vsop::evaluate(body, jd_later);
	CHECK(at_later.radius == Approx(1.05));
}
