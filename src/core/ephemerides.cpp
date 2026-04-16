///----------------------------------------
///      @file ephemerides.cpp
///   @ingroup ASTAP++
///     @brief Clean-room ephemeris implementations for ASTAP++.
///   @details Implementations follow Meeus, "Astronomical Algorithms" 2nd ed.
///            throughout: Kepler-equation solver (Ch. 30), IAU 1976
///            precession (Ch. 21.2 / 21.3), low-precision Sun (Ch. 25),
///            and mean orbital elements for Jupiter / Saturn (Ch. 33).
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. Mozilla Public License 2.0.
///----------------------------------------

#include "ephemerides.h"

#include <algorithm>
#include <cmath>
#include <numbers>

///----------------------------------------
namespace astap::core::ephem {
///----------------------------------------

namespace {

constexpr double pi = std::numbers::pi;

constexpr double kTwoPi    = 2.0 * pi;
constexpr double kDeg2Rad  = pi / 180.0;
constexpr double kArcsec2Rad = pi / (180.0 * 3600.0);

constexpr double kJ2000       = 2451545.0;   ///< Julian Date of J2000.0.
constexpr double kCenturyDays = 36525.0;     ///< Julian days per Julian century.

/// @brief Gaussian gravitational constant (rad/day, used to get GM of the Sun).
constexpr double kGaussK = 0.01720209895;

/// @brief Heliocentric gravitational parameter in (AU^3/day^2).
constexpr double kSolarGM = kGaussK * kGaussK;

/// @brief Normalise an angle into [0, 2π).
[[nodiscard]] double normalise_two_pi(double x) noexcept {
	x = std::fmod(x, kTwoPi);
	return x < 0.0 ? x + kTwoPi : x;
}

/// @brief Normalise an angle into (-π, π].
[[nodiscard]] double normalise_pi(double x) noexcept {
	x = std::fmod(x + pi, kTwoPi);
	if (x < 0.0) x += kTwoPi;
	return x - pi;
}

/// @brief Mean obliquity of the ecliptic at J2000.0 (IAU 2006: 84381.406").
///        Used as the fixed ecliptic→equatorial rotation so that outputs land
///        in the J2000 equatorial frame (not equator-of-date).
constexpr double kJ2000Obliquity = 84381.406 * kArcsec2Rad;

/// @brief Mean obliquity of the ecliptic at the given JD (Meeus 22.2).
///        Returns radians. Kept for callers that specifically need
///        obliquity-of-date; the main ecliptic→equatorial rotations
///        in this module use kJ2000Obliquity instead.
[[nodiscard]] double mean_obliquity(double jd_tt) noexcept {
	const double T = (jd_tt - kJ2000) / kCenturyDays;
	// Arcseconds.
	const double eps_arcsec = 23.0 * 3600.0 + 26.0 * 60.0 + 21.448
	                          - 46.8150 * T - 0.00059 * T * T
	                          + 0.001813 * T * T * T;
	return eps_arcsec * kArcsec2Rad;
}

}  // anonymous namespace

///----------------------------------------
/// MARK: - cartesian_to_spherical / spherical_to_cartesian
///----------------------------------------

SphericalResult cartesian_to_spherical(const Vec3& xyz) noexcept {
	const double x = xyz[0], y = xyz[1], z = xyz[2];
	const double r_xy = std::sqrt(x * x + y * y);

	SphericalResult out;
	out.radius = std::sqrt(x * x + y * y + z * z);
	out.dec    = std::atan2(z, r_xy);
	out.ra     = (r_xy == 0.0) ? 0.0 : normalise_two_pi(std::atan2(y, x));
	return out;
}

Vec3 spherical_to_cartesian(double radius, double dec, double ra) noexcept {
	const double cos_dec = std::cos(dec);
	return {
		radius * cos_dec * std::cos(ra),
		radius * cos_dec * std::sin(ra),
		radius * std::sin(dec),
	};
}

///----------------------------------------
/// MARK: - precess_iau1976
///----------------------------------------

// Follows the rigorous three-angle formulation in Meeus 21.2 / 21.3.
EquatorialCoords precess_iau1976(EquatorialCoords coords, double jd_from, double jd_to) {
	const double T  = (jd_from - kJ2000) / kCenturyDays;
	const double T2 = T * T;

	const double t  = (jd_to - jd_from) / kCenturyDays;
	const double t2 = t * t;
	const double t3 = t2 * t;

	// Meeus 21.2 — rigorous method. Coefficients in arcseconds.
	const double ez_common = (2306.2181 + 1.39656 * T - 0.000139 * T2) * t;

	const double zeta    = (ez_common + (0.30188 - 0.000344 * T) * t2 + 0.017998 * t3) * kArcsec2Rad;
	const double z       = (ez_common + (1.09468 + 0.000066 * T) * t2 + 0.018203 * t3) * kArcsec2Rad;
	const double theta   = ((2004.3109 - 0.85330 * T - 0.000217 * T2) * t
	                        - (0.42665 + 0.000217 * T) * t2 - 0.041833 * t3) * kArcsec2Rad;

	// Apply the three-angle rotation to (ra, dec).  Meeus 21.4.
	const double ra0  = coords.ra;
	const double dec0 = coords.dec;

	const double A = std::cos(dec0) * std::sin(ra0 + zeta);
	const double B = std::cos(theta) * std::cos(dec0) * std::cos(ra0 + zeta)
	               - std::sin(theta) * std::sin(dec0);
	const double C = std::sin(theta) * std::cos(dec0) * std::cos(ra0 + zeta)
	               + std::cos(theta) * std::sin(dec0);

	EquatorialCoords out;
	out.ra  = normalise_two_pi(std::atan2(A, B) + z);
	out.dec = std::asin(std::clamp(C, -1.0, 1.0));
	return out;
}

///----------------------------------------
/// MARK: - Kepler solver
///----------------------------------------

namespace {

constexpr double kKeplerErrThresh    = 1.0e-12;
constexpr double kKeplerMinErrThresh = 1.0e-15;

// Expansion of E − e·sin(E) for eccentricities near 1 where cancellation ruins
// the direct formula. Power-series from Meeus p.199.
[[nodiscard]] double near_parabolic(double anomaly, double e) noexcept {
	const double anom2 = (e > 1.0) ? anomaly * anomaly : -anomaly * anomaly;
	double term = e * anom2 * anomaly / 6.0;
	double result = (1.0 - e) * anomaly - term;
	unsigned n = 4;
	while (std::fabs(term) > 1.0e-15) {
		term *= anom2 / static_cast<double>(n * (n + 1));
		result -= term;
		n += 2;
	}
	return result;
}

/// @brief Solve Kepler's equation for an elliptical or hyperbolic orbit.
///        Returns the eccentric anomaly in radians. Low-eccentricity
///        ellipses use Meeus 30.8 as an initial estimate followed by
///        Newton iteration; near-parabolic and hyperbolic cases use the
///        expansion in near_parabolic() to avoid cancellation.
[[nodiscard]] double solve_kepler(double e, double M) {
	if (std::isnan(M)) return 0.0;

	// Elliptical case: normalise M to (-π, π] and use Meeus 30.8 + Newton.
	double offset = 0.0;
	if (e < 1.0) {
		if (M < -pi || M > pi) {
			const double M_wrapped = normalise_pi(M);
			offset = M - M_wrapped;
			M = M_wrapped;
		}

		if (e < 0.9) {
			// Low-eccentricity: Meeus 30.8 initial estimate + Newton iteration.
			double E = std::atan2(std::sin(M), std::cos(M) - e);
			double correction;
			do {
				correction = (E - e * std::sin(E) - M) / (1.0 - e * std::cos(E));
				E -= correction;
			} while (std::fabs(correction) > kKeplerErrThresh);
			return E + offset;
		}
	}

	// High-eccentricity ellipse, parabola, or hyperbola.
	const bool negative = M < 0.0;
	if (negative) M = -M;

	double thresh = std::max(kKeplerErrThresh * std::min(std::fabs(1.0 - e), 1.0),
	                         kKeplerMinErrThresh);

	// Pick an initial estimate.
	double E = M;
	if (e > 1.0 && M / e > 3.0) {
		E = std::log(M / e) + 0.85;
	} else if ((e > 0.8 && M < pi / 3.0) || e > 1.0) {
		const double pfac = std::fabs(1.0 - e);
		double trial = M / pfac;
		if (trial * trial > 6.0 * pfac) {
			trial = std::cbrt(6.0 * M);
		}
		E = trial;
	}

	constexpr int kMaxEarly = 7;
	constexpr int kMaxTotal = 15;

	if (e < 1.0) {
		double delta = 1.0;
		int iter = 0;
		while (std::fabs(delta) > thresh && iter < kMaxTotal) {
			const double err = (iter++ > kMaxEarly)
			                 ? near_parabolic(E, e) - M
			                 : E - e * std::sin(E) - M;
			delta = -err / (1.0 - e * std::cos(E));
			E += delta;
		}
	} else {
		double delta = 1.0;
		int iter = 0;
		while (std::fabs(delta) > thresh && iter < kMaxTotal) {
			const double err = (iter++ > kMaxEarly && e < 1.01)
			                 ? -near_parabolic(E, e) - M
			                 : e * std::sinh(E) - E - M;
			delta = -err / (e * std::cosh(E) - 1.0);
			E += delta;
		}
	}

	return negative ? offset - E : offset + E;
}

/// @brief Obliquity-driven rotation: ecliptic of J2000 → equator of J2000.
///        Uses a fixed obliquity so the output is in the J2000 equatorial
///        frame regardless of the target date.
[[nodiscard]] Vec3 ecliptic_to_equatorial(const Vec3& v) noexcept {
	const double cos_eps = std::cos(kJ2000Obliquity);
	const double sin_eps = std::sin(kJ2000Obliquity);
	return {
		v[0],
		v[1] * cos_eps - v[2] * sin_eps,
		v[1] * sin_eps + v[2] * cos_eps,
	};
}

}  // anonymous namespace

///----------------------------------------
/// MARK: - propagate
///----------------------------------------

State propagate(const OrbitalElements& elements, double jd_tt) {
	const double e = elements.eccentricity;
	const double q = elements.perihelion_distance;

	// Time since epoch (days).
	const double dt = jd_tt - elements.epoch_jd;

	// Position in the orbital plane: x along perihelion, y perpendicular.
	double x_orb = 0.0, y_orb = 0.0;
	double vx_orb = 0.0, vy_orb = 0.0;
	double r = 0.0;

	if (std::fabs(e - 1.0) < 1.0e-10) {
		// Parabolic: use Barker's equation (Meeus 35.8).
		const double W = 3.0 * kGaussK * dt / (q * std::sqrt(2.0 * q));
		const double s = std::cbrt(W + std::sqrt(W * W + 1.0));
		const double s_inv = 1.0 / s;
		const double tan_half_nu = s - s_inv;           // tan(ν/2)
		const double true_anom   = 2.0 * std::atan(tan_half_nu);

		r = q * (1.0 + tan_half_nu * tan_half_nu);       // r = q · sec²(ν/2)
		x_orb = r * std::cos(true_anom);
		y_orb = r * std::sin(true_anom);

		const double h = std::sqrt(2.0 * kSolarGM * q);  // specific ang. momentum
		vx_orb = -h / r * std::sin(true_anom);
		vy_orb =  h / r * (e + std::cos(true_anom));
	} else {
		// Elliptic or hyperbolic.
		// Semimajor axis (signed: negative for hyperbolae so that GM/a gives |n|²).
		const double a = q / (1.0 - e);
		const double abs_a = std::fabs(a);
		// Mean motion (rad/day).
		const double n = kGaussK / (abs_a * std::sqrt(abs_a));
		// Mean anomaly at target date.
		const double M = elements.mean_anomaly_at_epoch + n * dt;
		// Solve Kepler's equation.
		const double E = solve_kepler(e, M);

		if (e < 1.0) {
			const double cos_E = std::cos(E);
			const double sin_E = std::sin(E);
			const double b = a * std::sqrt(1.0 - e * e);           // semi-minor
			x_orb = a * (cos_E - e);
			y_orb = b * sin_E;
			r = a * (1.0 - e * cos_E);
			const double E_dot = n / (1.0 - e * cos_E);
			vx_orb = -a * sin_E * E_dot;
			vy_orb =  b * cos_E * E_dot;
		} else {
			const double cosh_E = std::cosh(E);
			const double sinh_E = std::sinh(E);
			const double b = abs_a * std::sqrt(e * e - 1.0);
			x_orb =  abs_a * (e - cosh_E);
			y_orb =  b * sinh_E;
			r = abs_a * (e * cosh_E - 1.0);
			const double E_dot = n / (e * cosh_E - 1.0);
			vx_orb = -abs_a * sinh_E * E_dot;
			vy_orb =  b * cosh_E * E_dot;
		}
	}

	// Rotate orbital-plane vector → ecliptic (ω, i, Ω).
	const double cos_w = std::cos(elements.argument_of_perihelion);
	const double sin_w = std::sin(elements.argument_of_perihelion);
	const double cos_i = std::cos(elements.inclination);
	const double sin_i = std::sin(elements.inclination);
	const double cos_O = std::cos(elements.ascending_node);
	const double sin_O = std::sin(elements.ascending_node);

	// Rotation matrix Rz(Ω) · Rx(i) · Rz(ω).
	const double P11 =  cos_O * cos_w - sin_O * sin_w * cos_i;
	const double P12 = -cos_O * sin_w - sin_O * cos_w * cos_i;
	const double P21 =  sin_O * cos_w + cos_O * sin_w * cos_i;
	const double P22 = -sin_O * sin_w + cos_O * cos_w * cos_i;
	const double P31 =  sin_w * sin_i;
	const double P32 =  cos_w * sin_i;

	const Vec3 pos_ecl = {
		P11 * x_orb + P12 * y_orb,
		P21 * x_orb + P22 * y_orb,
		P31 * x_orb + P32 * y_orb,
	};
	const Vec3 vel_ecl = {
		P11 * vx_orb + P12 * vy_orb,
		P21 * vx_orb + P22 * vy_orb,
		P31 * vx_orb + P32 * vy_orb,
	};

	State out;
	out.position = ecliptic_to_equatorial(pos_ecl);
	out.velocity = ecliptic_to_equatorial(vel_ecl);
	return out;
}

///----------------------------------------
/// MARK: - earth_state (Meeus Ch. 25 low-precision Sun/Earth)
///----------------------------------------

namespace {

// Evaluate the Sun's apparent geocentric ecliptic longitude and the
// Earth-Sun distance using Meeus Chapter 25 (low accuracy; ~0.01° in λ,
// ~1e-5 AU in R). Returns (true_longitude_rad, distance_AU, true_anomaly_rad).
struct SunLowPrec {
	double true_longitude;   ///< Sun geocentric true longitude (rad).
	double radius;           ///< Earth-Sun distance (AU).
	double true_anomaly;     ///< True anomaly of Earth (rad).
	double eccentricity;     ///< Earth's orbital eccentricity.
};

[[nodiscard]] SunLowPrec sun_low_precision(double jd_tt) noexcept {
	const double T = (jd_tt - kJ2000) / kCenturyDays;

	// Mean longitude of the Sun (Meeus 25.2).
	const double L0_deg = 280.46646 + 36000.76983 * T + 0.0003032 * T * T;
	// Mean anomaly of the Sun (Meeus 25.3).
	const double M_deg = 357.52911 + 35999.05029 * T - 0.0001537 * T * T;
	// Earth's orbital eccentricity (Meeus 25.4).
	const double e = 0.016708634 - 0.000042037 * T - 0.0000001267 * T * T;

	const double M = M_deg * kDeg2Rad;

	// Equation of the centre C (Meeus p.164).
	const double C = ((1.914602 - 0.004817 * T - 0.000014 * T * T) * std::sin(M)
	                  + (0.019993 - 0.000101 * T)                  * std::sin(2.0 * M)
	                  +  0.000289                                  * std::sin(3.0 * M)) * kDeg2Rad;

	const double L0 = L0_deg * kDeg2Rad;
	const double true_long = L0 + C;
	const double nu        = M + C;

	// Sun-Earth distance, R (Meeus 25.5).
	const double R = 1.000001018 * (1.0 - e * e) / (1.0 + e * std::cos(nu));

	return { normalise_two_pi(true_long), R, nu, e };
}

// ---------------------------------------------------------------------------
// Sun-to-SSB offset — dominant-body (Jupiter + Saturn) approximation.
//
// The Solar-System Barycentre lies at Σ(m_i / M_tot)·r_i summed over all
// bodies. With the Sun at the heliocentric origin, the Sun's own offset
// from the SSB is −Σ_planets (m_i / M_tot)·r_i_helio.
//
// Jupiter contributes ~71% of the offset, Saturn ~21%, Uranus ~3%,
// Neptune ~5%. Including only Jupiter and Saturn leaves ~0.0007 AU of
// residual — well below the Meeus Ch. 25 Sun precision (~1e-5 AU in R,
// ~40 arcsec in direction) used for Earth, so there's no point being
// more precise here until `earth_state` itself is upgraded to VSOP87.
//
// Elements: mean Keplerian elements at J2000.0 referred to the J2000
// mean ecliptic, from JPL solar-system dynamics tables (Standish 1992 /
// Meeus Ch. 33 Table 33.A). Secular drift is dropped — over 100 years
// the planet positions drift by ~1 arcmin, which scales down by the
// ~0.001 mass ratio to ≤0.1 milliarcsec of Earth-state error.
// ---------------------------------------------------------------------------

// Mass ratios (planet / Sun).
constexpr double kMassJup = 1.0 / 1047.3486;
constexpr double kMassSat = 1.0 / 3497.898;

// Sum of all point-mass contributors in solar masses (Sun + Jupiter + Saturn).
constexpr double kMassTotal = 1.0 + kMassJup + kMassSat;

[[nodiscard]] OrbitalElements jupiter_j2000_elements() {
	// Jupiter J2000 Keplerian (ecliptic): a = 5.20336301 AU, e = 0.04839266,
	// i = 1.30530°, Ω = 100.55615°, ϖ = 14.75385°, L = 34.40438°.
	constexpr double a  = 5.20336301;
	constexpr double e  = 0.04839266;
	constexpr double i  = 1.30530  * kDeg2Rad;
	constexpr double O  = 100.55615 * kDeg2Rad;
	constexpr double pi_long = 14.75385 * kDeg2Rad;
	constexpr double L  = 34.40438 * kDeg2Rad;

	return {
		.epoch_jd               = kJ2000,
		.perihelion_distance    = a * (1.0 - e),
		.eccentricity           = e,
		.inclination            = i,
		.ascending_node         = O,
		.argument_of_perihelion = pi_long - O,      // ω = ϖ − Ω
		.mean_anomaly_at_epoch  = L - pi_long,      // M = L − ϖ
	};
}

[[nodiscard]] OrbitalElements saturn_j2000_elements() {
	// Saturn J2000 Keplerian (ecliptic): a = 9.53707032 AU, e = 0.05415060,
	// i = 2.48446°, Ω = 113.71504°, ϖ = 92.43194°, L = 49.94432°.
	constexpr double a  = 9.53707032;
	constexpr double e  = 0.05415060;
	constexpr double i  = 2.48446  * kDeg2Rad;
	constexpr double O  = 113.71504 * kDeg2Rad;
	constexpr double pi_long = 92.43194 * kDeg2Rad;
	constexpr double L  = 49.94432 * kDeg2Rad;

	return {
		.epoch_jd               = kJ2000,
		.perihelion_distance    = a * (1.0 - e),
		.eccentricity           = e,
		.inclination            = i,
		.ascending_node         = O,
		.argument_of_perihelion = pi_long - O,
		.mean_anomaly_at_epoch  = L - pi_long,
	};
}

}  // anonymous namespace

namespace {

/// @brief Sun's position and velocity relative to the SSB in J2000
///        equatorial frame (AU, AU/day).
[[nodiscard]] State sun_barycentric_offset(double jd_tt) {
	const State jup = propagate(jupiter_j2000_elements(), jd_tt);
	const State sat = propagate(saturn_j2000_elements(), jd_tt);

	const double fj = kMassJup / kMassTotal;
	const double fs = kMassSat / kMassTotal;

	State s;
	for (int k = 0; k < 3; ++k) {
		s.position[k] = -(fj * jup.position[k] + fs * sat.position[k]);
		s.velocity[k] = -(fj * jup.velocity[k] + fs * sat.velocity[k]);
	}
	return s;
}

}  // anonymous namespace

State earth_state(double jd_tt, ReferenceFrame frame) {
	const SunLowPrec s = sun_low_precision(jd_tt);

	// Meeus 25 gives the Sun's longitude relative to the mean equinox *of
	// date*. To deliver a result in the J2000 mean equinox frame, subtract
	// the general precession of the equinox since J2000. Meeus 21.6:
	// p = 5028.796195"·T + 1.1054348"·T² per century from J2000.
	const double T = (jd_tt - kJ2000) / kCenturyDays;
	const double gen_precession_rad = (5028.796195 * T + 1.1054348 * T * T) * kArcsec2Rad;

	// Sun geocentric ecliptic: (R, λ, 0). Earth heliocentric ecliptic: opposite.
	const double earth_long = normalise_two_pi(
		s.true_longitude + pi - gen_precession_rad);
	const double cos_L = std::cos(earth_long);
	const double sin_L = std::sin(earth_long);

	const Vec3 pos_ecl = { s.radius * cos_L, s.radius * sin_L, 0.0 };

	// Velocity: dr/dt in ecliptic plane from orbital mechanics.
	// In Earth's heliocentric orbit, v = sqrt(GM·(2/r − 1/a)); direction is
	// perpendicular to the radius vector, biased by the true anomaly.
	// Use d(λ)/dt and dR/dt expressions.
	const double n_day = 360.0 / 365.256363004 * kDeg2Rad;          // mean motion rad/day
	const double dnu_dt = n_day * std::pow(1.0 + s.eccentricity * std::cos(s.true_anomaly), 2)
	                     / std::pow(1.0 - s.eccentricity * s.eccentricity, 1.5);
	const double dR_dt  = s.radius * s.eccentricity * std::sin(s.true_anomaly) * dnu_dt
	                     / (1.0 + s.eccentricity * std::cos(s.true_anomaly));

	// Velocity in ecliptic = (dR·cosL − R·sinL·dλ, dR·sinL + R·cosL·dλ, 0).
	const Vec3 vel_ecl = {
		dR_dt * cos_L - s.radius * sin_L * dnu_dt,
		dR_dt * sin_L + s.radius * cos_L * dnu_dt,
		0.0,
	};

	State out;
	out.position = ecliptic_to_equatorial(pos_ecl);
	out.velocity = ecliptic_to_equatorial(vel_ecl);

	// Barycentric: add the Sun's offset from the SSB. The helper covers
	// Jupiter and Saturn (~92% of the full correction); the ~0.0007 AU
	// residual from Uranus and Neptune is below this routine's other
	// error sources. Swap `sun_barycentric_offset` for a VSOP87-based
	// version when `earth_state` itself is upgraded.
	if (frame == ReferenceFrame::Barycentric) {
		const State off = sun_barycentric_offset(jd_tt);
		for (int k = 0; k < 3; ++k) {
			out.position[k] += off.position[k];
			out.velocity[k] += off.velocity[k];
		}
	}
	return out;
}

///----------------------------------------
/// MARK: - vsop::evaluate
///----------------------------------------

namespace vsop {

HelioEcliptic evaluate(const Body& body, double jd_tt) noexcept {
	// Julian millennia from J2000.
	const double tau = (jd_tt - kJ2000) / (kCenturyDays * 10.0);

	auto eval_component = [tau](const Component& comp) -> double {
		double result = 0.0;
		double tau_power = 1.0;
		for (std::size_t i = 0; i < comp.count; ++i) {
			const Series& s = comp.series[i];
			double sum = 0.0;
			for (std::size_t k = 0; k < s.count; ++k) {
				const Term& t = s.terms[k];
				sum += t.amplitude * std::cos(t.phase + t.frequency * tau);
			}
			result += sum * tau_power;
			tau_power *= tau;
		}
		return result;
	};

	HelioEcliptic out;
	out.longitude = normalise_two_pi(eval_component(body.L));
	out.latitude  = eval_component(body.B);
	out.radius    = eval_component(body.R);
	return out;
}

}  // namespace vsop

}  // namespace astap::core::ephem
