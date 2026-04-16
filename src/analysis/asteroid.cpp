///----------------------------------------
///      @file asteroid.cpp
///   @ingroup analysis
///     @brief Asteroid and comet ephemeris computation and catalog parsing.
///   @details Faithful port of the pure-algorithmic functions from
///            unit_asteroid.pas. All Pascal 1-based string indexing has been
///            converted to 0-based std::string_view::substr calls.
///            Original copyright (C) 2021 by Han Kleijn, www.hnsky.org.
///            Licensed under the Mozilla Public License, v. 2.0.
///    @author Created by John Stephen on 4/15/26.
/// @copyright Copyright 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "asteroid.h"

#include <cmath>
#include <numbers>
#include <string>

#include "../core/ephemerides.h"
#include "../core/hjd.h"

// julian_calc lives in the stacking module; forward-declare to avoid pulling
// the full stacking header into the analysis module.
namespace astap::stacking {
[[nodiscard]] double julian_calc(int yyyy, int mm, double dd,
                                 double hours, double minutes, double seconds);
} // namespace astap::stacking

///----------------------------------------
namespace astap::analysis {
///----------------------------------------

namespace {

/// @brief Light-time for 1 AU in seconds.
static constexpr double kTau = 499.004782;

/// @brief Light-time for 1 AU in days (= kTau / 86400).
///        Used with ephem::propagate's AU/day velocities for the
///        planetary-aberration correction.
static constexpr double kTauDays = kTau / 86400.0;

/// @brief Parse a substring as a double, returning 0.0 on failure.
[[nodiscard]] double safe_stod(std::string_view sv) noexcept {
	// stod requires a null-terminated string; build a small local copy.
	try {
		auto s = std::string(sv);
		return std::stod(s);
	} catch (...) {
		return 0.0;
	}
}

/// @brief Parse a substring as an integer, returning 0 on failure.
[[nodiscard]] int safe_stoi(std::string_view sv) noexcept {
	try {
		auto s = std::string(sv);
		return std::stoi(s);
	} catch (...) {
		return 0;
	}
}

} // anonymous namespace

///----------------------------------------
/// MARK: deltaT_calc
///----------------------------------------

[[nodiscard]] double deltaT_calc(double jd) noexcept {
	auto y = 2000.0 + (jd - 2451544.5) / 365.25;
	auto year = static_cast<int>(std::round(y));
	auto result = 0.0;

	if (year >= 2016 && year <= 2020) {
		// (71-68.3)/5 = 0.54
		auto t = y - 2016.0;
		result = 68.3 + t * 0.54; // seconds
	} else if (year >= 2021 && year <= 2024) {
		// (73-71)/4 = 0.5
		auto t = y - 2021.0;
		result = 71.0 + t * 0.5; // seconds
	} else if (year >= 2025 && year <= 2049) {
		auto t = y - 2000.0;
		result = 61.46 + t * (0.32217 + t * 0.005589); // seconds
	} else if (year >= 2050 && year <= 2149) {
		auto u = (y - 1820.0) / 100.0;
		auto t = 2150.0 - y;
		result = -20.0 + 32.0 * u * u - 0.5788 * t; // seconds
	} else if (year >= 2150 && year <= 2999) {
		// End of Espenak range
		auto u = (y - 1820.0) / 100.0;
		result = -20.0 + 32.0 * u * u; // seconds
	} else {
		result = 60.0; // seconds
	}

	// Convert seconds to days
	return result / (24.0 * 3600.0);
}

///----------------------------------------
/// MARK: parallax_xyz
///----------------------------------------

void parallax_xyz(double wtime, double latitude,
                  double& x, double& y, double& z) noexcept {
	// X, Y, Z in AU. Parallax can be 8.8 arcsec per AU distance.
	// See Meeus, Astronomical Algorithms, p. 78.
	static constexpr double kAE = 149597870.700; // AU in km (IAU 2012)
	static constexpr double kHeightAboveSea = 100.0; // metres
	static constexpr double kFlatteningEarth = 0.99664719; // Earth is not perfectly round

	// Reduced latitude
	auto u = std::atan(kFlatteningEarth * std::sin(latitude) / std::cos(latitude));

	auto sin_lat_corr = kFlatteningEarth * std::sin(u)
	                   + kHeightAboveSea * std::sin(latitude) / 6378140.0;
	auto cos_lat_corr = std::cos(u)
	                   + kHeightAboveSea * std::cos(latitude) / 6378140.0;

	// Observer position in AU
	auto x_observ = (6378.14 / kAE) * cos_lat_corr * std::cos(wtime);
	auto y_observ = (6378.14 / kAE) * cos_lat_corr * std::sin(wtime);
	auto z_observ = (6378.14 / kAE) * sin_lat_corr;

	x -= x_observ;
	y -= y_observ;
	z -= z_observ;
}

///----------------------------------------
/// MARK: minor_planet
///----------------------------------------

void minor_planet(bool sun_earth_vector, double julian,
                  int year, int month, double day,
                  double a_e, double a_or_q, double a_i,
                  double a_ohm, double a_w, double a_M,
                  double wtime, double site_lat,
                  double& ra3, double& dec3,
                  double& delta, double& sun_delta,
                  bool& outdated,
                  astap::core::ephem::Vec3& ph_earth,
                  astap::core::ephem::Vec3& ph_pln) {
	namespace ephem = astap::core::ephem;
	static constexpr auto pi = std::numbers::pi;

	// Compute heliocentric Earth position if not already done.
	auto vh_earth = ephem::Vec3{};
	if (!sun_earth_vector) {
		const auto es = ephem::earth_state(julian, ephem::ReferenceFrame::Heliocentric);
		ph_earth = es.position;
		vh_earth = es.velocity;
	}

	// Build the element set. The unified `propagate` API handles
	// asteroids/comets/parabolic/hyperbolic uniformly via eccentricity.
	// The legacy API distinguished asteroid (JFORM=2) vs comet (JFORM=3)
	// explicitly; we preserve the behaviour via the "comet sentinel" below.
	ephem::OrbitalElements elements{};
	elements.epoch_jd              = astap::stacking::julian_calc(year, month, day, 0, 0, 0);
	elements.perihelion_distance   = a_or_q;      // AU
	elements.eccentricity          = a_e;
	elements.inclination           = a_i   * pi / 180.0;
	elements.ascending_node        = a_ohm * pi / 180.0;
	elements.argument_of_perihelion = a_w  * pi / 180.0;

	if (a_M < 1e98) {
		// Asteroid — semimajor axis was passed as a_or_q, but `propagate`
		// wants perihelion distance q. Convert: q = a * (1 - e).
		elements.perihelion_distance   = a_or_q * (1.0 - a_e);
		elements.mean_anomaly_at_epoch = a_M * pi / 180.0;
		const auto days_since_epoch = std::abs(julian - elements.epoch_jd);
		outdated = days_since_epoch > 120.0;
	} else {
		// Comet — epoch is perihelion passage, mean anomaly at epoch is zero.
		elements.mean_anomaly_at_epoch = 0.0;
		outdated = false;
	}

	const auto state = ephem::propagate(elements, julian);
	const auto& pos = state.position;
	const auto& vel = state.velocity;

	// Geometric distance minor-planet to Earth (AU).
	const auto dx = pos[0] - ph_earth[0];
	const auto dy = pos[1] - ph_earth[1];
	const auto dz = pos[2] - ph_earth[2];
	const auto r_geo = std::sqrt(dx * dx + dy * dy + dz * dz);

	// Light time in days (velocity from propagate is AU/day, so the
	// correction tl_days × vel gives an AU-valued displacement).
	const auto tl_days = kTauDays * r_geo;

	// Correct position for planetary aberration using velocity components.
	// The Earth position is already aberration-corrected via earth_state.
	auto x_pln = dx - tl_days * vel[0];
	auto y_pln = dy - tl_days * vel[1];
	auto z_pln = dz - tl_days * vel[2];

	// Observer parallax correction (should be in equinox of date; small J2000 error).
	parallax_xyz(wtime, site_lat, x_pln, y_pln, z_pln);

	// Convert to polar coordinates -> RA, Dec.
	const auto sph = ephem::cartesian_to_spherical({x_pln, y_pln, z_pln});
	delta = sph.radius;
	dec3  = sph.dec;
	ra3   = sph.ra;

	// Store heliocentric planet position for illumination calculation.
	ph_pln[0] = pos[0];
	ph_pln[1] = pos[1];
	ph_pln[2] = pos[2];

	// Heliocentric distance of the minor planet.
	sun_delta = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]);

	// propagate() always succeeds (the Kepler solver converges for all
	// physical inputs), so no status return is produced.
	(void)vh_earth; // velocity is unused today; kept for future Earth-velocity needs
}

///----------------------------------------
/// MARK: illum2
///----------------------------------------

void illum2(double x, double y, double z,
            double xe, double ye, double ze,
            double& r_sp, double& r_ep,
            double& elong, double& phi, double& phase) noexcept {
	static constexpr auto pi = std::numbers::pi;

	// Minor planet geocentric position
	auto xp = x - xe;
	auto yp = y - ye;
	auto zp = z - ze;

	// Distances in the Sun-Earth-planet triangle
	r_sp = std::sqrt(x * x + y * y + z * z);       // Sun -- planet
	auto re = std::sqrt(xe * xe + ye * ye + ze * ze); // Sun -- Earth
	r_ep = std::sqrt(xp * xp + yp * yp + zp * zp);   // Earth -- planet

	// Elongation
	elong = (180.0 / pi) * std::acos((r_ep * r_ep + re * re - r_sp * r_sp)
	                                  / (2.0 * r_ep * re));

	// Phase angle
	auto c_phi = (r_ep * r_ep + r_sp * r_sp - re * re) / (2.0 * r_ep * r_sp);
	phi = (180.0 / pi) * std::acos(c_phi); // degrees

	// Phase fraction 0..100
	phase = 100.0 * 0.5 * (1.0 + c_phi);
}

///----------------------------------------
/// MARK: asteroid_magn_comp
///----------------------------------------

[[nodiscard]] double asteroid_magn_comp(double g, double b) noexcept {
	// Meeus, Astronomical Algorithms, formula 32.14.
	// g = slope parameter, b = phase angle (radians, Sun-asteroid-Earth).
	auto b2 = std::sin(b * 0.5) / std::cos(b * 0.5); // tan(b/2)

	// Power via exp(exponent * ln(base)):
	//   q1 = exp(-3.33 * (tan(b/2))^0.63)
	//   q2 = exp(-1.87 * (tan(b/2))^1.22)
	auto q1 = std::exp(-3.33 * std::exp(0.63 * std::log(b2 + 1e-8)));
	auto q2 = std::exp(-1.87 * std::exp(1.22 * std::log(b2 + 1e-8)));

	return -2.5 * std::log((1.0 - g) * q1 + g * q2) / std::log(10.0);
}

///----------------------------------------
/// MARK: convert_MPCORB_line
///----------------------------------------

void convert_MPCORB_line(std::string_view txt,
                         std::string& desn, std::string& name,
                         int& yy, int& mm,
                         double& dd, double& a_e, double& a_a,
                         double& a_i, double& a_ohm, double& a_w,
                         double& a_M, double& h, double& g) {
	desn.clear(); // assume failure

	if (txt.size() < 195) {
		return;
	}

	// Epoch in packed form. Pascal: txt[21] is 1-based -> 0-based index 20.
	// Ord(txt[21]) - 55 gives the century prefix: 'J'=74-55=19, 'K'=75-55=20.
	auto century_val = static_cast<int>(txt[20]) - 55;
	auto century_str = std::to_string(century_val);

	if (century_str != "19" && century_str != "20" && century_str != "21") {
		return;
	}

	// Name: Pascal copy(txt,167,28) -> 0-based substr(166, 28)
	name = std::string(txt.substr(166, 28));

	// Designation: Pascal copy(txt,1,7) -> 0-based substr(0, 7), trimmed
	auto desn_sv = txt.substr(0, 7);
	auto end_pos = desn_sv.find_last_not_of(' ');
	desn = (end_pos != std::string_view::npos)
	     ? std::string(desn_sv.substr(0, end_pos + 1))
	     : std::string(desn_sv);

	// H: Pascal copy(txt,8,5) -> 0-based substr(7, 5)
	h = safe_stod(txt.substr(7, 5));

	// G: Pascal copy(txt,14,6) -> 0-based substr(13, 6)
	g = safe_stod(txt.substr(13, 6));

	// Epoch year: century string + txt[22] + txt[23] -> indices 21, 22
	auto year_str = century_str + std::string(1, txt[21]) + std::string(1, txt[22]);
	yy = safe_stoi(year_str);

	// Epoch month: Pascal txt[24] -> index 23
	auto code_month = static_cast<int>(txt[23]);
	if (code_month < 65) {
		mm = code_month - 48; // '0'..'9' -> 0..9
	} else {
		mm = code_month - 55; // 'A'..'Z' -> 10..35
	}

	// Epoch day: Pascal txt[25] -> index 24
	auto code_day = static_cast<int>(txt[24]);
	if (code_day < 65) {
		dd = static_cast<double>(code_day - 48);
	} else {
		dd = static_cast<double>(code_day - 55);
	}

	// Mean anomaly: Pascal copy(txt,27,9) -> 0-based substr(26, 9)
	a_M = safe_stod(txt.substr(26, 9));

	// Argument of perihelion: Pascal copy(txt,38,9) -> 0-based substr(37, 9)
	a_w = safe_stod(txt.substr(37, 9));

	// Longitude of ascending node: Pascal copy(txt,49,9) -> 0-based substr(48, 9)
	a_ohm = safe_stod(txt.substr(48, 9));

	// Inclination: Pascal copy(txt,60,9) -> 0-based substr(59, 9)
	a_i = safe_stod(txt.substr(59, 9));

	// Eccentricity: Pascal copy(txt,71,9) -> 0-based substr(70, 9)
	a_e = safe_stod(txt.substr(70, 9));

	// Semimajor axis: Pascal copy(txt,93,11) -> 0-based substr(92, 11)
	a_a = safe_stod(txt.substr(92, 11));
}

///----------------------------------------
/// MARK: convert_comet_line
///----------------------------------------

void convert_comet_line(std::string_view txt,
                        std::string& desn, std::string& name,
                        int& yy, int& mm,
                        double& dd, double& ecc, double& q,
                        double& inc2, double& lan, double& aop,
                        double& M_anom, double& H, double& k) {
	desn.clear(); // assume failure

	if (txt.size() < 168) {
		return;
	}

	// Epoch year: Pascal copy(txt,15,4) -> 0-based substr(14, 4)
	yy = safe_stoi(txt.substr(14, 4));

	if (yy <= 1900 || yy >= 2200) {
		return;
	}

	// Epoch month: Pascal copy(txt,20,2) -> 0-based substr(19, 2)
	mm = safe_stoi(txt.substr(19, 2));

	// Epoch day: Pascal copy(txt,23,7) -> 0-based substr(22, 7)
	dd = safe_stod(txt.substr(22, 7));

	// Perihelion distance q: Pascal copy(txt,31,9) -> 0-based substr(30, 9)
	q = safe_stod(txt.substr(30, 9));

	// Eccentricity: Pascal copy(txt,41,9) -> 0-based substr(40, 9)
	ecc = safe_stod(txt.substr(40, 9));

	// Argument of perihelion: Pascal copy(txt,51,9) -> 0-based substr(50, 9)
	aop = safe_stod(txt.substr(50, 9));

	// Longitude of ascending node: Pascal copy(txt,61,9) -> 0-based substr(60, 9)
	lan = safe_stod(txt.substr(60, 9));

	// Inclination: Pascal copy(txt,71,9) -> 0-based substr(70, 9)
	inc2 = safe_stod(txt.substr(70, 9));

	// Mean anomaly sentinel: flags this as a comet
	M_anom = 1e99;

	// Absolute magnitude H: Pascal copy(txt,91,5) -> 0-based substr(90, 5)
	H = safe_stod(txt.substr(90, 5));

	// Slope/activity g -> k: Pascal copy(txt,97,4) -> 0-based substr(96, 4)
	auto g_val = safe_stod(txt.substr(96, 4));
	k = g_val * 2.5; // comet activity

	// Name: Pascal copy(txt,103,28) -> 0-based substr(102, 28)
	name = std::string(txt.substr(102, 28));

	// Designation: Pascal copy(txt,160,9) -> 0-based substr(159, 9)
	desn = std::string(txt.substr(159, 9));
}

} // namespace
