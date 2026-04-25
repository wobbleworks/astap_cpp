///----------------------------------------
///      @file asteroid_overlay.cpp
///   @ingroup ASTAP++
///     @brief Implementation of @ref scan_asteroids_in_field.
///    @author Created by John Stephen on 4/23/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "asteroid_overlay.h"

#include "asteroid.h"
#include "../core/ephemerides.h"
#include "../core/globals.h"
#include "../stacking/stack.h"

#include <cmath>
#include <fstream>
#include <numbers>
#include <string>

///----------------------------------------
namespace astap::analysis {
///----------------------------------------

namespace {

constexpr auto kPi = std::numbers::pi;

// Constants from Meeus (new editions, Ch. 11 / Ch. 83):
// sidereal time at 2000 Jan 1.5 UT and Earth's angular velocity per day.
constexpr auto kSiderealTime2000 = 280.46061837 * kPi / 180.0;
constexpr auto kEarthAngularVelocity = 2.0 * kPi * 1.00273790935;

// Positive-modulo wrap, matching Pascal's fnmodulo.
[[nodiscard]] double positive_mod(double x, double y) noexcept {
	auto r = std::fmod(x, y);
	if (r < 0) r += y;
	return r;
}

// Scan one catalog (asteroids or comets). Returns early when max_count is
// reached. Shares @p earth_cached across iterations to skip recomputing
// the Earth heliocentric vector for every line.
void scan_file(const std::filesystem::path& path, bool is_comet,
               const Header& head,
               double jd_tt, double fov_sqr,
               double wtime, double site_lat,
               double max_magnitude, int remaining,
               const std::function<void(const std::string&)>& log,
               std::vector<AsteroidDetection>& out) {
	std::ifstream in(path);
	if (!in.is_open()) {
		if (log) log("Cannot open " + path.string());
		return;
	}

	// Cached earth state + planet state. minor_planet writes to these;
	// passing sun_earth_vector=true on subsequent calls skips the earth
	// recompute but still produces correct planet positions.
	astap::core::ephem::Vec3 ph_earth{};
	astap::core::ephem::Vec3 ph_pln{};
	auto sun_earth_valid = false;

	auto line = std::string{};
	auto scanned = 0;
	while (remaining > 0 && std::getline(in, line)) {
		if (line.size() < 10) continue;

		auto desn = std::string{};
		auto name = std::string{};
		auto yy = 0, mm = 0;
		auto dd = 0.0, a_e = 0.0, a_or_q = 0.0, a_i = 0.0, a_ohm = 0.0,
		     a_w = 0.0, a_M = 0.0, H = 0.0, g_or_k = 0.0;

		if (is_comet) {
			convert_comet_line(line, desn, name, yy, mm, dd,
			                   a_e, a_or_q, a_i, a_ohm, a_w, a_M, H, g_or_k);
		} else {
			convert_MPCORB_line(line, desn, name, yy, mm, dd,
			                    a_e, a_or_q, a_i, a_ohm, a_w, a_M, H, g_or_k);
		}
		if (desn.empty() || a_or_q == 0.0) continue;
		--remaining;
		++scanned;

		auto ra = 0.0, dec = 0.0, delta = 0.0, sun_delta = 0.0;
		auto outdated = false;
		try {
			minor_planet(sun_earth_valid, jd_tt, yy, mm, dd, a_e, a_or_q,
			             a_i, a_ohm, a_w, a_M, wtime, site_lat,
			             ra, dec, delta, sun_delta, outdated,
			             ph_earth, ph_pln);
			sun_earth_valid = true;
		} catch (...) {
			continue;   // malformed element — skip this row
		}

		// FOV gate: Euclidean distance in (Δra·cos(dec), Δdec) space,
		// squared. Matches Pascal's sqr(dra*cos) + sqr(ddec) < sqr(fov).
		const auto cos_dec0 = std::cos(head.dec0);
		const auto dra = (ra - head.ra0) * cos_dec0;
		const auto ddec = dec - head.dec0;
		if (dra * dra + ddec * ddec >= fov_sqr) continue;

		// Magnitude: asteroid H-G vs comet H-k (see Meeus; comet mag
		// depends only on heliocentric distance to the power k).
		const auto log10 = [](double v) { return std::log(v) / std::log(10.0); };
		auto mag = 0.0;
		if (is_comet) {
			mag = H + 5.0 * log10(delta) + g_or_k * log10(sun_delta);
		} else {
			mag = H + 5.0 * log10(delta * sun_delta);
			// Phase correction. illum2 gives us the phase angle; asteroid_magn_comp
			// is additive. The phase angle must come from the same
			// (earth, planet) vectors minor_planet just computed.
			auto r_sp = 0.0, r_ep = 0.0, elong = 0.0, phi = 0.0, phase = 0.0;
			illum2(ph_pln[0], ph_pln[1], ph_pln[2],
			       ph_earth[0], ph_earth[1], ph_earth[2],
			       r_sp, r_ep, elong, phi, phase);
			mag += asteroid_magn_comp(g_or_k, phi * kPi / 180.0);
		}
		if (mag > max_magnitude) continue;

		auto det = AsteroidDetection{};
		det.desn      = std::move(desn);
		det.name      = std::move(name);
		det.ra        = ra;
		det.dec       = dec;
		det.magnitude = mag;
		det.is_comet  = is_comet;
		det.outdated  = outdated;
		out.push_back(std::move(det));

		if (log && (scanned % 10000) == 0) {
			log("Scanned " + std::to_string(scanned) + " entries…");
		}
	}
}

}  // namespace

std::vector<AsteroidDetection> scan_asteroids_in_field(
		const Header& head, const AsteroidScanOptions& opts) {
	auto out = std::vector<AsteroidDetection>{};

	if (head.naxis == 0 || head.cd1_1 == 0.0) {
		if (opts.log) opts.log("Plate-solve the image first.");
		return out;
	}

	// Populate jd_start / jd_mid / jd_end from DATE-OBS etc.
	astap::stacking::date_to_jd(head.date_obs, head.date_avg, head.exposure);
	if (astap::jd_start <= 2400000.0) {
		if (opts.log) opts.log("Cannot parse DATE-OBS / DATE-AVG.");
		return out;
	}

	const auto jd_tt = astap::jd_mid + deltaT_calc(astap::jd_mid);

	// Local sidereal time (radians). Mirrors Pascal unit_asteroid line ~722.
	// East-positive longitude convention per ESA SSA.
	const auto wtime = positive_mod(
		astap::site_long_radians + kSiderealTime2000 +
		(astap::jd_mid - 2451545.0) * kEarthAngularVelocity,
		2.0 * kPi);

	// Field-of-view radius (radians), squared for cheap gate.
	// 1.5x oversize matches Pascal: catches objects just off-frame in case
	// the WCS isn't pixel-perfect at the edges.
	const auto half_w_deg = 0.5 * head.width  * std::abs(head.cdelt1);
	const auto half_h_deg = 0.5 * head.height * std::abs(head.cdelt2);
	const auto fov_rad = 1.5 * std::sqrt(half_w_deg * half_w_deg +
	                                      half_h_deg * half_h_deg) * kPi / 180.0;
	const auto fov_sqr = fov_rad * fov_rad;

	if (!opts.mpcorb_path.empty()) {
		scan_file(opts.mpcorb_path, /*is_comet=*/false, head,
		          jd_tt, fov_sqr, wtime, astap::site_lat_radians,
		          opts.max_magnitude, opts.max_count, opts.log, out);
	}
	if (!opts.comets_path.empty()) {
		scan_file(opts.comets_path, /*is_comet=*/true, head,
		          jd_tt, fov_sqr, wtime, astap::site_lat_radians,
		          opts.max_magnitude, opts.max_count, opts.log, out);
	}

	if (opts.log) opts.log("Found " + std::to_string(out.size()) +
	                       " objects in field.");
	return out;
}

} // namespace
