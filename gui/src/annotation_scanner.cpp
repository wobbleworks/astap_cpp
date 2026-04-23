///----------------------------------------
///      @file annotation_scanner.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the deep-sky annotation scanner.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "annotation_scanner.h"

#include "../../src/analysis/constellations.h"
#include "../../src/core/wcs.h"
#include "../../src/reference/online_gaia.h"
#include "../../src/reference/star_database.h"
#include "../../src/reference/stars_wide_field.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <numbers>
#include <string>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

namespace {

constexpr double kPi = std::numbers::pi;

// Module-level catalog storage. Loaded once; queried per-solve.
std::vector<std::string> g_catalog;

struct CatalogEntry {
	double ra = 0.0;     // radians
	double dec = 0.0;
	std::string name;
	double lengthAm = 0.0; // major axis, arcminutes
	double widthAm = 0.0;  // minor axis, arcminutes
	double paDeg = 999.0;  // position angle, degrees
};

std::vector<CatalogEntry> g_entries;

// Parse one CSV line into a CatalogEntry. Returns false if the line is
// a header or malformed.
[[nodiscard]] bool parse_line(const std::string& line, CatalogEntry& out) {
	// Format: RA_encoded,DEC_encoded,name(s),length,width,PA
	// RA is in units of 0.1 arcsec (0..864000 → 0..24h)
	// Dec is in units of 1 arcsec (-324000..+324000 → -90..+90°)
	if (line.empty() || line[0] < '0' || line[0] > '9') {
		return false;
	}

	auto pos = std::size_t{0};
	auto next_field = [&]() -> std::string_view {
		const auto start = pos;
		const auto comma = line.find(',', pos);
		if (comma == std::string::npos) {
			pos = line.size();
			return {line.data() + start, line.size() - start};
		}
		pos = comma + 1;
		return {line.data() + start, comma - start};
	};

	// RA (encoded as integer × 0.1 arcsec, range 0..864000)
	const auto ra_str = next_field();
	char* end1 = nullptr;
	const auto ra_enc = std::strtod(std::string(ra_str).c_str(), &end1);
	if (end1 == std::string(ra_str).c_str()) {
		return false;
	}
	out.ra = ra_enc * 0.1 / 3600.0 * (kPi / 12.0);

	// Dec (encoded as integer × 1 arcsec, range -324000..+324000)
	const auto dec_str = next_field();
	char* end2 = nullptr;
	const auto dec_enc = std::strtod(std::string(dec_str).c_str(), &end2);
	if (end2 == std::string(dec_str).c_str()) {
		return false;
	}
	out.dec = dec_enc / 3600.0 * (kPi / 180.0);

	// Name(s) — slash-separated; take the first and replace underscores
	const auto names = next_field();
	const auto slash = names.find('/');
	auto name = (slash != std::string_view::npos) ? names.substr(0, slash) : names;
	out.name = std::string(name);
	for (auto& c : out.name) {
		if (c == '_') {
			c = ' ';
		}
	}

	// Length (0.1 arcmin)
	const auto len_str = next_field();
	if (!len_str.empty()) {
		out.lengthAm = std::strtod(std::string(len_str).c_str(), nullptr) * 0.1;
	}

	// Width (0.1 arcmin)
	const auto wid_str = next_field();
	if (!wid_str.empty()) {
		out.widthAm = std::strtod(std::string(wid_str).c_str(), nullptr) * 0.1;
	}

	// PA (degrees; 999 = unknown)
	const auto pa_str = next_field();
	if (!pa_str.empty()) {
		out.paDeg = std::strtod(std::string(pa_str).c_str(), nullptr);
	}

	return true;
}

} // namespace

bool load_deepsky_catalog(const std::filesystem::path& dir) {
	auto path = dir / "deep_sky.csv";
	std::ifstream in(path);
	if (!in.is_open()) {
		return false;
	}

	g_catalog.clear();
	g_entries.clear();

	std::string line;
	while (std::getline(in, line)) {
		CatalogEntry entry;
		if (parse_line(line, entry)) {
			g_entries.push_back(std::move(entry));
		}
	}

	return g_entries.size() > 100;
}

std::vector<AnnotationMarker> annotate_image(const astap::Header& head) {
	auto result = std::vector<AnnotationMarker>{};

	if (g_entries.empty() || head.cdelt2 == 0.0) {
		return result;
	}

	// Image centre RA/Dec
	const auto cx = (head.width + 1.0) * 0.5;
	const auto cy = (head.height + 1.0) * 0.5;
	double ra0 = 0.0;
	double dec0 = 0.0;
	astap::core::pixel_to_celestial(head, cx, cy, /*formalism=*/0, ra0, dec0);

	// Diagonal FOV radius in radians (generous, to catch edge objects).
	// cdelt2 is in degrees per pixel.
	const auto fovRad = std::sqrt(
		static_cast<double>(head.width) * head.width +
		static_cast<double>(head.height) * head.height)
		* std::abs(head.cdelt2) * (kPi / 180.0) * 0.6;

	// Pixel scale for converting arcminute axes to pixels.
	// cdelt2 is degrees/pixel → arcsec/pixel = cdelt2 × 3600.
	const auto arcsecPerPx = std::abs(head.cdelt2) * 3600.0;

	for (const auto& e : g_entries) {
		// Quick angular-distance rejection.
		double sep = 0.0;
		astap::core::ang_sep(e.ra, e.dec, ra0, dec0, sep);
		if (sep > fovRad) {
			continue;
		}

		// Project to FITS pixel coordinates.
		double px = 0.0;
		double py = 0.0;
		astap::core::celestial_to_pixel(head, e.ra, e.dec, px, py);

		// Reject if outside image bounds (with small margin for labels).
		if (px < -50 || px > head.width + 50 || py < -50 || py > head.height + 50) {
			continue;
		}

		AnnotationMarker m;
		m.name = QString::fromStdString(e.name);
		m.x = px;
		m.y = py;
		m.majorPx = (e.lengthAm * 60.0) / arcsecPerPx * 0.5;
		m.minorPx = (e.widthAm > 0 ? e.widthAm : e.lengthAm) * 60.0 / arcsecPerPx * 0.5;
		m.paDeg = (e.paDeg < 999) ? e.paDeg : 0.0;
		result.push_back(std::move(m));
	}

	return result;
}

ConstellationOverlay build_constellation_overlay(const astap::Header& head) {
	ConstellationOverlay out;

	if (head.cdelt2 == 0.0 || head.width <= 0 || head.height <= 0) {
		return out;
	}

	// Image centre for angular-distance pre-filter.
	const auto cx = (head.width + 1.0) * 0.5;
	const auto cy = (head.height + 1.0) * 0.5;
	double ra0 = 0.0;
	double dec0 = 0.0;
	astap::core::pixel_to_celestial(head, cx, cy, /*formalism=*/0, ra0, dec0);

	const auto fovRad = std::sqrt(
		static_cast<double>(head.width) * head.width +
		static_cast<double>(head.height) * head.height)
		* std::abs(head.cdelt2) * (kPi / 180.0) * 0.7;

	const auto margin = head.width * 0.1;

	auto projectStar = [&](const astap::analysis::ConstStar& s,
	                       double& px, double& py) -> bool {
		const auto ra  = s.ra  / 1000.0 * (kPi / 12.0);
		const auto dec = s.dec / 100.0  * (kPi / 180.0);
		double sep = 0.0;
		astap::core::ang_sep(ra, dec, ra0, dec0, sep);
		if (sep > fovRad) {
			return false;
		}
		astap::core::celestial_to_pixel(head, ra, dec, px, py);
		return (px > -margin && px < head.width + margin &&
		        py > -margin && py < head.height + margin);
	};

	// Stick-figure lines: dm=-2 starts a new strip, dm=-1 draws to next.
	using astap::analysis::kDrawModeStart;
	using astap::analysis::kDrawModeDraw;
	const auto& cons = astap::analysis::constellation;

	auto prevX = 0.0;
	auto prevY = 0.0;
	auto inStrip = false;

	for (std::size_t i = 0; i < cons.size(); ++i) {
		double px = 0.0;
		double py = 0.0;
		const auto visible = projectStar(cons[i], px, py);

		if (cons[i].dm == kDrawModeStart) {
			// Start a new strip — remember position, no line yet.
			prevX = px;
			prevY = py;
			inStrip = visible;
		} else if (cons[i].dm == kDrawModeDraw) {
			// Draw line from previous to current if either end is visible.
			if (inStrip || visible) {
				out.lines.push_back({prevX, prevY, px, py});
			}
			prevX = px;
			prevY = py;
			inStrip = visible;
		}
	}

	// Constellation name labels.
	const auto& pos = astap::analysis::constPos;
	const auto& names = astap::analysis::constShortName;

	for (std::size_t i = 0; i < pos.size(); ++i) {
		const auto ra  = pos[i][0] / 1000.0 * (kPi / 12.0);
		const auto dec = pos[i][1] / 100.0  * (kPi / 180.0);
		double sep = 0.0;
		astap::core::ang_sep(ra, dec, ra0, dec0, sep);
		if (sep > fovRad) {
			continue;
		}
		double px = 0.0;
		double py = 0.0;
		astap::core::celestial_to_pixel(head, ra, dec, px, py);
		if (px > -margin && px < head.width + margin &&
		    py > -margin && py < head.height + margin) {
			out.labels.push_back({
				QString::fromLatin1(names[i].data(),
					static_cast<qsizetype>(names[i].size())),
				px, py});
		}
	}

	return out;
}

///----------------------------------------
/// MARK: Star-catalog overlay
///
/// Walks the active catalog (local .1476/.290 tiles, wide-field w08, or
/// online Gaia — whichever `select_star_database` has picked) and projects
/// each entry into pixel space. Mirrors the iteration logic in
/// src/core/photometry_catalog.cpp, but emits overlay markers instead of
/// running HFD at each projected position.
///----------------------------------------

namespace {

constexpr auto kDeg2Rad  = kPi / 180.0;
constexpr auto kTile1476 = 5.142857143 * kDeg2Rad;
constexpr auto kTile290  = 9.53 * kDeg2Rad;

/// @brief Project (ra, dec) + append a marker if it lands in a reasonable
///        bbox around the frame (with 50-pixel slop).
void project_and_emit(const astap::Header& head,
                       double ra, double dec,
                       double magn_times_10, double bp_rp,
                       std::vector<CatalogStarMarker>& out) {
	auto fitsX = 0.0;
	auto fitsY = 0.0;
	astap::core::celestial_to_pixel(head, ra, dec, fitsX, fitsY);
	if (fitsX < -50.0 || fitsX > head.width  + 50.0
			|| fitsY < -50.0 || fitsY > head.height + 50.0) {
		return;
	}
	out.push_back({
		.x    = fitsX,
		.y    = fitsY,
		.magn = magn_times_10 * 0.1,   // catalog stores mag*10
		.bpRp = (bp_rp == 999.0) ? 999.0 : bp_rp * 0.1,
	});
}

}  // namespace

std::vector<CatalogStarMarker> scan_catalog_stars(
		const astap::Header& head, int max_markers) {
	std::vector<CatalogStarMarker> markers;
	if (head.naxis == 0 || head.cd1_1 == 0.0 || max_markers <= 0) {
		return markers;
	}

	auto telescope_ra  = 0.0;
	auto telescope_dec = 0.0;
	astap::core::pixel_to_celestial(head,
	                                 (head.width  + 1) * 0.5,
	                                 (head.height + 1) * 0.5,
	                                 /*formalism=*/1,
	                                 telescope_ra, telescope_dec);
	const auto fov_org = std::sqrt(
		(head.width  * head.cdelt1) * (head.width  * head.cdelt1) +
		(head.height * head.cdelt2) * (head.height * head.cdelt2)) * kDeg2Rad;

	const auto db_type = astap::reference::database_type;
	markers.reserve(static_cast<std::size_t>(std::min(max_markers, 2048)));

	if (db_type == 1476 || db_type == 290) {
		// Tile-based: select up to 4 cells covering the FOV, stream each.
		const auto tile_cap = (db_type == 1476) ? kTile1476 : kTile290;
		const auto fov = std::min(fov_org, tile_cap);

		auto area1 = 0, area2 = 0, area3 = 0, area4 = 0;
		auto f1 = 0.0, f2 = 0.0, f3 = 0.0, f4 = 0.0;
		astap::reference::find_areas(telescope_ra, telescope_dec, fov,
		                              area1, area2, area3, area4,
		                              f1, f2, f3, f4);
		const std::array<int, 4> areas{area1, area2, area3, area4};

		for (const auto area : areas) {
			if (area == 0 || static_cast<int>(markers.size()) >= max_markers) {
				continue;
			}
			if (!astap::reference::open_database(telescope_dec, area)) {
				continue;
			}
			auto ra = 0.0, dec = 0.0, mag_x10 = 0.0, bp_rp = 999.0;
			while (static_cast<int>(markers.size()) < max_markers
				&& astap::reference::readdatabase290(
					telescope_ra, telescope_dec, fov, ra, dec, mag_x10, bp_rp)) {
				project_and_emit(head, ra, dec, mag_x10, bp_rp, markers);
			}
			astap::reference::close_star_database();
		}
	}
	else if (db_type == 1) {
		// Wide-field (w08): iterate the flat triples buffer.
		if (astap::reference::wide_database != astap::reference::name_database) {
			if (!astap::reference::read_stars_wide_field()) {
				return markers;
			}
		}
		const auto& buf = astap::reference::wide_field_stars;
		const auto n = buf.size() / 3;
		const auto sep_limit = std::min(
			fov_org * 0.5 * 0.9 * (2.0 / std::sqrt(kPi)),
			kPi * 0.5);
		for (std::size_t i = 0; i < n && static_cast<int>(markers.size()) < max_markers; ++i) {
			const auto mag_x10 = static_cast<double>(buf[i * 3 + 0]);
			const auto ra      = static_cast<double>(buf[i * 3 + 1]);
			const auto dec     = static_cast<double>(buf[i * 3 + 2]);
			auto sep = 0.0;
			astap::core::ang_sep(ra, dec, telescope_ra, telescope_dec, sep);
			if (sep < sep_limit) {
				project_and_emit(head, ra, dec, mag_x10, /*bp_rp=*/999.0, markers);
			}
		}
	}
	else {
		// Online Gaia: caller is responsible for having populated the table.
		const auto& odb = astap::reference::online_database;
		if (odb[0].empty()) {
			return markers;
		}
		for (std::size_t i = 0; i < odb[0].size() && static_cast<int>(markers.size()) < max_markers; ++i) {
			const auto magf = odb[5][i];
			if (magf == 0.0) continue;
			project_and_emit(head, odb[0][i], odb[1][i],
				/*mag_x10=*/magf * 10.0, /*bp_rp=*/999.0, markers);
		}
	}

	// Brightest first so the viewer can draw big-to-small without overlap fights.
	std::sort(markers.begin(), markers.end(),
		[](const auto& a, const auto& b) { return a.magn < b.magn; });
	return markers;
}

///----------------------------------------
/// MARK: - Online overlays (VSX/VSP, Simbad, Vizier)
///----------------------------------------

namespace {

// Reject markers that fall well outside the frame. Mirrors the slop used by
// the deep-sky and catalog-star overlays above.
constexpr double kOverlaySlopPx = 50.0;

[[nodiscard]] bool inFrame(const astap::Header& head, double px, double py) noexcept {
	return px >= -kOverlaySlopPx && px <= head.width  + kOverlaySlopPx
	    && py >= -kOverlaySlopPx && py <= head.height + kOverlaySlopPx;
}

}  // namespace

std::vector<VarStarMarker> project_var_stars(const astap::Header& head) {
	auto out = std::vector<VarStarMarker>{};
	if (head.naxis == 0 || head.cd1_1 == 0.0) {
		return out;
	}

	out.reserve(astap::core::vsx.size() + astap::core::vsp.size());

	for (const auto& v : astap::core::vsx) {
		auto px = 0.0, py = 0.0;
		astap::core::celestial_to_pixel(head, v.ra, v.dec, px, py);
		if (!inFrame(head, px, py)) continue;
		auto m = VarStarMarker{};
		m.name = QString::fromStdString(v.name);
		// Display range "max → min" if both present, else whichever exists.
		if (!v.maxmag.empty() && !v.minmag.empty()) {
			m.magText = QString::fromStdString(v.maxmag + "-" + v.minmag);
		} else if (!v.maxmag.empty()) {
			m.magText = QString::fromStdString(v.maxmag);
		} else if (!v.minmag.empty()) {
			m.magText = QString::fromStdString(v.minmag);
		}
		m.x = px;
		m.y = py;
		m.isComparison = false;
		out.push_back(std::move(m));
	}

	for (const auto& c : astap::core::vsp) {
		auto px = 0.0, py = 0.0;
		astap::core::celestial_to_pixel(head, c.ra, c.dec, px, py);
		if (!inFrame(head, px, py)) continue;
		auto m = VarStarMarker{};
		m.name = QString::fromStdString(c.auid);
		// Prefer V, else B, else any band that's present (skip "?" placeholders).
		auto pick = [](const std::string& s) {
			return (!s.empty() && s != "?") ? QString::fromStdString(s) : QString{};
		};
		m.magText = pick(c.Vmag);
		if (m.magText.isEmpty()) m.magText = pick(c.Bmag);
		if (m.magText.isEmpty()) m.magText = pick(c.Rmag);
		m.x = px;
		m.y = py;
		m.isComparison = true;
		out.push_back(std::move(m));
	}
	return out;
}

std::vector<SimbadMarker> project_simbad(
		const std::vector<astap::core::SimbadObject>& objects,
		const astap::Header& head) {
	auto out = std::vector<SimbadMarker>{};
	if (head.naxis == 0 || head.cd1_1 == 0.0) {
		return out;
	}
	out.reserve(objects.size());

	// Simbad encodes ra_units = round(ra_hours * 864000 / 24)
	//   → ra_hours = ra_units * 24 / 864000 = ra_units / 36000
	//   → ra_rad   = ra_hours * pi / 12
	// And dec_units = round(sign * dec_deg * 324000 / 90)
	//   → dec_deg  = dec_units * 90 / 324000 = dec_units / 3600
	//   → dec_rad  = dec_deg * pi / 180
	for (const auto& o : objects) {
		const auto ra_rad  = static_cast<double>(o.ra_units)  / 36000.0 * (kPi / 12.0);
		const auto dec_rad = static_cast<double>(o.dec_units) / 3600.0  * (kPi / 180.0);
		auto px = 0.0, py = 0.0;
		astap::core::celestial_to_pixel(head, ra_rad, dec_rad, px, py);
		if (!inFrame(head, px, py)) continue;

		auto m = SimbadMarker{};
		m.name = QString::fromStdString(o.name);
		m.type = QString::fromStdString(o.type);
		m.x = px;
		m.y = py;
		m.magnitude = o.magnitude;
		m.sizeArcmin = o.size_arcmin;
		out.push_back(std::move(m));
	}
	return out;
}

std::vector<VizierMarker> project_vizier(
		const std::vector<astap::core::VizierObject>& rows,
		const astap::Header& head) {
	auto out = std::vector<VizierMarker>{};
	if (head.naxis == 0 || head.cd1_1 == 0.0) {
		return out;
	}
	out.reserve(rows.size());

	// Vizier encodes ra_units = round(ra_deg * 864000 / 360) = ra_deg * 2400
	//   → ra_rad  = ra_units / 2400 * pi / 180
	// dec_units = round(dec_deg * 324000 / 90) = dec_deg * 3600
	//   → dec_rad = ra_units / 3600 * pi / 180
	for (const auto& r : rows) {
		const auto ra_rad  = static_cast<double>(r.ra_units)  / 2400.0 * (kPi / 180.0);
		const auto dec_rad = static_cast<double>(r.dec_units) / 3600.0 * (kPi / 180.0);
		auto px = 0.0, py = 0.0;
		astap::core::celestial_to_pixel(head, ra_rad, dec_rad, px, py);
		if (!inFrame(head, px, py)) continue;
		out.push_back({.x = px, .y = py, .magnitude = r.magnitude});
	}

	// Brightest first — the painter draws big circles first so labels read
	// cleanly without smaller overlapping markers covering them.
	std::sort(out.begin(), out.end(),
		[](const auto& a, const auto& b) { return a.magnitude < b.magnitude; });
	return out;
}

} // namespace astap::gui
