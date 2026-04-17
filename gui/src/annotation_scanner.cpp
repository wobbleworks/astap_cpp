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

} // namespace astap::gui
