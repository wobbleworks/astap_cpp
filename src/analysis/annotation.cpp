/// @file annotation.cpp
/// Image annotation utilities — ported from unit_annotation.pas.

#include "annotation.h"

#include "../core/from_chars_fp.h"   // astap::core::from_chars (portable FP parse)
#include "../core/wcs.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace astap::analysis {

// ---------------------------------------------------------------------------
// annotation_to_array
// ---------------------------------------------------------------------------

void annotation_to_array(std::string_view text, bool transparent, int colour,
                         int size, int x, int y, ImageArray& img)
{
	if (img.empty() || img[0].empty() || img[0][0].empty()) return;

	const int h = static_cast<int>(img[0].size());
	const int w = static_cast<int>(img[0][0].size());

	const int len = static_cast<int>(text.size());
	for (int k = 0; k < len; ++k) {
		const int value = static_cast<unsigned char>(text[k]);
		if (value < 33 || value > 126) continue;

		const int glyph = value - 33;

		for (int j = 9 * size - 1; j >= 0; --j) {
			for (int i = 0; i < 5 * size; ++i) {
				const int x2 = x + i + k * 7 * size;
				const int y2 = y - j;

				if (x2 < 0 || x2 >= w || y2 < 0 || y2 >= h) continue;

				const uint8_t pixel = kFont5x9[glyph][j / size][i / size];
				if (transparent && pixel == 0) continue;

				img[0][y2][x2] = static_cast<float>(pixel * colour);
			}
		}
	}
}

// ---------------------------------------------------------------------------
// find_object
// ---------------------------------------------------------------------------

namespace {

/// Parse a slash-delimited name field into up to three names.
struct ParsedNames {
	std::string_view naam2;
	std::string_view naam3;
	std::string_view naam4;
};

ParsedNames parse_names(std::string_view field)
{
	// Trim leading spaces.
	while (!field.empty() && field.front() == ' ')
		field.remove_prefix(1);

	ParsedNames n;
	const auto slash1 = field.find('/');
	if (slash1 == std::string_view::npos) {
		n.naam2 = field;
		return n;
	}
	n.naam2 = field.substr(0, slash1);
	const auto rest = field.substr(slash1 + 1);
	const auto slash2 = rest.find('/');
	if (slash2 == std::string_view::npos) {
		n.naam3 = rest;
	} else {
		n.naam3 = rest.substr(0, slash2);
		n.naam4 = rest.substr(slash2 + 1);
	}
	return n;
}

/// Case-insensitive ASCII comparison.
bool iequals(std::string_view a, std::string_view b)
{
	if (a.size() != b.size()) return false;
	for (size_t i = 0; i < a.size(); ++i) {
		if (std::toupper(static_cast<unsigned char>(a[i])) !=
		    std::toupper(static_cast<unsigned char>(b[i])))
			return false;
	}
	return true;
}

/// Trim trailing spaces from a string_view.
std::string_view rtrim(std::string_view s)
{
	while (!s.empty() && s.back() == ' ')
		s.remove_suffix(1);
	return s;
}

/// Extract the next comma-delimited field from @p line starting at @p pos.
/// Advances @p pos past the delimiter.
std::string_view next_field(std::string_view line, size_t& pos)
{
	if (pos >= line.size()) return {};
	const auto start = pos;
	const auto comma = line.find(',', pos);
	if (comma == std::string_view::npos) {
		pos = line.size();
		return rtrim(line.substr(start));
	}
	pos = comma + 1;
	return rtrim(line.substr(start, comma - start));
}

/// Parse an integer from a string_view, returning 0 on failure.
int parse_int(std::string_view s)
{
	int val = 0;
	std::from_chars(s.data(), s.data() + s.size(), val);
	return val;
}

/// Parse a double from a string_view, returning 0.0 on failure.
double parse_double(std::string_view s)
{
	double val = 0.0;
	astap::core::from_chars(s.data(), s.data() + s.size(), val);
	return val;
}

}  // anonymous namespace

std::optional<DeepSkyMatch> find_object(std::string_view object_name,
                                        std::span<const std::string> database_lines)
{
	if (object_name.size() <= 1) return std::nullopt;

	// Database lines start at index 2 (skip header lines), matching the
	// Pascal linepos:=2 convention.
	for (size_t row = 2; row < database_lines.size(); ++row) {
		const std::string_view line = database_lines[row];
		size_t pos = 0;

		// Field 1: RA encoded as integer (unit = 0.1 arcsec of time).
		const auto ra_field = next_field(line, pos);
		const int ra_int = parse_int(ra_field);
		if (ra_field.empty()) continue;  // malformed line

		// Field 2: Dec encoded as integer (unit = 1 arcsec).
		const auto dec_field = next_field(line, pos);
		const int dec_int = parse_int(dec_field);

		const double ra  = ra_int  * std::numbers::pi * 2.0 / 864000.0;
		const double dec = dec_int * std::numbers::pi * 0.5 / 324000.0;

		// Field 3: name(s), slash-separated.
		const auto name_field = next_field(line, pos);
		const auto names = parse_names(name_field);

		if (!iequals(object_name, names.naam2) &&
		    !iequals(object_name, names.naam3) &&
		    !iequals(object_name, names.naam4))
			continue;

		// Field 4: length (major axis).
		const auto length_field = next_field(line, pos);
		const double length = parse_double(length_field);

		// Field 5: width (minor axis).
		const auto width_field = next_field(line, pos);
		const double width = parse_double(width_field);

		// Field 6: position angle.
		const auto pa_field = next_field(line, pos);
		double pa = 999.0;
		if (!pa_field.empty()) {
			const double v = parse_double(pa_field);
			// from_chars returns 0.0 on failure — but PA=0 is valid,
			// so we check whether parsing actually consumed input.
			double tmp = 0.0;
			auto [ptr, ec] = astap::core::from_chars(pa_field.data(),
			                                         pa_field.data() + pa_field.size(),
			                                         tmp);
			if (ec == std::errc{})
				pa = v;
		}

		// Build canonical name.
		DeepSkyMatch match;
		match.ra     = ra;
		match.dec    = dec;
		match.length = length;
		match.width  = width;
		match.pa     = pa;

		if (names.naam3.empty()) {
			match.name = std::string(names.naam2);
		} else {
			match.name.reserve(names.naam2.size() + 1 + names.naam3.size());
			match.name += names.naam2;
			match.name += '_';
			match.name += names.naam3;
		}

		return match;
	}

	return std::nullopt;
}

// ---------------------------------------------------------------------------
// DeepSkyDatabase — catalog loading
// ---------------------------------------------------------------------------

namespace {

/// On-disk file name and required header version for a catalog. An empty
/// required_version means the catalog carries no version gate (Pascal only
/// checks the magnitude-13 and -15 variable-star files).
struct CatalogSpec {
	std::string_view file_name;
	std::string_view required_version;
};

CatalogSpec catalog_spec(DeepSkyCatalog catalog)
{
	switch (catalog) {
		case DeepSkyCatalog::DeepSky:    return {"deep_sky.csv", {}};
		case DeepSkyCatalog::HyperLeda:  return {"hyperleda.csv", {}};
		case DeepSkyCatalog::Variable:   return {"variable_stars.csv", {}};
		case DeepSkyCatalog::Variable13: return {"variable_stars_13.csv", "V003"};
		case DeepSkyCatalog::Variable15: return {"variable_stars_15.csv", "V003"};
	}
	return {};  // unreachable; keeps the compiler happy
}

/// Read every line of a text file, dropping a trailing CR so CRLF catalogs
/// parse identically to LF ones. Returns false when the file cannot be opened.
bool read_lines(const std::filesystem::path& path, std::vector<std::string>& out)
{
	std::ifstream in(path, std::ios::binary);
	if (!in.is_open()) return false;

	out.clear();
	std::string line;
	while (std::getline(in, line)) {
		if (!line.empty() && line.back() == '\r')
			line.pop_back();
		out.push_back(std::move(line));
	}
	return true;
}

}  // anonymous namespace

CatalogLoadStatus DeepSkyDatabase::load(DeepSkyCatalog catalog,
                                        const std::filesystem::path& base_path)
{
	// Cache hit: the requested catalog is already resident (Pascal database_nr).
	if (resident_ == catalog) return CatalogLoadStatus::AlreadyLoaded;

	const CatalogSpec spec = catalog_spec(catalog);

	std::vector<std::string> lines;
	if (!read_lines(base_path / spec.file_name, lines)) {
		// Pascal clears deepstring and abandons the load on any file error.
		lines_.clear();
		resident_.reset();
		return CatalogLoadStatus::NotFound;
	}

	lines_ = std::move(lines);
	resident_ = catalog;

	// Versioned catalogs stay loaded even when out of date; the caller decides
	// how to warn (Pascal shows a message box but keeps the data).
	if (!spec.required_version.empty()) {
		const std::string_view header = lines_.empty()
		                              ? std::string_view{}
		                              : std::string_view{lines_.front()};
		if (header.substr(0, spec.required_version.size()) != spec.required_version)
			return CatalogLoadStatus::Outdated;
	}

	return CatalogLoadStatus::Loaded;
}

// ---------------------------------------------------------------------------
// read_deepsky
// ---------------------------------------------------------------------------

std::vector<DeepSkyObject> read_deepsky(std::span<const std::string> database_lines,
                                        DeepSkySearch mode,
                                        double telescope_ra, double telescope_dec,
                                        double fov)
{
	std::vector<DeepSkyObject> results;

	const double cos_dec0 = std::cos(telescope_dec);
	const double fov2 = fov * fov;

	// Records start at index 2, skipping the two header lines (Pascal linepos:=2).
	for (size_t row = 2; row < database_lines.size(); ++row) {
		const std::string_view line = database_lines[row];
		size_t pos = 0;

		// Field 1: RA as integer, unit = 0.1 s of time (RA 00:00 00.1 = 1).
		const auto ra_field = next_field(line, pos);
		int ra_int = 0;
		{
			auto [ptr, ec] = std::from_chars(ra_field.data(),
			                                 ra_field.data() + ra_field.size(), ra_int);
			(void)ptr;
			if (ec != std::errc{}) continue;  // Pascal: fout<>0 -> skip line
		}

		// Field 2: Dec as integer, unit = 1 arcsec (DEC 00:00 01 = 1).
		const auto dec_field = next_field(line, pos);
		int dec_int = 0;
		{
			auto [ptr, ec] = std::from_chars(dec_field.data(),
			                                 dec_field.data() + dec_field.size(), dec_int);
			(void)ptr;
			if (ec != std::errc{}) continue;
		}

		const double ra  = ra_int  * std::numbers::pi * 2.0 / 864000.0;
		const double dec = dec_int * std::numbers::pi * 0.5 / 324000.0;

		// Field gate: skip records outside the FOV disc (Pascal 'S' mode).
		if (mode != DeepSkySearch::FullDatabase) {
			double delta_ra = std::abs(ra - telescope_ra);
			if (delta_ra > std::numbers::pi)
				delta_ra = std::numbers::pi * 2.0 - delta_ra;
			const double a = delta_ra * cos_dec0;
			const double b = dec - telescope_dec;
			if (a * a + b * b > fov2) continue;
		}

		// Field 3: name(s), slash-separated.
		const auto names = parse_names(next_field(line, pos));

		// Field 4: length (major axis), floating point.
		const double length = parse_double(next_field(line, pos));

		// Field 5: width (minor axis), floating point.
		const double width = parse_double(next_field(line, pos));

		// Field 6: position angle; PA=0 is valid, so an unparseable field
		// (including empty) maps to 999 = unknown, matching Pascal.
		const auto pa_field = next_field(line, pos);
		double pa = 999.0;
		{
			double v = 0.0;
			auto [ptr, ec] = astap::core::from_chars(pa_field.data(),
			                                         pa_field.data() + pa_field.size(), v);
			(void)ptr;
			if (ec == std::errc{}) pa = v;
		}

		DeepSkyObject obj;
		obj.name   = std::string(names.naam2);
		obj.name2  = std::string(names.naam3);
		obj.name3  = std::string(names.naam4);
		obj.ra     = ra;
		obj.dec    = dec;
		obj.length = length;
		obj.width  = width;
		obj.pa     = pa;
		results.push_back(std::move(obj));
	}

	return results;
}

// ---------------------------------------------------------------------------
// plot_deepsky
// ---------------------------------------------------------------------------

std::vector<PlottedDeepSky> plot_deepsky(const Header& head,
                                         std::span<const std::string> database_lines,
                                         bool variable_catalog,
                                         int base_font_size,
                                         int formalism,
                                         bool flip_horizontal,
                                         bool flip_vertical)
{
	std::vector<PlottedDeepSky> results;

	// Need an image and an astrometric solution (Pascal: naxis<>0 and cd1_1<>0).
	if (head.naxis == 0 || head.cd1_1 == 0.0) return results;

	// Field centre from the image centre — robust to an off-centre CRPIX.
	double telescope_ra = 0.0;
	double telescope_dec = 0.0;
	astap::core::pixel_to_celestial(head,
	                                (head.width + 1) / 2.0, (head.height + 1) / 2.0,
	                                formalism, telescope_ra, telescope_dec);

	// Field of view with 50% extra, as a radius in radians.
	const double half_w = 0.5 * head.width * head.cdelt1;
	const double half_h = 0.5 * head.height * head.cdelt2;
	const double fov = 1.5 * std::sqrt(half_w * half_w + half_h * half_h)
	                 * std::numbers::pi / 180.0;

	// Determinant sign flags a mirror-flipped image.
	const double flipped =
	    (head.cd1_1 * head.cd2_2 - head.cd1_2 * head.cd2_1 > 0.0) ? -1.0 : 1.0;

	const auto objects = read_deepsky(database_lines, DeepSkySearch::WithinField,
	                                  telescope_ra, telescope_dec, fov);

	for (const auto& obj : objects) {
		double fits_x = 0.0;
		double fits_y = 0.0;
		astap::core::celestial_to_pixel(head, obj.ra, obj.dec, fits_x, fits_y);

		// FITS pixels are 1-based; the image array is 0-based.
		int x = static_cast<int>(std::lround(fits_x - 1.0));
		int y = static_cast<int>(std::lround(fits_y - 1.0));

		// Cull to the image plus a 25% margin (pre-flip coordinates).
		if (!(x > -0.25 * head.width && x <= 1.25 * head.width &&
		      y > -0.25 * head.height && y <= 1.25 * head.height))
			continue;

		// Object radius in pixels; length is in 0.1-arcminute catalog units.
		const double len = obj.length / (std::abs(head.cdelt2) * 60.0 * 10.0 * 2.0);

		// Size gate: drop tiny objects on wide fields; variables always pass.
		if (!(head.cdelt2 < 0.25 / 60.0 || len >= 1.0 || variable_catalog))
			continue;

		double gx_orientation = (obj.pa + head.crota2) * flipped;

		if (flip_horizontal) {
			x = (head.width - 1) - x;
			gx_orientation = -gx_orientation;
		}
		if (flip_vertical)
			gx_orientation = -gx_orientation;
		else
			y = (head.height - 1) - y;

		PlottedDeepSky p;
		p.ra = obj.ra;
		p.dec = obj.dec;
		p.x = x;
		p.y = y;
		p.radius = len;
		p.orientation = gx_orientation;

		// Compose a designation for identification (always, when named).
		if (!obj.name.empty()) {
			if (obj.name2.empty())
				p.name = obj.name;
			else if (obj.name3.empty())
				p.name = obj.name + '/' + obj.name2;
			else
				p.name = obj.name + '/' + obj.name2 + '/' + obj.name3;
			p.abbr = obj.name;  // primary designation (Pascal naam2)
			p.is_reference = (obj.name.front() == '0');
		}

		// A label is drawn only when the centre is on-image and named.
		if (x >= 0 && x <= head.width - 1 && y >= 0 && y <= head.height - 1 &&
		    !obj.name.empty()) {
			p.has_label = true;
			p.font_size = static_cast<int>(std::lround(std::min(
			    20.0, std::max(static_cast<double>(base_font_size), len / 2.0))));
		}

		// Marker shape: a zero minor axis collapses to a circle.
		double length1 = obj.length;
		double width1 = obj.width;
		double pa = obj.pa;
		if (width1 == 0.0) {
			width1 = length1;
			pa = 999.0;
		}

		p.axis_ratio = (length1 != 0.0) ? (width1 / length1) : 1.0;
		p.pen_width = static_cast<int>(
		    std::min(4.0, std::max(1.0, std::round(len / 70.0))));

		if (len <= 2.0)
			p.marker = DeepSkyMarker::Dots;
		else if (pa != 999.0)
			p.marker = DeepSkyMarker::Galaxy;
		else
			p.marker = DeepSkyMarker::Circle;

		results.push_back(std::move(p));
	}

	return results;
}

// ---------------------------------------------------------------------------
// extract_visible
// ---------------------------------------------------------------------------

std::vector<VariableInfo> extract_visible(std::span<const PlottedDeepSky> plotted,
                                          int source, std::size_t max_count)
{
	std::vector<VariableInfo> targets;

	for (const PlottedDeepSky& p : plotted) {
		// Only on-image, named objects become photometry targets (Pascal draws
		// the label and records the target under the same gate).
		if (!p.has_label) continue;
		if (targets.size() >= max_count) break;

		VariableInfo v;
		v.ra = p.ra;
		v.dec = p.dec;
		v.abbr = p.abbr;
		v.source = source;
		v.index = static_cast<int>(targets.size());
		targets.push_back(std::move(v));
	}

	return targets;
}

// ---------------------------------------------------------------------------
// layout_labels
// ---------------------------------------------------------------------------

std::vector<LabelBox> layout_labels(std::span<const LabelInput> labels,
                                    int image_width, int image_height)
{
	std::vector<LabelBox> boxes;
	boxes.reserve(labels.size());

	for (const auto& lab : labels) {
		const int x = lab.x;
		const int y = lab.y;
		const int th = lab.th;

		int x1 = x;
		int y1 = y;
		int x2 = x + lab.tw;
		int y2 = y + th;

		// Shift left when the text runs off the right edge.
		if (x1 <= image_width && x2 > image_width) {
			x1 -= (x2 - image_width);
			x2 = image_width;
		}

		if (!boxes.empty()) {
			bool overlap = false;
			do {
				overlap = false;
				for (const LabelBox& b : boxes) {
					const bool hit =
					    (x1 >= b.x1 && x1 <= b.x2 && y1 >= b.y1 && y1 <= b.y2) ||
					    (x2 >= b.x1 && x2 <= b.x2 && y1 >= b.y1 && y1 <= b.y2) ||
					    (x1 >= b.x1 && x1 <= b.x2 && y2 >= b.y1 && y2 <= b.y2) ||
					    (x2 >= b.x1 && x2 <= b.x2 && y2 >= b.y1 && y2 <= b.y2) ||
					    (b.x1 >= x1 && b.x1 <= x2 && b.y1 >= y1 && b.y1 <= y2) ||
					    (b.x2 >= x1 && b.x2 <= x2 && b.y2 >= y1 && b.y2 <= y2);
					if (hit) {
						overlap = true;
						// Try shifting one third of a line down.
						y1 += th / 3;
						y2 += th / 3;
						if (y2 >= image_height) {
							// No vertical space: fall back to the anchor.
							y1 = y;
							y2 = y + th;
							overlap = false;
						}
						break;  // Restart the scan from the first box.
					}
				}
			} while (overlap);
		}

		LabelBox box;
		box.x1 = x1;
		box.y1 = y1;
		box.x2 = x2;
		box.y2 = y2;
		box.connector = (y1 != y);
		boxes.push_back(box);
	}

	return boxes;
}

}  // namespace astap::analysis
