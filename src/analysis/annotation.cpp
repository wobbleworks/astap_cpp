/// @file annotation.cpp
/// Image annotation utilities — ported from unit_annotation.pas.

#include "annotation.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <optional>
#include <span>
#include <string>
#include <string_view>

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
	std::from_chars(s.data(), s.data() + s.size(), val);
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
			auto [ptr, ec] = std::from_chars(pa_field.data(),
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

}  // namespace astap::analysis
