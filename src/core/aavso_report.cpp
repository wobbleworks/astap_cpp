///----------------------------------------
///      @file aavso_report.cpp
///   @ingroup ASTAP++
///     @brief Implementation of @ref astap::core::format_aavso_report.
///    @author Created by John Stephen on 4/24/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "aavso_report.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <sstream>
#include <string>

///----------------------------------------
namespace astap::core {
///----------------------------------------

namespace {

// AAVSO line endings are CR/LF per the spec — they go through Windows tooling
// most of the time and the spec is explicit about it.
constexpr auto kCrlf = "\r\n";

[[nodiscard]] std::string fixed3(double v) {
	return std::format("{:.3f}", v);
}

[[nodiscard]] std::string fixed4(double v) {
	return std::format("{:.4f}", v);
}

[[nodiscard]] std::string fixed5(double v) {
	return std::format("{:.5f}", v);
}

[[nodiscard]] std::string str_or_na(double v, double na_threshold = 98.0) {
	if (v >= na_threshold) return "na";
	return std::format("{:.3f}", v);
}

// Build the AAVSO header block. Returns a string ending in "#\r\n#NAME...\r\n".
[[nodiscard]] std::string build_header(const AavsoOptions& opts) {
	const auto type_field = opts.baa_style ? "AAVSO EXT BAA V1.00" : "Extended";

	auto out = std::string{};
	out += "#TYPE=";        out += type_field;        out += kCrlf;
	out += "#OBSCODE=";     out += opts.observer_code; out += kCrlf;
	out += "#SOFTWARE=ASTAP++";
	if (!opts.software_version.empty()) {
		out += ", v";
		out += opts.software_version;
	}
	if (!opts.software_settings.empty()) {
		out += " (";
		out += opts.software_settings;
		out += ")";
	}
	out += kCrlf;

	// Convert tab to the literal word "tab" — that's what AAVSO ingests.
	const auto delim_label = (opts.delimiter == "\t") ? std::string("tab")
	                                                  : opts.delimiter;
	out += "#DELIM=";       out += delim_label;        out += kCrlf;
	out += "#DATE=";        out += (opts.hjd_date ? "HJD" : "JD"); out += kCrlf;
	out += "#OBSTYPE=CCD";  out += kCrlf;
	out += "#COMMENTS=";    out += opts.comments;      out += kCrlf;

	if (opts.baa_style) {
		out += "#LOCATION=";
		out += opts.site_lat;
		out += ' ';
		out += opts.site_long;
		out += ' ';
		out += opts.site_elev;
		out += kCrlf;
		out += "#TELESCOPE="; out += opts.telescope; out += kCrlf;
		out += "#CAMERA=";    out += opts.camera;    out += kCrlf;
	}

	out += "#";  out += kCrlf;

	// Column header line — uses the actual delimiter.
	const auto& d = opts.delimiter;
	out += "#NAME"; out += d; out += "DATE";    out += d; out += "MAG";
	out += d;      out += "MERR";   out += d;   out += "FILT";
	out += d;      out += "TRANS";  out += d;   out += "MTYPE";
	out += d;      out += "CNAME";  out += d;   out += "CMAG";
	out += d;      out += "KNAME";  out += d;   out += "KMAG";
	out += d;      out += "AIRMASS"; out += d;  out += "GROUP";
	out += d;      out += "CHART";  out += d;   out += "NOTES";
	out += kCrlf;

	return out;
}

[[nodiscard]] std::string build_row(const AavsoMeasurement& m,
                                    const AavsoOptions& opts) {
	const auto& d = opts.delimiter;

	// Magnitudes potentially get the comp-vs-catalog correction applied
	// when not in ensemble mode.
	auto var_mag   = m.var_magnitude;
	auto check_mag = m.check_magnitude;
	auto cname     = std::string{"ENSEMBLE"};
	auto cmag      = std::string{"na"};

	if (!opts.ensemble && m.comp_catalog_mag != 0.0 && m.comp_magnitude != 0.0) {
		// Pascal: magn_correction = documented - measured; var_mag += correction;
		// check_mag += correction. No flux conversion needed — magnitude deltas
		// hold for any ratio.
		const auto correction = m.comp_catalog_mag - m.comp_magnitude;
		var_mag   += correction;
		check_mag += correction;
		cname     = m.comp_name.empty() ? std::string{"COMP"} : m.comp_name;
		cmag      = fixed3(m.comp_catalog_mag);
	}

	// Apply the slope/colour-transform correction last (per Pascal:
	// var_magn := var_magn + delta_bv * magnitude_slope).
	var_mag += opts.delta_bv * opts.magnitude_slope;

	// Variable error: prefer the explicitly-supplied value; otherwise derive
	// 2/SNR (a conservative practical default — matches Pascal's err_by_snr).
	auto var_err = m.var_error;
	if (var_err <= 0.0 && m.snr > 0.0) {
		var_err = 2.0 / m.snr;
	}

	auto row = std::string{};
	row += m.variable_name;                            row += d;
	row += fixed5(m.jd);                               row += d;
	row += fixed3(var_mag);                            row += d;
	row += fixed4(var_err);                            row += d;
	row += m.filter_band;                              row += d;
	row += "NO";                                       row += d;  // TRANS
	row += "STD";                                      row += d;  // MTYPE
	row += cname;                                      row += d;
	row += cmag;                                       row += d;
	row += m.check_name;                               row += d;
	row += fixed3(check_mag);                          row += d;
	row += str_or_na(m.airmass);                       row += d;
	row += opts.group_id;                              row += d;
	row += opts.chart_id;                              row += d;
	row += opts.notes.empty() ? std::string{"na"} : opts.notes;
	row += kCrlf;
	return row;
}

}  // namespace

std::string clean_abbreviation(std::string s) {
	// Drop "<text>, σ=…" tail.
	if (auto pos = s.find(','); pos != std::string::npos) {
		s.erase(pos);
	}
	// Trim leading/trailing whitespace.
	const auto first = s.find_first_not_of(" \t\r\n");
	if (first == std::string::npos) {
		return {};
	}
	const auto last = s.find_last_not_of(" \t\r\n");
	s = s.substr(first, last - first + 1);
	// Underscores → spaces (matches Pascal stringreplace).
	std::ranges::replace(s, '_', ' ');
	return s;
}

std::string format_aavso_report(const AavsoMeasurement& m,
                                const AavsoOptions& opts) {
	auto out = build_header(opts);
	out += build_row(m, opts);
	return out;
}

std::string format_aavso_report(const std::vector<AavsoMeasurement>& measurements,
                                const AavsoOptions& opts) {
	auto out = build_header(opts);
	for (const auto& m : measurements) {
		out += build_row(m, opts);
	}
	return out;
}

} // namespace
