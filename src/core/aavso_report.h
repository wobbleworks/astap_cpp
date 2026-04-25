///----------------------------------------
///      @file aavso_report.h
///   @ingroup ASTAP++
///     @brief AAVSO Extended File Format report generator.
///   @details Ports the report-formatting half of Han Kleijn's
///            @c unit_aavso.pas. Contains:
///             - data structures for a per-frame photometry measurement,
///             - configuration knobs (observer code, filter, delimiter,
///               JD-vs-HJD, BAA-style header, ensemble vs. comp star),
///             - @ref format_aavso_report which produces the AAVSO submission
///               text from one or more measurements,
///             - @ref clean_abbreviation, which strips the
///               @c ", σ=0.012" tail the GUI appends to comp/check labels.
///
///            Multi-frame photometry data collection (Pascal's @c listview7)
///            is a separate concern; the dialog feeds this module @em already
///            measured rows.
///    @author Created by John Stephen on 4/24/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <string>
#include <string_view>
#include <vector>

///----------------------------------------
namespace astap::core {
///----------------------------------------

///----------------------------------------
/// @struct AavsoMeasurement
/// @brief One photometry observation. A single-frame report uses one of
///        these; a time-series report uses a vector.
/// @details Magnitudes here are the @em instrumental ones (head.mzero already
///          applied). When @ref AavsoOptions::ensemble is @c false and the
///          comp star's catalog magnitude is provided, the report formatter
///          applies the ensemble-vs-comp delta correction itself.
///----------------------------------------

struct AavsoMeasurement {
	std::string variable_name;        ///< @brief AAVSO variable designation (e.g. "SS Cyg").
	std::string check_name;           ///< @brief AAVSO check-star label (e.g. "000-BCP-306").
	std::string comp_name;            ///< @brief Comp-star label, or "ENSEMBLE".
	double      var_magnitude    = 0; ///< @brief Instrumental V mag of the variable.
	double      var_error        = 0; ///< @brief 1σ error in mag (0 = "use SNR-derived default").
	double      check_magnitude  = 0; ///< @brief Instrumental V mag of the check star.
	double      comp_magnitude   = 0; ///< @brief Instrumental V mag of the comp star (ignored in ensemble mode).
	double      comp_catalog_mag = 0; ///< @brief Documented (catalog) V mag of the comp star (ignored in ensemble mode).
	double      jd               = 0; ///< @brief Julian Date for the @c DATE column (caller chooses JD or HJD).
	double      airmass          = 99; ///< @brief Airmass; 99 = "na" in the report.
	double      snr              = 0; ///< @brief Variable's SNR. Used to derive @c var_error when the latter is 0.
	std::string filter_band      = "V"; ///< @brief Filter symbol the report should use ("V", "B", "R", "TR", etc.).
};

///----------------------------------------
/// @struct AavsoOptions
/// @brief Per-report configuration knobs.
///----------------------------------------

struct AavsoOptions {
	std::string observer_code;             ///< @brief AAVSO observer code (e.g. "JST01").
	std::string delimiter      = ",";      ///< @brief Field delimiter; pass @c "\t" for tab.
	bool        hjd_date       = false;    ///< @brief @c true → @c #DATE=HJD; @c false → @c #DATE=JD.
	bool        baa_style      = false;    ///< @brief @c true → @c #TYPE=AAVSO EXT BAA V1.00 plus @c #LOCATION/#TELESCOPE/#CAMERA.
	bool        ensemble       = false;    ///< @brief @c true → CNAME = "ENSEMBLE", CMAG = "na".
	double      delta_bv       = 0;        ///< @brief (B−V)var − (B−V)comp for the slope correction.
	double      magnitude_slope = 0;       ///< @brief 2nd-order transformation slope.

	// #SOFTWARE field components:
	std::string software_version;          ///< @brief e.g. "2026.04.24".
	std::string software_settings;         ///< @brief Free-form (catalog name, aperture, annulus).

	// BAA-style header extras (only emitted when @c baa_style is @c true):
	std::string site_lat;
	std::string site_long;
	std::string site_elev;
	std::string telescope;
	std::string camera;

	// Per-row defaults (rarely changed in practice):
	std::string chart_id  = "na";          ///< @brief KCHART column.
	std::string group_id  = "na";          ///< @brief GROUP column.
	std::string comments;                  ///< @brief @c #COMMENTS header value.
	std::string notes;                     ///< @brief NOTES column (free-form per-row, but uniform here).
};

///----------------------------------------
/// @brief Strip the GUI's @c ", σ=…" tail and replace underscores with spaces.
/// @details The Pascal port populates comp/check combos with strings like
///          @c "000-BCP-306, σ=0.012" or @c "SS_Cyg, σ=0.005". Both halves
///          need cleaning before the abbreviation lands in the AAVSO report.
/// @param s Combo-box text.
/// @return Cleaned designation (trailing comma-onwards trimmed, underscores → spaces).
///----------------------------------------

[[nodiscard]] std::string clean_abbreviation(std::string s);

///----------------------------------------
/// @brief Format a single-row AAVSO Extended File Format report.
/// @param m   The measurement.
/// @param opts Formatting options.
/// @return Newline-terminated report string ready for clipboard / file.
///----------------------------------------

[[nodiscard]] std::string format_aavso_report(const AavsoMeasurement& m,
                                              const AavsoOptions& opts);

///----------------------------------------
/// @brief Format a multi-row (time-series) AAVSO report.
/// @param measurements Per-frame measurements; emitted in input order.
/// @param opts         Formatting options applied to all rows.
/// @return Newline-terminated report string ready for clipboard / file.
///----------------------------------------

[[nodiscard]] std::string format_aavso_report(
    const std::vector<AavsoMeasurement>& measurements,
    const AavsoOptions& opts);

} // namespace
