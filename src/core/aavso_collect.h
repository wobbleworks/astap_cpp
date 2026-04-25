///----------------------------------------
///      @file aavso_collect.h
///   @ingroup ASTAP++
///     @brief Multi-frame photometry collector for AAVSO time-series reports.
///   @details Phase 5b of the @c unit_aavso.pas port. Given a list of
///            previously plate-solved + photometrically-calibrated FITS
///            files and the celestial coordinates of three reference
///            stars (variable, check, comp), this module loads each
///            frame in turn, projects the RA/Dec of each star into that
///            frame's pixel space (so the same physical star is measured
///            even if pointing drifted), runs aperture photometry, and
///            packs the results into a vector of @ref AavsoMeasurement
///            rows ready for @ref format_aavso_report.
///    @author Created by John Stephen on 4/25/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include "aavso_report.h"

#include <filesystem>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

///----------------------------------------
namespace astap::core {
///----------------------------------------

///----------------------------------------
/// @struct PhotometryStarSpec
/// @brief A reference star fixed by celestial coordinates so it can be
///        re-located in every frame regardless of pointing drift.
///----------------------------------------

struct PhotometryStarSpec {
	std::string name;     ///< @brief AAVSO designation or label.
	double      ra  = 0;  ///< @brief Right ascension, radians.
	double      dec = 0;  ///< @brief Declination, radians.
};

///----------------------------------------
/// @struct AavsoCollectOptions
/// @brief Inputs governing how the collector measures each frame.
///----------------------------------------

struct AavsoCollectOptions {
	PhotometryStarSpec variable;
	PhotometryStarSpec check;
	PhotometryStarSpec comp;             ///< @brief Ignored when @c ensemble is true.
	bool   ensemble                = true;
	double comp_catalog_magnitude  = 0.0; ///< @brief Required when @c ensemble is false.
	bool   hjd_date                = false;
	std::string filter_band        = "V";
	int    annulus_radius          = 14;  ///< @brief HFD search/aperture radius in pixels.

	/// @brief Optional progress callback: per-file message, never nullable.
	std::function<void(const std::string&)> log;
};

///----------------------------------------
/// @struct AavsoCollectResult
/// @brief Row-by-row outcome of @ref collect_aavso_measurements.
///----------------------------------------

struct AavsoCollectResult {
	std::vector<AavsoMeasurement> rows;          ///< @brief One per accepted frame.
	std::vector<std::string>      skipped_files; ///< @brief Files that failed to load / measure.
	int frames_total    = 0;
	int frames_accepted = 0;
};

///----------------------------------------
/// @brief Walk @p files, measure the three stars in each, and return
///        ready-to-report rows.
/// @details Loads each file via @c astap::core::load_fits. Skips files
///          without a plate solution (@c head.cd1_1 == 0) or without
///          @c MZERO. Projects each @ref PhotometryStarSpec to pixel
///          coords with @c celestial_to_pixel, runs @c HFD at that
///          location to get a refined centroid + flux + SNR, computes
///          @c mag = head.mzero − 2.5 log10(flux). Builds the per-row
///          @c AavsoMeasurement with the JD (or HJD if requested) from
///          the frame's own @c DATE-OBS / @c DATE-AVG.
/// @param files Paths to FITS frames in observation order.
/// @param opts  Star specs + options.
/// @return Per-row results plus a list of skipped files for diagnostics.
///----------------------------------------

[[nodiscard]] AavsoCollectResult collect_aavso_measurements(
    const std::vector<std::filesystem::path>& files,
    const AavsoCollectOptions& opts);

///----------------------------------------
/// @brief Find the documented AAVSO VSP magnitude of the comparison star
///        nearest @p ra / @p dec in the given filter band.
/// @details Walks the @c astap::core::vsp global cache (populated by
///          @c download_vsp) and returns the catalog magnitude of the
///          nearest entry whose distance is below @p max_arcsec. The
///          filter symbol selects which colour the AAVSO chart returned;
///          recognised values are @c "V", @c "B", @c "R", @c "SG",
///          @c "SR", @c "SI". Returns @c std::nullopt when the cache is
///          empty, no entry is within range, or the requested band is
///          unavailable for the matched entry.
/// @param filter_band The AAVSO filter symbol.
/// @param ra          Right ascension of the picked comp star (radians).
/// @param dec         Declination of the picked comp star (radians).
/// @param max_arcsec  Maximum acceptable separation (default 5″).
///----------------------------------------

[[nodiscard]] std::optional<double> find_vsp_magnitude(
    std::string_view filter_band, double ra, double dec,
    double max_arcsec = 5.0);

} // namespace
