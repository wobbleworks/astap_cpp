///----------------------------------------
///      @file asteroid_overlay.h
///   @ingroup ASTAP++
///     @brief Scan an MPCORB.DAT (and optional comet-elements) catalog for
///            objects whose ephemerides fall inside a solved image's field
///            of view.
///   @details Ported from the orchestrator half of Han Kleijn's
///            @c unit_asteroid.pas::plot_mpcorb. The per-object propagation,
///            parallax, and magnitude correction primitives already live
///            in @ref asteroid.h — this module glues them together and
///            yields (RA, Dec, mag, name) tuples the GUI can project.
///    @author Created by John Stephen on 4/23/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include "../types.h"

#include <filesystem>
#include <functional>
#include <string>
#include <vector>

///----------------------------------------
namespace astap::analysis {
///----------------------------------------

///----------------------------------------
/// @struct AsteroidDetection
/// @brief One minor planet or comet resolved inside the current frame.
///----------------------------------------

struct AsteroidDetection {
	std::string desn;      ///< Packed designation (short).
	std::string name;      ///< Full name (may be empty for unnumbered objects).
	double ra        = 0.0;    ///< Apparent right ascension, radians.
	double dec       = 0.0;    ///< Apparent declination, radians.
	double magnitude = 99.0;   ///< Apparent magnitude at @p jd_mid.
	bool   is_comet  = false;
	bool   outdated  = false;  ///< Elements are >120 days from observation epoch.
};

///----------------------------------------
/// @struct AsteroidScanOptions
/// @brief Inputs for @ref scan_asteroids_in_field.
///----------------------------------------

struct AsteroidScanOptions {
	/// @brief Path to an MPCORB.DAT (or compatible) asteroid catalog.
	///        Empty = skip asteroids.
	std::filesystem::path mpcorb_path;

	/// @brief Path to an MPC comet-elements file. Empty = skip comets.
	std::filesystem::path comets_path;

	/// @brief Upper bound on lines scanned per catalog (caps runtime on
	///        the full MPCORB.DAT, which has ~1.3M entries).
	int max_count = 1'500'000;

	/// @brief Drop objects fainter than this magnitude.
	double max_magnitude = 18.0;

	/// @brief Optional progress callback (called every ~10k lines).
	std::function<void(const std::string&)> log;
};

///----------------------------------------
/// @brief Scan an asteroid (and optional comet) catalog for objects whose
///        apparent position at the frame's @c DATE-AVG falls inside the
///        image's field of view.
/// @details Reads @c jd_mid / @c site_lat_radians / @c site_long_radians
///          from the engine globals (populated by
///          @c astap::stacking::date_to_jd and the FITS loader). The
///          caller must therefore make sure @p head has been loaded and
///          has a valid WCS solution and @c DATE-OBS.
/// @param head Solved FITS header (requires @c cdelt, @c ra0/dec0, etc.).
/// @param opts Catalog paths + magnitude / count limits + logger.
/// @return Detections whose apparent position is inside the FOV and
///         whose magnitude is at or below @p opts.max_magnitude.
///----------------------------------------

[[nodiscard]] std::vector<AsteroidDetection> scan_asteroids_in_field(
    const Header& head, const AsteroidScanOptions& opts);

} // namespace
