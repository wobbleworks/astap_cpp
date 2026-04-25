///----------------------------------------
///      @file annotation_scanner.h
///   @ingroup ASTAP++
///     @brief Finds deep-sky objects within a solved image's field of view.
///   @details Loads the ASTAP @c deep_sky.csv catalog and projects matching
///            objects into FITS pixel coordinates using the image's WCS
///            solution. Produces overlay markers that @ref ImageViewer can
///            render.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include "../../src/analysis/asteroid_overlay.h"
#include "../../src/core/online.h"
#include "../../src/types.h"

#include <QString>
#include <QPointF>

#include <filesystem>
#include <string>
#include <vector>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @struct AnnotationMarker
/// @brief A deep-sky object projected onto image pixel coordinates.
///----------------------------------------

struct AnnotationMarker {
	QString name;
	double x = 0.0;        // FITS pixel, 1-based
	double y = 0.0;
	double majorPx = 0.0;  // ellipse semi-major axis in pixels
	double minorPx = 0.0;  // ellipse semi-minor axis in pixels
	double paDeg = 0.0;    // position angle in degrees (N through E)
};

///----------------------------------------
/// @brief Load @c deep_sky.csv from the given directory into memory.
/// @param dir Folder containing the catalog file (typically database_path).
/// @return True if the catalog was loaded (at least a few hundred lines).
///----------------------------------------

[[nodiscard]] bool load_deepsky_catalog(const std::filesystem::path& dir);

///----------------------------------------
/// @brief Scan the loaded catalog for objects inside the current image's
///        solved field of view and return overlay markers.
/// @param head FITS header with a valid WCS solution.
/// @return Annotation markers for objects within the FOV.
///----------------------------------------

[[nodiscard]] std::vector<AnnotationMarker> annotate_image(const astap::Header& head);

///----------------------------------------
/// @struct ConstellationLine
/// @brief A single line segment of a constellation stick figure.
///----------------------------------------

struct ConstellationLine {
	double x1 = 0.0;  // FITS pixel, 1-based
	double y1 = 0.0;
	double x2 = 0.0;
	double y2 = 0.0;
};

///----------------------------------------
/// @struct ConstellationLabel
/// @brief A constellation name positioned in the image.
///----------------------------------------

struct ConstellationLabel {
	QString name;
	double x = 0.0;
	double y = 0.0;
};

///----------------------------------------
/// @struct ConstellationOverlay
/// @brief All constellation stick-figure lines and labels within the FOV.
///----------------------------------------

struct ConstellationOverlay {
	std::vector<ConstellationLine> lines;
	std::vector<ConstellationLabel> labels;
};

///----------------------------------------
/// @brief Project constellation stick figures and labels into pixel
///        coordinates for the current solved image.
/// @param head FITS header with a valid WCS solution.
/// @return Constellation lines and labels within the FOV.
///----------------------------------------

[[nodiscard]] ConstellationOverlay build_constellation_overlay(const astap::Header& head);

///----------------------------------------
/// @struct CatalogStarMarker
/// @brief A catalog star projected into pixel coordinates, with magnitude
///        and colour information suitable for overlay rendering.
///----------------------------------------

struct CatalogStarMarker {
	double x = 0.0;          ///< FITS pixel, 1-based
	double y = 0.0;
	double magn = 0.0;       ///< Gaia magnitude (e.g. 7.5)
	double bpRp = 999.0;     ///< Gaia Bp-Rp colour (sentinel 999 = unknown/mono catalog)
};

///----------------------------------------
/// @brief Scan the active star catalog for entries inside the solved image
///        and return overlay markers.
/// @details Uses whichever catalog `astap::reference::select_star_database`
///          has chosen (local .1476/.290 tiles, wide-field w08, or online Gaia).
///          Caller must ensure a catalog is loaded (e.g. via the plate
///          solver or photometric calibration). Results are sorted
///          brightest-first and capped at @p max_markers entries.
/// @param head FITS header with a valid WCS solution.
/// @param max_markers Upper bound on returned markers (typical: 1000).
/// @return Overlay markers for catalog stars inside the field of view.
///----------------------------------------

[[nodiscard]] std::vector<CatalogStarMarker> scan_catalog_stars(
    const astap::Header& head, int max_markers = 1000);

///----------------------------------------
/// @struct VarStarMarker
/// @brief Variable / comparison star projected into pixel coordinates.
///----------------------------------------

struct VarStarMarker {
	QString name;
	QString magText;     ///< @brief V mag (or B fallback) as text; empty = unknown
	double x = 0.0;      ///< @brief FITS pixel, 1-based
	double y = 0.0;
	bool isComparison = false;  ///< @brief true = AAVSO VSP comp/check; false = VSX variable
};

///----------------------------------------
/// @brief Project the cached AAVSO @c vsx (variables) and @c vsp
///        (comparison/check) catalogs into pixel space.
/// @details Caller must have populated the caches via
///          @c astap::core::download_vsx and @c astap::core::download_vsp
///          (or @c variable_star_annotation) first.
/// @param head FITS header with a valid WCS solution.
/// @return Variable + comparison-star markers inside the FOV.
///----------------------------------------

[[nodiscard]] std::vector<VarStarMarker> project_var_stars(
    const astap::Header& head);

///----------------------------------------
/// @struct SimbadMarker
/// @brief Simbad-resolved object projected into pixel coordinates.
///----------------------------------------

struct SimbadMarker {
	QString name;
	QString type;        ///< @brief Simbad maintype (e.g. "Galaxy", "OpC", "*")
	double x = 0.0;
	double y = 0.0;
	double magnitude = 0.0;  ///< @brief 0 = unknown
	double sizeArcmin = 0.0; ///< @brief 0 = unknown (single-object reply only)
};

///----------------------------------------
/// @brief Project parsed Simbad objects into pixel space.
/// @param objects Result of @c astap::core::plot_simbad on a fetched body.
/// @param head    FITS header with a valid WCS solution.
/// @return Markers inside the FOV.
///----------------------------------------

[[nodiscard]] std::vector<SimbadMarker> project_simbad(
    const std::vector<astap::core::SimbadObject>& objects,
    const astap::Header& head);

///----------------------------------------
/// @struct VizierMarker
/// @brief Vizier (Gaia) star projected into pixel coordinates.
///----------------------------------------

struct VizierMarker {
	double x = 0.0;
	double y = 0.0;
	double magnitude = 0.0;
};

///----------------------------------------
/// @brief Project parsed Vizier rows into pixel space.
/// @param rows  Result of @c astap::core::plot_vizier on a fetched body.
/// @param head  FITS header with a valid WCS solution.
/// @return Markers inside the FOV, brightest first.
///----------------------------------------

[[nodiscard]] std::vector<VizierMarker> project_vizier(
    const std::vector<astap::core::VizierObject>& rows,
    const astap::Header& head);

///----------------------------------------
/// @struct AsteroidMarker
/// @brief A minor planet or comet projected into pixel coordinates.
///----------------------------------------

struct AsteroidMarker {
	QString label;      ///< Designation + (optionally) magnitude.
	double x = 0.0;     ///< FITS pixel, 1-based
	double y = 0.0;
	bool isComet = false;
	bool outdated = false;
};

///----------------------------------------
/// @brief Project asteroid/comet detections into pixel space.
/// @param detections Result of @ref astap::analysis::scan_asteroids_in_field.
/// @param head FITS header with a valid WCS solution.
/// @return Markers inside the image bounds.
///----------------------------------------

[[nodiscard]] std::vector<AsteroidMarker> project_asteroids(
    const std::vector<astap::analysis::AsteroidDetection>& detections,
    const astap::Header& head);

} // namespace astap::gui
