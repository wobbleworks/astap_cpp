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

} // namespace astap::gui
