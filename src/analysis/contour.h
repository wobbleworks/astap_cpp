#pragma once

/// @file contour.h
/// Streak detection helpers — ported from unit_contour.pas.
///
/// Contains the pure-math functions for satellite-streak line fitting:
/// perpendicular point-to-line distance, linear least-squares trendline,
/// and an outlier-rejecting variant. The GUI-dependent `contour()` procedure
/// (Moore Neighbor Contour Tracing with canvas drawing) is intentionally
/// omitted; see the TODO below.
///
/// Copyright (C) 2023, 2024 by Han Kleijn, www.hnsky.org (original Pascal).
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this file,
/// You can obtain one at https://mozilla.org/MPL/2.0/.

#include "../types.h"

namespace astap::analysis {

using astap::StarList;

/// A detected satellite streak described by the line Y = slope * X + intercept.
struct Streak {
	double slope{};
	double intercept{};
};

/// Result of a trendline fit with outlier rejection.
struct TrendlineResult {
	double slope{};
	double intercept{};
	double sd{};          ///< Standard deviation of perpendicular distances.
};

/// Returns the perpendicular distance from point (@p fits_x, @p fits_y)
/// to the line Y = @p slope * X + @p intercept.
///
/// @see https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
[[nodiscard]] double line_distance(double fits_x, double fits_y,
                                   double slope, double intercept);

/// Ordinary least-squares linear fit: Y = slope * X + intercept.
///
/// @param xylist  StarList where [0] holds X values and [1] holds Y values.
/// @param len     Number of data points to use (may be less than xylist[0].size()).
/// @param[out] slope     Fitted slope.
/// @param[out] intercept Fitted intercept.
void trendline(const StarList& xylist, int len,
               double& slope, double& intercept);

/// Two-pass trendline: fits once, computes the standard deviation of
/// perpendicular distances, removes outliers beyond 1.5 * SD, then re-fits.
///
/// @param xylist  StarList where [0] holds X values and [1] holds Y values.
/// @param len     Number of data points to use.
/// @param[out] slope     Fitted slope (after outlier removal).
/// @param[out] intercept Fitted intercept (after outlier removal).
/// @param[out] sd        Standard deviation of perpendicular distances from
///                        the first (pre-rejection) fit.
void trendline_without_outliers(const StarList& xylist, int len,
                                double& slope, double& intercept, double& sd);

// TODO(astap-port): The `contour()` procedure from unit_contour.pas implements
// Moore Neighbor Contour Tracing with Gaussian blur, background estimation,
// and streak detection. It is deeply coupled to the Lazarus GUI
// (mainwindow.image1.Canvas drawing, application.processmessages, etc.) and
// needs a canvas-free refactor before it can be ported here.

}  // namespace astap::analysis
