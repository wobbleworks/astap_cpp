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

#include <vector>

namespace astap::analysis {

using astap::StarList;
using astap::ImageArray;

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

/// Core satellite-streak detector — the canvas-free heart of unit_contour.pas
/// contour(). Runs Moore-Neighbor contour tracing over @p img_bk (a mono,
/// already-blurred image), fills each traced blob's interior to measure its
/// area, and reports blobs long and thin enough to be a streak: surface
/// > 400 px, length > 200 px, length^2 / surface > 10, and a trendline fit
/// with sigma < 10. Depends on nothing outside this file, so it is directly
/// unit-testable.
///
/// @param img_bk          Mono, already-blurred image (channel 0 is used).
/// @param detection_level Pixel threshold (sigmafactor * noise + background).
/// @param binning         1 for mono/colour, 2 for a binned Bayer frame;
///                        scales the reported intercept/sigma back to full res.
/// @param detection_grid  Only start a trace from pixels on this coordinate
///                        grid (<= 0 traces from every pixel).
/// @return One @ref Streak per detected line, in full-resolution coordinates.
[[nodiscard]] std::vector<Streak> detect_streaks(const ImageArray& img_bk,
                                                 double detection_level,
                                                 int binning,
                                                 int detection_grid);

/// Full streak detection matching unit_contour.pas contour(): converts a
/// colour/Bayer @p img to mono, applies a Gaussian blur of radius @p blur,
/// estimates the background/noise, then runs @ref detect_streaks with a
/// threshold of @p sigmafactor * noise + background. Unlike the original — which
/// blurs the shared image buffer in place — this leaves @p img untouched.
///
/// @param detection_grid Grid spacing (full-res); divided by the binning used.
[[nodiscard]] std::vector<Streak> contour(const ImageArray& img,
                                          astap::Header& head,
                                          double blur, double sigmafactor,
                                          int detection_grid = 400);

}  // namespace astap::analysis
