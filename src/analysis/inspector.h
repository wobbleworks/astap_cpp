#pragma once

/// @file inspector.h
/// CCD inspector helpers — ported from unit_inspector_plot.pas.
///
/// Pure algorithmic functions for star-field analysis: star aspect-ratio
/// (elongation) measurement, HFD filtering via nearest-neighbour medians,
/// and angular modulo wrapping.  GUI / canvas code from the original unit
/// is intentionally omitted.
///
/// Copyright (C) 2018, 2022 by Han Kleijn, www.hnsky.org (original Pascal).
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this file,
/// You can obtain one at https://mozilla.org/MPL/2.0/.

#include "../types.h"

namespace astap::analysis {

using astap::ImageArray;
using astap::StarList;

/// Result of a single-star elongation measurement.
struct StarAspect {
	double aspect{999.0};     ///< Elongation ratio (max/min spread); 999 = failure.
	int    orientation{0};    ///< Angle of minimum spread in [0, 179] degrees.
};

/// Measure the elongation (aspect ratio) and orientation of a single star.
///
/// For each pixel within radius @p rs of the centroid (@p x1, @p y1) that is
/// brighter than 7 * @p sd_bg above @p star_bg, the function builds a
/// weighted-distance map and then sweeps 0-179 degrees to find the angle that
/// minimises the sum of |distance-to-line|.  The aspect ratio is max/min of
/// those sums.
///
/// @param img      Image buffer (only channel 0 is used).
/// @param x1       Star centroid X (subpixel).
/// @param y1       Star centroid Y (subpixel).
/// @param rs       Search radius in pixels (clamped to 51).
/// @param star_bg  Local background level around the star.
/// @param sd_bg    Standard deviation of the local background.
/// @return StarAspect with aspect = 999 on failure (fewer than 4 pixels or
///         aspect > 5).
[[nodiscard]] StarAspect measure_star_aspect(const ImageArray& img,
                                             double x1, double y1, int rs,
                                             double star_bg, double sd_bg);

/// Filter an array of HFD values using nearest-neighbour median smoothing.
///
/// For each star, the two spatially closest neighbours are found and the
/// median of the three HFD values replaces the original.  On return,
/// @p hfd_values[2][i] holds the filtered value.
///
/// @param hfd_values  StarList where [0] = X, [1] = Y, [2] = HFD * 100.
///                    Modified in-place.
/// @param nr          Number of stars to process (may be less than
///                    hfd_values[0].size()).
/// @param[out] mean       Arithmetic mean of the filtered HFD values.
/// @param[out] min_value  Minimum filtered HFD value.
/// @param[out] max_value  Maximum filtered HFD value.
void filter_hfd(StarList& hfd_values, int nr,
                float& mean, float& min_value, float& max_value);

/// Wrap an angle into the symmetric interval [-range/2, +range/2].
///
/// @param x     The angle (or value) to wrap.
/// @param range The total span of the interval (e.g. 360 for degrees,
///              2*pi for radians).
/// @return The wrapped value.
[[nodiscard]] constexpr double fnmodulo2(double x, double range) noexcept
{
	while (x < -range / 2.0)
		x += range;
	while (x > range / 2.0)
		x -= range;
	return x;
}

}  // namespace astap::analysis
