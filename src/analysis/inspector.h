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

#include <string>

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

/// Result of @ref ccd_inspector_analyse: the detected star list plus the
/// reported summary statistics. The GUI overlays of the original (HFD Voronoi/
/// contour maps, elongation vectors, on-image text) are not produced.
struct InspectorAnalysis {
	bool     ok{false};       ///< false when no image is loaded or < 10 useful stars found.
	int      nhfd{0};         ///< Number of useful stars.
	double   median_90{0.0};  ///< 90th-percentile HFD (or aspect ratio), in real units.
	float    mean{0.0f};      ///< Mean of the filtered values (valid when detype != ' ').
	float    min_value{0.0f}; ///< Minimum filtered value.
	float    max_value{0.0f}; ///< Maximum filtered value.
	StarList stars{};         ///< [0]=x, [1]=y, [2]=HFD*1000 (or aspect*1000), [3]=orientation.
};

/// Detect stars across the whole frame and summarise focus/tilt quality — the
/// computational core of unit_inspector_plot.pas CCDinspector_analyse. Uses the
/// same four-level detection retry ladder as the solver, measures each star's
/// HFD (and optionally its elongation), and reports the 90th-percentile value
/// ("10% of measurements are worse than this").
///
/// Precondition: the image histogram must already be current for @p img (the
/// original relies on the load-time histogram; call get_hist(0, img) first if
/// unsure). The map drawing, vector overlays and on-image text are omitted.
///
/// @param img    Image buffer (channel 0 used).
/// @param head   Image header (dimensions, datamax).
/// @param detype ' ' = plain HFD, 'V' = Voronoi map, '2' = contour map. Only
///               selects whether nearest-neighbour median smoothing
///               (@ref filter_hfd) is applied to the values; the maps
///               themselves are GUI-only.
/// @param aspect When true, stores each star's aspect ratio (elongation)
///               instead of its HFD, and skips the saturation rejection.
[[nodiscard]] InspectorAnalysis ccd_inspector_analyse(const ImageArray& img,
                                                      const astap::Header& head,
                                                      char detype, bool aspect);

/// Tilt classification derived from the best/worst region HFD medians.
struct TiltResult {
	double      tilt_hfd{};       ///< median_worst - median_best (HFD units).
	double      tilt_percent{};   ///< 100 * tilt_hfd / hfd_median.
	std::string classification;   ///< none / almost none / mild / moderate / severe / extreme.
};

/// Classify sensor tilt from the best and worst region HFD medians. Pure
/// (unit_inspector_plot.pas:575-587); thresholds 5/10/15/20/30 percent.
[[nodiscard]] TiltResult classify_tilt(double median_best, double median_worst,
                                       double hfd_median);

/// Full CCD-inspector result: focus (median HFD/FWHM), sensor tilt (best/worst
/// region HFD spread), and field curvature (outer-ring minus centre HFD). The
/// GUI octagon/vector overlay is not produced.
struct InspectorTilt {
	bool        ok{false};                 ///< false when no image or too few stars in the needed regions.
	int         nhfd{0};                   ///< Total detected stars.
	double      hfd_median{0.0};           ///< Median HFD over all stars.
	double      fwhm_median{0.0};          ///< Median FWHM over all stars.
	double      tilt_hfd{100.0};           ///< median_worst - median_best (100 = not computed).
	double      tilt_percent{0.0};         ///< 100 * tilt_hfd / hfd_median.
	std::string classification;            ///< Tilt severity label (empty if not computed).
	bool        has_off_axis{false};       ///< Whether off_axis_aberration is valid.
	double      off_axis_aberration{0.0};  ///< Outer-ring median HFD minus centre (region 22).
	// Per-region HFD medians (0 = not computed). Layout (FITS 1,1 = bottom-left):
	//   13 23 33 / 12 22 32 / 11 21 31, plus the outer ring.
	double      median_11{}, median_21{}, median_31{};
	double      median_12{}, median_22{}, median_32{};
	double      median_13{}, median_23{}, median_33{};
	double      median_outer_ring{};
};

/// Detect stars across the frame, bin them into a 3x3 region grid (or three
/// 120-degree sectors when @p triangle is set) and report focus, sensor tilt
/// and field curvature — the computational core of unit_inspector_plot.pas
/// CCDinspector. GUI drawing is omitted. Precondition: the histogram is current
/// for @p img (call get_hist(0, img) first if unsure).
///
/// @param snr_min        Minimum SNR for a usable star (typically 30).
/// @param triangle       Three-sector "screw" tilt mode instead of the 9-grid.
/// @param measuring_angle Sector reference angle (degrees) for triangle mode.
[[nodiscard]] InspectorTilt ccd_inspector(const ImageArray& img,
                                          const astap::Header& head,
                                          double snr_min, bool triangle,
                                          double measuring_angle);

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
