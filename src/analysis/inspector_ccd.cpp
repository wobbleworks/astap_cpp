/// @file inspector_ccd.cpp
/// CCD inspector analysis driver — the computational core of
/// unit_inspector_plot.pas CCDinspector_analyse. Detects stars across the whole
/// frame (four-level detection retry ladder, like the solver), measures each
/// star's HFD (or elongation), and reports the 90th-percentile value. Kept in
/// its own translation unit so the pure helpers in inspector.cpp
/// (measure_star_aspect, filter_hfd) stay unit-testable without linking the
/// HFD / background-estimation machinery.
///
/// The GUI output of the original (HFD Voronoi/contour maps, elongation vector
/// overlays, on-image text) is omitted; this returns the computed data only.
///
/// Copyright (C) 2018, 2022 by Han Kleijn, www.hnsky.org (original Pascal).
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this file,
/// You can obtain one at https://mozilla.org/MPL/2.0/.

#include "inspector.h"

#include "../core/photometry.h"   // HFD, get_background, HfdResult, HfdScratch
#include "../core/util.h"         // smedian

#include <algorithm>
#include <cmath>
#include <numbers>
#include <span>
#include <vector>

namespace astap::analysis {

InspectorAnalysis ccd_inspector_analyse(const ImageArray& img,
                                        const astap::Header& head,
                                        char detype, bool aspect) {
	InspectorAnalysis out;
	if (head.naxis == 0 || img.empty() || img[0].empty() || img[0][0].empty()) {
		return out;  // no image loaded
	}

	const int width  = head.width;
	const int height = head.height;
	constexpr int max_stars = 1000;

	// hfd_values rows: [0]=x, [1]=y, [2]=HFD*1000 (or aspect*1000), [3]=orientation.
	StarList hfd_values(4);

	astap::Background bck{};
	astap::core::get_background(0, img, /*calc_hist=*/false,
	                            /*calc_noise_level=*/true, bck);

	const double data_max = head.datamax_org - 1.0;

	// -1 = star-free, 1 = occupied by an already-counted star. Reset per retry.
	std::vector<std::vector<int>> img_sa(
		static_cast<std::size_t>(height), std::vector<int>(static_cast<std::size_t>(width), -1));

	astap::core::HfdScratch scratch{};
	int nhfd = 0;
	int retries = 3;  // up to four passes at progressively lower detection levels

	do {
		double detection_level = 0.0;
		if (retries == 3) {
			if (bck.star_level > 30 * bck.noise_level) detection_level = bck.star_level;
			else retries = 2;
		}
		if (retries == 2) {
			if (bck.star_level2 > 30 * bck.noise_level) detection_level = bck.star_level2;
			else retries = 1;
		}
		if (retries == 1) detection_level = 30 * bck.noise_level;
		if (retries == 0) detection_level = 7 * bck.noise_level;

		nhfd = 0;
		for (auto& row : hfd_values) row.clear();
		for (auto& row : img_sa) std::fill(row.begin(), row.end(), -1);

		for (int fitsY = 0; fitsY < height - 1; ++fitsY) {
			for (int fitsX = 0; fitsX < width - 1; ++fitsX) {
				if (img_sa[static_cast<std::size_t>(fitsY)][static_cast<std::size_t>(fitsX)] <= 0 &&
				    img[0][static_cast<std::size_t>(fitsY)][static_cast<std::size_t>(fitsX)] - bck.backgr
				        > detection_level) {
					astap::core::HfdResult res;
					astap::core::HFD(img, fitsX, fitsY, 14, 99.0, 0.0, 1.0, res, scratch);
					double hfd1 = res.hfd;
					if (hfd1 >= 1.3 && res.snr > 30 && hfd1 < 99) {
						// Mark a 5*HFD disk occupied (factor 5 for the inspector,
						// vs 3 in the solver) to prevent double detections.
						const int radius = static_cast<int>(std::lround(5.0 * hfd1));
						const int sqr_radius = radius * radius;
						const int xci = static_cast<int>(std::lround(res.xc));
						const int yci = static_cast<int>(std::lround(res.yc));
						for (int n = -radius; n <= radius; ++n) {
							for (int m = -radius; m <= radius; ++m) {
								const int jj = n + yci;
								const int ii = m + xci;
								if (jj >= 0 && ii >= 0 && jj < height && ii < width &&
								    (m * m + n * n) <= sqr_radius) {
									img_sa[static_cast<std::size_t>(jj)][static_cast<std::size_t>(ii)] = 1;
								}
							}
						}

						int orientation = 0;
						if (aspect) {
							const auto sa = measure_star_aspect(
								img, res.xc, res.yc,
								static_cast<int>(std::lround(hfd1 * 1.5)),
								scratch.star_bg, scratch.sd_bg);
							hfd1 = sa.aspect;  // hfd1 now carries the elongation
							orientation = sa.orientation;
						}

						if (hfd1 != 999) {
							const int cx = static_cast<int>(std::lround(res.xc));
							const int cy = static_cast<int>(std::lround(res.yc));
							bool not_saturated = true;
							if (cx - 1 >= 0 && cx + 1 < width && cy - 1 >= 0 && cy + 1 < height) {
								for (int dy = -1; dy <= 1 && not_saturated; ++dy) {
									for (int dx = -1; dx <= 1; ++dx) {
										if (img[0][static_cast<std::size_t>(cy + dy)]
										        [static_cast<std::size_t>(cx + dx)] >= data_max) {
											not_saturated = false;
											break;
										}
									}
								}
							}
							if (not_saturated || aspect) {
								hfd_values[0].push_back(res.xc);
								hfd_values[1].push_back(res.yc);
								hfd_values[2].push_back(hfd1 * 1000.0);
								hfd_values[3].push_back(static_cast<double>(orientation));
								++nhfd;
							}
						}
					}
				}
			}
		}
		--retries;
	} while (!(nhfd >= max_stars || retries < 0));

	if (nhfd < 10) {
		return out;  // ok stays false: "only N useful stars"
	}

	// For a contour/Voronoi map the values are smoothed by nearest-neighbour
	// median first; this also fills mean/min/max.
	if (detype != ' ') {
		filter_hfd(hfd_values, nhfd, out.mean, out.min_value, out.max_value);
	}

	// 90th percentile: "10% of the HFD measurements are worse than this".
	std::vector<double> hfds(hfd_values[2].begin(),
	                         hfd_values[2].begin() + nhfd);
	std::sort(hfds.begin(), hfds.end());
	const auto idx = static_cast<std::size_t>(std::lround((nhfd - 1) * 0.9));
	const double med = hfds[idx];

	out.ok = true;
	out.nhfd = nhfd;
	out.median_90 = med / 1000.0;
	out.stars = std::move(hfd_values);
	return out;
}

namespace {

/// SMedian over the first @p n entries of a region's HFD list.
[[nodiscard]] double region_median(const std::vector<double>& v, int n) {
	return astap::core::smedian(std::span<const double>(v.data(), static_cast<std::size_t>(n)), n);
}

}  // namespace

InspectorTilt ccd_inspector(const ImageArray& img, const astap::Header& head,
                            double snr_min, bool triangle,
                            double measuring_angle) {
	InspectorTilt out;
	if (head.naxis == 0 || img.empty() || img[0].empty() || img[0][0].empty()) {
		return out;
	}

	const int width  = head.width;
	const int height = head.height;
	constexpr int max_stars = 500;
	const double hfd_min = 0.8;  // ignore hot pixels (two-pixel minimum)

	astap::Background bck{};
	astap::core::get_background(0, img, /*calc_hist=*/false,
	                            /*calc_noise_level=*/true, bck);
	const double data_max = head.datamax_org - 1.0;

	std::vector<double> hfd_list, fwhm_list, star_x, star_y;
	std::vector<std::vector<int>> img_sa(
		static_cast<std::size_t>(height), std::vector<int>(static_cast<std::size_t>(width), -1));
	astap::core::HfdScratch scratch{};
	int nhfd = 0;
	int retries = 3;

	do {
		double detection_level = 0.0;
		if (retries == 3) {
			if (bck.star_level > 30 * bck.noise_level) detection_level = bck.star_level;
			else retries = 2;
		}
		if (retries == 2) {
			if (bck.star_level2 > 30 * bck.noise_level) detection_level = bck.star_level2;
			else retries = 1;
		}
		if (retries == 1) detection_level = 30 * bck.noise_level;
		if (retries == 0) detection_level = 7 * bck.noise_level;

		nhfd = 0;
		hfd_list.clear(); fwhm_list.clear(); star_x.clear(); star_y.clear();
		for (auto& row : img_sa) std::fill(row.begin(), row.end(), -1);

		for (int fitsY = 0; fitsY < height - 1; ++fitsY) {
			for (int fitsX = 0; fitsX < width - 1; ++fitsX) {
				if (img_sa[static_cast<std::size_t>(fitsY)][static_cast<std::size_t>(fitsX)] <= 0 &&
				    img[0][static_cast<std::size_t>(fitsY)][static_cast<std::size_t>(fitsX)] - bck.backgr
				        > detection_level) {
					astap::core::HfdResult res;
					astap::core::HFD(img, fitsX, fitsY, 14, 99.0, 0.0, 1.0, res, scratch);
					if (res.hfd <= 30.0 && res.snr > snr_min && res.hfd > hfd_min) {
						const int radius = static_cast<int>(std::lround(3.0 * res.hfd));
						const int sqr_radius = radius * radius;
						const int xci = static_cast<int>(std::lround(res.xc));
						const int yci = static_cast<int>(std::lround(res.yc));
						for (int n = -radius; n <= radius; ++n) {
							for (int m = -radius; m <= radius; ++m) {
								const int jj = n + yci;
								const int ii = m + xci;
								if (jj >= 0 && ii >= 0 && jj < height && ii < width &&
								    (m * m + n * n) <= sqr_radius) {
									img_sa[static_cast<std::size_t>(jj)][static_cast<std::size_t>(ii)] = 1;
								}
							}
						}

						bool not_saturated = true;
						if (xci - 1 >= 0 && xci + 1 < width && yci - 1 >= 0 && yci + 1 < height) {
							for (int dy = -1; dy <= 1 && not_saturated; ++dy) {
								for (int dx = -1; dx <= 1; ++dx) {
									if (img[0][static_cast<std::size_t>(yci + dy)]
									        [static_cast<std::size_t>(xci + dx)] >= data_max) {
										not_saturated = false;
										break;
									}
								}
							}
						}
						if (not_saturated) {
							hfd_list.push_back(res.hfd);
							fwhm_list.push_back(res.fwhm);
							star_x.push_back(res.xc);
							star_y.push_back(res.yc);
							++nhfd;
						}
					}
				}
			}
		}
		--retries;
	} while (!(nhfd >= max_stars || retries < 0));

	if (nhfd == 0) {
		return out;
	}

	out.nhfd = nhfd;
	out.hfd_median  = region_median(hfd_list, nhfd);
	out.fwhm_median = region_median(fwhm_list, nhfd);

	// Bin stars into the 3x3 grid (or three sectors) + outer ring.
	const int smx = width / 2;
	const int smy = height / 2;
	const double ring_thresh = 0.75 * 0.75 *
		(static_cast<double>(smx) * smx + static_cast<double>(smy) * smy);
	std::vector<double> r11, r21, r31, r12, r22, r32, r13, r23, r33, ring;

	double screw1 = 0, screw2 = 0, screw3 = 0;
	if (triangle) {
		screw1 = fnmodulo2(measuring_angle, 360.0);
		screw2 = fnmodulo2(measuring_angle + 120.0, 360.0);
		screw3 = fnmodulo2(measuring_angle - 120.0, 360.0);
	}

	for (int i = 0; i < nhfd; ++i) {
		const double hfd1 = hfd_list[static_cast<std::size_t>(i)];
		const int starX = static_cast<int>(std::lround(star_x[static_cast<std::size_t>(i)]));
		const int starY = static_cast<int>(std::lround(star_y[static_cast<std::size_t>(i)]));

		const double dx = starX - smx;
		const double dy = starY - smy;
		if (dx * dx + dy * dy > ring_thresh) {
			ring.push_back(hfd1);
		}

		if (!triangle) {
			const double w3 = width / 3.0, w23 = 2.0 * width / 3.0;
			const double h3 = height / 3.0, h23 = 2.0 * height / 3.0;
			if (starX < w3  && starY < h3)  r11.push_back(hfd1);
			if (starX > w23 && starY < h3)  r31.push_back(hfd1);
			if (starX > w23 && starY > h23) r33.push_back(hfd1);
			if (starX < w3  && starY > h23) r13.push_back(hfd1);
			if (starX > w3 && starX < w23 && starY > h23) r23.push_back(hfd1);
			if (starX < w3 && starY > h3 && starY < h23)  r12.push_back(hfd1);
			if (starX > w3 && starX < w23 && starY > h3 && starY < h23) r22.push_back(hfd1);
			if (starX > w23 && starY > h3 && starY < h23) r32.push_back(hfd1);
			if (starX > w3 && starX < w23 && starY < h3)  r21.push_back(hfd1);
		} else {
			const double xc = starX - smx;
			const double yc = starY - smy;
			// atan2(x, y): angle from the Y axis (x/y swapped, as in the original).
			const double theangle = std::atan2(xc, yc) * 180.0 / std::numbers::pi;
			// NOTE: faithful to unit_inspector_plot.pas, which sums x^2 twice here.
			const double sqrradius = xc * xc + xc * xc;
			const double theradius = std::sqrt(sqrradius);
			const double outer = 0.75 * 0.75 *
				(static_cast<double>(smx) * smx + static_cast<double>(smy) * smy);
			const double inner = 0.25 * 0.25 *
				(static_cast<double>(smx) * smx + static_cast<double>(smy) * smy);
			if (sqrradius <= outer) {
				if (sqrradius >= inner) {
					if (std::abs(fnmodulo2(theangle - screw1, 360.0)) < 30 && theradius < smy) r11.push_back(hfd1);
					if (std::abs(fnmodulo2(theangle - screw2, 360.0)) < 30 && theradius < smy) r21.push_back(hfd1);
					if (std::abs(fnmodulo2(theangle - screw3, 360.0)) < 30 && theradius < smy) r31.push_back(hfd1);
				} else {
					r22.push_back(hfd1);
				}
			}
		}
	}

	// Field curvature: outer-ring minus centre HFD.
	if (!r22.empty() && !ring.empty()) {
		out.median_22 = region_median(r22, static_cast<int>(r22.size()));
		out.median_outer_ring = region_median(ring, static_cast<int>(ring.size()));
		out.off_axis_aberration = out.median_outer_ring - out.median_22;
		out.has_off_axis = true;
	}

	// Sensor tilt: worst-minus-best region median.
	double median_best = 0.0, median_worst = 0.0;
	bool have_tilt = false;
	if (triangle && !r11.empty() && !r21.empty() && !r31.empty()) {
		out.median_11 = region_median(r11, static_cast<int>(r11.size()));
		out.median_21 = region_median(r21, static_cast<int>(r21.size()));
		out.median_31 = region_median(r31, static_cast<int>(r31.size()));
		median_best  = std::min({out.median_11, out.median_21, out.median_31});
		median_worst = std::max({out.median_11, out.median_21, out.median_31});
		have_tilt = true;
	} else if (!triangle && !r11.empty() && !r21.empty() && !r31.empty() &&
	           !r12.empty() && !r32.empty() && !r13.empty() && !r23.empty() && !r33.empty()) {
		out.median_11 = region_median(r11, static_cast<int>(r11.size()));
		out.median_21 = region_median(r21, static_cast<int>(r21.size()));
		out.median_31 = region_median(r31, static_cast<int>(r31.size()));
		out.median_12 = region_median(r12, static_cast<int>(r12.size()));
		out.median_32 = region_median(r32, static_cast<int>(r32.size()));
		out.median_13 = region_median(r13, static_cast<int>(r13.size()));
		out.median_23 = region_median(r23, static_cast<int>(r23.size()));
		out.median_33 = region_median(r33, static_cast<int>(r33.size()));
		// Tilt is measured from the four corners.
		median_best  = std::min({out.median_13, out.median_33, out.median_11, out.median_31});
		median_worst = std::max({out.median_13, out.median_33, out.median_11, out.median_31});
		have_tilt = true;
	}

	out.ok = true;
	if (have_tilt) {
		const auto tilt = classify_tilt(median_best, median_worst, out.hfd_median);
		out.tilt_hfd = tilt.tilt_hfd;
		out.tilt_percent = tilt.tilt_percent;
		out.classification = tilt.classification;
	}
	return out;
}

}  // namespace astap::analysis
