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

#include <algorithm>
#include <cmath>
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

}  // namespace astap::analysis
