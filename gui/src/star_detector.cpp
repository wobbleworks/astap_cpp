///----------------------------------------
///      @file star_detector.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the GUI-side star detector.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "star_detector.h"

#include "../../src/core/photometry.h"

#include <algorithm>
#include <cmath>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

DetectionResult detect_stars(const astap::ImageArray& img, double snr_min) {
	DetectionResult out;

	if (img.empty() || img[0].empty() || img[0][0].empty()) {
		return out;
	}

	const auto height = static_cast<int>(img[0].size());
	const auto width = static_cast<int>(img[0][0].size());

	// Background + noise + star thresholds via the engine primitive.
	astap::Background bck{};
	astap::core::get_background(/*colour=*/0, img,
		/*calc_hist=*/true, /*calc_noise_level=*/true, bck);

	out.background = bck.backgr;
	out.noise = bck.noise_level;
	out.starLevel = bck.star_level;
	out.starLevel2 = bck.star_level2;

	// Visited-pixel mask — mirrors analyse_image's img_sa exclusion zone.
	auto visited = std::vector<std::vector<std::uint8_t>>(
		height, std::vector<std::uint8_t>(width, 0));

	// Retry ladder from analyse_image: try the bright-star threshold first;
	// if that finds nothing, fall back to progressively deeper noise-floor
	// thresholds. Unlike analyse_image we always try the deepest fallback
	// (7·σ) if earlier attempts didn't hit — that matches the Pascal
	// behaviour more faithfully and surfaces stars on stretched images
	// whose bright-star threshold is dominated by gradient.
	const auto hfd_min = 0.8;

	// Ladder: (threshold value, retries-label for diagnostics)
	struct Level {
		double detection;
		int label;
	};
	auto levels = std::vector<Level>{};
	if (bck.star_level > 30 * bck.noise_level) {
		levels.push_back({bck.star_level, 3});
	}
	if (bck.star_level2 > 30 * bck.noise_level) {
		levels.push_back({bck.star_level2, 2});
	}
	levels.push_back({30 * bck.noise_level, 1});
	levels.push_back({ 7 * bck.noise_level, 0});

	astap::core::HfdScratch scratch{};

	for (const auto& level : levels) {
		out.detectionLevel = level.detection;
		out.retriesUsed = level.label;
		out.candidatesAbove = 0;
		out.candidatesRejected = 0;

		// Reset visited mask and star list for this attempt.
		for (auto& row : visited) {
			std::fill(row.begin(), row.end(), 0);
		}
		out.stars.clear();

		for (auto y = 0; y < height; ++y) {
			for (auto x = 0; x < width; ++x) {
				if (visited[y][x]) {
					continue;
				}
				if ((img[0][y][x] - bck.backgr) <= level.detection) {
					continue;
				}
				++out.candidatesAbove;

				astap::core::HfdResult r;
				astap::core::HFD(img, x, y, /*rs=*/14,
					/*aperture_small=*/99.0, /*adu_e=*/0.0,
					/*xbinning=*/1.0, r, scratch);

				if (r.hfd <= 30.0 && r.snr > snr_min && r.hfd > hfd_min) {
					DetectedStar s;
					// FITS convention: 1-based pixel coordinates.
					s.x = r.xc + 1.0;
					s.y = r.yc + 1.0;
					s.hfd = r.hfd;
					s.fwhm = r.fwhm;
					s.snr = r.snr;
					s.flux = r.flux;
					out.stars.push_back(s);

					// Mask a circular exclusion zone so we don't re-detect
					// the same star from a nearby bright pixel.
					const auto radius = static_cast<int>(std::round(3.0 * r.hfd));
					const auto sqrRadius = radius * radius;
					const auto cx = static_cast<int>(std::round(r.xc));
					const auto cy = static_cast<int>(std::round(r.yc));
					for (auto n = -radius; n <= radius; ++n) {
						for (auto m = -radius; m <= radius; ++m) {
							const auto j = n + cy;
							const auto i = m + cx;
							if (j >= 0 && i >= 0 && j < height && i < width
									&& (m * m + n * n) <= sqrRadius) {
								visited[j][i] = 1;
							}
						}
					}
				} else {
					++out.candidatesRejected;
				}
			}
		}

		if (!out.stars.empty()) {
			break;
		}
	}

	// Median HFD for convenience
	if (!out.stars.empty()) {
		auto hfds = std::vector<double>(out.stars.size());
		for (auto i = std::size_t{0}; i < out.stars.size(); ++i) {
			hfds[i] = out.stars[i].hfd;
		}
		std::ranges::nth_element(hfds, hfds.begin() + hfds.size() / 2);
		out.medianHfd = hfds[hfds.size() / 2];
	}

	return out;
}

} // namespace astap::gui
