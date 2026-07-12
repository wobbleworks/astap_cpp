/// @file contour_streak.cpp
/// Preprocessing wrapper for satellite-streak detection — the colour/Bayer →
/// mono conversion, Gaussian blur and background estimation from
/// unit_contour.pas contour(). The dependency-free tracing core lives in
/// contour.cpp (detect_streaks); this wrapper is split out so that core stays
/// unit-testable without linking the blur/background/crop machinery.
///
/// Copyright (C) 2023, 2024 by Han Kleijn, www.hnsky.org (original Pascal).
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this file,
/// You can obtain one at https://mozilla.org/MPL/2.0/.

#include "contour.h"

#include "../core/globals.h"                  // bayerpat
#include "../core/imaging.h"                  // get_hist
#include "../core/photometry.h"               // get_background
#include "../solving/astrometric_solving.h"   // binX1_crop, binX2_crop
#include "../stacking/gaussian_blur.h"        // gaussian_blur2

namespace astap::analysis {

std::vector<Streak> contour(const ImageArray& img, astap::Header& head,
                            double blur, double sigmafactor, int detection_grid) {
	int binning = 1;
	ImageArray img_bk;
	bool restore_his = false;

	if (head.naxis3 > 1) {
		// Colour image: crop to mono, no binning.
		astap::solving::binX1_crop(1, img, img_bk);
		astap::core::get_hist(0, img_bk);
		restore_his = true;
	} else if (!astap::bayerpat.empty()) {
		// Raw Bayer image: combine 2x2 pixels to mono.
		astap::solving::binX2_crop(1, img, img_bk);
		astap::core::get_hist(0, img_bk);
		restore_his = true;
		binning = 2;
	} else {
		// Mono: copy so the blur below does not mutate the caller's image
		// (the original shares the buffer and blurs it in place).
		img_bk = img;
	}

	astap::stacking::gaussian_blur2(img_bk, blur);
	astap::Background bck{};
	astap::core::get_background(0, img_bk, /*calc_hist=*/false,
	                            /*calc_noise_level=*/true, bck);

	const double detection_level = sigmafactor * bck.noise_level + bck.backgr;
	const int grid = detection_grid / binning;

	auto streaks = detect_streaks(img_bk, detection_level, binning, grid);

	if (restore_his) {
		astap::core::get_hist(0, img);  // restore the original image's histogram
	}
	return streaks;
}

}  // namespace astap::analysis
