///----------------------------------------
///      @file detect_stars.cpp
///   @ingroup ASTAP++/cli
///     @brief Dev harness that runs the stacker's alignment star-detection
///            path on a single FITS file and dumps each detected star.
///   @details Not meant for users — built to iterate on find_stars tuning
///            without launching the GUI. Invoke:
///                detect_stars <fits_path>
///    @author Created by John Stephen on 4/21/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "../src/core/globals.h"
#include "../src/core/image_io.h"
#include "../src/core/photometry.h"
#include "../src/solving/astrometric_solving.h"
#include "../src/solving/star_align.h"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
	if (argc < 2) {
		std::cerr << "usage: detect_stars <fits_path>\n";
		return 1;
	}

	astap::filename2 = argv[1];
	auto ok = astap::core::load_image(astap::filename2,
	                                   astap::img_loaded,
	                                   astap::head,
	                                   astap::memo1_lines,
	                                   /*re_center=*/false,
	                                   /*plot=*/false);
	if (!ok) {
		std::cerr << "failed to load: " << astap::filename2 << '\n';
		return 2;
	}

	std::cout << "loaded: " << astap::head.width << "x" << astap::head.height
	          << " naxis3=" << astap::head.naxis3
	          << " nrbits=" << astap::nrbits << '\n';

	astap::StarList starlist;
	std::string warning;
	const auto hfd_min = std::max(0.8, astap::hfd_min_setting);
	const auto hfd_max = astap::hfd_max_setting;
	const auto max_stars = astap::max_stars_setting;

	std::cout << "settings: hfd_min=" << hfd_min
	          << " hfd_max=" << hfd_max
	          << " max_stars=" << max_stars << '\n';

	astap::solving::bin_and_find_stars(
		astap::img_loaded, /*binning=*/1, /*cropping=*/1.0,
		hfd_min, hfd_max, max_stars, /*get_hist=*/true,
		starlist, warning);

	if (!warning.empty()) { std::cout << "warn: " << warning << '\n'; }

	const auto n = starlist.empty() ? 0 : static_cast<int>(starlist[0].size());
	std::cout << "detected " << n << " stars:\n";
	for (int i = 0; i < n; ++i) {
		std::printf("  [%2d] x=%7.2f y=%7.2f\n",
		            i, starlist[0][i], starlist[1][i]);
	}

	return 0;
}
