///----------------------------------------
///      @file star_detector.h
///   @ingroup ASTAP++
///     @brief GUI-side star detection using the engine's background /
///            HFD primitives.
///   @details Mirrors the detection loop from @c astap::stacking::analyse_image
///            but returns the full per-star list (not just a count and median)
///            so the GUI can render markers on top of the image. Keeps the
///            engine's API untouched.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include "../../src/types.h"

#include <vector>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @struct DetectedStar
/// @brief A single star located in the image, with its HFD / SNR / flux.
/// @details Pixel coordinates follow the FITS convention: 1-based, origin
///          at bottom-left, y increasing upward.
///----------------------------------------

struct DetectedStar {
	double x = 0.0;
	double y = 0.0;
	double hfd = 0.0;
	double fwhm = 0.0;
	double snr = 0.0;
	double flux = 0.0;
};

///----------------------------------------
/// @struct DetectionResult
/// @brief Bundle of per-star detections plus summary statistics.
///----------------------------------------

struct DetectionResult {
	std::vector<DetectedStar> stars;
	double medianHfd = 0.0;
	double background = 0.0;
	double noise = 0.0;
	double starLevel = 0.0;        // threshold for "bright" stars (HFD ≈ 2.25)
	double starLevel2 = 0.0;       // threshold for "dim" stars (HFD ≈ 4.5)
	double detectionLevel = 0.0;   // threshold actually used
	int retriesUsed = 0;           // 3→0 (3=star_level, 0=7·σ fallback)
	std::uint64_t candidatesAbove = 0;     // pixels that cleared detection_level
	std::uint64_t candidatesRejected = 0;  // candidates HFD rejected
};

///----------------------------------------
/// @brief Detect stars in @p img using the engine's background and HFD
///        primitives.
/// @details Thread-safe relative to the GUI: reads @p img, does not mutate
///          any global state beyond the detection scratch (histogram and
///          per-HFD scratch). Call from a worker thread.
/// @param img     Loaded image data.
/// @param snr_min Minimum signal-to-noise ratio to accept a star.
/// @return Detected stars plus summary statistics.
///----------------------------------------

[[nodiscard]] DetectionResult detect_stars(const astap::ImageArray& img,
                                           double snr_min = 5.0);

} // namespace astap::gui
