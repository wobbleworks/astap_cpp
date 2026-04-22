///----------------------------------------
///     @file photometry_calibration_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for `calibrate_flux` (src/core/photometry.cpp).
///  @details Exercises the end-to-end flux-calibration pipeline with an
///           injectable @ref StarSource: a synthetic image with Gaussian
///           stars at known pixel positions is paired with a catalog source
///           yielding matching (RA, Dec, magnitude) tuples, and the
///           resulting @c head.mzero is asserted against the expected
///           value.
///   @author Created by John Stephen on 4/19/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "core/photometry.h"
#include "core/fits.h"
#include "core/wcs.h"

#include <cmath>
#include <array>
#include <cstdint>
#include <numbers>
#include <string>
#include <vector>

using namespace astap::core;
using astap::Header;
using astap::ImageArray;

///----------------------------------------
/// MARK: Test-local stubs
///
/// photometry.cpp / wcs.cpp reference a handful of externals that live in
/// other translation units the test does not link in. Stub them here.
///----------------------------------------

namespace astap {

// Globals normally provided by globals.cpp — stubbed with harmless defaults.
std::vector<std::string> memo1_lines;
int                      max_stars_setting = 500;
int                      nrbits            = 16;
std::string              filename2;

}  // namespace astap

namespace astap::core {

// imaging.cpp-owned histogram state used by get_background. Provide
// zeroed stubs so photometry.cpp's background-dependent paths can link.
// Types MUST match imaging.h exactly — MSVC bakes the full type into the
// mangled symbol name, so `int[3][65536]` does not link against the
// `std::array<std::array<std::uint32_t, 65536>, 3>` reference.
std::array<std::array<std::uint32_t, 65536>, 3> histogram{};
std::array<int, 3>                              his_mean{};

// get_hist is part of imaging.cpp; calibrate_flux doesn't use it but
// photometry.cpp contains other functions that do (find_star_center etc.)
// that may get pulled in by the linker. Stub as a no-op.
void get_hist(int, const ImageArray&) {}

// WCS externals used by pixel_to_celestial / celestial_to_pixel.
void dsspos(double, double, double& ra, double& dec) { ra = 0; dec = 0; }
void EQU_GAL(double, double, double& l, double& b)   { l  = 0; b   = 0; }
bool calculate_az_alt_basic(double, double, double& az, double& alt) {
	az = alt = 0;
	return false;
}
double position_angle(double, double, double, double) { return 0.0; }

// get_background — histogram-aware overload from imaging.h. photometry.cpp's
// calibrate_photometry path uses it, but calibrate_flux does not; stub so
// the overload-resolution / linker are satisfied.
void get_background(int, const ImageArray&, bool, bool,
                    astap::Background& bck,
                    int (&)[3][65536], int (&)[3],
                    int, int,
                    const std::string&,
                    std::vector<std::string>&) {
	bck = astap::Background{};
}

// analyse_image: referenced by calibrate_photometry (not calibrate_flux).
void analyse_image(const ImageArray&, Header&, int, int,
                   int& hfd_counter, astap::Background& bck, double& hfd_med) {
	hfd_counter = 0;
	bck = astap::Background{};
	hfd_med = 0.0;
}

// plot_and_measure_stars: referenced by calibrate_photometry (not calibrate_flux).
void plot_and_measure_stars(const ImageArray&, std::vector<std::string>&,
                             Header&, bool, bool, bool) {}

// update_float / update_text: real impls live in fits.cpp which drags in the
// whole FITS loader. calibrate_flux only uses them to write human-readable
// header cards; noop stubs suffice (the test asserts head fields directly).
void update_float(std::vector<std::string>&, std::string_view, std::string_view,
                   bool, double) {}
void update_text(std::vector<std::string>&, std::string_view, std::string_view) {}

// Smedian (capital S): simple median helper, used for the limiting-magnitude
// branch. Provide a minimal impl so we do not need link_stubs.cpp.
double Smedian(std::vector<double>& list, int len) {
	if (len <= 0) return 0.0;
	if (len == 1) return list[0];
	std::sort(list.begin(), list.begin() + len);
	if (len % 2 == 1) {
		return list[(len - 1) / 2];
	}
	return (list[len / 2 - 1] + list[len / 2]) / 2.0;
}

}  // namespace astap::core

///----------------------------------------
/// MARK: Helpers
///----------------------------------------

/// @brief Build a simple TAN-projection Header with 1 arcsec/pixel scale.
[[nodiscard]] static Header make_header(int width, int height,
                                         double ra0_rad, double dec0_rad) {
	Header h{};
	h.width   = width;
	h.height  = height;
	h.naxis   = 2;
	h.naxis3  = 1;
	h.crpix1  = width  * 0.5;
	h.crpix2  = height * 0.5;
	constexpr double arcsec_deg = 1.0 / 3600.0;
	h.cdelt1  = -arcsec_deg;
	h.cdelt2  =  arcsec_deg;
	h.ra0     = ra0_rad;
	h.dec0    = dec0_rad;
	h.cd1_1   = h.cdelt1;
	h.cd1_2   = 0.0;
	h.cd2_1   = 0.0;
	h.cd2_2   = h.cdelt2;
	h.datamax_org = 65535.0;
	h.mzero_radius = 99.0;  // extended-objects mode; matches calibrate_flux callers
	return h;
}

/// @brief Construct a single-channel image with a constant background.
[[nodiscard]] static ImageArray make_image(int width, int height, float background) {
	ImageArray img;
	img.assign(1,
		std::vector<std::vector<float>>(height,
			std::vector<float>(width, background)));
	return img;
}

/// @brief Paint a Gaussian star centered at (cx, cy) with std-dev @p sigma
///        and peak amplitude @p peak on top of whatever is currently there.
static void paint_gaussian(ImageArray& img, double cx, double cy,
                           double sigma, double peak) {
	const auto height = static_cast<int>(img[0].size());
	const auto width  = static_cast<int>(img[0][0].size());
	// 5-sigma radius is more than enough.
	const auto r = static_cast<int>(std::ceil(5.0 * sigma));
	const auto inv2sig2 = 1.0 / (2.0 * sigma * sigma);

	for (int y = std::max(0, static_cast<int>(cy) - r);
	         y <= std::min(height - 1, static_cast<int>(cy) + r); ++y) {
		for (int x = std::max(0, static_cast<int>(cx) - r);
		         x <= std::min(width - 1, static_cast<int>(cx) + r); ++x) {
			const auto dx = x - cx;
			const auto dy = y - cy;
			const auto v  = peak * std::exp(-(dx*dx + dy*dy) * inv2sig2);
			img[0][y][x] = static_cast<float>(img[0][y][x] + v);
		}
	}
}

/// @brief Add a low-amplitude deterministic noise pattern so sd_bg > 0.
static void add_noise(ImageArray& img, float amplitude) {
	const auto height = static_cast<int>(img[0].size());
	const auto width  = static_cast<int>(img[0][0].size());
	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			// Deterministic pseudo-random pattern: XOR-hash of (x, y).
			const auto h = static_cast<std::uint32_t>((x * 0x9e3779b1u)
			                                        ^ (y * 0x85ebca77u));
			const auto rnd = static_cast<double>(h) / 4294967295.0;  // [0, 1]
			img[0][y][x] = static_cast<float>(img[0][y][x] + amplitude * (rnd - 0.5));
		}
	}
}

/// @brief One star for the injected catalog source.
struct TestStar {
	double px;          ///< @brief FITS pixel X (1-based).
	double py;          ///< @brief FITS pixel Y (1-based).
	double magn_times_10; ///< @brief Catalog magnitude * 10.
	double Bp_Rp{999.0};  ///< @brief Gaia Bp-Rp * 10; 999 = mono database.
};

/// @brief Build a @ref StarSource that yields the given stars in order,
///        converting each (pixel_x, pixel_y) into (RA, Dec) via the header.
[[nodiscard]] static StarSource make_source(const Header& head,
                                             std::vector<TestStar> stars) {
	auto index = std::make_shared<std::size_t>(0);
	auto data  = std::make_shared<std::vector<TestStar>>(std::move(stars));
	return [=, &head](CatalogStar& out) {
		if (*index >= data->size()) {
			return false;
		}
		const auto& s = (*data)[*index];
		double ra  = 0;
		double dec = 0;
		pixel_to_celestial(head, s.px, s.py, /*formalism=*/0, ra, dec);
		out.ra    = ra;
		out.dec   = dec;
		out.magn  = s.magn_times_10;
		out.Bp_Rp = s.Bp_Rp;
		++(*index);
		return true;
	};
}

///----------------------------------------
/// MARK: Baseline success
///----------------------------------------

TEST_CASE("calibrate_flux: 8 identical Gaussian stars produce a finite MZERO") {
	// Arrange: 500x500 image with uniform background, 8 Gaussian stars on a
	// 4x2 grid well inside the frame. All stars have identical pixel flux
	// and identical catalog magnitudes, so the per-star flux ratio should
	// be constant and the SEM tiny.
	constexpr int kW = 500;
	constexpr int kH = 500;
	constexpr float kBackground = 100.0f;
	constexpr double kSigma = 2.0;
	constexpr double kPeak  = 8000.0;

	auto head = make_header(kW, kH, 1.0, 0.3);
	auto img  = make_image(kW, kH, kBackground);
	add_noise(img, 10.0f);   // ensure sd_bg > 0

	std::vector<TestStar> stars;
	// Pixel positions for 8 stars (FITS 1-based).
	constexpr double positions[][2] = {
		{100,  100}, {200,  100}, {300,  100}, {400,  100},
		{100,  400}, {200,  400}, {300,  400}, {400,  400},
	};
	for (const auto& p : positions) {
		paint_gaussian(img, p[0] - 1.0, p[1] - 1.0, kSigma, kPeak);
		stars.push_back({.px = p[0], .py = p[1], .magn_times_10 = 100.0});
	}

	// Act.
	std::vector<std::string> memo;
	auto result = calibrate_flux(img, head, memo, make_source(head, stars),
	                              "V", /*annulus=*/14,
	                              /*aperture_setting=*/0.0,
	                              /*report_lim_magn=*/false);

	// Assert: detected all 8 stars, got a finite MZERO with small SEM.
	CHECK(result.success == true);
	CHECK(result.stars_measured == 8);
	CHECK(std::isfinite(head.mzero));
	CHECK(head.mzero > 10.0);
	CHECK(head.mzero < 30.0);
	CHECK(head.passband_database == "V");
	CHECK(result.standard_error_mean >= 0.0);
	CHECK(result.standard_error_mean < 0.2);   // identical stars → SEM tiny
}

///----------------------------------------
/// MARK: Flux-doubling shifts MZERO by +2.5 * log10(2)
///----------------------------------------

TEST_CASE("calibrate_flux: doubling pixel flux shifts MZERO by 2.5*log10(2)") {
	constexpr int kW = 500;
	constexpr int kH = 500;
	constexpr float kBackground = 100.0f;
	constexpr double kSigma = 2.0;

	auto run_with_peak = [&](double peak) {
		auto head = make_header(kW, kH, 1.0, 0.3);
		auto img  = make_image(kW, kH, kBackground);
		add_noise(img, 10.0f);
		std::vector<TestStar> stars;
		constexpr double positions[][2] = {
			{100,  100}, {200,  100}, {300,  100}, {400,  100},
			{100,  400}, {200,  400}, {300,  400}, {400,  400},
		};
		for (const auto& p : positions) {
			paint_gaussian(img, p[0] - 1.0, p[1] - 1.0, kSigma, peak);
			stars.push_back({.px = p[0], .py = p[1], .magn_times_10 = 100.0});
		}
		std::vector<std::string> memo;
		auto r = calibrate_flux(img, head, memo, make_source(head, stars),
		                        "V", 14, 0.0, false);
		return std::pair{r, head.mzero};
	};

	const auto [r1, mzero1] = run_with_peak(4000.0);
	const auto [r2, mzero2] = run_with_peak(8000.0);

	REQUIRE(r1.success);
	REQUIRE(r2.success);

	// Doubling flux → shift mzero by +2.5 * log10(2) ≈ +0.7526.
	constexpr double expected_shift = 2.5 * 0.30102999566398114;
	CHECK((mzero2 - mzero1) == doctest::Approx(expected_shift).epsilon(0.05));
}

///----------------------------------------
/// MARK: Rejection paths
///----------------------------------------

TEST_CASE("calibrate_flux: red stars (Bp_Rp > 12) are skipped") {
	constexpr int kW = 500;
	constexpr int kH = 500;
	constexpr double kSigma = 2.0;
	constexpr double kPeak  = 8000.0;

	auto head = make_header(kW, kH, 1.0, 0.3);
	auto img  = make_image(kW, kH, 100.0f);
	add_noise(img, 10.0f);

	// 8 stars painted, but one tagged as "too red" (Bp_Rp = 15 means 1.5 mag).
	std::vector<TestStar> stars;
	constexpr double positions[][2] = {
		{100,  100}, {200,  100}, {300,  100}, {400,  100},
		{100,  400}, {200,  400}, {300,  400}, {400,  400},
	};
	for (std::size_t i = 0; i < 8; ++i) {
		paint_gaussian(img, positions[i][0] - 1.0, positions[i][1] - 1.0,
		               kSigma, kPeak);
		stars.push_back({
			.px = positions[i][0], .py = positions[i][1],
			.magn_times_10 = 100.0,
			.Bp_Rp = (i == 0) ? 15.0 : 999.0,   // first star is too red
		});
	}

	std::vector<std::string> memo;
	auto result = calibrate_flux(img, head, memo, make_source(head, stars),
	                              "V", 14, 0.0, false);

	CHECK(result.success);
	CHECK(result.stars_measured == 7);  // red star skipped
}

TEST_CASE("calibrate_flux: Bp_Rp == -128 (unreliable Johnson-V) is skipped") {
	constexpr int kW = 500;
	constexpr int kH = 500;
	constexpr double kSigma = 2.0;
	constexpr double kPeak  = 8000.0;

	auto head = make_header(kW, kH, 1.0, 0.3);
	auto img  = make_image(kW, kH, 100.0f);
	add_noise(img, 10.0f);

	std::vector<TestStar> stars;
	constexpr double positions[][2] = {
		{100,  100}, {200,  100}, {300,  100}, {400,  100},
		{100,  400}, {200,  400}, {300,  400}, {400,  400},
	};
	for (std::size_t i = 0; i < 8; ++i) {
		paint_gaussian(img, positions[i][0] - 1.0, positions[i][1] - 1.0,
		               kSigma, kPeak);
		stars.push_back({
			.px = positions[i][0], .py = positions[i][1],
			.magn_times_10 = 100.0,
			.Bp_Rp = (i == 0) ? -128.0 : 999.0,
		});
	}

	std::vector<std::string> memo;
	auto result = calibrate_flux(img, head, memo, make_source(head, stars),
	                              "V", 14, 0.0, false);

	CHECK(result.success);
	CHECK(result.stars_measured == 7);
}

TEST_CASE("calibrate_flux: saturated star (center ≥ datamax - 1000) is rejected") {
	constexpr int kW = 500;
	constexpr int kH = 500;
	constexpr double kSigma = 2.0;
	constexpr double kPeak  = 8000.0;

	auto head = make_header(kW, kH, 1.0, 0.3);
	auto img  = make_image(kW, kH, 100.0f);
	add_noise(img, 10.0f);

	std::vector<TestStar> stars;
	constexpr double positions[][2] = {
		{100,  100}, {200,  100}, {300,  100}, {400,  100},
		{100,  400}, {200,  400}, {300,  400}, {400,  400},
	};
	for (std::size_t i = 0; i < 8; ++i) {
		paint_gaussian(img, positions[i][0] - 1.0, positions[i][1] - 1.0,
		               kSigma, kPeak);
		stars.push_back({.px = positions[i][0], .py = positions[i][1],
		                 .magn_times_10 = 100.0});
	}

	// Saturate the core of the first star: paint a 3x3 block at the peak.
	// datamax_org is 65535; cutoff is datamax_org - 1000 = 64535. Set to
	// 65000 to ensure the saturation check trips.
	for (int dy = -1; dy <= 1; ++dy) {
		for (int dx = -1; dx <= 1; ++dx) {
			img[0][99 + dy][99 + dx] = 65000.0f;
		}
	}

	std::vector<std::string> memo;
	auto result = calibrate_flux(img, head, memo, make_source(head, stars),
	                              "V", 14, 0.0, false);

	CHECK(result.success);
	CHECK(result.stars_measured == 7);  // saturated star rejected
}

///----------------------------------------
/// MARK: Failure paths
///----------------------------------------

TEST_CASE("calibrate_flux: fewer than 3 usable stars reports failure") {
	constexpr int kW = 500;
	constexpr int kH = 500;
	constexpr double kSigma = 2.0;
	constexpr double kPeak  = 8000.0;

	auto head = make_header(kW, kH, 1.0, 0.3);
	auto img  = make_image(kW, kH, 100.0f);
	add_noise(img, 10.0f);

	std::vector<TestStar> stars;
	constexpr double positions[][2] = {{100, 100}, {400, 400}};
	for (const auto& p : positions) {
		paint_gaussian(img, p[0] - 1.0, p[1] - 1.0, kSigma, kPeak);
		stars.push_back({.px = p[0], .py = p[1], .magn_times_10 = 100.0});
	}

	std::vector<std::string> memo;
	auto result = calibrate_flux(img, head, memo, make_source(head, stars),
	                              "V", 14, 0.0, false);

	CHECK(result.success == false);
	CHECK(result.stars_measured == 2);
	CHECK(head.mzero == 0.0);  // untouched on failure
}

TEST_CASE("calibrate_flux: empty image reports failure immediately") {
	auto head = make_header(500, 500, 1.0, 0.3);
	ImageArray img;   // empty

	std::vector<std::string> memo;
	auto result = calibrate_flux(img, head, memo,
	                              [](CatalogStar&) { return false; },
	                              "V", 14, 0.0, false);

	CHECK(result.success == false);
	CHECK(result.stars_measured == 0);
}

TEST_CASE("calibrate_flux: header without WCS reports failure") {
	auto head = make_header(500, 500, 1.0, 0.3);
	head.cd1_1 = 0.0;   // break WCS
	auto img  = make_image(500, 500, 100.0f);

	std::vector<std::string> memo;
	auto result = calibrate_flux(img, head, memo,
	                              [](CatalogStar&) { return false; },
	                              "V", 14, 0.0, false);

	CHECK(result.success == false);
	CHECK(result.stars_measured == 0);
}

///----------------------------------------
/// MARK: Passband recording
///----------------------------------------

TEST_CASE("calibrate_flux: passband_active is written to head.passband_database") {
	constexpr int kW = 500;
	constexpr int kH = 500;
	constexpr double kSigma = 2.0;
	constexpr double kPeak  = 8000.0;

	auto head = make_header(kW, kH, 1.0, 0.3);
	auto img  = make_image(kW, kH, 100.0f);
	add_noise(img, 10.0f);

	std::vector<TestStar> stars;
	constexpr double positions[][2] = {
		{100,  100}, {200,  100}, {300,  100}, {400,  100},
		{100,  400}, {200,  400}, {300,  400}, {400,  400},
	};
	for (const auto& p : positions) {
		paint_gaussian(img, p[0] - 1.0, p[1] - 1.0, kSigma, kPeak);
		stars.push_back({.px = p[0], .py = p[1], .magn_times_10 = 100.0});
	}

	std::vector<std::string> memo;
	auto result = calibrate_flux(img, head, memo, make_source(head, stars),
	                              "SG", 14, 0.0, false);

	REQUIRE(result.success);
	CHECK(head.passband_database == "SG");
}
