///----------------------------------------
///     @file util_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for astap::core utility helpers (src/core/util.cpp).
///  @details Covers sort/median/MAD, angle wrapping, HSV<->RGB conversion,
///           locale-independent numeric formatting, filename classifiers,
///           filename-metadata extractors, and the encrypt/decrypt pair.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "core/util.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

using namespace astap::core;

///----------------------------------------
/// MARK: Math — quicksort / smedian / mad_median / fnmodulo
///----------------------------------------

TEST_CASE("quicksort sorts inclusive range") {
	SUBCASE("shuffled array, full range") {
		std::vector<double> a{3, 1, 4, 1, 5, 9, 2, 6, 5, 3};
		quicksort(a, 0, static_cast<int>(a.size()) - 1);
		const std::vector<double> expected{1, 1, 2, 3, 3, 4, 5, 5, 6, 9};
		CHECK(a == expected);
	}

	SUBCASE("sub-range leaves outside untouched") {
		std::vector<double> a{9, 4, 2, 1, 3, 7};
		// Sort only indices 1..4 — element 0 (9) and 5 (7) must stay put.
		quicksort(a, 1, 4);
		CHECK(a[0] == 9.0);
		CHECK(a[5] == 7.0);
		// The middle run should now be non-decreasing.
		CHECK(a[1] <= a[2]);
		CHECK(a[2] <= a[3]);
		CHECK(a[3] <= a[4]);
	}

	SUBCASE("already-sorted input is a no-op") {
		std::vector<double> a{-2.5, 0.0, 1.1, 3.14, 42.0};
		const auto copy = a;
		quicksort(a, 0, static_cast<int>(a.size()) - 1);
		CHECK(a == copy);
	}
}

TEST_CASE("smedian") {
	SUBCASE("single element returns that element") {
		const double v = 42.0;
		CHECK(smedian({&v, 1}, 1) == 42.0);
	}

	SUBCASE("odd length >= 5 averages three central values") {
		const std::vector<double> v{1, 2, 3, 4, 5};
		// Sorted median is 3; avg of three centre (2, 3, 4) is also 3.
		CHECK(smedian(v, 5) == doctest::Approx(3.0));
	}

	SUBCASE("odd length >= 5, uneven distribution") {
		const std::vector<double> v{1, 2, 3, 4, 100};
		// Three centre values (2, 3, 4) → avg = 3.
		CHECK(smedian(v, 5) == doctest::Approx(3.0));
	}

	SUBCASE("even length averages two central values") {
		const std::vector<double> v{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
		CHECK(smedian(v, 10) == doctest::Approx(5.5));
	}

	SUBCASE("length == 0 returns NaN") {
		CHECK(std::isnan(smedian({}, 0)));
	}
}

TEST_CASE("mad_median does not modify input") {
	std::vector<double> v{1, 2, 3, 4, 5};
	const auto before = v;

	double mad     = -1.0;
	double median  = -1.0;
	mad_median(v, static_cast<int>(v.size()), mad, median);

	CHECK(v == before);
	CHECK(median == doctest::Approx(3.0));

	// Deviations from median 3: {2, 1, 0, 1, 2}.  After sorting (0,1,1,2,2),
	// ASTAP's median routine averages the three central values for odd
	// lengths >= 5 — (1+1+2)/3 = 4/3.
	CHECK(mad == doctest::Approx(4.0 / 3.0));
}

TEST_CASE("get_best_mean") {
	SUBCASE("empty input returns zeroed result") {
		const auto r = get_best_mean({}, 0);
		CHECK(r.count == 0);
		CHECK(r.mean == 0.0);
		CHECK(r.standard_error_mean == 0.0);
	}

	SUBCASE("single element is returned verbatim with SEM=0") {
		const std::vector<double> v{7.5};
		const auto r = get_best_mean(v, 1);
		CHECK(r.count == 1);
		CHECK(r.mean == doctest::Approx(7.5));
		CHECK(r.standard_error_mean == 0.0);
	}

	SUBCASE("two elements average with SEM=0") {
		const std::vector<double> v{2.0, 4.0};
		const auto r = get_best_mean(v, 2);
		CHECK(r.count == 2);
		CHECK(r.mean == doctest::Approx(3.0));
		CHECK(r.standard_error_mean == 0.0);
	}

	SUBCASE("tight cluster with one extreme outlier rejects the outlier") {
		// Nine samples at ~10, plus one wild 1000.
		const std::vector<double> v{10.0, 10.1,  9.9, 10.2,  9.8,
		                            10.0, 10.1,  9.9, 10.0, 1000.0};
		const auto r = get_best_mean(v, static_cast<int>(v.size()));

		CHECK(r.count == 9);                           // outlier rejected
		CHECK(r.mean == doctest::Approx(10.0).epsilon(0.02));
		CHECK(r.standard_error_mean > 0.0);            // must report non-zero SEM
		CHECK(r.standard_error_mean < 1.0);            // but bounded
	}

	SUBCASE("SEM shrinks as sample size grows (same spread)") {
		// Alternating ±0.5 around 100.0 — use even n so mean is exactly 100.
		auto sample = [](int n) {
			std::vector<double> v;
			v.reserve(n);
			for (int i = 0; i < n; ++i) {
				v.push_back(100.0 + ((i & 1) ? 0.5 : -0.5));
			}
			return v;
		};
		const auto small  = sample(6);
		const auto large  = sample(46);
		const auto rs = get_best_mean(small,  6);
		const auto rl = get_best_mean(large, 46);

		CHECK(rs.mean == doctest::Approx(100.0));
		CHECK(rl.mean == doctest::Approx(100.0));
		// SEM must strictly decrease as count grows (same sigma, larger N).
		CHECK(rl.standard_error_mean < rs.standard_error_mean);
	}

	SUBCASE("all-identical samples: degenerate MAD=0 case rejects everything") {
		// With MAD=0, sigma=0, so the strict `abs(x-median) < 1.5*sigma`
		// test never passes. This matches the Pascal original — callers
		// should not rely on this edge case returning a useful mean.
		const std::vector<double> v{5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
		const auto r = get_best_mean(v, 6);
		CHECK(r.count == 0);
		CHECK(r.mean == 0.0);
	}
}

TEST_CASE("fnmodulo wraps into [0, range)") {
	CHECK(fnmodulo(0.0,   360.0) == doctest::Approx(0.0));
	CHECK(fnmodulo(180.0, 360.0) == doctest::Approx(180.0));
	CHECK(fnmodulo(361.0, 360.0) == doctest::Approx(1.0));
	CHECK(fnmodulo(720.5, 360.0) == doctest::Approx(0.5));

	// Negative input wraps forward.
	CHECK(fnmodulo(-1.5,   360.0) == doctest::Approx(358.5));
	CHECK(fnmodulo(-360.0, 360.0) == doctest::Approx(0.0));
}

///----------------------------------------
/// MARK: Colour — HSV<->RGB / hue / COLORREF decomposition / Kelvin
///----------------------------------------

TEST_CASE("HSV<->RGB round-trips") {
	// Skip achromatic corner (s==0) since hue is undefined there.
	const std::array<std::array<float, 3>, 5> rgbs{{
		{1.0f,  0.0f,  0.0f},
		{0.0f,  1.0f,  0.0f},
		{0.0f,  0.0f,  1.0f},
		{0.5f,  0.25f, 0.75f},
		{0.9f,  0.8f,  0.1f},
	}};

	for (const auto& rgb : rgbs) {
		float h = 0, s = 0, v = 0;
		rgb_to_hsv(rgb[0], rgb[1], rgb[2], h, s, v);

		float r = 0, g = 0, b = 0;
		hsv_to_rgb(h, s, v, r, g, b);

		CHECK(r == doctest::Approx(rgb[0]).epsilon(1e-4));
		CHECK(g == doctest::Approx(rgb[1]).epsilon(1e-4));
		CHECK(b == doctest::Approx(rgb[2]).epsilon(1e-4));
	}
}

TEST_CASE("rgb_to_h primary colours") {
	// Red and 360 are interchangeable; allow either.
	const int red_h = rgb_to_h(1.0f, 0.0f, 0.0f);
	CHECK((red_h == 0 || red_h == 360));

	CHECK(rgb_to_h(0.0f, 1.0f, 0.0f) == 120);
	CHECK(rgb_to_h(0.0f, 0.0f, 1.0f) == 240);
}

TEST_CASE("intensity_rgb decomposes COLORREF") {
	// 0x00BBGGRR — low byte is Red.
	const auto rgb = intensity_rgb(0x00112233u);
	CHECK(rgb[0] == std::uint8_t{0x33}); // R
	CHECK(rgb[1] == std::uint8_t{0x22}); // G
	CHECK(rgb[2] == std::uint8_t{0x11}); // B

	const auto pureRed = intensity_rgb(0x000000FFu);
	CHECK(pureRed[0] == std::uint8_t{0xFF});
	CHECK(pureRed[1] == std::uint8_t{0});
	CHECK(pureRed[2] == std::uint8_t{0});
}

///----------------------------------------
/// MARK: Number / string formatting
///----------------------------------------

TEST_CASE("floattostr* emit fixed decimals with dot separator") {
	SUBCASE("six decimals") {
		const auto s = floattostr6(3.14159265358979);
		CHECK(s == "3.141593");
	}

	SUBCASE("four decimals") {
		CHECK(floattostr4(2.5) == "2.5000");
	}

	SUBCASE("two decimals") {
		CHECK(floattostr2(1.0 / 3.0) == "0.33");
	}

	SUBCASE("dot separator regardless of locale") {
		// Never a comma in the output, even if the host locale uses ",".
		const auto s = floattostr2(0.5);
		CHECK(s.find(',') == std::string::npos);
		CHECK(s.find('.') != std::string::npos);
	}
}

TEST_CASE("inttostr5 right-aligns in width 5") {
	CHECK(inttostr5(0)     == "    0");
	CHECK(inttostr5(42)    == "   42");
	CHECK(inttostr5(12345) == "12345");
}

TEST_CASE("strtoint2 is fault-tolerant") {
	CHECK(strtoint2("42",  99) == 42);
	CHECK(strtoint2("0",   99) == 0);
	CHECK(strtoint2("-7",  99) == -7);
	CHECK(strtoint2("abc", 99) == 99);
	CHECK(strtoint2("",    99) == 99);
}

TEST_CASE("strtofloat2 accepts dot or comma, trims whitespace") {
	CHECK(strtofloat2("3.14")   == doctest::Approx(3.14));
	CHECK(strtofloat2("3,14")   == doctest::Approx(3.14));
	CHECK(strtofloat2(" 2.5 ")  == doctest::Approx(2.5));
	CHECK(strtofloat2("not a number") == doctest::Approx(0.0));
	CHECK(strtofloat2("")       == doctest::Approx(0.0));
}

TEST_CASE("strtofloat1 accepts only the dot separator") {
	CHECK(strtofloat1("3.14") == doctest::Approx(3.14));
	// Comma should fail-parse to 0, not be accepted as a thousands sep.
	CHECK(strtofloat1("3,14") == doctest::Approx(0.0));
	CHECK(strtofloat1("")     == doctest::Approx(0.0));
}

///----------------------------------------
/// MARK: Filename classifiers
///----------------------------------------

TEST_CASE("fits_file_name recognises FITS-family extensions") {
	CHECK(fits_file_name("image.fit")     == true);
	CHECK(fits_file_name("image.fits")    == true);
	CHECK(fits_file_name("image.fts")     == true);
	CHECK(fits_file_name("image.FIT")     == true);    // case-insensitive
	CHECK(fits_file_name("archive.tgz")   == false);
	CHECK(fits_file_name("photo.png")     == false);
	CHECK(fits_file_name("")              == false);
}

TEST_CASE("tiff_file_name recognises TIFF-only extensions") {
	CHECK(tiff_file_name("image.tif")     == true);
	CHECK(tiff_file_name("image.tiff")    == true);
	CHECK(tiff_file_name("image.TIF")     == true);
	CHECK(tiff_file_name("image.fit")     == false);
	CHECK(tiff_file_name("image.png")     == false);
}

TEST_CASE("fits_tiff_file_name recognises both families") {
	CHECK(fits_tiff_file_name("image.fit")  == true);
	CHECK(fits_tiff_file_name("image.tif")  == true);
	CHECK(fits_tiff_file_name("image.jpg")  == false);
}

TEST_CASE("image_file_name recognises common raster formats") {
	CHECK(image_file_name("a.fit")  == true);
	CHECK(image_file_name("a.fits") == true);
	CHECK(image_file_name("a.tif")  == true);
	CHECK(image_file_name("a.png")  == true);
	CHECK(image_file_name("a.jpg")  == true);
	CHECK(image_file_name("a.jpeg") == true);
	CHECK(image_file_name("a.doc")  == false);
	CHECK(image_file_name("")       == false);

	// Matches the faithful port: BMP / PPM / PGM / PFM / XISF are handled
	// by dedicated loaders but are NOT recognised by image_file_name() —
	// the Pascal source gates them separately in load_image().
	CHECK(image_file_name("a.bmp")  == false);
	CHECK(image_file_name("a.xisf") == false);
}

///----------------------------------------
/// MARK: Filename metadata extraction
///----------------------------------------

TEST_CASE("extract_exposure_from_filename") {
	SUBCASE("'<digits>SEC' token") {
		CHECK(extract_exposure_from_filename("M31_300SEC_bin1.fit") == 300);
	}

	SUBCASE("no exposure token returns 0") {
		CHECK(extract_exposure_from_filename("session1.fit") == 0);
		CHECK(extract_exposure_from_filename("")              == 0);
	}
}

TEST_CASE("extract_temperature_from_filename") {
	SUBCASE("no temperature token returns the 999 sentinel") {
		CHECK(extract_temperature_from_filename("session.fit") == 999);
	}
}

///----------------------------------------
/// MARK: Light obfuscation (encrypt / decrypt)
///----------------------------------------

TEST_CASE("encrypt / decrypt round-trip") {
	const std::array<std::string, 5> inputs{
		"",
		"a",
		"password",
		"API-KEY-1234567890",
		"mixed 123 !@# abc",
	};
	for (const auto& s : inputs) {
		CAPTURE(s);
		CHECK(decrypt(encrypt(s)) == s);
	}
}
