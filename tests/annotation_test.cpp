///----------------------------------------
///     @file annotation_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for annotation utilities (src/analysis/annotation.cpp).
///  @details Covers the 5x9 bitmap font tables, text-to-image rendering,
///           the rotate() helper, equatorial_standard, and find_object.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright (C) 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "analysis/annotation.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <string>
#include <vector>

using namespace astap::analysis;
using astap::ImageArray;
using astap::Header;

///----------------------------------------
/// MARK: Test-local stubs
///
/// plot_deepsky links core/wcs.cpp for pixel/celestial projection, which
/// references a few externals (dsspos, EQU_GAL, calculate_az_alt_basic,
/// position_angle). Provide no-op stubs so the TAN projection path links in
/// isolation — same set the wcs_test target defines.
///----------------------------------------

namespace astap::core {

void dsspos(double, double, double& ra, double& dec) { ra = 0; dec = 0; }
void EQU_GAL(double, double, double& l, double& b)   { l  = 0; b   = 0; }
bool calculate_az_alt_basic(double, double, double& az, double& alt) {
	az = alt = 0;
	return false;
}
double position_angle(double, double, double, double) { return 0.0; }

}  // namespace astap::core

///----------------------------------------
/// MARK: Header helper
///----------------------------------------

/// @brief Build a 2000x2000 TAN header, 1 arcsec/px, centred on (ra0, dec0).
[[nodiscard]] static Header make_tan_header(double ra0_rad, double dec0_rad) {
	Header h{};
	h.width  = 2000;
	h.height = 2000;
	h.naxis  = 2;
	h.naxis3 = 1;
	h.crpix1 = 1000.0;
	h.crpix2 = 1000.0;
	constexpr double arcsec_deg = 1.0 / 3600.0;
	h.cdelt1 = -arcsec_deg;   // RA decreases with X, by convention
	h.cdelt2 =  arcsec_deg;
	h.ra0    = ra0_rad;
	h.dec0   = dec0_rad;
	h.crota1 = 0.0;
	h.crota2 = 0.0;
	h.cd1_1  = h.cdelt1;
	h.cd1_2  = 0.0;
	h.cd2_1  = 0.0;
	h.cd2_2  = h.cdelt2;
	return h;
}

///----------------------------------------
/// MARK: Font table dimensions
///----------------------------------------

TEST_CASE("kFont5x9: dimensions and value range") {
	// 94 characters (ASCII 33..126), each 9 rows of 5 columns.
	static_assert(sizeof(kFont5x9) == 94 * 9 * 5);

	for (int ch = 0; ch < 94; ++ch) {
		for (int row = 0; row < 9; ++row) {
			for (int col = 0; col < 5; ++col) {
				CHECK_MESSAGE((kFont5x9[ch][row][col] == 0 ||
				               kFont5x9[ch][row][col] == 1),
				              "invalid pixel at [" << ch << "][" << row << "][" << col << "]");
			}
		}
	}
}

///----------------------------------------
/// MARK: Font '!' character
///----------------------------------------

TEST_CASE("kFont5x9: '!' character shape") {
	// kFont5x9[0] = '!' (ASCII 33). Column 2 is set in rows 0-6, blank row 7,
	// dot at row 8 col 2.
	for (int row = 0; row <= 6; ++row) {
		CHECK(kFont5x9[0][row][2] == 1);
	}
	CHECK(kFont5x9[0][7][2] == 0);
	CHECK(kFont5x9[0][8][2] == 1);
}

///----------------------------------------
/// MARK: Font 'A' character
///----------------------------------------

TEST_CASE("kFont5x9: 'A' character row 0") {
	// kFont5x9[32] = 'A' (ASCII 65). Row 0 = {0,0,1,0,0}.
	CHECK(kFont5x9[32][0][0] == 0);
	CHECK(kFont5x9[32][0][1] == 0);
	CHECK(kFont5x9[32][0][2] == 1);
	CHECK(kFont5x9[32][0][3] == 0);
	CHECK(kFont5x9[32][0][4] == 0);
}

///----------------------------------------
/// MARK: annotation_to_array writes pixels
///----------------------------------------

TEST_CASE("annotation_to_array: writes pixels for 'A'") {
	ImageArray img(1);
	img[0].assign(50, std::vector<float>(50, 0.0f));

	annotation_to_array("A", false, 1000, 1, 5, 20, img);

	// Some pixels at the expected location should be non-zero.
	bool found_nonzero = false;
	for (int r = 10; r <= 25; ++r) {
		for (int c = 5; c < 12; ++c) {
			if (img[0][r][c] != 0.0f) {
				found_nonzero = true;
				break;
			}
		}
		if (found_nonzero) break;
	}
	CHECK(found_nonzero);
}

///----------------------------------------
/// MARK: annotation_to_array transparent mode
///----------------------------------------

TEST_CASE("annotation_to_array: transparent mode preserves background") {
	ImageArray img(1);
	img[0].assign(50, std::vector<float>(50, 42.0f));

	annotation_to_array("A", true, 1000, 1, 5, 20, img);

	// In transparent mode, pixels where font=0 should retain 42.0.
	// Find at least one background pixel in the glyph bounding box.
	bool found_preserved = false;
	for (int r = 12; r <= 20; ++r) {
		for (int c = 5; c < 10; ++c) {
			if (img[0][r][c] == 42.0f) {
				found_preserved = true;
				break;
			}
		}
		if (found_preserved) break;
	}
	CHECK(found_preserved);
}

///----------------------------------------
/// MARK: rotate identity
///----------------------------------------

TEST_CASE("rotate: identity at rot=0") {
	auto [x2, y2] = rotate(0.0, 1.0, 0.0);
	// x2 = x*sin(0) + y*cos(0) = 0, y2 = -x*cos(0) + y*sin(0) = -1
	CHECK(x2 == doctest::Approx(0.0).epsilon(1e-12));
	CHECK(y2 == doctest::Approx(-1.0).epsilon(1e-12));
}

TEST_CASE("rotate: rot=0 with (0, 1)") {
	auto [x2, y2] = rotate(0.0, 0.0, 1.0);
	// x2 = 0*sin(0) + 1*cos(0) = 1, y2 = -0*cos(0) + 1*sin(0) = 0
	CHECK(x2 == doctest::Approx(1.0).epsilon(1e-12));
	CHECK(y2 == doctest::Approx(0.0).epsilon(1e-12));
}

///----------------------------------------
/// MARK: equatorial_standard at reference point
///----------------------------------------

TEST_CASE("equatorial_standard: at reference point returns (0, 0)") {
	double ra0 = 1.5;
	double dec0 = 0.8;
	auto [xx, yy] = equatorial_standard(ra0, dec0, ra0, dec0, 1.0);
	CHECK(xx == doctest::Approx(0.0).epsilon(1e-8));
	CHECK(yy == doctest::Approx(0.0).epsilon(1e-8));
}

///----------------------------------------
/// MARK: find_object with mock database
///----------------------------------------

TEST_CASE("find_object: finds object in mock database") {
	// Minimal CSV database matching the expected format.
	// Lines 0 and 1 are headers (skipped). Line 2+ are data.
	// Fields: RA_encoded, Dec_encoded, names, length, width, PA
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"432000,324000, M31/NGC224,178.0,63.0,35"
	};

	auto result = find_object("M31", db);
	REQUIRE(result.has_value());
	CHECK(result->name == "M31_NGC224");
	CHECK(result->length == doctest::Approx(178.0));
	CHECK(result->width == doctest::Approx(63.0));
	CHECK(result->pa == doctest::Approx(35.0));
}

TEST_CASE("find_object: returns nullopt for missing object") {
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"432000,324000, M31/NGC224,178.0,63.0,35"
	};

	auto result = find_object("M42", db);
	CHECK_FALSE(result.has_value());
}

///----------------------------------------
/// MARK: read_deepsky field parsing
///----------------------------------------

TEST_CASE("read_deepsky: parses a record inside the field") {
	// RA 432000 -> pi rad; Dec 0 -> 0 rad.
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"432000,0, M31/NGC224,178.0,63.0,35"
	};

	auto objs = read_deepsky(db, DeepSkySearch::WithinField,
	                         std::numbers::pi, 0.0, 0.1);
	REQUIRE(objs.size() == 1);
	CHECK(objs[0].name == "M31");
	CHECK(objs[0].name2 == "NGC224");
	CHECK(objs[0].name3.empty());
	CHECK(objs[0].ra == doctest::Approx(std::numbers::pi).epsilon(1e-9));
	CHECK(objs[0].dec == doctest::Approx(0.0).epsilon(1e-12));
	CHECK(objs[0].length == doctest::Approx(178.0));
	CHECK(objs[0].width == doctest::Approx(63.0));
	CHECK(objs[0].pa == doctest::Approx(35.0));
}

///----------------------------------------
/// MARK: read_deepsky FOV gate
///----------------------------------------

TEST_CASE("read_deepsky: FOV gate excludes distant objects") {
	// First object at the field centre (pi, 0); second far away (0, pi/2).
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"432000,0, M31,178.0,63.0,35",
		"0,324000, FAR,10.0,10.0,0"
	};

	auto objs = read_deepsky(db, DeepSkySearch::WithinField,
	                         std::numbers::pi, 0.0, 0.1);
	REQUIRE(objs.size() == 1);
	CHECK(objs[0].name == "M31");
}

TEST_CASE("read_deepsky: FullDatabase ignores the FOV gate") {
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"432000,0, M31,178.0,63.0,35",
		"0,324000, FAR,10.0,10.0,0"
	};

	auto objs = read_deepsky(db, DeepSkySearch::FullDatabase,
	                         std::numbers::pi, 0.0, 0.1);
	CHECK(objs.size() == 2);
}

///----------------------------------------
/// MARK: read_deepsky RA wraparound
///----------------------------------------

TEST_CASE("read_deepsky: RA difference wraps across 0/2pi") {
	// RA 863999 -> ~2pi, i.e. ~7.3e-6 rad short of a full turn.  A field
	// centred at RA 0 must include it via the wrapped delta_ra fold.
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"863999,0, WRAP,10.0,5.0,0"
	};

	auto objs = read_deepsky(db, DeepSkySearch::WithinField,
	                         0.0, 0.0, 0.001);
	REQUIRE(objs.size() == 1);
	CHECK(objs[0].name == "WRAP");
}

///----------------------------------------
/// MARK: read_deepsky unknown position angle
///----------------------------------------

TEST_CASE("read_deepsky: unparseable PA maps to 999") {
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"432000,0, NOPA,10.0,5.0,"
	};

	auto objs = read_deepsky(db, DeepSkySearch::WithinField,
	                         std::numbers::pi, 0.0, 0.1);
	REQUIRE(objs.size() == 1);
	CHECK(objs[0].pa == doctest::Approx(999.0));
}

///----------------------------------------
/// MARK: read_deepsky skips header lines
///----------------------------------------

TEST_CASE("read_deepsky: skips the two header lines") {
	// The first two lines look like valid records but must be treated as
	// headers and skipped; only the third record is returned.
	std::vector<std::string> db = {
		"432000,0,HDR1,1,1,1",
		"432000,0,HDR2,1,1,1",
		"432000,0, REAL,10.0,5.0,35"
	};

	auto objs = read_deepsky(db, DeepSkySearch::WithinField,
	                         std::numbers::pi, 0.0, 0.1);
	REQUIRE(objs.size() == 1);
	CHECK(objs[0].name == "REAL");
}

///----------------------------------------
/// MARK: plot_deepsky projects an object at the field centre
///----------------------------------------

TEST_CASE("plot_deepsky: places a catalog object at the image centre") {
	const auto head = make_tan_header(1.0, 0.3);

	// One object encoded at (ra0, dec0) = (1.0, 0.3). RA unit = 2pi/864000,
	// Dec unit = (pi/2)/324000. length 60 (0.1'), width 30, PA 45.
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"137510,61879, NGC7000,60.0,30.0,45"
	};

	auto plotted = plot_deepsky(head, db);
	REQUIRE(plotted.size() == 1);

	const auto& p = plotted[0];
	CHECK(p.name == "NGC7000");
	CHECK(p.has_label);
	// Marker centre near the image centre; default vertical flip puts y at
	// (height-1) - ~999 = ~1000.
	CHECK(std::abs(p.x - 999) <= 2);
	CHECK(std::abs(p.y - 1000) <= 2);
	// len = 60 / (|cdelt2| * 60 * 10 * 2) = 60 / ((1/3600)*1200) = 180 px.
	CHECK(p.radius == doctest::Approx(180.0).epsilon(1e-3));
	CHECK(p.axis_ratio == doctest::Approx(0.5));
	CHECK(p.marker == DeepSkyMarker::Galaxy);
}

///----------------------------------------
/// MARK: plot_deepsky FOV exclusion
///----------------------------------------

TEST_CASE("plot_deepsky: excludes an object outside the field") {
	const auto head = make_tan_header(1.0, 0.3);

	// Same RA, but Dec offset by ~30 degrees — far outside the ~0.6 deg field.
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"137510,169879, FAR,60.0,30.0,45"
	};

	auto plotted = plot_deepsky(head, db);
	CHECK(plotted.empty());
}

///----------------------------------------
/// MARK: plot_deepsky requires a solution
///----------------------------------------

TEST_CASE("plot_deepsky: returns nothing without an astrometric solution") {
	Header head{};  // naxis == 0, cd1_1 == 0
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"137510,61879, NGC7000,60.0,30.0,45"
	};

	auto plotted = plot_deepsky(head, db);
	CHECK(plotted.empty());
}

///----------------------------------------
/// MARK: layout_labels single label unchanged
///----------------------------------------

TEST_CASE("layout_labels: a lone label keeps its anchor") {
	std::vector<LabelInput> labels = {{100, 100, 50, 20}};
	auto boxes = layout_labels(labels, 2000, 2000);

	REQUIRE(boxes.size() == 1);
	CHECK(boxes[0].x1 == 100);
	CHECK(boxes[0].y1 == 100);
	CHECK(boxes[0].x2 == 150);
	CHECK(boxes[0].y2 == 120);
	CHECK_FALSE(boxes[0].connector);
}

///----------------------------------------
/// MARK: layout_labels overlap shifts down + connector
///----------------------------------------

TEST_CASE("layout_labels: an overlapping label shifts down and gets a connector") {
	// Two identical anchors: the second must move clear of the first.
	std::vector<LabelInput> labels = {{100, 100, 50, 20}, {100, 100, 50, 20}};
	auto boxes = layout_labels(labels, 2000, 2000);

	REQUIRE(boxes.size() == 2);
	// First stays put.
	CHECK(boxes[0].y1 == 100);
	CHECK_FALSE(boxes[0].connector);
	// Second shifts down by th/3 (=6) until clear: 100->124.
	CHECK(boxes[1].y1 == 124);
	CHECK(boxes[1].connector);
}

///----------------------------------------
/// MARK: layout_labels right-edge clamp
///----------------------------------------

TEST_CASE("layout_labels: a label past the right edge shifts left") {
	std::vector<LabelInput> labels = {{1990, 100, 50, 20}};
	auto boxes = layout_labels(labels, 2000, 2000);

	REQUIRE(boxes.size() == 1);
	CHECK(boxes[0].x2 == 2000);
	CHECK(boxes[0].x1 == 1950);
	CHECK_FALSE(boxes[0].connector);
}

///----------------------------------------
/// MARK: layout_labels falls back to anchor with no vertical space
///----------------------------------------

TEST_CASE("layout_labels: no vertical space falls back to the anchor") {
	// Short image: shifting down would cross the bottom, so the second label
	// stays at its anchor (accepting overlap) with no connector.
	std::vector<LabelInput> labels = {{100, 100, 50, 20}, {100, 100, 50, 20}};
	auto boxes = layout_labels(labels, 2000, 130);

	REQUIRE(boxes.size() == 2);
	CHECK(boxes[1].y1 == 100);
	CHECK_FALSE(boxes[1].connector);
}

///----------------------------------------
/// MARK: DeepSkyDatabase catalog loading
///----------------------------------------

namespace {

/// Write @p contents to @p path, creating parent directories as needed.
static void write_file(const std::filesystem::path& path, std::string_view contents) {
	std::filesystem::create_directories(path.parent_path());
	std::ofstream out(path, std::ios::binary | std::ios::trunc);
	out.write(contents.data(), static_cast<std::streamsize>(contents.size()));
}

/// A unique scratch directory for one test case, cleared on entry.
static std::filesystem::path scratch_dir(std::string_view leaf) {
	auto dir = std::filesystem::temp_directory_path() / "astappp_deepsky_test" / leaf;
	std::error_code ec;
	std::filesystem::remove_all(dir, ec);
	std::filesystem::create_directories(dir);
	return dir;
}

}  // anonymous namespace

TEST_CASE("DeepSkyDatabase: loads a catalog and exposes its lines") {
	const auto dir = scratch_dir("load");
	write_file(dir / "deep_sky.csv",
	           "header line 1\nheader line 2\n432000,0, M31/NGC224,178.0,63.0,35\n");

	DeepSkyDatabase db;
	CHECK(db.load(DeepSkyCatalog::DeepSky, dir) == CatalogLoadStatus::Loaded);
	REQUIRE(db.resident().has_value());
	CHECK(db.resident().value() == DeepSkyCatalog::DeepSky);
	REQUIRE(db.lines().size() == 3);

	// The resident lines feed read_deepsky directly.
	auto objs = read_deepsky(db.lines(), DeepSkySearch::WithinField,
	                         std::numbers::pi, 0.0, 0.1);
	REQUIRE(objs.size() == 1);
	CHECK(objs[0].name == "M31");
}

TEST_CASE("DeepSkyDatabase: re-requesting the resident catalog is a no-op") {
	const auto dir = scratch_dir("cache");
	write_file(dir / "deep_sky.csv",
	           "header line 1\nheader line 2\n432000,0, M31,178.0,63.0,35\n");

	DeepSkyDatabase db;
	CHECK(db.load(DeepSkyCatalog::DeepSky, dir) == CatalogLoadStatus::Loaded);
	CHECK(db.load(DeepSkyCatalog::DeepSky, dir) == CatalogLoadStatus::AlreadyLoaded);
	CHECK(db.lines().size() == 3);
}

TEST_CASE("DeepSkyDatabase: switching catalogs reloads") {
	const auto dir = scratch_dir("switch");
	write_file(dir / "deep_sky.csv", "h1\nh2\n432000,0, M31,178.0,63.0,35\n");
	write_file(dir / "hyperleda.csv", "h1\nh2\n432000,0, PGC1,10.0,5.0,0\n432000,0, PGC2,9.0,4.0,0\n");

	DeepSkyDatabase db;
	CHECK(db.load(DeepSkyCatalog::DeepSky, dir) == CatalogLoadStatus::Loaded);
	CHECK(db.load(DeepSkyCatalog::HyperLeda, dir) == CatalogLoadStatus::Loaded);
	CHECK(db.resident().value() == DeepSkyCatalog::HyperLeda);
	CHECK(db.lines().size() == 4);
}

TEST_CASE("DeepSkyDatabase: a missing file reports NotFound and empties the database") {
	const auto dir = scratch_dir("missing");
	write_file(dir / "deep_sky.csv", "h1\nh2\n432000,0, M31,178.0,63.0,35\n");

	DeepSkyDatabase db;
	REQUIRE(db.load(DeepSkyCatalog::DeepSky, dir) == CatalogLoadStatus::Loaded);

	// variable_stars.csv was never written.
	CHECK(db.load(DeepSkyCatalog::Variable, dir) == CatalogLoadStatus::NotFound);
	CHECK_FALSE(db.resident().has_value());
	CHECK(db.lines().empty());
}

TEST_CASE("DeepSkyDatabase: a versioned catalog with the right header loads clean") {
	const auto dir = scratch_dir("v003_ok");
	write_file(dir / "variable_stars_13.csv",
	           "V003 variable stars\nheader 2\n432000,0, 000-BBB-001,5.0,5.0,0\n");

	DeepSkyDatabase db;
	CHECK(db.load(DeepSkyCatalog::Variable13, dir) == CatalogLoadStatus::Loaded);
	CHECK(db.resident().value() == DeepSkyCatalog::Variable13);
}

TEST_CASE("DeepSkyDatabase: a versioned catalog with an old header loads but reports Outdated") {
	const auto dir = scratch_dir("v003_old");
	write_file(dir / "variable_stars_15.csv",
	           "V001 old variable stars\nheader 2\n432000,0, 000-CCC-001,5.0,5.0,0\n");

	DeepSkyDatabase db;
	// Data still resident (Pascal keeps it and only warns), status flags the age.
	CHECK(db.load(DeepSkyCatalog::Variable15, dir) == CatalogLoadStatus::Outdated);
	REQUIRE(db.resident().has_value());
	CHECK(db.resident().value() == DeepSkyCatalog::Variable15);
	CHECK(db.lines().size() == 3);
}

///----------------------------------------
/// MARK: extract_visible photometry targets
///----------------------------------------

TEST_CASE("extract_visible: keeps only on-image, named objects, in order") {
	std::vector<PlottedDeepSky> plotted(3);
	// Labeled target.
	plotted[0].ra = 1.0; plotted[0].dec = 0.3;
	plotted[0].abbr = "NGC7000"; plotted[0].has_label = true;
	// Off-image (no label) — excluded.
	plotted[1].ra = 1.1; plotted[1].dec = 0.31;
	plotted[1].abbr = "OFF"; plotted[1].has_label = false;
	// Second labeled target.
	plotted[2].ra = 1.2; plotted[2].dec = 0.32;
	plotted[2].abbr = "IC5070"; plotted[2].has_label = true;

	auto targets = extract_visible(plotted);
	REQUIRE(targets.size() == 2);

	CHECK(targets[0].abbr == "NGC7000");
	CHECK(targets[0].ra == doctest::Approx(1.0));
	CHECK(targets[0].dec == doctest::Approx(0.3));
	CHECK(targets[0].source == 0);
	CHECK(targets[0].index == 0);

	CHECK(targets[1].abbr == "IC5070");
	CHECK(targets[1].index == 1);
}

TEST_CASE("extract_visible: honours the target cap") {
	std::vector<PlottedDeepSky> plotted(5);
	for (auto& p : plotted) {
		p.abbr = "X";
		p.has_label = true;
	}

	auto targets = extract_visible(plotted, 0, 3);
	CHECK(targets.size() == 3);
}

TEST_CASE("extract_visible: source tag is stored on every target") {
	std::vector<PlottedDeepSky> plotted(1);
	plotted[0].abbr = "V001";
	plotted[0].has_label = true;

	auto targets = extract_visible(plotted, 1);  // 1 = VSX, e.g.
	REQUIRE(targets.size() == 1);
	CHECK(targets[0].source == 1);
}

TEST_CASE("extract_visible: abbr is the primary designation, not the composed label") {
	const auto head = make_tan_header(1.0, 0.3);

	// A three-name object: composed label is naam2/naam3/naam4, but the
	// photometry abbr must be naam2 alone.
	std::vector<std::string> db = {
		"header line 1",
		"header line 2",
		"137510,61879, NGC7000/C20/LBN373,60.0,30.0,45"
	};

	auto plotted = plot_deepsky(head, db);
	REQUIRE(plotted.size() == 1);
	CHECK(plotted[0].name == "NGC7000/C20/LBN373");
	CHECK(plotted[0].abbr == "NGC7000");

	auto targets = extract_visible(plotted);
	REQUIRE(targets.size() == 1);
	CHECK(targets[0].abbr == "NGC7000");
}
