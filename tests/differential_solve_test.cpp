///----------------------------------------
///     @file differential_solve_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Differential plate-solve tests: ASTAP++ vs original astap_cli.
///  @details Drives both binaries through a full solve (`-f img -d db`) over a
///           shared image corpus and star database, then compares the WCS
///           solutions written to the `.ini` files. Rather than diffing the
///           CD matrix element-by-element (two different parameterisations can
///           describe the same mapping), it projects a grid of pixels through
///           each solution and measures the on-sky angular separation — one
///           number that folds in centre, scale, rotation, and skew.
///
///           Plate solving is not bit-deterministic across implementations
///           (slightly different star centroids and quad picks), so agreement
///           is expected at the sub-arcsecond-to-arcsecond level, not exact.
///
///           Gating: skips gracefully when the oracle, the star database, or a
///           corpus image is absent, so CI without the private assets stays
///           green.
///   @author Created by John Stephen on 7/11/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <optional>
#include <string>

#if defined(_WIN32)
  #include <io.h>
  #define popen  _popen
  #define pclose _pclose
#else
  #include <sys/wait.h>
#endif

namespace fs = std::filesystem;

///----------------------------------------
/// MARK: Configuration (injected at configure time)
///----------------------------------------

static const std::string kPortBin   = ASTAP_BIN_PATH;
static const std::string kOracleBin = ASTAP_ORACLE_BIN;
static const fs::path    kCorpusDir = ASTAP_CORPUS_DIR;
static const fs::path    kDbDir     = ASTAP_DB_DIR;

///----------------------------------------
/// MARK: Process + staging helpers
///----------------------------------------

[[nodiscard]] static std::string q(const std::string& s) { return "\"" + s + "\""; }

struct RunResult { int exit_code; std::string stdout_text; };

[[nodiscard]] static RunResult run(const std::string& cmd) {
	RunResult r{};
	FILE* pipe = ::popen((cmd + " 2>&1").c_str(), "r");
	if (!pipe) { r.exit_code = -1; return r; }
	char buf[1024];
	while (std::fgets(buf, sizeof(buf), pipe)) r.stdout_text.append(buf);
	const int raw = ::pclose(pipe);
#if defined(_WIN32)
	r.exit_code = raw;
#else
	r.exit_code = WIFEXITED(raw) ? WEXITSTATUS(raw) : -1;
#endif
	return r;
}

[[nodiscard]] static fs::path stage(const fs::path& src, const std::string& case_name) {
	auto dir = fs::temp_directory_path() / ("astap_solve_" + case_name);
	fs::remove_all(dir);
	fs::create_directories(dir);
	auto dst = dir / src.filename();
	fs::copy_file(src, dst, fs::copy_options::overwrite_existing);
	return dst;
}

[[nodiscard]] static bool assets_available() {
	if (!fs::exists(kOracleBin)) { MESSAGE("oracle absent (" << kOracleBin << ") — skipping"); return false; }
	if (!fs::exists(kDbDir))     { MESSAGE("star database absent (" << kDbDir << ") — skipping"); return false; }
	return true;
}

[[nodiscard]] static std::optional<fs::path> corpus_image(const std::string& name) {
	auto p = kCorpusDir / name;
	if (fs::exists(p)) return p;
	MESSAGE("corpus image absent (" << p.string() << ") — skipping");
	return std::nullopt;
}

///----------------------------------------
/// MARK: WCS parse + projection
///----------------------------------------

struct Wcs {
	bool   solved = false;
	double crpix1{}, crpix2{}, crval1{}, crval2{};
	double cd1_1{}, cd1_2{}, cd2_1{}, cd2_2{};
	double width{}, height{};   // from the DIMENSIONS= line, when present
};

/// @brief Parse the solver .ini (PLTSOLVD + WCS keys). std::stod tolerates the
///        Pascal `str()` layout (leading space, 16 decimals, E+003 exponent).
[[nodiscard]] static Wcs parse_ini(const fs::path& path) {
	Wcs w{};
	std::ifstream ifs(path);
	std::string line;
	auto val = [](const std::string& s) {
		try { return std::stod(s.substr(s.find('=') + 1)); } catch (...) { return 0.0; }
	};
	while (std::getline(ifs, line)) {
		if (line.rfind("PLTSOLVD=", 0) == 0) w.solved = (line.find("=T") != std::string::npos);
		else if (line.rfind("CRPIX1=", 0) == 0) w.crpix1 = val(line);
		else if (line.rfind("CRPIX2=", 0) == 0) w.crpix2 = val(line);
		else if (line.rfind("CRVAL1=", 0) == 0) w.crval1 = val(line);
		else if (line.rfind("CRVAL2=", 0) == 0) w.crval2 = val(line);
		else if (line.rfind("CD1_1=", 0) == 0)  w.cd1_1  = val(line);
		else if (line.rfind("CD1_2=", 0) == 0)  w.cd1_2  = val(line);
		else if (line.rfind("CD2_1=", 0) == 0)  w.cd2_1  = val(line);
		else if (line.rfind("CD2_2=", 0) == 0)  w.cd2_2  = val(line);
		else if (line.rfind("DIMENSIONS=", 0) == 0) {
			// Format: "DIMENSIONS=<w> x <h>".
			const auto rhs = line.substr(line.find('=') + 1);
			const auto xpos = rhs.find('x');
			if (xpos != std::string::npos) {
				try {
					w.width  = std::stod(rhs.substr(0, xpos));
					w.height = std::stod(rhs.substr(xpos + 1));
				} catch (...) {}
			}
		}
	}
	return w;
}

/// @brief Project a pixel through a gnomonic (TAN) WCS to (ra, dec) radians.
static void pixel_to_sky(const Wcs& w, double px, double py, double& ra, double& dec) {
	constexpr double d2r = std::numbers::pi / 180.0;
	// Intermediate world coordinates (degrees → radians).
	const double dx = px - w.crpix1;
	const double dy = py - w.crpix2;
	const double xi  = (w.cd1_1 * dx + w.cd1_2 * dy) * d2r;
	const double eta = (w.cd2_1 * dx + w.cd2_2 * dy) * d2r;
	const double ra0 = w.crval1 * d2r;
	const double dec0 = w.crval2 * d2r;
	// Inverse TAN projection.
	const double r = std::hypot(xi, eta);
	if (r == 0.0) { ra = ra0; dec = dec0; return; }
	const double c = std::atan(r);
	const double sc = std::sin(c), cc = std::cos(c);
	dec = std::asin(cc * std::sin(dec0) + eta * sc * std::cos(dec0) / r);
	ra  = ra0 + std::atan2(xi * sc,
	                       r * std::cos(dec0) * cc - eta * std::sin(dec0) * sc);
}

/// @brief Angular separation (arcsec) between two sky points given in radians.
[[nodiscard]] static double sep_arcsec(double ra1, double dec1, double ra2, double dec2) {
	const double sdd = std::sin((dec2 - dec1) / 2.0);
	const double sda = std::sin((ra2 - ra1) / 2.0);
	const double h = sdd * sdd + std::cos(dec1) * std::cos(dec2) * sda * sda;
	return 2.0 * std::asin(std::min(1.0, std::sqrt(h))) * 180.0 / std::numbers::pi * 3600.0;
}

/// @brief Max angular separation over a pixel grid (corners + centre + edges).
[[nodiscard]] static double max_grid_sep(const Wcs& a, const Wcs& b, double w, double h) {
	double worst = 0.0;
	for (double fy = 0.0; fy <= 1.0; fy += 0.5) {
		for (double fx = 0.0; fx <= 1.0; fx += 0.5) {
			const double px = 1.0 + fx * (w - 1.0);
			const double py = 1.0 + fy * (h - 1.0);
			double ra1, dec1, ra2, dec2;
			pixel_to_sky(a, px, py, ra1, dec1);
			pixel_to_sky(b, px, py, ra2, dec2);
			worst = std::max(worst, sep_arcsec(ra1, dec1, ra2, dec2));
		}
	}
	return worst;
}

/// @brief Pixel scale (arcsec/pixel) from the CD matrix determinant.
[[nodiscard]] static double scale_arcsec(const Wcs& w) {
	return std::sqrt(std::abs(w.cd1_1 * w.cd2_2 - w.cd1_2 * w.cd2_1)) * 3600.0;
}

///----------------------------------------
/// MARK: Case driver
///----------------------------------------

/// @param must_solve when true, both binaries are expected to solve; when
///        false, only PLTSOLVD agreement is asserted (either both solve or
///        both fail).
/// @param extra_args appended verbatim to both command lines (e.g. the
///        `-ra`/`-spd`/`-r` start-position hint). Identical for both binaries,
///        so any divergence is a port-vs-oracle behaviour difference.
/// @param expect_unsolved when true, additionally asserts that neither binary
///        solved — used by the bad-hint cases so the agreement check can't pass
///        vacuously by both spuriously solving.
static void solve_case(const std::string& name, const std::string& file,
                       bool must_solve, const std::string& extra_args = "",
                       bool expect_unsolved = false) {
	if (!assets_available()) return;
	auto img = corpus_image(file);
	if (!img) return;

	const auto o_in = stage(*img, name + "_oracle");
	const auto p_in = stage(*img, name + "_port");
	auto o_ini = o_in; o_ini.replace_extension(".ini");
	auto p_ini = p_in; p_ini.replace_extension(".ini");

	const std::string db = " -d " + q(kDbDir.string());
	const std::string extra = extra_args.empty() ? "" : (" " + extra_args);
	(void)run(q(kOracleBin) + " -f " + q(o_in.string()) + db + extra);
	(void)run(q(kPortBin)   + " -f " + q(p_in.string()) + db + extra);

	REQUIRE_MESSAGE(fs::exists(o_ini), "oracle wrote no .ini at " << o_ini.string());
	REQUIRE_MESSAGE(fs::exists(p_ini), "port wrote no .ini at " << p_ini.string());
	const Wcs o = parse_ini(o_ini);
	const Wcs p = parse_ini(p_ini);

	MESSAGE(name << ": PLTSOLVD oracle=" << o.solved << " port=" << p.solved);
	CHECK(o.solved == p.solved);
	if (expect_unsolved) {
		CHECK_MESSAGE(!o.solved, "oracle unexpectedly solved with a bad hint");
		CHECK_MESSAGE(!p.solved, "port unexpectedly solved with a bad hint");
	}

	if (!o.solved || !p.solved) {
		if (must_solve) { CHECK_MESSAGE(o.solved, "oracle failed to solve " << file);
		                  CHECK_MESSAGE(p.solved, "port failed to solve " << file); }
		return;
	}

	const double centre = sep_arcsec(o.crval1 * std::numbers::pi / 180.0,
	                                 o.crval2 * std::numbers::pi / 180.0,
	                                 p.crval1 * std::numbers::pi / 180.0,
	                                 p.crval2 * std::numbers::pi / 180.0);
	// Both WCS are sampled at the same pixels, so the exact image size only
	// affects which corner of the field is probed. Take it from the .ini.
	const double gw = p.width > 0 ? p.width : (o.width > 0 ? o.width : 1000.0);
	const double gh = p.height > 0 ? p.height : (o.height > 0 ? o.height : 1000.0);
	const double grid = max_grid_sep(o, p, gw, gh);
	const double so = scale_arcsec(o), sp = scale_arcsec(p);
	const double scale_rel = std::abs(so - sp) / so;

	MESSAGE(name << ": centre Δ=" << centre << "\"  max-grid Δ=" << grid
	             << "\"  scale oracle=" << so << " port=" << sp
	             << " (rel " << scale_rel << ")");

	// Cross-implementation solve noise is sub-arcsecond here; 5" leaves margin
	// while still catching a genuinely wrong solution (off by arcminutes+).
	CHECK(centre < 5.0);
	CHECK(grid < 5.0);
	CHECK(scale_rel < 1e-3);
}

///----------------------------------------
/// MARK: Tests
///----------------------------------------

TEST_CASE("solve: M31_a wide field agrees with oracle") {
	solve_case("M31_a", "M31_a.fits", /*must_solve=*/true);
}

TEST_CASE("solve: M31_b agrees with oracle") {
	solve_case("M31_b", "M31_b.fits", /*must_solve=*/true);
}

TEST_CASE("solve: M31_c agrees with oracle") {
	solve_case("M31_c", "M31_c.fits", /*must_solve=*/true);
}

TEST_CASE("solve: ngc_990 sparse field — both decline to solve") {
	// Only ~1 solver star in this 600x600 frame; both must reach the same
	// verdict (PLTSOLVD=F), not one solving while the other fails.
	solve_case("ngc990i", "ngc_990.i.unconv.fits", /*must_solve=*/false);
}

// M31 lies at RA ~00h42.7m (0.712 h) / Dec ~+41.27deg, i.e. south-pole
// distance (spd = 90 + dec) ~131.27deg. The two cases below drive the
// start-position hint (-ra hours, -spd degrees) with a bounded search radius
// and require identical verdicts from both binaries.

TEST_CASE("solve: M31_a with correct -ra/-spd hint agrees with oracle") {
	// A correct hint plus a tight radius must still solve. Guards the -spd
	// units wiring: before the fix the port stored (spd-90) as a raw degree
	// value in a radians field, seeding a nonsense declination that no search
	// radius recovers — the port would fail here while the oracle solved.
	solve_case("M31_a_hint", "M31_a.fits", /*must_solve=*/true,
	           "-ra 0.712 -spd 131.27 -r 10");
}

TEST_CASE("solve: M31_a with wrong -ra hint — both decline to solve") {
	// Correct declination hint but an RA 12h (180deg) away, inside a 5deg
	// search radius. Both binaries must fail identically. Guards the -ra
	// wiring: before the fix the port ignored -ra and kept the correct header
	// RA, so it would solve while the oracle (honouring the bad hint) failed.
	solve_case("M31_a_badra", "M31_a.fits", /*must_solve=*/false,
	           "-ra 12.0 -spd 131.27 -r 5", /*expect_unsolved=*/true);
}
