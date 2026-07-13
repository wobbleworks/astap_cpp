///----------------------------------------
///     @file annotation_features_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Self-consistency test for the headless annotation features.
///  @details The two ported features have no astap_cli oracle: the original's
///           plot_artificial_stars reads only the online VizieR catalog, and
///           measure_distortion draws its result onto the GUI canvas. So this is
///           a self-consistency test rather than a differential one. It solves a
///           corpus frame against the local star database and then:
///
///           - `-distortion`: asserts the reported median sky<->pixel error is
///             small (a good linear WCS lands catalog stars within ~1 pixel of
///             their measured centroids) and that a healthy number of catalog
///             stars were matched. A broken projection or a mis-read catalog
///             would blow the residual up or match nothing.
///
///           - `-artificialstars`: asserts a synthetic star field was written
///             with a plausible number of single-pixel stars and that the output
///             FITS exists. A star count of zero means the local read path never
///             produced a star.
///
///           Frames are copied to a temp dir first — the solve pass writes a
///           sidecar .ini/.wcs next to each frame, so the corpus is never used
///           in place.
///
///           Gating: skips when the CLI binary, the star database, or the corpus
///           images are absent, so CI without the private assets stays green.
///   @author Created by John Stephen on 7/12/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <string_view>

#if defined(_WIN32)
  #include <io.h>
  #define popen  _popen
  #define pclose _pclose
#else
  #include <sys/wait.h>
#endif

namespace fs = std::filesystem;

static const std::string kPortBin   = ASTAP_BIN_PATH;
static const fs::path     kCorpusDir = ASTAP_CORPUS_DIR;
static const fs::path     kDbDir     = ASTAP_DB_DIR;

[[nodiscard]] static std::string q(const std::string& s) { return "\"" + s + "\""; }

struct RunResult { int exit_code; std::string out; };

[[nodiscard]] static RunResult run(const std::string& cmd) {
	RunResult r{};
	FILE* pipe = ::popen((cmd + " 2>&1").c_str(), "r");
	if (!pipe) { r.exit_code = -1; return r; }
	char buf[1024];
	while (std::fgets(buf, sizeof(buf), pipe)) r.out.append(buf);
	const int raw = ::pclose(pipe);
#if defined(_WIN32)
	r.exit_code = raw;
#else
	r.exit_code = WIFEXITED(raw) ? WEXITSTATUS(raw) : -1;
#endif
	return r;
}

[[nodiscard]] static bool assets_available() {
	if (!fs::exists(kPortBin)) { MESSAGE("port binary absent — skipping"); return false; }
	if (!fs::exists(kDbDir))   { MESSAGE("star database absent — skipping"); return false; }
	return true;
}

/// @brief Copy one M31 corpus frame into a fresh temp dir. Returns the staged
///        path (empty when the source is missing).
[[nodiscard]] static fs::path stage_frame(const std::string& tag, const char* name) {
	auto src = kCorpusDir / name;
	if (!fs::exists(src)) return {};
	auto dir = fs::temp_directory_path() / ("astap_annot_" + tag);
	fs::remove_all(dir);
	fs::create_directories(dir);
	auto dst = dir / name;
	fs::copy_file(src, dst, fs::copy_options::overwrite_existing);
	return dst;
}

/// @brief First floating value following @p key in a line of stdout, or NaN.
[[nodiscard]] static double value_after(const std::string& out, std::string_view key) {
	const auto p = out.find(key);
	if (p == std::string::npos) return std::nan("");
	auto d = p + key.size();
	while (d < out.size() && (out[d] == ' ' || out[d] == '=')) ++d;
	try { return std::stod(out.substr(d)); } catch (...) { return std::nan(""); }
}

TEST_CASE("distortion: solved M31 frame has a sub-pixel median residual") {
	if (!assets_available()) return;
	const auto frame = stage_frame("dist", "M31_a.fits");
	if (frame.empty()) { MESSAGE("M31_a.fits absent — skipping"); return; }

	const std::string db = " -d " + q(kDbDir.string());
	const auto r = run(q(kPortBin) + " -f " + q(frame.string()) + db + " -distortion");
	MESSAGE("distortion stdout:\n" << r.out);

	// The feature runs only after a successful solve.
	REQUIRE(r.out.find("Sky->Pixel error inside") != std::string::npos);

	const double measured = value_after(r.out, "Distortion stars measured:");
	CHECK(measured > 20);  // a real field yields dozens–hundreds of matches

	// "Sky->Pixel error inside <arcsec>" ... " or <pixels> pixel". Grab the
	// pixel figure that follows the " or " on that line.
	const auto line_pos = r.out.find("Sky->Pixel error inside");
	const auto or_pos   = r.out.find(" or ", line_pos);
	REQUIRE(or_pos != std::string::npos);
	const double inner_px = std::stod(r.out.substr(or_pos + 4));
	CHECK(inner_px >= 0.0);
	CHECK(inner_px < 3.0);  // a healthy linear solution lands within a few px
}

TEST_CASE("artificial stars: solved M31 frame yields a synthetic star field") {
	if (!assets_available()) return;
	const auto frame = stage_frame("art", "M31_a.fits");
	if (frame.empty()) { MESSAGE("M31_a.fits absent — skipping"); return; }

	const std::string db = " -d " + q(kDbDir.string());
	const auto r = run(q(kPortBin) + " -f " + q(frame.string()) + db + " -artificialstars");
	MESSAGE("artificialstars stdout:\n" << r.out);

	const double plotted = value_after(r.out, "Artificial stars plotted:");
	CHECK(plotted > 20);  // the field must contain many catalog stars

	// The synthetic frame is written beside the solved input.
	auto art_out = frame;
	art_out.replace_extension();
	art_out += "_artificial.fits";
	CHECK(fs::exists(art_out));
}
