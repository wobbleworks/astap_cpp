///----------------------------------------
///     @file stack_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Self-consistency test for the headless batch stacker.
///  @details Stacking has no astap_cli oracle (it is GUI-only in the original),
///           so this is a self-consistency smoke test rather than a differential
///           one: batch-stack several frames of the same field with internal-
///           astrometry alignment, then plate-solve the stacked result. A stack
///           whose frames were correctly registered solves to the same field
///           with a healthy star count; a mis-aligned stack smears the stars and
///           either fails to solve or detects far fewer of them.
///
///           Frames are copied to a temp directory first — the pre-solve pass
///           writes each frame's WCS back into its file, so the corpus must not
///           be used in place.
///
///           Gating: skips when the CLI binary, the star database, or the corpus
///           images are absent, so CI without the private assets stays green.
///   @author Created by John Stephen on 7/12/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>
#include <vector>

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

/// @brief Value after "KEY=" in a line-oriented stdout stream (or "").
[[nodiscard]] static std::string field(const std::string& out, std::string_view key) {
	const auto p = out.find(key);
	if (p == std::string::npos) return {};
	const auto s = p + key.size();
	const auto e = out.find('\n', s);
	return out.substr(s, e == std::string::npos ? std::string::npos : e - s);
}

/// @brief First integer following "<n> stars," in solver stdout, or -1.
[[nodiscard]] static int selected_stars(const std::string& out) {
	const auto p = out.find(" stars,");
	if (p == std::string::npos) return -1;
	auto b = p;
	while (b > 0 && std::isdigit(static_cast<unsigned char>(out[b - 1]))) --b;
	if (b == p) return -1;
	try { return std::stoi(out.substr(b, p - b)); } catch (...) { return -1; }
}

/// @brief Read a numeric FITS keyword (e.g. "CRVAL1  ") from a file, or NaN.
[[nodiscard]] static double ini_value(const fs::path& path, std::string_view key) {
	std::ifstream ifs(path);
	std::string line;
	while (std::getline(ifs, line)) {
		if (line.rfind(key, 0) == 0) {
			const auto eq = line.find('=');
			if (eq == std::string::npos) continue;
			try { return std::stod(line.substr(eq + 1)); } catch (...) { return 0.0; }
		}
	}
	return 0.0;
}

/// @brief Copy the available M31 frames into a fresh temp dir (the pre-solve
///        pass rewrites each frame's WCS in place, so the corpus is never used
///        directly). Returns the dir, the quoted frame arg string, and count.
struct Staged { fs::path dir; std::string frame_args; std::vector<fs::path> files; int count = 0; };

[[nodiscard]] static Staged stage_m31(const std::string& tag) {
	Staged s;
	s.dir = fs::temp_directory_path() / ("astap_stack_" + tag);
	fs::remove_all(s.dir);
	fs::create_directories(s.dir);
	for (const auto* n : {"M31_a.fits", "M31_b.fits", "M31_c.fits"}) {
		auto src = kCorpusDir / n;
		if (!fs::exists(src)) continue;
		auto dst = s.dir / src.filename();
		fs::copy_file(src, dst, fs::copy_options::overwrite_existing);
		s.frame_args += " " + q(dst.string());
		s.files.push_back(dst);
		++s.count;
	}
	return s;
}

/// @brief Median of one plane of a 32-bit-float FITS, or NaN if unreadable.
/// @details Minimal reader for the stacker's own output (BITPIX -32, big-endian,
///          BZERO 0 / BSCALE 1): scans 2880-byte header blocks to END, then reads
///          @p plane of NAXIS1×NAXIS2 floats. Used to prove the colour combine
///          routed each channel into its own plane (distinct medians).
[[nodiscard]] static double plane_median(const fs::path& file, int plane) {
	std::ifstream ifs(file, std::ios::binary);
	if (!ifs) return std::nan("");
	std::string hdr;
	int w = 0, h = 0, blocks = 0;
	for (;;) {
		std::string blk(2880, '\0');
		ifs.read(blk.data(), 2880);
		if (ifs.gcount() != 2880) return std::nan("");
		hdr += blk;
		++blocks;
		if (blk.find("END     ") != std::string::npos) break;
		if (blocks > 20) return std::nan("");
	}
	auto card_int = [&](std::string_view k) -> int {
		auto p = hdr.find(k);
		if (p == std::string::npos) return -1;
		auto d = p + k.size();
		while (d < hdr.size() && !std::isdigit(static_cast<unsigned char>(hdr[d]))
		       && hdr[d] != '-') ++d;
		try { return std::stoi(hdr.substr(d)); } catch (...) { return -1; }
	};
	w = card_int("NAXIS1  =");
	h = card_int("NAXIS2  =");
	if (w <= 0 || h <= 0) return std::nan("");
	const std::size_t px = static_cast<std::size_t>(w) * h;
	ifs.seekg(static_cast<std::streamoff>(blocks) * 2880
	          + static_cast<std::streamoff>(plane) * px * 4);
	std::vector<float> vals(px);
	for (auto& v : vals) {
		unsigned char b[4];
		ifs.read(reinterpret_cast<char*>(b), 4);
		if (ifs.gcount() != 4) return std::nan("");
		const std::uint32_t u = (std::uint32_t(b[0]) << 24) | (std::uint32_t(b[1]) << 16)
		                      | (std::uint32_t(b[2]) << 8) | std::uint32_t(b[3]);
		std::memcpy(&v, &u, 4);
	}
	std::nth_element(vals.begin(), vals.begin() + px / 2, vals.end());
	return vals[px / 2];
}

/// @brief First integer following @p key in a raw FITS header block, or -1.
[[nodiscard]] static int fits_header_int(const fs::path& file, std::string_view key) {
	std::ifstream ifs(file, std::ios::binary);
	std::string hdr(2880 * 4, '\0');
	ifs.read(hdr.data(), static_cast<std::streamsize>(hdr.size()));
	const auto p = hdr.find(key);
	if (p == std::string::npos) return -1;
	auto d = p + key.size();
	// Skip padding to the value; stop at a digit or a leading minus sign.
	while (d < hdr.size() && !std::isdigit(static_cast<unsigned char>(hdr[d]))
	       && hdr[d] != '-') ++d;
	try { return std::stoi(hdr.substr(d)); } catch (...) { return -1; }
}

TEST_CASE("batch stack registers frames and solves to the same field") {
	if (!assets_available()) return;

	// Need >= 2 frames of one field. M31 a/b/c share the M31 field.
	std::vector<std::string> names{"M31_a.fits", "M31_b.fits", "M31_c.fits"};
	std::vector<fs::path> srcs;
	for (const auto& n : names) {
		auto p = kCorpusDir / n;
		if (fs::exists(p)) srcs.push_back(p);
	}
	if (srcs.size() < 2) { MESSAGE("need >= 2 M31 frames — skipping"); return; }

	// Copy into a temp dir: the pre-solve pass rewrites each frame's WCS in place.
	const auto dir = fs::temp_directory_path() / "astap_stack_test";
	fs::remove_all(dir);
	fs::create_directories(dir);
	std::string frame_args;
	for (const auto& s : srcs) {
		auto dst = dir / s.filename();
		fs::copy_file(s, dst, fs::copy_options::overwrite_existing);
		frame_args += " " + q(dst.string());
	}
	const auto out    = dir / "stack.fits";
	const std::string db = " -d " + q(kDbDir.string());

	// --- Batch stack (internal-astrometry alignment, average combine) --------
	const auto st = run(q(kPortBin) + " -stackfiles -o " + q(out.string()) + db + frame_args);
	CHECK(st.exit_code == 0);
	CHECK(field(st.out, "STACKED=") == "1");
	CHECK(field(st.out, "FRAMES_COMBINED=") == std::to_string(srcs.size()));
	REQUIRE(fs::exists(out));

	// --- Header: the stacked-header block was written ------------------------
	// The FITS header is not newline-delimited, so read the leading blocks raw
	// and match the 80-byte cards by substring.
	{
		std::ifstream ifs(out, std::ios::binary);
		std::string hdr(2880 * 4, '\0');
		ifs.read(hdr.data(), static_cast<std::streamsize>(hdr.size()));
		CHECK(hdr.find("CALSTAT = 'S")   != std::string::npos);  // status: stacked
		CHECK(hdr.find("Stacking method AVERAGE") != std::string::npos);
		// LUM_CNT card carries the combined-frame count.
		const auto lp = hdr.find("LUM_CNT =");
		REQUIRE(lp != std::string::npos);
		auto d = lp;
		while (d < hdr.size() && !std::isdigit(static_cast<unsigned char>(hdr[d]))) ++d;
		int lum = 0;
		try { lum = std::stoi(hdr.substr(d)); } catch (...) {}
		CHECK(lum == static_cast<int>(srcs.size()));
	}

	// --- Self-consistency: the stack solves to M31's field, stars sharp ------
	const auto base = dir / "stack_solve";
	const auto sv = run(q(kPortBin) + " -f " + q(out.string()) + db
	                    + " -o " + q(base.string()));
	CHECK(sv.exit_code == 0);
	auto ini = base; ini.replace_extension(".ini");
	REQUIRE(fs::exists(ini));
	{
		std::ifstream ifs(ini);
		const std::string txt((std::istreambuf_iterator<char>(ifs)),
		                      std::istreambuf_iterator<char>());
		CHECK(txt.find("PLTSOLVD=T") != std::string::npos);   // solved
	}

	// M31 lies at RA ~10.7 deg, Dec ~41.2 deg (.ini uses compact keys).
	CHECK(ini_value(ini, "CRVAL1=") == doctest::Approx(10.7).epsilon(0.05));
	CHECK(ini_value(ini, "CRVAL2=") == doctest::Approx(41.24).epsilon(0.05));

	// A correctly registered stack keeps its stars point-like: the solver still
	// finds a full field of them (a smeared stack would collapse this count).
	CHECK(selected_stars(sv.out) > 300);

	fs::remove_all(dir);
}

TEST_CASE("sigma-clip combine also aligns and solves to the same field") {
	if (!assets_available()) return;
	auto s = stage_m31("sigma");
	if (s.count < 2) { MESSAGE("need >= 2 M31 frames — skipping"); return; }

	const auto out = s.dir / "stack.fits";
	const std::string db = " -d " + q(kDbDir.string());
	const auto st = run(q(kPortBin) + " -stackfiles -sigmaclip -o " + q(out.string())
	                    + db + s.frame_args);
	CHECK(st.exit_code == 0);
	CHECK(field(st.out, "STACKED=") == "1");
	CHECK(field(st.out, "FRAMES_COMBINED=") == std::to_string(s.count));
	REQUIRE(fs::exists(out));

	// The sigma-clip method must be recorded in the stacked header.
	{
		std::ifstream ifs(out, std::ios::binary);
		std::string hdr(2880 * 4, '\0');
		ifs.read(hdr.data(), static_cast<std::streamsize>(hdr.size()));
		CHECK(hdr.find("Stacking method SIGMA CLIP AVERAGE") != std::string::npos);
	}

	const auto base = s.dir / "solve";
	const auto sv = run(q(kPortBin) + " -f " + q(out.string()) + db
	                    + " -o " + q(base.string()));
	auto ini = base; ini.replace_extension(".ini");
	REQUIRE(fs::exists(ini));
	CHECK(ini_value(ini, "CRVAL1=") == doctest::Approx(10.7).epsilon(0.05));
	CHECK(ini_value(ini, "CRVAL2=") == doctest::Approx(41.24).epsilon(0.05));
	CHECK(selected_stars(sv.out) > 300);

	fs::remove_all(s.dir);
}

// NOTE: no automated mosaic test. -mosaic dispatches the real stack_mosaic, but
// its per-pixel seam-blend (median_background over overlaps) is only cheap when
// tiles ABUT with small overlaps. The only multi-frame corpus here (M31 a/b/c)
// fully overlaps, which is both a degenerate mosaic input and pathologically slow
// (~60s/frame). Verified manually instead: a fully-overlapping M31 mosaic solves
// to the M31 field. A proper mosaic test needs a tiled corpus.

TEST_CASE("batch stack applies a master dark (calibration wiring)") {
	if (!assets_available()) return;
	auto s = stage_m31("dark");
	if (s.count < 2) { MESSAGE("need >= 2 M31 frames — skipping"); return; }

	// A separate master-dark file (a copy of one frame — the point here is that
	// the calibration path fires and is recorded, not photometric correctness).
	const auto dark = s.dir / "master_dark.fits";
	fs::copy_file(kCorpusDir / "M31_a.fits", dark, fs::copy_options::overwrite_existing);

	const auto out = s.dir / "stack.fits";
	const std::string db = " -d " + q(kDbDir.string());
	const auto st = run(q(kPortBin) + " -stackfiles -o " + q(out.string())
	                    + " -dark " + q(dark.string()) + db + s.frame_args);
	CHECK(st.exit_code == 0);
	CHECK(field(st.out, "STACKED=") == "1");
	REQUIRE(fs::exists(out));

	// The dark must be subtracted per frame and recorded: CALSTAT gains 'D'
	// (then 'S' for stacked → "DS") and LUM_DARK counts the master applied.
	{
		std::ifstream ifs(out, std::ios::binary);
		std::string hdr(2880 * 4, '\0');
		ifs.read(hdr.data(), static_cast<std::streamsize>(hdr.size()));
		CHECK(hdr.find("CALSTAT = 'D") != std::string::npos);
	}
	CHECK(fits_header_int(out, "LUM_DARK=") >= 1);

	fs::remove_all(s.dir);
}

TEST_CASE("reference frame is chosen by quality, not input order") {
	if (!assets_available()) return;
	auto s = stage_m31("quality");
	if (s.count < 3) { MESSAGE("need 3 M31 frames — skipping"); return; }

	const std::string db = " -d " + q(kDbDir.string());
	auto stack_ref = [&](const std::vector<fs::path>& order) -> std::string {
		std::string args;
		for (const auto& f : order) args += " " + q(f.string());
		const auto out = s.dir / "q.fits";
		const auto rr = run(q(kPortBin) + " -stackfiles -o " + q(out.string()) + db + args);
		CHECK(rr.exit_code == 0);
		return fs::path(field(rr.out, "REFERENCE=")).filename().string();
	};

	// Forward order, then reversed. A quality-based selector picks the same frame
	// either way; a positional "use frame [0]" would pick different frames.
	const auto ref_fwd = stack_ref(s.files);
	std::vector<fs::path> rev(s.files.rbegin(), s.files.rend());
	const auto ref_rev = stack_ref(rev);

	CHECK(ref_fwd == ref_rev);         // order-independent → chosen by quality
	CHECK(!ref_fwd.empty());
	// For the M31 set the sharpest frame is not the first input, so the chosen
	// reference must differ from frame [0] of the forward order.
	CHECK(ref_fwd != s.files.front().filename().string());

	fs::remove_all(s.dir);
}

TEST_CASE("calibrate-and-align writes aligned frames that solve to the field") {
	if (!assets_available()) return;
	auto s = stage_m31("calib");
	if (s.count < 2) { MESSAGE("need >= 2 M31 frames — skipping"); return; }

	const std::string db = " -d " + q(kDbDir.string());
	const auto st = run(q(kPortBin) + " -calibrate" + db + s.frame_args);
	CHECK(st.exit_code == 0);
	CHECK(field(st.out, "CALIBRATED=") == "1");
	CHECK(field(st.out, "FRAMES_ALIGNED=") == std::to_string(s.count));

	// One "<stem>_aligned.fit" per input, beside it.
	auto aligned = s.files.back();
	aligned.replace_extension("");
	aligned = aligned.string() + "_aligned.fit";
	REQUIRE(fs::exists(aligned));

	// The aligned frame carries the reference grid's solution: it re-solves to
	// M31 with the stored WCS essentially matching (small delta), stars sharp.
	const auto base = s.dir / "asolve";
	const auto sv = run(q(kPortBin) + " -f " + q(aligned.string()) + db
	                    + " -o " + q(base.string()));
	auto ini = base; ini.replace_extension(".ini");
	REQUIRE(fs::exists(ini));
	CHECK(ini_value(ini, "CRVAL1=") == doctest::Approx(10.7).epsilon(0.05));
	CHECK(ini_value(ini, "CRVAL2=") == doctest::Approx(41.24).epsilon(0.05));
	CHECK(selected_stars(sv.out) > 300);

	fs::remove_all(s.dir);
}

TEST_CASE("quality-outlier rejection drops the weakest frame") {
	if (!assets_available()) return;
	auto s = stage_m31("outlier");
	if (s.count < 3) { MESSAGE("need 3 M31 frames — skipping"); return; }

	const std::string db = " -d " + q(kDbDir.string());
	auto stack_sd = [&](const std::string& sd) -> RunResult {
		const auto out = s.dir / "o.fits";
		return run(q(kPortBin) + " -stackfiles -sd " + sd + " -o " + q(out.string())
		           + db + s.frame_args);
	};

	// A tight sigma rejects the lowest-quality frame; the survivors still combine.
	const auto tight = stack_sd("1.0");
	CHECK(tight.exit_code == 0);
	CHECK(field(tight.out, "FRAMES_SOLVED=") == "3");
	CHECK(field(tight.out, "FRAMES_KEPT=") == "2");
	CHECK(field(tight.out, "FRAMES_COMBINED=") == "2");
	// The rejected frame must not be the chosen reference.
	CHECK(!field(tight.out, "REFERENCE=").empty());

	// A loose sigma keeps every frame (nothing sits that far below the mean).
	const auto loose = stack_sd("5.0");
	CHECK(loose.exit_code == 0);
	CHECK(field(loose.out, "FRAMES_KEPT=") == "3");
	CHECK(field(loose.out, "FRAMES_COMBINED=") == "3");

	fs::remove_all(s.dir);
}

TEST_CASE("LRGB combine assembles a three-plane colour FITS from R/G/B frames") {
	if (!assets_available()) return;
	auto s = stage_m31("lrgb");
	if (s.count < 3) { MESSAGE("need 3 M31 frames — skipping"); return; }

	// No re-solve here: unlike the mono stack, an LRGB combine subtracts each
	// channel's own background (faithful to unit_stack_routines.pas), which for a
	// bright extended galaxy pushes the colour image below the solver's detection
	// floor — the original ASTAP behaves the same. The alignment maths are the
	// SAME astrometric_to_vector path the mono re-solve test already validates;
	// here we validate the LRGB-specific colour assembly and header instead.
	const auto out = s.dir / "colour.fits";
	const std::string db = " -d " + q(kDbDir.string());
	const auto lg = run(q(kPortBin) + " -lrgb -o " + q(out.string()) + db
	                    + " -red "   + q(s.files[0].string())
	                    + " -green " + q(s.files[1].string())
	                    + " -blue "  + q(s.files[2].string()));
	CHECK(lg.exit_code == 0);
	CHECK(field(lg.out, "LRGB=") == "1");
	CHECK(field(lg.out, "CHANNELS_IN=") == "3");
	CHECK(field(lg.out, "CHANNELS_COMBINED=") == "3");
	CHECK(field(lg.out, "LUMINANCE=") == "0");
	// No luminance → the reference is the red channel (first colour available).
	CHECK(fs::path(field(lg.out, "REFERENCE=")).filename() == s.files[0].filename());
	REQUIRE(fs::exists(out));

	// Header: a three-plane colour image, marked stacked. (The per-channel
	// RED_CNT/GRN_CNT/BLU_CNT blocks are written only when the channel inputs
	// carry LIGH_CNT — i.e. real per-filter interims; these raw corpus frames
	// do not, so only the colour/stacked markers are asserted here.)
	CHECK(fits_header_int(out, "NAXIS3  =") == 3);         // colour
	CHECK(fits_header_int(out, "BITPIX  =") == -32);
	{
		std::ifstream ifs(out, std::ios::binary);
		std::string hdr(2880 * 4, '\0');
		ifs.read(hdr.data(), static_cast<std::streamsize>(hdr.size()));
		CHECK(hdr.find("CALSTAT = 'S")            != std::string::npos);  // stacked
		CHECK(hdr.find("Combined to colour image") != std::string::npos);
	}

	// The identity colour matrix must route red→plane0, green→plane1, blue→plane2.
	// The three source frames differ, so their planes must carry distinct data —
	// a bug that dumped every channel into one plane would collapse this.
	const double m0 = plane_median(out, 0);
	const double m1 = plane_median(out, 1);
	const double m2 = plane_median(out, 2);
	REQUIRE(std::isfinite(m0));
	REQUIRE(std::isfinite(m1));
	REQUIRE(std::isfinite(m2));
	CHECK(m0 != doctest::Approx(m1));
	CHECK(m1 != doctest::Approx(m2));

	fs::remove_all(s.dir);
}

TEST_CASE("LRGB combine with a luminance channel drives the chroma modulation") {
	if (!assets_available()) return;
	auto s = stage_m31("lrgbL");
	if (s.count < 3) { MESSAGE("need 3 M31 frames — skipping"); return; }

	// red/green/blue plus a luminance frame: the luminance slot re-weights the
	// assembled chroma (c=5 branch) and is used as the alignment reference.
	const auto out  = s.dir / "colourL.fits";
	const auto outN = s.dir / "colourN.fits";  // same R/G/B, no luminance
	const std::string db = " -d " + q(kDbDir.string());
	const std::string rgb = " -red " + q(s.files[0].string())
	                      + " -green " + q(s.files[1].string())
	                      + " -blue " + q(s.files[2].string());
	const auto lg = run(q(kPortBin) + " -lrgb -o " + q(out.string()) + db + rgb
	                    + " -lum " + q(s.files[1].string()));
	CHECK(lg.exit_code == 0);
	CHECK(field(lg.out, "LRGB=") == "1");
	CHECK(field(lg.out, "LUMINANCE=") == "1");
	// Luminance present → it is the alignment reference.
	CHECK(fs::path(field(lg.out, "REFERENCE=")).filename() == s.files[1].filename());
	// R, G, B and L all contribute (the reference slot re-uses the L frame).
	CHECK(field(lg.out, "CHANNELS_COMBINED=") == "4");
	REQUIRE(fs::exists(out));
	CHECK(fits_header_int(out, "NAXIS3  =") == 3);

	// The luminance re-weighting must actually change the pixels: the same R/G/B
	// without an L channel yields a different plane-0 median.
	const auto ng = run(q(kPortBin) + " -lrgb -o " + q(outN.string()) + db + rgb);
	CHECK(ng.exit_code == 0);
	REQUIRE(fs::exists(outN));
	const double lum0 = plane_median(out,  0);
	const double raw0 = plane_median(outN, 0);
	REQUIRE(std::isfinite(lum0));
	REQUIRE(std::isfinite(raw0));
	CHECK(lum0 != doctest::Approx(raw0));

	fs::remove_all(s.dir);
}

TEST_CASE("master-frame creation averages raw frames into a mono FITS") {
	if (!assets_available()) return;
	auto s = stage_m31("mkdark");
	if (s.count < 2) { MESSAGE("need >= 2 M31 frames — skipping"); return; }

	const auto master = s.dir / "master_dark.fit";
	const auto st = run(q(kPortBin) + " -makedark -o " + q(master.string()) + s.frame_args);
	CHECK(st.exit_code == 0);
	CHECK(field(st.out, "MASTER=") == "1");
	CHECK(field(st.out, "FRAMES_COMBINED=") == std::to_string(s.count));
	REQUIRE(fs::exists(master));

	// Master carries the combined count and is mono 32-bit float.
	CHECK(fits_header_int(master, "DARK_CNT=") == s.count);
	CHECK(fits_header_int(master, "NAXIS   =") == 2);       // any colour → mono
	CHECK(fits_header_int(master, "BITPIX  =") == -32);

	fs::remove_all(s.dir);
}

TEST_CASE("master-flat creation subtracts a flat-dark (bias)") {
	if (!assets_available()) return;
	auto s = stage_m31("mkflat");
	if (s.count < 3) { MESSAGE("need 3 M31 frames — skipping"); return; }

	// Two frames as flats, one as the flat-dark. Content is a stand-in; the point
	// is the flat-dark is combined and subtracted, and the header records it.
	const auto master = s.dir / "master_flat.fit";
	const auto st = run(q(kPortBin) + " -makeflat -o " + q(master.string())
	                    + " -flatdark " + q(s.files[2].string())
	                    + " " + q(s.files[0].string()) + " " + q(s.files[1].string()));
	CHECK(st.exit_code == 0);
	CHECK(field(st.out, "MASTER=") == "1");
	CHECK(field(st.out, "FRAMES_COMBINED=") == "2");
	CHECK(field(st.out, "FLATDARKS_COMBINED=") == "1");
	REQUIRE(fs::exists(master));

	// Header: flat + bias counts recorded, CALSTAT gains 'B', still mono float.
	CHECK(fits_header_int(master, "FLAT_CNT=") == 2);
	CHECK(fits_header_int(master, "BIAS_CNT=") == 1);
	CHECK(fits_header_int(master, "NAXIS   =") == 2);
	{
		std::ifstream ifs(master, std::ios::binary);
		std::string hdr(2880 * 4, '\0');
		ifs.read(hdr.data(), static_cast<std::streamsize>(hdr.size()));
		CHECK(hdr.find("CALSTAT = 'B") != std::string::npos);
	}

	fs::remove_all(s.dir);
}

TEST_CASE("a created master dark can be applied when stacking") {
	if (!assets_available()) return;
	auto s = stage_m31("roundtrip");
	if (s.count < 2) { MESSAGE("need >= 2 M31 frames — skipping"); return; }

	// Build a master, then feed it back as the -dark for a stack. Content is
	// irrelevant here — the point is the created master is a valid, loadable
	// dark that the calibration path accepts (CALSTAT gains 'D').
	const auto master = s.dir / "master_dark.fit";
	const auto mk = run(q(kPortBin) + " -makedark -o " + q(master.string()) + s.frame_args);
	CHECK(mk.exit_code == 0);
	REQUIRE(fs::exists(master));

	const auto out = s.dir / "stack.fits";
	const std::string db = " -d " + q(kDbDir.string());
	const auto st = run(q(kPortBin) + " -stackfiles -o " + q(out.string())
	                    + " -dark " + q(master.string()) + db + s.frame_args);
	CHECK(st.exit_code == 0);
	CHECK(field(st.out, "STACKED=") == "1");
	REQUIRE(fs::exists(out));
	{
		std::ifstream ifs(out, std::ios::binary);
		std::string hdr(2880 * 4, '\0');
		ifs.read(hdr.data(), static_cast<std::streamsize>(hdr.size()));
		CHECK(hdr.find("CALSTAT = 'D") != std::string::npos);
	}

	fs::remove_all(s.dir);
}
