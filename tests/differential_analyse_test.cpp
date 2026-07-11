///----------------------------------------
///     @file differential_analyse_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Differential tests: ASTAP++ vs the original astap_cli.
///  @details Drives BOTH the ported `astap` binary and the reference
///           `astap_cli` oracle with identical `-analyse` / `-extract`
///           arguments over a shared image corpus, then compares the
///           results with tolerance- and field-aware comparators. The
///           original ASTAP ships no test suite of its own, so it serves
///           purely as a behavioural oracle here.
///
///           This is star-detection phase 1: `-analyse` (median HFD + star
///           count from stdout) and `-extract` (per-star CSV). No star
///           database is required — detection is image-only.
///
///           The oracle (astap_cli CLI-2025.07.14) is a near-but-not-exact
///           match for the ported source (ASTAP 2024.11.13, the newest that
///           is still archived), so tolerances are deliberately loose during
///           bring-up: the goal is to catch gross divergence and record the
///           deltas, not to pin sub-percent agreement. Every case emits its
///           numbers via MESSAGE so the comparison is visible even on pass.
///
///           Gating: the whole suite skips gracefully (passing, with a
///           MESSAGE) when the oracle binary or a corpus image is absent, so
///           CI without the private corpus stays green.
///   @author Created by John Stephen on 7/11/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

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

/// @brief The ported astap binary under test.
static const std::string kPortBin   = ASTAP_BIN_PATH;
/// @brief The reference oracle (original astap_cli).
static const std::string kOracleBin = ASTAP_ORACLE_BIN;
/// @brief Directory holding the shared image corpus.
static const fs::path    kCorpusDir = ASTAP_CORPUS_DIR;

///----------------------------------------
/// MARK: Process + parsing helpers
///----------------------------------------

/// @brief Shell-quote a path/argument.
[[nodiscard]] static std::string q(const std::string& s) { return "\"" + s + "\""; }

struct RunResult {
	int         exit_code;
	std::string stdout_text;
};

/// @brief Run a fully-composed shell command, capturing stdout+stderr.
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

/// @brief Extract the numeric value following `key=` from analyse stdout.
[[nodiscard]] static std::optional<double> value_after(const std::string& text,
                                                       const std::string& key) {
	const auto pos = text.find(key);
	if (pos == std::string::npos) return std::nullopt;
	try {
		return std::stod(text.substr(pos + key.size()));
	} catch (...) {
		return std::nullopt;
	}
}

struct Analysis {
	std::optional<double> hfd_median;
	std::optional<double> stars;
};

[[nodiscard]] static Analysis parse_analysis(const std::string& text) {
	return { value_after(text, "HFD_MEDIAN="), value_after(text, "STARS=") };
}

struct Star {
	double x{}, y{}, hfd{}, snr{}, flux{};
};

/// @brief Parse an extract CSV (`x,y,hfd,snr,flux[,ra,dec]`), skipping header.
[[nodiscard]] static std::vector<Star> parse_csv(const fs::path& path) {
	std::vector<Star> out;
	std::ifstream ifs(path);
	std::string line;
	while (std::getline(ifs, line)) {
		if (line.empty() || (line[0] < '0' || line[0] > '9')) continue;  // header / blanks
		std::vector<double> f;
		std::size_t start = 0;
		while (start <= line.size()) {
			const auto comma = line.find(',', start);
			const auto tok = line.substr(start, comma - start);
			if (!tok.empty()) { try { f.push_back(std::stod(tok)); } catch (...) {} }
			if (comma == std::string::npos) break;
			start = comma + 1;
		}
		if (f.size() >= 5) out.push_back({ f[0], f[1], f[2], f[3], f[4] });
	}
	return out;
}

/// @brief Copy a corpus image into a fresh per-case temp dir so both binaries
///        write their side-effect files (`.csv`) in isolation. Returns the
///        staged input path (its parent dir is unique to @p case_name).
[[nodiscard]] static fs::path stage(const fs::path& src, const std::string& case_name) {
	auto dir = fs::temp_directory_path() / ("astap_diff_" + case_name);
	fs::remove_all(dir);
	fs::create_directories(dir);
	auto dst = dir / src.filename();
	fs::copy_file(src, dst, fs::copy_options::overwrite_existing);
	return dst;
}

///----------------------------------------
/// MARK: Skip gating
///----------------------------------------

/// @brief True when the oracle is present; otherwise MESSAGE + skip.
[[nodiscard]] static bool oracle_available() {
	if (fs::exists(kOracleBin)) return true;
	MESSAGE("oracle binary absent (" << kOracleBin << ") — skipping differential case");
	return false;
}

/// @brief Resolve a corpus image, or MESSAGE + return nullopt to skip.
[[nodiscard]] static std::optional<fs::path> corpus_image(const std::string& name) {
	auto p = kCorpusDir / name;
	if (fs::exists(p)) return p;
	MESSAGE("corpus image absent (" << p.string() << ") — skipping differential case");
	return std::nullopt;
}

///----------------------------------------
/// MARK: Comparators
///----------------------------------------

/// @brief The port is expected to reproduce the oracle's star count exactly.
///        Always logs both values so a regression shows its magnitude.
static void expect_star_count_agrees(const std::string& label, double oracle, double port) {
	MESSAGE(label << ": STARS oracle=" << oracle << " port=" << port
	              << " Δ=" << std::abs(oracle - port));
	CHECK(oracle == port);
}

/// @brief Nearest-neighbour star-list match within `tol_px`. Reports match
///        rate and per-field medians; asserts a bring-up match-rate floor.
static void compare_star_lists(const std::string& label,
                               const std::vector<Star>& oracle,
                               const std::vector<Star>& port) {
	constexpr double tol_px = 2.0;
	std::vector<char> port_taken(port.size(), 0);
	int matched = 0;
	std::vector<double> dhfd, dflux_pct;
	for (const auto& o : oracle) {
		int best = -1;
		double best_d2 = tol_px * tol_px;
		for (std::size_t j = 0; j < port.size(); ++j) {
			if (port_taken[j]) continue;
			const double dx = o.x - port[j].x, dy = o.y - port[j].y;
			const double d2 = dx * dx + dy * dy;
			if (d2 <= best_d2) { best_d2 = d2; best = static_cast<int>(j); }
		}
		if (best >= 0) {
			port_taken[best] = 1;
			++matched;
			dhfd.push_back(std::abs(o.hfd - port[best].hfd));
			if (o.flux > 0.0)
				dflux_pct.push_back(std::abs(o.flux - port[best].flux) / o.flux * 100.0);
		}
	}
	const auto denom = std::max(oracle.size(), port.size());
	const double rate = denom ? static_cast<double>(matched) / denom : 1.0;
	const auto median = [](std::vector<double>& v) -> double {
		if (v.empty()) return 0.0;
		std::sort(v.begin(), v.end());
		return v[v.size() / 2];
	};
	const auto med_hfd = median(dhfd);
	const auto med_flux = median(dflux_pct);
	MESSAGE(label << ": oracle=" << oracle.size() << " port=" << port.size()
	              << " matched=" << matched
	              << " only-oracle=" << (oracle.size() - matched)
	              << " only-port=" << (port.size() - matched)
	              << " match-rate=" << rate
	              << " median|Δhfd|=" << med_hfd
	              << " median|Δflux|%=" << med_flux);
	// Exact equivalence: every star matched, no strays either side, and the
	// per-star HFD/flux identical to display precision.
	CHECK(oracle.size() == port.size());
	CHECK(static_cast<std::size_t>(matched) == oracle.size());
	CHECK(rate == doctest::Approx(1.0));
	CHECK(med_hfd == doctest::Approx(0.0));
	CHECK(med_flux == doctest::Approx(0.0));
}

///----------------------------------------
/// MARK: Case drivers
///----------------------------------------

static void analyse_case(const std::string& name, const std::string& file,
                         double snr, bool expect_rich,
                         const std::string& extra_args = "") {
	if (!oracle_available()) return;
	auto img = corpus_image(file);
	if (!img) return;

	const auto staged = stage(*img, name);
	const std::string common = " -f " + q(staged.string())
	                         + " -analyse " + std::to_string(snr)
	                         + (extra_args.empty() ? "" : " " + extra_args);
	const auto o = parse_analysis(run(q(kOracleBin) + common).stdout_text);
	const auto p = parse_analysis(run(q(kPortBin)   + common).stdout_text);

	REQUIRE(o.stars.has_value());
	REQUIRE(p.stars.has_value());
	expect_star_count_agrees(name, *o.stars, *p.stars);

	// A rich field must stay rich on both sides; a sparse field must stay
	// sparse. This catches a wrong detection *class* before any tolerance.
	if (expect_rich) {
		CHECK(*o.stars > 100);
		CHECK(*p.stars > 100);
	} else {
		CHECK(*o.stars < 50);
		CHECK(*p.stars < 50);
	}

	if (o.hfd_median && p.hfd_median && *o.hfd_median < 90.0 && *p.hfd_median < 90.0) {
		MESSAGE(name << ": HFD_MEDIAN oracle=" << *o.hfd_median
		             << " port=" << *p.hfd_median);
		// Both print to one decimal; 0.05 tolerance means same rounded value.
		CHECK(std::abs(*o.hfd_median - *p.hfd_median) <= 0.05);
	}
}

static void extract_case(const std::string& name, const std::string& file, double snr) {
	if (!oracle_available()) return;
	auto img = corpus_image(file);
	if (!img) return;

	const auto o_in = stage(*img, name + "_oracle");
	const auto p_in = stage(*img, name + "_port");
	auto o_csv = o_in; o_csv.replace_extension(".csv");
	auto p_csv = p_in; p_csv.replace_extension(".csv");

	(void)run(q(kOracleBin) + " -f " + q(o_in.string()) + " -extract " + std::to_string(snr));
	(void)run(q(kPortBin)   + " -f " + q(p_in.string()) + " -extract " + std::to_string(snr));

	REQUIRE_MESSAGE(fs::exists(o_csv), "oracle wrote no CSV at " << o_csv.string());
	REQUIRE_MESSAGE(fs::exists(p_csv), "port wrote no CSV at " << p_csv.string()
	                << " (is the stack.cpp CSV writer implemented?)");

	compare_star_lists(name, parse_csv(o_csv), parse_csv(p_csv));
}

///----------------------------------------
/// MARK: Tests
///----------------------------------------

TEST_CASE("oracle reports its version") {
	if (!oracle_available()) return;
	const auto r = run(q(kOracleBin) + " -h");
	MESSAGE("oracle: " << r.stdout_text.substr(0, r.stdout_text.find('\n')));
	CHECK(r.stdout_text.find("version") != std::string::npos);
}

TEST_CASE("analyse: M31 rich field agrees") {
	analyse_case("M31_analyse", "M31_a.fits", 20, /*expect_rich=*/true);
}

TEST_CASE("analyse: ngc_990 sparse field agrees") {
	analyse_case("ngc990i_analyse", "ngc_990.i.unconv.fits", 10, /*expect_rich=*/false);
}

TEST_CASE("analyse: -s override honored identically") {
	// The port must feed -s into analyse's max_stars exactly as the original
	// does (astap::max_stars_setting <- max_stars1). Counts stay equal.
	analyse_case("M31_analyse_s2000", "M31_a.fits", 20, /*expect_rich=*/true, "-s 2000");
}

TEST_CASE("extract: M31 per-star list agrees") {
	extract_case("M31_extract", "M31_a.fits", 20);
}

TEST_CASE("extract: ngc_990 per-star list agrees") {
	extract_case("ngc990i_extract", "ngc_990.i.unconv.fits", 10);
}
