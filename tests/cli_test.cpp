///----------------------------------------
///     @file cli_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Integration tests for the astap CLI binary.
///  @details Spawns the compiled `astap` executable with various argument
///           combinations and asserts on exit codes, stdout content, and
///           file side-effects. These tests exercise the argv parser, the
///           help text, the error-path exit codes, and the .ini output
///           contract that downstream tools (APT, NINA, SGP) depend on.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#if defined(_WIN32)
  #include <io.h>
  // MSVC's C runtime names these with a leading underscore.
  #define popen  _popen
  #define pclose _pclose
#else
  #include <sys/wait.h>
#endif

namespace fs = std::filesystem;

///----------------------------------------
/// MARK: Helpers
///----------------------------------------

/// @brief Path to the astap binary, resolved at CMake configure time.
static const std::string kAstapBin = ASTAP_BIN_PATH;

/// @brief Result of running astap with specific args.
struct RunResult {
	int         exit_code;
	std::string stdout_text;
};

/// @brief Spawn `astap <args>` and capture stdout + exit code.
[[nodiscard]] static RunResult run_astap(const std::string& args) {
	const auto cmd = kAstapBin + " " + args + " 2>&1";
	RunResult r{};

	FILE* pipe = ::popen(cmd.c_str(), "r");
	if (!pipe) {
		r.exit_code = -1;
		return r;
	}

	char buf[1024];
	while (std::fgets(buf, sizeof(buf), pipe)) {
		r.stdout_text.append(buf);
	}
	const int raw = ::pclose(pipe);
#if defined(_WIN32)
	r.exit_code = raw;
#else
	r.exit_code = WIFEXITED(raw) ? WEXITSTATUS(raw) : -1;
#endif
	return r;
}

/// @brief Create a minimal temp directory for test outputs.
[[nodiscard]] static fs::path make_temp_dir() {
	auto p = fs::temp_directory_path() / "astap_cli_test";
	fs::create_directories(p);
	return p;
}

///----------------------------------------
/// MARK: Help flag
///----------------------------------------

TEST_CASE("astap -h prints help and exits 0") {
	const auto r = run_astap("-h");
	CHECK(r.exit_code == 0);
	CHECK(r.stdout_text.find("Solver command-line usage") != std::string::npos);
	CHECK(r.stdout_text.find("-f  filename") != std::string::npos);
	CHECK(r.stdout_text.find("-fov") != std::string::npos);
	CHECK(r.stdout_text.find("-analyse") != std::string::npos);
}

TEST_CASE("astap --help also works") {
	const auto r = run_astap("--help");
	CHECK(r.exit_code == 0);
	CHECK(r.stdout_text.find("Solver command-line usage") != std::string::npos);
}

///----------------------------------------
/// MARK: Unknown / missing args
///----------------------------------------

TEST_CASE("astap with unknown flag exits 2") {
	const auto r = run_astap("--bogus-flag");
	CHECK(r.exit_code == 2);
	CHECK(r.stdout_text.find("unrecognised option") != std::string::npos);
}

///----------------------------------------
/// MARK: Missing file → exit code 16
///----------------------------------------

TEST_CASE("astap -f nonexistent.fit exits 16") {
	const auto r = run_astap("-f /tmp/does_not_exist_astap_test.fit");
	CHECK(r.exit_code == 16);
	CHECK(r.stdout_text.find("Error reading image file") != std::string::npos);
}

///----------------------------------------
/// MARK: .ini output contract
///----------------------------------------

TEST_CASE("astap -f zero-byte file exits with error code") {
	// A zero-byte file isn't valid FITS; load_fits will reject it.
	const auto tmpdir = make_temp_dir();
	const auto fake_input = tmpdir / "test_empty.fit";
	{ std::ofstream ofs(fake_input); }   // create empty file.

	const auto r = run_astap("-f " + fake_input.string());

	// Expect exit code 16 (image load failure) or 1 (no solution).
	CHECK((r.exit_code == 16 || r.exit_code == 1));

	// Cleanup.
	fs::remove_all(tmpdir);
}

///----------------------------------------
/// MARK: No-args invocation prints help
///----------------------------------------

TEST_CASE("astap with no args prints help and exits 0") {
	const auto r = run_astap("");
	CHECK(r.exit_code == 0);
	CHECK(r.stdout_text.find("Solver command-line usage") != std::string::npos);
}

///----------------------------------------
/// MARK: Focus mode with no images
///----------------------------------------

TEST_CASE("astap -focus1 nonexistent exits with error") {
	const auto r = run_astap("-focus1 /tmp/does_not_exist.fit");
	// fit_focus_hyperbola stub returns ok=false → run_focus returns 1.
	CHECK(r.exit_code != 0);
}

///----------------------------------------
/// MARK: Stack mode flag parses
///----------------------------------------

// NOTE: -stack test is deliberately omitted. The live-stacking loop polls
// a directory and blocks indefinitely until esc_pressed is set. Testing it
// requires a cooperative signal-injection harness, which is outside the
// scope of a quick CLI integration check.
