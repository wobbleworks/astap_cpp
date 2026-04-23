///----------------------------------------
///      @file astrometry_net.cpp
///   @ingroup ASTAP++
///     @brief Implementation of @ref astap::solving::astrometry_net.
///    @author Created by John Stephen on 4/23/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "astrometry_net.h"

#include "../core/platform.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

///----------------------------------------
namespace astap::solving {
///----------------------------------------

namespace fs = std::filesystem;

namespace {

// Wrap a token in double quotes, escaping embedded quotes. Good enough for
// the handful of paths we pass to solve-field — not a full POSIX shell quoter.
[[nodiscard]] std::string shell_quote(std::string_view s) {
	auto out = std::string{"\""};
	for (auto c : s) {
		if (c == '"') {
			out += "\\\"";
		} else {
			out.push_back(c);
		}
	}
	out += '"';
	return out;
}

// Split @p args on whitespace. Quotes aren't honoured because solve-field's
// own option values don't include spaces in any realistic user input.
[[nodiscard]] std::vector<std::string> tokenise_whitespace(std::string_view args) {
	auto out = std::vector<std::string>{};
	auto ss = std::istringstream{std::string{args}};
	auto tok = std::string{};
	while (ss >> tok) {
		out.emplace_back(std::move(tok));
	}
	return out;
}

#if !defined(_WIN32)

// POSIX launch: popen solve-field, stream stdout line-by-line to @p log.
[[nodiscard]] bool launch_posix(const AstrometryNetOptions& opts,
                                const MessageHook& log) {
	auto err = std::error_code{};

	// Resolve the solve-field binary. The Pascal UI lets the user point
	// at a directory containing solve-field, so accept both forms here.
	auto solve_field = fs::path{};
	if (opts.solve_field_path.empty()) {
		solve_field = "solve-field";  // rely on $PATH
	} else if (fs::is_directory(opts.solve_field_path, err)) {
		solve_field = opts.solve_field_path / "solve-field";
	} else {
		solve_field = opts.solve_field_path;
	}
	if (solve_field.has_parent_path() && !fs::exists(solve_field, err)) {
		if (log) log("Cannot find solve-field at: " + solve_field.string() +
		             " — install astrometry.net or correct the path.");
		return false;
	}

	auto cmd = std::string{};
	cmd += shell_quote(solve_field.string());
	cmd += ' ';
	cmd += shell_quote(opts.fits_path.string());
	cmd += " --overwrite --no-plots";
	for (const auto& t : tokenise_whitespace(opts.extra_args)) {
		cmd += ' ';
		cmd += shell_quote(t);
	}
	cmd += " 2>&1";  // merge stderr so progress lines reach the log

	if (log) log("Running: " + cmd);

	auto* pipe = ::popen(cmd.c_str(), "r");
	if (!pipe) {
		if (log) log("Failed to spawn solve-field.");
		return false;
	}
	char buf[4096];
	while (std::fgets(buf, sizeof(buf), pipe)) {
		if (log) {
			auto line = std::string{buf};
			while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) {
				line.pop_back();
			}
			log(line);
		}
	}
	const auto rc = ::pclose(pipe);
	if (rc != 0) {
		if (log) log("solve-field exited with non-zero status.");
		return false;
	}
	return true;
}

#else // _WIN32

// Windows launch: spawn a Cygwin or WSL bash.exe pointing at solve-field.
// Mirrors @c unit_astrometry_net.pas lines 181–199.
[[nodiscard]] bool launch_windows(const AstrometryNetOptions& opts,
                                  const MessageHook& log) {
	auto err = std::error_code{};

	if (opts.solve_field_path.empty()) {
		if (log) log("Windows: point 'solve-field' at a Cygwin or WSL "
		             "bash.exe that has astrometry.net installed.");
		return false;
	}
	if (!fs::exists(opts.solve_field_path, err)) {
		if (log) log("Cannot find bash.exe at: " +
		             opts.solve_field_path.string());
		return false;
	}

	const auto bash_path = opts.solve_field_path.string();
	const auto is_wsl = (bash_path.find("System32") != std::string::npos
	                   || bash_path.find("system32") != std::string::npos);

	// Build the trailing solve-field argument block, common to both
	// flavours. No per-token quoting — Pascal strips quotes before handing
	// to its ExecuteAndWait, trusting the user not to embed spaces in
	// their extra-args.
	auto args = std::string{"--overwrite --no-plots"};
	if (!opts.extra_args.empty()) {
		args += ' ';
		args += opts.extra_args;
	}

	// Both Cygwin and WSL accept forward slashes in Windows paths.
	auto fwd_path = opts.fits_path.string();
	std::replace(fwd_path.begin(), fwd_path.end(), '\\', '/');

	auto cmd = std::string{};
	if (!is_wsl) {
		// Cygwin: bash --login lets cygwin set up its environment,
		// then solve-field handles Windows-style paths natively.
		cmd = shell_quote(bash_path) + " --login solve-field \""
		    + fwd_path + "\" " + args;
	} else {
		// WSL: translate C:/foo/bar.fit -> /mnt/c/foo/bar.fit (lowercase
		// drive letter — WSL requires it).
		auto wsl_path = fwd_path;
		if (wsl_path.size() >= 2 && wsl_path[1] == ':') {
			const auto drive = static_cast<char>(
				std::tolower(static_cast<unsigned char>(wsl_path[0])));
			wsl_path = std::string{"/mnt/"} + drive + wsl_path.substr(2);
		}
		// Single-quote the inner path so bash -c sees it as one arg.
		cmd = shell_quote(bash_path) + " -c \"solve-field '"
		    + wsl_path + "' " + args + "\"";
	}

	if (opts.keep_console_open) {
		cmd = "cmd.exe /k " + cmd;
	}

	if (log) log("Running: " + cmd);

	// execute_and_wait uses CreateProcess + STARTF_USESHOWWINDOW to honour
	// the show-console flag. Output is not captured on Windows — the user
	// watches the console window (or the closed window, if show_console
	// is false and keep_console_open is false).
	const auto ok = astap::core::execute_and_wait(cmd, opts.show_console);
	if (!ok) {
		if (log) log("solve-field exited with non-zero status "
		             "(or failed to launch — check the bash.exe path).");
		return false;
	}
	return true;
}

#endif

}  // namespace

bool astrometry_net(const AstrometryNetOptions& opts, MessageHook log) {
	auto err = std::error_code{};
	if (!fs::exists(opts.fits_path, err)) {
		if (log) log("Input FITS does not exist: " + opts.fits_path.string());
		return false;
	}

#if defined(_WIN32)
	if (!launch_windows(opts, log)) {
		return false;
	}
#else
	if (!launch_posix(opts, log)) {
		return false;
	}
#endif

	// solve-field writes <basename>.new next to the input. If it's missing,
	// the solve failed even if the process exit code was zero.
	auto new_file = opts.fits_path;
	new_file.replace_extension(".new");
	if (!fs::exists(new_file, err)) {
		if (log) log("solve-field did not produce " + new_file.filename().string() +
		             " — astrometry could not solve this field.");
		return false;
	}

	// Mirror the Pascal behaviour: back up the original, swap .new in.
	auto bak = opts.fits_path;
	bak.replace_extension(".bak");
	if (fs::exists(bak, err)) {
		fs::remove(bak, err);
	}
	fs::rename(opts.fits_path, bak, err);
	if (err) {
		if (log) log("Failed to rename original to .bak: " + err.message());
		return false;
	}
	fs::rename(new_file, opts.fits_path, err);
	if (err) {
		if (log) log("Failed to move .new into place: " + err.message() +
		             ". Attempting to restore original.");
		auto restore_err = std::error_code{};
		fs::rename(bak, opts.fits_path, restore_err);
		return false;
	}

	if (opts.cleanup_tmp) {
		const auto dir = opts.fits_path.parent_path().string();
		if (!dir.empty()) {
			for (const auto* mask : {"*.wcs", "*.corr", "*.match",
			                          "*.rdls", "*.solved", "*.axy", "*.xyls"}) {
				astap::core::delete_files(dir, mask);
			}
		}
	}

	if (log) log("Solved. Original saved as " + bak.filename().string() + ".");
	return true;
}

} // namespace astap::solving
