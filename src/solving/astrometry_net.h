///----------------------------------------
///      @file astrometry_net.h
///   @ingroup ASTAP++
///     @brief Thin wrapper around a locally-installed @c solve-field from
///            astrometry.net, for use as a fallback solver.
///   @details Ported from Han Kleijn's @c unit_astrometry_net.pas. Spawns
///            @c solve-field on the given FITS file, streams its stdout to
///            a caller-supplied logger, and on success renames the
///            resulting @c <name>.new over the original (backing the
///            original up as @c <name>.bak — matching the Pascal
///            behaviour). Optionally removes the scratch @c .wcs /
///            @c .corr / @c .match / @c .rdls / @c .solved / @c .axy /
///            @c .xyls files that @c solve-field emits alongside.
///    @author Created by John Stephen on 4/23/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <filesystem>
#include <functional>
#include <string>

///----------------------------------------
namespace astap::solving {
///----------------------------------------

///----------------------------------------
/// @struct AstrometryNetOptions
/// @brief Parameters for a single @ref astrometry_net invocation.
///----------------------------------------

struct AstrometryNetOptions {
	/// @brief FITS file to solve. Must exist on disk.
	std::filesystem::path fits_path;

	/// @brief On POSIX: the full path to the @c solve-field binary or a
	///        directory containing it (empty = rely on @c $PATH). On
	///        Windows: the full path to a @c bash.exe — either Cygwin's
	///        (e.g. @c C:\cygwin\bin\bash.exe) or WSL's
	///        @c C:\Windows\System32\bash.exe. The engine auto-detects
	///        which is in use by looking for @c "System32" in the path.
	std::filesystem::path solve_field_path;

	/// @brief Extra command-line arguments, whitespace-separated. Common
	///        example: @c "--downsample 2 --objs 150". Appended verbatim
	///        after the always-on @c --overwrite and @c --no-plots flags.
	std::string extra_args;

	/// @brief When true, remove scratch @c .wcs / @c .corr / @c .match /
	///        @c .rdls / @c .solved / @c .axy / @c .xyls files left by
	///        @c solve-field alongside the FITS after a successful solve.
	bool cleanup_tmp = true;

	/// @brief Windows only: show the Cygwin/WSL console window during the
	///        solve. Ignored on POSIX (stdout is always streamed to the
	///        logger).
	bool show_console = true;

	/// @brief Windows only: wrap the launch in @c "cmd.exe /k ..." so the
	///        console window remains open after @c solve-field exits, for
	///        inspecting final messages. Ignored on POSIX.
	bool keep_console_open = false;
};

///----------------------------------------
/// @brief Callback for streaming @c solve-field stdout lines.
///        Empty function = silent.
///----------------------------------------

using MessageHook = std::function<void(const std::string&)>;

///----------------------------------------
/// @brief Run local @c solve-field on @p opts.fits_path.
/// @details Blocks until @c solve-field exits. On success renames
///          @c <name>.fits → @c <name>.bak and @c <name>.new →
///          @c <name>.fits, optionally removes scratch files.
///
///          Windows is not supported in this build — the Pascal original
///          branches between Cygwin and WSL there, neither of which is
///          wired up yet. Returns @c false with an explanatory log line
///          on Windows.
/// @param opts  Run configuration.
/// @param log   Line logger; each line from @c solve-field stdout is
///              forwarded here. May be empty.
/// @return @c true iff @c solve-field exited successfully AND the
///         expected @c .new file was found and renamed into place.
///----------------------------------------

[[nodiscard]] bool astrometry_net(const AstrometryNetOptions& opts,
                                  MessageHook log = {});

} // namespace astap::solving
