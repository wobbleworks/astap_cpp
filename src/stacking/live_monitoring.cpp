///----------------------------------------
///      @file live_monitoring.cpp
///   @ingroup ASTAP++
///     @brief Implementation of @ref LiveMonitorSession.
///    @author Ported from Han Kleijn's unit_live_monitoring.pas (ASTAP).
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "live_monitoring.h"

#include "stack.h"          // apply_dark_and_flat
#include "../core/demosaic.h"
#include "../core/fits.h"
#include "../core/globals.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <fstream>
#include <string_view>
#include <system_error>
#include <thread>

namespace fs = std::filesystem;

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

namespace {

/// Same acceptable-extensions list used by LiveStackSession. Kept file-local
/// to avoid a shared header for a one-line constant.
constexpr std::array<std::string_view, 21> kSupportedExts = {
	".fit",  ".fits", ".png", ".jpg", ".bmp", ".tif", ".tiff", ".xisf",
	".raw",  ".crw",  ".cr2", ".cr3", ".kdc", ".dcr", ".mrw",  ".arw",
	".nef",  ".nrw",  ".dng", ".orf", ".raf"
};

[[nodiscard]] bool ext_supported(const fs::path& p) {
	auto e = p.extension().string();
	std::transform(e.begin(), e.end(), e.begin(),
		[](unsigned char c) { return std::tolower(c); });
	return std::find(kSupportedExts.begin(), kSupportedExts.end(),
		std::string_view{e}) != kSupportedExts.end();
}

// Returns true only if the file is openable right now — filters partially-
// written files still in the capture software's hands.
[[nodiscard]] bool file_readable(const fs::path& p) {
	std::ifstream f(p, std::ios::binary);
	return f.good();
}

}  // namespace

LiveMonitorSession::LiveMonitorSession(fs::path watch_dir)
	: watch_dir_(std::move(watch_dir)) {
}

bool LiveMonitorSession::file_available(const fs::path& dir, fs::path& out_file) {
	auto ec = std::error_code{};
	if (!fs::is_directory(dir, ec)) {
		return false;
	}

	// Pick the OLDEST unprocessed file so captures are displayed in the
	// order they were written — matters when the user is actively capturing
	// and we poll less frequently than the cadence.
	auto best_time = fs::file_time_type::max();
	auto best_path = fs::path{};

	for (const auto& entry : fs::directory_iterator(dir, ec)) {
		if (ec) break;
		if (!entry.is_regular_file(ec)) continue;
		const auto& p = entry.path();
		if (!ext_supported(p)) continue;

		const auto t = entry.last_write_time(ec);
		if (ec) continue;
		if (t <= latest_time_) continue;         // already processed
		if (t < best_time) {
			best_time = t;
			best_path = p;
		}
	}

	if (best_path.empty()) return false;
	if (!file_readable(best_path)) return false;   // still being written

	out_file = best_path;
	latest_time_ = best_time;
	return true;
}

bool LiveMonitorSession::process_frame(const fs::path& filename) {
	if (!astap::core::load_fits(filename, /*light=*/true, /*load_data=*/true,
	                             /*update_memo=*/true, /*get_ext=*/0,
	                             astap::memo1_lines, astap::head,
	                             astap::img_loaded)) {
		emit_message("Error loading " + filename.string());
		return false;
	}

	// Apply calibration if a master dark/flat is loaded. This is idempotent
	// (the apply checks internal state) and mirrors LiveStackSession.
	(void)apply_dark_and_flat(astap::img_loaded, astap::head);

	// Auto-demosaic single-channel frames with a Bayer pattern.
	if (astap::head.naxis3 == 1 && !astap::bayerpat.empty()) {
		const auto pattern = astap::core::get_demosaic_pattern(
			2, astap::xbayroff, astap::ybayroff, astap::roworder);
		astap::core::demosaic_bayer(astap::img_loaded, astap::head, pattern,
			astap::core::DemosaicMethod::Bilinear);
	}

	++total_;
	astap::filename2 = filename.string();  // tell the GUI which file is shown
	emit_message("Loaded " + filename.filename().string() +
	             " (frame " + std::to_string(total_) + ").");
	emit_frame_loaded();
	return true;
}

void LiveMonitorSession::run() {
	astap::live_stacking.store(true);
	astap::pause_pressed.store(false);
	astap::esc_pressed.store(false);
	total_       = 0;
	latest_time_ = fs::file_time_type{};   // accept anything on first pass

	emit_message("Monitor started. Watching " + watch_dir_.string());

	auto waiting = false;
	while (!astap::esc_pressed.load()) {
		auto filename = fs::path{};
		const auto have_file = !astap::pause_pressed.load() &&
		                       file_available(watch_dir_, filename);

		if (!have_file) {
			if (!waiting) {
				emit_message(astap::pause_pressed.load()
					? "Paused."
					: "Waiting for files…");
				waiting = true;
			}
			std::this_thread::sleep_for(std::chrono::milliseconds(1000));
			continue;
		}
		waiting = false;

		(void)process_frame(filename);
	}

	astap::live_stacking.store(false);
	emit_message("Monitor stopped. Displayed " + std::to_string(total_) +
	             " frames.");
}

void LiveMonitorSession::emit_frame_loaded() {
	if (frame_hook_) frame_hook_(total_);
}

void LiveMonitorSession::emit_message(const std::string& msg) {
	if (message_hook_) message_hook_(msg);
}

} // namespace
