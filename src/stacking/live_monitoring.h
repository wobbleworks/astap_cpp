///----------------------------------------
///      @file live_monitoring.h
///   @ingroup ASTAP++
///     @brief Directory watcher that displays each newly-landed image as it
///            arrives, with optional dark/flat + OSC demosaic. No alignment,
///            no stacking — a minimal "live view" companion to
///            @ref astap::stacking::LiveStackSession.
///    @author Ported from Han Kleijn's unit_live_monitoring.pas (ASTAP).
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <filesystem>
#include <functional>
#include <string>

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

///----------------------------------------
/// @class LiveMonitorSession
/// @brief Polls a directory and loads each new image file into the engine
///        globals (@c astap::img_loaded, @c astap::head) so the GUI can
///        refresh its viewer. Pending master-dark / master-flat / OSC
///        Bayer settings apply automatically.
/// @details Uses modification-time tracking (no file rename) so the
///          capture software's output directory is left untouched. Picks
///          the oldest unprocessed file each poll, preserving FIFO order.
///----------------------------------------

class LiveMonitorSession final {
public:
	/// @brief Notified after each newly-loaded frame. @p total is the
	///        cumulative frame count since @ref run started.
	using FrameLoadedHook = std::function<void(int total)>;

	/// @brief Notified with a human-readable status line.
	using MessageHook = std::function<void(const std::string&)>;

	explicit LiveMonitorSession(std::filesystem::path watch_dir);

	void set_frame_loaded_hook(FrameLoadedHook hook) { frame_hook_   = std::move(hook); }
	void set_message_hook     (MessageHook     hook) { message_hook_ = std::move(hook); }

	/// @brief Main polling loop. Blocks until @c astap::esc_pressed is set.
	void run();

private:
	[[nodiscard]] bool process_frame(const std::filesystem::path& filename);

	[[nodiscard]] bool file_available(const std::filesystem::path& dir,
	                                  std::filesystem::path&       out_file);

	void emit_frame_loaded();
	void emit_message(const std::string& msg);

	std::filesystem::path  watch_dir_;
	std::filesystem::file_time_type latest_time_{};
	int                    total_ = 0;

	FrameLoadedHook frame_hook_;
	MessageHook     message_hook_;
};

} // namespace
