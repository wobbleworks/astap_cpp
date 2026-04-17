///----------------------------------------
///      @file live_stacking.h
///   @ingroup ASTAP++
///     @brief Live-stacking watcher: polls a directory for FITS/image files,
///            plate-solves and aligns each onto a reference, accumulates a
///            running average, and optionally exports a JPG snapshot.
///   @details Preserves the original polling-based loop (no inotify or
///            FSEvents). Heavy dependence on astap_main globals (image buffer,
///            FITS header, plate-solve state) is recorded with TODO markers
///            for the modules not yet ported.
///    @author Ported from Han Kleijn's unit_live_stacking.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <atomic>
#include <chrono>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

#include "../types.h"
#include "../core/globals.h"   // esc_pressed, pause_pressed, live_stacking

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

using astap::ImageArray;

///----------------------------------------
/// MARK: Public API
///----------------------------------------

///----------------------------------------
///   @brief Run the live-stacking polling loop for @p watch_dir.
/// @details Constructs a @ref LiveStackSession for the directory and runs its
///          polling loop until @c astap::esc_pressed is set.
///   @param watch_dir Directory to watch for incoming image files.
///----------------------------------------

void stack_live(const std::filesystem::path& watch_dir);

///----------------------------------------
/// MARK: LiveStackSession
///----------------------------------------

///----------------------------------------
///   @class LiveStackSession
///   @brief Per-run live-stacking state and polling loop.
/// @details Owns the per-run state that the original source kept as locals
///          inside @c stack_live and as unit-scope typed constants
///          (@c oldra0, @c olddec0, @c oldexposure, @c memo1_text). Exposed
///          for testability; normal callers should use @ref stack_live.
///----------------------------------------

class LiveStackSession final {
public:
    /// @brief Notified after each frame is processed (accepted or rejected).
    /// @param accepted Frames successfully aligned + accumulated.
    /// @param rejected Frames dropped (alignment failure, load failure, etc.).
    /// @param total Counter that resets on slew/exposure change.
    using FrameAddedHook = std::function<void(int accepted, int rejected, int total)>;

    /// @brief Notified with a human-readable status line.
    using MessageHook = std::function<void(const std::string&)>;

    ///----------------------------------------
    ///  @brief Construct a session that watches @p watch_dir.
    ///  @param watch_dir Directory to poll for incoming image files.
    ///----------------------------------------

    explicit LiveStackSession(std::filesystem::path watch_dir);

    /// @brief Install a hook called after every frame attempt.
    void set_frame_added_hook(FrameAddedHook hook) { frame_hook_ = std::move(hook); }

    /// @brief Install a hook called with log messages.
    void set_message_hook(MessageHook hook) { message_hook_ = std::move(hook); }

    ///----------------------------------------
    ///  @brief Main polling / stacking loop. Returns when @c astap::esc_pressed is set.
    ///----------------------------------------

    void run();

private:
    /// @brief Load + calibrate + align + accumulate a single frame.
    /// @return @c true if the frame was accepted into the stack.
    [[nodiscard]] bool process_frame(const std::filesystem::path& filename);

    /// @brief Forward to frame_hook_ if installed.
    void emit_frame_added();

    /// @brief Forward to message_hook_ if installed.
    void emit_message(const std::string& msg);

    /// @brief Reset per-stack accumulators.
    void reset_var() noexcept;
    
    /// @brief Update the running FITS header memo with current stacking metadata.
    void update_header();
    
    ///----------------------------------------
    ///      @brief Poll @p dir for the first ready image file.
    ///      @param dir Directory to scan.
    /// @param[out] out_file Receives the selected path on success.
    ///     @return @c true if a file was found and is readable.
    ///----------------------------------------
    
    [[nodiscard]] bool file_available(const std::filesystem::path& dir,
                                      std::filesystem::path& out_file) const;
                                      
    ///----------------------------------------
    ///  @brief Placeholder JPG export. TODO: wire a real JPEG encoder.
    ///  @param path Destination file path.
    ///  @param img Source image buffer.
    /// @return @c true if the placeholder file was written successfully.
    ///----------------------------------------
    
    [[nodiscard]] static bool save_as_jpg(const std::filesystem::path& path,
                                          const ImageArray& img);
                                          
    ///----------------------------------------
    ///  @brief Build a yyyymmdd_hhmmss timestamp for filename stamping.
    /// @return Formatted timestamp string.
    ///----------------------------------------
    
    [[nodiscard]] static std::string current_date_string();
    
    std::filesystem::path watch_dir_;
    
    // State formerly held as unit-scope typed constants.
    double old_ra0_      = 0.0;
    double old_dec0_     = -3.14159265358979323846 / 2.0;
    double old_exposure_ = 0.0;
    
    // Saved FITS header text from the first frame.
    std::string memo1_text_;
    
    // Per-stack accumulators.
    bool init_          = false;
    int  counter_       = 0;
    int  bad_counter_   = 0;
    int  total_counter_ = 0;
    int  spinner_       = 0;
    int  width_max_     = 0;
    int  height_max_    = 0;
    int  old_width_     = 0;
    int  old_height_    = 0;
    int  binning_       = 1;

    ImageArray img_average_;

    FrameAddedHook frame_hook_;
    MessageHook    message_hook_;
};
    
} // namespace
