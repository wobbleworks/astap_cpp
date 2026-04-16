///----------------------------------------
///      @file live_stacking.cpp
///   @ingroup ASTAP++
///     @brief Live-stacking implementation.
///   @details Polls a watch directory for newly captured frames, aligns each
///            onto the first accepted frame, and maintains a running average
///            image. Many collaborator modules (loader, dark/flat pipeline,
///            UI memos) are not yet ported; those integration points are
///            marked with TODO.
///    @author Ported from Han Kleijn's unit_live_stacking.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "live_stacking.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <format>
#include <fstream>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>

namespace fs = std::filesystem;

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

///----------------------------------------
/// MARK: File-local helpers
///----------------------------------------

namespace {

constexpr auto kPi = 3.14159265358979323846;

/// @brief Image extensions accepted by the watcher. Matched case-insensitively.
constexpr std::array<std::string_view, 21> kSupportedExts = {
    ".fit",  ".fits", ".png", ".jpg", ".bmp", ".tif", ".tiff", ".xisf",
    ".raw",  ".crw",  ".cr2", ".cr3", ".kdc", ".dcr", ".mrw",  ".arw",
    ".nef",  ".nrw",  ".dng", ".orf", ".raf"
};

///----------------------------------------
///  @brief Test whether a path's extension is in @ref kSupportedExts.
///  @param p Path to test.
/// @return @c true if the extension is supported.
///----------------------------------------

[[nodiscard]] bool ext_supported(const fs::path& p) {
    auto e = p.extension().string();
    std::transform(e.begin(), e.end(), e.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return std::find(kSupportedExts.begin(), kSupportedExts.end(),
                     std::string_view{e}) != kSupportedExts.end();
}

///----------------------------------------
///  @brief Test whether a file can currently be opened for reading.
/// @details Detects files still being written by the capture program.
///  @param p Path to test.
/// @return @c true if the file is openable.
///----------------------------------------

[[nodiscard]] bool file_readable(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return f.good();
}

///----------------------------------------
///  @brief Angular separation between two equatorial coordinates.
/// @details Re-implemented locally so this TU does not need to forward-declare
///          the star-align module.
///  @param ra1 First right ascension (radians).
///  @param dec1 First declination (radians).
///  @param ra2 Second right ascension (radians).
///  @param dec2 Second declination (radians).
/// @return Angular separation in radians.
///----------------------------------------

[[nodiscard]] double ang_sep(double ra1, double dec1, double ra2, double dec2) noexcept {
    const auto s = std::sin(dec1) * std::sin(dec2) +
                   std::cos(dec1) * std::cos(dec2) * std::cos(ra1 - ra2);
    return std::acos(std::clamp(s, -1.0, 1.0));
}
    
} // namespace

///----------------------------------------
/// MARK: External dependencies (TODO)
///----------------------------------------
// These live in astap_main / unit_stack / unit_star_align /
// unit_astrometric_solving in the original source and have not been ported
// yet. TODO markers below record the integration points so this file is
// self-contained syntactically.
//
// TODO(astap_main): port `head` (FITS header struct), `img_loaded`,
// `bayerpat`, `process_as_osc`, `sum_exp`, `sum_temp`, `jd_*`,
// `light_exposure`, `light_temperature`, `flat_filter`, `counterL`,
// `solution_vectorX/Y`, `nr_references*`, `quad_star_distances*`,
// `filename2`, etc.
//
// TODO(unit_stack): apply_dark_and_flat, analyse_listview,
// reset_solution_vectors, use_histogram, plot_fits, demosaic_bayer,
// test_bayer_matrix, report_binning, bin_and_find_stars, find_quads,
// find_offset_and_rotation, date_to_jd, JdToDate, memo2_message.
//
// TODO(ui): mainwindow.Memo1 / Image1, stackmenu1 controls
// (live_stacking1, live_stacking_pause1, write_jpeg1, etc.).

///----------------------------------------
/// MARK: LiveStackSession
///----------------------------------------

LiveStackSession::LiveStackSession(fs::path watch_dir) :
    
    // Initialize members
    watch_dir_(std::move(watch_dir)) {
}

void LiveStackSession::reset_var() noexcept {
    init_        = false;
    counter_     = 0;
    bad_counter_ = 0;
    // TODO(astap_main globals): also reset
    //   sum_exp = 0; sum_temp = 0;
    //   jd_sum = 0; jd_start_first = 1e99; jd_end_last = 0;
    //   light_exposure = 987654321; light_temperature = 987654321;
    //   flat_filter = "987654321";
}

bool LiveStackSession::file_available(const fs::path& dir,
                                      fs::path& out_file) const {
    // Bail out if the watch target is not a directory
    std::error_code ec;
    if (!fs::is_directory(dir, ec)) {
        return false;
    }
    
    // Scan for the first readable, supported image file
    for (const auto& entry : fs::directory_iterator(dir, ec)) {
        if (ec) {
            break;
        }
        if (!entry.is_regular_file(ec)) {
            continue;
        }
        const auto& p = entry.path();
        if (!ext_supported(p)) {
            continue;
        }
        if (!file_readable(p)) {
            continue;
        }
        out_file = p;
        return true;
    }
    return false;
}

std::string LiveStackSession::current_date_string() {
    using namespace std::chrono;
    const auto now = floor<seconds>(system_clock::now());
    // YYYYMMDD_HHMMSS — no separators inside the date or time portions.
    return std::format("{:%Y%m%d_%H%M%S}", now);
}

bool LiveStackSession::save_as_jpg(const fs::path& path,
                                   [[maybe_unused]] const ImageArray& img) {
    // TODO: real JPEG encoding. For now write a placeholder so callers can see
    // the snapshot was attempted.
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        return false;
    }
    constexpr auto kStub = std::string_view{
        "ASTAP++ live_stacking JPG export stub — TODO wire JPEG encoder.\n"};
    out.write(kStub.data(), static_cast<std::streamsize>(kStub.size()));
    return out.good();
}

void LiveStackSession::update_header() {
    // TODO: rewrite once the FITS header helpers (update_text/update_integer/
    // update_generic, JdToDate, head.*) are ported. The original mutates the
    // in-memory header text: COMMENT 1, HISTORY 1, EXPTIME, CALSTAT, DATE-OBS,
    // JD-AVG, DATE-AVG, LIGH_CNT, DARK_CNT, FLAT_CNT, BIAS_CNT.
}

void LiveStackSession::run() {
    // Enter the live-stacking state shared with the rest of the port
    astap::live_stacking.store(true);
    reset_var();
    astap::pause_pressed.store(false);
    astap::esc_pressed.store(false);
    total_counter_ = 0;
    spinner_       = 0;
    
    auto waiting = false;
    
    // TODO: prepare darks/flats — analyse_listview(...).
    // TODO: read colour_correction, hfd_min, max_stars from UI settings.
    
    while (!astap::esc_pressed.load()) {
        auto filename = fs::path{};
        const auto have_file =
            !astap::pause_pressed.load() && file_available(watch_dir_, filename);
            
        if (have_file) {
            waiting = false;
            auto transition_image = false;
            
            // TODO: load_image(filename, img_loaded, head, ...).
            // If load fails or esc was pressed, bail out.
            const auto loaded = false;  // placeholder
            if (astap::esc_pressed.load() || !loaded) {
                // TODO: memo2_message("Error loading file"); reset UI flags.
                astap::live_stacking.store(false);
                return;
            }
            
            // Detect mount slew via change in (ra0, dec0).
            // TODO: pull head.ra0/head.dec0 from real header.
            const auto ra0  = 0.0;
            const auto dec0 = 0.0;
            const auto distance = ang_sep(ra0, dec0, old_ra0_, old_dec0_);
            old_ra0_  = ra0;
            old_dec0_ = dec0;
            if (distance > (0.2 * kPi / 180.0)) {
                reset_var();
                if (total_counter_ != 0) {
                    transition_image = true;
                    // TODO: memo2_message("New telescope position ...").
                }
            } else {
                // TODO: detect head.exposure change vs old_exposure_; if
                // different, reset_var() and message.
            }
            // TODO: old_exposure_ = head.exposure;
            
            if (!transition_image) {
                if (!init_) {
                    // TODO: decide process_as_osc from head.naxis3, Xbinning,
                    // bayerpat and make_osc_color setting.
                    // TODO: memo1_text_ = mainwindow.Memo1.Text;
                }
                
                // TODO: apply_dark_and_flat(img_loaded, head);
                // TODO: memo2_message("Adding file: ...");
                if (astap::esc_pressed.load()) {
                    return;
                }
                
                if (!init_) {
                    // TODO: old_width_  = head.width;
                    //       old_height_ = head.height;
                } else {
                    // TODO: warn on size mismatch.
                }
                
                // TODO: if process_as_osc > 0, demosaic_bayer(img_loaded).
                
                if (!init_) {
                    // TODO: binning_ = report_binning(head.height);
                    // TODO: bin_and_find_stars(...);
                    // TODO: find_quads(starlist1, quad_star_distances1);
                    // TODO: setlength(img_average_, head.naxis3, H, W) and
                    //       zero-fill.
                }
                
                auto solution = true;
                if (init_) {
                    // TODO: bin_and_find_stars on second image, find_quads,
                    //       find_offset_and_rotation. On failure:
                    //         memo2_message("Not enough quad matches ...");
                    //         solution = false;
                } else {
                    // TODO: reset_solution_vectors(1);
                }
                init_ = true;
                
                if (solution) {
                    ++counter_;
                    ++total_counter_;
                    // TODO: sum_exp += head.exposure; sum_temp += ...;
                    // TODO: date_to_jd(...); update jd_start_first, jd_sum.
                    
                    // Affine alignment coefficients.
                    // TODO: pull from solution_vectorX/Y.
                    [[maybe_unused]] const auto aa = 1.0;
                    [[maybe_unused]] const auto bb = 0.0;
                    [[maybe_unused]] const auto cc = 0.0;
                    [[maybe_unused]] const auto dd = 0.0;
                    [[maybe_unused]] const auto ee = 1.0;
                    [[maybe_unused]] const auto ff = 0.0;
                    
                    // TODO: running-average accumulation into img_average_,
                    // either plain or with colour correction. Pseudocode:
                    //   for (y) for (x) {
                    //     int xn = round(aa*x + bb*y + cc);
                    //     int yn = round(dd*x + ee*y + ff);
                    //     if in-bounds:
                    //       for (col) img_average_[col][yn][xn] =
                    //         (img_average_[col][yn][xn] * (counter_-1)
                    //          + img_loaded[col][y][x]) / counter_;
                    //   }
                    
                    // TODO: head.cd1_1 = 0; head.height = height_max_;
                    //       head.width = width_max_;
                    //       img_loaded = img_average_;  // share buffer
                    // TODO: if (counter_ == 1) use_histogram(img_loaded, true);
                    // TODO: plot_fits(mainwindow.image1, false, false);
                    
                    // Snapshot exports.
                    // TODO: gate on write_jpeg setting.
                    [[maybe_unused]] const auto saved =
                        save_as_jpg(watch_dir_ / "stack.jpeg", img_average_);
                    // TODO: clipboard copy if interim_to_clipboard enabled.
                    // TODO: write log.txt if write_log enabled.
                } else {
                    ++bad_counter_;
                }
                // TODO: update files_live_stacked caption.
            }
            
            // Mark file as processed by renaming it. Appends "_@<date><ext>_"
            // or just "<ext>_" if already marked.
            const auto ext = filename.extension().string();
            const auto stem_path = filename.string().substr(
                0, filename.string().size() - ext.size());
            auto renamed = fs::path{};
            if (filename.string().find("_@") == std::string::npos) {
                renamed = stem_path + "_@" + current_date_string() + ext + "_";
            } else {
                renamed = stem_path + ext + "_";
            }
            std::error_code ec;
            fs::rename(filename, renamed, ec);
            if (ec) {
                std::fputc('\a', stderr);  // audible beep on failure
            }
        } else {
            // No new file (or paused). Show a simple status once, then sleep
            // and animate a tiny spinner.
            if (!waiting) {
                if (astap::pause_pressed.load() && counter_ > 0) {
                    // TODO: counterL = counter_; update_header();
                    // TODO: memo2_message("Live stack is suspended.");
                } else {
                    // TODO: memo2_message("Live stack is waiting for files.");
                }
            }
            waiting = true;
            
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
            
            // TODO: animate live-stacking caption ("  >", "> ", " > ").
            spinner_ = (spinner_ + 1) % 3;
        }
    }
    
    astap::live_stacking.store(false);
    // TODO: memo2_message("Live stack stopped. Save result if required");
    // TODO: counterL = counter_; if (counter_ > 0) update_header();
    memo1_text_.clear();
}

///----------------------------------------
/// MARK: Free-function entry point
///----------------------------------------

void stack_live(const fs::path& watch_dir) {
    auto session = LiveStackSession{watch_dir};
    session.run();
}
    
} // namespace
