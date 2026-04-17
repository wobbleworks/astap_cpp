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

#include "stack.h"
#include "stack_routines.h"
#include "../core/demosaic.h"
#include "../core/fits.h"
#include "../core/globals.h"
#include "../solving/astrometric_solving.h"
#include "../solving/star_align.h"

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
/// MARK: Engine imports
///----------------------------------------

using astap::solving::bin_and_find_stars;
using astap::solving::find_quads;
using astap::solving::find_offset_and_rotation;
using astap::solving::reset_solution_vectors;
using astap::solving::report_binning;
using astap::solving::quad_star_distances1;
using astap::solving::quad_star_distances2;

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

void LiveStackSession::emit_frame_added() {
    if (frame_hook_) {
        frame_hook_(counter_, bad_counter_, total_counter_);
    }
}

void LiveStackSession::emit_message(const std::string& msg) {
    if (message_hook_) {
        message_hook_(msg);
    }
}

bool LiveStackSession::process_frame(const fs::path& filename) {
    if (!astap::core::load_fits(filename, /*light=*/true, /*load_data=*/true,
                                /*update_memo=*/true, /*get_ext=*/0,
                                astap::memo1_lines, astap::head,
                                astap::img_loaded)) {
        emit_message("Error loading " + filename.string());
        return false;
    }

    // Detect mount slew via change in (ra0, dec0).
    const auto distance = ang_sep(astap::head.ra0, astap::head.dec0,
                                  old_ra0_, old_dec0_);
    old_ra0_  = astap::head.ra0;
    old_dec0_ = astap::head.dec0;
    if (distance > (0.2 * kPi / 180.0) && total_counter_ != 0) {
        emit_message("Mount slew detected — restarting stack.");
        reset_var();
    }

    // Exposure change resets the accumulator (mixing different exposures
    // would need weighted averaging; simpler to start fresh).
    if (total_counter_ != 0 && old_exposure_ != 0.0 &&
        std::abs(astap::head.exposure - old_exposure_) > 0.01) {
        emit_message("Exposure changed — restarting stack.");
        reset_var();
    }
    old_exposure_ = astap::head.exposure;

    (void)apply_dark_and_flat(astap::img_loaded, astap::head);

    // OSC demosaic for mono-with-bayer frames.
    if (astap::head.naxis3 == 1 && !astap::bayerpat.empty()) {
        const auto pattern = astap::core::get_demosaic_pattern(
            2, astap::xbayroff, astap::ybayroff, astap::roworder);
        astap::core::demosaic_bayer(astap::img_loaded, astap::head, pattern,
                                    astap::core::DemosaicMethod::Bilinear);
    }

    if (!init_) {
        astap::head_ref = astap::head;
        width_max_  = astap::head.width;
        height_max_ = astap::head.height;
        old_width_  = astap::head.width;
        old_height_ = astap::head.height;
        binning_    = report_binning(astap::head.height);
        img_average_.assign(astap::head.naxis3,
            std::vector<std::vector<float>>(astap::head.height,
                std::vector<float>(astap::head.width, 0.0f)));
    } else if (astap::head.width != old_width_ ||
               astap::head.height != old_height_) {
        emit_message("Size mismatch vs reference — skipping " +
                     filename.filename().string());
        ++bad_counter_;
        emit_frame_added();
        return false;
    }

    // Star detection.
    auto starlist = StarList{};
    auto warning = std::string{};
    bin_and_find_stars(astap::img_loaded, binning_, /*cropping=*/1.0,
                       /*hfd_min=*/std::max(0.8, astap::hfd_min_setting),
                       astap::max_stars_setting,
                       /*get_hist=*/true, starlist, warning);

    // Alignment.
    if (!init_) {
        find_quads(starlist, quad_star_distances1);
        reset_solution_vectors(1.0);
    } else {
        find_quads(starlist, quad_star_distances2);
        if (!find_offset_and_rotation(3, astap::quad_tolerance)) {
            emit_message("Not enough quad matches — skipping " +
                         filename.filename().string());
            ++bad_counter_;
            emit_frame_added();
            return false;
        }
    }
    init_ = true;

    ++counter_;
    ++total_counter_;

    // Running-average accumulation. calc_newx_newy maps (fitsX+1, fitsY+1)
    // through solution_vectorX/Y (vector-based branch) to reference pixel
    // space, leaving results in x_new_float / y_new_float (0-based).
    for (auto fitsY = 0; fitsY < astap::head.height; ++fitsY) {
        for (auto fitsX = 0; fitsX < astap::head.width; ++fitsX) {
            calc_newx_newy(/*vector_based=*/true,
                           static_cast<double>(fitsX + 1),
                           static_cast<double>(fitsY + 1));
            const auto xn = static_cast<int>(std::round(astap::x_new_float));
            const auto yn = static_cast<int>(std::round(astap::y_new_float));
            if (xn < 0 || xn >= width_max_ ||
                yn < 0 || yn >= height_max_) {
                continue;
            }
            for (auto col = 0; col < astap::head.naxis3; ++col) {
                auto& acc = img_average_[col][yn][xn];
                acc = (acc * (counter_ - 1) +
                       astap::img_loaded[col][fitsY][fitsX]) / counter_;
            }
        }
    }

    // Publish the running average to the viewer globals so the GUI can
    // refresh its canvas.
    astap::img_loaded = img_average_;
    astap::head = astap::head_ref;
    astap::head.light_count = counter_;

    emit_message("Added " + filename.filename().string() +
                 " — total " + std::to_string(counter_) + ".");
    emit_frame_added();
    return true;
}

void LiveStackSession::run() {
    astap::live_stacking.store(true);
    reset_var();
    astap::pause_pressed.store(false);
    astap::esc_pressed.store(false);
    total_counter_ = 0;
    spinner_       = 0;

    emit_message("Live stack started. Watching " + watch_dir_.string());

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
            spinner_ = (spinner_ + 1) % 3;
            continue;
        }
        waiting = false;

        (void)process_frame(filename);
        if (astap::esc_pressed.load()) {
            break;
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
            emit_message("Warning: could not rename " +
                         filename.filename().string());
        }
    }

    astap::live_stacking.store(false);
    emit_message("Live stack stopped. Accepted " +
                 std::to_string(counter_) + ", rejected " +
                 std::to_string(bad_counter_) + ".");
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
