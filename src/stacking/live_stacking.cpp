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
#include "../core/photometry.h"   // get_background (colour correction)
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

ColourCorrection colour_correction_factors(const ImageArray& img) {
    if (img.size() < 3) {
        return {};   // not a 3-colour image → identity, inactive (Pascal early exit)
    }
    // Measure each channel's background + star level. Red is measured last so the
    // module histogram is left describing red, matching the Pascal ordering.
    auto bR = astap::Background{};
    auto bG = astap::Background{};
    auto bB = astap::Background{};
    astap::core::get_background(1, img, /*calc_hist=*/true, /*calc_noise_level=*/true, bG);
    astap::core::get_background(2, img, /*calc_hist=*/true, /*calc_noise_level=*/true, bB);
    astap::core::get_background(0, img, /*calc_hist=*/true, /*calc_noise_level=*/true, bR);
    return colour_correction_from_backgrounds(bR, bG, bB);
}

void LiveStackSession::reset_var() noexcept {
    init_        = false;
    counter_     = 0;
    bad_counter_ = 0;
    // Reset the exposure / temperature / Julian-day accumulators that feed the
    // running stacked header (unit_live_stacking.pas:142-145).
    astap::sum_exp        = 0.0;
    astap::sum_temp       = 0.0;
    astap::jd_sum         = 0.0;
    astap::jd_start_first = 1e99;
    astap::jd_end_last    = 0.0;
    header_snapshot_.clear();
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
    if (counter_ <= 0) { return; }

    // Start from the first frame's header, then overwrite with the accumulated
    // stacking metadata (unit_live_stacking.pas:68-90).
    auto& memo = astap::memo1_lines;
    memo = header_snapshot_;

    astap::core::update_text(memo, "COMMENT 1",
        "  Written by Astrometric Stacking Program. www.hnsky.org");
    astap::core::update_text(memo, "HISTORY 1", "  Stacking method LIVE STACKING");
    astap::core::update_integer(memo, "EXPTIME =",
        " / Total luminance exposure time in seconds.      ",
        static_cast<int>(std::lround(astap::sum_exp)));
    astap::core::update_text(memo, "CALSTAT =", "'" + astap::head.calstat + "'");
    astap::core::update_text(memo, "DATE-OBS=",
        "'" + jd_to_date(astap::jd_start_first) + "'");
    const double jd_avg = astap::jd_sum / counter_;
    astap::core::update_generic(memo, "JD-AVG  ",
        std::format("{:.6f}", jd_avg),
        "Julian Day of the observation mid-point.       ");
    astap::head.date_avg = jd_to_date(jd_avg);
    astap::core::update_text(memo, "DATE-AVG=", "'" + astap::head.date_avg + "'");
    astap::core::update_integer(memo, "LIGH_CNT=",
        " / Light frames combined.                  ", counter_);
    astap::core::update_integer(memo, "DARK_CNT=",
        " / Darks used for luminance.               ", astap::head.dark_count);
    astap::core::update_integer(memo, "FLAT_CNT=",
        " / Flats used for luminance.               ", astap::head.flat_count);
    astap::core::update_integer(memo, "BIAS_CNT=",
        " / Flat-darks used for luminance.          ", astap::head.flatdark_count);
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

    // Detect mount slew via change in (ra0, dec0). A slew starts a fresh stack
    // AND skips this frame — the exposure spanning the slew may be trailed
    // (transition_image, unit_live_stacking.pas:181-219). A slew at the very
    // start (total_counter_ == 0) does not skip.
    bool transition_image = false;
    const auto distance = ang_sep(astap::head.ra0, astap::head.dec0,
                                  old_ra0_, old_dec0_);
    old_ra0_  = astap::head.ra0;
    old_dec0_ = astap::head.dec0;
    if (distance > (0.2 * kPi / 180.0) && total_counter_ != 0) {
        emit_message("Mount slew detected — restarting stack; skipping the "
                     "transition frame.");
        reset_var();
        transition_image = true;
    }

    // Exposure change resets the accumulator (mixing different exposures
    // would need weighted averaging; simpler to start fresh). Unlike a slew this
    // does not skip the frame (unit_live_stacking.pas:211-216).
    else if (total_counter_ != 0 && old_exposure_ != 0.0 &&
             std::abs(astap::head.exposure - old_exposure_) > 0.01) {
        emit_message("Exposure changed — restarting stack.");
        reset_var();
    }
    old_exposure_ = astap::head.exposure;

    if (transition_image) {
        return false;  // possibly slew-trailed — do not stack this frame
    }

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
        // Snapshot the first frame's header for the running stacked header
        // (Pascal saves memo1_text here, unit_live_stacking.pas:233).
        header_snapshot_ = astap::memo1_lines;
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
                       /*hfd_max=*/astap::hfd_max_setting,
                       astap::max_stars_setting,
                       /*get_hist=*/true, starlist, warning);

    // Alignment. Match the Pascal quad assignment (unit_live_stacking.pas:262,
    // 313): the reference (first) frame populates quad_star_distances2 and each
    // new source frame populates quad_star_distances1. find_offset_and_rotation
    // solves solution := quad2 -> quad1, so this assignment makes it map
    // reference -> source (Pascal's deliberate "inverse mapping"). The
    // accumulation loop below depends on that direction to gather each reference
    // pixel from the source pixel it came from.
    if (!init_) {
        find_quads(starlist, quad_star_distances2);
        reset_solution_vectors(1.0);
    } else {
        find_quads(starlist, quad_star_distances1);
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

    // On the first accepted (reference) frame, derive the colour-correction factors
    // from it (Pascal computes these at c==0 init, unit_live_stacking.pas:281-303).
    if (counter_ == 1 && colour_correction_) {
        cc_ = colour_correction_factors(astap::img_loaded);
        if (cc_.active) {
            emit_message("Colour correction: white-balancing green/blue to red.");
        }
    }

    // Accumulate exposure / temperature / Julian-day totals for the stacked
    // header (unit_live_stacking.pas:335-344). The live path tracks only the
    // start date and the midpoint sum (no jd_end_last, unlike the batch stackers).
    astap::sum_exp  += astap::head.exposure;
    astap::sum_temp += astap::head.set_temperature;
    date_to_jd(astap::head.date_obs, astap::head.date_avg, astap::head.exposure);
    astap::jd_start_first = std::min(astap::jd_start, astap::jd_start_first);
    astap::jd_sum        += astap::jd_mid;

    // Running-average accumulation via inverse mapping (unit_live_stacking.pas:353-410).
    // Iterate the reference-image grid; for each reference pixel find the source
    // pixel it came from. calc_newx_newy's vector-based branch applies
    // solution_vectorX/Y (reference -> source, per the quad assignment above),
    // taking a 1-based reference coordinate and leaving the 0-based source pixel
    // in x_new_float / y_new_float. Read that source pixel, write the reference
    // pixel. Gathering visits every reference pixel exactly once, so the running
    // mean is never double-written within a frame and no destination pixel is
    // left unwritten (zero) — the two defects of the previous forward scatter.
    for (auto fitsY = 0; fitsY < height_max_; ++fitsY) {
        for (auto fitsX = 0; fitsX < width_max_; ++fitsX) {
            calc_newx_newy(/*vector_based=*/true,
                           static_cast<double>(fitsX + 1),
                           static_cast<double>(fitsY + 1));
            const auto xn = static_cast<int>(std::round(astap::x_new_float));
            const auto yn = static_cast<int>(std::round(astap::y_new_float));
            // Bounds-check the SOURCE pixel against the source-frame dimensions
            // (Pascal width_maxS / height_maxS, unit_live_stacking.pas:338-339,360).
            if (xn < 0 || xn >= astap::head.width ||
                yn < 0 || yn >= astap::head.height) {
                continue;
            }
            for (auto col = 0; col < astap::head.naxis3; ++col) {
                auto& acc = img_average_[col][fitsY][fitsX];
                if (colour_correction_) {
                    // Pascal colour-correction path (unit_live_stacking.pas:381-406):
                    // skip zero (masked) source pixels, else (value + add)·multiply /
                    // largest, clamped >= 0, into the running average.
                    double dum = astap::img_loaded[col][yn][xn];
                    if (dum == 0.0) {
                        continue;   // 'if dum<>0' — leave the running average untouched
                    }
                    dum = (dum + cc_.add[col]) * cc_.multiply[col] / cc_.largest;
                    if (dum < 0.0) { dum = 0.0; }
                    acc = (acc * (counter_ - 1) + dum) / counter_;
                } else {
                    acc = (acc * (counter_ - 1) +
                           astap::img_loaded[col][yn][xn]) / counter_;
                }
            }
        }
    }

    // Publish the running average to the viewer globals so the GUI can
    // refresh its canvas.
    astap::img_loaded = img_average_;
    astap::head = astap::head_ref;
    astap::head.light_count = counter_;

    // Refresh the in-memory FITS header with the accumulated stacking metadata.
    update_header();

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
