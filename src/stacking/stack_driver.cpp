///----------------------------------------
///      @file stack_driver.cpp
///   @ingroup ASTAP++
///     @brief Headless batch image-stacking driver (internal-astrometry path).
///    @author Ported from Han Kleijn's unit_stack.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "stack_driver.h"

#include <cmath>
#include <format>
#include <string>
#include <vector>

#include "stack.h"            // update_solution_and_save, jd_to_date
#include "stack_routines.h"   // stack_average, stack_sigmaclip, FileToDo
#include "../core/fits.h"     // load_fits, save_fits, header-edit helpers
#include "../core/globals.h"  // head, img_loaded, memo1_lines, accumulators, flags
#include "../types.h"

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

// Defined in stack_routines.cpp (astap::stacking) → forwards to the core sink.
void memo2_message(std::string_view msg);

namespace {

///----------------------------------------
///  @brief Emit the ASTAP stacked-header keyword block for a mono stack.
///  @details Mirrors the monochrome branch of Tstackmenu1.stack_button1Click
///           (unit_stack.pas:13716-13804): it strips the per-frame interim
///           keywords and writes the aggregate EXPTIME / CALSTAT / PEDESTAL /
///           DATE-OBS/END/AVG / AIRMASS / JD-AVG / LUM_* block. The frame-loop
///           accumulators (jd_sum, jd_start_first, jd_end_last, airmass_sum) are
///           populated by the combiner; @p counterL is its combined-frame count.
///----------------------------------------

void write_stacked_header(StackMethod method, int counterL) {
    auto& memo = astap::memo1_lines;
    auto& head = astap::head;

    // Interim per-frame keywords are replaced by the aggregate ones below.
    astap::core::remove_key(memo, "DATE    ", false);
    astap::core::remove_key(memo, "EXPOSURE", false);
    astap::core::remove_key(memo, "CCD-TEMP", false);
    astap::core::remove_key(memo, "SET-TEMP", false);
    astap::core::remove_key(memo, "LIGH_CNT", false);
    astap::core::remove_key(memo, "DARK_CNT", false);
    astap::core::remove_key(memo, "FLAT_CNT", false);
    astap::core::remove_key(memo, "BIAS_CNT", false);

    astap::core::update_text(memo, "COMMENT 1", "  Written by ASTAP. www.hnsky.org");

    head.calstat += "S";  // status: stacked
    astap::core::update_integer(memo, "PEDESTAL=",
        " / Value added during calibration or stacking     ",
        static_cast<int>(std::lround(head.pedestal)));
    astap::core::update_text(memo, "CALSTAT =", "'" + head.calstat + "'");
    astap::core::add_text(memo, "ISSUES  =", "'" + head.issues + "'");

    head.date_obs = jd_to_date(astap::jd_start_first);
    astap::core::update_text(memo, "DATE-OBS=",
        "'" + head.date_obs + "'/ Date and time of the start of the observation.");
    astap::core::update_text(memo, "DATE-END=",
        "'" + jd_to_date(astap::jd_end_last) + "'/ Date and time of the end of the observation.");

    const int exposureL    = static_cast<int>(std::lround(head.exposure));
    const int temperatureL = head.set_temperature;

    if (counterL > 0) {
        astap::jd_mid  = astap::jd_sum / counterL;       // average mid-point
        astap::airmass = astap::airmass_sum / counterL;  // average airmass
        astap::core::update_generic(memo, "AIRMASS ",
            std::format("{:.4f}", astap::airmass),
            "Average relative optical path.                ");
        astap::core::update_generic(memo, "JD-AVG  ",
            std::format("{:.6f}", astap::jd_mid),
            "Julian Day of the observation mid-point.       ");
        head.date_avg = jd_to_date(astap::jd_mid);
        astap::core::update_text(memo, "DATE-AVG=", "'" + head.date_avg + "'");
        // Estimate DATE-OBS as the mid-point minus half the per-frame exposure.
        head.date_obs = jd_to_date(astap::jd_mid
                                   - head.exposure / (2.0 * 24.0 * 60.0 * 60.0));
    }

    astap::core::update_text(memo, "HISTORY 1",
        method == StackMethod::SigmaClip ? "  Stacking method SIGMA CLIP AVERAGE"
      : method == StackMethod::Mosaic    ? "  Stacking method MOSAIC"
                                         : "  Stacking method AVERAGE");
    astap::core::update_text(memo, "HISTORY 2", "  Processed as gray scale images.");

    astap::core::update_integer(memo, "EXPTIME =",
        " / Total exposure time in seconds.      ",
        static_cast<int>(std::lround(static_cast<double>(counterL) * exposureL)));
    astap::core::update_integer(memo, "SET-TEMP=",
        " / Average set temperature used for luminance.    ", temperatureL);
    astap::core::add_integer(memo, "LUM_EXP =",
        " / Average luminance exposure time.               ", exposureL);
    astap::core::add_integer(memo, "LUM_CNT =",
        " / Luminance images combined.                     ", counterL);
    astap::core::add_integer(memo, "LUM_DARK=",
        " / Darks used for luminance.                      ", head.dark_count);
    astap::core::add_integer(memo, "LUM_FLAT=",
        " / Flats used for luminance.                      ", head.flat_count);
    astap::core::add_integer(memo, "LUM_BIAS=",
        " / Flat-darks used for luminance.                 ", head.flatdark_count);
}

///----------------------------------------
///  @brief Mean-combine raw calibration frames into a mono image.
///  @details Loads each frame into the resident head/memo state, sums it (colour
///           input is reduced to mono via (R+G+B)/3), and divides by the count.
///           Frames whose dimensions differ from the first are skipped. Ports the
///           core of `average` (`unit_stack.pas:4814`); `temperature_avg` returns
///           the mean set-temperature over the combined frames.
///  @return Number of frames actually combined.
///----------------------------------------

[[nodiscard]] int average_frames(std::span<const std::filesystem::path> frames,
                                 ImageArray& out_img, double& temperature_avg) {
    out_img.clear();
    temperature_avg = 0.0;
    int count = 0, w = 0, h = 0;
    ImageArray tmp;
    for (const auto& path : frames) {
        if (!astap::core::load_fits(path, /*light=*/false, /*load_data=*/true,
                                    /*update_memo=*/true, /*get_ext=*/0,
                                    astap::memo1_lines, astap::head, tmp)) {
            memo2_message("Skipping unreadable frame: " + path.string());
            continue;
        }
        if (count == 0) {
            w = astap::head.width;
            h = astap::head.height;
            out_img.assign(1, std::vector<std::vector<float>>(
                                  h, std::vector<float>(w, 0.0f)));
        } else if (astap::head.width != w || astap::head.height != h) {
            // The Pascal test is an OR (a latent typo that risks an out-of-bounds
            // read); require both dimensions to match.
            memo2_message("Skipping frame with mismatched dimensions: " + path.string());
            continue;
        }
        const int nax3 = astap::head.naxis3;
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                out_img[0][y][x] += (nax3 == 3)
                    ? (tmp[0][y][x] + tmp[1][y][x] + tmp[2][y][x]) / 3.0f
                    : tmp[0][y][x];
            }
        }
        temperature_avg += astap::head.set_temperature;
        ++count;
    }
    if (count > 1) {
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                out_img[0][y][x] /= static_cast<float>(count);
            }
        }
    }
    if (count > 0) { temperature_avg /= count; }
    return count;
}

///----------------------------------------
///  @brief Emit the ASTAP stacked-header keyword block for a colour (LRGB) stack.
///  @details Mirrors the LRGB branch of Tstackmenu1.stack_button1Click
///           (unit_stack.pas:13717-13904): it strips the per-frame interim
///           keywords, marks the frame stacked/colour (NAXIS/NAXIS3=3, FILTER
///           wiped) and writes the per-channel LUM_/RED_/GRN_/BLU_ blocks plus
///           the total-exposure COMMENT lines from the counter globals the
///           combiner populated. Unlike the mono block there is no JD-AVG /
///           AIRMASS pass — that is written only for a mono (NAXIS3=1) stack.
///----------------------------------------

void write_lrgb_header() {
    auto& memo = astap::memo1_lines;
    auto& head = astap::head;

    // Interim per-frame keywords are replaced by the per-channel ones below.
    astap::core::remove_key(memo, "DATE    ", false);
    astap::core::remove_key(memo, "EXPOSURE", false);
    astap::core::remove_key(memo, "CCD-TEMP", false);
    astap::core::remove_key(memo, "SET-TEMP", false);
    astap::core::remove_key(memo, "LIGH_CNT", false);
    astap::core::remove_key(memo, "DARK_CNT", false);
    astap::core::remove_key(memo, "FLAT_CNT", false);
    astap::core::remove_key(memo, "BIAS_CNT", false);

    astap::core::update_text(memo, "COMMENT 1", "  Written by ASTAP. www.hnsky.org");

    head.calstat += "S";  // status: stacked
    astap::core::update_integer(memo, "PEDESTAL=",
        " / Value added during calibration or stacking     ",
        static_cast<int>(std::lround(head.pedestal)));
    astap::core::update_text(memo, "CALSTAT =", "'" + head.calstat + "'");
    astap::core::add_text(memo, "ISSUES  =", "'" + head.issues + "'");

    head.date_obs = jd_to_date(astap::jd_start_first);
    astap::core::update_text(memo, "DATE-OBS=",
        "'" + head.date_obs + "'/ Date and time of the start of the observation.");
    astap::core::update_text(memo, "DATE-END=",
        "'" + jd_to_date(astap::jd_end_last) + "'/ Date and time of the end of the observation.");

    astap::core::update_text(memo, "HISTORY 1", "  Stacking method AVERAGE");
    astap::core::update_text(memo, "HISTORY 2", "  Combined to colour image.");

    // Mark the header as a three-plane colour image and wipe the filter name
    // (the assembled RGB no longer belongs to a single filter).
    head.naxis  = 3;
    head.naxis3 = 3;
    astap::core::update_text(memo, "FILTER  =", "'        '");

    if (astap::counterL > 0) {
        astap::core::update_integer(memo, "EXPTIME =",
            " / Total exposure time in seconds.      ",
            static_cast<int>(std::lround(static_cast<double>(astap::counterL) * astap::exposureL)));
        astap::core::add_integer(memo, "LUM_EXP =", " / Luminance exposure time.                       ", astap::exposureL);
        astap::core::add_integer(memo, "LUM_CNT =", " / Luminance images combined.                     ", astap::counterL);
        astap::core::add_integer(memo, "LUM_DARK=", " / Darks used for luminance.                      ", astap::counterLdark);
        astap::core::add_integer(memo, "LUM_FLAT=", " / Flats used for luminance.                      ", astap::counterLflat);
        astap::core::add_integer(memo, "LUM_BIAS=", " / Flat-darks used for luminance.                 ", astap::counterLbias);
        astap::core::add_integer(memo, "LUM_TEMP=", " / Average set temperature used for luminance.    ", static_cast<int>(std::lround(astap::temperatureL)));
    }
    if (astap::counterR > 0) {
        astap::core::add_integer(memo, "RED_EXP =", " / Red exposure time.                             ", astap::exposureR);
        astap::core::add_integer(memo, "RED_CNT =", " / Red filter images combined.                    ", astap::counterR);
        astap::core::add_integer(memo, "RED_DARK=", " / Darks used for red.                            ", astap::counterRdark);
        astap::core::add_integer(memo, "RED_FLAT=", " / Flats used for red.                            ", astap::counterRflat);
        astap::core::add_integer(memo, "RED_BIAS=", " / Flat-darks used for red.                       ", astap::counterRbias);
        astap::core::add_integer(memo, "RED_TEMP=", " / Set temperature used for red.                  ", static_cast<int>(std::lround(astap::temperatureR)));
    }
    if (astap::counterG > 0) {
        astap::core::add_integer(memo, "GRN_EXP =", " / Green exposure time.                           ", astap::exposureG);
        astap::core::add_integer(memo, "GRN_CNT =", " / Green filter images combined.                  ", astap::counterG);
        astap::core::add_integer(memo, "GRN_DARK=", " / Darks used for green.                          ", astap::counterGdark);
        astap::core::add_integer(memo, "GRN_FLAT=", " / Flats used for green.                          ", astap::counterGflat);
        astap::core::add_integer(memo, "GRN_BIAS=", " / Flat-darks used for green.                     ", astap::counterGbias);
        astap::core::add_integer(memo, "GRN_TEMP=", " / Set temperature used for green.                ", static_cast<int>(std::lround(astap::temperatureG)));
    }
    if (astap::counterB > 0) {
        astap::core::add_integer(memo, "BLU_EXP =", " / Blue exposure time.                            ", astap::exposureB);
        astap::core::add_integer(memo, "BLU_CNT =", " / Blue filter images combined.                   ", astap::counterB);
        astap::core::add_integer(memo, "BLU_DARK=", " / Darks used for blue.                           ", astap::counterBdark);
        astap::core::add_integer(memo, "BLU_FLAT=", " / Flats used for blue.                           ", astap::counterBflat);
        astap::core::add_integer(memo, "BLU_BIAS=", " / Flat-darks used for blue.                      ", astap::counterBbias);
        astap::core::add_integer(memo, "BLU_TEMP=", " / Set temperature used for blue.                 ", static_cast<int>(std::lround(astap::temperatureB)));
    }

    if (astap::counterL > 0) astap::core::add_text(memo, "COMMENT 2",
        "  Total luminance exposure " + std::to_string(std::lround(static_cast<double>(astap::counterL) * astap::exposureL)));
    if (astap::counterR > 0) astap::core::add_text(memo, "COMMENT 3",
        "  Total red exposure " + std::to_string(std::lround(static_cast<double>(astap::counterR) * astap::exposureR)));
    if (astap::counterG > 0) astap::core::add_text(memo, "COMMENT 4",
        "  Total green exposure " + std::to_string(std::lround(static_cast<double>(astap::counterG) * astap::exposureG)));
    if (astap::counterB > 0) astap::core::add_text(memo, "COMMENT 5",
        "  Total blue exposure " + std::to_string(std::lround(static_cast<double>(astap::counterB) * astap::exposureB)));
}

///----------------------------------------
///  @brief Load the optional master dark/flat into the resident calibration
///         state so every combiner applies them per frame via apply_dark_and_flat.
///----------------------------------------

void load_masters(const std::filesystem::path& master_dark,
                  const std::filesystem::path& master_flat) {
    if (!master_dark.empty()) {
        MasterFrameInfo info{};
        if (set_master_dark(master_dark, info)) {
            memo2_message("Master dark loaded: " + master_dark.string());
        } else {
            memo2_message("Warning: could not load master dark " + master_dark.string());
        }
    }
    if (!master_flat.empty()) {
        MasterFrameInfo info{};
        if (set_master_flat(master_flat, info)) {
            memo2_message("Master flat loaded: " + master_flat.string());
        } else {
            memo2_message("Warning: could not load master flat " + master_flat.string());
        }
    }
}

} // namespace

///----------------------------------------
/// MARK: create_master
///----------------------------------------

MasterResult create_master(std::span<const std::filesystem::path> frames,
                           MasterKind kind,
                           const std::filesystem::path& output,
                           std::span<const std::filesystem::path> flatdark_frames) {
    MasterResult r;
    r.frames_input = static_cast<int>(frames.size());
    if (frames.empty()) { r.message = "no input frames"; return r; }

    // For a master flat, mean-combine the flat-darks/bias FIRST (into img_bias),
    // then the flats — so the resident header ends up carrying the flat's cards,
    // not the flat-dark's (order matters; unit_stack.pas:12428-12450).
    ImageArray img_bias;
    int bias_count = 0;
    const bool do_flatdark = (kind == MasterKind::Flat) && !flatdark_frames.empty();
    if (do_flatdark) {
        double bias_temp = 0.0;
        bias_count = average_frames(flatdark_frames, img_bias, bias_temp);
        if (bias_count == 0) {
            memo2_message("Warning: no flat-dark could be averaged; skipping bias subtraction.");
        }
    }

    ImageArray img;
    double temperature_avg = 0.0;
    const int count = average_frames(frames, img, temperature_avg);
    r.frames_combined = count;
    if (count == 0) { r.message = "no frame could be averaged"; return r; }

    // The averaged master is mono; the header carried the last frame's cards.
    auto& memo = astap::memo1_lines;
    auto& head = astap::head;
    head.naxis  = 2;
    head.naxis3 = 1;

    // Subtract the combined flat-dark from the combined flat (both mono).
    if (do_flatdark && bias_count > 0) {
        const int w = head.width, h = head.height;
        if (static_cast<int>(img_bias[0].size()) == h
            && static_cast<int>(img_bias[0][0].size()) == w) {
            for (int y = 0; y < h; ++y) {
                for (int x = 0; x < w; ++x) {
                    img[0][y][x] -= img_bias[0][y][x];
                }
            }
            head.flatdark_count = bias_count;
            head.calstat += "B";
            r.flatdarks_combined = bias_count;
            memo2_message("Applied the combined flat-dark on the combined flat.");
        } else {
            memo2_message("Warning: flat and flat-dark dimensions differ; "
                          "skipping bias subtraction.");
        }
    }

    const std::string_view count_key =
        kind == MasterKind::Dark ? "DARK_CNT="
      : kind == MasterKind::Flat ? "FLAT_CNT=" : "BIAS_CNT=";
    const std::string_view count_cmt =
        kind == MasterKind::Dark ? " / Number of dark images combined                "
      : kind == MasterKind::Flat ? " / Number of flat images combined                "
                                 : " / Number of bias images combined                ";

    astap::core::update_text(memo, "COMMENT 1", "  Written by ASTAP. www.hnsky.org");
    astap::core::update_integer(memo, count_key, count_cmt, count);
    astap::core::update_integer(memo, "CCD-TEMP=",
        " / Average sensor temperature (Celsius)           ",
        static_cast<int>(std::lround(temperature_avg)));
    if (r.flatdarks_combined > 0) {
        astap::core::update_integer(memo, "BIAS_CNT=",
            " / Number of flat-dark or bias images combined.   ", r.flatdarks_combined);
        astap::core::update_text(memo, "CALSTAT =", "'" + head.calstat + "'");
    }

    r.ok = astap::core::save_fits(img, memo, output, /*type1=*/-32, /*override2=*/true);
    r.output  = output;
    r.message = r.ok ? "ok" : "save failed";
    return r;
}

///----------------------------------------
/// MARK: stack_files
///----------------------------------------

StackResult stack_files(std::span<const std::filesystem::path> frames,
                        StackMethod method,
                        const std::filesystem::path& output,
                        const std::filesystem::path& master_dark,
                        const std::filesystem::path& master_flat,
                        double outlier_sigma) {
    StackResult r;
    r.frames_input = static_cast<int>(frames.size());
    if (frames.empty()) { r.message = "no input frames"; return r; }

    // The combiner applies these per frame via apply_dark_and_flat.
    load_masters(master_dark, master_flat);

    // ---- A + B: pre-solve pass -------------------------------------------
    // Load each frame; if it lacks a WCS, plate-solve it and persist the
    // solution back into the file. Frames that will not solve are dropped so
    // they never reach the combiner. While each frame is loaded, measure its
    // quality (star_detections / HFD^2) so the sharpest can be chosen as the
    // alignment reference. Mirrors unit_stack.pas:13129-13190 + 12821.
    struct Candidate { std::string name; double quality; int width; double hfd; };
    std::vector<Candidate> cands;
    cands.reserve(frames.size());
    for (const auto& path : frames) {
        if (!astap::core::load_fits(path, /*light=*/true, /*load_data=*/true,
                                    /*update_memo=*/true, /*get_ext=*/0,
                                    astap::memo1_lines, astap::head,
                                    astap::img_loaded)) {
            memo2_message("Skipping unreadable frame: " + path.string());
            continue;
        }
        bool solved = astap::head.cd1_1 != 0.0;
        if (!solved) {
            solved = update_solution_and_save(astap::img_loaded, astap::head,
                                              astap::memo1_lines, path);
        }
        if (!solved) {
            memo2_message("Dropping unsolved frame: " + path.string());
            continue;
        }
        int stars = 0;
        double hfd_median = 0.0;
        astap::Background bck{};
        analyse_image(astap::img_loaded, astap::head, /*snr_min=*/30.0,
                      /*report_type=*/0, stars, bck, hfd_median, nullptr);
        const double quality = (hfd_median > 0.0)
            ? static_cast<double>(stars) / (hfd_median * hfd_median) : 0.0;
        cands.push_back({path.string(), quality, astap::head.width, hfd_median});
    }
    r.frames_solved = static_cast<int>(cands.size());
    if (cands.empty()) { r.message = "no frame could be solved for alignment"; return r; }

    // ---- B': quality-outlier rejection (list_remove_outliers, :1600) ------
    // Drop frames whose quality sits more than outlier_sigma standard deviations
    // below the group mean (and any with a runaway HFD>90). Disabled when
    // outlier_sigma<=0. Runs before the reference pick so a rejected frame can
    // never become the reference.
    if (outlier_sigma > 0.0 && cands.size() > 2) {
        double sum = 0.0;
        int good = 0;
        for (const auto& c : cands) {
            if (c.hfd > 90.0) continue;
            sum += c.quality;
            ++good;
        }
        if (good > 0) {
            const double mean = sum / good;
            double var = 0.0;
            for (const auto& c : cands) {
                if (c.hfd > 90.0) continue;
                var += (mean - c.quality) * (mean - c.quality);
            }
            const double sd = std::sqrt(var / good);
            memo2_message("Analysing " + std::to_string(good)
                + " frames for outliers. Average quality (star_detections/sqr(hfd))="
                + std::format("{:.0f}", mean) + ", sigma=" + std::format("{:.1f}", sd));
            std::vector<Candidate> kept;
            kept.reserve(cands.size());
            for (auto& c : cands) {
                const bool reject = c.hfd > 90.0
                                 || (mean - c.quality) > outlier_sigma * sd;
                if (reject) {
                    memo2_message("Dropping low-quality outlier: " + c.name);
                } else {
                    kept.push_back(std::move(c));
                }
            }
            cands.swap(kept);
        }
    }
    r.frames_kept = static_cast<int>(cands.size());
    if (cands.empty()) { r.message = "all frames rejected as outliers"; return r; }

    // Put the best-quality frame first — the combiner takes the first frame as
    // the alignment reference (put_best_quality_on_top, unit_stack.pas:12821):
    // prefer a larger image, then the highest star_detections / HFD^2.
    {
        std::size_t best = 0;
        int    best_width   = cands[0].width;
        double best_quality = cands[0].quality;
        for (std::size_t i = 1; i < cands.size(); ++i) {
            if (cands[i].width > best_width
                || (cands[i].width == best_width && cands[i].quality > best_quality)) {
                best_width   = cands[i].width;
                best_quality = cands[i].quality;
                best         = i;
            }
        }
        if (best != 0) { std::swap(cands[0], cands[best]); }
        r.reference = cands[0].name;
        memo2_message("Reference image selected on quality "
                      "(star_detections/sqr(hfd)): " + cands[0].name);
    }

    std::vector<FileToDo> keep;
    keep.reserve(cands.size());
    for (const auto& c : cands) { keep.push_back(FileToDo{c.name, -1}); }

    // ---- C: select internal-astrometry alignment + dispatch --------------
    astap::use_astrometry_internal = true;
    astap::use_ephemeris_alignment = false;
    astap::use_manual_align        = false;
    const int process_as_osc = 0;
    int counter = 0;
    switch (method) {
        case StackMethod::Average:   stack_average(process_as_osc, keep, counter);   break;
        case StackMethod::SigmaClip: stack_sigmaclip(process_as_osc, keep, counter); break;
        case StackMethod::Mosaic:
            // Background-difference bound left at 0 (no tile bg-equalisation) —
            // the headless caller supplies pre-matched frames.
            stack_mosaic(process_as_osc, keep, /*max_dev_backgr=*/0.0, counter);
            break;
    }
    r.frames_combined = counter;
    if (counter == 0) { r.message = "combiner produced no output"; return r; }

    // ---- D: stacked-header keyword block + save --------------------------
    write_stacked_header(method, counter);
    r.ok = astap::core::save_fits(astap::img_loaded, astap::memo1_lines,
                                  output, /*type1=*/-32, /*override2=*/true);
    r.output  = output;
    r.message = r.ok ? "ok" : "save failed";
    return r;
}

///----------------------------------------
/// MARK: calibrate_align_files
///----------------------------------------

CalibrateResult calibrate_align_files(std::span<const std::filesystem::path> frames,
                                      const std::filesystem::path& master_dark,
                                      const std::filesystem::path& master_flat) {
    CalibrateResult r;
    r.frames_input = static_cast<int>(frames.size());
    if (frames.empty()) { r.message = "no input frames"; return r; }

    load_masters(master_dark, master_flat);

    // Pre-solve every frame and persist its WCS: calibration_and_alignment aligns
    // via astrometric_to_vector, which needs each frame to carry a solution.
    std::vector<FileToDo> keep;
    keep.reserve(frames.size());
    for (const auto& path : frames) {
        if (!astap::core::load_fits(path, /*light=*/true, /*load_data=*/true,
                                    /*update_memo=*/true, /*get_ext=*/0,
                                    astap::memo1_lines, astap::head,
                                    astap::img_loaded)) {
            memo2_message("Skipping unreadable frame: " + path.string());
            continue;
        }
        bool solved = astap::head.cd1_1 != 0.0;
        if (!solved) {
            solved = update_solution_and_save(astap::img_loaded, astap::head,
                                              astap::memo1_lines, path);
        }
        if (!solved) {
            memo2_message("Dropping unsolved frame: " + path.string());
            continue;
        }
        keep.push_back(FileToDo{path.string(), -1});
    }
    r.frames_solved = static_cast<int>(keep.size());
    if (keep.empty()) { r.message = "no frame could be solved for alignment"; return r; }

    astap::use_astrometry_internal = true;
    astap::use_ephemeris_alignment = false;
    astap::use_manual_align        = false;
    int counter = 0;
    calibration_and_alignment(/*process_as_osc=*/0, keep, counter);
    r.frames_aligned = counter;

    // The routine writes "<stem>_aligned.fit" beside each frame it kept; a frame
    // that failed alignment has had its name cleared in `keep`.
    for (const auto& f : keep) {
        if (f.name.empty()) { continue; }
        std::filesystem::path p(f.name);
        r.outputs.push_back(p.parent_path() / (p.stem().string() + "_aligned.fit"));
    }
    r.ok = counter > 0;
    r.message = r.ok ? "ok" : "no frame aligned";
    return r;
}

///----------------------------------------
/// MARK: combine_lrgb
///----------------------------------------

LrgbResult combine_lrgb(const LrgbInputs& inputs,
                        const std::filesystem::path& output) {
    LrgbResult r;
    const int colours = static_cast<int>(!inputs.red.empty())
                      + static_cast<int>(!inputs.green.empty())
                      + static_cast<int>(!inputs.blue.empty())
                      + static_cast<int>(!inputs.rgb.empty());
    r.channels_input = colours;
    r.has_luminance  = !inputs.luminance.empty();
    if (colours < 2) {
        r.message = "LRGB stacking needs at least two colour channels";
        return r;
    }

    // Reset the per-channel accumulators the combiner fills from each channel's
    // header (LIGH_CNT etc.) and the LRGB header block later reads. The GUI does
    // this once per object at unit_stack.pas:13288-13312; stack_LRGB itself does
    // not, so a stale count from an earlier stack would otherwise leak a bogus
    // block into a channel that was never supplied.
    astap::counterR = astap::counterG = astap::counterB = astap::counterL = 0;
    astap::counterRGB = 0;
    astap::counterRdark = astap::counterGdark = astap::counterBdark = astap::counterLdark = 0;
    astap::counterRflat = astap::counterGflat = astap::counterBflat = astap::counterLflat = 0;
    astap::counterRbias = astap::counterGbias = astap::counterBbias = astap::counterLbias = 0;
    astap::exposureR = astap::exposureG = astap::exposureB = astap::exposureL = astap::exposureRGB = 0;
    astap::temperatureR = astap::temperatureG = astap::temperatureB = astap::temperatureL = 0.0;

    // The reference (slot 0, used only to fix the alignment grid) is the
    // luminance if present, else the first available colour. Matches
    // unit_stack.pas:13596-13598.
    std::filesystem::path reference = inputs.luminance;
    if (reference.empty()) reference = inputs.red;
    if (reference.empty()) reference = inputs.green;
    if (reference.empty()) reference = inputs.blue;
    if (reference.empty()) reference = inputs.rgb;

    // files_to_process = [reference, R, G, B, RGB, L]; the combiner skips the
    // reference slot (c=0) when adding, so a colour re-used as the reference
    // still contributes through its own slot.
    std::vector<FileToDo> files(6);
    files[0] = FileToDo{reference.string(),         -1};
    files[1] = FileToDo{inputs.red.string(),        -1};
    files[2] = FileToDo{inputs.green.string(),      -1};
    files[3] = FileToDo{inputs.blue.string(),       -1};
    files[4] = FileToDo{inputs.rgb.string(),        -1};
    files[5] = FileToDo{inputs.luminance.string(),  -1};

    astap::use_astrometry_internal = true;
    astap::use_ephemeris_alignment = false;
    astap::use_manual_align        = false;

    int counter_colours = 0;
    stack_LRGB(files, counter_colours);
    r.channels_combined = counter_colours;
    r.reference         = reference;
    if (counter_colours == 0) {
        r.message = "colour combine produced no output";
        return r;
    }

    write_lrgb_header();
    r.ok = astap::core::save_fits(astap::img_loaded, astap::memo1_lines,
                                  output, /*type1=*/-32, /*override2=*/true);
    r.output  = output;
    r.message = r.ok ? "ok" : "save failed";
    return r;
}

} // namespace astap::stacking
