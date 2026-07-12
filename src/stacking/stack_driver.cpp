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

} // namespace

///----------------------------------------
/// MARK: stack_files
///----------------------------------------

StackResult stack_files(std::span<const std::filesystem::path> frames,
                        StackMethod method,
                        const std::filesystem::path& output,
                        const std::filesystem::path& master_dark,
                        const std::filesystem::path& master_flat) {
    StackResult r;
    r.frames_input = static_cast<int>(frames.size());
    if (frames.empty()) { r.message = "no input frames"; return r; }

    // Load the optional master calibration frames into the resident dark/flat
    // state. The combiner applies them per frame via apply_dark_and_flat.
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

    // ---- A + B: pre-solve pass -------------------------------------------
    // Load each frame; if it lacks a WCS, plate-solve it and persist the
    // solution back into the file. Frames that will not solve are dropped so
    // they never reach the combiner — the reference frame is the first survivor
    // and must carry WCS. Mirrors unit_stack.pas:13129-13190.
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

    // ---- C: select internal-astrometry alignment + dispatch --------------
    astap::use_astrometry_internal = true;
    astap::use_ephemeris_alignment = false;
    astap::use_manual_align        = false;
    const int process_as_osc = 0;
    int counter = 0;
    switch (method) {
        case StackMethod::Average:   stack_average(process_as_osc, keep, counter);   break;
        case StackMethod::SigmaClip: stack_sigmaclip(process_as_osc, keep, counter); break;
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

} // namespace astap::stacking
