///----------------------------------------
///      @file focus_fit.cpp
///   @ingroup ASTAP++
///     @brief Multi-frame focus-curve fit for the -focus1 / -focus2 CLI and
///            the GUI Focus dialog.
///   @details For each supplied FITS: load, read FOCUSPOS from the header,
///            run @ref astap::stacking::analyse_image to get the median HFD,
///            collect (position, hfd) pairs, and fit a hyperbola via
///            @ref astap::core::find_best_hyperbola_fit. The result includes
///            the predicted focuser position at best focus and the mean
///            residual.
///    @author Created by John Stephen on 4/21/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "focus_fit.h"

#include "../core/fits.h"
#include "../core/hyperbola.h"
#include "../core/image_io.h"
#include "../stacking/stack.h"

#include <string>
#include <vector>

///----------------------------------------
namespace astap::solving {
///----------------------------------------

FocusFitResult fit_focus_hyperbola(std::span<const std::filesystem::path> images,
                                    FocusFitSample* samples_out) {
    FocusFitResult out{};
    if (images.empty()) {
        out.message = "No images supplied.";
        return out;
    }

    auto points = std::vector<astap::core::FocusPoint>{};
    points.reserve(images.size());

    for (const auto& path : images) {
        auto head = astap::Header{};
        auto img  = astap::ImageArray{};
        auto memo = std::vector<std::string>{};
        auto fname = path.string();   // load_image takes std::string&
        if (!astap::core::load_image(fname, img, head, memo,
                                      /*re_center=*/false, /*plot=*/false)) {
            out.message = "Failed to load: " + path.string();
            return out;
        }

        auto bck = astap::Background{};
        auto star_counter = 0;
        auto hfd_median = 0.0;
        astap::stacking::analyse_image(img, head,
                                        /*snr_min=*/30.0,
                                        /*report_type=*/0,
                                        star_counter, bck, hfd_median);

        FocusFitSample sample{
            .path = path,
            .position = static_cast<double>(head.focus_pos),
            .hfd = hfd_median,
            .stars = star_counter,
        };
        if (samples_out != nullptr) {
            // Caller supplied a sink — append. (Using raw pointer+append so
            // the span-only overload stays simple for tests / CLI.)
        }

        // Skip frames with unusable HFD (analyse_image returns 99.0 when
        // detection found no stars).
        if (hfd_median <= 0.0 || hfd_median >= 98.0) {
            if (samples_out != nullptr) {
                // Persist the bad sample too so the GUI can show why it was
                // dropped. Caller is responsible for treating hfd ≥ 98 as
                // "rejected".
                *samples_out = sample;
                ++samples_out;
            }
            continue;
        }
        if (head.focus_pos == 0) {
            out.message = "FOCUSPOS keyword missing in " + path.string();
            if (samples_out != nullptr) {
                *samples_out = sample;
                ++samples_out;
            }
            continue;
        }

        points.push_back({.position = sample.position, .hfd = sample.hfd});
        if (samples_out != nullptr) {
            *samples_out = sample;
            ++samples_out;
        }
    }

    if (points.size() < 3) {
        out.message = "Need at least 3 images with FOCUSPOS + usable stars "
                       "(got " + std::to_string(points.size()) + ").";
        return out;
    }

    const auto fit = astap::core::find_best_hyperbola_fit(
        std::span<const astap::core::FocusPoint>{points});

    out.focus_best = fit.focus_position;
    out.lowest_error = fit.mean_error;
    out.a = fit.a;
    out.b = fit.b;
    out.samples_used = static_cast<int>(points.size());
    out.ok = true;
    out.message = "OK";
    return out;
}

} // namespace
