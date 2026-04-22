///----------------------------------------
///      @file focus_fit.h
///   @ingroup ASTAP++
///     @brief Multi-frame focus-curve fit.
///    @author Created by John Stephen on 4/21/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <filesystem>
#include <span>
#include <string>

///----------------------------------------
namespace astap::solving {
///----------------------------------------

///----------------------------------------
/// @brief Per-frame focus measurement emitted by @ref fit_focus_hyperbola.
///----------------------------------------

struct FocusFitSample {
    std::filesystem::path path;
    double position{};       ///< FOCUSPOS from the FITS header.
    double hfd{};            ///< Median HFD from analyse_image (≥98 means rejected).
    int    stars{};          ///< Number of stars contributing to the median.
};

///----------------------------------------
/// @brief Result of a multi-frame hyperbola focus-curve fit.
///----------------------------------------

struct FocusFitResult {
    bool   ok{false};        ///< True when a fit was produced.
    double focus_best{};     ///< Focuser position at best focus.
    double lowest_error{};   ///< Mean relative-HFD residual after convergence.
    double a{};              ///< Hyperbola a parameter (min HFD).
    double b{};              ///< Hyperbola b parameter (asymptote slope).
    int    samples_used{};   ///< How many frames contributed valid (position, hfd) pairs.
    std::string message;     ///< Human-readable status / failure reason.
};

///----------------------------------------
/// @brief Fit a hyperbola to the focus curve implied by a set of FITS frames.
/// @details For each image: loads the FITS, reads @c FOCUSPOS from the header,
///          runs @ref astap::stacking::analyse_image for the median HFD, and
///          accumulates (position, hfd) points for
///          @ref astap::core::find_best_hyperbola_fit. Frames with no FOCUSPOS
///          or no measurable stars are skipped; ≥ 3 usable frames are required
///          for a fit.
///  @param images List of FITS paths (focus sweep).
///  @param samples_out Optional contiguous buffer (size == images.size()) that
///                    receives per-frame measurements (position, hfd, stars)
///                    including rejected frames. Pass @c nullptr to skip.
/// @return Best-focus parameters + residuals, or @c ok=false with a message.
///----------------------------------------

[[nodiscard]] FocusFitResult fit_focus_hyperbola(
    std::span<const std::filesystem::path> images,
    FocusFitSample* samples_out = nullptr);

} // namespace
