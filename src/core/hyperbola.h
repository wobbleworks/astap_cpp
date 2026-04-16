///----------------------------------------
///      @file hyperbola.h
///   @ingroup ASTAP++
///     @brief Hyperbola modeling of star disk size (HFD) vs. focuser position.
///   @details Fits a vertical hyperbola to the V-curve of half-flux-diameter
///            measurements to find the optimal focus position. Used by the
///            -focus1 -focus2 CLI option.
///    @author Ported from Han Kleijn's unit_hyperbola.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <cstddef>
#include <span>

///----------------------------------------
namespace astap::core {
///----------------------------------------

/// MARK: Types

///----------------------------------------
/// @brief A single focus-curve sample: focuser position and measured HFD.
///----------------------------------------

struct FocusPoint {
    double position;
    double hfd;
};

///----------------------------------------
/// @brief Result of a hyperbola curve fit.
/// @details focus_position is the focuser position at the hyperbola minimum
///          (best focus). a and b are the hyperbola parameters where a is the
///          minimum HFD and the asymptote is y = +-x*a/b. mean_error is the
///          average relative HFD error per point after convergence.
///----------------------------------------

struct HyperbolaFit {
    double focus_position;
    double a;
    double b;
    double mean_error;
};

/// MARK: Fitting

///----------------------------------------
/// @brief Fit a hyperbola to focus-curve data and find the best focus position.
/// @details Uses a coarse sweep that is progressively refined until convergence.
/// @param data Focus-curve samples (focuser position vs. star HFD).
/// @return Best-fit hyperbola parameters including focus position.
///----------------------------------------

[[nodiscard]] HyperbolaFit find_best_hyperbola_fit(std::span<const FocusPoint> data);

///----------------------------------------
/// @brief Calculate HFD from focuser position using hyperbola parameters.
/// @param position Current focuser position.
/// @param perfect_focus Focuser position at the hyperbola minimum.
/// @param a Minimum HFD value at focus.
/// @param b Hyperbola shape parameter.
/// @return Predicted HFD at the given position.
///----------------------------------------

[[nodiscard]] double hfd_calc(double position, double perfect_focus, double a, double b) noexcept;

///----------------------------------------
/// @brief Calculate focuser steps from a given HFD to perfect focus.
/// @details Returns the positive branch (there are two solutions on either
///          side of the hyperbola).
/// @param hfd Current half-flux-diameter.
/// @param a Minimum HFD value at focus.
/// @param b Hyperbola shape parameter.
/// @return Number of focuser steps to reach perfect focus.
///----------------------------------------

[[nodiscard]] double steps_to_focus(double hfd, double a, double b) noexcept;
    
} // namespace
