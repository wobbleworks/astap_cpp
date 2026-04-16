///----------------------------------------
///      @file hyperbola.cpp
///   @ingroup ASTAP++
///     @brief Hyperbola modeling of star disk size (HFD) vs. focuser position.
///    @author Ported from Han Kleijn's unit_hyperbola.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "hyperbola.h"

#include <cmath>
#include <cstddef>
#include <span>

///----------------------------------------
namespace astap::core {
///----------------------------------------

/// MARK: - Anonymous Helpers

namespace {

///----------------------------------------
/// @brief Calculate the averaged error between measured V-curve and hyperbola.
/// @details The error per sample is outlier-limited: overshoot is scaled by the
///          measurement, undershoot by the simulation, saturating near 1.
/// @param data Focus-curve samples.
/// @param perfect_focus Focuser position at the hyperbola minimum.
/// @param a Minimum HFD value at focus.
/// @param b Hyperbola shape parameter.
/// @return Mean relative HFD error per point.
///----------------------------------------

[[nodiscard]] double mean_error_hyperbola(std::span<const FocusPoint> data,
                                          double perfect_focus,
                                          double a, double b) noexcept {
    auto total_error = 0.0;
    auto n = data.size();
    
    for (std::size_t i = 0; i < n; ++i) {
        auto hfd_simulation = hfd_calc(data[i].position, perfect_focus, a, b);
        
        // Smart error calculation which limits error for outliers
        auto error = hfd_simulation - data[i].hfd;
        if (error < 0.0) {
            total_error -= error / data[i].hfd;
        } else {
            total_error += error / hfd_simulation;
        }
    }
    
    return total_error / static_cast<double>(n);
}
    
}  // namespace

/// MARK: - HFD Calculation

double hfd_calc(double position, double perfect_focus, double a, double b) noexcept {
    // Vertical hyperbola: sqr(y/a) - sqr(x/b) = 1
    auto x = perfect_focus - position;
    return a * std::sqrt(1.0 + (x / b) * (x / b));
}

/// MARK: - Steps to Focus

double steps_to_focus(double hfd, double a, double b) noexcept {
    // Vertical hyperbola: sqr(y/a) - sqr(x/b) = 1, positive branch
    auto k = hfd / a;
    if (k < 1.0) {
        k = 1.0;
    }
    return b * std::sqrt(k * k - 1.0);
}

/// MARK: - Curve Fitting

HyperbolaFit find_best_hyperbola_fit(std::span<const FocusPoint> data) {
    auto lowest_error = 1e99;
    auto old_error = 1e99;
    auto n = data.size();
    
    // Find start values for the hyperbola loop
    auto highest_hfd = 0.0;
    auto lowest_hfd = 1e99;
    auto highest_hfd_position = 0.0;
    auto lowest_hfd_position = 0.0;
    
    for (std::size_t i = 0; i < n; ++i) {
        if (data[i].hfd > highest_hfd) {
            highest_hfd = data[i].hfd;
            highest_hfd_position = data[i].position;
        }
        if (data[i].hfd < lowest_hfd && data[i].hfd > 0.1) {
            lowest_hfd = data[i].hfd;
            lowest_hfd_position = data[i].position;
        }
    }
    
    if (highest_hfd_position < lowest_hfd_position) {
        // Go up always
        highest_hfd_position =
            (lowest_hfd_position - highest_hfd_position) + lowest_hfd_position;
    }
    
    // Get good starting values for a, b and p
    auto a = lowest_hfd;
    auto b = (highest_hfd_position - lowest_hfd_position) /
             std::sqrt(-1.0 + (highest_hfd / a) * (highest_hfd / a));
    auto p = lowest_hfd_position;
    
    auto iteration_cycles = 0;
    
    // Set starting test range
    auto a_range = a;
    auto b_range = b;
    auto p_range = (highest_hfd_position - lowest_hfd_position);
    
    do {
        auto p0 = p;
        auto b0 = b;
        auto a0 = a;
        
        // Reduce scan range by 50%
        a_range *= 0.5;
        b_range *= 0.5;
        p_range *= 0.5;
        
        // Position loop
        auto p1 = p0 - p_range;
        while (p1 <= p0 + p_range) {
            // a loop
            auto a1 = a0 - a_range;
            while (a1 <= a0 + a_range) {
                // b loop
                auto b1 = b0 - b_range;
                while (b1 <= b0 + b_range) {
                    auto error1 = mean_error_hyperbola(data, p1, a1, b1);
                    if (error1 < lowest_error) {
                        old_error = lowest_error;
                        lowest_error = error1;
                        a = a1;
                        b = b1;
                        p = p1;
                    }
                    b1 = b1 + b_range * 0.1;
                }
                a1 = a1 + a_range * 0.1;
            }
            p1 = p1 + p_range * 0.1;
        }
        ++iteration_cycles;
    } while (!((old_error - lowest_error < 1e-5)
               || (lowest_error <= 1e-5)
               || (iteration_cycles >= 30)));
               
    return HyperbolaFit{p, a, b, lowest_error};
}
    
} // namespace
