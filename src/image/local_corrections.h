///----------------------------------------
///      @file local_corrections.h
///   @ingroup ASTAP++
///     @brief Local pixel-level corrections: linear gradient removal and
///            dust-spot interpolation.
///   @details Ports two of @c astap_main.pas's small in-place image-fix
///            routines:
///             - @ref remove_linear_gradient subtracts a linear ramp running
///               from a "dark" image-space anchor to a "bright" anchor.
///             - @ref remove_dust_spot fills a user-selected box (optionally
///               clipped to its inscribed ellipse) with a bilinearly-interpolated
///               surface fitted to the four surrounding 20-pixel-margin patches.
///
///            Both operate on the canonical @c astap::ImageArray and read the
///            colour count from @c head.naxis3.
///    @author Ported from Han Kleijn's ASTAP (MPL-2.0) by John Stephen.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include "../types.h"

///----------------------------------------
namespace astap::image {
///----------------------------------------

///----------------------------------------
/// @brief Subtract a linear gradient between two anchor points from the image.
/// @details Around each anchor the routine samples a 41×41 patch and takes
///          its mode (most common pixel value) per channel — that defines the
///          "dark" and "bright" levels. For every pixel @c P, the orthogonal
///          projection @c p of @c P onto the line dark→bright is computed; the
///          subtraction is @c (bright_mode − dark_mode) · (a − p) / a where
///          @c a is the dark↔bright distance. Identical to @c gradient_removal1
///          in the Pascal source.
/// @param img      Image to modify in place.
/// @param head     Header (only @c width / @c height / @c naxis3 are read).
/// @param dark_x   Anchor X for the dark area (0-based pixel).
/// @param dark_y   Anchor Y for the dark area.
/// @param bright_x Anchor X for the bright area.
/// @param bright_y Anchor Y for the bright area.
/// @return @c true on success; @c false if the anchors are too close (< 100
///         pixels apart) — matches the Pascal guard.
///----------------------------------------

[[nodiscard]] bool remove_linear_gradient(ImageArray& img, const Header& head,
                                          int dark_x, int dark_y,
                                          int bright_x, int bright_y);

///----------------------------------------
/// @brief Fill a small image region with a bilinear estimate from its
///        surroundings, optionally clipping to an inscribed ellipse.
/// @details For each colour channel the routine samples the modes of the
///          four 20-pixel-thick patches just outside the box's corners and
///          builds a bilinearly-interpolated "expected" value for every pixel
///          inside the box. The difference between the expected value and the
///          actual pixel is Gaussian-blurred (radius 5) so the seam blends in,
///          then added back. When @p use_ellipse is @c true only pixels inside
///          the inscribed ellipse are touched; otherwise the entire box. Matches
///          @c dust_spot_removal1 in the Pascal source.
/// @param img         Image to modify in place.
/// @param head        Header (only @c width / @c height / @c naxis3 are read).
/// @param x1          Box bound 1 (X), 0-based pixel; order does not matter.
/// @param y1          Box bound 1 (Y).
/// @param x2          Box bound 2 (X).
/// @param y2          Box bound 2 (Y).
/// @param use_ellipse Restrict the fix to the inscribed ellipse.
/// @return @c true on success; @c false if the box is too small (Δx<3 OR Δy<3).
///----------------------------------------

[[nodiscard]] bool remove_dust_spot(ImageArray& img, const Header& head,
                                    int x1, int y1, int x2, int y2,
                                    bool use_ellipse);

} // namespace astap::image
