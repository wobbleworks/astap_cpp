///----------------------------------------
///      @file gaussian_blur.h
///   @ingroup ASTAP++
///     @brief In-place separable Gaussian blur for image buffers.
///   @details Convolves an @ref astap::ImageArray with a 1-D Gaussian kernel in
///            two passes (horizontal, then vertical). The kernel's effective
///            radius is trimmed based on the expected data range so very small
///            weights are skipped.
///    @author Ported from Han Kleijn's unit_stack.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include "../types.h"

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

using astap::ImageArray;

///----------------------------------------
///   @brief Apply a separable Gaussian blur to @p img in place.
/// @details Uses the Gaussian as a separable kernel: a horizontal convolution
///          followed by a vertical one. Pixels outside the image edge are
///          clamped to the nearest valid column / row. Does nothing when
///          @p radius is below @c 0.001 or when @p img is empty.
/// @param[in,out] img Image buffer indexed @c [channel][row][col].
///          Supports 1-, 2-, or 3-channel buffers.
///   @param radius Gaussian standard deviation in pixels.
///----------------------------------------

void gaussian_blur2(ImageArray& img, double radius);
	
} // namespace
