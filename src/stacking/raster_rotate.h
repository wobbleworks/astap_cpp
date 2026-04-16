///----------------------------------------
///      @file raster_rotate.h
///   @ingroup ASTAP++
///     @brief Area-weighted raster rotation of a multi-plane image about a chosen center.
///   @details Rotates an @ref ImageArray in place. Exact multiples of 90 degrees
///            take a lossless fast path (pure pixel reshuffle); other angles use
///            a per-destination-pixel flux-weighted average over the up-to-9
///            source pixels overlapping the rotated destination square. The
///            destination buffer may be larger than the source to avoid clipping.
///    @author Ported from Han Kleijn's unit_raster_rotate.pas (ASTAP),
///            itself based on code by Sudonull (2012):
///            https://sudonull.com/post/134233-Precise-rotation-of-the-bitmap-image-at-an-arbitrary-angle
///            MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / Sudonull / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include "../types.h"

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

using astap::ImageArray;

///----------------------------------------
///   @brief Rotate @p img about (@p CX, @p CY) by @p angle degrees.
/// @details Exact multiples of 90 degrees take a fast lossless path (pure pixel
///          reshuffle); other angles use a per-destination-pixel flux-weighted
///          average over the up-to-9 source pixels overlapping the rotated
///          destination square. The destination buffer may be larger than the
///          source to avoid clipping.
///   @param angle Rotation angle in degrees.
///   @param CX X coordinate of the center of rotation (source-pixel units).
///   @param CY Y coordinate of the center of rotation (source-pixel units).
/// @param[in,out] img Multi-plane image; replaced in place with the rotated result.
///----------------------------------------

void raster_rotate(double angle, double CX, double CY, ImageArray& img);
	
} // namespace
