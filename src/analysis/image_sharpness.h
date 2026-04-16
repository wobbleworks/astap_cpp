#pragma once

/// @file image_sharpness.h
/// Image sharpness measurement — ported from unit_image_sharpness.pas.
///
/// Measures image sharpness for astronomical images (Moon, Sun, stars)
/// via RMS of min/max differences in 4x4 pixel blocks (2x2 groups of
/// RGGB Bayer quads).  The result is inverted and scaled so that the
/// final value approximates a star HFD measurement: lower = sharper.

#include "../types.h"

namespace astap::analysis {

using astap::ImageArray;

/// Measure image sharpness.
///
/// Computes the root-mean-square of local contrast across 4x4 pixel
/// blocks, then inverts and scales the value so that it tracks the
/// half-flux diameter (HFD) convention: a *lower* return value means
/// a *sharper* image.
///
/// @param img  Single- or multi-channel image buffer (only channel 0
///             is used).
/// @return Approximate HFD-equivalent sharpness value.
double image_sharpness(const ImageArray& img);

}  // namespace astap::analysis
