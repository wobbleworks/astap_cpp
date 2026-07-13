///----------------------------------------
///      @file live_colour_correction.h
///   @ingroup ASTAP++
///     @brief Live-stacking per-channel colour correction (white-balance to red).
///   @details The pure white-balance math is split into this dependency-free unit
///            so the formula is unit-tested with hand-set backgrounds, without
///            pulling in the live-stacking / plate-solve machinery.
///    @author Ported from Han Kleijn's unit_live_stacking.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include "../types.h"   // astap::Background

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

///----------------------------------------
///   @brief Per-channel live colour-correction factors.
/// @details White-balances the green and blue channels to red (the reference):
///          each channel's source pixel becomes @c (value + add)·multiply / largest,
///          where @c add shifts its background and @c multiply matches its star level
///          to red, and @c largest renormalises so the peak stays near 65535. Red is
///          the reference (@c add=0, @c multiply=1). Ports @c auto_background_level +
///          the colour-correction path of @c unit_live_stacking.pas.
///----------------------------------------

struct ColourCorrection {
    bool   active = false;                  ///< True once derived from a 3-colour frame.
    double add[3] = {0.0, 0.0, 0.0};        ///< Per-channel additive offset (R/G/B).
    double multiply[3] = {1.0, 1.0, 1.0};   ///< Per-channel multiplier (R/G/B).
    double largest = 1.0;                   ///< Renormalisation divisor (>0).
};

///----------------------------------------
///   @brief Colour-correction math from three measured channel backgrounds.
/// @details Green/blue are balanced to @p red. A channel whose @c star_level is the
///          @c 1 "no signal" sentinel is left at identity. Ports
///          @c unit_stack.pas:11707-11730 + @c unit_live_stacking.pas:296-302.
///   @param red   Reference channel background + star level.
///   @param green Green channel background + star level.
///   @param blue  Blue channel background + star level.
///   @return Correction factors (always @c active).
///----------------------------------------

[[nodiscard]] ColourCorrection colour_correction_from_backgrounds(
    const astap::Background& red, const astap::Background& green,
    const astap::Background& blue);

} // namespace astap::stacking
