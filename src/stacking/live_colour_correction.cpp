///----------------------------------------
///      @file live_colour_correction.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the live colour-correction white-balance math.
///    @author Ported from Han Kleijn's unit_live_stacking.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "live_colour_correction.h"

#include <algorithm>

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

ColourCorrection colour_correction_from_backgrounds(
        const astap::Background& red, const astap::Background& green,
        const astap::Background& blue) {
    ColourCorrection cc;
    // White-balance green and blue to red (unit_stack.pas:11707-11730). star_level==1
    // is the "no signal" sentinel → leave that channel at identity (add 0, multiply 1).
    if (green.star_level != 1.0 && red.star_level != 1.0) {
        cc.add[1]      = red.backgr * (green.star_level / red.star_level) - green.backgr;
        cc.multiply[1] = red.star_level / green.star_level;
    }
    if (blue.star_level != 1.0 && red.star_level != 1.0) {
        cc.add[2]      = red.backgr * (blue.star_level / red.star_level) - blue.backgr;
        cc.multiply[2] = red.star_level / blue.star_level;
    }

    // Renormalise so the brightest channel's ceiling stays at 65535
    // (unit_live_stacking.pas:296-302).
    const double sR = (65535.0 + cc.add[0]) * cc.multiply[0] / 65535.0;
    const double sG = (65535.0 + cc.add[1]) * cc.multiply[1] / 65535.0;
    const double sB = (65535.0 + cc.add[2]) * cc.multiply[2] / 65535.0;
    // Red is never rewritten, so sR ≡ 1 and largest ≥ 1 always: this floor is
    // unreachable defensive code (the Pascal max has no such guard). Kept so a future
    // change that lets red's scale move cannot divide by a non-positive value.
    cc.largest = std::max({sR, sG, sB});
    if (cc.largest <= 0.0) {
        cc.largest = 1.0;
    }
    cc.active = true;
    return cc;
}

} // namespace astap::stacking
