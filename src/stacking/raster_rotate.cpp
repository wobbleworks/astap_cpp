///----------------------------------------
///      @file raster_rotate.cpp
///   @ingroup ASTAP++
///     @brief Area-weighted raster rotation implementation.
///   @details Rotates an image by gathering, for each destination pixel, the
///            flux contribution of up to 9 overlapping source pixels based on
///            the fractional area covered by the rotated destination square.
///    @author Ported from Han Kleijn's unit_raster_rotate.pas (ASTAP),
///            itself based on code by Sudonull (2012):
///            https://sudonull.com/post/134233-Precise-rotation-of-the-bitmap-image-at-an-arbitrary-angle
///            MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / Sudonull / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "raster_rotate.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <vector>

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

///----------------------------------------
/// MARK: File-local types and state
///----------------------------------------

namespace {

///----------------------------------------
/// @brief Extended-precision point used for bounding-box corners.
///----------------------------------------

struct ExPoint {
    long double x;
    long double y;
};

///----------------------------------------
///   @brief Single source-pixel contribution to a destination pixel.
/// @details @c weight is the fraction of the destination square covered
///          (range 0..1); @c x and @c y are source-pixel coordinates.
///----------------------------------------

struct PixelContribution {
    float weight;
    int   y;
    int   x;
};

/// @brief Up to 9 contributing source pixels per destination pixel.
using Flux = std::array<PixelContribution, 9>;

// File-local state shared between raster_rotate and its helper. These used to
// be unit-level globals in the original source and are scoped to the anonymous
// namespace here so callers don't see them.
int s_width   = 0;
int s_height  = 0;
int s_colours = 0;

double tgAd2       = 0.0;
double tgBd2       = 0.0;
double sinB        = 0.0;
double cosB        = 0.0;
double sinBucosBd2 = 0.0;
double cosAu2      = 0.0;
double cosBu2      = 0.0;
double triS        = 0.0;
double pdx         = 0.0;
double pdy         = 0.0;

///----------------------------------------
/// MARK: Numeric helpers
///----------------------------------------

///----------------------------------------
/// @brief Fractional part of @p x, carrying the sign of @p x (matches the
///        behavior of the original @c frac intrinsic).
///----------------------------------------

[[nodiscard]] inline double frac_signed(double x) noexcept {
    auto ipart = 0.0;
    return std::modf(x, &ipart);
}

///----------------------------------------
/// @brief Truncate @p x toward zero and return the result as an @c int.
///----------------------------------------

[[nodiscard]] inline int trunc_int(double x) noexcept {
    return static_cast<int>(std::trunc(x));
}

///----------------------------------------
///   @brief Angle in degrees from the +Y axis (12 o'clock) to the line from
///          (@p CentrX, @p CentrY) to (@p x, @p y), measured clockwise.
///   @param x Point X coordinate.
///   @param y Point Y coordinate.
///   @param CentrX Center X coordinate.
///   @param CentrY Center Y coordinate.
///   @param R Precomputed distance from center to point.
///  @return Angle in the range [0, 360) degrees, or 0 when @p R is zero.
///----------------------------------------

[[nodiscard]] double Point2Ang(double x, double y, double CentrX, double CentrY, double R) noexcept {
    if (R != 0.0) {
        const auto dx = x - CentrX;
        const auto dy = y - CentrY;
        if (dx < 0.0) {
            if (dy < 0.0) {
                return 360.0 + 180.0 * std::asin(dx / R) / M_PI;
            }
            else {
                return 180.0 - 180.0 * std::asin(dx / R) / M_PI;
            }
        }
        else {
            if (dy < 0.0) {
                return 180.0 * std::asin(dx / R) / M_PI;
            }
            else {
                return 180.0 - 180.0 * std::asin(dx / R) / M_PI;
            }
        }
    }
    return 0.0;
}

///----------------------------------------
///   @brief Inverse of @ref Point2Ang.
/// @details Returns the point at distance @p R and given @p angle (0..360
///          degrees) from origin, measured clockwise from +Y.
///   @param angle Angle in degrees.
///   @param R Distance from origin.
///  @return Offset from origin as an @ref ExPoint.
///----------------------------------------

[[nodiscard]] ExPoint Ang2Point(double angle, double R) noexcept {
    auto result = ExPoint{0.0L, 0.0L};
    const auto z = trunc_int(angle / 180.0);
    switch (z) {
        case 0:
            result.x = static_cast<long double>(R *  std::sin(M_PI * angle / 180.0));
            result.y = static_cast<long double>(R * -std::cos(M_PI * angle / 180.0));
            break;
        case 1: {
            angle = angle - 180.0 * z;
            result.x = static_cast<long double>(R * -std::sin(M_PI * angle / 180.0));
            result.y = static_cast<long double>(R *  std::cos(M_PI * angle / 180.0));
            break;
        }
        default:
            // Only z in {0,1} is meaningful; leave at zero otherwise.
            break;
    }
    return result;
}

///----------------------------------------
///   @brief Round @p x away from zero unless its fractional part is below
///          1e-6 in magnitude, in which case round toward zero.
/// @details Used to bound the rotated image without losing edge pixels.
///----------------------------------------

[[nodiscard]] int DoMax(double x) noexcept {
    if (x < 0.0) {
        if (frac_signed(std::fabs(x)) > 1e-6) {
            return trunc_int(x) - 1;
        }
        else {
            return trunc_int(x);
        }
    }
    else {
        if (frac_signed(x) > 1e-6) {
            return trunc_int(x) + 1;
        }
        else {
            return trunc_int(x);
        }
    }
}

///----------------------------------------
/// MARK: Source-pixel contribution computation
///----------------------------------------

///----------------------------------------
///   @brief For destination pixel (@p x, @p y), compute up to 9 source pixel
///          contributions (coordinate + fractional area).
/// @details If the rotated center of the destination pixel falls outside the
///          source image, the first entry's weight is set to -1.0f as a
///          sentinel. Otherwise each returned entry gives an overlapping
///          source-pixel coordinate and the fractional area covered.
///   @param x Destination X coordinate.
///   @param y Destination Y coordinate.
///   @param CX Center-of-rotation X coordinate.
///   @param CY Center-of-rotation Y coordinate.
///   @param angle Rotation angle in degrees.
///  @return Up to 9 @ref PixelContribution entries; unused slots are zeroed.
///----------------------------------------

[[nodiscard]] Flux calculate_relevant_source_pixels(int x, int y, double CX, double CY, double angle) {
    auto result = Flux{};
    for (auto& c : result) {
        c.weight = 0.0f;
        c.x = 0;
        c.y = 0;
    }
    
    const auto dx1 = x + 0.5 - CX;  // top-left -> center
    const auto dy1 = y + 0.5 - CY;
    const auto R   = std::sqrt(dx1 * dx1 + dy1 * dy1);
    
    auto ug = 0.0;
    if (R > 0.0) {
        // Inlined Point2Ang(x+0.5, y+0.5, CX, CY, R) + angle.
        if (dx1 < 0.0) {
            if (dy1 < 0.0) {
                ug = angle + 360.0 + 180.0 * std::asin(dx1 / R) / M_PI;
            }
            else {
                ug = angle + 180.0 - 180.0 * std::asin(dx1 / R) / M_PI;
            }
        }
        else {
            if (dy1 < 0.0) {
                ug = angle + 180.0 * std::asin(dx1 / R) / M_PI;
            }
            else {
                ug = angle + 180.0 - 180.0 * std::asin(dx1 / R) / M_PI;
            }
        }
    }
    else {
        ug = angle;
    }
    
    if (ug >= 360.0) {
        ug = ug - 360.0 * trunc_int(ug / 360.0);
    }
    
    // Inlined Ang2Point(ug, R) + (CX, CY).
    auto vcx = 0.0;
    auto vcy = 0.0;
    const auto half = trunc_int(ug / 180.0);
    if (half == 0) {
        vcx = CX + R *  std::sin(M_PI * ug / 180.0);
        vcy = CY - R *  std::cos(M_PI * ug / 180.0);
    }
    else {
        ug = ug - 180.0 * trunc_int(ug / 180.0);
        vcx = CX - R *  std::sin(M_PI * ug / 180.0);
        vcy = CY + R *  std::cos(M_PI * ug / 180.0);
    }
    
    if (vcx < 0.0 || vcy < 0.0 || vcx >= s_width || vcy >= s_height) {
        result[0].weight = -1.0f;  // out-of-image sentinel
        return result;
    }
    
    // Coordinates of the rotated destination square's corners.
    double x1, y1, x2, y2, x3, y3, x4, y4;
    if (vcx < 3.0 || vcy < 3.0) {
        // Shift by +10 so frac() arguments stay in the positive quadrant.
        x1 = vcx + pdx + 10.0;  y1 = vcy + pdy + 10.0;
        x2 = vcx - pdy + 10.0;  y2 = vcy + pdx + 10.0;
        x3 = vcx - pdx + 10.0;  y3 = vcy - pdy + 10.0;
        x4 = vcx + pdy + 10.0;  y4 = vcy - pdx + 10.0;
    }
    else {
        x1 = vcx + pdx;  y1 = vcy + pdy;
        x2 = vcx - pdy;  y2 = vcy + pdx;
        x3 = vcx - pdx;  y3 = vcy - pdy;
        x4 = vcx + pdy;  y4 = vcy - pdx;
    }
    
    const auto ix1 = trunc_int(x1);
    const auto iy1 = trunc_int(y1);
    const auto ix2 = trunc_int(x2);
    const auto iy2 = trunc_int(y2);
    const auto ix3 = trunc_int(x3);
    const auto iy3 = trunc_int(y3);
    const auto ix4 = trunc_int(x4);
    const auto iy4 = trunc_int(y4);
    
    // Case analysis over how the rotated square overlaps the integer grid.
    // All formulas below assume the corners are in the positive quadrant of
    // the coordinate system. Translated mechanically; do not simplify.
    double dx1f, dx3, dx4, dy2, dy3, dy4;
    double s1, s2, s3, s4, s5, s6, s7, s8, s9;
    
    if (iy3 == iy2) {
        if (iy4 == iy1) {
            if (ix4 == ix3) {
                if (ix1 == ix2) {
                    // option 1
                    dy2 = frac_signed(y2);
                    dx3 = 1.0 - frac_signed(x3);
                    dy3 = frac_signed(y3);
                    dx4 = 1.0 - frac_signed(x4);
                    dy4 = 1.0 - frac_signed(y4);
                    s3 = dx3 * dy3 + (dx3 * dx3 - dy3 * dy3) * tgAd2;
                    s4 = dx4 * dy4 + (dy4 * dy4 - dx4 * dx4) * tgAd2;
                    s1 = 0.5 + (dy4 - dy2) / cosAu2;
                    s2 = 1.0 - s1 - s3;
                    s1 = s1 - s4;
                    result[0] = {(float)s1, iy1, ix1};
                    result[1] = {(float)s2, iy2, ix2};
                    result[2] = {(float)s3, iy3, ix3};
                    result[3] = {(float)s4, iy4, ix4};
                }
                else {
                    // option 15
                    dx1f = frac_signed(x1);
                    dy2  = frac_signed(y2);
                    dy4  = 1.0 - frac_signed(y4);
                    s4 = dx1f / cosB + dy2 / sinB - 1.0;
                    s4 = s4 * s4 * triS;
                    result[0] = {(float)s4, iy2, ix1};
                    s2 = dx1f * dx1f / sinBucosBd2 - s4;
                    s1 = 0.5 + (dy4 - dy2) / cosAu2;
                    s3 = 1.0 - s1 - s4;
                    s1 = s1 - s2;
                    result[1] = {(float)s1, iy4, ix4};
                    result[2] = {(float)s2, iy1, ix1};
                    result[3] = {(float)s3, iy2, ix2};
                }
            }
            else if (ix1 == ix2) {
                // option 11
                dy2 = frac_signed(y2);
                dx3 = 1.0 - frac_signed(x3);
                dy4 = 1.0 - frac_signed(y4);
                s1 = dx3 / cosB + dy4 / sinB - 1.0;
                s1 = s1 * s1 * triS;
                result[0] = {(float)s1, iy4, ix3};
                s3 = dx3 * dx3 / sinBucosBd2 - s1;
                s2 = 0.5 + (dy4 - dy2) / cosAu2;
                s4 = 1.0 - s2 - s3;
                s2 = s2 - s1;
                result[1] = {(float)s2, iy4, ix4};
                result[2] = {(float)s3, iy3, ix3};
                result[3] = {(float)s4, iy2, ix2};
            }
            else {
                // option 4
                dx1f = frac_signed(x1);
                dy2  = frac_signed(y2);
                dx3  = 1.0 - frac_signed(x3);
                dy4  = 1.0 - frac_signed(y4);
                s1 = dx3 / cosB + dy4 / sinB - 1.0;
                if (s1 > 0.0) {
                    s1 = s1 * s1 * triS;
                    result[0] = {(float)s1, iy4, ix3};
                }
                else {
                    s1 = 0.0;
                }
                s6 = dx1f / cosB + dy2 / sinB - 1.0;
                if (s6 > 0.0) {
                    s6 = s6 * s6 * triS;
                    result[1] = {(float)s6, iy2, ix1};
                }
                else {
                    s6 = 0.0;
                }
                s4 = dx3 * dx3 / sinBucosBd2 - s1;
                s3 = dx1f * dx1f / sinBucosBd2 - s6;
                s2 = 0.5 + (dy4 - dy2) / cosAu2;
                s5 = 1.0 - s2 - s4 - s6;
                s2 = s2 - s1 - s3;
                result[2] = {(float)s2, iy4, ix4};
                result[3] = {(float)s4, iy3, ix3};
                result[4] = {(float)s5, iy2, ix2};
                result[5] = {(float)s3, iy1, ix1};
            }
        }
        else if (ix4 == ix2) {
            if (ix4 == ix1) {
                // option 10
                dx3 = 1.0 - frac_signed(x3);
                dy4 = 1.0 - frac_signed(y4);
                s1 = dx3 / cosB + dy4 / sinB - 1.0;
                s1 = s1 * s1 * triS;
                result[0] = {(float)s1, iy4, ix3};
                s2 = dy4 * dy4 / sinBucosBd2 - s1;
                s3 = dx3 * dx3 / sinBucosBd2 - s1;
                s4 = 1.0 - s1 - s2 - s3;
                result[1] = {(float)s2, iy4, ix4};
                result[2] = {(float)s3, iy3, ix3};
                result[3] = {(float)s4, iy3, ix4};
            }
            else if (ix3 == ix2) {
                // option 17
                dx1f = frac_signed(x1);
                dy4  = 1.0 - frac_signed(y4);
                s2 = dx1f / sinB + dy4 / cosB - 1.0;
                s2 = s2 * s2 * triS;
                result[0] = {(float)s2, iy4, ix1};
                s1 = dy4 * dy4 / sinBucosBd2 - s2;
                s4 = dx1f * dx1f / sinBucosBd2 - s2;
                s3 = 1.0 - s1 - s2 - s4;
                result[1] = {(float)s1, iy4, ix4};
                result[2] = {(float)s4, iy1, ix1};
                result[3] = {(float)s3, iy1, ix4};
            }
            else {
                // option 5
                dx1f = frac_signed(x1);
                dx3  = 1.0 - frac_signed(x3);
                dy4  = 1.0 - frac_signed(y4);
                s1 = dx3 / cosB + dy4 / sinB - 1.0;
                if (s1 > 0.0) {
                    s1 = s1 * s1 * triS;
                    result[0] = {(float)s1, iy4, ix3};
                }
                else {
                    s1 = 0.0;
                }
                s3 = dx1f / sinB + dy4 / cosB - 1.0;
                if (s3 > 0.0) {
                    s3 = s3 * s3 * triS;
                    result[1] = {(float)s3, iy4, ix1};
                }
                else {
                    s3 = 0.0;
                }
                s2 = dy4 * dy4 / sinBucosBd2;
                s4 = dx3 * dx3 / sinBucosBd2 - s1;
                s6 = dx1f * dx1f / sinBucosBd2 - s3;
                s5 = 1.0 - s2 - s4 - s6;
                s2 = s2 - s1 - s3;
                result[2] = {(float)s2, iy4, ix4};
                result[3] = {(float)s4, iy3, ix3};
                result[4] = {(float)s6, iy1, ix1};
                result[5] = {(float)s5, iy3, ix4};
            }
        }
        else if (ix4 == ix3) {
            // option 16
            dx1f = frac_signed(x1);
            dx3  = 1.0 - frac_signed(x3);
            dy4  = 1.0 - frac_signed(y4);
            s2 = dx1f / sinB + dy4 / cosB - 1.0;
            s2 = s2 * s2 * triS;
            result[0] = {(float)s2, iy4, ix1};
            s1 = dy4 * dy4 / sinBucosBd2 - s2;
            s4 = 0.5 + (dx1f - dx3) / cosAu2;
            s3 = 1.0 - s4 - s1;
            s4 = s4 - s2;
            result[1] = {(float)s1, iy4, ix4};
            result[2] = {(float)s3, iy3, ix3};
            result[3] = {(float)s4, iy1, ix1};
        }
        else {
            // option 12
            dx1f = frac_signed(x1);
            dx3  = 1.0 - frac_signed(x3);
            dy4  = 1.0 - frac_signed(y4);
            s1 = dx3 / cosB + dy4 / sinB - 1.0;
            s1 = s1 * s1 * triS;
            result[0] = {(float)s1, iy4, ix3};
            s2 = dy4 * dy4 / sinBucosBd2 - s1;
            s3 = 0.5 + (dx3 - dx1f) / cosBu2;
            s4 = 1.0 - s3 - s2;
            s3 = s3 - s1;
            result[1] = {(float)s2, iy4, ix4};
            result[2] = {(float)s3, iy3, ix3};
            result[3] = {(float)s4, iy1, ix1};
        }
    }
    else if (iy3 == iy4) {
        if (ix3 == ix2) {
            if (ix4 == ix1) {
                if (iy2 == iy1) {
                    // option 2
                    dx4  = frac_signed(x4);
                    dy2  = frac_signed(y2);
                    dy4  = 1.0 - frac_signed(y4);
                    dx1f = frac_signed(x1);
                    const auto dy1f = frac_signed(y1);
                    s4 = dx4 * dy4 + (dy4 * dy4 - dx4 * dx4) * tgBd2;
                    s1 = dx1f * dy1f + (dx1f * dx1f - dy1f * dy1f) * tgBd2;
                    s3 = 0.5 + (dy4 - dy2) / cosBu2;
                    s2 = 1.0 - s3 - s1;
                    s3 = s3 - s4;
                    result[0] = {(float)s1, iy1, ix1};
                    result[1] = {(float)s2, iy2, ix2};
                    result[2] = {(float)s3, iy3, ix3};
                    result[3] = {(float)s4, iy4, ix4};
                }
                else {
                    // option 13
                    dx1f = frac_signed(x1);
                    dy2  = frac_signed(y2);
                    dx3  = 1.0 - frac_signed(x3);
                    s4 = dx1f / cosB + dy2 / sinB - 1.0;
                    s4 = s4 * s4 * triS;
                    result[0] = {(float)s4, iy2, ix1};
                    s3 = dy2 * dy2 / sinBucosBd2 - s4;
                    s1 = 0.5 + (dx3 - dx1f) / cosBu2;
                    s2 = 1.0 - s1 - s4;
                    s1 = s1 - s3;
                    result[1] = {(float)s1, iy3, ix3};
                    result[2] = {(float)s3, iy2, ix2};
                    result[3] = {(float)s2, iy1, ix1};
                }
            }
            else if (iy2 == iy1) {
                // option 18
                dx1f = frac_signed(x1);
                dy2  = frac_signed(y2);
                dy4  = 1.0 - frac_signed(y4);
                s2 = dx1f / sinB + dy4 / cosB - 1.0;
                s2 = s2 * s2 * triS;
                result[0] = {(float)s2, iy4, ix1};
                s4 = dx1f * dx1f / sinBucosBd2 - s2;
                s1 = 0.5 + (dy4 - dy2) / cosBu2;
                s3 = 1.0 - s1 - s4;
                s1 = s1 - s2;
                result[1] = {(float)s1, iy4, ix4};
                result[2] = {(float)s3, iy2, ix2};
                result[3] = {(float)s4, iy1, ix1};
            }
            else {
                // option 14
                dx1f = frac_signed(x1);
                dy2  = frac_signed(y2);
                s4 = dx1f / cosB + dy2 / sinB - 1.0;
                s4 = s4 * s4 * triS;
                result[0] = {(float)s4, iy2, ix1};
                s2 = dx1f * dx1f / sinBucosBd2 - s4;
                s3 = dy2 * dy2 / sinBucosBd2 - s4;
                s1 = 1.0 - s2 - s3 - s4;
                result[1] = {(float)s2, iy1, ix1};
                result[2] = {(float)s3, iy2, ix2};
                result[3] = {(float)s1, iy1, ix2};
            }
        }
        else if (ix2 == ix1) {
            if (ix3 == ix4) {
                // option 9
                dx1f = frac_signed(x1);
                dy2  = frac_signed(y2);
                dx3  = 1.0 - frac_signed(x3);
                s3 = dx3 / sinB + dy2 / cosB - 1.0;
                s3 = s3 * s3 * triS;
                result[0] = {(float)s3, iy2, ix3};
                s4 = dy2 * dy2 / sinBucosBd2 - s3;
                s2 = 0.5 + (dx1f - dx3) / cosAu2;
                s1 = 1.0 - s2 - s3;
                s2 = s2 - s4;
                result[1] = {(float)s4, iy2, ix2};
                result[2] = {(float)s2, iy1, ix1};
                result[3] = {(float)s1, iy3, ix3};
            }
            else if (iy3 == iy1) {
                // option 8
                dy2 = frac_signed(y2);
                dx3 = 1.0 - frac_signed(x3);
                s3 = dx3 / sinB + dy2 / cosB - 1.0;
                s3 = s3 * s3 * triS;
                result[0] = {(float)s3, iy2, ix3};
                s1 = dx3 * dx3 / sinBucosBd2 - s3;
                s4 = dy2 * dy2 / sinBucosBd2 - s3;
                s2 = 1.0 - s1 - s3 - s4;
                result[1] = {(float)s1, iy3, ix3};
                result[2] = {(float)s4, iy2, ix2};
                result[3] = {(float)s2, iy3, ix2};
            }
            else {
                // option 7
                dy2 = frac_signed(y2);
                dx3 = 1.0 - frac_signed(x3);
                dy4 = 1.0 - frac_signed(y4);
                s3 = dx3 / sinB + dy2 / cosB - 1.0;
                s3 = s3 * s3 * triS;
                result[0] = {(float)s3, iy2, ix3};
                s1 = dx3 * dx3 / sinBucosBd2 - s3;
                s2 = 0.5 + (dy4 - dy2) / cosBu2;
                s4 = 1.0 - s2 - s3;
                s2 = s2 - s1;
                result[1] = {(float)s1, iy3, ix3};
                result[2] = {(float)s2, iy4, ix4};
                result[3] = {(float)s4, iy2, ix4};
            }
        }
        else if (iy3 == iy1) {
            // option 6
            dx1f = frac_signed(x1);
            dy2  = frac_signed(y2);
            dx3  = 1.0 - frac_signed(x3);
            s4 = dx3 / sinB + dy2 / cosB - 1.0;
            if (s4 > 0.0) {
                s4 = s4 * s4 * triS;
                result[0] = {(float)s4, iy2, ix3};
            }
            else {
                s4 = 0.0;
            }
            s6 = dx1f / cosB + dy2 / sinB - 1.0;
            if (s6 > 0.0) {
                s6 = s6 * s6 * triS;
                result[1] = {(float)s6, iy2, ix1};
            }
            else {
                s6 = 0.0;
            }
            s1 = dx3 * dx3 / sinBucosBd2 - s4;
            s3 = dx1f * dx1f / sinBucosBd2 - s6;
            s5 = dy2 * dy2 / sinBucosBd2;
            s2 = 1.0 - s1 - s3 - s5;
            s5 = s5 - s4 - s6;
            result[2] = {(float)s1, iy3, ix3};
            result[3] = {(float)s3, iy1, ix1};
            result[4] = {(float)s5, iy2, ix2};
            result[5] = {(float)s2, iy3, ix2};
        }
        else {
            // option 3
            dx1f = frac_signed(x1);
            dy2  = frac_signed(y2);
            dx3  = 1.0 - frac_signed(x3);
            dy4  = 1.0 - frac_signed(y4);
            s4 = dx3 / sinB + dy2 / cosB - 1.0;
            if (s4 > 0.0) {
                s4 = s4 * s4 * triS;
                result[0] = {(float)s4, iy2, ix3};
            }
            else {
                s4 = 0.0;
            }
            s3 = dx1f / sinB + dy4 / cosB - 1.0;
            if (s3 > 0.0) {
                s3 = s3 * s3 * triS;
                result[1] = {(float)s3, iy4, ix1};
            }
            else {
                s3 = 0.0;
            }
            s1 = dx3 * dx3 / sinBucosBd2 - s4;
            s6 = dx1f * dx1f / sinBucosBd2 - s3;
            s2 = 0.5 + (dy4 - dy2) / cosBu2;
            s5 = 1.0 - s2 - s4 - s6;
            s2 = s2 - s1 - s3;
            result[2] = {(float)s6, iy1, ix1};
            result[3] = {(float)s5, iy2, ix2};
            result[4] = {(float)s1, iy3, ix3};
            result[5] = {(float)s2, iy4, ix4};
        }
    }
    else if (ix4 == ix2) {
        if (ix3 == ix4) {
            // option 21
            dx1f = frac_signed(x1);
            dy2  = frac_signed(y2);
            dy4  = 1.0 - frac_signed(y4);
            s2 = dx1f / sinB + dy4 / cosB - 1.0;
            if (s2 > 0.0) {
                s2 = s2 * s2 * triS;
                result[0] = {(float)s2, iy4, ix1};
            }
            else {
                s2 = 0.0;
            }
            s6 = dx1f / cosB + dy2 / sinB - 1.0;
            if (s6 > 0.0) {
                s6 = s6 * s6 * triS;
                result[1] = {(float)s6, iy2, ix1};
            }
            else {
                s6 = 0.0;
            }
            s1 = dy4 * dy4 / sinBucosBd2 - s2;
            s5 = dy2 * dy2 / sinBucosBd2 - s6;
            s4 = dx1f * dx1f / sinBucosBd2;
            s3 = 1.0 - s4 - s1 - s5;
            s4 = s4 - s2 - s6;
            result[2] = {(float)s1, iy4, ix4};
            result[3] = {(float)s4, iy1, ix1};
            result[4] = {(float)s5, iy2, ix2};
            result[5] = {(float)s3, iy1, ix2};
        }
        else if (ix4 == ix1) {
            // option 2 (variant)
            dy2 = frac_signed(y2);
            dx3 = 1.0 - frac_signed(x3);
            dy4 = 1.0 - frac_signed(y4);
            s1 = dx3 / cosB + dy4 / sinB - 1.0;
            if (s1 > 0.0) {
                s1 = s1 * s1 * triS;
                result[0] = {(float)s1, iy4, ix3};
            }
            else {
                s1 = 0.0;
            }
            s5 = dx3 / sinB + dy2 / cosB - 1.0;
            if (s5 > 0.0) {
                s5 = s5 * s5 * triS;
                result[1] = {(float)s5, iy2, ix3};
            }
            else {
                s5 = 0.0;
            }
            s2 = dy4 * dy4 / sinBucosBd2 - s1;
            s6 = dy2 * dy2 / sinBucosBd2 - s5;
            s3 = dx3 * dx3 / sinBucosBd2;
            s4 = 1.0 - s3 - s2 - s6;
            s3 = s3 - s1 - s5;
            result[2] = {(float)s2, iy4, ix4};
            result[3] = {(float)s3, iy3, ix3};
            result[4] = {(float)s6, iy2, ix2};
            result[5] = {(float)s4, iy3, ix2};
        }
        else {
            // option 23
            dx1f = frac_signed(x1);
            dy2  = frac_signed(y2);
            dx3  = 1.0 - frac_signed(x3);
            dy4  = 1.0 - frac_signed(y4);
            s1 = dx3 / cosB + dy4 / sinB - 1.0;
            if (s1 > 0.0) {
                s1 = s1 * s1 * triS;
                result[0] = {(float)s1, iy4, ix3};
            }
            else {
                s1 = 0.0;
            }
            s3 = dx1f / sinB + dy4 / cosB - 1.0;
            if (s3 > 0.0) {
                s3 = s3 * s3 * triS;
                result[1] = {(float)s3, iy4, ix1};
            }
            else {
                s3 = 0.0;
            }
            s7 = dx3 / sinB + dy2 / cosB - 1.0;
            if (s7 > 0.0) {
                s7 = s7 * s7 * triS;
                result[2] = {(float)s7, iy2, ix3};
            }
            else {
                s7 = 0.0;
            }
            s9 = dx1f / cosB + dy2 / sinB - 1.0;
            if (s9 > 0.0) {
                s9 = s9 * s9 * triS;
                result[3] = {(float)s9, iy2, ix1};
            }
            else {
                s9 = 0.0;
            }
            s2 = dy4 * dy4 / sinBucosBd2 - s1 - s3;
            s8 = dy2 * dy2 / sinBucosBd2 - s7 - s9;
            s4 = dx3 * dx3 / sinBucosBd2;
            s6 = dx1f * dx1f / sinBucosBd2;
            s5 = 1.0 - s4 - s6 - s2 - s8;
            s4 = s4 - s1 - s7;
            s6 = s6 - s3 - s9;
            result[4] = {(float)s4, iy3, ix3};
            result[5] = {(float)s2, iy4, ix4};
            result[6] = {(float)s6, iy1, ix1};
            result[7] = {(float)s8, iy2, ix2};
            result[8] = {(float)s5, iy3, ix4};
        }
    }
    else if (ix3 == ix4) {
        // option 20
        dx1f = frac_signed(x1);
        dy2  = frac_signed(y2);
        dx3  = 1.0 - frac_signed(x3);
        dy4  = 1.0 - frac_signed(y4);
        s2 = dx1f / sinB + dy4 / cosB - 1.0;
        if (s2 > 0.0) {
            s2 = s2 * s2 * triS;
            result[0] = {(float)s2, iy4, ix1};
        }
        else {
            s2 = 0.0;
        }
        s5 = dx3 / sinB + dy2 / cosB - 1.0;
        if (s5 > 0.0) {
            s5 = s5 * s5 * triS;
            result[1] = {(float)s5, iy2, ix3};
        }
        else {
            s5 = 0.0;
        }
        s6 = dy2 * dy2 / sinBucosBd2 - s5;
        s1 = dy4 * dy4 / sinBucosBd2 - s2;
        s3 = 0.5 + (dx3 - dx1f) / cosAu2;
        s4 = 1.0 - s3 - s2 - s6;
        s3 = s3 - s1 - s5;
        result[2] = {(float)s1, iy4, ix4};
        result[3] = {(float)s3, iy3, ix3};
        result[4] = {(float)s4, iy1, ix1};
        result[5] = {(float)s6, iy2, ix2};
    }
    else {
        // option 19
        dx1f = frac_signed(x1);
        dy2  = frac_signed(y2);
        dx3  = 1.0 - frac_signed(x3);
        dy4  = 1.0 - frac_signed(y4);
        s1 = dx3 / cosB + dy4 / sinB - 1.0;
        if (s1 > 0.0) {
            s1 = s1 * s1 * triS;
            result[0] = {(float)s1, iy4, ix3};
        }
        else {
            s1 = 0.0;
        }
        s6 = dx1f / cosB + dy2 / sinB - 1.0;
        if (s6 > 0.0) {
            s6 = s6 * s6 * triS;
            result[1] = {(float)s6, iy2, ix1};
        }
        else {
            s6 = 0.0;
        }
        s2 = dy4 * dy4 / sinBucosBd2 - s1;
        s5 = dy2 * dy2 / sinBucosBd2 - s6;
        s3 = 0.5 + (dx3 - dx1f) / cosBu2;
        s4 = 1.0 - s3 - s2 - s6;
        s3 = s3 - s1 - s5;
        result[2] = {(float)s2, iy4, ix4};
        result[3] = {(float)s3, iy3, ix3};
        result[4] = {(float)s4, iy1, ix1};
        result[5] = {(float)s5, iy2, ix2};
    }
    
    if (vcx < 3.0 || vcy < 3.0) {
        // Reverse the +10 shift applied above so the coordinates land back
        // on the real source pixels.
        for (auto& c : result) {
            c.x -= 10;
            c.y -= 10;
        }
    }
    
    return result;
}
    
} // namespace

///----------------------------------------
/// MARK: Public API
///----------------------------------------

void raster_rotate(double angle, double CX, double CY, ImageArray& img) {
    // Reduce the incoming angle modulo 360.
    if (angle >= 360.0) {
        angle = angle - 360.0 * trunc_int(angle / 360.0);
    }
    
    if (angle == 0.0) {
        return;
    }
    
    // Capture source dimensions.
    s_colours = static_cast<int>(img.size());
    s_height  = s_colours > 0 ? static_cast<int>(img[0].size()) : 0;
    s_width   = (s_colours > 0 && s_height > 0)
                    ? static_cast<int>(img[0][0].size()) : 0;
                    
    // Normalize angle to [0, 360).
    angle = angle - 360.0 * trunc_int(angle / 360.0);
    if (angle < 0.0) {
        angle = 360.0 - std::fabs(angle);
    }
    
    const auto U90 = trunc_int(angle / 90.0);
    if ((angle - U90 * 90.0) == 0.0) {
        // Exact multiple of 90 deg: pure pixel reshuffle, lossless and fast.
        auto temp_img = ImageArray{};
        switch (U90) {
            case 1: { // 90 deg
                temp_img.assign(s_colours,
                                std::vector<std::vector<float>>(
                                    s_width, std::vector<float>(s_height, 0.0f)));
                for (int k = 0; k < s_colours; ++k) {
                    for (int j = 0; j < s_height; ++j) {
                        for (int i = 0; i < s_width; ++i) {
                            temp_img[k][i][s_height - 1 - j] = img[k][j][i];
                        }
                    }
                }
                break;
            }
            case 2: { // 180 deg
                temp_img.assign(s_colours,
                                std::vector<std::vector<float>>(
                                    s_height, std::vector<float>(s_width, 0.0f)));
                for (int k = 0; k < s_colours; ++k) {
                    for (int j = 0; j < s_height; ++j) {
                        for (int i = 0; i < s_width; ++i) {
                            temp_img[k][s_height - 1 - j][s_width - 1 - i] = img[k][j][i];
                        }
                    }
                }
                break;
            }
            case 3: { // 270 deg
                temp_img.assign(s_colours,
                                std::vector<std::vector<float>>(
                                    s_width, std::vector<float>(s_height, 0.0f)));
                for (int k = 0; k < s_colours; ++k) {
                    for (int j = 0; j < s_height; ++j) {
                        for (int i = 0; i < s_width; ++i) {
                            temp_img[k][s_width - 1 - i][j] = img[k][j][i];
                        }
                    }
                }
                break;
            }
            default:
                // U90 == 0 already excluded by the angle == 0 check above.
                break;
        }
        
        img.clear();
        img = std::move(temp_img);
        return;
    }
    
    // Compute the bounding box of the rotated source rectangle (each corner
    // is rotated individually so we don't lose any edge pixels via DoMax).
    double R, ug, dx, dy;
    ExPoint ep1, ep2, ep3, ep4;
    
    R  = std::sqrt(CX * CX + CY * CY);
    ug = Point2Ang(0.0, 0.0, CX, CY, R) + angle;
    if (ug >= 360.0) {
        ug = ug - 360.0 * trunc_int(ug / 360.0);
    }
    ep1 = Ang2Point(ug, R);
    ep1.x += static_cast<long double>(CX);
    ep1.y += static_cast<long double>(CY);
    
    dx = CX - s_width;
    dy = CY - s_height;
    R  = std::sqrt(dx * dx + CY * CY);
    ug = Point2Ang(s_width, 0.0, CX, CY, R) + angle;
    if (ug >= 360.0) {
        ug = ug - 360.0 * trunc_int(ug / 360.0);
    }
    ep2 = Ang2Point(ug, R);
    ep2.x += static_cast<long double>(CX);
    ep2.y += static_cast<long double>(CY);
    
    R  = std::sqrt(dx * dx + dy * dy);
    ug = Point2Ang(s_width, s_height, CX, CY, R) + angle;
    if (ug >= 360.0) {
        ug = ug - 360.0 * trunc_int(ug / 360.0);
    }
    ep3 = Ang2Point(ug, R);
    ep3.x += static_cast<long double>(CX);
    ep3.y += static_cast<long double>(CY);
    
    R  = std::sqrt(CX * CX + dy * dy);
    ug = Point2Ang(0.0, s_height, CX, CY, R) + angle;
    if (ug >= 360.0) {
        ug = ug - 360.0 * trunc_int(ug / 360.0);
    }
    ep4 = Ang2Point(ug, R);
    ep4.x += static_cast<long double>(CX);
    ep4.y += static_cast<long double>(CY);
    
    // Pick the destination-corner and size mapping based on the 90-deg quadrant.
    int yd, xd, rw, rh;
    switch (U90) {
        case 0:
            yd = DoMax(static_cast<double>(ep1.y));
            xd = DoMax(static_cast<double>(ep4.x));
            rw = std::abs(DoMax(static_cast<double>(ep2.x)) - xd);
            rh = std::abs(DoMax(static_cast<double>(ep3.y)) - yd);
            break;
        case 1:
            yd = DoMax(static_cast<double>(ep4.y));
            xd = DoMax(static_cast<double>(ep3.x));
            rw = std::abs(DoMax(static_cast<double>(ep1.x)) - xd);
            rh = std::abs(DoMax(static_cast<double>(ep2.y)) - yd);
            break;
        case 2:
            yd = DoMax(static_cast<double>(ep3.y));
            xd = DoMax(static_cast<double>(ep2.x));
            rw = std::abs(DoMax(static_cast<double>(ep4.x)) - xd);
            rh = std::abs(DoMax(static_cast<double>(ep1.y)) - yd);
            break;
        default:
            yd = DoMax(static_cast<double>(ep2.y));
            xd = DoMax(static_cast<double>(ep1.x));
            rw = std::abs(DoMax(static_cast<double>(ep3.x)) - xd);
            rh = std::abs(DoMax(static_cast<double>(ep4.y)) - yd);
            break;
    }
    
    // Allocate destination buffer (typically larger than source).
    auto temp_img = ImageArray{};
    temp_img.assign(s_colours,
                    std::vector<std::vector<float>>(
                        rh, std::vector<float>(rw, 0.0f)));
                        
    // Rotate backwards so we sample from the source for each destination pixel.
    angle = 360.0 - angle;
    ug = angle - trunc_int(angle / 90.0) * 90.0;  // 0..90 deg
    
    // Precomputed helpers for calculate_relevant_source_pixels.
    pdx          = std::sqrt(0.5)  * std::sin(M_PI * (ug + 45.0) / 180.0);
    pdy          = -std::sqrt(0.5) * std::cos(M_PI * (ug + 45.0) / 180.0);
    tgBd2        = std::tan(M_PI * (90.0 - ug) / 180.0) / 2.0;
    tgAd2        = std::tan(M_PI * ug / 180.0) / 2.0;
    sinB         = std::sin(M_PI * (90.0 - ug) / 180.0);
    cosB         = std::cos(M_PI * (90.0 - ug) / 180.0);
    cosAu2       = std::cos(M_PI * ug / 180.0) * 2.0;
    cosBu2       = std::cos(M_PI * (90.0 - ug) / 180.0) * 2.0;
    triS         = std::cos(M_PI * ug / 180.0) * std::sin(M_PI * ug / 180.0) / 2.0;
    sinBucosBd2  = std::sin(M_PI * (90.0 - ug) / 180.0)
                   * std::cos(M_PI * (90.0 - ug) / 180.0) * 2.0;
                   
    // Walk every destination pixel and gather flux from the source.
    for (int i = yd; i < yd + rh; ++i) {
        for (int j = xd; j < xd + rw; ++j) {
            const auto X = calculate_relevant_source_pixels(j, i, CX, CY, angle);
            if (X[0].weight == -1.0f) {
                continue;
            }
            
            auto red   = 0.0;
            auto green = 0.0;
            auto blue  = 0.0;
            for (int k = 0; k < 9; ++k) {
                if (X[k].weight > 1e-9f) {
                    if (X[k].y >= 0 && X[k].x >= 0
                        && X[k].y < s_height && X[k].x < s_width) {
                        red += static_cast<double>(X[k].weight) * img[0][X[k].y][X[k].x];
                        if (s_colours > 2) {
                            green += static_cast<double>(X[k].weight) * img[1][X[k].y][X[k].x];
                            blue  += static_cast<double>(X[k].weight) * img[2][X[k].y][X[k].x];
                        }
                    }
                }
            }
            temp_img[0][i - yd][j - xd] = static_cast<float>(red);
            if (s_colours > 2) {
                temp_img[1][i - yd][j - xd] = static_cast<float>(green);
                temp_img[2][i - yd][j - xd] = static_cast<float>(blue);
            }
        }
    }
    
    img.clear();
    img = std::move(temp_img);
}
    
} // namespace
