///----------------------------------------
///      @file demosaic.h
///   @ingroup ASTAP++
///     @brief Bayer / X-Trans demosaic family for raw sensor data.
///   @details Provides multiple demosaic algorithms (bilinear, astro-simple,
///            mono-preserving, colour-preserving, superpixel, and Fuji X-Trans)
///            plus a pattern resolver and dispatcher. Each variant allocates an
///            internal 3-channel buffer, fills it, then moves it into the
///            caller's ImageArray. The Header is updated to naxis=3 / naxis3=3
///            on success, except for preserve_colour_saturated_bayer which
///            operates on a single mosaic plane.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP), MPL-2.0.
///            C++ port by John Stephen.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <string>

#include "../types.h"

///----------------------------------------
namespace astap::core {
///----------------------------------------

/// MARK: - Bayer Pattern

///----------------------------------------
/// @brief Logical Bayer pattern enumeration.
/// @details Integer values match the convention used throughout the original ASTAP source.
///          Pattern 0-3 are standard Bayer layouts; 4 is Fuji X-Trans (6x6).
///----------------------------------------

enum class BayerPattern : int {
    GRBG   = 0,
    BGGR   = 1,
    RGGB   = 2,
    GBRG   = 3,
    XTrans = 4,
};

/// MARK: - Demosaic Variants

///----------------------------------------
/// @brief RGB bilinear interpolation demosaic.
/// @param img   Image array to demosaic (replaced in-place).
/// @param head  FITS header (updated: naxis=3, naxis3=3).
/// @param pattern Bayer pattern index (0-3).
///----------------------------------------

void demosaic_bilinear_interpolation(ImageArray& img, Header& head, int pattern);

///----------------------------------------
/// @brief Fuji X-Trans (6x6 cell) demosaic.
/// @param img  Image array to demosaic (replaced in-place).
/// @param head FITS header (updated: naxis=3, naxis3=3).
///----------------------------------------

void demosaic_x_trans(ImageArray& img, Header& head);

///----------------------------------------
/// @brief Han.k 2x2 spread demosaic, well-suited to oversampled astro frames.
/// @param img     Image array to demosaic (replaced in-place).
/// @param head    FITS header (updated: naxis=3, naxis3=3).
/// @param pattern Bayer pattern index (0-3).
///----------------------------------------

void demosaic_astrosimple(ImageArray& img, Header& head, int pattern);

///----------------------------------------
/// @brief Combined-spread variant of astrosimple (unused per original source).
/// @param img     Image array to demosaic (replaced in-place).
/// @param head    FITS header (updated: naxis=3, naxis3=3).
/// @param pattern Bayer pattern index (0-3).
///----------------------------------------

void demosaic_astrosimplebayercombined(ImageArray& img, Header& head, int pattern);

///----------------------------------------
/// @brief Mono-preserving bilinear interpolation demosaic.
/// @details Collapses pixels with steep magnitude slopes back to luminance.
/// @param img     Image array to demosaic (replaced in-place).
/// @param head    FITS header (updated: naxis=3, naxis3=3).
/// @param pattern Bayer pattern index (0-3).
///----------------------------------------

void demosaic_astroM_bilinear_interpolation(ImageArray& img, Header& head, int pattern);

///----------------------------------------
/// @brief Colour-preserving bilinear interpolation demosaic.
/// @details Detects saturation, reconstructs colour from a local circular
///          neighbourhood, and prevents purple stars.
/// @param img        Image array to demosaic (replaced in-place).
/// @param head       FITS header (updated: naxis=3, naxis3=3).
/// @param saturation Saturation threshold in ADU.
/// @param pattern    Bayer pattern index (0-3).
///----------------------------------------

void demosaic_astroC_bilinear_interpolation(ImageArray& img, Header& head,
                                            int saturation, int pattern);
                                            
///----------------------------------------
/// @brief 2x2 down-sampling super-pixel demosaic.
/// @param img     Image array to demosaic (replaced in-place, half dimensions).
/// @param head    FITS header (updated: width/2, height/2, naxis=3, naxis3=3).
/// @param pattern Bayer pattern index (0-3).
///----------------------------------------

void demosaic_superpixel(ImageArray& img, Header& head, int pattern);

///----------------------------------------
/// @brief In-place Bayer-mosaic saturation filler (single channel).
/// @param img  Image array to modify in-place.
/// @param head FITS header (read-only).
///----------------------------------------

void preserve_colour_saturated_bayer(ImageArray& img, const Header& head);

/// MARK: - Pattern Resolution and Dispatch

///----------------------------------------
/// @brief Resolve the effective Bayer pattern from a hint and FITS metadata.
/// @details Applies XBAYROFF / YBAYROFF / ROWORDER parity adjustments to the
///          pattern hint. Pass -1 to fall back to RGGB.
/// @param pattern_hint Raw resolved pattern (0-4) before offset flipping.
/// @param xbayroff     FITS XBAYROFF value.
/// @param ybayroff     FITS YBAYROFF value.
/// @param roworder     FITS ROWORDER string.
/// @return Final pattern index in {0..4}.
///----------------------------------------

[[nodiscard]] int get_demosaic_pattern(int pattern_hint,
                                      double xbayroff,
                                      double ybayroff,
                                      const std::string& roworder);
                                      
///----------------------------------------
/// @brief Demosaic method tag for the dispatcher.
///----------------------------------------

enum class DemosaicMethod {
    Bilinear,      ///< Default bilinear interpolation.
    AstroC,        ///< Colour-preserving variant.
    Simple,        ///< Han.k 2x2 spread.
    AstroM,        ///< Mono-preserving variant.
    Superpixel,    ///< 2x2 down-sample.
};

///----------------------------------------
/// @brief Dispatch to a demosaic variant based on the resolved pattern and method.
/// @details Pattern must already be the final 0-4 value (run through
///          get_demosaic_pattern first). X-Trans patterns always use the
///          X-Trans algorithm regardless of the method parameter.
/// @param img     Image array to demosaic (replaced in-place).
/// @param head    FITS header (updated on success).
/// @param pattern Final Bayer pattern index (0-4).
/// @param method  Demosaic algorithm to apply.
///----------------------------------------

void demosaic_bayer(ImageArray& img, Header& head, int pattern,
                    DemosaicMethod method);
                    
///----------------------------------------
/// @brief Top-level advanced demosaic wrapper.
/// @details Calls demosaic_bayer. Post-demosaic UI steps (auto-level,
///          colour-smooth, histogram refresh) are not yet ported.
/// @param img     Image array to demosaic (replaced in-place).
/// @param head    FITS header (updated on success).
/// @param pattern Final Bayer pattern index (0-4).
/// @param method  Demosaic algorithm to apply.
///----------------------------------------

void demosaic_advanced(ImageArray& img, Header& head, int pattern,
                       DemosaicMethod method);
    
} // namespace
