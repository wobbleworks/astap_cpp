///----------------------------------------
///      @file astrometric_solving.h
///   @ingroup ASTAP++
///     @brief Astrometric plate-solving entry point and related coordinate conversions.
///   @details Reads reference stars from a local star database, bins/crops the
///            image, finds stars in it, attempts to match its pattern against
///            the database, and computes a full WCS solution (optionally with
///            SIP polynomial distortion). Also provides coordinate-system
///            conversions (RA/DEC <-> CCD x/y) used by the rest of the app.
///    @author Ported from Han Kleijn's unit_astrometric_solving.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <string>
#include <vector>

#include "../types.h"
#include "star_align.h"

///----------------------------------------
namespace astap::solving {
///----------------------------------------

using astap::Header;

///----------------------------------------
/// MARK: Public API
///----------------------------------------

///----------------------------------------
/// @brief Main plate-solving entry point.
/// @details Finds a match between @p img and the star database and, on
///          success, fills the WCS fields in @p hd. @p log receives progress
///          and diagnostic messages.
/// @param img Image to solve.
/// @param hd Header to populate with the WCS solution.
/// @param log Receives human-readable progress messages.
/// @param get_hist Whether to include histogram metrics while measuring background.
/// @param check_patternfilter Whether to check the Bayer / OSC pattern filter.
/// @return @c true if a solution was found.
///----------------------------------------

[[nodiscard]] bool solve_image(ImageArray& img,
                               Header& hd,
                               std::vector<std::string>& log,
                               bool get_hist,
                               bool check_patternfilter);
                               
///----------------------------------------
///      @brief Bin/crop an image, measure the background, and detect stars.
///    @details On exit @p starlist3 holds the detected star (x, y) positions
///             in the coordinate system of the original (non-binned/cropped)
///             image. @p short_warning carries a brief advisory string suitable
///             for FITS headers (empty when none).
///      @param img Source image.
///      @param binning Downsampling factor (1, 2, 3, 4).
///      @param cropping Fractional crop (1.0 = no crop).
///      @param hfd_min Minimum half-flux diameter for detections.
///      @param hfd_max Maximum half-flux diameter for detections (10 matches the Pascal original; raise for DSS2 / wide-bloom plates).
///      @param max_stars Maximum number of stars to retain.
///      @param get_hist Whether to include histogram metrics while measuring background.
/// @param[out] starlist3 Detected star positions (2xN).
/// @param[out] short_warning Advisory message (may be empty).
///----------------------------------------

void bin_and_find_stars(const ImageArray& img,
                        int binning,
                        double cropping,
                        double hfd_min,
                        double hfd_max,
                        int max_stars,
                        bool get_hist,
                        StarList& starlist3,
                        std::string& short_warning);
                        
///----------------------------------------
///  @brief Select a downsampling factor for solving based on image height.
///  @param height Image height in pixels.
/// @return Binning factor to apply before solving.
///----------------------------------------

[[nodiscard]] int report_binning(double height);

///----------------------------------------
///  @brief Position angle of (ra1, dec1) as seen from (ra0, dec0).
/// @details Rigorous method (Meeus, @e Astronomical @e Algorithms, formula 46.5 / 48.5).
///  @param ra1 Target right ascension (radians).
///  @param dec1 Target declination (radians).
///  @param ra0 Reference right ascension (radians).
///  @param dec0 Reference declination (radians).
/// @return Position angle in radians.
///----------------------------------------

[[nodiscard]] double position_angle(double ra1, double dec1, double ra0, double dec0);

///----------------------------------------
///      @brief Transform equatorial coordinates (ra, dec) into CCD standard
///             coordinates (xx, yy) via a tangent-plane projection.
///      @param ra0 Projection centre right ascension (radians).
///      @param dec0 Projection centre declination (radians).
///      @param ra Target right ascension (radians).
///      @param dec Target declination (radians).
///      @param cdelt CCD scale in arcsec per pixel.
/// @param[out] xx Standard X coordinate.
/// @param[out] yy Standard Y coordinate.
///----------------------------------------

void equatorial_standard(double ra0, double dec0,
                         double ra, double dec,
                         double cdelt,
                         double& xx, double& yy);
                         
///----------------------------------------
///      @brief Read stars from the local star database around the telescope
///             pointing, covering @p search_field (radians).
///    @details On success, fills @p starlist (2xN: [0]=x, [1]=y, in standard
///             coordinates) and sets @p nrstars. Returns @c false on database
///             read failure.
///      @param telescope_ra Pointing right ascension (radians).
///      @param telescope_dec Pointing declination (radians).
///      @param search_field Search field size (radians).
///      @param database_type Database tile type (1, 290, 1476).
///      @param nrstars_required Target number of stars.
/// @param[out] starlist Star coordinates (2xN).
/// @param[out] nrstars Number of stars actually read.
///     @return @c true on success, @c false on database read failure.
///----------------------------------------

[[nodiscard]] bool read_stars(double telescope_ra,
                              double telescope_dec,
                              double search_field,
                              int database_type,
                              int nrstars_required,
                              StarList& starlist,
                              int& nrstars);
                              
///----------------------------------------
///      @brief Combine 2x2 pixel blocks into single mono samples, optionally cropping.
///      @param crop Fractional crop (1.0 = no crop).
///      @param img Source image.
/// @param[out] img2 Binned/cropped output (mono, single channel).
///----------------------------------------

void binX2_crop(double crop, const ImageArray& img, ImageArray& img2);

///----------------------------------------
///      @brief Crop image, make mono, no binning.
///      @param crop Fractional crop (1.0 = no crop).
///      @param img Source image.
/// @param[out] img2 Cropped output (mono, single channel).
///----------------------------------------

void binX1_crop(double crop, const ImageArray& img, ImageArray& img2);
                        
} // namespace
