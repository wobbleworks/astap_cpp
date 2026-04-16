///----------------------------------------
///      @file wcs.h
///   @ingroup ASTAP++
///     @brief World Coordinate System (WCS) and celestial-coordinate utilities.
///   @details Provides FITS pixel-to-celestial transforms, angular separation,
///            text parsers/formatters for RA/Dec, and FITS header I/O. SIP / DSS
///            distortion coefficients and related globals live in this namespace
///            as free variables.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP, MPL-2.0)
///            by John Stephen.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "../types.h"
#include "globals.h"

///----------------------------------------
namespace astap::core {
///----------------------------------------

using astap::Header;

/// MARK: Module-level state

/// @brief DSS subsampling factor consulted by @c pixel_to_celestial (formalism == 2).
extern int subsamp;

/// @brief Coordinate frame currently displayed by the host.
extern int coord_frame;

/// MARK: Coordinate transforms

///----------------------------------------
/// @brief Convert a FITS pixel position (1-based) to RA/Dec in radians.
/// @param head    FITS header with WCS solution.
/// @param fitsx   X pixel coordinate (1-based).
/// @param fitsy   Y pixel coordinate (1-based).
/// @param formalism Projection type: 0 = TAN, 1 = SIP, 2 = DSS.
/// @param[out] ra  Right ascension in radians.
/// @param[out] dec Declination in radians.
///----------------------------------------

void pixel_to_celestial(const Header& head,
                        double fitsx, double fitsy,
                        int formalism,
                        double& ra, double& dec);
                        
///----------------------------------------
/// @brief Convert RA/Dec (radians) to FITS pixel position (1-based).
/// @param head    FITS header with WCS solution.
/// @param ra      Right ascension in radians.
/// @param dec     Declination in radians.
/// @param[out] fitsX X pixel coordinate (1-based).
/// @param[out] fitsY Y pixel coordinate (1-based).
///----------------------------------------

void celestial_to_pixel(const Header& head,
                        double ra, double dec,
                        double& fitsX, double& fitsY);
                        
///----------------------------------------
/// @brief Convert CCD standard coordinates (tangent-plane) to equatorial.
/// @param ra0   Reference RA in radians.
/// @param dec0  Reference Dec in radians.
/// @param x     Tangent-plane X offset.
/// @param y     Tangent-plane Y offset.
/// @param cdelt Pixel scale in arcsec/px.
/// @param[out] ra  Computed RA in radians.
/// @param[out] dec Computed Dec in radians.
///----------------------------------------

void standard_equatorial2(double ra0, double dec0,
                          double x, double y, double cdelt,
                          double& ra, double& dec);
                          
///----------------------------------------
/// @brief Compute angular separation between two RA/Dec positions (Meeus 16.1).
/// @param ra1  First position RA in radians.
/// @param dec1 First position Dec in radians.
/// @param ra2  Second position RA in radians.
/// @param dec2 Second position Dec in radians.
/// @param[out] sep Angular separation in radians.
///----------------------------------------

void ang_sep(double ra1, double dec1, double ra2, double dec2, double& sep);

///----------------------------------------
/// @brief DSS astrometric position computation (external linkage).
/// @param x X pixel coordinate.
/// @param y Y pixel coordinate.
/// @param[out] ra  RA in radians.
/// @param[out] dec Dec in radians.
///----------------------------------------

void dsspos(double x, double y, double& ra, double& dec);

///----------------------------------------
/// @brief Convert equatorial coordinates to galactic.
/// @param ra  RA in radians.
/// @param dec Dec in radians.
/// @param[out] l Galactic longitude in radians.
/// @param[out] b Galactic latitude in radians.
///----------------------------------------

void EQU_GAL(double ra, double dec, double& l, double& b);

///----------------------------------------
/// @brief Compute basic azimuth and altitude from equatorial coordinates.
/// @param ra  RA in radians.
/// @param dec Dec in radians.
/// @param[out] az  Azimuth in radians.
/// @param[out] alt Altitude in radians.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool calculate_az_alt_basic(double ra, double dec, double& az, double& alt);

///----------------------------------------
/// @brief Compute position angle between two equatorial positions.
/// @param ra2  Second position RA in radians.
/// @param dec2 Second position Dec in radians.
/// @param ra1  First position RA in radians.
/// @param dec1 First position Dec in radians.
/// @return Position angle in radians.
///----------------------------------------

[[nodiscard]] double position_angle(double ra2, double dec2, double ra1, double dec1);

/// MARK: Higher-level helpers

///----------------------------------------
/// @brief Format a separation and position-angle pair between two FITS pixel positions.
/// @details Falls back to pixel distance if no astrometric solution (cdelt2 == 0).
/// @param head            FITS header with WCS solution.
/// @param polynomialIndex Formalism enum forwarded to @c pixel_to_celestial.
/// @param fitsx1          First position X pixel.
/// @param fitsy1          First position Y pixel.
/// @param fitsx2          Second position X pixel.
/// @param fitsy2          Second position Y pixel.
/// @param[out] seperation Formatted separation string.
/// @param[out] pa         Formatted position angle string.
///----------------------------------------

void ang_sep_two_positions(const Header& head, int polynomialIndex,
                           double fitsx1, double fitsy1,
                           double fitsx2, double fitsy2,
                           std::string& seperation, std::string& pa);
                           
/// MARK: Text parsers / formatters

///----------------------------------------
/// @brief Parse coordinate text in various formats into RA/Dec radians.
/// @details Accepts "12h34m56.7s +12d34m56.7s", Simbad, and other common forms.
///          When input is "C"/"c", computes the image centre via @c pixel_to_celestial.
/// @param data0           Input coordinate string.
/// @param head            FITS header (used for image-centre fallback).
/// @param polynomialIndex Formalism enum for pixel_to_celestial.
/// @param[out] ra4        Parsed RA in radians.
/// @param[out] dec4       Parsed Dec in radians.
/// @return True on successful parse.
///----------------------------------------

[[nodiscard]] bool decode_string(std::string data0,
                                 const Header& head, int polynomialIndex,
                                 double& ra4, double& dec4);
                                 
///----------------------------------------
/// @brief Parse RA text in any accepted format to radians.
/// @param inp      Input RA string.
/// @param[out] ra  Parsed RA in radians.
/// @param[out] errorRA True if parsing failed.
///----------------------------------------

void ra_text_to_radians(std::string inp, double& ra, bool& errorRA);

///----------------------------------------
/// @brief Parse Dec text in any accepted format to radians.
/// @param inp       Input Dec string.
/// @param[out] dec  Parsed Dec in radians.
/// @param[out] errorDEC True if parsing failed.
///----------------------------------------

void dec_text_to_radians(std::string inp, double& dec, bool& errorDEC);

///----------------------------------------
/// @brief Format a position with the requested separator (frame-aware).
/// @param sep Separator string between RA and Dec portions.
/// @param ra  RA in radians.
/// @param dec Dec in radians.
/// @return Formatted position string.
///----------------------------------------

[[nodiscard]] std::string position_to_string(std::string_view sep, double ra, double dec);

///----------------------------------------
/// @brief Format RA as "HH:MM.D" with the given separator.
/// @param rax RA in radians.
/// @param sep Separator between hours and minutes.
/// @return Formatted string.
///----------------------------------------

[[nodiscard]] std::string prepare_ra5(double rax, std::string_view sep);

///----------------------------------------
/// @brief Format Dec as "+/-DD:MM" with the given separator.
/// @param decx Dec in radians.
/// @param sep  Separator between degrees and minutes.
/// @return Formatted string.
///----------------------------------------

[[nodiscard]] std::string prepare_dec4(double decx, std::string_view sep);

///----------------------------------------
/// @brief Format RA as "HH MM SS" with the given separator.
/// @param rax RA in radians.
/// @param sep Separator between hours and minutes.
/// @return Formatted string.
///----------------------------------------

[[nodiscard]] std::string prepare_ra6(double rax, std::string_view sep);

///----------------------------------------
/// @brief Format RA as "HH MM SS.D" with the given separator.
/// @param rax RA in radians.
/// @param sep Separator between hours and minutes.
/// @return Formatted string.
///----------------------------------------

[[nodiscard]] std::string prepare_ra(double rax, std::string_view sep);

///----------------------------------------
/// @brief Format Dec as "+/-DD MM SS" with the given separator.
/// @param decx Dec in radians.
/// @param sep  Separator between degrees and minutes.
/// @return Formatted string.
///----------------------------------------

[[nodiscard]] std::string prepare_dec(double decx, std::string_view sep);

///----------------------------------------
/// @brief Format RA as "HH MM SS.DD" with the given separator.
/// @param rax RA in radians.
/// @param sep Separator between hours and minutes.
/// @return Formatted string.
///----------------------------------------

[[nodiscard]] std::string prepare_ra8(double rax, std::string_view sep);

///----------------------------------------
/// @brief Format Dec as "+/-DD MM SS.D" with the given separator.
/// @param decx Dec in radians.
/// @param sep  Separator between degrees and minutes.
/// @return Formatted string.
///----------------------------------------

[[nodiscard]] std::string prepare_dec2(double decx, std::string_view sep);

///----------------------------------------
/// @brief Format IAU designation in HHMMSS.s+DDMMSS form.
/// @param rax  RA in radians.
/// @param decx Dec in radians.
/// @return IAU designation string.
///----------------------------------------

[[nodiscard]] std::string prepare_IAU_designation(double rax, double decx);

/// MARK: File / time helpers

///----------------------------------------
/// @brief Convert a Julian Date to MPC date string ("YYYY MM DD.DDDDD").
/// @param jd Julian Date.
/// @return Formatted MPC date string.
///----------------------------------------

[[nodiscard]] std::string Jd_To_MPCDate(double jd);

///----------------------------------------
/// @brief Write FITS header lines to a file as a valid FITS header.
/// @details Writes 80-char records padded to a 2880-byte multiple.
/// @param filen       Output file path.
/// @param headerLines FITS header lines to write.
///----------------------------------------

void write_astronomy_wcs(const std::filesystem::path& filen,
                         const std::vector<std::string>& headerLines);
                        
} // namespace
