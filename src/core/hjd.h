///----------------------------------------
///      @file hjd.h
///   @ingroup ASTAP++
///     @brief Time, sky-geometry, and coordinate-conversion helpers.
///   @details Julian Date to Heliocentric Julian Date correction, equatorial to
///            galactic conversion, altitude/azimuth with atmospheric refraction,
///            air-mass, and atmospheric extinction estimation.
///    @author Ported from Han Kleijn's unit_hjd.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include "../types.h"

///----------------------------------------
namespace astap::core {
///----------------------------------------

/// MARK: Julian Date

///----------------------------------------
/// @brief Julian Day for the current UTC wall-clock moment.
///----------------------------------------

[[nodiscard]] double calc_jd_now();

///----------------------------------------
/// @brief Convert Julian Date to Heliocentric Julian Date for a target position.
/// @param jd Julian Date.
/// @param ra_object Right ascension of the target in radians.
/// @param dec_object Declination of the target in radians.
/// @return Heliocentric Julian Date.
///----------------------------------------

[[nodiscard]] double JD_to_HJD(double jd, double ra_object, double dec_object);

/// MARK: Coordinate Conversion

///----------------------------------------
/// @brief Equatorial (ra, dec in radians, J2000) to galactic (l, b in radians).
/// @param ra Right ascension in radians (J2000).
/// @param dec Declination in radians (J2000).
/// @param[out] l Galactic longitude in radians.
/// @param[out] b Galactic latitude in radians.
///----------------------------------------

void equ_gal(double ra, double dec, double& l, double& b) noexcept;

///----------------------------------------
/// @brief Cartesian (x, y, z) to polar (r, theta, phi).
/// @param x Cartesian x coordinate.
/// @param y Cartesian y coordinate.
/// @param z Cartesian z coordinate.
/// @param[out] r Radius.
/// @param[out] theta Elevation in [-pi/2, pi/2].
/// @param[out] phi Azimuth in [0, 2*pi).
///----------------------------------------

void polar2(double x, double y, double z,
            double& r, double& theta, double& phi) noexcept;
            
///----------------------------------------
/// @brief Azimuth/altitude to right ascension/declination (all radians).
/// @param az Azimuth in [0, 2*pi).
/// @param alt Altitude in [-pi/2, pi/2].
/// @param lat Observer latitude in radians.
/// @param longitude Observer longitude in radians.
/// @param t Local sidereal time in radians.
/// @param[out] ra Right ascension in radians.
/// @param[out] dcr Declination in radians.
///----------------------------------------

void az_ra(double az, double alt, double lat, double longitude, double t,
           double& ra, double& dcr) noexcept;
           
///----------------------------------------
/// @brief Right ascension/declination to azimuth/altitude (all radians).
/// @param ra Right ascension in radians.
/// @param de Declination in radians.
/// @param lat Observer latitude in radians.
/// @param longitude Observer longitude in radians.
/// @param t Local sidereal time in radians.
/// @param[out] azimuth2 Azimuth in radians.
/// @param[out] altitude2 Altitude in radians.
///----------------------------------------

void ra_az(double ra, double de, double lat, double longitude, double t,
           double& azimuth2, double& altitude2) noexcept;
           
/// MARK: Air-mass and Extinction

///----------------------------------------
/// @brief Pickering (2002) air-mass from apparent altitude in degrees.
/// @param h Apparent altitude in degrees.
/// @return Air-mass value, or 999 for altitudes at or below the horizon.
///----------------------------------------

[[nodiscard]] double airmass_calc(double h) noexcept;

///----------------------------------------
/// @brief Atmospheric extinction in magnitudes for the given air-mass.
/// @param airmass Air-mass (Schaefer 1992 / ICQ convention: 0.2811 mag per air-mass).
/// @return Extinction in magnitudes relative to vacuum.
///----------------------------------------

[[nodiscard]] double atmospheric_absorption(double airmass) noexcept;

/// MARK: Altitude and Refraction

///----------------------------------------
/// @brief Compute (az, alt) in degrees for the image's pointing centre.
/// @param calc_mode 0 = use header when present; 1 = force precession only;
///                  2 = force full apparent (aberration + nutation).
/// @param head FITS header with pointing and date information.
/// @param[out] az Azimuth in degrees.
/// @param[out] alt Altitude in degrees.
///----------------------------------------

void calculate_az_alt(int calc_mode, Header& head, double& az, double& alt);

///----------------------------------------
/// @brief Altitude/azimuth with atmospheric refraction correction.
/// @param lat Observer latitude in radians.
/// @param longitude Observer longitude in radians (positive east).
/// @param julian Julian day.
/// @param temperature Temperature in Celsius (>=100 treated as unknown).
/// @param pressure Barometric pressure in mbar.
/// @param ra3 Right ascension in radians (current equinox).
/// @param dec3 Declination in radians (current equinox).
/// @param[out] az Azimuth in degrees.
/// @param[out] alt Altitude in degrees (refraction-corrected).
///----------------------------------------

void altitude_and_refraction(double lat, double longitude, double julian,
                             double temperature, double pressure,
                             double ra3, double dec3,
                             double& az, double& alt);
                             
///----------------------------------------
/// @brief Basic precession + refraction calculation for mouse-pointer display.
/// @param ra Right ascension in radians (J2000).
/// @param dec Declination in radians (J2000).
/// @param[out] az Azimuth in radians.
/// @param[out] alt Altitude in radians (refraction-corrected).
/// @return False on error (e.g. site coordinates unavailable).
///----------------------------------------

[[nodiscard]] bool calculate_az_alt_basic(double ra, double dec, double& az, double& alt);

/// MARK: Site

///----------------------------------------
/// @brief Retrieve observation site latitude and longitude in radians from settings.
/// @param[out] site_lat_radians Site latitude in radians.
/// @param[out] site_long_radians Site longitude in radians.
/// @return False if the coordinate strings do not parse.
///----------------------------------------

[[nodiscard]] bool get_lat_long(double& site_lat_radians, double& site_long_radians);
           
} // namespace
