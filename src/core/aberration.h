///----------------------------------------
///      @file aberration.h
///   @ingroup ASTAP++
///     @brief Annual aberration and nutation corrections to equatorial coordinates.
///   @details Based on Montenbruck & Pfleger, "Astronomy on the Personal Computer".
///    @author Ported from Han Kleijn's unit_aberration.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

///----------------------------------------
namespace astap::core {
///----------------------------------------

///----------------------------------------
/// @brief Apply annual aberration in J2000 equinox.
/// @param julian_et Julian Ephemeris Time.
/// @param[in,out] ra Right ascension in radians.
/// @param[in,out] dec Declination in radians.
///----------------------------------------

void aberration_correction_equatorial(double julian_et, double& ra, double& dec);

///----------------------------------------
/// @brief Apply nutation to mean-equinox coordinates (M&P page 125).
/// @param julian_et Julian Ephemeris Time.
/// @param[in,out] ra Right ascension in radians.
/// @param[in,out] dec Declination in radians.
///----------------------------------------

void nutation_correction_equatorial(double julian_et, double& ra, double& dec);

///----------------------------------------
/// @brief Convert J2000 mean position to apparent position (no refraction).
/// @param jd Julian Date.
/// @param[in,out] ra Right ascension in radians.
/// @param[in,out] dec Declination in radians.
///----------------------------------------

void J2000_to_apparent(double jd, double& ra, double& dec);
	
} // namespace
