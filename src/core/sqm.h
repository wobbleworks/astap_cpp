///----------------------------------------
///      @file sqm.h
///   @ingroup ASTAP++
///     @brief Sky Quality Measurement (SQM) — zenith sky-surface-brightness and Bortle scale.
///   @details Calculates the zenith-equivalent sky-surface-brightness (magnitudes per
///            arcsec^2) from a photometrically-calibrated FITS frame and classifies
///            the result on the Bortle scale.
///    @author Ported from Han Kleijn's unit_sqm.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <string>

#include "../types.h"

///----------------------------------------
namespace astap::core {
///----------------------------------------

/// MARK: SQM Computation

///----------------------------------------
/// @brief Measure sky background and compute mag/arcsec^2 sky quality value.
/// @details Converts the photometric zero-point and pixel scale from the header
///          into a surface brightness, then applies an atmospheric-extinction
///          correction scaled to zenith.
/// @param head FITS header with MZERO, CDELT2, and pointing information.
/// @param get_backgr When true, recompute the sky background.
/// @param get_histogr When true, recompute the histogram.
/// @param[in,out] pedestal Dark pedestal value; zeroed when already subtracted or invalid.
/// @return True when a valid SQM value was produced and stored in @c sqmfloat.
///----------------------------------------

[[nodiscard]] bool calculate_sqm(Header& head, bool get_backgr, bool get_histogr, int& pedestal);

///----------------------------------------
/// @brief Map a mag/arcsec^2 value to the Bortle dark-sky classification string.
/// @param sqm Sky quality in magnitudes per arcsec^2.
/// @return Bortle classification string.
///----------------------------------------

[[nodiscard]] std::string bortle(double sqm) noexcept;

/// MARK: Globals

// Result values written by calculate_sqm and shared with the CLI / FITS-writer.

/// @brief Image-centre altitude as a formatted string.
extern std::string centalt;

/// @brief Image-centre altitude in degrees.
extern double altitudefloat;

/// @brief Last computed sky-quality value in mag/arcsec^2.
extern double sqmfloat;

// `airmass` is declared in globals.h (astap::airmass). SQM code writes to it
// too; no re-declaration here.

/// MARK: Forward Declarations

///----------------------------------------
/// @brief Pickering (2002) air-mass from apparent altitude in degrees.
/// @param altitude_deg Apparent altitude in degrees.
/// @return Air-mass value.
///----------------------------------------

[[nodiscard]] double airmass_calc(double altitude_deg) noexcept;

/// MARK: Site / Date Inputs

// Site and date inputs that must be populated before calculate_sqm runs.

/// @brief Site latitude (sexagesimal string).
extern std::string sitelat;

/// @brief Site longitude (sexagesimal string).
extern std::string sitelong;

/// @brief Persisted default latitude.
extern std::string lat_default;

/// @brief Persisted default longitude.
extern std::string long_default;

/// @brief Ambient temperature in Celsius.
extern double      temperature_c;

/// @brief Barometric pressure in hPa.
extern double      pressure_hpa;
	
} // namespace
