///----------------------------------------
///      @file online_gaia.h
///   @ingroup ASTAP++
///     @brief Online Gaia catalog query, parsing, and photometric passband conversion.
///   @details Queries the Gaia DR3 catalog via the VizieR TAP-like text endpoint,
///            parses the plain-text/CSV response, and populates an in-memory table
///            of star positions and magnitudes.  Also provides photometric-passband
///            transformations from Gaia G/BP/RP magnitudes to Johnson-Cousins
///            (B, V, R) and SDSS (g, r, i).
///    @author Ported from Han Kleijn's ASTAP (unit_online_gaia.pas). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <array>
#include <string>
#include <vector>

#include "../types.h"

///----------------------------------------
namespace astap::reference {
///----------------------------------------

/// MARK: Module-level state

/// Abstract HTTP client -- canonical definition lives in types.h.
using astap::IHttpClient;

/// In-memory star table.  Layout per star i:
///   [0][i] = RA  (radians)
///   [1][i] = Dec (radians)
///   [2][i] = G   magnitude
///   [3][i] = BP  magnitude
///   [4][i] = RP  magnitude
///   [5][i] = currently-active passband magnitude (starts as BP; overwritten
///            by @ref convert_magnitudes).
extern std::array<std::vector<double>, 6> online_database;

/// @brief Telescope pointing (radians) associated with the currently cached
///        @ref online_database; used by callers to decide whether the cache is
///        still valid for a new request.
extern double gaia_ra;

/// @brief Telescope declination (radians) paired with @ref gaia_ra.
extern double gaia_dec;

/// @brief Name of the passband currently stored in @c online_database[5].
///        Updated by @ref convert_magnitudes; consulted to avoid redundant
///        re-conversion.
extern std::string passband_active;

/// MARK: Public API

///----------------------------------------
///   @brief Download Gaia stars around a sky position and fill @ref online_database.
/// @details Queries VizieR for Gaia DR3 stars within a square box of side
///          @p search_field radians centred on (@p telescope_ra, @p telescope_dec),
///          limited to Gmag < @p magli.  On success sets @ref gaia_ra / @ref gaia_dec.
///   @param http           HTTP client used to perform the GET request.
///   @param telescope_ra   Right ascension of the pointing centre (radians).
///   @param telescope_dec  Declination of the pointing centre (radians).
///   @param search_field   Half-width of the search box (radians).
///   @param magli          Limiting G magnitude.
///  @return @c true on success; @c false on empty reply or HTTP failure.
///----------------------------------------

[[nodiscard]] bool read_stars_online(IHttpClient& http,
                                     double telescope_ra,
                                     double telescope_dec,
                                     double search_field,
                                     double magli);
                                     
///----------------------------------------
///   @brief Rewrite @c online_database[5] with magnitudes in the requested passband.
/// @details No-op when @p passband already equals @ref passband_active.
///   @param passband Target passband identifier (e.g. "BP", "V", "R", "B",
///                   "SG", "SR", "SI").
///----------------------------------------

void convert_magnitudes(const std::string& passband);

///----------------------------------------
///   @brief Photometric transformation from Gaia (G, BP, RP) to a named filter.
/// @details Supported filters: "BP" (pass-through), "V", "R", "B" (Johnson-Cousins),
///          "SG", "SR", "SI" (SDSS g/r/i).  Coefficients are from the Gaia
///          EDR3/DR3 calibration.
///   @param filter Filter name.
///   @param magG   Gaia G magnitude.
///   @param magBP  Gaia BP magnitude.
///   @param magRP  Gaia RP magnitude.
///  @return Transformed magnitude, or 0.0 if inputs are invalid or the colour
///          lies outside the polynomial's validity range.
///----------------------------------------

[[nodiscard]] double transform_gaia(const std::string& filter,
                                    double magG, double magBP, double magRP);
                                    
///----------------------------------------
///      @brief Look up the first star within 5 arcsec of a position and report
///             its transformed magnitudes.
///   @details All output parameters are set to 0 when no match is found.
///      @param ra  Right ascension to search (radians).
///      @param dec Declination to search (radians).
/// @param[out] b  Johnson-Cousins B magnitude.
/// @param[out] v  Johnson-Cousins V magnitude.
/// @param[out] r  Johnson-Cousins R magnitude.
/// @param[out] sg SDSS g magnitude.
/// @param[out] sr SDSS r magnitude.
/// @param[out] si SDSS i magnitude.
/// @param[out] g  Gaia G magnitude.
/// @param[out] bp Gaia BP magnitude.
/// @param[out] rp Gaia RP magnitude.
///----------------------------------------

void report_one_star_magnitudes(double ra, double dec,
                                double& b, double& v, double& r,
                                double& sg, double& sr, double& si,
                                double& g, double& bp, double& rp);
                                
} // namespace
