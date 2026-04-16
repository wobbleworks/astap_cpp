///----------------------------------------
///      @file asteroid.h
///   @ingroup analysis
///     @brief Asteroid and comet ephemeris computation and catalog parsing.
///   @details Ported from unit_asteroid.pas. Contains the pure algorithmic
///            functions for computing minor-planet/comet positions, Delta T,
///            observer parallax, illumination geometry, magnitude corrections,
///            and MPC catalog line parsing. GUI/form code is not ported.
///            Original copyright (C) 2021 by Han Kleijn, www.hnsky.org.
///            Licensed under the Mozilla Public License, v. 2.0.
///    @author Created by John Stephen on 4/15/26.
/// @copyright Copyright 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <string>
#include <string_view>

#include "../core/ephemerides.h"
#include "../types.h"

///----------------------------------------
namespace astap::analysis {
///----------------------------------------

///----------------------------------------
/// @brief Orbital elements for an asteroid or comet.
/// @details For comets, a_M is set to 1e99 as a sentinel (mean anomaly is
///          meaningless for a comet whose epoch is the perihelion date).
///----------------------------------------

struct Asteroid {
	int yy{};                ///< Epoch year.
	int mm{};                ///< Epoch month.
	double dd{};             ///< Epoch day (fractional).
	double a_e{};            ///< Orbital eccentricity.
	double a_or_q{};         ///< Semimajor axis (asteroid, AU) or perihelion distance (comet, AU).
	double a_i{};            ///< Inclination (degrees).
	double a_ohm{};          ///< Longitude of ascending node (degrees).
	double a_w{};            ///< Argument of perihelion (degrees).
	double a_M{};            ///< Mean anomaly (degrees); 1e99 = comet sentinel.
	double H{};              ///< Absolute magnitude.
	double a_g{};            ///< Slope parameter (asteroid) or activity coeff (comet).
	std::string desn;        ///< Designation, up to 9 chars.
	std::string name;        ///< Name, up to 28 chars.
};

///----------------------------------------
/// @brief Compute Delta T (TT - UT1) in days for a given Julian Date.
/// @details Piecewise polynomial approximation covering years 2016--2999.
///          Returns ~60 seconds (as days) for dates outside that range.
///          This is the real implementation replacing the stub in link_stubs.cpp.
/// @param jd Julian Date.
/// @return Delta T in days.
///----------------------------------------

[[nodiscard]] double deltaT_calc(double jd) noexcept;

///----------------------------------------
/// @brief Observer parallax correction, converting geocentric to topocentric.
/// @details Adjusts the Cartesian position (x, y, z) in AU by subtracting the
///          observer's offset from the Earth centre. Uses WGS84-like flattening
///          (0.99664719), equatorial radius 6378.14 km, and a fixed height of
///          100 m above sea level. See Meeus, Astronomical Algorithms, p. 78.
/// @param wtime Local apparent sidereal time (radians).
/// @param latitude Observer geodetic latitude (radians).
/// @param[in,out] x Cartesian X in AU, corrected in place.
/// @param[in,out] y Cartesian Y in AU, corrected in place.
/// @param[in,out] z Cartesian Z in AU, corrected in place.
///----------------------------------------

void parallax_xyz(double wtime, double latitude,
                  double& x, double& y, double& z) noexcept;

///----------------------------------------
/// @brief Compute the apparent RA/Dec of a minor planet or comet.
/// @details Propagates the orbit to the observation epoch, applies light-time
///          and parallax corrections, and converts to equatorial coordinates.
///          Uses ephem::earth_state for the Earth heliocentric vector (unless
///          sun_earth_vector is true, indicating it was already computed),
///          ephem::propagate for the body state, and
///          ephem::cartesian_to_spherical for the final RA/Dec.
/// @param sun_earth_vector If true, skip recomputing the Earth heliocentric
///        vector (assumes it was already set by a prior call).
/// @param julian Dynamical time (JD, TT).
/// @param year Epoch year of the orbital elements.
/// @param month Epoch month.
/// @param day Epoch day (fractional).
/// @param a_e Eccentricity.
/// @param a_or_q Semimajor axis (AU) or perihelion distance (AU).
/// @param a_i Inclination (degrees).
/// @param a_ohm Longitude of ascending node (degrees).
/// @param a_w Argument of perihelion (degrees).
/// @param a_M Mean anomaly (degrees); 1e99 flags a comet.
/// @param wtime Local apparent sidereal time for parallax (radians).
/// @param site_lat Observer latitude for parallax (radians).
/// @param[out] ra3 Right ascension (radians).
/// @param[out] dec3 Declination (radians).
/// @param[out] delta Geocentric distance (AU).
/// @param[out] sun_delta Heliocentric distance (AU).
/// @param[out] outdated True if the observation is >120 days from the elements epoch.
/// @param[out] ph_earth Heliocentric Earth position (AU), written on first call.
/// @param[out] ph_pln Heliocentric planet position (AU), for subsequent illum2 calls.
///----------------------------------------

void minor_planet(bool sun_earth_vector, double julian,
                  int year, int month, double day,
                  double a_e, double a_or_q, double a_i,
                  double a_ohm, double a_w, double a_M,
                  double wtime, double site_lat,
                  double& ra3, double& dec3,
                  double& delta, double& sun_delta,
                  bool& outdated,
                  astap::core::ephem::Vec3& ph_earth,
                  astap::core::ephem::Vec3& ph_pln);

///----------------------------------------
/// @brief Illumination geometry in the Sun-Earth-planet triangle.
/// @details Computes distances, elongation, phase angle, and phase fraction
///          from heliocentric positions of the planet (x, y, z) and the Earth
///          (xe, ye, ze).
/// @param x  Heliocentric X of minor planet (AU).
/// @param y  Heliocentric Y of minor planet (AU).
/// @param z  Heliocentric Z of minor planet (AU).
/// @param xe Heliocentric X of Earth (AU).
/// @param ye Heliocentric Y of Earth (AU).
/// @param ze Heliocentric Z of Earth (AU).
/// @param[out] r_sp Distance Sun--planet (AU).
/// @param[out] r_ep Distance Earth--planet (AU).
/// @param[out] elong Elongation (degrees).
/// @param[out] phi Phase angle (degrees).
/// @param[out] phase Phase fraction (0--100).
///----------------------------------------

void illum2(double x, double y, double z,
            double xe, double ye, double ze,
            double& r_sp, double& r_ep,
            double& elong, double& phi, double& phase) noexcept;

///----------------------------------------
/// @brief Phase-dependent magnitude correction for asteroids.
/// @details Implements Meeus, Astronomical Algorithms, formula 32.14.
///          Returns the magnitude offset due to the phase angle.
/// @param g Slope parameter (typically 0.15).
/// @param b Phase angle in radians (Sun--asteroid--Earth).
/// @return Magnitude correction (additive).
///----------------------------------------

[[nodiscard]] double asteroid_magn_comp(double g, double b) noexcept;

///----------------------------------------
/// @brief Parse one line of an MPCORB.DAT file into orbital elements.
/// @details Handles the packed-epoch encoding (K=20xx century, letter encoding
///          for month/day). On failure, desn is set to an empty string.
/// @param txt The raw MPCORB.DAT line (>= ~200 characters expected).
/// @param[out] desn Designation (up to 7 chars), empty on parse failure.
/// @param[out] name Object name (up to 28 chars).
/// @param[out] yy Epoch year.
/// @param[out] mm Epoch month.
/// @param[out] dd Epoch day.
/// @param[out] a_e Eccentricity.
/// @param[out] a_a Semimajor axis (AU).
/// @param[out] a_i Inclination (degrees).
/// @param[out] a_ohm Longitude of ascending node (degrees).
/// @param[out] a_w Argument of perihelion (degrees).
/// @param[out] a_M Mean anomaly (degrees).
/// @param[out] h Absolute magnitude H.
/// @param[out] g Slope parameter G.
///----------------------------------------

void convert_MPCORB_line(std::string_view txt,
                         std::string& desn, std::string& name,
                         int& yy, int& mm,
                         double& dd, double& a_e, double& a_a,
                         double& a_i, double& a_ohm, double& a_w,
                         double& a_M, double& h, double& g);

///----------------------------------------
/// @brief Parse one line of an MPC comet-elements file into orbital elements.
/// @details On failure, desn is set to an empty string.
/// @param txt The raw comet-elements line.
/// @param[out] desn Designation (up to 9 chars), empty on parse failure.
/// @param[out] name Object name (up to 28 chars).
/// @param[out] yy Epoch year.
/// @param[out] mm Epoch month.
/// @param[out] dd Epoch day (fractional).
/// @param[out] ecc Eccentricity.
/// @param[out] q Perihelion distance (AU).
/// @param[out] inc2 Inclination (degrees).
/// @param[out] lan Longitude of ascending node (degrees).
/// @param[out] aop Argument of perihelion (degrees).
/// @param[out] M_anom Mean anomaly sentinel (always set to 1e99).
/// @param[out] H Absolute magnitude.
/// @param[out] k Comet activity parameter (G * 2.5).
///----------------------------------------

void convert_comet_line(std::string_view txt,
                        std::string& desn, std::string& name,
                        int& yy, int& mm,
                        double& dd, double& ecc, double& q,
                        double& inc2, double& lan, double& aop,
                        double& M_anom, double& H, double& k);

} // namespace
