///----------------------------------------
///      @file ephemerides.h
///   @ingroup ASTAP++
///     @brief Clean-room ephemeris API for ASTAP++.
///   @details Provides heliocentric/barycentric Earth state, orbital-element
///            propagation for asteroids and comets, IAU 1976 precession, and
///            cartesian↔spherical coordinate conversion. Algorithms follow
///            Meeus, "Astronomical Algorithms" 2nd ed.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <array>
#include <cstddef>

///----------------------------------------
namespace astap::core::ephem {
///----------------------------------------

///----------------------------------------
/// MARK: - Value types
///----------------------------------------

/// @brief 3-element Cartesian vector. Units are context-dependent.
using Vec3 = std::array<double, 3>;

///----------------------------------------
/// @brief Heliocentric (or barycentric) state: position in AU, velocity in
///        AU/day, J2000 mean equator & equinox.
///----------------------------------------

struct State {
	Vec3 position{};
	Vec3 velocity{};
};

///----------------------------------------
/// @brief Spherical equatorial coordinates in radians.
///----------------------------------------

struct EquatorialCoords {
	double ra{};    ///< Right ascension, radians, normalised to [0, 2π).
	double dec{};   ///< Declination, radians.
};

///----------------------------------------
/// @brief Spherical conversion result: distance + equatorial pair.
///----------------------------------------

struct SphericalResult {
	double radius{};   ///< Distance (same units as input Cartesian).
	double dec{};      ///< Declination, radians.
	double ra{};       ///< Right ascension, radians, in [0, 2π).
};

///----------------------------------------
/// MARK: - Orbital elements
///----------------------------------------

///----------------------------------------
/// @brief Osculating orbital elements.
///  @details A single, unified element set that works for all orbital shapes:
///           elliptical (e < 1), parabolic (e == 1), and hyperbolic (e > 1).
///           The mean anomaly at epoch should be set to zero for a comet whose
///           epoch is the perihelion passage.
///----------------------------------------

struct OrbitalElements {
	double epoch_jd{};                   ///< Epoch of elements (Julian Date, TT).
	double perihelion_distance{};        ///< q (AU). Unambiguous for all eccentricities.
	double eccentricity{};               ///< e (≥ 0).
	double inclination{};                ///< i (radians).
	double ascending_node{};             ///< Ω, longitude of ascending node (radians).
	double argument_of_perihelion{};     ///< ω (radians).
	double mean_anomaly_at_epoch{};      ///< M at epoch (radians); zero for comets at perihelion.
};

///----------------------------------------
/// MARK: - Earth ephemeris
///----------------------------------------

/// @brief Which frame to report the Earth state in.
enum class ReferenceFrame {
	Heliocentric,   ///< Sun-to-Earth vector.
	Barycentric,    ///< Solar-System-Barycentre-to-Earth vector.
};

///----------------------------------------
///    @brief J2000 Earth position and velocity at a given epoch.
///  @details Uses the Meeus Ch. 25 low-precision Sun/Earth ephemeris
///           (~40 arcsec accuracy over 2000 ± 50 yr). For the barycentric
///           frame, the Sun-barycentre offset is approximated as zero,
///           adding up to ~2 solar-radii of error (≈ 0.009 AU). Both are
///           comfortably within ASTAP's astrometric needs, but for
///           sub-arcsecond work replace with the VSOP87 evaluator below
///           and a full barycentric correction.
///    @param jd_tt Julian Date, Terrestrial Time.
///    @param frame Heliocentric (default) or Barycentric.
///    @return Position (AU) and velocity (AU/day) in J2000 equatorial frame.
///----------------------------------------

[[nodiscard]] State earth_state(double jd_tt,
                                ReferenceFrame frame = ReferenceFrame::Heliocentric);

///----------------------------------------
/// MARK: - Orbit propagation
///----------------------------------------

///----------------------------------------
///    @brief Propagate osculating elements to a new epoch.
///  @details Solves Kepler's equation for the appropriate conic section
///           (elliptic, parabolic, or hyperbolic) and returns the
///           heliocentric state in the J2000 equatorial frame.
///    @param elements Osculating elements at some reference epoch.
///    @param jd_tt Target epoch, Julian Date TT.
///    @return Heliocentric state (AU, AU/day).
///----------------------------------------

[[nodiscard]] State propagate(const OrbitalElements& elements, double jd_tt);

///----------------------------------------
/// MARK: - Precession
///----------------------------------------

///----------------------------------------
///    @brief Precess equatorial coordinates between two epochs (IAU 1976 / FK5).
///  @details Uses the rigorous three-angle formulation from Meeus 21.2 /
///           21.3. Result RA is normalised to [0, 2π).
///    @param coords Source coordinates.
///    @param jd_from Source epoch (Julian Date).
///    @param jd_to Target epoch (Julian Date).
///    @return Precessed coordinates.
///----------------------------------------

[[nodiscard]] EquatorialCoords precess_iau1976(EquatorialCoords coords,
                                               double jd_from, double jd_to);

///----------------------------------------
/// MARK: - Coordinate conversions
///----------------------------------------

///----------------------------------------
///    @brief Convert a Cartesian vector to (radius, dec, ra).
///    @param xyz Cartesian vector.
///    @return Radius (same units), declination (rad), right ascension (rad, [0, 2π)).
///----------------------------------------

[[nodiscard]] SphericalResult cartesian_to_spherical(const Vec3& xyz) noexcept;

///----------------------------------------
///    @brief Convert (radius, dec, ra) to a Cartesian vector.
///----------------------------------------

[[nodiscard]] Vec3 spherical_to_cartesian(double radius, double dec, double ra) noexcept;

///----------------------------------------
/// MARK: - VSOP87-style evaluator (extension point)
///----------------------------------------

///----------------------------------------
/// @namespace astap::core::ephem::vsop
/// @brief A minimal VSOP87-style series evaluator. Data tables are out
///        of scope for this prototype — wire them in from Meeus
///        Appendix III or an IMCCE file when sub-arcsecond Earth/planet
///        positions are needed.
///----------------------------------------
namespace vsop {

/// @brief One VSOP series term: A·cos(B + C·τ).
struct Term {
	double amplitude;   ///< A
	double phase;       ///< B (radians)
	double frequency;   ///< C (radians / Julian millennium)
};

/// @brief A list of VSOP terms for one power of time.
struct Series {
	const Term* terms{};
	std::size_t count{};
};

/// @brief Multiple time-powers for one coordinate component (L, B, or R).
struct Component {
	const Series* series{};
	std::size_t count{};
};

/// @brief Complete VSOP87 data set for one body: longitude, latitude, radius.
struct Body {
	Component L;   ///< Ecliptic longitude.
	Component B;   ///< Ecliptic latitude.
	Component R;   ///< Heliocentric radius.
};

///----------------------------------------
/// @brief Ecliptic heliocentric coordinates produced by a VSOP evaluation.
///----------------------------------------

struct HelioEcliptic {
	double longitude;   ///< Ecliptic longitude (rad, normalised).
	double latitude;    ///< Ecliptic latitude (rad).
	double radius;      ///< Heliocentric radius (AU).
};

///----------------------------------------
///    @brief Evaluate a VSOP87 series for a given body at a given epoch.
///    @param body VSOP term tables for the body.
///    @param jd_tt Julian Date, Terrestrial Time.
///    @return Heliocentric ecliptic longitude, latitude, and radius.
///----------------------------------------

[[nodiscard]] HelioEcliptic evaluate(const Body& body, double jd_tt) noexcept;

}  // namespace vsop

}  // namespace astap::core::ephem
