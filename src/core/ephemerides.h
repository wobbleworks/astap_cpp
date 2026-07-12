///----------------------------------------
///      @file ephemerides.h
///   @ingroup ASTAP++
///     @brief Ephemeris API for ASTAP++.
///   @details Provides heliocentric/barycentric Earth state, major-planet
///            positions, orbital-element propagation for asteroids and comets,
///            IAU 1976 precession, and cartesian↔spherical conversion. The
///            ephemeris/propagation routines are faithful ports of SLALIB;
///            precession follows Meeus, "Astronomical Algorithms" 2nd ed.
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
///  @details Backed by SLALIB sla_EPV2 (VSOP-simplified; @ref sla::epv2),
///           sub-arcsecond over roughly 1900-2100. The barycentric frame uses
///           EPV2's own SSB→Sun series, so no separate barycentre offset is
///           applied.
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
///    @brief Heliocentric J2000-equatorial state of a major planet.
///  @details Backed by SLALIB sla_PLANET (Simon 1994); accuracy is
///           arcsecond-level over a few centuries around J2000.
///    @param np Planet number 1..8 (Mercury..Neptune).
///    @param jd_tt Julian Date, Terrestrial Time.
///    @return Position (AU) and velocity (AU/day), J2000 equatorial. A zero
///            state is returned for an illegal planet number or numerical error.
///----------------------------------------

[[nodiscard]] State planet_state(int np, double jd_tt);

///----------------------------------------
/// MARK: - SLALIB universal-variable propagation
///----------------------------------------

///----------------------------------------
/// @namespace astap::core::ephem::sla
/// @brief Faithful port of SLALIB's universal-variable orbit routines
///        (P.T.Wallace; via Han Kleijn's Pascal, from pyslalib). These solve
///        the two-body problem with Stumpff universal variables, valid for
///        every conic (ellipse, parabola, hyperbola) without special-casing.
///        `propagate` above is a thin wrapper over `planel`.
///
///        Conventions: angles radians, distances AU. Input elements are
///        referred to the J2000 ecliptic & equinox; the returned position and
///        velocity are J2000 **equatorial** (mean equator/equinox of J2000),
///        with position in AU and velocity in **AU/second**. `date`/`epoch`
///        may be MJD or full JD — only their difference is used.
///----------------------------------------
namespace sla {

/// @brief Universal orbital elements (SLALIB U array). 1..13 → 0..12 here.
using UElements = std::array<double, 13>;

/// @brief Position (AU) + velocity (AU/s): SLALIB PV array.
using PV = std::array<double, 6>;

/// @brief sla_UE2PV — position/velocity from universal elements.
/// @param u in/out: entries [11],[12] are refreshed for fast re-prediction.
/// @return JSTAT: 0 OK, −1 null radius vector, −2 failed to converge.
int ue2pv(double date, UElements& u, PV& pv);

/// @brief sla_PV2UE — universal elements from a position/velocity state.
/// @param pmass companion mass in solar masses (0 for a test particle).
/// @return JSTAT: 0 OK, −1 negative pmass, −2 too close, −3 too slow.
int pv2ue(const PV& pv, double date, double pmass, UElements& u);

/// @brief sla_EL2UE — universal elements from conventional osculating elements.
/// @param jform 1 = major planet (a, mean longitude L, daily motion dm),
///              2 = minor planet (a, mean anomaly M), 3 = comet (q).
/// @return JSTAT: 0 OK, −1 illegal jform, −2 illegal e, −3 illegal aorq,
///                −4 illegal dm, −5 numerical error.
int el2ue(double date, int jform, double epoch, double orbinc, double anode,
          double perih, double aorq, double e, double aorl, double dm,
          UElements& u);

/// @brief sla_PLANEL — heliocentric J2000-equatorial PV (AU, AU/s) from
///        conventional elements. Wrapper: el2ue then ue2pv.
int planel(double date, int jform, double epoch, double orbinc, double anode,
           double perih, double aorq, double e, double aorl, double dm, PV& pv);

/// @brief sla_PLANET — approximate heliocentric J2000-equatorial PV (AU, AU/s)
///        of a major planet from Simon et al. 1994 mean elements + periodic
///        terms. @p np is 1..8 (Mercury..Neptune; Pluto removed).
/// @return JSTAT: 0 OK, 1 date outside ~1000 yr of J2000 (warning, still
///                usable), −1 illegal @p np, −2 numerical error.
int planet(int np, double date, PV& pv);

/// @brief sla_EPV2 — J2000 Earth position (AU) and velocity (AU/day),
///        heliocentric (@p bary = false) or barycentric (@p bary = true), in the
///        J2000 equatorial (ICRS) frame. VSOP-simplified series (Simon/Moisson,
///        after Bretagnon); sub-arcsecond over ~1900-2100. @p date may be MJD or
///        full JD (only its offset from J2000 is used). Unlike the other `sla`
///        routines the outputs are AU/**day**, matching this module's State.
void epv2(double date, bool bary, Vec3& pe, Vec3& ve);

}  // namespace sla

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
