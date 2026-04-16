///----------------------------------------
///      @file aberration.cpp
///   @ingroup ASTAP++
///     @brief Annual aberration and nutation corrections to equatorial coordinates.
///    @author Ported from Han Kleijn's unit_aberration.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "aberration.h"

#include <cmath>
#include <numbers>

#include "ephemerides.h"

///----------------------------------------
namespace astap::core {
///----------------------------------------

// External dependencies provided by other translated units.
double deltaT_calc(double jd);

/// MARK: - Anonymous Helpers

namespace {

/// @brief Fractional part (truncation toward zero).
[[nodiscard]] inline double frac(double x) noexcept {
    return x - std::trunc(x);
}

///----------------------------------------
/// @brief Transformation of mean to true coordinates (IAU 1980, terms > 0.1").
/// @param T Julian centuries from J2000: (JD - 2451545.0) / 36525.0.
/// @param[in,out] X Cartesian X.
/// @param[in,out] Y Cartesian Y.
/// @param[in,out] Z Cartesian Z.
///----------------------------------------

void nutequ(double T, double& X, double& Y, double& Z) {
    constexpr auto ARC = 206264.8062;   // arcseconds per radian
    constexpr auto P2  = 6.283185307;   // 2*pi
    
    auto LS  = P2 * frac(0.993133 +   99.997306 * T);
    auto D   = P2 * frac(0.827362 + 1236.853087 * T);
    auto F   = P2 * frac(0.259089 + 1342.227826 * T);
    auto N   = P2 * frac(0.347346 -    5.372447 * T);
    auto EPS = 0.4090928 - 2.2696e-4 * T;
    
    auto DPSI = ( -17.200 * std::sin(N)
                  - 1.319 * std::sin(2*(F - D + N))
                  - 0.227 * std::sin(2*(F + N))
                  + 0.206 * std::sin(2*N)
                  + 0.143 * std::sin(LS) ) / ARC;
    auto DEPS = ( +  9.203 * std::cos(N)
                  +  0.574 * std::cos(2*(F - D + N))
                  +  0.098 * std::cos(2*(F + N))
                  -  0.090 * std::cos(2*N) ) / ARC;
                  
    auto C = DPSI * std::cos(EPS);
    auto S = DPSI * std::sin(EPS);
    auto DX = -(C*Y + S*Z);
    auto DY =  (C*X - DEPS*Z);
    auto DZ =  (S*X + DEPS*Y);
    X += DX;
    Y += DY;
    Z += DZ;
}

///----------------------------------------
/// @brief Polar (r, theta, phi) to cartesian (x, y, z).
/// @param R Radius.
/// @param THETA Elevation in [-pi/2, pi/2].
/// @param PHI Azimuth.
/// @param[out] X Cartesian X.
/// @param[out] Y Cartesian Y.
/// @param[out] Z Cartesian Z.
///----------------------------------------

void cart2(double R, double THETA, double PHI,
           double& X, double& Y, double& Z) noexcept {
    auto sin_theta = std::sin(THETA);
    auto cos_theta = std::cos(THETA);
    auto sin_phi   = std::sin(PHI);
    auto cos_phi   = std::cos(PHI);
    auto RCST = R * cos_theta;
    X = RCST * cos_phi;
    Y = RCST * sin_phi;
    Z = R * sin_theta;
}

///----------------------------------------
/// @brief Cartesian (x, y, z) to polar (r, theta, phi).
/// @param x Cartesian X.
/// @param y Cartesian Y.
/// @param z Cartesian Z.
/// @param[out] r Radius.
/// @param[out] theta Elevation in [-pi/2, pi/2].
/// @param[out] phi Azimuth in [0, 2*pi).
///----------------------------------------

void polar2(double x, double y, double z,
            double& r, double& theta, double& phi) noexcept {
    auto rho = x*x + y*y;
    r = std::sqrt(rho + z*z);
    phi = std::atan2(y, x);
    if (phi < 0.0) {
        phi += 2.0 * std::numbers::pi;
    }
    rho = std::sqrt(rho);
    theta = std::atan2(z, rho);
}
    
}  // namespace

/// MARK: - Nutation

void nutation_correction_equatorial(double julian_et, double& ra, double& dec) {
    auto x0 = 0.0;
    auto y0 = 0.0;
    auto z0 = 0.0;
    cart2(1.0, dec, ra, x0, y0, z0);
    
    // Add nutation
    nutequ((julian_et - 2451545.0) / 36525.0, x0, y0, z0);
    
    auto r = 0.0;
    polar2(x0, y0, z0, r, dec, ra);
}

/// MARK: - Aberration

void aberration_correction_equatorial(double julian_et, double& ra, double& dec) {
    auto x0 = 0.0;
    auto y0 = 0.0;
    auto z0 = 0.0;
    cart2(1.0, dec, ra, x0, y0, z0);
    
    // Barycentric Earth velocity (AU/day). NOTE: the new API's
    // barycentric frame currently coincides with heliocentric (sub-0.01
    // AU offset); error in aberration is < 0.001 arcsec — acceptable.
    const auto earth = ephem::earth_state(julian_et, ephem::ReferenceFrame::Barycentric);

    // Convert AU/day to fraction of speed of light (~1/173)
    constexpr auto AU_PER_DAY_TO_C = 0.00577552;
    x0 += earth.velocity[0] * AU_PER_DAY_TO_C;
    y0 += earth.velocity[1] * AU_PER_DAY_TO_C;
    z0 += earth.velocity[2] * AU_PER_DAY_TO_C;
    
    auto r = 0.0;
    polar2(x0, y0, z0, r, dec, ra);
}

/// MARK: - J2000 to Apparent

void J2000_to_apparent(double jd, double& ra, double& dec) {
    // Julian Day based on dynamic time (TT = UTC + delta_T)
    auto jde = jd + deltaT_calc(jd);
    
    // Aberration in J2000 equinox
    aberration_correction_equatorial(jde, ra, dec);
    
    // Precession from J2000 to Jnow
    const auto precessed = ephem::precess_iau1976({ra, dec}, 2451545.0, jde);
    ra = precessed.ra;
    dec = precessed.dec;
    
    // Nutation (mean equinox, M&P page 125)
    nutation_correction_equatorial(jde, ra, dec);
}
    
} // namespace
