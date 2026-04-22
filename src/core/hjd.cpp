///----------------------------------------
///      @file hjd.cpp
///   @ingroup ASTAP++
///     @brief Time, sky-geometry, and coordinate-conversion helpers.
///    @author Ported from Han Kleijn's unit_hjd.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "hjd.h"

#include <chrono>
#include <cmath>
#include <numbers>
#include <string>

#include "../core/aberration.h"
#include "../core/ephemerides.h"
#include "../core/util.h"
#include "../core/wcs.h"
#include "../stacking/stack.h"
#include "../core/globals.h"

///----------------------------------------
namespace astap::core {
///----------------------------------------

/// MARK: External State

// Cross-module state populated by the settings loader before these functions
// are called. Pull the canonical `astap::` globals (globals.cpp) into this
// TU; local astap::core:: placeholders are only used for the ones without a
// canonical home.

using astap::focus_temp;
using astap::pressure;
using astap::site_lat_radians;
using astap::site_long_radians;

extern double wtime2actual;
extern std::string sitelat;
extern std::string sitelong;
extern std::string lat_default;
extern std::string long_default;
extern std::string centaz;
extern std::string centalt;

// Host-supplied log hook.
void memo2_message(std::string_view text);

/// MARK: - Anonymous Helpers

namespace {

///----------------------------------------
/// @brief Low-precision solar coordinates (~1 arc-minute).
/// @param jd Julian Date.
/// @param[out] ra Solar right ascension in radians.
/// @param[out] dec Solar declination in radians.
///----------------------------------------

void sun(double jd, double& ra, double& dec) {
    using std::numbers::pi;
    constexpr auto cos_ecl = 0.9174369450421712;  // cos(23.43929111 deg)
    constexpr auto sin_ecl = 0.3977769691456674;  // sin(23.43929111 deg)
    
    auto t = (jd - 2451545.0) / 36525.0;
    auto m = 2.0 * pi * (0.993133 + 99.997361 * t
                          - std::floor(0.993133 + 99.997361 * t));
    auto dl = 6893.0 * std::sin(m) + 72.0 * std::sin(2.0 * m);
    
    auto angle_raw = 0.7859453 + m / (2.0 * pi) + (6191.2 * t + dl) / 1296.0e3;
    auto angle = angle_raw - std::floor(angle_raw);
    if (angle < 0) {
        angle += 1.0;
    }
    
    auto l = 2.0 * pi * angle;
    auto sin_l = std::sin(l);
    auto cos_l = std::cos(l);
    auto y   = cos_ecl * sin_l;
    auto z   = sin_ecl * sin_l;
    auto rho = std::sqrt(1.0 - z * z);
    dec = std::atan(z / rho);
    ra  = 2.0 * std::atan(y / (cos_l + rho));
    if (ra < 0) {
        ra += 2.0 * pi;
    }
}

///----------------------------------------
/// @brief Equatorial to horizontal. Azimuth counted S -> W -> N -> E -> S.
/// @param dec Declination in radians.
/// @param tau Hour angle in radians.
/// @param phi Observer latitude in radians.
/// @param[out] h Altitude in radians.
/// @param[out] az Azimuth in radians.
///----------------------------------------

void equhor2(double dec, double tau, double phi, double& h, double& az) {
    auto sin_phi = std::sin(phi);
    auto cos_phi = std::cos(phi);
    auto sin_dec = std::sin(dec);
    auto cos_dec = std::cos(dec);
    auto sin_tau = std::sin(tau);
    auto cos_tau = std::cos(tau);
    auto x = cos_dec * sin_phi * cos_tau - sin_dec * cos_phi;
    auto y = cos_dec * sin_tau;
    auto z = cos_dec * cos_phi * cos_tau + sin_dec * sin_phi;
    auto dummy = 0.0;
    polar2(x, y, z, dummy, h, az);
}

///----------------------------------------
/// @brief Atmospheric refraction correction (Meeus 1991, 15.4).
/// @param altitude_real Apparent altitude in degrees.
/// @param p Barometric pressure in mbar.
/// @param t Temperature in Celsius.
/// @return Refraction correction in degrees.
///----------------------------------------

[[nodiscard]] double atmospheric_refraction(double altitude_real, double p, double t) noexcept {
    using std::numbers::pi;
    auto hn = (altitude_real + 10.3 / (altitude_real + 5.11)) * pi / 180.0;
    return ((p / 1010.0) * 283.0 / (273.0 + t)) * (1.02 / 60.0)
           / (std::sin(hn) / std::cos(hn));
}
  
}  // namespace

/// MARK: - Julian Date

double calc_jd_now() {
    using namespace std::chrono;
    auto now = system_clock::now();
    auto dp  = floor<days>(now);
    auto ymd = year_month_day{dp};
    auto hms = hh_mm_ss{now - dp};
    
    auto yy     = static_cast<int>(ymd.year());
    auto mm     = static_cast<unsigned>(ymd.month());
    auto dd     = static_cast<unsigned>(ymd.day());
    auto hour   = hms.hours().count();
    auto minute = hms.minutes().count();
    auto sec    = static_cast<double>(hms.seconds().count())
                + duration_cast<duration<double>>(hms.subseconds()).count();
                
    return astap::stacking::julian_calc(yy, mm, dd, hour, minute, sec);
}

/// MARK: - Heliocentric Correction

double JD_to_HJD(double jd, double ra_object, double dec_object) {
    auto ra_sun = 0.0;
    auto dec_sun = 0.0;
    sun(jd, ra_sun, dec_sun);
    
    // Precess solar coordinates to J2000 equinox
    const auto sun_j2k = ephem::precess_iau1976({ra_sun, dec_sun}, jd, 2451545.0);
    ra_sun = sun_j2k.ra;
    dec_sun = sun_j2k.dec;
    
    auto sin_dec_object = std::sin(dec_object);
    auto cos_dec_object = std::cos(dec_object);
    auto sin_dec_sun    = std::sin(dec_sun);
    auto cos_dec_sun    = std::cos(dec_sun);
    
    // 500 sec assumed light-travel time Sun - Earth
    return jd - 500.0 * (1.0 / (24.0 * 3600.0))
              * (sin_dec_object * sin_dec_sun
                 + cos_dec_object * cos_dec_sun
                   * std::cos(ra_object - ra_sun));
}

/// MARK: - Coordinate Conversion

void equ_gal(double ra, double dec, double& l, double& b) noexcept {
    using std::numbers::pi;
    
    // North-galactic-pole constants from Liu et al. 2011 (A&A)
    constexpr auto pole_ra  = (12.0 + 51.0 / 60.0 + 26.27549 / 3600.0) * std::numbers::pi / 12.0;
    constexpr auto pole_dec = (27.0 +  7.0 / 60.0 + 41.7043  / 3600.0) * std::numbers::pi / 180.0;
    constexpr auto posangle = (122.93191857 - 90.0) * std::numbers::pi / 180.0;
    
    auto sin_pole_dec = std::sin(pole_dec);
    auto cos_pole_dec = std::cos(pole_dec);
    b = std::asin(std::cos(dec) * cos_pole_dec * std::cos(ra - pole_ra)
                  + std::sin(dec) * sin_pole_dec);
    auto sin_b = std::sin(b);
    l = std::atan2(std::sin(dec) - sin_b * sin_pole_dec,
                   std::cos(dec) * cos_pole_dec * std::sin(ra - pole_ra))
        + posangle;
}

void polar2(double x, double y, double z,
            double& r, double& theta, double& phi) noexcept {
    using std::numbers::pi;
    auto rho = x * x + y * y;
    r = std::sqrt(rho + z * z);
    phi = std::atan2(y, x);
    if (phi < 0) {
        phi += 2.0 * pi;
    }
    rho = std::sqrt(rho);
    theta = std::atan2(z, rho);
}

void ra_az(double ra, double de, double lat, double longitude, double t,
           double& azimuth2, double& altitude2) noexcept {
    using std::numbers::pi;
    equhor2(de, ra - longitude - t, lat, altitude2, azimuth2);
    azimuth2 = pi - azimuth2;
    if (azimuth2 < 0) {
        azimuth2 += 2.0 * pi;
    }
}

void az_ra(double az, double alt, double lat, double longitude, double t,
           double& ra, double& dcr) noexcept {
    using std::numbers::pi;
    equhor2(alt, az, lat, dcr, ra);
    ra = pi - ra + longitude + t;
    while (ra <  0) {
        ra += 2.0 * pi;
    }
    while (ra >= 2.0 * pi) {
        ra -= 2.0 * pi;
    }
}

/// MARK: - Altitude and Refraction

void altitude_and_refraction(double lat, double longitude, double julian,
                             double temperature, double pressure_mbar,
                             double ra3, double dec3,
                             double& az, double& alt) {
    using std::numbers::pi;
    constexpr auto siderealtime2000      = 280.46061837 * std::numbers::pi / 180.0;
    constexpr auto earth_angular_velocity = std::numbers::pi * 2.0 * 1.00273790935;
    
    auto w = fnmodulo(
        +longitude + siderealtime2000
        + (julian - 2451545.0) * earth_angular_velocity,
        2.0 * pi);
    ra_az(ra3, dec3, lat, 0.0, w, az, alt);
    
    // Correct altitude for refraction; fall back to 10 C when unknown
    if (temperature >= 100.0) {
        temperature = 10.0;
    }
    az  = az  * 180.0 / pi;
    alt = alt * 180.0 / pi;
    alt += atmospheric_refraction(alt, pressure_mbar, temperature);
}

bool get_lat_long(double& site_lat_radians_out, double& site_long_radians_out) {
    using std::numbers::pi;
    
    if (sitelat.empty()) {
        sitelat  = lat_default;
        sitelong = long_default;
    }
    
    auto errordecode = false;
    
    // Latitude: try plain decimal degrees first, else dd mm ss
    auto lat_deg = strtofloat2(sitelat);
    if (lat_deg != 0.0 || sitelat == "0" || sitelat == "0.0") {
        site_lat_radians_out = lat_deg * pi / 180.0;
    } else {
        dec_text_to_radians(sitelat, site_lat_radians_out, errordecode);
    }
    
    if (!errordecode) {
        auto long_deg = strtofloat2(sitelong);
        if (long_deg != 0.0 || sitelong == "0" || sitelong == "0.0") {
            site_long_radians_out = long_deg * pi / 180.0;
        } else {
            dec_text_to_radians(sitelong, site_long_radians_out, errordecode);
        }
    }
    
    if (errordecode) {
        memo2_message("Error decoding site longitude or latitude!");
        return false;
    }
    return true;
}

void calculate_az_alt(int calc_mode, Header& head_in, double& az, double& alt) {
    az  = strtofloat2(centaz);
    alt = strtofloat2(centalt);
    
    if (!(((alt == 0) || (calc_mode > 0)) && (head_in.cd1_1 != 0))) {
        return;
    }
    
    auto local_lat = 0.0;
    auto local_long = 0.0;
    if (!get_lat_long(local_lat, local_long)) {
        return;
    }
    
    astap::stacking::date_to_jd(head_in.date_obs, head_in.date_avg, head_in.exposure);
    if (astap::jd_mid <= 2400000.0) {
        memo2_message("Error decoding Julian day!");
        return;
    }
    
    // Duplicate to protect J2000 position
    auto ra  = head_in.ra0;
    auto dec = head_in.dec0;
    if (calc_mode == 2) {
        J2000_to_apparent(astap::jd_mid, ra, dec);
    } else {
        const auto precessed = ephem::precess_iau1976({ra, dec}, 2451545.0, astap::jd_mid);
        ra = precessed.ra;
        dec = precessed.dec;
    }
    
    altitude_and_refraction(local_lat, local_long, astap::jd_mid,
                            focus_temp, pressure, ra, dec, az, alt);
}

bool calculate_az_alt_basic(double ra, double dec, double& az, double& alt) {
    using std::numbers::pi;
    constexpr auto siderealtime2000      = 280.46061837 * std::numbers::pi / 180.0;
    constexpr auto earth_angular_velocity = std::numbers::pi * 2.0 * 1.00273790935;
    
    if (site_lat_radians > 998.0) {
        if (!get_lat_long(site_lat_radians, site_long_radians)) {
            return false;
        }
        astap::stacking::date_to_jd(astap::head.date_obs,
                                    astap::head.date_avg,
                                    astap::head.exposure);
        if (astap::jd_mid > 2400000.0) {
            wtime2actual = fnmodulo(
                +site_long_radians + siderealtime2000
                + (astap::jd_mid - 2451545.0) * earth_angular_velocity,
                2.0 * pi);
        } else {
            memo2_message("Error decoding Julian day!");
            return false;
        }
    }
    
    {
        const auto precessed = ephem::precess_iau1976({ra, dec}, 2451545.0, astap::jd_mid);
        ra = precessed.ra;
        dec = precessed.dec;
    }
    ra_az(ra, dec, site_lat_radians, 0.0, wtime2actual, az, alt);
    
    auto temperature = (focus_temp >= 100.0) ? 10.0 : focus_temp;
    alt = alt * 180.0 / pi;
    alt = (alt + atmospheric_refraction(alt, pressure, temperature)) * pi / 180.0;
    return true;
}

/// MARK: - Air-mass and Extinction

double airmass_calc(double h) noexcept {
    using std::numbers::pi;
    if (h >= 1.0e-7) {
        return 1.0 / std::sin((pi / 180.0)
                              * (h + 244.0 / (165.0 + 47.0 * std::pow(h, 1.1))));
    }
    return 999.0;
}

double atmospheric_absorption(double airmass) noexcept {
    // Schaefer (1992) coefficients; total ~= 0.2811 * X
    auto a_ozon = airmass * 0.016;
    auto a_ray  = airmass * 0.1451;
    auto a_aer  = airmass * 0.120;
    return a_ozon + a_ray + a_aer;
}
  
} // namespace
