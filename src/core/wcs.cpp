///----------------------------------------
///      @file wcs.cpp
///   @ingroup ASTAP++
///     @brief World Coordinate System (WCS) and celestial-coordinate utilities.
///   @details Implementation of FITS pixel-to-celestial transforms, angular
///            separation, text parsers/formatters, and FITS header I/O.
///            Order of operations is preserved for numerical parity with the
///            original ASTAP source.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP, MPL-2.0)
///            by John Stephen.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "wcs.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <format>
#include <fstream>
#include <numbers>
#include <string>
#include <string_view>
#include <vector>

///----------------------------------------
namespace astap::core {
///----------------------------------------

namespace {

constexpr auto kPi = std::numbers::pi;

/// MARK: Ported helper functions

// Simultaneous sine and cosine.
inline void sincos(double a, double& s, double& c) {
    s = std::sin(a);
    c = std::cos(a);
}

// Square of a value.
[[nodiscard]] inline double sqr(double v) noexcept {
    return v * v;
}

// atan2 wrapper.
[[nodiscard]] inline double arctan2(double y, double x) noexcept {
    return std::atan2(y, x);
}

// atan wrapper.
[[nodiscard]] inline double arctan(double v) noexcept {
    return std::atan(v);
}

// asin wrapper.
[[nodiscard]] inline double arcsin(double v) noexcept {
    return std::asin(v);
}

// acos wrapper.
[[nodiscard]] inline double arccos(double v) noexcept {
    return std::acos(v);
}

// Fractional part with sign.
[[nodiscard]] inline double frac(double v) noexcept {
    return v - std::trunc(v);
}

// Integer part toward zero.
[[nodiscard]] inline double pas_int(double v) noexcept {
    return std::trunc(v);
}

// Integer toward zero as long long.
[[nodiscard]] inline long long pas_trunc(double v) noexcept {
    return static_cast<long long>(std::trunc(v));
}

// At least 2 digits with leading '0' if needed.
[[nodiscard]] std::string leading_zero(long long w) {
    auto s = std::to_string(w);
    if (s.size() == 1) {
        s = "0" + s;
    }
    return s;
}

// 1-based position of sub in s, 0 if not found.
[[nodiscard]] int pos_str(std::string_view sub, std::string_view s) noexcept {
    auto p = s.find(sub);
    if (p == std::string_view::npos) {
        return 0;
    }
    return static_cast<int>(p) + 1;
}

// 1-based search from index start (1-based).
[[nodiscard]] int posex_str(std::string_view sub, std::string_view s, int start) noexcept {
    if (start < 1) {
        start = 1;
    }
    if (static_cast<size_t>(start - 1) > s.size()) {
        return 0;
    }
    auto p = s.find(sub, static_cast<size_t>(start - 1));
    if (p == std::string_view::npos) {
        return 0;
    }
    return static_cast<int>(p) + 1;
}

// 1-based substring copy.
[[nodiscard]] std::string pas_copy(const std::string& s, int start, int count) {
    if (start < 1) {
        start = 1;
    }
    auto i = static_cast<size_t>(start - 1);
    if (i >= s.size() || count <= 0) {
        return {};
    }
    auto n = std::min<size_t>(count, s.size() - i);
    return s.substr(i, n);
}

// 1-based in-place deletion.
void pas_delete(std::string& s, int start, int count) {
    if (start < 1 || count <= 0) {
        return;
    }
    auto i = static_cast<size_t>(start - 1);
    if (i >= s.size()) {
        return;
    }
    auto n = std::min<size_t>(count, s.size() - i);
    s.erase(i, n);
}

// ASCII uppercase.
[[nodiscard]] std::string uppercase(std::string s) {
    for (auto& c : s) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    }
    return s;
}

// Replace first occurrence only.
[[nodiscard]] std::string string_replace_first(std::string s, std::string_view from, std::string_view to) {
    auto p = s.find(from);
    if (p == std::string::npos) {
        return s;
    }
    s.replace(p, from.size(), to);
    return s;
}

// Replace all occurrences.
[[nodiscard]] std::string string_replace_all(std::string s, std::string_view from, std::string_view to) {
    if (from.empty()) {
        return s;
    }
    size_t p = 0;
    while ((p = s.find(from, p)) != std::string::npos) {
        s.replace(p, from.size(), to);
        p += to.size();
    }
    return s;
}

// Trim ASCII whitespace.
[[nodiscard]] std::string trim(std::string s) {
    size_t a = 0;
    while (a < s.size() && static_cast<unsigned char>(s[a]) <= 0x20) {
        ++a;
    }
    auto b = s.size();
    while (b > a && static_cast<unsigned char>(s[b - 1]) <= 0x20) {
        --b;
    }
    return s.substr(a, b - a);
}

// Parse a real number; err == 0 on success.
void pas_val(const std::string& s, double& x, int& err) {
    if (s.empty()) {
        x = 0;
        err = 1;
        return;
    }
    try {
        size_t consumed = 0;
        auto copy = s;
        x = std::stod(copy, &consumed);
        
        // Skip trailing whitespace; non-empty tail counts as failure.
        while (consumed < copy.size() && static_cast<unsigned char>(copy[consumed]) <= 0x20) {
            ++consumed;
        }
        err = (consumed == copy.size()) ? 0 : static_cast<int>(consumed) + 1;
    } catch (...) {
        x = 0;
        err = 1;
    }
}

// Fixed-point with specified decimal places.
[[nodiscard]] std::string floatstr_fixed(double v, int digits) {
    return std::format("{:.{}f}", v, digits);
}

// Integer formatting via std::format.
template <class T>
[[nodiscard]] std::string fmt_int(T v) {
    return std::format("{}", v);
}

// Build a string of n spaces.
[[nodiscard]] std::string spaces(size_t n) {
    return std::string(n, ' ');
}
    
}  // namespace

/// MARK: SIP / DSS globals

int a_order = 0;

double a_0_0 = 0, a_0_1 = 0, a_0_2 = 0, a_0_3 = 0;
double a_1_0 = 0, a_1_1 = 0, a_1_2 = 0;
double a_2_0 = 0, a_2_1 = 0;
double a_3_0 = 0;
double b_0_0 = 0, b_0_1 = 0, b_0_2 = 0, b_0_3 = 0;
double b_1_0 = 0, b_1_1 = 0, b_1_2 = 0;
double b_2_0 = 0, b_2_1 = 0;
double b_3_0 = 0;

double ap_0_0 = 0, ap_0_1 = 0, ap_0_2 = 0, ap_0_3 = 0;
double ap_1_0 = 0, ap_1_1 = 0, ap_1_2 = 0;
double ap_2_0 = 0, ap_2_1 = 0;
double ap_3_0 = 0;
double bp_0_0 = 0, bp_0_1 = 0, bp_0_2 = 0, bp_0_3 = 0;
double bp_1_0 = 0, bp_1_1 = 0, bp_1_2 = 0;
double bp_2_0 = 0, bp_2_1 = 0;
double bp_3_0 = 0;

bool sip = false;
int subsamp = 1;
int coord_frame = 0;

/// MARK: celestial_to_pixel / pixel_to_celestial

void celestial_to_pixel(const Header& head, double ra, double dec,
                        double& fitsX, double& fitsY) {
    auto SIN_dec = 0.0, COS_dec = 0.0;
    auto SIN_dec_ref = 0.0, COS_dec_ref = 0.0;
    auto SIN_delta_ra = 0.0, COS_delta_ra = 0.0;
    
    sincos(dec, SIN_dec, COS_dec);
    sincos(head.dec0, SIN_dec_ref, COS_dec_ref);
    sincos(ra - head.ra0, SIN_delta_ra, COS_delta_ra);
    
    auto H  = SIN_dec * SIN_dec_ref + COS_dec * COS_dec_ref * COS_delta_ra;
    auto xi = (COS_dec * SIN_delta_ra / H) * 180.0 / kPi;
    auto eta = ((SIN_dec * COS_dec_ref - COS_dec * SIN_dec_ref * COS_delta_ra) / H) * 180.0 / kPi;
    
    auto det = head.cd2_2 * head.cd1_1 - head.cd1_2 * head.cd2_1;
    
    auto u0 = -(head.cd1_2 * eta - head.cd2_2 * xi) / det;
    auto v0 = +(head.cd1_1 * eta - head.cd2_1 * xi) / det;
    
    if (sip) {
        fitsX = (head.crpix1 + u0
                 + ap_0_0 + ap_0_1 * v0 + ap_0_2 * v0 * v0 + ap_0_3 * v0 * v0 * v0
                 + ap_1_0 * u0 + ap_1_1 * u0 * v0 + ap_1_2 * u0 * v0 * v0
                 + ap_2_0 * u0 * u0 + ap_2_1 * u0 * u0 * v0
                 + ap_3_0 * u0 * u0 * u0);
        fitsY = (head.crpix2 + v0
                 + bp_0_0 + bp_0_1 * v0 + bp_0_2 * v0 * v0 + bp_0_3 * v0 * v0 * v0
                 + bp_1_0 * u0 + bp_1_1 * u0 * v0 + bp_1_2 * u0 * v0 * v0
                 + bp_2_0 * u0 * u0 + bp_2_1 * u0 * u0 * v0
                 + bp_3_0 * u0 * u0 * u0);
    } else {
        fitsX = head.crpix1 + u0;
        fitsY = head.crpix2 + v0;
    }
}

void pixel_to_celestial(const Header& head, double fitsx, double fitsy,
                        int formalism, double& ra, double& dec) {
    // Default for wrong index or head.cd1_1 == 0
    ra = 0;
    dec = 0;
    
    if (formalism == 2) {
        // DSS survey
        auto fits_unsampledX = subsamp * (fitsx - 0.5) + 0.5;
        auto fits_unsampledY = subsamp * (fitsy - 0.5) + 0.5;
        dsspos(fits_unsampledX, fits_unsampledY, ra, dec);
    } else if (head.cd1_1 != 0) {
        auto u2 = 0.0, v2 = 0.0;
        
        if ((formalism == 1) && (a_order >= 2)) {
            auto u = fitsx - head.crpix1;
            auto v = fitsy - head.crpix2;
            u2 = u + a_0_0 + a_0_1 * v + a_0_2 * v * v + a_0_3 * v * v * v
                 + a_1_0 * u + a_1_1 * u * v + a_1_2 * u * v * v
                 + a_2_0 * u * u + a_2_1 * u * u * v + a_3_0 * u * u * u;
            v2 = v + b_0_0 + b_0_1 * v + b_0_2 * v * v + b_0_3 * v * v * v
                 + b_1_0 * u + b_1_1 * u * v + b_1_2 * u * v * v
                 + b_2_0 * u * u + b_2_1 * u * u * v + b_3_0 * u * u * u;
        } else {
            u2 = fitsx - head.crpix1;
            v2 = fitsy - head.crpix2;
        }
        
        auto xi  = (head.cd1_1 * (u2) + head.cd1_2 * (v2)) * kPi / 180.0;
        auto eta = (head.cd2_1 * (u2) + head.cd2_2 * (v2)) * kPi / 180.0;
        
        auto sindec0 = 0.0, cosdec0 = 0.0;
        sincos(head.dec0, sindec0, cosdec0);
        auto delta = cosdec0 - eta * sindec0;
        ra = head.ra0 + arctan2(xi, delta);
        dec = arctan((sindec0 + eta * cosdec0) / std::sqrt(sqr(xi) + sqr(delta)));
        
        if (ra < 0) {
            ra = ra + kPi * 2;
        }
        if (ra > kPi * 2) {
            ra = ra - kPi * 2;
        }
    }
}

/// MARK: standard_equatorial2

void standard_equatorial2(double ra0, double dec0,
                          double x, double y, double cdelt,
                          double& ra, double& dec) {
    auto sin_dec0 = 0.0, cos_dec0 = 0.0;
    sincos(dec0, sin_dec0, cos_dec0);
    x = x * cdelt / (3600.0 * 180.0 / kPi);
    y = y * cdelt / (3600.0 * 180.0 / kPi);
    
    ra = ra0 + arctan2(-x, cos_dec0 - y * sin_dec0);
    if (ra > kPi * 2) {
        ra = ra - kPi * 2;
    }
    if (ra < 0) {
        ra = ra + kPi * 2;
    }
    dec = arcsin((sin_dec0 + y * cos_dec0) / std::sqrt(1.0 + x * x + y * y));
}

/// MARK: ang_sep / ang_sep_two_positions

void ang_sep(double ra1, double dec1, double ra2, double dec2, double& sep) {
    auto sin_dec1 = 0.0, cos_dec1 = 0.0, sin_dec2 = 0.0, cos_dec2 = 0.0;
    sincos(dec1, sin_dec1, cos_dec1);
    sincos(dec2, sin_dec2, cos_dec2);
    auto cos_sep = std::max(-1.0,
                       std::min(1.0,
                                sin_dec1 * sin_dec2
                                    + cos_dec1 * cos_dec2 * std::cos(ra1 - ra2)));
    sep = arccos(cos_sep);
}

void ang_sep_two_positions(const Header& head, int polynomialIndex,
                           double fitsx1, double fitsy1,
                           double fitsx2, double fitsy2,
                           std::string& seperation, std::string& pa) {
    if (head.cdelt2 != 0) {
        auto ra1 = 0.0, dec1 = 0.0, ra2 = 0.0, dec2 = 0.0, sep = 0.0;
        pixel_to_celestial(head, fitsx1, fitsy1, polynomialIndex, ra1, dec1);
        pixel_to_celestial(head, fitsx2, fitsy2, polynomialIndex, ra2, dec2);
        ang_sep(ra1, dec1, ra2, dec2, sep);
        
        // Convert to degrees
        sep = sep * 180.0 / kPi;
        if (sep < 1.0 / 60.0) {
            seperation = std::to_string(static_cast<long long>(std::lround(sep * 3600.0))) + "\"";
        } else if (sep < 1.0) {
            seperation = floatstr_fixed(sep * 60.0, 2) + "'";
        } else {
            seperation = floatstr_fixed(sep, 2) + "\xC2\xB0";  // UTF-8 degree
        }
        pa = floatstr_fixed(position_angle(ra2, dec2, ra1, dec1) * 180.0 / kPi, 0) + "\xC2\xB0";
    } else {
        seperation = floatstr_fixed(std::sqrt(sqr(fitsx2 - fitsx1) + sqr(fitsy2 - fitsy1)), 2)
                     + " pixels";
        pa = floatstr_fixed(arctan2(fitsy2 - fitsy1, fitsx2 - fitsx1) * 180.0 / kPi, 0)
             + "\xC2\xB0";
    }
}

/// MARK: Text -> radians

void ra_text_to_radians(std::string inp, double& ra, bool& errorRA) {
    auto error1 = 0, error2 = 0, error3 = 0;
    pas_val(inp, ra, error1);
    
    if (error1 != 0) {
        inp = uppercase(inp);
        auto degrees = pos_str("D", inp) > 0;
        inp = string_replace_all(inp, ",", ".");
        
        // Filter to digits, '.', '-', and spaces
        auto data = std::string{};
        for (size_t i = 0; i < inp.size(); ++i) {
            auto c = static_cast<unsigned char>(inp[i]);
            if (((c >= 48) && (c <= 57)) || (inp[i] == '.') || (inp[i] == '-')) {
                data += inp[i];
            } else {
                data += ' ';
            }
        }
        
        // Collapse double-spaces
        auto idx = 0;
        do {
            idx = pos_str("  ", data);
            if (idx > 0) {
                pas_delete(data, idx, 1);
            }
        } while (idx != 0);
        
        data = trim(data) + ' ';
        auto plusmin = (pos_str("-", data) > 0) ? -1.0 : 1.0;
        
        // Parse hours/minutes/seconds
        auto position1 = pos_str(" ", data);
        auto rah = 0.0, ram = 0.0, ras = 0.0;
        pas_val(pas_copy(data, 1, position1 - 1), rah, error1);
        if (degrees) {
            rah = rah * 24.0 / 360.0;
        }
        
        auto position2 = posex_str(" ", data, position1 + 1);
        if (position2 - position1 > 1) {
            pas_val(pas_copy(data, position1 + 1, position2 - position1 - 1), ram, error2);
            auto position3 = posex_str(" ", data, position2 + 1);
            if (position3 - position2 > 1) {
                pas_val(pas_copy(data, position2 + 1, position3 - position2 - 1), ras, error3);
            } else {
                ras = 0;
                error3 = 0;
            }
        } else {
            ram = 0;
            error2 = 0;
            ras = 0;
            error3 = 0;
        }
        
        ra = plusmin * (std::abs(rah) + ram / 60.0 + ras / 3600.0) * kPi / 12.0;
        errorRA = ((error1 != 0) || (error2 > 1) || (error3 != 0) || (ra > 2 * kPi));
    } else {
        errorRA = false;
        ra = ra * kPi / 12.0;
    }
}

void dec_text_to_radians(std::string inp, double& dec, bool& errorDEC) {
    auto error1 = 0, error2 = 0, error3 = 0;
    pas_val(inp, dec, error1);
    
    if (error1 != 0) {
        inp = string_replace_all(inp, ",", ".");
        
        // Filter to digits, '.', '-', and spaces
        auto data = std::string{};
        for (size_t i = 0; i < inp.size(); ++i) {
            auto c = static_cast<unsigned char>(inp[i]);
            if (((c >= 48) && (c <= 57)) || (inp[i] == '.') || (inp[i] == '-')) {
                data += inp[i];
            } else {
                data += ' ';
            }
        }
        
        // Collapse double-spaces
        auto idx = 0;
        do {
            idx = pos_str("  ", data);
            if (idx > 0) {
                pas_delete(data, idx, 1);
            }
        } while (idx != 0);
        
        data = trim(data) + ' ';
        auto plusmin = (pos_str("-", data) > 0) ? -1 : 1;
        
        // Parse degrees/minutes/seconds
        auto position1 = pos_str(" ", data);
        auto decd = 0.0, decm = 0.0, decs = 0.0;
        pas_val(pas_copy(data, 1, position1 - 1), decd, error1);
        
        auto position2 = posex_str(" ", data, position1 + 1);
        if (position2 - position1 > 1) {
            pas_val(pas_copy(data, position1 + 1, position2 - position1 - 1), decm, error2);
            auto position3 = posex_str(" ", data, position2 + 1);
            if (position3 - position2 > 1) {
                pas_val(pas_copy(data, position2 + 1, position3 - position2 - 1), decs, error3);
            } else {
                decs = 0;
                error3 = 0;
            }
        } else {
            decm = 0;
            error2 = 0;
            decs = 0;
            error3 = 0;
        }
        
        dec = plusmin * (std::abs(decd) + decm / 60.0 + decs / 3600.0) * kPi / 180.0;
        errorDEC = ((error1 != 0) || (error2 > 1) || (error3 != 0));
    } else {
        errorDEC = false;
        dec = dec * kPi / 180.0;
    }
}

/// MARK: decode_string

bool decode_string(std::string data0,
                   const Header& head, int polynomialIndex,
                   double& ra4, double& dec4) {
    auto error1 = false, error2 = false;
    
    data0 = uppercase(data0);
    auto degrees = pos_str("D", data0) > 0;
    data0 = string_replace_first(data0, "S.",   ".");
    data0 = string_replace_first(data0, "\".",  ".");
    data0 = string_replace_first(data0, "R.A.", "");
    data0 = string_replace_first(data0, "DECL.", "");
    
    if (data0 == "c" || data0 == "C") {
        // Pre-uppercase compare against lower-case "c" never matches; both
        // branches kept for exact numerical parity with the original.
        pixel_to_celestial(head,
                           (head.width + 1) / 2.0,
                           (head.height + 1) / 2.0,
                           polynomialIndex, ra4, dec4);
        error1 = false;
        error2 = false;
    } else {
        // Filter to digits, '.', '-', and spaces
        auto data1 = std::string{};
        for (size_t i = 0; i < data0.size(); ++i) {
            auto c = static_cast<unsigned char>(data0[i]);
            if (((c >= 48) && (c <= 57)) || (data0[i] == '.') || (data0[i] == '-')) {
                data1 += data0[i];
            } else {
                data1 += ' ';
            }
        }
        
        // Collapse double-spaces
        auto idx = 0;
        do {
            idx = pos_str("  ", data1);
            if (idx > 0) {
                pas_delete(data1, idx, 1);
            }
        } while (idx != 0);
        
        // Trim leading/trailing spaces
        while (!data1.empty() && data1.front() == ' ') {
            pas_delete(data1, 1, 1);
        }
        while (!data1.empty() && data1.back() == ' ') {
            pas_delete(data1, static_cast<int>(data1.size()), 1);
        }
        
        // Locate space-delimited fields
        auto pos1 = pos_str(" ", data1);
        if (pos1 == 0) {
            return false;
        }
        auto pos2 = posex_str(" ", data1, pos1 + 1);
        if (pos2 == 0) {
            pos2 = static_cast<int>(data1.size()) + 1;
        }
        auto pos3 = posex_str(" ", data1, pos2 + 1);
        if (pos3 == 0) {
            pos3 = static_cast<int>(data1.size()) + 1;
        }
        auto pos4 = posex_str(" ", data1, pos3 + 1);
        if (pos4 == 0) {
            pos4 = static_cast<int>(data1.size()) + 1;
        }
        auto pos5 = posex_str(" ", data1, pos4 + 1);
        if (pos5 == 0) {
            pos5 = static_cast<int>(data1.size()) + 1;
        }
        auto pos6 = posex_str(" ", data1, pos5 + 1);
        if (pos6 == 0) {
            pos6 = static_cast<int>(data1.size()) + 1;
        }
        
        // Split into RA and Dec text
        auto ra_text = std::string{};
        auto dec_text = std::string{};
        if (pos5 != pos6) {
            ra_text  = pas_copy(data1, 1, pos3);
            dec_text = pas_copy(data1, pos3 + 1, 99);
        } else if (pos3 != pos4) {
            ra_text  = pas_copy(data1, 1, pos2);
            dec_text = pas_copy(data1, pos2 + 1, 99);
        } else {
            ra_text = pas_copy(data1, 1, pos1);
            if (degrees) {
                ra_text = std::string("D") + ra_text;
            }
            dec_text = pas_copy(data1, pos1 + 1, 99);
        }
        
        ra_text_to_radians(ra_text, ra4, error1);
        dec_text_to_radians(dec_text, dec4, error2);
    }
    return ((error1 == false) && (error2 == false));
}

/// MARK: position_to_string

std::string position_to_string(std::string_view sep, double ra, double dec) {
    if (coord_frame == 0) {
        return prepare_ra8(ra, ": ") + std::string(sep)
             + prepare_dec2(dec, "\xC2\xB0 ");  // UTF-8 degree
    } else if (coord_frame == 1) {
        return floatstr_fixed(ra * 180.0 / kPi, 8) + "\xC2\xB0"
             + std::string(sep)
             + floatstr_fixed(dec * 180.0 / kPi, 8) + "\xC2\xB0";
    } else if (coord_frame == 2) {
        auto l = 0.0, b = 0.0;
        EQU_GAL(ra, dec, l, b);
        return floatstr_fixed(l * 180.0 / kPi, 3) + std::string(sep)
             + floatstr_fixed(b * 180.0 / kPi, 3) + " \xC2\xB0 gal";
    } else if (coord_frame == 3) {
        auto az = 0.0, alt = 0.0;
        if (calculate_az_alt_basic(ra, dec, az, alt) == false) {
            return {};
        }
        return floatstr_fixed(az * 180.0 / kPi, 1) + "\xC2\xB0"
             + std::string(sep)
             + floatstr_fixed(alt * 180.0 / kPi, 1) + "\xC2\xB0";
    }
    return "Error";
}

/// MARK: prepare_* family

std::string prepare_ra5(double rax, std::string_view sep) {
    rax = rax + kPi * 0.1 / (24.0 * 60.0);
    rax = rax * 12.0 / kPi;
    auto h  = static_cast<int>(std::trunc(rax));
    auto m  = static_cast<int>(std::trunc((rax - h) * 60.0));
    auto dm = static_cast<int>(std::trunc((rax - h - m / 60.0) * 600.0));
    auto b = std::format("{:2}", static_cast<int>(std::trunc(static_cast<double>(h))));
    auto dchar = static_cast<char>(dm + 48);
    return b + std::string(sep) + leading_zero(m) + "." + std::string(1, dchar);
}

std::string prepare_dec4(double decx, std::string_view sep) {
    auto sign = (decx < 0) ? '-' : '+';
    decx = std::abs(decx) + kPi / (360.0 * 60.0);
    decx = decx * 180.0 / kPi;
    auto g = static_cast<int>(std::trunc(decx));
    auto m = static_cast<int>(std::trunc((decx - g) * 60.0));
    auto b = std::format("{}", g);
    return std::string(1, sign) + b + std::string(sep) + leading_zero(m);
}

std::string prepare_ra6(double rax, std::string_view sep) {
    rax = rax + kPi / (24.0 * 60.0 * 60.0);
    rax = rax * 12.0 / kPi;
    auto h = static_cast<int>(std::trunc(rax));
    auto m = static_cast<int>(std::trunc((rax - h) * 60.0));
    auto s = static_cast<int>(std::trunc((rax - h - m / 60.0) * 3600.0));
    return leading_zero(h) + std::string(sep) + leading_zero(m) + " " + leading_zero(s);
}

std::string prepare_ra(double rax, std::string_view sep) {
    rax = rax + kPi * 0.1 / (24.0 * 60.0 * 60.0);
    rax = rax * 12.0 / kPi;
    auto h = static_cast<int>(std::trunc(rax));
    auto m = static_cast<int>(std::trunc((rax - h) * 60.0));
    auto s = static_cast<int>(std::trunc((rax - h - m / 60.0) * 3600.0));
    auto ds = static_cast<int>(std::trunc((rax - h - m / 60.0 - s / 3600.0) * 36000.0));
    auto dchar = static_cast<char>(ds + 48);
    return leading_zero(h) + std::string(sep) + leading_zero(m) + " "
         + leading_zero(s) + "." + std::string(1, dchar);
}

std::string prepare_dec(double decx, std::string_view sep) {
    auto sign = (decx < 0) ? '-' : '+';
    decx = std::abs(decx) + kPi / (360.0 * 60.0 * 60.0);
    decx = decx * 180.0 / kPi;
    auto g = static_cast<int>(std::trunc(decx));
    auto m = static_cast<int>(std::trunc((decx - g) * 60.0));
    auto s = static_cast<int>(std::trunc((decx - g - m / 60.0) * 3600.0));
    auto pad = std::min<size_t>(sep.size(), 2);
    return std::string(1, sign) + leading_zero(g) + std::string(sep) + leading_zero(m)
         + spaces(pad) + leading_zero(s);
}

std::string prepare_ra8(double rax, std::string_view sep) {
    rax = rax + kPi * 0.01 / (24.0 * 60.0 * 60.0);
    rax = rax * 12.0 / kPi;
    auto h = static_cast<int>(std::trunc(rax));
    auto m = static_cast<int>(std::trunc((rax - h) * 60.0));
    auto s = static_cast<int>(std::trunc((rax - h - m / 60.0) * 3600.0));
    auto ds = static_cast<int>(std::trunc((rax - h - m / 60.0 - s / 3600.0) * 360000.0));
    auto b = std::format("{:2}", h);
    auto pad = std::min<size_t>(sep.size(), 2);
    return b + std::string(sep) + leading_zero(m) + spaces(pad)
         + leading_zero(s) + "." + leading_zero(ds);
}

std::string prepare_dec2(double decx, std::string_view sep) {
    auto sign = (decx < 0) ? '-' : '+';
    decx = std::abs(decx) + kPi * 0.1 / (360.0 * 60.0 * 60.0);
    decx = decx * 180.0 / kPi;
    auto g = static_cast<int>(std::trunc(decx));
    auto m = static_cast<int>(std::trunc((decx - g) * 60.0));
    auto s = static_cast<int>(std::trunc((decx - g - m / 60.0) * 3600.0));
    auto ds = static_cast<int>(std::trunc((decx - g - m / 60.0 - s / 3600.0) * 36000.0));
    auto b = std::format("{:2}", g);
    auto ds2 = std::format("{:1}", ds);
    auto pad = std::min<size_t>(sep.size(), 2);
    return std::string(1, sign) + b + std::string(sep) + leading_zero(m)
         + spaces(pad) + leading_zero(s) + "." + ds2;
}

std::string prepare_IAU_designation(double rax, double decx) {
    rax = rax + kPi * 2.0 * 0.05 / (24.0 * 60.0 * 60.0);
    rax = rax * 12.0 / kPi;
    auto hh = static_cast<int>(std::trunc(rax));
    auto mm = static_cast<int>(std::trunc((rax - hh) * 60.0));
    auto ss = static_cast<int>(std::trunc((rax - hh - mm / 60.0) * 3600.0));
    auto ds = static_cast<int>(std::trunc((rax - hh - mm / 60.0 - ss / 3600.0) * 36000.0));
    
    auto sign = (decx < 0) ? '-' : '+';
    decx = std::abs(decx) + kPi * 2.0 * 0.5 / (360.0 * 60.0 * 60.0);
    decx = decx * 180.0 / kPi;
    auto g = static_cast<int>(std::trunc(decx));
    auto m = static_cast<int>(std::trunc((decx - g) * 60.0));
    auto s = static_cast<int>(std::trunc((decx - g - m / 60.0) * 3600.0));
    
    auto dchar = static_cast<char>(ds + 48);
    return leading_zero(hh) + leading_zero(mm) + leading_zero(ss) + "."
         + std::string(1, dchar) + std::string(1, sign)
         + leading_zero(g) + leading_zero(m) + leading_zero(s);
}

/// MARK: Jd_To_MPCDate

std::string Jd_To_MPCDate(double jd) {
    jd = jd + (0.5 / (24.0 * 3600.0));
    
    auto Z = std::trunc(jd + 0.5);
    auto F = frac(jd + 0.5);
    auto A = 0.0;
    
    if (Z < 2299160.5) {
        A = Z;
    } else {
        auto g = pas_int((Z - 1867216.25) / 36524.25);
        A = Z + 1 + g - std::trunc(g / 4.0);
    }
    auto B = A + 1524;
    auto C = std::trunc((B - 122.1) / 365.25);
    auto D = std::trunc(365.25 * C);
    auto E = std::trunc((B - D) / 30.6001);
    auto T = B - D - pas_int(30.6001 * E) + F;  // day of the month
    auto M = (E < 14) ? (E - 1) : (E - 13);
    auto J = (M > 2) ? (C - 4716) : (C - 4715);
    
    // Year: width 4, space pad
    auto year3 = std::format("{:4}", static_cast<long long>(std::trunc(J)));
    
    // Day: width 8, 5 decimals
    auto day_val = std::trunc(T) + F;
    auto day = std::format("{:8.5f}", day_val);
    if (!day.empty() && day[0] == ' ') {
        day[0] = '0';
    }
    
    return year3 + " " + leading_zero(static_cast<long long>(std::trunc(M))) + " " + day;
}

/// MARK: write_astronomy_wcs

void write_astronomy_wcs(const std::filesystem::path& filen,
                         const std::vector<std::string>& headerLines) {
    auto out = std::ofstream(filen, std::ios::out | std::ios::trunc);
    if (!out) {
        return;
    }
    
    try {
        const auto empthy_line = std::string(80, ' ');
        auto i = size_t{0};
        
        while (true) {
            if (i < headerLines.size()) {
                auto line0 = headerLines[i];
                while (line0.size() < 80) {
                    line0 += ' ';
                }
                if (line0.size() > 80) {
                    line0.resize(80);
                }
                out.write(line0.data(), 80);
            } else {
                out.write(empthy_line.data(), 80);
            }
            ++i;
            
            // Pad to 2880-byte FITS block boundary
            if (i >= headerLines.size() && frac(static_cast<double>(i) * 80.0 / 2880.0) == 0.0) {
                break;
            }
        }
    } catch (...) {
        // Swallow exceptions to match original behavior.
    }
}
    
} // namespace
