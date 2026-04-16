///----------------------------------------
///      @file online_gaia.cpp
///   @ingroup ASTAP++
///     @brief Online Gaia catalog query, parsing, and photometric passband conversion.
///    @author Ported from Han Kleijn's ASTAP (unit_online_gaia.pas). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "online_gaia.h"

#include <charconv>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>
#include <string_view>
#include <vector>

///----------------------------------------
namespace astap::reference {
///----------------------------------------

/// MARK: Module-level globals

std::array<std::vector<double>, 6> online_database{};
double gaia_ra = 0.0;
double gaia_dec = 0.0;
std::string passband_active{};

///----------------------------------------
/// MARK: File-local helpers
///----------------------------------------

namespace {

constexpr auto kPi = std::numbers::pi;

/// @brief Square of a double value.
[[nodiscard]] constexpr double sqr(double x) noexcept { return x * x; }

///----------------------------------------
/// @brief Split text on newlines, stripping trailing carriage returns.
/// @param text The input string to split.
/// @return A vector of individual lines.
///----------------------------------------

[[nodiscard]] std::vector<std::string> split_lines(const std::string& text) {
    auto lines = std::vector<std::string>{};
    auto start = std::size_t{0};
    for (auto i = std::size_t{0}; i <= text.size(); ++i) {
        if (i == text.size() || text[i] == '\n') {
            auto end = i;
            if (end > start && text[end - 1] == '\r') {
                --end;
            }
            lines.emplace_back(text.data() + start, end - start);
            start = i + 1;
        }
    }
    return lines;
}

///----------------------------------------
/// @brief Find the next ASCII space in @p s at or after @p from.
/// @param s    String view to search.
/// @param from Starting offset (0-based).
/// @return Position of the space, or @c std::string::npos if none.
///----------------------------------------

[[nodiscard]] std::size_t find_space(std::string_view s, std::size_t from) noexcept {
    return s.find(' ', from);
}

///----------------------------------------
/// @brief Parse a double from a string view, returning 0.0 on failure.
/// @param s Input string view (leading/trailing whitespace is trimmed).
/// @return Parsed value, or 0.0 on parse error or empty input.
///----------------------------------------

[[nodiscard]] double parse_double(std::string_view s) noexcept {
    // Trim leading whitespace
    while (!s.empty() && (s.front() == ' ' || s.front() == '\t')) {
        s.remove_prefix(1);
    }
    
    // Trim trailing whitespace
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r')) {
        s.remove_suffix(1);
    }
    
    if (s.empty()) {
        return 0.0;
    }
    
    auto out = 0.0;
    auto [p, ec] = std::from_chars(s.data(), s.data() + s.size(), out);
    if (ec != std::errc{}) {
        return 0.0;
    }
    return out;
}

///----------------------------------------
/// @brief Angular separation between two sky positions (radians).
/// @details Uses the Vincenty formula form, adequate for small separations.
/// @param ra1  Right ascension of the first position (radians).
/// @param dec1 Declination of the first position (radians).
/// @param ra2  Right ascension of the second position (radians).
/// @param dec2 Declination of the second position (radians).
/// @param[out] sep Angular separation (radians).
///----------------------------------------

void ang_sep(double ra1, double dec1, double ra2, double dec2, double& sep) noexcept {
    auto const dra = ra1 - ra2;
    auto const sd1 = std::sin(dec1);
    auto const sd2 = std::sin(dec2);
    auto const cd1 = std::cos(dec1);
    auto const cd2 = std::cos(dec2);
    auto const x   = sd1 * sd2 + cd1 * cd2 * std::cos(dra);
    
    // Clamp to [-1, 1] to guard against FP drift
    auto const c = std::max(-1.0, std::min(1.0, x));
    sep = std::acos(c);
}

///----------------------------------------
/// @brief Format a double as fixed-point with 2 fractional digits.
/// @param v Value to format.
/// @return Formatted string.
///----------------------------------------

[[nodiscard]] std::string fixed2(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.2f", v);
    return std::string(buf);
}

///----------------------------------------
/// @brief Format a double as fixed-point with 6 fractional digits.
/// @param v Value to format.
/// @return Formatted string.
///----------------------------------------

[[nodiscard]] std::string fixed6(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.6f", v);
    return std::string(buf);
}

///----------------------------------------
/// @brief Format a double as fixed-point with 10 fractional digits.
/// @param v Value to format.
/// @return Formatted string.
///----------------------------------------

[[nodiscard]] std::string fixed10(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.10f", v);
    return std::string(buf);
}

///----------------------------------------
/// @brief Parse the VizieR text response and fill @ref online_database.
/// @param slist Lines of the VizieR response.
///----------------------------------------

void extract_stars(const std::vector<std::string>& slist) {
    // Reset and size the six columns to an upper bound; shrink at the end
    auto const n = slist.size();
    for (auto& col : online_database) {
        col.assign(n, 0.0);
    }
    
    auto datalines = false;
    auto count2 = std::size_t{0};
    
    // Skip the VizieR header lines and iterate the body
    auto count = std::size_t{35};
    while (n >= 2 && count < n - 2) {
        auto const& regel = slist[count];
        ++count;
        
        if (datalines) {
            // Find four space-separated column boundaries
            auto const p1 = regel.find(' ');
            if (p1 == std::string::npos) {
                continue;
            }
            auto const p2 = find_space(regel, p1 + 3);
            if (p2 == std::string::npos) {
                continue;
            }
            auto const p3 = find_space(regel, p2 + 3);
            if (p3 == std::string::npos) {
                continue;
            }
            auto const p4 = find_space(regel, p3 + 3);
            
            // Guard: first character must be a digit
            if (regel.empty()) {
                continue;
            }
            auto const c0 = static_cast<unsigned char>(regel[0]);
            if (!(c0 >= '0' && c0 <= '9')) {
                continue;
            }
            
            auto const sv = std::string_view{regel};
            
            // RA (degrees -> radians)
            auto const ra = parse_double(sv.substr(0, p1));
            online_database[0][count2] = ra * kPi / 180.0;
            
            // Dec (degrees -> radians)
            auto const dec = parse_double(sv.substr(p1 + 1, p2 - p1 - 1));
            online_database[1][count2] = dec * kPi / 180.0;
            
            // G magnitude
            auto const g = parse_double(sv.substr(p2 + 1, p3 - p2 - 1));
            online_database[2][count2] = g;
            
            // BP magnitude
            auto const bp_end = (p4 == std::string::npos) ? sv.size() : p4;
            auto const bp = parse_double(sv.substr(p3 + 1, bp_end - p3 - 1));
            online_database[3][count2] = bp;
            
            // RP magnitude
            auto rp = 0.0;
            if (p4 != std::string::npos && p4 + 1 <= sv.size()) {
                rp = parse_double(sv.substr(p4 + 1));
            }
            online_database[4][count2] = rp;
            
            // Default the "active" column to BP until convert_magnitudes overwrites it
            online_database[5][count2] = bp;
            ++count2;
        } else if (regel.compare(0, 7, "RA_ICRS") == 0) {
            datalines = true;
            ++count; // skip the dashes line
        }
    }
    
    // Shrink columns to actual star count
    for (auto& col : online_database) {
        col.resize(count2);
    }
}
    
} // namespace

///----------------------------------------
/// MARK: transform_gaia
///----------------------------------------

[[nodiscard]] double transform_gaia(const std::string& filter,
                                    double magG, double magBP, double magRP) {
    // BP is a direct pass-through
    if (filter == "BP") {
        return magBP;
    }
    
    // Assume failure
    auto result = 0.0;
    
    if (magG != 0.0 && magBP != 0.0 && magRP != 0.0) {
        // Straylight quality check via flux ratio
        auto const Gflux  = std::pow(10.0, (22.0 - magG)  / 2.5);
        auto const BPflux = std::pow(10.0, (22.0 - magBP) / 2.5);
        auto const RPflux = std::pow(10.0, (22.0 - magRP) / 2.5);
        auto const c = (BPflux + RPflux) / Gflux;
        auto const straylight = (c > 4.0) && (magG > magBP);
        
        if (!straylight) {
            auto const BminR = magBP - magRP;
            
            // Coefficients taken verbatim from the Gaia EDR3/DR3 calibration;
            // do not edit without cross-checking upstream.
            if (filter == "V") {
                // Johnson-Cousins V
                if (BminR >= -0.5 && BminR <= 5.0) {
                    result = magG + 0.02704 - 0.01424 * BminR
                           + 0.2156 * sqr(BminR) - 0.01426 * sqr(BminR) * BminR;
                }
            } else if (filter == "R") {
                // Johnson-Cousins R
                if (BminR > 0.0 && BminR < 4.0) {
                    result = magG + 0.02275 - 0.3961 * BminR
                           + 0.1243 * sqr(BminR) + 0.01396 * sqr(BminR) * BminR
                           - 0.003775 * sqr(sqr(BminR));
                }
            } else if (filter == "B") {
                // Johnson-B via Tycho Bt/Vt
                if (BminR > -0.3 && BminR < 3.0) {
                    auto const Vt = magG + 0.01077 + 0.0682 * BminR
                                  + 0.2387 * sqr(BminR) - 0.02342 * sqr(BminR) * BminR;
                    auto const Bt = magG + 0.004288 + 0.8547 * BminR
                                  - 0.1244 * sqr(BminR) + 0.9085 * sqr(BminR) * BminR
                                  - 0.4843 * sqr(sqr(BminR))
                                  + 0.06814 * sqr(sqr(BminR)) * BminR;
                    auto const V  = magG + 0.02704 - 0.01424 * BminR
                                  + 0.2156 * sqr(BminR) - 0.01426 * sqr(BminR) * BminR;
                    result = V + 0.850 * (Bt - Vt);
                }
            } else if (filter == "SR") {
                // SDSS-r
                if (BminR > 0.0 && BminR < 3.0) {
                    result = magG + 0.09837 - 0.08592 * BminR
                           - 0.1907 * sqr(BminR) + 0.1701 * sqr(BminR) * BminR
                           - 0.02263 * sqr(sqr(BminR));
                }
            } else if (filter == "SI") {
                // SDSS-i
                if (BminR > 0.5 && BminR < 2.0) {
                    result = magG + 0.293 - 0.6404 * BminR
                           + 0.09609 * sqr(BminR) + 0.002104 * sqr(BminR) * BminR;
                }
            } else if (filter == "SG") {
                // SDSS-g
                if (BminR > 0.3 && BminR < 3.0) {
                    result = magG - 0.2199 + 0.6365 * BminR
                           + 0.1548 * sqr(BminR) - 0.0064 * sqr(BminR) * BminR;
                }
            }
        }
    }
    
    return result;
}

///----------------------------------------
/// MARK: convert_magnitudes
///----------------------------------------

void convert_magnitudes(const std::string& passband) {
    // Already active -- nothing to do
    if (passband == passband_active) {
        return;
    }
    
    // Transform every star's active-column magnitude
    auto const n = online_database[0].size();
    for (auto i = std::size_t{0}; i < n; ++i) {
        online_database[5][i] = transform_gaia(passband,
                                               online_database[2][i],
                                               online_database[3][i],
                                               online_database[4][i]);
    }
    passband_active = passband;
}

///----------------------------------------
/// MARK: report_one_star_magnitudes
///----------------------------------------

void report_one_star_magnitudes(double ra, double dec,
                                double& b, double& v, double& r,
                                double& sg, double& sr, double& si,
                                double& g, double& bp, double& rp) {
    // Zero all outputs
    b = 0; v = 0; r = 0; sg = 0; sr = 0; si = 0; g = 0; bp = 0; rp = 0;
    
    auto const n = online_database[0].size();
    if (n == 0) {
        return;
    }
    
    // 5 arcsec threshold in radians
    constexpr auto threshold = 5.0 * std::numbers::pi / (180.0 * 60.0 * 60.0);
    auto sep = 0.0;
    
    for (auto i = std::size_t{0}; i < n; ++i) {
        ang_sep(ra, dec, online_database[0][i], online_database[1][i], sep);
        if (sep < threshold) {
            auto const gM  = online_database[2][i];
            auto const bpM = online_database[3][i];
            auto const rpM = online_database[4][i];
            
            // Johnson-Cousins
            b  = transform_gaia("B",  gM, bpM, rpM);
            v  = transform_gaia("V",  gM, bpM, rpM);
            r  = transform_gaia("R",  gM, bpM, rpM);
            
            // SDSS
            sg = transform_gaia("SG", gM, bpM, rpM);
            sr = transform_gaia("SR", gM, bpM, rpM);
            si = transform_gaia("SI", gM, bpM, rpM);
            
            // Raw Gaia
            g  = gM;
            bp = bpM;
            rp = rpM;
            break;
        }
    }
}

///----------------------------------------
/// MARK: read_stars_online
///----------------------------------------

[[nodiscard]] bool read_stars_online(IHttpClient& http,
                                     double telescope_ra,
                                     double telescope_dec,
                                     double search_field,
                                     double magli) {
    // Format pointing and field for the URL
    auto const ra8  = fixed10(std::abs(telescope_ra  * 180.0 / kPi));
    auto const dec8 = fixed10(std::abs(telescope_dec * 180.0 / kPi));
    auto const sgn  = std::string{(telescope_dec >= 0.0) ? "%2B" : "%2D"};
    
    auto const field       = fixed6(search_field * 3600.0 * 180.0 / kPi);
    auto const window_size = std::string{"&-c.bs="} + field + "/" + field;
    auto const mag_lim     = fixed2(magli);
    
    // Build the VizieR query URL
    auto const url = std::string{
        "http://vizier.u-strasbg.fr/viz-bin/asu-txt?-source=I/355/Gaiadr3"
        "&-out=RA_ICRS,DE_ICRS,Gmag,BPmag,RPmag&-c="}
        + ra8 + sgn + dec8 + window_size
        + "&-out.max=200000&Gmag=<" + mag_lim;
        
    // Perform the HTTP request
    auto resp = http.get(url);
    if (!resp.has_value()) {
        return false;
    }
    
    // Parse the response
    auto const slist = split_lines(*resp);
    if (slist.size() <= 31) {
        return false;
    }
    
    // Update cache and extract stars
    gaia_ra  = telescope_ra;
    gaia_dec = telescope_dec;
    extract_stars(slist);
    return true;
}
    
} // namespace
