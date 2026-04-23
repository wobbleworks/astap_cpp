///----------------------------------------
///      @file online.cpp
///   @ingroup ASTAP++
///     @brief Online catalog query helpers for AAVSO VSP/VSX, Simbad, and Vizier.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "online.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

///----------------------------------------
namespace astap::core {
///----------------------------------------

/// MARK: - Globals

// Module-level cached arrays.
std::vector<Auid> vsp;
std::vector<VarStar> vsx;

namespace {

constexpr auto kPi = std::numbers::pi;

/// MARK: - String Helpers

// 1-based position search (0 = not found), preserving the original
// offset semantics for faithful porting.
[[nodiscard]] std::size_t posex(std::string_view needle, std::string_view hay, std::size_t offset_1based) {
    if (offset_1based == 0) {
        offset_1based = 1;
    }
    if (offset_1based > hay.size()) {
        return 0;
    }
    auto p = hay.find(needle, offset_1based - 1);
    if (p == std::string_view::npos) {
        return 0;
    }
    return p + 1;
}

[[nodiscard]] std::size_t pos1(std::string_view needle, std::string_view hay) {
    auto p = hay.find(needle);
    if (p == std::string_view::npos) {
        return 0;
    }
    return p + 1;
}

// 1-based substring extraction, clamped.
[[nodiscard]] std::string copy1(std::string_view s, std::size_t start_1based, std::size_t len) {
    if (start_1based == 0 || start_1based > s.size()) {
        return {};
    }
    auto i = start_1based - 1;
    return std::string(s.substr(i, std::min(len, s.size() - i)));
}

// 1-based char access; returns '\0' on out-of-bounds.
[[nodiscard]] char at1(std::string_view s, std::size_t i_1based) noexcept {
    if (i_1based == 0 || i_1based > s.size()) {
        return '\0';
    }
    return s[i_1based - 1];
}

[[nodiscard]] std::string trim(std::string_view s) {
    auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string_view::npos) {
        return {};
    }
    auto b = s.find_last_not_of(" \t\r\n");
    return std::string(s.substr(a, b - a + 1));
}

[[nodiscard]] std::string string_replace_all(std::string_view s, char from, char to) {
    auto out = std::string(s);
    std::ranges::replace(out, from, to);
    return out;
}

// Permissive double parse (leading/trailing spaces). Returns true on
// success and writes @p out. On failure leaves @p out alone.
[[nodiscard]] bool val_double(std::string_view s, double& out) {
    auto t = trim(s);
    if (t.empty()) {
        return false;
    }
    const auto* first = t.data();
    const auto* last = first + t.size();
    auto v = 0.0;
    auto [ptr, ec] = std::from_chars(first, last, v);
    if (ec != std::errc{}) {
        return false;
    }
    out = v;
    return true;
}

// Lenient float parser; returns 0 on failure.
[[nodiscard]] double strtofloat1(std::string_view s) {
    auto v = 0.0;
    (void)val_double(s, v);
    return v;
}

// Alternative lenient parser permitting comma decimal mark.
[[nodiscard]] double strtofloat2(std::string_view s) {
    auto t = std::string(s);
    std::ranges::replace(t, ',', '.');
    return strtofloat1(t);
}

// Fixed-precision float-to-string.
[[nodiscard]] std::string fixed(double v, int digits) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", digits, v);
    return std::string(buf);
}

[[nodiscard]] std::string floattostr6(double v) { return fixed(v, 6); }
[[nodiscard]] std::string floattostr4(double v) { return fixed(v, 4); }

/// MARK: - Coordinate Parsing

// Parse "HH MM SS.s" or "HH:MM:SS.s" into radians.
void ra_text_to_radians(std::string_view txt, double& ra_rad, bool& err) {
    err = true;
    auto t = std::string(txt);
    for (auto& c : t) {
        if (c == ':') {
            c = ' ';
        }
    }
    auto is = std::istringstream(t);
    auto h = 0.0;
    auto m = 0.0;
    auto s = 0.0;
    if (!(is >> h >> m >> s)) {
        return;
    }
    ra_rad = (h + m / 60.0 + s / 3600.0) * 15.0 * kPi / 180.0;
    err = false;
}

void dec_text_to_radians(std::string_view txt, double& dec_rad, bool& err) {
    err = true;
    auto t = std::string(txt);
    for (auto& c : t) {
        if (c == ':') {
            c = ' ';
        }
    }
    auto a = t.find_first_not_of(" \t");
    auto negative = (a != std::string::npos && t[a] == '-');
    if (a != std::string::npos && (t[a] == '+' || t[a] == '-')) {
        t[a] = ' ';
    }
    auto is = std::istringstream(t);
    auto d = 0.0;
    auto m = 0.0;
    auto s = 0.0;
    if (!(is >> d >> m >> s)) {
        return;
    }
    auto v = (std::abs(d) + m / 60.0 + s / 3600.0) * kPi / 180.0;
    dec_rad = negative ? -v : v;
    err = false;
}

// Great-circle separation in radians (haversine).
[[nodiscard]] double ang_sep(double ra1, double dec1, double ra2, double dec2) noexcept {
    auto dra = ra1 - ra2;
    auto cs = std::sin(dec1) * std::sin(dec2) +
              std::cos(dec1) * std::cos(dec2) * std::cos(dra);
    cs = std::clamp(cs, -1.0, 1.0);
    return std::acos(cs);
}

/// MARK: - Line Splitting

// Split text at line breaks.
[[nodiscard]] std::vector<std::string> split_lines(std::string_view text) {
    std::vector<std::string> out;
    std::size_t i = 0;
    while (i <= text.size()) {
        auto j = i;
        while (j < text.size() && text[j] != '\n' && text[j] != '\r') {
            ++j;
        }
        out.emplace_back(text.substr(i, j - i));
        if (j == text.size()) {
            break;
        }
        if (text[j] == '\r' && j + 1 < text.size() && text[j + 1] == '\n') {
            ++j;
        }
        i = j + 1;
    }
    return out;
}

// Pixel -> celestial stub. The full implementation lives elsewhere;
// TODO: replace with the real pixel_to_celestial once ported.
void pixel_to_celestial_stub([[maybe_unused]] const Header& head,
                             [[maybe_unused]] double x,
                             [[maybe_unused]] double y,
                             [[maybe_unused]] int formalism,
                             [[maybe_unused]] double& ra,
                             [[maybe_unused]] double& dec) {
    // Intentionally empty: leaves ra/dec untouched.
}
 
} // namespace

/// MARK: - download_vsp

bool download_vsp(IHttpClient& http, const Header& head, double limiting_mag) {
    // Compute field-of-view in arcminutes
    auto fov = static_cast<int>(std::round(
        std::sqrt(static_cast<double>(head.width) * head.width +
                  static_cast<double>(head.height) * head.height) *
        std::abs(head.cdelt2 * 60.0)));
        
    if (fov > 180) {
        // Clamp limiting magnitude for wide fields
        limiting_mag = std::min(limiting_mag, 12.0);
    }
    
    // Build the AAVSO VSP API URL
    auto url =
        "https://www.aavso.org/apps/vsp/api/chart/?format=json&ra=" +
        floattostr6(head.ra0 * 180.0 / kPi) +
        "&dec=" + floattostr6(head.dec0 * 180.0 / kPi) +
        "&fov=" + std::to_string(fov) +
        "&maglimit=" + floattostr4(limiting_mag);
        
    // Fetch the chart
    auto resp = http.get(url);
    if (!resp) {
        return false;
    }
    const auto& s = *resp;
    if (s.empty()) {
        return false;
    }
    if (s.size() < 256) {
        return false;  // no data for this field
    }
    
    vsp.clear();
    vsp.reserve(2000);
    
    auto count = std::size_t{0};
    auto j = std::size_t{150};  // skip header (original magic offset)
    
    while (count < 10000) {  // practical cap
        // Find the next auid field
        auto i = posex("\"auid\":\"", s, j);
        if (i == 0) {
            break;
        }
        i += std::string_view("\"auid\":\"").size();
        j = posex("\"", s, i);
        if (j == 0) {
            break;
        }
        
        auto entry = Auid{};
        entry.auid = copy1(s, i, j - i);
        
        // Parse RA
        i = posex("\"ra\":\"", s, j);
        if (i == 0) {
            break;
        }
        i += std::string_view("\"ra\":\"").size();
        j = posex("\"", s, i);
        if (j == 0) {
            break;
        }
        auto errRA = false;
        ra_text_to_radians(copy1(s, i, j - i), entry.ra, errRA);
        
        // Parse DEC
        i = posex("\"dec\":\"", s, j);
        if (i == 0) {
            break;
        }
        i += std::string_view("\"dec\":\"").size();
        j = posex("\"", s, i);
        if (j == 0) {
            break;
        }
        auto errDEC = false;
        dec_text_to_radians(copy1(s, i, j - i), entry.dec, errDEC);
        
        // Default magnitudes to "?" (unknown)
        entry.Bmag = "?";  entry.Berr = "";
        entry.Vmag = "?";  entry.Verr = "";
        entry.Rmag = "?";  entry.Rerr = "";
        entry.SGmag = "?"; entry.SGerr = "";
        entry.SRmag = "?"; entry.SRerr = "";
        entry.SImag = "?"; entry.SIerr = "";
        
        auto k = std::size_t{0};
        auto val_c = '\0';
        
        // Parse optional band magnitudes
        while (true) {
            val_c = at1(s, j);
            ++j;
            auto val2 = at1(s, j);
            
            auto read_band = [&](std::string& mag, std::string& err) {
                auto ii = posex("\"mag\":", s, j);
                if (ii == 0) {
                    return;
                }
                ii += std::string_view("\"mag\":").size();
                k = posex(",", s, ii);
                if (k == 0) {
                    return;
                }
                mag = copy1(s, ii, k - ii);
                
                ii = posex("error\":", s, k);
                if (ii == 0) {
                    return;
                }
                ii += std::string_view("error\":").size();
                k = posex("}", s, ii);
                if (k == 0) {
                    return;
                }
                err = copy1(s, ii, k - ii);
            };
            
            if (val_c == '"' && val2 == 'B') {
                read_band(entry.Bmag, entry.Berr);
            } else if (val_c == '"' && val2 == 'V') {
                read_band(entry.Vmag, entry.Verr);
            } else if (val_c == '"' && val2 == 'R') {
                read_band(entry.Rmag, entry.Rerr);
            }
            
            // Non-mutually-exclusive branches for Sloan bands
            if (val_c == 'S' && val2 == 'G') {
                read_band(entry.SGmag, entry.SGerr);
            }
            if (val_c == 'S' && val2 == 'R') {
                read_band(entry.SRmag, entry.SRerr);
            }
            if (val_c == 'S' && val2 == 'I') {
                read_band(entry.SImag, entry.SIerr);
            }
            
            // Recover when k moved forward
            if (j < k) {
                j = k;
            }
            
            if (val_c == ']' || j >= s.size()) {
                break;
            }
        }
        
        vsp.push_back(std::move(entry));
        ++count;
    }
    
    return true;
}

/// MARK: - download_vsx

bool download_vsx(IHttpClient& http,
                  const Header& head,
                  double limiting_mag,
                  double years_since_2000,
                  bool period_filter) {
    // Compute search radius in degrees
    auto radius = std::sqrt(
        static_cast<double>(head.width) * head.width +
        static_cast<double>(head.height) * head.height) *
        std::abs(head.cdelt2 / 2.0);
        
    // AAVSO requirement for wide fields
    if (radius > 3.0) {
        limiting_mag = std::min(12.0, limiting_mag);
    }
    
    // Build the AAVSO VSX API URL
    auto url =
        "https://www.aavso.org/vsx/index.php?view=api.list&ra=" +
        floattostr6(head.ra0 * 180.0 / kPi) +
        "&dec=" + floattostr6(head.dec0 * 180.0 / kPi) +
        "&radius=" + floattostr6(radius) +
        "&tomag=" + floattostr4(limiting_mag) +
        "&format=json";
        
    // Fetch the list
    auto resp = http.get(url);
    if (!resp) {
        return false;
    }
    const auto& s = *resp;
    if (s.empty()) {
        return false;
    }
    if (s.size() < 25) {
        return false;  // no stars in this field
    }
    
    vsx.clear();
    vsx.reserve(2000);
    
    auto count = std::size_t{0};
    auto j = std::size_t{25};
    
    while (count < 10000) {
        // Find the next Name field
        auto i = posex("\"Name\":\"", s, j);
        if (i == 0) {
            break;
        }
        i += std::string_view("\"Name\":\"").size();
        j = posex("\"", s, i);
        if (j == 0) {
            break;
        }
        
        auto entry = VarStar{};
        entry.name = string_replace_all(copy1(s, i, j - i), ' ', '_');
        
        // Parse RA2000
        i = posex("\"RA2000\":\"", s, j);
        if (i == 0) {
            break;
        }
        i += std::string_view("\"RA2000\":\"").size();
        j = posex("\"", s, i);
        if (j == 0) {
            break;
        }
        auto ra_deg = 0.0;
        (void)val_double(copy1(s, i, j - i), ra_deg);
        entry.ra = ra_deg * kPi / 180.0;
        
        // Parse Declination2000
        i = posex("\"Declination2000\":\"", s, j);
        if (i == 0) {
            break;
        }
        i += std::string_view("\"Declination2000\":\"").size();
        j = posex("\"", s, i);
        if (j == 0) {
            break;
        }
        auto dec_deg = 0.0;
        (void)val_double(copy1(s, i, j - i), dec_deg);
        entry.dec = dec_deg * kPi / 180.0;
        
        // Default fields
        entry.maxmag = "?";
        entry.minmag = "?";
        entry.period = "?";
        entry.category = "?";
        
        auto k = std::size_t{0};
        
        // Parse optional fields (MaxMag, MinMag, Category, Period, Epoch, proper motion)
        while (true) {
            ++j;
            auto c0 = at1(s, j);
            auto c1 = at1(s, j + 1);
            auto c2 = at1(s, j + 2);
            auto c3 = at1(s, j + 3);
            auto c12 = at1(s, j + 12);
            auto c13 = at1(s, j + 13);
            
            if (c0 == 'M' && c1 == 'a' && c2 == 'x') {
                auto ii = j + std::string_view("\"MaxMag:\"").size();
                k = posex("\"", s, ii);
                if (k != 0) {
                    entry.maxmag = copy1(s, ii, k - ii);
                }
            } else if (c0 == 'M' && c1 == 'i' && c2 == 'n') {
                auto ii = j + std::string_view("\"MinMag:\"").size();
                k = posex("\"", s, ii);
                if (k != 0) {
                    entry.minmag = copy1(s, ii, k - ii);
                }
            } else if (c0 == 'C' && c1 == 'a' && c2 == 't') {
                auto ii = j + std::string_view("\"Category:\"").size();
                k = posex("\"", s, ii);
                // Read only 3 chars (verbatim from original)
                entry.category = copy1(s, ii, 3);
            } else if (c0 == 'P' && c1 == 'e' && c2 == 'r') {
                auto ii = j + std::string_view("\"Period:\"").size();
                k = posex("\"", s, ii);
                if (k != 0) {
                    entry.period = copy1(s, ii, k - ii);
                }
            } else if (c0 == 'E' && c1 == 'p' && c2 == 'o') {
                auto ii = j + std::string_view("\"Epoch:\"").size();
                k = posex("\"", s, ii);
                if (k != 0) {
                    entry.epoch = copy1(s, ii, k - ii);
                }
            } else if (c0 == 'P' && c1 == 'r' && c2 == 'o' && c3 == 'p' &&
                       c12 == 'D' && c13 == 'e') {
                // ProperMotionDec
                auto ii = j + std::string_view("\"ProperMotionDec:\"").size();
                k = posex("\"", s, ii);
                if (k != 0) {
                    auto pmDec = 0.0;
                    if (val_double(copy1(s, ii, k - ii), pmDec)) {
                        entry.dec += pmDec * years_since_2000 /
                                     ((1000.0 * 3600.0) * 180.0 / kPi);
                    }
                }
            } else if (c0 == 'P' && c1 == 'r' && c2 == 'o' && c3 == 'p' &&
                       c12 == 'R' && c13 == 'A') {
                // ProperMotionRA
                auto ii = j + std::string_view("\"ProperMotionRA:\"").size();
                k = posex("\"", s, ii);
                if (k != 0) {
                    auto pmRA = 0.0;
                    if (val_double(copy1(s, ii, k - ii), pmRA)) {
                        entry.ra += pmRA * years_since_2000 /
                                    (std::cos(entry.dec) *
                                     (1000.0 * 3600.0) * 180.0 / kPi);
                    }
                }
            }
            
            if (j < k) {
                j = k;
            }
            if (at1(s, j) == '}' || j == s.size() - 1) {
                break;
            }
        }
        
        // Apply period filter if requested
        auto skip = false;
        if (period_filter) {
            auto var_period = strtofloat1(entry.period);
            if (var_period == 0 || var_period >= 3) {
                skip = true;
            }
        }
        if (!skip) {
            vsx.push_back(std::move(entry));
            ++count;
        }
    }
    
    return true;
}

/// MARK: - aavso_update_required

bool aavso_update_required(const Header& head) noexcept {
    if (vsx.empty()) {
        return true;
    }
    auto sep = ang_sep(vsx[0].ra, vsx[0].dec, head.ra0, head.dec0);
    auto threshold = head.width * head.cdelt2 * kPi / 180.0;
    return !(sep < threshold);
}

/// MARK: - variable_star_annotation

void variable_star_annotation(IHttpClient& http,
                              const Header& head,
                              int annotate_mode,
                              double years_since_2000,
                              [[maybe_unused]] bool extract_visible,
                              MessageHook log) {
    auto lim_magnitude = -99.0;

    switch (annotate_mode) {
        case 0: case 1: case 2: case 3:
            // Local variable-star catalog modes — depend on data files
            // (variable.csv / variable13.csv / variable15.csv) that ASTAP
            // distributes separately. Not ported.
            if (log) log("Local variable-star catalog modes are not available "
                         "in this build; use the online modes (4–15).");
            return;
        case 4: case 8: case 12:  lim_magnitude = 11; break;
        case 5: case 9: case 13:  lim_magnitude = 13; break;
        case 6: case 19: case 14: lim_magnitude = 15; break;
        case 7: case 11: case 15: lim_magnitude = 99; break;
        default: lim_magnitude = 99; break;
    }

    // Online flow: refresh both halves only when the cache is stale.
    if (aavso_update_required(head)) {
        auto period_filt = annotate_mode < 8;
        if (!download_vsx(http, head, lim_magnitude,
                          years_since_2000, period_filt)) {
            if (log) log("No VSX data!");
            return;
        }
        if (!download_vsp(http, head, lim_magnitude)) {
            if (log) log("No VSP data!");
            return;
        }
    }
    // The caches @c vsx and @c vsp are now populated. Drawing onto the
    // canvas is the GUI overlay layer's responsibility (see
    // gui/src/annotation_scanner: project_var_stars).
}

/// MARK: - annotation_position

void annotation_position(const Header& head,
                         const std::vector<std::string>& memo_lines,
                         const std::string& aname,
                         double& ra,
                         double& dec) {
    if (head.naxis == 0) {
        return;
    }
    if (head.cd1_1 == 0) {
        return;
    }
    
    // Scan backwards through memo lines
    for (auto it = memo_lines.rbegin(); it != memo_lines.rend(); ++it) {
        const auto& line = *it;
        if (line.compare(0, 8, "ANNOTATE") != 0) {
            continue;
        }
        
        // Extract the semicolon-separated payload bounded by the closing apostrophe
        auto end = posex("'", line, 20);
        if (end == 0 || end <= 12) {
            continue;
        }
        auto payload = copy1(line, 12, end - 12);
        
        // Split on semicolons
        auto parts = std::vector<std::string>{};
        auto start = std::size_t{0};
        for (auto i = std::size_t{0}; i <= payload.size(); ++i) {
            if (i == payload.size() || payload[i] == ';') {
                parts.emplace_back(payload.substr(start, i - start));
                start = i + 1;
            }
        }
        if (parts.size() < 6) {
            continue;
        }
        if (parts[5] != aname) {
            continue;
        }
        
        // Parse the bounding box coordinates
        auto x1 = static_cast<int>(std::lround(strtofloat2(parts[0])));
        auto y1 = static_cast<int>(std::lround(strtofloat2(parts[1])));
        auto x2 = static_cast<int>(std::lround(strtofloat2(parts[2])));
        auto y2 = static_cast<int>(std::lround(strtofloat2(parts[3])));
        
        // TODO: plumb the WCS formalism selector through; default 0.
        constexpr auto formalism = 0;
        pixel_to_celestial_stub(head, (x1 + x2) / 2.0, (y1 + y2) / 2.0,
                                formalism, ra, dec);
        return;
    }
}

/// MARK: - plot_simbad

std::vector<SimbadObject> plot_simbad(std::string_view info) {
    auto out = std::vector<SimbadObject>{};
    auto lines = split_lines(info);
    
    // Single-object scratch state that flows from the "Object" line to
    // "Coordinates" to "Flux" to "Angular".
    auto name = std::string{};
    auto type = std::string{};
    auto m = 0.0;
    auto rah = 0.0;
    auto ram = 0.0;
    auto ras = 0.0;
    auto decd = 0.0;
    auto decm = 0.0;
    auto decs = 0.0;
    auto sign = 1.0;
    
    auto read_position = [&](std::string_view regel, std::size_t start, std::size_t stop) {
        auto ra1 = posex(" ", regel, start + 1);
        auto ra2 = posex(" ", regel, ra1 + 1);
        auto ra3 = posex(" ", regel, ra2 + 1);
        auto dec1 = posex(" ", regel, ra3 + 3);
        auto dec2 = posex(" ", regel, dec1 + 1);
        if (stop == 0) {
            stop = posex(" ", regel, dec2 + 2);
        }
        
        rah = strtofloat1(copy1(regel, start + 1, ra1 - start - 1));
        ram = strtofloat1(copy1(regel, ra1 + 1, ra2 - ra1 - 1));
        ras = strtofloat1(copy1(regel, ra2 + 1, ra3 - ra2 - 1));
        decd = strtofloat1(trim(copy1(regel, ra3 + 1, dec1 - ra3 - 1)));
        decm = strtofloat1(copy1(regel, dec1 + 1, dec2 - dec1 - 1));
        decs = strtofloat1(trim(copy1(regel, dec2 + 1, stop - dec2 - 1)));
        sign = (pos1("-", copy1(regel, ra3 + 1, 3)) > 0) ? -1.0 : 1.0;
    };
    
    auto count = std::size_t{5};  // skip first part
    while (count + 1 <= lines.size()) {
        auto regel = std::string_view(lines[count]);
        ++count;
        
        // Single-object reply
        if (regel.substr(0, 6) == "Object") {
            auto minnen1 = pos1("---", regel);
            auto minnen2 = posex("---", regel, minnen1 + 4);
            if (minnen1 != 0 && minnen2 != 0) {
                name = string_replace_all(
                    trim(copy1(regel, 8, minnen1 - 2 - 8)), ' ', '_');
                type = trim(copy1(regel, minnen1 + 4, minnen2 - (minnen1 + 5)));
            }
            continue;
        }
        if (regel.substr(0, 16) == "Coordinates(ICRS") {
            read_position(regel, 36, 0);
            continue;
        }
        if (regel.substr(0, 4) == "Flux") {
            auto colour = copy1(regel, 6, 1);
            if (colour == "B" || colour == "V") {
                auto end_pos = posex(" ", regel, 11);
                (void)val_double(copy1(regel, 10, end_pos - 10), m);
            }
            continue;
        }
        if (regel.substr(0, 7) == "Angular") {
            auto size = 0.0;
            auto end_pos = posex(" ", regel, 11);
            (void)val_double(copy1(regel, 15, end_pos - 10), size);
            auto obj = SimbadObject{};
            obj.name = name;
            obj.type = type;
            obj.magnitude = m;
            obj.size_arcmin = size;
            obj.ra_units = std::llround((rah + ram / 60.0 + ras / 3600.0) * 864000.0 / 24.0);
            obj.dec_units = std::llround(sign * (std::abs(decd) + decm / 60.0 + decs / 3600.0) *
                                         324000.0 / 90.0);
            out.push_back(std::move(obj));
            break;  // single-object reply terminates here
        }
        
        // List reply
        if (regel.size() >= 130 && count >= 10) {
            auto p1 = pos1("|", regel);
            auto p2 = posex("|", regel, p1 + 1);
            auto p3 = posex("|", regel, p2 + 1);
            auto p4 = posex("|", regel, p3 + 1);
            auto p5 = posex("|", regel, p4 + 1);
            auto p6 = posex("|", regel, p5 + 1);
            auto p7 = posex("|", regel, p6 + 1);
            
            if (p7 > 0) {
                read_position(regel, p3, p4);
                m = 0;
                (void)val_double(trim(copy1(regel, p6 + 1, p7 - p6 - 1)), m);  // V
                if (m == 0) {
                    (void)val_double(trim(copy1(regel, p5 + 1, p6 - p5 - 1)), m);  // B fallback
                }
                
                auto typ_l = trim(copy1(regel, p2 + 1, p3 - p2 - 1));
                auto name_l = string_replace_all(
                    trim(copy1(regel, p1 + 1, p2 - p1 - 1)), ' ', '_');
                    
                auto obj = SimbadObject{};
                obj.name = name_l;
                obj.type = typ_l;
                obj.magnitude = m;
                obj.ra_units = std::llround((rah + ram / 60.0 + ras / 3600.0) * 864000.0 / 24.0);
                obj.dec_units = std::llround(sign * (std::abs(decd) + decm / 60.0 + decs / 3600.0) *
                                             324000.0 / 90.0);
                out.push_back(std::move(obj));
            }
        }
    }
    
    return out;
}

/// MARK: - plot_vizier

std::vector<VizierObject> plot_vizier(std::string_view info,
                                      std::string_view filter,
                                      TransformGaiaFn transform_gaia) {
    auto out = std::vector<VizierObject>{};
    auto lines = split_lines(info);
    
    auto datalines = false;
    auto count = std::size_t{35};  // skip first part
    
    while (count + 1 <= lines.size()) {
        auto regel = std::string_view(lines[count]);
        ++count;
        
        if (datalines && regel.size() > 10) {
            // Parse space-separated columns: RA, DEC, G, BP, RP
            auto p1 = pos1(" ", regel);
            auto p2 = posex(" ", regel, p1 + 3);  // double-space tolerated
            auto p3 = posex(" ", regel, p2 + 3);
            auto p4 = posex(" ", regel, p3 + 3);
            
            if (p3 > 0) {
                auto rad = strtofloat1(copy1(regel, 1, p1 - 1));
                auto decd = strtofloat1(copy1(regel, p1 + 1, p2 - p1 - 1));
                auto g = strtofloat1(copy1(regel, p2 + 1, p3 - p2 - 1));
                auto bp = strtofloat1(copy1(regel, p3 + 1, p4 - p3 - 1));
                auto rp = strtofloat1(copy1(regel, p4 + 1, 99));
                
                auto themagn = transform_gaia ? transform_gaia(filter, g, bp, rp) : g;
                if (themagn != 0) {
                    auto obj = VizierObject{};
                    obj.ra_units = std::llround(rad * 864000.0 / 360.0);
                    obj.dec_units = std::llround(decd * 324000.0 / 90.0);
                    obj.magnitude = themagn;
                    out.push_back(obj);
                }
            }
        } else if (regel.substr(0, 7) == "RA_ICRS") {
            datalines = true;
            ++count;  // skip dashed separator line
        }
    }
    
    return out;
}

/// MARK: - URL builders (Simbad / Vizier)

namespace {

// Pascal `str(x:3:10, s)` → fixed-width 10-decimal string.
[[nodiscard]] std::string fixed10(double v) { return fixed(v, 10); }

// Field-of-view in arcseconds for the full image, mirroring the Pascal
// `(stopX-startX)*head.cdelt2*3600` idiom but for the entire frame.
struct FieldArcSec {
    double w_as = 0.0;
    double h_as = 0.0;
};

[[nodiscard]] FieldArcSec field_arcsec(const Header& head) noexcept {
    const auto scale = std::abs(head.cdelt2) * 3600.0;
    return { head.width * scale, head.height * scale };
}

// Build the centre coordinate fragment used by both endpoints:
// "<RA degrees, 10dp>%2B<DEC degrees, 10dp>" for north dec, "%2D" for south.
// Matches Pascal lines 14966–14970.
[[nodiscard]] std::string centre_fragment(const Header& head) {
    const auto ra_deg  = std::abs(head.ra0)  * 180.0 / kPi;
    const auto dec_deg = std::abs(head.dec0) * 180.0 / kPi;
    const auto sgn = (head.dec0 >= 0.0) ? std::string("%2B") : std::string("%2D");
    return fixed10(ra_deg) + sgn + fixed10(dec_deg);
}

}  // namespace

std::string make_simbad_url(const Header& head,
                            SimbadQuery query,
                            std::string_view maintype) {
    const auto centre = centre_fragment(head);
    const auto fov = field_arcsec(head);

    if (query == SimbadQuery::SingleAt) {
        // sim-coo single-object resolver. Pascal line 15037.
        const auto ra_deg  = std::abs(head.ra0)  * 180.0 / kPi;
        const auto dec_deg = std::abs(head.dec0) * 180.0 / kPi;
        const auto sgn = (head.dec0 >= 0.0) ? std::string("%2B")
                                            : std::string("%2D");
        const auto radius = std::max(fov.w_as, fov.h_as) / 2.0;
        return std::string("https://simbad.u-strasbg.fr/simbad/sim-coo?")
             + "Radius.unit=arcsec&Radius=" + floattostr6(radius)
             + "&Coord=" + fixed10(ra_deg) + "d" + sgn + fixed10(dec_deg) + "d"
             + "&OutputMode=LIST&output.format=ASCII";
    }

    // sim-sam box query. Criteria differs by query kind.
    auto criteria = std::string{};
    switch (query) {
        case SimbadQuery::DeepSky:         criteria = "(maintype!=*)"; break;
        case SimbadQuery::DeepSkyFiltered: criteria = "(maintype="
                                              + std::string(maintype) + ")"; break;
        case SimbadQuery::Stars:           criteria = "(maintype=*)"; break;
        case SimbadQuery::SingleAt:        break;  // handled above
    }

    return std::string("https://simbad.u-strasbg.fr/simbad/sim-sam?")
         + "submit=submit+query&maxObject=1000&Criteria=" + criteria
         + "%26+region(box," + centre + ",+"
         + floattostr4(fov.w_as) + "s+" + floattostr4(fov.h_as) + "s)"
         + "&OutputMode=LIST&output.format=ASCII";
}

std::string make_vizier_gaia_url(const Header& head, double limiting_mag) {
    const auto centre = centre_fragment(head);
    const auto fov = field_arcsec(head);
    return std::string("https://vizier.u-strasbg.fr/viz-bin/asu-txt?")
         + "-source=I/355/Gaiadr3"
         + "&-out=RA_ICRS,DE_ICRS,Gmag,BPmag,RPmag"
         + "&-c=" + centre
         + "&-c.bs=" + floattostr6(fov.w_as) + "/" + floattostr6(fov.h_as)
         + "&-out.max=10000"
         + "&Gmag=<" + floattostr4(limiting_mag);
}

} // namespace
