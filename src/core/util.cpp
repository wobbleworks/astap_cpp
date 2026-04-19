///----------------------------------------
///      @file util.cpp
///   @ingroup ASTAP++
///     @brief Implementations of miscellaneous utility helpers declared in util.h.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP); MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "util.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <format>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

///----------------------------------------
namespace astap::core {
///----------------------------------------

namespace {

/// MARK: String helpers

[[nodiscard]] std::string ascii_upper(std::string_view s) {
    auto out = std::string{};
    out.reserve(s.size());
    for (auto c : s) {
        out.push_back(static_cast<char>(
            std::toupper(static_cast<unsigned char>(c))));
    }
    return out;
}

[[nodiscard]] std::string extract_ext(std::string_view path) {
    return std::filesystem::path(path).extension().string();
}

[[nodiscard]] std::string extract_filename(std::string_view path) {
    return std::filesystem::path(path).filename().string();
}

[[nodiscard]] std::string trim(std::string_view s) {
    auto is_space = [](unsigned char c) {
        return std::isspace(c) != 0;
    };
    auto b = std::size_t{0};
    auto e = s.size();
    while (b < e && is_space(static_cast<unsigned char>(s[b]))) {
        ++b;
    }
    while (e > b && is_space(static_cast<unsigned char>(s[e - 1]))) {
        --e;
    }
    return std::string(s.substr(b, e - b));
}

// Locale-independent floating-point parser (handles inf/nan, exponent).
[[nodiscard]] bool parse_double(std::string_view s, double& out) {
    auto trimmed = trim(s);
    if (trimmed.empty()) {
        return false;
    }
    auto* first = trimmed.data();
    auto* last  = first + trimmed.size();
    auto [ptr, ec] = std::from_chars(first, last, out);
    if (ec != std::errc{}) {
        return false;
    }
    return ptr == last;
}

// Format with N decimals, dot separator, no thousands grouping.
[[nodiscard]] std::string fixed_decimals(double x, int decimals) {
    return std::format("{:.{}f}", x, decimals);
}
 
}  // namespace

/// MARK: Math

void quicksort(std::span<double> a, int lo, int hi) {
    auto Lo = lo;
    auto Hi = hi;
    const auto pivot = a[(Lo + Hi) / 2];
    do {
        while (a[Lo] < pivot) {
            ++Lo;
        }
        while (a[Hi] > pivot) {
            --Hi;
        }
        if (Lo <= Hi) {
            std::swap(a[Lo], a[Hi]);
            ++Lo;
            --Hi;
        }
    } while (Lo <= Hi);
    if (Hi > lo) {
        quicksort(a, lo, Hi);
    }
    if (Lo < hi) {
        quicksort(a, Lo, hi);
    }
}

[[nodiscard]] double smedian(std::span<const double> list, int leng) {
    if (leng == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (leng == 1) {
        return list[0];
    }
    
    // Sort a working copy.
    auto work = std::vector<double>(list.begin(), list.begin() + leng);
    quicksort(std::span<double>{work}, 0, leng - 1);
    constexpr auto min_for_three_avg = 5;
    const auto mid = (leng - 1) / 2;
    if ((leng & 1) != 0) {
        if (leng < min_for_three_avg) {
            return work[mid];
        }
        return (work[mid - 1] + work[mid] + work[mid + 1]) / 3.0;
    }
    return (work[mid] + work[mid + 1]) / 2.0;
}

void mad_median(std::span<const double> list, int leng,
                double& mad, double& median) {
    auto work = std::vector<double>(list.begin(), list.begin() + leng);
    median = smedian(std::span<const double>{work}, leng);
    for (auto i = 0; i < leng; ++i) {
        work[i] = std::abs(list[i] - median);
    }
    mad = smedian(std::span<const double>{work}, leng);
}

BestMeanResult get_best_mean(std::span<const double> list, int leng) {
    // Fast paths that skip MAD altogether (matches the Pascal original).
    if (leng <= 0) {
        return {};
    }
    if (leng == 1) {
        return {.mean = list[0], .standard_error_mean = 0.0, .count = 1};
    }
    if (leng == 2) {
        return {.mean = (list[0] + list[1]) * 0.5,
                .standard_error_mean = 0.0,
                .count = 2};
    }

    // Compute robust sigma from the median absolute deviation.
    auto mad    = 0.0;
    auto median = 0.0;
    mad_median(list, leng, mad, median);
    const auto sigma = mad * 1.4826;

    // Sum only samples within 1.5 sigma of the median.
    auto sum   = 0.0;
    auto count = 0;
    for (auto i = 0; i < leng; ++i) {
        if (std::abs(list[i] - median) < 1.5 * sigma) {
            sum += list[i];
            ++count;
        }
    }

    if (count <= 0) {
        return {};
    }

    // Final mean and SEM. SEM uses sigma (not the in-range stdev) to match the
    // original — the reasoning is that sigma is a better estimator of the true
    // population spread than the clipped sample stdev.
    const auto mean = sum / static_cast<double>(count);
    const auto sem  = sigma / std::sqrt(static_cast<double>(count));
    return {.mean = mean, .standard_error_mean = sem, .count = count};
}

[[nodiscard]] double fnmodulo(double x, double range) noexcept {
    if (x >= 0.0 && x < range) {
        return x;
    }
    
    // Truncation toward zero to match original behaviour.
    auto r = x - range * std::trunc(x / range);
    if (r < 0.0) {
        r += range;
    }
    return r;
}

/// MARK: Colour

void hsv_to_rgb(float h, float s, float v, float& r, float& g, float& b) noexcept {
    if (s == 0.0f) {
        r = v;
        g = v;
        b = v;
        return;
    }
    
    // Compute chroma and intermediate value.
    const auto c  = v * s;
    const auto h2 = h / 60.0f;
    const auto h2mod2 =
        h2 - 2.0f * static_cast<float>(static_cast<int>(h2 / 2.0f));
    const auto x = c * (1.0f - std::abs(h2mod2 - 1.0f));
    
    if (h2 < 1.0f) {
        r = c; g = x; b = 0;
    } else if (h2 < 2.0f) {
        r = x; g = c; b = 0;
    } else if (h2 < 3.0f) {
        r = 0; g = c; b = x;
    } else if (h2 < 4.0f) {
        r = 0; g = x; b = c;
    } else if (h2 < 5.0f) {
        r = x; g = 0; b = c;
    } else {
        r = c; g = 0; b = x;
    }
    
    // Add the match value.
    const auto m = v - c;
    r += m;
    g += m;
    b += m;
}

void rgb_to_hsv(float r, float g, float b, float& h, float& s, float& v) noexcept {
    const auto rgbmax = std::max(r, std::max(g, b));
    const auto rgbmin = std::min(r, std::min(g, b));
    
    if (rgbmax == rgbmin) {
        h = 0.0f;
    } else {
        if (r == rgbmax) {
            h = 60.0f * (g - b) / (rgbmax - rgbmin);
        } else if (g == rgbmax) {
            h = 60.0f * (2.0f + (b - r) / (rgbmax - rgbmin));
        } else {
            h = 60.0f * (4.0f + (r - g) / (rgbmax - rgbmin));
        }
        if (h < 0.0f) {
            h += 360.0f;
        }
    }
    
    s = (rgbmax == 0.0f) ? 0.0f : (rgbmax - rgbmin) / rgbmax;
    v = rgbmax;
}

[[nodiscard]] int rgb_to_h(float r, float g, float b) noexcept {
    const auto channel_max = std::max(r, std::max(g, b));
    const auto channel_min = std::min(r, std::min(g, b));
    auto h = 0.0f;
    
    if (channel_max != channel_min) {
        const auto d = channel_max - channel_min;
        if (r == channel_max) {
            h = (g - b) / d;
        } else if (g == channel_max) {
            h = 2.0f + (b - r) / d;
        } else {
            h = 4.0f + (r - g) / d;
        }
        h *= 60.0f;
        if (h < 0.0f) {
            h += 360.0f;
        }
    }
    
    return static_cast<int>(std::lround(h));
}

[[nodiscard]] std::array<std::uint8_t, 3> intensity_rgb(std::uint32_t colour) noexcept {
    return {
        static_cast<std::uint8_t>(colour & 0xFF),
        static_cast<std::uint8_t>((colour >> 8) & 0xFF),
        static_cast<std::uint8_t>((colour >> 16) & 0xFF),
    };
}

[[nodiscard]] std::string rgb_kelvin(float red, float blue) {
    if (!(blue >= 18.0f && red >= 18.0f)) {
        return "";
    }
    
    const auto ratio = static_cast<double>(blue) / static_cast<double>(red);
    if (!(ratio > 0.04 && ratio < 1.55)) {
        return "";
    }
    
    // Fifth-order polynomial fit.
    const auto r2 = ratio * ratio;
    const auto r3 = r2 * ratio;
    const auto r4 = r3 * ratio;
    const auto r5 = r4 * ratio;
    const auto temperature =
          4817.4  * r5
        - 4194.2  * r4
        - 7126.7  * r3
        + 12922.0 * r2
        - 2082.2  * ratio
        + 2189.8;
    return std::format("{}K", static_cast<long long>(std::lround(temperature)));
}

/// MARK: Number / string formatting

[[nodiscard]] std::string floattostr8(double x) { return fixed_decimals(x, 8); }
[[nodiscard]] std::string floattostr6(double x) { return fixed_decimals(x, 6); }
[[nodiscard]] std::string floattostr4(double x) { return fixed_decimals(x, 4); }
[[nodiscard]] std::string floattostr2(double x) { return fixed_decimals(x, 2); }

[[nodiscard]] std::string floattostrE(double x) {
    return std::format("{:e}", x);
}

[[nodiscard]] std::string inttostr5(int x) {
    return std::format("{:5}", x);
}

[[nodiscard]] int strtoint2(std::string_view s, int default_value) {
    auto value = 0.0;
    if (!parse_double(s, value)) {
        return default_value;
    }
    return static_cast<int>(std::lround(value));
}

[[nodiscard]] double strtofloat2(std::string_view s) {
    auto buf = std::string(s);
    std::replace(buf.begin(), buf.end(), ',', '.');
    auto v = 0.0;
    if (!parse_double(buf, v)) {
        return 0.0;
    }
    return v;
}

[[nodiscard]] double strtofloat1(std::string_view s) {
    auto v = 0.0;
    if (!parse_double(s, v)) {
        return 0.0;
    }
    return v;
}

void addstring(std::string& target_line, double inp) {
    // Format with a 16-wide scientific field.
    auto s = std::format("{:>16e}", inp);
    
    // Ensure the target string is at least 30 characters long.
    if (target_line.size() < 30) {
        target_line.resize(30, ' ');
    }
    
    // Right-align result into columns 11..30 (1-based).
    const auto s_len = static_cast<int>(s.size());
    for (auto i = 11; i <= 30; ++i) {
        if (i + s_len <= 30) {
            target_line[i - 1] = ' ';
        } else {
            const auto src_index = (i + s_len - 30) - 1;
            target_line[i - 1] = s[src_index];
        }
    }
}

/// MARK: Filename helpers

[[nodiscard]] bool fits_file_name(std::string_view path) {
    const auto e = ascii_upper(extract_ext(path));
    return e == ".FIT" || e == ".FITS" || e == ".FTS" || e == ".WCS";
}

[[nodiscard]] bool fits_tiff_file_name(std::string_view path) {
    const auto e = ascii_upper(extract_ext(path));
    return e == ".FIT"  || e == ".FITS" || e == ".FTS"
        || e == ".TIF"  || e == ".TIFF" || e == ".WCS";
}

[[nodiscard]] bool tiff_file_name(std::string_view path) {
    const auto e = ascii_upper(extract_ext(path));
    return e == ".TIF" || e == ".TIFF";
}

[[nodiscard]] bool check_raw_file_extension(std::string_view ext) {
    const auto e = ascii_upper(ext);
    return e == ".RAW" || e == ".CRW" || e == ".CR2" || e == ".CR3"
        || e == ".KDC" || e == ".DCR" || e == ".MRW" || e == ".ARW"
        || e == ".NEF" || e == ".NRW" || e == ".DNG" || e == ".ORF"
        || e == ".PTX" || e == ".PEF" || e == ".RW2" || e == ".SRW"
        || e == ".RAF";
}

[[nodiscard]] bool image_file_name(std::string_view path) {
    const auto e = ascii_upper(extract_ext(path));
    if (e == ".FIT"  || e == ".FITS" || e == ".FTS"
     || e == ".JPG"  || e == ".JPEG" || e == ".TIF" || e == ".PNG") {
        return true;
    }
    return check_raw_file_extension(e);
}

[[nodiscard]] int extract_exposure_from_filename(std::string_view filename) {
    const auto upper = ascii_upper(extract_filename(filename));
    
    // Look for SEC or S_ token.
    auto find = [&](std::string_view needle) {
        return upper.find(needle);
    };
    
    auto pos = find("SEC");
    if (pos == std::string::npos) {
        pos = find("S_");
    }
    
    // Need at least two chars before the keyword.
    if (pos == std::string::npos || pos < 2) {
        return 0;
    }
    
    // Walk backwards collecting digits.
    auto i = static_cast<int>(pos);
    if (upper[i - 1] == ' ') {
        --i;
    }
    auto digits = std::string{};
    while (i >= 1) {
        const auto ch = upper[i - 1];
        if (ch >= '0' && ch <= '9') {
            digits.insert(digits.begin(), ch);
        } else {
            break;
        }
        --i;
    }
    
    if (digits.empty()) {
        return 0;
    }
    auto value = 0;
    auto [p, ec] = std::from_chars(
        digits.data(), digits.data() + digits.size(), value);
    if (ec != std::errc{}) {
        return 0;
    }
    return value;
}

[[nodiscard]] int extract_temperature_from_filename(std::string_view filename) {
    const auto upper = ascii_upper(extract_filename(filename));
    
    // Search for a digit followed by 'C'.
    auto pos = std::string::npos;
    for (auto d = '0'; d <= '9'; ++d) {
        const char needle[3] = { d, 'C', 0 };
        pos = upper.find(needle);
        if (pos != std::string::npos) {
            break;
        }
    }
    if (pos == std::string::npos) {
        return 999;
    }
    
    // Walk backwards collecting digits and a leading '-'.
    auto i = static_cast<int>(pos) + 1;
    auto digits = std::string{};
    while (i >= 1) {
        const auto ch = upper[i - 1];
        if ((ch >= '0' && ch <= '9') || ch == '-') {
            digits.insert(digits.begin(), ch);
        } else {
            break;
        }
        --i;
    }
    
    if (digits.empty()) {
        return 999;
    }
    auto value = 0;
    auto [p, ec] = std::from_chars(
        digits.data(), digits.data() + digits.size(), value);
    if (ec != std::errc{}) {
        return 999;
    }
    return value;
}

[[nodiscard]] std::string extract_objectname_from_filename(std::string_view filename) {
    const auto upper = ascii_upper(extract_filename(filename));
    
    auto prefix = std::string{};
    auto i = std::string::npos;
    
    auto try_prefix = [&](std::string_view needle) {
        if (!prefix.empty()) {
            return;
        }
        const auto p = upper.find(needle);
        if (p != std::string::npos) {
            prefix = std::string(needle);
            i = p + needle.size();
        }
    };
    
    // Try known catalogue prefixes.
    try_prefix("NGC");
    try_prefix("IC");
    try_prefix("SH2-");
    try_prefix("PGC");
    try_prefix("UGC");
    try_prefix("M");
    
    if (prefix.empty() || i == std::string::npos) {
        return "";
    }
    
    // Skip a single leading space, then accumulate digits.
    if (i < upper.size() && upper[i] == ' ') {
        ++i;
    }
    auto result = prefix;
    while (i < upper.size() && upper[i] >= '0' && upper[i] <= '9') {
        result.push_back(upper[i]);
        ++i;
    }
    return result;
}

/// MARK: Light obfuscation

[[nodiscard]] std::string encrypt(std::string_view inp) {
    auto out = std::string{"1"};
    out.reserve(inp.size() + 1);
    for (auto k = std::size_t{0}; k < inp.size(); ++k) {
        const auto i = static_cast<int>(k) + 1;
        const auto v = static_cast<unsigned char>(inp[k]) + i - 11;
        out.push_back(static_cast<char>(static_cast<unsigned char>(v & 0xFF)));
    }
    return out;
}

[[nodiscard]] std::string decrypt(std::string_view inp) {
    auto out = std::string{};
    if (inp.size() < 2 || inp[0] != '1') {
        return out;
    }
    out.reserve(inp.size() - 1);
    for (auto k = std::size_t{1}; k < inp.size(); ++k) {
        const auto i = static_cast<int>(k) + 1;
        const auto v = static_cast<unsigned char>(inp[k]) - i + 11 + 1;
        out.push_back(static_cast<char>(static_cast<unsigned char>(v & 0xFF)));
    }
    return out;
}
 
} // namespace
