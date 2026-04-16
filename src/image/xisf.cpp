///----------------------------------------
///      @file xisf.cpp
///   @ingroup ASTAP++
///     @brief XISF image loader implementation.
///    @author Ported from Han Kleijn's ASTAP by John Stephen.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "xisf.h"

#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <numbers>
#include <string>
#include <string_view>
#include <vector>

///----------------------------------------
namespace astap::image {
///----------------------------------------

namespace {

constexpr auto npos = std::string::npos;

/// MARK: - String helpers

///----------------------------------------
/// @brief Find @p sub in @p s starting at offset @p start (0-based).
///----------------------------------------

[[nodiscard]] inline std::size_t find_from(std::string_view s,
                                           std::string_view sub,
                                           std::size_t      start) noexcept {
    return s.find(sub, start);
}

///----------------------------------------
/// @brief Trim leading and trailing whitespace (space, tab, CR, LF).
///----------------------------------------

[[nodiscard]] inline std::string trim(std::string_view s) {
    constexpr auto ws = std::string_view{" \t\r\n"};
    const auto a = s.find_first_not_of(ws);
    if (a == npos) {
        return {};
    }
    const auto b = s.find_last_not_of(ws);
    return std::string{s.substr(a, b - a + 1)};
}

///----------------------------------------
/// @brief Parse an integer from @p s, returning @p def on failure.
///----------------------------------------

[[nodiscard]] inline int to_int_def(std::string_view s, int def) noexcept {
    auto value = def;
    
    // Skip leading whitespace
    std::size_t i = 0;
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) {
        ++i;
    }
    
    const auto* first = s.data() + i;
    const auto* last  = s.data() + s.size();
    const auto [ptr, ec] = std::from_chars(first, last, value);
    if (ec != std::errc()) {
        return def;
    }
    return value;
}

///----------------------------------------
/// @brief Parse a double from @p s into @p out, stripping whitespace and apostrophes.
/// @param s Input string to parse.
/// @param[out] out Parsed value on success.
/// @return True on success.
///----------------------------------------

[[nodiscard]] inline bool to_double(std::string_view s, double& out) noexcept {
    // Strip whitespace and apostrophes
    std::size_t a = 0;
    while (a < s.size() && (s[a] == ' ' || s[a] == '\t' || s[a] == '\'')) {
        ++a;
    }
    auto b = s.size();
    while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\'')) {
        --b;
    }
    if (b <= a) {
        return false;
    }
    
    const auto* first = s.data() + a;
    const auto* last  = s.data() + b;
    auto v = 0.0;
    const auto [ptr, ec] = std::from_chars(first, last, v);
    if (ec != std::errc()) {
        return false;
    }
    out = v;
    return true;
}

/// MARK: - XML keyword extraction

///----------------------------------------
/// @brief Extract a string keyword value from a FITS keyword XML attribute.
/// @details Searches for `<FITSKeyword name="KEY" value="THEVALUE" .../>` and
///          returns the text up to the next `"`. Surrounding apostrophes
///          (FITS-style quoted strings) are stripped.
///   @param aline The XML header string to search.
///   @param keyword The FITS keyword name to look for.
///  @return The extracted value, or an empty string if not found.
///----------------------------------------

[[nodiscard]] static std::string extract_string_keyword(std::string_view aline,
                                                        std::string_view keyword) {
    auto needle = std::string{};
    needle.reserve(keyword.size() + 9);
    needle.append(keyword);
    needle.append("\" value=\"");
    
    auto b = aline.find(needle);
    if (b == npos) {
        return {};
    }
    b += needle.size();
    
    auto c = aline.find('"', b);
    if (c == npos) {
        return {};
    }
    
    // Strip apostrophes
    while (b < c && aline[b] == '\'') {
        ++b;
    }
    while (c > b && aline[c - 1] == '\'') {
        --c;
    }
    return std::string{aline.substr(b, c - b)};
}

///----------------------------------------
/// @brief Extract a double keyword value from an XML header.
/// @details Leaves @p value untouched if the keyword is not found or unparseable.
///   @param aline The XML header string to search.
///   @param keyword The FITS keyword name to look for.
///   @param[out] value Updated with the parsed value on success.
///----------------------------------------

static void extract_double_keyword(std::string_view aline,
                                   std::string_view keyword,
                                   double&          value) {
    auto needle = std::string{};
    needle.reserve(keyword.size() + 9);
    needle.append(keyword);
    needle.append("\" value=\"");
    
    auto b = aline.find(needle);
    if (b == npos) {
        return;
    }
    b += needle.size();
    
    auto c = aline.find('"', b);
    if (c == npos) {
        return;
    }
    
    auto v = 0.0;
    if (to_double(aline.substr(b, c - b), v)) {
        value = v;
    }
}

/// MARK: - Byte-swap helpers

///----------------------------------------
/// @brief Load a little-endian value of type @p T from a byte pointer.
/// @tparam T The numeric type to load (uint16_t, uint32_t, float, double).
///   @param p Pointer to the raw bytes.
///  @return The loaded value, byte-swapped if the host is big-endian.
///----------------------------------------

template <typename T>
[[nodiscard]] inline T load_le(const std::uint8_t* p) noexcept {
    auto v = T{};
    std::memcpy(&v, p, sizeof(T));
    if constexpr (std::endian::native == std::endian::big) {
        if constexpr (sizeof(T) == 2) {
            auto u = std::uint16_t{};
            std::memcpy(&u, &v, 2);
            u = std::byteswap(u);
            std::memcpy(&v, &u, 2);
        } else if constexpr (sizeof(T) == 4) {
            auto u = std::uint32_t{};
            std::memcpy(&u, &v, 4);
            u = std::byteswap(u);
            std::memcpy(&v, &u, 4);
        } else if constexpr (sizeof(T) == 8) {
            auto u = std::uint64_t{};
            std::memcpy(&u, &v, 8);
            u = std::byteswap(u);
            std::memcpy(&v, &u, 8);
        }
    }
    return v;
}
    
}  // namespace

/// MARK: - load_xisf

bool load_xisf(const std::filesystem::path& filen,
               Header&                      head,
               ImageArray&                  img_loaded2,
               std::vector<std::string>&    log) {
    // Reset header
    head       = Header{};
    head.naxis = 0;
    
    auto in = std::ifstream{filen, std::ios::binary};
    if (!in) {
        log.emplace_back("Error loading file!");
        return false;
    }
    
    // Read 16-byte signature block
    auto sig = std::array<std::uint8_t, 16>{};
    in.read(reinterpret_cast<char*>(sig.data()),
            static_cast<std::streamsize>(sig.size()));
    if (!in) {
        log.emplace_back("Error");
        return false;
    }
    
    // Verify XISF magic bytes
    constexpr char xisf_magic[8] = {'X', 'I', 'S', 'F', '0', '1', '0', '0'};
    if (std::memcmp(sig.data(), xisf_magic, 8) != 0) {
        log.emplace_back(
            "Error loading XISF file!! Keyword XSIF100 not found.");
        return false;
    }
    
    // 32-bit little-endian header length at offset 8
    const auto header_length =
        static_cast<std::uint32_t>(sig[8]) |
        (static_cast<std::uint32_t>(sig[9]) << 8) |
        (static_cast<std::uint32_t>(sig[10]) << 16) |
        (static_cast<std::uint32_t>(sig[11]) << 24);
        
    auto reader_position = std::size_t{16};
    
    // Read the XML header
    auto aline = std::string(header_length, '\0');
    in.read(aline.data(), static_cast<std::streamsize>(header_length));
    if (!in) {
        log.emplace_back("Error reading XISF header.");
        return false;
    }
    reader_position += header_length;
    
    // Locate <Image ...>
    const auto start_image = aline.find("<Image ");
    if (start_image == npos) {
        log.emplace_back("Error!. No <Image element.");
        return false;
    }
    
    // Reject compressed files
    if (find_from(aline, "compression=", start_image) != npos) {
        log.emplace_back("Error, can not read compressed XISF files!!");
        return false;
    }
    
    // Parse geometry="W:H:C"
    {
        const auto a = find_from(aline, "geometry=", start_image);
        if (a != npos) {
            auto b = aline.find('"', a);
            if (b == npos) {
                log.emplace_back("Malformed geometry.");
                return false;
            }
            ++b;
            
            auto c = aline.find(':', b);
            if (c == npos) {
                log.emplace_back("Malformed geometry.");
                return false;
            }
            head.width = to_int_def(trim(aline.substr(b, c - b)), 0);
            
            b = c + 1;
            c = aline.find(':', b);
            if (c == npos) {
                log.emplace_back("Malformed geometry.");
                return false;
            }
            head.height = to_int_def(trim(aline.substr(b, c - b)), 0);
            
            b = c + 1;
            c = aline.find('"', b);
            if (c == npos) {
                log.emplace_back("Malformed geometry.");
                return false;
            }
            head.naxis3 = to_int_def(trim(aline.substr(b, c - b)), 0);
        }
    }
    
    head.naxis = (head.naxis3 > 1) ? 3 : 2;
    
    // Parse location="attachment:offset:size"
    auto attachment = std::size_t{0};
    {
        const auto a = find_from(aline, "location=\"attachment", start_image);
        if (a == npos) {
            log.emplace_back("Error!. Can not read this format, no attachment");
            head.naxis = 0;
            return false;
        }
        
        auto b = aline.find(':', a);
        if (b == npos) {
            log.emplace_back("Error!. Malformed attachment.");
            head.naxis = 0;
            return false;
        }
        ++b;
        
        auto c = aline.find(':', b);
        if (c == npos) {
            log.emplace_back("Error!. Malformed attachment.");
            head.naxis = 0;
            return false;
        }
        
        auto att_d = 0.0;
        if (!to_double(trim(aline.substr(b, c - b)), att_d)) {
            log.emplace_back("Error!. Can not read this format, no attachment");
            head.naxis = 0;
            return false;
        }
        attachment = static_cast<std::size_t>(att_d);
    }
    
    // Parse sampleFormat="..."
    auto nrbits = 0;
    {
        const auto a = find_from(aline, "sampleFormat=", start_image);
        if (a == npos) {
            log.emplace_back("Can not read this format.");
            head.naxis = 0;
            return false;
        }
        
        auto b = aline.find('"', a);
        if (b == npos) {
            log.emplace_back("Can not read this format.");
            head.naxis = 0;
            return false;
        }
        ++b;
        
        auto c = aline.find('"', b);
        if (c == npos) {
            log.emplace_back("Can not read this format.");
            head.naxis = 0;
            return false;
        }
        
        const auto fmt = trim(aline.substr(b, c - b));
        if (fmt == "Float32") {
            nrbits = -32;
        } else if (fmt == "UInt16") {
            nrbits = 16;
        } else if (fmt == "UInt8") {
            nrbits = 8;
        } else if (fmt == "Float64") {
            nrbits = -64;
        } else if (fmt == "UInt32") {
            nrbits = 32;
        } else {
            log.emplace_back("Can not read this format.");
            head.naxis = 0;
            return false;
        }
    }
    
    // Set data range based on bit depth
    if (nrbits == 8) {
        head.datamin_org = 0.0;
        head.datamax_org = 255.0;
    } else {
        head.datamin_org = 0.0;
        head.datamax_org = 65535.0;
    }
    
    // Log basic header lines
    log.emplace_back("SIMPLE  =                    T");
    {
        auto s = std::string{"BITPIX  = "};
        s += std::to_string(nrbits);
        log.emplace_back(s);
    }
    {
        auto s = std::string{"NAXIS   = "};
        s += std::to_string(head.naxis);
        log.emplace_back(s);
    }
    {
        auto s = std::string{"NAXIS1  = "};
        s += std::to_string(head.width);
        log.emplace_back(s);
    }
    {
        auto s = std::string{"NAXIS2  = "};
        s += std::to_string(head.height);
        log.emplace_back(s);
    }
    if (head.naxis3 > 1) {
        auto s = std::string{"NAXIS3  = "};
        s += std::to_string(head.naxis3);
        log.emplace_back(s);
    }
    
    // Extract FITS-style keywords from XML
    head.date_obs = extract_string_keyword(aline, "DATE-OBS");
    if (head.date_obs.empty()) {
        head.date_obs = extract_string_keyword(aline, "DATE");
    }
    
    head.filter_name = extract_string_keyword(aline, "FILTER");
    
    // Extract optional metadata keywords
    [[maybe_unused]] auto bayerpat = extract_string_keyword(aline, "BAYERPAT");
    [[maybe_unused]] auto sitelong = extract_string_keyword(aline, "SITELONG");
    if (sitelong.empty()) {
        sitelong = extract_string_keyword(aline, "LONG-OBS");
    }
    [[maybe_unused]] auto sitelat = extract_string_keyword(aline, "SITELAT");
    if (sitelat.empty()) {
        sitelat = extract_string_keyword(aline, "LAT-OBS");
    }
    
    auto xbayroff = 0.0;
    auto ybayroff = 0.0;
    extract_double_keyword(aline, "XBAYROFF", xbayroff);
    extract_double_keyword(aline, "YBAYROFF", ybayroff);
    [[maybe_unused]] auto roworder = extract_string_keyword(aline, "ROWORDER");
    
    // Extract WCS CD matrix keywords
    extract_double_keyword(aline, "CD1_1", head.cd1_1);
    if (head.cd1_1 != 0.0) {
        extract_double_keyword(aline, "CD1_2", head.cd1_2);
        extract_double_keyword(aline, "CD2_1", head.cd2_1);
        extract_double_keyword(aline, "CD2_2", head.cd2_2);
        extract_double_keyword(aline, "CRPIX1", head.crpix1);
        extract_double_keyword(aline, "CRPIX2", head.crpix2);
    }
    
    // Extract temperature keywords
    auto set_temp = 0.0;
    extract_double_keyword(aline, "CCD-TEMP", set_temp);
    extract_double_keyword(aline, "SET-TEMP", set_temp);
    head.set_temperature = static_cast<int>(std::lround(set_temp));
    
    // Extract exposure keywords
    extract_double_keyword(aline, "EXPTIME ", head.exposure);
    extract_double_keyword(aline, "EXPOSURE", head.exposure);
    
    // Extract CDELT/CROTA keywords
    extract_double_keyword(aline, "CDELT2", head.cdelt2);
    if (head.cdelt2 != 0.0) {
        extract_double_keyword(aline, "CROTA1", head.crota1);
        extract_double_keyword(aline, "CROTA2", head.crota2);
        extract_double_keyword(aline, "CDELT1", head.cdelt1);
    }
    
    // Extract focal length, pressure, and ambient temperature
    auto focallen = 0.0;
    auto pressure = 0.0;
    auto focus_temp = 0.0;
    extract_double_keyword(aline, "FOCALLEN", focallen);
    extract_double_keyword(aline, "PRESSURE", pressure);
    extract_double_keyword(aline, "AOCBAROM", pressure);
    extract_double_keyword(aline, "FOCUSTEM", focus_temp);
    extract_double_keyword(aline, "FOCTEMP",  focus_temp);
    extract_double_keyword(aline, "AMB-TEMP", focus_temp);
    extract_double_keyword(aline, "AOCAMBT",  focus_temp);
    
    // Extract pixel size
    extract_double_keyword(aline, "XPIXSZ", head.xpixsz);
    
    // Compute plate scale from focal length and pixel size if no CD matrix
    if (head.cd1_1 == 0.0) {
        if (focallen != 0.0 && head.xpixsz != 0.0) {
            head.cdelt2 =
                180.0 / (std::numbers::pi * 1000.0) * head.xpixsz / focallen;
        }
        if (head.cdelt2 == 0.0) {
            extract_double_keyword(aline, "SCALE", head.cdelt2);
            head.cdelt2 /= 3600.0;
        }
        if (head.cdelt2 == 0.0) {
            extract_double_keyword(aline, "SECPIX1", head.cdelt1);
            head.cdelt1 /= 3600.0;
        }
        if (head.cdelt2 == 0.0) {
            extract_double_keyword(aline, "SECPIX2", head.cdelt2);
            head.cdelt2 /= 3600.0;
        }
    }
    
    // Extract WCS reference coordinates
    extract_double_keyword(aline, "CRVAL1", head.ra0);
    extract_double_keyword(aline, "CRVAL2", head.dec0);
    
    // Fall back to RA/DEC mount keywords if CRVAL not set
    auto ra_mount = 999.0;
    auto dec_mount = 999.0;
    extract_double_keyword(aline, "RA",  ra_mount);
    extract_double_keyword(aline, "DEC", dec_mount);
    if (ra_mount < 999.0) {
        if (head.ra0  == 0.0) {
            head.ra0  = ra_mount;
        }
        if (head.dec0 == 0.0) {
            head.dec0 = dec_mount;
        }
    }
    
    // Convert degrees to radians
    head.ra0  *= std::numbers::pi / 180.0;
    head.dec0 *= std::numbers::pi / 180.0;
    
    // Pull every <FITSKeyword .../> tag into the log
    {
        auto d = start_image;
        for (;;) {
            const auto a = find_from(aline, "<FITSKeyword name=", d);
            if (a == npos) {
                break;
            }
            const auto e = find_from(aline, "/>", a + 1);
            if (e == npos) {
                break;
            }
            
            auto b = aline.find('"', a + 1);
            if (b == npos) {
                break;
            }
            ++b;
            auto c = aline.find('"', b);
            if (c == npos) {
                break;
            }
            const auto key = trim(aline.substr(b, c - b));
            
            auto value = std::string{};
            const auto av = find_from(aline, "value=", c);
            if (av != npos && av <= e) {
                b = aline.find('"', av + 1);
                if (b == npos) {
                    break;
                }
                ++b;
                c = aline.find('"', b);
                if (c == npos) {
                    break;
                }
                value = trim(aline.substr(b, c - b));
            }
            
            auto comment = std::string{};
            const auto ac = find_from(aline, "comment=", c);
            if (ac != npos && ac <= e) {
                b = aline.find('"', ac);
                if (b == npos) {
                    break;
                }
                ++b;
                c = aline.find('"', b);
                if (c == npos) {
                    break;
                }
                comment = trim(aline.substr(b, c - b));
            }
            
            auto rec = std::string{key};
            rec.append(" = ");
            rec.append(value);
            if (!comment.empty()) {
                rec.append(" / ");
                rec.append(comment);
            }
            log.push_back(std::move(rec));
            
            d = c;
        }
    }
    
    log.emplace_back("HISTORY Imported from XISF file by the ASTAP program");
    
    // WCS conversion between (cd1_1, cd1_2, cd2_1, cd2_2) and
    // (cdelt1, cdelt2, crota1, crota2) is handled by callers in
    // astap::solving to avoid a cyclic dependency.
    
    // Clamp rotation angles
    if (head.crota2 > 999.0) {
        head.crota2 = 0.0;
    }
    if (head.crota1 > 999.0) {
        head.crota1 = head.crota2;
    }
    
    // Skip zero-padding between XML header end and attachment
    if (attachment > reader_position) {
        in.seekg(static_cast<std::streamoff>(attachment), std::ios::beg);
        if (!in) {
            head.naxis = 0;
            return false;
        }
    }
    
    // Read image data
    constexpr auto bytes_per_pix_factor = 8;
    const auto bytes_per_pix = std::abs(nrbits) / bytes_per_pix_factor;
    const auto row_bytes =
        static_cast<std::size_t>(head.width) * static_cast<std::size_t>(bytes_per_pix);
    auto row_buf = std::vector<std::uint8_t>(row_bytes);
    
    // Allocate output image array
    img_loaded2.assign(
        static_cast<std::size_t>(head.naxis3),
        std::vector<std::vector<float>>(
            static_cast<std::size_t>(head.height),
            std::vector<float>(static_cast<std::size_t>(head.width), 0.0f)));
            
    // Read pixel data plane by plane, flipping row order (XISF is top-down,
    // FITS convention is bottom-up)
    for (auto k = 0; k < head.naxis3; ++k) {
        for (auto i = head.height - 1; i >= 0; --i) {
            in.read(reinterpret_cast<char*>(row_buf.data()),
                    static_cast<std::streamsize>(row_bytes));
            if (!in) {
                head.naxis = 0;
                return false;
            }
            
            for (auto j = 0; j < head.width; ++j) {
                const auto* p =
                    row_buf.data() + static_cast<std::size_t>(j) * bytes_per_pix;
                auto v = 0.0f;
                
                switch (nrbits) {
                    case 16: {
                        v = static_cast<float>(load_le<std::uint16_t>(p));
                        break;
                    }
                    case 8: {
                        v = static_cast<float>(*p);
                        break;
                    }
                    case -32: {
                        // Float32 in [0,1] scaled to [0,65535]
                        v = 65535.0f * load_le<float>(p);
                        break;
                    }
                    case -64: {
                        v = static_cast<float>(65535.0 * load_le<double>(p));
                        break;
                    }
                    case 32: {
                        // UInt32 scaled to float with 0..65535 range
                        v = static_cast<float>(load_le<std::uint32_t>(p)) /
                            65535.0f;
                        break;
                    }
                    default:
                        break;
                }
                
                img_loaded2[static_cast<std::size_t>(k)]
                           [static_cast<std::size_t>(i)]
                           [static_cast<std::size_t>(j)] = v;
            }
        }
    }
    
    return head.naxis != 0;
}
    
} // namespace
