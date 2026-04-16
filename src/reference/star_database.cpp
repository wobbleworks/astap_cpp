///----------------------------------------
///      @file star_database.cpp
///   @ingroup ASTAP++
///     @brief HNSKY star database reader implementation.
///    @author Ported from Han Kleijn's ASTAP. MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "star_database.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <ios>
#include <string>
#include <string_view>

///----------------------------------------
namespace astap::reference {
///----------------------------------------

/// MARK: File-local helpers

namespace {

/// MARK: Ring tables

/// @brief Number of cells in each declination ring for the 290-cell tiling.
/// @details From south pole (ring 1) to north pole (last ring). Sum = 290.
inline constexpr std::array<int, 18> ring_counts_290 = {
    1,  4,  8, 12, 16, 20, 24, 28, 32,
    32, 28, 24, 20, 16, 12,  8,  4,  1
};

/// @brief Number of cells in each declination ring for the 1476-cell tiling.
/// @details From south pole (ring 1) to north pole (last ring). Sum = 1476.
inline constexpr std::array<int, 36> ring_counts_1476 = {
    1,  3,  9, 15, 21, 27, 33, 38, 43, 48, 52, 56,
    60, 63, 65, 67, 68, 69, 69, 68, 67, 65, 63, 60,
    56, 52, 48, 43, 38, 33, 27, 21, 15,  9,  3,  1
};

static_assert([] {
    auto s = 0; for (auto n : ring_counts_290) { s += n; } return s;
}() == 290, "290 ring cells must total 290");

static_assert([] {
    auto s = 0; for (auto n : ring_counts_1476) { s += n; } return s;
}() == 1476, "1476 ring cells must total 1476");

/// MARK: Filename tables

/// @brief Compact compile-time filename for a star database cell file.
struct StarFilename {
    std::array<char, 10> data{};  ///< "RRCC.1476" + NUL = 10 bytes; ".290" is shorter
    std::uint8_t size = 0;
    
    [[nodiscard]] constexpr std::string_view view() const noexcept {
        return {data.data(), size};
    }
};

///----------------------------------------
///  @brief Build a StarFilename from ring, cell, and extension at compile time.
/// @param ring 1-based ring index.
/// @param cell 1-based cell index within the ring.
/// @param ext File extension without the dot (e.g. "290", "1476").
/// @return A StarFilename containing the formatted name.
///----------------------------------------

[[nodiscard]] consteval StarFilename make_name(int ring, int cell, std::string_view ext) {
    auto f = StarFilename{};
    f.data[0] = static_cast<char>('0' + (ring / 10));
    f.data[1] = static_cast<char>('0' + (ring % 10));
    f.data[2] = static_cast<char>('0' + (cell / 10));
    f.data[3] = static_cast<char>('0' + (cell % 10));
    f.data[4] = '.';
    auto n = std::uint8_t{5};
    for (auto c : ext) {
        f.data[n++] = c;
    }
    f.size = n;
    return f;
}

///----------------------------------------
///   @brief Build a full array of StarFilenames from ring counts at compile time.
/// @tparam NRings Number of declination rings.
/// @tparam Total Total number of cells across all rings.
/// @tparam Rings Array of per-ring cell counts.
/// @tparam Ext File extension string_view.
///  @return Array of Total StarFilenames.
///----------------------------------------

template <std::size_t NRings, std::size_t Total, std::array<int, NRings> Rings,
          std::string_view const& Ext>
[[nodiscard]] consteval std::array<StarFilename, Total> build_filenames() {
    auto out = std::array<StarFilename, Total>{};
    auto idx = std::size_t{0};
    for (auto r = std::size_t{0}; r < NRings; ++r) {
        for (auto c = 1; c <= Rings[r]; ++c) {
            out[idx++] = make_name(static_cast<int>(r + 1), c, Ext);
        }
    }
    return out;
}

inline constexpr std::string_view kExt290  = "290";
inline constexpr std::string_view kExt1476 = "1476";

inline constexpr auto filenames290 =
    build_filenames<18, 290, ring_counts_290, kExt290>();
    
inline constexpr auto filenames1476 =
    build_filenames<36, 1476, ring_counts_1476, kExt1476>();
    
// Spot-check a handful against the original so a drift in the ring tables
// would be caught at build time.
static_assert(filenames290[0].view()   == "0101.290");
static_assert(filenames290[1].view()   == "0201.290");
static_assert(filenames290[5].view()   == "0301.290");
static_assert(filenames290[289].view() == "1801.290");
static_assert(filenames1476[0].view()    == "0101.1476");
static_assert(filenames1476[4].view()    == "0301.1476");
static_assert(filenames1476[1475].view() == "3601.1476");

/// MARK: Declination ring boundaries (radians)

inline constexpr auto kPi    = 3.14159265358979323846;
inline constexpr auto kTwoPi = 2.0 * kPi;
inline constexpr auto kD2R   = kPi / 180.0;

inline constexpr std::array<double, 19> dec_boundaries_290 = {
    -90.0         * kD2R,
    -85.23224404  * kD2R,   // arcsin(1-1/289)
    -75.66348756  * kD2R,   // arcsin(1-(1+8)/289)
    -65.99286637  * kD2R,   // arcsin(1-(1+8+16)/289)
    -56.14497387  * kD2R,
    -46.03163067  * kD2R,
    -35.54307745  * kD2R,
    -24.53348115  * kD2R,
    -12.79440589  * kD2R,
     0.0,
     12.79440589  * kD2R,
     24.53348115  * kD2R,
     35.54307745  * kD2R,
     46.03163067  * kD2R,
     56.14497387  * kD2R,
     65.99286637  * kD2R,
     75.66348756  * kD2R,
     85.23224404  * kD2R,
     90.0         * kD2R,
};

inline constexpr std::array<double, 37> dec_boundaries_1476 = {
    -90.0         * kD2R,
    -87.42857143  * kD2R,
    -82.28571429  * kD2R,
    -77.14285714  * kD2R,
    -72.0         * kD2R,
    -66.85714286  * kD2R,
    -61.71428571  * kD2R,
    -56.57142857  * kD2R,
    -51.42857143  * kD2R,
    -46.28571429  * kD2R,
    -41.14285714  * kD2R,
    -36.0         * kD2R,
    -30.85714286  * kD2R,
    -25.71428571  * kD2R,
    -20.57142857  * kD2R,
    -15.42857143  * kD2R,
    -10.28571429  * kD2R,
    -5.142857143  * kD2R,
     0.0,
     5.142857143  * kD2R,
     10.28571429  * kD2R,
     15.42857143  * kD2R,
     20.57142857  * kD2R,
     25.71428571  * kD2R,
     30.85714286  * kD2R,
     36.0         * kD2R,
     41.14285714  * kD2R,
     46.28571429  * kD2R,
     51.42857143  * kD2R,
     56.57142857  * kD2R,
     61.71428571  * kD2R,
     66.85714286  * kD2R,
     72.0         * kD2R,
     77.14285714  * kD2R,
     82.28571429  * kD2R,
     87.42857143  * kD2R,
     90.0         * kD2R,
};

/// MARK: Module-local state

/// @brief Module-local file stream for the currently open cell file.
std::ifstream thefile_stars;

///----------------------------------------
/// @brief Fractional part of a non-negative value (floor-based).
/// @param x Input value (expected non-negative in context).
/// @return Fractional part of x.
///----------------------------------------

[[nodiscard]] inline double frac_pos(double x) noexcept {
    return x - std::floor(x);
}

///----------------------------------------
/// @brief Write a diagnostic message to stderr.
/// @details TODO: route through the application log once the hosting module
///          is ported.
/// @param msg Message to output.
///----------------------------------------

void memo2_message(std::string_view msg) {
    std::fwrite(msg.data(), 1, msg.size(), stderr);
    std::fputc('\n', stderr);
}

/// @brief Warning string set by select_star_database. TODO: migrate to the ported owner.
std::string warning_str;

/// MARK: Area lookup helpers

void area_and_boundaries(double ra1, double dec1,
                         int& area_nr,
                         double& spaceE, double& spaceW,
                         double& spaceN, double& spaceS);
                         
void area_and_boundaries1476(double ra1, double dec1,
                             int& area_nr,
                             double& spaceE, double& spaceW,
                             double& spaceN, double& spaceS);
                             
///----------------------------------------
///      @brief Compute area number and boundary offsets for a generic ring.
///   @details Shared between both 290 and 1476 tilings.
///      @param ra1 Right ascension in radians.
///      @param dec1 Declination in radians.
///      @param cos_dec Precomputed cos(dec1).
///      @param base Sum of cell counts in rings below this ring.
///      @param nra Number of RA cells in this ring.
///      @param dec_lo Southern boundary of this ring (radians).
///      @param dec_hi Northern boundary of this ring (radians).
/// @param[out] area_nr 1-based area index.
/// @param[out] spaceE Angular distance to the eastern boundary.
/// @param[out] spaceW Angular distance to the western boundary.
/// @param[out] spaceN Angular distance to the northern boundary.
/// @param[out] spaceS Angular distance to the southern boundary.
///----------------------------------------

inline void fill_generic_ring(double ra1, double dec1, double cos_dec,
                              int base, int nra, double dec_lo, double dec_hi,
                              int& area_nr,
                              double& spaceE, double& spaceW,
                              double& spaceN, double& spaceS) {
    const auto rot = ra1 * static_cast<double>(nra) / kTwoPi;
    area_nr = base + 1 + static_cast<int>(std::trunc(rot));
    spaceS = dec1 - dec_lo;
    spaceN = dec_hi - dec1;
    const auto step = kTwoPi / static_cast<double>(nra);
    spaceW = step * frac_pos(rot) * cos_dec;          // RA decreases to the west
    spaceE = step * (1.0 - frac_pos(rot)) * cos_dec;
}

void area_and_boundaries(double ra1, double dec1,
                         int& area_nr,
                         double& spaceE, double& spaceW,
                         double& spaceN, double& spaceS) {
    const auto cos_dec = std::cos(dec1);
    const auto& b = dec_boundaries_290;
    
    // Cumulative cell count from the south pole, mirroring the nested
    // if/else ladder in the original source (lines 2139-2317).
    if (dec1 > b[17]) {
        area_nr = 290;  // celestial north pole
        spaceS = dec1 - b[17];
        spaceN = b[18] - b[17];
        spaceW = kTwoPi;
        spaceE = kTwoPi;
        return;
    }
    
    if (dec1 > b[16]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4+8+12+16+20+24+28+32+32+28+24+20+16+12+8,  4, b[16], b[17], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[15]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4+8+12+16+20+24+28+32+32+28+24+20+16+12,    8, b[15], b[16], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[14]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4+8+12+16+20+24+28+32+32+28+24+20+16,      12, b[14], b[15], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[13]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4+8+12+16+20+24+28+32+32+28+24+20,         16, b[13], b[14], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[12]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4+8+12+16+20+24+28+32+32+28+24,            20, b[12], b[13], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[11]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4+8+12+16+20+24+28+32+32+28,               24, b[11], b[12], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[10]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4+8+12+16+20+24+28+32+32,                  28, b[10], b[11], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[ 9]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4+8+12+16+20+24+28+32,                     32, b[ 9], b[10], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[ 8]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4+8+12+16+20+24+28,                        32, b[ 8], b[ 9], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[ 7]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4+8+12+16+20+24,                           28, b[ 7], b[ 8], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[ 6]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4+8+12+16+20,                              24, b[ 6], b[ 7], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[ 5]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4+8+12+16,                                 20, b[ 5], b[ 6], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[ 4]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4+8+12,                                    16, b[ 4], b[ 5], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[ 3]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4+8,                                       12, b[ 3], b[ 4], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[ 2]) { fill_generic_ring(ra1, dec1, cos_dec, 1+4,                                          8, b[ 2], b[ 3], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    if (dec1 > b[ 1]) { fill_generic_ring(ra1, dec1, cos_dec, 1,                                            4, b[ 1], b[ 2], area_nr, spaceE, spaceW, spaceN, spaceS); return; }
    
    // South pole cap.
    area_nr = 1;
    spaceS = b[1] - b[0];
    spaceN = b[1] - dec1;
    spaceW = kTwoPi;
    spaceE = kTwoPi;
}

/// @brief Prefix sums of ring_counts_1476 (0-based boundary index -> cells below).
inline constexpr auto cumulative_1476 = [] {
    auto a = std::array<int, 37>{};
    a[0] = 0;
    for (auto i = std::size_t{0}; i < ring_counts_1476.size(); ++i) {
        a[i + 1] = a[i] + ring_counts_1476[i];
    }
    return a;
}();

static_assert(cumulative_1476[36] == 1476);

void area_and_boundaries1476(double ra1, double dec1,
                             int& area_nr,
                             double& spaceE, double& spaceW,
                             double& spaceN, double& spaceS) {
    const auto cos_dec = std::cos(dec1);
    const auto& b = dec_boundaries_1476;
    const auto& c = cumulative_1476;
    
    if (dec1 > b[35]) {
        area_nr = 1476;  // celestial north pole
        spaceS = dec1 - b[35];
        spaceN = b[36] - b[35];
        spaceW = kTwoPi;
        spaceE = kTwoPi;
        return;
    }
    
    // Rings 2..35 (indices 1..34 into ring_counts_1476, upper-boundary index
    // = ring index). ring_counts[r-1] is the number of RA cells in ring r.
    for (auto r = 35; r >= 2; --r) {
        if (dec1 > b[r - 1]) {
            const auto nra = ring_counts_1476[static_cast<std::size_t>(r - 1)];
            fill_generic_ring(ra1, dec1, cos_dec,
                              c[static_cast<std::size_t>(r - 1)],
                              nra, b[r - 1], b[r],
                              area_nr, spaceE, spaceW, spaceN, spaceS);
                              
            // Preserve a transcription quirk in the original source (line 2513):
            // in the ring just north of the equator (dec_boundaries1476[17])
            // spaceN uses dec_boundaries1476[19] rather than [18]. We replicate
            // that here to stay byte-identical with the original behaviour.
            if (r == 18) {
                spaceN = b[19] - dec1;
            }
            return;
        }
    }
    
    // South pole cap.
    area_nr = 1;
    spaceS = b[1] - b[0];
    spaceN = b[1] - dec1;
    spaceW = kTwoPi;
    spaceE = kTwoPi;
}
 
}  // namespace

/// MARK: Module-level state definitions

double                 cos_telescope_dec{};
std::array<char, 110>  database2{};
std::string            name_database{};
int                    cache_valid_pos{0};
int                    database_type{1476};
bool                   file_open{false};
double                 area2{1.0 * kD2R};
int                    old_area{9999999};
int                    cache_position{0};

// Defined here as a weak fallback so this TU links standalone; the real owner
// (astap main) will provide the definitive value at link time or by assignment
// after process start.
std::filesystem::path database_path{};

namespace detail {
int                              record_size{11};
std::vector<std::uint8_t>        cache_array{};
int                              cache_size{0};
std::int8_t                      dec9_storage{0};
std::array<std::uint8_t, 11>     buf2{};
}  // namespace detail

/// MARK: fill_cache

void fill_cache(int offset, int bytes) {
    if (bytes <= 0) {
        return;
    }
    
    thefile_stars.read(
        reinterpret_cast<char*>(detail::cache_array.data()) + offset,
        bytes);
}

/// MARK: Public API

void find_areas(double ra1, double dec1, double fov,
                int& area1, int& area2_out, int& area3, int& area4,
                double& frac1, double& frac2, double& frac3, double& frac4) {
    // Clamp FOV to just under the cell size so we never skip beyond a
    // neighbouring cell.
    if (database_type == 290) {
        fov = std::min(fov, 9.53      * kD2R);
    } else {
        fov = std::min(fov, 5.142857  * kD2R);
    }
    
    const auto fov_half = fov / 2.0;
    const auto dec_cornerN = dec1 + fov_half;  // >+pi/2 is still north pole
    const auto dec_cornerS = dec1 - fov_half;
    
    auto wrap_ra = [](double ra) {
        if (ra < 0.0) {
            ra += kTwoPi;
        } else if (ra >= kTwoPi) {
            ra -= kTwoPi;
        }
        return ra;
    };
    
    const auto ra_cornerWN = wrap_ra(ra1 - fov_half / std::cos(dec_cornerN));
    const auto ra_cornerEN = wrap_ra(ra1 + fov_half / std::cos(dec_cornerN));
    const auto ra_cornerWS = wrap_ra(ra1 - fov_half / std::cos(dec_cornerS));
    const auto ra_cornerES = wrap_ra(ra1 + fov_half / std::cos(dec_cornerS));
    
    auto spaceE = 0.0;
    auto spaceW = 0.0;
    auto spaceN = 0.0;
    auto spaceS = 0.0;
    auto lookup = (database_type == 290) ? &area_and_boundaries
                                          : &area_and_boundaries1476;
                                          
    // Corner 1: north-east
    lookup(ra_cornerEN, dec_cornerN, area1, spaceE, spaceW, spaceN, spaceS);
    frac1 = std::min(spaceW, fov) * std::min(spaceS, fov) / (fov * fov);
    
    // Corner 2: north-west
    lookup(ra_cornerWN, dec_cornerN, area2_out, spaceE, spaceW, spaceN, spaceS);
    frac2 = std::min(spaceE, fov) * std::min(spaceS, fov) / (fov * fov);
    
    // Corner 3: south-east
    lookup(ra_cornerES, dec_cornerS, area3, spaceE, spaceW, spaceN, spaceS);
    frac3 = std::min(spaceW, fov) * std::min(spaceN, fov) / (fov * fov);
    
    // Corner 4: south-west
    lookup(ra_cornerWS, dec_cornerS, area4, spaceE, spaceW, spaceN, spaceS);
    frac4 = std::min(spaceE, fov) * std::min(spaceN, fov) / (fov * fov);
    
    // De-duplicate overlapping corners that map to the same cell.
    if (area2_out == area1) { area2_out = 0; frac2 = 0; }
    if (area3 == area1)     { area3 = 0;     frac3 = 0; }
    if (area4 == area1)     { area4 = 0;     frac4 = 0; }
    if (area3 == area2_out) { area3 = 0;     frac3 = 0; }
    if (area4 == area2_out) { area4 = 0;     frac4 = 0; }
    if (area4 == area3)     { area4 = 0;     frac4 = 0; }
    
    // Drop coverage below 1%.
    if (frac1 < 0.01) { area1     = 0; frac1 = 0; }
    if (frac2 < 0.01) { area2_out = 0; frac2 = 0; }
    if (frac3 < 0.01) { area3     = 0; frac3 = 0; }
    if (frac4 < 0.01) { area4     = 0; frac4 = 0; }
}

[[nodiscard]] bool select_star_database(const std::string& database_in, double fov) {
    auto warning = false;
    database_type = 1476;
    
    // Lowercase the caller's choice.
    auto database = database_in;
    std::transform(database.begin(), database.end(), database.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                   
    auto typ = database.empty() ? 'a' : database[0];
    
    auto exists = [](const std::filesystem::path& p) {
        auto ec = std::error_code{};
        return std::filesystem::exists(p, ec);
    };
    
    auto has = [&](std::string_view suffix) {
        return exists(database_path / (database + std::string(suffix)));
    };
    
    if (typ != 'a') {
        if (typ == 'w') {
            if (has("_0101.001")) {
                name_database = database;
                database_type = 1;
                return true;
            }
        } else if (typ == 'd' || typ == 'v' || typ == 'h') {
            if (has("_0101.1476")) {
                name_database = database;
                return true;
            }
        }
        // No else: there is both v50 (1476) and v17 (290).
        if (typ == 'v' || typ == 'g') {
            if (has("_0101.290")) {
                name_database = database;
                database_type = 290;
                return true;
            }
        }
    }
    
    auto exists_prefix = [&](std::string_view prefix) {
        return exists(database_path / std::string(prefix));
    };
    
    if (fov > 20) {
        if (exists_prefix("w08_0101.001")) {
            name_database = "w08";
            database_type = 1;
            return true;
        }
        memo2_message("Could not find w08 star database. Will try with an other database.");
    }
    
    if (fov > 6 && exists_prefix("g05_0101.290"))      { name_database = "g05"; database_type = 290; }
    else if (fov > 6 && exists_prefix("v05_0101.290")) { name_database = "v05"; database_type = 290; }
    else if (fov > 6 && exists_prefix("v17_0101.290")) { name_database = "v17"; database_type = 290; warning = true; }
    else if (exists_prefix("d80_0101.1476"))           { name_database = "d80"; }
    else if (exists_prefix("d50_0101.1476"))           { name_database = "d50"; }
    else if (exists_prefix("v50_0101.1476"))           { name_database = "v50"; }
    else if (exists_prefix("d20_0101.1476"))           { name_database = "d20"; }
    else if (exists_prefix("d05_0101.1476"))           { name_database = "d05"; }
    else if (exists_prefix("g05_0101.290"))            { name_database = "g05"; database_type = 290; }
    else if (exists_prefix("v05_0101.290"))            { name_database = "v05"; database_type = 290; }
    else if (exists_prefix("h18_0101.1476"))           { name_database = "h18"; warning = true; }
    else if (exists_prefix("g18_0101.290"))            { name_database = "g18"; warning = true; database_type = 290; }
    else if (exists_prefix("h17_0101.1476"))           { name_database = "h17"; warning = true; }
    else if (exists_prefix("v17_0101.290"))            { name_database = "v17"; database_type = 290; warning = true; }
    else if (exists_prefix("g17_0101.290"))            { name_database = "g17"; database_type = 290; warning = true; }
    else {
        // TODO: surface this to the GUI via the caller.
        memo2_message("No star database found! Download and install one star database.");
        return false;
    }
    
    if (warning) {
        warning_str = "Old database!";
    }
    return true;
}

void close_star_database() noexcept {
    if (file_open) {
        thefile_stars.close();
        file_open = false;
    }
}

[[nodiscard]] bool open_database(double telescope_dec, int area290) {
    // Cache cos(telescope_dec) for readdatabase290.
    cos_telescope_dec = std::cos(telescope_dec);
    
    if (area290 != old_area || !file_open) {
        close_star_database();
        
        auto namefile = name_database + "_";
        if (database_type == 290) {
            // area290 is 1-based at the API boundary.
            namefile += std::string(filenames290[static_cast<std::size_t>(area290 - 1)].view());
        } else {
            namefile += std::string(filenames1476[static_cast<std::size_t>(area290 - 1)].view());
        }
        
        const auto full_path = database_path / namefile;
        thefile_stars.open(full_path, std::ios::binary | std::ios::in);
        if (!thefile_stars.is_open()) {
            return false;
        }
        file_open = true;
        
        cache_valid_pos = 0;  // new file
        
        // Read 110-byte header (10 rows x 11 chars).
        thefile_stars.read(database2.data(), 110);
        
        if (database2[109] == ' ') {
            detail::record_size = 11;  // default
        } else {
            detail::record_size = static_cast<int>(static_cast<unsigned char>(database2[109]));
        }
        
        // File size minus the 110-byte header == number of record bytes.
        const auto cur = thefile_stars.tellg();
        thefile_stars.seekg(0, std::ios::end);
        const auto end = thefile_stars.tellg();
        thefile_stars.seekg(cur, std::ios::beg);
        detail::cache_size = static_cast<int>(end) - 110;
        
        if (detail::cache_size > static_cast<int>(detail::cache_array.size())) {
            // Clear first to avoid copying the old contents during resize.
            detail::cache_array.clear();
            detail::cache_array.shrink_to_fit();
            detail::cache_array.resize(static_cast<std::size_t>(detail::cache_size));
        }
        
        old_area = area290;
    }
    // else: file for this area is already open; keep cache warm.
    
    cache_position = 0;
    return true;
}
 
} // namespace
