///----------------------------------------
///      @file star_database.h
///   @ingroup ASTAP++
///     @brief HNSKY star database reader (.290 and .1476 formats).
///   @details The .290 format divides the sky into 290 equal-area cells (except
///            for the poles, which are half-size) across 18 declination rings.
///            The .1476 format does the same across 36 rings of ~5.14 degrees
///            each. Each cell is a separate file named <dbname>_RRCC.<ext> where
///            RR is the 1-based ring index and CC is the 1-based cell-in-ring
///            index, both two digits, zero padded. Each file has a 110-byte
///            textual header followed by packed star records of 5, 6, 7, 10 or
///            11 bytes. ASTAP only supports record sizes 5 and 6 for solving.
///
///            --- Record layout ---
///            Each record is 5 or 6 bytes (packed, little-endian-ish byte
///            ordering with an explicit DEC9 cache):
///
///              hnskyhdr1476_5  (5 bytes):  ra7, ra8, ra9, dec7, dec8
///              hnskyhdr1476_6  (6 bytes):  ra7, ra8, ra9, dec7, dec8, B_R
///
///            RA is a 3-byte unsigned integer (ra7 | ra8<<8 | ra9<<16) giving
///            angle as ra = raw * (2*pi / (2^24 - 1)) radians.
///
///            DEC9 is the high byte of a 3-byte two's-complement DEC. Since
///            stars are sorted by magnitude within each cell and clustered by
///            DEC9, DEC9 is stored once in a preceding header-record
///            (ra_raw == 0xFFFFFF) with:
///              dec9 = dec7 - 128  (shortint, recovered from a biased byte)
///              mag  = (dec8 - 16)/10  magnitudes (raw byte stored as mag*10+16)
///            Subsequent normal records reuse the cached dec9:
///              dec = ((dec9<<16) | (dec8<<8) | dec7) * (pi/2 / (2^23-1))  rad
///
///            Record size 6 adds a signed Gaia Bp-Rp colour byte (B_R, *10).
///
///            --- Area selection ---
///            find_areas() returns up to four area indices (1..290 or 1..1476)
///            covering a telescope FOV, along with the fractional area coverage
///            of each. open_database opens one of these files and
///            select_star_database picks a DB by name or by FOV.
///            readdatabase290 streams one star at a time and is inlined because
///            it is called in the innermost loop (once per star).
///
///            --- Index convention ---
///            The original source uses 1-based indices for area numbers (1..290,
///            1..1476) and those indices leak into the public API (find_areas
///            outputs, open_database input). This port preserves the 1-based
///            area indexing at the API boundary (value 0 == "no area"). Internal
///            std::array lookups use `area - 1`.
///    @author Ported from Han Kleijn's ASTAP. MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

///----------------------------------------
namespace astap::reference {
///----------------------------------------

/// MARK: Packed on-disk star record types

#pragma pack(push, 1)

/// @brief Most compact star record (5 bytes). Used by Gaia .1476 databases.
struct hnskyhdr1476_5 {
    std::uint8_t ra7;
    std::uint8_t ra8;
    std::uint8_t ra9;
    std::uint8_t dec7;
    std::uint8_t dec8;
};

/// @brief 6-byte record; same layout as _5 but with an extra signed Gaia Bp-Rp byte.
struct hnskyhdr1476_6 {
    std::uint8_t ra7;
    std::uint8_t ra8;
    std::uint8_t ra9;
    std::uint8_t dec7;
    std::uint8_t dec8;
    std::int8_t  B_R;  ///< Gaia (Bp-Rp) * 10
};

#pragma pack(pop)

static_assert(sizeof(hnskyhdr1476_5) == 5, "hnskyhdr1476_5 must be packed to 5 bytes");
static_assert(sizeof(hnskyhdr1476_6) == 6, "hnskyhdr1476_6 must be packed to 6 bytes");

/// MARK: Module-level state

/// @brief cos(telescope_dec); must be set before each read series so readdatabase290
///        can cheaply reject stars outside the FOV.
extern double cos_telescope_dec;

/// @brief 110-byte database header text read from the current file.
extern std::array<char, 110> database2;

/// @brief Name (prefix) of the currently selected database, e.g. "g17", "d80".
extern std::string name_database;

/// @brief Byte offset within cache_array up to which the cache has been filled.
extern int cache_valid_pos;

/// @brief 290 or 1476 depending on which tiling the active database uses.
extern int database_type;

/// @brief True while a star database file is held open between calls.
extern bool file_open;

/// @brief Search half-size (radians). Default pi/180 (1 deg). Kept for API parity;
///        the reader uses field_diameter directly.
extern double area2;

/// @brief Sentinel for "no previously opened area".
extern int old_area;

/// @brief Current byte offset within cache_array for streaming reads.
extern int cache_position;

/// @brief Path to the directory containing star database cell files.
/// @details Owned by another (not-yet-ported) module and is resolved from
///          application settings. Declared here so select_star_database and
///          open_database can locate cell files. The defining TU should provide it.
extern std::filesystem::path database_path;

/// MARK: Public API

///----------------------------------------
///    @brief Select a star database by prefix (or auto-select by FOV).
///   @param database Database name prefix to try, or empty for auto-selection.
///   @param fov Field of view in radians for auto-selection heuristics.
///  @return False if no suitable database files are present under database_path.
///----------------------------------------

[[nodiscard]] bool select_star_database(const std::string& database, double fov);

///----------------------------------------
///    @brief Resolve the photometric passband from a FITS filter name + user-chosen
///           reference database.
/// @details Mirrors the Pascal original in @c unit_annotation.pas. When
///          @p reference_database contains "auto" or "Local", the passband
///          is inferred from @p filter_str using ASTAP's filter-name conventions
///          (empty/CV → BP, Sloan S* filters → SG/SR/SI, Johnson-Cousins
///          G/V → V, B → B, R → R; otherwise BP). When the user has explicitly
///          selected a non-auto reference (e.g. "Gaia BP", "Johnson V"), the
///          passband is parsed directly from @p reference_database using
///          priority-ordered substring matching: BP, V, B, R, SG, SR, SI.
///          Returns "??" if nothing matches.
///    @param filter_str Filter name from the FITS header (e.g. "G", "V", "Sloan g").
///    @param reference_database User-selected reference database string.
/// @return Two-character passband identifier ("BP", "V", "B", "R", "SG", "SR",
///         "SI", or "??").
///----------------------------------------

[[nodiscard]] std::string get_database_passband(std::string_view filter_str,
                                                 std::string_view reference_database);

///----------------------------------------
///      @brief For a telescope position and FOV, identify up to four area
///             indices covering the field of view.
///   @details Area indices are 1-based (1..290 or 1..1476 depending on
///            database_type). Zero indicates "unused".
///      @param ra1 Telescope right ascension in radians.
///      @param dec1 Telescope declination in radians.
///      @param fov Field of view in radians.
/// @param[out] area1 First area index.
/// @param[out] area2 Second area index (0 if unused).
/// @param[out] area3 Third area index (0 if unused).
/// @param[out] area4 Fourth area index (0 if unused).
/// @param[out] frac1 Fractional coverage of area1.
/// @param[out] frac2 Fractional coverage of area2.
/// @param[out] frac3 Fractional coverage of area3.
/// @param[out] frac4 Fractional coverage of area4.
///----------------------------------------

void find_areas(double ra1, double dec1, double fov,
                int& area1, int& area2, int& area3, int& area4,
                double& frac1, double& frac2, double& frac3, double& frac4);
                
///----------------------------------------
///    @brief Open (and cache) the file backing a specific area.
///   @param telescope_dec Telescope declination in radians (used to set cos_telescope_dec).
///   @param area290 1-based area index.
///  @return False if the file cannot be opened.
///----------------------------------------

[[nodiscard]] bool open_database(double telescope_dec, int area290);

///----------------------------------------
/// @brief Close any currently open database file.
///----------------------------------------

void close_star_database() noexcept;

///----------------------------------------
/// @brief Forward declaration of the inline hot-path reader.
///----------------------------------------

inline bool readdatabase290(double telescope_ra, double telescope_dec,
                            double field_diameter,
                            double& ra2, double& dec2, double& mag2, double& Bp_Rp);
                            
///----------------------------------------
///  @brief Non-inline cache-fill helper referenced from the inlined reader.
/// @details Defined in the .cpp because it owns the std::ifstream that backs
///          the cache. Declared here so the call site compiles cleanly.
///  @param offset Byte offset within cache_array to write into.
///  @param bytes Number of bytes to read from the open file stream.
///----------------------------------------

void fill_cache(int offset, int bytes);

/// MARK: Internal streaming state (detail)

namespace detail {

/// @brief Record size in bytes (5, 6, 7, 10 or 11) as read from database2[109].
extern int record_size;

/// @brief Whole-file cache of star records (after the 110-byte header).
extern std::vector<std::uint8_t> cache_array;

/// @brief Valid size of cache_array (== file size minus 110).
extern int cache_size;

/// @brief Most recently decoded DEC9 byte from a header-record (signed).
extern std::int8_t dec9_storage;

/// @brief One-record scratch buffer (max record size is 11).
extern std::array<std::uint8_t, 11> buf2;
    
}  // namespace detail

/// MARK: Inline hot-path reader

///----------------------------------------
///    @brief Stream one usable star at a time from the current cached file.
///   @details Returns false when the file is exhausted.
///
///            Preconditions:
///            - open_database() succeeded for the target area.
///            - cos_telescope_dec == cos(telescope_dec).
///
///            Outputs (on true):
///            - ra2, dec2: in radians
///            - mag2: magnitude * 10 (raw header-record byte minus 16)
///            - Bp_Rp: Gaia Bp-Rp * 10 for record size 6; unmodified otherwise
///
///            Only record sizes 5 and 6 are supported.
///      @param telescope_ra Telescope right ascension in radians.
///      @param telescope_dec Telescope declination in radians.
///      @param field_diameter Field of view diameter in radians.
/// @param[out] ra2 Decoded star right ascension in radians.
/// @param[out] dec2 Decoded star declination in radians.
/// @param[out] mag2 Decoded star magnitude * 10.
/// @param[out] Bp_Rp Decoded Gaia Bp-Rp * 10 (record size 6 only).
///    @return True if a star was decoded, false if the file is exhausted.
///----------------------------------------

inline bool readdatabase290(double telescope_ra, double telescope_dec,
                            double field_diameter,
                            double& ra2, double& dec2, double& mag2, double& Bp_Rp) {
    constexpr auto kPi        = 3.14159265358979323846;
    constexpr auto kTwoPi     = 2.0 * kPi;
    constexpr auto kRaScale   = kTwoPi / ((256.0 * 256.0 * 256.0) - 1.0);
    constexpr auto kDecScale  = (kPi * 0.5) / ((128.0 * 256.0 * 256.0) - 1.0);
    constexpr auto kBlockSize = 5 * 6 * 4 * 1024;  // multiple of 5 and 6
    
    // These live in detail:: so the non-inline .cpp owns their storage while
    // the body can still see them.
    using detail::record_size;
    using detail::cache_array;
    using detail::cache_size;
    using detail::dec9_storage;
    using detail::buf2;
    
    auto ra_raw = 0;
    auto delta_ra = 0.0;
    auto header_record = false;
    
    do {
        do {
            if (cache_position >= cache_size) {
                return false;  // end of file
            }
            
            if (cache_position >= cache_valid_pos) {
                // Top up the cache from the open file stream.
                auto block_to_read = std::min(cache_size, kBlockSize);
                block_to_read = std::min(block_to_read, cache_size - cache_valid_pos);
                // Streaming read is delegated to a non-inline helper that
                // owns the ifstream; see star_database.cpp.
                fill_cache(cache_valid_pos, block_to_read);
                cache_valid_pos += block_to_read;
            }
            
            // Copy one record out of the cache into the scratch buffer.
            std::copy_n(cache_array.begin() + cache_position, record_size, buf2.begin());
            cache_position += record_size;
            
            header_record = false;
            
            if (record_size == 5) {
                const auto* p5 = reinterpret_cast<const hnskyhdr1476_5*>(buf2.data());
                ra_raw = (p5->ra7) | (p5->ra8 << 8) | (p5->ra9 << 16);
                if (ra_raw == 0xFFFFFF) {
                    // Special magnitude/DEC9 header record.
                    mag2 = static_cast<double>(p5->dec8) - 16.0;
                    dec9_storage = static_cast<std::int8_t>(p5->dec7 - 128);
                    header_record = true;
                } else {
                    ra2 = static_cast<double>(ra_raw) * kRaScale;
                    const auto dec_raw = (static_cast<int>(dec9_storage) << 16)
                                       + (static_cast<int>(p5->dec8) << 8)
                                       +  static_cast<int>(p5->dec7);
                    dec2 = static_cast<double>(dec_raw) * kDecScale;
                }
            } else if (record_size == 6) {
                const auto* p6 = reinterpret_cast<const hnskyhdr1476_6*>(buf2.data());
                ra_raw = (p6->ra7) | (p6->ra8 << 8) | (p6->ra9 << 16);
                if (ra_raw == 0xFFFFFF) {
                    mag2 = static_cast<double>(p6->dec8) - 16.0;
                    dec9_storage = static_cast<std::int8_t>(p6->dec7 - 128);
                    header_record = true;
                } else {
                    ra2 = static_cast<double>(ra_raw) * kRaScale;
                    const auto dec_raw = (static_cast<int>(dec9_storage) << 16)
                                       + (static_cast<int>(p6->dec8) << 8)
                                       +  static_cast<int>(p6->dec7);
                    dec2 = static_cast<double>(dec_raw) * kDecScale;
                    Bp_Rp = static_cast<double>(p6->B_R);
                }
            }
            // Other record sizes fall through with header_record == false,
            // matching the original's empty case.
        } while (header_record);
        
        delta_ra = std::abs(ra2 - telescope_ra);
        if (delta_ra > kPi) {
            delta_ra = kTwoPi - delta_ra;
        }
    } while (!((delta_ra * cos_telescope_dec < field_diameter / 2.0)
               && (std::abs(dec2 - telescope_dec) < field_diameter / 2.0)));
               
    return true;
}
    
} // namespace
