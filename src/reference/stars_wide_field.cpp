///----------------------------------------
///      @file stars_wide_field.cpp
///   @ingroup ASTAP++
///     @brief Implementation of wide-field star catalog loading.
///    @author Ported from Han Kleijn's ASTAP (MPL-2.0).
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "stars_wide_field.h"

#include "star_database.h"

#include <cstdint>
#include <fstream>
#include <ios>

///----------------------------------------
namespace astap::reference {
///----------------------------------------

/// MARK: - Module state

std::vector<float> wide_field_stars;
std::string wide_database;

/// MARK: - Catalog loading

///----------------------------------------
/// @brief Load the wide-field catalog from disk into @c wide_field_stars.
/// @details File layout (little-endian, packed IEEE-754 singles):
///          - int32 count
///          - count * { float mag_x10; float ra_rad; float dec_rad; }
///          Stars are sorted bright-to-faint.
/// @return @c true on success, @c false on any I/O failure.
///----------------------------------------

[[nodiscard]] bool read_stars_wide_field() {
    // Build the catalog file path
    const auto path = database_path / (name_database + "_0101.001");
    
    // Open the binary file
    auto file = std::ifstream(path, std::ios::binary);
    if (!file) {
        return false;
    }
    
    // Read the star count
    auto count = std::int32_t{0};
    if (!file.read(reinterpret_cast<char*>(&count), sizeof(count))) {
        return false;
    }
    if (count <= 0) {
        return false;
    }
    
    // Allocate storage for mag/ra/dec triples
    constexpr auto floats_per_star = std::size_t{3};
    const auto total_floats = static_cast<std::size_t>(count) * floats_per_star;
    wide_field_stars.assign(total_floats, 0.0f);
    
    // Read the star data in one block
    const auto bytes = static_cast<std::streamsize>(total_floats * sizeof(float));
    if (!file.read(reinterpret_cast<char*>(wide_field_stars.data()), bytes)) {
        wide_field_stars.clear();
        return false;
    }
    
    // Record which database is now loaded
    wide_database = name_database;
    return true;
}
    
} // namespace astap::reference
