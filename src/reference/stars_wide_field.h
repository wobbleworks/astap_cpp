///----------------------------------------
///      @file stars_wide_field.h
///   @ingroup ASTAP++
///     @brief Runtime access to the wide-field star catalog (e.g. W08).
///   @details The catalog is stored on disk as
///            @c \<database_path\>/\<name_database\>_0101.001 and loaded into a
///            flat float buffer (mag, ra, dec triples) sorted bright-to-faint.
///    @author Ported from Han Kleijn's ASTAP (MPL-2.0).
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <string>
#include <vector>

///----------------------------------------
namespace astap::reference {
///----------------------------------------

/// MARK: - Wide-field star data

/// @brief Runtime flat buffer of stars.
/// @details Layout: mag1, ra1, dec1, mag2, ra2, dec2, ...
///          Magnitude is stored as mag * 10; RA and Dec are in radians.
extern std::vector<float> wide_field_stars;

/// @brief Name of the database currently loaded into @c wide_field_stars
///        (empty if nothing loaded).
extern std::string wide_database;

///----------------------------------------
/// @brief Load the wide-field catalog from disk into @c wide_field_stars.
/// @details On success sets @c wide_database to @c name_database.
/// @return @c true on success, @c false on any I/O failure (open, short read).
///----------------------------------------

[[nodiscard]] bool read_stars_wide_field();
	
} // namespace astap::reference
