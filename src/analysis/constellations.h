///----------------------------------------
///      @file constellations.h
///   @ingroup analysis
///     @brief Constellation stick-figure line data, label positions, and IAU abbreviations.
///   @details Ported from unit_constellations.pas (version 2001-10-21).
///            Original copyright (C) 1997, 2021 by Han Kleijn, www.hnsky.org.
///            Licensed under the Mozilla Public License, v. 2.0.
///    @author Created by John Stephen on 4/15/26.
/// @copyright Copyright 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <array>
#include <cstdint>
#include <string_view>

///----------------------------------------
namespace astap::analysis {
///----------------------------------------

/// @brief Drawing-mode sentinel: start a new line strip.
static constexpr int16_t kDrawModeStart = -2;

/// @brief Drawing-mode sentinel: draw a line to this star.
static constexpr int16_t kDrawModeDraw = -1;

///----------------------------------------
/// @brief A single vertex in a constellation stick-figure, carrying its Bayer designation.
///----------------------------------------

struct ConstStar {
	int16_t dm;          ///< Drawing mode: -2 = start new strip, -1 = draw to this point.
	int16_t ra;          ///< Right ascension, RA * 1000 (range 0..24000).
	int16_t dec;         ///< Declination, Dec * 100 (range -9000..9000).
	const char* bay;     ///< Bayer letter (UTF-8, e.g. "a", "b") or empty string.
};

/// @brief Number of entries in the constellation stick-figure array.
static constexpr size_t kConstellationLength = 603;

/// @brief Number of IAU constellations.
static constexpr size_t kConstellationCount = 89;

/// @brief Constellation stick-figure vertices (603 entries, index 0..602).
extern const std::array<ConstStar, kConstellationLength> constellation;

/// @brief Label positions for each constellation: [i][0] = RA*1000, [i][1] = Dec*100.
extern const std::array<std::array<int16_t, 2>, kConstellationCount> constPos;

/// @brief Three-letter IAU abbreviations for each constellation (89 entries, index 0..88).
extern const std::array<std::string_view, kConstellationCount> constShortName;

} // namespace
