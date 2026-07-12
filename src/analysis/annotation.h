#pragma once

/// @file annotation.h
/// Image annotation utilities — ported from unit_annotation.pas.
///
/// Provides a built-in 5x9 bitmap font, text-to-image rendering,
/// coordinate transforms, and deep-sky object lookup.  GUI-dependent
/// procedures (plot_deepsky, plot_vsx_vsp, measure_distortion, etc.)
/// are intentionally omitted.

#include "../types.h"

#include <cmath>
#include <cstdint>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace astap::analysis {

// ---------------------------------------------------------------------------
// VariableInfo — photometry-tab record (Pascal: tvariable_list)
// ---------------------------------------------------------------------------

/// Per-star metadata extracted during annotation for the photometry tab.
struct VariableInfo {
	double ra{};
	double dec{};
	std::string abbr;
	int source{};  // 0=local, 1=VSX, 2=VSP, 3=not in AAVSO
	int index{};
};

// ---------------------------------------------------------------------------
// 5x9 bitmap font (ASCII 33 '!' through 126 '~')
// ---------------------------------------------------------------------------

/// Built-in 5-wide, 9-tall bitmap font covering printable ASCII 33..126.
/// Index 0 corresponds to '!' (ASCII 33).  Each element is 0 or 1.
///
/// Pascal: packed array[33..126, 0..8, 0..4] of byte.
inline constexpr uint8_t kFont5x9[94][9][5] = {
	// --- declared in annotation.cpp via an explicit-specialisation trick ---
	// (kept here so the symbol is available to every TU that includes this
	//  header; the full 94-character table is defined below.)

	// [0] '!'
	{{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,0,0,0},{0,0,1,0,0}},
	// [1] '"'
	{{0,1,0,1,0},{0,1,0,1,0},{0,1,0,1,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
	// [2] '#'
	{{0,1,0,1,0},{0,1,0,1,0},{0,1,0,1,0},{1,1,1,1,1},{0,1,0,1,0},{1,1,1,1,1},{0,1,0,1,0},{0,1,0,1,0},{0,1,0,1,0}},
	// [3] '$'
	{{0,0,1,0,0},{0,0,1,0,0},{0,1,1,1,1},{1,0,1,0,0},{0,1,1,1,0},{0,0,1,0,1},{1,1,1,1,0},{0,0,1,0,0},{0,0,1,0,0}},
	// [4] '%'
	{{1,1,1,0,0},{1,0,1,0,0},{1,1,1,0,1},{0,0,0,1,0},{0,0,1,0,0},{0,1,0,0,0},{1,0,1,1,1},{0,0,1,0,1},{0,0,1,1,1}},
	// [5] '&'
	{{0,0,0,0,0},{0,1,1,0,0},{1,0,0,1,0},{1,0,0,1,0},{0,1,1,0,0},{1,0,1,0,0},{1,0,0,1,0},{1,0,0,1,1},{0,1,1,0,0}},
	// [6] '\''
	{{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
	// [7] '('
	{{0,0,0,1,0},{0,0,1,0,0},{0,0,1,0,0},{0,1,0,0,0},{0,1,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,0,1,0}},
	// [8] ')'
	{{0,1,0,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,1,0},{0,0,0,1,0},{0,0,1,0,0},{0,0,1,0,0},{0,1,0,0,0}},
	// [9] '*'
	{{0,0,0,0,0},{0,0,1,0,0},{1,0,1,0,1},{1,1,1,1,1},{0,1,1,1,0},{1,1,1,1,1},{1,0,1,0,1},{0,0,1,0,0},{0,0,0,0,0}},
	// [10] '+'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,1,0,0},{0,0,1,0,0},{1,1,1,1,1},{0,0,1,0,0},{0,0,1,0,0},{0,0,0,0,0},{0,0,0,0,0}},
	// [11] ','
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,1,0,0,0}},
	// [12] '-'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{1,1,1,1,1},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
	// [13] '.'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,1,1,0,0},{0,1,1,0,0}},
	// [14] '/'
	{{0,0,0,0,0},{0,0,0,1,0},{0,0,0,1,0},{0,0,1,0,0},{0,0,1,0,0},{0,1,0,0,0},{0,1,0,0,0},{1,0,0,0,0},{1,0,0,0,0}},
	// [15] '0'
	{{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{0,1,1,1,0}},
	// [16] '1'
	{{0,0,1,0,0},{0,1,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,1,1,1,0}},
	// [17] '2'
	{{0,1,1,1,0},{1,0,0,0,1},{0,0,0,0,1},{0,0,0,0,1},{0,0,0,1,0},{0,0,1,0,0},{0,1,0,0,0},{1,0,0,0,0},{1,1,1,1,1}},
	// [18] '3'
	{{1,1,1,1,0},{0,0,0,0,1},{0,0,0,0,1},{0,0,0,0,1},{0,1,1,1,1},{0,0,0,0,1},{0,0,0,0,1},{0,0,0,0,1},{1,1,1,1,0}},
	// [19] '4'
	{{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,1,1,1,1},{0,0,0,0,1},{0,0,0,0,1},{0,0,0,0,1},{0,0,0,0,1}},
	// [20] '5'
	{{1,1,1,1,1},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,1,1,1,0},{0,0,0,0,1},{0,0,0,0,1},{1,0,0,0,1},{0,1,1,1,0}},
	// [21] '6'
	{{0,1,1,1,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{0,1,1,1,0}},
	// [22] '7'
	{{1,1,1,1,1},{0,0,0,0,1},{0,0,0,0,1},{0,0,0,1,0},{0,0,0,1,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0}},
	// [23] '8'
	{{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{0,1,1,1,0}},
	// [24] '9'
	{{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{0,1,1,1,1},{0,0,0,0,1},{0,0,0,0,1},{0,0,0,1,0},{0,1,1,0,0}},
	// [25] ':'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,1,1,0},{0,0,1,1,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,1,1,0},{0,0,1,1,0}},
	// [26] ';'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,1,1,0},{0,0,1,1,0},{0,0,0,0,0},{0,0,1,1,0},{0,0,1,1,0},{0,0,0,1,0},{0,1,1,1,0}},
	// [27] '<'
	{{0,0,0,0,1},{0,0,0,1,0},{0,0,1,0,0},{0,1,0,0,0},{1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1}},
	// [28] '='
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{1,1,1,1,1},{0,0,0,0,0},{0,0,0,0,0},{1,1,1,1,1},{0,0,0,0,0},{0,0,0,0,0}},
	// [29] '>'
	{{1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1},{0,0,0,1,0},{0,0,1,0,0},{0,1,0,0,0},{1,0,0,0,0}},
	// [30] '?'
	{{1,1,1,1,0},{1,0,0,0,1},{0,0,0,0,1},{0,0,0,1,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,0,0,0},{0,0,1,0,0},{0,0,1,0,0}},
	// [31] '@'
	{{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,1,1,1},{1,0,1,0,1},{1,0,1,1,1},{1,0,0,0,0},{1,0,0,0,1},{0,1,1,1,0}},
	// [32] 'A'
	{{0,0,1,0,0},{0,1,0,1,0},{0,1,0,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,1,1,1,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1}},
	// [33] 'B'
	{{1,1,1,1,0},{0,1,0,0,1},{0,1,0,0,1},{0,1,0,0,1},{0,1,1,1,0},{0,1,0,0,1},{0,1,0,0,1},{0,1,0,0,1},{1,1,1,1,0}},
	// [34] 'C'
	{{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,1},{0,1,1,1,0}},
	// [35] 'D'
	{{1,1,1,1,0},{0,1,0,0,1},{0,1,0,0,1},{0,1,0,0,1},{0,1,0,0,1},{0,1,0,0,1},{0,1,0,0,1},{0,1,0,0,1},{1,1,1,1,0}},
	// [36] 'E'
	{{1,1,1,1,1},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,1,1,1,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,1,1,1,1}},
	// [37] 'F'
	{{1,1,1,1,1},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,1,1,1,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0}},
	// [38] 'G'
	{{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,1,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,1,1},{0,1,1,0,1}},
	// [39] 'H'
	{{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,1,1,1,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1}},
	// [40] 'I'
	{{1,1,1,1,1},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{1,1,1,1,1}},
	// [41] 'J'
	{{0,0,0,1,1},{0,0,0,0,1},{0,0,0,0,1},{0,0,0,0,1},{0,0,0,0,1},{0,0,0,0,1},{0,0,0,0,1},{1,0,0,0,1},{0,1,1,1,0}},
	// [42] 'K'
	{{1,0,0,0,1},{1,0,0,0,1},{1,0,0,1,0},{1,0,1,0,0},{1,1,0,0,0},{1,0,1,0,0},{1,0,0,1,0},{1,0,0,0,1},{1,0,0,0,1}},
	// [43] 'L'
	{{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,1,1,1,1}},
	// [44] 'M'
	{{1,0,0,0,1},{1,1,0,1,1},{1,1,0,1,1},{1,0,1,0,1},{1,0,1,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1}},
	// [45] 'N'
	{{1,0,0,0,1},{1,1,0,0,1},{1,1,0,0,1},{1,0,1,0,1},{1,0,1,0,1},{1,0,0,1,1},{1,0,0,1,1},{1,0,0,0,1},{1,0,0,0,1}},
	// [46] 'O'
	{{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{0,1,1,1,0}},
	// [47] 'P'
	{{1,1,1,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,1,1,1,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0}},
	// [48] 'Q'
	{{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,1,0,1},{0,1,1,1,0},{0,0,0,1,0},{0,0,0,0,1}},
	// [49] 'R'
	{{1,1,1,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,1,1,1,0},{1,1,0,0,0},{1,0,1,0,0},{1,0,0,1,0},{1,0,0,0,1}},
	// [50] 'S'
	{{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1},{1,0,0,0,1},{0,1,1,1,0}},
	// [51] 'T'
	{{1,1,1,1,1},{1,0,1,0,1},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,1,1,1,0}},
	// [52] 'U'
	{{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{0,1,1,1,0}},
	// [53] 'V'
	{{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{0,1,0,1,0},{0,1,0,1,0},{0,0,1,0,0}},
	// [54] 'W'
	{{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,1,0,1},{1,0,1,0,1},{1,1,0,1,1},{1,1,0,1,1},{1,0,0,0,1}},
	// [55] 'X'
	{{1,0,0,0,1},{1,0,0,0,1},{0,1,0,1,0},{0,1,0,1,0},{0,0,1,0,0},{0,1,0,1,0},{0,1,0,1,0},{1,0,0,0,1},{1,0,0,0,1}},
	// [56] 'Y'
	{{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{0,1,0,1,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0}},
	// [57] 'Z'
	{{1,1,1,1,1},{0,0,0,0,1},{0,0,0,0,1},{0,0,0,1,0},{0,0,1,0,0},{0,1,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,1,1,1,1}},
	// [58] '['
	{{0,1,1,1,1},{0,1,0,0,0},{0,1,0,0,0},{0,1,0,0,0},{0,1,0,0,0},{0,1,0,0,0},{0,1,0,0,0},{0,1,0,0,0},{0,1,1,1,1}},
	// [59] '\\'
	{{0,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{0,1,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,1,0}},
	// [60] ']'
	{{1,1,1,1,0},{0,0,0,1,0},{0,0,0,1,0},{0,0,0,1,0},{0,0,0,1,0},{0,0,0,1,0},{0,0,0,1,0},{0,0,0,1,0},{1,1,1,1,0}},
	// [61] '^'
	{{0,0,1,0,0},{0,1,0,1,0},{1,0,0,0,1},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
	// [62] '_'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{1,1,1,1,1}},
	// [63] '`'
	{{0,0,1,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
	// [64] 'a'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,1,1,1,0},{0,0,0,0,1},{0,1,1,1,1},{1,0,0,0,1},{1,0,0,1,1},{0,1,1,0,1}},
	// [65] 'b'
	{{0,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,1,1,0},{1,1,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,1,1,1,0}},
	// [66] 'c'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,1,1,1,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,1},{0,1,1,1,0}},
	// [67] 'd'
	{{0,0,0,0,0},{0,0,0,0,1},{0,0,0,0,1},{0,0,0,0,1},{1,1,1,1,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,1,1},{0,1,1,0,1}},
	// [68] 'e'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,1,1,1,0},{1,0,0,0,0},{0,1,1,1,1}},
	// [69] 'f'
	{{0,0,1,1,0},{0,1,0,0,1},{0,1,0,0,0},{0,1,0,0,0},{1,1,1,1,0},{0,1,0,0,0},{0,1,0,0,0},{0,1,0,0,0},{0,1,0,0,0}},
	// [70] 'g'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,1,1,1,1},{1,0,0,0,1},{1,0,0,1,1},{0,1,1,0,1},{0,0,0,0,1},{1,1,1,1,0}},
	// [71] 'h'
	{{0,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,1,1,0},{1,1,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1}},
	// [72] 'i'
	{{0,0,0,0,0},{0,1,1,0,0},{0,0,0,0,0},{0,1,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{1,1,1,1,1}},
	// [73] 'j'
	{{0,0,0,0,0},{0,0,0,1,1},{0,0,0,0,0},{0,0,0,1,1},{0,0,0,0,1},{0,0,0,0,1},{0,0,0,0,1},{1,0,0,0,1},{0,1,1,1,0}},
	// [74] 'k'
	{{0,0,0,0,0},{1,0,0,0,0},{1,0,0,0,0},{1,0,0,0,1},{1,0,0,1,0},{1,1,1,0,0},{1,0,1,0,0},{1,0,0,1,0},{1,0,0,0,1}},
	// [75] 'l'
	{{0,0,0,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,1,1}},
	// [76] 'm'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{1,1,0,1,0},{1,0,1,0,1},{1,0,1,0,1},{1,0,1,0,1},{1,0,0,0,1},{1,0,0,0,1}},
	// [77] 'n'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{1,0,1,1,0},{1,1,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1}},
	// [78] 'o'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,1,1,1,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{0,1,1,1,0}},
	// [79] 'p'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{1,0,1,1,0},{1,1,0,0,1},{1,0,0,0,1},{1,1,1,1,0},{1,0,0,0,0},{1,0,0,0,0}},
	// [80] 'q'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,1,1,1,1},{1,0,0,0,1},{1,0,0,1,1},{0,1,1,0,1},{0,0,0,0,1},{0,0,0,0,1}},
	// [81] 'r'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{1,1,0,1,1},{0,1,1,0,1},{0,1,0,0,1},{0,1,0,0,0},{0,1,0,0,0},{1,1,1,0,0}},
	// [82] 's'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,1,1,1,1},{1,0,0,0,0},{0,1,1,1,0},{0,0,0,0,1},{0,0,0,0,1},{1,1,1,1,0}},
	// [83] 't'
	{{0,0,0,0,0},{0,1,0,0,0},{0,1,0,0,0},{1,1,1,1,0},{0,1,0,0,0},{0,1,0,0,0},{0,1,0,0,0},{0,1,0,0,1},{0,0,1,1,0}},
	// [84] 'u'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,1,1},{0,1,1,0,1}},
	// [85] 'v'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,0,1},{0,1,0,1,0},{0,1,0,1,0},{0,0,1,0,0}},
	// [86] 'w'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,1,0,1},{1,0,1,0,1},{1,1,0,1,1},{1,0,0,0,1}},
	// [87] 'x'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{1,0,0,0,1},{0,1,0,1,0},{0,0,1,0,0},{0,1,0,1,0},{1,0,0,0,1},{1,0,0,0,1}},
	// [88] 'y'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{1,0,0,0,1},{1,0,0,0,1},{1,0,0,1,1},{0,1,1,0,1},{0,0,0,0,1},{1,1,1,1,0}},
	// [89] 'z'
	{{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{1,1,1,1,1},{0,0,0,1,0},{0,0,1,0,0},{0,1,0,0,0},{1,0,0,0,0},{1,1,1,1,1}},
	// [90] '{'
	{{0,0,1,1,0},{0,1,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{1,1,1,0,0},{0,0,1,0,0},{0,1,0,0,0},{0,1,0,0,0},{0,0,1,1,0}},
	// [91] '|'
	{{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0},{0,0,1,0,0}},
	// [92] '}'
	{{0,1,1,0,0},{0,0,0,1,0},{0,0,0,1,0},{0,0,1,0,0},{0,0,1,1,1},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,1,0},{0,1,1,0,0}},
	// [93] '~'
	{{0,0,0,0,0},{0,0,0,0,0},{0,1,0,0,1},{1,0,1,0,1},{1,0,1,0,1},{1,0,0,1,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
};

// ---------------------------------------------------------------------------
// Text rendering
// ---------------------------------------------------------------------------

/// Render a text string into an image buffer using the built-in 5x9 font.
///
/// Each glyph is scaled by @p size (1 = native 5x9, 2 = 10x18, etc.).
/// Characters outside ASCII 33..126 are treated as spaces (advance only).
///
/// @param text        String to render.
/// @param transparent If true, only "on" pixels are written (background
///                    is left unchanged).  If false, "off" pixels are
///                    written as 0.
/// @param colour      Pixel value written for "on" font pixels.
/// @param size        Integer scale factor (>=1).
/// @param x           Starting x-coordinate (column) in @p img.
/// @param y           Starting y-coordinate (row) in @p img.
/// @param img         Destination image buffer (channel 0 is written).
void annotation_to_array(std::string_view text, bool transparent, int colour,
                         int size, int x, int y, ImageArray& img);

// ---------------------------------------------------------------------------
// Coordinate utilities
// ---------------------------------------------------------------------------

/// Rotate a 2-D vector by @p rot radians (counter-clockwise, seen from
/// the y-axis).
///
/// @return {x2, y2} where x2 = x*sin(rot) + y*cos(rot),
///                        y2 = -x*cos(rot) + y*sin(rot).
[[nodiscard]]
inline std::pair<double, double> rotate(double rot, double x, double y)
{
	const double s = std::sin(rot);
	const double c = std::cos(rot);
	return {x * s + y * c, -x * c + y * s};
}

/// Convert equatorial (RA, Dec) to tangent-plane CCD pixel offsets
/// relative to the optical axis (ra0, dec0).
///
/// Uses the gnomonic (tangent-plane) projection.  The output xx, yy are
/// in pixel units, scaled by @p cdelt (arcsec per pixel).
///
/// @param ra0   Right ascension of the optical axis (radians).
/// @param dec0  Declination of the optical axis (radians).
/// @param ra    Right ascension of the target (radians).
/// @param dec   Declination of the target (radians).
/// @param cdelt CCD scale in arcseconds per pixel.
/// @return {xx, yy} pixel offsets from the optical axis.
[[nodiscard]]
inline std::pair<double, double> equatorial_standard(double ra0, double dec0,
                                                     double ra, double dec,
                                                     double cdelt)
{
	const double sin_dec0 = std::sin(dec0);
	const double cos_dec0 = std::cos(dec0);
	const double sin_dec  = std::sin(dec);
	const double cos_dec  = std::cos(dec);
	const double delta_ra = ra - ra0;
	const double sin_dra  = std::sin(delta_ra);
	const double cos_dra  = std::cos(delta_ra);

	const double dv = (cos_dec0 * cos_dec * cos_dra + sin_dec0 * sin_dec)
	                / (3600.0 * 180.0 / std::numbers::pi) * cdelt;

	const double xx = -cos_dec * sin_dra / dv;
	const double yy = -(sin_dec0 * cos_dec * cos_dra - cos_dec0 * sin_dec) / dv;
	return {xx, yy};
}

// ---------------------------------------------------------------------------
// Deep-sky object lookup
// ---------------------------------------------------------------------------

/// Result of a successful deep-sky object search.
struct DeepSkyMatch {
	std::string name;   ///< Canonical name(s), joined with '_'.
	double ra{};        ///< Right ascension (radians).
	double dec{};       ///< Declination (radians).
	double length{};    ///< Major axis (arcminutes).
	double width{};     ///< Minor axis (arcminutes).
	double pa{999.0};   ///< Position angle (degrees); 999 = unknown.
};

/// Search the deep-sky CSV database for an object by name.
///
/// The database file (deep_sky.csv) must have been loaded into
/// @p database_lines before calling this function.  Each line is a
/// comma-separated record: RA(encoded), Dec(encoded), names(/‑separated),
/// length, width, PA.
///
/// @param object_name  Name to search for (case-insensitive).
/// @param database_lines  Pre-loaded lines from deep_sky.csv.
/// @return A DeepSkyMatch on success, or std::nullopt if not found.
[[nodiscard]]
std::optional<DeepSkyMatch> find_object(std::string_view object_name,
                                        std::span<const std::string> database_lines);

// ---------------------------------------------------------------------------
// Deep-sky field iteration
// ---------------------------------------------------------------------------

/// A single deep-sky record returned by read_deepsky.
///
/// The three name fields mirror the Pascal naam2/naam3/naam4 globals: a record
/// may carry one to three slash-separated designations, kept separate so the
/// caller can build a display label however it likes.
struct DeepSkyObject {
	std::string name;    ///< Primary designation (Pascal naam2).
	std::string name2;   ///< First alternate designation (naam3); empty if none.
	std::string name3;   ///< Second alternate designation (naam4); empty if none.
	double ra{};         ///< Right ascension (radians).
	double dec{};        ///< Declination (radians).
	double length{};     ///< Major axis, catalog units (0.1 arcminute).
	double width{};      ///< Minor axis, catalog units (0.1 arcminute).
	double pa{999.0};    ///< Position angle (degrees); 999 = unknown.
};

/// Search scope for read_deepsky.
enum class DeepSkySearch {
	WithinField,   ///< Return only objects inside the FOV disc (Pascal 'S').
	FullDatabase,  ///< Return every record, ignoring the FOV gate (Pascal 'T').
};

/// Collect deep-sky records from a loaded catalog that fall within a field.
///
/// Mirrors unit_annotation.pas read_deepsky, but returns the full result set in
/// one pass instead of streaming a single record per call through a global
/// cursor.  The caller supplies the catalog lines — deep_sky.csv, hyperleda.csv
/// and the variable-star CSVs all share this format — and the first two lines
/// are treated as headers and skipped.
///
/// The field gate matches Pascal exactly: the RA offset is wrapped
/// (delta_ra > pi is folded to 2*pi - delta_ra) and the small-angle metric is
/// sqr(delta_ra * cos(dec0)) + sqr(dec - dec0), compared against sqr(fov).
///
/// @param database_lines  Pre-loaded catalog lines (headers first, then records).
/// @param mode            WithinField applies the FOV gate; FullDatabase skips it.
/// @param telescope_ra    Field-centre right ascension (radians).
/// @param telescope_dec   Field-centre declination (radians).
/// @param fov             Field radius (radians) used by the WithinField gate.
/// @return All matching records, in catalog order.
[[nodiscard]]
std::vector<DeepSkyObject> read_deepsky(std::span<const std::string> database_lines,
                                        DeepSkySearch mode,
                                        double telescope_ra, double telescope_dec,
                                        double fov);

// ---------------------------------------------------------------------------
// Deep-sky overlay (headless projection of a catalog onto a solved image)
// ---------------------------------------------------------------------------

/// Marker shape chosen for a plotted deep-sky object.
enum class DeepSkyMarker {
	Dots,    ///< len <= 2: too small to outline, drawn as four corner dots.
	Circle,  ///< len > 2 with unknown/degenerate PA: plain circle.
	Galaxy,  ///< len > 2 with a known PA: oriented oval.
};

/// One catalog object projected onto a solved image by plot_deepsky.
///
/// Positions are in image pixels (0-based), after the default FITS->image
/// vertical flip and any requested horizontal/vertical flips — i.e. exactly
/// where the marker and label would be drawn.
struct PlottedDeepSky {
	double ra{};             ///< Right ascension (radians).
	double dec{};            ///< Declination (radians).
	int x{};                 ///< Marker centre X.
	int y{};                 ///< Marker centre Y.
	DeepSkyMarker marker{};  ///< Chosen marker shape.
	double radius{};         ///< Object radius "len" in pixels.
	double orientation{};    ///< Oval orientation in degrees (Galaxy only).
	double axis_ratio{};     ///< Minor/major axis ratio (width/length).
	int pen_width{};         ///< Outline pen width, min(4, max(1, round(len/70))).
	std::string name;        ///< Composed designation ('/'-joined); empty if none.
	bool has_label{};        ///< True when the centre is on-image and named.
	int font_size{};         ///< Label point size, round(min(20, max(base, len/2))).
	bool is_reference{};     ///< naam2 begins with '0' (AAVSO comparison star).
};

/// Project a deep-sky catalog onto a solved image, headless.
///
/// Mirrors unit_annotation.pas plot_deepsky but performs no drawing: it returns
/// the placement of every catalog object that falls on (or just off) the image,
/// so a caller can render or diff the result.  The field centre and radius come
/// from @p head exactly as Pascal (image centre via pixel_to_celestial,
/// fov = 1.5 * half-diagonal); read_deepsky then supplies candidates, which are
/// projected with celestial_to_pixel, culled to a 25%-margin box, size-gated,
/// flipped, and assigned a marker shape.
///
/// The label-collision layout is intentionally not applied here (it depends on
/// font metrics the caller owns); feed the has_label entries to layout_labels.
///
/// @param head              Solved FITS header (needs naxis and cd1_1 set).
/// @param database_lines    Loaded catalog lines (deep_sky/hyperleda/variable).
/// @param variable_catalog  True for a variable-star catalog (Pascal db 3..5):
///                          skips the large-FOV size gate.
/// @param base_font_size    Minimum label point size (Pascal font_size arg).
/// @param formalism         Projection formalism forwarded to pixel_to_celestial.
/// @param flip_horizontal   Mirror X (Pascal flip_horizontal1 checkbox).
/// @param flip_vertical     Suppress the default vertical flip (flip_vertical1).
/// @return Placements in catalog order.
[[nodiscard]]
std::vector<PlottedDeepSky> plot_deepsky(const Header& head,
                                         std::span<const std::string> database_lines,
                                         bool variable_catalog = false,
                                         int base_font_size = 8,
                                         int formalism = 0,
                                         bool flip_horizontal = false,
                                         bool flip_vertical = false);

/// A label to be placed by layout_labels: its anchor and text extent.
struct LabelInput {
	int x{};   ///< Anchor X (marker centre).
	int y{};   ///< Anchor Y (marker centre).
	int tw{};  ///< Text width in pixels.
	int th{};  ///< Text height in pixels.
};

/// A resolved label box produced by layout_labels.
struct LabelBox {
	int x1{};          ///< Left.
	int y1{};          ///< Top (after any downward shift).
	int x2{};          ///< Right.
	int y2{};          ///< Bottom.
	bool connector{};  ///< True when the label moved and needs a leader line.
};

/// Resolve non-overlapping label positions, matching plot_deepsky exactly.
///
/// Labels are placed in order; each starts at its anchor, is clamped left off
/// the right image edge, then shifted down by th/3 until it clears every
/// previously placed box.  If it reaches the bottom edge it falls back to the
/// anchor (with no leader line).  A label whose top moved from its anchor is
/// flagged so the caller can draw a connector to the marker.
///
/// @param labels        Labels in draw order.
/// @param image_width   Image width in pixels.
/// @param image_height  Image height in pixels.
/// @return One box per input label, in the same order.
[[nodiscard]]
std::vector<LabelBox> layout_labels(std::span<const LabelInput> labels,
                                    int image_width, int image_height);

}  // namespace astap::analysis
