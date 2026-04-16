///----------------------------------------
///      @file util.h
///   @ingroup ASTAP++
///     @brief Miscellaneous utility helpers: sorting, colour-space conversions,
///            locale-independent formatting, filename classifiers, and light
///            obfuscation.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP); MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>

///----------------------------------------
namespace astap::core {
///----------------------------------------

/// MARK: Math

///----------------------------------------
/// @brief In-place quicksort over elements with indices in [lo, hi] (inclusive).
/// @param a   Span of doubles to sort.
/// @param lo  First index (inclusive).
/// @param hi  Last index (inclusive).
///----------------------------------------

void quicksort(std::span<double> a, int lo, int hi);

///----------------------------------------
/// @brief Median of the first @p leng doubles.
/// @details Sorts a working copy.  For odd lengths >= 5 returns the average of
///          the three central values.  For even lengths returns the average of
///          the two central values.  Returns NaN when @p leng == 0.
/// @param list Source data (not modified).
/// @param leng Number of elements to consider.
/// @return Median value.
///----------------------------------------

[[nodiscard]] double smedian(std::span<const double> list, int leng);

///----------------------------------------
/// @brief Median absolute deviation and median, computed without modifying @p list.
/// @param      list   Source data.
/// @param      leng   Number of elements to consider.
/// @param[out] mad    Median absolute deviation.
/// @param[out] median Median value.
///----------------------------------------

void mad_median(std::span<const double> list, int leng,
                double& mad, double& median);
                
///----------------------------------------
/// @brief Floating-point modulo helper for wrapping angles into [0, range).
/// @param x     Value to wrap.
/// @param range Modulus (must be positive).
/// @return Wrapped value in [0, range).
///----------------------------------------

[[nodiscard]] double fnmodulo(double x, double range) noexcept;

/// MARK: Colour

///----------------------------------------
/// @brief Convert HSV to RGB.
/// @param h         Hue in [0, 360).
/// @param s         Saturation in [0, 1].
/// @param v         Value in [0, 1].
/// @param[out] r    Red in [0, 1].
/// @param[out] g    Green in [0, 1].
/// @param[out] b    Blue in [0, 1].
///----------------------------------------

void hsv_to_rgb(float h, float s, float v, float& r, float& g, float& b) noexcept;

///----------------------------------------
/// @brief Convert RGB to HSV.
/// @param r         Red in [0, 1].
/// @param g         Green in [0, 1].
/// @param b         Blue in [0, 1].
/// @param[out] h    Hue in [0, 360).
/// @param[out] s    Saturation in [0, 1].
/// @param[out] v    Value in [0, 1].
///----------------------------------------

void rgb_to_hsv(float r, float g, float b, float& h, float& s, float& v) noexcept;

///----------------------------------------
/// @brief Hue in [0, 360] for given RGB values, rounded to nearest integer.
/// @param r Red in [0, 1].
/// @param g Green in [0, 1].
/// @param b Blue in [0, 1].
/// @return Hue in degrees.
///----------------------------------------

[[nodiscard]] int rgb_to_h(float r, float g, float b) noexcept;

///----------------------------------------
/// @brief Decompose a packed COLORREF (0x00BBGGRR) into R/G/B byte components.
/// @param colour Packed colour value.
/// @return Array of {R, G, B} bytes.
///----------------------------------------

[[nodiscard]] std::array<std::uint8_t, 3> intensity_rgb(std::uint32_t colour) noexcept;

///----------------------------------------
/// @brief Estimate colour temperature in Kelvin from red and blue intensities.
/// @details Returns a formatted string like "5500K", or "" when the input is
///          saturated/dark or the ratio falls outside the valid range.
/// @param red  Red channel intensity.
/// @param blue Blue channel intensity.
/// @return Formatted temperature string, or empty on failure.
///----------------------------------------

[[nodiscard]] std::string rgb_kelvin(float red, float blue);

/// MARK: Number / string formatting

///----------------------------------------
/// @brief Format double with 8 decimal places.
///----------------------------------------

[[nodiscard]] std::string floattostr8(double x);

///----------------------------------------
/// @brief Format double with 6 decimal places.
///----------------------------------------

[[nodiscard]] std::string floattostr6(double x);

///----------------------------------------
/// @brief Format double with 4 decimal places.
///----------------------------------------

[[nodiscard]] std::string floattostr4(double x);

///----------------------------------------
/// @brief Format double with 2 decimal places.
///----------------------------------------

[[nodiscard]] std::string floattostr2(double x);

///----------------------------------------
/// @brief Format double in scientific notation.
///----------------------------------------

[[nodiscard]] std::string floattostrE(double x);

///----------------------------------------
/// @brief Integer formatted right-aligned in a field of width 5.
///----------------------------------------

[[nodiscard]] std::string inttostr5(int x);

///----------------------------------------
/// @brief Fault-tolerant string to integer; returns @p default_value on failure.
/// @param s             Input string.
/// @param default_value Fallback value.
/// @return Parsed integer or @p default_value.
///----------------------------------------

[[nodiscard]] int strtoint2(std::string_view s, int default_value);

///----------------------------------------
/// @brief Fault-tolerant string to double, accepts both ',' and '.' separators.
/// @param s Input string.
/// @return Parsed value, or 0 on failure.
///----------------------------------------

[[nodiscard]] double strtofloat2(std::string_view s);

///----------------------------------------
/// @brief Fault-tolerant string to double, dot separator only.
/// @param s Input string.
/// @return Parsed value, or 0 on failure.
///----------------------------------------

[[nodiscard]] double strtofloat1(std::string_view s);

///----------------------------------------
/// @brief Write @p inp right-aligned into columns 11..30 of @p target_line.
/// @param[in,out] target_line String to modify (grown if necessary).
/// @param         inp         Value to format.
///----------------------------------------

void addstring(std::string& target_line, double inp);

/// MARK: Filename helpers

///----------------------------------------
/// @brief Test whether path has a FITS extension (.fit/.fits/.fts/.wcs).
///----------------------------------------

[[nodiscard]] bool fits_file_name(std::string_view path);

///----------------------------------------
/// @brief Test whether path has a FITS or TIFF extension.
///----------------------------------------

[[nodiscard]] bool fits_tiff_file_name(std::string_view path);

///----------------------------------------
/// @brief Test whether path has a TIFF extension (.tif/.tiff).
///----------------------------------------

[[nodiscard]] bool tiff_file_name(std::string_view path);

///----------------------------------------
/// @brief Test whether the extension denotes a raw camera format.
///----------------------------------------

[[nodiscard]] bool check_raw_file_extension(std::string_view ext);

///----------------------------------------
/// @brief Test whether path has any readable image extension.
///----------------------------------------

[[nodiscard]] bool image_file_name(std::string_view path);

///----------------------------------------
/// @brief Recover an exposure time (seconds) from an embedded filename token.
/// @details Looks for "<digits>SEC" or "<digits>S_" patterns.
/// @param filename Input filename.
/// @return Exposure in seconds, or 0 when no match is found.
///----------------------------------------

[[nodiscard]] int extract_exposure_from_filename(std::string_view filename);

///----------------------------------------
/// @brief Recover a sensor temperature (degrees C) from a filename token.
/// @param filename Input filename.
/// @return Temperature in degrees C, or 999 (sentinel for unknown).
///----------------------------------------

[[nodiscard]] int extract_temperature_from_filename(std::string_view filename);

///----------------------------------------
/// @brief Recover an object designation (e.g. "M31", "NGC7000") from a filename.
/// @param filename Input filename.
/// @return Object name, or "" when no recognised prefix is found.
///----------------------------------------

[[nodiscard]] std::string extract_objectname_from_filename(std::string_view filename);

/// MARK: Light obfuscation

///----------------------------------------
/// @brief Light obfuscation for stored credentials (NOT cryptographic).
/// @param inp Plaintext input.
/// @return Obfuscated string.
///----------------------------------------

[[nodiscard]] std::string encrypt(std::string_view inp);

///----------------------------------------
/// @brief Reverse the obfuscation applied by encrypt().
/// @param inp Obfuscated input.
/// @return Plaintext string, or "" on invalid input.
///----------------------------------------

[[nodiscard]] std::string decrypt(std::string_view inp);
                
} // namespace
