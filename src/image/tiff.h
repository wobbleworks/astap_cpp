///----------------------------------------
///      @file tiff.h
///   @ingroup ASTAP++
///     @brief Uncompressed TIFF writer for 16/32-bit grayscale and 48/96-bit RGB.
///   @details Based originally on 8-bit routines from bit2tiff.pas by Wolfgang Krug.
///            Heavily modified for 16-bit integer and 32-bit float gray and colour.
///    @author Ported from Han Kleijn's ASTAP. MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <filesystem>
#include <string_view>

#include "../types.h"

///----------------------------------------
namespace astap::image {
///----------------------------------------

using astap::ImageArray;

/// MARK: - 16-bit grayscale

///----------------------------------------
/// @brief Save image as a 16-bit grayscale TIFF.
/// @details Not used in ASTAP; preserved for completeness.
///  @param img        Source image array (single channel used).
///  @param filename   Output path; extension is forced to .tif.
///  @param description TIFF ImageDescription tag content.
///  @param flip_h     Mirror horizontally.
///  @param flip_v     Mirror vertically.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool save_tiff_16(const ImageArray& img,
                                const std::filesystem::path& filename,
                                std::string_view description,
                                bool flip_h,
                                bool flip_v);
                                
/// MARK: - 32-bit float grayscale

///----------------------------------------
/// @brief Save image as a 32-bit float grayscale TIFF.
///  @param img        Source image array (single channel used).
///  @param filename   Output path; extension is forced to .tif.
///  @param description TIFF ImageDescription tag content.
///  @param flip_h     Mirror horizontally.
///  @param flip_v     Mirror vertically.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool save_tiff_32(const ImageArray& img,
                                const std::filesystem::path& filename,
                                std::string_view description,
                                bool flip_h,
                                bool flip_v);
                                
/// MARK: - 48-bit (3x16) colour

///----------------------------------------
/// @brief Save image as a 48-bit (3x16) colour TIFF.
/// @details Not used in ASTAP; preserved for completeness.
///  @param img        Source image array (three channels required).
///  @param filename   Output path; extension is forced to .tif.
///  @param description TIFF ImageDescription tag content.
///  @param flip_h     Mirror horizontally.
///  @param flip_v     Mirror vertically.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool save_tiff_48(const ImageArray& img,
                                const std::filesystem::path& filename,
                                std::string_view description,
                                bool flip_h,
                                bool flip_v);
                                
/// MARK: - 96-bit (3x32 float) colour

///----------------------------------------
/// @brief Save image as a 96-bit (3x32 float) colour TIFF.
///  @param img        Source image array (three channels required).
///  @param filename   Output path; extension is forced to .tif.
///  @param description TIFF ImageDescription tag content.
///  @param flip_h     Mirror horizontally.
///  @param flip_v     Mirror vertically.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool save_tiff_96(const ImageArray& img,
                                const std::filesystem::path& filename,
                                std::string_view description,
                                bool flip_h,
                                bool flip_v);
                                
} // namespace
