///----------------------------------------
///      @file xisf.h
///   @ingroup ASTAP++
///     @brief XISF image loader for uncompressed PixInsight files.
///   @details Reads an uncompressed XISF file (PixInsight's XML-based FITS-like
///            astronomical image format) with 8, 16, 32, -32, or -64 bit pixel
///            data. Parses the embedded XML header, extracts pixel data, and
///            populates a FITS-style Header struct.
///
///            XISF format reference:
///              http://pixinsight.com/doc/docs/XISF-1.0-spec/XISF-1.0-spec.html
///    @author Ported from Han Kleijn's ASTAP by John Stephen.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <filesystem>
#include <string>
#include <vector>

#include "../types.h"

///----------------------------------------
namespace astap::image {
///----------------------------------------

using astap::ImageArray;
using astap::Header;

///----------------------------------------
///   @brief Load an uncompressed XISF file into a FITS-style header and pixel array.
/// @details On success fills @p head with FITS-style metadata, populates
///          @p img_loaded2 with pixel data (scaled to the 0..65535 convention
///          used elsewhere in ASTAP), and appends status/header lines to @p log.
///          On failure returns false and @p head.naxis is set to 0.
///   @param filen Path to the XISF file.
///   @param[out] head FITS-style header to populate.
///   @param[out] img_loaded2 Pixel data array to fill.
///   @param[out] log Status and header lines appended here.
///  @return True on success, false on any failure.
///----------------------------------------

[[nodiscard]] bool load_xisf(const std::filesystem::path& filen,
                             Header&                      head,
                             ImageArray&                  img_loaded2,
                             std::vector<std::string>&    log);
                             
} // namespace

