///----------------------------------------
///      @file image_io.h
///   @ingroup ASTAP++
///     @brief Non-FITS image I/O for ASTAP++ (PPM/PGM/PFM, TIFF/PNG/JPEG/BMP, RAW).
///   @details Covers PPM/PGM/PFM, TIFF/PNG/JPEG/BMP loading, RAW conversion via
///            DCRAW/LibRaw shell-out, the dispatcher that picks the right loader
///            by file extension, and the recent-files list helper.
///            FITS I/O lives in fits.h. Format-specific encoders live in
///            ../image/{tiff,xisf,...}.h.
///    @author Ported from Han Kleijn's ASTAP; MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "../types.h"

///----------------------------------------
namespace astap::core {
///----------------------------------------

using astap::Header;
using astap::ImageArray;

/// MARK: Image decoder backend

///----------------------------------------
/// @class IImageDecoder
/// @brief Abstract backend for raster and RAW image decoding.
/// @details ASTAP++ doesn't ship with a PNG / JPEG / TIFF / RAW decoder. The
///          library code operates against this interface; the host application
///          supplies a concrete implementation at startup (typically backed by
///          stb_image, libpng+libjpeg+libtiff, or platform-native APIs for
///          raster formats, and LibRaw/DCRAW for camera RAWs).
///
///          On success, fill @c out_pixels with decoded samples in
///          channel-separated, row-major order matching astap::ImageArray
///          ([channel][row][col], channel count in [1, 3]), set @c out_width,
///          @c out_height, @c out_channels, and return true.
///
///          Populate @c out_bits_per_sample with 8, 16, or 32; the loader uses
///          this to scale samples to the 0..65535 convention.
///
///          Return false on unsupported format, corrupt file, or I/O error.
///          The library will log to @c memo and continue with a synthesised
///          minimal header.
///----------------------------------------

struct IImageDecoder {
    struct DecodedImage {
        int               width             = 0;
        int               height            = 0;
        int               channels          = 0;  // 1 (mono) or 3 (RGB)
        int               bits_per_sample   = 8;  // 8, 16, or 32 (float)
        ImageArray        pixels;                 // [channel][row][col]
    };
    
    ///----------------------------------------
    /// @brief Decode a raster image (TIFF / PNG / JPEG / BMP).
    /// @param path File path to decode.
    /// @param ext_upper Upper-cased file extension (e.g. ".PNG", ".JPG").
    /// @param[out] out Decoded image output.
    /// @param[out] error_out Error description on failure.
    /// @return True on success.
    ///----------------------------------------
    
    virtual bool decode_raster(const std::filesystem::path& path,
                               std::string_view             ext_upper,
                               DecodedImage&                out,
                               std::string&                 error_out) = 0;
                               
    ///----------------------------------------
    /// @brief Decode a camera RAW file.
    /// @details Typically shells out to DCRAW or LibRaw. Callers pass
    ///          @p save_intermediate as true to produce a .pgm / .ppm alongside.
    /// @param path File path to decode.
    /// @param save_intermediate True to keep the intermediate file.
    /// @param[out] out Decoded image output.
    /// @param[out] error_out Error description on failure.
    /// @return True on success.
    ///----------------------------------------
    
    virtual bool decode_raw(const std::filesystem::path& path,
                            bool                         save_intermediate,
                            DecodedImage&                out,
                            std::string&                 error_out) = 0;
                            
    virtual ~IImageDecoder() = default;
};

/// MARK: Decoder management

///----------------------------------------
/// @brief Get the currently installed image decoder.
/// @return Reference to the active decoder (null decoder by default).
///----------------------------------------

[[nodiscard]] IImageDecoder& image_decoder() noexcept;

///----------------------------------------
/// @brief Install a decoder implementation.
/// @details Ownership is shared so the library doesn't force any particular
///          lifetime story on the host.
/// @param impl Shared pointer to the decoder; nullptr restores the null decoder.
///----------------------------------------

void set_image_decoder(std::shared_ptr<IImageDecoder> impl);

/// MARK: PPM / PGM / PFM

///----------------------------------------
/// @brief Load a PPM (P6), PGM (P5), or PFM (PF/Pf) image.
/// @details Faithful port of the hand-rolled scanner that also recognises a
///          custom DCRAW/LibRaw comment dialect (#EXPTIME=, #TIMESTAMP=, etc.).
///          Output pixels are scaled to the 0..65535 convention.
/// @param filen Path to the image file.
/// @param[out] head Parsed header fields.
/// @param[out] img_loaded2 Decoded image data.
/// @param[out] memo Synthesised FITS header cards.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool load_ppm_pgm_pfm(const std::filesystem::path& filen,
                                    Header&                      head,
                                    ImageArray&                  img_loaded2,
                                    std::vector<std::string>&    memo);
                                    
/// MARK: TIFF / PNG / JPEG / BMP

///----------------------------------------
/// @brief Load a TIFF, PNG, JPEG, or BMP file.
/// @details Dispatches on file extension. The underlying decode is delegated
///          to the installed @ref IImageDecoder.
/// @param filen Path to the image file.
/// @param light True for full header population.
/// @param[out] head Parsed header fields.
/// @param[out] img Decoded image data.
/// @param[out] memo Synthesised FITS header cards.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool load_tiff_png_jpeg(const std::filesystem::path& filen,
                                      bool                         light,
                                      Header&                      head,
                                      ImageArray&                  img,
                                      std::vector<std::string>&    memo);
                                      
/// MARK: TIFF save

///----------------------------------------
/// @brief Atomic 16-bit TIFF save via temporary file and rename.
/// @details Writes to <stem>.tmp, then deletes the destination and renames the
///          temp file into place, ensuring the destination is never lost if the
///          encode fails halfway.
/// @param img Source image data.
/// @param memo Header cards to embed as TIFF ImageDescription.
/// @param filen2 Destination file path.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool save_tiff16_secure(const ImageArray&            img,
                                      const std::vector<std::string>& memo,
                                      const std::filesystem::path& filen2);
                                      
/// MARK: RAW conversion

///----------------------------------------
/// @brief Convert a RAW camera file to FITS via DCRAW or LibRaw.
/// @details On success @p filename3 is updated to the new file path (".fits").
///          @p loadfile controls whether @p img / @p head are populated;
///          @p savefile controls whether the intermediate is kept.
/// @param loadfile True to populate @p img and @p head.
/// @param savefile True to keep the intermediate file.
/// @param[in,out] filename3 File path; updated to the output path on success.
/// @param[out] head Parsed header fields.
/// @param[out] img Decoded image data.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool convert_raw(bool                   loadfile,
                               bool                   savefile,
                               std::string&           filename3,
                               Header&                head,
                               ImageArray&            img);
                               
///----------------------------------------
/// @brief Convert any supported non-FITS file to FITS in place.
/// @details Writes <stem>.fits alongside the source and updates @p filen to
///          the new path.
/// @param[in,out] filen File path; updated to the .fits path on success.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool convert_to_fits(std::string& filen);

/// MARK: Master dispatcher

///----------------------------------------
/// @brief Load any supported image format, dispatching by file extension.
/// @details Picks FITS, XISF, TIFF/PNG/JPEG, PPM/PGM/PFM, RAW, or .fz
///          (CFITSIO) based on file extension and loads accordingly.
/// @param[in,out] filename2 File path; may be updated for format conversions.
/// @param[out] img Decoded image data.
/// @param[out] head Parsed header fields.
/// @param[out] memo Header card vector.
/// @param re_center True to centre the image after load (GUI hint).
/// @param plot True to re-render the image after load (GUI hint).
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool load_image(std::string&                 filename2,
                              ImageArray&                  img,
                              Header&                      head,
                              std::vector<std::string>&    memo,
                              bool                         re_center,
                              bool                         plot);
                              
/// MARK: Recent files

///----------------------------------------
/// @brief Append a path to the recent-files list.
/// @details Hoists an existing entry to position 0 and caps the list at 8 entries.
/// @param f File path to add.
/// @param[in,out] recent_files The recent-files list.
///----------------------------------------

void add_recent_file(std::string_view              f,
                     std::vector<std::string>&     recent_files);
    
} // namespace
