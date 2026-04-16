///----------------------------------------
///      @file avi.h
///   @ingroup ASTAP++
///     @brief Streaming writer for uncompressed RGB24 AVI (Audio Video Interleave) files.
///   @details Hand-rolled RIFF/AVI container with a single uncompressed video stream.
///            Produces files that are bit-identical to the original implementation
///            (modulo the pixel source abstraction).
///    @author Ported from Han Kleijn's ASTAP (MPL-2.0) by John Stephen.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>

#include "../types.h"

///----------------------------------------
namespace astap::image {
///----------------------------------------

/// @brief Abstract pixel source (canonical definition in types.h).
using astap::PixelSource;

///----------------------------------------
/// MARK: AviWriter
///----------------------------------------

///----------------------------------------
///  @class AviWriter
///  @brief Stateful streaming writer for an uncompressed RGB24 AVI file.
/// @details Usage:
///          @code
///          astap::image::AviWriter w;
///          if (!w.open("out.avi", "25", nrframes, width, height)) return;
///          for (int i = 0; i < nrframes; ++i)
///              w.write_frame(src, 0, 0, width, height);
///          w.close(nrframes);
///          @endcode
///----------------------------------------

class AviWriter final {
public:
    AviWriter() = default;
    ~AviWriter();
    
    AviWriter(const AviWriter&)            = delete;
    AviWriter& operator=(const AviWriter&) = delete;
    AviWriter(AviWriter&&)                 = delete;
    AviWriter& operator=(AviWriter&&)      = delete;
    
    ///----------------------------------------
    ///   @brief Create file and write the RIFF / hdrl / strl / movi headers.
    /// @details @p frame_rate is a decimal string (e.g. "25", "29.97"); parsed
    ///          with std::stod.
    ///   @param path       Output file path.
    ///   @param frame_rate Frames per second as a decimal string.
    ///   @param nr_frames  Total number of frames that will be written.
    ///   @param width      Frame width in pixels.
    ///   @param height     Frame height in pixels.
    ///  @return @c false on failure (file open, bad frame_rate, etc.).
    ///----------------------------------------
    
    [[nodiscard]] bool open(const std::filesystem::path& path,
                            std::string_view             frame_rate,
                            int                          nr_frames,
                            int                          width,
                            int                          height);
                            
    ///----------------------------------------
    ///   @brief Append one frame to the AVI file.
    /// @details Reads a (w x h) window at (x, y) from @p src and emits a
    ///          '00db' chunk with BGR24 scanlines, bottom-up, padded per
    ///          line to a 4-byte multiple.
    ///   @param src Pixel source to read from.
    ///   @param x   Left edge of the source window.
    ///   @param y   Top edge of the source window.
    ///   @param w   Width of the source window in pixels.
    ///   @param h   Height of the source window in pixels.
    ///  @return @c false on write failure.
    ///----------------------------------------
    
    [[nodiscard]] bool write_frame(const PixelSource& src, int x, int y, int w, int h);
    
    ///----------------------------------------
    /// @brief Write the idx1 index and close the file.
    /// @param nr_frames Must match the value passed to open().
    ///----------------------------------------
    
    void close(int nr_frames);
    
    ///----------------------------------------
    /// @brief Whether the underlying file stream is currently open.
    ///----------------------------------------
    
    [[nodiscard]] bool is_open() const noexcept { return _file.is_open(); }
    
private:
    std::ofstream _file;
    int           _extra = 0;          ///< Zero-byte pad per scanline (0..3).
    std::uint32_t _sizeimage = 0;      ///< Bytes per frame including padding.
};

///----------------------------------------
/// MARK: Legacy free-function API
///----------------------------------------

///----------------------------------------
/// @brief Legacy shim — opens an AVI via a module-level singleton writer.
/// @deprecated Use astap::image::AviWriter instead.
///----------------------------------------

[[deprecated("Use astap::image::AviWriter")]]
bool write_avi_head(const std::filesystem::path& filen,
                    std::string_view             frame_rate,
                    int                          nrframes,
                    int                          w,
                    int                          h);
                    
///----------------------------------------
/// @brief Legacy shim — writes one frame via a module-level singleton writer.
/// @deprecated Use astap::image::AviWriter instead.
///----------------------------------------

[[deprecated("Use astap::image::AviWriter")]]
bool write_avi_frame(const PixelSource& src, int x, int y, int w, int h);

///----------------------------------------
/// @brief Legacy shim — closes the AVI via a module-level singleton writer.
/// @deprecated Use astap::image::AviWriter instead.
///----------------------------------------

[[deprecated("Use astap::image::AviWriter")]]
void close_the_avi(int nrframes);
    
} // namespace
