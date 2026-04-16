///----------------------------------------
///      @file yuv4mpeg2.h
///   @ingroup ASTAP++
///     @brief Stateful streaming writer for YUV4MPEG2 (y4m) uncompressed video.
///   @details Writes an ASCII stream header followed by concatenated YUV 4:4:4
///            (or mono) frames.  See https://wiki.multimedia.cx/index.php/YUV4MPEG2
///    @author Ported from Han Kleijn's ASTAP. MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string_view>

#include "../types.h"

///----------------------------------------
namespace astap::image {
///----------------------------------------

using astap::PixelSource;

///----------------------------------------
/// @class Y4mWriter
///  @brief Streaming writer for a single y4m file.
///----------------------------------------

class Y4mWriter final {
public:
    Y4mWriter() = default;
    ~Y4mWriter();
    
    Y4mWriter(const Y4mWriter&) = delete;
    Y4mWriter& operator=(const Y4mWriter&) = delete;
    Y4mWriter(Y4mWriter&&) = default;
    Y4mWriter& operator=(Y4mWriter&&) = default;
    
    ///----------------------------------------
    ///   @brief Open/create the file and write the y4m stream header.
    /// @details @p framerate is the numerator as ASCII (denominator is fixed at 1),
    ///          e.g. "30" yields "F30:1".
    ///   @param path      Output file path.
    ///   @param framerate  Frame-rate numerator as ASCII text.
    ///   @param colour     True for 4:4:4 colour, false for mono.
    ///   @param width      Frame width in pixels.
    ///   @param height     Frame height in pixels.
    ///  @return False on I/O failure.
    ///----------------------------------------
    
    [[nodiscard]] bool open(const std::filesystem::path& path,
                            std::string_view framerate,
                            bool colour,
                            int width,
                            int height);
                            
    ///----------------------------------------
    ///   @brief Write a single frame preceded by the "FRAME\n" marker.
    /// @details Reads pixels from @p src over the rectangle [x, x+w) x [y, y+h).
    ///          When @p colour is true a 4:4:4 Y/U/V triple-plane frame is
    ///          emitted; otherwise only the Y (luma) plane is written.
    ///   @param colour True for colour, false for mono.
    ///   @param src    Pixel source providing packed RGB pixels.
    ///   @param x      Left edge of the source rectangle.
    ///   @param y      Top edge of the source rectangle.
    ///   @param w      Width of the source rectangle.
    ///   @param h      Height of the source rectangle.
    ///  @return False on I/O failure.
    ///----------------------------------------
    
    [[nodiscard]] bool write_frame(bool colour,
                                   const PixelSource& src,
                                   int x,
                                   int y,
                                   int w,
                                   int h);
                                   
    ///----------------------------------------
    /// @brief Close the underlying stream.  Safe to call multiple times.
    ///----------------------------------------
    
    void close() noexcept;
    
    ///----------------------------------------
    /// @brief Whether the output stream is currently open.
    ///----------------------------------------
    
    [[nodiscard]] bool is_open() const noexcept { return _out.is_open(); }
    
private:
    std::ofstream _out;
    int _w{};
    int _h{};
};
    
} // namespace
