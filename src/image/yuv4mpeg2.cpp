///----------------------------------------
///      @file yuv4mpeg2.cpp
///   @ingroup ASTAP++
///     @brief Implementation of @ref Y4mWriter — YUV4MPEG2 streaming output.
///    @author Ported from Han Kleijn's ASTAP. MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "yuv4mpeg2.h"

#include <cstdint>
#include <format>
#include <ios>
#include <string>
#include <vector>

///----------------------------------------
namespace astap::image {
///----------------------------------------

/// MARK: - BT.601 colour-space helpers

namespace {

/// @brief BT.601 full-swing RGB-to-Y (luma) conversion.
[[nodiscard]] constexpr auto to_y(std::uint8_t r, std::uint8_t g, std::uint8_t b) noexcept
    -> std::uint8_t {
    // trunc(R*77/256 + G*150/256 + B*29/256)
    auto v = (r * 77 + g * 150 + b * 29) / 256;
    return static_cast<std::uint8_t>(v);
}

/// @brief BT.601 full-swing RGB-to-U (Cb) conversion.
[[nodiscard]] constexpr auto to_u(std::uint8_t r, std::uint8_t g, std::uint8_t b) noexcept
    -> std::uint8_t {
    auto num = -43 * static_cast<int>(r)
               + -84 * static_cast<int>(g)
               + 127 * static_cast<int>(b);
    // Integer division in C++ truncates toward zero (C++11+), matching trunc().
    auto v = num / 256 + 128;
    return static_cast<std::uint8_t>(v);
}

/// @brief BT.601 full-swing RGB-to-V (Cr) conversion.
[[nodiscard]] constexpr auto to_v(std::uint8_t r, std::uint8_t g, std::uint8_t b) noexcept
    -> std::uint8_t {
    auto num = 127 * static_cast<int>(r)
               + -106 * static_cast<int>(g)
               + -21 * static_cast<int>(b);
    auto v = num / 256 + 128;
    return static_cast<std::uint8_t>(v);
}

/// @brief Unpack a 0x00RRGGBB pixel into its red, green, and blue components.
constexpr void unpack_rgb(std::uint32_t pixel,
                          std::uint8_t& r,
                          std::uint8_t& g,
                          std::uint8_t& b) noexcept {
    r = static_cast<std::uint8_t>((pixel >> 16) & 0xFF);
    g = static_cast<std::uint8_t>((pixel >> 8) & 0xFF);
    b = static_cast<std::uint8_t>(pixel & 0xFF);
}
    
} // namespace

/// MARK: - Y4mWriter

Y4mWriter::~Y4mWriter() {
    close();
}

bool Y4mWriter::open(const std::filesystem::path& path,
                     std::string_view framerate,
                     bool colour,
                     int width,
                     int height) {
    // Close any previously open stream
    close();
    
    // Open the output file in binary mode
    _out.open(path, std::ios::out | std::ios::binary | std::ios::trunc);
    if (!_out.is_open()) {
        return false;
    }
    
    _w = width;
    _h = height;
    
    // Write the stream header:
    // 'YUV4MPEG2 W{w} H{h} F{fr}:1 Ip A0:0 C{444|mono}\n'
    auto header = std::format(
        "YUV4MPEG2 W{} H{} F{}:1 Ip A0:0 {}\n",
        width,
        height,
        framerate,
        colour ? "C444" : "Cmono");
        
    _out.write(header.data(), static_cast<std::streamsize>(header.size()));
    if (!_out.good()) {
        _out.close();
        return false;
    }
    return true;
}

bool Y4mWriter::write_frame(bool colour,
                            const PixelSource& src,
                            int x,
                            int y,
                            int w,
                            int h) {
    if (!_out.is_open()) {
        return false;
    }
    
    // Write per-frame marker
    constexpr auto frame_marker = "FRAME\n";
    _out.write(frame_marker, 6);
    if (!_out.good()) {
        return false;
    }
    
    // 4:4:4 frames: full Y plane, then U plane, then V plane.
    // For mono only the Y plane is written (matches the Cmono header).
    constexpr auto mono_steps = 0;
    constexpr auto colour_steps = 2;
    auto steps = colour ? colour_steps : mono_steps;
    
    auto row = std::vector<std::uint8_t>(static_cast<std::size_t>(w));
    
    for (auto k = 0; k <= steps; ++k) {
        for (auto yy = y; yy < y + h; ++yy) {
            for (auto xx = x; xx < x + w; ++xx) {
                auto r = std::uint8_t{};
                auto g = std::uint8_t{};
                auto b = std::uint8_t{};
                unpack_rgb(src.pixel(xx, yy), r, g, b);
                
                auto sample = std::uint8_t{};
                if (k == 0) {
                    sample = to_y(r, g, b);
                } else if (k == 1) {
                    sample = to_u(r, g, b);
                } else {
                    sample = to_v(r, g, b);
                }
                row[static_cast<std::size_t>(xx - x)] = sample;
            }
            _out.write(reinterpret_cast<const char*>(row.data()),
                       static_cast<std::streamsize>(row.size()));
            if (!_out.good()) {
                return false;
            }
        }
    }
    return true;
}

void Y4mWriter::close() noexcept {
    if (_out.is_open()) {
        _out.close();
    }
}
    
} // namespace
