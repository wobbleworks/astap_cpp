///----------------------------------------
///      @file avi.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the uncompressed RGB24 AVI streaming writer.
///   @details Byte layout notes (all little-endian):
///            - Header (hdrl LIST): sizeof = 0xC8 (200 bytes).
///            - Stream header (strl LIST + strh + strf): sizeof = 0x84 (132 bytes).
///            - movi LIST header: sizeof = 0x0C (12 bytes).
///            - Per-frame: frame_start (8) + sizeimage bytes of BGR24 + padding.
///            - idx1 header (8) + nrframes * 16-byte index entries.
///    @author Ported from Han Kleijn's ASTAP (MPL-2.0) by John Stephen.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "avi.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>

///----------------------------------------
namespace astap::image {
///----------------------------------------

namespace {

///----------------------------------------
/// MARK: FOURCC helpers
///----------------------------------------

/// @brief Build a little-endian FOURCC from 4 chars.
[[nodiscard]] constexpr std::uint32_t fourcc(char a, char b, char c, char d) noexcept {
    return  static_cast<std::uint32_t>(static_cast<std::uint8_t>(a))
         | (static_cast<std::uint32_t>(static_cast<std::uint8_t>(b)) << 8)
         | (static_cast<std::uint32_t>(static_cast<std::uint8_t>(c)) << 16)
         | (static_cast<std::uint32_t>(static_cast<std::uint8_t>(d)) << 24);
}

constexpr auto FCC_RIFF = fourcc('R','I','F','F');
constexpr auto FCC_AVI  = fourcc('A','V','I',' ');
constexpr auto FCC_LIST = fourcc('L','I','S','T');
constexpr auto FCC_HDRL = fourcc('h','d','r','l');
constexpr auto FCC_AVIH = fourcc('a','v','i','h');
constexpr auto FCC_STRL = fourcc('s','t','r','l');
constexpr auto FCC_STRH = fourcc('s','t','r','h');
constexpr auto FCC_VIDS = fourcc('v','i','d','s');
constexpr auto FCC_STRF = fourcc('s','t','r','f');
constexpr auto FCC_MOVI = fourcc('m','o','v','i');
constexpr auto FCC_00DB = fourcc('0','0','d','b');
constexpr auto FCC_IDX1 = fourcc('i','d','x','1');

constexpr auto NRCOLORS = 3; // RGB24 — mono is not standardised in AVI.

///----------------------------------------
/// MARK: Little-endian emit helpers
///----------------------------------------

/// @brief Write a 16-bit unsigned integer in little-endian order.
void put_u16(std::ofstream& f, std::uint16_t v) {
    auto b = std::array<std::uint8_t, 2>{
        static_cast<std::uint8_t>( v       & 0xFF),
        static_cast<std::uint8_t>((v >> 8) & 0xFF),
    };
    f.write(reinterpret_cast<const char*>(b.data()), 2);
}

/// @brief Write a 32-bit unsigned integer in little-endian order.
void put_u32(std::ofstream& f, std::uint32_t v) {
    auto b = std::array<std::uint8_t, 4>{
        static_cast<std::uint8_t>( v        & 0xFF),
        static_cast<std::uint8_t>((v >> 8)  & 0xFF),
        static_cast<std::uint8_t>((v >> 16) & 0xFF),
        static_cast<std::uint8_t>((v >> 24) & 0xFF),
    };
    f.write(reinterpret_cast<const char*>(b.data()), 4);
}

///----------------------------------------
/// MARK: Main AVI header (hdrl LIST + avih)
///----------------------------------------

constexpr auto HEADER_SIZE = std::uint32_t{0xC8};

struct AviHeader {
    std::uint32_t riff      = FCC_RIFF;
    std::uint32_t riffsize  = 0;
    std::uint32_t avi       = FCC_AVI;
    std::uint32_t list      = FCC_LIST;
    std::uint32_t lsize     = 0x000000C0;
    std::uint32_t hdrL      = FCC_HDRL;
    std::uint32_t fcc       = FCC_AVIH;
    std::uint32_t cb        = 0x00000038;   // 14 * 4 = 56
    std::uint32_t dwMicroSecPerFrame    = 0x000F4240; // 1 sec default
    std::uint32_t dwMaxBytesPerSec      = 0;
    std::uint32_t dwPaddingGranularity  = 0;
    std::uint32_t dwFlags               = 0x00000010;
    std::uint32_t dwTotalFrames         = 0;
    std::uint32_t dwInitialFrames       = 0;
    std::uint32_t dwStreams             = 1;
    std::uint32_t dwSuggestedBufferSize = 0;
    std::uint32_t dwWidth               = 16;
    std::uint32_t dwHeight              = 8;
    std::uint32_t dwReserved1 = 0, dwReserved2 = 0, dwReserved3 = 0, dwReserved4 = 0;
    
    void write(std::ofstream& f) const {
        put_u32(f, riff);      put_u32(f, riffsize);  put_u32(f, avi);
        put_u32(f, list);      put_u32(f, lsize);     put_u32(f, hdrL);
        put_u32(f, fcc);       put_u32(f, cb);
        put_u32(f, dwMicroSecPerFrame);    put_u32(f, dwMaxBytesPerSec);
        put_u32(f, dwPaddingGranularity);  put_u32(f, dwFlags);
        put_u32(f, dwTotalFrames);         put_u32(f, dwInitialFrames);
        put_u32(f, dwStreams);             put_u32(f, dwSuggestedBufferSize);
        put_u32(f, dwWidth);               put_u32(f, dwHeight);
        put_u32(f, dwReserved1);           put_u32(f, dwReserved2);
        put_u32(f, dwReserved3);           put_u32(f, dwReserved4);
    }
};

///----------------------------------------
/// MARK: Stream header (strl LIST + strh + strf)
///----------------------------------------

/// On-disk size of the stream header block (field-validated against the
/// original source which sets the LIST size to $74 = 116, i.e. total - 12).
constexpr auto STREAMHEADER_ON_DISK_BYTES = std::uint32_t{
      12   // LIST + size + strl
    + 8    // strh + hsize
    + 52   // fccType..dwSampleSize (13 dwords)
    + 4    // wPriority + wLanguage
    + 8    // rcframe*
    + 12   // strf + Ssize + fsize
    + 8    // width + height
    + 4    // planes + bitcount
    + 24   // compression..nr_important_colours
};

static_assert(STREAMHEADER_ON_DISK_BYTES == 132,
              "stream header size derivation");
              
struct StreamHeader {
    std::uint32_t list        = FCC_LIST;
    std::uint32_t size        = 0x74;
    std::uint32_t strl        = FCC_STRL;
    
    std::uint32_t strh        = FCC_STRH;
    std::uint32_t hsize       = 0x00000038; // 56
    
    std::uint32_t fccType     = FCC_VIDS;
    std::uint32_t fccHandler  = 0;
    std::uint32_t dwFlags     = 0;
    std::uint16_t wPriority   = 0;
    std::uint16_t wLanguage   = 0;
    std::uint32_t dwInitialFrames = 0;
    std::uint32_t dwScale     = 1;
    std::uint32_t dwRate      = 1;
    std::uint32_t dwStart     = 0;
    std::uint32_t dwLength    = 0;
    std::uint32_t dwSuggestedBufferSize = 0;
    std::uint32_t dwQuality   = 0;
    std::uint32_t dwSampleSize = 0;
    std::uint16_t rcframew1   = 0;
    std::uint16_t rcframeh1   = 0;
    std::uint16_t rcframew2   = 200;
    std::uint16_t rcframeh2   = 100;
    
    std::uint32_t strf        = FCC_STRF;
    std::uint32_t Ssize       = 40;
    std::uint32_t fsize       = 40;
    std::uint32_t width       = 200;
    std::uint32_t height      = 100;
    std::uint16_t planes      = 1;
    std::uint16_t bitcount    = 24;
    std::uint32_t compression = 0;
    std::uint32_t sizeimage   = 200 * 100 + 200;
    std::uint32_t pixels_per_meterH = 0x0EC4;
    std::uint32_t pixels_per_meterV = 0x0EC4;
    std::uint32_t nr_colours_used        = 0;
    std::uint32_t nr_important_colours   = 0;
    
    void write(std::ofstream& f) const {
        put_u32(f, list);  put_u32(f, size); put_u32(f, strl);
        put_u32(f, strh);  put_u32(f, hsize);
        put_u32(f, fccType);      put_u32(f, fccHandler);   put_u32(f, dwFlags);
        put_u16(f, wPriority);    put_u16(f, wLanguage);
        put_u32(f, dwInitialFrames);
        put_u32(f, dwScale);      put_u32(f, dwRate);       put_u32(f, dwStart);
        put_u32(f, dwLength);     put_u32(f, dwSuggestedBufferSize);
        put_u32(f, dwQuality);    put_u32(f, dwSampleSize);
        put_u16(f, rcframew1);    put_u16(f, rcframeh1);
        put_u16(f, rcframew2);    put_u16(f, rcframeh2);
        put_u32(f, strf);         put_u32(f, Ssize);        put_u32(f, fsize);
        put_u32(f, width);        put_u32(f, height);
        put_u16(f, planes);       put_u16(f, bitcount);
        put_u32(f, compression);  put_u32(f, sizeimage);
        put_u32(f, pixels_per_meterH); put_u32(f, pixels_per_meterV);
        put_u32(f, nr_colours_used);   put_u32(f, nr_important_colours);
    }
};

///----------------------------------------
/// MARK: movi LIST header
///----------------------------------------

constexpr auto MOVIHEADER_SIZE = std::uint32_t{12};

struct MoviHeader {
    std::uint32_t list = FCC_LIST;
    std::uint32_t size = 0;
    std::uint32_t movi = FCC_MOVI;
    
    void write(std::ofstream& f) const {
        put_u32(f, list); put_u32(f, size); put_u32(f, movi);
    }
};

///----------------------------------------
/// MARK: Per-frame chunk header
///----------------------------------------

constexpr auto FRAMESTART_SIZE = std::uint32_t{8};

struct FrameStart {
    std::uint32_t db = FCC_00DB;
    std::uint32_t x  = 0;
    
    void write(std::ofstream& f) const { put_u32(f, db); put_u32(f, x); }
};

///----------------------------------------
/// MARK: idx1 index records
///----------------------------------------

struct IndexStart {
    std::uint32_t idx1 = FCC_IDX1;
    std::uint32_t size = 0;
    
    void write(std::ofstream& f) const { put_u32(f, idx1); put_u32(f, size); }
};

constexpr auto INDEX_ENTRY_SIZE = std::uint32_t{16};

struct IndexEntry {
    std::uint32_t db       = FCC_00DB;
    std::uint32_t x        = 0x10;
    std::uint32_t position = 0;
    std::uint32_t size     = 0;
    
    void write(std::ofstream& f) const {
        put_u32(f, db); put_u32(f, x); put_u32(f, position); put_u32(f, size);
    }
};

///----------------------------------------
/// MARK: Legacy singleton
///----------------------------------------

/// @brief Returns the module-level singleton used by deprecated free functions.
AviWriter& legacy_writer() {
    static AviWriter w;
    return w;
}
    
} // namespace

///----------------------------------------
/// MARK: AviWriter
///----------------------------------------

AviWriter::~AviWriter() {
    if (_file.is_open()) {
        _file.close();
    }
}

bool AviWriter::open(const std::filesystem::path& path,
                     std::string_view             frame_rate,
                     int                          nr_frames,
                     int                          width,
                     int                          height) {
    // Parse frame rate
    auto fps = 0.0;
    try {
        fps = std::stod(std::string(frame_rate));
    } catch (const std::exception&) {
        return false;
    }
    fps = std::max(fps, 0.00001);
    
    // Each written scanline must be a multiple of 4 bytes; compute padding
    _extra = (width * NRCOLORS) % 4;
    if (_extra != 0) {
        _extra = 4 - _extra;
    }
    
    // Build headers with runtime-known values
    auto head = AviHeader{};
    auto sh   = StreamHeader{};
    auto mh   = MoviHeader{};
    
    head.dwWidth       = static_cast<std::uint32_t>(width);
    head.dwHeight      = static_cast<std::uint32_t>(height);
    head.dwTotalFrames = static_cast<std::uint32_t>(nr_frames);
    head.dwMicroSecPerFrame =
        static_cast<std::uint32_t>(std::llround(1'000'000.0 / fps));
        
    sh.bitcount  = static_cast<std::uint16_t>(8 * NRCOLORS);
    sh.width     = static_cast<std::uint32_t>(width);
    sh.height    = static_cast<std::uint32_t>(height);
    sh.rcframew2 = static_cast<std::uint16_t>(width);
    sh.rcframeh2 = static_cast<std::uint16_t>(height);
    
    sh.sizeimage = static_cast<std::uint32_t>(
        width * height * (sh.bitcount / 8) + height * _extra);
    sh.dwLength  = head.dwTotalFrames;
    sh.dwSuggestedBufferSize   = sh.sizeimage;
    head.dwSuggestedBufferSize = sh.sizeimage;
    
    mh.size = 4 /* length dword 'movi' */ +
              (sh.sizeimage + FRAMESTART_SIZE) * head.dwTotalFrames;
              
    // riffsize = (sizeof(head)-8) + sizeof(streamhead) + sizeof(movihead)
    //          + (frame_start + sizeimage + index_entry) * nrframes + idx1 header
    head.riffsize =
        (HEADER_SIZE - 8)
        + STREAMHEADER_ON_DISK_BYTES
        + MOVIHEADER_SIZE
        + (FRAMESTART_SIZE + sh.sizeimage + INDEX_ENTRY_SIZE) * head.dwTotalFrames
        + 8 /* idx1 header = fourcc + size dword */;
        
    // Prototype frame start (written here only to mirror original layout)
    auto frame_start_proto = FrameStart{};
    frame_start_proto.x = sh.sizeimage;
    
    // Open the file
    _file.open(path, std::ios::binary | std::ios::out | std::ios::trunc);
    if (!_file.is_open()) {
        return false;
    }
    
    // Write headers
    head.write(_file);
    sh.write(_file);
    mh.write(_file);
    
    _sizeimage = sh.sizeimage;
    return static_cast<bool>(_file);
}

bool AviWriter::write_frame(const PixelSource& src,
                            int x, int y, int w, int h) {
    if (!_file.is_open()) {
        return false;
    }
    
    // Write frame chunk header
    auto frame_start = FrameStart{};
    frame_start.x = _sizeimage;
    frame_start.write(_file);
    
    auto row = std::vector<std::uint8_t>(static_cast<std::size_t>(NRCOLORS * w));
    constexpr auto zero = std::uint32_t{0};
    
    // Walk bottom-up as AVI DIBs expect
    for (auto yy = y + h - 1; yy >= y; --yy) {
        for (auto xx = x; xx < x + w; ++xx) {
            const auto rgb = src.pixel(xx, yy);
            const auto R = static_cast<std::uint8_t>((rgb >> 16) & 0xFF);
            const auto G = static_cast<std::uint8_t>((rgb >>  8) & 0xFF);
            const auto B = static_cast<std::uint8_t>( rgb        & 0xFF);
            const auto i = NRCOLORS * (xx - x);
            row[i    ] = B;
            row[i + 1] = G;
            row[i + 2] = R;
        }
        
        // Write scanline
        _file.write(reinterpret_cast<const char*>(row.data()),
                    static_cast<std::streamsize>(row.size()));
                    
        // Write padding bytes
        if (_extra > 0) {
            _file.write(reinterpret_cast<const char*>(&zero),
                        static_cast<std::streamsize>(_extra));
        }
        
        if (!_file) {
            return false;
        }
    }
    
    return static_cast<bool>(_file);
}

void AviWriter::close(int nr_frames) {
    if (!_file.is_open()) {
        return;
    }
    
    // Write idx1 header
    auto index_start = IndexStart{};
    index_start.size = static_cast<std::uint32_t>(nr_frames) * 0x10u;
    index_start.write(_file);
    
    // Write index entries
    auto indx = IndexEntry{};
    indx.position = 0x4;
    indx.size     = _sizeimage;
    for (auto i = 1; i <= nr_frames; ++i) {
        indx.write(_file);
        indx.position += FRAMESTART_SIZE + _sizeimage;
    }
    
    _file.close();
    _sizeimage = 0;
    _extra     = 0;
}

///----------------------------------------
/// MARK: Legacy free-function shims
///----------------------------------------

bool write_avi_head(const std::filesystem::path& filen,
                    std::string_view             frame_rate,
                    int                          nrframes,
                    int                          w,
                    int                          h) {
    return legacy_writer().open(filen, frame_rate, nrframes, w, h);
}

bool write_avi_frame(const PixelSource& src, int x, int y, int w, int h) {
    return legacy_writer().write_frame(src, x, y, w, h);
}

void close_the_avi(int nrframes) {
    legacy_writer().close(nrframes);
}
    
} // namespace
