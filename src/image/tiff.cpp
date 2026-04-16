///----------------------------------------
///      @file tiff.cpp
///   @ingroup ASTAP++
///     @brief Uncompressed TIFF writer for 16/32-bit grayscale and 48/96-bit RGB.
///   @details Byte-accurate translation: the IFD layout, tag ordering, and file
///            layout are preserved exactly as in the original. All multi-byte
///            fields are written little-endian (TIFF classic "Intel byte order").
///            This matches the host byte order on all target platforms, but is
///            done explicitly for portability.
///    @author Ported from Han Kleijn's ASTAP. MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "tiff.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

///----------------------------------------
namespace astap::image {
///----------------------------------------

namespace {

/// @brief Shared buffer width in bytes.
/// @details Large enough for a single scanline at any reasonable width up to the
///          limits of this writer (120000 bytes = 10000 pixels of 3x32-bit RGB).
constexpr auto kBufWide = std::size_t{1024 * 120};

/// @brief IFD directory entry, packed to 12 bytes as required by the TIFF spec.
#pragma pack(push, 1)
struct TDirEntry {
    std::uint16_t tag;
    std::uint16_t type;
    std::uint32_t count;
    std::uint32_t value;
};
#pragma pack(pop)
static_assert(sizeof(TDirEntry) == 12, "TIFF IFD entry must be 12 bytes");

/// @brief Software identifier written into tag 0x0131.
/// @details The trailing NUL is preserved (some readers require it).
constexpr char kSoftwareName[] = "ASTAP";
constexpr auto kSoftwareNameLen = std::size_t{sizeof(kSoftwareName)};  // includes NUL

/// @brief TIFF file header: little-endian, version 42, placeholder IFD offset.
/// @details The IFD offset at bytes [4..7] is patched per-file.
constexpr std::array<std::uint8_t, 8> kTifHeaderTemplate = {
    0x49, 0x49,              // Intel byte order ("II")
    0x2a, 0x00,              // TIFF version (42)
    0x08, 0x00, 0x00, 0x00,  // Pointer to the first directory (patched later)
};

/// @brief Build number-of-entries header preceding each IFD (little-endian uint16).
[[nodiscard]] constexpr std::array<std::uint8_t, 2> make_nodirs(std::uint8_t n) noexcept {
    return {n, 0x00};
}

// IFD templates. Values marked "set later" are patched before writing.
// Indices in comments match the original source for easy cross-reference.

constexpr auto kSize16 = std::size_t{16};
using Dir16 = std::array<TDirEntry, kSize16>;
constexpr Dir16 kDirectoryBW16Template = {{
    {0x00FE, 0x0004, 0x00000001, 0x00000000},  // 0  NewSubFile
    {0x0100, 0x0003, 0x00000001, 0x00000000},  // 1  ImageWidth      (later)
    {0x0101, 0x0003, 0x00000001, 0x00000000},  // 2  ImageLength     (later)
    {0x0102, 0x0003, 0x00000001, 0x00000010},  // 3  BitsPerSample = 16
    {0x0103, 0x0003, 0x00000001, 0x00000001},  // 4  Compression = none
    {0x0106, 0x0003, 0x00000001, 0x00000001},  // 5  Photometric = BlackIsZero
    {0x010E, 0x0002, 0x0000000A, 0x00000000},  // 6  ImageDescription (later)
    {0x0111, 0x0004, 0x00000001, 0x00000000},  // 7  StripOffsets     (later)
    {0x0115, 0x0003, 0x00000001, 0x00000001},  // 8  SamplesPerPixel = 1
    {0x0116, 0x0004, 0x00000001, 0x00000000},  // 9  RowsPerStrip     (later)
    {0x0117, 0x0004, 0x00000001, 0x00000000},  // 10 StripByteCounts  (later)
    {0x011A, 0x0005, 0x00000001, 0x00000000},  // 11 X-Resolution     (later)
    {0x011B, 0x0005, 0x00000001, 0x00000000},  // 12 Y-Resolution     (later)
    {0x0128, 0x0003, 0x00000001, 0x00000002},  // 13 ResolutionUnit = inch
    {0x0131, 0x0002, 0x0000000A, 0x00000000},  // 14 Software         (later)
    {0x0153, 0x0003, 0x00000001, 0x00000001},  // 15 SampleFormat = uint
}};

constexpr auto kSize32 = std::size_t{16};
using Dir32 = std::array<TDirEntry, kSize32>;
constexpr Dir32 kDirectoryBW32Template = {{
    {0x00FE, 0x0004, 0x00000001, 0x00000000},  // 0  NewSubFile
    {0x0100, 0x0003, 0x00000001, 0x00000000},  // 1  ImageWidth      (later)
    {0x0101, 0x0003, 0x00000001, 0x00000000},  // 2  ImageLength     (later)
    {0x0102, 0x0003, 0x00000001, 0x00000020},  // 3  BitsPerSample = 32
    {0x0103, 0x0003, 0x00000001, 0x00000001},  // 4  Compression = none
    {0x0106, 0x0003, 0x00000001, 0x00000001},  // 5  Photometric = BlackIsZero
    {0x010E, 0x0002, 0x0000000A, 0x00000000},  // 6  ImageDescription (later)
    {0x0111, 0x0004, 0x00000001, 0x00000000},  // 7  StripOffsets     (later)
    {0x0115, 0x0003, 0x00000001, 0x00000001},  // 8  SamplesPerPixel = 1
    {0x0116, 0x0004, 0x00000001, 0x00000000},  // 9  RowsPerStrip     (later)
    {0x0117, 0x0004, 0x00000001, 0x00000000},  // 10 StripByteCounts  (later)
    {0x011A, 0x0005, 0x00000001, 0x00000000},  // 11 X-Resolution     (later)
    {0x011B, 0x0005, 0x00000001, 0x00000000},  // 12 Y-Resolution     (later)
    {0x0128, 0x0003, 0x00000001, 0x00000002},  // 13 ResolutionUnit = inch
    {0x0131, 0x0002, 0x0000000A, 0x00000000},  // 14 Software         (later)
    {0x0153, 0x0003, 0x00000001, 0x00000003},  // 15 SampleFormat = float
}};

constexpr auto kSize48 = std::size_t{17};
using Dir48 = std::array<TDirEntry, kSize48>;
constexpr Dir48 kDirectoryRGB48Template = {{
    {0x00FE, 0x0004, 0x00000001, 0x00000000},  // 0  NewSubFile
    {0x0100, 0x0003, 0x00000001, 0x00000000},  // 1  ImageWidth      (later)
    {0x0101, 0x0003, 0x00000001, 0x00000000},  // 2  ImageLength     (later)
    {0x0102, 0x0003, 0x00000003, 0x00000000},  // 3  BitsPerSample  -> offset (later)
    {0x0103, 0x0003, 0x00000001, 0x00000001},  // 4  Compression = none
    {0x0106, 0x0003, 0x00000001, 0x00000002},  // 5  Photometric = RGB
    {0x010E, 0x0002, 0x0000000A, 0x00000000},  // 6  ImageDescription (later)
    {0x0111, 0x0004, 0x00000001, 0x00000000},  // 7  StripOffsets     (later)
    {0x0115, 0x0003, 0x00000001, 0x00000003},  // 8  SamplesPerPixel = 3
    {0x0116, 0x0004, 0x00000001, 0x00000000},  // 9  RowsPerStrip     (later)
    {0x0117, 0x0004, 0x00000001, 0x00000000},  // 10 StripByteCounts  (later)
    {0x011A, 0x0005, 0x00000001, 0x00000000},  // 11 X-Resolution     (later)
    {0x011B, 0x0005, 0x00000001, 0x00000000},  // 12 Y-Resolution     (later)
    {0x011C, 0x0003, 0x00000001, 0x00000001},  // 13 PlanarConfig = chunky
    {0x0128, 0x0003, 0x00000001, 0x00000002},  // 14 ResolutionUnit = inch
    {0x0131, 0x0002, 0x0000000A, 0x00000000},  // 15 Software         (later)
    {0x0153, 0x0003, 0x00000001, 0x00000001},  // 16 SampleFormat = uint
}};

constexpr auto kSize96 = std::size_t{17};
using Dir96 = std::array<TDirEntry, kSize96>;
constexpr Dir96 kDirectoryRGB96Template = {{
    {0x00FE, 0x0004, 0x00000001, 0x00000000},  // 0  NewSubFile
    {0x0100, 0x0003, 0x00000001, 0x00000000},  // 1  ImageWidth      (later)
    {0x0101, 0x0003, 0x00000001, 0x00000000},  // 2  ImageLength     (later)
    {0x0102, 0x0003, 0x00000003, 0x00000000},  // 3  BitsPerSample  -> offset (later)
    {0x0103, 0x0003, 0x00000001, 0x00000001},  // 4  Compression = none
    {0x0106, 0x0003, 0x00000001, 0x00000002},  // 5  Photometric = RGB
    {0x010E, 0x0002, 0x0000000A, 0x00000000},  // 6  ImageDescription (later)
    {0x0111, 0x0004, 0x00000001, 0x00000000},  // 7  StripOffsets     (later)
    {0x0115, 0x0003, 0x00000001, 0x00000003},  // 8  SamplesPerPixel = 3
    {0x0116, 0x0004, 0x00000001, 0x00000000},  // 9  RowsPerStrip     (later)
    {0x0117, 0x0004, 0x00000001, 0x00000000},  // 10 StripByteCounts  (later)
    {0x011A, 0x0005, 0x00000001, 0x00000000},  // 11 X-Resolution     (later)
    {0x011B, 0x0005, 0x00000001, 0x00000000},  // 12 Y-Resolution     (later)
    {0x011C, 0x0003, 0x00000001, 0x00000001},  // 13 PlanarConfig = chunky
    {0x0128, 0x0003, 0x00000001, 0x00000002},  // 14 ResolutionUnit = inch
    {0x0131, 0x0002, 0x0000000A, 0x00000000},  // 15 Software         (later)
    {0x0153, 0x0003, 0x00000001, 0x00000003},  // 16 SampleFormat = float
}};

constexpr std::array<std::uint8_t, 4> kNullString = {0x00, 0x00, 0x00, 0x00};

/// @brief X/Y-Resolution value: rational 877/10 = 87.7 dpi.
constexpr std::array<std::uint8_t, 8> kXResValue = {
    0x6D, 0x03, 0x00, 0x00, 0x0A, 0x00, 0x00, 0x00};
constexpr std::array<std::uint8_t, 8> kYResValue = {
    0x6D, 0x03, 0x00, 0x00, 0x0A, 0x00, 0x00, 0x00};
    
/// @brief BitsPerSample arrays for RGB variants.
constexpr std::array<std::uint16_t, 3> kBitsPerSample48 = {0x0010, 0x0010, 0x0010};
constexpr std::array<std::uint16_t, 3> kBitsPerSample96 = {0x0020, 0x0020, 0x0020};

/// @brief Write a contiguous byte sequence.
/// @return False on stream error.
template <class T, std::size_t N>
[[nodiscard]] bool write_bytes(std::ofstream& f, const std::array<T, N>& arr) {
    f.write(std::bit_cast<const char*>(arr.data()),
            static_cast<std::streamsize>(sizeof(T) * N));
    return static_cast<bool>(f);
}

/// @brief Write raw bytes from an arbitrary pointer.
/// @return False on stream error.
[[nodiscard]] bool write_raw(std::ofstream& f, const void* data, std::size_t size) {
    f.write(static_cast<const char*>(data), static_cast<std::streamsize>(size));
    return static_cast<bool>(f);
}

/// @brief Write an IFD directory array.
/// @return False on stream error.
template <class Dir>
[[nodiscard]] bool write_dir(std::ofstream& f, const Dir& d) {
    return write_raw(f, d.data(), sizeof(TDirEntry) * d.size());
}

/// @brief Patch a little-endian uint32 into a 4-byte slot of the header.
void patch_le32(std::array<std::uint8_t, 8>& header, std::size_t offset,
                std::uint32_t value) noexcept {
    header[offset + 0] = static_cast<std::uint8_t>(value & 0xFF);
    header[offset + 1] = static_cast<std::uint8_t>((value >> 8) & 0xFF);
    header[offset + 2] = static_cast<std::uint8_t>((value >> 16) & 0xFF);
    header[offset + 3] = static_cast<std::uint8_t>((value >> 24) & 0xFF);
}

/// @brief Compose description bytes as "<description>\0" -- some readers require the trailing NUL.
[[nodiscard]] std::string with_nul(std::string_view description) {
    auto s = std::string{};
    s.reserve(description.size() + 1);
    s.append(description);
    s.push_back('\0');
    return s;
}

/// @brief Adjust extension to ".tif".
[[nodiscard]] std::filesystem::path with_tif_ext(const std::filesystem::path& p) {
    auto out = p;
    out.replace_extension(".tif");
    return out;
}

/// @brief Open for binary write (overwrites).
/// @return An open stream, or a closed stream on failure.
[[nodiscard]] std::ofstream open_create(const std::filesystem::path& p) {
    auto f = std::ofstream{p, std::ios::binary | std::ios::trunc};
    return f;
}

/// @brief Quantize a single-precision sample to uint16 with saturate + round.
[[nodiscard]] std::uint16_t to_u16(float sample) noexcept {
    auto dum = static_cast<double>(sample);
    if (dum > 65535.0) {
        dum = 65535.0;
    }
    if (dum < 0.0) {
        dum = 0.0;
    }
    return static_cast<std::uint16_t>(std::lround(dum));
}
    
}  // namespace

/// MARK: - save_tiff_16

bool save_tiff_16(const ImageArray& img,
                  const std::filesystem::path& filename,
                  std::string_view description,
                  bool flip_h,
                  bool flip_v) {
    if (img.empty() || img[0].empty() || img[0][0].empty()) {
        return false;
    }
    
    const auto path = with_tif_ext(filename);
    
    // Callers must handle the "overwrite?" decision before invoking.
    const auto width2 = static_cast<int>(img[0][0].size());
    const auto height2 = static_cast<int>(img[0].size());
    const auto desc = with_nul(description);
    
    auto f = open_create(path);
    if (!f) {
        return false;
    }
    
    // Prepare per-file header and IFD copy.
    auto header = kTifHeaderTemplate;
    auto dir = kDirectoryBW16Template;
    
    dir[1].value = static_cast<std::uint32_t>(width2);
    dir[2].value = static_cast<std::uint32_t>(height2);
    dir[9].value = static_cast<std::uint32_t>(height2);
    dir[10].value = static_cast<std::uint32_t>(2 * width2 * height2);
    dir[6].count = static_cast<std::uint32_t>(desc.size());
    dir[14].count = static_cast<std::uint32_t>(kSoftwareNameLen);
    
    // Layout: [Header][XRes][YRes][Descr][Software][IFD]...
    const auto offset_dir =
        static_cast<std::uint32_t>(kTifHeaderTemplate.size() +
                                   kXResValue.size() + kYResValue.size() +
                                   desc.size() + kSoftwareNameLen);
    patch_le32(header, 4, offset_dir);
    
    if (!write_bytes(f, header)) {
        return false;
    }
    
    const auto offset_x_res = static_cast<std::uint32_t>(f.tellp());
    if (!write_bytes(f, kXResValue)) {
        return false;
    }
    
    const auto offset_y_res = static_cast<std::uint32_t>(f.tellp());
    if (!write_bytes(f, kYResValue)) {
        return false;
    }
    
    const auto offset_descrip = static_cast<std::uint32_t>(f.tellp());
    if (!write_raw(f, desc.data(), desc.size())) {
        return false;
    }
    
    const auto offset_software = static_cast<std::uint32_t>(f.tellp());
    if (!write_raw(f, kSoftwareName, kSoftwareNameLen)) {
        return false;
    }
    
    const auto nodirs = make_nodirs(static_cast<std::uint8_t>(kSize16));
    const auto offset_strip =
        offset_dir + static_cast<std::uint32_t>(nodirs.size()) +
        static_cast<std::uint32_t>(sizeof(TDirEntry) * kSize16) +
        static_cast<std::uint32_t>(kNullString.size());
        
    dir[7].value = offset_strip;
    dir[11].value = offset_x_res;
    dir[12].value = offset_y_res;
    dir[14].value = offset_software;
    dir[6].value = offset_descrip;
    
    if (!write_bytes(f, nodirs)) {
        return false;
    }
    if (!write_dir(f, dir)) {
        return false;
    }
    if (!write_bytes(f, kNullString)) {
        return false;
    }
    
    // Write scanlines. kBufWide accommodates width2 * 2 bytes up to 61440 px.
    auto buf = std::array<std::uint8_t, kBufWide>{};
    for (auto i = 0; i < height2; ++i) {
        const auto k = flip_v ? i : (height2 - 1 - i);
        for (auto j = 0; j < width2; ++j) {
            const auto m = flip_h ? (width2 - 1 - j) : j;
            const auto v = to_u16(img[0][k][m]);
            buf[static_cast<std::size_t>(m) * 2 + 0] =
                static_cast<std::uint8_t>(v & 0xFF);
            buf[static_cast<std::size_t>(m) * 2 + 1] =
                static_cast<std::uint8_t>((v >> 8) & 0xFF);
        }
        if (!write_raw(f, buf.data(), static_cast<std::size_t>(width2) * 2)) {
            return false;
        }
    }
    
    f.flush();
    return static_cast<bool>(f);
}

/// MARK: - save_tiff_32

bool save_tiff_32(const ImageArray& img,
                  const std::filesystem::path& filename,
                  std::string_view description,
                  bool flip_h,
                  bool flip_v) {
    if (img.empty() || img[0].empty() || img[0][0].empty()) {
        return false;
    }
    
    const auto path = with_tif_ext(filename);
    const auto width2 = static_cast<int>(img[0][0].size());
    const auto height2 = static_cast<int>(img[0].size());
    const auto desc = with_nul(description);
    
    auto f = open_create(path);
    if (!f) {
        return false;
    }
    
    auto header = kTifHeaderTemplate;
    auto dir = kDirectoryBW32Template;
    
    dir[1].value = static_cast<std::uint32_t>(width2);
    dir[2].value = static_cast<std::uint32_t>(height2);
    dir[9].value = static_cast<std::uint32_t>(height2);
    dir[10].value = static_cast<std::uint32_t>(4 * width2 * height2);
    dir[6].count = static_cast<std::uint32_t>(desc.size());
    dir[14].count = static_cast<std::uint32_t>(kSoftwareNameLen);
    
    const auto offset_dir =
        static_cast<std::uint32_t>(kTifHeaderTemplate.size() +
                                   kXResValue.size() + kYResValue.size() +
                                   desc.size() + kSoftwareNameLen);
    patch_le32(header, 4, offset_dir);
    
    if (!write_bytes(f, header)) {
        return false;
    }
    
    const auto offset_x_res = static_cast<std::uint32_t>(f.tellp());
    if (!write_bytes(f, kXResValue)) {
        return false;
    }
    
    const auto offset_y_res = static_cast<std::uint32_t>(f.tellp());
    if (!write_bytes(f, kYResValue)) {
        return false;
    }
    
    const auto offset_descrip = static_cast<std::uint32_t>(f.tellp());
    if (!write_raw(f, desc.data(), desc.size())) {
        return false;
    }
    
    const auto offset_software = static_cast<std::uint32_t>(f.tellp());
    if (!write_raw(f, kSoftwareName, kSoftwareNameLen)) {
        return false;
    }
    
    const auto nodirs = make_nodirs(static_cast<std::uint8_t>(kSize32));
    const auto offset_strip =
        offset_dir + static_cast<std::uint32_t>(nodirs.size()) +
        static_cast<std::uint32_t>(sizeof(TDirEntry) * kSize32) +
        static_cast<std::uint32_t>(kNullString.size());
        
    dir[7].value = offset_strip;
    dir[11].value = offset_x_res;
    dir[12].value = offset_y_res;
    dir[14].value = offset_software;
    dir[6].value = offset_descrip;
    
    if (!write_bytes(f, nodirs)) {
        return false;
    }
    if (!write_dir(f, dir)) {
        return false;
    }
    if (!write_bytes(f, kNullString)) {
        return false;
    }
    
    // Write scanlines as float, normalised to [0,1].
    static_assert(sizeof(float) == 4);
    auto buf = std::array<float, kBufWide / 4>{};
    for (auto i = 0; i < height2; ++i) {
        const auto k = flip_v ? i : (height2 - 1 - i);
        for (auto j = 0; j < width2; ++j) {
            const auto m = flip_h ? (width2 - 1 - j) : j;
            buf[static_cast<std::size_t>(m)] = img[0][k][m] / 65535.0f;
        }
        if (!write_raw(f, buf.data(), static_cast<std::size_t>(width2) * 4)) {
            return false;
        }
    }
    
    f.flush();
    return static_cast<bool>(f);
}

/// MARK: - save_tiff_48

bool save_tiff_48(const ImageArray& img,
                  const std::filesystem::path& filename,
                  std::string_view description,
                  bool flip_h,
                  bool flip_v) {
    if (img.size() < 3 || img[0].empty() || img[0][0].empty()) {
        return false;
    }
    
    const auto path = with_tif_ext(filename);
    const auto width2 = static_cast<int>(img[0][0].size());
    const auto height2 = static_cast<int>(img[0].size());
    const auto desc = with_nul(description);
    
    auto f = open_create(path);
    if (!f) {
        return false;
    }
    
    auto header = kTifHeaderTemplate;
    auto dir = kDirectoryRGB48Template;
    
    dir[1].value = static_cast<std::uint32_t>(width2);
    dir[2].value = static_cast<std::uint32_t>(height2);
    dir[9].value = static_cast<std::uint32_t>(height2);
    dir[10].value = static_cast<std::uint32_t>(2 * 3 * width2 * height2);
    dir[6].count = static_cast<std::uint32_t>(desc.size());
    dir[15].count = static_cast<std::uint32_t>(kSoftwareNameLen);
    
    // Layout: [Header][XRes][YRes][BitsPerSample][Descr][Software][IFD]...
    const auto offset_dir = static_cast<std::uint32_t>(
        kTifHeaderTemplate.size() + kXResValue.size() + kYResValue.size() +
        sizeof(kBitsPerSample48) + desc.size() + kSoftwareNameLen);
    patch_le32(header, 4, offset_dir);
    
    if (!write_bytes(f, header)) {
        return false;
    }
    
    const auto offset_x_res = static_cast<std::uint32_t>(f.tellp());
    if (!write_bytes(f, kXResValue)) {
        return false;
    }
    
    const auto offset_y_res = static_cast<std::uint32_t>(f.tellp());
    if (!write_bytes(f, kYResValue)) {
        return false;
    }
    
    const auto offset_bps = static_cast<std::uint32_t>(f.tellp());
    if (!write_bytes(f, kBitsPerSample48)) {
        return false;
    }
    
    const auto offset_descrip = static_cast<std::uint32_t>(f.tellp());
    if (!write_raw(f, desc.data(), desc.size())) {
        return false;
    }
    
    const auto offset_software = static_cast<std::uint32_t>(f.tellp());
    if (!write_raw(f, kSoftwareName, kSoftwareNameLen)) {
        return false;
    }
    
    const auto nodirs = make_nodirs(static_cast<std::uint8_t>(kSize48));
    const auto offset_strip =
        offset_dir + static_cast<std::uint32_t>(nodirs.size()) +
        static_cast<std::uint32_t>(sizeof(TDirEntry) * kSize48) +
        static_cast<std::uint32_t>(kNullString.size());
        
    dir[3].value = offset_bps;
    dir[7].value = offset_strip;
    dir[11].value = offset_x_res;
    dir[12].value = offset_y_res;
    dir[15].value = offset_software;
    dir[6].value = offset_descrip;
    
    if (!write_bytes(f, nodirs)) {
        return false;
    }
    if (!write_dir(f, dir)) {
        return false;
    }
    if (!write_bytes(f, kNullString)) {
        return false;
    }
    
    // Write interleaved RGB scanlines.
    auto buf = std::array<std::uint8_t, kBufWide>{};
    for (auto i = 0; i < height2; ++i) {
        const auto k = flip_v ? i : (height2 - 1 - i);
        for (auto j = 0; j < width2; ++j) {
            const auto m = flip_h ? (width2 - 1 - j) : j;
            const auto base = static_cast<std::size_t>(m) * 6;
            for (auto c = 0; c < 3; ++c) {
                const auto v = to_u16(img[c][k][m]);
                buf[base + c * 2 + 0] = static_cast<std::uint8_t>(v & 0xFF);
                buf[base + c * 2 + 1] = static_cast<std::uint8_t>((v >> 8) & 0xFF);
            }
        }
        if (!write_raw(f, buf.data(), static_cast<std::size_t>(width2) * 6)) {
            return false;
        }
    }
    
    // Write a duplicate IFD at end-of-file (matches original for bit-for-bit compatibility).
    if (!write_bytes(f, nodirs)) {
        return false;
    }
    if (!write_dir(f, dir)) {
        return false;
    }
    if (!write_bytes(f, kNullString)) {
        return false;
    }
    
    f.flush();
    return static_cast<bool>(f);
}

/// MARK: - save_tiff_96

bool save_tiff_96(const ImageArray& img,
                  const std::filesystem::path& filename,
                  std::string_view description,
                  bool flip_h,
                  bool flip_v) {
    if (img.size() < 3 || img[0].empty() || img[0][0].empty()) {
        return false;
    }
    
    const auto path = with_tif_ext(filename);
    const auto width2 = static_cast<int>(img[0][0].size());
    const auto height2 = static_cast<int>(img[0].size());
    const auto desc = with_nul(description);
    
    auto f = open_create(path);
    if (!f) {
        return false;
    }
    
    auto header = kTifHeaderTemplate;
    auto dir = kDirectoryRGB96Template;
    
    dir[1].value = static_cast<std::uint32_t>(width2);
    dir[2].value = static_cast<std::uint32_t>(height2);
    dir[9].value = static_cast<std::uint32_t>(height2);
    dir[10].value = static_cast<std::uint32_t>(4 * 3 * width2 * height2);
    dir[6].count = static_cast<std::uint32_t>(desc.size());
    dir[15].count = static_cast<std::uint32_t>(kSoftwareNameLen);
    
    const auto offset_dir = static_cast<std::uint32_t>(
        kTifHeaderTemplate.size() + kXResValue.size() + kYResValue.size() +
        sizeof(kBitsPerSample96) + desc.size() + kSoftwareNameLen);
    patch_le32(header, 4, offset_dir);
    
    if (!write_bytes(f, header)) {
        return false;
    }
    
    const auto offset_x_res = static_cast<std::uint32_t>(f.tellp());
    if (!write_bytes(f, kXResValue)) {
        return false;
    }
    
    const auto offset_y_res = static_cast<std::uint32_t>(f.tellp());
    if (!write_bytes(f, kYResValue)) {
        return false;
    }
    
    const auto offset_bps = static_cast<std::uint32_t>(f.tellp());
    if (!write_bytes(f, kBitsPerSample96)) {
        return false;
    }
    
    const auto offset_descrip = static_cast<std::uint32_t>(f.tellp());
    if (!write_raw(f, desc.data(), desc.size())) {
        return false;
    }
    
    const auto offset_software = static_cast<std::uint32_t>(f.tellp());
    if (!write_raw(f, kSoftwareName, kSoftwareNameLen)) {
        return false;
    }
    
    const auto nodirs = make_nodirs(static_cast<std::uint8_t>(kSize96));
    const auto offset_strip =
        offset_dir + static_cast<std::uint32_t>(nodirs.size()) +
        static_cast<std::uint32_t>(sizeof(TDirEntry) * kSize96) +
        static_cast<std::uint32_t>(kNullString.size());
        
    dir[3].value = offset_bps;
    dir[7].value = offset_strip;
    dir[11].value = offset_x_res;
    dir[12].value = offset_y_res;
    dir[15].value = offset_software;
    dir[6].value = offset_descrip;
    
    if (!write_bytes(f, nodirs)) {
        return false;
    }
    if (!write_dir(f, dir)) {
        return false;
    }
    if (!write_bytes(f, kNullString)) {
        return false;
    }
    
    // Write interleaved 3x float scanlines, each sample normalised to [0,1]
    // via bit_cast (avoids UB from type-punning).
    static_assert(sizeof(float) == 4);
    auto buf = std::array<std::uint8_t, kBufWide>{};
    for (auto i = 0; i < height2; ++i) {
        const auto k = flip_v ? i : (height2 - 1 - i);
        for (auto j = 0; j < width2; ++j) {
            const auto m = flip_h ? (width2 - 1 - j) : j;
            const auto base = static_cast<std::size_t>(m) * 12;
            for (auto c = 0; c < 3; ++c) {
                const auto sample = img[c][k][m] / 65535.0f;
                const auto raw = std::bit_cast<std::array<std::uint8_t, 4>>(sample);
                buf[base + c * 4 + 0] = raw[0];
                buf[base + c * 4 + 1] = raw[1];
                buf[base + c * 4 + 2] = raw[2];
                buf[base + c * 4 + 3] = raw[3];
            }
        }
        if (!write_raw(f, buf.data(), static_cast<std::size_t>(width2) * 12)) {
            return false;
        }
    }
    
    f.flush();
    return static_cast<bool>(f);
}
    
} // namespace
