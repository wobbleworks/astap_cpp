///----------------------------------------
///      @file image_io.cpp
///   @ingroup ASTAP++
///     @brief Non-FITS image I/O -- implementation.
///   @details See image_io.h for the public-API doc-comments.
///    @author Ported from Han Kleijn's ASTAP; MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "image_io.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include "globals.h"
#include "../image/tiff.h"
#include "../image/xisf.h"

#include <memory>

///----------------------------------------
namespace astap::core {
///----------------------------------------

using astap::nrbits;
using astap::extend_type;
using astap::instrum;
using astap::bck;
using astap::cwhite;
using astap::head1;
using astap::unsaved_import;
using astap::application_path;
using astap::raw_conversion_program_index;
using astap::esc_pressed;
using astap::commandline_execution;
using astap::head;
using astap::img_loaded;
using astap::recent_files;

// ---------------------------------------------------------------------------
// Null decoder: the factory default. Refuses all requests with a helpful
// error so failures surface at decode time (not at link time).
// Hosts override this by calling set_image_decoder() at startup.
// ---------------------------------------------------------------------------
namespace {

struct NullDecoder : IImageDecoder {
    bool decode_raster(const std::filesystem::path& /*path*/,
                       std::string_view             /*ext_upper*/,
                       DecodedImage&                /*out*/,
                       std::string&                 error_out) override {
        error_out = "no image decoder installed; call astap::core::set_image_decoder()";
        return false;
    }
    bool decode_raw(const std::filesystem::path& /*path*/,
                    bool                         /*save_intermediate*/,
                    DecodedImage&                /*out*/,
                    std::string&                 error_out) override {
        error_out = "no RAW decoder installed; call astap::core::set_image_decoder()";
        return false;
    }
};

// Storage for the currently-installed decoder. Initialised lazily so a null
// implementation always exists, even during static init (call order is
// otherwise unspecified across translation units).
std::shared_ptr<IImageDecoder>& decoder_slot() {
    static std::shared_ptr<IImageDecoder> s_decoder = std::make_shared<NullDecoder>();
    return s_decoder;
}
    
}  // anonymous namespace

IImageDecoder& image_decoder() noexcept {
    return *decoder_slot();
}

void set_image_decoder(std::shared_ptr<IImageDecoder> impl) {
    decoder_slot() = impl ? std::move(impl) : std::make_shared<NullDecoder>();
}

// ---------------------------------------------------------------------------
// Forward declarations of helpers that live in sister translation units.
// These are populated as the surrounding modules are ported; until then the
// stubs in their respective .cpp files keep the link viable.
// ---------------------------------------------------------------------------

// FITS I/O (fits.h, being ported in parallel).
bool load_fits(const std::filesystem::path& filename,
               bool                         light,
               bool                         load_data,
               bool                         update_memo,
               int                          hdu_index,
               std::vector<std::string>&    memo,
               Header&                      head,
               ImageArray&                  img);
               
bool save_fits(const ImageArray&               img,
               const std::vector<std::string>& memo,
               const std::filesystem::path&    filename,
               int                             nrbits,
               bool                            overwrite);
               
// .fz CFITSIO unpack — wraps `funpack` shell-out. Lives in fits.cpp.
bool unpack_cfitsio(std::string& filename);

// Extension-classification helpers (util.h).
bool check_raw_file_extension(std::string_view ext);

// Header-text helpers (header_utils.h, also being ported).
void reset_fits_global_variables(bool light, Header& head);
void update_integer(std::vector<std::string>& memo,
                    std::string_view          key,
                    std::string_view          comment,
                    long long                 value);
void update_float(std::vector<std::string>& memo,
                  std::string_view          key,
                  std::string_view          comment,
                  bool                      exponent,
                  double                    value);
void update_text(std::vector<std::string>& memo,
                 std::string_view          key,
                 std::string_view          value);
void add_text(std::vector<std::string>& memo,
              std::string_view          key,
              std::string_view          value);
void add_long_comment(std::vector<std::string>& memo, std::string_view text);
void read_keys_memo(bool light, Header& head, std::vector<std::string>& memo);

// Filename heuristics — used to back-fill metadata when the file lacks it.
double      extract_exposure_from_filename(std::string_view filename);
int         extract_temperature_from_filename(std::string_view filename);
std::string extract_objectname_from_filename(std::string_view filename);

// Calendar conversion. JD -> "YYYY-MM-DDThh:mm:ss".
std::string jd_to_date(double jd);
double      file_age_jd(const std::filesystem::path& path);

// Platform (platform.cpp). Show-console is the "show window" flag.
bool execute_and_wait(const std::string& cmd, bool show_console);

// Recent-file menu refresh (GUI layer). No-op in headless builds.
void update_recent_file_menu();

// GUI / status notifications. In headless builds these are no-ops; the
// concrete implementations live alongside the respective UIs.
void show_error_label(std::string_view text);
void memo2_message(std::string_view text);

// Shared globals (nrbits, extend_type, instrum, bck, cwhite, head1,
// unsaved_import, application_path, raw_conversion_program_index,
// esc_pressed, commandline_execution) live in ../core/globals.h which is
// included via this module's header.
//
// `memox` is a scratch memo pointer used by convert_to_fits to inject lines
// into the caller's memo without touching mainwindow.Memo1; it stays local.
extern std::vector<std::string>* memox;

// ---------------------------------------------------------------------------
// Internal helpers.
// ---------------------------------------------------------------------------
namespace {

constexpr std::size_t kReaderBufferSize = 0x60000;  // 393216 bytes.

inline std::string to_upper(std::string_view s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    }
    return out;
}

inline std::string ext_upper(const std::filesystem::path& p) {
    return to_upper(p.extension().string());
}

inline std::string change_ext(const std::filesystem::path& p, std::string_view new_ext) {
    auto q = p;
    q.replace_extension(new_ext);
    return q.string();
}

// Tolerant integer parser — returns 0 on failure and writes err non-zero.
inline long val_int(std::string_view s, int& err) {
    err = 0;
    if (s.empty()) { err = 1; return 0; }
    char* end = nullptr;
    auto str = std::string(s);
    long v = std::strtol(str.c_str(), &end, 10);
    if (end == str.c_str() || *end != '\0') err = 1;
    return v;
}

inline double val_double(std::string_view s, int& err) {
    err = 0;
    if (s.empty()) { err = 1; return 0.0; }
    char* end = nullptr;
    auto str = std::string(s);
    double v = std::strtod(str.c_str(), &end);
    if (end == str.c_str() || *end != '\0') err = 1;
    return v;
}

// Swap the two bytes of a 16-bit value.
inline std::uint16_t swap16(std::uint16_t v) {
    return static_cast<std::uint16_t>((v << 8) | (v >> 8));
}

// Tolerant float parser that accepts either '.' or ',' as decimal separator.
inline double strtofloat2(std::string_view s) {
    std::string copy(s);
    for (char& c : copy) if (c == ',') c = '.';
    int err = 0;
    return val_double(copy, err);
}

// Convert a unix timestamp (seconds since 1970-01-01) to JD.
inline double unix_to_jd(long long secs) {
    return 2440587.5 + static_cast<double>(secs) / 86400.0;
}

// Simple stream reader wrapping ifstream with position tracking.
struct ReaderStream {
    std::ifstream f;
    int           position{0};
    
    explicit ReaderStream(const std::filesystem::path& path)
        : f(path, std::ios::binary) {
        f.rdbuf()->pubsetbuf(nullptr, 0);  // unbuffered; we manage chunking.
    }
    
    bool ok() const { return f.good() || f.eof(); }
    bool is_open() const { return f.is_open(); }
    
    bool read_byte(std::uint8_t& out) {
        char c{};
        if (!f.read(&c, 1)) return false;
        out = static_cast<std::uint8_t>(c);
        ++position;
        return true;
    }
    
    std::size_t read(void* buf, std::size_t n) {
        f.read(static_cast<char*>(buf), static_cast<std::streamsize>(n));
        auto got = static_cast<std::size_t>(f.gcount());
        position += static_cast<int>(got);
        return got;
    }
};
    
}  // namespace

// ===========================================================================
// load_PPM_PGM_PFM
// ===========================================================================
bool load_ppm_pgm_pfm(const std::filesystem::path& filen,
                      Header&                      head,
                      ImageArray&                  img_loaded2,
                      std::vector<std::string>&    memo) {
    head.naxis = 0;
    
    ReaderStream reader(filen);
    if (!reader.is_open()) {
        show_error_label("Error, accessing the file!");
        return false;
    }
    
    memo.clear();
    
    reset_fits_global_variables(true, head);
    
    // ---- Magic number: "P5\n", "P6\n", "PF\n", "Pf\n" -----------------------
    std::array<std::uint8_t, 3> magic{};
    for (int i = 0; i < 3; ++i) {
        if (!reader.read_byte(magic[i])) {
            show_error_label("Error loading PGM/PPM/PFM file!! Keyword P5, P6, PF. Pf not found.");
            return false;
        }
    }
    
    bool color7 = false;
    bool pfm    = false;
    if      (magic[0] == 'P' && magic[1] == '5' && magic[2] == 0x0A) { color7 = false; }
    else if (magic[0] == 'P' && magic[1] == '6' && magic[2] == 0x0A) { color7 = true;  }
    else if (magic[0] == 'P' && magic[1] == 'F' && magic[2] == 0x0A) { color7 = true;  pfm = true; }
    else if (magic[0] == 'P' && magic[1] == 'f' && magic[2] == 0x0A) { color7 = false; pfm = true; }
    else {
        show_error_label("Error loading PGM/PPM/PFM file!! Keyword P5, P6, PF. Pf not found.");
        return false;
    }
    
    // ---- ASCII header lines: width / height / bits, with #comment lines ----
    std::string w1, h1, bits;
    int header_field = 0;  // 1=width, 2=height, 3=bits
    
    auto try_strtoint = [](std::string_view s) -> long {
        int err = 0;
        return val_int(s, err);
    };
    
    bool expdet      = false;
    bool timedet     = false;
    bool isodet      = false;
    bool instdet     = false;
    bool ccdtempdet  = false;
    
    while (header_field < 3 && reader.position <= 200) {
        bool        comment = false;
        std::string aline;
        std::string comm;
        
        for (;;) {
            std::uint8_t ch{};
            if (!reader.read_byte(ch)) break;
            
            if (ch == '#') comment = true;
            
            if (comment) {
                bool is_delim = (ch == ';' || ch == '#' || ch == ' ' || ch == 0x0A);
                if (!is_delim) {
                    comm.push_back(static_cast<char>(ch));
                } else {
                    if (expdet)     { head.exposure        = strtofloat2(comm); expdet = false; }
                    if (isodet)     { head.gain            = comm;              isodet = false; }
                    if (instdet)    { instrum              = comm;              instdet = false; }
                    if (ccdtempdet) { head.set_temperature = static_cast<int>(std::lround(strtofloat2(comm))); ccdtempdet = false; }
                    if (timedet) {
                        const double jd2 = unix_to_jd(try_strtoint(comm));
                        head.date_obs = jd_to_date(jd2);
                        timedet = false;
                    }
                    comm.clear();
                }
                if      (comm == "EXPTIME=")   { expdet     = true; comm.clear(); }
                else if (comm == "TIMESTAMP=") { timedet    = true; comm.clear(); }
                else if (comm == "ISOSPEED=")  { isodet     = true; comm.clear(); }
                else if (comm == "MODEL=")     { instdet    = true; comm.clear(); }
                else if (comm == "CCD-TEMP=")  { ccdtempdet = true; comm.clear(); }
            } else if (ch > 32) {
                aline.push_back(static_cast<char>(ch));
            }
            
            if (ch == 0x0A) comment = false;
            
            // Loop termination mirrors:
            //   until ((comment=false) and (ord(ch)<=32)) or (reader_position>200)
            if ((!comment && ch <= 32) || reader.position > 200) break;
        }
        
        if (aline.size() > 1) {
            ++header_field;
            if      (header_field == 1) w1   = aline;
            else if (header_field == 2) h1   = aline;
            else                        bits = aline;
        }
    }
    
    int err = 0, err2 = 0, err3 = 0;
    head.width  = static_cast<int>(val_int(w1, err));
    head.height = static_cast<int>(val_int(h1, err2));
    const double range = val_double(bits, err3);
    nrbits = static_cast<int>(std::lround(range));
    
    if (pfm) {
        nrbits = -32;
        head.datamax_org = 0xFFFF;
    } else if (nrbits == 65535) {
        nrbits = 16; head.datamax_org = 0xFFFF;
    } else if (nrbits == 255) {
        nrbits = 8;  head.datamax_org = 0xFF;
    } else {
        err3 = 999;
    }
    
    if (err != 0 || err2 != 0 || err3 != 0) {
        show_error_label("Incompatible PPM/PGM/PFM file !!");
        head.naxis = 0;
        return false;
    }
    
    head.datamin_org = 0;
    bck.backgr = head.datamin_org;
    cwhite     = head.datamax_org;
    
    int package = 0;  // bytes per pixel
    if (color7) {
        package      = static_cast<int>(std::lround(std::abs(nrbits) * 3.0 / 8.0));
        head.naxis3  = 3;
        head.naxis   = 3;
    } else {
        package      = static_cast<int>(std::lround(std::abs(nrbits) / 8.0));
        head.naxis3  = 1;
        head.naxis   = 2;
    }
    
    // Allocate one row's worth on demand.
    img_loaded2.assign(head.naxis3,
                       std::vector<std::vector<float>>(head.height,
                                                       std::vector<float>(head.width, 0.0f)));
                                                       
    std::vector<std::uint8_t> row(static_cast<std::size_t>(head.width) * package);
    
    for (int i = 0; i < head.height; ++i) {
        const auto got = reader.read(row.data(), row.size());
        if (got < row.size()) {
            // Silently tolerate short reads (mirrors original behavior).
        }
        
        for (int j = 0; j < head.width; ++j) {
            if (!color7) {
                if (nrbits == 8) {
                    img_loaded2[0][i][j] = row[j];
                } else if (nrbits == 16) {
                    std::uint16_t be{};
                    std::memcpy(&be, &row[j * 2], 2);
                    img_loaded2[0][i][j] = swap16(be);
                } else if (pfm) {
                    std::uint32_t raw{};
                    std::memcpy(&raw, &row[j * 4], 4);
                    if (range < 0) {
                        // Little-endian floats — payload already native.
                        float v;
                        std::memcpy(&v, &raw, 4);
                        img_loaded2[0][i][j] = static_cast<float>(v * 65535.0 / (-range));
                    } else {
                        // Big-endian floats — byteswap and reinterpret.
                        const std::uint32_t le = std::byteswap(raw);
                        float v;
                        std::memcpy(&v, &le, 4);
                        img_loaded2[0][i][j] = static_cast<float>(v * 65535.0 / range);
                    }
                }
            } else {
                if (nrbits == 8) {
                    img_loaded2[0][i][j] = row[j * 3 + 0];
                    img_loaded2[1][i][j] = row[j * 3 + 1];
                    img_loaded2[2][i][j] = row[j * 3 + 2];
                } else if (nrbits == 16) {
                    std::uint16_t r{}, g{}, b{};
                    std::memcpy(&r, &row[j * 6 + 0], 2);
                    std::memcpy(&g, &row[j * 6 + 2], 2);
                    std::memcpy(&b, &row[j * 6 + 4], 2);
                    img_loaded2[0][i][j] = swap16(r);
                    img_loaded2[1][i][j] = swap16(g);
                    img_loaded2[2][i][j] = swap16(b);
                } else if (pfm) {
                    auto load_float = [&](int byte_off) -> float {
                        std::uint32_t raw{};
                        std::memcpy(&raw, &row[byte_off], 4);
                        if (range < 0) {
                            float v;
                            std::memcpy(&v, &raw, 4);
                            return static_cast<float>(v * 65535.0 / (-range));
                        }
                        const std::uint32_t le = std::byteswap(raw);
                        float v;
                        std::memcpy(&v, &le, 4);
                        return static_cast<float>(v * 65535.0 / range);
                    };
                    img_loaded2[0][i][j] = load_float(j * 12 + 0);
                    img_loaded2[1][i][j] = load_float(j * 12 + 4);
                    img_loaded2[2][i][j] = load_float(j * 12 + 8);
                }
            }
        }
    }
    
    unsaved_import = false;
    
    // ---- Synthesize a minimal FITS header in memo --------------------------
    for (int j = 0; j <= 10; ++j) {
        if (j == 5 && head.naxis3 == 1) continue;  // skip NAXIS3 for mono.
        if (j < static_cast<int>(head1.size())) memo.push_back(head1[j]);
    }
    if (27 < static_cast<int>(head1.size())) memo.push_back(head1[27]);  // END.
    
    update_integer(memo, "BITPIX  =", " / Bits per entry                                 ", nrbits);
    update_integer(memo, "NAXIS   =", " / Number of dimensions                           ", head.naxis);
    update_integer(memo, "NAXIS1  =", " / length of x axis                               ", head.width);
    update_integer(memo, "NAXIS2  =", " / length of y axis                               ", head.height);
    if (head.naxis3 != 1) {
        update_integer(memo, "NAXIS3  =", " / length of z axis (mostly colors)               ", head.naxis3);
    }
    update_integer(memo, "DATAMIN =", " / Minimum data value                             ", 0);
    update_integer(memo, "DATAMAX =", " / Maximum data value                           ",
                   static_cast<long long>(std::lround(head.datamax_org)));
                   
    if (head.exposure != 0) {
        update_float(memo, "EXPTIME =", " / duration of exposure in seconds                ", false, head.exposure);
    }
    if (!head.gain.empty()) {
        int e = 0;
        update_integer(memo, "GAIN    =", " / iso speed                                      ", val_int(head.gain, e));
    }
    
    if (!head.date_obs.empty()) update_text(memo, "DATE-OBS=", "'" + head.date_obs + "'");
    if (!instrum.empty())       update_text(memo, "INSTRUME=", "'" + instrum + "'");
    
    update_text(memo, "BAYERPAT=", "'T'                  / Unknown Bayer color pattern                  ");
    update_text(memo, "COMMENT 1", "  Written by ASTAP, Astrometric STAcking Program. www.hnsky.org");
    
    return true;
}

// ===========================================================================
// load_TIFFPNGJPEG
//
// The actual pixel decode is delegated to the installed IImageDecoder;
// dispatch + header population is faithful to the original.
// ===========================================================================
bool load_tiff_png_jpeg(const std::filesystem::path& filen,
                        bool                         light,
                        Header&                      head,
                        ImageArray&                  img,
                        std::vector<std::string>&    memo) {
    head.naxis = 0;
    
    [[maybe_unused]] auto tiff = false;
    [[maybe_unused]] auto png  = false;
    [[maybe_unused]] auto jpeg = false;
    auto bmp  = false;
    auto saved_header = false;
    
    const std::string ext = ext_upper(filen);
    
    if      (ext == ".TIF"  || ext == ".TIFF") tiff = true;
    else if (ext == ".PNG")                    png  = true;
    else if (ext == ".JPG"  || ext == ".JPEG") jpeg = true;
    else if (ext == ".BMP")                    bmp  = true;
    else                                       return false;
    
    // Ask the installed decoder to hand us a DecodedImage. The null default
    // returns false + an error string; hosts call set_image_decoder() with
    // their own stb_image / libpng / libjpeg / libtiff-backed implementation
    // at startup.
    IImageDecoder::DecodedImage decoded;
    std::string                 decode_err;
    if (!image_decoder().decode_raster(filen, ext, decoded, decode_err)) {
        memo.push_back("Image decode failed: " + decode_err);
        return false;
    }
    
    int  image_width  = decoded.width;
    int  image_height = decoded.height;
    bool colour       = (decoded.channels >= 3) || bmp;
    std::string descrip;  // TIFF ImageDescription — see note below.
    reset_fits_global_variables(light, head);
    
    if (!colour) { head.naxis = 2; head.naxis3 = 1; }
    else         { head.naxis = 3; head.naxis3 = 3; }
    
    memo.clear();
    
    extend_type      = 0;
    // Normalise to the 0..65535 convention used elsewhere in ASTAP.
    // 8-bit images are stored as 16-bit (scaled x256); 16-bit as-is.
    nrbits           = (decoded.bits_per_sample == 32) ? -32 : 16;
    head.datamin_org = 0;
    head.datamax_org = (nrbits == -32) ? 1.0 : 0xFFFF;
    bck.backgr       = head.datamin_org;
    cwhite           = head.datamax_org;
    
    head.width  = image_width;
    head.height = image_height;
    img.assign(head.naxis3,
               std::vector<std::vector<float>>(head.height,
                                               std::vector<float>(head.width, 0.0f)));
                                               
    // Copy decoded samples into `img`, flipping vertically to match the
    // FITS convention (origin bottom-left vs. raster top-left).
    const int  src_ch   = decoded.channels;
    const auto scale    = (decoded.bits_per_sample == 8)  ? 256.0f
                        : (decoded.bits_per_sample == 16) ? 1.0f
                        : 65535.0f;   // float sources renormalised to 0..65535.
    for (int c = 0; c < head.naxis3; ++c) {
        const int sc = (src_ch == 1) ? 0 : std::min(c, src_ch - 1);
        for (int y = 0; y < image_height; ++y) {
            const int src_y = image_height - 1 - y;
            const auto& src_row = decoded.pixels[sc][src_y];
            auto&       dst_row = img[c][y];
            for (int x = 0; x < image_width; ++x) {
                dst_row[x] = src_row[x] * scale;
            }
        }
    }
    
    // TIFF ImageDescription: some decoders surface tag metadata through an
    // implementation-specific channel. If the host's decoder populates
    // `decoded.pixels` metadata via a side channel it can set `descrip`
    // directly; left empty here, parallels the original when the tag is absent.
    
    if (descrip.size() >= 6 && descrip.compare(0, 6, "SIMPLE") == 0) {
        // FITS header embedded in TIFF description — split on \n into memo.
        std::stringstream ss(descrip);
        std::string       line;
        while (std::getline(ss, line)) memo.push_back(line);
        read_keys_memo(light, head, memo);
        saved_header = true;
    } else {
        for (int j = 0; j <= 10; ++j) {
            if (j == 5 && head.naxis3 == 1) continue;
            if (j < static_cast<int>(head1.size())) memo.push_back(head1[j]);
        }
        if (27 < static_cast<int>(head1.size())) memo.push_back(head1[27]);
        if (!descrip.empty()) add_long_comment(memo, descrip);
    }
    
    update_integer(memo, "BITPIX  =", " / Bits per entry                                 ", nrbits);
    update_integer(memo, "NAXIS   =", " / Number of dimensions                           ", head.naxis);
    update_integer(memo, "NAXIS1  =", " / length of x axis                               ", head.width);
    update_integer(memo, "NAXIS2  =", " / length of y axis                               ", head.height);
    update_integer(memo, "DATAMIN =", " / Minimum data value                             ", 0);
    update_integer(memo, "DATAMAX =", " / Maximum data value                             ",
                   static_cast<long long>(std::lround(head.datamax_org)));
                   
    if (!saved_header) {
        const double jd2 = file_age_jd(filen);
        head.date_obs = jd_to_date(jd2);
        update_text(memo, "DATE-OBS=", "'" + head.date_obs + "'");
    }
    
    update_text(memo, "COMMENT 1", "  Written by ASTAP, Astrometric STAcking Program. www.hnsky.org");
    
    unsaved_import = true;
    return true;
}

// ===========================================================================
// save_tiff16_secure
// ===========================================================================
bool save_tiff16_secure(const ImageArray&               img,
                        const std::vector<std::string>& memo,
                        const std::filesystem::path&    filen2) {
    auto tmp = filen2; tmp.replace_extension(".tmp");
    
    // Build the TIFF "ImageDescription" payload from the memo.
    std::string description;
    for (const auto& line : memo) { description += line; description.push_back('\n'); }
    
    if (!astap::image::save_tiff_16(img, tmp, description, false, false)) {
        return false;
    }
    
    std::error_code ec;
    std::filesystem::remove(filen2, ec);  // Ignored if missing.
    if (ec && std::filesystem::exists(filen2)) return false;
    
    std::filesystem::rename(tmp, filen2, ec);
    return !ec;
}

// ===========================================================================
// convert_raw
// ===========================================================================
bool convert_raw(bool         loadfile,
                 bool         savefile,
                 std::string& filename3,
                 Header&      head,
                 ImageArray&  img) {
    bool        result = true;
    const int   conv_index = raw_conversion_program_index;
    std::string filename4;
    std::string commando;
    std::string param;
    
    namespace fs = std::filesystem;
    
    // ---- LibRaw branch (conv_index <= 1) ----------------------------------
    if (conv_index <= 1) {
        param  = (conv_index == 1) ? "-i" : "-f";
        
#if defined(_WIN32)
        const auto exe = application_path / "unprocessed_raw.exe";
        if (!fs::exists(exe)) {
            result = false;
        } else {
            // Pass the long path through; let the platform layer handle quoting.
            const std::string cmd = exe.string() + " " + param + " \"" + filename3 + "\"";
            execute_and_wait(cmd, false);
            filename4 = filename3 + ".fits";
        }
#elif defined(__APPLE__)
        const auto exe = application_path / "unprocessed_raw";
        if (!fs::exists(exe)) {
            result = false;
        } else {
            execute_and_wait(exe.string() + " " + param + " \"" + filename3 + "\"", false);
            filename4 = filename3 + ".fits";
        }
#else  // Linux
        const auto local_exe = application_path / "unprocessed_raw-astap";
        if (!fs::exists(local_exe)) {
            if (!fs::exists("/usr/lib/libraw/unprocessed_raw")) {
                if (!fs::exists("/usr/bin/unprocessed_raw")) {
                    result = false;
                } else {
                    execute_and_wait("/usr/bin/unprocessed_raw \"" + filename3 + "\"", false);
                    filename4 = filename3 + ".pgm";
                }
            } else {
                execute_and_wait("/usr/lib/libraw/unprocessed_raw \"" + filename3 + "\"", false);
                filename4 = filename3 + ".pgm";
            }
        } else {
            execute_and_wait(local_exe.string() + " " + param + " \"" + filename3 + "\"", false);
            filename4 = filename3 + ".fits";
        }
#endif
    }
    
    // ---- DCRAW branch (conv_index == 2) -----------------------------------
    if (conv_index == 2) {
        if (ext_upper(filename3) == ".CR3") return false;  // dcraw can't process CR3.
        commando = "-D -4 -t 0";
        
#if defined(_WIN32)
        const auto exe = application_path / "dcraw.exe";
        if (!fs::exists(exe)) {
            result = false;
        } else {
            execute_and_wait(exe.string() + " " + commando + " \"" + filename3 + "\"", false);
        }
#elif defined(__APPLE__)
        const auto exe = application_path / "dcraw";
        if (!fs::exists(exe)) {
            result = false;
        } else {
            execute_and_wait(exe.string() + " " + commando + " \"" + filename3 + "\"", false);
        }
#else  // Linux
        const auto local_exe = application_path / "dcraw-astap";
        if (!fs::exists(local_exe)) {
            if (!fs::exists("/usr/bin/dcraw-astap")) {
                if (!fs::exists("/usr/local/bin/dcraw-astap")) {
                    if (!fs::exists("/usr/bin/dcraw")) {
                        if (!fs::exists("/usr/local/bin/dcraw")) {
                            result = false;
                        } else {
                            execute_and_wait("/usr/local/bin/dcraw " + commando + " \"" + filename3 + "\"", false);
                        }
                    } else {
                        execute_and_wait("/usr/bin/dcraw " + commando + " \"" + filename3 + "\"", false);
                    }
                } else {
                    execute_and_wait("/usr/local/bin/dcraw-astap " + commando + " \"" + filename3 + "\"", false);
                }
            } else {
                execute_and_wait("/usr/bin/dcraw-astap " + commando + " \"" + filename3 + "\"", false);
            }
        } else {
            execute_and_wait(local_exe.string() + " " + commando + " \"" + filename3 + "\"", false);
        }
#endif
        
        if (!result) {
            memo2_message("DCRAW executable not found! Will try unprocessed_raw as alternative.");
        } else {
            filename4 = change_ext(filename3, ".pgm");
        }
    }
    
    // ---- No conversion program found --------------------------------------
    if (!result) {
        // Delegate error reporting to the host UI hook.
        if (conv_index == 2) {
            memo2_message("Could not find dcraw executable.");
        }
        if (conv_index <= 1) {
            memo2_message("Could not find unprocessed_raw (LibRaw) executable.");
        }
        return false;
    }
    
    // ---- Post-processing the produced intermediate ------------------------
    // TODO: mainwindow_memo1_lines pointer was a GUI-owned memo. When no
    // GUI is attached the fall-through to memox handles the memo sink.
    extern std::vector<std::string>* mainwindow_memo1_lines;
    
    if (ext_upper(filename4) == ".PGM") {
        std::vector<std::string>& memo_target = mainwindow_memo1_lines ? *mainwindow_memo1_lines : *memox;
        
        if (load_ppm_pgm_pfm(filename4, head, img, memo_target)) {
            std::error_code ec;
            std::filesystem::remove(filename4, ec);  // delete temp pgm.
            filename4 = change_ext(filename4, ".fits");
            
            if (head.date_obs.empty()) {
                const double jd2 = file_age_jd(filename3);
                head.date_obs = jd_to_date(jd2);
                update_text(memo_target, "DATE-OBS=", "'" + head.date_obs + "'");
            }
            update_text(memo_target, "BAYERPAT=", "'\?\?\?\?'");
            add_text(memo_target, "HISTORY  ", std::string("Converted from ") + filename3);
            result = true;
        } else {
            result = false;
        }
        
        if (savefile && conv_index == 2 && result) {
            head.set_temperature = extract_temperature_from_filename(filename4);
            update_text(memo_target,
                        "OBJECT  =",
                        "'" + extract_objectname_from_filename(filename4) + "'");
            result = save_fits(img, memo_target, filename4, 16, true);
        }
        if (!loadfile) img.clear();
    } else {
        // FITS produced directly by the patched unprocessed_raw.
        if (loadfile) {
            std::vector<std::string>& memo_target = mainwindow_memo1_lines ? *mainwindow_memo1_lines : *memox;
            result = load_fits(filename4, true, true, true, 0, memo_target, head, img);
            if (result && !savefile) {
                std::error_code ec;
                std::filesystem::remove(filename4, ec);
                filename4 = change_ext(filename3, ".fits");
            }
        }
    }
    
    if (result) filename3 = filename4;
    return result;
}

// ===========================================================================
// convert_to_fits
// ===========================================================================
bool convert_to_fits(std::string& filen) {
    const std::string ext = ext_upper(filen);
    bool result = false;
    
    Header     headX;
    ImageArray img_temp;
    
    // `head` and `img_loaded` reach us via the using-decls at the top of
    // this TU (canonical in ../core/globals.h).
    
    if (check_raw_file_extension(ext)) {
        result = convert_raw(false, true, filen, head, img_loaded);
    } else if (ext == ".FZ") {
        result = unpack_cfitsio(filen);
    } else {
        if (ext == ".PPM" || ext == ".PGM" || ext == ".PFM" || ext == ".PBM") {
            result = load_ppm_pgm_pfm(filen, headX, img_temp, *memox);
        } else if (ext == ".XISF") {
            result = astap::image::load_xisf(filen, head, img_loaded, *memox);
        } else if (ext == ".JPG" || ext == ".JPEG" || ext == ".PNG" || ext == ".TIF" || ext == ".TIFF") {
            result = load_tiff_png_jpeg(filen, true, headX, img_temp, *memox);
        }
        
        if (result) {
            if (head.exposure == 0) {
                head.exposure = extract_exposure_from_filename(filen);
                update_text(*memox, "OBJECT  =",
                            "'" + extract_objectname_from_filename(filen) + "'");
                head.set_temperature = extract_temperature_from_filename(filen);
            }
            filen  = change_ext(filen, ".fits");
            result = save_fits(img_loaded, *memox, filen, nrbits, false);
        }
    }
    return result;
}

// ===========================================================================
// load_image
// ===========================================================================
bool load_image(std::string&                 filename2,
                ImageArray&                  img,
                Header&                      head,
                std::vector<std::string>&    memo,
                bool                         re_center,
                bool                         plot) {
    // GUI-side bookkeeping (caption / shape markers / updown reset) is
    // delegated to the host; we only retain the data-layer behavior.
    std::string filename_org;
    if (plot) filename_org = filename2;
    
    const std::string ext1 = ext_upper(filename2);
    bool result = false;
    
    if (ext1 == ".FIT" || ext1 == ".FITS" || ext1 == ".FTS" || ext1 == ".NEW" ||
        ext1 == ".WCS" || ext1 == ".AXY"  || ext1 == ".XYLS"|| ext1 == ".GSC" ||
        ext1 == ".BAK") {
        result = load_fits(filename2, true, true, true, 0, memo, head, img);
        if (head.naxis < 2) result = false;
    } else if (ext1 == ".FZ") {
        if (unpack_cfitsio(filename2)) {
            result = load_fits(filename2, true, true, true, 0, memo, head, img);
        }
    } else if (check_raw_file_extension(ext1)) {
        result = convert_raw(true, false, filename2, head, img);
        if (result) filename2 = change_ext(filename2, ".fits");
    } else if (ext1 == ".PPM" || ext1 == ".PGM" || ext1 == ".PFM" || ext1 == ".PBM") {
        result = load_ppm_pgm_pfm(filename2, head, img, memo);
    } else if (ext1 == ".XISF") {
        result = astap::image::load_xisf(filename2, head, img, memo);
    } else if (ext1 == ".TIF" || ext1 == ".TIFF" || ext1 == ".PNG" ||
               ext1 == ".JPG" || ext1 == ".JPEG" || ext1 == ".BMP") {
        result = load_tiff_png_jpeg(filename2, true, head, img, memo);
    }
    
    if (!result) {
        // Delegated to GUI layer (update_menu).
        return false;
    }
    
    if (plot) {
        // Demosaic / histogram / replot all live in the GUI/analysis layer;
        // the only data-layer side effect we keep is the recent-files entry.
        // TODO: hook demosaic_advanced + use_histogram + plot_fits when those
        // analysis modules are wired up.
        add_recent_file(filename_org, recent_files);
    }
    
    if (!commandline_execution) {
        // TODO: backup-buffer management lives in the undo module.
    }
    return true;
}

// ===========================================================================
// add_recent_file
// ===========================================================================
void add_recent_file(std::string_view              f,
                     std::vector<std::string>&     recent_files) {
    // Hoist any existing entry to the front.
    for (auto it = recent_files.begin(); it != recent_files.end(); ++it) {
        if (*it == f) { recent_files.erase(it); break; }
    }
    recent_files.insert(recent_files.begin(), std::string(f));
    if (recent_files.size() > 8) recent_files.resize(8);
    update_recent_file_menu();
}
    
} // namespace
