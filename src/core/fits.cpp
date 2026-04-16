///----------------------------------------
///      @file fits.cpp
///   @ingroup ASTAP++
///     @brief FITS I/O and header-edit primitives -- implementation.
///   @details Faithful port. Quirks of the original (e.g. validate_double
///            scanning past position 30 to tolerate CFITSIO violations of
///            FITS standard 4, BZERO special-cased for MaximDL's
///            signed-overflow bug) are preserved.
///    @author Ported from Han Kleijn's ASTAP; MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "fits.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstring>
#include <fstream>
#include <numbers>
#include <string>
#include <string_view>
#include <vector>

#include "ephemerides.h"   // ephem::precess_iau1976

///----------------------------------------
namespace astap::core {
///----------------------------------------

namespace {

// ----- TODO: living in core/util.h ------------------------------------------
// These helpers (floattostr4, str(x:20,s), strtoint2, etc.) should move to
// util.h; for now stub them inline so this file compiles standalone for review.
std::string floattostr20(double x) {
    // Right-aligned width-20 default float format.
    char buf[64];
    int n = std::snprintf(buf, sizeof(buf), "%20.10g", x);
    if (n < 0) return std::string{};
    return std::string(buf, static_cast<std::size_t>(n));
}
std::string inttostr20(long long x) {
    char buf[64];
    int n = std::snprintf(buf, sizeof(buf), "%20lld", x);
    if (n < 0) return std::string{};
    return std::string(buf, static_cast<std::size_t>(n));
}
std::string trim(std::string_view s) {
    auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string_view::npos) return std::string{};
    auto last = s.find_last_not_of(" \t\r\n");
    return std::string(s.substr(first, last - first + 1));
}

// ----- TODO: extern globals from astap_main ---------------------------------
// The C++ port will expose these from a global state header; declare here so
// the body of the port type-checks. Default-init only; the real definitions
// live elsewhere.
double ra_mount = 999, dec_mount = 999;
double ra_radians = 0, dec_radians = 0;
double focallen = 0;
int subsamp = 1;
int xbayroff = 0, ybayroff = 0;
std::string roworder;
double a_order = 0, ap_order = 0;
double a_0_0=0,a_0_1=0,a_0_2=0,a_0_3=0,a_1_0=0,a_1_1=0,a_1_2=0,a_2_0=0,a_2_1=0,a_3_0=0;
double b_0_0=0,b_0_1=0,b_0_2=0,b_0_3=0,b_1_0=0,b_1_1=0,b_1_2=0,b_2_0=0,b_2_1=0,b_3_0=0;
double ap_0_0=0,ap_0_1=0,ap_0_2=0,ap_0_3=0,ap_1_0=0,ap_1_1=0,ap_1_2=0,ap_2_0=0,ap_2_1=0,ap_3_0=0;
double bp_0_0=0,bp_0_1=0,bp_0_2=0,bp_0_3=0,bp_1_0=0,bp_1_1=0,bp_1_2=0,bp_2_0=0,bp_2_1=0,bp_3_0=0;
std::string centalt, centaz;
std::array<double, 20> x_coeff{}, y_coeff{};
std::array<double, 20> ppo_coeff{};
std::string telescop, instrum, origin, object_name;
std::string sitelat, sitelong, siteelev;
double focus_temp = 999;
int focus_pos = 0;
double pressure = 1010;
double airmass = 0;
bool annotated = false;
double site_lat_radians = 999;
std::string sqm_value;
double equinox = 2000;
std::string imagetype;
std::string bayerpat;
int nrbits = 0;
int extend_type = 0;  // 0 image / 1 image-ext / 2 ascii table / 3 bintable
bool last_extension = true;
bool unsaved_import = false;
double bandpass = 0;
double plate_ra = 0, plate_dec = 0;
int dec_sign = 1;
double x_pixel_size = 0, y_pixel_size = 0;
int x_pixel_offset = 0, y_pixel_offset = 0;
double cwhite = 0;
astap::Background bck{};
bool sip = false;
std::string filename2;
bool esc_pressed = false;
std::array<std::string, 1> head1{};  // placeholder for the saved primary card
std::string sqm_key = "SQM     ";  // 8 chars; adjustable in the original
astap::Header head{};               // global "current head" (not the one passed in)

// ----- TODO: GUI helpers replaced with stubs --------------------------------
void memo2_message(const std::string& /*msg*/) {
    // TODO: route to a logging sink; the original wrote to mainwindow.memo2.
}

// ----- low-level helpers ----------------------------------------------------

constexpr auto kPi = std::numbers::pi_v<double>;

template <class T>
T swap_endian(T v) {
    static_assert(std::is_trivially_copyable_v<T>);
    auto bytes = std::bit_cast<std::array<std::uint8_t, sizeof(T)>>(v);
    std::reverse(bytes.begin(), bytes.end());
    return std::bit_cast<T>(bytes);
}

// Swap the two bytes of a 16-bit value.
std::uint16_t swap16(std::uint16_t v) {
    return static_cast<std::uint16_t>((v >> 8) | (v << 8));
}

// JdToDate -- TODO: living in core/time.h. Return ISO-ish stub.
std::string jd_to_date(double jd) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "JD%.6f", jd);
    return buf;
}

// Precession lives in core/ephemerides.cpp as ephem::precess_iau1976.

// Format a card line padded to 80 characters (truncates if longer).
std::string pad80(std::string_view s) {
    std::string out(s);
    if (out.size() < kFitsCardSize) out.append(kFitsCardSize - out.size(), ' ');
    else if (out.size() > kFitsCardSize) out.resize(kFitsCardSize);
    return out;
}

// pos1: 1-based index of needle in haystack, 0 if not found.
std::size_t pos1(std::string_view needle, std::string_view hay) {
    auto p = hay.find(needle);
    return p == std::string_view::npos ? 0 : p + 1;
}

// Position of the line ending after a 1-based offset (returns 1-based pos, or
// hay.size()+1 if there is no line end).
std::size_t line_end_after(std::string_view hay, std::size_t one_based) {
    auto p = hay.find('\n', one_based - 1);
    if (p == std::string_view::npos) return hay.size() + 1;
    return p + 1;  // 1-based
}
    
}  // namespace

// ---------------------------------------------------------------------------
// reset_fits_global_variables
// ---------------------------------------------------------------------------

void reset_fits_global_variables(bool light, astap::Header& head_out) noexcept {
    if (light) {
        head_out.crota2 = 999;
        head_out.crota1 = 999;
        head_out.ra0 = 0;
        head_out.dec0 = 0;
        ra_mount = 999;
        dec_mount = 999;
        head_out.cdelt1 = 0;
        head_out.cdelt2 = 0;
        head_out.xpixsz = 0;
        head_out.ypixsz = 0;
        focallen = 0;
        subsamp = 1;
        head_out.cd1_1 = 0;
        head_out.cd1_2 = 0;
        head_out.cd2_1 = 0;
        head_out.cd2_2 = 0;
        xbayroff = 0;
        ybayroff = 0;
        roworder.clear();
        
        a_order = 0;
        ap_order = 0;
        a_0_0=a_0_1=a_0_2=a_0_3=a_1_0=a_1_1=a_1_2=a_2_0=a_2_1=a_3_0=0;
        b_0_0=b_0_1=b_0_2=b_0_3=b_1_0=b_1_1=b_1_2=b_2_0=b_2_1=b_3_0=0;
        ap_0_0=ap_0_1=ap_0_2=ap_0_3=ap_1_0=ap_1_1=ap_1_2=ap_2_0=ap_2_1=ap_3_0=0;
        bp_0_0=bp_0_1=bp_0_2=bp_0_3=bp_1_0=bp_1_1=bp_1_2=bp_2_0=bp_2_1=bp_3_0=0;
        
        centalt.clear();
        centaz.clear();
        
        x_coeff[0] = 0;
        y_coeff[0] = 0;
        
        head_out.xbinning = 1;
        head_out.ybinning = 1;
        head_out.mzero = 0;
        head_out.mzero_radius = 99;
        head_out.magn_limit = 0;
        head_out.pedestal = 0;
        
        telescop.clear(); instrum.clear(); origin.clear(); object_name.clear();
        sitelat.clear(); sitelong.clear(); siteelev.clear();
        
        focus_temp = 999;
        focus_pos = 0;
        pressure = 1010;
        airmass = 0;
        annotated = false;
        site_lat_radians = 999;
        
        sqm_value.clear();
        equinox = 2000;
    }
    
    head_out.date_obs.clear();
    head_out.date_avg.clear();
    head_out.calstat.clear();
    head_out.filter_name = "CV";
    head_out.naxis = -1;
    head_out.naxis3 = 1;
    head_out.datamin_org = 0;
    imagetype.clear();
    head_out.exposure = 0;
    head_out.set_temperature = 999;
    head_out.gain.clear();
    head_out.egain.clear();
    head_out.passband_database.clear();
    bayerpat.clear();
    head_out.issues.clear();
}

// ---------------------------------------------------------------------------
// load_fits
// ---------------------------------------------------------------------------

namespace {

// Nested helpers: validate_double / get_string / get_as_string.
// `header` is the 2880-byte block, `i` the offset of the start of the current
// 80-byte card.
double validate_double(const std::array<char, 2881>& header, int i, int& err) {
    // Read positions i+10..i+30, stop at '/'.
    std::string t;
    int r = i + 10;
    while (header[r] != '/' && r <= i + 30) {
        if (header[r] != ' ') t.push_back(header[r]);
        ++r;
    }
    if (t.empty()) { err = 1; return 0.0; }
    try {
        std::size_t consumed = 0;
        double v = std::stod(t, &consumed);
        err = (consumed == t.size()) ? 0 : 1;
        return v;
    } catch (...) {
        err = 1;
        return 0.0;
    }
}

std::string get_string(const std::array<char, 2881>& header, int i) {
    // Read between single quotes starting near position i+11.
    std::string out;
    int r = i + 11;
    while (header[r - 1] != '\'' && r < i + 77) ++r;
    while (header[r] != '\'' && r < i + 79) {
        out.push_back(header[r]);
        ++r;
    }
    // trim
    auto start = out.find_first_not_of(' ');
    if (start == std::string::npos) return std::string{};
    auto end = out.find_last_not_of(' ');
    return out.substr(start, end - start + 1);
}

std::string get_as_string(const std::array<char, 2881>& header, int i) {
    // Universal: handles both string-quoted and bare numeric values.
    std::string out;
    char c = header[i + 10];
    if (c != '\'') out.push_back(c);
    int r = i + 11;
    while (header[r] != '\'' && r < i + 30) {
        out.push_back(header[r]);
        ++r;
    }
    return out;
}

bool starts_with(const std::array<char, 2881>& h, int off, std::string_view s) {
    for (std::size_t k = 0; k < s.size(); ++k)
        if (h[off + static_cast<int>(k)] != s[k]) return false;
    return true;
}
    
}  // namespace

bool load_fits(const std::filesystem::path& filen,
               bool light,
               bool load_data,
               bool update_memo,
               int get_ext,
               std::vector<std::string>& memo,
               astap::Header& head_out,
               astap::ImageArray& img_loaded2) {
    bool result = false;
    img_loaded2.clear();
    
    // TODO: mainwindow.caption := ExtractFileName(filen) -- GUI only.
    
    // TODO: tiff_file_name(filen) -> delegate to load_TIFFPNGJPEG.
    // Astro-TIFF short-circuit is skipped here; TIFF path lives in
    // image/tiff.cpp.
    
    std::ifstream the_file(filen, std::ios::binary);
    if (!the_file) {
        // TODO: mainwindow.error_label1 -- GUI only.
        return false;
    }
    
    the_file.seekg(0, std::ios::end);
    const auto file_size = static_cast<std::int64_t>(the_file.tellg());
    the_file.seekg(0, std::ios::beg);
    
    if (update_memo) memo.clear();
    
    reset_fits_global_variables(light, head_out);
    
    if (get_ext == 0) extend_type = 0;
    int naxis1 = 0;
    long long bzero = 0;  // MaximDL overflow workaround
    float bscale = 1.0f;
    double ccd_temperature = 999.0;
    float measured_max = 0;
    double pc1_1 = 0, pc1_2 = 0, pc2_1 = 0, pc2_2 = 0;
    
    int header_count = 0;
    bool bintable = false, asciitable = false;
    bool simple = false, image = false;
    int reader_position = 0;
    
    std::array<char, 2881> header{};
    
    bool end_record = false;
    
    // Outer loop over 2880-byte header blocks until END found.
    do {
        int i = 0;
        
        // Inner loop: read blocks until we land on the requested extension.
        for (;;) {
            the_file.read(header.data(), kFitsBlockSize);
            if (!the_file) {
                // TODO: error_label
                return false;
            }
            reader_position += static_cast<int>(kFitsBlockSize);
            
            if (reader_position == static_cast<int>(kFitsBlockSize) &&
                starts_with(header, 0, "SIMPLE ")) {
                simple = true;
                image = true;
            }
            if (!simple) {
                // TODO: error_label, "SIMPLE not found".
                return false;
            }
            
            if (header_count < get_ext && starts_with(header, 0, "XTENSION=")) {
                ++header_count;
                image = starts_with(header, 11, "IMAGE ");
                bintable = starts_with(header, 11, "BINTAB");
                asciitable = starts_with(header, 11, "TABLE ");
                std::string ext_name = get_string(header, 0);
                if (ext_name.find("BINTABLE") != std::string::npos) extend_type = 3;
                else if (ext_name == "TABLE") extend_type = 2;
                else extend_type = 1;
            }
            
            if (simple && header_count >= get_ext) break;
        }
        
        // Walk the 36 cards in this 2880 block.
        std::vector<std::string> tform, ttype, tunit;
        std::vector<int> tform_nr, tbcol;
        int tfields = 0;
        
        do {
            if (load_data) {
                std::string aline(&header[i], kFitsCardSize);
                if (update_memo) memo.push_back(aline);
            }
            
            int err = 0;
            auto vd = [&](int& e) { return validate_double(header, i, e); };
            
            // NAXIS / NAXIS1 / NAXIS2 / NAXIS3
            if (starts_with(header, i, "NAXIS")) {
                if (header[i + 5] == ' ')
                    head_out.naxis = static_cast<int>(std::lround(vd(err)));
                else if (header[i + 5] == '1') {
                    naxis1 = static_cast<int>(std::lround(vd(err)));
                    head_out.width = naxis1;
                } else if (header[i + 5] == '2')
                    head_out.height = static_cast<int>(std::lround(vd(err)));
                else if (header[i + 5] == '3') {
                    head_out.naxis3 = static_cast<int>(std::lround(vd(err)));
                    if (head_out.naxis == 3 && naxis1 == 3) {
                        head_out.width = head_out.height;
                        head_out.height = head_out.naxis3;
                        head_out.naxis3 = 1;
                    }
                    if (head_out.naxis3 > 3) {
                        head_out.naxis3 = 1;
                        memo2_message(
                            "Warning more then three colours. Displayed only the first one.");
                    }
                }
            }
            
            if (image) {
                // BITPIX, BZERO, BSCALE, BAYERPAT, BIAS_CNT
                if (header[i] == 'B') {
                    if (starts_with(header, i + 1, "AYERPA"))
                        bayerpat = get_string(header, i);
                    else if (starts_with(header, i + 1, "ITPIX"))
                        nrbits = static_cast<int>(std::lround(vd(err)));
                    else if (starts_with(header, i + 1, "ZERO")) {
                        double tempval = vd(err);
                        if (tempval > 2147483647.0) bzero = -2147483648LL;
                        else bzero = static_cast<long long>(std::llround(tempval));
                    } else if (starts_with(header, i + 1, "SCAL"))
                        bscale = static_cast<float>(vd(err));
                    else if (starts_with(header, i + 1, "IAS_CNT"))
                        head_out.flatdark_count = static_cast<int>(std::lround(vd(err)));
                }
                
                if (header[i] == 'C') {
                    if (starts_with(header, i + 1, "ALSTA"))
                        head_out.calstat = get_string(header, i);
                    else if (starts_with(header, i + 1, "CD-TEM"))
                        ccd_temperature = vd(err);
                }
                
                if (header[i] == 'E') {
                    if (starts_with(header, i + 1, "GAIN"))
                        head_out.egain = trim(get_as_string(header, i));
                    else if (starts_with(header, i + 1, "XP")) {
                        if (starts_with(header, i + 3, "OSUR"))
                            head_out.exposure = vd(err);
                        if (starts_with(header, i + 3, "TIME "))
                            head_out.exposure = vd(err);
                    }
                }
                
                if (starts_with(header, i, "SET-TEM")) {
                    try {
                        head_out.set_temperature =
                            static_cast<int>(std::lround(vd(err)));
                    } catch (...) {}
                }
                
                if (header[i] == 'I') {
                    if (starts_with(header, i + 1, "MAGETY"))
                        imagetype = get_string(header, i);
                    if (starts_with(header, i + 1, "SSUES"))
                        head_out.issues = get_string(header, i);
                }
                
                if (header[i] == 'F') {
                    if (starts_with(header, i + 1, "ILTER "))
                        head_out.filter_name = get_string(header, i);
                    else if (starts_with(header, i + 1, "LAT_CNT"))
                        head_out.flat_count = static_cast<int>(std::lround(vd(err)));
                }
                
                if (starts_with(header, i, "XBINNI"))
                    head_out.xbinning = vd(err);
                if (starts_with(header, i, "YBINNI"))
                    head_out.ybinning = vd(err);
                    
                if (starts_with(header, i, "GAIN "))
                    head_out.gain = trim(get_as_string(header, i));
                if (starts_with(header, i, "ISOSP"))
                    if (head_out.gain.empty())
                        head_out.gain = trim(get_as_string(header, i));
                        
                if (starts_with(header, i, "LIGH_CNT"))
                    head_out.light_count = static_cast<int>(std::lround(vd(err)));
                    
                if (starts_with(header, i, "TIME-OB")) {
                    if (head_out.date_obs.size() == 10)
                        head_out.date_obs += "T" + get_string(header, i);
                }
                
                if (header[i] == 'J' && header[i + 1] == 'D') {
                    if (header[i + 2] == ' ' && header[i + 3] == ' ' &&
                        header[i + 4] == ' ') {
                        if (head_out.date_obs.empty()) {
                            double jd2 = vd(err);
                            head_out.date_obs = jd_to_date(jd2);
                        }
                    } else if (header[i + 2] == '-' && header[i + 3] == 'A' &&
                               header[i + 4] == 'G') {
                        if (head_out.date_avg.empty()) {
                            double jd2 = vd(err);
                            head_out.date_avg = jd_to_date(jd2);
                        }
                    }
                }
                
                if (header[i] == 'D' && header[i + 1] == 'A') {
                    if (header[i + 2] == 'T' && header[i + 3] == 'E' &&
                        header[i + 4] == '-') {
                        if (header[i + 5] == 'O' && header[i + 6] == 'B')
                            head_out.date_obs = get_string(header, i);
                        else if (header[i + 5] == 'A' && header[i + 6] == 'V')
                            head_out.date_avg = get_string(header, i);
                    } else if (starts_with(header, i + 2, "RK_CNT"))
                        head_out.dark_count = static_cast<int>(std::lround(vd(err)));
                }
                
                if (light) {
                    // The full WCS / SIP / mount-conditions block. Heavy and
                    // repetitive; kept inline so callers don't have to index
                    // by keyword name.
                    if (header[i] == 'A') {
                        if (starts_with(header, i + 1, "MB-TEM"))
                            focus_temp = vd(err);
                        else if (starts_with(header, i + 1, "OC")) {
                            if (starts_with(header, i + 3, "BARO"))
                                pressure = vd(err);
                            else if (starts_with(header, i + 3, "AMBT"))
                                focus_temp = vd(err);
                        } else if (starts_with(header, i + 1, "NNOTAT"))
                            annotated = true;
                            
                        if (header[i + 1] == 'M' && header[i + 2] == 'D') {
                            // AMDX / AMDY index
                            if (header[i + 3] == 'X' || header[i + 3] == 'Y') {
                                std::string s;
                                s.push_back(header[i + 4]);
                                if (header[i + 5] != ' ') s.push_back(header[i + 5]);
                                int nr = std::atoi(s.c_str());
                                double v = vd(err);
                                if (nr >= 1 && nr <= 20) {
                                    if (header[i + 3] == 'X') x_coeff[nr - 1] = v;
                                    else y_coeff[nr - 1] = v;
                                }
                            }
                        } else if (starts_with(header, i + 1, "IRMAS"))
                            airmass = vd(err);
                            
                        if (header[i + 1] == '_') {
                            // SIP A_*
                            if (starts_with(header, i + 2, "ORD"))
                                a_order = static_cast<int>(std::lround(vd(err)));
                            #define SIP_AB(name, ch1, ch2, ch3) \
                                if (header[i+2]==(ch1) && header[i+3]==(ch2) && header[i+4]==(ch3)) name = vd(err);
                            SIP_AB(a_0_0, '0','_','0')
                            SIP_AB(a_0_1, '0','_','1')
                            SIP_AB(a_0_2, '0','_','2')
                            SIP_AB(a_0_3, '0','_','3')
                            SIP_AB(a_1_0, '1','_','0')
                            SIP_AB(a_1_1, '1','_','1')
                            SIP_AB(a_1_2, '1','_','2')
                            SIP_AB(a_2_0, '2','_','0')
                            SIP_AB(a_2_1, '2','_','1')
                            SIP_AB(a_3_0, '3','_','0')
                            #undef SIP_AB
                        }
                        if (header[i + 1] == 'P' && header[i + 2] == '_') {
                            if (starts_with(header, i + 3, "ORD"))
                                ap_order = static_cast<int>(std::lround(vd(err)));
                            #define SIP_AP(name, ch1, ch2, ch3) \
                                if (header[i+3]==(ch1) && header[i+4]==(ch2) && header[i+5]==(ch3)) name = vd(err);
                            SIP_AP(ap_0_0,'0','_','0')
                            SIP_AP(ap_0_1,'0','_','1')
                            SIP_AP(ap_0_2,'0','_','2')
                            SIP_AP(ap_0_3,'0','_','3')
                            SIP_AP(ap_1_0,'1','_','0')
                            SIP_AP(ap_1_1,'1','_','1')
                            SIP_AP(ap_1_2,'1','_','2')
                            SIP_AP(ap_2_0,'2','_','0')
                            SIP_AP(ap_2_1,'2','_','1')
                            SIP_AP(ap_3_0,'3','_','0')
                            #undef SIP_AP
                        }
                    }
                    
                    if (header[i] == 'B') {
                        if (starts_with(header, i + 1, "ANDPAS")) {
                            bandpass = vd(err);
                            if (bandpass == 35 || bandpass == 8)
                                head_out.filter_name = "red";
                            else if (bandpass == 18 || bandpass == 7)
                                head_out.filter_name = "blue";
                            else
                                head_out.filter_name = std::to_string(bandpass);
                        }
                        if (header[i + 1] == '_') {
                            #define SIP_B(name, ch1, ch2, ch3) \
                                if (header[i+2]==(ch1) && header[i+3]==(ch2) && header[i+4]==(ch3)) name = vd(err);
                            SIP_B(b_0_0,'0','_','0')
                            SIP_B(b_0_1,'0','_','1')
                            SIP_B(b_0_2,'0','_','2')
                            SIP_B(b_0_3,'0','_','3')
                            SIP_B(b_1_0,'1','_','o')   // Original has lowercase 'o' here -- preserved.
                            SIP_B(b_1_1,'1','_','1')
                            SIP_B(b_1_2,'1','_','2')
                            SIP_B(b_2_0,'2','_','0')
                            SIP_B(b_2_1,'2','_','1')
                            SIP_B(b_3_0,'3','_','0')
                            #undef SIP_B
                        }
                        if (header[i + 1] == 'P' && header[i + 2] == '_') {
                            #define SIP_BP(name, ch1, ch2, ch3) \
                                if (header[i+3]==(ch1) && header[i+4]==(ch2) && header[i+5]==(ch3)) name = vd(err);
                            SIP_BP(bp_0_0,'0','_','0')
                            SIP_BP(bp_0_1,'0','_','1')
                            SIP_BP(bp_0_2,'0','_','2')
                            SIP_BP(bp_0_3,'0','_','3')
                            SIP_BP(bp_1_0,'1','_','0')
                            SIP_BP(bp_1_1,'1','_','1')
                            SIP_BP(bp_1_2,'1','_','2')
                            SIP_BP(bp_2_0,'2','_','0')
                            SIP_BP(bp_2_1,'2','_','1')
                            SIP_BP(bp_3_0,'3','_','0')
                            #undef SIP_BP
                        }
                    }
                    
                    if (header[i] == 'C') {
                        if (header[i + 1] == 'R') {
                            if (starts_with(header, i + 2, "OTA")) {
                                if (header[i + 5] == '2') head_out.crota2 = vd(err);
                                else if (header[i + 5] == '1') head_out.crota1 = vd(err);
                            } else if (starts_with(header, i + 2, "PIX")) {
                                if (header[i + 5] == '1') head_out.crpix1 = vd(err);
                                else if (header[i + 5] == '2') head_out.crpix2 = vd(err);
                            }
                        } else if (starts_with(header, i + 1, "DELT")) {
                            if (header[i + 5] == '1') head_out.cdelt1 = vd(err);
                            else if (header[i + 5] == '2') head_out.cdelt2 = vd(err);
                        }
                        if (starts_with(header, i + 1, "RVAL")) {
                            if (header[i + 5] == '1') head_out.ra0 = vd(err) * kPi / 180;
                            if (header[i + 5] == '2') head_out.dec0 = vd(err) * kPi / 180;
                        } else if (header[i + 1] == 'D') {
                            if (header[i + 2] == '1' && header[i + 3] == '_' &&
                                header[i + 4] == '1') head_out.cd1_1 = vd(err);
                            if (header[i + 2] == '1' && header[i + 3] == '_' &&
                                header[i + 4] == '2') head_out.cd1_2 = vd(err);
                            if (header[i + 2] == '2' && header[i + 3] == '_' &&
                                header[i + 4] == '1') head_out.cd2_1 = vd(err);
                            if (header[i + 2] == '2' && header[i + 3] == '_' &&
                                header[i + 4] == '2') head_out.cd2_2 = vd(err);
                        } else if (starts_with(header, i + 1, "ENT")) {
                            if (starts_with(header, i + 4, "ALT"))
                                centalt = get_as_string(header, i);
                            else if (starts_with(header, i + 4, "AZ"))
                                centaz = get_as_string(header, i);
                        }
                        if (starts_with(header, i + 1, "NPIX")) {
                            if (header[i + 5] == '1')
                                x_pixel_offset = static_cast<int>(std::lround(vd(err)));
                            else if (header[i + 5] == '2')
                                y_pixel_offset = static_cast<int>(std::lround(vd(err)));
                        }
                    }
                    
                    if (header[i] == 'D' && header[i + 1] == 'E' &&
                        header[i + 2] == 'C' && header[i + 3] == ' ') {
                        int derr = 0;
                        double tempval = validate_double(header, i, derr) * kPi / 180;
                        if (derr == 0) {
                            dec_mount = tempval;
                            if (head_out.dec0 == 0) head_out.dec0 = tempval;
                        }
                    }
                    
                    if (header[i] == 'E') {
                        if (starts_with(header, i + 1, "QUINOX"))
                            equinox = vd(err);
                        if (starts_with(header, i + 1, "XTEND")) {
                            if (get_as_string(header, i).find('T') != std::string::npos)
                                last_extension = false;
                        }
                    }
                    
                    if (starts_with(header, i, "SECPIX") ||
                        starts_with(header, i, "SCALE ") ||
                        starts_with(header, i, "PIXSCA")) {
                        if (head_out.cdelt2 == 0) {
                            head_out.cdelt2 = vd(err) / 3600.0;
                            head_out.cdelt1 = head_out.cdelt2;
                        }
                    }
                    
                    if (starts_with(header, i, "FOC")) {
                        if (starts_with(header, i + 3, "ALL"))
                            focallen = vd(err);
                        else if (starts_with(header, i + 3, "USPO") ||
                                 starts_with(header, i + 3, "POS "))
                            try {
                                focus_pos = static_cast<int>(std::lround(vd(err)));
                            } catch (...) {}
                        else if (starts_with(header, i + 3, "USTE") ||
                                 starts_with(header, i + 3, "TEMP"))
                            focus_temp = vd(err);
                    }
                    
                    if (starts_with(header, i, "INSTRUM"))
                        instrum = get_string(header, i);
                        
                    if (starts_with(header, i, "MZERO")) {
                        if (header[i + 5] == 'R') head_out.mzero = vd(err);
                        if (header[i + 5] == 'A') head_out.mzero_radius = vd(err);
                        if (header[i + 5] == 'P')
                            head_out.passband_database = get_string(header, i);
                    }
                    
                    if (header[i] == 'O') {
                        if (starts_with(header, i + 1, "BS")) {
                            if (starts_with(header, i + 3, "LAT") ||
                                (header[i + 3] == '-' && starts_with(header, i + 4, "LA")))
                                sitelat = get_as_string(header, i);
                            if (starts_with(header, i + 3, "LON") ||
                                (header[i + 3] == '-' && starts_with(header, i + 4, "LO")))
                                sitelong = get_as_string(header, i);
                            if (starts_with(header, i + 3, "GEO-")) {
                                if (header[i + 7] == 'B') sitelat = get_as_string(header, i);
                                else if (header[i + 7] == 'L') sitelong = get_as_string(header, i);
                            }
                        }
                        if (starts_with(header, i + 1, "RIGIN"))
                            origin = get_string(header, i);
                        if (starts_with(header, i + 1, "BJ")) {
                            if (header[i + 3] == 'C' && header[i + 4] == 'T') {
                                if (header[i + 5] == 'R' && header[i + 6] == 'A' &&
                                    ra_mount >= 999) {
                                    // TODO: mainwindow.ra1.text := ...
                                    ra_mount = ra_radians;
                                } else if (header[i + 5] == 'D' && header[i + 6] == 'E' &&
                                           dec_mount >= 999) {
                                    // TODO: mainwindow.dec1.text := ...
                                    dec_mount = dec_radians;
                                } else if (header[i + 5] == 'A' && header[i + 6] == 'L' &&
                                           centalt.empty())
                                    centalt = get_as_string(header, i);
                                else if (header[i + 5] == 'A' && header[i + 6] == 'Z' &&
                                         centaz.empty())
                                    centaz = get_as_string(header, i);
                            } else if (header[i + 3] == 'E' && header[i + 4] == 'C' &&
                                       header[i + 5] == 'T')
                                object_name = get_string(header, i);
                        }
                    }
                    
                    if (header[i] == 'P') {
                        if (starts_with(header, i + 1, "RESSUR"))
                            pressure = vd(err);
                        if (starts_with(header, i + 1, "EDESTAL"))
                            head_out.pedestal = std::abs(vd(err));
                        if (header[i + 1] == 'L' && header[i + 2] == 'T') {
                            if (header[i + 3] == 'R' && header[i + 4] == 'A') {
                                if (header[i + 5] == 'H')
                                    plate_ra = vd(err) * kPi / 12;
                                if (header[i + 5] == 'M')
                                    plate_ra += vd(err) * kPi / (60 * 12);
                                if (header[i + 5] == 'S')
                                    plate_ra += vd(err) * kPi / (60 * 60 * 12);
                            } else if (header[i + 3] == 'D' && header[i + 4] == 'E') {
                                if (header[i + 7] == 'N')
                                    dec_sign = (header[i + 11] == '-') ? -1 : +1;
                                if (header[i + 6] == 'D')
                                    plate_dec = vd(err) * kPi / 180;
                                if (header[i + 6] == 'M')
                                    plate_dec += vd(err) * kPi / (60 * 180);
                                if (header[i + 6] == 'S')
                                    plate_dec = dec_sign * (plate_dec + vd(err) * kPi / (60 * 60 * 180));
                            }
                        }
                        if (header[i + 1] == 'P' && header[i + 2] == 'O') {
                            if (header[i + 3] == '3') ppo_coeff[2] = vd(err);
                            if (header[i + 3] == '6') ppo_coeff[5] = vd(err);
                        }
                        if (header[i + 1] == 'C') {
                            if (header[i + 2] == '1' && header[i + 3] == '_' && header[i + 4] == '1') pc1_1 = vd(err);
                            if (header[i + 2] == '1' && header[i + 3] == '_' && header[i + 4] == '2') pc1_2 = vd(err);
                            if (header[i + 2] == '2' && header[i + 3] == '_' && header[i + 4] == '1') pc2_1 = vd(err);
                            if (header[i + 2] == '2' && header[i + 3] == '_' && header[i + 4] == '2') pc2_2 = vd(err);
                        }
                    }
                    
                    if (header[i] == 'R') {
                        if (header[i + 1] == 'A' && header[i + 2] == ' ') {
                            int derr = 0;
                            double tempval = validate_double(header, i, derr) * kPi / 180;
                            if (derr == 0) {
                                ra_mount = tempval;
                                if (head_out.ra0 == 0) head_out.ra0 = tempval;
                            }
                        }
                        if (starts_with(header, i + 1, "OWORDE"))
                            roworder = get_string(header, i);
                    }
                    
                    if (header[i] == 'S') {
                        if (starts_with(header, i + 1, "ITE")) {
                            if (starts_with(header, i + 4, "LAT")) sitelat = get_as_string(header, i);
                            if (starts_with(header, i + 4, "LON")) sitelong = get_as_string(header, i);
                            if (starts_with(header, i + 4, "ELE")) siteelev = get_as_string(header, i);
                        }
                        if (starts_with(header, i + 1, "UBSAM"))
                            subsamp = static_cast<int>(std::lround(vd(err)));
                    }
                    
                    if (starts_with(header, i, "TELESCO"))
                        telescop = get_string(header, i);
                        
                    // Adjustable SQM-style keyword.
                    if (sqm_key.size() == 8 &&
                        std::memcmp(&header[i], sqm_key.data(), 8) == 0) {
                        sqm_value = trim(get_as_string(header, i));
                    }
                    
                    if (header[i] == 'X') {
                        if (starts_with(header, i + 1, "PIXEL"))
                            x_pixel_size = vd(err);
                        if (starts_with(header, i + 1, "PIXSZ"))
                            head_out.xpixsz = vd(err);
                        if (starts_with(header, i + 1, "BAYROF"))
                            xbayroff = static_cast<int>(std::lround(vd(err)));
                    }
                    if (header[i] == 'Y') {
                        if (starts_with(header, i + 1, "PIXEL"))
                            y_pixel_size = vd(err);
                        if (starts_with(header, i + 1, "PIXSZ"))
                            head_out.ypixsz = vd(err);
                        if (starts_with(header, i + 1, "BAYROF"))
                            ybayroff = static_cast<int>(std::lround(vd(err)));
                    }
                }  // light
            }    // image header
            else {
                // ---- table header --------------------------------------
                // BINTABLE / TABLE; build column descriptors so we can later
                // dump the table to the (skipped) memo3 sink.
                if (starts_with(header, i, "TFIELDS")) {
                    tfields = static_cast<int>(std::lround(vd(err)));
                    ttype.assign(tfields, "");
                    tform.assign(tfields, "");
                    tform_nr.assign(tfields, 0);
                    tbcol.assign(tfields, 0);
                    tunit.assign(tfields, "");
                }
                if (starts_with(header, i, "ZCMPTY")) last_extension = true;
                
                auto read_index = [&]() -> int {
                    std::string number;
                    number.push_back(header[i + 5]);
                    number.push_back(header[i + 6]);
                    number.push_back(header[i + 7]);
                    return std::atoi(number.c_str()) - 1;
                };
                
                if (starts_with(header, i, "TFORM")) {
                    int idx = read_index();
                    if (idx >= 0 && idx < tfields) {
                        tform[idx] = get_string(header, i);
                        // letter-prefix logic compressed: pick first matching
                        // type letter, parse the leading count.
                        for (char letter : std::string("EDLXBIJKA")) {
                            auto p = tform[idx].find(letter);
                            if (p != std::string::npos) {
                                std::string aline = trim(tform[idx]);
                                std::string count_s = aline.substr(0, p);
                                tform[idx] = std::string(1, letter);
                                tform_nr[idx] = std::max(1, std::atoi(("0" + count_s).c_str()));
                                break;
                            }
                        }
                    }
                }
                if (starts_with(header, i, "TBCOL")) {
                    int idx = read_index();
                    if (idx >= 0 && idx < tfields)
                        tbcol[idx] = static_cast<int>(std::lround(vd(err)));
                }
                if (starts_with(header, i, "TTYPE")) {
                    int idx = read_index();
                    if (idx >= 0 && idx < tfields) ttype[idx] = get_string(header, i);
                }
                if (starts_with(header, i, "TUNIT")) {
                    int idx = read_index();
                    if (idx >= 0 && idx < tfields) tunit[idx] = get_string(header, i);
                }
            }
            
            end_record =
                header[i] == 'E' && header[i + 1] == 'N' &&
                header[i + 2] == 'D' && header[i + 3] == ' ';
            i += 80;
        } while (i < static_cast<int>(kFitsBlockSize) && !end_record);
    } while (!end_record);
    
    if (head_out.naxis < 2) {
        if (head_out.naxis == 0) result = true;  // wcs file
        // TODO: mainwindow.image1.visible := false
        image = false;
    }
    
    if (image) {
        // RGB FITS NAXIS=3 + NAXIS1=3 -> handled as 2D 24-bit packed.
        if (head_out.naxis == 3 && naxis1 == 3) {
            nrbits = 24;
            head_out.naxis3 = 3;
        }
        
        if (light) {
            if (head_out.cd1_1 != 0 &&
                (head_out.cdelt1 == 0 || head_out.crota2 >= 999))
                new_to_old_WCS(head_out);
            else if (head_out.cd1_1 == 0 && head_out.cdelt2 != 0) {
                if (pc1_1 != 0) {
                    head_out.cd1_1 = pc1_1 * head_out.cdelt1;
                    head_out.cd1_2 = pc1_2 * head_out.cdelt1;
                    head_out.cd2_1 = pc2_1 * head_out.cdelt2;
                    head_out.cd2_2 = pc2_2 * head_out.cdelt2;
                    new_to_old_WCS(head_out);
                } else if (head_out.crota2 < 999) {
                    if (head_out.crota1 == 999) head_out.crota1 = head_out.crota2;
                    old_to_new_WCS(head_out);
                }
            }
            
            if (head_out.cd1_1 == 0 && head_out.cdelt2 == 0) {
                if (focallen != 0 && head_out.xpixsz != 0)
                    head_out.cdelt2 =
                        180.0 / (kPi * 1000) * head_out.xpixsz / focallen;
            }
            
            sip = (ap_order > 0);
            // TODO: mainwindow.Polynomial1.itemindex := ...
            
            if (head_out.ra0 != 0 || head_out.dec0 != 0 || equinox != 2000) {
                if (equinox != 2000) {
                    double jd_obs = (equinox - 2000) * 365.25 + 2451545;
                    {
                        const auto p = ephem::precess_iau1976({head_out.ra0, head_out.dec0}, jd_obs, 2451545);
                        head_out.ra0 = p.ra;
                        head_out.dec0 = p.dec;
                    }
                    if (dec_mount < 999) {
                        const auto p = ephem::precess_iau1976({ra_mount, dec_mount}, jd_obs, 2451545);
                        ra_mount = p.ra;
                        dec_mount = p.dec;
                    }
                }
                // TODO: mainwindow.ra1/dec1 update.
            }
        }
        
        if (head_out.set_temperature == 999)
            head_out.set_temperature = static_cast<int>(std::lround(ccd_temperature));
            
        unsaved_import = false;
        
        if (!load_data) {
            return true;
        }
        
        // ---- read image data --------------------------------------------
        const int bytes_per_pixel = std::abs(nrbits) / 8;
        const int line_bytes = head_out.width * bytes_per_pixel;
        if (line_bytes > static_cast<int>(kFitsBufWide)) {
            // TODO: GUI textout. Abort.
            return false;
        }
        
        img_loaded2.assign(
            head_out.naxis3,
            std::vector<std::vector<float>>(
                head_out.height, std::vector<float>(head_out.width, 0.0f)));
                
        std::vector<std::uint8_t> row_buf(line_bytes);
        
        auto read_row = [&]() -> bool {
            the_file.read(reinterpret_cast<char*>(row_buf.data()), line_bytes);
            return static_cast<bool>(the_file);
        };
        
        if (nrbits == 16) {
            for (int k = 0; k < head_out.naxis3; ++k)
                for (int j = 0; j < head_out.height; ++j) {
                    if (!read_row()) head_out.naxis = 0;
                    for (int x = 0; x < head_out.width; ++x) {
                        std::uint16_t w =
                            static_cast<std::uint16_t>(row_buf[x * 2]) << 8 |
                            static_cast<std::uint16_t>(row_buf[x * 2 + 1]);
                        std::int16_t s = static_cast<std::int16_t>(w);
                        float col = s * bscale + static_cast<float>(bzero);
                        img_loaded2[k][j][x] = col;
                        if (col > measured_max) measured_max = col;
                    }
                }
        } else if (nrbits == -32) {
            for (int k = 0; k < head_out.naxis3; ++k)
                for (int j = 0; j < head_out.height; ++j) {
                    if (!read_row()) head_out.naxis = 0;
                    for (int x = 0; x < head_out.width; ++x) {
                        std::uint32_t lw;
                        std::memcpy(&lw, &row_buf[x * 4], 4);
                        lw = swap_endian(lw);
                        float v = std::bit_cast<float>(lw);
                        float col = v * bscale + static_cast<float>(bzero);
                        if (std::isnan(col)) col = measured_max;
                        img_loaded2[k][j][x] = col;
                        if (col > measured_max) measured_max = col;
                    }
                }
        } else if (nrbits == 8) {
            for (int k = 0; k < head_out.naxis3; ++k)
                for (int j = 0; j < head_out.height; ++j) {
                    if (!read_row()) head_out.naxis = 0;
                    for (int x = 0; x < head_out.width; ++x)
                        img_loaded2[k][j][x] = row_buf[x] * bscale + static_cast<float>(bzero);
                }
        } else if (nrbits == 24) {
            for (int j = 0; j < head_out.height; ++j) {
                if (!read_row()) head_out.naxis = 0;
                for (int x = 0; x < head_out.width; ++x) {
                    img_loaded2[0][j][x] = row_buf[x * 3 + 0];
                    img_loaded2[1][j][x] = row_buf[x * 3 + 1];
                    img_loaded2[2][j][x] = row_buf[x * 3 + 2];
                }
            }
        } else if (nrbits == 32) {
            for (int k = 0; k < head_out.naxis3; ++k)
                for (int j = 0; j < head_out.height; ++j) {
                    if (!read_row()) head_out.naxis = 0;
                    for (int x = 0; x < head_out.width; ++x) {
                        std::uint32_t lw;
                        std::memcpy(&lw, &row_buf[x * 4], 4);
                        lw = swap_endian(lw);
                        std::int32_t i32 = static_cast<std::int32_t>(lw);
                        float col = i32 * bscale + static_cast<float>(bzero);
                        img_loaded2[k][j][x] = col;
                        if (col > measured_max) measured_max = col;
                    }
                }
        } else if (nrbits == -64) {
            for (int k = 0; k < head_out.naxis3; ++k)
                for (int j = 0; j < head_out.height; ++j) {
                    if (!read_row()) head_out.naxis = 0;
                    for (int x = 0; x < head_out.width; ++x) {
                        std::uint64_t qw;
                        std::memcpy(&qw, &row_buf[x * 8], 8);
                        qw = swap_endian(qw);
                        double v = std::bit_cast<double>(qw);
                        float col = static_cast<float>(v * bscale + bzero);
                        img_loaded2[k][j][x] = col;
                        if (col > measured_max) measured_max = col;
                    }
                }
        }
        
        // Rescale floating-point images outside the 0..65535 range.
        if (nrbits <= -32 || nrbits == 32) {
            float scalefactor = 1.0f;
            if (measured_max > 0)
                if (measured_max <= 1.0f * 1.5f || measured_max > 65535.0f * 1.5f)
                    scalefactor = 65535.0f / measured_max;
            if (scalefactor != 1.0f) {
                for (int k = 0; k < head_out.naxis3; ++k)
                    for (int j = 0; j < head_out.height; ++j)
                        for (int x = 0; x < head_out.width; ++x)
                            img_loaded2[k][j][x] *= scalefactor;
                head_out.datamax_org = 65535;
            } else
                head_out.datamax_org = measured_max;
        } else if (nrbits == 8)
            head_out.datamax_org = 255;
        else if (nrbits == 24) {
            head_out.datamax_org = 255;
            nrbits = 8;
        } else
            head_out.datamax_org = measured_max;
            
        bck.backgr = head_out.datamin_org;
        cwhite = head_out.datamax_org;
        
        result = (head_out.naxis != 0);
        reader_position += head_out.width * head_out.height * (std::abs(nrbits) / 8);
    } else if (head_out.naxis == 2 && (bintable || asciitable)) {
        // TODO: read table -- routes to mainwindow.Memo3 in the original.
        // Skipped: table content is GUI-side. Still advance reader_position.
        if (bintable) extend_type = 3;
        if (asciitable) extend_type = 2;
        reader_position += head_out.width * head_out.height;
    }
    
    if (!last_extension) {
        if (file_size - reader_position > static_cast<std::int64_t>(kFitsBlockSize)) {
            // TODO: mainwindow.Memo3 / pagecontrol1 -- GUI hooks.
            last_extension = false;
        } else {
            last_extension = true;
        }
    }
    // TODO: mainwindow.tabsheet1.caption := ...
    
    return result;
}

// ---------------------------------------------------------------------------
// read_keys_memo
// ---------------------------------------------------------------------------

namespace {
// Read a 20-char float starting at memo[index] position 11 (0-based 10).
double read_float_card(const std::string& card) {
    if (card.size() < 11) return 0;
    std::string sub = card.substr(10, std::min<std::size_t>(20, card.size() - 10));
    try { return std::stod(sub); } catch (...) { return 0; }
}
int read_integer_card(const std::string& card) {
    return static_cast<int>(std::lround(read_float_card(card)));
}
std::string read_string_card(const std::string& card) {
    if (card.size() < 11) return std::string{};
    std::string sub = card.substr(10, std::min<std::size_t>(70, card.size() - 10));
    auto p1 = sub.find('\'');
    if (p1 == std::string::npos) {
        // Allow reading bare floats/integers as strings.
        return trim(sub.substr(0, std::min<std::size_t>(11, sub.size())));
    }
    auto p2 = sub.find('\'', p1 + 1);
    if (p2 == std::string::npos) p2 = 21;
    return trim(sub.substr(p1 + 1, p2 - p1 - 1));
}
}  // namespace

void read_keys_memo(bool light, astap::Header& head_out,
                    std::vector<std::string>& memo) {
    int count1 = static_cast<int>(memo.size()) - 1 - 1;
    double ccd_temperature = 999;
    annotated = false;
    
    auto ensure_default = [&](std::size_t pos, std::string_view key,
                              std::string_view dflt) {
        if (pos > memo.size()) return;
        if (pos == 0) return;
        if (memo[pos - 1].substr(0, key.size()) != key) {
            memo.insert(memo.begin() + pos, std::string(dflt));
            ++count1;
        }
    };
    
    int index = 1;
    while (index <= count1) {
        const std::string& line = memo[index];
        std::string key = line.substr(0, std::min<std::size_t>(9, line.size()));
        
        if (index == 1) ensure_default(1, "BITPIX  =",
            "BITPIX  =                   16 / Bits per entry                                 ");
        if (index == 2) ensure_default(2, "NAXIS   =",
            "NAXIS   =                    2 / Number of dimensions                           ");
        if (index == 3) ensure_default(3, "NAXIS1  =",
            "NAXIS1  =                  100 / length of x axis                               ");
        if (index == 4) ensure_default(4, "NAXIS2  =",
            "NAXIS2  =                  100 / length of y axis                               ");
        if (index == 5 && head_out.naxis3 > 1)
            ensure_default(5, "NAXIS3  =",
                "NAXIS3  =                    3 / length of z axis (mostly colors)               ");
                
        if      (key == "CD1_1   =") head_out.cd1_1 = read_float_card(line);
        else if (key == "CD1_2   =") head_out.cd1_2 = read_float_card(line);
        else if (key == "CD2_1   =") head_out.cd2_1 = read_float_card(line);
        else if (key == "CD2_2   =") head_out.cd2_2 = read_float_card(line);
        else if (key == "CRPIX1  =") head_out.crpix1 = read_float_card(line);
        else if (key == "CRPIX2  =") head_out.crpix2 = read_float_card(line);
        else if (key == "CRVAL1  =") head_out.ra0 = read_float_card(line) * kPi / 180;
        else if (key == "CRVAL2  =") head_out.dec0 = read_float_card(line) * kPi / 180;
        else if (key == "RA      =") {
            ra_mount = read_float_card(line) * kPi / 180;
            if (head_out.ra0 == 0) head_out.ra0 = ra_mount;
        } else if (key == "DEC     =") {
            dec_mount = read_float_card(line) * kPi / 180;
            if (head_out.dec0 == 0) head_out.dec0 = dec_mount;
        } else if (key == "OBJCTRA =" && ra_mount >= 999) {
            // TODO: mainwindow.ra1.text := read_string_card(line);
            ra_mount = ra_radians;
        } else if (key == "OBJCTDEC=" && dec_mount >= 999) {
            // TODO: mainwindow.dec1.text := read_string_card(line);
            dec_mount = dec_radians;
        } else if (key == "OBJECT  =") object_name = read_string_card(line);
        else if (key == "EXPOSURE=" || key == "EXPTIME =") head_out.exposure = read_float_card(line);
        else if (key == "XBINNING=") head_out.xbinning = read_integer_card(line);
        else if (key == "YBINNING=") head_out.ybinning = read_integer_card(line);
        else if (key == "FOCALLEN=") focallen = read_float_card(line);
        else if (key == "XPIXSZ  =") head_out.xpixsz = read_float_card(line);
        else if (key == "YPIXSZ  =") head_out.ypixsz = read_float_card(line);
        else if (key == "CDELT1  =") head_out.cdelt1 = read_float_card(line);
        else if (key == "CDELT2  =") head_out.cdelt2 = read_float_card(line);
        else if (key == "EQUINOX =") equinox = read_float_card(line);
        else if (key == "SECPIX2 =" || key == "PIXSCALE=" || key == "SCALE   =") {
            if (head_out.cdelt2 == 0)
                head_out.cdelt2 = read_float_card(line) / 3600;
        }
        else if (key == "GAIN    =") head_out.gain = read_string_card(line).substr(0, 5);
        else if (key == "EGAIN   =") head_out.egain = read_string_card(line).substr(0, 5);
        else if (key == "CCD-TEMP=") ccd_temperature = std::lround(read_float_card(line));
        else if (key == "SET-TEMP=") head_out.set_temperature = static_cast<int>(std::lround(read_float_card(line)));
        else if (key == "LIGH_CNT=") head_out.light_count = read_integer_card(line);
        else if (key == "DARK_CNT=") head_out.dark_count = read_integer_card(line);
        else if (key == "FLAT_CNT=") head_out.flat_count = read_integer_card(line);
        else if (key == "BIAS_CNT=") head_out.flatdark_count = read_integer_card(line);
        else if (key == "PEDESTAL=") head_out.pedestal = std::lround(read_float_card(line));
        else if (key == "CALSTAT =") head_out.calstat = read_string_card(line);
        else if (key == "FILTER  =") head_out.filter_name = read_string_card(line);
        else if (key == "ISSUES  =") head_out.issues = read_string_card(line);
        else if (key == "DATE-OBS=") head_out.date_obs = read_string_card(line);
        else if (key == "BAYERPAT=") bayerpat = read_string_card(line);
        else if (key == "ROWORDER=") roworder = read_string_card(line);
        else if (key == "XBAYROFF=") xbayroff = static_cast<int>(std::lround(read_float_card(line)));
        else if (key == "YBAYROFF=") ybayroff = static_cast<int>(std::lround(read_float_card(line)));
        else if (key == "PRESSURE=") pressure = std::lround(read_float_card(line));
        else if (key == "AIRMASS =") airmass = std::lround(read_float_card(line));
        else if (key == "AOCBAROM=") pressure = std::lround(read_float_card(line));
        else if (key == "FOCUSTEM=" || key == "FOCTEMP =" ||
                 key == "AMB-TEMP=" || key == "AOCAMBT =")
            focus_temp = std::lround(read_float_card(line));
        else if (key == "MZEROR  =") head_out.mzero = read_float_card(line);
        else if (key == "MZEROAPT=") head_out.mzero_radius = read_float_card(line);
        else if (key == "MZEROPAS=") head_out.passband_database = read_string_card(line);
        else if (key == "ANNOTATE=") annotated = true;
        else if (key == "TELESCOP=") telescop = read_string_card(line);
        else if (key == "INSTRUME=") instrum = read_string_card(line);
        else if (key == "CENTALT =") centalt = read_string_card(line);
        else if (key == "SITELAT =") sitelat = read_string_card(line);
        else if (key == "SITELONG=") sitelong = read_string_card(line);
        
        if (key == sqm_key + "=") sqm_value = read_string_card(line);
        
        ++index;
    }
    
    if (light && (head_out.ra0 != 0 || head_out.dec0 != 0)) {
        if (equinox != 2000) {
            double jd_obs = (equinox - 2000) * 365.25 + 2451545;
            {
                const auto p = ephem::precess_iau1976({head_out.ra0, head_out.dec0}, jd_obs, 2451545);
                head_out.ra0 = p.ra;
                head_out.dec0 = p.dec;
            }
            if (dec_mount < 999) {
                const auto p = ephem::precess_iau1976({ra_mount, dec_mount}, jd_obs, 2451545);
                ra_mount = p.ra;
                dec_mount = p.dec;
            }
        }
        // TODO: mainwindow.ra1.text / dec1.text update.
    }
    
    if (head_out.cd1_1 != 0 &&
        (head_out.cdelt1 == 0 || head_out.crota2 >= 999))
        new_to_old_WCS(head_out);
    else if (head_out.cd1_1 == 0 && head_out.crota2 < 999 && head_out.cdelt2 != 0)
        old_to_new_WCS(head_out);
        
    if (head_out.cd1_1 == 0 && head_out.cdelt2 == 0)
        if (focallen != 0 && head_out.xpixsz != 0)
            head_out.cdelt2 = 180.0 / (kPi * 1000) * head_out.xpixsz / focallen;
            
    if (head_out.crota2 > 999) head_out.crota2 = 0;
    if (head_out.crota1 > 999) head_out.crota1 = head_out.crota2;
    
    if (head_out.set_temperature == 999)
        head_out.set_temperature = static_cast<int>(std::lround(ccd_temperature));
}

// ---------------------------------------------------------------------------
// Header-edit primitives
// ---------------------------------------------------------------------------

namespace {

// Helper used by update_float / update_integer: replace cards in-place
// using the legacy "memo as one big text buffer" model. We materialise the
// text, edit, then split back into cards.
std::string memo_to_text(const std::vector<std::string>& memo) {
    std::string out;
    out.reserve(memo.size() * (kFitsCardSize + 1));
    for (auto const& l : memo) {
        out.append(l);
        out.push_back('\n');
    }
    return out;
}
void text_to_memo(const std::string& buf, std::vector<std::string>& memo) {
    memo.clear();
    std::size_t start = 0;
    for (std::size_t i = 0; i < buf.size(); ++i) {
        if (buf[i] == '\n') {
            memo.emplace_back(buf, start, i - start);
            start = i + 1;
        }
    }
    if (start < buf.size()) memo.emplace_back(buf.substr(start));
}
    
}  // namespace

void update_float(std::vector<std::string>& memo,
                  std::string_view inpt, std::string_view comment1,
                  bool preserve_comment, double x) {
    std::string s = floattostr20(x);
    std::string text = memo_to_text(memo);
    std::size_t cnt = pos1(inpt, text);
    if (cnt > 0) {
        std::size_t line_end = line_end_after(text, cnt + 1);
        std::string aline = text.substr(cnt - 1, line_end - cnt);
        
        if (preserve_comment && aline.size() >= 32 && aline[31] == '/') {
            // delete pos 11..30 (1-based) -> 0-based [10..30)
            if (aline.size() > 10) aline.erase(10, std::min<std::size_t>(20, aline.size() - 10));
            aline.insert(10, s);
        } else {
            if (aline.size() > 10) aline.erase(10);
            aline += s;
            aline += comment1;
        }
        
        // Splice back with space-padding so the line keeps its width.
        for (std::size_t k = 0; k < line_end - cnt; ++k) {
            if (k < aline.size())
                text[cnt - 1 + k] = aline[k];
            else
                text[cnt - 1 + k] = ' ';
        }
        text_to_memo(text, memo);
        return;
    }
    // not found, append before the END line
    if (memo.empty()) memo.emplace_back("END");
    memo.insert(memo.end() - 1, std::string(inpt) + " " + s + std::string(comment1));
}

void update_integer(std::vector<std::string>& memo,
                    std::string_view inpt, std::string_view comment1, int x) {
    std::string s = inttostr20(x);
    std::string text = memo_to_text(memo);
    std::size_t cnt = pos1(inpt, text);
    if (cnt > 0) {
        std::size_t line_end = line_end_after(text, cnt + 1);
        std::string aline = text.substr(cnt - 1, line_end - cnt);
        if (aline.size() > 10) aline.erase(10);
        aline += s;
        aline += comment1;
        for (std::size_t k = 0; k < line_end - cnt; ++k) {
            if (k < aline.size()) text[cnt - 1 + k] = aline[k];
            else text[cnt - 1 + k] = ' ';
        }
        text_to_memo(text, memo);
        return;
    }
    if (memo.empty()) memo.emplace_back("END");
    auto built = std::string(inpt) + " " + s + std::string(comment1);
    if (inpt == "NAXIS1  =" && memo.size() > 3)
        memo.insert(memo.begin() + 3, built);
    else if (inpt == "NAXIS2  =" && memo.size() > 4)
        memo.insert(memo.begin() + 4, built);
    else if (inpt == "NAXIS3  =" && memo.size() > 5)
        memo.insert(memo.begin() + 5, built);
    else
        memo.insert(memo.end() - 1, built);
}

void add_integer(std::vector<std::string>& memo,
                 std::string_view inpt, std::string_view comment1, int x) {
    std::string s = inttostr20(x);
    if (memo.empty()) memo.emplace_back("END");
    memo.insert(memo.end() - 1, std::string(inpt) + " " + s + std::string(comment1));
}

void update_generic(std::vector<std::string>& memo,
                    std::string_view message_key,
                    std::string_view message_value,
                    std::string_view message_comment) {
    std::string mk(message_key);
    std::string mv(message_value);
    
    if (mk.find("HISTORY") == std::string::npos &&
        mk.find("COMMENT") == std::string::npos) {
        while (mv.size() < 20) mv.insert(mv.begin(), ' ');
        while (mk.size() < 8) mk.push_back(' ');
        
        for (int c = static_cast<int>(memo.size()) - 1; c >= 0; --c) {
            if (memo[c].find(mk) != std::string::npos) {
                memo[c] = mk + "= " + mv + " / " + std::string(message_comment);
                return;
            }
        }
        if (memo.empty()) memo.emplace_back("END");
        memo.insert(memo.end() - 1,
                    mk + "= " + mv + " / " + std::string(message_comment));
    } else {
        if (memo.empty()) memo.emplace_back("END");
        memo.insert(memo.end() - 1,
                    mk + " " + mv + std::string(message_comment));
    }
}

void update_text(std::vector<std::string>& memo,
                 std::string_view inpt, std::string_view comment1) {
    for (int c = static_cast<int>(memo.size()) - 1; c >= 0; --c) {
        if (memo[c].find(inpt) != std::string::npos) {
            memo[c] = std::string(inpt) + " " + std::string(comment1);
            return;
        }
    }
    if (memo.empty()) memo.emplace_back("END");
    memo.insert(memo.end() - 1, std::string(inpt) + " " + std::string(comment1));
}

void update_longstr(std::vector<std::string>& memo,
                    std::string_view inpt, std::string_view thestr) {
    // Delete existing keyword and any continuation lines.
    for (int c = static_cast<int>(memo.size()) - 1; c >= 0; --c) {
        if (memo[c].find(inpt) != std::string::npos) {
            std::size_t idx = static_cast<std::size_t>(c);
            memo.erase(memo.begin() + idx);
            while (idx < memo.size() &&
                   memo[idx].find("CONTINUE=") != std::string::npos)
                memo.erase(memo.begin() + idx);
        }
    }
    
    if (memo.empty()) memo.emplace_back("END");
    
    const std::size_t m = thestr.size();
    if (m > 68) {
        memo.insert(memo.end() - 1,
                    std::string(inpt) + " '" + std::string(thestr.substr(0, 67)) + "&'");
        std::size_t k = 67;  // 1-based 68 -> 0-based 67
        for (;;) {
            std::string amp = (m - k > 67) ? "&" : "";
            std::size_t take = std::min<std::size_t>(67, m - k);
            memo.insert(memo.end() - 1,
                        "CONTINUE= '" + std::string(thestr.substr(k, take)) + amp + "'");
            k += 67;
            if (k >= m) break;
        }
    } else {
        memo.insert(memo.end() - 1,
                    std::string(inpt) + " '" + std::string(thestr) + "'");
    }
}

void add_text(std::vector<std::string>& memo,
              std::string_view inpt, std::string_view comment1) {
    if (memo.empty()) memo.emplace_back("END");
    std::string c(comment1);
    if (c.size() > 79 - inpt.size()) c.resize(79 - inpt.size());
    memo.insert(memo.end() - 1, std::string(inpt) + " " + c);
}

void add_long_comment(std::vector<std::string>& memo,
                      std::string_view descrip) {
    if (memo.empty()) memo.emplace_back("END");
    std::size_t i = 0;
    const std::size_t j = descrip.size();
    while (i < j) {
        std::size_t take = std::min<std::size_t>(72, j - i);
        memo.insert(memo.end() - 1,
                    "COMMENT " + std::string(descrip.substr(i, take)));
        i += 72;
    }
}

void remove_key(std::vector<std::string>& memo, std::string_view inpt, bool all) {
    for (int c = static_cast<int>(memo.size()) - 1; c >= 0; --c) {
        if (memo[c].find(inpt) != std::string::npos) {
            memo.erase(memo.begin() + c);
            if (!all) return;
        }
    }
}

void remove_solution(std::vector<std::string>& memo, bool keep_wcs,
                     astap::Header& head_out) {
    auto rm = [&](std::string_view k) { remove_key(memo, k, false); };
    
    if (!keep_wcs) {
        head_out.cd1_1 = 0;
        rm("CD1_1   ="); rm("CD1_2   ="); rm("CD2_1   ="); rm("CD2_2   =");
    }
    
    a_order = 0;
    rm("A_ORDER ="); rm("A_0_0   ="); rm("A_0_1   ="); rm("A_0_2   =");
    rm("A_0_3   ="); rm("A_1_0   ="); rm("A_1_1   ="); rm("A_1_2   =");
    rm("A_2_0   ="); rm("A_2_1   ="); rm("A_3_0   =");
    
    rm("B_ORDER ="); rm("B_0_0   ="); rm("B_0_1   ="); rm("B_0_2   =");
    rm("B_0_3   ="); rm("B_1_0   ="); rm("B_1_1   ="); rm("B_1_2   =");
    rm("B_2_0   ="); rm("B_2_1   ="); rm("B_3_0   =");
    
    rm("AP_ORDER="); rm("AP_0_0  ="); rm("AP_0_1  ="); rm("AP_0_2  =");
    rm("AP_0_3  ="); rm("AP_1_0  ="); rm("AP_1_1  ="); rm("AP_1_2  =");
    rm("AP_2_0  ="); rm("AP_2_1  ="); rm("AP_3_0  =");
    
    rm("BP_ORDER="); rm("BP_0_0  ="); rm("BP_0_1  ="); rm("BP_0_2  =");
    rm("BP_0_3  ="); rm("BP_1_0  ="); rm("BP_1_1  ="); rm("BP_1_2  =");
    rm("BP_2_0  ="); rm("BP_2_1  ="); rm("BP_3_0  =");
    
    rm("CROTA1  ="); rm("CROTA2  =");
}

// ---------------------------------------------------------------------------
// WCS conversion
// ---------------------------------------------------------------------------

void old_to_new_WCS(astap::Header& head_out) noexcept {
    head_out.cd1_1 = +head_out.cdelt1 * std::cos(head_out.crota1 * kPi / 180);
    head_out.cd2_1 = +head_out.cdelt1 * std::sin(head_out.crota1 * kPi / 180);
    head_out.cd1_2 = -head_out.cdelt2 * std::sin(head_out.crota2 * kPi / 180);
    head_out.cd2_2 = +head_out.cdelt2 * std::cos(head_out.crota2 * kPi / 180);
}

void new_to_old_WCS(astap::Header& head_out) noexcept {
    double crota_1, crota_2;
    
    if (head_out.cd2_1 > 0)
        crota_1 = std::atan2(-head_out.cd2_1, -head_out.cd1_1);
    else if (head_out.cd2_1 < 0)
        crota_1 = std::atan2(+head_out.cd2_1, +head_out.cd1_1);
    else
        crota_1 = 0;
        
    if (head_out.cd1_2 > 0)
        crota_2 = std::atan2(-head_out.cd1_2, head_out.cd2_2);
    else if (head_out.cd1_2 < 0)
        crota_2 = std::atan2(head_out.cd1_2, -head_out.cd2_2);
    else
        crota_2 = 0;
        
    if (std::abs(head_out.cd1_1) > std::abs(head_out.cd2_1)) {
        head_out.cdelt1 = +head_out.cd1_1 / std::cos(crota_1);
        head_out.cdelt2 = +head_out.cd2_2 / std::cos(crota_2);
    } else {
        head_out.cdelt1 = +head_out.cd2_1 / std::sin(crota_1);
        head_out.cdelt2 = -head_out.cd1_2 / std::sin(crota_2);
    }
    
    head_out.crota1 = crota_1 * 180 / kPi;
    head_out.crota2 = crota_2 * 180 / kPi;
    
    if (head_out.cdelt2 < 0) {
        if (head_out.crota2 < 0) {
            head_out.crota2 += 180;
            head_out.crota1 += 180;
        } else {
            head_out.crota2 -= 180;
            head_out.crota1 -= 180;
        }
        head_out.cdelt2 = -head_out.cdelt2;
        head_out.cdelt1 = -head_out.cdelt1;
    }
}

// ---------------------------------------------------------------------------
// save_fits
// ---------------------------------------------------------------------------

namespace {
// Float -> big-endian uint32.
std::uint32_t int_ieee4_reverse(float x) {
    return swap_endian(std::bit_cast<std::uint32_t>(x));
}
}  // namespace

bool save_fits(const astap::ImageArray& img,
               std::vector<std::string>& memo,
               const std::filesystem::path& filen2,
               int type1,
               [[maybe_unused]] bool override2) {
    if (img.empty()) {
        memo2_message("Error,  no image");
        return false;
    }
    auto colours5   = static_cast<int>(img.size());
    auto height5    = static_cast<int>(img[0].size());
    auto width5     = static_cast<int>(img[0].empty() ? 0 : img[0][0].size());
    auto dimensions = (colours5 == 1) ? 2 : 3;
    
    if (type1 == 24 && colours5 < 3) {
        // TODO: messagebox
        return false;
    }
    
    // TODO: override2 / overwrite confirmation -- GUI only.
    // TODO: extend_type==1 multi-extension handling -- GUI only.
    
    filename2 = filen2.string();
    progress_indicator(0, "");
    
    std::ofstream fs(filen2, std::ios::binary | std::ios::trunc);
    if (!fs) return false;
    
    int bzero2 = 0;
    
    if (type1 != 24) {
        update_integer(memo, "BITPIX  =",
            " / Bits per entry                                 ", type1);
        update_integer(memo, "NAXIS   =",
            " / Number of dimensions                           ", dimensions);
        update_integer(memo, "NAXIS1  =",
            " / length of x axis                               ", width5);
        update_integer(memo, "NAXIS2  =",
            " / length of y axis                               ", height5);
        if (colours5 != 1)
            update_integer(memo, "NAXIS3  =",
                " / length of z axis (mostly colors)               ", colours5);
        else
            remove_key(memo, "NAXIS3  ", false);
            
        bzero2 = (type1 == 16) ? 32768 : 0;
        
        update_integer(memo, "BZERO   =",
            " / physical_value = BZERO + BSCALE * array_value  ", bzero2);
        update_integer(memo, "BSCALE  =",
            " / physical_value = BZERO + BSCALE * array_value  ", 1);
            
        if (type1 != 8) {
            update_integer(memo, "DATAMIN =",
                " / Minimum data value                             ",
                static_cast<int>(std::lround(head.datamin_org)));
            update_integer(memo, "DATAMAX =",
                " / Maximum data value                             ",
                static_cast<int>(std::lround(head.datamax_org)));
            update_integer(memo, "CBLACK  =",
                " / Black point used for displaying image.         ",
                static_cast<int>(std::lround(bck.backgr)));
            update_integer(memo, "CWHITE  =",
                " / White point used for displaying the image.     ",
                static_cast<int>(std::lround(cwhite)));
        } else {
            update_integer(memo, "DATAMIN =",
                " / Minimum data value                             ", 0);
            update_integer(memo, "DATAMAX =",
                " / Maximum data value                             ", 255);
        }
    } else {
        // Special 8-bit RGB packed via NAXIS1=3.
        update_integer(memo, "BITPIX  =",
            " / Bits per entry                                 ", 8);
        update_integer(memo, "NAXIS   =",
            " / Number of dimensions                           ", dimensions);
        update_integer(memo, "NAXIS1  =",
            " / length of x axis                               ", 3);
        update_integer(memo, "NAXIS2  =",
            " / length of y axis                               ", width5);
        update_integer(memo, "NAXIS3  =",
            " / length of z axis (mostly colors)               ", height5);
        update_integer(memo, "DATAMIN =",
            " / Minimum data value                             ", 0);
        update_integer(memo, "DATAMAX =",
            " / Maximum data value                             ", 255);
        update_integer(memo, "BZERO   =",
            " / physical_value = BZERO + BSCALE * array_value  ", 0);
        update_integer(memo, "BSCALE  =",
            " / physical_value = BZERO + BSCALE * array_value  ", 1);
    }
    
    // Write header in 80-byte cards, padded out to a 2880-byte block.
    std::string empty_line(kFitsCardSize, ' ');
    std::size_t i = 0;
    do {
        if (i < memo.size()) {
            std::string padded = pad80(memo[i]);
            fs.write(padded.data(), kFitsCardSize);
        } else {
            fs.write(empty_line.data(), kFitsCardSize);
        }
        ++i;
    } while (!(i >= memo.size() && (i * kFitsCardSize) % kFitsBlockSize == 0));
    
    // Pixel data.
    int progressC = 0;
    auto bump_progress = [&]() {
        ++progressC;
        int pct = static_cast<int>(std::lround(progressC * 100.0 / (colours5 * height5)));
        if (pct % 5 == 0) progress_indicator(pct, "");
    };
    
    if (type1 == 8) {
        // TODO: minimum/maximum from mainwindow.minimum1/maximum1 sliders.
        const int minimum = 0;
        const int maximum = 255;
        std::vector<std::uint8_t> row(width5);
        for (int k = 0; k < colours5; ++k)
            for (int j = 0; j < height5; ++j) {
                bump_progress();
                for (int x = 0; x < width5; ++x) {
                    float dd = img[k][j][x];
                    int dum = static_cast<int>(std::lround(
                        (dd - minimum) * 255.0 / (maximum - minimum)));
                    if (dum < 0) dum = 0;
                    if (dum > 255) dum = 255;
                    row[x] = static_cast<std::uint8_t>(dum);
                }
                fs.write(reinterpret_cast<const char*>(row.data()), width5);
            }
    } else if (type1 == 24) {
        const int minimum = 0;
        const int maximum = 255;
        std::vector<std::uint8_t> row(width5 * 3);
        for (int j = 0; j < height5; ++j) {
            bump_progress();
            for (int x = 0; x < width5; ++x)
                for (int k = 0; k < 3; ++k) {
                    float dd = img[k][j][x];
                    int dum = static_cast<int>(std::lround(
                        (dd - minimum) * 255.0 / (maximum - minimum)));
                    if (dum < 0) dum = 0;
                    if (dum > 255) dum = 255;
                    row[x * 3 + k] = static_cast<std::uint8_t>(dum);
                }
            fs.write(reinterpret_cast<const char*>(row.data()), width5 * 3);
        }
    } else if (type1 == 16) {
        std::vector<std::uint16_t> row(width5);
        for (int k = 0; k < colours5; ++k)
            for (int j = 0; j < height5; ++j) {
                bump_progress();
                for (int x = 0; x < width5; ++x) {
                    int v = std::max(0, std::min(65535,
                                static_cast<int>(std::lround(img[k][j][x]))));
                    int dum = (v - bzero2) & 0xFFFF;
                    row[x] = swap16(static_cast<std::uint16_t>(dum));
                }
                fs.write(reinterpret_cast<const char*>(row.data()), width5 * 2);
            }
    } else if (type1 == -32) {
        std::vector<std::uint32_t> row(width5);
        for (int k = 0; k < colours5; ++k)
            for (int j = 0; j < height5; ++j) {
                bump_progress();
                for (int x = 0; x < width5; ++x)
                    row[x] = int_ieee4_reverse(img[k][j][x]);
                fs.write(reinterpret_cast<const char*>(row.data()), width5 * 4);
            }
    }
    
    // Pad data section to the next 2880-byte block.
    std::int64_t pos = fs.tellp();
    std::int64_t remain = (kFitsBlockSize - (pos % kFitsBlockSize)) % kFitsBlockSize;
    if (remain > 0) {
        std::vector<char> zeros(remain, 0);
        fs.write(zeros.data(), remain);
    }
    
    unsaved_import = false;
    progress_indicator(-100, "");
    return static_cast<bool>(fs);
}

// ---------------------------------------------------------------------------
// savefits_update_header
// ---------------------------------------------------------------------------

bool savefits_update_header(std::vector<std::string>& memo,
                            const std::filesystem::path& filen2) {
    auto filename_tmp = filen2;
    filename_tmp.replace_extension(".tmp");
    
    std::ifstream src(filen2, std::ios::binary);
    std::ofstream dst(filename_tmp, std::ios::binary | std::ios::trunc);
    if (!src || !dst) return false;
    
    src.seekg(0, std::ios::end);
    const std::int64_t file_size = src.tellg();
    src.seekg(0, std::ios::beg);
    
    // Find END card in source so we can skip its header.
    std::int64_t reader_position = 0;
    bool endfound = false;
    std::array<char, 2881> header{};
    int i = 0;
    while (!endfound && i < static_cast<int>(sizeof(header)) - 16) {
        src.read(&header[i], 80);
        if (!src) return false;
        reader_position += 80;
        endfound = (header[i] == 'E' && header[i + 1] == 'N' &&
                    header[i + 2] == 'D' && header[i + 3] == ' ');
        if (!endfound) i += 80;
    }
    if (!endfound) {
        memo2_message("Abort, error reading source FITS file!!");
        return false;
    }
    
    // Skip the rest of the current 2880 block.
    double fract = (reader_position % kFitsBlockSize) / static_cast<double>(kFitsBlockSize);
    if (fract != 0) {
        int skip = static_cast<int>(std::lround((1.0 - fract) * kFitsBlockSize));
        std::vector<char> dump(skip);
        src.read(dump.data(), skip);
        reader_position += skip;
    }
    
    // Write the new header from the memo.
    std::string empty_line(kFitsCardSize, ' ');
    std::size_t k = 0;
    do {
        if (k < memo.size())
            dst.write(pad80(memo[k]).data(), kFitsCardSize);
        else
            dst.write(empty_line.data(), kFitsCardSize);
        ++k;
    } while (!(k >= memo.size() && (k * kFitsCardSize) % kFitsBlockSize == 0));
    
    // Stream the rest of the source data block by block.
    std::vector<char> buf(kFitsBufWide);
    while (reader_position < file_size) {
        std::int64_t to_read = std::min<std::int64_t>(buf.size(), file_size - reader_position);
        src.read(buf.data(), to_read);
        reader_position += to_read;
        dst.write(buf.data(), to_read);
    }
    
    src.close();
    dst.close();
    
    std::error_code ec;
    std::filesystem::remove(filen2, ec);
    if (ec) return false;
    std::filesystem::rename(filename_tmp, filen2, ec);
    return !ec;
}

// ---------------------------------------------------------------------------
// unpack_cfitsio / pack_cfitsio
// ---------------------------------------------------------------------------

bool unpack_cfitsio(std::filesystem::path& filename3) {
    const std::string commando = "-D";
    
    if (filename3.string().find('(') != std::string::npos) {
        // funpack chokes on '(' in filenames.
        std::string fname = filename3.filename().string();
        std::replace(fname.begin(), fname.end(), '(', '_');
        auto newpath = filename3.parent_path() / fname;
        std::error_code ec;
        std::filesystem::rename(filename3, newpath, ec);
        if (!ec) filename3 = newpath;
        if (filename3.string().find('(') != std::string::npos) {
            memo2_message("Error!. Can not process a path with the \"(\" character");
            return false;
        }
    }
    
    // TODO: living in core/platform.cpp -- use the right funpack path per OS.
    std::string cmd = "funpack " + commando + " \"" + filename3.string() + "\"";
    if (!execute_and_wait(cmd, false)) return false;
    
    auto s = filename3.string();
    auto p = s.find(".fz");
    if (p != std::string::npos) s.erase(p, 3);
    filename3 = s;
    return true;
}

bool pack_cfitsio(const std::filesystem::path& filename3) {
    // TODO: living in core/platform.cpp -- use the right fpack path per OS.
    std::string cmd = "fpack \"" + filename3.string() + "\"";
    return execute_and_wait(cmd, false);
}

// ---------------------------------------------------------------------------
// rice_encoding (experimental)
// ---------------------------------------------------------------------------

bool rice_encoding(const std::vector<std::uint16_t>& inp,
                   int k, int bitdepth,
                   std::vector<std::uint16_t>& outp,
                   int& compressed_size) {
    const std::size_t len = inp.size();
    outp.assign(len, 0);
    
    std::uint32_t m = 1u << k;
    --m;  // bitmask for remainder extraction
    
    compressed_size = 0;
    std::uint32_t bitpointer = 0;
    
    auto put_bit = [&](std::uint8_t value) {
        std::uint32_t y = bitpointer / static_cast<std::uint32_t>(bitdepth);
        std::uint32_t rest = bitpointer - y * bitdepth;
        if (static_cast<int>(y) > compressed_size) {
            ++compressed_size;
            outp[y] = 0;
        }
        outp[y] = static_cast<std::uint16_t>(outp[y] | (value << rest));
        ++bitpointer;
    };
    
    bool result = true;
    for (std::size_t j = 0; j < len; ++j) {
        std::uint32_t q = inp[j] >> k;
        std::uint32_t r = inp[j] & m;
        for (std::uint32_t ii = 0; ii < q; ++ii) put_bit(1);
        put_bit(0);
        for (int ii = k - 1; ii >= 0; --ii)
            put_bit(static_cast<std::uint8_t>((r >> ii) & 1));
        if (compressed_size >= static_cast<int>(len) - 1) {
            result = false;
            break;
        }
    }
    ++compressed_size;
    return result;
}

// ---------------------------------------------------------------------------
// binX2X3_file
// ---------------------------------------------------------------------------

bool binX2X3_file(int binfactor) {
    // TODO: load_fits + bin_X2X3X4 + save_fits / save_tiff16. The bin_X2X3X4
    // implementation lives in image/binning.cpp; once it lands replace this
    // stub with the real wiring.
    astap::Header headX{};
    astap::ImageArray img_temp;
    std::vector<std::string> memox;
    
    if (!load_fits(filename2, /*light=*/true, /*load_data=*/true,
                   /*update_memo=*/true, /*get_ext=*/0,
                   memox, headX, img_temp))
        return false;
        
    // TODO: bin_X2X3X4(img_temp, headX, memox, binfactor);
    remove_key(memox, "BAYERPAT=", false);
    
    // TODO: fits_file_name(filename2) check
    std::filesystem::path src(filename2);
    std::string suffix =
        (binfactor == 2) ? "_bin2x2.fit" : "_bin3x3.fit";
    auto out = src;
    out.replace_extension(suffix);
    filename2 = out.string();
    
    return save_fits(img_temp, memox, out, nrbits, true);
}

// ---------------------------------------------------------------------------
// FITS_BMP
// ---------------------------------------------------------------------------

void FITS_BMP(const std::filesystem::path& filen) {
    [[maybe_unused]] astap::ImageArray img;
    auto memo = std::vector<std::string>{};
    if (!load_fits(filen, /*light=*/true, /*load_data=*/true,
                   /*update_memo=*/true, /*get_ext=*/0, memo, head, img)) {
        return;
    }
    
    // TODO: use_histogram(img, true) -- lives in image/histogram.cpp
    // TODO: bitmap rendering lives in the GUI layer. Skipped here.
}
    
} // namespace
