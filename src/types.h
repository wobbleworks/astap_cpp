#pragma once

// Shared types for the ASTAP++ port. Consolidates aliases and structs that
// were redeclared in per-module headers during the unit-by-unit translations
// from the Object Pascal source. Module-specific types (e.g. IQuadCanvas,
// AviWriter, Y4mWriter) stay in their own headers.

#include <array>
#include <cstdint>
#include <expected>
#include <string>
#include <vector>

namespace astap {

// 3D image buffer indexed [channel][row][col], where channel is in [0, naxis3).
// Pascal: image_array = array of array of array of Single.
using ImageArray = std::vector<std::vector<std::vector<float>>>;

// Variable-width 2-D table used for star lists, quad distance tables, etc.
// Pascal: star_list = array of array of double.
using StarList = std::vector<std::vector<double>>;

// Dense row-major matrix. Pascal: TMatrix = array of array of Double.
using Matrix = std::vector<std::vector<double>>;

// Dynamic 1-D vector of doubles. Pascal: TVector = array of Double.
using Vector = std::vector<double>;

// 3-element affine image-alignment solution (a, b, c in x' = a*x + b*y + c).
using SolutionVector = std::array<double, 3>;

// 2-D star coordinate pair. Pascal: s_star = record x, y: Double; end;
struct SStar {
    double x;
    double y;
};

using StarArray = std::vector<SStar>;

// Cubic bivariate polynomial coefficients produced by calc_trans_cubic —
// maps positions in one star array onto positions in another.
struct TransCoeffs {
    double x00, x10, x01, x20, x11, x02, x30, x21, x12, x03;
    double y00, y10, y01, y20, y11, y02, y30, y21, y12, y03;
};

// Image-background / noise statistics. Pascal: Tbackground.
struct Background {
    double backgr{};
    double star_level{};
    double star_level2{};
    double noise_level{};
};

// A file queued for stacking. Pascal: TfileToDo.
struct FileToDo {
    std::string name;
    int listview_index{};
};

// FITS-style header metadata. Mirrors Pascal Theader from astap_main.pas
// (fields 542-584). Kept flat and default-initialised; callers populate
// what they know.
struct Header {
    int width{};
    int height{};
    int naxis{};
    int naxis3{};

    double crpix1{};
    double crpix2{};
    double cdelt1{};
    double cdelt2{};
    double ra0{};
    double dec0{};
    double crota1{};
    double crota2{};
    double cd1_1{};
    double cd1_2{};
    double cd2_1{};
    double cd2_2{};

    double exposure{};
    double datamin_org{};
    double datamax_org{};
    double xbinning{};
    double ybinning{};
    double xpixsz{};
    double ypixsz{};
    double mzero{};
    double mzero_radius{};
    double magn_limit{};
    double pedestal{};

    int set_temperature{};
    int dark_count{};
    int light_count{};
    int flat_count{};
    int flatdark_count{};
    int focus_pos{};              // FOCUSPOS keyword: focuser stepper position.

    std::string egain;
    std::string gain;
    std::string date_obs;
    std::string date_avg;
    std::string calstat;
    std::string filter_name;
    std::string passband_database;
    std::string airmass;
    std::string issues;
};

// Abstract pixel source used by the AVI and YUV4MPEG2 writers. Decouples
// the encoders from whichever bitmap backend the host is using.
// Return value is packed 0x00RRGGBB.
struct PixelSource {
    virtual std::uint32_t pixel(int x, int y) const = 0;
    virtual ~PixelSource() = default;
};

// Abstract HTTP GET client used by the online-catalog modules (Gaia, VSP,
// VSX, Simbad, Vizier). Host code supplies the concrete implementation
// (libcurl, cpp-httplib, platform HTTP, etc.); the library code never
// depends on any specific HTTP stack. On success returns the full response
// body; on failure returns a human-readable error string.
struct IHttpClient {
    virtual std::expected<std::string, std::string> get(const std::string& url) = 0;
    virtual ~IHttpClient() = default;
};

}  // namespace astap
