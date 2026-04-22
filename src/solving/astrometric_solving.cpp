///----------------------------------------
///      @file astrometric_solving.cpp
///   @ingroup ASTAP++
///     @brief Astrometric plate-solving entry point for ASTAP.
///    @author Ported from Han Kleijn's unit_astrometric_solving.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "astrometric_solving.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <format>
#include <numbers>
#include <string>
#include <vector>

#include "../core/fits.h"
#include "../core/globals.h"
#include "../core/photometry.h"
#include "../core/util.h"
#include "../core/wcs.h"
#include "../reference/star_database.h"
#include "../reference/stars_wide_field.h"
#include "../stacking/stack.h"
#include "calc_trans_cubic.h"
#include "star_align.h"

///----------------------------------------
namespace astap::solving {
///----------------------------------------

// Bring the shared globals into unqualified scope inside this namespace so
// call-sites can keep using their bare names.
using astap::mag2;
using astap::commandline_execution;
using astap::solve_show_log;
using astap::fov_specified;
using astap::esc_pressed;
using astap::errorlevel;
using astap::warning_str;
using astap::filename2;
using astap::astap_version;
using astap::ra_radians;
using astap::dec_radians;
using astap::ra_mount;
using astap::dec_mount;
using astap::sip;
using astap::a_order;
using astap::b_order;
using astap::ap_order;
using astap::bp_order;
using astap::a_0_0; using astap::a_0_1; using astap::a_0_2; using astap::a_0_3;
using astap::a_1_0; using astap::a_1_1; using astap::a_1_2;
using astap::a_2_0; using astap::a_2_1;
using astap::a_3_0;
using astap::b_0_0; using astap::b_0_1; using astap::b_0_2; using astap::b_0_3;
using astap::b_1_0; using astap::b_1_1; using astap::b_1_2;
using astap::b_2_0; using astap::b_2_1;
using astap::b_3_0;
using astap::ap_0_0; using astap::ap_0_1; using astap::ap_0_2; using astap::ap_0_3;
using astap::ap_1_0; using astap::ap_1_1; using astap::ap_1_2;
using astap::ap_2_0; using astap::ap_2_1;
using astap::ap_3_0;
using astap::bp_0_0; using astap::bp_0_1; using astap::bp_0_2; using astap::bp_0_3;
using astap::bp_1_0; using astap::bp_1_1; using astap::bp_1_2;
using astap::bp_2_0; using astap::bp_2_1;
using astap::bp_3_0;
using astap::bck;
using astap::reference::name_database;
using astap::reference::wide_database;
using astap::reference::wide_field_stars;
using astap::reference::cos_telescope_dec;
using astap::reference::database_type;
using astap::reference::database_path;

///----------------------------------------
/// MARK: StackMenuConfig + small helpers
///----------------------------------------

// Solver configuration — sourced from the stacking UI in the original,
// populated from engine globals via populate_stackcfg_from_globals() at
// the top of solve_image.
struct StackMenuConfig {
    std::string star_database;
    std::string radius_search_deg;
    std::string search_fov_deg;
    std::string max_stars;
    std::string quad_tolerance;
    std::string min_star_size;
    int  downsample_for_solving_index = 0;
    bool force_oversize = false;
    bool use_triples    = false;
    bool add_sip        = false;
};

///----------------------------------------
/// MARK: Delegates to real implementations
///----------------------------------------

// Each function preserves the local signature expected by solve_image's
// call sites while forwarding to the real implementation in the correct
// namespace (astap::core, astap::reference, astap::stacking).
StackMenuConfig stackcfg;

inline void find_areas(double ra, double dec, double fov,
                       int& a1, int& a2, int& a3, int& a4,
                       double& f1, double& f2, double& f3, double& f4) {
    astap::reference::find_areas(ra, dec, fov, a1, a2, a3, a4, f1, f2, f3, f4);
}

inline bool open_database(double telescope_dec, int area290) {
    return astap::reference::open_database(telescope_dec, area290);
}

inline bool readdatabase290(double ra, double dec, double fov,
                            double& ra2, double& dec2, double& mag2v, double& bp_rp) {
    return astap::reference::readdatabase290(ra, dec, fov, ra2, dec2, mag2v, bp_rp);
}

inline void read_stars_wide_field() {
    (void)astap::reference::read_stars_wide_field();
}

inline void ang_sep(double ra1, double dec1, double ra2, double dec2, double& sep) {
    astap::core::ang_sep(ra1, dec1, ra2, dec2, sep);
}

inline double fnmodulo(double x, double range) {
    return astap::core::fnmodulo(x, range);
}

inline void get_background(int colour, const ImageArray& img,
                           bool calc_hist, bool calc_noise, astap::Background& b) {
    astap::core::get_background(colour, img, calc_hist, calc_noise, b);
}

inline void check_pattern_filter(ImageArray& img) {
    astap::stacking::check_pattern_filter(img);
}

inline bool select_star_database(const std::string& requested, double fov) {
    return astap::reference::select_star_database(requested, fov);
}

inline void plot_stars_used_for_solving(const StarList&, const StarList&,
                                        const Header&, double, double) {}

inline void update_integer(std::vector<std::string>& m, const std::string& k,
                           const std::string& c, int v) {
    astap::core::update_integer(m, k, c, v);
}
inline void update_float(std::vector<std::string>& m, const std::string& k,
                         const std::string& c, bool p, double v) {
    astap::core::update_float(m, k, c, p, v);
}
inline void update_text(std::vector<std::string>& m, const std::string& k,
                        const std::string& v) {
    astap::core::update_text(m, k, v);
}
inline void update_longstr(std::vector<std::string>& m, const std::string& k,
                           const std::string& v) {
    astap::core::update_longstr(m, k, v);
}
inline void remove_key(std::vector<std::string>& m, const std::string& k, bool all) {
    astap::core::remove_key(m, k, all);
}

void populate_stackcfg_from_globals() {
    stackcfg.star_database = astap::reference::name_database;
    stackcfg.radius_search_deg = std::format("{}", astap::search_radius_deg);
    stackcfg.search_fov_deg = std::format("{}", astap::search_fov_deg);
    stackcfg.max_stars = std::format("{}", astap::max_stars_setting);
    stackcfg.quad_tolerance = std::format("{}", astap::quad_tolerance);
    stackcfg.min_star_size = std::format("{}", astap::min_star_size_arcsec);
    stackcfg.downsample_for_solving_index = astap::downsample_setting;
    stackcfg.force_oversize = astap::force_oversize;
    stackcfg.add_sip = astap::add_sip;
}

namespace {

///----------------------------------------
/// MARK: Small helpers
///----------------------------------------

constexpr auto kPi = std::numbers::pi;

inline void sincos_d(double a, double& s, double& c) {
    s = std::sin(a);
    c = std::cos(a);
}

inline void memo2_message(std::vector<std::string>& log, const std::string& msg) {
    log.push_back(msg);
}

///----------------------------------------
/// @brief Convert an angular distance (radians) to a human-readable string.
/// @details @p dist selects the unit ("\"", "'", or "d"); @p inp is the value actually printed.
/// @param dist Angular magnitude used to pick the unit.
/// @param inp Value to print in the selected unit.
/// @return Formatted distance string.
///----------------------------------------

[[nodiscard]] std::string distance_to_string(double dist, double inp) {
    const auto abs_d = std::abs(dist);
    if (abs_d < kPi / (180.0 * 60.0)) {
        return std::format("{:.1f}\"", inp * 3600.0 * 180.0 / kPi);
    }
    if (abs_d < kPi / 180.0) {
        return std::format("{:.1f}'", inp * 60.0 * 180.0 / kPi);
    }
    return std::format("{:.1f}d", inp * 180.0 / kPi);
}

///----------------------------------------
/// @brief Inverse of @ref equatorial_standard: CCD standard coordinates -> (ra, dec).
///----------------------------------------

inline void standard_equatorial(double ra0, double dec0,
                                double x, double y,
                                double cdelt,
                                double& ra, double& dec) {
    double sin_dec0 = 0.0;
    double cos_dec0 = 0.0;
    sincos_d(dec0, sin_dec0, cos_dec0);
    x *= cdelt / (3600.0 * 180.0 / kPi);
    y *= cdelt / (3600.0 * 180.0 / kPi);
    
    const auto delta = cos_dec0 - y * sin_dec0;
    // atan2 handles the celestial pole correctly.
    ra  = ra0 + std::atan2(-x, delta);
    dec = std::atan((sin_dec0 + y * cos_dec0) / std::sqrt(x * x + delta * delta));
    
    // Keep RA inside [0, 2*pi) to avoid confusing later direction detection.
    if (ra > kPi * 2.0) {
        ra -= kPi * 2.0;
    }
    if (ra < 0.0) {
        ra += kPi * 2.0;
    }
}

///----------------------------------------
/// MARK: Binning / cropping helpers (private)
///----------------------------------------

void binX3_crop(double crop, const ImageArray& img, ImageArray& img2) {
    const auto nrcolors = static_cast<int>(img.size());
    const auto width5   = static_cast<int>(img[0][0].size());
    const auto height5  = static_cast<int>(img[0].size());
    
    // 1/3 size, cropped.
    const auto w = static_cast<int>(crop * width5  / 3);
    const auto h = static_cast<int>(crop * height5 / 3);
    
    img2.assign(1, std::vector<std::vector<float>>(h, std::vector<float>(w, 0.0f)));
    
    const auto shiftX = static_cast<int>(std::lround(width5  * (1.0 - crop) / 2.0));
    const auto shiftY = static_cast<int>(std::lround(height5 * (1.0 - crop) / 2.0));
    
    for (int fitsY = 0; fitsY < h; ++fitsY) {
        for (int fitsX = 0; fitsX < w; ++fitsX) {
            auto val = 0.0f;
            for (int k = 0; k < nrcolors; ++k) {
                val += (img[k][shiftY + fitsY * 3    ][shiftX + fitsX * 3    ] +
                        img[k][shiftY + fitsY * 3    ][shiftX + fitsX * 3 + 1] +
                        img[k][shiftY + fitsY * 3    ][shiftX + fitsX * 3 + 2] +
                        img[k][shiftY + fitsY * 3 + 1][shiftX + fitsX * 3    ] +
                        img[k][shiftY + fitsY * 3 + 1][shiftX + fitsX * 3 + 1] +
                        img[k][shiftY + fitsY * 3 + 1][shiftX + fitsX * 3 + 2] +
                        img[k][shiftY + fitsY * 3 + 2][shiftX + fitsX * 3    ] +
                        img[k][shiftY + fitsY * 3 + 2][shiftX + fitsX * 3 + 1] +
                        img[k][shiftY + fitsY * 3 + 2][shiftX + fitsX * 3 + 2]) / 9.0f;
            }
            img2[0][fitsY][fitsX] = val / static_cast<float>(nrcolors);
        }
    }
}

void binX4_crop(double crop, const ImageArray& img, ImageArray& img2) {
    const auto nrcolors = static_cast<int>(img.size());
    const auto width5   = static_cast<int>(img[0][0].size());
    const auto height5  = static_cast<int>(img[0].size());
    
    // 1/4 size, cropped.
    const auto w = static_cast<int>(crop * width5  / 4);
    const auto h = static_cast<int>(crop * height5 / 4);
    
    img2.assign(1, std::vector<std::vector<float>>(h, std::vector<float>(w, 0.0f)));
    
    const auto shiftX = static_cast<int>(std::lround(width5  * (1.0 - crop) / 2.0));
    const auto shiftY = static_cast<int>(std::lround(height5 * (1.0 - crop) / 2.0));
    
    for (int fitsY = 0; fitsY < h; ++fitsY) {
        for (int fitsX = 0; fitsX < w; ++fitsX) {
            auto val = 0.0f;
            for (int k = 0; k < nrcolors; ++k) {
                auto sum = 0.0f;
                for (int dy = 0; dy < 4; ++dy) {
                    for (int dx = 0; dx < 4; ++dx) {
                        sum += img[k][shiftY + fitsY * 4 + dy][shiftX + fitsX * 4 + dx];
                    }
                }
                val += sum / 16.0f;
            }
            img2[0][fitsY][fitsX] = val / static_cast<float>(nrcolors);
        }
    }
}

///----------------------------------------
///      @brief Build a grid of (x, y) positions equally spread over the image,
///             given relative to the image centre.
///    @details Used for SIP grid generation.
///      @param width2 Image width.
///      @param height2 Image height.
///      @param nrpoints Grid divisions per axis; produces @c nrpoints^2 points.
/// @param[out] grid_list Populated with @c nrpoints*nrpoints positions.
///----------------------------------------

void create_grid_list(int width2, int height2, int nrpoints, StarArray& grid_list) {
    const auto middleX = width2  / 2.0;
    const auto middleY = height2 / 2.0;
    grid_list.resize(static_cast<std::size_t>(nrpoints) * nrpoints);
    auto counter = 0;
    for (int y = 0; y < nrpoints; ++y) {
        for (int x = 0; x < nrpoints; ++x) {
            grid_list[counter].x = -middleX + x * width2  / double(nrpoints - 1);
            grid_list[counter].y = -middleY + y * height2 / double(nrpoints - 1);
            ++counter;
        }
    }
}

///----------------------------------------
///  @brief Refine a first-order solve into cubic SIP forward + inverse polynomials.
/// @details After the base solve, refines the mapping between measured and
///          reference stars into cubic SIP polynomials, and stashes the
///          coefficients into the module globals for later FITS serialisation.
///
///          Algorithm:
///          1. Solve image with 1st-order solver.
///          2. Get (x, y) of detected stars  -> @c stars_measured.
///          3. Get (x, y) of reference stars -> @c stars_reference.
///          4. Shift measured coords so [0, 0] is at (CRPIX1, CRPIX2).
///          5. Convert reference coords to the same system as measured coords.
///          6. Now both sets match except for distortion, origin at CRPIX.
///          7. pixel_to_sky: @c calc_trans_cubic(stars_measured, stars_reference).
///          8. sky_to_pixel: @c calc_trans_cubic(stars_reference, stars_measured).
///  @param hd Header used to read the first-order CD matrix and CRPIX.
///  @param memo Log sink.
///  @param ra_database Reference RA used during the search.
///  @param dec_database Reference DEC used during the search.
/// @return @c true on success.
///----------------------------------------

[[nodiscard]] bool add_sip(const Header& hd,
                           std::vector<std::string>& memo,
                           double ra_database, double dec_database) {
    const auto len = b_Xrefpositions.size();
    if (len < 20) {
        memo2_message(memo, "Not enough quads for calculating SIP.");
        return false;
    }
    
    auto stars_measured  = StarArray(len);
    auto stars_reference = StarArray(len);
    
    // Hoisted out of the loop.
    double sin_dec_ref = 0.0;
    double cos_dec_ref = 0.0;
    sincos_d(hd.dec0, sin_dec_ref, cos_dec_ref);
    
    for (std::size_t i = 0; i < len; ++i) {
        // Measured stars: centre at CRPIX1/CRPIX2 in FITS range 1..width.
        stars_measured[i].x = 1.0 + A_XYpositions[0][i] - hd.crpix1;
        stars_measured[i].y = 1.0 + A_XYpositions[1][i] - hd.crpix2;
        
        double ra_t = 0.0;
        double dec_t = 0.0;
        standard_equatorial(ra_database, dec_database,
                            b_Xrefpositions[i],
                            b_Yrefpositions[i],
                            1.0,
                            ra_t, dec_t);
                            
        // (RA, DEC) -> image x, y (FITS range 1..max).
        double sin_dec_t = 0.0;
        double cos_dec_t = 0.0;
        sincos_d(dec_t, sin_dec_t, cos_dec_t);
        
        const auto delta_ra = ra_t - hd.ra0;
        double sin_delta_ra = 0.0;
        double cos_delta_ra = 0.0;
        sincos_d(delta_ra, sin_delta_ra, cos_delta_ra);
        
        const auto H = sin_dec_t * sin_dec_ref + cos_dec_t * cos_dec_ref * cos_delta_ra;
        const auto dRA  = (cos_dec_t * sin_delta_ra / H) * 180.0 / kPi;
        const auto dDEC = ((sin_dec_t * cos_dec_ref - cos_dec_t * sin_dec_ref * cos_delta_ra) / H) * 180.0 / kPi;
        
        const auto det = hd.cd2_2 * hd.cd1_1 - hd.cd1_2 * hd.cd2_1;
        stars_reference[i].x = -(hd.cd1_2 * dDEC - hd.cd2_2 * dRA) / det;
        stars_reference[i].y = +(hd.cd1_1 * dDEC - hd.cd2_1 * dRA) / det;
    }
    
    // sky -> pixel.
    auto r1 = calc_trans_cubic(stars_reference, stars_measured);
    if (!r1) {
        memo2_message(memo, r1.error());
        return false;
    }
    const auto trans_sky_to_pixel = *r1;
    
    // Third-order inverse.
    ap_order = 3;
    ap_0_0 = trans_sky_to_pixel.x00;
    ap_0_1 = trans_sky_to_pixel.x01;
    ap_0_2 = trans_sky_to_pixel.x02;
    ap_0_3 = trans_sky_to_pixel.x03;
    ap_1_0 = -1.0 + trans_sky_to_pixel.x10;
    ap_1_1 = trans_sky_to_pixel.x11;
    ap_1_2 = trans_sky_to_pixel.x12;
    ap_2_0 = trans_sky_to_pixel.x20;
    ap_2_1 = trans_sky_to_pixel.x21;
    ap_3_0 = trans_sky_to_pixel.x30;
    
    bp_0_0 = trans_sky_to_pixel.y00;
    bp_0_1 = -1.0 + trans_sky_to_pixel.y01;
    bp_0_2 = trans_sky_to_pixel.y02;
    bp_0_3 = trans_sky_to_pixel.y03;
    bp_1_0 = trans_sky_to_pixel.y10;
    bp_1_1 = trans_sky_to_pixel.y11;
    bp_1_2 = trans_sky_to_pixel.y12;
    bp_2_0 = trans_sky_to_pixel.y20;
    bp_2_1 = trans_sky_to_pixel.y21;
    bp_3_0 = trans_sky_to_pixel.y30;
    
    // Inverse: pixel -> sky. Swap the arrays; valid because the offset is
    // small in this regime.
    auto r2 = calc_trans_cubic(stars_measured, stars_reference);
    if (!r2) {
        memo2_message(memo, r2.error());
        return false;
    }
    const auto trans_pixel_to_sky = *r2;
    
    sip = true;
    // SIP definitions: https://irsa.ipac.caltech.edu/data/SPITZER/docs/files/spitzer/shupeADASS.pdf
    
    // Pixel -> sky.
    a_order = 3;
    a_0_0 = trans_pixel_to_sky.x00;
    a_0_1 = trans_pixel_to_sky.x01;
    a_0_2 = trans_pixel_to_sky.x02;
    a_0_3 = trans_pixel_to_sky.x03;
    a_1_0 = -1.0 + trans_pixel_to_sky.x10;
    a_1_1 = trans_pixel_to_sky.x11;
    a_1_2 = trans_pixel_to_sky.x12;
    a_2_0 = trans_pixel_to_sky.x20;
    a_2_1 = trans_pixel_to_sky.x21;
    a_3_0 = trans_pixel_to_sky.x30;
    
    b_0_0 = trans_pixel_to_sky.y00;
    b_0_1 = -1.0 + trans_pixel_to_sky.y01;
    b_0_2 = trans_pixel_to_sky.y02;
    b_0_3 = trans_pixel_to_sky.y03;
    b_1_0 = trans_pixel_to_sky.y10;
    b_1_1 = trans_pixel_to_sky.y11;
    b_1_2 = trans_pixel_to_sky.y12;
    b_2_0 = trans_pixel_to_sky.y20;
    b_2_1 = trans_pixel_to_sky.y21;
    b_3_0 = trans_pixel_to_sky.y30;
    
    // Serialise into output memo/header.
    update_integer(memo, "A_ORDER =", " / Polynomial order, axis 1. Pixel to Sky         ", 3);
    update_float  (memo, "A_0_0   =", " / SIP coefficient                                ", false, a_0_0);
    update_float  (memo, "A_1_0   =", " / SIP coefficient                                ", false, a_1_0);
    update_float  (memo, "A_0_1   =", " / SIP coefficient                                ", false, a_0_1);
    update_float  (memo, "A_2_0   =", " / SIP coefficient                                ", false, a_2_0);
    update_float  (memo, "A_1_1   =", " / SIP coefficient                                ", false, a_1_1);
    update_float  (memo, "A_0_2   =", " / SIP coefficient                                ", false, a_0_2);
    update_float  (memo, "A_3_0   =", " / SIP coefficient                                ", false, a_3_0);
    update_float  (memo, "A_2_1   =", " / SIP coefficient                                ", false, a_2_1);
    update_float  (memo, "A_1_2   =", " / SIP coefficient                                ", false, a_1_2);
    update_float  (memo, "A_0_3   =", " / SIP coefficient                                ", false, a_0_3);
    
    update_integer(memo, "B_ORDER =", " / Polynomial order, axis 2. Pixel to sky.        ", 3);
    update_float  (memo, "B_0_0   =", " / SIP coefficient                                ", false, b_0_0);
    update_float  (memo, "B_0_1   =", " / SIP coefficient                                ", false, b_0_1);
    update_float  (memo, "B_1_0   =", " / SIP coefficient                                ", false, b_1_0);
    update_float  (memo, "B_2_0   =", " / SIP coefficient                                ", false, b_2_0);
    update_float  (memo, "B_1_1   =", " / SIP coefficient                                ", false, b_1_1);
    update_float  (memo, "B_0_2   =", " / SIP coefficient                                ", false, b_0_2);
    update_float  (memo, "B_3_0   =", " / SIP coefficient                                ", false, b_3_0);
    update_float  (memo, "B_2_1   =", " / SIP coefficient                                ", false, b_2_1);
    update_float  (memo, "B_1_2   =", " / SIP coefficient                                ", false, b_1_2);
    update_float  (memo, "B_0_3   =", " / SIP coefficient                                ", false, b_0_3);
    
    update_integer(memo, "AP_ORDER=", " / Inv polynomial order, axis 1. Sky to pixel.      ", 3);
    update_float  (memo, "AP_0_0  =", " / SIP coefficient                                ", false, ap_0_0);
    update_float  (memo, "AP_1_0  =", " / SIP coefficient                                ", false, ap_1_0);
    update_float  (memo, "AP_0_1  =", " / SIP coefficient                                ", false, ap_0_1);
    update_float  (memo, "AP_2_0  =", " / SIP coefficient                                ", false, ap_2_0);
    update_float  (memo, "AP_1_1  =", " / SIP coefficient                                ", false, ap_1_1);
    update_float  (memo, "AP_0_2  =", " / SIP coefficient                                ", false, ap_0_2);
    update_float  (memo, "AP_3_0  =", " / SIP coefficient                                ", false, ap_3_0);
    update_float  (memo, "AP_2_1  =", " / SIP coefficient                                ", false, ap_2_1);
    update_float  (memo, "AP_1_2  =", " / SIP coefficient                                ", false, ap_1_2);
    update_float  (memo, "AP_0_3  =", " / SIP coefficient                                ", false, ap_0_3);
    
    update_integer(memo, "BP_ORDER=", " / Inv polynomial order, axis 2. Sky to pixel.    ", 3);
    update_float  (memo, "BP_0_0  =", " / SIP coefficient                                ", false, bp_0_0);
    update_float  (memo, "BP_1_0  =", " / SIP coefficient                                ", false, bp_1_0);
    update_float  (memo, "BP_0_1  =", " / SIP coefficient                                ", false, bp_0_1);
    update_float  (memo, "BP_2_0  =", " / SIP coefficient                                ", false, bp_2_0);
    update_float  (memo, "BP_1_1  =", " / SIP coefficient                                ", false, bp_1_1);
    update_float  (memo, "BP_0_2  =", " / SIP coefficient                                ", false, bp_0_2);
    update_float  (memo, "BP_3_0  =", " / SIP coefficient                                ", false, bp_3_0);
    update_float  (memo, "BP_2_1  =", " / SIP coefficient                                ", false, bp_2_1);
    update_float  (memo, "BP_1_2  =", " / SIP coefficient                                ", false, bp_1_2);
    update_float  (memo, "BP_0_3  =", " / SIP coefficient                                ", false, bp_0_3);
    
    return true;
}

///----------------------------------------
/// @brief Fixed-digit double-to-string convenience wrapper.
///----------------------------------------

[[nodiscard]] inline std::string ff_fixed(double v, int digits) {
    return std::format("{:.{}f}", v, digits);
}

[[nodiscard]] inline double try_stod(const std::string& s, double fallback) {
    try {
        return std::stod(s);
    }
    catch (...) {
        return fallback;
    }
}

[[nodiscard]] inline int try_stoi(const std::string& s, int fallback) {
    try {
        return std::stoi(s);
    }
    catch (...) {
        return fallback;
    }
}

[[nodiscard]] inline std::string prepare_ra(double r, const std::string& s)   { return astap::core::prepare_ra(r, s); }
[[nodiscard]] inline std::string prepare_dec(double d, const std::string& s)  { return astap::core::prepare_dec(d, s); }
[[nodiscard]] inline std::string prepare_ra8(double r, const std::string& s)  { return astap::core::prepare_ra8(r, s); }
[[nodiscard]] inline std::string prepare_dec2(double d, const std::string& s) { return astap::core::prepare_dec2(d, s); }

} // namespace

///----------------------------------------
/// MARK: Public API
///----------------------------------------

double position_angle(double ra1, double dec1, double ra0, double dec0) {
    // Meeus, Astronomical Algorithms (1991 46.5 / 1998 48.5).
    double sinDeltaRa = 0.0;
    double cosDeltaRa = 0.0;
    double sinDec0    = 0.0;
    double cosDec0    = 0.0;
    double sinDec1    = 0.0;
    double cosDec1    = 0.0;
    sincos_d(ra1 - ra0, sinDeltaRa, cosDeltaRa);
    sincos_d(dec0,      sinDec0,    cosDec0);
    sincos_d(dec1,      sinDec1,    cosDec1);
    return std::atan2(cosDec1 * sinDeltaRa,
                      sinDec1 * cosDec0 - cosDec1 * sinDec0 * cosDeltaRa);
}

void equatorial_standard(double ra0, double dec0,
                         double ra, double dec,
                         double cdelt,
                         double& xx, double& yy) {
    double sin_dec0    = 0.0;
    double cos_dec0    = 0.0;
    double sin_dec     = 0.0;
    double cos_dec     = 0.0;
    double sin_deltaRA = 0.0;
    double cos_deltaRA = 0.0;
    sincos_d(dec0,     sin_dec0,    cos_dec0);
    sincos_d(dec,      sin_dec,     cos_dec);
    sincos_d(ra - ra0, sin_deltaRA, cos_deltaRA);
    
    // Factor to convert standard coordinates to CCD pixels.
    const auto dv = (cos_dec0 * cos_dec * cos_deltaRA + sin_dec0 * sin_dec)
                  * cdelt / (3600.0 * 180.0 / kPi);
                  
    // Tangent of the angle in RA / DEC.
    xx = -cos_dec * sin_deltaRA / dv;
    yy = -(sin_dec0 * cos_dec * cos_deltaRA - cos_dec0 * sin_dec) / dv;
}

bool read_stars(double telescope_ra, double telescope_dec,
                double search_field,
                int database_type_in, int nrstars_required,
                StarList& starlist, int& nrstars) {
    nrstars = 0;
    
    // Initialise to avoid NaN when no record is read.
    auto ra2   = 0.0;
    auto dec2  = 0.0;
    auto bp_rp = 0.0;
    auto sep   = 0.0;
    
    // Layout: 2 rows, nrstars_required columns.
    starlist.assign(2, std::vector<double>(static_cast<std::size_t>(nrstars_required), 0.0));
    
    if (database_type_in > 1) {
        // 1476 or 290 tile databases.
        // Treat the search field as straddling up to four tiles. Allocate
        // each area's share of the requested stars proportionally.
        auto area1 = 0;
        auto area2 = 0;
        auto area3 = 0;
        auto area4 = 0;
        auto frac1 = 0.0;
        auto frac2 = 0.0;
        auto frac3 = 0.0;
        auto frac4 = 0.0;
        find_areas(telescope_ra, telescope_dec, search_field,
                   area1, area2, area3, area4,
                   frac1, frac2, frac3, frac4);
                   
        auto read_tile = [&](int area, double cumulative_frac) -> bool {
            if (area == 0) {
                return true;
            }
            if (!open_database(telescope_dec, area)) {
                return false;
            }
            const auto nrstars_required2 = std::min(
                nrstars_required,
                static_cast<int>(nrstars_required * cumulative_frac));
            while (nrstars < nrstars_required2 &&
                   readdatabase290(telescope_ra, telescope_dec, search_field,
                                   ra2, dec2, mag2, bp_rp)) {
                equatorial_standard(telescope_ra, telescope_dec,
                                    ra2, dec2, 1.0,
                                    starlist[0][nrstars], starlist[1][nrstars]);
                ++nrstars;
            }
            return true;
        };
        
        if (!read_tile(area1, frac1)) {
            return false;
        }
        if (!read_tile(area2, frac1 + frac2)) {
            return false;
        }
        if (!read_tile(area3, frac1 + frac2 + frac3)) {
            return false;
        }
        if (!read_tile(area4, frac1 + frac2 + frac3 + frac4)) {
            return false;
        }
    }
    else {
        // Wide field database (database_type == 1).
        if (wide_database != name_database) {
            read_stars_wide_field();
        }
        auto count = 0;
        cos_telescope_dec = std::cos(telescope_dec);
        const auto nmax = static_cast<int>(wide_field_stars.size() / 3);
        while (nrstars < nrstars_required && count < nmax) {
            ra2  = wide_field_stars[count * 3 + 1];
            dec2 = wide_field_stars[count * 3 + 2];
            ang_sep(ra2, dec2, telescope_ra, telescope_dec, sep);
            // 2/sqrt(pi) adapts circle search to square, 0.9 is a fudge for
            // trees/houses/dark corners; sep<pi/2 bounds equatorial_standard.
            if (sep < search_field * 0.5 * 0.9 * (2.0 / std::sqrt(kPi)) &&
                sep < kPi / 2.0) {
                equatorial_standard(telescope_ra, telescope_dec,
                                    ra2, dec2, 1.0,
                                    starlist[0][nrstars], starlist[1][nrstars]);
                ++nrstars;
            }
            ++count;
        }
        if (count > 0) {
            // Faintest magnitude actually used.
            mag2 = wide_field_stars[(count - 1) * 3];
        }
    }
    
    if (nrstars < nrstars_required) {
        // Trim trailing allocated columns.
        starlist[0].resize(nrstars);
        starlist[1].resize(nrstars);
    }
    return true;
}

void binX1_crop(double crop, const ImageArray& img, ImageArray& img2) {
    const auto nrcolors = static_cast<int>(img.size());
    const auto width5   = static_cast<int>(img[0][0].size());
    const auto height5  = static_cast<int>(img[0].size());
    
    const auto w = static_cast<int>(crop * width5);
    const auto h = static_cast<int>(crop * height5);
    
    img2.assign(1, std::vector<std::vector<float>>(h, std::vector<float>(w, 0.0f)));
    
    const auto shiftX = static_cast<int>(std::lround(width5  * (1.0 - crop) / 2.0));
    const auto shiftY = static_cast<int>(std::lround(height5 * (1.0 - crop) / 2.0));
    
    for (int fitsY = 0; fitsY < h; ++fitsY) {
        for (int fitsX = 0; fitsX < w; ++fitsX) {
            auto val = 0.0f;
            for (int k = 0; k < nrcolors; ++k) {
                val += img[k][shiftY + fitsY][shiftX + fitsX];
            }
            img2[0][fitsY][fitsX] = val / static_cast<float>(nrcolors);
        }
    }
}

void binX2_crop(double crop, const ImageArray& img, ImageArray& img2) {
    const auto nrcolors = static_cast<int>(img.size());
    const auto width5   = static_cast<int>(img[0][0].size());
    const auto height5  = static_cast<int>(img[0].size());
    
    // Half size and cropped. trunc avoids an exception for odd-width images
    // like the 1391-px-wide M27 test.
    const auto w = static_cast<int>(crop * width5  / 2);
    const auto h = static_cast<int>(crop * height5 / 2);
    
    img2.assign(1, std::vector<std::vector<float>>(h, std::vector<float>(w, 0.0f)));
    
    const auto shiftX = static_cast<int>(std::lround(width5  * (1.0 - crop) / 2.0));
    const auto shiftY = static_cast<int>(std::lround(height5 * (1.0 - crop) / 2.0));
    
    for (int fitsY = 0; fitsY < h; ++fitsY) {
        for (int fitsX = 0; fitsX < w; ++fitsX) {
            auto val = 0.0f;
            for (int k = 0; k < nrcolors; ++k) {
                val += (img[k][shiftY + fitsY * 2    ][shiftX + fitsX * 2    ] +
                        img[k][shiftY + fitsY * 2 + 1][shiftX + fitsX * 2    ] +
                        img[k][shiftY + fitsY * 2    ][shiftX + fitsX * 2 + 1] +
                        img[k][shiftY + fitsY * 2 + 1][shiftX + fitsX * 2 + 1]) / 4.0f;
            }
            img2[0][fitsY][fitsX] = val / static_cast<float>(nrcolors);
        }
    }
}

void bin_and_find_stars(const ImageArray& img,
                        int binning, double cropping, double hfd_min,
                        double hfd_max,
                        int max_stars, bool get_hist,
                        StarList& starlist3, std::string& short_warning) {
    short_warning.clear();
    
    const auto width5  = static_cast<int>(img[0][0].size());
    const auto height5 = static_cast<int>(img[0].size());
    
    // Local discard; user-facing messages are emitted from solve_image.
    auto dummy_log = std::vector<std::string>{};
    
    if (binning > 1 || cropping < 1.0) {
        if (binning > 1) {
            memo2_message(dummy_log,
                "Creating grayscale x " + std::to_string(binning) +
                " binning image for solving or star alignment.");
        }
        if (cropping != 1.0) {
            memo2_message(dummy_log, "Cropping image x " + ff_fixed(cropping, 2));
        }
        
        auto img_binned = ImageArray{};
        if (binning == 2) {
            binX2_crop(cropping, img, img_binned);
        }
        else if (binning == 3) {
            binX3_crop(cropping, img, img_binned);
        }
        else if (binning == 4) {
            binX4_crop(cropping, img, img_binned);
        }
        else if (binning == 1) {
            binX1_crop(cropping, img, img_binned);
        }
        
        get_background(0, img_binned, true, true, bck);
        find_stars(img_binned, hfd_min, hfd_max, max_stars, bck, starlist3);
        
        if (static_cast<int>(img_binned[0].size()) < 960) {
            short_warning = "Warning, remaining image dimensions too low! ";
            memo2_message(dummy_log,
                "Warning, remaining image dimensions too low! Try to REDUCE OR "
                "REMOVE DOWNSAMPLING. Set this option in stack menu, tab alignment.");
        }
        img_binned.clear();
        
        const auto nrstars = static_cast<int>(starlist3[0].size());
        for (int i = 0; i < nrstars; ++i) {
            // Correct star positions for binning/cropping. After e.g. 2x2
            // binning, position [3.5, 3.5] -> [1, 1]; multiplying back by 2
            // lands at [3, 3], a half-pixel shift. General rule:
            //   x := (binning-1)*0.5 + binning*x + shiftX
            starlist3[0][i] = (binning - 1) * 0.5 + starlist3[0][i] * binning
                            + (width5  * (1.0 - cropping) / 2.0);
            starlist3[1][i] = (binning - 1) * 0.5 + starlist3[1][i] * binning
                            + (height5 * (1.0 - cropping) / 2.0);
        }
    }
    else {
        if (height5 > 2500) {
            short_warning = "Warning, increase downsampling!! ";
            memo2_message(dummy_log,
                "Info: DOWNSAMPLING IS RECOMMENDED FOR LARGE IMAGES. "
                "Set this option in stack menu, tab alignment.");
        }
        else if (height5 < 960) {
            short_warning = "Warning, small image dimensions! ";
            memo2_message(dummy_log, "Warning, small image dimensions!");
        }
        
        get_background(0, img, get_hist, true, bck);
        find_stars(img, hfd_min, hfd_max, max_stars, bck, starlist3);
    }
}

int report_binning(double height) {
    auto result = stackcfg.downsample_for_solving_index;
    // 0 = Auto, -1 = none.
    if (result <= 0) {
        result = (height > 2500) ? 2 : 1;
    }
    return result;
}

///----------------------------------------
/// MARK: Main plate solver
///----------------------------------------

bool solve_image(ImageArray& img, Header& hd,
                 std::vector<std::string>& memo,
                 bool get_hist, bool check_patternfilter) {
    // Sync solver settings from engine globals (set by CLI arg parser or
    // GUI dialog) into the file-local StackMenuConfig that the solver reads.
    populate_stackcfg_from_globals();

    // Sync RA/Dec search centre from the header's position hint.
    ra_radians = hd.ra0;
    dec_radians = hd.dec0;

    // Verbose per-position logging. Controlled by the global flag which the
    // GUI could expose as a checkbox in the Solve dialog later.
    solve_show_log = astap::commandline_log;

    esc_pressed = false;
    warning_str.clear();
    const auto startTick = std::chrono::steady_clock::now();
    auto popup_warningG05 = std::string{};

    if (check_patternfilter) {
        check_pattern_filter(img);
        get_hist = true;
    }
    
    // Cap CLI-specified quad tolerance.
    auto quad_tolerance = try_stod(stackcfg.quad_tolerance, 0.0);
    quad_tolerance = std::min(quad_tolerance, 0.01);
    
    auto max_stars = try_stoi(stackcfg.max_stars, 500);
    const auto use_triples = stackcfg.use_triples;
    
    const auto ra_start  = ra_radians;
    const auto dec_start = dec_radians;
    
    auto fov_org = 0.0;
    if (!fov_specified && hd.cdelt2 != 0.0) {
        fov_org = std::min(180.0, hd.height * std::abs(hd.cdelt2));
    }
    else {
        fov_org = std::min(180.0, try_stod(stackcfg.search_fov_deg, 0.0));
    }
    
    if (!select_star_database(stackcfg.star_database, fov_org)) {
        errorlevel = 32;
        return false;
    }
    {
        auto upper = name_database;
        std::transform(upper.begin(), upper.end(), upper.begin(),
                       [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
        memo2_message(memo, "Using star database " + upper);
    }
    
    if (fov_org > 30 && database_type != 1) {
        warning_str += "Very large FOV, use W08 database! ";
    }
    else if (fov_org > 6 && database_type == 1476) {
        warning_str += "Large FOV, use G05 (or V05) database! ";
    }
    if (!warning_str.empty()) {
        memo2_message(memo, warning_str);
    }
    popup_warningG05 = "\n" + warning_str;
    
    auto max_fov = 0.0;
    if (database_type == 1476) {
        // 1476 tile size.
        max_fov = 5.142857143;
    }
    else if (database_type == 290) {
        // 290 tile size.
        max_fov = 9.53;
    }
    else {
        max_fov = 180.0;
    }
    
    if (max_stars == 0) {
        max_stars = 500;
    }
    
    auto database_density = 0;
    {
        // Parse two chars starting at offset 2 (1-based) of @c name_database.
        auto parsed = 0;
        auto ok = false;
        if (name_database.size() >= 3) {
            try {
                parsed = std::stoi(name_database.substr(1, 2));
                ok = true;
            }
            catch (...) {
                ok = false;
            }
        }
        if (!ok || parsed == 17 || parsed == 18) {
            // Old V17, G17, G18, H17, H18.
            database_density = 9999;
        }
        else {
            database_density = parsed * 100;
        }
    }
    
    const auto min_star_size_arcsec = try_stod(stackcfg.min_star_size, 0.0);
    const auto autoFOV              = (fov_org == 0.0);
    auto fov_min                    = 0.0;
    
    auto starlist1 = StarList{};
    auto starlist2 = StarList{};
    auto nrstars              = 0;
    auto nrstars_required     = 0;
    auto nrstars_required2    = 0;
    auto database_stars       = 0;
    auto binning              = 0;
    auto nr_quads             = 0;
    auto minimum_quads        = 3;
    auto match_nr             = 0;
    auto ra_database          = 0.0;
    auto dec_database         = 0.0;
    auto cropping             = 1.0;
    auto fov2                 = 0.0;
    auto hfd_min              = 0.0;
    auto oversize             = 1.0;
    auto oversize2            = 1.0;
    auto radius               = 0.0;
    auto centerX              = 0.0;
    auto centerY              = 0.0;
    auto correctionX          = 0.0;
    auto correctionY          = 0.0;
    auto sep_search           = 0.0;
    auto solution             = false;
    auto go_ahead             = false;
    auto yes_use_triples      = false;
    auto quads_str            = std::string{" quads"};
    auto warning_downsample   = std::string{};
    
    // AutoFOV loop.
    do {
        if (autoFOV) {
            if (fov_org == 0.0) {
                if (database_type != 1) {
                    fov_org = 9.5;
                    fov_min = 0.38;
                }
                else {
                    fov_org = 90.0;
                    fov_min = 12.0;
                }
            }
            else {
                fov_org /= 1.5;
            }
            memo2_message(memo, "Trying FOV: " + ff_fixed(fov_org, 1));
        }
        
        if (fov_org > max_fov) {
            cropping = max_fov / fov_org;
            fov2     = max_fov;
        }
        else {
            cropping = 1.0;
            fov2     = fov_org;
        }
        
        // limit = density * surface_of_full_image.
        auto limit = static_cast<int>(std::lround(
            database_density * fov2 * fov2 * hd.width / hd.height));
        if (limit < max_stars) {
            max_stars = limit;
            memo2_message(memo, "Database limit for this FOV is " +
                                std::to_string(max_stars) + " stars.");
        }
        
        binning = report_binning(hd.height * cropping);
        hfd_min = std::max(0.8,
            min_star_size_arcsec / (binning * fov_org * 3600.0 / hd.height));
            
        bin_and_find_stars(img, binning, cropping, hfd_min, /*hfd_max=*/10.0,
                           max_stars, get_hist, starlist2, warning_downsample);
        nrstars = static_cast<int>(starlist2[0].size());
        
        if (hd.xpixsz != 0.0 && hd.ypixsz != 0.0 &&
            std::abs(hd.xpixsz - hd.ypixsz) > 0.1) {
            // Non-square pixels (extremely rare: e.g. QHY6).
            memo2_message(memo, "Rare none square pixels specified.");
            const auto pixel_aspect_ratio = hd.xpixsz / hd.ypixsz;
            for (int i = 0; i < nrstars; ++i) {
                starlist2[0][i] = hd.width / 2.0 +
                                  (starlist2[0][i] - hd.width / 2.0) * pixel_aspect_ratio;
            }
        }
        
        // Informational strings; would populate tray popup.
        [[maybe_unused]] const auto popup_warningSample =
            warning_downsample.empty() ? std::string{} : "\n" + warning_downsample;
        [[maybe_unused]] const auto keep_popup_warningG05 = &popup_warningG05;
        
        nrstars_required = static_cast<int>(std::lround(
            nrstars * (static_cast<double>(hd.height) / hd.width)));
            
        solution = false;
        // Bare minimum for three quads.
        go_ahead = (nrstars >= 6);
        
        if (go_ahead) {
            yes_use_triples = (nrstars < 30) && use_triples;
            if (yes_use_triples) {
                find_triples_using_quads(starlist2, quad_star_distances2);
                quad_tolerance = 0.002;
                quads_str = " triples";
                if (solve_show_log) {
                    memo2_message(memo,
                        "For triples the hash code tolerance is forced to " +
                        ff_fixed(quad_tolerance, 6) + ".");
                }
            }
            else {
                find_quads(starlist2, quad_star_distances2);
                quads_str = " quads";
            }
            
            nr_quads = static_cast<int>(quad_star_distances2[0].size());
            go_ahead = (nr_quads >= 3);
            
            // Grow search area up to 2x for low quad counts to guarantee
            // every image quad is covered as we step through the sky.
            if (nr_quads < 25) {
                oversize = 2.0;
            }
            else if (nr_quads > 100) {
                oversize = 1.0;
            }
            else {
                oversize = 2.0 * std::sqrt(25.0 / nr_quads);
            }
            
            if (stackcfg.force_oversize) {
                oversize = 2.0;
            }
            
            oversize = std::min(oversize, max_fov / fov2);
            radius = try_stod(stackcfg.radius_search_deg, 0.0);
            // Harder threshold for dense images.
            minimum_quads = 3 + nr_quads / 100;
        }
        else {
            memo2_message(memo, "Only " + std::to_string(nrstars) +
                                " stars found in image. Abort");
            errorlevel = 2;
        }
        
        if (go_ahead) {
            const auto search_field_init = fov2 * (kPi / 180.0);
            auto STEP_SIZE = search_field_init;
            auto max_distance = 0;
            
            if (database_type == 1) {
                // Small steps for reliable wide-field solving.
                STEP_SIZE *= 0.1;
                max_distance = static_cast<int>(std::lround(radius / (0.1 * fov2 + 0.00001)));
                memo2_message(memo, "Wide field, making small steps for reliable solving.");
            }
            else {
                max_distance = static_cast<int>(std::lround(radius / (fov2 + 0.00001)));
            }
            
            memo2_message(memo,
                std::to_string(nrstars) + " stars, " +
                std::to_string(nr_quads) + quads_str + " selected in the image. " +
                std::to_string(nrstars_required) + " database stars, " +
                std::to_string(static_cast<int>(std::lround(
                    double(nr_quads) * nrstars_required / nrstars))) +
                " database" + quads_str + " required for the " +
                ff_fixed(oversize * fov2, 2) + " degree square search window. " +
                "Step size " + ff_fixed(fov2, 2) + " degree.");
                
            match_nr = 0;
            // Maximum-accuracy loop: refit from the centre after first hit.
            do {
                auto count     = 0;
                auto distance  = 0.0;
                auto spiral_x  = 0;
                auto spiral_y  = 0;
                auto spiral_dx = 0;
                auto spiral_dy = -1;
                auto spiral_t  = 0;
                
                // Squared-spiral search.
                do {
                    if (count != 0) {
                        // Walk [0,0], [1,0], [1,1], [0,1], [-1,1], [-1,0],
                        // [-1,-1], [0,-1], [1,-1], [2,-1], [2,0], ...
                        if ((spiral_x == spiral_y) ||
                            (spiral_x <  0 && spiral_x == -spiral_y) ||
                            (spiral_x >  0 && spiral_x == 1 - spiral_y)) {
                            spiral_t  = spiral_dx;
                            spiral_dx = -spiral_dy;
                            spiral_dy = spiral_t;
                        }
                        spiral_x += spiral_dx;
                        spiral_y += spiral_dy;
                    }
                    
                    dec_database = STEP_SIZE * spiral_y + dec_radians;
                    auto flip = 0.0;
                    if (dec_database > +kPi / 2.0) {
                        // Crossed the pole.
                        dec_database = kPi - dec_database;
                        flip = kPi;
                    }
                    else if (dec_database < -kPi / 2.0) {
                        dec_database = -kPi - dec_database;
                        flip = kPi;
                    }
                    
                    const auto extra = (dec_database > 0.0) ? STEP_SIZE / 2.0
                                                            : -STEP_SIZE / 2.0;
                    const auto ra_database_offset = STEP_SIZE * spiral_x
                                                  / std::cos(dec_database - extra);
                                                  
                    if (ra_database_offset <= +kPi / 2.0 + STEP_SIZE / 2.0 &&
                        ra_database_offset >= -kPi / 2.0) {
                        ra_database = fnmodulo(flip + ra_radians + ra_database_offset,
                                               2.0 * kPi);
                        auto seperation = 0.0;
                        ang_sep(ra_database, dec_database, ra_radians, dec_radians,
                                seperation);
                                
                        if (seperation <= radius * kPi / 180.0 + STEP_SIZE / 2.0) {
                            if (seperation * 180.0 / kPi > distance + fov_org) {
                                distance = seperation * 180.0 / kPi;
                                [[maybe_unused]] const auto distancestr =
                                    std::to_string(static_cast<int>(std::lround(distance))) + "d";
                            }
                            
                            if (match_nr == 0) {
                                oversize2 = oversize;
                            }
                            else {
                                // Second solve: use full image, but limit to one tile.
                                const auto full = std::sqrt(
                                    double(hd.width) * hd.width /
                                    (double(hd.height) * hd.height) + 1.0);
                                oversize2 = std::min(max_fov / fov2,
                                                     std::max(oversize, full));
                            }
                            nrstars_required2 = static_cast<int>(std::lround(
                                nrstars_required * oversize2 * oversize2));
                                
                            if (!read_stars(ra_database, dec_database,
                                            search_field_init * oversize2,
                                            database_type, nrstars_required2,
                                            starlist1, database_stars)) {
                                memo2_message(memo,
                                    "No star database found at " + database_path.string() +
                                    " ! Download and install one star database.");
                                errorlevel = 33;
                                return false;
                            }
                            
                            if (yes_use_triples) {
                                find_triples_using_quads(starlist1, quad_star_distances1);
                            }
                            else {
                                find_quads(starlist1, quad_star_distances1);
                            }
                            
                            if (solve_show_log) {
                                memo2_message(memo,
                                    "Search " + std::to_string(count) +
                                    ", [" + std::to_string(spiral_x) + "," +
                                    std::to_string(spiral_y) + "]\tposition: \t" +
                                    prepare_ra(ra_database, ": ") + "\t" +
                                    prepare_dec(dec_database, " deg ") +
                                    "\t Down to magn " + ff_fixed(mag2 / 10.0, 1) +
                                    "\t " + std::to_string(database_stars) + " database stars" +
                                    "\t " + std::to_string(quad_star_distances1[0].size()) +
                                    " database quads to compare.");
                            }
                            
                            solution = find_offset_and_rotation(minimum_quads, quad_tolerance);

                            if (solve_show_log && count <= 3) {
                                memo2_message(memo,
                                    "  match diag: q1=" + std::to_string(diag_nrquads1) +
                                    " q2=" + std::to_string(diag_nrquads2) +
                                    " pass1=" + std::to_string(diag_pass1_matches) +
                                    " pass2=" + std::to_string(diag_pass2_matches));
                            }

                            if (esc_pressed) {
                                return false;
                            }
                        }
                    }
                    ++count;
                } while (!solution && spiral_x <= max_distance);
                
                if (solution) {
                    centerX = (hd.width  - 1) / 2.0;
                    centerY = (hd.height - 1) / 2.0;
                    
                    // Plug solved centre equatorial position back into the
                    // ra/dec start-search to refit on the second pass.
                    standard_equatorial(ra_database, dec_database,
                        solution_vectorX[0] * centerX +
                        solution_vectorX[1] * centerY +
                        solution_vectorX[2],
                        solution_vectorY[0] * centerX +
                        solution_vectorY[1] * centerY +
                        solution_vectorY[2],
                        1.0,
                        ra_radians, dec_radians);
                    ++match_nr;
                }
                else {
                    // Second solve unlikely, but reset anyway.
                    match_nr = 0;
                }
            } while (solution && match_nr < 2);
        }
    } while (autoFOV && !solution && fov2 > fov_min);
    
    if (solution) {
        hd.ra0    = ra_radians;
        hd.dec0   = dec_radians;
        // FITS convention: 1-based.
        hd.crpix1 = centerX + 1;
        hd.crpix2 = centerY + 1;
        
        ang_sep(ra_radians, dec_radians, ra_start, dec_start, sep_search);
        
        memo2_message(memo,
            std::to_string(nr_references) + " of " + std::to_string(nr_references2) +
            quads_str + " selected matching within " + ff_fixed(quad_tolerance, 6) +
            " tolerance.  Solution[\"] x:=" + ff_fixed(solution_vectorX[0], 6) + "x+ " +
            ff_fixed(solution_vectorX[1], 6) + "y+ " + ff_fixed(solution_vectorX[2], 6) +
            ",  y:=" + ff_fixed(solution_vectorY[0], 6) + "x+ " +
            ff_fixed(solution_vectorY[1], 6) + "y+ " + ff_fixed(solution_vectorY[2], 6));
            
        // New-in-2023 method for correct rotation near the celestial pole.
        auto flipped_image = 0.0;
        if (solution_vectorX[0] * solution_vectorY[1] -
            solution_vectorX[1] * solution_vectorY[0] > 0.0) {
            // Flipped vertically OR horizontally (not both).
            flipped_image = -1.0;
        }
        else {
            flipped_image = +1.0;
        }
        
        // Position +1 pixel along CRPIX2 direction.
        auto ra7  = 0.0;
        auto dec7 = 0.0;
        standard_equatorial(ra_database, dec_database,
            solution_vectorX[0] * centerX + solution_vectorX[1] * (centerY + 1) +
            solution_vectorX[2],
            solution_vectorY[0] * centerX + solution_vectorY[1] * (centerY + 1) +
            solution_vectorY[2],
            1.0, ra7, dec7);
            
        auto crota2 = -position_angle(ra7, dec7, hd.ra0, hd.dec0);
        
        // Position flipped_image pixels along CRPIX1 direction.
        standard_equatorial(ra_database, dec_database,
            solution_vectorX[0] * (centerX + flipped_image) +
            solution_vectorX[1] * centerY + solution_vectorX[2],
            solution_vectorY[0] * (centerX + flipped_image) +
            solution_vectorY[1] * centerY + solution_vectorY[2],
            1.0, ra7, dec7);
            
        auto crota1 = kPi / 2.0 - position_angle(ra7, dec7, hd.ra0, hd.dec0);
        // Wrap into [-pi, +pi].
        if (crota1 > kPi) {
            crota1 -= 2.0 * kPi;
        }
        
        hd.cdelt1 = flipped_image * std::sqrt(
            solution_vectorX[0] * solution_vectorX[0] +
            solution_vectorX[1] * solution_vectorX[1]) / 3600.0;
        hd.cdelt2 = std::sqrt(
            solution_vectorY[0] * solution_vectorY[0] +
            solution_vectorY[1] * solution_vectorY[1]) / 3600.0;
            
        hd.cd1_1 = +hd.cdelt1 * std::cos(crota1);
        hd.cd1_2 = -hd.cdelt1 * std::sin(crota1) * flipped_image;
        hd.cd2_1 = +hd.cdelt2 * std::sin(crota2) * flipped_image;
        hd.cd2_2 = +hd.cdelt2 * std::cos(crota2);
        
        hd.crota2 = crota2 * 180.0 / kPi;
        hd.crota1 = crota1 * 180.0 / kPi;
        
        const auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - startTick).count();
        const auto solved_in = " Solved in " +
            ff_fixed(std::round(elapsed_ms / 100.0) / 10.0, 1) + " sec.";
            
        const auto offset_found = distance_to_string(sep_search, sep_search) + ".";
        
        auto mount_info_str   = std::string{};
        auto mount_offset_str = std::string{};
        if (ra_mount < 99) {
            // Total mount error (used only for scaling the distance unit).
            const auto mount_ra_sep  = kPi * (std::fmod((ra_mount - ra_radians) / kPi, 1.0))
                                       * std::cos((dec_mount + dec_radians) * 0.5);
            const auto mount_dec_sep = dec_mount - dec_radians;
            const auto mount_sep = std::sqrt(
                mount_ra_sep * mount_ra_sep + mount_dec_sep * mount_dec_sep);
                
            const auto ra_offset_str  = distance_to_string(mount_sep, mount_ra_sep);
            const auto dec_offset_str = distance_to_string(mount_sep, mount_dec_sep);
            mount_offset_str = " Mount offset RA=" + ra_offset_str +
                               ", DEC=" + dec_offset_str;
            mount_info_str = " Mount dRA=" + ra_offset_str +
                             ",  dDEC=" + dec_offset_str + ". \t";
        }
        
        memo2_message(memo,
            "Solution found: " + prepare_ra8(hd.ra0, ": ") + "\t" +
            prepare_dec2(hd.dec0, " deg ") + "\t" + solved_in + "\t delta was " +
            offset_found + "\t" + mount_info_str +
            " Used stars down to magnitude: " + ff_fixed(mag2 / 10.0, 1));
            
        // SIP distortion (optional). About 50 ms is spent on header update.
        if (stackcfg.add_sip && add_sip(hd, memo, ra_database, dec_database)) {
            update_text(memo, "CTYPE1  =",
                        "'RA---TAN-SIP'       / TAN (gnomic) projection + SIP distortions      ");
            update_text(memo, "CTYPE2  =",
                        "'DEC--TAN-SIP'       / TAN (gnomic) projection + SIP distortions      ");
        }
        else {
            update_text(memo, "CTYPE1  =",
                        "'RA---TAN'           / first parameter RA,    projection TANgential   ");
            update_text(memo, "CTYPE2  =",
                        "'DEC--TAN'           / second parameter DEC,  projection TANgential   ");
        }
        
        update_text(memo, "CUNIT1  =",
                    "'deg     '           / Unit of coordinates                            ");
        update_text(memo, "EQUINOX =",
                    "              2000.0 / Equinox of coordinates                         ");
                    
        update_float(memo, "CRPIX1  =", " / X of reference pixel                           ", false, hd.crpix1);
        update_float(memo, "CRPIX2  =", " / Y of reference pixel                           ", false, hd.crpix2);
        update_float(memo, "CRVAL1  =", " / RA of reference pixel (deg)                    ", false, hd.ra0  * 180.0 / kPi);
        update_float(memo, "CRVAL2  =", " / DEC of reference pixel (deg)                   ", false, hd.dec0 * 180.0 / kPi);
        update_float(memo, "CDELT1  =", " / X pixel size (deg)                             ", false, hd.cdelt1);
        update_float(memo, "CDELT2  =", " / Y pixel size (deg)                             ", false, hd.cdelt2);
        update_float(memo, "CROTA1  =", " / Image twist X axis (deg)                       ", false, hd.crota1);
        update_float(memo, "CROTA2  =", " / Image twist Y axis (deg) E of N if not flipped.", false, hd.crota2);
        update_float(memo, "CD1_1   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, hd.cd1_1);
        update_float(memo, "CD1_2   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, hd.cd1_2);
        update_float(memo, "CD2_1   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, hd.cd2_1);
        update_float(memo, "CD2_2   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, hd.cd2_2);
        update_text (memo, "PLTSOLVD=",
            "                   T / Astrometric solved by ASTAP v" + astap_version + ".       ");
        update_text (memo, "COMMENT 7", solved_in + " Offset " + offset_found + mount_offset_str);
        
        if (solve_show_log) {
            equatorial_standard(ra_database, dec_database, hd.ra0, hd.dec0, 1.0,
                                correctionX, correctionY);
            plot_stars_used_for_solving(starlist1, starlist2, hd,
                                        correctionX, correctionY);
            memo2_message(memo,
                "See viewer image for image stars used (red) and database star used (yellow)");
        }
        
        if (fov_org > 1.05 * (hd.height * hd.cdelt2) ||
            fov_org < 0.95 * (hd.height * hd.cdelt2)) {
            // hd.cdelt2 is always positive in ASTAP; no abs() needed.
            auto suggest_str = std::string{};
            if (hd.xpixsz != 0.0) {
                suggest_str = "Warning inexact scale! Set FOV=" +
                              ff_fixed(hd.height * hd.cdelt2, 2) +
                              "d or scale=" + ff_fixed(hd.cdelt2 * 3600.0, 1) +
                              "\"/pix or FL=" +
                              std::to_string(static_cast<int>(std::lround(
                                  180.0 / (kPi * 1000.0) * hd.xpixsz / hd.cdelt2))) +
                              "mm ";
            }
            else {
                suggest_str = "Warning inexact scale! Set FOV=" +
                              ff_fixed(hd.height * hd.cdelt2, 2) +
                              "d or scale=" + ff_fixed(hd.cdelt2 * 3600.0, 1) +
                              "\"/pix ";
            }
            memo2_message(memo, suggest_str);
            warning_str = suggest_str + warning_str;
        }
    }
    else {
        memo2_message(memo, "No solution found!  :(");
        update_text(memo, "PLTSOLVD=",
                    "                   F / No plate solution found.   ");
        remove_key(memo, "COMMENT 7", false);
    }
    
    // Bring in last warning from the AutoFOV loop.
    warning_str += warning_downsample;
    if (!warning_str.empty()) {
        update_longstr(memo, "WARNING =", warning_str);
    }
    return solution;
}
    
} // namespace
