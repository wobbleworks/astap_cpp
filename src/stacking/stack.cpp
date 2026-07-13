///----------------------------------------
///      @file stack.cpp
///   @ingroup ASTAP++
///     @brief Non-GUI algorithmic implementations for the image stacking module.
///    @author Ported from Han Kleijn's unit_stack.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "stack.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <numbers>
#include <span>
#include <string>
#include <vector>

#include "../core/fits.h"
#include "../core/globals.h"
#include "../core/image_io.h"      // save_tiff16_secure (atomic memo→tiff writer)
#include "../core/photometry.h"
#include "../core/util.h"          // smedian, fits_file_name
#include "../core/wcs.h"           // pixel_to_celestial (the real WCS projector)
#include "../solving/astrometric_solving.h"  // the real internal plate solver

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

using astap::head;
using astap::head_ref;
using astap::img_loaded;
using astap::jd_start;
using astap::jd_mid;
using astap::jd_end;

// Defined in stack_routines.cpp (astap::stacking) → forwards to the core sink.
void memo2_message(std::string_view msg);

///----------------------------------------
/// MARK: File-local helpers
///----------------------------------------

namespace {

// Background / histogram / star detection — delegates to the ported
// implementations in astap::core::photometry. Wrappers keep the call sites
// in analyse_image / plot_stars untouched while the rest of the port moves
// to the explicit-parameter style.
inline void get_background(int color,
                           const ImageArray& img,
                           bool measure_hist,
                           bool calc_noise,
                           Background& bck) {
    astap::core::get_background(color, img, measure_hist, calc_noise, bck);
}

// Persistent HFD scratch — star_bg / sd_bg / r_aperture carry over between
// calls in the original source.
astap::core::HfdScratch& hfd_scratch() {
    thread_local astap::core::HfdScratch s;
    return s;
}

inline void hfd_measure(const ImageArray& img,
                        int x, int y,
                        int annulus, int aperture,
                        double adu_e,
                        double& hfd1, double& fwhm, double& snr,
                        double& flux, double& xc, double& yc) {
    astap::core::HfdResult r;
    const auto xbinning = (head.xbinning > 0) ? head.xbinning : 1.0;
    astap::core::HFD(img, x, y, annulus,
                     static_cast<double>(aperture), adu_e, xbinning,
                     r, hfd_scratch());
    hfd1 = r.hfd;
    fwhm = r.fwhm;
    snr  = r.snr;
    flux = r.flux;
    xc   = r.xc;
    yc   = r.yc;
}

void pixel_to_celestial(const Header& hd, double px, double py,
                        int formalism, double& ra, double& decl) {
    // The real projector (TAN / SIP / DSS) lives in core/wcs.cpp; forward to it
    // so the analyse CSV gets true RA/Dec instead of zeros.
    astap::core::pixel_to_celestial(hd, px, py, formalism, ra, decl);
}

/// @brief Catmull-Rom cubic spline through p1..p2 (p0/p3 shape the ends).
///        Tension 0.5. Returns the interpolated value at fraction @p t in [0,1].
[[nodiscard]] double catmull_rom(double p0, double p1, double p2, double p3,
                                 double t) {
    const auto t2 = t * t;
    const auto t3 = t2 * t;
    const auto a0 = 0.5 * (-p0 + 3.0 * p1 - 3.0 * p2 + p3);
    const auto a1 = 0.5 * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3);
    const auto a2 = 0.5 * (-p0 + p2);
    const auto a3 = 0.5 * 2.0 * p1;
    return a0 * t3 + a1 * t2 + a2 * t + a3;
}

/// @brief Bicubic (Catmull-Rom) sub-pixel sample. Matches unit_interpolate.pas:
///        needs a 4x4 neighbourhood; the outer 2-pixel border samples as black.
void bicubic_interpolate(const ImageArray& img, double x, double y,
                         std::array<float, 3>& pixel) {
    const auto width5    = static_cast<int>(img[0][0].size());
    const auto height5   = static_cast<int>(img[0].size());
    const auto nrcolours = static_cast<int>(img.size());

    const auto int_x = static_cast<int>(std::floor(x));
    const auto int_y = static_cast<int>(std::floor(y));
    const auto dx = x - int_x;
    const auto dy = y - int_y;

    if (int_x <= 1 || int_x >= width5 - 2 || int_y <= 1 || int_y >= height5 - 2) {
        pixel = {0.0f, 0.0f, 0.0f};
        return;
    }

    for (int col = 0; col < nrcolours; ++col) {
        std::array<double, 4> row{};
        for (int i = -1; i <= 2; ++i) {
            row[i + 1] = catmull_rom(img[col][int_y + i][int_x - 1],
                                     img[col][int_y + i][int_x],
                                     img[col][int_y + i][int_x + 1],
                                     img[col][int_y + i][int_x + 2],
                                     dx);
        }
        pixel[col] = static_cast<float>(
            catmull_rom(row[0], row[1], row[2], row[3], dy));
    }
}

// The calibration pipeline maintains these as globals.
// Stubbed here; real values come from astap_main / live stacking.
// NOTE: head_ref is the canonical global in core/globals.h, brought into
// scope via `using astap::head_ref` above.
ImageArray  img_dark;
ImageArray  img_flat;
Header      head_dark;
Header      head_flat;
double      dark_norm_value  = 0.0;
int         last_light_jd    = 0;
int         process_as_osc   = 0;

// Defined further down in this TU; used by the library lookup below.
[[nodiscard]] bool load_master_frame(const std::filesystem::path& path,
                                     ImageArray& img_out, Header& head_out,
                                     int& count_field);

///----------------------------------------
/// MARK: Master calibration library (per-frame hot-reload)
///----------------------------------------

// The candidate libraries the per-frame selector chooses from (pre-analysed
// header metadata). Empty ⇒ the single-resident-master path (set_master_dark) is
// used instead.
std::vector<MasterMetadata> dark_library;
std::vector<MasterMetadata> flat_library;
std::string last_dark_loaded;
std::string last_flat_loaded;
MasterMatchOptions master_match;

// Case-insensitive string equality (Pascal AnsiCompareText). Used by the public
// master_flat_issue below (accessible throughout the TU).
[[nodiscard]] bool iequals(std::string_view a, std::string_view b) {
    if (a.size() != b.size()) { return false; }
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::tolower(static_cast<unsigned char>(a[i]))
            != std::tolower(static_cast<unsigned char>(b[i]))) {
            return false;
        }
    }
    return true;
}

void build_library(std::span<const std::filesystem::path> masters,
                   std::vector<MasterMetadata>& library) {
    library.clear();
    for (const auto& path : masters) {
        MasterMetadata cand;
        if (analyse_master(path, cand)) {
            library.push_back(std::move(cand));
        } else {
            memo2_message("Skipping unreadable master: " + path.string());
        }
    }
}

///----------------------------------------
///  @brief Select + load the closest-matching master dark for one light frame.
///  @details Ports load_master_dark (unit_stack.pas:12028): among the candidate
///           library, keep those matching the enabled classify criteria and the
///           mandatory width/height, pick the one with the closest JD, and reload
///           only when it differs from the last (the hot-reload cache). Emits the
///           compatibility warnings. A no-op when no library is installed, so the
///           single-resident-master path (set_master_dark) is preserved.
///----------------------------------------

void load_master_dark(int jd, Header& hd) {
    if (dark_library.empty()) {
        return;  // single-resident-master path (set_master_dark) — nothing to do
    }

    const int         light_exposure    = static_cast<int>(std::lrint(hd.exposure));
    const int         light_temperature = hd.set_temperature;
    const std::string dark_gain = hd.egain.empty() ? hd.gain : hd.egain;

    std::string filen;
    double day_offset = 99999999.0;
    for (const auto& cand : dark_library) {
        if (!master_dark_issue(hd, cand, master_match).empty()) {
            continue;  // incompatible with this light
        }
        const double d = std::abs(cand.jd - jd);
        if (d < day_offset) {
            filen      = cand.path.string();
            day_offset = d;
        }
    }

    if (filen.empty()) {
        memo2_message("Warning, could not find a suitable dark for "
            + std::to_string(light_exposure) + " sec, temperature "
            + std::to_string(light_temperature) + " deg, gain " + dark_gain
            + " and width " + std::to_string(hd.width)
            + ". De-classify temperature/exposure or add correct darks.");
        head_dark.dark_count = 0;
        return;
    }

    if (head_ref.dark_count == 0 || filen != last_dark_loaded) {
        memo2_message("Loading master dark file " + filen);
        if (!load_master_frame(filen, img_dark, head_dark, head_dark.dark_count)) {
            memo2_message("Error loading master dark " + filen);
            head_ref.dark_count = 0;
            return;
        }
        // Compatibility warnings vs the light (even when it was the best match).
        if (std::lrint(head_dark.exposure) != 0
            && std::lrint(hd.exposure) != std::lrint(head_dark.exposure)) {
            memo2_message("Warning: dark exposure time differs from the light.");
            head_dark.issues += "D_exposure|";
        }
        if (head_dark.set_temperature != 999
            && hd.set_temperature != head_dark.set_temperature) {
            memo2_message("Warning: dark sensor temperature differs from the light.");
            head_dark.issues += "D_temperature|";
        }
        if (!head_dark.gain.empty() && hd.gain != head_dark.gain) {
            memo2_message("Warning: dark gain differs from the light.");
            head_dark.issues += "D_gain|";
        }
        if (!head_dark.egain.empty() && hd.egain != head_dark.egain) {
            memo2_message("Warning: dark egain differs from the light.");
            head_dark.issues += "D_egain|";
        }
        last_dark_loaded = filen;
        if (head_dark.dark_count == 0) { head_dark.dark_count = 1; }
    }
}

///----------------------------------------
///  @brief Select + load the closest-matching master flat for one light frame.
///  @details Ports load_master_flat (unit_stack.pas:12132): match on filter name
///           (when classified) + mandatory dimensions, reject a candidate that is
///           itself a calibrated light/dark (D or F in CALSTAT), pick closest JD,
///           reload only on change, and warn if the flat has no flat-dark/bias.
///           A no-op when no library is installed (set_master_flat preserved).
///----------------------------------------

void load_master_flat(int jd, Header& hd) {
    if (flat_library.empty()) {
        return;
    }

    std::string filen;
    double day_offset = 99999999.0;
    for (const auto& cand : flat_library) {
        if (!master_flat_issue(hd, cand, master_match).empty()) {
            continue;  // incompatible with this light
        }
        const double d = std::abs(cand.jd - jd);
        if (d < day_offset) {
            filen      = cand.path.string();
            day_offset = d;
        }
    }

    if (filen.empty()) {
        memo2_message("Warning, could not find a suitable flat for \""
            + hd.filter_name + "\". De-classify flat filter or add correct flats.");
        head_flat.flat_count = 0;
        return;
    }

    if (head_ref.flat_count == 0 || filen != last_flat_loaded) {
        memo2_message("Loading master flat file " + filen);
        if (!load_master_frame(filen, img_flat, head_flat, head_flat.flat_count)) {
            memo2_message("Error loading master flat " + filen);
            head_flat.flat_count = 0;
            return;
        }
        last_flat_loaded = filen;
        if (head_flat.calstat.find('B') == std::string::npos) {
            memo2_message("Warning: flat not calibrated with a flat-dark/bias "
                          "(CALSTAT/BIAS_CNT).");
        }
        if (head_flat.flat_count == 0) { head_flat.flat_count = 1; }
    }
}

///----------------------------------------
///  @brief Parse an integer from a 1-based substring range.
///  @param s Source string.
///  @param start1 1-based start position.
///  @param len Character count.
/// @param[out] ok Set to @c true on success.
/// @return Parsed value, or 0 on failure.
///----------------------------------------

[[nodiscard]] int parse_int_range(std::string_view s, int start1, int len, bool& ok) {
    ok = false;
    // Require only that the field STARTS within the string; substr clamps a field
    // that runs past the end (Pascal `copy` semantics). This accepts the standard
    // 19-char "YYYY-MM-DDThh:mm:ss" FITS date whose trailing field is shorter than
    // the fixed-width read.
    if (start1 < 1 || start1 - 1 >= static_cast<int>(s.size())) {
        return 0;
    }
    try {
        auto v = std::stoi(std::string(s.substr(start1 - 1, len)));
        ok = true;
        return v;
    } catch (...) {
        return 0;
    }
}

///----------------------------------------
///  @brief Parse a double from a 1-based substring range.
///  @param s Source string.
///  @param start1 1-based start position.
///  @param len Character count.
/// @param[out] ok Set to @c true on success.
/// @return Parsed value, or 0.0 on failure.
///----------------------------------------

[[nodiscard]] double parse_double_range(std::string_view s, int start1, int len, bool& ok) {
    ok = false;
    // Clamp a field that runs past the end rather than rejecting it (Pascal `copy`
    // semantics) so a 19-char FITS date with a 2-char seconds field still parses.
    if (start1 < 1 || start1 - 1 >= static_cast<int>(s.size())) {
        return 0.0;
    }
    try {
        auto v = std::stod(std::string(s.substr(start1 - 1, len)));
        ok = true;
        return v;
    } catch (...) {
        return 0.0;
    }
}

///----------------------------------------
/// @brief Format an integer as a zero-padded two-digit string.
/// @param v Value to format.
/// @return Two-character string.
///----------------------------------------

[[nodiscard]] std::string leading_zero(int v) {
    char buf[8];
    std::snprintf(buf, sizeof(buf), "%02d", v);
    return buf;
}

///----------------------------------------
/// @brief Floating-point modulo that returns a non-negative result.
/// @param a Dividend.
/// @param b Divisor.
/// @return Non-negative remainder.
///----------------------------------------

[[nodiscard]] double fnmodulo(double a, double b) noexcept {
    auto r = std::fmod(a, b);
    if (r < 0.0) {
        r += b;
    }
    return r;
}

///----------------------------------------
/// @brief Clamp @p v into [@p lo, @p hi].
///----------------------------------------

template <typename T>
[[nodiscard]] constexpr T clamp3(T lo, T v, T hi) noexcept {
    return std::min(hi, std::max(lo, v));
}

///----------------------------------------
/// @brief Quicksort-backed median for double arrays.
/// @param v Sample buffer (mutated by nth_element).
/// @param count Number of valid elements at the front of @p v.
/// @return Median value, or 0.0 if @p count <= 0.
///----------------------------------------

[[nodiscard]] double s_median(std::vector<double>& v, int count) {
    if (count <= 0) {
        return 0.0;
    }
    auto mid = v.begin() + count / 2;
    std::nth_element(v.begin(), mid, v.begin() + count);
    return *mid;
}
 
}  // namespace

///----------------------------------------
/// MARK: box_blur
///----------------------------------------

void box_blur(int colors, int range, ImageArray& img) {
    const auto col = static_cast<int>(img.size());
    if (col == 0) {
        return;
    }
    const auto h = static_cast<int>(img[0].size());
    const auto w = static_cast<int>(img[0][0].size());
    
    auto minimum = 0;
    auto maximum = 0;
    if (range == 2) {
        minimum = 0;
        maximum = 1;
    } else if (range == 3) {
        minimum = -1;
        maximum = 1;
    } else if (range == 4) {
        minimum = -1;
        maximum = 2;
    } else {
        minimum = -(range / 2);
        maximum = range / 2;
    }
    
    auto img_temp2 = ImageArray(col,
        std::vector<std::vector<float>>(h, std::vector<float>(w, 0.0f)));
        
    for (auto k = 0; k < col; ++k) {
        for (auto fitsY = 0; fitsY < h; ++fitsY) {
            for (auto fitsX = 0; fitsX < w; ++fitsX) {
                auto value = 0.0;
                auto counter = 0;
                for (auto i = minimum; i <= maximum; ++i) {
                    for (auto j = minimum; j <= maximum; ++j) {
                        auto x1 = fitsX + i;
                        auto y1 = fitsY + j;
                        if (x1 >= 0 && x1 <= w - 1 && y1 >= 0 && y1 <= h - 1) {
                            auto value2 = img[k][y1][x1];
                            if (value2 != 0.0f) {
                                value += value2;
                                ++counter;
                            }
                        }
                    }
                }
                img_temp2[k][fitsY][fitsX] = counter != 0
                    ? static_cast<float>(value / counter)
                    : 0.0f;
            }
        }
    }
    
    if (colors == 1 && col == 3) {
        // Collapse to mono. The original loop writes the averaged value
        // inside an inner `for k` loop but only reads channel 0 of
        // img_temp2. Preserved verbatim to keep output identical.
        for (auto fitsY = 0; fitsY < h; ++fitsY) {
            for (auto fitsX = 0; fitsX < w; ++fitsX) {
                for (auto k = 0; k < col; ++k) {
                    img[0][fitsY][fitsX] =
                        (img_temp2[0][fitsY][fitsX] + img_temp2[1][fitsY][fitsX]
                         + img_temp2[2][fitsY][fitsX]) / 3.0f;
                }
            }
        }
    } else {
        img = std::move(img_temp2);
    }
    
    // TODO: the original also sets head.naxis3 := colors.
    // The global `head` is not plumbed into this module yet; the
    // caller should update naxis3 after calling box_blur.
}

///----------------------------------------
/// MARK: check_pattern_filter
///----------------------------------------

void check_pattern_filter(ImageArray& img) {
    const auto col = static_cast<int>(img.size());
    if (col == 0) {
        return;
    }
    const auto h = static_cast<int>(img[0].size());
    const auto w = static_cast<int>(img[0][0].size());
    
    if (col > 1) {
        memo2_message("Skipping normalise filter. This filter works only for raw OSC images!");
        return;
    }
    memo2_message("Normalise raw OSC image by applying check pattern filter.");
    
    auto value1 = 0.0;
    auto value2 = 0.0;
    auto value3 = 0.0;
    auto value4 = 0.0;
    auto counter1 = 0;
    auto counter2 = 0;
    auto counter3 = 0;
    auto counter4 = 0;
    
    // Accumulate means for each Bayer cell within the central 50%
    for (auto fitsY = h / 4; fitsY <= (h * 3) / 4; ++fitsY) {
        for (auto fitsX = w / 4; fitsX <= (w * 3) / 4; ++fitsX) {
            auto oddX = (fitsX & 1) != 0;
            auto oddY = (fitsY & 1) != 0;
            if (!oddX && !oddY) {
                value1 += img[0][fitsY][fitsX];
                ++counter1;
            } else if (oddX && !oddY) {
                value2 += img[0][fitsY][fitsX];
                ++counter2;
            } else if (!oddX && oddY) {
                value3 += img[0][fitsY][fitsX];
                ++counter3;
            } else {
                value4 += img[0][fitsY][fitsX];
                ++counter4;
            }
        }
    }
    
    if (counter1 == 0 || counter2 == 0 || counter3 == 0 || counter4 == 0) {
        return;
    }
    
    // Compute normalisation factors
    value1 /= counter1;
    value2 /= counter2;
    value3 /= counter3;
    value4 /= counter4;
    auto maxval = std::max(std::max(value1, value2), std::max(value3, value4));
    value1 = maxval / value1;
    value2 = maxval / value2;
    value3 = maxval / value3;
    value4 = maxval / value4;
    
    // Apply normalisation to the full image
    for (auto fitsY = 0; fitsY < h; ++fitsY) {
        for (auto fitsX = 0; fitsX < w; ++fitsX) {
            auto oddX = (fitsX & 1) != 0;
            auto oddY = (fitsY & 1) != 0;
            if (value1 != 1.0 && !oddX && !oddY) {
                img[0][fitsY][fitsX] =
                    static_cast<float>(std::round(img[0][fitsY][fitsX] * value1));
            } else if (value2 != 1.0 && oddX && !oddY) {
                img[0][fitsY][fitsX] =
                    static_cast<float>(std::round(img[0][fitsY][fitsX] * value2));
            } else if (value3 != 1.0 && !oddX && oddY) {
                img[0][fitsY][fitsX] =
                    static_cast<float>(std::round(img[0][fitsY][fitsX] * value3));
            } else if (value4 != 1.0 && oddX && oddY) {
                img[0][fitsY][fitsX] =
                    static_cast<float>(std::round(img[0][fitsY][fitsX] * value4));
            }
        }
    }
}

///----------------------------------------
/// MARK: black_spot_filter
///----------------------------------------

void black_spot_filter(ImageArray& img) {
    const auto col = static_cast<int>(img.size());
    if (col == 0) {
        return;
    }
    const auto h = static_cast<int>(img[0].size());
    const auto w = static_cast<int>(img[0][0].size());
    
    auto img_temp2 = ImageArray(col,
        std::vector<std::vector<float>>(h, std::vector<float>(w, 0.0f)));
        
    for (auto k = 0; k < col; ++k) {
        // Detect non-black borders per channel
        auto left = -1;
        do {
            ++left;
        } while (!(img[k][h / 2][left] != 0.0f || left >= w - 1));
        
        auto right = w;
        do {
            --right;
        } while (!(img[k][h / 2][right] != 0.0f || right <= 0));
        
        auto bottom = -1;
        do {
            ++bottom;
        } while (!(img[k][bottom][w / 2] != 0.0f || bottom >= h - 1));
        
        auto top = h;
        do {
            --top;
        } while (!(img[k][top][w / 2] != 0.0f || top <= 0));
        
        // Replace interior zero-pixels with the mean of non-zero neighbours
        for (auto fitsY = 0; fitsY < h; ++fitsY) {
            for (auto fitsX = 0; fitsX < w; ++fitsX) {
                auto value = static_cast<double>(img[k][fitsY][fitsX]);
                if (value == 0.0
                    && fitsX >= left && fitsX <= right
                    && fitsY >= bottom && fitsY <= top) {
                    auto range = 1;
                    auto counter = 0;
                    auto acc = 0.0;
                    do {
                        counter = 0;
                        acc = value;
                        for (auto i = -range; i <= range; ++i) {
                            for (auto j = -range; j <= range; ++j) {
                                if (std::abs(i) == range || std::abs(j) == range) {
                                    auto x1 = fitsX + i;
                                    auto y1 = fitsY + j;
                                    if (x1 >= left && x1 <= right
                                        && y1 >= bottom && y1 <= top) {
                                        auto v2 = img[k][y1][x1];
                                        if (v2 != 0.0f) {
                                            acc += v2;
                                            ++counter;
                                        }
                                    }
                                }
                            }
                        }
                        if (counter != 0) {
                            value = acc / counter;
                        } else {
                            ++range;
                        }
                    } while (!(counter != 0 || range >= 10));
                }
                img_temp2[k][fitsY][fitsX] = static_cast<float>(value);
            }
        }
    }
    
    img = std::move(img_temp2);
}

///----------------------------------------
/// MARK: update_solution_and_save
///----------------------------------------

[[nodiscard]] bool update_solution_and_save(const ImageArray& img,
                                            Header& hd,
                                            std::vector<std::string>& memo,
                                            const std::filesystem::path& filename) {
    // Plate-solve the already-loaded interim frame with the internal solver,
    // then persist the WCS back into the frame's own header so the alignment
    // vector — and any later pass that reloads the file — can read it.
    // solve_image bins/crops in place, so hand it a copy of the caller's image.
    auto img_copy = img;
    if (!astap::solving::solve_image(img_copy, hd, memo, /*get_hist=*/true,
                                     /*check_patternfilter=*/false)) {
        return false;
    }

    // FITS keeps its pixels and rewrites the header in place; anything else
    // (TIFF, PNG, JPG) is written out as a 16-bit TIFF carrying the header.
    const bool ok = astap::core::fits_file_name(filename.string())
        ? astap::core::savefits_update_header(memo, filename)
        : astap::core::save_tiff16_secure(img_copy, memo, filename);
    if (!ok) {
        memo2_message("Write error !!" + filename.string());
    }
    return ok;
}

///----------------------------------------
/// MARK: apply_dark_and_flat
///----------------------------------------

[[nodiscard]] bool apply_dark_and_flat(ImageArray& img, Header& hd) {
    auto result = false;
    date_to_jd(hd.date_obs, hd.date_avg, hd.exposure);
    
    // Dark subtraction
    if (hd.calstat.find('D') != std::string::npos) {
        memo2_message("Skipping dark calibration, already applied. See header keyword CALSTAT");
    } else {
        load_master_dark(static_cast<int>(std::round(jd_start)), hd);
        if (head_dark.dark_count > 0 && !img_dark.empty() && !img_dark[0].empty()) {
            dark_norm_value = 0.0;
            for (auto fitsY = (hd.width / 2) - 7; fitsY <= (hd.width / 2) + 8; ++fitsY) {
                for (auto fitsX = (hd.height / 2) - 7; fitsX <= (hd.height / 2) + 8; ++fitsX) {
                    dark_norm_value += img_dark[0][fitsY][fitsX];
                }
            }
            dark_norm_value /= (16.0 * 16.0);
            
            // Subtract dark from all channels
            for (auto fitsY = 0; fitsY < hd.height; ++fitsY) {
                for (auto fitsX = 0; fitsX < hd.width; ++fitsX) {
                    auto value = img_dark[0][fitsY][fitsX];
                    for (auto k = 0; k < hd.naxis3; ++k) {
                        img[k][fitsY][fitsX] -= value;
                    }
                }
            }
            
            head_ref.calstat      = "D";
            head_ref.dark_count   = head_dark.dark_count;
            head_ref.datamax_org  = hd.datamax_org - dark_norm_value;
            
            hd.calstat      = "D";
            hd.dark_count   = head_dark.dark_count;
            hd.datamax_org  = hd.datamax_org - dark_norm_value;
            result = true;
        }
    }
    
    // Flat division
    if (hd.calstat.find('F') != std::string::npos) {
        memo2_message("Skipping flat calibration, already applied. See header keyword CALSTAT");
    } else {
        load_master_flat(static_cast<int>(std::round(jd_start)), hd);
        last_light_jd = static_cast<int>(std::round(jd_start));
        
        if (head_flat.flat_count != 0 && !img_flat.empty() && !img_flat[0].empty()) {
            auto flat_norm_value = 0.0;
            auto flatNorm11 = 0.0;
            auto flatNorm12 = 0.0;
            auto flatNorm21 = 0.0;
            auto flatNorm22 = 0.0;
            
            // Compute flat normalisation values from the centre 16x16
            for (auto fitsY = (head_flat.width / 2) - 7; fitsY <= (head_flat.width / 2) + 8; ++fitsY) {
                for (auto fitsX = (head_flat.height / 2) - 7; fitsX <= (head_flat.height / 2) + 8; ++fitsX) {
                    auto value = static_cast<double>(img_flat[0][fitsY][fitsX]);
                    flat_norm_value += value;
                    auto oddX = (fitsX & 1) != 0;
                    auto oddY = (fitsY & 1) != 0;
                    if (oddX) {
                        if (oddY) {
                            flatNorm11 += value;
                        } else {
                            flatNorm12 += value;
                        }
                    } else {
                        if (oddY) {
                            flatNorm21 += value;
                        } else {
                            flatNorm22 += value;
                        }
                    }
                }
            }
            flat_norm_value /= (16.0 * 16.0);
            
            if (process_as_osc > 0) {
                // OSC flat: per-Bayer-cell normalisation
                auto hi = std::max(std::max(flatNorm11, flatNorm12),
                                   std::max(flatNorm21, flatNorm22));
                auto lo = std::min(std::min(flatNorm11, flatNorm12),
                                   std::min(flatNorm21, flatNorm22));
                if (lo > 0.0 && hi / lo > 2.0) {
                    memo2_message("Warning flat pixel colour values differ too much. "
                                  "Use white light for OSC flats!! Will compensate accordingly.");
                }
                
                flatNorm11 /= (8.0 * 8.0);
                flatNorm12 /= (8.0 * 8.0);
                flatNorm21 /= (8.0 * 8.0);
                flatNorm22 /= (8.0 * 8.0);
                
                for (auto fitsY = 0; fitsY < hd.height; ++fitsY) {
                    for (auto fitsX = 0; fitsX < hd.width; ++fitsX) {
                        auto oddX = (fitsX & 1) != 0;
                        auto oddY = (fitsY & 1) != 0;
                        auto denom = static_cast<double>(img_flat[0][fitsY][fitsX]) + 0.001;
                        auto flat_factor = 0.0;
                        if (oddX) {
                            flat_factor = oddY ? flatNorm11 / denom : flatNorm12 / denom;
                        } else {
                            flat_factor = oddY ? flatNorm21 / denom : flatNorm22 / denom;
                        }
                        flat_factor = clamp3<double>(-4.0, flat_factor, 4.0);
                        img[0][fitsY][fitsX] =
                            static_cast<float>(img[0][fitsY][fitsX] * flat_factor);
                    }
                }
            } else {
                // Mono flat: uniform normalisation
                for (auto k = 0; k < hd.naxis3; ++k) {
                    for (auto fitsY = 0; fitsY < hd.height; ++fitsY) {
                        for (auto fitsX = 0; fitsX < hd.width; ++fitsX) {
                            auto flat_factor = flat_norm_value
                                / (img_flat[0][fitsY][fitsX] + 0.001);
                            flat_factor = clamp3<double>(-4.0, flat_factor, 4.0);
                            img[k][fitsY][fitsX] =
                                static_cast<float>(img[k][fitsY][fitsX] * flat_factor);
                        }
                    }
                }
            }
            
            head_ref.calstat        += "F" + head_flat.calstat;
            head_ref.flat_count      = head_flat.flat_count;
            head_ref.flatdark_count  = head_flat.flatdark_count;
            head_ref.issues          = head_dark.issues + head_flat.issues;
            
            hd.calstat        += "F" + head_flat.calstat;
            hd.flat_count      = head_flat.flat_count;
            hd.flatdark_count  = head_flat.flatdark_count;
            result = true;
        }
    }
    return result;
}

///----------------------------------------
/// MARK: Master calibration setters
///----------------------------------------

namespace {

bool load_master_frame(const std::filesystem::path& path,
                       ImageArray& img_out,
                       Header& head_out,
                       int& count_field) {
    img_out.clear();
    head_out = Header{};
    count_field = 0;
    if (path.empty()) {
        return true;
    }
    auto memo = std::vector<std::string>{};
    if (!astap::core::load_fits(path, /*light=*/false, /*load_data=*/true,
                                /*update_memo=*/true, /*get_ext=*/0,
                                memo, head_out, img_out)) {
        img_out.clear();
        head_out = Header{};
        return false;
    }
    // The CALSTAT / NLIGHT counts from the master's header already live in
    // head_out; if missing, treat this as a single-frame master so the
    // apply_dark_and_flat predicate (count > 0) still fires.
    if (count_field == 0) {
        count_field = 1;
    }
    return true;
}

}  // namespace

namespace {

void fill_info(const Header& h, MasterFrameInfo& out) {
    out.width    = h.width;
    out.height   = h.height;
    out.naxis3   = h.naxis3;
    out.exposure = h.exposure;
    out.calstat  = h.calstat;
}

}  // namespace

bool set_master_dark(const std::filesystem::path& path, MasterFrameInfo& info) {
    if (!load_master_frame(path, img_dark, head_dark, head_dark.dark_count)) {
        return false;
    }
    fill_info(head_dark, info);
    return true;
}

bool set_master_flat(const std::filesystem::path& path, MasterFrameInfo& info) {
    if (!load_master_frame(path, img_flat, head_flat, head_flat.flat_count)) {
        return false;
    }
    fill_info(head_flat, info);
    return true;
}

void set_master_match_options(const MasterMatchOptions& options) {
    master_match = options;
}

void set_dark_library(std::span<const std::filesystem::path> masters) {
    build_library(masters, dark_library);
    last_dark_loaded.clear();
    if (dark_library.empty()) {
        // No candidates ⇒ no dark. Also drop the resident master so a prior run's
        // dark can't leak into a later library-less run in the same process (the
        // per-frame load_master_dark early-returns when the library is empty and
        // would otherwise leave stale img_dark/head_dark for apply_dark_and_flat
        // to reuse). Pascal resets dark_count to 0 when nothing matches
        // (unit_stack.pas:12126). The single-resident-master path (set_master_dark,
        // used by the CLI batch driver) never calls this, so it is unaffected.
        img_dark.clear();
        head_dark = Header{};
    }
}

void set_flat_library(std::span<const std::filesystem::path> masters) {
    build_library(masters, flat_library);
    last_flat_loaded.clear();
    if (flat_library.empty()) {
        img_flat.clear();
        head_flat = Header{};
    }
}

bool analyse_master(const std::filesystem::path& path, MasterMetadata& out) {
    auto memo = std::vector<std::string>{};
    auto img  = ImageArray{};
    auto head = Header{};
    if (!astap::core::load_fits(path, /*light=*/false, /*load_data=*/false,
                                /*update_memo=*/true, /*get_ext=*/0,
                                memo, head, img)) {
        return false;
    }
    out.path            = path;
    out.exposure        = static_cast<int>(std::lrint(head.exposure));
    out.exposure_seconds = head.exposure;
    out.set_temperature = head.set_temperature;
    out.gain            = head.egain.empty() ? head.gain : head.egain;
    out.xbinning        = head.xbinning;
    out.width           = head.width;
    out.height          = head.height;
    out.date_obs        = head.date_obs;
    out.filter_name     = head.filter_name;
    out.calstat         = head.calstat;
    // JD from DATE-OBS (date_to_jd writes the jd_start global).
    date_to_jd(head.date_obs, head.date_avg, head.exposure);
    out.jd = astap::jd_start;
    return true;
}

std::string master_dark_issue(const Header& light, const MasterMetadata& candidate,
                              const MasterMatchOptions& options) {
    // Same gate order + operators as the Pascal (unit_stack.pas:12049-12064); the
    // returned string is the reason ASTAP writes to its Compatibility column
    // (:12068-12080), showing the light's required value. Empty ⇒ eligible.
    const int light_exposure = static_cast<int>(std::lrint(light.exposure));
    if (options.classify_exposure && light_exposure != candidate.exposure) {
        return "exposure<>" + std::to_string(light_exposure);
    }
    if (options.classify_temperature
        && std::abs(light.set_temperature - candidate.set_temperature) > options.delta_temp) {
        return "temperature<>" + std::to_string(light.set_temperature);
    }
    const std::string light_gain = light.egain.empty() ? light.gain : light.egain;
    if (options.classify_gain && light_gain != candidate.gain) {
        return "gain<>" + light_gain;
    }
    if (light.width != candidate.width)   { return "width<>"  + std::to_string(light.width); }
    if (light.height != candidate.height) { return "height<>" + std::to_string(light.height); }
    return "";
}

std::string master_flat_issue(const Header& light, const MasterMetadata& candidate,
                              const MasterMatchOptions& options) {
    // Pascal order (unit_stack.pas:12145-12163): filter (non-fatal in the GUI, but
    // still excludes the candidate) → width → height → not-a-flat. The port surfaces
    // the filter mismatch as a reason string too (more informative; selection is
    // unchanged). A D/F in CALSTAT marks a calibrated light/dark, not a flat; the
    // Pascal `pos('F',imagetype)` override is omitted (port Header has no IMAGETYP).
    if (options.classify_filter && !iequals(light.filter_name, candidate.filter_name)) {
        return "filter<>" + light.filter_name;
    }
    if (light.width != candidate.width)   { return "width<>"  + std::to_string(light.width); }
    if (light.height != candidate.height) { return "height<>" + std::to_string(light.height); }
    if (candidate.calstat.find('D') != std::string::npos
        || candidate.calstat.find('F') != std::string::npos) {
        return "Calibration, not a flat!";
    }
    return "";
}

///----------------------------------------
/// MARK: smart_colour_smooth
///----------------------------------------

void smart_colour_smooth(ImageArray& img,
                         double wide,
                         double sd,
                         bool preserve_r_nebula,
                         bool measurehist) {
    if (img.size() < 3) {
        return;
    }
    
    const auto width5  = static_cast<int>(img[0][0].size());
    const auto height5 = static_cast<int>(img[0].size());
    
    auto img_temp2 = ImageArray(3,
        std::vector<std::vector<float>>(height5, std::vector<float>(width5, 0.0f)));
        
    auto step = static_cast<int>(std::round(wide)) / 2;
    
    // Measure background and noise per channel
    auto bckR = Background{};
    auto bckG = Background{};
    auto bckB = Background{};
    get_background(0, img, measurehist, true, bckR);
    auto bgR = bckR.backgr;
    get_background(1, img, measurehist, true, bckG);
    auto bgG = bckG.backgr;
    get_background(2, img, measurehist, true, bckB);
    auto bgB = bckB.backgr;
    
    auto star_level = 30.0 * std::max(bckR.noise_level,
                                      std::max(bckG.noise_level, bckB.noise_level));
    auto bg = (bgR + bgG + bgB) / 3.0;
    
    for (auto fitsY = 0; fitsY < height5; ++fitsY) {
        for (auto fitsX = 0; fitsX < width5; ++fitsX) {
            auto red = 0.0;
            auto green = 0.0;
            auto blue = 0.0;
            auto count = 0;
            auto peak = 0.0;
            auto bgR2 = 65535.0;
            auto bgG2 = 65535.0;
            auto bgB2 = 65535.0;
            
            auto r2 = static_cast<float>(img[0][fitsY][fitsX] - bgR);
            auto g2 = static_cast<float>(img[1][fitsY][fitsX] - bgG);
            auto b2 = static_cast<float>(img[2][fitsY][fitsX] - bgB);
            
            // Only process pixels above the noise threshold
            if (r2 > sd * bckR.noise_level
                || g2 > sd * bckG.noise_level
                || b2 > sd * bckB.noise_level) {
                for (auto y = -step; y <= step; ++y) {
                    for (auto x = -step; x <= step; ++x) {
                        auto x2 = fitsX + x;
                        auto y2 = fitsY + y;
                        if (x2 >= 0 && x2 < width5 && y2 >= 0 && y2 < height5) {
                            auto sqr_dist = x * x + y * y;
                            if (sqr_dist <= step * step) {
                                auto r = static_cast<double>(img[0][y2][x2]);
                                auto g = static_cast<double>(img[1][y2][x2]);
                                auto b = static_cast<double>(img[2][y2][x2]);
                                if (r > peak) {
                                    peak = r;
                                }
                                if (g > peak) {
                                    peak = g;
                                }
                                if (b > peak) {
                                    peak = b;
                                }
                                if (r < bgR2) {
                                    bgR2 = r;
                                }
                                if (g < bgG2) {
                                    bgG2 = g;
                                }
                                if (b < bgB2) {
                                    bgB2 = b;
                                }
                                if (r < 60000 && g < 60000 && b < 60000) {
                                    if (r - bgR > 0) {
                                        red   += r - bgR;
                                    }
                                    if (g - bgG > 0) {
                                        green += g - bgG;
                                    }
                                    if (b - bgB > 0) {
                                        blue  += b - bgB;
                                    }
                                    ++count;
                                }
                            }
                        }
                    }
                }
            }
            
            auto copydata = true;
            if (count >= 1) {
                red   /= count;
                green /= count;
                blue  /= count;
                
                if (peak > star_level) {
                    auto highest_colour = std::max({static_cast<double>(r2),
                                                    static_cast<double>(g2),
                                                    static_cast<double>(b2)});
                    auto red_nebula = false;
                    if (preserve_r_nebula) {
                        red_nebula = (highest_colour == r2)
                            && ((r2 - (bgR2 - bgR)) < 150.0)
                            && ((bgR2 - bgR) > 3.0 * bckR.noise_level);
                    }
                    if (!red_nebula) {
                        // Correct green channel relative to red/blue
                        if (red < blue * 1.06) {
                            green = std::max(green, 0.6604 * red + 0.3215 * blue);
                        }
                        auto flux = static_cast<double>(r2) + g2 + b2;
                        auto rgb  = red + green + blue + 0.00001;
                        auto strongest_colour_local =
                            std::max({red, green, blue});
                        auto top = bg + strongest_colour_local * (flux / rgb);
                        if (top >= 65534.99) {
                            flux -= (top - 65534.99) * rgb / strongest_colour_local;
                        }
                        auto lumr = flux / rgb;
                        img_temp2[0][fitsY][fitsX] = static_cast<float>(bg + red   * lumr);
                        img_temp2[1][fitsY][fitsX] = static_cast<float>(bg + green * lumr);
                        img_temp2[2][fitsY][fitsX] = static_cast<float>(bg + blue  * lumr);
                        copydata = false;
                    }
                }
            }
            if (copydata) {
                img_temp2[0][fitsY][fitsX] = static_cast<float>(std::max(0.0, bg + r2));
                img_temp2[1][fitsY][fitsX] = static_cast<float>(std::max(0.0, bg + g2));
                img_temp2[2][fitsY][fitsX] = static_cast<float>(std::max(0.0, bg + b2));
            }
        }
    }
    
    img = std::move(img_temp2);
}

///----------------------------------------
/// MARK: green_purple_filter
///----------------------------------------

void green_purple_filter(ImageArray& img) {
    if (img.size() < 3) {
        return;
    }
    const auto h = static_cast<int>(img[0].size());
    const auto w = static_cast<int>(img[0][0].size());
    
    for (auto fitsY = 0; fitsY < h; ++fitsY) {
        for (auto fitsX = 0; fitsX < w; ++fitsX) {
            auto r2 = static_cast<double>(img[0][fitsY][fitsX]);
            auto g2 = static_cast<double>(img[1][fitsY][fitsX]);
            auto b2 = static_cast<double>(img[2][fitsY][fitsX]);
            
            // Green dominant: push green down to match the weaker of R/B
            if (g2 > r2 && g2 > b2) {
                auto lum = r2 + g2 + b2;
                auto ratio = 0.0;
                if (r2 >= b2) {
                    ratio = std::min(r2 / std::max(b2, 0.001), 30.0);
                    r2 = lum * ratio / (ratio + ratio + 1);
                    g2 = r2;
                    b2 = lum * 1.0 / (ratio + ratio + 1);
                } else {
                    ratio = std::min(b2 / std::max(r2, 0.001), 30.0);
                    b2 = lum * ratio / (ratio + ratio + 1);
                    g2 = b2;
                    r2 = lum * 1.0 / (ratio + ratio + 1);
                }
                img[0][fitsY][fitsX] = static_cast<float>(r2);
                img[1][fitsY][fitsX] = static_cast<float>(g2);
                img[2][fitsY][fitsX] = static_cast<float>(b2);
            }
            
            // Green deficit: push green up to match the weaker of R/B
            if (g2 < r2 && g2 < b2) {
                auto lum = r2 + g2 + b2;
                auto ratio = 0.0;
                if (r2 >= b2) {
                    ratio = std::min(r2 / std::max(b2, 0.001), 30.0);
                    r2 = lum / (1 + 1 + ratio);
                    g2 = r2;
                    b2 = lum * ratio / (1 + 1 + ratio);
                } else {
                    ratio = std::min(b2 / std::max(r2, 0.001), 30.0);
                    b2 = lum / (1 + 1 + ratio);
                    g2 = b2;
                    r2 = lum * ratio / (1 + 1 + ratio);
                }
                img[0][fitsY][fitsX] = static_cast<float>(r2);
                img[1][fitsY][fitsX] = static_cast<float>(g2);
                img[2][fitsY][fitsX] = static_cast<float>(b2);
            }
        }
    }
}

///----------------------------------------
/// MARK: julian_calc
///----------------------------------------

[[nodiscard]] double julian_calc(int yyyy, int mm, double dd,
                                 double hours, double minutes, double seconds) {
    auto Y = 0;
    auto M = 0;
    if (mm > 2) {
        Y = yyyy;
        M = mm;
    } else {
        Y = yyyy - 1;
        M = mm + 12;
    }
    
    auto DD = dd + hours / 24.0 + minutes / (24.0 * 60.0)
            + seconds / (24.0 * 60.0 * 60.0);
            
    auto B = 0.0;
    if ((yyyy + mm / 100.0 + DD / 10000.0) < 1582.10149999) {
        B = 0.0;
    } else {
        auto A = std::trunc(Y / 100.0);
        B = 2.0 - A + std::trunc(A / 4.0);
    }
    
    auto xx = (Y < 0) ? 0.75 : 0.0;
    return std::trunc(365.25 * Y - xx) + std::trunc(30.6001 * (M + 1))
         + DD + B + 1720994.5;
}

///----------------------------------------
/// MARK: date_to_jd
///----------------------------------------

void date_to_jd(std::string_view date_obs, std::string_view date_avg, double exp) {
    jd_start = 0;
    jd_end   = 0;
    jd_mid   = 0;
    
    // Parse a FITS date string into components
    auto parse = [](std::string_view s,
                    int& yy, int& mm, int& dd, int& hh, int& mn, double& ss) -> bool {
        auto ok = false;
        ss = parse_double_range(s, 18, 7, ok);
        if (!ok) {
            return false;
        }
        mn = parse_int_range(s, 15, 2, ok);
        if (!ok) {
            return false;
        }
        hh = parse_int_range(s, 12, 2, ok);
        if (!ok) {
            return false;
        }
        dd = parse_int_range(s, 9, 2, ok);
        if (!ok) {
            return false;
        }
        mm = parse_int_range(s, 6, 2, ok);
        if (!ok) {
            return false;
        }
        yy = parse_int_range(s, 1, 4, ok);
        if (!ok) {
            return false;
        }
        return true;
    };
    
    auto yy = 0;
    auto mm = 0;
    auto dd = 0;
    auto hh = 0;
    auto mn = 0;
    auto ss = 0.0;
    if (!date_obs.empty()) {
        if (!parse(date_obs, yy, mm, dd, hh, mn, ss)) {
            return;
        }
        jd_start = julian_calc(yy, mm, dd, hh, mn, ss);
        jd_end   = jd_start + exp / (24.0 * 3600.0);
    }
    
    if (date_avg.empty()) {
        jd_mid = jd_start + exp / (2.0 * 24.0 * 3600.0);
    } else {
        if (!parse(date_avg, yy, mm, dd, hh, mn, ss)) {
            return;
        }
        jd_mid = julian_calc(yy, mm, dd, hh, mn, ss);
    }
}

///----------------------------------------
/// MARK: jd_to_date
///----------------------------------------

[[nodiscard]] std::string jd_to_date(double jd) {
    if (std::abs(jd) > 1461.0 * 10000.0) {
        return "Error, JD outside allowed range!";
    }
    jd = jd + (0.5 / (24.0 * 3600.0));
    
    auto Z = std::trunc(jd + 0.5);
    auto F = (jd + 0.5) - std::trunc(jd + 0.5);
    auto A = 0.0;
    if (Z < 2299160.5) {
        A = Z;
    } else {
        auto g = std::trunc((Z - 1867216.25) / 36524.25);
        A = Z + 1 + g - std::trunc(g / 4.0);
    }
    auto B = A + 1524.0 + (1461.0 * 10000.0);
    auto C = std::trunc((B - 122.1) / 365.25);
    auto D = std::trunc(365.25 * C);
    auto E = std::trunc((B - D) / 30.6001);
    auto T = B - D - std::trunc(30.6001 * E) + F;
    auto M = (E < 14) ? E - 1 : E - 13;
    auto J = (M > 2) ? C - 4716 : C - 4715;
    J = J - 4 * 10000;
    
    F = fnmodulo(F, 1.0);
    auto HH = static_cast<int>(std::trunc(F * 24.0));
    auto MM = static_cast<int>(std::trunc((F - HH / 24.0) * (24.0 * 60.0)));
    auto SS = static_cast<int>(std::trunc((F - HH / 24.0 - MM / (24.0 * 60.0))
                                          * (24.0 * 3600.0)));
                                          
    char year4[16];
    std::snprintf(year4, sizeof(year4), "%4d", static_cast<int>(std::trunc(J)));
    
    return std::string(year4) + "-"
         + leading_zero(static_cast<int>(std::trunc(M))) + "-"
         + leading_zero(static_cast<int>(std::trunc(T))) + "T"
         + leading_zero(HH) + ":"
         + leading_zero(MM) + ":"
         + leading_zero(SS);
}

///----------------------------------------
/// MARK: resize_img_loaded
///----------------------------------------

void resize_img_loaded(double ratio) {
    if (ratio <= 0.0 || std::abs(ratio - 1.0) < 1e-6 || img_loaded.empty()) {
        return;
    }
    const auto w1 = head.width;
    const auto h1 = head.height;
    const auto w2 = static_cast<int>(std::lround(ratio * w1));
    const auto h2 = static_cast<int>(std::lround(ratio * h1));
    if (w2 < 5) {
        return;  // matches the original's guard against degenerate targets
    }

    // Bicubic (Catmull-Rom) resample into a fresh buffer, then swap in. Matches
    // the original resize_img_loaded (unit_stack.pas), which samples img_loaded
    // at (fitsX/ratio, fitsY/ratio); the outer border samples as black.
    auto out = ImageArray{};
    out.assign(head.naxis3,
        std::vector<std::vector<float>>(h2, std::vector<float>(w2, 0.0f)));

    std::array<float, 3> colour{};
    for (auto fitsY = 0; fitsY < h2; ++fitsY) {
        const auto y = fitsY / ratio;
        for (auto fitsX = 0; fitsX < w2; ++fitsX) {
            const auto x = fitsX / ratio;
            bicubic_interpolate(img_loaded, x, y, colour);
            for (auto col = 0; col < head.naxis3; ++col) {
                out[col][fitsY][fitsX] = colour[col];
            }
        }
    }
    img_loaded = std::move(out);

    // Update WCS so the scaled frame still points at the right sky.
    // Pixel count scales by ratio; scale per pixel scales by 1/ratio;
    // reference pixel index scales by ratio.
    head.width  = w2;
    head.height = h2;
    head.cdelt1 /= ratio;
    head.cdelt2 /= ratio;
    head.crpix1 *= ratio;
    head.crpix2 *= ratio;
    head.cd1_1  /= ratio;
    head.cd1_2  /= ratio;
    head.cd2_1  /= ratio;
    head.cd2_2  /= ratio;
}

///----------------------------------------
/// MARK: median_background
///----------------------------------------

[[nodiscard]] double median_background(const ImageArray& img,
                                       int color, int sizeX, int sizeY, int x, int y) {
    // Round up to odd dimensions
    if ((sizeX / 2) * 2 == sizeX) {
        sizeX += 1;
    }
    if ((sizeY / 2) * 2 == sizeY) {
        sizeY += 1;
    }
    auto size2 = sizeX * sizeY;
    auto pixArray = std::vector<double>(size2);
    auto stepX = sizeX / 2;
    auto stepY = sizeY / 2;
    auto count = 0;
    const auto w = img.empty() || img[0].empty() ? 0
                 : static_cast<int>(img[0][0].size());
    const auto h = img.empty() ? 0 : static_cast<int>(img[0].size());
    
    for (auto j = y - stepY; j <= y + stepY; ++j) {
        for (auto i = x - stepX; i <= x + stepX; ++i) {
            if (i >= 0 && i < w && j >= 0 && j < h) {
                auto value = static_cast<double>(img[color][j][i]);
                if (value != 0.0) {
                    pixArray[count++] = value;
                }
            }
        }
    }
    
    if (count == 0) {
        return 0.0;
    }
    return s_median(pixArray, count);
}

///----------------------------------------
/// MARK: analyse_image
///----------------------------------------

void analyse_image(const ImageArray& img,
                   const Header& head,
                   double snr_min,
                   int report_type,
                   int& star_counter,
                   Background& bck,
                   double& hfd_median,
                   std::string* csv_out) {
    star_counter = 0;
    hfd_median   = 99.0;
    
    if (img.empty() || img[0].empty()) {
        return;
    }
    
    const auto width5  = static_cast<int>(img[0][0].size());
    const auto height5 = static_cast<int>(img[0].size());
    
    // max_stars mirrors the original's stackmenu1.max_stars1 field, which the
    // CLI -s option populates (astap::max_stars_setting, default 500). It caps
    // the retry loop that lowers the detection level until enough stars appear.
    const auto max_stars = astap::max_stars_setting;
    // hfd_min matches the original's max(0.8, min_star_size_stacking): the
    // stacking field is not command-line settable, so it stays at the
    // two-pixel floor of 0.8 that rejects hot pixels.
    const auto hfd_min = 0.8;
    
    auto len = 1000;
    auto hfd_list = std::vector<double>(len);
    
    get_background(0, img, true, true, bck);
    auto detection_level = bck.star_level;
    
    // TODO: ap_order is a global in unit_astrometric_solving; assume 0.
    constexpr auto ap_order = 0;
    constexpr auto formalism = ap_order > 0 ? 1 : 0;
    
    auto retries = 3;
    
    // Low-dynamic-range images (8-bit, or a datamax at/below 255) get a
    // min_background of 0; everything else 8.
    const auto min_background =
        (astap::nrbits == 8 || head.datamax_org <= 255) ? 0.0 : 8.0;
    
    auto startext = std::string{};
    
    if (bck.backgr < 60000 && bck.backgr > min_background) {
        auto img_sa = ImageArray(1,
            std::vector<std::vector<float>>(height5,
                std::vector<float>(width5, -1.0f)));
                
        do {
            if (retries == 3) {
                if (bck.star_level > 30 * bck.noise_level) {
                    detection_level = bck.star_level;
                } else {
                    retries = 2;
                }
            }
            if (retries == 2) {
                if (bck.star_level2 > 30 * bck.noise_level) {
                    detection_level = bck.star_level2;
                } else {
                    retries = 1;
                }
            }
            if (retries == 1) {
                detection_level = 30 * bck.noise_level;
            }
            if (retries == 0) {
                detection_level = 7 * bck.noise_level;
            }
            
            star_counter = 0;
            
            if (report_type > 0) {
                startext = "x,y,hfd,snr,flux,ra[0..360],dec[0..360]\n";
            }
            
            // Reset visited-star mask
            for (auto& row : img_sa[0]) {
                std::fill(row.begin(), row.end(), -1.0f);
            }
            
            for (auto fitsY = 0; fitsY < height5; ++fitsY) {
                for (auto fitsX = 0; fitsX < width5; ++fitsX) {
                    if (img_sa[0][fitsY][fitsX] <= 0.0f
                        && (img[0][fitsY][fitsX] - bck.backgr > detection_level)) {
                        auto hfd1 = 0.0;
                        auto fwhm = 0.0;
                        auto snr = 0.0;
                        auto flux = 0.0;
                        auto xc = 0.0;
                        auto yc = 0.0;
                        hfd_measure(img, fitsX, fitsY, 14, 99, 0.0,
                                    hfd1, fwhm, snr, flux, xc, yc);
                        // Accept as a star when HFD is within the ceiling, the
                        // SNR clears the threshold, and HFD exceeds hfd_min (the
                        // two-pixel floor that rejects hot pixels). The ceiling
                        // of 30 keeps extended sources — soft/large stars and
                        // compact galaxy cores — in the count.
                        if (hfd1 <= 30.0 && snr > snr_min && hfd1 > hfd_min) {
                            const auto xci = static_cast<int>(std::round(xc));
                            const auto yci = static_cast<int>(std::round(yc));

                            if (star_counter >= len) {
                                len += 1000;
                                hfd_list.resize(len);
                            }
                            hfd_list[star_counter++] = hfd1;

                            // Mark the whole circular star area as occupied to
                            // prevent double detections. Anything between 2.5 and
                            // 3.5 * HFD performs the same; a real star PSF has
                            // wider wings than a Gaussian predicts.
                            const auto radius = static_cast<int>(std::round(3.0 * hfd1));
                            const auto sqr_radius = radius * radius;
                            for (auto n = -radius; n <= radius; ++n) {
                                for (auto m = -radius; m <= radius; ++m) {
                                    auto j = n + yci;
                                    auto i = m + xci;
                                    if (j >= 0 && i >= 0 && j < height5 && i < width5
                                        && (m * m + n * n) <= sqr_radius) {
                                        img_sa[0][j][i] = 1.0f;
                                    }
                                }
                            }
                            
                            // Build CSV report line
                            if (report_type > 0) {
                                char buf[256];
                                if (head.cd1_1 == 0.0) {
                                    std::snprintf(buf, sizeof(buf),
                                        "%.4f,%.4f,%.4f,%d,%d\n",
                                        xc + 1.0, yc + 1.0, hfd1,
                                        static_cast<int>(std::round(snr)),
                                        static_cast<int>(std::round(flux)));
                                } else {
                                    auto ra = 0.0;
                                    auto decl = 0.0;
                                    pixel_to_celestial(head, xc + 1.0, yc + 1.0,
                                                       formalism, ra, decl);
                                    std::snprintf(buf, sizeof(buf),
                                        "%.4f,%.4f,%.4f,%d,%d,%.8f,%.8f\n",
                                        xc + 1.0, yc + 1.0, hfd1,
                                        static_cast<int>(std::round(snr)),
                                        static_cast<int>(std::round(flux)),
                                        ra  * 180.0 / std::numbers::pi,
                                        decl * 180.0 / std::numbers::pi);
                                }
                                startext += buf;
                            }
                        }
                    }
                }
            }
            
            --retries;
        } while (!(star_counter >= max_stars || retries < 0));
        
        // Compute median HFD. Uses the CCDciel-style median (averages the
        // middle two/three samples), matching the original analyse_image's
        // SMedian — not the plain nth-element s_median used for pixel stacks.
        if (star_counter > 0 && report_type <= 1) {
            hfd_median = astap::core::smedian(std::span<const double>{hfd_list},
                                              star_counter);
        } else {
            hfd_median = 99.0;
        }
    }
    
    if (report_type > 0 && csv_out != nullptr) {
        // Hand the assembled CSV back to the caller (CLI / test harness),
        // which owns the output path. Mirrors the original Pascal, which
        // wrote ChangeFileExt(filename2, '.csv').
        *csv_out = std::move(startext);
    }
}

namespace {

/// @brief Most common pixel value in a window — the modal background. Matches
///        astap_main.pas mode(): ignores black (< 1) and near-saturated
///        (>= max1 - 10) pixels. @p greylevels receives the count of distinct
///        occupied bins; @p ellipse restricts the window to its inscribed
///        ellipse (apply_most_common passes false).
[[nodiscard]] int mode(const ImageArray& img, bool ellipse, int colorm,
                       int xmin, int xmax, int ymin, int ymax, int max1,
                       int& greylevels) {
    const auto height3 = static_cast<int>(img[0].size());
    const auto width3  = static_cast<int>(img[0][0].size());

    max1 -= 10;  // do not measure saturated pixels
    if (xmin < 0) xmin = 0;
    if (xmax > width3 - 1) xmax = width3 - 1;
    if (ymin < 0) ymin = 0;
    if (ymax > height3 - 1) ymax = height3 - 1;
    greylevels = 0;
    if (max1 < 1) {
        return 0;  // degenerate (very low datamax): empty histogram
    }

    auto histogram = std::vector<int>(static_cast<std::size_t>(max1) + 1, 0);

    const auto centerX = (xmax + xmin) / 2.0;
    const auto centerY = (ymax + ymin) / 2.0;
    const auto a = (xmax - xmin - 1) / 2.0;
    const auto b = (ymax - ymin - 1) / 2.0;

    for (auto i = ymin; i <= ymax; ++i) {
        for (auto j = xmin; j <= xmax; ++j) {
            if (!ellipse ||
                ((j - centerX) * (j - centerX) / (a * a) +
                 (i - centerY) * (i - centerY) / (b * b) < 1.0)) {
                const auto val = static_cast<int>(std::lround(img[colorm][i][j]));
                if (val >= 1 && val < max1) {
                    ++histogram[static_cast<std::size_t>(val)];
                }
            }
        }
    }

    auto result = 0;
    auto value_count = 0;
    for (auto i = 1; i <= max1; ++i) {  // most common but ignore 0
        const auto val = histogram[static_cast<std::size_t>(i)];
        if (val != 0) ++greylevels;
        if (val > value_count) {
            value_count = val;
            result = i;
        }
    }
    return result;
}

}  // namespace

///----------------------------------------
/// MARK: apply_most_common
///----------------------------------------

void apply_most_common(const ImageArray& sourc, ImageArray& dest,
                       double datamax, int radius) {
    auto diameter = radius * 2;
    if (diameter <= 0 || sourc.empty() || sourc[0].empty()) {
        return;
    }
    const auto colors3 = static_cast<int>(sourc.size());
    const auto height3 = static_cast<int>(sourc[0].size());
    const auto width3  = static_cast<int>(sourc[0][0].size());

    for (auto k = 0; k < colors3; ++k) {
        for (auto fitsY = 0;
             fitsY <= static_cast<int>(std::round((height3 - 1.0) / diameter)); ++fitsY) {
            for (auto fitsX = 0;
                 fitsX <= static_cast<int>(std::round((width3 - 1.0) / diameter)); ++fitsX) {
                auto x = fitsX * diameter;
                auto y = fitsY * diameter;
                auto greylevels = 0;
                auto most_common = mode(sourc, /*ellipse=*/false, k,
                    x - radius, x + radius - 1,
                    y - radius, y + radius - 1,
                    static_cast<int>(datamax), greylevels);
                    
                // Fill the tile with the modal value
                for (auto i = -radius; i < radius; ++i) {
                    for (auto j = -radius; j < radius; ++j) {
                        auto x2 = x + i;
                        auto y2 = y + j;
                        if (x2 >= 0 && x2 < width3 && y2 >= 0 && y2 < height3) {
                            dest[k][y2][x2] = static_cast<float>(most_common);
                        }
                    }
                }
            }
        }
    }
}

///----------------------------------------
/// MARK: report_results (stub)
///----------------------------------------

void report_results([[maybe_unused]] std::string_view object_to_process,
                    [[maybe_unused]] std::string_view stack_info,
                    [[maybe_unused]] int object_counter,
                    [[maybe_unused]] int colorinfo) {
    // TODO: GUI-only side effects (populate stackmenu1.ListView5).
}

///----------------------------------------
/// MARK: apply_factors (stub)
///----------------------------------------

void apply_factors() {
    // TODO: reimplement as
    //   apply_factors(ImageArray&, float mulR, float mulG, float mulB,
    //                  float addR, float addG, float addB, bool accept_zero)
    // once the GUI-facing wrapper is rebuilt.
}

///----------------------------------------
/// MARK: remove_special_chars
///----------------------------------------

[[nodiscard]] std::string remove_special_chars(std::string_view s) {
    static constexpr char invalid[] = {'.', '\\', '/', '*', '"', ':', '|', '<', '>'};
    auto out = std::string{};
    out.reserve(s.size());
    for (auto c : s) {
        auto bad = false;
        for (auto k : invalid) {
            if (c == k) {
                bad = true;
                break;
            }
        }
        if (!bad) {
            out.push_back(c);
        }
    }
    return out;
}
 
} // namespace
