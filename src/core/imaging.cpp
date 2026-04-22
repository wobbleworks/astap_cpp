///----------------------------------------
///      @file imaging.cpp
///   @ingroup ASTAP++
///     @brief Image-manipulation utilities implementation.
///   @details Binning, mono conversion, rotation, duplication, coordinate
///            flipping, raw Bayer extraction, histogram computation, and
///            luminance/saturation stretch.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP, MPL-2.0)
///            by John Stephen.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "imaging.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "globals.h"

///----------------------------------------
namespace astap::core {
///----------------------------------------

using astap::roworder;
using astap::filename2;

/// MARK: Forward declarations

// Header-update memo helpers from util.h / fits.h.
void update_integer(std::vector<std::string>& memo,
                    const std::string& keyword,
                    const std::string& comment,
                    int value);
void update_float(std::vector<std::string>& memo,
                  const std::string& keyword,
                  const std::string& comment,
                  bool exponential,
                  double value);
void update_text(std::vector<std::string>& memo,
                 const std::string& keyword,
                 const std::string& value);
void remove_key(std::vector<std::string>& memo,
                const std::string& keyword,
                bool all);
void add_text(std::vector<std::string>& memo,
              const std::string& keyword,
              const std::string& text);
              
// FITS I/O from core/fits.h.
bool load_fits(const std::filesystem::path& filename,
               bool light,
               bool update_memo,
               int hdu_index,
               std::vector<std::string>& memo,
               Header& head,
               ImageArray& img);
bool save_fits(const ImageArray& img,
               const std::vector<std::string>& memo,
               const std::filesystem::path& filename,
               int nrbits,
               bool overwrite);
               
// Bayer-pattern detection and old-style WCS conversion.
int get_demosaic_pattern();
void old_to_new_WCS(Header& head);

// HSV <-> RGB conversion used by stretch_image.
void RGB2HSV(float r, float g, float b, float& h, float& s, float& v);
void HSV2RGB(float h, float s, float v, float& r, float& g, float& b);

// Accurate raster rotation with canvas expansion.
void raster_rotate(double angle, double cx, double cy, ImageArray& img);

// GUI-side stubs.
void memo2_message(const std::string& msg);
void remove_solution(bool keep_wcs);

/// MARK: Module-level state

std::array<std::array<std::uint32_t, 65536>, 3> histogram{};
std::uint32_t his_total_red = 0;
std::array<int, 3> his_mean{};
std::array<float, 32769> stretch_c{};
int hist_range = 255;
// Canonical `cwhite` lives in `astap::` (globals.cpp). A duplicate definition
// used to live here as `astap::core::cwhite`, silently shadowing the loader's
// writes; resolved via `using astap::cwhite;` below.
using astap::cwhite;
bool stretch_on = false;
float saturation_factor = 1.0F;

namespace {

// Always non-negative modulo.
[[nodiscard]] double fnmodulo(double x, double m) noexcept {
    auto r = std::fmod(x, m);
    if (r < 0.0) {
        r += m;
    }
    return r;
}

// Truncate toward zero.
[[nodiscard]] inline int trunc_to_int(double v) noexcept {
    return static_cast<int>(v);
}

// Allocate a [c][h][w] image, zero-initialised.
[[nodiscard]] ImageArray make_image(int c, int h, int w) {
    return ImageArray(static_cast<std::size_t>(c),
                      std::vector<std::vector<float>>(
                          static_cast<std::size_t>(h),
                          std::vector<float>(static_cast<std::size_t>(w), 0.0F)));
}

// Replace file extension with a forced suffix.
[[nodiscard]] std::filesystem::path change_file_ext(const std::filesystem::path& p,
                                                    const std::string& new_ext) {
    auto out = p;
    out.replace_extension();
    out += new_ext;
    return out;
}

// Check whether s contains sub.
[[nodiscard]] inline bool contains(const std::string& s, const std::string& sub) noexcept {
    return s.find(sub) != std::string::npos;
}
 
}  // namespace

/// MARK: bin_X2X3X4

void bin_X2X3X4(ImageArray& img,
                Header& head,
                std::vector<std::string>& memo,
                int binfactor,
                const std::string& filename2_in) {
    binfactor = std::min(4, binfactor);
    const auto w = trunc_to_int(static_cast<double>(head.width) / binfactor);
    const auto h = trunc_to_int(static_cast<double>(head.height) / binfactor);
    auto img_temp2 = make_image(head.naxis3, h, w);
    
    if (binfactor == 2) {
        for (auto k = 0; k < head.naxis3; ++k) {
            for (auto fitsY = 0; fitsY < h; ++fitsY) {
                for (auto fitsX = 0; fitsX < w; ++fitsX) {
                    img_temp2[k][fitsY][fitsX] =
                        (img[k][fitsY * 2][fitsX * 2] +
                         img[k][fitsY * 2 + 1][fitsX * 2] +
                         img[k][fitsY * 2][fitsX * 2 + 1] +
                         img[k][fitsY * 2 + 1][fitsX * 2 + 1]) / 4.0F;
                }
            }
        }
    } else if (binfactor == 3) {
        for (auto k = 0; k < head.naxis3; ++k) {
            for (auto fitsY = 0; fitsY < h; ++fitsY) {
                for (auto fitsX = 0; fitsX < w; ++fitsX) {
                    img_temp2[k][fitsY][fitsX] =
                        (img[k][fitsY * 3    ][fitsX * 3    ] +
                         img[k][fitsY * 3    ][fitsX * 3 + 1] +
                         img[k][fitsY * 3    ][fitsX * 3 + 2] +
                         img[k][fitsY * 3 + 1][fitsX * 3    ] +
                         img[k][fitsY * 3 + 1][fitsX * 3 + 1] +
                         img[k][fitsY * 3 + 1][fitsX * 3 + 2] +
                         img[k][fitsY * 3 + 2][fitsX * 3    ] +
                         img[k][fitsY * 3 + 2][fitsX * 3 + 1] +
                         img[k][fitsY * 3 + 2][fitsX * 3 + 2]) / 9.0F;
                }
            }
        }
    } else {
        // Bin 4x4
        for (auto k = 0; k < head.naxis3; ++k) {
            for (auto fitsY = 0; fitsY < h; ++fitsY) {
                for (auto fitsX = 0; fitsX < w; ++fitsX) {
                    auto sum = 0.0F;
                    for (auto dy = 0; dy < 4; ++dy) {
                        for (auto dx = 0; dx < 4; ++dx) {
                            sum += img[k][fitsY * 4 + dy][fitsX * 4 + dx];
                        }
                    }
                    img_temp2[k][fitsY][fitsX] = sum / 16.0F;
                }
            }
        }
    }
    
    img = std::move(img_temp2);
    head.width = w;
    head.height = h;
    
    // Update FITS header entries
    update_integer(memo, "NAXIS1  =", " / length of x axis                               ", head.width);
    update_integer(memo, "NAXIS2  =", " / length of y axis                               ", head.height);
    
    if (head.crpix1 != 0) {
        head.crpix1 /= binfactor;
        update_float(memo, "CRPIX1  =", " / X of reference pixel                           ", false, head.crpix1);
    }
    if (head.crpix2 != 0) {
        head.crpix2 /= binfactor;
        update_float(memo, "CRPIX2  =", " / Y of reference pixel                           ", false, head.crpix2);
    }
    if (head.cdelt1 != 0) {
        head.cdelt1 *= binfactor;
        update_float(memo, "CDELT1  =", " / X pixel size (deg)                             ", false, head.cdelt1);
    }
    if (head.cdelt2 != 0) {
        head.cdelt2 *= binfactor;
        update_float(memo, "CDELT2  =", " / Y pixel size (deg)                             ", false, head.cdelt2);
    }
    if (head.cd1_1 != 0) {
        head.cd1_1 *= binfactor;
        head.cd1_2 *= binfactor;
        head.cd2_1 *= binfactor;
        head.cd2_2 *= binfactor;
        update_float(memo, "CD1_1   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, head.cd1_1);
        update_float(memo, "CD1_2   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, head.cd1_2);
        update_float(memo, "CD2_1   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, head.cd2_1);
        update_float(memo, "CD2_2   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, head.cd2_2);
    }
    
    head.xbinning *= binfactor;
    head.ybinning *= binfactor;
    update_integer(memo, "XBINNING=", " / Binning factor in width                         ",
                   static_cast<int>(std::lround(head.xbinning)));
    update_integer(memo, "YBINNING=", " / Binning factor in height                        ",
                   static_cast<int>(std::lround(head.ybinning)));
                   
    if (head.xpixsz != 0) {
        head.xpixsz *= binfactor;
        head.ypixsz *= binfactor;
        update_float(memo, "XPIXSZ  =", " / Pixel width in microns (after binning)          ", false, head.xpixsz);
        update_float(memo, "YPIXSZ  =", " / Pixel height in microns (after binning)         ", false, head.ypixsz);
        update_float(memo, "PIXSIZE1=", " / Pixel width in microns (after binning)          ", false, head.xpixsz);
        update_float(memo, "PIXSIZE2=", " / Pixel height in microns (after binning)         ", false, head.ypixsz);
    }
    
    const auto fact = std::to_string(binfactor);
    add_text(memo, "HISTORY   ", "BIN" + fact + "x" + fact + " version of " + filename2_in);
}

/// MARK: convert_mono

void convert_mono(ImageArray& img, Header& head) {
    if (head.naxis3 < 3) {
        return;
    }
    
    memo2_message("Converting to mono.");
    auto img_temp = make_image(1, head.height, head.width);
    
    for (auto fitsY = 0; fitsY < head.height; ++fitsY) {
        for (auto fitsX = 0; fitsX < head.width; ++fitsX) {
            img_temp[0][fitsY][fitsX] =
                (img[0][fitsY][fitsX] + img[1][fitsY][fitsX] + img[2][fitsY][fitsX]) / 3.0F;
        }
    }
    
    head.naxis = 2;
    head.naxis3 = 1;
    img = std::move(img_temp);
}

/// MARK: rotate_arbitrary

void rotate_arbitrary(double angle,
                      double flipped_view,
                      double flipped_image,
                      ImageArray& img,
                      Header& head,
                      std::vector<std::string>& memo) {
    const auto centerxs = head.width / 2.0;
    const auto centerys = head.height / 2.0;
    
    raster_rotate(flipped_view * angle, centerxs, centerys, img);
    
    head.width = static_cast<int>(img[0][0].size());
    head.height = static_cast<int>(img[0].size());
    
    // Update dimensions
    update_integer(memo, "NAXIS1  =", " / length of x axis                               ", head.width);
    update_integer(memo, "NAXIS2  =", " / length of y axis                               ", head.height);
    
    if (head.cd1_1 != 0) {
        if ((head.crpix1 != 0.5 + centerxs) || (head.crpix2 != 0.5 + centerys)) {
            // Reference is not center; drop the solution.
            remove_solution(true);
        }
        head.crota2 = fnmodulo(head.crota2 + angle * flipped_image * flipped_view, 360.0);
        head.crota1 = fnmodulo(head.crota1 + angle * flipped_image * flipped_view, 360.0);
        head.crpix1 = head.width / 2.0;
        head.crpix2 = head.height / 2.0;
        old_to_new_WCS(head);
        
        // Update CD matrix
        update_float(memo, "CD1_1   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, head.cd1_1);
        update_float(memo, "CD1_2   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, head.cd1_2);
        update_float(memo, "CD2_1   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, head.cd2_1);
        update_float(memo, "CD2_2   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, head.cd2_2);
        
        // Update reference pixel
        update_float(memo, "CRPIX1  =", " / X of reference pixel                           ", false, head.crpix1);
        update_float(memo, "CRPIX2  =", " / Y of reference pixel                           ", false, head.crpix2);
        
        // Update rotation angles
        update_float(memo, "CROTA1  =", " / Image twist X axis (deg)                       ", false, head.crota1);
        update_float(memo, "CROTA2  =", " / Image twist Y axis (deg) E of N if not flipped.", false, head.crota2);
    }
    
    remove_key(memo, "ANNOTATE", true);
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.2f", angle);
    add_text(memo, "HISTORY   ", std::string("Rotated CCW by angle ") + buf);
}

/// MARK: duplicate

ImageArray duplicate(const ImageArray& img) {
    // std::vector copy already deep-copies; this matches the original helper's
    // intent (force a real allocation rather than aliasing).
    return img;
}

/// MARK: flip

void flip(int x1, int y1, int& x2, int& y2,
          const Header& head, bool flip_horizontal, bool flip_vertical) {
    if (flip_horizontal) {
        x2 = (head.width - 1) - x1;
    } else {
        x2 = x1;
    }
    
    if (!flip_vertical) {
        y2 = (head.height - 1) - y1;
    } else {
        y2 = y1;
    }
}

/// MARK: extract_raw_colour_to_file

std::string extract_raw_colour_to_file(const std::filesystem::path& filename,
                                       const std::string& filtern,
                                       int xp,
                                       int yp,
                                       std::vector<std::string>& memo) {
    auto head = Header{};
    auto img = ImageArray{};
    if (!load_fits(filename, /*light=*/true, /*update_memo=*/true, /*hdu=*/0, memo, head, img)) {
        return {};
    }
    
    // Validate: must be a single-channel raw image without prior extraction
    if (!(!contains(head.filter_name, "TR") &&
          !contains(head.filter_name, "TG") &&
          !contains(head.filter_name, "TB") &&
          head.naxis3 == 1)) {
        if (head.naxis3 > 1) {
            memo2_message("Skipped COLOUR image " + filename.string() +
                          ", Raw red, green or blue pixel extraction is only possible for raw images.");
        } else {
            memo2_message("Skipped image " + filename.string() +
                          ", FILTER indicates earlier extraction!");
        }
        return {};
    }
    
    constexpr auto ratio = 0.5;
    const auto w = trunc_to_int(head.width / 2.0);
    const auto h = trunc_to_int(head.height / 2.0);
    auto img_temp11 = make_image(1, h, w);
    
    // Determine Bayer cell offsets from the demosaic pattern
    const auto pattern = get_demosaic_pattern();
    auto get_green = false;
    auto xp2 = 0;
    auto yp2 = 0;
    
    if (filtern == "TR") {
        switch (pattern) {
            case 0: xp = 2; yp = 1; break;  // GRBG
            case 1: xp = 2; yp = 2; break;  // BGGR
            case 2: xp = 1; yp = 1; break;  // RGGB
            case 3: xp = 1; yp = 2; break;  // GBRG
            default: break;
        }
    } else if (filtern == "TB") {
        switch (pattern) {
            case 0: xp = 1; yp = 2; break;
            case 1: xp = 1; yp = 1; break;
            case 2: xp = 2; yp = 2; break;
            case 3: xp = 2; yp = 1; break;
            default: break;
        }
    } else if (filtern == "TG") {
        // Green averages two pixels
        get_green = true;
        switch (pattern) {
            case 0: xp = 1; yp = 1; xp2 = 2; yp2 = 2; break;  // GRBG
            case 1: xp = 2; yp = 1; xp2 = 1; yp2 = 2; break;  // BGGR
            case 2: xp = 2; yp = 1; xp2 = 1; yp2 = 2; break;  // RGGB
            case 3: xp = 1; yp = 1; xp2 = 2; yp2 = 2; break;  // GBRG
            default: break;
        }
    }
    
    // Adjust pattern label for BOTTOM-UP row order
    auto pattern2 = pattern;
    if (contains(roworder, "BOT")) {
        switch (pattern) {
            case 0: pattern2 = 1; break;
            case 1: pattern2 = 0; break;
            case 2: pattern2 = 3; break;
            case 3: pattern2 = 2; break;
            default: break;
        }
    }
    
    // Log the extraction
    const auto filtern_letter = filtern.size() >= 2 ? std::string(1, filtern[1]) : std::string();
    switch (pattern2) {
        case 0: memo2_message("GRBG => " + filtern_letter); break;
        case 1: memo2_message("BGGR => " + filtern_letter); break;
        case 2: memo2_message("RGGB => " + filtern_letter); break;
        case 3: memo2_message("GBRG => " + filtern_letter); break;
        default: break;
    }
    
    // Extract the Bayer plane
    for (auto fitsY = 0; fitsY < h; ++fitsY) {
        for (auto fitsX = 0; fitsX < w; ++fitsX) {
            auto val = img[0][fitsY * 2 + yp - 1][fitsX * 2 + xp - 1];
            if (get_green) {
                val = (val + img[0][fitsY * 2 + yp2 - 1][fitsX * 2 + xp2 - 1]) / 2.0F;
            }
            img_temp11[0][fitsY][fitsX] = val;
        }
    }
    
    head.width = w;
    head.height = h;
    
    // Update header dimensions
    update_integer(memo, "NAXIS1  =", " / length of x axis                               ", head.width);
    update_integer(memo, "NAXIS2  =", " / length of y axis                               ", head.height);
    update_integer(memo, "NAXIS1  =", " / length of x axis                               ", head.width);
    update_integer(memo, "NAXIS2  =", " / length of y axis                               ", head.height);
    
    // Update reference pixel
    if (head.crpix1 != 0) {
        head.crpix1 *= ratio;
        update_float(memo, "CRPIX1  =", " / X of reference pixel                           ", false, head.crpix1);
    }
    if (head.crpix2 != 0) {
        head.crpix2 *= ratio;
        update_float(memo, "CRPIX2  =", " / Y of reference pixel                           ", false, head.crpix2);
    }
    
    // Update pixel scale
    if (head.cdelt1 != 0) {
        head.cdelt1 /= ratio;
        update_float(memo, "CDELT1  =", " / X pixel size (deg)                             ", false, head.cdelt1);
    }
    if (head.cdelt2 != 0) {
        head.cdelt2 /= ratio;
        update_float(memo, "CDELT2  =", " / Y pixel size (deg)                             ", false, head.cdelt2);
    }
    
    // Update CD matrix
    if (head.cd1_1 != 0) {
        head.cd1_1 /= ratio;
        head.cd1_2 /= ratio;
        head.cd2_1 /= ratio;
        head.cd2_2 /= ratio;
        update_float(memo, "CD1_1   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, head.cd1_1);
        update_float(memo, "CD1_2   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, head.cd1_2);
        update_float(memo, "CD2_1   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, head.cd2_1);
        update_float(memo, "CD2_2   =", " / CD matrix to convert (x,y) to (Ra, Dec)        ", false, head.cd2_2);
    }
    
    // Update binning
    head.xbinning /= ratio;
    head.ybinning /= ratio;
    update_float(memo, "XBINNING=", " / Binning factor in width                         ", false, head.xbinning);
    update_float(memo, "YBINNING=", " / Binning factor in height                        ", false, head.ybinning);
    
    // Update pixel size
    if (head.xpixsz != 0) {
        head.xpixsz /= ratio;
        head.ypixsz /= ratio;
        update_float(memo, "XPIXSZ  =", " / Pixel width in microns (after binning)          ", false, head.xpixsz);
        update_float(memo, "YPIXSZ  =", " / Pixel height in microns (after binning)         ", false, head.ypixsz);
        update_float(memo, "PIXSIZE1=", " / Pixel width in microns (after binning)          ", false, head.xpixsz);
        update_float(memo, "PIXSIZE2=", " / Pixel height in microns (after binning)         ", false, head.ypixsz);
    }
    
    // Update filter metadata
    add_text(memo, "HISTORY   ", "One raw colour extracted.");
    remove_key(memo, "BAYERPAT=", false);
    auto padded = "'" + filtern + "'" + std::string(20, ' ');
    padded = padded.substr(0, 21) + "/ Filter name";
    update_text(memo, "FILTER  =", padded);
    
    // Save the extracted plane
    img = std::move(img_temp11);
    auto out = change_file_ext(filename, std::string("_") + filtern + ".fit");
    if (!save_fits(img, memo, out, /*nrbits=*/16, /*overwrite=*/true)) {
        return {};
    }
    return out.string();
}

/// MARK: split_raw

std::size_t split_raw(const std::vector<std::filesystem::path>& filenames,
                      int xp,
                      int yp,
                      const std::string& filtern,
                      std::vector<std::string>& memo) {
    // TODO host hook: backup_img / restore_img + Screen.Cursor.
    auto ok = size_t{0};
    for (const auto& f : filenames) {
        // TODO host hook: Application.ProcessMessages, esc_pressed.
        if (!extract_raw_colour_to_file(f, filtern, xp, yp, memo).empty()) {
            ++ok;
        }
    }
    return ok;
}

/// MARK: get_hist

void get_hist(int colour, const ImageArray& img) {
    if (static_cast<std::size_t>(colour) + 1 > img.size()) {
        colour = 0;
    }
    
    // Clear the histogram for this channel
    for (auto i = 0; i <= 65535; ++i) {
        histogram[colour][i] = 0;
    }
    
    auto his_total = std::uint32_t{0};
    auto total_value = 0.0;
    auto count = 1;  // avoid divide by zero
    const auto width5 = static_cast<int>(img[0][0].size());
    const auto height5 = static_cast<int>(img[0].size());
    
    // Exclude LibRaw border
    const auto offsetW = trunc_to_int(width5 * 0.042);
    const auto offsetH = trunc_to_int(height5 * 0.015);
    
    // Accumulate histogram
    for (auto i = offsetH; i <= height5 - 1 - offsetH; ++i) {
        for (auto j = offsetW; j <= width5 - 1 - offsetW; ++j) {
            const auto col = static_cast<int>(std::lround(img[colour][i][j]));
            if ((col >= 1) && (col < 65000)) {
                ++histogram[colour][col];
                ++his_total;
                total_value += col;
                ++count;
            }
        }
    }
    
    if (colour == 0) {
        his_total_red = his_total;
    }
    his_mean[colour] = static_cast<int>(std::lround(total_value / count));
}

/// MARK: use_histogram

void use_histogram(const ImageArray& img,
                   bool update_hist,
                   int range_index,
                   const Header& head,
                   int& minm,
                   int& maxm) {
    const auto number_colors = static_cast<int>(img.size());
    
    // Recompute histograms if requested
    if (update_hist) {
        get_hist(0, img);
        if (number_colors > 1) {
            get_hist(1, img);
        }
        if (number_colors > 2) {
            get_hist(2, img);
        }
    }
    
    const auto max_range = static_cast<int>(
        std::lround(std::min(head.datamax_org, 65535.0)));
        
    auto above_R = 0.001;
    minm = 0;
    maxm = 0;
    
    // Determine range based on index
    switch (range_index) {
        case -1:
        case 0:
        case 1:
            above_R = 0.001;
            break;
        case 2:
        case 3:
            above_R = 0.003;
            break;
        case 4:
        case 5:
            above_R = 0.01;
            break;
        case 6:
        case 7:
            minm = static_cast<int>(std::lround(head.datamin_org));
            maxm = static_cast<int>(std::lround(head.datamax_org));
            break;
        case 8:
            minm = static_cast<int>(std::lround(max_range * 0.95));
            maxm = max_range;
            break;
        case 9:
            minm = 0;
            maxm = 65535;
            break;
        default:
            break;
    }
    
    if (range_index <= 5) {
        // Auto-detect mode
        minm = 0;
        maxm = 0;
        auto above = 0.0;
        auto histo_peak_position = 0;
        auto histo_peakR = std::int64_t{-99999999};
        
        // Find the peak position
        for (auto i = 1; i <= max_range - 1; ++i) {
            if (static_cast<std::int64_t>(histogram[0][i]) > histo_peakR) {
                histo_peakR = histogram[0][i];
                histo_peak_position = i;
            }
        }
        
        // Find minimum from peak downward
        auto i = histo_peak_position;
        while (minm == 0 && i > 0) {
            --i;
            if (histogram[0][i] < 0.1 * histogram[0][histo_peak_position]) {
                minm = i;
            }
        }
        
        // Find maximum from top downward
        i = max_range;
        while (maxm == 0 && i > minm + 1) {
            --i;
            above += histogram[0][i];
            if (above > above_R * his_total_red) {
                maxm = i;
            }
        }
    }
    
    hist_range = static_cast<int>(std::lround(
        std::min(65535.0, std::min(head.datamax_org, 2.0 * maxm))));
        
    // TODO host: paint mainwindow.histogram1 (logarithmic-scaled per-channel bars).
}

/// MARK: stretch_image

ImageArray stretch_image(const ImageArray& img, const Background& bck) {
    const auto colours5 = static_cast<int>(img.size());
    const auto height5 = static_cast<int>(img[0].size());
    const auto width5 = static_cast<int>(img[0][0].size());
    
    auto result = make_image(colours5, height5, width5);
    const auto sat_factor = saturation_factor;
    const auto denom = (cwhite - bck.backgr);
    
    for (auto fitsY = 0; fitsY < height5; ++fitsY) {
        for (auto fitsX = 0; fitsX < width5; ++fitsX) {
            if (colours5 == 3) {
                auto col_r = img[0][fitsY][fitsX];
                auto col_g = img[1][fitsY][fitsX];
                auto col_b = img[2][fitsY][fitsX];
                
                auto colrr = static_cast<float>((col_r - bck.backgr) / denom);
                auto colgg = static_cast<float>((col_g - bck.backgr) / denom);
                auto colbb = static_cast<float>((col_b - bck.backgr) / denom);
                
                // Apply saturation adjustment
                if (sat_factor != 1.0F) {
                    auto h = 0.0F, s = 0.0F, v = 0.0F;
                    RGB2HSV(colrr, colgg, colbb, h, s, v);
                    HSV2RGB(h, s * sat_factor, v, colrr, colgg, colbb);
                }
                
                // Clamp to minimum
                if (colrr <= 1e-11F) {
                    colrr = 1e-11F;
                }
                if (colgg <= 1e-11F) {
                    colgg = 1e-11F;
                }
                if (colbb <= 1e-11F) {
                    colbb = 1e-11F;
                }
                
                // Normalize if any channel exceeds 1.0
                auto largest = colrr;
                if (colgg > largest) {
                    largest = colgg;
                }
                if (colbb > largest) {
                    largest = colbb;
                }
                if (largest > 1.0F) {
                    colrr /= largest;
                    colgg /= largest;
                    colbb /= largest;
                    largest = 1.0F;
                }
                
                // Apply stretch or linear scale
                if (stretch_on) {
                    const auto luminance = (colrr + colgg + colbb) / 3.0F;
                    const auto idx = std::min(32768, std::max(0,
                        static_cast<int>(32768.0F * luminance)));
                    const auto luminance_stretched = stretch_c[idx];
                    auto factor = luminance_stretched / luminance;
                    if (factor * largest > 1.0F) {
                        factor = 1.0F / largest;
                    }
                    col_r = std::round(colrr * factor * 65535.0F);
                    col_g = std::round(colgg * factor * 65535.0F);
                    col_b = std::round(colbb * factor * 65535.0F);
                } else {
                    col_r = std::round(65535.0F * colrr);
                    col_g = std::round(65535.0F * colgg);
                    col_b = std::round(65535.0F * colbb);
                }
                
                result[0][fitsY][fitsX] = col_r;
                result[1][fitsY][fitsX] = col_g;
                result[2][fitsY][fitsX] = col_b;
            } else {
                // Mono
                auto col_r = img[0][fitsY][fitsX];
                auto colrr = static_cast<float>((col_r - bck.backgr) / denom);
                
                // Clamp
                if (colrr <= 1e-11F) {
                    colrr = 1e-11F;
                }
                if (colrr > 1.0F) {
                    colrr = 1.0F;
                }
                
                // Apply stretch or linear scale
                if (stretch_on) {
                    const auto idx = std::min(32768, std::max(0,
                        static_cast<int>(32768.0F * colrr)));
                    col_r = std::round(65535.0F * stretch_c[idx]);
                } else {
                    col_r = std::round(65535.0F * colrr);
                }
                
                result[0][fitsY][fitsX] = col_r;
            }
        }
    }
    return result;
}
 
} // namespace
