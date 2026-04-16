///----------------------------------------
///      @file photometry.cpp
///   @ingroup ASTAP++
///     @brief Photometry and star-measurement implementation.
///   @details Faithful C++23 port of the photometry primitives from
///            astap_main.pas. Loop structure, numeric ordering, and even the
///            redundant double-zero-init in HFD are preserved so the math
///            matches the original output.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP), MPL-2.0.
///            C++ port by John Stephen.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "photometry.h"

#include "globals.h"
#include "imaging.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <format>
#include <numbers>
#include <string>
#include <vector>

///----------------------------------------
namespace astap::core {
///----------------------------------------

// Forward declarations of helpers ported in other modules.
double Smedian(std::vector<double>& data, int n);

// GUI-side hooks the original source called inline. Provide weak no-op stubs
// so this translation unit links standalone.
[[maybe_unused]] static void memo2_message([[maybe_unused]] const std::string& msg) {}

// Forward declarations from sibling modules. Stubs here so this file is
// self-contained at the type level; the linker resolves to the real
// implementations once the rest of ASTAP++ is wired up. get_hist comes
// from core/imaging.h (included above).
void analyse_image(const ImageArray& img, Header& head, int snr_min,
                   int report_mode,
                   int& hfd_counter, Background& bck, double& hfd_med);
void plot_and_measure_stars(const ImageArray& img,
                            std::vector<std::string>& memo, Header& head,
                            bool calibration, bool plot_stars,
                            bool report_lim_magnitude);
                            
namespace {

constexpr auto kPi = std::numbers::pi;

///----------------------------------------
/// @brief Fractional part of a double.
///----------------------------------------

[[nodiscard]] inline double frac(double x) noexcept {
    auto ip = 0.0;
    return std::modf(x, &ip);
}
 
}  // namespace

/// MARK: - HFD

void HFD(const ImageArray& img,
         int x1, int y1, int rs,
         double aperture_small, double adu_e,
         double xbinning,
         HfdResult& result,
         HfdScratch& scratch) {
    constexpr auto kMaxRi = 74;
    constexpr auto kSamplePoints = 5;
    
    const auto width5  = static_cast<int>(img[0][0].size());
    const auto height5 = static_cast<int>(img[0].size());
    
    // Local sub-pixel sampler
    auto value_subpixel = [&](double xs, double ys) -> double {
        const auto x_trunc = static_cast<int>(xs);
        const auto y_trunc = static_cast<int>(ys);
        if (x_trunc <= 0 || x_trunc >= (width5 - 2) ||
            y_trunc <= 0 || y_trunc >= (height5 - 2)) {
            return 0.0;
        }
        const auto x_frac = frac(xs);
        const auto y_frac = frac(ys);
        auto r = (img[0][y_trunc    ][x_trunc    ]) * (1 - x_frac) * (1 - y_frac);
        r     += (img[0][y_trunc    ][x_trunc + 1]) * (    x_frac) * (1 - y_frac);
        r     += (img[0][y_trunc + 1][x_trunc    ]) * (1 - x_frac) * (    y_frac);
        r     += (img[0][y_trunc + 1][x_trunc + 1]) * (    x_frac) * (    y_frac);
        return r;
    };
    
    auto annulus_width = (aperture_small < 99) ? 3 : 1;
    auto r1_square = rs * rs;
    auto r2 = rs + annulus_width;
    auto r2_square = r2 * r2;
    
    if ((x1 - r2 <= 0) || (x1 + r2 >= width5 - 1) ||
        (y1 - r2 <= 0) || (y1 + r2 >= height5 - 1)) {
        result.hfd = 999;
        result.snr = 0;
        return;
    }
    
    auto valmax = 0.0;
    result.hfd = 999;
    result.snr = 0;
    
    auto r_aperture = scratch.r_aperture;
    auto star_bg = scratch.star_bg;
    auto sd_bg   = scratch.sd_bg;
    auto xc = result.xc;
    auto yc = result.yc;
    auto boxed = false;
    
    // Outer try/catch mirrors the original catch-all exception handler
    try {
        // Compute annulus background
        auto counter = 0;
        auto background = std::array<double, 1001>{};
        for (auto i = -r2; i <= r2; ++i) {
            for (auto j = -r2; j <= r2; ++j) {
                const auto distance = i * i + j * j;
                if ((distance > r1_square) && (distance <= r2_square)) {
                    background[counter] = img[0][y1 + j][x1 + i];
                    ++counter;
                }
            }
        }
        
        auto bg_vec = std::vector<double>(background.begin(), background.begin() + counter);
        star_bg = Smedian(bg_vec, counter);
        for (auto i = 0; i < counter; ++i) {
            bg_vec[i] = std::abs(bg_vec[i] - star_bg);
        }
        auto mad_bg = Smedian(bg_vec, counter);
        sd_bg = mad_bg * 1.4826;
        sd_bg = std::max(sd_bg, 1.0);
        
        auto SumVal = 0.0;
        auto SumValX = 0.0;
        auto SumValY = 0.0;
        auto signal_counter = 0;
        
        // Iterative centroid refinement
        do {
            SumVal = 0;
            SumValX = 0;
            SumValY = 0;
            signal_counter = 0;
            
            for (auto i = -rs; i <= rs; ++i) {
                for (auto j = -rs; j <= rs; ++j) {
                    auto val = (img[0][y1 + j][x1 + i]) - star_bg;
                    if (val > 3.0 * sd_bg) {
                        SumVal  += val;
                        SumValX += val * i;
                        SumValY += val * j;
                        ++signal_counter;
                    }
                }
            }
            
            if (SumVal <= 12 * sd_bg) {
                scratch.star_bg = star_bg;
                scratch.sd_bg = sd_bg;
                return;  // hfd already 999
            }
            
            const auto Xg = SumValX / SumVal;
            const auto Yg = SumValY / SumVal;
            xc = x1 + Xg;
            yc = y1 + Yg;
            
            if ((xc - rs < 0) || (xc + rs > width5 - 1) ||
                (yc - rs < 0) || (yc + rs > height5 - 1)) {
                scratch.star_bg = star_bg;
                scratch.sd_bg = sd_bg;
                result.xc = xc;
                result.yc = yc;
                return;
            }
            
            boxed = (signal_counter >= (2.0 / 9.0) * (rs + rs + 1) * (rs + rs + 1));
            
            if (!boxed) {
                if (rs > 4) {
                    rs -= 2;
                } else {
                    rs -= 1;
                }
            }
            
            if (signal_counter <= 1) {
                scratch.star_bg = star_bg;
                scratch.sd_bg = sd_bg;
                result.xc = xc;
                result.yc = yc;
                return;
            }
        } while (!(boxed || (rs <= 1)));
        
        // Add some space
        rs += 2;
        
        // Build signal histogram from centre of gravity
        auto distance_histogram = std::array<int, kMaxRi + 1>{};
        for (auto i = 0; i <= rs; ++i) {
            distance_histogram[i] = 0;
        }
        for (auto i = -rs; i <= rs; ++i) {
            for (auto j = -rs; j <= rs; ++j) {
                const auto distance = static_cast<int>(std::lround(std::sqrt(
                    static_cast<double>(i * i + j * j))));
                if (distance <= rs) {
                    const auto val = value_subpixel(xc + i, yc + j) - star_bg;
                    if (val > 3.0 * sd_bg) {
                        distance_histogram[distance] += 1;
                        if (val > valmax) {
                            valmax = val;
                        }
                    }
                }
            }
        }
        
        // Determine aperture radius from the distance histogram
        r_aperture = -1;
        auto distance_top_value = 0;
        auto HistStart = false;
        auto illuminated_pixels = 0;
        do {
            ++r_aperture;
            illuminated_pixels += distance_histogram[r_aperture];
            if (distance_histogram[r_aperture] > 0) {
                HistStart = true;
            }
            if (distance_top_value < distance_histogram[r_aperture]) {
                distance_top_value = distance_histogram[r_aperture];
            }
        } while (!((r_aperture >= rs) ||
                   (HistStart && (distance_histogram[r_aperture] <= 0.1 * distance_top_value))));
                   
        if (r_aperture >= rs) {
            scratch.star_bg = star_bg;
            scratch.sd_bg = sd_bg;
            scratch.r_aperture = r_aperture;
            result.xc = xc;
            result.yc = yc;
            return;
        }
        
        if ((r_aperture > 2) &&
            (illuminated_pixels < 0.35 * (r_aperture + r_aperture - 2) *
                                          (r_aperture + r_aperture - 2))) {
            scratch.star_bg = star_bg;
            scratch.sd_bg = sd_bg;
            scratch.r_aperture = r_aperture;
            result.xc = xc;
            result.yc = yc;
            return;
        }
    } catch (...) {
        // Original catch-all -- swallow
    }
    
    // Get HFD (the intentional double-zero-init matches the original)
    auto SumVal = 0.0;
    auto Sumval_small = 0.0;
    auto SumValR = 0.0;
    auto pixel_counter = 0.0;
    SumVal = 0;
    Sumval_small = 0;
    SumValR = 0;
    pixel_counter = 0;
    
    auto flux = 0.0;
    auto radius = 0.0;
    if (r_aperture < aperture_small) {
        for (auto i = -r_aperture; i <= r_aperture; ++i) {
            for (auto j = -r_aperture; j <= r_aperture; ++j) {
                auto Val = value_subpixel(xc + i, yc + j) - star_bg;
                if (Val >= valmax * 0.5) {
                    pixel_counter += 1;
                }
                const auto r = std::sqrt(static_cast<double>(i * i + j * j));
                SumVal  += Val;
                SumValR += Val * r;
            }
        }
        flux = std::max(SumVal, 0.00001);
        radius = r_aperture;
    } else {
        for (auto i = -r_aperture * kSamplePoints; i <= r_aperture * kSamplePoints; ++i) {
            for (auto j = -r_aperture * kSamplePoints; j <= r_aperture * kSamplePoints; ++j) {
                const auto dx = static_cast<double>(i) / kSamplePoints;
                const auto dy = static_cast<double>(j) / kSamplePoints;
                auto Val = value_subpixel(xc + dx, yc + dy) - star_bg;
                if (Val >= valmax * 0.5) {
                    pixel_counter += 1.0 / (kSamplePoints * kSamplePoints);
                }
                Val = Val / (kSamplePoints * kSamplePoints);
                const auto r = std::sqrt(dx * dx + dy * dy);
                if (r <= aperture_small) {
                    Sumval_small += Val;
                }
                SumVal  += Val;
                SumValR += Val * r;
            }
        }
        flux = std::max(Sumval_small, 0.00001);
        radius = aperture_small;
    }
    
    auto star_fwhm = 2 * std::sqrt(pixel_counter / kPi);
    auto hfd1 = 2 * SumValR / flux;
    hfd1 = std::max(0.7, hfd1);
    
    if (adu_e != 0) {
        flux  = flux  * adu_e * (xbinning * xbinning);
        sd_bg = sd_bg * adu_e *  xbinning;
    }
    
    auto snr = 0.0;
    if (flux >= 1) {
        snr = flux / std::sqrt(flux + radius * radius * kPi * sd_bg * sd_bg);
    } else {
        snr = 0;
    }
    
    result.hfd  = hfd1;
    result.fwhm = star_fwhm;
    result.snr  = snr;
    result.flux = flux;
    result.xc   = xc;
    result.yc   = yc;
    
    scratch.star_bg    = star_bg;
    scratch.sd_bg      = sd_bg;
    scratch.r_aperture = r_aperture;
}

/// MARK: - find_star_center

void find_star_center(const ImageArray& img,
                      int box, int x1, int y1,
                      int head_width, int head_height,
                      double& xc, double& yc) {
    // Local sub-pixel sampler
    auto value_subpixel = [&](double xs, double ys) -> double {
        const auto x_trunc = static_cast<int>(xs);
        const auto y_trunc = static_cast<int>(ys);
        if (x_trunc <= 0 || x_trunc >= (head_width - 2) ||
            y_trunc <= 0 || y_trunc >= (head_height - 2)) {
            return 0.0;
        }
        const auto x_frac = frac(xs);
        const auto y_frac = frac(ys);
        auto r = (img[0][y_trunc    ][x_trunc    ]) * (1 - x_frac) * (1 - y_frac);
        r     += (img[0][y_trunc    ][x_trunc + 1]) * (    x_frac) * (1 - y_frac);
        r     += (img[0][y_trunc + 1][x_trunc    ]) * (1 - x_frac) * (    y_frac);
        r     += (img[0][y_trunc + 1][x_trunc + 1]) * (    x_frac) * (    y_frac);
        return r;
    };
    
    const auto w = static_cast<int>(img[0][0].size());
    const auto h = static_cast<int>(img[0].size());
    
    // Bounds check
    if (!((x1 >= box) && (x1 < w - box) && (y1 >= box) && (y1 < h - box))) {
        xc = x1;
        yc = y1;
        return;
    }
    
    xc = x1;
    yc = y1;
    
    // Two iterations of weighted centroid refinement
    for (auto k = 1; k <= 2; ++k) {
        auto value = -99999.0;
        const auto xc_r = static_cast<int>(std::lround(xc));
        const auto yc_r = static_cast<int>(std::lround(yc));
        
        // Find peak value in the box
        for (auto i = xc_r - box; i <= xc_r + box; ++i) {
            for (auto j = yc_r - box; j <= yc_r + box; ++j) {
                const auto val = img[0][j][i];
                if (val > value) {
                    value = val;
                }
            }
        }
        
        // Weighted centroid
        auto SumVal = 0.0;
        auto SumValX = 0.0;
        auto SumValY = 0.0;
        for (auto i = -box; i <= box; ++i) {
            for (auto j = -box; j <= box; ++j) {
                auto val = value_subpixel(xc + i, yc + j) - value / 2;
                if (val > 0) {
                    val = val * val;
                }
                SumVal  += val;
                SumValX += val * i;
                SumValY += val * j;
            }
        }
        const auto Xg = SumValX / SumVal;
        const auto Yg = SumValY / SumVal;
        xc = xc + Xg;
        yc = yc + Yg;
    }
}

/// MARK: - measure_magnitudes

void measure_magnitudes(int annulus_rad,
                        int x1, int y1, int x2, int y2,
                        bool deep,
                        const ImageArray& img_loaded,
                        Header& head,
                        Background& bck,
                        double min_star_size_stacking,
                        StarList& stars) {
    // 5 rows, initial column allocation 5000
    stars.assign(5, std::vector<double>(5000, 0.0));
    
    // Local star-occupancy map
    auto img_sa = ImageArray(1,
        std::vector<std::vector<float>>(head.height,
            std::vector<float>(head.width, -1.0f)));
            
    get_background(0, img_loaded, /*calc_hist=*/false, /*calc_noise_level=*/true, bck);
                   
    auto detection_level = deep ? 5 * bck.noise_level : bck.star_level;
    auto hfd_min = std::max(0.8, min_star_size_stacking);
    
    auto nrstars = 0;
    
    auto saturation_level = (head.calstat.empty()) ? 64000.0f : 60000.0f;
    saturation_level = std::min(static_cast<float>(head.datamax_org - 1),
                                saturation_level);
                                
    // Initialise occupancy map
    for (auto fitsY = 0; fitsY < head.height; ++fitsY) {
        for (auto fitsX = 0; fitsX < head.width; ++fitsX) {
            img_sa[0][fitsY][fitsX] = -1.0f;
        }
    }
    
    auto scratch = HfdScratch{};
    auto res = HfdResult{};
    
    // Scan the region for stars
    for (auto fitsY = y1; fitsY < y2; ++fitsY) {
        for (auto fitsX = x1; fitsX < x2; ++fitsX) {
            if ((img_sa[0][fitsY][fitsX] <= 0) &&
                (img_loaded[0][fitsY][fitsX] - bck.backgr > detection_level)) {
                HFD(img_loaded, fitsX, fitsY, annulus_rad, head.mzero_radius,
                    /*adu_e*/0, head.xbinning, res, scratch);
                    
                if ((res.hfd < 15) && (res.hfd >= hfd_min) && (res.snr > 10)) {
                    // Mark the star's footprint in the occupancy map
                    auto radius = static_cast<int>(std::lround(3.0 * res.hfd));
                    auto sqr_radius = radius * radius;
                    auto xci = static_cast<int>(std::lround(res.xc));
                    auto yci = static_cast<int>(std::lround(res.yc));
                    for (auto n = -radius; n <= radius; ++n) {
                        for (auto m = -radius; m <= radius; ++m) {
                            auto j = n + yci;
                            auto i = m + xci;
                            if (j >= 0 && i >= 0 &&
                                j < head.height && i < head.width &&
                                (m * m + n * n) <= sqr_radius) {
                                img_sa[0][j][i] = 1.0f;
                            }
                        }
                    }
                    
                    // Check saturation in a 3x3 neighbourhood
                    const auto rxc = static_cast<int>(std::lround(res.xc));
                    const auto ryc = static_cast<int>(std::lround(res.yc));
                    if ((img_loaded[0][ryc    ][rxc    ] < saturation_level) &&
                        (img_loaded[0][ryc    ][rxc - 1] < saturation_level) &&
                        (img_loaded[0][ryc    ][rxc + 1] < saturation_level) &&
                        (img_loaded[0][ryc - 1][rxc    ] < saturation_level) &&
                        (img_loaded[0][ryc + 1][rxc    ] < saturation_level) &&
                        (img_loaded[0][ryc - 1][rxc - 1] < saturation_level) &&
                        (img_loaded[0][ryc + 1][rxc - 1] < saturation_level) &&
                        (img_loaded[0][ryc - 1][rxc + 1] < saturation_level) &&
                        (img_loaded[0][ryc + 1][rxc + 1] < saturation_level)) {
                        ++nrstars;
                        
                        // Grow the star list if needed
                        if (nrstars >= static_cast<int>(stars[0].size())) {
                            for (auto& row : stars) {
                                row.resize(nrstars + 5000, 0.0);
                            }
                        }
                        stars[0][nrstars - 1] = res.xc;
                        stars[1][nrstars - 1] = res.yc;
                        stars[2][nrstars - 1] = res.hfd;
                        stars[3][nrstars - 1] = res.flux;
                        stars[4][nrstars - 1] = res.snr;
                    }
                }
            }
        }
    }
    
    // Trim to actual count
    for (auto& row : stars) {
        row.resize(nrstars);
    }
}

/// MARK: - calibrate_photometry

void calibrate_photometry(const ImageArray& img,
                          std::vector<std::string>& memo,
                          Header& head,
                          [[maybe_unused]] Background& bck,
                          bool update,
                          double aperture_ratio_setting,
                          double annulus_radius_setting,
                          double& aperture_ratio,
                          int& annulus_radius,
                          std::string& passband_active) {
    if ((head.naxis == 0) || (head.cd1_1 == 0)) {
        return;
    }
    
    const auto apert = aperture_ratio_setting;
    
    if (update || (head.mzero == 0) || (aperture_ratio != apert) ||
        (passband_active != head.passband_database)) {
        memo.push_back("Photometric calibration of the measured stellar flux.");
        annulus_radius = 14;
        head.mzero_radius = 99;
        
        aperture_ratio = apert;
        
        if (apert != 0) {
            auto hfd_counter = 0;
            auto hfd_med = 0.0;
            auto bck_local = Background{};
            analyse_image(img, head, 30, 0, hfd_counter, bck_local, hfd_med);
            
            if (hfd_med != 0) {
                char buf[128];
                std::snprintf(buf, sizeof(buf),
                              "Median HFD is %g. Aperture and annulus will be adapted accordingly.",
                              hfd_med);
                memo.push_back(buf);
                head.mzero_radius = hfd_med * apert / 2;
                const auto annul = annulus_radius_setting;
                annulus_radius = std::min(50,
                    static_cast<int>(std::lround(hfd_med * annul / 2)) - 1);
            }
        } else {
            memo.push_back("To increase the accuracy of point sources magnitudes "
                           "set a smaller aperture diameter in tab \"photometry\".");
        }
        
        plot_and_measure_stars(img, memo, head, /*calibration*/true,
                               /*plot_stars*/false,
                               /*report_lim_magnitude*/true);
    }
}

/// MARK: - test_star_spectrum

float test_star_spectrum(float r, float g, float b) noexcept {
    if ((b < (0x12 / 255.0f) * r) || (b > (255.0f / 0xAD) * r)) {
        return 1;
    }
    
    if ((r <= 1) || (g <= 1) || (b <= 1)) {
        return 0;
    }
    
    const auto RdivG = r / g;
    return std::abs((b / g) - (0.6427f * RdivG * RdivG - 2.868f * RdivG + 3.3035f));
}

/// MARK: - measure_hotpixels

void measure_hotpixels(int x1, int y1, int x2, int y2, int col,
                       double sd, double mean,
                       const ImageArray& img,
                       double& hotpixel_perc, double& hotpixel_adu) {
    const auto w = static_cast<int>(img[0][0].size());
    const auto h = static_cast<int>(img[0].size());
    
    x1 = std::max(x1, 0);
    x2 = std::min(x2, w - 1);
    y1 = std::max(y1, 0);
    y2 = std::min(y2, h - 1);
    
    if ((y1 > y2) || (x1 > x2)) {
        return;
    }
    
    hotpixel_adu = 0;
    auto counter = 0;
    auto counter2 = 0;
    for (auto j = y1; j <= y2; ++j) {
        for (auto i = x1; i <= x2; ++i) {
            const auto value = std::abs(img[col][j][i] - mean);
            if (value <= 3 * sd) {
                ++counter;
            } else {
                hotpixel_adu += value * value;
                ++counter2;
            }
        }
    }
    
    if (counter2 > 0) {
        hotpixel_adu = std::sqrt(hotpixel_adu / counter2);
    } else {
        hotpixel_adu = 0;
    }
    
    hotpixel_perc = 0.997 - static_cast<double>(counter) / (counter + counter2);
}

/// MARK: - local_sd

void local_sd(int x1, int y1, int x2, int y2, int col,
              const ImageArray& img,
              double& sd, double& mean, int& iterations) {
    const auto w = static_cast<int>(img[0][0].size());
    const auto h = static_cast<int>(img[0].size());
    
    x1 = std::max(x1, 0);
    x2 = std::min(x2, w - 1);
    y1 = std::max(y1, 0);
    y2 = std::min(y2, h - 1);
    
    sd = 99999;
    mean = 0;
    if ((y1 > y2) || (x1 > x2)) {
        return;
    }
    
    iterations = 0;
    auto sd_old = 0.0;
    do {
        // Compute mean
        auto counter = 0;
        auto meanx = 0.0;
        for (auto j = y1; j <= y2; ++j) {
            for (auto i = x1; i <= x2; ++i) {
                const auto value = img[col][j][i];
                if ((iterations == 0) || (std::abs(value - mean) <= 3 * sd)) {
                    ++counter;
                    meanx += value;
                }
            }
        }
        if (counter != 0) {
            mean = meanx / counter;
        }
        
        // Compute standard deviation
        sd_old = sd;
        counter = 0;
        for (auto j = y1; j <= y2; ++j) {
            for (auto i = x1; i <= x2; ++i) {
                const auto value = img[col][j][i] - mean;
                if ((value < mean) &&
                    ((iterations == 0) || (std::abs(value) <= 3 * sd_old))) {
                    sd = sd + value * value;
                    ++counter;
                }
            }
        }
        if (counter != 0) {
            sd = std::sqrt(sd / counter);
        }
        ++iterations;
    } while (!(((sd_old - sd) < 0.03 * sd) || (iterations >= 7)));
}

/// MARK: - mode

int mode(const ImageArray& img, bool ellipse,
         int colorm, int xmin, int xmax, int ymin, int ymax,
         int max1, int& greylevels) {
    const auto height3 = static_cast<int>(img[0].size());
    const auto width3  = static_cast<int>(img[0][0].size());
    
    max1 -= 10;
    if (xmin < 0) {
        xmin = 0;
    }
    if (xmax > width3 - 1) {
        xmax = width3 - 1;
    }
    if (ymin < 0) {
        ymin = 0;
    }
    if (ymax > height3 - 1) {
        ymax = height3 - 1;
    }
    
    auto histogram = std::vector<int>(max1 + 1, 0);
    
    const auto centerX = (xmax + xmin) / 2.0;
    const auto centerY = (ymax + ymin) / 2.0;
    const auto a = (xmax - xmin - 1) / 2.0;
    const auto b = (ymax - ymin - 1) / 2.0;
    
    for (auto i = ymin; i <= ymax; ++i) {
        for (auto j = xmin; j <= xmax; ++j) {
            if (!ellipse ||
                ((j - centerX) * (j - centerX) / (a * a) +
                 (i - centerY) * (i - centerY) / (b * b) < 1)) {
                const auto val = static_cast<int>(std::lround(img[colorm][i][j]));
                if ((val >= 1) && (val < max1)) {
                    histogram[val] += 1;
                }
            }
        }
    }
    
    auto result = 0;
    greylevels = 0;
    auto value_count = 0;
    for (auto i = 1; i <= max1; ++i) {
        const auto val = histogram[i];
        if (val != 0) {
            ++greylevels;
        }
        if (val > value_count) {
            value_count = val;
            result = i;
        }
    }
    
    return result;
}

/// MARK: - get_negative_noise_level

double get_negative_noise_level(const ImageArray& img,
                                int colorm,
                                int xmin, int xmax, int ymin, int ymax,
                                double common_level,
                                int head_width, int head_height) {
    if (xmin < 0) {
        xmin = 0;
    }
    if (xmax > head_width  - 1) {
        xmax = head_width  - 1;
    }
    if (ymin < 0) {
        ymin = 0;
    }
    if (ymax > head_height - 1) {
        ymax = head_height - 1;
    }
    
    auto result = 0.0;
    auto count_neg = 0;
    for (auto i = ymin; i <= ymax; ++i) {
        for (auto j = xmin; j <= xmax; ++j) {
            const auto col = static_cast<int>(std::lround(img[colorm][i][j]));
            if ((col >= 1) && (col <= common_level)) {
                ++count_neg;
                result += (col - common_level) * (col - common_level);
            }
        }
    }
    
    if (count_neg >= 1) {
        return std::sqrt(result / count_neg);
    }
    
    return 0;
}

/// MARK: - get_background

void get_background(int colour,
                    const ImageArray& img,
                    bool calc_hist,
                    bool calc_noise_level,
                    Background& back) {
    if (calc_hist) {
        // Populate the module-level histogram + his_mean for this channel.
        get_hist(colour, img);
    }

    back.backgr = img[0][0][0];

    // Find the most common pixel value (mode)
    auto pixels = std::uint32_t{0};
    auto max_range = his_mean[colour];
    for (auto i = 1; i <= max_range; ++i) {
        if (histogram[colour][i] > pixels) {
            pixels = histogram[colour][i];
            back.backgr = i;
        }
    }

    // Fall back to mean if it is much higher than the mode
    if (his_mean[colour] > 1.5 * back.backgr) {
        memo1_lines.push_back(std::format(
            "{}, will use mean value {} as background rather then most common value {}",
            filename2, his_mean[colour],
            static_cast<int>(std::lround(back.backgr))));
        back.backgr = his_mean[colour];
    }

    if (calc_noise_level) {
        const auto width5  = static_cast<int>(img[0][0].size());
        const auto height5 = static_cast<int>(img[0].size());
        auto stepsize = static_cast<int>(std::lround(height5 / 71.0));
        if ((stepsize % 2) == 0) {
            stepsize += 1;
        }

        // Iterative sigma-clipped noise estimation
        auto sd = 99999.0;
        auto iterations = 0;
        auto sd_old = 0.0;
        do {
            auto fitsX = 15;
            auto counter = 1;
            sd_old = sd;
            while (fitsX <= width5 - 1 - 15) {
                auto fitsY = 15;
                while (fitsY <= height5 - 1 - 15) {
                    const auto value = img[colour][fitsY][fitsX];
                    if ((value < back.backgr * 2) && (value != 0)) {
                        if ((iterations == 0) ||
                            (std::abs(value - back.backgr) <= 3 * sd_old)) {
                            sd = sd + (value - back.backgr) * (value - back.backgr);
                            ++counter;
                        }
                    }
                    fitsY += stepsize;
                }
                fitsX += stepsize;
            }
            sd = std::sqrt(sd / counter);
            ++iterations;
        } while (!(((sd_old - sd) < 0.05 * sd) || (iterations >= 7)));
        back.noise_level = sd;

        // Compute star detection thresholds from the histogram
        if ((nrbits == 8) || (nrbits == 24)) {
            max_range = 255;
        } else {
            max_range = 65001;
        }

        back.star_level = 0;
        back.star_level2 = 0;
        auto i = max_range;
        const auto factor  =  6LL * max_stars_setting;
        const auto factor2 = 24LL * max_stars_setting;
        auto above = 0LL;
        while ((back.star_level == 0) && (i > back.backgr + 1)) {
            --i;
            above += histogram[colour][i];
            if (above >= factor) {
                back.star_level = i;
            }
        }
        while ((back.star_level2 == 0) && (i > back.backgr + 1)) {
            --i;
            above += histogram[colour][i];
            if (above >= factor2) {
                back.star_level2 = i;
            }
        }
        back.star_level  = std::max(std::max(3.5 * sd, 1.0),
                                    back.star_level  - back.backgr - 1);
        back.star_level2 = std::max(std::max(3.5 * sd, 1.0),
                                    back.star_level2 - back.backgr - 1);
    }
}

/// MARK: - retrieve_ADU_to_e_unbinned

double retrieve_ADU_to_e_unbinned(const std::string& head_egain,
                                  double egain_extra_factor,
                                  double egain_default,
                                  bool noise_in_electron) {
    if ((egain_extra_factor != 0) && noise_in_electron) {
        auto egain = 0.0;
        try {
            egain = std::stod(head_egain);
        } catch (...) {
            egain = 0;
        }
        if (egain == 0) {
            egain = egain_default;
        }
        return egain / egain_extra_factor;
    }
    
    return 0;
}

/// MARK: - noise_to_electrons

std::string noise_to_electrons(double adu_e, double sd) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.1f", sd);
    auto result = std::string(buf);
    if (adu_e != 0) {
        result += " e-";
    }
    return result;
}
 
} // namespace
