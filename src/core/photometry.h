///----------------------------------------
///      @file photometry.h
///   @ingroup ASTAP++
///     @brief Photometry and star-measurement primitives.
///   @details HFD/FWHM/SNR/flux measurement, centroid refinement, background
///            statistics, and helpers used by the photometric calibration
///            pipeline. Exposes routines as free functions with previously-global
///            scratch values threaded through reference parameters.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP), MPL-2.0.
///            C++ port by John Stephen.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <string>
#include <vector>

#include "../types.h"

///----------------------------------------
namespace astap::core {
///----------------------------------------

/// MARK: - HFD Types

///----------------------------------------
/// @brief Aggregate result from the HFD measurement procedure.
/// @details Bundles the six output values that the original returned as
///          separate out-parameters.
///----------------------------------------

struct HfdResult {
    double hfd{999.0};
    double fwhm{0.0};
    double snr{0.0};
    double flux{0.0};
    double xc{0.0};
    double yc{0.0};
};

///----------------------------------------
/// @brief Scratch state shared between HFD and downstream callers.
/// @details Mirrors the unit-level globals star_bg, sd_bg, and r_aperture
///          from the original source.
///----------------------------------------

struct HfdScratch {
    double star_bg{0.0};
    double sd_bg{0.0};
    int    r_aperture{-1};
};

/// MARK: - Core Measurement

///----------------------------------------
/// @brief Measure the half-flux diameter, FWHM, SNR, and flux of a star.
/// @param img            Read-only image data.
/// @param x1             Initial X centre estimate.
/// @param y1             Initial Y centre estimate.
/// @param rs             Search radius (reduced internally if not boxed).
/// @param aperture_small Small-aperture threshold.
/// @param adu_e          ADU-to-electron conversion factor (0 to skip).
/// @param xbinning       X binning factor.
/// @param[out] result    Measured HFD, FWHM, SNR, flux, and refined centre.
/// @param[in,out] scratch  Shared star_bg / sd_bg / r_aperture state.
///----------------------------------------

void HFD(const ImageArray& img,
         int x1, int y1, int rs,
         double aperture_small, double adu_e,
         double xbinning,
         HfdResult& result,
         HfdScratch& scratch);
         
///----------------------------------------
/// @brief Refine a star centre position using iterative centroiding.
/// @param img         Read-only image data.
/// @param box         Half-size of the search box.
/// @param x1          Initial X position.
/// @param y1          Initial Y position.
/// @param head_width  Image width from header.
/// @param head_height Image height from header.
/// @param[out] xc     Refined X centre.
/// @param[out] yc     Refined Y centre.
///----------------------------------------

void find_star_center(const ImageArray& img,
                      int box, int x1, int y1,
                      int head_width, int head_height,
                      double& xc, double& yc);
                      
///----------------------------------------
/// @brief Detect and measure stars in a sub-region for magnitude calibration.
/// @param annulus_rad           Annulus radius for HFD measurement.
/// @param x1                   Left bound of the search region.
/// @param y1                   Top bound of the search region.
/// @param x2                   Right bound of the search region.
/// @param y2                   Bottom bound of the search region.
/// @param deep                 Use deeper detection threshold (5x noise).
/// @param img_loaded           Read-only loaded image.
/// @param head                 FITS header.
/// @param bck                  Background statistics (updated).
/// @param min_star_size_stacking Minimum HFD for valid stars.
/// @param[out] stars           Detected star list (5 rows: x, y, hfd, flux, snr).
///----------------------------------------

void measure_magnitudes(int annulus_rad,
                        int x1, int y1, int x2, int y2,
                        bool deep,
                        const ImageArray& img_loaded,
                        Header& head,
                        Background& bck,
                        double min_star_size_stacking,
                        StarList& stars);
                        
///----------------------------------------
/// @brief Run photometric calibration on the image.
/// @param img                    Read-only image data.
/// @param memo                   Log message output.
/// @param head                   FITS header (updated with mzero etc.).
/// @param bck                    Background statistics.
/// @param update                 Force recalibration.
/// @param aperture_ratio_setting Aperture ratio from settings (0 = max).
/// @param annulus_radius_setting Annulus radius from settings.
/// @param[in,out] aperture_ratio Current aperture ratio state.
/// @param[out] annulus_radius    Computed annulus radius.
/// @param[in,out] passband_active Active passband string.
///----------------------------------------

void calibrate_photometry(const ImageArray& img,
                          std::vector<std::string>& memo,
                          Header& head,
                          Background& bck,
                          bool update,
                          double aperture_ratio_setting,
                          double annulus_radius_setting,
                          double& aperture_ratio,
                          int& annulus_radius,
                          std::string& passband_active);
                          
///----------------------------------------
/// @brief Test how well an RGB triplet matches a stellar blackbody spectrum.
/// @param r Red channel value.
/// @param g Green channel value.
/// @param b Blue channel value.
/// @return Deviation from the expected spectrum (0 = perfect match, 1 = out of range).
///----------------------------------------

[[nodiscard]] float test_star_spectrum(float r, float g, float b) noexcept;

///----------------------------------------
/// @brief Measure hot-pixel statistics in a sub-region.
/// @param x1   Left bound.
/// @param y1   Top bound.
/// @param x2   Right bound.
/// @param y2   Bottom bound.
/// @param col  Colour channel index.
/// @param sd   Standard deviation of the background.
/// @param mean Background mean.
/// @param img  Read-only image data.
/// @param[out] hotpixel_perc Fraction of hot pixels above 3-sigma.
/// @param[out] hotpixel_adu  RMS ADU of hot-pixel deviations.
///----------------------------------------

void measure_hotpixels(int x1, int y1, int x2, int y2, int col,
                       double sd, double mean,
                       const ImageArray& img,
                       double& hotpixel_perc, double& hotpixel_adu);
                       
///----------------------------------------
/// @brief Compute iterative sigma-clipped local standard deviation and mean.
/// @param x1   Left bound.
/// @param y1   Top bound.
/// @param x2   Right bound.
/// @param y2   Bottom bound.
/// @param col  Colour channel index.
/// @param img  Read-only image data.
/// @param[out] sd         Computed standard deviation.
/// @param[out] mean       Computed mean.
/// @param[out] iterations Number of sigma-clip iterations performed.
///----------------------------------------

void local_sd(int x1, int y1, int x2, int y2, int col,
              const ImageArray& img,
              double& sd, double& mean, int& iterations);
              
///----------------------------------------
/// @brief Compute the histogram mode (most common pixel value) of a sub-region.
/// @param img        Read-only image data.
/// @param ellipse    If true, restrict sampling to an elliptical region.
/// @param colorm     Colour channel index.
/// @param xmin       Left bound.
/// @param xmax       Right bound.
/// @param ymin       Top bound.
/// @param ymax       Bottom bound.
/// @param max1       Upper histogram bin limit.
/// @param[out] greylevels Number of distinct non-zero grey levels found.
/// @return The most frequent pixel value (the mode).
///----------------------------------------

[[nodiscard]] int mode(const ImageArray& img, bool ellipse,
                       int colorm, int xmin, int xmax, int ymin, int ymax,
                       int max1, int& greylevels);
                       
///----------------------------------------
/// @brief Estimate noise from the negative side of the pixel distribution.
/// @param img          Read-only image data.
/// @param colorm       Colour channel index.
/// @param xmin         Left bound.
/// @param xmax         Right bound.
/// @param ymin         Top bound.
/// @param ymax         Bottom bound.
/// @param common_level Mode or background level.
/// @param head_width   Image width.
/// @param head_height  Image height.
/// @return RMS of negative-side deviations, or 0 if no qualifying pixels.
///----------------------------------------

[[nodiscard]] double get_negative_noise_level(const ImageArray& img,
                                              int colorm,
                                              int xmin, int xmax, int ymin, int ymax,
                                              double common_level,
                                              int head_width, int head_height);
                                              
///----------------------------------------
/// @brief Compute background level, noise, and star-detection thresholds.
/// @param colour            Colour channel index.
/// @param img               Read-only image data.
/// @param calc_hist         Whether to recompute the histogram first.
/// @param calc_noise_level  Whether to compute noise and star levels.
/// @param[out] back         Background statistics output.
/// @param histogram         Shared 3-channel histogram array.
/// @param his_mean          Shared per-channel mean values.
/// @param nrbits            Bit depth of the image (8, 16, 24).
/// @param max_stars_setting Maximum stars tuning knob.
/// @param filename2         Filename for log messages.
/// @param memo              Log message output.
///----------------------------------------

void get_background(int colour,
                    const ImageArray& img,
                    bool calc_hist,
                    bool calc_noise_level,
                    Background& back,
                    int (&histogram)[3][65536],
                    int (&his_mean)[3],
                    int nrbits,
                    int max_stars_setting,
                    const std::string& filename2,
                    std::vector<std::string>& memo);
                    
///----------------------------------------
/// @brief Retrieve the ADU-to-electron conversion factor from FITS metadata.
/// @param head_egain          EGAIN string from the FITS header.
/// @param egain_extra_factor  Extra gain factor from settings.
/// @param egain_default       Default gain if parsing fails.
/// @param noise_in_electron   Whether noise units are in electrons.
/// @return Conversion factor, or 0 if not applicable.
///----------------------------------------

[[nodiscard]] double retrieve_ADU_to_e_unbinned(const std::string& head_egain,
                                                double egain_extra_factor,
                                                double egain_default,
                                                bool noise_in_electron);
                                                
///----------------------------------------
/// @brief Format a noise value with optional electron units.
/// @param adu_e ADU-to-electron factor (0 = ADU units only).
/// @param sd    Standard deviation value.
/// @return Formatted string like "12.3" or "12.3 e-".
///----------------------------------------

[[nodiscard]] std::string noise_to_electrons(double adu_e, double sd);
    
} // namespace
