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

#include <functional>
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
    /// @brief Count of pixels > 3σ above local background inside the
    ///        centroid-refinement box. Exposed so callers can reject
    ///        extended sources (galaxies / nebulae) whose fill fraction is
    ///        much higher than a point source's.
    int    signal_counter{0};
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
                        
/// MARK: - Photometric flux calibration

///----------------------------------------
/// @brief One catalog star yielded by a @ref StarSource.
/// @details Units match the Pascal convention used by the .1476/.290 readers:
///          magnitude is stored as @c magn = mag * 10, and Gaia Bp-Rp is
///          @c Bp_Rp = (Bp - Rp) * 10. Sentinels are 999 (no colour recorded)
///          and -128 (unreliable Johnson-V, skip).
///----------------------------------------

struct CatalogStar {
    double ra{};      ///< @brief Right ascension in radians.
    double dec{};     ///< @brief Declination in radians.
    double magn{};    ///< @brief Magnitude * 10.
    double Bp_Rp{};   ///< @brief Gaia Bp-Rp * 10, or 999 (mono database) / -128 (unreliable).
};

///----------------------------------------
/// @brief Catalog-star iterator protocol.
/// @details Callers pass a source that yields one star per invocation, writing
///          into @p out and returning @c true. When the catalog is exhausted
///          the source returns @c false. Three concrete sources exist in the
///          wider codebase: local .1476/.290 files, the wide-field @c w08
///          buffer, and the online Gaia query cache.
///----------------------------------------

using StarSource = std::function<bool(CatalogStar& out)>;

///----------------------------------------
/// @brief Result of a call to @ref calibrate_flux.
///----------------------------------------

struct FluxCalibrationResult {
    bool        success{false};                ///< @brief True if @p head.mzero was set.
    int         stars_measured{0};             ///< @brief Count of stars that survived all filters.
    double      standard_error_mean{0.0};      ///< @brief SEM of the log-ratio mean, in magnitudes.
    std::string message;                       ///< @brief Human-readable summary of the result.
};

///----------------------------------------
/// @brief Photometric flux calibration: measure catalog stars and compute MZERO.
/// @details C++23 port of the flux-calibration half of
///          @c plot_and_measure_stars in @c unit_annotation.pas, refactored so
///          the catalog iterator is injectable. For each star yielded by
///          @p next_star the routine projects the RA/Dec onto the image,
///          runs @ref HFD at the predicted pixel, rejects saturated and
///          low-SNR candidates, and accumulates @c flux_ratio =
///          @c flux / 10^((21 - mag) / 2.5). After the loop, a MAD-based
///          robust mean of the flux ratios yields
///          @c head.mzero = 21 + 2.5 log10(mean_flux_ratio). At least three
///          usable stars are required for success.
///
///          When @p report_lim_magn is @c true the routine also writes
///          @c head.magn_limit from the per-star HFD * SD products.
///
///          This routine updates the following @p head fields on success:
///          @c mzero, @c passband_database, and (if @p report_lim_magn)
///          @c magn_limit. It does @em not touch @c mzero_radius — the caller
///          must set that to either the point-source aperture radius or 99
///          (extended objects) before calling.
///
///          No GUI drawing is performed; the annotation-overlay half of the
///          Pascal original lives in a separate routine (future work).
///    @param img             Read-only image data (single channel used).
///    @param[in,out] head    FITS header (updated on success).
///    @param[in,out] memo    Log sink for progress messages and MZERO card updates.
///    @param next_star       Injectable catalog iterator.
///    @param passband_active Passband identifier written to @c head.passband_database.
///    @param annulus_radius  HFD annulus radius in pixels (typically 14).
///    @param aperture_setting Aperture diameter setting in HFD multiples; 0 means
///                           "max" (extended objects), which forces a large
///                           virtual aperture for the limiting-magnitude calc.
///    @param report_lim_magn When @c true, compute and store @c head.magn_limit.
/// @return @ref FluxCalibrationResult carrying success flag, star count, SEM,
///         and a summary message.
///----------------------------------------

[[nodiscard]] FluxCalibrationResult calibrate_flux(
    const ImageArray& img,
    Header& head,
    std::vector<std::string>& memo,
    const StarSource& next_star,
    std::string_view passband_active,
    int annulus_radius,
    double aperture_setting,
    bool report_lim_magn);

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
/// @details Uses the module-level @c histogram, @c his_mean, @c nrbits,
///          @c max_stars_setting, @c filename2 and @c memo1_lines globals.
/// @param colour            Colour channel index.
/// @param img               Read-only image data.
/// @param calc_hist         Whether to recompute the histogram first.
/// @param calc_noise_level  Whether to compute noise and star levels.
/// @param[out] back         Background statistics output.
///----------------------------------------

void get_background(int colour,
                    const ImageArray& img,
                    bool calc_hist,
                    bool calc_noise_level,
                    Background& back);
                    
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
