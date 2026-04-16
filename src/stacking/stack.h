///----------------------------------------
///      @file stack.h
///   @ingroup ASTAP++
///     @brief Non-GUI algorithmic exports from the image stacking module.
///   @details Ports the non-GUI procedures and functions originally declared
///            in unit_stack.pas plus a handful of file-local algorithmic
///            helpers. GUI event handlers, list-view updaters and memo/tab
///            side effects are deliberately skipped.
///    @author Ported from Han Kleijn's unit_stack.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "../types.h"
#include "../core/globals.h"
#include "../solving/astrometric_solving.h"

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

///----------------------------------------
/// MARK: Local aliases
///----------------------------------------

/// [channel][row][column] float image buffer.
using astap::ImageArray;

/// FITS/WCS header struct.
using astap::Header;

/// File-to-process record.
using astap::FileToDo;

/// Background/noise measurement record.
using astap::Background;

///----------------------------------------
/// MARK: Non-GUI algorithmic exports
///----------------------------------------

///----------------------------------------
///   @brief Blur by averaging values of neighbouring pixels, skipping zero
///          samples (dead pixels or unfilled stack regions).
/// @details @p range selects the kernel footprint: 2 -> 2x2, 3 -> 3x3,
///          4 -> 4x4, N -> NxN centred on the current pixel. @p colors is
///          the requested output channel count (1 collapses RGB to a mono
///          average).
///   @param colors Output channel count.
///   @param range Kernel footprint selector.
/// @param[in,out] img Image buffer to blur in place.
///----------------------------------------

void box_blur(int colors, int range, ImageArray& img);

///----------------------------------------
///  @brief Normalise the four cells of a Bayer 2x2 tile so a raw OSC flat
///         does not introduce colour shifts when a non-white flat source
///         was used.
/// @details Mono/colour (naxis3 > 1) inputs are skipped.
/// @param[in,out] img Image buffer to normalise in place.
///----------------------------------------

void check_pattern_filter(ImageArray& img);

///----------------------------------------
///  @brief Replace interior zero pixels (dead pixels / saturation markers)
///         with the mean of their non-zero neighbours.
/// @details Uses an expanding square search up to radius 10. Leaves the
///          outer black borders (detected per channel) alone.
/// @param[in,out] img Image buffer to filter in place.
///----------------------------------------

void black_spot_filter(ImageArray& img);

///----------------------------------------
///  @brief Run the internal plate solver on an already-loaded image, then
///         save an updated FITS/TIFF alongside the source.
///  @param img Image to solve.
/// @param[in,out] hd Header updated with solution metadata.
/// @param[in,out] memo Log lines.
/// @return @c true iff both solving and saving succeeded.
///----------------------------------------

[[nodiscard]] bool update_solution_and_save(const ImageArray& img,
                                            Header& hd,
                                            std::vector<std::string>& memo);
                                            
///----------------------------------------
///  @brief Apply master dark and master flat to @p img if a matching
///         calibration set exists and has not already been applied.
/// @details Updates @c hd.calstat, counts, pedestal, etc.
/// @param[in,out] img Image buffer to calibrate in place.
/// @param[in,out] hd Header updated with calibration metadata.
/// @return @c true iff any calibration was applied.
///----------------------------------------

[[nodiscard]] bool apply_dark_and_flat(ImageArray& img, Header& hd);

///----------------------------------------
///   @brief Bright-star colour smoothing: combines colour ratios over a
///          @p wide x @p wide aperture while preserving luminance.
/// @param[in,out] img Image buffer to smooth in place.
///    @param wide Aperture size in pixels.
///    @param sd Detection sigma above noise.
///    @param preserve_r_nebula Skip pixels dominated by red nebulosity.
///    @param measurehist Force a histogram re-measurement.
///----------------------------------------

void smart_colour_smooth(ImageArray& img,
                         double wide,
                         double sd,
                         bool preserve_r_nebula,
                         bool measurehist);
                         
///----------------------------------------
///  @brief Rebalance RGB ratios so pixels where green dominates both
///         others (or where green is weaker than both) are pushed toward
///         a natural colour.
/// @details Used to remove green/purple casts from Hubble-palette data.
/// @param[in,out] img Image buffer to filter in place.
///----------------------------------------

void green_purple_filter(ImageArray& img);

///----------------------------------------
///  @brief Parse FITS DATE-OBS / DATE-AVG strings and the exposure time
///         into the module-level globals @c jd_start, @c jd_end and
///         @c jd_mid.
///  @param date_obs FITS DATE-OBS value.
///  @param date_avg FITS DATE-AVG value (may be empty).
///  @param exp Exposure time in seconds.
///----------------------------------------

void date_to_jd(std::string_view date_obs,
                std::string_view date_avg,
                double exp);
                
///----------------------------------------
///  @brief Convert a Julian Date back to an ISO-like date string
///         (@c YYYY-MM-DDThh:mm:ss).
/// @details Meeus, @e Astronomical @e Algorithms, chapter 7.
///  @param jd Julian Date.
/// @return ISO-like date string, or a diagnostic if @p jd is out of range.
///----------------------------------------

[[nodiscard]] std::string jd_to_date(double jd);

///----------------------------------------
///  @brief Free-ratio bicubic resize of the globally loaded image
///         @c img_loaded and matching WCS-related header fields.
///  @param ratio Scale factor applied to width and height.
///----------------------------------------

void resize_img_loaded(double ratio);

///----------------------------------------
///   @brief Median ADU value of a @p sizeX x @p sizeY window centred at
///          (@p x, @p y) in channel @p color.
/// @details Window sizes are rounded up to an odd number. Zero-valued
///          pixels are skipped. Used by several background-equalisation
///          tools.
///   @param img Source image.
///   @param color Channel index.
///   @param sizeX Window width.
///   @param sizeY Window height.
///   @param x Window centre column.
///   @param y Window centre row.
///  @return Median of the non-zero samples in the window, or 0.0 if none.
///----------------------------------------

[[nodiscard]] double median_background(const ImageArray& img,
                                       int color,
                                       int sizeX,
                                       int sizeY,
                                       int x,
                                       int y);
                                       
///----------------------------------------
///      @brief Scan @p img for stars above the noise floor, returning the
///             total star count, the measured background block and the
///             median HFD.
///    @details @p report_type: 0 = report only, 1 = report + write .csv,
///             2 = .csv only.
///      @param img Source image.
///      @param head Image header (for WCS lookups).
///      @param snr_min Minimum signal-to-noise threshold.
///      @param report_type Reporting mode selector.
/// @param[out] star_counter Number of accepted stars.
/// @param[out] bck Measured background block.
/// @param[out] hfd_median Median half-flux diameter.
///----------------------------------------

void analyse_image(const ImageArray& img,
                   const Header& head,
                   double snr_min,
                   int report_type,
                   int& star_counter,
                   Background& bck,
                   double& hfd_median);
                   
///----------------------------------------
///   @brief Apply a "most common value" (mode) filter.
/// @details Partitions the source image into @p radius * 2 tiles and
///          replaces every pixel in a tile with the tile's modal value.
///   @param sourc Source image.
/// @param[out] dest Destination image.
///   @param datamax Maximum expected data value.
///   @param radius Tile radius.
///----------------------------------------

void apply_most_common(const ImageArray& sourc,
                       ImageArray& dest,
                       double datamax,
                       int radius);
                       
///----------------------------------------
///  @brief Julian day for the given calendar date/time.
/// @details Gregorian switchover at 1582-10-15. Meeus 7.1.
///  @param yyyy Calendar year.
///  @param mm Calendar month (1-12).
///  @param dd Day of month (fractional).
///  @param hours Hour component.
///  @param minutes Minute component.
///  @param seconds Second component.
/// @return Julian Date.
///----------------------------------------

[[nodiscard]] double julian_calc(int yyyy,
                                 int mm,
                                 double dd,
                                 double hours,
                                 double minutes,
                                 double seconds);
                                 
///----------------------------------------
///  @brief Strip characters that are illegal in a filename on
///         Windows/macOS: @c { '.', '\\', '/', '*', '"', ':', '|', '<', '>' }.
///  @param s Input string.
/// @return Input with illegal characters removed.
///----------------------------------------

[[nodiscard]] std::string remove_special_chars(std::string_view s);

/// @note @c jd_start / @c jd_end / @c jd_mid are written by @ref date_to_jd
///       and consumed by the stacking/photometry pipeline. Canonical
///       declarations live in @c ../core/globals.h.
                
} // namespace
