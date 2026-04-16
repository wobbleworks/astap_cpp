///----------------------------------------
///      @file imaging.h
///   @ingroup ASTAP++
///     @brief Image-manipulation utilities: binning, mono conversion, rotation,
///            duplication, coordinate flipping, Bayer extraction, histogram,
///            and luminance/saturation stretch.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP, MPL-2.0)
///            by John Stephen.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "../types.h"

///----------------------------------------
namespace astap::core {
///----------------------------------------

using astap::Background;
using astap::Header;
using astap::ImageArray;

/// MARK: Module-level histogram state

/// @brief Histogram bins: three channels, 0..65535 inclusive.
extern std::array<std::array<std::uint32_t, 65536>, 3> histogram;

/// @brief Total number of red samples inside the histogram window.
extern std::uint32_t his_total_red;

/// @brief Mean intensity per channel measured during @c get_hist.
extern std::array<int, 3> his_mean;

/// @brief Stretch lookup table indexed 0..32768; values in 0..1.
extern std::array<float, 32769> stretch_c;

/// @brief Histogram x-axis range set by @c use_histogram for the viewer widget.
extern int hist_range;

/// @brief Current white point used by @c stretch_image.
extern double cwhite;

/// @brief Whether @c stretch_image should apply the stretch_c LUT.
extern bool stretch_on;

/// @brief Saturation factor used by @c stretch_image.
extern float saturation_factor;

/// MARK: Image manipulation

///----------------------------------------
/// @brief Bin an image in-place by 2x2, 3x3, or 4x4.
/// @details Updates header NAXIS1/2, CRPIXn, CDELTn, the CD matrix,
///          X/YBINNING, X/YPIXSZ and appends HISTORY entries via memo.
/// @param img        Image to bin (modified in-place).
/// @param head       FITS header (modified in-place).
/// @param memo       Header memo lines (receives HISTORY updates).
/// @param binfactor  Bin factor, clamped to [2, 4].
/// @param filename2  Original filename for HISTORY annotation.
///----------------------------------------

void bin_X2X3X4(ImageArray& img,
                Header& head,
                std::vector<std::string>& memo,
                int binfactor,
                const std::string& filename2);
                
///----------------------------------------
/// @brief Average three colour channels into a single mono channel.
/// @details No-op if the image is already mono. Updates head.naxis / head.naxis3.
/// @param img  Image to convert (modified in-place).
/// @param head FITS header (modified in-place).
///----------------------------------------

void convert_mono(ImageArray& img, Header& head);

///----------------------------------------
/// @brief Rotate an image by an arbitrary angle CCW, expanding the canvas.
/// @details Updates head WCS (CRPIX/CRotA/CD-matrix) accordingly.
/// @param angle         Rotation angle in degrees CCW.
/// @param flipped_view  View flip factor (+1 or -1).
/// @param flipped_image Image flip factor (+1 or -1).
/// @param img           Image to rotate (modified in-place).
/// @param head          FITS header (modified in-place).
/// @param memo          Header memo lines (receives HISTORY updates).
///----------------------------------------

void rotate_arbitrary(double angle,
                      double flipped_view,
                      double flipped_image,
                      ImageArray& img,
                      Header& head,
                      std::vector<std::string>& memo);
                      
///----------------------------------------
/// @brief Deep-copy an image.
/// @param img Source image.
/// @return Independent copy of the image.
///----------------------------------------

[[nodiscard]] ImageArray duplicate(const ImageArray& img);

///----------------------------------------
/// @brief Flip array/screen coordinates respecting horizontal/vertical flags.
/// @param x1              Input X coordinate.
/// @param y1              Input Y coordinate.
/// @param[out] x2         Output X coordinate.
/// @param[out] y2         Output Y coordinate.
/// @param head            FITS header (provides image dimensions).
/// @param flip_horizontal Whether to flip horizontally.
/// @param flip_vertical   Whether to flip vertically.
///----------------------------------------

void flip(int x1,
          int y1,
          int& x2,
          int& y2,
          const Header& head,
          bool flip_horizontal,
          bool flip_vertical);
          
///----------------------------------------
/// @brief Extract a single Bayer plane and write it to a sibling FITS file.
/// @details The suffix _<filtern>.fit is appended to the filename.
/// @param filename Input FITS file path.
/// @param filtern  Filter name: "TR", "TG", or "TB".
/// @param xp       X start position in the 2x2 Bayer cell.
/// @param yp       Y start position in the 2x2 Bayer cell.
/// @param memo     Header memo lines (modified).
/// @return Output filename, or empty string on failure.
///----------------------------------------

[[nodiscard]] std::string extract_raw_colour_to_file(const std::filesystem::path& filename,
                                                     const std::string& filtern,
                                                     int xp,
                                                     int yp,
                                                     std::vector<std::string>& memo);
                                                     
///----------------------------------------
/// @brief Extract Bayer planes from multiple files.
/// @param filenames List of input FITS file paths.
/// @param xp        X start position in the 2x2 Bayer cell.
/// @param yp        Y start position in the 2x2 Bayer cell.
/// @param filtern   Filter name: "TR", "TG", or "TB".
/// @param memo      Header memo lines (modified).
/// @return Number of files successfully written.
///----------------------------------------

[[nodiscard]] std::size_t split_raw(const std::vector<std::filesystem::path>& filenames,
                                    int xp,
                                    int yp,
                                    const std::string& filtern,
                                    std::vector<std::string>& memo);
                                    
///----------------------------------------
/// @brief Recompute the histogram for one channel (0=R, 1=G, 2=B).
/// @details Falls back to channel 0 when the requested channel does not exist.
/// @param colour Channel index (0, 1, or 2).
/// @param img    Source image.
///----------------------------------------

void get_hist(int colour, const ImageArray& img);

///----------------------------------------
/// @brief Recompute histograms and derive auto-stretch min/max and hist_range.
/// @param img          Source image.
/// @param update_hist  Whether to recompute histograms first.
/// @param range_index  Range selector (0..9; -1 same as 0).
/// @param head         FITS header.
/// @param[out] minm    Computed minimum value.
/// @param[out] maxm    Computed maximum value.
///----------------------------------------

void use_histogram(const ImageArray& img,
                   bool update_hist,
                   int range_index,
                   const Header& head,
                   int& minm,
                   int& maxm);
                   
///----------------------------------------
/// @brief Apply the current stretch to an image and return 16-bit-scaled result.
/// @param img Source image.
/// @param bck Background statistics.
/// @return Stretched image.
///----------------------------------------

[[nodiscard]] ImageArray stretch_image(const ImageArray& img, const Background& bck);
          
} // namespace
