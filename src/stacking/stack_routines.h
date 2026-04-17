///----------------------------------------
///      @file stack_routines.h
///   @ingroup ASTAP++
///     @brief Core stacking algorithms: LRGB, simple average, mosaic, sigma-clip,
///            comet (ephemeris-aligned) and calibration-only passes, plus the
///            support helpers used by the stacking pipeline.
///   @details Provides the public stacking entry points along with the
///            coordinate-transform helper @ref calc_newx_newy that maps a pixel
///            of the current image into the reference-image grid, either via a
///            first-order vector solution or the full astrometric solution
///            (with optional SIP distortion correction).
///    @author Ported from Han Kleijn's unit_stack_routines.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <cmath>
#include <functional>
#include <span>
#include <string>
#include <vector>

#include "../types.h"
#include "../core/globals.h"
#include "../solving/astrometric_solving.h"
#include "../solving/star_align.h"

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

///----------------------------------------
/// MARK: Shared type aliases
///----------------------------------------

using astap::FileToDo;
using astap::ImageArray;
using astap::Header;
using astap::SolutionVector;

// Pull the star-align solution vectors into this namespace so the inline
// @ref calc_newx_newy below can reach them without qualification.
using astap::solving::solution_vectorX;
using astap::solving::solution_vectorY;

///----------------------------------------
/// MARK: UI callback sinks
///----------------------------------------

/// Message line emitted to the solver/stacker log.
using MemoSink = std::function<void(const std::string&)>;

/// Progress update: value in [0,100], label is a short status tag.
using ProgressSink = std::function<void(double, const std::string&)>;

///----------------------------------------
/// @brief Install a sink for solver/stacker log messages.
/// @details Pass a null function to uninstall. Thread-safe only if the sink
///          itself is — sinks are called from whatever thread the stacker
///          runs on.
///----------------------------------------

void set_memo2_sink(MemoSink sink);

///----------------------------------------
/// @brief Install a sink for progress updates.
///----------------------------------------

void set_progress_sink(ProgressSink sink);

///----------------------------------------
/// MARK: Public API
///----------------------------------------

///----------------------------------------
///   @brief Stack an LRGB set.
/// @details @p files_to_process is expected to contain the reference, R, G, B,
///          RGB and L slots in that order.
///   @param files_to_process List of files to stack; slot 0 is the reference.
/// @param[out] counter Number of frames successfully added.
///----------------------------------------

void stack_LRGB(std::span<FileToDo> files_to_process, int& counter);

///----------------------------------------
///   @brief Simple weighted-average stacker.
///   @param process_as_osc Non-zero to demosaic one-shot-colour inputs.
///   @param files_to_process Input frames (first entry is the reference).
/// @param[out] counter Number of frames successfully added.
///----------------------------------------

void stack_average(int process_as_osc,
                   std::span<FileToDo> files_to_process,
                   int& counter);
                   
///----------------------------------------
///   @brief Mosaic / tile-mode stacker.
/// @details Combines overlapping tiles into a single mosaic image.
///   @param process_as_osc Non-zero to demosaic one-shot-colour inputs.
///   @param files_to_process Input tiles.
///   @param max_dev_backgr Bound on the per-channel background correction
///                         between overlapping tiles.
/// @param[out] counter Number of tiles successfully added.
///----------------------------------------

void stack_mosaic(int process_as_osc,
                  std::span<FileToDo> files_to_process,
                  double max_dev_backgr,
                  int& counter);
                  
///----------------------------------------
///   @brief Sigma-clip average stacker.
/// @details Three-pass algorithm: average, variance, outlier rejection.
///   @param process_as_osc Non-zero to demosaic one-shot-colour inputs.
///   @param files_to_process Input frames.
/// @param[out] counter Number of frames successfully added.
///----------------------------------------

void stack_sigmaclip(int process_as_osc,
                     std::span<FileToDo> files_to_process,
                     int& counter);
                     
///----------------------------------------
///   @brief Calibration and alignment only.
/// @details Writes @c *_aligned.fit for every input; does not combine.
///   @param process_as_osc Non-zero to demosaic one-shot-colour inputs.
///   @param files_to_process Input frames.
/// @param[out] counter Number of frames successfully aligned.
///----------------------------------------

void calibration_and_alignment(int process_as_osc,
                               std::span<FileToDo> files_to_process,
                               int& counter);
                               
///----------------------------------------
///   @brief Ephemeris-based comet stacker.
/// @details Aligns the comet across frames while rejecting drifting-star
///          contributions.
///   @param process_as_osc Non-zero to demosaic one-shot-colour inputs.
///   @param files_to_process Input frames.
/// @param[out] counter Number of frames successfully combined.
///----------------------------------------

void stack_comet(int process_as_osc,
                 std::span<FileToDo> files_to_process,
                 int& counter);
                 
///----------------------------------------
/// @brief Convert the astrometric solution of the current image into a
///        first-order vector solution stored in @ref solution_vectorX and
///        @ref solution_vectorY.
///----------------------------------------

void astrometric_to_vector();

///----------------------------------------
/// @brief Pre-compute @c sin / @c cos of @c head.dec0 and @c head_ref.dec0.
///----------------------------------------

void initialise_calc_sincos_dec0();

///----------------------------------------
///  @brief Statistical test for the presence of a Bayer matrix.
/// @details Approximate; executes in roughly 1 ms for a 3040x2016 image.
///  @param img Image buffer to test.
/// @return @c true if a Bayer pattern is likely present.
///----------------------------------------

[[nodiscard]] bool test_bayer_matrix(const ImageArray& img);

///----------------------------------------
/// MARK: Inline helpers
///----------------------------------------

// External dependency forward-declared for use by the inline helper below.
// The full prototype and real implementation live in the astrometric-solving
// module currently being ported.
// TODO: remove once those are available via a real header.
void sincos(double angle, double& s, double& c);

///----------------------------------------
///   @brief Map a pixel of the current image onto the reference image.
/// @details Applies either the vector solution or the full astrometric
///          solution (with optional SIP distortion) to convert an input pixel
///          @c (fitsXfloat, fitsYfloat) to output @ref x_new_float and
///          @ref y_new_float. Inputs are 1-based (FITS convention); outputs
///          are 0-based (array indices, range @c 0..width-1 / @c 0..height-1).
///   @param vector_based When @c true, apply the affine vector solution;
///                       otherwise apply the full astrometric solution.
///   @param fitsXfloat FITS X coordinate (1-based).
///   @param fitsYfloat FITS Y coordinate (1-based).
///----------------------------------------

inline void calc_newx_newy(bool vector_based,
                           double fitsXfloat,
                           double fitsYfloat) {
    if (vector_based) {
        // x_new := a*X + b*Y + c (results in 0..width-1 array index range).
        x_new_float = solution_vectorX[0] * (fitsXfloat - 1.0)
                    + solution_vectorX[1] * (fitsYfloat - 1.0)
                    + solution_vectorX[2];
        y_new_float = solution_vectorY[0] * (fitsXfloat - 1.0)
                    + solution_vectorY[1] * (fitsYfloat - 1.0)
                    + solution_vectorY[2];
        return;
    }
    
    // Astrometric-based correction.
    // Step 6: (x,y) -> (RA,DEC) for the image to be added.
    auto u0 = fitsXfloat - head.crpix1;
    auto v0 = fitsYfloat - head.crpix2;
    
    auto u = 0.0;
    auto v = 0.0;
    if (a_order >= 2) {
        // SIP forward correction up to third order.
        u = u0 + a_0_0 + a_0_1 * v0 + a_0_2 * v0 * v0 + a_0_3 * v0 * v0 * v0
            + a_1_0 * u0 + a_1_1 * u0 * v0 + a_1_2 * u0 * v0 * v0
            + a_2_0 * u0 * u0 + a_2_1 * u0 * u0 * v0
            + a_3_0 * u0 * u0 * u0;
        v = v0 + b_0_0 + b_0_1 * v0 + b_0_2 * v0 * v0 + b_0_3 * v0 * v0 * v0
            + b_1_0 * u0 + b_1_1 * u0 * v0 + b_1_2 * u0 * v0 * v0
            + b_2_0 * u0 * u0 + b_2_1 * u0 * u0 * v0
            + b_3_0 * u0 * u0 * u0;
    } else {
        u = u0;
        v = v0;
    }
    
    constexpr auto kPi = 3.14159265358979323846;
    
    auto dRa  = (head.cd1_1 * u + head.cd1_2 * v) * kPi / 180.0;
    auto dDec = (head.cd2_1 * u + head.cd2_2 * v) * kPi / 180.0;
    const auto delta   = COS_dec0 - dDec * SIN_dec0;
    const auto gamma   = std::sqrt(dRa * dRa + delta * delta);
    const auto ra_new  = head.ra0 + std::atan(dRa / delta);
    const auto dec_new = std::atan((SIN_dec0 + dDec * COS_dec0) / gamma);
    
    // Step 5: (RA,DEC) -> (x,y) of reference image.
    auto sin_dec_new = 0.0;
    auto cos_dec_new = 0.0;
    sincos(dec_new, sin_dec_new, cos_dec_new);
    
    const auto delta_ra = ra_new - head_ref.ra0;
    auto sin_delta_ra = 0.0;
    auto cos_delta_ra = 0.0;
    sincos(delta_ra, sin_delta_ra, cos_delta_ra);
    
    const auto H = sin_dec_new * SIN_dec_ref
                 + cos_dec_new * COS_dec_ref * cos_delta_ra;
    dRa  = (cos_dec_new * sin_delta_ra / H) * 180.0 / kPi;
    dDec = ((sin_dec_new * COS_dec_ref
             - cos_dec_new * SIN_dec_ref * cos_delta_ra) / H) * 180.0 / kPi;
             
    const auto det = head_ref.cd2_2 * head_ref.cd1_1
                   - head_ref.cd1_2 * head_ref.cd2_1;
                   
    u0 = -(head_ref.cd1_2 * dDec - head_ref.cd2_2 * dRa) / det;
    v0 = +(head_ref.cd1_1 * dDec - head_ref.cd2_1 * dRa) / det;
    
    if (ap_order >= 2) {
        // SIP inverse correction up to third order.
        x_new_float = (head_ref.crpix1 + u0
                       + ap_0_1 * v0 + ap_0_2 * v0 * v0 + ap_0_3 * v0 * v0 * v0
                       + ap_1_0 * u0 + ap_1_1 * u0 * v0 + ap_1_2 * u0 * v0 * v0
                       + ap_2_0 * u0 * u0 + ap_2_1 * u0 * u0 * v0
                       + ap_3_0 * u0 * u0 * u0) - 1.0;
        y_new_float = (head_ref.crpix2 + v0
                       + bp_0_1 * v0 + bp_0_2 * v0 * v0 + bp_0_3 * v0 * v0 * v0
                       + bp_1_0 * u0 + bp_1_1 * u0 * v0 + bp_1_2 * u0 * v0 * v0
                       + bp_2_0 * u0 * u0 + bp_2_1 * u0 * u0 * v0
                       + bp_3_0 * u0 * u0 * u0) - 1.0;
    } else {
        x_new_float = (head_ref.crpix1 + u0) - 1.0;  // 0..width-1
        y_new_float = (head_ref.crpix2 + v0) - 1.0;
    }
}
    
} // namespace
