///----------------------------------------
///      @file stack_driver.h
///   @ingroup ASTAP++
///     @brief Headless batch image-stacking driver (internal-astrometry path).
///   @details Reproduces the essential flow of the GUI stack button for the
///            internal-astrometry alignment mode: pre-solve every input frame,
///            drop any that will not solve, align + combine the survivors with
///            the chosen combiner, and write a stacked FITS carrying the ASTAP
///            stacked-header keyword block. Monochrome Average / Sigma-clip
///            only for now; LRGB, OSC, master-calibration creation, mosaic and
///            comet remain future work.
///    @author Ported from Han Kleijn's unit_stack.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <filesystem>
#include <span>
#include <string>

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

/// @brief Frame-combination method.
enum class StackMethod { Average, SigmaClip };

///----------------------------------------
///  @brief Outcome of a batch stack.
///  @details @c frames_input is what was handed in, @c frames_solved is how many
///           survived the pre-solve pass (carry a usable WCS), and
///           @c frames_combined is how many the combiner actually merged.
///----------------------------------------

struct StackResult {
    bool                  ok = false;
    int                   frames_input = 0;
    int                   frames_solved = 0;
    int                   frames_combined = 0;
    std::filesystem::path output;
    std::string           message;
};

///----------------------------------------
///  @brief Pre-solve, align (internal astrometry), combine and save a batch.
///  @details Each frame is loaded and, if it has no WCS yet, plate-solved with
///           the internal solver and its solution persisted back into the file;
///           frames that will not solve are dropped (the reference frame is the
///           first survivor and must carry WCS). The survivors are then aligned
///           and combined with @p method, and the result written to @p output
///           as a 32-bit float FITS with the stacked-header keyword block.
///  @param frames Input light-frame paths (must be non-empty).
///  @param method Average or SigmaClip.
///  @param output Destination path for the stacked FITS.
///  @return Result with per-stage counts and an @c ok flag.
///----------------------------------------

[[nodiscard]] StackResult stack_files(std::span<const std::filesystem::path> frames,
                                      StackMethod method,
                                      const std::filesystem::path& output);

} // namespace
