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
    std::filesystem::path reference;   // frame chosen as the alignment reference
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
///  @param master_dark Optional master dark subtracted from every frame (empty
///                     to skip). @param master_flat Optional master flat divided
///                     out of every frame (empty to skip). Both are applied
///                     per-frame by the combiner via @c apply_dark_and_flat.
///  @return Result with per-stage counts and an @c ok flag.
///----------------------------------------

[[nodiscard]] StackResult stack_files(std::span<const std::filesystem::path> frames,
                                      StackMethod method,
                                      const std::filesystem::path& output,
                                      const std::filesystem::path& master_dark = {},
                                      const std::filesystem::path& master_flat = {});

///----------------------------------------
///  @brief Kind of master calibration frame to create.
///  @details Selects only the count keyword written (DARK_CNT / FLAT_CNT /
///           BIAS_CNT); the combine itself is identical (a plain mean). A master
///           flat needs no pre-normalisation — @c apply_dark_and_flat normalises
///           it at application time from its centre.
///----------------------------------------

enum class MasterKind { Dark, Flat, Bias };

///----------------------------------------
///  @brief Outcome of a master-frame creation.
///----------------------------------------

struct MasterResult {
    bool                  ok = false;
    int                   frames_input = 0;
    int                   frames_combined = 0;
    std::filesystem::path output;
    std::string           message;
};

///----------------------------------------
///  @brief Mean-combine raw calibration frames into a mono master FITS.
///  @details Loads each frame, converts any colour input to mono, averages them,
///           and writes a 32-bit-float mono FITS carrying the @p kind count
///           keyword plus the average CCD temperature. Frames whose dimensions
///           differ from the first are skipped. Matches the core of the original
///           `replace_by_master_dark`/`replace_by_master_flat` (`unit_stack.pas`);
///           the GUI's exposure/temperature/gain classification is intentionally
///           left to the caller (headless: supply the frames you want combined).
///  @param frames Raw calibration-frame paths (must be non-empty).
///  @param kind Dark, Flat or Bias — selects the count keyword written.
///  @param output Destination path for the master FITS.
///  @return Result with the combined-frame count and an @c ok flag.
///----------------------------------------

[[nodiscard]] MasterResult create_master(std::span<const std::filesystem::path> frames,
                                         MasterKind kind,
                                         const std::filesystem::path& output);

} // namespace
