///----------------------------------------
///      @file stack_driver.h
///   @ingroup ASTAP++
///     @brief Headless batch image-stacking driver (internal-astrometry path).
///   @details Reproduces the essential flow of the GUI stack button for the
///            internal-astrometry alignment mode: pre-solve every input frame,
///            drop any that will not solve, align + combine the survivors with
///            the chosen combiner, and write a stacked FITS carrying the ASTAP
///            stacked-header keyword block. Monochrome Average / Sigma-clip,
///            master-calibration creation and LRGB colour combine are provided;
///            OSC, mosaic and comet remain future work.
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

///----------------------------------------
///  @brief Per-channel inputs for an LRGB colour combine.
///  @details Each entry is a single, already-calibrated mono FITS for that
///           filter — typically a per-filter interim produced by @ref stack_files
///           (or a single exposure). @c rgb is an optional pre-combined
///           one-shot-colour frame added on top, and @c luminance an optional
///           lightness channel that modulates the chroma of the assembled RGB.
///           Supply at least two of @c red / @c green / @c blue (the Pascal
///           "minimum of two colours" rule); empty paths are skipped. The GUI's
///           filter-name classification that groups raw lights into these
///           channels is intentionally left to the caller (headless: hand in
///           the per-channel frames you want combined).
///----------------------------------------

struct LrgbInputs {
    std::filesystem::path red;
    std::filesystem::path green;
    std::filesystem::path blue;
    std::filesystem::path luminance;   // optional lightness channel
    std::filesystem::path rgb;         // optional pre-combined colour frame
};

///----------------------------------------
///  @brief Outcome of an LRGB colour combine.
///  @details @c channels_input counts the non-empty colour paths handed in,
///           @c channels_combined how many the combiner actually merged, and
///           @c reference is the frame used to fix the alignment grid (the
///           luminance if present, else the first colour available).
///----------------------------------------

struct LrgbResult {
    bool                  ok = false;
    int                   channels_input = 0;
    int                   channels_combined = 0;
    bool                  has_luminance = false;
    std::filesystem::path output;
    std::filesystem::path reference;
    std::string           message;
};

///----------------------------------------
///  @brief Assemble a colour (LRGB) FITS from per-channel mono frames.
///  @details Builds the reference/R/G/B/RGB/L slot array, aligns every channel
///           to the reference with the internal astrometric solver (each channel
///           is plate-solved and its solution persisted if it lacks a WCS),
///           combines them with an identity colour matrix, optionally modulates
///           the chroma by the luminance, and writes a 3-plane 32-bit-float FITS
///           carrying the LRGB stacked-header keyword block. Ports the
///           internal-astrometry path of @c stack_LRGB + the colour branch of
///           Tstackmenu1.stack_button1Click (`unit_stack.pas`).
///  @param inputs Per-channel frame paths (at least two colours required).
///  @param output Destination path for the colour FITS.
///  @return Result with per-channel counts and an @c ok flag.
///----------------------------------------

[[nodiscard]] LrgbResult combine_lrgb(const LrgbInputs& inputs,
                                      const std::filesystem::path& output);

} // namespace
