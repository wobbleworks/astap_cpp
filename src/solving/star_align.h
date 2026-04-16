///----------------------------------------
///      @file star_align.h
///   @ingroup ASTAP++
///     @brief Star-alignment routines: detect stars, build quad/triple patterns,
///            match them across images, and compute an affine solution vector.
///   @details Detect stars in an image, build quad/triple star-pattern
///            descriptors, match them between a reference and a target image,
///            and compute an affine solution vector (scale/rotation/offset) via
///            a least-squares fit using GIVENS rotations.
///    @author Ported from Han Kleijn's unit_star_align.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <string>
#include <vector>

#include "../types.h"
#include "calc_trans_cubic.h"  // SStar, StarArray (shared solver types)

///----------------------------------------
namespace astap::solving {
///----------------------------------------

/// 2xN / 8xN / 10xN double table used for star lists and quad descriptors.
using astap::StarList;

/// 3-element affine-solution vector.
using astap::SolutionVector;

/// 3-D float image buffer.
using astap::ImageArray;

///----------------------------------------
/// MARK: Module-level globals
///----------------------------------------

///----------------------------------------
/// @brief Quad descriptors for the reference (database) image.
/// @details 8xN: rows 0..5 hold sorted distances (row 0 = largest absolute
///          length, rows 1..5 are ratios relative to row 0), and rows 6..7
///          hold the quad centroid.
///----------------------------------------

extern int diag_nrquads1, diag_nrquads2, diag_pass1_matches, diag_pass2_matches;
extern StarList quad_star_distances1;

///----------------------------------------
/// @brief Quad descriptors for the target image. Same layout as @ref quad_star_distances1.
///----------------------------------------

extern StarList quad_star_distances2;

///----------------------------------------
/// @brief Least-squares A matrix, 3xN with rows (x, y, 1).
///----------------------------------------

extern StarList A_XYpositions;

///----------------------------------------
/// @brief Least-squares b vector: matching reference x positions.
///----------------------------------------

extern std::vector<double> b_Xrefpositions;

///----------------------------------------
/// @brief Least-squares b vector: matching reference y positions.
///----------------------------------------

extern std::vector<double> b_Yrefpositions;

///----------------------------------------
/// @brief Number of accepted quad correspondences after outlier rejection.
///----------------------------------------

extern int nr_references;

///----------------------------------------
/// @brief Number of quad correspondences before outlier rejection.
///----------------------------------------

extern int nr_references2;

///----------------------------------------
/// @brief Affine solution vector for the X axis: x' = [0]*x + [1]*y + [2].
///----------------------------------------

extern SolutionVector solution_vectorX;

///----------------------------------------
/// @brief Affine solution vector for the Y axis: y' = [0]*x + [1]*y + [2].
///----------------------------------------

extern SolutionVector solution_vectorY;

///----------------------------------------
/// @brief Per-channel black-level offsets used during stacking.
///----------------------------------------

extern SolutionVector solution_cblack;

///----------------------------------------
/// MARK: Canvas abstraction
///----------------------------------------

///----------------------------------------
///   @class IQuadCanvas
///   @brief Abstract drawing backend for @ref display_quads.
/// @details Decouples the module from any particular GUI toolkit.
///----------------------------------------

struct IQuadCanvas {
    virtual void draw_line(int x1, int y1, int x2, int y2) = 0;
    virtual ~IQuadCanvas() = default;
};

///----------------------------------------
/// MARK: Public API
///----------------------------------------

///----------------------------------------
///      @brief Find stars in @p img (hfd > hfd_min, snr > 10) and append their
///             (x, y) positions to @p starlist1 as a 2xN StarList.
///    @details If more than @p max_stars are found, the brightest @p max_stars
///             are kept.
///      @param img Source image buffer.
///      @param hfd_min Minimum half-flux diameter for a detection to be kept.
///      @param max_stars Maximum number of stars retained.
/// @param[out] starlist1 Output star list (2xN).
///----------------------------------------

void find_stars(const ImageArray& img, double hfd_min, int max_stars,
                const Background& bck, StarList& starlist1);

///----------------------------------------
///      @brief Build quad descriptors (largest length + 5 ratios + centroid)
///             from the supplied 2xN star list.
///      @param starlist Input star list; may be sorted in place.
/// @param[out] quad_star_distances Output 8xK quad descriptor table.
///----------------------------------------

void find_quads(StarList& starlist, StarList& quad_star_distances);

///----------------------------------------
///      @brief Same as @ref find_quads, but each detected quad is expanded into
///             four triples (123, 124, 134, 234).
///    @details Helps for low star counts (<30) where the brightest four stars
///             may differ between images.
///      @param starlist Input star list; may be sorted in place.
/// @param[out] quad_star_distances Output 8xK triple descriptor table.
///----------------------------------------

void find_triples_using_quads(StarList& starlist, StarList& quad_star_distances);

///----------------------------------------
///      @brief Build quads and return them as a 10xK list (8 star xy coords + xy centroid).
///    @details FOR DISPLAY ONLY. Used by @ref display_quads.
///      @param starlist Input star list.
/// @param[out] starlistquads Output 10xK table.
///----------------------------------------

void find_quads_xy(const StarList& starlist, StarList& starlistquads);

///----------------------------------------
/// @brief Match reference quads against image quads and solve the affine fit.
/// @details Populates the global least-squares inputs
///          (@ref A_XYpositions, @ref b_Xrefpositions, @ref b_Yrefpositions)
///          and solves for @ref solution_vectorX / @ref solution_vectorY.
///          On failure the solution vectors are reset to near-zero.
///   @param minimum_quads Required quad count (6 for stacking, 3 for plate-solving).
///   @param tolerance Maximum allowed deviation for ratio matching.
///  @return @c true on success, @c false otherwise.
///----------------------------------------

[[nodiscard]] bool find_offset_and_rotation(int minimum_quads, double tolerance);

///----------------------------------------
/// @brief Reset the affine solution vectors to a scaled identity.
/// @param factor Placed on the diagonal (normally 1.0; 0.001 is used as a "nullify" sentinel).
///----------------------------------------

void reset_solution_vectors(double factor) noexcept;

///----------------------------------------
/// @brief Draw the quads in @p starlistquads (10xK, produced by @ref find_quads_xy).
/// @param starlistquads Quad data (10xK).
/// @param canvas Drawing backend.
///----------------------------------------

void display_quads(const StarList& starlistquads, IQuadCanvas& canvas);

///----------------------------------------
///  @brief Human-readable rendering of @ref solution_vectorX / @ref solution_vectorY.
/// @return Formatted single-line summary of the current affine solution.
///----------------------------------------

[[nodiscard]] std::string solution_str();
    
} // namespace
