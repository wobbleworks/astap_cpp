///----------------------------------------
///      @file calc_trans_cubic.h
///   @ingroup ASTAP++
///     @brief Third-order transfer function between two sets of matching star positions.
///   @details Produces ten transfer coefficients for the X axis and ten for the
///            Y axis via a least-squares fit solved with Gauss-Jordan
///            elimination (partial pivoting).
///    @author Ported from Han Kleijn's unit_calc_trans_cubic.pas (ASTAP),
///            itself derived from Michael Richmond's Match package; the
///            ten-coefficient cubic extension is due to Cecile Melis (Siril).
///            MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / Michael Richmond / Cecile Melis / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <expected>
#include <string>

#include "../types.h"

///----------------------------------------
namespace astap::solving {
///----------------------------------------

using astap::SStar;
using astap::StarArray;
using astap::Matrix;
using astap::Vector;
using astap::TransCoeffs;

///----------------------------------------
///   @brief Solve for the 20 cubic transfer coefficients that map star
///          positions in @p stars_reference onto @p stars_distorted.
/// @details Uses a least-squares fit: builds two 10x10 normal-equation systems
///          (one per output axis) and solves each with Gauss-Jordan elimination.
///          Requires at least 10 matched pairs.
///   @param stars_reference Matched reference positions.
///   @param stars_distorted Matched distorted/target positions.
///  @return Twenty coefficients on success, or a diagnostic error string if
///          too few points were supplied or the normal-equation matrix is
///          singular / ill-conditioned.
///----------------------------------------

[[nodiscard]] std::expected<TransCoeffs, std::string>
calc_trans_cubic(const StarArray& stars_reference, const StarArray& stars_distorted);
	
} // namespace
