///----------------------------------------
///      @file calc_trans_cubic.cpp
///   @ingroup ASTAP++
///     @brief Cubic transfer-function solver.
///    @author Ported from Han Kleijn's unit_calc_trans_cubic.pas (ASTAP),
///            itself derived from Michael Richmond's Match package; the
///            ten-coefficient cubic extension is due to Cecile Melis (Siril).
///            MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / Michael Richmond / Cecile Melis / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "calc_trans_cubic.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <string>
#include <utility>

///----------------------------------------
namespace astap::solving {
///----------------------------------------

///----------------------------------------
/// MARK: File-local helpers
///----------------------------------------

namespace {

using SolutionVector10 = std::array<double, 10>;

constexpr auto matrixTol = 1e-100;

///----------------------------------------
///  @brief Partial-pivot step of Gauss-Jordan elimination.
/// @details Given a square @p matrix (num x num) and a @p vector of @p num
///          elements, find the largest value in the matrix at or below the
///          given @p row (scaled by each row's largest absolute element, in
///          @p biggest_val). If that largest value is not already in @p row,
///          swap rows so it is - and swap the corresponding elements of
///          @p vector and @p biggest_val to match.
///  @param matrix Coefficient matrix (mutated).
///  @param num Matrix dimension.
///  @param vector Right-hand side (mutated).
///  @param biggest_val Per-row scale factors (mutated).
///  @param row Current pivot row.
///----------------------------------------

void gauss_pivot(Matrix& matrix,
                 int num,
                 SolutionVector10& vector,
                 Vector& biggest_val,
                 int row) {
    auto pivot_row = row;
    auto big = std::fabs(matrix[row][row] / biggest_val[row]);
    
    for (int i = row + 1; i < num; ++i) {
        const auto other_big = std::fabs(matrix[i][row] / biggest_val[i]);
        if (other_big > big) {
            big = other_big;
            pivot_row = i;
        }
    }
    
    if (pivot_row != row) {
        for (int col = row; col < num; ++col) {
            std::swap(matrix[pivot_row][col], matrix[row][col]);
        }
        std::swap(vector[pivot_row], vector[row]);
        std::swap(biggest_val[pivot_row], biggest_val[row]);
    }
}

///----------------------------------------
///   @brief Solve @c matrix*x=vector by Gaussian elimination with partial
///          pivoting and back-substitution.
/// @details See Chapra & Canale, @e Numerical @e Methods @e for @e Engineers
///          (McGraw-Hill 1998), chapter 9, for the algorithm. On success
///          overwrites @p vector with the solution and returns @c true;
///          otherwise fills @p err_mess and returns @c false.
///   @param matrix Coefficient matrix (mutated during elimination).
///   @param num Matrix dimension.
///   @param vector Right-hand side; overwritten with the solution.
///   @param err_mess Receives a diagnostic string on failure.
///  @return @c true on success, @c false if the matrix is singular or ill-conditioned.
///----------------------------------------

[[nodiscard]] bool gauss_matrix(Matrix& matrix,
                                int num,
                                SolutionVector10& vector,
                                std::string& err_mess) {
    err_mess.clear();
    
    auto biggest_val    = Vector(static_cast<std::size_t>(num), 0.0);
    auto solution_vector = Vector(static_cast<std::size_t>(num), 0.0);
    
    // Step 1: record the largest absolute value in each row; we use these to scale pivots.
    for (int i = 0; i < num; ++i) {
        biggest_val[i] = std::fabs(matrix[i][0]);
        for (int j = 1; j < num; ++j) {
            if (std::fabs(matrix[i][j]) > biggest_val[i]) {
                biggest_val[i] = std::fabs(matrix[i][j]);
            }
        }
        if (biggest_val[i] == 0.0) {
            err_mess = "Gauss_matrix: biggest val in row is zero";
            return false;
        }
    }
    
    // Step 2: forward elimination to upper-triangular form.
    for (int i = 0; i < num - 1; ++i) {
        gauss_pivot(matrix, num, vector, biggest_val, i);
        
        if (std::fabs(matrix[i][i] / biggest_val[i]) < matrixTol) {
            err_mess = "Gauss_matrix error: row has a too tiny value";
            return false;
        }
        
        for (int j = i + 1; j < num; ++j) {
            const auto factor = matrix[j][i] / matrix[i][i];
            for (int k = i + 1; k < num; ++k) {
                matrix[j][k] -= factor * matrix[i][k];
            }
            vector[j] -= factor * vector[i];
        }
    }
    
    if (std::fabs(matrix[num - 1][num - 1] / biggest_val[num - 1]) < matrixTol) {
        err_mess = "Gauss_matrix error: last row has a too tiny value";
        return false;
    }
    
    // Step 3: back-substitution.
    solution_vector[num - 1] = vector[num - 1] / matrix[num - 1][num - 1];
    for (int i = num - 2; i >= 0; --i) {
        auto sum = 0.0;
        for (int j = i + 1; j < num; ++j) {
            sum += matrix[i][j] * solution_vector[j];
        }
        solution_vector[i] = (vector[i] - sum) / matrix[i][i];
    }
    
    // Step 4: copy solution back into the caller's vector.
    for (int i = 0; i < num; ++i) {
        vector[i] = solution_vector[i];
    }
    return true;
}
    
} // namespace

///----------------------------------------
/// MARK: Public API
///----------------------------------------

std::expected<TransCoeffs, std::string>
calc_trans_cubic(const StarArray& stars_reference, const StarArray& stars_distorted) {
    // In the variable names below, a trailing '1' refers to a coordinate of
    // star s1 (the reference, appearing on both sides of the matrix equation)
    // and a '2' refers to star s2 (the distorted/target, appearing only on
    // the left-hand side).
    
    if (stars_reference.size() < 10) {
        return std::unexpected(std::string{"Calc_Trans_Cubic: Not enough equations."});
    }
    
    auto matrix = Matrix(10, Vector(10, 0.0));
    auto vector = SolutionVector10{};
    
    // Sums that make up the elements of the normal-equation matrix and its
    // two right-hand sides (one per output axis).
    auto sum = 0.0;
    auto sumx1 = 0.0, sumy1 = 0.0;
    auto sumx1sq = 0.0, sumx1y1 = 0.0, sumy1sq = 0.0;
    auto sumx1cu = 0.0, sumx1sqy1 = 0.0, sumx1y1sq = 0.0, sumy1cu = 0.0;
    auto sumx1qu = 0.0, sumx1cuy1 = 0.0, sumx1sqy1sq = 0.0, sumx1y1cu = 0.0, sumy1qu = 0.0;
    auto sumx1pe = 0.0, sumx1quy1 = 0.0, sumx1cuy1sq = 0.0, sumx1sqy1cu = 0.0, sumx1y1qu = 0.0, sumy1pe = 0.0;
    auto sumx1he = 0.0, sumx1pey1 = 0.0, sumx1quy1sq = 0.0, sumx1cuy1cu = 0.0,
         sumx1sqy1qu = 0.0, sumx1y1pe = 0.0, sumy1he = 0.0;
         
    auto sumx2 = 0.0, sumx2x1 = 0.0, sumx2y1 = 0.0;
    auto sumx2x1sq = 0.0, sumx2x1y1 = 0.0, sumx2y1sq = 0.0;
    auto sumx2x1cu = 0.0, sumx2x1sqy1 = 0.0, sumx2x1y1sq = 0.0, sumx2y1cu = 0.0;
    
    auto sumy2 = 0.0, sumy2x1 = 0.0, sumy2y1 = 0.0;
    auto sumy2x1sq = 0.0, sumy2x1y1 = 0.0, sumy2y1sq = 0.0;
    auto sumy2x1cu = 0.0, sumy2x1sqy1 = 0.0, sumy2x1y1sq = 0.0, sumy2y1cu = 0.0;
    
    // Take the minimum length of the two arrays (should be equal in normal
    // operation; defensive against mismatched input).
    const auto n = std::min(stars_reference.size(), stars_distorted.size());
    for (std::size_t i = 0; i < n; ++i) {
        const auto rx = stars_reference[i].x;
        const auto ry = stars_reference[i].y;
        const auto dx = stars_distorted[i].x;
        const auto dy = stars_distorted[i].y;
        
        sumx2       += dx;
        sumx2x1     += dx * rx;
        sumx2y1     += dx * ry;
        sumx2x1sq   += dx * rx * rx;
        sumx2x1y1   += dx * rx * ry;
        sumx2y1sq   += dx * ry * ry;
        sumx2x1cu   += dx * rx * rx * rx;
        sumx2x1sqy1 += dx * rx * rx * ry;
        sumx2x1y1sq += dx * rx * ry * ry;
        sumx2y1cu   += dx * ry * ry * ry;
        
        sumy2       += dy;
        sumy2x1     += dy * rx;
        sumy2y1     += dy * ry;
        sumy2x1sq   += dy * rx * rx;
        sumy2x1y1   += dy * rx * ry;
        sumy2y1sq   += dy * ry * ry;
        sumy2x1cu   += dy * rx * rx * rx;
        sumy2x1sqy1 += dy * rx * rx * ry;
        sumy2x1y1sq += dy * rx * ry * ry;
        sumy2y1cu   += dy * ry * ry * ry;
        
        // Matrix elements (moments of the reference coordinates).
        sum         += 1.0;
        sumx1       += rx;
        sumy1       += ry;
        sumx1sq     += rx * rx;
        sumx1y1     += rx * ry;
        sumy1sq     += ry * ry;
        sumx1cu     += rx * rx * rx;
        sumx1sqy1   += rx * rx * ry;
        sumx1y1sq   += rx * ry * ry;
        sumy1cu     += ry * ry * ry;
        sumx1qu     += rx * rx * rx * rx;
        sumx1cuy1   += rx * rx * rx * ry;
        sumx1sqy1sq += rx * rx * ry * ry;
        sumx1y1cu   += rx * ry * ry * ry;
        sumy1qu     += ry * ry * ry * ry;
        sumx1pe     += rx * rx * rx * rx * rx;
        sumx1quy1   += rx * rx * rx * rx * ry;
        sumx1cuy1sq += rx * rx * rx * ry * ry;
        sumx1sqy1cu += rx * rx * ry * ry * ry;
        sumx1y1qu   += rx * ry * ry * ry * ry;
        sumy1pe     += ry * ry * ry * ry * ry;
        sumx1he     += rx * rx * rx * rx * rx * rx;
        sumx1pey1   += rx * rx * rx * rx * rx * ry;
        sumx1quy1sq += rx * rx * rx * rx * ry * ry;
        sumx1cuy1cu += rx * rx * rx * ry * ry * ry;
        sumx1sqy1qu += rx * rx * ry * ry * ry * ry;
        sumx1y1pe   += rx * ry * ry * ry * ry * ry;
        sumy1he     += ry * ry * ry * ry * ry * ry;
    }
    
    // Fill the lower triangle of the symmetric normal-equation matrix, then
    // mirror it across the diagonal. Done twice (once per axis) because
    // gauss_matrix overwrites @c matrix in place during elimination.
    auto fill_lower_and_transpose = [&]() {
        // column 0
        matrix[0][0] = sum;
        matrix[1][0] = sumx1;
        matrix[2][0] = sumy1;
        matrix[3][0] = sumx1sq;
        matrix[4][0] = sumx1y1;
        matrix[5][0] = sumy1sq;
        matrix[6][0] = sumx1cu;
        matrix[7][0] = sumx1sqy1;
        matrix[8][0] = sumx1y1sq;
        matrix[9][0] = sumy1cu;
        // column 1
        matrix[1][1] = sumx1sq;
        matrix[2][1] = sumx1y1;
        matrix[3][1] = sumx1cu;
        matrix[4][1] = sumx1sqy1;
        matrix[5][1] = sumx1y1sq;
        matrix[6][1] = sumx1qu;
        matrix[7][1] = sumx1cuy1;
        matrix[8][1] = sumx1sqy1sq;
        matrix[9][1] = sumx1y1cu;
        // column 2
        matrix[2][2] = sumy1sq;
        matrix[3][2] = sumx1sqy1;
        matrix[4][2] = sumx1y1sq;
        matrix[5][2] = sumy1cu;
        matrix[6][2] = sumx1cuy1;
        matrix[7][2] = sumx1sqy1sq;
        matrix[8][2] = sumx1y1cu;
        matrix[9][2] = sumy1qu;
        // column 3
        matrix[3][3] = sumx1qu;
        matrix[4][3] = sumx1cuy1;
        matrix[5][3] = sumx1sqy1sq;
        matrix[6][3] = sumx1pe;
        matrix[7][3] = sumx1quy1;
        matrix[8][3] = sumx1cuy1sq;
        matrix[9][3] = sumx1sqy1cu;
        // column 4
        matrix[4][4] = sumx1sqy1sq;
        matrix[5][4] = sumx1y1cu;
        matrix[6][4] = sumx1quy1;
        matrix[7][4] = sumx1cuy1sq;
        matrix[8][4] = sumx1sqy1cu;
        matrix[9][4] = sumx1y1qu;
        // column 5
        matrix[5][5] = sumy1qu;
        matrix[6][5] = sumx1cuy1sq;
        matrix[7][5] = sumx1sqy1cu;
        matrix[8][5] = sumx1y1qu;
        matrix[9][5] = sumy1pe;
        // column 6
        matrix[6][6] = sumx1he;
        matrix[7][6] = sumx1pey1;
        matrix[8][6] = sumx1quy1sq;
        matrix[9][6] = sumx1cuy1cu;
        // column 7
        matrix[7][7] = sumx1quy1sq;
        matrix[8][7] = sumx1cuy1cu;
        matrix[9][7] = sumx1sqy1qu;
        // column 8
        matrix[8][8] = sumx1sqy1qu;
        matrix[9][8] = sumx1y1pe;
        // column 9
        matrix[9][9] = sumy1he;
        
        for (int r = 0; r <= 8; ++r) {
            for (int c = r + 1; c <= 9; ++c) {
                matrix[r][c] = matrix[c][r];
            }
        }
    };
    
    // Solve for x-axis coefficients (x00..x03).
    fill_lower_and_transpose();
    vector[0] = sumx2;
    vector[1] = sumx2x1;
    vector[2] = sumx2y1;
    vector[3] = sumx2x1sq;
    vector[4] = sumx2x1y1;
    vector[5] = sumx2y1sq;
    vector[6] = sumx2x1cu;
    vector[7] = sumx2x1sqy1;
    vector[8] = sumx2x1y1sq;
    vector[9] = sumx2y1cu;
    
    auto err_mess = std::string{};
    if (!gauss_matrix(matrix, 10, vector, err_mess)) {
        err_mess += ", Calc_trans_cubic: can not solve for coeffs A,B,C,D,E,F,G,H,I,J";
        return std::unexpected(std::move(err_mess));
    }
    
    auto trans = TransCoeffs{};
    trans.x00 = vector[0];
    trans.x10 = vector[1];
    trans.x01 = vector[2];
    trans.x20 = vector[3];
    trans.x11 = vector[4];
    trans.x02 = vector[5];
    trans.x30 = vector[6];
    trans.x21 = vector[7];
    trans.x12 = vector[8];
    trans.x03 = vector[9];
    
    // Solve for y-axis coefficients (y00..y03). The matrix was mutated by the
    // previous elimination, so rebuild it.
    fill_lower_and_transpose();
    vector[0] = sumy2;
    vector[1] = sumy2x1;
    vector[2] = sumy2y1;
    vector[3] = sumy2x1sq;
    vector[4] = sumy2x1y1;
    vector[5] = sumy2y1sq;
    vector[6] = sumy2x1cu;
    vector[7] = sumy2x1sqy1;
    vector[8] = sumy2x1y1sq;
    vector[9] = sumy2y1cu;
    
    if (!gauss_matrix(matrix, 10, vector, err_mess)) {
        err_mess += ", Calc_trans_cubic: Can not solve for coeffs y00,y10,y01,y20,y11,y02,y30,y21,y12,y03";
        return std::unexpected(std::move(err_mess));
    }
    
    trans.y00 = vector[0];
    trans.y10 = vector[1];
    trans.y01 = vector[2];
    trans.y20 = vector[3];
    trans.y11 = vector[4];
    trans.y02 = vector[5];
    trans.y30 = vector[6];
    trans.y21 = vector[7];
    trans.y12 = vector[8];
    trans.y03 = vector[9];
    
    return trans;
}
    
} // namespace
