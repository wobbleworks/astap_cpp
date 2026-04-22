///----------------------------------------
///      @file star_align.cpp
///   @ingroup ASTAP++
///     @brief Star-alignment implementation.
///    @author Ported from Han Kleijn's unit_star_align.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "star_align.h"

#include "../core/globals.h"
#include "../core/photometry.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <format>
#include <limits>
#include <random>
#include <ranges>
#include <string>
#include <vector>

///----------------------------------------
namespace astap::solving {
///----------------------------------------

///----------------------------------------
/// MARK: Global state
///----------------------------------------

StarList            quad_star_distances1;
StarList            quad_star_distances2;
int                 diag_nrquads1 = 0;
int                 diag_nrquads2 = 0;
int                 diag_pass1_matches = 0;
int                 diag_pass2_matches = 0;
StarList            A_XYpositions;
std::vector<double> b_Xrefpositions;
std::vector<double> b_Yrefpositions;
int                 nr_references  = 0;
int                 nr_references2 = 0;
SolutionVector      solution_vectorX{};
SolutionVector      solution_vectorY{};
SolutionVector      solution_cblack{};

///----------------------------------------
/// MARK: File-local helpers
///----------------------------------------

namespace {

/// @brief Fixed-point format with 5 fractional digits and dot separator.
[[nodiscard]] std::string float_to_str6(double x) {
    return std::format("{:.5f}", x);
}

/// @brief Allocate a row-major 2D star-list array of shape [rows][cols].
void set_length_2d(StarList& a, std::size_t rows, std::size_t cols) {
    a.assign(rows, std::vector<double>(cols, 0.0));
}

/// @brief Resize both dimensions, preserving as much existing content as fits.
void resize_2d_preserve(StarList& a, std::size_t rows, std::size_t cols) {
    a.resize(rows);
    for (auto& row : a) {
        row.resize(cols);
    }
}

/// @brief Quick sort by row 0 (ascending), swapping rows 0 and 1 together.
void quicksort_starlist(StarList& A, int iLo, int iHi) {
    auto Lo = iLo;
    auto Hi = iHi;
    const auto pivot = A[0][(Lo + Hi) / 2];
    
    do {
        while (A[0][Lo] < pivot) {
            ++Lo;
        }
        while (A[0][Hi] > pivot) {
            --Hi;
        }
        if (Lo <= Hi) {
            std::swap(A[0][Lo], A[0][Hi]);
            std::swap(A[1][Lo], A[1][Hi]);
            ++Lo;
            --Hi;
        }
    } while (Lo <= Hi);
    
    if (Hi > iLo) {
        quicksort_starlist(A, iLo, Hi);
    }
    if (Lo < iHi) {
        quicksort_starlist(A, Lo, iHi);
    }
}

///----------------------------------------
///   @brief Solve the overdetermined system @c A*x=b by least squares using
///          GIVENS rotations.
/// @details See Montenbruck & Pfleger, @e Astronomy @e on @e the @e personal
///          @e computer. @p A_matrix is [nr_columns][nr_equations]; @p x_matrix
///          has length nr_columns (3).
///
///          Zero the lower triangle column by column: pick p, q with p^2+q^2=1
///          so that rotating rows (j, i) annihilates A[j,i], keeping the
///          residual norm invariant. Then back-substitute.
///   @param A_matrix Coefficient matrix (column-major: [column][equation]).
///   @param b_matrix Right-hand side; taken by value so it can be mutated.
///   @param x_matrix Output solution vector.
///  @return @c true on success, @c false if the system is singular.
///----------------------------------------

[[nodiscard]] bool lsq_fit(const StarList& A_matrix,
                           std::vector<double> b_matrix,
                           SolutionVector& x_matrix) {
    constexpr auto tiny = 1e-10;
    
    const auto nr_columns   = static_cast<int>(A_matrix.size());       // 3
    const auto nr_equations = static_cast<int>(A_matrix[0].size());
    
    // Duplicate A so the caller's data is untouched.
    auto temp_matrix = A_matrix;
    
    for (int j = 0; j < nr_columns; ++j) {
        for (int i = j + 1; i < nr_equations; ++i) {
            if (temp_matrix[j][i] == 0.0) {
                continue;
            }
            
            // Determine the Givens coefficients p, q.
            double p = 0.0;
            double q = 0.0;
            if (std::fabs(temp_matrix[j][j]) < tiny * std::fabs(temp_matrix[j][i])) {
                p = 0.0;
                q = 1.0;
                temp_matrix[j][j] = -temp_matrix[j][i];
                temp_matrix[j][i] = 0.0;
            }
            else {
                // h = +-sqrt(A11^2 + A21^2); p = A11/h; q = -A21/h
                auto h = std::sqrt(temp_matrix[j][j] * temp_matrix[j][j] +
                                   temp_matrix[j][i] * temp_matrix[j][i]);
                if (temp_matrix[j][j] < 0.0) {
                    h = -h;
                }
                p = temp_matrix[j][j] / h;
                q = -temp_matrix[j][i] / h;
                temp_matrix[j][j] = h;
                temp_matrix[j][i] = 0.0;
            }
            
            // Apply the same rotation to the rest of the row.
            for (int k = j + 1; k < nr_columns; ++k) {
                const auto h = p * temp_matrix[k][j] - q * temp_matrix[k][i];
                temp_matrix[k][i] = q * temp_matrix[k][j] + p * temp_matrix[k][i];
                temp_matrix[k][j] = h;
            }
            const auto h = p * b_matrix[j] - q * b_matrix[i];
            b_matrix[i] = q * b_matrix[j] + p * b_matrix[i];
            b_matrix[j] = h;
        }
    }
    
    // Zero the output before back-substitution.
    x_matrix.fill(0.0);
    
    // Back-substitute upper-triangular system.
    for (int i = nr_columns - 1; i >= 0; --i) {
        auto h = b_matrix[i];
        for (int k = i + 1; k < nr_columns; ++k) {
            h -= temp_matrix[k][i] * x_matrix[k];
        }
        if (std::fabs(temp_matrix[i][i]) > 1e-30) {
            x_matrix[i] = h / temp_matrix[i][i];
        }
        else {
            return false;  // singular: fail rather than divide by zero
        }
    }
    return true;
}

///----------------------------------------
/// @brief Keep the brightest @p nr_stars_required stars from @p starlist1.
/// @details Uses an SNR histogram with 201 bins scaled against @p highest_snr.
///   @param nr_stars_required Target count.
///   @param highest_snr Maximum SNR observed (for histogram scaling).
///   @param snr_list SNR value per star in @p starlist1.
///   @param starlist1 Star list filtered in place.
///----------------------------------------

void get_brightest_stars(int nr_stars_required,
                         double highest_snr,
                         const std::vector<double>& snr_list,
                         StarList& starlist1) {
    constexpr auto range = 200;
    auto snr_histogram = std::array<int, range + 1>{};
    
    // Populate SNR histogram.
    for (const auto& snr : snr_list) {
        auto snr_scaled = static_cast<int>(snr * range / highest_snr);
        if (snr_scaled < 0) {
            snr_scaled = 0;
        }
        if (snr_scaled > range) {
            snr_scaled = range;
        }
        ++snr_histogram[snr_scaled];
    }
    
    // Walk bins from brightest until we have enough stars.
    auto count = 0;
    auto i = range + 1;
    do {
        --i;
        count += snr_histogram[i];
    } while (!(i <= 0 || count >= nr_stars_required));
    
    const auto snr_required = highest_snr * i / range;
    
    // Compact kept stars in place.
    count = 0;
    const auto nrstars = static_cast<int>(starlist1[0].size());
    for (int k = 0; k < nrstars; ++k) {
        if (snr_list[k] >= snr_required) {
            starlist1[0][count] = starlist1[0][k];
            starlist1[1][count] = starlist1[1][k];
            ++count;
        }
    }
    resize_2d_preserve(starlist1, 2, static_cast<std::size_t>(count));
}

///----------------------------------------
///  @brief Core of @ref find_offset_and_rotation.
/// @details Cross-matches quad descriptors, filters on the median of the
///          longest-length ratio, and populates the LSQ system.
///  @param minimum_count Required count of accepted pairs.
///  @param quad_tolerance Ratio-matching tolerance.
/// @return @c true on success, @c false if too few matches survive.
///----------------------------------------

[[nodiscard]] bool find_fit(int minimum_count, double quad_tolerance) {
    const auto nrquads1 = static_cast<int>(quad_star_distances1[0].size());
    const auto nrquads2 = static_cast<int>(quad_star_distances2[0].size());
    diag_nrquads1 = nrquads1;
    diag_nrquads2 = nrquads2;
    diag_pass1_matches = 0;
    diag_pass2_matches = 0;

    if (nrquads1 < minimum_count || nrquads2 < minimum_count) {
        nr_references = 0;
        return false;
    }
    
    // Pass 1: all quad pairs whose 5 length-ratios are within tolerance.
    auto matchlist2 = std::vector<std::array<int, 2>>{};
    matchlist2.reserve(1000);
    nr_references2 = 0;
    
    for (int i = 0; i < nrquads1; ++i) {
        for (int j = 0; j < nrquads2; ++j) {
            if (std::fabs(quad_star_distances1[1][i] - quad_star_distances2[1][j]) <= quad_tolerance &&
                std::fabs(quad_star_distances1[2][i] - quad_star_distances2[2][j]) <= quad_tolerance &&
                std::fabs(quad_star_distances1[3][i] - quad_star_distances2[3][j]) <= quad_tolerance &&
                std::fabs(quad_star_distances1[4][i] - quad_star_distances2[4][j]) <= quad_tolerance &&
                std::fabs(quad_star_distances1[5][i] - quad_star_distances2[5][j]) <= quad_tolerance) {
                matchlist2.push_back({i, j});
                ++nr_references2;
            }
        }
    }
    
    diag_pass1_matches = nr_references2;
    if (nr_references2 < minimum_count) {
        nr_references = 0;
        return false;
    }
    
    // Median of the longest-length ratio across candidate matches.
    auto ratios = std::vector<double>(static_cast<std::size_t>(nr_references2));
    for (int k = 0; k < nr_references2; ++k) {
        ratios[k] = quad_star_distances1[0][matchlist2[k][0]] /
                    quad_star_distances2[0][matchlist2[k][1]];
    }
    auto sorted_ratios = ratios;
    std::ranges::nth_element(sorted_ratios, sorted_ratios.begin() + sorted_ratios.size() / 2);
    const auto median_ratio = sorted_ratios[sorted_ratios.size() / 2];
    
    // Pass 2: drop outliers in longest-length ratio.
    nr_references = 0;
    auto matchlist1 = std::vector<std::array<int, 2>>{};
    matchlist1.reserve(1000);
    for (int k = 0; k < nr_references2; ++k) {
        if (std::fabs(median_ratio - ratios[k]) <= quad_tolerance * median_ratio) {
            matchlist1.push_back(matchlist2[k]);
            ++nr_references;
        }
        // TODO: external dep (logger) - log dropped outliers.
    }
    
    diag_pass2_matches = nr_references;
    if (nr_references < 3) {
        matchlist1.clear();
        matchlist2.clear();
        return false;
    }
    
    // Fill the 3xN LSQ system: rows (x, y, 1). b_* are the ref centroid coords.
    set_length_2d(A_XYpositions, 3, static_cast<std::size_t>(nr_references));
    b_Xrefpositions.assign(static_cast<std::size_t>(nr_references), 0.0);
    b_Yrefpositions.assign(static_cast<std::size_t>(nr_references), 0.0);
    
    for (int k = 0; k < nr_references; ++k) {
        A_XYpositions[0][k] = quad_star_distances2[6][matchlist1[k][1]];
        A_XYpositions[1][k] = quad_star_distances2[7][matchlist1[k][1]];
        A_XYpositions[2][k] = 1.0;
        b_Xrefpositions[k]  = quad_star_distances1[6][matchlist1[k][0]];
        b_Yrefpositions[k]  = quad_star_distances1[7][matchlist1[k][0]];
    }
    return true;
}
 
} // namespace

///----------------------------------------
/// MARK: Public API
///----------------------------------------

std::string solution_str() {
    return std::format(
        "Solution[px] x:={}x+ {}y+ {},  y:={}x+ {}y+ {}",
        float_to_str6(solution_vectorX[0]),
        float_to_str6(solution_vectorX[1]),
        float_to_str6(solution_vectorX[2]),
        float_to_str6(solution_vectorY[0]),
        float_to_str6(solution_vectorY[1]),
        float_to_str6(solution_vectorY[2]));
}

void reset_solution_vectors(double factor) noexcept {
    // Identity (scaled) on the diagonal.
    solution_vectorX[0] = factor;
    solution_vectorX[1] = 0.0;
    solution_vectorX[2] = 0.0;
    
    solution_vectorY[0] = 0.0;
    solution_vectorY[1] = factor;
    solution_vectorY[2] = 0.0;
}

void display_quads(const StarList& starlistquads, IQuadCanvas& canvas) {
    // TODO: external dep - original guards on head.naxis==0 (file loaded?),
    // sets pen thickness from head.height / mainwindow.image1.height, and uses
    // mainwindow.flip_horizontal1/flip_vertical1 for the axis origin. These
    // must be provided by the concrete IQuadCanvas implementation now.
    
    if (starlistquads.empty() || starlistquads[0].empty()) {
        return;
    }
    
    const auto nrquads = static_cast<int>(starlistquads[0].size()) - 1;
    
    // Matches the original "coin-flip" pen colour / solid-vs-dotted style.
    // A real canvas would apply these per-line; here we just iterate.
    [[maybe_unused]] auto rng = std::mt19937{12345};
    [[maybe_unused]] auto colour_dist = std::uniform_int_distribution<int>{0, 0x9F9F9F};
    
    // The canvas is expected to apply any flipping in its own coordinate system.
    constexpr auto x = 0;
    constexpr auto y = 0;
    constexpr auto flipx = 1;
    constexpr auto flipy = 1;
    
    for (int i = 0; i <= nrquads; ++i) {
        const auto px = [&](int row) { return x + flipx * static_cast<int>(std::lround(starlistquads[row][i])); };
        const auto py = [&](int row) { return y + flipy * static_cast<int>(std::lround(starlistquads[row][i])); };
        
        const auto x1 = px(0);
        const auto y1 = py(1);
        const auto x2 = px(2);
        const auto y2 = py(3);
        const auto x3 = px(4);
        const auto y3 = py(5);
        const auto x4 = px(6);
        const auto y4 = py(7);
        
        canvas.draw_line(x1, y1, x2, y2);  // star1 -> star2
        canvas.draw_line(x2, y2, x3, y3);  // star2 -> star3
        canvas.draw_line(x3, y3, x1, y1);  // star3 -> star1
        canvas.draw_line(x1, y1, x4, y4);  // star1 -> star4
        canvas.draw_line(x4, y4, x3, y3);  // star4 -> star3
        canvas.draw_line(x4, y4, x2, y2);  // star4 -> star2
    }
    // TODO: external dep - memo2_message('<n> quads found.')
}

void find_quads(StarList& starlist, StarList& quad_star_distances) {
    const auto nrstars = static_cast<int>(starlist[0].size());
    
    if (nrstars < 4) {
        set_length_2d(quad_star_distances, 8, 0);
        return;
    }
    
    // Tolerance band ~ 2x average stellar spacing if stars are uniform.
    auto tolerance = 1;  // pre-filter off by default
    if (nrstars >= 150) {
        quicksort_starlist(starlist, 0, nrstars - 1);  // sort in X
        tolerance = static_cast<int>(std::lround(0.5 * std::sqrt(static_cast<double>(nrstars))));
    }
    
    auto nrquads = 0;
    set_length_2d(quad_star_distances, 8, static_cast<std::size_t>(nrstars));
    
    auto j_distance1 = 0;
    auto j_distance2 = 0;
    auto j_distance3 = 0;
    
    for (int i = 0; i < nrstars; ++i) {
        auto distance1 = 1e99;
        auto distance2 = 1e99;
        auto distance3 = 1e99;
        
        const auto band   = nrstars / tolerance;
        const auto Sstart = std::max(0, i - band);
        const auto Send   = std::min(nrstars - 1, i + band);
        
        // Track the three nearest neighbours of star i.
        for (int j = Sstart; j <= Send; ++j) {
            if (j == i) {
                continue;
            }
            const auto dy    = starlist[1][j] - starlist[1][i];
            const auto disty = dy * dy;
            if (disty < distance3) {
                const auto dx = starlist[0][j] - starlist[0][i];
                const auto distance = dx * dx + disty;
                if (distance > 1.0) {  // not an identical star
                    if (distance < distance1) {
                        distance3 = distance2; j_distance3 = j_distance2;
                        distance2 = distance1; j_distance2 = j_distance1;
                        distance1 = distance;  j_distance1 = j;
                    }
                    else if (distance < distance2) {
                        distance3 = distance2; j_distance3 = j_distance2;
                        distance2 = distance;  j_distance2 = j;
                    }
                    else if (distance < distance3) {
                        distance3 = distance;  j_distance3 = j;
                    }
                }
            }
        }
        
        // Gather the four stars of the quad.
        const auto x1 = starlist[0][i];
        const auto y1 = starlist[1][i];
        const auto x2 = starlist[0][j_distance1];
        const auto y2 = starlist[1][j_distance1];
        const auto x3 = starlist[0][j_distance2];
        const auto y3 = starlist[1][j_distance2];
        const auto x4 = starlist[0][j_distance3];
        const auto y4 = starlist[1][j_distance3];
        
        const auto xt = (x1 + x2 + x3 + x4) / 4.0;
        const auto yt = (y1 + y2 + y3 + y4) / 4.0;
        
        // Skip quads whose centroid coincides with an earlier one.
        auto identical = false;
        for (int k = 0; k < nrquads; ++k) {
            if (std::fabs(xt - quad_star_distances[6][k]) < 1.0 &&
                std::fabs(yt - quad_star_distances[7][k]) < 1.0) {
                identical = true;
                break;
            }
        }
        if (identical) {
            continue;
        }
        
        // Compute the six inter-star distances.
        auto dist1 = std::sqrt(distance1);
        auto dist2 = std::sqrt(distance2);
        auto dist3 = std::sqrt(distance3);
        auto dist4 = std::sqrt((x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3));
        auto dist5 = std::sqrt((x2 - x4) * (x2 - x4) + (y2 - y4) * (y2 - y4));
        auto dist6 = std::sqrt((x3 - x4) * (x3 - x4) + (y3 - y4) * (y3 - y4));
        
        // 5-pass bubble sort (descending). Keeps the original sort order.
        for (int j = 1; j <= 5; ++j) {
            if (dist6 > dist5) {
                std::swap(dist5, dist6);
            }
            if (dist5 > dist4) {
                std::swap(dist4, dist5);
            }
            if (dist4 > dist3) {
                std::swap(dist3, dist4);
            }
            if (dist3 > dist2) {
                std::swap(dist2, dist3);
            }
            if (dist2 > dist1) {
                std::swap(dist1, dist2);
            }
        }
        
        // Store the descriptor row.
        quad_star_distances[0][nrquads] = dist1;         // largest
        quad_star_distances[1][nrquads] = dist2 / dist1; // scaled ratios
        quad_star_distances[2][nrquads] = dist3 / dist1;
        quad_star_distances[3][nrquads] = dist4 / dist1;
        quad_star_distances[4][nrquads] = dist5 / dist1;
        quad_star_distances[5][nrquads] = dist6 / dist1;
        quad_star_distances[6][nrquads] = xt;
        quad_star_distances[7][nrquads] = yt;
        ++nrquads;
    }
    resize_2d_preserve(quad_star_distances, 8, static_cast<std::size_t>(nrquads));
}

void find_triples_using_quads(StarList& starlist, StarList& quad_star_distances) {
    const auto nrstars = static_cast<int>(starlist[0].size());
    
    if (nrstars < 4) {
        set_length_2d(quad_star_distances, 8, 0);
        return;
    }
    
    auto tolerance = 1;
    if (nrstars >= 150) {
        quicksort_starlist(starlist, 0, nrstars - 1);
        tolerance = static_cast<int>(std::lround(0.5 * std::sqrt(static_cast<double>(nrstars))));
    }
    
    auto nrquads     = 0;  // triples-as-quads count
    auto nrrealquads = 0;  // distinct underlying quads
    set_length_2d(quad_star_distances, 8, static_cast<std::size_t>(nrstars * 4));
    StarList quad_centers;
    set_length_2d(quad_centers, 2, static_cast<std::size_t>(nrstars));
    
    auto j_distance1 = 0;
    auto j_distance2 = 0;
    auto j_distance3 = 0;
    
    // Closure-equivalent: add a new (xt, yt, sorted d1/d2/d3) triple record,
    // skipping duplicates by centroid.
    auto get_triple = [&](double x1, double y1, double x2, double y2,
                          double x3, double y3,
                          double d1, double d2, double d3) {
        const auto xt = (x1 + x2 + x3) / 3.0;
        const auto yt = (y1 + y2 + y3) / 3.0;
        
        for (int k = 0; k < nrquads; ++k) {
            if (std::fabs(xt - quad_star_distances[6][k]) < 1.0 &&
                std::fabs(yt - quad_star_distances[7][k]) < 1.0) {
                return;
            }
        }
        // Sort three distances descending (two passes).
        for (int j = 1; j <= 2; ++j) {
            if (d3 > d2) {
                std::swap(d2, d3);
            }
            if (d2 > d1) {
                std::swap(d1, d2);
            }
        }
        quad_star_distances[0][nrquads] = d1;
        quad_star_distances[1][nrquads] = d2 / d1;
        quad_star_distances[2][nrquads] = d3 / d1;
        quad_star_distances[3][nrquads] = 0.0;
        quad_star_distances[4][nrquads] = 0.0;
        quad_star_distances[5][nrquads] = 0.0;
        quad_star_distances[6][nrquads] = xt;
        quad_star_distances[7][nrquads] = yt;
        ++nrquads;
    };
    
    for (int i = 0; i < nrstars; ++i) {
        auto distance1 = 1e99;
        auto distance2 = 1e99;
        auto distance3 = 1e99;
        
        const auto band   = nrstars / tolerance;
        const auto Sstart = std::max(0, i - band);
        const auto Send   = std::min(nrstars - 1, i + band);
        
        for (int j = Sstart; j <= Send; ++j) {
            if (j == i) {
                continue;
            }
            const auto dy    = starlist[1][j] - starlist[1][i];
            const auto disty = dy * dy;
            if (disty < distance3) {
                const auto dx = starlist[0][j] - starlist[0][i];
                const auto distance = dx * dx + disty;
                if (distance > 1.0) {
                    if (distance < distance1) {
                        distance3 = distance2; j_distance3 = j_distance2;
                        distance2 = distance1; j_distance2 = j_distance1;
                        distance1 = distance;  j_distance1 = j;
                    }
                    else if (distance < distance2) {
                        distance3 = distance2; j_distance3 = j_distance2;
                        distance2 = distance;  j_distance2 = j;
                    }
                    else if (distance < distance3) {
                        distance3 = distance;  j_distance3 = j;
                    }
                }
            }
        }
        
        const auto x1a = starlist[0][i];
        const auto y1a = starlist[1][i];
        const auto x2a = starlist[0][j_distance1];
        const auto y2a = starlist[1][j_distance1];
        const auto x3a = starlist[0][j_distance2];
        const auto y3a = starlist[1][j_distance2];
        const auto x4a = starlist[0][j_distance3];
        const auto y4a = starlist[1][j_distance3];
        
        const auto xt = (x1a + x2a + x3a + x4a) / 4.0;
        const auto yt = (y1a + y2a + y3a + y4a) / 4.0;
        
        auto identical = false;
        for (int k = 0; k < nrrealquads; ++k) {
            if (std::fabs(xt - quad_centers[0][k]) < 1.0 &&
                std::fabs(yt - quad_centers[1][k]) < 1.0) {
                identical = true;
                break;
            }
        }
        if (identical) {
            continue;
        }
        
        quad_centers[0][nrrealquads] = xt;
        quad_centers[1][nrrealquads] = yt;
        ++nrrealquads;
        
        const auto dist12 = std::sqrt(distance1);
        const auto dist13 = std::sqrt(distance2);
        const auto dist14 = std::sqrt(distance3);
        const auto dist23 = std::sqrt((x2a - x3a) * (x2a - x3a) + (y2a - y3a) * (y2a - y3a));
        const auto dist24 = std::sqrt((x2a - x4a) * (x2a - x4a) + (y2a - y4a) * (y2a - y4a));
        const auto dist34 = std::sqrt((x3a - x4a) * (x3a - x4a) + (y3a - y4a) * (y3a - y4a));
        
        get_triple(x1a, y1a, x2a, y2a, x3a, y3a, dist12, dist23, dist13); // 123
        get_triple(x1a, y1a, x2a, y2a, x4a, y4a, dist12, dist24, dist14); // 124
        get_triple(x1a, y1a, x3a, y3a, x4a, y4a, dist13, dist34, dist14); // 134
        get_triple(x2a, y2a, x3a, y3a, x4a, y4a, dist23, dist34, dist24); // 234
    }
    quad_centers.clear();
    resize_2d_preserve(quad_star_distances, 8, static_cast<std::size_t>(nrquads));
}

void find_quads_xy(const StarList& starlist, StarList& starlistquads) {
    const auto nrstars_min_one = static_cast<int>(starlist[0].size()) - 1;
    
    if (nrstars_min_one < 3) {
        set_length_2d(starlistquads, 10, 0);
        return;
    }
    
    auto nrquads = 0;
    set_length_2d(starlistquads, 10, static_cast<std::size_t>(nrstars_min_one));
    
    auto j_distance1 = 0;
    auto j_distance2 = 0;
    auto j_distance3 = 0;
    
    for (int i = 0; i <= nrstars_min_one; ++i) {
        auto distance1 = 1e99;
        auto distance2 = 1e99;
        auto distance3 = 1e99;
        
        for (int j = 0; j <= nrstars_min_one; ++j) {
            if (j == i) {
                continue;
            }
            const auto dx = starlist[0][j] - starlist[0][i];
            const auto dy = starlist[1][j] - starlist[1][i];
            const auto distance = dx * dx + dy * dy;
            if (distance < distance1) {
                distance3 = distance2; j_distance3 = j_distance2;
                distance2 = distance1; j_distance2 = j_distance1;
                distance1 = distance;  j_distance1 = j;
            }
            else if (distance < distance2) {
                distance3 = distance2; j_distance3 = j_distance2;
                distance2 = distance;  j_distance2 = j;
            }
            else if (distance < distance3) {
                distance3 = distance;  j_distance3 = j;
            }
        }
        
        const auto x1 = starlist[0][i];
        const auto y1 = starlist[1][i];
        const auto x2 = starlist[0][j_distance1];
        const auto y2 = starlist[1][j_distance1];
        const auto x3 = starlist[0][j_distance2];
        const auto y3 = starlist[1][j_distance2];
        const auto x4 = starlist[0][j_distance3];
        const auto y4 = starlist[1][j_distance3];
        
        const auto xt = (x1 + x2 + x3 + x4) / 4.0;
        const auto yt = (y1 + y2 + y3 + y4) / 4.0;
        
        auto identical = false;
        for (int k = 0; k < nrquads; ++k) {
            if (std::fabs(xt - starlistquads[8][k]) < 1.0 &&
                std::fabs(yt - starlistquads[9][k]) < 1.0) {
                identical = true;
                break;
            }
        }
        if (identical) {
            continue;
        }
        
        starlistquads[0][nrquads] = x1;
        starlistquads[1][nrquads] = y1;
        starlistquads[2][nrquads] = x2;
        starlistquads[3][nrquads] = y2;
        starlistquads[4][nrquads] = x3;
        starlistquads[5][nrquads] = y3;
        starlistquads[6][nrquads] = x4;
        starlistquads[7][nrquads] = y4;
        starlistquads[8][nrquads] = xt;
        starlistquads[9][nrquads] = yt;
        ++nrquads;
    }
    resize_2d_preserve(starlistquads, 10, static_cast<std::size_t>(nrquads));
}

void find_stars(const ImageArray& img, double hfd_min, double hfd_max,
                int max_stars,
                const Background& bck, StarList& starlist1) {
    constexpr auto buffersize = 5000;

    if (img.empty() || img[0].empty() || img[0][0].empty()) {
        set_length_2d(starlist1, 2, 0);
        return;
    }

    const auto width2  = static_cast<int>(img[0][0].size());
    const auto height2 = static_cast<int>(img[0].size());

    set_length_2d(starlist1, 2, buffersize);
    auto snr_list = std::vector<double>(buffersize, 0.0);

    ImageArray img_sa;
    img_sa.assign(1, std::vector<std::vector<float>>(
        height2, std::vector<float>(width2, -1.0f)));

    auto retries     = 3;
    auto nrstars     = 0;
    auto highest_snr = 0.0;

    // HFD search radius (aperture-growth cap). 14 matches the Pascal original
    // and handles star HFDs up to ~10 comfortably. Extended sources with
    // HFD > 14 land in the aperture_cap exit and get rejected — which is
    // what we want (they are galaxies or nebulosity, not stars).
    constexpr auto rs = 14;

    astap::core::HfdScratch scratch{};

    do {
        auto detection_level = 0.0;
        if (retries == 3) {
            if (bck.star_level > 30 * bck.noise_level) {
                detection_level = bck.star_level;
            } else {
                retries = 2;
            }
        }
        if (retries == 2) {
            if (bck.star_level2 > 30 * bck.noise_level) {
                detection_level = bck.star_level2;
            } else {
                retries = 1;
            }
        }
        if (retries == 1) {
            detection_level = 30 * bck.noise_level;
        }
        if (retries == 0) {
            detection_level = 7 * bck.noise_level;
        }
        if (retries == -1) {
            // Last-resort floor for low-contrast images (e.g. DSS2 photographic
            // plates) where 7*noise still rejects genuine faint stars. SNR > 10
            // and HFD gates remain the final filters.
            detection_level = 3 * bck.noise_level;
        }

        highest_snr = 0.0;
        nrstars     = 0;

        for (auto& row : img_sa[0]) {
            std::fill(row.begin(), row.end(), -1.0f);
        }

        for (int fitsY = 0; fitsY < height2 - 1; ++fitsY) {
            for (int fitsX = 0; fitsX < width2 - 1; ++fitsX) {
                if (img_sa[0][fitsY][fitsX] > 0.0f) {
                    continue;
                }
                if ((img[0][fitsY][fitsX] - bck.backgr) <= detection_level) {
                    continue;
                }

                astap::core::HfdResult r;
                astap::core::HFD(img, fitsX, fitsY, rs,
                    /*aperture_small=*/99.0, /*adu_e=*/0.0,
                    /*xbinning=*/1.0, r, scratch);
                if (r.hfd <= hfd_max && r.snr > 10.0 && r.hfd > hfd_min) {
                    // 1.5*HFD disk is wide enough to mask a star's own core
                    // (preventing re-detection) but narrow enough not to
                    // suppress a legitimate neighbour. Pascal uses 3*HFD
                    // which merges distinct stars ~10px apart.
                    const auto radius     = static_cast<int>(std::lround(1.5 * r.hfd));
                    const auto sqr_radius = radius * radius;
                    const auto xci        = static_cast<int>(std::lround(r.xc));
                    const auto yci        = static_cast<int>(std::lround(r.yc));

                    // If HFD's refined centroid landed inside an existing
                    // mask, we just re-detected a source we already recorded
                    // (typical for wide extended sources like galaxies where
                    // a halo pixel outside the first mask triggers HFD and
                    // its centroid drifts back to the galaxy core). Skip.
                    if (xci >= 0 && xci < width2 && yci >= 0 && yci < height2
                            && img_sa[0][yci][xci] > 0.0f) {
                        continue;
                    }

                    for (int n = -radius; n <= radius; ++n) {
                        for (int m = -radius; m <= radius; ++m) {
                            const auto j = n + yci;
                            const auto i = m + xci;
                            if (j >= 0 && i >= 0 && j < height2 && i < width2
                                    && (m * m + n * n) <= sqr_radius) {
                                img_sa[0][j][i] = 1.0f;
                            }
                        }
                    }

                    ++nrstars;
                    if (nrstars >= static_cast<int>(starlist1[0].size())) {
                        resize_2d_preserve(starlist1, 2,
                            static_cast<std::size_t>(nrstars + buffersize));
                        snr_list.resize(
                            static_cast<std::size_t>(nrstars + buffersize), 0.0);
                    }
                    starlist1[0][nrstars - 1] = r.xc;
                    starlist1[1][nrstars - 1] = r.yc;
                    snr_list[nrstars - 1]     = r.snr;
                    if (r.snr > highest_snr) {
                        highest_snr = r.snr;
                    }
                }
            }
        }

        --retries;
    } while (!(nrstars >= max_stars || retries < -1));

    img_sa.clear();

    resize_2d_preserve(starlist1, 2, static_cast<std::size_t>(nrstars));
    snr_list.resize(static_cast<std::size_t>(nrstars));

    if (nrstars > max_stars) {
        get_brightest_stars(max_stars, highest_snr, snr_list, starlist1);
    }
}

bool find_offset_and_rotation(int minimum_quads, double tolerance) {
    tolerance = std::min(tolerance, 0.008);  // clamp
    
    // Require >= 3 quads (3 centroid correspondences).
    if (!find_fit(minimum_quads, tolerance)) {
        reset_solution_vectors(0.001);  // nullify
        return false;
    }
    
    // Solve b_X = A * solution_vectorX and b_Y = A * solution_vectorY.
    if (!lsq_fit(A_XYpositions, b_Xrefpositions, solution_vectorX)) {
        reset_solution_vectors(0.001);
        return false;
    }
    if (!lsq_fit(A_XYpositions, b_Yrefpositions, solution_vectorY)) {
        reset_solution_vectors(0.001);
        return false;
    }
    
    // Sanity check: the two axes should have roughly equal scale.
    const auto xy_sqr_ratio =
        (solution_vectorX[0] * solution_vectorX[0] +
         solution_vectorX[1] * solution_vectorX[1]) /
        (1e-8 + solution_vectorY[0] * solution_vectorY[0] +
                solution_vectorY[1] * solution_vectorY[1]);
    if (xy_sqr_ratio < 0.9 || xy_sqr_ratio > 1.1) {
        reset_solution_vectors(0.001);
        // TODO: external dep - memo2_message('Solution skipped on XY ratio: ...')
        return false;
    }
    return true;
}
 
} // namespace
