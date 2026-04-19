///----------------------------------------
///      @file photometry_catalog.cpp
///   @ingroup ASTAP++
///     @brief Glue that pairs a local HNSKY star catalog with @ref calibrate_flux.
///   @details Provides the port's @c plot_and_measure_stars — the driver the
///            Pascal original uses to wire an image's WCS, a star catalog
///            selection, and the per-star HFD/flux measurement together to
///            produce @c head.mzero. This port currently supports the local
///            .1476 / .290 database path only; wide-field (w08) and online
///            Gaia catalog paths are TODO.
///
///            The star-overlay half of the Pascal's @c plot_and_measure_stars
///            (annotations drawn on the GUI canvas) is intentionally not
///            ported here — that lives in the GUI layer.
///    @author Ported from Han Kleijn's @c unit_annotation.pas (ASTAP); MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "fits.h"
#include "photometry.h"
#include "wcs.h"
#include "../reference/star_database.h"

#include <array>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

///----------------------------------------
namespace astap::core {
///----------------------------------------

namespace {

constexpr auto kPi         = std::numbers::pi;
constexpr auto kDeg2Rad    = kPi / 180.0;
constexpr auto kTile1476   = 5.142857143 * kDeg2Rad;  ///< .1476 cell diameter
constexpr auto kTile290    =  9.53 * kDeg2Rad;        ///< .290 cell diameter

// Record-size cap for fluxpoint calibration — fainter stars carry too much
// noise to contribute usefully. Scaled by image area against the Pascal's
// reference 2328*1760 frame.
constexpr auto kStarsPerReferencePixel = 730.0 / (2328.0 * 1760.0);

///----------------------------------------
/// @brief Decide the passband identifier for a local database name.
/// @details The first character of the database prefix encodes the passband
///          (v/V → Johnson-V, everything else → Gaia BP). Matches the Pascal
///          @c uppercase(copy(name_database,1,1))='V' branch.
///----------------------------------------

[[nodiscard]] std::string passband_for_database(std::string_view name) {
    if (!name.empty() && (name.front() == 'v' || name.front() == 'V')) {
        return "V";
    }
    return "BP";
}

}  // namespace

/// MARK: - plot_and_measure_stars

void plot_and_measure_stars(const ImageArray& img,
                            std::vector<std::string>& memo, Header& head,
                            bool flux_calibration, bool plot_stars,
                            bool report_lim_magnitude) {
    // Guard — no image or no WCS means nothing to do.
    if (img.empty() || img[0].empty() || head.naxis == 0 || head.cd1_1 == 0) {
        return;
    }

    // The star-annotation overlay path (plot_stars=true) is handled by the
    // GUI layer's own render pipeline and is not ported here.
    if (plot_stars) {
        memo.push_back("plot_and_measure_stars: star overlay rendering is provided by the GUI layer.");
    }

    if (!flux_calibration) {
        return;
    }

    // Pointing at image centre. FITS pixels are 1-based so the centre of an
    // N-pixel axis is (N+1)/2.
    auto telescope_ra  = 0.0;
    auto telescope_dec = 0.0;
    pixel_to_celestial(head,
                        (head.width  + 1) * 0.5,
                        (head.height + 1) * 0.5,
                        /*formalism=*/1,
                        telescope_ra, telescope_dec);

    // Frame diagonal FOV in radians — used to decide which catalog tiles to open.
    const auto fov_org = std::sqrt(
        (head.width  * head.cdelt1) * (head.width  * head.cdelt1) +
        (head.height * head.cdelt2) * (head.height * head.cdelt2)) * kDeg2Rad;

    // Target star count. Cap at the brightest set to keep noise down.
    auto max_nr_stars = static_cast<int>(std::lround(
        static_cast<double>(head.width) * head.height * kStarsPerReferencePixel));

    // Select a local catalog sized to the frame. The auto-selection inside
    // select_star_database considers database_path's installed DBs and the FOV.
    if (!astap::reference::select_star_database("",
            head.height * std::abs(head.cdelt2))) {
        memo.push_back("No local star database available for photometric calibration.");
        return;
    }
    memo.push_back("Using star database " + astap::reference::name_database);

    // The wide-field (type 1) and online-Gaia (type 0) paths are not yet ported
    // in this translation unit; only the tile-based .1476 / .290 formats.
    if (astap::reference::database_type != 1476 &&
        astap::reference::database_type != 290) {
        memo.push_back("Photometric calibration currently requires a .1476 or .290 "
                       "database; wide-field and online catalogs are not yet supported.");
        return;
    }

    const auto passband = passband_for_database(astap::reference::name_database);

    // FOV clamp: a single database read pass must stay inside a tile so the
    // reader never skips over a cell.
    const auto tile_cap = (astap::reference::database_type == 1476)
                             ? kTile1476 : kTile290;
    const auto fov = std::min(fov_org, tile_cap);
    if (fov_org > fov) {
        max_nr_stars = static_cast<int>(
            max_nr_stars * (fov * fov) / (fov_org * fov_org));
    }

    // Identify up to four cell files covering the FOV.
    auto area1 = 0, area2 = 0, area3 = 0, area4 = 0;
    auto frac1 = 0.0, frac2 = 0.0, frac3 = 0.0, frac4 = 0.0;
    astap::reference::find_areas(telescope_ra, telescope_dec, fov,
                                  area1, area2, area3, area4,
                                  frac1, frac2, frac3, frac4);

    const auto areas = std::array<int, 4>{area1, area2, area3, area4};
    const auto fracs = std::array<double, 4>{frac1, frac2, frac3, frac4};

    // Build a StarSource that lazily walks the four tiles. The cumulative-
    // fraction cap mirrors the Pascal's per-area `nrstars_required2` check.
    struct CatalogIterator {
        std::array<int, 4>    areas;
        std::array<double, 4> fracs;
        double telescope_ra;
        double telescope_dec;
        double fov;
        int    max_nr_stars;
        int    star_total_counter{0};
        int    area_index{0};
        double accumulated_frac{0.0};
        bool   area_open{false};

        bool next(CatalogStar& out) {
            while (true) {
                // Advance to the next area when the current one is exhausted or
                // past its allotted share of stars.
                while (area_index < 4) {
                    const auto current_area = areas[area_index];
                    if (current_area == 0) {
                        ++area_index;
                        continue;
                    }
                    const auto cumulative_cap =
                        accumulated_frac + fracs[area_index];
                    if (star_total_counter >=
                        static_cast<int>(max_nr_stars * cumulative_cap)) {
                        accumulated_frac = cumulative_cap;
                        ++area_index;
                        if (area_open) {
                            astap::reference::close_star_database();
                            area_open = false;
                        }
                        continue;
                    }
                    break;
                }
                if (area_index >= 4) {
                    if (area_open) {
                        astap::reference::close_star_database();
                        area_open = false;
                    }
                    return false;
                }

                // Open the area's cell file on first touch.
                if (!area_open) {
                    if (!astap::reference::open_database(telescope_dec,
                                                         areas[area_index])) {
                        ++area_index;
                        continue;
                    }
                    area_open = true;
                }

                // Pull the next star out of the current cell.
                auto ra  = 0.0;
                auto dec = 0.0;
                auto mag = 0.0;
                auto bp_rp = 999.0;
                if (astap::reference::readdatabase290(
                        telescope_ra, telescope_dec, fov,
                        ra, dec, mag, bp_rp)) {
                    ++star_total_counter;
                    out.ra    = ra;
                    out.dec   = dec;
                    out.magn  = mag;
                    out.Bp_Rp = bp_rp;
                    return true;
                }

                // End of this area file — advance and keep walking.
                accumulated_frac += fracs[area_index];
                ++area_index;
                astap::reference::close_star_database();
                area_open = false;
            }
        }
    };

    auto iter = std::make_shared<CatalogIterator>(CatalogIterator{
        .areas = areas,
        .fracs = fracs,
        .telescope_ra  = telescope_ra,
        .telescope_dec = telescope_dec,
        .fov = fov,
        .max_nr_stars  = max_nr_stars,
    });

    // Extended-objects mode for SQM: aperture_setting=0 means the limiting-
    // magnitude calc uses a very large virtual aperture. Callers that want
    // point-source aperture calibration set these before calling.
    head.mzero_radius = 99.0;
    constexpr auto kAnnulusRadius = 14;
    constexpr auto kAperture = 0.0;

    [[maybe_unused]] const auto result = calibrate_flux(
        img, head, memo,
        [iter](CatalogStar& out) { return iter->next(out); },
        passband, kAnnulusRadius, kAperture, report_lim_magnitude);

    // Close any cell that stayed open (e.g. calibrate_flux stopped early
    // before the StarSource reported exhaustion).
    astap::reference::close_star_database();
}

/// MARK: - non-const-ref overload
///
/// A few callers (see @c unit_sqm.pas line 69 in the original) pass a mutable
/// @c image_array even though the function does not modify it. Provide the
/// overload so the linker resolves both call sites.

void plot_and_measure_stars(ImageArray& img,
                            std::vector<std::string>& memo, Header& head,
                            bool flux_calibration, bool plot_stars,
                            bool report_lim_magnitude) {
    plot_and_measure_stars(static_cast<const ImageArray&>(img), memo, head,
                           flux_calibration, plot_stars, report_lim_magnitude);
}

} // namespace
