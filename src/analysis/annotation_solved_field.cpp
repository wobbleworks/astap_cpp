/// @file annotation_solved_field.cpp
/// Solved-field annotation features that read the local star database —
/// ported from unit_annotation.pas.
///
/// Split out from annotation.cpp so the lean font / deep-sky unit tests need
/// not link the star-database and photometry (HFD) subsystems these two
/// features pull in. Declarations live in annotation.h alongside the rest.
///
///   - plot_artificial_stars (unit_annotation.pas:2711) — writes each catalog
///     star as a single magnitude-valued pixel, for supernova / minor-planet
///     reference frames. The original's tile-based database branch was
///     commented out in favour of the online VizieR path; this port restores
///     the local branch and drops the out-of-scope online branch.
///   - measure_distortion (unit_annotation.pas:2465) — the headless residual
///     computation only: HFD-measures each catalog star and reports the median
///     sky<->pixel astrometric error. The GUI distortion-vector overlay and the
///     scale bars are dropped.

#include "annotation.h"

#include "../core/photometry.h"
#include "../core/util.h"
#include "../core/wcs.h"
#include "../reference/star_database.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <vector>

namespace astap::analysis {

// ---------------------------------------------------------------------------
// Local star-catalog iteration (shared by plot_artificial_stars and
// measure_distortion)
// ---------------------------------------------------------------------------

namespace {

constexpr double kDeg2Rad  = std::numbers::pi / 180.0;
constexpr double kTile1476 = 5.142857143 * kDeg2Rad;  // .1476 cell diameter
constexpr double kTile290  = 9.53 * kDeg2Rad;         // .290 cell diameter

/// Visit every local-catalog star within the frame's field of view, matching
/// the solver / photometric-calibration read path (find_areas -> open_database
/// -> readdatabase290), frac-apportioned across up to four database tiles.
///
/// @p visit receives (ra, dec, magn*10) per star and returns true to continue.
/// Iteration also stops once @p max_stars stars have been visited (apportioned
/// per tile by the find_areas fractions; pass a huge value for "no cap").
/// Requires select_star_database to have set database_type to a tile format
/// (1476 or 290); other formats visit nothing. Returns the number of stars
/// visited.
///
/// The FOV is clamped to a single cell diameter so one read pass never skips a
/// tile — the same clamp the original applies to its database branch.
template <class Visit>
int read_local_catalog(double telescope_ra, double telescope_dec,
                       double fov_org, double max_stars, Visit&& visit)
{
	namespace ref = astap::reference;

	const int db_type = ref::database_type;
	if (db_type != 1476 && db_type != 290) return 0;

	const double tile_cap = (db_type == 1476) ? kTile1476 : kTile290;
	const double fov = std::min(fov_org, tile_cap);

	int area1 = 0, area2 = 0, area3 = 0, area4 = 0;
	double frac1 = 0.0, frac2 = 0.0, frac3 = 0.0, frac4 = 0.0;
	ref::find_areas(telescope_ra, telescope_dec, fov,
	                area1, area2, area3, area4, frac1, frac2, frac3, frac4);

	const std::array<int, 4> areas{area1, area2, area3, area4};
	const std::array<double, 4> fracs{frac1, frac2, frac3, frac4};

	int star_total_counter = 0;
	double accumulated_frac = 0.0;

	for (int a = 0; a < 4; ++a) {
		if (areas[a] == 0) { accumulated_frac += fracs[a]; continue; }

		const double cap = max_stars * (accumulated_frac + fracs[a]);
		if (!ref::open_database(telescope_dec, areas[a])) {
			accumulated_frac += fracs[a];
			continue;
		}

		double ra = 0.0, dec = 0.0, mag = 0.0, bp_rp = 999.0;
		bool stopped = false;
		while (star_total_counter < cap &&
		       ref::readdatabase290(telescope_ra, telescope_dec, fov,
		                            ra, dec, mag, bp_rp)) {
			++star_total_counter;
			if (!visit(ra, dec, mag)) { stopped = true; break; }
		}
		ref::close_star_database();
		accumulated_frac += fracs[a];
		if (stopped) break;
	}

	return star_total_counter;
}

}  // anonymous namespace

// ---------------------------------------------------------------------------
// plot_artificial_stars
// ---------------------------------------------------------------------------

ArtificialStarsResult plot_artificial_stars(ImageArray& img, const Header& head,
                                            const std::string& star_database)
{
	ArtificialStarsResult result;

	// Need a single-plane buffer and an astrometric solution.
	if (img.empty() || img[0].empty() || img[0][0].empty()
	    || head.naxis == 0 || head.cd1_1 == 0.0)
		return result;

	// Field centre from the image centre (FITS pixels are 1-based, so the centre
	// of an N-pixel axis is (N+1)/2). Robust to an off-centre CRPIX.
	double telescope_ra = 0.0;
	double telescope_dec = 0.0;
	astap::core::pixel_to_celestial(head, (head.width + 1) / 2.0,
	                                (head.height + 1) / 2.0, 0,
	                                telescope_ra, telescope_dec);

	// Frame diagonal FOV in radians, 0% extra (Pascal fov_org).
	const double fov_org = std::sqrt(
	    (head.width * head.cdelt1) * (head.width * head.cdelt1) +
	    (head.height * head.cdelt2) * (head.height * head.cdelt2)) * kDeg2Rad;

	// Select a local catalog sized to the frame height, exactly like the
	// photometric-calibration path.
	if (!astap::reference::select_star_database(star_database,
	        head.height * std::abs(head.cdelt2)))
		return result;

	// Go one magnitude fainter than the solve's limit; the reader reports
	// magnitude*10, so compare against (magn_limit + 1) * 10 (Pascal m_limit*10).
	// head.magn_limit is populated only by photometric calibration, not by a
	// plain solve; when it is uncalibrated (<= 0) there is no meaningful faint
	// cut, so the whole (already magnitude-bounded) local catalog is plotted.
	const bool apply_limit = head.magn_limit > 0.0;
	const double m_limit_x10 = (head.magn_limit + 1.0) * 10.0;

	const int width = head.width;
	const int height = head.height;
	int plotted = 0;

	// No star budget — plot every catalog star in the field (Pascal plots the
	// whole online result set). The original's tile branch stops at the first
	// star fainter than the limit; because tile stars are DEC-clustered rather
	// than globally magnitude-sorted, this port instead skips fainter stars and
	// keeps reading, which is what the online branch effectively does.
	(void)read_local_catalog(
	    telescope_ra, telescope_dec, fov_org,
	    std::numeric_limits<double>::max(),
	    [&](double ra2, double dec2, double mag2) {
		    if (apply_limit && mag2 > m_limit_x10) return true;  // fainter than limit: skip

		    double fits_x = 0.0;
		    double fits_y = 0.0;
		    astap::core::celestial_to_pixel(head, ra2, dec2, fits_x, fits_y);
		    const int x = static_cast<int>(std::lround(fits_x - 1.0));
		    const int y = static_cast<int>(std::lround(fits_y - 1.0));

		    if (x >= 0 && x <= width - 1 && y >= 0 && y <= height - 1) {
			    float& pixel = img[0][y][x];
			    if (pixel >= 1000.0f)
				    pixel = static_cast<float>(mag2);  // empty: place the star
			    else
				    // Flux-combine two overlapping stars' magnitudes.
				    pixel = static_cast<float>(
				        -2.5 * std::log10(std::pow(10.0, -0.4 * pixel) +
				                          std::pow(10.0, -0.4 * mag2)));
			    ++plotted;
		    }
		    return true;
	    });

	result.star_count = plotted;
	result.ok = plotted > 0;
	return result;
}

// ---------------------------------------------------------------------------
// measure_distortion
// ---------------------------------------------------------------------------

DistortionResult measure_distortion(const ImageArray& img, const Header& head,
                                    int formalism, const std::string& star_database)
{
	DistortionResult result;

	if (img.empty() || img[0].empty() || img[0][0].empty()
	    || head.naxis == 0 || head.cd1_1 == 0.0)
		return result;

	// Field centre from the image centre (works even for an off-centre CRPIX).
	double telescope_ra = 0.0;
	double telescope_dec = 0.0;
	astap::core::pixel_to_celestial(head, (head.width + 1) / 2.0,
	                                (head.height + 1) / 2.0, 0,
	                                telescope_ra, telescope_dec);

	// Star budget scaled from Pascal's 1216-in-2328x1760 reference frame.
	const int max_nr_stars = static_cast<int>(std::lround(
	    static_cast<double>(head.width) * head.height *
	    (1216.0 / (2328.0 * 1760.0))));
	if (max_nr_stars <= 0) return result;

	// Residual arrays, split by radius from CRPIX: inner (within 0.25*height),
	// outer (beyond 0.5*height); sky->pixel (pixel separation) and pixel->sky
	// (angular separation of the re-projected centroid).
	std::vector<double> errors_sky_pixel1(max_nr_stars, 0.0);
	std::vector<double> errors_sky_pixel2(max_nr_stars, 0.0);
	std::vector<double> errors_pixel_sky1(max_nr_stars, 0.0);
	std::vector<double> errors_pixel_sky2(max_nr_stars, 0.0);

	if (!astap::reference::select_star_database(star_database,
	        head.height * std::abs(head.cdelt2)))
		return result;

	result.distortion_data.reserve(max_nr_stars);

	// Diagonal FOV, 0% extra (Pascal fov_org).
	const double fov_org = std::sqrt(
	    (head.width * head.cdelt1) * (head.width * head.cdelt1) +
	    (head.height * head.cdelt2) * (head.height * head.cdelt2)) * kDeg2Rad;

	astap::core::HfdScratch scratch;

	const int width = head.width;
	const int height = head.height;
	const double crpix1 = head.crpix1;
	const double crpix2 = head.crpix2;
	const double r_inner2 = (0.25 * height) * (0.25 * height);
	const double r_outer2 = (0.5 * height) * (0.5 * height);

	read_local_catalog(
	    telescope_ra, telescope_dec, fov_org,
	    static_cast<double>(max_nr_stars),
	    [&](double ra2, double dec2, double /*mag2*/) {
		    // Project the catalog position onto the image (FITS 1-based -> 0-based).
		    double fits_x = 0.0;
		    double fits_y = 0.0;
		    astap::core::celestial_to_pixel(head, ra2, dec2, fits_x, fits_y);
		    const double x = fits_x - 1.0;
		    const double y = fits_y - 1.0;

		    // Within the frame plus a 50-pixel overlap.
		    if (!(x > -50 && x <= width + 50 && y > -50 && y <= height + 50))
			    return true;

		    astap::core::HfdResult hfd;
		    astap::core::HFD(img,
		                     static_cast<int>(std::lround(x)),
		                     static_cast<int>(std::lround(y)),
		                     14 /*annulus_radius*/, 99 /*flux_aperture*/,
		                     0 /*adu_e*/, head.xbinning, hfd, scratch);

		    // Star detected in the image?
		    if (!(hfd.hfd < 15 && hfd.hfd >= 0.8 && hfd.snr > 10))
			    return true;

		    const double dx = x - crpix1;
		    const double dy = y - crpix2;
		    const double r2 = dx * dx + dy * dy;

		    // Inner band residuals.
		    if (r2 < r_inner2 && result.stars_inner < max_nr_stars) {
			    errors_sky_pixel1[result.stars_inner] =
			        std::sqrt((x - hfd.xc) * (x - hfd.xc) +
			                  (y - hfd.yc) * (y - hfd.yc));
			    ++result.stars_inner;

			    if (result.stars_inner_ps < max_nr_stars) {
				    double ra3 = 0.0, dec3 = 0.0;
				    astap::core::pixel_to_celestial(head, hfd.xc + 1, hfd.yc + 1,
				                                    formalism, ra3, dec3);
				    double sep = 0.0;
				    astap::core::ang_sep(ra3, dec3, ra2, dec2, sep);
				    errors_pixel_sky1[result.stars_inner_ps] = sep;
				    ++result.stars_inner_ps;
			    }
		    }
		    // Outer band residuals.
		    if (r2 > r_outer2 && result.stars_outer < max_nr_stars) {
			    errors_sky_pixel2[result.stars_outer] =
			        std::sqrt((x - hfd.xc) * (x - hfd.xc) +
			                  (y - hfd.yc) * (y - hfd.yc));
			    ++result.stars_outer;

			    if (result.stars_outer_ps < max_nr_stars) {
				    double ra3 = 0.0, dec3 = 0.0;
				    astap::core::pixel_to_celestial(head, hfd.xc + 1, hfd.yc + 1,
				                                    formalism, ra3, dec3);
				    double sep = 0.0;
				    astap::core::ang_sep(ra3, dec3, ra2, dec2, sep);
				    errors_pixel_sky2[result.stars_outer_ps] = sep;
				    ++result.stars_outer_ps;
			    }
		    }

		    if (result.stars_measured < max_nr_stars) {
			    result.distortion_data.push_back(
			        DistortionSample{x, y, hfd.xc, hfd.yc});
			    ++result.stars_measured;
		    }
		    return true;
	    });

	// Median residuals. The pixel->sky counters (stars_inner_ps / _outer_ps)
	// track the sky->pixel counters in lockstep because their arrays are the
	// same length, so the medians match the original's smedian(..., sub_counter).
	result.sky_pixel_error_inner =
	    astap::core::smedian(errors_sky_pixel1, result.stars_inner);
	result.sky_pixel_error_outer =
	    astap::core::smedian(errors_sky_pixel2, result.stars_outer);
	result.pixel_sky_error_inner =
	    astap::core::smedian(errors_pixel_sky1, result.stars_inner_ps);
	result.pixel_sky_error_outer =
	    astap::core::smedian(errors_pixel_sky2, result.stars_outer_ps);

	return result;
}

}  // namespace astap::analysis
