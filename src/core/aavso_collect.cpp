///----------------------------------------
///      @file aavso_collect.cpp
///   @ingroup ASTAP++
///     @brief Implementation of @ref collect_aavso_measurements.
///    @author Created by John Stephen on 4/25/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "aavso_collect.h"

#include "fits.h"
#include "globals.h"
#include "hjd.h"
#include "photometry.h"
#include "wcs.h"
#include "../stacking/stack.h"

#include <cmath>

///----------------------------------------
namespace astap::core {
///----------------------------------------

namespace {

// Run HFD at a star's projected pixel position. Returns true when a real
// star was found (flux > 0) and writes mag/snr through the out-params.
[[nodiscard]] bool measure_star(const ImageArray& img, const Header& head,
                                double ra_rad, double dec_rad,
                                int annulus_rad,
                                double& mag_out, double& snr_out) {
	auto px = 0.0, py = 0.0;
	celestial_to_pixel(head, ra_rad, dec_rad, px, py);
	const auto x = static_cast<int>(std::lround(px - 1.0));
	const auto y = static_cast<int>(std::lround(py - 1.0));
	if (x < 0 || y < 0 || x >= head.width || y >= head.height) {
		return false;
	}

	auto scratch = HfdScratch{};
	auto r       = HfdResult{};
	HFD(img, x, y, annulus_rad, /*aperture_small=*/99.0, /*adu_e=*/0.0,
	    /*xbinning=*/1.0, r, scratch);
	if (r.flux <= 0.0) return false;

	mag_out = head.mzero - 2.5 * std::log10(r.flux);
	snr_out = r.snr;
	return true;
}

}  // namespace

AavsoCollectResult collect_aavso_measurements(
		const std::vector<std::filesystem::path>& files,
		const AavsoCollectOptions& opts) {
	auto result = AavsoCollectResult{};
	result.frames_total = static_cast<int>(files.size());

	auto log = [&](const std::string& msg) { if (opts.log) opts.log(msg); };

	for (const auto& path : files) {
		auto memo = std::vector<std::string>{};
		auto head = Header{};
		auto img  = ImageArray{};

		if (!load_fits(path, /*light=*/true, /*load_data=*/true,
		               /*update_memo=*/true, /*get_ext=*/0,
		               memo, head, img)) {
			log("Skipping (load failed): " + path.string());
			result.skipped_files.push_back(path.string());
			continue;
		}
		if (head.cd1_1 == 0.0) {
			log("Skipping (no plate solution): " + path.string());
			result.skipped_files.push_back(path.string());
			continue;
		}
		if (head.mzero == 0.0) {
			log("Skipping (no MZERO — run Photometric Calibration first): "
			    + path.string());
			result.skipped_files.push_back(path.string());
			continue;
		}

		auto var_mag = 0.0, var_snr = 0.0;
		auto chk_mag = 0.0, chk_snr = 0.0;
		auto cmp_mag = 0.0, cmp_snr = 0.0;

		if (!measure_star(img, head, opts.variable.ra, opts.variable.dec,
		                  opts.annulus_radius, var_mag, var_snr)) {
			log("Skipping (variable not measurable): " + path.string());
			result.skipped_files.push_back(path.string());
			continue;
		}
		if (!measure_star(img, head, opts.check.ra, opts.check.dec,
		                  opts.annulus_radius, chk_mag, chk_snr)) {
			log("Skipping (check not measurable): " + path.string());
			result.skipped_files.push_back(path.string());
			continue;
		}
		if (!opts.ensemble) {
			if (!measure_star(img, head, opts.comp.ra, opts.comp.dec,
			                  opts.annulus_radius, cmp_mag, cmp_snr)) {
				log("Skipping (comp not measurable): " + path.string());
				result.skipped_files.push_back(path.string());
				continue;
			}
		}

		// Per-frame JD / HJD from this frame's own DATE-OBS / DATE-AVG.
		// date_to_jd writes to module globals jd_start/jd_mid/jd_end.
		astap::stacking::date_to_jd(head.date_obs, head.date_avg, head.exposure);
		auto jd = astap::jd_mid;
		if (opts.hjd_date && jd > 2400000.0) {
			jd = JD_to_HJD(jd, head.ra0, head.dec0);
		}

		// Airmass from the FITS header field (string). Best-effort parse;
		// 99 means "na".
		auto airmass = 99.0;
		try {
			if (!head.airmass.empty()) {
				airmass = std::stod(head.airmass);
			}
		} catch (...) { airmass = 99.0; }

		auto m = AavsoMeasurement{};
		m.variable_name   = opts.variable.name;
		m.check_name      = opts.check.name;
		m.var_magnitude   = var_mag;
		m.check_magnitude = chk_mag;
		m.snr             = var_snr;
		m.jd              = jd;
		m.airmass         = airmass;
		m.filter_band     = opts.filter_band;
		if (!opts.ensemble) {
			m.comp_name        = opts.comp.name;
			m.comp_magnitude   = cmp_mag;
			m.comp_catalog_mag = opts.comp_catalog_magnitude;
		}

		result.rows.push_back(std::move(m));
		++result.frames_accepted;
		log("Measured " + path.filename().string() +
		    "  var=" + std::to_string(var_mag) +
		    "  chk=" + std::to_string(chk_mag));
	}

	return result;
}

} // namespace
