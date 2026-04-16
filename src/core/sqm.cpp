///----------------------------------------
///      @file sqm.cpp
///   @ingroup ASTAP++
///     @brief Sky Quality Measurement implementation.
///    @author Ported from Han Kleijn's unit_sqm.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "sqm.h"

#include <cmath>
#include <string>

#include "globals.h"
#include "imaging.h"
#include "photometry.h"
#include "wcs.h"
#include "../types.h"

///----------------------------------------
namespace astap::core {
///----------------------------------------

/// MARK: Globals

std::string centalt;
double      altitudefloat = 0.0;
double      sqmfloat      = 0.0;
double      airmass       = 0.0;

std::string sitelat;
std::string sitelong;
std::string lat_default;
std::string long_default;
double      temperature_c = 10.0;
double      pressure_hpa  = 1010.0;

/// MARK: External Dependencies

// TODO(astap-port): lives in the stacking module.
void bin_X2X3X4(ImageArray& img, Header& head,
                std::vector<std::string>& memo, int factor);
void backup_img();
void restore_img();

// TODO(astap-port): lives in the main application module.
extern double annulus_radius;
void plot_and_measure_stars(ImageArray& img,
                            std::vector<std::string>& memo,
                            Header& head,
                            bool calibration,
                            bool plot_stars,
                            bool report_lim_magnitude);
                            
void calculate_az_alt(int mode, Header& head, double& az, double& alt);
double atmospheric_absorption(double airmass);

void memo2_message(const std::string& s);

extern std::string bayerpat;

/// MARK: - SQM Computation

bool calculate_sqm(Header& headx, bool get_bk, bool get_his, int& pedestal2) {
    auto backup_made = false;
    
    // form_exist was form_sqm1 <> nil. The CLI path has no form.
    constexpr auto form_exist = false;
    
    [[maybe_unused]] auto bayer = (!bayerpat.empty() && headx.xbinning == 1);
    
    if constexpr (form_exist) {
        // Intentionally unreachable in the CLI port.
    }
    
    // Re-run photometric calibration for extended sources if MZERO was
    // either missing or calibrated for point sources (mzero_radius != 99).
    if (headx.mzero == 0.0 || headx.mzero_radius != 99.0) {
        annulus_radius = 14.0;
        headx.mzero_radius = 99.0;
        plot_and_measure_stars(img_loaded, memo1_lines, headx,
                               /*calibration=*/true,
                               /*plot_stars=*/false,
                               /*report_lim_magnitude=*/false);
    }
    
    auto result = false;
    
    if (headx.mzero > 0.0) {
        if (get_bk) {
            get_background(/*colour=*/0, img_loaded,
                           /*calc_hist=*/get_his,
                           /*calc_noise_level=*/false,
                           bck);
        }
        
        if (headx.calstat.find('D') != std::string::npos) {
            if (pedestal2 > 0) {
                memo2_message("Dark already applied! Pedestal value ignored.");
                pedestal2 = 0;
            }
        } else if (pedestal2 == 0) {
            memo2_message("Pedestal value missing!");
            warning_str += "Pedestal value missing!";
        }
        
        if (static_cast<double>(pedestal2) >= bck.backgr) {
            memo2_message("Too high pedestal value!");
            warning_str += "Too high pedestal value!";
            pedestal2 = 0;
        }
        
        // Flux per arcsec^2 -> mag/arcsec^2
        auto pixel_arcsec = headx.cdelt2 * 3600.0;
        auto flux_per_arcsec2 =
            (bck.backgr - static_cast<double>(pedestal2) - headx.pedestal)
            / (pixel_arcsec * pixel_arcsec);
        sqmfloat = headx.mzero
                   - std::log(flux_per_arcsec2) * 2.5 / std::log(10.0);
                   
        auto az = 0.0;
        calculate_az_alt(/*mode=*/1, head, az, altitudefloat);
        
        if (altitudefloat > 0.0) {
            auto airm = airmass_calc(altitudefloat);
            airmass = airm;
            
            // Zenith correction: subtract 0.28 from the empirical absorption at airmass = 1
            auto correction = atmospheric_absorption(airm) - 0.28;
            sqmfloat += correction;
            result = true;
        } else {
            memo2_message("Negative altitude calculated!");
            warning_str += "Negative altitude calculated!";
        }
    } else {
        memo2_message("MZERO calibration failure!");
        warning_str += "MZERO calibration failure!";
    }
    
    if (backup_made) {
        restore_img();
        backup_made = false;
    }
    
    return result;
}

/// MARK: - Bortle Classification

std::string bortle(double sqm) noexcept {
    if (sqm > 21.99) {
        return "Bortle 1, excellent dark-sky site";
    }
    if (sqm > 21.89) {
        return "Bortle 2, truly dark site";
    }
    if (sqm > 21.69) {
        return "Bortle 3, dark rural sky";
    }
    if (sqm > 21.25) {
        return "Bortle 4, rural sky";
    }
    if (sqm > 20.49) {
        return "Bortle 4.5, rural/suburban sky";
    }
    if (sqm > 19.50) {
        return "Bortle 5, suburban sky";
    }
    if (sqm > 18.94) {
        return "Bortle 6, bright suburban sky";
    }
    if (sqm > 18.38) {
        return "Bortle=7, suburban/urban sky";
    }
    if (sqm > 17.80) {
        return "Bortle 8, city sky";
    }
    return "Bortle 9, inner-city sky";
}
    
} // namespace
