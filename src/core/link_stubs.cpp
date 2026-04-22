///----------------------------------------
///      @file link_stubs.cpp
///   @ingroup ASTAP++
///     @brief Link-time stubs for symbols the port references but doesn't
///            yet define with the exact signature the caller wants.
///   @details Lets the CLI binary link while pieces of the port are still
///            evolving (signature mismatches between caller expectations
///            and the real definition, plus GUI helpers that stay stubbed
///            in the headless build). None of these are intended to do
///            anything useful at runtime — the call sites are all gated
///            behind TODO markers in the stacking/solving pipelines.
///    @author Created by John Stephen on 4/15/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "../types.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace astap::core {

// -----------------------------------------------------------------------------
// Module-level externs declared by hjd.cpp / astrometric_solving.cpp / etc.
// These would normally be populated by the settings loader.
// -----------------------------------------------------------------------------

// sqm.cpp defines: centalt, sitelat, sitelong, lat_default, long_default.
// Canonicals in astap:: (globals.cpp): bayerpat, pressure, focus_temp,
// site_lat_radians, site_long_radians. Only the truly-orphan symbols remain
// here as astap::core placeholders.
std::string centaz;
double      wtime2actual      = 0.0;

std::vector<std::string>  recent_files_local;
std::vector<std::string>* mainwindow_memo1_lines = nullptr;
std::vector<std::string>* memox                  = nullptr;

// -----------------------------------------------------------------------------
// Signature-mismatch bridges. Each of these has a "real" implementation
// elsewhere in the port, but one or more call sites expect a different
// signature than what the canonical function provides. Rather than rewrite
// the call sites (which would ripple), we provide an adapter that ignores
// the caller's arguments and returns a sentinel.
// -----------------------------------------------------------------------------

void HSV2RGB(float, float, float, float& r, float& g, float& b) { r = g = b = 0.0f; }
void RGB2HSV(float, float, float, float& h, float& s, float& v) { h = s = v = 0.0f; }

// fits.cpp loader expects this `path&`-based signature; the image_io.cpp
// caller wants the `string` variant. Provide the string-taking overload.
bool load_fits(const std::filesystem::path& /*filen*/,
               bool /*light*/,
               bool /*load_data*/,
               int  /*get_ext*/,
               std::vector<std::string>& /*memo*/,
               Header& /*head*/,
               ImageArray& /*img*/) {
    return false;
}

bool save_fits(const ImageArray& /*img*/,
               const std::vector<std::string>& /*memo*/,
               const std::filesystem::path& /*filen*/,
               int /*type1*/,
               bool /*override2*/) {
    return false;
}

void update_integer(std::vector<std::string>& /*memo*/,
                    const std::string& /*key*/,
                    const std::string& /*comment*/,
                    int /*value*/) {}
void add_text(std::vector<std::string>& /*memo*/,
              const std::string& /*key*/,
              const std::string& /*value*/) {}

void progress_indicator(double /*percent*/, const std::string& /*info*/) {}

Background get_hist_bck() { return {}; }
void get_hist(int /*colour*/, const ImageArray& /*img*/,
              int (&/*out*/)[3][65536],
              int (&/*mean*/)[3],
              int /*nrbits*/) {}

void remove_solution(bool /*keep_wcs*/) {}
void show_error_label(std::string_view /*text*/) {}
void update_recent_file_menu() {}

int get_demosaic_pattern() { return 0; }

// -----------------------------------------------------------------------------
// Unported Pascal units — supplied as no-op stubs so the linker is happy.
// -----------------------------------------------------------------------------

// plot_and_measure_stars is ported in core/photometry_catalog.cpp.

void dsspos(double, double, double& ra, double& dec)  { ra = 0; dec = 0; }
void EQU_GAL(double, double, double& l, double& b)    { l  = 0; b   = 0; }

// -----------------------------------------------------------------------------
// Smedian with capital S — used by HFD and the stacking pipeline. Ported
// from astap_main.pas. Sorts @p list in place (caller tolerates that) and
// returns the 3-value average around the middle for odd n > 3, the 2-value
// average for even n, and the single element for n <= 1.
// -----------------------------------------------------------------------------
double Smedian(std::vector<double>& list, int len) {
    if (len <= 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (len == 1) {
        return list[0];
    }
    std::sort(list.begin(), list.begin() + len);
    const auto mid = (len - 1) / 2;
    if ((len % 2) == 1) {
        if (len <= 3) {
            return list[mid];
        }
        return (list[mid - 1] + list[mid] + list[mid + 1]) / 3.0;
    }
    return (list[mid] + list[mid + 1]) / 2.0;
}

// -----------------------------------------------------------------------------
// Additional signature-mismatch bridges discovered at link time.
// -----------------------------------------------------------------------------

void restore_img() {}

void update_text(std::vector<std::string>& /*memo*/,
                 const std::string& /*key*/,
                 const std::string& /*value*/) {}

void update_float(std::vector<std::string>& /*memo*/,
                  const std::string& /*key*/,
                  const std::string& /*comment*/,
                  bool /*exponential*/,
                  double /*value*/) {}

// update_integer overload taking string_view keys + long long (used by PPM/PGM
// loader). The real `update_integer` is in fits.cpp with `const std::string&`.
void update_integer(std::vector<std::string>& /*memo*/,
                    std::string_view /*key*/,
                    std::string_view /*comment*/,
                    long long /*value*/) {}

// analyse_image in astap::core (differs from the astap::stacking version).
void analyse_image(const ImageArray& /*img*/,
                   Header& /*head*/,
                   int /*report_type*/,
                   int /*max_stars*/,
                   int& star_counter,
                   Background& bck,
                   double& hfd_median) {
    star_counter = 0;
    bck = Background{};
    hfd_median = 0.0;
}

// memo2_message: two call-signatures exist in the port. Provide both.
void memo2_message(std::string_view /*msg*/) {}
void memo2_message(const std::string& /*msg*/) {}

// raster_rotate overload with (angle, cx, cy, img) — real one is in
// astap::stacking with the same signature but different namespace.
void raster_rotate(double /*angle*/, double /*cx*/, double /*cy*/,
                   ImageArray& /*img*/) {}

double annulus_radius = 4.0;  // Pascal default; sqm.cpp reads this.

// Position angle of a second point as seen from a first (Meeus 48.5).
double position_angle(double /*ra1*/, double /*dec1*/,
                      double /*ra0*/, double /*dec0*/) {
    return 0.0;
}

// Unpacker for .fz compressed FITS — in the real port this shells out to
// funpack. Stub always says "no".
bool unpack_cfitsio(std::string& /*filename*/) { return false; }

// deltaT (TT - UT1) in days. ~70 s for modern epochs.
double deltaT_calc(double /*jd*/) { return 70.0 / 86400.0; }

// File's JD based on filesystem mtime. Stub returns J2000.
double file_age_jd(const std::filesystem::path& /*p*/) { return 2451545.0; }

// Julian Date -> ISO-ish date string.
std::string jd_to_date(double /*jd*/) { return "2000-01-01T12:00:00"; }

// FITS key remover overload used by image_io.cpp (const std::string&).
void remove_key(std::vector<std::string>& /*memo*/,
                const std::string& /*key*/,
                bool /*all*/) {}

}  // namespace astap::core

// -----------------------------------------------------------------------------
// The anonymous-namespace prepare_ra / prepare_dec in astrometric_solving.cpp
// need TU-local definitions; those stubs live inside that file. Nothing to
// add here for them.
// -----------------------------------------------------------------------------
