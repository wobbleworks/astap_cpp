///----------------------------------------
///      @file globals.h
///   @ingroup ASTAP++
///     @brief Cross-module globals for the ASTAP++ port.
///   @details During the per-unit translations each .cpp forward-declared the
///            globals it touched with a "TODO: consolidate" note.  This file is
///            that consolidation -- the single canonical site for variables the
///            original source kept at astap_main unit-scope and that multiple
///            ported modules now need to see.
///
///            Module-owned state stays in the owning module's header:
///              reference/star_database.h   -> database cache, name_database,
///                                             database_path, database_type
///              reference/stars_wide_field.h-> wide_field_stars, wide_database
///              reference/online_gaia.h     -> online_database, gaia_ra/dec,
///                                             passband_active
///              core/online.h               -> vsp, vsx
///              core/imaging.h              -> histogram, stretch_c,
///                                             saturation_factor...
///              solving/star_align.h        -> solution_vectorX/Y/cblack,
///                                             quad_star_distances1/2,
///                                             A_XYpositions,
///                                             b_Xref/Yrefpositions,
///                                             nr_references(2)
///              stacking/live_stacking.h    -> (drops its atomics; they live
///                                             here now)
///
///            Everything else that crosses module boundaries is declared here
///            and defined once in globals.cpp.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP); MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <array>
#include <atomic>
#include <filesystem>
#include <string>
#include <vector>

#include "../types.h"

///----------------------------------------
namespace astap {
///----------------------------------------

/// MARK: Process control

extern std::atomic<bool> esc_pressed;      ///< @brief User requested abort (Ctrl-C).
extern std::atomic<bool> pause_pressed;    ///< @brief Live-stacking pause flag.
extern std::atomic<bool> live_stacking;    ///< @brief Live-stack loop is active.
extern bool        commandline_execution;  ///< @brief Running headless.
extern bool        commandline_log;        ///< @brief -log specified.
extern int         errorlevel;             ///< @brief Exit code.
extern std::string filename2;              ///< @brief Current image filename.
extern std::string cmdline;                ///< @brief Full argv joined into one string.
extern std::string warning_str;            ///< @brief Cumulative warnings.
extern std::string astap_version;          ///< @brief Version string "YYYY.MM.DD".
extern bool        solve_show_log;         ///< @brief Emit solver progress to memo2.
extern bool        fov_specified;          ///< @brief -fov was supplied.

/// MARK: Current image state

extern Header     head;                    ///< @brief FITS header of img_loaded.
extern Header     head_ref;                ///< @brief Reference-image header (stacking).
extern ImageArray img_loaded;              ///< @brief Current image buffer.
extern int        nrbits;                  ///< @brief BITPIX of last load.
extern int        extend_type;             ///< @brief 0=none, 1=image, 2=ASCII, 3=bin.
extern std::string instrum;                ///< @brief INSTRUME header value.
extern Background bck;                     ///< @brief Last-measured background stats.
extern double     cwhite;                  ///< @brief Display white-point.
extern std::string bayerpat;               ///< @brief Bayer pattern string.
extern std::string object_name;            ///< @brief OBJECT header.
extern bool       unsaved_import;          ///< @brief True if not on disk as FITS.
extern std::string roworder;               ///< @brief Top-down vs bottom-up.
extern int        xbayroff, ybayroff;      ///< @brief Bayer offsets.

/// MARK: Mount / pointing

extern double ra_radians;                  ///< @brief Parsed from RA text field.
extern double dec_radians;                 ///< @brief Parsed from Dec text field.
extern double ra_mount;                    ///< @brief Last reported by mount.
extern double dec_mount;                   ///< @brief Last reported by mount.

/// MARK: SIP distortion coefficients

extern bool sip;                           ///< @brief SIP distortion enabled.
extern int  a_order, b_order;              ///< @brief Forward polynomial order.
extern int  ap_order, bp_order;            ///< @brief Inverse polynomial order.

extern double a_0_0, a_0_1, a_0_2, a_0_3;
extern double a_1_0, a_1_1, a_1_2;
extern double a_2_0, a_2_1;
extern double a_3_0;
extern double b_0_0, b_0_1, b_0_2, b_0_3;
extern double b_1_0, b_1_1, b_1_2;
extern double b_2_0, b_2_1;
extern double b_3_0;

extern double ap_0_0, ap_0_1, ap_0_2, ap_0_3;
extern double ap_1_0, ap_1_1, ap_1_2;
extern double ap_2_0, ap_2_1;
extern double ap_3_0;
extern double bp_0_0, bp_0_1, bp_0_2, bp_0_3;
extern double bp_1_0, bp_1_1, bp_1_2;
extern double bp_2_0, bp_2_1;
extern double bp_3_0;

/// MARK: Solver cached kinematics and outputs

extern double SIN_dec0, COS_dec0;          ///< @brief Cached sin/cos of head.dec0.
extern double SIN_dec_ref, COS_dec_ref;    ///< @brief Same for head_ref.dec0.
extern double x_new_float, y_new_float;    ///< @brief calc_newx_newy out-params.
extern double referenceX, referenceY;      ///< @brief Manual-alignment anchor pixel.
extern double pedestal_s;                  ///< @brief Target background for normalise.
extern double mag2;                        ///< @brief Faintest mag from read_stars.
extern std::string solution_str;           ///< @brief Last solve summary string.

/// MARK: Stacking accumulator state

extern double sum_exp, sum_temp;
extern double airmass_sum, airmass;
extern double jd_sum;
extern double jd_start_first, jd_end_last;
extern double jd_start, jd_mid, jd_end;
extern double jd_mid_reference;
extern int    images_selected;

extern int counterR,   counterRdark,   counterRflat,   counterRbias;
extern int counterG,   counterGdark,   counterGflat,   counterGbias;
extern int counterB,   counterBdark,   counterBflat,   counterBbias;
extern int counterL,   counterLdark,   counterLflat,   counterLbias;
extern int counterRGB, counterRGBdark, counterRGBflat, counterRGBbias;
extern int exposureR, exposureG, exposureB, exposureL, exposureRGB;
extern double temperatureR, temperatureG, temperatureB;
extern double temperatureL, temperatureRGB;

/// MARK: Solver settings

extern int    max_stars_setting;           ///< @brief Max stars for solving.
extern int    downsample_setting;          ///< @brief Downsample for solving (0 = auto).
extern double quad_tolerance;              ///< @brief Quad match tolerance.
extern double min_star_size_arcsec;        ///< @brief Minimum star size in arcseconds.
extern double search_radius_deg;           ///< @brief Search radius in degrees.
extern double search_fov_deg;              ///< @brief Search FOV in degrees (0 = auto).
extern bool   force_oversize;              ///< @brief Force oversize solving.
extern bool   add_sip;                     ///< @brief Add SIP distortion to solution.
extern bool   check_pattern_filter;        ///< @brief Pattern filter check enabled.

/// MARK: Stacking settings

extern bool   use_manual_align;            ///< @brief Align via manual-marked reference star.
extern bool   use_ephemeris_alignment;     ///< @brief Align via ephemeris (comet/asteroid).
extern bool   use_astrometry_internal;     ///< @brief Align via per-frame plate solve.
extern bool   skip_alignment;              ///< @brief Use identity transform per frame (assumes pre-aligned inputs).
extern double hfd_min_setting;             ///< @brief Minimum star HFD for stacking (0 = auto).
extern double hfd_max_setting;             ///< @brief Maximum star HFD for stacking (default 10; raise for DSS2 plates).
extern double sigma_clip_factor;           ///< @brief Stdev factor for sigma-clip rejection.

/// MARK: Memos

extern std::vector<std::string> memo1_lines;   ///< @brief FITS header memo.
extern std::vector<std::string> memo2_lines;   ///< @brief Solver log memo.

/// MARK: Filesystem / executable environment

extern std::filesystem::path application_path;    ///< @brief Directory of running binary.
extern std::vector<std::string> recent_files;
extern int raw_conversion_program_index;          ///< @brief 0=libraw -f, 1=libraw -i, 2=dcraw.
extern std::vector<std::string> head1;            ///< @brief Canonical FITS header template.
	
} // namespace
