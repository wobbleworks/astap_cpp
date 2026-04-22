///----------------------------------------
///      @file globals.cpp
///   @ingroup ASTAP++
///     @brief Definitions for cross-module globals declared in globals.h.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP); MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "globals.h"

///----------------------------------------
namespace astap {
///----------------------------------------

/// MARK: Process control

std::atomic<bool> esc_pressed{false};
std::atomic<bool> pause_pressed{false};
std::atomic<bool> live_stacking{false};
bool        commandline_execution = false;
bool        commandline_log       = false;
int         errorlevel            = 0;
std::string filename2;
std::string cmdline;
std::string warning_str;
std::string astap_version         = "0.0.0";
bool        solve_show_log        = false;
bool        fov_specified         = false;

/// MARK: Current image state

Header     head{};
Header     head_ref{};
ImageArray img_loaded;
int        nrbits          = 0;
int        extend_type     = 0;
std::string instrum;
Background bck{};
double     cwhite          = 65535.0;
std::string bayerpat;
std::string object_name;
bool       unsaved_import  = false;
std::string roworder;
int        xbayroff        = 0;
int        ybayroff        = 0;

/// MARK: Mount / pointing

double ra_radians  = 0.0;
double dec_radians = 0.0;
double ra_mount    = 0.0;
double dec_mount   = 0.0;

/// MARK: SIP distortion coefficients

bool sip       = false;
int  a_order   = 0;
int  b_order   = 0;
int  ap_order  = 0;
int  bp_order  = 0;

double a_0_0 = 0.0, a_0_1 = 0.0, a_0_2 = 0.0, a_0_3 = 0.0;
double a_1_0 = 0.0, a_1_1 = 0.0, a_1_2 = 0.0;
double a_2_0 = 0.0, a_2_1 = 0.0;
double a_3_0 = 0.0;
double b_0_0 = 0.0, b_0_1 = 0.0, b_0_2 = 0.0, b_0_3 = 0.0;
double b_1_0 = 0.0, b_1_1 = 0.0, b_1_2 = 0.0;
double b_2_0 = 0.0, b_2_1 = 0.0;
double b_3_0 = 0.0;

double ap_0_0 = 0.0, ap_0_1 = 0.0, ap_0_2 = 0.0, ap_0_3 = 0.0;
double ap_1_0 = 0.0, ap_1_1 = 0.0, ap_1_2 = 0.0;
double ap_2_0 = 0.0, ap_2_1 = 0.0;
double ap_3_0 = 0.0;
double bp_0_0 = 0.0, bp_0_1 = 0.0, bp_0_2 = 0.0, bp_0_3 = 0.0;
double bp_1_0 = 0.0, bp_1_1 = 0.0, bp_1_2 = 0.0;
double bp_2_0 = 0.0, bp_2_1 = 0.0;
double bp_3_0 = 0.0;

/// MARK: Solver cached kinematics

double SIN_dec0     = 0.0;
double COS_dec0     = 1.0;
double SIN_dec_ref  = 0.0;
double COS_dec_ref  = 1.0;
double x_new_float  = 0.0;
double y_new_float  = 0.0;
double referenceX   = 0.0;
double referenceY   = 0.0;
double pedestal_s   = 0.0;
double mag2         = 0.0;
std::string solution_str;

/// MARK: Stacking accumulators

double sum_exp         = 0.0;
double sum_temp        = 0.0;
double airmass_sum     = 0.0;
double airmass         = 0.0;
double jd_sum          = 0.0;
double jd_start_first  = 0.0;
double jd_end_last     = 0.0;
double jd_start        = 0.0;
double jd_mid          = 0.0;
double jd_end          = 0.0;
double jd_mid_reference = 0.0;
int    images_selected = 0;

int counterR   = 0, counterRdark   = 0, counterRflat   = 0, counterRbias   = 0;
int counterG   = 0, counterGdark   = 0, counterGflat   = 0, counterGbias   = 0;
int counterB   = 0, counterBdark   = 0, counterBflat   = 0, counterBbias   = 0;
int counterL   = 0, counterLdark   = 0, counterLflat   = 0, counterLbias   = 0;
int counterRGB = 0, counterRGBdark = 0, counterRGBflat = 0, counterRGBbias = 0;
int exposureR = 0, exposureG = 0, exposureB = 0, exposureL = 0, exposureRGB = 0;
double temperatureR   = 0.0, temperatureG = 0.0, temperatureB = 0.0;
double temperatureL   = 0.0, temperatureRGB = 0.0;

/// MARK: Solver settings

int    max_stars_setting     = 500;
int    downsample_setting    = 0;      // 0 = auto
double quad_tolerance        = 0.007;
double min_star_size_arcsec  = 1.5;
double search_radius_deg     = 180.0;
double search_fov_deg        = 0.0;    // 0 = auto from header
bool   force_oversize        = false;
bool   add_sip               = false;
bool   check_pattern_filter  = false;

/// MARK: Stacking settings

bool   use_manual_align         = false;
bool   use_ephemeris_alignment  = false;
bool   use_astrometry_internal  = false;
bool   skip_alignment           = false;
double hfd_min_setting          = 0.0;    // 0 = auto (engine defaults to 0.8)
double hfd_max_setting          = 10.0;   // Matches Pascal. Real stars fit well inside; extended sources (galaxies) are correctly excluded.
double sigma_clip_factor        = 2.0;

/// MARK: Memos

std::vector<std::string> memo1_lines;
std::vector<std::string> memo2_lines;


/// MARK: Filesystem / executable environment

std::filesystem::path application_path;
std::vector<std::string> recent_files;
int raw_conversion_program_index = 0;
std::vector<std::string> head1 = {
	/*  0 */ "SIMPLE  =                    T / FITS header                                    ",
	/*  1 */ "BITPIX  =                    8 / Bits per entry                                 ",
	/*  2 */ "NAXIS   =                    2 / Number of dimensions                           ",
	/*  3 */ "NAXIS1  =                  100 / length of x axis                               ",
	/*  4 */ "NAXIS2  =                  100 / length of y axis                               ",
	/*  5 */ "NAXIS3  =                    3 / length of z axis (mostly colors)               ",
	/*  6 */ "EQUINOX =               2000.0 / Equinox of coordinates                         ",
	/*  7 */ "DATAMIN =                    0 / Minimum data value                             ",
	/*  8 */ "DATAMAX =                  255 / Maximum data value                             ",
	/*  9 */ "BZERO   =                  0.0 / physical_value = BZERO + BSCALE * array_value  ",
	/* 10 */ "BSCALE  =                  1.0 / physical_value = BZERO + BSCALE * array_value  ",
	/* 11 */ "CTYPE1  = 'RA---TAN'           / first parameter RA  ,  projection TANgential   ",
	/* 12 */ "CTYPE2  = 'DEC--TAN'           / second parameter DEC,  projection TANgential   ",
	/* 13 */ "CUNIT1  = 'deg     '           / Unit of coordinates                            ",
	/* 14 */ "CRPIX1  =                  0.0 / X of reference pixel                           ",
	/* 15 */ "CRPIX2  =                  0.0 / Y of reference pixel                           ",
	/* 16 */ "CRVAL1  =                  0.0 / RA of reference pixel (deg)                    ",
	/* 17 */ "CRVAL2  =                  0.0 / DEC of reference pixel (deg)                   ",
	/* 18 */ "CDELT1  =                  0.0 / X pixel size (deg)                             ",
	/* 19 */ "CDELT2  =                  0.0 / Y pixel size (deg)                             ",
	/* 20 */ "CROTA1  =                  0.0 / Image twist X axis(deg)                        ",
	/* 21 */ "CROTA2  =                  0.0 / Image twist Y axis deg) E of N if not flipped  ",
	/* 22 */ "CD1_1   =                  0.0 / CD matrix to convert (x,y) to (Ra, Dec)        ",
	/* 23 */ "CD1_2   =                  0.0 / CD matrix to convert (x,y) to (Ra, Dec)        ",
	/* 24 */ "CD2_1   =                  0.0 / CD matrix to convert (x,y) to (Ra, Dec)        ",
	/* 25 */ "CD2_2   =                  0.0 / CD matrix to convert (x,y) to (Ra, Dec)        ",
	/* 26 */ "PLTSOLVD=                    T / ASTAP from hnsky.org                           ",
	/* 27 */ "END                                                                             ",
	/* 28 */ "                                                                                ",
};
	
} // namespace
