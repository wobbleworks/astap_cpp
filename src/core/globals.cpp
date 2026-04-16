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

/// MARK: Memos

std::vector<std::string> memo1_lines;
std::vector<std::string> memo2_lines;

/// MARK: Filesystem / executable environment

std::filesystem::path application_path;
std::vector<std::string> recent_files;
int raw_conversion_program_index = 0;
std::vector<std::string> head1;
	
} // namespace
