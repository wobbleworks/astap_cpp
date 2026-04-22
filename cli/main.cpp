// ASTAP++ command-line front-end.
//
// Ported from the CLI path inside Tmainwindow.FormShow in astap_main.pas
// (lines ~13440-13810). Maintains argument-name compatibility with the
// original so existing ASTAP callers (APT, CCDCiel, NINA, SGP, Ekos) work
// unchanged.
//
// Dispatch:
//   1. PlateSolve2-style positional comma-separated argv   -> platesolve2
//   2. -h / -help                                          -> print_help
//   3. -focus1 ... -focus2 ... (+ optional -focus3, ...)   -> focus fit
//   4. -analyse / -extract / -extract2                     -> analysis-only
//   5. -f <file>                                           -> solve
//   6. -stack <path>                                       -> live-stack
//   7. positional first argv                               -> image view
//      (CLI viewer not implemented here; requires GUI)
//
// Exit codes (from Pascal source lines 13739-13748):
//    0 no errors
//    1 no solution
//    2 not enough stars detected
//   16 error reading image file
//   32 no star database found
//   33 error reading star database
//   34 error updating input file
//
// MPL-2.0. See readme.txt in the original Pascal source.

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "../src/types.h"
#include "../src/core/globals.h"         // all cross-module runtime state

// -----------------------------------------------------------------------------
// Backing-module headers. Many of these are faithful-but-stubbed ports; the
// CLI scaffolds the real dispatch and will start working as the remaining
// TODOs in those modules get filled in.
// -----------------------------------------------------------------------------
#include "../src/core/fits.h"            // load_fits, save_fits, savefits_update_header, header-edit
#include "../src/core/image_io.h"        // load_image, add_recent_file
#include "../src/core/wcs.h"             // write_astronomy_wcs
#include "../src/core/platform.h"        // execute_and_wait, wait_ms, log_to_file, IniFile
#include "../src/core/util.h"            // fits_file_name, tiff_file_name, strtofloat2
#include "../src/core/imaging.h"         // bin_X2X3X4
#include "../src/reference/star_database.h"       // database_path, name_database
#include "../src/solving/astrometric_solving.h"   // solve_image
#include "../src/stacking/stack.h"       // analyse_image, Background
#include "../src/stacking/live_stacking.h"        // stack_live
#include "../src/image/tiff.h"           // save_tiff_16
#include "stb_image_decoder.h"           // install_stb_image_decoder()

// Forward decls for functionality not yet consolidated into a ported header.
// When their owning module lands, drop these and include it instead.
// Include the modules that back the CLI's extra features.
#include "../src/core/sqm.h"            // calculate_sqm, bortle, altitudefloat, sqmfloat, airmass, centalt
#include "../src/core/hjd.h"            // airmass_calc
#include "../src/core/hyperbola.h"      // find_best_hyperbola_fit

namespace astap::core {
// Pascal: save_annotated_jpg (astap_main.pas:13040) — renders annotations
// on top of the loaded image and saves a JPEG. GUI-coupled in the original;
// a headless port needs a JPEG encoder (stb_image_write or libjpeg).
// Declared here as a weak stub so -annotate dispatch still compiles; the
// stub returns false without producing a file.
inline bool save_annotated_jpg(const std::filesystem::path& /*output_base*/) {
    std::cerr << "astap: -annotate requires a JPEG encoder adapter "
                 "(not yet installed); skipping.\n";
    return false;
}

// Pascal: platesolve2_command (astap_main.pas:13102) — detects and processes
// the positional PlateSolve2 argv convention. Declared here with the
// (argc, argv) shape the CLI expects; the backing implementation lives in
// platform.cpp and takes no arguments (reads sys argv directly).
inline bool platesolve2_command(int /*argc*/, char* /*argv*/[]) {
    return platesolve2_command();
}

// Pascal: do_stretching — builds the gamma LUT used by the display stretch.
// GUI-only in the original; no-op for CLI.
inline void do_stretching() noexcept {}

// Pascal: use_histogram (astap_main.pas:7946). Reachable via imaging.h; the
// CLI only needs the side-effect of populating the histogram arrays for
// subsequent plot_annotations. Here, no-op — the plot path is GUI only.
inline void use_histogram(const astap::ImageArray& /*img*/, bool /*update_hist*/) {}
}  // namespace astap::core

namespace astap::solving {
// CLI-level wrapper: loads each image, measures HFD, fits the hyperbola.
// Backed by astap::core::find_best_hyperbola_fit (pure math, in core/hyperbola.h).
struct FocusFitResult {
    double focus_best;
    double lowest_error;
    bool   ok;
};
FocusFitResult fit_focus_hyperbola(std::span<const std::filesystem::path> images);
}  // namespace astap::solving

// All cross-module globals (esc_pressed, commandline_execution, filename2,
// head, img_loaded, memo1_lines, memo2_lines, errorlevel, and the solver
// settings `max_stars_setting`, `quad_tolerance`, `search_fov_deg`, etc.)
// live in src/core/globals.{h,cpp}. `database_path` and `name_database`
// come from src/reference/star_database.h (astap::reference namespace).

// -----------------------------------------------------------------------------
// Argument parsing. Mirrors Lazarus TApplication.hasOption / GetOptionValue
// enough that Pascal-style "-key value" and "--key=value" both work.
// -----------------------------------------------------------------------------
namespace {

struct Args {
    // Flags / options (only set when explicitly passed).
    bool                       help = false;
    bool                       debug = false;
    bool                       file_specified = false;
    bool                       focus_request = false;
    bool                       analyse_specified = false;
    bool                       extract_specified = false;
    bool                       extract2_specified = false;
    bool                       wcs_output = false;
    bool                       update_input = false;
    bool                       annotate = false;
    bool                       sip = false;
    bool                       check = false;
    bool                       log = false;
    bool                       stack_mode = false;
    std::optional<std::string> file;
    std::optional<double>      radius_deg;
    std::optional<double>      fov_deg;
    std::optional<double>      ra_hours;
    std::optional<double>      spd_deg;       // south pole distance
    std::optional<int>         max_stars;
    std::optional<double>      tolerance;
    std::optional<double>      min_star_size;
    std::optional<int>         downsample;    // 0..4
    std::optional<std::string> database_path; // -d
    std::optional<std::string> database_name; // -D
    std::optional<std::string> output;        // -o base
    std::optional<std::string> speed;         // "auto" | "slow"
    std::optional<double>      analyse_snr;
    std::optional<double>      extract_snr;
    std::optional<double>      extract2_snr;
    std::optional<int>         tofits_binning;
    std::optional<int>         sqm_pedestal;
    std::optional<std::string> stack_path;
    std::vector<std::filesystem::path> focus_files;
    std::string                raw_cmdline;
    // Positional argv remainder (Pascal treats non-option argv[1] as filename).
    std::vector<std::string>   positional;
};

struct ArgParseError { std::string what; };

// Check whether `raw` (already stripped of leading '-' or '--') is an exact
// match for `name`, or starts with `name=`. Returns the value after '=' if
// any, or an empty optional otherwise.
std::optional<std::string_view> match_key(std::string_view raw,
                                          std::string_view name) {
    if (raw == name) return std::string_view{};
    if (raw.starts_with(name) && raw.size() > name.size()
        && raw[name.size()] == '=') {
        return raw.substr(name.size() + 1);
    }
    return std::nullopt;
}

// Fetch the value for `key`. Lazarus accepts both "-key value" (space
// separated) and "--key=value". `advance_to` is updated to the next unused
// argv index.
std::optional<std::string> take_value(std::span<char* const> argv,
                                      std::size_t            i,
                                      std::string_view       key,
                                      std::size_t&           advance_to) {
    std::string_view arg = argv[i];
    std::string_view raw = arg;
    if (raw.starts_with("--")) raw.remove_prefix(2);
    else if (raw.starts_with('-')) raw.remove_prefix(1);
    else return std::nullopt;

    auto matched = match_key(raw, key);
    if (!matched) return std::nullopt;

    if (!matched->empty()) { advance_to = i + 1; return std::string{*matched}; }
    if (i + 1 >= argv.size()) {
        throw ArgParseError{std::format("option -{} requires a value", key)};
    }
    advance_to = i + 2;
    return std::string{argv[i + 1]};
}

bool take_flag(std::string_view arg, std::string_view key) {
    std::string_view raw = arg;
    if (raw.starts_with("--")) raw.remove_prefix(2);
    else if (raw.starts_with('-')) raw.remove_prefix(1);
    else return false;
    return raw == key;
}

double parse_double(std::string_view s) {
    return astap::core::strtofloat2(std::string{s});
}

int parse_int(std::string_view s, int dflt = 0) {
    return astap::core::strtoint2(std::string{s}, dflt);
}

Args parse_args(int argc, char* argv[]) {
    Args out;

    // Record the raw command-line as a single string (for -log + FITS
    // history COMMENT), matching Pascal's `cmdline` global.
    for (int i = 0; i < argc; ++i) {
        if (i) out.raw_cmdline.push_back(' ');
        out.raw_cmdline.append(argv[i]);
    }

    std::span<char* const> av(argv + 1, argv + argc);
    std::size_t next = 0;
    for (std::size_t i = 0; i < av.size(); i = next) {
        next = i + 1;
        std::string_view arg = av[i];

        auto try_value = [&](std::string_view key, auto setter) -> bool {
            if (auto v = take_value(av, i, key, next); v) { setter(*v); return true; }
            return false;
        };

        if (take_flag(arg, "h") || take_flag(arg, "help")) { out.help = true; continue; }
        if (take_flag(arg, "debug"))    { out.debug = true; continue; }
        if (take_flag(arg, "check"))    { out.check = true; continue; }
        if (take_flag(arg, "sip"))      { out.sip = true; continue; }
        if (take_flag(arg, "wcs"))      { out.wcs_output = true; continue; }
        if (take_flag(arg, "log"))      { out.log = true; continue; }
        if (take_flag(arg, "update"))   { out.update_input = true; continue; }
        if (take_flag(arg, "annotate")) { out.annotate = true; continue; }

        if (try_value("f",      [&](auto v){ out.file = v; out.file_specified = true; })) continue;
        if (try_value("r",      [&](auto v){ out.radius_deg   = parse_double(v); }))     continue;
        if (try_value("fov",    [&](auto v){ out.fov_deg      = parse_double(v); }))     continue;
        if (try_value("ra",     [&](auto v){ out.ra_hours     = parse_double(v); }))     continue;
        if (try_value("spd",    [&](auto v){ out.spd_deg      = parse_double(v); }))     continue;
        if (try_value("s",      [&](auto v){ out.max_stars    = parse_int(v, 500); }))   continue;
        if (try_value("t",      [&](auto v){ out.tolerance    = parse_double(v); }))     continue;
        if (try_value("m",      [&](auto v){ out.min_star_size= parse_double(v); }))     continue;
        if (try_value("z",      [&](auto v){ out.downsample   = parse_int(v, 0); }))     continue;
        if (try_value("d",      [&](auto v){ out.database_path= v; }))                   continue;
        if (try_value("D",      [&](auto v){ out.database_name= v; }))                   continue;
        if (try_value("o",      [&](auto v){ out.output       = v; }))                   continue;
        if (try_value("speed",  [&](auto v){ out.speed        = v; }))                   continue;
        if (try_value("analyse",[&](auto v){ out.analyse_snr  = parse_double(v); out.analyse_specified  = true; })) continue;
        if (try_value("extract",[&](auto v){ out.extract_snr  = parse_double(v); out.extract_specified  = true; })) continue;
        if (try_value("extract2",[&](auto v){out.extract2_snr = parse_double(v); out.extract2_specified = true; })) continue;
        if (try_value("tofits", [&](auto v){ out.tofits_binning = parse_int(v, 1); }))   continue;
        if (try_value("sqm",    [&](auto v){ out.sqm_pedestal = parse_int(v, 0); }))     continue;
        if (try_value("stack",  [&](auto v){ out.stack_path   = v; out.stack_mode = true; })) continue;

        // -focusN  (N = 1, 2, 3, ...)
        if (arg.starts_with("-focus") || arg.starts_with("--focus")) {
            std::string_view raw = arg;
            if (raw.starts_with("--")) raw.remove_prefix(2);
            else raw.remove_prefix(1);
            std::string_view numstr = raw.substr(5);
            // Digits-only check; any non-digit falls through to "unknown".
            bool digits_only = !numstr.empty()
                && std::all_of(numstr.begin(), numstr.end(),
                               [](char c){ return c >= '0' && c <= '9'; });
            if (digits_only) {
                if (i + 1 >= av.size()) {
                    throw ArgParseError{std::format("option {} requires a filename", arg)};
                }
                out.focus_files.emplace_back(av[i + 1]);
                out.focus_request = true;
                next = i + 2;
                continue;
            }
        }

        // Unknown -flag: Pascal's TApplication silently accepts; we pass
        // through as positional so the caller sees where we stopped.
        if (arg.starts_with('-')) {
            throw ArgParseError{std::format("unrecognised option: {}", arg)};
        }
        out.positional.emplace_back(arg);
    }

    return out;
}

void print_help(std::ostream& os) {
    os <<
        "Solver command-line usage:\n"
        "-f  filename\n"
        "-r  radius_area_to_search[degrees]\n"
        "-fov height_field[degrees]\n"
        "-ra  center_right_ascension[hours]\n"
        "-spd center_south_pole_distance[degrees]\n"
        "-s  max_number_of_stars {typical 500}\n"
        "-t  tolerance\n"
        "-m  minimum_star_size[\"]\n"
        "-z  downsample_factor[0,1,2,3,4] {Downsample prior to solving. 0 is auto}\n"
        "\n"
        "-check  {Apply check pattern filter prior to solving. Use for raw OSC images only when binning is 1x1}\n"
        "-d  path {Specify a path to the star database}\n"
        "-D  abbreviation {Specify a star database [d80,d50,..]}\n"
        "-o  file {Name the output files with this base path & file name}\n"
        "-sip     {Add SIP (Simple Image Polynomial) coefficients}\n"
        "-speed mode[auto/slow] {Slow forces a larger database area for more overlap}\n"
        "-wcs  {Write a .wcs file in similar format as Astrometry.net. Else text style.}\n"
        "-log   {Write the solver log to file}\n"
        "-update  {update the FITS/TIFF header with the found solution. Jpeg/png will be written as fits}\n"
        "\n"
        "Analyse options:\n"
        "-analyse snr_min {Analyse only and report median HFD and number of stars used}\n"
        "-extract snr_min {As -analyse but additionally export info of all detectable stars to a .csv file}\n"
        "-extract2 snr_min {Solve and export info of all detectable stars to a .csv file including ra, dec.}\n"
        "\n"
        "Extra options:\n"
        "-annotate  {Produce deepsky annotated jpg file}\n"
        "-debug  {Show GUI and stop prior to solving}  (GUI disabled in this build)\n"
        "-tofits  binning[1,2,3,4]  {Make new fits file from PNG/JPG file input}\n"
        "-sqm pedestal  {add measured sqm, centalt, airmass values to the solution}\n"
        "-focus1 file1.fit -focus2 file2.fit ....  {Find best focus using files and hyperbola curve fitting.\n"
        "                                           Errorlevel is focuspos*1E4 + rem.error*1E3}\n"
        "-stack  path {startup with live stack tab and path selected}\n"
        "\n"
        "Preference will be given to the command-line values.\n"
        "CSV files are written with a dot as decimal separator.\n"
        "Solver result will be written to filename.ini and filename.wcs.\n";
}

// Apply CLI-supplied overrides to the module-level settings consumed by the
// solver. Matches the Pascal `hasoption(...)` block at lines 13536-13567.
void apply_solver_overrides(const Args& a) {
    if (a.fov_deg)       { astap::search_fov_deg = *a.fov_deg; astap::fov_specified = true; }
    if (a.radius_deg)    { astap::search_radius_deg = *a.radius_deg; }
    if (a.ra_hours)      { /* TODO: push into Header/head.ra0 once wired */ }
    if (a.spd_deg)       { astap::head.dec0 = *a.spd_deg - 90.0; }
    if (a.downsample)    { astap::downsample_setting = *a.downsample; }
    if (a.max_stars)     { astap::max_stars_setting = *a.max_stars; }
    if (a.tolerance)     { astap::quad_tolerance = *a.tolerance; }
    if (a.min_star_size) { astap::min_star_size_arcsec = *a.min_star_size; }
    if (a.sip)           { astap::add_sip = true; }
    if (a.check)         { astap::check_pattern_filter = true; }
    if (a.speed)         { astap::force_oversize = (*a.speed == "slow"); }
    if (a.database_path) {
        std::filesystem::path p(*a.database_path);
        if (!p.empty() && p.native().back() != std::filesystem::path::preferred_separator) {
            p /= "";   // ensure trailing separator (matches Pascal DirectorySeparator)
        }
        astap::reference::database_path = p;
    }
    if (a.database_name) { astap::reference::name_database = *a.database_name; }
}

// -----------------------------------------------------------------------------
// Dispatch branches
// -----------------------------------------------------------------------------

int run_analyse(const Args& a,
                const std::filesystem::path& output_base,
                bool extract_csv) {
    using namespace astap;
    double snr_min = a.analyse_snr.value_or(a.extract_snr.value_or(0.0));
    if (snr_min == 0.0) snr_min = 30.0;

    int report_type = extract_csv ? 2 /* csv only */ : 0 /* report only */;
    int    star_counter = 0;
    double hfd_median   = 0.0;
    astap::Background bck{};

    astap::stacking::analyse_image(astap::img_loaded,
                                   astap::head,
                                   snr_min,
                                   report_type,
                                   star_counter,
                                   bck,
                                   hfd_median);

    std::cout << std::format("HFD_MEDIAN={:.1f}\n", hfd_median);
    std::cout << std::format("STARS={}\n", star_counter);

#if defined(_WIN32)
    // Packed errorlevel matches the Pascal Windows branch: hfd_median*100*1e6 + stars.
    long long packed = static_cast<long long>(std::round(hfd_median * 100.0)) * 1'000'000LL
                     + star_counter;
    if (packed > std::numeric_limits<int>::max()) packed = std::numeric_limits<int>::max();
    return static_cast<int>(packed);
#else
    return astap::errorlevel;
#endif
}

int run_focus(const Args& a) {
    auto r = astap::solving::fit_focus_hyperbola(a.focus_files);
    if (!r.ok) return 1;

    std::cout << std::format("FOCUS={:.1f}\n", r.focus_best);
    std::cout << std::format("ERROR_MIN={:.5f}\n", r.lowest_error);

#if defined(_WIN32)
    int f = static_cast<int>(std::round(r.focus_best));
    int e = std::min(9999, static_cast<int>(std::round(r.lowest_error * 1000.0)));
    return f * 10000 + e;
#else
    return astap::errorlevel;
#endif
}

int run_stack(const Args& a) {
    // TODO: Pascal opens the live-stacking tab and calls stack_images1Click(nil).
    // For a headless CLI this should drive astap::stacking::stack_live(path)
    // directly. That function is a blocking polling loop that terminates when
    // esc_pressed becomes true. Hooking Ctrl-C to a sigaction that sets
    // esc_pressed = true is the next step; for now we just call through.
    astap::stacking::stack_live(std::filesystem::path{*a.stack_path});
    return astap::errorlevel;
}

int run_solve(const Args& a, const std::filesystem::path& output_base) {
    using namespace astap;

    bool solved = astap::solving::solve_image(astap::img_loaded,
                                              astap::head,
                                              astap::memo1_lines,
                                              /*get_hist=*/true,
                                              astap::check_pattern_filter);
    if (!solved) {
        (void)astap::core::write_ini(output_base, false);
        if (astap::errorlevel == 0) astap::errorlevel = 1;
        return astap::errorlevel;
    }

    // -sqm: compute sky background magnitude.
    if (a.sqm_pedestal) {
        int ped = *a.sqm_pedestal;
        if (astap::core::calculate_sqm(astap::head, /*get_backgr=*/false,
                                       /*get_histogr=*/false, ped)) {
            if (astap::core::centalt.empty()) {
                astap::core::centalt = std::format("{:.2f}", astap::core::altitudefloat);
                astap::core::update_text(astap::memo1_lines, "CENTALT =",
                    "'" + astap::core::centalt + "'              / [deg] Nominal altitude of center of image    ");
                astap::core::update_text(astap::memo1_lines, "OBJCTALT=",
                    "'" + astap::core::centalt + "'              / [deg] Nominal altitude of center of image    ");
            }
            astap::core::update_text(astap::memo1_lines, "SQM     = ",
                std::format("{:.2f}               / Sky background [magn/arcsec^2]",
                            astap::core::sqmfloat));
            astap::core::update_text(astap::memo1_lines, "COMMENT SQM",
                std::format(", used {} as pedestal value", ped));
        } else {
            astap::core::update_text(astap::memo1_lines, "SQM     =",
                "'Error calculating SQM value! Check in the SQM menu (ctrl+Q) first.'");
        }

        if (astap::airmass == 0.0) {
            astap::airmass = astap::core::airmass_calc(astap::core::altitudefloat);
            astap::core::update_generic(astap::memo1_lines, "AIRMASS ",
                std::format("{:.4f}", astap::airmass),
                "Relative optical path.                        ");
        }
    }

    (void)astap::core::write_ini(output_base, true);
    astap::core::add_long_comment(astap::memo1_lines, "cmdline:" + a.raw_cmdline);

    // -update: write the found solution back to the input file.
    if (a.update_input) {
        bool wresult = false;
        std::filesystem::path in(*a.file);
        if (astap::core::fits_file_name(in.extension().string())) {
            wresult = astap::core::savefits_update_header(astap::memo1_lines, in);
        } else if (astap::core::tiff_file_name(in.extension().string())) {
            wresult = astap::image::save_tiff_16(astap::img_loaded, in,
                                                 /*description=*/"", false, false);
        } else {
            auto fits_out = in;
            fits_out.replace_extension(".fits");
            wresult = astap::core::save_fits(astap::img_loaded,
                                             astap::memo1_lines, fits_out,
                                             16, /*override=*/true);
        }
        if (!wresult) {
            std::cerr << "Error updating input file\n";
            astap::errorlevel = 34;
        }
    }

    // Trim NAXIS1/2, force NAXIS=0 for the WCS/header text file.
    astap::core::remove_key(astap::memo1_lines, "NAXIS1  =", true);
    astap::core::remove_key(astap::memo1_lines, "NAXIS2  =", true);
    astap::core::update_integer(astap::memo1_lines, "NAXIS   =",
        " / Minimal header                                 ", 0);

    auto wcs_path = output_base;
    wcs_path.replace_extension(".wcs");
    if (a.wcs_output) {
        astap::core::write_astronomy_wcs(wcs_path, astap::memo1_lines);
    } else {
        try {
            std::ofstream ofs(wcs_path);
            for (auto const& line : astap::memo1_lines) ofs << line << '\n';
        } catch (...) {
            // Matches Pascal's swallowed exception (sometimes APT locks file).
        }
    }

    // -annotate: requires a GUI-rendered histogram + stretch + JPG; stubbed.
    if (a.annotate) {
        astap::core::use_histogram(astap::img_loaded, /*update_hist=*/false);
        astap::core::save_annotated_jpg(output_base);
    }

    // -tofits: convert non-FITS input with optional binning.
    if (a.tofits_binning) {
        std::filesystem::path in(*a.file);
        if (!astap::core::fits_file_name(in.extension().string())) {
            int bin = *a.tofits_binning;
            if (bin > 1) astap::core::bin_X2X3X4(astap::img_loaded, astap::head,
                                                  astap::memo1_lines, bin,
                                                  astap::filename2);
            astap::core::use_histogram(astap::img_loaded, false);
            auto fit_out = output_base;
            fit_out.replace_extension(".fit");
            (void)astap::core::save_fits(astap::img_loaded, astap::memo1_lines, fit_out, 8, true);
        }
    }

    // -extract2: export full star catalog with celestial coords after solving.
    if (a.extract2_specified) {
        double snr_min = a.extract2_snr.value_or(30.0);
        if (snr_min == 0.0) snr_min = 30.0;
        int    dummy_count = 0;
        double dummy_hfd   = 0.0;
        astap::Background bck{};
        astap::stacking::analyse_image(astap::img_loaded, astap::head, snr_min,
                                       /*report_type=*/2, dummy_count, bck, dummy_hfd);
    }

    if (astap::commandline_log) {
        auto logp = output_base;
        logp.replace_extension(".log");
        std::ofstream logf(logp);
        for (auto const& line : astap::memo2_lines) logf << line << '\n';
    }

    return astap::errorlevel;
}

}  // namespace

// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------
// Load each image, measure median HFD at the centre with
// astap::core::HFD, then fit the focus-hyperbola via the ported pure-math
// primitive. This is the CLI-level glue; the pure math lives in
// astap::core::find_best_hyperbola_fit.
namespace astap::solving {
FocusFitResult fit_focus_hyperbola(std::span<const std::filesystem::path> images) {
    FocusFitResult r{0.0, 0.0, false};
    // TODO: extract focus_position from each frame (FOCUSPOS/FITS header) and
    //       measure_hfd via astap::core::HFD. Pending core::HFD wire-up the
    //       fit is deferred. The current stub reports ok=false so the CLI
    //       exits with a non-zero errorlevel instead of emitting bogus data.
    (void)images;
    return r;
}
}  // namespace astap::solving

int main(int argc, char* argv[]) {
    // Install the stb_image-backed raster decoder as the default. Hosts that
    // want richer format coverage (TIFF, RAW, ...) can replace this by
    // calling astap::core::set_image_decoder() themselves before dispatch.
    astap::cli::install_stb_image_decoder();

    // Wire Ctrl-C / SIGTERM to astap::esc_pressed so long-running loops
    // (stack_live, solver sweeps) can quit cleanly.
    (void)astap::core::install_cancel_signal_handler();

    // 1) PlateSolve2-style positional argv. platesolve2_command() inspects
    //    argv itself and returns true if it consumed the command (writing
    //    its own .apm output and setting errorlevel before returning).
    if (astap::core::platesolve2_command(argc, argv)) {
        return astap::errorlevel;
    }

    // 2) Parse the ASTAP-native option set.
    Args a;
    try {
        a = parse_args(argc, argv);
    } catch (const ArgParseError& e) {
        std::cerr << "astap: " << e.what << "\nTry -h for help.\n";
        return 2;
    }
    astap::cmdline = a.raw_cmdline;

    if (a.help) { print_help(std::cout); return 0; }

    // 3) -stack live-stacking mode short-circuits everything else.
    if (a.stack_mode) { return run_stack(a); }

    // Anything below requires -f, -debug, -focus1+, or a positional filename.
    if (!a.file_specified && !a.focus_request && a.positional.empty()) {
        // Pascal behaviour: no useful args -> just builds the gamma curve and
        // exits (for GUI reuse). In the CLI this means: print help briefly.
        print_help(std::cout);
        return 0;
    }

    astap::commandline_execution = true;
    astap::commandline_log = (a.debug || a.log);
    if (astap::commandline_log) astap::memo2_lines.push_back(astap::cmdline);

    // 4) Load the input image (if specified).
    if (a.file_specified) {
        astap::filename2 = *a.file;
        bool loaded = astap::core::load_image(astap::filename2,
                                              astap::img_loaded,
                                              astap::head,
                                              astap::memo1_lines,
                                              /*re_center=*/false,
                                              /*plot=*/false);
        if (!loaded) {
            astap::errorlevel = 16;
            std::cerr << "Error reading image file: " << astap::filename2 << '\n';
            return astap::errorlevel;
        }
    }

    apply_solver_overrides(a);

    // 5) Focus-fit mode.
    if (a.focus_request) { return run_focus(a); }

    // 6) Determine the output base path (-o overrides the input filename).
    std::filesystem::path output_base =
        a.output.has_value() ? std::filesystem::path{*a.output}
                             : std::filesystem::path{astap::filename2};

    // 7) Analysis-only branches (-analyse, -extract). No solve, no .wcs.
    if (a.file_specified && (a.analyse_specified || a.extract_specified)) {
        return run_analyse(a, output_base, a.extract_specified);
    }

    // 8) Full solve. -debug would normally reveal the GUI here; this CLI
    //    build has no GUI, so debug falls through to a plain solve.
    if (a.debug) {
        std::cerr << "-debug: GUI not available in this build; running a "
                     "plain solve instead.\n";
    }

    return run_solve(a, output_base);
}
