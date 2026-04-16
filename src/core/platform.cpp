///----------------------------------------
///      @file platform.cpp
///   @ingroup ASTAP++
///     @brief Platform integration and settings persistence for ASTAP++.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "platform.h"

#include "globals.h"   // esc_pressed, cmdline, warning_str, errorlevel, head

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <numbers>
#include <sstream>
#include <system_error>
#include <thread>

#if defined(_WIN32)
  #include <windows.h>
#else
  #include <spawn.h>
  #include <sys/wait.h>
  #include <unistd.h>
  extern char** environ;
#endif

///----------------------------------------
namespace astap::core {
///----------------------------------------

namespace fs = std::filesystem;

/// MARK: - Process Launch Helpers

namespace {

// Tokenise a command into argv-style arguments, respecting double-quoted spans.
// Good enough for ASTAP's use -- it only shells out to known executables with
// paths the host assembles. Not a full POSIX shell parser.
[[nodiscard]] std::vector<std::string> tokenise_command(std::string_view cmd) {
    auto out = std::vector<std::string>{};
    auto current = std::string{};
    auto in_quote = false;
    for (auto c : cmd) {
        if (c == '"') {
            in_quote = !in_quote;
            continue;
        }
        if (!in_quote && (c == ' ' || c == '\t')) {
            if (!current.empty()) {
                out.emplace_back(std::move(current));
                current.clear();
            }
            continue;
        }
        current.push_back(c);
    }
    if (!current.empty()) {
        out.emplace_back(std::move(current));
    }
    return out;
}
 
}  // anonymous namespace

bool execute_and_wait(const std::string& cmd, bool show_console) {
#if defined(_WIN32)
    // Windows: CreateProcess with STARTUPINFO; show_console controls visibility.
    STARTUPINFOA si{};
    si.cb = sizeof(si);
    si.dwFlags = STARTF_USESHOWWINDOW;
    si.wShowWindow = show_console ? SW_SHOWMINNOACTIVE : SW_HIDE;
    
    PROCESS_INFORMATION pi{};
    // CreateProcess mutates its command-line argument, so copy to a buffer.
    auto buf = cmd;
    auto ok = ::CreateProcessA(nullptr, buf.data(), nullptr, nullptr, FALSE,
                               0, nullptr, nullptr, &si, &pi);
    if (!ok) {
        return false;
    }
    
    ::WaitForSingleObject(pi.hProcess, INFINITE);
    auto code = DWORD{1};
    ::GetExitCodeProcess(pi.hProcess, &code);
    ::CloseHandle(pi.hProcess);
    ::CloseHandle(pi.hThread);
    return code == 0;
#else
    // POSIX: posix_spawn + waitpid. show_console is ignored --
    // UNIX processes don't have a hidden-window concept.
    [[maybe_unused]] auto console_flag = show_console;
    
    auto tokens = tokenise_command(cmd);
    if (tokens.empty()) {
        return false;
    }
    
    auto argv = std::vector<char*>{};
    argv.reserve(tokens.size() + 1);
    for (auto& t : tokens) {
        argv.push_back(t.data());
    }
    argv.push_back(nullptr);
    
    auto pid = pid_t{0};
    auto spawn_rc = ::posix_spawnp(&pid, argv[0], nullptr, nullptr,
                                   argv.data(), environ);
    if (spawn_rc != 0) {
        errno = spawn_rc;
        return false;
    }
    
    auto status = 0;
    while (true) {
        auto w = ::waitpid(pid, &status, 0);
        if (w == -1) {
            if (errno == EINTR) {
                continue;  // e.g. Ctrl-C in parent; retry.
            }
            return false;
        }
        break;
    }
    if (WIFEXITED(status)) {
        return WEXITSTATUS(status) == 0;
    }
    if (WIFSIGNALED(status)) {
        return false;
    }
    return false;
#endif
}

bool execute_unix(const std::string& executable,
                  const std::vector<std::string>& params,
                  bool show_output,
                  std::string& captured_output) {
    captured_output.clear();
    
    // Build the shell command
    auto cmd = std::ostringstream{};
    cmd << '"' << executable << '"';
    for (const auto& p : params) {
        cmd << ' ' << '"' << p << '"';
    }
    if (show_output) {
        cmd << " 2>&1";
    }
    
    // Execute via popen
#if defined(_WIN32)
    auto* pipe = ::_popen(cmd.str().c_str(), "r");
#else
    auto* pipe = ::popen(cmd.str().c_str(), "r");
#endif
    if (!pipe) {
        return false;
    }
    
    // Read output
    char buffer[4096];
    while (std::fgets(buffer, sizeof(buffer), pipe)) {
        if (show_output) {
            captured_output.append(buffer);
        }
    }
    
    // Close and check exit code
#if defined(_WIN32)
    const auto rc = ::_pclose(pipe);
#else
    const auto rc = ::pclose(pipe);
#endif
    return rc == 0;
}

int execute_unix2(const std::string& cmd) {
#if ASTAP_NO_SYSTEM
	(void)cmd;
	return -1;
#else
	return std::system(cmd.c_str());
#endif
}

std::string get_short_path(const std::string& long_path) {
#if defined(_WIN32)
    // TODO: round-trip through GetShortPathNameW once the wide-string story
    // for the rest of the port is settled.
    return long_path;
#else
    // On non-Windows there is no notion of short paths.
    return long_path;
#endif
}

/// MARK: - Filesystem Helpers

void delete_files(const std::string& lpath, const std::string& file_spec) {
    // Translate the mask (e.g. "*.wcs") into a simple suffix match.
    std::error_code ec;
    if (!fs::is_directory(lpath, ec)) {
        return;
    }
    
    // Strip any leading '*' so "*.wcs" becomes ".wcs"
    auto suffix = file_spec;
    if (!suffix.empty() && suffix.front() == '*') {
        suffix.erase(0, 1);
    }
    
    for (const auto& entry : fs::directory_iterator(lpath, ec)) {
        if (ec) {
            break;
        }
        if (!entry.is_regular_file()) {
            continue;
        }
        const auto name = entry.path().filename().string();
        if (suffix.empty() ||
            (name.size() >= suffix.size() &&
             name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)) {
            auto rmec = std::error_code{};
            fs::remove(entry.path(), rmec);
        }
    }
}

std::int64_t textfile_size(const std::string& name) noexcept {
    auto ec = std::error_code{};
    const auto sz = fs::file_size(name, ec);
    if (ec) {
        return 0;
    }
    return static_cast<std::int64_t>(sz);
}

void wait(std::chrono::milliseconds duration) {
    // Sleep without GUI event polling.
    std::this_thread::sleep_for(duration);
}

void log_to_file(const std::string& logf, const std::string& mess) {
    auto f = std::ofstream(logf, std::ios::app);
    if (f) {
        f << mess << '\n';
    }
}

void log_to_file2(const std::string& logf, const std::string& mess) {
    // Truncating write.
    auto f = std::ofstream(logf, std::ios::trunc);
    if (f) {
        f << mess << '\n';
    }
}

/// MARK: - Progress / Status Sink

void progress_indicator(double percent, std::string_view info) {
    // Headless sink: emit to stderr with a carriage return so successive
    // updates overwrite the same line when the stream is a TTY.
    if (percent < 0.0) {
        // Sentinel used to reset ("done") the indicator.
        std::cerr << '\r' << std::string(60, ' ') << '\r' << std::flush;
        return;
    }
    const auto pct = std::clamp(percent, 0.0, 100.0);
    std::cerr << std::format("\r[{:5.1f}%] {}", pct, info) << std::flush;
    if (pct >= 100.0) {
        std::cerr << '\n';
    }
}

/// MARK: - PlateSolve2

bool platesolve2_command() {
    // TODO: faithful port. The implementation relies heavily on GUI singletons
    // and global state. See the detailed roadmap in the original porting notes.
    return false;
}

/// MARK: - IniFile

namespace {

[[nodiscard]] std::string trim(std::string s) {
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

[[nodiscard]] bool to_bool(const std::string& v, bool def) {
    if (v.empty()) {
        return def;
    }
    // The original writes "1"/"0" by default but accepts true/false too.
    if (v == "1" || v == "true" || v == "TRUE" || v == "True") {
        return true;
    }
    if (v == "0" || v == "false" || v == "FALSE" || v == "False") {
        return false;
    }
    return def;
}
 
}  // anonymous namespace

IniFile::IniFile(const std::string& path) {
    (void)load(path);
}

bool IniFile::load(const std::string& path) {
    _sections.clear();
    auto in = std::ifstream(path);
    if (!in) {
        return false;
    }
    
    auto line = std::string{};
    auto current_section = std::string{};
    while (std::getline(in, line)) {
        // Strip CR for CRLF files.
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        auto t = trim(line);
        if (t.empty() || t.front() == ';' || t.front() == '#') {
            continue;
        }
        
        if (t.front() == '[' && t.back() == ']') {
            current_section = t.substr(1, t.size() - 2);
            // Ensure section exists
            [[maybe_unused]] auto& sec = _sections[current_section];
            continue;
        }
        
        const auto eq = t.find('=');
        if (eq == std::string::npos) {
            continue;
        }
        auto key = trim(t.substr(0, eq));
        auto val = trim(t.substr(eq + 1));
        _sections[current_section][key] = val;
    }
    return true;
}

bool IniFile::save(const std::string& path) const {
    auto out = std::ofstream(path, std::ios::trunc);
    if (!out) {
        return false;
    }
    for (const auto& [section, kv] : _sections) {
        out << '[' << section << "]\n";
        for (const auto& [key, value] : kv) {
            out << key << '=' << value << '\n';
        }
        out << '\n';
    }
    return true;
}

void IniFile::clear() { _sections.clear(); }

bool IniFile::has(const std::string& section, const std::string& key) const noexcept {
    auto sit = _sections.find(section);
    if (sit == _sections.end()) {
        return false;
    }
    return sit->second.find(key) != sit->second.end();
}

std::string IniFile::read_string(const std::string& section,
                                 const std::string& key,
                                 const std::string& def) const {
    auto sit = _sections.find(section);
    if (sit == _sections.end()) {
        return def;
    }
    auto kit = sit->second.find(key);
    if (kit == sit->second.end()) {
        return def;
    }
    return kit->second;
}

int IniFile::read_int(const std::string& section, const std::string& key, int def) const {
    const auto v = read_string(section, key, {});
    if (v.empty()) {
        return def;
    }
    try {
        return std::stoi(v);
    } catch (...) {
        return def;
    }
}

bool IniFile::read_bool(const std::string& section, const std::string& key, bool def) const {
    return to_bool(read_string(section, key, {}), def);
}

double IniFile::read_float(const std::string& section, const std::string& key, double def) const {
    const auto v = read_string(section, key, {});
    if (v.empty()) {
        return def;
    }
    try {
        return std::stod(v);
    } catch (...) {
        return def;
    }
}

void IniFile::write_string(const std::string& section, const std::string& key, const std::string& value) {
    _sections[section][key] = value;
}

void IniFile::write_int(const std::string& section, const std::string& key, int value) {
    _sections[section][key] = std::to_string(value);
}

void IniFile::write_bool(const std::string& section, const std::string& key, bool value) {
    _sections[section][key] = value ? "1" : "0";
}

void IniFile::write_float(const std::string& section, const std::string& key, double value) {
    auto oss = std::ostringstream{};
    oss << value;
    _sections[section][key] = oss.str();
}

/// MARK: - load_settings / save_settings

bool load_settings(const std::string& path) {
    auto ini = IniFile{};
    if (!ini.load(path)) {
        return false;
    }
    
    // Sentinel: if window_left is absent, treat as a fresh install.
    if (!ini.has("main", "window_left")) {
        // TODO: reset window position to (0, 0).
        return false;
    }
    
    // ----- [main] window geometry ------------------------------------------
    // TODO: restore window geometry from ini "main" section.
    
    // ----- [main] font ------------------------------------------------------
    // TODO: font_color, font_size, font_name2, font_style, font_charset, pedestal
    
    // ----- [main] stretch / range / saturation -----------------------------
    // TODO: minimum_position, maximum_position, range, saturation_factor,
    //       polynomial, thumbnails_width, thumbnails_height
    
    // ----- [main] viewer toggles -------------------------------------------
    // TODO: inversemousewheel, fliphorizontal, flipvertical, annotations,
    //       north_east, star_profile, mount_position, constellations,
    //       grid, grid_az, pos_date, freetxt, f_text
    
    // ----- [main] photometry display ---------------------------------------
    // TODO: noise_e, egain_d, egain_ext, add_marker, preview_demosaic,
    //       s_overwrite, maintain_date, add_lim_magn
    
    // ----- [main] textboxes / paths ----------------------------------------
    // TODO: marker_position, ra, dec, gamma, last_file, export_index,
    //       anno_magn, cal_batch
    
    // ----- [ast] asteroids/comets ------------------------------------------
    // TODO: mpcorb_path, cometels_path, maxcount_asteroid, maxmag_asteroid,
    //       font_follows_diameter, showfullnames, showmagnitude, add_date,
    //       lat_default, long_default, annotation_color, annotation_diameter,
    //       add_annotations
    
    // ----- [anet] astrometry.net -------------------------------------------
    // TODO: astrometry_extra_options, show_console, cygwin_path
    
    // ----- [sqm] -----------------------------------------------------------
    // TODO: sqm_applyDF = ini.read_bool("sqm", "apply_df", false);
    
    // ----- [main] recent files ---------------------------------------------
    // TODO: recent_files population
    
    // ----- [stack] window geometry / tabs ----------------------------------
    // TODO: stackmenu geometry, mosaic_crop, stack_method, box_blur_factor,
    //       stack_tab, demosaic_method2, conv_progr
    
    // ----- [stack] OSC color -----------------------------------------------
    // TODO: osc_color_convert, osc_al, osc_cs, osc_pr, osc_cw, osc_sd
    
    // ----- [stack] smoothing -----------------------------------------------
    // TODO: smooth_dia, smooth_stars
    
    // ----- [stack] keys ----------------------------------------------------
    // TODO: sqm_key, centaz_key
    
    // ----- [stack] LRGB -----------------------------------------------------
    // TODO: lrgb_al, green_fl, lrgb_cs, lrgb_pr, lrgb_sm,
    //       lrgb_smd, lrgb_sms, lrgb_sw, lrgb_sd
    
    // ----- [stack] mosaic / classify ---------------------------------------
    // TODO: ignore_header_solution, equalise_background, merge_overlap,
    //       limit_back_corr, classify_*, add_time, copy_sett, uncheck_outliers,
    //       blur_factor
    
    // ----- [stack] alignment method ----------------------------------------
    // TODO: align_method "1".."4"
    
    // ----- [stack] solver / detection --------------------------------------
    // TODO: write_log, align_blink, time_stamp, force_slow, use_triples, sip,
    //       star_database, solve_search_field, radius_search, quad_tolerance,
    //       maximum_stars, min_star_size, min_star_size_stacking,
    //       manual_centering, downsample, sd_factor, most_common_filter_radius,
    //       extract_background_box_size, dark_areas_box_size,
    //       ring_equalise_factor, gradient_filter_factor
    
    // ----- [stack] bayer + per-filter strings ------------------------------
    // TODO: bayer_pat, red_filter1/2, green_filter1/2, blue_filter1/2,
    //       luminance_filter1/2
    
    // ----- [stack] colour mixing -------------------------------------------
    // TODO: rr/rg/rb/gr/gg/gb/br/bg/bb factors, filter_add, add_value,
    //       multiply, smart_smooth_width, star_level_colouring,
    //       filter_artificial_colouring, resize_factor
    
    // ----- [stack] photometry tab ------------------------------------------
    // TODO: nr_stars_p, flux_aperture, annulus_radius, font_size_p,
    //       annotate_m, reference_d, measure_all, ign_saturation
    
    // ----- [stack] decolour / noise / hue ----------------------------------
    // TODO: sigma_decolour, sd_factor_list, noisefilter_blur, noisefilter_sd,
    //       hue_fuzziness, saturation_tolerance, blend, sample_size,
    //       usm_amount, usm_radius, usm_thresh
    
    // ----- [stack] mount / video / contour ---------------------------------
    // TODO: wcs, video_index, frame_rate, contour_gaus, contour_sd,
    //       contour_grid, groupsize
    
    // ----- [aavso] photometry report ---------------------------------------
    // TODO: obscode, delim_pos, baa_style, sort_alphabetically, hjd_date,
    //       ensemble_database, pfilter, slope, vsp_stars
    
    // ----- [live] live stacking / monitor ----------------------------------
    // TODO: live_stack_dir, monitor_dir, write_jpeg, to_clipboard,
    //       live_inspect, monitor_df
    
    // ----- [insp] inspector window -----------------------------------------
    // TODO: insp_left, insp_top, insp_angle, contour, voronoi, values,
    //       vectors, 3corners, extra_stars, insp_binning, insp_grid, insp_grad
    
    // ----- [files] listview contents ---------------------------------------
    // TODO: port listview population once the listview model lands.
    
    // TODO: restore stack menu visibility.
    
    return true;
}

void save_settings(const std::string& path) {
    auto ini = IniFile{};
    ini.clear();  // start from a blank file
    
    // ----- [main] window geometry ------------------------------------------
    // TODO: ini.write_int("main", "window_left/top/height/width", ...);
    
    // ----- [main] font ------------------------------------------------------
    // TODO: font_color, font_size, font_name2, font_style, font_charset, pedestal
    
    // ----- [main] stretch / range / saturation -----------------------------
    // TODO: minimum_position, maximum_position, range, saturation_factor,
    //       polynomial, thumbnails_width, thumbnails_height
    
    // ----- [main] viewer toggles -------------------------------------------
    // TODO: inversemousewheel, fliphorizontal, flipvertical, annotations,
    //       north_east, star_profile, mount_position, constellations,
    //       grid, grid_az, pos_date, freetxt, f_text
    
    // ----- [main] photometry display ---------------------------------------
    // TODO: noise_e, egain_d, egain_ext, add_marker, preview_demosaic,
    //       s_overwrite, maintain_date, add_lim_magn, ra, dec, gamma,
    //       marker_position, last_file, export_index, anno_magn, cal_batch
    
    // ----- [ast] asteroids/comets ------------------------------------------
    // TODO: mpcorb_path, cometels_path, maxcount, maxmag, font_follows,
    //       showfullnames, showmagnitude, add_date, lat/long_default,
    //       annotation_color, annotation_diameter, add_annotations
    
    // ----- [anet] / [sqm] --------------------------------------------------
    // TODO: cygwin_path, show_console, astrometry_extra_options, apply_df
    
    // ----- [main] recent files ---------------------------------------------
    // TODO: recent_files persistence
    
    // ----- [stack] window geometry / tabs ----------------------------------
    // TODO: stackmenu_visible, stackmenu geometry, stack_method, mosaic_crop,
    //       box_blur_factor, stack_tab, bayer_pat, demosaic_method2, conv_progr
    
    // ----- [stack] OSC ------------------------------------------------------
    // NB: load uses key "osc_cw" but save uses "osc_sw" -- preserve this
    //     asymmetry verbatim so existing .cfg files remain compatible.
    // TODO: osc_color_convert, osc_al, osc_cs, osc_pr, osc_sw, osc_sd,
    //       smooth_dia, smooth_stars
    
    // ----- [stack] keys ----------------------------------------------------
    // TODO: sqm_key + "*" sentinel
    
    // ----- [stack] LRGB ----------------------------------------------------
    // TODO: lrgb_al, green_fl, lrgb_cs, lrgb_pr, lrgb_sm,
    //       lrgb_smd, lrgb_sms, lrgb_sw, lrgb_sd
    
    // ----- [stack] mosaic / classify ---------------------------------------
    // TODO: ignore_header_solution, equalise_background, merge_overlap,
    //       limit_back_corr, classify_*, add_time, copy_sett, uncheck_outliers,
    //       write_log, align_blink, time_stamp, force_slow, use_triples, sip
    
    // ----- [stack] alignment method ----------------------------------------
    // TODO: write the chosen "1".."4"
    
    // ----- [stack] solver / detection --------------------------------------
    // TODO: star_database, solve_search_field, radius_search, quad_tolerance,
    //       maximum_stars, min_star_size, min_star_size_stacking,
    //       manual_centering, downsample, sd_factor, blur_factor,
    //       most_common_filter_radius, extract_background_box_size,
    //       dark_areas_box_size, ring_equalise_factor, gradient_filter_factor
    
    // ----- [stack] colour mixing -------------------------------------------
    // TODO: red_filter1/2, green_filter1/2, blue_filter1/2, luminance_filter1/2,
    //       rr/rg/rb/gr/gg/gb/br/bg/bb factors, filter_add, add_value,
    //       multiply, smart_smooth_width, star_level_colouring,
    //       filter_artificial_colouring, resize_factor
    
    // ----- [stack] photometry tab ------------------------------------------
    // TODO: nr_stars_p, flux_aperture, annulus_radius, font_size_p,
    //       annotate_m, reference_d, measure_all, ign_saturation
    
    // ----- [stack] decolour / noise / hue ----------------------------------
    // TODO: sigma_decolour, sd_factor_list, noisefilter_blur, noisefilter_sd,
    //       hue_fuzziness, saturation_tolerance, blend,
    //       usm_amount, usm_radius, usm_thresh, sample_size
    
    // ----- [stack] mount / video / contour ---------------------------------
    // TODO: wcs, video_index, frame_rate, contour_gaus, contour_sd,
    //       contour_grid, groupsize
    
    // ----- [aavso] ---------------------------------------------------------
    // TODO: obscode, delim_pos, baa_style, sort_alphabetically, hjd_date,
    //       ensemble, pfilter, slope, vsp_stars
    
    // ----- [live] ----------------------------------------------------------
    // TODO: live_stack_dir, monitor_dir, write_jpeg, to_clipboard,
    //       live_inspect, monitor_df
    
    // ----- [insp] ----------------------------------------------------------
    // TODO: insp_left, insp_top, insp_angle, contour, voronoi, values,
    //       vectors, 3corners, extra_stars, insp_binning, insp_grid, insp_grad
    
    // ----- [files] listview contents ---------------------------------------
    // TODO: port once the listview model is available.
    
    (void)ini.save(path);
}

/// MARK: - write_ini

bool write_ini(const std::filesystem::path& filen,
               bool                         solution,
               double                       hfd_median,
               int                          hfd_counter) {
    auto path = filen;
    path.replace_extension(".ini");
    auto f = std::ofstream(path, std::ios::trunc);
    if (!f) {
        return false;
    }
    
    // Scientific notation formatter matching the original floattostrE
    auto e = [](double v) { return std::format("{:e}", v); };
    
    const auto& head = astap::head;
    
    if (solution) {
        f << "PLTSOLVD=T\n";
        f << "CRPIX1=" << e(head.crpix1) << '\n';
        f << "CRPIX2=" << e(head.crpix2) << '\n';
        f << "CRVAL1=" << e(head.ra0  * 180.0 / std::numbers::pi) << '\n';
        f << "CRVAL2=" << e(head.dec0 * 180.0 / std::numbers::pi) << '\n';
        f << "CDELT1=" << e(head.cdelt1) << '\n';
        f << "CDELT2=" << e(head.cdelt2) << '\n';
        f << "CROTA1=" << e(head.crota1) << '\n';
        f << "CROTA2=" << e(head.crota2) << '\n';
        f << "CD1_1="  << e(head.cd1_1)  << '\n';
        f << "CD1_2="  << e(head.cd1_2)  << '\n';
        f << "CD2_1="  << e(head.cd2_1)  << '\n';
        f << "CD2_2="  << e(head.cd2_2)  << '\n';
        
        // sqmfloat lives in core/sqm.cpp if the -sqm path was exercised.
        extern double sqmfloat;
        if (sqmfloat > 0.0) {
            f << "SQM=" << e(sqmfloat) << '\n';
        }
        
        if (hfd_median > 0.0) {
            f << "HFD=" << e(hfd_median) << '\n';
        }
        if (hfd_counter > 0) {
            f << "STARS=" << e(static_cast<double>(hfd_counter)) << '\n';
        }
    } else {
        f << "PLTSOLVD=F\n";
    }
    
    f << "CMDLINE=" << astap::cmdline << '\n';
    f << "DIMENSIONS=" << head.width << " x " << head.height << '\n';
    
    switch (astap::errorlevel) {
        case 2:  f << "ERROR=Not enough stars.\n"; break;
        case 16: f << "ERROR=Error reading image file.\n"; break;
        case 32: f << "ERROR=No star database found.\n"; break;
        case 33: f << "ERROR=Error reading star database.\n"; break;
        default: break;
    }
    
    if (!astap::warning_str.empty()) {
        f << "WARNING=" << astap::warning_str << '\n';
    }
    return static_cast<bool>(f);
}

/// MARK: - Cancellation Signal Handler

namespace {

#if defined(_WIN32)
BOOL WINAPI win_console_ctrl_handler(DWORD ctrl_type) {
    switch (ctrl_type) {
        case CTRL_C_EVENT:
        case CTRL_BREAK_EVENT:
        case CTRL_CLOSE_EVENT:
            astap::esc_pressed.store(true, std::memory_order_relaxed);
            return TRUE;  // handled; don't terminate immediately.
        default:
            return FALSE;
    }
}
#else
extern "C" void posix_cancel_signal([[maybe_unused]] int sig) {
    astap::esc_pressed.store(true, std::memory_order_relaxed);
}
#endif
 
}  // anonymous namespace

bool install_cancel_signal_handler() noexcept {
#if defined(_WIN32)
    return ::SetConsoleCtrlHandler(&win_console_ctrl_handler, TRUE) != 0;
#else
    struct sigaction sa{};
    sa.sa_handler = &posix_cancel_signal;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags   = SA_RESTART;
    const auto ok_int  = ::sigaction(SIGINT,  &sa, nullptr) == 0;
    const auto ok_term = ::sigaction(SIGTERM, &sa, nullptr) == 0;
    return ok_int && ok_term;
#endif
}
 
} // namespace
