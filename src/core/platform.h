///----------------------------------------
///      @file platform.h
///   @ingroup ASTAP++
///     @brief Platform integration and settings persistence for ASTAP++.
///   @details Process-launch helpers, filesystem utilities, INI-backed settings
///            persistence, solver result output, and cancellation signal handling.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include "../types.h"

#include <chrono>
#include <cstdint>
#include <filesystem>
#include <map>
#include <string>
#include <string_view>
#include <vector>

///----------------------------------------
namespace astap::core {
///----------------------------------------

/// MARK: - Process Launch Helpers

///----------------------------------------
/// @brief Spawn @p cmd and block until the child exits.
/// @param cmd          Command line to execute.
/// @param show_console On Windows, controls SW_SHOWMINNOACTIVE vs SW_HIDE;
///                     ignored on POSIX.
/// @return @c true if the child exited with code 0.
///----------------------------------------

[[nodiscard]] bool execute_and_wait(const std::string& cmd, bool show_console);

///----------------------------------------
/// @brief Spawn @p executable with @p params, capturing stdout (merged with stderr).
/// @param executable      Path to the executable.
/// @param params          Arguments to pass.
/// @param show_output     Whether to capture output at all.
/// @param[out] captured_output  Collected stdout text.
/// @return @c true if the child exited with code 0.
///----------------------------------------

[[nodiscard]] bool execute_unix(const std::string& executable,
                                const std::vector<std::string>& params,
                                bool show_output,
                                std::string& captured_output);
                                
///----------------------------------------
/// @brief Fire-and-forget shell exec.
/// @param cmd  Shell command to run.
/// @return Raw exit code from @c std::system.
///----------------------------------------

[[nodiscard]] int execute_unix2(const std::string& cmd);

///----------------------------------------
/// @brief Return the short-path form of @p long_path (Windows only).
/// @details On non-Windows platforms the input is returned unchanged.
/// @param long_path  Full filesystem path.
/// @return Short path on Windows, or @p long_path on other platforms.
///----------------------------------------

[[nodiscard]] std::string get_short_path(const std::string& long_path);

/// MARK: - Filesystem Helpers

///----------------------------------------
/// @brief Delete every file in @p lpath matching @p file_spec (e.g. "*.wcs").
/// @param lpath      Directory to scan.
/// @param file_spec  Glob-style mask (leading-"*" suffix match).
///----------------------------------------

void delete_files(const std::string& lpath, const std::string& file_spec);

///----------------------------------------
/// @brief Return the size of @p name in bytes, or 0 if it does not exist.
/// @param name  Filesystem path to query.
/// @return File size in bytes, or 0 on error.
///----------------------------------------

[[nodiscard]] std::int64_t textfile_size(const std::string& name) noexcept;

///----------------------------------------
/// @brief Sleep for @p duration, without GUI event polling.
/// @param duration  How long to sleep (default 500 ms).
///----------------------------------------

void wait(std::chrono::milliseconds duration = std::chrono::milliseconds{500});

///----------------------------------------
/// @brief Append a single line to @p logf (creates the file if missing).
/// @param logf  Log file path.
/// @param mess  Message to append.
///----------------------------------------

void log_to_file(const std::string& logf, const std::string& mess);

///----------------------------------------
/// @brief Truncate @p logf and write a single line.
/// @param logf  Log file path.
/// @param mess  Message to write.
///----------------------------------------

void log_to_file2(const std::string& logf, const std::string& mess);

/// MARK: - Progress / Status Sink

///----------------------------------------
/// @brief Emit a progress update to stderr.
/// @details A negative @p percent resets the indicator. Values >= 100 emit
///          a trailing newline.
/// @param percent  Progress percentage (0-100, or negative to reset).
/// @param info     Descriptive text for the status line.
///----------------------------------------

void progress_indicator(double percent, std::string_view info);

/// MARK: - PlateSolve2

///----------------------------------------
/// @brief PlateSolve2 command-line entry point.
/// @details Reads sys argv, runs the solver, and writes the .apm file that
///          PlateSolve2 callers (SGP, APT) expect. Currently stubbed.
/// @return @c true on a successful solve.
///----------------------------------------

[[nodiscard]] bool platesolve2_command();

/// MARK: - INI-Backed Settings

///----------------------------------------
/// @class IniFile
/// @brief Tiny hand-rolled INI store for flat section.key=value files.
/// @details Flat section.key=value store with typed read/write helpers that
///          preserve the original semantics (default-on-missing, bool/int/float
///          formatting, etc.).
///----------------------------------------

class IniFile {
public:
    IniFile() = default;
    explicit IniFile(const std::string& path);
    
    ///----------------------------------------
    /// @brief Load settings from an INI file on disk.
    /// @param path  Filesystem path to the INI file.
    /// @return @c true if the file was readable.
    ///----------------------------------------
    
    [[nodiscard]] bool load(const std::string& path);
    
    ///----------------------------------------
    /// @brief Save all settings to an INI file on disk.
    /// @param path  Filesystem path to write.
    /// @return @c true on successful write.
    ///----------------------------------------
    
    [[nodiscard]] bool save(const std::string& path) const;
    
    /// @brief Remove all sections and keys.
    void clear();
    
    ///----------------------------------------
    /// @brief Check whether a key exists.
    /// @param section  INI section name.
    /// @param key      Key within the section.
    /// @return @c true if the key is present.
    ///----------------------------------------
    
    [[nodiscard]] bool has(const std::string& section, const std::string& key) const noexcept;
    
    ///----------------------------------------
    /// @brief Read a string value, returning @p def if missing.
    /// @param section  INI section name.
    /// @param key      Key within the section.
    /// @param def      Default value if key is absent.
    /// @return The stored string, or @p def.
    ///----------------------------------------
    
    [[nodiscard]] std::string read_string(const std::string& section,
                                          const std::string& key,
                                          const std::string& def = {}) const;
                                          
    ///----------------------------------------
    /// @brief Read an integer value, returning @p def if missing or unparseable.
    /// @param section  INI section name.
    /// @param key      Key within the section.
    /// @param def      Default value.
    /// @return The stored integer, or @p def.
    ///----------------------------------------
    
    [[nodiscard]] int read_int(const std::string& section,
                               const std::string& key,
                               int def = 0) const;
                               
    ///----------------------------------------
    /// @brief Read a boolean value, returning @p def if missing or unparseable.
    /// @param section  INI section name.
    /// @param key      Key within the section.
    /// @param def      Default value.
    /// @return The stored boolean, or @p def.
    ///----------------------------------------
    
    [[nodiscard]] bool read_bool(const std::string& section,
                                 const std::string& key,
                                 bool def = false) const;
                                 
    ///----------------------------------------
    /// @brief Read a floating-point value, returning @p def if missing or unparseable.
    /// @param section  INI section name.
    /// @param key      Key within the section.
    /// @param def      Default value.
    /// @return The stored double, or @p def.
    ///----------------------------------------
    
    [[nodiscard]] double read_float(const std::string& section,
                                    const std::string& key,
                                    double def = 0.0) const;
                                    
    /// @brief Write a string value.
    void write_string(const std::string& section, const std::string& key, const std::string& value);
    
    /// @brief Write an integer value.
    void write_int(const std::string& section, const std::string& key, int value);
    
    /// @brief Write a boolean value.
    void write_bool(const std::string& section, const std::string& key, bool value);
    
    /// @brief Write a floating-point value.
    void write_float(const std::string& section, const std::string& key, double value);
    
private:
    // section -> (key -> value). std::map keeps deterministic on-disk order.
    std::map<std::string, std::map<std::string, std::string>> _sections;
};

/// MARK: - Settings Load/Save

///----------------------------------------
/// @brief Load every persisted setting from @p path.
/// @details Returns @c true if the file was readable AND contained a valid
///          window_left key (the sentinel used to detect a fresh install).
/// @param path  Filesystem path to the settings INI file.
/// @return @c true on successful load with valid sentinel key.
///----------------------------------------

[[nodiscard]] bool load_settings(const std::string& path);

///----------------------------------------
/// @brief Persist every setting to @p path (overwrites existing contents).
/// @param path  Filesystem path to write.
///----------------------------------------

void save_settings(const std::string& path);

/// MARK: - Solver Result Output

///----------------------------------------
/// @brief Write the solver's .ini result file alongside @p filen.
/// @details Extension is replaced with .ini. When @p solution is @c true,
///          the WCS fields from the global header plus hfd/sqm/stars are
///          emitted; when @c false only a PLTSOLVD=F marker + cmdline +
///          error string are written.
/// @param filen       Base filepath (extension replaced with .ini).
/// @param solution    Whether a valid plate solution was found.
/// @param hfd_median  Median HFD value (0 to omit).
/// @param hfd_counter Number of stars used for HFD (0 to omit).
/// @return @c true on successful write.
///----------------------------------------

[[nodiscard]] bool write_ini(const std::filesystem::path& filen,
                             bool                         solution,
                             double                       hfd_median = 0.0,
                             int                          hfd_counter = 0);
                             
/// MARK: - Cancellation

///----------------------------------------
/// @brief Install a signal handler for SIGINT/SIGTERM that sets
///        @c astap::esc_pressed to @c true.
/// @details Idempotent; subsequent calls replace the previous handler with
///          the same behaviour.
/// @return @c true on success.
///----------------------------------------

[[nodiscard]] bool install_cancel_signal_handler() noexcept;
    
} // namespace
