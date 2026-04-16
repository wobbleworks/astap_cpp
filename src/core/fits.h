///----------------------------------------
///      @file fits.h
///   @ingroup ASTAP++
///     @brief FITS I/O and header-edit primitives for in-memory card images.
///   @details This module owns the binary read/write paths for FITS images, the
///            in-memory header-line ("memo") manipulation primitives that the
///            rest of ASTAP uses to stamp WCS and instrumentation keywords into
///            a header before saving, and a few related helpers (WCS form
///            conversion, cfitsio shell-out wrappers, an experimental Rice
///            encoder, and the FITS-to-BMP convenience wrapper).
///
///            "memo" is a std::vector<std::string> where each entry is one
///            80-char card image of a FITS header. The header-edit primitives
///            operate on this vector, NOT on streams. Actual disk writes happen
///            in save_fits / savefits_update_header which serialise the memo
///            into 2880-byte blocks of space-padded ASCII.
///    @author Ported from Han Kleijn's ASTAP; MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include <cstdint>
#include <expected>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "../types.h"

///----------------------------------------
namespace astap::core {
///----------------------------------------

/// MARK: Constants

constexpr auto kFitsBlockSize    = std::size_t{2880};
constexpr auto kFitsCardSize     = std::size_t{80};
constexpr auto kFitsCardsPerBlock = kFitsBlockSize / kFitsCardSize;  // 36
constexpr auto kFitsBufWide      = std::size_t{1024 * 120};

/// MARK: Forward declarations (stubs in sister translation units)

///----------------------------------------
/// @brief Run a child process and wait for completion.
/// @param cmd Shell command to execute.
/// @param show_console Whether to show a console window (platform-specific).
/// @return True on successful execution.
///----------------------------------------

[[nodiscard]] bool execute_and_wait(const std::string& cmd, bool show_console);

///----------------------------------------
/// @brief Append a path to the recent-files list.
/// @param path File path to add.
///----------------------------------------

void add_recent_file(const std::string& path);

///----------------------------------------
/// @brief Delete files matching a pattern from a directory.
/// @param dir Directory to scan.
/// @param spec Glob/wildcard specification.
///----------------------------------------

void delete_files(const std::filesystem::path& dir, std::string_view spec);

///----------------------------------------
/// @brief Pump the message loop for a short delay.
/// @param ms Delay in milliseconds (default 500).
///----------------------------------------

void wait_ms(float ms = 500.0f);

///----------------------------------------
/// @brief Update the progress indicator.
/// @param pct Percentage 0..100 (-100 resets, -101 marks failure).
/// @param info Optional status string.
///----------------------------------------

void progress_indicator(double pct, const std::string& info);

/// MARK: Header loading

///----------------------------------------
/// @brief Load a FITS image (or Astro-TIFF passthrough).
/// @details When @p light is true the loader populates WCS / mount /
///          acquisition fields; otherwise only the minimum needed to treat the
///          file as a dark/flat is read. When @p load_data is false only the
///          header is parsed. When @p update_memo is true the memo is cleared
///          and refilled with the parsed header card lines. @p get_ext chooses
///          which header to land on (0 = primary, 1 = first extension, ...).
/// @param filen Path to the FITS file.
/// @param light True for full WCS / mount / acquisition parsing.
/// @param load_data True to read pixel data, false for header-only.
/// @param update_memo True to clear and refill @p memo with card lines.
/// @param get_ext HDU index (0 = primary).
/// @param[in,out] memo In-memory header card vector.
/// @param[out] head Parsed FITS header struct.
/// @param[out] img_loaded2 Decoded image data.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool load_fits(const std::filesystem::path& filen,
                             bool light,
                             bool load_data,
                             bool update_memo,
                             int get_ext,
                             std::vector<std::string>& memo,
                             astap::Header& head,
                             astap::ImageArray& img_loaded2);
                             
/// MARK: Header saving

///----------------------------------------
/// @brief Save an image array to a FITS file.
/// @details @p type1 selects bit depth: 8, 16, -32, or 24 (special: 8-bit
///          interleaved RGB stored as NAXIS1=3). When @p override2 is false, any
///          existing @p filen2 triggers an overwrite confirmation.
/// @param img Source image data.
/// @param[in,out] memo Header card vector (may be updated with BITPIX, etc.).
/// @param filen2 Destination file path.
/// @param type1 Bit depth selector.
/// @param override2 True to silently overwrite.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool save_fits(const astap::ImageArray& img,
                             std::vector<std::string>& memo,
                             const std::filesystem::path& filen2,
                             int type1,
                             bool override2);
                             
///----------------------------------------
/// @brief Rewrite an existing FITS file's header in-place.
/// @details The data section is copied verbatim, so the new memo MUST occupy
///          the same number of 2880-byte blocks as the original.
/// @param[in,out] memo New header card vector.
/// @param filen2 Path to the FITS file to update.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool savefits_update_header(std::vector<std::string>& memo,
                                         const std::filesystem::path& filen2);
                                         
/// MARK: State reset

///----------------------------------------
/// @brief Reset the global FITS state to defaults.
/// @param light When true, performs the wider reset path.
/// @param[out] head Header struct to reinitialise.
///----------------------------------------

void reset_fits_global_variables(bool light, astap::Header& head) noexcept;

/// MARK: TIFF-side header decoding

///----------------------------------------
/// @brief Decode a FITS-style header from a TIFF description tag.
/// @param light True for full WCS parsing.
/// @param[in,out] head Header struct to populate.
/// @param[in,out] memo Card lines parsed from the description.
///----------------------------------------

void read_keys_memo(bool light, astap::Header& head,
                    std::vector<std::string>& memo);
                    
/// MARK: Header-edit primitives

///----------------------------------------
/// @brief Update or insert a floating-point keyword in the memo.
/// @param[in,out] memo Header card vector.
/// @param inpt Keyword with trailing '=' (e.g. "BITPIX  =").
/// @param comment1 FITS comment string.
/// @param preserve_comment True to keep existing comment text.
/// @param x Value to write.
///----------------------------------------

void update_float(std::vector<std::string>& memo,
                  std::string_view inpt, std::string_view comment1,
                  bool preserve_comment, double x);
                  
///----------------------------------------
/// @brief Update or insert an integer keyword in the memo.
/// @param[in,out] memo Header card vector.
/// @param inpt Keyword with trailing '='.
/// @param comment1 FITS comment string.
/// @param x Value to write.
///----------------------------------------

void update_integer(std::vector<std::string>& memo,
                    std::string_view inpt, std::string_view comment1,
                    int x);
                    
///----------------------------------------
/// @brief Insert an integer keyword (always appends, never updates).
/// @param[in,out] memo Header card vector.
/// @param inpt Keyword with trailing '='.
/// @param comment1 FITS comment string.
/// @param x Value to write.
///----------------------------------------

void add_integer(std::vector<std::string>& memo,
                 std::string_view inpt, std::string_view comment1, int x);
                 
///----------------------------------------
/// @brief Update or insert a generic keyword (HISTORY/COMMENT aware).
/// @param[in,out] memo Header card vector.
/// @param message_key Keyword name.
/// @param message_value Value string.
/// @param message_comment Comment string.
///----------------------------------------

void update_generic(std::vector<std::string>& memo,
                    std::string_view message_key,
                    std::string_view message_value,
                    std::string_view message_comment);
                    
///----------------------------------------
/// @brief Update or insert a text keyword.
/// @param[in,out] memo Header card vector.
/// @param inpt Keyword with trailing '='.
/// @param comment1 Value and comment string.
///----------------------------------------

void update_text(std::vector<std::string>& memo,
                 std::string_view inpt, std::string_view comment1);
                 
///----------------------------------------
/// @brief Update or insert a long-string keyword with CONTINUE cards.
/// @param[in,out] memo Header card vector.
/// @param inpt Keyword with trailing '='.
/// @param thestr String value (may exceed 68 characters).
///----------------------------------------

void update_longstr(std::vector<std::string>& memo,
                    std::string_view inpt, std::string_view thestr);
                    
///----------------------------------------
/// @brief Insert a text keyword (always appends, never updates).
/// @param[in,out] memo Header card vector.
/// @param inpt Keyword with trailing '='.
/// @param comment1 Value and comment string.
///----------------------------------------

void add_text(std::vector<std::string>& memo,
              std::string_view inpt, std::string_view comment1);
              
///----------------------------------------
/// @brief Insert a long COMMENT keyword, wrapping at 72 characters.
/// @param[in,out] memo Header card vector.
/// @param descrip Description text to wrap across COMMENT cards.
///----------------------------------------

void add_long_comment(std::vector<std::string>& memo,
                      std::string_view descrip);
                      
///----------------------------------------
/// @brief Remove a keyword from the memo.
/// @param[in,out] memo Header card vector.
/// @param inpt Keyword to search for.
/// @param all True to remove all occurrences, false for first only.
///----------------------------------------

void remove_key(std::vector<std::string>& memo,
                std::string_view inpt, bool all);
                
///----------------------------------------
/// @brief Strip WCS / SIP solution keywords from the memo.
/// @param[in,out] memo Header card vector.
/// @param keep_wcs True to preserve the CDx_y matrix.
/// @param[in,out] head Header struct (cd1_1 zeroed when removing WCS).
///----------------------------------------

void remove_solution(std::vector<std::string>& memo, bool keep_wcs,
                     astap::Header& head);
                     
/// MARK: WCS form conversion

///----------------------------------------
/// @brief Convert old-style (CDELT/CROTA) WCS to new-style (CD matrix).
/// @param[in,out] head Header with cdelt/crota fields to read and CD fields to write.
///----------------------------------------

void old_to_new_WCS(astap::Header& head) noexcept;

///----------------------------------------
/// @brief Convert new-style (CD matrix) WCS to old-style (CDELT/CROTA).
/// @param[in,out] head Header with CD fields to read and cdelt/crota fields to write.
///----------------------------------------

void new_to_old_WCS(astap::Header& head) noexcept;

/// MARK: cfitsio shell-out wrappers

///----------------------------------------
/// @brief Unpack a .fz compressed FITS file via funpack.
/// @param[in,out] filename3 Path to the .fz file; updated to the unpacked path on success.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool unpack_cfitsio(std::filesystem::path& filename3);

///----------------------------------------
/// @brief Compress a FITS file via fpack.
/// @param filename3 Path to the FITS file to compress.
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool pack_cfitsio(const std::filesystem::path& filename3);

/// MARK: Rice encoder

///----------------------------------------
/// @brief Experimental Rice encoder for 16-bit data.
/// @details Encodes @p inp using parameter @p k and given @p bitdepth.
/// @param inp Input 16-bit samples.
/// @param k Rice coding parameter.
/// @param bitdepth Bit depth of the output words.
/// @param[out] outp Encoded output buffer (same size as input).
/// @param[out] compressed_size Number of output words used.
/// @return False if the encoded stream would exceed the input length.
///----------------------------------------

[[nodiscard]] bool rice_encoding(const std::vector<std::uint16_t>& inp,
                                 int k, int bitdepth,
                                 std::vector<std::uint16_t>& outp,
                                 int& compressed_size);
                                 
/// MARK: Binning

///----------------------------------------
/// @brief Re-bin the current file by a factor of 2, 3, or 4.
/// @param binfactor Bin factor (2, 3, or 4).
/// @return True on success.
///----------------------------------------

[[nodiscard]] bool binX2X3_file(int binfactor);

/// MARK: FITS to BMP

///----------------------------------------
/// @brief Render the current FITS file as a BMP alongside it.
/// @param filen Path to the source FITS file.
///----------------------------------------

void FITS_BMP(const std::filesystem::path& filen);
              
} // namespace
