///----------------------------------------
///      @file from_chars_fp.h
///   @ingroup ASTAP++
///     @brief Portable floating-point std::from_chars.
///    @author John Stephen. Mozilla Public License 2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///----------------------------------------

#pragma once

#include <cerrno>
#include <charconv>
#include <cstddef>
#include <cstdlib>
#include <string>

///----------------------------------------
namespace astap::core {
///----------------------------------------

///----------------------------------------
/// @brief Portable floating-point @c std::from_chars.
/// @details Apple's libc++ (through Xcode 16.x) leaves the floating-point
///          @c std::from_chars overloads deleted, so only the integer ones resolve
///          and a float/double call fails to compile on that toolchain (while it
///          builds on newer libc++, libstdc++ and MSVC). This strtod-based helper
///          honours the same @c [first,last) + @c from_chars_result contract on every
///          toolchain, matching the general-format quirks (rejects leading whitespace
///          and a leading '+', accepts '-', reports the past-the-number position and
///          @c errc). Locale-independent for '.'-decimal input. Integer parsing keeps
///          using @c std::from_chars, which every supported libc++ implements.
///
///          Header-only @c inline so any translation unit can call it without a shared
///          definition, and declared in its own minimal header so callers do not drag
///          in the rest of util.h (several units keep file-local parse helpers whose
///          names would otherwise clash).
/// @param first Start of the character range.
/// @param last  One past the end of the range.
/// @param[out] value Parsed value on success (untouched on failure).
/// @return @c from_chars_result with @c ptr past the consumed characters and @c ec.
///----------------------------------------

[[nodiscard]] inline std::from_chars_result from_chars(const char* first, const char* last,
                                                       double& value) noexcept {
    if (first == last) {
        return {first, std::errc::invalid_argument};
    }
    // std::from_chars (general format) rejects leading whitespace and a leading '+';
    // strtod would accept them, so reject up front to preserve the exact contract.
    const char c = *first;
    if (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f'
        || c == '\v' || c == '+') {
        return {first, std::errc::invalid_argument};
    }
    // strtod needs a null-terminated buffer; copy the (small) field.
    const std::string buf(first, last);
    const char* base = buf.c_str();
    char* end = nullptr;
    errno = 0;
    const double v = std::strtod(base, &end);
    if (end == base) {
        return {first, std::errc::invalid_argument};
    }
    const auto consumed = static_cast<std::size_t>(end - base);
    if (errno == ERANGE) {
        return {first + consumed, std::errc::result_out_of_range};
    }
    value = v;
    return {first + consumed, std::errc{}};
}

} // namespace astap::core
