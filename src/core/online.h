///----------------------------------------
///      @file online.h
///   @ingroup ASTAP++
///     @brief Online catalog query helpers for AAVSO VSP/VSX, Simbad, and Vizier.
///   @details Covers AAVSO VSP/VSX comparison-and-variable-star fetching,
///            variable-star annotation orchestration, FITS header ANNOTATE
///            lookup, and Simbad / Vizier fetch+parse. The GUI drawing half
///            is intentionally left out: the goal is data retrieval and
///            parsed in-memory arrays.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#pragma once

#include "../types.h"

#include <expected>
#include <functional>
#include <string>
#include <string_view>
#include <vector>

///----------------------------------------
namespace astap::core {
///----------------------------------------

using astap::IHttpClient;

///----------------------------------------
/// @brief Logger callback used by online flows for user-visible status
///        messages (e.g. "No VSX data!"). Empty by default = silent.
///----------------------------------------

using MessageHook = std::function<void(const std::string&)>;

/// MARK: - Data Structures

///----------------------------------------
/// @brief AAVSO VSP comparison-star record.
/// @details RA/DEC stored in radians; magnitudes are kept as text because
///          the source transports them verbatim into chart annotation strings
///          (the "?" placeholder must round-trip).
///----------------------------------------

struct Auid {
    std::string auid;
    double ra{};   ///< @brief Radians.
    double dec{};  ///< @brief Radians.
    std::string Vmag;
    std::string Verr;
    std::string Bmag;
    std::string Berr;
    std::string Rmag;
    std::string Rerr;
    std::string SGmag;
    std::string SGerr;
    std::string SRmag;
    std::string SRerr;
    std::string SImag;
    std::string SIerr;
};

///----------------------------------------
/// @brief AAVSO VSX variable-star record.
/// @details RA/DEC in radians, proper-motion corrected.
///----------------------------------------

struct VarStar {
    std::string name;
    double ra{};   ///< @brief Radians, proper-motion corrected.
    double dec{};  ///< @brief Radians, proper-motion corrected.
    std::string maxmag;
    std::string minmag;
    std::string period;
    std::string epoch;
    std::string category;
};

/// In-memory catalog cache. Kept as globals so callers can be ported
/// one at a time.
extern std::vector<Auid> vsp;
extern std::vector<VarStar> vsx;

/// MARK: - AAVSO Fetch

///----------------------------------------
/// @brief Fetch the AAVSO VSP comparison-star chart and (re)populate @c vsp.
/// @param http   HTTP client abstraction.
/// @param head   FITS header describing the field of view.
/// @param limiting_mag  Faintest magnitude to request.
/// @return @c true if the chart was downloaded and parsed; @c false on
///         network error or empty field.
///----------------------------------------

[[nodiscard]] bool download_vsp(IHttpClient& http, const Header& head, double limiting_mag);

///----------------------------------------
/// @brief Fetch the AAVSO VSX variable-star list and (re)populate @c vsx.
/// @param http             HTTP client abstraction.
/// @param head             FITS header describing the field of view.
/// @param limiting_mag     Faintest magnitude to request.
/// @param years_since_2000 Epoch offset for proper-motion correction.
/// @param period_filter    When @c true, drops entries with period == 0 or
///                         period >= 3.
/// @return @c true on success; @c false on network error or empty field.
///----------------------------------------

[[nodiscard]] bool download_vsx(IHttpClient& http,
                                const Header& head,
                                double limiting_mag,
                                double years_since_2000,
                                bool period_filter);
                                
///----------------------------------------
/// @brief Check whether the cached @c vsx is stale relative to @p head.
/// @param head  FITS header describing the current field of view.
/// @return @c true when the cache is empty or the first entry is farther
///         from the image centre than the image half-width.
///----------------------------------------

[[nodiscard]] bool aavso_update_required(const Header& head) noexcept;

///----------------------------------------
/// @brief Drive the variable-star annotation flow.
/// @details Decides between local and online catalogs based on
///          @p annotate_mode, calls download_vsx / download_vsp as needed.
///          On success the global @c vsx and @c vsp caches are populated
///          and the GUI overlay layer is responsible for projecting them.
///          The local-database modes (0..3) require catalog files that
///          ASTAP distributes separately; they are no-ops here.
/// @param http             HTTP client abstraction.
/// @param head             FITS header describing the field of view.
/// @param annotate_mode    Annotation mode selector (0-15).
/// @param years_since_2000 Epoch offset for proper-motion correction.
/// @param extract_visible  Whether to extract visible objects.
/// @param log              Optional logger; empty function = silent.
///----------------------------------------

void variable_star_annotation(IHttpClient& http,
                              const Header& head,
                              int annotate_mode,
                              double years_since_2000,
                              bool extract_visible,
                              MessageHook log = {});
                              
/// MARK: - FITS Annotation Lookup

///----------------------------------------
/// @brief Look up a named ANNOTATE record in the FITS header memo.
/// @details Walks @p memo_lines looking for an ANNOTATE record whose label
///          matches @p aname, then converts the centre pixel of that
///          annotation to celestial coordinates.
/// @param head       FITS header (WCS).
/// @param memo_lines The FITS header memo lines.
/// @param aname      Annotation name to search for.
/// @param[out] ra    Right ascension result (untouched if not found).
/// @param[out] dec   Declination result (untouched if not found).
///----------------------------------------

void annotation_position(const Header& head,
                         const std::vector<std::string>& memo_lines,
                         const std::string& aname,
                         double& ra,
                         double& dec);
                         
/// MARK: - Simbad

///----------------------------------------
/// @brief Which Simbad query to build with @ref make_simbad_url.
///----------------------------------------

enum class SimbadQuery {
    DeepSky,           ///< @brief All non-stellar objects in the field box.
    DeepSkyFiltered,   ///< @brief Objects matching @p maintype (e.g. "Galaxy").
    Stars,             ///< @brief Stars only.
    SingleAt,          ///< @brief Resolve one object at the field centre.
};

///----------------------------------------
/// @brief Build a Simbad sim-sam / sim-coo URL for the field of view in
///        @p head.
/// @param head      FITS header (uses ra0/dec0/width/height/cdelt1/cdelt2).
/// @param query     Which query to build.
/// @param maintype  Simbad maintype filter (only used for DeepSkyFiltered).
/// @return URL string suitable for an @c IHttpClient::get call.
///----------------------------------------

[[nodiscard]] std::string make_simbad_url(const Header& head,
                                          SimbadQuery query,
                                          std::string_view maintype = {});

///----------------------------------------
/// @brief Parsed Simbad object as produced by @ref plot_simbad.
///----------------------------------------

struct SimbadObject {
    long long ra_units{};   ///< @brief round(ra_hours * 864000 / 24).
    long long dec_units{};  ///< @brief round(sign * dec_deg * 324000 / 90).
    std::string name;
    std::string type;
    double magnitude{};     ///< @brief 0 means unknown.
    double size_arcmin{};   ///< @brief 0 means unknown (single-object branch only).
};

///----------------------------------------
/// @brief Parse a Simbad ASCII response.
/// @details Handles both the single-object reply and the pipe-delimited
///          list reply. Drawing onto the canvas is stubbed.
/// @param info  Raw Simbad response text.
/// @return Parsed Simbad objects.
///----------------------------------------

[[nodiscard]] std::vector<SimbadObject> plot_simbad(std::string_view info);

/// MARK: - Vizier

///----------------------------------------
/// @brief Build a Vizier asu-txt URL fetching Gaia DR3 rows for the field of
///        view in @p head.
/// @param head             FITS header (uses ra0/dec0/width/height/cdelt1/cdelt2).
/// @param limiting_mag     Faintest Gaia G magnitude to request.
/// @return URL string suitable for an @c IHttpClient::get call.
///----------------------------------------

[[nodiscard]] std::string make_vizier_gaia_url(const Header& head,
                                               double limiting_mag);

///----------------------------------------
/// @brief Parsed Vizier (Gaia) row.
///----------------------------------------

struct VizierObject {
    long long ra_units{};   ///< @brief round(ra_deg * 864000 / 360).
    long long dec_units{};  ///< @brief round(dec_deg * 324000 / 90).
    double magnitude{};     ///< @brief Already passed through transform_gaia.
};

///----------------------------------------
/// @brief Callback that converts (G, BP, RP) into the requested filter band.
/// @details Returning 0 marks the row as filtered out.
///----------------------------------------

using TransformGaiaFn = double(*)(std::string_view filter,
                                  double g, double bp, double rp);
                                  
///----------------------------------------
/// @brief Parse a Vizier asu-txt response.
/// @details Skips header until the "RA_ICRS" line plus the dashed separator,
///          then reads RA/DE/G/BP/RP rows. Drawing onto the canvas is stubbed.
/// @param info            Raw Vizier response text.
/// @param filter          Filter band identifier.
/// @param transform_gaia  Colour-transform callback.
/// @return Parsed Vizier rows.
///----------------------------------------

[[nodiscard]] std::vector<VizierObject> plot_vizier(std::string_view info,
                                                    std::string_view filter,
                                                    TransformGaiaFn transform_gaia);
    
} // namespace

