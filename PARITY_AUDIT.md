# ASTAP++ vs ASTAP — Full Parity Audit

Differential source audit of the C++23 port against the Pascal original
(`/Users/john/Documents/ASTAP`). Seven functional clusters read on both sides in
full. Scope = CLI + computational library; the Qt GUI is out of scope (the
original's LCL GUI has no parity obligation).

Goal restated: the port should match the original library's behavior exactly.

---

## Cross-cutting themes

1. **Stubbed computational units that silently return zero.** Several routines
   compile and are called on live paths but return `0`/`(0,0)`/`{0,0,0}`,
   producing wrong output with no error: `dsspos` (DSS WCS), `pixel_to_celestial`
   (analyse CSV RA/DEC), `bicubic_interpolate`, `mode`/`apply_most_common`, and
   the photometry-side `analyse_image`.

2. **The solver's own star detection was never given the five analyse-path
   fixes.** `star_align.cpp` still carries the pre-reconciliation constants
   (1.5·HFD mask, extra retry pass, SMedian-as-nth_element). The differential
   analyse fixes touched `stack.cpp`; the `-solve` detection path is a *separate*
   routine and is still divergent.

3. **Unwired config globals.** Values are read from config/globals that nothing
   ever populates, so whole modes are permanently off: `use_triples`,
   `process_as_osc`, and the `-ra` CLI hint.

4. **Whole feature families absent.** Master-calibration creation/selection,
   the CCD-inspector suite, point-source photometry/`MZERO` calibration,
   variable-star & spatial deep-sky overlays, satellite-streak `contour`, and
   the planetary/high-precision ephemeris.

5. **Clean-room reductions that trade accuracy or precision** (may be
   intentional): SLALIB ephemeris → low-precision Meeus; bicubic → bilinear;
   `extended` (80-bit) → `float`/`double`; `round` half-to-even → `lround`.

---

## Tier 1 — Active defects on the solve / CLI path (fix first)

> **STATUS: DONE (2026-07-11).** All items below fixed and verified against the
> oracle; `run_tests.sh` green (26/26), with two new differential solve cases
> guarding the `-ra`/`-spd` hint wiring (correct-hint solve + wrong-RA decline).
> The `bin_and_find_stars` `hfd_max` param is benign — the solve path passes
> `10.0`, matching Pascal's hardcoded `hfd1<=10`.

These affect the code paths the differential harness runs today.

- **[HIGH] Solver star-detection mask radius 1.5·HFD vs Pascal 3.0·HFD** —
  `star_align.cpp:836` vs `unit_star_align.pas:905`. Same bug already fixed in
  analyse; still live in the solve path. Changes star list → quads → every match.
- **[HIGH] Solver adds an extra 3·noise retry pass** — `star_align.cpp:804-809,880`
  vs `unit_star_align.pas:870-939`. Port drops to a lower detection floor than
  Pascal ever uses on sparse fields.
- **[HIGH] Solver adds a "centroid already masked" rejection Pascal lacks** —
  `star_align.cpp:846-849` vs `unit_star_align.pas:894-930`.
- **[HIGH] DB star-count floor `max(30, …)` in the port** — `astrometric_solving.cpp:1028-1030`
  vs `unit_astrometric_solving.pas:902-907`. Pascal lets `max_stars` drop below 30
  on tight fields; the port never does, changing requested density/faintest mag.
- **[HIGH] `find_fit` median = single `nth_element`, not Pascal SMedian 2/3-avg** —
  `star_align.cpp:290-292` vs `unit_star_align.pas:710` (`smedian`). Shifts which
  quad matches survive outlier rejection. (Same SMedian issue fixed for HFD.)
- **[HIGH] `use_triples` permanently dead** — read at `astrometric_solving.cpp:886`
  but never populated (`:173-183`), so the low-star-count (`nrstars<30`) triples
  fallback can never fire. `platform.cpp:530,634` confirm it's an unparsed TODO.
- **[HIGH] `-ra` parsed then discarded** — `main.cpp:340` (`/* TODO: push into
  head.ra0 */`) while `-spd` *is* applied (`:341`). A blind solve given both RA and
  Dec hints honors only Dec, corrupting the seed. Primary mount/hint path, dead.
- **[MED] `bin_and_find_stars` carries an extra `hfd_max` parameter** absent in
  Pascal (`unit_astrometric_solving.pas:479`) — signature divergence in the shared
  detection theme; verify it doesn't alter acceptance.

## Tier 2 — Stubs that silently produce zero/wrong output

> **STATUS: 5 of 6 DONE (2026-07-11).** `dsspos` (real DSS plate model, in
> fits.cpp next to its coefficients), `pixel_to_celestial` (stack.cpp now
> forwards to the real wcs.cpp projector → analyse CSV RA/Dec are real),
> photometry `analyse_image` (re-pointed to `astap::stacking::analyse_image`,
> mirroring the Pascal cross-unit call), `mode`/`apply_most_common` (real modal
> filter), and `bicubic_interpolate` + `resize_img_loaded` (real Catmull-Rom,
> matching the original). `run_tests.sh` 26/26. **Internal-astrometry stacking
> deferred to Tier 3:** `update_solution_and_save` hard-stubs `filename2` to an
> empty string (needs threading from astap_main) and sits inside the incomplete
> stacking pipeline, so wiring `solve_image`/save stubs now is net-negative and
> untestable — do it with the Tier 3 stacking work.

- **[HIGH] DSS plate model `dsspos` → `(0,0)`** — `link_stubs.cpp:105` vs whole
  `unit_dss.pas` (13-term polynomial + gnomonic). Called live at `wcs.cpp:304`
  for `formalism==2`; every DSS-survey plate pixel resolves to RA=0/Dec=0.
- **[HIGH] Photometry `analyse_image` resolves to the no-op stub** —
  `photometry.cpp:46-48` forward-decl binds to `link_stubs.cpp:156-166` (wrong
  namespace/signature vs the real `stack.cpp:1185`). `hfd_med` is always 0, so
  `calibrate_photometry`'s aperture/annulus sizing (`photometry.cpp:559`) never runs.
- **[HIGH] `pixel_to_celestial` → ra=dec=0** — `stack.cpp:84`. `analyse_image`
  CSV RA/DEC columns emit `0,0` (also breaks any internal-astrometry consumer).
- **[HIGH] Internal-astrometry stacking non-functional** — `solve_image`,
  `fits_file_name`, `savefits_update_header`, `save_tiff16` all stub to `false`
  (`stack.cpp:93-120,511`); every `use_astrometry_internal` alignment aborts.
- **[MED] `mode()` / `apply_most_common` → 0** — `stack.cpp:1380-1384` fills the
  image with zeros vs real `mode()` (`unit_stack.pas:3438,3472`).
- **[MED] `bicubic_interpolate` → {0,0,0} TODO; resize falls back to bilinear** —
  `stack.cpp:122,1086` vs `unit_interpolate.pas:107` (Catmull-Rom is the *used*
  method; port inverted the used/unused pair). Resized/rescaled output diverges.

## Tier 3 — Missing feature families (roadmap / scope decisions)

Calibration & stacking:
- **[HIGH] Master-frame CREATION absent** — no averaging of raw darks/flats/
  bias/flat-darks into masters (`unit_stack.pas:4595,4664,11235,11385`). Port only
  loads a single pre-existing master (`stack.cpp:716,724`).
- **[HIGH] `load_master_dark`/`load_master_flat` empty stubs** — `stack.cpp:142-151`
  vs `unit_stack.pas:11064,11168`. No exposure/temp/gain/dimension matching, no
  closest-JD selection, no compatibility warnings.
- **[HIGH] OSC flat path dead** — `process_as_osc` never set (`stack.cpp:140,608`);
  colour flats always take the mono-division branch → OSC mis-calibration.
- **[HIGH] Reference-frame quality selection absent** — no `put_best_quality_on_top`
  (`unit_stack.pas:11839`); port always uses frame [0].
- **[HIGH] "Calibration only" method (4 & 5) missing** — `unit_stack.pas:11740`.
- **[HIGH] Live stacking gaps** — colour correction, JD/EXPTIME/temp header
  accumulation (`update_header` stub), and master dark/flat analysis all absent
  (`unit_live_stacking.pas:165,278-389,331-336,170-171`).
- **[MED] Per-frame quality outlier rejection absent** (`list_remove_outliers`,
  `unit_stack.pas:1495`); mosaic bg-equalisation/overlap-merge/crop/SIP hardcoded
  off (`stack_routines.cpp:1204-1213,1334-1335,153`); LRGB colour-mix matrix
  hardwired to identity (`stack_routines.cpp:501-505`).
- **[MED] Richer `unit_monitoring.pas` not ported** — `report_delta` guidance and
  `monitor_action` dispatch missing; monitoring picks OLDEST vs Pascal NEWEST
  (`live_monitoring.cpp:69-95`) and always calibrates (no gate).

Analysis & photometry:
- **[HIGH] Point-source photometry / flux calibration absent** — the `MZERO`/
  `MZEROR`/`MZEROAPT`/`LIM_MAGN` pipeline (`plot_and_measure_stars`,
  `unit_annotation.pas:1927`) has no equivalent; `photometry_catalog.cpp:330-336`
  hardcodes extended-object mode (`mzero_radius=99`, `aperture=0`).
- **[HIGH] Online-Gaia calibration flow not ported** — `photometry_catalog.cpp:284-325`
  never fetches/converts; requires pre-populated list (`unit_annotation.pas:2069-2098`).
- **[HIGH] VSX/VSP variable-star & spatial deep-sky overlays absent** — only the
  name-search path is ported (`annotation.cpp:142`); FOV-culled overlay selection
  (`plot_vsx_vsp`, `read_deepsky('S')`, `extract_visible`) is not.
- **[HIGH] Entire CCD-inspector suite absent** — tilt/curvature/off-axis
  (`CCDinspector`, `CCDinspector_analyse`), Voronoi/contour HFD maps, background
  grids (`unit_inspector_plot.pas:130-714,775-897,1026-1185,1361-1512`).
- **[HIGH] Satellite-streak `contour()` detection absent** — TODO stub
  (`contour.h:64-68`) vs `unit_contour.pas:175-478`.
- **[MED] Annotation catalog loaders + multi-DB selection missing**
  (`load_deep/variable/variable_13/variable_15/hyperleda`, `unit_annotation.pas:1046-1131`);
  distortion measurement and artificial-star generation absent.

Ephemerides (see parity-vs-intent below):
- **[HIGH] High-precision Earth ephemeris `sla_EPV2` → 2-body Meeus Sun**
  (`ephemerides.cpp:484,365` vs `unit_ephemerides.pas:1885`).
- **[HIGH] Full planetary theory `sla_PLANET` (Mercury–Neptune) absent** — only a
  hard-coded 2-body Jupiter+Saturn stub (`ephemerides.cpp:419,440`).
- **[MED] Universal-variable propagation → Meeus conic branches**
  (`ephemerides.cpp:254` vs `unit_ephemerides.pas:1776`); no JSTAT status;
  asteroid a→q conversion couples to it (`asteroid.cpp:173`).

## Tier 4 — Output fields & numeric parity

- **[MED] Solve `.ini` omits `HFD=` and `STARS=`** — `main.cpp:488` passes counts
  defaulted to 0 (`platform.h:305-306`); Pascal emits both (`astap_main.pas:13081-13082`).
- **[MED] Unknown/extra CLI flag throws → exit code 2** ("not enough stars")
  (`main.cpp:284-285,621-624`); Pascal silently ignores unknown options. Misleads
  NINA/APT/Ekos passing benign flags. The in-code comment contradicts the `throw`.
- **[MED] `-extract2` no longer forces SIP** (`astap_main.pas:13600`) and only runs
  after a *successful* solve (`main.cpp:554`; Pascal runs it either way).
- **[MED] AAVSO report** — MERR drops the `1.4826·MAD` check-star floor
  (`aavso_report.cpp:127-130`); NOTES reduced to static text; `#SOFTWARE` writes
  `ASTAP++` not `ASTAP` (`aavso_report.cpp:51`); local comp-magnitude parse absent.
- **[MED] Saturation cutoff off-by-one** — port `datamax_org-1000` vs Pascal
  effective `-1001` (`photometry.cpp:606` vs `unit_annotation.pas:1977,2047`).
- **[MED] Tile-catalog gating counts out-of-frame/red stars** —
  `photometry_catalog.cpp:198` reaches the cap sooner than `unit_annotation.pas:1940-1947`.
- **[MED] `calibration_and_alignment` omits PEDESTAL/CRPIX/DARK_CNT/CALSTAT
   provenance headers** (`stack_routines.cpp:2604-2626`).
- **[LOW rollup] Numeric/precision parity:** `round` half-even vs `lround` 1-LSB
  (`tiff.cpp:235`, `xisf.cpp:476`, `online.cpp`); raster-rotate `float` weights &
  `long double` corners vs Pascal `extended` (~1-ULP, matters for byte-exact
  diffing); XBINNING/XPIXSZ header rescale on resize omitted; `box_blur` leaves
  `head.naxis3` stale; near-edge & +X-boundary rejections one pixel stricter.
- **[LOW rollup] CLI parsing edge cases:** `take_value` can swallow a following
  flag; `-analyse`/`-extract` with no trailing value throws instead of defaulting
  to 30; `-sip n` (explicit disable) unsupported; help text omits the DB path line;
  `-debug` falls through to a solve instead of exiting.
- **[LOW rollup] Format/RAW coverage:** live stacking omits ptx/pef/rw2/srw RAWs;
  no `download_file` (fetch-to-disk); `demosaic_advanced` post-processing TODO
  (auto-bg-level, apply_factors, smart_colour_smooth); XISF loader defers inline
  WCS derivation to callers.

## Parity-vs-intent — needs a decision, not a blind fix

These may be deliberate re-imaginings (per the re-imagine-don't-transliterate
doctrine) rather than defects. Decide target before touching:

- Ephemeris clean-room (SLALIB → Meeus, no planetary theory). Trades sub-km /
  sub-arcsec SLALIB accuracy for low-precision series. If parity is required for
  aberration/annotation, this is a large re-port; if "good enough" is acceptable,
  document the tolerance.
- Bicubic → bilinear resample inversion.
- `ASTAP++` vs `ASTAP` in the AAVSO `#SOFTWARE` header (submission identity).
- `extended`/80-bit reductions to `double`/`float` on Apple ARM (no 80-bit).

## Verified faithful (no disparity) — do not "fix"

Star DB decode (sizes 5/6, byte-identical), all DB format selection, Gaia
photometric transforms & VizieR URL; FITS read/write (BITPIX set, BZERO/BSCALE,
float rescale, `.fz`, NAXIS3 color); TIFF/AVI/YUV4MPEG2 writers & XISF decode
(byte-exact); quad geometry, `lsq_fit`, cubic/SIP `gauss` solver, CD-matrix/CROTA,
pole-crossing, spiral search, second max-accuracy solve loop; `HFD`,
`find_star_center`, `get_background`, `measure_magnitudes`, hot-pixel, noise/
electron conversions, flux→MZERO & limiting-mag formulas; SQM & HJD in full;
gaussian blur, raster-rotate geometry (all cases), sigma-clip/comet combine, SIP
`calc_newx_newy`, dark/flat arithmetic; constellation figures, aberration/
nutation/precession, hyperbola focus-fit, image-sharpness metric.
