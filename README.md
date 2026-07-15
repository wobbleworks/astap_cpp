# ASTAP++

A C++23 port of [ASTAP](https://www.hnsky.org/astap.htm) — the
Astrometric STAcking Program by Han Kleijn. ASTAP solves astronomical
images against star catalogs, stacks exposures, and measures photometry.
ASTAP++ brings the core engine to modern C++ as a portable static
library with a command-line front end and an optional Qt 6 desktop GUI.
This library has been built and is tested against Linux, Windows, and
all Apple platforms.

## Status

Port is functionally complete. The static library, the `astap` CLI, and
the Qt 6 desktop GUI build and run on Linux, Windows, and macOS; a
22-suite doctest harness covers the pure-math modules (solving, stacking,
photometry, ephemerides, WCS, SIP, annotation, inspector, etc.). The
original Lazarus GUI is reimplemented as a new Qt front end rather than
ported line-for-line.

## Capabilities

- **Plate solving** — astrometric calibration via quad-matching, affine
  fits, and SIP polynomial distortion.
- **Image stacking** — LRGB, average, sigma-clip, mosaic, comet.
- **Bayer demosaicing** — nine variants including bilinear,
  astroSimple, astroC, astroM, superpixel, and X-Trans.
- **Photometry** — HFD, sky-quality measurement, Bortle scale, airmass.
- **Ephemerides** — Earth heliocentric and barycentric state, unified
  orbital-element propagator for asteroids and comets (elliptic /
  parabolic / hyperbolic), IAU 1976 precession, nutation, and annual
  aberration. Implemented from Meeus, *Astronomical Algorithms* (Ch. 21,
  25, 30, 33).
- **WCS** — pixel↔celestial transforms with tangent-plane projection
  and SIP distortion.
- **Catalog support** — Gaia DR3 extracts (D05/D20/D50/V50/G05), deep-sky
  (30k objects), HyperLeda, GCVS variable stars.

## Building

### Static library + CLI (Linux, Windows, macOS)

```sh
cmake -S . -B build
cmake --build build -j
./build/astap -h
```

Builds the `astap_lib` static library and the `astap` command-line
binary with any C++23 toolchain (GCC, Clang, or MSVC).

### Qt 6 desktop GUI (optional)

```sh
cmake -S . -B build -DASTAP_BUILD_GUI=ON
cmake --build build -j
```

Requires Qt 6. The GUI builds on Linux, Windows, and macOS.

### Apple xcframework (optional)

For embedding the library in iOS, macOS, watchOS, or visionOS apps:

```sh
./build_frameworks
```

Produces `lib/libastap.xcframework` with slices for iOS, iOS Simulator,
visionOS, visionOS Simulator, watchOS, watchOS Simulator, and macOS
(arm64 + x86_64 universal).

### Running the test suite

```sh
./run_tests.sh
```

## Requirements

- C++23 compiler (GCC 13+, Clang 16+, MSVC 19.37+, or Apple Clang 17)
- CMake 3.20+
- [doctest](https://github.com/doctest/doctest) as a sibling directory
  (`../doctest/doctest/doctest.h`) for tests
- Optional: [Qt 6](https://www.qt.io/) for the desktop GUI
- Optional: Xcode Command Line Tools, only for the Apple xcframework build
- Optional: [stb](https://github.com/nothings/stb) as a sibling
  directory (`../stb/stb_image.h`) so the CLI picks up PNG / JPEG / TGA
  decoding

## License

**Mozilla Public License 2.0** — see `LICENSE`. Database credits are in
`ACKNOWLEDGEMENTS`.

## Credits

- Original ASTAP: Han Kleijn, www.hnsky.org
- C++ port: John Stephen / wobbleworks.com
