# ASTAP++

A C++23 port of [ASTAP](https://www.hnsky.org/astap.htm) — the
Astrometric STAcking Program by Han Kleijn. ASTAP solves astronomical
images against star catalogs, stacks exposures, and measures photometry.
ASTAP++ brings the core engine to modern C++ with a headless CLI and a
static library suitable for embedding in iOS, macOS, watchOS, and
visionOS applications.

## Status

Port is functionally complete. The CLI binary builds and runs; a 22-suite
doctest harness covers the pure-math modules (solving, stacking,
photometry, ephemerides, WCS, SIP, annotation, inspector, etc.). The
GUI layer from the original Lazarus program is not part of this port.

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

### CLI binary (native platform)

```sh
cmake -S . -B build
cmake --build build -j
./build/astap -h
```

### Static library + xcframework (Apple platforms)

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

- C++23 compiler (tested with Apple Clang 17)
- CMake 3.20+
- Xcode Command Line Tools (for Apple framework builds)
- [doctest](https://github.com/doctest/doctest) at
  `../third_party/doctest/doctest/doctest.h` (for tests only)

## License

**Mozilla Public License 2.0** — see `LICENSE`. Database credits are in
`ACKNOWLEDGEMENTS`.

## Credits

- Original ASTAP: Han Kleijn, www.hnsky.org
- C++ port: John Stephen / wobbleworks.com
