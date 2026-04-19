##----------------------------------------------------------------------------
##      @file CMakeLists.txt
##   @ingroup ASTAP++/tests
##     @brief Doctest-based unit tests for the ASTAP++ library modules.
##   @details Each *_test.cpp becomes its own executable, registered with
##            CTest.  The subject-under-test .cpp files are compiled into
##            each test binary directly (no library artifact required) so
##            individual tests stay independent of the broader build.
##    @author Created by John Stephen on 4/15/26.
## @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
##----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.20)
project(ASTAPpp_tests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

# ---- Paths ------------------------------------------------------------------

# Expect a sibling `doctest` directory next to this repo:
#   <parent>/astap_cpp/
#   <parent>/doctest/doctest/doctest.h
set(ASTAP_DOCTEST_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/../../doctest")
set(ASTAP_SRC             "${CMAKE_CURRENT_SOURCE_DIR}/../src")

if(NOT EXISTS "${ASTAP_DOCTEST_HEADERS}/doctest/doctest.h")
	message(FATAL_ERROR "doctest header not found at ${ASTAP_DOCTEST_HEADERS}/doctest/doctest.h")
endif()

# ---- Helpers ----------------------------------------------------------------

# astap_add_test(NAME [extra_source1 extra_source2 ...])
#   Builds tests/${NAME}.cpp plus any extra sources into an executable called
#   ${NAME}, adds the doctest headers + ASTAP src/ include path, and
#   registers the result with ctest.
function(astap_add_test NAME)
	add_executable(${NAME} ${NAME}.cpp ${ARGN})
	target_include_directories(${NAME}
		PRIVATE
			"${ASTAP_DOCTEST_HEADERS}"
			"${ASTAP_SRC}"
	)
	target_compile_features(${NAME} PRIVATE cxx_std_23)
	add_test(NAME ${NAME} COMMAND ${NAME})
endfunction()

# ---- Tests ------------------------------------------------------------------

astap_add_test(util_test
	"${ASTAP_SRC}/core/util.cpp"
)

astap_add_test(calc_trans_cubic_test
	"${ASTAP_SRC}/solving/calc_trans_cubic.cpp"
)

astap_add_test(hyperbola_test
	"${ASTAP_SRC}/core/hyperbola.cpp"
)

astap_add_test(aberration_test
	"${ASTAP_SRC}/core/aberration.cpp"
	# ephem::earth_state and ephem::precess_iau1976 are stubbed inside the
	# test to isolate the aberration/nutation math. Link the real
	# ephemerides.cpp only if those stubs are removed.
)

astap_add_test(ephemerides_test
	"${ASTAP_SRC}/core/ephemerides.cpp"
)

astap_add_test(hjd_test
	"${ASTAP_SRC}/core/hjd.cpp"
	"${ASTAP_SRC}/core/util.cpp"         # strtofloat2, fnmodulo
	"${ASTAP_SRC}/core/globals.cpp"      # astap::head, astap::jd_mid, ...
	"${ASTAP_SRC}/core/ephemerides.cpp"  # precess_iau1976
)

astap_add_test(gaussian_blur_test
	"${ASTAP_SRC}/stacking/gaussian_blur.cpp"
)

astap_add_test(raster_rotate_test
	"${ASTAP_SRC}/stacking/raster_rotate.cpp"
)

astap_add_test(online_gaia_test
	"${ASTAP_SRC}/reference/online_gaia.cpp"
)

astap_add_test(demosaic_test
	"${ASTAP_SRC}/core/demosaic.cpp"
)

astap_add_test(wcs_test
	"${ASTAP_SRC}/core/wcs.cpp"
	"${ASTAP_SRC}/core/util.cpp"         # strtofloat2 (used by ra/dec_text_to_radians)
)

astap_add_test(star_align_test
	"${ASTAP_SRC}/solving/star_align.cpp"
)

astap_add_test(star_database_test
	"${ASTAP_SRC}/reference/star_database.cpp"
)

astap_add_test(sqm_test
	"${ASTAP_SRC}/core/sqm.cpp"
	"${ASTAP_SRC}/core/globals.cpp"
	"${ASTAP_SRC}/core/link_stubs.cpp"
	"${ASTAP_SRC}/core/util.cpp"          # for any util helpers calculate_sqm touches
	"${ASTAP_SRC}/core/hjd.cpp"           # airmass_calc + calculate_az_alt
	"${ASTAP_SRC}/core/imaging.cpp"       # histogram arrays
	"${ASTAP_SRC}/core/ephemerides.cpp"  # precess_iau1976 (via hjd.cpp)
)

astap_add_test(photometry_calibration_test
	"${ASTAP_SRC}/core/photometry.cpp"
	"${ASTAP_SRC}/core/wcs.cpp"
	"${ASTAP_SRC}/core/util.cpp"          # get_best_mean, smedian, mad_median
)

astap_add_test(imaging_test
	"${ASTAP_SRC}/core/imaging.cpp"
	"${ASTAP_SRC}/core/globals.cpp"
	"${ASTAP_SRC}/core/link_stubs.cpp"
	"${ASTAP_SRC}/core/util.cpp"
	"${ASTAP_SRC}/stacking/raster_rotate.cpp"   # rotate_arbitrary delegates here
)

astap_add_test(image_sharpness_test
    "${ASTAP_SRC}/analysis/image_sharpness.cpp"
)

astap_add_test(contour_test
    "${ASTAP_SRC}/analysis/contour.cpp"
)

astap_add_test(constellations_test
    "${ASTAP_SRC}/analysis/constellations.cpp"
)

astap_add_test(asteroid_test
    "${ASTAP_SRC}/analysis/asteroid.cpp"
)

astap_add_test(annotation_test
    "${ASTAP_SRC}/analysis/annotation.cpp"
)

astap_add_test(inspector_test
    "${ASTAP_SRC}/analysis/inspector.cpp"
)

# ---- CLI integration test (spawns the astap binary) -------------------------

# Only register the CLI test if the binary has been built. Pass its absolute
# path via a compile definition so the test TU can find it without hardcoding.
set(ASTAP_CLI_BIN "${CMAKE_CURRENT_SOURCE_DIR}/../build/astap")
if(EXISTS "${ASTAP_CLI_BIN}")
	add_executable(cli_test cli_test.cpp)
	target_include_directories(cli_test PRIVATE "${ASTAP_DOCTEST_HEADERS}")
	target_compile_definitions(cli_test PRIVATE
		ASTAP_BIN_PATH="${ASTAP_CLI_BIN}")
	target_compile_features(cli_test PRIVATE cxx_std_23)
	add_test(NAME cli_test COMMAND cli_test)
else()
	message(STATUS "CLI binary not found at ${ASTAP_CLI_BIN} — skipping cli_test. "
	               "Build the top-level project first.")
endif()
