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
# Default: sibling `doctest` directory alongside the repo root.
# CI (and anyone with doctest elsewhere) can override via
# `-DASTAP_DOCTEST_HEADERS=/path/to/doctest`.
if(NOT DEFINED ASTAP_DOCTEST_HEADERS OR ASTAP_DOCTEST_HEADERS STREQUAL "")
	set(ASTAP_DOCTEST_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/../../doctest")
endif()
set(ASTAP_SRC             "${CMAKE_CURRENT_SOURCE_DIR}/../src")

if(NOT EXISTS "${ASTAP_DOCTEST_HEADERS}/doctest/doctest.h")
	message(FATAL_ERROR "doctest header not found at ${ASTAP_DOCTEST_HEADERS}/doctest/doctest.h "
		"(override with -DASTAP_DOCTEST_HEADERS=<path>)")
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

astap_add_test(live_stacking_test
	"${ASTAP_SRC}/stacking/live_colour_correction.cpp"
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
	"${ASTAP_SRC}/core/util.cpp"          # smedian (find_fit median)
)

astap_add_test(star_database_test
	"${ASTAP_SRC}/reference/star_database.cpp"
	"memo2_test_stub.cpp"
)

astap_add_test(sqm_test
	"${ASTAP_SRC}/core/sqm.cpp"
	"${ASTAP_SRC}/core/globals.cpp"
	"${ASTAP_SRC}/core/link_stubs.cpp"
	"${ASTAP_SRC}/core/util.cpp"          # for any util helpers calculate_sqm touches
	"${ASTAP_SRC}/core/hjd.cpp"           # airmass_calc + calculate_az_alt
	"${ASTAP_SRC}/core/imaging.cpp"       # histogram arrays
	"${ASTAP_SRC}/core/ephemerides.cpp"  # precess_iau1976 (via hjd.cpp)
	"memo2_test_stub.cpp"                 # no-op memo2_message
)

astap_add_test(photometry_calibration_test
	"${ASTAP_SRC}/core/photometry.cpp"
	"${ASTAP_SRC}/core/wcs.cpp"
	"${ASTAP_SRC}/core/util.cpp"          # get_best_mean, smedian, mad_median
	"memo2_test_stub.cpp"                 # no-op memo2_message
)

astap_add_test(imaging_test
	"${ASTAP_SRC}/core/imaging.cpp"
	"${ASTAP_SRC}/core/globals.cpp"
	"${ASTAP_SRC}/core/link_stubs.cpp"
	"${ASTAP_SRC}/core/util.cpp"
	"${ASTAP_SRC}/stacking/raster_rotate.cpp"   # rotate_arbitrary delegates here
	"memo2_test_stub.cpp"                        # no-op memo2_message
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
    "${ASTAP_SRC}/core/wcs.cpp"           # plot_deepsky uses pixel/celestial projection
    "${ASTAP_SRC}/core/util.cpp"          # strtofloat2 (via wcs ra/dec text parsers)
)

astap_add_test(inspector_test
    "${ASTAP_SRC}/analysis/inspector.cpp"
)

astap_add_test(aavso_report_test
    "${ASTAP_SRC}/core/aavso_report.cpp"
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

# ---- Differential test (ASTAP++ vs original astap_cli oracle) ---------------

# Compares the ported binary against the reference oracle over a shared image
# corpus. Registered only when the port binary exists; at runtime it skips
# individual cases (passing) when the oracle or a corpus image is missing, so
# CI without the private corpus stays green.
#
#   ASTAP_ORACLE_BIN  reference astap_cli   (default: tests/oracle/astap_cli)
#   ASTAP_CORPUS_DIR  shared image corpus   (default: ~/Desktop/ASTAP Tests)
if(NOT DEFINED ASTAP_ORACLE_BIN OR ASTAP_ORACLE_BIN STREQUAL "")
	set(ASTAP_ORACLE_BIN "${CMAKE_CURRENT_SOURCE_DIR}/oracle/astap_cli")
endif()
if(NOT DEFINED ASTAP_CORPUS_DIR OR ASTAP_CORPUS_DIR STREQUAL "")
	set(ASTAP_CORPUS_DIR "$ENV{HOME}/Desktop/ASTAP Tests")
endif()

if(EXISTS "${ASTAP_CLI_BIN}")
	add_executable(differential_analyse_test differential_analyse_test.cpp)
	target_include_directories(differential_analyse_test PRIVATE "${ASTAP_DOCTEST_HEADERS}")
	target_compile_definitions(differential_analyse_test PRIVATE
		ASTAP_BIN_PATH="${ASTAP_CLI_BIN}"
		ASTAP_ORACLE_BIN="${ASTAP_ORACLE_BIN}"
		ASTAP_CORPUS_DIR="${ASTAP_CORPUS_DIR}")
	target_compile_features(differential_analyse_test PRIVATE cxx_std_23)
	add_test(NAME differential_analyse_test COMMAND differential_analyse_test)
	message(STATUS "differential_analyse_test: oracle=${ASTAP_ORACLE_BIN} corpus=${ASTAP_CORPUS_DIR}")
else()
	message(STATUS "CLI binary not found — skipping differential_analyse_test.")
endif()

# ---- Differential solve test (needs the shared star database) ---------------

#   ASTAP_DB_DIR  shared star database (default: /usr/local/opt/astap)
if(NOT DEFINED ASTAP_DB_DIR OR ASTAP_DB_DIR STREQUAL "")
	set(ASTAP_DB_DIR "/usr/local/opt/astap")
endif()

if(EXISTS "${ASTAP_CLI_BIN}")
	add_executable(differential_solve_test differential_solve_test.cpp)
	target_include_directories(differential_solve_test PRIVATE "${ASTAP_DOCTEST_HEADERS}")
	target_compile_definitions(differential_solve_test PRIVATE
		ASTAP_BIN_PATH="${ASTAP_CLI_BIN}"
		ASTAP_ORACLE_BIN="${ASTAP_ORACLE_BIN}"
		ASTAP_CORPUS_DIR="${ASTAP_CORPUS_DIR}"
		ASTAP_DB_DIR="${ASTAP_DB_DIR}")
	target_compile_features(differential_solve_test PRIVATE cxx_std_23)
	add_test(NAME differential_solve_test COMMAND differential_solve_test)
	message(STATUS "differential_solve_test: db=${ASTAP_DB_DIR}")
else()
	message(STATUS "CLI binary not found — skipping differential_solve_test.")
endif()

# ---- Batch-stacker self-consistency test (spawns the astap binary) ----------

if(EXISTS "${ASTAP_CLI_BIN}")
	add_executable(stack_test stack_test.cpp)
	target_include_directories(stack_test PRIVATE "${ASTAP_DOCTEST_HEADERS}")
	target_compile_definitions(stack_test PRIVATE
		ASTAP_BIN_PATH="${ASTAP_CLI_BIN}"
		ASTAP_CORPUS_DIR="${ASTAP_CORPUS_DIR}"
		ASTAP_DB_DIR="${ASTAP_DB_DIR}")
	target_compile_features(stack_test PRIVATE cxx_std_23)
	add_test(NAME stack_test COMMAND stack_test)
	message(STATUS "stack_test: db=${ASTAP_DB_DIR} corpus=${ASTAP_CORPUS_DIR}")
else()
	message(STATUS "CLI binary not found — skipping stack_test.")
endif()

# ---- Artificial-stars / distortion self-consistency test --------------------

# Solves a corpus frame with the local star database, then exercises the two
# headless annotation features through the CLI (-artificialstars, -distortion).
# No astap_cli oracle exists (both are GUI-only / online-only in the original),
# so this is a self-consistency test. Skips when the binary, database, or corpus
# images are absent.
if(EXISTS "${ASTAP_CLI_BIN}")
	add_executable(annotation_features_test annotation_features_test.cpp)
	target_include_directories(annotation_features_test PRIVATE "${ASTAP_DOCTEST_HEADERS}")
	target_compile_definitions(annotation_features_test PRIVATE
		ASTAP_BIN_PATH="${ASTAP_CLI_BIN}"
		ASTAP_CORPUS_DIR="${ASTAP_CORPUS_DIR}"
		ASTAP_DB_DIR="${ASTAP_DB_DIR}")
	target_compile_features(annotation_features_test PRIVATE cxx_std_23)
	add_test(NAME annotation_features_test COMMAND annotation_features_test)
	message(STATUS "annotation_features_test: db=${ASTAP_DB_DIR} corpus=${ASTAP_CORPUS_DIR}")
else()
	message(STATUS "CLI binary not found — skipping annotation_features_test.")
endif()
