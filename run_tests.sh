#!/usr/bin/env bash
#
# Build and run the ASTAP++ test suites.
#
# Usage:
#   ./run_tests.sh          # configure + build + run all tests
#   ./run_tests.sh --build  # build only (skip configure if already done)
#   ./run_tests.sh --run    # run only (skip build)
#
# Exit code is 0 when all tests pass, non-zero otherwise.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TESTS_SRC="${SCRIPT_DIR}/tests"
TESTS_BUILD="${TESTS_SRC}/build"
CLI_BUILD="${SCRIPT_DIR}/build"

configure() {
    echo "==> Configuring tests..."
    cmake -S "${TESTS_SRC}" -B "${TESTS_BUILD}"
}

build_cli() {
    echo "==> Building CLI binary..."
    cmake -S "${SCRIPT_DIR}" -B "${CLI_BUILD}"
    cmake --build "${CLI_BUILD}" -j
}

build_tests() {
    echo "==> Building test suites..."
    cmake --build "${TESTS_BUILD}" -j
}

run() {
    echo "==> Running tests..."
    ctest --test-dir "${TESTS_BUILD}" --output-on-failure
}

case "${1:-all}" in
    --build)
        build_cli
        configure
        build_tests
        ;;
    --run)
        run
        ;;
    *)
        build_cli
        configure
        build_tests
        run
        ;;
esac
