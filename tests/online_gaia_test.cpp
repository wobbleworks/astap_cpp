///----------------------------------------
///     @file online_gaia_test.cpp
///   @ingroup ASTAP++/tests
///    @brief Unit tests for the Gaia photometric transformations
///           (src/reference/online_gaia.cpp).
///  @details Verifies the pure-math `transform_gaia` polynomial for every
///           supported passband, pass-through behaviour for "BP", out-of-
///           range rejection, and `convert_magnitudes` cache side effect.
///   @author Created by John Stephen on 4/15/26.
///@copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "reference/online_gaia.h"

#include <cmath>
#include <string>

using namespace astap::reference;

///----------------------------------------
/// MARK: transform_gaia
///----------------------------------------

TEST_CASE("transform_gaia pass-through for BP") {
	// For "BP" the function should return the Gaia BP magnitude unchanged.
	CHECK(transform_gaia("BP", 12.5, 13.0, 12.1) == doctest::Approx(13.0));
	CHECK(transform_gaia("BP", 0.0,   7.5,  0.0) == doctest::Approx(7.5));
}

TEST_CASE("transform_gaia Johnson V is near G magnitude for solar-type stars") {
	// For a G2V star (BP - RP ≈ 0.82), V-G correction is small (~0.05 mag).
	// Our reference point is G=10.0, BP-RP=0.82 → V ≈ G + small offset.
	const double V = transform_gaia("V", /*G=*/10.0, /*BP=*/10.41, /*RP=*/9.59);

	CHECK(V > 0.0);                         // produced a valid result
	CHECK(std::abs(V - 10.0) < 0.2);        // close to G by less than 0.2 mag
}

TEST_CASE("transform_gaia returns zero for colours outside the valid range") {
	// The Pascal polynomials are validity-gated around |BP-RP| <= ~5. A
	// wildly blue or red star should fall outside and return 0.
	CHECK(transform_gaia("V", 10.0, 0.0,  20.0) == doctest::Approx(0.0));
	CHECK(transform_gaia("V", 10.0, 20.0,  0.0) == doctest::Approx(0.0));
}

TEST_CASE("transform_gaia Johnson B > V for a red-ish star") {
	// Red star: BP-RP = 1.5 (K-type). B should be fainter than V.
	const double B = transform_gaia("B", 12.0, 12.6, 11.1);
	const double V = transform_gaia("V", 12.0, 12.6, 11.1);

	// Both defined and B > V (larger magnitude = fainter at red end).
	CHECK(B > 0.0);
	CHECK(V > 0.0);
	CHECK(B > V);
}

TEST_CASE("transform_gaia supports SDSS passbands SG, SR, SI") {
	const double sg = transform_gaia("SG", 11.0, 11.3, 10.6);
	const double sr = transform_gaia("SR", 11.0, 11.3, 10.6);
	const double si = transform_gaia("SI", 11.0, 11.3, 10.6);

	CHECK(sg > 0.0);
	CHECK(sr > 0.0);
	CHECK(si > 0.0);
	// SDSS passbands nest from blue (g) through red (i); for an average
	// star colour, sg > sr > si in magnitude (redder = fainter at short
	// wavelengths).
	CHECK(sg > sr);
	CHECK(sr > si);
}

TEST_CASE("transform_gaia unknown filter returns zero") {
	CHECK(transform_gaia("NOPE", 10.0, 10.5, 9.5) == doctest::Approx(0.0));
	CHECK(transform_gaia("",     10.0, 10.5, 9.5) == doctest::Approx(0.0));
}

///----------------------------------------
/// MARK: convert_magnitudes / passband_active
///----------------------------------------

TEST_CASE("convert_magnitudes is a no-op when the passband is already active") {
	// Seed the cache: pretend online_database was last converted to V.
	passband_active = "V";
	// Also seed some dummy data in the online_database so we can verify
	// it's unchanged by a same-band request.
	online_database[5] = {1.0, 2.0, 3.0};

	convert_magnitudes("V");

	CHECK(passband_active == "V");
	CHECK(online_database[5] == std::vector<double>{1.0, 2.0, 3.0});
}

TEST_CASE("convert_magnitudes updates passband_active on a real switch") {
	// The Pascal function only writes the derived magnitudes when the
	// call is actually changing the cached passband. Set up a fake
	// transformation and confirm the cache tag moves.
	passband_active = "";
	online_database[0] = {0.0};   // RA
	online_database[1] = {0.0};   // Dec
	online_database[2] = {10.0};  // G
	online_database[3] = {10.3};  // BP
	online_database[4] = {9.7};   // RP
	online_database[5] = {0.0};   // active

	convert_magnitudes("V");
	CHECK(passband_active == "V");

	convert_magnitudes("BP");
	CHECK(passband_active == "BP");

	// Subsequent same-band request is a no-op — passband_active stays.
	convert_magnitudes("BP");
	CHECK(passband_active == "BP");
}
