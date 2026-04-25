///----------------------------------------
///      @file aavso_report_test.cpp
///   @ingroup ASTAP++/tests
///     @brief Doctest coverage for @ref astap::core::format_aavso_report.
///    @author Created by John Stephen on 4/24/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "core/aavso_report.h"

using astap::core::AavsoMeasurement;
using astap::core::AavsoOptions;
using astap::core::clean_abbreviation;
using astap::core::format_aavso_report;

TEST_CASE("clean_abbreviation strips sigma tail and underscores") {
	CHECK(clean_abbreviation("000-BCP-306") == "000-BCP-306");
	CHECK(clean_abbreviation("000-BCP-306, σ=0.012") == "000-BCP-306");
	CHECK(clean_abbreviation("SS_Cyg") == "SS Cyg");
	CHECK(clean_abbreviation("SS_Cyg, σ=0.005") == "SS Cyg");
	CHECK(clean_abbreviation("  trim me  ") == "trim me");
	CHECK(clean_abbreviation("") == "");
	CHECK(clean_abbreviation(", σ=0.1") == "");
}

TEST_CASE("single-row Extended-format report assembles correctly") {
	auto m = AavsoMeasurement{};
	m.variable_name   = "SS Cyg";
	m.check_name      = "000-BCP-306";
	m.comp_name       = "ENSEMBLE";
	m.var_magnitude   = 9.123;
	m.var_error       = 0.0123;
	m.check_magnitude = 11.456;
	m.jd              = 2461000.50000;
	m.airmass         = 1.42;
	m.snr             = 100;
	m.filter_band     = "V";

	auto opts = AavsoOptions{};
	opts.observer_code     = "JST01";
	opts.delimiter         = ",";
	opts.ensemble          = true;
	opts.software_version  = "2026.04.24";
	opts.software_settings = "G05, aperture=1.5 HFD";

	const auto report = format_aavso_report(m, opts);

	// Header sanity
	CHECK(report.find("#TYPE=Extended\r\n") != std::string::npos);
	CHECK(report.find("#OBSCODE=JST01\r\n") != std::string::npos);
	CHECK(report.find("#SOFTWARE=ASTAP++, v2026.04.24 (G05, aperture=1.5 HFD)\r\n")
	      != std::string::npos);
	CHECK(report.find("#DELIM=,\r\n") != std::string::npos);
	CHECK(report.find("#DATE=JD\r\n") != std::string::npos);
	CHECK(report.find("#OBSTYPE=CCD\r\n") != std::string::npos);
	// No BAA header lines when baa_style is false
	CHECK(report.find("#LOCATION=") == std::string::npos);
	CHECK(report.find("#TELESCOPE=") == std::string::npos);

	// Row sanity (ensemble mode → CNAME=ENSEMBLE, CMAG=na)
	CHECK(report.find("SS Cyg,2461000.50000,9.123,0.0123,V,NO,STD,"
	                  "ENSEMBLE,na,000-BCP-306,11.456,1.420,na,na,na\r\n")
	      != std::string::npos);
}

TEST_CASE("ensemble=false applies comp/catalog magnitude correction") {
	auto m = AavsoMeasurement{};
	m.variable_name      = "SS Cyg";
	m.check_name         = "000-BCP-142";
	m.comp_name          = "000-BCP-306";
	m.var_magnitude      = 9.000;       // instrumental
	m.check_magnitude    = 11.000;      // instrumental
	m.comp_magnitude     = 10.500;      // instrumental
	m.comp_catalog_mag   = 10.700;      // documented
	m.var_error          = 0.0100;
	m.jd                 = 2461001.0;
	m.airmass            = 1.0;
	m.filter_band        = "V";

	auto opts = AavsoOptions{};
	opts.observer_code = "JST01";
	opts.delimiter     = ",";
	opts.ensemble      = false;

	// correction = 10.700 - 10.500 = +0.200; var_mag = 9.200; check_mag = 11.200
	const auto report = format_aavso_report(m, opts);
	CHECK(report.find("SS Cyg,2461001.00000,9.200,0.0100,V,NO,STD,"
	                  "000-BCP-306,10.700,000-BCP-142,11.200,1.000,na,na,na\r\n")
	      != std::string::npos);
}

TEST_CASE("delta_bv * magnitude_slope is added to variable mag") {
	auto m = AavsoMeasurement{};
	m.variable_name = "X";
	m.var_magnitude = 10.000;
	m.var_error     = 0.01;
	m.jd            = 2461000.0;

	auto opts = AavsoOptions{};
	opts.delimiter      = ",";
	opts.ensemble       = true;
	opts.delta_bv       = 0.5;
	opts.magnitude_slope = 0.04;     // expected delta = 0.5 * 0.04 = +0.020
	const auto report = format_aavso_report(m, opts);
	CHECK(report.find(",10.020,") != std::string::npos);
}

TEST_CASE("HJD date and tab delimiter") {
	auto m = AavsoMeasurement{};
	m.variable_name = "Y";
	m.var_magnitude = 12.0;
	m.jd            = 2461002.5;
	m.snr           = 50;            // → derived var_err = 0.04

	auto opts = AavsoOptions{};
	opts.delimiter = "\t";
	opts.hjd_date  = true;
	opts.ensemble  = true;

	const auto report = format_aavso_report(m, opts);
	CHECK(report.find("#DELIM=tab\r\n") != std::string::npos);
	CHECK(report.find("#DATE=HJD\r\n") != std::string::npos);
	CHECK(report.find("Y\t2461002.50000\t12.000\t0.0400\t") != std::string::npos);
}

TEST_CASE("BAA-style header includes location/telescope/camera") {
	auto m = AavsoMeasurement{};
	m.variable_name = "Z";
	m.var_magnitude = 8.0;
	m.var_error     = 0.01;

	auto opts = AavsoOptions{};
	opts.delimiter  = ",";
	opts.baa_style  = true;
	opts.ensemble   = true;
	opts.site_lat   = "+50.123";
	opts.site_long  = "-1.234";
	opts.site_elev  = "75";
	opts.telescope  = "C8 EdgeHD";
	opts.camera     = "ASI2600MC";

	const auto report = format_aavso_report(m, opts);
	CHECK(report.find("#TYPE=AAVSO EXT BAA V1.00\r\n") != std::string::npos);
	CHECK(report.find("#LOCATION=+50.123 -1.234 75\r\n") != std::string::npos);
	CHECK(report.find("#TELESCOPE=C8 EdgeHD\r\n") != std::string::npos);
	CHECK(report.find("#CAMERA=ASI2600MC\r\n") != std::string::npos);
}

TEST_CASE("multi-frame report keeps shared header and emits one row per measurement") {
	auto rows = std::vector<AavsoMeasurement>{};
	auto base = AavsoMeasurement{};
	base.variable_name = "M";
	base.var_error     = 0.01;
	base.filter_band   = "V";
	for (auto i = 0; i < 3; ++i) {
		base.jd            = 2461000.0 + i;
		base.var_magnitude = 9.0 + i * 0.001;
		rows.push_back(base);
	}

	auto opts = AavsoOptions{};
	opts.delimiter = ",";
	opts.ensemble  = true;

	const auto report = format_aavso_report(rows, opts);

	// One header
	const auto header_pos = report.find("#TYPE=");
	CHECK(header_pos != std::string::npos);
	CHECK(report.find("#TYPE=", header_pos + 1) == std::string::npos);

	// Three rows
	CHECK(report.find("M,2461000.00000,9.000,") != std::string::npos);
	CHECK(report.find("M,2461001.00000,9.001,") != std::string::npos);
	CHECK(report.find("M,2461002.00000,9.002,") != std::string::npos);
}
