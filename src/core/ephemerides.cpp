///----------------------------------------
///      @file ephemerides.cpp
///   @ingroup ASTAP++
///     @brief Ephemeris implementations for ASTAP++.
///   @details Orbit propagation, planet positions and the Earth ephemeris are
///            faithful ports of SLALIB (sla_UE2PV/PV2UE/EL2UE/PLANEL,
///            sla_PLANET Simon-1994, sla_EPV2 VSOP-simplified). IAU 1976
///            precession follows Meeus Ch. 21.2 / 21.3.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. Mozilla Public License 2.0.
///----------------------------------------

#include "ephemerides.h"
#include "ephemerides_epv2_data.h"

#include <algorithm>
#include <cmath>
#include <numbers>

///----------------------------------------
namespace astap::core::ephem {
///----------------------------------------

namespace {

constexpr double pi = std::numbers::pi;

constexpr double kTwoPi    = 2.0 * pi;
constexpr double kDeg2Rad  = pi / 180.0;
constexpr double kArcsec2Rad = pi / (180.0 * 3600.0);

constexpr double kJ2000       = 2451545.0;   ///< Julian Date of J2000.0.
constexpr double kCenturyDays = 36525.0;     ///< Julian days per Julian century.

/// @brief Gaussian gravitational constant (rad/day, used to get GM of the Sun).
constexpr double kGaussK = 0.01720209895;

/// @brief Heliocentric gravitational parameter in (AU^3/day^2).
constexpr double kSolarGM = kGaussK * kGaussK;

/// @brief Normalise an angle into [0, 2π).
[[nodiscard]] double normalise_two_pi(double x) noexcept {
	x = std::fmod(x, kTwoPi);
	return x < 0.0 ? x + kTwoPi : x;
}

/// @brief Normalise an angle into (-π, π].
[[nodiscard]] double normalise_pi(double x) noexcept {
	x = std::fmod(x + pi, kTwoPi);
	if (x < 0.0) x += kTwoPi;
	return x - pi;
}

/// @brief Mean obliquity of the ecliptic at J2000.0 (IAU 2006: 84381.406").
///        Used as the fixed ecliptic→equatorial rotation so that outputs land
///        in the J2000 equatorial frame (not equator-of-date).
constexpr double kJ2000Obliquity = 84381.406 * kArcsec2Rad;

/// @brief Mean obliquity of the ecliptic at the given JD (Meeus 22.2).
///        Returns radians. Kept for callers that specifically need
///        obliquity-of-date; the main ecliptic→equatorial rotations
///        in this module use kJ2000Obliquity instead.
[[nodiscard]] double mean_obliquity(double jd_tt) noexcept {
	const double T = (jd_tt - kJ2000) / kCenturyDays;
	// Arcseconds.
	const double eps_arcsec = 23.0 * 3600.0 + 26.0 * 60.0 + 21.448
	                          - 46.8150 * T - 0.00059 * T * T
	                          + 0.001813 * T * T * T;
	return eps_arcsec * kArcsec2Rad;
}

}  // anonymous namespace

///----------------------------------------
/// MARK: - cartesian_to_spherical / spherical_to_cartesian
///----------------------------------------

SphericalResult cartesian_to_spherical(const Vec3& xyz) noexcept {
	const double x = xyz[0], y = xyz[1], z = xyz[2];
	const double r_xy = std::sqrt(x * x + y * y);

	SphericalResult out;
	out.radius = std::sqrt(x * x + y * y + z * z);
	out.dec    = std::atan2(z, r_xy);
	out.ra     = (r_xy == 0.0) ? 0.0 : normalise_two_pi(std::atan2(y, x));
	return out;
}

Vec3 spherical_to_cartesian(double radius, double dec, double ra) noexcept {
	const double cos_dec = std::cos(dec);
	return {
		radius * cos_dec * std::cos(ra),
		radius * cos_dec * std::sin(ra),
		radius * std::sin(dec),
	};
}

///----------------------------------------
/// MARK: - precess_iau1976
///----------------------------------------

// Follows the rigorous three-angle formulation in Meeus 21.2 / 21.3.
EquatorialCoords precess_iau1976(EquatorialCoords coords, double jd_from, double jd_to) {
	const double T  = (jd_from - kJ2000) / kCenturyDays;
	const double T2 = T * T;

	const double t  = (jd_to - jd_from) / kCenturyDays;
	const double t2 = t * t;
	const double t3 = t2 * t;

	// Meeus 21.2 — rigorous method. Coefficients in arcseconds.
	const double ez_common = (2306.2181 + 1.39656 * T - 0.000139 * T2) * t;

	const double zeta    = (ez_common + (0.30188 - 0.000344 * T) * t2 + 0.017998 * t3) * kArcsec2Rad;
	const double z       = (ez_common + (1.09468 + 0.000066 * T) * t2 + 0.018203 * t3) * kArcsec2Rad;
	const double theta   = ((2004.3109 - 0.85330 * T - 0.000217 * T2) * t
	                        - (0.42665 + 0.000217 * T) * t2 - 0.041833 * t3) * kArcsec2Rad;

	// Apply the three-angle rotation to (ra, dec).  Meeus 21.4.
	const double ra0  = coords.ra;
	const double dec0 = coords.dec;

	const double A = std::cos(dec0) * std::sin(ra0 + zeta);
	const double B = std::cos(theta) * std::cos(dec0) * std::cos(ra0 + zeta)
	               - std::sin(theta) * std::sin(dec0);
	const double C = std::sin(theta) * std::cos(dec0) * std::cos(ra0 + zeta)
	               + std::cos(theta) * std::sin(dec0);

	EquatorialCoords out;
	out.ra  = normalise_two_pi(std::atan2(A, B) + z);
	out.dec = std::asin(std::clamp(C, -1.0, 1.0));
	return out;
}

///----------------------------------------
/// MARK: - SLALIB universal-variable routines
///----------------------------------------

// Faithful transcription of the SLALIB Fortran (via Han Kleijn's Pascal port,
// unit_ephemerides.pas). Pascal 1-based indices [1..13] map to [0..12] here.
namespace sla {

namespace {
constexpr double kGCon = 0.01720209895;      // Gaussian gravitational constant.
constexpr double kCD2S = kGCon / 86400.0;    // Canonical days → seconds.
}  // namespace

int ue2pv(double date, UElements& u, PV& pv) {
	constexpr double kTest = 1.0e-13;
	constexpr int kNitMax  = 25;

	// Unpack the parameters.
	const double CM     = u[0];
	const double ALPHA  = u[1];
	const double T0     = u[2];
	double P0[3] = { u[3], u[4], u[5] };
	double V0[3] = { u[6], u[7], u[8] };
	const double R0     = u[9];
	const double SIGMA0 = u[10];
	const double T      = u[11];
	double PSI          = u[12];

	// Approximately update the universal eccentric anomaly.
	PSI = PSI + (date - T) * kGCon / R0;
	// Time from reference epoch to date (canonical days).
	const double DT = (date - T0) * kGCon;

	int NIT = 1;
	double W = 1.0;
	double TOL = 0.0;
	double S0 = 0.0, S1 = 0.0, S2 = 0.0, S3 = 0.0, R = 0.0;
	double FF = 0.0, FLAST = 0.0, PLAST = 0.0;

	while (std::fabs(W) > TOL) {
		// Form half angles until BETA small enough.
		int N = 0;
		double PSJ  = PSI;
		double PSJ2 = PSJ * PSJ;
		double BETA = ALPHA * PSJ2;
		while (std::fabs(BETA) > 0.7) {
			N += 1;
			BETA /= 4.0;
			PSJ  /= 2.0;
			PSJ2 /= 4.0;
		}
		// Universal variables S0..S3 by nested series.
		S3 = PSJ * PSJ2 * ((((((BETA / 210 + 1)
		                       * BETA / 156 + 1)
		                       * BETA / 110 + 1)
		                       * BETA / 72 + 1)
		                       * BETA / 42 + 1)
		                       * BETA / 20 + 1) / 6;
		S2 = PSJ2 * ((((((BETA / 182 + 1)
		                   * BETA / 132 + 1)
		                   * BETA / 90 + 1)
		                   * BETA / 56 + 1)
		                   * BETA / 30 + 1)
		                   * BETA / 12 + 1) / 2;
		S1 = PSJ + ALPHA * S3;
		S0 = 1 + ALPHA * S2;
		// Undo the angle-halving.
		TOL = kTest;
		while (N > 0) {
			S3 = 2 * (S0 * S3 + PSJ * S2);
			S2 = 2 * S1 * S1;
			S1 = 2 * S0 * S1;
			S0 = 2 * S0 * S0 - 1;
			PSJ = PSJ + PSJ;
			TOL = TOL + TOL;
			N = N - 1;
		}
		// F and F' at the current psi.
		FF = R0 * S1 + SIGMA0 * S2 + CM * S3 - DT;
		R  = R0 * S0 + SIGMA0 * S1 + CM * S2;
		if (NIT == 1) FLAST = FF;
		if (FF * FLAST < 0.0) {
			// Sign change: secant method.
			W = FF * (PLAST - PSI) / (FLAST - FF);
		} else {
			// No sign change: Newton-Raphson.
			if (R == 0.0) return -1;   // null radius vector
			W = FF / R;
		}
		PLAST = PSI;
		FLAST = FF;
		PSI = PSI - W;
		if (NIT > kNitMax) return -2;  // failed to converge
		NIT += 1;
	}

	// Project position/velocity (velocity scaled to AU/s).
	W = CM * S2;
	const double F  = 1 - W / R0;
	const double G  = DT - CM * S3;
	const double FD = -CM * S1 / (R0 * R);
	const double GD = 1 - W / R;
	for (int i = 0; i < 3; ++i) {
		pv[static_cast<std::size_t>(i)]     = P0[i] * F + V0[i] * G;
		pv[static_cast<std::size_t>(i + 3)] = kCD2S * (P0[i] * FD + V0[i] * GD);
	}
	u[11] = date;  // update for fast re-prediction
	u[12] = PSI;
	return 0;
}

int pv2ue(const PV& pv, double date, double pmass, UElements& u) {
	constexpr double kRMin = 1.0e-3;   // minimum distance (AU)
	constexpr double kVMin = 1.0e-3;   // minimum speed (AU per canonical day)

	if (pmass < 0.0) return -1;
	const double CM = 1 + pmass;
	// Unpack, expressing velocity in AU per canonical day.
	const double X  = pv[0], Y = pv[1], Z = pv[2];
	const double XD = pv[3] / kCD2S, YD = pv[4] / kCD2S, ZD = pv[5] / kCD2S;
	const double R  = std::sqrt(X * X + Y * Y + Z * Z);
	const double V2 = XD * XD + YD * YD + ZD * ZD;
	const double V  = std::sqrt(V2);
	if (R < kRMin) return -2;   // too close
	if (V < kVMin) return -3;   // too slow
	const double ALPHA = V2 - 2 * CM / R;      // total energy
	const double RDV   = X * XD + Y * YD + Z * ZD;
	u = { CM, ALPHA, date, X, Y, Z, XD, YD, ZD, R, RDV, date, 0.0 };
	return 0;
}

int el2ue(double date, int jform, double epoch, double orbinc, double anode,
          double perih, double aorq, double e, double aorl, double dm,
          UElements& u) {
	// Sin/cos of the J2000 mean obliquity (IAU 1976).
	constexpr double SE = 0.3977771559319137;
	constexpr double CE = 0.9174820620691818;

	if (jform < 1 || jform > 3) return -1;
	if (e < 0.0 || e > 10.0 || (e >= 1.0 && jform != 3)) return -2;
	if (aorq <= 0.0) return -3;
	if (jform == 1 && dm <= 0.0) return -4;

	// Transform elements into standard form (perihelion epoch/arg/distance/mass).
	double PHT = 0.0, ARGPH = 0.0, Q = 0.0, CM = 0.0, W = 0.0;
	if (jform == 1) {
		PHT   = epoch - (aorl - perih) / dm;
		ARGPH = perih - anode;
		Q     = aorq * (1 - e);
		W     = dm / kGCon;
		CM    = W * W * aorq * aorq * aorq;
	} else if (jform == 2) {
		PHT   = epoch - aorl * std::sqrt(aorq * aorq * aorq) / kGCon;
		ARGPH = perih;
		Q     = aorq * (1 - e);
		CM    = 1;
	} else {
		PHT   = epoch;
		ARGPH = perih;
		Q     = aorq;
		CM    = 1;
	}
	// Universal variable alpha (∝ total energy) and speed at perihelion.
	const double ALPHA = CM * (e - 1) / Q;
	const double PHS   = std::sqrt(ALPHA + 2 * CM / Q);

	// Rotate the perihelion position/velocity into J2000 equatorial frame
	// (Euler: ω about z, i about x, Ω about z, then obliquity ε about x).
	const double SW = std::sin(ARGPH), CW = std::cos(ARGPH);
	const double SI = std::sin(orbinc), CI = std::cos(orbinc);
	const double SO = std::sin(anode),  CO = std::cos(anode);

	// Position at perihelion (AU).
	double X = Q * CW;
	double Y = Q * SW;
	double Z = Y * SI;
	Y = Y * CI;
	const double PX = X * CO - Y * SO;
	Y = X * SO + Y * CO;
	const double PY = Y * CE - Z * SE;
	const double PZ = Y * SE + Z * CE;

	// Velocity at perihelion (AU per canonical day).
	X = -PHS * SW;
	Y = PHS * CW;
	Z = Y * SI;
	Y = Y * CI;
	const double VX = X * CO - Y * SO;
	Y = X * SO + Y * CO;
	const double VY = Y * CE - Z * SE;
	const double VZ = Y * SE + Z * CE;

	// Time from perihelion to date (canonical days).
	const double DT = (date - PHT) * kGCon;
	// First approximation to the universal eccentric anomaly, PSI, blending
	// the circle (FC) and parabola (FP) values.
	const double FC = DT / Q;
	W = std::pow(3 * DT + std::sqrt(9 * DT * DT + 8 * Q * Q * Q), 1.0 / 3.0);
	const double FP = W - 2 * Q / W;
	const double PSI = (1 - e) * FC + e * FP;

	UElements ul = { CM, ALPHA, PHT, PX, PY, PZ, VX, VY, VZ, Q, 0.0, date, PSI };

	// Predict PV at the epoch of osculation, then convert to universal elements.
	PV pv{};
	if (ue2pv(date, ul, pv) != 0) return -5;
	if (pv2ue(pv, date, CM - 1, u) != 0) return -5;
	return 0;
}

int planel(double date, int jform, double epoch, double orbinc, double anode,
           double perih, double aorq, double e, double aorl, double dm, PV& pv) {
	UElements u{};
	int j = el2ue(date, jform, epoch, orbinc, anode, perih, aorq, e, aorl, dm, u);
	if (j == 0) {
		j = ue2pv(date, u, pv);
		if (j != 0) j = -5;
	}
	return j;
}

namespace {

// Simon et al. 1994 mean Keplerian elements (rows = planets 1..8 =
// Mercury..Neptune, 0-based here) and the periodic-perturbation tables, from
// SLALIB sla_PLANET via unit_ephemerides.pas. Generated by a parser, not hand-
// copied. Columns of A/E are AU/dimensionless; DLM/PI/DINC/OMEGA hold
// (degrees, arcsec/millennium, arcsec/millennium^2).
constexpr std::array<double, 8> kAMAS = { 6023600, 408523.5, 328900.5, 3098710, 1047.355, 3498.5, 22869, 19314 };

constexpr std::array<std::array<double, 3>, 8> kA = {{
	{{ 0.3870983098, 0, 0 }},
	{{ 0.7233298200, 0, 0 }},
	{{ 1.0000010178, 0, 0 }},
	{{ 1.5236793419, 3E-10, 0 }},
	{{ 5.2026032092, 19132E-10, -39E-10 }},
	{{ 9.5549091915, -0.0000213896, 444E-10 }},
	{{ 19.2184460618, -3716E-10, 979E-10 }},
	{{ 30.1103868694, -16635E-10, 686E-10 }},
}};

constexpr std::array<std::array<double, 3>, 8> kDLM = {{
	{{ 252.25090552, 5381016286.88982, -1.92789 }},
	{{ 181.97980085, 2106641364.33548, 0.59381 }},
	{{ 100.46645683, 1295977422.83429, -2.04411 }},
	{{ 355.43299958, 689050774.93988, 0.94264 }},
	{{ 34.35151874, 109256603.77991, -30.60378 }},
	{{ 50.07744430, 43996098.55732, 75.61614 }},
	{{ 314.05500511, 15424811.93933, -1.75083 }},
	{{ 304.34866548, 7865503.20744, 0.21103 }},
}};

constexpr std::array<std::array<double, 3>, 8> kE = {{
	{{ 0.2056317526, 0.0002040653, -28349E-10 }},
	{{ 0.0067719164, -0.0004776521, 98127E-10 }},
	{{ 0.0167086342, -0.0004203654, -0.0000126734 }},
	{{ 0.0934006477, 0.0009048438, -80641E-10 }},
	{{ 0.0484979255, 0.0016322542, -0.0000471366 }},
	{{ 0.0555481426, -0.0034664062, -0.0000643639 }},
	{{ 0.0463812221, -0.0002729293, 0.0000078913 }},
	{{ 0.0094557470, 0.0000603263, 0 }},
}};

constexpr std::array<std::array<double, 3>, 8> kPI = {{
	{{ 77.45611904, 5719.11590, -4.83016 }},
	{{ 131.56370300, 175.48640, -498.48184 }},
	{{ 102.93734808, 11612.35290, 53.27577 }},
	{{ 336.06023395, 15980.45908, -62.32800 }},
	{{ 14.33120687, 7758.75163, 259.95938 }},
	{{ 93.05723748, 20395.49439, 190.25952 }},
	{{ 173.00529106, 3215.56238, -34.09288 }},
	{{ 48.12027554, 1050.71912, 27.39717 }},
}};

constexpr std::array<std::array<double, 3>, 8> kDINC = {{
	{{ 7.00498625, -214.25629, 0.28977 }},
	{{ 3.39466189, -30.84437, -11.67836 }},
	{{ 0, 469.97289, -3.35053 }},
	{{ 1.84972648, -293.31722, -8.11830 }},
	{{ 1.30326698, -71.55890, 11.95297 }},
	{{ 2.48887878, 91.85195, -17.66225 }},
	{{ 0.77319689, -60.72723, 1.25759 }},
	{{ 1.76995259, 8.12333, 0.08135 }},
}};

constexpr std::array<std::array<double, 3>, 8> kOMEGA = {{
	{{ 48.33089304, -4515.21727, -31.79892 }},
	{{ 76.67992019, -10008.48154, -51.32614 }},
	{{ 174.87317577, -8679.27034, 15.34191 }},
	{{ 49.55809321, -10620.90088, -230.57416 }},
	{{ 100.46440702, 6362.03561, 326.52178 }},
	{{ 113.66550252, -9240.19942, -66.23743 }},
	{{ 74.00595701, 2669.15033, 145.93964 }},
	{{ 131.78405702, -221.94322, -0.78728 }},
}};

constexpr std::array<std::array<double, 9>, 8> kDKP = {{
	{{ 69613, 75645, 88306, 59899, 15746, 71087, 142173, 3086, 0 }},
	{{ 21863, 32794, 26934, 10931, 26250, 43725, 53867, 28939, 0 }},
	{{ 16002, 21863, 32004, 10931, 14529, 16368, 15318, 32794, 0 }},
	{{ 6345, 7818, 15636, 7077, 8184, 14163, 1107, 4872, 0 }},
	{{ 1760, 1454, 1167, 880, 287, 2640, 19, 2047, 1454 }},
	{{ 574, 0, 880, 287, 19, 1760, 1167, 306, 574 }},
	{{ 204, 0, 177, 1265, 4, 385, 200, 208, 204 }},
	{{ 0, 102, 106, 4, 98, 1367, 487, 204, 0 }},
}};

constexpr std::array<std::array<double, 9>, 8> kCA = {{
	{{ 4, -13, 11, -9, -9, -3, -1, 4, 0 }},
	{{ -156, 59, -42, 6, 19, -20, -10, -12, 0 }},
	{{ 64, -152, 62, -8, 32, -41, 19, -11, 0 }},
	{{ 124, 621, -145, 208, 54, -57, 30, 15, 0 }},
	{{ -23437, -2634, 6601, 6259, -1507, -1821, 2620, -2115, -1489 }},
	{{ 62911, -119919, 79336, 17814, -24241, 12068, 8306, -4893, 8902 }},
	{{ 389061, -262125, -44088, 8387, -22976, -2093, -615, -9720, 6633 }},
	{{ -412235, -157046, -31430, 37817, -9740, -13, -7449, 9644, 0 }},
}};

constexpr std::array<std::array<double, 9>, 8> kSA = {{
	{{ -29, -1, 9, 6, -6, 5, 4, 0, 0 }},
	{{ -48, -125, -26, -37, 18, -13, -20, -2, 0 }},
	{{ -150, -46, 68, 54, 14, 24, -28, 22, 0 }},
	{{ -621, 532, -694, -20, 192, -94, 71, -73, 0 }},
	{{ -14614, -19828, -5869, 1881, -4372, -2255, 782, 930, 913 }},
	{{ 139737, 0, 24667, 51123, -5102, 7429, -4095, -1976, -9566 }},
	{{ -138081, 0, 37205, -49039, -41901, -33872, -27037, -12474, 18797 }},
	{{ 0, 28492, 133236, 69654, 52322, -49577, -26430, -3593, 0 }},
}};

constexpr std::array<std::array<double, 10>, 8> kDKQ = {{
	{{ 3086, 15746, 69613, 59899, 75645, 88306, 12661, 2658, 0, 0 }},
	{{ 21863, 32794, 10931, 73, 4387, 26934, 1473, 2157, 0, 0 }},
	{{ 10, 16002, 21863, 10931, 1473, 32004, 4387, 73, 0, 0 }},
	{{ 10, 6345, 7818, 1107, 15636, 7077, 8184, 532, 10, 0 }},
	{{ 19, 1760, 1454, 287, 1167, 880, 574, 2640, 19, 1454 }},
	{{ 19, 574, 287, 306, 1760, 12, 31, 38, 19, 574 }},
	{{ 4, 204, 177, 8, 31, 200, 1265, 102, 4, 204 }},
	{{ 4, 102, 106, 8, 98, 1367, 487, 204, 4, 102 }},
}};

constexpr std::array<std::array<double, 10>, 8> kCLO = {{
	{{ 21, -95, -157, 41, -5, 42, 23, 30, 0, 0 }},
	{{ -160, -313, -235, 60, -74, -76, -27, 34, 0, 0 }},
	{{ -325, -322, -79, 232, -52, 97, 55, -41, 0, 0 }},
	{{ 2268, -979, 802, 602, -668, -33, 345, 201, -55, 0 }},
	{{ 7610, -4997, -7689, -5841, -2617, 1115, -748, -607, 6074, 354 }},
	{{ -18549, 30125, 20012, -730, 824, 23, 1289, -352, -14767, -2062 }},
	{{ -135245, -14594, 4197, -4030, -5630, -2898, 2540, -306, 2939, 1986 }},
	{{ 89948, 2103, 8963, 2695, 3682, 1648, 866, -154, -1963, -283 }},
}};

constexpr std::array<std::array<double, 10>, 8> kSLO = {{
	{{ -342, 136, -23, 62, 66, -52, -33, 17, 0, 0 }},
	{{ 524, -149, -35, 117, 151, 122, -71, -62, 0, 0 }},
	{{ -105, -137, 258, 35, -116, -88, -112, -80, 0, 0 }},
	{{ 854, -205, -936, -240, 140, -341, -97, -232, 536, 0 }},
	{{ -56980, 8016, 1012, 1448, -3024, -3710, 318, 503, 3767, 577 }},
	{{ 138606, -13478, -4964, 1441, -1319, -1482, 427, 1236, -9167, -1918 }},
	{{ 71234, -41116, 5334, -4935, -1848, 66, 434, -1748, 3780, -701 }},
	{{ -47645, 11647, 2166, 3194, 679, 0, -244, -419, -2531, 48 }},
}};

}  // namespace

int planet(int np, double date, PV& pv) {
	if (np < 1 || np > 8) {
		pv = PV{};
		return -1;
	}
	const int n = np - 1;   // 0-based table row

	// Time: Julian millennia since J2000.
	const double T = (date - kJ2000) / 365250.0;
	const int jstat = (std::abs(T) <= 1.0) ? 0 : 1;   // warn for remote epochs

	// Mean elements (first column of the angle tables is in degrees → arcsec).
	double da        = kA[n][0] + (kA[n][1] + kA[n][2] * T) * T;
	double dl        = (3600.0 * kDLM[n][0] + (kDLM[n][1] + kDLM[n][2] * T) * T) * kArcsec2Rad;
	const double de  = kE[n][0] + (kE[n][1] + kE[n][2] * T) * T;
	const double dpe = normalise_two_pi((3600.0 * kPI[n][0] + (kPI[n][1] + kPI[n][2] * T) * T) * kArcsec2Rad);
	const double di  = (3600.0 * kDINC[n][0] + (kDINC[n][1] + kDINC[n][2] * T) * T) * kArcsec2Rad;
	const double dom = normalise_two_pi((3600.0 * kOMEGA[n][0] + (kOMEGA[n][1] + kOMEGA[n][2] * T) * T) * kArcsec2Rad);

	// Periodic perturbation terms (added to a and mean longitude).
	const double dmu = 0.35953620 * T;
	for (int jj = 0; jj < 8; ++jj) {
		const double arga = kDKP[n][jj] * dmu;
		const double argl = kDKQ[n][jj] * dmu;
		da += (kCA[n][jj]  * std::cos(arga) + kSA[n][jj]  * std::sin(arga)) * 1e-7;
		dl += (kCLO[n][jj] * std::cos(argl) + kSLO[n][jj] * std::sin(argl)) * 1e-7;
	}
	{
		const double arga = kDKP[n][8] * dmu;
		da += T * (kCA[n][8] * std::cos(arga) + kSA[n][8] * std::sin(arga)) * 1e-7;
	}
	for (int jj = 8; jj < 10; ++jj) {
		const double argl = kDKQ[n][jj] * dmu;
		dl += T * (kCLO[n][jj] * std::cos(argl) + kSLO[n][jj] * std::sin(argl)) * 1e-7;
	}
	dl = normalise_two_pi(dl);

	// Daily motion, then predict via the universal-variable machinery. JFORM=1
	// takes (a, mean longitude, daily motion); epoch == date so no propagation.
	const double dm = kGaussK * std::sqrt((1.0 + 1.0 / kAMAS[n]) / (da * da * da));
	const int j = planel(date, 1, date, di, dom, dpe, da, de, dl, dm, pv);
	if (j < 0) {
		return -2;
	}
	return jstat;
}

namespace {

// SLALIB EPV2 constants.
constexpr double kEpvDjy = 365.25;   // Julian year in days (SLALIB DJY)

// Ecliptic → ICRS rotation (SLALIB, empirically fit to DE405 over 1900-2100).
// AM11 = 1, AM31 = 0 implicitly.
constexpr double kAM12 = +0.000000211284;
constexpr double kAM13 = -0.000000091603;
constexpr double kAM21 = -0.000000230286;
constexpr double kAM22 = +0.917482137087;
constexpr double kAM23 = -0.397776982902;
constexpr double kAM32 = +0.397776982902;
constexpr double kAM33 = +0.917482137087;

// Accumulate one Cartesian axis over a table's T^0/T^1/T^2 series, matching
// unit_ephemerides.pas sla_EPV2. @p t is Julian years since J2000.
template <std::size_t M0, std::size_t M1, std::size_t M2>
void epv2_axis(int k, double t, double t2,
               const std::array<std::array<std::array<double, 3>, M0>, 3>& e0, int n0,
               const std::array<std::array<std::array<double, 3>, M1>, 3>& e1, int n1,
               const std::array<std::array<std::array<double, 3>, M2>, 3>& e2, int n2,
               double& xyz, double& xyzd) {
	const auto kk = static_cast<std::size_t>(k);
	for (int j = 0; j < n0; ++j) {
		const auto& term = e0[kk][static_cast<std::size_t>(j)];
		const double p = term[1] + term[2] * t;
		xyz  += term[0] * std::cos(p);
		xyzd -= term[0] * term[2] * std::sin(p);
	}
	for (int j = 0; j < n1; ++j) {
		const auto& term = e1[kk][static_cast<std::size_t>(j)];
		const double ct = term[2] * t;
		const double cp = std::cos(term[1] + ct);
		xyz  += term[0] * t * cp;
		xyzd += term[0] * (cp - ct * std::sin(term[1] + ct));
	}
	for (int j = 0; j < n2; ++j) {
		const auto& term = e2[kk][static_cast<std::size_t>(j)];
		const double ct = term[2] * t;
		const double cp = std::cos(term[1] + ct);
		xyz  += term[0] * t2 * cp;
		xyzd += term[0] * t * (2.0 * cp - ct * std::sin(term[1] + ct));
	}
}

}  // namespace

void epv2(double date, bool bary, Vec3& pe, Vec3& ve) {
	using namespace epv2_data;

	const double t  = (date - kJ2000) / kEpvDjy;   // Julian years since J2000
	const double t2 = t * t;

	Vec3 hp{}, hv{}, bp{}, bv{};
	for (int k = 0; k < 3; ++k) {
		const auto kk = static_cast<std::size_t>(k);
		double xyz  = 0.0;
		double xyzd = 0.0;

		// Sun-to-Earth vector.
		epv2_axis(k, t, t2, kE0, kNE0[kk], kE1, kNE1[kk], kE2, kNE2[kk], xyz, xyzd);
		hp[kk] = xyz;
		hv[kk] = xyzd / kEpvDjy;   // AU/day

		if (bary) {
			// SSB-to-Sun vector, accumulated onto the Sun-to-Earth vector.
			epv2_axis(k, t, t2, kS0, kNS0[kk], kS1, kNS1[kk], kS2, kNS2[kk], xyz, xyzd);
			bp[kk] = xyz;
			bv[kk] = xyzd / kEpvDjy;
		}
	}

	const Vec3& p = bary ? bp : hp;
	const Vec3& v = bary ? bv : hv;

	// Rotate ecliptic → ICRS (J2000 equatorial).
	pe = { p[0] + kAM12 * p[1] + kAM13 * p[2],
	       kAM21 * p[0] + kAM22 * p[1] + kAM23 * p[2],
	       kAM32 * p[1] + kAM33 * p[2] };
	ve = { v[0] + kAM12 * v[1] + kAM13 * v[2],
	       kAM21 * v[0] + kAM22 * v[1] + kAM23 * v[2],
	       kAM32 * v[1] + kAM33 * v[2] };
}

}  // namespace sla

///----------------------------------------
/// MARK: - propagate
///----------------------------------------

State propagate(const OrbitalElements& elements, double jd_tt) {
	const double e = elements.eccentricity;
	const double q = elements.perihelion_distance;

	// Choose the SLALIB element form. Minor-planet (JFORM=2) needs a bound
	// orbit with a real mean anomaly; everything else — comets, parabolae and
	// hyperbolae, and the "at perihelion at epoch" case (M==0) — uses the comet
	// form (JFORM=3), whose epoch is the perihelion passage. This mirrors the
	// asteroid-vs-comet split the Pascal caller made.
	int jform;
	double aorq, aorl;
	if (e < 1.0 && elements.mean_anomaly_at_epoch != 0.0) {
		jform = 2;
		aorq  = q / (1.0 - e);                       // semimajor axis a
		aorl  = elements.mean_anomaly_at_epoch;      // mean anomaly M
	} else {
		jform = 3;
		aorq  = q;                                   // perihelion distance q
		aorl  = 0.0;
	}

	// DATE and EPOCH share the JD offset, which cancels in every difference
	// SLALIB forms, so full JD is passed directly.
	sla::PV pv{};
	const int jstat = sla::planel(jd_tt, jform, elements.epoch_jd,
		elements.inclination, elements.ascending_node,
		elements.argument_of_perihelion, aorq, e, aorl, /*dm=*/0.0, pv);

	State out;
	if (jstat != 0) return out;   // failed: return a zero state

	// SLALIB returns J2000 equatorial position (AU) and velocity (AU/s);
	// the port's contract is AU/day.
	for (int k = 0; k < 3; ++k) {
		out.position[static_cast<std::size_t>(k)] = pv[static_cast<std::size_t>(k)];
		out.velocity[static_cast<std::size_t>(k)] = pv[static_cast<std::size_t>(k + 3)] * 86400.0;
	}
	return out;
}

///----------------------------------------
/// MARK: - planet_state
///----------------------------------------

State planet_state(int np, double jd_tt) {
	sla::PV pv{};
	const int jstat = sla::planet(np, jd_tt, pv);
	State out;
	if (jstat < 0) {
		return out;   // illegal planet number or numerical error → zero state
	}
	// SLALIB PV: J2000 equatorial, AU and AU/s → port contract AU and AU/day.
	for (int k = 0; k < 3; ++k) {
		out.position[static_cast<std::size_t>(k)] = pv[static_cast<std::size_t>(k)];
		out.velocity[static_cast<std::size_t>(k)] = pv[static_cast<std::size_t>(k + 3)] * 86400.0;
	}
	return out;
}

///----------------------------------------
/// MARK: - earth_state (SLALIB sla_EPV2, VSOP-simplified)
///----------------------------------------

State earth_state(double jd_tt, ReferenceFrame frame) {
	// SLALIB sla_EPV2 gives the J2000-equatorial Earth state directly for both
	// frames: heliocentric (Sun→Earth) and barycentric (SSB→Earth, via its
	// SSB→Sun terms). No separate Sun-barycentre offset is needed.
	State out;
	sla::epv2(jd_tt, frame == ReferenceFrame::Barycentric,
	          out.position, out.velocity);
	return out;
}

///----------------------------------------
/// MARK: - vsop::evaluate
///----------------------------------------

namespace vsop {

HelioEcliptic evaluate(const Body& body, double jd_tt) noexcept {
	// Julian millennia from J2000.
	const double tau = (jd_tt - kJ2000) / (kCenturyDays * 10.0);

	auto eval_component = [tau](const Component& comp) -> double {
		double result = 0.0;
		double tau_power = 1.0;
		for (std::size_t i = 0; i < comp.count; ++i) {
			const Series& s = comp.series[i];
			double sum = 0.0;
			for (std::size_t k = 0; k < s.count; ++k) {
				const Term& t = s.terms[k];
				sum += t.amplitude * std::cos(t.phase + t.frequency * tau);
			}
			result += sum * tau_power;
			tau_power *= tau;
		}
		return result;
	};

	HelioEcliptic out;
	out.longitude = normalise_two_pi(eval_component(body.L));
	out.latitude  = eval_component(body.B);
	out.radius    = eval_component(body.R);
	return out;
}

}  // namespace vsop

}  // namespace astap::core::ephem
