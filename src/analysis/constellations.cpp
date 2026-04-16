///----------------------------------------
///      @file constellations.cpp
///   @ingroup analysis
///     @brief Constellation stick-figure line data, label positions, and IAU abbreviations.
///   @details Ported from unit_constellations.pas (version 2001-10-21).
///            Original copyright (C) 1997, 2021 by Han Kleijn, www.hnsky.org.
///            Licensed under the Mozilla Public License, v. 2.0.
///    @author Created by John Stephen on 4/15/26.
/// @copyright Copyright 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "constellations.h"

///----------------------------------------
namespace astap::analysis {
///----------------------------------------

///----------------------------------------
/// MARK: Constellation stick-figure vertices
///----------------------------------------

const std::array<ConstStar, kConstellationLength> constellation = {{
	// Andromeda (0)
	{-2,   140,  2909, "\xCE\xB1"},  // Alpha And
	{-1,   655,  3086, "\xCE\xB4"},  // Delta And
	{-1,  1162,  3562, "\xCE\xB2"},  // Beta And
	{-1,  2065,  4233, "\xCE\xB3"},  // Gamma 1 And
	{-2,  1162,  3562, "\xCE\xB2"},  // Beta And
	{-1,   946,  3850, "\xCE\xBC"},  // Mu And
	{-1,   830,  4108, "\xCE\xBD"},  // Nu And

	// Antlia (1)
	{-2, 10945, -3714, "\xCE\xB9"},  // Iota Ant
	{-1, 10453, -3107, "\xCE\xB1"},  // Alpha Ant
	{-1,  9487, -3595, "\xCE\xB5"},  // Epsilon Ant

	// Apus (2)
	{-2, 14798, -7904, "\xCE\xB1"},  // Alpha Aps
	{-1, 16558, -7890, "\xCE\xB3"},  // Gamma Aps
	{-1, 16718, -7752, "\xCE\xB2"},  // Beta Aps

	// Aquarius (3)
	{-2, 22877,  -758, "\xCE\xBB"},  // Lambda Aqr
	{-1, 22589,   -12, "\xCE\xB7"},  // Eta Aqr
	{-1, 22480,    -2, "\xCE\xB6"},  // Zeta 1 Aqr
	{-1, 22361,  -139, "\xCE\xB3"},  // Gamma Aqr
	{-1, 22096,   -32, "\xCE\xB1"},  // Alpha Aqr
	{-1, 21526,  -557, "\xCE\xB2"},  // Beta Aqr
	{-1, 20795,  -950, "\xCE\xB5"},  // Epsilon Aqr
	{-2, 22361,  -139, "\xCE\xB3"},  // Gamma Aqr
	{-1, 22281,  -778, "\xCE\xB8"},  // Theta Aqr
	{-1, 22911, -1582, "\xCE\xB4"},  // Delta Aqr

	// Aquila (4)
	{-2, 20189,   -82, "\xCE\xB8"},  // Theta Aql
	{-1, 19922,   641, "\xCE\xB2"},  // Beta Aql
	{-1, 19846,   887, "\xCE\xB1"},  // Alpha Aql
	{-1, 19771,  1061, "\xCE\xB3"},  // Gamma Aql
	{-1, 19090,  1386, "\xCE\xB6"},  // Zeta Aql
	{-1, 18994,  1507, "\xCE\xB5"},  // Epsilon Aql
	{-2, 19846,   887, "\xCE\xB1"},  // Alpha Aql
	{-1, 19425,   311, "\xCE\xB4"},  // Delta Aql
	{-1, 19104,  -488, "\xCE\xBB"},  // Lambda Aql

	// Ara (5)
	{-2, 18110, -5009, "\xCE\xB8"},  // Theta Ara
	{-1, 17531, -4988, "\xCE\xB1"},  // Alpha Ara
	{-1, 17422, -5553, "\xCE\xB2"},  // Beta Ara
	{-1, 17423, -5638, "\xCE\xB3"},  // Gamma Ara
	{-1, 17518, -6068, "\xCE\xB4"},  // Delta Ara
	{-2, 17422, -5553, "\xCE\xB2"},  // Beta Ara
	{-1, 16977, -5599, "\xCE\xB6"},  // Zeta Ara
	{-1, 16830, -5904, "\xCE\xB7"},  // Eta Ara

	// Aries (6)
	{-2,  2833,  2726, "41"},        // 41 Ari
	{-1,  2120,  2346, "\xCE\xB1"},  // Alpha Ari
	{-1,  1911,  2081, "\xCE\xB2"},  // Beta Ari
	{-1,  1892,  1930, "\xCE\xB3"},  // Gamma 1 Ari B

	// Auriga (7)
	{-2,  5278,  4600, "\xCE\xB1"},  // Alpha Aur
	{-1,  5033,  4382, "\xCE\xB5"},  // Epsilon Aur
	{-1,  4950,  3317, "\xCE\xB9"},  // Iota Aur
	{-1,  5438,  2861, "\xCE\xB2"},  // Beta Tau (shared)
	{-1,  5995,  3721, "\xCE\xB8"},  // Theta Aur
	{-1,  5992,  4495, "\xCE\xB2"},  // Beta Aur
	{-1,  5278,  4600, "\xCE\xB1"},  // Alpha Aur

	// Bootes (8)
	{-2, 13911,  1840, "\xCE\xB7"},  // Eta Boo
	{-1, 14261,  1918, "\xCE\xB1"},  // Alpha Boo
	{-1, 14531,  3037, "\xCF\x81"},  // Rho Boo
	{-1, 14535,  3831, "\xCE\xB3"},  // Gamma Boo
	{-1, 15032,  4039, "\xCE\xB2"},  // Beta Boo
	{-1, 15258,  3331, "\xCE\xB4"},  // Delta Boo
	{-1, 14750,  2707, "\xCE\xB5"},  // Epsilon Boo
	{-1, 14261,  1918, "\xCE\xB1"},  // Alpha Boo
	{-1, 14686,  1373, "\xCE\xB6"},  // Zeta Boo

	// Caelum (9)
	{-2,  4514, -4495, "\xCE\xB4"},  // Delta Cae
	{-1,  4676, -4186, "1"},         // 1 Cae
	{-1,  4701, -3714, "\xCE\xB2"},  // Beta Cae
	{-1,  5073, -3548, "\xCE\xB3"},  // Gamma Cae

	// Camelopardalis (10)
	{-2, 12821,  8341, ""},          // SAO2102 Cam
	{-1,  7001,  7698, ""},          // 6022 Cam
	{-1,  6314,  6932, ""},          // SAO13788 Cam
	{-1,  4901,  6634, "\xCE\xB1"},  // Alpha Cam
	{-1,  5057,  6044, "\xCE\xB2"},  // Beta Cam
	{-1,  4955,  5375, "7"},         // 7 Cam
	{-2,  3484,  5994, "CS"},        // CS Cam
	{-1,  3825,  6553, "BE"},        // BE Cam
	{-1,  3839,  7133, "\xCE\xB3"},  // Gamma Cam
	{-1,  6314,  6932, ""},          // SAO13788 Cam

	// Cancer (11)
	{-2,  8975,  1186, "\xCE\xB1"},  // Alpha Cnc
	{-1,  8745,  1815, "\xCE\xB4"},  // Delta Cnc
	{-1,  8721,  2147, "\xCE\xB3"},  // Gamma Cnc
	{-1,  8778,  2876, "\xCE\xB9"},  // Iota Cnc
	{-2,  8275,   919, "\xCE\xB2"},  // Beta Cnc
	{-1,  8745,  1815, "\xCE\xB4"},  // Delta Cnc
	{-1,  8204,  1765, "\xCE\xB6"},  // Zeta Cnc

	// Canes Venatici (12)
	{-2, 12933,  3831, "\xCE\xB1"},  // Alpha 1 CVn
	{-1, 12562,  4136, "\xCE\xB2"},  // Beta CVn

	// Canis Major (13)
	{-2,  7402, -2930, "\xCE\xB7"},  // Eta CMa
	{-1,  7140, -2639, "\xCE\xB4"},  // Delta CMa
	{-1,  6752, -1672, "\xCE\xB1"},  // Alpha CMa
	{-1,  6378, -1796, "\xCE\xB2"},  // Beta CMa
	{-2,  7063, -1563, "\xCE\xB3"},  // Gamma CMa
	{-1,  6936, -1705, "\xCE\xB9"},  // Iota CMa
	{-1,  6752, -1672, "\xCE\xB1"},  // Alpha CMa
	{-2,  6977, -2897, "\xCE\xB5"},  // Epsilon CMa
	{-1,  7140, -2639, "\xCE\xB4"},  // Delta CMa

	// Canis Minor (14)
	{-2,  7655,   523, "\xCE\xB1"},  // Alpha CMi
	{-1,  7453,   829, "\xCE\xB2"},  // Beta CMi
	{-1,  7469,   893, "\xCE\xB3"},  // Gamma CMi

	// Capricornus (15)
	{-2, 20294, -1251, "\xCE\xB1"},  // Alpha 1 Cap
	{-1, 20301, -1255, "\xCE\xB1" "2"}, // Alpha2 Cap
	{-1, 20350, -1478, "\xCE\xB2"},  // Beta 1 Cap
	{-1, 20768, -2527, "\xCF\x88"},  // Psi Cap
	{-1, 20864, -2692, "\xCF\x89"},  // Omega Cap
	{-1, 21444, -2241, "\xCE\xB6"},  // Zeta Cap
	{-1, 21784, -1613, "\xCE\xB4"},  // Delta Cap
	{-1, 21668, -1666, "\xCE\xB3"},  // Gamma Cap
	{-1, 21371, -1683, "\xCE\xB9"},  // Iota Cap
	{-1, 21099, -1723, "\xCE\xB8"},  // Theta Cap
	{-1, 20350, -1478, "\xCE\xB2"},  // Beta 1 Cap

	// Carina (16)
	{-2,  6399, -5270, "\xCE\xB1"},  // Alpha Car
	{-1,  7946, -5298, "\xCF\x87"},  // Chi Car
	{-1,  8375, -5951, "\xCE\xB5"},  // Epsilon Car
	{-1,  9183, -5897, "a"},         // a Car
	{-1,  9285, -5928, "\xCE\xB9"},  // Iota Car
	{-1, 10285, -6133, "q"},         // q Car
	{-1, 10716, -6439, "\xCE\xB8"},  // Theta Car
	{-1, 10229, -7004, "\xCF\x89"},  // Omega Car
	{-1,  9220, -6972, "\xCE\xB2"},  // Beta Car
	{-1,  9785, -6507, "\xCF\x85"},  // Upsilon Car
	{-1,  9285, -5928, "\xCE\xB9"},  // Iota Car

	// Cassiopeia (17)
	{-2,  1907,  6367, "\xCE\xB5"},  // Epsilon Cas
	{-1,  1430,  6024, "\xCE\xB4"},  // Delta Cas
	{-1,   945,  6072, "\xCE\xB3"},  // Gamma Cas
	{-1,   675,  5654, "\xCE\xB1"},  // Alpha Cas
	{-1,   153,  5915, "\xCE\xB2"},  // Beta Cas

	// Centaurus (18)
	{-2, 14660, -6084, "\xCE\xB1"},  // Alpha 1 Cen
	{-1, 14064, -6037, "\xCE\xB2"},  // Beta Cen
	{-1, 13665, -5347, "\xCE\xB5"},  // Epsilon Cen
	{-1, 13926, -4729, "\xCE\xB6"},  // Zeta Cen
	{-1, 14592, -4216, "\xCE\xB7"},  // Eta Cen
	{-1, 14111, -3637, "\xCE\xB8"},  // Theta Cen
	{-1, 13343, -3671, "\xCE\xB9"},  // Iota Cen
	{-1, 12692, -4896, "\xCE\xB3"},  // Gamma Cen
	{-2, 12139, -5072, "\xCE\xB4"},  // Delta Cen
	{-1, 12692, -4896, "\xCE\xB3"},  // Gamma Cen
	{-1, 13665, -5347, "\xCE\xB5"},  // Epsilon Cen

	// Cepheus (19)
	{-2, 21478,  7056, "\xCE\xB2"},  // Beta Cep
	{-1, 21310,  6259, "\xCE\xB1"},  // Alpha Cep
	{-1, 22486,  5842, "\xCE\xB4"},  // Delta Cep
	{-1, 22828,  6620, "\xCE\xB9"},  // Iota Cep
	{-1, 23656,  7763, "\xCE\xB3"},  // Gamma Cep
	{-1, 21478,  7056, "\xCE\xB2"},  // Beta Cep
	{-1, 22828,  6620, "\xCE\xB9"},  // Iota Cep

	// Cetus (20)
	{-2,  3038,   409, "\xCE\xB1"},  // Alpha Cet
	{-1,  2722,   324, "\xCE\xB3"},  // Gamma Cet
	{-1,  2658,    33, "\xCE\xB4"},  // Delta Cet
	{-1,  2322,  -298, "\xCE\xBF"},  // Omicron Cet
	{-1,  1858, -1034, "\xCE\xB6"},  // Zeta Cet
	{-1,  1734, -1594, "\xCF\x84"},  // Tau Cet
	{-1,   726, -1799, "\xCE\xB2"},  // Beta Cet
	{-1,   324,  -882, "\xCE\xB9"},  // Iota Cet
	{-2,  1858, -1034, "\xCE\xB6"},  // Zeta Cet
	{-1,  1400,  -818, "\xCE\xB8"},  // Theta Cet
	{-1,  1143, -1018, "\xCE\xB7"},  // Eta Cet
	{-1,   726, -1799, "\xCE\xB2"},  // Beta Cet

	// Chamaeleon (21)
	{-2,  8309, -7692, "\xCE\xB1"},  // Alpha Cha
	{-1, 10591, -7861, "\xCE\xB3"},  // Gamma Cha
	{-1, 12306, -7931, "\xCE\xB2"},  // Beta Cha
	{-1, 10763, -8054, "\xCE\xB4" "2"}, // Delta2 Cha
	{-1,  8344, -7748, "\xCE\xB8"},  // Theta Cha
	{-1,  8309, -7692, "\xCE\xB1"},  // Alpha Cha

	// Circinus (22)
	{-2, 15292, -5880, "\xCE\xB2"},  // Beta Cir
	{-1, 14708, -6498, "\xCE\xB1"},  // Alpha Cir
	{-1, 15390, -5932, "\xCE\xB3"},  // Gamma Cir

	// Columba (23)
	{-2,  6369, -3344, "\xCE\xB4"},  // Delta Col
	{-1,  5849, -3577, "\xCE\xB2"},  // Beta Col
	{-1,  5661, -3407, "\xCE\xB1"},  // Alpha Col
	{-1,  5520, -3547, "\xCE\xB5"},  // Epsilon Col
	{-1,  5849, -3577, "\xCE\xB2"},  // Beta Col

	// Coma Berenices (24)
	{-2, 13166,  1753, "\xCE\xB1"},  // Alpha Com
	{-1, 13198,  2788, "\xCE\xB2"},  // Beta Com
	{-1, 12449,  2827, "\xCE\xB3"},  // Gamma Com

	// Corona Australis (25)
	{-2, 19107, -3706, "\xCE\xB3"},  // Gamma CrA
	{-1, 19158, -3790, "\xCE\xB1"},  // Alpha CrA
	{-1, 19167, -3934, "\xCE\xB2"},  // Beta CrA

	// Corona Borealis (26)
	{-2, 15960,  2688, "\xCE\xB5"},  // Epsilon CrB
	{-1, 15712,  2630, "\xCE\xB3"},  // Gamma CrB
	{-1, 15578,  2671, "\xCE\xB1"},  // Alpha CrB
	{-1, 15464,  2911, "\xCE\xB2"},  // Beta CrB
	{-1, 15549,  3136, "\xCE\xB8"},  // Theta CrB

	// Corvus (27)
	{-2, 12140, -2473, "\xCE\xB1"},  // Alpha Crv
	{-1, 12169, -2262, "\xCE\xB5"},  // Epsilon Crv
	{-1, 12263, -1754, "\xCE\xB3"},  // Gamma Crv
	{-1, 12498, -1652, "\xCE\xB4"},  // Delta Crv
	{-1, 12573, -2340, "\xCE\xB2"},  // Beta Crv
	{-1, 12169, -2262, "\xCE\xB5"},  // Epsilon Crv

	// Crater (28)
	{-2, 10996, -1830, "\xCE\xB1"},  // Alpha Crt
	{-1, 11194, -2283, "\xCE\xB2"},  // Beta Crt
	{-1, 11415, -1768, "\xCE\xB3"},  // Gamma Crt
	{-1, 11322, -1478, "\xCE\xB4"},  // Delta Crt
	{-1, 10996, -1830, "\xCE\xB1"},  // Alpha Crt

	// Crux (29)
	{-2, 12443, -6310, "\xCE\xB1"},  // Alpha 1 Cru
	{-1, 12519, -5711, "\xCE\xB3"},  // Gamma Cru a
	{-2, 12795, -5969, "\xCE\xB2"},  // Beta Crux
	{-1, 12252, -5875, "\xCE\xB4"},  // Delta Cru

	// Cygnus (30)
	{-2, 20691,  4528, "\xCE\xB1"},  // Alpha Cyg
	{-1, 20370,  4026, "\xCE\xB3"},  // Gamma Cyg
	{-1, 19938,  3508, "\xCE\xB7"},  // Eta Cyg
	{-1, 19843,  3291, "\xCF\x87"},  // Chi Cyg
	{-1, 19513,  2797, "\xCE\xB2"},  // Beta Cyg
	{-2, 21216,  3023, "\xCE\xB6"},  // Zeta Cyg
	{-1, 20770,  3397, "\xCE\xB5"},  // Epsilon Cyg
	{-1, 20370,  4026, "\xCE\xB3"},  // Gamma Cyg
	{-1, 19750,  4513, "\xCE\xB4"},  // Delta Cyg
	{-1, 19495,  5173, "\xCE\xB9"},  // Iota Cyg
	{-1, 19285,  5337, "\xCE\xBA"},  // Kappa Cyg

	// Delphinus (31)
	{-2, 20554,  1130, "\xCE\xB5"},  // Epsilon Del
	{-1, 20588,  1467, "\xCE\xB6"},  // Zeta Del
	{-1, 20661,  1591, "\xCE\xB1"},  // Alpha Del
	{-1, 20777,  1612, "\xCE\xB3"},  // Gamma 1 Del
	{-1, 20724,  1507, "\xCE\xB4"},  // Delta Del
	{-1, 20626,  1460, "\xCE\xB2"},  // Beta Del
	{-1, 20588,  1467, "\xCE\xB6"},  // Zeta Del

	// Dorado (32)
	{-2,  4267, -5149, "\xCE\xB3"},  // Gamma Dor
	{-1,  4567, -5505, "\xCE\xB1"},  // Alpha Dor
	{-1,  5560, -6249, "\xCE\xB2"},  // Beta Dor
	{-1,  5746, -6574, "\xCE\xB4"},  // Delta Dor

	// Draco (33)
	{-2, 12558,  6979, "\xCE\xBA"},  // Kappa Dra
	{-1, 14073,  6438, "\xCE\xB1"},  // Alpha Dra
	{-1, 15415,  5897, "\xCE\xB9"},  // Iota Dra
	{-1, 16031,  5857, "\xCE\xB8"},  // Theta Dra
	{-1, 16400,  6151, "\xCE\xB7"},  // Eta Dra
	{-1, 17146,  6571, "\xCE\xB6"},  // Zeta Dra
	{-1, 19803,  7027, "\xCE\xB5"},  // Epsilon Dra
	{-1, 19209,  6766, "\xCE\xB4"},  // Delta Dra
	{-1, 17536,  5518, "\xCE\xBD"},  // Nu 1 Dra
	{-1, 17507,  5230, "\xCE\xB2"},  // Beta Dra
	{-1, 17943,  5149, "\xCE\xB3"},  // Gamma Dra
	{-1, 17536,  5518, "\xCE\xBD"},  // Nu 1 Dra

	// Equuleus (34)
	{-2, 21264,   525, "\xCE\xB1"},  // Alpha Equ
	{-1, 21241,  1001, "\xCE\xB4"},  // Delta Equ
	{-1, 21172,  1013, "\xCE\xB3"},  // Gamma Equ
	{-1, 20985,   429, "\xCE\xB5"},  // Epsilon Equ
	{-1, 21264,   525, "\xCE\xB1"},  // Alpha Equ

	// Eridanus (35)
	{-2,  5131,  -509, "\xCE\xB2"},  // Beta Eri
	{-1,  4605,  -335, "\xCE\xBD"},  // Nu Eri
	{-1,  3967, -1351, "\xCE\xB3"},  // Gamma Eri
	{-1,  3721,  -976, "\xCE\xB4"},  // Delta Eri
	{-1,  3549,  -946, "\xCE\xB5"},  // Epsilon Eri
	{-1,  2940,  -890, "\xCE\xB7"},  // Eta Eri
	{-1,  4298, -3380, "41"},        // 41 Eri
	{-1,  2971, -4030, "\xCE\xB8"},  // Theta 1 Eri
	{-1,  2275, -5150, "\xCF\x86"},  // Phi Eri
	{-1,  1933, -5161, "\xCF\x87"},  // Chi Eri
	{-1,  1629, -5724, "\xCE\xB1"},  // Alpha Eri

	// Fornax (36)
	{-2,  3704, -3194, "\xCE\xB4"},  // Delta For
	{-1,  3201, -2899, "\xCE\xB1"},  // Alpha For
	{-1,  2818, -3241, "\xCE\xB2"},  // Beta For
	{-1,  2075, -2930, "\xCE\xBD"},  // Nu For

	// Gemini (37)
	{-2,  6629,  1640, "\xCE\xB3"},  // Gamma Gem
	{-1,  7068,  2057, "\xCE\xB6"},  // Zeta Gem
	{-1,  7335,  2198, "\xCE\xB4"},  // Delta Gem
	{-1,  7755,  2803, "\xCE\xB2"},  // Beta Gem
	{-1,  7577,  3189, "\xCE\xB1"},  // Alpha Gem
	{-1,  6732,  2513, "\xCE\xB5"},  // Epsilon Gem
	{-1,  6383,  2251, "\xCE\xBC"},  // Mu Gem
	{-1,  6248,  2251, "\xCE\xB7"},  // Eta Gem

	// Grus (38)
	{-2, 21899, -3737, "\xCE\xB3"},  // Gamma Gru
	{-1, 22488, -4350, "\xCE\xB4"},  // Delta 1 Gru
	{-1, 22496, -4375, "\xCE\xB4" "2"}, // Delta2 Gru
	{-1, 22711, -4688, "\xCE\xB2"},  // Beta Gru
	{-1, 22809, -5132, "\xCE\xB5"},  // Epsilon Gru
	{-1, 23015, -5275, "\xCE\xB6"},  // Zeta Gru
	{-2, 22137, -4696, "\xCE\xB1"},  // Alpha Gru
	{-1, 22488, -4350, "\xCE\xB4"},  // Delta 1 Gru

	// Hercules (39)
	{-2, 17938,  3725, "\xCE\xB8"},  // Theta Her
	{-1, 17251,  3681, "\xCF\x80"},  // Pi Her
	{-1, 17005,  3093, "\xCE\xB5"},  // Epsilon Her
	{-1, 17251,  2484, "\xCE\xB4"},  // Delta Her
	{-1, 17244,  1439, "\xCE\xB1"},  // Alpha 1 Her
	{-2, 17251,  3681, "\xCF\x80"},  // Pi Her
	{-1, 16715,  3892, "\xCE\xB7"},  // Eta Her
	{-1, 16688,  3160, "\xCE\xB6"},  // Zeta Her
	{-1, 16504,  2149, "\xCE\xB2"},  // Beta Her
	{-1, 16365,  1915, "\xCE\xB3"},  // Gamma Her
	{-2, 17005,  3093, "\xCE\xB5"},  // Epsilon Her
	{-1, 16688,  3160, "\xCE\xB6"},  // Zeta Her
	{-2, 16715,  3892, "\xCE\xB7"},  // Eta Her
	{-1, 16568,  4244, "\xCF\x83"},  // Sigma Her
	{-1, 16329,  4631, "\xCF\x84"},  // Tau Her

	// Horologium (40)
	{-2,  4233, -4229, "\xCE\xB1"},  // Alpha Hor
	{-1,  2623, -5254, "\xCE\xB7"},  // Eta Hor
	{-1,  2980, -6407, "\xCE\xB2"},  // Beta Hor

	// Hydra (41)
	{-2, 14106, -2668, "\xCF\x80"},  // Pi Hya
	{-1, 13495, -2328, "R"},         // R Hya
	{-1, 13315, -2317, "\xCE\xB3"},  // Gamma Hya
	{-1, 11882, -3391, "\xCE\xB2"},  // Beta Hya
	{-1, 11550, -3186, "\xCE\xBE"},  // Xi Hya
	{-1, 10827, -1619, "\xCE\xBD"},  // Nu Hya
	{-1, 10176, -1235, "\xCE\xBB"},  // Lambda Hya
	{-1,  9858, -1485, "\xCF\x85"},  // Upsilon 1 Hya
	{-1,  9460,  -866, "\xCE\xB1"},  // Alpha Hya
	{-1,  9664,  -114, "\xCE\xB9"},  // Iota Hya
	{-1,  9239,   231, "\xCE\xB8"},  // Theta Hya
	{-1,  8923,   595, "\xCE\xB6"},  // Zeta Hya
	{-1,  8780,   642, "\xCE\xB5"},  // Epsilon Hya
	{-1,  8628,   570, "\xCE\xB4"},  // Delta Hya
	{-1,  8720,   340, "\xCE\xB7"},  // Eta Hya
	{-1,  8923,   595, "\xCE\xB6"},  // Zeta Hya

	// Hydrus (42)
	{-2,  1980, -6157, "\xCE\xB1"},  // Alpha Hyi
	{-1,   429, -7725, "\xCE\xB2"},  // Beta Hyi
	{-1,  3787, -7424, "\xCE\xB3"},  // Gamma Hyi
	{-1,  1980, -6157, "\xCE\xB1"},  // Alpha Hyi

	// Indus (43)
	{-2, 20626, -4729, "\xCE\xB1"},  // Alpha Ind
	{-1, 21331, -5345, "\xCE\xB8"},  // Theta Ind
	{-1, 20913, -5845, "\xCE\xB2"},  // Beta Ind
	{-2, 21331, -5345, "\xCE\xB8"},  // Theta Ind
	{-1, 21965, -5499, "\xCE\xB4"},  // Delta Ind

	// Lacerta (44)
	{-2, 22393,  5223, "\xCE\xB2"},  // Beta Lac
	{-1, 22522,  5028, "\xCE\xB1"},  // Alpha Lac
	{-1, 22266,  3775, "1"},         // 1 Lac

	// Leo (45)
	{-2,  9764,  2377, "\xCE\xB5"},  // Epsilon Leo
	{-1, 10278,  2342, "\xCE\xB6"},  // Zeta Leo
	{-1, 10333,  1984, "\xCE\xB3"},  // Gamma Leo
	{-1, 10122,  1676, "\xCE\xB7"},  // Eta Leo
	{-1, 10140,  1197, "\xCE\xB1"},  // Alpha Leo
	{-1, 11237,  1543, "\xCE\xB8"},  // Theta Leo
	{-1, 11818,  1457, "\xCE\xB2"},  // Beta Leo
	{-1, 11235,  2052, "\xCE\xB4"},  // Delta Leo
	{-1, 10333,  1984, "\xCE\xB3"},  // Gamma Leo

	// Leo Minor (46)
	{-2, 10889,  3421, "46"},        // 46 LMi
	{-1, 10465,  3671, "\xCE\xB2"},  // Beta LMi
	{-1, 10124,  3524, "21"},        // 21 LMi
	{-1,  9570,  3640, "10"},        // 10 LMi

	// Lepus (47)
	{-2,  5940, -1417, "\xCE\xB7"},  // Eta Lep
	{-1,  5855, -2088, "\xCE\xB4"},  // Delta Lep
	{-1,  5741, -2245, "\xCE\xB3"},  // Gamma Lep
	{-1,  5471, -2076, "\xCE\xB2"},  // Beta Lep
	{-1,  5091, -2237, "\xCE\xB5"},  // Epsilon Lep
	{-1,  5216, -1621, "\xCE\xBC"},  // Mu Lep
	{-1,  5545, -1782, "\xCE\xB1"},  // Alpha Lep
	{-1,  5855, -2088, "\xCE\xB4"},  // Delta Lep
	{-1,  5783, -1482, "\xCE\xB6"},  // Zeta Lep
	{-2,  5545, -1782, "\xCE\xB1"},  // Alpha Lep
	{-1,  5471, -2076, "\xCE\xB2"},  // Beta Lep

	// Libra (48)
	{-2, 15592, -1479, "\xCE\xB3"},  // Gamma Lib
	{-1, 15283,  -938, "\xCE\xB2"},  // Beta Lib
	{-1, 14848, -1604, "\xCE\xB1" "2"}, // Alpha2 Lib

	// Lupus (49)
	{-2, 14699, -4739, "\xCE\xB1"},  // Alpha Lup
	{-1, 14976, -4313, "\xCE\xB2"},  // Beta Lup
	{-1, 15356, -4065, "\xCE\xB4"},  // Delta Lup
	{-1, 16002, -3840, "\xCE\xB7"},  // Eta Lup
	{-1, 15586, -4117, "\xCE\xB3"},  // Gamma Lup
	{-1, 15378, -4469, "\xCE\xB5"},  // Epsilon Lup
	{-1, 15199, -4874, "\xCE\xBA"},  // Kappa 1 Lup
	{-1, 15205, -5210, "\xCE\xB6"},  // Zeta Lup
	{-1, 14699, -4739, "\xCE\xB1"},  // Alpha Lup

	// Lynx (50)
	{-2,  9351,  3439, "\xCE\xB1"},  // Alpha Lyn
	{-1,  9314,  3680, "38"},        // 38 Lyn
	{-1,  9011,  4178, "10"},        // 10 Lyn
	{-1,  8381,  4319, "31"},        // 31 Lyn
	{-1,  7445,  4921, "21"},        // 21 Lyn

	// Lyra (51)
	{-2, 18616,  3878, "\xCE\xB1"},  // Alpha Lyr
	{-1, 18746,  3761, "\xCE\xB6"},  // Zeta 1 Lyr
	{-1, 18835,  3336, "\xCE\xB2"},  // Beta Lyr
	{-1, 18982,  3269, "\xCE\xB3"},  // Gamma Lyr
	{-1, 18908,  3690, "\xCE\xB4" "2"}, // Delta2 Lyr
	{-1, 18746,  3761, "\xCE\xB6"},  // Zeta 1 Lyr

	// Mensa (52)
	{-2,  6171, -7475, "\xCE\xB1"},  // Alpha Men
	{-1,  5531, -7634, "\xCE\xB3"},  // Gamma Men
	{-1,  4920, -7494, "\xCE\xB7"},  // Eta Men
	{-1,  5045, -7131, "\xCE\xB2"},  // Beta Men

	// Microscopium (53)
	{-2, 20833, -3378, "\xCE\xB1"},  // Alpha Mic
	{-1, 21022, -3226, "\xCE\xB3"},  // Gamma Mic
	{-1, 21299, -3217, "\xCE\xB5"},  // Epsilon Mic

	// Monoceros (54)
	{-2,  8143,  -298, "\xCE\xB6"},  // Zeta Mon
	{-1,  7687,  -955, "\xCE\xB1"},  // Alpha Mon
	{-1,  7198,   -49, "\xCE\xB4"},  // Delta Mon
	{-2,  6396,   459, "\xCE\xB5"},  // Epsilon Mon
	{-1,  7198,   -49, "\xCE\xB4"},  // Delta Mon
	{-1,  6480,  -703, "\xCE\xB2"},  // Beta Mon
	{-1,  6248,  -627, "\xCE\xB3"},  // Gamma Mon

	// Musca (55)
	{-2, 11760, -6673, "\xCE\xBB"},  // Lambda Mus
	{-1, 12293, -6796, "\xCE\xB5"},  // Epsilon Mus
	{-1, 12620, -6914, "\xCE\xB1"},  // Alpha Mus
	{-1, 12771, -6811, "\xCE\xB2"},  // Beta Mus
	{-1, 13038, -7155, "\xCE\xB4"},  // Delta Mus
	{-1, 12541, -7213, "\xCE\xB3"},  // Gamma Mus
	{-1, 12620, -6914, "\xCE\xB1"},  // Alpha Mus

	// Norma (56)
	{-2, 16331, -5016, "\xCE\xB3" "2"}, // Gamma2 Nor
	{-1, 16054, -4923, "\xCE\xB7"},  // Eta Nor

	// Octans (57)
	{-2, 22768, -8138, "\xCE\xB2"},  // Beta Oct
	{-1, 14449, -8367, "\xCE\xB4"},  // Delta Oct
	{-1, 21691, -7739, "\xCE\xBD"},  // Nu Oct
	{-1, 22768, -8138, "\xCE\xB2"},  // Beta Oct

	// Ophiuchus (58)
	{-2, 17367, -2500, "\xCE\xB8"},  // Theta Oph
	{-1, 17173, -1572, "\xCE\xB7"},  // Eta Oph
	{-1, 17725,   457, "\xCE\xB2"},  // Beta Oph
	{-1, 17582,  1256, "\xCE\xB1"},  // Alpha Oph
	{-1, 16961,   938, "\xCE\xBA"},  // Kappa Oph
	{-1, 16239,  -369, "\xCE\xB4"},  // Delta Oph
	{-1, 16305,  -469, "\xCE\xB5"},  // Epsilon Oph
	{-1, 16619, -1057, "\xCE\xB6"},  // Zeta Oph
	{-1, 17173, -1572, "\xCE\xB7"},  // Eta Oph

	// Orion (59)
	{-2,  5679,  -194, "\xCE\xB6"},  // Zeta Ori
	{-1,  5920,   741, "\xCE\xB1"},  // Alpha Ori
	{-1,  5586,   993, "\xCE\xBB"},  // Lambda Ori
	{-1,  5419,   635, "\xCE\xB3"},  // Gamma Ori
	{-1,  5533,   -28, "\xCE\xB4"},  // Delta Ori
	{-1,  5604,  -120, "\xCE\xB5"},  // Epsilon Ori
	{-1,  5679,  -194, "\xCE\xB6"},  // Zeta Ori
	{-1,  5796,  -967, "\xCE\xBA"},  // Kappa Ori
	{-1,  5242,  -820, "\xCE\xB2"},  // Beta Ori
	{-1,  5533,   -28, "\xCE\xB4"},  // Delta Ori

	// Pavo (60)
	{-2, 21441, -6537, "\xCF\x84"},  // Tau Pav
	{-1, 20749, -6620, "\xCE\xB2"},  // Beta Pav
	{-1, 20010, -7291, "\xCE\xB5"},  // Epsilon Pav
	{-1, 18717, -7143, "\xCE\xB6"},  // Zeta Pav
	{-1, 17762, -6472, "\xCE\xB7"},  // Eta Pav
	{-1, 18387, -6149, "\xCE\xBE"},  // Xi Pav
	{-1, 20145, -6618, "\xCE\xB4"},  // Delta Pav
	{-1, 20749, -6620, "\xCE\xB2"},  // Beta Pav

	// Pegasus (61)
	{-2, 22717,  3022, "\xCE\xB7"},  // Eta Peg
	{-1, 23063,  2808, "\xCE\xB2"},  // Beta Peg
	{-1,   140,  2909, "\xCE\xB1"},  // Alpha And
	{-1,   221,  1518, "\xCE\xB3"},  // Gamma Peg
	{-1, 23079,  1521, "\xCE\xB1"},  // Alpha Peg
	{-1, 22691,  1083, "\xCE\xB6"},  // Zeta Peg
	{-1, 22170,   620, "\xCE\xB8"},  // Theta Peg
	{-1, 21736,   988, "\xCE\xB5"},  // Epsilon Peg
	{-2, 23063,  2808, "\xCE\xB2"},  // Beta Peg
	{-1, 23079,  1521, "\xCE\xB1"},  // Alpha Peg

	// Perseus (62)
	{-2,  2845,  5590, "\xCE\xB7"},  // Eta Per
	{-1,  3080,  5351, "\xCE\xB3"},  // Gamma Per
	{-1,  3405,  4986, "\xCE\xB1"},  // Alpha Per
	{-1,  3715,  4779, "\xCE\xB4"},  // Delta Per
	{-1,  3964,  4001, "\xCE\xB5"},  // Epsilon Per
	{-1,  3902,  3188, "\xCE\xB6"},  // Zeta Per
	{-2,  3405,  4986, "\xCE\xB1"},  // Alpha Per
	{-1,  3136,  4096, "\xCE\xB2"},  // Beta Per
	{-1,  3086,  3884, "\xCF\x81"},  // Rho Per

	// Phoenix (63)
	{-2,   157, -4575, "\xCE\xB5"},  // Epsilon Phe
	{-1,   437, -4368, "\xCE\xBA"},  // Kappa Phe
	{-1,  1101, -4672, "\xCE\xB2"},  // Beta Phe
	{-1,  1473, -4332, "\xCE\xB3"},  // Gamma Phe
	{-1,  1521, -4907, "\xCE\xB4"},  // Delta Phe
	{-1,  1140, -5525, "\xCE\xB6"},  // Zeta Phe
	{-1,  1101, -4672, "\xCE\xB2"},  // Beta Phe

	// Pictor (64)
	{-2,  6803, -6194, "\xCE\xB1"},  // Alpha Pic
	{-1,  5830, -5617, "\xCE\xB3"},  // Gamma Pic
	{-1,  5788, -5107, "\xCE\xB2"},  // Beta Pic

	// Pisces (65)
	{-2,  1525,  1535, "\xCE\xB7"},  // Eta Psc
	{-1,  1757,   916, "\xCE\xBF"},  // Omicron Psc
	{-1,  2034,   276, "\xCE\xB1"},  // Alpha Psc
	{-1,  1049,   789, "\xCE\xB5"},  // Epsilon Psc
	{-1,   811,   759, "\xCE\xB4"},  // Delta Psc
	{-1, 23989,   686, "\xCF\x89"},  // Omega Psc
	{-1, 23666,   563, "\xCE\xB9"},  // Iota Psc
	{-1, 23466,   638, "\xCE\xB8"},  // Theta Psc
	{-1, 23065,   382, "\xCE\xB2"},  // Beta Psc
	{-1, 23286,   328, "\xCE\xB3"},  // Gamma Psc
	{-1, 23666,   563, "\xCE\xB9"},  // Iota Psc

	// Piscis Austrinus (66)
	{-2, 21749, -3303, "\xCE\xB9"},  // Iota PsA
	{-1, 22525, -3235, "\xCE\xB2"},  // Beta PsA
	{-1, 22875, -3288, "\xCE\xB3"},  // Gamma PsA
	{-1, 22932, -3254, "\xCE\xB4"},  // Delta PsA
	{-1, 22961, -2962, "\xCE\xB1"},  // Alpha PsA
	{-1, 22678, -2704, "\xCE\xB5"},  // Epsilon PsA
	{-1, 22525, -3235, "\xCE\xB2"},  // Beta PsA

	// Puppis (67)
	{-2,  8126, -2430, "\xCF\x81"},  // Rho Pup
	{-1,  7822, -2486, "\xCE\xBE"},  // Xi Pup
	{-1,  8060, -4000, "\xCE\xB6"},  // Zeta Pup
	{-1,  7286, -3710, "\xCF\x80"},  // Pi Pup
	{-1,  7487, -4330, "\xCF\x83"},  // Sigma Pup
	{-1,  7226, -4464, "L2"},        // L2 Pup
	{-1,  6629, -4320, "\xCE\xBD"},  // Nu Pup
	{-1,  6832, -5061, "\xCF\x84"},  // Tau Pup
	{-1,  7226, -4464, "L2"},        // L2 Pup
	{-2,  8060, -4000, "\xCE\xB6"},  // Zeta Pup
	{-1,  7487, -4330, "\xCF\x83"},  // Sigma Pup

	// Pyxis (68)
	{-2,  8668, -3531, "\xCE\xB2"},  // Beta Pyx
	{-1,  8727, -3319, "\xCE\xB1"},  // Alpha Pyx
	{-1,  8842, -2771, "\xCE\xB3"},  // Gamma Pyx

	// Reticulum (69)
	{-2,  4240, -6247, "\xCE\xB1"},  // Alpha Ret
	{-1,  3737, -6481, "\xCE\xB2"},  // Beta Ret
	{-1,  4015, -6216, "\xCE\xB3"},  // Gamma Ret
	{-1,  3979, -6140, "4"},         // 4 Ret
	{-1,  4275, -5930, "\xCE\xB5"},  // Epsilon Ret
	{-1,  4240, -6247, "\xCE\xB1"},  // Alpha Ret

	// Sagitta (70)
	{-2, 19668,  1801, "\xCE\xB1"},  // Alpha Sge
	{-1, 19790,  1853, "\xCE\xB4"},  // Delta Sge
	{-1, 19684,  1748, "\xCE\xB2"},  // Beta Sge
	{-2, 19790,  1853, "\xCE\xB4"},  // Delta Sge
	{-1, 19979,  1949, "\xCE\xB3"},  // Gamma Sge

	// Sagittarius (71)
	{-2, 19387, -4480, "\xCE\xB2" "2"}, // Beta2 Sgr
	{-1, 19377, -4446, "\xCE\xB2"},  // Beta 1 Sgr
	{-1, 19398, -4062, "\xCE\xB1"},  // Alpha Sgr
	{-1, 19044, -2988, "\xCE\xB6"},  // Zeta Sgr
	{-1, 18921, -2630, "\xCF\x83"},  // Sigma Sgr
	{-1, 18466, -2542, "\xCE\xBB"},  // Lambda Sgr
	{-1, 18350, -2983, "\xCE\xB4"},  // Delta Sgr
	{-1, 18403, -3438, "\xCE\xB5"},  // Epsilon Sgr
	{-1, 18294, -3676, "\xCE\xB7"},  // Eta Sgr
	{-2, 18097, -3042, "\xCE\xB3"},  // Gamma Sgr
	{-1, 18350, -2983, "\xCE\xB4"},  // Delta Sgr
	{-2, 18229, -2106, "\xCE\xBC"},  // Mu Sgr
	{-1, 18466, -2542, "\xCE\xBB"},  // Lambda Sgr
	{-2, 18962, -2111, "\xCE\xBE""2"}, // Xi2 Sgr
	{-1, 18921, -2630, "\xCF\x83"},  // Sigma Sgr

	// Scorpius (72)
	{-2, 17560, -3710, "\xCE\xBB"},  // Lambda Sco
	{-1, 17708, -3903, "\xCE\xBA"},  // Kappa Sco
	{-1, 17622, -4300, "\xCE\xB8"},  // Theta Sco
	{-1, 17203, -4324, "\xCE\xB7"},  // Eta Sco
	{-1, 16910, -4236, "\xCE\xB6" "2"}, // Zeta2 Sco
	{-1, 16864, -3805, "\xCE\xBC"},  // Mu 1 Sco
	{-1, 16836, -3429, "\xCE\xB5"},  // Epsilon Sco
	{-1, 16490, -2643, "\xCE\xB1"},  // Alpha Sco
	{-1, 16091, -1981, "\xCE\xB2"},  // Beta 1 Sco
	{-1, 16006, -2262, "\xCE\xB4"},  // Delta Sco
	{-1, 16490, -2643, "\xCE\xB1"},  // Alpha Sco

	// Sculptor (73)
	{-2,   977, -2936, "\xCE\xB1"},  // Alpha Scl
	{-1, 23815, -2813, "\xCE\xB4"},  // Delta Scl
	{-1, 23314, -3253, "\xCE\xB3"},  // Gamma Scl
	{-1, 23550, -3782, "\xCE\xB2"},  // Beta Scl

	// Scutum (74)
	{-2, 18786,  -475, "\xCE\xB2"},  // Beta Sct
	{-1, 18587,  -824, "\xCE\xB1"},  // Alpha Sct
	{-1, 18487, -1457, "\xCE\xB3"},  // Gamma Sct

	// Serpens (76) - Caput
	{-2, 15941,  1566, "\xCE\xB3"},  // Gamma Ser
	{-1, 15770,  1542, "\xCE\xB2"},  // Beta Ser
	{-1, 15580,  1054, "\xCE\xB4"},  // Delta Ser
	{-1, 15738,   643, "\xCE\xB1"},  // Alpha Ser
	{-1, 15847,   448, "\xCE\xB5"},  // Epsilon Ser
	{-1, 15827,  -343, "\xCE\xBC"},  // Mu Ser
	{-1, 16239,  -369, "\xCE\xB4"},  // Delta Oph (shared)

	// Serpens (76) - Cauda
	{-2, 17173, -1572, "\xCE\xB7"},  // Eta Oph (shared)
	{-1, 17626, -1540, "\xCE\xBE"},  // Xi Ser
	{-1, 18355,  -290, "\xCE\xB7"},  // Eta Ser
	{-1, 18937,   420, "\xCE\xB8"},  // Theta 1 Ser

	// Sextans (77)
	{-2, 10505,   -64, "\xCE\xB2"},  // Beta Sex
	{-1, 10132,   -37, "\xCE\xB1"},  // Alpha Sex
	{-1,  9875,  -811, "\xCE\xB3"},  // Gamma Sex

	// Taurus (78)
	{-2,  5627,  2114, "\xCE\xB6"},  // Zeta Tau
	{-1,  4599,  1651, "\xCE\xB1"},  // Alpha Tau
	{-1,  5438,  2861, "\xCE\xB2"},  // Beta Tau
	{-2,  4599,  1651, "\xCE\xB1"},  // Alpha Tau
	{-1,  4478,  1587, "\xCE\xB8" "2"}, // Theta2 Tau
	{-1,  4330,  1563, "\xCE\xB3"},  // Gamma Tau
	{-1,  4382,  1754, "\xCE\xB4"},  // Delta 1 Tau
	{-1,  4477,  1918, "\xCE\xB5"},  // Epsilon Tau
	{-2,  3791,  2411, "\xCE\xB7"},  // Eta Tau
	{-1,  4330,  1563, "\xCE\xB3"},  // Gamma Tau
	{-1,  4011,  1249, "\xCE\xBB"},  // Lambda Tau

	// Telescopium (79)
	{-2, 18187, -4595, "\xCE\xB5"},  // Epsilon Tel
	{-1, 18450, -4597, "\xCE\xB1"},  // Alpha Tel
	{-1, 18481, -4907, "\xCE\xB6"},  // Zeta Tel

	// Triangulum (80)
	{-2,  1885,  2958, "\xCE\xB1"},  // Alpha Tri
	{-1,  2159,  3499, "\xCE\xB2"},  // Beta Tri
	{-1,  2289,  3385, "\xCE\xB3"},  // Gamma Tri
	{-1,  1885,  2958, "\xCE\xB1"},  // Alpha Tri

	// Triangulum Australe (81)
	{-2, 16811, -6903, "\xCE\xB1"},  // Alpha TrA
	{-1, 15919, -6343, "\xCE\xB2"},  // Beta TrA
	{-1, 15315, -6868, "\xCE\xB3"},  // Gamma TrA
	{-1, 16811, -6903, "\xCE\xB1"},  // Alpha TrA

	// Tucana (82)
	{-2, 22308, -6026, "\xCE\xB1"},  // Alpha Tuc
	{-1, 23290, -5824, "\xCE\xB3"},  // Gamma Tuc
	{-1,   526, -6296, "\xCE\xB2"},  // Beta 1 Tuc
	{-1,   335, -6488, "\xCE\xB6"},  // Zeta Tuc
	{-1, 22308, -6026, "\xCE\xB1"},  // Alpha Tuc

	// Ursa Major (83)
	{-2, 13792,  4931, "\xCE\xB7"},  // Eta UMa
	{-1, 13399,  5493, "\xCE\xB6"},  // Zeta UMa
	{-1, 12900,  5596, "\xCE\xB5"},  // Epsilon UMa
	{-1, 12257,  5703, "\xCE\xB4"},  // Delta UMa
	{-1, 11062,  6175, "\xCE\xB1"},  // Alpha UMa
	{-1, 11031,  5638, "\xCE\xB2"},  // Beta UMa
	{-1, 11897,  5369, "\xCE\xB3"},  // Gamma UMa
	{-1, 12257,  5703, "\xCE\xB4"},  // Delta UMa

	// Ursa Minor (84)
	{-2,  2531,  8926, "\xCE\xB1"},  // Alpha UMi
	{-1, 17537,  8659, "\xCE\xB4"},  // Delta UMi
	{-1, 16766,  8204, "\xCE\xB5"},  // Epsilon UMi
	{-1, 15734,  7779, "\xCE\xB6"},  // Zeta UMi
	{-1, 14845,  7416, "\xCE\xB2"},  // Beta UMi
	{-1, 15345,  7183, "\xCE\xB3"},  // Gamma UMi
	{-1, 16292,  7576, "\xCE\xB7"},  // Eta UMi
	{-1, 15734,  7779, "\xCE\xB6"},  // Zeta UMi

	// Vela (85)
	{-2, 10779, -4942, "\xCE\xBC"},  // Mu Vel
	{-1,  9948, -5457, "\xCF\x86"},  // Phi Vel
	{-1,  9369, -5501, "\xCE\xBA"},  // Kappa Vel
	{-1,  8745, -5471, "\xCE\xB4"},  // Delta Vel
	{-1,  8158, -4735, "\xCE\xB3"},  // Gamma 1 Vel
	{-1,  9133, -4343, "\xCE\xBB"},  // Lambda Vel
	{-1,  9512, -4047, "\xCF\x88"},  // Psi Vel
	{-2,  9133, -4343, "\xCE\xBB"},  // Lambda Vel
	{-1,  9369, -5501, "\xCE\xBA"},  // Kappa Vel

	// Virgo (86)
	{-2, 14718,  -566, "\xCE\xBC"},  // Mu Vir
	{-1, 14267,  -600, "\xCE\xB9"},  // Iota Vir
	{-1, 13420, -1116, "\xCE\xB1"},  // Alpha Vir
	{-1, 13166,  -554, "\xCE\xB8"},  // Theta Vir
	{-1, 12694,  -145, "\xCE\xB3"},  // Gamma Vir
	{-1, 12332,   -67, "\xCE\xB7"},  // Eta Vir
	{-1, 11845,   176, "\xCE\xB2"},  // Beta Vir
	{-2, 12694,  -145, "\xCE\xB3"},  // Gamma Vir
	{-1, 12927,   340, "\xCE\xB4"},  // Delta Vir
	{-1, 13578,   -60, "\xCE\xB6"},  // Zeta Vir
	{-1, 13420, -1116, "\xCE\xB1"},  // Alpha Vir
	{-2, 13036,  1096, "\xCE\xB5"},  // Epsilon Vir
	{-1, 12927,   340, "\xCE\xB4"},  // Delta Vir

	// Volans (87)
	{-2,  9041, -6640, "\xCE\xB1"},  // Alpha Vol
	{-1,  8132, -6862, "\xCE\xB5"},  // Epsilon Vol
	{-1,  7281, -6796, "\xCE\xB4"},  // Delta Vol
	{-1,  7145, -7050, "\xCE\xB3"},  // Gamma 1 Vol
	{-1,  7697, -7261, "\xCE\xB6"},  // Zeta Vol
	{-1,  8132, -6862, "\xCE\xB5"},  // Epsilon Vol

	// Vulpecula (88)
	{-2, 20018,  2775, "15"},        // 15 Vul
	{-1, 19891,  2408, "13"},        // 13 Vul
	{-1, 19478,  2467, "\xCE\xB1"},  // Alpha Vul
	{-1, 19270,  2139, "1"},         // 1 Vul
}};

///----------------------------------------
/// MARK: Constellation label positions
///----------------------------------------

const std::array<std::array<int16_t, 2>, kConstellationCount> constPos = {{
	{  564,  3925},  // Andromeda
	{10118, -3365},  // Antlia
	{16000, -8000},  // Apus
	{22697, -1053},  // Aquarius
	{19690,   337},  // Aquila
	{17231, -5189},  // Ara
	{ 2676,  2257},  // Aries
	{ 5500,  4240},  // Auriga
	{14687,  3233},  // Bootes
	{ 4721, -3883},  // Caelum
	{ 6151,  7196},  // Camelopardalis
	{ 8497,  2356},  // Cancer
	{13020,  4235},  // Canes Venatici
	{ 6830, -2269},  // Canis Major
	{ 7624,   676},  // Canis Minor
	{21048, -1965},  // Capricornus
	{ 7761, -5717},  // Carina
	{  870,  6030},  // Cassiopeia
	{12950, -4400},  // Centaurus
	{22417,  7256},  // Cepheus
	{ 1709,  -664},  // Cetus
	{12000, -8000},  // Chamaeleon
	{14527, -6770},  // Circinus
	{ 5705, -3708},  // Columba
	{12748,  2265},  // Coma Berenices
	{18655, -4125},  // Corona Australis
	{15880,  3263},  // Corona Borealis
	{12387, -1836},  // Corvus
	{11348, -1325},  // Crater
	{12600, -6070},  // Crux
	{20598,  4958},  // Cygnus
	{20663,  1210},  // Delphinus
	{ 5333, -6399},  // Dorado
	{17945,  6606},  // Draco
	{21251,   794},  // Equuleus
	{ 3876, -1701},  // Eridanus
	{ 2766, -2694},  // Fornax
	{ 7300,  2600},  // Gemini
	{22457, -4586},  // Grus
	{17425,  3123},  // Hercules
	{ 3212, -5200},  // Horologium
	{ 9136, -1132},  // Hydra
	{ 2589, -7208},  // Hydrus
	{21138, -5268},  // Indus
	{22514,  4667},  // Lacerta
	{10600,  1800},  // Leo
	{10316,  3324},  // Leo Minor
	{ 5436, -1935},  // Lepus
	{15187, -1545},  // Libra
	{15376, -4228},  // Lupus
	{ 7734,  4783},  // Lynx
	{18908,  4065},  // Lyra
	{ 5500, -7999},  // Mensa
	{20942, -3620},  // Microscopium
	{ 6962,  -500},  // Monoceros
	{12460, -6987},  // Musca
	{16042, -5229},  // Norma
	{22173, -8473},  // Octans
	{17037,  -265},  // Ophiuchus
	{ 5660,   500},  // Orion
	{19160, -6514},  // Pavo
	{22617,  1965},  // Pegasus
	{ 3514,  4489},  // Perseus
	{  732, -4823},  // Phoenix
	{ 5381, -5163},  // Pictor
	{  891,  1548},  // Pisces
	{22410, -3143},  // Piscis Austrinus
	{ 7873, -3239},  // Puppis
	{ 8891, -2921},  // Pyxis
	{ 3899, -6049},  // Reticulum
	{19667,  1700},  // Sagitta
	{19385, -2911},  // Sagittarius
	{16865, -3567},  // Scorpius
	{  500, -3500},  // Sculptor
	{18651, -1011},  // Scutum
	{15731,  1085},  // Serpens Caput
	{17958, -1352},  // Serpens Cauda
	{10102,  -187},  // Sextans
	{ 4095,  1734},  // Taurus
	{19244, -5154},  // Telescopium
	{ 2043,  3234},  // Triangulum
	{16124, -6590},  // Triangulum Australe
	{23828, -6406},  // Tucana
	{10263,  5748},  // Ursa Major
	{15000,  7600},  // Ursa Minor
	{ 9337, -4851},  // Vela
	{13343,  -349},  // Virgo
	{ 7659, -6939},  // Volans
	{20367,  2503},  // Vulpecula
}};

///----------------------------------------
/// MARK: Constellation short names (IAU abbreviations)
///----------------------------------------

const std::array<std::string_view, kConstellationCount> constShortName = {{
	"And",  // Andromeda
	"Ant",  // Antlia
	"Aps",  // Apus
	"Aqr",  // Aquarius
	"Aql",  // Aquila
	"Ara",  // Ara
	"Ari",  // Aries
	"Aur",  // Auriga
	"Boo",  // Bootes
	"Cae",  // Caelum
	"Cam",  // Camelopardalis
	"Cnc",  // Cancer
	"CVn",  // Canes Venatici
	"CMa",  // Canis Major
	"CMi",  // Canis Minor
	"Cap",  // Capricornus
	"Car",  // Carina
	"Cas",  // Cassiopeia
	"Cen",  // Centaurus
	"Cep",  // Cepheus
	"Cet",  // Cetus
	"Cha",  // Chamaeleon
	"Cir",  // Circinus
	"Col",  // Columba
	"Com",  // Coma Berenices
	"CrA",  // Corona Australis
	"CrB",  // Corona Borealis
	"Crv",  // Corvus
	"Crt",  // Crater
	"Cru",  // Crux
	"Cyg",  // Cygnus
	"Del",  // Delphinus
	"Dor",  // Dorado
	"Dra",  // Draco
	"Equ",  // Equuleus
	"Eri",  // Eridanus
	"For",  // Fornax
	"Gem",  // Gemini
	"Gru",  // Grus
	"Her",  // Hercules
	"Hor",  // Horologium
	"Hya",  // Hydra
	"Hyi",  // Hydrus
	"Ind",  // Indus
	"Lac",  // Lacerta
	"Leo",  // Leo
	"LMi",  // Leo Minor
	"Lep",  // Lepus
	"Lib",  // Libra
	"Lup",  // Lupus
	"Lyn",  // Lynx
	"Lyr",  // Lyra
	"Men",  // Mensa
	"Mic",  // Microscopium
	"Mon",  // Monoceros
	"Mus",  // Musca
	"Nor",  // Norma
	"Oct",  // Octans
	"Oph",  // Ophiuchus
	"Ori",  // Orion
	"Pav",  // Pavo
	"Peg",  // Pegasus
	"Per",  // Perseus
	"Phe",  // Phoenix
	"Pic",  // Pictor
	"Psc",  // Pisces
	"PsA",  // Piscis Austrinus
	"Pup",  // Puppis
	"Pyx",  // Pyxis
	"Ret",  // Reticulum
	"Sge",  // Sagitta
	"Sgr",  // Sagittarius
	"Sco",  // Scorpius
	"Scl",  // Sculptor
	"Sct",  // Scutum
	"Ser",  // Serpens Caput
	"Ser",  // Serpens Cauda
	"Sex",  // Sextans
	"Tau",  // Taurus
	"Tel",  // Telescopium
	"Tri",  // Triangulum
	"TrA",  // Triangulum Australe
	"Tuc",  // Tucana
	"UMa",  // Ursa Major
	"UMi",  // Ursa Minor
	"Vel",  // Vela
	"Vir",  // Virgo
	"Vol",  // Volans
	"Vul",  // Vulpecula
}};

} // namespace
