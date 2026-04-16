/// @file image_sharpness.cpp
/// Image sharpness measurement — ported from unit_image_sharpness.pas.

#include "image_sharpness.h"

#include <algorithm>
#include <cmath>

namespace astap::analysis {

double image_sharpness(const ImageArray& img)
{
	const auto h = static_cast<int>(img[0].size());      // height
	const auto w = static_cast<int>(img[0][0].size());    // width

	double result  = 0.0;
	double average = 0.0;

	// Process 4x4 pixel blocks (two 2x2 RGGB quads per axis).
	for (int i = 0; i <= (h - 4) / 4; ++i) {
		for (int j = 0; j <= (w - 4) / 4; ++j) {
			// Bottom-left 2x2 quad sum.
			const double v1_lo = img[0][i * 4    ][j * 4    ]
			                   + img[0][i * 4 + 1][j * 4    ]
			                   + img[0][i * 4    ][j * 4 + 1]
			                   + img[0][i * 4 + 1][j * 4 + 1];

			// Bottom-right 2x2 quad sum (offset +2 in y).
			const double v2_lo = img[0][i * 4 + 2][j * 4    ]
			                   + img[0][i * 4 + 3][j * 4    ]
			                   + img[0][i * 4 + 2][j * 4 + 1]
			                   + img[0][i * 4 + 3][j * 4 + 1];

			const double max_a = std::max(v1_lo, v2_lo);
			const double min_a = std::min(v1_lo, v2_lo);

			// Top-left 2x2 quad sum (offset +2 in x).
			const double v1_hi = img[0][i * 4    ][j * 4 + 2]
			                   + img[0][i * 4 + 1][j * 4 + 2]
			                   + img[0][i * 4    ][j * 4 + 3]
			                   + img[0][i * 4 + 1][j * 4 + 3];

			// Top-right 2x2 quad sum (offset +2 in both axes).
			const double v2_hi = img[0][i * 4 + 2][j * 4 + 2]
			                   + img[0][i * 4 + 3][j * 4 + 2]
			                   + img[0][i * 4 + 2][j * 4 + 3]
			                   + img[0][i * 4 + 3][j * 4 + 3];

			const double max_b = std::max(v1_hi, v2_hi);
			const double min_b = std::min(v1_hi, v2_hi);

			const double minimum = std::min(min_a, min_b);
			const double maximum = std::max(max_a, max_b);

			const double diff = maximum - minimum;
			result  += diff * diff;
			average += maximum + minimum;
		}
	}

	const double area = static_cast<double>(w) * h;
	result  = std::sqrt(result / area);
	average = average / area;

	// Invert and scale to approximate HFD: lower = sharper.
	return 4.0 * average / (result + 1e-18);
}

}  // namespace astap::analysis
