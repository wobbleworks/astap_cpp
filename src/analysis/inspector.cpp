/// @file inspector.cpp
/// CCD inspector helpers — ported from unit_inspector_plot.pas.

#include "inspector.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>

namespace astap::analysis {

// ── measure_star_aspect ─────────────────────────────────────────────

StarAspect measure_star_aspect(const ImageArray& img,
                               double x1, double y1, int rs,
                               double star_bg, double sd_bg)
{
	static constexpr int kMaxRadius = 51;
	static constexpr double kPi    = std::numbers::pi;

	const auto h = static_cast<int>(img[0].size());      // height
	const auto w = static_cast<int>(img[0][0].size());    // width

	rs = std::min(rs, kMaxRadius);

	// Bounds check: the search box must fit inside the image.
	if (x1 - rs < 0 || x1 + rs >= w || y1 - rs <= 0 || y1 + rs >= h)
		return {};  // aspect = 999 (failure)

	// Bilinear sub-pixel interpolation on channel 0.
	auto value_subpixel = [&](double sx, double sy) -> double {
		const auto x_trunc = static_cast<int>(sx);
		const auto y_trunc = static_cast<int>(sy);
		if (x_trunc <= 0 || x_trunc >= w - 2 || y_trunc <= 0 || y_trunc >= h - 2)
			return 0.0;
		const double x_frac = sx - x_trunc;
		const double y_frac = sy - y_trunc;
		return img[0][y_trunc    ][x_trunc    ] * (1.0 - x_frac) * (1.0 - y_frac)
		     + img[0][y_trunc    ][x_trunc + 1] * (       x_frac) * (1.0 - y_frac)
		     + img[0][y_trunc + 1][x_trunc    ] * (1.0 - x_frac) * (       y_frac)
		     + img[0][y_trunc + 1][x_trunc + 1] * (       x_frac) * (       y_frac);
	};

	// Build weighted-distance map.  data[i][j] = sqrt(val) * r, where val
	// is the background-subtracted pixel value and r is the distance from
	// the centroid, but only where val > 7 * sd_bg.
	constexpr int kSide = 2 * kMaxRadius + 1;
	std::array<std::array<double, kSide>, kSide> data{};

	const double threshold = 7.0 * sd_bg;
	int pixel_counter = 0;

	for (int i = -rs; i <= rs; ++i) {
		for (int j = -rs; j <= rs; ++j) {
			const double val = value_subpixel(x1 + i, y1 + j) - star_bg;
			if (val > threshold) {
				const double r = std::sqrt(static_cast<double>(i * i + j * j));
				data[static_cast<std::size_t>(i + kMaxRadius)]
				    [static_cast<std::size_t>(j + kMaxRadius)] = std::sqrt(val) * r;
				++pixel_counter;
			}
			// else: zero-initialised by std::array value-init
		}
	}

	if (pixel_counter < 4)
		return {};  // not enough pixels

	// Sweep 0-179 degrees.  For each angle, compute the sum of
	// |distance-to-line| weighted by the data map.  The line passes through
	// the centroid at the given angle.  Because the data already includes
	// the radial distance r, the perpendicular distance is
	// data[i][j] * |sin(angle - atan2(j, i))|.
	double themax = 0.0;
	double themin = 1e99;
	int orientation_min = 0;

	for (int angle = 0; angle < 180; ++angle) {
		const double angle_rad = angle * kPi / 180.0;
		double distance = 0.0;

		for (int i = -rs; i <= rs; ++i) {
			for (int j = -rs; j <= rs; ++j) {
				const double d = data[static_cast<std::size_t>(i + kMaxRadius)]
				                     [static_cast<std::size_t>(j + kMaxRadius)];
				if (d > 0.0) {
					const double g = std::atan2(static_cast<double>(j),
					                            static_cast<double>(i));
					const double delta_angle = angle_rad - g;
					distance += d * std::abs(std::sin(delta_angle));
				}
			}
		}

		if (distance > themax)
			themax = distance;
		if (distance < themin) {
			themin = distance;
			orientation_min = angle;
		}
	}

	const double aspect = themax / (themin + 1e-5);
	if (aspect > 5.0)
		return {};  // failure

	return {.aspect = aspect, .orientation = orientation_min};
}

// ── filter_hfd ──────────────────────────────────────────────────────

void filter_hfd(StarList& hfd_values, int nr,
                float& mean, float& min_value, float& max_value)
{
	max_value = 0.0f;
	min_value = 65535.0f;
	mean      = 0.0f;

	for (int i = 0; i < nr; ++i) {
		double closest_dist    = 999999.0;
		double second_closest  = 999999.0;
		int nr_closest         = 0;
		int nr_second_closest  = 0;

		for (int j = 0; j < nr; ++j) {
			if (i == j)
				continue;
			const double dx = hfd_values[0][i] - hfd_values[0][j];
			const double dy = hfd_values[1][i] - hfd_values[1][j];
			const double dist_sqr = dx * dx + dy * dy;

			if (dist_sqr < closest_dist) {
				second_closest    = closest_dist;
				nr_second_closest = nr_closest;
				closest_dist      = dist_sqr;
				nr_closest        = j;
			} else if (dist_sqr < second_closest) {
				second_closest    = dist_sqr;
				nr_second_closest = j;
			}
		}

		// Median of the three HFD values (current star + two nearest).
		double a = hfd_values[2][i];
		double b = hfd_values[2][nr_closest];
		double c = hfd_values[2][nr_second_closest];

		// Sort so that a >= b >= c.
		if (a < b) std::swap(a, b);
		if (a < c) std::swap(a, c);
		if (b < c) std::swap(b, c);

		hfd_values[2][i] = b;  // median

		if (b > max_value) max_value = static_cast<float>(b);
		if (b < min_value) min_value = static_cast<float>(b);
		mean += static_cast<float>(b);
	}

	if (nr > 0)
		mean /= static_cast<float>(nr);
}

}  // namespace astap::analysis
