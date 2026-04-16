/// @file contour.cpp
/// Streak detection helpers — faithful port of the pure-math functions from
/// ASTAP's unit_contour.pas.
///
/// Copyright (C) 2023, 2024 by Han Kleijn, www.hnsky.org (original Pascal).
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this file,
/// You can obtain one at https://mozilla.org/MPL/2.0/.

#include "contour.h"

#include <cmath>
#include <vector>

namespace astap::analysis {

double line_distance(double fits_x, double fits_y,
                     double slope, double intercept) {
	// Line: y = slope*x + intercept  =>  0 = slope*x - y + intercept
	// Distance = |slope*x - y + intercept| / sqrt(slope^2 + 1)
	// See https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
	return std::abs(slope * fits_x - fits_y + intercept)
	       / std::sqrt(slope * slope + 1.0);
}

void trendline(const StarList& xylist, int len,
               double& slope, double& intercept) {
	// Standard least-squares linear regression: Y = slope * X + intercept.
	// Idea from https://stackoverflow.com/questions/43224/how-do-i-calculate-a-trendline-for-a-graph
	double sum_x  = 0.0;
	double sum_x2 = 0.0;
	double sum_y  = 0.0;
	double sum_xy = 0.0;

	for (int i = 0; i < len; ++i) {
		const double x = xylist[0][static_cast<std::size_t>(i)];
		const double y = xylist[1][static_cast<std::size_t>(i)];
		sum_x  += x;
		sum_x2 += x * x;
		sum_y  += y;
		sum_xy += x * y;
	}

	const double n = static_cast<double>(len);
	slope     = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
	intercept = (sum_y - slope * sum_x) / n;
}

void trendline_without_outliers(const StarList& xylist, int len,
                                double& slope, double& intercept, double& sd) {
	// Pass 1: ordinary least-squares fit.
	trendline(xylist, len, slope, intercept);

	// Compute standard deviation of perpendicular distances.
	// Distance^2 = (slope*x - y + intercept)^2 / (slope^2 + 1).
	const double denom = slope * slope + 1.0;
	double sum_sq = 0.0;
	for (int i = 0; i < len; ++i) {
		const double x = xylist[0][static_cast<std::size_t>(i)];
		const double y = xylist[1][static_cast<std::size_t>(i)];
		const double residual = slope * x - y + intercept;
		sum_sq += (residual * residual) / denom;
	}
	sd = std::sqrt(sum_sq / static_cast<double>(len));

	// Pass 2: rebuild data excluding outliers beyond 1.5 * sd.
	StarList filtered(2);
	filtered[0].reserve(static_cast<std::size_t>(len));
	filtered[1].reserve(static_cast<std::size_t>(len));

	for (int i = 0; i < len; ++i) {
		const double x = xylist[0][static_cast<std::size_t>(i)];
		const double y = xylist[1][static_cast<std::size_t>(i)];
		const double error = std::abs(y - (slope * x + intercept));
		if (error < 1.5 * sd) {
			filtered[0].push_back(x);
			filtered[1].push_back(y);
		}
	}

	// Re-fit on the filtered data.
	trendline(filtered, static_cast<int>(filtered[0].size()), slope, intercept);
}

}  // namespace astap::analysis
