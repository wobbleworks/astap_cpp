/// @file contour.cpp
/// Streak detection helpers — faithful port of the pure-math functions from
/// ASTAP's unit_contour.pas.
///
/// Copyright (C) 2023, 2024 by Han Kleijn, www.hnsky.org (original Pascal).
/// This Source Code Form is subject to the terms of the Mozilla Public
/// License, v. 2.0. If a copy of the MPL was not distributed with this file,
/// You can obtain one at https://mozilla.org/MPL/2.0/.

#include "contour.h"

#include <algorithm>
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

std::vector<Streak> detect_streaks(const ImageArray& img_bk,
                                   double detection_level,
                                   int binning, int detection_grid) {
	std::vector<Streak> streaks;
	if (img_bk.empty() || img_bk[0].empty() || img_bk[0][0].empty()) {
		return streaks;
	}
	const int hh = static_cast<int>(img_bk[0].size());
	const int ww = static_cast<int>(img_bk[0][0].size());

	// img_sa: -1 = untested/star-free, >= 0 = traced boundary or filled interior.
	std::vector<std::vector<int>> img_sa(
		static_cast<std::size_t>(hh), std::vector<int>(static_cast<std::size_t>(ww), -1));

	// Moore-Neighbor tracing tables (unit_contour.pas).
	static constexpr int newdirection[8] = {-1, 0, 0, +1, +1, +2, +2, -1};
	static constexpr int directions[8][2] = {
		{-1, -1}, {-1, 0}, {-1, +1}, {0, +1}, {+1, +1}, {+1, 0}, {+1, -1}, {0, -1}};

	auto img_protected = [&](int xx, int yy) -> bool {
		if (xx >= 0 && xx < ww - 1 && yy >= 0 && yy < hh - 1) {
			return img_bk[0][static_cast<std::size_t>(yy)][static_cast<std::size_t>(xx)]
			       > detection_level;
		}
		return false;
	};

	// Boundary points traced for the current blob.
	std::vector<int> cx, cy;
	cx.reserve(static_cast<std::size_t>(4 * ww));
	cy.reserve(static_cast<std::size_t>(4 * ww));

	auto find_contour = [&](int fx, int fy) {
		int direction = 1;  // north=0, west=1, south=2, east=3
		const int startX = fx;
		const int startY = fy;
		int counter = 0;
		cx.clear();
		cy.clear();

		while (true) {
			bool detection = false;
			for (int i = 0; i < 8; ++i) {
				const int j = (i + direction * 2) & 0x7;
				if (img_protected(fx + directions[j][0], fy + directions[j][1])) {
					fx += directions[j][0];
					fy += directions[j][1];
					detection = true;
					direction += newdirection[i];
					break;
				}
			}
			if (!detection) {
				break;
			}
			cx.push_back(fx);
			cy.push_back(fy);

			img_sa[static_cast<std::size_t>(fy)][static_cast<std::size_t>(fx)] += 1;
			if (img_sa[static_cast<std::size_t>(fy)][static_cast<std::size_t>(fx)] > 1) {
				break;  // looping locally
			}
			++counter;
			if ((fx == startX && fy == startY) || counter > 4 * ww) {
				break;
			}
		}

		const int counterC = static_cast<int>(cx.size());

		// Mark the interior of the contour (span-fill per shared row) and
		// measure its surface area plus bounding box.
		double surface = 0.0;
		int maxX = 0, minX = 999999, maxY = 0, minY = 999999;
		for (int i = 0; i < counterC; ++i) {
			minX = std::min(cx[static_cast<std::size_t>(i)], minX);
			maxX = std::max(cx[static_cast<std::size_t>(i)], maxX);
			minY = std::min(cy[static_cast<std::size_t>(i)], minY);
			maxY = std::max(cy[static_cast<std::size_t>(i)], maxY);
			for (int j = 0; j < counterC; ++j) {
				if (cy[static_cast<std::size_t>(i)] == cy[static_cast<std::size_t>(j)]) {
					const int k0 = std::min(cx[static_cast<std::size_t>(i)],
					                        cx[static_cast<std::size_t>(j)]);
					const int k1 = std::max(cx[static_cast<std::size_t>(i)],
					                        cx[static_cast<std::size_t>(j)]);
					const auto row = static_cast<std::size_t>(cy[static_cast<std::size_t>(i)]);
					for (int k = k0; k <= k1; ++k) {
						if (img_sa[row][static_cast<std::size_t>(k)] < 0) {
							surface += 1.0;
							img_sa[row][static_cast<std::size_t>(k)] = 1;
						}
					}
				}
			}
		}

		if (surface > 200 * 2) {
			const double maxleng = std::sqrt(
				static_cast<double>(maxY - minY) * (maxY - minY) +
				static_cast<double>(maxX - minX) * (maxX - minX));
			if (maxleng > 200 && maxleng * maxleng / surface > 10) {
				StarList contour_pts(2);
				contour_pts[0].resize(static_cast<std::size_t>(counterC));
				contour_pts[1].resize(static_cast<std::size_t>(counterC));
				for (int i = 0; i < counterC; ++i) {
					contour_pts[0][static_cast<std::size_t>(i)] = cx[static_cast<std::size_t>(i)];
					contour_pts[1][static_cast<std::size_t>(i)] = cy[static_cast<std::size_t>(i)];
				}
				double slope = 0.0, intercept = 0.0, sd = 0.0;
				trendline_without_outliers(contour_pts, counterC, slope, intercept, sd);
				intercept *= binning;  // back to full-resolution coordinates
				sd *= binning;
				if (sd < 10) {  // a real line: sigma ~ line thickness + a nearby star
					streaks.push_back(Streak{slope, intercept});
				}
			}
		}
	};

	for (int fitsY = 0; fitsY < hh; ++fitsY) {
		for (int fitsX = 0; fitsX < ww; ++fitsX) {
			// Sample along an overlay grid of vertical + horizontal lines.
			if (detection_grid <= 0 || (fitsX % detection_grid == 0) ||
			    (fitsY % detection_grid == 0)) {
				if (img_sa[static_cast<std::size_t>(fitsY)][static_cast<std::size_t>(fitsX)] < 0 &&
				    img_bk[0][static_cast<std::size_t>(fitsY)][static_cast<std::size_t>(fitsX)]
				        > detection_level) {
					find_contour(fitsX, fitsY);
				}
			}
		}
	}
	return streaks;
}

}  // namespace astap::analysis
