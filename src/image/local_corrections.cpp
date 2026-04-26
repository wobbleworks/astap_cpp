///----------------------------------------
///      @file local_corrections.cpp
///   @ingroup ASTAP++
///     @brief Implementation of @ref remove_linear_gradient and
///            @ref remove_dust_spot.
///    @author Ported from Han Kleijn's ASTAP (MPL-2.0) by John Stephen.
/// @copyright Copyright (C) Han Kleijn / John Stephen. Mozilla Public License 2.0.
///----------------------------------------

#include "local_corrections.h"

#include "../core/photometry.h"
#include "../stacking/gaussian_blur.h"

#include <algorithm>
#include <cmath>

///----------------------------------------
namespace astap::image {
///----------------------------------------

namespace {

constexpr auto kPatchHalfSize = 20;
constexpr auto kModeCeiling   = 65535;
constexpr auto kDustModeCeil  = 32000;
constexpr auto kBlurRadius    = 5.0;

[[nodiscard]] double mode_at(const ImageArray& img, int channel,
                             int xc, int yc, int half, int ceiling) {
	auto greylevels = 0;
	return astap::core::mode(img, /*ellipse=*/false, channel,
	                          xc - half, xc + half,
	                          yc - half, yc + half,
	                          ceiling, greylevels);
}

}  // namespace

bool remove_linear_gradient(ImageArray& img, const Header& head,
                            int dark_x, int dark_y,
                            int bright_x, int bright_y) {
	if (head.width <= 0 || head.height <= 0) return false;
	if (std::abs(bright_x - dark_x) <= 100 &&
	    std::abs(bright_y - dark_y) <= 100) {
		return false;  // Pascal's "OR" guard: at least one axis must be far apart.
	}

	const auto channels = std::min<int>(head.naxis3, static_cast<int>(img.size()));
	if (channels <= 0) return false;

	auto dark_mode   = std::array<double, 3>{};
	auto bright_mode = std::array<double, 3>{};
	for (auto c = 0; c < channels && c < 3; ++c) {
		dark_mode[c]   = mode_at(img, c, dark_x,   dark_y,   kPatchHalfSize, kModeCeiling);
		bright_mode[c] = mode_at(img, c, bright_x, bright_y, kPatchHalfSize, kModeCeiling);
	}

	const auto a = std::sqrt(
		static_cast<double>(bright_x - dark_x) * (bright_x - dark_x) +
		static_cast<double>(bright_y - dark_y) * (bright_y - dark_y));
	if (a == 0.0) return false;

	for (auto y = 0; y < head.height; ++y) {
		for (auto x = 0; x < head.width; ++x) {
			const auto bx = static_cast<double>(x - dark_x);
			const auto by = static_cast<double>(y - dark_y);
			const auto cx = static_cast<double>(x - bright_x);
			const auto cy = static_cast<double>(y - bright_y);
			const auto b  = std::sqrt(bx * bx + by * by);
			const auto cc = std::sqrt(cx * cx + cy * cy);
			// Orthogonal projection of P onto the dark→bright line, expressed
			// as distance along the line measured from the dark anchor.
			const auto p = -((b * b - a * a - cc * cc) / (2.0 * a));
			const auto t = (a - p) / a;
			for (auto c = 0; c < channels && c < 3; ++c) {
				const auto delta = (bright_mode[c] - dark_mode[c]) * t;
				img[c][y][x] = static_cast<float>(img[c][y][x] - delta);
			}
		}
	}
	return true;
}

bool remove_dust_spot(ImageArray& img, const Header& head,
                      int x1, int y1, int x2, int y2,
                      bool use_ellipse) {
	if (head.width <= 0 || head.height <= 0) return false;
	auto sx = std::min(x1, x2);
	auto ex = std::max(x1, x2);
	auto sy = std::min(y1, y2);
	auto ey = std::max(y1, y2);
	if ((ex - sx) <= 3 || (ey - sy) <= 3) return false;

	// Clip to image bounds plus a kPatchHalfSize margin so the corner patches
	// stay valid. (mode() itself clamps but giving it a bad rect would produce
	// garbage modes.)
	sx = std::clamp(sx, kPatchHalfSize, head.width  - kPatchHalfSize - 1);
	ex = std::clamp(ex, kPatchHalfSize, head.width  - kPatchHalfSize - 1);
	sy = std::clamp(sy, kPatchHalfSize, head.height - kPatchHalfSize - 1);
	ey = std::clamp(ey, kPatchHalfSize, head.height - kPatchHalfSize - 1);
	if (ex - sx <= 3 || ey - sy <= 3) return false;

	const auto cx = (sx + ex) / 2.0;
	const auto cy = (sy + ey) / 2.0;
	const auto rx = (ex - 1 - sx) / 2.0;
	const auto ry = (ey - 1 - sy) / 2.0;
	if (rx <= 0.0 || ry <= 0.0) return false;

	const auto w = ex - sx;
	const auto h = ey - sy;
	const auto channels = std::min<int>(head.naxis3, static_cast<int>(img.size()));
	if (channels <= 0) return false;

	// Reusable scratch buffer matching gaussian_blur2's signature.
	auto delta = ImageArray{};
	delta.assign(1, std::vector<std::vector<float>>(
		h, std::vector<float>(w, 0.0f)));

	for (auto k = 0; k < channels; ++k) {
		// Sample the four corner-region modes just outside the box.
		const auto m_lb = mode_at(img, k, sx - kPatchHalfSize / 2, sy - kPatchHalfSize / 2,
		                          kPatchHalfSize / 2, kDustModeCeil);
		const auto m_lt = mode_at(img, k, sx - kPatchHalfSize / 2, ey + kPatchHalfSize / 2,
		                          kPatchHalfSize / 2, kDustModeCeil);
		const auto m_rb = mode_at(img, k, ex + kPatchHalfSize / 2, sy - kPatchHalfSize / 2,
		                          kPatchHalfSize / 2, kDustModeCeil);
		const auto m_rt = mode_at(img, k, ex + kPatchHalfSize / 2, ey + kPatchHalfSize / 2,
		                          kPatchHalfSize / 2, kDustModeCeil);

		const auto box_w = static_cast<double>(ex - sx);
		const auto box_h = static_cast<double>(ey - sy);

		for (auto y = sy; y < ey; ++y) {
			for (auto x = sx; x < ex; ++x) {
				const auto fx = (ex - x) / box_w;
				const auto fy = (ey - y) / box_h;
				const auto line_bottom = m_lb * fx + m_rb * (1.0 - fx);
				const auto line_top    = m_lt * fx + m_rt * (1.0 - fx);
				const auto expected    = line_bottom * fy + line_top * (1.0 - fy);
				// The Pascal uses max(0, expected − pixel) so light contamination
				// (a brighter star inside the patch) is left alone.
				const auto raw = expected - img[k][y][x];
				delta[0][y - sy][x - sx] = static_cast<float>(std::max(0.0, raw));
			}
		}

		astap::stacking::gaussian_blur2(delta, kBlurRadius);

		for (auto y = sy; y < ey; ++y) {
			for (auto x = sx; x < ex; ++x) {
				if (use_ellipse) {
					const auto u = (x - cx) / rx;
					const auto v = (y - cy) / ry;
					if (u * u + v * v >= 1.0) continue;
				}
				img[k][y][x] = static_cast<float>(
					img[k][y][x] + delta[0][y - sy][x - sx]);
			}
		}
	}
	return true;
}

} // namespace astap::image
