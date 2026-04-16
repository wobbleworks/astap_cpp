///----------------------------------------
///      @file gaussian_blur.cpp
///   @ingroup ASTAP++
///     @brief Separable Gaussian blur implementation.
///    @author Ported from Han Kleijn's unit_stack.pas (ASTAP). MPL-2.0.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "gaussian_blur.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

///----------------------------------------
namespace astap::stacking {
///----------------------------------------

///----------------------------------------
/// MARK: File-local helpers
///----------------------------------------

namespace {

constexpr auto kMaxKernelSize = 100;
constexpr auto kDataGranularity = 1.0;
constexpr auto kMaxData = 255.0;

///----------------------------------------
///   @brief 1-D Gaussian convolution kernel.
/// @details Weights are stored 0-indexed at @c weights[kMaxKernelSize + offset]
///          so that valid offsets in @c [-size, +size] map to array slots
///          @c [kMaxKernelSize - size, kMaxKernelSize + size].
///----------------------------------------

struct TKernel {
    int size = 1;
    std::array<float, 2 * kMaxKernelSize + 1> weights{};
    
    [[nodiscard]] float& weight_at(int offset) noexcept {
        return weights[static_cast<std::size_t>(kMaxKernelSize + offset)];
    }
    
    [[nodiscard]] float weight_at(int offset) const noexcept {
        return weights[static_cast<std::size_t>(kMaxKernelSize + offset)];
    }
};

///----------------------------------------
///   @brief Build a normalized Gaussian kernel with standard deviation @p radius.
/// @details Fills the full @c 2*kMaxKernelSize+1 window with the unnormalized
///          Gaussian, normalizes to unit sum, then trims the kernel size until
///          the discarded tail exceeds @c data_granularity / (2 * max_data).
///          The surviving entries are re-normalized so they sum to exactly 1.
/// @param[out] K Kernel to populate.
///   @param radius Gaussian standard deviation in pixels.
///   @param max_data Maximum expected sample value (used to pick a trim threshold).
///   @param data_granularity Smallest meaningful sample increment.
///----------------------------------------

void make_gaussian_kernel(TKernel& K, double radius,
                          double max_data, double data_granularity) {
    // Fill all weights with the unnormalized Gaussian.
    for (int j = -kMaxKernelSize; j <= kMaxKernelSize; ++j) {
        const auto t = static_cast<double>(j) / radius;
        K.weight_at(j) = static_cast<float>(std::exp(-t * t / 2.0));
    }
    
    // Normalize so the full kernel sums to 1.
    auto sum = 0.0;
    for (int j = -kMaxKernelSize; j <= kMaxKernelSize; ++j) {
        sum += K.weight_at(j);
    }
    for (int j = -kMaxKernelSize; j <= kMaxKernelSize; ++j) {
        K.weight_at(j) = static_cast<float>(K.weight_at(j) / sum);
    }
    
    // Trim insignificant tails: shrink size until the discarded mass exceeds
    // delta = data_granularity / (2 * max_data).
    auto kernel_size = kMaxKernelSize;
    const auto delta = data_granularity / (2.0 * max_data);
    auto trimmed = 0.0;
    while (trimmed < delta && kernel_size > 1) {
        trimmed += 2.0 * K.weight_at(kernel_size);
        --kernel_size;
    }
    K.size = kernel_size;
    
    // Re-normalize the surviving entries so they sum to exactly 1.
    sum = 0.0;
    for (int j = -K.size; j <= K.size; ++j) {
        sum += K.weight_at(j);
    }
    for (int j = -K.size; j <= K.size; ++j) {
        K.weight_at(j) = static_cast<float>(K.weight_at(j) / sum);
    }
}

///----------------------------------------
///   @brief Horizontal convolution pass.
/// @details Convolves @p src rows into @p dst using kernel @p K, clamping out-of-
///          range columns to the image edge. @p h and @p w are inclusive maxima
///          (height-1, width-1).
///   @param src Source image buffer.
/// @param[out] dst Destination image buffer (same shape as @p src).
///   @param K Gaussian kernel.
///   @param h Largest valid row index (height - 1).
///   @param w Largest valid column index (width - 1).
///   @param colors Number of active channels (1, 2, or 3).
///----------------------------------------

void blur_h(const ImageArray& src, ImageArray& dst, const TKernel& K,
            int h, int w, int colors) {
    for (int i = 0; i <= h; ++i) {
        for (int j = 0; j <= w; ++j) {
            auto valr = 0.0f;
            auto valg = 0.0f;
            auto valb = 0.0f;
            for (int jx = -K.size; jx <= K.size; ++jx) {
                auto x = j + jx;
                if (x < 0) {
                    x = 0;
                }
                if (x > w) {
                    x = w;
                }
                
                const auto weight = K.weight_at(jx);
                valr += src[0][i][x] * weight;
                if (colors >= 2) {
                    valg += src[1][i][x] * weight;
                }
                if (colors >= 3) {
                    valb += src[2][i][x] * weight;
                }
            }
            dst[0][i][j] = valr;
            if (colors >= 2) {
                dst[1][i][j] = valg;
            }
            if (colors >= 3) {
                dst[2][i][j] = valb;
            }
        }
    }
}

///----------------------------------------
///   @brief Vertical convolution pass.
/// @details Convolves @p src columns into @p dst using kernel @p K, clamping
///          out-of-range rows to the image edge.
///   @param src Source image buffer.
/// @param[out] dst Destination image buffer (same shape as @p src).
///   @param K Gaussian kernel.
///   @param h Largest valid row index (height - 1).
///   @param w Largest valid column index (width - 1).
///   @param colors Number of active channels (1, 2, or 3).
///----------------------------------------

void blur_v(const ImageArray& src, ImageArray& dst, const TKernel& K,
            int h, int w, int colors) {
    for (int i = 0; i <= h; ++i) {
        for (int j = 0; j <= w; ++j) {
            auto valr = 0.0f;
            auto valg = 0.0f;
            auto valb = 0.0f;
            for (int iy = -K.size; iy <= K.size; ++iy) {
                auto y = i + iy;
                if (y < 0) {
                    y = 0;
                }
                if (y > h) {
                    y = h;
                }
                
                const auto weight = K.weight_at(iy);
                valr += src[0][y][j] * weight;
                if (colors >= 2) {
                    valg += src[1][y][j] * weight;
                }
                if (colors >= 3) {
                    valb += src[2][y][j] * weight;
                }
            }
            dst[0][i][j] = valr;
            if (colors >= 2) {
                dst[1][i][j] = valg;
            }
            if (colors >= 3) {
                dst[2][i][j] = valb;
            }
        }
    }
}
    
} // namespace

///----------------------------------------
/// MARK: Public API
///----------------------------------------

void gaussian_blur2(ImageArray& img, double radius) {
    // Skip tiny radii to avoid division-by-zero when building the kernel.
    if (radius < 0.001) {
        return;
    }
    if (img.empty() || img[0].empty() || img[0][0].empty()) {
        return;
    }
    
    // Build the trimmed, normalized 1-D Gaussian kernel.
    TKernel K;
    make_gaussian_kernel(K, radius, kMaxData, kDataGranularity);
    
    // Grab the buffer shape.
    const auto colors = static_cast<int>(img.size());
    const auto h = static_cast<int>(img[0].size());
    const auto w = static_cast<int>(img[0][0].size());
    
    // Scratch buffer matching img's shape: [colors][h][w].
    auto tmp = ImageArray(colors,
                          std::vector<std::vector<float>>(
                              h, std::vector<float>(w, 0.0f)));
                              
    // Separable convolution: horizontal into scratch, then vertical back into img.
    blur_h(img, tmp, K, h - 1, w - 1, colors);
    blur_v(tmp, img, K, h - 1, w - 1, colors);
}
    
} // namespace
