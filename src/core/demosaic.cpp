///----------------------------------------
///      @file demosaic.cpp
///   @ingroup ASTAP++
///     @brief Bayer / X-Trans demosaic implementation.
///   @details Faithful port of astap_main.pas demosaic procedures (~6707-7651).
///            Loop structure and numeric ordering are preserved so the math
///            matches the original output. Header mutations (naxis, naxis3,
///            width, height) match the originals.
///    @author Ported from Han Kleijn's astap_main.pas (ASTAP), MPL-2.0.
///            C++ port by John Stephen.
/// @copyright Copyright (C) Han Kleijn / John Stephen.
///            Mozilla Public License 2.0.
///----------------------------------------

#include "demosaic.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <utility>

///----------------------------------------
namespace astap::core {
///----------------------------------------

namespace {

/// MARK: - Helpers

///----------------------------------------
/// @brief Test whether an integer is odd.
///----------------------------------------

[[nodiscard]] constexpr bool is_odd(int v) noexcept {
    return (v & 1) != 0;
}

///----------------------------------------
/// @brief Allocate a fresh 3 x H x W zero-initialised image buffer.
///----------------------------------------

[[nodiscard]] ImageArray make_rgb_buffer(int height, int width) {
    return ImageArray(3,
                      std::vector<std::vector<float>>(
                          static_cast<std::size_t>(height),
                          std::vector<float>(static_cast<std::size_t>(width), 0.0f)));
}

///----------------------------------------
/// @brief Look up Bayer pattern offsets.
/// @param pattern  Bayer pattern index (0-3).
/// @param[out] offsetx  X offset for the pattern.
/// @param[out] offsety  Y offset for the pattern.
/// @return True if the pattern is in 0-3; false otherwise.
///----------------------------------------

[[nodiscard]] bool pattern_offsets(int pattern, int& offsetx, int& offsety) noexcept {
    switch (pattern) {
        case 0: offsetx = 0; offsety = 0; return true;  // GRBG
        case 1: offsetx = 0; offsety = 1; return true;  // BGGR
        case 2: offsetx = 1; offsety = 0; return true;  // RGGB
        case 3: offsetx = 1; offsety = 1; return true;  // GBRG
        default: return false;
    }
}

/// Saturated-pixel marker used by the astroC algorithm (0xFFFFFF).
constexpr auto kSaturatedMarker = 16777215.0f;
    
}  // namespace

/// MARK: - demosaic_bilinear_interpolation

void demosaic_bilinear_interpolation(ImageArray& img, Header& head, int pattern) {
    auto offsetx = 0;
    auto offsety = 0;
    if (!pattern_offsets(pattern, offsetx, offsety)) {
        return;
    }
    
    const auto H = head.height;
    const auto W = head.width;
    auto tmp = make_rgb_buffer(H, W);
    
    for (auto y = 1; y <= H - 2; ++y) {
        for (auto x = 1; x <= W - 2; ++x) {
            const auto green_even = is_odd(x + 1 + offsetx) && is_odd(y + 1 + offsety);
            const auto green_odd  = is_odd(x + offsetx)     && is_odd(y + offsety);
            const auto red        = is_odd(x + offsetx)     && is_odd(y + 1 + offsety);
            const auto blue       = is_odd(x + 1 + offsetx) && is_odd(y + offsety);
            
            if (green_odd) {
                tmp[0][y][x] = (img[0][y - 1][x]     + img[0][y + 1][x]    ) / 2.0f;
                tmp[1][y][x] =  img[0][y][x];
                tmp[2][y][x] = (img[0][y][x - 1]     + img[0][y][x + 1]    ) / 2.0f;
            } else if (green_even) {
                tmp[0][y][x] = (img[0][y][x - 1]     + img[0][y][x + 1]    ) / 2.0f;
                tmp[1][y][x] =  img[0][y][x];
                tmp[2][y][x] = (img[0][y - 1][x]     + img[0][y + 1][x]    ) / 2.0f;
            } else if (red) {
                tmp[0][y][x] =  img[0][y][x];
                tmp[1][y][x] = (img[0][y][x - 1]     + img[0][y][x + 1]    +
                                img[0][y - 1][x]     + img[0][y + 1][x]    ) / 4.0f;
                tmp[2][y][x] = (img[0][y - 1][x - 1] + img[0][y + 1][x - 1] +
                                img[0][y - 1][x + 1] + img[0][y + 1][x + 1]) / 4.0f;
            } else if (blue) {
                tmp[0][y][x] = (img[0][y - 1][x - 1] + img[0][y + 1][x - 1] +
                                img[0][y - 1][x + 1] + img[0][y + 1][x + 1]) / 4.0f;
                tmp[1][y][x] = (img[0][y][x - 1]     + img[0][y][x + 1]    +
                                img[0][y - 1][x]     + img[0][y + 1][x]    ) / 4.0f;
                tmp[2][y][x] =  img[0][y][x];
            }
        }
    }
    
    img = std::move(tmp);
    head.naxis3 = 3;
    head.naxis = 3;
}

/// MARK: - demosaic_x_trans

void demosaic_x_trans(ImageArray& img, Header& head) {
    const auto H = head.height;
    const auto W = head.width;
    auto tmp = make_rgb_buffer(H, W);
    
    for (auto y = 2; y <= H - 2; ++y) {
        for (auto x = 2; x <= W - 2; ++x) {
            const auto x2 = x - 1;
            const auto y2 = y - 1;
            const auto xpos = 1 + x2 - (x2 / 3) * 3;
            const auto ypos = 1 + y2 - (y2 / 3) * 3;
            const auto xpos6 = 1 + x2 - (x2 / 6) * 6;
            const auto ypos6 = 1 + y2 - (y2 / 6) * 6;
            
            auto red = 0.0f;
            auto blue = 0.0f;
            
            // Greens -- five orientations within the 3x3 cell
            if (xpos == 1 && ypos == 1) {
                red          = img[0][y + 1][x];
                tmp[1][y][x] = img[0][y][x];
                blue         = img[0][y][x + 1];
            } else if (xpos == 3 && ypos == 1) {
                red          = img[0][y + 1][x];
                tmp[1][y][x] = img[0][y][x];
                blue         = img[0][y][x - 1];
            } else if (xpos == 2 && ypos == 2) {
                red          = img[0][y][x + 1];
                tmp[1][y][x] = img[0][y][x];
                blue         = img[0][y + 1][x];
            } else if (xpos == 1 && ypos == 3) {
                red          = img[0][y - 1][x];
                tmp[1][y][x] = img[0][y][x];
                blue         = img[0][y][x + 1];
            } else if (xpos == 3 && ypos == 3) {
                red          = img[0][y - 1][x];
                tmp[1][y][x] = img[0][y][x];
                blue         = img[0][y][x - 1];
            } else if (xpos == 2 && ypos == 1) {
                red          = img[0][y - 1][x];
                tmp[1][y][x] = img[0][y][x + 1];
                blue         = img[0][y][x];
            } else if (xpos == 2 && ypos == 3) {
                red          = img[0][y + 1][x];
                tmp[1][y][x] = img[0][y][x + 1];
                blue         = img[0][y][x];
            } else if (xpos == 1 && ypos == 2) {
                red          = img[0][y][x];
                tmp[1][y][x] = img[0][y][x + 1];
                blue         = img[0][y][x - 1];
            } else if (xpos == 3 && ypos == 2) {
                red          = img[0][y][x];
                tmp[1][y][x] = img[0][y + 1][x];
                blue         = img[0][y][x + 1];
            }
            
            // Fix red/blue swap by 6x6 quadrant
            if (xpos6 <= 3 && ypos6 <= 3) {
                tmp[0][y][x] = red;
                tmp[2][y][x] = blue;
            } else if (xpos6 > 3 && ypos6 <= 3) {
                tmp[0][y][x] = blue;
                tmp[2][y][x] = red;
            } else if (xpos6 <= 3 && ypos6 > 3) {
                tmp[0][y][x] = blue;
                tmp[2][y][x] = red;
            } else if (xpos6 > 3 && ypos6 > 3) {
                tmp[0][y][x] = red;
                tmp[2][y][x] = blue;
            }
        }
    }
    
    img = std::move(tmp);
    head.naxis3 = 3;
    head.naxis = 3;
}

/// MARK: - demosaic_astrosimple

void demosaic_astrosimple(ImageArray& img, Header& head, int pattern) {
    auto offsetx = 0;
    auto offsety = 0;
    if (!pattern_offsets(pattern, offsetx, offsety)) {
        return;
    }
    
    const auto H = head.height;
    const auto W = head.width;
    auto tmp = make_rgb_buffer(H, W);
    
    for (auto y = 0; y <= H - 2; ++y) {
        for (auto x = 0; x <= W - 2; ++x) {
            const auto green_even = is_odd(x + 1 + offsetx) && is_odd(y + 1 + offsety);
            const auto green_odd  = is_odd(x + offsetx)     && is_odd(y + offsety);
            const auto red        = is_odd(x + offsetx)     && is_odd(y + 1 + offsety);
            const auto blue       = is_odd(x + 1 + offsetx) && is_odd(y + offsety);
            
            auto value = img[0][y][x];
            
            if (green_odd || green_even) {
                value /= 2.0f;
                tmp[1][y    ][x    ] += value;
                tmp[1][y + 1][x    ] += value;
                tmp[1][y    ][x + 1] += value;
                tmp[1][y + 1][x + 1] += value;
            } else if (red) {
                tmp[0][y    ][x    ] = value;
                tmp[0][y    ][x + 1] = value;
                tmp[0][y + 1][x    ] = value;
                tmp[0][y + 1][x + 1] = value;
            } else if (blue) {
                tmp[2][y    ][x    ] = value;
                tmp[2][y    ][x + 1] = value;
                tmp[2][y + 1][x    ] = value;
                tmp[2][y + 1][x + 1] = value;
            }
        }
    }
    
    img = std::move(tmp);
    head.naxis3 = 3;
    head.naxis = 3;
}

/// MARK: - demosaic_astrosimplebayercombined

void demosaic_astrosimplebayercombined(ImageArray& img, Header& head, int pattern) {
    auto offsetx = 0;
    auto offsety = 0;
    if (!pattern_offsets(pattern, offsetx, offsety)) {
        return;
    }
    
    const auto H = head.height;
    const auto W = head.width;
    auto tmp = make_rgb_buffer(H, W);
    
    for (auto y = 0; y <= H - 2; ++y) {
        for (auto x = 0; x <= W - 2; ++x) {
            const auto green_even = is_odd(x + 1 + offsetx) && is_odd(y + 1 + offsety);
            const auto green_odd  = is_odd(x + offsetx)     && is_odd(y + offsety);
            const auto red        = is_odd(x + offsetx)     && is_odd(y + 1 + offsety);
            const auto blue       = is_odd(x + 1 + offsetx) && is_odd(y + offsety);
            
            auto value = img[0][y][x];
            
            if (green_even) {
                // Bounds guard
                if (x == 0 || y == 0) {
                    continue;
                }
                value /= 2.0f;
                tmp[1][y    ][x    ] += value;
                tmp[1][y    ][x - 1] += value;
                tmp[1][y - 1][x    ] += value;
                tmp[1][y - 1][x - 1] += value;
            } else if (green_odd) {
                value /= 2.0f;
                tmp[1][y    ][x    ] += value;
                tmp[1][y    ][x + 1] += value;
                tmp[1][y + 1][x    ] += value;
                tmp[1][y + 1][x + 1] += value;
            } else if (red) {
                if (y == 0) {
                    continue;
                }
                tmp[0][y    ][x    ] = value;
                tmp[0][y    ][x + 1] = value;
                tmp[0][y - 1][x    ] = value;
                tmp[0][y - 1][x + 1] = value;
            } else if (blue) {
                if (x == 0) {
                    continue;
                }
                tmp[2][y    ][x    ] = value;
                tmp[2][y    ][x - 1] = value;
                tmp[2][y + 1][x    ] = value;
                tmp[2][y + 1][x - 1] = value;
            }
        }
    }
    
    img = std::move(tmp);
    head.naxis3 = 3;
    head.naxis = 3;
}

/// MARK: - demosaic_astroM_bilinear_interpolation

void demosaic_astroM_bilinear_interpolation(ImageArray& img, Header& head, int pattern) {
    auto offsetx = 0;
    auto offsety = 0;
    if (!pattern_offsets(pattern, offsetx, offsety)) {
        return;
    }
    
    const auto H = head.height;
    const auto W = head.width;
    auto tmp = make_rgb_buffer(H, W);
    
    // Mean background sampled on a coarse grid (replicates verbatim)
    auto count = 0;
    auto bg = 0.0f;
    {
        const auto yend = (H - 10) / 100;
        const auto xend = (W - 10) / 100;
        for (auto y = 10; y <= yend; ++y) {
            for (auto x = 10; x <= xend; ++x) {
                bg += img[0][y][x] +
                      img[0][y    ][x + 1] +
                      img[0][y + 1][x    ] +
                      img[0][y + 1][x + 1];
                count += 4;
            }
        }
    }
    
    // Gracefully handle zero-count (prevents divide-by-zero)
    if (count == 0) {
        bg = 0.0f;
    } else {
        bg /= static_cast<float>(count);
    }
    
    const auto signal  = 0.5f * bg;
    const auto signal2 = signal / 1.67f;
    
    for (auto y = 1; y <= H - 2; ++y) {
        for (auto x = 1; x <= W - 2; ++x) {
            const auto green_even = is_odd(x + 1 + offsetx) && is_odd(y + 1 + offsety);
            const auto green_odd  = is_odd(x + offsetx)     && is_odd(y + offsety);
            const auto red        = is_odd(x + offsetx)     && is_odd(y + 1 + offsety);
            const auto blue       = is_odd(x + 1 + offsetx) && is_odd(y + offsety);
            
            auto a1 = 0.0f, a2 = 0.0f, a3 = 0.0f, a4 = 0.0f;
            auto a5 = 0.0f, a6 = 0.0f, a7 = 0.0f, a8 = 0.0f;
            auto average1 = 0.0f, average2 = 0.0f, average3 = 0.0f, luminance = 0.0f;
            
            if (green_odd) {
                a1 = img[0][y - 1][x];
                a2 = img[0][y + 1][x];
                average1 = (a1 + a2) / 2.0f;
                average2 = img[0][y][x];
                a3 = img[0][y][x - 1];
                a4 = img[0][y][x + 1];
                average3 = (a3 + a4) / 2.0f;
                
                if (a1 > average1 + signal || a2 > average1 + signal ||
                    a3 > average2 + signal || a4 > average2 + signal) {
                    luminance = (average1 + average2 + average3) / 3.0f;
                    tmp[0][y][x] = luminance;
                    tmp[1][y][x] = luminance;
                    tmp[2][y][x] = luminance;
                } else {
                    tmp[0][y][x] = average1;
                    tmp[1][y][x] = average2;
                    tmp[2][y][x] = average3;
                }
            } else if (green_even) {
                a1 = img[0][y][x - 1];
                a2 = img[0][y][x + 1];
                average1 = (a1 + a2) / 2.0f;
                average2 = img[0][y][x];
                a3 = img[0][y - 1][x];
                a4 = img[0][y + 1][x];
                average3 = (a3 + a4) / 2.0f;
                
                if (a1 > average1 + signal || a2 > average1 + signal ||
                    a3 > average2 + signal || a4 > average2 + signal) {
                    luminance = (average1 + average2 + average3) / 3.0f;
                    tmp[0][y][x] = luminance;
                    tmp[1][y][x] = luminance;
                    tmp[2][y][x] = luminance;
                } else {
                    tmp[0][y][x] = average1;
                    tmp[1][y][x] = average2;
                    tmp[2][y][x] = average3;
                }
            } else if (red) {
                average1 = img[0][y][x];
                
                a1 = img[0][y][x - 1];
                a2 = img[0][y][x + 1];
                a3 = img[0][y - 1][x];
                a4 = img[0][y + 1][x];
                average2 = (a1 + a2 + a3 + a4) / 4.0f;
                
                a5 = img[0][y - 1][x - 1];
                a6 = img[0][y + 1][x - 1];
                a7 = img[0][y - 1][x + 1];
                a8 = img[0][y + 1][x + 1];
                average3 = (a5 + a6 + a7 + a8) / 4.0f;
                
                if (a1 > average2 + signal2 || a2 > average2 + signal2 ||
                    a3 > average2 + signal2 || a4 > average2 + signal2 ||
                    a5 > average3 + signal2 || a6 > average3 + signal2 ||
                    a7 > average3 + signal2 || a8 > average3 + signal2) {
                    luminance = (average1 + average2 + average3) / 3.0f;
                    tmp[0][y][x] = luminance;
                    tmp[1][y][x] = luminance;
                    tmp[2][y][x] = luminance;
                } else {
                    tmp[0][y][x] = average1;
                    tmp[1][y][x] = average2;
                    tmp[2][y][x] = average3;
                }
            } else if (blue) {
                // average1 is overwritten below (matches original dead code)
                average1 = img[0][y][x];
                
                a1 = img[0][y - 1][x - 1];
                a2 = img[0][y + 1][x - 1];
                a3 = img[0][y - 1][x + 1];
                a4 = img[0][y + 1][x + 1];
                average1 = (a1 + a2 + a3 + a4) / 4.0f;
                
                a5 = img[0][y][x - 1];
                a6 = img[0][y][x + 1];
                a7 = img[0][y - 1][x];
                a8 = img[0][y + 1][x];
                average2 = (a5 + a6 + a7 + a8) / 4.0f;
                
                average3 = img[0][y][x];
                
                if (a1 > average1 + signal2 || a2 > average1 + signal2 ||
                    a3 > average1 + signal2 || a4 > average1 + signal2 ||
                    a5 > average2 + signal2 || a6 > average2 + signal2 ||
                    a7 > average2 + signal2 || a8 > average2 + signal2) {
                    luminance = (average1 + average2 + average3) / 3.0f;
                    tmp[0][y][x] = luminance;
                    tmp[1][y][x] = luminance;
                    tmp[2][y][x] = luminance;
                } else {
                    tmp[0][y][x] = average1;
                    tmp[1][y][x] = average2;
                    tmp[2][y][x] = average3;
                }
            }
        }
    }
    
    img = std::move(tmp);
    head.naxis3 = 3;
    head.naxis = 3;
}

/// MARK: - demosaic_astroC_bilinear_interpolation

void demosaic_astroC_bilinear_interpolation(ImageArray& img, Header& head,
                                            int saturation, int pattern) {
    auto offsetx = 0;
    auto offsety = 0;
    if (!pattern_offsets(pattern, offsetx, offsety)) {
        return;
    }
    
    constexpr auto step = 5;
    const auto H = head.height;
    const auto W = head.width;
    const auto sat = static_cast<float>(saturation);
    
    auto tmp = make_rgb_buffer(H, W);
    
    auto bg = 0.0;
    auto counter = 0;
    
    for (auto y = 1; y <= H - 2; ++y) {
        for (auto x = 1; x <= W - 2; ++x) {
            const auto green_even = is_odd(x + 1 + offsetx) && is_odd(y + 1 + offsety);
            const auto green_odd  = is_odd(x + offsetx)     && is_odd(y + offsety);
            const auto red        = is_odd(x + offsetx)     && is_odd(y + 1 + offsety);
            const auto blue       = is_odd(x + 1 + offsetx) && is_odd(y + offsety);
            
            auto a1 = 0.0f, a2 = 0.0f, a3 = 0.0f, a4 = 0.0f;
            auto a5 = 0.0f, a6 = 0.0f, a7 = 0.0f, a8 = 0.0f;
            auto average1 = 0.0f, average2 = 0.0f, average3 = 0.0f;
            
            if (green_odd) {
                a1 = img[0][y - 1][x];
                a2 = img[0][y + 1][x];
                average1 = (a1 + a2) / 2.0f;
                average2 = img[0][y][x];
                a3 = img[0][y][x - 1];
                a4 = img[0][y][x + 1];
                average3 = (a3 + a4) / 2.0f;
                
                if (a1 > sat || a2 > sat || a3 > sat || a4 > sat) {
                    tmp[0][y][x] = (average1 + average2 + average3) / 3.0f;
                    tmp[1][y][x] = kSaturatedMarker;
                } else {
                    tmp[0][y][x] = average1;
                    tmp[1][y][x] = average2;
                    tmp[2][y][x] = average3;
                }
            } else if (green_even) {
                a1 = img[0][y][x - 1];
                a2 = img[0][y][x + 1];
                average1 = (a1 + a2) / 2.0f;
                average2 = img[0][y][x];
                a3 = img[0][y - 1][x];
                a4 = img[0][y + 1][x];
                average3 = (a3 + a4) / 2.0f;
                
                if (a1 > sat || a2 > sat || a3 > sat || a4 > sat) {
                    tmp[0][y][x] = (average1 + average2 + average3) / 3.0f;
                    tmp[1][y][x] = kSaturatedMarker;
                } else {
                    tmp[0][y][x] = average1;
                    tmp[1][y][x] = average2;
                    tmp[2][y][x] = average3;
                }
            } else if (red) {
                average1 = img[0][y][x];
                
                a1 = img[0][y][x - 1];
                a2 = img[0][y][x + 1];
                a3 = img[0][y - 1][x];
                a4 = img[0][y + 1][x];
                average2 = (a1 + a2 + a3 + a4) / 4.0f;
                
                a5 = img[0][y - 1][x - 1];
                a6 = img[0][y + 1][x - 1];
                a7 = img[0][y - 1][x + 1];
                a8 = img[0][y + 1][x + 1];
                average3 = (a5 + a6 + a7 + a8) / 4.0f;
                
                if (a1 > sat || a2 > sat || a3 > sat || a4 > sat ||
                    a5 > sat || a6 > sat || a7 > sat || a8 > sat) {
                    tmp[0][y][x] = (average1 + average2 + average3) / 3.0f;
                    tmp[1][y][x] = kSaturatedMarker;
                } else {
                    tmp[0][y][x] = average1;
                    tmp[1][y][x] = average2;
                    tmp[2][y][x] = average3;
                    bg += static_cast<double>(average1) +
                          static_cast<double>(average2) +
                          static_cast<double>(average3);
                    counter += 3;
                }
            } else if (blue) {
                average1 = img[0][y][x];
                
                a1 = img[0][y - 1][x - 1];
                a2 = img[0][y + 1][x - 1];
                a3 = img[0][y - 1][x + 1];
                a4 = img[0][y + 1][x + 1];
                average1 = (a1 + a2 + a3 + a4) / 4.0f;
                
                a5 = img[0][y][x - 1];
                a6 = img[0][y][x + 1];
                a7 = img[0][y - 1][x];
                a8 = img[0][y + 1][x];
                average2 = (a5 + a6 + a7 + a8) / 4.0f;
                
                average3 = img[0][y][x];
                
                if (a1 > sat || a2 > sat || a3 > sat || a4 > sat ||
                    a5 > sat || a6 > sat || a7 > sat || a8 > sat) {
                    tmp[0][y][x] = (average1 + average2 + average3) / 3.0f;
                    tmp[1][y][x] = kSaturatedMarker;
                } else {
                    tmp[0][y][x] = average1;
                    tmp[1][y][x] = average2;
                    tmp[2][y][x] = average3;
                }
            }
        }
    }
    
    // Move the demosaiced buffer into img, then run the saturation fix-up
    // sweep directly on img (matches original two-pass design)
    img = std::move(tmp);
    
    [[maybe_unused]] auto sat_counter = 0;
    
    if (counter > 0) {
        const auto bg_avg = bg / static_cast<double>(counter);
        const auto bgf = static_cast<float>(bg_avg);
        
        for (auto fitsY = 0; fitsY < H; ++fitsY) {
            for (auto fitsX = 0; fitsX < W; ++fitsX) {
                if (img[1][fitsY][fitsX] != kSaturatedMarker) {
                    continue;
                }
                
                auto colred = 0.0f;
                auto colgreen = 0.0f;
                auto colblue = 0.0f;
                auto local_count = 0;
                ++sat_counter;
                auto luminance = img[0][fitsY][fitsX];
                luminance -= bgf;
                
                for (auto dy = -step; dy <= step; ++dy) {
                    for (auto dx = -step; dx <= step; ++dx) {
                        const auto x2 = fitsX + dx;
                        const auto y2 = fitsY + dy;
                        if (x2 < 0 || x2 >= W || y2 < 0 || y2 >= H) {
                            continue;
                        }
                        const auto sqr_dist = static_cast<double>(dx) * dx +
                                              static_cast<double>(dy) * dy;
                        if (sqr_dist > static_cast<double>(step) * step) {
                            continue;
                        }
                        
                        const auto g = img[1][y2][x2];
                        if (g == kSaturatedMarker) {
                            continue;
                        }
                        const auto r = img[0][y2][x2];
                        const auto b = img[2][y2][x2];
                        
                        if ((r - bgf) > 0.0f) {
                            colred += (r - bgf);
                        }
                        if ((g - bgf) > 0.0f) {
                            colgreen += (g - bgf);
                        }
                        if ((b - bgf) > 0.0f) {
                            colblue += (b - bgf);
                        }
                        ++local_count;
                    }
                }
                
                if (local_count >= 1) {
                    colred   /= static_cast<float>(local_count);
                    colgreen /= static_cast<float>(local_count);
                    colblue  /= static_cast<float>(local_count);
                    
                    // Prevent purple stars
                    const auto lowest = (colred > colblue) ? colblue : colred;
                    if (colgreen < lowest) {
                        colgreen = lowest;
                    }
                    
                    const auto rgb = (colred + colgreen + colblue + 0.00001f) / 3.0f;
                    img[0][fitsY][fitsX] = bgf + luminance * colred   / rgb;
                    img[1][fitsY][fitsX] = bgf + luminance * colgreen / rgb;
                    img[2][fitsY][fitsX] = bgf + luminance * colblue  / rgb;
                } else {
                    // No unsaturated neighbours found; replicate luminance
                    img[1][fitsY][fitsX] = img[0][fitsY][fitsX];
                    img[2][fitsY][fitsX] = img[0][fitsY][fitsX];
                }
            }
        }
    } else {
        // Fully saturated image -- fill all channels
        const auto satf = static_cast<float>(saturation);
        for (auto fitsY = 0; fitsY < H; ++fitsY) {
            for (auto fitsX = 0; fitsX < W; ++fitsX) {
                img[0][fitsY][fitsX] = satf;
                img[1][fitsY][fitsX] = satf;
                img[2][fitsY][fitsX] = satf;
            }
        }
    }
    
    head.naxis3 = 3;
    head.naxis = 3;
}

/// MARK: - demosaic_superpixel

void demosaic_superpixel(ImageArray& img, Header& head, int pattern) {
    const auto w = head.width / 2;
    const auto h = head.height / 2;
    auto tmp = make_rgb_buffer(h, w);
    
    auto fill = [&](auto kernel) {
        for (auto y = 0; y < h; ++y) {
            for (auto x = 0; x < w; ++x) {
                kernel(x, y);
            }
        }
    };
    
    if (pattern == 0) {  // GRBG
        fill([&](int x, int y) {
            const auto x2 = x + x;
            const auto y2 = y + y;
            tmp[0][y][x] =  img[0][y2    ][x2 + 1];
            tmp[1][y][x] = (img[0][y2    ][x2    ] + img[0][y2 + 1][x2 + 1]) / 2.0f;
            tmp[2][y][x] =  img[0][y2 + 1][x2    ];
        });
    } else if (pattern == 1) {  // BGGR
        fill([&](int x, int y) {
            const auto x2 = x + x;
            const auto y2 = y + y;
            tmp[0][y][x] =  img[0][y2 + 1][x2 + 1];
            tmp[1][y][x] = (img[0][y2    ][x2 + 1] + img[0][y2 + 1][x2    ]) / 2.0f;
            tmp[2][y][x] =  img[0][y2    ][x2    ];
        });
    } else if (pattern == 2) {  // RGGB
        fill([&](int x, int y) {
            const auto x2 = x + x;
            const auto y2 = y + y;
            tmp[0][y][x] =  img[0][y2    ][x2    ];
            tmp[1][y][x] = (img[0][y2    ][x2 + 1] + img[0][y2 + 1][x2    ]) / 2.0f;
            tmp[2][y][x] =  img[0][y2 + 1][x2 + 1];
        });
    } else if (pattern == 3) {  // GBRG
        fill([&](int x, int y) {
            const auto x2 = x + x;
            const auto y2 = y + y;
            tmp[0][y][x] =  img[0][y2 + 1][x2    ];
            tmp[1][y][x] = (img[0][y2    ][x2    ] + img[0][y2 + 1][x2 + 1]) / 2.0f;
            tmp[2][y][x] =  img[0][y2    ][x2 + 1];
        });
    }
    
    img = std::move(tmp);
    head.width = w;
    head.height = h;
    head.naxis3 = 3;
    head.naxis = 3;
}

/// MARK: - preserve_colour_saturated_bayer

void preserve_colour_saturated_bayer(ImageArray& img, const Header& head) {
    const auto w = head.width / 2;
    const auto h = head.height / 2;
    
    for (auto fitsY = 0; fitsY < h; ++fitsY) {
        for (auto fitsX = 1; fitsX < w; ++fitsX) {
            const auto yy0 = fitsY * 2;
            const auto yy1 = fitsY * 2 + 1;
            const auto xx0 = fitsX * 2;
            const auto xx1 = fitsX * 2 + 1;
            const auto px0 = (fitsX - 1) * 2;
            const auto px1 = (fitsX - 1) * 2 + 1;
            if (img[0][yy0][xx0] > 65500.0f ||
                img[0][yy0][xx1] > 65500.0f ||
                img[0][yy1][xx0] > 65500.0f ||
                img[0][yy1][xx1] > 65500.0f) {
                img[0][yy0][xx0] = img[0][yy0][px0];
                img[0][yy0][xx1] = img[0][yy0][px1];
                img[0][yy1][xx0] = img[0][yy1][px0];
                img[0][yy1][xx1] = img[0][yy1][px1];
            }
        }
    }
}

/// MARK: - get_demosaic_pattern

int get_demosaic_pattern(int pattern_hint,
                         double xbayroff,
                         double ybayroff,
                         const std::string& roworder) {
    auto result = (pattern_hint < 0 || pattern_hint > 4) ? 2 : pattern_hint;
    
    // X-Trans is not affected by the Bayer-offset parity flips
    if (result == 4) {
        return result;
    }
    
    // ROWORDER correction
    auto ybayroff2 = ybayroff;
    if (roworder.find("BOT") != std::string::npos) {
        ybayroff2 += 1.0;
    }
    
    if (is_odd(static_cast<int>(std::lround(xbayroff)))) {
        switch (result) {
            case 2: result = 0; break;
            case 0: result = 2; break;
            case 1: result = 3; break;
            case 3: result = 1; break;
            default: break;
        }
    }
    
    if (is_odd(static_cast<int>(std::lround(ybayroff2)))) {
        switch (result) {
            case 1: result = 0; break;
            case 0: result = 1; break;
            case 3: result = 2; break;
            case 2: result = 3; break;
            default: break;
        }
    }
    
    return result;
}

/// MARK: - demosaic_bayer

void demosaic_bayer(ImageArray& img, Header& head, int pattern,
                    DemosaicMethod method) {
    if (pattern == 4) {
        demosaic_x_trans(img, head);
        return;
    }
    
    switch (method) {
        case DemosaicMethod::AstroC: {
            auto sat = 0;
            if (head.datamax_org > 16384.0) {
                sat = 65535 / 2;  // 16-bit
            } else if (head.datamax_org > 4096.0) {
                sat = 16383 / 2;  // 14-bit
            } else {
                sat = 4095  / 2;  // 12-bit
            }
            demosaic_astroC_bilinear_interpolation(img, head, sat, pattern);
            break;
        }
        case DemosaicMethod::Simple:
            demosaic_astrosimple(img, head, pattern);
            break;
        case DemosaicMethod::AstroM:
            demosaic_astroM_bilinear_interpolation(img, head, pattern);
            break;
        case DemosaicMethod::Superpixel:
            demosaic_superpixel(img, head, pattern);
            break;
        case DemosaicMethod::Bilinear:
        default:
            demosaic_bilinear_interpolation(img, head, pattern);
            break;
    }
}

/// MARK: - demosaic_advanced

void demosaic_advanced(ImageArray& img, Header& head, int pattern,
                       DemosaicMethod method) {
    // Dispatches to the selected algorithm
    demosaic_bayer(img, head, pattern, method);
    
    // TODO(astap++ port): once the UI/auto-level/colour-smooth pipeline is
    // translated, restore the post-demosaic steps (auto background level,
    // apply_factors, smart_colour_smooth, histogram refresh).
}
    
} // namespace
