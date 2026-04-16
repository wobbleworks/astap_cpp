#pragma once

// Ready-to-use IImageDecoder backed by stb_image.h.
//
// stb_image.h is a single-header public-domain image decoder supporting
// JPEG, PNG, BMP, TGA, PSD, GIF, HDR, PIC, and PNM. This adapter is header-
// only; it builds and installs as the active decoder when you call
// install_stb_image_decoder() from main().
//
// Usage:
//   1. Drop stb_image.h into the include path (https://github.com/nothings/stb).
//   2. In ONE translation unit (typically cli/stb_image_decoder.cpp) define
//      STB_IMAGE_IMPLEMENTATION before including stb_image.h.
//   3. Call astap::cli::install_stb_image_decoder() early in main().
//
// TIFF is not covered by stb_image. If you need TIFF, combine this adapter
// with a libtiff-backed IImageDecoder that delegates .TIF/.TIFF to libtiff
// and everything else to this one (see IImageDecoder interface for the
// composition pattern).

#include <memory>

#include "../src/core/image_io.h"

namespace astap::cli {

// Creates and installs the stb_image-backed decoder as the active
// IImageDecoder via astap::core::set_image_decoder(). Returns a strong
// reference so callers can re-install later if they want to swap it out.
std::shared_ptr<astap::core::IImageDecoder> install_stb_image_decoder();

}  // namespace astap::cli
