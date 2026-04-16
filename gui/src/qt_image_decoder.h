///----------------------------------------
///      @file qt_image_decoder.h
///   @ingroup ASTAP++
///     @brief Qt-backed @ref astap::core::IImageDecoder for the GUI build.
///   @details Decodes TIFF / PNG / JPEG / BMP (and anything else Qt's image
///            plugins recognise) via @c QImageReader. No extra third-party
///            deps beyond Qt itself. RAW camera files still require LibRaw
///            or DCRAW and are not handled here.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include "../../src/core/image_io.h"

#include <memory>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @brief Construct a Qt-backed image decoder and install it as the active
///        @ref astap::core::IImageDecoder.
/// @return Strong reference to the decoder so the caller can keep it alive
///         or swap it out later.
///----------------------------------------

[[nodiscard]] std::shared_ptr<astap::core::IImageDecoder> install_qt_image_decoder();

} // namespace astap::gui
