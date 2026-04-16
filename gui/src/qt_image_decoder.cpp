///----------------------------------------
///      @file qt_image_decoder.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the Qt-backed image decoder.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "qt_image_decoder.h"

#include <QImage>
#include <QImageReader>
#include <QString>

#include <cstdint>
#include <vector>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

namespace {

///----------------------------------------
/// @class QtImageDecoder
/// @brief @ref astap::core::IImageDecoder implementation backed by QImageReader.
/// @details Stores raw sample values (0..255 for 8-bit, 0..65535 for 16-bit)
///          per the contract documented on @ref astap::core::IImageDecoder.
///          The upstream loader rescales 8-bit samples by 256 to fit the
///          0..65535 convention.
///----------------------------------------

class QtImageDecoder final : public astap::core::IImageDecoder {

public:
	bool decode_raster(const std::filesystem::path& path,
	                   std::string_view             /*ext_upper*/,
	                   DecodedImage&                out,
	                   std::string&                 error_out) override {

		// Decode the file via Qt's image plugins
		QImageReader reader(QString::fromStdString(path.string()));
		reader.setAutoTransform(true);
		auto img = reader.read();
		if (img.isNull()) {
			error_out = reader.errorString().toStdString();
			return false;
		}

		// Decide whether to preserve >8-bit samples
		const auto wantsHi = img.depth() > 32
			|| img.format() == QImage::Format_Grayscale16
			|| img.format() == QImage::Format_RGBA64
			|| img.format() == QImage::Format_RGBA64_Premultiplied;

		const auto isMono = img.isGrayscale();
		const auto channels = isMono ? 1 : 3;

		if (wantsHi) {
			return convertHi(img, isMono, channels, out, error_out);
		}
		return convertLo(img, isMono, channels, out, error_out);
	}

	bool decode_raw(const std::filesystem::path& /*path*/,
	                bool                         /*save_intermediate*/,
	                DecodedImage&                /*out*/,
	                std::string&                 error_out) override {
		// Not supported by Qt; require a separate LibRaw/DCRAW adapter.
		error_out = "Qt-backed decoder does not support camera RAW files.";
		return false;
	}

private:
	[[nodiscard]] static bool convertLo(QImage&          img,
	                                    bool             isMono,
	                                    int              channels,
	                                    DecodedImage&    out,
	                                    std::string&     error_out) {

		// Force a known 8-bit layout so scanLine has a fixed stride
		const auto target = isMono ? QImage::Format_Grayscale8 : QImage::Format_RGB888;
		if (img.format() != target) {
			img = img.convertToFormat(target);
		}
		if (img.isNull()) {
			error_out = "QImage::convertToFormat failed (8-bit)";
			return false;
		}

		const auto w = img.width();
		const auto h = img.height();
		out.bits_per_sample = 8;
		out.width = w;
		out.height = h;
		out.channels = channels;
		out.pixels.assign(channels,
			std::vector<std::vector<float>>(h, std::vector<float>(w, 0.0f)));

		// Pull samples scanline by scanline
		for (int y = 0; y < h; ++y) {
			const auto* row = img.constScanLine(y);
			if (isMono) {
				auto& dst = out.pixels[0][y];
				for (int x = 0; x < w; ++x) {
					dst[x] = static_cast<float>(row[x]);
				}
			} else {
				for (int x = 0; x < w; ++x) {
					out.pixels[0][y][x] = static_cast<float>(row[3 * x + 0]);
					out.pixels[1][y][x] = static_cast<float>(row[3 * x + 1]);
					out.pixels[2][y][x] = static_cast<float>(row[3 * x + 2]);
				}
			}
		}
		return true;
	}

	[[nodiscard]] static bool convertHi(QImage&          img,
	                                    bool             isMono,
	                                    int              channels,
	                                    DecodedImage&    out,
	                                    std::string&     error_out) {

		// Force a known 16-bit layout so scanLine has a fixed stride
		const auto target = isMono ? QImage::Format_Grayscale16 : QImage::Format_RGBA64;
		if (img.format() != target) {
			img = img.convertToFormat(target);
		}
		if (img.isNull()) {
			error_out = "QImage::convertToFormat failed (16-bit)";
			return false;
		}

		const auto w = img.width();
		const auto h = img.height();
		out.bits_per_sample = 16;
		out.width = w;
		out.height = h;
		out.channels = channels;
		out.pixels.assign(channels,
			std::vector<std::vector<float>>(h, std::vector<float>(w, 0.0f)));

		// Pull samples scanline by scanline (16 bits per channel)
		for (int y = 0; y < h; ++y) {
			const auto* row = reinterpret_cast<const std::uint16_t*>(img.constScanLine(y));
			if (isMono) {
				auto& dst = out.pixels[0][y];
				for (int x = 0; x < w; ++x) {
					dst[x] = static_cast<float>(row[x]);
				}
			} else {
				// RGBA64 = 4 × uint16 per pixel; ignore alpha
				for (int x = 0; x < w; ++x) {
					out.pixels[0][y][x] = static_cast<float>(row[4 * x + 0]);
					out.pixels[1][y][x] = static_cast<float>(row[4 * x + 1]);
					out.pixels[2][y][x] = static_cast<float>(row[4 * x + 2]);
				}
			}
		}
		return true;
	}
};

} // namespace

std::shared_ptr<astap::core::IImageDecoder> install_qt_image_decoder() {
	auto decoder = std::make_shared<QtImageDecoder>();
	astap::core::set_image_decoder(decoder);
	return decoder;
}

} // namespace astap::gui
