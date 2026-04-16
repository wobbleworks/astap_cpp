///----------------------------------------
///      @file image_viewer.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the central image canvas.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "image_viewer.h"

#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>
#include <QResizeEvent>
#include <QWheelEvent>

#include "../../src/core/photometry.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

namespace {

constexpr double kZoomStep = 1.25;
constexpr double kMinZoom = 0.02;
constexpr double kMaxZoom = 64.0;

// Pixel-value range that the working data is normalised to (load_image
// stores 8-bit sources scaled ×256 and 16-bit / float sources renormalised
// to 0..65535; floats from -32 land in 0..1 and are scaled at decode time).
constexpr float kFullRange = 65535.0f;

// log1p denominator used for the logarithmic stretch path.
constexpr float kLogDenom = 6.5512f;  // log1p(700) — matches a 700:1 stretch knee.

[[nodiscard]] constexpr std::uint8_t toByte(float t) noexcept {
	if (t <= 0.0f) {
		return 0;
	}
	if (t >= 1.0f) {
		return 255;
	}
	return static_cast<std::uint8_t>(t * 255.0f + 0.5f);
}

} // namespace

///----------------------------------------
/// MARK: ImageViewer
///----------------------------------------

ImageViewer::ImageViewer(QWidget* parent) :
	QWidget(parent) {

	// Solid background, opaque for fast paint
	setAttribute(Qt::WA_OpaquePaintEvent);
	setAutoFillBackground(false);
	setMouseTracking(true);
	setFocusPolicy(Qt::StrongFocus);

	// Coalesce rapid stretch / flip changes (e.g. slider drags) into a
	// single render. Without this, a full image re-render on every slider
	// tick blocks the GUI thread long enough that the slider thumb can't
	// repaint cleanly and trails artifacts.
	_renderTimer.setSingleShot(true);
	_renderTimer.setInterval(40);
	connect(&_renderTimer, &QTimer::timeout, this, [this]() {
		renderImage();
		update();
	});
}

void ImageViewer::setImage(astap::ImageArray image, astap::Header header) {
	// A pending coalesced render would just redo what we're about to do.
	_renderTimer.stop();

	// Markers from a previous image no longer apply.
	_stars.clear();
	_annotations.clear();

	// Take ownership
	_image = std::move(image);
	_header = header;

	// Compute statistics + histogram + reset stretch defaults
	computeStatistics();

	// Render to QImage cache and reset view
	renderImage();
	fitToWindow();

	emit imageLoaded();
	emit stretchChanged();
	emit flipChanged();
}

void ImageViewer::setStarMarkers(std::vector<DetectedStar> stars) {
	_stars = std::move(stars);
	update();
}

void ImageViewer::clearStarMarkers() {
	_stars.clear();
	update();
}

void ImageViewer::setAnnotations(std::vector<AnnotationMarker> annotations) {
	_annotations = std::move(annotations);
	update();
}

void ImageViewer::clearAnnotations() {
	_annotations.clear();
	update();
}

void ImageViewer::clear() {
	// Drop everything
	_image.clear();
	_header = {};
	_rendered = QImage();
	_histogram.clear();
	_stats.clear();
	_stars.clear();
	_sampleMin = _sampleMax = 0.0f;
	_dataLo = 0.0f;
	_dataHi = kFullRange;
	_stretchLo = 0.0f;
	_stretchHi = kFullRange;
	_zoom = 1.0;
	_fitMode = true;
	_pan = {0.0, 0.0};
	update();
	emit imageLoaded();
}

void ImageViewer::setStretchRange(float lo, float hi) {
	// Reject invalid order (caller's responsibility, but be defensive)
	if (hi <= lo) {
		hi = lo + 1.0f;
	}
	if (lo == _stretchLo && hi == _stretchHi) {
		return;
	}
	_stretchLo = lo;
	_stretchHi = hi;
	scheduleRender();
	emit stretchChanged();
}

void ImageViewer::setLogarithmicStretch(bool on) {
	if (on == _logStretch) {
		return;
	}
	_logStretch = on;
	scheduleRender();
	emit stretchChanged();
}

void ImageViewer::setFlipHorizontal(bool on) {
	if (on == _flipH) {
		return;
	}
	_flipH = on;
	scheduleRender();
	emit flipChanged();
}

void ImageViewer::setFlipVertical(bool on) {
	if (on == _flipV) {
		return;
	}
	_flipV = on;
	scheduleRender();
	emit flipChanged();
}

void ImageViewer::autoStretch() {
	// Re-apply the load-time defaults
	setStretchRange(_dataLo, _dataHi);
}

void ImageViewer::zoomIn() {
	// Zoom around viewport centre
	const auto centre = QPointF(width() * 0.5, height() * 0.5);
	applyZoomAtViewport(_zoom * kZoomStep, centre);
}

void ImageViewer::zoomOut() {
	// Zoom around viewport centre
	const auto centre = QPointF(width() * 0.5, height() * 0.5);
	applyZoomAtViewport(_zoom / kZoomStep, centre);
}

void ImageViewer::fitToWindow() {
	// Restore fit-to-window mode (the actual zoom is computed in paintEvent)
	_fitMode = true;
	_pan = {0.0, 0.0};
	update();
}

void ImageViewer::actualSize() {
	// 1:1 pixel-for-pixel
	_fitMode = false;
	_zoom = 1.0;
	_pan = {0.0, 0.0};
	update();
}

void ImageViewer::paintEvent(QPaintEvent* event) {
	QPainter painter(this);

	// Background
	painter.fillRect(event->rect(), Qt::black);

	if (_rendered.isNull()) {
		// No image — show a hint
		painter.setPen(Qt::lightGray);
		painter.drawText(rect(), Qt::AlignCenter, tr("No image loaded"));
		return;
	}

	// Compute the effective zoom
	auto effectiveZoom = _zoom;
	if (_fitMode) {
		const auto sx = static_cast<double>(width()) / _rendered.width();
		const auto sy = static_cast<double>(height()) / _rendered.height();
		effectiveZoom = std::min(sx, sy);
	}

	// Compute destination rect, centred (plus pan offset)
	const auto dstW = _rendered.width() * effectiveZoom;
	const auto dstH = _rendered.height() * effectiveZoom;
	const auto cx = (width() - dstW) * 0.5 + _pan.x();
	const auto cy = (height() - dstH) * 0.5 + _pan.y();
	const auto dst = QRectF(cx, cy, dstW, dstH);

	// Smooth interpolation only when zoomed out; nearest when zoomed in,
	// so the user can see actual pixels.
	painter.setRenderHint(QPainter::SmoothPixmapTransform, effectiveZoom < 1.0);
	painter.drawImage(dst, _rendered);

	// Star overlays — draw circles at detected star positions with radius
	// proportional to HFD, scaled to the current zoom so markers track the
	// image as the user pans / zooms.
	if (!_stars.empty() && _header.width > 0 && _header.height > 0) {
		painter.setRenderHint(QPainter::Antialiasing, true);
		painter.setBrush(Qt::NoBrush);
		QPen pen(QColor(80, 255, 120, 220));
		pen.setWidthF(1.2);
		painter.setPen(pen);

		const auto w = _header.width;
		const auto h = _header.height;

		for (const auto& star : _stars) {
			// star.x, star.y are 1-based FITS coords (y up from bottom).
			// Reverse those to image-row / image-column, apply user flips.
			const auto imgCol = star.x - 1.0;
			const auto imgRow = (h - 1) - (star.y - 1.0);
			const auto srcX = _flipH ? (w - 1 - imgCol) : imgCol;
			const auto srcY = _flipV ? ((h - 1) - imgRow) : imgRow;

			// Map to viewport pixels using the same transform as the image.
			const auto vx = dst.left() + (srcX + 0.5) * effectiveZoom;
			const auto vy = dst.top() + (srcY + 0.5) * effectiveZoom;
			const auto rPx = std::max(3.0, 1.5 * star.hfd * effectiveZoom);
			painter.drawEllipse(QPointF(vx, vy), rPx, rPx);
		}
	}

	// Deep-sky annotation overlay — names + ellipses for cataloged objects.
	if (!_annotations.empty() && _header.width > 0 && _header.height > 0) {
		painter.setRenderHint(QPainter::Antialiasing, true);

		QPen ellipsePen(QColor(255, 200, 60, 200));
		ellipsePen.setWidthF(1.0);
		QFont labelFont = painter.font();
		labelFont.setPointSizeF(std::max(7.0, 9.0 * std::min(1.0, effectiveZoom)));
		painter.setFont(labelFont);

		const auto w = _header.width;
		const auto h = _header.height;

		for (const auto& ann : _annotations) {
			const auto imgCol = ann.x - 1.0;
			const auto imgRow = (h - 1) - (ann.y - 1.0);
			const auto srcX = _flipH ? (w - 1 - imgCol) : imgCol;
			const auto srcY = _flipV ? ((h - 1) - imgRow) : imgRow;

			const auto vx = dst.left() + (srcX + 0.5) * effectiveZoom;
			const auto vy = dst.top() + (srcY + 0.5) * effectiveZoom;

			// Ellipse (if the object has a measurable size)
			const auto majorPx = ann.majorPx * effectiveZoom;
			const auto minorPx = ann.minorPx * effectiveZoom;
			if (majorPx > 3.0) {
				painter.setPen(ellipsePen);
				painter.setBrush(Qt::NoBrush);
				painter.save();
				painter.translate(vx, vy);
				painter.rotate(-ann.paDeg);
				painter.drawEllipse(QPointF(0, 0), majorPx, minorPx);
				painter.restore();
			}

			// Label — offset slightly above the ellipse
			painter.setPen(QColor(255, 200, 60, 230));
			const auto labelOffset = std::max(8.0, majorPx + 4.0);
			painter.drawText(QPointF(vx + 4, vy - labelOffset), ann.name);
		}
	}
}

void ImageViewer::wheelEvent(QWheelEvent* event) {
	if (_rendered.isNull()) {
		event->ignore();
		return;
	}

	// Each notch is 120 units (Qt convention).
	const auto delta = event->angleDelta().y();
	if (delta == 0) {
		event->ignore();
		return;
	}

	// Take the current effective zoom as starting point, then leave fit mode.
	auto effectiveZoom = _zoom;
	if (_fitMode) {
		const auto sx = static_cast<double>(width()) / _rendered.width();
		const auto sy = static_cast<double>(height()) / _rendered.height();
		effectiveZoom = std::min(sx, sy);
		_fitMode = false;
		_zoom = effectiveZoom;
	}

	// Step proportionally to the wheel delta
	const auto steps = delta / 120.0;
	const auto factor = std::pow(kZoomStep, steps);
	applyZoomAtViewport(_zoom * factor, event->position());
	event->accept();
}

void ImageViewer::mousePressEvent(QMouseEvent* event) {
	if (event->button() == Qt::LeftButton && hasImage()) {
		// Begin pan
		_panning = true;
		_panAnchor = event->pos();
		_panAnchorOffset = _pan;
		setCursor(Qt::ClosedHandCursor);
	}
}

void ImageViewer::mouseMoveEvent(QMouseEvent* event) {
	if (_panning) {
		// Translate the image by the cursor delta
		_pan = _panAnchorOffset + (event->pos() - _panAnchor);
		_fitMode = false;  // any pan exits fit mode
		clampPan();
		update();
	}

	// Always advertise the cursor position so the host can show pixel /
	// celestial coordinates.
	bool inImage = false;
	const auto imgPos = viewportToImage(event->position(), inImage);
	emit cursorMoved(imgPos, inImage);
}

void ImageViewer::mouseReleaseEvent(QMouseEvent* event) {
	if (event->button() == Qt::LeftButton && _panning) {
		// End pan
		_panning = false;
		unsetCursor();
	}
}

void ImageViewer::resizeEvent(QResizeEvent* /*event*/) {
	// In fit mode the painter recomputes scale; in zoom mode keep pan in range.
	if (!_fitMode) {
		clampPan();
	}
}

void ImageViewer::scheduleRender() {
	// Single-shot: rapid back-to-back changes coalesce into one render.
	if (!_renderTimer.isActive()) {
		_renderTimer.start();
	}
}

void ImageViewer::computeStatistics() {
	_histogram.clear();
	_stats.clear();
	_sampleMin = _sampleMax = 0.0f;
	_dataLo = 0.0f;
	_dataHi = kFullRange;
	_stretchLo = 0.0f;
	_stretchHi = kFullRange;

	if (_image.empty() || _header.width <= 0 || _header.height <= 0) {
		return;
	}

	const auto numChannels = static_cast<int>(_image.size());
	_histogram.assign(numChannels, HistogramBins{});
	_stats.assign(numChannels, ChannelStats{});

	auto combinedLo = std::numeric_limits<float>::infinity();
	auto combinedHi = -std::numeric_limits<float>::infinity();

	// Single pass per channel: histogram + sum + sum-of-squares + min/max.
	for (int c = 0; c < numChannels; ++c) {
		auto& bins = _histogram[c];
		auto& s = _stats[c];

		double sum = 0.0;
		double sumSq = 0.0;
		std::uint64_t count = 0;
		auto mn = std::numeric_limits<float>::infinity();
		auto mx = -std::numeric_limits<float>::infinity();

		for (const auto& row : _image[c]) {
			for (const auto v : row) {
				if (v < mn) {
					mn = v;
				}
				if (v > mx) {
					mx = v;
				}
				sum += v;
				sumSq += static_cast<double>(v) * v;
				++count;

				// Bin index in [0, 255] mapping the working 0..kFullRange domain.
				const auto binF = (v / kFullRange) * 255.0f;
				const auto bin = static_cast<int>(std::clamp(binF, 0.0f, 255.0f));
				++bins[static_cast<std::size_t>(bin)];
			}
		}

		if (count == 0) {
			continue;
		}

		s.count = count;
		s.min = mn;
		s.max = mx;
		s.mean = static_cast<float>(sum / count);

		// σ = sqrt(E[X²] - (E[X])²), clamped to avoid tiny-negative rounding.
		const auto meanD = sum / count;
		const auto variance = std::max(0.0, sumSq / count - meanD * meanD);
		s.stddev = static_cast<float>(std::sqrt(variance));

		// Median — cumulative walk over the histogram, linearly interpolated
		// inside the containing bin for better resolution than bin-centre.
		const auto half = count / 2;
		std::uint64_t running = 0;
		for (int i = 0; i < 256; ++i) {
			const auto inBin = bins[static_cast<std::size_t>(i)];
			if (running + inBin >= half) {
				const auto frac = inBin > 0
					? static_cast<double>(half - running) / inBin
					: 0.5;
				s.median = static_cast<float>((i + frac) * (kFullRange / 256.0));
				break;
			}
			running += inBin;
		}

		// Background, noise, star-level: delegate to the engine's
		// get_background which uses the full 65536-bin histogram, iterative
		// sigma-clipped noise estimation, and proper star-level thresholds.
		astap::Background bck{};
		astap::core::get_background(c, _image,
			/*calc_hist=*/true, /*calc_noise_level=*/true, bck);
		s.background = static_cast<float>(bck.backgr);
		s.noise = static_cast<float>(bck.noise_level);
		s.starLevel = static_cast<float>(bck.star_level);

		if (mn < combinedLo) {
			combinedLo = mn;
		}
		if (mx > combinedHi) {
			combinedHi = mx;
		}
	}

	if (!std::isfinite(combinedLo) || !std::isfinite(combinedHi)
			|| combinedHi <= combinedLo) {
		_sampleMin = 0.0f;
		_sampleMax = kFullRange;
	} else {
		_sampleMin = combinedLo;
		_sampleMax = combinedHi;
	}

	// Default stretch: background + sigma-clipped noise from the engine.
	// Uses the channel with the highest noise (typically the noisiest
	// filter, often red) so no channel gets crushed.
	auto bestChannel = 0;
	for (int c = 1; c < numChannels; ++c) {
		if (_stats[c].noise > _stats[bestChannel].noise) {
			bestChannel = c;
		}
	}
	const auto& ref = _stats[bestChannel];

	float lo = ref.background - 1.0f * ref.noise;
	float hi = ref.background + 7.0f * ref.noise;
	// Guard rails for degenerate inputs (flat fields, all-black frames).
	lo = std::clamp(lo, 0.0f, kFullRange);
	hi = std::clamp(hi, 0.0f, kFullRange);
	if (hi <= lo) {
		lo = _sampleMin;
		hi = _sampleMax;
		if (hi <= lo) {
			hi = lo + 1.0f;
		}
	}

	_dataLo = lo;
	_dataHi = hi;
	_stretchLo = lo;
	_stretchHi = hi;
}

void ImageViewer::renderImage() {
	if (_image.empty() || _header.width <= 0 || _header.height <= 0) {
		_rendered = QImage();
		return;
	}

	const auto w = _header.width;
	const auto h = _header.height;
	const auto numChannels = static_cast<int>(_image.size());

	// Stretch parameters derived once per render
	const auto lo = _stretchLo;
	const auto span = std::max(_stretchHi - _stretchLo, 1.0e-6f);
	const auto invSpan = 1.0f / span;
	const auto invLog = 1.0f / kLogDenom;

	// Render to RGB888. Channel 0 is the gray fallback for mono sources.
	_rendered = QImage(w, h, QImage::Format_RGB888);

	const auto& r = _image[0];
	const auto& g = (numChannels >= 2) ? _image[1] : _image[0];
	const auto& b = (numChannels >= 3) ? _image[2] : _image[0];

	auto stretch = [lo, invSpan, invLog, this](float v) -> std::uint8_t {
		auto t = (v - lo) * invSpan;
		if (t <= 0.0f) {
			return 0;
		}
		if (t >= 1.0f) {
			return 255;
		}
		if (_logStretch) {
			t = std::log1p(t * 700.0f) * invLog;
			if (t > 1.0f) {
				t = 1.0f;
			}
		}
		return toByte(t);
	};

	// FITS counts y from the bottom; QImage from the top. Apply user flips
	// on top of the FITS convention so the default presentation is "up = up".
	for (int y = 0; y < h; ++y) {
		const auto srcY = _flipV ? y : (h - 1) - y;
		auto* dst = _rendered.scanLine(y);
		const auto& rrow = r[srcY];
		const auto& grow = g[srcY];
		const auto& brow = b[srcY];
		for (int x = 0; x < w; ++x) {
			const auto srcX = _flipH ? (w - 1) - x : x;
			dst[3 * x + 0] = stretch(rrow[srcX]);
			dst[3 * x + 1] = stretch(grow[srcX]);
			dst[3 * x + 2] = stretch(brow[srcX]);
		}
	}
}

QPointF ImageViewer::viewportToImage(QPointF viewportPoint, bool& inImage) const {
	if (_rendered.isNull()) {
		inImage = false;
		return {0.0, 0.0};
	}

	// Mirror paintEvent's transform to find the image rect on screen.
	auto effectiveZoom = _zoom;
	if (_fitMode) {
		const auto sx = static_cast<double>(width()) / _rendered.width();
		const auto sy = static_cast<double>(height()) / _rendered.height();
		effectiveZoom = std::min(sx, sy);
	}
	const auto dstW = _rendered.width() * effectiveZoom;
	const auto dstH = _rendered.height() * effectiveZoom;
	const auto cx = (width() - dstW) * 0.5 + _pan.x();
	const auto cy = (height() - dstH) * 0.5 + _pan.y();

	// Render-space (top-left origin) image coords.
	const auto rx = (viewportPoint.x() - cx) / effectiveZoom;
	const auto ry = (viewportPoint.y() - cy) / effectiveZoom;

	const auto w = _rendered.width();
	const auto h = _rendered.height();
	inImage = (rx >= 0.0 && rx < w && ry >= 0.0 && ry < h);

	// Reverse the per-axis flips applied in renderImage.
	const auto srcRX = _flipH ? (w - 1 - rx) : rx;
	const auto srcRY = _flipV ? ry : (h - 1 - ry);

	// FITS pixels are 1-based, origin bottom-left, y up.
	return {srcRX + 1.0, srcRY + 1.0};
}

void ImageViewer::clampPan() {
	if (_rendered.isNull()) {
		return;
	}

	// Allow panning until the image edge meets the viewport edge.
	const auto effectiveZoom = _zoom;
	const auto dstW = _rendered.width() * effectiveZoom;
	const auto dstH = _rendered.height() * effectiveZoom;

	const auto maxX = std::max(0.0, (dstW - width()) * 0.5);
	const auto maxY = std::max(0.0, (dstH - height()) * 0.5);

	_pan.setX(std::clamp(_pan.x(), -maxX, maxX));
	_pan.setY(std::clamp(_pan.y(), -maxY, maxY));
}

void ImageViewer::applyZoomAtViewport(double newZoom, QPointF viewportPoint) {
	if (_rendered.isNull()) {
		return;
	}

	const auto clamped = std::clamp(newZoom, kMinZoom, kMaxZoom);

	// Convert viewport point → image-pixel coordinates at the *current* zoom.
	const auto cx = width() * 0.5 + _pan.x();
	const auto cy = height() * 0.5 + _pan.y();
	const auto imgX = (viewportPoint.x() - cx) / _zoom + _rendered.width() * 0.5;
	const auto imgY = (viewportPoint.y() - cy) / _zoom + _rendered.height() * 0.5;

	// Apply new zoom and recompute pan so the same image pixel stays under
	// the cursor.
	_zoom = clamped;
	_fitMode = false;
	const auto newCx = viewportPoint.x() - (imgX - _rendered.width() * 0.5) * _zoom;
	const auto newCy = viewportPoint.y() - (imgY - _rendered.height() * 0.5) * _zoom;
	_pan.setX(newCx - width() * 0.5);
	_pan.setY(newCy - height() * 0.5);

	clampPan();
	update();
}

} // namespace astap::gui
