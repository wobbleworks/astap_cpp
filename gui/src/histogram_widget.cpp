///----------------------------------------
///      @file histogram_widget.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the histogram display widget.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "histogram_widget.h"

#include <QPaintEvent>
#include <QPainter>
#include <QPainterPath>

#include <algorithm>
#include <cmath>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

namespace {

constexpr float kFullRange = 65535.0f;

// Overlay colours; alpha lets the channels combine visually without
// completely covering each other.
const QColor kColorMono{200, 200, 200, 180};
const QColor kColorR{220, 70, 70, 160};
const QColor kColorG{70, 200, 90, 160};
const QColor kColorB{90, 130, 220, 160};

const QColor kColorBackground{20, 20, 24};
const QColor kColorAxis{120, 120, 120};
const QColor kColorBlackPoint{255, 220, 60, 220};
const QColor kColorWhitePoint{60, 220, 255, 220};

} // namespace

///----------------------------------------
/// MARK: HistogramWidget
///----------------------------------------

HistogramWidget::HistogramWidget(QWidget* parent) :
	QWidget(parent) {

	// Solid background so log-axis paint doesn't see prior frames.
	setAttribute(Qt::WA_OpaquePaintEvent);
	setAutoFillBackground(false);
}

void HistogramWidget::attachViewer(ImageViewer* viewer) {
	if (_viewer == viewer) {
		return;
	}

	// Drop old connections (Qt::UniqueConnection makes the wires below idempotent
	// if the user reattaches the same viewer).
	if (_viewer) {
		disconnect(_viewer, nullptr, this, nullptr);
	}

	_viewer = viewer;
	if (_viewer) {
		connect(_viewer, &ImageViewer::imageLoaded,
			this, &HistogramWidget::onImageLoaded, Qt::UniqueConnection);
		connect(_viewer, &ImageViewer::stretchChanged,
			this, &HistogramWidget::onStretchChanged, Qt::UniqueConnection);
	}
	update();
}

void HistogramWidget::paintEvent(QPaintEvent* event) {
	QPainter painter(this);

	// Fill background
	painter.fillRect(event->rect(), kColorBackground);

	if (!_viewer || !_viewer->hasImage()) {
		// Empty hint
		painter.setPen(Qt::darkGray);
		painter.drawText(rect(), Qt::AlignCenter, tr("No image"));
		return;
	}

	const auto& bins = _viewer->histogram();
	if (bins.empty()) {
		return;
	}

	// Plot area: leave room for axis baseline.
	const auto plotRect = rect().adjusted(2, 2, -2, -10);
	const auto plotW = plotRect.width();
	const auto plotH = plotRect.height();
	if (plotW <= 0 || plotH <= 0) {
		return;
	}

	// Find the highest bin count across all channels for log-scale normalisation.
	std::uint32_t maxBin = 1;
	for (const auto& chBins : bins) {
		for (const auto count : chBins) {
			if (count > maxBin) {
				maxBin = count;
			}
		}
	}
	const auto logMax = std::log1p(static_cast<float>(maxBin));

	// Per-channel polyline drawn as a filled curve.
	auto drawChannel = [&](const ImageViewer::HistogramBins& chBins, const QColor& colour) {
		QPainterPath path;
		path.moveTo(plotRect.left(), plotRect.bottom());
		for (int i = 0; i < 256; ++i) {
			const auto x = plotRect.left() + (i / 255.0) * plotW;
			const auto t = std::log1p(static_cast<float>(chBins[static_cast<std::size_t>(i)])) / logMax;
			const auto y = plotRect.bottom() - t * plotH;
			path.lineTo(x, y);
		}
		path.lineTo(plotRect.right(), plotRect.bottom());
		path.closeSubpath();
		painter.fillPath(path, colour);
	};

	if (bins.size() == 1) {
		drawChannel(bins[0], kColorMono);
	} else {
		const auto end = std::min<std::size_t>(bins.size(), 3);
		const QColor cols[3] = {kColorR, kColorG, kColorB};
		for (std::size_t c = 0; c < end; ++c) {
			drawChannel(bins[c], cols[c]);
		}
	}

	// Axis baseline
	painter.setPen(kColorAxis);
	painter.drawLine(plotRect.bottomLeft(), plotRect.bottomRight());

	// Stretch endpoints — vertical guides at the lo/hi positions.
	auto markerX = [&](float v) {
		const auto t = std::clamp(v / kFullRange, 0.0f, 1.0f);
		return plotRect.left() + t * plotW;
	};

	painter.setPen(QPen(kColorBlackPoint, 1.0));
	const auto xLo = markerX(_viewer->stretchLo());
	painter.drawLine(QPointF(xLo, plotRect.top()), QPointF(xLo, plotRect.bottom()));

	painter.setPen(QPen(kColorWhitePoint, 1.0));
	const auto xHi = markerX(_viewer->stretchHi());
	painter.drawLine(QPointF(xHi, plotRect.top()), QPointF(xHi, plotRect.bottom()));
}

void HistogramWidget::onImageLoaded() {
	update();
}

void HistogramWidget::onStretchChanged() {
	update();
}

} // namespace astap::gui
