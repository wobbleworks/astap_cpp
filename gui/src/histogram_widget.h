///----------------------------------------
///      @file histogram_widget.h
///   @ingroup ASTAP++
///     @brief Histogram display widget for the controls panel.
///   @details Plots a 256-bin per-channel histogram on a log Y-scale, with
///            vertical guides for the current black/white stretch points.
///            Read-only display; the sliders in @ref ControlsPanel drive
///            stretch range changes.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include "image_viewer.h"

#include <QPointer>
#include <QWidget>

#include <vector>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @class HistogramWidget
/// @brief Per-channel histogram chart with stretch-range overlay.
///----------------------------------------

class HistogramWidget final : public QWidget {
	Q_OBJECT

public:
	explicit HistogramWidget(QWidget* parent = nullptr);
	~HistogramWidget() override = default;

	/// @brief Bind the widget to a viewer; subsequent updates pull from it.
	void attachViewer(ImageViewer* viewer);

	[[nodiscard]] QSize minimumSizeHint() const override { return {180, 90}; }
	[[nodiscard]] QSize sizeHint() const override { return {220, 110}; }

protected:
	void paintEvent(QPaintEvent* event) override;

private slots:
	void onImageLoaded();
	void onStretchChanged();

private:
	QPointer<ImageViewer> _viewer;
};

} // namespace astap::gui
