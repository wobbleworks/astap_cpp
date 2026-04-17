///----------------------------------------
///      @file image_inspector_dialog.h
///   @ingroup ASTAP++
///     @brief Image Inspector: per-cell HFD / FWHM / star-count grid for
///            diagnosing focus tilt and field aberrations.
///    @author Created by John Stephen on 4/17/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include "star_detector.h"

#include <QDialog>
#include <QFutureWatcher>

#include <memory>

class QLabel;
class QTableWidget;

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @class ImageInspectorDialog
/// @brief Modeless dialog showing a 3x3 grid of star-shape statistics.
/// @details Runs detect_stars() on a worker thread, bins detections into
///          a 3x3 image-space grid, then displays per-cell star count,
///          median HFD, and median FWHM. Cells are colour-coded by HFD
///          relative to the sharpest cell, which makes tilt / backspace
///          / collimation issues visually obvious.
///----------------------------------------

class ImageInspectorDialog final : public QDialog {
	Q_OBJECT

public:
	explicit ImageInspectorDialog(QWidget* parent = nullptr);
	~ImageInspectorDialog() override;

	/// @brief Kick off detection on the currently loaded image.
	void analyseCurrentImage();

private slots:
	void onDetectionFinished();

private:
	void buildUi();
	void populateTable(const DetectionResult& result);

	QTableWidget* _table = nullptr;
	QLabel* _summary = nullptr;
	QLabel* _status = nullptr;

	std::unique_ptr<QFutureWatcher<DetectionResult>> _watcher;
	int _imageWidth = 0;
	int _imageHeight = 0;
};

} // namespace astap::gui
