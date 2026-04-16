///----------------------------------------
///      @file controls_panel.h
///   @ingroup ASTAP++
///     @brief Right-hand display controls (histogram, stretch, flip).
///   @details Hosts a @ref HistogramWidget and the slider/spin/checkbox
///            controls that drive the attached @ref ImageViewer's stretch
///            range, stretch mode, and flip flags. The panel observes the
///            viewer for external changes (loaded image, autoStretch from
///            menu) and refreshes its widgets to match.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <QPointer>
#include <QWidget>

class QCheckBox;
class QPushButton;
class QSlider;
class QSpinBox;
class QTableWidget;

///----------------------------------------
namespace astap::gui {
///----------------------------------------

class HistogramWidget;
class ImageViewer;

///----------------------------------------
/// @class ControlsPanel
/// @brief Composite widget exposing histogram + stretch/flip controls.
///----------------------------------------

class ControlsPanel final : public QWidget {
	Q_OBJECT

public:
	explicit ControlsPanel(QWidget* parent = nullptr);
	~ControlsPanel() override = default;

	/// @brief Bind the panel to a viewer; controls drive its display state.
	void attachViewer(ImageViewer* viewer);

	[[nodiscard]] HistogramWidget* histogramWidget() const noexcept { return _histogram; }

private slots:
	void onLoSliderChanged(int v);
	void onHiSliderChanged(int v);
	void onLoSpinChanged(int v);
	void onHiSpinChanged(int v);
	void onLogToggled(bool on);
	void onAutoClicked();
	void onFlipHToggled(bool on);
	void onFlipVToggled(bool on);

	void onViewerLoaded();
	void onViewerStretchChanged();
	void onViewerFlipChanged();

private:
	void syncFromViewer();
	void pushStretchToViewer();
	void refreshStats();

	QPointer<ImageViewer> _viewer;
	HistogramWidget* _histogram = nullptr;

	QSlider* _loSlider = nullptr;
	QSlider* _hiSlider = nullptr;
	QSpinBox* _loSpin = nullptr;
	QSpinBox* _hiSpin = nullptr;
	QPushButton* _autoButton = nullptr;
	QCheckBox* _logCheck = nullptr;
	QCheckBox* _flipHCheck = nullptr;
	QCheckBox* _flipVCheck = nullptr;

	QTableWidget* _statsTable = nullptr;

	bool _suppress = false;  // re-entrancy guard during programmatic updates
};

} // namespace astap::gui
