///----------------------------------------
///      @file stack_window.h
///   @ingroup ASTAP++
///     @brief Image stacking window.
///   @details Phase 5b: Lights / Calibration / Settings tabs. Drives the
///            engine's stack_average / stack_sigmaclip / stack_LRGB passes
///            with configurable alignment mode and master dark / flat.
///            Result is pushed to the main viewer.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include "../../src/types.h"

#include <QWidget>

class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QLineEdit;
class QListWidget;
class QProgressBar;
class QPushButton;
class QSpinBox;
class QTabWidget;

///----------------------------------------
namespace astap::gui {
///----------------------------------------

class ImageViewer;

///----------------------------------------
/// @class StackWindow
/// @brief Non-modal stacking window with file list and stack controls.
///----------------------------------------

class StackWindow final : public QWidget {
	Q_OBJECT

public:
	explicit StackWindow(QWidget* parent = nullptr);
	~StackWindow() override = default;

	/// @brief Set the viewer that receives the stacked result.
	void setViewer(ImageViewer* viewer) { _viewer = viewer; }

signals:
	/// @brief Emitted after a successful stack so the main window can
	///        update title / status.
	void stackCompleted(int frameCount);

private slots:
	void addFiles();
	void removeSelected();
	void clearList();
	void startStack();
	void browseMasterDark();
	void browseMasterFlat();
	void clearMasterDark();
	void clearMasterFlat();

private:
	void buildLightsTab();
	void buildCalibrationTab();
	void buildSettingsTab();
	void applySettingsToEngine();

	QTabWidget* _tabs = nullptr;

	// Lights tab
	QListWidget* _fileList = nullptr;
	QPushButton* _addButton = nullptr;
	QPushButton* _removeButton = nullptr;
	QPushButton* _clearButton = nullptr;

	// Calibration tab
	QLineEdit* _darkPath = nullptr;
	QLabel*    _darkStatus = nullptr;
	QPushButton* _darkBrowseButton = nullptr;
	QPushButton* _darkClearButton = nullptr;
	QLineEdit* _flatPath = nullptr;
	QLabel*    _flatStatus = nullptr;
	QPushButton* _flatBrowseButton = nullptr;
	QPushButton* _flatClearButton = nullptr;

	// Settings tab
	QComboBox* _methodCombo = nullptr;
	QComboBox* _alignmentCombo = nullptr;
	QDoubleSpinBox* _sigmaFactor = nullptr;
	QSpinBox* _maxStars = nullptr;

	// Footer
	QPushButton* _stackButton = nullptr;
	QProgressBar* _progress = nullptr;

	ImageViewer* _viewer = nullptr;
};

} // namespace astap::gui
