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

#include <QTemporaryDir>
#include <QWidget>

#include <memory>

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QLineEdit;
class QProgressBar;
class QPushButton;
class QSpinBox;
class QTabWidget;
class QTableWidget;

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
	// Dark / flat candidate-library management (mirrors the Pascal dark/flat tabs).
	void addDarkFiles();
	void addFlatFiles();
	void removeDarks();
	void removeFlats();
	void clearDarks();
	void clearFlats();

private:
	void buildLightsTab();
	void buildDarksTab();
	void buildFlatsTab();
	void buildSettingsTab();
	QWidget* buildClassifyBar();
	void applySettingsToEngine();
	void applyLibrariesToEngine();
	void addMasterFiles(QTableWidget* table, const QString& title);
	void removeSelectedRows(QTableWidget* table);
	void refreshCompatibility();

	QTabWidget* _tabs = nullptr;

	// Lights tab
	QTableWidget* _fileTable = nullptr;
	QPushButton* _addButton = nullptr;
	QPushButton* _removeButton = nullptr;
	QPushButton* _clearButton = nullptr;

	// Darks / Flats candidate tables (the master library the engine matches from).
	QTableWidget* _darkTable = nullptr;
	QTableWidget* _flatTable = nullptr;

	// Shared "Classify by" controls (Pascal classify_dark_* / classify_flat_filter
	// + delta_temp).
	QCheckBox* _classDarkExp = nullptr;
	QCheckBox* _classDarkTemp = nullptr;
	QCheckBox* _classDarkGain = nullptr;
	QCheckBox* _classFlatFilter = nullptr;
	QSpinBox*  _deltaTemp = nullptr;

	// Settings tab
	QComboBox* _methodCombo = nullptr;
	QComboBox* _alignmentCombo = nullptr;
	QDoubleSpinBox* _sigmaFactor = nullptr;
	QSpinBox* _maxStars = nullptr;

	// Footer
	QPushButton* _stackButton = nullptr;
	QProgressBar* _progress = nullptr;
	QLabel* _phaseLabel = nullptr;   ///< Current phase for multi-phase runs.

	ImageViewer* _viewer = nullptr;

	/// Temp dir for LRGB auto-chain's intermediate channel masters. Lazy
	/// initialised on first use; cleaned up when the window is destroyed.
	std::unique_ptr<QTemporaryDir> _tempDir;
};

} // namespace astap::gui
