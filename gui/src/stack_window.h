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
#include "../../src/stacking/stack.h"   // MasterMetadata cache

#include <QHash>
#include <QString>
#include <QTemporaryDir>
#include <QWidget>

#include <memory>
#include <vector>

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
	// Dark / flat / flat-dark candidate-library management (mirrors the Pascal tabs).
	void addDarkFiles();
	void addFlatFiles();
	void addFlatDarkFiles();
	void removeDarks();
	void removeFlats();
	void removeFlatDarks();
	void clearDarks();
	void clearFlats();
	void clearFlatDarks();
	// "Replace check-marked by … master" — combine checked raw frames into master(s)
	// and swap them into the table (Pascal replace_by_master_dark/flat1Click).
	void replaceDarksByMaster();
	void replaceFlatsByMaster();

private:
	void buildLightsTab();
	void buildDarksTab();
	void buildFlatsTab();
	void buildFlatDarksTab();
	void buildSettingsTab();
	QWidget* buildClassifyBar();
	void applySettingsToEngine();
	void applyLibrariesToEngine();
	void addMasterFiles(QTableWidget* table, const QString& title);
	bool insertMasterRow(QTableWidget* table, const QString& path);
	void replaceRows(QTableWidget* table,
	                 const std::vector<QString>& consumed,
	                 const std::vector<QString>& created);
	void removeSelectedRows(QTableWidget* table);
	void refreshCompatibility();
	void updateCompatibility(QTableWidget* table, bool isFlat,
	                         const astap::Header& light,
	                         const astap::stacking::MasterMatchOptions& opts,
	                         bool haveLight);

	QTabWidget* _tabs = nullptr;

	// Lights tab
	QTableWidget* _fileTable = nullptr;
	QPushButton* _addButton = nullptr;
	QPushButton* _removeButton = nullptr;
	QPushButton* _clearButton = nullptr;

	// Darks / Flats candidate tables (the master library the engine matches from).
	// The Flat darks table holds the raw flat-dark/bias frames subtracted during
	// master-flat creation (Pascal listview4).
	QTableWidget* _darkTable = nullptr;
	QTableWidget* _flatTable = nullptr;
	QTableWidget* _flatDarkTable = nullptr;

	// Pre-analysed metadata cache (path → header fields) so the Compatibility
	// column doesn't re-read every candidate header on each classify change.
	QHash<QString, astap::stacking::MasterMetadata> _masterMeta;

	// Shared "Classify by" controls (Pascal classify_groupbox1). Light filter /
	// object classify the lights; the dark_* / flat_filter gate master matching.
	QCheckBox* _lightFilter = nullptr;   // activates LRGB (classify lights by filter)
	QCheckBox* _lightObject = nullptr;   // stack different objects in one run
	QCheckBox* _classDarkExp = nullptr;
	QCheckBox* _classDarkTemp = nullptr;
	QCheckBox* _classDarkGain = nullptr;
	QCheckBox* _classFlatFilter = nullptr;
	QSpinBox*  _deltaTemp = nullptr;

	// Per-tab "… for master creation" grouping options (Pascal classify_dark_date1 /
	// classify_flat_date1 / classify_flat_duration1). These govern how raw frames are
	// grouped into masters by the Replace buttons; they do not affect run-time master
	// selection (which always uses closest-JD).
	QCheckBox* _classDarkDate = nullptr;
	QCheckBox* _classFlatDate = nullptr;
	QCheckBox* _classFlatDuration = nullptr;

	// Settings tab
	QComboBox* _methodCombo = nullptr;
	QComboBox* _alignmentCombo = nullptr;
	QDoubleSpinBox* _sigmaFactor = nullptr;
	QSpinBox* _maxStars = nullptr;
	// LRGB colour-mix matrix (Pascal rr1..bb1). Rows = input channel R/G/B,
	// columns = output channel R/G/B; identity = a straight RGB combine.
	QDoubleSpinBox* _colorMix[3][3] = {};

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
