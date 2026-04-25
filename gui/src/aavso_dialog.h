///----------------------------------------
///      @file aavso_dialog.h
///   @ingroup ASTAP++
///     @brief Modeless dialog for generating a single-frame AAVSO Extended
///            File Format report from the currently-loaded image.
///   @details Phase 5a of the Pascal `unit_aavso.pas` port. The user picks
///            three stars on the image (variable, check, comp), enters
///            naming + observation metadata, and the dialog formats an
///            AAVSO report and copies it to the clipboard or saves to file.
///
///            Multi-frame time-series workflow (Pascal's @c listview7)
///            lives in Phase 5b; this dialog is intentionally one row.
///    @author Created by John Stephen on 4/24/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include "image_viewer.h"

#include <QDialog>
#include <QFutureWatcher>
#include <QPointF>
#include <QString>

#include <array>
#include <memory>
#include <optional>

namespace astap::gui { class LightCurveView; }

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLabel;
class QLineEdit;
class QListWidget;
class QPushButton;
class QTableWidget;
class QTabWidget;
class QPlainTextEdit;

namespace astap::core { struct AavsoCollectResult; }

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @class AavsoDialog
/// @brief Modeless single-frame AAVSO photometry report generator.
///----------------------------------------

class AavsoDialog final : public QDialog {
	Q_OBJECT

public:
	explicit AavsoDialog(QWidget* parent, ImageViewer* viewer);
	~AavsoDialog() override;

private slots:
	void onPickButtonClicked();
	void onPicked(QPointF imagePos, int role);
	void onEnsembleToggled(bool ensemble);
	void onGenerate();   // primary: copy to clipboard
	void onSaveAs();     // save to .txt
	void onFetchVspMagnitude();

	// Multi-frame slots
	void onAddFiles();
	void onRemoveSelectedFiles();
	void onClearFiles();
	void onMeasureAllFrames();
	void onMeasureFinished();
	void onCopyMultiReport();
	void onSaveMultiReport();

private:
	struct StarRow {
		ImageViewer::PickRole role = ImageViewer::PickRole::None;
		QLineEdit* nameEdit = nullptr;
		QPushButton* pickButton = nullptr;
		QLabel* magLabel = nullptr;
		std::optional<QPointF> position;        // 1-based FITS pixel
		double instrumentalMag = 0.0;            // computed at pick time
		double snr = 0.0;
		double ra  = 0.0;                        // celestial coords of pick (rad)
		double dec = 0.0;
	};

	void measureAtPick(StarRow& row, QPointF imagePos);
	[[nodiscard]] bool validateForReport(QString& whyOut) const;
	[[nodiscard]] QString buildReport() const;
	void setStatus(const QString& msg, bool ok);

	void loadSettings();
	void saveSettings() const;

	// Multi-frame helpers
	void rebuildResultsTable();
	void appendMultiLog(const QString& line);

	ImageViewer* _viewer = nullptr;

	std::array<StarRow, 3> _rows{};   // 0=Variable, 1=Check, 2=Comp

	// Comp-only fields
	QDoubleSpinBox* _compCatalogMag = nullptr;
	QPushButton* _fetchVspBtn = nullptr;
	QCheckBox* _ensembleCheck = nullptr;

	// Common metadata
	QLineEdit* _observerCode = nullptr;
	QComboBox* _filter = nullptr;
	QComboBox* _delimiter = nullptr;
	QCheckBox* _hjd = nullptr;
	QCheckBox* _baaStyle = nullptr;
	QDoubleSpinBox* _deltaBv = nullptr;
	QDoubleSpinBox* _magnitudeSlope = nullptr;

	QLabel* _status = nullptr;
	QPushButton* _generateBtn = nullptr;
	QPushButton* _saveAsBtn = nullptr;

	// Multi-frame UI
	QTabWidget* _tabs = nullptr;
	QListWidget* _fileList = nullptr;
	QPushButton* _addFilesBtn = nullptr;
	QPushButton* _removeFilesBtn = nullptr;
	QPushButton* _clearFilesBtn = nullptr;
	QPushButton* _measureAllBtn = nullptr;
	QTableWidget* _resultsTable = nullptr;
	LightCurveView* _lightCurve = nullptr;
	QPlainTextEdit* _multiLog = nullptr;
	QPushButton* _copyMultiBtn = nullptr;
	QPushButton* _saveMultiBtn = nullptr;
	QLabel* _multiStatus = nullptr;

	std::unique_ptr<QFutureWatcher<astap::core::AavsoCollectResult>> _measureWatcher;
	std::unique_ptr<astap::core::AavsoCollectResult> _measureResult;
};

} // namespace astap::gui
