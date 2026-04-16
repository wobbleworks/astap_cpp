///----------------------------------------
///      @file solve_dialog.h
///   @ingroup ASTAP++
///     @brief Modal dialog collecting plate-solver parameters.
///   @details On accept, writes the user's choices into the engine's solver
///            settings globals and Header position hints. The caller then
///            kicks off @ref astap::solving::solve_image on a worker thread.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include "../../src/types.h"

#include <QDialog>

class QCheckBox;
class QComboBox;
class QDoubleSpinBox;
class QLineEdit;
class QPushButton;
class QSpinBox;

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @class SolveDialog
/// @brief Collects search hints and solver parameters before a plate solve.
///----------------------------------------

class SolveDialog final : public QDialog {
	Q_OBJECT

public:
	explicit SolveDialog(QWidget* parent = nullptr);
	~SolveDialog() override = default;

	/// @brief Pre-fill fields from the currently loaded header (if it carries
	///        a previous WCS solution, scale + centre come from there).
	void prefillFromHeader(const astap::Header& head);

protected:
	void accept() override;

private slots:
	void browseDatabase();

private:
	void buildLayout();
	void loadFromGlobals();

	QLineEdit* _databasePath = nullptr;
	QPushButton* _databaseBrowse = nullptr;
	QComboBox* _databaseName = nullptr;

	QDoubleSpinBox* _fovDeg = nullptr;
	QLineEdit* _ra = nullptr;
	QLineEdit* _dec = nullptr;
	QDoubleSpinBox* _searchRadiusDeg = nullptr;

	QSpinBox* _maxStars = nullptr;
	QComboBox* _downsample = nullptr;
	QDoubleSpinBox* _minStarArcsec = nullptr;
	QDoubleSpinBox* _quadTolerance = nullptr;

	QCheckBox* _forceOversize = nullptr;
	QCheckBox* _addSip = nullptr;
	QCheckBox* _patternFilter = nullptr;
};

} // namespace astap::gui
