///----------------------------------------
///      @file photometry_dialog.h
///   @ingroup ASTAP++
///     @brief Photometry dialog: runs a photometric flux-calibration pass
///            and reports MZERO, passband, star count, standard error, and
///            limiting magnitude.
///    @author Created by John Stephen on 4/19/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <QDialog>

class QCheckBox;
class QDoubleSpinBox;
class QLabel;
class QPushButton;
class QSpinBox;
class QTextEdit;

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @class PhotometryDialog
/// @brief Modeless dialog for running a photometric calibration pass.
/// @details Lets the user pick an aperture-ratio (in HFD multiples) and an
///          annulus radius, then calls @c astap::core::calibrate_photometry
///          on the current image. Displays MZERO, passband, stars used, the
///          standard error, and the limiting magnitude (if reported).
///          Settings persist via QSettings.
///----------------------------------------

class PhotometryDialog final : public QDialog {
	Q_OBJECT

public:
	explicit PhotometryDialog(QWidget* parent = nullptr);
	~PhotometryDialog() override = default;

	/// @brief Pre-fill fields from persisted defaults.
	void prefillFromSettings();

private slots:
	void runCalibration();
	void onExtendedObjectsToggled(bool on);

private:
	void buildLayout();
	void saveToSettings() const;
	void displayResult();

	QDoubleSpinBox* _apertureRatio = nullptr;   ///< 0 means "Max" (extended)
	QCheckBox*      _extendedMode  = nullptr;
	QSpinBox*       _annulusRadius = nullptr;

	QPushButton* _runButton = nullptr;

	QLabel* _mzeroValue       = nullptr;
	QLabel* _passbandValue    = nullptr;
	QLabel* _starsValue       = nullptr;
	QLabel* _semValue         = nullptr;
	QLabel* _limMagnValue     = nullptr;
	QLabel* _apertureValue    = nullptr;
	QTextEdit* _log           = nullptr;
};

} // namespace astap::gui
