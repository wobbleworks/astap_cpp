///----------------------------------------
///      @file sqm_dialog.h
///   @ingroup ASTAP++
///     @brief Sky Quality Meter dialog: runs photometric flux calibration
///            and reports the zenith-equivalent sky-surface-brightness plus
///            Bortle classification.
///    @author Created by John Stephen on 4/19/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <QDialog>

class QDoubleSpinBox;
class QLabel;
class QLineEdit;
class QPushButton;
class QSpinBox;
class QTextEdit;

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @class SqmDialog
/// @brief Modeless dialog for the Sky Quality Meter feature.
/// @details Collects site latitude/longitude, temperature, pressure, and
///          dark-pedestal. On Compute, writes those into the engine globals
///          (@c astap::core::sitelat etc.) and calls
///          @c astap::core::calculate_sqm. Displays SQM value,
///          Bortle classification, image-centre altitude, and airmass.
///          Lat/long defaults are persisted via QSettings.
///----------------------------------------

class SqmDialog final : public QDialog {
	Q_OBJECT

public:
	explicit SqmDialog(QWidget* parent = nullptr);
	~SqmDialog() override = default;

	/// @brief Pre-fill fields from persisted defaults and current header.
	void prefillFromSettings();

private slots:
	void computeSqm();

private:
	void buildLayout();
	void saveToSettings() const;
	void applyInputsToEngine() const;
	void displayResult(bool ok);

	QLineEdit*      _latitude  = nullptr;
	QLineEdit*      _longitude = nullptr;
	QDoubleSpinBox* _temperature = nullptr;
	QDoubleSpinBox* _pressure    = nullptr;
	QSpinBox*       _pedestal    = nullptr;

	QPushButton* _computeButton = nullptr;

	QLabel*    _sqmValue      = nullptr;
	QLabel*    _bortleValue   = nullptr;
	QLabel*    _altitudeValue = nullptr;
	QLabel*    _airmassValue  = nullptr;
	QLabel*    _passbandValue = nullptr;
	QLabel*    _starsValue    = nullptr;
	QTextEdit* _log           = nullptr;
};

} // namespace astap::gui
