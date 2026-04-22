///----------------------------------------
///      @file sqm_dialog.cpp
///   @ingroup ASTAP++
///     @brief SqmDialog implementation.
///    @author Created by John Stephen on 4/19/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "sqm_dialog.h"

#include "../../src/core/globals.h"
#include "../../src/core/sqm.h"

#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFontDatabase>
#include <QFormLayout>
#include <QGroupBox>
#include <QGuiApplication>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSettings>
#include <QSpinBox>
#include <QTextEdit>
#include <QVBoxLayout>

#include <cmath>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

namespace {

constexpr auto kKeyLatitude    = "sqm/latitude";
constexpr auto kKeyLongitude   = "sqm/longitude";
constexpr auto kKeyTemperature = "sqm/temperature_c";
constexpr auto kKeyPressure    = "sqm/pressure_hpa";
constexpr auto kKeyPedestal    = "sqm/pedestal";

/// @brief Read a batch of messages appended to @c astap::memo1_lines since the
///        given count; returns them joined by newline.
[[nodiscard]] QString new_memo_lines(std::size_t since) {
	QString out;
	for (auto i = since; i < astap::memo1_lines.size(); ++i) {
		if (!out.isEmpty()) {
			out.append('\n');
		}
		out.append(QString::fromStdString(astap::memo1_lines[i]));
	}
	return out;
}

}  // namespace

/// MARK: - Construction

SqmDialog::SqmDialog(QWidget* parent) :
	QDialog(parent) {

	setWindowTitle(tr("Sky Quality Meter"));
	setModal(false);
	resize(520, 520);
	buildLayout();
	prefillFromSettings();
}

/// MARK: - UI

void SqmDialog::buildLayout() {
	auto* root = new QVBoxLayout(this);

	// Site + atmosphere inputs.
	auto* siteGroup = new QGroupBox(tr("Observing site"), this);
	auto* siteForm  = new QFormLayout(siteGroup);

	_latitude = new QLineEdit(siteGroup);
	_latitude->setPlaceholderText(tr("e.g. 52 22 15 N"));
	siteForm->addRow(tr("Latitude"), _latitude);

	_longitude = new QLineEdit(siteGroup);
	_longitude->setPlaceholderText(tr("e.g. 4 53 50 E"));
	siteForm->addRow(tr("Longitude"), _longitude);

	_temperature = new QDoubleSpinBox(siteGroup);
	_temperature->setRange(-50.0, 60.0);
	_temperature->setDecimals(1);
	_temperature->setSuffix(tr(" °C"));
	_temperature->setValue(10.0);
	siteForm->addRow(tr("Temperature"), _temperature);

	_pressure = new QDoubleSpinBox(siteGroup);
	_pressure->setRange(500.0, 1100.0);
	_pressure->setDecimals(0);
	_pressure->setSuffix(tr(" hPa"));
	_pressure->setValue(1010.0);
	siteForm->addRow(tr("Pressure"), _pressure);

	_pedestal = new QSpinBox(siteGroup);
	_pedestal->setRange(0, 65535);
	_pedestal->setValue(0);
	_pedestal->setToolTip(tr("Dark-frame pedestal ADU already applied by calibration."));
	siteForm->addRow(tr("Pedestal (ADU)"), _pedestal);

	root->addWidget(siteGroup);

	// Compute button.
	auto* actionRow = new QHBoxLayout();
	_computeButton = new QPushButton(tr("Compute"), this);
	_computeButton->setDefault(true);
	connect(_computeButton, &QPushButton::clicked, this, &SqmDialog::computeSqm);
	actionRow->addStretch();
	actionRow->addWidget(_computeButton);
	root->addLayout(actionRow);

	// Result grid.
	auto* resultGroup = new QGroupBox(tr("Result"), this);
	auto* resultForm  = new QFormLayout(resultGroup);

	const auto valueFont = QFontDatabase::systemFont(QFontDatabase::FixedFont);
	auto make_value = [&](const QString& placeholder) {
		auto* lbl = new QLabel(placeholder, resultGroup);
		lbl->setFont(valueFont);
		lbl->setTextInteractionFlags(Qt::TextSelectableByMouse);
		return lbl;
	};
	_sqmValue      = make_value(tr("—"));
	_bortleValue   = make_value(tr("—"));
	_altitudeValue = make_value(tr("—"));
	_airmassValue  = make_value(tr("—"));
	_passbandValue = make_value(tr("—"));
	_starsValue    = make_value(tr("—"));

	resultForm->addRow(tr("SQM"),          _sqmValue);
	resultForm->addRow(tr("Bortle class"), _bortleValue);
	resultForm->addRow(tr("Altitude"),     _altitudeValue);
	resultForm->addRow(tr("Airmass"),      _airmassValue);
	resultForm->addRow(tr("Passband"),     _passbandValue);
	resultForm->addRow(tr("Stars used"),   _starsValue);

	root->addWidget(resultGroup);

	// Log area.
	_log = new QTextEdit(this);
	_log->setReadOnly(true);
	_log->setFont(valueFont);
	_log->setMinimumHeight(100);
	root->addWidget(_log, 1);

	// Close button.
	auto* buttons = new QDialogButtonBox(QDialogButtonBox::Close, this);
	connect(buttons, &QDialogButtonBox::rejected, this, &SqmDialog::reject);
	root->addWidget(buttons);
}

/// MARK: - Settings persistence

void SqmDialog::prefillFromSettings() {
	QSettings settings;
	_latitude  ->setText(settings.value(kKeyLatitude,   "").toString());
	_longitude ->setText(settings.value(kKeyLongitude,  "").toString());
	_temperature->setValue(settings.value(kKeyTemperature, 10.0).toDouble());
	_pressure   ->setValue(settings.value(kKeyPressure,  1010.0).toDouble());
	_pedestal   ->setValue(settings.value(kKeyPedestal,      0).toInt());
}

void SqmDialog::saveToSettings() const {
	QSettings settings;
	settings.setValue(kKeyLatitude,    _latitude->text());
	settings.setValue(kKeyLongitude,   _longitude->text());
	settings.setValue(kKeyTemperature, _temperature->value());
	settings.setValue(kKeyPressure,    _pressure->value());
	settings.setValue(kKeyPedestal,    _pedestal->value());
}

/// MARK: - Engine wire-up

void SqmDialog::applyInputsToEngine() const {
	// These globals feed into calculate_az_alt / airmass_calc inside
	// calculate_sqm. Sitelat/long are sexagesimal text; the engine parses
	// them on demand.
	astap::core::sitelat      = _latitude->text().toStdString();
	astap::core::sitelong     = _longitude->text().toStdString();
	astap::core::temperature_c = _temperature->value();
	astap::core::pressure_hpa  = _pressure->value();
}

/// MARK: - Compute

void SqmDialog::computeSqm() {
	// Guard — need an image with a WCS solution to compute altitude and to
	// perform photometric flux calibration.
	if (astap::head.naxis == 0 || astap::head.cd1_1 == 0.0) {
		_log->setPlainText(tr("No plate-solved image is loaded. "
		                      "Open and solve an image first."));
		return;
	}

	saveToSettings();
	applyInputsToEngine();

	_computeButton->setEnabled(false);
	QGuiApplication::setOverrideCursor(Qt::WaitCursor);
	const auto memoStart = astap::memo1_lines.size();

	auto pedestal = _pedestal->value();
	const auto ok = astap::core::calculate_sqm(
		astap::head, /*get_backgr=*/true, /*get_histogr=*/true,
		pedestal);

	QGuiApplication::restoreOverrideCursor();
	_computeButton->setEnabled(true);

	displayResult(ok);

	// Append any new messages from the engine's memo.
	const auto newLog = new_memo_lines(memoStart);
	if (!newLog.isEmpty()) {
		_log->setPlainText(newLog);
	}
}

/// MARK: - Display

void SqmDialog::displayResult(bool ok) {
	if (!ok) {
		_sqmValue     ->setText(tr("— (calibration failed)"));
		_bortleValue  ->setText(tr("—"));
		_altitudeValue->setText(tr("—"));
		_airmassValue ->setText(tr("—"));
		_passbandValue->setText(tr("—"));
		_starsValue   ->setText(tr("—"));
		return;
	}

	_sqmValue->setText(QString::number(astap::core::sqmfloat, 'f', 2) +
	                   tr(" mag/arcsec²"));
	_bortleValue->setText(QString::fromStdString(
		astap::core::bortle(astap::core::sqmfloat)));

	// Altitude / airmass come from calculate_az_alt / airmass_calc inside
	// calculate_sqm.
	const auto alt = astap::core::altitudefloat;
	_altitudeValue->setText(QString::number(alt, 'f', 2) + tr(" °"));

	const auto airmass = astap::core::airmass;
	if (airmass > 0.0 && std::isfinite(airmass)) {
		_airmassValue->setText(QString::number(airmass, 'f', 2));
	} else {
		_airmassValue->setText(tr("—"));
	}

	_passbandValue->setText(QString::fromStdString(astap::head.passband_database));

	// stars_measured isn't exposed through the SQM path — best-effort parse
	// from the engine's most recent log line.
	_starsValue->setText(tr("see log"));
}

} // namespace astap::gui
