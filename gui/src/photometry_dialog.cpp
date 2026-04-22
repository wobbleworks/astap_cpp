///----------------------------------------
///      @file photometry_dialog.cpp
///   @ingroup ASTAP++
///     @brief PhotometryDialog implementation.
///    @author Created by John Stephen on 4/19/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "photometry_dialog.h"

#include "../../src/core/globals.h"
#include "../../src/core/photometry.h"

#include <QCheckBox>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFontDatabase>
#include <QFormLayout>
#include <QGroupBox>
#include <QGuiApplication>
#include <QHBoxLayout>
#include <QLabel>
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

constexpr auto kKeyApertureRatio = "photometry/apertureRatio";
constexpr auto kKeyExtended      = "photometry/extended";
constexpr auto kKeyAnnulus       = "photometry/annulusRadius";

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

PhotometryDialog::PhotometryDialog(QWidget* parent) :
	QDialog(parent) {

	setWindowTitle(tr("Photometric Calibration"));
	setModal(false);
	resize(540, 520);
	buildLayout();
	prefillFromSettings();
}

/// MARK: - UI

void PhotometryDialog::buildLayout() {
	auto* root = new QVBoxLayout(this);

	// Aperture and annulus settings.
	auto* apertureGroup = new QGroupBox(tr("Aperture"), this);
	auto* apertureForm  = new QFormLayout(apertureGroup);

	_apertureRatio = new QDoubleSpinBox(apertureGroup);
	_apertureRatio->setRange(0.5, 10.0);
	_apertureRatio->setSingleStep(0.1);
	_apertureRatio->setDecimals(1);
	_apertureRatio->setValue(2.0);
	_apertureRatio->setToolTip(tr("Aperture diameter as a multiple of the "
	                               "median star HFD. Smaller is better for "
	                               "point sources; larger for extended."));
	apertureForm->addRow(tr("Aperture (× HFD)"), _apertureRatio);

	_extendedMode = new QCheckBox(tr("Extended objects (use Max aperture)"),
	                               apertureGroup);
	_extendedMode->setToolTip(tr("When enabled, the aperture radius is set to "
	                              "99 px and the per-star aperture is grown to "
	                              "cover the full source. Required for galaxy "
	                              "and nebula SQM."));
	connect(_extendedMode, &QCheckBox::toggled,
	        this, &PhotometryDialog::onExtendedObjectsToggled);
	apertureForm->addRow(QString(), _extendedMode);

	_annulusRadius = new QSpinBox(apertureGroup);
	_annulusRadius->setRange(5, 50);
	_annulusRadius->setValue(14);
	_annulusRadius->setSuffix(tr(" px"));
	_annulusRadius->setToolTip(tr("Background-annulus radius. 14 px is the "
	                               "engine default."));
	apertureForm->addRow(tr("Annulus radius"), _annulusRadius);

	root->addWidget(apertureGroup);

	// Action row.
	auto* actionRow = new QHBoxLayout();
	_runButton = new QPushButton(tr("Run calibration"), this);
	_runButton->setDefault(true);
	connect(_runButton, &QPushButton::clicked,
	        this, &PhotometryDialog::runCalibration);
	actionRow->addStretch();
	actionRow->addWidget(_runButton);
	root->addLayout(actionRow);

	// Results.
	auto* resultGroup = new QGroupBox(tr("Result"), this);
	auto* resultForm  = new QFormLayout(resultGroup);

	const auto valueFont = QFontDatabase::systemFont(QFontDatabase::FixedFont);
	auto make_value = [&](const QString& placeholder) {
		auto* lbl = new QLabel(placeholder, resultGroup);
		lbl->setFont(valueFont);
		lbl->setTextInteractionFlags(Qt::TextSelectableByMouse);
		return lbl;
	};
	_mzeroValue    = make_value(tr("—"));
	_passbandValue = make_value(tr("—"));
	_starsValue    = make_value(tr("—"));
	_semValue      = make_value(tr("—"));
	_limMagnValue  = make_value(tr("—"));
	_apertureValue = make_value(tr("—"));

	resultForm->addRow(tr("MZERO"),              _mzeroValue);
	resultForm->addRow(tr("Passband"),           _passbandValue);
	resultForm->addRow(tr("Stars used"),         _starsValue);
	resultForm->addRow(tr("Standard error"),     _semValue);
	resultForm->addRow(tr("Limiting magnitude"), _limMagnValue);
	resultForm->addRow(tr("Aperture radius"),    _apertureValue);

	root->addWidget(resultGroup);

	// Log area.
	_log = new QTextEdit(this);
	_log->setReadOnly(true);
	_log->setFont(valueFont);
	_log->setMinimumHeight(120);
	root->addWidget(_log, 1);

	// Close button.
	auto* buttons = new QDialogButtonBox(QDialogButtonBox::Close, this);
	connect(buttons, &QDialogButtonBox::rejected, this, &PhotometryDialog::reject);
	root->addWidget(buttons);

	onExtendedObjectsToggled(_extendedMode->isChecked());
}

void PhotometryDialog::onExtendedObjectsToggled(bool on) {
	// Aperture ratio is ignored when extended-objects mode is on — the
	// engine uses a very-wide virtual aperture in that case.
	_apertureRatio->setEnabled(!on);
}

/// MARK: - Settings persistence

void PhotometryDialog::prefillFromSettings() {
	QSettings settings;
	_apertureRatio->setValue(settings.value(kKeyApertureRatio, 2.0).toDouble());
	_extendedMode ->setChecked(settings.value(kKeyExtended,   false).toBool());
	_annulusRadius->setValue(settings.value(kKeyAnnulus,          14).toInt());
	onExtendedObjectsToggled(_extendedMode->isChecked());
}

void PhotometryDialog::saveToSettings() const {
	QSettings settings;
	settings.setValue(kKeyApertureRatio, _apertureRatio->value());
	settings.setValue(kKeyExtended,      _extendedMode->isChecked());
	settings.setValue(kKeyAnnulus,       _annulusRadius->value());
}

/// MARK: - Calibration

void PhotometryDialog::runCalibration() {
	if (astap::head.naxis == 0 || astap::head.cd1_1 == 0.0) {
		_log->setPlainText(tr("No plate-solved image is loaded. "
		                      "Open and solve an image first."));
		return;
	}

	saveToSettings();

	// Force a fresh calibration regardless of the prior MZERO state — the
	// user just pressed Run, so they want the current settings applied.
	astap::head.mzero = 0.0;

	const auto apert = _extendedMode->isChecked() ? 0.0
	                                               : _apertureRatio->value();
	auto aperture_state = -1.0;  // sentinel so calibrate_photometry re-runs
	auto annulus_radius_state = 0;
	auto passband_state = std::string{};

	_runButton->setEnabled(false);
	QGuiApplication::setOverrideCursor(Qt::WaitCursor);
	const auto memoStart = astap::memo1_lines.size();

	// annulus_radius_setting is a multiplier-of-HFD in the original; the
	// engine computes the true annulus from median HFD when apert != 0.
	// For simplicity, expose the final annulus radius directly and pass a
	// unit multiplier when in extended mode.
	astap::core::calibrate_photometry(
		astap::img_loaded,
		astap::memo1_lines,
		astap::head,
		astap::bck,
		/*update=*/true,
		/*aperture_ratio_setting=*/apert,
		/*annulus_radius_setting=*/static_cast<double>(_annulusRadius->value()),
		aperture_state,
		annulus_radius_state,
		passband_state);

	QGuiApplication::restoreOverrideCursor();
	_runButton->setEnabled(true);

	displayResult();
	const auto newLog = new_memo_lines(memoStart);
	if (!newLog.isEmpty()) {
		_log->setPlainText(newLog);
	}
}

/// MARK: - Display

void PhotometryDialog::displayResult() {
	if (astap::head.mzero == 0.0) {
		_mzeroValue   ->setText(tr("— (calibration failed)"));
		_passbandValue->setText(tr("—"));
		_starsValue   ->setText(tr("—"));
		_semValue     ->setText(tr("—"));
		_limMagnValue ->setText(tr("—"));
		_apertureValue->setText(tr("—"));
		return;
	}

	_mzeroValue->setText(QString::number(astap::head.mzero, 'f', 3));
	_passbandValue->setText(QString::fromStdString(astap::head.passband_database));
	_apertureValue->setText(
		QString::number(astap::head.mzero_radius, 'f', 2) + tr(" px"));

	if (astap::head.magn_limit > 0.0 && std::isfinite(astap::head.magn_limit)) {
		_limMagnValue->setText(QString::number(astap::head.magn_limit, 'f', 2));
	} else {
		_limMagnValue->setText(tr("—"));
	}

	// Stars-used and SEM are not pushed into head; they appear only in the
	// engine log. Defer to the log pane for those values.
	_starsValue->setText(tr("see log"));
	_semValue  ->setText(tr("see log"));
}

} // namespace astap::gui
