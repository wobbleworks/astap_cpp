///----------------------------------------
///      @file solve_dialog.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the plate-solve parameter dialog.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "solve_dialog.h"

#include "../../src/core/globals.h"
#include "../../src/core/wcs.h"
#include "../../src/reference/star_database.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QSpinBox>
#include <QVBoxLayout>

#include <cmath>
#include <numbers>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

namespace {

constexpr double kRad2Deg = 180.0 / std::numbers::pi;
constexpr double kRad2Hours = 12.0 / std::numbers::pi;

// Standard ASTAP database abbreviations seen in the original GUI.
const char* const kDatabaseNames[] = {"d05", "d20", "d50", "d80", "v50", "g05"};

} // namespace

///----------------------------------------
/// MARK: SolveDialog
///----------------------------------------

SolveDialog::SolveDialog(QWidget* parent) :
	QDialog(parent) {

	// Title and modal behaviour
	setWindowTitle(tr("Plate Solve"));
	setModal(true);

	// Build widgets and pre-load from current solver settings
	buildLayout();
	loadFromGlobals();
}

void SolveDialog::buildLayout() {
	auto* root = new QVBoxLayout(this);

	// --- Database group -----------------------------------------------------
	auto* dbGroup = new QGroupBox(tr("Star database"), this);
	auto* dbLayout = new QFormLayout(dbGroup);

	auto* dbRow = new QHBoxLayout();
	_databasePath = new QLineEdit(dbGroup);
	_databasePath->setPlaceholderText(tr("Folder containing the catalog files"));
	_databaseBrowse = new QPushButton(tr("Browse…"), dbGroup);
	dbRow->addWidget(_databasePath, 1);
	dbRow->addWidget(_databaseBrowse);
	dbLayout->addRow(tr("Path:"), dbRow);

	_databaseName = new QComboBox(dbGroup);
	for (const auto* name : kDatabaseNames) {
		_databaseName->addItem(QString::fromLatin1(name));
	}
	_databaseName->setEditable(true);
	dbLayout->addRow(tr("Name:"), _databaseName);

	root->addWidget(dbGroup);

	// --- Field hint group ---------------------------------------------------
	auto* fieldGroup = new QGroupBox(tr("Field hint (optional)"), this);
	auto* fieldLayout = new QFormLayout(fieldGroup);

	_fovDeg = new QDoubleSpinBox(fieldGroup);
	_fovDeg->setRange(0.0, 180.0);
	_fovDeg->setSuffix(tr(" °"));
	_fovDeg->setDecimals(4);
	_fovDeg->setSpecialValueText(tr("auto"));
	fieldLayout->addRow(tr("FOV height:"), _fovDeg);

	_ra = new QLineEdit(fieldGroup);
	_ra->setPlaceholderText(tr("e.g. 12 34 56 or 12.5827 (hours)"));
	fieldLayout->addRow(tr("RA centre:"), _ra);

	_dec = new QLineEdit(fieldGroup);
	_dec->setPlaceholderText(tr("e.g. +12 34 56 or 12.5828 (degrees)"));
	fieldLayout->addRow(tr("Dec centre:"), _dec);

	_searchRadiusDeg = new QDoubleSpinBox(fieldGroup);
	_searchRadiusDeg->setRange(0.0, 180.0);
	_searchRadiusDeg->setSuffix(tr(" °"));
	_searchRadiusDeg->setDecimals(2);
	fieldLayout->addRow(tr("Search radius:"), _searchRadiusDeg);

	root->addWidget(fieldGroup);

	// --- Solver tuning group -----------------------------------------------
	auto* solverGroup = new QGroupBox(tr("Solver"), this);
	auto* solverLayout = new QFormLayout(solverGroup);

	_maxStars = new QSpinBox(solverGroup);
	_maxStars->setRange(50, 100000);
	_maxStars->setSingleStep(50);
	solverLayout->addRow(tr("Max stars:"), _maxStars);

	_downsample = new QComboBox(solverGroup);
	_downsample->addItem(tr("auto"), 0);
	_downsample->addItem(tr("1× (none)"), 1);
	_downsample->addItem(tr("2×"), 2);
	_downsample->addItem(tr("3×"), 3);
	_downsample->addItem(tr("4×"), 4);
	solverLayout->addRow(tr("Downsample:"), _downsample);

	_minStarArcsec = new QDoubleSpinBox(solverGroup);
	_minStarArcsec->setRange(0.0, 60.0);
	_minStarArcsec->setSingleStep(0.1);
	_minStarArcsec->setDecimals(2);
	_minStarArcsec->setSuffix(tr(" \""));
	solverLayout->addRow(tr("Min star size:"), _minStarArcsec);

	_quadTolerance = new QDoubleSpinBox(solverGroup);
	_quadTolerance->setRange(0.0001, 0.1);
	_quadTolerance->setSingleStep(0.001);
	_quadTolerance->setDecimals(4);
	solverLayout->addRow(tr("Quad tolerance:"), _quadTolerance);

	_forceOversize = new QCheckBox(tr("Force oversize search (slow)"), solverGroup);
	_addSip = new QCheckBox(tr("Add SIP distortion to solution"), solverGroup);
	_patternFilter = new QCheckBox(tr("Apply OSC pattern filter"), solverGroup);
	solverLayout->addRow(_forceOversize);
	solverLayout->addRow(_addSip);
	solverLayout->addRow(_patternFilter);

	root->addWidget(solverGroup);

	// --- Buttons ------------------------------------------------------------
	auto* buttons = new QDialogButtonBox(
		QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
	buttons->button(QDialogButtonBox::Ok)->setText(tr("Solve"));
	root->addWidget(buttons);

	// --- Wires --------------------------------------------------------------
	connect(_databaseBrowse, &QPushButton::clicked, this, &SolveDialog::browseDatabase);
	connect(buttons, &QDialogButtonBox::accepted, this, &QDialog::accept);
	connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
}

void SolveDialog::loadFromGlobals() {
	// Database location
	_databasePath->setText(QString::fromStdString(
		astap::reference::database_path.string()));
	const auto current = QString::fromStdString(astap::reference::name_database);
	const auto idx = _databaseName->findText(current);
	if (idx >= 0) {
		_databaseName->setCurrentIndex(idx);
	} else if (!current.isEmpty()) {
		_databaseName->setCurrentText(current);
	}

	// Field hint defaults
	_fovDeg->setValue(astap::fov_specified ? astap::search_fov_deg : 0.0);
	_searchRadiusDeg->setValue(astap::search_radius_deg);

	// Solver tuning
	_maxStars->setValue(astap::max_stars_setting);
	{
		const auto idx2 = _downsample->findData(astap::downsample_setting);
		_downsample->setCurrentIndex(idx2 >= 0 ? idx2 : 0);
	}
	_minStarArcsec->setValue(astap::min_star_size_arcsec);
	_quadTolerance->setValue(astap::quad_tolerance);
	_forceOversize->setChecked(astap::force_oversize);
	_addSip->setChecked(astap::add_sip);
	_patternFilter->setChecked(astap::check_pattern_filter);
}

void SolveDialog::prefillFromHeader(const astap::Header& head) {
	// If the loaded image already carries a WCS solution, surface it as the
	// hint so a re-solve uses the existing centre.
	if (head.cdelt2 != 0.0) {
		// cdelt2 is degrees/pixel (engine convention; see the note in
		// project CLAUDE). Multiplying by height yields degrees directly —
		// no /3600 factor.
		const auto fov = std::abs(head.cdelt2) * head.height;
		if (fov > 0.0) {
			_fovDeg->setValue(fov);
		}
	}
	if (head.ra0 != 0.0 || head.dec0 != 0.0) {
		// Format using the engine's own helpers for consistency.
		_ra->setText(QString::fromStdString(astap::core::prepare_ra(head.ra0, ":")));
		_dec->setText(QString::fromStdString(astap::core::prepare_dec(head.dec0, ":")));
	}
}

void SolveDialog::accept() {
	// Validate database path
	const auto dbPath = _databasePath->text().trimmed();
	if (dbPath.isEmpty()) {
		QMessageBox::warning(this, windowTitle(),
			tr("Please choose a star database folder."));
		return;
	}

	// Parse RA / Dec text via the engine's parsers (handles many formats).
	auto parseRA = [this](double& out) -> bool {
		const auto txt = _ra->text().trimmed();
		if (txt.isEmpty()) {
			out = 0.0;
			return true;
		}
		bool err = false;
		auto inp = txt.toStdString();
		astap::core::ra_text_to_radians(inp, out, err);
		return !err;
	};
	auto parseDec = [this](double& out) -> bool {
		const auto txt = _dec->text().trimmed();
		if (txt.isEmpty()) {
			out = 0.0;
			return true;
		}
		bool err = false;
		auto inp = txt.toStdString();
		astap::core::dec_text_to_radians(inp, out, err);
		return !err;
	};

	double ra = 0.0;
	double dec = 0.0;
	if (!parseRA(ra)) {
		QMessageBox::warning(this, windowTitle(),
			tr("Could not parse the RA value."));
		return;
	}
	if (!parseDec(dec)) {
		QMessageBox::warning(this, windowTitle(),
			tr("Could not parse the Dec value."));
		return;
	}

	// Push everything into engine globals — solve_image reads these.
	std::filesystem::path dbp(dbPath.toStdString());
	if (!dbp.empty() && dbp.native().back() != std::filesystem::path::preferred_separator) {
		dbp /= "";
	}
	astap::reference::database_path = dbp;
	astap::reference::name_database = _databaseName->currentText().toStdString();

	const auto fov = _fovDeg->value();
	astap::search_fov_deg = fov;
	astap::fov_specified = (fov > 0.0);
	astap::search_radius_deg = _searchRadiusDeg->value();

	astap::max_stars_setting = _maxStars->value();
	astap::downsample_setting = _downsample->currentData().toInt();
	astap::min_star_size_arcsec = _minStarArcsec->value();
	astap::quad_tolerance = _quadTolerance->value();
	astap::force_oversize = _forceOversize->isChecked();
	astap::add_sip = _addSip->isChecked();
	astap::check_pattern_filter = _patternFilter->isChecked();

	// Position hint goes into the Header (solve_image picks it up from there).
	if (!_ra->text().trimmed().isEmpty()) {
		astap::head.ra0 = ra;
	}
	if (!_dec->text().trimmed().isEmpty()) {
		astap::head.dec0 = dec;
	}

	QDialog::accept();
}

void SolveDialog::browseDatabase() {
	const auto picked = QFileDialog::getExistingDirectory(
		this, tr("Star database folder"), _databasePath->text());
	if (!picked.isEmpty()) {
		_databasePath->setText(picked);
	}
}

} // namespace astap::gui
