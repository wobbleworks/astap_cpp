///----------------------------------------
///      @file preferences_dialog.cpp
///   @ingroup ASTAP++
///     @brief PreferencesDialog implementation.
///    @author Created by John Stephen on 4/22/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "preferences_dialog.h"

#include "../../src/core/globals.h"
#include "../../src/reference/star_database.h"

#include <QComboBox>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QSettings>
#include <QTabWidget>
#include <QVBoxLayout>
#include <QWidget>

#include <filesystem>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

namespace {

// QSettings keys — mirror the rest of the app.
constexpr auto kKeyDatabasePath = "solver/databasePath";
constexpr auto kKeyDatabaseName = "solver/databaseName";
constexpr auto kKeySqmLatitude  = "sqm/latitude";
constexpr auto kKeySqmLongitude = "sqm/longitude";
constexpr auto kKeySqmTemperature = "sqm/temperature_c";
constexpr auto kKeySqmPressure    = "sqm/pressure_hpa";
constexpr auto kKeySiteElevation  = "site/elevation";
constexpr auto kKeySiteTelescope  = "site/telescope";
constexpr auto kKeySiteCamera     = "site/camera";
constexpr auto kKeyRecentFiles    = "files/recent";
constexpr auto kKeyGeometry       = "window/geometry";
constexpr auto kKeyWindowState    = "window/state";
constexpr auto kKeySplitter       = "window/splitter";

}  // namespace

PreferencesDialog::PreferencesDialog(QWidget* parent) :
	QDialog(parent) {

	setWindowTitle(tr("Preferences"));
	setModal(true);
	resize(520, 420);
	buildLayout();
	loadFromSettings();
}

/// MARK: - UI

void PreferencesDialog::buildLayout() {
	auto* root = new QVBoxLayout(this);
	auto* tabs = new QTabWidget(this);
	root->addWidget(tabs, 1);

	// --- Catalog tab ---
	{
		auto* page = new QWidget(tabs);
		auto* form = new QFormLayout(page);

		auto* pathRow = new QHBoxLayout();
		_databasePath = new QLineEdit(page);
		_databasePath->setPlaceholderText(tr("Directory containing g17, d80, "
			"d50, d05 star databases"));
		_browseButton = new QPushButton(tr("Browse…"), page);
		pathRow->addWidget(_databasePath, 1);
		pathRow->addWidget(_browseButton);
		form->addRow(tr("Star database path"), pathRow);

		_databaseName = new QComboBox(page);
		_databaseName->setEditable(true);
		_databaseName->addItems({
			QString(), tr("auto"), "g17", "g18", "v17", "v50",
			"d80", "d50", "d20", "d05",
		});
		_databaseName->setToolTip(tr(
			"Preferred database prefix for the solver. Leave blank for auto-select "
			"based on FOV."));
		form->addRow(tr("Preferred database"), _databaseName);

		connect(_browseButton, &QPushButton::clicked,
		        this, &PreferencesDialog::browseDatabase);

		tabs->addTab(page, tr("Catalog"));
	}

	// --- Site tab ---
	{
		auto* page = new QWidget(tabs);
		auto* form = new QFormLayout(page);
		auto* intro = new QLabel(tr(
			"Default observing-site values used by the Sky Quality Meter "
			"dialog when it first opens."), page);
		intro->setWordWrap(true);
		intro->setStyleSheet("color: gray;");
		form->addRow(intro);

		_latitude = new QLineEdit(page);
		_latitude->setPlaceholderText(tr("e.g. 52 22 15 N"));
		form->addRow(tr("Latitude"), _latitude);

		_longitude = new QLineEdit(page);
		_longitude->setPlaceholderText(tr("e.g. 4 53 50 E"));
		form->addRow(tr("Longitude"), _longitude);

		_elevation = new QLineEdit(page);
		_elevation->setPlaceholderText(tr("e.g. 25m"));
		_elevation->setToolTip(tr(
			"Site elevation; used for the BAA-style AAVSO #LOCATION line."));
		form->addRow(tr("Elevation"), _elevation);

		_temperature = new QDoubleSpinBox(page);
		_temperature->setRange(-50.0, 60.0);
		_temperature->setDecimals(1);
		_temperature->setSuffix(tr(" °C"));
		_temperature->setValue(10.0);
		form->addRow(tr("Temperature"), _temperature);

		_pressure = new QDoubleSpinBox(page);
		_pressure->setRange(500.0, 1100.0);
		_pressure->setDecimals(0);
		_pressure->setSuffix(tr(" hPa"));
		_pressure->setValue(1010.0);
		form->addRow(tr("Pressure"), _pressure);

		_telescope = new QLineEdit(page);
		_telescope->setPlaceholderText(tr("e.g. SCT 8\" f/10"));
		_telescope->setToolTip(tr(
			"Telescope description used by the BAA-style AAVSO #TELESCOPE line."));
		form->addRow(tr("Telescope"), _telescope);

		_camera = new QLineEdit(page);
		_camera->setPlaceholderText(tr("e.g. ASI2600MM"));
		_camera->setToolTip(tr(
			"Camera description used by the BAA-style AAVSO #CAMERA line. "
			"If left blank the FITS INSTRUME field is used as a fallback."));
		form->addRow(tr("Camera"), _camera);

		tabs->addTab(page, tr("Site"));
	}

	// --- Files tab ---
	{
		auto* page = new QWidget(tabs);
		auto* col  = new QVBoxLayout(page);

		_clearRecentsButton = new QPushButton(tr("Clear recent-files list"), page);
		col->addWidget(_clearRecentsButton);
		connect(_clearRecentsButton, &QPushButton::clicked,
		        this, &PreferencesDialog::clearRecents);

		_resetLayoutButton = new QPushButton(
			tr("Reset window geometry and layout"), page);
		col->addWidget(_resetLayoutButton);
		connect(_resetLayoutButton, &QPushButton::clicked,
		        this, &PreferencesDialog::resetWindowLayout);

		auto* note = new QLabel(tr(
			"Window-geometry changes take effect the next time ASTAP launches."),
			page);
		note->setWordWrap(true);
		note->setStyleSheet("color: gray;");
		col->addWidget(note);

		col->addStretch(1);
		tabs->addTab(page, tr("Files"));
	}

	auto* buttons = new QDialogButtonBox(
		QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
	connect(buttons, &QDialogButtonBox::accepted, this, &PreferencesDialog::accept);
	connect(buttons, &QDialogButtonBox::rejected, this, &PreferencesDialog::reject);
	root->addWidget(buttons);
}

/// MARK: - Settings I/O

void PreferencesDialog::loadFromSettings() {
	QSettings s;
	_databasePath->setText(s.value(kKeyDatabasePath).toString());
	const auto name = s.value(kKeyDatabaseName).toString();
	if (!name.isEmpty()) {
		const auto idx = _databaseName->findText(name);
		if (idx >= 0) {
			_databaseName->setCurrentIndex(idx);
		} else {
			_databaseName->setEditText(name);
		}
	}
	_latitude   ->setText(s.value(kKeySqmLatitude).toString());
	_longitude  ->setText(s.value(kKeySqmLongitude).toString());
	_elevation  ->setText(s.value(kKeySiteElevation).toString());
	_temperature->setValue(s.value(kKeySqmTemperature, 10.0).toDouble());
	_pressure   ->setValue(s.value(kKeySqmPressure,  1010.0).toDouble());
	_telescope  ->setText(s.value(kKeySiteTelescope).toString());
	_camera     ->setText(s.value(kKeySiteCamera).toString());
}

void PreferencesDialog::accept() {
	QSettings s;
	s.setValue(kKeyDatabasePath, _databasePath->text());
	s.setValue(kKeyDatabaseName, _databaseName->currentText());
	s.setValue(kKeySqmLatitude,    _latitude->text());
	s.setValue(kKeySqmLongitude,   _longitude->text());
	s.setValue(kKeySiteElevation,  _elevation->text());
	s.setValue(kKeySqmTemperature, _temperature->value());
	s.setValue(kKeySqmPressure,    _pressure->value());
	s.setValue(kKeySiteTelescope,  _telescope->text());
	s.setValue(kKeySiteCamera,     _camera->text());

	// Propagate catalog path to the engine so the next operation uses it
	// without restart.
	const auto dbPath = _databasePath->text();
	if (!dbPath.isEmpty()) {
		astap::reference::database_path =
			std::filesystem::path(dbPath.toStdString());
	}

	QDialog::accept();
}

/// MARK: - Actions

void PreferencesDialog::browseDatabase() {
	const auto initial = _databasePath->text();
	const auto dir = QFileDialog::getExistingDirectory(
		this, tr("Star database directory"), initial);
	if (!dir.isEmpty()) {
		_databasePath->setText(dir);
	}
}

void PreferencesDialog::clearRecents() {
	const auto response = QMessageBox::question(
		this, tr("Clear recent files"),
		tr("Clear the recent-files list?"),
		QMessageBox::Yes | QMessageBox::No);
	if (response != QMessageBox::Yes) {
		return;
	}
	QSettings s;
	s.remove(kKeyRecentFiles);
	QMessageBox::information(this, tr("Recent files"),
		tr("Recent-files list cleared. The File → Open Recent menu will update "
		   "the next time the main window's menu is rebuilt."));
}

void PreferencesDialog::resetWindowLayout() {
	const auto response = QMessageBox::question(
		this, tr("Reset window layout"),
		tr("Reset window size, position, and panel state to defaults on next launch?"),
		QMessageBox::Yes | QMessageBox::No);
	if (response != QMessageBox::Yes) {
		return;
	}
	QSettings s;
	s.remove(kKeyGeometry);
	s.remove(kKeyWindowState);
	s.remove(kKeySplitter);
	QMessageBox::information(this, tr("Window layout"),
		tr("Window layout will reset on the next launch."));
}

} // namespace astap::gui
