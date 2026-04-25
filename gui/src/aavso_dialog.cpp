///----------------------------------------
///      @file aavso_dialog.cpp
///   @ingroup ASTAP++
///     @brief Implementation of @ref AavsoDialog.
///    @author Created by John Stephen on 4/24/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "aavso_dialog.h"

#include "../../src/core/aavso_collect.h"
#include "../../src/core/aavso_report.h"
#include "../../src/core/fits.h"
#include "../../src/core/globals.h"
#include "../../src/core/hjd.h"
#include "../../src/core/photometry.h"
#include "../../src/core/sqm.h"
#include "../../src/core/wcs.h"
#include "../../src/stacking/stack.h"

#include <QApplication>
#include <QCheckBox>
#include <QClipboard>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QFrame>
#include <QFuture>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QMetaObject>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSettings>
#include <QStringList>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QTabWidget>
#include <QVBoxLayout>
#include <QtConcurrent>

#include <cmath>
#include <fstream>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

namespace {
constexpr auto kKeyObserver = "aavso/observer_code";
constexpr auto kKeyFilter   = "aavso/filter";
constexpr auto kKeyDelim    = "aavso/delimiter";
constexpr auto kKeyHjd      = "aavso/hjd";
constexpr auto kKeyBaa      = "aavso/baa_style";
constexpr auto kKeyEnsemble = "aavso/ensemble";
constexpr auto kKeyDeltaBv  = "aavso/delta_bv";
constexpr auto kKeySlope    = "aavso/magnitude_slope";

[[nodiscard]] const char* roleLetter(ImageViewer::PickRole r) {
	switch (r) {
		case ImageViewer::PickRole::Variable: return "V";
		case ImageViewer::PickRole::Check:    return "K";
		case ImageViewer::PickRole::Comp:     return "C";
		default:                              return "?";
	}
}

[[nodiscard]] QString roleName(ImageViewer::PickRole r) {
	switch (r) {
		case ImageViewer::PickRole::Variable: return AavsoDialog::tr("Variable");
		case ImageViewer::PickRole::Check:    return AavsoDialog::tr("Check");
		case ImageViewer::PickRole::Comp:     return AavsoDialog::tr("Comp");
		default:                              return {};
	}
}
}  // namespace

AavsoDialog::AavsoDialog(QWidget* parent, ImageViewer* viewer)
	: QDialog(parent), _viewer(viewer) {
	setWindowTitle(tr("AAVSO Photometry Report"));
	setWindowFlag(Qt::WindowContextHelpButtonHint, false);
	resize(720, 640);

	auto* root = new QVBoxLayout(this);

	// Tab widget — single-frame and multi-frame share the metadata block
	// below but use distinct primary inputs and outputs.
	_tabs = new QTabWidget(this);
	root->addWidget(_tabs, /*stretch=*/1);

	// ---- Single-frame tab ------------------------------------------------
	auto* singleTab = new QWidget(this);
	auto* singleLayout = new QVBoxLayout(singleTab);

	// Three star-pick rows.
	_rows[0].role = ImageViewer::PickRole::Variable;
	_rows[1].role = ImageViewer::PickRole::Check;
	_rows[2].role = ImageViewer::PickRole::Comp;

	auto* starsGroup = new QGroupBox(tr("Stars (click on the image)"), singleTab);
	auto* starsLayout = new QFormLayout(starsGroup);
	for (auto& r : _rows) {
		r.nameEdit = new QLineEdit(singleTab);
		r.nameEdit->setPlaceholderText(roleName(r.role) + tr(" designation (e.g. SS Cyg or 000-BCP-306)"));

		r.pickButton = new QPushButton(tr("Pick on image…"), singleTab);
		r.pickButton->setProperty("role", static_cast<int>(r.role));
		connect(r.pickButton, &QPushButton::clicked,
		        this, &AavsoDialog::onPickButtonClicked);

		r.magLabel = new QLabel(tr("—"), singleTab);
		r.magLabel->setMinimumWidth(120);

		auto* row = new QHBoxLayout{};
		row->addWidget(r.nameEdit, /*stretch=*/2);
		row->addWidget(r.pickButton);
		row->addWidget(r.magLabel, /*stretch=*/1);
		starsLayout->addRow(roleName(r.role) + ":", row);
	}
	singleLayout->addWidget(starsGroup);

	// Comp-specific extras.
	auto* compRow = new QHBoxLayout{};
	_compCatalogMag = new QDoubleSpinBox(singleTab);
	_compCatalogMag->setRange(-30.0, 30.0);
	_compCatalogMag->setDecimals(3);
	_compCatalogMag->setSingleStep(0.1);
	_compCatalogMag->setValue(0.0);
	_compCatalogMag->setSpecialValueText(tr("(none)"));
	_ensembleCheck = new QCheckBox(tr("Ensemble (CNAME=ENSEMBLE)"), singleTab);
	connect(_ensembleCheck, &QCheckBox::toggled,
	        this, &AavsoDialog::onEnsembleToggled);
	compRow->addWidget(new QLabel(tr("Comp catalog magnitude:"), singleTab));
	compRow->addWidget(_compCatalogMag, /*stretch=*/1);
	compRow->addWidget(_ensembleCheck);
	singleLayout->addLayout(compRow);
	singleLayout->addStretch(1);
	_tabs->addTab(singleTab, tr("Single frame"));

	// ---- Multi-frame tab -------------------------------------------------
	auto* multiTab = new QWidget(this);
	auto* multiLayout = new QVBoxLayout(multiTab);

	multiLayout->addWidget(new QLabel(tr(
		"Pick the variable, check, and comp stars on the Single-frame tab "
		"first — their celestial coordinates carry into multi-frame mode. "
		"Then add additional FITS files (each must already be plate-solved "
		"AND photometrically calibrated) and click Measure."), multiTab));

	auto* fileRow = new QHBoxLayout{};
	_fileList = new QListWidget(multiTab);
	_fileList->setSelectionMode(QAbstractItemView::ExtendedSelection);
	_addFilesBtn    = new QPushButton(tr("Add files…"), multiTab);
	_removeFilesBtn = new QPushButton(tr("Remove selected"), multiTab);
	_clearFilesBtn  = new QPushButton(tr("Clear"), multiTab);
	_measureAllBtn  = new QPushButton(tr("Measure all frames"), multiTab);
	connect(_addFilesBtn,    &QPushButton::clicked, this, &AavsoDialog::onAddFiles);
	connect(_removeFilesBtn, &QPushButton::clicked, this, &AavsoDialog::onRemoveSelectedFiles);
	connect(_clearFilesBtn,  &QPushButton::clicked, this, &AavsoDialog::onClearFiles);
	connect(_measureAllBtn,  &QPushButton::clicked, this, &AavsoDialog::onMeasureAllFrames);

	auto* fileButtons = new QVBoxLayout{};
	fileButtons->addWidget(_addFilesBtn);
	fileButtons->addWidget(_removeFilesBtn);
	fileButtons->addWidget(_clearFilesBtn);
	fileButtons->addStretch(1);
	fileButtons->addWidget(_measureAllBtn);

	fileRow->addWidget(_fileList, /*stretch=*/1);
	fileRow->addLayout(fileButtons);
	multiLayout->addLayout(fileRow, /*stretch=*/1);

	_resultsTable = new QTableWidget(multiTab);
	_resultsTable->setColumnCount(7);
	_resultsTable->setHorizontalHeaderLabels(
		{tr("File"), tr("JD"), tr("Var mag"), tr("Check mag"),
		 tr("Comp mag"), tr("Var SNR"), tr("Airmass")});
	_resultsTable->horizontalHeader()->setStretchLastSection(true);
	_resultsTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
	_resultsTable->setSelectionBehavior(QAbstractItemView::SelectRows);
	multiLayout->addWidget(_resultsTable, /*stretch=*/1);

	_multiLog = new QPlainTextEdit(multiTab);
	_multiLog->setReadOnly(true);
	_multiLog->setMaximumBlockCount(1000);
	_multiLog->setMaximumHeight(120);
	_multiLog->setLineWrapMode(QPlainTextEdit::NoWrap);
	multiLayout->addWidget(_multiLog);

	_multiStatus = new QLabel(multiTab);
	_multiStatus->setWordWrap(true);
	multiLayout->addWidget(_multiStatus);

	auto* multiActions = new QHBoxLayout{};
	_copyMultiBtn = new QPushButton(tr("Copy multi-row report to clipboard"), multiTab);
	_saveMultiBtn = new QPushButton(tr("Save multi-row report as…"), multiTab);
	connect(_copyMultiBtn, &QPushButton::clicked, this, &AavsoDialog::onCopyMultiReport);
	connect(_saveMultiBtn, &QPushButton::clicked, this, &AavsoDialog::onSaveMultiReport);
	multiActions->addWidget(_copyMultiBtn);
	multiActions->addWidget(_saveMultiBtn);
	multiActions->addStretch(1);
	multiLayout->addLayout(multiActions);

	_tabs->addTab(multiTab, tr("Multi-frame (time series)"));

	// ---- Shared metadata block (under the tabs) --------------------------
	auto* meta = new QFormLayout{};

	_observerCode = new QLineEdit(this);
	_observerCode->setPlaceholderText(tr("e.g. JST01"));
	meta->addRow(tr("Observer code:"), _observerCode);

	_filter = new QComboBox(this);
	for (const auto* f : {"V", "B", "R", "I", "U", "TR", "TG", "TB", "CV", "SG", "SR", "SI"}) {
		_filter->addItem(QString::fromLatin1(f));
	}
	meta->addRow(tr("Filter:"), _filter);

	_delimiter = new QComboBox(this);
	_delimiter->addItem(tr(", (comma)"), QStringLiteral(","));
	_delimiter->addItem(tr("tab"),       QStringLiteral("\t"));
	_delimiter->addItem(tr("|"),         QStringLiteral("|"));
	_delimiter->addItem(tr(";"),         QStringLiteral(";"));
	meta->addRow(tr("Delimiter:"), _delimiter);

	_hjd = new QCheckBox(tr("Use HJD instead of JD"), this);
	meta->addRow(QString{}, _hjd);

	_baaStyle = new QCheckBox(tr("BAA-style header (location / telescope / camera)"), this);
	meta->addRow(QString{}, _baaStyle);

	_deltaBv = new QDoubleSpinBox(this);
	_deltaBv->setRange(-5.0, 5.0);
	_deltaBv->setDecimals(3);
	_deltaBv->setSingleStep(0.05);
	meta->addRow(tr("Δ(B−V):"), _deltaBv);

	_magnitudeSlope = new QDoubleSpinBox(this);
	_magnitudeSlope->setRange(-5.0, 5.0);
	_magnitudeSlope->setDecimals(3);
	_magnitudeSlope->setSingleStep(0.01);
	meta->addRow(tr("Magnitude slope:"), _magnitudeSlope);

	root->addLayout(meta);

	// Status.
	_status = new QLabel(this);
	_status->setWordWrap(true);
	_status->setFrameShape(QFrame::StyledPanel);
	_status->setMinimumHeight(40);
	root->addWidget(_status);

	// Buttons. The single-frame report buttons live here (multi-frame has
	// its own buttons inside the tab).
	auto* buttons = new QDialogButtonBox(this);
	_generateBtn = buttons->addButton(tr("Copy single-frame report to clipboard"),
		QDialogButtonBox::AcceptRole);
	_saveAsBtn = buttons->addButton(tr("Save single-frame report as…"),
		QDialogButtonBox::ActionRole);
	auto* closeBtn = buttons->addButton(QDialogButtonBox::Close);
	connect(_generateBtn, &QPushButton::clicked, this, &AavsoDialog::onGenerate);
	connect(_saveAsBtn,   &QPushButton::clicked, this, &AavsoDialog::onSaveAs);
	connect(closeBtn,     &QPushButton::clicked, this, &QDialog::close);
	root->addWidget(buttons);

	// Hook up the viewer's pick signal once.
	if (_viewer) {
		connect(_viewer, &ImageViewer::picked,
		        this, &AavsoDialog::onPicked);
	}

	loadSettings();
	onEnsembleToggled(_ensembleCheck->isChecked());
	setStatus(tr("Ready. Pick the variable, check, and comp stars on the image."), true);
}

AavsoDialog::~AavsoDialog() {
	saveSettings();
	if (_viewer && _viewer->pickMode() != ImageViewer::PickRole::None) {
		_viewer->setPickMode(ImageViewer::PickRole::None);
	}
}

void AavsoDialog::onPickButtonClicked() {
	auto* btn = qobject_cast<QPushButton*>(sender());
	if (!btn || !_viewer) return;

	const auto role = static_cast<ImageViewer::PickRole>(btn->property("role").toInt());
	_viewer->setPickMode(role);
	setStatus(tr("Click the %1 star on the image…").arg(roleName(role)), true);
}

void AavsoDialog::onPicked(QPointF imagePos, int roleInt) {
	const auto role = static_cast<ImageViewer::PickRole>(roleInt);
	for (auto& r : _rows) {
		if (r.role == role) {
			measureAtPick(r, imagePos);
			return;
		}
	}
}

void AavsoDialog::measureAtPick(StarRow& row, QPointF imagePos) {
	row.position = imagePos;
	_viewer->setPickMarker(row.role, imagePos);

	// Run HFD at the picked location to get a refined centroid + flux.
	auto scratch = astap::core::HfdScratch{};
	auto result  = astap::core::HfdResult{};
	const auto px = static_cast<int>(std::lround(imagePos.x() - 1.0));
	const auto py = static_cast<int>(std::lround(imagePos.y() - 1.0));
	astap::core::HFD(astap::img_loaded, px, py, /*rs=*/14,
	                  /*aperture_small=*/99.0, /*adu_e=*/0.0,
	                  /*xbinning=*/1.0, result, scratch);

	if (result.flux <= 0.0) {
		row.magLabel->setText(tr("no star found"));
		row.magLabel->setStyleSheet("color: #c0392b;");
		row.instrumentalMag = 0.0;
		row.snr = 0.0;
		setStatus(tr("Could not measure a star at that pixel — try clicking closer to a real source."), false);
		return;
	}

	if (astap::head.mzero == 0.0) {
		row.magLabel->setText(tr("flux=%1 (no MZERO)").arg(result.flux, 0, 'f', 0));
		row.magLabel->setStyleSheet("color: #b35d00;");
		row.instrumentalMag = 0.0;
		row.snr = result.snr;
		setStatus(tr("Run Image → Photometric Calibration first so MZERO is "
		             "set; then re-pick the star."), false);
		return;
	}

	const auto mag = astap::head.mzero - 2.5 * std::log10(result.flux);
	row.instrumentalMag = mag;
	row.snr = result.snr;
	// Update marker to refined centroid (HFD reports xc/yc 0-based; FITS is 1-based).
	const auto refined = QPointF(result.xc + 1.0, result.yc + 1.0);
	row.position = refined;
	_viewer->setPickMarker(row.role, refined);

	// Capture celestial coords so the multi-frame collector can re-locate
	// this physical star in every frame regardless of pointing drift.
	auto ra = 0.0, dec = 0.0;
	astap::core::pixel_to_celestial(astap::head, refined.x(), refined.y(),
	                                 /*formalism=*/0, ra, dec);
	row.ra  = ra;
	row.dec = dec;

	row.magLabel->setText(tr("mag=%1  HFD=%2  SNR=%3")
		.arg(mag, 0, 'f', 3)
		.arg(result.hfd, 0, 'f', 2)
		.arg(result.snr, 0, 'f', 1));
	row.magLabel->setStyleSheet("color: #1e7a1e;");
	setStatus(tr("%1 star measured: mag=%2.")
		.arg(roleName(row.role))
		.arg(mag, 0, 'f', 3), true);
}

void AavsoDialog::onEnsembleToggled(bool ensemble) {
	auto& compRow = _rows[2];
	compRow.nameEdit->setEnabled(!ensemble);
	compRow.pickButton->setEnabled(!ensemble);
	_compCatalogMag->setEnabled(!ensemble);
}

bool AavsoDialog::validateForReport(QString& why) const {
	if (astap::head.naxis == 0 || astap::head.cd1_1 == 0.0) {
		why = tr("Image is not plate-solved.");
		return false;
	}
	if (astap::head.mzero == 0.0) {
		why = tr("MZERO not set. Run Image → Photometric Calibration first.");
		return false;
	}
	if (_observerCode->text().trimmed().isEmpty()) {
		why = tr("Observer code is required.");
		return false;
	}
	if (_rows[0].nameEdit->text().trimmed().isEmpty()) {
		why = tr("Variable name is required.");
		return false;
	}
	if (!_rows[0].position) {
		why = tr("Variable star not picked.");
		return false;
	}
	if (!_rows[1].position) {
		why = tr("Check star not picked.");
		return false;
	}
	if (!_ensembleCheck->isChecked() && !_rows[2].position) {
		why = tr("Comp star not picked (or enable Ensemble mode).");
		return false;
	}
	if (!_ensembleCheck->isChecked() && _compCatalogMag->value() == _compCatalogMag->minimum()) {
		// Special-value text was shown — no real catalog magnitude.
		why = tr("Enter the comp star's documented (catalog) magnitude.");
		return false;
	}
	return true;
}

QString AavsoDialog::buildReport() const {
	auto m = astap::core::AavsoMeasurement{};
	m.variable_name   = astap::core::clean_abbreviation(
		_rows[0].nameEdit->text().toStdString());
	m.check_name      = astap::core::clean_abbreviation(
		_rows[1].nameEdit->text().toStdString());
	m.var_magnitude   = _rows[0].instrumentalMag;
	m.check_magnitude = _rows[1].instrumentalMag;
	m.snr             = _rows[0].snr;

	if (!_ensembleCheck->isChecked()) {
		m.comp_name        = astap::core::clean_abbreviation(
			_rows[2].nameEdit->text().toStdString());
		m.comp_magnitude   = _rows[2].instrumentalMag;
		m.comp_catalog_mag = _compCatalogMag->value();
	}

	m.filter_band = _filter->currentText().toStdString();

	// JD or HJD from the FITS header timestamp.
	astap::stacking::date_to_jd(astap::head.date_obs, astap::head.date_avg,
	                             astap::head.exposure);
	auto jd = astap::jd_mid;
	if (_hjd->isChecked() && jd > 2400000.0) {
		jd = astap::core::JD_to_HJD(jd, astap::head.ra0, astap::head.dec0);
	}
	m.jd = jd;
	m.airmass = (astap::airmass > 0.0 && astap::airmass < 99.0)
	          ? astap::airmass : 99.0;

	auto opts = astap::core::AavsoOptions{};
	opts.observer_code = _observerCode->text().trimmed().toStdString();
	opts.delimiter     = _delimiter->currentData().toString().toStdString();
	opts.hjd_date      = _hjd->isChecked();
	opts.baa_style     = _baaStyle->isChecked();
	opts.ensemble      = _ensembleCheck->isChecked();
	opts.delta_bv      = _deltaBv->value();
	opts.magnitude_slope = _magnitudeSlope->value();
	opts.software_version = std::string(astap::astap_version);

	if (opts.baa_style) {
		// Site coordinates come from the SQM globals (settable via the SQM
		// dialog or the FITS header). Site elevation, telescope name, and
		// camera name are not exposed as cross-module globals yet — leave
		// blank for now; the user can edit the saved report or set them
		// up via Preferences in a follow-up.
		opts.site_lat   = astap::core::sitelat;
		opts.site_long  = astap::core::sitelong;
		opts.telescope  = "";
		opts.camera     = astap::instrum;
	}

	const auto report = astap::core::format_aavso_report(m, opts);
	return QString::fromStdString(report);
}

void AavsoDialog::onGenerate() {
	saveSettings();

	auto why = QString{};
	if (!validateForReport(why)) {
		setStatus(why, false);
		return;
	}
	const auto text = buildReport();
	QGuiApplication::clipboard()->setText(text);
	setStatus(tr("Report copied to clipboard (%1 bytes).").arg(text.size()), true);
}

void AavsoDialog::onSaveAs() {
	saveSettings();

	auto why = QString{};
	if (!validateForReport(why)) {
		setStatus(why, false);
		return;
	}

	const auto text = buildReport();

	auto basename = QFileInfo(QString::fromStdString(astap::filename2)).completeBaseName();
	if (basename.isEmpty()) basename = QStringLiteral("aavso_report");
	const auto suggested = basename + QStringLiteral("_aavso.txt");

	const auto path = QFileDialog::getSaveFileName(this, tr("Save AAVSO report"),
		suggested, tr("Text files (*.txt);;All files (*)"));
	if (path.isEmpty()) return;

	std::ofstream out(path.toStdString(), std::ios::binary);
	if (!out) {
		setStatus(tr("Failed to write %1").arg(path), false);
		return;
	}
	out.write(text.toStdString().data(), static_cast<std::streamsize>(text.size()));
	if (!out) {
		setStatus(tr("Write error on %1").arg(path), false);
		return;
	}
	setStatus(tr("Report saved to %1").arg(path), true);
}

void AavsoDialog::setStatus(const QString& msg, bool ok) {
	_status->setText(msg);
	_status->setStyleSheet(ok ? "color: #1e7a1e;" : "color: #c0392b;");
}

void AavsoDialog::loadSettings() {
	QSettings settings;
	_observerCode->setText(settings.value(kKeyObserver).toString());
	if (const auto idx = _filter->findText(
			settings.value(kKeyFilter, "V").toString()); idx >= 0) {
		_filter->setCurrentIndex(idx);
	}
	if (const auto idx = _delimiter->findData(
			settings.value(kKeyDelim, ",").toString()); idx >= 0) {
		_delimiter->setCurrentIndex(idx);
	}
	_hjd->setChecked     (settings.value(kKeyHjd, false).toBool());
	_baaStyle->setChecked(settings.value(kKeyBaa, false).toBool());
	_ensembleCheck->setChecked(settings.value(kKeyEnsemble, true).toBool());
	_deltaBv->setValue       (settings.value(kKeyDeltaBv, 0.0).toDouble());
	_magnitudeSlope->setValue(settings.value(kKeySlope,   0.0).toDouble());
}

void AavsoDialog::saveSettings() const {
	QSettings settings;
	settings.setValue(kKeyObserver, _observerCode->text());
	settings.setValue(kKeyFilter,   _filter->currentText());
	settings.setValue(kKeyDelim,    _delimiter->currentData().toString());
	settings.setValue(kKeyHjd,      _hjd->isChecked());
	settings.setValue(kKeyBaa,      _baaStyle->isChecked());
	settings.setValue(kKeyEnsemble, _ensembleCheck->isChecked());
	settings.setValue(kKeyDeltaBv,  _deltaBv->value());
	settings.setValue(kKeySlope,    _magnitudeSlope->value());
}

// ---- Multi-frame ----------------------------------------------------------

void AavsoDialog::onAddFiles() {
	QSettings settings;
	const auto lastDir = settings.value("aavso/last_multi_dir").toString();
	const auto picked = QFileDialog::getOpenFileNames(this,
		tr("Add FITS frames"), lastDir,
		tr("FITS images (*.fit *.fits *.fts);;All files (*)"));
	if (picked.isEmpty()) return;
	settings.setValue("aavso/last_multi_dir",
		QFileInfo(picked.first()).absolutePath());
	for (const auto& p : picked) {
		_fileList->addItem(p);
	}
	_multiStatus->setText(tr("%1 file(s) queued.").arg(_fileList->count()));
}

void AavsoDialog::onRemoveSelectedFiles() {
	const auto selected = _fileList->selectedItems();
	for (auto* item : selected) {
		delete _fileList->takeItem(_fileList->row(item));
	}
}

void AavsoDialog::onClearFiles() {
	_fileList->clear();
	_resultsTable->setRowCount(0);
	_multiLog->clear();
	_measureResult.reset();
}

void AavsoDialog::onMeasureAllFrames() {
	if (_fileList->count() == 0) {
		_multiStatus->setText(tr("Add FITS files first."));
		return;
	}
	if (!_rows[0].position) {
		_multiStatus->setText(tr("Pick the variable star on the Single-frame "
		                          "tab first (its RA/Dec is reused per frame)."));
		return;
	}
	if (!_rows[1].position) {
		_multiStatus->setText(tr("Pick the check star on the Single-frame tab first."));
		return;
	}
	if (!_ensembleCheck->isChecked() && !_rows[2].position) {
		_multiStatus->setText(tr("Pick the comp star on the Single-frame tab "
		                          "first (or enable Ensemble mode)."));
		return;
	}
	if (_measureWatcher) return;  // already running

	saveSettings();
	_multiLog->clear();
	_resultsTable->setRowCount(0);
	_measureAllBtn->setEnabled(false);
	_multiStatus->setText(tr("Measuring %1 frame(s)…").arg(_fileList->count()));

	auto opts = astap::core::AavsoCollectOptions{};
	opts.variable.name = _rows[0].nameEdit->text().toStdString();
	opts.variable.ra   = _rows[0].ra;
	opts.variable.dec  = _rows[0].dec;
	opts.check.name    = _rows[1].nameEdit->text().toStdString();
	opts.check.ra      = _rows[1].ra;
	opts.check.dec     = _rows[1].dec;
	opts.ensemble      = _ensembleCheck->isChecked();
	if (!opts.ensemble) {
		opts.comp.name              = _rows[2].nameEdit->text().toStdString();
		opts.comp.ra                = _rows[2].ra;
		opts.comp.dec               = _rows[2].dec;
		opts.comp_catalog_magnitude = _compCatalogMag->value();
	}
	opts.hjd_date    = _hjd->isChecked();
	opts.filter_band = _filter->currentText().toStdString();

	auto files = std::vector<std::filesystem::path>{};
	files.reserve(_fileList->count());
	for (auto i = 0; i < _fileList->count(); ++i) {
		files.emplace_back(_fileList->item(i)->text().toStdString());
	}

	// Marshal log lines back to the GUI thread.
	opts.log = [this](const std::string& s) {
		QMetaObject::invokeMethod(this, "appendMultiLog", Qt::QueuedConnection,
			Q_ARG(QString, QString::fromStdString(s)));
	};

	_measureWatcher = std::make_unique<QFutureWatcher<astap::core::AavsoCollectResult>>();
	connect(_measureWatcher.get(),
	        &QFutureWatcher<astap::core::AavsoCollectResult>::finished,
	        this, &AavsoDialog::onMeasureFinished);
	auto future = QtConcurrent::run(
		[files = std::move(files), opts = std::move(opts)]() {
			return astap::core::collect_aavso_measurements(files, opts);
		});
	_measureWatcher->setFuture(future);
}

void AavsoDialog::onMeasureFinished() {
	auto result = _measureWatcher->result();
	_measureWatcher.reset();
	_measureAllBtn->setEnabled(true);

	_measureResult = std::make_unique<astap::core::AavsoCollectResult>(std::move(result));
	rebuildResultsTable();
	_multiStatus->setText(tr("Measured %1 of %2 frame(s); %3 skipped.")
		.arg(_measureResult->frames_accepted)
		.arg(_measureResult->frames_total)
		.arg(static_cast<int>(_measureResult->skipped_files.size())));
}

void AavsoDialog::appendMultiLog(const QString& line) {
	_multiLog->appendPlainText(line);
}

void AavsoDialog::rebuildResultsTable() {
	_resultsTable->setRowCount(0);
	if (!_measureResult) return;

	const auto& rows = _measureResult->rows;
	_resultsTable->setRowCount(static_cast<int>(rows.size()));
	for (auto i = 0; i < static_cast<int>(rows.size()); ++i) {
		const auto& m = rows[i];
		const auto file = (i < _fileList->count())
			? QFileInfo(_fileList->item(i)->text()).fileName()
			: QStringLiteral("?");
		_resultsTable->setItem(i, 0, new QTableWidgetItem(file));
		_resultsTable->setItem(i, 1, new QTableWidgetItem(
			QString::number(m.jd, 'f', 5)));
		_resultsTable->setItem(i, 2, new QTableWidgetItem(
			QString::number(m.var_magnitude, 'f', 3)));
		_resultsTable->setItem(i, 3, new QTableWidgetItem(
			QString::number(m.check_magnitude, 'f', 3)));
		_resultsTable->setItem(i, 4, new QTableWidgetItem(
			m.comp_magnitude == 0.0 ? QStringLiteral("—")
			                         : QString::number(m.comp_magnitude, 'f', 3)));
		_resultsTable->setItem(i, 5, new QTableWidgetItem(
			QString::number(m.snr, 'f', 1)));
		_resultsTable->setItem(i, 6, new QTableWidgetItem(
			m.airmass >= 99.0 ? QStringLiteral("na")
			                  : QString::number(m.airmass, 'f', 3)));
	}
	_resultsTable->resizeColumnsToContents();
}

void AavsoDialog::onCopyMultiReport() {
	if (!_measureResult || _measureResult->rows.empty()) {
		_multiStatus->setText(tr("Run Measure first."));
		return;
	}
	if (_observerCode->text().trimmed().isEmpty()) {
		_multiStatus->setText(tr("Observer code is required."));
		return;
	}
	saveSettings();

	auto opts = astap::core::AavsoOptions{};
	opts.observer_code = _observerCode->text().trimmed().toStdString();
	opts.delimiter     = _delimiter->currentData().toString().toStdString();
	opts.hjd_date      = _hjd->isChecked();
	opts.baa_style     = _baaStyle->isChecked();
	opts.ensemble      = _ensembleCheck->isChecked();
	opts.delta_bv      = _deltaBv->value();
	opts.magnitude_slope = _magnitudeSlope->value();
	opts.software_version = std::string(astap::astap_version);
	if (opts.baa_style) {
		opts.site_lat   = astap::core::sitelat;
		opts.site_long  = astap::core::sitelong;
		opts.camera     = astap::instrum;
	}

	const auto report = QString::fromStdString(
		astap::core::format_aavso_report(_measureResult->rows, opts));
	QGuiApplication::clipboard()->setText(report);
	_multiStatus->setText(tr("Multi-row report copied to clipboard "
	                          "(%1 rows, %2 bytes).")
		.arg(_measureResult->rows.size())
		.arg(report.size()));
}

void AavsoDialog::onSaveMultiReport() {
	if (!_measureResult || _measureResult->rows.empty()) {
		_multiStatus->setText(tr("Run Measure first."));
		return;
	}
	if (_observerCode->text().trimmed().isEmpty()) {
		_multiStatus->setText(tr("Observer code is required."));
		return;
	}
	saveSettings();

	auto opts = astap::core::AavsoOptions{};
	opts.observer_code = _observerCode->text().trimmed().toStdString();
	opts.delimiter     = _delimiter->currentData().toString().toStdString();
	opts.hjd_date      = _hjd->isChecked();
	opts.baa_style     = _baaStyle->isChecked();
	opts.ensemble      = _ensembleCheck->isChecked();
	opts.delta_bv      = _deltaBv->value();
	opts.magnitude_slope = _magnitudeSlope->value();
	opts.software_version = std::string(astap::astap_version);
	if (opts.baa_style) {
		opts.site_lat  = astap::core::sitelat;
		opts.site_long = astap::core::sitelong;
		opts.camera    = astap::instrum;
	}

	const auto suggested = QStringLiteral("aavso_timeseries.txt");
	const auto path = QFileDialog::getSaveFileName(this,
		tr("Save AAVSO multi-row report"),
		suggested, tr("Text files (*.txt);;All files (*)"));
	if (path.isEmpty()) return;

	const auto text = astap::core::format_aavso_report(_measureResult->rows, opts);
	std::ofstream out(path.toStdString(), std::ios::binary);
	if (!out) {
		_multiStatus->setText(tr("Failed to write %1").arg(path));
		return;
	}
	out.write(text.data(), static_cast<std::streamsize>(text.size()));
	_multiStatus->setText(tr("Saved %1 row(s) to %2")
		.arg(_measureResult->rows.size()).arg(path));
}

} // namespace astap::gui
