///----------------------------------------
///      @file stack_window.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the stacking window.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "stack_window.h"
#include "image_viewer.h"

#include "../../src/core/fits.h"
#include "../../src/core/globals.h"
#include "../../src/core/image_io.h"
#include "../../src/stacking/stack.h"
#include "../../src/stacking/stack_routines.h"

#include <QComboBox>
#include <QCoreApplication>
#include <QDebug>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QProgressBar>
#include <QPushButton>
#include <QSettings>
#include <QSpinBox>
#include <QRegularExpression>
#include <QSet>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QTabWidget>
#include <QVBoxLayout>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <functional>
#include <vector>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

namespace {

enum MethodIndex { kMethodAverage = 0, kMethodSigmaClip = 1, kMethodLRGB = 2 };
enum AlignmentIndex {
	kAlignStarMatch = 0,
	kAlignAstrometric = 1,
	kAlignManual = 2,
	kAlignNone = 3
};

// Channel role for LRGB stacking. "Light" is untagged / treated as a mono
// frame by Average and Sigma-clip methods.
enum Channel {
	kChanLight = 0,
	kChanL = 1,
	kChanR = 2,
	kChanG = 3,
	kChanB = 4,
	kChanRGB = 5,
};

QStringList channelLabels() {
	return {
		QObject::tr("Light"),
		QStringLiteral("L"),
		QStringLiteral("R"),
		QStringLiteral("G"),
		QStringLiteral("B"),
		QStringLiteral("RGB"),
	};
}

// Read FITS header only (skip pixel data) to get width × height. Returns
// area = w * h, or 0 on failure.
long long probe_area(const QString& path) {
	astap::Header h;
	astap::ImageArray img;
	auto memo = std::vector<std::string>{};
	if (!astap::core::load_fits(std::filesystem::path(path.toStdString()),
	                            /*light=*/true, /*load_data=*/false,
	                            /*update_memo=*/false, /*get_ext=*/0,
	                            memo, h, img)) {
		return 0;
	}
	return static_cast<long long>(h.width) * h.height;
}

// Guess a channel from a filename. Matches "_L_", "_R_", "_G_", "_B_",
// "_RGB_", "_Lum_", etc. Falls back to Light.
Channel guess_channel(const QString& path) {
	const auto name = QFileInfo(path).completeBaseName().toLower();
	// Check RGB first so "_rgb_" doesn't match as just R.
	const auto tokens = name.split(QRegularExpression("[^a-z0-9]+"),
	                               Qt::SkipEmptyParts);
	for (const auto& t : tokens) {
		if (t == "rgb") return kChanRGB;
		if (t == "l" || t == "lum" || t == "luminance") return kChanL;
		if (t == "r" || t == "red")   return kChanR;
		if (t == "g" || t == "green") return kChanG;
		if (t == "b" || t == "blue")  return kChanB;
	}
	return kChanLight;
}

}  // namespace

StackWindow::StackWindow(QWidget* parent) :
	QWidget(parent, Qt::Window) {

	setWindowTitle(tr("Stack"));
	resize(520, 600);

	auto* root = new QVBoxLayout(this);

	_tabs = new QTabWidget(this);
	root->addWidget(_tabs, 1);

	buildLightsTab();
	buildCalibrationTab();
	buildSettingsTab();

	_progress = new QProgressBar(this);
	_progress->setRange(0, 100);
	_progress->setValue(0);
	_progress->setTextVisible(true);
	root->addWidget(_progress);

	auto* stackRow = new QHBoxLayout();
	_stackButton = new QPushButton(tr("Stack"), this);
	_stackButton->setDefault(true);
	stackRow->addStretch(1);
	stackRow->addWidget(_stackButton);
	root->addLayout(stackRow);

	connect(_stackButton, &QPushButton::clicked, this, &StackWindow::startStack);

	hydrateCalibrationFromSettings();
}

void StackWindow::hydrateCalibrationFromSettings() {
	QSettings settings;
	const auto darkPath = settings.value("calibration/masterDark").toString();
	if (!darkPath.isEmpty() && QFileInfo::exists(darkPath)) {
		astap::stacking::MasterFrameInfo info;
		if (astap::stacking::set_master_dark(
		        std::filesystem::path(darkPath.toStdString()), info)) {
			_darkPath->setText(darkPath);
			_darkStatus->setText(tr("Loaded: %1 × %2, exp %3s")
				.arg(info.width).arg(info.height)
				.arg(info.exposure, 0, 'f', 1));
		} else {
			settings.remove("calibration/masterDark");
		}
	}
	const auto flatPath = settings.value("calibration/masterFlat").toString();
	if (!flatPath.isEmpty() && QFileInfo::exists(flatPath)) {
		astap::stacking::MasterFrameInfo info;
		if (astap::stacking::set_master_flat(
		        std::filesystem::path(flatPath.toStdString()), info)) {
			_flatPath->setText(flatPath);
			_flatStatus->setText(tr("Loaded: %1 × %2")
				.arg(info.width).arg(info.height));
		} else {
			settings.remove("calibration/masterFlat");
		}
	}
}

void StackWindow::buildLightsTab() {
	auto* page = new QWidget(_tabs);
	auto* layout = new QVBoxLayout(page);

	_fileTable = new QTableWidget(0, 2, page);
	_fileTable->setHorizontalHeaderLabels({tr("File"), tr("Channel")});
	_fileTable->horizontalHeader()->setSectionResizeMode(
		0, QHeaderView::Stretch);
	_fileTable->horizontalHeader()->setSectionResizeMode(
		1, QHeaderView::ResizeToContents);
	_fileTable->verticalHeader()->setVisible(false);
	_fileTable->setSelectionBehavior(QAbstractItemView::SelectRows);
	_fileTable->setSelectionMode(QAbstractItemView::ExtendedSelection);
	_fileTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
	layout->addWidget(_fileTable, 1);

	auto* buttonRow = new QHBoxLayout();
	_addButton = new QPushButton(tr("Add…"), page);
	_removeButton = new QPushButton(tr("Remove"), page);
	_clearButton = new QPushButton(tr("Clear"), page);
	buttonRow->addWidget(_addButton);
	buttonRow->addWidget(_removeButton);
	buttonRow->addWidget(_clearButton);
	buttonRow->addStretch(1);
	layout->addLayout(buttonRow);

	connect(_addButton, &QPushButton::clicked, this, &StackWindow::addFiles);
	connect(_removeButton, &QPushButton::clicked, this, &StackWindow::removeSelected);
	connect(_clearButton, &QPushButton::clicked, this, &StackWindow::clearList);

	_tabs->addTab(page, tr("Lights"));
}

void StackWindow::buildCalibrationTab() {
	auto* page = new QWidget(_tabs);
	auto* layout = new QVBoxLayout(page);

	auto makeRow = [&](const QString& title,
	                   QLineEdit*& pathEdit, QLabel*& statusLabel,
	                   QPushButton*& browseBtn, QPushButton*& clearBtn) {
		auto* group = new QGroupBox(title, page);
		auto* groupLayout = new QVBoxLayout(group);
		auto* row = new QHBoxLayout();
		pathEdit = new QLineEdit(group);
		pathEdit->setReadOnly(true);
		pathEdit->setPlaceholderText(tr("No file loaded"));
		browseBtn = new QPushButton(tr("Browse…"), group);
		clearBtn = new QPushButton(tr("Clear"), group);
		row->addWidget(pathEdit, 1);
		row->addWidget(browseBtn);
		row->addWidget(clearBtn);
		groupLayout->addLayout(row);
		statusLabel = new QLabel(tr("—"), group);
		statusLabel->setStyleSheet("color: gray;");
		groupLayout->addWidget(statusLabel);
		layout->addWidget(group);
	};

	makeRow(tr("Master Dark"), _darkPath, _darkStatus,
	        _darkBrowseButton, _darkClearButton);
	makeRow(tr("Master Flat"), _flatPath, _flatStatus,
	        _flatBrowseButton, _flatClearButton);

	auto* biasNote = new QLabel(
		tr("Bias subtraction is not yet implemented; load a bias-subtracted "
		   "master dark instead."), page);
	biasNote->setWordWrap(true);
	biasNote->setStyleSheet("color: gray; font-style: italic;");
	layout->addWidget(biasNote);
	layout->addStretch(1);

	connect(_darkBrowseButton, &QPushButton::clicked,
	        this, &StackWindow::browseMasterDark);
	connect(_darkClearButton, &QPushButton::clicked,
	        this, &StackWindow::clearMasterDark);
	connect(_flatBrowseButton, &QPushButton::clicked,
	        this, &StackWindow::browseMasterFlat);
	connect(_flatClearButton, &QPushButton::clicked,
	        this, &StackWindow::clearMasterFlat);

	_tabs->addTab(page, tr("Calibration"));
}

void StackWindow::buildSettingsTab() {
	auto* page = new QWidget(_tabs);
	auto* form = new QFormLayout(page);

	_methodCombo = new QComboBox(page);
	_methodCombo->addItem(tr("Average (weighted)"), kMethodAverage);
	_methodCombo->addItem(tr("Sigma-clip average"), kMethodSigmaClip);
	_methodCombo->addItem(tr("LRGB combine"), kMethodLRGB);
	form->addRow(tr("Method:"), _methodCombo);

	_alignmentCombo = new QComboBox(page);
	_alignmentCombo->addItem(tr("Star-match (quads)"), kAlignStarMatch);
	_alignmentCombo->addItem(tr("Astrometric (per-frame solve)"), kAlignAstrometric);
	_alignmentCombo->addItem(tr("Manual reference star"), kAlignManual);
	_alignmentCombo->addItem(tr("None (pre-aligned)"), kAlignNone);
	form->addRow(tr("Alignment:"), _alignmentCombo);

	_sigmaFactor = new QDoubleSpinBox(page);
	_sigmaFactor->setRange(0.5, 10.0);
	_sigmaFactor->setSingleStep(0.1);
	_sigmaFactor->setValue(2.0);
	_sigmaFactor->setDecimals(1);
	_sigmaFactor->setSuffix(tr(" σ"));
	form->addRow(tr("Sigma-clip threshold:"), _sigmaFactor);

	_maxStars = new QSpinBox(page);
	_maxStars->setRange(50, 5000);
	_maxStars->setSingleStep(50);
	_maxStars->setValue(500);
	form->addRow(tr("Max stars for matching:"), _maxStars);

	_tabs->addTab(page, tr("Settings"));
}

void StackWindow::applySettingsToEngine() {
	const auto align = _alignmentCombo->currentData().toInt();
	astap::use_manual_align        = (align == kAlignManual);
	astap::use_ephemeris_alignment = false;  // no ephemeris UI yet
	astap::use_astrometry_internal = (align == kAlignAstrometric);
	astap::skip_alignment          = (align == kAlignNone);

	astap::sigma_clip_factor = _sigmaFactor->value();
	astap::max_stars_setting = _maxStars->value();
}

void StackWindow::addFiles() {
	QSettings settings;
	const auto lastDir = settings.value("files/lastStackDir").toString();

	const auto paths = QFileDialog::getOpenFileNames(
		this,
		tr("Add light frames"),
		lastDir,
		tr("FITS images (*.fit *.fits *.fts *.new);;"
		   "All images (*.fit *.fits *.fts *.new "
		               "*.png *.jpg *.jpeg *.bmp *.tif *.tiff);;"
		   "All files (*)"));
	if (paths.isEmpty()) {
		return;
	}

	settings.setValue("files/lastStackDir",
		QFileInfo(paths.first()).absolutePath());

	// Find already-present paths to skip duplicates.
	auto existing = QSet<QString>{};
	for (int r = 0; r < _fileTable->rowCount(); ++r) {
		existing.insert(_fileTable->item(r, 0)->text());
	}

	const auto labels = channelLabels();
	for (const auto& p : paths) {
		if (existing.contains(p)) {
			continue;
		}
		const auto row = _fileTable->rowCount();
		_fileTable->insertRow(row);
		auto* pathItem = new QTableWidgetItem(p);
		pathItem->setToolTip(p);
		_fileTable->setItem(row, 0, pathItem);

		auto* combo = new QComboBox(_fileTable);
		combo->addItems(labels);
		combo->setCurrentIndex(guess_channel(p));
		_fileTable->setCellWidget(row, 1, combo);
	}
}

void StackWindow::removeSelected() {
	auto rows = QList<int>{};
	for (const auto* item : _fileTable->selectedItems()) {
		if (item->column() == 0) {
			rows.append(item->row());
		}
	}
	std::sort(rows.begin(), rows.end(), std::greater<int>());
	for (const auto r : rows) {
		_fileTable->removeRow(r);
	}
}

void StackWindow::clearList() {
	_fileTable->setRowCount(0);
}

void StackWindow::browseMasterDark() {
	QSettings settings;
	const auto lastDir = settings.value("files/lastCalDir").toString();
	const auto path = QFileDialog::getOpenFileName(
		this, tr("Select master dark"), lastDir,
		tr("FITS images (*.fit *.fits *.fts);;All files (*)"));
	if (path.isEmpty()) {
		return;
	}
	settings.setValue("files/lastCalDir", QFileInfo(path).absolutePath());

	astap::stacking::MasterFrameInfo info;
	if (!astap::stacking::set_master_dark(
	        std::filesystem::path(path.toStdString()), info)) {
		QMessageBox::warning(this, tr("Master Dark"),
			tr("Failed to load: %1").arg(path));
		return;
	}
	_darkPath->setText(path);
	_darkStatus->setText(tr("Loaded: %1 × %2, exp %3s")
		.arg(info.width).arg(info.height)
		.arg(info.exposure, 0, 'f', 1));
	settings.setValue("calibration/masterDark", path);
}

void StackWindow::browseMasterFlat() {
	QSettings settings;
	const auto lastDir = settings.value("files/lastCalDir").toString();
	const auto path = QFileDialog::getOpenFileName(
		this, tr("Select master flat"), lastDir,
		tr("FITS images (*.fit *.fits *.fts);;All files (*)"));
	if (path.isEmpty()) {
		return;
	}
	settings.setValue("files/lastCalDir", QFileInfo(path).absolutePath());

	astap::stacking::MasterFrameInfo info;
	if (!astap::stacking::set_master_flat(
	        std::filesystem::path(path.toStdString()), info)) {
		QMessageBox::warning(this, tr("Master Flat"),
			tr("Failed to load: %1").arg(path));
		return;
	}
	_flatPath->setText(path);
	_flatStatus->setText(tr("Loaded: %1 × %2")
		.arg(info.width).arg(info.height));
	settings.setValue("calibration/masterFlat", path);
}

void StackWindow::clearMasterDark() {
	astap::stacking::MasterFrameInfo info;
	(void)astap::stacking::set_master_dark({}, info);
	_darkPath->clear();
	_darkStatus->setText(tr("—"));
	QSettings().remove("calibration/masterDark");
}

void StackWindow::clearMasterFlat() {
	astap::stacking::MasterFrameInfo info;
	(void)astap::stacking::set_master_flat({}, info);
	_flatPath->clear();
	_flatStatus->setText(tr("—"));
	QSettings().remove("calibration/masterFlat");
}

void StackWindow::startStack() {
	const auto count = _fileTable->rowCount();
	if (count < 1) {
		QMessageBox::information(this, tr("Stack"),
			tr("Add at least one frame to stack."));
		return;
	}

	applySettingsToEngine();

	_stackButton->setEnabled(false);
	_progress->setValue(0);

	// Gather (path, channel) rows from the table. Pre-probe FITS headers
	// for pixel dimensions so we can pick the largest frame as the stacking
	// reference — the engine upsamples smaller frames to match it.
	struct Row { QString path; int channel; long long area; };
	auto rows = std::vector<Row>{};
	rows.reserve(count);
	for (int i = 0; i < count; ++i) {
		const auto path = _fileTable->item(i, 0)->text();
		auto* combo = qobject_cast<QComboBox*>(
			_fileTable->cellWidget(i, 1));
		const auto chan = combo ? combo->currentIndex() : kChanLight;
		rows.push_back({path, chan, probe_area(path)});
	}

	astap::stacking::set_progress_sink(
		[this](double value, const std::string& /*label*/) {
			_progress->setValue(static_cast<int>(std::round(value)));
			QCoreApplication::processEvents();
		});
	astap::stacking::set_memo2_sink([](const std::string& msg) {
		qDebug().noquote() << "[stack]" << QString::fromStdString(msg);
	});

	auto counter = 0;
	const auto method = _methodCombo->currentData().toInt();
	const auto osc = 0;  // TODO: OSC/Bayer toggle

	if (method == kMethodLRGB) {
		// Assemble the 6-slot LRGB span: ref, R, G, B, RGB, L.
		// R, G, B are required. L and RGB are optional (engine skips empty
		// slots).
		auto firstOf = [&](Channel c) -> QString {
			for (const auto& r : rows) {
				if (r.channel == c) return r.path;
			}
			return {};
		};
		const auto lPath = firstOf(kChanL);
		const auto rPath = firstOf(kChanR);
		const auto gPath = firstOf(kChanG);
		const auto bPath = firstOf(kChanB);
		const auto rgbPath = firstOf(kChanRGB);

		auto missing = QStringList{};
		if (rPath.isEmpty()) missing << "R";
		if (gPath.isEmpty()) missing << "G";
		if (bPath.isEmpty()) missing << "B";
		if (!missing.isEmpty()) {
			QMessageBox::warning(this, tr("LRGB"),
				tr("LRGB combine needs at least one file tagged for each of "
				   "R, G, B. Missing: %1.").arg(missing.join(", ")));
			astap::stacking::set_progress_sink(nullptr);
			astap::stacking::set_memo2_sink(nullptr);
			_stackButton->setEnabled(true);
			return;
		}

		// Engine expects [ref, R, G, B, RGB, L]. Pick the ref as the
		// largest-dim tagged file so the engine upsamples smaller channels
		// to match. Falls back to L when present at equal size.
		auto areaOf = [&](const QString& p) -> long long {
			for (const auto& r : rows) {
				if (r.path == p) return r.area;
			}
			return 0;
		};
		auto refPath = !lPath.isEmpty() ? lPath : rPath;
		auto refArea = areaOf(refPath);
		for (const auto* p : {&rPath, &gPath, &bPath, &lPath}) {
			if (!p->isEmpty() && areaOf(*p) > refArea) {
				refPath = *p;
				refArea = areaOf(*p);
			}
		}
		auto files = std::vector<astap::FileToDo>{};
		files.push_back({refPath.toStdString(), 0});
		files.push_back({rPath.toStdString(), 0});
		files.push_back({gPath.toStdString(), 0});
		files.push_back({bPath.toStdString(), 0});
		files.push_back({rgbPath.toStdString(), 0});  // may be empty — engine skips
		files.push_back({lPath.toStdString(), 0});    // may be empty — engine skips

		astap::stacking::stack_LRGB(std::span<astap::FileToDo>(files), counter);
	} else {
		// Average / Sigma-clip: feed everything untagged by channel. Put
		// the largest-dim frame first so it becomes the engine's reference.
		auto ordered = rows;
		std::stable_sort(ordered.begin(), ordered.end(),
			[](const Row& a, const Row& b) { return a.area > b.area; });
		auto files = std::vector<astap::FileToDo>{};
		files.reserve(ordered.size());
		for (int i = 0; i < static_cast<int>(ordered.size()); ++i) {
			files.push_back({ordered[i].path.toStdString(), i});
		}
		const auto span = std::span<astap::FileToDo>(files);
		if (method == kMethodSigmaClip) {
			astap::stacking::stack_sigmaclip(osc, span, counter);
		} else {
			astap::stacking::stack_average(osc, span, counter);
		}
	}

	astap::stacking::set_progress_sink(nullptr);
	astap::stacking::set_memo2_sink(nullptr);

	if (counter < 2 || astap::img_loaded.empty()) {
		QMessageBox::warning(this, tr("Stack"),
			tr("Could not stack enough frames (%1 of %2).")
				.arg(counter).arg(count));
		_progress->setValue(0);
		_stackButton->setEnabled(true);
		return;
	}

	astap::head.light_count = counter;
	astap::filename2 = "Stacked_" + std::to_string(counter) + "_frames";

	if (_viewer) {
		_viewer->setImage(astap::img_loaded, astap::head);
	}

	_progress->setValue(100);
	_stackButton->setEnabled(true);

	emit stackCompleted(counter);
}

} // namespace astap::gui
