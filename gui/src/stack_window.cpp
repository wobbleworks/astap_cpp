///----------------------------------------
///      @file stack_window.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the stacking window.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "stack_window.h"
#include "image_viewer.h"

#include "../../src/core/globals.h"
#include "../../src/core/image_io.h"
#include "../../src/stacking/stack.h"
#include "../../src/stacking/stack_routines.h"

#include <QComboBox>
#include <QCoreApplication>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QProgressBar>
#include <QPushButton>
#include <QSettings>
#include <QSpinBox>
#include <QTabWidget>
#include <QVBoxLayout>

#include <cmath>
#include <filesystem>
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
}

void StackWindow::buildLightsTab() {
	auto* page = new QWidget(_tabs);
	auto* layout = new QVBoxLayout(page);

	_fileList = new QListWidget(page);
	_fileList->setSelectionMode(QAbstractItemView::ExtendedSelection);
	_fileList->setDragDropMode(QAbstractItemView::NoDragDrop);
	layout->addWidget(_fileList, 1);

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

	for (const auto& p : paths) {
		const auto items = _fileList->findItems(p, Qt::MatchExactly);
		if (items.isEmpty()) {
			_fileList->addItem(p);
		}
	}
}

void StackWindow::removeSelected() {
	const auto selected = _fileList->selectedItems();
	for (auto* item : selected) {
		delete _fileList->takeItem(_fileList->row(item));
	}
}

void StackWindow::clearList() {
	_fileList->clear();
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
}

void StackWindow::clearMasterDark() {
	astap::stacking::MasterFrameInfo info;
	(void)astap::stacking::set_master_dark({}, info);
	_darkPath->clear();
	_darkStatus->setText(tr("—"));
}

void StackWindow::clearMasterFlat() {
	astap::stacking::MasterFrameInfo info;
	(void)astap::stacking::set_master_flat({}, info);
	_flatPath->clear();
	_flatStatus->setText(tr("—"));
}

void StackWindow::startStack() {
	const auto count = _fileList->count();
	if (count < 2) {
		QMessageBox::information(this, tr("Stack"),
			tr("Add at least 2 light frames to stack."));
		return;
	}

	applySettingsToEngine();

	_stackButton->setEnabled(false);
	_progress->setValue(0);

	auto files = std::vector<astap::FileToDo>{};
	files.reserve(count);
	for (int i = 0; i < count; ++i) {
		files.push_back({_fileList->item(i)->text().toStdString(), i});
	}

	astap::stacking::set_progress_sink(
		[this](double value, const std::string& /*label*/) {
			_progress->setValue(static_cast<int>(std::round(value)));
			QCoreApplication::processEvents();
		});
	astap::stacking::set_memo2_sink([](const std::string& /*msg*/) {
		// TODO: route to a Stack log pane; swallow for now.
	});

	auto counter = 0;
	const auto method = _methodCombo->currentData().toInt();
	const auto span = std::span<astap::FileToDo>(files);
	const auto osc = 0;  // TODO: OSC/Bayer toggle

	switch (method) {
	case kMethodSigmaClip:
		astap::stacking::stack_sigmaclip(osc, span, counter);
		break;
	case kMethodLRGB:
		astap::stacking::stack_LRGB(span, counter);
		break;
	case kMethodAverage:
	default:
		astap::stacking::stack_average(osc, span, counter);
		break;
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
