///----------------------------------------
///      @file focus_dialog.cpp
///   @ingroup ASTAP++
///     @brief FocusDialog implementation.
///    @author Created by John Stephen on 4/21/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "focus_dialog.h"

#include "../../src/solving/focus_fit.h"

#include <QDialogButtonBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFontDatabase>
#include <QFormLayout>
#include <QGuiApplication>
#include <QHBoxLayout>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QSettings>
#include <QTextEdit>
#include <QVBoxLayout>

#include <filesystem>
#include <vector>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

/// MARK: - Construction

FocusDialog::FocusDialog(QWidget* parent) :
	QDialog(parent) {

	setWindowTitle(tr("Focuser V-Curve Fit"));
	setModal(false);
	resize(560, 560);
	buildLayout();
}

/// MARK: - UI

void FocusDialog::buildLayout() {
	auto* root = new QVBoxLayout(this);

	auto* intro = new QLabel(tr(
		"Add 3+ FITS frames taken across a focuser sweep. Each must have a "
		"FOCUSPOS keyword in its header. The dialog measures median HFD per "
		"frame and fits a hyperbola to predict the best-focus position."), this);
	intro->setWordWrap(true);
	root->addWidget(intro);

	_fileList = new QListWidget(this);
	_fileList->setSelectionMode(QAbstractItemView::ExtendedSelection);
	_fileList->setWordWrap(false);
	root->addWidget(_fileList, 1);

	auto* fileRow = new QHBoxLayout();
	_addButton    = new QPushButton(tr("Add…"),    this);
	_removeButton = new QPushButton(tr("Remove"),  this);
	_clearButton  = new QPushButton(tr("Clear"),   this);
	fileRow->addWidget(_addButton);
	fileRow->addWidget(_removeButton);
	fileRow->addWidget(_clearButton);
	fileRow->addStretch(1);
	_runButton = new QPushButton(tr("Fit"), this);
	_runButton->setDefault(true);
	fileRow->addWidget(_runButton);
	root->addLayout(fileRow);

	connect(_addButton,    &QPushButton::clicked, this, &FocusDialog::addFiles);
	connect(_removeButton, &QPushButton::clicked, this, &FocusDialog::removeSelected);
	connect(_clearButton,  &QPushButton::clicked, this, &FocusDialog::clearList);
	connect(_runButton,    &QPushButton::clicked, this, &FocusDialog::runFit);

	// Result summary.
	auto* form = new QFormLayout();
	const auto valueFont = QFontDatabase::systemFont(QFontDatabase::FixedFont);
	auto make_value = [&](const QString& placeholder) {
		auto* lbl = new QLabel(placeholder, this);
		lbl->setFont(valueFont);
		lbl->setTextInteractionFlags(Qt::TextSelectableByMouse);
		return lbl;
	};
	_bestFocusLabel = make_value(tr("—"));
	_residualLabel  = make_value(tr("—"));
	_samplesLabel   = make_value(tr("—"));
	form->addRow(tr("Best focus (steps)"), _bestFocusLabel);
	form->addRow(tr("Mean residual"),      _residualLabel);
	form->addRow(tr("Samples used"),       _samplesLabel);
	root->addLayout(form);

	// Per-frame log.
	_log = new QTextEdit(this);
	_log->setReadOnly(true);
	_log->setFont(valueFont);
	_log->setMinimumHeight(120);
	root->addWidget(_log, 1);

	auto* buttons = new QDialogButtonBox(QDialogButtonBox::Close, this);
	connect(buttons, &QDialogButtonBox::rejected, this, &FocusDialog::reject);
	root->addWidget(buttons);
}

/// MARK: - File list

void FocusDialog::addFiles() {
	QSettings settings;
	const auto lastDir = settings.value("files/lastFocusDir").toString();
	const auto paths = QFileDialog::getOpenFileNames(
		this, tr("Focus sweep frames"), lastDir,
		tr("FITS images (*.fit *.fits *.fts);;All files (*)"));
	if (paths.isEmpty()) {
		return;
	}
	settings.setValue("files/lastFocusDir",
		QFileInfo(paths.first()).absolutePath());

	// Skip duplicates by comparing canonical paths (Qt::UserRole).
	auto existing = QSet<QString>{};
	for (auto i = 0; i < _fileList->count(); ++i) {
		existing.insert(_fileList->item(i)->data(Qt::UserRole).toString());
	}
	for (const auto& p : paths) {
		if (existing.contains(p)) {
			continue;
		}
		auto* item = new QListWidgetItem(QFileInfo(p).fileName(), _fileList);
		item->setData(Qt::UserRole, p);
		item->setToolTip(p);
	}
}

void FocusDialog::removeSelected() {
	const auto items = _fileList->selectedItems();
	for (auto* item : items) {
		delete _fileList->takeItem(_fileList->row(item));
	}
}

void FocusDialog::clearList() {
	_fileList->clear();
}

/// MARK: - Fit

void FocusDialog::runFit() {
	const auto n = _fileList->count();
	if (n < 3) {
		_log->setPlainText(tr("Need at least 3 frames for a hyperbola fit."));
		return;
	}

	auto paths = std::vector<std::filesystem::path>{};
	paths.reserve(n);
	for (auto i = 0; i < n; ++i) {
		paths.emplace_back(_fileList->item(i)->data(Qt::UserRole)
			.toString().toStdString());
	}

	auto samples = std::vector<astap::solving::FocusFitSample>(paths.size());

	_runButton->setEnabled(false);
	QGuiApplication::setOverrideCursor(Qt::WaitCursor);
	const auto result = astap::solving::fit_focus_hyperbola(
		std::span<const std::filesystem::path>{paths},
		samples.data());
	QGuiApplication::restoreOverrideCursor();
	_runButton->setEnabled(true);

	// Per-frame listing first, so the log always shows something.
	QString log;
	for (auto i = 0; i < n; ++i) {
		const auto& s = samples[i];
		const auto name = _fileList->item(i)->text();
		if (s.hfd <= 0.0 || s.hfd >= 98.0) {
			log += QString("  %1: HFD=— (no stars)\n").arg(name);
		} else if (s.position == 0.0) {
			log += QString("  %1: FOCUSPOS missing, HFD=%2\n")
				.arg(name).arg(s.hfd, 0, 'f', 2);
		} else {
			log += QString("  %1: pos=%2  HFD=%3  stars=%4\n")
				.arg(name)
				.arg(s.position, 0, 'f', 0)
				.arg(s.hfd,      0, 'f', 2)
				.arg(s.stars);
		}
	}
	if (!result.ok) {
		log += QString("\nFit: FAILED — %1").arg(QString::fromStdString(result.message));
	}
	_log->setPlainText(log);

	if (!result.ok) {
		_bestFocusLabel->setText(tr("— (fit failed)"));
		_residualLabel ->setText(tr("—"));
		_samplesLabel  ->setText(tr("—"));
		return;
	}

	_bestFocusLabel->setText(QString::number(result.focus_best, 'f', 1));
	_residualLabel ->setText(QString::number(result.lowest_error, 'f', 4)
		+ tr(" (relative)"));
	_samplesLabel  ->setText(QString::number(result.samples_used));
}

} // namespace astap::gui
