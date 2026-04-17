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

#include <QCoreApplication>
#include <QFileDialog>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QListWidget>
#include <QMessageBox>
#include <QProgressBar>
#include <QPushButton>
#include <QSettings>
#include <QVBoxLayout>

#include <cmath>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

StackWindow::StackWindow(QWidget* parent) :
	QWidget(parent, Qt::Window) {

	setWindowTitle(tr("Stack"));
	resize(480, 520);

	auto* root = new QVBoxLayout(this);

	// Lights group
	auto* lightsGroup = new QGroupBox(tr("Light frames"), this);
	auto* lightsLayout = new QVBoxLayout(lightsGroup);

	_fileList = new QListWidget(lightsGroup);
	_fileList->setSelectionMode(QAbstractItemView::ExtendedSelection);
	_fileList->setDragDropMode(QAbstractItemView::NoDragDrop);
	lightsLayout->addWidget(_fileList, 1);

	auto* buttonRow = new QHBoxLayout();
	_addButton = new QPushButton(tr("Add…"), lightsGroup);
	_removeButton = new QPushButton(tr("Remove"), lightsGroup);
	_clearButton = new QPushButton(tr("Clear"), lightsGroup);
	buttonRow->addWidget(_addButton);
	buttonRow->addWidget(_removeButton);
	buttonRow->addWidget(_clearButton);
	buttonRow->addStretch(1);
	lightsLayout->addLayout(buttonRow);

	root->addWidget(lightsGroup, 1);

	// Progress + Stack
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

	// Wires
	connect(_addButton, &QPushButton::clicked, this, &StackWindow::addFiles);
	connect(_removeButton, &QPushButton::clicked, this, &StackWindow::removeSelected);
	connect(_clearButton, &QPushButton::clicked, this, &StackWindow::clearList);
	connect(_stackButton, &QPushButton::clicked, this, &StackWindow::startStack);
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
		// Avoid duplicates
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

void StackWindow::startStack() {
	const auto count = _fileList->count();
	if (count < 2) {
		QMessageBox::information(this, tr("Stack"),
			tr("Add at least 2 light frames to stack."));
		return;
	}

	_stackButton->setEnabled(false);
	_progress->setValue(0);

	// Simple average stacking. Phase 5a: no calibration, no alignment.
	// Load each frame, accumulate into a double buffer, divide at the end.
	astap::ImageArray accumulator;
	astap::Header refHead{};
	auto framesLoaded = 0;

	for (int i = 0; i < count; ++i) {
		const auto path = _fileList->item(i)->text().toStdString();

		// Load frame via the engine's universal loader.
		astap::ImageArray frame;
		astap::Header frameHead{};
		auto memo = std::vector<std::string>{};
		auto pathCopy = path;
		const auto ok = astap::core::load_image(
			pathCopy, frame, frameHead, memo,
			/*re_center=*/false, /*plot=*/false);

		if (!ok || frame.empty()) {
			continue;
		}

		if (framesLoaded == 0) {
			// First frame: initialise accumulator with matching dimensions.
			refHead = frameHead;
			const auto nch = static_cast<int>(frame.size());
			const auto h = frameHead.height;
			const auto w = frameHead.width;
			accumulator.assign(nch,
				std::vector<std::vector<float>>(h,
					std::vector<float>(w, 0.0f)));
		}

		// Accumulate (must match dimensions of the first frame).
		const auto nch = std::min(accumulator.size(), frame.size());
		const auto h = std::min(accumulator[0].size(), frame[0].size());
		for (std::size_t c = 0; c < nch; ++c) {
			for (std::size_t y = 0; y < h; ++y) {
				const auto w = std::min(accumulator[c][y].size(),
					frame[c][y].size());
				for (std::size_t x = 0; x < w; ++x) {
					accumulator[c][y][x] += frame[c][y][x];
				}
			}
		}

		++framesLoaded;
		_progress->setValue(static_cast<int>(
			(i + 1) * 100.0 / count));
		QCoreApplication::processEvents();
	}

	if (framesLoaded < 2) {
		QMessageBox::warning(this, tr("Stack"),
			tr("Could not load enough frames (%1 of %2).")
				.arg(framesLoaded).arg(count));
		_stackButton->setEnabled(true);
		return;
	}

	// Divide by frame count to get average.
	const auto divisor = static_cast<float>(framesLoaded);
	for (auto& ch : accumulator) {
		for (auto& row : ch) {
			for (auto& v : row) {
				v /= divisor;
			}
		}
	}

	// Push the result into the engine globals and the viewer.
	astap::img_loaded = std::move(accumulator);
	astap::head = refHead;
	astap::head.light_count = framesLoaded;
	astap::filename2 = "Stacked_" + std::to_string(framesLoaded) + "_frames";

	if (_viewer) {
		_viewer->setImage(astap::img_loaded, astap::head);
	}

	_progress->setValue(100);
	_stackButton->setEnabled(true);

	emit stackCompleted(framesLoaded);
}

} // namespace astap::gui
