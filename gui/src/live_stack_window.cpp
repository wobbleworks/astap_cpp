///----------------------------------------
///      @file live_stack_window.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the Live Stacking window.
///    @author Created by John Stephen on 4/17/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "live_stack_window.h"
#include "image_viewer.h"

#include "../../src/core/globals.h"
#include "../../src/stacking/live_monitoring.h"
#include "../../src/stacking/live_stacking.h"
#include "../../src/stacking/stack.h"

#include <QCloseEvent>
#include <QComboBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSettings>
#include <QVBoxLayout>

#include <filesystem>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

// MARK: - LiveStackWorker

LiveStackWorker::LiveStackWorker(QString watchDir, Mode mode, QObject* parent) :
	QThread(parent),
	_watchDir(std::move(watchDir)),
	_mode(mode) {
}

LiveStackWorker::~LiveStackWorker() {
	requestStop();
	wait();
}

void LiveStackWorker::requestStop() {
	// Both session types poll astap::esc_pressed and exit when set.
	// Setting it here is safe to do repeatedly.
	astap::esc_pressed.store(true);
}

void LiveStackWorker::run() {
	const auto dir = std::filesystem::path(_watchDir.toStdString());

	auto messageHook = [this](const std::string& msg) {
		emit message(QString::fromStdString(msg));
	};

	if (_mode == Mode::MonitorOnly) {
		_monitorSession = std::make_unique<astap::stacking::LiveMonitorSession>(dir);
		_monitorSession->set_message_hook(messageHook);
		_monitorSession->set_frame_loaded_hook([this](int total) {
			// Map to the shared signal: accepted == total, rejected == 0.
			// The window's onWorkerFrame path refreshes the viewer when
			// `accepted > 0`, which is what we want here too.
			emit frameProcessed(total, 0, total);
		});
		_monitorSession->run();
		_monitorSession.reset();
	} else {
		_stackSession = std::make_unique<astap::stacking::LiveStackSession>(dir);
		_stackSession->set_message_hook(messageHook);
		_stackSession->set_frame_added_hook(
			[this](int accepted, int rejected, int total) {
				emit frameProcessed(accepted, rejected, total);
			});
		_stackSession->run();
		_stackSession.reset();
	}
}

// MARK: - LiveStackWindow

LiveStackWindow::LiveStackWindow(QWidget* parent) :
	QWidget(parent, Qt::Window) {

	setWindowTitle(tr("Live Stacking"));
	resize(560, 480);

	auto* root = new QVBoxLayout(this);

	// Folder row
	auto* folderRow = new QHBoxLayout();
	folderRow->addWidget(new QLabel(tr("Watch folder:"), this));
	_folderEdit = new QLineEdit(this);
	_folderEdit->setPlaceholderText(tr("Pick a folder to watch for incoming frames"));
	{
		QSettings settings;
		_folderEdit->setText(settings.value("liveStack/lastDir").toString());
	}
	_browseButton = new QPushButton(tr("Browse…"), this);
	folderRow->addWidget(_folderEdit, 1);
	folderRow->addWidget(_browseButton);
	root->addLayout(folderRow);

	// Mode selector: live-stack (align + average) vs. monitor-only (just
	// display each new frame as it lands).
	auto* modeRow = new QHBoxLayout();
	modeRow->addWidget(new QLabel(tr("Mode:"), this));
	_modeCombo = new QComboBox(this);
	_modeCombo->addItem(tr("Live stack (align + average)"),
		static_cast<int>(LiveStackWorker::Mode::LiveStack));
	_modeCombo->addItem(tr("Monitor only (display each new frame)"),
		static_cast<int>(LiveStackWorker::Mode::MonitorOnly));
	{
		QSettings settings;
		const auto persisted = settings.value("liveStack/mode",
			static_cast<int>(LiveStackWorker::Mode::LiveStack)).toInt();
		const auto idx = _modeCombo->findData(persisted);
		if (idx >= 0) _modeCombo->setCurrentIndex(idx);
	}
	modeRow->addWidget(_modeCombo, 1);
	root->addLayout(modeRow);

	// Calibration group — Dark + Flat pickers. Writes the same engine
	// globals as the batch Stack window's Calibration tab, so state is
	// shared between the two windows.
	auto makeCalRow = [&](const QString& title,
	                     QLineEdit*& pathEdit, QLabel*& statusLabel,
	                     QPushButton*& browseBtn, QPushButton*& clearBtn) {
		auto* group = new QGroupBox(title, this);
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
		root->addWidget(group);
	};
	makeCalRow(tr("Master Dark"), _darkPath, _darkStatus,
	           _darkBrowse, _darkClear);
	makeCalRow(tr("Master Flat"), _flatPath, _flatStatus,
	           _flatBrowse, _flatClear);

	// Control row
	auto* controlRow = new QHBoxLayout();
	_startButton = new QPushButton(tr("Start"), this);
	_startButton->setDefault(true);
	_pauseButton = new QPushButton(tr("Pause"), this);
	_pauseButton->setEnabled(false);
	_stopButton = new QPushButton(tr("Stop"), this);
	_stopButton->setEnabled(false);
	controlRow->addWidget(_startButton);
	controlRow->addWidget(_pauseButton);
	controlRow->addWidget(_stopButton);
	controlRow->addStretch(1);
	root->addLayout(controlRow);

	// Counters
	_countersLabel = new QLabel(tr("Idle."), this);
	_countersLabel->setStyleSheet("font-weight: bold;");
	root->addWidget(_countersLabel);

	// Log
	_log = new QPlainTextEdit(this);
	_log->setReadOnly(true);
	_log->setMaximumBlockCount(500);
	root->addWidget(_log, 1);

	connect(_browseButton, &QPushButton::clicked,
	        this, &LiveStackWindow::browseFolder);
	connect(_startButton, &QPushButton::clicked,
	        this, &LiveStackWindow::startStacking);
	connect(_pauseButton, &QPushButton::clicked,
	        this, &LiveStackWindow::togglePause);
	connect(_stopButton, &QPushButton::clicked,
	        this, &LiveStackWindow::stopStacking);
	connect(_darkBrowse, &QPushButton::clicked,
	        this, &LiveStackWindow::browseMasterDark);
	connect(_darkClear, &QPushButton::clicked,
	        this, &LiveStackWindow::clearMasterDark);
	connect(_flatBrowse, &QPushButton::clicked,
	        this, &LiveStackWindow::browseMasterFlat);
	connect(_flatClear, &QPushButton::clicked,
	        this, &LiveStackWindow::clearMasterFlat);

	hydrateCalibrationFromSettings();
}

void LiveStackWindow::hydrateCalibrationFromSettings() {
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

LiveStackWindow::~LiveStackWindow() {
	if (_worker) {
		_worker->requestStop();
		_worker->wait();
	}
}

void LiveStackWindow::closeEvent(QCloseEvent* event) {
	if (_worker && _worker->isRunning()) {
		const auto answer = QMessageBox::question(this, tr("Live Stacking"),
			tr("Stop live stacking before closing?"),
			QMessageBox::Yes | QMessageBox::Cancel);
		if (answer != QMessageBox::Yes) {
			event->ignore();
			return;
		}
		stopStacking();
	}
	event->accept();
}

void LiveStackWindow::browseMasterDark() {
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

void LiveStackWindow::browseMasterFlat() {
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

void LiveStackWindow::clearMasterDark() {
	astap::stacking::MasterFrameInfo info;
	(void)astap::stacking::set_master_dark({}, info);
	_darkPath->clear();
	_darkStatus->setText(tr("—"));
	QSettings().remove("calibration/masterDark");
}

void LiveStackWindow::clearMasterFlat() {
	astap::stacking::MasterFrameInfo info;
	(void)astap::stacking::set_master_flat({}, info);
	_flatPath->clear();
	_flatStatus->setText(tr("—"));
	QSettings().remove("calibration/masterFlat");
}

void LiveStackWindow::browseFolder() {
	QSettings settings;
	const auto start = _folderEdit->text().isEmpty()
		? settings.value("liveStack/lastDir").toString()
		: _folderEdit->text();
	const auto folder = QFileDialog::getExistingDirectory(
		this, tr("Pick watch folder"), start);
	if (folder.isEmpty()) {
		return;
	}
	_folderEdit->setText(folder);
	settings.setValue("liveStack/lastDir", folder);
}

void LiveStackWindow::startStacking() {
	const auto folder = _folderEdit->text().trimmed();
	if (folder.isEmpty() || !QFileInfo(folder).isDir()) {
		QMessageBox::warning(this, tr("Live Stacking"),
			tr("Pick a valid folder first."));
		return;
	}
	if (_worker && _worker->isRunning()) {
		return;
	}

	_log->clear();
	_paused = false;
	astap::pause_pressed.store(false);

	const auto mode = static_cast<LiveStackWorker::Mode>(
		_modeCombo->currentData().toInt());
	QSettings().setValue("liveStack/mode", static_cast<int>(mode));

	_worker = new LiveStackWorker(folder, mode, this);
	connect(_worker, &LiveStackWorker::message,
	        this, &LiveStackWindow::onWorkerMessage,
	        Qt::QueuedConnection);
	connect(_worker, &LiveStackWorker::frameProcessed,
	        this, &LiveStackWindow::onWorkerFrame,
	        Qt::QueuedConnection);
	connect(_worker, &QThread::finished,
	        this, &LiveStackWindow::onWorkerFinished);
	_worker->start();
	setRunningUi(true);
}

void LiveStackWindow::togglePause() {
	_paused = !_paused;
	astap::pause_pressed.store(_paused);
	_pauseButton->setText(_paused ? tr("Resume") : tr("Pause"));
}

void LiveStackWindow::stopStacking() {
	if (!_worker) {
		return;
	}
	_worker->requestStop();
	_stopButton->setEnabled(false);
	_pauseButton->setEnabled(false);
}

void LiveStackWindow::onWorkerMessage(QString text) {
	_log->appendPlainText(text);
}

void LiveStackWindow::onWorkerFrame(int accepted, int rejected, int total) {
	_countersLabel->setText(tr("Accepted: %1 · Rejected: %2 · Total: %3")
		.arg(accepted).arg(rejected).arg(total));

	// Refresh main viewer with the running stack written to the engine
	// globals by the session. Only refresh on accepted frames — rejected
	// ones leave img_loaded unchanged.
	if (_viewer && accepted > 0 && accepted != 0 &&
	    !astap::img_loaded.empty()) {
		_viewer->setImage(astap::img_loaded, astap::head);
		emit stackUpdated(accepted);
	}
}

void LiveStackWindow::onWorkerFinished() {
	setRunningUi(false);
	if (_worker) {
		_worker->deleteLater();
		_worker = nullptr;
	}
	// Clear the engine's abort flag so other engine calls (solve / batch
	// stack) are not tripped by our leftover signal.
	astap::esc_pressed.store(false);
	astap::pause_pressed.store(false);
}

void LiveStackWindow::setRunningUi(bool running) {
	_startButton->setEnabled(!running);
	_pauseButton->setEnabled(running);
	_stopButton->setEnabled(running);
	_browseButton->setEnabled(!running);
	_folderEdit->setEnabled(!running);
	_modeCombo->setEnabled(!running);
	_darkBrowse->setEnabled(!running);
	_darkClear->setEnabled(!running);
	_flatBrowse->setEnabled(!running);
	_flatClear->setEnabled(!running);
	if (!running) {
		_pauseButton->setText(tr("Pause"));
	}
}

} // namespace astap::gui
