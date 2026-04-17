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
#include "../../src/stacking/live_stacking.h"

#include <QCloseEvent>
#include <QFileDialog>
#include <QFileInfo>
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

LiveStackWorker::LiveStackWorker(QString watchDir, QObject* parent) :
	QThread(parent),
	_watchDir(std::move(watchDir)) {
}

LiveStackWorker::~LiveStackWorker() {
	requestStop();
	wait();
}

void LiveStackWorker::requestStop() {
	// The engine loop polls astap::esc_pressed and exits when set. Setting
	// it here is safe to do repeatedly.
	astap::esc_pressed.store(true);
}

void LiveStackWorker::run() {
	_session = std::make_unique<astap::stacking::LiveStackSession>(
		std::filesystem::path(_watchDir.toStdString()));

	_session->set_message_hook([this](const std::string& msg) {
		emit message(QString::fromStdString(msg));
	});
	_session->set_frame_added_hook(
		[this](int accepted, int rejected, int total) {
			emit frameProcessed(accepted, rejected, total);
		});

	_session->run();
	_session.reset();
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

	_worker = new LiveStackWorker(folder, this);
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
	if (!running) {
		_pauseButton->setText(tr("Pause"));
	}
}

} // namespace astap::gui
