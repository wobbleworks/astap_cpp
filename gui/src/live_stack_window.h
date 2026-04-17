///----------------------------------------
///      @file live_stack_window.h
///   @ingroup ASTAP++
///     @brief Live-stacking window: watches a folder, aligns and averages
///            new frames as they land, pushes the running stack to the main
///            viewer.
///    @author Created by John Stephen on 4/17/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <QThread>
#include <QWidget>

#include <atomic>
#include <memory>

class QLabel;
class QLineEdit;
class QPlainTextEdit;
class QPushButton;

namespace astap::stacking {
class LiveStackSession;
}

///----------------------------------------
namespace astap::gui {
///----------------------------------------

class ImageViewer;

///----------------------------------------
/// @class LiveStackWorker
/// @brief QThread that runs a LiveStackSession and emits Qt signals for
///        each session-level event (frame, message, finished).
///----------------------------------------

class LiveStackWorker final : public QThread {
	Q_OBJECT

public:
	explicit LiveStackWorker(QString watchDir, QObject* parent = nullptr);
	~LiveStackWorker() override;

	void requestStop();

signals:
	void message(QString text);
	void frameProcessed(int accepted, int rejected, int total);

protected:
	void run() override;

private:
	QString _watchDir;
	std::unique_ptr<astap::stacking::LiveStackSession> _session;
};

///----------------------------------------
/// @class LiveStackWindow
/// @brief Non-modal window with folder picker + Start/Pause/Stop + counters
///        + live log.
///----------------------------------------

class LiveStackWindow final : public QWidget {
	Q_OBJECT

public:
	explicit LiveStackWindow(QWidget* parent = nullptr);
	~LiveStackWindow() override;

	void setViewer(ImageViewer* viewer) { _viewer = viewer; }

signals:
	void stackUpdated(int counter);

protected:
	void closeEvent(QCloseEvent* event) override;

private slots:
	void browseFolder();
	void startStacking();
	void togglePause();
	void stopStacking();
	void onWorkerMessage(QString text);
	void onWorkerFrame(int accepted, int rejected, int total);
	void onWorkerFinished();

private:
	void setRunningUi(bool running);

	QLineEdit* _folderEdit = nullptr;
	QPushButton* _browseButton = nullptr;
	QPushButton* _startButton = nullptr;
	QPushButton* _pauseButton = nullptr;
	QPushButton* _stopButton = nullptr;
	QLabel* _countersLabel = nullptr;
	QPlainTextEdit* _log = nullptr;

	LiveStackWorker* _worker = nullptr;
	ImageViewer* _viewer = nullptr;
	bool _paused = false;
};

} // namespace astap::gui
