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

class QComboBox;
class QLabel;
class QLineEdit;
class QPlainTextEdit;
class QPushButton;

namespace astap::stacking {
class LiveStackSession;
class LiveMonitorSession;
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
	enum class Mode { LiveStack, MonitorOnly };

	LiveStackWorker(QString watchDir, Mode mode, QObject* parent = nullptr);
	~LiveStackWorker() override;

	void requestStop();

signals:
	/// @brief Emitted with each log line from the session.
	void message(QString text);

	/// @brief Emitted after a frame attempt.
	///        LiveStack mode: (accepted, rejected, total).
	///        MonitorOnly mode: (total, 0, total) — every frame counts.
	void frameProcessed(int accepted, int rejected, int total);

protected:
	void run() override;

private:
	QString _watchDir;
	Mode    _mode;
	std::unique_ptr<astap::stacking::LiveStackSession>   _stackSession;
	std::unique_ptr<astap::stacking::LiveMonitorSession> _monitorSession;
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
	void browseMasterDark();
	void browseMasterFlat();
	void clearMasterDark();
	void clearMasterFlat();
	void onWorkerMessage(QString text);
	void onWorkerFrame(int accepted, int rejected, int total);
	void onWorkerFinished();

private:
	void setRunningUi(bool running);
	void hydrateCalibrationFromSettings();

	QLineEdit* _folderEdit = nullptr;
	QPushButton* _browseButton = nullptr;
	QComboBox* _modeCombo = nullptr;
	QPushButton* _startButton = nullptr;
	QPushButton* _pauseButton = nullptr;
	QPushButton* _stopButton = nullptr;

	QLineEdit* _darkPath = nullptr;
	QLabel* _darkStatus = nullptr;
	QPushButton* _darkBrowse = nullptr;
	QPushButton* _darkClear = nullptr;
	QLineEdit* _flatPath = nullptr;
	QLabel* _flatStatus = nullptr;
	QPushButton* _flatBrowse = nullptr;
	QPushButton* _flatClear = nullptr;

	QLabel* _countersLabel = nullptr;
	QPlainTextEdit* _log = nullptr;

	LiveStackWorker* _worker = nullptr;
	ImageViewer* _viewer = nullptr;
	bool _paused = false;
};

} // namespace astap::gui
