///----------------------------------------
///      @file main_window.h
///   @ingroup ASTAP++
///     @brief Top-level viewer window for the ASTAP++ Qt GUI.
///   @details Hosts the central @ref ImageViewer + @ref ControlsPanel
///            splitter, the menu bar, and the status bar (cursor pixel /
///            celestial position, plus image metadata). Drives plate-solve
///            jobs through @c QtConcurrent and refreshes WCS readouts on
///            completion.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <QFutureWatcher>
#include <QMainWindow>
#include <QPointF>
#include <memory>

class QLabel;
class QMenu;
class QCloseEvent;
class QDragEnterEvent;
class QDropEvent;

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @class MainWindow
/// @brief Main viewer window. Hosts the image canvas, menu, and status bar.
///----------------------------------------

class MainWindow final : public QMainWindow {
	Q_OBJECT

public:
	explicit MainWindow(QWidget* parent = nullptr);
	~MainWindow() override;

private slots:
	void openFile();
	void showAbout();
	void solveImage();
	void solveWithAstrometryNet();
	void onSolveFinished();
	void onCursorMoved(QPointF imagePos, bool inImage);
	void showSolverLog();
	void openRecent();
	void analyseStars();
	void onAnalyseFinished();
	void annotateDeepSky();
	void overlayCatalogStars();
	void overlayVariableStars();
	void overlaySimbadObjects();
	void overlayVizierGaia();
	void overlayAsteroids();
	void saveFile();
	void saveFileAs();
	void inspectImage();
	void openPhotometryDialog();
	void openAavsoDialog();
	void openSqmDialog();
	void openFocusDialog();
	void openPreferences();
	void startGradientRemoval();
	void startDustSpotRemoval();
	void onSelectionMade(QPointF start, QPointF end, int mode);

protected:
	void closeEvent(QCloseEvent* event) override;
	void dragEnterEvent(QDragEnterEvent* event) override;
	void dropEvent(QDropEvent* event) override;

private:
	void loadImageAt(const QString& path);
	void rebuildRecentMenu();
	void rememberRecent(const QString& path);

	void saveAppSettings() const;
	void restoreAppSettings();

	void updateCursorReadout(QPointF imagePos, bool inImage);
	void updateWcsReadout();
	void ensureLogWindow();

	std::unique_ptr<Ui::MainWindow> _ui;
	class AavsoDialog* _aavsoDialog = nullptr;
	class AstrometryNetDialog* _astrometryNetDialog = nullptr;
	class LogWindow* _logWindow = nullptr;
	class StackWindow* _stackWindow = nullptr;
	class ImageInspectorDialog* _inspectorDialog = nullptr;
	class LiveStackWindow* _liveStackWindow = nullptr;
	class PhotometryDialog* _photometryDialog = nullptr;
	class SqmDialog* _sqmDialog = nullptr;
	class FocusDialog* _focusDialog = nullptr;
	class QMenu* _recentMenu = nullptr;

	// Permanent status-bar widgets (cursor pixel, cursor celestial, image dims).
	QLabel* _cursorPixelLabel = nullptr;
	QLabel* _cursorCelestialLabel = nullptr;
	QLabel* _imageInfoLabel = nullptr;

	// Background solver job, if one is running.
	std::unique_ptr<QFutureWatcher<bool>> _solveWatcher;

	// Background star-detection job, if one is running.
	struct DetectionJob;
	std::unique_ptr<DetectionJob> _analyseJob;
};

} // namespace astap::gui
