///----------------------------------------
///      @file main_window.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the top-level viewer window.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "main_window.h"
#include "controls_panel.h"
#include "image_viewer.h"
#include "log_window.h"
#include "solve_dialog.h"
#include "star_detector.h"
#include "ui_main_window.h"

#include "../../src/core/globals.h"
#include "../../src/core/image_io.h"
#include "../../src/core/wcs.h"
#include "../../src/reference/star_database.h"
#include "../../src/solving/astrometric_solving.h"

#include <QAction>
#include <QCloseEvent>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QFileDialog>
#include <QFileInfo>
#include <QFuture>
#include <QGuiApplication>
#include <QLabel>
#include <QMenu>
#include <QMessageBox>
#include <QMimeData>
#include <QSettings>
#include <QSignalBlocker>
#include <QStatusBar>
#include <QtConcurrent>

#include <cmath>
#include <numbers>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

namespace {

constexpr double kRad2Deg = 180.0 / std::numbers::pi;

// QSettings keys kept in one place so typos fail loudly.
constexpr auto kKeyGeometry       = "window/geometry";
constexpr auto kKeyWindowState    = "window/state";
constexpr auto kKeySplitter       = "window/splitter";
constexpr auto kKeyControlsShown  = "window/controlsShown";
constexpr auto kKeyRecentFiles    = "files/recent";
constexpr auto kKeyLastOpenDir    = "files/lastOpenDir";
constexpr auto kKeyDatabasePath   = "solver/databasePath";
constexpr auto kKeyDatabaseName   = "solver/databaseName";

constexpr int kMaxRecentFiles = 8;

} // namespace

struct MainWindow::DetectionJob {
	QFutureWatcher<DetectionResult> watcher;
};

MainWindow::MainWindow(QWidget* parent) :
	QMainWindow(parent),
	_ui(std::make_unique<Ui::MainWindow>()) {

	// Build the form
	_ui->setupUi(this);

	// Default splitter sizes: image takes most, controls panel ~280px.
	_ui->centralSplitter->setSizes({900, 280});

	// Bind controls panel to the viewer
	_ui->controlsPanel->attachViewer(_ui->imageViewer);

	// Permanent readouts in the status bar (right-justified).
	_imageInfoLabel = new QLabel(this);
	_cursorPixelLabel = new QLabel(this);
	_cursorCelestialLabel = new QLabel(this);
	_cursorCelestialLabel->setMinimumWidth(180);
	statusBar()->addPermanentWidget(_imageInfoLabel);
	statusBar()->addPermanentWidget(_cursorPixelLabel);
	statusBar()->addPermanentWidget(_cursorCelestialLabel);

	// Dynamic "Open Recent" submenu
	_recentMenu = new QMenu(tr("Open &Recent"), this);
	_ui->menuFile->insertMenu(_ui->actionQuit, _recentMenu);
	_ui->menuFile->insertSeparator(_ui->actionQuit);

	// Accept dropped files
	setAcceptDrops(true);

	// Menu wires
	connect(_ui->actionOpen, &QAction::triggered, this, &MainWindow::openFile);
	connect(_ui->actionQuit, &QAction::triggered, this, &QWidget::close);
	connect(_ui->actionAbout, &QAction::triggered, this, &MainWindow::showAbout);
	connect(_ui->actionSolve, &QAction::triggered, this, &MainWindow::solveImage);
	connect(_ui->actionShowSolverLog, &QAction::triggered, this, &MainWindow::showSolverLog);
	connect(_ui->actionAnalyse, &QAction::triggered, this, &MainWindow::analyseStars);
	connect(_ui->actionClearMarkers, &QAction::triggered,
		_ui->imageViewer, &ImageViewer::clearStarMarkers);

	connect(_ui->actionZoomIn, &QAction::triggered,
		_ui->imageViewer, &ImageViewer::zoomIn);
	connect(_ui->actionZoomOut, &QAction::triggered,
		_ui->imageViewer, &ImageViewer::zoomOut);
	connect(_ui->actionFitToWindow, &QAction::triggered,
		_ui->imageViewer, &ImageViewer::fitToWindow);
	connect(_ui->actionActualSize, &QAction::triggered,
		_ui->imageViewer, &ImageViewer::actualSize);
	connect(_ui->actionAutoStretch, &QAction::triggered,
		_ui->imageViewer, &ImageViewer::autoStretch);

	connect(_ui->actionFlipHorizontal, &QAction::toggled,
		_ui->imageViewer, &ImageViewer::setFlipHorizontal);
	connect(_ui->actionFlipVertical, &QAction::toggled,
		_ui->imageViewer, &ImageViewer::setFlipVertical);

	connect(_ui->imageViewer, &ImageViewer::flipChanged, this, [this]() {
		const QSignalBlocker bh(_ui->actionFlipHorizontal);
		const QSignalBlocker bv(_ui->actionFlipVertical);
		_ui->actionFlipHorizontal->setChecked(_ui->imageViewer->flipHorizontal());
		_ui->actionFlipVertical->setChecked(_ui->imageViewer->flipVertical());
	});

	connect(_ui->imageViewer, &ImageViewer::cursorMoved,
		this, &MainWindow::onCursorMoved);
	connect(_ui->imageViewer, &ImageViewer::imageLoaded,
		this, &MainWindow::updateWcsReadout);

	connect(_ui->actionToggleControls, &QAction::toggled,
		_ui->controlsPanel, &QWidget::setVisible);

	_ui->actionSolve->setEnabled(false);
	_ui->actionAnalyse->setEnabled(false);
	connect(_ui->imageViewer, &ImageViewer::imageLoaded, this, [this]() {
		const auto has = _ui->imageViewer->hasImage();
		_ui->actionSolve->setEnabled(has);
		_ui->actionAnalyse->setEnabled(has);
	});

	// Persisted state (geometry, splitter, recent files, database path...)
	restoreAppSettings();
	rebuildRecentMenu();

	statusBar()->showMessage(tr("Ready"));
}

MainWindow::~MainWindow() = default;

void MainWindow::closeEvent(QCloseEvent* event) {
	// Capture current layout before shutdown.
	saveAppSettings();
	QMainWindow::closeEvent(event);
}

void MainWindow::dragEnterEvent(QDragEnterEvent* event) {
	// Accept any local-file drop; dropEvent validates the extension.
	if (event->mimeData()->hasUrls()) {
		for (const auto& url : event->mimeData()->urls()) {
			if (url.isLocalFile()) {
				event->acceptProposedAction();
				return;
			}
		}
	}
}

void MainWindow::dropEvent(QDropEvent* event) {
	// Open the first dropped local file; ignore the rest.
	if (!event->mimeData()->hasUrls()) {
		return;
	}
	for (const auto& url : event->mimeData()->urls()) {
		if (url.isLocalFile()) {
			loadImageAt(url.toLocalFile());
			event->acceptProposedAction();
			return;
		}
	}
}

void MainWindow::openFile() {
	// Start the picker in the last-used directory.
	QSettings settings;
	const auto lastDir = settings.value(kKeyLastOpenDir).toString();

	const auto path = QFileDialog::getOpenFileName(
		this,
		tr("Open image"),
		lastDir,
		tr("FITS images (*.fit *.fits *.fts *.new *.fz);;"
		   "All images (*.fit *.fits *.fts *.new *.fz "
		                "*.png *.jpg *.jpeg *.bmp *.tif *.tiff "
		                "*.ppm *.pgm *.pfm *.xisf);;"
		   "All files (*)"));
	if (path.isEmpty()) {
		return;
	}
	loadImageAt(path);
}

void MainWindow::loadImageAt(const QString& path) {
	astap::filename2 = path.toStdString();

	QGuiApplication::setOverrideCursor(Qt::WaitCursor);
	const bool ok = astap::core::load_image(
		astap::filename2,
		astap::img_loaded,
		astap::head,
		astap::memo1_lines,
		/*re_center=*/false,
		/*plot=*/false);
	QGuiApplication::restoreOverrideCursor();

	if (!ok) {
		statusBar()->showMessage(tr("Failed to load: %1").arg(path));
		QMessageBox::warning(this, tr("Open"),
			tr("Failed to load %1.").arg(path));
		return;
	}

	// Hand the loaded buffer to the viewer.
	_ui->imageViewer->setImage(astap::img_loaded, astap::head);

	// Remember this location for the next open dialog + recent-files menu.
	QSettings settings;
	settings.setValue(kKeyLastOpenDir, QFileInfo(path).absolutePath());
	rememberRecent(path);

	setWindowTitle(tr("%1 — ASTAP").arg(QFileInfo(path).fileName()));
	_imageInfoLabel->setText(
		tr("%1 × %2  ·  %3 ch")
			.arg(astap::head.width)
			.arg(astap::head.height)
			.arg(astap::head.naxis3 > 0 ? astap::head.naxis3 : 1));
	statusBar()->showMessage(path);

	// Refresh the header/log window if it's already open.
	if (_logWindow && _logWindow->isVisible()) {
		_logWindow->refreshFromEngine();
	}
}

void MainWindow::openRecent() {
	auto* sender = qobject_cast<QAction*>(QObject::sender());
	if (!sender) {
		return;
	}
	loadImageAt(sender->data().toString());
}

void MainWindow::rebuildRecentMenu() {
	_recentMenu->clear();

	QSettings settings;
	const auto list = settings.value(kKeyRecentFiles).toStringList();
	if (list.isEmpty()) {
		auto* placeholder = _recentMenu->addAction(tr("(none)"));
		placeholder->setEnabled(false);
		return;
	}

	for (int i = 0; i < list.size(); ++i) {
		const auto& path = list[i];
		// Numbered shortcuts on the first nine entries
		const auto label = (i < 9)
			? QStringLiteral("&%1  %2").arg(i + 1).arg(QFileInfo(path).fileName())
			: QFileInfo(path).fileName();
		auto* action = _recentMenu->addAction(label);
		action->setData(path);
		action->setToolTip(path);
		connect(action, &QAction::triggered, this, &MainWindow::openRecent);
	}

	_recentMenu->addSeparator();
	auto* clear = _recentMenu->addAction(tr("Clear menu"));
	connect(clear, &QAction::triggered, this, [this]() {
		QSettings settings;
		settings.remove(kKeyRecentFiles);
		rebuildRecentMenu();
	});
}

void MainWindow::rememberRecent(const QString& path) {
	QSettings settings;
	auto list = settings.value(kKeyRecentFiles).toStringList();

	// Hoist any existing entry so the most-recent is always at index 0
	list.removeAll(path);
	list.prepend(path);
	while (list.size() > kMaxRecentFiles) {
		list.removeLast();
	}
	settings.setValue(kKeyRecentFiles, list);

	rebuildRecentMenu();
}

void MainWindow::solveImage() {
	if (!_ui->imageViewer->hasImage() || _solveWatcher) {
		return;
	}

	SolveDialog dialog(this);
	dialog.prefillFromHeader(astap::head);
	if (dialog.exec() != QDialog::Accepted) {
		return;
	}

	// Persist the chosen database location across sessions.
	QSettings settings;
	settings.setValue(kKeyDatabasePath,
		QString::fromStdString(astap::reference::database_path.string()));
	settings.setValue(kKeyDatabaseName,
		QString::fromStdString(astap::reference::name_database));

	_ui->actionSolve->setEnabled(false);
	statusBar()->showMessage(tr("Solving…"));
	QGuiApplication::setOverrideCursor(Qt::BusyCursor);

	_solveWatcher = std::make_unique<QFutureWatcher<bool>>();
	connect(_solveWatcher.get(), &QFutureWatcher<bool>::finished,
		this, &MainWindow::onSolveFinished);

	auto future = QtConcurrent::run([]() -> bool {
		return astap::solving::solve_image(
			astap::img_loaded,
			astap::head,
			astap::memo1_lines,
			/*get_hist=*/true,
			astap::check_pattern_filter);
	});
	_solveWatcher->setFuture(future);
}

void MainWindow::onSolveFinished() {
	const bool ok = _solveWatcher && _solveWatcher->result();

	QGuiApplication::restoreOverrideCursor();
	_ui->actionSolve->setEnabled(_ui->imageViewer->hasImage());

	_solveWatcher.reset();

	if (!ok) {
		statusBar()->showMessage(tr("Solve failed"));
		ensureLogWindow();
		_logWindow->refreshFromEngine();
		_logWindow->show();
		_logWindow->raise();
		QMessageBox::warning(this, tr("Plate solve"),
			tr("No solution found. The solver log has been opened."));
		return;
	}

	if (_logWindow && _logWindow->isVisible()) {
		_logWindow->refreshFromEngine();
	}

	updateWcsReadout();

	const auto cdelt = std::abs(astap::head.cdelt2) * 3600.0;
	statusBar()->showMessage(tr("Solved at %1 arcsec/px").arg(cdelt, 0, 'f', 2));
}

void MainWindow::analyseStars() {
	if (!_ui->imageViewer->hasImage() || _analyseJob) {
		return;
	}

	_ui->actionAnalyse->setEnabled(false);
	statusBar()->showMessage(tr("Detecting stars…"));
	QGuiApplication::setOverrideCursor(Qt::BusyCursor);

	// Copy the image data so the worker can operate independently of the
	// main thread. ImageArray is cheap to move but we need a strong copy
	// because the worker outlives this call.
	auto imgCopy = astap::img_loaded;

	_analyseJob = std::make_unique<DetectionJob>();
	connect(&_analyseJob->watcher, &QFutureWatcher<DetectionResult>::finished,
		this, &MainWindow::onAnalyseFinished);

	auto future = QtConcurrent::run(
		[img = std::move(imgCopy)]() -> DetectionResult {
			return detect_stars(img, /*snr_min=*/10.0);
		});
	_analyseJob->watcher.setFuture(future);
}

void MainWindow::onAnalyseFinished() {
	if (!_analyseJob) {
		return;
	}

	const auto result = _analyseJob->watcher.result();

	QGuiApplication::restoreOverrideCursor();
	_analyseJob.reset();
	_ui->actionAnalyse->setEnabled(_ui->imageViewer->hasImage());

	_ui->imageViewer->setStarMarkers(result.stars);

	if (result.stars.empty()) {
		statusBar()->showMessage(
			tr("No stars detected  ·  bgnd=%1  σ=%2  threshold=%3 (retry %4)  "
			   "candidates=%5 rejected=%6")
				.arg(result.background, 0, 'f', 0)
				.arg(result.noise, 0, 'f', 2)
				.arg(result.detectionLevel, 0, 'f', 0)
				.arg(result.retriesUsed)
				.arg(result.candidatesAbove)
				.arg(result.candidatesRejected));
	} else {
		statusBar()->showMessage(
			tr("%1 stars  ·  median HFD = %2  ·  bgnd = %3  ·  σ = %4  "
			   "(retry %5, threshold %6)")
				.arg(result.stars.size())
				.arg(result.medianHfd, 0, 'f', 2)
				.arg(result.background, 0, 'f', 0)
				.arg(result.noise, 0, 'f', 2)
				.arg(result.retriesUsed)
				.arg(result.detectionLevel, 0, 'f', 0));
	}
}

void MainWindow::onCursorMoved(QPointF imagePos, bool inImage) {
	updateCursorReadout(imagePos, inImage);
}

void MainWindow::updateCursorReadout(QPointF imagePos, bool inImage) {
	if (!inImage) {
		_cursorPixelLabel->clear();
		_cursorCelestialLabel->clear();
		return;
	}

	_cursorPixelLabel->setText(
		tr("(%1, %2)")
			.arg(imagePos.x(), 0, 'f', 1)
			.arg(imagePos.y(), 0, 'f', 1));

	if (astap::head.cdelt2 != 0.0) {
		double ra = 0.0;
		double dec = 0.0;
		astap::core::pixel_to_celestial(astap::head, imagePos.x(), imagePos.y(),
			/*formalism=*/0, ra, dec);
		const auto raStr = QString::fromStdString(astap::core::prepare_ra(ra, ":"));
		const auto decStr = QString::fromStdString(astap::core::prepare_dec(dec, ":"));
		_cursorCelestialLabel->setText(tr("%1  %2").arg(raStr, decStr));
	} else {
		_cursorCelestialLabel->clear();
	}
}

void MainWindow::updateWcsReadout() {
	if (astap::head.cdelt2 == 0.0) {
		_cursorCelestialLabel->clear();
	}
}

void MainWindow::showSolverLog() {
	ensureLogWindow();
	_logWindow->refreshFromEngine();
	_logWindow->show();
	_logWindow->raise();
	_logWindow->activateWindow();
}

void MainWindow::ensureLogWindow() {
	if (!_logWindow) {
		_logWindow = new LogWindow(this);
	}
}

void MainWindow::saveAppSettings() const {
	QSettings settings;
	settings.setValue(kKeyGeometry, saveGeometry());
	settings.setValue(kKeyWindowState, saveState());
	settings.setValue(kKeySplitter, _ui->centralSplitter->saveState());
	settings.setValue(kKeyControlsShown, _ui->controlsPanel->isVisible());
}

void MainWindow::restoreAppSettings() {
	QSettings settings;

	// Window geometry / dock layout / splitter sizes
	const auto geom = settings.value(kKeyGeometry).toByteArray();
	if (!geom.isEmpty()) {
		restoreGeometry(geom);
	}
	const auto state = settings.value(kKeyWindowState).toByteArray();
	if (!state.isEmpty()) {
		restoreState(state);
	}
	const auto splitter = settings.value(kKeySplitter).toByteArray();
	if (!splitter.isEmpty()) {
		_ui->centralSplitter->restoreState(splitter);
	}

	const auto controlsShown = settings.value(kKeyControlsShown, true).toBool();
	_ui->controlsPanel->setVisible(controlsShown);
	_ui->actionToggleControls->setChecked(controlsShown);

	// Seed the engine's solver globals from persisted settings so the
	// Plate Solve dialog shows the user's previous choices.
	const auto dbPath = settings.value(kKeyDatabasePath).toString();
	if (!dbPath.isEmpty()) {
		astap::reference::database_path = std::filesystem::path(dbPath.toStdString());
	}
	const auto dbName = settings.value(kKeyDatabaseName).toString();
	if (!dbName.isEmpty()) {
		astap::reference::name_database = dbName.toStdString();
	}
}

void MainWindow::showAbout() {
	QMessageBox::about(
		this,
		tr("About ASTAP"),
		tr("<h3>ASTAP</h3>"
		   "<p>Astrometric STAcking Program — C++23 / Qt 6 port.</p>"
		   "<p>Original ASTAP by Han Kleijn (www.hnsky.org).</p>"
		   "<p>C++ port by John Stephen / wobbleworks.com.</p>"));
}

} // namespace astap::gui
