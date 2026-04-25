///----------------------------------------
///      @file main_window.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the top-level viewer window.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "main_window.h"
#include "aavso_dialog.h"
#include "annotation_scanner.h"
#include "astrometry_net_dialog.h"
#include "controls_panel.h"
#include "image_viewer.h"
#include "image_inspector_dialog.h"
#include "live_stack_window.h"
#include "log_window.h"
#include "focus_dialog.h"
#include "photometry_dialog.h"
#include "preferences_dialog.h"
#include "qt_http_client.h"
#include "solve_dialog.h"
#include "sqm_dialog.h"
#include "stack_window.h"
#include "star_detector.h"
#include "ui_main_window.h"

#include "../../src/core/fits.h"
#include "../../src/core/globals.h"
#include "../../src/core/image_io.h"
#include "../../src/analysis/asteroid_overlay.h"
#include "../../src/core/online.h"
#include "../../src/core/wcs.h"
#include "../../src/image/tiff.h"
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
	connect(_ui->actionSave, &QAction::triggered, this, &MainWindow::saveFile);
	connect(_ui->actionSaveAs, &QAction::triggered, this, &MainWindow::saveFileAs);
	connect(_ui->actionQuit, &QAction::triggered, this, &QWidget::close);
	connect(_ui->actionAbout, &QAction::triggered, this, &MainWindow::showAbout);
	connect(_ui->actionSolve, &QAction::triggered, this, &MainWindow::solveImage);
	connect(_ui->actionSolveAstrometryNet, &QAction::triggered,
		this, &MainWindow::solveWithAstrometryNet);
	connect(_ui->actionShowSolverLog, &QAction::triggered, this, &MainWindow::showSolverLog);
	connect(_ui->actionAnalyse, &QAction::triggered, this, &MainWindow::analyseStars);
	connect(_ui->actionInspect, &QAction::triggered, this, &MainWindow::inspectImage);
	connect(_ui->actionAnnotate, &QAction::triggered, this, &MainWindow::annotateDeepSky);
	connect(_ui->actionCatalogStars, &QAction::triggered, this, &MainWindow::overlayCatalogStars);
	connect(_ui->actionVariableStars, &QAction::triggered, this, &MainWindow::overlayVariableStars);
	connect(_ui->actionSimbadObjects, &QAction::triggered, this, &MainWindow::overlaySimbadObjects);
	connect(_ui->actionVizierGaia, &QAction::triggered, this, &MainWindow::overlayVizierGaia);
	connect(_ui->actionAsteroids, &QAction::triggered, this, &MainWindow::overlayAsteroids);
	connect(_ui->actionPhotometry, &QAction::triggered, this, &MainWindow::openPhotometryDialog);
	connect(_ui->actionAavsoReport, &QAction::triggered, this, &MainWindow::openAavsoDialog);
	connect(_ui->actionSqm, &QAction::triggered, this, &MainWindow::openSqmDialog);
	connect(_ui->actionFocus, &QAction::triggered, this, &MainWindow::openFocusDialog);
	connect(_ui->actionPreferences, &QAction::triggered, this, &MainWindow::openPreferences);
	connect(_ui->actionStack, &QAction::triggered, this, [this]() {
		if (!_stackWindow) {
			_stackWindow = new StackWindow(this);
			_stackWindow->setViewer(_ui->imageViewer);
			connect(_stackWindow, &StackWindow::stackCompleted,
				this, [this](int n) {
					setWindowTitle(tr("Stacked %1 frames — ASTAP").arg(n));
					_imageInfoLabel->setText(
						tr("%1 × %2  ·  %3 ch  ·  %4 frames")
							.arg(astap::head.width)
							.arg(astap::head.height)
							.arg(astap::head.naxis3 > 0 ? astap::head.naxis3 : 1)
							.arg(n));
					statusBar()->showMessage(
						tr("Stacked %1 frames (simple average)").arg(n));
				});
		}
		_stackWindow->show();
		_stackWindow->raise();
		_stackWindow->activateWindow();
	});
	connect(_ui->actionLiveStack, &QAction::triggered, this, [this]() {
		if (!_liveStackWindow) {
			_liveStackWindow = new LiveStackWindow(this);
			_liveStackWindow->setViewer(_ui->imageViewer);
			connect(_liveStackWindow, &LiveStackWindow::stackUpdated,
				this, [this](int n) {
					setWindowTitle(tr("Live — %1 frames — ASTAP").arg(n));
					statusBar()->showMessage(
						tr("Live stack: %1 frames").arg(n));
				});
		}
		_liveStackWindow->show();
		_liveStackWindow->raise();
		_liveStackWindow->activateWindow();
	});
	connect(_ui->actionClearMarkers, &QAction::triggered, this, [this]() {
		_ui->imageViewer->clearStarMarkers();
		_ui->imageViewer->clearAnnotations();
		_ui->imageViewer->clearConstellations();
		_ui->imageViewer->clearCatalogStars();
		_ui->imageViewer->clearVarStars();
		_ui->imageViewer->clearSimbadObjects();
		_ui->imageViewer->clearVizierStars();
		_ui->imageViewer->clearAsteroids();
	});

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
	_ui->actionAnnotate->setEnabled(false);
	_ui->actionSave->setEnabled(false);
	_ui->actionSaveAs->setEnabled(false);
	connect(_ui->imageViewer, &ImageViewer::imageLoaded, this, [this]() {
		const auto has = _ui->imageViewer->hasImage();
		_ui->actionSolve->setEnabled(has);
		_ui->actionAnalyse->setEnabled(has);
		_ui->actionAnnotate->setEnabled(has && astap::head.cdelt2 != 0.0);
		_ui->actionSave->setEnabled(has);
		_ui->actionSaveAs->setEnabled(has);
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

void MainWindow::solveWithAstrometryNet() {
	if (!_ui->imageViewer->hasImage()) {
		QMessageBox::information(this, tr("Solve with astrometry.net"),
			tr("Open a FITS file first."));
		return;
	}
	if (astap::filename2.empty()) {
		QMessageBox::information(this, tr("Solve with astrometry.net"),
			tr("The current image has no filename on disk "
			   "(it may have come from a drag-dropped buffer). "
			   "Save it first, then re-open it from disk."));
		return;
	}

	if (!_astrometryNetDialog) {
		_astrometryNetDialog = new AstrometryNetDialog(this);
		connect(_astrometryNetDialog, &AstrometryNetDialog::solved,
			this, [this](const QString& path) {
				// Reload from disk — the FITS now carries solve-field's WCS.
				loadImageAt(path);
				statusBar()->showMessage(tr("Solved with astrometry.net"));
			});
	}
	_astrometryNetDialog->setFitsPath(
		QString::fromStdString(astap::filename2));
	_astrometryNetDialog->show();
	_astrometryNetDialog->raise();
	_astrometryNetDialog->activateWindow();
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
	_ui->actionAnnotate->setEnabled(true);

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

// Produce a clean FITS header from memo1_lines. The solver pushes log
// strings after END, and its update_float/update_text may insert WCS
// cards after END too (because log lines moved the "last line" past END).
// This helper keeps every line that looks like a valid FITS card and
// discards solver log debris, then ensures END is the very last entry.
static void sanitize_memo_for_save(std::vector<std::string>& memo) {
	auto isFitsCard = [](const std::string& line) -> bool {
		if (line.empty()) {
			return false;
		}
		// Keyword = value cards have '=' at column 9.
		if (line.size() >= 10 && line[8] == '=') {
			return true;
		}
		// COMMENT, HISTORY, and blank cards.
		if (line.starts_with("COMMENT") || line.starts_with("HISTORY")) {
			return true;
		}
		// All-space padding.
		if (line.find_first_not_of(' ') == std::string::npos) {
			return false;
		}
		return false;
	};

	// Collect valid FITS cards from the entire memo (before AND after END),
	// skipping END itself and any non-FITS log lines.
	auto clean = std::vector<std::string>{};
	clean.reserve(memo.size());
	for (auto& line : memo) {
		if (line.starts_with("END")) {
			continue;
		}
		if (isFitsCard(line)) {
			clean.push_back(std::move(line));
		}
	}
	clean.emplace_back(
		"END                                                                             ");
	memo = std::move(clean);
}

void MainWindow::saveFile() {
	if (astap::filename2.empty()) {
		saveFileAs();
		return;
	}

	const auto path = std::filesystem::path(astap::filename2);
	const auto ext = path.extension().string();
	const auto isFits = (ext == ".fit" || ext == ".fits" || ext == ".fts"
		|| ext == ".FIT" || ext == ".FITS" || ext == ".FTS");

	if (isFits) {
		auto memo = astap::memo1_lines;
		sanitize_memo_for_save(memo);
		const auto ok = astap::core::savefits_update_header(memo, path);
		if (ok) {
			statusBar()->showMessage(tr("Header saved to %1").arg(
				QString::fromStdString(astap::filename2)));
		} else {
			QMessageBox::warning(this, tr("Save"),
				tr("Failed to update the header of %1.").arg(
					QString::fromStdString(astap::filename2)));
		}
	} else {
		saveFileAs();
	}
}

namespace {

/// @brief Clamp a float sample to [0, 65535] and round to 16-bit integer.
[[nodiscard]] std::uint16_t clamp_u16(float v) noexcept {
	if (!(v > 0.0f)) {
		return 0;
	}
	if (v >= 65535.0f) {
		return 65535;
	}
	return static_cast<std::uint16_t>(v + 0.5f);
}

/// @brief Build a 16-bit QImage from the raw (un-stretched) image data.
/// @details Uses @c QImage::Format_Grayscale16 for single-channel sources and
///          @c QImage::Format_RGBX64 for three-channel. Bayer frames are
///          exported via the grayscale path. Pixel values are clamped to
///          [0, 65535] to fit the 16-bit container — this preserves full
///          precision for 16-bit inputs and loses nothing of consequence for
///          stacked float results (the science value is linear even after
///          clamp, since very bright stars are already saturated anyway).
[[nodiscard]] QImage build_raw_16bit_image(const astap::ImageArray& img) {
	if (img.empty() || img[0].empty() || img[0][0].empty()) {
		return {};
	}
	const auto height = static_cast<int>(img[0].size());
	const auto width  = static_cast<int>(img[0][0].size());
	const auto channels = static_cast<int>(img.size());

	if (channels >= 3) {
		// 16-bit-per-channel RGBA; alpha set to 65535 (opaque).
		QImage out(width, height, QImage::Format_RGBX64);
		for (int y = 0; y < height; ++y) {
			auto* row = reinterpret_cast<quint16*>(out.scanLine(y));
			for (int x = 0; x < width; ++x) {
				row[x * 4 + 0] = clamp_u16(img[0][y][x]);
				row[x * 4 + 1] = clamp_u16(img[1][y][x]);
				row[x * 4 + 2] = clamp_u16(img[2][y][x]);
				row[x * 4 + 3] = 65535;
			}
		}
		return out;
	}

	QImage out(width, height, QImage::Format_Grayscale16);
	for (int y = 0; y < height; ++y) {
		auto* row = reinterpret_cast<quint16*>(out.scanLine(y));
		for (int x = 0; x < width; ++x) {
			row[x] = clamp_u16(img[0][y][x]);
		}
	}
	return out;
}

}  // namespace

void MainWindow::saveFileAs() {
	// Multi-format save. Filter order matches the extension-dispatch order
	// below so the dialog's "save as type" preselects the right handler.
	auto selectedFilter = QString{};
	const auto path = QFileDialog::getSaveFileName(
		this,
		tr("Save As"),
		QString(),
		tr("FITS images (*.fits *.fit *.fts);;"
		   "TIFF images (*.tif *.tiff);;"
		   "PNG images (*.png);;"
		   "JPEG images (*.jpg *.jpeg);;"
		   "All files (*)"),
		&selectedFilter);
	if (path.isEmpty()) {
		return;
	}

	// Resolve format from the chosen extension. Append one from the filter if
	// the user did not type one.
	auto suffix = QFileInfo(path).suffix().toLower();
	auto outPath = path;
	if (suffix.isEmpty()) {
		if      (selectedFilter.contains("TIFF")) { suffix = "tif";  outPath += ".tif";  }
		else if (selectedFilter.contains("PNG"))  { suffix = "png";  outPath += ".png";  }
		else if (selectedFilter.contains("JPEG")) { suffix = "jpg";  outPath += ".jpg";  }
		else                                      { suffix = "fits"; outPath += ".fits"; }
	}

	auto ok = false;

	if (suffix == "fits" || suffix == "fit" || suffix == "fts") {
		auto memo = astap::memo1_lines;
		sanitize_memo_for_save(memo);
		ok = astap::core::save_fits(
			astap::img_loaded,
			memo,
			std::filesystem::path(outPath.toStdString()),
			/*type1=*/16,
			/*override2=*/true);
	}
	else if (suffix == "tif" || suffix == "tiff") {
		// 16-bit grayscale or 48-bit RGB depending on channel count.
		const auto channels = astap::img_loaded.size();
		const auto stdPath = std::filesystem::path(outPath.toStdString());
		const auto desc = "Saved by ASTAP";
		if (channels >= 3) {
			ok = astap::image::save_tiff_48(astap::img_loaded, stdPath,
				desc, /*flip_h=*/false, /*flip_v=*/false);
		} else {
			ok = astap::image::save_tiff_16(astap::img_loaded, stdPath,
				desc, /*flip_h=*/false, /*flip_v=*/false);
		}
	}
	else if (suffix == "png") {
		// Raw 16-bit export — preserves scientific value. Qt writes
		// 16-bit PNG when the source QImage is Format_Grayscale16 /
		// Format_RGBX64.
		const auto img = build_raw_16bit_image(astap::img_loaded);
		if (!img.isNull()) {
			ok = img.save(outPath, "PNG");
		}
	}
	else if (suffix == "jpg" || suffix == "jpeg") {
		// JPEG is 8-bit only, so use the viewer's stretched render —
		// that's what the user actually wants to share.
		const auto& rendered = _ui->imageViewer->renderedImage();
		if (!rendered.isNull()) {
			ok = rendered.save(outPath, "JPEG", /*quality=*/92);
		}
	}
	else {
		QMessageBox::warning(this, tr("Save As"),
			tr("Unsupported extension: .%1").arg(suffix));
		return;
	}

	if (ok) {
		astap::filename2 = outPath.toStdString();
		statusBar()->showMessage(tr("Saved to %1").arg(outPath));
	} else {
		QMessageBox::warning(this, tr("Save As"),
			tr("Failed to save %1.").arg(outPath));
	}
}

void MainWindow::inspectImage() {
	if (!_ui->imageViewer->hasImage()) {
		QMessageBox::information(this, tr("Inspect Image"),
			tr("Open an image first."));
		return;
	}
	if (!_inspectorDialog) {
		_inspectorDialog = new ImageInspectorDialog(this);
	}
	_inspectorDialog->show();
	_inspectorDialog->raise();
	_inspectorDialog->activateWindow();
	_inspectorDialog->analyseCurrentImage();
}

void MainWindow::openPhotometryDialog() {
	if (!_ui->imageViewer->hasImage()) {
		QMessageBox::information(this, tr("Photometric Calibration"),
			tr("Open an image first."));
		return;
	}
	if (astap::head.cd1_1 == 0.0) {
		QMessageBox::information(this, tr("Photometric Calibration"),
			tr("Plate-solve the image first (Image → Plate Solve)."));
		return;
	}
	if (!_photometryDialog) {
		_photometryDialog = new PhotometryDialog(this);
	}
	_photometryDialog->prefillFromSettings();
	_photometryDialog->show();
	_photometryDialog->raise();
	_photometryDialog->activateWindow();
}

void MainWindow::openAavsoDialog() {
	if (!_ui->imageViewer->hasImage()) {
		QMessageBox::information(this, tr("AAVSO Report"),
			tr("Open an image first."));
		return;
	}
	if (astap::head.cd1_1 == 0.0) {
		QMessageBox::information(this, tr("AAVSO Report"),
			tr("Plate-solve the image first (Image → Plate Solve)."));
		return;
	}
	if (astap::head.mzero == 0.0) {
		const auto answer = QMessageBox::question(this, tr("AAVSO Report"),
			tr("Photometric calibration has not been run on this image, "
			   "so MZERO is unset and magnitudes can't be computed.\n\n"
			   "Open the AAVSO dialog anyway?"),
			QMessageBox::Yes | QMessageBox::No);
		if (answer != QMessageBox::Yes) return;
	}
	if (!_aavsoDialog) {
		_aavsoDialog = new AavsoDialog(this, _ui->imageViewer);
	}
	_aavsoDialog->show();
	_aavsoDialog->raise();
	_aavsoDialog->activateWindow();
}

void MainWindow::openSqmDialog() {
	if (!_ui->imageViewer->hasImage()) {
		QMessageBox::information(this, tr("Sky Quality Meter"),
			tr("Open an image first."));
		return;
	}
	if (astap::head.cd1_1 == 0.0) {
		QMessageBox::information(this, tr("Sky Quality Meter"),
			tr("Plate-solve the image first (Image → Plate Solve)."));
		return;
	}
	if (!_sqmDialog) {
		_sqmDialog = new SqmDialog(this);
	}
	_sqmDialog->prefillFromSettings();
	_sqmDialog->show();
	_sqmDialog->raise();
	_sqmDialog->activateWindow();
}

void MainWindow::openFocusDialog() {
	// Unlike SQM/photometry, Focus doesn't need a loaded image — the user
	// picks the sweep's frames inside the dialog.
	if (!_focusDialog) {
		_focusDialog = new FocusDialog(this);
	}
	_focusDialog->show();
	_focusDialog->raise();
	_focusDialog->activateWindow();
}

void MainWindow::openPreferences() {
	PreferencesDialog dlg(this);
	if (dlg.exec() == QDialog::Accepted) {
		// Refresh the recent-files menu in case it was cleared.
		rebuildRecentMenu();
	}
}

void MainWindow::annotateDeepSky() {
	if (astap::head.cdelt2 == 0.0) {
		QMessageBox::information(this, tr("Annotate"),
			tr("Plate-solve the image first (Image → Plate Solve)."));
		return;
	}

	// Ensure the catalog is loaded. Try the database path first, then
	// the directory containing the loaded image (some installs put the
	// CSV alongside the database files).
	if (!load_deepsky_catalog(astap::reference::database_path)) {
		const auto imgDir = std::filesystem::path(astap::filename2).parent_path();
		if (!load_deepsky_catalog(imgDir)) {
			QMessageBox::warning(this, tr("Annotate"),
				tr("Could not load deep_sky.csv from the database path (%1).\n\n"
				   "Make sure the ASTAP catalog files are installed.")
					.arg(QString::fromStdString(
						astap::reference::database_path.string())));
			return;
		}
	}

	QGuiApplication::setOverrideCursor(Qt::WaitCursor);
	const auto markers = annotate_image(astap::head);
	const auto consOverlay = build_constellation_overlay(astap::head);
	QGuiApplication::restoreOverrideCursor();

	_ui->imageViewer->setAnnotations(markers);
	_ui->imageViewer->setConstellations(consOverlay);
	statusBar()->showMessage(
		tr("%1 deep-sky objects, %2 constellation lines")
			.arg(markers.size())
			.arg(consOverlay.lines.size()));
}

void MainWindow::overlayCatalogStars() {
	if (astap::head.cdelt2 == 0.0) {
		QMessageBox::information(this, tr("Catalog Star Overlay"),
			tr("Plate-solve the image first (Image → Plate Solve)."));
		return;
	}

	// Toggle behaviour: if an overlay is already showing, clear it.
	if (!_ui->imageViewer->catalogStars().empty()) {
		_ui->imageViewer->clearCatalogStars();
		statusBar()->showMessage(tr("Catalog star overlay cleared"));
		return;
	}

	// The plate-solve flow leaves `astap::reference::database_type` and the
	// selected catalog name populated; scan_catalog_stars uses them.
	if (astap::reference::name_database.empty()) {
		if (!astap::reference::select_star_database("",
				astap::head.height * std::abs(astap::head.cdelt2))) {
			QMessageBox::warning(this, tr("Catalog Star Overlay"),
				tr("No star catalog is available. Install the ASTAP catalog "
				   "(e.g. d05 or d50) in the database path."));
			return;
		}
	}

	QGuiApplication::setOverrideCursor(Qt::WaitCursor);
	auto markers = scan_catalog_stars(astap::head, /*max_markers=*/1000);
	QGuiApplication::restoreOverrideCursor();

	if (markers.empty()) {
		QMessageBox::information(this, tr("Catalog Star Overlay"),
			tr("No catalog stars found in this field. Verify the plate solve "
			   "and that the matching catalog covers this sky area."));
		return;
	}

	const auto n = markers.size();
	_ui->imageViewer->setCatalogStars(std::move(markers));
	statusBar()->showMessage(tr("%1 catalog stars overlaid").arg(n));
}

void MainWindow::overlayVariableStars() {
	if (astap::head.cdelt2 == 0.0) {
		QMessageBox::information(this, tr("Variable Stars"),
			tr("Plate-solve the image first (Image → Plate Solve)."));
		return;
	}
	if (!_ui->imageViewer->varStars().empty()) {
		_ui->imageViewer->clearVarStars();
		statusBar()->showMessage(tr("Variable-star overlay cleared"));
		return;
	}

	QGuiApplication::setOverrideCursor(Qt::WaitCursor);
	auto http = QtHttpClient{};
	auto loggedMsg = QString{};
	auto log = [&](const std::string& s) {
		if (loggedMsg.isEmpty()) loggedMsg = QString::fromStdString(s);
	};
	// Mode 6 → online VSX/VSP, limiting magnitude 15 (a sensible default for
	// most amateur fields). 4 = m11, 5 = m13, 6 = m15, 7 = unlimited.
	astap::core::variable_star_annotation(http, astap::head, /*mode=*/6,
	                                      /*years_since_2000=*/0.0,
	                                      /*extract_visible=*/false, log);
	auto markers = project_var_stars(astap::head);
	QGuiApplication::restoreOverrideCursor();

	if (markers.empty()) {
		QMessageBox::information(this, tr("Variable Stars"),
			loggedMsg.isEmpty()
				? tr("No VSX/VSP entries found in this field.")
				: loggedMsg);
		return;
	}

	const auto n = markers.size();
	_ui->imageViewer->setVarStars(std::move(markers));
	statusBar()->showMessage(tr("%1 VSX/VSP markers overlaid").arg(n));
}

void MainWindow::overlaySimbadObjects() {
	if (astap::head.cdelt2 == 0.0) {
		QMessageBox::information(this, tr("Simbad Objects"),
			tr("Plate-solve the image first (Image → Plate Solve)."));
		return;
	}
	if (!_ui->imageViewer->simbadObjects().empty()) {
		_ui->imageViewer->clearSimbadObjects();
		statusBar()->showMessage(tr("Simbad overlay cleared"));
		return;
	}

	const auto url = astap::core::make_simbad_url(astap::head,
		astap::core::SimbadQuery::DeepSky);

	QGuiApplication::setOverrideCursor(Qt::WaitCursor);
	auto http = QtHttpClient{};
	auto reply = http.get(url);
	QGuiApplication::restoreOverrideCursor();

	if (!reply) {
		QMessageBox::warning(this, tr("Simbad Objects"),
			tr("Simbad request failed: %1").arg(QString::fromStdString(reply.error())));
		return;
	}
	const auto objects = astap::core::plot_simbad(*reply);
	auto markers = project_simbad(objects, astap::head);
	if (markers.empty()) {
		QMessageBox::information(this, tr("Simbad Objects"),
			tr("Simbad returned no deep-sky objects in this field."));
		return;
	}
	const auto n = markers.size();
	_ui->imageViewer->setSimbadObjects(std::move(markers));
	statusBar()->showMessage(tr("%1 Simbad objects overlaid").arg(n));
}

void MainWindow::overlayVizierGaia() {
	if (astap::head.cdelt2 == 0.0) {
		QMessageBox::information(this, tr("Vizier Gaia"),
			tr("Plate-solve the image first (Image → Plate Solve)."));
		return;
	}
	if (!_ui->imageViewer->vizierStars().empty()) {
		_ui->imageViewer->clearVizierStars();
		statusBar()->showMessage(tr("Vizier overlay cleared"));
		return;
	}

	// Default limiting Gaia G magnitude — matches the VSX default tier.
	constexpr auto kLimitMag = 15.0;
	const auto url = astap::core::make_vizier_gaia_url(astap::head, kLimitMag);

	QGuiApplication::setOverrideCursor(Qt::WaitCursor);
	auto http = QtHttpClient{};
	auto reply = http.get(url);
	QGuiApplication::restoreOverrideCursor();

	if (!reply) {
		QMessageBox::warning(this, tr("Vizier Gaia"),
			tr("Vizier request failed: %1").arg(QString::fromStdString(reply.error())));
		return;
	}
	// No filter transform — keep raw G magnitudes for now (filter selection
	// can be exposed in a follow-up dialog along with the magnitude limit).
	const auto rows = astap::core::plot_vizier(*reply, "G", nullptr);
	auto markers = project_vizier(rows, astap::head);
	if (markers.empty()) {
		QMessageBox::information(this, tr("Vizier Gaia"),
			tr("Vizier returned no Gaia rows in this field."));
		return;
	}
	const auto n = markers.size();
	_ui->imageViewer->setVizierStars(std::move(markers));
	statusBar()->showMessage(tr("%1 Vizier Gaia stars overlaid").arg(n));
}

void MainWindow::overlayAsteroids() {
	if (astap::head.cdelt2 == 0.0) {
		QMessageBox::information(this, tr("Asteroids & Comets"),
			tr("Plate-solve the image first (Image → Plate Solve)."));
		return;
	}
	if (!_ui->imageViewer->asteroids().empty()) {
		_ui->imageViewer->clearAsteroids();
		statusBar()->showMessage(tr("Asteroid overlay cleared"));
		return;
	}

	QSettings settings;
	auto mpcorbPath  = settings.value("asteroid/mpcorb_path").toString();
	auto cometsPath  = settings.value("asteroid/comets_path").toString();

	// Prompt if we don't have an MPCORB path yet.
	if (mpcorbPath.isEmpty() || !QFileInfo::exists(mpcorbPath)) {
		const auto picked = QFileDialog::getOpenFileName(this,
			tr("Locate MPCORB.DAT"),
			QString{},
			tr("Asteroid catalogs (*.DAT *.dat);;All files (*)"));
		if (picked.isEmpty()) return;
		mpcorbPath = picked;
		settings.setValue("asteroid/mpcorb_path", mpcorbPath);
	}

	astap::analysis::AsteroidScanOptions opts;
	opts.mpcorb_path = std::filesystem::path(mpcorbPath.toStdString());
	if (!cometsPath.isEmpty() && QFileInfo::exists(cometsPath)) {
		opts.comets_path = std::filesystem::path(cometsPath.toStdString());
	}
	opts.max_magnitude = settings.value("asteroid/max_magnitude", 18.0).toDouble();
	opts.max_count     = settings.value("asteroid/max_count", 1500000).toInt();

	QGuiApplication::setOverrideCursor(Qt::WaitCursor);
	const auto detections = astap::analysis::scan_asteroids_in_field(astap::head, opts);
	auto markers = project_asteroids(detections, astap::head);
	QGuiApplication::restoreOverrideCursor();

	if (markers.empty()) {
		QMessageBox::information(this, tr("Asteroids & Comets"),
			tr("No asteroids or comets brighter than %1 mag fall in this field.\n\n"
			   "(If you expected hits, verify the MPCORB.DAT path under "
			   "Settings → Preferences, or increase the magnitude limit in "
			   "asteroid/max_magnitude.)")
				.arg(opts.max_magnitude, 0, 'f', 1));
		return;
	}

	const auto n = markers.size();
	_ui->imageViewer->setAsteroids(std::move(markers));
	statusBar()->showMessage(tr("%1 asteroid/comet markers overlaid").arg(n));
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
	const auto version = QString::fromStdString(astap::astap_version);
	const auto qtVer   = QString::fromLatin1(qVersion());
	const auto dbPath  = QString::fromStdString(
		astap::reference::database_path.string());
	const auto dbName  = QString::fromStdString(
		astap::reference::name_database);

	const auto html = tr(
		"<h3>ASTAP</h3>"
		"<p><b>Astrometric STAcking Program</b> — C++23 / Qt 6 desktop port.</p>"
		"<table cellpadding='3'>"
		"<tr><td><b>Engine</b></td><td>%1</td></tr>"
		"<tr><td><b>Qt</b></td><td>%2</td></tr>"
		"<tr><td><b>Catalog</b></td><td>%3%4</td></tr>"
		"</table>"
		"<p>Original ASTAP by Han Kleijn "
		"(<a href=\"https://www.hnsky.org\">hnsky.org</a>).<br>"
		"C++ port by John Stephen / wobbleworks.com.</p>"
		"<p>Licensed under the Mozilla Public License 2.0.</p>")
		.arg(version.isEmpty() ? tr("unknown") : version,
		     qtVer,
		     dbName.isEmpty() ? tr("<i>none selected</i>") : dbName,
		     dbPath.isEmpty() ? QString() : tr(" &nbsp;(%1)").arg(dbPath));

	QMessageBox box(this);
	box.setWindowTitle(tr("About ASTAP"));
	box.setTextFormat(Qt::RichText);
	box.setTextInteractionFlags(Qt::TextBrowserInteraction);
	box.setText(html);
	box.setIconPixmap({});   // No default icon — the HTML is self-contained.
	box.addButton(QMessageBox::Close);
	box.exec();
}

} // namespace astap::gui
