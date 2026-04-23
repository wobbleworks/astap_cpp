///----------------------------------------
///      @file astrometry_net_dialog.cpp
///   @ingroup ASTAP++
///     @brief Implementation of @ref AstrometryNetDialog.
///    @author Created by John Stephen on 4/23/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "astrometry_net_dialog.h"

#include "../../src/solving/astrometry_net.h"

#include <QCheckBox>
#include <QDialogButtonBox>
#include <QFileDialog>
#include <QFormLayout>
#include <QFuture>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMetaObject>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSettings>
#include <QVBoxLayout>
#include <QtConcurrent>

#include <filesystem>
#include <string>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

namespace {
constexpr auto kKeyPath        = "astrometry_net/solve_field_path";
constexpr auto kKeyArgs        = "astrometry_net/extra_args";
constexpr auto kKeyCleanup     = "astrometry_net/cleanup_tmp";
constexpr auto kKeyShowConsole = "astrometry_net/show_console";
constexpr auto kKeyKeepConsole = "astrometry_net/keep_console_open";
}

AstrometryNetDialog::AstrometryNetDialog(QWidget* parent)
	: QDialog(parent) {
	setWindowTitle(tr("Solve with local astrometry.net"));
	setWindowFlag(Qt::WindowContextHelpButtonHint, false);
	resize(640, 480);

	auto* root = new QVBoxLayout(this);

	auto* form = new QFormLayout{};
	root->addLayout(form);

	// solve-field path + browse. On Windows the target is bash.exe
	// (Cygwin or WSL); on POSIX it's the solve-field binary or its dir.
	auto* pathRow = new QHBoxLayout{};
	_pathEdit = new QLineEdit(this);
#if defined(Q_OS_WIN)
	_pathEdit->setPlaceholderText(tr(
		"Path to Cygwin or WSL bash.exe (with astrometry.net installed)"));
#else
	_pathEdit->setPlaceholderText(tr(
		"Path to solve-field binary or directory containing it"));
#endif
	_browseButton = new QPushButton(tr("Browse…"), this);
	connect(_browseButton, &QPushButton::clicked,
	        this, &AstrometryNetDialog::browse);
	pathRow->addWidget(_pathEdit, /*stretch=*/1);
	pathRow->addWidget(_browseButton);
#if defined(Q_OS_WIN)
	form->addRow(tr("bash.exe:"), pathRow);
#else
	form->addRow(tr("solve-field:"), pathRow);
#endif

	// Extra arguments
	_extraArgsEdit = new QLineEdit(this);
	_extraArgsEdit->setPlaceholderText(tr("e.g. --downsample 2 --objs 150"));
	form->addRow(tr("Extra arguments:"), _extraArgsEdit);

	// Cleanup checkbox
	_cleanupCheck = new QCheckBox(
		tr("Remove scratch .wcs/.corr/.match/… files after a successful solve"),
		this);
	form->addRow(QString{}, _cleanupCheck);

#if defined(Q_OS_WIN)
	// Windows-only: console-visibility controls. On POSIX we stream
	// solve-field's output into the log pane below instead.
	_showConsoleCheck = new QCheckBox(tr("Show console window during solve"), this);
	form->addRow(QString{}, _showConsoleCheck);
	_keepConsoleCheck = new QCheckBox(
		tr("Keep console window open after solve finishes"), this);
	form->addRow(QString{}, _keepConsoleCheck);
#endif

	// Log area
	root->addWidget(new QLabel(tr("Output:"), this));
	_log = new QPlainTextEdit(this);
	_log->setReadOnly(true);
	_log->setLineWrapMode(QPlainTextEdit::NoWrap);
	{
		auto font = _log->font();
		font.setFamily(QStringLiteral("Menlo"));
		font.setStyleHint(QFont::Monospace);
		_log->setFont(font);
	}
	root->addWidget(_log, /*stretch=*/1);

	// Buttons
	auto* buttons = new QDialogButtonBox(this);
	_solveButton = buttons->addButton(tr("Solve"), QDialogButtonBox::AcceptRole);
	_closeButton = buttons->addButton(QDialogButtonBox::Close);
	connect(_solveButton, &QPushButton::clicked,
	        this, &AstrometryNetDialog::runSolve);
	connect(_closeButton, &QPushButton::clicked, this, &QDialog::close);
	root->addWidget(buttons);

	loadSettings();
}

AstrometryNetDialog::~AstrometryNetDialog() {
	if (_watcher) {
		_watcher->waitForFinished();
	}
}

void AstrometryNetDialog::setFitsPath(const QString& path) {
	_fitsPath = path;
	auto base = QFileInfo(path).fileName();
	if (base.isEmpty()) {
		setWindowTitle(tr("Solve with local astrometry.net"));
	} else {
		setWindowTitle(tr("Solve with local astrometry.net — %1").arg(base));
	}
}

void AstrometryNetDialog::browse() {
	auto hint = _pathEdit->text().trimmed();
	if (hint.isEmpty()) {
		hint = QStringLiteral("/usr/local/bin");
	}
	// Accept either a directory or an executable file — the engine handles
	// both forms.
	const auto picked = QFileDialog::getOpenFileName(this,
		tr("Locate solve-field"), hint);
	if (!picked.isEmpty()) {
		_pathEdit->setText(picked);
	}
}

void AstrometryNetDialog::runSolve() {
	if (_fitsPath.isEmpty()) {
		appendLog(tr("No image loaded."));
		return;
	}
	if (_watcher) {
		return;  // already running
	}

	saveSettings();
	_log->clear();
	_solveButton->setEnabled(false);

	astap::solving::AstrometryNetOptions opts;
	opts.fits_path        = std::filesystem::path(_fitsPath.toStdString());
	opts.solve_field_path = std::filesystem::path(_pathEdit->text().trimmed().toStdString());
	opts.extra_args       = _extraArgsEdit->text().toStdString();
	opts.cleanup_tmp      = _cleanupCheck->isChecked();
	if (_showConsoleCheck) opts.show_console      = _showConsoleCheck->isChecked();
	if (_keepConsoleCheck) opts.keep_console_open = _keepConsoleCheck->isChecked();

	// Marshal log lines from the worker thread back to the GUI thread via
	// a queued invokeMethod call on ourselves.
	auto logFn = [this](const std::string& line) {
		QMetaObject::invokeMethod(this, "appendLog", Qt::QueuedConnection,
			Q_ARG(QString, QString::fromStdString(line)));
	};

	_watcher = std::make_unique<QFutureWatcher<bool>>();
	connect(_watcher.get(), &QFutureWatcher<bool>::finished,
	        this, &AstrometryNetDialog::onFinished);
	auto future = QtConcurrent::run(
		[opts = std::move(opts), logFn = std::move(logFn)]() -> bool {
			return astap::solving::astrometry_net(opts, logFn);
		});
	_watcher->setFuture(future);
}

void AstrometryNetDialog::onFinished() {
	const auto ok = _watcher && _watcher->result();
	_watcher.reset();
	_solveButton->setEnabled(true);

	if (ok) {
		appendLog(tr("─── solve-field succeeded ───"));
		emit solved(_fitsPath);
	} else {
		appendLog(tr("─── solve-field failed ───"));
	}
}

void AstrometryNetDialog::appendLog(const QString& line) {
	_log->appendPlainText(line);
}

void AstrometryNetDialog::loadSettings() {
	QSettings settings;
	_pathEdit->setText(settings.value(kKeyPath).toString());
	_extraArgsEdit->setText(
		settings.value(kKeyArgs, QStringLiteral("--downsample 2")).toString());
	_cleanupCheck->setChecked(settings.value(kKeyCleanup, true).toBool());
	if (_showConsoleCheck) {
		_showConsoleCheck->setChecked(settings.value(kKeyShowConsole, true).toBool());
	}
	if (_keepConsoleCheck) {
		_keepConsoleCheck->setChecked(settings.value(kKeyKeepConsole, false).toBool());
	}
}

void AstrometryNetDialog::saveSettings() const {
	QSettings settings;
	settings.setValue(kKeyPath,    _pathEdit->text().trimmed());
	settings.setValue(kKeyArgs,    _extraArgsEdit->text());
	settings.setValue(kKeyCleanup, _cleanupCheck->isChecked());
	if (_showConsoleCheck) {
		settings.setValue(kKeyShowConsole, _showConsoleCheck->isChecked());
	}
	if (_keepConsoleCheck) {
		settings.setValue(kKeyKeepConsole, _keepConsoleCheck->isChecked());
	}
}

} // namespace astap::gui
