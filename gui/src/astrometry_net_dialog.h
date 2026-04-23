///----------------------------------------
///      @file astrometry_net_dialog.h
///   @ingroup ASTAP++
///     @brief Modeless dialog that drives a local astrometry.net
///            @c solve-field run as a fallback to the built-in solver.
///    @author Created by John Stephen on 4/23/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <QDialog>
#include <QFutureWatcher>
#include <QString>

#include <memory>

class QCheckBox;
class QLineEdit;
class QPlainTextEdit;
class QPushButton;

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @class AstrometryNetDialog
/// @brief Modeless dialog for running the local astrometry.net
///        @c solve-field on the currently-loaded FITS.
/// @details Emits @ref solved when the run succeeds so the host can
///          reload the FITS (which now carries the new WCS headers).
///----------------------------------------

class AstrometryNetDialog final : public QDialog {
	Q_OBJECT

public:
	explicit AstrometryNetDialog(QWidget* parent = nullptr);
	~AstrometryNetDialog() override;

	/// @brief Set the FITS path to solve. Prefill the title label.
	void setFitsPath(const QString& path);

signals:
	/// @brief Emitted after a successful solve. @p path is the FITS file
	///        whose header was just updated (the caller should reload it).
	void solved(const QString& path);

private slots:
	void browse();
	void runSolve();
	void onFinished();
	void appendLog(const QString& line);

private:
	void loadSettings();
	void saveSettings() const;

	QLineEdit* _pathEdit = nullptr;
	QPushButton* _browseButton = nullptr;
	QLineEdit* _extraArgsEdit = nullptr;
	QCheckBox* _cleanupCheck = nullptr;
	QCheckBox* _showConsoleCheck = nullptr;   // Windows only (null elsewhere)
	QCheckBox* _keepConsoleCheck = nullptr;   // Windows only (null elsewhere)
	QPlainTextEdit* _log = nullptr;
	QPushButton* _solveButton = nullptr;
	QPushButton* _closeButton = nullptr;

	QString _fitsPath;
	std::unique_ptr<QFutureWatcher<bool>> _watcher;
};

} // namespace astap::gui
