///----------------------------------------
///      @file focus_dialog.h
///   @ingroup ASTAP++
///     @brief Focuser V-curve dialog: load a sweep of FITS files, fit a
///            hyperbola to their HFDs, and report the best-focus position.
///    @author Created by John Stephen on 4/21/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <QDialog>

class QLabel;
class QListWidget;
class QPushButton;
class QTextEdit;

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @class FocusDialog
/// @brief Modeless dialog for focuser V-curve analysis.
/// @details User drops in a set of FITS frames taken across a focuser sweep;
///          the dialog reads FOCUSPOS from each, measures median HFD via
///          analyse_image, and fits a hyperbola to the curve. Reports the
///          best-focus position plus a per-frame (position, HFD, stars)
///          listing.
///----------------------------------------

class FocusDialog final : public QDialog {
	Q_OBJECT

public:
	explicit FocusDialog(QWidget* parent = nullptr);
	~FocusDialog() override = default;

private slots:
	void addFiles();
	void removeSelected();
	void clearList();
	void runFit();

private:
	void buildLayout();

	QListWidget* _fileList = nullptr;
	QPushButton* _addButton = nullptr;
	QPushButton* _removeButton = nullptr;
	QPushButton* _clearButton = nullptr;
	QPushButton* _runButton = nullptr;

	QLabel*    _bestFocusLabel = nullptr;
	QLabel*    _residualLabel = nullptr;
	QLabel*    _samplesLabel = nullptr;
	QTextEdit* _log = nullptr;
};

} // namespace astap::gui
