///----------------------------------------
///      @file stack_window.h
///   @ingroup ASTAP++
///     @brief Image stacking window.
///   @details Phase 5a MVP: Lights tab with a file list, add/remove/clear
///            buttons, and a Stack button that simple-averages the loaded
///            frames without calibration or alignment. Result is pushed
///            to the main viewer.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include "../../src/types.h"

#include <QWidget>

class QListWidget;
class QProgressBar;
class QPushButton;

///----------------------------------------
namespace astap::gui {
///----------------------------------------

class ImageViewer;

///----------------------------------------
/// @class StackWindow
/// @brief Non-modal stacking window with file list and stack controls.
///----------------------------------------

class StackWindow final : public QWidget {
	Q_OBJECT

public:
	explicit StackWindow(QWidget* parent = nullptr);
	~StackWindow() override = default;

	/// @brief Set the viewer that receives the stacked result.
	void setViewer(ImageViewer* viewer) { _viewer = viewer; }

signals:
	/// @brief Emitted after a successful stack so the main window can
	///        update title / status.
	void stackCompleted(int frameCount);

private slots:
	void addFiles();
	void removeSelected();
	void clearList();
	void startStack();

private:
	QListWidget* _fileList = nullptr;
	QPushButton* _addButton = nullptr;
	QPushButton* _removeButton = nullptr;
	QPushButton* _clearButton = nullptr;
	QPushButton* _stackButton = nullptr;
	QProgressBar* _progress = nullptr;

	ImageViewer* _viewer = nullptr;
};

} // namespace astap::gui
