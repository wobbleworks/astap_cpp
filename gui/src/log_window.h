///----------------------------------------
///      @file log_window.h
///   @ingroup ASTAP++
///     @brief Non-modal window displaying the latest solver log.
///   @details Shows the contents of @c astap::memo1_lines (populated by
///            @ref astap::solving::solve_image). Read-only; a Refresh
///            button re-pulls from the global so the user can reopen
///            after a new solve without reconstructing the window.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#pragma once

#include <QPointer>
#include <QWidget>

class QPlainTextEdit;

///----------------------------------------
namespace astap::gui {
///----------------------------------------

///----------------------------------------
/// @class LogWindow
/// @brief Read-only viewer for the latest plate-solver log.
///----------------------------------------

class LogWindow final : public QWidget {
	Q_OBJECT

public:
	explicit LogWindow(QWidget* parent = nullptr);
	~LogWindow() override = default;

	/// @brief Pull the current @c astap::memo1_lines into the text view.
	void refreshFromEngine();

private:
	QPlainTextEdit* _text = nullptr;
};

} // namespace astap::gui
