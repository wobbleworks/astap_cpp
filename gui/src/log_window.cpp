///----------------------------------------
///      @file log_window.cpp
///   @ingroup ASTAP++
///     @brief Implementation of the solver log viewer.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "log_window.h"

#include "../../src/core/globals.h"

#include <QFontDatabase>
#include <QHBoxLayout>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QVBoxLayout>

///----------------------------------------
namespace astap::gui {
///----------------------------------------

LogWindow::LogWindow(QWidget* parent) :
	QWidget(parent, Qt::Window) {

	setWindowTitle(tr("FITS header / solver log"));
	resize(720, 420);

	// Text area: read-only monospace
	_text = new QPlainTextEdit(this);
	_text->setReadOnly(true);
	_text->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
	_text->setLineWrapMode(QPlainTextEdit::NoWrap);

	// Refresh / Close buttons
	auto* refresh = new QPushButton(tr("Refresh"), this);
	auto* close = new QPushButton(tr("Close"), this);
	connect(refresh, &QPushButton::clicked, this, &LogWindow::refreshFromEngine);
	connect(close, &QPushButton::clicked, this, &QWidget::close);

	auto* buttons = new QHBoxLayout();
	buttons->addWidget(refresh);
	buttons->addStretch(1);
	buttons->addWidget(close);

	auto* root = new QVBoxLayout(this);
	root->addWidget(_text, 1);
	root->addLayout(buttons);

	refreshFromEngine();
}

void LogWindow::refreshFromEngine() {
	// Concatenate the engine's log vector into the text widget.
	QString joined;
	joined.reserve(static_cast<qsizetype>(astap::memo1_lines.size()) * 64);
	for (const auto& line : astap::memo1_lines) {
		joined += QString::fromStdString(line);
		joined += QLatin1Char('\n');
	}
	_text->setPlainText(joined);

	// Scroll to the bottom so the most recent messages are visible.
	auto cursor = _text->textCursor();
	cursor.movePosition(QTextCursor::End);
	_text->setTextCursor(cursor);
}

} // namespace astap::gui
