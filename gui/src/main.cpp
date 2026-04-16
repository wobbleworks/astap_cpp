///----------------------------------------
///      @file main.cpp
///   @ingroup ASTAP++
///     @brief Entry point for the ASTAP++ Qt 6 desktop application.
///    @author Created by John Stephen on 4/16/26.
/// @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
///----------------------------------------

#include "main_window.h"
#include "qt_image_decoder.h"

#include <QApplication>

int main(int argc, char* argv[]) {
	QApplication app(argc, argv);
	QCoreApplication::setOrganizationName("wobbleworks");
	QCoreApplication::setOrganizationDomain("wobbleworks.com");
	QCoreApplication::setApplicationName("ASTAP");

	// Install Qt-backed raster decoder so load_image can handle non-FITS
	// inputs (PNG / JPEG / BMP / TIFF, etc.). FITS is decoded directly by
	// the engine and doesn't need a decoder.
	(void)astap::gui::install_qt_image_decoder();

	astap::gui::MainWindow window;
	window.show();
	return QApplication::exec();
}
