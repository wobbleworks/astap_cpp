##----------------------------------------------------------------------------
##      @file CMakeLists.txt
##   @ingroup ASTAP++
##     @brief Build for the ASTAP++ Qt 6 desktop GUI.
##   @details Builds a Qt 6 / C++23 desktop application linking against
##            astap_lib. Qt is located via (1) sibling third_party/qt/
##            populated by build_qt, (2) sibling ../qt/ when ASTAP++ is
##            cloned standalone, or (3) system Qt via standard find_package.
##    @author Created by John Stephen on 4/16/26.
## @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
##----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.20)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# ---- Locate Qt --------------------------------------------------------------

if(NOT DEFINED Qt6_DIR AND NOT DEFINED CMAKE_PREFIX_PATH)
	if(APPLE)
		set(_qt_host "mac")
	elseif(WIN32)
		set(_qt_host "windows")
	else()
		set(_qt_host "linux")
	endif()

	# Search third_party/qt (when ASTAP++ is a submodule of wobbleworks)
	# and ../qt (when ASTAP++ is cloned standalone with a sibling qt dir).
	# aqtinstall lays out as <root>/<host>/<version>/<arch>/, where <arch>
	# is e.g. clang_64, macos, linux_gcc_64, win64_msvc2022_64. We probe
	# every candidate and pick whichever exposes lib/cmake/Qt6/.
	foreach(_root
			"${CMAKE_CURRENT_SOURCE_DIR}/../../qt"
			"${CMAKE_CURRENT_SOURCE_DIR}/../qt")
		if(IS_DIRECTORY "${_root}/${_qt_host}")
			file(GLOB _qt_arches LIST_DIRECTORIES TRUE
				"${_root}/${_qt_host}/*/*")
			foreach(_a ${_qt_arches})
				if(IS_DIRECTORY "${_a}/lib/cmake/Qt6")
					list(APPEND CMAKE_PREFIX_PATH "${_a}")
				endif()
			endforeach()
		endif()
	endforeach()
endif()

find_package(Qt6 6.4 REQUIRED COMPONENTS Core Gui Widgets Concurrent Network)

qt_standard_project_setup()

# ---- Sources ----------------------------------------------------------------

set(GUI_SOURCES
	src/main.cpp
	src/main_window.cpp
	src/main_window.h
	src/main_window.ui
	src/image_viewer.cpp
	src/image_viewer.h
	src/histogram_widget.cpp
	src/histogram_widget.h
	src/controls_panel.cpp
	src/controls_panel.h
	src/solve_dialog.cpp
	src/solve_dialog.h
	src/log_window.cpp
	src/log_window.h
	src/star_detector.cpp
	src/star_detector.h
	src/image_inspector_dialog.cpp
	src/image_inspector_dialog.h
	src/live_stack_window.cpp
	src/live_stack_window.h
	src/annotation_scanner.cpp
	src/annotation_scanner.h
	src/stack_window.cpp
	src/stack_window.h
	src/photometry_dialog.cpp
	src/photometry_dialog.h
	src/sqm_dialog.cpp
	src/sqm_dialog.h
	src/focus_dialog.cpp
	src/focus_dialog.h
	src/preferences_dialog.cpp
	src/preferences_dialog.h
	src/qt_image_decoder.cpp
	src/qt_image_decoder.h
	src/qt_http_client.cpp
	src/qt_http_client.h
)

# ---- Target -----------------------------------------------------------------

qt_add_executable(astap_gui ${GUI_SOURCES})

target_include_directories(astap_gui PRIVATE
	"${CMAKE_CURRENT_SOURCE_DIR}/src"
)

target_link_libraries(astap_gui PRIVATE
	astap_lib
	Qt6::Core
	Qt6::Gui
	Qt6::Widgets
	Qt6::Concurrent
	Qt6::Network
)

target_compile_features(astap_gui PRIVATE cxx_std_23)

set_target_properties(astap_gui PROPERTIES
	OUTPUT_NAME "ASTAP"
	MACOSX_BUNDLE TRUE
	WIN32_EXECUTABLE TRUE
)
