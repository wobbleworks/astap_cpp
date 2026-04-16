##----------------------------------------------------------------------------
##      @file CMakeLists.txt
##   @ingroup ASTAP++
##     @brief Top-level build for the ASTAP++ CLI binary.
##   @details Compiles every .cpp under src/{core,image,reference,solving,stacking}
##            plus cli/{main,stb_image_decoder}.cpp into a single `astap`
##            executable. Tests are kept as their own CMake project at tests/.
##    @author Created by John Stephen on 4/15/26.
## @copyright Copyright © 2026 wobbleworks.com. All rights reserved.
##----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.20)
project(astap LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---- Static library ---------------------------------------------------------

file(GLOB_RECURSE LIB_SOURCES CONFIGURE_DEPENDS
	"${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_library(astap_lib STATIC ${LIB_SOURCES})
target_include_directories(astap_lib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_compile_features(astap_lib PUBLIC cxx_std_23)

if(CMAKE_SYSTEM_NAME MATCHES "iOS|watchOS|visionOS")
	target_compile_definitions(astap_lib PRIVATE ASTAP_NO_SYSTEM=1)
endif()

install(TARGETS astap_lib ARCHIVE DESTINATION lib)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/"
	DESTINATION include
	FILES_MATCHING PATTERN "*.h"
)

# ---- CLI executable (optional) ----------------------------------------------

option(ASTAP_BUILD_CLI "Build the astap command-line binary" ON)

if(ASTAP_BUILD_CLI)
	file(GLOB CLI_SOURCES CONFIGURE_DEPENDS
		"${CMAKE_CURRENT_SOURCE_DIR}/cli/*.cpp"
	)

	add_executable(astap ${CLI_SOURCES})
	target_link_libraries(astap PRIVATE astap_lib)
	target_include_directories(astap PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/cli")

	# Optional: sibling `stb` directory for stb_image.h
	if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../stb")
		target_include_directories(astap PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../stb")
	endif()
endif()

# ---- GUI executable (optional, Qt 6 desktop) --------------------------------

option(ASTAP_BUILD_GUI "Build the Qt 6 desktop GUI" OFF)

if(ASTAP_BUILD_GUI)
	add_subdirectory(gui)
endif()
