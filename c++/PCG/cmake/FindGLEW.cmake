#
# Find GLEW
#
# Try to find GLEW library.
# This module defines the following variables:
# - GLEW_INCLUDE_DIRS
# - GLEW_LIBRARIES
# - GLEW_FOUND
#
# The following variables can be set as arguments for the module.
# - GLEW_ROOT_DIR : Root library directory of GLEW

# Additional modules
include(FindPackageHandleStandardArgs)

if (WIN32)
	# Find include files
	find_path(
		GLEW_INCLUDE_DIR
		NAMES GL/glew.h
		PATHS
		$ENV{PROGRAMFILES}/include
		${GLEW_ROOT_DIR}/include
		DOC "The directory where GLEW/GLEW.h resides")

	# Find library files
	find_library(
		GLEW_LIBRARY
		NAMES glew32s
		PATHS
		$ENV{PROGRAMFILES}/lib
		${GLEW_ROOT_DIR}/lib
		${GLEW_ROOT_DIR}/lib/Release/x64/
		${GLEW_ROOT_DIR}/lib/Debug/x64/
		${GLEW_ROOT_DIR}/src
        ${GLEW_ROOT_DIR}/src/Release
        ${GLEW_ROOT_DIR}/src/Debug)
else()
	# Find include files
	find_path(
		GLEW_INCLUDE_DIR
		NAMES GL/glew.h
		PATHS
		/usr/include
		/usr/local/include
		/sw/include
		/opt/local/include
		${GLEW_ROOT_DIR}/include
		DOC "The directory where GL/GLEW.h resides")

	# Find library files
	# Try to use static libraries
	find_library(
		GLEW_LIBRARY
		NAMES glew
		PATHS
		/usr/lib64
		/usr/lib
		/usr/local/lib64
		/usr/local/lib
		/sw/lib
		/opt/local/lib
		${GLEW_ROOT_DIR}/lib
		${GLEW_ROOT_DIR}/src
		DOC "The GLEW library")
endif()

# Handle REQUIRD argument, define *_FOUND variable
find_package_handle_standard_args(GLEW DEFAULT_MSG GLEW_INCLUDE_DIR GLEW_LIBRARY)

# Define GLEW_LIBRARIES and GLEW_INCLUDE_DIRS
if (GLEW_FOUND)
	set(GLEW_LIBRARIES ${OPENGL_LIBRARIES} ${GLEW_LIBRARY})
	set(GLEW_INCLUDE_DIRS ${GLEW_INCLUDE_DIR})
endif()

# Hide some variables
mark_as_advanced(GLEW_INCLUDE_DIR GLEW_LIBRARY)