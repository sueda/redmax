# Find pardiso
#
# Find the pardiso library
# 
# if you need to add a custom library search path, do it via CMAKE_PREFIX_PATH 
# 
# This module defines
#  PARDISO_LIBRARIES, the libraries needed to use PARDISO.
#  PARDISO_FOUND, If false, do not try to use PARDISO.
#  PARDISO_INCLUDE_PREFIX, include prefix for PARDISO

# Additional modules
include(FindPackageHandleStandardArgs)

IF(NOT PARDISO_ROOT_DIR)
 	SET(PARDISO_ROOT_DIR "$ENV{PARDISO_ROOT_DIR}")
ENDIF()

if (WIN32)
	# Find library files
	find_library(
		PARDISO_LIBRARY
		NAMES pardiso600-WIN-X86-64
		PATHS
		$ENV{PROGRAMFILES}/lib
		${PARDISO_ROOT_DIR}
		${PARDISO_ROOT_DIR}/lib
		${PARDISO_ROOT_DIR}/src
        ${PARDISO_ROOT_DIR}/src/Release
        ${PARDISO_ROOT_DIR}/src/Debug)
elseif(APPLE)	
	# Find library files
	# Try to use static libraries
	find_library(
		PARDISO_LIBRARY
		NAMES pardiso600-MACOS-X86-64
		PATHS
		/usr/lib64
		/usr/lib
		/usr/local/lib64
		/usr/local/lib
		/sw/lib
		/opt/local/lib
		${PARDISO_ROOT_DIR}
		${PARDISO_ROOT_DIR}/lib
		${PARDISO_ROOT_DIR}/src
		DOC "The PARDISO library")
else()
	# Find library files
	# Try to use static libraries
	find_library(
		PARDISO_LIBRARY
		NAMES pardiso600-MACOS-X86-64
		PATHS
		/usr/lib64
		/usr/lib
		/usr/local/lib64
		/usr/local/lib
		/sw/lib
		/opt/local/lib
		${PARDISO_ROOT_DIR}
		${PARDISO_ROOT_DIR}/lib
		${PARDISO_ROOT_DIR}/src
		DOC "The PARDISO library")
endif()

# Handle REQUIRD argument, define *_FOUND variable
find_package_handle_standard_args(PARDISO DEFAULT_MSG PARDISO_LIBRARY)
if (PARDISO_FOUND)
	set(PARDISO_LIBRARIES ${PARDISO_LIBRARY})	
endif()

# Hide some variables
mark_as_advanced (PARDISO_LIBRARY)