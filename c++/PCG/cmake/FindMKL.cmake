#
# Find MKL
#
# Try to find MKL library.
# This module defines the following variables:
# - MKL_INCLUDE_DIRS
# - MKL_LIBRARIES
# - MKL_FOUND
#
# The following variables can be set as arguments for the module.
# - MKL_ROOT_DIR : Root library directory of MKL

# Additional modules
include(FindPackageHandleStandardArgs)

if (WIN32)
	# Find include files
	find_path(
		MKL_INCLUDE_DIR
		NAMES mkl.h
		PATHS
		$ENV{PROGRAMFILES}/include/
		${MKL_ROOT_DIR}/include/
		DOC "The directory where MKL.h resides")

	# Find library files
	find_library(
		MKL_LIBRARY
		NAMES mkl_core mkl_intel_lp64 mkl_sequential
		PATHS
		$ENV{PROGRAMFILES}/lib/
		${MKL_ROOT_DIR}/lib/intel64_win/
		${MKL_ROOT_DIR}/src/)
else()
	# Find include files
	find_path(
		MKL_INCLUDE_DIR
		NAMES mkl.h
		PATHS
		/usr/include/
		/usr/local/include/
		/sw/include/
		/opt/local/include/
		${MKL_ROOT_DIR}/include/
		DOC "The directory where MKL.h resides")

	# Find library files
	# Try to use static libraries
	find_library(
		MKL_LIBRARY
		NAMES mkl_core mkl_intel_lp64 mkl_sequential
		PATHS
		/usr/lib64/
		/usr/lib/
		/usr/local/lib64/
		/usr/local/lib/
		/sw/lib/
		/opt/local/lib/
		${MKL_ROOT_DIR}/lib/
		${MKL_ROOT_DIR}/lib/intel64/
		${MKL_ROOT_DIR}/src/
		DOC "The MKL library")
endif()

# Handle REQUIRD argument, define *_FOUND variable
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_INCLUDE_DIR MKL_LIBRARY)

# Define MKL_LIBRARIES and MKL_INCLUDE_DIRS
if (MKL_FOUND)
	set(MKL_LIBRARIES ${OPENGL_LIBRARIES} ${MKL_LIBRARY})
	set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})
endif()

# Hide some variables
mark_as_advanced(MKL_INCLUDE_DIR MKL_LIBRARY)