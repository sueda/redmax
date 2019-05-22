# Find jsoncpp
#
# Find the jsoncpp includes and library
# 
# if you need to add a custom library search path, do it via CMAKE_PREFIX_PATH 
# 
# This module defines
#  JSONCPP_INCLUDE_DIRS, where to find header, etc.
#  JSONCPP_LIBRARIES, the libraries needed to use jsoncpp.
#  JSONCPP_FOUND, If false, do not try to use jsoncpp.
#  JSONCPP_INCLUDE_PREFIX, include prefix for jsoncpp
# References:
# - https://github.com/majorx234/FindJsoncpp.cmake/blob/master/FindJsoncpp.cmake

# Additional modules
include(FindPackageHandleStandardArgs)

IF(NOT JSONCPP_ROOT_DIR)
 	SET(JSONCPP_ROOT_DIR "$ENV{JSONCPP_ROOT_DIR}")
ENDIF()

if (WIN32)
	# Find include files
	find_path(
		JSONCPP_INCLUDE_DIR
		NAMES json/json.h
		PATHS
		$ENV{PROGRAMFILES}/include/
		${JSONCPP_ROOT_DIR}/include/
		DOC "The directory where jsoncpp/json/json.h or json/json.h resides")

	# Find library files
	find_library(
		JSONCPP_LIBRARY
		NAMES jsoncpp
		PATHS
		$ENV{PROGRAMFILES}/lib/
		${JSONCPP_ROOT_DIR}/lib/
		${JSONCPP_ROOT_DIR}/src/
        ${JSONCPP_ROOT_DIR}/src/Release/
        ${JSONCPP_ROOT_DIR}/src/Debug/)
else()
	# Find include files
	find_path(
		JSONCPP_INCLUDE_DIR
		NAMES json/json.h
		PATHS
		${JSONCPP_ROOT_DIR}/include/
		DOC "The directory where jsoncpp/json/json.h or json/json.h resides")

	# Find library files
	# Try to use static libraries
	find_library(
		JSONCPP_LIBRARY
		NAMES jsoncpp libjsoncpp
		PATHS
		/usr/lib64/
		/usr/lib/
		/usr/local/lib64/
		/usr/local/lib/
		/sw/lib/
		/opt/local/lib/
		${JSONCPP_ROOT_DIR}/lib/
		${JSONCPP_ROOT_DIR}/src/
		DOC "The JSONCPP library")
endif()

# Handle REQUIRD argument, define *_FOUND variable
find_package_handle_standard_args(jsoncpp DEFAULT_MSG JSONCPP_INCLUDE_DIR JSONCPP_LIBRARY)
if (JSONCPP_FOUND)
	set(JSONCPP_INCLUDE_DIRS ${JSONCPP_INCLUDE_DIR})
	set(JSONCPP_LIBRARIES ${JSONCPP_LIBRARY})	
endif()

# Hide some variables
mark_as_advanced (JSONCPP_INCLUDE_DIR JSONCPP_LIBRARY)
