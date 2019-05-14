#
# Find GLFW
#
# Try to find GLFW library.
# This module defines the following variables:
# - GLFW_INCLUDE_DIRS
# - GLFW_LIBRARIES
# - GLFW_FOUND
#
# The following variables can be set as arguments for the module.
# - GLFW_ROOT_DIR : Root library directory of GLFW
# - GLFW_USE_STATIC_LIBS : Specifies to use static version of GLFW library (Windows only)
#
# References:
# - https://github.com/progschj/OpenGL-Examples/blob/master/cmake_modules/FindGLFW.cmake
# - https://bitbucket.org/Ident8/cegui-mk2-opengl3-renderer/src/befd47200265/cmake/FindGLFW.cmake
# - https://github.com/lighttransport/nanogi/blob/master/cmake/FindGLFW.cmake
# - https://github.com/LevSmirnov1997/vse/blob/d11b7d30f4451e2c990d76fe88ef72599963fe94/al/cmake/FindGLFW.cmake

# Additional modules
include(FindPackageHandleStandardArgs)

if (WIN32)
	# Find include files
	find_path(
		GLFW_INCLUDE_DIR
		NAMES GLFW/glfw3.h
		PATHS
		$ENV{PROGRAMFILES}/include
		${GLFW_ROOT_DIR}/include
		DOC "The directory where GLFW/glfw.h resides")

	# Find library files
	find_library(
		GLFW_LIBRARY
		NAMES glfw3
		PATHS
		$ENV{PROGRAMFILES}/lib
		${GLFW_ROOT_DIR}/lib
		${GLFW_ROOT_DIR}/src
        ${GLFW_ROOT_DIR}/src/Release
        ${GLFW_ROOT_DIR}/src/Debug)

else()
	# Find include files
	find_path(
		GLFW_INCLUDE_DIR
		NAMES GLFW/glfw3.h
		PATHS
		/usr/include
		/usr/local/include
		/sw/include
		/opt/local/include
		${GLFW_ROOT_DIR}/include
		DOC "The directory where GL/glfw.h resides")

	# Find library files
	# Try to use static libraries
	find_library(
		GLFW_LIBRARY
		NAMES glfw3
		PATHS
		/usr/lib64
		/usr/lib
		/usr/local/lib64
		/usr/local/lib
		/sw/lib
		/opt/local/lib
		${GLFW_ROOT_DIR}/lib
		${GLFW_ROOT_DIR}/src
		DOC "The GLFW library")
endif()

# Handle REQUIRD argument, define *_FOUND variable
find_package_handle_standard_args(GLFW DEFAULT_MSG GLFW_INCLUDE_DIR GLFW_LIBRARY)

# Define GLFW_LIBRARIES and GLFW_INCLUDE_DIRS
if (GLFW_FOUND)
	set(GLFW_LIBRARIES ${OPENGL_LIBRARIES} ${GLFW_LIBRARY})
	set(GLFW_INCLUDE_DIRS ${GLFW_INCLUDE_DIR})
endif()

# Hide some variables
mark_as_advanced(GLFW_INCLUDE_DIR GLFW_LIBRARY)