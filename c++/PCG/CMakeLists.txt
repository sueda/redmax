CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
SET(CMAKE_CXX_STANDARD 14)

# Name of the project
PROJECT(RedMaxPCG)

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# find_package(dlib REQUIRED)
# find_package(jsoncpp CONFIG REQUIRED)
# find_package(intel-mkl REQUIRED)

# FOR LAB MACHINE
# SET(DEF_DIR_GLM "D:\\lib_c++\\glm-0.9.8.5\\glm")
# SET(DEF_DIR_GLFW "D:\\lib_c++\\glfw-3.2.1\\glfw-3.2.1")
# SET(DEF_DIR_GLEW "D:\\lib_c++\\glew-2.1.0-win32\\glew-2.1.0")
# SET(DEF_DIR_EIGEN "D:\\lib_c++\\eigen-eigen-5a0156e40feb\\eigen-eigen-5a0156e40feb")
# SET(DEF_DIR_MOSEK "D:\\lib_c++\\mosek\\8\\tools\\platform\\win64x86")
# SET(DEF_DIR_LIBIGL_INCLUDE "D:\\lib_c++\\libigl-master\\libigl-master\\include")
# SET(DEF_DIR_PARDISO "D:\\lib_c++\\PARDISO")
# SET(DEF_DIR_JSONCPP "D:\\lib_c++\\vcpkg\\packages\\jsoncpp_x64-windows")

# FOR MAC
SET(DEF_DIR_GLM "/Users/yingwang/Libraries/glm-0.9.9.3")
SET(DEF_DIR_GLFW "/Users/yingwang/Libraries/glfw-master")
SET(DEF_DIR_GLEW "/Users/yingwang/Libraries/glew-2.1.0")
SET(DEF_DIR_LIBIGL_INCLUDE "/Users/yingwang/Libraries/libigl-2.0.0/include")
SET(DEF_DIR_PARDISO "/Users/yingwang/Libraries/pardiso")
SET(DEF_DIR_MKL "/opt/intel/mkl")
SET(DEF_DIR_JSONCPP "/Users/yingwang/Libraries/jsoncpp-0.10.7")

# Is this the solution?
# Override with `cmake -DSOL=ON ..`
OPTION(SOL "Solution" OFF)

# Use glob to get the list of all source files.
# We don't really need to include header and resource files to build, but it's
# nice to have them also show up in IDEs.
IF(${SOL})
  FILE(GLOB_RECURSE SOURCES "src0/*.cpp")
  FILE(GLOB_RECURSE HEADERS "src0/*.h")
ELSE()
  FILE(GLOB_RECURSE SOURCES "src/*.cpp" "src/*.cc")
  FILE(GLOB_RECURSE HEADERS "src/*.h" "src/*.hh" "src/*.hpp")
ENDIF()
FILE(GLOB_RECURSE GLSL "resources/*.glsl")

#INCLUDE_DIRECTORIES(src/qp_interface)

# Set the executable.
ADD_EXECUTABLE(${CMAKE_PROJECT_NAME} ${SOURCES} ${HEADERS} ${GLSL})

# Get the EIGEN environment variable. (header-only)
find_package (Eigen3 3.3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

# Get the GLM environment variable. (header-only)
SET(GLM_INCLUDE_DIR "$ENV{GLM_INCLUDE_DIR}")
IF(NOT GLM_INCLUDE_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable GLM_INCLUDE_DIR to the root directory of your GLM installation.")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for GLM in ${DEF_DIR_GLM}")
    IF(IS_DIRECTORY ${DEF_DIR_GLM})
      MESSAGE(STATUS "Found!")
      SET(GLM_INCLUDE_DIR ${DEF_DIR_GLM})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSE()
  	MESSAGE(STATUS "Looking for GLM in ${DEF_DIR_GLM}")
  	IF(IS_DIRECTORY ${DEF_DIR_GLM})
      MESSAGE(STATUS "Found!")
      SET(GLM_INCLUDE_DIR ${DEF_DIR_GLM})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ENDIF()
ENDIF()
INCLUDE_DIRECTORIES(${GLM_INCLUDE_DIR})

# Get the GLFW environment variable. There should be a CMakeLists.txt in the 
# specified directory.
SET(GLFW_DIR "$ENV{GLFW_DIR}")
IF(NOT GLFW_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable GLFW_DIR to the root directory of your GLFW installation.")
  # IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for GLFW in ${DEF_DIR_GLFW}")
    IF(IS_DIRECTORY ${DEF_DIR_GLFW})
      MESSAGE(STATUS "Found!")
      SET(GLFW_DIR ${DEF_DIR_GLFW})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  # ELSE()
  #   MESSAGE(FATAL_ERROR ${ERR_MSG})
  # ENDIF()
ENDIF()
OPTION(GLFW_BUILD_EXAMPLES "GLFW_BUILD_EXAMPLES" OFF)
OPTION(GLFW_BUILD_TESTS "GLFW_BUILD_TESTS" OFF)
OPTION(GLFW_BUILD_DOCS "GLFW_BUILD_DOCS" OFF)
IF(CMAKE_BUILD_TYPE MATCHES Release)
  ADD_SUBDIRECTORY(${GLFW_DIR} ${GLFW_DIR}/release)
ELSE()
  ADD_SUBDIRECTORY(${GLFW_DIR} ${GLFW_DIR}/debug)
ENDIF()
INCLUDE_DIRECTORIES(${GLFW_DIR}/include)
TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} glfw ${GLFW_LIBRARIES})

# Get the GLEW environment variable.
SET(GLEW_DIR "$ENV{GLEW_DIR}")
IF(NOT GLEW_DIR)
  # The environment variable was not set
  SET(ERR_MSG "Please point the environment variable GLEW_DIR to the root directory of your GLEW installation.")
  # IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for GLEW in ${DEF_DIR_GLEW}")
    IF(IS_DIRECTORY ${DEF_DIR_GLEW})
      MESSAGE(STATUS "Found!")
      SET(GLEW_DIR ${DEF_DIR_GLEW})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  # ELSE()
  #   MESSAGE(FATAL_ERROR ${ERR_MSG})
  # ENDIF()
ENDIF()
INCLUDE_DIRECTORIES(${GLEW_DIR}/include)
IF(WIN32)
  # With prebuilt binaries
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${GLEW_DIR}/lib/Release/x64/glew32s.lib)
ELSE()
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${GLEW_DIR}/lib/libGLEW.a)
ENDIF()

# # For MOSEK functionality
# SET(MOSEK_DIR "$ENV{MOSEK_DIR}")
# IF(NOT MOSEK_DIR)
  # # The environment variable was not set
  # SET(ERR_MSG "Please point the environment variable MOSEK_DIR to the root directory of your MOSEK installation.")
  # IF(WIN32)
    # # On Windows, try the default location
    # MESSAGE(STATUS "Looking for Mosek in ${MOSEK_DIR}")
    # IF(IS_DIRECTORY ${DEF_DIR_MOSEK})
      # MESSAGE(STATUS "Found!")
      # SET(MOSEK_DIR ${DEF_DIR_MOSEK})
    # ELSE()
      # MESSAGE(FATAL_ERROR ${ERR_MSG})
    # ENDIF()
  # ELSE()
    # MESSAGE(FATAL_ERROR ${ERR_MSG})
  # ENDIF()
# ENDIF()
# INCLUDE_DIRECTORIES(${MOSEK_DIR}/h)
# TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${MOSEK_DIR}/bin/mosek64_8_1.lib)

# # TBB stuff
# SET(TBB_INCLUDE_DIR "$ENV{TBB_INCLUDE_DIR}")
# IF(NOT TBB_INCLUDE_DIR)
  # # The environment variable was not set
  # SET(ERR_MSG "Please point the environment variable TBB_INCLUDE_DIR to the root directory of your TBB installation.")
  # IF(WIN32)
    # # On Windows, try the default location
    # MESSAGE(STATUS "Looking for TBB in ${TBB_INCLUDE_DIR}")
    # IF(IS_DIRECTORY ${DEF_DIR_TBB})
      # MESSAGE(STATUS "Found!")
      # SET(TBB_INCLUDE_DIR ${DEF_DIR_TBB})
    # ELSE()
      # MESSAGE(FATAL_ERROR ${ERR_MSG})
    # ENDIF()
  # ELSE()
    # MESSAGE(FATAL_ERROR ${ERR_MSG})
  # ENDIF()
# ENDIF()
# TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${TBB_INCLUDE_DIR}/lib/intel64/vc14/tbb.lib)
# TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${TBB_INCLUDE_DIR}/lib/intel64/vc14/tbb_debug.lib)
# INCLUDE_DIRECTORIES(${TBB_INCLUDE_DIR}/include/tbb)

# IGL stuff
#SET(LIBIGL_INCLUDE_DIR "$ENV{LIBIGL_INCLUDE_DIR}")
#IF(NOT LIBIGL_INCLUDE_DIR)
#  # The environment variable was not set
#  SET(ERR_MSG "Please point the environment variable LIBIGL_INCLUDE_DIR to the root directory of your IGL installation.")
#  IF(WIN32)
#    # On Windows, try the default location
#    MESSAGE(STATUS "Looking for IGL in ${LIBIGL_INCLUDE_DIR}")
#    IF(IS_DIRECTORY ${DEF_DIR_LIBIGL_INCLUDE})
#      MESSAGE(STATUS "Found!")
#      SET(LIBIGL_INCLUDE_DIR ${DEF_DIR_LIBIGL_INCLUDE})
#    ELSE()
#      MESSAGE(FATAL_ERROR ${ERR_MSG})
#    ENDIF()
#  ELSE()
#    MESSAGE(FATAL_ERROR ${ERR_MSG})
#  ENDIF()
#ENDIF()
#INCLUDE_DIRECTORIES(${LIBIGL_INCLUDE_DIR})

# #Include nlopt-master
# SET(NLOPT_INCLUDE_DIR "$ENV{NLOPT_INCLUDE_DIR}")
# IF(NOT NLOPT_INCLUDE_DIR)
 # # The environment variable was not set
 # SET(ERR_MSG "Please point the environment variable NLOPT_INCLUDE_DIR to the root directory of your nlopt installation.")
  # IF(WIN32)
    # # On Windows, try the default location
    # MESSAGE(STATUS "Looking for nlopt in ${DEF_DIR_NLOPT}")
    # IF(IS_DIRECTORY ${DEF_DIR_NLOPT})
      # MESSAGE(STATUS "Found!")
      # SET(NLOPT_INCLUDE_DIR ${DEF_DIR_NLOPT})
    # ELSE()
      # MESSAGE(FATAL_ERROR ${ERR_MSG})
    # ENDIF()
  # ELSE()
    # MESSAGE(FATAL_ERROR ${ERR_MSG})
  # ENDIF()
# ENDIF()
# INCLUDE_DIRECTORIES(${NLOPT_INCLUDE_DIR})
# TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${NLOPT_INCLUDE_DIR}/libnlopt-0.lib)

#Include PARDISO
SET(PARDISO_INCLUDE_DIR "$ENV{PARDISO_INCLUDE_DIR}")
IF(NOT PARDISO_INCLUDE_DIR)
 # The environment variable was not set
 SET(ERR_MSG "Please point the environment variable PARDISO_INCLUDE_DIR to the root directory of your PARDISO installation.")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for PARDISO in ${DEF_DIR_PARDISO}")
    IF(IS_DIRECTORY ${DEF_DIR_PARDISO})
      MESSAGE(STATUS "Found!")
      SET(PARDISO_INCLUDE_DIR ${DEF_DIR_PARDISO})
      TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${PARDISO_INCLUDE_DIR}/libpardiso600-WIN-X86-64.lib)
    ENDIF()
  ELSEIF(APPLE)
  	MESSAGE(STATUS "Looking for PARDISO in ${DEF_DIR_PARDISO}")
    IF(IS_DIRECTORY ${DEF_DIR_PARDISO})
      MESSAGE(STATUS "Found!")
      SET(PARDISO_INCLUDE_DIR ${DEF_DIR_PARDISO})
      TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} -L${PARDISO_INCLUDE_DIR}/ libpardiso600-MACOS-X86-64.dylib)
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR ${ERR_MSG})
  ENDIF()
ENDIF()

#Include JSONCPP
SET(JSONCPP_INCLUDE_DIR "$ENV{JSONCPP_INCLUDE_DIR}")
IF(NOT JSONCPP_INCLUDE_DIR)
 # The environment variable was not set
 SET(ERR_MSG "Please point the environment variable JSONCPP_INCLUDE_DIR to the root directory of your JSONCPP installation.")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for JSONCPP in ${DEF_DIR_JSONCPP}")
    IF(IS_DIRECTORY ${DEF_DIR_JSONCPP})
      MESSAGE(STATUS "Found!")
      SET(JSONCPP_INCLUDE_DIR ${DEF_DIR_JSONCPP})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSEIF(APPLE)
    MESSAGE(STATUS "Looking for JSONCPP in ${DEF_DIR_JSONCPP}")
    IF(IS_DIRECTORY ${DEF_DIR_JSONCPP})
      MESSAGE(STATUS "Found!")
      SET(JSONCPP_INCLUDE_DIR ${DEF_DIR_JSONCPP})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR ${ERR_MSG})
  ENDIF()
ENDIF()
INCLUDE_DIRECTORIES(${JSONCPP_INCLUDE_DIR}/include)
IF(WIN32)
	INCLUDE_DIRECTORIES(${JSONCPP_INCLUDE_DIR}/include)
	TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${JSONCPP_INCLUDE_DIR}/lib/jsoncpp.lib)
ELSEIF(APPLE)
	INCLUDE_DIRECTORIES(${JSONCPP_INCLUDE_DIR}/include)
	TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} -L${JSONCPP_INCLUDE_DIR}/lib/ libjsoncpp.a)
ENDIF()

# Get the MKL environment variable.
SET(MKL_DIR "$ENV{MKL_DIR}")
IF(NOT MKL_DIR)
  # The environment variable was not set
  #SET(ERR_MSG "Please point the environment variable MKL_DIR to the root directory of your MKL installation")
  IF(WIN32)
    # On Windows, try the default location
    MESSAGE(STATUS "Looking for MKL in ${DEF_DIR_MKL}")
    IF(IS_DIRECTORY ${DEF_DIR_MKL})
      MESSAGE(STATUS "Found!")
      SET(MKL_DIR ${DEF_DIR_MKL})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
    INCLUDE_DIRECTORIES(${MKL_DIR}/include)
    TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${MKL_DIR}/lib/intel64_win/mkl_core.lib ${MKL_DIR}/lib/intel64_win/mkl_intel_lp64.lib ${MKL_DIR}/lib/intel64_win/mkl_sequential.lib)
  ELSEIF(APPLE)
    MESSAGE(STATUS "Looking for MKL in ${DEF_DIR_MKL}")
    IF(IS_DIRECTORY ${DEF_DIR_MKL})
      MESSAGE(STATUS "Found!")
      SET(MKL_DIR ${DEF_DIR_MKL})
    ELSE()
      MESSAGE(FATAL_ERROR ${ERR_MSG})
    ENDIF()
    INCLUDE_DIRECTORIES(${MKL_DIR}/include)
    TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} -L${MKL_DIR}/lib/ libmkl_core.a libmkl_intel_lp64.a libmkl_sequential.a)
  ELSE()
    MESSAGE(FATAL_ERROR ${ERR_MSG})
  ENDIF()
ENDIF()


# # CMAES stuff
# SET(CMAES_INCLUDE_DIR "$ENV{CMAES_INCLUDE_DIR}")
# IF(NOT CMAES_INCLUDE_DIR)
  # # The environment variable was not set
  # SET(ERR_MSG "Please point the environment variable CMAES_INCLUDE_DIR to the root directory of your CMA installation.")
  # IF(WIN32)
    # # On Windows, try the default location
    # MESSAGE(STATUS "Looking for CMAES in ${CMAES_INCLUDE_DIR}")
    # IF(IS_DIRECTORY ${DEF_DIR_CMAES})
      # MESSAGE(STATUS "Found!")
      # SET(CMAES_INCLUDE_DIR ${DEF_DIR_CMAES})
    # ELSE()
      # MESSAGE(FATAL_ERROR ${ERR_MSG})
    # ENDIF()
  # ELSE()
    # MESSAGE(FATAL_ERROR ${ERR_MSG})
  # ENDIF()
# ENDIF()
# INCLUDE_DIRECTORIES(${CMAES_INCLUDE_DIR}/cma-es)


#Include liblbfgs-master
#SET(LINLBFGS_INCLUDE_DIR "$ENV{LINLBFGS_INCLUDE_DIR}")
#IF(NOT LINLBFGS_INCLUDE_DIR)
# # The environment variable was not set
# SET(ERR_MSG "Please point the environment variable LINLBFGS_INCLUDE_DIR to the root directory of your liblbfgs installation.")
#  IF(WIN32)
#    # On Windows, try the default location
#    MESSAGE(STATUS "Looking for liblbfgs in ${DEF_DIR_LIBLBFGS}")
#    IF(IS_DIRECTORY ${DEF_DIR_LIBLBFGS})
#      MESSAGE(STATUS "Found!")
#      SET(LINLBFGS_INCLUDE_DIR ${DEF_DIR_LIBLBFGS})
#    ELSE()
#      MESSAGE(FATAL_ERROR ${ERR_MSG})
#    ENDIF()
#  ELSE()
#    MESSAGE(FATAL_ERROR ${ERR_MSG})
#  ENDIF()
#ENDIF()
#INCLUDE_DIRECTORIES(${LINLBFGS_INCLUDE_DIR}/include)
#TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${LINLBFGS_INCLUDE_DIR}/x64/Debug/lbfgs_debug.lib)
#TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${LINLBFGS_INCLUDE_DIR}/x64/Release/lbfgs.lib)

#Include gurobi
#SET(GUROBI_INCLUDE_DIR "$ENV{GUROBI_INCLUDE_DIR}")
#IF(NOT GUROBI_INCLUDE_DIR)
# # The environment variable was not set
# SET(ERR_MSG "Please point the environment variable GUROBI_INCLUDE_DIR to the root directory of your gurobi installation.")
#  IF(WIN32)
#    # On Windows, try the default location
#    MESSAGE(STATUS "Looking for gurobi in ${DEF_DIR_GUROBI}")
#    IF(IS_DIRECTORY ${DEF_DIR_GUROBI})
#      MESSAGE(STATUS "Found!")
#      SET(GUROBI_INCLUDE_DIR ${DEF_DIR_GUROBI})
#    ELSE()
#      MESSAGE(FATAL_ERROR ${ERR_MSG})
#    ENDIF()
#  ELSE()
#    MESSAGE(FATAL_ERROR ${ERR_MSG})
#  ENDIF()
#ENDIF()
#INCLUDE_DIRECTORIES(${GUROBI_INCLUDE_DIR}/include)
#TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} ${GUROBI_INCLUDE_DIR}/lib/gurobi_c++md2017.lib)


#Include yixuan-LBFGSapp
#SET(YIXUAN_INCLUDE_DIR "$ENV{YIXUAN_INCLUDE_DIR}")
#IF(NOT YIXUAN_INCLUDE_DIR)
# # The environment variable was not set
# SET(ERR_MSG "Please point the environment variable YIXUAN_INCLUDE_DIR to the root directory of your yixuan installation.")
#  IF(WIN32)
#    # On Windows, try the default location
#    MESSAGE(STATUS "Looking for yixuan in ${DEF_DIR_YIXUAN}")
#    IF(IS_DIRECTORY ${DEF_DIR_YIXUAN})
#      MESSAGE(STATUS "Found!")
#      SET(YIXUAN_INCLUDE_DIR ${DEF_DIR_YIXUAN})
#    ELSE()
#      MESSAGE(FATAL_ERROR ${ERR_MSG})
#    ENDIF()
#  ELSE()
#    MESSAGE(FATAL_ERROR ${ERR_MSG})
#  ENDIF()
#ENDIF()
#INCLUDE_DIRECTORIES(${YIXUAN_INCLUDE_DIR}/include)

#Include dlib
#SET(DLIB_INCLUDE_DIR "$ENV{DLIB_INCLUDE_DIR}")
#IF(NOT DLIB_INCLUDE_DIR)
# # The environment variable was not set
# SET(ERR_MSG "Please point the environment variable DLIB_INCLUDE_DIR to the root directory of your dlib installation.")
#  IF(WIN32)
#    # On Windows, try the default location
#    MESSAGE(STATUS "Looking for dlib in ${DEF_DIR_DLIB}")
#    IF(IS_DIRECTORY ${DEF_DIR_DLIB})
#      MESSAGE(STATUS "Found!")
#      SET(DLIB_INCLUDE_DIR ${DEF_DIR_DLIB})
#    ELSE()
#      MESSAGE(FATAL_ERROR ${ERR_MSG})
#    ENDIF()
#  ELSE()
#    MESSAGE(FATAL_ERROR ${ERR_MSG})
#  ENDIF()
#ENDIF()
#ADD_SUBDIRECTORY(${DLIB_INCLUDE_DIR}/dlib dlib_build)
#TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} dlib::dlib)


# OS specific options and libraries
IF(WIN32)
  # c++11 is enabled by default.
  # -Wall produces way too many warnings.
  # -pedantic is not supported.
  # Disable warning 4996.
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4996")
  TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} opengl32.lib)
ELSE()
  # Enable all pedantic warnings.
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic")
  IF(APPLE)
    # Add required frameworks for GLFW.
    TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} "-framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo")
  ELSE()
    #Link the Linux OpenGL library
    TARGET_LINK_LIBRARIES(${CMAKE_PROJECT_NAME} "GL")
  ENDIF()
ENDIF()
