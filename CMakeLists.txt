cmake_minimum_required(VERSION 2.8)

set(gmx_clusterByFeatures_MAJOR_VERSION 1)
set(gmx_clusterByFeatures_MINOR_VERSION 0)
set(gmx_clusterByFeatures_PATCH_VERSION 0)
set(gmx_clusterByFeatures_VERSION
  ${gmx_clusterByFeatures_MAJOR_VERSION}.${gmx_clusterByFeatures_MINOR_VERSION}.${gmx_clusterByFeatures_PATCH_VERSION})

## Added for ECLIPSE ####
# SET(CMAKE_VERBOSE_MAKEFILE ON)

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
	set(DEF_INSTALL_CMAKE_DIR CMake)
else()
set(DEF_INSTALL_CMAKE_DIR lib/CMake/gmx_clusterByFeatures)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
	set(var INSTALL_${p}_DIR)
	if(NOT IS_ABSOLUTE "${${var}}")
		set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
	endif()
endforeach()

# Cmake modules/macros are in a subdirectory to keep this file cleaner
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

##############################################################
################### PYTHON CHECKING #########################
##############################################################
macro(check_python_module module)
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c" "import ${module}; print(${module}.__version__)" RESULT_VARIABLE FOUND_MODULE
                    OUTPUT_VARIABLE MODULE_VERSION
                    ERROR_QUIET)
    if (NOT ${FOUND_MODULE} EQUAL 0)
            message(FATAL_ERROR "\n ${module} Python module not found.... \n Please install ${module} by following command:\n sudo pip3 install ${module} \n or \n sudo pip3 install ${module} \n")
    else (NOT ${FOUND_MODULE} EQUAL 0)
        string(STRIP ${MODULE_VERSION} MODULE_VERSION)
        message(STATUS "${module} (version ${MODULE_VERSION}) Python module found.")
    endif(NOT ${FOUND_MODULE} EQUAL 0)
endmacro(check_python_module module)

set(Python_ADDITIONAL_VERSIONS 3.4 3.5 3.6 3.7)
find_package(PythonInterp)
if (NOT PYTHONINTERP_FOUND)
        message(FATAL_ERROR "Python3 not found... Please install python3. ")
endif(NOT PYTHONINTERP_FOUND)
find_package(PythonLibs)
if (NOT PythonLibs_FOUND)
        message(FATAL_ERROR "Python developement library not found. On Ubuntu-like OS install python3-dev. On fedora/centos/RHEL install python3-devel.")
endif(NOT PythonLibs_FOUND)
check_python_module("sklearn")
check_python_module("matplotlib")
include_directories( ${PYTHON_INCLUDE_DIRS} )


##############################################################
################### CHECKING GROMACS PATH ####################
##############################################################
find_package(GROMACS COMPONENTS "libgromacs")
if (NOT GROMACS_FOUND)
	message(STATUS "Continue checking for other Gromacs version library file")
endif(NOT GROMACS_FOUND)

if (NOT GROMACS_FOUND)
	message(FATAL_ERROR "\nGromacs path is not defined.....\nUSE:\n -DGMX_PATH=/path/to/gromacs\n  or  \n export CMAKE_PREFIX_PATH=\"/path/to/gromacs\"\n" )
endif(NOT GROMACS_FOUND)

message(STATUS "Found Gromacs include directory: " ${GROMACS_INCLUDE_DIRS})

if(${GROMACS_VERSION} VERSION_LESS "5.2")
	message(FATAL_ERROR "Gromacs-${GROMACS_VERSION} not supported. Please use Gromacs-2016 and above versions.\n")
endif(${GROMACS_VERSION} VERSION_LESS "5.2")

include_directories( ${GROMACS_INCLUDE_DIRS} )
add_definitions( ${GROMACS_DEFINITIONS} )
##############################################################

##############################################################
# checking for C++11 or newer
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(FATAL_ERROR "C++11 complier not found. Please use a newer compiler, or a newer standard library.")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")


##############################################################
# Add pybind11
set(PYBIND11_CPP_STANDARD -std=c++11)
add_subdirectory(pybind11)

##############################################################
# Python module compilation and installation path
if (INPLACE)
    set(PYTHON_INSTDIR "${CMAKE_SOURCE_DIR}")
    message(STATUS "gmx_clusterByFeatures will be built in place and executed from here")
else (INPLACE)
    # Determine correct installation directory for Python bindings
    # From http://plplot.svn.sourceforge.net/viewvc/plplot/trunk/cmake/modules/python.cmake?revision=11014
    execute_process(
        COMMAND
        ${PYTHON_EXECUTABLE} -c "import site; print(site.getsitepackages()[0])"
        OUTPUT_VARIABLE PYTHON_INSTDIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
endif (INPLACE)
message(STATUS "Python Install directory ${PYTHON_INSTDIR}")
set(gmx_clusterByFeatures_PYTHON_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/gmx_clusterByFeatures")
configure_file("gmx_clusterByFeatures.py.in" "gmx_clusterByFeatures.py" @ONLY)

# Change python-code file to C hex-code
execute_process(COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/pyCode2Hex.py ${CMAKE_SOURCE_DIR}/src/cluster.py ${CMAKE_SOURCE_DIR}/src/cluster.pyhex)


# Add sub-directories
add_subdirectory(src)

if(NOT INPLACE)
    INSTALL(DIRECTORY "${CMAKE_SOURCE_DIR}/gmx_clusterByFeatures" DESTINATION ${PYTHON_INSTDIR}
    COMPONENT python
    PATTERN ".svn" EXCLUDE
    PATTERN "test_data" EXCLUDE
    PATTERN "testData" EXCLUDE
    PATTERN "test_list*" EXCLUDE
    PATTERN "CMake*" EXCLUDE
    PATTERN "Basement" EXCLUDE
    PATTERN "UnitTest*" EXCLUDE
    PATTERN "__pycache__*" EXCLUDE
    PATTERN ".localhistory*" EXCLUDE
    )
endif(NOT INPLACE)

INSTALL(PROGRAMS "${CMAKE_BINARY_DIR}/gmx_clusterByFeatures.py" DESTINATION "${INSTALL_BIN_DIR}" RENAME "gmx_clusterByFeatures" COMPONENT bin)
