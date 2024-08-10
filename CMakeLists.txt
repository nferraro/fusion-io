cmake_minimum_required(VERSION 3.24)
project(fusion_io VERSION 1.0.0 LANGUAGES C CXX Fortran)

enable_testing()
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(FindPkgConfig)
pkg_check_modules(HDF5 REQUIRED IMPORTED_TARGET hdf5)

option(FUSIONIO_ENABLE_PYTHON "Enable Python" ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "-fdefault-real-8")
endif()

add_subdirectory(m3dc1_lib)
add_subdirectory(fusion_io)

option(FUSIONIO_ENABLE_TRACE "Enable Trace" ON)
if (FUSIONIO_ENABLE_TRACE)
    add_subdirectory(trace)
endif()

option(BUILD_TESTING "Build for CTest" OFF)
if (BUILD_TESTING)
    add_subdirectory(examples)
endif()