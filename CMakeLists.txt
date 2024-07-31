cmake_minimum_required(VERSION 3.20)
project(fusion_io VERSION 1.0.0 LANGUAGES C CXX Fortran)

enable_testing()
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(FindPkgConfig)

pkg_check_modules(HDF5 REQUIRED IMPORTED_TARGET hdf5)

set(CMAKE_Fortran_FLAGS "-fdefault-real-8 -fimplicit-none -fPIC")
set(CMAKE_CXX_FLAGS "-fPIC")
set(CMAKE_C_FLAGS "-fPIC")

add_subdirectory(m3dc1_lib)
add_subdirectory(fusion_io)
add_subdirectory(trace)

option(BUILD_TESTING "Build for CTest" OFF)
if (BUILD_TESTING)
    add_subdirectory(examples)
endif()