cmake_minimum_required(VERSION 3.20)
project(fusion_io VERSION 1.0.0 LANGUAGES C CXX Fortran)
enable_testing()
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

find_package(HDF5 REQUIRED)

add_subdirectory(m3dc1_lib)
add_subdirectory(fusion_io)
add_subdirectory(trace)

option(BUILD_TESTING "Build for CTest" OFF)
if (BUILD_TESTING)
    add_subdirectory(examples)
endif()