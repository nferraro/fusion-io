cmake_minimum_required(VERSION 3.12)

# set(CMAKE_DISABLE_SOURCE_CHANGES ON)
# set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(m3dc1 LANGUAGES CXX Fortran C)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# enable_language(Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Source files
set(SOURCES
    m3dc1_mesh.cpp
    m3dc1_field.cpp
    m3dc1_coord_map.cpp
    m3dc1_stell_field.cpp
    m3dc1_timeslice.cpp
    m3dc1_file.cpp
)

find_package(HDF5 REQUIRED)

# Object library for common object files
# add_library(m3dc1_objs OBJECT ${SOURCES})
# target_link_libraries(m3dc1_objs HDF5::HDF5)

# Static library
# add_library(m3dc1_static STATIC $<TARGET_OBJECTS:m3dc1_objs>)
add_library(m3dc1_static STATIC ${SOURCES})
set_target_properties(m3dc1_static PROPERTIES OUTPUT_NAME "m3dc1")

# Shared library
# add_library(m3dc1_shared SHARED $<TARGET_OBJECTS:m3dc1_objs>)
add_library(m3dc1_shared SHARED ${SOURCES})
set_target_properties(m3dc1_shared PROPERTIES OUTPUT_NAME "m3dc1")

# Fortran static library
add_library(m3dc1_fortran STATIC m3dc1_fortran.cpp)

# Linking libraries
target_link_libraries(m3dc1_static PRIVATE HDF5::HDF5)
target_link_libraries(m3dc1_shared PRIVATE HDF5::HDF5)
target_link_libraries(m3dc1_fortran PRIVATE m3dc1 HDF5::HDF5)

# # Fortran executable
# add_executable(m3dc1_fortran_test test.f90)
# target_link_libraries(m3dc1_fortran_test PRIVATE m3dc1_fortran m3dc1_static)
#
# Install targets
install(TARGETS m3dc1_static m3dc1_shared m3dc1_fortran
        ARCHIVE DESTINATION "lib"
        LIBRARY DESTINATION "lib")
