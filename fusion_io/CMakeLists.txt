cmake_minimum_required(VERSION 3.12)

# set(CMAKE_DISABLE_SOURCE_CHANGES ON)
# set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(fusion_io LANGUAGES CXX Fortran C)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# enable_language(Fortran)

set(INCLUDE_DIRS
    "../m3dc1_lib"
)

# Add include directories
include_directories(${INCLUDE_DIRS})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Set compiler options for Fortran
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DFORTRAN")

# Source files
set(SOURCES
    interpolate.cpp
    options.cpp
    fusion_io_species.cpp
    fusion_io_series.cpp
    fusion_io_field.cpp
    fio_operations.cpp
    compound_field.cpp
    isosurface.cpp
    m3dc1_source.cpp
    m3dc1_field.cpp
    mars_source.cpp
    mars_field.cpp
    geqdsk_source.cpp
    geqdsk_field.cpp
    gato_source.cpp
    gato_field.cpp
    gpec_source.cpp
    gpec_field.cpp
    fusion_io.cpp
    c_interface.cpp
    fusion_io_fortran.F90
)

set_property(
	SOURCE fusion_io_fortran.F90
    PROPERTY COMPILE_DEFINITIONS
	FORTRAN
)

# set(HDF5_USE_STATIC_LIBRARIES TRUE)
find_package(HDF5 REQUIRED)

# Library targets
add_library(fusionio_static STATIC ${SOURCES})
set_target_properties(fusionio_static PROPERTIES OUTPUT_NAME "fusionio")

# Shared library
add_library(fusionio_shared SHARED ${SOURCES})
set_target_properties(fusionio_shared PROPERTIES OUTPUT_NAME "fusionio")

# set_target_properties(fusionio PROPERTIES OUTPUT_NAME "fusionio")

# Link libraries
target_link_libraries(fusionio_static PUBLIC m3dc1_static HDF5::HDF5)
target_link_libraries(fusionio_shared PUBLIC m3dc1_shared HDF5::HDF5)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  target_compile_options(fusionio_static PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-s real64>)
  target_compile_options(fusionio_shared PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-s real64>)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(fusionio_static PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-real-8>)
  target_compile_options(fusionio_shared PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-real-8>)
else()
  target_compile_options(fusionio_static PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-r8>)
  target_compile_options(fusionio_shared PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-r8>)
endif()


# Python binding
if(ENABLE_PYTHON)
  find_package(Python COMPONENTS Development)
  add_library(fio_py SHARED fio_setup.py python_interface.cpp)
  set_target_properties(fio_py PROPERTIES PREFIX "")
  target_include_directories(fio_py PRIVATE ${Python_INCLUDE_DIRS})
  target_link_libraries(fio_py PRIVATE fusionio_shared m3dc1_shared ${Python_LIBRARIES})
endif()

# Install targets
install(TARGETS fusionio_static fusionio_shared
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

# Install headers and Python scripts
install(FILES fusion_io_defs.h fusion_io_c.h
        DESTINATION include)
if(ENABLE_PYTHON)
  if(NOT DEFINED PYTHON_MODULE_INSTALL_PATH)
    set(PYTHON_MODULE_INSTALL_PATH lib)
  endif()
  install(FILES fpy.py
          DESTINATION ${PYTHON_MODULE_INSTALL_PATH})
  install(TARGETS fio_py
          DESTINATION ${PYTHON_MODULE_INSTALL_PATH})
endif()
#
# # Installation script modifications
# # Note: Actual modification of scripts (e.g., using 'sed') should be done outside of CMake or via add_custom_command
#
# # Conda specific installation
# # Handle Conda installation specifics here or via a custom target
