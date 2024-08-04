cmake_minimum_required(VERSION 3.10)

# Set the project name
project(trace LANGUAGES CXX C)

# Set the source directory
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
   message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

# # Include platform-specific options
# include("${FIO_ROOT}/install/make.inc.${FIO_ARCH}")

# # Set the libraries and include directories
# set(LIBS
#     "${FIO_ROOT}/fusion_io/_${FIO_ARCH}"
#     "${FIO_ROOT}/m3dc1_lib/_${FIO_ARCH}"
#     ${HDF5_LIBS}
# )

set(INCLUDE_DIRS
    "../fusion_io"
    "../m3dc1_lib"
)

# Add include directories
include_directories(${INCLUDE_DIRS})

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Set the source files
set(SOURCES
    trace_integrator.cpp
    main.cpp
)

find_package(HDF5 REQUIRED)

# Create the trace executable
add_executable(trace ${SOURCES})
target_link_libraries(trace PUBLIC fusionio_static m3dc1_static HDF5::HDF5 MPI::MPI_C)

# Installation
install(TARGETS trace
    RUNTIME DESTINATION "bin"
)

# # Define the installation directory
# install(DIRECTORY DESTINATION "${FIO_INSTALL_DIR}/bin")

# # Include common footer
# include("${FIO_ROOT}/install/common_footer.mk")
