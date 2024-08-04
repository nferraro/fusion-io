cmake_minimum_required(VERSION 3.8)
project(fusion-io LANGUAGES Fortran C CXX)

cmake_host_system_information(RESULT HOST QUERY HOSTNAME)
message(STATUS "Configuring to build on: ${HOST}")

OPTION(ENABLE_TRACE "Build trace program" ON)
OPTION(ENABLE_PYTHON "Build Python binding" ON)

message(STATUS "Configuring m3dc1_lib")
add_subdirectory(m3dc1_lib)
message(STATUS "Configuring fusion_io")
add_subdirectory(fusion_io)
if(ENABLE_TRACE)
  message(STATUS "Configuring trace")
  add_subdirectory(trace)
endif()

