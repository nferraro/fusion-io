cmake_minimum_required(VERSION 3.20)
project(fusion_io C CXX Fortran)

find_package(HDF5 REQUIRED)

add_subdirectory(m3dc1_lib)
add_subdirectory(fusion_io)
add_subdirectory(trace)
add_subdirectory(examples)