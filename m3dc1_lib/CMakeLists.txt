add_library(fusionio_m3dc1
        m3dc1_mesh.cpp
        m3dc1_field.cpp
        m3dc1_coord_map.cpp
        m3dc1_stell_field.cpp
        m3dc1_timeslice.cpp
        m3dc1_file.cpp
)
target_include_directories(fusionio_m3dc1 PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(fusionio_m3dc1 PUBLIC PkgConfig::HDF5)
target_compile_options(fusionio_m3dc1 PRIVATE -fPIC)
add_library(fusionio::m3dc1 ALIAS fusionio_m3dc1)

add_library(fusionio_m3dc1_fortran
            m3dc1_fortran.cpp)
target_link_libraries(fusionio_m3dc1_fortran PUBLIC fusionio::m3dc1)
add_library(fusionio::m3dc1_fortran ALIAS fusionio_m3dc1_fortran)

# tests
add_executable(m3dc1_fortran_test test.f90)
target_link_libraries(m3dc1_fortran_test PRIVATE fusionio::m3dc1_fortran)
target_compile_options(m3dc1_fortran_test PRIVATE -fdefault-real-8 -fimplicit-none)
