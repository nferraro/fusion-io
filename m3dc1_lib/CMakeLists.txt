set(FUSIONIO_M3DC1_HEADERS
        m3dc1_coord_map.h
        m3dc1_field.h
        m3dc1_file.h
        m3dc1_mesh.h
        m3dc1_stell_field.h
        m3dc1_timeslice.h
)
add_library(fusionio_m3dc1
        m3dc1_mesh.cpp
        m3dc1_field.cpp
        m3dc1_coord_map.cpp
        m3dc1_stell_field.cpp
        m3dc1_timeslice.cpp
        m3dc1_file.cpp
)
set_target_properties(fusionio_m3dc1 PROPERTIES PUBLIC_HEADER "${FUSIONIO_M3DC1_HEADERS}")
target_include_directories(fusionio_m3dc1 PUBLIC 
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(fusionio_m3dc1 PUBLIC PkgConfig::HDF5)
add_library(fusionio::m3dc1 ALIAS fusionio_m3dc1)

add_library(fusionio_m3dc1_fortran
            m3dc1_fortran.cpp)
target_link_libraries(fusionio_m3dc1_fortran PUBLIC fusionio::m3dc1)
add_library(fusionio::m3dc1_fortran ALIAS fusionio_m3dc1_fortran)


install(TARGETS fusionio_m3dc1
        EXPORT fusionio_m3dc1-target)
install(EXPORT fusionio_m3dc1-target
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fusionio
        NAMESPACE fusionio::
        FILE fusionio_m3dc1-target.cmake
        COMPONENT fusionio_m3dc1-target)

# tests
add_executable(m3dc1_fortran_test test.f90)
target_link_libraries(m3dc1_fortran_test PRIVATE fusionio::m3dc1_fortran)
set_property(TARGET m3dc1_fortran_test PROPERTY LINKER_LANGUAGE Fortran)
