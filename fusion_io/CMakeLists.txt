set(FUSIONIO_HEADERS
        coil_source.h
        compound_field.h
        diiid_coils.h
        fio_operations.h
        fusion_io_c.h
        fusion_io_defs.h
        fusion_io_field.h
        fusion_io_series.h
        fusion_io_source.h
        fusion_io_species.h
        fusion_io.h
        gato_field.h
        gato_source.h
        geqdsk_field.h
        geqdsk_source.h
        gpec_field.h
        gpec_source.h
        interpolate.h
        isosurface.h
        m3dc1_field.h
        m3dc1_source.h
        mars_field.h
        mars_source.h
        options.h
)

add_library(fusionio_fusionio
        interpolate.cpp
        options.cpp
        fusion_io_species.cpp
        fusion_io_field.cpp
        fusion_io_series.cpp
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

set_target_properties(fusionio_fusionio PROPERTIES 
        PUBLIC_HEADER "${FUSIONIO_HEADERS}"
        Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(fusionio_fusionio PUBLIC 
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:$<TARGET_PROPERTY:fusionio_fusionio,Fortran_MODULE_DIRECTORY>>
)
target_link_libraries(fusionio_fusionio PUBLIC fusionio_m3dc1)
target_compile_definitions(fusionio_fusionio PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>: FORTRAN>)
add_library(fusionio::fusionio ALIAS fusionio_fusionio)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/fusion_io.mod TYPE INCLUDE)


add_library(fusionio_push push_interface.f90)
target_link_libraries(fusionio_push PUBLIC fusionio::fusionio)
set_target_properties(fusionio_push PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(fusionio_push PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:$<TARGET_PROPERTY:fusionio_push,Fortran_MODULE_DIRECTORY>>
)
target_compile_definitions(fusionio_push PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>: FORTRAN>)
add_library(fusionio::push ALIAS fusionio_push)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/fio_push.mod TYPE INCLUDE)


if (FUSIONIO_ENABLE_PYTHON)
        find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
        add_library(fio_py SHARED python_interface.cpp)
        target_link_libraries(fio_py PUBLIC fusionio::fusionio Python3::Python)

        # note this will fail on windows. If building on windows is a priority
        # we can configure a setup.py and use add_custom_target to build the module
        set_target_properties(fio_py PROPERTIES PREFIX "" SUFFIX ".so")

        install(FILES fpy.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

        if(DEFINED PYTHON_MODULE_INSTALL_PATH)
                install(FILES fpy.py DESTINATION ${PYTHON_MODULE_INSTALL_PATH})
                install(TARGETS fio_py DESTINATION ${PYTHON_MODULE_INSTALL_PATH})
        endif()
endif()

install(TARGETS fusionio_fusionio 
        EXPORT fusionio-target)
install(EXPORT fusionio-target
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fusionio
        NAMESPACE fusionio::
        FILE fusionio-target.cmake
        COMPONENT fusionio-target)
install(FILES fusionio-config.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fusionio)
write_basic_package_version_file(
        "${PROJECT_BINARY_DIR}/fusionio-config-version.cmake"
        COMPATIBILITY AnyNewerVersion)
install(FILES
        "${PROJECT_BINARY_DIR}/fusionio-config-version.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fusionio)