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
        m3dc1_numpy_source.h
        m3dc1_source.h
        mars_field.h
        mars_source.h
        nimrod_field.h
        nimrod_source.h
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
        m3dc1_numpy_source.cpp
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
        nimrod_source.cpp
        nimrod_field.cpp
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

if (FUSIONIO_ENABLE_NIMROD)
    target_compile_definitions(fusionio_fusionio PUBLIC FUSIONIO_ENABLE_NIMROD)
    target_link_libraries(fusionio_fusionio PUBLIC ${NIMEVAL_LIBRARY} ${NIMLIB_LIBRARY} ${FFTW3_LIBRARIES})
endif()
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

        # Build fio_py.so by compiling all sources directly into the shared
        # library rather than linking against static archives.  This ensures
        # that RTTI typeinfo and vtable symbols (e.g. for m3dc1_numpy_source,
        # needed by dynamic_cast in python_interface.cpp) are present in
        # fio_py.so at runtime.

        # Collect sources from fusionio_fusionio.
        get_target_property(_fusionio_sources fusionio_fusionio SOURCES)
        get_target_property(_fusionio_src_dir fusionio_fusionio SOURCE_DIR)

        # Collect sources from fusionio_m3dc1 (m3dc1_lib).
        get_target_property(_m3dc1_sources fusionio_m3dc1 SOURCES)
        get_target_property(_m3dc1_src_dir fusionio_m3dc1 SOURCE_DIR)

        # Resolve relative source paths to absolute paths.
        set(_fio_py_sources python_interface.cpp)

        foreach(_src ${_fusionio_sources})
            if(IS_ABSOLUTE "${_src}")
                list(APPEND _fio_py_sources "${_src}")
            else()
                list(APPEND _fio_py_sources "${_fusionio_src_dir}/${_src}")
            endif()
        endforeach()

        foreach(_src ${_m3dc1_sources})
            if(IS_ABSOLUTE "${_src}")
                list(APPEND _fio_py_sources "${_src}")
            else()
                list(APPEND _fio_py_sources "${_m3dc1_src_dir}/${_src}")
            endif()
        endforeach()

        # Remove Fortran sources — they cannot be compiled into a Python .so.
        list(FILTER _fio_py_sources EXCLUDE REGEX "\\.F90$|\\.f90$|\\.f$|\\.F$")

        # Deduplicate (m3dc1_field.cpp etc. may appear in both lists).
        list(REMOVE_DUPLICATES _fio_py_sources)

        add_library(fio_py SHARED ${_fio_py_sources})

        # Inherit include directories from both upstream targets.
        target_include_directories(fio_py PRIVATE
            $<TARGET_PROPERTY:fusionio_fusionio,INCLUDE_DIRECTORIES>
            $<TARGET_PROPERTY:fusionio_m3dc1,INCLUDE_DIRECTORIES>
        )

        # Inherit compile definitions from fusionio_fusionio.
        target_compile_definitions(fio_py PRIVATE
            $<TARGET_PROPERTY:fusionio_fusionio,COMPILE_DEFINITIONS>
        )

        # With Intel oneAPI (icpx), RTTI typeinfo symbols can be hidden even
        # when sources are compiled directly into the .so.  Force default
        # symbol visibility, enable RTTI explicitly, export all symbols to the
        # dynamic symbol table, and link libstdc++ explicitly so that the C++
        # runtime type_info infrastructure is fully resolved at load time.
        # These flags are harmless no-ops on GCC/Clang but are required on
        # Intel oneAPI to avoid "undefined symbol: _ZTI..." errors at import.
        if(CMAKE_CXX_COMPILER_ID MATCHES "Intel" OR
           CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
            target_compile_options(fio_py PRIVATE
                $<$<COMPILE_LANGUAGE:CXX>:-fvisibility=default>
                $<$<COMPILE_LANGUAGE:CXX>:-frtti>
            )
            target_link_options(fio_py PRIVATE
                -Wl,--export-dynamic
            )
            target_link_libraries(fio_py PRIVATE stdc++)
        endif()

        # Collect the transitive link libraries needed by fusionio_m3dc1
        # (e.g. HDF5) so fio_py.so can resolve those symbols.  We do NOT
        # link against fusionio_fusionio or fusionio_m3dc1 themselves because
        # their sources are compiled directly into fio_py above.
        get_target_property(_m3dc1_link_libs fusionio_m3dc1 LINK_LIBRARIES)

        # Also collect any extra link libraries from fusionio_fusionio
        # (e.g. NIMROD libraries) that are not already covered.
        get_target_property(_fusionio_link_libs fusionio_fusionio LINK_LIBRARIES)

        # Filter out the fusionio targets themselves from the link list to
        # avoid circular / duplicate symbol issues.
        set(_fio_py_link_libs)
        foreach(_lib ${_m3dc1_link_libs} ${_fusionio_link_libs})
            if(NOT "${_lib}" MATCHES "fusionio_fusionio|fusionio_m3dc1|fusionio::fusionio|fusionio::m3dc1")
                list(APPEND _fio_py_link_libs "${_lib}")
            endif()
        endforeach()
        list(REMOVE_DUPLICATES _fio_py_link_libs)

        target_link_libraries(fio_py PRIVATE ${_fio_py_link_libs} Python3::Python)

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
