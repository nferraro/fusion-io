file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_executable(example example.cpp)
target_link_libraries(example PRIVATE fusionio::fusionio)

add_executable(example_c example_c.c)
target_link_libraries(example_c PRIVATE fusionio::fusionio)

add_executable(example_fortran example_fortran.f90)
target_link_libraries(example_fortran PRIVATE fusionio::fusionio)
target_compile_options(example_fortran PRIVATE -fdefault-real-8 -fimplicit-none)

add_executable(example_push example_push.f90)
target_link_libraries(example_push PRIVATE fusionio::push)
target_compile_options(example_push PRIVATE -fdefault-real-8 -fimplicit-none)


add_test(cpp_example_m3dc1 ./example m3dc1)
add_test(cpp_example_geqdsk ./example geqdsk)
add_test(cpp_example_gpec ./example gpec)

add_test(c_example ./example_c)
add_test(fortran_example ./example_fortran)