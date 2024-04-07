
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
