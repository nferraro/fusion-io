
add_library(fusionio_trace
    trace_integrator.cpp
)
target_include_directories(fusionio_trace PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(fusionio_trace PUBLIC fusionio::fusionio)
target_include_directories(fusionio_trace PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_library(fusionio::trace ALIAS fusionio_trace)

find_package(MPI REQUIRED)
add_executable(trace
    main.cpp
)
target_link_libraries(trace PRIVATE fusionio::trace MPI::MPI_CXX)

find_package(LAPACK REQUIRED)
add_executable(plot_dr plot_dr.c)
target_link_libraries(plot_dr PRIVATE LAPACK::LAPACK)

install(TARGETS trace plot_dr
   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
