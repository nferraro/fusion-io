include(CMakeFindDependencyMacro)
find_dependency(HDF5)

include(${CMAKE_CURRENT_LIST_DIR}/fusionio_m3dc1-target.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/fusionio-target.cmake)

check_required_components(fusionio)