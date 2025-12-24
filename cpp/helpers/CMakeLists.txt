find_package(xtensor-blas REQUIRED)

# create a petsird_helpers library target
# Currently it is just empty, but allows creating target properties which are transitive
add_library(petsird_helpers INTERFACE)
target_link_libraries(petsird_helpers INTERFACE petsird xtensor-blas)
target_include_directories(petsird_helpers
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

add_executable(petsird_generator petsird_generator.cpp)
target_link_libraries(petsird_generator petsird_helpers)

add_executable(petsird_analysis petsird_analysis.cpp)
target_link_libraries(petsird_analysis petsird_helpers)
