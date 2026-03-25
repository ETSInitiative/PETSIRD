find_package(xtensor-blas REQUIRED)

# create a petsird_helpers library target
# Currently it is just empty, but allows creating target properties which are transitive
add_library(petsird_helpers INTERFACE)
# Future change for when we have sources:
# add_library(petsird_helpers ${PETSIRD_HELPERS_SOURCES})

target_include_directories(petsird_helpers
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(petsird_helpers
  INTERFACE
    petsird
    xtensor-blas
)

add_executable(petsird_generator petsird_generator.cpp)
target_link_libraries(petsird_generator PRIVATE petsird_generated petsird_helpers)

add_executable(petsird_analysis petsird_analysis.cpp)
target_link_libraries(petsird_analysis PRIVATE petsird_generated petsird_helpers)

