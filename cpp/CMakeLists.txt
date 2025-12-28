cmake_minimum_required(VERSION 3.12.0) # older would work, but could give warnings on policy CMP0074
project(petsird VERSION 0.7.2 LANGUAGES CXX)
include(GNUInstallDirs)
set(PETSIRD_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/PETSIRD")

include(CMakePackageConfigHelpers)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(WIN32)
  add_compile_options(/W3 /WX)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  message(STATUS "ccache found, so we will use it.")
endif()

# assume that yardl was run with the following in the _package.yml
#   sourcesOutputDir: ../cpp/generated/petsird

add_subdirectory(generated/petsird)
add_subdirectory(helpers)
# create a petsird library target
# Currently it is just empty, but allows creating target properties which are transitive
add_library(petsird)
add_library(PETSIRD::petsird ALIAS petsird)

target_link_libraries(petsird
  PUBLIC
    petsird_generated
    petsird_helpers
)

target_compile_features(petsird PUBLIC cxx_std_17)

target_include_directories(petsird
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/generated>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS petsird petsird_generated petsird_helpers
  EXPORT petsirdTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/generated/petsird/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/petsird
)

install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/helpers/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT petsirdTargets
  NAMESPACE PETSIRD::
  DESTINATION ${PETSIRD_CMAKE_DIR}
)


configure_package_config_file(
  cmake/PETSIRDConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/PETSIRDConfig.cmake
  INSTALL_DESTINATION ${PETSIRD_CMAKE_DIR}
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/PETSIRDConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/PETSIRDConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/PETSIRDConfigVersion.cmake
  DESTINATION ${PETSIRD_CMAKE_DIR}
)
