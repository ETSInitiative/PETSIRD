#  Copyright (C) 2025 UMCG
#  Copyright (C) 2023-2026 University College London
#  SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.12.0)
cmake_policy(SET CMP0074 NEW) # find_package() uses <PackageName>_ROOT variables

project(PETSIRD VERSION 0.13.3 LANGUAGES CXX)
include(GNUInstallDirs)
set(PETSIRD_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/PETSIRD-${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}")

include(CMakePackageConfigHelpers)

# Not required, as this is set in the yardl-generated file
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# make sure that shared libraries can link against this static one
set(CMAKE_POSITION_INDEPENDENT_CODE True)

if(WIN32)
  add_compile_options(/W2)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  message(STATUS "ccache found, so we will use it.")
endif()

option(USE_HDF5 "Enable PETSIRD HDF5 support" ON)
# set internal option in the generated CMakeLists.txt from our own
# (as we can't get access to the value of the internal one)
if (CMAKE_VERSION VERSION_LESS 4.2)
  set(PETSIRD_GENERATED_USE_HDF5  ${USE_HDF5} CACHE BOOL "Use HDF5")
else()
  set(CACHE{PETSIRD_GENERATED_USE_HDF5} TYPE BOOL HELP "Use HDF5" VALUE ${USE_HDF5})
endif()

# assume that yardl was run with the following in the _package.yml
#   sourcesOutputDir: ../cpp/generated/petsird
add_subdirectory(generated/petsird)

mark_as_advanced(PETSIRD_GENERATED_USE_HDF5)

target_include_directories(petsird_generated
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/generated>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# petsird_generated is an OBJECT library. Those are a bit hard
# to handle for external projects, so create a "real" library as well.
# This is the one people should use.
add_library(petsird STATIC)
target_link_libraries(petsird PUBLIC petsird_generated)

add_subdirectory(helpers)

install(TARGETS petsird petsird_generated petsird_helpers
  EXPORT PETSIRDTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(EXPORT PETSIRDTargets
  NAMESPACE PETSIRD::
  DESTINATION ${PETSIRD_CMAKE_DIR}
)

install(
  DIRECTORY generated/petsird/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/petsird
)

install(
  DIRECTORY helpers/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
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
