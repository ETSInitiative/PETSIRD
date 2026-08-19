#  Copyright (C) 2026 University College London
#  SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.12.0)
cmake_policy(SET CMP0074 NEW) # find_package() uses <PackageName>_ROOT variables

project(PETSIRDtest VERSION 0.1 LANGUAGES CXX)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  message(STATUS "ccache found, so we will use it.")
endif()

find_package(PETSIRD CONFIG REQUIRED)

add_executable(minimal_test minimal_test.cpp)
target_link_libraries(minimal_test PETSIRD::petsird)

add_executable(minimal_test_helpers minimal_test_helpers.cpp)
target_link_libraries(minimal_test_helpers   PETSIRD::petsird_helpers)
