cmake_minimum_required(VERSION 3.12.0) # older would work, but could give warnings on policy CMP0074
project(petsird VERSION 0.7.2)

set(CMAKE_CXX_STANDARD 17)

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
# create a petsird library target
# Currently it is just empty, but allows creating target properties which are transitive
add_library(petsird)
target_link_libraries(petsird PUBLIC petsird_generated)
target_include_directories(petsird PUBLIC "${PROJECT_SOURCE_DIR}/generated")

add_subdirectory(helpers)
