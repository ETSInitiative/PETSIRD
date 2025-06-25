cmake_minimum_required(VERSION 3.12.0) # older would work, but could give warnings on policy CMP0074
project(petsird VERSION 0.7.0)

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

add_executable(petsird_generator petsird_generator.cpp)
target_link_libraries(petsird_generator petsird_generated)
target_include_directories(petsird_generator PRIVATE "${PROJECT_SOURCE_DIR}/helpers/include")
# needed for generated PETSIRD files
target_include_directories(petsird_generator PRIVATE "${PROJECT_SOURCE_DIR}/")

add_executable(petsird_analysis petsird_analysis.cpp)
target_link_libraries(petsird_analysis petsird_generated)
target_include_directories(petsird_analysis PRIVATE "${PROJECT_SOURCE_DIR}/helpers/include")
# needed for generated PETSIRD files
target_include_directories(petsird_analysis PRIVATE "${PROJECT_SOURCE_DIR}/")

add_subdirectory(generated)
