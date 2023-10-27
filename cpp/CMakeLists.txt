cmake_minimum_required(VERSION 3.12.0) # older would work, but could give warnings on policy CMP0074
project(prd VERSION 0.1.1)

set(CMAKE_CXX_STANDARD 17)

if(WIN32)
  add_compile_options(/W3 /WX)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

add_executable(prd_generator prd_generator.cpp)
target_link_libraries(prd_generator prd_generated)

add_executable(prd_analysis prd_analysis.cpp)
target_link_libraries(prd_analysis prd_generated)

add_subdirectory(generated)
