cmake_minimum_required(VERSION 3.0.0)
project(prd VERSION 0.0.1)

set(CMAKE_CXX_STANDARD 17)

if(WIN32)
  add_compile_options(/W3 /WX)
else()
  add_compile_options(-Wall -Wextra -pedantic -Werror)
endif()

add_executable(prd_generator prd_generator.cpp)
target_link_libraries(prd_generator prd_generated)

add_subdirectory(generated)
