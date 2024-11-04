cmake_minimum_required(VERSION 3.12.0) # older would work, but could give warnings on policy CMP0074
project(petsird_test VERSION 0.2.0)

set(CMAKE_CXX_STANDARD 17)

if(WIN32)
  add_compile_options(/W3 /WX)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

find_package(petsird REQUIRED)
