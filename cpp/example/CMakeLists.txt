cmake_minimum_required(VERSION 3.12.0)
project(PETSIRDTest VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
find_package(PETSIRD REQUIRED CONFIG)
add_executable(PETSIRD_test PETSIRD_test.cpp)
target_link_libraries(PETSIRD_test PUBLIC PETSIRD::petsird)
add_executable(PETSIRD_test_helpers PETSIRD_test_helpers.cpp)
target_link_libraries(PETSIRD_test_helpers PUBLIC PETSIRD::helpers)

# install(TARGETS PETSIRD_test DESTINATION bin)
