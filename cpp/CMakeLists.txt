cmake_minimum_required(VERSION 3.12.0) # older would work, but could give warnings on policy CMP0074
project(petsird VERSION 0.2.0)

set(CMAKE_CXX_STANDARD 17)

if(WIN32)
  add_compile_options(/W3 /WX)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

add_executable(petsird_generator petsird_generator.cpp)
target_link_libraries(petsird_generator petsird_generated)

add_executable(petsird_analysis petsird_analysis.cpp)
target_link_libraries(petsird_analysis petsird_generated)

add_subdirectory(generated)

# PETSIRD Library (libpetsird)
add_library(petsird SHARED $<TARGET_OBJECTS:petsird_generated>)
target_link_libraries(petsird petsird_generated)
install(TARGETS petsird DESTINATION lib)

# Set the target include directory where headers will be installed
set(INSTALL_INCLUDE_DIR "include/petsird")

# Collect all header files in the current directory and subdirectories
file(GLOB_RECURSE HEADER_FILES "*.h")

# Specify the installation directory for header files
install(FILES ${HEADER_FILES}
    DESTINATION ${INSTALL_INCLUDE_DIR}
)
