
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
